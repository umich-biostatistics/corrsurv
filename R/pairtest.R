

#' Title for matched.survtest function. <10 words
#' 
#' Short description for function. What does it do? What does the output include?
#' Is this the main function or a helper?
#' 
#' @param x1 insert description and data type, default value if needed
#' @param delta1 insert description and data type, default value if needed
#' @param x2 insert description and data type, default value if needed
#' @param delta2 insert description and data type, default value if needed
#' @param n the integer number of correlated individuals and we assume they are
#' sorted as in the paper, so that correlated observations come first up to the nth 
#' observation and then uncorrelated observations follow.
#' @param tm insert description and data type, default value if needed
#' @param weights string, options are "left" and "right". Option "left" will give pf and 
#' yls weighting using upperlimit similar to logrank. "right" will use weights 
#' that are similar, but right continuous.
#' @param maxtau insert description and data type, default value if needed
#' @param fudge.factor insert description and data type, default value if needed
#'
#' @return insert description and data type and what it contains
#'
#' @seealso 
#'
#' @references Murray, Susan. Nonparametric Rank-Based Methods for Group
#' Sequential Monitoring of Paired Censored Survival Data. 2000. Biometrics,
#' 56, pp. 984-990.
#'
#' @examples
#' data(pairdata)
#' eyeresults = pairtest(x1 = pairdata$x1, delta1 = pairdata$delta1, 
#'                       x2 = pairdata$x2, delta2 = pairdata$delta2, n = 3711)
#' print(eyeresults)
#' 
#' @export 
#'

pairtest = function(x1, delta1, x2, delta2, n, tm = sort(unique(c(0,x1,x2))), 
                    weights = "left", maxtau = 100000, fudge.factor = 0)  {

  #weights=left will give pf and yls weighting using upperlimit similar to logrank
  #weights=right will use weights that are similar, but right continuous
  
  maxtau = max(tm)

  n1 = length(x1)
  n2 = length(x2)
  nstar = (n1*n2)/n

  nx = n  #number of matched pairs
  x1pair = x1[1:n]
  x2pair = x2[1:n]
  delta1pair = delta1[1:n]
  delta2pair = delta2[1:n]

  nt = length(tm)

  #dN marginal terms and Y marginal terms need to use all data:

  failframe1 = as.data.frame(cbind(x1, delta1))
  failinfo1 = getriskdat(failframe1, tm)
  failframe2 = as.data.frame(cbind(x2, delta2))
  failinfo2 = getriskdat(failframe2, tm)
  dN1 = failinfo1$dN
  Y1 = failinfo1$Y
  dN2 = failinfo2$dN
  Y2 = failinfo2$Y
  dbarN = dN1 + dN2
  Ybar = Y1 + Y2
  haz1 = dN1/Y1
  haz1[dN1 == 0] = 0
  haz2 = dN2/Y2
  haz2[dN2 == 0] = 0
  hazbar = dbarN/Ybar
  hazbar[is.na(hazbar)] = 0

  #still need to calculate Kaplan-Meier curves: can use dN stuff
  haz1complement = 1 - haz1
  haz2complement = 1 - haz2
  hazbarcomplement = 1 - hazbar
 
  #check these later, but are probably correct Kaplan Meier estimates
  km1 = cumprod(haz1complement)
  km2 = cumprod(haz2complement)
  km.pooled = cumprod(hazbarcomplement)
  km.left.cont = c(1, km.pooled[1:nt-1])

  junk = 1 - delta1
  censorframe1 = as.data.frame(cbind(x1, junk))
  censorinfo1 = getriskdat(censorframe1, tm)
  junk = 1 - delta2    #just fixed this (was absent)
  censorframe2 = as.data.frame(cbind(x2, junk))
  censorinfo2 = getriskdat(censorframe2, tm)
 
  h.left.cont1 = c(1, censorinfo1$kmsurv[1:nt-1])
  h.left.cont2 = c(1, censorinfo2$kmsurv[1:nt-1])
 
  pr1dNmat = matrix(0, nrow = nx, ncol = nt)
  pr2dNmat = matrix(0, nrow = nx, ncol = nt)
  pr1Ymat = matrix(0, nrow = nx, ncol = nt)
  pr2Ymat = matrix(0, nrow = nx, ncol = nt)

  #This is only among the paired indiv's

  dNinfo = 
    .Fortran("bivsurvmats", as.integer(nx), as.integer(nt), as.double(tm),
	           as.double(x1pair), as.double(delta1pair), as.double(x2pair), as.double(delta2pair),
	           matrix(as.double(pr1dNmat), nrow = nx, ncol = nt),
	           matrix(as.double(pr2dNmat), nrow = nx, ncol = nt),
	           matrix(as.double(pr1Ymat), nrow = nx, ncol = nt),
	           matrix(as.double(pr2Ymat), nrow = nx, ncol = nt))

  dN12 = t(dNinfo[[8]]) %*% dNinfo[[9]]  #joint dN for groups 1 and 2
  dN1.2 = t(dNinfo[[8]]) %*% dNinfo[[11]]  #conditional dN 1 given at risk in 2
  #maybe transpose the one below?
  dN2.1 = t(dNinfo[[10]]) %*% dNinfo[[9]]  #conditional dN 2 given at risk in 1
  #dbarN=dN1+dN2
  #Y1.small=rep(0,nt)
  #Y1.small=apply(dNinfo[[10]],2,sum)
  #Y2.small=rep(0,nt)
  #Y2.small=apply(dNinfo[[11]],2,sum)
  Y12 = t(dNinfo[[10]]) %*% dNinfo[[11]]  #joint Y for groups 1 and 2

  J = rep(1, nt)                     
  J[Y1*Y2 == 0] = 0

  #correct to here
  Jmat = J %*% t(J)
  #Now I want to calculate G12(u,v)=term1-term2-term3+term4 and G12.pooled(u,v)
  term1 = nstar*dN12 / (Y1 %*% t(Y2))
  term1[Jmat == 0] = 0

  term1.p = (1/nx)*dN12 / ( (km.left.cont*h.left.cont1) %*% t(km.left.cont*h.left.cont2) )
  term1.p[Jmat == 0] = 0

  #temp=dN1.2(u|v)dN2(v)
  temp = t(apply(dN1.2, 1, "*", dN2))
  term2 = nstar*temp / (Y1 %*% t(Y2^2))
  term2[Jmat == 0] = 0

  temp.p = t(apply(dN1.2, 1, "*", dbarN))
  term2.p = (1/nx)*temp.p / ( (km.left.cont*h.left.cont1) %*% t(km.left.cont*h.left.cont2*Ybar) )
  term2.p[Jmat == 0] = 0

  #temp=dN2.1(v|u)dN1(u)
  temp = dN2.1 * dN1
  term3 = nstar*temp / ((Y1^2) %*% t(Y2))
  term3[Jmat == 0] = 0

  temp.p = dN2.1*dbarN
  term3.p = (1/nx)*temp.p / ( (km.left.cont*h.left.cont1*Ybar) %*% t(km.left.cont*h.left.cont2) )
  term3.p[Jmat == 0] = 0

  #temp=dN1(u)dN2(v)
  temp = (dN1 %*% t(dN2))
  term4 = (nstar*Y12*temp) / ((Y1^2) %*% t(Y2^2))
  term4[Jmat == 0] = 0

  temp.p = (dbarN %*% t(dbarN))
  term4.p = (1/nx)*temp.p*Y12 / ( (km.left.cont*h.left.cont1*Ybar) %*% 
	          t(km.left.cont*h.left.cont2*Ybar))
  term4.p[Jmat == 0] = 0



  G12 = term1 - term2 - term3 + term4
  G12.p = term1.p - term2.p - term3.p + term4.p

  #J[Y1.small*Y2.small==0]=0

  weight1 = J*(Y1/n1)*(Y2/n2)*( (n1+n2)/(Y1+Y2) )
  #(1/nx)*(Y1*Y2)/(Y1+Y2)
  weight1[is.na(weight1)] = 0
  weight2 = J*(Y1*Y2) / (n1*n2)

  w1 = weight1 %*% t(weight1)    #create matrices of weights across u,v times
  w2 = weight2 %*% t(weight2)

  #logrank covariance integrand,gehan covariance integrand

  lgrk.cov.integrand = w1*G12
  gehan.cov.integrand = w2*G12

  lgrk.cov.integrand.p = w1*G12.p
  gehan.cov.integrand.p = w2*G12.p

  #covariance piece

  logrank.cov.part = sum(c(lgrk.cov.integrand))
  gehan.cov.part = sum(c(gehan.cov.integrand))

  logrank.cov.part.p = sum(c(lgrk.cov.integrand.p))
  gehan.cov.part.p = sum(c(gehan.cov.integrand.p))

  group1.lgrk.integrand = n1*(weight1^2)*dN1 / (Y1^2)
  group1.lgrk.integrand[dN1 == 0] = 0
  group1.lgrk = sum(group1.lgrk.integrand)

  group1.lgrk.integrand.p = (weight1^2)*dbarN / (km.left.cont*h.left.cont1*Ybar)
  group1.lgrk.integrand.p[is.na(group1.lgrk.integrand.p)] = 0
  group1.lgrk.p = sum(group1.lgrk.integrand.p)

  group2.lgrk.integrand = n2*(weight1^2)*dN2 / (Y2^2)
  group2.lgrk.integrand[dN2 == 0] = 0
  group2.lgrk = sum(group2.lgrk.integrand)

  group2.lgrk.integrand.p = (weight1^2)*dbarN / (km.left.cont*h.left.cont2*Ybar)
  group2.lgrk.integrand.p[is.na(group2.lgrk.integrand.p)] = 0
  group2.lgrk.p = sum(group2.lgrk.integrand.p)

  pi1 = n1/(n1+n2)
  pi2 = n2/(n1+n2)
  theta = (2*n)/(n1+n2)
  lgrk.type.var = pi2*group1.lgrk + pi1*group2.lgrk - theta*logrank.cov.part
  lgrk.type.var.p = pi2*group1.lgrk.p + pi1*group2.lgrk.p - theta*logrank.cov.part.p
  lgrk.nopair.type.var = pi2*group1.lgrk + pi1*group2.lgrk
  lgrk.nopair.type.var.p = pi2*group1.lgrk.p + pi1*group2.lgrk.p

  group1.gehan.integrand = n1*(weight2^2)*dN1 / (Y1^2)
  group1.gehan.integrand[dN1 == 0] = 0
  group1.gehan = sum(group1.gehan.integrand)

  group1.gehan.integrand.p = (weight2^2)*dbarN / (km.left.cont*h.left.cont1*Ybar)
  group1.gehan.integrand.p[is.na(group1.gehan.integrand.p)] = 0
  group1.gehan.p = sum(group1.gehan.integrand.p)

  group2.gehan.integrand = n2*(weight2^2)*dN2 / (Y2^2)
  group2.gehan.integrand[dN2 == 0] = 0
  group2.gehan = sum(group2.gehan.integrand)

  group2.gehan.integrand.p = (weight2^2)*dbarN / (km.left.cont*h.left.cont2*Ybar)
  group2.gehan.integrand.p[is.na(group2.gehan.integrand.p)] = 0
  group2.gehan.p = sum(group2.gehan.integrand.p)

  gehan.type.var = pi2*group1.gehan + pi1*group2.gehan - theta*gehan.cov.part
  gehan.type.var.p = pi2*group1.gehan.p + pi1*group2.gehan.p - theta*gehan.cov.part.p
  gehan.nopair.type.var = pi2*group1.gehan + pi1*group2.gehan
  gehan.nopair.type.var.p = pi2*group1.gehan.p + pi1*group2.gehan.p

  #don't forget the numerator of the test statistic:

  temp = (n1*n2) / (n1+n2)
  num.lgrk = sqrt(temp)*sum(weight1*(haz1 - haz2))
  num.gehan = sqrt(temp)*sum(weight2*(haz1 - haz2))

  #test statistics for logrank and gehan version
  lgrk.type = num.lgrk / sqrt(lgrk.type.var)
  gehan.type = num.gehan / sqrt(gehan.type.var)
  lgrk.type.p = num.lgrk / sqrt(lgrk.type.var.p)
  gehan.type.p = num.gehan / sqrt(gehan.type.var.p)

  lgrk.nopair.type = num.lgrk / sqrt(lgrk.nopair.type.var)
  gehan.nopair.type = num.gehan / sqrt(gehan.nopair.type.var)
  lgrk.nopair.type.p = num.lgrk / sqrt(lgrk.nopair.type.var.p)
  gehan.nopair.type.p = num.gehan / sqrt(gehan.nopair.type.var.p)

  #From here we can create pepe-fleming statistics
  #fairly easily

  if (weights=="left") {
	  if(min(J) == 0) {
		  upperlimhint = min(tm[J == 0])
		  upperlimmark = c(1:nt)[tm == upperlimhint] - 1
		  upperlim = tm[upperlimmark] 
		  } else {
		    upperlim = max(tm) 
		  }

      pfweight1 = rep(1, nt)
      pfweight1[J == 0] = 0
      #for integration purposes you want the weight=0 at the upper limit
      #weights for the pf and yls tests always place weights with integrals
      #and the weight function is visually taking on step-function
      #values until the last time where someone is at risk
      pfweight1[tm == upperlim] = 0
      #Also program PF recommended weight
      pfweight2 = (J*h.left.cont1*h.left.cont2) /                         #left-continuous version
                  (pi1*h.left.cont1 + pi2*h.left.cont2)
      pfweight2[J*h.left.cont1*h.left.cont2 == 0] = 0
      pfweight2[tm == upperlim] = 0
  } else {
	  if(min(J*censorinfo1$kmsurv*censorinfo2$kmsurv) == 0) {
		  upperlimhint = min(tm[J*censorinfo1$kmsurv*censorinfo2$kmsurv == 0])
      upperlimmark = c(1:nt)[tm == upperlimhint] - 1
      upperlim = tm[upperlimmark] 
	  } else { 
      upperlim=max(tm) 
    }
    pfweight1 = rep(1, nt)
    pfweight1[J*censorinfo1$kmsurv*censorinfo2$kmsurv == 0] = 0  #right-continuous version
    pfweight1[tm == upperlim] = 0
    pfweight2 = (J*censorinfo1$kmsurv*censorinfo2$kmsurv) / 
 	              (pi1*censorinfo1$kmsurv + pi2*censorinfo2$kmsurv)      #right-continuous version
    pfweight2[J*censorinfo1$kmsurv*censorinfo2$kmsurv == 0] = 0 
    pfweight2[tm == upperlim] = 0 
  }	

  A1.weight1 = A.1look.revised(km1,tm, maxtau, pfweight1)		
  A2.weight1 = A.1look.revised(km2, tm, maxtau, pfweight1)
  Abar.weight1 = A.1look.revised(km.pooled, tm, maxtau, pfweight1)

  A1.weight2 = A.1look.revised(km1, tm, maxtau, pfweight2)
  A2.weight2 = A.1look.revised(km2, tm, maxtau, pfweight2)
  Abar.weight2 = A.1look.revised(km.pooled, tm, maxtau, pfweight2)

  Amat1 = A1.weight1 %*% t(A2.weight1)    #create matrices of A1(u)A2(v) across u,v times
  Amat2 = A1.weight2 %*% t(A2.weight2)
  Abarmat1 = Abar.weight1 %*% t(Abar.weight1)
  Abarmat2 = Abar.weight2 %*% t(Abar.weight2)

  #yls covariance integrand, pf recommended weight covariance integrand

  yls.cov.integrand = Amat1*G12
  pf.cov.integrand = Amat2*G12
  yls.cov.integrand.p = Abarmat1*G12.p
  pf.cov.integrand.p = Abarmat2*G12.p
 
  #covariance piece
 
  yls.cov.part = sum(c(yls.cov.integrand))
  pf.cov.part = sum(c(pf.cov.integrand))

  yls.cov.part.p = sum(c(yls.cov.integrand.p))
  pf.cov.part.p = sum(c(pf.cov.integrand.p))

  group1.yls.integrand = n1*(A1.weight1^2)*dN1 / (Y1^2)
  group1.yls.integrand[J == 0] = 0
  group1.yls = sum(group1.yls.integrand)
 
  group1.yls.integrand.p = (Abar.weight1^2)*dbarN / 
                           (km.left.cont*h.left.cont1*Ybar)
  group1.yls.integrand.p[J == 0] = 0
  group1.yls.p = sum(group1.yls.integrand.p)
 
  group2.yls.integrand = n2*(A2.weight1^2)*dN2 / (Y2^2)
  group2.yls.integrand[J == 0] = 0
  group2.yls = sum(group2.yls.integrand)
 
  group2.yls.integrand.p = (Abar.weight1^2)*dbarN / 
                           (km.left.cont*h.left.cont2*Ybar)
  group2.yls.integrand.p[is.na(group2.yls.integrand.p)] = 0
  group2.yls.p = sum(group2.yls.integrand.p)
 
  yls.type.var = pi2*group1.yls + pi1*group2.yls - theta*yls.cov.part
  yls.type.var.p = pi2*group1.yls.p + pi1*group2.yls.p - theta*yls.cov.part.p
  yls.nopair.type.var = pi2*group1.yls + pi1*group2.yls
  yls.nopair.type.var.p = pi2*group1.yls.p + pi1*group2.yls.p

  group1.pf.integrand = n1*(A1.weight2^2)*dN1 / (Y1^2)
  group1.pf.integrand[J == 0] = 0
  group1.pf = sum(group1.pf.integrand)
 
  group1.pf.integrand.p = (Abar.weight2^2)*dbarN / 
                          (km.left.cont*h.left.cont1*Ybar)
  group1.pf.integrand.p[J == 0] = 0
  group1.pf.p = sum(group1.pf.integrand.p)

  group2.pf.integrand = n2*(A2.weight2^2)*dN2 / (Y2^2)
  group2.pf.integrand[J == 0] = 0
  group2.pf = sum(group2.pf.integrand)
 
  group2.pf.integrand.p = (Abar.weight2^2)*dbarN / 
                          (km.left.cont*h.left.cont2*Ybar)
  group2.pf.integrand.p[J == 0] = 0
  group2.pf.p = sum(group2.pf.integrand.p)
 
  pf.type.var = pi1*group1.pf + pi2*group2.pf - theta*pf.cov.part
  pf.type.var.p = pi1*group1.pf.p + pi2*group2.pf.p - theta*pf.cov.part.p
  pf.nopair.type.var = pi1*group1.pf + pi2*group2.pf
  pf.nopair.type.var.p = pi1*group1.pf.p + pi2*group2.pf.p

  yls.mean1 = truncmean(km1, tm, maxtau, pfweight1)
  #should be same as A1.weight1[1]  (check this)
  yls.mean2 = truncmean(km2, tm, maxtau, pfweight1)
  #should be same as A2.weight1[1] 
  yls.num = sqrt(temp)*(yls.mean1 - yls.mean2)

  pf.mean1 = truncmean(km1, tm, maxtau, pfweight2)
  pf.mean2 = truncmean(km2, tm, maxtau, pfweight2)

  pf.num = sqrt(temp)*(pf.mean1 - pf.mean2)

  #test statistics for yls and pf version
  yls.type = yls.num / sqrt(yls.type.var)
  pf.type = pf.num / sqrt(pf.type.var)
  yls.type.p = yls.num / sqrt(yls.type.var.p)
  pf.type.p = pf.num / sqrt(pf.type.var.p)

  yls.nopair.type = yls.num / sqrt(yls.nopair.type.var)
  pf.nopair.type = pf.num / sqrt(pf.nopair.type.var)
  yls.nopair.type.p = yls.num / sqrt(yls.nopair.type.var.p)
  pf.nopair.type.p = pf.num / sqrt(pf.nopair.type.var.p)


  answer = 
    list(
      n = nx,
      n1 = n1,
      n2 = n2,
      theta = theta,
      pi1 = pi1,
      pi2 = pi2,
      #dN1=dN1,dN2=dN2,dN12=dN12,dN1.2=dN1.2,dN2.1=dN2.1,
      Y1 = Y1,
      Y2 = Y2,
      #Y12=Y12,G12=G12,
      logrank.cov.part = logrank.cov.part,
      gehan.cov.part = gehan.cov.part,
      group1.lgrk = group1.lgrk,
      group2.lgrk = group2.lgrk,
      lgrk.type.var = lgrk.type.var,
      gehan.cov.part = gehan.cov.part,
      group1.gehan = group1.gehan,
      group2.gehan = group2.gehan,
      gehan.type.var = gehan.type.var,
      #haz1=haz1,haz2=haz2,
      num.lgrk = num.lgrk,
      num.gehan = num.gehan,
      lgrk.stat = lgrk.type,
      gehan.stat = gehan.type,
      km1 = km1,
      km2 = km2,
      pfweight1 = pfweight1,
      pfweight2 = pfweight2,
      upperlim = upperlim,
      yls.cov.part = yls.cov.part,
      pf.cov.part = pf.cov.part,
      group1.yls = group1.yls,
      group2.yls = group2.yls,
      yls.type.var = yls.type.var,
      group1.pf = group1.pf,
      group2.pf = group2.pf,
      pf.type.var = pf.type.var,
      yls.mean1 = yls.mean1,
      yls.mean2 = yls.mean2,
      pf.mean1 = pf.mean1,
      pf.mean2 = pf.mean2,
      yls.stat = yls.type,
      pf.stat = pf.type,
      yls.num = yls.num,
      pf.num = pf.num,
      lgrk.stat.p = lgrk.type.p,
      gehan.stat.p = gehan.type.p,
      lgrk.type.var.p = lgrk.type.var.p,
      gehan.type.var.p = gehan.type.var.p,
      logrank.cov.part.p = logrank.cov.part.p,
      gehan.cov.part.p = gehan.cov.part.p,
      group1.lgrk.p = group1.lgrk.p,
      group2.lgrk.p = group2.lgrk.p,
      group1.gehan.p = group1.gehan.p,
      group2.gehan.p = group2.gehan.p,
      yls.cov.part.p = yls.cov.part.p,
      pf.cov.part.p = pf.cov.part.p,
      group1.yls.p = group1.yls.p,
      group2.yls.p = group2.yls.p,
      yls.type.var.p = yls.type.var.p,
      group1.pf.p = group1.pf.p,
      group2.pf.p = group2.pf.p,
      pf.type.var.p = pf.type.var.p,
      yls.stat.p = yls.type.p,
      pf.stat.p = pf.type.p,
      lgrk.nopair.type.var = lgrk.nopair.type.var,
      lgrk.nopair.type.var.p = lgrk.nopair.type.var.p,
      gehan.nopair.type.var = gehan.nopair.type.var,
      gehan.nopair.type.var.p = gehan.nopair.type.var.p,
      lgrk.nopair.stat = lgrk.nopair.type,
      gehan.nopair.stat = gehan.nopair.type,
      lgrk.nopair.stat.p = lgrk.nopair.type.p,
      gehan.nopair.stat.p = gehan.nopair.type.p,
      yls.nopair.type.var = yls.nopair.type.var,
      yls.nopair.type.var.p = yls.nopair.type.var.p,
      pf.nopair.type.var = pf.nopair.type.var,
      pf.nopair.type.var.p = pf.nopair.type.var.p,
      yls.nopair.stat = yls.nopair.type,
      pf.nopair.stat = pf.nopair.type,
      yls.nopair.stat.p = yls.nopair.type.p,
      pf.nopair.stat.p = pf.nopair.type.p
    )
  
  attr(answer, "class") = "pairtest"
  
  return(answer)
}


#' Print a pairtest object
#'
#' @param object an object of class 'pairtest'

print.pairtest = function(object) {
  
  integration_upper_limit = object$upperlim
  
  logrank_statistic = object$lgrk.stat.p
  logrank_statistic_p = dtail(logrank_statistic)
  
  gehan_statistic = object$gehan.stat.p
  gehan_statistic_p = dtail(gehan_statistic)
  
  yls_statistic = object$yls.stat.p
  yls_statistic_p = dtail(yls_statistic)
  
  pepe_flem_statistic = object$pf.stat.p
  pepe_flem_statistic_p = dtail(pepe_flem_statistic)
  
  logrank_assuming_indep = object$lgrk.nopair.stat.p
  logrank_assuming_indep_p = dtail(logrank_assuming_indep)
  
  gehan_assuming_indep = object$gehan.nopair.stat.p
  gehan_assuming_indep_p = dtail(gehan_assuming_indep)
  
  yls_assuming_indep = object$yls.nopair.stat.p
  yls_assuming_indep_p = dtail(yls_assuming_indep)
  
  pf_assuming_indep = object$ pf.nopair.stat.p
  pf_assuming_indep_p = dtail(pf_assuming_indep)
  
  n1 = object$n1
  n2 = object$n2
  nstuff = ((n1*n2) / (n1 + n2))^(.5)
  
  yls.diff = object$yls.num / nstuff
  
  yls.upper = yls.diff + 1.96*sqrt(object$yls.type.var) / nstuff
  yls.lower = yls.diff - 1.96*sqrt(object$yls.type.var) / nstuff

  yls.upper2 = yls.diff + 1.96*sqrt(object$yls.nopair.type.var) / nstuff
  yls.lower2 = yls.diff - 1.96*sqrt(object$yls.nopair.type.var) / nstuff
  
  cat("\n")
  cat(" ************  Paired sample test results  ************")
  cat("\n")
  cat("\n")
  cat(" Reference Paper: Murray, Susan. Nonparametric Rank-Based Methods for Group")
  cat("\n")
  cat(" Sequential Monitoring of Paired Censored Survival Data. 2000. Biometrics,")
  cat("\n")
  cat(" 56, pp. 984-990.")
  cat("\n")
  cat("\n")
  cat(paste(" Upper limit of integration is", integration_upper_limit))
  cat("\n")
  cat("\n")
  cat("  Logrank statistic: ")
  cat("\n")
  cat(paste("   Logrank estimate = ", logrank_statistic, ", ", "p-value = ", 
            logrank_statistic_p, sep = ""))
  cat("\n")
  cat("\n")
  cat("  Gehan statistic: ")
  cat("\n")
  cat(paste("   Gehan estimate = ", gehan_statistic, ", ", "p-value = ", 
            gehan_statistic_p, sep = ""))
  cat("\n")
  cat("\n")
  cat("  Pepe and Fleming statistic: ")
  cat("\n")
  cat(paste("   PF estimate = ", pepe_flem_statistic, ", ", "p-value = ", 
            pepe_flem_statistic_p, sep = ""))
  cat("\n")
  cat("\n")
  cat("  Years of Life Statistic: ")
  cat("\n")
  cat(paste("   YLS estimate = ", yls_statistic, ", ", "p-value = ", 
            yls_statistic_p, sep = ""))
  cat("\n")
  cat("\n")
  cat(paste("   YLS area between curves = ", yls.diff, "\n   ", "  95% CI for YLS area between curves = (", 
            yls.lower, ", ", yls.upper, ")", sep = ""))
  cat("\n")
  cat("\n")
  cat("\n")
  cat("\n")
  cat(" *****  Results assuming independence (included for comparison)  *****")
  cat("\n")
  cat("\n")
  cat("  Logrank statistic: ")
  cat("\n")
  cat(paste("   Logrank estimate = ", logrank_assuming_indep, ", ", "p-value = ", 
            logrank_assuming_indep_p, sep = ""))
  cat("\n")
  cat("\n")
  cat("  Gehan statistic: ")
  cat("\n")
  cat(paste("   Gehan estimate = ", gehan_assuming_indep, ", ", "p-value = ", 
            gehan_assuming_indep_p, sep = ""))
  cat("\n")
  cat("\n")
  cat("  Pepe and Fleming statistic: ")
  cat("\n")
  cat(paste("   PF estimate = ", pf_assuming_indep, ", ", "p-value = ", 
            pf_assuming_indep_p, sep = ""))
  cat("\n")
  cat("\n")
  cat("  Years of Life Statistic: ")
  cat("\n")
  cat(paste("   YLS estimate = ", yls_assuming_indep, ", ", "p-value = ", 
            yls_assuming_indep_p, sep = ""))
  cat("\n")
  cat("\n")
  cat(paste("   YLS area between curves = ", yls.diff, "\n", "    95% CI for YLS area between curves = (", 
            yls.lower2, ", ", yls.upper2, ")", sep = ""))
  cat("\n")
  cat("\n")
  
}


#' Choose appropriate tail for two-sided test
#' 
#' For internal use only.
#' 
#' @param stat test statistic
#' 
#' @return numeric. p-value for the two-sided test

dtail = function(stat) {
  if(stat <= 0) { return(2*(pnorm(stat, lower.tail = T))) }
  else { return(2*(pnorm(stat, lower.tail = F))) }
}

