
getriskdat = function(time.fail.frame, extralooks = c(NA,NA)) {
  
  # If no specific time points are given, calculate estimates at all 
  # unique time points (includes both failure and censoring times) in the data.

  time = time.fail.frame[,1]
  fail = time.fail.frame[,2]
  tmpts = sort(unique(c(0, time, extralooks)))
  nt = length(tmpts)
  nx = length(time)
  outmat = matrix(0, nt, 4)

  temp = surv.info(time, fail, tmpts, outmat, nx, nt)

  list(
    Y = temp$Y,
    Ylag = temp$Ylag,
    dN = temp$dN,
    kmsurv = temp$kmsurv,
    uniq.times = tmpts
  )
}

## Routine which calls fortran routine for doing all
## calculations for at risk info, etc...

surv.info = function(time, fail, tmpts, outmat, nx, nt) {

  dead = as.vector(sort(unique(time[fail == 1])))
  n = length(dead)
  i = sort.list(time)
  x = time[i]

  if (n > 0) {
    storage.mode(outmat) = "double"
    surv = rep(0, nt)
    delta = fail[i]
    Y = as.double((rep(0, n)))
    dN = as.double((rep(0, n)))
    Ylag = as.double((rep(0, n)))
    chaz = as.double((rep(0, n)))

    temp = .Fortran("survinfo",
                    as.integer(n),
                    as.integer(nx),
                    as.integer(nt),
                    as.double(x),
                    as.double(delta),
                    as.double(dead),
                    as.double(tmpts),
                    ans2 = as.double(surv),
                    ans = outmat,
                    Y = as.double(Y),
 	                  dN = as.double(dN),
	                  as.double(Ylag),
	                  chaz = as.double(chaz),
                    PACKAGE = "pairsurv")

    Y = temp$ans[,3]
    dN = temp$ans[,2]
    NA.est = temp$ans[,1]
    NA.var = temp$ans[,4]
    Ylag = Y - dN
    NA.var[Ylag == 0] = NA
    kmsurv = temp$ans2
  } else {

    Y = rep(0, nt)
    temp = .Fortran("yallcens",
    	              as.integer(nx),
    	              as.integer(nt),
                    as.double(x),
                    as.double(tmpts),
                    Y = as.double(Y),
                    PACKAGE = "pairsurv")
    Y = temp$Y
    dN = rep(0, nt)
    NA.est = rep(0, nt)
    NA.var = rep(0, nt)
    Ylag = Y
    kmsurv = rep(1, nt)
  }  #end else

  list(
    Y = Y,
    dN = dN,
    Ylag = Ylag,
    NA.est = NA.est,
    NA.var = NA.var,
    kmsurv = kmsurv
  )
}

