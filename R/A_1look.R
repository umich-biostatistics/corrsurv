

#This program calculates the value of Ai1 when it is passed the output
#from SBayes1 and an upperlimit to integrate to


A.1look.long = function(survb, uniqtimes, upperlim, weight) {
  
  #should check to make sure that upperlim is legal
  #before continuing.  Ie the parameter passed to this
  #program should come from min(upperlim,SBayes1$upperlim)
  
  survb=survb[uniqtimes<=upperlim]
  times=uniqtimes[uniqtimes<=upperlim]
  weight=weight[uniqtimes<=upperlim]
  luniqtim=length(uniqtimes)
  #delta(y) vector over ALL time values
  deltay=as.matrix(c((times[2:luniqtim]-times[1:luniqtim-1]),0))
  f.to.integrate=survb[1:luniqtim]*weight*deltay
  f.reverse=rev(f.to.integrate)
  A.reverse=cumsum(f.reverse)
  A=rev(A.reverse)
  
  answer = 
    list(
      A=A,
      deltay=deltay,
      times=times,
      survb=survb,
      f.to.integrate=f.to.integrate,
      f.reverse=f.reverse,
      A.reverse=A.reverse
    )
  
  return(answer)
}



A.1look = function(survb, uniqtimes, upperlim, weight) {
  
  #should check to make sure that upperlim is legal
  #before continuing.  Ie the parameter passed to this
  #program should come from min(upperlim,SBayes1$upperlim)
  #This program might not be correct, might require luniqtim
  #=length(times)  check later
  
  survb=survb[uniqtimes<=upperlim]
  times=uniqtimes[uniqtimes<=upperlim]
  weight=weight[uniqtimes<=upperlim]
  luniqtim=length(uniqtimes)
  #delta(y) vector over ALL time values
  deltay=as.matrix(c((times[2:luniqtim]-times[1:luniqtim-1]),0))
  f.to.integrate=survb[1:luniqtim]*weight*deltay
  f.reverse=rev(f.to.integrate)
  A.reverse=cumsum(f.reverse)
  A=rev(A.reverse)
  
  return(A)
}



A.1look.revised = function(survb, uniqtimes, upperlim, weight) {
  
  #should check to make sure that upperlim is legal
  #before continuing.  Ie the parameter passed to this
  #program should come from min(upperlim,SBayes1$upperlim)
  
  survb=survb[uniqtimes<=upperlim]
  times=uniqtimes[uniqtimes<=upperlim]
  weight=weight[uniqtimes<=upperlim]
  luniqtim=length(times)
  deltay=as.matrix(c((times[2:luniqtim]-times[1:luniqtim-1]),0))
  f.to.integrate=survb*weight*deltay
  f.reverse=rev(f.to.integrate)
  A.reverse=cumsum(f.reverse)
  A=rev(A.reverse)
  
  return(A)
}



truncmean = function(survb, uniqtimes, upperlim, weight) {
  
  #this was debugged on 6/16
  survb=survb[uniqtimes<=upperlim]
  times=uniqtimes[uniqtimes<=upperlim]
  weight=weight[uniqtimes<=upperlim]
  luniqtim=length(uniqtimes)
  deltay=as.matrix(c((times[2:luniqtim]-times[1:luniqtim-1]),0))
  f.to.integrate=survb[1:luniqtim]*weight*deltay
  truncmean=sum(f.to.integrate)
  
  return(truncmean)
}



A.with.Jweight = function(survb, uniqtimes, Jweight) {

  #should check to make sure that Jweight is correct
  #before continuing.  The Jweight, if correct,
  #should be set to zero from the nth term on where
  #n is the upperlimit of integration, also
  #if the length of Jweight is exactly the same as survb,
  #then the last item needs to be zero

  luniqtim=length(uniqtimes)

  #first check dimensions
  if ( Jweight[luniqtim] != 0 ) {
	  print(c("error in A.with.Jweight:1"))
  }
  
    #survb=survb[uniqtimes<=upperlim]
    #times=uniqtimes[uniqtimes<=upperlim]
    #weight=Jweight[uniqtimes<=upperlim]
  
    deltay=as.matrix(c((uniqtimes[2:luniqtim]-uniqtimes[1:luniqtim-1]),0))
    f.to.integrate=survb*Jweight*deltay
    f.reverse=rev(f.to.integrate)
    A.reverse=cumsum(f.reverse)
    A=rev(A.reverse)
    
    return(A)
}

