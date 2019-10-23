      subroutine survinfo(n,nx,nt,x,delta,d,tm,sv,outm,Y,dN,Ylag,chaz)

      integer n,nt,i,nx,j,k,lastY
      double precision d(n),x(nx),delta(nx),tm(nt),outm(nt,4),sv(nt)
      double precision Y(n),dN(n),old,temp1,temp2,Ylag(n),tempsv
      double precision chaz(n),oldch 

c     this fortran program calculates the Nelson-Aalen estimator,
c     the Kaplan-Meier survival estimator, and reports the number
c     at risk Y(t) as well as dN(t)

c     x: sorted observed failure times
c     delta: corresponding failure indicator
c     d: vector of sorted observed death times
c     n: number of observed death times (length of d)
c     nx: number of event times (length of x)
c     nt: length of unique tmpts of interest
c     tm: vector of sorted tmpts of interest
c     sv: survival curve
c     outm: outmatrix

      old=0.d0
      oldch=0.d0
      lastY=nx 
 
c     Create dN(t), Y(t) at each unique death time
c     chaz is value of nelson estimator at each death time
 
      do 5 i=1,n
c     (ie, over unique death times)
          temp1=0.d0
          temp2=0.d0
          Y(i)=0.d0
          dN(i)=0.d0
 
          do 7 j=1,nx
c         (ie, over number of event times)
             if(x(j).ge.d(i)) temp1=temp1+1.d0
c            (above counts the number at risk Y(t)
             if(x(j).le.d(i)) temp2=temp2+delta(j)
c            (above creates the counting process N(t)
 7       continue

c     small sample adjustment - if no one at risk, set Y=1.

         Y(i) = max(temp1,1.d0)

         dN(i)=temp2-old
         old=old+dN(i)
         Ylag(i)=max((Y(i)-dN(i)),1.d0)
c        (Ylag is number at risk - number of deaths)
         chaz(i)=oldch+dN(i)/Y(i)
         oldch=chaz(i)
 
 5       continue
 
c     set KM(0)=1
 
      sv(1)=1.d0

c     set Y(0)= total number of event times 
c     calculate KM, NA, greenwood's formula for nelson hazard estimator, and
c     collect  Y(t) and dN(t) for all time points requested
 
      outm(1,3)=nx
      outm(1,4)=0.0d0
      
      do 10 i=1,nt
c     (ie, over all unique times requested) 
         tempsv=1.d0
         outm(i,1)=0.d0
         outm(i,2)=0.d0

 
         do 20 j=1,n
c        (ie, over unique death times)

            if (d(j).le.tm(i)) then
 
c     calculate product limit estimator
 
               tempsv=tempsv*(1.d0-dN(j)/Y(j))

c     calculate NA(x) 
 
               outm(i,1)=outm(i,1)+dN(j)/Y(j)

c     computes the greenwood version of the variance of NA
 
               outm(i,4)=outm(i,4)+dN(j)/(Y(j)*Ylag(j))
 
            endif

c     spreads dN(t) over all the unique times requested
 
            if (d(j).eq.tm(i)) outm(i,2)=dN(j)

 20         continue
 
c     computes KM estimate
 
            if (i.gt.1) sv(i)=tempsv

c     computes number at risk at all times requested

	 if (i.lt.nt) then 
c       (ie, we don't want to have an error for i+1 = nt +1)
  
	 outm(i+1,3)=lastY

         do 30 k=1,nx-1
c       (ie, over all but the last event time)

              if (x(k).eq.tm(i)) then

                 outm(i+1,3)=max((lastY-1.0d0),0.0d0)  
                 lastY=outm(i+1,3)   

              endif

 30         continue
         
          endif 
          if (tm(i).gt.x(nx)) outm(i,3) = 0.0d0

 10         continue

          if (tm(nt).gt.x(nx)) outm(ny,3) = 0.0d0
 
            return
            end
 



c      program test
c
c      integer  nt, nx 
c      double precision x, delta, d, Y, dN, Ylag, chaz
c      double precision sv(5), outm(5,4), tm(5)
c     
c      data (tm(i),i=1,5)/0.0,0.02982376,0.03871822,5.417064,8.0/
c
cc     initialize, read in data
c
c      n=1
c      nt=5
c      nx=1
c      x=5.417064d0
c      delta=1.0d0
c      d=5.417064d0
c      Y=0.0d0
c      dN=0.0d0
c      Ylag=0.0d0
c      chaz=0.0d0
c      do 40 i=1,nt
c         sv(i)=0.0d0
c         outm(i,1)=0.0d0
c         outm(i,2)=0.0d0 
c         outm(i,3)=0.0d0 
c         outm(i,4)=0.0d0
c 40    continue 
cc       tm=(/ 0.00000000d0, 0.02982376d0, 0.03871822d0, 5.417064d0, 8.0d0 /)
c     
c      call survinfo(n,nx,nt,x,delta,d,tm,sv,outm,Y,dN,Ylag,chaz)
c      print *, (outm(i,3), i=1,5) 
c
c      end 
