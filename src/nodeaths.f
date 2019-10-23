      subroutine yallcens(nx,nt,x,tm,Y)
 
      integer nt,i,nx,k,lastY
      double precision x(nx),tm(nt), Y(nt)
 
c     this fortran program reports the number
c     at risk Y(t) when only censoring is present 
c     in the data  (ie, no observed failures)
 
c     x: sorted observed failure times
c     nx: number of event times (length of x)
c     nt: length of unique tmpts of interest
c     tm: vector of sorted tmpts of interest

      Y(1)=nx 
      lastY=nx

      do 5 i=1,nt-1
c     (ie, over all unique times requested)
 
      Y(i+1)=lastY

         do 10 k=1,nx-1
c       (ie, over all but the last event time)
 
              if (x(k).eq.tm(i)) then
 
                 Y(i+1)=max((lastY-1.0d0),0.0d0)
                 lastY=Y(i+1)
 
              endif
 
 10       continue
 
      if (tm(i).gt.x(nx)) Y(i) = 0.0d0
 
 5    continue

      if (tm(nt).gt.x(nx)) Y(nt) = 0.0d0 

      return
      end
 
 
c        program test
c
c      integer, parameter::  nt=350, nx=21 
c      double precision x(nx), tm(nt), Y(nt)
c	
c  initialize, read in data
c
c        open(11,file='testdat',status='old')
c        read(11, *) (x(i), i=1,nx)
c	read(11, *) (tm(i), i=1,nt)

	

