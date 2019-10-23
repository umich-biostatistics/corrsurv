      subroutine bivsurvmats (nx,nt,tm,x1,delta1,x2,delta2,
     1     pr1dNmat,pr2dNmat,pr1Ymat,pr2Ymat)

      integer nx,nt,i,j
      double precision x1(nx),delta1(nx),x2(nx),delta2(nx),tm(nt)
      double precision pr1dNmat(nx,nt),pr2dNmat(nx,nt),
     1     pr1Ymat(nx,nt),pr2Ymat(nx,nt)
c     nx: number of event times (common to both groups, length of x1)
c     nt: length of unique tmpts of interest
c     tm: vector of sorted tmpts of interest
c     x1: observed event times for pair member 1, sorted by pair id
c     delta1: corresponding failure indicator for event time 1
c     x2: observed event times for pair member 2, sorted by pair id
c     delta2: corresponding failure indicator for event time 2
c     pr1dNmat: nx by nt matrix where row are across individuals
c               and columns are over unique tmpts of interest
c	        collecting I(X1=u,delta=1)  (for 1st member of pairs)
c     pr2dNmat: nx by nt matrix where row are across individuals
c               and columns are over unique tmpts of interest
c               collecting I(X2=v,delta=1)  (for 2nd member of pairs)
c     pr1Ymat: nx by nt matrix where row are across individuals
c               and columns are over unique tmpts of interest
c               collecting I(X1>=u) (for 1st member of pairs)
c     pr2Ymat: nx by nt matrix where row are across individuals
c               and columns are over unique tmpts of interest
c               collecting I(X2>=u) (for 2nd member of pairs)


      do 7 i=1,nx
c     (ie, over unique individuals)

           do 9 j=1,nt
c              (ie, over all unique timepoints of interest)


c	          fill in values of pr1dNmat

                  pr1dNmat(i,j)=0.0d0 
        	  if( (x1(i).eq.tm(j))  .and. (delta1(i).eq.1) ) 
     1        		pr1dNmat(i,j)=1.0d0

c                 fill in values of pr2dNmat

	          pr2dNmat(i,j)=0.0d0 
                  if( (x2(i).eq.tm(j))  .and. (delta2(i).eq.1) ) 
     1		        pr2dNmat(i,j)=1.0d0

c	          fill in values of pr1Ymat

		  pr1Ymat(i,j)=0.0d0
		  if( x1(i).ge.tm(j) ) pr1Ymat(i,j)=1.0d0

c                 fill in values of pr2Ymat

		  pr2Ymat(i,j)=0.0d0
                  if( x2(i).ge.tm(j) ) pr2Ymat(i,j)=1.0d0

 9             continue

 7        continue




            return
            end

