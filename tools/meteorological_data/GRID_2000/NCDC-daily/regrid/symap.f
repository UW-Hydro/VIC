c       Symap routine
c
c       PURPOSE: perform interpolation from station data to cell centers
c	INPUT:   total number of cells
c                current cellnumber
c                number of nearest neighbors
c                inverse distance from each nearest neighbors to cell center 
c                (theta):  1 minus cos(c) for nearest neighbors
c                one timestep's worth of input data (all stations)
c
c       OUTPUT:  a single interpolated value for the cellnumber.
c 
c       Note:  If the number of nearest neighbors is 1, then the
c       interpolation routine (Symap algorithm) is bypassed.
c
c       AW-051100:  will use up to the 'max_neigh' nearest neighbors
c                   in order to find 'min_neigh' valid data values
c 

      subroutine symap(numcells,cellnum,mn_neigh,mx_neigh,
     &                   nnarrpointer,invdist,theta,inputdataarray,
     &                   outputdatapt,numstations,void,timestep,
     &                   tpsw,cellelev,statelev,stopflag,numb_used,
     &                   mult,outformsw)

      include 'array_dims.f'

c      implicit none
      real chgtemp
      external chgtemp
      
c     Declare shared variables
      character tpsw*2
      integer numcells,cellnum,mx_neigh,timestep,numstations,mn_neigh
      integer nnarrpointer(DNCELL,DMAXNEIGH)
      integer stopflag, outformsw
      real invdist(DNCELL,DMAXNEIGH)
      real theta(DNCELL, DMAXNEIGH, DMAXNEIGH)
      real inputdataarray(DNSTATS)
      real outputdatapt
      real mult, tempmult
      real cellelev(DNCELL),statelev(DNSTATS)
      real void
      integer numb_used(DNCELL)

c     Declare local variables   
      integer i,j,k
      real H,sumweight
      integer reduced_nn  !number, per cell & timestep, of neighbors needed
      integer nncount     !number of valid cells used (should equal mn_neigh)

c     -- following declarations work for HP-UX compiler; otherwise
c        use dynamic arrays for these three
      integer exclude(DMAXNEIGH)
      real T(DMAXNEIGH)
      real weight(DMAXNEIGH)

c     -- START CALCULATIONS ----
      do i=1,mx_neigh
	T(i)=0
      end do	

c     first find how far into array you need to go to get enough values
      i=0
      nncount=0

      do while (nncount .lt. mn_neigh .and. i .lt. mx_neigh)
        i = i+1
        exclude(i)=0    ! void value default
  	if(inputdataarray(nnarrpointer(cellnum,i)).ne.void) then
	  exclude(i)=1  ! means good value
          nncount = nncount+1     !number of actual good values
        end if
      end do
      numb_used(cellnum) = i   !this number of cells is needed
      reduced_nn = numb_used(cellnum)  !to truncate the calculations below

      if(reduced_nn .lt. 1) then  !GET OUT, something very wrong
        print *, 'reduced_nn=',reduced_nn, ': BIG PROBLEM'
        stop
      end if 

      do i=1, reduced_nn
	if(exclude(i).eq.0) then  ! skip over T calc
	  go to 99
	end if
	do j=1, reduced_nn !i.e., nearest neighbor used count
          if((inputdataarray(nnarrpointer(cellnum,j)).ne.void).and.
     $      (i.ne.j)) then
	    T(i)=T(i) + invdist(cellnum,j)*theta(cellnum,i,j)
	  end if
	end do	
 99	continue
      end do

c     Calculate H
      H=0
      do k=1,reduced_nn
        H=H+invdist(cellnum, k)*exclude(k)
      end do

c     If H=0 then there are no good data points.  Print message, exit

      if(H.eq.0) then        
        write(*,*) 'There are no near neighbors with valid data',
     &             ' for comp cell',cellnum,'at timestep',timestep
	write(*,*) 'Void value written if no stop specified'
	outputdatapt=void
        stopflag=1
	return      ! skip the datapoint calculation
      end if	

c     Calculate final weights and the sum of weights
      sumweight=0
      do k=1,reduced_nn
        weight(k)=(invdist(cellnum,k)**2)*(H+T(k))*exclude(k)
        sumweight=sumweight + weight(k)	
      end do

c     Calculate interpolated data point
      outputdatapt=0
      if(tpsw.eq.'p')then
	do k=1,reduced_nn
	  outputdatapt=outputdatapt + weight(k)* 
     &    inputdataarray(nnarrpointer(cellnum,k))
	end do
      else
	do k=1,reduced_nn
	  outputdatapt=outputdatapt + weight(k)* 
     $    (inputdataarray(nnarrpointer(cellnum,k))+
     $    chgtemp(statelev(nnarrpointer(cellnum,k)),cellelev(cellnum)))
	end do
      end if

      outputdatapt=outputdatapt/sumweight

c     check for crazy values

      if(tpsw.eq.'p' .and. 
     &  (outputdatapt .gt. 350 .or. outputdatapt .lt. -350)) then
        print*, 'possible bad value at cell ',cellnum,
     &          ' timestep ', timestep
        stopflag = 1
        print*, 'cellnum', cellnum, ' data=',outputdatapt
        print*, 'some debugging data:'
        print*, 'sumweight=',sumweight
        do i=1,k-1
          print*,'pt=',i,' datapt=',outputdatapt,' weight=',weight(i),
     &       ' inputdata=',inputdataarray(nnarrpointer(cellnum,i)),
     &       ' normalized weight=',weight(i)/sumweight
        end do
        print*, 'tpsw=',tpsw
      else if (tpsw.eq.'t' .and.
     &  (outputdatapt.gt.100.or.outputdatapt.lt.-100)) then
        print*, 'possible bad value at cell ',cellnum,
     &          ' timestep ', timestep
        stopflag = 1
        print*, 'cellnum', cellnum, ' data=',outputdatapt
        print*, 'some debugging data:'
        print*, 'sumweight=',sumweight
        do i=1,k-1
          print*,'pt=',i,' datapt=',outputdatapt,' weight=',weight(i),
     &       ' inputdata=',inputdataarray(nnarrpointer(cellnum,i)),
     &       ' chgtemp=',
     &     chgtemp(statelev(nnarrpointer(cellnum,i)),cellelev(cellnum)),
     &     ' stn_elev=',statelev(nnarrpointer(cellnum,i)), 
     &     ' cellelev=',cellelev(cellnum), 
     &     ' normalized weight=',weight(i)/sumweight
        end do
        print*, 'tpsw=',tpsw
      end if
c check to be sure binary data is small enough for multiplier
      if (outformsw.eq.2) then
         if ((outputdatapt*mult).ge.32767 .or. 
     &        (outputdatapt*mult).le.-32767) then
            tempmult=32767.0/outputdatapt
            print*, 'multiplier too high at timestep ', timestep
            print*, 'cellnum ', cellnum, ' data= ',outputdatapt
            print*, 'reduce multiplier below ',tempmult
            print*, 'if value is not an error '
            stopflag = 1
         end if
      end if
      return
      end




