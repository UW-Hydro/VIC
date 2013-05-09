      program regrid
c     implicit none

c     This program regrids raw data.  If the temperature sw
c     is set to "t", temperature data is lapsed to the cell
c     elevation before interpolation.  See documentation.

c     Code revised to correct error in calling function latlong 3/9/98
c     in subroutine near_neighbor

c     Code revised to correct error in loading array cellelev() 
c     Line 147 in grd.f  4/17/98

c     Corrected variable type error that caused incorrect read of cell
c     elevations on meter. 10/12/98

c     Corrected code to print no data error messages when using 1 nearneighbor

c     Code revised 10/12/98 to read arc info style mask file as option,
c     write 4-byte binary data,or netCDF format, stop when encountering 
c     NODATA as an option,

c     Note, if MALLOC is used in this code, it is advised that the code
c     is not optimised above level 2.

c     AW-051100:  rewritten to use variable nearest neighbors, depending
c                 on min & max number to use specified by user in argument 
c                 file (regrid.args). 
c                 o  maximum is really just to keep array sizes down -- 
c                    could be set at something like 40 or 50.  
c                 o  minimum should be set around 4, which I believe is
c                    the recommended minimum for the Symap algorithm.
c
c                 I ran this on HP-UX, didn't try on other OSs
c
c     NOTE: the parameter 'printarray' below can be used to 
c           output an array of number of nearest neighbors used,
c           each timestep (row in file) for each cell (column).  
c           the file is large, so the default is no such output.

c     Modified 1/2001 by EdM to remove dynamic memory allocation
c     all array dimensions are defined in array_dims.f
      
c ---------- DECLARATIONS ---------------------------------

      include 'array_dims.f'
      integer ic, is, ir, in, run,ix    !file unit labels
      integer printarray   ! 1=print nearest neighbor array file
                           ! 0=don't print it ("./cellarray.neighbors")
      parameter(printarray=0)

      parameter(run=16, ic=10, is=11, ir=12, io=14,in=15,ix=26)

      integer rc(2)      !rc(1)=rows; rc(2)=cols
      real crnr(2)       !crnr(1)=ylleft; crnr(2)=xlleft
      real cell, void    !resolution; void value
      integer imiss      !arc no-data value (should be real?)

c     counter for output used for writing direct access binary
      integer tcount
      character cjunk*80     !dummy character string
      integer stopflag       !flag to stop program if bad data found
      integer ncell, nstats  !number of computational cells & stations

c     min and max # of nearest neighbors to use in interpolation
      integer min_neigh, max_neigh

      integer i,n,timestep,intwout*2
      real out      !block variable for the output of "symap"

      integer comp(DRC1*DRC2)   !mask variable
      real readmask(DRC1*DRC2)  !temp variable for reading values from mask
      integer neigh(DNCELL,DMAXNEIGH)  !array indices used to rank neighbors
      real wout(DNCELL)         ! output array
      real xcoord(DNSTATS)      !latitudes for input stns
      real ycoord(DNSTATS)      !longitudes for input stns
      real ndist(DNCELL,DMAXNEIGH)    !inverse distances from cell to stns.
      real theta(DNCELL,DMAXNEIGH,DMAXNEIGH)    !weights=(1-cos(c) associated
c                                        !with angle c between 2 nearest neigh
      real timeser(DNSTATS)            !one timesteps worth of station data
      real statelev(DNSTATS)           ! station elevations
      real cellelev(DRC1*DRC2)         ! cell elevations
      integer num_used(DRC1*DRC2)      ! number of stations used per cell
      integer voidcount  !if stopflag=1, use at end of output to tell
                         ! user how many voids are in output

c     Declare command line or run file arguments.
      character tpsw*2
      character infofile*80
      character maskfile*80
      character inputfile*80
      character outputfile*80
      character nneighfile*80
      integer stoprunsw
      integer outformsw
      real mult

c ^^^^^^^^^^ END DECLARATIONS ^^^^^^^^^^^^^^^

c vvvvvvvvvvv  BASIC INPUT vvvvvvvvvvvvvvvvvvv

c     get run information from file
      open(run,file='regrid.runfile',status='old')
      read(run,*) tpsw
      read(run,*) min_neigh, max_neigh      
      read(run,*) infofile
      read(run,*) maskfile
      read(run,*) inputfile
      read(run,*) outputfile
      read(run,*) stoprunsw
      read(run,*) outformsw
      read(run,*) mult
      close(run)
      nneighfile='cellarray.neighbors'
      print*,nneighfile

c check values
      if(max_neigh .le. min_neigh .or. max_neigh .lt. 1 .or. 
     &   min_neigh .lt. 1) then
        print*,'check your max and min nearest neighbor values:'
        print*,'    max = ',max_neigh, '   min = ', min_neigh
        print*,'this combination may give poor results'
        stop
      end if
      if(max_neigh .gt. DMAXNEIGH) then
        print*,'specified maximum neighbors: ',max_neigh
        print*,'exceeds max value defined in array_dim.f: ',DMAXNEIGH
        stop
      end if

      print*, tpsw
      print*,infofile
      print*,maskfile
      print*,inputfile
      print*,outputfile
      print*,stoprunsw
      print*,outformsw
      print*,mult
      
      void=-99      ! define void val 

c     Open mask file
      open(ic,file=maskfile,status='old')

c     read arcinfo style mask file
      read(ic,*) cjunk, rc(2)
      read(ic,*) cjunk, rc(1)
      read(ic,*) cjunk, crnr(2)
      read(ic,*) cjunk, crnr(1)
      read(ic,*) cjunk, cell
      read(ic,*) cjunk, xmiss
      imiss=int(xmiss)

      print*, crnr(1), crnr(2), rc(1), rc(2), cell, imiss
      if(rc(1).eq.0 .or. rc(2).eq.0) then
         print*, 'Error in reading computational grid header.'
         stop
      endif
c     check read values against specified values for array dimensions
      if(rc(1).gt.DRC1 .or. rc(2).gt.DRC2) then
        print*,'rows and cols in mask file: ',rc(1),rc(2)
        print*,'exceeds max values in array_dim.f: ',DRC1,DRC2
        stop
      end if


c     Read the mask and calculate the number of cells
c     and load the cell elevations to an array.
      ncell=0
      read(ic,*) (readmask(i), i=1, rc(1)*rc(2))
      
c     load integer values to comp
      do i=1,rc(1)*rc(2)
         comp(i)=int(readmask(i))
         end do
c     load comp with zero or one and load cellelev with real readmask
      do i=1, rc(1)*rc(2)
         if(comp(i).ne.imiss) then
            ncell=ncell+1
            cellelev(ncell)=readmask(i)
            comp(i)=1
         else
            comp(i)=0
         endif
      end do
      close(ic)
      print*, 'Processing ', ncell, ' active grid cells'

c     Open info file and get max_neigh and nstats
      open(in,file=infofile,status='old')

      read(in,*) nstats

c     check ncell and nstats against specified values for dimensions
      if(ncell.gt.DNCELL .or. nstats.gt.DNSTATS) then
        print*,'numbers of cells or stations: ',ncell,nstats
        print*,'exceeds max values in array_dim.f: ',DNCELL,DNSTATS
        stop
      end if

c     Check to make sure max_neigh not > nstats
      if(max_neigh.gt.nstats)then
      print*, 'Maximum number of interpolation points cannot be'
      print*,'greater than number of stations. Revise info file.'
      stop
      end if

c     Read in lat lon elev data for stations
      do i=1,nstats
        read(in,*) xcoord(i),ycoord(i),statelev(i)
      end do
      close(in)
c     Write out lat lon of input grid
         print*, 'Lat Lon:'
      do i=1,nstats
         print*, xcoord(i),ycoord(i),statelev(i)
      end do

c     Open input data file 
      open(ir,file=inputfile,status='old')      

c     Call subroutine to load theta(x,x,x),ndist(x,x),neigh(x,x)
c     for all computational cells
      print*, 'calculating distances for',ncell,' cells to',nstats,
     &        ' stations...'
      call near_neigh(comp,cell,crnr,ncell,imiss,max_neigh,nstats,
     $     theta,ndist,neigh,rc,xcoord,ycoord)

c  ----------   Final gridding --------------------------------

c     Open output files

      if (printarray .eq. 1) then
        open(ix,file=nneighfile,status='unknown')
      end if

      if(outformsw.eq.0) then
        open(io,file=outputfile)
      else if(outformsw.eq.1)then
        open(io,file=outputfile,access='direct',form='unformatted',
     &       recl=4, status='unknown')
      else if(outformsw.eq.2)then
        open(io,file=outputfile,access='direct',form='unformatted',
     &       recl=2, status='unknown')
      else
        print*,"netCDF option not installed in this version."
        print*,"No output written."
        stop
      end if

c     Read in one timestep from input data file.
c     When end of file is reached, end program.
      timestep=0
      stopflag=0
      tcount=0
      voidcount=0

c %%%%%%%%%%%%%%%%%%%%% TOP OF TIMESERIES LOOP %%%%%%%%%%%%%%%%%%%%%%

 5    read(ir,*,end=199)(timeser(i),i=1,nstats)

        timestep=timestep+1
      
        print*, 'Processing timestep',timestep

c       Interpolate for all cells in comp grid
c       one timestep at a time and write to outfile.
      
        do n=1, ncell     ! active cell loop------
c       Call interpolation subroutine
          call symap(ncell,n,min_neigh,max_neigh,neigh,
     &               ndist,theta,timeser,out,nstats,void,
     &               timestep,tpsw,cellelev,statelev,
     &               stopflag,num_used,mult,outformsw)

          wout(n)=out
          if ( stopflag.eq.1 .and. stoprunsw.eq.1 ) then
            print*, 'Run stopped due to no data at timestep=', 
     &              timestep
            stop
          else if (stopflag .eq.1 .and. stoprunsw .eq. 0) then
            voidcount = voidcount+1
          end if
        end do  ! end of cell loop --------------

c       write out ascii file describing n_neigh used per cell
        if (printarray .eq. 1) then
          write(ix,'(10000I4)'), (num_used(n), n=1, ncell)
        end if

c       write out text of binary according to outformsw
        if (outformsw.eq.0) then
          write(io,'(10000f7.2)'), (wout(n), n=1, ncell)
        else if(outformsw.eq.1) then
          do n=1,ncell
            write(io,rec=(n+ncell*tcount)), wout(n)
          end do 
          tcount=tcount+1
        else if(outformsw.eq.2) then
          do n=1,ncell
            intwout=int(wout(n)*mult)
            write(io,rec=(n+ncell*tcount)), intwout
          end do 
          tcount=tcount+1
        else
          continue
        end if

      goto 5  ! %%%%%%%%%%%%% END OF TIMESTEP LOOP %%%%%%%%%%%%%%%%

c     Bail out instruction for end of input file
 199  continue
      close(ir)
      close(io)
      if (printarray .eq. 1) then
        close(ix)
      end if

      if (voidcount .gt. 0) then
        print*, 'IMPORTANT: regrid wrote ', voidcount,
     &          'voids in the data.'
        print*, 'void value (hardwired) is: ',void
        print*, 'try increasing the max neighbors parameter'
      end if

      print*, 'DONE!'

      stop
      end
      
