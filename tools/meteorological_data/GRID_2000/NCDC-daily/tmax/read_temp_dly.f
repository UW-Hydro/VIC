c     File:           read_temp_dly.f
c     Modified:       09.28.98 by G.M.O.D
c                     Looping rewired to reduce memory overheads.
c                     Code generally cleaned.
c                     Check data and information files are consistent.

c     Programmers:    Greg O'Donnell 1997/2000
c                     Bernt Viggo Matheussen 1998
c                     Univeristy of Washington
c                     Dept of Civil Engineering
c                     Wilcox Hall, Box 352700
c                     Seattle, Washington 98105
c                     tempbvm@ce.washington.edu

      program read_temp_dly

c     This program reads the output from the script preproc_tmax/tmin.scr
c     and formats the daily temperature so the regrid program can read them
c     Only the *.fmt (cmb_prcp.fmt) file is output
      

c     number of stations
      integer nstat

c     ****************************************************************
c     CHANGE NUMBER OF STATIONS HERE
c        number of stations <= nstat
      parameter( in=10, io=11, inf=12, nstat=5000 )
c     ****************************************************************

      integer isaleap
      external isaleap

      character statfile*72      
      character infil*72, outfil*72, info*72
      character*23 name, name_old, stat_name

      real intomm
      parameter( intomm = 25.4 )
      real temp(nstat,366)
      integer id(2)
      logical done

c     leap year always included in input data
      
      print*, ' '
      print*, 'Be sure that the void number used is -99'
      print*, ' '    

      print*, 'Input file timeseries data file: '
      read(*,'(a)') infil
      open(in, file=infil, status='old')

      print*, 'Input file information data file: '
      read(*,'(a)') info
      open(inf, file=info, status='old')
c     read number of stations to be processed and check bounds
      read(inf,*) lines
      if(lines .gt. nstat)then
         print*, 'Increase dimension of nstat to:', lines
         stop
      endif


      print*, 'Output file for temperature time-series: '
      read(*,'(a)') outfil
      open(io, file=outfil, status='unknown')

c      all data will be reference to this date
      print*, 'Start and End years for data:       [INCLUSIVE]'
      read*, istrt_yr, iend_yr

c     on 1st loop prevent accessing array(0) 
      name_old='zzzzzzzzzzzzzzzzzzzz'
      id_old=-999
      nyr=istrt_yr
      iyr_old=iend_yr
      icount=0
      done=.true.
c     5    read(in,'(a22,i5)',end=199) name, iyr
 5    read(in,'(i6,1x,a23,1x,i4)',end=199) id(1), name, iyr
c     station has changed
c      if(name.ne.name_old) then
      if(name.ne.name_old .or. id(1).ne.id_old) then
         icount=icount+1
c     write station list out on 1st loop
         if(nyr.eq.istrt_yr)then
c     check filename in info and timseries files consistent
c            read(inf,'(33x,a22)') stat_name
            read(inf,'(31x,a23,1x,i6)') stat_name, id(2)
c            call name_comp(stat_name,name)
            call name_comp(stat_name,name,id)
         endif
         name_old=name
         id_old=id(1)
c     previous station data not available
         if(.not. done)then
            do i=1, 366
c               print*, icount, nyr
               temp(icount-1,i)=-99
            end do
         endif
         done=.false.
      endif
      
c     required year of data
      if(iyr.eq.nyr)then
         read(in,*) (temp(icount,i), i=1, 366)
         do i=1, 366
c     convert to celcius
            if(temp(icount,i).ne.-99)then
              temp(icount,i)=5*(temp(icount,i)-32)/9.0
	       if(temp(icount,i).gt.50 .or. temp(icount,i).lt.-40)
     $               print*, 'Extreme temp (celcius):  ', temp(icount,i)
            endif
         end do
         done=.true.
      else
c     skip this year
         read(in,*)
      endif
      
      
      iyr_old=iyr

      goto 5

 199  continue


c     missing data for last station
      if(.not.done)then
         do i=1, 366
            temp(icount,i)=-99
         end do
      end if

c     write 1-yr of data for all stations

      LP = isaleap(nyr)
                    
      do i=1, 366
c     negate this is you want - write out the daily values
c     except for feb 29th - unless leap year
         if(LP.eq.0 .and. i.eq.60) goto 99
         write(io,'(9000f7.2)')(max(-99,temp(n,i)), n=1, icount)
 99   end do

      print*, 'Processing', nyr
      nyr=nyr+1
c     check to see if all required years output
      if(nyr.le.iend_yr)then
         iyr_old=iend_yr
         icount=0
         done=.true.
         rewind(in)
         goto 5
      endif

      print*, 'Number of stations processed: ', icount

      stop
      end


****************************************************************
*** SUBROUTINE
****************************************************************

      subroutine name_comp( name1, name2, id )

c     Check that string name2 occurs in name1.
c     name2 must be truncated, of the form ### TEXT_STRING

      character name1*23, name2*23
      integer find, id(2)
      logical id_match

      id_match=.false.

      find=index(name1,name2)
      id_tmp=mod(id(2),10000*(id(2)/10000))
      if(id(1).eq.id_tmp)id_match=.true.
      if(find.ne.1 .or. .not.id_match)then
         print*, 'Information and timeseries filenames do not match'
         print'(i6,1x,a23,1x,i6,1x,a23,1x,i4)',  
     $        id(1), name1, id_tmp, name2, find
         stop
      else
         print*,  name1, ' ', name2
      endif

      return
      end

****************************************************************
*** FUNCTIONS
****************************************************************

      integer function isaleap( iyr )

c     return 1 if a leap yr else 0

      if( (mod(iyr,4) .eq. 0 .and. mod(iyr,100) .ne.0) 
     $                       .or. mod(iyr,400) .eq. 0) then
         isaleap = 1
      else
         isaleap = 0
      endif

      end
