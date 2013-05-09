c     File:           prec_tob_adjust.f

c     Programmers:    Greg O'Donnell 1997/2000
c                     Ed Maurer 2000
c                     Univeristy of Washington
c                     Dept of Civil Engineering
c                     Wilcox Hall, Box 352700
c                     Seattle, Washington 98105

      program prec_tob_adjust

c     This program reads the output from the read_prec_dly.f program
c     and adjusts the daily precipitation to account for the time of
c     observation for each gauge. This prepares a new file in the same
c     format as the .fmt file for input to the regridding program.
      

c     number of stations
      integer nstat

c     ****************************************************************
c     CHANGE NUMBER OF STATIONS HERE :- nstat
c        number of stations <= nstat
      parameter( in=10, io=11, inf=12, it=14, nstat=5000, ntime=40000)
c     ****************************************************************
c     CHANGE UPPER LIMIT ON ALLOWED PRECIPITATION HERE
c     maxprec is the upper plausible limit on allowed precipitation
c     which is in inches/day at this point
      real maxprec
      parameter( maxprec=25 )
c     ****************************************************************
      real void_no
      parameter( void_no=-99 )
c     ****************************************************************

      integer isaleap
      external isaleap
      integer countday
      external countday

      character infil*72, outfil*72, info*72, gtime*72
      character*23 name, name_old, stat_name, off_name(ntime)
      character off_sdate(ntime)*8,off_edate(ntime)*8,off_hrmin(ntime)*4
      integer id(2), off_id(ntime)
      integer obs
      parameter(obs=10)
      real off_frac(nstat,obs)
      real thres
      integer nobs(nstat), date(nstat,obs,2)
      data nobs / nstat*0 /

      real prec(nstat,366), str_p(nstat), stor_frac
      integer ipres(nstat)
      data ipres / nstat*1 /
      logical done
      
      print*, 'Input file:formatted file-output from read_prec_dly: '
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

      print*, 'Input gauge time information data file: '
      read(*,'(a)') gtime
      open(it, file=gtime, status='old')
      call get_gauge_times(it, off_name, 
     $     off_id, off_edate, off_sdate, off_hrmin, ntime, ngt)
      close(it)

      print*, 'Input threshold, mm/d, below which no adjustment: '
      read*, thresh

      print*, 'Output file for formatted precip time-series: '
      read(*,'(a)') outfil
      open(io, file=outfil, status='unknown')

c      all data will be reference to this date
      print*, 'Start year and month for data: '
      read*, istrt_yr, istrt_mo
      print*, 'End year and month for data: '
      read*, iend_yr, iend_mo

      nyr=istrt_yr-1
      nyears=iend_yr-istrt_yr+1

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     begin loop for adjusting data

      do yr=1,nyears
         nyr=nyr+1
c        nday=365+isaleap(nyr)

         if(nyr.eq.istrt_yr) then
            ist_day = countday(istrt_mo,isaleap(nyr))
         else
            ist_day=1
         endif

         if(nyr.eq.iend_yr) then
            nday = countday(iend_mo+1,isaleap(nyr))-1
         else
            nday=365+isaleap(nyr)
         endif
c     read in one year of lines of data file
         do j=ist_day,nday
            read(in,*), (prec(icount,j), icount=1,lines)
         end do

c     check filename in info and timseries and compute offset
         if(nyr.eq.istrt_yr) then
         do icount=1,lines
            read(inf,'(31x,a23,1x,i6)') name, id(2)
            call time_offset(name,id(2),ntime,nstat,obs,ngt,off_frac,
     $           off_name,off_id,off_sdate, off_edate, off_hrmin, 
     $           istrt_yr,nobs,date,icount)
c            print*,'icount,id, nobs(n),off_frac(1): ',icount,id(2),
c     $           nobs(icount),off_frac(icount,1)
         end do
         endif
      
 199  continue
cccccccc
c     for each station, adjust current year
cccccccc
      do n=1, lines
         if(nobs(n).eq.0)then
            continue
c     only one station or no change in offset
         elseif(nobs(n).eq.1 .or. date(n,ipres(n),1).ne.nyr)then
            do i=ist_day, nday-1
c       begin added for threshold
             temp_prec=0.0
             stor_frac=off_frac(n,ipres(n))
             if(prec(n,i).lt.thresh.and.off_frac(n,ipres(n)).le.0.5)then
                if(prec(n,i+1).gt.thresh) then
                   temp_prec=prec(n,i+1)*off_frac(n,ipres(n))
                endif
                off_frac(n,ipres(n))=0.0
             endif
             if(prec(n,i+1).lt.thresh.and.
     $            off_frac(n,ipres(n)).gt.0.5)then
                if(prec(n,i).gt.thresh) then
                   temp_prec=prec(n,i)*(1-off_frac(n,ipres(n)))
                endif
                off_frac(n,ipres(n))=1.0
             endif
c       end added for threshold
               if(prec(n,i).ne.void_no)then
                  if(prec(n,i+1).ne.void_no .or.
     $                 off_frac(n,ipres(n)).eq.0 )then
                     prec(n,i)=prec(n,i)*(1-off_frac(n,ipres(n)))+
     $                   prec(n,i+1)*off_frac(n,ipres(n))+temp_prec
                  else
c                     prec(n,i)=prec(n,i)*(1-off_frac(n,1))
                     prec(n,i)=void_no
                  endif
               endif
c       add line for threshold -- return frac to original value
              off_frac(n,ipres(n))=stor_frac
            end do

c     either more than one station listing or offset changes this year
c     do year in two portions
         else
            if(ist_day.lt. date(n,ipres(n),2)) then
            do i=ist_day, min(date(n,ipres(n),2),nday)
c       begin added for threshold
             temp_prec=0.0
             stor_frac=off_frac(n,ipres(n))
             if(prec(n,i).lt.thresh.and.off_frac(n,ipres(n)).le.0.5)then
                if(prec(n,i+1).gt.thresh) then
                   temp_prec=prec(n,i+1)*off_frac(n,ipres(n))
                endif
                off_frac(n,ipres(n))=0.0
             endif
             if(prec(n,i+1).lt.thresh.and.
     $            off_frac(n,ipres(n)).gt.0.5)then
                if(prec(n,i).gt.thresh) then
                   temp_prec=prec(n,i)*(1-off_frac(n,ipres(n)))
                endif
                off_frac(n,ipres(n))=1.0
             endif
c       end added for threshold

               if(prec(n,i).ne.void_no)then
                  if(prec(n,i+1).ne.void_no .or.
     $                 off_frac(n,ipres(n)).eq.0 )then
                     prec(n,i)=prec(n,i)*(1-off_frac(n,ipres(n)))+
     $                    prec(n,i+1)*off_frac(n,ipres(n))+temp_prec
                  else
c                     prec(n,i)=prec(n,i)*(1-off_frac(n,1))
                     prec(n,i)=void_no
                  endif
               endif
c       add line for threshold -- return frac to original value
              off_frac(n,ipres(n))=stor_frac
            end do
            endif
            ipres(n)=ipres(n)+1

            if(nday.gt. date(n,ipres(n)-1,2)) then
            do i= max(ist_day,date(n,ipres(n)-1,2)+1), nday-1
c           do i=date(n,ipres(n)-1,2)+1, nday-1
c       begin added for threshold
             temp_prec=0.0
             stor_frac=off_frac(n,ipres(n))
             if(prec(n,i).lt.thresh.and.off_frac(n,ipres(n)).le.0.5)then
                if(prec(n,i+1).gt.thresh) then
                   temp_prec=prec(n,i+1)*off_frac(n,ipres(n))
                endif
                off_frac(n,ipres(n))=0.0
             endif
             if(prec(n,i+1).lt.thresh.and.
     $            off_frac(n,ipres(n)).gt.0.5)then
                if(prec(n,i).gt.thresh) then
                   temp_prec=prec(n,i)*(1-off_frac(n,ipres(n)))
                endif
                off_frac(n,ipres(n))=1.0
             endif
c       end added for threshold
               if(prec(n,i).ne.void_no)then
                  if(prec(n,i+1).ne.void_no .or.
     $                 off_frac(n,ipres(n)).eq.0 )then
                     prec(n,i)=prec(n,i)*(1-off_frac(n,ipres(n)))+
     $                    prec(n,i+1)*off_frac(n,ipres(n))+temp_prec
                  else
c                    prec(n,i)=prec(n,i)*(1-off_frac(n,1))
                    prec(n,i)=void_no
                  endif
               endif
c       add line for threshold -- return frac to original value
              off_frac(n,ipres(n))=stor_frac
            end do
         endif
         endif
      end do

c     write dec 31 of previous yr
      if(nyr.gt.istrt_yr)then
         do n=1, lines
c       begin added for threshold
             temp_prec=0.0
             stor_frac=off_frac(n,ipres(n))
             if(str_p(n).lt.thresh.and.off_frac(n,ipres(n)).le.0.5)then
                if(prec(n,1).gt.thresh) then
                   temp_prec=prec(n,1)*off_frac(n,ipres(n))
                endif
                off_frac(n,ipres(n))=0.0
             endif
             if(prec(n,1).lt.thresh.and.
     $            off_frac(n,ipres(n)).gt.0.5)then
                off_frac(n,ipres(n))=1.0
             endif
c       end added for threshold
            if(str_p(n).ne.void_no .and. (prec(n,1).ne.void_no
     $           .or.off_frac(n,ipres(n)).eq.0.0)) then
              str_p(n)=str_p(n)+prec(n,1)*off_frac(n,ipres(n))+temp_prec
            endif
c       add line for threshold -- return frac to original value
            off_frac(n,ipres(n))=stor_frac
         end do
         write(io,'(9000f7.2)'), ((str_p(n)), n=1, lines)
      endif

c     (re)initialize str_p -- store dec 31 for next year
      if(nyr.ge.istrt_yr)then
         do n=1, lines
c       begin added for threshold
            stor_frac=off_frac(n,ipres(n))
            if(prec(n,nday).lt.thresh) then
               if(off_frac(n,ipres(n)).le.0.5)then
                  off_frac(n,ipres(n))=0.0
               else
                  off_frac(n,ipres(n))=1.0
               endif
            endif
c       end added for threshold
            if(nobs(n).eq.0)then
               str_p(n)=prec(n,nday)
            else if(prec(n,nday).ne.void_no)then
               str_p(n)=prec(n,nday)*(1-off_frac(n,ipres(n)))
            else
               str_p(n)=void_no
            endif
c       add line for threshold -- return frac to original value
            off_frac(n,ipres(n))=stor_frac
         end do
      endif

      do i=ist_day, nday-1
            write(io,'(9000f7.2)'), 
     $           (prec(n,i), n=1, lines)
            str_p(n)=0
      end do
      
      print*, 'Processing ', nyr, nday

      end do

c     write final values for dec31
      write(io,'(9000f7.2)'), 
     $     (str_p(n), n=1, lines)
      

      stop
      end

************************************************************************
      subroutine time_offset(name,id,ntime,nstat,obs,ngt,off_frac,
     $     off_name,off_id, off_sdate, off_edate, off_hrmin, istrt_yr, 
     $     nobs,date,istat)
************************************************************************

c     find time offset for the station
c     store as yr, julian day

      integer obs
      character*23 name, off_name(ntime)
      character*8 off_sdate(ntime), off_edate(ntime)
      character*4 off_hrmin(ntime)
      integer off_id(ntime)

      integer nobs(nstat),date(nstat,obs,2)
      real off_frac(nstat,obs)
      
      integer iloc(20)

c     find how many entries in library there are for current station
      icount=0
      do i=1, ngt
         if(id.eq.off_id(i))then
            icount=icount+1
            iloc(icount)=i
         endif
      end do

c     set up offset times
c     #1 no stations available
      if(icount.eq.0)then
         nobs(istat)=0
c     #2 one station
      else if(icount.eq.1)then
         if(off_hrmin(iloc(1)).ne.'    ')then
            nobs(istat)=1
            off_frac(istat,1)=get_frac(off_hrmin(iloc(1)))
            if(off_frac(istat,1).eq.0)nobs(istat)=0
         else
             nobs(istat)=0
          endif
      else
c     #3 multi-stations
         do i=1, icount-1
            read(off_edate(iloc(i)),'(i4,i2,i2)')iyr,imnth,idy
            if(iyr.ge.istrt_yr)then
               if( off_hrmin(iloc(i)).ne.'    ' .and.
     $              off_hrmin(iloc(i)).ne.off_hrmin(iloc(i+1)))then
                  nobs(istat)=nobs(istat)+1
                  off_frac(istat,nobs(istat))=
     $                 get_frac(off_hrmin(iloc(i)))
c     calculate swap over date
                  read(off_edate(iloc(i)),'(i4)') iyr
                  date(istat,nobs(istat),1)=iyr
                  date(istat,nobs(istat),2)=julday(off_edate(iloc(i)))
               endif
            endif
         enddo

c     determine whether last record required
         if(off_hrmin(iloc(icount)).ne.'    ')then
            nobs(istat)=nobs(istat)+1
            off_frac(istat,nobs(istat))=
     $           get_frac(off_hrmin(iloc(icount)))
            date(istat,nobs(istat),1)=9999
         elseif(nobs(istat).gt.1)then
            date(istat,nobs(istat),1)=9999
         endif

c     no stations had times
         if(nobs(istat).eq.0)then
            nobs(istat)=0
         endif

c         print*, off_hrmin(iloc(icount)), iyr, iloc(icount),nobs(istat)
      endif

c      do i=1, nobs(istat)
c         print*, nobs(istat),date(istat,i,1),
c     $        date(istat,i,2),off_frac(istat,i)
c      end do

c      if(nobs(istat).eq.0)print*, 'No time'

      return
      end

****************

      subroutine get_gauge_times(it, off_name, 
     $     off_id, off_edate, off_sdate, off_hrmin, ntime, ngt)

c     process the gauge time file

      character*4 VOID
      parameter(VOID='2400')

      character*4 off_hrmin(ntime), tmp_hrmin
      character*8 off_sdate(ntime), off_edate(ntime)
      character*23 off_name(ntime)
      integer off_id(ntime)

      ngt=1
 15   read(it,901,end=299)  off_id(ngt), off_sdate(ngt), off_edate(ngt),
     $     off_name(ngt), off_hrmin(ngt)
      tmp_hrmin=off_hrmin(ngt)

c     see coop_top.doc for codes
      if(tmp_hrmin.eq.'  SR')then
         off_hrmin(ngt)='0700'
      elseif(tmp_hrmin.eq.'  SS')then
         off_hrmin(ngt)='1900'
      elseif(tmp_hrmin.eq.'0000')then
         off_hrmin(ngt)='2400'
      elseif(tmp_hrmin.eq.' MID')then
         off_hrmin(ngt)='2400'
      elseif(tmp_hrmin.eq.'3030')then
         off_hrmin(ngt)='    '
      elseif(tmp_hrmin.eq.'7777')then
         off_hrmin(ngt)='    '
      elseif(tmp_hrmin.eq.'8888')then
         off_hrmin(ngt)='    '
      elseif(tmp_hrmin.eq.'9999')then
         off_hrmin(ngt)='    '
      elseif(ichar(tmp_hrmin(4:4)).gt.58)then
         off_hrmin(ngt)='    '
      endif

      ngt=ngt+1
      goto 15

 299  ngt=ngt-1

      print*, 'Number of gauge times: ', ngt

      return
 901  format(i6,1x,a8,1x,a8,1x,a23,8x,a4)
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

**********************************************************
      real function get_frac(hrmin)

c     return fraction of precip which belongs to previous day

      character*4 hrmin
      integer hr,min

      read(hrmin,'(2i2)') hr, min
      get_frac=(24.0-(real(hr)+min/60.0))/24.00
      
      return
      end

**********************************************************
      integer function julday( cdate )

      integer iyr, imnth, idy
      character*8 cdate

      integer mnth(12)

      data mnth / 0, 31, 59, 90, 120, 151, 181, 212, 243,
     $     273, 304, 334 /

      read(cdate,'(i4,i2,i2)'), iyr, imnth, idy

      julday=mnth(imnth)+idy

      if(isaleap(iyr).eq.1 .and. imnth.gt.2)julday=julday+1

      return
      end
***********************************************************
c     returns julian day of first day of month passed to function
c     for month 13 it returns 365 or 366

      integer function countday( imonth, ileap )

      integer ileap
      integer imonth

      integer mnth(13)

      data mnth / 0, 31, 59, 90, 120, 151, 181, 212, 243,
     $     273, 304, 334, 365 /

      countday=mnth(imonth)+1
      if(imonth.gt.2) countday=countday+ileap

      return
      end
