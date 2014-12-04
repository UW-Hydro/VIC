module ndepStreamMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: ndepStreamMod
! 
! !DESCRIPTION: 
! Contains methods for reading in nitrogen deposition data file
! Also includes functions for dynamic ndep file handling and 
! interpolation.
!
! !USES
  use shr_kind_mod, only: r8 => shr_kind_r8, CL => shr_kind_cl

! !PUBLIC TYPES:
  implicit none 
  save

! !PUBLIC MEMBER FUNCTIONS:
  public :: ndep_init   ! position datasets for dynamic ndep
  public :: ndep_interp ! interpolates between two years of ndep file data
!
!EOP

! ! PRIVATE TYPES

  type, private :: shr_strdata_type
    real(r8), pointer :: syr(:)           ! file year
    real(r8), pointer :: slat(:)          ! file lat
    real(r8), pointer :: slon(:)          ! file lon
    real(r8), pointer :: ndep(:,:,:)      ! N deposition
  end type shr_strdata_type

  type(shr_strdata_type), private  :: sdat         ! input data stream

   integer, private   :: nyrs      ! number of years in file
   integer, private   :: nlat      ! number of latitudes in file
   integer, private   :: nlon      ! number of longitudes in file

!=======================================================================
contains
!=======================================================================

  subroutine ndep_init(begg, endg)

   !----------------------------------------------------------------------- 
   ! Initialize data stream information.  
   !----------------------------------------------------------------------- 
   ! Uses:
   use clmtype, only : clm_a2l
   ! arguments
   implicit none
   integer :: begg, endg
   ! local variables
   integer            :: t, j, i, n   ! indices
   !-----------------------------------------------------------------------

   open(10, file='fndep_clm_hist_simyr1849-2006_1.9x2.5_c100428')
   read(10, '(I4)') nyrs
   read(10, '(I4)') nlat
   read(10, '(I4)') nlon

   allocate(sdat%syr(1:nyrs))
   allocate(sdat%slat(1:nlat))
   allocate(sdat%slon(1:nlon))
   allocate(sdat%ndep(1:nyrs,1:nlat,1:nlon))

   t = 1
   j = 1
   i = 1
   do n = 1, nyrs * nlat * nlon
     read(10, '(I4,F18.13,F7.3,F9.6)') sdat%syr(t), sdat%slat(j), &
	sdat%slon(i), sdat%ndep(t,j,i)
     i = i + 1
     if(i > nlon) then
       j = j + 1
       i = 1
     endif
     if(j > nlat) then
       t = t + 1
       j = 1
       i = 1
     endif
   enddo
   close(10)

   allocate(clm_a2l%forc_ndep(begg:endg))

 end subroutine ndep_init
  
!================================================================

 subroutine ndep_interp( year, jday, begg, endg )

   !-----------------------------------------------------------------------
   use clmtype         , only : clm3, clm_a2l
   use clm_time_manager, only : get_days_per_year
   use clm_varcon      , only : secspday

   implicit none
   ! arguments
   integer :: year    ! year
   integer :: jday    ! day in year
   integer :: begg, endg ! bounds
   ! Local variables
   integer :: g, igx, igy, igt, j
   integer :: dayspyr ! days per year
   integer :: midyear ! mid-year day
   !-----------------------------------------------------------------------

   igt = year - sdat%syr(1) + 1

   dayspyr = get_days_per_year(year)
   do g = begg,endg
     igy = 0
     if(clm3%g%lat(g) < sdat%slat(1) + 0.5 * (sdat%slat(1) + sdat%slat(2))) &
	then
       igy = 1
     else
       do j = 2, nlat - 1
         if(clm3%g%lat(g) > sdat%slat(j) - 0.5 * (sdat%slat(j - 1) + &
		sdat%slat(j)) .and. clm3%g%lat(g) < sdat%slat(j) + 0.5 * &
		(sdat%slat(j) + sdat%slat(j + 1))) then
           igy = j
         endif
       enddo
       if(igy == 0) then
         igy = nlat
       endif
     endif

     igx = 0
     if(clm3%g%lon(g) + 360. < sdat%slon(1) + 0.5 * (sdat%slon(1) + &
	sdat%slon(2)) .or. clm3%g%lon(g) + 360. > sdat%slon(nlon) + 0.5 * &
	(sdat%slon(nlon - 1) + sdat%slon(nlon))) then
       igx = 1
     else
       do j = 2, nlon - 1
         if(clm3%g%lon(g) + 360. > sdat%slon(j) - 0.5 * (sdat%slon(j - 1) + &
		sdat%slon(j)) .and. clm3%g%lon(g) + 360. < sdat%slon(j) + 0.5 &
		* (sdat%slon(j) + sdat%slon(j + 1))) then
           igx = j
         endif
       enddo
       if(igx == 0) then
         igx = nlon
       endif
     endif

      if(mod(year, 4) == 0) then
        midyear = 167
      else
        midyear = 166
      endif
      if(jday <= midyear) then
        if(year > sdat%syr(1)) then
          clm_a2l%forc_ndep(g) = sdat%ndep(igt-1,igy,igx) + &
		(sdat%ndep(igt-1,igy,igx) + sdat%ndep(igt,igy,igx)) / &
		(secspday * dayspyr)
        else
          clm_a2l%forc_ndep(g) = sdat%ndep(igt,igy,igx)
        endif
      else
        if(mod(year, 4) == 0) then
          dayspyr = dayspyr - 1
        endif
        if(mod(year + 1, 4) == 0) then
          dayspyr = dayspyr + 1
        endif
        if(year < sdat%syr(nyrs)) then
          clm_a2l%forc_ndep(g) = sdat%ndep(igt,igy,igx) + &
		(sdat%ndep(igt,igy,igx) + sdat%ndep(igt + 1,igy,igx)) / &
		(secspday * dayspyr)
        else
          clm_a2l%forc_ndep(g) = sdat%ndep(igt,igy,igx)
        endif
      endif
   end do
   
 end subroutine ndep_interp

end module ndepStreamMod

