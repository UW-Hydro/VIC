module clm_time_manager

   use shr_kind_mod, only: r8 => shr_kind_r8
   use clm_varcon  , only: isecspday

   implicit none
   private

! Public methods

   public :: get_curr_calday, &! return calendar day at end of current timestep
      get_days_per_year       ! return the days per year for current year

! Private module data

contains

real function get_curr_calday(day_in_year, hour, offset)

! Return calendar day at end of current timestep with optional offset.
! Calendar day 1.0 = 0Z on Jan 1.

   implicit none
! Arguments
   integer, intent(in)           :: day_in_year ! Julian day in year
   integer, intent(in)           :: hour    ! Beginning of current hour
   integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
                                            ! Positive for future times, negative 
                                            ! for previous times.
! Return value

!-----------------------------------------------------------------------------------------

   get_curr_calday = real(day_in_year) + (real(hour) / 24.)

   if (present(offset)) then
     get_curr_calday = get_curr_calday + real(offset / 86400)
   else
     get_curr_calday = get_curr_calday + (1. / 24.)
   end if

   return
end function get_curr_calday

!=========================================================================================

real function get_days_per_year( year )

  ! Get the number of days per year for currrent year
  implicit none
  !
  ! Arguments
  integer, intent(in)           :: year    ! Current year

  if(mod(year, 4) == 0) then
    get_days_per_year = 366._r8
  else
    get_days_per_year = 365._r8
  end if
  
  return
end function get_days_per_year

end module clm_time_manager
