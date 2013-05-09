*************************************
      integer function isaleap(year)
*************************************
*     determine whether present year is a leap year
 
      implicit none
 
      integer year
 
      if( (mod(year,4).eq.0 .and. mod(year,100).ne.0)
     $     .or. mod(year,400).eq.0) then
         isaleap=1
      else
         isaleap=0
      endif
 
      end
