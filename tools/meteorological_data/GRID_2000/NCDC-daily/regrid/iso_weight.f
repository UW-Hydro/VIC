      real function iso_weight( a, b, c )

*     calculate the isolation weight
*     returns cos A, where A is the angle about the prism centre
*     linking two rain gages

      implicit none

      real pi
      parameter( pi=3.141593 )

      real radius
      parameter( radius=6378.0 )
      
      real circum
      parameter( circum=2.0*pi*radius )

      real dtor
      parameter( dtor=pi/180.0 )

      real a,b,c

*     first convert lengths into degrees - note lengths are great arcs
*     and convert to radians
      a=(a/circum)*360.0*dtor
      b=(b/circum)*360.0*dtor
      c=(c/circum)*360.0*dtor
      iso_weight=(cos(a)-(cos(b)*cos(c))) / MAX(1.0E-04,(sin(b)*sin(c)))
*     slight adjustments required for floating-point errors
      if(abs(iso_weight).gt.1.0)then
         iso_weight=nint(iso_weight)
      end if
      return
      end
