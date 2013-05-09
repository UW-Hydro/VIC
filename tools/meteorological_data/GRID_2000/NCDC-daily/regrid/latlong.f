      real function latlong(lat1,long1,lat2,long2)

*     find distance between pair of lat long coordinates 
*     south and west are negative
      
      implicit none

      real lat1, long1, lat2, long2

      real pi
      parameter( pi=3.141593 )
      real radius
      parameter( radius=6378.0 )
      real dtor
      parameter( dtor=2.0*pi/360.0 )
      
      real theta1, theta2, phi1, phi2, temp
      real term1, term2, term3
      
      theta1 = dtor*long1
      phi1 = dtor*lat1
      theta2 = dtor*long2
      phi2 = dtor*lat2
      term1 = cos(phi1)*cos(theta1)*cos(phi2)*cos(theta2)
      term2 = cos(phi1)*sin(theta1)*cos(phi2)*sin(theta2)
      term3 = sin(phi1)*sin(phi2)
      temp = term1+term2+term3
      temp=min(temp,1.0)
      latlong = radius*acos(temp)
      return
      end
