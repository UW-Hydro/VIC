      real function chgtemp(stationelev,cellelev)
      
      implicit none
c     This function takes a station elevation(m) and a
c     cell elevation(m), and calculates a correction to 
c     the temperature obs at the station elevation with
c     reference to the cell elevation.  A change of 0.0065
c     degrees C per meter of elevation is assumed.

      real stationelev,cellelev

      chgtemp=0.0065*(stationelev-cellelev)
      return
      end
