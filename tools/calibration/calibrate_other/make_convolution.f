
      SUBROUTINE MAKE_CONVOLUTION
     & (NCOL, NROW, NOB, PMAX, DAYS, CATCHIJ, 
     &  BASE, RUNO, FLOW, KE, UH_DAY, UH_S, FRACTION, FACTOR_SUM,
     $     XC, YC, SIZE, DPREC, INPATH,ICOL,NDAY)

      IMPLICIT NONE

      INTEGER     N, I, J, DAYS, NDAY, II, JJ  
      INTEGER     NCOL,NROW,ICOL,NOB,PMAX,KE,UH_DAY
      INTEGER     CATCHIJ(PMAX,2)
      REAL        UH_S(PMAX,KE+UH_DAY-1)
      REAL        BASE(DAYS), RUNO(DAYS), FLOW(DAYS) 
      REAL        FRACTION(NCOL,NROW)

      REAL        PI, RERD, FACTOR, FACTOR_SUM

      PARAMETER   (RERD  = 6371229.0)

      CHARACTER*20 LOC
      CHARACTER*72 INPATH

      INTEGER DPREC, CLEN

      REAL        JLOC, ILOC
      REAL        XC, YC, SIZE
      REAL        AREA, AREA_SUM

      REAL        STORAGE, K_CONST
      REAL        DUM1,DUM2
  
      INTEGER     IDAY, IMONTH, IYEAR

C     *** 0 <= K_CONST = 1.0 
C *** K_CONST smaller 1.0 makes it a simple linear storage
    
      K_CONST = 1.0

      PI = ATAN(1.0) * 4.0

      AREA_SUM   = 0.0
      FACTOR_SUM = 0.0

      DO I = 1,NDAY
         FLOW(I) = 0.0
      END DO

      DO N = 1,NOB
         STORAGE = 0.0
         DO I = 1,NDAY
            RUNO(I) = 0.0
            BASE(I) = 0.0
         END DO
         II = CATCHIJ(N,1)
         JJ = CATCHIJ(N,2)
         
c     the grid has been flipped left to right
c     find the revised cooordinates

         ILOC=XC + (ICOL-II)*SIZE + SIZE/2.0
         JLOC=YC + JJ*SIZE - SIZE/2.0

C        CONVERSIONFACTOR for mm/day to ft**3/sec


         AREA =  RERD**2*ABS(SIZE)*PI/180*
     &        ABS(SIN((JLOC-SIZE/2.0)*PI/180)-
     $        SIN((JLOC+SIZE/2.0)*PI/180))

         
         AREA_SUM = AREA_SUM + AREA

c        WRITE(*,*) N, ILOC, JLOC

         FACTOR = FRACTION(II,JJ)*35.315*AREA/(86400.0*1000.0)

         FACTOR_SUM = FACTOR_SUM + FACTOR
         
         call create_vic_names(jloc,iloc,loc,clen,dprec)

c     print*, INPATH(1:INDEX(INPATH,' ')-1)//LOC(1:CLEN)
         OPEN(20,FILE=INPATH(1:(INDEX(INPATH,' ')-1))//
     $        LOC(1:CLEN),
     $        STATUS='OLD',ERR=9001)

         DO I = 1,NDAY
            READ(20,*,END=9001,ERR=9001) IYEAR, IMONTH, IDAY, 
     $        RUNO(I), BASE(I)
         END DO

         DO I = 1,NDAY
            RUNO(I) = RUNO(I) * FACTOR
            BASE(I) = BASE(I) * FACTOR
         END DO
         DO I = 1,NDAY
            DO J = 1,KE+UH_DAY-1
               IF ((I-J+1) .GE. 1) THEN
                  FLOW(I) = FLOW(I)+UH_S(N,J)*(BASE(I-J+1)+RUNO(I-J+1))
               END IF
            END DO
         END DO
         CLOSE(20)
      END DO
      RETURN
 9001 WRITE(*,*) 'Error reading time-series data, ',
     $     'insufficient data or missing input file',
     $     INPATH(1:INDEX(INPATH,' ')-1)//LOC(1:CLEN)
      END

