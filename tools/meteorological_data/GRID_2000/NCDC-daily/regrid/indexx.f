
      subroutine indexx(nn,n,arrdis,arrindx)

      include 'array_dims.f'

      INTEGER nn,n,arrindx(DNSTATS)
      REAL arrdis(DNSTATS)
      INTEGER i,j,temp
      

c  Initialize array arrindx(n)

      do 11 j=1,n
        arrindx(j)=j
11    continue
  
c  Float data cell numbers with lowest distances 
c  to top of indx array. Do only enough to get 
c  required nearest neighbors. 

      do 13 i=1,nn
         do 12 j=n,2,-1
            if(arrdis(arrindx(j)).lt.arrdis(arrindx(j-1))) then
               temp=arrindx(j-1)
               arrindx(j-1)=arrindx(j)
               arrindx(j)=temp
            endif

12 	continue

c       print*,"bubblesort",i     
c       print*,arrdis(arrindx(1)),arrdis(arrindx(2))
c       print*,arrindx(1),arrindx(2)

13     continue

        RETURN
        END
   




