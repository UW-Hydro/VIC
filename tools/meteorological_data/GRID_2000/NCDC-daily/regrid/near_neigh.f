      subroutine near_neigh(comp,cell,crnr,ncell,imiss,
     $nneigh,nstats,theta,ndist,neigh,rc,xcoord,ycoord)
      include 'array_dims.f'

c     Declare functions latlong and iso_weight
      real latlong
      external latlong

      real iso_weight
      external iso_weight
      
c     Declare shared block variables
c     Set array limits 

      integer ncell, nneigh, nstats
      integer rc(2), comp(DRC1*DRC2), imiss
      real cell, crnr(2)
c     max_neigh is the same as nneigh
      real xcoord(DNSTATS), ycoord(DNSTATS)
      real ndist(DNCELL,DMAXNEIGH), theta(DNCELL,DMAXNEIGH,DMAXNEIGH)
      integer neigh(DNCELL,DMAXNEIGH)

c     local variables
      integer i,j,k,m,n,cellnum
      real xloc, yloc
      
c     Declare variables used in subroutine indexx

      real dis(DNSTATS)
      integer indx(DNSTATS)

c     Change crnr(1) to the top right of grid 
      crnr(1)=crnr(1)+rc(1)*cell

c     Main loop. Loop over all rows(j)and columns(i).
c     Consult the mask variable comp.  If comp=1:
c     process that cell.  Otherwise check next cell.

      cellnum=0      
         do j= 1, rc(1)
            do i=1, rc(2) 
               if(comp((j-1)*rc(2)+i).eq.1)then
                  cellnum=cellnum+1

c     find center lat/lon coords of computational cell
c     using grid resolution variable "cell" and the
c     row (j) and column (i) of the comp cell 
c     in the mask file.
 
                yloc=crnr(2)+i*cell-cell/2.0
                xloc=crnr(1)-(j-1)*cell-cell/2.0

 
c     For the current cell, calculate the distances
c     to all station pts, and load to array "dis()".
                do n=1, nstats
2                   dis(n)=latlong(xcoord(n),ycoord(n),xloc,yloc)
      
                end do

c     Get the array pointer for the nearest neighbors
c     using the bubble sort subroutine indexx.
c     The array pointers are returned in array indx().
c     indx(1) stores the position in the input data array 
c     of the closest neighbor,for example.

                call indexx(nneigh,nstats,dis,indx)

c     Load the 2-dim arrays neigh()and ndist()
                do n=1, nneigh
                   neigh(cellnum,n)=indx(n)
                   ndist(cellnum,n)=1.0/(MAX(dis(indx(n)),1.0E-6))

                end do
c      print*, neigh(cellnum,1),neigh(cellnum,2)      
c     Get the isolation weight (1-cos(c)) for the current
c     cell and all combinations of nearest neighbors and
c     load to the 3-dim array theta().
c     The loop is a little wasteful in that theta(x,y,z)=
c     theta(x,z,y), but the coding makes later processing 
c     very straight-forward.

      do k=1, nneigh
         do m=1, nneigh
               a=latlong(xcoord(indx(k)),ycoord(indx(k)),
     $                     xcoord(indx(m)),ycoord(indx(m)))
               b=dis(indx(m))
               c=dis(indx(k))
               theta(cellnum,k,m)= max(1-iso_weight(a,b,c),1.0e-6)
          end do
        end do
c     End of if statement on comp variable status
              endif

c     End of i and j do loops
          end do
      end do
      return
      end
