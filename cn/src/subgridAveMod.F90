module subgridAveMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: subgridAveMod
!
! !DESCRIPTION:
! Utilities to perfrom subgrid averaging
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clm_varcon, only : spval

! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: p2c_1d   ! Perfrom an average from pfts to columns
  public :: p2c_2d   ! Perfrom an average from pfts to columns
  public :: p2g_1d   ! Perfrom an average from pfts to gridcells
  public :: p2g_2d   ! Perfrom an average from pfts to gridcells
  public :: c2g_1d   ! Perfrom an average from columns to gridcells
  public :: c2g_2d   ! Perfrom an average from columns to gridcells

! PRIVATE MEMBER FUNCTIONS

!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
! 1/31/13: Adopted for use in VIC by Michael Brunke
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2c_1d
!
! !INTERFACE:
  subroutine p2c_1d (lbp, ubp, lbc, ubc, parr, carr, p2c_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from pfts to columns.
! Averaging is only done for points that are not equal to "spval".
!
! !USES:
    use clm_varpar, only : max_pft_per_col
    use pftvarcon , only: npcropmin
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbp, ubp              ! beginning and ending pft
    integer , intent(in)  :: lbc, ubc              ! beginning and ending column
    real(r8), intent(in)  :: parr(lbp:ubp)         ! pft array
    real(r8), intent(out) :: carr(lbc:ubc)         ! column array
    character(len=*), intent(in) :: p2c_scale_type ! scale type
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
! Adapted by Michael Brunke for use in VIC 1/31/13
! Added ignoring crops in PFT loops (Michael Brunke) 2/7/14.
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: pi,p,c,index           ! indices
    real(r8) :: scale_p2c(lbp:ubp)     ! scale factor for column->landunit mapping
    logical  :: found                  ! temporary for error check
    real(r8) :: sumwt(lbc:ubc)         ! sum of weights
    real(r8), pointer :: wtcol(:)      ! weight of pft relative to column
    integer , pointer :: pcolumn(:)    ! column index of corresponding pft
    integer , pointer :: npfts(:)      ! number of pfts in column
    integer , pointer :: pfti(:)       ! initial pft index in column
    integer , pointer :: ivt(:)        ! pft vegetation index
!------------------------------------------------------------------------

    wtcol    => clm3%g%c%p%wtcol
    pcolumn  => clm3%g%c%p%column
    npfts    => clm3%g%c%npfts
    pfti     => clm3%g%c%pfti
    ivt      => clm3%g%c%p%itype

    if (p2c_scale_type == 'unity') then
       do p = lbp,ubp
          scale_p2c(p) = 1.0_r8
       end do
    else
       write(*,*)'p2c_1d error: scale type ',p2c_scale_type,' not supported'
       do p = lbp,ubp
          scale_p2c(p) = 1.0_r8
       end do
    end if

    carr(lbc:ubc) = spval
    sumwt(lbc:ubc) = 0._r8
    do p = lbp,ubp
       if (wtcol(p) /= 0._r8 .and. ivt(p) < npcropmin) then
          if (parr(p) /= spval) then
             c = pcolumn(p)
             if (sumwt(c) == 0._r8) carr(c) = 0._r8
             carr(c) = carr(c) + parr(p) * scale_p2c(p) * wtcol(p)
             sumwt(c) = sumwt(c) + wtcol(p)
          end if
       end if
    end do
    found = .false.
    do c = lbc,ubc
       if (sumwt(c) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = c
       else if (sumwt(c) /= 0._r8) then
          carr(c) = carr(c)/sumwt(c)
       end if
    end do
!    if (found) then
!       write(*,*)'p2c error: sumwt is greater than 1.0 at c= ',index
!       call endrun()
!    end if

  end subroutine p2c_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2c_2d
!
! !INTERFACE:
  subroutine p2c_2d (lbp, ubp, lbc, ubc, num2d, parr, carr, p2c_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from landunits to gridcells.
! Averaging is only done for points that are not equal to "spval".
!
! !USES:
    use clm_varpar, only : max_pft_per_col
    use pftvarcon , only: npcropmin
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbp, ubp              ! beginning and ending pft
    integer , intent(in)  :: lbc, ubc              ! beginning and ending column
    integer , intent(in)  :: num2d                 ! size of second dimension
    real(r8), intent(in)  :: parr(lbp:ubp,num2d)   ! pft array
    real(r8), intent(out) :: carr(lbp:ubp,num2d)   ! column array
    character(len=*), intent(in) :: p2c_scale_type ! scale type
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
! Adapted by Michael Brunke for use in VIC 1/31/13
! Added ignoring crops in PFT loops (Michael Brunke) 2/7/14.
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: j,pi,p,c,index         ! indices
    real(r8) :: scale_p2c(lbp:ubp)     ! scale factor for column->landunit mapping
    logical  :: found                  ! temporary for error check
    real(r8) :: sumwt(lbc:ubc)         ! sum of weights
    real(r8), pointer :: wtcol(:)      ! weight of pft relative to column
    integer , pointer :: pcolumn(:)    ! column index of corresponding pft
    integer , pointer :: npfts(:)      ! number of pfts in column
    integer , pointer :: pfti(:)       ! initial pft index in column
    integer , pointer :: ivt(:)        ! pft vegetation index
!------------------------------------------------------------------------

    wtcol    => clm3%g%c%p%wtcol
    pcolumn  => clm3%g%c%p%column
    npfts    => clm3%g%c%npfts
    pfti     => clm3%g%c%pfti
    ivt      => clm3%g%c%p%itype

    if (p2c_scale_type == 'unity') then
       do p = lbp,ubp
          scale_p2c(p) = 1.0_r8
       end do
    else
       write(*,*)'p2c_2d error: scale type ',p2c_scale_type,' not supported'
       do p = lbp,ubp
          scale_p2c(p) = 1.0_r8
       end do
    end if

    carr(:,:) = spval
    do j = 1,num2d
       sumwt(:) = 0._r8
       do p = lbp,ubp
          if (wtcol(p) /= 0._r8 .and. ivt(p) < npcropmin) then
             if (parr(p,j) /= spval) then
                c = pcolumn(p)
                if (sumwt(c) == 0._r8) carr(c,j) = 0._r8
                carr(c,j) = carr(c,j) + parr(p,j) * scale_p2c(p) * wtcol(p)
                sumwt(c) = sumwt(c) + wtcol(p)
             end if
          end if
       end do
       found = .false.
       do c = lbc,ubc
          if (sumwt(c) > 1.0_r8 + 1.e-6_r8) then
             found = .true.
             index = c
          else if (sumwt(c) /= 0._r8) then
             carr(c,j) = carr(c,j)/sumwt(c)
          end if
       end do
!       if (found) then
!          write(*,*)'p2c_2d error: sumwt is greater than 1.0 at c= ',index,' lev= ',j
!          call endrun()
!       end if
    end do 
  end subroutine p2c_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2g_1d
!
! !INTERFACE:
  subroutine p2g_1d(lbp, ubp, lbc, ubc, lbg, ubg, parr, garr, &
       p2c_scale_type, c2g_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from pfts to gridcells.
! Averaging is only done for points that are not equal to "spval".
!
! !USES:
    use clm_varpar, only : max_pft_per_gcell
    use pftvarcon , only: npcropmin
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbp, ubp            ! beginning and ending pft indices
    integer , intent(in)  :: lbc, ubc            ! beginning and ending column indices
    integer , intent(in)  :: lbg, ubg            ! beginning and ending gridcell indices
    real(r8), intent(in)  :: parr(lbp:ubp)       ! input pft array
    real(r8), intent(out) :: garr(lbg:ubg)       ! output gridcell array
    character(len=*), intent(in) :: p2c_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: c2g_scale_type ! scale factor type for averaging
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
! 1/31/13: Adapted for use in VIC by Michael Brunke
! Added ignoring crops in PFT loops (Michael Brunke) 2/7/14.
!
!  !LOCAL VARIABLES:
!EOP
    integer  :: pi,p,c,l,g,index       ! indices
    logical  :: found                  ! temporary for error check
    real(r8) :: scale_p2c(lbp:ubp)     ! scale factor
    real(r8) :: scale_c2g(lbc:ubc)     ! scale factor
    real(r8) :: sumwt(lbg:ubg)         ! sum of weights
    real(r8), pointer :: wtgcell(:)    ! weight of pfts relative to gridcells
    integer , pointer :: pcolumn(:)    ! column of corresponding pft
    integer , pointer :: plandunit(:)  ! landunit of corresponding pft
    integer , pointer :: pgridcell(:)  ! gridcell of corresponding pft
    integer , pointer :: npfts(:)      ! number of pfts in gridcell
    integer , pointer :: pfti(:)       ! initial pft index in gridcell
    integer , pointer :: ivt(:)        ! pft vegetation index
    integer , pointer :: ctype(:)      ! column type
    integer , pointer :: clandunit(:)  ! landunit of corresponding column
    integer , pointer :: ltype(:)      ! landunit type
    real(r8), pointer :: canyon_hwr(:) ! urban canyon height to width ratio
!------------------------------------------------------------------------

    ctype      => clm3%g%c%itype
    wtgcell   => clm3%g%c%p%wtgcell
    pcolumn   => clm3%g%c%p%column
    pgridcell => clm3%g%c%p%gridcell
    npfts     => clm3%g%npfts
    pfti      => clm3%g%pfti
    ivt      => clm3%g%c%p%itype

    if (c2g_scale_type == 'unity') then
       do l = lbc,ubc
          scale_c2g(l) = 1.0_r8
       end do
    else
       write(*,*)'p2g_1d error: scale type ',p2c_scale_type,' not supported'
       do c = lbc,ubc
          scale_c2g(l) = 1.0_r8
       end do
    end if

    if (p2c_scale_type == 'unity') then
       do p = lbp,ubp
          scale_p2c(p) = 1.0_r8
       end do
!    else
!       write(*,*)'p2g_1d error: scale type ',c2g_scale_type,' not supported'
!       call endrun()
    end if

    garr(:) = spval
    sumwt(:) = 0._r8
    do p = lbp,ubp
       if (wtgcell(p) /= 0._r8 .and. ivt(p) < npcropmin) then
          c = pcolumn(p)
          if (parr(p) /= spval) then
             g = pgridcell(p)
             if (sumwt(g) == 0._r8) garr(g) = 0._r8
             garr(g) = garr(g) + parr(p) * scale_p2c(p) * scale_c2g(l) * wtgcell(p)
             sumwt(g) = sumwt(g) + wtgcell(p)
          end if
       end if
    end do
    found = .false.
    do g = lbg, ubg
       if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = g
       else if (sumwt(g) /= 0._r8) then
          garr(g) = garr(g)/sumwt(g)
       end if
    end do
!    if (found) then
!       write(*,*)'p2g_1d error: sumwt is greater than 1.0 at g= ',index
!       call endrun()
!    end if

  end subroutine p2g_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2g_2d
!
! !INTERFACE:
  subroutine p2g_2d(lbp, ubp, lbc, ubc, lbg, ubg, num2d, &
       parr, garr, p2c_scale_type, c2g_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from pfts to gridcells.
! Averaging is only done for points that are not equal to "spval".
!
! !USES:
    use clm_varpar, only : max_pft_per_gcell
    use pftvarcon , only: npcropmin
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbp, ubp              ! beginning and ending pft indices
    integer , intent(in)  :: lbc, ubc              ! beginning and ending column indices
    integer , intent(in)  :: lbg, ubg              ! beginning and ending gridcell indices
    integer , intent(in)  :: num2d                 ! size of second dimension
    real(r8), intent(in)  :: parr(lbp:ubp,num2d)   ! input pft array
    real(r8), intent(out) :: garr(lbg:ubg,num2d)   ! output gridcell array
    character(len=*), intent(in) :: p2c_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: c2g_scale_type ! scale factor type for averaging
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
! 1/31/13: Adapted for use in VIC by Michael Brunke
! Added ignoring crops in PFT loops (Michael Brunke) 2/7/14.
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: j,pi,p,c,l,g,index     ! indices
    logical  :: found                  ! temporary for error check
    real(r8) :: scale_p2c(lbp:ubp)     ! scale factor
    real(r8) :: scale_c2g(lbc:ubc)     ! scale factor
    real(r8) :: sumwt(lbg:ubg)         ! sum of weights
    real(r8), pointer :: wtgcell(:)    ! weight of pfts relative to gridcells
    integer , pointer :: pcolumn(:)    ! column of corresponding pft
    integer , pointer :: plandunit(:)  ! landunit of corresponding pft
    integer , pointer :: pgridcell(:)  ! gridcell of corresponding pft
    integer , pointer :: npfts(:)      ! number of pfts in gridcell
    integer , pointer :: pfti(:)       ! initial pft index in gridcell
    integer , pointer :: ivt(:)        ! pft vegetation index
    integer , pointer :: clandunit(:)  ! landunit of corresponding column
    integer , pointer :: ctype(:)      ! column type
    integer , pointer :: ltype(:)      ! landunit type
    real(r8), pointer :: canyon_hwr(:) ! urban canyon height to width ratio
!------------------------------------------------------------------------

    ctype        => clm3%g%c%itype
    wtgcell      => clm3%g%c%p%wtgcell
    pcolumn      => clm3%g%c%p%column
    pgridcell    => clm3%g%c%p%gridcell
    npfts        => clm3%g%npfts
    pfti         => clm3%g%pfti
    ivt          => clm3%g%c%p%itype

    if (c2g_scale_type == 'unity') then
       do c = lbc,ubc
          scale_c2g(l) = 1.0_r8
       end do
    else
       write(*,*)'p2g_2d error: scale type ',c2g_scale_type,' not supported'
       do c = lbc,ubc
          scale_c2g(l) = 1.0_r8
       end do
    end if

    if (p2c_scale_type == 'unity') then
       do p = lbp,ubp
          scale_p2c(p) = 1.0_r8
       end do
    else
       write(*,*)'p2g_2d error: scale type ',p2c_scale_type,' not supported'
       do p = lbp,ubp
          scale_p2c(p) = 1.0_r8
       end do
    end if

    garr(:,:) = spval
    do j = 1,num2d
       sumwt(:) = 0._r8
       do p = lbp,ubp
          if (wtgcell(p) /= 0._r8 .and. ivt(p) < npcropmin) then
             c = pcolumn(p)
             if (parr(p,j) /= spval .and. scale_c2g(c) /= spval) then
                g = pgridcell(p)
                if (sumwt(g) == 0._r8) garr(g,j) = 0._r8
                garr(g,j) = garr(g,j) + parr(p,j) * scale_p2c(p) * scale_c2g(c) * wtgcell(p)
                sumwt(g) = sumwt(g) + wtgcell(p)
             end if
          end if
       end do
       found = .false.
       do g = lbg, ubg
          if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
             found = .true.
             index = g
          else if (sumwt(g) /= 0._r8) then
             garr(g,j) = garr(g,j)/sumwt(g)
          end if
       end do
!       if (found) then
!          write(*,*)'p2g_2d error: sumwt gt 1.0 at g/sumwt = ',index,sumwt(index)
!          call endrun()
!       end if
    end do

  end subroutine p2g_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: c2g_1d
!
! !INTERFACE:
  subroutine c2g_1d(lbc, ubc, lbg, ubg, carr, garr, c2g_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from columns to gridcells.
! Averaging is only done for points that are not equal to "spval".
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbc, ubc              ! beginning and ending column indices
    integer , intent(in)  :: lbg, ubg              ! beginning and ending landunit indices
    real(r8), intent(in)  :: carr(lbc:ubc)         ! input column array
    real(r8), intent(out) :: garr(lbg:ubg)         ! output gridcell array
    character(len=*), intent(in) :: c2g_scale_type ! scale factor type for averaging
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: ci,c,l,g,index         ! indices
    integer  :: max_col_per_gcell      ! max columns per gridcell; on the fly
    logical  :: found                  ! temporary for error check
    real(r8) :: scale_c2g(lbc:ubc)     ! scale factor
    real(r8) :: sumwt(lbg:ubg)         ! sum of weights
    real(r8), pointer :: wtgcell(:)    ! weight of columns relative to gridcells
    integer , pointer :: clandunit(:)  ! landunit of corresponding column
    integer , pointer :: cgridcell(:)  ! gridcell of corresponding column
    integer , pointer :: ncolumns(:)   ! number of columns in gridcell
    integer , pointer :: coli(:)       ! initial column index in gridcell
    integer , pointer :: ctype(:)      ! column type
!------------------------------------------------------------------------

    ctype      => clm3%g%c%itype
    wtgcell    => clm3%g%c%wtgcell
    clandunit  => clm3%g%c%landunit
    cgridcell  => clm3%g%c%gridcell
    ncolumns   => clm3%g%ncolumns
    coli       => clm3%g%coli

    if (c2g_scale_type == 'unity') then
       do c = lbc,ubc
          scale_c2g(c) = 1.0_r8
       end do
    else
       write(*,*)'c2g_1d error: scale type ',c2g_scale_type,' not supported'
       do c = lbc,ubc
          scale_c2g(c) = 1.0_r8
       end do
    end if

    garr(:) = spval
    sumwt(:) = 0._r8
    do c = lbc,ubc
       if ( wtgcell(c) /= 0._r8) then
          if (carr(c) /= spval .and. scale_c2g(c) /= spval) then
             g = cgridcell(c)
             if (sumwt(g) == 0._r8) garr(g) = 0._r8
             garr(g) = garr(g) + carr(c) * scale_c2g(c) * wtgcell(c)
             sumwt(g) = sumwt(g) + wtgcell(c)
          end if
       end if
    end do
    found = .false.
    do g = lbg, ubg
       if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
          found = .true.
          index = g
       else if (sumwt(g) /= 0._r8) then
          garr(g) = garr(g)/sumwt(g)
       end if
    end do
!    if (found) then
!       write(*,*)'c2g_1d error: sumwt is greater than 1.0 at g= ',index
!       call endrun()
!    end if

  end subroutine c2g_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: c2g_2d
!
! !INTERFACE:
  subroutine c2g_2d(lbc, ubc, lbg, ubg, num2d, carr, garr, &
       c2g_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from columns to gridcells.
! Averaging is only done for points that are not equal to "spval".
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbc, ubc              ! beginning and ending column indices
    integer , intent(in)  :: lbg, ubg              ! beginning and ending gridcell indices
    integer , intent(in)  :: num2d                 ! size of second dimension
    real(r8), intent(in)  :: carr(lbc:ubc,num2d)   ! input column array
    real(r8), intent(out) :: garr(lbg:ubg,num2d)   ! output gridcell array
    character(len=*), intent(in) :: c2g_scale_type ! scale factor type for averaging
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: j,ci,c,g,l,index       ! indices
    integer  :: max_col_per_gcell      ! max columns per gridcell; on the fly
    logical  :: found                  ! temporary for error check
    real(r8) :: scale_c2g(lbc:ubc)     ! scale factor
    real(r8) :: sumwt(lbg:ubg)         ! sum of weights
    real(r8), pointer :: wtgcell(:)    ! weight of columns relative to gridcells
    integer , pointer :: clandunit(:)  ! landunit of corresponding column
    integer , pointer :: cgridcell(:)  ! gridcell of corresponding column
    integer , pointer :: ncolumns(:)   ! number of columns in gridcell
    integer , pointer :: coli(:)       ! initial column index in gridcell
    integer , pointer :: ctype(:)      ! column type
    integer , pointer :: ltype(:)      ! landunit type
!------------------------------------------------------------------------

    ctype      => clm3%g%c%itype
    wtgcell    => clm3%g%c%wtgcell
    cgridcell  => clm3%g%c%gridcell
    ncolumns   => clm3%g%ncolumns
    coli       => clm3%g%coli

    if (c2g_scale_type == 'unity') then
       do c = lbc,ubc
          scale_c2g(c) = 1.0_r8
       end do
    else
       write(*,*)'c2g_2d error: scale type ',c2g_scale_type,' not supported'
       do c = lbc,ubc
          scale_c2g(c) = 1.0_r8
       end do
    end if

    garr(:,:) = spval
    do j = 1,num2d
       sumwt(:) = 0._r8
       do c = lbc,ubc
          if (wtgcell(c) /= 0._r8) then
             if (carr(c,j) /= spval .and. scale_c2g(c) /= spval) then
                g = cgridcell(c)
                if (sumwt(g) == 0._r8) garr(g,j) = 0._r8
                garr(g,j) = garr(g,j) + carr(c,j) * scale_c2g(c) * wtgcell(c)
                sumwt(g) = sumwt(g) + wtgcell(c)
             end if
          end if
       end do
       found = .false.
       do g = lbg, ubg
          if (sumwt(g) > 1.0_r8 + 1.e-6_r8) then
             found = .true.
             index = g
          else if (sumwt(g) /= 0._r8) then
             garr(g,j) = garr(g,j)/sumwt(g)
          end if
       end do
!       if (found) then
!          write(*,*)'c2g_2d error: sumwt is greater than 1.0 at g= ',index
!          call endrun()
!       end if
    end do

  end subroutine c2g_2d

end module subgridAveMod
