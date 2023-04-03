 module mem_tuv
   USE ModTuv, ONLY: nw,kz,ks

  implicit none

  integer,parameter :: firstBio = 10
  integer,parameter :: lastBio= 24
  integer,parameter :: nbio=lastBio-firstBio+1
  character(len=47),parameter, dimension(nbio) :: nameBio=(/ &
    'PAR, 400-700 nm, umol m-2 s-1                 ', & !01 10
    'Exponential decay, 14 nm/10                   ', & !02 11
    'DNA damage, in vitro (Setlow, 1974)           ', & !03 12
    'SCUP-mice (de Gruijl et al., 1993)            ', & !04 13
    'SCUP-human (de Gruijl and van der Leun, 1994) ', & !05 14
    'CIE human erythema (McKinlay and Diffey, 1987)', & !06 15
    'UV index (WMO, 1994)                          ', & !07 16
    'Erythema, humans (Anders et al., 1995)        ', & !08 17
    'Occupational TLV (ACGIH, 1992)                ', & !09 18
    'Phytoplankton (Boucher et al., 1994)          ', & !10 19
    'Phytoplankton, phaeo (Cullen et al., 1992)    ', & !11 20
    'Phytoplankton, proro (Cullen et al., 1992)    ', & !12 21
    'Cataract, pig (Oriowo et al., 2001)           ', & !13 22
    'Plant damage (Caldwell, 1971)                 ', & !14 23
    'Plant damage (Flint & Caldwell, 2003)         '  & !15 24
    /)
    !         1         2         3         4
    !12345678901234567890123456789012345678901234567890
    character(len=47),parameter, dimension(nbio) :: nameVar=(/ &
      'P47UM', & !01 10
      'EXD14', & !02 11
      'DNADV', & !03 12
      'SCUPM', & !04 13
      'SCUPH', & !05 14
      'CIEHE', & !06 15
      'UVIND', & !07 16
      'ERYTH', & !08 17
      'OCTLV', & !09 18
      'PHYTO', & !10 19
      'PHYPH', & !11 20
      'PHYPR', & !12 21
      'CAPIG', & !13 22
      'PLDCW', & !14 23
      'PLDFL'  & !15 24
      /)




  type mtuv
  REAL,POINTER,DIMENSION(:,:,:,:) :: g_tauaer
  REAL,POINTER,DIMENSION(:,:,:,:) :: g_wol
  REAL,POINTER,DIMENSION(:,:,:,:) :: g_gol
  REAL,POINTER,DIMENSION(:,:,:,:) :: g_taucld
  REAL,POINTER,DIMENSION(:,:,:,:) :: g_wcld
  REAL,POINTER,DIMENSION(:,:,:,:) :: g_gcld
  end type
  type (mtuv), allocatable, dimension(:) :: carma_tuv,carma_tuvm

  type tuvbio
    REAL,POINTER,DIMENSION(:,:) :: rateUV
  end type tuvbio
  type(tuvbio), allocatable, dimension(:,:) :: tuv_bio,tuv_biom

  INTEGER,DIMENSION(:),allocatable :: tuv2carma


 contains

  subroutine alloc_carma_tuv(carma_tuv, n1, n2, n3, n4, ng)
    use mem_globrad, only : ntotal, nlayer
    !    use mem_grid, only : ngrids,nnxp,nnyp

    type (mtuv), intent(inout) :: carma_tuv
    integer, intent(in) :: n1, n2, n3, n4, ng

    integer i

    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(alloc_carma_tuv)**"
    !print *,'LFR->test, nlayer: ',ntotal,nlayer; CALL flush(6)

    !allocate (carma_tuv(ngrids),stat=ierr)
    ! if (ierr/=0) then
    !     write(c0,"(i8)") ierr
    !     call fatal_error(h//"Allocating carma_tuv(ngrids) fails with stat="//trim(adjustl(c0)))
    ! endif

    !do i=1,ngrids
!     print*,'i=',i,ntotal,nlayer,nnxp(i),nnyp(i); call flush(6)

!     allocate (carma_tuv%g_tauaer(ntotal,nlayer,nnxp(i),nnyp(i)),stat=ierr)
     allocate (carma_tuv%g_tauaer(n1,n2,n3,n4),stat=ierr)
     if (ierr/=0) then
         write(c0,"(i8)") ierr
         call fatal_error(h//"Allocating carma_tuv%g_tauaer fails with stat="//trim(adjustl(c0)))
     endif

!     allocate (carma_tuv%g_wol   (ntotal,nlayer,nnxp(i),nnyp(i)),stat=ierr)
     allocate (carma_tuv%g_wol   (n1,n2,n3,n4),stat=ierr)
     if (ierr/=0) then
         write(c0,"(i8)") ierr
         call fatal_error(h//"Allocating carma_tuv%g_wol fails with stat="//trim(adjustl(c0)))
     endif

!     allocate (carma_tuv%g_gol   (ntotal,nlayer,nnxp(i),nnyp(i)),stat=ierr)
     allocate (carma_tuv%g_gol   (n1,n2,n3,n4),stat=ierr)
     if (ierr/=0) then
         write(c0,"(i8)") ierr
         call fatal_error(h//"Allocating carma_tuv%g_gol fails with stat="//trim(adjustl(c0)))
     endif

!     allocate (carma_tuv%g_taucld(ntotal,nlayer,nnxp(i),nnyp(i)),stat=ierr)
     allocate (carma_tuv%g_taucld(n1,n2,n3,n4),stat=ierr)
     if (ierr/=0) then
         write(c0,"(i8)") ierr
         call fatal_error(h//"Allocating carma_tuv%g_taucld fails with stat="//trim(adjustl(c0)))
     endif

!     allocate (carma_tuv%g_wcld  (ntotal,nlayer,nnxp(i),nnyp(i)),stat=ierr)
     allocate (carma_tuv%g_wcld  (n1,n2,n3,n4),stat=ierr)
     if (ierr/=0) then
         write(c0,"(i8)") ierr
         call fatal_error(h//"Allocating carma_tuv%g_wcld fails with stat="//trim(adjustl(c0)))
     endif

!     allocate (carma_tuv%g_gcld  (ntotal,nlayer,nnxp(i),nnyp(i)),stat=ierr)
     allocate (carma_tuv%g_gcld  (n1,n2,n3,n4),stat=ierr)
     if (ierr/=0) then
         write(c0,"(i8)") ierr
         call fatal_error(h//"Allocating carma_tuv%g_gcld fails with stat="//trim(adjustl(c0)))
     endif

!--(BRAMS-5.0-INI)---------------------------------------------------
     carma_tuv%g_tauaer = 0.
     carma_tuv%g_wol    = 0.
     carma_tuv%g_gol    = 0.
     carma_tuv%g_taucld = 0.
     carma_tuv%g_wcld   = 0.
     carma_tuv%g_gcld   = 0.
!--(BRAMS-5.0-FIM)---------------------------------------------------

  end subroutine alloc_carma_tuv

  subroutine alloc_tuv_bio(tuv_bio,n3,n4)
    integer, intent(in) :: n3,n4
    type (tuvbio), intent(inout) :: tuv_bio

    character(len=8) :: c0
    character(len=*), parameter :: h="**(alloc_tuv_bio)**"
    integer :: k,ierr

    allocate (tuv_bio%rateUV(n3,n4),stat=ierr)
    if (ierr/=0) then
      write(c0,"(i8)") ierr
      call fatal_error(h//"Allocating tuv_bio(ng,ks)%rateUV fails with stat="//trim(adjustl(c0)))
    endif
    tuv_bio%rateUV=0.0

  end subroutine alloc_tuv_bio

  subroutine nullify_tuv_bio(tuv_bio)
    implicit none
    type (tuvbio), intent(inout) :: tuv_bio

    if(associated(tuv_bio%rateUV)) nullify(tuv_bio%rateUV)

  end subroutine nullify_tuv_bio

  subroutine dealloc_tuv_bio(tuv_bio)

    implicit none
    type (tuvbio), intent(inout) :: tuv_bio

    if (associated(tuv_bio%rateUV))    deallocate (tuv_bio%rateUV   )

  end subroutine dealloc_tuv_bio

  subroutine nullify_carma_tuv(carma_tuv)

    implicit none
    type (mtuv), intent(inout) :: carma_tuv

    if (associated(carma_tuv%g_tauaer)) nullify (carma_tuv%g_tauaer)
    if (associated(carma_tuv%g_wol   )) nullify (carma_tuv%g_wol   )
    if (associated(carma_tuv%g_gol   )) nullify (carma_tuv%g_gol   )
    if (associated(carma_tuv%g_taucld)) nullify (carma_tuv%g_taucld)
    if (associated(carma_tuv%g_wcld  )) nullify (carma_tuv%g_wcld  )
    if (associated(carma_tuv%g_gcld  )) nullify (carma_tuv%g_gcld  )

  end subroutine nullify_carma_tuv

  subroutine filltab_tuv_bio(tuv_bio, tuv_biom, imean, n1, n2, n3, ng)

    use var_tables, only: InsertVTab
    use mem_stilt, only: iexev
    implicit none

    type (tuvbio), intent(in) :: tuv_bio, tuv_biom
    integer, intent(in)           :: imean, n1, n2, n3, ng
    integer, parameter :: i8 = selected_int_kind(14) !Kind for 64-bits Integer Numbers
    integer(kind=i8) :: npts
    real, pointer :: var, varm
    integer :: k

    ! Fill pointers to arrays into variable tables

    npts = n2*n3

    if (associated(tuv_bio%rateUV)) then
      call InsertVTab (tuv_bio%rateUV(:,:),tuv_biom%rateUV(:,:)  &
         ,ng, npts, imean,  &
         nameVar(n1)//' :2:hist:anal:mpti:mpt3:mpt2')
    end if

  end subroutine filltab_tuv_bio

end module mem_tuv
