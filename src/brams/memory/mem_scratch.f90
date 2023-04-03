!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_scratch


  !
  ! JP: uses excessive memory due to scr1, scr2 and scr3, that
  ! are dimensioned to store the full domain, not the mpi decomposed one
  !

  use grid_dims
  use node_mod, only: mchnum, master_num

  type scratch_vars
     real, pointer :: scr1(:)
     real, pointer :: vt3da(:)
     real, pointer :: vt3db(:)
     real, pointer :: vt3dc(:)
     real, pointer :: vt3dd(:)
     real, pointer :: vt3de(:)
     real, pointer :: vt3df(:)
     real, pointer :: vt3dg(:)
     real, pointer :: vt3dh(:)
     real, pointer :: vt3di(:)
     real, pointer :: vt3dj(:)
     real, pointer :: vt3dk(:)
     real, pointer :: vt3dl(:)
     real, pointer :: vt3dm(:)
     real, pointer :: vt3dn(:)
     real, pointer :: vt3do(:)
     real, pointer :: vt3dp(:)
!--(DMK-CCATT-INI)----------------------------------------------------------------
     real, pointer :: vt3dq(:)
!--(DMK-CCATT-FIM)----------------------------------------------------------------
     real, pointer :: vt2da(:)
     real, pointer :: vt2db(:)
     real, pointer :: vt2dc(:)
     real, pointer :: vt2dd(:)
     real, pointer :: vt2de(:)
     real, pointer :: vt2df(:)
  end type scratch_vars

  type (scratch_vars) :: scratch

  !-------------------------------------------------------------------
  REAL, ALLOCATABLE :: vctr1(:), vctr2(:), vctr3(:),          &
       vctr4(:),  vctr5(:),  vctr6(:),  vctr7(:),  vctr8(:),  &
       vctr9(:),  vctr10(:), vctr11(:), vctr12(:), vctr13(:), &
       vctr14(:), vctr15(:), vctr16(:), vctr17(:), vctr18(:), &
       vctr19(:), vctr20(:), vctr21(:), vctr22(:), vctr23(:), &
       vctr24(:), vctr25(:), vctr26(:), vctr27(:), vctr28(:), &
       vctr29(:), vctr30(:), vctr31(:), vctr32(:), vctr33(:), &
       vctr34(:), vctr35(:), vctr36(:), vctr37(:), vctr38(:), &
       vctr39(:), vctr40(:), vctr41(:)

  !---------------------------------------------------------------

contains

  subroutine alloc_scratch(nmzp, nmxp, nmyp, nnzp, nnxp, nnyp,  &
       maxgrds, ngrids, nzg, nzs, npatch, proc_type, maxnxp, maxnyp, maxnzp)

    ! FOR CATT
    use mem_aerad,   only: nwave      !intent(in)
    use mem_radiate, only: ilwrtyp, & !intent(in)
         iswrtyp                      !intent(in)

!--(DMK-CCATT-INI)-----------------------------------------------------
!    use catt_start, only: CATT        ! intent(in)
!--(DMK-CCATT-FIM)-----------------------------------------------------

    implicit none

    integer, intent(in ) :: maxgrds
    integer, intent(in ) :: nmzp(maxgrds)
    integer, intent(in ) :: nmxp(maxgrds)
    integer, intent(in ) :: nmyp(maxgrds)
    integer, intent(in ) :: nnzp(maxgrds)
    integer, intent(in ) :: nnxp(maxgrds)
    integer, intent(in ) :: nnyp(maxgrds)
    integer, intent(in ) :: ngrids
    integer, intent(in ) :: nzg
    integer, intent(in ) :: nzs
    integer, intent(in ) :: npatch
    integer, intent(in ) :: proc_type
    integer, intent(out) :: maxnxp
    integer, intent(out) :: maxnyp
    integer, intent(out) :: maxnzp

    ! Local Variables:
    integer :: ng,ntpts,ntpts1,ntpts2,ntptsx !,maxnxp,maxnyp,maxnzp
    ! For CATT
    integer :: ntpts_catt
    !         Find the maximum number of grid points needed for any grid.
    !           The max points in each direction are passed back for use by
    !           various nesting things.
    character(len=*), parameter :: h="**(alloc_scratch)**"

    maxnxp = 0
    maxnyp = 0
    maxnzp = 0
    ntpts=0
    ntpts2=0
    ntpts_catt = 0 ! CATT
    do ng=1,ngrids
       maxnxp = max(maxnxp,nnxp(ng))
       maxnyp = max(maxnyp,nnyp(ng))
       maxnzp = max(maxnzp,nnzp(ng))
       ntpts=max( nmxp(ng)*nmyp(ng)*nmzp(ng),ntpts )
       ntpts2=max( nmxp(ng)*nmyp(ng),ntpts2 )
    enddo
    ! scr1 and scr2 needs to be the max of a passed field
    ntptsx=max(maxnxp*maxnyp*maxnzp,ntpts2*nzg*npatch,ntpts2*nzs*npatch,maxnzp*40)+1000

    ! For CARMA
    !if (CATT == 1) then
    if (ilwrtyp==4 .or. iswrtyp==4) then
       ntpts_catt = max(ntptsx,ntpts_catt,(maxnxp*maxnyp*nwave))+1000
    else
       ntpts_catt = ntptsx
    endif

    if(proc_type==1) then
       ntpts1=1
    else
       ntpts1=ntpts
    endif

    ! Allocate arrays based on options (if necessary).
    !-scr1 and scr2 need to be allocated to full domain (even on compute nodes)
    ! to max(nx)*max(ny)*max(nz)
    !-do not need all these arrays if it is a master process in a parallel run,
    ! so just allocate some to 1 word.


    ! For CATT
    allocate (scratch%scr1 (ntpts_catt))
    scratch%scr1 = 0.
    allocate (scratch%vt3da(ntpts))
    scratch%vt3da = 0.
    allocate (scratch%vt3db(ntpts))
    scratch%vt3db = 0.
    allocate (scratch%vt3dc(ntpts1))
    scratch%vt3dc = 0.
    allocate (scratch%vt3dd(ntpts1))
    scratch%vt3dd = 0.
    allocate (scratch%vt3de(ntpts1))
    scratch%vt3de = 0.
    allocate (scratch%vt3df(ntpts1))
    scratch%vt3df = 0.
    allocate (scratch%vt3dg(ntpts1))
    scratch%vt3dg = 0.
    allocate (scratch%vt3dh(ntpts1))
    scratch%vt3dh = 0.
    allocate (scratch%vt3di(ntpts1))
    scratch%vt3di = 0.
    allocate (scratch%vt3dj(ntpts1))
    scratch%vt3dj = 0.
    allocate (scratch%vt3dk(ntpts1))
    scratch%vt3dk = 0.
    allocate (scratch%vt3dl(ntpts1))
    scratch%vt3dl = 0.
    allocate (scratch%vt3dm(ntpts1))
    scratch%vt3dm = 0.
    allocate (scratch%vt3dn(ntpts1))
    scratch%vt3dn = 0.
    allocate (scratch%vt3do(ntpts1))
    scratch%vt3do = 0.
    allocate (scratch%vt3dp(ntpts1))
    scratch%vt3dp = 0.

!--(DMK-CCATT-INI)----------------------------------------------------------------
    allocate (scratch%vt3dq(ntpts1))
    scratch%vt3dq = 0. 
!--(DMK-CCATT-FIM)----------------------------------------------------------------
    
    allocate (scratch%vt2da(ntpts2))
    scratch%vt2da = 0.
    allocate (scratch%vt2db(ntpts2))
    scratch%vt2db = 0.
    allocate (scratch%vt2dc(ntpts2))
    scratch%vt2dc = 0.
    allocate (scratch%vt2dd(ntpts2))
    scratch%vt2dd = 0.
    allocate (scratch%vt2de(ntpts2))
    scratch%vt2de = 0.
    allocate (scratch%vt2df(ntpts2))
    scratch%vt2df = 0.
    return
  end subroutine alloc_scratch

  !---------------------------------------------------------------

  subroutine nullify_scratch()

    implicit none

    ! Deallocate all scratch arrays

    if (associated(scratch%scr1 ))  nullify (scratch%scr1 )
!!$    if (associated(scratch%scr2 ))  nullify (scratch%scr2 )
    if (associated(scratch%vt3da))  nullify (scratch%vt3da)
    if (associated(scratch%vt3db))  nullify (scratch%vt3db)
    if (associated(scratch%vt3dc))  nullify (scratch%vt3dc)
    if (associated(scratch%vt3dd))  nullify (scratch%vt3dd)
    if (associated(scratch%vt3de))  nullify (scratch%vt3de)
    if (associated(scratch%vt3df))  nullify (scratch%vt3df)
    if (associated(scratch%vt3dg))  nullify (scratch%vt3dg)
    if (associated(scratch%vt3dh))  nullify (scratch%vt3dh)
    if (associated(scratch%vt3di))  nullify (scratch%vt3di)
    if (associated(scratch%vt3dj))  nullify (scratch%vt3dj)
    if (associated(scratch%vt3dk))  nullify (scratch%vt3dk)
    if (associated(scratch%vt3dl))  nullify (scratch%vt3dl)
    if (associated(scratch%vt3dm))  nullify (scratch%vt3dm)
    if (associated(scratch%vt3dn))  nullify (scratch%vt3dn)
    if (associated(scratch%vt3do))  nullify (scratch%vt3do)
    if (associated(scratch%vt3dp))  nullify (scratch%vt3dp)
    
!--(DMK-CCATT-INI)----------------------------------------------------------------
    if (associated(scratch%vt3dq))  nullify (scratch%vt3dq)
!--(DMK-CCATT-FIM)----------------------------------------------------------------
    
    if (associated(scratch%vt2da))  nullify (scratch%vt2da)
    if (associated(scratch%vt2db))  nullify (scratch%vt2db)
    if (associated(scratch%vt2dc))  nullify (scratch%vt2dc)
    if (associated(scratch%vt2dd))  nullify (scratch%vt2dd)
    if (associated(scratch%vt2de))  nullify (scratch%vt2de)
    if (associated(scratch%vt2df))  nullify (scratch%vt2df)


    return
  end subroutine nullify_scratch
  !---------------------------------------------------------------

  subroutine dealloc_scratch()

    implicit none

    ! Deallocate all scratch arrays

    if (associated(scratch%scr1 ))  deallocate (scratch%scr1 )
!!$    if (associated(scratch%scr2 ))  deallocate (scratch%scr2 )
    if (associated(scratch%vt3da))  deallocate (scratch%vt3da)
    if (associated(scratch%vt3db))  deallocate (scratch%vt3db)
    if (associated(scratch%vt3dc))  deallocate (scratch%vt3dc)
    if (associated(scratch%vt3dd))  deallocate (scratch%vt3dd)
    if (associated(scratch%vt3de))  deallocate (scratch%vt3de)
    if (associated(scratch%vt3df))  deallocate (scratch%vt3df)
    if (associated(scratch%vt3dg))  deallocate (scratch%vt3dg)
    if (associated(scratch%vt3dh))  deallocate (scratch%vt3dh)
    if (associated(scratch%vt3di))  deallocate (scratch%vt3di)
    if (associated(scratch%vt3dj))  deallocate (scratch%vt3dj)
    if (associated(scratch%vt3dk))  deallocate (scratch%vt3dk)
    if (associated(scratch%vt3dl))  deallocate (scratch%vt3dl)
    if (associated(scratch%vt3dm))  deallocate (scratch%vt3dm)
    if (associated(scratch%vt3dn))  deallocate (scratch%vt3dn)
    if (associated(scratch%vt3do))  deallocate (scratch%vt3do)
    if (associated(scratch%vt3dp))  deallocate (scratch%vt3dp)
    
!--(DMK-CCATT-INI)----------------------------------------------------------------
    if (associated(scratch%vt3dq))  deallocate (scratch%vt3dq)
!--(DMK-CCATT-FIM)---------------------------------------------------------------- 
    
    if (associated(scratch%vt2da))  deallocate (scratch%vt2da)
    if (associated(scratch%vt2db))  deallocate (scratch%vt2db)
    if (associated(scratch%vt2dc))  deallocate (scratch%vt2dc)
    if (associated(scratch%vt2dd))  deallocate (scratch%vt2dd)
    if (associated(scratch%vt2de))  deallocate (scratch%vt2de)
    if (associated(scratch%vt2df))  deallocate (scratch%vt2df)

    return
  end subroutine dealloc_scratch

  !--------------------------------------------------------------------

  SUBROUTINE createVctr(ngrids,nnxp, nnyp, nnzp)

    IMPLICIT NONE
    ! Arguments:
    INTEGER, INTENT(IN) :: ngrids
!TO    INTEGER, INTENT(IN) :: nnxp(maxgrds) ! From RAMSIN
!TO    INTEGER, INTENT(IN) :: nnyp(maxgrds) ! From RAMSIN
!TO    INTEGER, INTENT(IN) :: nnzp(maxgrds) ! From RAMSIN
    INTEGER, INTENT(IN),TARGET :: nnxp(ngrids) ! From RAMSIN
    INTEGER, INTENT(IN),TARGET :: nnyp(ngrids) ! From RAMSIN
    INTEGER, INTENT(IN),TARGET :: nnzp(ngrids) ! From RAMSIN
    ! Local variables:
    INTEGER :: ierr, maxx, maxy, maxz, maxxyz

    maxx   = maxval(nnxp(1:ngrids))
    maxy   = maxval(nnyp(1:ngrids))
    maxz   = maxval(nnzp(1:ngrids))
    maxxyz = max(maxx, maxy, maxz) + 2

    ALLOCATE(vctr1(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr1 (createVctr)")
    ALLOCATE(vctr2(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr2 (createVctr)")
    ALLOCATE(vctr3(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr3 (createVctr)")
    ALLOCATE(vctr4(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr4 (createVctr)")
    ALLOCATE(vctr5(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr5 (createVctr)")
    ALLOCATE(vctr6(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr6 (createVctr)")
    ALLOCATE(vctr7(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr7 (createVctr)")
    ALLOCATE(vctr8(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr8 (createVctr)")
    ALLOCATE(vctr9(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr9 (createVctr)")
    ALLOCATE(vctr10(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr10 (createVctr)")
    ALLOCATE(vctr11(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr11 (createVctr)")
    ALLOCATE(vctr12(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr12 (createVctr)")
    ALLOCATE(vctr13(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr13 (createVctr)")
    ALLOCATE(vctr14(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr14 (createVctr)")
    ALLOCATE(vctr15(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr15 (createVctr)")
    ALLOCATE(vctr16(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr16 (createVctr)")
    ALLOCATE(vctr17(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr17 (createVctr)")
    ALLOCATE(vctr18(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr18 (createVctr)")
    ALLOCATE(vctr19(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr19 (createVctr)")
    ALLOCATE(vctr20(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr20 (createVctr)")
    ALLOCATE(vctr21(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr21 (createVctr)")
    ALLOCATE(vctr22(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr22 (createVctr)")
    ALLOCATE(vctr23(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr23 (createVctr)")
    ALLOCATE(vctr24(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr24 (createVctr)")
    ALLOCATE(vctr25(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr25 (createVctr)")
    ALLOCATE(vctr26(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr26 (createVctr)")
    ALLOCATE(vctr27(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr27 (createVctr)")
    ALLOCATE(vctr28(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr28 (createVctr)")
    ALLOCATE(vctr29(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr29 (createVctr)")
    ALLOCATE(vctr30(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr30 (createVctr)")
    ALLOCATE(vctr31(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr31 (createVctr)")
    ALLOCATE(vctr32(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr32 (createVctr)")
    ALLOCATE(vctr33(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr33 (createVctr)")
    ALLOCATE(vctr34(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr34 (createVctr)")
    ALLOCATE(vctr35(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr35 (createVctr)")
    ALLOCATE(vctr36(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr36 (createVctr)")
    ALLOCATE(vctr37(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr37 (createVctr)")
    ALLOCATE(vctr38(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr38 (createVctr)")
    ALLOCATE(vctr39(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr39 (createVctr)")
    ALLOCATE(vctr40(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr40 (createVctr)")
    ALLOCATE(vctr41(maxxyz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating vctr41 (createVctr)")

!!$    ALLOCATE(ivctr(maxxyz), STAT=ierr)
!!$    IF (ierr/=0) CALL fatal_error("ERROR allocating ivctr (createVctr)")

!--(DMK-LFR NEC-SX6)----------------------------------------------
    vctr1  = 0.
    vctr2  = 0.
    vctr3  = 0.
    vctr4  = 0.
    vctr5  = 0.
    vctr6  = 0.
    vctr7  = 0.     
    vctr9  = 0.
    vctr10 = 0.
    vctr11 = 0.
    vctr12 = 0.
    vctr13 = 0.
    vctr14 = 0.
    vctr15 = 0.
    vctr16 = 0.
    vctr17 = 0.
    vctr19 = 0.
    vctr20 = 0.
    vctr21 = 0.
    vctr22 = 0.
    vctr23 = 0.
    vctr24 = 0.
    vctr25 = 0.
    vctr26 = 0.
    vctr27 = 0.
    vctr28 = 0.
    vctr29 = 0.
    vctr30 = 0.
    vctr31 = 0.
    vctr32 = 0.
    vctr33 = 0.
    vctr34 = 0.
    vctr35 = 0. 
    vctr36 = 0.
    vctr37 = 0.
    vctr38 = 0.
    vctr39 = 0.
    vctr40 = 0.
    vctr41 = 0.
!--(DMK-LFR NEC-SX6)----------------------------------------------



  END SUBROUTINE createVctr

  !--------------------------------------------------------------------

  SUBROUTINE destroyVctr()
    IMPLICIT NONE
    ! Local variables:
    INTEGER :: ierr

    DEALLOCATE(vctr1, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr1 (destroyVctr)")
    DEALLOCATE(vctr2, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr2 (destroyVctr)")
    DEALLOCATE(vctr3, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr3 (destroyVctr)")
    DEALLOCATE(vctr4, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr4 (destroyVctr)")
    DEALLOCATE(vctr5, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr5 (destroyVctr)")
    DEALLOCATE(vctr6, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr6 (destroyVctr)")
    DEALLOCATE(vctr7, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr7 (destroyVctr)")
    DEALLOCATE(vctr8, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr8 (destroyVctr)")
    DEALLOCATE(vctr9, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr9 (destroyVctr)")
    DEALLOCATE(vctr10, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr10 (destroyVctr)")
    DEALLOCATE(vctr11, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr11 (destroyVctr)")
    DEALLOCATE(vctr12, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr12 (destroyVctr)")
    DEALLOCATE(vctr13, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr13 (destroyVctr)")
    DEALLOCATE(vctr14, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr14 (destroyVctr)")
    DEALLOCATE(vctr15, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr15 (destroyVctr)")
    DEALLOCATE(vctr16, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr16 (destroyVctr)")
    DEALLOCATE(vctr17, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr17 (destroyVctr)")
    DEALLOCATE(vctr18, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr18 (destroyVctr)")
    DEALLOCATE(vctr19, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr19 (destroyVctr)")
    DEALLOCATE(vctr20, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr20 (destroyVctr)")
    DEALLOCATE(vctr21, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr21 (destroyVctr)")
    DEALLOCATE(vctr22, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr22 (destroyVctr)")
    DEALLOCATE(vctr23, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr23 (destroyVctr)")
    DEALLOCATE(vctr24, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr24 (destroyVctr)")
    DEALLOCATE(vctr25, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr25 (destroyVctr)")
    DEALLOCATE(vctr26, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr26 (destroyVctr)")
    DEALLOCATE(vctr27, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr27 (destroyVctr)")
    DEALLOCATE(vctr28, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr28 (destroyVctr)")
    DEALLOCATE(vctr29, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr29 (destroyVctr)")
    DEALLOCATE(vctr30, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr30 (destroyVctr)")
    DEALLOCATE(vctr31, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr31 (destroyVctr)")
    DEALLOCATE(vctr32, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr32 (destroyVctr)")
    DEALLOCATE(vctr33, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr33 (destroyVctr)")
    DEALLOCATE(vctr34, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr34 (destroyVctr)")
    DEALLOCATE(vctr35, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr35 (destroyVctr)")
    DEALLOCATE(vctr36, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr36 (destroyVctr)")
    DEALLOCATE(vctr37, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr37 (destroyVctr)")
    DEALLOCATE(vctr38, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr38 (destroyVctr)")
    DEALLOCATE(vctr39, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr39 (destroyVctr)")
    DEALLOCATE(vctr40, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr40 (destroyVctr)")
    DEALLOCATE(vctr41, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating vctr41 (destroyVctr)")

  END SUBROUTINE destroyVctr

end module mem_scratch
