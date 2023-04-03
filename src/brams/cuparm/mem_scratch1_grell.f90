! Module necessary to Grell Cumulus param.
! Scratch variables - module 1

MODULE mem_scratch1_grell

     !srf- feb-05-2002 : Variables for cumulus transport scheme
     !                   adapted in july-15-2002 for 5.x version
     !
     !ngrids_cp = numero de grades onde a parametrizacao de cumulus e' usada
     !
     implicit none

     logical :: scratch1_grell_alloc=.false.
     integer :: iruncon = 0

    !4d dependence (mgmxp,mgmyp,maxiens,ngrids_cp)
     INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:) ::  &
          ierr4d,                             & ! ierr4d=0 there is convection
          jmin4d,                             & ! ETL
          kdet4d,                             & ! Detrainemnt level
          k224d,                              & !
          kbcon4d,                            & ! Cloud base
          ktop4d,                             & ! Cloud top
          kstabi4d,                           & ! Cloud base
          kstabm4d,                           & ! Cloud top
          kpbl4d,                             &! PBL height
	  klcl4d

     REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: &
          xmb4d,                              & ! updraft mass flux
          pwav4d,                             &
          edt4d,                              &
	  cprr4d,                             &
	  sigma4d                              ! MFLX_DOWN/MFLX_UP


     !5d dependece  (mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp)
     REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) ::   &
          enup5d,                             & ! entrain up   rate
          endn5d,                             & ! entrain down rate
          deup5d,                             & ! detrain up   rate
          dedn5d,                             & ! detrain down rate
           zup5d,                             & ! norm mass flux up
           zdn5d,                             & ! norm mass flux down
          zcup5d,                             & ! z level
          pcup5d,                             & !p level
          prup5d,                             &
          prdn5d,                             &
         clwup5d,                             &
           tup5d,                             &
!-only for GF scheme
          up_massentr5d,                      &
	  up_massdetr5d,                      &
	  dd_massentr5d,                      &
	  dd_massdetr5d,                      &
          conv_cld_fr5d
!-only for GF scheme


  !END TYPE scratch1_grell_vars

CONTAINS

  SUBROUTINE alloc_scratch1_grell !(scratch1_grell)
    use dump, only: &
      dumpMessage

    USE mem_grell_param, ONLY : mgmxp,  & ! INTENT(IN)
         mgmyp,                         & ! INTENT(IN)
         mgmzp,                         & ! INTENT(IN)
         maxiens,                       & ! INTENT(IN)
         ngrids_cp                        ! INTENT(IN)

    !use mem_cuparm     , only: NNQPARM   ! INTENT(IN)
    !use shcu_vars_const, only: NNSHCU    ! INTENT(IN)
    use ccatt_start, only: ccatt
    IMPLICIT NONE

    include "constants.f90"
    character(len=*), parameter :: header='**(alloc_scratch1_grell)**'

    !-srf this array must be allocated anyway because it will be used to stor
    !-srf the sub-grid scale cloud liquid water to feed the radiation schemes.
    ALLOCATE(clwup5d(mgmzp, mgmxp, mgmyp, maxiens, ngrids_cp)); clwup5d=0.0
    ALLOCATE(ierr4d (mgmxp, mgmyp, maxiens, ngrids_cp))       ; ierr4d =-999.
    ALLOCATE(xmb4d  (mgmxp, mgmyp, maxiens, ngrids_cp))       ; xmb4d  =0.0
    ALLOCATE(zup5d  (mgmzp, mgmxp, mgmyp, maxiens, ngrids_cp)); zup5d  =0.0

    !- allocate anyway because of STILT as well.
    !IF(CCATT == 0 ) return

    if(scratch1_grell_alloc) then
       print *,'ERROR: scratch1_grell already allocated'
       print *,'Routine: alloc_scratch1_grell File: mem_scratch1_grell.f90'
       !call fatal_error("Fatal error at SCRATCH_GRELL allocation.") !stop
       iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion,c_fatal, &
          "error at SCRATCH_GRELL allocation")

    end if

    ALLOCATE (  jmin4d(mgmxp, mgmyp, maxiens, ngrids_cp));  jmin4d=0.0
    ALLOCATE (  kdet4d(mgmxp, mgmyp, maxiens, ngrids_cp));  kdet4d=0.0
    ALLOCATE (   k224d(mgmxp, mgmyp, maxiens, ngrids_cp));   k224d=0.0
    ALLOCATE ( kbcon4d(mgmxp, mgmyp, maxiens, ngrids_cp)); kbcon4d=0.0
    ALLOCATE (  ktop4d(mgmxp, mgmyp, maxiens, ngrids_cp));  ktop4d=0.0
    ALLOCATE (kstabi4d(mgmxp, mgmyp, maxiens, ngrids_cp));kstabi4d=0.0
    ALLOCATE (kstabm4d(mgmxp, mgmyp, maxiens, ngrids_cp));kstabm4d=0.0
    ALLOCATE (  kpbl4d(mgmxp, mgmyp, maxiens, ngrids_cp));  kpbl4d=0.0
    ALLOCATE (  klcl4d(mgmxp, mgmyp, maxiens, ngrids_cp));  klcl4d=0.0

    ALLOCATE (edt4d(mgmxp, mgmyp, maxiens, ngrids_cp));edt4d=0.0
    ALLOCATE (pwav4d(mgmxp, mgmyp, maxiens, ngrids_cp));pwav4d=0.0
    ALLOCATE (sigma4d(mgmxp, mgmyp, maxiens, ngrids_cp));sigma4d=0.0
    ALLOCATE (cprr4d(mgmxp, mgmyp, maxiens, ngrids_cp));cprr4d=0.0


!-only for GD scheme
    ALLOCATE (enup5d(mgmzp, mgmxp, mgmyp, maxiens, ngrids_cp));  enup5d=0.0
    ALLOCATE (endn5d(mgmzp, mgmxp, mgmyp, maxiens, ngrids_cp));  endn5d=0.0
    ALLOCATE (deup5d(mgmzp, mgmxp, mgmyp, maxiens, ngrids_cp));  deup5d=0.0
    ALLOCATE (dedn5d(mgmzp, mgmxp, mgmyp, maxiens, ngrids_cp));  dedn5d=0.0
    ALLOCATE (zcup5d(mgmzp, mgmxp, mgmyp, maxiens, ngrids_cp));  zcup5d=0.0

!-only for GF scheme
    ALLOCATE (pcup5d       (mgmzp, mgmxp, mgmyp, maxiens, ngrids_cp)); pcup5d=0.0
    ALLOCATE (up_massentr5d(mgmzp, mgmxp, mgmyp, maxiens, ngrids_cp)); up_massentr5d=0.0
    ALLOCATE (up_massdetr5d(mgmzp, mgmxp, mgmyp, maxiens, ngrids_cp)); up_massdetr5d=0.0
    ALLOCATE (dd_massentr5d(mgmzp, mgmxp, mgmyp, maxiens, ngrids_cp)); dd_massentr5d=0.0
    ALLOCATE (dd_massdetr5d(mgmzp, mgmxp, mgmyp, maxiens, ngrids_cp)); dd_massdetr5d=0.0
    ALLOCATE (conv_cld_fr5d(mgmzp, mgmxp, mgmyp, maxiens, ngrids_cp)); conv_cld_fr5d=0.0


!-for all schemes
    ALLOCATE ( zdn5d(mgmzp, mgmxp, mgmyp, maxiens, ngrids_cp));   zdn5d=0.0
    ALLOCATE (prup5d(mgmzp, mgmxp, mgmyp, maxiens, ngrids_cp));  prup5d=0.0
    ALLOCATE (prdn5d(mgmzp, mgmxp, mgmyp, maxiens, ngrids_cp));  prdn5d=0.0
    ALLOCATE ( tup5d(mgmzp, mgmxp, mgmyp, maxiens, ngrids_cp));   tup5d=0.0




    scratch1_grell_alloc=.true.

  END SUBROUTINE alloc_scratch1_grell

  SUBROUTINE dealloc_scratch1_grell !(scratch1_grell)

    IMPLICIT NONE

    DEALLOCATE (ierr4d)
    DEALLOCATE (jmin4d)
    DEALLOCATE (kdet4d)
    DEALLOCATE (k224d)
    DEALLOCATE (kbcon4d)
    DEALLOCATE (ktop4d)
    DEALLOCATE (kstabi4d)
    DEALLOCATE (kstabm4d)
    DEALLOCATE (kpbl4d)

    DEALLOCATE (xmb4d)
    DEALLOCATE (edt4d)
    DEALLOCATE (pwav4d)

    DEALLOCATE (zcup5d)
    DEALLOCATE (pcup5d)

    DEALLOCATE (enup5d)
    DEALLOCATE (endn5d)
    DEALLOCATE (deup5d)
    DEALLOCATE (dedn5d)
    DEALLOCATE (zup5d)
    DEALLOCATE (zdn5d)

    DEALLOCATE ( prup5d)
    DEALLOCATE (clwup5d)
    DEALLOCATE (  tup5d)

    DEALLOCATE (up_massentr5d)
    DEALLOCATE (up_massdetr5d)
    DEALLOCATE (dd_massentr5d)
    DEALLOCATE (dd_massdetr5d)
    DEALLOCATE (conv_cld_fr5d)
    DEALLOCATE (prdn5d)
  END SUBROUTINE dealloc_scratch1_grell

!!$  SUBROUTINE filltab_scratch1_grell()
!!$
!!$    USE var_tables
!!$
!!$    IMPLICIT NONE
!!$
!!$    ! Can't think of anything to do here...
!!$
!!$    RETURN
!!$  END SUBROUTINE filltab_scratch1_grell

! subroutine zero_scratch1_grell()
!
!   implicit none
!
!    zcup5d=0.
!    pcup5d=0.
!    prup5d=0.
!    clwup5d=0.
!    tup5d=0.
!
!  end subroutine zero_scratch1_grell




END MODULE mem_scratch1_grell
