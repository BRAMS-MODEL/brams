! Module necessary to Grell Cumulus param.
! Scratch variables - module 2

MODULE mem_scratch2_grell

  !TYPE scratch2_grell_vars

     !                   adapted in july-15-2002 for 5.x version
     !
     ! 3d dependence (mgmxp,mgmyp,ensdim)
     REAL, ALLOCATABLE, DIMENSION(:,:,:) :: massfln

     ! 2d dependence (mgmxp,mgmzp)
     REAL, ALLOCATABLE, DIMENSION(:,:) :: T,      &
          Q,                                  &
          P,                                  &
          PO,                                 &
          TN,                                 &
          QO,                                 &
          OUTT,                               &
          OUTQ,                               &
          outqc,                              &
          US_Grell,                           & ! Substitui US original
          VS_Grell,                           & ! Substitui VS original
          omeg
     ! 2d dependence (mgmxp, mgmyp)
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: KDT, &
          iact_gr,                            &
          iact_old_gr
     REAL, ALLOCATABLE, DIMENSION(:,:) :: xland,  &
          tkeg,                                   &
          rcpg,                                   &
          massflx

     ! 1d dependence (mgmxp)
     REAL, ALLOCATABLE, DIMENSION(:) :: mconv,    &
          umean,                              &
          vmean,                              &
          pmean,                              &
          direction
     REAL, ALLOCATABLE, DIMENSION(:) ::       &
           AA0,                               &
          PRET,                               &
          PSUR,                               &
          TER11,                              &
          !srf-20082006 for training
          glatg,                              &
	  glong 

     INTEGER, ALLOCATABLE, DIMENSION(:) :: KDET

  !END TYPE scratch2_grell_vars

CONTAINS

  SUBROUTINE alloc_scratch2_grell !(scratch2_grell)

    USE mem_grell_param, ONLY : mgmxp,  & ! INTENT(IN)
         mgmyp,                         & ! INTENT(IN)
         mgmzp,                         & ! INTENT(IN)
         ensdim                           ! INTENT(IN)
    USE node_mod, ONLY : mynum, &   ! INTENT(IN)
         MXP,                   &   ! INTENT(IN)
         MYP,                   &   ! INTENT(IN)
         MZP,                   &   ! INTENT(IN)
         IA,                    &   ! INTENT(IN)
         IZ,                    &   ! INTENT(IN)
         JA,                    &   ! INTENT(IN)
         JZ,                    &   ! INTENT(IN)
         I0,                    &   ! INTENT(IN)
         J0                         ! INTENT(IN)

    IMPLICIT NONE
    !TYPE (scratch2_grell_vars) :: scratch2_grell

    INTEGER :: i

    ALLOCATE (massfln(mgmxp, mgmyp, ensdim)) ; massfln =0.0

    ALLOCATE (mconv    (mgmxp))              ;mconv     = 0.0
    ALLOCATE (umean    (mgmxp)) 	     ;umean     = 0.0
    ALLOCATE (vmean    (mgmxp)) 	     ;vmean     = 0.0
    ALLOCATE (pmean    (mgmxp)) 	     ;pmean     = 0.0
    ALLOCATE (direction(mgmxp)) 	     ;direction = 0.0
    ALLOCATE (PRET     (mgmxp)) 	     ;PRET      = 0.0
    ALLOCATE (PSUR     (mgmxp)) 	     ;PSUR      = 0.0
    ALLOCATE (TER11    (mgmxp)) 	     ;TER11     = 0.0
!srf 20082006				     ;		= 0.0
    ALLOCATE (glatg    (mgmxp)) 	     ;glatg	= 0.0
    ALLOCATE (glong    (mgmxp)) 	     ;glong	= 0.0
!					     ;		= 0.0
    ALLOCATE (AA0      (mgmxp)) 	     ;AA0 	= 0.0
    ALLOCATE (KDET     (mgmxp)) 	     ;KDET	= 0.0
!					     ;		= 0.0
    ALLOCATE (T          (mgmxp, mgmzp))     ;T 	= 0.0 
    ALLOCATE (Q          (mgmxp, mgmzp))     ;Q 	= 0.0 
    ALLOCATE (P          (mgmxp, mgmzp))     ;P 	= 0.0 
    ALLOCATE (PO         (mgmxp, mgmzp))     ;PO	= 0.0 
    ALLOCATE (TN         (mgmxp, mgmzp))     ;TN	= 0.0 
    ALLOCATE (QO         (mgmxp, mgmzp))     ;QO	= 0.0 
    ALLOCATE (OUTT       (mgmxp, mgmzp))     ;OUTT	= 0.0 
    ALLOCATE (OUTQ       (mgmxp, mgmzp))     ;OUTQ	= 0.0 
    ALLOCATE (outqc      (mgmxp, mgmzp))     ;outqc	= 0.0 
    ALLOCATE (US_Grell   (mgmxp, mgmzp))     ;US_Grell  = 0.0 
    ALLOCATE (VS_Grell   (mgmxp, mgmzp))     ;VS_Grell  = 0.0 
    ALLOCATE (omeg       (mgmxp, mgmzp))     ;omeg	= 0.0 
    ALLOCATE (KDT        (mgmxp, mgmyp))     ;KDT	= 0.0 
    ALLOCATE (xland      (mgmxp, mgmyp))     ;xland	= 0.0 
    ALLOCATE (massflx    (mgmxp, mgmyp))     ;massflx	= 0.0 
    ALLOCATE (iact_gr    (mgmxp, mgmyp))     ;iact_gr	=0.0
    ALLOCATE (iact_old_gr(mgmxp, mgmyp))     ;iact_old_gr=0.0

    ALLOCATE (tkeg       (mgmxp, mgmzp))  ;tkeg=0.0
    ALLOCATE (rcpg       (mgmxp, mgmzp))  ;rcpg=0.0

    RETURN
  END SUBROUTINE alloc_scratch2_grell

  SUBROUTINE dealloc_scratch2_grell !(scratch2_grell)

    IMPLICIT NONE
    !TYPE (scratch2_grell_vars) :: scratch2_grell

    DEALLOCATE (massfln)

    DEALLOCATE (mconv)
    DEALLOCATE (umean)
    DEALLOCATE (vmean)
    DEALLOCATE (pmean)
    DEALLOCATE (direction)

    DEALLOCATE (T)
    DEALLOCATE (Q)
    DEALLOCATE (P)
    DEALLOCATE (PO)
    DEALLOCATE (TN)
    DEALLOCATE (QO)
    DEALLOCATE (OUTT)
    DEALLOCATE (OUTQ)
    DEALLOCATE (outqc)
    DEALLOCATE (PRET)
    DEALLOCATE (PSUR)
    DEALLOCATE (TER11)
    DEALLOCATE (US_Grell)
    DEALLOCATE (VS_Grell)
    DEALLOCATE (omeg)
    DEALLOCATE (AA0)

    DEALLOCATE (KDET)
    DEALLOCATE (KDT)

    DEALLOCATE (xland)
    DEALLOCATE (massflx)
    DEALLOCATE (iact_gr)
    DEALLOCATE (iact_old_gr)

    DEALLOCATE (tkeg)
    DEALLOCATE (rcpg)

    DEALLOCATE (glatg) 
    DEALLOCATE (glong) 

    RETURN
  END SUBROUTINE dealloc_scratch2_grell

!  SUBROUTINE zero_scratch2_grell()
!	massfln=0.
!	T=0.
!	Q=0.
!	P=0.
!	PO=0.
!	TN=0.
!	QO=0.
!	OUTT=0.
!	OUTQ=0.
!	outqc=0.
!	US_Grell=0.
!	VS_Grell=0.
!	omeg=0.
!	xland=0.
!	massflx=0.
!	mconv=0.
!	umean=0.
!	vmean=0.
!	pmean=0.
!	direction=0.
!	AA0=0.
!	PRET=0.
!	PSUR=0.
!	TER11=0.
!       KDT=0
!       iact_gr=0
!       iact_old_gr=0
!       KDET=0
!       tkeg=0.
!       rcpg=0.
!       glatg=0.
!       glong=0.
! END SUBROUTINE zero_scratch2_grell

END MODULE mem_scratch2_grell
