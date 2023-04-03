! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in Veg3 Ecosystem Demography
! *****************************COPYRIGHT****************************************

MODULE veg3_parm_mod

IMPLICIT NONE

!Set up object for veg3 control

TYPE :: veg3_ctrl_type
  INTEGER :: land_pts,nsurft,npft,nnpft,soil,triffid_period,nstep_trif,nmasst
  REAL    :: timestep,frac_min
END TYPE veg3_ctrl_type

!Set up objects containing everything we need for litter calculation

TYPE :: litter_parm_type
  REAL, ALLOCATABLE ::                                                        &
    g_wood(:),              & ! Turnover rate for woody biomass (/s).
    g_root(:),              & ! Turnover rate for root biomass (/s).
    g_leaf(:)                 ! Turnover rate for leaf biomass (/s).
END TYPE litter_parm_type

TYPE :: red_parm_type
  INTEGER, ALLOCATABLE ::                                                     &
    veg_type(:),                                                              &
    mclass(:)
  REAL, ALLOCATABLE ::                                                        &
    phi_g(:),                                                                 &
    phi_h(:),                                                                 &
    phi_a(:),                                                                 &
    alpha(:),                                                                 &
    frac_min(:),                                                              &
    mass0(:),                                                                 &
    height0(:),                                                               &
    crwn_area0(:),                                                            &
    mort_base(:),                                                             &
    mclass_geom_mult(:) 
END TYPE red_parm_type

TYPE(veg3_ctrl_type)   :: veg3_ctrl
TYPE(litter_parm_type) :: litter_parms
TYPE(red_parm_type)    :: red_parms

! Switches
LOGICAL :: l_veg3 = .FALSE.

!Private by default
PRIVATE

!Expose routines
PUBLIC :: veg3_parm_init,veg3_parm_allocate

!Expose data
PUBLIC :: veg3_ctrl, litter_parms, red_parms, l_veg3

!Expose data structures
PUBLIC :: veg3_ctrl_type, litter_parm_type, red_parm_type

!Allow external code to read but not write
PROTECTED :: l_veg3, litter_parms, veg3_ctrl, red_parms

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='VEG3_PARM_MOD'

CONTAINS
!-------------------------------------------------------------------------------

SUBROUTINE veg3_parm_allocate(land_pts,nsurft,nnpft,npft)

IMPLICIT NONE
INTEGER, INTENT(IN) :: land_pts, nsurft, nnpft, npft

!End of Header

! Set up litter_parms object

ALLOCATE(litter_parms%g_wood(nnpft))
ALLOCATE(litter_parms%g_leaf(nnpft))
ALLOCATE(litter_parms%g_root(nnpft))

! Allocate red_data_type
ALLOCATE(red_parms%veg_type           (nnpft))
ALLOCATE(red_parms%mclass             (nnpft))
ALLOCATE(red_parms%phi_g              (nnpft))
ALLOCATE(red_parms%phi_h              (nnpft))
ALLOCATE(red_parms%phi_a              (nnpft))
ALLOCATE(red_parms%alpha              (nnpft))
ALLOCATE(red_parms%frac_min           (nnpft))
ALLOCATE(red_parms%mass0              (nnpft))
ALLOCATE(red_parms%height0            (nnpft))
ALLOCATE(red_parms%crwn_area0         (nnpft))
ALLOCATE(red_parms%mort_base          (nnpft))
ALLOCATE(red_parms%mclass_geom_mult   (nnpft))

RETURN
END SUBROUTINE veg3_parm_allocate

!-------------------------------------------------------------------------------
SUBROUTINE veg3_set_parms(land_pts,nsurft,nnpft,npft,nmasst)

!Source parms, etc from io modules - these should be available on and offline

USE pftparm,                  ONLY: g_leaf_0
USE trif,                     ONLY: g_root, g_wood
! Above only allocated if triffid on - needs to be addressed
USE jules_vegetation_mod,     ONLY: l_triffid, triffid_period,frac_min

!Get the timestep length
#if ! defined(UM_JULES)
USE model_time_mod,           ONLY: timestep_len => timestep
#endif
USE timestep_mod,             ONLY: timestep
USE conversions_mod,          ONLY: rsec_per_day

USE jules_surface_types_mod,  ONLY: soil

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, nsurft, nnpft, npft,nmasst

INTEGER :: l,n

!End of header

IF (l_veg3 .AND. l_triffid) THEN

  ! Allocate memory for veg3 objects
  !  CALL veg3_parm_allocate(land_pts,nsurft,nnpft,npft)

  ! Setup parameters/control

  !Get the timestep length
  !Convert JULES to real number - important for call later
  !This way call occurs on same timestep on or offline - addresses issue in
  !current setup such that REAL is always used

  veg3_ctrl%timestep = REAL(timestep)
  veg3_ctrl%triffid_period = triffid_period
  veg3_ctrl%nstep_trif = INT(rsec_per_day * veg3_ctrl%triffid_period          &
    / veg3_ctrl%timestep)

  veg3_ctrl%land_pts = land_pts
  veg3_ctrl%nsurft   = nsurft
  veg3_ctrl%nnpft    = nnpft
  veg3_ctrl%npft     = npft
  veg3_ctrl%nmasst   = nmasst
  veg3_ctrl%soil     = soil
  veg3_ctrl%frac_min = frac_min

  ! Convert to s-1 from 360d-1
  litter_parms%g_wood(:) = g_wood(1:nnpft)   / rsec_per_day / 360.0
  litter_parms%g_leaf(:) = g_leaf_0(1:nnpft) / rsec_per_day / 360.0
  litter_parms%g_root(:) = g_root(1:nnpft)   / rsec_per_day / 360.0

  !RED 
  red_parms%veg_type(:)         = 0
  red_parms%mclass(:)           = 0
  red_parms%phi_g(:)            = 0.0
  red_parms%phi_h(:)            = 0.0
  red_parms%phi_a(:)            = 0.0
  red_parms%alpha(:)            = 0.0
  red_parms%frac_min(:)         = 0.0
  red_parms%mass0(:)            = 0.0
  red_parms%height0(:)          = 0.0
  red_parms%crwn_area0(:)       = 0.0
  red_parms%mort_base(:)        = 0.0
  red_parms%mclass_geom_mult(:) = 0.0
  red_parms%frac_min(:)         = frac_min

END IF

RETURN
END SUBROUTINE veg3_set_parms

!-------------------------------------------------------------------------------

SUBROUTINE veg3_parm_init(land_pts,nsurft,nnpft,npft,nmasst)

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, nnpft, npft, nsurft, nmasst

!End of header

CALL veg3_set_parms(land_pts,nsurft,nnpft,npft,nmasst)

RETURN
END SUBROUTINE veg3_parm_init


END MODULE veg3_parm_mod
