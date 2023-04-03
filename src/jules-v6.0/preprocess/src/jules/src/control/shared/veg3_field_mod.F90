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

MODULE veg3_field_mod

USE veg3_parm_mod, ONLY: red_parms, l_veg3

USE um_types, ONLY: real_jlslsm

!Use at module level
USE ancil_info,    ONLY: ainfo_type
USE prognostics,    ONLY: progs_type

IMPLICIT NONE

TYPE :: veg_state_type
  REAL, ALLOCATABLE ::                                                        &
      leafC(:,:),                                                             &
      rootC(:,:),                                                             &
      woodC(:,:),                                                             &
      lai_bal(:,:),                                                           &
      frac(:,:),                                                              &
      vegCpft(:,:),                                                           &
      vegC(:),                                                                &
      npp_acc(:,:),                                                           &
      canht(:,:)            ! Copy of prognostic
END TYPE veg_state_type

TYPE :: red_state_type
  REAL, ALLOCATABLE ::                                                        &
    mass_mass(:,:),                                                           &
    ht_mass(:,:),                                                             &
    crwn_area_mass(:,:),                                                      &
    g_mass_scale(:,:),                                                        &
    mclass_geom_mult(:),                                                      &
    plantNumDensity(:,:,:),                              & ! Copy of prognostic
    mort(:,:,:)                 
END TYPE red_state_type

TYPE(veg_state_type)   :: veg_state
TYPE(red_state_type)   :: red_state

!Private by default
PRIVATE

!Expose routines
PUBLIC :: veg3_field_init, veg3_field_allocate

!Expose data
PUBLIC :: veg_state, red_state

!Expose data structures
PUBLIC :: veg_state_type,red_state_type

!Allow external code to read but not write
!PROTECTED :: 

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='VEG3_FIELD_MOD'

CONTAINS
!-------------------------------------------------------------------------------

SUBROUTINE veg3_field_allocate(land_pts,nsurft,nnpft,npft,nmasst)

IMPLICIT NONE
INTEGER, INTENT(IN) :: land_pts, nsurft, nnpft, npft, nmasst

!End of Header

ALLOCATE(veg_state%leafC        ( land_pts, nnpft) )
ALLOCATE(veg_state%rootC        ( land_pts, nnpft) )
ALLOCATE(veg_state%woodC        ( land_pts, nnpft) )
ALLOCATE(veg_state%vegCpft      ( land_pts, nnpft) )
ALLOCATE(veg_state%lai_bal      ( land_pts, nnpft) )
ALLOCATE(veg_state%canht        ( land_pts, nnpft) )
ALLOCATE(veg_state%npp_acc      ( land_pts, nnpft) )

ALLOCATE(veg_state%frac         ( land_pts, nsurft) )
ALLOCATE(veg_state%vegC         ( land_pts) )

!Initialise
veg_state%leafC(:,:)          = 0.0
veg_state%rootC(:,:)          = 0.0
veg_state%woodC(:,:)          = 0.0
veg_state%vegCpft(:,:)        = 0.0
veg_state%lai_bal(:,:)        = 0.0
veg_state%canht(:,:)          = 0.0
veg_state%frac(:,:)           = 0.0
veg_state%vegC(:)             = 0.0
veg_state%npp_acc(:,:)        = 0.0


! RED

! Allocate red_data_type
ALLOCATE(red_state%mass_mass          (nnpft, nmasst ))
ALLOCATE(red_state%ht_mass            (nnpft, nmasst ))
ALLOCATE(red_state%crwn_area_mass     (nnpft, nmasst ))
ALLOCATE(red_state%g_mass_scale       (nnpft, nmasst ))
ALLOCATE(red_state%plantNumDensity    (land_pts, nnpft, nmasst ))
ALLOCATE(red_state%mort               (land_pts, nnpft, nmasst ))

! Initialise red_data_type

red_state%mass_mass(:,:)          = 0.0
red_state%ht_mass(:,:)            = 0.0
red_state%crwn_area_mass(:,:)     = 0.0
red_state%g_mass_scale(:,:)       = 0.0
red_state%plantNumDensity(:,:,:)  = 0.0
red_state%mort(:,:,:)             = 0.0

RETURN
END SUBROUTINE veg3_field_allocate

!-------------------------------------------------------------------------------
SUBROUTINE veg3_set_fields(land_pts,nsurft,nnpft,npft,ainfo,progs)

!Source parms, etc from io modules - these should be available on and offline

USE pftparm,                  ONLY: g_leaf_0
USE trif,                     ONLY: g_root, g_wood
! Above only allocated if triffid on - needs to be addressed
USE jules_vegetation_mod,     ONLY: l_triffid, triffid_period,frac_min

!Get the timestep length
USE model_time_mod,           ONLY: timestep_len
USE timestep_mod,             ONLY: timestep
USE conversions_mod,          ONLY: rsec_per_day

USE jules_surface_types_mod,  ONLY: soil

!Functions
USE calc_c_comps_triffid_mod, ONLY: calc_c_comps_triffid
USE gridbox_mean_mod,         ONLY: pfttiles_to_gbm

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, nsurft, nnpft, npft

TYPE(ainfo_type), INTENT(IN) :: ainfo
TYPE(progs_type), INTENT(IN) :: progs

INTEGER :: l,n,k

!End of header

IF (l_veg3 .AND. l_triffid) THEN

  ! Allocate memory for veg3 objects
  !  CALL veg3_allocate(land_pts,nsurft,nnpft,npft)

  ! Set up fields
  DO n = 1, nnpft
    DO l = 1, land_pts
      ! Take a copy of the prognostic to store in object
      veg_state%canht(l,n) = progs%canht_pft(l,n)

      CALL calc_c_comps_triffid(n,                                            &
                    veg_state%canht(l,n),         &       !In
                    veg_state%lai_bal(l,n),       &       !Out
                    veg_state%leafC(l,n),         &       !Out
                    veg_state%rootC(l,n),         &       !Out
                    veg_state%woodC(l,n),         &       !Out
                    veg_state%vegCpft(l,n))               !Out

    END DO
  END DO

  DO n = 1, nsurft
    DO l = 1, land_pts
      ! Now do fractions
      veg_state%frac(l,n) = ainfo%frac_surft(l,n)
    END DO
  END DO

  ! Aggregate vegC to gridbox for use as diagnostic
  veg_state%vegC = pfttiles_to_gbm(veg_state%vegCpft, ainfo)

END IF

RETURN
END SUBROUTINE veg3_set_fields

!-------------------------------------------------------------------------------
SUBROUTINE veg3_red_set_fields(land_pts,nsurft,nnpft,npft,nmasst)

!Source parms, etc from io modules - these should be available on and offline

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, nsurft, nnpft, npft, nmasst

INTEGER :: l,n,k

!End of header

! Setup Allometry
DO k = 1,nmasst
  DO n = 1,nnpft
    IF (k == 1) THEN
      red_state%mass_mass(n,k) = red_parms%mass0(n)
      red_state%ht_mass(n,k) = red_parms%height0(n)
      red_state%crwn_area_mass(n,k) = red_parms%crwn_area0(n)
      red_state%g_mass_scale(n,k) = 1.0
    ELSE IF (k > 1 .AND. k <= red_parms%mclass(n)) THEN
      red_state%mass_mass(n,k) = red_state%mass_mass(n,k-1) * red_parms%mclass_geom_mult(n)
      red_state%ht_mass(n,k) = red_parms%height0(n) * (red_state%mass_mass(n,k) &
        / red_parms%mass0(n))** red_parms%phi_h(n)
      red_state%crwn_area_mass(n,k) = red_parms%crwn_area0(n) * (red_state%mass_mass(n,k) &
        / red_parms%mass0(n))** red_parms%phi_a(n)
      red_state%g_mass_scale(n,k) = (red_state%mass_mass(n,k)                 &
        / red_parms%mass0(n))** red_parms%phi_g(n)
    END IF
  END DO
END DO

RETURN
END SUBROUTINE veg3_red_set_fields

!-------------------------------------------------------------------------------

SUBROUTINE veg3_field_init(land_pts,nsurft,nnpft,npft,nmasst,ainfo,progs)

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, nnpft, npft, nsurft,nmasst

TYPE(ainfo_type), INTENT(IN) :: ainfo
TYPE(progs_type), INTENT(IN) :: progs

!End of header

CALL veg3_set_fields(land_pts,nsurft,nnpft,npft,ainfo,progs)
CALL veg3_red_set_fields(land_pts,nsurft,nnpft,npft,nmasst)

RETURN
END SUBROUTINE veg3_field_init


END MODULE veg3_field_mod
