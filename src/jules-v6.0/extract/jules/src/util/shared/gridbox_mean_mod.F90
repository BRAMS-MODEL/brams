MODULE gridbox_mean_mod

! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

PRIVATE
PUBLIC surftiles_to_gbm, soiltiles_to_gbm, pfttiles_to_gbm

CONTAINS

FUNCTION surftiles_to_gbm(tile_data, ainfo, tile_mask) RESULT(gbm_data)

USE ancil_info, ONLY: land_pts, nsurft, surft_pts,                            &
                      ainfo_type

USE jules_surface_mod, ONLY: l_aggregate

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes data on land points * nsurft and calculates a gridbox mean value
!   for each land point. If mask is given, only the tiles for which mask
!   is true are included
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
REAL, INTENT(IN) :: tile_data(land_pts,nsurft)  ! Per tile data

TYPE(ainfo_type), INTENT(IN) :: ainfo

LOGICAL, INTENT(IN), OPTIONAL :: tile_mask(nsurft)
                                  ! T - include tile in calculation of gbm
                                  ! F - do not include tile in calculation

! Return type
REAL :: gbm_data(land_pts)

! Work variables
INTEGER :: i, p, t  ! Index variables

LOGICAL :: tile_mask_local(nsurft)  ! Local version of tile_mask that is
                                    ! always present

!-------------------------------------------------------------------------------
! Set mask.
IF (PRESENT( tile_mask )) THEN
  tile_mask_local(:) = tile_mask(:)
ELSE
  tile_mask_local(:) = .TRUE.
END IF

! Initialise the average.
gbm_data(:) = 0.0

IF ( l_aggregate ) THEN
  ! If l_aggregate is .TRUE., then all tile variables are essentially gridbox
  ! means already
  gbm_data(:) = tile_data(:,1)
ELSE
  ! Otherwise, we can just use ainfo%frac_surft, since nsurft=ntype
  DO t = 1,nsurft
    IF ( tile_mask_local(t) ) THEN
      DO i = 1,surft_pts(t)
        p = ainfo%surft_index(i,t)
        gbm_data(p) = gbm_data(p) + ainfo%frac_surft(p,t) * tile_data(p,t)
      END DO
    END IF
  END DO
END IF

RETURN

END FUNCTION surftiles_to_gbm

!-------------------------------------------------------------------------------

FUNCTION soiltiles_to_gbm(tile_data, ainfo, tile_mask) RESULT(gbm_data)

USE ancil_info, ONLY: land_pts, nsoilt, soilt_pts,                            &
                      ainfo_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes data on land points * nsoilt and calculates a gridbox mean value
!   for each land point. If mask is given, only the tiles for which mask
!   is true are included
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
REAL, INTENT(IN) :: tile_data(land_pts,nsoilt)  ! Per tile data

TYPE(ainfo_type), INTENT(IN) :: ainfo

LOGICAL, INTENT(IN), OPTIONAL :: tile_mask(nsoilt)
                                  ! T - include tile in calculation of gbm
                                  ! F - do not include tile in calculation

! Return type
REAL :: gbm_data(land_pts)

! Work variables
INTEGER :: i, p, t  ! Index variables

LOGICAL :: tile_mask_local(nsoilt)  ! Local version of tile_mask that is
                                    ! always present

!-------------------------------------------------------------------------------
! Set mask.
IF (PRESENT( tile_mask )) THEN
  tile_mask_local(:) = tile_mask(:)
ELSE
  tile_mask_local(:) = .TRUE.
END IF

! Initialise the average
gbm_data(:) = 0.0

! 1 soil tile is a common use case, so allow for a straight copy
IF ( nsoilt == 1) THEN
  gbm_data(:) = tile_data(:,1)
ELSE
  ! Otherwise add up the contributions
  DO t = 1,nsoilt
    IF ( tile_mask_local(t) ) THEN
      DO i = 1,soilt_pts(t)
        p = ainfo%soilt_index(i,t)
        gbm_data(p) = gbm_data(p) + (ainfo%frac_soilt(p,t) * tile_data(p,t))
      END DO
    END IF
  END DO
END IF

RETURN
END FUNCTION soiltiles_to_gbm

!-----------------------------------------------------------------------------

FUNCTION pfttiles_to_gbm(tile_data, ainfo, tile_mask) RESULT(gbm_data)

USE ancil_info, ONLY: land_pts, nsurft, surft_pts,                            &
                      ainfo_type

USE jules_surface_mod, ONLY: l_aggregate
USE jules_surface_types_mod,      ONLY: npft

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes data on land points * npft and calculates a gridbox mean value
!   for each land point. If mask is given, only the tiles for which mask
!   is true are included
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
REAL, INTENT(IN) :: tile_data(land_pts,npft)  ! Per tile data

TYPE(ainfo_type), INTENT(IN) :: ainfo

LOGICAL, INTENT(IN), OPTIONAL :: tile_mask(npft)
                                  ! T - include tile in calculation of gbm
                                  ! F - do not include tile in calculation

! Return type
REAL :: gbm_data(land_pts)

! Work variables
INTEGER :: i, p, t  ! Index variables

LOGICAL :: tile_mask_local(npft)  ! Local version of tile_mask that is
                                    ! always present

!-------------------------------------------------------------------------------
! Set mask.
IF (PRESENT( tile_mask )) THEN
  tile_mask_local(:) = tile_mask(:)
ELSE
  tile_mask_local(:) = .TRUE.
END IF

! Initialise the average.
gbm_data(:) = 0.0

IF ( l_aggregate ) THEN
  ! If l_aggregate is .TRUE., then all tile variables are essentially gridbox
  ! means already
  gbm_data(:) = tile_data(:,1)
ELSE
  ! Otherwise, we can just use ainfo%frac_surft, since npft is the first part of nsurft
  DO t = 1,npft
    IF ( tile_mask_local(t) ) THEN
      DO i = 1,surft_pts(t)
        p = ainfo%surft_index(i,t)
        gbm_data(p) = gbm_data(p) + ainfo%frac_surft(p,t) * tile_data(p,t)
      END DO
    END IF
  END DO
END IF

RETURN

END FUNCTION pfttiles_to_gbm

END MODULE gridbox_mean_mod
