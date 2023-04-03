! *****************************COPYRIGHT*******************************
! (c) Crown copyright, Met Office, All Rights Reserved.
! Please refer to file $UMDIR/vn$VN/copyright.txt for further details
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE elevate_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ELEVATE_MOD'

CONTAINS

!     SUBROUTINE ELEVATE ------------------------------------------

!     Purpose:
!     Calculate temperature and humidity at a given elevation above the
!     mean gridbox surface

!     ------------------------------------------------------------------
SUBROUTINE elevate (                                                          &
 land_pts,nsurft,surft_pts,land_index,surft_index,                            &
 tl_1,qw_1,pstar,surf_hgt,l_elev_absolute_height,z_land,                      &
 t_elev,q_elev)

USE atm_fields_bounds_mod, ONLY: tdims
USE theta_field_sizes, ONLY: t_i_length, t_j_length

USE planet_constants_mod, ONLY:                                               &
  cp, c_virtual,                                                              &
  grcp,                                                                       &
    ! Dry adiabatic lapse rate
  g, r, repsilon
USE dewpnt_mod, ONLY: dewpnt

USE water_constants_mod, ONLY: lc
USE qsat_mod, ONLY: qsat

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                        &
 land_pts                                                                     &
                       ! IN No of land points being processed.
,nsurft                                                                       &
                       ! IN Number of land tiles per land point.
,land_index(land_pts)                                                         &
                       ! IN Index of land points.
,surft_index(land_pts,nsurft)                                                 &
                       ! IN Index of tile points.
,surft_pts(nsurft)      ! IN Number of tile points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
 tl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                       ! IN Liquid/frozen water temperature for
!                            !    lowest atmospheric layer (K).
,qw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                       ! IN Total water content of lowest
!                            !    atmospheric layer (kg per kg air).
,pstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                       ! IN Surface pressure (Pascals).
,z_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                       ! IN Gridbox mean height
,surf_hgt(land_pts,nsurft)
                       ! IN Height of elevated tile above
!                            !        mean gridbox surface (m)

LOGICAL, INTENT(IN) ::                                                        &
 l_elev_absolute_height(nsurft)
                       ! IN switch for whether surf_hgt is an offset
                       ! to the gridbox mean or an absolute height

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
 t_elev(land_pts,nsurft)                                                      &
                          ! OUT Temperature at elevated height (k)
,q_elev(land_pts,nsurft)  ! OUT Specific humidity at elevated
!                               !     height (kg per kg air)


! Local variables
REAL(KIND=real_jlslsm) ::                                                     &
 tdew(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                        ! Dew point temperature for input
!                             ! specific humidity (K)
,tv                                                                           &
                        ! Virtual temeprature for TDEW
,salr                                                                         &
                        ! Saturated adiabatic lapse rate for
!                             ! input specific humidity
,z                                                                            &
                        ! Height at which air becomes saturated
,delevation
                        ! Height to adjust climate by

! Scalars
INTEGER ::                                                                    &
 i,j                                                                          &
                     ! Horizontal field index.
,k                                                                            &
                     ! Tile field index.
,l                                                                            &
                     ! Land point field index.
,n                   ! Tile index loop counter

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ELEVATE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate the dew point temperature
CALL dewpnt(qw_1,pstar,tl_1,t_i_length * t_j_length,tdew)

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(n,l,i,j,k,delevation,tv,salr,z)          &
!$OMP SHARED(nsurft,land_pts,t_elev,q_elev,surft_pts,surft_index,land_index,  &
!$OMP        t_i_length,l_elev_absolute_height,surf_hgt,z_land,tl_1,grcp,     &
!$OMP        tdew,qw_1,c_virtual,r,g,repsilon,cp,pstar)

! Initialise output arrays (for clarity and safety).
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
DO n = 1,nsurft
  DO l = 1,land_pts
    t_elev(l,n) = 0.0
    q_elev(l,n) = 0.0
  END DO
END DO
!$OMP END DO

DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
  DO k = 1,surft_pts(n)
    l = surft_index(k,n)
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length

    IF (l_elev_absolute_height(n)) THEN
      delevation = surf_hgt(l,n) - z_land(i,j)
    ELSE
      delevation = surf_hgt(l,n)
    END IF

    t_elev(l,n) = tl_1(i,j) - delevation * grcp

    IF (t_elev(l,n) <  tdew(i,j)) THEN
      ! Temperature following a dry adiabat is less than dew point temperature
      ! Therefore need to follow a saturated adiabate from height of
      ! dew point temperature

      tv = tdew(i,j) * (1.0 + c_virtual * qw_1(i,j))
      salr = g  * (1.0 + lc * qw_1(i,j) / (r * tv * (1.0 - qw_1(i,j)))) /     &
             ( cp + lc**2.0 * qw_1(i,j) * repsilon /                          &
                      (r * tv**2.0 * (1.0 - qw_1(i,j))))

      z = (tl_1(i,j) - tdew(i,j)) / grcp

      t_elev(l,n) = tdew(i,j) - (delevation - z) * salr

      CALL qsat(q_elev(l,n),t_elev(l,n),pstar(i,j))

    ELSE
      ! Temperature follows a dry adiabatic lapse rate and humidity remains constant

      q_elev(l,n) = qw_1(i,j)

    END IF

  END DO
!$OMP END DO NOWAIT
END DO

!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE elevate

END MODULE elevate_mod
