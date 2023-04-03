! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE sf_resist_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SF_RESIST_MOD'

CONTAINS
!    L  SUBROUTINE SF_RESIST----------------------------------------------
!
!    Purpose: Calculate surface moisture flux resistance factors.
!
!    Arguments --------------------------------------------------------
SUBROUTINE sf_resist (                                                        &
 land_pts,surft_pts,land_index,surft_index,                                   &
 canopy,catch,ch,dq,epdt,flake,gc,gc_stom_surft,snow_surft,vshr,              &
 fraca,resfs,resft,gc_irr,resfs_irr,resfs_stom,l_et_stom,l_et_stom_surft)

USE atm_fields_bounds_mod, ONLY: tdims
USE theta_field_sizes, ONLY: t_i_length

USE jules_irrig_mod, ONLY: l_irrig_dmd

USE jules_snow_mod, ONLY: maskd, frac_snow_subl_melt
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                        &
 land_pts                                                                     &
                       ! IN Total number of land points.
,surft_pts                                                                    &
                       ! IN Number of tile points.
,land_index(land_pts)                                                         &
                       ! IN Index of land points.
,surft_index(land_pts) ! IN Index of tile points.

LOGICAL, INTENT(IN) ::                                                        &
 l_et_stom                                                                    &
                        !IN flag for calculating transpiration for gridbox
,l_et_stom_surft
                        !IN flag for calculating transpiration on tiles


REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
 canopy(land_pts)                                                             &
                     ! IN Surface water (kg per sq metre).  F642.
,catch(land_pts)                                                              &
                     ! IN Surface capacity (max. surface water)
!                          !    (kg per sq metre).  F6416.
,ch(land_pts)                                                                 &
                     ! IN Transport coefficient for heat and
!                          !    moisture transport
,dq(land_pts)                                                                 &
                     ! IN Sp humidity difference between surface
!                          !    and lowest atmospheric level (Q1 - Q*).
,epdt(land_pts)                                                               &
                     ! IN "Potential" Evaporation * Timestep.
!                          !    Dummy variable for first call to routine
,flake(land_pts)                                                              &
                     ! IN Lake fraction.
,gc(land_pts)                                                                 &
                     ! IN Interactive canopy conductance
!                          !    to evaporation (m/s)
,gc_stom_surft(land_pts)                                                      &
                     ! IN canopy conductance (excluding soil) (m/s)
,snow_surft(land_pts)                                                         &
                     ! IN Lying snow amount (kg per sq metre).
,vshr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                     ! IN Magnitude of surface-to-lowest-level
!                          !     windshear

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
 fraca(land_pts)                                                              &
                     ! OUT Fraction of surface moisture flux with
!                          !     only aerodynamic resistance.
,resfs(land_pts)                                                              &
                     ! OUT Combined soil, stomatal and aerodynamic
!                          !     resistance factor for fraction 1-FRACA.
,resft(land_pts)                                                              &
                     ! OUT Total resistance factor
!                          !     FRACA+(1-FRACA)*RESFS.
,resfs_stom(land_pts)
                     ! OUT Combined stomatal and aerodynamic
!                          !     resistance factor for fraction 1-FRACA.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
! Needed for irrigation
 gc_irr(land_pts)
!                    ! IN Interactive canopy conductance
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
 resfs_irr(land_pts)
!                    ! OUT Combined soil, stomatal and aerodynamic
!                          !     resistance factor for fraction 1-FRACA.

! Workspace -----------------------------------------------------------
INTEGER ::                                                                    &
 i,j                                                                          &
             ! Horizontal field index.
,k                                                                            &
             ! Tile field index.
,l           ! Land field index.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SF_RESIST'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
!     Evaporation over land surfaces without snow is limited by
!     soil moisture availability and stomatal resistance.
!     Set FRACA (= fA in the documentation) according to P243.68,
!     and RESFS (= fS) according to P243.75 and P243.61.
!-----------------------------------------------------------------------

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(l,k,j,i)                                                        &
!$OMP SHARED(surft_pts,surft_index,land_index,t_i_length,fraca,dq,snow_surft, &
!$OMP        catch,frac_snow_subl_melt,maskd,resfs,gc,ch,vshr,l_et_stom,      &
!$OMP        l_et_stom_surft,resfs_stom,gc_stom_surft,flake,l_irrig_dmd,      &
!$OMP        resfs_irr,gc_irr,canopy,epdt,resft)
DO k = 1,surft_pts
  l = surft_index(k)
  j=(land_index(l) - 1) / t_i_length + 1
  i = land_index(l) - (j-1) * t_i_length

  !-----------------------------------------------------------------------
  ! Calculate the fraction of the flux with only aerodynamic resistance
  ! (canopy evaporation).
  ! Set to 1 for negative moisture flux or snow-covered land
  ! (no surface/stomatal resistance to condensation).
  !-----------------------------------------------------------------------
  fraca(l) = 1.0
  IF (dq(l) <  0.0 .AND. snow_surft(l) <= 0.0) fraca(l) = 0.0
  IF (dq(l) <  0.0 .AND. snow_surft(l) <= 0.0 .AND. catch(l) >  0.0)          &
    fraca(l) = canopy(l) / ( epdt(l) + catch(l) )
  IF (snow_surft(l) > 0.0 .AND. frac_snow_subl_melt == 1) THEN
    fraca(l) = 1.0 - EXP(-maskd * snow_surft(l))
  END IF
  fraca(l) = MIN(fraca(l),1.0)

  !-----------------------------------------------------------------------
  ! Calculate resistance factors for transpiration from vegetation tiles
  ! and bare soil evaporation from soil tiles.
  !-----------------------------------------------------------------------
  resfs(l) = gc(l) / ( gc(l) + ch(l) * vshr(i,j) )
  IF (l_et_stom .OR. l_et_stom_surft) THEN
    resfs_stom(l) = gc_stom_surft(l) / ( gc_stom_surft(l) + ch(l) * vshr(i,j) )
  END IF
  resft(l) = flake(l) + (1.0 - flake(l)) *                                    &
                        ( fraca(l) + (1.0 - fraca(l)) * resfs(l) )
  IF (l_irrig_dmd) THEN
    resfs_irr(l) = gc_irr(l) / ( gc_irr(l) + ch(l) * vshr(i,j) )
  END IF

END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE sf_resist
END MODULE sf_resist_mod
