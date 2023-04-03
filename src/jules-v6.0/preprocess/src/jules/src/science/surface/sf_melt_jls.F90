! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE sf_melt_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SF_MELT_MOD'

CONTAINS
! SUBROUTINE SF_MELT----------------------------------------------------
!
! Purpose : Calculates surface melting (snow and sea-ice) and increments
!           surface fluxes to satisfy energy balance.
!           Sub-surface snowmelt is calculated and snowdepth incremented
!           by melt and sublimation in P251.
!-----------------------------------------------------------------------
SUBROUTINE sf_melt (                                                          &
 points,pts_index                                                             &
,surft_index,surft_pts,fld_sea                                                &
,alpha1,ashtf_prime,dtrdz_1                                                   &
,resft,rhokh_1,tile_frac,timestep,r_gamma                                     &
,ei_surft,fqw_1,ftl_1,fqw_surft,ftl_surft                                     &
,tstar_surft,snow_surft,snowdepth                                             &
,melt_surft                                                                   &
 )

USE atm_fields_bounds_mod, ONLY: tdims
USE theta_field_sizes, ONLY: t_i_length

USE planet_constants_mod, ONLY: cp

USE jules_snow_mod, ONLY: maskd, rho_snow_const,                              &
                          frac_snow_subl_melt, l_snowdep_surf

USE water_constants_mod, ONLY:                                                &
 lc, lf, rho_water, tm

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER ::                                                                    &
 points                                                                       &
                      ! IN Total number of points.
,pts_index(points)                                                            &
                      ! IN Index of points.
,surft_index(points)                                                          &
!                           ! IN Index of tile points.
,surft_pts             ! IN Number of tile points.

REAL(KIND=real_jlslsm) ::                                                     &
fld_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
!                           ! IN Fraction of land or sea.
,alpha1(points)                                                               &
!                           ! IN Gradients of saturated specific
!                           !    humidity with respect to temp.
!                           !    between the bottom model layer
!                           !    and surface.
,ashtf_prime(points)                                                          &
!                           ! IN Adjusted SEB coefficient
,dtrdz_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
!                           ! IN -g.dt/dp for surface layer
,resft(points)                                                                &
                       !IN Resistance factor.
,rhokh_1(points)                                                              &
!                           ! IN Surface exchange coefficient.
,tile_frac(points)                                                            &
!                           ! IN Tile fractions.
,timestep                                                                     &
                      ! IN Timestep (sec).
,r_gamma              ! IN implicit weight in level 1

REAL(KIND=real_jlslsm) ::                                                     &
 ei_surft(points)                                                             &
!                           ! INOUT Sublimation for tile (kg/m2/s)
,fqw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
!                           ! INOUT GBM surface moisture flux (kg/m2/s).
,ftl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
!                           ! INOUT GBM surface sens. heat flux (W/m2).
,fqw_surft(points)                                                            &
!                           ! INOUT FQW for tile.
,ftl_surft(points)                                                            &
!                           ! INOUT FTL for tile.
,tstar_surft(points)                                                          &
!                           ! INOUT Tile surface temperatures (K).
,snow_surft(points)                                                           &
!                           ! INOUT Lying snow on tile (kg/m2).
,snowdepth(points)
!                           ! INOUT Depth of snow on tile (m).

REAL(KIND=real_jlslsm) ::                                                     &
 melt_surft(points)
!                           ! OUT Surface snowmelt on tiles (kg/m2/s).

REAL(KIND=real_jlslsm) ::                                                     &
 dfqw                                                                         &
                      ! Moisture flux increment.
,dftl                                                                         &
                      ! Sensible heat flux increment.
,dtstar                                                                       &
                      ! Surface temperature increment.
,lcmelt                                                                       &
                      ! Temporary in melt calculations.
,lsmelt                                                                       &
                      ! Temporary in melt calculations.
,rhokh1_prime                                                                 &
                      ! Modified forward time-weighted
!                           ! transfer coefficient.
,snow_max                                                                     &
                      ! Snow available for melting.
,snow_density
                      ! Density of snow on input
INTEGER ::                                                                    &
 i,j                                                                          &
                      ! Loop counter - full horizontal field.
,k                                                                            &
                      ! Loop counter - tile field
,l
                      ! Loop counter - land field.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SF_MELT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

melt_surft(:) = 0.0

!-----------------------------------------------------------------------
!  Melt snow on tile if TSTAR_SURFT is greater than TM.
!-----------------------------------------------------------------------

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(SHARED)                                                         &
!$OMP PRIVATE(k,j,i,l,snow_density,snow_max,rhokh1_prime,lcmelt,lsmelt,dtstar,&
!$OMP         dftl,dfqw)
DO k = 1,surft_pts
  l = surft_index(k)
  j=(pts_index(l) - 1) / t_i_length + 1
  i = pts_index(l) - (j-1) * t_i_length
  !
  IF (snowdepth(l) > SQRT(TINY(1.0))) THEN
    snow_density = snow_surft(l) / snowdepth(l)
  ELSE
    snow_density = rho_snow_const
  END IF
  !
  ! The l_snowdep_surf if-test below is temporary to preserve
  ! bit comparison. It should be removed when this switch becomes
  ! default .TRUE.
  !
  IF (l_snowdep_surf) THEN
    snow_density = MAX(rho_snow_const,snow_density)
    snow_density = MIN(rho_water     ,snow_density)
  ELSE
    snow_density = MAX(   1.0,snow_density)
    snow_density = MIN(1000.0,snow_density)
  END IF
  !
  snow_max = MAX( 0.0, snow_surft(l) - ei_surft(l) * timestep )
  IF ( snow_max >  0.0 .AND. tstar_surft(l) >  tm ) THEN
    rhokh1_prime = 1.0 / ( 1.0 / rhokh_1(l)                                   &
                       + r_gamma * dtrdz_1(i,j) )
    lcmelt = (cp + lc * alpha1(l) * resft(l)) * rhokh1_prime                  &
             + ashtf_prime(l)
    lsmelt = lcmelt + lf * alpha1(l) * rhokh1_prime
    IF (frac_snow_subl_melt == 1) THEN
      dtstar = - MIN( (tstar_surft(l) - tm) *                                 &
               (1.0 - EXP(-maskd * snow_max / snow_density)),                 &
               lf * snow_max / (lcmelt * timestep) )
    ELSE
      dtstar = - MIN( tstar_surft(l) - tm ,                                   &
                      lf * snow_max / (lcmelt * timestep) )
    END IF
    tstar_surft(l) = tstar_surft(l) + dtstar
    melt_surft(l) = - lsmelt * dtstar / lf
    dftl = cp * rhokh1_prime * dtstar
    dfqw = alpha1(l) * resft(l) * rhokh1_prime * dtstar
    ftl_surft(l) = ftl_surft(l) + dftl
    fqw_surft(l) = fqw_surft(l) + dfqw
    ei_surft(l) = ei_surft(l) + dfqw
    !-----------------------------------------------------------------------
    !  Update gridbox-mean quantities
    !-----------------------------------------------------------------------
    dftl = tile_frac(l) * dftl
    dfqw = tile_frac(l) * dfqw
    ftl_1(i,j) = ftl_1(i,j) + fld_sea(i,j) * dftl
    fqw_1(i,j) = fqw_1(i,j) + fld_sea(i,j) * dfqw
  END IF
END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE sf_melt
END MODULE sf_melt_mod
