! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!     SUBROUTINE SF_FLUX -----------------------------------------------
! Description:
! Subroutines SF_FLUX to calculate explicit surface fluxes of
! heat and moisture
!-----------------------------------------------------------------------
MODULE sf_flux_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SF_FLUX_MOD'

CONTAINS
SUBROUTINE sf_flux (                                                          &
 points,surft_pts,fld_sea,pts_index,surft_index,                              &
 nsnow,n,canhc,dzsurf,hcons,qs1_elev,qstar,q_elev,radnet,resft,               &
 rhokh_1,l_soil_point,snowdepth,tile_frac,timestep,                           &
 t_elev,ts1_elev,tstar,vfrac,rhokh_can,                                       &
 z0h,z0m_eff,zdt,z1_tq,lh0,emis_surft,emis_soil,                              &
 salinityfactor,anthrop_heat,scaling_urban,l_vegdrag,                         &
 fqw_1_gb,ftl_1_gb,                                                           &
 alpha1,ashtf_prime,fqw_1,epot,ftl_1,rhokm_1,vshr,tau_1,dtstar,               &
 sea_point, sf_diag                                                           &
 )

USE atm_fields_bounds_mod, ONLY: tdims
USE theta_field_sizes, ONLY: t_i_length

USE csigma, ONLY: sbcon
USE water_constants_mod, ONLY: lc, tm
USE planet_constants_mod, ONLY: cp, r, c_virtual,epsil=>repsilon
USE jules_sea_seaice_mod, ONLY: l_use_dtstar_sea
USE jules_snow_mod, ONLY: snow_hcon
USE jules_surface_mod, ONLY: grcp,ls
USE jules_surface_types_mod, ONLY: urban_roof
USE jules_vegetation_mod, ONLY: l_vegcan_soilfx
USE switches_urban, ONLY: l_moruses_storage
USE jules_surface_mod, ONLY: l_aggregate, l_epot_corr
USE jules_science_fixes_mod, ONLY: l_fix_moruses_roof_rad_coupling
USE sf_diags_mod, ONLY: strnewsfdiag

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                        &
 points                                                                       &
                           ! IN Total number of points.
,surft_pts                                                                    &
                           ! IN Number of tile points.
,pts_index(points)                                                            &
                           ! IN Index of points.
,surft_index(points)                                                          &
                           ! IN Index of tile points.
,nsnow(points)                                                                &
                           ! IN Number of snow layers
,n                         ! IN Tile number.
                           ! For sea and sea-ice this = 0

LOGICAL, INTENT(IN) ::                                                        &
 l_soil_point(points)                                                         &
                      ! IN Boolean to test for soil points
,l_vegdrag
                      ! IN Option for vegetation canopy drag scheme.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
 fld_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                           ! IN Fraction of land or sea
,canhc(points)                                                                &
                           ! IN Areal heat capacity of canopy (J/K/m2).
,dzsurf(points)                                                               &
                           ! IN Surface layer thickness (m).
,qs1_elev(points)                                                             &
                           ! IN Sat. specific humidity qsat(T_ELEV,PSTAR
,qstar(points)                                                                &
                           ! IN Surface qsat.
,q_elev(points)                                                               &
                           ! IN Total water content of lowest
                           !    atmospheric layer (kg per kg air).
,radnet(points)                                                               &
                           ! IN Net surface radiation (W/m2) positive
                           !    downwards
,resft(points)                                                                &
                           ! IN Total resistance factor.
,rhokh_1(points)                                                              &
                           ! IN Surface exchange coefficient.
,rhokm_1(points)                                                              &
!                          ! IN Surface momentum exchange coefficient
,snowdepth(points)                                                            &
                           ! IN Snow depth (on ground) (m)
,tile_frac(points)                                                            &
                           ! IN Tile fraction.
,timestep                                                                     &
                           ! IN Timestep (s).
,vshr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                           ! IN Magnitude of surface-to-lowest-level
!                          !    wind shear.
,t_elev(points)                                                               &
                           ! IN Liquid/frozen water temperature for
                           !     lowest atmospheric layer (K).
,ts1_elev(points)                                                             &
                           ! IN Temperature of surface layer (K).
,tstar(points)                                                                &
                           ! IN Surface temperature (K).
,vfrac(points)                                                                &
                           ! IN Fractional canopy coverage.
,rhokh_can(points)                                                            &
                           ! IN Exchange coefficient for canopy air
                           !     to surface
,z0h(points)                                                                  &
                           ! IN Roughness length for heat and moisture
,z0m_eff(points)                                                              &
                           ! IN Effective roughness length for momentum
,zdt(points)                                                                  &
                           ! IN Difference between the canopy height and
                           !    displacement height (m)
,z1_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                           ! IN Height of lowest atmospheric level (m).
,emis_surft(points)                                                           &
                           ! IN Emissivity for land tiles
,emis_soil(points)                                                            &
                           ! IN Emissivity of underlying soil
,lh0                                                                          &
                           ! IN Latent heat for snow free surface
                           !    =LS for sea-ice, =LC otherwise
,salinityfactor                                                               &
                           ! IN Factor allowing for the effect of the
                           !    salinity of sea water on the
                           !    evaporative flux.
,anthrop_heat(points)                                                         &
                           ! IN Anthropogenic contribution to surface
                           !    heat flux (W/m2). Zero except for
                           !    urban and L_ANTHROP_HEAT=.true.
                           !    or for urban_canyon & urban_roof when
                           !    l_urban2t=.true.
,scaling_urban(points)
                           ! IN MORUSES: ground heat flux scaling;
                           ! canyon tile only coupled to soil.
                           ! This equals 1.0 except for urban tiles when
                           ! MORUSES is used.

REAL(KIND=real_jlslsm) ::                                                     &
 hcons(points)
                           ! IN Soil thermal conductivity (W/m/K).

REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
 fqw_1_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                &
                          ! INOUT GBM surface flux of
                          !       QW (kg/m2/s).
,ftl_1_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                          ! INOUT GBM surface flux of TL.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
 ashtf_prime(points)                                                          &
                           ! OUT Adjusted SEB coefficient
,alpha1(points)                                                               &
                           ! OUT Gradient of saturated specific humidity
                           !     with respect to temperature between the
                           !     bottom model layer and the surface.
,fqw_1(points)                                                                &
                           ! OUT Local surface flux of QW (kg/m2/s).
,epot(points)                                                                 &
                           ! OUT
,ftl_1(points)                                                                &
                           ! OUT Local surface flux of TL.
,tau_1(points)                                                                &
!                          ! OUT Local surface momentum flux
,dtstar(points)            ! OUT Change in TSTAR over timestep

REAL(KIND=real_jlslsm) ::                                                     &
 sea_point                 ! =1.0 IF SEA POINT, =0.0 OTHERWISE

! Diagnostics
TYPE (strnewsfdiag), INTENT(INOUT) :: sf_diag


! Workspace
REAL(KIND=real_jlslsm) ::                                                     &
 ashtf(points)                                                                &
                           ! Coefficient to calculate surface
                           ! heat flux into soil (W/m2/K).
,dtstar_pot(points)                                                           &
                           ! Change in TSTAR over timestep that is
                           ! appropriate for the potential evaporation
,surf_ht_flux              ! Flux of heat from surface to sub-surface

! Scalars
INTEGER ::                                                                    &
 i,j                                                                          &
                           ! Horizontal field index.
,k                                                                            &
                           ! Tile field index.
,l                         ! Points field index.

REAL(KIND=real_jlslsm) ::                                                     &
 ds_ratio                                                                     &
                           ! 2 * snowdepth / depth of top soil layer.
,d_t                                                                          &
                           ! Temporary in calculation of alpha1.
,lh                        ! Latent heat (J/K/kg).

REAL(KIND=real_jlslsm) :: lambda
                           ! Attenuation factor for influence of soil
                           ! temperature

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SF_FLUX'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! This allows fluxes to be corrected for sea points
IF (l_use_dtstar_sea) THEN
  sea_point = 0.0
END IF

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(SHARED)                                                         &
!$OMP PRIVATE(l,k,d_t,ds_ratio,lambda,lh,surf_ht_flux,i,j)

!-----------------------------------------------------------------------
!!  0 initialise
!-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
DO l = 1,points
  ashtf(l)  = 0.0
  alpha1(l) = 0.0
  ftl_1(l)  = 0.0
  epot(l)   = 0.0
END DO
!$OMP END DO

!-----------------------------------------------------------------------
!!  1 Calculate gradient of saturated specific humidity for use in
!!    calculation of surface fluxes
!-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
DO k = 1,surft_pts
  l = surft_index(k)
  d_t = tstar(l) - t_elev(l)
  IF (d_t > 0.05 .OR. d_t < -0.05) THEN
    alpha1(l) = (qstar(l) - qs1_elev(l)) / d_t
  ELSE IF (t_elev(l) > tm) THEN
    alpha1(l) = epsil * lc * qs1_elev(l) *                                    &
                (1.0 + c_virtual * qs1_elev(l)) /                             &
                ( r * t_elev(l) * t_elev(l))
  ELSE
    alpha1(l) = epsil * ls * qs1_elev(l) *                                    &
                (1.0 + c_virtual * qs1_elev(l)) /                             &
                ( r * t_elev(l) * t_elev(l))
  END IF
END DO
!$OMP END DO

! MORUSES: Uncouple the roof for perfect insulation. hcons should only be zero
! to change the "conductive" coupling to "uncoupled" otherwise it is radiatively
! coupled.
IF ( .NOT. l_aggregate .AND. l_moruses_storage .AND. n == urban_roof ) THEN
  IF ( l_fix_moruses_roof_rad_coupling ) THEN
!$OMP DO SCHEDULE(STATIC)
    DO l = 1,points
      IF (vfrac(l) == 0.0) THEN
        hcons(l) = 0.0
      END IF
    END DO
!$OMP END DO
  ELSE
!$OMP DO SCHEDULE(STATIC)
    DO l = 1, points
      hcons(l) = 0.0
    END DO
!$OMP END DO
  END IF
END IF

!$OMP DO SCHEDULE(STATIC)
DO k = 1,surft_pts
  l = surft_index(k)
  ! Except when n == urban_canyon when MORUSES is used
  ! scaling_urban(l) = 1.0
  ashtf(l) = 2.0 * hcons(l) / dzsurf(l)
  ashtf(l) = ashtf(l) * scaling_urban(l)

  IF (snowdepth(l) > 0.0 .AND. l_soil_point(l) .AND. nsnow(l) == 0 ) THEN
    IF ( l_moruses_storage .AND. n == urban_roof ) THEN
      ! This required as HCONS(L) = 0 in this case.
      ashtf(l) =  0.0
    ELSE
      ds_ratio = 2.0 * snowdepth(l) / dzsurf(l)
      IF (ds_ratio <= 1.0) THEN
        ashtf(l) =  ashtf(l) /                                                &
                    (1.0 + ds_ratio * (hcons(l) / snow_hcon - 1.0))
      ELSE
        ashtf(l) =  ashtf(l) * snow_hcon / hcons(l)
      END IF
    END IF
  END IF
END DO
!$OMP END DO

! If conduction in the soil beneath the vegetative canopy is not
! explicitly considered, the attentuation facor may be set equal
! to 1 everywhere.
lambda = 1.0

!$OMP DO SCHEDULE(STATIC)
DO k = 1,surft_pts
  l = surft_index(k)
  j=(pts_index(l) - 1) / t_i_length + 1
  i = pts_index(l) - (j-1) * t_i_length

  ! Calculate the attenuation factor if different from 1.
  IF (l_vegcan_soilfx)                                                        &
    lambda = 2.0 * hcons(l) / dzsurf(l) /                                     &
             ( 2.0 * hcons(l) / dzsurf(l) + rhokh_can(l) +                    &
             4.0 * emis_soil(l) * emis_surft(l) * sbcon *                     &
             tstar(l)**3 )

  lh = lh0
  IF (snowdepth(l) > 0.0) lh = ls

  IF (l_vegdrag) THEN
    ftl_1(l) = rhokh_1(l) * (tstar(l) - t_elev(l) -                           &
                    grcp * (z1_tq(i,j) + zdt(l) - z0h(l)))
  ELSE
    ftl_1(l) = rhokh_1(l) * (tstar(l) - t_elev(l) -                           &
                    grcp * (z1_tq(i,j) + z0m_eff(l) - z0h(l)))
  END IF
  epot(l) = rhokh_1(l) * (salinityfactor * qstar(l) - q_elev(l))
  fqw_1(l) = resft(l) * epot(l)

  surf_ht_flux = ((1.0 - vfrac(l)) * ashtf(l) +                               &
                    vfrac(l) * rhokh_can(l) * lambda ) *                      &
                                   (tstar(l) - ts1_elev(l)) +                 &
                 vfrac(l) * emis_soil(l) * emis_surft(l) * sbcon *            &
                  lambda * (tstar(l)**4.0 - ts1_elev(l)**4.0)


  ashtf_prime(l) = 4.0 * (1.0 + lambda * emis_soil(l) * vfrac(l)) *           &
                          emis_surft(l) * sbcon * tstar(l)**3.0 +             &
                          lambda * vfrac(l) * rhokh_can(l) +                  &
                          (1.0 - vfrac(l)) * ashtf(l) + canhc(l) / timestep

  dtstar(l) = (radnet(l) + anthrop_heat(l) - cp * ftl_1(l) -                  &
                         lh * fqw_1(l) - surf_ht_flux)  /                     &
               ( rhokh_1(l) * (cp + lh * alpha1(l) * resft(l)) +              &
                    ashtf_prime(l) )

  ! Correction to surface fluxes due to change in surface temperature
  ftl_1(l) = ftl_1(l) + rhokh_1(l) * dtstar(l) * (1.0 - sea_point)
  fqw_1(l) = fqw_1(l) + resft(l) * rhokh_1(l) * alpha1(l) *                   &
                        dtstar(l) * (1.0 - sea_point)

  IF (l_epot_corr) THEN
    dtstar_pot(l) = (radnet(l) + anthrop_heat(l) - cp * ftl_1(l) -            &
                         lh * epot(l) - surf_ht_flux)  /                      &
               ( rhokh_1(l) * (cp + lh * alpha1(l)) + ashtf_prime(l) )
    epot(l) = epot(l) + rhokh_1(l) * alpha1(l) * dtstar_pot(l)
  ELSE
    epot(l) = epot(l) + rhokh_1(l) * alpha1(l) * dtstar(l)
  END IF

  ftl_1_gb(i,j) = ftl_1_gb(i,j) + fld_sea(i,j) * tile_frac(l) * ftl_1(l)
  fqw_1_gb(i,j) = fqw_1_gb(i,j) + fld_sea(i,j) * tile_frac(l) * fqw_1(l)

END DO
!$OMP END DO

! Calculate surface momentum flux
IF (sf_diag%l_tau_surft .OR. sf_diag%l_tau_1) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 1,surft_pts
    l = surft_index(k)
    j=(pts_index(l) - 1) / t_i_length + 1
    i = pts_index(l) - (j-1) * t_i_length
    tau_1(l) = rhokm_1(l) * vshr(i,j)
  END DO
!$OMP END DO
END IF

IF (sf_diag%l_tau_1) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 1,surft_pts
    l = surft_index(k)
    j=(pts_index(l) - 1) / t_i_length + 1
    i = pts_index(l) - (j-1) * t_i_length
    sf_diag%tau_1(i,j) = sf_diag%tau_1(i,j) +                                 &
                    fld_sea(i,j) * tile_frac(l) * tau_1(l)
  END DO
!$OMP END DO
END IF

!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE sf_flux
END MODULE sf_flux_mod
