#if !defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine SOILCARB_LAYERS --------------------------------------------------------
!
! Purpose : Updates carbon and nitrogen contents of the soil.
!
!Note that triffid is not compatible with soil tiling at this time, so all
!_soilt variables have their soil tile index hard-coded to 1
!
! ----------------------------------------------------------------------------
MODULE soilcarb_layers_mod

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SOILCARB_LAYERS_MOD'

CONTAINS

SUBROUTINE soilcarb_layers (land_pts, trif_pts, trif_index, forw, r_gamma,    &
                            lit_c, lit_c_t, lit_n_t, resp_frac,               &
                            resp_s_pot, cs,                                   &
                            ns_gb, neg_n, implicit_resp_correction,           &
                            burnt_soil, isunfrozen,                           &
                           !New arguments replacing USE statements
                           !prognostics
                           ns_pool_gb, n_inorg_soilt_lyrs, n_inorg_avail_pft, &
                           t_soil_soilt_acc,                                  &
                           !trif_vars_mod
                           burnt_carbon_dpm, g_burn_gb, burnt_carbon_rpm,     &
                           minl_n_gb, minl_n_pot_gb, immob_n_gb, immob_n_pot_gb, &
                           fn_gb, resp_s_diag_gb, resp_s_pot_diag_gb,         &
                           dpm_ratio_gb, n_gas_gb, resp_s_to_atmos_gb,        &
                           !p_s_parms
                           sthu_soilt)



USE jules_surface_types_mod, ONLY: npft
USE jules_soil_biogeochem_mod, ONLY: bio_hum_cn, tau_lit

USE jules_vegetation_mod, ONLY: l_nitrogen

USE jules_soil_mod, ONLY: cs_min, dzsoil, sm_levels
USE veg_param, ONLY: litc_norm
USE pftparm, ONLY: rootd_ft
USE ancil_info, ONLY: dim_cslayer, nsoilt, dim_cs1
USE jules_vegetation_mod, ONLY: l_trif_fire
USE root_frac_mod, ONLY: root_frac
USE soilcarb_mix_mod, ONLY: soilcarb_mix
USE dpm_rpm_mod, ONLY: dpm_rpm
USE decay_mod, ONLY: decay

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN).
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN)::                                                         &
  land_pts,                                                                   &
    ! Total number of land points.
  trif_pts
    ! Number of points on which TRIFFID may operate.

INTEGER, INTENT(IN) ::                                                        &
  trif_index(land_pts)
    ! Indices of land points on which TRIFFID may operate.

REAL(KIND=real_jlslsm), INTENT(IN)::                                          &
  forw,                                                                       &
    ! Forward timestep weighting.
  r_gamma
    ! Inverse timestep (/360days).

REAL(KIND=real_jlslsm), INTENT(IN)::                                          &
  lit_c(land_pts,npft),                                                       &
    ! Carbon Litter (kg C/m2/360days).
  lit_c_t(land_pts),                                                          &
    ! Total carbon litter (kg C/m2/360days).
  lit_n_t(land_pts),                                                          &
    ! Total nitrogen litter (kg N/m2/360days).
  resp_frac(land_pts,dim_cslayer),                                            &
    ! The fraction of RESP_S (soil respiration) that forms new soil C.
    ! This is the fraction that is NOT released to the atmosphere.
  isunfrozen(land_pts, dim_cslayer)
    ! Frozen status of soil

!-----------------------------------------------------------------------------
! Arguments with INTENT(INOUT).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(INOUT)::                                       &
  resp_s_pot(land_pts,dim_cslayer,5),                                         &
    ! Soil respiration (kg C/m2/360days).
  cs(land_pts,dim_cslayer,4)
    ! Soil carbon (kg C/m2).
    ! The 4 soil C pools are DPM, RPM, biomass and humus.

!-----------------------------------------------------------------------------
! Arguments with INTENT(OUT).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT)::                                         &
  ns_gb(land_pts,dim_cslayer),                                                &
    ! Total Soil N (kg N/m2).
  neg_n(land_pts),                                                            &
    ! Negative N required to prevent ns<0 (kg N).
  implicit_resp_correction(land_pts),                                         &
    ! Respiration carried to next triffid timestep to account for applying
    ! minimum soil carbon constraint (kg m-2).
  burnt_soil(land_pts)
    ! Burnt C in RPM and DPM pools, for cnsrv diagnostics (kg/m2/360days).

!New arguments replacing USE statements
!prognostics
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: ns_pool_gb(land_pts,dim_cslayer,dim_cs1)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: n_inorg_soilt_lyrs(land_pts,nsoilt,dim_cslayer)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: n_inorg_avail_pft(land_pts,npft,dim_cslayer)
REAL(KIND=real_jlslsm), INTENT(IN):: t_soil_soilt_acc(land_pts,nsoilt,sm_levels)

!trif_vars_mod
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: burnt_carbon_dpm(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: g_burn_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: burnt_carbon_rpm(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: minl_n_gb(land_pts,dim_cslayer,dim_cs1+1)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: minl_n_pot_gb(land_pts,dim_cslayer,dim_cs1+1)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: immob_n_gb(land_pts,dim_cslayer,dim_cs1+1)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: immob_n_pot_gb(land_pts,dim_cslayer,dim_cs1+1)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: fn_gb(land_pts,dim_cslayer)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: resp_s_diag_gb(land_pts,dim_cslayer,dim_cs1+1)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: resp_s_pot_diag_gb(land_pts,dim_cslayer,dim_cs1+1)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: dpm_ratio_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: n_gas_gb(land_pts,dim_cslayer)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: resp_s_to_atmos_gb(land_pts,dim_cslayer)

!p_s_parms
REAL(KIND=real_jlslsm), INTENT(IN) :: sthu_soilt(land_pts,nsoilt,sm_levels)

!-----------------------------------------------------------------------------
! Local scalar parameters.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), PARAMETER :: lit_cn    = 300.0
    ! Maximum-allowed C:N for soil litter pools.
REAL(KIND=real_jlslsm), PARAMETER :: nminl_gas = 0.01
    ! Fraction of net mineralisation of N that is lost as gas.

REAL(KIND=real_jlslsm) , PARAMETER ::                                         &
  ccdpm_min = 0.8,                                                            &
  ccdpm_max = 1.0,                                                            &
    ! Decomposable Plant Material burns between 80 to 100 %.
  ccrpm_min = 0.0,                                                            &
  ccrpm_max = 0.2
    ! Resistant Plant Material burns between 0 to 20 %.
    ! These values are also set in inferno_mod to calculate emitted_carbon_DPM
    ! and emitted_carbon_RPM, and are also set in soilcarb.

!-----------------------------------------------------------------------------
! Local variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  l,t,n,i  ! Loop counters

REAL(KIND=real_jlslsm) ::                                                     &
  resp_frac_mult46(land_pts,dim_cslayer),                                     &
    ! Respired fraction of RESP_S multiplied by 0.46.
  resp_frac_mult54(land_pts,dim_cslayer),                                     &
    ! Respired fraction of RESP_S multiplied by 0.54.
  lit_frac(dim_cslayer),                                                      &
    ! Litter fraction into each layer.
  dcs(land_pts,dim_cslayer,4),                                                &
    ! Increment to the soil carbon (kg C/m2).
  dpc_dcs(land_pts,dim_cslayer,4),                                            &
    ! Rate of change of PC with soil carbon (/360days).
  pc(land_pts,dim_cslayer,4),                                                 &
    ! Net carbon accumulation in the soil (kg C/m2/360days).
  pn(land_pts,dim_cslayer,4),                                                 &
    ! Net nitrogen accumulation in the soil (kg N/m2/360days).
  cn(land_pts,dim_cslayer,5),                                                 &
    ! C:N ratios of pools.
  mix_s(land_pts,dim_cslayer-1,4),                                            &
    ! Diffusion rate (m^2/360days).
  resp_s(land_pts,dim_cslayer,5),                                             &
    ! Soil respiration (kg C/m2/360days).
  f_root_pft(npft,dim_cslayer),                                               &
    ! Root fraction in each soil layer.
  f_root_pft_dz(npft,dim_cslayer),                                            &
    ! Exponential root fraction in each soil layer.
  cs_min_lit_cn,                                                              &
    ! cs_min/lit_cn: speeding up calculations.
  cs_min_bio_hum_cn,                                                          &
    ! cs_min/bio_hum_cn: speeding up calculations
  mix_term(land_pts,dim_cslayer,4)
    ! Mixing term for calculating mixing (kg/m^2/360days)

! Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SOILCARB_LAYERS'

!-----------------------------------------------------------------------------
!end of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Note: The 4 soil carbon pools are  1 decomposable plant material,
! 2 resistant plant material, 3 biomass, 4 humus.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Initialisation.
!-----------------------------------------------------------------------------
burnt_soil(:)     = 0.0
fn_gb(:,:)        = 1.0
n_gas_gb(:,:)     = 0.0

!Calcuate constants
cs_min_lit_cn     = cs_min / lit_cn
cs_min_bio_hum_cn = cs_min / bio_hum_cn

!Calculate root profiles for updating plant available N
DO i = 1,npft
  CALL root_frac(i, dim_cslayer, dzsoil, rootd_ft(i), f_root_pft(i,:))
  f_root_pft_dz(i,:) = f_root_pft(i,:) / dzsoil(:) /                          &
                       (f_root_pft(i,1) / dzsoil(1)) 
END DO

DO t = 1,trif_pts
  l = trif_index(t)

  DO n = 1,dim_cslayer

    !See eq. 72 & 73 in Clark et al. (2011); doi:10.5194/gmd-4-701-2011
    resp_frac_mult46(l,n) = 0.46 * resp_frac(l,n)
    resp_frac_mult54(l,n) = 0.54 * resp_frac(l,n)

    resp_s_pot(l,n,5) = resp_s_pot(l,n,1) + resp_s_pot(l,n,2)                 &
                      + resp_s_pot(l,n,3) + resp_s_pot(l,n,4)

    cn(l,n,1) = cs(l,n,1) / ns_pool_gb(l,n,1) ! C:N ratio of DPM - Prognostic
    cn(l,n,2) = cs(l,n,2) / ns_pool_gb(l,n,2) ! C:N ratio of DPM - Prognostic
    cn(l,n,3) = bio_hum_cn                    ! C:N ratio of BIO - Fixed
    cn(l,n,4) = bio_hum_cn                    ! C:N ratio of HUM - Fixed

    ! N pool BIO+HUM are diagnostic from equivalent cSoil pool
    ns_pool_gb(l,n,3) = cs(l,n,3) / bio_hum_cn
    ns_pool_gb(l,n,4) = cs(l,n,4) / bio_hum_cn

    ! These 4 lines are needed if cs is not initialised correctly
    cs(l,n,1) = MIN(MAX(cs(l,n,1),1.0e-6),1000.0)
    cs(l,n,2) = MIN(MAX(cs(l,n,2),1.0e-6),1000.0)
    cs(l,n,3) = MIN(MAX(cs(l,n,3),1.0e-6),1000.0)
    cs(l,n,4) = MIN(MAX(cs(l,n,4),1.0e-6),1000.0)

    cn(l,n,1) = MIN(MAX(cn(l,n,1),1.0e-6),1000.0)
    cn(l,n,2) = MIN(MAX(cn(l,n,2),1.0e-6),1000.0)
    cn(l,n,3) = MIN(MAX(cn(l,n,3),1.0e-6),1000.0)
    cn(l,n,4) = MIN(MAX(cn(l,n,4),1.0e-6),1000.0)

    cn(l,n,5) = SUM(cs(l,n,:)) / SUM(ns_pool_gb(l,n,:))

    minl_n_pot_gb(l,n,1) = resp_s_pot(l,n,1) / cn(l,n,1)
    minl_n_pot_gb(l,n,2) = resp_s_pot(l,n,2) / cn(l,n,2)
    minl_n_pot_gb(l,n,3) = resp_s_pot(l,n,3) / cn(l,n,3)
    minl_n_pot_gb(l,n,4) = resp_s_pot(l,n,4) / cn(l,n,4)

    minl_n_pot_gb(l,n,5) = minl_n_pot_gb(l,n,1) + minl_n_pot_gb(l,n,2) +      &
                           minl_n_pot_gb(l,n,3) + minl_n_pot_gb(l,n,4)

    immob_n_pot_gb(l,n,1) = (resp_frac_mult46(l,n)                            &
                             * resp_s_pot(l,n,1) / cn(l,n,3))                 &
                            + (resp_frac_mult54(l,n) * resp_s_pot(l,n,1)      &
                               / cn(l,n,4))
    immob_n_pot_gb(l,n,2) = (resp_frac_mult46(l,n)                            &
                             * resp_s_pot(l,n,2) / cn(l,n,3))                 &
                            + (resp_frac_mult54(l,n) * resp_s_pot(l,n,2)      &
                               / cn(l,n,4))
    immob_n_pot_gb(l,n,3) = (resp_frac_mult46(l,n)                            &
                             * resp_s_pot(l,n,3) / cn(l,n,3))                 &
                            + (resp_frac_mult54(l,n) * resp_s_pot(l,n,3)      &
                               / cn(l,n,4))
    immob_n_pot_gb(l,n,4) = (resp_frac_mult46(l,n)                            &
                             * resp_s_pot(l,n,4) / cn(l,n,3))                 &
                            + (resp_frac_mult54(l,n) * resp_s_pot(l,n,4)      &
                               / cn(l,n,4))

    immob_n_pot_gb(l,n,5) = immob_n_pot_gb(l,n,1) + immob_n_pot_gb(l,n,2)     &
                          + immob_n_pot_gb(l,n,3) + immob_n_pot_gb(l,n,4)

    !-------------------------------------------------------------------------
    ! If soil N demand greater than N_available then limit decay through fn_gb
    !-------------------------------------------------------------------------
    IF ( l_nitrogen .AND.                                                     &
         (immob_n_pot_gb(l,n,5) - minl_n_pot_gb(l,n,5)) / r_gamma             &
         > n_inorg_soilt_lyrs(l,1,n) * isunfrozen(l,n) ) THEN

      fn_gb(l,n) = ( ((minl_n_pot_gb(l,n,3) + minl_n_pot_gb(l,n,4)            &
                   - immob_n_pot_gb(l,n,3) - immob_n_pot_gb(l,n,4)) / r_gamma)&
                   + n_inorg_soilt_lyrs(l,1,n) * isunfrozen(l,n) ) /          &
                   ( (immob_n_pot_gb(l,n,1) - minl_n_pot_gb(l,n,1)) / r_gamma &
                   + (immob_n_pot_gb(l,n,2) - minl_n_pot_gb(l,n,2))           &
                   / r_gamma )
      fn_gb(l,n) = MIN(MAX(fn_gb(l,n), 0.0), 1.0)
    END IF

    minl_n_gb(l,n,1)  = fn_gb(l,n) * minl_n_pot_gb(l,n,1)
    minl_n_gb(l,n,2)  = fn_gb(l,n) * minl_n_pot_gb(l,n,2)
    minl_n_gb(l,n,3)  = minl_n_pot_gb(l,n,3)
    minl_n_gb(l,n,4)  = minl_n_pot_gb(l,n,4)

    immob_n_gb(l,n,1) = fn_gb(l,n) * immob_n_pot_gb(l,n,1)
    immob_n_gb(l,n,2) = fn_gb(l,n) * immob_n_pot_gb(l,n,2)
    immob_n_gb(l,n,3) = immob_n_pot_gb(l,n,3)
    immob_n_gb(l,n,4) = immob_n_pot_gb(l,n,4)

    resp_s(l,n,1)     = fn_gb(l,n) * resp_s_pot(l,n,1)
    resp_s(l,n,2)     = fn_gb(l,n) * resp_s_pot(l,n,2)
    resp_s(l,n,3)     = resp_s_pot(l,n,3)
    resp_s(l,n,4)     = resp_s_pot(l,n,4)

    resp_s(l,n,5)     = resp_s(l,n,1) + resp_s(l,n,2) +                       &
                        resp_s(l,n,3) + resp_s(l,n,4)
    minl_n_gb(l,n,5)  = minl_n_gb(l,n,1) + minl_n_gb(l,n,2) +                 &
                        minl_n_gb(l,n,3) + minl_n_gb(l,n,4)
    immob_n_gb(l,n,5) = immob_n_gb(l,n,1) + immob_n_gb(l,n,2) +               &
                        immob_n_gb(l,n,3) + immob_n_gb(l,n,4)

  END DO !end loop for dim_cslayer
END DO  !  trif_pts

! calculate DPM:RPM ratio of input litter Carbon
CALL dpm_rpm(land_pts, trif_pts, trif_index, lit_c, dpm_ratio_gb)

!-----------------------------------------------------------------------------
! Calculate vertical profile of litter inputs.
!-----------------------------------------------------------------------------
lit_frac(1) = dzsoil(1) *  EXP( -tau_lit * 0.5 * dzsoil(1) ) / litc_norm
DO n = 2,dim_cslayer
  lit_frac(n) = dzsoil(n) * EXP( -tau_lit *                                   &
                ( SUM(dzsoil(1:n-1)) + 0.5 * dzsoil(n) )  ) / litc_norm
END DO

!---------------------------------------------------------------------------
! Calculate mixing terms to soil nitrogen consistent with the soil carbon
! that will mix between layers on this timestep.
! This follows the soil carbon (cs) mixing (see below)
!---------------------------------------------------------------------------
CALL soilcarb_mix(land_pts, trif_pts, trif_index,                             &
                    cs(:,:,1:4) / cn(:,:,1:4),                                &
                    t_soil_soilt_acc, mix_term, mix_s)
DO t = 1,trif_pts
  l = trif_index(t)
  !---------------------------------------------------------------------------
  ! Diagnose the net local nitrogen flux into the soil
  !---------------------------------------------------------------------------
  DO n = 1,dim_cslayer
    pn(l,n,1) = dpm_ratio_gb(l) * lit_n_t(l) * lit_frac(n)                    &
               - minl_n_gb(l,n,1)
    pn(l,n,2) = (1.0 - dpm_ratio_gb(l)) * lit_n_t(l) * lit_frac(n)            &
                - minl_n_gb(l,n,2)
    pn(l,n,3) = 0.46 * immob_n_gb(l,n,5) - minl_n_gb(l,n,3)
    pn(l,n,4) = 0.54 * immob_n_gb(l,n,5) - minl_n_gb(l,n,4)
  END DO

  !---------------------------------------------------------------------------
  ! Apply mixing term to pn
  !---------------------------------------------------------------------------
  DO i = 1,4  !soil carbon pools
    DO n = 1,dim_cslayer
      pn(l,n,i) = pn(l,n,i) + mix_term(l,n,i)
    END DO
  END DO ! soil carbon pools

  !---------------------------------------------------------------------------
  ! Update soil nitrogen pools
  !---------------------------------------------------------------------------
  DO n = 1,dim_cslayer
    ns_pool_gb(l,n,1) = ns_pool_gb(l,n,1) + pn(l,n,1) / r_gamma
    ns_pool_gb(l,n,2) = ns_pool_gb(l,n,2) + pn(l,n,2) / r_gamma
    ns_pool_gb(l,n,3) = ns_pool_gb(l,n,3) + pn(l,n,3) / r_gamma
    ns_pool_gb(l,n,4) = ns_pool_gb(l,n,4) + pn(l,n,4) / r_gamma

    neg_n(l) = 0.0
    IF (ns_pool_gb(l,n,1) <  cs_min_lit_cn) THEN
      neg_n(l)  = neg_n(l) + ns_pool_gb(l,n,1) - cs_min_lit_cn
    END IF
    IF (ns_pool_gb(l,n,2) <  cs_min_lit_cn) THEN
      neg_n(l)  = neg_n(l) + ns_pool_gb(l,n,2) - cs_min_lit_cn
    END IF
    IF (ns_pool_gb(l,n,3) <  cs_min_bio_hum_cn) THEN
      neg_n(l) = neg_n(l) + ns_pool_gb(l,n,3) - cs_min_bio_hum_cn
    END IF
    IF (ns_pool_gb(l,n,4) <  cs_min_bio_hum_cn) THEN
      neg_n(l) = neg_n(l) + ns_pool_gb(l,n,4) - cs_min_bio_hum_cn
    END IF

    ns_pool_gb(l,n,1) = MAX(ns_pool_gb(l,n,1),cs_min_lit_cn)
    ns_pool_gb(l,n,2) = MAX(ns_pool_gb(l,n,2),cs_min_lit_cn)
    ns_pool_gb(l,n,3) = MAX(ns_pool_gb(l,n,3),cs_min_bio_hum_cn)
    ns_pool_gb(l,n,4) = MAX(ns_pool_gb(l,n,4),cs_min_bio_hum_cn)

    ! calculate mineralised gas emissions
    ! let gas be negative if required
    n_gas_gb(l,n) = nminl_gas * MAX( minl_n_gb(l,n,5) - immob_n_gb(l,n,5),    &
                                     0.0)

    !-------------------------------------------------------------------------
    ! Update inorganic N
    !-------------------------------------------------------------------------
    n_inorg_soilt_lyrs(l,1,n) = n_inorg_soilt_lyrs(l,1,n)                     &
                                + (minl_n_gb(l,n,5) - immob_n_gb(l,n,5)       &
                                   - n_gas_gb(l,n)) / r_gamma
    ns_gb(l,n) = ns_pool_gb(l,n,1) + ns_pool_gb(l,n,2) +                      &
                 ns_pool_gb(l,n,3) + ns_pool_gb(l,n,4)
    DO i = 1,npft !Update plant available inorganic nitrogen.
      n_inorg_avail_pft(l,i,n) = n_inorg_avail_pft(l,i,n)                     &
                                 + MAX( (f_root_pft_dz(i,n)                   &
                                         * (minl_n_gb(l,n,5)                  &
                                             - immob_n_gb(l,n,5)              &
                                             - n_gas_gb(l,n)) / r_gamma),     &
                                         n_inorg_avail_pft(l,i,n) * (-1.0) )
    END DO !npft
    
    !add on neg_n after modifying n_inorg so you dont break the N balance
    n_gas_gb(l,n) = n_gas_gb(l,n) + neg_n(l) * r_gamma
  
  END DO !dim_cslayer
END DO

!---------------------------------------------------------------------------
! Calculate a mixing term
! (These are diffusion terms: first term on right hand side of
! equations 10-13 in Burke et al 2017. Link to paper:
! https://www.geosci-model-dev.net/10/959/2017/gmd-10-959-2017.pdf)
!---------------------------------------------------------------------------
CALL soilcarb_mix(land_pts, trif_pts, trif_index, cs(:,:,1:4),                &
                  t_soil_soilt_acc, mix_term, mix_s)

DO t = 1,trif_pts
  l = trif_index(t)
  DO n = 1,dim_cslayer

    !-------------------------------------------------------------------------
    ! Diagnose the gridbox mean soil-to-atmosphere respiration carbon flux
    ! [kg m-2 (360 days)-1]
    !-------------------------------------------------------------------------
    resp_s_to_atmos_gb(l,n) = (1.0 - resp_frac(l,n)) * resp_s(l,n,5)

    !-------------------------------------------------------------------------
    ! Diagnose the net local carbon flux into the soil
    !-------------------------------------------------------------------------
    pc(l,n,1) = (dpm_ratio_gb(l) * lit_c_t(l) * lit_frac(n)) - resp_s(l,n,1)
    pc(l,n,2) = ((1.0 - dpm_ratio_gb(l)) * lit_c_t(l) * lit_frac(n))          &
                - resp_s(l,n,2)
    pc(l,n,3) = (resp_frac_mult46(l,n) * resp_s(l,n,5)) - resp_s(l,n,3)
    pc(l,n,4) = (resp_frac_mult54(l,n) * resp_s(l,n,5)) - resp_s(l,n,4)

    !-------------------------------------------------------------------------
    ! Variables required for the implicit and equilibrium calculations
    !-------------------------------------------------------------------------
    dpc_dcs(l,n,1) = resp_s(l,n,1) / cs(l,n,1)
    dpc_dcs(l,n,2) = resp_s(l,n,2) / cs(l,n,2)
    dpc_dcs(l,n,3) = resp_s(l,n,3) / cs(l,n,3)
    dpc_dcs(l,n,4) = resp_s(l,n,4) / cs(l,n,4)

    !-------------------------------------------------------------------------
    ! Save current value of soil carbon
    !-------------------------------------------------------------------------
    dcs(l,n,1) = cs(l,n,1)
    dcs(l,n,2) = cs(l,n,2)
    dcs(l,n,3) = cs(l,n,3)
    dcs(l,n,4) = cs(l,n,4)

  END DO !end soil layers

  !---------------------------------------------------------------------------
  ! Apply mixing term to soil carbon
  !---------------------------------------------------------------------------
  DO i = 1,4  !soil carbon pools
    DO n = 1,dim_cslayer
      pc(l,n,i) = pc(l,n,i) + mix_term(l,n,i)
    END DO
  END DO ! soil carbon pools
END DO ! trif_pts

!-----------------------------------------------------------------------------
! Update soil carbon
!-----------------------------------------------------------------------------
DO n = 1,dim_cslayer
  CALL decay(land_pts, trif_pts, trif_index, forw, r_gamma,                   &
             dpc_dcs(:,n,:), pc(:,n,:), cs(:,n,:))
END DO

!-----------------------------------------------------------------------------
! Remove burnt litter from top soil layer only
!-----------------------------------------------------------------------------
IF (l_trif_fire) THEN
  DO t = 1,trif_pts
    l = trif_index(t)
    burnt_carbon_dpm(l) = MAX(g_burn_gb(l) *                                  &
                          (cs(l,1,1) * (ccdpm_min + (ccdpm_max - ccdpm_min)   &
                          * (1.0 - (sthu_soilt(l,1,1))))) ,0.0)
    burnt_carbon_rpm(l) = MAX(g_burn_gb(l) *                                  &
                          (cs(l,1,2) * (ccrpm_min + (ccrpm_max - ccrpm_min)   &
                          * (1.0 - (sthu_soilt(l,1,1))))) ,0.0)
    cs(l,1,1)     = cs(l,1,1) - burnt_carbon_dpm(l) / r_gamma
    cs(l,1,2)     = cs(l,1,2) - burnt_carbon_rpm(l) / r_gamma
    burnt_soil(l) = burnt_carbon_dpm(l) + burnt_carbon_rpm(l)
  END DO
END IF
!-----------------------------------------------------------------------------
! Apply implicit correction to the soil respiration rate.
!-----------------------------------------------------------------------------
DO t = 1,trif_pts
  l = trif_index(t)
  DO n = 1,dim_cslayer

    dcs(l,n,1) = cs(l,n,1) - dcs(l,n,1)
    dcs(l,n,2) = cs(l,n,2) - dcs(l,n,2)
    dcs(l,n,3) = cs(l,n,3) - dcs(l,n,3)
    dcs(l,n,4) = cs(l,n,4) - dcs(l,n,4)

    resp_s(l,n,1) = resp_s(l,n,1) + forw * dpc_dcs(l,n,1) * dcs(l,n,1)
    resp_s(l,n,2) = resp_s(l,n,2) + forw * dpc_dcs(l,n,2) * dcs(l,n,2)
    resp_s(l,n,3) = resp_s(l,n,3) + forw * dpc_dcs(l,n,3) * dcs(l,n,3)
    resp_s(l,n,4) = resp_s(l,n,4) + forw * dpc_dcs(l,n,4) * dcs(l,n,4)

  END DO

  implicit_resp_correction(l) = ( r_gamma * SUM(dcs(l,1:dim_cslayer,1:4))     &
                                  - SUM(pc(l,1:dim_cslayer,1:4)) ) / r_gamma
END DO

! Sum total respiration.
DO t = 1,trif_pts
  l = trif_index(t)
  DO n = 1,dim_cslayer
    resp_s(l,n,5)             = resp_s(l,n,1) + resp_s(l,n,2) +               &
                                resp_s(l,n,3) + resp_s(l,n,4)

    resp_s_diag_gb(l,n,:)     = resp_s(l,n,:)
    resp_s_pot_diag_gb(l,n,:) = resp_s_pot(l,n,:)

    resp_s_pot(l,n,1) = resp_s(l,n,1)
    resp_s_pot(l,n,2) = resp_s(l,n,2)
    resp_s_pot(l,n,3) = resp_s(l,n,3)
    resp_s_pot(l,n,4) = resp_s(l,n,4)

    resp_s_pot(l,n,5) = resp_s_pot(l,n,1) + resp_s_pot(l,n,2)  +              &
                        resp_s_pot(l,n,3) + resp_s_pot(l,n,4)

  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE soilcarb_layers
END MODULE soilcarb_layers_mod
#endif
