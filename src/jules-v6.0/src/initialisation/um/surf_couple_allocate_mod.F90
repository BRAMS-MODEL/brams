#if defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Initialise model from JULES namelist data
!
MODULE surf_couple_allocate_mod

IMPLICIT NONE

! Description:
!  init_jules_from_namelists initialises the relevant JULES arrays
!  from data read in via namelists.
!
! Method:
!  Arrays are allocated in the first call, then initialised from namelist data.
!  Should only be called after all the JULES namelists have been read in.
!
! Code Owner: Please refer to ModuleLeaders.txt and UM file CodeOwners.txt
!
! Code Description:
!  Language: FORTRAN 95.
!  This code is written to UMDP3 v8.5 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                       &
  ModuleName='SURF_COUPLE_ALLOCATE_MOD'

CONTAINS

SUBROUTINE surf_couple_allocate (land_field, ntiles, sm_levels,               &
                                      nice, nice_use)

USE atm_fields_bounds_mod,    ONLY: tdims, vdims, udims
USE theta_field_sizes,        ONLY: t_i_length, t_j_length,                   &
                                    u_i_length, u_j_length,                   &
                                    v_i_length, v_j_length
USE ancil_info,               ONLY: dim_cs1,                                  &
                                    land_pts, nsurft
USE jules_soil_mod,           ONLY: sm_levels_jules => sm_levels

USE allocate_jules_arrays_mod, ONLY: allocate_jules_arrays
USE jules_surface_types_mod,    ONLY: npft, nnvg, soil
USE jules_surface_mod,          ONLY: l_aggregate

USE jules_vegetation_mod, ONLY: l_triffid, l_phenol, can_rad_mod

USE pftparm_io, ONLY:                                                         &
! namelist variables:
  c3_io,           orient_io,        a_wl_io,                                 &
  a_ws_io,         albsnc_max_io,    albsnc_min_io,                           &
  albsnf_maxu_io,  albsnf_max_io,    albsnf_maxl_io,                          &
  alpha_io,        alniru_io,        alnir_io,                                &
  alnirl_io,       alparu_io,        alpar_io,                                &
  alparl_io,       b_wl_io,          catch0_io,                               &
  dcatch_dlai_io,  dgl_dm_io,        dgl_dt_io,                               &
  dqcrit_io,       dz0v_dh_io,       z0v_io,                                  &
  eta_sl_io,       fd_io,            fsmc_of_io,                              &
  fsmc_p0_io,      f0_io,                                                     &
  g_leaf_0_io,     glmin_io,         infil_f_io,                              &
  kext_io,         kpar_io,          neff_io,                                 &
  nl0_io,          nr_nl_io,         ns_nl_io,                                &
  omegau_io,       omega_io,         omegal_io,                               &
  omniru_io,       omnir_io,         omnirl_io,                               &
  r_grow_io,       rootd_ft_io,      sigl_io,                                 &
  tleaf_of_io,     tlow_io,          tupp_io,                                 &
  emis_pft_io,     z0hm_pft_io,      z0hm_classic_pft_io,                     &
  dust_veg_scj_io, fl_o3_ct_io,      dfp_dcuo_io,                             &
  ci_st_io,        gpp_st_io,        ief_io,                                  &
  tef_io,          mef_io,           aef_io,                                  &
  avg_ba_io,       ccleaf_max_io,    ccleaf_min_io,                           &
  ccwood_max_io,   ccwood_min_io,    fef_bc_io,                               &
  fef_ch4_io,      fef_co2_io,       fef_co_io,                               &
  fef_nox_io,      fef_oc_io,        fef_so2_io,                              &
  lma_io,          nmass_io,         q10_leaf_io,                             &
  vsl_io,          vint_io,          kn_io,                                   &
  lai_alb_lim_io,  hw_sw_io,         nr_io,                                   &
  nsw_io,          can_struct_a_io,  gsoil_f_io

USE pftparm, ONLY:                                                            &
  c3,              orient,           a_wl,                                    &
  a_ws,            albsnc_max,       albsnc_min,                              &
  albsnf_maxu,     albsnf_max,       albsnf_maxl,                             &
  alpha,           alniru,           alnir,                                   &
  alnirl,          alparu,           alpar,                                   &
  alparl,          b_wl,             catch0,                                  &
  dcatch_dlai,     dgl_dm,           dgl_dt,                                  &
  dqcrit,          dz0v_dh,          z0v,                                     &
  eta_sl,          fd,               fsmc_of,                                 &
  fsmc_p0,         f0,                                                        &
  g_leaf_0,        glmin,            infil_f,                                 &
  kext,            kpar,             neff,                                    &
  nl0,             nr_nl,            ns_nl,                                   &
  omegau,          omega,            omegal,                                  &
  omniru,          omnir,            omnirl,                                  &
  r_grow,          rootd_ft,         sigl,                                    &
  tleaf_of,        tlow,             tupp,                                    &
  emis_pft,                                                                   &
  dust_veg_scj,    fl_o3_ct,         dfp_dcuo,                                &
  ci_st,           gpp_st,           ief,                                     &
  tef,             mef,              aef,                                     &
  avg_ba,          ccleaf_max,       ccleaf_min,                              &
  ccwood_max,      ccwood_min,       fef_bc,                                  &
  fef_ch4,         fef_co2,          fef_co,                                  &
  fef_nox,         fef_oc,           fef_so2,                                 &
  lma,             nmass,            q10_leaf,                                &
  vsl,             vint,             kn,                                      &
  lai_alb_lim,     hw_sw,            nr,                                      &
  nsw,             can_struct_a,     gsoil_f

USE nvegparm_io, ONLY:                                                        &
! namelist variables:
  albsnc_nvg_io,   albsnf_nvgu_io,  albsnf_nvg_io,                            &
  albsnf_nvgl_io,  catch_nvg_io,    gs_nvg_io,                                &
  infil_nvg_io,    z0_nvg_io,       ch_nvg_io,                                &
  vf_nvg_io,       emis_nvg_io,     z0hm_nvg_io,                              &
  z0hm_classic_nvg_io

USE nvegparm, ONLY:                                                           &
  albsnc_nvg,      albsnf_nvgu,     albsnf_nvg,                               &
  albsnf_nvgl,     catch_nvg,       gs_nvg,                                   &
  infil_nvg,       z0_nvg,          ch_nvg,                                   &
  vf_nvg,          emis_nvg

! For use with jules_pftparm and jules_nvegparm
USE c_z0h_z0m, ONLY:                                                          &
  z0h_z0m,         z0h_z0m_classic

USE trif_io, ONLY:                                                            &
! namelist variables:
  crop_io,         g_area_io,        g_grow_io,                               &
  g_root_io,       g_wood_io,        lai_max_io,                              &
  lai_min_io,      alloc_fast_io,    alloc_med_io,                            &
  alloc_slow_io,   dpm_rpm_ratio_io, retran_l_io,                             &
  retran_r_io

USE trif, ONLY:                                                               &
  crop,            g_area,           g_grow,                                  &
  g_root,          g_wood,           lai_max,                                 &
  lai_min,         alloc_fast,       alloc_med,                               &
  alloc_slow,      dpm_rpm_ratio,    retran_l,                                &
  retran_r

USE jules_irrig_mod, ONLY:                                                    &
   irrtiles, irrigtiles, nirrtile, l_irrig_dmd

USE c_elevate, ONLY:                                                          &
! namelist variables:
  surf_hgt_io

USE dust_parameters_mod, ONLY: z0_soil

! This check must be run once we have access to sm_levels
USE jules_soil_mod, ONLY: check_jules_soil

!Pointer association functions
USE crop_vars_mod, ONLY: crop_vars_assoc
USE p_s_parms, ONLY: psparms_assoc
USE top_pdm, ONLY: top_pdm_assoc
USE fire_vars_mod, ONLY: fire_vars_assoc
USE ancil_info, ONLY: ancil_info_assoc
USE trif_vars_mod, ONLY: trif_vars_assoc
USE soil_ecosse_vars_mod, ONLY: soil_ecosse_vars_assoc
USE aero, ONLY: aero_assoc
USE urban_param_mod, ONLY:urban_param_assoc
USE prognostics, ONLY: prognostics_assoc
USE trifctl, ONLY: trifctl_assoc
USE coastal, ONLY: coastal_assoc
USE jules_vars_mod, ONLY: jules_vars_assoc

!USE in instances of the JULES types
USE atm_fields_mod, ONLY: crop_vars, crop_vars_data,                          &
                          psparms, psparms_data,                              &
                          toppdm, top_pdm_data,                               &
                          fire_vars, fire_vars_data,                          &
                          ainfo, ainfo_data,                                  &
                          trif_vars, trif_vars_data,                          &
                          soil_ecosse_vars_data, soilecosse,                  &
                          aero_data, aerotype,                                &
                          urban_param_data, urban_param,                      &
                          progs_data, progs,                                  &
                          toppdm, top_pdm_data,                               &
                          trifctltype, trifctl_data,                          &
                          coastal_data, coast,                                &
                          jules_vars_data, jules_vars

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: land_field   ! Number of land points
INTEGER, INTENT(IN) :: ntiles       ! Number of surface tiles
INTEGER, INTENT(IN) :: sm_levels    ! Number of soil layers
INTEGER, INTENT(IN) :: nice         ! Number of sea ice categories
INTEGER, INTENT(IN) :: nice_use     ! Number of sea ice categories used
                                    ! fully in surface exch and radiation

! Local variables
INTEGER             :: i            ! Looper

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER   :: RoutineName='SURF_COUPLE_ALLOCATE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!Set sizes of dimensions

! Compute lengths in i and j direction. This is the earliest place that they
! are needed. They will be kept in the module from here onward.
t_i_length = tdims%i_end - tdims%i_start + 1
t_j_length = tdims%j_end - tdims%j_start + 1

u_i_length = udims%i_end - udims%i_start + 1
u_j_length = udims%j_end - udims%j_start + 1

v_i_length = vdims%i_end - vdims%i_start + 1
v_j_length = vdims%j_end - vdims%j_start + 1

!We need to set the sizes of a few dimensions. These are also exist in
!atm_step_local and set by the same logic in atm_step_phys_init. However,
!This happens after the call to this routine.
!Pending a better approach to JULES memory management we'll repeat the logic
!here
IF (l_triffid) THEN
  dim_cs1 = 4
ELSE
  dim_cs1 = 1
END IF

!Copy dimension sizes from UM modules to JULES ones.
!Rather ugly but removes the need for ifdefs elsewhere
land_pts        = land_field
nsurft          = ntiles
sm_levels_jules = sm_levels

! Allocate JULES arrays
CALL allocate_jules_arrays(crop_vars_data,psparms_data,top_pdm_data,          &
                           fire_vars_data,ainfo_data,trif_vars_data,          &
                           soil_ecosse_vars_data, aero_data, urban_param_data,&
                           progs_data,trifctl_data, coastal_data,             &
                           jules_vars_data)

!Associate the JULES pointer types to the data types.
!Doing this here reflects when the same is done in standalone JULES
!in init.F90
CALL crop_vars_assoc(crop_vars, crop_vars_data)
CALL psparms_assoc(psparms,psparms_data)
CALL top_pdm_assoc(toppdm, top_pdm_data)
CALL fire_vars_assoc(fire_vars, fire_vars_data)
CALL ancil_info_assoc(ainfo, ainfo_data)
CALL trif_vars_assoc(trif_vars, trif_vars_data)
CALL soil_ecosse_vars_assoc(soilecosse, soil_ecosse_vars_data)
CALL aero_assoc(aerotype, aero_data)
CALL urban_param_assoc(urban_param, urban_param_data)
CALL prognostics_assoc(progs, progs_data)
CALL trifctl_assoc(trifctltype, trifctl_data)
CALL coastal_assoc(coast, coastal_data)
CALL jules_vars_assoc(jules_vars, jules_vars_data)

! Begin initialising arrays from:
! jules_pftparm
c3(:)           = c3_io(1:npft)
orient(:)       = orient_io(1:npft)
a_wl(:)         = a_wl_io(1:npft)
a_ws(:)         = a_ws_io(1:npft)
albsnc_max(:)   = albsnc_max_io(1:npft)
albsnc_min(:)   = albsnc_min_io(1:npft)
albsnf_maxu(:)  = albsnf_maxu_io(1:npft)
albsnf_max(:)   = albsnf_max_io(1:npft)
albsnf_maxl(:)  = albsnf_maxl_io(1:npft)
alpha(:)        = alpha_io(1:npft)
alniru(:)       = alniru_io(1:npft)
alnir(:)        = alnir_io(1:npft)
alnirl(:)       = alnirl_io(1:npft)
alparu(:)       = alparu_io(1:npft)
alpar(:)        = alpar_io(1:npft)
alparl(:)       = alparl_io(1:npft)
b_wl(:)         = b_wl_io(1:npft)
catch0(:)       = catch0_io(1:npft)
dcatch_dlai(:)  = dcatch_dlai_io(1:npft)
dgl_dm(:)       = dgl_dm_io(1:npft)
dgl_dt(:)       = dgl_dt_io(1:npft)
dqcrit(:)       = dqcrit_io(1:npft)
dz0v_dh(:)      = dz0v_dh_io(1:npft)
z0v(:)          = z0v_io(1:npft)
eta_sl(:)       = eta_sl_io(1:npft)
fd(:)           = fd_io(1:npft)
fsmc_of(:)      = fsmc_of_io(1:npft)
fsmc_p0(:)      = fsmc_p0_io(1:npft)
f0(:)           = f0_io(1:npft)
g_leaf_0(:)     = g_leaf_0_io(1:npft)
glmin(:)        = glmin_io(1:npft)
infil_f(:)      = infil_f_io(1:npft)
kext(:)         = kext_io(1:npft)
kpar(:)         = kpar_io(1:npft)
can_struct_a(:) = can_struct_a_io(1:npft)
gsoil_f(:)      = gsoil_f_io(1:npft)
lai_alb_lim(:)  = lai_alb_lim_io(1:npft)
neff(:)         = neff_io(1:npft)
nl0(:)          = nl0_io(1:npft)
nr_nl(:)        = nr_nl_io(1:npft)
ns_nl(:)        = ns_nl_io(1:npft)
omegau(:)       = omegau_io(1:npft)
omega(:)        = omega_io(1:npft)
omegal(:)       = omegal_io(1:npft)
omniru(:)       = omniru_io(1:npft)
omnir(:)        = omnir_io(1:npft)
omnirl(:)       = omnirl_io(1:npft)
r_grow(:)       = r_grow_io(1:npft)
rootd_ft(:)     = rootd_ft_io(1:npft)
sigl(:)         = sigl_io(1:npft)
tleaf_of(:)     = tleaf_of_io(1:npft)
tlow(:)         = tlow_io(1:npft)
tupp(:)         = tupp_io(1:npft)
emis_pft(:)     = emis_pft_io(1:npft)
dust_veg_scj(:) = dust_veg_scj_io(1:npft)
fl_o3_ct(:)     = fl_o3_ct_io(1:npft)
dfp_dcuo(:)     = dfp_dcuo_io(1:npft)
ci_st(:)        = ci_st_io(1:npft)
gpp_st(:)       = gpp_st_io(1:npft)
ief(:)          = ief_io(1:npft)
tef(:)          = tef_io(1:npft)
mef(:)          = mef_io(1:npft)
aef(:)          = aef_io(1:npft)

! INFERNO parameters
avg_ba(:)       = avg_ba_io(1:npft)
ccleaf_max(:)   = ccleaf_max_io(1:npft)
ccleaf_min(:)   = ccleaf_min_io(1:npft)
ccwood_max(:)   = ccwood_max_io(1:npft)
ccwood_min(:)   = ccwood_min_io(1:npft)
fef_bc(:)       = fef_bc_io(1:npft)
fef_ch4(:)      = fef_ch4_io(1:npft)
fef_co2(:)      = fef_co2_io(1:npft)
fef_co(:)       = fef_co_io(1:npft)
fef_nox(:)      = fef_nox_io(1:npft)
fef_oc(:)       = fef_oc_io(1:npft)
fef_so2(:)      = fef_so2_io(1:npft)

! Trait physiology parameters
lma(:)          = lma_io(1:npft)
nmass(:)        = nmass_io(1:npft)
vsl(:)          = vsl_io(1:npft)
vint(:)         = vint_io(1:npft)
kn(:)           = kn_io(1:npft)
q10_leaf(:)     = q10_leaf_io(1:npft)
hw_sw(:)        = hw_sw_io(1:npft)
nr(:)           = nr_io(1:npft)
nsw(:)          = nsw_io(1:npft)

z0h_z0m(1:npft)         = z0hm_pft_io(1:npft)
z0h_z0m_classic(1:npft) = z0hm_classic_pft_io(1:npft)

! jules_nvegparm
albsnc_nvg(:) = albsnc_nvg_io(1:nnvg)
albsnf_nvgu(:)= albsnf_nvgu_io(1:nnvg)
albsnf_nvg(:) = albsnf_nvg_io(1:nnvg)
albsnf_nvgl(:)= albsnf_nvgl_io(1:nnvg)
catch_nvg(:)  = catch_nvg_io(1:nnvg)
gs_nvg(:)     = gs_nvg_io(1:nnvg)
infil_nvg(:)  = infil_nvg_io(1:nnvg)
z0_nvg(:)     = z0_nvg_io(1:nnvg)
ch_nvg(:)     = ch_nvg_io(1:nnvg)
vf_nvg(:)     = vf_nvg_io(1:nnvg)
emis_nvg(:)   = emis_nvg_io(1:nnvg)
z0h_z0m(npft+1:npft + nnvg)         = z0hm_nvg_io(1:nnvg)
z0h_z0m_classic(npft+1:npft + nnvg) = z0hm_classic_nvg_io(1:nnvg)

! jules_triffid
! Space only allocated if at least phenology is enabled
IF ( l_triffid .OR. l_phenol ) THEN
  crop(:)    = crop_io(1:npft)
  g_area(:)  = g_area_io(1:npft)
  g_grow(:)  = g_grow_io(1:npft)
  g_root(:)  = g_root_io(1:npft)
  g_wood(:)  = g_wood_io(1:npft)
  lai_max(:) = lai_max_io(1:npft)
  lai_min(:) = lai_min_io(1:npft)
  alloc_fast(:) = alloc_fast_io(1:npft)
  alloc_med(:)  = alloc_med_io(1:npft)
  alloc_slow(:) = alloc_slow_io(1:npft)
  dpm_rpm_ratio(:) = dpm_rpm_ratio_io(1:npft)
  retran_l(:)      = retran_l_io(1:npft)
  retran_r(:)      = retran_r_io(1:npft)
END IF

! jules_irrig
IF ( l_irrig_dmd ) THEN
  ! irrtiles is now allocated with SIZE(npft) instead of SIZE(npft_max).
  ! It should really be allocated with SIZE(nirrtile) though, this will 
  ! be changed through ticket jules:#1065
  irrtiles(1:nirrtile) = irrigtiles(1:nirrtile)
END IF

! Now we have access to sm_levels from the dump, we can check the JULES_SOIL
! namelist for consistency
CALL check_jules_soil(sm_levels)


! jules_surf_param
! set diffuse fraction to 0.4 for CanRadMod 6
jules_vars%diff_frac(:) = 0.0
IF (can_rad_mod == 6) jules_vars%diff_frac(:) = 0.4

! jules_elevate
! Use same height for all land points.
IF ( l_aggregate ) THEN
  jules_vars%surf_hgt_surft(:,:) = 0.0
ELSE
  DO i = 1,ntiles
    jules_vars%surf_hgt_surft(:,i) = surf_hgt_io(i)
  END DO
END IF

  !----------------------------------------------------------------------
  ! Set bare soil roughness for use in 1 tile dust scheme
  !----------------------------------------------------------------------
  ! fix to enable CRUNs for two bin dust with 1 surface tile
  ! this mirrors code in sparm. Currently this is a simple real number
  ! in the future this may become an array and thus will need to be
  ! updated in line with the sparm code.

z0_soil = z0_nvg(soil - npft)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE surf_couple_allocate

END MODULE surf_couple_allocate_mod
#endif
