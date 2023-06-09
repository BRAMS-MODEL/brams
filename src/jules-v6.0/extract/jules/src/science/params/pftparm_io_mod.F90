! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This module contains variables used for reading in pftparm data
! and initialisations

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.


MODULE pftparm_io

USE max_dimensions, ONLY:                                                     &
  npft_max
USE missing_data_mod, ONLY: imdi, rmdi
USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!---------------------------------------------------------------------
! Set up variables to use in IO (a fixed size version of each array
! in pftparm that we want to initialise).
!---------------------------------------------------------------------
INTEGER ::                                                                    &
  c3_io(npft_max) = imdi,                                                     &
  orient_io(npft_max) = imdi

REAL(KIND=real_jlslsm) ::                                                     &
  a_wl_io(npft_max) = rmdi,                                                   &
  a_ws_io(npft_max) = rmdi,                                                   &
  act_jmax_io(npft_max) = rmdi,                                               &
  act_vcmax_io(npft_max) = rmdi,                                              &
  albsnc_max_io(npft_max) = rmdi,                                             &
  albsnc_min_io(npft_max) = rmdi,                                             &
  albsnf_maxu_io(npft_max) = rmdi,                                            &
  albsnf_max_io(npft_max) = rmdi,                                             &
  albsnf_maxl_io(npft_max) = rmdi,                                            &
  alpha_io(npft_max) = rmdi,                                                  &
  alpha_elec_io(npft_max) = rmdi,                                             &
  alniru_io(npft_max) = rmdi,                                                 &
  alnir_io(npft_max) = rmdi,                                                  &
  alnirl_io(npft_max) = rmdi,                                                 &
  alparu_io(npft_max) = rmdi,                                                 &
  alpar_io(npft_max) = rmdi,                                                  &
  alparl_io(npft_max) = rmdi,                                                 &
  b_wl_io(npft_max) = rmdi,                                                   &
  catch0_io(npft_max) = rmdi,                                                 &
  dcatch_dlai_io(npft_max) = rmdi,                                            &
  deact_jmax_io(npft_max) = rmdi,                                             &
  deact_vcmax_io(npft_max) = rmdi,                                            &
  dgl_dm_io(npft_max) = rmdi,                                                 &
  dgl_dt_io(npft_max) = rmdi,                                                 &
  dqcrit_io(npft_max) = rmdi,                                                 &
  ds_jmax_io(npft_max) = rmdi,                                                &
  ds_vcmax_io(npft_max) = rmdi,                                               &
  dz0v_dh_io(npft_max) = rmdi,                                                &
  z0v_io(npft_max) = rmdi,                                                    &
  eta_sl_io(npft_max) = rmdi,                                                 &
  fd_io(npft_max) = rmdi,                                                     &
  fsmc_of_io(npft_max) = rmdi,                                                &
  fsmc_p0_io(npft_max) = rmdi,                                                &
  f0_io(npft_max) = rmdi,                                                     &
  g1_stomata_io(npft_max) = rmdi,                                             &
  g_leaf_0_io(npft_max) = rmdi,                                               &
  glmin_io(npft_max) = rmdi,                                                  &
  infil_f_io(npft_max) = rmdi,                                                &
  gsoil_f_io(npft_max) = rmdi,                                                &
  jv25_ratio_io(npft_max) = rmdi,                                             &
  kext_io(npft_max) = rmdi,                                                   &
  kpar_io(npft_max) = rmdi,                                                   &
  lai_alb_lim_io(npft_max) = rmdi,                                            &
  neff_io(npft_max) = rmdi,                                                   &
  nl0_io(npft_max) = rmdi,                                                    &
  nr_nl_io(npft_max) = rmdi,                                                  &
  ns_nl_io(npft_max) = rmdi,                                                  &
  nsw_io(npft_max) = rmdi,                                                    &
  nr_io(npft_max) = rmdi,                                                     &
  hw_sw_io(npft_max) = rmdi,                                                  &
  can_struct_a_io(npft_max) = rmdi,                                           &
  omegau_io(npft_max) = rmdi,                                                 &
  omega_io(npft_max) = rmdi,                                                  &
  omegal_io(npft_max) = rmdi,                                                 &
  omniru_io(npft_max) = rmdi,                                                 &
  omnir_io(npft_max) = rmdi,                                                  &
  omnirl_io(npft_max) = rmdi,                                                 &
  r_grow_io(npft_max) = rmdi,                                                 &
  rootd_ft_io(npft_max) = rmdi,                                               &
  sigl_io(npft_max) = rmdi,                                                   &
  tleaf_of_io(npft_max) = rmdi,                                               &
  tlow_io(npft_max) = rmdi,                                                   &
  tupp_io(npft_max) = rmdi,                                                   &
  emis_pft_io(npft_max) = rmdi,                                               &
  z0hm_pft_io(npft_max) = rmdi,                                               &
  z0hm_classic_pft_io(npft_max) = rmdi,                                       &
  dust_veg_scj_io(npft_max) = rmdi,                                           &
  fl_o3_ct_io(npft_max) = rmdi,                                               &
  dfp_dcuo_io(npft_max) = rmdi,                                               &
  ci_st_io(npft_max) = rmdi,                                                  &
  gpp_st_io(npft_max) = rmdi,                                                 &
  ief_io(npft_max) = rmdi,                                                    &
  tef_io(npft_max) = rmdi,                                                    &
  mef_io(npft_max) = rmdi,                                                    &
  aef_io(npft_max) = rmdi,                                                    &
  fef_co2_io(npft_max) = rmdi,                                                &
  fef_co_io(npft_max) = rmdi,                                                 &
  fef_ch4_io(npft_max) = rmdi,                                                &
  fef_nox_io(npft_max) = rmdi,                                                &
  fef_so2_io(npft_max) = rmdi,                                                &
  fef_oc_io(npft_max) = rmdi,                                                 &
  fef_bc_io(npft_max) = rmdi,                                                 &
  ccleaf_min_io(npft_max) = rmdi,                                             &
  ccleaf_max_io(npft_max) = rmdi,                                             &
  ccwood_min_io(npft_max) = rmdi,                                             &
  ccwood_max_io(npft_max) = rmdi,                                             &
  avg_ba_io(npft_max) = rmdi,                                                 &
  lma_io(npft_max) = rmdi,                                                    &
  nmass_io(npft_max) = rmdi,                                                  &
  vsl_io(npft_max) = rmdi,                                                    &
  vint_io(npft_max) = rmdi,                                                   &
  kn_io(npft_max) = rmdi,                                                     &
  knl_io(npft_max) = rmdi,                                                    &
  q10_leaf_io(npft_max) = rmdi,                                               &
  fire_mort_io(npft_max) = rmdi

#if !defined(UM_JULES)
INTEGER ::                                                                    &
  fsmc_mod_io(npft_max) = imdi

REAL(KIND=real_jlslsm) ::                                                     &
  canht_ft_io(npft_max) = rmdi,                                               &
  lai_io(npft_max) = rmdi,                                                    &
  psi_close_io(npft_max) = rmdi,                                              &
  psi_open_io(npft_max) = rmdi
#endif

!---------------------------------------------------------------------
! Set up a namelist for reading and writing these arrays
!---------------------------------------------------------------------
NAMELIST  / jules_pftparm/                                                    &
#if !defined(UM_JULES)
                         fsmc_mod_io, canht_ft_io, lai_io,                    &
                         psi_close_io, psi_open_io,                           &
#endif
                         c3_io,orient_io,a_wl_io,a_ws_io,                     &
                         act_jmax_io,act_vcmax_io,                            &
                         albsnc_max_io,albsnc_min_io,albsnf_maxu_io,          &
                         albsnf_max_io,albsnf_maxl_io,                        &
                         alpha_io,alpha_elec_io,alniru_io,alnir_io,           &
                         alnirl_io,alparu_io,alpar_io,alparl_io,b_wl_io,      &
                         catch0_io,dcatch_dlai_io,deact_jmax_io,              &
                         deact_vcmax_io,dgl_dm_io,                            &
                         dgl_dt_io,dqcrit_io,ds_jmax_io,                      &
                         ds_vcmax_io,dz0v_dh_io,z0v_io,eta_sl_io,             &
                         fd_io,fsmc_of_io,fsmc_p0_io,f0_io,g1_stomata_io,     &
                         g_leaf_0_io, glmin_io,infil_f_io,jv25_ratio_io,      &
                         kext_io,kpar_io,lai_alb_lim_io,                      &
                         neff_io,nl0_io,nr_nl_io,ns_nl_io,                    &
                         nr_io,nsw_io,hw_sw_io,                               &
                         can_struct_a_io,gsoil_f_io,                          &
                         omegau_io,omega_io,omegal_io,omniru_io,              &
                         omnir_io,omnirl_io,r_grow_io,rootd_ft_io,            &
                         sigl_io,tleaf_of_io,tlow_io,tupp_io,                 &
                         emis_pft_io,z0hm_pft_io,z0hm_classic_pft_io,         &
                         dust_veg_scj_io,fl_o3_ct_io,dfp_dcuo_io,             &
                         ci_st_io,gpp_st_io,                                  &
                         ief_io,tef_io,mef_io,aef_io,                         &
                         fef_co2_io, fef_co_io, fef_ch4_io,                   &
                         fef_nox_io, fef_so2_io, fef_oc_io,                   &
                         fef_bc_io,ccleaf_min_io, ccleaf_max_io,              &
                         ccwood_min_io, ccwood_max_io, avg_ba_io,             &
                         lma_io,nmass_io,vsl_io,vint_io,kn_io,knl_io,         &
                         q10_leaf_io, fire_mort_io
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PFTPARM_IO'

CONTAINS

SUBROUTINE print_nlist_jules_pftparm()
USE jules_print_mgr, ONLY: jules_print
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer

CALL jules_print('pftparm_io',                                                &
    'Contents of namelist jules_pftparm')

#if !defined(UM_JULES)
WRITE(lineBuffer,*)' canht_ft_io = ',canht_ft_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' lai_io = ',lai_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' fsmc_mod_io = ',fsmc_mod_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' psi_close_io = ',psi_close_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' psi_open_io = ',psi_open_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' fef_co2_io = ',fef_co2_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' fef_co_io = ',fef_co_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' fef_ch4_io = ',fef_ch4_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' fef_nox_io = ',fef_nox_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' fef_so2_io = ',fef_so2_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' fef_oc_io = ',fef_oc_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' fef_bc_io = ',fef_bc_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' ccleaf_min_io = ',ccleaf_min_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' ccleaf_max_io = ',ccleaf_max_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' ccwood_min_io = ',ccwood_min_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' ccwood_max_io = ',ccwood_max_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' avg_ba_io = ',avg_ba_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' fire_mort_io = ',fire_mort_io
CALL jules_print('pftparm_io',lineBuffer)
#endif
WRITE(lineBuffer,*)' c3_io = ',c3_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' orient_io = ',orient_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' a_wl_io = ',a_wl_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' a_ws_io = ',a_ws_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' act_jmax_io = ',act_jmax_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' act_vcmax_io = ',act_vcmax_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' albsnc_max_io = ',albsnc_max_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' albsnc_min_io = ',albsnc_min_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' albsnf_max_io = ',albsnf_max_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' alpha_elec_io = ',alpha_elec_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' alpha_io = ',alpha_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' alnir_io = ',alnir_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' alpar_io = ',alpar_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' b_wl_io = ',b_wl_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' catch0_io = ',catch0_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' dcatch_dlai_io = ',dcatch_dlai_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' deact_vcmax_io = ',deact_jmax_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' deact_vcmax_io = ',deact_vcmax_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' dgl_dm_io = ',dgl_dm_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' dgl_dt_io = ',dgl_dt_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' dqcrit_io = ',dqcrit_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' ds_jmax_io = ',ds_jmax_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' ds_vcmax_io = ',ds_vcmax_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' dz0v_dh_io = ',dz0v_dh_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' z0v_io = ',z0v_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' eta_sl_io = ',eta_sl_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' fd_io = ',fd_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' fsmc_of_io = ',fsmc_of_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' fsmc_p0_io = ',fsmc_p0_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' f0_io = ',f0_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' g1_stomata_io = ',g1_stomata_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' g_leaf_0_io = ',g_leaf_0_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' glmin_io = ',glmin_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' infil_f_io = ',infil_f_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' gsoil_f_io = ',gsoil_f_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' jv25_ratio_io = ',jv25_ratio_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' kext_io = ',kext_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' kpar_io = ',kpar_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' lai_alb_lim_io = ',lai_alb_lim_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' neff_io = ',neff_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' nl0_io = ',nl0_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' nr_nl_io = ',nr_nl_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' ns_nl_io = ',ns_nl_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' nsw_io = ',nsw_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' nr_io = ',nr_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' hw_sw_io = ',hw_sw_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' can_struct_a_io = ',can_struct_a_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' omega_io = ',omega_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' omnir_io = ',omnir_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' r_grow_io = ',r_grow_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' rootd_ft_io = ',rootd_ft_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' sigl_io = ',sigl_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' tleaf_of_io = ',tleaf_of_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' tlow_io = ',tlow_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' tupp_io = ',tupp_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' emis_pft_io = ',emis_pft_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' z0hm_pft_io = ',z0hm_pft_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' dust_veg_scj_io = ',dust_veg_scj_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' fl_o3_ct_io = ',fl_o3_ct_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' dfp_dcuo_io = ',dfp_dcuo_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' ci_st_io = ',ci_st_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' gpp_st_io = ',gpp_st_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' ief_io = ',ief_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' tef_io = ',tef_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' mef_io = ',mef_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' aef_io = ',aef_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' avg_ba_io = ',avg_ba_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' ccleaf_max_io = ',ccleaf_max_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' ccleaf_min_io = ',ccleaf_min_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' ccwood_max_io = ',ccwood_max_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' ccwood_min_io = ',ccwood_min_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' fef_bc_io = ',fef_bc_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' fef_ch4_io = ',fef_ch4_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' fef_co2_io = ',fef_co2_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' fef_co_io = ',fef_co_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' fef_nox_io = ',fef_nox_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' fef_oc_io = ',fef_oc_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' fef_so2_io = ',fef_so2_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' lma_io = ',lma_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' nmass_io = ',nmass_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' vsl_io = ',vsl_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' vint_io = ',vint_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' kn_io = ',kn_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' knl_io = ',knl_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' q10_leaf_io = ',q10_leaf_io
CALL jules_print('pftparm_io',lineBuffer)
WRITE(lineBuffer,*)' fire_mort_io = ',fire_mort_io
CALL jules_print('pftparm_io',lineBuffer)

CALL jules_print('pftparm_io',                                                &
    '- - - - - - end of namelist - - - - - -')

END SUBROUTINE print_nlist_jules_pftparm

#if defined(UM_JULES)
SUBROUTINE read_nml_jules_pftparm (unitnumber)

! Description:
!  Read the JULES_PFTPARM namelist

USE setup_namelist, ONLY: setup_nml_type
USE check_iostat_mod, ONLY:  check_iostat
USE UM_parcore,       ONLY:  mype
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: unitnumber
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_JULES_PFTPARM'
INTEGER(KIND=jpim), PARAMETER          :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER          :: zhook_out = 1
CHARACTER(LEN=errormessagelength) :: iomessage

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_int = 2 * npft_max
INTEGER, PARAMETER :: n_real = 95 * npft_max

TYPE my_namelist
  SEQUENCE
  INTEGER :: c3_io(npft_max)
  INTEGER :: orient_io(npft_max)
  REAL(KIND=real_jlslsm) :: a_wl_io(npft_max)
  REAL(KIND=real_jlslsm) :: a_ws_io(npft_max)
  REAL(KIND=real_jlslsm) :: act_jmax_io(npft_max)
  REAL(KIND=real_jlslsm) :: act_vcmax_io(npft_max)
  REAL(KIND=real_jlslsm) :: albsnc_max_io(npft_max)
  REAL(KIND=real_jlslsm) :: albsnc_min_io(npft_max)
  REAL(KIND=real_jlslsm) :: albsnf_maxu_io(npft_max)
  REAL(KIND=real_jlslsm) :: albsnf_max_io(npft_max)
  REAL(KIND=real_jlslsm) :: albsnf_maxl_io(npft_max)
  REAL(KIND=real_jlslsm) :: alpha_elec_io(npft_max)
  REAL(KIND=real_jlslsm) :: alpha_io(npft_max)
  REAL(KIND=real_jlslsm) :: alniru_io(npft_max)
  REAL(KIND=real_jlslsm) :: alnir_io(npft_max)
  REAL(KIND=real_jlslsm) :: alnirl_io(npft_max)
  REAL(KIND=real_jlslsm) :: alparu_io(npft_max)
  REAL(KIND=real_jlslsm) :: alpar_io(npft_max)
  REAL(KIND=real_jlslsm) :: alparl_io(npft_max)
  REAL(KIND=real_jlslsm) :: b_wl_io(npft_max)
  REAL(KIND=real_jlslsm) :: catch0_io(npft_max)
  REAL(KIND=real_jlslsm) :: dcatch_dlai_io(npft_max)
  REAL(KIND=real_jlslsm) :: deact_jmax_io(npft_max)
  REAL(KIND=real_jlslsm) :: deact_vcmax_io(npft_max)
  REAL(KIND=real_jlslsm) :: dgl_dm_io(npft_max)
  REAL(KIND=real_jlslsm) :: dgl_dt_io(npft_max)
  REAL(KIND=real_jlslsm) :: dqcrit_io(npft_max)
  REAL(KIND=real_jlslsm) :: ds_jmax_io(npft_max)
  REAL(KIND=real_jlslsm) :: ds_vcmax_io(npft_max)
  REAL(KIND=real_jlslsm) :: dz0v_dh_io(npft_max)
  REAL(KIND=real_jlslsm) :: z0v_io(npft_max)
  REAL(KIND=real_jlslsm) :: eta_sl_io(npft_max)
  REAL(KIND=real_jlslsm) :: fd_io(npft_max)
  REAL(KIND=real_jlslsm) :: fsmc_of_io(npft_max)
  REAL(KIND=real_jlslsm) :: fsmc_p0_io(npft_max)
  REAL(KIND=real_jlslsm) :: f0_io(npft_max)
  REAL(KIND=real_jlslsm) :: g1_stomata_io(npft_max)
  REAL(KIND=real_jlslsm) :: g_leaf_0_io(npft_max)
  REAL(KIND=real_jlslsm) :: glmin_io(npft_max)
  REAL(KIND=real_jlslsm) :: infil_f_io(npft_max)
  REAL(KIND=real_jlslsm) :: gsoil_f_io(npft_max)
  REAL(KIND=real_jlslsm) :: jv25_ratio_io(npft_max)
  REAL(KIND=real_jlslsm) :: kext_io(npft_max)
  REAL(KIND=real_jlslsm) :: kpar_io(npft_max)
  REAL(KIND=real_jlslsm) :: lai_alb_lim_io(npft_max)
  REAL(KIND=real_jlslsm) :: neff_io(npft_max)
  REAL(KIND=real_jlslsm) :: nl0_io(npft_max)
  REAL(KIND=real_jlslsm) :: nr_nl_io(npft_max)
  REAL(KIND=real_jlslsm) :: ns_nl_io(npft_max)
  REAL(KIND=real_jlslsm) :: nsw_io(npft_max)
  REAL(KIND=real_jlslsm) :: nr_io(npft_max)
  REAL(KIND=real_jlslsm) :: hw_sw_io(npft_max)
  REAL(KIND=real_jlslsm) :: can_struct_a_io(npft_max)
  REAL(KIND=real_jlslsm) :: omegau_io(npft_max)
  REAL(KIND=real_jlslsm) :: omega_io(npft_max)
  REAL(KIND=real_jlslsm) :: omegal_io(npft_max)
  REAL(KIND=real_jlslsm) :: omniru_io(npft_max)
  REAL(KIND=real_jlslsm) :: omnir_io(npft_max)
  REAL(KIND=real_jlslsm) :: omnirl_io(npft_max)
  REAL(KIND=real_jlslsm) :: r_grow_io(npft_max)
  REAL(KIND=real_jlslsm) :: rootd_ft_io(npft_max)
  REAL(KIND=real_jlslsm) :: sigl_io(npft_max)
  REAL(KIND=real_jlslsm) :: tleaf_of_io(npft_max)
  REAL(KIND=real_jlslsm) :: tlow_io(npft_max)
  REAL(KIND=real_jlslsm) :: tupp_io(npft_max)
  REAL(KIND=real_jlslsm) :: emis_pft_io(npft_max)
  REAL(KIND=real_jlslsm) :: z0hm_pft_io(npft_max)
  REAL(KIND=real_jlslsm) :: z0hm_classic_pft_io(npft_max)
  REAL(KIND=real_jlslsm) :: dust_veg_scj_io(npft_max)
  REAL(KIND=real_jlslsm) :: fl_o3_ct_io(npft_max)
  REAL(KIND=real_jlslsm) :: dfp_dcuo_io(npft_max)
  REAL(KIND=real_jlslsm) :: ief_io(npft_max)
  REAL(KIND=real_jlslsm) :: tef_io(npft_max)
  REAL(KIND=real_jlslsm) :: mef_io(npft_max)
  REAL(KIND=real_jlslsm) :: aef_io(npft_max)
  REAL(KIND=real_jlslsm) :: ci_st_io(npft_max)
  REAL(KIND=real_jlslsm) :: gpp_st_io(npft_max)
  REAL(KIND=real_jlslsm) :: fef_co2_io(npft_max)
  REAL(KIND=real_jlslsm) :: fef_co_io(npft_max)
  REAL(KIND=real_jlslsm) :: fef_ch4_io(npft_max)
  REAL(KIND=real_jlslsm) :: fef_nox_io(npft_max)
  REAL(KIND=real_jlslsm) :: fef_so2_io(npft_max)
  REAL(KIND=real_jlslsm) :: fef_oc_io(npft_max)
  REAL(KIND=real_jlslsm) :: fef_bc_io(npft_max)
  REAL(KIND=real_jlslsm) :: ccleaf_min_io(npft_max)
  REAL(KIND=real_jlslsm) :: ccleaf_max_io(npft_max)
  REAL(KIND=real_jlslsm) :: ccwood_min_io(npft_max)
  REAL(KIND=real_jlslsm) :: ccwood_max_io(npft_max)
  REAL(KIND=real_jlslsm) :: avg_ba_io(npft_max)
  REAL(KIND=real_jlslsm) :: lma_io(npft_max)
  REAL(KIND=real_jlslsm) :: nmass_io(npft_max)
  REAL(KIND=real_jlslsm) :: vsl_io(npft_max)
  REAL(KIND=real_jlslsm) :: vint_io(npft_max)
  REAL(KIND=real_jlslsm) :: kn_io(npft_max)
  REAL(KIND=real_jlslsm) :: knl_io(npft_max)
  REAL(KIND=real_jlslsm) :: q10_leaf_io(npft_max)
  REAL(KIND=real_jlslsm) :: fire_mort_io(npft_max)
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in = n_int,              &
                    n_real_in = n_real)

IF (mype == 0) THEN

  READ (UNIT = unitnumber, NML = jules_pftparm, IOSTAT = errorstatus,         &
        IOMSG = iomessage)
  CALL check_iostat(errorstatus, "namelist jules_pftparm",                    &
           iomessage)

  my_nml % c3_io          = c3_io
  my_nml % orient_io      = orient_io
  my_nml % a_wl_io        = a_wl_io
  my_nml % a_ws_io        = a_ws_io
  my_nml % act_jmax_io    = act_jmax_io
  my_nml % act_vcmax_io   = act_vcmax_io
  my_nml % albsnc_max_io  = albsnc_max_io
  my_nml % albsnc_min_io  = albsnc_min_io
  my_nml % albsnf_maxu_io = albsnf_maxu_io
  my_nml % albsnf_max_io  = albsnf_max_io
  my_nml % albsnf_maxl_io = albsnf_maxl_io
  my_nml % alpha_elec_io  = alpha_elec_io
  my_nml % alpha_io       = alpha_io
  my_nml % alniru_io      = alniru_io
  my_nml % alnir_io       = alnir_io
  my_nml % alnirl_io      = alnirl_io
  my_nml % alparu_io      = alparu_io
  my_nml % alpar_io       = alpar_io
  my_nml % alparl_io      = alparl_io
  my_nml % b_wl_io        = b_wl_io
  my_nml % catch0_io      = catch0_io
  my_nml % dcatch_dlai_io = dcatch_dlai_io
  my_nml % deact_jmax_io  = deact_jmax_io
  my_nml % deact_vcmax_io = deact_vcmax_io
  my_nml % dgl_dm_io      = dgl_dm_io
  my_nml % dgl_dt_io      = dgl_dt_io
  my_nml % dqcrit_io      = dqcrit_io
  my_nml % ds_jmax_io     = ds_jmax_io
  my_nml % ds_vcmax_io    = ds_vcmax_io
  my_nml % dz0v_dh_io     = dz0v_dh_io
  my_nml % z0v_io         = z0v_io
  my_nml % eta_sl_io      = eta_sl_io
  my_nml % fd_io          = fd_io
  my_nml % fsmc_of_io     = fsmc_of_io
  my_nml % fsmc_p0_io     = fsmc_p0_io
  my_nml % f0_io          = f0_io
  my_nml % g1_stomata_io  = g1_stomata_io
  my_nml % g_leaf_0_io    = g_leaf_0_io
  my_nml % glmin_io       = glmin_io
  my_nml % infil_f_io     = infil_f_io
  my_nml % gsoil_f_io     = gsoil_f_io
  my_nml % jv25_ratio_io  = jv25_ratio_io
  my_nml % kext_io        = kext_io
  my_nml % kpar_io        = kpar_io
  my_nml % lai_alb_lim_io = lai_alb_lim_io
  my_nml % neff_io        = neff_io
  my_nml % nl0_io         = nl0_io
  my_nml % nr_nl_io       = nr_nl_io
  my_nml % ns_nl_io       = ns_nl_io
  my_nml % nsw_io         = nsw_io
  my_nml % nr_io          = nr_io
  my_nml % hw_sw_io       = hw_sw_io
  my_nml % can_struct_a_io = can_struct_a_io
  my_nml % omegau_io      = omegau_io
  my_nml % omega_io       = omega_io
  my_nml % omegal_io      = omegal_io
  my_nml % omniru_io      = omniru_io
  my_nml % omnir_io       = omnir_io
  my_nml % omnirl_io      = omnirl_io
  my_nml % r_grow_io      = r_grow_io
  my_nml % rootd_ft_io    = rootd_ft_io
  my_nml % sigl_io        = sigl_io
  my_nml % tleaf_of_io    = tleaf_of_io
  my_nml % tlow_io        = tlow_io
  my_nml % tupp_io        = tupp_io
  my_nml % emis_pft_io    = emis_pft_io
  my_nml % z0hm_pft_io    = z0hm_pft_io
  my_nml % z0hm_classic_pft_io = z0hm_classic_pft_io
  my_nml % dust_veg_scj_io = dust_veg_scj_io
  my_nml % fl_o3_ct_io    = fl_o3_ct_io
  my_nml % dfp_dcuo_io    = dfp_dcuo_io
  my_nml % ief_io         = ief_io
  my_nml % tef_io         = tef_io
  my_nml % mef_io         = mef_io
  my_nml % aef_io         = aef_io
  my_nml % ci_st_io       = ci_st_io
  my_nml % gpp_st_io      = gpp_st_io
  my_nml % fef_co2_io     = fef_co2_io
  my_nml % fef_co_io      = fef_co_io
  my_nml % fef_ch4_io     = fef_ch4_io
  my_nml % fef_nox_io     = fef_nox_io
  my_nml % fef_so2_io     = fef_so2_io
  my_nml % fef_oc_io      = fef_oc_io
  my_nml % fef_bc_io      = fef_bc_io
  my_nml % ccleaf_min_io  = ccleaf_min_io
  my_nml % ccleaf_max_io  = ccleaf_max_io
  my_nml % ccwood_min_io  = ccwood_min_io
  my_nml % ccwood_max_io  = ccwood_max_io
  my_nml % avg_ba_io      = avg_ba_io
  my_nml % lma_io         = lma_io
  my_nml % nmass_io       = nmass_io
  my_nml % vsl_io         = vsl_io
  my_nml % vint_io        = vint_io
  my_nml % kn_io          = kn_io
  my_nml % knl_io         = knl_io
  my_nml % q10_leaf_io    = q10_leaf_io
  my_nml % fire_mort_io   = fire_mort_io
END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  c3_io           = my_nml % c3_io
  orient_io       = my_nml % orient_io
  a_wl_io         = my_nml % a_wl_io
  a_ws_io         = my_nml % a_ws_io
  act_jmax_io     = my_nml % act_jmax_io
  act_vcmax_io    = my_nml % act_vcmax_io
  albsnc_max_io   = my_nml % albsnc_max_io
  albsnc_min_io   = my_nml % albsnc_min_io
  albsnf_maxu_io  = my_nml % albsnf_maxu_io
  albsnf_max_io   = my_nml % albsnf_max_io
  albsnf_maxl_io  = my_nml % albsnf_maxl_io
  alpha_elec_io   = my_nml % alpha_elec_io
  alpha_io        = my_nml % alpha_io
  alniru_io       = my_nml % alniru_io
  alnir_io        = my_nml % alnir_io
  alnirl_io       = my_nml % alnirl_io
  alparu_io       = my_nml % alparu_io
  alpar_io        = my_nml % alpar_io
  alparl_io       = my_nml % alparl_io
  b_wl_io         = my_nml % b_wl_io
  catch0_io       = my_nml % catch0_io
  dcatch_dlai_io  = my_nml % dcatch_dlai_io
  deact_jmax_io   = my_nml % deact_jmax_io
  deact_vcmax_io  = my_nml % deact_vcmax_io
  dgl_dm_io       = my_nml % dgl_dm_io
  dgl_dt_io       = my_nml % dgl_dt_io
  dqcrit_io       = my_nml % dqcrit_io
  ds_jmax_io      = my_nml % ds_jmax_io
  ds_vcmax_io     = my_nml % ds_vcmax_io
  dz0v_dh_io      = my_nml % dz0v_dh_io
  z0v_io          = my_nml % z0v_io
  eta_sl_io       = my_nml % eta_sl_io
  fd_io           = my_nml % fd_io
  fsmc_of_io      = my_nml % fsmc_of_io
  fsmc_p0_io      = my_nml % fsmc_p0_io
  f0_io           = my_nml % f0_io
  g1_stomata_io   = my_nml % g1_stomata_io
  g_leaf_0_io     = my_nml % g_leaf_0_io
  glmin_io        = my_nml % glmin_io
  infil_f_io      = my_nml % infil_f_io
  gsoil_f_io      = my_nml % gsoil_f_io
  jv25_ratio_io   = my_nml % jv25_ratio_io
  kext_io         = my_nml % kext_io
  kpar_io         = my_nml % kpar_io
  lai_alb_lim_io  = my_nml % lai_alb_lim_io
  neff_io         = my_nml % neff_io
  nl0_io          = my_nml % nl0_io
  nr_nl_io        = my_nml % nr_nl_io
  ns_nl_io        = my_nml % ns_nl_io
  nsw_io          = my_nml % nsw_io
  nr_io           = my_nml % nr_io
  hw_sw_io        = my_nml % hw_sw_io
  can_struct_a_io = my_nml % can_struct_a_io
  omegau_io       = my_nml % omegau_io
  omega_io        = my_nml % omega_io
  omegal_io       = my_nml % omegal_io
  omniru_io       = my_nml % omniru_io
  omnir_io        = my_nml % omnir_io
  omnirl_io       = my_nml % omnirl_io
  r_grow_io       = my_nml % r_grow_io
  rootd_ft_io     = my_nml % rootd_ft_io
  sigl_io         = my_nml % sigl_io
  tleaf_of_io     = my_nml % tleaf_of_io
  tlow_io         = my_nml % tlow_io
  tupp_io         = my_nml % tupp_io
  emis_pft_io     = my_nml % emis_pft_io
  z0hm_pft_io     = my_nml % z0hm_pft_io
  z0hm_classic_pft_io = my_nml % z0hm_classic_pft_io
  dust_veg_scj_io = my_nml % dust_veg_scj_io
  fl_o3_ct_io     = my_nml % fl_o3_ct_io
  dfp_dcuo_io     = my_nml % dfp_dcuo_io
  ief_io          = my_nml % ief_io
  tef_io          = my_nml % tef_io
  mef_io          = my_nml % mef_io
  aef_io          = my_nml % aef_io
  ci_st_io        = my_nml % ci_st_io
  gpp_st_io       = my_nml % gpp_st_io
  fef_co2_io      = my_nml % fef_co2_io
  fef_co_io       = my_nml % fef_co_io
  fef_ch4_io      = my_nml % fef_ch4_io
  fef_nox_io      = my_nml % fef_nox_io
  fef_so2_io      = my_nml % fef_so2_io
  fef_oc_io       = my_nml % fef_oc_io
  fef_bc_io       = my_nml % fef_bc_io
  ccleaf_min_io   = my_nml % ccleaf_min_io
  ccleaf_max_io   = my_nml % ccleaf_max_io
  ccwood_min_io   = my_nml % ccwood_min_io
  ccwood_max_io   = my_nml % ccwood_max_io
  avg_ba_io       = my_nml % avg_ba_io
  lma_io          = my_nml % lma_io
  nmass_io        = my_nml % nmass_io
  vsl_io          = my_nml % vsl_io
  vint_io         = my_nml % vint_io
  kn_io           = my_nml % kn_io
  knl_io          = my_nml % knl_io
  q10_leaf_io     = my_nml % q10_leaf_io
  fire_mort_io    = my_nml % fire_mort_io
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_nml_jules_pftparm
#endif

END MODULE pftparm_io
