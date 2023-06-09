! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module holds surface parameters for each Plant Functional Type (but
! not parameters that are only used by TRIFFID).



! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

MODULE pftparm

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Radiation and albedo parameters.
!-----------------------------------------------------------------------------
INTEGER, ALLOCATABLE ::                                                       &
 orient(:)                                                                    &
                 ! Flag for leaf orientation: 1 for horizontal,
!                    0 for spherical.
,fsmc_mod(:)
                   ! Flag for whether water stress is calculated from
                   ! available water in layers weighted by root fraction (0)
                   ! or
                   ! whether water stress calculated from available
                   ! water in root zone (1)

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
 albsnc_max(:)                                                                &
                 ! Snow-covered albedo for large LAI.
,albsnc_min(:)                                                                &
                 ! Snow-covered albedo for zero LAI.
,albsnf_maxu(:)                                                               &
                 ! Max Snow-free albedo (max LAI) when scaled to obs
,albsnf_max(:)                                                                &
                 ! Snow-free albedo for large LAI.
,albsnf_maxl(:)                                                               &
                 ! Min Snow-free albedo (max LAI) when scaled to obs
,alniru(:)                                                                    &
                 ! upper limit on alnir, when scaled to albedo obs
,alnir(:)                                                                     &
                 ! Leaf reflection coefficient for near infra-red.
,alnirl(:)                                                                    &
                 ! lower limit on alnir, when scaled to albedo obs
,alparu(:)                                                                    &
                 ! upper limit on alpar, when scaled to albedo obs
,alpar(:)                                                                     &
                 ! Leaf reflection coefficient for PAR.
,alparl(:)                                                                    &
                 ! lower limit on alpar, when scaled to albedo obs
,kext(:)                                                                      &
                 ! Light extinction coefficient - used to
!                    calculate weightings for soil and veg.
  ,kpar(:)                                                                    &
                   ! PAR Extinction coefficient
!                    (m2 leaf/m2 ground)
  ,lai_alb_lim(:)                                                             &
!                  ! Lower limit on permitted LAI in albedo
  ,omegau(:)                                                                  &
                   ! upper limit on omega, when scaled to albedo obs
  ,omega(:)                                                                   &
                   ! Leaf scattering coefficient for PAR.
  ,omegal(:)                                                                  &
                   ! lower limit on omega, when scaled to albedo obs
  ,omniru(:)                                                                  &
                   ! upper limit on omnir, when scaled to albedo obs
  ,omnir(:)                                                                   &
                   ! Leaf scattering coefficient for near infra-red.
  ,omnirl(:)
                   ! lower limit on omnir, when scaled to albedo obs

!-----------------------------------------------------------------------------
! Parameters for phoyosynthesis and respiration.
!-----------------------------------------------------------------------------
INTEGER, ALLOCATABLE ::                                                       &
 c3(:)           ! Flag for C3 types: 1 for C3 Plants,
!                    0 for C4 Plants.

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
 alpha(:)                                                                     &
                 ! Quantum efficiency of photosynthesis
!                    (mol CO2/mol PAR photons).
  ,dqcrit(:)                                                                  &
                   ! Critical humidity deficit (kg H2O/kg air), used with the
                   ! Jacobs closure.
  ,fd(:)                                                                      &
                   ! Dark respiration coefficient.
  ,f0(:)                                                                      &
                   ! CI/CA for DQ = 0, used with the Jacobs closure.
  ,g1_stomata(:)                                                              &
                   ! Parameter g1 for the Medlyn et al. (2011) model of
                   ! stomatal conductance (kPa**0.5) - see Eqn.11 of
                   ! doi: 10.1111/j.1365-2486.2012.02790.x.
  ,neff(:)                                                                    &
                  ! Constant relating VCMAX and leaf N (mol/m2/s)
!                   from Schulze et al. 1994
!                   (AMAX = 0.4e-3 * NL  - assuming dry matter is
!                   40% carbon by mass)
!                   and Jacobs 1994:
!                   C3 : VCMAX = 2 * AMAX ;
!                   C4 : VCMAX = AMAX  ..
  ,nl0(:)                                                                     &
                   ! Top leaf nitrogen concentration
!                    (kg N/kg C).
  ,nr_nl(:)                                                                   &
                   ! Ratio of root nitrogen concentration to
!                    leaf nitrogen concentration.
  ,ns_nl(:)                                                                   &
                   ! Ratio of stem nitrogen concentration to
!                    leaf nitrogen concentration.
  ,nr(:)                                                                      &
                   ! Root nitrogen concentration (kg N/kg C)
  ,nsw(:)                                                                     &
                   ! Stem nitrogen concentration (kg N/kg C)
  ,hw_sw(:)                                                                   &
                   ! Heart:Stemwood Ratio (kg N/kg N)
  ,can_struct_a(:)                                                            &
                   ! Pinty canopy structure factor corresponding to overhead
                   ! sun (dimensionless)
  ,r_grow(:)                                                                  &
                   ! Growth respiration fraction.
  ,tlow(:)                                                                    &
                   ! Lower temperature for photosynthesis (deg C).
  ,tupp(:)                                                                    &
                   ! Upper temperature for photosynthesis (deg C).
  ,lma(:)                                                                     &
                   ! Leaf mass by area (1/SLA), (kg leaf/ m2)
  ,nmass(:)                                                                   &
                   ! Leaf nitrogen dry weight (g N/g leaf)
  ,vsl(:)                                                                     &
                   ! Slope of the Narea to Vcmax relationship
                   ! from Kattge et al. (2009)
  ,vint(:)                                                                    &
                   ! Y intercept of the Narea to Vcmax relationship
                   ! from Kattge et al. (2009)
  ,kn(:)                                                                      &
                   ! Exponential for N profile in canopy, used with
                   ! can_rad_mod=4, 5 (decay is a function of layers).
  ,knl(:)                                                                     &
                   ! Decay coefficient for N profile in canopy, used with
                   ! can_rad_mod=6 (decay is a function of LAI).
  ,q10_leaf(:)                                                                &
                   ! Factor for leaf respiration.
  !---------------------------------------------------------------------------
  ! Parameters for the Farquhar photosynthesis model.
  ! Jmax is the potential rate of electron transport.
  ! Vcmax is the maximum rate of carboxylation of Rubisco.
  !---------------------------------------------------------------------------
  ,act_jmax(:)                                                                &
                   ! Activation energy for temperature response of Jmax
                   ! (J mol-1).
  ,act_vcmax(:)                                                               &
                   ! Activation energy for temperature response of Vcmax
                   ! (J mol-1).
  ,alpha_elec(:)                                                              &
                   ! Quantum yield of electron transport
                   ! (mol electrons/mol PAR photons).
  ,deact_jmax(:)                                                              &
                   ! Deactivation energy for temperature response of Jmax
                   ! (J mol-1). This describes the rate of decrease
                   ! above the optimum temperature.
  ,deact_vcmax(:)                                                             &
                   ! Deactivation energy for temperature response of Vcmax
                   ! and Jmax (J mol-1). This describes the rate of decrease
                   ! above the optimum temperature.
  ,ds_jmax(:)                                                                 &
                   ! Entropy factor for temperature reponse of Jmax
                   ! (J mol-1 K-1).
  ,ds_vcmax(:)                                                                &
                   ! Entropy factor for temperature reponse of Vcmax
                   ! (J mol-1 K-1).
  ,jv25_ratio(:)
                   ! Ratio of Jmax to Vcmax at 25 deg C
                   ! (mol electrons mol-1 CO2).

!-----------------------------------------------------------------------------
! Allometric and other parameters.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
 a_wl(:)                                                                      &
                 ! Allometric coefficient relating the target
!                    woody biomass to the leaf area index
!                    (kg C/m2)
  ,a_ws(:)                                                                    &
                   ! Woody biomass as a multiple of live
!                    stem biomass.
  ,b_wl(:)                                                                    &
                   ! Allometric exponent relating the target
!                    woody biomass to the leaf area index.
  ,eta_sl(:)                                                                  &
                   ! Live stemwood coefficient (kg C m-1 (m2 leaf)-1)
  ,sigl(:)
                   ! Specific density of leaf carbon
!                    (kg C/m2 leaf).

!-----------------------------------------------------------------------------
! Phenology parameters.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
 g_leaf_0(:)                                                                  &
                 ! Minimum turnover rate for leaves (/360days).

,dgl_dm(:)                                                                    &
                 ! Rate of change of leaf turnover rate with
!                    moisture availability.
  ,fsmc_of(:)                                                                 &
                   ! Moisture availability below which leaves
!                    are dropped.
  ,dgl_dt(:)                                                                  &
                   ! Rate of change of leaf turnover rate with
!                    temperature (/K)
  ,tleaf_of(:)
                   ! Temperature below which leaves are
!                    dropped (K)

!-----------------------------------------------------------------------------
! Parameters for hydrological, thermal and other "physical" characteristics.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
 catch0(:)                                                                    &
                 ! Minimum canopy capacity (kg/m2).
,dcatch_dlai(:)                                                               &
                 ! Rate of change of canopy capacity with LAI.
,infil_f(:)                                                                   &
                 ! Infiltration enhancement factor.
,gsoil_f(:)                                                                   &
                 ! Soil evaporation enhancement factor (no units).
,glmin(:)                                                                     &
                 ! Minimum leaf conductance for H2O (m/s).
,dz0v_dh(:)                                                                   &
                 ! Rate of change of vegetation roughness
!                    length with height.
,z0v(:)                                                                       &
                 ! Specified vegetation roughness length.
!                    used with l_spec_veg_z0 = .true.
  ,rootd_ft(:)                                                                &
                   ! e-folding depth (m) of the root density.
  ,psi_close(:)                                                               &
                   ! soil matric potential (Pa) below which soil moisture
                   ! stress factor fsmc is zero. Should be negative.
  ,psi_open(:)                                                                &
                   ! soil matric potential (Pa) above which soil moisture
                   ! stress factor fsmc is one. Should be negative.
  ,fsmc_p0(:)                                                                 &
                   ! parameter in calculation of the
                   ! soil moisture at which the plant begins to experience
                   ! water stress
  ,emis_pft(:)                                                                &
                   !  Surface emissivity
  ,dust_veg_scj(:) ! Dust emission scaling factor for  each PFT



!-----------------------------------------------------------------------------
! Parameters for ozone damage
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
 fl_o3_ct(:)                                                                  &
                 ! Critical flux of O3 to vegetation (nmol/m2/s).
,dfp_dcuo(:)     ! Plant type specific O3 sensitivity parameter
                 ! (nmol-1 m2 s).

!-----------------------------------------------------------------------------
! Parameters for BVOC emissions
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
 ci_st(:),                                                                    &
                 ! Internal CO2 partial pressure (Pa)
                 !   at standard conditions
 gpp_st(:),                                                                   &
                 ! Gross primary productivity (KgC/m2/s)
                 !   at standard conditions
 ief(:),                                                                      &
                 ! Isoprene Emission Factor (ugC/g/h)
                 ! See Pacifico et al., (2011) Atm. Chem. Phys.
 tef(:),                                                                      &
                 ! (Mono-)Terpene Emission Factor (ugC/g/h)
 mef(:),                                                                      &
                 ! Methanol Emission Factor (ugC/g/h)
 aef(:)
                 ! Acetone Emission Factor (ugC/g/h)

!-----------------------------------------------------------------------------
! Parameters for INFERNO combustion
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
 ccleaf_min(:),                                                               &
                 ! Leaf minimum combustion completeness (kg/kg)
 ccleaf_max(:),                                                               &
                 ! Leaf maximum combustion completeness (kg/kg)
 ccwood_min(:),                                                               &
                 ! Wood (or Stem) minimum combustion completeness (kg/kg)
 ccwood_max(:),                                                               &
                 ! Wood (or Stem) maximum combustion completeness (kg/kg)
 avg_ba(:),                                                                   &
                 ! Average PFT Burnt Area per fire
 fire_mort(:)
                 ! Fire mortality per PFT

!-----------------------------------------------------------------------------
! Parameters for INFERNO emissions
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
 fef_co2(:),                                                                  &
                 ! Fire CO2 Emission Factor (g/kg)
                 ! See Thonicke et al., (2005,2010)
 fef_co(:),                                                                   &
                 ! Fire CO Emission Factor (g/kg)
 fef_ch4(:),                                                                  &
                 ! Fire CH4 Emission Factor (g/kg)
 fef_nox(:),                                                                  &
                 ! Fire NOx Emission Factor (g/kg)
 fef_so2(:),                                                                  &
                 ! Fire SO2 Emission Factor (g/kg)
 fef_oc(:),                                                                   &
                 ! Fire OC Emission Factor (g/kg)
 fef_bc(:)
                 ! Fire BC Emission Factor (g/kg)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PFTPARM'

CONTAINS
 
SUBROUTINE pftparm_alloc(npft)

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: npft

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PFTPARM_ALLOC'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!  ====pftparm module common====
! Veg surface type variables.
ALLOCATE( act_jmax(npft))
ALLOCATE( act_vcmax(npft))
ALLOCATE( albsnc_max(npft))
ALLOCATE( albsnc_min(npft))
ALLOCATE( albsnf_maxu(npft))
ALLOCATE( albsnf_max(npft))
ALLOCATE( albsnf_maxl(npft))
ALLOCATE( alpha_elec(npft))
ALLOCATE( alpha(npft))
ALLOCATE( alniru(npft))
ALLOCATE( alnir(npft))
ALLOCATE( alnirl(npft))
ALLOCATE( alparu(npft))
ALLOCATE( alpar(npft))
ALLOCATE( alparl(npft))
ALLOCATE( a_wl(npft))
ALLOCATE( a_ws(npft))
ALLOCATE( b_wl(npft))
ALLOCATE( catch0(npft))
ALLOCATE( c3(npft))
ALLOCATE( dcatch_dlai(npft))
ALLOCATE( deact_jmax(npft))
ALLOCATE( deact_vcmax(npft))
ALLOCATE( dgl_dm(npft))
ALLOCATE( dgl_dt(npft))
ALLOCATE( dqcrit(npft))
ALLOCATE( ds_jmax(npft))
ALLOCATE( ds_vcmax(npft))
ALLOCATE( dz0v_dh(npft))
ALLOCATE( z0v(npft))
ALLOCATE( emis_pft(npft))
ALLOCATE( eta_sl(npft))
ALLOCATE( fd(npft))
ALLOCATE( fsmc_of(npft))
ALLOCATE( f0(npft))
ALLOCATE( g1_stomata(npft))
ALLOCATE( glmin(npft))
ALLOCATE( g_leaf_0(npft))
ALLOCATE( infil_f(npft))
ALLOCATE( gsoil_f(npft))
ALLOCATE( jv25_ratio(npft))
ALLOCATE( kext(npft))
ALLOCATE( kpar(npft))
ALLOCATE( lai_alb_lim(npft))
ALLOCATE( neff(npft))
ALLOCATE( nl0(npft))
ALLOCATE( nr_nl(npft))
ALLOCATE( ns_nl(npft))
ALLOCATE( nsw(npft))
ALLOCATE( nr(npft))
ALLOCATE( hw_sw(npft))
ALLOCATE( can_struct_a(npft))
ALLOCATE( omegau(npft))
ALLOCATE( omega(npft))
ALLOCATE( omegal(npft))
ALLOCATE( omniru(npft))
ALLOCATE( omnir(npft))
ALLOCATE( omnirl(npft))
ALLOCATE( orient(npft))
ALLOCATE( r_grow(npft))
ALLOCATE( rootd_ft(npft))
ALLOCATE( psi_close(npft))
ALLOCATE( psi_open(npft))
ALLOCATE( fsmc_p0(npft))
ALLOCATE( fsmc_mod(npft))
ALLOCATE( sigl(npft))
ALLOCATE( tleaf_of(npft))
ALLOCATE( tlow(npft))
ALLOCATE( tupp(npft))
ALLOCATE( lma(npft))
ALLOCATE( nmass(npft))
ALLOCATE( vsl(npft))
ALLOCATE( vint(npft))
ALLOCATE( kn(npft))
ALLOCATE( knl(npft))
ALLOCATE( q10_leaf(npft))

act_jmax(:)     = 0.0
act_vcmax(:)    = 0.0
deact_jmax(:)   = 0.0
deact_vcmax(:)  = 0.0
albsnc_max(:)   = 0.0
albsnc_min(:)   = 0.0
albsnf_maxu(:)  = 0.0
albsnf_max(:)   = 0.0
albsnf_maxl(:)  = 0.0
alpha_elec(:)   = 0.0
alpha(:)        = 0.0
alniru(:)       = 0.0
alnir(:)        = 0.0
alnirl(:)       = 0.0
alparu(:)       = 0.0
alpar(:)        = 0.0
alparl(:)       = 0.0
a_wl(:)         = 0.0
a_ws(:)         = 0.0
b_wl(:)         = 0.0
catch0(:)       = 0.0
c3(:)           = 0
dcatch_dlai(:)  = 0.0
dgl_dm(:)       = 0.0
dgl_dt(:)       = 0.0
dqcrit(:)       = 0.0
ds_jmax(:)      = 0.0
ds_vcmax(:)     = 0.0
dz0v_dh(:)      = 0.0
z0v(:)          = 0.0
emis_pft(:)     = 0.0
eta_sl(:)       = 0.0
fd(:)           = 0.0
fsmc_of(:)      = 0.0
f0(:)           = 0.0
g1_stomata(:)   = 0.0
glmin(:)        = 0.0
g_leaf_0(:)     = 0.0
infil_f(:)      = 0.0
gsoil_f(:)      = 0.0
jv25_ratio(:)   = 0.0
kext(:)         = 0.0
kpar(:)         = 0.0
lai_alb_lim(:)  = 0.0
neff(:)         = 0.0
nl0(:)          = 0.0
nr_nl(:)        = 0.0
ns_nl(:)        = 0.0
nsw(:)          = 0.0
nr(:)           = 0.0
hw_sw(:)        = 0.0
can_struct_a(:) = 0.0
omegau(:)       = 0.0
omega(:)        = 0.0
omegal(:)       = 0.0
omniru(:)       = 0.0
omnir(:)        = 0.0
omnirl(:)       = 0.0
orient(:)       = 0
r_grow(:)       = 0.0
rootd_ft(:)     = 0.0
psi_close(:)    = 0.0
psi_open(:)     = 0.0
fsmc_p0(:)      = 0.0
fsmc_mod(:)     = 0
sigl(:)         = 0.0
tleaf_of(:)     = 0.0
tlow(:)         = 0.0
tupp(:)         = 0.0
lma(:)          = 0.0
nmass(:)        = 0.0
vsl(:)          = 0.0
vint(:)         = 0.0
kn(:)           = 0.0
knl(:)          = 0.0
q10_leaf(:)     = 0.0

! Ozone damage parameters
ALLOCATE( fl_o3_ct(npft))
ALLOCATE( dfp_dcuo(npft))

fl_o3_ct(:) = 0.0
dfp_dcuo(:) = 0.0

! BVOC emission parameters
ALLOCATE( ci_st(npft))
ALLOCATE( gpp_st(npft))
ALLOCATE( ief(npft))
ALLOCATE( tef(npft))
ALLOCATE( mef(npft))
ALLOCATE( aef(npft))

ci_st(:)  = 0.0
gpp_st(:) = 0.0
ief(:)    = 0.0
tef(:)    = 0.0
mef(:)    = 0.0
aef(:)    = 0.0

! INFERNO combustion parameters
ALLOCATE( ccleaf_min(npft))
ALLOCATE( ccleaf_max(npft))
ALLOCATE( ccwood_min(npft))
ALLOCATE( ccwood_max(npft))
ALLOCATE( avg_ba(npft))
ALLOCATE( fire_mort(npft))

ccleaf_min(:) = 0.0
ccleaf_max(:) = 0.0
ccwood_min(:) = 0.0
ccwood_max(:) = 0.0
avg_ba(:)     = 0.0
fire_mort(:)  = 0.0

! INFERNO emission parameters
ALLOCATE( fef_co2(npft))
ALLOCATE( fef_co(npft))
ALLOCATE( fef_ch4(npft))
ALLOCATE( fef_nox(npft))
ALLOCATE( fef_so2(npft))
ALLOCATE( fef_oc(npft))
ALLOCATE( fef_bc(npft))

fef_co2(:) = 0.0
fef_co(:)  = 0.0
fef_ch4(:) = 0.0
fef_nox(:) = 0.0
fef_so2(:) = 0.0
fef_oc(:)  = 0.0
fef_bc(:)  = 0.0

!UM only
#if defined(UM_JULES)
ALLOCATE( dust_veg_scj(npft))
dust_veg_scj(:) = 0.0
#endif

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE pftparm_alloc

END MODULE pftparm
