!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT**************************************

MODULE jules_soil_biogeochem_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Contains soil biogeochemistry options and parameters, and a namelist for
!   setting them.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Public scope by default.

!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------
! Parameters identifying alternative soil biogeochemistry models.
! (The "bgc" is for biogeochemistry!)
! These should be >0 and unique.
INTEGER, PARAMETER ::                                                         &
  soil_model_1pool = 1,                                                       &
    ! A 1-pool model of soil carbon turnover in which the pool is not
    ! prognostic (not updated). Historically this was only used to
    ! calculate soil respiration when the TRIFFID vegetation model was not
    ! used.
  soil_model_rothc = 2,                                                       &
    ! The RothC (4 pool) model of soil carbon. Historically this was
    ! bundled with the TRIFFID vegetation model, with l_triffid=.TRUE.
    ! effectively implying RothC was used.
  soil_model_ecosse = 3
    ! ECOSSE model of soil carbon and nitrogen.

! Parameters identifying alternative wetland methane substrate models.
! These should be >0 and unique.
INTEGER, PARAMETER ::                                                         &
  ch4_substrate_soil = 1,                                                     &
    ! Soil carbon provides the substrate for wetland methane emissions.
  ch4_substrate_npp = 2,                                                      &
    ! NPP provides the substrate for wetland methane emissions.
  ch4_substrate_soil_resp = 3
    ! Soil respiration provides the substrate for wetland methane emissions.

!-----------------------------------------------------------------------------
! Module variables
!-----------------------------------------------------------------------------

! Items set in namelist jules_soil_biogeochem.

!-----------------------------------------------------------------------------
! Switches that control what soil model is used.
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  soil_bgc_model = soil_model_1pool
    ! Indicates choice of soil model.
    ! Valid values are given by the soil_model_* parameters.

!-----------------------------------------------------------------------------
! Namelist variables used by both 1pool and RothC models.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                     &
  q10_soil = 2.0
    ! Q10 factor for soil respiration.

LOGICAL ::                                                                    &
  l_layeredC = .FALSE.,                                                       &
    ! Switch to select layered soil carbon model.
    ! Note this is not used with the ECOSSE model.
    ! .TRUE.  = use layered model
    ! .FALSE. = no layers (bulk pool)
  l_q10 = .TRUE.,                                                             &
    ! Switch for temperature function for soil respiration.
    ! .TRUE.  = use Q10 formulation
    ! .FALSE. = use RothC formulation
! Switch for bug fix.
    l_soil_resp_lev2 = .FALSE.,                                               &
      ! Switch used to control the soil tempoerature and moisture used
      ! un the soil respiration calculation.
      ! .TRUE.  means use total (frozen+unfrozen) soil moisture.
      ! .FALSE. means use unfrozen soil moisture.
      ! Depending on l_layeredC, l_soil_resp_lev2 can affect the layer from
      ! which the temperature and moisture are taken for respiration.
      ! If l_layeredC=.TRUE.: use T and moisture for each layer.
      ! If l_layeredC=.FALSE.:
      !   l_soil_resp_lev2=T means uses T and moisture from layer 2
      !   l_soil_resp_lev2=F means uses T and moisture from layer 1
    l_ch4_interactive = .FALSE.,                                              &
        ! Switch to couple methane release into the carbon cycle
        ! CH4 flux will be removed from soil carbon pools.
        ! Must have l_ch4_tlayered = .TRUE.
    l_ch4_tlayered = .FALSE.,                                                 &
        ! Switch to calculate CH4 according to layered soil temperature
        ! instead of top 1m average.
    l_ch4_microbe = .FALSE.
        ! Switch to use microbial methane production scheme

INTEGER ::                                                                    &
  ch4_substrate = ch4_substrate_soil,                                         &
    ! Indicates choice of methane substratel model.
    ! Valid values are given by the ch4_substrate_* parameters.
  dim_ch4layer = 1
    ! If methane is calculated from individual soil layer temperatures, this
    ! will be equal to sm_levels (depends on l_ch4_tlayered)

!-----------------------------------------------------------------------------
! Namelist variables used only by the 1-pool model.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                     &
  kaps = 0.5e-8
    ! Specific soil respiration rate at 25 degC and optimum soil moisture
    ! (s-1). Only used for 1-pool model.

!-----------------------------------------------------------------------------
! Namelist variables used only by the RothC model.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                     &
  bio_hum_CN = 10.0,                                                          &
    ! Soil Bio and Hum CN ratio parameter
  sorp = 10.0,                                                                &
    ! Soil inorganic N factor in leaching.
  N_inorg_turnover = 1.0,                                                     &
    ! Inorganic N turnover rate (per 360 days).
  diff_n_pft = 100.0,                                                         &
    ! Inorganic N diffusion in soil (determines how quickly it reaches the
    ! roots after the roots uptake from the soil around them)
    ! per 360 days. Should be quicker than the turnover rate of inorganic
    ! N hence choice of value (100 vs 1).
  tau_resp = 2.0
    ! Parameter controlling decay of respiration with depth (m-1)

REAL(KIND=real_jlslsm) ::                                                     &
  kaps_roth(4) = (/ 3.22e-7, 9.65e-9, 2.12e-8, 6.43e-10 /)
    ! Specific soil respiration rate for RothC (s-1).

!-----------------------------------------------------------------------------
! Namelist variables that are potentially used by more than one soil model.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                     &
  tau_lit = 5.0
    ! Parameter controlling the decay of litter inputs with depth (m-1).

!-----------------------------------------------------------------------------
! Namelist variables used in the CH4 Emission Scheme
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                     &
  t0_ch4 = 273.15,                                                            &
      ! Reference temperature for the Q10 function in the CH4 calculations
      ! (T_0 in equations 74 and 75 in Clark et al. 2010)
  const_ch4_cs = 7.41e-12,                                                    &
      ! Scale factor for CH4 emissions when soil carbon is the substrate
      ! (k in equation 74 in Clark et al. 2010)
      ! In the case of UM simulations const_ch4_cs should be set to specific
      ! values depending on vegation model. These are currently set using the
      ! rose app-upgrade functionality
      ! The UM values are:
      !       IF (l_triffid==false)  const_ch4_cs = 5.41e-12
      !       IF (l_triffid==true)   const_ch4_cs = 5.41e-10
  const_ch4_npp  = 9.99e-3,                                                   &
      ! Scale factor for CH4 emissions when NPP is the substrate
      ! (k in equation 74 in Clark et al. 2010)
  const_ch4_resps  = 4.36e-3,                                                 &
      ! Scale factor for CH4 emissions when soil resp. is the substrate
      ! (k in equation 74 in Clark et al. 2010)
  q10_ch4_cs = 3.7,                                                           &
      ! Q10 factor for CH4 emissions when soil carbon is the substrate
      ! ( Q10_CH4(T_0) in equation 75 in Clark et al. 2010 )
  q10_ch4_npp  = 1.5,                                                         &
      ! Q10 factor for CH4 emissions when NPP is the substrate
      ! ( Q10_CH4(T_0) in equation 75 in Clark et al. 2010 )
  q10_ch4_resps  = 1.5,                                                       &
      ! Q10 factor for CH4 emissions when soil resp. is the substrate
      ! ( Q10_CH4(T_0) in equation 75 in Clark et al. 2010 )
  k2_ch4      = 0.01,                                                         &
      ! Baseline methanogenic respiration rate (hr-1)
  kd_ch4      = 0.0003,                                                       &
      ! Baseline methanogenic death/turnover rate (hr-1)
  rho_ch4     = 47.0,                                                         &
      ! Factor in substrate limitation function (related to half saturation of
      ! substrate for methanogenic respiration) ( (mgC/m3)-1 )
  q10_mic_ch4 = 4.3,                                                          &
      ! Q10 factor for methanogens
  cue_ch4     = 0.03,                                                         &
      ! Carbon use efficiency of methanogenic growth
  mu_ch4      = 0.00042,                                                      &
      ! Threshold growth rate below which methanogens die (hr-1)
  alpha_ch4 = 0.001,                                                          &
      ! Ratio between maintenance and growth respiration rates for methanogens
  frz_ch4 = 0.5,                                                              &
      ! Factor to reduce CH4 substrate production when soil is sufficiently
      ! frozen (only in microbial scheme)
  tau_ch4 = 6.5,                                                              &
      ! Parameter controlling decay of methane oxidation with depth (m-1)
  ch4_cpow = 1.0,                                                             &
      ! Methane (l_ch4_microbe=false) or substrate (l_ch4_microbe=true)
      ! production dependence on soil carbon goes like cs**ch4_cpow
  ev_ch4 = 5.0,                                                               &
      ! Timescale over which methanogenic traits adapt to temperature change
      ! (yr)
  q10_ev_ch4 = 2.2
      ! Q10 for temperature response of methanogenic traits under adaptation

!-----------------------------------------------------------------------------
! Namelist definition
!-----------------------------------------------------------------------------
NAMELIST  / jules_soil_biogeochem/                                            &
! Shared
    soil_bgc_model, ch4_substrate, kaps, kaps_roth, q10_soil, sorp,           &
    n_inorg_turnover, diff_n_pft, tau_resp, tau_lit, bio_hum_CN, l_layeredC,  &
    l_q10, l_soil_resp_lev2, l_ch4_interactive, l_ch4_tlayered, l_ch4_microbe,&
    t0_ch4, const_ch4_cs, const_ch4_npp, const_ch4_resps, q10_ch4_cs,         &
    q10_ch4_npp, q10_ch4_resps, tau_ch4, ch4_cpow, k2_ch4, kd_ch4, rho_ch4,   &
    q10_mic_ch4, cue_ch4, mu_ch4, alpha_ch4, frz_ch4, ev_ch4, q10_ev_ch4

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                       &
  ModuleName = 'JULES_SOIL_BIOGEOCHEM_MOD'

CONTAINS

!#############################################################################

SUBROUTINE check_jules_soil_biogeochem()

USE jules_surface_mod, ONLY:                                                  &
  ! imported scalars
  l_aggregate

USE jules_vegetation_mod, ONLY:                                               &
  ! imported scalars
  l_triffid, l_trif_fire


USE ereport_mod, ONLY: ereport

!-----------------------------------------------------------------------------
! Description:
!   Checks JULES_SOIL_BIOGEOCHEM namelist for consistency.
!-----------------------------------------------------------------------------

IMPLICIT NONE

! Local scalar parameters.
INTEGER :: errorstatus

CHARACTER(LEN=*), PARAMETER ::                                                &
   RoutineName = 'check_jules_soil_biogeochem'   ! Name of this procedure.

! Check that a valid soil model is selected.
SELECT CASE ( soil_bgc_model )
CASE ( soil_model_1pool, soil_model_rothc, soil_model_ecosse )
  !  Acceptable values.
CASE DEFAULT
  errorstatus = 100
  CALL ereport( TRIM(routineName), errorstatus,                               &
                "Invalid value for soil model" )
END SELECT

! Check that a suitable soil model is used with TRIFFID.
IF ( l_triffid ) THEN
  SELECT CASE ( soil_bgc_model )
  CASE ( soil_model_ecosse, soil_model_rothc )
    ! These are OK.
  CASE ( soil_model_1pool )
    errorstatus = 100
    CALL ereport(TRIM(RoutineName), errorstatus,                              &
                 'TRIFFID needs a prognostic soil model - use RothC.')
  END SELECT
END IF

SELECT CASE ( ch4_substrate )
CASE (  ch4_substrate_npp, ch4_substrate_soil, ch4_substrate_soil_resp )
  ! Acceptable values, nothing to do.
CASE DEFAULT
  errorstatus = 100
  CALL ereport(TRIM(RoutineName), errorstatus,                                &
               "Invalid value for ch4_substrate" )
END SELECT

IF ( l_ch4_interactive .AND. .NOT. l_ch4_tlayered ) THEN
  errorstatus = 100
  CALL ereport( TRIM(RoutineName), errorstatus,                               &
      'To couple CH4 to soil carbon (l_ch4_interactive) you must use' //      &
      'the layered soil temperature calculation (l_ch4_tlayered)' )
END IF

! Check that certain soil models are only used with a vegetation model.
SELECT CASE ( soil_bgc_model )
CASE ( soil_model_ecosse, soil_model_rothc )
  IF ( .NOT. l_triffid ) THEN
    errorstatus = 100
    CALL ereport( RoutineName, errorstatus,                                   &
                  'RothC and ECOSSE soil models need a veg model ' //         &
                  '(TRIFFID). Set l_triffid=T.' )
  END IF
END SELECT

! Check that certain soil models are not used with the aggregate surface
! scheme. At present this follows from needing TRIFFID, but test again.
SELECT CASE ( soil_bgc_model )
CASE ( soil_model_ecosse, soil_model_rothc )
  IF ( l_aggregate ) THEN
    errorstatus = 100
    CALL ereport( RoutineName, errorstatus,                                   &
                  'RothC and ECOSSE soil models cannot be used with ' //      &
                  'the aggregated surface scheme (l_aggregate = true)' )
  END IF
END SELECT

! If using the single-pool C model make sure l_q10=T.
! This is done anyway in MICROBE, so might as well do here so that it is
! reported to the user.
IF ( soil_bgc_model == soil_model_1pool ) l_q10 = .TRUE.

! Check that l_layeredC=T is only used with 1-pool and RothC models.
! In particular, layers are set up differently with ECOSSE.
IF ( l_layeredc ) THEN
  SELECT CASE ( soil_bgc_model )
  CASE ( soil_model_1pool, soil_model_rothc )
    ! Fine - nothing more to do.
  CASE DEFAULT
    errorstatus = 100
    CALL ereport(TRIM(RoutineName), errorstatus,                              &
                 'l_layeredc should be FALSE with this soil model.')
  END SELECT
END IF

! Check that ECOSSE is not used with l_trif_fire=T. This combination
! would require that further code is added to ECOSSE.
IF ( soil_bgc_model == soil_model_ecosse .AND. l_trif_fire ) THEN
  errorstatus = 100
  CALL ereport( RoutineName, errorstatus,                                     &
                'ECOSSE cannot be used with l_trif_fire=T' )
END IF

#if defined(UM_JULES)
! UM-only code.

! Until layered soil C is allowed with the UM, ensure it is not selected!
! If layering is requested, reset switch and issue a warning.
IF ( l_layeredC ) THEN
  l_layeredC = .FALSE.
  errorstatus=-100
  CALL ereport(TRIM(RoutineName), errorstatus,                                &
               'l_layeredC reset to FALSE until layering is ' //              &
               'allowed in the UM.')
END IF

! Currently ECOSSE is not allowed in the UM.
IF ( soil_bgc_model == soil_model_ecosse ) THEN
  errorstatus = 100
  CALL ereport( RoutineName, errorstatus,                                     &
                "ECOSSE soil model not allowed with UM." )
END IF
#endif

END SUBROUTINE check_jules_soil_biogeochem

!#############################################################################

SUBROUTINE print_nlist_jules_soil_biogeochem()

USE jules_print_mgr, ONLY: jules_print

IMPLICIT NONE

CHARACTER(LEN=50000) :: lineBuffer

CALL jules_print('jules_soil_biogeochem_mod',                                 &
                 'Contents of namelist jules_soil_biogeochem')

WRITE(lineBuffer,*) ' soil_bgc_model = ', soil_bgc_model
CALL jules_print('jules_soil_biogeochem_mod',lineBuffer)

WRITE(lineBuffer, *) ' l_layeredC = ', l_layeredC
CALL jules_print('jules_soil_biogeochem_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_q10 = ',l_q10
CALL jules_print('jules_soil_biogeochem_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_soil_resp_lev2 = ',l_soil_resp_lev2
CALL jules_print('jules_soil_biogeochem_mod',lineBuffer)

WRITE(lineBuffer, *) ' q10_soil = ', q10_soil
CALL jules_print('jules_soil_biogeochem_mod', lineBuffer)

WRITE(lineBuffer, *) ' kaps = ', kaps
CALL jules_print('jules_soil_biogeochem_mod', lineBuffer)

WRITE(lineBuffer, *) ' kaps_roth = ', kaps_roth
CALL jules_print('jules_soil_biogeochem_mod', lineBuffer)

WRITE(lineBuffer, *) ' sorp = ', sorp
CALL jules_print('jules_soil_biogeochem_mod', lineBuffer)

WRITE(lineBuffer, *) ' bio_hum_CN = ', bio_hum_CN
CALL jules_print('jules_soil_biogeochem_mod', lineBuffer)

WRITE(lineBuffer, *) ' n_inorg_turnover = ', n_inorg_turnover
CALL jules_print('jules_soil_biogeochem_mod', lineBuffer)

WRITE(lineBuffer, *) ' tau_resp = ', tau_resp
CALL jules_print('jules_soil_biogeochem_mod', lineBuffer)

WRITE(lineBuffer, *) ' tau_lit = ', tau_lit
CALL jules_print('jules_soil_biogeochem_mod', lineBuffer)

WRITE(lineBuffer, *) ' diff_n_pft = ', diff_n_pft
CALL jules_print('jules_soil_biogeochem_mod', lineBuffer)

WRITE(lineBuffer, *) ' l_ch4_interactive = ', l_ch4_interactive
CALL jules_print('jules_soil_biogeochem_mod',lineBuffer)

WRITE(lineBuffer, *) ' l_ch4_tlayered = ', l_ch4_tlayered
CALL jules_print('jules_soil_biogeochem_mod',lineBuffer)

WRITE(lineBuffer, *) ' tau_ch4 = ', tau_ch4
CALL jules_print('jules_soil_biogeochem_mod',lineBuffer)

WRITE(lineBuffer, *) ' ch4_cpow = ', ch4_cpow
CALL jules_print('jules_soil_biogeochem_mod',lineBuffer)

WRITE(lineBuffer, *) ' l_ch4_microbe = ', l_ch4_microbe
CALL jules_print('jules_soil_biogeochem_mod',lineBuffer)

WRITE(lineBuffer, *) ' k2_ch4 = ', k2_ch4
CALL jules_print('jules_soil_biogeochem_mod',lineBuffer)

WRITE(lineBuffer, *) ' kd_ch4 = ', kd_ch4
CALL jules_print('jules_soil_biogeochem_mod',lineBuffer)

WRITE(lineBuffer, *) ' rho_ch4 = ', rho_ch4
CALL jules_print('jules_soil_biogeochem_mod',lineBuffer)

WRITE(lineBuffer, *) ' q10_mic_ch4 = ', q10_mic_ch4
CALL jules_print('jules_soil_biogeochem_mod',lineBuffer)

WRITE(lineBuffer, *) ' cue_ch4 = ', cue_ch4
CALL jules_print('jules_soil_biogeochem_mod',lineBuffer)

WRITE(lineBuffer, *) ' mu_ch4 = ', mu_ch4
CALL jules_print('jules_soil_biogeochem_mod',lineBuffer)

WRITE(lineBuffer, *) ' frz_ch4 = ', frz_ch4
CALL jules_print('jules_soil_biogeochem_mod',lineBuffer)

WRITE(lineBuffer, *) ' alpha_ch4 = ', alpha_ch4
CALL jules_print('jules_soil_biogeochem_mod',lineBuffer)

WRITE(lineBuffer, *) ' ev_ch4 = ', ev_ch4
CALL jules_print('jules_soil_biogeochem_mod',lineBuffer)

WRITE(lineBuffer, *) ' q10_ev_ch4 = ', q10_ev_ch4
CALL jules_print('jules_soil_biogeochem_mod',lineBuffer)

WRITE(lineBuffer, *) ' ch4_substrate = ', ch4_substrate
CALL jules_print('jules_soil_biogeochem_mod',lineBuffer)

WRITE(lineBuffer, *) '  t0_ch4 = ', t0_ch4
CALL jules_print('jules_soil_biogeochem_mod', lineBuffer)

WRITE(lineBuffer, *) '  const_ch4_cs = ', const_ch4_cs
CALL jules_print('jules_soil_biogeochem_mod', lineBuffer)

WRITE(lineBuffer, *) '  const_ch4_npp = ', const_ch4_npp
CALL jules_print('jules_soil_biogeochem_mod', lineBuffer)

WRITE(lineBuffer, *) '  const_ch4_resps = ', const_ch4_resps
CALL jules_print('jules_soil_biogeochem_mod', lineBuffer)

WRITE(lineBuffer, *) '  q10_ch4_cs = ', q10_ch4_cs
CALL jules_print('jules_soil_biogeochem_mod', lineBuffer)

WRITE(lineBuffer, *) '  q10_ch4_npp = ', q10_ch4_npp
CALL jules_print('jules_soil_biogeochem_mod', lineBuffer)

WRITE(lineBuffer, *) 'q10_ch4_resps = ', q10_ch4_resps
CALL jules_print('jules_soil_biogeochem_mod', lineBuffer)

CALL jules_print('jules_soil_biogeochem_mod',                                 &
    '- - - - - - end of namelist - - - - - -')

END SUBROUTINE print_nlist_jules_soil_biogeochem

!#############################################################################

#if defined(UM_JULES) && !defined(LFRIC)

SUBROUTINE read_nml_jules_soil_biogeochem (unitnumber)

! Description:
!  Read the JULES_SOIL_BIOGEOCHEM namelist

USE setup_namelist,   ONLY: setup_nml_type
USE check_iostat_mod, ONLY: check_iostat
USE UM_parcore,       ONLY: mype


USE parkind1,         ONLY: jprb, jpim
USE yomhook,          ONLY: lhook, dr_hook

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: unitnumber
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName=                                   &
                                'READ_NML_JULES_SOIL_BIOGEOCHEM'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=errormessagelength) :: iomessage

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 2
INTEGER, PARAMETER :: n_real = 27 + 4
INTEGER, PARAMETER :: n_log = 6

TYPE my_namelist
  SEQUENCE
  INTEGER :: soil_bgc_model
  INTEGER :: ch4_substrate
  REAL(KIND=real_jlslsm) :: q10_soil
  REAL(KIND=real_jlslsm) :: kaps
  REAL(KIND=real_jlslsm) :: kaps_roth(4)
  REAL(KIND=real_jlslsm) :: sorp
  REAL(KIND=real_jlslsm) :: bio_hum_cn
  REAL(KIND=real_jlslsm) :: n_inorg_turnover
  REAL(KIND=real_jlslsm) :: tau_resp
  REAL(KIND=real_jlslsm) :: tau_lit
  REAL(KIND=real_jlslsm) :: tau_ch4
  REAL(KIND=real_jlslsm) :: ch4_cpow
  REAL(KIND=real_jlslsm) :: diff_n_pft
  REAL(KIND=real_jlslsm) :: t0_ch4
  REAL(KIND=real_jlslsm) :: const_ch4_cs
  REAL(KIND=real_jlslsm) :: const_ch4_npp
  REAL(KIND=real_jlslsm) :: const_ch4_resps
  REAL(KIND=real_jlslsm) :: q10_ch4_cs
  REAL(KIND=real_jlslsm) :: q10_ch4_npp
  REAL(KIND=real_jlslsm) :: q10_ch4_resps
  REAL(KIND=real_jlslsm) :: k2_ch4
  REAL(KIND=real_jlslsm) :: kd_ch4
  REAL(KIND=real_jlslsm) :: rho_ch4
  REAL(KIND=real_jlslsm) :: q10_mic_ch4
  REAL(KIND=real_jlslsm) :: cue_ch4
  REAL(KIND=real_jlslsm) :: mu_ch4
  REAL(KIND=real_jlslsm) :: frz_ch4
  REAL(KIND=real_jlslsm) :: alpha_ch4
  REAL(KIND=real_jlslsm) :: ev_ch4
  REAL(KIND=real_jlslsm) :: q10_ev_ch4
  LOGICAL :: l_layeredC
  LOGICAL :: l_q10
  LOGICAL :: l_soil_resp_lev2
  LOGICAL :: l_ch4_interactive
  LOGICAL :: l_ch4_tlayered
  LOGICAL :: l_ch4_microbe
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                         &
                        zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in = n_int,              &
                    n_real_in = n_real, n_log_in = n_log)

IF (mype == 0) THEN

  READ (UNIT = unitnumber, NML = jules_soil_biogeochem,                       &
        IOSTAT = errorstatus, IOMSG = iomessage)
  CALL check_iostat(errorstatus, "namelist jules_soil_biogeochem",            &
                    iomessage)

  my_nml % soil_bgc_model    = soil_bgc_model
  my_nml % ch4_substrate     = ch4_substrate
  my_nml % q10_soil          = q10_soil
  my_nml % kaps              = kaps
  my_nml % kaps_roth         = kaps_roth
  my_nml % sorp              = sorp
  my_nml % bio_hum_cn        = bio_hum_cn
  my_nml % n_inorg_turnover  = n_inorg_turnover
  my_nml % tau_resp          = tau_resp
  my_nml % tau_lit           = tau_lit
  my_nml % tau_ch4           = tau_ch4
  my_nml % ch4_cpow          = ch4_cpow
  my_nml % diff_n_pft        = diff_n_pft
  my_nml % l_layeredC        = l_layeredC
  my_nml % l_q10             = l_q10
  my_nml % l_soil_resp_lev2  = l_soil_resp_lev2
  my_nml % l_ch4_interactive = l_ch4_interactive
  my_nml % l_ch4_tlayered    = l_ch4_tlayered
  my_nml % l_ch4_microbe     = l_ch4_microbe
  my_nml % t0_ch4           = t0_ch4
  my_nml % const_ch4_cs     = const_ch4_cs
  my_nml % const_ch4_npp    = const_ch4_npp
  my_nml % const_ch4_resps  = const_ch4_resps
  my_nml % q10_ch4_cs       = q10_ch4_cs
  my_nml % q10_ch4_npp      = q10_ch4_npp
  my_nml % q10_ch4_resps    = q10_ch4_resps
  my_nml % k2_ch4           = k2_ch4
  my_nml % kd_ch4           = kd_ch4
  my_nml % rho_ch4          = rho_ch4
  my_nml % q10_mic_ch4      = q10_mic_ch4
  my_nml % cue_ch4          = cue_ch4
  my_nml % mu_ch4           = mu_ch4
  my_nml % frz_ch4          = frz_ch4
  my_nml % alpha_ch4        = alpha_ch4
  my_nml % ev_ch4           = ev_ch4
  my_nml % q10_ev_ch4       = q10_ev_ch4

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN
  soil_bgc_model    = my_nml % soil_bgc_model
  ch4_substrate     = my_nml % ch4_substrate
  q10_soil          = my_nml % q10_soil
  kaps              = my_nml % kaps
  kaps_roth         = my_nml % kaps_roth
  sorp              = my_nml % sorp
  bio_hum_CN        = my_nml % bio_hum_CN
  n_inorg_turnover  = my_nml % n_inorg_turnover
  tau_resp          = my_nml % tau_resp
  tau_lit           = my_nml % tau_lit
  tau_ch4           = my_nml % tau_ch4
  ch4_cpow          = my_nml % ch4_cpow
  diff_n_pft        = my_nml % diff_n_pft
  l_layeredC        = my_nml % l_layeredC
  l_q10             = my_nml % l_q10
  l_soil_resp_lev2  = my_nml % l_soil_resp_lev2
  l_ch4_interactive = my_nml % l_ch4_interactive
  l_ch4_tlayered    = my_nml % l_ch4_tlayered
  l_ch4_microbe     = my_nml % l_ch4_microbe
  t0_ch4           = my_nml % t0_ch4
  const_ch4_cs     = my_nml % const_ch4_cs
  const_ch4_npp    = my_nml % const_ch4_npp
  const_ch4_resps  = my_nml % const_ch4_resps
  q10_ch4_cs       = my_nml % q10_ch4_cs
  q10_ch4_npp      = my_nml % q10_ch4_npp
  q10_ch4_resps    = my_nml % q10_ch4_resps
  k2_ch4           = my_nml % k2_ch4
  kd_ch4           = my_nml % kd_ch4
  rho_ch4          = my_nml % rho_ch4
  q10_mic_ch4      = my_nml % q10_mic_ch4
  cue_ch4          = my_nml % cue_ch4
  mu_ch4           = my_nml % mu_ch4
  frz_ch4          = my_nml % frz_ch4
  alpha_ch4        = my_nml % alpha_ch4
  ev_ch4           = my_nml % ev_ch4
  q10_ev_ch4       = my_nml % q10_ev_ch4
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                         &
                        zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_nml_jules_soil_biogeochem

#endif

END MODULE jules_soil_biogeochem_mod
