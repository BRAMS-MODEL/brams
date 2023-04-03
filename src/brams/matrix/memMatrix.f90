MODULE memMatrix
   use ModNamelistFile, only: namelistFile
   USE aer1_list      , ONLY: nmodes_aer=>nmodes,matrix_level,N_matrix_level
   
 
   integer :: aerosol
   real :: aerfrq
   integer,PARAMETER :: mech=N_matrix_level

   type matrixType
         !Inputs to matrix
         DOUBLE PRECISION,dimension(:,:)  ,pointer  :: emis_mas      !< mass emission rates [ug/m^3/s]
         DOUBLE PRECISION,dimension(:)    ,pointer  :: tk            !< absolute temperature [k]
         DOUBLE PRECISION,dimension(:)    ,pointer  :: rh            !< relative humidity [0-1]
         DOUBLE PRECISION,dimension(:)    ,pointer  :: pres          !< ambient pressure [pa]
         DOUBLE PRECISION,dimension(:)    ,pointer  :: aqso4rate     !< in-cloud so4 production rate [ug/m^3/s]
         DOUBLE PRECISION,dimension(:)    ,pointer  :: wupdraft      !< cloud updraft velocity [m/s]
         DOUBLE PRECision,DIMENSION(:,:)  ,pointer  :: dens_mode_dry !< Average mode density calculated from component concentrations [g/cm^3]
         DOUBLE PRECision,DIMENSION(:)    ,pointer  :: dn0           !
         DOUBLE PRECision,DIMENSION(:)    ,pointer  :: pi0           !
         DOUBLE PRECision,DIMENSION(:)    ,pointer  :: pp            !
         DOUBLE PRECision,DIMENSION(:)    ,pointer  :: ustar         !
         DOUBLE PRECision,DIMENSION(:)    ,pointer  :: tstar         !
         DOUBLE PRECision                           :: zi            !
         !Inputs/Outputs to/from matrix
         DOUBLE PRECISION,dimension(:,:)  ,pointer  :: aero          !< aerosol conc. [ug/m^3] or [#/m^3]
         DOUBLE PRECISION,dimension(:,:)  ,pointer  :: gas           !< gas-phase conc. [ug/m^3]
         !Outputs from matrix
         DOUBLE PRECision,DIMENSION(:,:)  ,pointer  :: nAct          !< Activated number concentration for each mode [#/m^3]
         DOUBLE PRECision,DIMENSION(:,:)  ,pointer  :: vddep_aero    !< Dry deposition velocity [m/s]
         DOUBLE PRECision,DIMENSION(:,:)  ,pointer  :: diam          !< Dry diameter of average mass for each mode [m]
         DOUBLE PRECISION,dimension(:,:,:),pointer  :: diag          !< budget or tendency diagnostics [ug/m^3/s] or [#/m^3/s]
         DOUBLE PRECision,DIMENSION(:,:)  ,pointer  :: ac            !< minimum dry radius to activate for each mode [um].
         DOUBLE PRECision,DIMENSION(:,:)  ,pointer  :: aer2mp_eff     !< save aerosol properties for uphysics
	 
   END type matrixType
   TYPE(matrixType),dimension(:),allocatable :: matrixVar !< Type used to call Matrix

   !- for 2 or 3 time levels integration schemes (RK 2 and 3)
   type matrixType_rk
         DOUBLE PRECISION,dimension(:,:)  ,pointer  :: aero  !< aerosol conc. [ug/m^3] or [#/m^3]
         DOUBLE PRECISION,dimension(:,:)  ,pointer  :: gas   !< gas-phase conc. [ug/m^3]
   END type matrixType_rk
   TYPE(matrixType_rk),dimension(:),allocatable :: matrixVarRK!< Type used to store intermediate
    

   INTEGER :: NumberOfColunms		     !< Total of columns in processor
   INTEGER, ALLOCATABLE :: iPos(:)	     !< i position of each column position
   INTEGER, ALLOCATABLE :: jPos(:)	     !< j position of each column position
   INTEGER, ALLOCATABLE :: nColumn(:,:)      !< Position of each (i,j) point in column strip
   integer,PARAMETER :: naerovars=51
   integer,PARAMETER :: nextra=3
   !-------------------------------------------------------------------------------------------------------------------------
   !     2. Set the number of quadrature points per mode (1-2); must use NPOINTS=1 for the present.
   !-------------------------------------------------------------------------------------------------------------------------
      integer, parameter :: nPoints=1
   !srf integer :: nmodes  ! define in aer1_list
   integer,PARAMETER :: nAeroBox=naerovars+nextra
   integer,PARAMETER :: nWeights=16*npoints
   INTEGER, PARAMETER :: nmodes=16 !internal for Matrix

     !from TRAMP_param
!----------------------------------------------------------------------------------------------------------------------
!
!  GEOMETRIC AND SCIENTIFIC CONSTANTS; DERIVED CONVERSION FACTORS.
!
!----------------------------------------------------------------------------------------------------------------------
   real(8), parameter :: pi          = 3.141592653589793d+00
   real(8), parameter :: pi6         = pi/6.0d+00
   real(8), parameter :: avo         = 6.0221367d+23            ! avogadro's number     [#/mole]
   real(8), parameter :: rgas_si     = 8.3145d+00               ! univers. gas constant [j/mol/k]
   real(8), parameter :: mw_h2so4    = 98.07948d+00             ! molar mass of h2so4   [g/mole]
   real(8), parameter :: mw_so4      = 96.06360d+00             ! molar mass of so4=    [g/mole]
   real(8), parameter :: mw_nh3      = 17.03056d+00             ! molar mass of nh3     [g/mole]
   real(8), parameter :: mw_h2o      = 18.01528d+00             ! molar mass of h2o     [g/mole]
   real(8), parameter :: mw_nh42so4  = 132.1406d+00             ! molar mass of nh42so4 [g/mole]
   real(8), parameter :: ugm3_ncm3   = 1.0d-12 * avo / mw_so4   ! [ugso4/m^3] to [#/cm^3]
   real(8), parameter :: convnh3     = rgas_si / mw_nh3         ! used in [ug nh3/m^3] to [ppmv]
   real(8), parameter :: rho_nh42so4 = 1.77d+00                 ! density of dry (nh4)2so4 [g/cm^3]
   real(8), parameter :: rho_h2so4   = 1.84d+00                 ! density of pure h2so4    [g/cm^3] - crc
   real(8), parameter :: rho_h2o     = 1.00d+00                 ! density of pure h2so4    [g/cm^3] - crc
   real(8), parameter :: densp       = 1.40d+00                 ! default ambient particle density [g/cm^3]
   real(8), parameter :: mw_no3      = 62.0049d+00              ! molar mass of no3     [g/mole]
   real(8), parameter :: mw_hno3     = 63.01d+00                ! molar mass of hno3    [g/mole]
   real(8), parameter :: mw_nh4      = 18.03851d+00             ! molar mass of nh4     [g/mole]
   !-------------------------------------------------------------------------------------------------------------------

   ! CONV_DP_TO_MASS converts Dp^3 [m^3] to particle mass [ug].
   ! CONV_MASS_TO_DP converts particle mass [ug] to Dp^3 [m^3].
   ! CONV_VOL_TO_DP_FAC converts particle volume [ug] to Dp^3 [m^3].
   ! The factor 1.0D+12 converts [m^3] to [cm^3] and [g] to [ug].
   !-------------------------------------------------------------------------------------------------------------------
   real(8), parameter :: conv_dp_to_mass    = 1.0d+12 * pi6 * densp
   real(8), parameter :: conv_mass_to_dp    = 1.0d+00 / conv_dp_to_mass
   real(8), parameter :: conv_vol_to_dp_fac = 1.0d-12 / pi6
   !-------------------------------------------------------------------------------------------------------------------
   ! Miniumum and maximum values of average mode diameters. [m]
   ! These are needed when number and/or mass concentrations are very small.
   !-------------------------------------------------------------------------------------------------------------------
   real(8), parameter :: dpmin_global =  0.001d-06   ! [m] -  1 nm
   real(8), parameter :: dpmax_global = 20.000d-06   ! [m] - 20 um
!----------------------------------------------------------------------------------------------------------------------
!
!  MODEL PARAMETERS AND VARIABLES THAT MAY NEED TO BE SET BY THE USER.
!
!----------------------------------------------------------------------------------------------------------------------
   integer, parameter :: aunit1            = 90 ! logical unit # - log file of module
   integer, parameter :: aunit2            = 91 ! logical unit # - test of coag. coef.
   integer, parameter :: nemis_spcs        = 10 ! number of emissions variables
   integer, parameter :: ndiag_aero        = 15 ! number of aerosol diagnostics collected
   integer, parameter :: kij_ndgs_set      = 31 ! default value=81; if no_microphysics=.true., set to 3 to save storage
   integer, parameter :: imtr_method       =  1 ! =1 no cut of pdf, =2 fixed-dp cut, =3 variable-dp cut as in cmaq
   integer, parameter :: activation_scheme =  3 ! =1 uses typical solubility only, =2 detailed multimodal activation
                                                !KML =3 calculate only the hygroscopicity parameter and leave it to the microphysics.
   integer, parameter :: update_kij        =  1 ! =0 use time-independent coagulation coefficients from lookup tables
                                                ! =1 use time-  dependent coagulation coefficients from lookup tables
   logical, parameter :: write_log   = .false.  ! write matrix log to unit aunit1: default setting is .false.
   logical, parameter :: mass_adj    = .true.   ! enforce precise mass conservation: default setting is .true.
   logical, parameter :: cpu_stats   = .false.  ! timer for sections of the matrix code
   logical, parameter :: update_dp   = .true.   ! update particle diameters at each time step: default is .true.
   logical, parameter :: update_diam = .true.   ! update particle diameters at each time step in global diam array
   logical, parameter :: update_vdep = .false.  ! update particle diameters at each time step in global array
   logical, parameter :: set_intermodal_transfer  = .true.  ! do akk -> acc transfer if mode akk is defined
   logical, parameter :: no_microphysics          = .false. ! no-microphysics option
   logical, parameter :: no_microphysics_w_thermo = .true.  ! do gas-particle partitioning when no_microphysics=.true.
   logical, parameter :: do_npf                   = .true.  ! include secondary particle formation
   logical, parameter :: discrete_eval_option     = .false. ! for evaluation with results from the discrete pdf model
   logical, parameter :: activation_comparison    = .false. ! for comparison of aerosol activation w/ all 8 mechanisms
   !-------------------------------------------------------------------------------------------------------------------
   ! AQSO4RATE_MIN is the min, aqueous SO4 production rate for call of the activation routine.
   ! The default value is 4.43D-11 [ug/m^3/s], equivalent to 1000 [molecules/cm3/h] = 1.6D-07 [ugSO4/m^3/h].
   !-------------------------------------------------------------------------------------------------------------------
   real(8), parameter :: aqso4rate_min  = 4.43d-11
   !-------------------------------------------------------------------------------------------------------------------
   ! The Maximum Inorganic Volume Fraction (MIVF) in modes DD1, DD2, BC1, and BC2.
   !-------------------------------------------------------------------------------------------------------------------
   real(8), parameter :: mivf_ddd = 0.05d+00   !
   real(8), parameter :: mivf_bc1 = 0.05d+00   !- these two are from mzj 2002, "analysis ..."
   real(8), parameter :: mivf_bc2 = 0.20d+00   !/
   !-------------------------------------------------------------------------------------------------------------------
   ! Scale factor for the geometric mean diameters of the emissions lognormals.
   !-------------------------------------------------------------------------------------------------------------------
   real(8), parameter :: scale_emis_diam = 1.0d+00
   !-------------------------------------------------------------------------------------------------------------------
   ! The minimum value of a number concentration or mass concentration leaving MATRIX.
   !-------------------------------------------------------------------------------------------------------------------
   real(8), parameter :: minconc   = 1.0d-15   ! [ug/m^3] and [#/m^3]
   real(8), parameter :: tinynumer = 1.0d-30   !
   real(8), parameter :: tinydenom = 1.0d-30   !
   !-------------------------------------------------------------------------------------------------------------------
   ! ZHEIGHT(I) is the mid-level height of model vertical layer I in [km].
   ! It is used in the calculation of ionization rates in the ion-ion
   ! recombination nucleation scheme, and in computing the pre-calculated
   ! factor in the condensational sink.
   !
   ! Set ZHEIGHT to global-average values typical of the vertical structure of the host GCM.
   !-------------------------------------------------------------------------------------------------------------------
   integer, parameter :: nlays = 42          ! number of model vertical layers
   real(8) :: zheight(nlays)                 ! typical (cmaq) mid-layer heights [km]
!LFc M 20 model
!LFc      DATA ZHEIGHT/ 0.007D+00,  0.024D+00,  0.052D+00,  0.100D+00,  0.210D+00,
!LFc     &              0.390D+00,  0.640D+00,  0.950D+00,  1.300D+00,  1.740D+00,
!LFc     &              2.260D+00,  2.810D+00,  3.390D+00,  4.000D+00,  4.700D+00,
!LFc     &              5.400D+00,  6.200D+00,  7.200D+00,  8.400D+00,  10.00D+00,
!LFc     &              12.40D+00 /
!LFc F 40 model
   DATA ZHEIGHT/0.01D+00,0.02D+00,0.04D+00,0.06D+00,0.09D+00,1.2D+00,1.5D+00,2.D+00,2.4D+00,3.D+00 &
               ,3.5D+00,4.D+00,6.7D+00,7.4D+00,8.1D+00,8.5D+00,9.D+00,10.D+00,11.D+00,12.D+00 &
               ,13.D+00,14.D+00,15.D+00,16.D+00,18.D+00,19.D+00,21.D+00,24.D+00,28.D+00,32.D+00 &
               ,36.D+00,40.D+00,44.D+00,48.D+00,53.D+00,58.D+00,62.D+00,66.D+00,72.D+00,80.D+00 &
               ,85.D+00,90.D+00/
!----------------------------------------------------------------------------------------------------------------------
!  efault values for the box model - do not change these.
!----------------------------------------------------------------------------------------------------------------------
!  INTEGER, PARAMETER :: NLAYS = 21          ! number of model vertical layers
!  REAL(8) :: ZHEIGHT(NLAYS)                 ! typical (CMAQ) mid-layer heights [km]
!  DATA ZHEIGHT/ 0.007D+00,  0.024D+00,  0.052D+00,  0.100D+00,  0.210D+00,
!                0.390D+00,  0.640D+00,  0.950D+00,  1.300D+00,  1.740D+00,
!                2.260D+00,  2.810D+00,  3.390D+00,  4.000D+00,  4.700D+00,
!                5.400D+00,  6.200D+00,  7.200D+00,  8.400D+00,  10.00D+00,
!                12.40D+00 /
!----------------------------------------------------------------------------------------------------------------------
   ! Characteristic lognormal parameters for each mode: DG_XXX[um],SG_XXX[1].
   !
   ! The Dg values from Easter et al., 2004 (E04) are the "diagnosed number" values, not the "emitted" values.
   !-------------------------------------------------------------------------------------------------------------------
   real(8), parameter :: dg_akk = 0.026d+00      ! e04, table 2, aitken       mode
   real(8), parameter :: dg_acc = 0.110d+00      ! e04, table 2, accumulation mode
   !-------------------------------------------------------------------------------------------------------------------
   ! DLW: 021507: These dust and sea salt Dg values were calculated approximate emissions sizes used at GISS.
   !
   ! The number mean diameters for the four GISS dust size classes were 0.46, 2.94, 5.88, and 11.77 micrometers.
   ! Size classes 1 and 2 were averaged with a 10:1 ratio, and the average converted to a lognormal geometric
   ! mean diameter for an assumed geometric standard deviation of 1.8.
   ! Likewise, size classes 3 and 4 were averaged with a 10:1 ratio, and the average converted to a lognormal
   ! geometric mean diameter for an assumed geometric standard deviation of 1.8.
   !
   ! The number mean diameters for the two GISS sea salt size classes were 0.44 and 5.0 micrometers, and were
   ! converted to lognormal geometric mean diameters for an assumed geometric standard deviation of 1.8 for the
   ! smaller (accumulation) size class and a standard deviation of 2.0 for the larger (coarse) size class.
    !-------------------------------------------------------------------------------------------------------------------
   real(8), parameter :: dg_dd1 = 0.580d+00 *2.     ! set to match giss dust emissions for average of sizes 1 & 2
   real(8), parameter :: dg_dd2 = 1.000d+00 *2.     ! set to match giss dust emissions for average of sizes 3 & 4
   real(8), parameter :: dg_ds1 = 0.580d+00 *2.     ! set to match giss dust emissions for average of sizes 1 & 2
   real(8), parameter :: dg_ds2 = 1.00d+00 *2.     ! set to match giss dust emissions for average of sizes 3 & 4
   real(8), parameter :: dg_ssa = 0.06d+00 *2.     ! set to match giss sea salt emissions
!   real(8), parameter :: dg_ssa = 0.370d+00 *2.     ! set to match giss sea salt emissions
   real(8), parameter :: dg_ssc = 1.d+00 *2.     ! set to match giss sea salt emissions
   real(8), parameter :: dg_sss = 0.690d+00 *2.     ! 10:1 average of modes ssa and ssc
   !-------------------------------------------------------------------------------------------------------------------
   ! DLW: 021507: End of emissions sizes used at GISS.
   !-------------------------------------------------------------------------------------------------------------------
   real(8), parameter :: dg_occ = 0.050d+00      ! geo. avg. of e04, table 2, aitken and accumulation modes
   real(8), parameter :: dg_bc1 = 0.050d+00      ! geo. avg. of e04, table 2, aitken and accumulation modes
   real(8), parameter :: dg_bc2 = 0.100d+00      ! 1.0164 * geo. avg. of e04, table 2, akk and acc modes, w/  5% shell
   real(8), parameter :: dg_bc3 = 0.100d+00      ! 1.0627 * geo. avg. of e04, table 2, akk and acc modes, w/ 20% shell
   real(8), parameter :: dg_dbc = 0.330d+00      ! geo. avg. of e04, table 2, accumulation and coarse dust modes
   real(8), parameter :: dg_boc = 0.100d+00      ! assuming additive volumes for bc1 and occ
   real(8), parameter :: dg_bcs = 0.070d+00      ! assuming add. vol. for bc1 and akk (acc) and greater weight for akk
   real(8), parameter :: dg_ocs = 0.070d+00      ! assuming add. vol. for bc1 and akk (acc) and greater weight for akk
   real(8), parameter :: dg_mxx = 0.300d+00      ! value is midrange considering all modes
   real(8), parameter :: sg_akk = 1.600d+00      ! e04, table 2, aitken       mode
   real(8), parameter :: sg_acc = 1.800d+00      ! e04, table 2, accumulation mode
   real(8), parameter :: sg_dd1 = 1.800d+00      ! e04, table 2, accumulation mode
   real(8), parameter :: sg_dd2 = 1.800d+00      ! e04, table 2, coarse       mode
   real(8), parameter :: sg_ds1 = 1.800d+00      ! e04, table 2, accumulation mode
   real(8), parameter :: sg_ds2 = 1.800d+00      ! e04, table 2, coarse       mode
   real(8), parameter :: sg_ssa = 1.800d+00      ! e04, table 2, accumulation mode
   real(8), parameter :: sg_ssc = 2.000d+00      ! e04, table 2, coarse       mode
   real(8), parameter :: sg_sss = 2.000d+00      ! same as ssc
   real(8), parameter :: sg_occ = 1.800d+00      ! e04, table 2, accumulation mode
   real(8), parameter :: sg_bc1 = 1.800d+00      ! e04, table 2, accumulation mode
   real(8), parameter :: sg_bc2 = 1.800d+00      ! e04, table 2, accumulation mode
   real(8), parameter :: sg_bc3 = 1.800d+00      ! e04, table 2, accumulation mode
   real(8), parameter :: sg_dbc = 1.800d+00      ! same as parent modes
   real(8), parameter :: sg_boc = 1.800d+00      ! same as parent modes
   real(8), parameter :: sg_bcs = 1.800d+00      ! same as parent modes
   real(8), parameter :: sg_ocs = 1.800d+00      ! same as parent modes
   real(8), parameter :: sg_mxx = 2.000d+00      ! likely a broad mode
   !-------------------------------------------------------------------------------------------------------------------
   ! Lognormal parameters for emissions into each mode: DG_XXX_EMIS[um],SG_XXX_EMIS[1].
   !
   ! These are used to convert mass emission rates to number emission rates.
   ! The Dg values from Easter et al., 2004 (E04) are the "emitted" values, not the "diagnosed number" values.
   ! All modes are assigned a value, even if they do not receive primary particles.
   !-------------------------------------------------------------------------------------------------------------------
   real(8), parameter :: dg_akk_emis = 0.013d+00      ! e04, table 2, aitken       mode
   real(8), parameter :: dg_acc_emis = 0.068d+00      ! e04, table 2, accumulation mode
   !-------------------------------------------------------------------------------------------------------------------
   ! DLW: 021507: These dust and sea salt Dg values were calculated approximate emissions sizes used at GISS.
   !
   ! The number mean diameters for the four GISS dust size classes were 0.46, 2.94, 5.88, and 11.77 micrometers.
   ! Size classes 1 and 2 were averaged with a 10:1 ratio, and the average converted to a lognormal geometric
   ! mean diameter for an assumed geometric standard deviation of 1.8.
   ! Likewise, size classes 3 and 4 were averaged with a 10:1 ratio, and the average converted to a lognormal
   ! geometric mean diameter for an assumed geometric standard deviation of 1.8.
   !
   ! The number mean diameters for the two GISS sea salt size classes were 0.44 and 5.0 micrometers, and were
   ! converted to lognormal geometric mean diameters for an assumed geometric standard deviation of 1.8 for the
   ! smaller (accumulation) size class and a standard deviation of 2.0 for the larger (coarse) size class.
   !-------------------------------------------------------------------------------------------------------------------
   real(8), parameter :: dg_dd1_emis = 0.580d+00 *2.     ! set to match giss dust emissions for average of sizes 1 & 2
   real(8), parameter :: dg_dd2_emis = 1.000d+00 *2.     ! set to match giss dust emissions for average of sizes 3 & 4
   real(8), parameter :: dg_ds1_emis = 0.580d+00 *2.     ! set to match giss dust emissions for average of sizes 1 & 2
   real(8), parameter :: dg_ds2_emis = 1.00d+00 *2.     ! set to match giss dust emissions for average of sizes 3 & 4
   real(8), parameter :: dg_ssa_emis = 0.060d+00 *2     ! set to match giss sea salt emissions
!   real(8), parameter :: dg_ssa_emis = 0.370d+00 *2     ! set to match giss sea salt emissions
   real(8), parameter :: dg_ssc_emis = 1.000d+00 *2.     ! set to match giss sea salt emissions
   real(8), parameter :: dg_sss_emis = 0.690d+00 *2.     ! 10:1 average of modes ssa and ssc
   !-------------------------------------------------------------------------------------------------------------------
   ! DLW: 021507: End of emissions sizes used at GISS.
   !-------------------------------------------------------------------------------------------------------------------
   real(8), parameter :: dg_occ_emis = 0.050d+00      ! geometric average of the akk and acc values in e04
   real(8), parameter :: dg_bc1_emis = 0.050d+00      ! geometric average of the akk and acc values in e04
   real(8), parameter :: dg_bc2_emis = 0.100d+00      !   currently no emissions into this mode
   real(8), parameter :: dg_bc3_emis = 0.100d+00      !   currently no emissions into this mode
   real(8), parameter :: dg_dbc_emis = 0.300d+00      !   currently no emissions into this mode
   real(8), parameter :: dg_boc_emis = 0.100d+00      ! assuming additive volumes for bc1 and occ
   real(8), parameter :: dg_bcs_emis = 0.140d+00      !   currently no emissions into this mode
   real(8), parameter :: dg_ocs_emis = 0.140d+00      !   currently no emissions into this mode
   real(8), parameter :: dg_mxx_emis = 0.500d+00      !   currently no emissions into this mode
   real(8), parameter :: sg_akk_emis = 1.600d+00      ! e04, table 2, aitken       mode
   real(8), parameter :: sg_acc_emis = 1.800d+00      ! e04, table 2, accumulation mode
   real(8), parameter :: sg_dd1_emis = 1.800d+00      ! e04, table 2, accumulation mode
   real(8), parameter :: sg_dd2_emis = 1.800d+00      ! e04, table 2, coarse       mode
   real(8), parameter :: sg_ds1_emis = 1.800d+00      !   currently no emissions into this mode
   real(8), parameter :: sg_ds2_emis = 1.800d+00      !   currently no emissions into this mode
   real(8), parameter :: sg_ssa_emis = 1.800d+00      ! e04, table 2, accumulation mode
   real(8), parameter :: sg_ssc_emis = 2.000d+00      ! e04, table 2, coarse       mode
   real(8), parameter :: sg_sss_emis = 2.000d+00      ! same as ssc
   real(8), parameter :: sg_occ_emis = 1.800d+00      ! e04, table 2, accumulation mode
   real(8), parameter :: sg_bc1_emis = 1.800d+00      ! e04, table 2, accumulation mode
   real(8), parameter :: sg_bc2_emis = 1.800d+00      !   currently no emissions into this mode
   real(8), parameter :: sg_bc3_emis = 1.800d+00      !   currently no emissions into this mode
   real(8), parameter :: sg_dbc_emis = 1.800d+00      !   currently no emissions into this mode
   real(8), parameter :: sg_boc_emis = 1.800d+00      ! same as bc1 and occ modes
   real(8), parameter :: sg_bcs_emis = 1.800d+00      !   currently no emissions into this mode
   real(8), parameter :: sg_ocs_emis = 1.800d+00      !   currently no emissions into this mode
   real(8), parameter :: sg_mxx_emis = 2.000d+00      !   currently no emissions into this mode
   !-------------------------------------------------------------------------------------------------------------------
   ! KAPPAI_XXX is the activating fraction for mode XXX.
   !-------------------------------------------------------------------------------------------------------------------
   real(8), parameter :: kappai_akk = 0.0d+00
   real(8), parameter :: kappai_acc = 1.0d+00
   real(8), parameter :: kappai_dd1 = 0.0d+00
   real(8), parameter :: kappai_dd2 = 0.0d+00
   real(8), parameter :: kappai_ds1 = 1.0d+00
   real(8), parameter :: kappai_ds2 = 1.0d+00
   real(8), parameter :: kappai_ssa = 1.0d+00
   real(8), parameter :: kappai_ssc = 1.0d+00
   real(8), parameter :: kappai_sss = 1.0d+00
   real(8), parameter :: kappai_occ = 0.7d+00
   real(8), parameter :: kappai_bc1 = 0.0d+00
   real(8), parameter :: kappai_bc2 = 1.0d+00
   real(8), parameter :: kappai_bc3 = 1.0d+00
   real(8), parameter :: kappai_dbc = 0.0d+00
   real(8), parameter :: kappai_boc = 0.5d+00
   real(8), parameter :: kappai_bcs = 1.0d+00
   real(8), parameter :: kappai_ocs = 1.0d+00
   real(8), parameter :: kappai_mxx = 1.0d+00

!KML-------------------------------------------------------------------------------------------------------------------
! KAPPAPK_XXX is the hygroscopicity parameter for mode XXX.
!-------------------------------------------------------------------------------------------------------------------
   real(8), parameter :: kappapk_SULF = 0.90d+00
   real(8), parameter :: kappapk_BCAR = 0.00d+00
   real(8), parameter :: kappapk_OCAR = 0.10d+00 
   real(8), parameter :: kappapk_DUST = 0.00d+00
   real(8), parameter :: kappapk_SEAS = 1.28d+00 

!-------------------------------------------------------------------------------------------------------------------
! Solubility per mode
!-------------------------------------------------------------------------------------------------------------------
   real(8), parameter :: solu_akk = 1.0d+00
   real(8), parameter :: solu_acc = 1.0d+00
   real(8), parameter :: solu_dd1 = 0.5d+00
   real(8), parameter :: solu_dd2 = 0.5d+00
   real(8), parameter :: solu_ds1 = 1.0d+00
   real(8), parameter :: solu_ds2 = 1.0d+00
   real(8), parameter :: solu_ssa = 1.0d+00
   real(8), parameter :: solu_ssc = 1.0d+00
   real(8), parameter :: solu_sss = 1.0d+00
   real(8), parameter :: solu_occ = 0.4d+00
   real(8), parameter :: solu_bc1 = 0.4d+00
   real(8), parameter :: solu_bc2 = 0.8d+00
   real(8), parameter :: solu_bc3 = 1.0d+00
   real(8), parameter :: solu_dbc = 0.0d+00
   real(8), parameter :: solu_boc = 0.6d+00
   real(8), parameter :: solu_bcs = 1.0d+00
   real(8), parameter :: solu_ocs = 1.0d+00
   real(8), parameter :: solu_mxx = 1.0d+00
!----------------------------------------------------------------------------------------------------------------------
!
!  MODEL PARAMETERS AND VARIABLES THAT PROBABLY DO NOT NEED TO BE CHANGED.
!
!----------------------------------------------------------------------------------------------------------------------
   integer, parameter :: ngases     = 3      ! number of gas-phase species
   integer, parameter :: nmass_spcs = 5      ! total number of mass species
   integer, parameter :: gas_h2so4  = 1      !\
   integer, parameter :: gas_hno3   = 2      !-indices in the gas array
   integer, parameter :: gas_nh3    = 3      !/
   integer, parameter :: prod_index_sulf = 1 ! sulf index in prod_index(:,:)
   integer, parameter :: prod_index_bcar = 2 ! bcar index in prod_index(:,:)
   integer, parameter :: prod_index_ocar = 3 ! ocar index in prod_index(:,:)
   integer, parameter :: prod_index_dust = 4 ! dust index in prod_index(:,:)
   integer, parameter :: prod_index_seas = 5 ! seas index in prod_index(:,:)
   !-------------------------------------------------------------------------------------------------------------------
   ! EMIS_DENS_XXXX is the dry particle density of emitted species XXXX.
   !
   ! These are used only for deriving number emission rates from mass emission rates.
   !
   ! For sulfate, the emissions are treated as pure dry ammonium sulfate.
   ! For sea salt, the emissions are treated as pure dry sodium chloride.
   !
   ! For the sulfate density, volume is converted to ammonium sulfate mass using
   ! the density of ammonium sulfate, which is then converted to sulfate (only) mass.
   !-------------------------------------------------------------------------------------------------------------------
   real(8), parameter :: emis_dens_sulf = 1.77d+00  & ! [gso4/cm^3] - nh42so4
                                        * mw_so4 / mw_nh42so4
   real(8), parameter :: emis_dens_bcar = 1.70d+00   ! [g/cm^3] - ghan et al. (2001) - mirage
   real(8), parameter :: emis_dens_ocar = 1.00d+00   ! [g/cm^3] - ghan et al. (2001) - mirage
   real(8), parameter :: emis_dens_dust = 2.60d+00   ! [g/cm^3] - ghan et al. (2001) - mirage
   real(8), parameter :: emis_dens_seas = 2.165d+00  ! [g/cm^3] - nacl
   real(8), parameter :: emis_dens_bocc = 0.50d+00 &  ! [g/cm^3] - average
                                       * ( emis_dens_bcar + emis_dens_ocar )
   real, dimension(nemis_spcs) :: emis_dens = (/  emis_dens_sulf, &
                 emis_dens_sulf, emis_dens_bcar, emis_dens_ocar, &
                 emis_dens_dust, emis_dens_seas, emis_dens_seas, &
                 emis_dens_bocc, emis_dens_bocc, emis_dens_dust /)
   !-------------------------------------------------------------------------------------------------------------------
   ! The aerosol chemical species are SO4, BC, OC, mineral dust, and sea salt.
   ! Nitrate, ammonium and water are not included here.
   !-------------------------------------------------------------------------------------------------------------------
   character(len=4) :: chem_spc_name(nmass_spcs) &
                   = (/'SULF','BCAR','OCAR','DUST','SEAS'/)
   !srf-kml CHARACTER(LEN=3) :: gas_name(3)=(/'NO3','NH4','H2O'/)
   !-------------------------------------------------------------------------------------------------------------------
   ! The Maximum Inorganic Mass Ratio (MIMR) in modes DD1, DD2, BC1, and BC2.
   !
   ! The above MIVF values are converted to MIMR values for computational efficiency.
   ! The volume of inorganic coating is converted to mass using the default ambient aerosol density.
   !-------------------------------------------------------------------------------------------------------------------
   real(8), parameter :: mimr_ddd = ( densp / emis_dens_dust ) &
                                 * mivf_ddd / ( 1.0d+00 - mivf_ddd )
   real(8), parameter :: mimr_bc1 = ( densp / emis_dens_bcar ) &
                                 * mivf_bc1 / ( 1.0d+00 - mivf_bc1 )
   real(8), parameter :: mimr_bc2 = ( densp / emis_dens_bcar ) &
                                 * mivf_bc2 / ( 1.0d+00 - mivf_bc2 )
!----------------------------------------------------------------------------------------------------------------------
!
!  VARIABLES HELD IN THE MANNER OF A COMMON BLOCK.
!
!----------------------------------------------------------------------------------------------------------------------
   integer, save :: ixxx, iyyy, ilay    ! current grid cell indices
   logical, save :: include_bc3         ! true if mechanism includes mode bc3; false otherwise
!----------------------------------------------------------------------------------------------------------------------
!  Aerosol mode names (and numbers) that might appear in one or more mechanisms.
!  These mode number only pertain to this set of all possible modes, and are not the mode numbers
!  used for any specific mechanism.
!----------------------------------------------------------------------------------------------------------------------
   integer, parameter :: nmodes_max=18
   character(len=3) :: mname(nmodes_max)
   ! Mode #     1     2     3     4     5     6     7     8     9    ! # to identify the mode in MODES1, etc. below.
   data mname/'AKK','ACC','DD1','DS1','DD2','DS2','SSA','SSC','SSS', &
              'OCC','BC1','BC2','BC3','OCS','DBC','BOC','BCS','MXX'/
   ! Mode #     10    11    12    13    14    15    16    17    18   ! # to identify the mode in MODES1, etc. below.
!----------------------------------------------------------------------------------------------------------------------
!  Aerosol species defined for each mode in each mechanism.
!----------------------------------------------------------------------------------------------------------------------
   integer, save :: mspcs(nmass_spcs,nmodes_max)
   data mspcs(1,1:nmodes_max)/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/   ! SULF: =0 no sulfate, =1 has sulfate
   data mspcs(2,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,1/   ! BCAR: =0 no BC     , =1 has BC
   data mspcs(3,1:nmodes_max)/0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,1/   ! OCAR: =0 no OC     , =1 has OC
   data mspcs(4,1:nmodes_max)/0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,0,0,1/   ! DUST: =0 no dust   , =1 has dust
   data mspcs(5,1:nmodes_max)/0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1/   ! SEAS: =0 no seasalt, =1 has seasalt
!----------------------------------------------------------------------------------------------------------------------
!  Aerosol modes used for each mechanism.
!----------------------------------------------------------------------------------------------------------------------
   integer, parameter :: nm1=16,nm2=16,nm3=13,nm4=10
   integer, parameter :: nm5=14,nm6=14,nm7=11,nm8=8
   integer :: modes1(nm1),modes2(nm2),modes3(nm3),modes4(nm4)
   integer :: modes5(nm5),modes6(nm6),modes7(nm7),modes8(nm8)
   data modes1/ 1, 2, 3, 4, 5, 6, 7, 8,10,11,12,13,15,16,17,18/
   data modes2/ 1, 2, 3, 4, 5, 6, 7, 8,10,11,12,14,15,16,17,18/
   data modes3/ 1, 2, 3, 4, 5, 6, 7, 8,10,11,12,16,18/
   data modes4/ 2, 3, 4, 5, 6, 9,10,11,12,18/
   data modes5/ 1, 2, 3, 4, 7, 8,10,11,12,13,15,16,17,18/
   data modes6/ 1, 2, 3, 4, 7, 8,10,11,12,14,15,16,17,18/
   data modes7/ 1, 2, 3, 4, 7, 8,10,11,12,16,18/
   data modes8/ 2, 3, 4, 9,10,11,12,18/

   !Comments from Luiz Flavio:
   !To understood: modes8 (mech 8) has 2=ACC, 3=DD1, 4=DS1, 9=SSS, 10=OCC, 11=BC1, 12=BC2 and 18=MXX
   !If You put this configuration inside an specie, e.g., SEAS (see mspcs) You got 9 and 18. (SEAS=5)
   ! at the same way, for DUST (4), You got 3, 4, 18, i.e., DD1, DS1 and MXX respectively.
   !Please, see that modes 1 (mech 1)  has all the modes, 1 to 18, but You can see that not all species has all modes ON.
   !

!----------------------------------------------------------------------------------------------------------------------
!  LFR-Aerosol sources used for each mechanism.
!----------------------------------------------------------------------------------------------------------------------
   integer, save :: msource(nmass_spcs,nmodes_max)
   data msource(1,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! SULF: =0 no source, =1 has source
   data msource(2,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! BCAR
   data msource(3,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! OCAR
   data msource(4,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! DUST
   data msource(5,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! SEAS

!----------------------------------------------------------------------------------------------------------------------
!  LFR-Aerosol dry deposition used for each mechanism.
!----------------------------------------------------------------------------------------------------------------------
   integer, save :: mdrydep(nmass_spcs,nmodes_max)
   data mdrydep(1,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! SULF: =0 no drydep, =1 has drydep
   data mdrydep(2,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! BCAR
   data mdrydep(3,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! OCAR
   data mdrydep(4,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! DUST
   data mdrydep(5,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! SEAS

!----------------------------------------------------------------------------------------------------------------------
!  LFR-Aerosol wet deposition used for each mechanism.
!----------------------------------------------------------------------------------------------------------------------
   integer, save :: mwetdep(nmass_spcs,nmodes_max)
   data mwetdep(1,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! SULF: =0 no wetdep, =1 has wetdep
   data mwetdep(2,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! BCAR
   data mwetdep(3,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! OCAR
   data mwetdep(4,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! DUST
   data mwetdep(5,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! SEAS

!----------------------------------------------------------------------------------------------------------------------
!  LFR-Aerosol fdda used for each mechanism.
!----------------------------------------------------------------------------------------------------------------------
   integer, save :: mfdda(nmass_spcs,nmodes_max)
   data mfdda(1,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! SULF: =0 no fdda, =1 has fdda
   data mfdda(2,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! BCAR
   data mfdda(3,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! OCAR
   data mfdda(4,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! DUST
   data mfdda(5,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! SEAS

!----------------------------------------------------------------------------------------------------------------------
!  LFR-Aerosol offline emission used for each mechanism.
!----------------------------------------------------------------------------------------------------------------------
   integer, save :: moffem(nmass_spcs,nmodes_max)
   data moffem(1,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! SULF: =0 no offline emission, =1 has offline emission
   data moffem(2,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! BCAR
   data moffem(3,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! OCAR
   data moffem(4,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! DUST
   data moffem(5,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! SEAS

!----------------------------------------------------------------------------------------------------------------------
!  LFR-Aerosol transport used for each mechanism.
!----------------------------------------------------------------------------------------------------------------------
   integer, save :: mtransp(nmass_spcs,nmodes_max)
   data mtransp(1,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! SULF: =0 no transport, =1 has transport
   data mtransp(2,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! BCAR
   data mtransp(3,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! OCAR
   data mtransp(4,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! DUST
   data mtransp(5,1:nmodes_max)/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/   ! SEAS


!----------------------------------------------------------------------------------------------------------------------
!  Indices of the AERO array. There are 78 possible indices.
!----------------------------------------------------------------------------------------------------------------------
   integer       :: mass_no3=1, mass_nh4=2, mass_h2o=3
   integer, save :: numb_akk_1, numb_akk_2, mass_akk_sulf,  &
                    numb_acc_1, numb_acc_2, mass_acc_sulf,  &
                    numb_dd1_1, numb_dd1_2, mass_dd1_sulf, mass_dd1_dust,  &
                    numb_ds1_1, numb_ds1_2, mass_ds1_sulf, mass_ds1_dust,  &
                    numb_dd2_1, numb_dd2_2, mass_dd2_sulf, mass_dd2_dust,  &
                    numb_ds2_1, numb_ds2_2, mass_ds2_sulf, mass_ds2_dust,  &
                    numb_ssa_1, numb_ssa_2, mass_ssa_sulf, mass_ssa_seas,  &
                    numb_ssc_1, numb_ssc_2, mass_ssc_sulf, mass_ssc_seas,  &
                    numb_sss_1, numb_sss_2, mass_sss_sulf, mass_sss_seas,  &
                    numb_occ_1, numb_occ_2, mass_occ_sulf, mass_occ_ocar,  &
                    numb_bc1_1, numb_bc1_2, mass_bc1_sulf, mass_bc1_bcar,  &
                    numb_bc2_1, numb_bc2_2, mass_bc2_sulf, mass_bc2_bcar,  &
                    numb_bc3_1, numb_bc3_2, mass_bc3_sulf, mass_bc3_bcar,  &
                    numb_ocs_1, numb_ocs_2, mass_ocs_sulf, mass_ocs_ocar,  &
                    numb_dbc_1, numb_dbc_2, mass_dbc_sulf, mass_dbc_bcar, mass_dbc_dust,  &
                    numb_boc_1, numb_boc_2, mass_boc_sulf, mass_boc_bcar, mass_boc_ocar,  &
                    numb_bcs_1, numb_bcs_2, mass_bcs_sulf, mass_bcs_bcar,  &
                    numb_mxx_1, numb_mxx_2, mass_mxx_sulf, mass_mxx_bcar, mass_mxx_ocar,  &
                    mass_mxx_dust, mass_mxx_seas
                                                                    !   set to unity.




   !From tramp_config


!-------------------------------------------------------------------------------------------------------------------------
!     3. Select modes to undergo condensational growth for the desired mechanism (1-8). (Ignore other mechanisms.)
!        ICONDn(I)=1, condensational growth done; ICONDn(I)=0, condensational growth not done.
!        Ordinarily, all modes would undergo condenational growth.
!-------------------------------------------------------------------------------------------------------------------------
      integer, save, dimension(nm1) :: icond1=(/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/)  ! mechanism 1
      integer, save, dimension(nm2) :: icond2=(/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/)  ! mechanism 2
      integer, save, dimension(nm3) :: icond3=(/1,1,1,1,1,1,1,1,1,1,1,1,1/)        ! mechanism 3
      integer, save, dimension(nm4) :: icond4=(/1,1,1,1,1,1,1,1,1,1/)              ! mechanism 4
      integer, save, dimension(nm5) :: icond5=(/1,1,1,1,1,1,1,1,1,1,1,1,1,1/)      ! mechanism 5
      integer, save, dimension(nm6) :: icond6=(/1,1,1,1,1,1,1,1,1,1,1,1,1,1/)      ! mechanism 6
      integer, save, dimension(nm7) :: icond7=(/1,1,1,1,1,1,1,1,1,1,1/)            ! mechanism 7
      integer, save, dimension(nm8) :: icond8=(/1,1,1,1,1,1,1,1/)                  ! mechanism 8
!-------------------------------------------------------------------------------------------------------------------------------------
!     4. Optionally edit the table of coagulation interactions.
!        The donor modes may not be modified. Each receptor mode must contain all species present in either donor mode.
!        Entering 'OFF' for the receptor mode name disables coagulation between the two donor modes.
!-------------------------------------------------------------------------------------------------------------------------------------
!
!     Mechanism 1.
!
!     FIRST MODE               AKK   ACC   DD1   DS1   DD2   DS2   SSA   SSC   OCC   BC1   BC2   BC3   DBC   BOC   BCS   MXX    SECOND
!                                                                                                                                MODE
!-------------------------------------------------------------------------------------------------------------------------------------
      character(len=3) :: citable1(nm1,nm1)
      data citable1(1:nm1, 1)/'AKK','ACC','DD1','DS1','DD2','DS2','SSA','SSC','OCC','BCS','BCS','BCS','DBC','BOC','BCS','MXX'/ ! AKK
      data citable1(1:nm1, 2)/'ACC','ACC','DD1','DS1','DD2','DS2','SSA','SSC','OCC','BCS','BCS','BCS','DBC','BOC','BCS','MXX'/ ! ACC
      data citable1(1:nm1, 3)/'DD1','DD1','DD1','DD1','DD2','DD2','MXX','MXX','MXX','DBC','DBC','DBC','DBC','MXX','DBC','MXX'/ ! DD1
      data citable1(1:nm1, 4)/'DS1','DS1','DD1','DS1','DD2','DS2','MXX','MXX','MXX','DBC','DBC','DBC','DBC','MXX','DBC','MXX'/ ! DS1
      data citable1(1:nm1, 5)/'DD2','DD2','DD2','DD2','DD2','DD2','MXX','MXX','MXX','DBC','DBC','DBC','DBC','MXX','DBC','MXX'/ ! DD2
      data citable1(1:nm1, 6)/'DS2','DS2','DD2','DS2','DD2','DS2','MXX','MXX','MXX','DBC','DBC','DBC','DBC','MXX','DBC','MXX'/ ! DS2
      data citable1(1:nm1, 7)/'SSA','SSA','MXX','MXX','MXX','MXX','SSA','SSC','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! SSA
      data citable1(1:nm1, 8)/'SSC','SSC','MXX','MXX','MXX','MXX','SSC','SSC','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! SSC
      data citable1(1:nm1, 9)/'OCC','OCC','MXX','MXX','MXX','MXX','MXX','MXX','OCC','BOC','BOC','BOC','MXX','BOC','BOC','MXX'/ ! OCC
      data citable1(1:nm1,10)/'BCS','BCS','DBC','DBC','DBC','DBC','MXX','MXX','BOC','BC1','BC1','BC1','DBC','BOC','BCS','MXX'/ ! BC1
      data citable1(1:nm1,11)/'BCS','BCS','DBC','DBC','DBC','DBC','MXX','MXX','BOC','BC1','BC2','BC2','DBC','BOC','BCS','MXX'/ ! BC2
      data citable1(1:nm1,12)/'BCS','BCS','DBC','DBC','DBC','DBC','MXX','MXX','BOC','BC1','BC2','BC3','DBC','BOC','BCS','MXX'/ ! BC3
      data citable1(1:nm1,13)/'DBC','DBC','DBC','DBC','DBC','DBC','MXX','MXX','MXX','DBC','DBC','DBC','DBC','MXX','DBC','MXX'/ ! DBC
      data citable1(1:nm1,14)/'BOC','BOC','MXX','MXX','MXX','MXX','MXX','MXX','BOC','BOC','BOC','BOC','MXX','BOC','BOC','MXX'/ ! BOC
      data citable1(1:nm1,15)/'BCS','BCS','DBC','DBC','DBC','DBC','MXX','MXX','BOC','BCS','BCS','BCS','DBC','BOC','BCS','MXX'/ ! BCS
      data citable1(1:nm1,16)/'MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! MXX
!-------------------------------------------------------------------------------------------------------------------------------------
!
!     Mechanism 2.
!
!     FIRST MODE               AKK   ACC   DD1   DS1   DD2   DS2   SSA   SSC   OCC   BC1   BC2   OCS   DBC   BOC   BCS   MXX    SECOND
!                                                                                                                                MODE
!-------------------------------------------------------------------------------------------------------------------------------------
      character(len=3) :: citable2(nm2,nm2)
      data citable2(1:nm2, 1)/'AKK','ACC','DD1','DS1','DD2','DS2','SSA','SSC','OCS','BCS','BCS','OCS','DBC','BOC','BCS','MXX'/ ! AKK
      data citable2(1:nm2, 2)/'ACC','ACC','DD1','DS1','DD2','DS2','SSA','SSC','OCS','BCS','BCS','OCS','DBC','BOC','BCS','MXX'/ ! ACC
      data citable2(1:nm2, 3)/'DD1','DD1','DD1','DD1','DD2','DD2','MXX','MXX','MXX','DBC','DBC','MXX','DBC','MXX','DBC','MXX'/ ! DD1
      data citable2(1:nm2, 4)/'DS1','DS1','DD1','DS1','DD2','DS2','MXX','MXX','MXX','DBC','DBC','MXX','DBC','MXX','DBC','MXX'/ ! DS1
      data citable2(1:nm2, 5)/'DD2','DD2','DD2','DD2','DD2','DD2','MXX','MXX','MXX','DBC','DBC','MXX','DBC','MXX','DBC','MXX'/ ! DD2
      data citable2(1:nm2, 6)/'DS2','DS2','DD2','DS2','DD2','DS2','MXX','MXX','MXX','DBC','DBC','MXX','DBC','MXX','DBC','MXX'/ ! DS2
      data citable2(1:nm2, 7)/'SSA','SSA','MXX','MXX','MXX','MXX','SSA','SSC','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! SSA
      data citable2(1:nm2, 8)/'SSC','SSC','MXX','MXX','MXX','MXX','SSC','SSC','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! SSC
      data citable2(1:nm2, 9)/'OCS','OCS','MXX','MXX','MXX','MXX','MXX','MXX','OCC','BOC','BOC','OCS','MXX','BOC','BOC','MXX'/ ! OCC
      data citable2(1:nm2,10)/'BCS','BCS','DBC','DBC','DBC','DBC','MXX','MXX','BOC','BC1','BC1','BOC','DBC','BOC','BCS','MXX'/ ! BC1
      data citable2(1:nm2,11)/'BCS','BCS','DBC','DBC','DBC','DBC','MXX','MXX','BOC','BC1','BC2','BOC','DBC','BOC','BCS','MXX'/ ! BC2
      data citable2(1:nm2,12)/'OCS','OCS','MXX','MXX','MXX','MXX','MXX','MXX','OCS','BOC','BOC','OCS','MXX','BOC','BOC','MXX'/ ! OCS
      data citable2(1:nm2,13)/'DBC','DBC','DBC','DBC','DBC','DBC','MXX','MXX','MXX','DBC','DBC','MXX','DBC','MXX','DBC','MXX'/ ! DBC
      data citable2(1:nm2,14)/'BOC','BOC','MXX','MXX','MXX','MXX','MXX','MXX','BOC','BOC','BOC','BOC','MXX','BOC','BOC','MXX'/ ! BOC
      data citable2(1:nm2,15)/'BCS','BCS','DBC','DBC','DBC','DBC','MXX','MXX','BOC','BCS','BCS','BOC','DBC','BOC','BCS','MXX'/ ! BCS
      data citable2(1:nm2,16)/'MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! MXX
!-------------------------------------------------------------------------------------------------------------------------------------
!
!     Mechanism 3.
!
!     FIRST MODE               AKK   ACC   DD1   DS1   DD2   DS2   SSA   SSC   OCC   BC1   BC2   BOC   MXX    SECOND
!                                                                                                              MODE
!-------------------------------------------------------------------------------------------------------------------------------------
      character(len=3) :: citable3(nm3,nm3)
      data citable3(1:nm3, 1)/'AKK','ACC','DD1','DS1','DD2','DS2','SSA','SSC','OCC','BC1','BC2','BOC','MXX'/ ! AKK
      data citable3(1:nm3, 2)/'ACC','ACC','DD1','DS1','DD2','DS2','SSA','SSC','OCC','BC1','BC2','BOC','MXX'/ ! ACC
      data citable3(1:nm3, 3)/'DD1','DD1','DD1','DD1','DD2','DD2','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! DD1
      data citable3(1:nm3, 4)/'DS1','DS1','DD1','DS1','DD2','DS2','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! DS1
      data citable3(1:nm3, 5)/'DD2','DD2','DD2','DD2','DD2','DD2','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! DD2
      data citable3(1:nm3, 6)/'DS2','DS2','DD2','DS2','DD2','DS2','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! DS2
      data citable3(1:nm3, 7)/'SSA','SSA','MXX','MXX','MXX','MXX','SSA','SSC','MXX','MXX','MXX','MXX','MXX'/ ! SSA
      data citable3(1:nm3, 8)/'SSC','SSC','MXX','MXX','MXX','MXX','SSC','SSC','MXX','MXX','MXX','MXX','MXX'/ ! SSC
      data citable3(1:nm3, 9)/'OCC','OCC','MXX','MXX','MXX','MXX','MXX','MXX','OCC','BOC','BOC','BOC','MXX'/ ! OCC
      data citable3(1:nm3,10)/'BC1','BC1','MXX','MXX','MXX','MXX','MXX','MXX','BOC','BC1','BC1','BOC','MXX'/ ! BC1
      data citable3(1:nm3,11)/'BC2','BC2','MXX','MXX','MXX','MXX','MXX','MXX','BOC','BC1','BC2','BOC','MXX'/ ! BC2
      data citable3(1:nm3,12)/'BOC','BOC','MXX','MXX','MXX','MXX','MXX','MXX','BOC','BOC','BOC','BOC','MXX'/ ! BOC
      data citable3(1:nm3,13)/'MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! MXX
!-------------------------------------------------------------------------------------------------------------------------------------
!
!     Mechanism 4.
!
!     FIRST MODE               ACC   DD1   DS1   DD2   DS2   SSS   OCC   BC1   BC2   MXX    SECOND
!                                                                                            MODE
!-------------------------------------------------------------------------------------------------------------------------------------
      character(len=3) :: citable4(nm4,nm4)
      data citable4(1:nm4, 1)/'ACC','DD1','DS1','DD2','DS2','SSS','OCC','BC1','BC2','MXX'/ ! ACC
      data citable4(1:nm4, 2)/'DD1','DD1','DD1','DD2','DD2','MXX','MXX','MXX','MXX','MXX'/ ! DD1
      data citable4(1:nm4, 3)/'DS1','DD1','DS1','DD2','DS2','MXX','MXX','MXX','MXX','MXX'/ ! DS1
      data citable4(1:nm4, 4)/'DD2','DD2','DD2','DD2','DD2','MXX','MXX','MXX','MXX','MXX'/ ! DD2
      data citable4(1:nm4, 5)/'DS2','DD2','DS2','DD2','DS2','MXX','MXX','MXX','MXX','MXX'/ ! DS2
      data citable4(1:nm4, 6)/'SSS','MXX','MXX','MXX','MXX','SSS','MXX','MXX','MXX','MXX'/ ! SSS
      data citable4(1:nm4, 7)/'OCC','MXX','MXX','MXX','MXX','MXX','OCC','MXX','MXX','MXX'/ ! OCC
      data citable4(1:nm4, 8)/'BC1','MXX','MXX','MXX','MXX','MXX','MXX','BC1','BC1','MXX'/ ! BC1
      data citable4(1:nm4, 9)/'BC2','MXX','MXX','MXX','MXX','MXX','MXX','BC1','BC2','MXX'/ ! BC2
      data citable4(1:nm4,10)/'MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! MXX
!-------------------------------------------------------------------------------------------------------------------------
!
!     Mechanism 5.
!
!     FIRST MODE               AKK   ACC   DD1   DS1   SSA   SSC   OCC   BC1   BC2   BC3   DBC   BOC   BCS   MXX    SECOND
!                                                                                                                    MODE
!-------------------------------------------------------------------------------------------------------------------------
      character(len=3) :: citable5(nm5,nm5)
      data citable5(1:nm5, 1)/'AKK','ACC','DD1','DS1','SSA','SSC','OCC','BCS','BCS','BCS','DBC','BOC','BCS','MXX'/ ! AKK
      data citable5(1:nm5, 2)/'ACC','ACC','DD1','DS1','SSA','SSC','OCC','BCS','BCS','BCS','DBC','BOC','BCS','MXX'/ ! ACC
      data citable5(1:nm5, 3)/'DD1','DD1','DD1','DD1','MXX','MXX','MXX','DBC','DBC','DBC','DBC','MXX','DBC','MXX'/ ! DD1
      data citable5(1:nm5, 4)/'DS1','DS1','DD1','DS1','MXX','MXX','MXX','DBC','DBC','DBC','DBC','MXX','DBC','MXX'/ ! DS1
      data citable5(1:nm5, 5)/'SSA','SSA','MXX','MXX','SSA','SSC','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! SSA
      data citable5(1:nm5, 6)/'SSC','SSC','MXX','MXX','SSC','SSC','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! SSC
      data citable5(1:nm5, 7)/'OCC','OCC','MXX','MXX','MXX','MXX','OCC','BOC','BOC','BOC','MXX','BOC','BOC','MXX'/ ! OCC
      data citable5(1:nm5, 8)/'BCS','BCS','DBC','DBC','MXX','MXX','BOC','BC1','BC1','BC1','DBC','BOC','BCS','MXX'/ ! BC1
      data citable5(1:nm5, 9)/'BCS','BCS','DBC','DBC','MXX','MXX','BOC','BC1','BC2','BC2','DBC','BOC','BCS','MXX'/ ! BC2
      data citable5(1:nm5,10)/'BCS','BCS','DBC','DBC','MXX','MXX','BOC','BC1','BC2','BC3','DBC','BOC','BCS','MXX'/ ! BC3
      data citable5(1:nm5,11)/'DBC','DBC','DBC','DBC','MXX','MXX','MXX','DBC','DBC','DBC','DBC','MXX','DBC','MXX'/ ! DBC
      data citable5(1:nm5,12)/'BOC','BOC','MXX','MXX','MXX','MXX','BOC','BOC','BOC','BOC','MXX','BOC','BOC','MXX'/ ! BOC
      data citable5(1:nm5,13)/'BCS','BCS','DBC','DBC','MXX','MXX','BOC','BCS','BCS','BCS','DBC','BOC','BCS','MXX'/ ! BCS
      data citable5(1:nm5,14)/'MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! MXX
!-------------------------------------------------------------------------------------------------------------------------
!
!     Mechanism 6.
!
!     FIRST MODE               AKK   ACC   DD1   DS1   SSA   SSC   OCC   BC1   BC2   OCS   DBC   BOC   BCS   MXX    SECOND
!                                                                                                                    MODE
!-------------------------------------------------------------------------------------------------------------------------
      character(len=3) :: citable6(nm6,nm6)
      data citable6(1:nm6, 1)/'AKK','ACC','DD1','DS1','SSA','SSC','OCS','BCS','BCS','OCS','DBC','BOC','BCS','MXX'/ ! AKK
      data citable6(1:nm6, 2)/'ACC','ACC','DD1','DS1','SSA','SSC','OCS','BCS','BCS','OCS','DBC','BOC','BCS','MXX'/ ! ACC
      data citable6(1:nm6, 3)/'DD1','DD1','DD1','DD1','MXX','MXX','MXX','DBC','DBC','MXX','DBC','MXX','DBC','MXX'/ ! DD1
      data citable6(1:nm6, 4)/'DS1','DS1','DD1','DS1','MXX','MXX','MXX','DBC','DBC','MXX','DBC','MXX','DBC','MXX'/ ! DS2
      data citable6(1:nm6, 5)/'SSA','SSA','MXX','MXX','SSA','SSC','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! SSA
      data citable6(1:nm6, 6)/'SSC','SSC','MXX','MXX','SSC','SSC','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! SSC
      data citable6(1:nm6, 7)/'OCS','OCS','MXX','MXX','MXX','MXX','OCC','BOC','BOC','OCS','MXX','BOC','BOC','MXX'/ ! OCC
      data citable6(1:nm6, 8)/'BCS','BCS','DBC','DBC','MXX','MXX','BOC','BC1','BC1','BOC','DBC','BOC','BCS','MXX'/ ! BC1
      data citable6(1:nm6, 9)/'BCS','BCS','DBC','DBC','MXX','MXX','BOC','BC1','BC2','BOC','DBC','BOC','BCS','MXX'/ ! BC2
      data citable6(1:nm6,10)/'OCS','OCS','MXX','MXX','MXX','MXX','OCS','BOC','BOC','OCS','MXX','MXX','MXX','MXX'/ ! OCS
      data citable6(1:nm6,11)/'DBC','DBC','DBC','DBC','MXX','MXX','MXX','DBC','DBC','MXX','DBC','MXX','MXX','MXX'/ ! DBC
      data citable6(1:nm6,12)/'BOC','BOC','MXX','MXX','MXX','MXX','BOC','BOC','BOC','MXX','MXX','BOC','MXX','MXX'/ ! BOC
      data citable6(1:nm6,13)/'BCS','BCS','DBC','DBC','MXX','MXX','BOC','BCS','BCS','MXX','MXX','MXX','BCS','MXX'/ ! BCS
      data citable6(1:nm6,14)/'MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! MXX
!-------------------------------------------------------------------------------------------------------------------------
!
!     Mechanism 7.
!
!     FIRST MODE               AKK   ACC   DD1   DS1   SSA   SSC   OCC   BC1   BC2   BOC  MXX     SECOND
!                                                                                                  MODE
!-------------------------------------------------------------------------------------------------------------------------
      character(len=3) :: citable7(nm7,nm7)
      data citable7(1:nm7, 1)/'AKK','ACC','DD1','DS1','SSA','SSC','OCC','BC1','BC2','BOC','MXX'/ ! AKK
      data citable7(1:nm7, 2)/'ACC','ACC','DD1','DS1','SSA','SSC','OCC','BC1','BC2','BOC','MXX'/ ! ACC
      data citable7(1:nm7, 3)/'DD1','DD1','DD1','DD1','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! DD1
      data citable7(1:nm7, 4)/'DS1','DS1','DD1','DS1','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! DS1
      data citable7(1:nm7, 5)/'SSA','SSA','MXX','MXX','SSA','SSC','MXX','MXX','MXX','MXX','MXX'/ ! SSA
      data citable7(1:nm7, 6)/'SSC','SSC','MXX','MXX','SSC','SSC','MXX','MXX','MXX','MXX','MXX'/ ! SSC
      data citable7(1:nm7, 7)/'OCC','OCC','MXX','MXX','MXX','MXX','OCC','BOC','BOC','BOC','MXX'/ ! OCC
      data citable7(1:nm7, 8)/'BC1','BC1','MXX','MXX','MXX','MXX','BOC','BC1','BC1','BOC','MXX'/ ! BC1
      data citable7(1:nm7, 9)/'BC2','BC2','MXX','MXX','MXX','MXX','BOC','BC1','BC2','BOC','MXX'/ ! BC2
      data citable7(1:nm7,10)/'BOC','BOC','MXX','MXX','MXX','MXX','BOC','BOC','BOC','BOC','MXX'/ ! BOC
      data citable7(1:nm7,11)/'MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! MXX
!-------------------------------------------------------------------------------------------------------------------------
!
!     Mechanism 8.
!
!     FIRST MODE               ACC   DD1   DS1   SSS   OCC   BC1   BC2   MXX    SECOND
!                                                                                MODE
!-------------------------------------------------------------------------------------------------------------------------
      character(len=3) :: citable8(nm8,nm8)
      data citable8(1:nm8, 1)/'ACC','DD1','DS1','SSS','OCC','BC1','BC2','MXX'/ ! ACC
      data citable8(1:nm8, 2)/'DD1','DD1','DD1','MXX','MXX','MXX','MXX','MXX'/ ! DD1
      data citable8(1:nm8, 3)/'DS1','DD1','DS1','MXX','MXX','MXX','MXX','MXX'/ ! DS1
      data citable8(1:nm8, 4)/'SSS','MXX','MXX','SSS','MXX','MXX','MXX','MXX'/ ! SSS
      data citable8(1:nm8, 5)/'OCC','MXX','MXX','MXX','OCC','MXX','MXX','MXX'/ ! OCC
      data citable8(1:nm8, 6)/'BC1','MXX','MXX','MXX','MXX','BC1','BC1','MXX'/ ! BC1
      data citable8(1:nm8, 7)/'BC2','MXX','MXX','MXX','MXX','BC1','BC2','MXX'/ ! BC2
      data citable8(1:nm8, 8)/'MXX','MXX','MXX','MXX','MXX','MXX','MXX','MXX'/ ! MXX








   !From TrampSetup

      integer, parameter :: n_dp_condtable = 1000              ! [1] number of tabulated particle diameters
      real(8), parameter :: dp_condtable_min = dpmin_global    ! [m] minimum ambient particle diameter
      real(8), parameter :: dp_condtable_max = dpmax_global    ! [m] maximum ambient particle diameter

      integer :: mode_numb_akk, mode_numb_acc, mode_numb_dd1, mode_numb_dd2
      integer :: mode_numb_ds1, mode_numb_ds2, mode_numb_ssa, mode_numb_ssc
      integer :: mode_numb_sss, mode_numb_occ, mode_numb_bc1, mode_numb_bc2
      integer :: mode_numb_bc3, mode_numb_dbc, mode_numb_boc, mode_numb_bcs
      integer :: mode_numb_ocs, mode_numb_mxx
      !-------------------------------------------------------------------------
      ! nmodes_xxxx is the number of modes containing species xxxx.
      ! nmodes_seas is the number of modes containing sea salt, not the
      !             number of sea salt modes.
      ! number_of_seasalt_modes is the number of sea salt modes, either 1 (sss)
      !                         or 2 (ssa and ssc).
      !-------------------------------------------------------------------------
      integer :: nmodes_sulf, nmodes_bcar, nmodes_ocar, nmodes_dust, nmodes_seas
      integer :: number_of_seasalt_modes
      integer, allocatable :: numb_map(:)                           ! [1]
      integer, allocatable :: mass_map(:,:)                ! [1]
      integer, allocatable :: prod_index(:,:)              ! [1]
      integer, allocatable :: sulf_map(:)                     ! [1]
      integer, allocatable :: bcar_map(:)                     ! [1]
      integer, allocatable :: ocar_map(:)                     ! [1]
      integer, allocatable :: dust_map(:)                     ! [1]
      integer, allocatable :: seas_map(:)                     ! [1]
      integer, allocatable :: seas_mode_map(:)                ! [1]
      integer, allocatable :: seas_mode_mass_map(:)           ! [1]
      integer, allocatable :: mode_numb_seas(:)               ! [1]
      integer, allocatable :: nm(:)                                 ! [1]
     ! character(len=4), save :: nm_spc_name(:,:)
      character(len=4), allocatable  :: nm_spc_name(:,:)

      type giklq_type
        integer :: n
        integer, allocatable :: l(:), k(:)
        integer, allocatable :: qq(:)
      end type giklq_type

      type dikl_type
        integer :: i
        integer :: k
        integer :: l
      end type dikl_type

      type (giklq_type), allocatable :: giklq_control(:)
      type (dikl_type), allocatable :: dikl_control(:)
      integer, save :: ndikl

      integer, allocatable :: giklq(:,:,:,:) ! [1]
      integer, allocatable :: dikl (:,:,:)            ! [1]
      integer, allocatable :: dij  (:,:)                     ! [1]
      real(8), allocatable :: xdij (:,:) ! [1]
      real(8), allocatable :: recip_part_mass(:)                  ! [1/ug]
      real(8), allocatable :: kci_coef_dp     (:,:)             ! [m^3/s]
      real(8), allocatable :: kci_coef_dp_aeq1(:,:)             ! [m^3/s]
      real(8), allocatable :: theta_poly (:)                        ! [1]
      real(8), allocatable :: dp0(:)      ! default mode diameters of average mass [m]
                                    !   calculated from the dgn0 and sig0 values
      real(8), allocatable :: dp0_emis(:) ! mode diameters of average mass [m] for emissions,
                                          !   calculated from the dgn0_emis and sig0_emis values
      !-------------------------------------------------------------------------
      ! diffcoef_m2s(i) is the diffusivity of h2so4 in air [m^2/s] for layer l.
      !-------------------------------------------------------------------------
      real(8), allocatable :: diffcoef_m2s(:)
      !-------------------------------------------------------------------------
      ! kappai(i) is the activating fraction for mode i.
      !-------------------------------------------------------------------------
      real(8), allocatable :: kappai(:)

      !KML-------------------------------------------------------------------------
      ! kappapk(i) is the hygroscopicity parameter
      !-------------------------------------------------------------------------
      real(8), allocatable :: kappapk(:)

      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      ! denspi(i) is the default particle density for mode i.
      ! dens_comp(i) is the density of chemical component i.
      ! recip_dens_comp(i) is the reciprocal density of chemical component i.
      !-------------------------------------------------------------------------
      real(8), allocatable :: denspi(:)           ! [g/cm^3]
      real(8), allocatable :: dens_comp(:)        ! [g/cm^3]
      real(8), allocatable :: recip_dens_comp(:)  ! [cm^3/g]
      !-------------------------------------------------------------------------
      ! characteristic lognormal parameters for each mode: dgn0 [um], sig0 [1].
      !-------------------------------------------------------------------------
      real(8), allocatable :: dgn0(:)
      real(8), allocatable :: sig0(:)
      real(8), allocatable :: lnsig0(:)  ! ln( sig0 )
      !-------------------------------------------------------------------------
      ! lognormal parameters for emissions into each mode:
      !   dgn0_emis [um], sig0_emis [1].
      !   these are used to convert mass emission rates to number emission rates.
      !-------------------------------------------------------------------------
      real(8), allocatable :: dgn0_emis(:)
      real(8), allocatable :: sig0_emis(:)
      !-------------------------------------------------------------------------
      ! conv_dpam_to_dgn(i) converts the diameter of average mass to the
      ! geometric mean diameter of the number size distribution for mode i,
      ! based on an assumed standard deviation for each mode.
      !-------------------------------------------------------------------------
      real(8), allocatable :: conv_dpam_to_dgn(:)
      !-------------------------------------------------------------------------
      ! emis_mode_map and emis_spcs_map have elements corresponding to
      !   the aerosol types (in this order): akk(=1), acc(=2), bcc(=8), occ(=7),
      !               dd1(=3), ssa(=5), ssc(=6), boc(bc=8), boc(oc=9), dd2(=10).
      ! emis_mode_map(j) is mode number receiving the emissions held
      !                  in emis_mass(j).
      ! emis_spcs_map(j) is the chemical species number (1-5) of the chemical
      !                  species held in emis_mass(j).
      !-------------------------------------------------------------------------
      integer                        :: emis_mode_map(nemis_spcs)
      integer, dimension(nemis_spcs) :: emis_spcs_map = (/1,1,2,3,4,5,5,2,3,4/)
      !-------------------------------------------------------------------------
      ! the dimensions of these arrays depends upon mechanism.
      ! seas_map(i) is the mean mass per particle for sea salt mode i.
      !-------------------------------------------------------------------------
      character(len= 3), allocatable :: mode_name(:)
      integer          , allocatable :: mode_spcs(:,:)
      character(len=16), allocatable :: aero_spcs(:)
      integer          , allocatable :: icond(:)
      real             , allocatable :: recip_seas_mpp(:)         ! [1/ug]
      character(len=3) , allocatable :: citable(:,:)
      logical                        :: intermodal_transfer
      !-------------------------------------------------------------------------------------------------------------------
      ! diameter of average mass, averaged over all modes, used in the kk02 parameterization in subr. npfrate.
      !-------------------------------------------------------------------------------------------------------------------
      real(8) :: avg_dp_of_avg_mass_meters = 150.0d-09    ! [m] initial value; updated in subr. matrix.
      !-------------------------------------------------------------------------------------------------------------------
      ! variables for the lookup table for condensational growth.
      !-------------------------------------------------------------------------------------------------------------------
      real(8) :: dp_condtable(n_dp_condtable)                  ! [m] tabulated particle ambient diameters.
      real(8) :: xln_scale_dp                            ! [1] ln of table diameter ratio.
      real(8) :: kci_dp_condtable     (n_dp_condtable,nlays)  ! [m^3/s] tabulated condensational growth factors.
      real(8) :: kci_dp_condtable_aeq1(n_dp_condtable,nlays)  ! [m^3/s] tabulated condensational growth factors


	                                                          !   for the mass accommodation coefficient

      !From aero_coag
      !-------------------------------------------------------------------------------------------------------------------------------
      ! mode-average coagulation coefficients for mode-mode (i-j) interactions.
      !-------------------------------------------------------------------------------------------------------------------------------
      integer, parameter :: kij_ndgs  = kij_ndgs_set ! number of geo. mean diameters
      integer, parameter :: kij_nsgs  =      3       ! number of geo. std. deviations
      real(4), parameter :: kij_temp1 =    325.0     ! [k]  \
      real(4), parameter :: kij_temp2 =    260.0     ! [k]  - must have t1 > t2 > t3
      real(4), parameter :: kij_temp3 =    200.0     ! [k]  /
      real(4), parameter :: kij_pres1 = 101325.0     ! [pa] \
      real(4), parameter :: kij_pres2 =  10132.50    ! [pa] - must have p1 > p2 > p3
      real(4), parameter :: kij_pres3 =   1013.250   ! [pa] /
      real(4), parameter :: kij_sigm1 =      1.6     ! [1]
      real(4), parameter :: kij_sigm2 =      1.8     ! [1]
      real(4), parameter :: kij_sigm3 =      2.0     ! [1]
      !-------------------------------------------------------------------------------------------------------------------------------
      ! kij_dgmin must be smaller than dpmin_global / smax, where smax = exp[1.5(ln sigma_max)^2] where sigma_max is the largest
      ! lognormal geometric standard deviation occurring for any mode. currently set for sigma_max = 2.0.
      !-------------------------------------------------------------------------------------------------------------------------------
      real(4), parameter :: kij_dgmin = 1.0d+06 * dpmin_global / 2.1d+00           ! [um] if any mode has sigma>2.0, must modify.
      real(4), parameter :: kij_dgmax =    100.0000                                ! [um]
      real(4) :: k0ij_temp1pres1(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)	! [m^3/s]
      real(4) :: k0ij_temp1pres2(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)	! [m^3/s]
      real(4) :: k0ij_temp1pres3(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)	! [m^3/s]
      real(4) :: k0ij_temp2pres1(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)	! [m^3/s]
      real(4) :: k0ij_temp2pres2(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)	! [m^3/s]
      real(4) :: k0ij_temp2pres3(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)	! [m^3/s]
      real(4) :: k0ij_temp3pres1(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)	! [m^3/s]
      real(4) :: k0ij_temp3pres2(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)	! [m^3/s]
      real(4) :: k0ij_temp3pres3(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)	! [m^3/s]
      real(4) :: k3ij_temp1pres1(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)	! [m^3/s]
      real(4) :: k3ij_temp1pres2(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)	! [m^3/s]
      real(4) :: k3ij_temp1pres3(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)	! [m^3/s]
      real(4) :: k3ij_temp2pres1(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)	! [m^3/s]
      real(4) :: k3ij_temp2pres2(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)	! [m^3/s]
      real(4) :: k3ij_temp2pres3(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)	! [m^3/s]
      real(4) :: k3ij_temp3pres1(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)	! [m^3/s]
      real(4) :: k3ij_temp3pres2(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)	! [m^3/s]
      real(4) :: k3ij_temp3pres3(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)	! [m^3/s]
      !-----------------------------------------------------------------------------------------------------------------
      ! kij_diameters contains the values of dg used in building the lookup tables.
      ! kij_sigmas    contains the values of sigmag used in building the lookup tables.
      ! index_sigg(i) is the index of kij_sigmas to obtain the sigmag value for mode i.
      !-----------------------------------------------------------------------------------------------------------------
      real(4)                      :: kij_diameters(kij_ndgs)                             ! [um]
      real(4), dimension(kij_nsgs) :: kij_sigmas = (/ kij_sigm1, kij_sigm2, kij_sigm3 /)  ! [1]
      integer, allocatable         :: index_sigg(:)                                ! [1]
      !-------------------------------------------------------------------------------------------------------------------------------
      ! constant coagulation coefficients for each mode-mode (i-j) interaction,
      ! based upon characteristic sizes for each mode.
      ! these are not mode-averaged and are no longer in use.
      !-------------------------------------------------------------------------------------------------------------------------------
      real(4), allocatable :: kij(:,:)  ! [m^3/s]
      !From tramp_diam
      real(8),allocatable :: diam_histogram(:,:,:)    ! [1]
      real*8, allocatable, dimension(:,:,:,:)     :: diam       ![m](i,j,l,nmodes)
      !From Tramp_init
      real(8), allocatable :: aero_in(:)
      real(8), allocatable :: gas_in(:)

CONTAINS

   !> @brief To allocate the matrix with cache and memory blocks
   !! @author Luiz Flavio
   !! @date Oct/2011
   subroutine allocateMatrix(ia,iz,ja,jz,m1,m2,m3,nzg)
      Implicit none

      INTEGER, INTENT(IN) :: ia !< i initial point
      INTEGER, INTENT(IN) :: iz !< i final point
      INTEGER, INTENT(IN) :: ja !< j initial point
      INTEGER, INTENT(IN) :: jz !< j final point
      INTEGER, INTENT(IN) :: m1 !< z size
      INTEGER, INTENT(IN) :: m2 !< i size
      INTEGER, INTENT(IN) :: m3 !< j size
      INTEGER, INTENT(IN) :: nzg

      INTEGER :: noc,istep,i1,i2

      NumberOfColunms=(iz-ia+1)*(jz-ja+1)        !Total of columns in this processor athmosphere
      ALLOCATE(matrixVar(NumberOfColunms))
      DO noc=1,numberOfColunms
         ALLOCATE(matrixVar(noc)%aero     (m1,naerobox))
         ALLOCATE(matrixVar(noc)%gas      (m1,ngases))
         ALLOCATE(matrixVar(noc)%emis_mas (m1,nemis_spcs))
         ALLOCATE(matrixVar(noc)%tk       (m1))
         ALLOCATE(matrixVar(noc)%rh       (m1))
         ALLOCATE(matrixVar(noc)%pres     (m1))
         ALLOCATE(matrixVar(noc)%aqso4rate(m1))
         ALLOCATE(matrixVar(noc)%wupdraft (m1))
         ALLOCATE(matrixVar(noc)%diag     (m1,ndiag_aero,naerobox))
         ALLOCATE(matrixVar(noc)%ac       (m1,nmodes))
         ALLOCATE(matrixVar(noc)%dens_mode_dry (m1,nmodes))
         ALLOCATE(matrixVar(noc)%nAct     ( m1,nmodes))
         ALLOCATE(matrixVar(noc)%vddep_aero(nmodes,2))
         ALLOCATE(matrixVar(noc)%diam( m1,nmodes))
         ALLOCATE(matrixVar(noc)%dn0           (m1))
         ALLOCATE(matrixVar(noc)%pi0           (m1))
         ALLOCATE(matrixVar(noc)%pp            (m1))
         ALLOCATE(matrixVar(noc)%ustar         (1))
         ALLOCATE(matrixVar(noc)%tstar         (1))
         ALLOCATE(matrixVar(noc)%aer2mp_eff(m1,4))
      END DO
      ALLOCATE(iPos(numberOfColunms))
      ALLOCATE(jPos(numberOfColunms))
      ALLOCATE(nColumn((m2),(m3)))
!srf- for RK2 and RK3
      ALLOCATE(matrixVarRK(2))
      DO istep=1,2
       ALLOCATE(matrixVarRK(istep)%aero(m1,naerobox))
       ALLOCATE(matrixVarRK(istep)%gas (m1,ngases)  )
      ENDDO
      !srf- for RK and RK3

   END subroutine allocateMatrix

   subroutine allocateMatrixSetup()
      allocate(numb_map(nweights))
      allocate(mass_map(nweights,nmass_spcs))
      allocate(prod_index(nweights,nmass_spcs))
      allocate(nm_spc_name(nweights,nmass_spcs))
      allocate(giklq(nweights,nweights,nweights,nmass_spcs))
      allocate(dikl (nweights,nweights,nweights))
      allocate(dij  (nweights,nweights))
      allocate( xdij (nweights,nweights))
      allocate(recip_part_mass(nemis_spcs))
      allocate(kci_coef_dp     (nweights,nlays))
      allocate(kci_coef_dp_aeq1(nweights,nlays))
      allocate(theta_poly (nweights))
      allocate(dp0(nweights))
      allocate(dp0_emis(nweights))
      allocate( diffcoef_m2s(nlays))
      allocate(kappai(nweights))
!KML
      allocate(kappapk(nweights))

      allocate(denspi(nweights))
      allocate(dens_comp(nweights))
      allocate(recip_dens_comp(nweights))
      allocate(dgn0(nweights))
      allocate(sig0(nweights))
      allocate(lnsig0(nweights))
      allocate(dgn0_emis(nweights))
      allocate(sig0_emis(nweights))
      allocate(conv_dpam_to_dgn(nweights))
      allocate(kij(nweights,nweights))
      allocate(index_sigg(nweights))
      allocate(diam_histogram(nmodes,n_dp_condtable,2))
      allocate(aero_in(naerobox))
      allocate(gas_in(ngases))
      allocate(nm(nmodes))
      gas_in = 0.0d-30
   end subroutine allocateMatrixSetup


   !> @brief To deallocate the matrix to release memory for another time
   !! @author Luiz Flavio
   !! @date Oct/2011
   subroutine deAllocateMatrix
      IMPLICIT NONE

      INTEGER :: noc,istep
      deALLOCATE(iPos)
      deALLOCATE(jPos)
      deALLOCATE(nColumn)
      DO noc=1,NumberOfColunms
         DEALLOCATE(matrixVar(noc)%aero)
         deALLOCATE(matrixVar(noc)%gas)
         deALLOCATE(matrixVar(noc)%emis_mas)
         deALLOCATE(matrixVar(noc)%tk)
         deALLOCATE(matrixVar(noc)%rh)
         deALLOCATE(matrixVar(noc)%pres)
         deALLOCATE(matrixVar(noc)%aqso4rate)
         deALLOCATE(matrixVar(noc)%wupdraft)
         deALLOCATE(matrixVar(noc)%diag)
         DEALLOCATE(matrixVar(noc)%ac)
         DEALLOCATE(matrixVar(noc)%dens_mode_dry)
         DEALLOCATE(matrixVar(noc)%nAct     )
         DEALLOCATE(matrixVar(noc)%vddep_aero)
         DEALLOCATE(matrixVar(noc)%diam)
         DEALLOCATE(matrixVar(noc)%dn0)
         DEALLOCATE(matrixVar(noc)%pi0)
         DEALLOCATE(matrixVar(noc)%pp)
         DEALLOCATE(matrixVar(noc)%ustar)
         DEALLOCATE(matrixVar(noc)%tstar)
         DEALLOCATE(matrixVar(noc)%aer2mp_eff)
      END DO
      DEALLOCATE(matrixVar)
!srf- for RK2 and RK3
      DO istep=1,2
       DEALLOCATE(matrixVarRK(istep)%aero)
       DEALLOCATE(matrixVarRK(istep)%gas )
      ENDDO
      DEALLOCATE(matrixVarRK)
!srf- for RK and RK3
   END subroutine deAllocateMatrix

   !LFR  subroutine StoreNamelistFileAtMatrix(oneNamelistFile)
      !LFR  type(namelistFile), pointer :: oneNamelistFile

      !LFR  aerosol=oneNamelistFile%aerosol
      !LFR  aerfrq=oneNamelistFile%aerfrq
      !LFR  mech=oneNamelistFile%mech

   !LFR  END subroutine StoreNamelistFileAtMatrix


END MODULE memMatrix
