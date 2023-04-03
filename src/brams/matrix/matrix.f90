!!>
!!-----------------------------------------------------------------------------------------------------------------------
!!
!!@brief     This is the top-level routine of the MATRIX aerosol microphysical module.
!!@author    Susanne Bauer/Doug Wright
!!
!!     Here are some items to be familiar with in setting up the call to MATRIX.
!!
!!     1. The GAS array: 
!!
!!        GAS(1) contains the current H2SO4(g) concentration expressed as
!!        ug SO4=/m^3, including the H2SO4(g) produced by gas-phase chemistry
!!        during the current time step.
!!     
!!        GAS(2) contains the current HNO3(g) concentration expressed as
!!        ug NO3-/m^3, including the HNO3(g) produced by gas-phase chemistry
!!        during the current time step.
!!     
!!        GAS(3) contains the current NH3(g) concentration expressed as
!!        ug NH4+/m^3. Any sources of NH3(g) during the time step should  
!!        probably be added in before the call to MATRIX.
!!
!!     2. H2SO4 is always represented as SO4= (MW=96g/mol) rather than 
!!        H2SO4 (MW=98g/mol) (except in the calculation of the mean thermal
!!        velocity of an H2SO4 molecule). Likewise, HNO3 is always 
!!        represented as NO3- and NH3 as NH4+.
!!
!!     3. In this model, sea salt contains NaCl and potentially non-sea salt 
!!        sulfate, nitrate, ammonium, and water, depending upon conditions.
!!        The model sea salt species MASS_SSA_SEAS (accumulation mode),
!!        MASS_SSC_SEAS (coarse mode), and MASS_MXX_SEAS (sea salt in 
!!        a mixed mode) do not include the non-sea salt sulfate, nitrate,
!!        ammonium, and water, are treated as pure NaCl, and the density
!!        and MW of NaCl are used with these species.
!!
!!     4. EMIS_MASS(J) contains the mass emissions rates for species J. 
!!        The species and units are:
!!
!!          J               SPECIES                      UNITS
!!        -----   -----------------------------    ------------------
!!          1       Aitken-mode sulfate              [ug SO4 /m^3/s]
!!          2       accumulation-mode sulfate        [ug SO4 /m^3/s]
!!          3       black carbon                     [ug BC  /m^3/s]
!!          4       organic carbon                   [ug OC  /m^3/s]
!!          5       mineral dust - size 1            [ug DUST/m^3/s]
!!          6       accumulation-mode sea salt       [ug NaCl/m^3/s]
!!          7       coarse-mode sea salt             [ug NaCl/m^3/s]
!!          8       BC in the mixed BC-OC mode       [ug BC  /m^3/s]
!!          9       OC in the mixed BC-OC mode       [ug OC  /m^3/s]
!!         10       mineral dust - size 2            [ug DUST/m^3/s]
!!
!!-----------------------------------------------------------------------------------------------------------------------
SUBROUTINE Matrix(var,tstep,m1,noc,mech,testCase,do_Setup_Kci_dynj)
!LFR      USE Aero_config, only:  
      ! Module AERO_SETUP uses module AERO_PARAM, so the statement USE AERO_PARAM is not needed here.
!srf
      USE memMatrix      , ONLY: nmodes
 
      USE Aero_setup, ONLY: &
                            Setup_Config,              & !subroutine
                            Setup_Species_Maps,        & !subroutine
                            Setup_Aero_Mass_Map,       & !subroutine
                            Setup_Coag_Tensors,        & !subroutine
                            Setup_Dp0,                 & !subroutine
                            Setup_Emis,                & !subroutine
                            Setup_Kci_dyn                !subroutine


      USE Aero_subs
                         
      USE Aero_coag,   ONLY: &
                            setup_kij_diameters,       & !subroutine
                            setup_kij_tables,          & !subroutine
                            get_kbarnij                  !subroutine
      USE Aero_npf,    ONLY: &
                            dnu,                       & !intent()
                            npfrate,                   & !subroutine
                            setup_npfmass,             & !subroutine
                            steady_state_h2so4           !subroutine   
      USE Aero_actv,   ONLY: &
                            getactfrac                   !subroutine
      USE Aero_depv,   ONLY: &
                            get_aero_depv               !subroutine
      USE memMatrix,   only: &
                            matrixType, &
			    naerovars, &
			    nextra, &
!srf			    nmodes, &
			    nweights, &
			    naerobox, &
                            dgn0,                      & !intent()
                            dp0,                       & !intent()
                            intermodal_transfer,       & !intent()
                            lnsig0,                    & !intent()
                            denspi,                    & !intent()
                            numb_map,                  & !intent()
                            seas_mode_map,             & !intent()
                            seas_mode_mass_map,        & !intent()
                            recip_seas_mpp,            & !intent()
                            number_of_seasalt_modes,   & !intent()
                            sulf_map,                  & !intent()
                            seas_mode_map,             & !intent()
                            dust_map,                  & !intent()
                            seas_map,                  & !intent()
                            nm,                        & !intent()
                            prod_index,                & !intent()
                            mass_map,                  & !intent()
                            nmodes_seas,               & !intent()
                            mode_numb_seas,            & !intent()
                            recip_dens_comp,           & !intent()
                            dp_condtable_min,          & !intent()
                            xln_scale_dp,              & !intent()
                            n_dp_condtable,            & !intent()
                            theta_poly,                & !intent()
                            kci_dp_condtable,          & !intent()
                            kci_dp_condtable_aeq1,     & !intent()
                            avg_dp_of_avg_mass_meters, & !intent()
                            conv_dpam_to_dgn,          & !intent()
                            emis_mode_map,             & !intent()
                            recip_part_mass,           & !intent()
                            kci_coef_dp,               & !intent()
                            kci_coef_dp_aeq1,          & !intent()               
                            ndikl,                     & !intent()
                            dikl_control,              & !intent()
                            xdij,                      & !intent()
                            emis_spcs_map,             & !intent()
                            giklq_control,             & !intent()
                            dij,                       & !intent()
                            kappai,                    & !intent()
			    kappapk_SULF,kappapk_BCAR,kappapk_OCAR,kappapk_DUST,kappapk_SEAS, & !KML
			    emis_dens_sulf,emis_dens_bcar,emis_dens_ocar,emis_dens_dust,emis_dens_seas,  & !KML
			    mode_name,                 & !KML
                            sig0,                      & !intent()
                            diam_histogram,            & !intent()
                            aero_spcs      
!      USE Aero_Isrpia, only : InitIsrpia
      IMPLICIT NONE

! Arguments.
 
      INTEGER, INTENT(IN)    :: m1            !< Size of each column
      INTEGER, INTENT(IN)    :: noc           !< Number of column
      INTEGER, INTENT(IN)    :: mech          !< Mechanism
      real(8), intent(in)    :: tstep         !< model physics time step [s]
      TYPE(matrixType),INTENT(INOUT) :: var   !< MatrixType (see: memMatrix)
      LOGICAL, INTENT(IN) :: testCase         !< Used for test of patterns
      LOGICAL, INTENT(IN) :: do_Setup_Kci_dynj!< Used for call/update setup_kci,dynj routine

!Local Parameters
      real(8), parameter :: kcmin = 1.0d-08 ! [1/s] minimum condensational sink - see notes of 10-18-06
      real(8), parameter :: dpakk0          = dnu     ! min. diameter of average mass for mode akk [m]
      real(8), parameter :: fnum_max        = 0.5d+00 ! max. value of fnum [1]
      real(8), parameter :: akk_minnum_imtr = 1.0d+06 ! min. akk number conc. to enable imtr [#/m^3]
      real(8), parameter :: dpcut_imtr  = sqrt( dg_akk * dg_acc )   ! fixed diameter for intermodal transfer [um] 
                                                                    !this is the formula of easter et al. 2004. (= 0.053 um)
      real(8), parameter :: lnsg_akk    = log( sg_akk )
      real(8), parameter :: lnsg_acc    = log( sg_acc )
      real(8), parameter :: xnum_factor = 1.0d+00 / ( sqrt( 2.0d+00 ) * lnsg_akk ) ! factor in xnum expression [1]
      real(8), parameter :: x3_term     = 3.0d+00 * lnsg_akk / sqrt( 2.0d+00 ) ! term   in x3   expression [1]
      integer, parameter :: imtr_exp = 4              ! exponent for akk --> acc intermodal transfer
      real(8), parameter :: piq_thresh = 1.0d-08  ! [1] threshold in number/mass conc. solver
      real(8), parameter :: xh2so4_nucl_min_ncm3 = 1.00d+03 ! min. [h2so4] to enter nucleation calculations [#/cm^3]  
      real(8), parameter :: xh2so4_nucl_min = xh2so4_nucl_min_ncm3 * mw_so4 * 1.0d+12 / avo  ! convert to [ugso4/m^3] 
      real(8), parameter :: xntau = 2.0d+00                 ! number of time constants in the current time step  
      ! These variables are used as default values in computing dry deposition velocities.
      real(8), parameter :: temp_ddep = 288.15d+00      ! temperature [ k ]
      real(8), parameter :: rhoa_ddep =  1.225d+00      ! air density [ kg/m^3 ]
      real(8), parameter :: lamb_ddep = 6.6332d-08      ! mean free path of air [ m ]
      real(8), parameter :: dvis_ddep = 1.7894d-05      ! dynamic viscosity of air [ kg/(m s) ]
      real(8), parameter :: wstr_ddep =    1.0d+00      ! convective velocity scale [ m/s ]
      real(8), parameter :: ustr_ddep =    0.5d+00      ! friction velocity [ m/s ]
      real(8), parameter :: raer_ddep =    5.0d+00      ! aerodynamic resistance [ s/m ]
      logical, parameter :: get_dep_vel_only = .false.  ! for early return after getting dep. velocities
      real(8), parameter :: n_min_diam_histogram = 1.0d+04  ! min. # conc. for count in diam_histogram [#/m^3] 
      INTEGER, PARAMETER :: first_layer=2
      REAL,PARAMETER :: cpi = 1. / 1004.,ep13=1./3.,g = 9.80
      real(8), parameter :: erfconst = -4.0d+00 / pi

! Local variables.
      integer :: i,j,k,l,q,qq              ! indices
      integer :: index_dp, index_dp_dry    ! index for condensation factor lookup table
      integer :: ibranch                   ! scratch debugging variable [1]
      real(8) :: bi(nweights)              ! number conc. coefficients [1/s]
      real(8) :: ci(nweights)              ! number conc. coefficients [#/m^3/s]
      real(8) :: ri(nweights)              ! number conc. coefficients [#/m^3/s]
      real(8) :: ni(nweights)              ! number concentrations [#/m^3]
      real(8) :: mjq(nweights,nmass_spcs)  ! mjq(j,q) is the avg. mass/particle of species q (=1-5) for mode j
      real(8) :: piq(nweights,nmass_spcs)  ! production terms for mass conc. [ug/m^3/s]
      real(8) :: fi(nweights)              ! loss coefficients for mass conc. [1/s]
      real(8) :: tot_sulf                  ! total sulfate conc.  [ug/m^3]
      real(8) :: tot_dust                  ! total dust conc.     [ug/m^3]
      real(8) :: tot_seas                  ! total sea salt conc. [ug/m^3]
      real(8) :: ssh2o                     ! total sea salt h2o   [ug/m^3]
      real(8) :: ssh2o_per_ssmass          ! total sea salt h2o / total sea-salt dry mass
      real(8) :: ftmp, voltmp, voltmp_dry  ! scratch variables for computing mode mean diameters 
      real(8) :: tot_mass    (nmodes)      ! total ambient mass conc. for each mode [ug/m^3]
      real(8) :: tot_mass_dry(nmodes)      ! total dry     mass conc. for each mode [ug/m^3]
      real(8) :: h2o_mass(nmodes)          ! water mass mass conc. for each mode [ug/m^3]
      real(8) :: mass_comp(nmodes,8)       ! mass conc. of each component for each mode [ug/m^3]
      real(8) :: dens_mode    (nmodes)     ! average mode density calculated from component concentrations [g/cm^3]
      real(8) :: optot_no3nh4h2o_to_sulf   ! 1 + ( total no3+nh4+h2o / total so4 )
      real(8) :: aero_water_actual         ! actual aerosol h2o conc. [ug/m^3]
      real(8) :: aero_water_wet            ! wet    aerosol h2o conc. [ug/m^3]
      real(8) :: rhd                       ! deliquescence   rh [0-1]
      real(8) :: rhc                       ! crystallization rh [0-1]
      real(8) :: dgn(nweights)             ! ambient geometric mean diameter of the number distribution for each mode [m]
      real(8) :: dgn_dry(nweights)         ! geometric mean dry diameter of the number distribution for each mode [m]
      real(8) :: dp(nweights)              ! ambient diameter of average mass of the number distribution for each mode [m]
      real(8) :: dp_dry(nweights)          ! dry diameter of average mass of the number distribution for each mode [m]
      real(8) :: p_emis_numb(nmodes)       ! number emission rates [#/m^3/s]       
      real(8) :: spcmass1(nmass_spcs+2)    ! initial total mass conc. of each model species [ug/m^3]
      real(8) :: spcmass2(nmass_spcs+2)    ! final   total mass conc. of each model species [ug/m^3]
      real(8) :: fseassulf                 ! fraction of sulfate assigned to sea salt in mode ssc (coarse mode)
      real(8) :: factor                    ! scratch variable in mass conc. equation solver
      real(8) :: y0,y,a,b,c,delta,r1,r2    ! scratch variables in the number concentration analytic solution
      real(8) :: gamma,gexpdt,expdt        ! scratch variables in the number concentration analytic solution
!KML
      real(8) :: kappapk(nmodes)
      real(8) :: kappapk_eff(4),vol_total(4),ni_total(4),diam_eff(4)
      real(8) :: mi5_sum1(nmodes)
      real(8) :: mi5_sum2   
      ! For the condensational sink, condensational growth, and cou. 

      !----------------------------------------------------------------------------------------------------------------
      ! This minimum value for the condensational sink, 1.0D-08 [1/s], is calculated as 
      ! 
      !   KC = 2 * PI * D * Beta * Dpbar * N = 2 * 3.14 * (1.0D-05 m^2/s) * (1) * (1.0D-07 m) * (1.0D+04 #/m^3)
      !      = 6.28D-08 1/s --> 1.0D-08 1/s  
      !----------------------------------------------------------------------------------------------------------------
      real(8) :: kc                         ! total condensational sink with arbitrary mass accommodation coef. [1/s]
      real(8) :: kc_aeq1                    ! total condensational sink with unity     mass accommodation coef. [1/s]
      real(8) :: pq_growth                  ! avg. cond. rate of h2so4 (as so4) over the time step [ugso4/m^3/s]
      real(8) :: tot_h2so4_loss             ! net loss of h2so4 (as so4) [ugso4/m^3]

      ! Mode-average coagulation coefficients. 
      real(8) :: kbar0ij(nweights,nweights)     ! mode-average coagulation coefficients [m^3/s] for number
      real(8) :: kbar3ij(nweights,nweights)     ! mode-average coagulation coefficients [m^3/s] for mass 

      ! These variables are used in intermodal transfer.
      real(8)            :: dpakk                     ! diameter of average mass for mode akk [m]
      real(8)            :: dpacc                     ! diameter of average mass for mode acc [m] 
      real(8)            :: fnum                      ! fraction of akk number transferred over the time step [1]
      real(8)            :: f3                        ! fraction of akk mass   transferred over the time step [1]
      real(8)            :: xnum, x3                  ! error function complement arguments [1]
      real(8)            :: dgn_akk_imtr              ! geo. mean diam. of the akk mode number distribution [m]
      real(8)            :: dgn_acc_imtr              ! geo. mean diam. of the acc mode number distribution [m]
      real(8)            :: del_numb                  ! number conc. transferred from akk to acc [ #/m^3]
      real(8)            :: del_mass                  ! mass   conc. transferred from akk to acc [ug/m^3]

      ! These variables and parameters are used in nucleation and new particle formation.
      integer :: klq
      integer :: ikl
      integer :: ih2so4_path               ! index for branch in the [h2so4] calc. for nucl and gr [1]
      real(8) :: xh2so4_ss                 ! steady-state h2so4 (as so4) conc. [ugso4/m^3]
      real(8) :: xh2so4_ss_wnpf            ! steady-state h2so4 (as so4) conc. including npf [ugso4/m^3]
      real(8) :: xh2so4_nucl               ! h2so4 (as so4) conc. used in nucleation and gr calculation [ugso4/m^3]
      real(8) :: xh2so4_init               ! h2so4 at the top of the time step [ugso4/m^3]
      real(8) :: so4rate                   ! h2so4(g) production rate in/as [ugso4/m^3/s]
      real(8) :: xnh3                      ! ammonia concentration [ppmv]
      real(8) :: dndt                      ! number production rate [1/m^3/s]     from npf
      real(8) :: dmdt_so4                  ! so4 mass product. rate [ugso4/m^3/s] from npf
                                                            ! required to invoke the steady-state approx. [1] 
      ! These variables are used in aerosol activation.
      real(8) :: rsum_activ        ! used in pcloud_i,q terms
      real(8) :: nsol    (nmodes)  ! number concentration of soluble particles for each mode [#/m^3] 
      real(8) :: fracactn(nmodes)  ! activated fraction of number concentration for each mode [1]
      real(8) :: fracactm(nmodes)  ! activated fraction of mass   concentration for each mode [1]
      real(8) :: mact    (nmodes)  ! activated mass   concentration for each mode [ug/m^3] 
      real(8) :: mi5     (nmodes,nmass_spcs)  ! mass conc. of each species for each mode [ug/m^3] 
      real(8) :: pis,temps,dens,amu,v_air,xlm,ustar,rhoa,wptp,wstar,tstar

      ! These variables are used in computing the DIAG or DIAM_HISTOGRAM diagnostic arrays. 
      real(8) :: diagtmp1(ndiag_aero,naerobox)          ! scratch work array [ug/m^3/s] or [#/m^3/s]
      real(8) :: aerotmp1(naerobox)                     ! input   aerosol concentrations [ug/m^3] or [#/m^3]
      real(8) :: aerotmp2(naerobox)                     ! scratch aerosol concentrations [ug/m^3] or [#/m^3]
      real(8) :: piqtmp(nweights,nmass_spcs)            ! scratch work array for mass production terms [ug/m^3/s]
     
      real(8) :: vi5(nmodes)                            ! vol [m^3] 
      real(8) :: number_in                              ! 
      !----------------------------------------------------------------------------------------------------------------
      ! Error function statement function derived from code in CMAQ v4.4. 
      ! 
      ! The CMAQ source is Meng, Z., and J.H.Seinfeld (1994) On the source of the submicrometer droplet mode of 
      ! urban and regional aerosols. Aerosol Sci. and Technology, 20:253-265, who cite
      ! Reasearch & Education Association (REA), (1991) Handbook of Mathematical, Scientific, 
      ! and Engineering Formulas, Tables, Functions, Graphs, Transforms: REA, Piscataway, NJ. p. 493.
      !----------------------------------------------------------------------------------------------------------------
      character(LEN=50) :: ca,cb,cc
      integer :: col
      real(8)            :: xx, erf, erfc                          ! error function, error function complement
      erf(xx) = sqrt(1.0d+00 - exp( erfconst * xx * xx ) )
      erfc(xx) = 1.0d+00 - erf(xx)

      !-Variables to prepare to dry-deposition velocity
      pis = (var%pp(1) + var%pp(first_layer) + var%pi0(1) + var%pi0(first_layer) ) * .5 * cpi
      temps= var%tk(first_layer)
      dens= (var%dn0(1) + var%dn0(first_layer)) * .5
      amu= 1.8325e-5*(416.16/(temps+120.))*(temps/296.16) * sqrt(temps/296.16)
      v_air = sqrt( 7.3102e+2 * temps)
      xlm = 2.* amu /(dens*v_air)
      ustar=var%ustar(1)
      tstar=var%tstar(1)
      rhoa=var%dn0(first_layer)
      wptp  = - ustar * tstar ! sensible heat flux
      wstar = ( max (0., g* var%Zi*  wptp/temps) )**ep13
      !print*,"1";call flush(6)
      
      
      !Getting the dry-deposition velocity for current column
      IF(testCase) THEN
         call Get_Aero_Depv(nmodes,temp_ddep,rhoa_ddep,lamb_ddep,dvis_ddep, &
                             wstr_ddep,ustr_ddep,raer_ddep,dgn0,lnsig0,denspi,var%vddep_aero)     
      ELSE
      !print*,nmodes,var%tk(first_layer),rhoa,xlm,amu, &
      !                        wstr_ddep,ustar,raer_ddep,dgn0,lnsig0,denspi,var%vddep_aero
!			      call flush(6)
         call Get_Aero_Depv(nmodes,var%tk(first_layer),rhoa,xlm,amu, &
                              wstr_ddep,ustar,raer_ddep,dgn0,lnsig0,denspi,var%vddep_aero)
      END IF
      
      !print*,"vdep=",var%vddep_aero
      !call flush(6)
      !stop 333
      !calculate the coefficients that multiply the number
      !concentrations, or the number concentrations times the particle diameters,
      !to obtain the condensational sink for each mode or quadrature point.
      !srf
      if(do_Setup_Kci_dynj) CALL Setup_Kci_dyn(var%pres(:),var%tk(:),m1) 
!----------------------------------------------------------------------------------------------------------------------
!     Begin execution.
!----------------------------------------------------------------------------------------------------------------------
      !print*,"1";call flush(6)
      
      DO col=1,m1 !Begin loop in column

         ilay=col
         ! these variables are used in intermodal transfer.        
         if( intermodal_transfer .and. numb_akk_1 .eq. 0 ) then
            write(*,*)'intermodal_transfer must be set to .false. for this mechanism.'
            stop
         endif
         IF(write_log) THEN
            if( intermodal_transfer ) then
               write(aunit1,'(/a/)') 'intermodal transfer (akk->acc) is turned on.'
               write(aunit1,'(a,5f12.6)')'dpcut_imtr, xnum_factor, x3_term, lnsg_akk, lnsg_acc = ', &
                                                      dpcut_imtr, xnum_factor, x3_term, lnsg_akk, lnsg_acc
            else
               write(aunit1,'(/a/)') 'intermodal transfer (akk->acc) is turned off.'
            endif
         END IF
         if ( .not. no_microphysics ) then
            if( update_Kij .eq. 0 ) then
               dgn(:)     = dgn0(:)
               dgn_dry(:) = dgn0(:)
               call Get_kBarnIj(0,var%tk(col),var%pres(col),dgn,kbar0ij,kbar3ij)
            endif
         endif      
   
         !----------------------------------------------------------------------------------------------------------------
         ! the mass-only, no-microphysics option.
         ! add in emissions, h2so4(g), so4(aq), and do gas-particle partitioning only.
         !---------------------------------------------------------------------------------------------------------------- 
         if( no_microphysics ) then
            var%emis_mas(col,:) = 0.0d+00 
            !!!LFR ????? call Aero_NoMicrophysics(aero,gas,emis_mass,tstep,var%tk(col),var%rh(col),var%pres(col),var%aqso4rate(col))
            PRINT *,'LFR -> No microphysics...where this routine? I did find it!'
            STOP
            !return
         endif
   
         !----------------------------------------------------------------------------------------------------------------
         ! zero entire diagnostics array. 
         !----------------------------------------------------------------------------------------------------------------
         diagtmp1(:,:) = 0.0d+00 
         aerotmp1(:)   = var%aero(col,:)   ! input concentrations
         
         !----------------------------------------------------------------------------------------------------------------
         ! get initial total mass conc. for each model species [ug/m^3].
         ! these are (optionally) used to enforce precise mass conservation.
         !----------------------------------------------------------------------------------------------------------------
         if( mass_adj ) call SpcMasses(var%aero(col,:),var%gas(col,:),spcmass1)
   
         !----------------------------------------------------------------------------------------------------------------
         ! load the number concentrations [#/m^3] into the local work array.
         !----------------------------------------------------------------------------------------------------------------
         ni(:) = var%aero(col, numb_map(:) ) + tinydenom
       
         !----------------------------------------------------------------------------------------------------------------
         ! for the sea salt modes, use characteristic mean particle masses
         ! to derive the number concentration from the mass concentration.
         !
         ! there may be two sea salt modes, ssa and ssc, or just one, sss.
         !
         ! seas_mode_map(i)      is the mode number of sea salt mode i.
         ! seas_mode_mass_map(i) is the location in aero of the sea salt (nacl)
         !                         mass concentration for sea salt mode i.
         ! recip_seas_mpp(i)     is the reciprocal of the mean nacl mass 
         !                         per sea salt particle for sea salt mode i.
         !----------------------------------------------------------------------------------------------------------------
         ni(seas_mode_map(:)) = var%aero(col, seas_mode_mass_map(:) ) * recip_seas_mpp(:) + tinydenom !ACrescentei tinydenom

!DO i=1,nmodes
!   write (*,FMT='(A,I1,A,E20.8)') 'LFR=DBG ni(',i,')= ',ni(i)
!END DO
  
         !----------------------------------------------------------------------------------------------------------------
         ! partition the sea salt sulfate between the ssa and ssc modes, if both
         ! of these modes are present. in this case, all sea salt sulfate is
         ! passed to this routine in mass species mass_ssa_sulf.
         !----------------------------------------------------------------------------------------------------------------
         ! write(*,*)'number_of_seasalt_modes = ', number_of_seasalt_modes
         if( number_of_seasalt_modes .eq. 2 ) then
         fseassulf = var%aero(col, mass_ssc_seas ) / ( var%aero(col, mass_ssc_seas ) + var%aero(col, mass_ssa_seas ) + tinydenom )
         var%aero(col, mass_ssc_sulf ) = var%aero(col, mass_ssa_sulf ) * fseassulf
         var%aero(col, mass_ssa_sulf ) = var%aero(col, mass_ssa_sulf ) * ( 1.0d+00 - fseassulf )
         endif
   
         !----------------------------------------------------------------------------------------------------------------
         ! get the total sulfate, total mineral dust, and total sea salt (nacl)
         ! mass concentrations summed over all quadrature points/modes. [ug/m^3]
         !----------------------------------------------------------------------------------------------------------------
         tot_sulf = sum( var%aero(col,sulf_map(:)) ) + tinynumer   ! insure nonzero value
         tot_dust = sum( var%aero(col,dust_map(:)) )
         tot_seas = sum( var%aero(col,seas_map(:)) ) + tinynumer   ! insure nonzero value
   
         if( write_log ) write(aunit1,'(/a,f12.4/)') 'tot_sulf = ', tot_sulf
   
   
         !----------------------------------------------------------------------------------------------------------------
         ! call the aerosol thermodynamic module to determine the bulk gas-particle 
         ! partitioning of the inorganic species and the water content
         ! associated with those species. also determine the water content
         ! associated with the nacl of sea salt.
         !
         ! note that aero(mass_h2o) includes the sea-salt associated water 
         ! when it is passed to this routine and passed to aero_thermo.
         ! upon return from aero_thermo, aero(mass_h2o) contains only the 
         ! non-sea salt-associated water. the sea salt-associated water is in ssh2o.
         !----------------------------------------------------------------------------------------------------------------
         aero_water_actual = var%aero(col,mass_h2o)       ! actual tracked aerosol water conc.
         call Aero_Thermo(tot_sulf,var%aero(col,mass_no3),var%aero(col,mass_nh4),var%aero(col,mass_h2o),var%gas(col,gas_nh3), &
                        var%gas(col,gas_hno3),tot_dust,tot_seas,ssh2o,var%tk(col),var%rh(col),var%pres(col),rhd,rhc)
         aero_water_wet = var%aero(col,mass_h2o) + ssh2o  ! total metastable aerosol water conc.
   
         !----------------------------------------------------------------------------------------------------------------
         ! adjust the aerosol water concentration for hysteresis.
         !
         ! if the aerosol water concentration is less than half its metastable (wet)
         ! concentration, treat the aerosol as dry. otherwise, the values
         ! in aero(mass_h2o) and ssh2o remain at their wet values. since the water
         ! associated with sea salt is not tracked separately, this treatment of
         ! hysteresis is based on the total aerosol water including that of sea
         ! salt, although the rhd governing this is that obtained (or set) for the
         ! non-sea salt inorganics.
         !
         ! this is done only for rh between the crystallization and deliquescence
         ! rhs of the non-sea salt inorganics. 
         !
         ! currently, rhc = 80%, as for ammonium sulfate (ghan et al. 2001), and 
         !            rhd = 35%, as for ammonium sulfate (ghan et al. 2001), both set in subr. aero_thermo. 
         !----------------------------------------------------------------------------------------------------------------
         if ( var%rh(col) .gt. rhc  .and.  var%rh(col) .lt. rhd ) then
         if ( aero_water_wet .gt. 0.0d+00 ) then
            if ( aero_water_actual/aero_water_wet .lt. 0.5d+00 ) then  ! dry aerosol
               var%aero(col,mass_h2o) = 0.0d+00  ! zero the non-sea salt-associated water.
               ssh2o = 0.0d+00           ! zero the     sea salt-associated water.
            else                                                       ! wet (metastable) aerosol
               ! leave aero(mass_h2o) and ssh2o at their metastable (wet) values. 
            endif
         endif
         elseif ( var%rh(col) .le. rhc ) then     ! insure that the aerosol is dry below the rhc.
         var%aero(col,mass_h2o) = 0.0d+00
         ssh2o = 0.0d+00               ! since the rhc of nacl is set to 45%, and rhc
         endif                           ! is currently set to 35% (ammonium sulfate),
                                       ! ssh2o should already be zero here. 
   
   
         !----------------------------------------------------------------------------------------------------------------
         ! get the mass of species q in mode (or quadrature point) i for the
         !   principal aerosol-phase chemical species.
         !   these species are sulfate, bc, oc, dust, and sea salt.
         !
         ! mass_map(i,q) is the location in aero(:) of the qth mass in mode i.
         ! mode i has nm(i) mass species defined for it, and nm(i) varies between 1 and nmass_spcs (=5).
         ! the second index of prod_index has nmass_spcs (=5) values:
         !   1=sulfate, 2=bc, 3=oc, 4=dust, 5=sea salt.
         ! the second index of mjq also has these same nmass_spcs (=5) values.
         !----------------------------------------------------------------------------------------------------------------
         if( write_log ) write(aunit1,'(/a/)')'i,q,mjq(i,qq),aero(mass_map(i,q)),ni(i),mass_map(i,q)'
         do i=1, nweights                               ! loop over modes (quadrature points)
         if ( ni(i) .gt. 1.0d-10 ) then        
            mjq(i,:) = 1.0d-15                         ! [ug/particle] - all species in mode i
            do q=1, nm(i)                              ! loop over species defined for mode i
               qq = prod_index(i,q)                     ! qq is in the range 1 to 5.
               mjq(i,qq) = var%aero(col,mass_map(i,q)) / ni(i)  ! [ug/particle]
               mjq(i,qq) = max( mjq(i,qq), 1.0d-15)     ! [ug/particle] - all species in mode i
               ! write(30,'(2i4,3d15.5,i4)') i,q,mjq(i,qq),aero(mass_map(i,q)),ni(i),mass_map(i,q)
               if( write_log ) write(aunit1,'(2i4,3d15.5,i4)') i,q,mjq(i,qq),var%aero(col,mass_map(i,q)),ni(i),mass_map(i,q)
               !print*,i,q,mass_map(i,q)
	    enddo
         else
            mjq(i,:) = 1.0d-15                         ! [ug/particle] - all species in mode i
         endif
         enddo
        
         !----------------------------------------------------------------------------------------------------------------
         ! get the ambient diameter of average mass for each mode.
         !----------------------------------------------------------------------------------------------------------------
         if ( update_dp ) then
         tot_sulf = sum( var%aero(col,sulf_map(:)) ) + tinynumer   ! insure nonzero value
         mass_comp(:,:)  = 0.0d+00    ! [ug/m^3] mass of each component in each mode
         tot_mass    (:) = 0.0d+00    ! [ug/m^3] total mass in each mode
         tot_mass_dry(:) = 0.0d+00    ! [ug/m^3] total mass in each mode
         do i=1, nweights
            do q=1, nm(i)
               qq = prod_index(i,q)    ! qq is in the range 1 to 5.
               !----------------------------------------------------------------------------------------------------------
               ! sulfate, bc, oc, dust, and sea salt concentrations.
               !----------------------------------------------------------------------------------------------------------
               mass_comp(i,qq) = var%aero(col,mass_map(i,q))   ! [ug/m^3]
            enddo
            !------------------------------------------------------------------------------------------------------------
            ! nitrate, ammonium, and (non-sea salt) water concentrations.
            !------------------------------------------------------------------------------------------------------------
            ftmp = var%aero(col,sulf_map(i)) / tot_sulf
            mass_comp(i,6) = ftmp * var%aero(col,1)           ! [ug/m^3] nitrate
            mass_comp(i,7) = ftmp * var%aero(col,2)           ! [ug/m^3] ammonium
            mass_comp(i,8) = ftmp * var%aero(col,3)           ! [ug/m^3] water
            ! write(*,'(i5,6d15.5)')i,ftmp,var%aero(col,1:3), tot_sulf, var%aero(col,sulf_map(i))
         enddo
         !--------------------------------------------------------------------------------------------------------------
         ! distribute the sea-salt associated water over modes in proportion to the sea salt mass present.
         !--------------------------------------------------------------------------------------------------------------
   
         ssh2o_per_ssmass = ssh2o / tot_seas                                  ! [ugh2o/ugnacl]
         do j=1, nmodes_seas    ! loop over all modes containing sea salt
            mass_comp(mode_numb_seas(j),8) = mass_comp(mode_numb_seas(j),8) &
                                          + ssh2o_per_ssmass * var%aero(col,seas_map(j)) ! [ug/m^3]
         enddo
         do i=1, nweights
            do j=1, 7                                                  ! all components except water 
               tot_mass_dry(i) = tot_mass_dry(i) + mass_comp(i,j)       ! [ug/m^3]
            enddo
            tot_mass(i) = tot_mass_dry(i) + mass_comp(i,8)             ! add in water [ug/m^3]
            !------------------------------------------------------------------------------------------------------------
            ! if(i.eq. 3) write(*,'(i4,a6,9f12.6 )') i, mode_name(i), tot_mass(i), mass_comp(i,:)
            ! if(i.eq. 8) write(*,'(i4,a6,9f12.6 )') i, mode_name(i), tot_mass(i), mass_comp(i,:)
            ! if(i.eq.16) write(*,'(i4,a6,9f12.6/)') i, mode_name(i), tot_mass(i), mass_comp(i,:)
            !------------------------------------------------------------------------------------------------------------
         enddo      
   
         do i=1, nweights
            !------------------------------------------------------------------------------------------------------------
            ! get the ambient and dry diameter of average mass for each mode (quadrature point).
            !------------------------------------------------------------------------------------------------------------
            voltmp_dry = 1.0d-30
            do j=1, 7                                                        ! all components except water 
               voltmp_dry = voltmp_dry + mass_comp(i,j) * recip_dens_comp(j)  ! mode dry volume conc. [10^6 cm^3/m^3]
            enddo
            voltmp = voltmp_dry + mass_comp(i,8) * recip_dens_comp(8)        ! mode ambient volume conc. [10^6 cm^3/m^3]
            dens_mode    (i) = tot_mass    (i) / voltmp                      ! mode ambient density [g/cm^3] 
            var%dens_mode_dry(col,i) = tot_mass_dry(i) / voltmp_dry                  ! mode dry     density [g/cm^3] 
            dp    (i) = ( conv_vol_to_dp_fac * voltmp     / ni(i) )**0.333333333333333  ! [m]
            dp_dry(i) = ( conv_vol_to_dp_fac * voltmp_dry / ni(i) )**0.333333333333333  ! [m]
            if (ni(i) .lt. 1.) then
               dp(i)     = dp0(i)
               dp_dry(i) = dp0(i)
            endif   
            dp    (i) = min( max( dp    (i), dpmin_global ), dpmax_global )
            dp_dry(i) = min( max( dp_dry(i), dpmin_global ), dpmax_global )
            var%diam(col,i) = dp(i)   ! [m] - store for use outside this routine.
            !------------------------------------------------------------------------------------------------------------
            ! update values of kci_coef_dp(i) for the current diameter of average mass for each mode.
            ! theta_poly(i) prevents excessive condensation due to treating the mode as monodisperse.
            !   see okuyama et al., studies in binary nucleation: the dibutylphthalate/dioctylphthalate system,
            !   j. chem. phys. 89, p. 6442, 1988. 
            !------------------------------------------------------------------------------------------------------------
            index_dp = nint( log( dp(i) / dp_condtable_min ) / xln_scale_dp ) + 1
            index_dp = max( min( index_dp, n_dp_condtable ), 1 )
            if( ni(i) .gt. n_min_diam_histogram ) then
               index_dp_dry = nint( log( dp_dry(i) / dp_condtable_min ) / xln_scale_dp ) + 1
               index_dp_dry = max( min( index_dp_dry, n_dp_condtable ), 1 )
               diam_histogram(i,index_dp,    1) = diam_histogram(i,index_dp,    1) + 1.0d+00 
               diam_histogram(i,index_dp_dry,2) = diam_histogram(i,index_dp_dry,2) + 1.0d+00 
            endif
            kci_coef_dp     (i,ilay) = theta_poly(i) * kci_dp_condtable     (index_dp,ilay)
            kci_coef_dp_aeq1(i,ilay) = theta_poly(i) * kci_dp_condtable_aeq1(index_dp,ilay)
         enddo
         if( write_log ) then
            write(aunit1,'(/a8,4a25/)') 'i','dp0(i) [um]','dp [um]','dens_mode [g/cm^3]','dens_mode_dry [g/cm^3]'
            do i=1, nweights
               write(aunit1,'(i8,4d25.6)') i, dp0(i)*1.0d+06, dp(i)*1.0d+06, dens_mode(i), var%dens_mode_dry(col,i)
            enddo
         endif
         avg_dp_of_avg_mass_meters = sum( dp(:)*ni(:) ) / sum( ni(:) )  ! [m] used in aero_npf, kk02 gamma expression 
         else
            do i=1,nweights
            var%diam(col,i) = dp(i)   ! [m] - store for use outside this routine.
            enddo
         endif
         
         !----------------------------------------------------------------------------------------------------------------
         ! update the ambient and dry geometric mean diameters. 
         !----------------------------------------------------------------------------------------------------------------
         dgn(:)     = 1.0d+06 * dp    (:) * conv_dpam_to_dgn(:)   ! [um]
         dgn_dry(:) = 1.0d+06 * dp_dry(:) * conv_dpam_to_dgn(:)   ! [um]
         !----------------------------------------------------------------------------------------------------------------
         ! do i=1, nweights
         !   write(*,'(i6,3f18.8)') i, 1.0d+06*dp_dry(i), conv_dpam_to_dgn(i), dgn_dry(i)
         ! enddo
         ! write(*,'(a)') 
         !----------------------------------------------------------------------------------------------------------------
   
         !----------------------------------------------------------------------------------------------------------------
         ! update the dry deposition velocities if desired. 
         !----------------------------------------------------------------------------------------------------------------
         !if( update_vdep .and. k==first_layer ) then 
         !   call Get_Aero_Depv(nmodes,temp_ddep,rhoa_ddep,lamb_ddep,dvis_ddep, &
         !                     wstr_ddep,ustr_ddep,raer_ddep,dgn,lnsig0,denspi,var%vddep_aero)
         !   if( get_dep_vel_only ) return
         !endif
   
   
         !----------------------------------------------------------------------------------------------------------------
         ! update the coagulation coefficients, if desired. [m^3/s]
         !----------------------------------------------------------------------------------------------------------------
         select case( update_kij )
         case ( 0 )
         case ( 1 )
         call Get_kBarnIj(1,var%tk(col),var%pres(col),dgn,kbar0ij,kbar3ij)
         !--------------------------------------------------------------------------------------------------------------
         ! kbar0ij(:,:) = 1.0d-14
         ! kbar3ij(:,:) = 1.0d-14
         ! kbar3ij(:,:) = kbar0ij(:,:)
         !--------------------------------------------------------------------------------------------------------------
         end select
   
   
         !----------------------------------------------------------------------------------------------------------------
         !
         ! get the production and loss terms for the number concentrations.
         !
         !----------------------------------------------------------------------------------------------------------------
         ! emissions terms. 
         !
         ! p_emis_numb is in [#/m^3/s].
         ! emis_mode_map(j) is the mode number receiving the emissions species j.
         ! emis_mass is the mass emissions rate [ug/m^3/s].
         ! recip_part_mass is the reciprocal of the mean mass
         !   of an emitted particle for this emissions species [1/ug].
         !----------------------------------------------------------------------------------------------------------------
         p_emis_numb(:) = 0.0d+00
         do j=1, nemis_spcs
         p_emis_numb(emis_mode_map(j)) = p_emis_numb(emis_mode_map(j)) &
                                       + var%emis_mas(col,j)*recip_part_mass(j)
         enddo 
   
         diagtmp1(1,numb_map(:)) = p_emis_numb(:)  
         ! mass emissions are already treated before
         var%emis_mas(col,:) = 0.0d+00 
   
         if( write_log ) then
         write(aunit1,'(/a4,a30/)')'i','p_emis_numb(i) [ #/m^3/s]'
         do i=1, nmodes
            ! write(aunit1,'(i4,2d30.6)') i, p_emis_numb(i)
         enddo
         endif
   
         !----------------------------------------------------------------------------------------------------------------
         ! get the secondary-particle formation rate dndt [#/m^3/s] and the 
         !   secondary-particle sulfate mass production rate dmdt_so4 [ugso4/m^3/s].
         !   the new particle formation rate dndt is smaller than the nucleation 
         !   rate j and is (for most parameterizations) derived from the nucleation rate. 
         !
         ! for the calculation of dndt and dmdt_so4, it is assumed that all of the
         !   current h2so4 concentration was produced by gas-phase oxidation of 
         !   so2 during the current time step. so4rate is the average production rate 
         !   of h2so4 (represented as so4, mw=96g/mol) over the time step based on this
         !   assumption. dmdt_so4 and dndt are limited by the available h2so4. 
         !
         ! when appropriate, the nucleation rate is calculated using a steady-state 
         !   h2so4 concentration. this is done to avoid the spuriously high nucleation 
         !   rates that would result if the current h2so4 concentration, which has 
         !   accumulated without loss over the time step, were used. when the  
         !   condensational sink (kc) is small enough such that steady-state will not
         !   be reached during the time step, an estimate of the h2so4 concentration
         !   at the mid-point of the time step is used in nucleation and condensation
         !   calculations. 
         !
         ! the condensational sink kc is the first-order rate constant for the
         !   loss of h2so4 due to condensation. this is summed over all modes.
         !   for the kerminen and kulmala (2002) parameterization for the conversion
         !   of the nucleation rate to the new particle formation rate, the 
         !   condensational sink kc_aeq1 should obtained with the mass accommodation 
         !   coefficient set of unity. 
         !----------------------------------------------------------------------------------------------------------------
         xh2so4_init = var%gas(col, gas_h2so4 ) + tinynumer             ! input h2so4 concentration [ugso4/m^3]
         xh2so4_nucl =  xh2so4_nucl_min !tinynumer ! xh2so4_nucl_min              ! for the case  xh2so4_init .lt. xh2so4_nucl_min
         kc = sum( kci_coef_dp(:,ilay)*ni(:) )                  ! total condensational sink [1/s] for any value of the
         kc = max( kc, kcmin )                                  !   mass accommodation coefficient [1/s]
         if( do_npf .and. xh2so4_init .gt. xh2so4_nucl_min ) then
         kc_aeq1 = sum( kci_coef_dp_aeq1(:,ilay)*ni(:) )      ! total condensational sink for the
         kc_aeq1 = max( kc_aeq1, kcmin )                      !   mass accommodation coefficient set to unity  [1/s]
         xnh3 = convnh3 * var%gas(col, gas_nh3 ) * var%tk(col) / var%pres(col)          ! nh3 concentration; from [ug/m^3] to [ppmv]
         so4rate = xh2so4_init / tstep                        ! average h2so4 production rate [ugso4/m^3/s]
         if(kc*tstep .ge. xntau ) then                        ! invoke steady-state assumption
            ih2so4_path = 1
            xh2so4_ss = min( so4rate/kc, xh2so4_init )         ! steady-state h2so4 [ugso4/m^3]
            call Steady_State_H2SO4(var%pres(col),var%rh(col),var%tk(col),xh2so4_ss,so4rate,xnh3,kc,tstep,xh2so4_ss_wnpf)
            xh2so4_nucl = xh2so4_ss_wnpf                       ! [h2so4] for nucl., gr, and cond. calculation [ugso4/m^3]
         else
            ih2so4_path = 2
            xh2so4_nucl = so4rate / ( (2.0d+00/tstep) + kc )   ! use [h2so4] at mid-time step [ugso4/m^3]
         endif
         call NpfRate(var%pres(col),var%rh(col),var%tk(col),xh2so4_nucl,so4rate,xnh3,kc_aeq1,dndt,dmdt_so4,0)
         ! write(34,90010) var%tk(col),var%rh(col),ugm3_ncm3*xh2so4_init,ugm3_ncm3*xh2so4_nucl,kc,xnh3,1.0d-06*dndt,ih2so4_path
         if( write_log ) then
            write(aunit1,'(/a)')'new particle formation'                          
            write(aunit1,'(/a)')'pres,var%rh(col),var%tk(col),xh2so4_init(ug/m3),xh2so4_nucl(#/cm3),so4rate,xnh3,kc,dndt,dmdt_so4'
            write(aunit1,90003) var%pres(col),var%rh(col),var%tk(col),xh2so4_init,xh2so4_nucl*ugm3_ncm3,&
	                        so4rate,xnh3,kc,dndt,dmdt_so4
            write(aunit1,'(/a)')'var%tk(col),var%rh(col),h2so4(#.cm^3),h2so4_nucl(#/cm3),kc,nh3,dndt(#/cm3/s),ih2so4_path'      
            write(aunit1,90010) var%tk(col),var%rh(col),ugm3_ncm3*xh2so4_init,ugm3_ncm3*xh2so4_nucl,kc,&
	                        xnh3,1.0d-06*dndt,ih2so4_path
         endif
         else
         dndt     = 0.0d+00
         dmdt_so4 = 0.0d+00
         endif
   
   
         !----------------------------------------------------------------------------------------------------------------
         ! get the b_i loss       terms due to intermodal coagulation. [1/s]
         ! get the r_i production terms due to intermodal coagulation. [#/m^3/s]
         ! get the c_i terms, which include all source terms.          [#/m^3/s]
         ! for the c_i terms, the secondary particle formation term dndt must be 
         !   added in after coupling to condensation below.
         ! the a_i terms for intramodal coagulation are directly computed 
         !   from the coagulation coefficients when the number equations 
         !   are integrated.
         ! if dikl(i,k,l) = 0, then modes k and l to not coagulate to form mode i.
         ! dij(i,j) is unity if coagulation of mode i with mode j results
         !   in the removal of particles from mode i, and zero otherwise. 
         !----------------------------------------------------------------------------------------------------------------
         bi(:) = 0.0d+00
         ri(:) = 0.0d+00
         do ikl = 1, ndikl
         i = dikl_control(ikl)%i
         k = dikl_control(ikl)%k
         l = dikl_control(ikl)%l
         ri(i) = ri(i) + kbar0ij(k,l) * ni(k) * ni(l)
         end do
         do k=1, nweights
         do i = 1, nweights
            bi(i) = bi(i) + kbar0ij(i,k)*ni(k)*xdij(i,k) ! coagulation of modes i and j removes
         end do
         end do
         ci(:) = ri(:) + p_emis_numb(:)    ! all source terms for number concentration 
         diagtmp1(3,numb_map(:)) = ri(:)  
   
         if( write_log ) then
         write(aunit1,'(/6a20/)') 'i','bi [1/s]','ci [#/m^3/s]','ri [#/m^3/s]','ni [#/m^3]','kii[m^3/s]'
         do i=1, nweights
            write(aunit1,90004) i, bi(i), ci(i), ri(i), ni(i), kbar0ij(i,i)
         enddo
         endif
   
   
         !----------------------------------------------------------------------------------------------------------------
         !
         ! get the production and loss terms for the mass concentrations.
         !
         !----------------------------------------------------------------------------------------------------------------
         ! the production terms piq(nweights,nmass_spcs) are in [ug/m^3/s].
         !   the first  index runs over all modes or quadrature points.
         !   the second index runs over the nmass_spcs (=5) principal mass
         !     species: sulfate, bc, oc, dust, sea salt. 
         !
         ! get the emissions production terms (pemis_i,q in the manuscript)
         !   in [ug/m^3/s] and put them in the total production rate array piq
         !   (the p_i,q in the manuscript).
         !
         ! emis_mode_map and emis_spcs_map have elements corresponding to 
         !   the aerosol types (in this order): akk(=1), acc(=2), bcc(=8), occ(=7),
         !               dd1(=3), ssa(=5), ssc(=6), boc(bc=8), boc(oc=9), dd2(=10).
         ! emis_mode_map(j) is mode number receiving the emissions held 
         !                  in var%emis_mas(col,j).
         ! emis_spcs_map(j) is the chemical species number (1-5) of the species
         !                  held in var%emis_mas(col,j).
         ! currently, emis_spcs_map = (/1,1,2,3,4,5,5,2,3,4/)
         !----------------------------------------------------------------------------------------------------------------
         piq(:,:) = 0.0d+00          ! zero all production terms for this time step.
         do j=1, nemis_spcs          ! loop over the emitted species. 
         piq( emis_mode_map(j), emis_spcs_map(j) ) = &
         piq( emis_mode_map(j), emis_spcs_map(j) ) + var%emis_mas(col,j)
         enddo
         do i=1, nmodes
         do q=1, nm(i)
            diagtmp1(8,mass_map(i,q)) = piq(i,prod_index(i,q))  
         enddo
         enddo      
   
   
         if( write_log ) then
         write(aunit1,'(/a/)') 'emissions only: i, piq(i,:)'
         do i=1, nmodes
            ! write(aunit1,'(i5,5d14.4)') i, piq(i,:)
         enddo
         endif
   
   
         !-------------------------------------------------------------------------------------------------------------------
         ! get the pcoag_i,q terms for the production of mass species q in mode or quadrature point i 
         !   due to coagulation between modes or quadrature points k and l.
         ! these are at the same time mapped to the piq array.
         !
         ! qq = prod_index(i,q) is the principal species index that is the second index of the piq and mjq arrays.
         ! it points to chemical species chem_spc_name(q) where
         !   chem_spc_name(nmass_spcs) = (/'sulf','bcar','ocar','dust','seas'/)
         !-------------------------------------------------------------------------------------------------------------------
         ! write(36,'(a)')'i,q,k,l,piq(i,qq),kbar3ij(k,l)*ni(k)*ni(l)*(mjq(k,q)+mjq(l,q))'
         ! write(36,'(a)')'i,q,k,l,kbar3ij(k,l),ni(k),ni(l),mjq(k,q),mjq(l,q)'
         !-------------------------------------------------------------------------------------------------------------------
   
         piqtmp(:,:) = 0.0d+00
   
         !-----------------------------------------------------------------------------------------------------------------
         ! sum over all k-l interactions and identify which interactions produce species q in mode i.
         !
         ! there are three cases where mass of q is transferred to mode i, if the q concentration 
         ! in the donor mode(s) k and/or l is non-zero.
         !
         !   1. mode i is the same as mode k but different from mode l, and q-mass from mode l is 
         !      transferred to mode i.
         !   2. mode i is the same as mode l but different from mode k, and q-mass from mode k is
         !      transferred to mode i.
         !   3. mode i is different from either mode k or mode l, and q-mass from both mode k and mode l is 
         !      transferred to mode i.
         !
         ! the k-l double sum are symmetric in k-l, so we can sum over either the 
         ! superdiagonal or subdiagonal part of the 'matrix', but not both. note that 
         ! kbar3ij(k,l) is not symmetric in k-l, so kbar3ij(k,l) and kbar3ij(l,k) are not interchangeable.
         !---------------------------------------------------------------------------------------------------------------
   
         do i=1, nweights          ! loop over all modes or quadrature points receiving mass species q
         do klq = 1, giklq_control(i)%n
            k = giklq_control(i)%k(klq)
            l = giklq_control(i)%l(klq)
            qq = giklq_control(i)%qq(klq)
   
            piqtmp(i,qq) = piqtmp(i,qq) + &
                  ni(k)*ni(l) *kbar3ij(l,k) * mjq(l,qq)
         end do
         end do
   
            !-------------------------------------------------------------------------------------------------------------
            ! if( i.eq.15 .and. k.le.2 .and. l.eq.10 ) then
            !   write(36,'(4i3,3d15.6,i6)') i,q,k,l,piqtmp(i,qq),ni(k),ni(l),ibranch
            !   write(36,'(12x,5d15.6)')             kbar3ij(k,l),kbar3ij(l,k),mjq(k,qq),mjq(l,qq)
            !   write(36,'(12x,5d15.6)')             kbar3ij(k,l)*mjq(k,qq),kbar3ij(l,k)*mjq(l,qq)
            !   write(36,'(12x,5d15.6)')             kbar3ij(k,l)*mjq(k,qq)+kbar3ij(l,k)*mjq(l,qq)
            !   write(36,'(12x,5d15.6)')ni(k)*ni(l)*(kbar3ij(k,l)*mjq(k,qq)+kbar3ij(l,k)*mjq(l,qq))
            ! endif
            !-------------------------------------------------------------------------------------------------------------
         do i=1, nmodes
         do q=1, nm(i)
            diagtmp1( 8,mass_map(i,q)) = piq   (i,prod_index(i,q))  ! due to mass emissions
            diagtmp1(10,mass_map(i,q)) = piqtmp(i,prod_index(i,q))  ! due to coagulation
         enddo
         enddo      
         piq(:,:) = piq(:,:) + piqtmp(:,:)
   
         if( write_log ) then
         write(aunit1,'(/a/)') 'i, q, qq, piq(i,qq) [ug/m^3/s] - after coagulation mass terms'
         do i=1, nweights
         do q=1, nm(i)
            qq = prod_index(i,q)
            ! write(aunit1,90000) i, q, qq, piq(i,qq)
         enddo
         enddo
         endif
   
         !----------------------------------------------------------------------------------------------------------------
         ! get the f_i terms for the loss of species q in quadrature point i due
         !   to intermodal coagulation of i with other quadrature points. 
         ! as all species in quadrature point i are lost through this coagulation,
         !   the species index q is not needed.
         ! dij(i,j) is unity if coagulation of mode i with mode j results
         !   in the removal of particles from mode i, and zero otherwise. 
         !   dij(i,j) is not symmetric in i-j.
         ! kbar3ij(i,j) is not symmetric in i-j, and the first index i refers
         !   to the 'donor' mode in the coagulation interaction. since mode i
         !   is here losing mass due to coagulation, mode i is the 'donor' mode.
         !----------------------------------------------------------------------------------------------------------------
         ! write(35,'(/a/)') 'dij(i,j), i, j, kbar3ij(i,j), ni(j), fi(i)'
         !----------------------------------------------------------------------------------------------------------------
         fi(:) = 0.0d+00
         do i=1, nweights
         do j=1, nweights    ! for each mode i, we must sum over all modes j.
            if( dij(i,j).gt.0 ) fi(i) = fi(i) + kbar3ij(i,j)*ni(j)
            !------------------------------------------------------------------------------------------------------------
            ! write(35,'(3i4,3d15.7)') dij(i,j), i, j, kbar3ij(i,j), ni(j), fi(i)
            !------------------------------------------------------------------------------------------------------------
         enddo
         enddo
         do i=1, nmodes
         do q=1, nm(i)
            diagtmp1(13,mass_map(i,q)) = - fi(i) * var%aero(col, mass_map(i,q) )  
         enddo
         enddo
   
         if( write_log ) then
         write(aunit1,'(/a/)') 'i, fi [1/s], ni [#/m^3]'
         do i=1, nweights
            ! write(aunit1,90001) i, fi(i), ni(i)
         enddo
         endif
   
   
         !----------------------------------------------------------------------------------------------------------------
         ! get the pcloud_i,q terms due to the in-cloud production of sulfate
         ! for each mode or quadrature point in [ugso4/m^3/s].
         !
         ! aqso4rate is the in-cloud sulfate production rate [ug/m^3/s].
         ! 
         ! kappai(i) is the soluble (or activating) fraction for mode i, as specified in aero_param.f.
         ! 
         ! if all solule particles in mode i are assumed to activate, then
         ! the product kappai(i) * ni(i) * rsum_activ is the fraction of the
         ! in-cloud sulfate going into mode i.
         ! 
         ! if a droplet activation calculation is done for each mode, then 
         ! the product var%nact(col,i) * rsum_activ is the fraction of the
         ! in-cloud sulfate going into mode i. the geometric mean radii passed to
         ! getactfrac is for the dry aerosol. 
         !----------------------------------------------------------------------------------------------------------------
         piqtmp(:,:) = 0.0d+00
         IF( var%aqso4rate(col) .gt. aqso4rate_min ) THEN  
         
	  IF(     activation_scheme .eq. 1 ) then  
            nsol(:) = kappai(:)*ni(:)
            rsum_activ = 1.0d+00 / sum( nsol(:) + tinydenom ) 
            piqtmp(:,prod_index_sulf) = ( kappai(:)*ni(:)*rsum_activ ) * var%aqso4rate(col)
            !nactv(ixxx,iyyy,ilay,:) = nsol(:)       ! [#/m] - store for use outside this routine.
            !------------------------------------------------------------------------------------------------------------
            ! write(40,'(/a,f15.6/)')'total number soluble (#/cm^3) = ', 1.0d-06/rsum_activ 
            ! do i=1, nmodes
            !   write(40,'(i3,f5.2,3f12.4,a5)')i,kappai(i),1.0d-06*ni(i),1.0d-06*nsol(i),nsol(i)*rsum_activ,mode_name(i)
            ! enddo  
            !------------------------------------------------------------------------------------------------------------
          ELSEIF( activation_scheme .eq. 2 ) then
            do i=1, nmodes
               mi5(i,:) = mjq(i,:) * ni(i)     ! mass conc. of each component in mode i [ug/m^3]
            enddo  
            call GetActFrac(nmodes,ni,mi5,0.5d+00*dgn_dry,sig0,var%tk(col),var%pres(col),var%wupdraft(col), &
                            var%ac(col,:),fracactn,fracactm,var%nact(col,:),mact)
            rsum_activ = 1.0d+00 / sum ( var%nact(col,:) + tinydenom ) 
            piqtmp(:,prod_index_sulf) = ( var%nact(col,:)*rsum_activ ) * var%aqso4rate(col)
            !LFR ATENCAO  nactv(ixxx,iyyy,ilay,:) = var%nact(col,:)       ! [#/m] - store for use outside this routine.
            !------------------------------------------------------------------------------------------------------------
            ! write(40,'(/a,f15.6/)')'total number activated (#/cm^3) = ', 1.0d-06/rsum_activ 
            ! do i=1, nmodes
            !   write(40,90009)i,fracactn(i),1.0d-06*ni(i),1.0d-06*var%nact(col,i),1.0d-06*var%nact(col,i),var%nact(col,i)*rsum_activ,
            !&                       mode_name(i),dgn(i),dgn_dry(i)
            ! enddo  
            !------------------------------------------------------------------------------------------------------------

          !-KML
          ELSEIF( activation_scheme .eq. 3 ) then
	   kappapk   (:)= 0.0
	   vi5       (:)= 0.0
	   
           do i=1, nmodes
		do j=1,nmass_spcs
		  if ( mjq(i,j) .lt.  1.0D-014 ) then
		     mi5(i,j) = mjq(i,j) * tinynumer
                  else
		     mi5(i,j) = mjq(i,j) * ni(i)     ! mass conc. of each component in mode i [ug/m^3]
		  endif
		enddo
		 !- partial volume for each mode
		 vi5(i) = vi5(i) + (4.0*PI/3.0)* dp_dry(i)**3
	    enddo
	    
            !- get kappa parameter for each mode	    
	    do i=1, nmodes
                               
		IF ( mode_name(i) .eq. 'AKK') kappapk(i) = kappapk_SULF         ! kappapk - hygroscopicity parameter
                IF ( mode_name(i) .eq. 'ACC') kappapk(i) = kappapk_SULF
                IF ( mode_name(i) .eq. 'DD1') kappapk(i) = kappapk_DUST
                
		IF ( mode_name(i) .eq. 'DD1' .OR. mode_name(i) .eq.'DD2' &
         .OR. mode_name(i) .eq.'DS1' .OR. mode_name(i) .eq.'DS2') then
		  kappapk(i) = (kappapk_DUST*mi5(i,4)/emis_dens_dust+  &
		                kappapk_SULF*mi5(i,1)/emis_dens_sulf)/vi5(i)
                ENDIF
		
		IF ( mode_name(i) .eq. 'SSA' .OR. mode_name(i) .eq.'SSC' .OR. mode_name(i) .eq.'SSS') then
		  kappapk(i) = (kappapk_SEAS*mi5(i,5)/emis_dens_seas+  &
		                kappapk_SULF*mi5(i,1)/emis_dens_sulf)/vi5(i)
                ENDIF


		IF ( mode_name(i) .eq. 'OCC' .OR. mode_name(i) .eq. 'OCS') then
		  kappapk(i) = (kappapk_OCAR*mi5(i,3)/emis_dens_ocar+  &
		                kappapk_SULF*mi5(i,1)/emis_dens_sulf)/vi5(i)
		ENDIF

		IF ( mode_name(i) .eq. 'BC1' .OR. mode_name(i) .eq.'BC2' &
          .OR. mode_name(i) .eq.'BC3'.OR. mode_name(i) .eq.'BCS') then
		  kappapk(i) = (kappapk_BCAR*mi5(i,2)/emis_dens_bcar+  &
		                kappapk_SULF*mi5(i,1)/emis_dens_sulf)/vi5(i)
                ENDIF

		IF ( mode_name(i) .eq. 'DBC') then
		  kappapk(i) = (kappapk_BCAR*mi5(i,2)/emis_dens_bcar+ &
		                kappapk_DUST*mi5(i,4)/emis_dens_dust+ &
				kappapk_SULF*mi5(i,1)/emis_dens_sulf)/vi5(i)
		ENDIF
		
		IF ( mode_name(i) .eq. 'BOC') then
		  kappapk(i) = (kappapk_BCAR*mi5(i,2)/emis_dens_bcar+  &
		                kappapk_OCAR*mi5(i,3)/emis_dens_ocar+  &
				kappapk_SULF*mi5(i,1)/emis_dens_sulf)/vi5(i)
                ENDIF
		
                IF ( mode_name(i) .eq. 'MXX') then         
		  kappapk(i) = (kappapk_BCAR*mi5(i,2)/emis_dens_bcar + &
		                kappapk_OCAR*mi5(i,3)/emis_dens_ocar + &
				kappapk_DUST*mi5(i,4)/emis_dens_dust + &
                                kappapk_SEAS*mi5(i,5)/emis_dens_seas + &
				kappapk_SULF*mi5(i,1)/emis_dens_sulf)/vi5(i) 
                ENDIF
            enddo
           
	    number_in     = 0.0	    
	    kappapk_eff(:)= 0.0	   
	    diam_eff   (:)= 0.0	   
	    ni_total   (:)= 0.0	   
            vol_total  (:)= 0.0	
	    DO i=1, nmodes
	     
	      ! -- ice nuclei particles
	      IF(kappapk(i) < 0.1 .AND. dp_dry(i) >= 0.1) THEN
	       
	        number_in = number_in + ni(i)
	     
	      ELSE
	      ! -- cloud water nuclei particles
	     
	        IF    ( (kappapk(i) >= 0.1 .AND. kappapk(i) <= 0.2)) THEN
	       
	         IF(dp_dry(i) <= 0.1) THEN
		  !- effective values for nucleation mode with low hygroscopicity
		  vol_total  (1) = vol_total  (1) + vi5(i)
		  ni_total   (1) = ni_total   (1) + ni(i)
		  kappapk_eff(1) = kappapk_eff(1) + kappapk(i)*vi5(i)
		  diam_eff   (1) = diam_eff   (1) + dp_dry(i)*ni(i)
		  
	         ELSEIF( dp_dry (i) >0.1 .AND. dp_dry (i) < 2.0) THEN
	        
		  !- effective values for accumulation mode with low hygroscopicity
		  vol_total  (3) = vol_total  (3) + vi5(i)
		  ni_total   (3) = ni_total   (3) + ni(i)
		  kappapk_eff(3) = kappapk_eff(3) + kappapk(i)*vi5(i)
		  diam_eff   (3) = diam_eff   (3) + dp_dry(i)*ni(i)
		  
		 ENDIF
	     
	        ELSEIF ( (kappapk(i) > 0.2) ) THEN
	       
	         IF(dp_dry(i) <=  0.1) THEN
	          
		  !- effective values for nucleation mode with high hygroscopicity
		  vol_total  (2) = vol_total  (2) + vi5(i)
		  ni_total   (2) = ni_total   (2) + ni(i)
		  kappapk_eff(2) = kappapk_eff(2) + kappapk(i)*vi5(i)
		  diam_eff   (2) = diam_eff   (2) + dp_dry(i)*ni(i)

	         ELSEIF( dp_dry (i) >0.1 .AND. dp_dry (i) <= 2.0) THEN
		
		  !- effective values for accumulation mode with high hygroscopicity
		  vol_total  (4) = vol_total  (4) + vi5(i)
		  ni_total   (4) = ni_total   (4) + ni(i)
		  kappapk_eff(4) = kappapk_eff(4) + kappapk(i)*vi5(i)
		  diam_eff   (4) = diam_eff   (4) + dp_dry(i)*ni(i)
	        
		 ENDIF
	       ENDIF
	      ENDIF
            ENDDO

	    DO i = 1, 4
               kappapk_eff(i)= kappapk_eff(i)/(vol_total(i)+tinynumer)
               diam_eff   (i)= diam_eff   (i)/(ni_total (i)+tinynumer)
            ENDDO
            
	    !- for now we are sending only the accumulation mode to the microphysics combining
	    !- low and high hygroscopicity
	    !-- effective kappa
	    var%aer2mp_eff(col,1) = (kappapk_eff(3)*vol_total(3)+ kappapk_eff(4)*vol_total(4)) &
	                          / (vol_total(3)+ vol_total(4) + tinynumer)
	    
	    !--total number of water friendly aerosol
	    var%aer2mp_eff(col,3) = ni_total(3) + ni_total(4) + tinynumer
	    
	    !-- effective diameter
	    var%aer2mp_eff(col,2) = ( diam_eff(3)*ni_total(3) + diam_eff(4)*ni_total(4) ) &
	                          /   var%aer2mp_eff(col,3)
	   	      
	    !--total number of ice friendly aerosol
	    var%aer2mp_eff(col,4) = number_in
	       
          ENDIF 

         ELSE                                       ! put all in-cloud sulfate into the accumulation mode.
           if( numb_akk_1 .ne. 0 ) then             ! mode akk exists and mode acc has mode number = 2.
            piqtmp(2,prod_index_sulf) = var%aqso4rate(col)  ! [ugso4/m^3/s]
           else                                     ! mode akk does not exist and mode acc has mode number = 1.
            piqtmp(1,prod_index_sulf) = var%aqso4rate(col)  ! [ugso4/m^3/s]
           endif      
         ENDIF
	 
         piq(:,prod_index_sulf) = piq(:,prod_index_sulf) + piqtmp(:,prod_index_sulf) 
         diagtmp1(12,mass_map(:,prod_index_sulf)) = piqtmp(:,prod_index_sulf)  
   
         !----------------------------------------------------------------------------------------------------------------
         ! get the pgrowth_i,q terms due to condensation and gas-particle 
         !   mass transfer for each mode (or quadrature point) and add them to the piq array. 
         !
         ! the net loss of h2so4 due to both secondary particle formation
         !   and condensation should not exceed the current h2so4 concentration.
         !   this is enforced by rescaling the two h2so4 consumption rates in a way 
         !   that preserves the relative magnitudes of these two loss processes. 
         !   the net condensation rate pq_growth is calculated from xh2so4_nucl, not
         !   the total accumulated h2so4 in xh2so4_init, for balance with the 
         !   treatment of new particle formation above. 
         !
         ! the expression in parentheses on the rhs of piq is that for h_i, 
         !
         !   h_i = kci_coef_dp(i,ilay) * ni(i) / kc
         !
         !   the ratio of the condensational sink of mode or quadrature point i 
         !   to the total condensational sink.
         !
         ! at this point, xh2so4_init = var%gas(col, gas_h2so4 ) + tinynumer.
         !
         ! the h2so4 concentration var%gas(col, gas_h2so4 ) is updated here.
         !----------------------------------------------------------------------------------------------------------------
         pq_growth = xh2so4_nucl * ( 1.0d+00 - exp(-kc*tstep) ) / tstep            ! [ugso4/m^3/s]
         tot_h2so4_loss = ( dmdt_so4 + pq_growth ) * tstep                         ! [ugso4/m^3]
         if ( tot_h2so4_loss .gt. xh2so4_init ) then                               ! xh2so4_init=var%gas(col,gas_h2so4)+tinynumer
         dmdt_so4  = dmdt_so4  * ( xh2so4_init / ( tot_h2so4_loss + tinydenom ) )! [ugso4/m^3/s]
         dndt      = dndt      * ( xh2so4_init / ( tot_h2so4_loss + tinydenom ) )! [  #  /m^3/s]
         pq_growth = pq_growth * ( xh2so4_init / ( tot_h2so4_loss + tinydenom ) )! [ugso4/m^3/s]
         var%gas(col, gas_h2so4 ) = tinynumer
         else
         var%gas(col, gas_h2so4 ) = var%gas(col, gas_h2so4 ) - tot_h2so4_loss + tinynumer
         endif
         diagtmp1(2,numb_map(1))               = dndt  
         diagtmp1(9,sulf_map(prod_index_sulf)) = dmdt_so4  
         
         ! write(35,'(8d13.5)')piq(:,prod_index_sulf), kci_coef_dp(:,ilay),ni(:),kc
         
         piqtmp(:,prod_index_sulf) = ( kci_coef_dp(:,ilay)*ni(:)/kc ) * pq_growth
         piq   (:,prod_index_sulf) = piq(:,prod_index_sulf) + piqtmp(:,prod_index_sulf) 
   
         diagtmp1(11,mass_map(:,prod_index_sulf)) = piqtmp(:,prod_index_sulf)  
         !-----------------------------------------------------------------------------------------------------------------
         ! add the secondary particle formation term dndt [#/m^3/s] for the number concentration term.
         ! add the secondary particle formation term dmdt_so4 [ug/m^3/s], calculated above, 
         !   to the total mass production rate array piq for the akk mode.
         ! if the akk mode is absent, the secondary particle formation terms still go into mode 1.
         !-----------------------------------------------------------------------------------------------------------------
         ci(1) = ci(1) + dndt                                        ! add secondary particle formation number term
         piq(1,prod_index_sulf) = piq(1,prod_index_sulf) + dmdt_so4  ! add secondary particle formation mass   term
   
         if( write_log ) then
         write(aunit1,'(/a,5x,3d15.8)')'xh2so4_init, xh2so4_nucl, pq_growth = ', xh2so4_init, xh2so4_nucl, pq_growth
         write(aunit1,*)'piq(1,prod_index_sulf) = ', piq(1,prod_index_sulf)
         endif
   
   
         !----------------------------------------------------------------------------------------------------------------
         !
         ! solve for the updated number concentrations (quadrature weights). 
         !
         !----------------------------------------------------------------------------------------------------------------
         do i=1, nweights
         y0 = ni(i)                                  ! initial number concentration [#/m^3/s].
         a = 0.5d+00*kbar0ij(i,i)
         b = bi(i)   
         c = ci(i)   
         if( c .gt. 1.0d-30 ) then
            delta = sqrt( b * b + 4.0d+00 * a * c )    
            r1 = 2.0d+00 * a * c / ( b + delta )       
            r2 = - 0.5d0 * ( b + delta )             
            gamma =  - ( r1 - a * y0 ) / ( r2 - a * y0 )
            gexpdt = gamma * exp( - delta * tstep )
            y = (  r1 + r2 * gexpdt ) / ( a * ( 1.0d+00 + gexpdt ) )
         else                                        ! when c = 0.0d+00, as we assume c is not negative.
            expdt = exp( - b * tstep )
            if( 1.0d+00-expdt .gt. piq_thresh ) then                     ! if( expdt .lt. 1.0d+00 ) then
               y = b * y0 * expdt / ( b + a * y0 * ( 1.0d+00 - expdt ) )
            else
               y = y0 / ( 1.0d+00 + a * y0 * tstep )
            endif
         endif
         var%aero(col, numb_map(i) ) = max ( y, minconc )   ! update the output array, not the work array ni(i)
         diagtmp1(4,numb_map(i)) = -b*y0  
         diagtmp1(5,numb_map(i)) = -a*y0*y0
         ! write(*,'(4d15.5)') a,b,c,y0  
         enddo
   
         if( write_log ) then
         write(aunit1,'(/a6,2a30/)') 'i','ni(i) [#/m^3]','var%aero(col, numb_map(i) ) [#/m^3]'
         do i=1, nweights
            ! write(aunit1,90006) i, ni(i), aero ( numb_map(i) )
         enddo
         endif
   
   
         !----------------------------------------------------------------------------------------------------------------
         !
         ! solve for the updated mass concentrations. 
         !
         !----------------------------------------------------------------------------------------------------------------
         ! update the sulfate, bc, oc, dust, and sea salt concentrations.
         !
         ! mass_map(i,q) is the location in aero(:) of the qth mass in mode i.
         ! mode i has nm(i) mass species defined for it, and nm(i) varies between
         !   1 and nmass_spcs (=5).
         !
         ! the second index of prod_index has nmass_spcs (=5) values:
         !   1=sulfate, 2=bc, 3=oc, 4=dust, 5=sea salt.
         !   prod_index(i,q) is the location in array piq(i,q) of chemical species
         !   chem_spc_name(q) for mode (quadrature weight) i.
         !
         ! if ( 1.0d+00-expdt .gt. piq_thresh ) --> the first (loss) term in aero update may be significant.
         !----------------------------------------------------------------------------------------------------------------
         !LFR-> see eq (10)
         do i=1, nweights
         expdt = exp( - fi(i) * tstep ) 
         if ( 1.0d+00-expdt .gt. piq_thresh ) then 
            factor = ( 1.0d+00 - expdt ) / fi(i)
            do q=1, nm(i)
               var%aero(col,mass_map(i,q)) = var%aero(col,mass_map(i,q)) * expdt + piq(i,prod_index(i,q)) * factor
            enddo
         else
            do q=1, nm(i)
               var%aero(col,mass_map(i,q)) = var%aero(col,mass_map(i,q)) + piq(i,prod_index(i,q)) * tstep
            enddo
         endif
         enddo
   
   
         !----------------------------------------------------------------------------------------------------------------
         ! update the nitrate, ammonium and aerosol water concentrations.
         !
         ! call the thermodynamic module again to determine the bulk gas-particle 
         ! partitioning of the inorganic species and the water content
         ! associated with those species. also determine the water content
         ! associated with the nacl of sea salt.
         !
         ! note that aero(mass_h2o) must include the sea-salt associated water 
         !   upon exit from this routine. 
         ! upon return from aero_thermo, aero(mass_h2o) contains only the 
         !   non-sea salt-associated water. the sea salt-associated water is in ssh2o.
         !----------------------------------------------------------------------------------------------------------------
         tot_sulf = sum( var%aero(col,sulf_map(:)) )
         tot_dust = sum( var%aero(col,dust_map(:)) )
         tot_seas = sum( var%aero(col,seas_map(:)) )
         if( write_log ) write(aunit1,'(/a,f12.4/)') 'tot_sulf = ', tot_sulf
         aero_water_actual = var%aero(col,mass_h2o) + ssh2o
         call Aero_Thermo(tot_sulf,var%aero(col,mass_no3),var%aero(col,mass_nh4),var%aero(col,mass_h2o),var%gas(col,gas_nh3), &
                        var%gas(col,gas_hno3),tot_dust,tot_seas,ssh2o,var%tk(col),var%rh(col),var%pres(col),rhd,rhc)
         aero_water_wet = var%aero(col,mass_h2o) + ssh2o
   
         !----------------------------------------------------------------------------------------------------------------
         ! again, adjust the water concentration for hysteresis.
         !----------------------------------------------------------------------------------------------------------------
         if ( var%rh(col) .gt. rhc  .and.  var%rh(col) .lt. rhd ) then
         if ( aero_water_wet .gt. 0.0d+00 ) then
            if ( aero_water_actual/aero_water_wet .lt. 0.5d+00 ) then
               var%aero(col,mass_h2o) = 0.0d+00  ! zero the non-sea salt-associated water.
               ssh2o = 0.0d+00           ! zero the     sea salt-associated water.
            endif
         endif
         elseif ( var%rh(col) .le. rhc ) then  
         var%aero(col,mass_h2o) = 0.0d+00
         ssh2o = 0.0d+00
         endif
         !----------------------------------------------------------------------------------------------------------------
         ! the total aerosol water must exit the routine in aero(mass_h2o).
         !----------------------------------------------------------------------------------------------------------------
         var%aero(col,mass_h2o) = var%aero(col,mass_h2o) + ssh2o         ! total aerosol water now.
   
   
         !----------------------------------------------------------------------------------------------------------------
         ! transfer mass and number concentrations from dd1 to ds1, dd2 to ds2,
         !   bc1 to bc2, and bc2 to bc3, if the appropriate modes are defined 
         !   in the present configuration.
         !
         ! this transfer is based upon the volume fraction of inorganic constituents.
         !   for computational efficiency, the maximum inorganic volume fraction (mivf)
         !   for a mode is transformed into the maximum inorganic mass ratio (mimr), 
         !   and the mimr are precomputed parameters.
         !
         ! in the following if statements, the first (and lengthy) quantity is the current
         !   inorganic-to-dust or inorganic-to-bc mass ratio. for example, in the second if statement below,  
         !
         !     aero(mass_dd1_sulf)*optot_no3nh4h2o_to_sulf 
         !
         !   is the total mass concentration of inorganic coating 
         !   (sulfate + nitrate + ammonium + water) in mode dd1, and 
         !
         !     aero(mass_dd1_sulf)*optot_no3nh4h2o_to_sulf/aero(mass_dd1_dust) is the 
         !
         !   inorganic-to-dust mass ratio for mode dd1.
         !
         ! there is no need to explicitly transfer nitrate, ammonium, or aerosol water. 
         !
         ! optot_no3nh4h2o_to_sulf is the total no3+nh4+h2o mass per unit mass so4, plus 1.
         !----------------------------------------------------------------------------------------------------------------
         aerotmp2(:) = var%aero(col,:) 
         var%aero(col,mass_dd1_dust) = max( var%aero(col,mass_dd1_dust), tinynumer )
         var%aero(col,mass_bc1_bcar) = max( var%aero(col,mass_bc1_bcar), tinynumer )
         var%aero(col,mass_bc2_bcar) = max( var%aero(col,mass_bc2_bcar), tinynumer )
         if( mass_dd2_dust .gt. 0 ) var%aero(col,mass_dd2_dust) = max( var%aero(col,mass_dd2_dust), tinynumer )   ! if mode dd2 exists. 
   
         tot_sulf = sum( var%aero(col,sulf_map(:)) ) + tinynumer
         optot_no3nh4h2o_to_sulf = 1.0d+00 + sum(var%aero(col,1:3)) / tot_sulf
   
         if( var%aero(col,mass_dd1_sulf)*optot_no3nh4h2o_to_sulf/var%aero(col,mass_dd1_dust) .gt. mimr_ddd ) then
         !--------------------------------------------------------------------------------------------------------------
         ! transfer mode dd1 to mode ds1.
         !--------------------------------------------------------------------------------------------------------------
         var%aero(col,mass_ds1_sulf) = var%aero(col,mass_ds1_sulf) + var%aero(col,mass_dd1_sulf)
         var%aero(col,mass_ds1_dust) = var%aero(col,mass_ds1_dust) + var%aero(col,mass_dd1_dust)
         var%aero(col,numb_ds1_1   ) = var%aero(col,numb_ds1_1   ) + var%aero(col,numb_dd1_1   )
         var%aero(col,mass_dd1_sulf) = tinynumer
         var%aero(col,mass_dd1_dust) = tinynumer
         var%aero(col,numb_dd1_1   ) = tinynumer
         endif
   
         if( mass_dd2_dust .gt. 0.0d+00 ) then
         if( var%aero(col,mass_dd2_sulf)*optot_no3nh4h2o_to_sulf/var%aero(col,mass_dd2_dust) .gt. mimr_ddd ) then
            !------------------------------------------------------------------------------------------------------------
            ! transfer mode dd2 to mode ds2.
            !------------------------------------------------------------------------------------------------------------
            var%aero(col,mass_ds2_sulf) = var%aero(col,mass_ds2_sulf) + var%aero(col,mass_dd2_sulf)
            var%aero(col,mass_ds2_dust) = var%aero(col,mass_ds2_dust) + var%aero(col,mass_dd2_dust)
            var%aero(col,numb_ds2_1   ) = var%aero(col,numb_ds2_1   ) + var%aero(col,numb_dd2_1   )
            var%aero(col,mass_dd2_sulf) = tinynumer
            var%aero(col,mass_dd2_dust) = tinynumer
            var%aero(col,numb_dd2_1   ) = tinynumer
         endif
         endif
   
         if( var%aero(col,mass_bc1_sulf)*optot_no3nh4h2o_to_sulf/var%aero(col,mass_bc1_bcar) .gt. mimr_bc1 ) then
         !--------------------------------------------------------------------------------------------------------------
         ! transfer mode bc1 to mode bc2.
         !--------------------------------------------------------------------------------------------------------------
         var%aero(col,mass_bc2_sulf) = var%aero(col,mass_bc2_sulf) + var%aero(col,mass_bc1_sulf)
         var%aero(col,mass_bc2_bcar) = var%aero(col,mass_bc2_bcar) + var%aero(col,mass_bc1_bcar)
         var%aero(col,numb_bc2_1   ) = var%aero(col,numb_bc2_1   ) + var%aero(col,numb_bc1_1   )
         var%aero(col,mass_bc1_sulf) = tinynumer
         var%aero(col,mass_bc1_bcar) = tinynumer
         var%aero(col,numb_bc1_1   ) = tinynumer
         endif
   
         if( var%aero(col,mass_bc2_sulf)*optot_no3nh4h2o_to_sulf/var%aero(col,mass_bc2_bcar) .gt. mimr_bc2 ) then
         !--------------------------------------------------------------------------------------------------------------
         ! transfer mode bc2 to mode bc3.
         !--------------------------------------------------------------------------------------------------------------
         if( include_bc3 ) then
            var%aero(col,mass_bc3_sulf) = var%aero(col,mass_bc3_sulf) + var%aero(col,mass_bc2_sulf)
            var%aero(col,mass_bc3_bcar) = var%aero(col,mass_bc3_bcar) + var%aero(col,mass_bc2_bcar)
            var%aero(col,numb_bc3_1   ) = var%aero(col,numb_bc3_1   ) + var%aero(col,numb_bc2_1   )
            var%aero(col,mass_bc2_sulf) = tinynumer
            var%aero(col,mass_bc2_bcar) = tinynumer
            var%aero(col,numb_bc2_1   ) = tinynumer
         endif
         endif
   
         if( write_log ) then
         write(aunit1,'(/a,2f15.6)')'mivf_ddd, mimr_ddd = ', mivf_ddd, mimr_ddd 
         write(aunit1,'( a,2f15.6)')'mivf_bc1, mimr_bc1 = ', mivf_bc1, mimr_bc1 
         write(aunit1,'( a,2f15.6)')'mivf_bc2, mimr_bc2 = ', mivf_bc2, mimr_bc2 
         write(aunit1,'(/a,2d15.6)')'tot_sulf, optot_no3nh4h2o_to_sulf = ', tot_sulf, optot_no3nh4h2o_to_sulf
         endif
   
         !----------------------------------------------------------------------------------------------------------------
         ! intermodal transfer from the aitken (akk) mode to the accumulation (acc) mode.
         !
         ! this is not a physical process, but a reclassification of particles. see binkowski and roselle (2003).
         !
         ! aero(mass_axx_sulf)*optot_no3nh4h2o_to_sulf is the total mass in mode axx (x=k, c) [ug/m^3]. 
         !   division by aero(numb_axx_1) converts to mean mass per particle [ug].
         !   multiplication by conv_mass_to_dp converts to dp^3 [m^3].
         !   taking the cube root yields the diameter of average mass for the mode [m].
         !
         ! fnum is the fraction of the akk mode mass and number concentrations transferred to the acc mode
         !      during this time step. binkowski and roselle (2003) limit fnum to a maximum of 0.5 for
         !      numerical stability, and the same is done here.
         !----------------------------------------------------------------------------------------------------------------
         if( intermodal_transfer .and. var%aero(col, numb_map(1) ) .gt. akk_minnum_imtr ) then
         dpakk = ( conv_mass_to_dp * var%aero(col, mass_akk_sulf ) &              ! diameter of average mass for akk [m]
                  * optot_no3nh4h2o_to_sulf / var%aero(col, numb_akk_1 ) )**0.33333333333
         dpacc = ( conv_mass_to_dp * var%aero(col, mass_acc_sulf ) &              ! diameter of average mass for acc [m]
                  * optot_no3nh4h2o_to_sulf / var%aero(col, numb_acc_1 ) )**0.33333333333
         if(     imtr_method .eq. 1 ) then
            !------------------------------------------------------------------------------------------------------------
            ! calculate the fraction transferred based on the relative difference 
            !   in mass mean diameters of the akk and acc modes. 
            !------------------------------------------------------------------------------------------------------------
            if( dpakk .ge. dpakk0 ) then                                     ! [m]
               fnum = ( ( dpakk - dpakk0 ) / ( dpacc - dpakk0 ) )**imtr_exp   ! fraction transferred from akk to acc
               fnum = max( min( fnum, fnum_max ), 0.0d+00 )                   ! limit transfer in a single transfer
            else
               fnum = 0.0d+00
            endif
            f3 = fnum
            ! write(34,'(7d12.4)')fnum,f3
         elseif( imtr_method .eq. 2 ) then
            !------------------------------------------------------------------------------------------------------------
            ! calculate the fraction transferred based on a fixed 
            !   threshold diameter dpcut_imtr. 
            !------------------------------------------------------------------------------------------------------------
            dgn_akk_imtr = 1.0d+06 * dpakk * conv_dpam_to_dgn(1)           ! [um]
            xnum = xnum_factor * log( dpcut_imtr / dgn_akk_imtr )          ! [1]
            xnum = max( xnum, x3_term )                                    ! limit for stability as in bs2003
            x3 = xnum - x3_term                                            ! [1]
            fnum = 0.5d+00 * erfc( xnum )                                  ! number fraction transferred from akk to acc
            f3   = 0.5d+00 * erfc( x3   )                                  ! mass   fraction transferred from akk to acc
            ! write(34,'(9d12.4)')dgn_akk_imtr,dpcut_imtr,dpakk*1.0d+06,var%aero(col,numb_akk_1),var%aero(col,mass_akk_sulf),fnum,f3
         elseif( imtr_method .eq. 3 ) then
            !------------------------------------------------------------------------------------------------------------
            ! calculate the fraction transferred based on the 
            !   diameter of intersection of the akk and acc modes. 
            !------------------------------------------------------------------------------------------------------------
            dgn_akk_imtr = 1.0d+06 * dpakk * conv_dpam_to_dgn(1)           ! [um]
            dgn_acc_imtr = 1.0d+06 * dpacc * conv_dpam_to_dgn(2)           ! [um]
            if( var%aero(col, numb_acc_1 ) .gt. 1.0d+06 ) then                     ! mode acc not essentially empty
               xnum = getxnum(var%aero(col,numb_akk_1),var%aero(col,numb_acc_1),dgn_akk_imtr,dgn_acc_imtr,lnsg_akk,lnsg_acc)   ! [1]
            else                                                           ! mode acc essentially empty - use method 2
               xnum = xnum_factor * log( dpcut_imtr / dgn_akk_imtr )        ! [1]
            endif
            xnum = max( xnum, x3_term )                                    ! limit for stability as in bs2003
            x3 = xnum - x3_term                                            ! [1]
            fnum = 0.5d+00 * erfc( xnum )                                  ! number fraction transferred from akk to acc
            f3   = 0.5d+00 * erfc( x3   )                                  ! mass   fraction transferred from akk to acc
            ! write(34,'(9d12.4)')dgn_akk_imtr,dgn_acc_imtr,var%aero(col,numb_akk_1),var%aero(col,numb_acc_1),fnum,f3
         endif
         !---srf
	 !f3  =f3*tstep/3600.0D+00
	 !fnum=fnum*tstep/3600.0D+00
	 !---srf
	 del_mass = var%aero(col, mass_akk_sulf ) * f3                            ! mass   concentration transferred [ug/m^3]
         del_numb = var%aero(col, numb_akk_1    ) * fnum                          ! number concentration transferred [# /m^3]
         var%aero(col, mass_akk_sulf ) = var%aero(col, mass_akk_sulf ) - del_mass         ! update akk mass   concentration  [ug/m^3]
         var%aero(col, numb_akk_1    ) = var%aero(col, numb_akk_1    ) - del_numb         ! update akk number concentration  [# /m^3]
         var%aero(col, mass_acc_sulf ) = var%aero(col, mass_acc_sulf ) + del_mass         ! update acc mass   concentration  [ug/m^3]
         var%aero(col, numb_acc_1    ) = var%aero(col, numb_acc_1    ) + del_numb         ! update acc number concentration  [# /m^3]
         if( write_log ) then
            write(aunit1,'(/a5,9d11.3)') 'imtr:',dpakk,dpacc,dpakk0,fnum,f3, &
                                                del_mass,del_numb,var%aero(col, numb_akk_1 ), var%aero(col, numb_acc_1 )
         endif
         endif
   
         diagtmp1( 6,:) = ( var%aero(col,:) - aerotmp2(:) ) / tstep
         diagtmp1(14,:) = ( var%aero(col,:) - aerotmp2(:) ) / tstep
         
         !----------------------------------------------------------------------------------------------------------------
         ! put all sea salt sulfate back into the accumulation mode (ssa) if the
         !   mechanism uses both the ssa and ssc modes.
         !----------------------------------------------------------------------------------------------------------------
         if ( number_of_seasalt_modes .eq. 2 ) then
         var%aero(col, mass_ssa_sulf ) = var%aero(col, mass_ssa_sulf ) +  var%aero(col, mass_ssc_sulf )
         var%aero(col, mass_ssc_sulf ) = tinynumer
         endif
   
         !----------------------------------------------------------------------------------------------------------------
         ! get final total mass concentration for each model species [ug/m^3].
         ! adjust final aerosol and gas-phase species by rescaling to enforce mass 
         !   conservation to machine precision.
         ! precise mass conservation to machine precision is not necessarily
         !   conserved otherwise due to formulation of the model equations
         !   in terms of production and loss terms that may not precisely
         !   cancel on some occasions that they should due to their distinct
         !   treatment.
         !--------------------- -------------------------------------------------------------------------------------------
         ! write(*,*)'tot_sulf 1 = ',sum( var%aero(col, sulf_map(:)) )
         if( mass_adj ) then
         call SpcMasses(var%aero(col,:),var%gas(col,:),spcmass2)
         call MassAdj(var%aero(col,:),var%gas(col,:),spcmass1,spcmass2,var%emis_mas(col,:), &
                      var%aqso4rate(col),tstep)
         endif
         ! write(*,*)'tot_sulf 2 = ',sum( var%aero(col, sulf_map(:)) )
   
         !----------------------------------------------------------------------------------------------------------------
         ! limit low mass or number concentrations.
         !----------------------------------------------------------------------------------------------------------------
         var%aero(col,:) = max( var%aero(col,:), minconc )
   
         !----------------------------------------------------------------------------------------------------------------
         ! budget diagnostics.
         !----------------------------------------------------------------------------------------------------------------
         diagtmp1(7,:)          = ( var%aero(col,:) - aerotmp1(:) ) / tstep   ! actual total differences for both mass and number.
         var%diag(col,:,:)              = 0.0d+00 
         var%diag(col,1: 7,numb_map(:)) = diagtmp1(1: 7,numb_map(:))  ! save diagnostics 1-7  for the number concentrations only.  
         var%diag(col,8:14,:)           = diagtmp1(8:14,:)            ! save diagnostics 8-15 for the mass   concentrations only.  
         var%diag(col,15,:)             = diagtmp1(7,:)               ! ... cont'd ...                      
         var%diag(col,8:15,numb_map(:)) = 0.0d+00                     ! zero diagnostics 8-15 for the number concentrations.
   
         !----------------------------------------------------------------------------------------------------------------
         ! rescale each term of the budget for species i to get the correct sum. 
         !----------------------------------------------------------------------------------------------------------------
         !LFR  c      do i=1, naerobox
         !LFR  c        aerotmp2(i) = 0.0d+00
         !LFR  c        do k=1, 6 
         !LFR  c         aerotmp2(i) = aerotmp2(i) + var%diag(col,k,i) + var%diag(col,k+7,i)
         !LFR  c        enddo 
         !LFR  c        aerotmp2(i) = aerotmp2(i) + var%diag(col,14,i)
         !LFR  c         do j=1, 6
         !LFR  c          var%diag(col,j,i) = var%diag(col,j,i) * ( diagtmp1(7,i) / ( aerotmp2(i) + tinydenom ) )
         !LFR  c        enddo 
         !LFR  c        do j=8, ndiag_aero-1
         !LFR  c          var%diag(col,j,i) = var%diag(col,j,i) * ( diagtmp1(7,i) / ( aerotmp2(i) + tinydenom ) )
         !LFR  c        enddo 
         !LFR  c      enddo 
   
         !LFR
         !Saving the radius os each particle and the density in aer1_list for use in drydep
         !CALL copyDensRadius(ac,dens_mode_dry,nmodes)
   
   
         90000 format(3i6,d15.5)
         90001 format(1i6,2d15.5)
         90002 format(1i6,d15.5)
         90003 format(f9.1,f6.3,f7.2,7d13.3)
         90004 format(i20,5d20.6)
         90005 format(i6,3d30.6)
         90006 format(i6,2d30.6)
         90007 format(f7.2,f5.2,7d13.5,f8.5)
         90008 format(2i5,4d15.5)
         90009 format(i3,f11.7,2f13.3,d11.3,f8.3,a5,2f8.3)
         90010 format(f7.2,f5.2,4d13.5,f16.8,i3)
      
      END DO   !end loop in column
      
END SUBROUTINE Matrix

SUBROUTINE spcmasses(aero,gas,spcmass)
   !----------------------------------------------------------------------------------------------------------------------
   !     Routine to calculate the total mass concentration of each model species:
   !     SULF, BCAR, OCAR, DUST, SEAS, NO3, NH4. Aerosol water is not treated.
   !----------------------------------------------------------------------------------------------------------------------
   use memMatrix, only: sulf_map, bcar_map, ocar_map, dust_map, seas_map, &
                   naerobox,ngases,nmass_spcs,gas_h2so4,gas_hno3,gas_nh3, &
		   mass_no3,mass_nh4
   IMPLICIT NONE
   REAL(8), INTENT(IN)                   :: aero(naerobox)
   REAL(8), INTENT(IN)                   :: gas(ngases)
   REAL(8), INTENT(OUT)                  :: spcmass(nmass_spcs+2)

   
   
   spcmass(1) = sum( aero( sulf_map(:) ) ) + gas( gas_h2so4 )
   spcmass(2) = sum( aero( bcar_map(:) ) )
   spcmass(3) = sum( aero( ocar_map(:) ) )
   spcmass(4) = sum( aero( dust_map(:) ) )
   spcmass(5) = sum( aero( seas_map(:) ) )
   spcmass(6) = aero( mass_no3 )           + gas( gas_hno3 )
   spcmass(7) = aero( mass_nh4 )           + gas( gas_nh3  )
   
   
   RETURN
END SUBROUTINE spcmasses
