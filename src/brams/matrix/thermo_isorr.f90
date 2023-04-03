      subroutine Aero_Thermo(aso4,ano3,anh4,ah2o,gnh3,ghno3,tot_dust, &
                            seas,ssh2o,tk,rh,pres,rhd,rhc)
!@sum
!@+     this routine sets up for and calls the thermodynamic module for aerosol
!@+     gas-particle partitioning.
!@+
!@+     this version of aero_thermo is for use with the isorropia thermodynamic module.
!@auth susanne bauer/doug wright

      use memMatrix, only: write_log, aunit1
      use Isorropia_Module, ONLY: isoropia
      IMPLICIT NONE

      !------------------------------------------------------------------------------------------------------
      ! arguments.
      !------------------------------------------------------------------------------------------------------
      real(8), intent(inout) :: aso4      ! aerosol sulfate       [ug/m^3]
      real(8), intent(inout) :: ano3      ! aerosol nitrate       [ug/m^3]
      real(8), intent(inout) :: anh4      ! aerosol ammonium      [ug/m^3]
      real(8), intent(inout) :: ah2o      ! aerosol water         [ug/m^3]
      real(8), intent(inout) :: gnh3      ! gas-phase ammonia     [ugnh4/m^3] as ammonium (mw)
      real(8), intent(inout) :: ghno3     ! gas-phase nitric acid [ugno3/m^3] as nitrate  (mw)
      real(8), intent(in)    :: tot_dust  ! total dust(sol+insol) [ug/m^3]
      real(8), intent(in)    :: seas      ! sea salt (nacl)       [ug/m^3]
      real(8), intent(out)   :: ssh2o     ! sea salt assoc. h2o   [ug/m^3]
      real(8), intent(in)    :: tk        ! absolute temperature  [k]          
      real(8), intent(in)    :: rh        ! relative humidity     [0-1]
      real(8), intent(in)    :: pres      ! ambient pressure      [pa]  
      real(8), intent(out)   :: rhd       ! rh of deliquescence   [0-1]
      real(8), intent(out)   :: rhc       ! rh of crystallization [0-1]

      !------------------------------------------------------------------------------------------------------
      ! input to isoropia.
      !------------------------------------------------------------------------------------------------------
      real(8) :: wi(5)        ! [moles/m^3]
      real(8) :: rhi          ! [0.0-1.0]
      real(8) :: tempi        ! [k]
      real(8) :: cntrl(2)     ! [1] control variables

      !------------------------------------------------------------------------------------------------------
      ! output from isoropia.
      !------------------------------------------------------------------------------------------------------
      real(8) :: wt(5)        ! [moles/m^3]
      real(8) :: gas(3)       ! [moles/m^3]
      real(8) :: aerliq(12)   ! [moles/m^3]
      real(8) :: aersld(9)    ! [moles/m^3]
      real(8) :: other(6)     ! 
      character(len=15) :: scasi = '               '
   
      !------------------------------------------------------------------------------------------------------
      ! parameters. double-precision molecular weights [g/mol] and their reciprocals.
      !------------------------------------------------------------------------------------------------------
      real(8), parameter :: mw_anh4   = 18.03850d+00  ! [g/mol]
      real(8), parameter :: mw_gnh3   = mw_anh4       ! [g/mol] nh3  is passed as equivalent conc. of nh4+
      real(8), parameter :: mw_ano3   = 62.00494d+00  ! [g/mol]
      real(8), parameter :: mw_ghno3  = mw_ano3       ! [g/mol] hno3 is passed as equivalent conc. of no3-
      real(8), parameter :: mw_aso4   = 96.0636d+00   ! [g/mol]
      real(8), parameter :: mw_na     = 22.989768d+00 ! [g/mol]
      real(8), parameter :: mw_cl     = 35.4527d+00   ! [g/mol]
      real(8), parameter :: mw_nacl   = 58.442468d+00 ! [g/mol]
      real(8), parameter :: mw_h2o    = 18.01528d+00  ! [g/mol]
      real(8), parameter :: rmw_na    = 1.0d-06 / mw_na    ! [mol/g]
      real(8), parameter :: rmw_aso4  = 1.0d-06 / mw_aso4  ! [mol/g]
      real(8), parameter :: rmw_anh4  = 1.0d-06 / mw_anh4  ! [mol/g]
      real(8), parameter :: rmw_gnh3  = 1.0d-06 / mw_gnh3  ! [mol/g]
      real(8), parameter :: rmw_ano3  = 1.0d-06 / mw_ano3  ! [mol/g]
      real(8), parameter :: rmw_ghno3 = 1.0d-06 / mw_ghno3 ! [mol/g]
      real(8), parameter :: rmw_cl    = 1.0d-06 / mw_cl    ! [mol/g]
      real(8), parameter :: cmw_na    = 1.0d+06 * mw_na    ! [ug/mol]
      real(8), parameter :: cmw_aso4  = 1.0d+06 * mw_aso4  ! [ug/mol]
      real(8), parameter :: cmw_anh4  = 1.0d+06 * mw_anh4  ! [ug/mol]
      real(8), parameter :: cmw_gnh3  = 1.0d+06 * mw_gnh3  ! [ug/mol]
      real(8), parameter :: cmw_ano3  = 1.0d+06 * mw_ano3  ! [ug/mol]
      real(8), parameter :: cmw_ghno3 = 1.0d+06 * mw_ghno3 ! [ug/mol]
      real(8), parameter :: cmw_cl    = 1.0d+06 * mw_cl    ! [ug/mol]
      real(8), parameter :: cmw_h2o   = 1.0d+06 * mw_h2o   ! [g/mol]

      !------------------------------------------------------------------------------------------------------
      ! fraction of sea salt (nacl) mass that is na, and is cl.
      !------------------------------------------------------------------------------------------------------
      real(8), parameter :: rat_na = mw_na / ( mw_na + mw_cl )  ! [1] 
      real(8), parameter :: rat_cl = mw_cl / ( mw_na + mw_cl )  ! [1] 

      !------------------------------------------------------------------------------------------------------
      ! other parameters.
      !------------------------------------------------------------------------------------------------------
      real(8), parameter :: dh2o   = 1.000d+00   ! density of water [g/cm^3]
      real(8), parameter :: dnacl  = 2.165d+00   ! density of nacl  [g/cm^3]
      real(8), parameter :: css    = 1.08d+00    ! for sea salt ...
      real(8), parameter :: bss    = 1.2d+00     ! for sea salt ...              
      real(8), parameter :: ssh2oa = (css*css*css*bss-1.0d+00)*dh2o/dnacl
      real(8), parameter :: ssh2ob = (css*css*css            )*dh2o/dnacl
      real(8), parameter :: rhmax  = 0.995d+00   ! [0-1]
      real(8), parameter :: rhmin  = 0.010d+00   ! [0-1]   
      real(8)            :: h                    ! local rh, with rhmin < h < rhmax

      !------------------------------------------------------------------------------------------------------
      ! call for the bulk non-sea salt inorganic aerosol.
      !------------------------------------------------------------------------------------------------------
      if ( write_log ) then
        WRITE(aunit1,'(/A,3F12.3/)') 'ISORROPIA: TK[K], RH[0-1], PRES[Pa] = ', tk, rh, pres
        WRITE(aunit1,'(A4,7A14  )') '   ','ASO4','ANO3','ANH4','AH2O', 'GNH3','GHNO3','TOT_DUST'
        WRITE(aunit1,'(A4,7E14.5)') 'TOP',aso4,ano3,anh4,ah2o,gnh3,ghno3,tot_dust
      ENDIF

      h = max( min( rh, rhmax ), rhmin )
!     wi(1) = rat_na*seas*rmw_na                       ! from [ug/m^3] to [mol/m^3]
      wi(1) = 0.0d+00                                  ! na
      wi(2) =        aso4*rmw_aso4                     ! from [ug/m^3] to [mol/m^3]
      wi(3) =        anh4*rmw_anh4 +  gnh3*rmw_gnh3    ! from [ug/m^3] to [mol/m^3]
      wi(4) =        ano3*rmw_ano3 + ghno3*rmw_ghno3   ! from [ug/m^3] to [mol/m^3]
!     wi(5) = rat_cl*seas*rmw_cl                       ! from [ug/m^3] to [mol/m^3]
      wi(5) = 0.0d+00                                  ! cl
      cntrl(1) = 0.0d+00  ! forward problem: wi contains the gas+aerosol concentrations
      cntrl(2) = 1.0d+00  ! 0 (solid & liquid phases), 1 (liquid only, metastable)

      wt(:)     = 0.0d+00
      gas(:)    = 0.0d+00
      aerliq(:) = 0.0d+00
      aersld(:) = 0.0d+00
      other(:)  = 0.0d+00

      call isoropia ( wi, h, tk, cntrl, wt, gas, aerliq, aersld, scasi, other )

      gnh3  = max( gas(1)*cmw_gnh3,  0.0d+00 )    ! from [mol/m^3] to [ug/m^3]
      ghno3 = max( gas(2)*cmw_ghno3, 0.0d+00 )    ! from [mol/m^3] to [ug/m^3]
      aso4  = wt(2)*cmw_aso4                      ! from [mol/m^3] to [ug/m^3]
      anh4  = wt(3)*cmw_anh4 - gnh3               ! from [mol/m^3] to [ug/m^3]
      ano3  = wt(4)*cmw_ano3 - ghno3              ! from [mol/m^3] to [ug/m^3]
      ah2o  = aerliq(8)*cmw_h2o                   ! from [mol/m^3] to [ug/m^3]
      anh4  = max( anh4, 0.0d+00 )                ! [ug/m^3]
      ano3  = max( ano3, 0.0d+00 )                ! [ug/m^3]

      rhd   = 0.80d+00                            ! rhd = 0.80 for ammonium sulfate (ghan et al., 2001).
      rhc   = 0.35d+00                            ! rhc = 0.35 for ammonium sulfate (ghan et al., 2001).

      if ( write_log ) then
        write(aunit1,'(A4,7E14.5)') 'END',aso4,ano3,anh4,ah2o,gnh3,ghno3,tot_dust
        WRITE(aunit1,'(A4,7F14.5)') 'RHD',rhd
      ENDIF 

      !-------------------------------------------------------------------------
      ! Get the sea salt-associated water (only).
      !
      ! A simple parameterization provided by E. Lewis is used.
      !-------------------------------------------------------------------------
      if ( write_log ) then
        write(aunit1,'(A4,3A12  )') '   ','SEAS','SSH2O'           
        WRITE(aunit1,'(A4,3F12.5)') 'TOP' ,seas
      ENDIF

      if ( h .gt. 0.45d+00 ) then     ! ... then we are above the crystallization rh of nacl
        ssh2o = seas * ( ssh2oa + ssh2ob / ( 1.0d+00 - h ) )
      else
        ssh2o = 0.0d+00
      endif

      if ( write_log ) then
        write(aunit1,'(A4,3F12.5)') 'END',seas,ssh2o
      ENDIF

      RETURN
      END SUBROUTINE Aero_Thermo


