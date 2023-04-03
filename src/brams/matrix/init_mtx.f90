      module aero_init  
!----------------------------------------------------------------------------------------------------------------------
!
!@sum     defines initial values for the aerosol and gas-phase species for the
!@+     stand-alone version of the matrix microphysical module.
!@auth  susanne bauer/doug wright
!----------------------------------------------------------------------------------------------------------------------
      use memMatrix,  only: conv_dp_to_mass, minconc, densp, aunit1, write_log, &
                            nweights
      use aero_setup
      use aero_discrete, only: discrete_init, ispca, ispcb, ispcc
      implicit none

      contains


      subroutine init_aero( icset, temp, pres )
!----------------------------------------------------------------------------------------------------------------------
!     routine to set initial concentrations for various test cases.
!----------------------------------------------------------------------------------------------------------------------
      implicit none
      integer :: i, index, index1, index2
      integer :: icset           ! identifier of the desired set of initial concentrations
      real(8) :: temp            ! ambient temperature [k]
      real(8) :: pres            ! ambient pressure [pa]
      real(8) :: mpp(nweights)   ! mass per particle in mode i [ug]
      real(8), parameter :: default_n_akk = 1.0d+10  ! 1.0d+10
      real(8), parameter :: default_n_acc = 1.0d+09  ! 1.0d+09
      real(8), parameter :: default_n_dd1 = 1.0d+07  ! 1.0d+07
      real(8), parameter :: default_n_dd2 = 1.0d+06  ! 1.0d+06
      real(8), parameter :: default_n_ssa = 1.0d+08  ! 1.0d+08
      real(8), parameter :: default_n_ssc = 1.0d+05  ! 1.0d+05
      real(8), parameter :: default_n_sss = 1.0d+05  ! 1.0d+05
      real(8), parameter :: default_n_occ = 1.0d+08  ! 1.0d+08
      real(8), parameter :: default_n_bc1 = 1.0d+09  ! 1.0d+08

      ! varibles for the discrete pdf model.

      real(8) :: na,      nb,      nc       ! number concentrations for modes a, b, and c [#/m^3]
      real(8) :: dga,     dgb,     dgc      ! geo. mean diameters   for modes a, b, and c [um]
      real(8) :: sigmaga, sigmagb, sigmagc  ! geo. std. deviations  for modes a, b, and c [1]
      real(8) :: massa,   massb,   massc    ! mass concentrations   for modes a, b, and c [ug/m^3]
      real(8) :: maform,  mbform,  mcform   ! formula mass concentrations for modes a, b, and c [ug/m^3]
      real(8), parameter :: dg_default     = 0.08d+00       ! [um]
      real(8), parameter :: sigmag_default = 1.80d+00       ! [1]


      aero_in(:) = minconc

      !----------------------------------------------------------------------------------------------------------------
      ! use the default mode lognormal parameters and default particle densities for each mode to get a mean mass
      !   per particle to obtain initial mass concentrations from initial number concentrations.
      !
      ! denspi(:) contains the default density for mode i based on its principal chemical component.
      !   parameter conv_dp_to_mass contains the default ambient aerosol density densp, which is
      !   divided out in favor of denspi(:).
      !----------------------------------------------------------------------------------------------------------------
      mpp(:) = denspi(:) * ( conv_dp_to_mass / densp ) &
             * ( 1.0d-06 * dgn0(:) )**3 * exp( 4.5d+00 * ( log( sig0(:) ) )**2 ) ! [ug/particle]

      if ( write_log ) then
        write(aunit1,'(/A/)') '   I    MPP(I) [ug]  <--  Initial mean mass per particle in subr. INIT_AERO'
        do i=1, nweights
          write(aunit1,'(i4,d15.5,f12.5)') i, mpp(i), denspi(i)
        enddo
      endif

      !----------------------------------------------------------------------------------------------------------------
      ! set number concentrations for each mode.
      !----------------------------------------------------------------------------------------------------------------
      select case ( icset )

      case( 0 )

      case( 1 )

        if ( numb_akk_1 .gt. 0 ) aero_in( numb_akk_1 ) = default_n_akk

      case( 2 )

        if ( numb_akk_1 .gt. 0 ) aero_in( numb_akk_1 ) = default_n_akk
        aero_in( numb_acc_1 ) = default_n_acc

      case( 3 )

        if ( numb_akk_1 .gt. 0 ) aero_in( numb_akk_1 ) = default_n_akk
        aero_in( numb_acc_1 ) = default_n_acc
        aero_in( numb_bc1_1 ) = default_n_bc1
!        if ( numb_dd2_1 .gt. 0 ) aero_in( numb_dd2_1 ) = default_n_dd2

      case( 4 )

        if ( numb_akk_1 .gt. 0 ) aero_in( numb_akk_1 ) = default_n_akk
        aero_in( numb_acc_1 ) = default_n_acc
        aero_in( numb_dd1_1 ) = default_n_dd1
        if ( numb_dd2_1 .gt. 0 ) aero_in( numb_dd2_1 ) = default_n_dd2
        if ( numb_ssa_1 .gt. 0 ) aero_in( numb_ssa_1 ) = default_n_ssa
        if ( numb_sss_1 .gt. 0 ) aero_in( numb_sss_1 ) = default_n_sss

      case( 5 )

        if ( numb_akk_1 .gt. 0 ) aero_in( numb_akk_1 ) = default_n_akk
        aero_in( numb_acc_1 ) = default_n_acc
        aero_in( numb_dd1_1 ) = default_n_dd1
        if ( numb_dd2_1 .gt. 0 ) aero_in( numb_dd2_1 ) = default_n_dd2
        if ( numb_ssa_1 .gt. 0 ) aero_in( numb_ssa_1 ) = default_n_ssa
        if ( numb_ssc_1 .gt. 0 ) aero_in( numb_ssc_1 ) = default_n_ssc
        if ( numb_sss_1 .gt. 0 ) aero_in( numb_sss_1 ) = default_n_sss

      case( 6 )

        if ( numb_akk_1 .gt. 0 ) aero_in( numb_akk_1 ) = default_n_akk
        aero_in( numb_acc_1 ) = default_n_acc
        aero_in( numb_dd1_1 ) = default_n_dd1
        if ( numb_dd2_1 .gt. 0 ) aero_in( numb_dd2_1 ) = default_n_dd2
        if ( numb_ssa_1 .gt. 0 ) aero_in( numb_ssa_1 ) = default_n_ssa
        if ( numb_ssc_1 .gt. 0 ) aero_in( numb_ssc_1 ) = default_n_ssc
        if ( numb_sss_1 .gt. 0 ) aero_in( numb_sss_1 ) = default_n_sss
        aero_in( numb_occ_1 ) = default_n_occ

      case( 7 )

        if ( numb_akk_1 .gt. 0 ) aero_in( numb_akk_1 ) = default_n_akk
        aero_in( numb_acc_1 ) = default_n_acc
        aero_in( numb_dd1_1 ) = default_n_dd1
        if ( numb_dd2_1 .gt. 0 ) aero_in( numb_dd2_1 ) = default_n_dd2
        if ( numb_ssa_1 .gt. 0 ) aero_in( numb_ssa_1 ) = default_n_ssa
        if ( numb_ssc_1 .gt. 0 ) aero_in( numb_ssc_1 ) = default_n_ssc
        if ( numb_sss_1 .gt. 0 ) aero_in( numb_sss_1 ) = default_n_sss
        aero_in( numb_occ_1 ) = default_n_occ
        aero_in( numb_bc1_1 ) = default_n_bc1

      case( 8 )

        if ( numb_akk_1 .gt. 0 ) aero_in( numb_akk_1 ) = default_n_akk
        aero_in( numb_acc_1 ) = default_n_acc * 1.0d+01
        aero_in( numb_dd1_1 ) = default_n_dd1 * 1.0d+01
        if ( numb_dd2_1 .gt. 0 ) aero_in( numb_dd2_1 ) = default_n_dd2 * 1.0d+01
        if ( numb_ssa_1 .gt. 0 ) aero_in( numb_ssa_1 ) = default_n_ssa
        if ( numb_ssc_1 .gt. 0 ) aero_in( numb_ssc_1 ) = default_n_ssc
        if ( numb_sss_1 .gt. 0 ) aero_in( numb_sss_1 ) = default_n_sss
        aero_in( numb_occ_1 ) = default_n_occ * 1.0d+01
        aero_in( numb_bc1_1 ) = default_n_bc1 * 1.0d+01

      case( 9 )

        if ( numb_akk_1 .gt. 0 ) aero_in( numb_akk_1 ) = default_n_akk
        aero_in( numb_acc_1 ) = default_n_acc
        aero_in( numb_dd1_1 ) = default_n_dd1
        if ( numb_dd2_1 .gt. 0 ) aero_in( numb_dd2_1 ) = default_n_dd2
        if ( numb_ssa_1 .gt. 0 ) aero_in( numb_ssa_1 ) = default_n_ssa
        if ( numb_ssc_1 .gt. 0 ) aero_in( numb_ssc_1 ) = default_n_ssc
        if ( numb_sss_1 .gt. 0 ) aero_in( numb_sss_1 ) = default_n_sss
        aero_in( numb_occ_1 ) = default_n_occ
        aero_in( numb_bc1_1 ) = default_n_bc1

      !----------------------------------------------------------------------------------------------------------------
      ! case 10 is for testing the discrete model of the pdf. mzj 2005, figure 15.2.
      !----------------------------------------------------------------------------------------------------------------
      case( 10 )

        ispca   = 2
        ispcb   = 2
        ispcc   = 2
        na      = 1.0d+12
        dga     = dg_default
        sigmaga = sigmag_default
        nb      = 0.0d+00
        dgb     = dg_default
        sigmagb = sigmag_default
        nc      = 0.0d+00
        dgc     = dg_default
        sigmagc = sigmag_default
        call discrete_init(icset,na,dga,sigmaga,nb,dgb,sigmagb,nc,dgc,sigmagc,temp,pres,massa,massb,massc)
        aero_in( mass_akk_sulf ) = aero_in( numb_map(1) ) * mpp(1)
        return

      !----------------------------------------------------------------------------------------------------------------
      ! case 11 is for testing the discrete model of the pdf. mzj 2005, figure 15.3.
      !----------------------------------------------------------------------------------------------------------------
      case( 11 )

        ispca   = 1
        ispcb   = 1
        ispcc   = 1
        na      = 1.0d+11
        dga     = dg_default
        sigmaga = sigmag_default
        nb      = 0.0d+00
        dgb     = dg_default
        sigmagb = sigmag_default
        nc      = 0.0d+00
        dgc     = dg_default
        sigmagc = sigmag_default
        call discrete_init(icset,na,dga,sigmaga,nb,dgb,sigmagb,nc,dgc,sigmagc,temp,pres,massa,massb,massc)
        aero_in( mass_akk_sulf ) = aero_in( numb_map(1) ) * mpp(1)
        return

      !----------------------------------------------------------------------------------------------------------------
      ! case 12 is for comparison with results for the discrete model of the pdf. modes akk and acc only.
      !----------------------------------------------------------------------------------------------------------------
      case( 12 )

        aero_in( numb_akk_1 ) = 1.0d+10     ! [#/m^3]
        aero_in( numb_acc_1 ) = 1.0d+09     ! [#/m^3]
        ispca   = 2
        ispcb   = 2
        ispcc   = 2
        index1  = 1
        index2  = 2
        na      = aero_in( numb_akk_1 )
        dga     = dgn0(index1)
        sigmaga = sig0(index1)
        nb      = aero_in( numb_acc_1 )
        dgb     = dgn0(index2)
        sigmagb = sig0(index2)
        nc      = 0.0d+00
        dgc     = dg_default
        sigmagc = sigmag_default
        call discrete_init(icset,na,dga,sigmaga,nb,dgb,sigmagb,nc,dgc,sigmagc,temp,pres,massa,massb,massc)
        maform = aero_in( numb_map(index1) ) * mpp(index1)
        WRITE(36,'(A,2D20.12)')'MASSA from discrete pdf model = ', massa
        WRITE(36,'(A,2D20.12)')'MASSA from analytic formula   = ', maform
        aero_in( mass_akk_sulf ) = massa
        mbform = aero_in( numb_map(index2) ) * mpp(index2)
        WRITE(36,'(A,2D20.12)')'MASSB from discrete pdf model = ', massb
        WRITE(36,'(A,2D20.12)')'MASSB from analytic formula   = ', mbform
        aero_in( mass_acc_sulf ) = massb
        return

      !----------------------------------------------------------------------------------------------------------------
      ! case 13 is for comparison with results for the discrete model of the pdf. acc + bc1 --> bcs.
      !----------------------------------------------------------------------------------------------------------------
      case( 13 )

        aero_in( numb_acc_1 ) = 1.0d+09     ! [#/m^3]
        aero_in( numb_bc1_1 ) = 1.0d+09     ! [#/m^3]
        ispca   = 2
        ispcb   = 2
        ispcc   = 2
        index1  = 2
        index2  = 10
        na      = aero_in( numb_acc_1 )
        dga     = dgn0(index1)
        sigmaga = sig0(index1)
        nb      = aero_in( numb_bc1_1 )
        dgb     = dgn0(index2)
        sigmagb = sig0(index2)
        nc      = 0.0d+00
        dgc     = dg_default
        sigmagc = sigmag_default
        call discrete_init(icset,na,dga,sigmaga,nb,dgb,sigmagb,nc,dgc,sigmagc,temp,pres,massa,massb,massc)
        maform = aero_in( numb_map(index1) ) * mpp(index1)
        WRITE(36,'(A,2D20.12)')'MASSA from discrete pdf model = ', massa
        WRITE(36,'(A,2D20.12)')'MASSA from analytic formula   = ', maform
        aero_in( mass_acc_sulf ) = massa
        mbform = aero_in( numb_map(index2) ) * mpp(index2)
        WRITE(36,'(A,2D20.12)')'MASSB from discrete pdf model = ', massb
        WRITE(36,'(A,2D20.12)')'MASSB from analytic formula   = ', mbform
        AERO_IN( MASS_BC1_BCAR ) = MASSB
        RETURN

      !----------------------------------------------------------------------------------------------------------------
      ! Case 14 is for comparison with results for the discrete model of the PDF. DD2 + OCC --> MXX.
      !----------------------------------------------------------------------------------------------------------------
      case( 14 )

        if ( numb_dd2_1 .gt. 0 ) then
          aero_in( numb_dd2_1 ) = 1.0d+07     ! [#/m^3]
        else
          write(*,*)'cannot use icset = 14 with this mechanism - must have mode dd2 to use this icset.'
          stop
        endif
        aero_in( numb_occ_1 ) = 1.0d+09     ! [#/m^3]
        ispca   = 2
        ispcb   = 2
        ispcc   = 2
        index1  = 5
        index2  = 9
        na      = aero_in( numb_dd2_1 )
        dga     = dgn0(index1)
        sigmaga = sig0(index1)
        nb      = aero_in( numb_occ_1 )
        dgb     = dgn0(index2)
        sigmagb = sig0(index2)
        nc      = 0.0d+00
        dgc     = dg_default
        sigmagc = sigmag_default
        call discrete_init(icset,na,dga,sigmaga,nb,dgb,sigmagb,nc,dgc,sigmagc,temp,pres,massa,massb,massc)
        maform = aero_in( numb_map(index1) ) * mpp(index1)
        WRITE(36,'(A,2D20.12)')'MASSA from discrete pdf model = ', massa
        WRITE(36,'(A,2D20.12)')'MASSA from analytic formula   = ', maform
        aero_in( mass_dd2_dust ) = massa
        mbform = aero_in( numb_map(index2) ) * mpp(index2)
        WRITE(36,'(A,2D20.12)')'MASSB from discrete pdf model = ', massb
        WRITE(36,'(A,2D20.12)')'MASSB from analytic formula   = ', mbform
        aero_in( mass_occ_ocar ) = massb
        RETURN

      !----------------------------------------------------------------------------------------------------------------
      ! Case 15 is for comparison with results for the discrete model of the PDF. AKK + BC1 --> BCS.
      !----------------------------------------------------------------------------------------------------------------
      CASE( 15 )

        aero_in( numb_akk_1 ) = 1.0d+10     ! [#/m^3]
        aero_in( numb_bc1_1 ) = 1.0d+09     ! [#/m^3]
        ispca   = 2
        ispcb   = 2
        ispcc   = 2
        index1  = 1
        index2  = 10
        na      = aero_in( numb_akk_1 )
        dga     = dgn0(index1)
        sigmaga = sig0(index1)
        nb      = aero_in( numb_bc1_1 )
        dgb     = dgn0(index2)
        sigmagb = sig0(index2)
        nc      = 0.0d+00
        dgc     = dg_default
        sigmagc = sigmag_default
        call discrete_init(icset,na,dga,sigmaga,nb,dgb,sigmagb,nc,dgc,sigmagc,temp,pres,massa,massb,massc)
        maform = aero_in( numb_map(index1) ) * mpp(index1)
        WRITE(36,'(A,2D20.12)')'MASSA from discrete pdf model = ', massa
        WRITE(36,'(A,2D20.12)')'MASSA from analytic formula   = ', maform
        aero_in( mass_akk_sulf ) = massa
        mbform = aero_in( numb_map(index2) ) * mpp(index2)
        WRITE(36,'(A,2D20.12)')'MASSB from discrete pdf model = ', massb
        WRITE(36,'(A,2D20.12)')'MASSB from analytic formula   = ', mbform
        aero_in( mass_bc1_bcar ) = massb
        return

      !----------------------------------------------------------------------------------------------------------------
      ! case 16 is for comparison with results for the discrete model of the pdf. dd1 + occ --> mxx.
      !----------------------------------------------------------------------------------------------------------------
      case( 16 )

        aero_in( numb_dd1_1 ) = 1.0d+09     ! [#/m^3]
        aero_in( numb_occ_1 ) = 1.0d+09     ! [#/m^3]
        ispca   = 2
        ispcb   = 2
        ispcc   = 2
        index1  = 3
        index2  = 9
        na      = aero_in( numb_dd1_1 )
        dga     = dgn0(index1)
        sigmaga = sig0(index1)
        nb      = aero_in( numb_occ_1 )
        dgb     = dgn0(index2)
        sigmagb = sig0(index2)
        nc      = 0.0d+00
        dgc     = dg_default
        sigmagc = sigmag_default
        call discrete_init(icset,na,dga,sigmaga,nb,dgb,sigmagb,nc,dgc,sigmagc,temp,pres,massa,massb,massc)
        maform = aero_in( numb_map(index1) ) * mpp(index1)
        WRITE(36,'(A,2D20.12)')'MASSA from discrete pdf model = ', massa
        WRITE(36,'(A,2D20.12)')'MASSA from analytic formula   = ', maform
        aero_in( mass_dd1_dust ) = massa
        mbform = aero_in( numb_map(index2) ) * mpp(index2)
        WRITE(36,'(A,2D20.12)')'MASSB from discrete pdf model = ', massb
        WRITE(36,'(A,2D20.12)')'MASSB from analytic formula   = ', mbform
        aero_in( mass_occ_ocar ) = massb
        RETURN

      !----------------------------------------------------------------------------------------------------------------
      ! Case 17 is for comparison with results for the discrete model of the PDF. DD1 + SSC --> MXX.
      !----------------------------------------------------------------------------------------------------------------
      CASE( 17 )

        aero_in( numb_dd1_1 ) = 1.00d+09     ! [#/m^3]
        aero_in( numb_ssc_1 ) = 0.40d+06     ! [#/m^3]
        ispca   = 2
        ispcb   = 2
        ispcc   = 2
        index1  = 3
        index2  = 8
        na      = aero_in( numb_dd1_1 )
        dga     = dgn0(index1)
        sigmaga = sig0(index1)
        nb      = aero_in( numb_ssc_1 )
        dgb     = dgn0(index2)
        sigmagb = sig0(index2)
        nc      = 0.0d+00
        dgc     = dg_default
        sigmagc = sigmag_default
        call discrete_init(icset,na,dga,sigmaga,nb,dgb,sigmagb,nc,dgc,sigmagc,temp,pres,massa,massb,massc)
        maform = aero_in( numb_map(index1) ) * mpp(index1)
        WRITE(36,'(A,2D20.12)')'MASSA from discrete pdf model = ', massa
        WRITE(36,'(A,2D20.12)')'MASSA from analytic formula   = ', maform
        aero_in( mass_dd1_dust ) = massa
        mbform = aero_in( numb_map(index2) ) * mpp(index2)
        WRITE(36,'(A,2D20.12)')'MASSB from discrete pdf model = ', massb
        WRITE(36,'(A,2D20.12)')'MASSB from analytic formula   = ', mbform
        aero_in( mass_ssc_seas ) = massb
        return

      !----------------------------------------------------------------------------------------------------------------
      ! case 18 is for comparison with the discrete model of the pdf. mode akk. intramodal coagulation only.
      !----------------------------------------------------------------------------------------------------------------
      case( 18 )

        if ( numb_akk_1 .gt. 0 ) aero_in( numb_akk_1 ) = default_n_akk
        ispca   = 2
        ispcb   = 2
        ispcc   = 2
        index1  = 1
        index2  = 2
        na      = aero_in( numb_akk_1 )
        dga     = dgn0(index1)
        sigmaga = sig0(index1)
        nb      = aero_in( numb_acc_1 )
        dgb     = dgn0(index2)
        sigmagb = sig0(index2)
        nc      = 0.0d+00
        dgc     = dg_default
        sigmagc = sigmag_default
        call discrete_init(icset,na,dga,sigmaga,nb,dgb,sigmagb,nc,dgc,sigmagc,temp,pres,massa,massb,massc)
        maform = aero_in( numb_map(index1) ) * mpp(index1)
        WRITE(36,'(A,2D20.12)')'MASSA from discrete pdf model = ', massa
        WRITE(36,'(A,2D20.12)')'MASSA from analytic formula   = ', maform
        aero_in( mass_akk_sulf ) = massa
        mbform = aero_in( numb_map(index2) ) * mpp(index2)
        WRITE(36,'(A,2D20.12)')'MASSB from discrete pdf model = ', massb
        WRITE(36,'(A,2D20.12)')'MASSB from analytic formula   = ', mbform
        aero_in( mass_acc_sulf ) = massb
        return

      !----------------------------------------------------------------------------------------------------------------
      ! case 19 is for comparison with results for the discrete model of the pdf. bc1 + occ --> boc.
      !----------------------------------------------------------------------------------------------------------------
      case( 19 )

        aero_in( numb_occ_1 ) = 1.0d+09     ! [#/m^3]
        aero_in( numb_bc1_1 ) = 1.0d+09     ! [#/m^3]
        ispca   = 2
        ispcb   = 2
        ispcc   = 2
        index1  = 9
        index2  = 10
        na      = aero_in( numb_occ_1 )
        dga     = dgn0(index1)
        sigmaga = sig0(index1)
        nb      = aero_in( numb_bc1_1 )
        dgb     = dgn0(index2)
        sigmagb = sig0(index2)
        nc      = 0.0d+00
        dgc     = dg_default
        sigmagc = sigmag_default
        call discrete_init(icset,na,dga,sigmaga,nb,dgb,sigmagb,nc,dgc,sigmagc,temp,pres,massa,massb,massc)
        maform = aero_in( numb_map(index1) ) * mpp(index1)
        WRITE(36,'(A,2D20.12)')'MASSA from discrete pdf model = ', massa
        WRITE(36,'(A,2D20.12)')'MASSA from analytic formula   = ', maform
        aero_in( mass_occ_ocar ) = massa
        mbform = aero_in( numb_map(index2) ) * mpp(index2)
        WRITE(36,'(A,2D20.12)')'MASSB from discrete pdf model = ', massb
        WRITE(36,'(A,2D20.12)')'MASSB from analytic formula   = ', mbform
        aero_in( mass_bc1_bcar ) = massb
        return

      !----------------------------------------------------------------------------------------------------------------
      ! case 20 is for comparison with results for the discrete model of the pdf. acc + dd1 --> dd1. (no ds1 transfer)
      !----------------------------------------------------------------------------------------------------------------
      case( 20 )

        aero_in( numb_acc_1 ) = 1.0d+09     ! [#/m^3]
        aero_in( numb_dd1_1 ) = 1.0d+09     ! [#/m^3]
        ispca   = 2
        ispcb   = 2
        ispcc   = 2
        index1  = 2
        index2  = 3
        na      = aero_in( numb_acc_1 )
        dga     = dgn0(index1)
        sigmaga = sig0(index1)
        nb      = aero_in( numb_dd1_1 )
        dgb     = dgn0(index2)
        sigmagb = sig0(index2)
        nc      = 0.0d+00
        dgc     = dg_default
        sigmagc = sigmag_default
        call discrete_init(icset,na,dga,sigmaga,nb,dgb,sigmagb,nc,dgc,sigmagc,temp,pres,massa,massb,massc)
        maform = aero_in( numb_map(index1) ) * mpp(index1)
        WRITE(36,'(A,2D20.12)')'MASSA from discrete pdf model = ', massa
        WRITE(36,'(A,2D20.12)')'MASSA from analytic formula   = ', maform
        aero_in( mass_acc_sulf ) = massa
        mbform = aero_in( numb_map(index2) ) * mpp(index2)
        WRITE(36,'(A,2D20.12)')'MASSB from discrete pdf model = ', massb
        WRITE(36,'(A,2D20.12)')'MASSB from analytic formula   = ', mbform
        aero_in( mass_dd1_dust ) = massb
        return

      !----------------------------------------------------------------------------------------------------------------
      ! case 21 is for illustration of aerosol activation.
      !----------------------------------------------------------------------------------------------------------------
      case( 21 )

        if ( numb_akk_1 .gt. 0 ) then
          aero_in( numb_akk_1 ) = 1.0d+09 * 2.0d+00
          aero_in( numb_acc_1 ) = 1.0d+08 * 2.0d+00
        else
          aero_in( numb_acc_1 ) = 1.1d+09 * 2.0d+00
        endif
        if ( numb_dd2_1 .gt. 0 ) then
          aero_in( numb_dd1_1 ) = 1.0d+08 * 0.1d+00
          aero_in( numb_dd2_1 ) = 1.0d+07 * 0.1d+00
        else
          aero_in( numb_dd1_1 ) = 1.1d+08 * 0.1d+00
        endif
        if ( numb_ssa_1 .gt. 0 ) aero_in( numb_ssa_1 ) = 1.000d+08 * 0.2d+00
        if ( numb_ssc_1 .gt. 0 ) aero_in( numb_ssc_1 ) = 0.004d+08 * 0.2d+00
        if ( numb_sss_1 .gt. 0 ) aero_in( numb_sss_1 ) = 1.004d+08 * 0.2d+00
        aero_in( numb_occ_1 ) = 1.0d+08 * 1.0d+01
        aero_in( numb_bc1_1 ) = 1.0d+08 * 1.0d+01

      case default

        WRITE(*,*)'BAD VALUE OF ICSET IN SUBR. INIT_AERO'
        STOP

      END SELECT

      !----------------------------------------------------------------------------------------------------------------
      ! Calculate masses for all primary modes.
      !----------------------------------------------------------------------------------------------------------------
      do i=1, nmodes

        ! set the defining chemical species for each mode to receive the initial concentrations.

        index = 0
        if( mode_name(i) .EQ. 'AKK' ) index = mass_akk_sulf
        if( mode_name(i) .EQ. 'ACC' ) index = mass_acc_sulf
        if( mode_name(i) .EQ. 'DD1' ) index = mass_dd1_dust
        if( mode_name(i) .EQ. 'DD2' ) index = mass_dd2_dust
        if( mode_name(i) .EQ. 'SSA' ) index = mass_ssa_seas
        if( mode_name(i) .EQ. 'SSC' ) index = mass_ssc_seas
        if( mode_name(i) .EQ. 'SSS' ) index = mass_sss_seas
        if( mode_name(i) .EQ. 'OCC' ) index = mass_occ_ocar
        if( mode_name(i) .EQ. 'BC1' ) index = mass_bc1_bcar

        if( index .gt. 0 ) then

          ! write(*,*)'index = ', index
          aero_in( index ) = aero_in( numb_map(i) ) * mpp(i)

          ! for some ic sets, initial sulfate concentrations are also given to these modes.

          if( icset.eq.9 .and. mode_name(i).ne.'ssc' ) aero_in( sulf_map(i) ) = aero_in( numb_map(i) ) * mpp(i)

        endif
      enddo

      ! set ammonium assuming ammonium sulfate.

      if( icset .eq. 21 ) aero_in( 2 ) = 0.3755532793d+00 * sum( aero_in(sulf_map(:)) )

      return
      end subroutine Init_Aero


      END MODULE Aero_Init
