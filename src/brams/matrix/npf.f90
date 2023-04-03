!>-------------------------------------------------------------------------------------------------------------------
!!
!!@brief     this module contains all sub-programs to calculate nucleation and 
!!           new particle formation rates.
!!@author  susanne bauer/doug wright
!-------------------------------------------------------------------------------------------------------------------
module Aero_Npf
      use memMatrix, ONLY: &
                           pi6,                         & !intent()
                           ugm3_ncm3,                   & !intent()
                           ilay,                        & !intent()
                           zheight,                     & !intent()
                           write_log,                   & !intent()
                           aunit1,                      & !intent()
                           rho_nh42so4,                 & !intent()
                           mw_so4,                      & !intent()
                           mw_nh42so4,                  & !intent()
                           rho_h2so4,                   & !intent()
                           mw_h2so4,                    & !intent()
                           mw_nh3,                      & !intent()
                           mw_h2o,                      & !intent()
                           pi,                          & !intent()
                           avo, &                         !intent()     
                           diffcoef_m2s,                & !intent()
                           avg_dp_of_avg_mass_meters      !intent()
      IMPLICIT NONE
      !-------------------------------------------------------------------------------------------------------------
      ! Select the nucleation scheme: inuc = 1, JVM 
      !                               inuc = 2, VEHKAMAKI
      !                               inuc = 3, NAPARI
      !                               inuc = 4, EISELE AND MCMURRY, 1997
      !-------------------------------------------------------------------------------------------------------------
      integer, parameter :: inuc = 3
      logical, parameter :: include_ion_ion = .true.      ! include turco scheme

      !-------------------------------------------------------------------------------------------------------------
      ! flag for use of the kerminem and kulmala (2002) parameterization for
      ! conversion of a nucleation rate to a particle formation rate at a
      ! larger user-selected size.
      !-------------------------------------------------------------------------------------------------------------
      logical, parameter :: kk02 = .true.
      logical, parameter :: write_f_kk02 = .false.

      !-------------------------------------------------------------------------------------------------------------
      ! new particle parameters.
      !-------------------------------------------------------------------------------------------------------------
      real(8), parameter :: dstar_nm =  1.0d+00   ! diameter of critical nucleus [nm]
      real(8), parameter :: dnpf_nm  =  3.0d+00   ! diameter of e&m(1997) new particles [nm]
      real(8), parameter :: dnu_nm   =  3.0d+00   ! diameter of a new particle [nm]
      real(8), parameter :: rnu_nm = dnu_nm*0.5d+00 ! radius of a new particle [m]
      real(8), parameter :: dnu = dnu_nm*1.0d-09  ! diameter of a new particle [m]
      real(8), parameter :: vnu = pi6*dnu*dnu*dnu ! volume of a new particle   [m^3]
      !-------------------------------------------------------------------------------------------------------------
      ! for conversion from 1-nm particles to dnu_nm-nm particles.
      !-------------------------------------------------------------------------------------------------------------
      real(8), parameter :: conv_nucl_to_npf = 1.0d+06  &  ! [#/cm^3/s] to [#/m^3/s]
               * dstar_nm * dstar_nm * dstar_nm / ( dnu_nm * dnu_nm * dnu_nm )
      !-------------------------------------------------------------------------------------------------------------
      ! for conversion from dnpf_nm-nm particles to dnu_nm-nm particles.
      !-------------------------------------------------------------------------------------------------------------
      real(8), parameter :: conv_emnpf_to_npf = 1.0d+06 &  ! [#/cm^3/s] to [#/m^3/s]
               * dnpf_nm  * dnpf_nm  * dnpf_nm  / ( dnu_nm * dnu_nm * dnu_nm )

      !-------------------------------------------------------------------------------------------------------------
      ! limits on the nucleation rate applied to all parameterizations,
      ! except that the upper limit is not applied to the vehkamaki et al. 2002
      ! parameterization for binary h2so4-h2o nucleation.
      !-------------------------------------------------------------------------------------------------------------
      real(8), parameter :: j_lower = 1.0d-07    ! [#/cm^3/s]
      real(8), parameter :: j_upper = 1.0d+07    ! [#/cm^3/s]

      !-------------------------------------------------------------------------------------------------------------
      ! npfmass(i,n) is the sulfate mass [ugso4] per new particle of volume vnu for i=nint[rh(%)].
      !
      ! particle composition is presently divided into three regimes:
      !   (1) ammonium sulfate   - second index set to 2
      !   (2) ammonium bisulfate - second index set to 1
      !   (3) sulfuric acid      - second index set to 0
      !-------------------------------------------------------------------------------------------------------------
      real(8), save :: npfmass(0:100,0:2)   ! [ugso4/particle]
      integer, save :: npfmass_regime = 2   ! second index to npfmass
      contains
 

!>-------------------------------------------------------------------------------------------------------------------
!!     dlw 2006.           
!!     routine to calculate the rate of production of new particles and the 
!!     corresponding rate of production of particulate sulfate mass.
!!-------------------------------------------------------------------------------------------------------------------
      subroutine NpfRate(prs,rh,temp,xh2so4,so4rate,xnh3,kc,dndt,dmdt_so4,icall)
      implicit none

      ! input arguments.

      real(8), intent(in)   :: prs       !< pressure [pa]
      real(8), intent(in)   :: rh        !< fractional relative humidity [1]
      real(8), intent(in)   :: temp      !< ambient temperature [k]
      real(8), intent(in)   :: xh2so4    !< sulfuric acid (as so4) concentration [ugso4/m^3]
      real(8), intent(in)   :: so4rate   !< gas-phase h2so4 (as so4) production rate [ugso4/m^3 s]
      real(8), intent(in)   :: xnh3      !< ammonia mixing ratio [ppmv]
      real(8), intent(in)   :: kc        !< condensational sink [1/s]
      integer, intent(in)   :: icall     !< flag signalling type of call

      ! output arguments.

      real(8), intent(out)  :: dndt      ! particle number production rate [m^-3 s^-1]
      real(8), intent(out)  :: dmdt_so4  ! so4 mass production rate        [ugso4/m^3 s]

      ! scratch local variables.

      real(8) :: jtot        ! total nucleation rate [cm^-3 s^-1]
      real(8) :: jgas        ! homogeneous nucleation rate [cm^-3 s^-1]
      real(8) :: jion        ! ion-ion recombination nucleation rate [cm^-3 s^-1]
      real(8) :: so4mass     ! mass of so4 per new particle [ugso4]
      real(8) :: h2so4_tmp   ! [h2so4] scratch variable [molecules/cm^3]
      real(8) :: nh3_tmp     ! [nh3]   scratch variable [ppt]                 
      real(8) :: jp_to_j     ! ratio of the particle formation rate to the nucleation rate [1]

      !-------------------------------------------------------------------------------------------------------------
      ! in the call to nucleation_rate, rh is converted to [%], the sulfuric
      ! acid (so4) concentration is converted from ugso4/m^3 to molecules/cm^3,
      ! and ammonia is converted from ppm to ppt.
      !-------------------------------------------------------------------------------------------------------------
      h2so4_tmp = ugm3_ncm3 * max ( xh2so4, 1.0d-30 )     ! [molecule/cm^3]
      nh3_tmp   =   1.0d+06 * max ( xnh3  , 1.0d-30 )     ! [ppt]

      call Nucleation_Rate( temp, 1.0d+02*rh, h2so4_tmp, nh3_tmp, &
                            zheight(ilay), jtot, jgas, jion )

      ! write(78,*)'temp, 100.0d+00*rh, h2so4_tmp, nh3_tmp, zheight(ilay), jtot, jgas, jion'
      ! write(78,*) temp, 100.0d+00*rh, h2so4_tmp, nh3_tmp, zheight(ilay), jtot, jgas, jion  

      !-------------------------------------------------------------------------------------------------------------
      ! convert the nucleation rate in [#/cm^3/s] into a new particle formation rate for 
      !   particles of diameter dnu_nm in [#/m^3/s].
      !
      ! there are two options, either the kerminen and kulmala (2002) parameterization, 
      !   simply reduce the nucleation rate based upon mass conservation as 
      !   1-nm diameter particles to converted to particles at the selected size. 
      !
      ! if the eisele and mcmurry (1997) curves are used, the conversion is
      !   from 3.0 nm (rather than 1 nm) to the selected size. 
      !-------------------------------------------------------------------------------------------------------------
      if( inuc.eq.1 .or. inuc.eq.2 .or. inuc.eq.3 ) then
        if( kk02 ) then
          call F_Kk02( prs, temp, 1.0d+02*rh, h2so4_tmp, nh3_tmp, kc, dnu_nm, dstar_nm, jp_to_j )
          dndt = ( 1.0d+06 * jtot ) * jp_to_j ! convert [#/cm^3/s] to [#/m^3/s].
        else
          dndt = conv_nucl_to_npf * jtot      ! conv_nucl_to_npf includes the [#/cm^3/s] to [#/m^3/s] conversion. 
        endif
      elseif( inuc.eq.4 ) then    ! eisele and mcmurry (1997) curves are used.
        if( kk02 ) then
          call F_Kk02( prs, temp, 1.0d+02*rh, h2so4_tmp, nh3_tmp, kc, dnu_nm, dnpf_nm, jp_to_j )
          dndt = ( 1.0d+06 * jtot ) * jp_to_j ! convert [#/cm^3/s] to [#/m^3/s].
        else
          dndt = conv_emnpf_to_npf * jtot     ! conv_emnpf_to_npf includes the [#/cm^3/s] to [#/m^3/s] conversion.
        endif
      endif
      ! write(34,'(a,4d13.5)')'in npfrate: jtot, jp_to_j, dndt, h2so4 = ', jtot, jp_to_j, dndt, h2so4_tmp

      !-------------------------------------------------------------------------------------------------------------
      ! calculate mass production rate [ugso4/m^3/s)], limited by the
      ! production rate by the production rate of h2so4.
      ! adjust the number production rate if necessary.
      !
      ! npfmass(i,n) is the sulfate mass [ugso4] per new particle of volume vnu for i=nint[rh(%)].
      !
      ! particle composition is presently divided into three regimes:--> 
      !   (1) ammonium sulfate   --> second index set to 2 -->  npfmass_regime = 2
      !   (2) ammonium bisulfate --> second index set to 1 -->  npfmass_regime = 1
      !   (3) sulfuric acid      --> second index set to 0 -->  npfmass_regime = 0
      !-------------------------------------------------------------------------------------------------------------
      so4mass = npfmass( nint( 100.0d+00*rh ), npfmass_regime )   ! [ugso4]
      dmdt_so4 = so4mass * dndt                                   ! [ugso4/m^3/s]
      if( icall .gt. 0 ) return                                   ! do not impose mass limitation for this call
      if ( dmdt_so4 .gt. so4rate ) then                           ! [ugso4/m^3/s]
        ! if ( dmdt_so4 .gt. 1.0d-10 ) write(34,*)'npfrate mass-limit imposed: dmdt_so4, so4rate=',dmdt_so4,so4rate
        dmdt_so4 = so4rate                                        ! [ugso4/m^3/s]
        dndt = dmdt_so4 / so4mass                                 ! [   # /m^3/s]
      endif

      return
      end subroutine NpfRate


!>-------------------------------------------------------------------------------------------------------------------
!!     dlw 2006.
!!     routine to pre-calculate the so4 mass in a single new particle as
!!     a function of rh and new particle volume, the array npfmass(0:100,0:2).
!!-------------------------------------------------------------------------------------------------------------------
      subroutine Setup_NpfMass
      implicit none
      integer :: i 
      real(8) :: h, xi_ke_sulfate, xi_ke_h2so4
      !-------------------------------------------------------------------------------------------------------------
      ! conv_v_to_so4_xxx is the mass of sulfate (mw=96g/mol) in [ugso4] for either xxx=sulfate or xxx=h2so4.
      !   in a dry particle of volume vnu [m^3]. checked on 7-27-06.
      !   in conv_v_to_so4, the 1.0d+12 converts [gso4/cm^3] to [ugso4/m^3].
      !
      ! rhc_nh42so4_rnu_nm is the crystallization rh for a 5-nm (radius) dry particle (e. lewis).
      !-------------------------------------------------------------------------------------------------------------
      real(8), parameter :: conv_v_to_so4_sulfate = rho_nh42so4 * 1.0d+12 * mw_so4 / mw_nh42so4
      real(8), parameter :: conv_v_to_so4_h2so4   = rho_h2so4   * 1.0d+12 * mw_so4 / mw_h2so4
      real(8), parameter :: rhc_nh42so4_dnu_03nm  = 0.00d+00
      real(8), parameter :: rhc_nh42so4_dnu_10nm  = 0.50d+00
      real(8), parameter :: rhc_nh42so4_dnu_20nm  = 0.50d+00
      real(8), parameter :: rhc_h2so4_dnu_03nm    = 0.10d+00
      real(8), parameter :: rhc_h2so4_dnu_10nm    = 0.05d+00
      real(8), parameter :: rhc_h2so4_dnu_20nm    = 0.05d+00

      !-------------------------------------------------------------------------------------------------------------
      ! npfmass(i,j) is the sulfate mass [ugso4] per new particle of volume vnu for i=nint[rh(%)].
      !
      ! particle composition is presently divided into three regimes:
      !   (1) ammonium sulfate   - second index j = 2
      !   (2) ammonium bisulfate - second index j = 1
      !   (3) sulfuric acid      - second index j = 0
      ! 
      ! the kelvin effect is taken into account.
      !
      ! xi_ke is the radius ratio r_ambient/r_dry, including the kelvin effect, and is used to convert
      !   ambient particle volume to dry particle volume.
      !
      ! conv_v_to_so4 converts dry volume ammonium sulfate [m^3] to dry mass sulfate (mw=96g/mol) [ugso4].
      !-------------------------------------------------------------------------------------------------------------
      if( write_log ) write(aunit1,90)
      do i=0, 100
        h = min( 1.0d-02 * dble(i), 0.999d+00 )
        !-----------------------------------------------------------------------------------------------------------
        ! ammonium sulfate. 
        ! ammonium bisulfate. the sulfate value is also used because of lack of data for bisulfate solutions.
        !-----------------------------------------------------------------------------------------------------------
        if    ( dnu_nm .eq.  3.0d+00 ) then
          if ( h .ge. rhc_nh42so4_dnu_03nm ) then     ! the particle is wet.
            xi_ke_sulfate = 1.0d+00 + 0.2d+00*h       ! linear fit over the entire rh range.
          else
            xi_ke_sulfate = 1.0d+00 + 0.2d+00*h       ! linear fit over the entire rh range.
          endif
        elseif( dnu_nm .eq. 10.0d+00 ) then
          if ( h .ge. rhc_nh42so4_dnu_10nm ) then     ! the particle is wet.
            xi_ke_sulfate = 0.677d+00 + h*(1.816d+00 + h*( -2.345d+00 + h*1.296d+00 ) )
          else                                        ! the particle is dry. 
            xi_ke_sulfate = 1.0d+00 + 0.16075d+00*(h/rhc_nh42so4_dnu_10nm)
          endif
        elseif( dnu_nm .eq. 20.0d+00 ) then
          if ( h .ge. rhc_nh42so4_dnu_20nm ) then     ! the particle is wet.
            xi_ke_sulfate = 0.175d+00 + h*(4.532d+00 + h*( -6.894d+00 + h*3.856d+00 ) )
          else                                        ! the particle is dry. 
            xi_ke_sulfate = 1.0d+00 + 0.1995d+00*(h/rhc_nh42so4_dnu_20nm)
          endif
        else
          WRITE(*,*)'Bad value of DNU_NM in subr. SETUP_NPFMASS: DNU_NM = ', DNU_NM
          stop
        endif
        npfmass(i,2) = ( vnu / xi_ke_sulfate**3 ) * conv_v_to_so4_sulfate
        npfmass(i,1) = npfmass(i,2)
        !-----------------------------------------------------------------------------------------------------------
        ! sulfuric acid.
        !-----------------------------------------------------------------------------------------------------------
        if    ( dnu_nm .eq.  3.0d+00 ) then
          if ( h .ge. rhc_h2so4_dnu_03nm ) then                 
            xi_ke_h2so4 = 1.14d+00 + h*(0.464d+00 + h*( -0.336d+00 + h*0.189d+00 ) )
          else
            xi_ke_h2so4 = 1.0d+00 + 0.179d+00*(h/rhc_h2so4_dnu_03nm) ! linear fit over this range.
          endif
        elseif( dnu_nm .eq. 10.0d+00 ) then
          if ( h .ge. rhc_h2so4_dnu_10nm ) then                
            xi_ke_h2so4 = 1.14d+00 + h*(0.765d+00 + h*( -0.850d+00 + h*0.745d+00 ) )
          else                                    
            xi_ke_h2so4 = 1.0d+00 + 0.211d+00*(h/rhc_h2so4_dnu_10nm) ! linear fit over this range.
          endif
        elseif( dnu_nm .eq. 20.0d+00 ) then
          if ( h .ge. rhc_h2so4_dnu_20nm ) then                
            xi_ke_h2so4 = 1.113d+00 + h*(1.190d+00 + h*( -2.001d+00 + h*1.750d+00 ) )
          else                                    
            xi_ke_h2so4 = 1.0d+00 + 0.168d+00*(h/rhc_h2so4_dnu_20nm) ! linear fit over this range.
          endif
        else
          write(*,*)'Bad value of DNU_NM in subr. SETUP_NPFMASS: DNU_NM = ', DNU_NM
          stop
        endif
        npfmass(i,0) = ( vnu / xi_ke_h2so4**3 ) * conv_v_to_so4_h2so4
        ! if( write_log ) write(aunit1,'(i6,2f25.6,3d15.6)') i, xi_ke_h2so4, xi_ke_sulfate, npfmass(i,0:2)
      enddo

90    FORMAT(/,' RH[%]','   R_ambient/R_dry[H2SO4]',' R_ambient/R_dry[SULFATE]', &
               '   NPFMASS(I,0:2)[ugSO4]')
      return
      end subroutine Setup_NpfMass


!>-----------------------------------------------------------------------------------------------------------------------
!!     lsc/dlw 2005-2006:
!!-----------------------------------------------------------------------------------------------------------------------
      subroutine Nucleation_Rate (t,rh,na,mb,z,jtot,jgas,jion)
      implicit none

      ! input arguments.

      real(8), intent( in ) :: t      !< temperature [k]
      real(8), intent( in ) :: rh     !< relative humidity [%]
      real(8), intent( in ) :: na     !< h2so4 concentration [molecules/cm^3]
      real(8), intent( in ) :: mb     !< nh3 concentration [ppt]
      real(8), intent( in ) :: z      !< height above the earth's surface [km]

      ! output arguments.

      real(8) :: jgas !< homogeneous nucleation rate [#/cm^3/s]
      real(8) :: jion !< ion-ion     nucleation rate [#/cm^3/s]
      real(8) :: jtot !< total       nucleation rate [#/cm^3/s]


      select case (inuc)
     
      case (1)     ! inuc=1: jvm binary nucleation scheme

        call Nucl_Jvm      (na,t,rh,jgas)

      case (2)     ! inuc=2: vehkamaki binary nucleation scheme

        call Nucl_Vehkamaki(na,t,rh,jgas)

      case (3)     ! inuc=3: napari ternary nucleation scheme

        if(mb.lt.0.1d+00) then
          call Nucl_Vehkamaki(na,t,rh,jgas)
        else
          call Nucl_Napari(na,mb,t,rh,jgas)
        endif
      
      case (4)     ! inuc=4: lines from plots in eisele and mcmurry, 1997.

        call Nucl_Eisele_Mcmurry(na,jgas)

      end select

      if ( include_ion_ion ) then
        call Nucl_Turco (na,z,jion)    ! ion-ion recombination scheme
      else
        jion = 0.0d+00
      endif

      if( inuc.eq.1 .or. inuc.eq.3 .or. inuc.eq.4 ) then
        jtot = min( max( jgas+jion, j_lower ), j_upper )
      elseif( inuc.eq.2 ) then   ! higher upper limit for vehkamaki et al. 2002
        jtot = min( max( jgas+jion, j_lower ), 1.0d+10 ) 
      endif

      return
      end subroutine Nucleation_Rate


!>-----------------------------------------------------------------------------------------------------------------------
!!     lsc/dlw 2005-2006:
!!-----------------------------------------------------------------------------------------------------------------------
      subroutine Nucl_Napari(na,mb,t,rh,j)
      implicit none

      ! input arguments.

      real(8) :: na   !< h2so4 concentration [molecules/cm^3]
      real(8) :: mb   !< nh3 concentration [ppt]
      real(8) :: t    !< temperature [k]
      real(8) :: rh   !< relative humidity [%]

      ! output arguments.

      real(8) :: j    !< nucleation rate [#/cm^3/s]

      ! local variables.

      logical :: valid_input

      ! parameters.

      real(8), parameter :: pi = 3.141592653589793d+00
      real(8), parameter :: avo = 6.0221367d+23
      real(8), parameter :: mwh2so4 = 98.07948d+00
      real(8), parameter :: mwnh3   = 17.03356d+00
      real(8), parameter :: mwh2o   = 18.01528d+00

      mb = min ( mb, 1.00d+02 )   ! cap at the maximum value for which the parameterization is valid.
      na = min ( na, 1.00d+09 )   ! cap at the maximum value for which the parameterization is valid.

      ! check the conditions of validity for input parameters.

      valid_input = .true.
      if     ( t  .lt. 240.00d+00  .or.  t  .gt. 300.00d+00 ) then
        valid_input = .false.
      elseif ( rh .lt.   0.50d+00  .or.  rh .gt.  95.00d+00 ) then
        valid_input = .false.
      elseif ( na .lt.   1.00d+04  .or.  na .gt.   1.00d+09 ) then    ! upper limit should have no effect.
        valid_input = .false.
      elseif ( mb .lt.   1.00d-01  .or.  mb .gt.   1.00d+02 ) then    ! upper limit should have no effect.
        valid_input = .false.
      endif

      if ( .not. valid_input ) then
        j = j_lower
        return
      endif

      j = j_napari(na,mb,t,rh)
      if ( j .lt. 1.0d-05 ) then      
        j = j_lower
        return
      elseif (  j .gt. 1.0d+06 ) then  
        j = j_upper
        return
      endif

      return
      end subroutine Nucl_Napari


!>-----------------------------------------------------------------------------------------------------------------------
!!     lsc/dlw 2005-2006:
!!-----------------------------------------------------------------------------------------------------------------------
      real(8) function Nh2so4_Napari(j,t)
      implicit none

      ! arguments.

      real(8) :: j    !< nucleation rate [#/cm^3/s]
      real(8) :: t    !< k
  
      ! scratch local variables.

      real(8) :: lnj

      lnj = log ( j )

      nh2so4_napari = 38.1645d+00 + 0.77410d+00*lnj +           &
                       0.00298879d+00*lnj*lnj-0.357605d+00*t -  &
                       0.00366358d+00*t*lnj + 0.0008553d+00*t*t
      nh2so4_napari = max (nh2so4_napari, 1.0d-30)

      return
      end function Nh2so4_Napari


!>-----------------------------------------------------------------------------------------------------------------------
!!     lsc/dlw 2005-2006:
!!-----------------------------------------------------------------------------------------------------------------------
      real(8) function Nnh3_Napari(j,t)
      implicit none

      ! arguments.

      real(8) :: j    !< nucleation rate [#/cm^3/s]
      real(8) :: t    !< k
  
      ! scratch local variables.

      real(8) :: lnj

      lnj = log ( j )

      nnh3_napari =   26.8982d+00 + 0.682905d+00*lnj +         &
                      0.0035752d+00*lnj*lnj -0.265748d+00*t -  &
                      0.00341895d+00*t*lnj + 0.000673454d+00*t*t
      nnh3_napari = max (nnh3_napari, 1.0d-30)

      return
      end function Nnh3_Napari


!>-----------------------------------------------------------------------------------------------------------------------
!!     lsc/dlw 2005-2006:
!!-----------------------------------------------------------------------------------------------------------------------
      real(8) function Ntot_Napari(j,t)
      implicit none

      ! arguments.

      real(8) :: t    !< [k]
      real(8) :: j    !< nucleation rate [#/cm^3/s]

      ! scratch local variables.

      real(8),parameter :: a = 79.3484d+00        ! coefficients
      real(8),parameter :: b = 1.7384d+00         ! coefficients
      real(8),parameter :: c = 0.00711403d+00     ! coefficients
      real(8),parameter :: d = -0.74493d+00       ! coefficients
      real(8),parameter :: e = -0.008202608d+00   ! coefficients
      real(8),parameter :: f = 0.0017855d+00      ! coefficients
      real(8) :: lnj

      lnj  = log ( j )

      ntot_napari = a + b*lnj + c*lnj*lnj + d*t + e*t*lnj + f*t*t 
      ntot_napari = max (ntot_napari, 1.0d-30)

      return
      end function Ntot_Napari


!>-----------------------------------------------------------------------------------------------------------------------
!!     lsc/dlw 2005-2006:
!!-----------------------------------------------------------------------------------------------------------------------
      real(8) function Rstar_Napari(j,t)
      implicit none

      ! arguments.

      real(8) :: t    !< temperature [k]
      real(8) :: j    !< nucleation rate [#/cm^3/s]

      ! scratch local variables.
      real(8) :: lnj    

      lnj = log(j)
      rstar_napari = 0.141027d+00 - 0.00122625d+00*lnj -       &
                     7.82211d-06*lnj*lnj - 0.00156727d+00*t -  &
                     0.00003076d+00*t*lnj + 0.0000108375d+00*t*t

      return
      end function Rstar_Napari


!>-----------------------------------------------------------------------------------------------------------------------
!!     lsc/dlw 2005-2006:
!!-----------------------------------------------------------------------------------------------------------------------
      real(8) function J_Napari(na,mb,t,rh)
      implicit none

      ! arguments.

      real(8) :: na   !< [molecules/cm^3]
      real(8) :: t    !< [k]
      real(8) :: rh   !< [%]
      real(8) :: mb   !< [ppt]

      ! scratch local variables.

      real(8) :: x,y,ex,ry,z,w,lnj
      real(8), dimension(20) :: f
      integer :: j

      ! parameters.

      real(8), dimension(20,0:3) :: a

      data a( 1,0:3)/  -0.355297000d+00,  -0.338449000d+02,   0.345360000d+00,  -0.824007000d-03/
      data a( 2,0:3)/   0.313735000d+01,  -0.772861000d+00,   0.561204000d-02,  -0.974576000d-05/
      data a( 3,0:3)/   0.190359000d+02,  -0.170957000d+00,   0.479808000d-03,  -0.414699000d-06/
      data a( 4,0:3)/   0.107605000d+01,   0.148932000d+01,  -0.796052000d-02,   0.761229000d-05/
      data a( 5,0:3)/   0.609160000d+01,  -0.125378000d+01,   0.939836000d-02,  -0.174927000d-04/
      data a( 6,0:3)/   0.311760000d+00,   0.164009000d+01,  -0.343852000d-02,  -0.109753000d-04/
      data a( 7,0:3)/  -0.200738000d-01,  -0.752115000d+00,   0.525813000d-02,  -0.898038000d-05/
      data a( 8,0:3)/   0.165536000d+00,   0.326623000d+01,  -0.489703000d-01,   0.146967000d-03/
      data a( 9,0:3)/   0.652645000d+01,  -0.258002000d+00,   0.143456000d-02,  -0.202036000d-05/
      data a(10,0:3)/   0.368024000d+01,  -0.204098000d+00,   0.106259000d-02,  -0.126560000d-05/
      data a(11,0:3)/  -0.665140000d-01,  -0.782382000d+01,   0.122938000d-01,   0.618554000d-04/
      data a(12,0:3)/   0.658740000d+00,   0.190542000d+00,  -0.165718000d-02,   0.341744000d-05/
      data a(13,0:3)/   0.599321000d-01,   0.596475000d+01,  -0.362432000d-01,   0.493370000d-04/
      data a(14,0:3)/  -0.732731000d+00,  -0.184179000d-01,   0.147186000d-03,  -0.237711000d-06/
      data a(15,0:3)/   0.728429000d+00,   0.364736000d+01,  -0.274220000d-01,   0.493478000d-04/
      data a(16,0:3)/   0.413016000d+02,  -0.357520000d+00,   0.904383000d-03,  -0.573788000d-06/
      data a(17,0:3)/  -0.160336000d+00,   0.889881000d-02,  -0.539514000d-04,   0.839522000d-07/
      data a(18,0:3)/   0.857868000d+01,  -0.112358000d+00,   0.472626000d-03,  -0.648365000d-06/
      data a(19,0:3)/   0.530167000d-01,  -0.198815000d+01,   0.157827000d-01,  -0.293564000d-04/
      data a(20,0:3)/  -0.232736000d+01,   0.234646000d-01,  -0.765190000d-04,   0.804590000d-07/

      ! statement function.

      real(8) :: z0,z1,z2,z3,zt
      z(z0,z1,z2,z3,zt) = z0 + z1*zt + z2*zt*zt + z3*zt*zt*zt
    
      x  = log ( rh * 0.01d+00 )
      ex = rh * 0.01d+00
      y  = log ( na )
      ry = 1.00d+00/y
      w  = log ( mb ) 

      do j=1,20
        f(j) = z(a(j,0),a(j,1),a(j,2),a(j,3),t)
      enddo
 
      lnj= -84.7551d+00 + f(1)*ry + f(2)*y + f(3)*y*y + f(4)*w         &
                        + f(5)*w*w + f(6)*ex + f(7)*x + f(8)*w*ry      &
                        + f(9)*w*y + f(10)*ex*y + f(11)*ex*ry          &
                        + f(12)*ex*w + f(13)*x*ry + f(14)*x*w          &
                        + f(15)*w*w*ry + f(16)*y*w*w + f(17)*y*y*w     &
                        + f(18)*ex*w*w + f(19)*ex*w*ry + f(20)*y*y*w*w

      j_napari = exp ( lnj )

      return
      end function J_Napari


!>-----------------------------------------------------------------------------------------------------------------------
!!     lsc/dlw 2005-2006:
!!     dlw:062005: created and checked for j value for 82 points in 
!!                 vehkamaki et al. 2002.
!!-----------------------------------------------------------------------------------------------------------------------
      subroutine Nucl_Vehkamaki(na,t,rh,j)
      implicit none

      ! input arguments.

      real(8) :: na   !< h2so4 concentration [molecules/cm^3]
      real(8) :: t    !< temperature [k]
      real(8) :: rh   !< relative humidity [%]

      ! output arguments.

      real(8) :: j    !< nucleation rate [#/cm^3/s]

      ! local variables.

      real(8) :: x    ! mole fraction h2so4
      real(8) :: n    ! number of molecules in the critical nucleus
      logical :: valid_input

      ! parameters.

      real(8), parameter :: pi = 3.141592653589793d+00
      real(8), parameter :: avo = 6.0221367d+23
      real(8), parameter :: mwh2so4 = 98.07948d+00
      real(8), parameter :: mwh2o   = 18.01528d+00


      ! check for conditions of validity for input parameters.
 
      na = min ( na, 1.00d+11 )  ! cap at the maximum value valid for the parameterization. 

      valid_input = .true.
      if     ( t  .lt. 190.00d+00  .or.  t  .gt. 305.15d+00 ) then
        valid_input = .false.
      elseif ( rh .lt.   0.01d+00  .or.  rh .gt. 100.000001d+00 ) then
        valid_input = .false.
      elseif ( na .lt.   1.00d+04  .or.  na .gt.   1.00d+11 ) then     ! upper limit should have no effect. 
        valid_input = .false.
      endif
      if ( .not. valid_input ) then
        j = j_lower
        return
      endif

      x = xstar_vehkamaki(na,t,rh)

      j = j_vehkamaki(na,t,rh,x)
      if ( j .lt. 1.0d-07 ) then 
        j = j_lower
        return
      elseif ( j .gt. 1.0d+10 ) then 
        j = 1.0d+10                     ! maximum value valid for this parameterization.
        return
      endif

      ! properties of the critical nucleus.

      n = ntot_vehkamaki(na,t,rh,x)
      if ( n .lt. 4.0d+00 ) then        ! check for condition on ntot.
        j = j_lower
        return
      endif

      return
      end subroutine Nucl_Vehkamaki


!>-----------------------------------------------------------------------------------------------------------------------
!!     lsc/dlw 2005-2006:
!!-----------------------------------------------------------------------------------------------------------------------
      real(8) function Xstar_Vehkamaki(na,t,rh)
      implicit none

      ! arguments.

      real(8) :: na   !< molecules/cm^3
      real(8) :: t    !< k
      real(8) :: rh   !< %

      ! scratch local variables.

      real(8) :: x

      x = log ( rh * 0.01d+00 )

      xstar_vehkamaki =   0.740997d+00       - 0.00266379d+00   *t        &
             + ( 0.0000504022d+00*t - 0.00349998d+00      ) * log ( na )  &
             + ( 0.00201048d+00     - 0.000183289d+00  *t ) * x           &
             + ( 0.00157407d+00     - 0.0000179059d+00 *t ) * x*x         &
             + ( 0.000184403d+00    - 1.50345d-06      *t ) * x**3

      return
      end function Xstar_Vehkamaki


!>-----------------------------------------------------------------------------------------------------------------------
!!     lsc/dlw 2005-2006:
!!-----------------------------------------------------------------------------------------------------------------------
      real(8) function Ntot_Vehkamaki(na,t,rh,xs)
      implicit none

      ! arguments.

      real(8) :: na   !< [molecules/cm^3]
      real(8) :: t    !< [k]
      real(8) :: rh   !< [%]
      real(8) :: xs   !< mole fraction h2so4 in the critical nucleus

      ! scratch local variables.

      real(8) :: a,b,c,d,e,f,g,h,i,j    ! coefficients
      real(8) :: x,y,rx

      ! parameters.

      real(8), parameter :: a1 = - 0.00295413d+00
      real(8), parameter :: a2 = - 0.0976834d+00
      real(8), parameter :: a3 = + 0.00102485d+00
      real(8), parameter :: a4 = - 2.18646d-06
      real(8), parameter :: a5 = - 0.101717d+00

      real(8), parameter :: b1 = - 0.00205064d+00
      real(8), parameter :: b2 = - 0.00758504d+00
      real(8), parameter :: b3 = + 0.000192654d+00
      real(8), parameter :: b4 = - 6.7043d-07 
      real(8), parameter :: b5 = - 0.255774d+00

      real(8), parameter :: c1 = + 0.00322308d+00
      real(8), parameter :: c2 = + 0.000852637d+00
      real(8), parameter :: c3 = - 0.0000154757d+00
      real(8), parameter :: c4 = + 5.66661d-08
      real(8), parameter :: c5 = + 0.0338444d+00

      real(8), parameter :: d1 = + 0.0474323d+00 
      real(8), parameter :: d2 = - 0.000625104d+00
      real(8), parameter :: d3 = + 2.65066d-06 
      real(8), parameter :: d4 = - 3.67471d-09
      real(8), parameter :: d5 = - 0.000267251d+00

      real(8), parameter :: e1 = - 0.0125211d+00
      real(8), parameter :: e2 = + 0.00580655d+00
      real(8), parameter :: e3 = - 0.000101674d+00
      real(8), parameter :: e4 = + 2.88195d-07
      real(8), parameter :: e5 = + 0.0942243d+00

      real(8), parameter :: f1 = - 0.038546d+00
      real(8), parameter :: f2 = - 0.000672316d+00
      real(8), parameter :: f3 = + 2.60288d-06    
      real(8), parameter :: f4 = + 1.19416d-08
      real(8), parameter :: f5 = - 0.00851515d+00

      real(8), parameter :: g1 = - 0.0183749d+00
      real(8), parameter :: g2 = + 0.000172072d+00
      real(8), parameter :: g3 = - 3.71766d-07    
      real(8), parameter :: g4 = - 5.14875d-10
      real(8), parameter :: g5 = + 0.00026866d+00

      real(8), parameter :: h1 = - 0.0619974d+00
      real(8), parameter :: h2 = + 0.000906958d+00
      real(8), parameter :: h3 = - 9.11728d-07    
      real(8), parameter :: h4 = - 5.36796d-09
      real(8), parameter :: h5 = - 0.00774234d+00

      real(8), parameter :: i1 = + 0.0121827d+00
      real(8), parameter :: i2 = - 0.00010655d+00
      real(8), parameter :: i3 = + 2.5346d-07      
      real(8), parameter :: i4 = - 3.63519d-10 
      real(8), parameter :: i5 = + 0.000610065d+00

      real(8), parameter :: j1 = + 0.000320184d+00
      real(8), parameter :: j2 = - 0.0000174762d+00 
      real(8), parameter :: j3 = + 6.06504d-08      
      real(8), parameter :: j4 = - 1.42177d-11
      real(8), parameter :: j5 = + 0.000135751d+00

      ! statement function.

      real(8) :: z
      real(8) :: z1,z2,z3,z4,z5,zt,zx
      z(z1,z2,z3,z4,z5,zt,zx) = z1 + z2*zt + z3*zt*zt + z4*zt**3 + z5*zx

      rx = 1.0d+00/xs
      x  = log ( rh * 0.01d+00 )
      y  = log ( na )

      a = z(a1,a2,a3,a4,a5,t,rx)
      b = z(b1,b2,b3,b4,b5,t,rx)
      c = z(c1,c2,c3,c4,c5,t,rx)
      d = z(d1,d2,d3,d4,d5,t,rx)
      e = z(e1,e2,e3,e4,e5,t,rx)
      f = z(f1,f2,f3,f4,f5,t,rx)
      g = z(g1,g2,g3,g4,g5,t,rx)
      h = z(h1,h2,h3,h4,h5,t,rx)
      i = z(i1,i2,i3,i4,i5,t,rx)
      j = z(j1,j2,j3,j4,j5,t,rx)

      ntot_vehkamaki = exp ( a + b*x + c*x*x + d*x**3 + e*y + f*x*y + &
                             g*x*x*y + h*y*y + i*x*y*y + j*y**3 )

      return
      end function Ntot_Vehkamaki


!>-----------------------------------------------------------------------------------------------------------------------
!!     lsc/dlw 2005-2006:
!!-----------------------------------------------------------------------------------------------------------------------
      real(8) function Rstar_Vehkamaki(xs,ntot)
      implicit none

      ! arguments.

      real(8) :: ntot  !< # molecules in the critical nucleus
      real(8) :: xs    !< mole fraction h2so4

      rstar_vehkamaki = exp( -1.6524245d+00 + 0.42316402d+00*xs &
                                           + 0.3346648d+00*log( ntot ) )

      return
      end function Rstar_Vehkamaki


!>-----------------------------------------------------------------------------------------------------------------------
!!     lsc/dlw 2005-2006:
!!-----------------------------------------------------------------------------------------------------------------------
      real(8) function J_Vehkamaki(na,t,rh,xs)
      implicit none

      ! arguments.

      real(8) :: na   !< [molecules/cm^3]
      real(8) :: t    !< [k]
      real(8) :: rh   !< [%]
      real(8) :: xs   !< mole fraction h2so4 in the critical nucleus

      ! scratch local variables.

      real(8) :: a,b,c,d,e,f,g,h,i,j    ! coefficients
      real(8) :: x,y,rx

      ! parameters.

      real(8), parameter :: a1 = + 0.14309d+00
      real(8), parameter :: a2 = + 2.21956d+00
      real(8), parameter :: a3 = - 0.0273911d+00
      real(8), parameter :: a4 = + 0.0000722811d+00
      real(8), parameter :: a5 = + 5.91822d+00

      real(8), parameter :: b1 = + 0.117489d+00
      real(8), parameter :: b2 = + 0.462532d+00
      real(8), parameter :: b3 = - 0.0118059d+00
      real(8), parameter :: b4 = + 0.0000404196d+00
      real(8), parameter :: b5 = + 15.7963d+00

      real(8), parameter :: c1 = - 0.215554d+00
      real(8), parameter :: c2 = - 0.0810269d+00
      real(8), parameter :: c3 = + 0.00143581d+00
      real(8), parameter :: c4 = - 4.7758d-06
      real(8), parameter :: c5 = - 2.91297d+00

      real(8), parameter :: d1 = - 3.58856d+00 
      real(8), parameter :: d2 = + 0.049508d+00
      real(8), parameter :: d3 = - 0.00021382d+00
      real(8), parameter :: d4 = + 3.10801d-07
      real(8), parameter :: d5 = - 0.0293333d+00

      real(8), parameter :: e1 = + 1.14598d+00
      real(8), parameter :: e2 = - 0.600796d+00
      real(8), parameter :: e3 = + 0.00864245d+00  
      real(8), parameter :: e4 = - 0.0000228947d+00
      real(8), parameter :: e5 = - 8.44985d+00  

      real(8), parameter :: f1 = + 2.15855d+00
      real(8), parameter :: f2 = + 0.0808121d+00
      real(8), parameter :: f3 = - 0.000407382d+00
      real(8), parameter :: f4 = - 4.01957d-07
      real(8), parameter :: f5 = + 0.721326d+00

      real(8), parameter :: g1 = + 1.6241d+00
      real(8), parameter :: g2 = - 0.0160106d+00
      real(8), parameter :: g3 = + 0.0000377124d+00
      real(8), parameter :: g4 = + 3.21794d-08
      real(8), parameter :: g5 = - 0.0113255d+00

      real(8), parameter :: h1 = + 9.71682d+00
      real(8), parameter :: h2 = - 0.115048d+00
      real(8), parameter :: h3 = + 0.000157098d+00
      real(8), parameter :: h4 = + 4.00914d-07 
      real(8), parameter :: h5 = + 0.71186d+00

      real(8), parameter :: i1 = - 1.05611d+00
      real(8), parameter :: i2 = + 0.00903378d+00
      real(8), parameter :: i3 = - 0.0000198417d+00
      real(8), parameter :: i4 = + 2.46048d-08 
      real(8), parameter :: i5 = - 0.0579087d+00

      real(8), parameter :: j1 = - 0.148712d+00
      real(8), parameter :: j2 = + 0.00283508d+00 
      real(8), parameter :: j3 = - 9.24619d-06     
      real(8), parameter :: j4 = + 5.00427d-09
      real(8), parameter :: j5 = - 0.0127081d+00   

      ! statement function.

      real(8) :: z
      real(8) :: z1,z2,z3,z4,z5,zt,zx
      z(z1,z2,z3,z4,z5,zt,zx) = z1 + z2*zt + z3*zt*zt + z4*zt**3 + z5*zx

      rx = 1.0d+00/xs
      x  = log ( rh * 0.01d+00 )
      y  = log ( na )

      a = z(a1,a2,a3,a4,a5,t,rx)
      b = z(b1,b2,b3,b4,b5,t,rx)
      c = z(c1,c2,c3,c4,c5,t,rx)
      d = z(d1,d2,d3,d4,d5,t,rx)
      e = z(e1,e2,e3,e4,e5,t,rx)
      f = z(f1,f2,f3,f4,f5,t,rx)
      g = z(g1,g2,g3,g4,g5,t,rx)
      h = z(h1,h2,h3,h4,h5,t,rx)
      i = z(i1,i2,i3,i4,i5,t,rx)
      j = z(j1,j2,j3,j4,j5,t,rx)

      j_vehkamaki = exp ( a + b*x + c*x*x + d*x**3 + e*y + f*x*y + &
                          g*x*x*y + h*y*y + i*x*y*y + j*y**3 )

      return
      end function J_Vehkamaki


!>-----------------------------------------------------------------------------------------------------------------------
!!     lsc/dlw 2005-2006:
!!-----------------------------------------------------------------------------------------------------------------------
      subroutine Nucl_Jvm(na,t,rh,j)
      implicit none

      ! input arguments.

      real(8) :: na   !< h2so4 concentration [molecules/cm^3]
      real(8) :: t    !< temperature [k]
      real(8) :: rh   !< relative humidity [%]

      ! output arguments.

      real(8) :: j    ! nucleation rate [#/cm^3/s]

      ! local variables.

      real(8), dimension(0:2) :: b    ! coefficient of fitted curve function
      logical :: valid_input
      real(8) :: log10j, xx, xlog10na

      na = min ( na, 1.00d+14 ) 

      ! check for conditions of validity for input parameters.

      call Jvm_Coef (t,rh,b)

      xx = -b(1) / (2.0d+00 * b(2))
      xlog10na = log10(na)

      valid_input = .true.
      if(xlog10na .gt. xx) then
        valid_input = .false.
      elseif ( na .lt. 1.00d+04  .or.  na .gt. 1.00d+14 ) then     ! upper limit should have no effect. 
        valid_input = .false.
      endif

      if ( .not. valid_input ) then
        j = j_lower
        return
      endif

      ! calculate the nucleation rate.

      log10j = b(0) + b(1)*xlog10na + b(2)*xlog10na*xlog10na

      j = 10.0**( log10j )

      if ( j .lt. 1.0d-03 ) then 
        j = j_lower 
        return
      elseif (  j .gt. 1.0d+05 ) then 
        j = j_upper  
        return
      endif

      return
      end subroutine Nucl_Jvm


!-----------------------------------------------------------------------------------------------------------------------
!     lsc/dlw 2005-2006:
!-----------------------------------------------------------------------------------------------------------------------
      subroutine Jvm_Coef(t,rh,b)
      implicit none

      ! input arguments.

      real(8) :: t    !< [k]
      real(8) :: rh   !< [%]

      ! output arguments.

      real(8), dimension(0:2) :: b !< ?

      ! scratch local variables.

      real(8) :: x,y,z
      real(8), dimension(0:4) :: bb0,bb1,bb2

      ! parameters.

      real(8), dimension(0:4,0:4) :: c0,c1,c2 
      integer :: j

      data c0(0,0:4)/  -0.493166930d+03,  -0.746060630d+03,   0.427001110d+04,   0.403205190d+04,  -0.103305460d+05/
      data c0(1,0:4)/   0.227147650d+04,   0.296173760d+04,  -0.345957800d+05,  -0.435346380d+05,   0.583845900d+05/
      data c0(2,0:4)/  -0.596070430d+04,  -0.921088590d+04,   0.919397290d+05,   0.131516870d+06,  -0.130410870d+06/
      data c0(3,0:4)/   0.673050820d+04,   0.112160080d+05,  -0.101584900d+06,  -0.151340760d+06,   0.133574610d+06/
      data c0(4,0:4)/  -0.261458200d+04,  -0.444305610d+04,   0.394068830d+05,   0.594771880d+05,  -0.496408160d+05/
      data c1(0,0:4)/   0.100665720d+03,   0.108313420d+03,  -0.113683880d+04,  -0.888430520d+03,   0.302914500d+04/
      data c1(1,0:4)/  -0.504046310d+03,  -0.474794750d+03,   0.904264540d+04,   0.909745930d+04,  -0.201185820d+05/
      data c1(2,0:4)/   0.137492990d+04,   0.161927870d+04,  -0.243558850d+05,  -0.269705410d+05,   0.498079760d+05/
      data c1(3,0:4)/  -0.156449310d+04,  -0.200427440d+04,   0.270531380d+05,   0.307133920d+05,  -0.531638600d+05/
      data c1(4,0:4)/   0.604879820d+03,   0.788620250d+03,  -0.104754690d+05,  -0.119805370d+05,   0.200426110d+05/
      data c2(0,0:4)/  -0.516074080d+01,  -0.398706400d+01,   0.720177880d+02,   0.509391250d+02,  -0.211700440d+03/
      data c2(1,0:4)/   0.284462970d+02,   0.212336580d+02,  -0.573120500d+03,  -0.497927540d+03,   0.152757040d+04/
      data c2(2,0:4)/  -0.795352710d+02,  -0.764306230d+02,   0.155886490d+04,   0.143625300d+04,  -0.392583270d+04/
      data c2(3,0:4)/   0.904137910d+02,   0.939499970d+02,  -0.173349330d+04,  -0.160777790d+04,   0.422689150d+04/
      data c2(4,0:4)/  -0.345602640d+02,  -0.362758330d+02,   0.667993290d+03,   0.619425480d+03,  -0.159032770d+04/

      ! statement function.

      real(8) :: z0,z1,z2,z3,z4,v
      z(z0,z1,z2,z3,z4,v) = z0+z1*v+z2*v*v+z3*v*v*v+z4*v*v*v*v

      x  = ( t - 273.15d+00 ) * 0.01d+00
      y  = rh * 0.01d+00

      do j = 0, 4
        bb0(j) = z(c0(j,0),c0(j,1),c0(j,2),c0(j,3),c0(j,4),x)
      enddo
      do j = 0, 4
        bb1(j) = z(c1(j,0),c1(j,1),c1(j,2),c1(j,3),c1(j,4),x)
      enddo
      do j = 0, 4
        bb2(j) = z(c2(j,0),c2(j,1),c2(j,2),c2(j,3),c2(j,4),x)
      enddo

      b(0) = z(bb0(0),bb0(1),bb0(2),bb0(3),bb0(4),y)
      b(1) = z(bb1(0),bb1(1),bb1(2),bb1(3),bb1(4),y)
      b(2) = z(bb2(0),bb2(1),bb2(2),bb2(3),bb2(4),y)

      return
      end subroutine Jvm_Coef


!>-----------------------------------------------------------------------------------------------------------------------
!!     lsc/dlw 2005-2006:
!!-----------------------------------------------------------------------------------------------------------------------
      subroutine Nucl_Turco(na,z,j)
      implicit none

      ! input arguments.

      real(8) :: na   !< h2so4 concentration [molecules/cm^3]
      real(8) :: z    !< geographical height [km]

      ! output arguments.

      real(8) :: j    !< nucleation rate [#/cm^3/s]

      j = j_turco(na,z)

      return
      end subroutine Nucl_Turco


!>-----------------------------------------------------------------------------------------------------------------------
!!     lsc/dlw 2005-2006:
!!-----------------------------------------------------------------------------------------------------------------------
      real(8) function J_Turco(na,z)
      implicit none

      ! arguments.

      real(8) :: na   !< h2so4 concentration [molecules/cm^3]
      real(8) :: z    !< geographical height [km]

      ! local variable
      real(8) :: q    ! local ionization rate [/cm^3/s]

      ! parameters.

      real(8), parameter :: na0 = 5.00d+06
      real(8), parameter :: f0 = 1.00d-03
      integer, parameter :: nstar = 3

 
      if(z .lt. 2.25) then
        q = 2.00d+00
      elseif(z .ge. 2.25d+00 .and. z .lt. 11.00d+00) then
        q = 3.10d+00 * (z - 2.25d+00) + 2.00d+00
      elseif ( z .ge. 11.00d+00) then
        q = 29.00d+00 + (z - 11.00d+00)
      endif

      j_turco = q * f0 * ( na / na0 )**nstar
      j_turco = min( q, j_turco )

      return
      end function J_Turco


!>-----------------------------------------------------------------------------------------------------------------------
!!     dlw: 7-27-06: these expressions were derived from figure 7 of 
!!          eisele, f. l., and mcmurry, p. h. (1997). 
!!          recent progress in understanding particle nucleation and growth,
!!          phil. trans. r. soc. lond. b, 352, 191-201.
!!-----------------------------------------------------------------------------------------------------------------------
      subroutine Nucl_Eisele_Mcmurry(na,j)
      implicit none

      ! input arguments.

      real(8) :: na   !< h2so4 concentration [molecules/cm^3]

      ! output arguments.

      real(8) :: j    !< nucleation rate [#/cm^3/s]

      ! parameters.

      real(8), parameter :: k1 = 5.8d-13
      real(8), parameter :: k2 = 3.5d-15
      real(8), parameter :: k3 = 3.7d-14
      logical, parameter :: lower_curve = .false.
      logical, parameter :: upper_curve = .false.
      logical, parameter :: avgul_curve = .true.

      ! local variables.

      ! integer :: i
      ! real(8) :: s

      if(lower_curve) j = k1*na
      if(upper_curve) j = k2*na*na
      if(avgul_curve) j = k3*na**1.5
       
!-----------------------------------------------------------------------------------------------------------------------
!     plot the expressions.
!-----------------------------------------------------------------------------------------------------------------------
!     do i=1,5
!       s = 10.0**(4+i-1)
!       write(78,91 )s, k1*s, k3*s**1.5, k2*s*s
!     enddo
!-----------------------------------------------------------------------------------------------------------------------

91    format(4d15.5)
      return
      end subroutine Nucl_Eisele_Mcmurry


!>-----------------------------------------------------------------------------------------------------------------------
!!     dlw, 8-1-06.
!!     routine to calculate the ratio of the new particle formation rate at
!!     a user-specified diameter to the nucleation rate.
!!-----------------------------------------------------------------------------------------------------------------------
      subroutine F_Kk02(prs,tk,rhpercent,napercm3,nh3ppt,kc,dnpf_nm,dnuc_nm,jp_to_j)
      implicit none

      ! input arguments.

      real(8) :: prs       !< ambient pressure               [pa]
      real(8) :: tk        !< ambient temperature            [k]
      real(8) :: rhpercent !< relative humidity              [%]
      real(8) :: napercm3  !< sulfuric acid concentration    [#/cm^3]
      real(8) :: nh3ppt    !< ammonia concentration          [ppt]         
      real(8) :: kc        !< condensational sink            [1/s]         
      real(8) :: dnpf_nm   !< user-selected diameter for new particles    [nm]         
      real(8) :: dnuc_nm   !< initial diameter before growth to size dnpf [nm]         

      ! output arguments.

      real(8) :: jp_to_j   !< ratio of new particle formation rate at diameter dnp
                           !< to the nucleation rate. [1]

      ! local scratch variables.

      real(8) :: m_effective  ! effective molar mass of h2so4 including the water
                              ! and ammonia that instantaneously equilibrate with
                              ! the particle [g/mol]
      real(8) :: gr_inorg     ! growth rate in the free-molecular regime due to 
                              ! inorganic species [nm/h]
      real(8) :: x_nh3        ! # of nh3 molecules condensing per h2so4 molecule
                              ! - lies in the range 0-2 [1]
      real(8) :: x_h2o        ! # of h2o molecules condensing per h2so4 molecule [1]
      real(8) :: gamma        ! eq.(22) of kk2002. [nm^2 m^2 h^-1]
      real(8) :: dmean_nm     ! number mean diameter over all modes [nm]
      real(8) :: cs_prime     ! condensational sink espressed as  s [m^-2]
      real(8) :: eta          ! eta parameter of the kk2002 model   [nm]

      !-----------------------------------------------------------------------------------------------------------------
      ! the growth rate due to h2so4+h2o+nh3 condensation can be expressed as
      !
      !     gr(nm/h) = gr_const * t^0.5 * meff(g/mol) * c(#/cm^3) 
      !
      ! where meff is the effective molar mass of h2so4 with water and ammonia
      ! instantaneously equilibrating to the particle. see notes of 7-19-05.
      ! free-molecular growth is assumed.
      ! 
      ! alpha is the mass accommodation coefficient[1]. see eq.(20) of kk02:
      !   alpha should be unity here, even if less than that elsewhere.  
      ! 
      ! dstar_density is the density of a critical nucleus [g/cm^3].
      !
      ! the mean thermal speed of a h2so4 molecule [m/s] is cmean_h2so4 * t^0.5,
      !   with t the absolute temperature [k].
      !
      ! avo is avogadro's number [#/mol]. 
      !
      ! the factor 3.6d+12 converts [m/s] to [nm/h].
      !-----------------------------------------------------------------------------------------------------------------
      real(8), parameter :: alpha = 1.0d+00                           ! [1]
      real(8), parameter :: dstar_density = 1.6d+00                   ! [g/cm^3]
      real(8), parameter :: cmean_h2so4 = 14.692171525317413d+00      ! [m/s/k^0.5]
      real(8), parameter :: gr_const = 0.5d+00 * alpha * cmean_h2so4 * 3.6d+12 &
                                     / ( avo * dstar_density )
      real(8), parameter :: fourpi = 4.0d+00 * pi 

      !-----------------------------------------------------------------------------------------------------------------
      ! the proportionality factor gamma for eq.22 of kk2002 have this 
      ! constant extracted.
      !-----------------------------------------------------------------------------------------------------------------
      real(8), parameter :: gamma_prime = 11.79348270d+00


      !-----------------------------------------------------------------------------------------------------------------
      ! get the number mean diameter over all modes [nm].
      !-----------------------------------------------------------------------------------------------------------------
      dmean_nm = 1.0d+09 * avg_dp_of_avg_mass_meters   ! stored in aero_param.f.

      !-----------------------------------------------------------------------------------------------------------------
      ! get gamma from eq.(22) in kk2002 [nm^2 m^2 /h].
      !
      ! dstar_density is in [g/cm^3], with the conversion from [kg/m^3] in
      ! kk02 eq.(22) folded into the parameter gamma_prime.
      !-----------------------------------------------------------------------------------------------------------------
      gamma = gamma_prime * dnuc_nm**0.2 * dnpf_nm**0.075 * dmean_nm**0.048  &
                          * dstar_density**(-0.33)        * tk**(-0.75)

      !-----------------------------------------------------------------------------------------------------------------
      ! get the effective molecular weight of condensing h2so4 [g/mol].
      !-----------------------------------------------------------------------------------------------------------------
      call Effective_Mw(prs,tk,rhpercent,napercm3,nh3ppt,x_nh3,x_h2o)
      m_effective = mw_h2so4 + x_nh3*mw_nh3 + x_h2o*mw_h2o

      !-----------------------------------------------------------------------------------------------------------------
      ! get the growth rate due to condensation of inorganics [nm/h].
      !-----------------------------------------------------------------------------------------------------------------
      gr_inorg = gr_const * sqrt(tk) * m_effective * napercm3 

      !-----------------------------------------------------------------------------------------------------------------
      ! get the condensational sink (cs') in [m^-2].
      !-----------------------------------------------------------------------------------------------------------------
      cs_prime = kc / ( fourpi * diffcoef_m2s(ilay) ) 

      !-----------------------------------------------------------------------------------------------------------------
      ! get eta of the kk2002 formulation [nm].
      !-----------------------------------------------------------------------------------------------------------------
      eta = gamma * cs_prime / max( gr_inorg, 1.0d-30 )
      ! write(34,*)'eta,cs_prime,gr_inorg,napercm3 = ',eta,cs_prime,gr_inorg,napercm3 

      if(eta .lt. 56.0d+00) then
        jp_to_j = exp( eta * ( (1.0d+00/dnpf_nm) - (1.0d+00/dnuc_nm) ) )
      else 
        jp_to_j = 1.0d-17
      endif
      ! write(34,*)'jp_to_j = ', jp_to_j 

      !-----------------------------------------------------------------------------------------------------------------
      ! dlw, 8-1-06: all of these variables were fine.
      !-----------------------------------------------------------------------------------------------------------------
      if( write_f_kk02 ) then
        write(aunit1,'(/A/)')   'IN SUBROUTINE F_KK02'
        write(aunit1,*)         'DSTAR_DENSITY,TK,DMEAN_NM,GAMMA,DIFFCOEF_M2S(ILAY)'
        write(aunit1,'(5d15.5)') dstar_density,tk,dmean_nm,gamma,diffcoef_m2s(ilay)
        write(aunit1,*)         'X_H2O,X_NH3,M_EFFECTIVE,DNUC_NM,DNPF_NM'
        write(aunit1,'(5d15.5)') x_h2o,x_nh3,m_effective,dnuc_nm,dnpf_nm 
        write(aunit1,*)         'CS_PRIME,GR_INORG,ETA,JP_TO_J'
        write(aunit1,'(5d15.5)') cs_prime,gr_inorg,eta,jp_to_j 
        write(aunit1,'(A)')     '  '
      ENDIF
      !-----------------------------------------------------------------------------------------------------------------

      return
      end subroutine F_Kk02


!>-----------------------------------------------------------------------------------------------------------------------
!!     routine to calculate x_nh3 and x_h2o as defined below.
!!-----------------------------------------------------------------------------------------------------------------------
      subroutine Effective_Mw(airprs,tk,rhpercent,napercm3,nh3ppt,x_nh3,x_h2o)
      implicit none

      ! input arguments.

      real(8) :: airprs    !< ambient pressure               [pa]
      real(8) :: tk        !< ambient temperature            [k]
      real(8) :: rhpercent !< relative humidity              [%]
      real(8) :: napercm3  !< sulfuric acid concentration    [#/cm^3]
      real(8) :: nh3ppt    !< ammonia mixing ratio           [ppt]

      ! output arguments.

      real(8) :: x_nh3     !< # of nh3 molecules condensing per h2so4 [#]
      real(8) :: x_h2o     !< # of h2o molecules condensing per h2so4 [#]

      ! parameters.

!     real(8), parameter :: dnu_nm   = 10.0d+00        ! diameter of new particles [nm]
      real(8), parameter :: rrnu_nm  = 2.0d+00/dnu_nm  ! recipocal radius ... [nm^-1]
      real(8), parameter :: sqrt_mwh2so4_over_mwnh3 = 2.400d+00
      real(8), parameter :: a0_h2so4 =  0.300d+00
      real(8), parameter :: a1_h2so4 = -0.608d+00
      real(8), parameter :: a2_h2so4 =  0.701d+00
      real(8), parameter :: a3_h2so4 = -0.392d+00
      real(8), parameter :: a0_biso4 =  0.866d+00
      real(8), parameter :: a1_biso4 = -1.965d+00
      real(8), parameter :: a2_biso4 =  1.897d+00
      real(8), parameter :: a3_biso4 = -0.800d+00
      real(8), parameter :: a0_amso4 =  0.951d+00
      real(8), parameter :: a1_amso4 = -2.459d+00
      real(8), parameter :: a2_amso4 =  2.650d+00
      real(8), parameter :: a3_amso4 = -1.142d+00
!     logical, parameter :: include_kelvin_effect = .true.
      logical, parameter :: include_kelvin_effect = .false.

      ! local scratch variables.

      real(8)       :: nh3percm3      ! ammonia number concentration [# cm^-3]
      real(8)       :: mfs            ! mole fraction sulfur in the particle:
                                      !   ( mol s ) / ( mol s + mol h2o )
      real(8)       :: rhf            ! (fractional rh) * exp factor for kelvin effect
      real(8)       :: a              ! scratch variable for kelvin effect

      napercm3 = max ( napercm3, 1.0d-30 )
      nh3ppt   = max ( nh3ppt,   1.0d-30 )

      !-----------------------------------------------------------------------------------------------------------------
      ! convert nh3 from ppt to molecules cm^-3.
      !-----------------------------------------------------------------------------------------------------------------
      nh3percm3 = 7.25d+04 * nh3ppt * airprs / tk

      !-----------------------------------------------------------------------------------------------------------------
      ! see notes of 7-08-05 for the x_nh3 calculation.
      ! the square root arises from a ratio of the thermal velocities for 
      !   the two molecules. 
      !-----------------------------------------------------------------------------------------------------------------
      x_nh3 = min ( sqrt_mwh2so4_over_mwnh3 * nh3percm3 / napercm3, 2.0d+00 )

      !-----------------------------------------------------------------------------------------------------------------
      ! now use the x_nh3 value to classify the new particles as either
      !
      !  (1) acidic            - treat as pure sulfuric acid
      !  (2) half-neutralized  - treat as pure ammonium bisulfate
      !  (3) fully-neutralized - treat as pure ammonium sulfate
      !
      ! then compute the water uptake for the selected composition.
      !-----------------------------------------------------------------------------------------------------------------
      if ( x_nh3 .lt. 0.5d+00 ) then

        ! acidic case.

        npfmass_regime = 0                          ! second index to npfmass(:,:)
        if ( include_kelvin_effect ) then
          a = 1.2d+00 - 0.0072d+00*(tk-273.15d+00)
          rhf = 0.01d+00 * rhpercent * exp ( -a * rrnu_nm )
        else
          rhf = 0.01d+00 * rhpercent
        endif
        mfs = a0_h2so4 + a1_h2so4*rhf + a2_h2so4*rhf*rhf + a3_h2so4*rhf**3

      elseif ( x_nh3 .ge. 0.5d+00 .and. x_nh3 .le. 1.5d+00 ) then

        ! half-neutralized case.

        npfmass_regime = 1                          ! second index to npfmass(:,:)
        if ( include_kelvin_effect ) then
          rhf = 0.01d+00 * rhpercent
          a = ( 1.2d+00 - 0.0072d+00*(tk-273.15d+00) ) &
            * ( 1.0d+00 - 0.17d+00*(1.0d+00-rhf))      &
            * ( 1.0d+00 + 0.95d+00*(1.0d+00-rhf))
          rhf = rhf * exp ( -a * rrnu_nm )
        else
          rhf = 0.01d+00 * rhpercent
        endif
        mfs = a0_biso4 + a1_biso4*rhf + a2_biso4*rhf*rhf + a3_biso4*rhf**3

      else

        ! fully-neutralized case.

        npfmass_regime = 2                          ! second index to npfmass(:,:)
        if ( include_kelvin_effect ) then
          rhf = 0.01d+00 * rhpercent
          a = ( 1.2d+00 - 0.0072d+00*(tk-273.15d+00) ) &
            * ( 1.0d+00 - 0.17d+00*(1.0d+00-rhf))      &
            * ( 1.0d+00 + 0.95d+00*(1.0d+00-rhf))
          rhf = rhf * exp ( -a * rrnu_nm )
        else
          rhf = 0.01d+00 * rhpercent
        endif
        mfs = a0_amso4 + a1_amso4*rhf + a2_amso4*rhf*rhf + a3_amso4*rhf**3

      endif

      if ( mfs .gt. 0.0d+00 ) then   ! the mole fraction s in the particle should exceed zero.
        x_h2o = (1.0d+00 - mfs ) / mfs
      else                           ! this case should not occur.
        x_h2o = 10.0d+00
      endif

      if( write_log ) write(aunit1,'(/A,F6.2,I6/)') 'X_NH3, NPFMASS_REGIME = ', x_nh3, npfmass_regime
      return
      end subroutine Effective_Mw


!>------------------------------------------------------------------------------------------------------------------
!!     101706, dlw: routine to estimate the steady-state concentration of sulfuric acid
!!                  including the consumption of h2so4 by new particle formation during the current time step.
!!------------------------------------------------------------------------------------------------------------------
      subroutine Steady_State_H2so4(prs,rh,temp,xh2so4_ss,so4rate,xnh3,kc,dt,xh2so4_ss_wnpf)
      implicit none

      ! input arguments.

      real(8), intent(in)   :: prs              !< pressure [pa]
      real(8), intent(in)   :: rh               !< fractional relative humidity [1]
      real(8), intent(in)   :: temp             !< ambient temperature [k]
      real(8), intent(in)   :: xh2so4_ss        !< initial steady-state [h2so4] (as so4) 
                                                !<   neglecting new particle formation [ugso4/m^3]
      real(8), intent(in)   :: so4rate          !< gas-phase h2so4 (as so4) production rate [ugso4/m^3 s]
      real(8), intent(in)   :: xnh3             !< ammonia mixing ratio [ppmv]
      real(8), intent(in)   :: kc               !< condensational sink due to pre-existing aerosol [1/s]
      real(8), intent(in)   :: dt               !< model physics time step [s]

      ! output arguments.

      real(8), intent(out)  :: xh2so4_ss_wnpf   !< steady-state [h2so4] including new particle formation [ugso4/m^3]

      ! scratch local variables.

      integer :: i                              ! loop counter
      real(8) :: dndt                           ! new particle formation rate [particles/m^3/s]
      real(8) :: dmdt_so4                       ! npf mass production rates [ugso4/m^3/s]
      real(8) :: fx                             ! steady-state equation is fx = 0. [ugso4/m^3/s]
      integer, parameter :: itmax = 100         ! loop limit for development code
      integer, parameter :: icallnpfrate = 1    ! =0 impose mass limitation, >0 do not impose mass limitation
      real(8), parameter :: xh2so4_thres_ncm3 = 1.001d+04  ! if [h2so4] is below this npf can be neglected.[#/cm^3] 
      real(8), parameter :: xh2so4_thres = xh2so4_thres_ncm3 * mw_so4 * 1.0d+12 / avo   ! converted to [ugso4/m^3] 
      real(8), parameter :: eps_xh2so4_ncm3   = 1.00d+00                                ! tiny [h2so4] [#/cm^3] 
      real(8), parameter :: eps_xh2so4 = eps_xh2so4_ncm3 * mw_so4 * 1.0d+12 / avo       ! convert to [ugso4/m^3] 
      real(8), parameter :: reduction_factor = 1.2d+00   ! factor by which [h2so4]ss is reduced each iteration [1] 
      logical, parameter :: early_return = .false.       ! flag for no-operation early exit 

      if( xh2so4_ss.lt.xh2so4_thres .or. inuc.ne.3 .or. early_return ) then     ! inuc=3 is the napari et al. 
        xh2so4_ss_wnpf = xh2so4_ss                                              ! [ugso4/m^3]
        return
      endif
      xh2so4_ss_wnpf = xh2so4_ss + eps_xh2so4                                   ! [ugso4/m^3]
      call NpfRate(prs,rh,temp,xh2so4_ss_wnpf,so4rate,xnh3,kc,dndt,dmdt_so4,icallnpfrate) 
      fx = so4rate - kc*xh2so4_ss_wnpf - dmdt_so4                               ! evaluate function [ugso4/m^3/s]
      if( fx .gt. 0.0d+00 ) then
        ! write(34,*)'fx(xmax) .gt. 0.0d+00 in steady_state_h2so4: fx = ', fx
        return
      endif
      ! write(34,'(a,i5,5d13.4)')'i,x_ss,x_ss_wnpf,p,dmdt_so4,fx=',0,xh2so4_ss,xh2so4_ss_wnpf,so4rate,dmdt_so4,fx
!------------------------------------------------------------------------------------------------------------------
!     reduce the steady-state h2so4 until fx changes sign from negative to positive. 
!     then the current value of xh2so4_ss_wnpf is within a factor of reduction_factor 
!     of the actual steady-state value. 
!------------------------------------------------------------------------------------------------------------------
      do i=1, itmax
        xh2so4_ss_wnpf = xh2so4_ss_wnpf / reduction_factor                      ! [ugso4/m^3]
        call NpfRate(prs,rh,temp,xh2so4_ss_wnpf,so4rate,xnh3,kc,dndt,dmdt_so4,icallnpfrate) 
        fx = so4rate - kc*xh2so4_ss_wnpf - dmdt_so4                             ! evaluate function [ugso4/m^3/s]
        ! write(34,'(a,i5,5d13.4)')'i,x_ss,x_ss_wnpf,p,dmdt_so4,fx=',i,xh2so4_ss,xh2so4_ss_wnpf,so4rate,dmdt_so4,fx
        if( fx .ge. 0.0d+00 ) exit
        if( xh2so4_ss_wnpf .lt. eps_xh2so4 ) return                             ! [ugso4/m^3]
      enddo
      return
      end subroutine Steady_State_H2so4


      end module aero_npf


