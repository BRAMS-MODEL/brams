   module Aero_Depv
      use memMatrix,  only: &
                              nlays, &   !intent()
                              ixxx, &    !intent()
                              iyyy, &    !intent()
                              ilay , &      !intent()
                              nmodes     !intent()
!      use Amp_Aerosol, only: &
!                              vddep_aero !intent()
!-------------------------------------------------------------------------------------------------------------------------
!     The array VDDEP_AERO(X,Y,Z,I,1) contains current values for the dry deposition velocities 
!     for aerosol number concentrations for mode I. 
!     The array VDDEP_AERO(X,Y,Z,I,2) contains current values for the dry deposition velocities 
!     for aerosol mass   concentrations for mode I. 
!     Values in VDDEP_AERO are saved in subr. MATRIX at each time step. 
!-------------------------------------------------------------------------------------------------------------------------
      
      CONTAINS

!>----------------------------------------------------------------------------------------------------------------------
!!     calculate deposition velocity for aitken, accumulation, and coarse modes.
!!     reference: binkowski f. s., and u. shankar, the regional particulate model 
!!     1. model description and preliminary results. j. geophys. res., 100, d12, 26191-26209, 1995.
!!
!!     12-18-06, dlw: derived from the cmaq routine aero_depv.f originally coded by f. s. binkowski. :: 
!!----------------------------------------------------------------------------------------------------------------------
      subroutine Get_Aero_Depv(n,tk,rhoa,xlm,amu,wstar,ustar,ra,dgn_ddep,xls_ddep,den_ddep,vddep_aero)
      implicit none

      ! arguments. 

      integer :: n                    !< number of modes [1]
      real(8) :: tk                   !< air temperature [k]
      real(8) :: rhoa                 !< air density  [kg/m^3]
      real(8) :: xlm                  !< atmospheric mean free path [m]
      real(8) :: amu                  !< atmospheric dynamic viscosity [kg/(m s)]
      real(8) :: wstar                !< convective velocity scale [m/s]
      real(8) :: ustar                !< friction velocity [m/s]
      real(8) :: ra                   !< aerodynamic resistance [s/m]
      real(8) :: dgn_ddep(n)          !< geo. mean diameter    for each mode [um] 
      real(8) :: xls_ddep(n)          !< ln(geo. std. dev.)    for each mode [1]
      real(8) :: den_ddep(n)          !< avg. particle density for each mode [g/cm^3]
      REAL(8),INTENT(OUt)  :: vddep_aero(n,2) !< Dry deposition velocities [m/s]

      ! local variables. 

      integer :: i      
      real(8) :: dgn_m                ! geo. mean diameter    [m] 
      real(8) :: den_kgm3             ! avg. particle density [kg/m^3]
      real(8) :: vdep(2)              ! deposition velocities [m/s]
	!print*,"get dpv=",n;call flush(6)


      do i=1, n
         !print*,"get dpv1=",i,n;call flush(6)
        dgn_m    = dgn_ddep(i) * 1.0d-06    ! convert from [um] to [m]
        den_kgm3 = den_ddep(i) * 1.0d+03    ! convert from [g/cm^3] to [kg/m^3]
        call GetDep_V( tk, rhoa, xlm, amu, wstar, ustar, ra, dgn_m, xls_ddep(i), den_kgm3, vdep )
!       vdep(:) = min( vdep(:), 10.0d+00 )  ! cap at 10 [m/s] = 1000 [cm/s]; should have no effect 
         !print*,"get dpv2=",i,n;call flush(6)
        vddep_aero(i,1) = vdep(1) ! for deposition of number [m/s]
        vddep_aero(i,2) = vdep(2) ! for deposition of mass   [m/s]
	!print*,"get dpv=",i,vddep_aero(i,1),vddep_aero(i,2);call flush(6)
      enddo

      return
      end subroutine Get_Aero_Depv


!>----------------------------------------------------------------------------------------------------------------------
!!     calculate deposition velocity for aitken, accumulation, and coarse modes.
!!     reference: binkowski f. s., and u. shankar, the regional particulate model 
!!     1. model description and preliminary results. j. geophys. res., 100, d12, 26191-26209, 1995.
!!
!!     12-18-06, dlw: derived from the cmaq routine aero_depv.f originally coded by f. s. binkowski. :: 
!!----------------------------------------------------------------------------------------------------------------------
      subroutine GetDep_V( blkta, blkdens, xlm, amu, blkwstar, blkustar, blkra, dgacc, xxlsgac, pdensac, vdep )
      implicit none

      ! arguments. 

      real(8) :: blkta    !< air temperature [k]
      real(8) :: blkdens  !< air density  [kg/m^3]
      real(8) :: xlm      !< atmospheric mean free path [m]
      real(8) :: amu      !< atmospheric dynamic viscosity [kg/(m s)]
      real(8) :: blkwstar !< convective velocity scale [m/s]
      real(8) :: blkustar !< friction velocity [m/s]
      real(8) :: blkra    !< aerodynamic resistance [s/m]
      real(8) :: dgacc    !< geo. mean diamter [m]
      real(8) :: xxlsgac  !< ln(geo. std. dev.) [1]
      real(8) :: pdensac  !< average particle density [kg/m^3]
      real(8) :: vdep(2)  !< deposition  velocity [ m/s ]

      ! local variables.

      real(8) :: knacc    ! modal knudsen [1]
      real(8) :: dchat0a  ! modal particle diffusivity for number [m^2/s]
      real(8) :: dchat3a  ! modal particle diffusivity for mass   [m^2/s]
      real(8) :: vghat0a  ! modal sedimentation velocity for number [m/s]
      real(8) :: vghat3a  ! modal sedimentation velocity for mass   [m/s]
      real(8) :: dconst1, dconst1a
      real(8) :: dconst2, dconst3a
      real(8) :: sc0a     ! schmidt numbers for number
      real(8) :: sc3a     ! schmidt numbers for 3rd moment
      real(8) :: st0a     ! stokes numbers for number
      real(8) :: st3a     ! stokes numbers for 3rd moment
      real(8) :: rd0a     ! canopy resistance for number
      real(8) :: rd3a     ! canopy resisteance for 3rd moment
      real(8) :: utscale  ! scratch function of ustar and wstar
      real(8) :: nu       ! kinematic viscosity [ m**2 s**-1 ]
      real(8) :: ustfac   ! scratch function of ustar, nu, and grav

      real(8), parameter :: bhat     = 1.246d+00             ! constant from cunningham slip correction
      real(8), parameter :: pi       = 3.141592653589793d+00        
      real(8), parameter :: pi6      = pi / 6.0d+00
      real(8), parameter :: threepi  = 3.0d+00 * pi
      real(8), parameter :: one3     = 1.0d+00 / 3.0d+00
      real(8), parameter :: two3     = 2.0d+00 / 3.0d+00
      real(8), parameter :: avo      = 6.0221367d+23         ! avogadro's constant  [1/mol]
      real(8), parameter :: rgasuniv = 8.314510d+00          ! universal gas const  [j/mol/k]
      real(8), parameter :: boltz    = rgasuniv / avo        ! boltzmann's constant [j/k]
      !----------------------------------------------------------------------------------------------------------------
      ! value is the mean of polar and equatorial values. crc handbook (76th ed) page 14-6. (fsb)
      !----------------------------------------------------------------------------------------------------------------
      real(8), parameter :: grav        = 9.80622d+00        ! mean gravitational accel [m/s^2]
      real(8), parameter :: dgacc_max   = 1.0d-06            ! 1.0 um, min. value for elimination of impaction term [m]
      real(8), parameter :: xxlsgac_max = 0.6931472d+00      ! ln(2),  min. value for elimination of impaction term [1]

      ! scratch variables for standard deviations.

      real(8) :: l2sgac           
      real(8) :: eac1             
      real(8) :: esac04    
      real(8) :: esac08    
      real(8) :: esac16    
      real(8) :: esac20     
      real(8) :: esac28       
      real(8) :: esac32    
      real(8) :: esac64    

         !print*,"1=";call flush(6)

      knacc  = 2.0d+00 * xlm / dgacc
      l2sgac = xxlsgac * xxlsgac
      eac1   = exp( 0.125d+00 * l2sgac )
      esac04 = eac1**4
      esac08 = esac04 * esac04
      esac16 = esac08 * esac08
      esac20 = esac16 * esac04
      esac28 = esac20 * esac08
      esac32 = esac16 * esac16
      esac64 = esac32 * esac32

      dconst1  = boltz * blkta / ( threepi * amu )
      dconst1a = dconst1 / dgacc
      dconst2  = grav / ( 18.0d+00 * amu )
      dconst3a = dconst2 * pdensac * dgacc*dgacc

      dchat0a = dconst1a * ( esac04  + bhat * knacc * esac16  )
      dchat3a = dconst1a * ( ( 1.0d+00 / esac20 ) + bhat * knacc / esac32 )
      vghat0a = dconst3a * ( esac16  + bhat * knacc * esac04  )
      vghat3a = dconst3a * ( esac64  + bhat * knacc * esac28  )

      nu      = amu / blkdens
      ustfac  = blkustar * blkustar / ( grav * nu )
      utscale = blkustar + 0.24d+00 * blkwstar * blkwstar / blkustar
         !print*,"3=";call flush(6)

      sc0a = nu / dchat0a
      st0a = max ( vghat0a * ustfac, 0.01d+00 )
      if( dgacc .lt. dgacc_max ) then                        ! not a coarse mode ...
        rd0a = 1.0d+00 / ( utscale * ( sc0a**( -two3 ) + 10.0**( -3.0d+00 / st0a ) ) )
      else
        rd0a = 1.0d+00 / ( utscale * ( sc0a**( -two3 ) ) )   ! eliminate impaction term for coarse modes as in cmaq.
      endif
      vdep(1) = vghat0a + 1.0d+00 / ( blkra + rd0a + rd0a * blkra * vghat0a )  ! for deposition of number.

      sc3a = nu / dchat3a
      st3a = max( vghat3a * ustfac , 0.01d+00 )
      if( dgacc .lt. dgacc_max  ) then                       ! not a coarse mode ...
        rd3a = 1.0d+00 / ( utscale * ( sc3a**( -two3 ) + 10.0**( -3.0d+00 / st3a ) ) )
      else
        rd3a = 1.0d+00 / ( utscale * ( sc3a**( -two3 ) ) )   ! eliminate impaction term for coarse modes as in cmaq.
      endif
      vdep(2) = vghat3a + 1.0d+00 / ( blkra + rd3a + rd3a * blkra * vghat3a )  ! for deposition of mass.

         !print*,"4=";call flush(6)
      return
      end subroutine GetDep_V
      
      
   end module Aero_Depv
