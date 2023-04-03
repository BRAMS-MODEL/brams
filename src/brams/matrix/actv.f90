   module Aero_Actv
!      use amp_aerosol, only: nactv
      use memMatrix,  only: &
                           nlays, & !intent()
                           aunit1, &   !intent()
                           nmodes   !intent()
!-------------------------------------------------------------------------------------------------------------------------
!@sum     The array NACTV(X,Y,Z,I) contains current values of the number of aerosol particles 
!@+       activated in clouds for each mode I for use outside of the MATRIX microphysical module.
!@+       Values in NACTV are saved in subr. MATRIX at each time step. 
!
!@auth    Susanne Bauer/Doug Wright
!
!-------------------------------------------------------------------------------------------------------------------------
      
      real(8), parameter :: dens_sulf = 1.77d+03    ! [kg/m^3] nh42so4
      real(8), parameter :: dens_bcar = 1.70d+03    ! [kg/m^3] ghan et al. (2001) - mirage
      real(8), parameter :: dens_ocar = 1.00d+03    ! [kg/m^3] ghan et al. (2001) - mirage
      real(8), parameter :: dens_dust = 2.60d+03    ! [kg/m^3] ghan et al. (2001) - mirage
      real(8), parameter :: dens_seas = 2.165d+03   ! [kg/m^3] nacl, ghan et al. (2001) used 1.90d+03

      CONTAINS


!>----------------------------------------------------------------------------------------------------------------------
!!    12-12-06, DLW: Routine to set up the call to subr. ACTFRAC_MAT to calculate the 
!!                   activated fraction of the number and mass concentrations, 
!!                    as well as the number and mass concentrations activated 
!!                    for each of NMODEX modes. The minimum dry radius for activation 
!!                    for each mode is also returned. 
!!
!!     Each mode is assumed to potentially contains 5 chemical species:
!!         (1) sulfate 
!!         (2) BC 
!!         (3) OC
!!         (4) mineral dust
!!         (5) sea salt 
!!
!!     The aerosol activation parameterizations are described in 
!!
!!         1. Abdul-Razzak et al.   1998, JGR, vol.103, p.6123-6131.
!!         2. Abdul-Razzak and Ghan 2000, JGR, vol.105, p.6837-6844. 
!!
!!     and values for many of the required parameters were taken from 
!!
!!         3. Ghan et al. 2001, JGR vol 106, p.5295-5316.
!!
!!     With the density of sea salt set to the value used in ref. 3 (1900 kg/m^3), this routine 
!!     yields values for the hygroscopicity parameters Bi in agreement with ref. 3. 
!!----------------------------------------------------------------------------------------------------------------------
      subroutine GetActFrac(nmodex,xnap,xmap5,rg,sigmag,tkelvin,ptot,wupdraft, &
                            ac,fracactn,fracactm,nact,mact)
      IMPLICIT NONE

      integer, parameter :: ncomps = 5

      ! arguments.
      
      integer :: nmodex               !< number of modes [1]      
      real(8) :: xnap(nmodex)         !< number concentration for each mode [#/m^3]
      real(8) :: xmap5(nmodex,ncomps) !< mass concentration of each of the 5 species for each mode [ug/m^3]
      real(8) :: rg(nmodex)           !< geometric mean dry radius for each mode [um]
      real(8) :: sigmag(nmodex)       !< geometric standard deviation for each mode [um]
      real(8) :: tkelvin              !< absolute temperature [k]
      real(8) :: ptot                 !< ambient pressure [pa]
      real(8) :: wupdraft             !< updraft velocity [m/s]
      real(8) :: ac(nmodex)           !< minimum dry radius for activation for each mode [um]
      real(8) :: fracactn(nmodex)     !< activating fraction of number conc. for each mode [1]
      real(8) :: fracactm(nmodex)     !< activating fraction of mass   conc. for each mode [1]
      real(8) :: nact(nmodex)         !< activating number concentration for each mode [#/m^3]
      real(8) :: mact(nmodex)         !< activating mass   concentration for each mode [ug/m^3]

      ! local variables. 

      integer :: i, j                 ! loop counters 
      real(8) :: xmap(nmodex)         ! total mass concentration for each mode [ug/m^3]
      real(8) :: bibar(nmodex)        ! hygroscopicity parameter for each mode [1]

      ! variables for mode-average hygroscopicity parameters. 
      
      real(8)       :: xr  (nmodex,ncomps)  ! mass fraction for component j in mode i [1]   
      real(8), save :: xnu (ncomps)         ! # of ions formed per formula unit solute for component j in mode i [1]
      real(8), save :: xphi(ncomps)         ! osmotic coefficient for component j in mode i [1]
      real(8), save :: xmw (ncomps)         ! molecular weight for component j in mode i [kg/mol]
      real(8), save :: xrho(ncomps)         ! density of component j in mode i [kg/m^3]
      real(8), save :: xeps(ncomps)         ! soluble fraction of component j in mode i [1]

      real(8) :: sumnumer, sumdenom         ! scratch variables 

      real(8), parameter :: nion_sulf = 3.00d+00    ! [1]
      real(8), parameter :: nion_bcar = 1.00d+00    ! [1]
      real(8), parameter :: nion_ocar = 1.00d+00    ! [1]
      real(8), parameter :: nion_dust = 2.30d+00    ! [1]
      real(8), parameter :: nion_seas = 2.00d+00    ! [1] nacl

      real(8), parameter :: xphi_sulf = 0.70d+00    ! [1]
      real(8), parameter :: xphi_bcar = 1.00d+00    ! [1]
      real(8), parameter :: xphi_ocar = 1.00d+00    ! [1]
      real(8), parameter :: xphi_dust = 1.00d+00    ! [1]
      real(8), parameter :: xphi_seas = 1.00d+00    ! [1] nacl

      real(8), parameter :: molw_sulf = 132.0d-03   ! [kg/mol]
      real(8), parameter :: molw_bcar = 100.0d-03   ! [kg/mol]
      real(8), parameter :: molw_ocar = 100.0d-03   ! [kg/mol]
      real(8), parameter :: molw_dust = 100.0d-03   ! [kg/mol]
      real(8), parameter :: molw_seas = 58.44d-03   ! [kg/m^3] nacl

      real(8), parameter :: xeps_sulf = 1.00d+00    ! [1]
      real(8), parameter :: xeps_bcar = 1.67d-06    ! [1]
      real(8), parameter :: xeps_ocar = 0.78d+00    ! [1]
      real(8), parameter :: xeps_dust = 0.13d+00    ! [1]
      real(8), parameter :: xeps_seas = 1.00d+00    ! [1] nacl

      real(8), parameter :: wmolmass = 18.01528d-03 ! molar mass of h2o     [kg/mol]
      real(8), parameter :: denh2o   =  1.00d+03    ! density of water [kg/m^3]

      logical, save :: firstime = .true.
      
      if( firstime ) then
        firstime = .false.
        xnu (1) = nion_sulf
        xnu (2) = nion_bcar
        xnu (3) = nion_ocar
        xnu (4) = nion_dust
        xnu (5) = nion_seas
        xphi(1) = xphi_sulf
        xphi(2) = xphi_bcar
        xphi(3) = xphi_ocar
        xphi(4) = xphi_dust
        xphi(5) = xphi_seas
        xmw (1) = molw_sulf
        xmw (2) = molw_bcar
        xmw (3) = molw_ocar
        xmw (4) = molw_dust
        xmw (5) = molw_seas
        xrho(1) = dens_sulf
        xrho(2) = dens_bcar
        xrho(3) = dens_ocar
        xrho(4) = dens_dust
        xrho(5) = dens_seas
        xeps(1) = xeps_sulf
        xeps(2) = xeps_bcar
        xeps(3) = xeps_ocar
        xeps(4) = xeps_dust
        xeps(5) = xeps_seas
      endif

      !--------------------------------------------------------------------------------------------------------------
      ! calculate the mass fraction component j for each mode i. 
      !--------------------------------------------------------------------------------------------------------------
      do i=1, nmodex
        xmap(i) = 0.0d+00
        do j=1, ncomps
          xmap(i) = xmap(i) + xmap5(i,j)
        enddo
        xr(i,:) = xmap5(i,:) / max( xmap(i), 1.0d-30 )   
        ! write(*,'(i4,5f12.6)') i,xr(i,:)
      enddo

      !--------------------------------------------------------------------------------------------------------------
      ! calculate the hygroscopicity parameter for each mode. 
      !--------------------------------------------------------------------------------------------------------------
      do i=1, nmodex
        sumnumer = 0.0d+00
        sumdenom = 0.0d+00
        do j=1, ncomps
          sumnumer = sumnumer + xr(i,j)*xnu(j)*xphi(j)*xeps(j)/xmw(j)     ! [mol/kg] 
          sumdenom = sumdenom + xr(i,j)/xrho(j)                           ! [m^3/kg] 
!write (*,FMT='(A,2(I3.3,1X),5(E20.8,1X))') 'LFR=DBG GetActFrac:',i,j, &
!        sumdenom,xmap5(i,j),xmap(i),xrho(j),0.0
        enddo

        bibar(i) = ( wmolmass*sumnumer ) / ( denh2o*sumdenom )            ! [1] 
      enddo
      ! write(*,'(8d15.6)') bibar(:)

      !--------------------------------------------------------------------------------------------------------------
      ! calculate the droplet activation parameters for each mode. 
      !--------------------------------------------------------------------------------------------------------------
      call ActFrac_Mat(nmodex,xnap,xmap,rg,sigmag,bibar,tkelvin,ptot,wupdraft, &
                       ac,fracactn,fracactm,nact,mact)

      do i=1, nmodex
        if(xnap(i) .lt. 1.0d-06 ) fracactn(i) = 1.0d-30
      enddo

      return 
      end subroutine GetActFrac


!>----------------------------------------------------------------------------------------------------------------------
!!     12-12-06, DLW: Routine to calculate the activated fraction of the number 
!!                    and mass concentrations, as well as the number and mass 
!!                    concentrations activated for each of NMODEX modes. The 
!!                    minimum dry radius for activation for each mode is also returned. 
!!
!!     The aerosol activation parameterizations are described in 
!!
!!         1. Abdul-Razzak et al.   1998, JGR, vol.103, p.6123-6131.
!!         2. Abdul-Razzak and Ghan 2000, JGR, vol.105, p.6837-6844. 
!! 
!!     This routine is for the multiple-aerosol type parameterization. 
!!----------------------------------------------------------------------------------------------------------------------
      subroutine ActFrac_Mat(nmodex,xnap,xmap,rg,sigmag,bibar,tkelvin,ptot,wupdraft, &
                             ac,fracactn,fracactm,nact,mact)
!!LFR  cccONLY FOR GCM      USE DOMAIN_DECOMP_ATM,only: am_i_root
      IMPLICIT NONE

      ! Arguments.
      
      integer :: nmodex            !< number of modes [1]      
      real(8) :: xnap(nmodex)      !< number concentration for each mode [#/m^3]
      real(8) :: xmap(nmodex)      !< mass   concentration for each mode [ug/m^3]
      real(8) :: rg(nmodex)        !< geometric mean radius for each mode [um]
      real(8) :: sigmag(nmodex)    !< geometric standard deviation for each mode [um]
      real(8) :: bibar(nmodex)     !< hygroscopicity parameter for each mode [1]
      real(8) :: tkelvin           !< absolute temperature [k]
      real(8) :: ptot              !< ambient pressure [pa]
      real(8) :: wupdraft          !< updraft velocity [m/s]
      real(8) :: ac(nmodex)        !< minimum dry radius for activation for each mode [um]
      real(8) :: fracactn(nmodex)  !< activating fraction of number conc. for each mode [1]
      real(8) :: fracactm(nmodex)  !< activating fraction of mass   conc. for each mode [1]
      real(8) :: nact(nmodex)      !< activating number concentration for each mode [#/m^3]
      real(8) :: mact(nmodex)      !< activating mass   concentration for each mode [ug/m^3]

      ! parameters.
      
      real(8), parameter :: pi            = 3.141592653589793d+00
      real(8), parameter :: twopi         = 2.0d+00 * pi
      real(8), parameter :: sqrt2         = 1.414213562d+00
      real(8), parameter :: threesqrt2by2 = 1.5d+00 * sqrt2

      real(8), parameter :: avgnum   = 6.0221367d+23       ! [1/mol]
      real(8), parameter :: rgasjmol = 8.31451d+00         ! [j/mol/k]
      real(8), parameter :: wmolmass = 18.01528d-03        ! molar mass of h2o     [kg/mol]
      real(8), parameter :: amolmass = 28.966d-03          ! molar mass of air     [kg/mol]
      real(8), parameter :: asmolmss = 132.1406d-03        ! molar mass of nh42so4 [kg/mol]
      real(8), parameter :: denh2o   = 1.00d+03            ! density of water [kg/m^3]
      real(8), parameter :: denamsul = 1.77d+03            ! density of pure ammonium sulfate [kg/m^3]
      real(8), parameter :: xnuamsul = 3.00d+00            ! # of ions formed when the salt is dissolved in water [1]
      real(8), parameter :: phiamsul = 1.000d+00           ! osmotic coefficient value in a-r 1998. [1] 
      real(8), parameter :: gravity  = 9.81d+00            ! grav. accel. at the earth's surface [m/s/s] 
      real(8), parameter :: heatvap  = 40.66d+03/wmolmass  ! latent heat of vap. for water and tnbp [j/kg] 
      real(8), parameter :: cpair    = 1006.0d+00          ! heat capacity of air [j/kg/k] 
      real(8), parameter :: t0dij    = 273.15d+00          ! reference temp. for dv [k] 
      real(8), parameter :: p0dij    = 101325.0d+00        ! reference pressure for dv [pa] 
      real(8), parameter :: dijh2o0  = 0.211d-04           ! reference value of dv [m^2/s] (p&k,2nd ed., p.503)
      !----------------------------------------------------------------------------------------------------------------    
      ! real(8), parameter :: t0dij    = 283.15d+00          ! reference temp. for dv [k] 
      ! real(8), parameter :: p0dij    = 80000.0d+00         ! reference pressure for dv [pa] 
      ! real(8), parameter :: dijh2o0  = 0.300d-04           ! reference value of dv [m^2/s] (p&k,2nd ed., p.503)
      !----------------------------------------------------------------------------------------------------------------
      real(8), parameter :: deltav   = 1.096d-07           ! vapor jump length [m]  
      real(8), parameter :: deltat   = 2.160d-07           ! thermal jump length [m]  
      real(8), parameter :: alphac   = 1.000d+00           ! condensation mass accommodation coefficient [1]  
      real(8), parameter :: alphat   = 0.960d+00           ! thermal accommodation coefficient [1]  

      ! local variables. 

      integer            :: i                              ! loop counter 
      real(8)            :: dv                             ! diffusion coefficient for water [m^2/s] 
      real(8)            :: dvprime                        ! modified diffusion coefficient for water [m^2/s] 
      real(8)            :: dumw, duma                     ! scratch variables [s/m] 
      real(8)            :: wpe                            ! saturation vapor pressure of water [pa]  
      real(8)            :: surten                         ! surface tension of air-water interface [j/m^2] 
      real(8)            :: xka                            ! thermal conductivity of air [j/m/s/k]  
      real(8)            :: xkaprime                       ! modified thermal conductivity of air [j/m/s/k]  
      real(8)            :: eta(nmodex)                    ! model parameter [1]  
      real(8)            :: zeta                           ! model parameter [1]  
      real(8)            :: xlogsigm(nmodex)               ! ln(sigmag) [1]   
      real(8)            :: a                              ! [m]
      real(8)            :: g                              ! [m^2/s]   
      real(8)            :: rdrp                           ! [m]   
      real(8)            :: f1                             ! [1]   
      real(8)            :: f2                             ! [1]
      real(8)            :: alpha                          ! [1/m]
      real(8)            :: gamma                          ! [m^3/kg]   
      real(8)            :: sm(nmodex)                     ! [1]   
      real(8)            :: dum                            ! [1/m]    
      real(8)            :: u                              ! argument to error function [1]
      real(8)            :: erf                            ! error function [1], but not declared in an f90 module 
      real(8)            :: smax                           ! maximum supersaturation [1]
      
!----------------------------------------------------------------------------------------------------------------------
!     rdrp is the radius value used in eqs.(17) & (18) and was adjusted to yield eta and zeta 
!     values close to those given in a-z et al. 1998 figure 5. 
!----------------------------------------------------------------------------------------------------------------------
      rdrp = 0.105d-06   ! [m] tuned to approximate the results in figures 1-5 in a-z et al. 1998.  
!----------------------------------------------------------------------------------------------------------------------
!     these variables are common to all modes and need only be computed once. 
!----------------------------------------------------------------------------------------------------------------------
      dv = dijh2o0*(p0dij/ptot)*(tkelvin/t0dij)**1.94d+00                 ! [m^2/s] (p&k,2nd ed., p.503)
      surten = 76.10d-03 - 0.155d-03 * (tkelvin-273.15d+00)               ! [j/m^2] 
      wpe = exp( 77.34491296d+00 - 7235.424651d+00/tkelvin - 8.2d+00*log(tkelvin) + tkelvin*5.7113d-03 )  ! [pa] 
      dumw = sqrt(twopi*wmolmass/rgasjmol/tkelvin)                        ! [s/m] 
      dvprime = dv / ( (rdrp/(rdrp+deltav)) + (dv*dumw/(rdrp*alphac)) )   ! [m^2/s] - eq. (17) 
      xka = (5.69d+00+0.017d+00*(tkelvin-273.15d+00))*418.4d-05           ! [j/m/s/k] (0.0238 j/m/s/k at 273.15 k)
      duma = sqrt(twopi*amolmass/rgasjmol/tkelvin)                        ! [s/m]
      xkaprime = xka / ( ( rdrp/(rdrp+deltat) ) + ( xka*duma/(rdrp*alphat*denh2o*cpair) ) )   ! [j/m/s/k]
      g = 1.0d+00 / ( (denh2o*rgasjmol*tkelvin) / (wpe*dvprime*wmolmass) &
                      + ( (heatvap*denh2o) / (xkaprime*tkelvin) ) &
                      * ( (heatvap*wmolmass) / (rgasjmol*tkelvin) - 1.0d+00 ) )               ! [m^2/s]
      a = (2.0d+00*surten*wmolmass)/(denh2o*rgasjmol*tkelvin)                                 ! [m] 
      alpha = (gravity/(rgasjmol*tkelvin))*((wmolmass*heatvap)/(cpair*tkelvin) - amolmass)    ! [1/m] 
      gamma = (rgasjmol*tkelvin)/(wpe*wmolmass) &
            + (wmolmass*heatvap*heatvap)/(cpair*ptot*amolmass*tkelvin)                        ! [m^3/kg]
      dum = sqrt(alpha*wupdraft/g)                  ! [1/m] 
      zeta = 2.d+00*a*dum/3.d+00                    ! [1] 
      !write(*,*) 'DBG:LFR dum,alpha,wupdraft,g: ',dum,alpha,wupdraft,g
      !write(*,*) 'DBG:LFR gravity,rgasjmol,tkelvin,wmolmass,heatvap,cpair,amolmass: ', &
      !gravity,rgasjmol,tkelvin,wmolmass,heatvap,cpair,amolmass
!----------------------------------------------------------------------------------------------------------------
      ! write(1,'(a27,4d15.5)')'surten,wpe,a            =',surten,wpe,a
      ! write(1,'(a27,4d15.5)')'xka,xkaprime,dv,dvprime =',xka,xkaprime,dv,dvprime
      ! write(1,'(a27,4d15.5)')'alpha,gamma,g, zeta     =',alpha,gamma,g,zeta
!----------------------------------------------------------------------------------------------------------------------
!     these variables must be computed for each mode. 
!----------------------------------------------------------------------------------------------------------------------
      xlogsigm(:) = log(sigmag(:))                                                    ! [1] 
      smax = 0.0d+00      
                                                          ! [1]
      do i=1, nmodex

        sm(i) = ( 2.0d+00/sqrt(bibar(i)) ) * ( a/(3.0d-06*rg(i)) )**1.5d+00           ! [1] 
        eta(i) = dum**3 / (twopi*denh2o*gamma*xnap(i))                                ! [1] 

        !--------------------------------------------------------------------------------------------------------------
        ! write(1,'(a27,i4,4d15.5)')'i,eta(i),sm(i) =',i,eta(i),sm(i)
        !--------------------------------------------------------------------------------------------------------------
        f1 = 0.5d+00 * exp(2.50d+00 * xlogsigm(i)**2)                                 ! [1] 
        f2 = 1.0d+00 +     0.25d+00 * xlogsigm(i)                                     ! [1] 
        smax = smax + (   f1*(  zeta  / eta(i)              )**1.50d+00 &
                        + f2*(sm(i)**2/(eta(i)+3.0d+00*zeta))**0.75d+00 ) / sm(i)**2  ! [1] - eq. (6)
      enddo 
      smax = 1.0d+00 / sqrt(smax)                                                     ! [1]
      do i=1, nmodex

        ac(i)       = rg(i) * ( sm(i) / smax )**0.66666666666666667d+00               ! [um]

        u           = log(ac(i)/rg(i)) / ( sqrt2 * xlogsigm(i) )                      ! [1]
        !write(*,*) 'DBG:LFR i,rg,xlog= ',i,u,ac(i),rg(i),( sqrt2 * xlogsigm(i) )       
        fracactn(i) = 0.5d+00 * (1.0d+00 - erf(u))                                    ! [1]
        fracactm(i) = 0.5d+00 * (1.0d+00 - erf(u - threesqrt2by2*xlogsigm(i) ) )      ! [1]
        nact(i)     = fracactn(i) * xnap(i)                                           ! [#/m^3]
        mact(i)     = fracactm(i) * xmap(i)                                           ! [ug/m^3]
        !--------------------------------------------------------------------------------------------------------------
        !if(fracactn(i) .gt. 0.9999999d+00 ) then
        !   write(*,*)i,ac(i),u,fracactn(i),xnap(i)
        !  print*,' susa',i,ac(i),u,fracactn(i),xnap(i)
        ! stop
        ! endif
        !--------------------------------------------------------------------------------------------------------------
      enddo 

      return
      end subroutine ActFrac_Mat


!>-----------------------------------------------------------------------------------------------------------------------
!!     see numerical recipes, w. press et al., 2nd edition.
!!-----------------------------------------------------------------------------------------------------------------------
      subroutine GcfMatrix(gammcf,a,x,gln)

      implicit none
      integer, parameter :: itmax=10000
      real(8), parameter :: eps=3.0d-07
      real(8), parameter :: fpmin=1.0d-30
      real(8) :: a,gammcf,gln,x
      integer :: i
      real(8) :: an,b,c,d,del,h
      !real(8) :: gammln   ! function names not declared in an f90 module 
      gln=gammln(a)
      b=x+1.0d+00-a
      c=1.0d+00/fpmin
      d=1.0d+00/b
      h=d
      do i=1,itmax
        an=-i*(i-a)
        b=b+2.0d+00
        d=an*d+b
        if(abs(d).lt.fpmin)d=fpmin
        c=b+an/c
        if(abs(c).lt.fpmin)c=fpmin
        d=1.0d+00/d
        del=d*c
        h=h*del
        if(abs(del-1.0d+00).lt.eps)goto 1
      enddo
      write(*,*)'AERO_ACTV: SUBROUTINE GCF: A TOO LARGE, ITMAX TOO SMALL', gammcf,a,x,gln
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      end subroutine GcfMatrix


!>-----------------------------------------------------------------------------------------------------------------------
!!     see numerical recipes, w. press et al., 2nd edition.
!!-----------------------------------------------------------------------------------------------------------------------
      subroutine Gser(gamser,a,x,gln)

      implicit none
      integer, parameter :: itmax=10000  ! was itmax=100   in press et al. 
      real(8), parameter :: eps=3.0d-09  ! was eps=3.0d-07 in press et al.
      real(8) :: a,gamser,gln,x
      integer :: n
      real(8) :: ap,del,sum
      !real(8) :: gammln   ! function names not declared in an f90 module 
      gln=gammln(a)
      if(x.le.0.d+00)then
        if(x.lt.0.)stop 'aero_actv: subroutine gser: x < 0 in gser'
        gamser=0.d+00
        return
      endif
      ap=a
      sum=1.d+00/a
      del=sum
      do n=1,itmax
        ap=ap+1.d+00
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*eps)goto 1
      enddo
      write(*,*)'aero_actv: subroutine gser: a too large, itmax too small'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      end subroutine Gser


!>-----------------------------------------------------------------------------------------------------------------------
!!     see numerical recipes, w. press et al., 2nd edition.
!!-----------------------------------------------------------------------------------------------------------------------
      double precision function GammLn(xx)

      implicit none
      real(8) :: xx
      integer j
      double precision ser,stp,tmp,x,y,cof(6)
      save cof,stp
      data cof,stp/76.18009172947146d0,-86.50532032941677d0,         &
      24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
      -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
      enddo
      gammln=tmp+log(stp*ser/x)
      return
      end function GammLn 


!>-----------------------------------------------------------------------------------------------------------------------
!!     see numerical recipes, w. press et al., 2nd edition.
!!-----------------------------------------------------------------------------------------------------------------------
      double precision function Erf(x)
      implicit none
      real(8) :: x
!u    uses gammp
      !LFR  real(8) :: gammp   ! function names not declared in an f90 module 
      erf = 0.d0
      if(x.lt.0.0d+00)then
        erf=-gammp(0.5d0,x**2)
      else
        erf= gammp(0.5d0,x**2)
      endif
      return
      end function Erf


!>-----------------------------------------------------------------------------------------------------------------------
!!     see numerical recipes, w. press et al., 2nd edition.
!!-----------------------------------------------------------------------------------------------------------------------
      double precision function GammP(a,x)
      implicit none
      real(8) :: a,x
      real(8) :: gammcf,gamser,gln
      if(x.lt.0.0d+00.or.a.le.0.0d+00)then
        write(*,*)'aero_actv: function gammp: bad arguments'
      endif

      if(x.lt.a+1.0d+00)then
        call Gser(gamser,a,x,gln)
        gammp=gamser
      else
        call GcfMatrix(gammcf,a,x,gln)
        gammp=1.0d+00-gammcf
      endif
      return
      end function GammP


      end module Aero_Actv
      
