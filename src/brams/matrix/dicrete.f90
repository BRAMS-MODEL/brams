   module Aero_Discrete
      use memMatrix, only: densp, pi6
      use Aero_Coag,  only: total_coag_coef, brownian_coag_coef
      implicit none
      integer, parameter :: nbins = 200
      integer, parameter :: nspcs =   2
      real(8) :: apdf(nbins,1+nspcs)       ! mode a: number conc. [#/m^3] and mass concs. [ug/m^3] for each bin. 
      real(8) :: bpdf(nbins,1+nspcs)       ! mode b: number conc. [#/m^3] and mass concs. [ug/m^3] for each bin. 
      real(8) :: cpdf(nbins,1+nspcs)       ! mode c: number conc. [#/m^3] and mass concs. [ug/m^3] for each bin. 
      real(8) :: dgrid(nbins)              ! fixed diameter        grid [um]
      real(8) :: vgrid(nbins)              ! fixed volume/particle grid [um^3/particle]
      real(8) :: mgrid(nbins)              ! fixed mass/particle   grid [ug/particle]
      real(8) :: d3grid(nbins)             ! fixed diameter-cubed  grid [um^3/particle]
      real(8) :: dlower(nbins)             ! lower boundary fixed diameter grid [um]
      real(8) :: dupper(nbins)             ! upper boundary fixed diameter grid [um]
      real(8), parameter :: dmin     =  0.001d+00       ! smallest particle diameter of the discrete grid [um]
      real(8), parameter :: dmax     = 20.000d+00       ! largest  particle diameter of the discrete grid [um]
      real(8), save      :: rdmin    =  0.000d+00       ! reciprocal of dmin to optimize coagulation [1/um]
      real(8), save      :: rdlogdsc =  0.000d+00       ! reciprocal of log10 of the grid spacing [1]
      real(8), parameter :: xkb      =  1.3806505d-23   ! [j/k] http://en.wikipedia.org/wiki/boltzmann_constant
      !-----------------------------------------------------------------------------------------------------------------
      ! ispca, ispcb, ispcc: 1=sulf, 2=bcar, 3=ocar, 4=dust, 5=seas
      !-----------------------------------------------------------------------------------------------------------------
      integer, save :: ispca               ! index of the chemical species for mode a: set in aero_init
      integer, save :: ispcb               ! index of the chemical species for mode b: set in aero_init
      integer, save :: ispcc               ! index of the chemical species for mode c: set in aero_init
      integer, parameter :: itcoag = 10    ! number of subdivisions of the time step for integration of coagulation 
      real(8) :: dtcoag = 0.0d+00
      real(8) :: kij_discrete(nbins,nbins) ! fixed coagulation coefficients for the discrete grid [m^3/s/particle]
      real(8) :: frac_lo (nbins,nbins)     ! used in discrete_intracoag [1]
      real(8) :: frac_hi (nbins,nbins)     ! used in discrete_intracoag [1]
      integer :: nlo_grid(nbins,nbins)     ! used in discrete_intracoag [1]
      integer :: nhi_grid(nbins,nbins)     ! used in discrete_intracoag [1]
      

      CONTAINS


      subroutine Discrete_Init(icset,na,dga,sigmaga,nb,dgb,sigmagb,nc,dgc,sigmagc,temp,pres,massa,massb,massc)
!-----------------------------------------------------------------------------------------------------------------------
!@auth    susanne bauer/doug wright
!-----------------------------------------------------------------------------------------------------------------------
      implicit none

      ! arguments.
 
      integer, intent(   in) :: icset                      ! identifies test case [1]
      real(8), intent(   in) :: na,      nb,      nc       ! number concentrations for modes a, b, and c [#/m^3]
      real(8), intent(   in) :: dga,     dgb,     dgc      ! geo. mean diameters   for modes a, b, and c [um]
      real(8), intent(   in) :: sigmaga, sigmagb, sigmagc  ! geo. std. deviations  for modes a, b, and c [1]
      real(8), intent(   in) :: temp                       ! ambient temperature [k]
      real(8), intent(   in) :: pres                       ! ambient pressure    [pa]
      real(8), intent(  out) :: massa,   massb,   massc    ! mass concentrations   for modes a, b, and c [ug/m^3]

      ! local variables.

      integer :: i, j, nhi, nlo
      real(8) :: scale, fa, fb, fc, wa(nbins), wb(nbins), wc(nbins), wtota, wtotb, wtotc, dnewcub
      real(8) :: dminl, dmaxl, etaa, vp_pdf, mfrac(nbins,nspcs)

      ! for the call to the coagulation coefficient routines

      real(4) :: di, dj    ! ambient particle diameters [um]
      real(4) :: betaij    ! coagulation coefficients [m^2/s/particle]
      real(4) :: templ     ! ambient temperature [k]
      real(4) :: presl     ! ambient pressure [pa]

      real(8), parameter :: onethird = 1.0d+00 / 3.0d+00 


      apdf(:,:) = 0.0d+00
      bpdf(:,:) = 0.0d+00
      cpdf(:,:) = 0.0d+00
      
      ! set up discrete grid.

      select case( icset )
      case ( 10 )
        dmaxl = 0.600d+00
        dminl = 0.006d+00
      case ( 11 )
        dmaxl = 1.000d+00
        dminl = 0.001d+00
      case ( 12 )
        dmaxl = dmax
        dminl = dmin
      case ( 13 )
        dmaxl = dmax
        dminl = dmin
      case ( 14 )
        dmaxl = 100.000d+00
        dminl =   0.001d+00
      case ( 15 )
        dmaxl = dmax
        dminl = dmin
      case ( 16 )
        dmaxl = dmax
        dminl = dmin
      case ( 17 )
        dmaxl = 200.000d+00
        dminl =   0.001d+00
      case ( 18 )
        dmaxl = 8.000d+00
        dminl = 0.001d+00
      case ( 19 )
        dmaxl = dmax
        dminl = dmin
      case ( 20 )
        dmaxl = dmax
        dminl = dmin
      case default
        WRITE(*,*)'Should not reach CASE DEFAULT in subr. discrete_init (A).'
        STOP
      end select
      scale    = ( dmaxl / dminl )**(1.0d+00/real(nbins-1))
      rdlogdsc = 1.0d+00 / log10( scale )
      rdmin    = 1.0d+00 / dminl
      ! WRITE(*,'(I4,5D14.5)') ICSET,DMINL,DMAXL,SCALE,RDLOGDSC,RDMIN
      WRITE(36,'(/A,I6,/)') 'ICSET = ', icset
      WRITE(36,*)'I, DLOWER(I) [um], DGRID(I) [um], DUPPER(I) [um], VGRID(I) [um^3/particle], MGRID(I) [ug/particle]'
      do i=1, nbins
        dgrid(i)  = dminl * scale**(i-1)                  ! [um]
        dlower(i) = dgrid(i) / scale**0.5d+00             ! [um]
        dupper(i) = dgrid(i) * scale**0.5d+00             ! [um]
        d3grid(i) = dgrid(i)**3                           ! [um^3/particle]
        vgrid(i)  = pi6 * dgrid(i)**3                     ! [um^3/particle]
        mgrid(i)  = 1.0d-06 * densp * vgrid(i)            ! [ug/particle]
        write(36,'(i4,3f12.7,2d16.7)') i, dlower(i), dgrid(i), dupper(i), vgrid(i), mgrid(i)
      enddo

      wtota = 0.0d+00
      wtotb = 0.0d+00
      wtotc = 0.0d+00
      fa = 0.0d+00
      fb = 0.0d+00
      fc = 0.0d+00
      write(36,'(/A,I6,/)') 'ICSET = ', icset
      WRITE(36,*)'I, DGRID(I) [um], FLN [/m^3/um], W(I) [#/m^3], WTOT [#/m^3]'
      write(36,*)'i, dgrid(i), fa, fb, fc, wa(i), wb(i), wc(i), wtota, wtotb, wtotc'
      select case( icset )
      case ( 10 )
        wa(:) = 0.0d+00
        wb(:) = 0.0d+00
        wc(:) = 0.0d+00
        wa(1) = na
        wtota = sum( wa(:) ) 
        wtotb = sum( wb(:) ) 
        wtotc = sum( wc(:) )
        do i=1, nbins 
          write(36,'(i4,f10.5,9d14.6)') i, dgrid(i), 0d0, 0d0, 0d0, wa(i), wb(i), wc(i), wtota, wtotb, wtotc
        enddo
      case ( 11 )
        vp_pdf = 0.0005236d+00
        wb(:) = 0.0d+00
        wc(:) = 0.0d+00
        wtotb = 0.0d+00
        wtotc = 0.0d+00
        do i=1, nbins
          wa(i) = na * ( (pi6*(dupper(i)**3-dlower(i)**3))/vp_pdf ) * exp( -vgrid(i)/vp_pdf ) 
          wtota = wtota + wa(i) 
          write(36,'(i4,f10.5,9d14.6)') i, dgrid(i), 0d0, 0d0, 0d0, wa(i), wb(i), wc(i), wtota, wtotb, wtotc
        enddo
      case ( 12, 13, 14, 15, 16, 17, 18, 19, 20 )
        do i=1, nbins
          fa = na * fln( dgrid(i), dga, sigmaga )
          fb = nb * fln( dgrid(i), dgb, sigmagb )
          fc = nc * fln( dgrid(i), dgc, sigmagc )
          wa(i) = fa * ( dupper(i) - dlower(i) )
          wb(i) = fb * ( dupper(i) - dlower(i) )
          wc(i) = fc * ( dupper(i) - dlower(i) )
          wtota = wtota + wa(i) 
          wtotb = wtotb + wb(i) 
          wtotc = wtotc + wc(i) 
          write(36,'(i4,f10.5,9d14.6)') i, dgrid(i), fa, fb, fc, wa(i), wb(i), wc(i), wtota, wtotb, wtotc
        enddo
      case default
        WRITE(*,*)'Should not reach CASE DEFAULT in subr. discrete_init (B).'
        STOP
      END SELECT

      wa(:) = wa(:) * ( na / max( wtota, 1.0d-30) )  ! renormalize to the precise input number concentration.
      wb(:) = wb(:) * ( nb / max( wtotb, 1.0d-30) )  ! renormalize to the precise input number concentration.
      wc(:) = wc(:) * ( nc / max( wtotc, 1.0d-30) )  ! renormalize to the precise input number concentration.

      write(36,*)'Check sums on WA(:), WB(:), WC(:)'
      write(36,'(3d22.14)') na,        nb,        nc
      write(36,'(3d22.14)') sum(wa(:)),sum(wb(:)),sum(wc(:))

      ! get total mass concentration summed over all bins. 

      massa = sum( wa(:) * mgrid(:) ) 
      massb = sum( wb(:) * mgrid(:) ) 
      massc = sum( wc(:) * mgrid(:) ) 

      ! initialize the main discrete grids for number and mass concentrations.
      ! write(*,*)ispca,ispcb,ispcc

      apdf(:,1) = wa(:)                                   ! [#/m^3]
      bpdf(:,1) = wb(:)                                   ! [#/m^3]
      cpdf(:,1) = wc(:)                                   ! [#/m^3]
      mfrac(:,:) = 0.0d+00                                ! [1]
      select case( icset )
      case( 10, 11 )                                  ! one-component for all modes 
        do i=1, nbins
          mfrac(i,2) = real(i) / real(nbins)              ! [1]
          mfrac(i,1) = 1.0d+00 - mfrac(i,2)               ! [1]
        enddo      
        do i=1, nbins
          do j=1, ispca
            apdf(i,1+j) = wa(i) * mgrid(i) * mfrac(i,j)   ! [ug/m^3]
            bpdf(i,1+j) = wb(i) * mgrid(i) * mfrac(i,j)   ! [ug/m^3]
            cpdf(i,1+j) = wc(i) * mgrid(i) * mfrac(i,j)   ! [ug/m^3]
          enddo
        enddo
      case( 12 )
        apdf(:,2:nspcs+1) = 0.0d+00                       ! [ug/m^3]
        bpdf(:,2:nspcs+1) = 0.0d+00                       ! [ug/m^3]
        cpdf(:,2:nspcs+1) = 0.0d+00                       ! [ug/m^3]
        apdf(:,2) = wa(:) * mgrid(:)                      ! [ug/m^3] pure sulfate
        bpdf(:,2) = wb(:) * mgrid(:)                      ! [ug/m^3] pure bc
      case( 13, 14, 15, 16, 17, 18, 19, 20 )
        apdf(:,2:nspcs+1) = 0.0d+00                       ! [ug/m^3]
        bpdf(:,2:nspcs+1) = 0.0d+00                       ! [ug/m^3]
        cpdf(:,2:nspcs+1) = 0.0d+00                       ! [ug/m^3]
        apdf(:,2) = wa(:) * mgrid(:)                      ! [ug/m^3] pure sulfate
        bpdf(:,3) = wb(:) * mgrid(:)                      ! [ug/m^3] pure bc
      case default
        WRITE(*,*)'Should not reach CASE DEFAULT in subr. discrete_init (C).'
        stop
      end select

      ! calculate particle partitioning between adjacient grid points to optimize the coagulation integration.

      frac_lo(:,:) = 0.0d+00
      frac_hi(:,:) = 0.0d+00
      nlo_grid(:,:) = 1
      nhi_grid(:,:) = 1
      do i=1, nbins
         do j=1, nbins
            dnewcub = d3grid(i) + d3grid(j)
            nlo  = int( log10( rdmin * dnewcub**onethird ) * rdlogdsc ) + 1
            nlo  = min( nbins-1, nlo )
            nhi  = nlo + 1
            nhi  = min( nbins, nhi )
            nlo_grid(i,j) = nlo
            nhi_grid(i,j) = nhi
            if(     1.000000001d+00*dnewcub .lt. d3grid(nlo) ) then
               frac_lo(i,j) = 1.0d+00
               frac_hi(i,j) = 1.0d+00 - frac_lo(i,j)
               if( nlo .eq. 1 .or. nlo .eq. nbins-1 ) CYCLE
               if( nhi .eq. 2 .or. nhi .eq. nbins   ) CYCLE
               write(*,90) nlo, nhi, d3grid(nlo), dnewcub, d3grid(nhi)
               stop
            elseif( 0.999999999d+00*dnewcub .gt. d3grid(nhi) ) then
               frac_lo(i,j) = 0.0d+00
               frac_hi(i,j) = 1.0d+00 - frac_lo(i,j)
               if( nlo .eq. 1 .or. nlo .eq. nbins-1 ) CYCLE
               if( nhi .eq. 2 .or. nhi .eq. nbins   ) CYCLE
               write(*,90) nlo, nhi, d3grid(nlo), dnewcub, d3grid(nhi)
               stop
            else 
               frac_lo(i,j) = ( d3grid(nhi) - dnewcub ) / ( d3grid(nhi) - d3grid(nlo) )
               frac_hi(i,j) = 1.0d+00 - frac_lo(i,j)
            endif
        ! WRITE(36,'(4I5,2D15.5)') I,J,NLO,NHI,FRAC_LO(I,J), FRAC_HI(I,J)
         ENDDO
      ENDDO

      ! calculate the fixed grid coagulation coefficients.

      select case( icset )
      case ( 10, 11 )    
        !---------------------------------------------------------------------------------------------------------------
        ! etaa: dynamic viscosity of air [kg/m/s], jacobson, 1999, eq.(4.55)
        !---------------------------------------------------------------------------------------------------------------
        etaa = 1.832d-05*(416.16d+00/(temp+120.0d+00))*(temp/296.16d+00)**1.5d+00
        kij_discrete(:,:) = 8.0d+00 * xkb * temp / ( 3.0d+00 * etaa )
      case( 12, 13, 14, 15, 16, 17, 18, 19, 20 )
        kij_discrete(:,:) = 0.0d+00
        templ = real( temp )                                         ! convert to single precision               
        presl = real( pres )                                         ! convert to single precision
        do i=1, nbins
            do j=i, nbins
               di = real( dgrid(i) )                                      ! convert to single precision
               dj = real( dgrid(j) )                                      ! convert to single precision
               !call brownian_coag_coef( di, dj, templ, presl, betaij )    ! all variables single precision
               call Total_Coag_Coef   ( di, dj, templ, presl, betaij )    ! all variables single precision
               kij_discrete(i,j) = real( betaij )                         ! convert to double precision
               kij_discrete(j,i) = real( betaij )                         ! convert to double precision
            enddo
        enddo
      case default
        WRITE(*,*)'Should not reach CASE DEFAULT in subr. discrete_init (D).'
        stop
        !---------------------------------------------------------------------------------------------------------------
        ! etaa = 1.832d-05*(416.16d+00/(temp+120.0d+00))*(temp/296.16d+00)**1.5d+00
        ! kij_discrete(:,:) = 8.0d+00 * xkb * temp / ( 3.0d+00 * etaa )
        ! kij_discrete(:,:) = 1.0d-14
        !---------------------------------------------------------------------------------------------------------------
      end select

      write(37,'(f8.2,11d13.5)') 0d0, &
                               1.0d-06*sum(apdf(:,1)), 1.0d-06*sum(bpdf(:,1)), 1.0d-06*sum(cpdf(:,1)), &
                               sum( apdf(:,2) ),  sum( bpdf(:,2) ),  sum( cpdf(:,2) ), &
                               sum( apdf(:,3) ),  sum( bpdf(:,3) ),  sum( cpdf(:,3) ), &
                               sum( apdf(:,2) ) + sum( bpdf(:,2) ) + sum( cpdf(:,2) ), &
                               sum( apdf(:,3) ) + sum( bpdf(:,3) ) + sum( cpdf(:,3) )

      call Discrete_Out(38,icset,0d0)      
      
90    FORMAT('DISCRETE_INIT: Bad NLO or NHI: NLO, NHI, D3GRID(NLO), DNEWCUB, D3GRID(NHI) = ',/,2I10,3D22.14)
      RETURN
      END SUBROUTINE Discrete_Init


      subroutine Discrete( icset, timeh, tstep )
!-----------------------------------------------------------------------------------------------------------------------
!     10-30-06, dlw: discrete model of the pdf used to obtain results for evaluation of matrix. 
!-----------------------------------------------------------------------------------------------------------------------

      ! arguments.
 
      integer             :: icset                 ! label for set of initial conditions [1] 
      real(8), intent(in) :: timeh                 ! model time [h]
      real(8), intent(in) :: tstep                 ! model physics time step [s]

      ! local variables.

      integer :: it
      logical, save      :: firstime           = .true.
      logical, parameter :: discrete_coag_flag = .true.

!----------------------------------------------------------------------------------------------------------------------
!     begin execution.
!----------------------------------------------------------------------------------------------------------------------
      if ( firstime ) then
        firstime = .false.
        dtcoag = tstep / real( itcoag )
      endif                           


      if( discrete_coag_flag ) then
        select case( icset )
        case ( 10, 11, 18 )
          do it=1, itcoag
            call discrete_intracoag    ( apdf )
          enddo
        case ( 12, 20 )
          do it=1, itcoag
            call discrete_intracoag    ( apdf )
            call discrete_intracoag    ( bpdf )
            call discrete_intercoag_ab ( apdf, bpdf )
          enddo
        case ( 13, 14, 15, 16, 17, 19 )
          do it=1, itcoag
            call discrete_intracoag    ( apdf )
            call discrete_intracoag    ( bpdf )
            call discrete_intracoag    ( cpdf )
            call discrete_intercoag_ab ( apdf, cpdf )
            call discrete_intercoag_ab ( bpdf, cpdf )
            call discrete_intercoag_abc( apdf, bpdf, cpdf )
          enddo
        case default
          WRITE(*,*)'Should not reach CASE DEFAULT in subr. discrete.'
          stop
        end select
      endif

      write(37,'(f8.2,11d13.5)') timeh+tstep/3.6d+03, &
                               1.0d-06*sum(apdf(:,1)), 1.0d-06*sum(bpdf(:,1)), 1.0d-06*sum(cpdf(:,1)), &
                               sum( apdf(:,2) ),  sum( bpdf(:,2) ),  sum( cpdf(:,2) ), &
                               sum( apdf(:,3) ),  sum( bpdf(:,3) ),  sum( cpdf(:,3) ), &
                               sum( apdf(:,2) ) + sum( bpdf(:,2) ) + sum( cpdf(:,2) ), &
                               sum( apdf(:,3) ) + sum( bpdf(:,3) ) + sum( cpdf(:,3) )

      return
      end subroutine Discrete


      subroutine Discrete_Intracoag( pdf )
!-----------------------------------------------------------------------------------------------------------------------
!     10-30-06, dlw: discrete model of the pdf used to obtain results for evaluation of matrix. 
!-----------------------------------------------------------------------------------------------------------------------

      ! arguments.
 
      real(8), intent(inout) :: pdf(nbins,1+nspcs)    ! working array for discrete variables [#/m^3], [ug/m^3]

      ! local variables.

      integer :: i, j, nlo, nhi 
      real(8) :: w(nbins), delw(nbins), dw, delm(nbins,nspcs), mp(nbins,nspcs), dm(nspcs)

      w(:)      = pdf(:,1)   ! load number concentrations in work array [#/m^3]
      delw(:)   = 0.0d+00
      delm(:,:) = 0.0d+00
      do i=1, nbins
        mp(i,:) = pdf(i,2:nspcs+1) / max( w(i), 1.0d-30 )   ! load mass per particle for each bin and species [ug]
      enddo
      do i=1, nbins
         do j=i, nbins
            if(i.eq.j) then
               dw        = 0.5d+00*kij_discrete(i,j)*w(i)*w(j)*dtcoag
               dm(:)     = 2.0d+00 * dw * mp(i,:)
               delw(i)   = delw(i)   - 2.0d+00*dw
               delm(i,:) = delm(i,:) - dm(:)
            else
               dw = kij_discrete(i,j)*w(i)*w(j)*dtcoag
               dm(:)     = dw * ( mp(i,:) + mp(j,:) )
               delw(i)   = delw(i)   - dw
               delw(j)   = delw(j)   - dw
               delm(i,:) = delm(i,:) - dw*mp(i,:)
               delm(j,:) = delm(j,:) - dw*mp(j,:)
            endif
            ! write(*,'(2i4,7d12.3)')i,j,dw,wlo,whi,w(i),w(j),frac_lo(i,j),frac_hi(i,j)
            nlo         = nlo_grid(i,j)
            nhi         = nhi_grid(i,j)
            delw(nlo)   = delw(nlo)   + dw    * frac_lo(i,j)
            delw(nhi)   = delw(nhi)   + dw    * frac_hi(i,j)
            delm(nlo,:) = delm(nlo,:) + dm(:) * frac_lo(i,j)
            delm(nhi,:) = delm(nhi,:) + dm(:) * frac_hi(i,j)
         enddo
      enddo
      pdf(:,1)         = w(:)             + delw(:)          ! update output array [#/m^3]
      pdf(:,2:nspcs+1) = pdf(:,2:nspcs+1) + delm(:,:)        ! update output array [ug/m^3]
      do i=1, nspcs+1
        pdf(:,i) = max( pdf(:,i), 0.0d-30 )
      enddo


      return
      end subroutine Discrete_Intracoag


      subroutine Discrete_Intercoag_Ab( pdfa, pdfb )
!-----------------------------------------------------------------------------------------------------------------------
!     10-30-06, dlw: discrete model of the pdf used to obtain results for evaluation of matrix. 
!
!     pdfa coagulates with pdfb to produce additional particles in pdfb. 
!-----------------------------------------------------------------------------------------------------------------------

      ! arguments.

      real(8), intent(inout) :: pdfa(nbins,1+nspcs)    ! working array for discrete variables [#/m^3], [ug/m^3]
      real(8), intent(inout) :: pdfb(nbins,1+nspcs)    ! working array for discrete variables [#/m^3], [ug/m^3]

      ! local variables.

      integer :: i, j, nlo, nhi 
      real(8) :: dw, dm(nspcs)
      real(8) :: wa(nbins), delwa(nbins), delma(nbins,nspcs), mpa(nbins,nspcs), dma(nspcs)
      real(8) :: wb(nbins), delwb(nbins), delmb(nbins,nspcs), mpb(nbins,nspcs), dmb(nspcs)

      wa(:)      = pdfa(:,1)             ! load number concentrations in work array [#/m^3]
      wb(:)      = pdfb(:,1)             ! load number concentrations in work array [#/m^3]
      delwa(:)   = 0.0d+00
      delwb(:)   = 0.0d+00
      delma(:,:) = 0.0d+00
      delmb(:,:) = 0.0d+00
      do i=1, nbins
        mpa(i,:) = pdfa(i,2:nspcs+1) / max( wa(i), 1.0d-30 )   ! load mass per particle for each bin and species [ug]
        mpb(i,:) = pdfb(i,2:nspcs+1) / max( wb(i), 1.0d-30 )   ! load mass per particle for each bin and species [ug]
      enddo
      do i=1, nbins
         do j=1, nbins
            dw           = kij_discrete(i,j)*wa(i)*wb(j)*dtcoag
            dma(:)       = dw * mpa(i,:)
            dmb(:)       = dw * mpb(j,:)
            dm(:)        = dma(:) + dmb(:)
            delwa(i)     = delwa(i)   - dw
            delwb(j)     = delwb(j)   - dw
            delma(i,:)   = delma(i,:) - dma(:)
            delmb(j,:)   = delmb(j,:) - dmb(:)
            nlo          = nlo_grid(i,j)
            nhi          = nhi_grid(i,j)
            delwb(nlo)   = delwb(nlo)   + dw    * frac_lo(i,j)
            delwb(nhi)   = delwb(nhi)   + dw    * frac_hi(i,j)
            delmb(nlo,:) = delmb(nlo,:) + dm(:) * frac_lo(i,j)
            delmb(nhi,:) = delmb(nhi,:) + dm(:) * frac_hi(i,j)
         enddo
      enddo
      pdfa(:,1)         = wa(:)             + delwa(:)          ! update output array [#/m^3]
      pdfb(:,1)         = wb(:)             + delwb(:)          ! update output array [#/m^3]
      pdfa(:,2:nspcs+1) = pdfa(:,2:nspcs+1) + delma(:,:)        ! update output array [ug/m^3]
      pdfb(:,2:nspcs+1) = pdfb(:,2:nspcs+1) + delmb(:,:)        ! update output array [ug/m^3]

      do i=1, nspcs+1
        pdfa(:,i) = max( pdfa(:,i), 0.0d-30 )
        pdfb(:,i) = max( pdfb(:,i), 0.0d-30 )
      enddo

      return
      end subroutine Discrete_Intercoag_Ab


      subroutine Discrete_Intercoag_Abc( pdfa, pdfb, pdfc )
!-----------------------------------------------------------------------------------------------------------------------
!     10-30-06, dlw: discrete model of the pdf used to obtain results for evaluation of matrix. 
!
!     pdfa coagulates with pdfb to produce additional particles in pdfc. 
!     either pdfa or pdfb may be identical with pdfc, but pdfa and pdfb cannot be identical.
!
!     if pdfa is not pdfc, and pdfb is not pdfc --> icase = 0
!     if pdfa is. pdfc                          --> icase = 1
!     if pdfb is. pdfc                          --> icase = 2
!-----------------------------------------------------------------------------------------------------------------------

      ! arguments.

      real(8), intent(inout) :: pdfa(nbins,1+nspcs)    ! working array for discrete variables [#/m^3], [ug/m^3]
      real(8), intent(inout) :: pdfb(nbins,1+nspcs)    ! working array for discrete variables [#/m^3], [ug/m^3]
      real(8), intent(inout) :: pdfc(nbins,1+nspcs)    ! working array for discrete variables [#/m^3], [ug/m^3]

      ! local variables.

      integer :: i, j, nlo, nhi 
      real(8) :: dw, dm(nspcs)
      real(8) :: wa(nbins), delwa(nbins), delma(nbins,nspcs), mpa(nbins,nspcs), dma(nspcs)
      real(8) :: wb(nbins), delwb(nbins), delmb(nbins,nspcs), mpb(nbins,nspcs), dmb(nspcs)
      real(8) ::            delwc(nbins), delmc(nbins,nspcs)

      wa(:)      = pdfa(:,1)             ! load number concentrations in work array [#/m^3]
      wb(:)      = pdfb(:,1)             ! load number concentrations in work array [#/m^3]
      delwa(:)   = 0.0d+00
      delwb(:)   = 0.0d+00
      delwc(:)   = 0.0d+00
      delma(:,:) = 0.0d+00
      delmb(:,:) = 0.0d+00
      delmc(:,:) = 0.0d+00
      do i=1, nbins
        mpa(i,:) = pdfa(i,2:nspcs+1) / max( wa(i), 1.0d-30 )   ! load mass per particle for each bin and species [ug]
        mpb(i,:) = pdfb(i,2:nspcs+1) / max( wb(i), 1.0d-30 )   ! load mass per particle for each bin and species [ug]
        ! write(39,'(i6,8d15.5)') i, mpb(i,:), pdfb(i,2:nspcs+1), wb(i)
      enddo
      ! write(39,*)'  '
      do i=1, nbins
         do j=1, nbins
            dw           = kij_discrete(i,j)*wa(i)*wb(j)*dtcoag
!         if( dw .gt. min( wa(i), wb(j) ) ) then
!         write(*,*)'dw .gt. min( wa(i), wb(j) ): i,j, dw, wa, wb = ', i, j, dw, wa(i), wb(j)
!         stop
!       endif
            dma(:)       = dw * mpa(i,:)
            dmb(:)       = dw * mpb(j,:)
            dm(:)        = dma(:) + dmb(:)
            delwa(i)     = delwa(i)   - dw
            delwb(j)     = delwb(j)   - dw
            delma(i,:)   = delma(i,:) - dma(:)
            delmb(j,:)   = delmb(j,:) - dmb(:)
            nlo          = nlo_grid(i,j)
            nhi          = nhi_grid(i,j)
            delwc(nlo)   = delwc(nlo)   + dw    * frac_lo(i,j)
            delwc(nhi)   = delwc(nhi)   + dw    * frac_hi(i,j)
            delmc(nlo,:) = delmc(nlo,:) + dm(:) * frac_lo(i,j)
            delmc(nhi,:) = delmc(nhi,:) + dm(:) * frac_hi(i,j)
         enddo
      enddo
      pdfa(:,1)         = wa(:)             + delwa(:)          ! update output array [#/m^3]
      pdfb(:,1)         = wb(:)             + delwb(:)          ! update output array [#/m^3]
      pdfc(:,1)         = pdfc(:,1)         + delwc(:)          ! update output array [#/m^3]
      pdfa(:,2:nspcs+1) = pdfa(:,2:nspcs+1) + delma(:,:)        ! update output array [ug/m^3]
      pdfb(:,2:nspcs+1) = pdfb(:,2:nspcs+1) + delmb(:,:)        ! update output array [ug/m^3]
      pdfc(:,2:nspcs+1) = pdfc(:,2:nspcs+1) + delmc(:,:)        ! update output array [ug/m^3]

      do i=1, nspcs+1
        pdfa(:,i) = max( pdfa(:,i), 0.0d-30 )
        pdfb(:,i) = max( pdfb(:,i), 0.0d-30 )
      enddo

      return
      end subroutine Discrete_Intercoag_Abc


      subroutine Discrete_Out(iunit,icset,timeh)
!-----------------------------------------------------------------------------------------------------------------------
!     11-01-06, dlw
!-----------------------------------------------------------------------------------------------------------------------
      implicit none

      ! arguments.
 
      integer, intent(   in) :: iunit                 ! output logical unit number [1]
      integer, intent(   in) :: icset                 ! identifies test case [1]
      real(8), intent(   in) :: timeh                 ! model time [h]

      ! local variables.

      integer :: i
      real(8) :: adndlogd,  bdndlogd,  cdndlogd
      real(8) :: adm1dlogd, bdm1dlogd, cdm1dlogd
      real(8) :: adm2dlogd, bdm2dlogd, cdm2dlogd

      do i=1, nbins
        adndlogd  = apdf(i,1) * rdlogdsc * 1.0d-06           ! convert from [#/m^3] to [#/cm^3]
        bdndlogd  = bpdf(i,1) * rdlogdsc * 1.0d-06           ! convert from [#/m^3] to [#/cm^3]
        cdndlogd  = cpdf(i,1) * rdlogdsc * 1.0d-06           ! convert from [#/m^3] to [#/cm^3]
        adndlogd  = max( adndlogd, 1.0d-30 )
        bdndlogd  = max( bdndlogd, 1.0d-30 )
        cdndlogd  = max( cdndlogd, 1.0d-30 )
        adm1dlogd = apdf(i,2) * rdlogdsc                     ! [ug/m^3]
        bdm1dlogd = bpdf(i,2) * rdlogdsc                     ! [ug/m^3]
        cdm1dlogd = cpdf(i,2) * rdlogdsc                     ! [ug/m^3]
        adm1dlogd = max( adm1dlogd, 1.0d-30 )
        bdm1dlogd = max( bdm1dlogd, 1.0d-30 )
        cdm1dlogd = max( cdm1dlogd, 1.0d-30 )
        adm2dlogd = apdf(i,3) * rdlogdsc                     ! [ug/m^3]
        bdm2dlogd = bpdf(i,3) * rdlogdsc                     ! [ug/m^3]
        cdm2dlogd = cpdf(i,3) * rdlogdsc                     ! [ug/m^3]
        adm2dlogd = max( adm2dlogd, 1.0d-30 )
        bdm2dlogd = max( bdm2dlogd, 1.0d-30 )
        cdm2dlogd = max( cdm2dlogd, 1.0d-30 )
        write(iunit,91) i, dgrid(i), adndlogd,  bdndlogd,  cdndlogd,  &
                                     adm1dlogd, bdm1dlogd, cdm1dlogd, &
                                     adm2dlogd, bdm2dlogd, cdm2dlogd
         
      enddo
      
91    format(i5,10d14.6)
      return
      end subroutine Discrete_Out


      real(8) function Fln(x,xg,sigmag)
      real(8) :: x      ! particle radius or diameter [any units]
      real(8) :: xg     ! geometric mean radius or diameter [any units]
      real(8) :: sigmag ! geometric standard deviation [monodisperse = 1.0]
      real(8), parameter :: sqrttwopi = 2.506628275d+00
      fln = exp(-0.5d+00*(log(x/xg)/log(sigmag))**2) / (x*log(sigmag)*sqrttwopi)
      return
      end function Fln

 
      end module Aero_Discrete

