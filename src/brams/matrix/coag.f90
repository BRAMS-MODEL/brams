!>-------------------------------------------------------------------------------------------------------------------------------------
!!
!!@brief     This module contains all routines for coagulation coefficients.
!!@author    Susanne Bauer/Doug Wright
!!
!!     The routines are included in this file in the following order.
!!     Only routine KBARNIJ is called after initialization is complete.
!!      
!!           setup_kij_diameters
!!           setup_kij
!!           setup_kij_tables
!!           build_kij_tables
!!           get_knij
!!           get_kbarnij
!!           brownian_coag_coef
!!           cbde_coag_coef
!!           gravcoll_coag_coef
!!           turb_coag_coef
!!           total_coag_coef
!!           test_coag_coef
!!
!!-------------------------------------------------------------------------------------------------------------------------------------
   MODULE Aero_Coag
   
   use memMatrix
!LFR>       use Aero_Param,  only: &
!LFR>                            aunit1,       & !intent()
!LFR>                            aunit2,       & !intent()
!LFR>                            write_log,    & !intent()
!LFR>                            kij_ndgs_set, & !intent()
!LFR>                            dpmin_global    !intent()
!LFR>       use Aero_Config, only: &
!LFR>                            nweights        !intent()
!LFR>       use Aero_Setup,  only: &
!LFR>                            citable         !intent()
      implicit none



      CONTAINS


!>-----------------------------------------------------------------------------------------------------------------------
!!     routine to define the geometric mean diameters dg of the lookup tables.
!!-----------------------------------------------------------------------------------------------------------------------
      subroutine Setup_Kij_Diameters
      implicit none
      integer :: i
      real(4) :: e, scale
      e = 1.0 / real( kij_ndgs - 1 )
      scale = ( kij_dgmax / kij_dgmin )**e
      do i=1, kij_ndgs
        kij_diameters(i) = kij_dgmin * scale**(i-1)             ! [um]
        ! write(*,'(i6,f16.6)') i, kij_diameters(i)
      enddo
      return
      end subroutine Setup_Kij_Diameters


!>-----------------------------------------------------------------------------------------------------------------------
!!     Routine to setup the KIJ array of coagulation coefficients [m^3/s].
!!     These are not currently used.
!!
!!     The KIJ(NWEIGHTS,NWEIGHTS) are constant coagulation coefficients
!!     for each mode-mode interaction, based upon characteristic sizes
!!     for each mode. A uniform and constant temperature and pressure are used. 
!!-----------------------------------------------------------------------------------------------------------------------
      subroutine Setup_Kij
      use Aero_Setup, only: &
                           dp0  ! [m] default value of diameter of average mass - intent()
      implicit none                              !     for the assumed lognormal for each mode      
      integer :: i, j
      real(4) :: dpi, dpj                        ! [um]
      real(4) :: betaij                          ! [m^3/s]
      real(4), parameter :: set_temp = 288.15    ! [k]
      real(4), parameter :: set_pres = 101325.0  ! {pa]

      betaij = 0.0
      if( write_log ) write(aunit1,'(/5A17/)') 'I','J','DPI[um]','DPJ[um]','KIJ(I,J)[m^3/s]'
      do i=1, nweights
         do j=1, nweights
            dpi = dp0(i)*1.0e+06
            dpj = dp0(j)*1.0e+06
!       call brownian_coag_coef( dpi, dpj, set_temp, set_pres, betaij )
            call Total_Coag_Coef   ( dpi, dpj, set_temp, set_pres, betaij )
            kij(i,j) = betaij  
            kij(j,i) = kij(i,j)
            ! if( write_log ) write(aunit1,90000) i, j, dpi, dpj, kij(i,j)
         enddo
      enddo

90000 format(2i17,2f17.6,d17.5)
      return
      end subroutine Setup_Kij


!>-----------------------------------------------------------------------------------------------------------------------
!!     routine to setup tables of mode-average coagulation coefficients [m^3/s].
!!     several temperatures and pressures are currently used. 
!!-----------------------------------------------------------------------------------------------------------------------
      subroutine Setup_Kij_Tables
      implicit none
      integer :: i
      real(4) :: k0ij_table(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)        ! [m^3/s]
      real(4) :: k3ij_table(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)        ! [m^3/s]

      !-------------------------------------------------------------------------
      ! build the table for each choice of temperature and pressure. 
      !-------------------------------------------------------------------------
      call Build_Kij_Tables(kij_temp1,kij_pres1,k0ij_table,k3ij_table)  ! t1, p1
      k0ij_temp1pres1(:,:,:,:) = k0ij_table(:,:,:,:)      ! [m^3/s]
      k3ij_temp1pres1(:,:,:,:) = k3ij_table(:,:,:,:)      ! [m^3/s]
      call Build_Kij_Tables(kij_temp1,kij_pres2,k0ij_table,k3ij_table)  ! t1, p2
      k0ij_temp1pres2(:,:,:,:) = k0ij_table(:,:,:,:)      ! [m^3/s]
      k3ij_temp1pres2(:,:,:,:) = k3ij_table(:,:,:,:)      ! [m^3/s]
      call Build_Kij_Tables(kij_temp1,kij_pres3,k0ij_table,k3ij_table)  ! t1, p3
      k0ij_temp1pres3(:,:,:,:) = k0ij_table(:,:,:,:)      ! [m^3/s]
      k3ij_temp1pres3(:,:,:,:) = k3ij_table(:,:,:,:)      ! [m^3/s]
      call Build_Kij_Tables(kij_temp2,kij_pres1,k0ij_table,k3ij_table)  ! t2, p1
      k0ij_temp2pres1(:,:,:,:) = k0ij_table(:,:,:,:)      ! [m^3/s]
      k3ij_temp2pres1(:,:,:,:) = k3ij_table(:,:,:,:)      ! [m^3/s]
      call Build_Kij_Tables(kij_temp2,kij_pres2,k0ij_table,k3ij_table)  ! t2, p2
      k0ij_temp2pres2(:,:,:,:) = k0ij_table(:,:,:,:)      ! [m^3/s]
      k3ij_temp2pres2(:,:,:,:) = k3ij_table(:,:,:,:)      ! [m^3/s]
      call Build_Kij_Tables(kij_temp2,kij_pres3,k0ij_table,k3ij_table)  ! t2, p3
      k0ij_temp2pres3(:,:,:,:) = k0ij_table(:,:,:,:)      ! [m^3/s]
      k3ij_temp2pres3(:,:,:,:) = k3ij_table(:,:,:,:)      ! [m^3/s]
      call Build_Kij_Tables(kij_temp3,kij_pres1,k0ij_table,k3ij_table)  ! t3, p1
      k0ij_temp3pres1(:,:,:,:) = k0ij_table(:,:,:,:)      ! [m^3/s]
      k3ij_temp3pres1(:,:,:,:) = k3ij_table(:,:,:,:)      ! [m^3/s]
      call Build_Kij_Tables(kij_temp3,kij_pres2,k0ij_table,k3ij_table)  ! t3, p2
      k0ij_temp3pres2(:,:,:,:) = k0ij_table(:,:,:,:)      ! [m^3/s]
      k3ij_temp3pres2(:,:,:,:) = k3ij_table(:,:,:,:)      ! [m^3/s]
      call Build_Kij_Tables(kij_temp3,kij_pres3,k0ij_table,k3ij_table)  ! t3, p3
      k0ij_temp3pres3(:,:,:,:) = k0ij_table(:,:,:,:)      ! [m^3/s]
      k3ij_temp3pres3(:,:,:,:) = k3ij_table(:,:,:,:)      ! [m^3/s]

      return
      end subroutine Setup_Kij_Tables


!>-------------------------------------------------------------------------------------------------------------------------------------
!!     routine to setup a table of mode-average coagulation coefficients [m^3/s]
!!       for a given temperature and pressure.
!!-------------------------------------------------------------------------------------------------------------------------------------
      subroutine Build_Kij_Tables( temp, pres, k0ij_table, k3ij_table )
      implicit none
      integer :: i, j, k, l
      real(4) :: temp                                             !< [k]
      real(4) :: pres                                             !< [pa]
      real(4) :: k0ij_table(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)  !< [m^3/s]
      real(4) :: k3ij_table(kij_ndgs,kij_nsgs,kij_ndgs,kij_nsgs)  !< [m^3/s]
      real(4) :: k0ij                                             !< [m^3/s]
      real(4) :: k3ij                                             !< [m^3/s]

      do i=1, kij_ndgs
         do j=1, kij_ndgs
            do k=1, kij_nsgs
               do l=1, kij_nsgs
                  call Get_Knij(temp,pres,kij_diameters(i),kij_sigmas(k), &
                                kij_diameters(j),kij_sigmas(l),k0ij,k3ij)
                  k0ij_table(i,k,j,l) = k0ij    ! [m^3/s]
                  k3ij_table(i,k,j,l) = k3ij    ! [m^3/s]
               enddo
            enddo
         enddo
      enddo

      return
      end subroutine Build_Kij_Tables


!>-------------------------------------------------------------------------------------------------------------------------------------
!!     Routine to calculate the mode-average coagulation coefficients 
!!     K0IJ and K3IJ [m^3/s] for the coagulation of mode I, having lognormal
!!     parameters DGI and SIGGI, with mode J, having lognormal parameters
!!     DGJ and SIGGJ, at a given temperature and pressure.
!!
!!     The integrals necessary to evaluate the mode-average coagulation
!!     coefficients are evaluated using n-point quadrature from 2n moments
!!     calculated from the lognormal parameters for each mode.
!!-------------------------------------------------------------------------------------------------------------------------------------
      subroutine Get_Knij( temp, pres, dgi, siggi, dgj, siggj, k0ij, k3ij )
      IMPLICIT NONE
      
      ! Arguments.

      real(4),INTENT(IN) :: temp   !< [k]  ambient temperature
      real(4),INTENT(IN) :: pres   !< [pa] ambient pressure
      real(4),INTENT(IN) :: dgi    !< [um] geometric mean diameters
      real(4),INTENT(IN) :: dgj    !< [um] geometric mean diameters
      real(4),INTENT(IN) :: siggi  !< [1]  geometric standard deviations
      real(4),INTENT(IN) :: siggj  !< [1]  geometric standard deviations
      real(4),INTENT(OUT) :: k0ij   !< [m^3/s] mode-average coagulation coefficients
      real(4),INTENT(OUT) :: k3ij   !< [m^3/s] mode-average coagulation coefficients

      ! local variables.

      integer, parameter :: npoints = 4         ! number of quadrature points for each mode

      integer :: i, j, k, l                     ! loop indices
      real(8) :: sgi, sgj                       ! [1] function of geometric standard deviation
      real(8) :: uki(2*npoints)                 ! [um^k] normalized moments for mode i
      real(8) :: ukj(2*npoints)                 ! [um^k] normalized moments for mode j
      real(8) :: xi(npoints), xj(npoints)       ! [um] quadrature abscissas for modes i and j
      real(8) :: wi(npoints), wj(npoints)       ! [1] normalized quadrature weights for modes i and j
      real(8) :: zfi, zfj                       ! [1] flags for failed quadrature inversion
      real(8) :: k0ij_tmp, k3ij_tmp             ! [m^3/s] double precision accumulators for k0ij and k3ij
      real(4) :: di, dj                         ! [um] single precision particle diameters for modes i and j
      real(4) :: betaij                         ! [m^3/s] single precision coagulation coefficient


      !write(89,fmt='("Inputs:",6(D15.5,1X))') temp,pres,dgi,siggi,dgj,siggj
      !write(89,fmt='(A)') '.....................................................'

      sgi = exp( 0.5d+00 * ( log( dble(siggi) )**2 ) )
      sgj = exp( 0.5d+00 * ( log( dble(siggj) )**2 ) )
      !-------------------------------------------------------------------------
      ! write(*,*)         'dgi,siggi,sgi,dgj,siggj,sgj'
      ! write(*,'(6f13.8)') dgi,siggi,sgi,dgj,siggj,sgj 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      ! compute the first 2*npoints diameter moments for each lognormal mode.
      !-------------------------------------------------------------------------
      do l=1, 2*npoints
        k = l-1
        uki(l) = dble(dgi)**k * sgi**(k*k)
        ukj(l) = dble(dgj)**k * sgj**(k*k)
        !-----------------------------------------------------------------------
        ! write(*,'(i6,2d15.5)') k, uki(l), ukj(l)
        !-----------------------------------------------------------------------
      enddo

      !-------------------------------------------------------------------------
      ! get the quadrature abscissas and weights for modes i and j.
      !-------------------------------------------------------------------------
      call Gauss_matrix(npoints,uki,xi,wi,zfi)
      if( zfi .gt. 1.0d-15 ) then
        WRITE(*,*)'Failed quadrature for mode I in subr. GET_KNIJ'
        WRITE(*,*)'ABSCISSAS = ', xi(:)
        WRITE(*,*)'WEIGHTS   = ', wi(:)
        STOP
      ENDIF
      
      !write(89,fmt='("npoints=",I3.3)') npoints
      !write(89,fmt='("uki=",8(D15.5,1X))') uki
      !write(89,fmt='(" xi=",4(D15.5,1X))') xi
      !write(89,fmt='(" wi=",4(D15.5,1X))') wi
      !write(89,fmt='("zfi=",4(D15.5,1X))') zfi
      
      call Gauss_Matrix(npoints,ukj,xj,wj,zfj)
      if( zfj .gt. 1.0d-15 ) then
        WRITE(*,*)'Failed quadrature for mode J in subr. GET_KNIJ'
        WRITE(*,*)'ABSCISSAS = ', xj(:)
        write(*,*)'weights   = ', wj(:)
        STOP
      ENDIF
      !write(89,fmt='("ukj=",8(D15.5,1X))') ukj
      !write(89,fmt='(" xj=",4(D15.5,1X))') xj
      !write(89,fmt='(" wj=",4(D15.5,1X))') wj
      !write(89,fmt='("zfj=",4(D15.5,1X))') zfj
      !-------------------------------------------------------------------------
      ! Write the abscissas and weights for modes I and J.
      !-------------------------------------------------------------------------
      ! DO L=1, NPOINTS
      !   WRITE(*,'(I6,4D15.5)') L, XI(L), WI(L), XJ(L), WJ(L)
      ! ENDDO
      !-------------------------------------------------------------------------

      k0ij_tmp = 0.0d+00
      k3ij_tmp = 0.0d+00

      do i=1, npoints
         do j=1, npoints
            di = real(xi(i))                                           ! convert to single precision
            dj = real(xj(j))                                           ! convert to single precision
!               call brownian_coag_coef( di, dj, temp, pres, betaij )      ! all variables single precision
            call Total_Coag_Coef   ( di, dj, temp, pres, betaij )      ! all variables single precision
            k0ij_tmp = k0ij_tmp + dble(betaij)*wi(i)*wj(j)             ! all factors are real(8)
            k3ij_tmp = k3ij_tmp + dble(betaij)*wi(i)*wj(j)*xi(i)**3    ! all factors are real(8)
         enddo
      enddo

      k0ij = real(k0ij_tmp)                                        ! divide by u0i*u0j=   1.0*1.0
      k3ij = real(k3ij_tmp/uki(4))                                 ! divide by u3i*u0j=uki(4)*1.0

      !write(89,fmt='("k0ij=",D15.5)') k0ij
      !write(89,fmt='("k3ij=",D15.5)') k3ij
      !write(89,fmt='(A)') '======================================================'
      !call flush(6)

!-------------------------------------------------------------------------------
!     write(aunit2,'(a,4f9.4,2e13.5)') 'dgi, siggi, dgj, siggj, k0ij, k3ij = ',
!    &                                  dgi, siggi, dgj, siggj, k0ij, k3ij
!-------------------------------------------------------------------------------
            
      return
      end subroutine Get_Knij


!>-------------------------------------------------------------------------------------------------------------------------------------
!!     Routine to setup tables of mode-average coagulation coefficients
!!     KBAR0IJ and KBAR3IJ [m^3/s] for arbitrary temperature and pressure.
!!     Modes are assumed to be lognormal, and Sigmag values are set to
!!     constants for each mode, and Dg values are derived from the current
!!     value of the diameter of average mass for each mode, before being passed
!!     to this routine.
!!-------------------------------------------------------------------------------------------------------------------------------------
      subroutine Get_kBarnIj( iupdate, tk, pres, diam, kbar0ij, kbar3ij )
      use Aero_Setup, only: &
                        sig0  ! [um], [1], default lognormal parameters - intent()
      implicit none                                     !            for each mode

      ! arguments.

      integer :: iupdate                                !< [1]  control flag
      real(8) :: tk                                     !< [k]  ambient temperature
      real(8) :: pres                                   !< [pa] ambient pressure
      real(8) :: diam(nweights)                         !< [um] geo. mean diameter for each mode
      real(8) :: kbar0ij(nweights,nweights)             !< [m^3/s] 0th mode-average coag. coef.
      real(8) :: kbar3ij(nweights,nweights)             !< [m^3/s] 3th mode-average coag. coef.

      ! local variables.

      integer       :: index_diami, index_diamj, index_diamip1, index_diamjp1 
      integer       :: i, j, itrange, iprange                     
      !------------------------------------------------------------------------------------------------------------
      ! index_sigg(i) is the index of kij_sigmas to obtain sigmag for mode i.
      !------------------------------------------------------------------------------------------------------------
      !integer, save :: index_sigg(nweights)                               ! [1]
      real(4), save :: deltalndg                                          ! [1] table spacing in ln(dg)
      real(4)       :: kbar0ij_ll, kbar0ij_lu, kbar0ij_ul, kbar0ij_uu     ! [m^3/s] for bilinear interpolation
      real(4)       :: kbar3ij_ll, kbar3ij_lu, kbar3ij_ul, kbar3ij_uu     ! [m^3/s] for bilinear interpolation
      real(4)       :: xinterpi, xinterpj, xinterpt, xinterpp             ! [1] for bilinear interpolation
      real(4)       :: tmp0, tmp3                                         ! [m^3/s] scratch variables
      real(4)       :: tpinterp_ll, tpinterp_lu, tpinterp_ul, tpinterp_uu ! [1] scratch variables
      real(4)       :: tuse                                               ! [k]  ambient temperature local variable 
      real(4)       :: puse                                               ! [pa] ambient pressure    local variable 
      real(4)       :: kbar0ij_ll_ll, kbar0ij_lu_ll, kbar0ij_ul_ll, kbar0ij_uu_ll  ! [m^3/s] for bilinear interp.
      real(4)       :: kbar3ij_ll_ll, kbar3ij_lu_ll, kbar3ij_ul_ll, kbar3ij_uu_ll  ! [m^3/s] for bilinear interp.
      real(4)       :: kbar0ij_ll_lu, kbar0ij_lu_lu, kbar0ij_ul_lu, kbar0ij_uu_lu  ! [m^3/s] for bilinear interp.
      real(4)       :: kbar3ij_ll_lu, kbar3ij_lu_lu, kbar3ij_ul_lu, kbar3ij_uu_lu  ! [m^3/s] for bilinear interp.
      real(4)       :: kbar0ij_ll_ul, kbar0ij_lu_ul, kbar0ij_ul_ul, kbar0ij_uu_ul  ! [m^3/s] for bilinear interp.
      real(4)       :: kbar3ij_ll_ul, kbar3ij_lu_ul, kbar3ij_ul_ul, kbar3ij_uu_ul  ! [m^3/s] for bilinear interp.
      real(4)       :: kbar0ij_ll_uu, kbar0ij_lu_uu, kbar0ij_ul_uu, kbar0ij_uu_uu  ! [m^3/s] for bilinear interp.
      real(4)       :: kbar3ij_ll_uu, kbar3ij_lu_uu, kbar3ij_ul_uu, kbar3ij_uu_uu  ! [m^3/s] for bilinear interp.
      logical, save :: firstime = .true.
      logical       :: flag

      integer :: arrindex_diam(nweights)
      real(4) :: arrxinterp(nweights)

      if( firstime ) then
        firstime = .false.
        !------------------------------------------------------------------------------------------------------------
        ! deltalng is the table spacing in ln(dg) and needed for interpolation.
        !------------------------------------------------------------------------------------------------------------
        deltalndg = log( kij_dgmax / kij_dgmin ) / real( kij_ndgs - 1 )  ! to interpolate in dg
        !------------------------------------------------------------------------------------------------------------
        ! to efficiently identify the sigmag value assigned to each mode.
        !------------------------------------------------------------------------------------------------------------
        index_sigg(:) = 1
        do i=1, nweights 
          do j=1, kij_nsgs
            if( abs( sig0(i) - dble( kij_sigmas(j) ) ) .lt. 1.0d-03 ) index_sigg(i) = j
          enddo
          ! write(aunit2,'(i6,f12.6,i6,f12.6)') i, sig0(i), index_sigg(i), diam(i)
        enddo
      endif

      !--------------------------------------------------------------------------------------------------------------
      ! temporary code to check the new value of kij_dgmin. 
      !--------------------------------------------------------------------------------------------------------------
      ! if( minval( diam(:) ) .lt. kij_dgmin ) then
      !  write(*,*)'minval( diam(:) ) .lt. kij_dgmin  in subr. get_kbarnij. :', minval( diam(:) ) 
      !  stop
      ! endif
      !--------------------------------------------------------------------------------------------------------------

      ! precompute common elements
      do i=1, nweights
        !----------------------------------------------------------------------------------------------------------
        ! for mode i, get the lower and upper bounding table diameters and the interpolation variable xinterpi.
        !----------------------------------------------------------------------------------------------------------
        index_diami = ( log( diam(i) / kij_dgmin ) / deltalndg ) + 1
        index_diami = min( max( index_diami, 1 ), kij_ndgs-1 )
        arrindex_diam(i) = index_diami

        xinterpi = log( diam(i) / kij_diameters(index_diami) ) / deltalndg 
        arrxinterp(i) = xinterpi
      end do
      !--------------------------------------------------------------------------------------------------------------
      ! iupdate .eq. 0: the mode-average coagulation coefficients are held constant throughout the simulation.
      !--------------------------------------------------------------------------------------------------------------
      if( iupdate .eq. 0 ) then
        flag = .false.
        do i=1, nweights
          !----------------------------------------------------------------------------------------------------------
          ! for mode i, get the lower and upper bounding table diameters and the interpolation variable xinterpi.
          !----------------------------------------------------------------------------------------------------------
          index_diami = arrindex_diam(i)
          index_diamip1 = index_diami+1
          if(diam(i) .lt. kij_diameters(index_diami  )) flag = .true.
          if(diam(i) .gt. kij_diameters(index_diamip1)) flag = .true.  

          if( flag ) then
            WRITE(*,*)'Problem in GET_KBARNIJ for IUPDATE = 0'
            write(*,'(2i6,8f11.5)')i,kij_diameters(index_diami),diam(i),kij_diameters(index_diamip1)
            stop
          endif
!--------------------------------------------------------------------------------------------------------------------
!           the lower and upper table diameters bounding diam(i)), and the interpolation variables
!           xinterpi were checked and found correct.
!--------------------------------------------------------------------------------------------------------------------
!           write(aunit2,'(2i6,8f11.5)')i,kij_diameters(index_diami),diam(i),kij_diameters(index_diamip1),xinterpi
!--------------------------------------------------------------------------------------------------------------------
        end do

        do i=1, nweights
          !----------------------------------------------------------------------------------------------------------
          ! for mode i, get the lower and upper bounding table diameters and the interpolation variable xinterpi.
          !----------------------------------------------------------------------------------------------------------
          index_diami = arrindex_diam(i)
          index_diamip1 = index_diami+1
          xinterpi = arrxinterp(i)

          do j=1, nweights
            !--------------------------------------------------------------------------------------------------------
            ! for mode j, get the lower and upper bounding table diameters and the interpolation variable xinterpj.
            !--------------------------------------------------------------------------------------------------------
            index_diamj = arrindex_diam(j)
            index_diamjp1 = index_diamj+1
            xinterpj = arrxinterp(j)
            !--------------------------------------------------------------------------------------------------------
            ! for each of the four points needed for bilinear interpolation, get the mode-average coagulation
            ! coefficients at the selected temperature and pressure.
            !--------------------------------------------------------------------------------------------------------
            kbar0ij_ll = k0ij_temp2pres1( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
            kbar0ij_lu = k0ij_temp2pres1( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
            kbar0ij_ul = k0ij_temp2pres1( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
            kbar0ij_uu = k0ij_temp2pres1( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )
            kbar3ij_ll = k3ij_temp2pres1( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
            kbar3ij_lu = k3ij_temp2pres1( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
            kbar3ij_ul = k3ij_temp2pres1( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
            kbar3ij_uu = k3ij_temp2pres1( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )
            !--------------------------------------------------------------------------------------------------------
            ! interpolate in dg(i) and dg(j) for modes i and j.
            !
            ! when diam(i) = kij_diameters(index_diami), the lower i-mode dg value, xinterpi = 0.0, so
            ! kbar0ij_ll and kbar0ij_lu should be multiplied by (1.0 - xinterpi ) = 1.0.
            !--------------------------------------------------------------------------------------------------------
            tmp0         = kbar0ij_ll*(1.0-xinterpi)*(1.0-xinterpj) &
                         + kbar0ij_lu*(1.0-xinterpi)*(    xinterpj) &
                         + kbar0ij_ul*(    xinterpi)*(1.0-xinterpj) &
                         + kbar0ij_uu*(    xinterpi)*(    xinterpj)
            tmp3         = kbar3ij_ll*(1.0-xinterpi)*(1.0-xinterpj) &
                         + kbar3ij_lu*(1.0-xinterpi)*(    xinterpj) &
                         + kbar3ij_ul*(    xinterpi)*(1.0-xinterpj) &
                         + kbar3ij_uu*(    xinterpi)*(    xinterpj)
            kbar0ij(i,j) = dble( tmp0 )
            kbar3ij(i,j) = dble( tmp3 )
            !--------------------------------------------------------------------------------------------------------
            ! for narrow distributions, kbar0ij and kbar3ij should be nearly equal, and were found to be so
            ! with sigmag = 1.1 for all modes.
            !--------------------------------------------------------------------------------------------------------
            ! write(aunit2,'(2i6,2e15.5)')i,j,kbar0ij(i,j),kbar3ij(i,j)
            !--------------------------------------------------------------------------------------------------------
          enddo
        enddo
      !--------------------------------------------------------------------------------------------------------------
      ! iupdate .eq. 1: the mode-average coagulation coefficients are updated at each time step.
      !--------------------------------------------------------------------------------------------------------------
      elseif( iupdate .eq. 1 ) then
        tuse = min( max( real(tk),   kij_temp3 ), kij_temp1 )   ! tmin=kij_temp3, tmax=kij_temp1  
        puse = min( max( real(pres), kij_pres3 ), kij_pres1 )   ! pmin=kij_pres3, pmax=kij_pres1    
        if( tuse .gt. kij_temp2 ) then                          ! tmiddle value=kij_temp2
          itrange = 12                                          ! use temperatures 1 and 2
          xinterpt = ( tuse - kij_temp2 ) / ( kij_temp1 - kij_temp2 ) 
        else
          itrange = 23                                          ! use temperatures 2 and 3
          xinterpt = ( tuse - kij_temp3 ) / ( kij_temp2 - kij_temp3 ) 
        endif
        if( puse .gt. kij_pres2 ) then                          ! pmiddle value=kij_pres2
          iprange = 12                                          ! use pressures 1 and 2
          xinterpp = ( puse - kij_pres2 ) / ( kij_pres1 - kij_pres2 ) 
        else
          iprange = 23                                          ! use pressures 2 and 3
          xinterpp = ( puse - kij_pres3 ) / ( kij_pres2 - kij_pres3 ) 
        endif
        tpinterp_ll = (1.0-xinterpt)*(1.0-xinterpp)             ! all weight at t_lower, p_lower
        tpinterp_lu = (1.0-xinterpt)*(    xinterpp)             ! all weight at t_lower, p_upper
        tpinterp_ul = (    xinterpt)*(1.0-xinterpp)             ! all weight at t_upper, p_lower
        tpinterp_uu = 1.0-tpinterp_ll-tpinterp_lu-tpinterp_ul   ! all weight at t_upper, p_upper
!--------------------------------------------------------------------------------------------------------------------
!       write(aunit2,'(/a/)')'new step'
!       write(aunit2,'(a40,2f15.6    )')'tuse, puse = ', tuse, puse
!       write(aunit2,'(a40,2i4,2f13.7)')'itrange, iprange, xinterpt, xinterpp = ',
!    &                                   itrange, iprange, xinterpt, xinterpp
!       write(aunit2,'(a40,4f12.6    )')'tpinterp_ll, tpinterp_lu, tpinterp_ul, tpinterp_uu = ', 
!    &                                   tpinterp_ll, tpinterp_lu, tpinterp_ul, tpinterp_uu
!--------------------------------------------------------------------------------------------------------------------
        do i=1, nweights
          !----------------------------------------------------------------------------------------------------------
          ! for mode i, get the lower and upper bounding table diameters and the interpolation variable xinterpi.
          !----------------------------------------------------------------------------------------------------------
          index_diami = arrindex_diam(i)
          index_diamip1 = index_diami+1
          xinterpi = arrxinterp(i)
          do j=1, nweights

            if (citable(i,j) == 'off') then ! turn off coagulation between selected modes.
              kbar0ij(i,j) = 1.0d-30
              kbar3ij(i,j) = 1.0d-30
              cycle
            end if

            !--------------------------------------------------------------------------------------------------------
            ! for mode j, get the lower and upper bounding table diameters and the interpolation variable xinterpj.
            !--------------------------------------------------------------------------------------------------------
            index_diamj = arrindex_diam(j)
            index_diamjp1 = index_diamj+1
            xinterpj = arrxinterp(j)
            !--------------------------------------------------------------------------------------------------------
            ! for each of the four points needed for bilinear interpolation in dg(i) and dg(j), get the 
            ! mode-average coagulation coefficients at each of the four temperature-pressure points.
            !
            ! in kbar0ij_ab_cd, a indicates the upper or lower value of dg(i)
            !                   b indicates the upper or lower value of dg(j)
            !                   c indicates the upper or lower value of temperature: t1 > t2 > t3
            !                   d indicates the upper or lower value of pressure:    p1 > p2 > p3
            !--------------------------------------------------------------------------------------------------------
            if( itrange .eq. 12 .and. iprange .eq. 12 ) then

              kbar0ij_ll_ll = k0ij_temp2pres2( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_lu_ll = k0ij_temp2pres2( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar0ij_ul_ll = k0ij_temp2pres2( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_uu_ll = k0ij_temp2pres2( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ll_ll = k3ij_temp2pres2( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_lu_ll = k3ij_temp2pres2( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ul_ll = k3ij_temp2pres2( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_uu_ll = k3ij_temp2pres2( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )

              kbar0ij_ll_lu = k0ij_temp2pres1( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_lu_lu = k0ij_temp2pres1( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar0ij_ul_lu = k0ij_temp2pres1( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_uu_lu = k0ij_temp2pres1( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ll_lu = k3ij_temp2pres1( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_lu_lu = k3ij_temp2pres1( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ul_lu = k3ij_temp2pres1( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_uu_lu = k3ij_temp2pres1( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )

              kbar0ij_ll_ul = k0ij_temp1pres2( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_lu_ul = k0ij_temp1pres2( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar0ij_ul_ul = k0ij_temp1pres2( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_uu_ul = k0ij_temp1pres2( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ll_ul = k3ij_temp1pres2( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_lu_ul = k3ij_temp1pres2( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ul_ul = k3ij_temp1pres2( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_uu_ul = k3ij_temp1pres2( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )

              kbar0ij_ll_uu = k0ij_temp1pres1( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_lu_uu = k0ij_temp1pres1( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar0ij_ul_uu = k0ij_temp1pres1( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_uu_uu = k0ij_temp1pres1( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ll_uu = k3ij_temp1pres1( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_lu_uu = k3ij_temp1pres1( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ul_uu = k3ij_temp1pres1( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_uu_uu = k3ij_temp1pres1( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )

            elseif( itrange .eq. 12 .and. iprange .eq. 23 ) then

              kbar0ij_ll_ll = k0ij_temp2pres3( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_lu_ll = k0ij_temp2pres3( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar0ij_ul_ll = k0ij_temp2pres3( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_uu_ll = k0ij_temp2pres3( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ll_ll = k3ij_temp2pres3( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_lu_ll = k3ij_temp2pres3( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ul_ll = k3ij_temp2pres3( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_uu_ll = k3ij_temp2pres3( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )

              kbar0ij_ll_lu = k0ij_temp2pres2( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_lu_lu = k0ij_temp2pres2( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar0ij_ul_lu = k0ij_temp2pres2( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_uu_lu = k0ij_temp2pres2( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ll_lu = k3ij_temp2pres2( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_lu_lu = k3ij_temp2pres2( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ul_lu = k3ij_temp2pres2( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_uu_lu = k3ij_temp2pres2( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )

              kbar0ij_ll_ul = k0ij_temp1pres3( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_lu_ul = k0ij_temp1pres3( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar0ij_ul_ul = k0ij_temp1pres3( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_uu_ul = k0ij_temp1pres3( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ll_ul = k3ij_temp1pres3( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_lu_ul = k3ij_temp1pres3( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ul_ul = k3ij_temp1pres3( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_uu_ul = k3ij_temp1pres3( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )

              kbar0ij_ll_uu = k0ij_temp1pres2( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_lu_uu = k0ij_temp1pres2( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar0ij_ul_uu = k0ij_temp1pres2( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_uu_uu = k0ij_temp1pres2( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ll_uu = k3ij_temp1pres2( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_lu_uu = k3ij_temp1pres2( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ul_uu = k3ij_temp1pres2( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_uu_uu = k3ij_temp1pres2( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )

            elseif( itrange .eq. 23 .and. iprange .eq. 12 ) then

              kbar0ij_ll_ll = k0ij_temp3pres2( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_lu_ll = k0ij_temp3pres2( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar0ij_ul_ll = k0ij_temp3pres2( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_uu_ll = k0ij_temp3pres2( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ll_ll = k3ij_temp3pres2( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_lu_ll = k3ij_temp3pres2( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ul_ll = k3ij_temp3pres2( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_uu_ll = k3ij_temp3pres2( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )

              kbar0ij_ll_lu = k0ij_temp3pres1( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_lu_lu = k0ij_temp3pres1( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar0ij_ul_lu = k0ij_temp3pres1( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_uu_lu = k0ij_temp3pres1( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ll_lu = k3ij_temp3pres1( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_lu_lu = k3ij_temp3pres1( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ul_lu = k3ij_temp3pres1( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_uu_lu = k3ij_temp3pres1( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )

              kbar0ij_ll_ul = k0ij_temp2pres2( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_lu_ul = k0ij_temp2pres2( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar0ij_ul_ul = k0ij_temp2pres2( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_uu_ul = k0ij_temp2pres2( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ll_ul = k3ij_temp2pres2( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_lu_ul = k3ij_temp2pres2( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ul_ul = k3ij_temp2pres2( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_uu_ul = k3ij_temp2pres2( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )

              kbar0ij_ll_uu = k0ij_temp2pres1( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_lu_uu = k0ij_temp2pres1( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar0ij_ul_uu = k0ij_temp2pres1( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_uu_uu = k0ij_temp2pres1( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ll_uu = k3ij_temp2pres1( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_lu_uu = k3ij_temp2pres1( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ul_uu = k3ij_temp2pres1( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_uu_uu = k3ij_temp2pres1( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )

            elseif( itrange .eq. 23 .and. iprange .eq. 23 ) then

              kbar0ij_ll_ll = k0ij_temp3pres3( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_lu_ll = k0ij_temp3pres3( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar0ij_ul_ll = k0ij_temp3pres3( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_uu_ll = k0ij_temp3pres3( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ll_ll = k3ij_temp3pres3( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_lu_ll = k3ij_temp3pres3( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ul_ll = k3ij_temp3pres3( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_uu_ll = k3ij_temp3pres3( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )

              kbar0ij_ll_lu = k0ij_temp3pres2( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_lu_lu = k0ij_temp3pres2( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar0ij_ul_lu = k0ij_temp3pres2( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_uu_lu = k0ij_temp3pres2( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ll_lu = k3ij_temp3pres2( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_lu_lu = k3ij_temp3pres2( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ul_lu = k3ij_temp3pres2( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_uu_lu = k3ij_temp3pres2( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )

              kbar0ij_ll_ul = k0ij_temp2pres3( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_lu_ul = k0ij_temp2pres3( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar0ij_ul_ul = k0ij_temp2pres3( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_uu_ul = k0ij_temp2pres3( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ll_ul = k3ij_temp2pres3( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_lu_ul = k3ij_temp2pres3( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ul_ul = k3ij_temp2pres3( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_uu_ul = k3ij_temp2pres3( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )

              kbar0ij_ll_uu = k0ij_temp2pres2( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_lu_uu = k0ij_temp2pres2( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar0ij_ul_uu = k0ij_temp2pres2( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar0ij_uu_uu = k0ij_temp2pres2( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ll_uu = k3ij_temp2pres2( index_diami,  index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_lu_uu = k3ij_temp2pres2( index_diami,  index_sigg(i),index_diamjp1,index_sigg(j) )
              kbar3ij_ul_uu = k3ij_temp2pres2( index_diamip1,index_sigg(i),index_diamj,  index_sigg(j) )
              kbar3ij_uu_uu = k3ij_temp2pres2( index_diamip1,index_sigg(i),index_diamjp1,index_sigg(j) )

            else

              WRITE(*,*)'Error in GET_KBARNIJ: ITRANGE, IPRANGE = ', itrange, iprange
              stop

            endif
            !--------------------------------------------------------------------------------------------------------
            ! interpolate in t and p for each of the four points needed for the dg(i) and dg(j) interpolation.
            !--------------------------------------------------------------------------------------------------------
            kbar0ij_ll = tpinterp_ll*kbar0ij_ll_ll + tpinterp_lu*kbar0ij_ll_lu &
                       + tpinterp_ul*kbar0ij_ll_ul + tpinterp_uu*kbar0ij_ll_uu
            kbar0ij_lu = tpinterp_ll*kbar0ij_lu_ll + tpinterp_lu*kbar0ij_lu_lu &
                       + tpinterp_ul*kbar0ij_lu_ul + tpinterp_uu*kbar0ij_lu_uu
            kbar0ij_ul = tpinterp_ll*kbar0ij_ul_ll + tpinterp_lu*kbar0ij_ul_lu &
                       + tpinterp_ul*kbar0ij_ul_ul + tpinterp_uu*kbar0ij_ul_uu
            kbar0ij_uu = tpinterp_ll*kbar0ij_uu_ll + tpinterp_lu*kbar0ij_uu_lu &
                       + tpinterp_ul*kbar0ij_uu_ul + tpinterp_uu*kbar0ij_uu_uu

            kbar3ij_ll = tpinterp_ll*kbar3ij_ll_ll + tpinterp_lu*kbar3ij_ll_lu &
                       + tpinterp_ul*kbar3ij_ll_ul + tpinterp_uu*kbar3ij_ll_uu
            kbar3ij_lu = tpinterp_ll*kbar3ij_lu_ll + tpinterp_lu*kbar3ij_lu_lu &
                       + tpinterp_ul*kbar3ij_lu_ul + tpinterp_uu*kbar3ij_lu_uu
            kbar3ij_ul = tpinterp_ll*kbar3ij_ul_ll + tpinterp_lu*kbar3ij_ul_lu &
                       + tpinterp_ul*kbar3ij_ul_ul + tpinterp_uu*kbar3ij_ul_uu
            kbar3ij_uu = tpinterp_ll*kbar3ij_uu_ll + tpinterp_lu*kbar3ij_uu_lu &
                       + tpinterp_ul*kbar3ij_uu_ul + tpinterp_uu*kbar3ij_uu_uu
!--------------------------------------------------------------------------------------------------------------------
!           write(aunit2,'(a40,4e13.5)')'kbar0ij_ll, kbar0ij_lu, kbar0ij_ul, kbar0ij_uu = ',
!    &                                   kbar0ij_ll, kbar0ij_lu, kbar0ij_ul, kbar0ij_uu
!           write(aunit2,'(a40,4e13.5)')'kbar3ij_ll, kbar3ij_lu, kbar3ij_ul, kbar3ij_uu = ',
!    &                                   kbar3ij_ll, kbar3ij_lu, kbar3ij_ul, kbar3ij_uu
!--------------------------------------------------------------------------------------------------------------------
            ! interpolate in dg(i) and dg(j) for modes i and j.
            !--------------------------------------------------------------------------------------------------------
            tmp0         = kbar0ij_ll*(1.0-xinterpi)*(1.0-xinterpj) &
                         + kbar0ij_lu*(1.0-xinterpi)*(    xinterpj) &
                         + kbar0ij_ul*(    xinterpi)*(1.0-xinterpj) &
                         + kbar0ij_uu*(    xinterpi)*(    xinterpj)
            tmp3         = kbar3ij_ll*(1.0-xinterpi)*(1.0-xinterpj) &
                         + kbar3ij_lu*(1.0-xinterpi)*(    xinterpj) &
                         + kbar3ij_ul*(    xinterpi)*(1.0-xinterpj) &
                         + kbar3ij_uu*(    xinterpi)*(    xinterpj)
            kbar0ij(i,j) = dble( tmp0 )
            kbar3ij(i,j) = dble( tmp3 )
            !--------------------------------------------------------------------------------------------------------
            ! for narrow distributions, kbar0ij and kbar3ij should be nearly equal, and were found to be so
            ! with sigmag = 1.1 for all modes.
            !--------------------------------------------------------------------------------------------------------
            ! write(aunit2,'(a40,2i6,4e15.5)')'i,j,kbar0ij,kbar3ij=',i,j,kbar0ij(i,j),kbar3ij(i,j)
            !--------------------------------------------------------------------------------------------------------
          enddo
        enddo
      endif

      return
      end subroutine Get_kBarnIj


!>-------------------------------------------------------------------------------------------------------------------------------------
!!     routine to calculate the brownian coagulation coefficient betaij 
!!     for particles of diameters di and dj at ambient temperature tempk
!!     and pressure pres, using the interpolation formula of fuchs (1964),
!!     as given in jacobson 1999, p.446, eq.(16.28) 
!!              or jacobson 2005, p.509, eq.(15.33). 
!!-------------------------------------------------------------------------------------------------------------------------------------
      subroutine Brownian_Coag_Coef( di, dj, tempk, pres, betaij )
      implicit none

      real(4) :: di          !< particle diameters [um]
      real(4) :: dj          !< particle diameters [um]
      real(4) :: tempk       !< ambient temperature [k]
      real(4) :: pres        !< ambient pressure [pa]
      real(4) :: betaij      !< coagulation coefficient [m^3/s/particle]

      real(4), parameter :: pi = 3.141592653589793
      real(4), parameter :: fourpi = 4.0 * pi
      real(4), parameter :: sixpi  = 6.0 * pi
      real(4), parameter :: fourpiby3 = 4.0 * pi / 3.0
      real(4), parameter :: kb = 1.38065e-23     ! boltzmann's constant [j/k]
      real(4), parameter :: na = 6.0221367e+23   ! avogadro's number [#/mole]
      real(4), parameter :: mwair = 28.9628e-03  ! molecular weight of air in [kg/mole]
      real(4), parameter :: m1air = mwair/na     ! average mass of one air molecule [kg]
      real(4), parameter :: vaconst = 8.0*kb/pi/m1air  ! used in eq.(16.21) for vabar
      real(4), parameter :: aprime = 1.249       ! for cunningham slip-flow correction
      real(4), parameter :: bprime = 0.42        !   as given in eq.(16.25);
      real(4), parameter :: cprime = 0.87        !   kasten (1968) values used here
      real(4), parameter :: vpconst = 8.0*kb/pi  ! used in eq.(16.27) for vpibar [j/k]
      real(4), parameter :: densp_kgm3 = 1.4e+03 ! particle density [kg/m^3]

      real(4) :: ri, rj          ! particle radii [m]
      real(4) :: dpi, dpj        ! particle diffusion coefficients [m^2/s]
      real(4) :: gi, gj          ! cunningham slip-flow corrections for particles i, j [1]
      real(4) :: etaa            ! dynamic viscosity of air [kg/m/s]
      real(4) :: vabar           ! mean thermal velocity of an air molecule [m/s]
      real(4) :: rhoa            ! density of air [kg/m^3]
      real(4) :: lambdaa         ! mean free path of air [m]
      real(4) :: knai, knaj      ! knudsen numbers in air for particles i, j [1]
      real(4) :: vpibar, vpjbar  ! mean thermal velocity  for particles i, j [m/s]
      real(4) :: mpi, mpj        ! particle masses [kg]                       
      real(4) :: voli, volj      ! particle volumes [m^3]
      real(4) :: lambdapi        ! mean free path for particle i [m]
      real(4) :: lambdapj        ! mean free path for particle j [m]
      real(4) :: deltai, deltaj  ! see eq.(16.29), p. 446 of mzj 1999
      real(4) :: denom1, denom2  ! scratch variables
!     real(4) :: dum1, dum2      ! scratch variables

      ri = 0.5e-06 * di           ! convert from [um] to [m]
      rj = 0.5e-06 * dj           ! convert from [um] to [m]
      voli = fourpiby3 * ri**3    ! volume of particle i [m^3]
      volj = fourpiby3 * rj**3    ! volume of particle j [m^3]
      !-------------------------------------------------------------------------
      ! dynamic viscosity of air [kg/m/s], jacobson, 2005, eq.(4.54).
      ! mean thermal velocity of an air molecule [m/s].
      ! density of air [kg/m^3] from the ideal gas law.
      ! mean free path of an air molecule [m], jacobson, 2005, eq.(15.24).
      !-------------------------------------------------------------------------
      etaa = 1.8325e-05 * (416.16/(tempk + 120.0)) * (tempk/296.16)**1.5
      vabar = sqrt( vaconst * tempk ) 
      rhoa = pres * m1air / ( kb * tempk ) 
      lambdaa = 2.0 * etaa / ( rhoa * vabar ) 
      knai = lambdaa / ri
      knaj = lambdaa / rj
!-------------------------------------------------------------------------------
!     dum1 = 0.0
!     dum2 = 0.0
!     write(*,     '(12d11.3)') dum1, dum2
!-------------------------------------------------------------------------------
!     dum1 = knai * ( aprime + bprime * exp(-cprime/knai) )
!     dum2 = knaj * ( aprime + bprime * exp(-cprime/knaj) )
!-------------------------------------------------------------------------------
!     write(*,     '(12d11.3)') dum1, dum2
!-------------------------------------------------------------------------------
!     gi = 1.0 + dum1
!     gj = 1.0 + dum2
!-------------------------------------------------------------------------------
!     write(*,     '(12d11.3)') di,dj,knai,knaj,gi,gj
!-------------------------------------------------------------------------------
      gi = 1.0 + knai * ( aprime + bprime * exp(-cprime/knai) )
      gj = 1.0 + knaj * ( aprime + bprime * exp(-cprime/knaj) )
      dpi = kb * tempk * gi / ( sixpi * ri * etaa ) 
      dpj = kb * tempk * gj / ( sixpi * rj * etaa ) 
      mpi = voli * densp_kgm3
      mpj = volj * densp_kgm3
      vpibar = sqrt ( vpconst * tempk / mpi )
      vpjbar = sqrt ( vpconst * tempk / mpj )
      lambdapi = 2.0 * dpi / ( pi * vpibar ) 
      lambdapj = 2.0 * dpj / ( pi * vpjbar ) 
      deltai = (2.0*ri + lambdapi)**3 - (4.0*ri*ri + lambdapi*lambdapi)**1.5
      deltai = deltai / ( 6.0*ri*lambdapi ) - 2.0*ri  
      deltaj = (2.0*rj + lambdapj)**3 - (4.0*rj*rj + lambdapj*lambdapj)**1.5
      deltaj = deltaj / ( 6.0*rj*lambdapj ) - 2.0*rj  
      denom1 = (ri+rj) / ( (ri+rj) + sqrt(deltai*deltai + deltaj*deltaj) )
      denom2 = 4.0*(dpi+dpj) / ( (ri+rj) * sqrt(vpibar*vpibar + vpjbar*vpjbar) )
      betaij = fourpi * ( ri + rj ) * ( dpi + dpj ) / ( denom1 + denom2 )
!-------------------------------------------------------------------------------
!     write(aunit2,'(12d11.3)') pres,rhoa,lambdaa,di,dj,knai,knaj,gi,gj,dpi,dpj,betaij
!     write(*,     '(12d11.3)') pres,rhoa,lambdaa,di,dj,knai,knaj,gi,gj,dpi,dpj,betaij
!     write(aunit2,*)'  '
!     write(aunit2,*)'di,dj,tempk,pres'
!     write(aunit2,*) di,dj,tempk,pres
!     write(aunit2,*)'  '
!     write(aunit2,*)'etaa,vabar,rhoa,lambdaa'
!     write(aunit2,*) etaa,vabar,rhoa,lambdaa 
!     write(aunit2,*)'  '
!     write(aunit2,*)'knai,knaj,gi,gj,dpi,dpj'
!     write(aunit2,*) knai,knaj,gi,gj,dpi,dpj 
!     write(aunit2,*)'  '
!     write(aunit2,*)'mpi,mpj,voli,volj,vpibar,vpjbar'
!     write(aunit2,*) mpi,mpj,voli,volj,vpibar,vpjbar 
!     write(aunit2,*)'  '
!     write(aunit2,*)'lambdapi,lambdapj,deltai,deltaj,denom1,denom2'
!     write(aunit2,*) lambdapi,lambdapj,deltai,deltaj,denom1,denom2 
!     write(aunit2,*)'  '
!     write(aunit2,*)'betaij=',betaij
!-------------------------------------------------------------------------------
      return
      end subroutine Brownian_Coag_Coef


!>-------------------------------------------------------------------------------------------------------------------------------------
!!     routine to calculate the convective brownian diffusion enhancement coagulation coefficient.
!!     see mzj, 2005, p. 510, eq.15.35.
!!-------------------------------------------------------------------------------------------------------------------------------------
      subroutine Cbde_Coag_Coef( di, dj, tempk, pres, kdeij )
      implicit none
      real(4) :: di      !< particle diameters [um]
      real(4) :: dj      !< particle diameters [um]
      real(4) :: tempk   !< ambient temperature [k]
      real(4) :: pres    !< ambient pressure [pa]
      real(4) :: kdeij   !< coagulation coefficient [m^3/s/particle]
      
      real(4), parameter :: pi = 3.141592653589793
      real(4), parameter :: sixpi  = 6.0 * pi
      real(4), parameter :: kb = 1.38065e-23     ! boltzmann's constant [j/k]
      real(4), parameter :: na = 6.0221367e+23   ! avogadro's number [#/mole]
      real(4), parameter :: mwair = 28.9628e-03  ! molecular weight of air in [kg/mole]
      real(4), parameter :: m1air = mwair/na     ! average mass of one air molecule [kg]
      real(4), parameter :: vaconst = 8.0*kb/pi/m1air  ! used in eq.(16.21) for vabar
      real(4), parameter :: aprime = 1.249       ! for cunningham slip-flow correction
      real(4), parameter :: bprime = 0.42        !   as given in eq.(16.25);
      real(4), parameter :: cprime = 0.87        !   kasten (1968) values used here
      real(4), parameter :: densp_kgm3 = 1.4e+03 ! particle density [kg/m^3]
      real(4), parameter :: ggrav = 9.81         ! gravitational acceleration [m/s^2]
      real(4) :: ri, rj                          ! particle radii [um]
      real(4) :: dpi                             ! particle diffusion coefficients [m^2/s]
      real(4) :: gi, gj                          ! cunningham slip-flow corrections for particles i, j [1]
      real(4) :: etaa                            ! dynamic viscosity of air [kg/m/s]
      real(4) :: vabar                           ! mean thermal velocity of an air molecule [m/s]
      real(4) :: rhoa                            ! density of air [kg/m^3]
      real(4) :: lambdaa                         ! mean free path of air [m]
      real(4) :: knai, knaj                      ! knudsen numbers in air for particles i, j [1]
      real(4) :: nua                             ! kinematic viscosity of air [m^2/s]
      real(4) :: vfj                             ! terminal fall speed [m/s]
      real(4) :: rej                             ! particle reynolds number [1]
      real(4) :: scpi                            ! particle schmidt  number [1]
      real(4) :: kbij                            ! brownian diffusion coefficient [m^3/s]

      !-----------------------------------------------------------------------------------------------------------------
      ! rj must be larger than or equal to ri for this coagulation mechanism.
      ! make rj greater than or equal ri for all pairs of particles. 
      !-----------------------------------------------------------------------------------------------------------------
      if( dj .ge. di ) then 
        ri = 0.5e-06 * di                                                     ! convert from [um] to [m]
        rj = 0.5e-06 * dj                                                     ! convert from [um] to [m]
      else
        rj = 0.5e-06 * di                                                     ! convert from [um] to [m]
        ri = 0.5e-06 * dj                                                     ! convert from [um] to [m]
      endif
      !-----------------------------------------------------------------------------------------------------------------
      ! dynamic viscosity of air [kg/m/s], jacobson, 2005, eq.(4.54).
      ! mean thermal velocity of an air molecule [m/s].
      ! density of air [kg/m^3] from the ideal gas law.
      ! mean free path of an air molecule [m], jacobson, 2005, eq.(15.24).
      !-----------------------------------------------------------------------------------------------------------------
      etaa = 1.8325e-05 * (416.16/(tempk + 120.0)) * (tempk/296.16)**1.5      ! [kg/m/s]
      vabar = sqrt( vaconst * tempk )                                         ! [m/s]
      rhoa = pres * m1air / ( kb * tempk )                                    ! [kg/m^3]
      lambdaa = 2.0 * etaa / ( rhoa * vabar )                                 ! [m]
      knai = lambdaa / ri                                                     ! [1]
      knaj = lambdaa / rj                                                     ! [1]
      gi = 1.0 + knai * ( aprime + bprime * exp(-cprime/knai) )               ! cunningham slip-flow correction [1]
      gj = 1.0 + knaj * ( aprime + bprime * exp(-cprime/knaj) )               ! cunningham slip-flow correction [1]
      dpi = kb * tempk * gi / ( sixpi * ri * etaa )                           ! particle diffusion coefficient  [m^2/s]
      nua = etaa / rhoa                                                       ! kinematic viscosity of air [m^2/s]
      vfj = 2.0 * rj*rj*(densp_kgm3-rhoa)*ggrav*gj/(9.0*etaa)                 ! fall velocity [m/s] 
      rej = 2.0 * rj * vfj / nua                                              ! particle reynolds number [1] 
      scpi = nua / dpi                                                        ! particle schmidt  number [1] 
      call brownian_coag_coef ( di, dj, tempk, pres, kbij )                   ! brownian coag. coef. [m^3/s]
      if( rej .le. 1.0 ) kdeij = 0.45 * kbij * rej**0.3333333 * scpi**0.3333333   ! cbde coag. coef. [m^3/s]
      if( rej .gt. 1.0 ) kdeij = 0.45 * kbij * rej**0.5000000 * scpi**0.3333333   ! cbde coag. coef. [m^3/s]
      ! write(*,*) 'vfj = ', vfj
      return
      end subroutine Cbde_Coag_Coef


!>-------------------------------------------------------------------------------------------------------------------------------------
!!     routine to calculate the gravitational collection coagulation coefficient. 
!!     see mzj, 2005, p. 510, eq.15.37.
!!-------------------------------------------------------------------------------------------------------------------------------------
      subroutine Gravcoll_Coag_Coef( di, dj, tempk, pres, kgcij )
      implicit none
      real(4) :: di     !< particle diameters [um]
      real(4) :: dj     !< particle diameters [um]
      real(4) :: tempk  !< ambient temperature [k]
      real(4) :: pres   !< ambient pressure [pa]
      real(4) :: kgcij  !< coagulation coefficient [m^3/s/particle]
      
      real(4), parameter :: pi = 3.141592653589793
      real(4), parameter :: sixpi  = 6.0 * pi
      real(4), parameter :: kb = 1.38065e-23     ! boltzmann's constant [j/k]
      real(4), parameter :: na = 6.0221367e+23   ! avogadro's number [#/mole]
      real(4), parameter :: mwair = 28.9628e-03  ! molecular weight of air in [kg/mole]
      real(4), parameter :: m1air = mwair/na     ! average mass of one air molecule [kg]
      real(4), parameter :: vaconst = 8.0*kb/pi/m1air  ! used in eq.(16.21) for vabar
      real(4), parameter :: aprime = 1.249       ! for cunningham slip-flow correction
      real(4), parameter :: bprime = 0.42        !   as given in eq.(16.25);
      real(4), parameter :: cprime = 0.87        !   kasten (1968) values used here
      real(4), parameter :: densp_kgm3 = 1.4e+03 ! particle density [kg/m^3]
      real(4), parameter :: ggrav = 9.81         ! gravitational acceleration [m/s^2]
      real(4) :: ri, rj                          ! particle radii [um]
      real(4) :: gi, gj                          ! cunningham slip-flow corrections for particles i, j [1]
      real(4) :: etaa                            ! dynamic viscosity of air [kg/m/s]
      real(4) :: vabar                           ! mean thermal velocity of an air molecule [m/s]
      real(4) :: rhoa                            ! density of air [kg/m^3]
      real(4) :: lambdaa                         ! mean free path of air [m]
      real(4) :: knai, knaj                      ! knudsen numbers in air for particles i, j [1]
      real(4) :: vfi, vfj                        ! terminal fall speed [m/s]
      real(4) :: ecollij                         ! collision efficiency [1]
      real(4) :: p                               ! p = min(ri,rj)/max(ri,rj) from pruppacher and klett 1980, p.377.

      ri = 0.5e-06 * di                                                       ! convert from [um] to [m]
      rj = 0.5e-06 * dj                                                       ! convert from [um] to [m]
      !-----------------------------------------------------------------------------------------------------------------
      ! dynamic viscosity of air [kg/m/s], jacobson, 2005, eq.(4.54).
      ! mean thermal velocity of an air molecule [m/s].
      ! density of air [kg/m^3] from the ideal gas law.
      ! mean free path of an air molecule [m], jacobson, 2005, eq.(15.24).
      !-----------------------------------------------------------------------------------------------------------------
      etaa = 1.8325e-05 * (416.16/(tempk + 120.0)) * (tempk/296.16)**1.5      ! [kg/m/s]
      vabar = sqrt( vaconst * tempk )                                         ! [m/s]
      rhoa = pres * m1air / ( kb * tempk )                                    ! [kg/m^3]
      lambdaa = 2.0 * etaa / ( rhoa * vabar )                                 ! [m]
      knai = lambdaa / ri                                                     ! [1]
      knaj = lambdaa / rj                                                     ! [1]
      gi = 1.0 + knai * ( aprime + bprime * exp(-cprime/knai) )               ! cunningham slip-flow correction [1]
      gj = 1.0 + knaj * ( aprime + bprime * exp(-cprime/knaj) )               ! cunningham slip-flow correction [1]
      vfi = 2.0 * ri*ri*(densp_kgm3-rhoa)*ggrav*gi/(9.0*etaa)                 ! fall velocity [m/s] 
      vfj = 2.0 * rj*rj*(densp_kgm3-rhoa)*ggrav*gj/(9.0*etaa)                 ! fall velocity [m/s] 
!-------------------------------------------------------------------------------------------------------------------------
!     taking the collision efficiency from pruppacher and klett, 1980, eq.12-78.
!-------------------------------------------------------------------------------------------------------------------------
      if( rj .gt. ri ) then
        p = ri / rj
      else
        p = rj / ri
      endif
      ecollij = 0.5 * ( p / ( 1.0 + p ) )**2                                  ! collision efficiency [1], eq.12-78
      kgcij = ecollij * pi * ( ri + rj )**2 * abs( vfi - vfj )                ! grav. collection coag. coef. [m^3/s]
!-------------------------------------------------------------------------------------------------------------------------
!     write(*,'(/,a    )')'kgcij, ri, rj, vfi, vfj, ecollij, stij, rej, evij, eaij'
!     write(*,'(10d12.4)') kgcij, ri, rj, vfi, vfj, ecollij, stij, rej, evij, eaij
!     write(*,*) 'nua = ', nua
!-------------------------------------------------------------------------------------------------------------------------
      return
      end subroutine Gravcoll_Coag_Coef


!>-------------------------------------------------------------------------------------------------------------------------------------
!!     routine to calculate the turbulent inertial coagulation coefficient. see mzj, 2005, p. 511, eq.15.40.
!!     routine to calculate the turbulent shear    coagulation coefficient. see mzj, 2005, p. 511, eq.15.41.
!!-------------------------------------------------------------------------------------------------------------------------------------
      subroutine Turb_Coag_Coef( di, dj, tempk, pres, ktiij, ktsij )
      implicit none
      real(4) :: di     !< particle diameters [um]
      real(4) :: dj     !< particle diameters [um]
      real(4) :: tempk  !< ambient temperature [k]
      real(4) :: pres   !< ambient pressure [pa]
      real(4) :: ktiij  !< turbulent inertial coagulation coefficient [m^3/s/particle]
      real(4) :: ktsij  !< turbulent shear    coagulation coefficient [m^3/s/particle]
      
      real(4), parameter :: pi = 3.141592653589793
      real(4), parameter :: kb = 1.38065e-23     ! boltzmann's constant [j/k]
      real(4), parameter :: na = 6.0221367e+23   ! avogadro's number [#/mole]
      real(4), parameter :: mwair = 28.9628e-03  ! molecular weight of air in [kg/mole]
      real(4), parameter :: m1air = mwair/na     ! average mass of one air molecule [kg]
      real(4), parameter :: vaconst = 8.0*kb/pi/m1air  ! used in eq.(16.21) for vabar
      real(4), parameter :: aprime = 1.249       ! for cunningham slip-flow correction
      real(4), parameter :: bprime = 0.42        !   as given in eq.(16.25);
      real(4), parameter :: cprime = 0.87        !   kasten (1968) values used here
      real(4), parameter :: densp_kgm3 = 1.4e+03 ! particle density [kg/m^3]
      real(4), parameter :: ggrav = 9.81         ! gravitational acceleration [m/s^2]
      real(4), parameter :: ed = 5.0e-04         ! dissipation rate of turbulent kinetic energy per gram of medium [m^2/s^3]
                                                 ! typical clear-sky value from mzj (2005, p.511,p.236) of 5.0e-05 [m^2/s^3]
                                                 !   from pruppacher and klett (1997).
      real(4) :: ri, rj                          ! particle radii [um]
      real(4) :: gi, gj                          ! cunningham slip-flow corrections for particles i, j [1]
      real(4) :: etaa                            ! dynamic viscosity of air [kg/m/s]
      real(4) :: vabar                           ! mean thermal velocity of an air molecule [m/s]
      real(4) :: rhoa                            ! density of air [kg/m^3]
      real(4) :: lambdaa                         ! mean free path of air [m]
      real(4) :: knai, knaj                      ! knudsen numbers in air for particles i, j [1]
      real(4) :: nua                             ! kinematic viscosity of air [m^2/s]
      real(4) :: vfi, vfj                        ! terminal fall speed [m/s]

      ri = 0.5e-06 * di                                                       ! convert from [um] to [m]
      rj = 0.5e-06 * dj                                                       ! convert from [um] to [m]
      !-----------------------------------------------------------------------------------------------------------------
      ! dynamic viscosity of air [kg/m/s], jacobson, 2005, eq.(4.54).
      ! mean thermal velocity of an air molecule [m/s].
      ! density of air [kg/m^3] from the ideal gas law.
      ! mean free path of an air molecule [m], jacobson, 2005, eq.(15.24).
      !-----------------------------------------------------------------------------------------------------------------
      etaa = 1.8325e-05 * (416.16/(tempk + 120.0)) * (tempk/296.16)**1.5      ! [kg/m/s]
      vabar = sqrt( vaconst * tempk )                                         ! [m/s]
      rhoa = pres * m1air / ( kb * tempk )                                    ! [kg/m^3]
      lambdaa = 2.0 * etaa / ( rhoa * vabar )                                 ! [m]
      knai = lambdaa / ri                                                     ! [1]
      knaj = lambdaa / rj                                                     ! [1]
      gi = 1.0 + knai * ( aprime + bprime * exp(-cprime/knai) )               ! cunningham slip-flow correction [1]
      gj = 1.0 + knaj * ( aprime + bprime * exp(-cprime/knaj) )               ! cunningham slip-flow correction [1]
      vfi = 2.0 * ri*ri*(densp_kgm3-rhoa)*ggrav*gi/(9.0*etaa)                 ! fall velocity [m/s] 
      vfj = 2.0 * rj*rj*(densp_kgm3-rhoa)*ggrav*gj/(9.0*etaa)                 ! fall velocity [m/s] 
      nua = etaa / rhoa                                                       ! kinematic viscosity of air [m^2/s]
      ktiij = (pi*ed**0.75/(ggrav*nua**0.25))*(ri+rj)**2*abs(vfi-vfj)         ! turbulent inertial coag. coef. [m^3/s]
      ktsij = sqrt(8.0*pi*ed/(15.0*nua))*(ri+rj)**3                           ! turbulent shear    coag. coef. [m^3/s]
      ! write(*,'(10d12.4)') ktiij, ktsij, ri, rj, vfi, vfj
      return
      end subroutine Turb_Coag_Coef


!>-------------------------------------------------------------------------------------------------------------------------------------
!!     routine to calculate the total coagulation coefficient. 
!!-------------------------------------------------------------------------------------------------------------------------------------
      subroutine total_coag_coef( di, dj, tempk, pres, ktotij )
      implicit none
      real(4) :: di     !< particle diameters [um]
      real(4) :: dj     !< particle diameters [um]
      real(4) :: tempk  !< ambient temperature [k]
      real(4) :: pres   !< ambient pressure [pa]
      real(4) :: ktotij !< total  kernel [m^3/s/particle]

      real(4) :: kbij                            ! brownian                                  kernel [m^3/s/particle]
      real(4) :: kdeij                           ! convective brownian diffusion enhancement kernel [m^3/s/particle]
      real(4) :: kgcij                           ! gravitational collection                  kernel [m^3/s/particle]
      real(4) :: ktiij                           ! turbulent inertial                        kernel [m^3/s/particle]
      real(4) :: ktsij                           ! turbulent shear                           kernel [m^3/s/particle]
      call Brownian_Coag_Coef( di, dj, tempk, pres, kbij  )
      call Cbde_Coag_Coef    ( di, dj, tempk, pres, kdeij )
      call Gravcoll_Coag_Coef( di, dj, tempk, pres, kgcij )
      call Turb_Coag_Coef    ( di, dj, tempk, pres, ktiij, ktsij )
      ktotij = kbij + kdeij + kgcij + ktiij + ktsij
      return
      end subroutine Total_Coag_Coef


!>-----------------------------------------------------------------------------------------------------------------------
!!     routine to test the various routines for coagulation coefficients.
!!-----------------------------------------------------------------------------------------------------------------------
      subroutine Test_Coag_Coef
      implicit none
      integer, parameter :: nsizes =   81         ! was 501 for mzj figure 15.7
      real(4), parameter :: dlower =  0.0030      ! [um]
      real(4), parameter :: dupper = 30.0000      ! [um]
      real(4) :: d(nsizes)                        ! [um]
      real(4) :: drat,dtest1,dtest2,beta          ! scratch variables
      real(4) :: tempk                            ! temperature [k]
      real(4) :: pres                             ! pressure [pa]
      real(4) :: kbij                             ! brownian                                  kernel [m^3/s/particle]
      real(4) :: kdeij                            ! convective brownian diffusion enhancement kernel [m^3/s/particle]
      real(4) :: kgcij                            ! gravitational collection                  kernel [m^3/s/particle]
      real(4) :: ktiij                            ! turbulent inertial                        kernel [m^3/s/particle]
      real(4) :: ktsij                            ! turbulent shear                           kernel [m^3/s/particle]
      real(4) :: ktotij, ktotij_tmp               ! total                                     kernel [m^3/s/particle]
      integer :: i,j,n
      integer, parameter :: size_number = 1
      logical, parameter :: mzj_figure    = .false.
      logical, parameter :: tp_dependence = .true.
      logical, parameter :: write_table   = .false.

      open(aunit2,file='kij.out',status='replace')

      ! write(aunit2,*)'n,d(n)'
      drat = (dupper/dlower)**(1.0/real(nsizes-1))
      do n=1,nsizes
        d(n) = dlower*(drat**(n-1))
        ! write(aunit2,90000)n,d(n)
      enddo
!---------------------------------------------------------------------------------------------------------------------------
!     for comparison with figure 16.4, p.448 of mzj 1999.
!     for comparison with figure 15.7, p.512 of mzj 2005.
!     dlw, 082906: checked against fig.15.7 of mzj 2005; gives close agreement.
!---------------------------------------------------------------------------------------------------------------------------
      if( mzj_figure ) then
        write(aunit2,'(/A/)')'    N   RTEST[um]    R(N)[um]  KBIJ[cm3/s] KDEIJ[cm3/s] KGCIJ[cm3/s] KTIIJ[cm3/s] KTSIJ[cm3/s]'
        pres = 101325.0
        tempk = 298.0  
        if( size_number .eq. 1 ) dtest1 =  0.02   ! for r= 0.01 um test
        if( size_number .eq. 2 ) dtest1 = 20.00   ! for r=10.00 um test
        if( size_number .eq. 3 ) dtest1 =  2.00   ! for r= 1.00 um test
        do n=1,nsizes
          call Brownian_Coag_Coef( dtest1, d(n), tempk, pres, kbij  )
          call Cbde_Coag_Coef    ( dtest1, d(n), tempk, pres, kdeij )
          call Gravcoll_Coag_Coef( d(n), dtest1, tempk, pres, kgcij )
          call Turb_Coag_Coef    ( dtest1, d(n), tempk, pres, ktiij, ktsij )
          call Total_Coag_Coef   ( dtest1, d(n), tempk, pres, ktotij )
          ktotij_tmp = kbij + kdeij + kgcij + ktiij + ktsij                     ! for check of routine total_coag_coef.
          write(aunit2,90001)n,0.5*dtest1,0.5*d(n),1.0e+06*kbij,1.0e+06*kdeij,1.0e+06*kgcij,1.0e+06*ktiij,1.0e+06*ktsij, &
                                                   1.0e+06*ktotij, 1.0e+06*ktotij_tmp 
        enddo
      endif
!---------------------------------------------------------------------------------------------------------------------------
!     for examination of the temperature and pressure dependence.
!---------------------------------------------------------------------------------------------------------------------------
      if( tp_dependence ) then
        write(aunit2,*)'N,DTEST1,DTEST2,BROWNIAN_BETA,TEMPK,PRES'
        pres = 101325.0
        tempk = 288.0
        dtest1 =  0.003
        dtest2 = 30.000
        do n=1,14
          call Brownian_Coag_Coef(dtest1,dtest2,tempk,pres,beta)
          write(aunit2,90008)n,dtest1,dtest2,1.0e+06*beta,tempk,pres
          ! tempk = tempk + 10.0
          pres = 0.70*pres
        enddo
      endif
!---------------------------------------------------------------------------------------------------------------------------
!     table of coagulation coefficients. 
!---------------------------------------------------------------------------------------------------------------------------
      if( write_table ) then
        write(aunit2,'(2A)')'I, J, R(I)[um], R(J)[um],',' KTOTIJ[cm^3/s]'
        pres = 101325.0   ! for examination of the pressure    dependence: *0.1,  *0.01
        tempk = 288.00    ! for examination of the temperature dependence: 200.0, 325.0 
        do j=1,nsizes
        do i=1,nsizes
          call Total_Coag_Coef ( d(i), d(j), tempk, pres, ktotij )
          write(aunit2,90002) i, j, 0.5*d(i), 0.5*d(j), 1.0e+06*ktotij
        enddo
        enddo
      endif
!---------------------------------------------------------------------------------------------------------------------------
      close(aunit2)
90000 format(i5,f15.6)
90001 format(i4,2f10.5,7e13.4)
90002 format(2i5,2f12.6,2e16.6)
90008 format(i5,2f12.6,f16.6,2f12.2)
      return
      end subroutine Test_Coag_Coef


      END MODULE Aero_Coag

