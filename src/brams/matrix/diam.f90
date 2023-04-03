  MODULE Aero_Diam
      use memMatrix,  only: &
                              nlays,  &         !intent()
                              aunit1, &         !intent()
                              write_log, &         !intent()
                              nmodes, &         !intent()
                              mech, &              !intent()
                              n_dp_condtable, & !intent()
                              dp_condtable,   & !intent()
                              mode_name,      &
                              diam     ,      &         !intent()
			      diam_histogram
!-------------------------------------------------------------------------------------------------------------------------
!     The array DIAM(x,y,z) contains current values of some measure of average ambient mode diameter for each mode  
!       for use outside of the MATRIX microphysical module where it is calculated. 
!       Values in DIAM are saved at the top of the subr. MATRIX before microphysical evolution 
!       for the current time step is done. 
!
!     The current measure of particle diameter is the diameter of average mass:
! 
!        DIAM(x,y,z) = [ (6/pi) * (Mi/Ni) * (1/D) ]^(1/3)
!
!     with Mi the total mass concentration (including water) in mode i, Ni the number concentration in mode i, and
!     D is the particle density. 
!-------------------------------------------------------------------------------------------------------------------------
      implicit none
      
      CONTAINS

!>-------------------------------------------------------------------------------------------------------------------------
!!@brief     routine to initialize diam before model time stepping is begun.
!!     this routine should be called only once. 
!!     the diameter of average mass is the cube root of the normalized third diameter moment: (m3/m0)**(1/3) = dg*sg**3
!!@author susanne bauer/doug wright
!-------------------------------------------------------------------------------------------------------------------------
      subroutine Setup_Diam
      use memMatrix, only: &
                           dgn0, & !intent()
                           sig0    !intent()
      integer :: i
      real(8) :: sg, d_number_mean
      if( write_log ) write(aunit1,'(/A/)') 'I,DGN0(I),SIG0(I),SG,DIAM(1,1,1,I)*1.0D+06,D_NUMBER_MEAN'
      do i=1, nmodes
        sg = exp( 0.5d+00*( log(sig0(i)) )**2 )
        diam(:,:,:,i) = 1.0d-06 * dgn0(i)*sg**3           ! convert from [um] to [m]
        d_number_mean =           dgn0(i)*sg             
        if( write_log ) write(aunit1,90) i,dgn0(i),sig0(i),sg,diam(1,1,1,i)*1.0d+06,d_number_mean
      enddo
      write(aunit1,'(a)') '  '

      ! zero histogram for diam values.
      
      diam_histogram(:,:,:) = 0.0d+00 

90    format(i5,6f12.5)
      return
      end subroutine Setup_Diam


!>-------------------------------------------------------------------------------------------------------------------------
!!     call this routine at the end of all time stepping.
!!-------------------------------------------------------------------------------------------------------------------------
      subroutine Write_Diam_Histogram
      integer :: i, n
      integer, parameter :: outunit  = 93
      integer, parameter :: nbinssum =  5
      real(8) :: sum_hist(nmodes),sum_hist_dry(nmodes),d,f16,f18,f20

      select case( mech ) 
      case( 1 )
        OPEN(OUTUNIT,STATUS='REPLACE',FILE='mech1_diam_histogram.plt')
      CASE( 2 )
        OPEN(OUTUNIT,STATUS='REPLACE',FILE='mech2_diam_histogram.plt')
      CASE( 3 )
        OPEN(OUTUNIT,STATUS='REPLACE',FILE='mech3_diam_histogram.plt')
      CASE( 4 )
        OPEN(OUTUNIT,STATUS='REPLACE',FILE='mech4_diam_histogram.plt')
      CASE( 5 )
        OPEN(OUTUNIT,STATUS='REPLACE',FILE='mech5_diam_histogram.plt')
      CASE( 6 )
        OPEN(OUTUNIT,STATUS='REPLACE',FILE='mech6_diam_histogram.plt')
      CASE( 7 )
        OPEN(OUTUNIT,STATUS='REPLACE',FILE='mech7_diam_histogram.plt')
      CASE( 8 )
        OPEN(OUTUNIT,STATUS='REPLACE',FILE='mech8_diam_histogram.plt')
      CASE DEFAULT
        OPEN(OUTUNIT,STATUS='REPLACE',FILE='mechx_diam_histogram.plt')
      END SELECT 

      !-------------------------------------------------------------------------------------------------------------------
      ! f16, f18, and f20 convert the diameter of average mass (d) to the number mean diameter for 
      ! geometric standard deviations of 1.6, 1.8, and 2.0, respectively. 
      !-------------------------------------------------------------------------------------------------------------------
      f16 = 1.0d+00 / ( exp( ( log(1.6d+00) )**2 ) ) 
      f18 = 1.0d+00 / ( exp( ( log(1.8d+00) )**2 ) ) 
      f20 = 1.0d+00 / ( exp( ( log(2.0d+00) )**2 ) ) 

      !-------------------------------------------------------------------------------------------------------------------
      ! for each mode, normalize the histogram to unity. 
      !-------------------------------------------------------------------------------------------------------------------
      do i=1, nmodes
        sum_hist(i) = 0.0d+00
        do n=1, n_dp_condtable
          sum_hist(i) = sum_hist(i) + diam_histogram(i,n,1)
        enddo
        diam_histogram(i,:,1) = diam_histogram(i,:,1) / ( sum_hist(i) + 1.0d-30 ) 
        diam_histogram(i,:,2) = diam_histogram(i,:,2) / ( sum_hist(i) + 1.0d-30 ) 
        write(outunit,'(A,I3,1X,A,F15.1)') 'Total count for mode I = ',I,' is ', sum_hist(i)
      ENDDO

      write(outunit,91) '         DAM','       DNM16','       DNM18','       DNM20', &
                                       '       DGN16','       DGN18','       DGN20',mode_name(:),mode_name(:)
      sum_hist    (:) = 0.0d+00
      sum_hist_dry(:) = 0.0d+00
      d = 1.0d+00
      n = 0
      do i=1, n_dp_condtable
        sum_hist    (1:nmodes) = sum_hist    (1:nmodes) + diam_histogram(1:nmodes,i,1) + 1.0d-20
        sum_hist_dry(1:nmodes) = sum_hist_dry(1:nmodes) + diam_histogram(1:nmodes,i,2) + 1.0d-20
        d = d*dp_condtable(i)
        if( mod( i, nbinssum ) .eq. 0 ) then
          n = n + 1
          d = 1.0d+06 * d**(1.0d+00/real(nbinssum))
          write(outunit,90)n,i,d,d*f16,d*f18,d*f20,d*f16**1.5,d*f18**1.5,d*f20**1.5,sum_hist(1:nmodes),sum_hist_dry(1:nmodes)
          sum_hist    (:) = 0.0d+00
          sum_hist_dry(:) = 0.0d+00
          d = 1.0d+00
        endif
      enddo
      close(outunit)
90    format(2i5,7f12.7,32f14.10)
91    format(10x,7a12,32a14)
      return
      end subroutine Write_Diam_Histogram

   END MODULE Aero_Diam
 
