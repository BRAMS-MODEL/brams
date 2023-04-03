!>
!!@brief     This module contains various aerosol microphysical routines.
!!@author    Susanne Bauer/Doug Wright
!!----------------------------------------------------------------------------------------------------------------------
   MODULE Aero_Subs
   
      USE memMatrix, ONLY: &
         naerobox, &
         ngases ,&
         nmass_spcs, &
         nemis_spcs, &
         tinydenom,          & !intent()
         mass_no3,           & !intent()
         mass_nh4,           & !intent()
         gas_h2so4,          & !intent()
         gas_hno3,           & !intent()
         gas_nh3,            & !intent()
         write_log,          & !intent()
         ngases,             & !intent()
         nmass_spcs,         & !intent()
         nemis_spcs,         & !intent()
         conv_mass_to_dp,    & !intent()
         no_microphysics,    & !intent()
         update_kij,         & !intent()
         numb_akk_1,         & !intent()
         dg_akk,             & !intent()
         dg_acc,             & !intent()
         sg_akk,             & !intent()
         sg_acc,             & !intent()
         aunit1,             & !intent()
         ixxx,               & !intent()
         iyyy,               & !intent()
         mass_adj,           & !intent()
         mass_ssc_seas,      & !intent()
         mass_ssa_seas,      & !intent()
         mass_ssc_sulf,      & !intent()
         mass_ssa_sulf,      & !intent()
         tinynumer,          & !intent()
         mass_h2o,           & !intent()
         update_dp,          & !intent()
         conv_vol_to_dp_fac, & !intent()
         dpmin_global,       & !intent()
         dpmax_global,       & !intent()
         ilay,               & !intent()
         update_vdep,        & !intent()
         do_npf,             & !intent()
         convnh3,            & !intent()
         ugm3_ncm3,          & !intent()
         aqso4rate_min,      & !intent()
         activation_scheme,  & !intent()
         prod_index_sulf,    & !intent()
         minconc,            & !intent()
         mass_dd1_dust,      & !intent()
         mass_bc1_bcar,      & !intent()
         mass_bc2_bcar,      & !intent()
         mass_dd2_dust,      & !intent()
         mass_dd1_sulf,      & !intent()
         mimr_ddd,           & !intent()
         mass_ds1_sulf,      & !intent()
         mass_ds1_dust,      & !intent()
         numb_ds1_1,         & !intent()
         numb_dd1_1,         & !intent()
         mass_dd2_sulf,      & !intent()
         mass_ds2_sulf,      & !intent()
         mass_ds2_dust,      & !intent()
         numb_ds2_1,         & !intent()
         numb_dd2_1,         & !intent()
         mass_bc1_sulf,      & !intent()
         mimr_bc1,           & !intent()
         mass_bc2_sulf,      & !intent()
         numb_bc2_1,         & !intent()
         numb_bc1_1,         & !intent()
         mimr_bc2,           & !intent()
         include_bc3,        & !intent()
         mass_bc3_sulf,      & !intent()
         mass_bc3_bcar,      & !intent()
         numb_bc3_1,         & !intent()
         mivf_ddd,           & !intent()
         mivf_bc1,           & !intent()
         mivf_bc2,           & !intent()
         mass_akk_sulf,      & !intent()
         mass_acc_sulf,      & !intent()
         numb_acc_1,         & !intent()
         imtr_method,        & !intent()
         ndiag_aero,         & !intent()
         numb_acc_1,         & !intent()
         imtr_method,        & !intent()
         mw_so4,             & !intent()
         avo,                & !intent()
         pi,                 & !intent()
         ndiag_aero            !intent()
         
         

!LFR>       use Aero_Param, ONLY: &  
!LFR>                             tinydenom,          & !intent()  
!LFR>                             mass_no3,           & !intent()  
!LFR>                             mass_nh4,           & !intent()  
!LFR>                             gas_h2so4,          & !intent()  
!LFR>                             gas_hno3,           & !intent()  
!LFR>                             gas_nh3,            & !intent()  
!LFR>                             write_log,          & !intent()  
!LFR>                             ngases,             & !intent()  
!LFR>                             nmass_spcs,         & !intent()  
!LFR>                             nemis_spcs,         & !intent()        
!LFR>                             conv_mass_to_dp,    & !intent()  
!LFR>                             no_microphysics,    & !intent()  
!LFR>                             update_kij,         & !intent()  
!LFR>                             numb_akk_1,         & !intent()  
!LFR>                             dg_akk,             & !intent()  
!LFR>                             dg_acc,             & !intent()  
!LFR>                             sg_akk,             & !intent()  
!LFR>                             sg_acc,             & !intent()  
!LFR>                             aunit1,             & !intent()  
!LFR>                             ixxx,               & !intent()  
!LFR>                             iyyy,               & !intent()  
!LFR>                             mass_adj,           & !intent()  
!LFR>                             mass_ssc_seas,      & !intent()  
!LFR>                             mass_ssa_seas,      & !intent()  
!LFR>                             mass_ssc_sulf,      & !intent()  
!LFR>                             mass_ssa_sulf,      & !intent()  
!LFR>                             tinynumer,          & !intent()  
!LFR>                             mass_h2o,           & !intent()  
!LFR>                             update_dp,          & !intent()  
!LFR>                             conv_vol_to_dp_fac, & !intent()  
!LFR>                             dpmin_global,       & !intent()  
!LFR>                             dpmax_global,       & !intent()  
!LFR>                             ilay,               & !intent()  
!LFR>                             update_vdep,        & !intent()  
!LFR>                             do_npf,             & !intent()  
!LFR>                             convnh3,            & !intent()  
!LFR>                             ugm3_ncm3,          & !intent()  
!LFR>                             aqso4rate_min,      & !intent()  
!LFR>                             activation_scheme,  & !intent()  
!LFR>                             prod_index_sulf,    & !intent()  
!LFR>                             minconc,            & !intent()  
!LFR>                             mass_dd1_dust,      & !intent()  
!LFR>                             mass_bc1_bcar,      & !intent()  
!LFR>                             mass_bc2_bcar,      & !intent()  
!LFR>                             mass_dd2_dust,      & !intent()  
!LFR>                             mass_dd1_sulf,      & !intent()  
!LFR>                             mimr_ddd,           & !intent()  
!LFR>                             mass_ds1_sulf,      & !intent()  
!LFR>                             mass_ds1_dust,      & !intent()  
!LFR>                             numb_ds1_1,         & !intent()  
!LFR>                             numb_dd1_1,         & !intent()  
!LFR>                             mass_dd2_sulf,      & !intent()  
!LFR>                             mass_ds2_sulf,      & !intent()  
!LFR>                             mass_ds2_dust,      & !intent()  
!LFR>                             numb_ds2_1,         & !intent()  
!LFR>                             numb_dd2_1,         & !intent()  
!LFR>                             mass_bc1_sulf,      & !intent()  
!LFR>                             mimr_bc1,           & !intent()  
!LFR>                             mass_bc2_sulf,      & !intent()  
!LFR>                             numb_bc2_1,         & !intent()  
!LFR>                             numb_bc1_1,         & !intent()  
!LFR>                             mimr_bc2,           & !intent()  
!LFR>                             include_bc3,        & !intent()  
!LFR>                             mass_bc3_sulf,      & !intent()  
!LFR>                             mass_bc3_bcar,      & !intent()  
!LFR>                             numb_bc3_1,         & !intent()  
!LFR>                             mivf_ddd,           & !intent()  
!LFR>                             mivf_bc1,           & !intent()  
!LFR>                             mivf_bc2,           & !intent()  
!LFR>                             mass_akk_sulf,      & !intent()  
!LFR>                             mass_acc_sulf,      & !intent()  
!LFR>                             numb_acc_1,         & !intent()  
!LFR>                             imtr_method,        & !intent()  
!LFR>                             ndiag_aero,         & !intent()  
!LFR>                             numb_acc_1,         & !intent()  
!LFR>                             imtr_method,        & !intent()  
!LFR>                             mw_so4,             & !intent()  
!LFR>                             avo,                & !intent()  
!LFR>                             pi,                 & !intent()  
!LFR>                             ndiag_aero            !intent()
!LFR>                
!LFR>       use Aero_Config, ONLY: &
!LFR>                             naerobox,           & !intent(in)
!LFR>                             nmodes,             & !intent()
!LFR>                             nweights,           & !intent()
!LFR>                             nmodes,             & !intent()
!LFR>                             nweights              !intent()

      implicit none

      contains

!>
!!----------------------------------------------------------------------------------------------------------------------
!!     this routine rescales all aerosol and gas-phase species to enforce
!!     mass conservation to machine precision.
!!----------------------------------------------------------------------------------------------------------------------
      subroutine MassAdj(aero,gas,spcmass1,spcmass2,emis_mass,aqso4rate,tstep)
      use Aero_Setup, only: &
                           sulf_map, &  !intent()
                           bcar_map, &  !intent()
                           ocar_map, &  !intent()
                           dust_map, &  !intent()
                           seas_map     !intent()
      IMPLICIT NONE

      ! Arguments.
 
      real(8), intent(inout) :: aero(naerobox)         !< aerosol conc. [ug/m^3] or [#/m^3]
      real(8), intent(inout) :: gas(ngases)            !< gas-phase conc. [ug/m^3]
      real(8), intent(in)    :: spcmass1(nmass_spcs+2) !< initial total mass spc. conc. [ug/m^3]
      real(8), intent(inout) :: spcmass2(nmass_spcs+2) !< final   total mass spc. conc. [ug/m^3]
      real(8), intent(in)    :: emis_mass(nemis_spcs)  !< mass emission rates [ug/m^3/s]
      real(8), intent(in)    :: aqso4rate              !< in-cloud so4 production rate [ug/m^3/s]
      real(8), intent(in)    :: tstep                  !< model physics time step [s]

      ! local variables.

      integer       :: i
      real(8)       :: scale(nmass_spcs+2)             ! scale factor for mass adjustment
      real(8), save :: scalemax = 1.0d-80
      real(8), save :: scalemin = 1.0d+80
      
      !----------------------------------------------------------------------------------------------------------------
      ! get the precise mass conc. that should exist at the end of the time
      ! step, divided by the actual mass conc. at the end of the time step.
      !----------------------------------------------------------------------------------------------------------------
      spcmass2(:) = spcmass2(:) + tinydenom 
      scale(1) = ( spcmass1(1) + ( aqso4rate + emis_mass(1) + emis_mass(2)  ) * tstep ) / spcmass2(1) 
      scale(2) = ( spcmass1(2) + (             emis_mass(3) + emis_mass(8)  ) * tstep ) / spcmass2(2) 
      scale(3) = ( spcmass1(3) + (             emis_mass(4) + emis_mass(9)  ) * tstep ) / spcmass2(3) 
      scale(4) = ( spcmass1(4) + (             emis_mass(5) + emis_mass(10) ) * tstep ) / spcmass2(4) 
      scale(5) = ( spcmass1(5) + (             emis_mass(6) + emis_mass(7)  ) * tstep ) / spcmass2(5) 
      scale(6) = ( spcmass1(6)                                                        ) / spcmass2(6) 
      scale(7) = ( spcmass1(7)                                                        ) / spcmass2(7) 

      ! write(*,'(7f14.9)') scale(:)
      ! write(*,'(7e14.6)') spcmass1(6), spcmass2(6), spcmass1(7), spcmass2(7)
      !----------------------------------------------------------------------------------------------------------------

      aero( sulf_map(:) ) = aero( sulf_map(:) ) * scale(1)
      aero( bcar_map(:) ) = aero( bcar_map(:) ) * scale(2)
      aero( ocar_map(:) ) = aero( ocar_map(:) ) * scale(3)
      aero( dust_map(:) ) = aero( dust_map(:) ) * scale(4)
      aero( seas_map(:) ) = aero( seas_map(:) ) * scale(5)
      aero( mass_no3    ) = aero( mass_no3    ) * scale(6)
      aero( mass_nh4    ) = aero( mass_nh4    ) * scale(7)
      gas ( gas_h2so4   ) = gas ( gas_h2so4   ) * scale(1)
      gas ( gas_hno3    ) = gas ( gas_hno3    ) * scale(6)
      gas ( gas_nh3     ) = gas ( gas_nh3     ) * scale(7)
       
      
!----------------------------------------------------------------------------------------------------------------------
!     track the maximum and minimum scale factors required.
!----------------------------------------------------------------------------------------------------------------------
      if( write_log ) then  
        write(31,90000) spcmass1(:)
        write(31,90000) spcmass2(:)
        write(31,90000) scale(:)
        write(31,*) '  '
        write(32,90000) scale(:)
        do i=1, nmass_spcs+2
          if(     scale(i) .gt. scalemax ) then
            scalemax = scale(i)
            write(33,90001) scale(i), scalemax, scalemin, i, spcmass1(i), spcmass2(i)
          elseif( scale(i) .lt. scalemin ) then
            scalemin = scale(i)
            write(33,90001) scale(i), scalemax, scalemin, i, spcmass1(i), spcmass2(i)
          endif
        enddo
      endif

90000 format(7d15.6)
90001 format(3d15.6,i6,2d15.6)
      return
      end subroutine MassAdj


!>----------------------------------------------------------------------------------------------------------------------
!!     Routine to calculate the total mass concentration of each model species:
!!     SULF, BCAR, OCAR, DUST, SEAS, NO3, NH4. Aerosol water is not treated. 
!!----------------------------------------------------------------------------------------------------------------------
      SUBROUTINE SpcMasses(aero,gas,spcmass)
      use aero_setup, only: &
                           sulf_map, &  !intent()
                           bcar_map, &  !intent()
                           ocar_map, &  !intent()
                           dust_map, &  !intent()
                           seas_map     !intent()
      implicit none
      real(8) :: aero(naerobox)
      real(8) :: gas(ngases)    
      real(8) :: spcmass(nmass_spcs+2)
      spcmass(1) = sum( aero( sulf_map(:) ) ) + gas( gas_h2so4 )
      spcmass(2) = sum( aero( bcar_map(:) ) )
      spcmass(3) = sum( aero( ocar_map(:) ) )
      spcmass(4) = sum( aero( dust_map(:) ) )
      spcmass(5) = sum( aero( seas_map(:) ) )
      spcmass(6) = aero( mass_no3 )           + gas( gas_hno3 )
      spcmass(7) = aero( mass_nh4 )           + gas( gas_nh3  )
      return
      end subroutine SpcMasses


!>---------------------------------------------------------------------------------------------------------------------
!! DLW, 102306: derived from function GETAF of CMAQ v4.4.
!!
!! GETXNUM = ln( Dij / Dgi ) / ( sqrt(2) * ln(Sgi) ), where
!!
!!      Dij is the diameter of intersection,
!!      Dgi is the median diameter of the smaller size mode, and
!!      Sgi is the geometric standard deviation of smaller mode.
!!
!! A quadratic equation is solved to obtain GETXNUM, following the method of Press et al. 1992.
!!  
!! REFERENCES:
!!
!!  1. Binkowski, F.S. and S.J. Roselle, Models-3 Community Multiscale Air Quality (CMAQ) 
!!     model aerosol component 1: Model Description.  J. Geophys. Res., Vol 108, No D6, 4183
!!     doi:10.1029/2001JD001409, 2003.
!!  2. Press, W.H., S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery, Numerical Recipes in 
!!     Fortran 77 - 2nd Edition. Cambridge University Press, 1992.
!!----------------------------------------------------------------------------------------------------------------------
      REAL(8) FUNCTION GetXNum(ni,nj,dgni,dgnj,xlsgi,xlsgj)
      IMPLICIT NONE
      
      ! Arguments.
      
      real(8) :: ni         !< Aitken       mode number concentration [#/m^3] 
      real(8) :: nj         !< accumulation mode number concentration [#/m^3] 
      real(8) :: dgni       !< Aitken       mode geo. mean diameter [um] 
      real(8) :: dgnj       !< accumulation mode geo. mean diameter [um]
      real(8) :: xlsgi      !< Aitken       mode ln(geo. std. dev.) [1]
      real(8) :: xlsgj      !< accumulation mode ln(geo. std. dev.) [1]

      ! Local variables. 

      real(8) :: aa, bb, cc, disc, qq, alfa, l, yji
      real(8), parameter :: sqrt2 = 1.414213562d+00

      alfa = xlsgi / xlsgj
      yji = log( dgnj / dgni ) / ( sqrt2 * xlsgi )
      l = log( alfa * nj / ni)

      ! calculate quadratic equation coefficients & discriminant.
      
      aa = 1.0d+00 - alfa * alfa
      bb = 2.0d+00 * yji * alfa * alfa
      cc = l - yji * yji * alfa * alfa
      disc = bb * bb - 4.0d+00 * aa * cc

      ! if roots are imaginary, return a negative getaf value so that no imtr takes place.
      
      if( disc .lt. 0.0d+00 ) then
        getxnum = - 5.0d+00         ! error in intersection
        return
      endif
      
      ! equation 5.6.4 of press et al. 1992.
      
      qq = -0.5d+00 * ( bb + sign( 1.0d+00, bb ) * sqrt(disc) )

      ! return solution of the quadratic equation that corresponds to a
      ! diameter of intersection lying between the median diameters of the 2 modes.
      
      getxnum = cc / qq       ! see equation 5.6.5 of press et al.
      
      ! write(*,*)'getxnum = ', getxnum
      return
      end function GetXNum


      END MODULE Aero_Subs
 
