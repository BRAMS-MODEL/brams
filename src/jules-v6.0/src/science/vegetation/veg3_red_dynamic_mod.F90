
! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in Veg3 Ecosystem Demography
! *****************************COPYRIGHT****************************************

MODULE veg3_red_dynamic_mod

IMPLICIT NONE

PRIVATE
PUBLIC :: veg3_red_dynamic

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='VEG3_RED_DYNAMIC_MOD'

CONTAINS

!-----------------------------------------------------------------------------
SUBROUTINE veg3_red_dynamic(                                                  &
                !IN Control vars
                dt,veg_index_pts,veg_index,veg3_ctrl,land_pts,                &
                nnpft,nmasst,                                                 &
                !IN red_parms
                red_parms,                                                    &
                !IN fields
                growth,mort_add,                                              &
                !IN state 
                veg_state,red_state,                                          &
                ! OUT Fields
                demographic_lit                                               &
                !OUT Diagnostics
                )

!Only get the data structures - the data comes through the calling tree
USE veg3_parm_mod,ONLY:  red_parm_type, veg3_ctrl_type
USE veg3_field_mod,ONLY:  veg_state_type, red_state_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Objects with INTENT IN
!-----------------------------------------------------------------------------
TYPE(red_state_type)  :: red_state
TYPE(red_parm_type)   :: red_parms
TYPE(veg_state_type)  :: veg_state
TYPE(veg3_ctrl_type)  :: veg3_ctrl

!----------------------------------------------------------------------------
! Integers with INTENT IN
!----------------------------------------------------------------------------
INTEGER, INTENT(IN) :: land_pts,nnpft,veg_index(land_pts),veg_index_pts,nmasst    

!----------------------------------------------------------------------------
! Reals with INTENT IN
!----------------------------------------------------------------------------
REAL, INTENT(IN)   ::                                                         &
growth(land_pts,nnpft),                                                       &
mort_add(land_pts,nnpft,nmasst),                                              &
dt

!-----------------------------------------------------------------------------
! Reals with INTENT OUT
!-----------------------------------------------------------------------------
REAL, INTENT(OUT)      ::                                                     &
demographic_lit(land_pts,nnpft)

!-----------------------------------------------------------------------------
!Local Vars
!-----------------------------------------------------------------------------
INTEGER                ::l,n,k

REAL                   ::                                                     &
G_tot(land_pts,nnpft),                                                        &
!  The total PFT carbon assimilate for the gridbox. (kgC/m2/year)
G_seed(land_pts,nnpft),                                                       &
              !  Total gridbox carbon assimilate devoted to recruitment. (kgC/m2/year)
G_str(land_pts,nnpft),                                                        &
              !  Total gridbox carbon assimilate devoted to vegetation structural growth. (kgC/m2/year)
g0(land_pts,nnpft),                                                           &
              !  Boundary growth for an individual member of the smallest mass cohort. (kgC/year)
g_mass(land_pts,nnpft,nmasst),                                                &
              !  Individual growth across the mass cohorts. (kgC/year)
plantNumDensity_g_sum(land_pts,nnpft),                                        &
              !  Summation of the relative cohort contribution towards the total PFT assimilate (/m2)
frac_shade(land_pts,3),                                                       &
              !  Compeitive shading of seedlings in each Plant Functional Group - Tree's, shrubs and grasses. (-)
dplantNumDensity_dt(land_pts,nnpft,nmasst),                                   &
              !  Net rate of change of population density within each mass cohort. (m2/year) 
flux_in(land_pts,nnpft,nmasst),                                               &
              ! Rate of change of population growing into a mass cohort. (m2/year) 
flux_out(land_pts,nnpft,nmasst),                                              &
              ! Rate of change of population growing out of a mass cohort. (m2/year)
frac_check(land_pts,nnpft)    
             !  The difference between the minimum vegetation fraction and the updated fraction. (-)

!End of headers

! Initialise vars
demographic_lit(:,:) = 0.0
!G_tot(:,:) = 0.
G_seed(:,:) = 0.0
G_str(:,:) = 0.0
g0(:,:) = 0.0
g_mass(:,:,:) = 0.0
plantNumDensity_g_sum(:,:) = 0.0
frac_shade(:,:) = 0.0
frac_check(:,:) = 0.0
dplantNumDensity_dt(:,:,:) = 0.0
flux_in(:,:,:) = 0.0
flux_out(:,:,:) = 0.0

! Loop finds the growth boundary condition, the fraction of shade for seedlings
DO k = 1,nmasst
  DO n = 1,nnpft
    DO l = 1,land_pts
      IF (k < red_parms%mclass(n)) THEN
        ! Sum product of the number density and the allometric scaling 
        plantNumDensity_g_sum(l,n) = plantNumDensity_g_sum(l,n)               &
          + red_state%plantNumDensity(l,n,k) * red_state%g_mass_scale(n,k)
      ELSE IF (k == red_parms%mclass(n)) THEN
        plantNumDensity_g_sum(l,n) = plantNumDensity_g_sum(l,n)               &
          + red_state%plantNumDensity(l,n,k) * red_state%g_mass_scale(n,k)
        ! Partition the growth

        IF (growth(l,n) >= 0.0) THEN
          G_tot(l,n) = veg_state%frac(l,n) * growth(l,n)
          G_seed(l,n) = red_parms%alpha(n) * G_tot(l,n)
          G_str(l,n) = (1.0 - red_parms%alpha(n)) * G_tot(l,n)
        ELSE
          G_tot(l,n) = 0.0
          G_seed(l,n) = 0.0
          G_str(l,n) = 0.0         
          demographic_lit(l,n) = demographic_lit(l,n)                         &
            - veg_state%frac(l,n) * growth(l,n)
        END IF
        
        IF (plantNumDensity_g_sum(l,n) > 0) THEN
          g0(l,n) = G_str(l,n) / plantNumDensity_g_sum(l,n)
        END IF
        
        ! Add the shaded regions per PFT hierarchy
        frac_shade(l,red_parms%veg_type(n)) =                                 &
          frac_shade(l,red_parms%veg_type(n))                                 &
          + veg_state%frac(l,n)
      END IF
    END DO
  END DO
END DO 

! Loop updates the demographic distribution.
DO k = 1,nmasst
  DO n = 1,nnpft
    DO l = 1,land_pts
      g_mass(l,n,k) = g0(l,n) * red_state%g_mass_scale(n,k)
      red_state%mort(l,n,k) = red_parms%mort_base(n) + mort_add(l,n,k)
      IF (k == 1) THEN
        ! Seedling flux
        flux_in(l,n,k) = G_seed(l,n) * (1.0 - frac_shade(l,red_parms%veg_type(n)))
        demographic_lit(l,n) = demographic_lit(l,n) +                         &
          frac_shade(l,red_parms%veg_type(n)) * G_seed(l,n)
      ELSE
        ! Flux into mass class 
        flux_in(l,n,k) = flux_out(l,n,k-1)
      END IF
      IF (k == red_parms%mclass(n)) THEN
        ! Truncate growth at the top mass class
        flux_out(l,n,k) = 0.0
        demographic_lit(l,n) = demographic_lit(l,n)                           &
          + red_state%plantNumDensity(l,n,k) * g_mass(l,n,k)
      ELSE
        flux_out(l,n,k) = red_state%plantNumDensity(l,n,k)                    &
        * g_mass(l,n,k) / (red_state%mass_mass(n,k+1) - red_state%mass_mass(n,k))
      END IF
        ! Update the number density

      IF (k <= red_parms%mclass(n)) THEN
        dplantNumDensity_dt(l,n,k) = flux_in(l,n,k) - flux_out(l,n,k)         &
          - red_state%mort(l,n,k) * red_state%plantNumDensity(l,n,k)
        demographic_lit(l,n) = demographic_lit(l,n) + red_state%mort(l,n,k)   &
          * red_state%plantNumDensity(l,n,k) * red_state%mass_mass(n,k)
        IF (ABS(dplantNumDensity_dt(l,n,k) * dt) >=                           &
          red_state%plantNumDensity(l,n,k)) THEN
          red_state%plantNumDensity(l,n,k) = 0.0
        ELSE
          red_state%plantNumDensity(l,n,k) = red_state%plantNumDensity(l,n,k) &
            + dplantNumDensity_dt(l,n,k) * dt
        END IF
        frac_check(l,n) = frac_check(l,n) + red_state%plantNumDensity(l,n,k)  &
          * red_state%crwn_area_mass(n,k)
        ! If the resultant vegetation fraction is less than the minimum       
        ! fraction, add trees to the lowest mass class to make up the difference
        IF (k == red_parms%mclass(n) .AND. frac_check(l,n)                    &
          < red_parms%frac_min(n)) THEN
          red_state%plantNumDensity(l,n,1) = (red_parms%frac_min(n)           &
            - frac_check(l,n)) / red_state%crwn_area_mass(n,1) 
          frac_check(l,n) = red_parms%frac_min(n)
        END IF
      END IF
    END DO
  END DO
END DO

! for now - sumproduct across the mass classes to output back to JULES
DO k = 1,nmasst
  DO n = 1,nnpft
    DO l = 1,land_pts
      IF (k == 1) THEN
        veg_state%frac(l,n) = 0.0
        veg_state%vegCpft(l,n) = 0.0
      END IF
      veg_state%frac(l,n) = veg_state%frac(l,n)                               &
        + red_state%plantNumDensity(l,n,k) * red_state%crwn_area_mass(n,k)
      veg_state%vegCpft(l,n) = veg_state%vegCpft(l,n)                         &
        + red_state%plantNumDensity(l,n,k) * red_state%mass_mass(n,k)
      IF (k == red_parms%mclass(n)) THEN
        veg_state%vegCpft(l,n) = veg_state%vegCpft(l,n) / veg_state%frac(l,n)
      END IF
    END DO
  END DO
END DO

END SUBROUTINE veg3_red_dynamic

!-----------------------------------------------------------------------------

END MODULE veg3_red_dynamic_mod
