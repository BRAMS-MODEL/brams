!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237] 
!******************************COPYRIGHT**************************************
SUBROUTINE delta_temp(                                                        &
  n_olevs,f_ocean,kappa,lambda_l,lambda_o,mu,q,dtemp_l,dtemp_o                &
)

USE jules_print_mgr, ONLY:                                                    &
  jules_message,                                                              &
  jules_print
IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   This subroutine calls a simple two-box thermal model in order to
!   estimate the global mean land surface temperature.
!   The code uses an implicit solver.
!
! Code Owner: Please refer to ModuleLeaders.txt
!             This file belongs in IMOGEN
! Written by: C. Huntingford, 11/09/98.
!
! Code Description:
!   Language: Fortran 90.
!   
!-----------------------------------------------------------------------------

INTEGER ::                                                                    &
  n_olevs  !IN Number of ocean thermal layers

REAL ::                                                                       &
  q,                                                                          &
           !IN Increase in global radiative forcing (W/m2)
  lambda_l,                                                                   &
           !IN Inverse climate sensitivity over land (W/K/m2)
  lambda_o,                                                                   &
           !IN Inverse climate sensitivity over ocean (W/K/m2)
  mu,                                                                         &
           !IN Ratio of land-to-ocean temperature anomalies (.)
  f_ocean,                                                                    &
           !IN Fractional coverage of planet with ocean
  kappa    !IN Ocean eddy diffusivity (W/m/K)

INTEGER ::                                                                    &
  iter_per_year,                                                              &
           !WORK Iterations per year
  timesteps,                                                                  &
           !WORK Number of iterations per call to DELTA_TEMP
  i,j      !WORK Looping parameter
PARAMETER(iter_per_year = 20)

REAL ::                                                                       &
  dtemp_l,                                                                    &
           !OUT Land mean temperature anomaly
  dtemp_o(n_olevs)
           !IN/OUT Land mean temperature anomaly

REAL ::                                                                       &
  rhocp,                                                                      &
           !WORK Rho times cp for the ocean (J/K/m3)
  flux_top,                                                                   &
           !WORK Flux into the top of the ocean (W/m2)
  flux_bottom
           !WORK Flux into the top of the ocean (W/m2)

REAL ::                                                                       &
  dtime,                                                                      &
           !WORK Timestep (s)
  dz(1:n_olevs),                                                              &
           !WORK Distance step (m)
  sec_year !WORK Seconds in a year (s)

REAL ::                                                                       &
  r1(1:n_olevs),                                                              &
           !WORK Mesh ratio
  r2(1:n_olevs),                                                              &
           !WORK Mesh ratio
  lambda_old,                                                                 &
           !WORK Inhomogenous component in implicit scheme (/s)
  lambda_new,                                                                 &
           !WORK Inhomogenous component in implicit scheme (/s)
  p        !WORK ``Dirchlet'' part of mixed boundary
           !     condition at ocean surface (/m)

REAL ::                                                                       &
  u_old(1:n_olevs),                                                           &
           !WORK Old ocean temperature anomalies
  u_new(1:n_olevs)
           !WORK New ocean temperature anomalies

REAL ::                                                                       &
  factor_dep,                                                                 &
           !WORK Factor by which layers are changed
  depth,                                                                      &
           !WORK Cumulated depth of layers (m)
  dz_top   !WORK Depth of the top layer (m)

PARAMETER(rhocp = 4.04e6)
PARAMETER(sec_year = 3.1536e7)


!-----------------------------------------------------------------------
! Set up parameters for the model
!-----------------------------------------------------------------------

dtime = sec_year / REAL(iter_per_year)
timesteps = iter_per_year


! Here the variable depths are calculated, along with the "R" values

! The mesh is as follows
!
!      -------------            U(1)
!                  (This is the top layer, called TEM(0) at points)
!          R(1)
!
!      ------------             U(2)
!
!          R(2)
!
!
!           |
!           |
!          \ /
!
!      ------------             U(N_OLEVS) = 0.0

dz_top = 2.0
depth = 0.0

DO i = 1,n_olevs
  factor_dep = 2.71828**(depth / 500.0)
  dz(i) = dz_top * factor_dep
  depth = depth + dz(i)
END DO

! Now arrange that the lowest depth (here DEPTH = U(N_OLEVS+1)) is at 5000m
dz_top = dz_top * (5000.0 / depth)
DO i = 1,n_olevs
  dz(i) = dz(i) * (5000.0 / depth)
END DO

r1(1) = (kappa / rhocp) * (dtime / (dz_top * dz_top))
r2(1) = 0.0
DO i = 2,n_olevs
  r1(i) = (kappa / rhocp) * (dtime / (dz(i-1) * (dz(i) + dz(i-1))))
  r2(i) = (kappa / rhocp) * (dtime / (dz(i) * (dz(i) + dz(i-1))))
END DO

! Reset the new values of U_OLD for the start of this run.
DO i = 1,n_olevs
  u_old(i) = dtemp_o(i)
END DO

DO i = 1,timesteps

  lambda_old = -q / (kappa * f_ocean)
  lambda_new = -q / (kappa * f_ocean)

  p = ((1.0 - f_ocean) * lambda_l * mu) / (f_ocean * kappa)                   &
    + (lambda_o / kappa)

  CALL invert(                                                                &
    u_old,u_new,p,lambda_old,lambda_new,r1,r2,dz_top,n_olevs                  &
  )

  ! Now check that the flux out of the bottom is not too large.
  flux_top =  - 0.5 * (kappa * lambda_old * dtime)                            &
              - 0.5 * (kappa * lambda_new * dtime)                            &
              - 0.5 * (p * kappa * u_old(1) * dtime)                          &
              - 0.5 * (p * kappa * u_new(1) * dtime)
  flux_bottom = (dtime / (2.0 * dz(n_olevs))) * kappa                         &
              * (u_old(n_olevs) + u_new(n_olevs))

  IF (ABS(flux_bottom / (flux_top+0.0000001)) > 0.01) THEN
    WRITE(jules_message,*) 'Flux at bottom of the ocean is greater ' //       &
               'than 0.01 of top'
    CALL jules_print('delta_temp',jules_message)
  END IF


  ! Set calculated values at now ``old'' values
  DO j = 1,n_olevs
    u_old(j) = u_new(j)
  END DO

END DO

! At end of this model run, reset values of LAMBDA_L and LAMBDA_O

DO j = 1,n_olevs
  dtemp_o(j) = u_new(j)
END DO
dtemp_l = mu * dtemp_o(1)

RETURN

END SUBROUTINE delta_temp
