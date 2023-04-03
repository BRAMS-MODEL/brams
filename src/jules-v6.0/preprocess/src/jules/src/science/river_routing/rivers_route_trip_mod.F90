! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!     Science routines for calculating river flow routing
!     using the TRIP model
!     see Oki et al 1999 J.Met.Soc.Japan, 77, 235-255.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

MODULE rivers_route_trip_mod

CONTAINS

!###############################################################################
! subroutine rivers_route_trip
!
!-----------------------------------------------------------------------------
! Description:
!   Calculates river outflow (kg s-1) and updates channel storage for the
!   TRIP river routing model.
!
! Method:
!   See Oki et al. 1999, J.Met.Soc.Japan, 77, 235-255.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

SUBROUTINE rivers_route_trip( sfc_runoff, sub_sfc_runoff, outflow, baseflow,  &
                              riverout_rgrid )

USE jules_rivers_mod, ONLY:                                                   &
!  imported scalars with intent(in)
     np_rivers,nstep_rivers,nseqmax,rivers_meander,rivers_speed               &
!  imported arrays with intent(in)
    ,rivers_next_rp,rivers_seq_rp,rivers_sto_rp                               &
    ,rivers_boxareas_rp, rivers_lat_rp, rivers_lon_rp, rivers_dir_rp

USE rivers_utils, ONLY:                                                       &
!  imported procedures
     get_rivers_len_rp

USE timestep_mod, ONLY:                                                       &
   timestep

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-------------------------------------------------------------------------------

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
!  arrays with intent(in)
     sfc_runoff(np_rivers)                                                    &
     !  average rate of surface runoff since last rivers call (kg m-2 s-1)
     ,sub_sfc_runoff(np_rivers)
     !  average rate of sub-surface runoff since last call (kg m-2 s-1)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
!  arrays with intent(out)
     outflow(np_rivers)                                                       &
     !  rate of channel surface flow leaving gridbox (kg m-2 s-1)
     ,baseflow(np_rivers)                                                     &
!  rate of channel base flow leaving gridbox (kg m-2 s-1)
    ,riverout_rgrid(np_rivers)
! River outflow into the ocean on river grid (kg s-1)

INTEGER ::                                                                    &
!  local scalars (work/loop counters)
     ip,iseq

REAL(KIND=real_jlslsm) ::                                                     &
!  local scalars
     coeff                                                                    &
     !  coefficient in the routing model (s-1)
     ,exp_coeffdt                                                             &
     !  working variable exp[c*dt]
     ,dt                                                                      &
     !  timestep of routing model (s)
     ,store_old
     !  channel storage (kg)

REAL(KIND=real_jlslsm) ::                                                     &
!  local arrays
    inflow(np_rivers)                                                         &
    !  rate of channel flow entering gridbox (kg s-1)
    ,riverslength(np_rivers)
    !  distance between gridpoints (m)

!-------------------------------------------------------------------------------

dt = REAL(nstep_rivers) * timestep

! Initialise inflow with total runoff generated over each gridbox.
! Convert runoff to ks s-1 flux
DO ip = 1, np_rivers
  inflow(ip) = ( sfc_runoff(ip) + sub_sfc_runoff(ip) ) *                      &
                                  rivers_boxareas_rp(ip)

  ! Initialise outflow (for clarity at non-rivers points).
  outflow(ip) = 0.0
  baseflow(ip) = 0.0
  riverout_rgrid(ip) = 0.0
END DO

!-------------------------------------------------------------------------------
! Calculate distance between grid points
!-------------------------------------------------------------------------------
CALL get_rivers_len_rp( np_rivers, rivers_next_rp,                            &
                        rivers_lat_rp, rivers_lon_rp, riverslength )

!-------------------------------------------------------------------------------
! Loop over rivers points with valid flow direction
!-------------------------------------------------------------------------------

DO iseq = 1, nseqmax

  DO ip = 1,np_rivers

    !   Get index (location in rivers vector) of the point to consider.
    IF (NINT(rivers_seq_rp(ip)) == iseq) THEN

      !-------------------------------------------------------------------------------
      !   Calculate the coefficient "c" of the model.
      !   c=u/(d*r), where u is effective flow speed,
      !   d is distance between gridpoints, and r is meander ratio.
      !-------------------------------------------------------------------------------
      coeff = rivers_speed / ( riverslength(ip) * rivers_meander )
      exp_coeffdt = EXP(-(coeff * dt))

      !-------------------------------------------------------------------------------
      !   Save value of channel storage at start of timestep.
      !-------------------------------------------------------------------------------
      store_old = rivers_sto_rp(ip)

      !-------------------------------------------------------------------------------
      !   Calculate channel storage at end of timestep.
      !   Eqn.4 of Oki et al, 1999, J.Met.Soc.Japan, 77, 235-255.
      !-------------------------------------------------------------------------------
      rivers_sto_rp(ip) = store_old * exp_coeffdt                             &
                          + ( 1.0 - exp_coeffdt ) * inflow(ip) / coeff

      !-------------------------------------------------------------------------------
      !   Calculate outflow as inflow minus change in storage.
      !-------------------------------------------------------------------------------
      outflow(ip) = inflow(ip) + (store_old - rivers_sto_rp(ip)) / dt

      !-------------------------------------------------------------------------------
      !   Add outflow to inflow of next downstream point.
      !-------------------------------------------------------------------------------
      IF ( rivers_next_rp(ip) > 0 ) THEN
        !     Get location in grid of next downstream point.
        inflow(rivers_next_rp(ip)) = inflow(rivers_next_rp(ip)) + outflow(ip)
      END IF

    END IF

  END DO !points
END DO !river sequences

! End main TRIP algorithm routine

!-------------------------------------------------------------------------------
! Add runoff on TRIP sea points to the river flow variable.
! Regridding may have resulted in some land runoff appearing in TRIP sea
! gridboxes, so we need to account for this water. Even without regridding,
! a land point that is sea in TRIP will not have been added to flow.
!-------------------------------------------------------------------------------
!  WHERE( .NOT. rivers_mask(:,:) ) outflow(:,:) = inflow(:,:)

! Catch all outflow going into the sea and send it to riverout_rgrid
DO ip = 1,np_rivers
  IF ( rivers_dir_rp(ip) == 9 ) THEN
    riverout_rgrid(ip) = outflow(ip)
  END IF

  ! Return flows in flux density units kg/m2/s
  outflow(ip) = outflow(ip) / rivers_boxareas_rp(ip)
END DO

END SUBROUTINE rivers_route_trip

!###############################################################################

END MODULE rivers_route_trip_mod

