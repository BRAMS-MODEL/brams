! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!     Science routines for calculating river flow routing
!     using the RFM kinematic wave model
!     see Bell et al. 2007 Hydrol. Earth Sys. Sci. 11. 532-549
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

MODULE rivers_route_rfm_mod

CONTAINS

!###############################################################################
! subroutine rivers_route_rfm
!
!-----------------------------------------------------------------------------
! Description:
!   Perform the kinematic wave routing of surface and sub-surface runoff
!   Calculates river outflow (kg m-2 s-1) and baseflow (kg m-2 s-1) for the
!   RFM kinematic wave river routing model.
!
! Method:
!   See Bell et al. 2007 Hydrol. Earth Sys. Sci. 11. 532-549
!   This subroutine routes surface and subsurface runoff
!   for the RCM. The routing procedure is performed on the
!   RCM grid, so no regridding is needed. The routing
!   procedure is currently based on the Kinematic Wave
!   routing model used at CEH Wallingford.
!
! Author: V.A.Bell, CEH Wallingford, 21.08.03
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!
! Code History:
!    Modified from river2a.f by sjd 13/05/05
!    Modified from UM routine riv_rout-river2a.F90 by hl 13/04/14
!    Updated to route runoff on riv_pts vector only by hl 24/04/14
!
!  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 5.5:
! VERSION  DATE
!   6.0   12/09/03  Change DEF from A20 to A26. D. Robinson
!   6.0   12.09.03  Routing code added. V.A.Bell
!   x.x   02/05/12  Additional implementation and developments S. Dadson
!   x.x   06/01/15  Formal implementation within JULES code base H. Lewis
!
! NOTE ON UNITS:
!   This routine, based on Simon Dadson's work, includes a modification to the
!   'standard' RFM routines which assume a regular x/y-grid based
!   implementation to account for potential variable grid box areas (e.g. from
!   lat/lon grid).
!   Stores are calculated in units of m x m2 rather than mm, and flows are
!   initially calculated in units of m3/s.
!   For consistency with other routines (for now), the output is converted
!   again to a flux density kg/m2/s

SUBROUTINE rivers_route_rfm( sfc_runoff, sub_sfc_runoff,                      &
                              outflow, baseflow, riverout_rgrid )

USE jules_rivers_mod, ONLY:                                                   &
!  imported scalar parameters
       nstep_rivers, np_rivers,rivers_dlat, rivers_first, rivers_dx           &
       ,runoff_factor, a_thresh, cland, criver, cbland, cbriver               &
       ,retl, retr                                                            &
!  imported arrays with intent(in)
       ,rfm_flowobs1_rp                                                       &
!  imported arrays with intent(inout)
       ,rivers_next_rp, rfm_iarea_rp, rfm_land_rp                             &
       ,rfm_substore_rp, rfm_surfstore_rp, rfm_flowin_rp, rfm_bflowin_rp      &
       ,rfm_rivflow_rp, rfm_baseflow_rp, rivers_boxareas_rp

USE timestep_mod, ONLY: timestep

USE planet_constants_mod, ONLY: planet_radius  ! the Earth's radius (m)

USE conversions_mod, ONLY: pi_over_180

USE water_constants_mod, ONLY: rho_water

USE jules_print_mgr, ONLY:                                                    &
   jules_message,                                                             &
   jules_print

!-------------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

! IN Arguments

! Array arguments with intent(in)
REAL(KIND=real_jlslsm), INTENT(IN) :: sfc_runoff(np_rivers)
       !  average rate of surface runoff since last rivers call (kg m-2 s-1)
REAL(KIND=real_jlslsm), INTENT(IN) :: sub_sfc_runoff(np_rivers)
       !  average rate of sub-surface runoff since last call (kg m-2 s-1)

REAL(KIND=real_jlslsm), INTENT(OUT) :: outflow(np_rivers)
       !  rate of channel surface flow leaving gridbox (kg m-2 s-1)
REAL(KIND=real_jlslsm), INTENT(OUT) :: baseflow(np_rivers)
       !  rate of channel base flow leaving gridbox (kg m-2 s-1)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: riverout_rgrid(np_rivers)
       ! River outflow into the ocean on river grid (kg s-1)

! internal variables
INTEGER ::                                                                    &
     landtype                                                                 &
       !  local for land type
     ,landtype_next                                                           &
       !  local for land type of downstream point
     ,rn                                                                      &
       !  local co-ords of downstream point
     ,ip
       !  co-ordinate counters in do loops

REAL(KIND=real_jlslsm) ::                                                     &
   landtheta, rivertheta                                                      &
       !  surface wave speed factors
   ,sublandtheta, subrivertheta                                               &
       !  sub-surface wave speed factors
   ,returnflow                                                                &
       !  returnflow (m3 per timestep)
   ,flowobs1_m3s                                                              &
       !  initial river flow [m3/s]
   ,dt                                                                        &
       !  river routing model timestep (s)
   ,dx
       !  distance between midpoints of neighbouring cells (m)

REAL(KIND=real_jlslsm) ::                                                     &
   substore_n(np_rivers)                                                      &
       !   subsurface store at next timestep (m3 per timestep)
   ,surfstore_n(np_rivers)                                                    &
       !   surface store at next timestep (m3 per timestep)
   ,flowin_n(np_rivers)                                                       &
       !   surface lateral inflow next time (m3 per timestep)
   ,bflowin_n(np_rivers)                                                      &
       !   sub-surface lateral inflow next time (m3 per timestep)
   ,surf_roff(np_rivers)                                                      &
       !   INTERNAL surf_runoff (m3 per timestep)
   ,sub_surf_roff(np_rivers)
       !   INTERNAL sub_surf_runoff (m3 per timestep)

INTEGER :: nland, nsea, nriv

!-------------------------------------------------------------------------------
! Set up rivers parameters
!-------------------------------------------------------------------------------

! rivers model timestep(s)
dt = REAL(nstep_rivers) * timestep

! horizontal gridsize (m)
IF (rivers_dx <= 0) THEN
  dx = planet_radius * (ABS(rivers_dlat) * pi_over_180)
ELSE
  dx = rivers_dx
END IF

! Wave speed factors (dimensionless)
rivertheta    = criver  * dt / dx
landtheta     = cland   * dt / dx
sublandtheta  = cbland  * dt / dx
subrivertheta = cbriver * dt / dx

! Check condition for numerical stability
IF (landtheta > 1.0 .OR. sublandtheta > 1.0 .OR.                              &
    rivertheta > 1.0 .OR. subrivertheta > 1.0) THEN

  WRITE(jules_message,*)'WARNING: rivers_route_rfm: ' //                      &
                        'Finite difference method will be unstable in RFM,'   &
                        // ' setting thetas to zero'
  CALL jules_print('rivers_route_rfm',jules_message)

  rivertheta    = 0.0
  landtheta     = 0.0
  sublandtheta  = 0.0
  subrivertheta = 0.0
END IF

!-------------------------------------------------------------------------------
! Initialise variables at first timestep
!-------------------------------------------------------------------------------
nland = 0
nsea  = 0
nriv  = 0

IF (rivers_first) THEN
  ! From the cumulative catchment areas dataset (rfm_iarea), determine
  ! which grid cells are land or river
  !------------------------------------------------------------------
  ! rfm_land_rp() is set at the first entry, and the values must be
  ! retained for subsequent timesteps
  !-------------------------------------------------------------------
  rivers_first = .FALSE.

  DO ip = 1, np_rivers

    IF (rfm_iarea_rp(ip) < 0) THEN
      rfm_land_rp(ip) = 0     !Sea
      nsea = nsea + 1
    ELSE IF (rfm_iarea_rp(ip) > a_thresh) THEN
      rfm_land_rp(ip) = 1     !river
      nriv = nriv + 1
    ELSE
      rfm_land_rp(ip) = 2     !land
      nland = nland + 1
    END IF
    !      rarea(ip) = REAL(rfm_iarea_rp(ip)) + 1.0
          ! include current point, so add 1

          ! set to sea if top or left edge drains outside model domain
    IF (rivers_next_rp(ip) < 1 .OR. rivers_next_rp(ip) > np_rivers) THEN
      rfm_land_rp(ip) = 0  !Sea
    END IF

    ! Initialise surface and sub-surface stores using flow observations if available
    IF ( rfm_flowobs1_rp(ip) > 0.0 ) THEN
      flowobs1_m3s = rfm_flowobs1_rp(ip) * rivers_boxareas_rp(ip) / rho_water
      rfm_surfstore_rp(ip) = flowobs1_m3s * dt / rivertheta
      rfm_substore_rp(ip) = flowobs1_m3s * dt / subrivertheta
    END IF
  END DO
END IF   ! end riversFirst

!-------------------------------------------------------------------------------
! Processing for each timestep
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Convert runoff from (kg m-2 s-1) to (m3 per gridcell per timestep)
!-------------------------------------------------------------------------------
DO ip = 1,np_rivers
  IF ( sfc_runoff(ip) >= 0.0 .AND. sub_sfc_runoff(ip) >= 0.0 ) THEN

    surf_roff(ip) = runoff_factor * sfc_runoff(ip) *                          &
                          dt * rivers_boxareas_rp(ip) / rho_water
    sub_surf_roff(ip) = runoff_factor * sub_sfc_runoff(ip) *                  &
                          dt * rivers_boxareas_rp(ip) / rho_water
  ELSE
    ! ignore no data (-1.00e20) and other negative values
    surf_roff(ip) = 0.0
    sub_surf_roff(ip) = 0.0
  END IF

  !-------------------------------------------------------------------------------
  ! Initialise accumulated inflows and stores for the next timestep
  !-------------------------------------------------------------------------------
  flowin_n(ip)    = 0.0
  bflowin_n(ip)   = 0.0
  surfstore_n(ip) = 0.0
  substore_n(ip)  = 0.0
END DO

!-------------------------------------------------------------------------------
! Rivers runoff using simple kinematic wave model (see Lewis et al. 2018:Appx.B)
!-------------------------------------------------------------------------------
DO ip = 1,np_rivers

  rn = rivers_next_rp(ip)
  landtype = NINT( rfm_land_rp(ip) )

  IF (landtype == 2) THEN  !land

    ! land surface
    surfstore_n(ip) = (1.0 - landtheta) * rfm_surfstore_rp(ip) +              &
                             rfm_flowin_rp(ip) + surf_roff(ip)

    ! land subsurface
    substore_n(ip) = (1.0 - sublandtheta) * rfm_substore_rp(ip) +             &
                             rfm_bflowin_rp(ip) + sub_surf_roff(ip)

    ! return flow
    IF (retl > 0) THEN ! changed to allow for -ve retr/l (sjd: 4/6/09)
      returnflow = MAX( ABS( substore_n(ip) * retl ), 0.0 )
    ELSE
      returnflow = -1.0 * MAX( ABS( surfstore_n(ip) * retl ), 0.0 )
    END IF

    substore_n(ip)  = substore_n(ip)  - returnflow
    surfstore_n(ip) = surfstore_n(ip) + returnflow

    IF ( rn > 0 .AND. rn /= ip ) THEN
      flowin_n(rn)  = flowin_n(rn)  + landtheta    * rfm_surfstore_rp(ip)
      bflowin_n(rn) = bflowin_n(rn) + sublandtheta * rfm_substore_rp(ip)
    END IF

    rfm_rivflow_rp(ip)  = rfm_surfstore_rp(ip) * (landtheta    / dt)
    rfm_baseflow_rp(ip) = rfm_substore_rp(ip)  * (sublandtheta / dt)

  ELSE IF (landtype == 1) THEN  !river

    ! river subsurface
    substore_n(ip) = (1.0 - subrivertheta) * rfm_substore_rp(ip) +            &
                          rfm_bflowin_rp(ip) + sub_surf_roff(ip)

    ! river surface
    surfstore_n(ip) = (1.0 - rivertheta) * rfm_surfstore_rp(ip) +             &
                            rfm_flowin_rp(ip) + surf_roff(ip)

    ! return flow
    IF (retr > 0) THEN ! changed to allow for -ve retr/l (sjd: 4/6/09)
      returnflow = MAX( ABS( substore_n(ip) * retr ), 0.0 )
    ELSE
      returnflow = -1.0 * MAX( ABS( surfstore_n(ip) * retr ), 0.0 )
    END IF
    substore_n(ip)  = substore_n(ip)  - returnflow
    surfstore_n(ip) = surfstore_n(ip) + returnflow

    IF ( rn > 0 .AND. rn /= ip ) THEN
      flowin_n(rn)  = flowin_n(rn)  + rivertheta    * rfm_surfstore_rp(ip)
      bflowin_n(rn) = bflowin_n(rn) + subrivertheta * rfm_substore_rp(ip)
    END IF

    rfm_rivflow_rp(ip)  = rfm_surfstore_rp(ip) * (rivertheta    / dt)
    rfm_baseflow_rp(ip) = rfm_substore_rp(ip)  * (subrivertheta / dt)

  END IF ! land or river

  ! accumulate flow into the sea if the next point is sea
  ! (assume it's equal to river flow in the adjacent land pt)
  IF ( rn > 0 ) THEN
    landtype_next = NINT( rfm_land_rp(rn) )
    IF ( landtype_next == 0) THEN
      rfm_rivflow_rp(rn) = rfm_rivflow_rp(rn) + rfm_rivflow_rp(ip)
      ! Calculate flow out to sea
      riverout_rgrid(rn) = riverout_rgrid(rn) +                               &
                         ( sfc_runoff(ip) + sub_sfc_runoff(ip) ) *            &
                           rivers_boxareas_rp(ip)
    END IF
  END IF
END DO !end of rivers loop, ip

!-------------------------------------------------------------------------------
! Housekeeping for next timestep
!-------------------------------------------------------------------------------
DO ip = 1,np_rivers
  ! keep inflows for next timestep
  rfm_flowin_rp(ip)  = flowin_n(ip)
  rfm_bflowin_rp(ip) = bflowin_n(ip)

  ! keep rivers stores for next timestep (m3)
  rfm_surfstore_rp(ip) = surfstore_n(ip)
  rfm_substore_rp(ip)  = substore_n(ip)

  !-------------------------------------------------------------------------------
  ! Return flows in flux density units kg/m2/s
  !-------------------------------------------------------------------------------
  ! HL: N.B. Add option to select what output units required for river flow
  outflow(ip)  = rfm_rivflow_rp(ip)  * rho_water / rivers_boxareas_rp(ip)
  baseflow(ip) = rfm_baseflow_rp(ip) * rho_water / rivers_boxareas_rp(ip)
END DO

END SUBROUTINE rivers_route_rfm

!###############################################################################
!###############################################################################

END MODULE rivers_route_rfm_mod
