!###############################################################################
!###############################################################################

MODULE rivers_utils

USE conversions_mod, ONLY: pi_over_180

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RIVERS_UTILS'

! The following constants are defined in the JULES standalone version of
! planet_constants_mod, but not in the corresponding UM parameters file (yet!)

! planet equatorial radius 'a' (m)
! from Appendix A of Oki and Sud, 1998, Earth Interactions, Vol.2, Paper 1.
REAL(KIND=real_jlslsm), PARAMETER :: planet_eq_radius    = 6378136.0

! eccentricity of planet spheroid
!         (spheroid - lines of constant latitude are circular,
!         meridional cross-section is an ellipse). From Appendix A of
!                 Oki and Sud, 1998, Earth Interactions, Vol.2, Paper 1.
REAL(KIND=real_jlslsm), PARAMETER :: eccen               = 0.08181974

REAL(KIND=real_jlslsm), PARAMETER :: eccensq             = 0.00669447
                                                         ! eccen*eccen

!-----------------------------------------------------------------------------
! Description:
!   Contains river routing utility functions for standalone running
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

CONTAINS

!###############################################################################
! subroutine get_rivers_len
!     Driver routine that calls procedures that calculates distance
!     between gridpoints.

SUBROUTINE get_rivers_len_rp( np_rivers, rivers_next, lat, lon, length )

USE jules_rivers_mod, ONLY:                                                   &
!  imported scalars with intent(in)
     rivers_dlon

IMPLICIT NONE

! Array arguments with intent(in)

INTEGER, INTENT(IN) ::  np_rivers,                                            &
                           ! size of river routing vector
                        rivers_next(np_rivers)
                          !  location of the next downstream point

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
                    lat(np_rivers),                                           &
                          !  latitude of points on rivers grid (degrees)
                    lon(np_rivers)
                          !  longitude of points on rivers grid (degrees)

! Array arguments with intent(out)

REAL(KIND=real_jlslsm), INTENT(OUT) :: length(np_rivers)
                          !  distance between each gridpoint and downstream
                          !  gridpoint (m)

! Local scalar variables.
INTEGER ::  ip,jp         !  work

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GET_RIVERS_LEN_RP'
!------------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Only calculate length at rivers grid points with valid flow direction
DO ip = 1,np_rivers
  !   Get coords of this point and any point immediately downstream
  jp = rivers_next(ip)

  IF ( jp > 0 ) THEN
    !     There is an outflow direction, and the downstream point is on the grid
    length(ip) = givelength( lat(ip), lon(ip), lat(jp), lon(jp) )

  ELSE
    !     Outflow to sea, or no outflow (pit or depression).
    !     This is relatively common, hence treatment separately from flow across
    !     grid edge. Use distance to an adjacent point at the same latitude
    length(ip) = givelength( lat(ip), lon(ip), lat(ip), lon(ip) + rivers_dlon )
  END IF
END DO ! np_rivers

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE get_rivers_len_rp
!###############################################################################

!###############################################################################
! subroutine get_rivers_len
!     Driver routine that calls procedures that calculates distance
!     between gridpoints.

SUBROUTINE get_rivers_len( x, xnext, y, ynext, lat, lon, length )

USE jules_rivers_mod, ONLY:                                                   &
!  imported scalars with intent(in)
     np_rivers,nx_rivers,ny_rivers,rivers_dlat,rivers_dlon                    &
!  imported arrays with intent(in)
    ,flow_dir_delta

IMPLICIT NONE

! Array arguments with intent(in)

INTEGER, INTENT(IN) ::                                                        &
                       x(np_rivers),                                          &
                          !  x index on the rivers grid
                       xnext(np_rivers),                                      &
                          !  x location of the next downstream point
                          !      <0 indicates next point is off grid
                          !      and gives the direction of the next point,
                          !      in terms of index in flow_dir_delta
                       y(np_rivers),                                          &
                          !  y index on the rivers grid
                       ynext(np_rivers)
                          !  y location of the next downstream point
                          !      <0 interpreted as for xnext.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
                    lat(nx_rivers,ny_rivers),                                 &
                          !  latitude of points on rivers grid (degrees)
                    lon(nx_rivers,ny_rivers)
                          !  longitude of points on rivers grid (degrees)

! Array arguments with intent(out)

REAL(KIND=real_jlslsm), INTENT(OUT) :: length(nx_rivers,ny_rivers)
                          !  distance between each gridpoint and downstream
                          !  gridpoint (m)

! Local scalar variables.
INTEGER ::  ip,ix,iy,jx,jy          !  work
REAL(KIND=real_jlslsm) ::  dx,dy    !  work

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GET_RIVERS_LEN'
!------------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Only calculate length at rivers grid points with valid flow direction
DO ip = 1,np_rivers
  !   Get coords of this point and any point immediately downstream
  ix = x(ip)
  iy = y(ip)
  jx = xnext(ip)
  jy = ynext(ip)

  IF ( jx > 0 ) THEN
    !     There is an outflow direction, and the downstream point is on the grid
    length(ix,iy) = givelength( lat(ix,iy), lon(ix,iy),                       &
                                lat(jx,jy), lon(jx,jy) )
  ELSE IF ( jx == -9 ) THEN
    !       Outflow to sea, or no outflow (pit or depression).
    !       This is relatively common, hence treatment separately from flow across
    !       grid edge. Use distance to an adjacent point at same latitude.

    length(ix,iy) = givelength( lat(ix,iy), lon(ix,iy),                       &
                                lat(ix,iy), lon(ix,iy) + rivers_dlon )
  ELSE
    !       Outflow across the edge of the grid. We would still like to calculate
    !       a length, so that, as far as possible, a regional implementation gives
    !       the same results as a global implementation, even for the edge points.
    !       NB This works for outflow from grid, but if a regional run is missing
    !       an inflow from outside the grid, answers will generally be different
    !       in global and regional implementations.
    !       The lat/lon of the downstream point is not available in these cases,
    !       so we ASSUME a "regular" lat/lon grid.
    !              jx is -1 * index in flowDirSet.

    !       Get number of gridboxes in each direction to the downstream location
    dx = REAL( flow_dir_delta(ABS(jx),1) )
    dy = REAL( flow_dir_delta(ABS(jx),2) )

    length(ix,iy) = givelength( lat(ix,iy), lon(ix,iy),                       &
                                lat(ix,iy) + dy * rivers_dlat,                &
                                lon(ix,iy) + dx * rivers_dlon )
  END IF
END DO ! np_rivers

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE get_rivers_len
!###############################################################################

!###############################################################################
! function rivers_earth_area
!
! Internal procedure in module earth_utils.
!     Calculates area of Earth surface between two lines of latitude and two of
!     longitude, assuming that the Earth is a spheriod.
!     Uses equation A6 of Oki and Sud, 1998, Earth Interactions, Vol.2, Paper 1.

FUNCTION rivers_earth_area( lat1d,lat2d,lon1d,lon2d ) RESULT( area )

IMPLICIT NONE

! Scalar function result

REAL(KIND=real_jlslsm) ::  area
                !  area (m2) (actually in units of earth_radius**2)

! Scalar arguments with intent(in)

REAL(KIND=real_jlslsm), INTENT(IN) :: lat1d
                           !  latitude of southern edge of strip (degrees)
REAL(KIND=real_jlslsm), INTENT(IN) :: lat2d
                           !  latitude of northern edge of strip (degrees)
REAL(KIND=real_jlslsm), INTENT(IN) :: lon1d
                           !  longitude of western edge of strip (degrees)
REAL(KIND=real_jlslsm), INTENT(IN) :: lon2d
                           !  longitude of eastern edge of strip (degrees)

! Local scalars

REAL(KIND=real_jlslsm) ::  esinlat1
                   !  product [eccen * SIN(lat1)], with lat1 in radians
REAL(KIND=real_jlslsm) ::  esinlat2
                   !  product [eccen * SIN(lat2)], with lat2 in radians
REAL(KIND=real_jlslsm) ::  val1, val2, val3      ! work variables

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RIVERS_EARTH_AREA'
!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate SIN(latitudes)
esinlat1 = eccen * SIN(lat1d * pi_over_180)
esinlat2 = eccen * SIN(lat2d * pi_over_180)

! Evaluate terms at each of lat1 and lat2
val1 = 0.5 * esinlat1 / ( 1.0 - (esinlat1 * esinlat1) )
val1 = val1 + 0.25 * LOG( ABS( (1.0 + esinlat1) / (1.0 - esinlat1) ) )

val2 = 0.5 * esinlat2 / ( 1.0 - (esinlat2 * esinlat2) )
val2 = val2 + 0.25 * LOG( ABS( (1.0 + esinlat2) / (1.0 - esinlat2) ) )

val3 = pi_over_180 * planet_eq_radius * planet_eq_radius *                    &
            (1.0 - eccensq) / eccen

area = ABS(lon2d - lon1d) * val3 * ( val2 - val1 )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION rivers_earth_area
!###############################################################################

!###############################################################################
! function givelength
!
! Internal procedure in module earth_utils.
!      Calculates distance on surface of Earth between two locations,
!      assuming Earth is a spheroid.
!      Based on function giveLen by Taikan Oki (26/August/1996).

FUNCTION givelength( lat1,lon1,lat2,lon2 ) RESULT ( length )

  !  Scalar function result
REAL(KIND=real_jlslsm) ::  length   !  distance (m)

!  Scalars with intent (in)
REAL(KIND=real_jlslsm), INTENT(IN) ::  lat1
                            ! latitude of southern edge of strip (degrees)
REAL(KIND=real_jlslsm), INTENT(IN) ::  lat2
                            ! latitude of northern edge of strip (degrees)
REAL(KIND=real_jlslsm), INTENT(IN) ::  lon1
                            ! longitude of western edge of strip (degrees)
REAL(KIND=real_jlslsm), INTENT(IN) ::  lon2
                            ! longitude of eastern edge of strip (degrees)

!  Local scalars
REAL(KIND=real_jlslsm) ::  dlat   !  difference in latitude (degrees)
REAL(KIND=real_jlslsm) ::  dlon   !  difference in longitude (degrees)
REAL(KIND=real_jlslsm) ::  dx     !  work
REAL(KIND=real_jlslsm) ::  dy     !  work
REAL(KIND=real_jlslsm) ::  lat    !  work: latitude (degrees)
REAL(KIND=real_jlslsm) ::  radius
                !  equivalent radius of Earth at given latitude (m)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GIVELENGTH'
!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

dlon = ABS( lon2 - lon1 )
IF ( dlon >= 180.0 ) THEN
  dlon = 360.0 - dlon
END IF
dlat = ABS( lat2 - lat1 )

IF ( dlon < EPSILON( dlon ) ) THEN
  !   Constant longitude.
  lat = ( lat1 + lat2 ) * 0.5
  length = givelatlength(lat) * dlat
ELSE IF ( dlat < EPSILON(dlat) ) THEN
  !   Constant latitude.
  length = givelonlength( lat1 ) * dlon
ELSE
  !   Both lat and lon change.
  !   Use equation A8 of Oki and Sud, 1998, Earth Interactions, Vol.2, Paper 1
  lat = ( lat1 + lat2 ) * 0.5
  radius = giveearthradius( lat )
  dx = givelonlength( lat ) * dlon / radius
  dy = givelatlength( lat ) * dlat / radius
  length = ACOS( COS(dx) * COS(dy) ) * radius
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION givelength

!###############################################################################

!###############################################################################
! function givelatlength
!
! Internal procedure in module earth_utils.
!     Calculates the distance (km) along the surface of the Earth between two
!     points per 1 degree of latitude and at the same longitude.
!     Based on function givelat by Taikan Oki (23/April/1996).
!     See EqnA2 of Oki and Sud, 1998, Earth Interactions, Vol.2, Paper 1.

FUNCTION givelatlength( lat )  RESULT (latlen)

IMPLICIT NONE

! Scalar function result
REAL(KIND=real_jlslsm) ::  latlen               !  distance between points (m)
                              !  [i.e. same units as planet_eq_radius]

! Scalars with intent (in)
REAL(KIND=real_jlslsm), INTENT(IN) ::  lat
                              !  effective (e.g. average) latitude (degrees)

! Local scalars
REAL(KIND=real_jlslsm) :: esinlat  !  product [eccen * SIN(lat_r)]

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GIVELATLENGTH'
!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

esinlat = eccen * SIN ( lat * pi_over_180 )

latlen = pi_over_180 * planet_eq_radius * ( 1 - eccensq)                      &
         / (SQRT(1.0 - ( esinlat * esinlat )))**3

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION givelatlength
!###############################################################################

!###############################################################################
! function givelonlength
!
! Internal procedure in module earth_utils.
!    Calculates the distance (km) along the surface of the Earth between two
!    points per degree of longitude and at the same latitude.
!    Based on function givelat by Taikan Oki (23/April/1996).
!    See EqnA2 of Oki and Sud, 1998, Earth Interactions, Vol.2, Paper 1.

FUNCTION givelonlength( lat )  RESULT (lonlen)

IMPLICIT NONE

! Scalar function result
REAL(KIND=real_jlslsm) ::  lonlen                !  distance between points (m)
                               !  [i.e. same units as planet_eq_radius]

! Scalars with intent (in)
REAL(KIND=real_jlslsm), INTENT(IN) ::  lat       !  latitude (degrees)

! Local scalars
REAL(KIND=real_jlslsm) ::  esinlat               !  product [eccen * SIN(lat_r)]

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GIVELONLENGTH'
!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

esinlat = eccen * SIN ( lat * pi_over_180 )

lonlen = pi_over_180 * planet_eq_radius * COS(lat * pi_over_180)              &
         / SQRT(1.0 - (esinlat * esinlat))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION givelonlength
!###############################################################################

!###############################################################################
! function giveearthradius
!
! Internal procedure in module earth_utils.
!    Calculates equivalent radius of Earth at a given latitude, assuming the
!    Earth is a spheroid.
!    From equations A10 and A11 in Appendix A of Oki and Sud, 1998, Earth
!    Interactions, Vol.2, Paper 1. Based on function giverade by Taikan Oki
!    (26/August/1996).

FUNCTION giveearthradius( lat ) RESULT( radius )

IMPLICIT NONE

! Scalar function result
REAL(KIND=real_jlslsm) ::  radius            !  equivalent radius (m)
                           !  [i.e. same units as planet_eq_radius]

! Scalars with intent (in)
REAL(KIND=real_jlslsm), INTENT(IN) ::  lat   !  latitude (degrees)

! Local scalars
REAL(KIND=real_jlslsm) ::  esinlat           !  product  eccen * SIN(lat_r)
REAL(KIND=real_jlslsm) ::  esinlat2          !  product [eccen * SIN(lat_r)]^2
REAL(KIND=real_jlslsm) ::  rn                !  work

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GIVEEARTHRADIUS'
!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

esinlat  = eccen * SIN(lat * pi_over_180)
esinlat2 = esinlat * esinlat

rn = planet_eq_radius / SQRT(1.0 - esinlat2)
radius = rn * SQRT(1.0-2.0 * eccen * esinlat +                                &
                   eccensq * esinlat2)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION giveearthradius
!###############################################################################

END MODULE rivers_utils
!###############################################################################
