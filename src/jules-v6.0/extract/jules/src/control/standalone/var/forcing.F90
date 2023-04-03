#if !defined(UM_JULES)
! Module containing all of the driving (atmospheric forcing) variables.

MODULE forcing

IMPLICIT NONE

!-------------------------------------------------------------------------------

! The forcing variables.
REAL, ALLOCATABLE ::                                                          &
  qw_1_ij(:,:),                                                               &
                      !  Total water content (Kg/Kg)
  tl_1_ij(:,:),                                                               &
                      ! Ice/liquid water temperature (k)
  u_0_ij(:,:),                                                                &
                      ! W'ly component of surface current (m/s)
  v_0_ij(:,:),                                                                &
                      ! S'ly component of surface current (m/s)
  u_1_ij(:,:),                                                                &
                      ! W'ly wind component (m/s)
  v_1_ij(:,:),                                                                &
                      ! S'ly wind component (m/s)
  pstar_ij(:,:),                                                              &
                      ! Surface pressure (Pascals)
  ls_rain_ij(:,:),                                                            &
                      ! Large-scale rain (kg/m2/s)
  con_rain_ij(:,:),                                                           &
                      ! Convective rain (kg/m2/s)
  ls_snow_ij(:,:),                                                            &
                      ! Large-scale snowfall (kg/m2/s)
  con_snow_ij(:,:),                                                           &
                      ! Convective snowfall (kg/m2/s)
  sw_down_ij(:,:),                                                            &
                      ! Surface downward SW radiation (W/m2)
  lw_down_ij(:,:),                                                            &
                      ! Surface downward LW radiation (W/m2)
  diurnal_temperature_range_ij(:,:)
                      ! diurnal temperature range (K), used when
                      !l_dailydisagg=T

! Variables that aid in the calculation of the actual forcing variables
REAL, ALLOCATABLE ::                                                          &
  diff_rad_ij(:,:)       ! Input diffuse radiation (W/m2)

!-------------------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FORCING'

CONTAINS
  
SUBROUTINE forcing_alloc(t_i_length,t_j_length)

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: t_i_length,t_j_length

!Local variables
INTEGER :: temp_size, temp_tiles, temp_layers

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FORCING_ALLOC'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE( qw_1_ij(t_i_length,t_j_length))
ALLOCATE( tl_1_ij(t_i_length,t_j_length))
ALLOCATE( u_0_ij(t_i_length,t_j_length))
ALLOCATE( v_0_ij(t_i_length,t_j_length))
ALLOCATE( u_1_ij(t_i_length,t_j_length))
ALLOCATE( v_1_ij(t_i_length,t_j_length))
ALLOCATE( pstar_ij(t_i_length,t_j_length))
ALLOCATE( ls_rain_ij(t_i_length,t_j_length))
ALLOCATE( con_rain_ij(t_i_length,t_j_length))
ALLOCATE( ls_snow_ij(t_i_length,t_j_length))
ALLOCATE( con_snow_ij(t_i_length,t_j_length))
ALLOCATE( sw_down_ij(t_i_length,t_j_length))
ALLOCATE( lw_down_ij(t_i_length,t_j_length))
ALLOCATE( diff_rad_ij(t_i_length,t_j_length))
ALLOCATE( diurnal_temperature_range_ij(t_i_length,t_j_length))

qw_1_ij(:,:)                       = 0.0
tl_1_ij(:,:)                       = 0.0
u_0_ij(:,:)                        = 0.0
v_0_ij(:,:)                        = 0.0
u_1_ij(:,:)                        = 0.0
v_1_ij(:,:)                        = 0.0
pstar_ij(:,:)                      = 0.0
ls_rain_ij(:,:)                    = 0.0
con_rain_ij(:,:)                   = 0.0
ls_snow_ij(:,:)                    = 0.0
con_snow_ij(:,:)                   = 0.0
sw_down_ij(:,:)                    = 0.0
lw_down_ij(:,:)                    = 0.0
diff_rad_ij(:,:)                   = 0.0
diurnal_temperature_range_ij(:,:)  = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE forcing_alloc

END MODULE forcing
#endif
