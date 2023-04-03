! File containing variables for ozone implementation

MODULE ozone_vars

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

! Ozone forcing
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  o3_gb(:)        ! Surface ozone concentration (ppb).


! Ozone diagnostics
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  flux_o3_pft(:,:)                                                            &
               ! Flux of O3 to stomata (nmol O3/m2/s).
 ,fo3_pft(:,:)  ! Ozone exposure factor.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='OZONE_VARS'

CONTAINS
 
SUBROUTINE ozone_vars_alloc(land_pts,npft)

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: land_pts,npft

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='OZONE_VARS_ALLOC'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE( o3_gb(land_pts))
ALLOCATE( flux_o3_pft(land_pts,npft))
ALLOCATE( fo3_pft(land_pts,npft))

o3_gb(:)         = 0.0
flux_o3_pft(:,:) = 0.0
fo3_pft(:,:)     = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ozone_vars_alloc

END MODULE ozone_vars
