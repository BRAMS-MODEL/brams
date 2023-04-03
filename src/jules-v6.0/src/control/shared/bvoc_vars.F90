! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

! Module containing BVOC diagnostics

MODULE bvoc_vars

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  isoprene_gb(:),                                                             &
          ! Gridbox mean isoprene emission flux (kgC/m2/s)
  isoprene_pft(:,:),                                                          &
          ! Isoprene emission flux on PFTs (kgC/m2/s)
  terpene_gb(:),                                                              &
          ! Gridbox mean (mono-)terpene emission flux (kgC/m2/s)
  terpene_pft(:,:),                                                           &
          ! (Mono-)Terpene emission flux on PFTs (kgC/m2/s)
  methanol_gb(:),                                                             &
          ! Gridbox mean methanol emission flux (kgC/m2/s)
  methanol_pft(:,:),                                                          &
          ! Methanol emission flux on PFTs (kgC/m2/s)
  acetone_gb(:),                                                              &
          ! Gridbox mean acetone emission flux (kgC/m2/s)
  acetone_pft(:,:)
          ! Acetone emission flux on PFTs (kgC/m2/s)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='BVOC_VARS'

CONTAINS

SUBROUTINE bvocvars_alloc(land_pts,npft)

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: land_pts,npft

!Local variables
INTEGER :: temp_size, temp_tiles, temp_layers

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BVOCVARS_ALLOC'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!  ====bvoc_vars module common====
ALLOCATE( isoprene_gb(land_pts))
ALLOCATE( isoprene_pft(land_pts,npft))
ALLOCATE( terpene_gb(land_pts))
ALLOCATE( terpene_pft(land_pts,npft))
ALLOCATE( methanol_gb(land_pts))
ALLOCATE( methanol_pft(land_pts,npft))
ALLOCATE( acetone_gb(land_pts))
ALLOCATE( acetone_pft(land_pts,npft))

isoprene_gb(:)    = 0.0
isoprene_pft(:,:) = 0.0
terpene_gb(:)     = 0.0
terpene_pft(:,:)  = 0.0
methanol_gb(:)    = 0.0
methanol_pft(:,:) = 0.0
acetone_gb(:)     = 0.0
acetone_pft(:,:)  = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE bvocvars_alloc
END MODULE bvoc_vars
