! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module holding parameter arrays for non-vegetation surface types.


! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.


MODULE nvegparm

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
 albsnc_nvg(:)                                                                &
                  ! Snow-covered albedo.
,albsnf_nvgu(:)                                                               &
                  ! Max Snow-free albedo, when scaled to obs
,albsnf_nvg(:)                                                                &
                  ! Snow-free albedo.
,albsnf_nvgl(:)                                                               &
                  ! Min Snow-free albedo, when scaled to obs
,catch_nvg(:)                                                                 &
                  ! Canopy capacity for water (kg/m2).
,gs_nvg(:)                                                                    &
                  ! Surface conductance (m/s).
,infil_nvg(:)                                                                 &
                  ! Infiltration enhancement factor.
,z0_nvg(:)                                                                    &
                  ! Roughness length (m).
,ch_nvg(:)                                                                    &
                  ! "Canopy" heat capacity (J/K/m2)
,vf_nvg(:)                                                                    &
                  ! Fractional "canopy" coverage
,emis_nvg(:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NVEGPARM'

CONTAINS
  
SUBROUTINE nvegparm_alloc(nnvg)

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: nnvg

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NVEGPARM_ALLOC'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE( albsnc_nvg(nnvg))
ALLOCATE( albsnf_nvgu(nnvg))
ALLOCATE( albsnf_nvg(nnvg))
ALLOCATE( albsnf_nvgl(nnvg))
ALLOCATE( catch_nvg(nnvg))
ALLOCATE( emis_nvg(nnvg))
ALLOCATE( gs_nvg(nnvg))
ALLOCATE( infil_nvg(nnvg))
ALLOCATE( z0_nvg(nnvg))
ALLOCATE( ch_nvg(nnvg))
ALLOCATE( vf_nvg(nnvg))

albsnc_nvg(:)  = 0.0
albsnf_nvgu(:) = 0.0
albsnf_nvg(:)  = 0.0
albsnf_nvgl(:) = 0.0
catch_nvg(:)   = 0.0
emis_nvg(:)    = 0.0
gs_nvg(:)      = 0.0
infil_nvg(:)   = 0.0
z0_nvg(:)      = 0.0
ch_nvg(:)      = 0.0
vf_nvg(:)      = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE nvegparm_alloc

END MODULE nvegparm
