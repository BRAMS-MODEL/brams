
MODULE init_cable_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE init_cable_grid()

USE cable_types_mod,          ONLY: l_tile_pts
USE ancil_info,               ONLY: land_pts
USE jules_fields_mod,         ONLY: ainfo
USE jules_surface_types_mod,  ONLY: ntype

IMPLICIT NONE

!------------------------------------------------------------------------------
! Description:
!   Initialises the JULES/CABLE grid array, which aligns JULES grid points
!   with CABLE land points
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!------------------------------------------------------------------------------

INTEGER :: i, j

ALLOCATE(l_tile_pts(land_pts, ntype))

l_tile_pts(:,:) = .FALSE.

DO j = 1, ntype
  DO i = 1, land_pts
    IF ( ainfo%frac_surft(i,j)  >   0.0 ) THEN
      l_tile_pts(i,j) = .TRUE.
    END IF
  END DO
END DO

RETURN

END SUBROUTINE init_cable_grid


SUBROUTINE init_cable_veg()

USE cable_types_mod,         ONLY: veg, vegin, mp, l_tile_pts
USE jules_fields_mod,        ONLY: ainfo

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   init_cables veg parameters using values read from namelist
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

INTEGER :: h

veg%iveg   = PACK(ainfo%surft_index, l_tile_pts)

! Prescribe parameters for current gridcell based on veg/soil type
! (which may have loaded from default value file or met file):
DO h = 1, mp          ! over each patch in current grid
  veg%taul(h,1)   = vegin%taul(1,veg%iveg(h))
  veg%taul(h,2)   = vegin%taul(2,veg%iveg(h))
  veg%refl(h,1)   = vegin%refl(1,veg%iveg(h))
  veg%refl(h,2)   = vegin%refl(2,veg%iveg(h))
END DO ! over each veg patch in land point

END SUBROUTINE init_cable_veg

END MODULE init_cable_mod

