!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


MODULE ref_sounding

  use ModNamelistFile, only: namelistFile

   USE grid_dims, ONLY: &
        nzpmax,         & ! INTENT(IN)
        maxgrds           ! INTENT(IN)

   PRIVATE

   public :: StoreNamelistFileAtRef_sounding

   INTEGER, PUBLIC             :: iref
   INTEGER, PUBLIC             :: jref
   REAL, PUBLIC                :: topref

!!$   real, public                :: u01dn(nzpmax,maxgrds)
!!$   real, public                :: v01dn(nzpmax,maxgrds)
!!$   real, public                :: pi01dn(nzpmax,maxgrds)
!!$   real, public                :: th01dn(nzpmax,maxgrds)
!!$   real, public                :: dn01dn(nzpmax,maxgrds)
!!$   real, public                :: rt01dn(nzpmax,maxgrds)
   REAL, ALLOCATABLE, PUBLIC      :: u01dn(:,:)
   REAL, ALLOCATABLE, PUBLIC      :: v01dn(:,:)
   REAL, ALLOCATABLE, PUBLIC, target :: pi01dn(:,:)
   REAL, ALLOCATABLE, PUBLIC, target :: th01dn(:,:)
   REAL, ALLOCATABLE, PUBLIC       :: dn01dn(:,:)
   REAL, ALLOCATABLE, PUBLIC      :: rt01dn(:,:)

   INTEGER, PARAMETER, PUBLIC  :: maxsndg=200
   INTEGER, PUBLIC             :: ipsflg ! from RAMSIN
   INTEGER, PUBLIC             :: itsflg ! from RAMSIN
   INTEGER, PUBLIC             :: irtsflg ! from RAMSIN
   INTEGER, PUBLIC             :: iusflg ! from RAMSIN
   INTEGER, PUBLIC             :: nsndg
   REAL, PUBLIC                :: us(maxsndg) ! from RAMSIN
   REAL, PUBLIC                :: vs(maxsndg) ! from RAMSIN
   REAL, PUBLIC                :: ts(maxsndg) ! from RAMSIN
   REAL, PUBLIC                :: thds(maxsndg)
   REAL, PUBLIC                :: ps(maxsndg) ! from RAMSIN
   REAL, PUBLIC                :: hs(maxsndg) ! from RAMSIN
   REAL, PUBLIC                :: rts(maxsndg) ! from RAMSIN

   PUBLIC :: createRefSounding, destroyRefSounding

 CONTAINS

   SUBROUTINE createRefSounding(ngrids, nnzp)
     IMPLICIT NONE
     ! Arguments:
     INTEGER, INTENT(in) :: ngrids
     INTEGER, INTENT(IN) :: nnzp(maxgrds) ! From RAMSIN
     ! Local Variables:
     INTEGER :: ierr, maxz

     maxz = MAXVAL(nnzp(1:ngrids))
     ALLOCATE(u01dn(maxz,ngrids), STAT=ierr)
     IF (ierr/=0) CALL fatal_error&
          ("ERROR allocating u01dn (createRefSounding)")
     ALLOCATE(v01dn(maxz,ngrids), STAT=ierr)
     IF (ierr/=0) CALL fatal_error&
          ("ERROR allocating v01dn (createRefSounding)")
     ALLOCATE(pi01dn(maxz,ngrids), STAT=ierr)
     IF (ierr/=0) CALL fatal_error&
          ("ERROR allocating pi01dn (createRefSounding)")
     ALLOCATE(th01dn(maxz,ngrids), STAT=ierr)
     IF (ierr/=0) CALL fatal_error&
          ("ERROR allocating th01dn (createRefSounding)")
     ALLOCATE(dn01dn(maxz,ngrids), STAT=ierr)
     IF (ierr/=0) CALL fatal_error&
          ("ERROR allocating dn01dn (createRefSounding)")
     ALLOCATE(rt01dn(maxz,ngrids), STAT=ierr)
     IF (ierr/=0) CALL fatal_error&
          ("ERROR allocating rt01dn (createRefSounding)")

   END SUBROUTINE createRefSounding

   ! ************************************************************

   SUBROUTINE destroyRefSounding()
     IMPLICIT NONE
     ! Local Variables:
     INTEGER :: ierr

     DEALLOCATE(u01dn, STAT=ierr)
     IF (ierr/=0) CALL fatal_error&
          ("ERROR deallocating u01dn (destroyRefSounding)")
     DEALLOCATE(v01dn, STAT=ierr)
     IF (ierr/=0) CALL fatal_error&
          ("ERROR deallocating v01dn (destroyRefSounding)")
     DEALLOCATE(pi01dn, STAT=ierr)
     IF (ierr/=0) CALL fatal_error&
          ("ERROR deallocating pi01dn (destroyRefSounding)")
     DEALLOCATE(th01dn, STAT=ierr)
     IF (ierr/=0) CALL fatal_error&
          ("ERROR deallocating th01dn (destroyRefSounding)")
     DEALLOCATE(dn01dn, STAT=ierr)
     IF (ierr/=0) CALL fatal_error&
          ("ERROR deallocating dn01dn (destroyRefSounding)")
     DEALLOCATE(rt01dn, STAT=ierr)
     IF (ierr/=0) CALL fatal_error&
          ("ERROR deallocating rt01dn (destroyRefSounding)")
     
   END SUBROUTINE destroyRefSounding

  subroutine StoreNamelistFileAtRef_sounding(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    hs = oneNamelistFile%hs
    ipsflg = oneNamelistFile%ipsflg
    irtsflg = oneNamelistFile%irtsflg
    itsflg = oneNamelistFile%itsflg
    iusflg = oneNamelistFile%iusflg
    ps = oneNamelistFile%ps
    rts = oneNamelistFile%rts
    ts = oneNamelistFile%ts
    us = oneNamelistFile%us
    vs = oneNamelistFile%vs
  end subroutine StoreNamelistFileAtRef_sounding
END MODULE ref_sounding
