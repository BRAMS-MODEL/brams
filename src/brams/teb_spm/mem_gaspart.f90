MODULE mem_gaspart

  TYPE gaspart_vars

     ! Variables to be dimensioned by (nzp,nxp,nyp)
     REAL, POINTER :: PCO(:,:,:)
     REAL, POINTER :: PNO(:,:,:)
     REAL, POINTER :: PNO2(:,:,:)
     REAL, POINTER :: PPM25(:,:,:)
     REAL, POINTER :: PVOC(:,:,:)
     REAL, POINTER :: PSO2(:,:,:)
     REAL, POINTER :: PROO(:,:,:)
     REAL, POINTER :: PSO4(:,:,:)
     REAL, POINTER :: PAER(:,:,:)
     REAL, POINTER :: PO3(:,:,:)
     REAL, POINTER :: PRHCO(:,:,:)
     REAL, POINTER :: PHO2(:,:,:)
     REAL, POINTER :: PO3P(:,:,:)
     REAL, POINTER :: PO1D(:,:,:)
     REAL, POINTER :: PHO(:,:,:)
     REAL, POINTER :: GASR(:,:,:)
     REAL, POINTER :: PEOXID(:,:,:)

     ! Variables to be dimensioned by (nzp,nxp)
     REAL, POINTER :: fusog(:,:)

     REAL, POINTER :: PCOT(:)
     REAL, POINTER :: PNOT(:)
     REAL, POINTER :: PNO2T(:)
     REAL, POINTER :: PPM25T(:)
     REAL, POINTER :: PVOCT(:)
     REAL, POINTER :: PSO2T(:)
     REAL, POINTER :: PSO4T(:)
     REAL, POINTER :: PAERT(:) 
     REAL, POINTER :: PO3T(:)
     REAL, POINTER :: PRHCOT(:)
     REAL, POINTER :: PHO2T(:)
     REAL, POINTER :: PO3PT(:)
     REAL, POINTER :: PO1DT(:)
     REAL, POINTER :: PHOT(:)
     REAL, POINTER :: PROOT(:)

  END TYPE gaspart_vars

  TYPE (gaspart_vars), ALLOCATABLE, TARGET :: gaspart_g(:), gaspartm_g(:)

CONTAINS

  SUBROUTINE alloc_gaspart(gaspart, n1, n2, n3)

    USE mem_emiss, ONLY: ichemi

    IMPLICIT NONE
    ! Arguments:
    TYPE (gaspart_vars), INTENT(INOUT) :: gaspart
    INTEGER, INTENT(in) :: n1, n2, n3

    ALLOCATE (gaspart%pco(n1,n2,n3))
    ALLOCATE (gaspart%pno(n1,n2,n3))
    ALLOCATE (gaspart%pno2(n1,n2,n3))
    ALLOCATE (gaspart%ppm25(n1,n2,n3))
    ALLOCATE (gaspart%pso2(n1,n2,n3))
    ALLOCATE (gaspart%pvoc(n1,n2,n3))
    ALLOCATE (gaspart%gasr(n1,n2,n3))
    ALLOCATE (gaspart%pso4(n1,n2,n3))
    ALLOCATE (gaspart%paer(n1,n2,n3))
    ALLOCATE (gaspart%PEOXID(n1,n2,n3))
    ALLOCATE (gaspart%fusog(n2,n3))

    if (ichemi==1) then
       ALLOCATE (gaspart%po3   (n1,n2,n3))
       ALLOCATE (gaspart%prhco (n1,n2,n3))
       ALLOCATE (gaspart%pho2  (n1,n2,n3))
       ALLOCATE (gaspart%po3p  (n1,n2,n3))
       ALLOCATE (gaspart%po1d  (n1,n2,n3))
       ALLOCATE (gaspart%pho   (n1,n2,n3))
       ALLOCATE (gaspart%proo  (n1,n2,n3))
    endif

  END SUBROUTINE alloc_gaspart



  SUBROUTINE zero_gaspart(gaspart, n1, n2, n3)

    USE mem_emiss, ONLY: ichemi

    IMPLICIT NONE
    ! Arguments:
    TYPE (gaspart_vars), INTENT(OUT) :: gaspart
    INTEGER, INTENT(in) :: n1, n2, n3

    gaspart%pco    = 0.0
    gaspart%pno    = 0.0
    gaspart%pno2   = 0.0
    gaspart%ppm25  = 0.0
    gaspart%pso2   = 0.0
    gaspart%pvoc   = 0.0
    gaspart%gasr   = 0.0
    gaspart%pso4   = 0.0
    gaspart%paer   = 0.0
    gaspart%PEOXID = 0.0
    gaspart%fusog  = 0.0

    if (ichemi==1) then
       gaspart%po3   = 0.0
       gaspart%prhco = 0.0
       gaspart%pho2  = 0.0
       gaspart%po3p  = 0.0
       gaspart%po1d  = 0.0
       gaspart%pho   = 0.0
       gaspart%proo  = 0.0
    endif

  END SUBROUTINE zero_gaspart


  SUBROUTINE nullify_gaspart(gaspart)
    USE mem_emiss, ONLY: ichemi

    IMPLICIT NONE
    ! Arguments:
    TYPE (gaspart_vars), INTENT(INOUT) :: gaspart

    IF (ASSOCIATED(gaspart%pco))    NULLIFY (gaspart%pco)
    IF (ASSOCIATED(gaspart%pno))    NULLIFY (gaspart%pno)
    IF (ASSOCIATED(gaspart%pno2))   NULLIFY (gaspart%pno2)
    IF (ASSOCIATED(gaspart%ppm25))  NULLIFY (gaspart%ppm25)
    IF (ASSOCIATED(gaspart%pvoc))   NULLIFY (gaspart%pvoc)
    IF (ASSOCIATED(gaspart%pso2))   NULLIFY (gaspart%pso2)
    IF (ASSOCIATED(gaspart%pso4))   NULLIFY (gaspart%pso4)
    IF (ASSOCIATED(gaspart%paer))   NULLIFY (gaspart%paer)
    IF (ASSOCIATED(gaspart%PEOXID))   NULLIFY (gaspart%PEOXID)
    IF (ASSOCIATED(gaspart%gasr))   NULLIFY (gaspart%gasr)
    IF (ASSOCIATED(gaspart%fusog))   NULLIFY (gaspart%fusog)

    if (ichemi==1) then
       IF (ASSOCIATED(gaspart%po3  ))  NULLIFY (gaspart%po3  )
       IF (ASSOCIATED(gaspart%prhco))  NULLIFY (gaspart%prhco)
       IF (ASSOCIATED(gaspart%pho2 ))  NULLIFY (gaspart%pho2 )
       IF (ASSOCIATED(gaspart%po3p ))  NULLIFY (gaspart%po3p )
       IF (ASSOCIATED(gaspart%po1d ))  NULLIFY (gaspart%po1d )
       IF (ASSOCIATED(gaspart%pho  ))  NULLIFY (gaspart%pho  )
       IF (ASSOCIATED(gaspart%proo ))  NULLIFY (gaspart%proo )
    endif

  END SUBROUTINE nullify_gaspart



  SUBROUTINE dealloc_gaspart(gaspart)

    USE mem_emiss, ONLY: ichemi

    IMPLICIT NONE

    ! Arguments
    TYPE (gaspart_vars), INTENT(INOUT) :: gaspart

    IF (ASSOCIATED(gaspart%pco))    DEALLOCATE (gaspart%pco)
    IF (ASSOCIATED(gaspart%pno))    DEALLOCATE (gaspart%pno)
    IF (ASSOCIATED(gaspart%pno2))   DEALLOCATE (gaspart%pno2)
    IF (ASSOCIATED(gaspart%ppm25))  DEALLOCATE (gaspart%ppm25)
    IF (ASSOCIATED(gaspart%pvoc))   DEALLOCATE (gaspart%pvoc)
    IF (ASSOCIATED(gaspart%pso2))   DEALLOCATE (gaspart%pso2)
    IF (ASSOCIATED(gaspart%pso4))   DEALLOCATE (gaspart%pso4)
    IF (ASSOCIATED(gaspart%paer))   DEALLOCATE (gaspart%paer)
    IF (ASSOCIATED(gaspart%PEOXID))   DEALLOCATE (gaspart%PEOXID)
    IF (ASSOCIATED(gaspart%gasr))   DEALLOCATE (gaspart%gasr)
    IF (ASSOCIATED(gaspart%fusog))   DEALLOCATE (gaspart%fusog)

    if (ichemi==1) then
       IF (ASSOCIATED(gaspart%po3  ))  DEALLOCATE (gaspart%po3  )
       IF (ASSOCIATED(gaspart%prhco))  DEALLOCATE (gaspart%prhco)
       IF (ASSOCIATED(gaspart%pho2 ))  DEALLOCATE (gaspart%pho2 )
       IF (ASSOCIATED(gaspart%po3p ))  DEALLOCATE (gaspart%po3p )
       IF (ASSOCIATED(gaspart%po1d ))  DEALLOCATE (gaspart%po1d )
       IF (ASSOCIATED(gaspart%pho  ))  DEALLOCATE (gaspart%pho  )
       IF (ASSOCIATED(gaspart%proo ))  DEALLOCATE (gaspart%proo )
    endif

  END SUBROUTINE dealloc_gaspart


  SUBROUTINE filltab_gaspart(gaspart, gaspartm, imean, n1, n2, n3, ng)

    USE var_tables
    USE mem_emiss, ONLY: ichemi

    IMPLICIT NONE
    include "i8.h"
    ! Arguments:
    TYPE (gaspart_vars), INTENT(IN) :: gaspart, gaspartm
    INTEGER, INTENT(IN) :: imean, n1, n2, n3, ng
    ! Local Variables:
    INTEGER(kind=i8)  :: npts
    REAL, POINTER :: var, varm

    ! Fill pointers to arrays into variable tables

    npts=n2*n3

    IF (ASSOCIATED(gaspart%fusog))  &
         CALL InsertVTab (gaspart%FUSOG,gaspartm%FUSOG &
         ,ng, npts, imean,  &
         'FUSOG:2:hist:anal:mpti:mpt3:mpt1')

    npts=n1*n2*n3

    IF (ASSOCIATED(gaspart%pco))  &
         CALL InsertVTab (gaspart%PCO,gaspartm%PCO&
         ,ng, npts, imean,  &
         'PCO:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%pno))  &
         CALL InsertVTab (gaspart%PNO,gaspartm%PNO&
         ,ng, npts, imean,  &
         'PNO:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%pno2))  &
         CALL InsertVTab (gaspart%PNO2,gaspartm%PNO2&
         ,ng, npts, imean,  &
         'PNO2:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%ppm25))  &
         CALL InsertVTab (gaspart%PPM25,gaspartm%PPM25&
         ,ng, npts, imean,  &
         'PPM25:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%pvoc))  &
         CALL InsertVTab (gaspart%PVOC,gaspartm%PVOC&
         ,ng, npts, imean,  &
         'PVOC:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%pso2))  &
         CALL InsertVTab (gaspart%PSO2,gaspartm%PSO2&
         ,ng, npts, imean,  &
         'PSO2:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%pso4))  &
         CALL InsertVTab (gaspart%PSO4,gaspartm%PSO4&
         ,ng, npts, imean,  &
         'PSO4:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%paer))  &
         CALL InsertVTab (gaspart%PAER,gaspartm%PAER&
         ,ng, npts, imean,  &
         'PAER:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%PEOXID))  &
         CALL InsertVTab (gaspart%PEOXID,gaspartm%PEOXID&
         ,ng, npts, imean,  &
         'PEOXID:3:hist:anal:mpti:mpt3:mpt1')
    IF (ASSOCIATED(gaspart%gasr))  &
         CALL InsertVTab (gaspart%GASR,gaspartm%GASR&
         ,ng, npts, imean,  &
         'GASR:3:mpti:mpt3:mpt1')

    if (ichemi==1) then
       IF (ASSOCIATED(gaspart%po3))  &
            CALL InsertVTab (gaspart%PO3,gaspartm%PO3&
            ,ng, npts, imean,  &
            'PO3:3:hist:anal:mpti:mpt3:mpt1')
       IF (ASSOCIATED(gaspart%prhco))  &
            CALL InsertVTab (gaspart%PRHCO,gaspartm%PRHCO&
            ,ng, npts, imean,  &
            'PRHCO:3:hist:anal:mpti:mpt3:mpt1')
       IF (ASSOCIATED(gaspart%pho2))  &
            CALL InsertVTab (gaspart%PHO2,gaspartm%PHO2&
            ,ng, npts, imean,  &
            'PHO2:3:hist:anal:mpti:mpt3:mpt1')
       IF (ASSOCIATED(gaspart%po3p))  &
            CALL InsertVTab (gaspart%PO3P,gaspartm%PO3P&
            ,ng, npts, imean,  &
            'PO3P:3:hist:anal:mpti:mpt3:mpt1')
       IF (ASSOCIATED(gaspart%po1d))  &
            CALL InsertVTab (gaspart%PO1D,gaspartm%PO1D&
            ,ng, npts, imean,  &
            'PO1D:3:hist:anal:mpti:mpt3:mpt1')
       IF (ASSOCIATED(gaspart%pho))  &
            CALL InsertVTab (gaspart%PHO,gaspartm%PHO&
            ,ng, npts, imean,  &
            'PHO:3:hist:anal:mpti:mpt3:mpt1')
       IF (ASSOCIATED(gaspart%proo))  &
            CALL InsertVTab (gaspart%PROO,gaspartm%PROO&
            ,ng, npts, imean,  &
            'PROO:3:hist:anal:mpti:mpt3:mpt1')
    endif

  END SUBROUTINE filltab_gaspart

END MODULE mem_gaspart
