MODULE solve_sparse

  IMPLICIT NONE

  PRIVATE

  TYPE solve_param
     LOGICAL :: initialized
     INTEGER :: size
     INTEGER :: non_zero
  END TYPE solve_param

  TYPE(solve_param),PRIVATE :: Lmatrix
  INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:) :: element
  INTEGER,ALLOCATABLE,DIMENSION(:) :: i_pos
  INTEGER,ALLOCATABLE,DIMENSION(:) :: j_pos
  INTEGER(KIND=8) :: matrixId

  DOUBLE PRECISION,PRIVATE,PARAMETER :: zero=1.0D-39

  INTEGER, PUBLIC :: status
  LOGICAL, PUBLIC :: get_non_zeros=.FALSE.

  PUBLIC :: Prepare,      & ! Subroutine
            Solve_linear, & ! Subroutine
            Eliminate_0,  & ! Subroutine
            reserve         ! Subroutine

CONTAINS

  !LFR  !=================================
  !LFR  SUBROUTINE Alloc_matrix(ijk,sizeOfMatrix)
  !LFR  !=================================
  !LFR    IMPLICIT NONE
  !LFR    INTEGER,INTENT(IN) :: ijk,sizeOfMatrix
  !LFR    
  !LFR    ALLOCATE(DLMat(ijk,sizeOfMatrix,sizeOfMatrix))
  !LFR    ALLOCATE(DLB1(ijk,sizeOfMatrix))
  !LFR    ALLOCATE(DLB2(ijk,sizeOfMatrix))
  !LFR    Lmatrix%size=sizeOfMatrix
  !LFR
  !LFR  END SUBROUTINE Alloc_matrix


  !=================================
  SUBROUTINE Prepare(sizeOfMatrix)
  !=================================
      
    INTEGER, INTENT(IN) :: sizeOfMatrix
    
    INTEGER(kind=8), EXTERNAL :: sfCreate_solve
    INTEGER :: error
    INTEGER,PARAMETER :: complex_t=0


    !create matrix
    matrixId = sfCreate_solve(sizeOfMatrix,complex_t,error)

    Lmatrix%initialized=.TRUE.
    IF(ALLOCATED(i_pos)) DEALLOCATE (i_pos)
    ALLOCATE(i_pos(sizeofmatrix*sizeofmatrix))

    IF(ALLOCATED(j_pos)) DEALLOCATE (j_pos)
    ALLOCATE(j_pos(sizeOfMatrix*sizeOfMatrix))


  END SUBROUTINE Prepare

  !================================================
  SUBROUTINE Eliminate_0(matrixElements,sizeOfMatrix)
  !================================================

    INTEGER                                               , INTENT(IN) :: sizeOfMatrix
    DOUBLE PRECISION,DIMENSION(sizeOfMatrix,sizeOfMatrix) , INTENT(IN) :: matrixElements

    INTEGER :: i,j

    IF(.NOT. Lmatrix%initialized) THEN 
       PRINT *,'ERROR - In Module Solve - Line 68 - solve.f90'
       PRINT *,'Matrix not initialized - please verify your code'
       CALL flush(6)
       STOP
    END IF

    Lmatrix%non_zero=0
    DO j=1,sizeOfMatrix
       DO i=1,sizeOfMatrix
          IF(matrixElements(i,j)<=zero .AND. matrixElements(i,j)>=(-1)*Zero) CYCLE
          Lmatrix%non_zero=Lmatrix%non_zero+1
          i_pos(Lmatrix%non_zero)=i
          j_pos(Lmatrix%non_zero)=j
       END DO
    END DO
    IF(ALLOCATED(element)) DEALLOCATE (element)
    ALLOCATE(element(lmatrix%non_zero))

  END SUBROUTINE Eliminate_0

  !==============================================
  SUBROUTINE Reserve()
  !=============================================

    INTEGER :: i

    INTEGER(kind=8), EXTERNAL :: sfGetElement

    !if(ijk ==1) then
    DO i=1,Lmatrix%non_zero
       element(i)=sfGetElement(matrixId,i_pos(i),j_pos(i))
    END DO
    !endif
    !clear matrix
    CALL sfZero(matrixId)


  END SUBROUTINE Reserve


  !=====================================================
  SUBROUTINE Load_matrix(matrixElements,sizeOfMatrix)
  !=====================================================

    INTEGER                                               , INTENT(IN) :: sizeOfMatrix
    DOUBLE PRECISION,DIMENSION(sizeOfMatrix,sizeOfMatrix) , INTENT(IN) :: matrixElements

    INTEGER :: i

    DO i=1,Lmatrix%non_zero
       CALL sfAdd1Real(element(i),matrixElements(i_pos(i),j_pos(i)))
    END DO

  END SUBROUTINE Load_matrix

  !call sfprint(matrix, .false., .false., .true.)


  !=============================================================================
  SUBROUTINE Fact_solve_matrix(matrixIndependent,matrixSolution,sizeOfMatrix)
  !=============================================================================

    INTEGER                                               , INTENT(IN)    :: sizeOfMatrix
    DOUBLE PRECISION,DIMENSION(sizeOfMatrix)              , INTENT(IN)    :: matrixIndependent
    DOUBLE PRECISION,DIMENSION(sizeOfMatrix,sizeOfMatrix) , INTENT(OUT) :: matrixSolution

    INTEGER, EXTERNAL :: sfFactor
    INTEGER :: error

    INTEGER, PARAMETER :: spOKAY       =0
    INTEGER, PARAMETER :: spSMALL_PIVOT=1
    INTEGER, PARAMETER :: spZERO_DIAG  =2
    INTEGER, PARAMETER :: spSINGULAR   =3
    INTEGER, PARAMETER :: spNO_MEMORY  =4
    INTEGER, PARAMETER :: spPANIC      =5

    !factor matrix
    error = sfFactor(matrixId)
    status=error
    IF(error/=0) THEN
       PRINT *,'Error detected - Factor Matrix - sfFactor',matrixId
       PRINT *,'See below:'
       SELECT CASE (error)
       CASE (spSMALL_PIVOT)
          PRINT *, ' WARNING(Non Fatal) When reordering the matrix, no element was found '
          PRINT *, ' which satisfies the absolute threshold criteria. The largest element'
          PRINT *, ' in the matrix was chosen as pivot.'
          CALL flush(6)
       CASE (spZERO_DIAG)
          PRINT *, 'FATAL ERROR.  A zero was encountered on the diagonal the matrix.'
          PRINT *, 'This does not necessarily imply that the matrix is singular. When'
          PRINT *, 'this error occurs, the matrix should be reconstructed and factored'
          PRINT *, 'using spOrderAndFactor()'
          STOP 444
          !	    CALL flush(6)
          !	    RETURN
       CASE (spSINGULAR)
          PRINT *, 'FATAL ERROR.  Matrix is singular, so no unique solution exists.'
          RETURN
       CASE(spNO_MEMORY)
          PRINT *, 'FATAL ERROR.  Not enough memory is available to '
          PRINT *, 'handle the matrix.'
          CALL flush(6)
          RETURN
       CASE(spPANIC)
          PRINT *, 'FATAL ERROR. Routines are not prepared to  handle the matrix'
          PRINT *, 'that has been requested.  This may occur when the matrix is '
          PRINT *, 'specified to be real and the routines are not compiled for real'
          PRINT *, 'matrices, or when the matrix is specified to be complex and the'
          PRINT *, 'routines are not compiled to handle complex matrices.'
          CALL flush(6)
          RETURN
       CASE DEFAULT
          PRINT *, 'WARNING: Error not defined. Your code must presents some disfunction'
          CALL flush(6)
          RETURN
       END SELECT
    END IF
    !PRINT *, 'Fact OK'

    !solve matrix
    CALL sfSolve(matrixId,matrixIndependent,matrixSolution)

    !PRINT *, 'Solve OK' 

  END SUBROUTINE Fact_solve_matrix

  SUBROUTINE   Solve_linear(matrixSize,matrixElements,matrixIndependent,matrixSolution)!,ijk)

    INTEGER         , INTENT(IN)                                     :: matrixSize!,ijk
    DOUBLE PRECISION, INTENT(IN)    , DIMENSION(matrixSize,matrixSize) :: matrixElements
    DOUBLE PRECISION, INTENT(IN)    , DIMENSION(matrixSize)            :: matrixIndependent
    DOUBLE PRECISION, INTENT(OUT) , DIMENSION(matrixSize)            :: matrixSolution

    INTEGER :: opt=1
    !if(ijk==1)
    !if(opt == 0) then
    !  CALL Eliminate_0(matrixElements,matrixSize)
    !
    !  CALL Reserve(ijk)
    !else 	

    !clear matrix
    CALL sfZero(matrixId)
    !endif


    ! IN: matrixElements, matrixSize
    CALL Load_matrix(matrixElements,matrixSize)

    ! IN: matrixIndependent, matrixSize
    ! OUT: matrixSolution
    CALL Fact_solve_matrix(matrixIndependent,matrixSolution,matrixSize)

  END SUBROUTINE Solve_linear


END MODULE solve_sparse
