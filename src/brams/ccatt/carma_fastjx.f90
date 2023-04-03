MODULE carma_fastjx
 
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: do3
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: daer
  LOGICAL :: carma_fastjx_allocated
  INTEGER, PARAMETER :: na = 4  ! number of aerosols

  CONTAINS
  
  SUBROUTINE alloc_carma_fastjx(nz,nx,ny)
     
     INTEGER,INTENT(IN) :: nx !max Number of points in i direction
     INTEGER,INTENT(IN) :: ny !max Number of points in j direction
     INTEGER,INTENT(IN) :: nz !max Number of points in k direction 
     IF(carma_fastjx_allocated) THEN
        PRINT *, 'carma_fastjx already allocated'
	STOP
     END IF
     carma_fastjx_allocated=.true.
     ALLOCATE(daer(na,nz,nx,ny));daer=0.0
     
     !-srf-06-05-2008 using harrington dataset 
     !ALLOCATE(do3(nz,nx,ny));do3=0.0
  
  END SUBROUTINE alloc_carma_fastjx
    
END MODULE carma_fastjx
