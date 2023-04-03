MODULE GridMod

   IMPLICIT NONE
   PRIVATE
   TYPE grid_type
      INTEGER :: nx
      INTEGER :: ny
      INTEGER :: nz
      INTEGER, POINTER, DIMENSION(:) :: processors
      INTEGER :: lastProcessor   
   END TYPE grid_type
   TYPE(grid_type),DIMENSION(:),ALLOCATABLE :: grid
   INTEGER :: maxGrids
   INTEGER, PARAMETER :: error=0
   
   PUBLIC createGrid,destroyGrid,setGrid,getNx,getNy,getNz,addProcessor
   PUBLIC getProcessor,getLastProcessor,removeProcessor,findProcessor,getMaxGrids
   
   CONTAINS
   
   INTEGER FUNCTION createGrid(ngrids)
      INTEGER,INTENT(IN) :: ngrids
      
      INTEGER :: i
      
      IF(allocated(grid)) THEN
         createGrid=error
      ELSE
         ALLOCATE(grid(ngrids))
         maxGrids=ngrids
         DO i=1,ngrids
            grid(i)%nx=-1
            grid(i)%ny=-1
            grid(i)%nz=-1
	    grid(i)%lastProcessor=0
	    IF(associated(grid(i)%processors)) nullify(grid(i)%processors)
         END DO	 
         createGrid=-1
      END IF
   END FUNCTION createGrid
   
   INTEGER FUNCTION addProcessor(np,ngrid)
      INTEGER, INTENT(IN) :: np
      INTEGER, INTENT(IN) ::ngrid
      
      INTEGER,ALLOCATABLE,DIMENSION(:) :: auxProcessors
      INTEGER :: i
      
      IF(.not. allocated(grid) .or. ngrid>maxGrids) THEN
         addProcessor=error
      ELSE
         IF(.not. associated(grid(ngrid)%processors)) THEN
	    allocate(grid(ngrid)%processors(1))
	    grid(ngrid)%lastProcessor=grid(ngrid)%lastProcessor+1
	    grid(ngrid)%processors(grid(ngrid)%lastProcessor)=np
	    addProcessor=-1	    
	 ELSE
	    allocate(auxProcessors(grid(ngrid)%lastProcessor))
	    DO i=1,grid(ngrid)%lastProcessor
	       auxProcessors(i)=grid(ngrid)%processors(i)
	    END DO
	    deallocate(grid(ngrid)%processors)
	    grid(ngrid)%lastProcessor=grid(ngrid)%lastProcessor+1
	    ALLOCATE(grid(ngrid)%processors(grid(ngrid)%lastProcessor))
	    DO i=1,grid(ngrid)%lastProcessor-1
	       grid(ngrid)%processors(i)=auxProcessors(i)
	    END DO
	    deallocate(auxProcessors)
	    grid(ngrid)%processors(grid(ngrid)%lastProcessor)=np
	    addProcessor=-1
	 END IF
      END IF
   END FUNCTION addProcessor

   INTEGER FUNCTION removeProcessor(position,ngrid)
      INTEGER, INTENT(IN) :: position,ngrid

      INTEGER,ALLOCATABLE,DIMENSION(:) :: auxProcessors
      INTEGER :: i,k
      
      IF(.not. allocated(grid) .or. .not. associated(grid(ngrid)%processors) .or. &
          position<0 .or. position>grid(ngrid)%lastProcessor .or. ngrid>maxGrids) THEN
         removeProcessor=error
      ELSE
         ALLOCATE(auxProcessors(grid(ngrid)%lastProcessor))
	 DO i=1,grid(ngrid)%lastProcessor
	    auxProcessors(i)=grid(ngrid)%processors(i)
	 END DO
	 deallocate(grid(ngrid)%processors)
	 grid(ngrid)%lastProcessor=grid(ngrid)%lastProcessor-1
	 ALLOCATE(grid(ngrid)%processors(grid(ngrid)%lastProcessor))
	 k=0
	 DO i=1,grid(ngrid)%lastProcessor+1	 
	    IF(i/=position) THEN
	       k=k+1
	       PRINT *, 'Copiando: ',auxProcessors(i),', posicao: ',position
	       grid(ngrid)%processors(k)=auxProcessors(i)
	    END IF
         END DO
         removeProcessor=-1
	 DEALLOCATE(auxProcessors)
      END IF
      
   END FUNCTION removeProcessor
   
   INTEGER FUNCTION findProcessor(np,ngrid)   
      INTEGER,INTENT(IN) :: np,ngrid
      
      INTEGER :: i

      IF(.not. allocated(grid) .or. .not. associated(grid(ngrid)%processors) .or. ngrid>maxGrids) THEN
         findProcessor=-1
      ELSE
         findProcessor=0
	 DO i=1,grid(ngrid)%lastProcessor
	    IF(np==grid(ngrid)%processors(i)) THEN
	       findProcessor=i
	       EXIT
	    END IF
	 END DO
     END IF
   END FUNCTION findProcessor
      
      
   INTEGER FUNCTION getLastProcessor(ngrid)
      INTEGER,INTENT(IN) :: ngrid
      
      IF(.not. allocated(grid)) THEN
         getLastProcessor=-1
      ELSE
         getLastProcessor=grid(ngrid)%lastProcessor
      END IF   
   END FUNCTION getLastProcessor

   INTEGER FUNCTION getProcessor(ngrid,position)
      INTEGER, INTENT(IN) :: ngrid
      INTEGER, INTENT(IN) :: position
      
      IF(position<0 .or. position>grid(ngrid)%lastProcessor .or. .not. allocated(grid) &
         .or. ngrid>maxGrids) THEN
         getProcessor=-1
      ELSE
         getProcessor=grid(ngrid)%processors(position)
      END IF
   END FUNCTION getProcessor
   
   INTEGER FUNCTION destroyGrid()
      
      IF(.not. allocated(grid)) THEN
         destroyGrid=error
      ELSE
         DEALLOCATE(grid)
         maxGrids=0
         destroyGrid=-1
      END IF
      
   END FUNCTION destroyGrid
   
   INTEGER FUNCTION setGrid(ngrid,n1,n2,n3)

      INTEGER,INTENT(IN) :: ngrid,n1,n2,n3
      IF(.not. allocated(grid) .or. ngrid>maxGrids) THEN
         setGrid=error
      ELSE
         grid(ngrid)%nx=n1
         grid(ngrid)%ny=n2
         grid(ngrid)%nz=n3
         setGrid=-1
      END IF
      
   END FUNCTION setGrid

   INTEGER FUNCTION getNx(ngrid)
      INTEGER,INTENT(IN) :: ngrid

      IF(.not. allocated(grid) .or. ngrid>maxGrids) THEN
         getNx=-1
      ELSE
         getNx=grid(ngrid)%nx
      END IF
      
   END FUNCTION  getNx    

   INTEGER FUNCTION getNy(ngrid)
      INTEGER,INTENT(IN) :: ngrid

      IF(.not. allocated(grid) .or. ngrid>maxGrids) THEN
         getNy=-1
      ELSE
         getNy=grid(ngrid)%ny
      END IF
      
   END FUNCTION  getNy    

   INTEGER FUNCTION getNz(ngrid)
      INTEGER,INTENT(IN) :: ngrid

      IF(.not. allocated(grid) .or. ngrid>maxGrids) THEN
         getNz=-1
      ELSE
         getNz=grid(ngrid)%nz
      END IF
      
   END FUNCTION  getNz   
   
   INTEGER FUNCTION getMaxGrids()
   
     getMaxGrids=maxGrids
     
   END FUNCTION getMaxGrids
       

END MODULE GridMod
