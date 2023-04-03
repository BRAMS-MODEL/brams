MODULE ProcessorMod
   IMPLICIT NONE
   
   PRIVATE
   
   TYPE processor_Type
      INTEGER :: number
      INTEGER,POINTER,DIMENSION(:) :: grids
      INTEGER :: lastGrid
      !North=1,south=2,east=3,west=4
      INTEGER,DIMENSION(4) :: neighbour !
   END TYPE processor_Type
   TYPE(processor_Type),DIMENSION(:,:),ALLOCATABLE :: processor
   INTEGER :: maxProcessors
   INTEGER,PARAMETER :: pNorth=1,pSouth=2,pEast=3,pWest=4
   INTEGER,PARAMETER :: error=0
   
   PUBLIC createProcessor,getNumber,destroyProcessor,setNumber!,addGrid,getGrid,getLastGrid
   PUBLIC pNorth,pSouth,pEast,pWest,setNeighbour,getPNeighbour!findGrids, removeGrids,
   PUBLIC getMaxProcessors,dumpNeighbour
   
   CONTAINS
   
   INTEGER FUNCTION createProcessor(nprocessors, ngrids)
      INTEGER,INTENT(IN) :: nprocessors, ngrids
      
      INTEGER :: i,j, n
      
      IF(allocated(processor)) THEN
         createProcessor=error
      ELSE
         ALLOCATE(processor(nprocessors, ngrids))
         maxProcessors=nprocessors
	DO n = 1, ngrids
	 DO i=1,nprocessors 
	    processor(i,n)%number=-1
	    processor(i,n)%lastGrid=0
	    DO j=1,4
	        processor(i,n)%neighbour(j)=0
	    END DO
	    IF(associated(processor(i,n)%grids)) nullify(processor(i,n)%grids)
	 END DO
	END DO
         createProcessor=-1
      END IF
   END FUNCTION createProcessor

   INTEGER FUNCTION destroyProcessor()
      
      IF(.not. allocated(Processor)) THEN
         destroyProcessor=error
      ELSE
         DEALLOCATE(Processor)
         maxProcessors=0
         destroyProcessor=-1
      END IF
      
   END FUNCTION destroyProcessor

   INTEGER FUNCTION setNeighbour(nProcessor,neig,direction,ngrid)
      INTEGER, INTENT(IN) :: nProcessor,neig,direction,ngrid
      
      IF(.not. allocated(processor) .or. nProcessor>maxProcessors .or. neig>maxProcessors .or. direction<1 .or. direction>4) THEN
         setNeighbour=error
      ELSE
         processor(nProcessor,ngrid)%neighbour(direction)=neig
	 setNeighbour=-1
      END IF
   END FUNCTION setNeighbour
   
   INTEGER FUNCTION getPNeighbour(nProcessor,direction,ngrid)
      INTEGER, INTENT(IN) :: nProcessor,direction,ngrid
      
      IF(.not. allocated(processor) .or. nProcessor>maxProcessors .or. direction<pNorth .or. direction>pWest) THEN
         getPNeighbour=-1
      ELSE
         getPNeighbour=processor(nProcessor,ngrid)%neighbour(direction)
      END IF
   END FUNCTION getPNeighbour
   
   INTEGER FUNCTION getNumber(nprocessor,ngrid)
      INTEGER,INTENT(IN) :: nprocessor,ngrid
      
      IF(.not. allocated(processor) .or. nProcessor>maxProcessors) THEN
        getNumber=-1
      ELSE
        getNumber=processor(nprocessor,ngrid)%number
      END IF
   
   END FUNCTION getNumber

   INTEGER FUNCTION setNumber(nprocessor,num,ngrid)
      INTEGER,INTENT(IN) :: nprocessor,num,ngrid
      
      IF(.not. allocated(processor) .or. nProcessor>maxProcessors) THEN
        setNumber=error
      ELSE
        processor(nprocessor,ngrid)%number=num
	setNumber=-1
      END IF
   
   END FUNCTION setNumber
   
!LFR>  INTEGER FUNCTION addGrid(ngrid,nProcessor)
!LFR>      INTEGER, INTENT(IN) :: nProcessor
!LFR>      INTEGER, INTENT(IN) ::ngrid
!LFR>      
!LFR>      INTEGER,ALLOCATABLE,DIMENSION(:) :: auxGrids !!!!!!!!
!LFR>      INTEGER :: i
!LFR>     
!LFR>      IF(.not. allocated(processor) .or. nProcessor>maxProcessors) THEN
!LFR>  	addGrid=error
!LFR>      ELSE
!LFR>  	IF(.not. associated(processor(nProcessor)%grids)) THEN
!LFR>  	   allocate(processor(nProcessor)%grids(1))
!LFR>  	   processor(nProcessor)%lastGrid=processor(nProcessor)%lastGrid+1
!LFR>  	   processor(nProcessor)%grids(processor(nProcessor)%lastGrid)=ngrid
!LFR>  	   addGrid=-1	   
!LFR>  	ELSE
!LFR>  	   allocate(auxGrids(processor(nProcessor)%lastGrid))
!LFR>  	   DO i=1,processor(nProcessor)%lastGrid
!LFR>  		      auxgrids(i)=processor(nProcessor)%grids(i)
!LFR>  	   END DO
!LFR>  	   deallocate(processor(nProcessor)%grids)
!LFR>  	   processor(nProcessor)%lastGrid=processor(nProcessor)%lastGrid+1
!LFR>  	   ALLOCATE(processor(nProcessor)%grids(processor(nProcessor)%lastGrid))
!LFR>  	   DO i=1,processor(nProcessor)%lastGrid-1
!LFR>  	      processor(nProcessor)%grids(i)=auxgrids(i)
!LFR>  	   END DO
!LFR>  	   deallocate(auxGrids)
!LFR>  	   processor(nProcessor)%grids(processor(nProcessor)%lastGrid)=nGrid
!LFR>  	   addGrid=-1
!LFR>  	END IF
!LFR>      END IF
!LFR>   END FUNCTION addGrid
!LFR> 
!LFR>   INTEGER FUNCTION getLastGrid(nProcessor)
!LFR>      INTEGER,INTENT(IN) :: nProcessor
!LFR>      
!LFR>      IF(.not. allocated(processor) .or. nProcessor<0 .or. nProcessor>maxProcessors) THEN
!LFR>  	getLastGrid=-1
!LFR>      ELSE
!LFR>  	getLastGrid=processor(nProcessor)%lastGrid
!LFR>      END IF   
!LFR>   END FUNCTION getLastGrid
!LFR> 
!LFR>   INTEGER FUNCTION getGrid(nprocessor,position)
!LFR>      INTEGER, INTENT(IN) :: nprocessor
!LFR>      INTEGER, INTENT(IN) :: position
!LFR>      
!LFR>      IF(position<0 .or. position>processor(nprocessor)%lastGrid .or. .not. allocated(processor) &
!LFR>  	.or. nprocessor>maxProcessors) THEN
!LFR>  	getGrid=-1
!LFR>      ELSE
!LFR>  	getGrid=processor(nprocessor)%grids(position)
!LFR>      END IF
!LFR>   END FUNCTION getGrid
!LFR>      
!LFR>   INTEGER FUNCTION findGrids(ng,nprocessor)   
!LFR>      INTEGER,INTENT(IN) :: ng,nProcessor
!LFR>      
!LFR>      INTEGER :: i
!LFR> 
!LFR>      IF(.not. allocated(processor) .or. .not. associated(processor(nProcessor)%grids) .or. &
!LFR>  	 nProcessor>maxProcessors) THEN
!LFR>  	findGrids=-1
!LFR>      ELSE
!LFR>  	findGrids=0
!LFR>  	DO i=1,processor(nProcessor)%lastGrid
!LFR>  	   IF(ng==processor(nProcessor)%grids(i)) THEN
!LFR>  	      findGrids=i
!LFR>  	      EXIT
!LFR>  	   END IF
!LFR>  	END DO
!LFR>     END IF
!LFR>   END FUNCTION findGrids
!LFR> 
!LFR>   INTEGER FUNCTION removeGrids(position,nProcessor)
!LFR>      INTEGER, INTENT(IN) :: position,nProcessor
!LFR> 
!LFR>      INTEGER,ALLOCATABLE,DIMENSION(:) :: auxgrids
!LFR>      INTEGER :: i,k
!LFR>      
!LFR>      IF(.not. allocated(processor) .or. .not. associated(processor(nprocessor)%grids) .or. &
!LFR>  	 position<0 .or. position>processor(nprocessor)%lastgrid .or. nprocessor>maxProcessors) THEN
!LFR>  	 removegrids=error
!LFR>      ELSE
!LFR>  	ALLOCATE(auxgrids(processor(nprocessor)%lastgrid))
!LFR>  	DO i=1,processor(nprocessor)%lastgrid
!LFR>  	   auxgrids(i)=processor(nprocessor)%grids(i)
!LFR>  	END DO
!LFR>  	deallocate(processor(nprocessor)%grids)
!LFR>  	processor(nprocessor)%lastgrid=processor(nprocessor)%lastgrid-1
!LFR>  	ALLOCATE(processor(nprocessor)%grids(processor(nprocessor)%lastgrid))
!LFR>  	k=0
!LFR>  	DO i=1,processor(nprocessor)%lastgrid+1 	
!LFR>  	   IF(i/=position) THEN
!LFR>  	      k=k+1
!LFR>  	      processor(nprocessor)%grids(k)=auxgrids(i)
!LFR>  	   END IF
!LFR>  	END DO
!LFR>  	removegrids=-1
!LFR>  	DEALLOCATE(auxgrids)
!LFR>      END IF
!LFR>      
!LFR>   END FUNCTION removeGrids

   INTEGER FUNCTION getMaxProcessors()
     
     getMaxProcessors=maxProcessors
   
   END FUNCTION getMaxProcessors
   
   SUBROUTINE dumpNeighbour(ngrids)
      INTEGER, INTENT(IN) :: ngrids
     
      INTEGER :: p, n
      
      DO n = 1, ngrids
 	     WRITE (*,FMT='(A)') '=========================================='
      	WRITE (*,FMT='(A5,1X,5(A3,1X))') 'Proc.',' N ',' S ',' E ',' W '
      	DO p=1,maxProcessors
          WRITE (*,FMT='(I5.5,1X,5(I3.3,1X))') p,processor(p,n)%neighbour(1), &
	               processor(p,n)%neighbour(2), &
	               processor(p,n)%neighbour(3), &
	               processor(p,n)%neighbour(4)
     	 END DO
      END DO	         
      
   END SUBROUTINE dumpNeighbour
   
END MODULE ProcessorMod
