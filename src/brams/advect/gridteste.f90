PROGRAM testeGrid
   USE GridMod, ONLY : createGrid,destroyGrid,setGrid,getNx,getNy,getNz,addProcessor, &
                       getProcessor,getLastProcessor,removeProcessor,findProcessor
   USE ProcessorMod, ONLY : createProcessor,getNumber,destroyProcessor,setNumber,addGrid, &
                            getGrid,getLastGrid,findGrids, removeGrids
   
   USE BoundaryMod, ONLY: createBoundary,destroyBoundary,addBoundary,getNeighbour, &
                          getfirstLine,getLastLine,getFirstColumn,getLastColumn, &
			  getFirstLevel,getLastLevel,getidTag,getNumberOfPoints, &
			  isToSend
   USE mapMod, ONLY: createmap,destroymap,setMap,setMapProcessor,getMapProcessor,getMaxX,getMaxY
   
   IMPLICIT NONE
   INTEGER :: n_grids,i,j,n_proc,k
   
   n_grids=4
   n_proc=10
   
   IF(createGrid(n_grids)==0) PRINT *, 'Error on create grid'
   !IF(createGrid(n_grids)==0) PRINT *, 'Error on create grid'
   DO i=1,n_grids
      IF(setGrid(i,3+i,4+i,5+i)==0) PRINT *, 'Error on set grid'
   END DO
   DO i=1,6!n_grids
      PRINT *, i,getNx(i),getNy(i),getNz(i)
   END DO  
   PRINT *, 'Testing add';CALL flush(6)
   IF(addProcessor(5,1)==0) PRINT *, 'Error adding processor'
   IF(addProcessor(6,1)==0) PRINT *, 'Error adding processor'
   PRINT *, 'List of processors'
   DO i=1,n_grids
      DO j=1,getLastProcessor(i)
         PRINT *,i,j,getProcessor(i,j)
      END DO
   END DO
   IF(removeProcessor(1,1)==0) PRINT *, 'Error removing processor'
   PRINT *, 'New List of processors'
   DO i=1,n_grids
      DO j=1,getLastProcessor(i)
         PRINT *,i,j,getProcessor(i,j)
      END DO
   END DO
   IF(addProcessor(1,1)==0) PRINT *, 'Error adding processor'
   IF(addProcessor(2,1)==0) PRINT *, 'Error adding processor'
   IF(addProcessor(8,1)==0) PRINT *, 'Error adding processor'
   IF(addProcessor(3,1)==0) PRINT *, 'Error adding processor'
   PRINT *, 'Newest List of processors'
   DO i=1,n_grids
      DO j=1,getLastProcessor(i)
         PRINT *,i,j,getProcessor(i,j)
      END DO
   END DO
   PRINT *,findProcessor(8,1)
   PRINT *, 'Creating Processor'
   IF(createProcessor(n_proc)==0) PRINT *, 'Error on create processor'
   PRINT *, 'Setting Processor number'
   DO i=1,n_proc
      IF(setNumber(i,i)==0) PRINT *, 'Error on setting processor'
   END DO
   PRINT *, 'getting Processor number'
   DO i=1,n_proc
      PRINT *,i,getNumber(i)
   END DO
   PRINT *, 'Adding Processor grids'
   IF(addgrid(1,1)==0) PRINT *, 'Error adding grid'
   IF(addgrid(2,1)==0) PRINT *, 'Error adding grid'
   IF(addgrid(8,1)==0) PRINT *, 'Error adding grid'
   IF(addgrid(3,1)==0) PRINT *, 'Error adding grid'
   PRINT *, 'List of grids'
   DO i=1,n_proc
      DO j=1,getLastGrid(i)
         PRINT *,i,j,getGrid(i,j)
      END DO
   END DO
   PRINT *,findGrids(8,1),findGrids(9,1)
   IF(removeGrids(2,1)==0) PRINT *, 'Error removing grid'
   PRINT *, 'New List of grids'
   DO i=1,n_proc
      DO j=1,getLastGrid(i)
         PRINT *,i,j,getGrid(i,j)
      END DO
   END DO
   IF(createBoundary(n_grids,n_proc)==0) PRINT *,'Error on create boundary'
   IF (addBoundary(nGrid=1,nProcessor=1,localNeighbour=2,localFirstLine=1,localLastLine=10, &
                                localFirstColumn=9,localLastColumn=10,localFirstLevel=1, &
				localLastLevel=33,localToSend=.true.)==0) &
				PRINT *, 'Error adding boundary'
   IF (addBoundary(nGrid=1,nProcessor=1,localNeighbour=2,localFirstLine=1,localLastLine=10, &
                                localFirstColumn=9,localLastColumn=10,localFirstLevel=1, &
				localLastLevel=33,localToSend=.false.)==0) &
				PRINT *, 'Error adding boundary'
   IF (addBoundary(nGrid=1,nProcessor=1,localNeighbour=3,localFirstLine=1,localLastLine=10, &
                                localFirstColumn=9,localLastColumn=10, &
				localFirstLevel=2,localLastLevel=31,localToSend=.true.)==0) &
				PRINT *, 'Error adding boundary'
   IF (addBoundary(nGrid=1,nProcessor=1,localNeighbour=3,localFirstLine=1,localLastLine=10, &
                                localFirstColumn=9,localLastColumn=10,localFirstLevel=2 &
				,localLastLevel=31,localToSend=.false.)==0) &
				PRINT *, 'Error adding boundary'

   DO i=1,1
      DO j=1,1
         DO k=1,4
	    WRITE (*,FMT='(3(I2.2,1X),7(I3.3,1X),I6.6,1X,I8.8,1X,L1)') i,j,k,getNeighbour(nGrid=j,nProcessor=i,nBound=k), &
	                  getFirstLine(nGrid=j,nProcessor=i,nBound=k), &
	                  getLastLine(nGrid=j,nProcessor=i,nBound=k), &
	                  getFirstColumn(nGrid=j,nProcessor=i,nBound=k), &
	                  getLastColumn(nGrid=j,nProcessor=i,nBound=k), &
	                  getFirstLevel(nGrid=j,nProcessor=i,nBound=k), &
	                  getLastLevel(nGrid=j,nProcessor=i,nBound=k), &
	                  getIdTag(nGrid=j,nProcessor=i,nBound=k), &
			  getNumberOfPoints(nGrid=j,nProcessor=i,nBound=k), &
			  isToSend(nGrid=j,nProcessor=i,nBound=k)
	 END DO
      END DO
   END DO
   
   IF(createMap(n_grids)==0) PRINT *, 'Error on create map'
   DO i=1,n_grids
      IF(setMap(i,10+i,15+i)==0) PRINT *, 'Error on set map'
   END DO
   DO i=1,n_grids
      DO j=1,getMaxX(i)
         DO k=1,getMaxY(i)
            IF(setMapProcessor(i,j,k,1)==0) PRINT *, 'Error on set mapprocessor'
	 END DO
      END DO
   END DO
   DO i=1,1
      DO j=1,getMaxX(i)
         DO k=1,getMaxY(i)
            PRINT *, i,j,k,getMapProcessor(i,j,k),getNeighbour(i,getMapProcessor(i,j,k),4)
	 END DO
      END DO
  END DO
   
   
END PROGRAM testeGrid
