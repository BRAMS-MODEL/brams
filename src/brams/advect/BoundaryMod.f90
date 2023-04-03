MODULE BoundaryMod
   IMPLICIT NONE
   PRIVATE
   
   TYPE bound_data
      INTEGER :: neighbour
      INTEGER :: firstLine
      INTEGER :: lastLine
      INTEGER :: firstColumn
      INTEGER :: lastColumn
      INTEGER :: firstLevel
      INTEGER :: lastLevel
      INTEGER :: itag
      INTEGER :: operation !1= send, 2=recv
      INTEGER :: axys !1 = x , 2 =y
   END TYPE bound_data

   TYPE bound_type
      TYPE(bound_data),POINTER,DIMENSION(:) :: boundLimits
      INTEGER :: lastBoundary
      INTEGER :: m2
      INTEGER :: m3
      INTEGER :: ia
      INTEGER :: ja
      INTEGER :: iz
      INTEGER :: jz
   END TYPE bound_type
   TYPE(bound_type),ALLOCATABLE,DIMENSION(:,:) :: boundary !dimension: ngrid,nprocessors
   
   INTEGER :: maxGrids,maxProcessors
   INTEGER,PARAMETER :: error=0
   INTEGER,PARAMETER :: welldone=-1
   
   PUBLIC createBoundary,destroyBoundary,addBoundary,getfirstLine,getLastLine,getFirstColumn
   PUBLIC getLastColumn,getFirstLevel,getLastLevel,dumpBounds,bound_data,getLastBoundary
   PUBLIC getNeighbour,getSize,getTag,getOperation,getAxys,adjust,getlimits,setLimits
      
   CONTAINS 
   
   INTEGER FUNCTION createBoundary(nGrids,nProcessors)
      INTEGER,INTENT(IN) :: nGrids,nProcessors
      
      INTEGER :: i,j
      
      IF(allocated(boundary)) THEN
         createBoundary=error
      ELSE
         ALLOCATE(boundary(nGrids,nProcessors))
         maxProcessors=nprocessors
	 maxGrids=nGrids
	 DO i=1,nprocessors 
	    DO j=1,nGrids
	       boundary(j,i)%lastBoundary=0
	       IF(associated(boundary(j,i)%boundLimits))    nullify(boundary(j,i)%boundLimits)
	    END DO
	 END DO
         createBoundary=-1
      END IF

   END FUNCTION createBoundary   
   
   INTEGER FUNCTION destroyBoundary()  
      INTEGER :: i,j
   
      IF(.not. allocated(boundary)) THEN
         destroyBoundary=error
      ELSE
	 DO i=1,maxprocessors 
	    DO j=1,maxGrids
	       boundary(j,i)%lastBoundary=0
	       IF(associated(boundary(j,i)%boundLimits))    nullify(boundary(j,i)%boundLimits)
	    END DO
	 END DO
         DEALLOCATE(boundary)
	 destroyBoundary=-1
      END IF   
   END FUNCTION destroyBoundary
   

   SUBROUTINE setLimits(ngrid,np,ia,iz,ja,jz)
      INTEGER,INTENT(IN) :: ngrid,np,ia,iz,ja,jz
      
      boundary(ngrid,np)%ia=ia
      boundary(ngrid,np)%iz=iz
      boundary(ngrid,np)%ja=ja
      boundary(ngrid,np)%jz=jz
            
   END SUBROUTINE setLimits
   
   SUBROUTINE adjust(incx,incy,ngrid,nmachs,m2,m3)
      
      INTEGER, INTENT(IN) :: ngrid,nmachs
      INTEGER,INTENT(IN),DIMENSION(maxProcessors) :: incx,incy,m2,m3
      
      INTEGER :: np,nb,ia,iz,ja,jz
      DO np=1,maxProcessors
         IF(nmachs>1) THEN
	    boundary(ngrid,np)%m2=m2(np)
	    boundary(ngrid,np)%m3=m3(np)
         ELSE
	    boundary(ngrid,np)%m2=boundary(ngrid,np)%iz-boundary(ngrid,np)%ia+3
	    boundary(ngrid,np)%m3=boundary(ngrid,np)%jz-boundary(ngrid,np)%ja+3
	 END IF
     END DO 
      DO np=1,maxProcessors
	  boundary(ngrid,np)%ia=boundary(ngrid,np)%ia-incx(np)
	  boundary(ngrid,np)%iz=boundary(ngrid,np)%iz-incx(np)
	  boundary(ngrid,np)%ja=boundary(ngrid,np)%ja-incy(np)
	  boundary(ngrid,np)%jz=boundary(ngrid,np)%jz-incy(np)
	  DO nb=1,boundary(ngrid,np)%lastBoundary
	     boundary(ngrid,np)%boundLimits(nb)%firstLine=boundary(ngrid,np)%boundLimits(nb)%firstLine-incx(np)
	     boundary(ngrid,np)%boundLimits(nb)%lastLine=boundary(ngrid,np)%boundLimits(nb)%lastLine-incx(np)
	     boundary(ngrid,np)%boundLimits(nb)%firstColumn=boundary(ngrid,np)%boundLimits(nb)%firstColumn-incy(np)
	     boundary(ngrid,np)%boundLimits(nb)%lastColumn=boundary(ngrid,np)%boundLimits(nb)%lastColumn-incy(np)
         END DO
     END DO 
      
   END SUBROUTINE adjust
   
   INTEGER FUNCTION getLAstBoundary(ngrid,nprocessor)
      INTEGER,INTENT(IN) :: ngrid,nProcessor
      IF(.not. allocated(boundary) .or. nProcessor>maxProcessors .or. nGrid>maxGrids) THEN
         getLAstBoundary=0
      ELSE
         getLAstBoundary=boundary(ngrid,nProcessor)%lastBoundary
      END IF
   END FUNCTION getLAstBoundary
   
   
   SUBROUTINE getlimits(ngrid,nprocessor,m2,m3,ia,iz,ja,jz)
      INTEGER,INTENT(IN) :: ngrid,nProcessor
      INTEGER,INTENT(OUT) :: m2,m3,ia,iz,ja,jz
      IF(.not. allocated(boundary) .or. nProcessor>maxProcessors .or. nGrid>maxGrids) THEN
         m2=0
	 m3=0
	 ia=0
	 iz=0
	 ja=0
	 jz=0
      ELSE
         m2=boundary(ngrid,nprocessor)%m2
         m3=boundary(ngrid,nprocessor)%m3
	 ia=boundary(ngrid,nprocessor)%ia
	 iz=boundary(ngrid,nprocessor)%iz
	 ja=boundary(ngrid,nprocessor)%ja
	 jz=boundary(ngrid,nprocessor)%jz
      END IF
      
   END SUBROUTINE getLImits

   INTEGER FUNCTION addBoundary(nGrid,nProcessor,neighbour,localFirstLine,localLastLine, &
                                localFirstColumn,localLastColumn,localFirstLevel,localLastLevel, &
				tag,oper,ax)

      INTEGER,INTENT(IN) :: localFirstLine,localLastLine,neighbour
      INTEGER,INTENT(IN) :: ngrid,nProcessor,localFirstColumn,localLastColumn
      INTEGER,INTENT(IN) :: localFirstLevel,localLastLevel,tag,oper,ax
      
      INTEGER :: i
      
      TYPE(bound_data) :: auxBoundary(boundary(ngrid,nProcessor)%lastBoundary+1)

      
      IF(.not. allocated(boundary) .or. nProcessor>maxProcessors .or. nGrid>maxGrids) THEN
         addBoundary=error
      ELSE
         IF(.not. associated(boundary(ngrid,nProcessor)%boundLimits)) THEN
	    allocate(boundary(ngrid,nProcessor)%boundLimits(1))
	    boundary(ngrid,nProcessor)%boundLimits(1)%neighbour=neighbour
	    boundary(ngrid,nProcessor)%boundLimits(1)%firstLine=localFirstLine
	    boundary(ngrid,nProcessor)%boundLimits(1)%lastLine=localLastLine
	    boundary(ngrid,nProcessor)%boundLimits(1)%firstColumn=localFirstColumn
	    boundary(ngrid,nProcessor)%boundLimits(1)%lastColumn=localLastColumn
	    boundary(ngrid,nProcessor)%boundLimits(1)%firstLevel=localFirstLevel
	    boundary(ngrid,nProcessor)%boundLimits(1)%lastLevel=localLastLevel
	    boundary(ngrid,nProcessor)%boundLimits(1)%itag=tag
	    boundary(ngrid,nProcessor)%boundLimits(1)%operation=oper
	    boundary(ngrid,nProcessor)%boundLimits(1)%axys=ax
	    boundary(ngrid,nProcessor)%lastBoundary=1

	 ELSE
            DO i=1,boundary(ngrid,nProcessor)%lastBoundary
	       auxBoundary(i)= boundary(ngrid,nProcessor)%boundLimits(i)
	    END DO
	    deallocate(boundary(ngrid,nProcessor)%boundLimits)
	    boundary(ngrid,nProcessor)%lastBoundary=boundary(ngrid,nProcessor)%lastBoundary+1
	    allocate(boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary))
            DO i=1,boundary(ngrid,nProcessor)%lastBoundary-1
	      boundary(ngrid,nProcessor)%boundLimits(i)=auxBoundary(i)
	    END DO
	    boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary)%neighbour=neighbour
	    boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary)%firstLine=localFirstLine
	    boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary)%lastLine=localLastLine
	    boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary)%firstColumn=localFirstColumn
	    boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary)%lastColumn=localLastColumn	    
	    boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary)%firstLevel=localFirstLevel
	    boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary)%lastLevel=localLastLevel	    
	    boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary)%itag=tag	    
	    boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary)%operation=oper	
	    boundary(ngrid,nProcessor)%boundLimits(boundary(ngrid,nProcessor)%lastBoundary)%axys=ax	
	 END IF
	 addBoundary=-1
      END IF	      
   END FUNCTION addBoundary
   
   INTEGER FUNCTION getNeighbour(nGrid,nProcessor,nBound)
       
       INTEGER,INTENT(IN) :: nGrid,nProcessor,nBound
       if (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
            nBound>boundary(ngrid,nProcessor)%lastBoundary) THEN
            getNeighbour=0
       ELSE
            getNeighbour=boundary(ngrid,nProcessor)%boundLimits(nBound)%neighbour
       END IF
       
   END FUNCTION getNeighbour
   
   INTEGER FUNCTION getAxys(nGrid,nProcessor,nBound)
       
       INTEGER,INTENT(IN) :: nGrid,nProcessor,nBound
       if (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
            nBound>boundary(ngrid,nProcessor)%lastBoundary) THEN
            getAxys=0
       ELSE
            getAxys=boundary(ngrid,nProcessor)%boundLimits(nBound)%axys
       END IF
       
   END FUNCTION getAxys

   INTEGER FUNCTION getOperation(nGrid,nProcessor,nBound)
       
       INTEGER,INTENT(IN) :: nGrid,nProcessor,nBound
       if (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
            nBound>boundary(ngrid,nProcessor)%lastBoundary) THEN
            getOperation=-1
       ELSE
            getOperation=boundary(ngrid,nProcessor)%boundLimits(nBound)%Operation
       END IF
       
   END FUNCTION getOperation

   INTEGER FUNCTION getTag(nGrid,nProcessor,nBound)
       
       INTEGER,INTENT(IN) :: nGrid,nProcessor,nBound
       if (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
            nBound>boundary(ngrid,nProcessor)%lastBoundary) THEN
            getTag=-1
       ELSE
            getTag=boundary(ngrid,nProcessor)%boundLimits(nBound)%itag
       END IF
       
   END FUNCTION getTag

   INTEGER FUNCTION getfirstLine(nGrid,nProcessor,nBound)
       INTEGER,INTENT(IN) :: nGrid,nProcessor,nBound
       
       if (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
            nBound>boundary(ngrid,nProcessor)%lastBoundary) THEN
	  getFirstLine=-1
       else
          getFirstLine=boundary(ngrid,nProcessor)%boundLimits(nBound)%firstLine
       end if
       
   END FUNCTION getFirstLine

   INTEGER FUNCTION getLastLine(nGrid,nProcessor,nBound)
       INTEGER,INTENT(IN) :: nGrid,nProcessor,nBound
       
       if (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
            nBound>boundary(ngrid,nProcessor)%lastBoundary) THEN
	  getLastLine=-1
       else
          getLastLine=boundary(ngrid,nProcessor)%boundLimits(nBound)%LastLine
       end if
       
   END FUNCTION getLastLine

   INTEGER FUNCTION getFirstColumn(nGrid,nProcessor,nBound)
       INTEGER,INTENT(IN) :: nGrid,nProcessor,nBound
       
       if (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
            nBound>boundary(ngrid,nProcessor)%lastBoundary) THEN
	  getFirstColumn=-1
       else
          getFirstColumn=boundary(ngrid,nProcessor)%boundLimits(nBound)%FirstColumn
       end if
       
   END FUNCTION getFirstColumn

   INTEGER FUNCTION getLastColumn(nGrid,nProcessor,nBound)
       INTEGER,INTENT(IN) :: nGrid,nProcessor,nBound
       
       if (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
            nBound>boundary(ngrid,nProcessor)%lastBoundary) THEN
	  getLastColumn=-1
       else
          getLastColumn=boundary(ngrid,nProcessor)%boundLimits(nBound)%LastColumn
       end if
       
   END FUNCTION getLastColumn

   INTEGER FUNCTION getFirstLevel(nGrid,nProcessor,nBound)
       INTEGER,INTENT(IN) :: nGrid,nProcessor,nBound
       
       if (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
            nBound>boundary(ngrid,nProcessor)%lastBoundary) THEN
	  getFirstLevel=-1
       else
          getFirstLevel=boundary(ngrid,nProcessor)%boundLimits(nBound)%FirstLevel
       end if
       
   END FUNCTION getFirstLevel

   INTEGER FUNCTION getLastLevel(nGrid,nProcessor,nBound)
       INTEGER,INTENT(IN) :: nGrid,nProcessor,nBound
       
       if (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
            nBound>boundary(ngrid,nProcessor)%lastBoundary) THEN
	  getLastLevel=-1
       else
          getLastLevel=boundary(ngrid,nProcessor)%boundLimits(nBound)%LastLevel
       end if
       
   END FUNCTION getLastLevel

   SUBROUTINE dumpBounds()
      INTEGER :: ng,np,nb
      
      DO ng=1,maxGrids
         WRITE (*,FMT='(A,I3.3,A)') '================== Grid ',ng,'=================='; CALL flush(6)
	 WRITE (*,FMT='(8(A3,1X),1X,A5,1X,A8,1X,A1,1X,A1)') 'Npr','Nei','ibg','ien','jbg','jen','zbg','zen','itag ','  Total ','O','D'
         DO np=1,maxProcessors
            DO nb=1,boundary(ng,np)%lastBoundary
	       WRITE (*,FMT='(8(I5,1X),1X,I5.5,1X,I8,1X,I1,1X,I1)') np, &
	                          boundary(ng,np)%boundLimits(nb)%neighbour, &
	                          boundary(ng,np)%boundLimits(nb)%firstLine, &
	                          boundary(ng,np)%boundLimits(nb)%lastLine, &
	                          boundary(ng,np)%boundLimits(nb)%firstColumn, &
	                          boundary(ng,np)%boundLimits(nb)%lastColumn, &	    
	                          boundary(ng,np)%boundLimits(nb)%firstLevel, &
	                          boundary(ng,np)%boundLimits(nb)%lastLevel, &
				  boundary(ng,np)%boundLimits(nb)%itag, &
				  getsize(ng,np,nb),getOperation(ng,np,nb), &
				  getAxys(ng,np,nb) ; CALL flush(6)	    
            END DO
	 END DO
      END DO
      
   END SUBROUTINE dumpBounds

   INTEGER FUNCTION getSize(nGrid,nProcessor,nBound)
       INTEGER,INTENT(IN) :: nGrid,nProcessor,nBound
       
       IF (.not. allocated(boundary) .or. nGrid>maxGrids .or. nProcessor>maxProcessors .or. &
            nBound>boundary(ngrid,nProcessor)%lastBoundary) THEN
	  getSize=0
       ELSE
         getSize=((boundary(ngrid,nProcessor)% boundLimits(nBound)%lastLine - &
	           boundary(ngrid,nProcessor)% boundLimits(nBound)%firstLine +1) * &
		  (boundary(ngrid,nProcessor)% boundLimits(nBound)%lastColumn - &
		  boundary(ngrid,nProcessor)% boundLimits(nBound)%firstColumn +1) * &
		  (boundary(ngrid,nProcessor)% boundLimits(nBound)%lastlevel - &
         	  boundary(ngrid,nProcessor)% boundLimits(nBound)%firstlevel +1))
       END IF
       
   END FUNCTION getSize

END MODULE BoundaryMod
