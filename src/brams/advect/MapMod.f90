MODULE mapMod

  PRIVATE
  
  TYPE map_type
     INTEGER,POINTER,DIMENSION(:,:) :: processor     
     INTEGER :: maxX
     INTEGER :: maxY
  END TYPE map_type
  TYPE(map_type),DIMENSION(:),ALLOCATABLE :: map
  INTEGER,PARAMETER :: error=0
  INTEGER :: maxMaps
  
  PUBLIC createmap,destroymap,setMap,setMapProcessor,getMapProcessor,getMaxX,getMaxY,dumpMap
  
  CONTAINS
  
  INTEGER FUNCTION createmap(ngrids)
      INTEGER,INTENT(IN) :: ngrids
      
      INTEGER :: i
      
      IF(allocated(map)) THEN
         createmap=error
      ELSE
         ALLOCATE(map(ngrids))
         maxmaps=ngrids
	 DO i=1,ngrids 
	    IF(associated(map(i)%processor)) nullify(map(i)%processor)
	 END DO
         createmap=-1
      END IF
   END FUNCTION createmap

   INTEGER FUNCTION destroymap()
      
      IF(.not. allocated(map)) THEN
         destroymap=error
      ELSE
         DEALLOCATE(map)
         maxmaps=0
         destroymap=-1
      END IF
      
   END FUNCTION destroymap

   INTEGER FUNCTION setMap(ngrid,x,y)
      INTEGER,INTENT(IN) :: ngrid,x,y
      
      IF(.not. allocated(map) .or. associated(map(ngrid)%processor)) THEN
         setMap=error
      ELSE
         allocate(map(ngrid)%processor(x,y))
	 map(ngrid)%maxX=x
	 map(ngrid)%maxY=y
         setMap=-1
	 map(ngrid)%processor=0
      END IF
      
   END FUNCTION setMap
   
   INTEGER FUNCTION getMaxX(ngrid)
      INTEGER, INTENT(IN) :: ngrid
      IF(.not. allocated(map)) THEN
         getMaxX=-1
      ELSE
        getMaxX= map(ngrid)%maxX
      END IF      
   END FUNCTION getMaxX

   INTEGER FUNCTION getMaxY(ngrid)
      INTEGER, INTENT(IN) :: ngrid
      IF(.not. allocated(map)) THEN
         getMaxY=-1
      ELSE
        getMaxY= map(ngrid)%maxY
      END IF      
   END FUNCTION getMaxY

   INTEGER FUNCTION setMapProcessor(ngrid,xpos,ypos,nprocessor)
      INTEGER, INTENT(IN) :: ngrid,xpos,ypos,nprocessor
      
      IF(.not. allocated(map) .or. .not. associated(map(ngrid)%processor) .or. xpos<0 .or. xpos>map(ngrid)%maxX &
          .or. ypos<0 .or. ypos>map(ngrid)%maxY) THEN
         setMapProcessor=error
      ELSE
         map(ngrid)%processor(xpos,ypos)=nprocessor
	 setMapProcessor=-1
      END IF

   END FUNCTION setMapProcessor
   
   INTEGER FUNCTION getMapProcessor(ngrid,xpos,ypos)
      INTEGER :: ngrid,xpos,ypos
      
      IF(.not. allocated(map) .or. .not. associated(map(ngrid)%processor) .or. xpos<0 .or. xpos>map(ngrid)%maxX &
         .or. ypos<0 .or. ypos>map(ngrid)%maxY) THEN
	 getMapProcessor=-1
      ELSE
         getMapProcessor=map(ngrid)%processor(xpos,ypos)
      END IF
   END FUNCTION getMapProcessor
   
   SUBROUTINE dumpMap()
     INTEGER :: nm,i,j
     !- not showing the proc map at log file
     return
     DO nm=1,maxMaps
        WRITE (*,FMT='("---------------------------- Map No ",I3.3," ---------------------------")') nm; CALL flush(6)
	WRITE (*,*) ''
	WRITE (*,FMT='(A,"  ",$)') 'I=    '
	DO i=1,map(nm)%maxX
	   WRITE(*,FMT='(I3.3,A,$)') i,' ';CALL flush(6)
	END DO
	WRITE (*,*) ''
	WRITE (*,*) ''
	DO j=1,map(nm)%maxY
	   WRITE(*,FMT='(A,I3.3,"  ",$)') 'J= ',j
	   DO i=1,map(nm)%maxX
	       WRITE(*,FMT='(I3.3,A,$)') map(nm)%processor(i,j),' ';CALL flush(6)
	   END DO
	   WRITE(*,FMT='(A)') '';CALL flush(6)
        END DO
	WRITE (*,FMT='(A)') '----------------------------------------------------------------------------------'
        CALL flush(6)
    END DO
	      
     
   END SUBROUTINE dumpMap

END MODULE mapMod
