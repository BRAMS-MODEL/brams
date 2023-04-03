MODULE AdvectData
   IMPLICIT NONE
   PRIVATE
   PUBLIC InitAdvect
   !PUBLIC showMessages
   
  CONTAINS

  SUBROUTINE showNeighbour(ngrids)
     USE ProcessorMod, ONLY: dumpNeighbour
     INTEGER, INTENT(IN) :: ngrids
     CALL dumpNeighbour(ngrids)
  END SUBROUTINE showNeighbour

  SUBROUTINE showMap()
     USE MapMod, ONLY: dumpMap
     CALL dumpMap()
  END SUBROUTINE showMap
  
  SUBROUTINE showBounds()
     USE BoundaryMod, ONLY: dumpBounds
     CALL dumpBounds()
  END SUBROUTINE showBounds

  SUBROUTINE showMessages(fout,ngrids)
     USE advMessageMod, ONLY: dumpMessages
     INTEGER, INTENT(IN) :: fout,ngrids
     CALL dumpMessages(fout,ngrids)
  END SUBROUTINE showMessages

  SUBROUTINE InitGrid(ngrids,nnxp,nnyp,nnzp)
     USE GridMod, ONLY: createGrid,setGrid
     USE errorMod, ONLY: printError
     
     INTEGER, INTENT(IN) :: ngrids
     INTEGER, INTENT(IN), DIMENSION(ngrids) :: nnxp,nnyp,nnzp
     INTEGER :: i,istat
     
     IF(creategrid(ngrids)==0) istat=printError(1,'Creating grid on Initialization',6)
     DO i=1,ngrids
       IF(setGrid(i,nnxp(i),nnyp(i),nnzp(i))==0) istat=printError(1,'Setting grid on Initialization',6)
     END DO
     
  END SUBROUTINE InitGrid
  
  SUBROUTINE InitAdvect(ngrids,nmachs,mynum,ghostzone,nnxp,nnyp,nnzp,ixb,ixe,iyb,iye)
     USE advMessageMod, ONLY: dumpMessages
     
     INTEGER, INTENT(IN) :: ngrids,nmachs,mynum,ghostZone
     INTEGER, INTENT(IN), OPTIONAL, DIMENSION(ngrids) :: nnxp,nnyp,nnzp
     INTEGER, INTENT(IN), OPTIONAL, DIMENSION(nmachs,ngrids) :: ixb,ixe,iyb,iye

    IF (mynum==1) THEN
       CALL InitGrid(ngrids,nnxp,nnyp,nnzp)
       CALL InitMap(nmachs,ngrids,ixb,ixe,iyb,iye,nnxp,nnyp,nnzp)
       CALL showMap()
       CALL InitBounds(ixb,ixe,iyb,iye,nmachs,ngrids,ghostZone)
       !CALL showBounds()
       !CALL showNeighbour(ngrids)
    END IF
    CALL sendInit(nmachs,mynum,ngrids)
    !IF(myNum/=0) CALL dumpMessages(70+mynum,ngrids)   
  END SUBROUTINE InitAdvect
  
  SUBROUTINE InitMap(nmachs,ngrids,ixb,ixe,iyb,iye,nnxp,nnyp,nnzp)
     USE MapMod, ONLY: createmap,setMap,setMapProcessor
     USE ProcessorMod, ONLY: createProcessor
     USE GridMod, ONLY: addProcessor,findProcessor
     USE errorMod, ONLY: printError
     
     INTEGER, INTENT(IN) :: nmachs,ngrids
     INTEGER, INTENT(IN), DIMENSION(nmachs,ngrids) :: ixb,ixe,iyb,iye
     INTEGER, INTENT(IN), DIMENSION(ngrids) :: nnxp,nnyp,nnzp
     INTEGER :: i,j,ngrid,nproc,istat
     
     IF(createProcessor(nmachs,ngrids)==0) istat=printError(1,'Creating processor on Initialization',0)
     IF(createmap(ngrids)==0) istat=printError(1,'Creating map on Initialization',0)
     DO i=1,ngrids
  	IF(setMap(i,nnxp(i),nnyp(i))==0) istat=printError(1,'Setting map on Initialization',0)
     END DO
     DO ngrid=1,ngrids
  	DO nproc=1,nmachs
  	   DO i=ixb(nproc,ngrid),ixe(nproc,ngrid)
  	      DO j=iyb(nproc,ngrid),iye(nproc,ngrid)
  		 !Putting processor in especific point of grid
 		 IF(setMapProcessor(ngrid,i,j,nproc)==0) istat=printError(1,'Attributing processor to map on Initialization',6)
  		 !Adding  a grid to processor
  		 !IF(findGrids(ngrid,nproc)==0) istat=addGrid(ngrid,nproc)
  		 !Adding processor to grid
  		 !IF(findProcessor(nproc,ngrid)==0) istat=addProcessor(nproc,ngrid)  
  	      END DO
  	   END DO
  	END DO
     END DO
     
  END SUBROUTINE InitMap
  
  SUBROUTINE InitBounds(ixb,ixe,iyb,iye,nmachs,ngrids,ghostZone)
     USE GridMod, ONLY: getMaxGrids,getNx,getNy,getNz
     USE MapMod, ONLY: getMapProcessor
     USE ProcessorMod, ONLY: pNorth,pSouth,pEast,pWest,setNeighbour,getPNeighbour,getMaxProcessors
     USE BoundaryMod, ONLY: createBoundary,destroyBoundary,addBoundary,adjust,setLimits,dumpBounds
     USE errorMod, ONLY: printError
     use dump

     include "constants.f90"

     INTEGER, INTENT(IN) :: nmachs,ngrids,ghostZone
     INTEGER, INTENT(IN), DIMENSION(nmachs,ngrids) :: ixb,ixe,iyb,iye
     INTEGER :: ngrid,i,istep,j,ini,fin,istat,p1,p2,ii,nextJ,jj,jl,il,yini,yfin,xini,xfin,gz,gzi
     INTEGER :: interval,err
     LOGICAL,ALLOCATABLE,DIMENSION(:) :: foundp
     INTEGER :: tagCount
     INTEGER :: incx(nmachs),incy(nmachs)
     INTEGER :: m2(nmachs),m3(nmachs)
     LOGICAL :: ThereAreNeighbour(4,nmachs)
     
     INTEGER,PARAMETER :: pSend=1,pRecv=2,we=1,ew=2,ns=3,sn=4
 
  !   +------------------------------------------------------------------------+
  !   | 	   N					  ^		       |
  !   | 	   |					  |		       |
  !   |        W <-.-> E   ------> x direction = i	  | y direction = j    |
  !   | 	   |					  |		       |
  !   | 	   S					  |		       |
  !   +------------------------------------------------------------------------+
  
  
     IF(createBoundary(getMaxGrids(),getMaxProcessors())==0) istat=printError(1,'Creating processor on Initialization',6)
     ALLOCATE(foundp(getMaxProcessors()))
     foundp=.false.
     tagCount=0
     ThereAreNeighbour=.false.
     gz=ghostZone
     gzi=gz
     DO ngrid=1,ngrids
        DO p1=1,nmachs
	   interval=ixe(p1,ngrid)-ixb(p1,ngrid)+2
	   gz=min(interval,gz)
	   interval=iye(p1,ngrid)-iyb(p1,ngrid)+2
	   gz=min(interval,gz)	   	   
	END DO
     END DO
     if(gzi/=gz) err=dumpMessage(c_tty,c_yes,'InitBounds','ADVC_MNT',c_notice,'Ghostzone required is different of ghostzone possible:',(/gzi,gz/),'I2')
     DO ngrid=1,ngrids
       DO p1=1,nmachs
          !write (88,fmt='(A,I4.4,1X,A,4(I4.4,1X))') '----- ',p1,'-----',ixb(p1,ngrid),ixe(p1,ngrid),iyb(p1,ngrid),iye(p1,ngrid)
	  !write (88,fmt='(2(A4,1X))') '  i ',' nb '
          do i=ixb(p1,ngrid),ixe(p1,ngrid)
	     p2=getMapProcessor(ngrid,i,iye(p1,ngrid)+1)
	     !write (88,fmt='(2(I4.4,1X))') i,p2
	     IF(p2/=0) THEN 
	        ThereAreNeighbour(pSouth,p1)=.true.
	        ThereAreNeighbour(pNorth,p2)=.true.
	     END IF
	  end do
	  !write (88,fmt='(2(A4,1X))') '  j ',' nb '
          do j=iyb(p1,ngrid),iye(p1,ngrid)
	     p2=getMapProcessor(ngrid,ixe(p1,ngrid)+1,j)
	     !write (88,fmt='(2(I4.4,1X))') j,p2
	     IF(p2/=0) THEN 
	        ThereAreNeighbour(pEast,p1)=.true.
	        ThereAreNeighbour(pWest,p2)=.true.
	     END IF
	  end do	  
       END DO 
       !write (88,fmt='(A3,1X,4(A1,1X))') 'Proc','N','S','L','W'
       !DO p1=1,nmachs
       !    write (88,fmt='(I3.3,1X,4(L1,1X))') p1,(ThereAreNeighbour(i,p1),i=1,4)
       !END DO
       !CALL flush(88)
           
       DO p1=1,nmachs
	   incx(p1)=max(ixb(p1,ngrid)-gz-1,0)
	   incy(p1)=max(iyb(p1,ngrid)-gz-1,0)

	   p2=getMapProcessor(ngrid,ixe(p1,ngrid)+1,iyb(p1,ngrid))
	   !IF(p1==275 .or. p1==259) WRITE (*,fmt='("B1: ",4(I3.3,1X))') p1,p2,ixe(p1,ngrid)+1,iyb(p1,ngrid); call flush(6)
	   IF(p2/=0) THEN
	      istat=setNeighbour(p1,p2,pEast,ngrid)
	      istat=setNeighbour(p2,p1,pWest,ngrid)	    
           END IF

	   p2=getMapProcessor(ngrid,ixb(p1,ngrid),iye(p1,ngrid)+1)
	   !IF(p1==275 .or. p1==259) WRITE (*,fmt='("B2: ",4(I3.3,1X))') p1,p2,ixb(p1,ngrid),iye(p1,ngrid)+1; call flush(6)
	   IF(p2/=0) THEN
	      istat=setNeighbour(p1,p2,pSouth,ngrid)
	      istat=setNeighbour(p2,p1,pNorth,ngrid)
	   END IF
      END DO

      DO p1=1,nmachs
!	   WRITE (*,FMT='("Debug Neig: ",7(I3.3,1X))') p1,getPNeighbour(p1,1),getPNeighbour(p1,2), &
!	                         getPNeighbour(p1,3),getPNeighbour(p1,4),incx(p1),incy(p1); CALL flush(6)
           !right (east) and left (west) boundaries
!         IF(getPNeighbour(p1,1,ngrid)/=0) THEN
         IF(ThereAreNeighbour(1,p1)) THEN
              yini=iyb(p1,ngrid)-gz
	 ELSE
	      yini=iyb(p1,ngrid)-1
	 END IF
	   
         IF(ThereAreNeighbour(2,p1)) THEN
!         IF(getPNeighbour(p1,2,ngrid)/=0) THEN
              yfin=iye(p1,ngrid)+gz
	 ELSE
	      yfin=iye(p1,ngrid)+1
	 END IF

         IF(ThereAreNeighbour(4,p1)) THEN
!         IF(getPNeighbour(p1,4,ngrid)/=0) THEN
              xini=ixb(p1,ngrid)-gz
	 ELSE
	      xini=ixb(p1,ngrid)-1
	 END IF
	   
         IF(ThereAreNeighbour(3,p1)) THEN
         !IF(getPNeighbour(p1,3,ngrid)/=0) THEN
              xfin=ixe(p1,ngrid)+gz
	 ELSE
	      xfin=ixe(p1,ngrid)+1
	 END IF

	   m2(p1)=xfin-xini+1
	   m3(p1)=yfin-yini+1
!LFR> 	   IF(p1==275 .or. p1==274) WRITE (*,FMT='("1.Proc,M2,M3: ",7(I4.4,1X))') &
!LFR> 	          p1,m2(p1),m3(p1),getPNeighbour(p1,1,ngrid),getPNeighbour(p1,2,ngrid),yini,yfin; call flush(6)
	END DO

    DO p1=1,nmachs
           !right (east) and left (west) boundaries
           IF(getPNeighbour(p1,1,ngrid)/=0) THEN
              yini=iyb(p1,ngrid)-gz
	   ELSE
	      yini=iyb(p1,ngrid)-1
	   END IF
	   !! if(yini .eq. -1)then ; print*, 'y', p1, yini, iyb(p1,ngrid), getPNeighbour(p1,1) ; stop ; endif!remove
	   p2=getMapProcessor(ngrid,ixe(p1,ngrid)+1,iyb(p1,ngrid))
	   IF(p2/=0) THEN
	      IF(iye(p2,ngrid)>iye(p1,ngrid)) THEN
                 IF(getPNeighbour(p2,2,ngrid)/=0) THEN
	            yfin=iye(p2,ngrid)+gz
		 ELSE
		    yfin=iye(p2,ngrid)+1
		 END IF
                 !P1->P2
	         tagCount=tagCount+1
	         istat=addBoundary(nGrid,p1,p2,ixe(p1,ngrid)+1-gz,ixe(p1,ngrid), &
		       yini,yfin,1,getNz(ngrid),tagCount,pSend,we)
                 istat=addBoundary(nGrid,p2,p1,ixe(p1,ngrid)+1-gz,ixe(p1,ngrid), &
		       yini,yfin,1,getNz(ngrid),tagCount,pRecv,we)		       
!	         WRITE (*,FMT='("a->",7(I2.2,1X))') ngrid,p1,p2,ixe(p1,ngrid)-gz,ixe(p1,ngrid),yini,yfin; CALL flush(6)
!	         WRITE (*,FMT='("a<-",7(I2.2,1X))') ngrid,p2,p1,ixe(p1,ngrid)-gz,ixe(p1,ngrid),yini,yfin; CALL flush(6)
                 !P2->P1
	         tagCount=tagCount+1
                 istat=addBoundary(nGrid,p2,p1,ixe(p1,ngrid)+1,ixe(p1,ngrid)+gz, &
		       yini,yfin,1,getNz(ngrid),tagCount,pSend,ew)
	         istat=addBoundary(nGrid,p1,p2,ixe(p1,ngrid)+1,ixe(p1,ngrid)+gz, &
		       yini,yfin,1,getNz(ngrid),tagCount,pRecv,ew)		       
!	         WRITE (*,FMT='("b->",7(I2.2,1X))') ngrid,p1,p2,ixe(p1,ngrid)+1,ixe(p1,ngrid)+gz,yini,yfin; CALL flush(6)
!	         WRITE (*,FMT='("b<-",7(I2.2,1X))') ngrid,p2,p1,ixe(p1,ngrid)+1,ixe(p1,ngrid)+gz,yini,yfin; CALL flush(6)
	      ELSE
	         tagCount=tagCount+1
                 IF(getPNeighbour(p1,2,ngrid)/=0) THEN
	            yfin=iye(p1,ngrid)+gz
                 ELSE
		    yfin=iye(p1,ngrid)+1
		 END IF
                 istat=addBoundary(nGrid,p1,p2,ixe(p1,ngrid)+1-gz,ixe(p1,ngrid), &
		       yini,yfin,1,getNz(ngrid),tagCount,pSend,we)
                 istat=addBoundary(nGrid,p2,p1,ixe(p1,ngrid)+1-gz,ixe(p1,ngrid), &
		       yini,yfin,1,getNz(ngrid),tagCount,pRecv,we)
!	         WRITE (*,FMT='("c->",7(I2.2,1X))') ngrid,p1,p2,ixe(p1,ngrid)+1-gz,ixe(p1,ngrid),yini,yfin; CALL flush(6)
!	         WRITE (*,FMT='("c<-",7(I2.2,1X))') ngrid,p2,p1,ixe(p1,ngrid)+1-gz,ixe(p1,ngrid),yini,yfin; CALL flush(6)
		 tagCount=tagCount+1     
                 istat=addBoundary(nGrid,p2,p1,ixe(p1,ngrid)+1,ixe(p1,ngrid)+gz, &
		       yini,yfin,1,getNz(ngrid),tagCount,pSend,ew)
                 istat=addBoundary(nGrid,p1,p2,ixe(p1,ngrid)+1,ixe(p1,ngrid)+gz, &
		       yini,yfin,1,getNz(ngrid),tagCount,pRecv,ew)
!	         WRITE (*,FMT='("d->",7(I2.2,1X))') ngrid,p1,p2,ixe(p1,ngrid)+1,ixe(p1,ngrid)+gz,yini,yfin; CALL flush(6)
!	         WRITE (*,FMT='("d<-",7(I2.2,1X))') ngrid,p2,p1,ixe(p1,ngrid)+1,ixe(p1,ngrid)+gz,yini,yfin; CALL flush(6)
	      END IF
           END IF

           !upper (north) and lower (south) boundaries
	    IF(getPNeighbour(p1,4,ngrid)/=0) THEN
	      xini=ixb(p1,ngrid)-gz
	   ELSE
	      xini=ixb(p1,ngrid)-1
	   END IF  
	   !! if(xini .eq. -1)then ; print*, 'x', p1, xini, ixb(p1,ngrid), getPNeighbour(p1,4) ; stop ; endif!remove
	   p2=getMapProcessor(ngrid,ixb(p1,ngrid),iye(p1,ngrid)+1)
	   IF(p2/=0) THEN
	      IF(ixe(p2,ngrid)<ixe(p1,ngrid)) THEN
	        tagCount=tagCount+1
                IF(getPNeighbour(p2,3,ngrid)/=0) THEN
	           xfin=ixe(p2,ngrid)+gz
		ELSE
		   xfin=ixe(p2,ngrid)+1
		END IF
                istat=addBoundary(nGrid,p1,p2,xini,xfin,iye(p1,ngrid)+1-gz,iye(p1,ngrid), &
		      1,getNz(ngrid),tagCount,pSend,ns)
                istat=addBoundary(nGrid,p2,p1,xini,xfin,iye(p1,ngrid)+1-gz, &
		      iye(p1,ngrid),1,getNz(ngrid),tagCount,pRecv,ns)
		tagCount=tagCount+1
                istat=addBoundary(nGrid,p2,p1,xini,xfin,iye(p1,ngrid)+1, &
		      iye(p1,ngrid)+gz,1,getNz(ngrid),tagCount,pSend,sn)
                istat=addBoundary(nGrid,p1,p2,xini,xfin,iye(p1,ngrid)+1,iye(p1,ngrid)+gz, &
		      1,getNz(ngrid),tagCount,pRecv,sn)
	        p2=getMapProcessor(ngrid,ixe(p2,ngrid)+1,iye(p1,ngrid)+1)
                IF(getPNeighbour(p2,4,ngrid)/=0) THEN
	           xini=ixb(p2,ngrid)-gz
		ELSE
	           xini=ixb(p2,ngrid)+1
		END IF
                IF(getPNeighbour(p1,3,ngrid)/=0) THEN
	           xfin=ixe(p1,ngrid)+gz
		ELSE
	           xfin=ixe(p1,ngrid)+1
		END IF
	        tagCount=tagCount+1
                istat=addBoundary(nGrid,p1,p2,xini,xfin,iye(p1,ngrid)+1-gz, &
		       iye(p1,ngrid),1,getNz(ngrid),tagCount,pSend,ns)
                istat=addBoundary(nGrid,p2,p1,xini,xfin,iye(p1,ngrid)+1-gz, &
		       iye(p1,ngrid),1,getNz(ngrid),tagCount,pRecv,ns)
		tagCount=tagCount+1       
                istat=addBoundary(nGrid,p2,p1,xini,xfin,iye(p1,ngrid)+1, &
		       iye(p1,ngrid)+gz,1,getNz(ngrid),tagCount,pSend,sn)
                istat=addBoundary(nGrid,p1,p2,xini,xfin,iye(p1,ngrid)+1, &
		       iye(p1,ngrid)+gz,1,getNz(ngrid),tagCount,pRecv,sn)
	      ELSE
                IF(getPNeighbour(p1,3,ngrid)/=0) THEN
	           xfin=ixe(p1,ngrid)+gz
		ELSE
		   xfin=ixe(p1,ngrid)+1
		END IF
	        tagCount=tagCount+1
                istat=addBoundary(nGrid,p1,p2,xini,xfin,iye(p1,ngrid)+1-gz, &
		      iye(p1,ngrid),1,getNz(ngrid),tagCount,pSend,ns)
                istat=addBoundary(nGrid,p2,p1,xini,xfin,iye(p1,ngrid)+1-gz, &
		      iye(p1,ngrid),1,getNz(ngrid),tagCount,pRecv,ns)
		tagCount=tagCount+1
                istat=addBoundary(nGrid,p2,p1,xini,xfin,iye(p1,ngrid)+1, &
		      iye(p1,ngrid)+gz,1,getNz(ngrid),tagCount,pSend,sn)
                istat=addBoundary(nGrid,p1,p2,xini,xfin,iye(p1,ngrid)+1, &
		      iye(p1,ngrid)+gz,1,getNz(ngrid),tagCount,pRecv,sn)
	      END IF
	  END IF
	END DO
        
	   DO p1=1,nmachs
	      !WRITE (*,'("DB: ",7(I3.3,1X))') &
	      !     p1,ixb(p1,ngrid),ixe(p1,ngrid),iyb(p1,ngrid),iye(p1,ngrid),incx(p1),incy(p1); CALL flush(6)
 	      CALL SetLimits(ngrid,p1,ixb(p1,ngrid),ixe(p1,ngrid),iyb(p1,ngrid),iye(p1,ngrid))
       END DO
	   !CALL dumpBounds() !debug

	   CALL adjust(incx,incy,ngrid,nmachs,m2,m3)
    END DO
   
  END SUBROUTINE InitBounds
   
  SUBROUTINE sendInit(nmachs,mynum, ngrids)
    USE ProcessorMod, ONLY: getPNeighbour,setNeighbour
    USE BoundaryMod, ONLY: bound_data,getLAstBoundary,getNeighbour,getfirstLine,getLastLine, &
                           getFirstColumn,getLastColumn,getFirstLevel,getLastLevel, &
			   getTag,getSize,getOperation,getAxys,getLimits
    USE advMessageMod, ONLY: SendMessageI,RecvMessageI,createMessage,allocMessages, &
                             newM2,newM3,newIa,newIz,newJa,newJz,SendMessageJ,RecvMessageJ, &
			     totalSendI,totalSendJ,TotalRecvI,totalRecvJ			   
    include "mpif.h"
    
    INTEGER, INTENT(IN) :: nmachs,mynum,ngrids
    INTEGER :: np,dir,istat,ierr,soc,i,ng,lb,nb,isi,iri,gz
    INTEGER :: socl,iia,iiz,ija,ijz,isj,irj,lbl
    INTEGER :: nSEndI,nsendJ,nRecvI,nRecvJ
    LOGICAL :: ing(4)
    INTEGER :: handle(nmachs)
    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER,ALLOCATABLE,DIMENSION(:) :: datacomm,datacomm_local
    INTEGER :: l_mynum
    
    INTEGER :: info(9)
    INTEGER :: newM2l(ngrids),newM3l(ngrids),newIal(ngrids)
    INTEGER :: newIzl(ngrids),newJal(ngrids),newJzl(ngrids)
    
    
    
    l_mynum=mynum-1
    CALL createMessage(ngrids)
     IF(l_mynum==0) THEN  !IF I am the master
	DO np=0,nmachs-1
	   DO ng=1,ngrids
	      lb=getLAstBoundary(ng,np+1)
	      soc=lb*11
              IF(np==0) socl=soc
	      info(1)=np+1
	      info(2)=ng
	      info(3)=lb
	      CALL getLimits(ng,np+1,info(4),info(5),info(6),info(7),info(8),info(9))
              IF(np/=0) & !LFR-5.0  
                  CALL MPI_send(info,9, MPI_INTEGER,np,20000+np, &
                           MPI_COMM_WORLD, ierr)
	      ALLOCATE(dataComm(soc))
              IF(np==0) ALLOCATE(dataComm_local(socl))               
	      i=1
	      DO nb=1,lb
               !print *, '-------------------------------'
                 !print *,'np,lb: ',np,nb;print *, '-------------------------------'
               	 dataComm(i)=getNeighbour(ng,np+1,nb)!-1 ;print *,'Neighbour: ',dataComm(i)
               	 dataComm(i+1)=getfirstLine(ng,np+1,nb) !;print *,'firstLine:',dataComm(i+1)
               	 dataComm(i+2)=getLAstLine(ng,np+1,nb) !;print *,'LAstLine:',dataComm(i+2)
               	 dataComm(i+3)=getfirstColumn(ng,np+1,nb) !;print *,'firstColumn:',dataComm(i+3)
               	 dataComm(i+4)=getLastColumn(ng,np+1,nb) !;print *,'LastColumn:',dataComm(i+4)
               	 dataComm(i+5)=getfirstLevel(ng,np+1,nb) !;print *,'firstLevel:',dataComm(i+5)
               	 dataComm(i+6)=getLastLevel(ng,np+1,nb) !;print *,'LastLevel:',dataComm(i+6)
               	 dataComm(i+7)=getSize(ng,np+1,nb) !;print *,'Size:',dataComm(i+7)
               	 dataComm(i+8)=getTag(ng,np+1,nb) !;print *,'Tag:',dataComm(i+8)
               	 dataComm(i+9)=getOperation(ng,np+1,nb) !;print *,'Operation:',dataComm(i+9)
               	 dataComm(i+10)=getAxys(ng,np+1,nb) !;print *,'Axys:',dataComm(i+10); CALL flush(6)
                 IF(np==0) THEN
                     dataComm_local(i)=getNeighbour(ng,np+1,nb)!-1! ;print *,'Neighbour: ',dataComm_local(i)
                     dataComm_local(i+1)=getfirstLine(ng,np+1,nb)! ;print *,'firstLine:',dataComm_local(i+1)
                     dataComm_local(i+2)=getLAstLine(ng,np+1,nb)! ;print *,'LAstLine:',dataComm_local(i+2)
                     dataComm_local(i+3)=getfirstColumn(ng,np+1,nb)! ;print *,'firstColumn:',dataComm_local(i+3)
                     dataComm_local(i+4)=getLastColumn(ng,np+1,nb)! ;print *,'LastColumn:',dataComm_local(i+4)
                     dataComm_local(i+5)=getfirstLevel(ng,np+1,nb)! ;print *,'firstLevel:',dataComm_local(i+5)
                     dataComm_local(i+6)=getLastLevel(ng,np+1,nb)! ;print *,'LastLevel:',dataComm_local(i+6)
                     dataComm_local(i+7)=getSize(ng,np+1,nb)! ;print *,'Size:',dataComm_local(i+7)
                     dataComm_local(i+8)=getTag(ng,np+1,nb)! ;print *,'Tag:',dataComm_local(i+8)
                     dataComm_local(i+9)=getOperation(ng,np+1,nb)! ;print *,'Operation:',dataComm_local(i+9)
                     dataComm_local(i+10)=getAxys(ng,np+1,nb)! ;print *,'Axys:',dataComm_local(i+10); CALL flush(6)
                 END IF
		 i=i+11
	     END DO
             IF(np/=0) & !LFR-5.0
               CALL MPI_send(dataComm,soc, MPI_INTEGER,np,30000+np, &
                           MPI_COMM_WORLD, ierr)
              IF(np==0) THEN !copia para variaveis locais ao inves de enviar
                  lbl=info(3)
                  newM2l(ng)=info(4)
                  newM3l(ng)=info(5)
                  newIal(ng)=info(6)
                  newIzl(ng)=info(7)
                  newJal(ng)=info(8)
                  newJzl(ng)=info(9)
              end if
	     DEALLOCATE(dataComm)
	  END DO
       END DO    
   END IF
!LFR-5.0     ELSE               !Else i am a slave
       DO ng=1,ngrids
         IF(l_mynum/=0) THEN  !LFR-5.0  
            CALL MPI_recv(info,9, MPI_INTEGER,0,20000+l_mynum, & !LFR - 5.0
               MPI_COMM_WORLD, status, ierr)
	      !IF(mynum==275) Write(*,fmt='("Recv info:",6(I4.4,1X))') &
	      !            info(4),info(5),info(6),info(7),info(8),info(9)
	  newM2(ng)=info(4)
	  newM3(ng)=info(5)
	  newIa(ng)=info(6)
	  newIz(ng)=info(7)
	  newJa(ng)=info(8)
	  newJz(ng)=info(9)
	  soc=info(3)*11 !Number of boundaries
          lb=info(3)
         ELSE
            newM2(ng)=newM2l(ng)
            newM3(ng)=newM3l(ng)
            newIa(ng)=newIal(ng)
            newIz(ng)=newIzl(ng)
            newJa(ng)=newJal(ng)
            newJz(ng)=newJzl(ng)
            lb=lbl
         END IF
         !print *,'mynum,lb: ',mynum,lb
!srf
! if(newIz(ng)+1 .gt. newM2(ng)) then
!   newM3(ng)=max(newM3(ng),newjz(ng)+3)
   !print*,'MNTERROR X=',newIz(ng),newM2(ng),mynum
  ! call flush(6)
   !stop 1001
! endif
! if(newJz(ng)+1 .gt. newM3(ng)) then
  ! print*,'MNTERROR Y=',newjz(ng),newM3(ng),mynum
   !call flush(6)
   !stop 1002
! endif

         ALLOCATE(datacomm(soc))

         IF(l_mynum/=0) THEN !LFR-5.0
            CALL MPI_recv(datacomm,soc, MPI_INTEGER,0,30000+l_mynum, & !LFR - 5.0
               MPI_COMM_WORLD, status, ierr)
         ELSE !Se eh o processador 1 apenas copia os dados pois eh local
            DO i=1,socl
               datacomm(i)=datacomm_local(i)
            END DO
         END IF

          !Counting the messages by types
	  nRecvI=0;nRecvJ=0;nSendI=0;nSendJ=0
	  i=1
	  DO nb=1,lb
             IF(dataComm(i+9)==1) THEN !Send
	        IF(dataComm(i+10)==1 .or. dataComm(i+10)==2) THEN !x axys
		   nsendI=nsendI+1
		ELSE ! y axys
		   nsendJ=nsendJ+1
		END IF
             ELSE
	        IF(dataComm(i+10)==1 .or. dataComm(i+10)==2) THEN !x axys
		   nRecvI=nRecvI+1
		ELSE ! y axys
		   nRecvJ=nRecvJ+1
		END IF
	     END IF
	     i=i+11
	  END DO

	  CALL allocMessages(ng,nSendI,.true.,.true.) !IstoSend, IsI
	  CALL allocMessages(ng,nRecvI,.false.,.true.) !IstoRecv, IsI
	  CALL allocMessages(ng,nSEndJ,.true.,.false.) !IstoSend, IsJ
	  CALL allocMessages(ng,nRecvJ,.false.,.false.) !IstoREcv, IsJ
	  i=1
	  isi=0;isj=0
	  iri=0;irj=0
	  DO nb=1,lb
             IF(dataComm(i+9)==1) THEN !SEnd
	        IF(dataComm(i+10)==1 .or. dataComm(i+10)==2) THEN !x axys
	           isi=isi+1
                   SendMessageI(ng)%Proc(isi)= dataComm(i+0)
                   SendMessageI(ng)%ia(isi)=	  dataComm(i+1)
                   SendMessageI(ng)%iz(isi)=	  dataComm(i+2)
                   SendMessageI(ng)%ja(isi)=	  dataComm(i+3)
                   SendMessageI(ng)%jz(isi)=	  dataComm(i+4)
                   SendMessageI(ng)%za(isi)=	  dataComm(i+5)
                   SendMessageI(ng)%zz(isi)=	  dataComm(i+6)
                   SendMessageI(ng)%mSize(isi)=dataComm(i+7)
                   SendMessageI(ng)%tag(isi)=  dataComm(i+8)	
	           SendMessageI(ng)%mSize(isi)=(SendMessageI(ng)%iz(isi)- &
	                                SendMessageI(ng)%ia(isi)+1)* &
					(SendMessageI(ng)%jz(isi) - &
					 SendMessageI(ng)%ja(isi)+1) * &
					 ( SendMessageI(ng)%zz(isi) - &
					 SendMessageI(ng)%za(isi)+1)
                   SendMessageI(ng)%start(isi)=totalsendI(ng)+1
		   totalsendI(ng)=totalsendI(ng)+SendMessageI(ng)%mSize(isi)
                   !print *,'Me, totalsendI: ',mynum,totalsendI(ng); call flush(6)
		ELSE

		   isj=isj+1
                   SendMessageJ(ng)%Proc(isj)= dataComm(i+0)
                   SendMessageJ(ng)%ia(isj)=	  dataComm(i+1)
                   SendMessageJ(ng)%iz(isj)=	  dataComm(i+2)
                   SendMessageJ(ng)%ja(isj)=	  dataComm(i+3)
                   SendMessageJ(ng)%jz(isj)=	  dataComm(i+4)
                   SendMessageJ(ng)%za(isj)=	  dataComm(i+5)
                   SendMessageJ(ng)%zz(isj)=	  dataComm(i+6)
                   SendMessageJ(ng)%mSize(isj)=dataComm(i+7)
                   SendMessageJ(ng)%tag(isj)=  dataComm(i+8)	
	           SendMessageJ(ng)%mSize(isj)=(SendMessageJ(ng)%iz(isj)- &
	                                SendMessageJ(ng)%ia(isj)+1)* &
					(SendMessageJ(ng)%jz(isj) - &
					 SendMessageJ(ng)%ja(isj)+1) * &
					 ( SendMessageJ(ng)%zz(isj) - &
					 SendMessageJ(ng)%za(isj)+1)
                   SendMessageJ(ng)%start(isj)=totalsendJ(ng)+1
		   totalsendJ(ng)=totalsendJ(ng)+SendMessageJ(ng)%mSize(isj)		
                   !print *,'Me, totalsendJ: ',mynum,totalsendJ(ng); call flush(6)
	        END IF
	     ELSE IF(dataComm(i+9)==2) THEN !Recv
	        IF(dataComm(i+10)==1 .or. dataComm(i+10)==2) THEN !x axys
	           iri=iri+1
 	           RecvMessageI(ng)%Proc(iri)= dataComm(i+0)
	           RecvMessageI(ng)%ia(iri)=	dataComm(i+1)
	           RecvMessageI(ng)%iz(iri)=	dataComm(i+2)
	           RecvMessageI(ng)%ja(iri)=	dataComm(i+3)
	           RecvMessageI(ng)%jz(iri)=	dataComm(i+4)
	           RecvMessageI(ng)%za(iri)=	     dataComm(i+5)
	           RecvMessageI(ng)%zz(iri)=	     dataComm(i+6)
	           RecvMessageI(ng)%mSize(iri)=dataComm(i+7)
	           RecvMessageI(ng)%tag(iri)=  dataComm(i+8)			   
	           RecvMessageI(ng)%mSize(iri)=(RecvMessageI(ng)%iz(iri)- &
				   RecvMessageI(ng)%ia(iri)+1)* &
				   (RecvMessageI(ng)%jz(iri) - &
				    RecvMessageI(ng)%ja(iri)+1) * &
				    ( RecvMessageI(ng)%zz(iri) - &
				    RecvMessageI(ng)%za(iri)+1)
                   RecvMessageI(ng)%start(iri)=totalRecvI(ng)+1
		   totalRecvI(ng)=totalRecvI(ng)+RecvMessageI(ng)%mSize(iri)
                   !print *,'Me, totalRecvI: ',mynum,totalRecvI(ng); call flush(6)
	        ELSE
		   irj=irj+1
 	           RecvMessageJ(ng)%Proc(irj)= dataComm(i+0)
	           RecvMessageJ(ng)%ia(irj)=	dataComm(i+1)
	           RecvMessageJ(ng)%iz(irj)=	dataComm(i+2)
	           RecvMessageJ(ng)%ja(irj)=	dataComm(i+3)
	           RecvMessageJ(ng)%jz(irj)=	dataComm(i+4)
	           RecvMessageJ(ng)%za(irj)=	     dataComm(i+5)
	           RecvMessageJ(ng)%zz(irj)=	     dataComm(i+6)
	           RecvMessageJ(ng)%mSize(irj)=dataComm(i+7)
	           RecvMessageJ(ng)%tag(irj)=  dataComm(i+8)			   
	           RecvMessageJ(ng)%mSize(irj)=(RecvMessageJ(ng)%iz(irj)- &
				   RecvMessageJ(ng)%ia(irj)+1)* &
				   (RecvMessageJ(ng)%jz(irj) - &
				    RecvMessageJ(ng)%ja(irj)+1) * &
				    ( RecvMessageJ(ng)%zz(irj) - &
				    RecvMessageJ(ng)%za(irj)+1)
		
                   RecvMessageJ(ng)%start(irj)=totalRecvJ(ng)+1
		   totalRecvJ(ng)=totalRecvJ(ng)+RecvMessageJ(ng)%mSize(irj)
                   !print *,'Me, totalRecvJ: ',mynum,totalrecvJ(ng); call flush(6)
		END IF   
	     END IF
	     i=i+11
	  END DO
	  DEALLOCATE(datacomm)
       END DO
!LFR-5.0 END IF

  END SUBROUTINE sendInit
  
        
END MODULE AdvectData
