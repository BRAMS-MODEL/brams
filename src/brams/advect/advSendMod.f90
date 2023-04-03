MODULE advMessageMod
   IMPLICIT NONE
   !PRIVATE
   
   TYPE t_message
     INTEGER,POINTER,DIMENSION(:)  :: Proc
     INTEGER,POINTER,DIMENSION(:)  :: ia
     INTEGER,POINTER,DIMENSION(:)  :: iz
     INTEGER,POINTER,DIMENSION(:)  :: ja
     INTEGER,POINTER,DIMENSION(:)  :: jz
     INTEGER,POINTER,DIMENSION(:)  :: za
     INTEGER,POINTER,DIMENSION(:)  :: zz
     INTEGER,POINTER,DIMENSION(:)  :: mSize
     INTEGER,POINTER,DIMENSION(:)  :: tag
     INTEGER,POINTER,DIMENSION(:)  :: start
   END TYPE t_message
   TYPE(t_message),ALLOCATABLE,DIMENSION(:) :: SendMessageI !ngrids,nrecv
   TYPE(t_message),ALLOCATABLE,DIMENSION(:) :: RecvMessageI !ngrids,nsend
   TYPE(t_message),ALLOCATABLE,DIMENSION(:) :: SendMessageJ !ngrids,nrecv
   TYPE(t_message),ALLOCATABLE,DIMENSION(:) :: RecvMessageJ !ngrids,nsend
   
   INTEGER,ALLOCATABLE,DIMENSION(:) :: nsendI,nrecvI,nsendJ,nrecvJ,totalSendI,totalSendJ,TotalRecvI,totalRecvJ
   INTEGER,ALLOCATABLE,DIMENSION(:) :: newIa,newJa,newIz,newJz,newM2,newM3 !ngrids
   LOGICAL :: isAllocated
   INTEGER :: TotalMessages
  
   PUBLIC :: nsendI,nrecvI,SendMessageI,RecvMessageI,totalSendI,totalSendJ,TotalRecvI,totalRecvJ
   PUBLIC :: nsendJ,nrecvJ,SendMessageJ,RecvMessageJ
   PUBLIC :: createMessage,allocMessages
   PUBLIC :: dumpMessages,newM2,newM3,newIa,newJa,newIz,newJz

   CONTAINS
   
   SUBROUTINE createMessage(ngrids)
   
     INTEGER,INTENT(IN) :: ngrids
     INTEGER :: i   
  
     IF(isAllocated) THEN
        PRINT *,'Error: Messages already allocated!'
	STOP
     END IF
!     PRINT *,'Allocating messages';call flush(6)
     ALLOCATE(nsendI(ngrids))
     ALLOCATE(nrecvI(ngrids))
     ALLOCATE(nsendJ(ngrids))
     ALLOCATE(nrecvJ(ngrids))
     ALLOCATE(SendMessageI(ngrids))
     ALLOCATE(RecvMessageI(ngrids))
     ALLOCATE(SendMessageJ(ngrids))
     ALLOCATE(RecvMessageJ(ngrids))
     ALLOCATE(newM2(ngrids),newM3(ngrids),newIa(ngrids),newIz(ngrids),newJa(ngrids),newJz(ngrids))
     ALLOCATE(totalSendI(ngrids),totalSendJ(ngrids),TotalRecvI(ngrids),totalRecvJ(ngrids))
     totalSendI=0
     totalSendJ=0
     TotalRecvI=0
     totalRecvJ=0
!     PRINT *,'Allocated messages';call flush(6)
     DO i=1,ngrids
        IF(ASSOCIATED(SendMessageI(i)%Proc))  NULLIFY(SendMessageI(i)%Proc)
        IF(ASSOCIATED(SendMessageI(i)%ia))    NULLIFY(SendMessageI(i)%ia)    
        IF(ASSOCIATED(SendMessageI(i)%iz))	 NULLIFY(SendMessageI(i)%iz)
        IF(ASSOCIATED(SendMessageI(i)%ja))	 NULLIFY(SendMessageI(i)%ja)
        IF(ASSOCIATED(SendMessageI(i)%jz))	 NULLIFY(SendMessageI(i)%jz)
        IF(ASSOCIATED(SendMessageI(i)%za))	 NULLIFY(SendMessageI(i)%za)
        IF(ASSOCIATED(SendMessageI(i)%zz))	 NULLIFY(SendMessageI(i)%zz)
        IF(ASSOCIATED(SendMessageI(i)%mSize)) NULLIFY(SendMessageI(i)%mSize)
        IF(ASSOCIATED(SendMessageI(i)%tag))	 NULLIFY(SendMessageI(i)%tag)
        IF(ASSOCIATED(SendMessageI(i)%start))	 NULLIFY(SendMessageI(i)%start)
        IF(ASSOCIATED(REcvMessageI(i)%Proc))  NULLIFY(REcvMessageI(i)%Proc)
        IF(ASSOCIATED(REcvMessageI(i)%ia))    NULLIFY(REcvMessageI(i)%ia)    
        IF(ASSOCIATED(REcvMessageI(i)%iz))	 NULLIFY(REcvMessageI(i)%iz)
        IF(ASSOCIATED(REcvMessageI(i)%ja))	 NULLIFY(REcvMessageI(i)%ja)
        IF(ASSOCIATED(REcvMessageI(i)%jz))	 NULLIFY(REcvMessageI(i)%jz)
        IF(ASSOCIATED(REcvMessageI(i)%za))	 NULLIFY(REcvMessageI(i)%za)
        IF(ASSOCIATED(REcvMessageI(i)%zz))	 NULLIFY(REcvMessageI(i)%zz)
        IF(ASSOCIATED(REcvMessageI(i)%mSize)) NULLIFY(REcvMessageI(i)%mSize)
        IF(ASSOCIATED(REcvMessageI(i)%tag))	 NULLIFY(REcvMessageI(i)%tag)
        IF(ASSOCIATED(REcvMessageI(i)%start))	 NULLIFY(REcvMessageI(i)%start)
        IF(ASSOCIATED(SendMessageJ(i)%Proc))  NULLIFY(SendMessageJ(i)%Proc)
        IF(ASSOCIATED(SendMessageJ(i)%ia))    NULLIFY(SendMessageJ(i)%ia)    
        IF(ASSOCIATED(SendMessageJ(i)%iz))	 NULLIFY(SendMessageJ(i)%iz)
        IF(ASSOCIATED(SendMessageJ(i)%ja))	 NULLIFY(SendMessageJ(i)%ja)
        IF(ASSOCIATED(SendMessageJ(i)%jz))	 NULLIFY(SendMessageJ(i)%jz)
        IF(ASSOCIATED(SendMessageJ(i)%za))	 NULLIFY(SendMessageJ(i)%za)
        IF(ASSOCIATED(SendMessageJ(i)%zz))	 NULLIFY(SendMessageJ(i)%zz)
        IF(ASSOCIATED(SendMessageJ(i)%mSize)) NULLIFY(SendMessageJ(i)%mSize)
        IF(ASSOCIATED(SendMessageJ(i)%tag))	 NULLIFY(SendMessageJ(i)%tag)
        IF(ASSOCIATED(SendMessageJ(i)%start))	 NULLIFY(SendMessageJ(i)%start)
        IF(ASSOCIATED(REcvMessageJ(i)%Proc))  NULLIFY(REcvMessageJ(i)%Proc)
        IF(ASSOCIATED(REcvMessageJ(i)%ia))    NULLIFY(REcvMessageJ(i)%ia)    
        IF(ASSOCIATED(REcvMessageJ(i)%iz))	 NULLIFY(REcvMessageJ(i)%iz)
        IF(ASSOCIATED(REcvMessageJ(i)%ja))	 NULLIFY(REcvMessageJ(i)%ja)
        IF(ASSOCIATED(REcvMessageJ(i)%jz))	 NULLIFY(REcvMessageJ(i)%jz)
        IF(ASSOCIATED(REcvMessageJ(i)%za))	 NULLIFY(REcvMessageJ(i)%za)
        IF(ASSOCIATED(REcvMessageJ(i)%zz))	 NULLIFY(REcvMessageJ(i)%zz)
        IF(ASSOCIATED(REcvMessageJ(i)%mSize)) NULLIFY(REcvMessageJ(i)%mSize)
        IF(ASSOCIATED(REcvMessageJ(i)%tag))	 NULLIFY(REcvMessageJ(i)%tag)
        IF(ASSOCIATED(REcvMessageJ(i)%start))	 NULLIFY(REcvMessageJ(i)%start)
     END DO
!     PRINT *,'Nullifyed messages';call flush(6)
     
   END SUBROUTINE createMessage

   SUBROUTINE allocMessages(ngrid,nmess,is2Send,isI)
      INTEGER, INTENT(IN) :: ngrid,nmess
      LOGICAL, INTENT(IN) :: is2Send,isI
      
      IF(is2Send) THEN
         IF(isI) THEN
            nsendI(ngrid)=nmess
            ALLOCATE(SendMessageI(ngrid)%Proc(nmess)) 
            ALLOCATE(SendMessageI(ngrid)%ia(nmess))   
            ALLOCATE(SendMessageI(ngrid)%iz(nmess))   
            ALLOCATE(SendMessageI(ngrid)%ja(nmess))   
            ALLOCATE(SendMessageI(ngrid)%jz(nmess))   
            ALLOCATE(SendMessageI(ngrid)%za(nmess))   
            ALLOCATE(SendMessageI(ngrid)%zz(nmess))   
            ALLOCATE(SendMessageI(ngrid)%mSize(nmess))
            ALLOCATE(SendMessageI(ngrid)%tag(nmess))   
            ALLOCATE(SendMessageI(ngrid)%start(nmess))   
        ELSE
            nsendJ(ngrid)=nmess
            ALLOCATE(SendMessageJ(ngrid)%Proc(nmess)) 
            ALLOCATE(SendMessageJ(ngrid)%ia(nmess))   
            ALLOCATE(SendMessageJ(ngrid)%iz(nmess))   
            ALLOCATE(SendMessageJ(ngrid)%ja(nmess))   
            ALLOCATE(SendMessageJ(ngrid)%jz(nmess))   
            ALLOCATE(SendMessageJ(ngrid)%za(nmess))   
            ALLOCATE(SendMessageJ(ngrid)%zz(nmess))   
            ALLOCATE(SendMessageJ(ngrid)%mSize(nmess))
            ALLOCATE(SendMessageJ(ngrid)%tag(nmess))   
            ALLOCATE(SendMessageJ(ngrid)%start(nmess))   
	END IF
      ELSE
         IF(isI) THEN
            nrecvI(ngrid)=nmess
            ALLOCATE(RecvMessageI(ngrid)%Proc(nmess)) 
            ALLOCATE(RecvMessageI(ngrid)%ia(nmess))   
            ALLOCATE(RecvMessageI(ngrid)%iz(nmess))   
            ALLOCATE(RecvMessageI(ngrid)%ja(nmess))   
            ALLOCATE(RecvMessageI(ngrid)%jz(nmess))   
            ALLOCATE(RecvMessageI(ngrid)%za(nmess))   
            ALLOCATE(RecvMessageI(ngrid)%zz(nmess))   
            ALLOCATE(RecvMessageI(ngrid)%mSize(nmess))
            ALLOCATE(RecvMessageI(ngrid)%tag(nmess))   
            ALLOCATE(RecvMessageI(ngrid)%start(nmess))   
         ELSE
            nrecvJ(ngrid)=nmess
	    ALLOCATE(RecvMessageJ(ngrid)%Proc(nmess)) 
            ALLOCATE(RecvMessageJ(ngrid)%ia(nmess))   
            ALLOCATE(RecvMessageJ(ngrid)%iz(nmess))   
            ALLOCATE(RecvMessageJ(ngrid)%ja(nmess))   
            ALLOCATE(RecvMessageJ(ngrid)%jz(nmess))   
            ALLOCATE(RecvMessageJ(ngrid)%za(nmess))   
            ALLOCATE(RecvMessageJ(ngrid)%zz(nmess))   
            ALLOCATE(RecvMessageJ(ngrid)%mSize(nmess))
            ALLOCATE(RecvMessageJ(ngrid)%tag(nmess))   
            ALLOCATE(RecvMessageJ(ngrid)%start(nmess))   
	 END IF
      END IF  

   END SUBROUTINE allocMessages

  SUBROUTINE dumpMessages(fout,ngrids)
       INTEGER, INTENT(IN) :: fout,ngrids
       
       INTEGER :: ng,m
       DO ng=1,ngrids
	  WRITE (fout,FMT='(A,I3.3,A)') '============================ Grid ',ng,'=============================='

	  WRITE (fout,FMT='("#Sends I: ",I5.5," #Sends J: ",I5.5," #Recs I:		   ",I5.5," #Recs J: ",I5.5 )') &
	                nsendI(ng),nsendJ(ng),nrecvI(ng),nrecvJ(ng)	  

	  WRITE (fout,FMT='(A)') '=================================================================================='
	  IF(nSendI(ng)>0) THEN
	  WRITE (fout,FMT='(A,I8.8,A)') '==== Send data ==== total ',totalSendI(ng),' ==== I direction ===='
	  DO m=1,nSendI(ng)
	     WRITE(fout,FMT='("#",I3.3," - Send     to  ",I5.5,1X,6(I3.3,1X),I8.8,1X,I5.5,1X,I5.5)') & 
			      m,SendMessageI(ng)%Proc(m), &
			      SendMessageI(ng)%ia(m) , &    
			      SendMessageI(ng)%iz(m) , &  
			      SendMessageI(ng)%ja(m) , &  
			      SendMessageI(ng)%jz(m) , &  
			      SendMessageI(ng)%za(m) , &  
			      SendMessageI(ng)%zz(m) , &  
			      SendMessageI(ng)%mSize(m), &
			      SendMessageI(ng)%tag(m), &
			      SendMessageI(ng)%start(m)
         END DO
         END IF	 
         IF(nsendJ(ng)>0) THEN
	 WRITE (fout,FMT='(A,I8.8,A)') '==== Send data ==== total ',totalSendJ(ng),' ==== J direction ===='	  
	 DO m=1,nSendJ(ng)
	     WRITE(fout,FMT='("#",I3.3," - Send     to  ",I5.5,1X,6(I3.3,1X),I8.8,1X,I5.5,1X,I5.5)') & 
			      m,SendMessageJ(ng)%Proc(m), &
			      SendMessageJ(ng)%ia(m) , &    
			      SendMessageJ(ng)%iz(m) , &  
			      SendMessageJ(ng)%ja(m) , &  
			      SendMessageJ(ng)%jz(m) , &  
			      SendMessageJ(ng)%za(m) , &  
			      SendMessageJ(ng)%zz(m) , &  
			      SendMessageJ(ng)%mSize(m), &
			      SendMessageJ(ng)%tag(m), &
			      SendMessageJ(ng)%start(m)
         END DO
         END IF
	 IF(nRecvI(ng)>0) THEN
	 WRITE (fout,FMT='(A,I8.8,A)') '==== Recv data ==== total ',totalRecvI(ng),' ==== I direction ===='	  
	  DO m=1,nRecvI(ng)
	     WRITE(fout,FMT='("#",I3.3," - Receive from ",I5.5,1X,6(I3.3,1X),I8.8,1X,I5.5,1X,I5.5)') & 
			      m,RecvMessageI(ng)%Proc(m), &
			      RecvMessageI(ng)%ia(m) , &    
			      RecvMessageI(ng)%iz(m) , &  
			      RecvMessageI(ng)%ja(m) , &  
			      RecvMessageI(ng)%jz(m) , &  
			      RecvMessageI(ng)%za(m) , &  
			      RecvMessageI(ng)%zz(m) , &  
			      RecvMessageI(ng)%mSize(m), &
			      RecvMessageI(ng)%tag(m), &
			      RecvMessageI(ng)%start(m)
         END DO
         END IF
	 IF(nRecvJ(ng)>0) THEN
	 WRITE (fout,FMT='(A,I8.8,A)') '==== Recv data ==== total ',totalRecvJ(ng),' ==== J direction ===='	  
	  DO m=1,nRecvJ(ng)
	     WRITE(fout,FMT='("#",I3.3," - Receive from ",I5.5,1X,6(I3.3,1X),I8.8,1X,I5.5,1X,I5.5)') & 
			      m,RecvMessageJ(ng)%Proc(m), &
			      RecvMessageJ(ng)%ia(m) , &    
			      RecvMessageJ(ng)%iz(m) , &  
			      RecvMessageJ(ng)%ja(m) , &  
			      RecvMessageJ(ng)%jz(m) , &  
			      RecvMessageJ(ng)%za(m) , &  
			      RecvMessageJ(ng)%zz(m) , &  
			      RecvMessageJ(ng)%mSize(m), &
			      RecvMessageJ(ng)%tag(m), &
			      RecvMessageJ(ng)%start(m)
         END DO
         END IF
      END DO
      CALL flush(fout)
   END SUBROUTINE dumpMessages
  
END MODULE advMessageMod
