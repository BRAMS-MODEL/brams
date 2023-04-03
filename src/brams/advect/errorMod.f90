MODULE ErrorMod
   PRIVATE
   
   TYPE type_error
      CHARACTER(LEN=80) :: message
      INTEGER :: class  !1=warning, 2=fatal
      LOGICAL :: action !T=stop, F=Continue
   END TYPE type_error
   TYPE(type_error), ALLOCATABLE,DIMENSION(:) :: error
   
   PUBLIC CreateError,DestroyError,setError,printError
   
   CONTAINS
   
  INTEGER FUNCTION CreateError(size)
     INTEGER, INTENT(IN) :: size
 !    
 !    IF(allocated(error)) THEN
         CreateError=0
 !    ELSE
 !       allocate(error(size))
 !    END IF
 !    
  END FUNCTION CreateError
 ! 
  INTEGER FUNCTION DestroyError()
 !    
 !    IF(.not. allocated(error)) THEN
         DestroyError=0
 !    ELSE
 !       deallocate(error)
 !    END IF
 !    
  END FUNCTION DestroyError     
!
  INTEGER FUNCTION setError(vnumber,vmessage,vclass, vaction)
     INTEGER, INTENT(IN) :: vnumber
     INTEGER, INTENT(IN) :: vclass
     LOGICAL, INTENT(IN) :: vaction
     CHARACTER(LEN=*), INTENT(IN) :: vmessage
 !    
 !    IF(.not. allocated(error) .or. vnumber>sizeOf(errors)) THEN
        setError=0
 !    ELSE
 !       error(vnumber)%message=vmessage
!	 error(vnumber)%class=vclass
!	 error(vnumber)%action=vaction
 !       setError=-1
 !    END IF
  END FUNCTION setError
 ! 
  FUNCTION printError(vnumber,newMessage, filenum)
     INTEGER, INTENT(IN) :: vnumber
     INTEGER, INTENT(IN) :: filenum
     CHARACTER(LEN=*), INTENT(IN) :: newMessage
 !    
 !    INTEGER :: fN
 !    
 !    IF(filenum==0) THEN 
 !       fN=6
 !    ELSE
 !       fN=filenum
 !    END IF
 !    
 !    IF(vnumber>sizeOf(errors) .or. vnumber<=0) THEN
 !       WRITE (fN,*) '============= Unespecified Error ============'; CALL flush(fN)
!	 WRITE (fN,*) trim(newMessage); CALL flush(fN)
!	 WRITE (fN,*) '========================================='; CALL flush(fN)
 !    ELSE
 !       IF(error(vnumber)%class==1) THEN
!	    WRITE (fN,*) '=============== Warning =================' ; CALL flush(fN)
!	    WRITE (fN,*) error(vnumber)%message; CALL flush(fN)
!	    WRITE (fN,*) trim(newMessage); CALL flush(fN)
!	    WRITE (fN,*) '========================================='; CALL flush(fN)
!	 ELSE
!	    WRITE (fN,*) '============= Fatal Error ===============' ; CALL flush(fN)
!	    WRITE (fN,*) error(vnumber)%message; CALL flush(fN)
!	    WRITE (fN,*) trim(newMessage); CALL flush(fN)
!	    WRITE (fN,*) '========================================='; CALL flush(fN)
!	    !STOP
!	 END IF
 !   END IF
!
      printError=0
  END FUNCTION printError        
   
END MODULE ErrorMod
