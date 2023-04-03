MODULE mem_spack

  ! USE Spack_utils,only: NOB_REAL=>NOB,maxblock_size

  USE Spack_utils,ONLY: &
       maxblock_size      ! IN

  USE chem1_list, ONLY: &
       nspecies,        & ! PARAMETER
       nr,              & ! PARAMETER
       nr_photo,        & ! PARAMETER
       PhotojMethod       ! PARAMETER

  IMPLICIT NONE

  TYPE, PUBLIC :: spack_type

     !3d Real
     DOUBLE PRECISION,POINTER,DIMENSION(:,:,:) :: DLdrdc

     !2d Real
    DOUBLE PRECISION,POINTER,DIMENSION(:,:)   :: sc_p_new
    DOUBLE PRECISION,POINTER,DIMENSION(:,:)   :: sc_p_4
    DOUBLE PRECISION,POINTER,DIMENSION(:,:)   :: DLr	 
    DOUBLE PRECISION,POINTER,DIMENSION(:,:)   :: DLr3	 


     DOUBLE PRECISION,POINTER,DIMENSION(:,:)   :: jphoto 
     DOUBLE PRECISION,POINTER,DIMENSION(:,:)   :: rk    
     DOUBLE PRECISION,POINTER,DIMENSION(:,:)   :: w    
     DOUBLE PRECISION,POINTER,DIMENSION(:,:)   :: sc_p

     !1D Real
     DOUBLE PRECISION,POINTER,DIMENSION(:)     :: temp	
     DOUBLE PRECISION,POINTER,DIMENSION(:)     :: press 
     DOUBLE PRECISION,POINTER,DIMENSION(:)     :: cosz	  
     DOUBLE PRECISION,POINTER,DIMENSION(:)     :: att	  
     DOUBLE PRECISION,POINTER,DIMENSION(:)     :: vapp 
     DOUBLE PRECISION,POINTER,DIMENSION(:)     :: volmol 
     DOUBLE PRECISION,POINTER,DIMENSION(:)     :: volmol_i 
     DOUBLE PRECISION,POINTER,DIMENSION(:)     :: xlw 
     DOUBLE PRECISION,POINTER,DIMENSION(:)     :: err    

  END TYPE spack_type

  TYPE, PUBLIC :: spack_type_2d
     DOUBLE PRECISION,POINTER,DIMENSION(:,:)   :: DLmat 
     DOUBLE PRECISION,POINTER,DIMENSION(:)     :: DLb1	 
     DOUBLE PRECISION,POINTER,DIMENSION(:)     :: DLb2	 
     DOUBLE PRECISION,POINTER,DIMENSION(:)     :: DLb3	 
     DOUBLE PRECISION,POINTER,DIMENSION(:)     :: DLb4	 
     DOUBLE PRECISION,POINTER,DIMENSION(:)     :: DLk1    
     DOUBLE PRECISION,POINTER,DIMENSION(:)     :: DLk2    
     DOUBLE PRECISION,POINTER,DIMENSION(:)     :: DLk3    
     DOUBLE PRECISION,POINTER,DIMENSION(:)     :: DLk4    

  END TYPE spack_type_2d

  PRIVATE

  DOUBLE PRECISION, PARAMETER ::  rtols=1.D-3 ! 1e-2 means two digits
  DOUBLE PRECISION, PARAMETER ::  atols=1.D+7 ! Jacobson (1998, SMVGEAR) range 1.e3-1.e7! 1.D0

  DOUBLE PRECISION   , PUBLIC, DIMENSION(nspecies)        :: ATOL, RTOL 
  TYPE(spack_type)   , PUBLIC, ALLOCATABLE,DIMENSION(:)   :: spack
  TYPE(spack_type_2d), PUBLIC, ALLOCATABLE,DIMENSION(:,:) :: spack_2d
  LOGICAL            , PUBLIC                             :: spack_alloc = .FALSE.

  PUBLIC :: alloc_spack


CONTAINS
  !========================================================================

  SUBROUTINE alloc_spack(CHEMISTRY)
    
    implicit none
    !integer, intent(in) :: NOB_MEM
    INTEGER i,ii,NOB,n
    integer, intent(in) :: CHEMISTRY

    IF(spack_alloc) THEN
       PRINT *,'ERROR: spack_alloc already allocated'
       PRINT *,'Routine: spack_alloc File: mem_spack.f90'
       STOP
    END IF

    !  IF(NOB_MEM==0) NOB=NOB_REAL
    !  IF(NOB_MEM==1) NOB=1


    DO n=1,nspecies
       ATOL(n) = atols
       RTOL(n) = rtols
    END DO

    !- allocating Spaces to copy structure

    nob=1 ! to save memory, the scratch will be re-used
    ALLOCATE(spack(nob))

    !- maxblock_size is the maximum block size including all grids
    !- all grids/blocks will share the same array/memory area
    ALLOCATE(spack_2d(maxblock_size,nob))

    DO i=1,nob

    !- 3D Variables
    ALLOCATE(spack(i)%DLdrdc  (1:maxblock_size,nspecies,nspecies));spack(i)%DLdrdc  = 0.0D0	 

    !- 2D Variables
    ALLOCATE(spack(i)%jphoto  (1:maxblock_size,nr_photo))	  ;spack(i)%jphoto  = 0.0D0
    
    ALLOCATE(spack(i)%rk      (1:maxblock_size,nr))		  ;spack(i)%rk      = 0.0d0
    ALLOCATE(spack(i)%w       (1:maxblock_size,nr))		  ;spack(i)%w	    = 0.0d0
    ALLOCATE(spack(i)%sc_p    (1:maxblock_size,nspecies))	  ;spack(i)%sc_p    = 0.0D0
    ALLOCATE(spack(i)%sc_p_new(1:maxblock_size,nspecies))	  ;spack(i)%sc_p_new= 0.0D0
  
    ALLOCATE(spack(i)%DLr     (1:maxblock_size,nspecies))	  ;spack(i)%DLr     = 0.0D0

    !- for RODAS 3 only for version 1
    !IF( CHEMISTRY == 4) then	   
    !  ALLOCATE(spack(i)%DLr3  (1:maxblock_size,nspecies))    ;spack(i)%DLr3	= 0.0d0
    !  ALLOCATE(spack(i)%sc_p_4 (1:maxblock_size,nspecies))   ;spack(i)%sc_p_4  = 0.0d0
    !ENDIF
    
    !- 1D variables
    ALLOCATE(spack(i)%temp	       (1:maxblock_size))  ;spack(i)%temp      = 0.0D0  
    ALLOCATE(spack(i)%press	       (1:maxblock_size))  ;spack(i)%press     = 0.0D0 
    ALLOCATE(spack(i)%cosz	       (1:maxblock_size))  ;spack(i)%cosz      = 0.0d0
    ALLOCATE(spack(i)%att	       (1:maxblock_size))  ;spack(i)%att      = 0.0d0
    ALLOCATE(spack(i)%vapp	       (1:maxblock_size))  ;spack(i)%vapp      = 0.0d0
    ALLOCATE(spack(i)%volmol	       (1:maxblock_size))  ;spack(i)%volmol    = 0.0D0
    ALLOCATE(spack(i)%volmol_i         (1:maxblock_size))  ;spack(i)%volmol_i  = 0.0d0
    ALLOCATE(spack(i)%xlw	       (1:maxblock_size))  ;spack(i)%xlw       = 0.0d0
    ALLOCATE(spack(i)%err	       (1:maxblock_size))  ;spack(i)%err       = 0.0D0
    
     DO ii=1,maxblock_size
      ALLOCATE(spack_2d(ii,i)%DLmat(nspecies, nspecies ))    ;spack_2d(ii,i)%DLmat= 0.0d0
      ALLOCATE(spack_2d(ii,i)% DLb1(nspecies))  	     ;spack_2d(ii,i)% DLb1= 0.0D0
      ALLOCATE(spack_2d(ii,i)% DLb2(nspecies))  	     ;spack_2d(ii,i)% DLb2= 0.0d0
      ALLOCATE(spack_2d(ii,i)% DLk1(nspecies))  	     ;spack_2d(ii,i)% DLk1= 0.0d0
      ALLOCATE(spack_2d(ii,i)% DLk2(nspecies))  	     ;spack_2d(ii,i)% DLk2= 0.0D0
      !- for RODAS 3 only
      IF( CHEMISTRY == 4) then
    	ALLOCATE(spack_2d(ii,i)% DLb3(nspecies)) ;spack_2d(ii,i)% DLb3=   0.0D0 
    	ALLOCATE(spack_2d(ii,i)% DLb4(nspecies)) ;spack_2d(ii,i)% DLb4=   0.0d0
    	ALLOCATE(spack_2d(ii,i)% DLk3(nspecies)) ;spack_2d(ii,i)% DLk3=   0.0d0
    	ALLOCATE(spack_2d(ii,i)% DLk4(nspecies)) ;spack_2d(ii,i)% DLk4=   0.0D0
      ENDIF

      ENDDO
    ENDDO

    spack_alloc=.TRUE.

  END SUBROUTINE alloc_spack
  !-----------------------------------------------------------------
  !subroutine dealloc_spack(NOB_MEM)
  !implicit none
  !integer, intent(in) :: NOB_MEM
  !integer i,ii,NOB
  !  !IF(NOB_MEM==0) NOB=NOB_REAL
  !  !IF(NOB_MEM==1) NOB=1
  !
  !  do i=1,NOB_REAL
  !
  !    !if (associated(spack(i)%DLmat   )) deallocate(spack(i)%DLmat   )  	 
  !    !if (associated(spack(i)%DLmatLu )) deallocate(spack(i)%DLmatLu )      
  !    if (associated(spack(i)%DLdrdc  )) deallocate(spack(i)%DLdrdc  )	    
  !
  !    !2D Variables
  !
  !    if  (associated(spack(i)%jphoto))  deallocate (spack(i)%jphoto)
  !    if  (associated(spack(i)%rk    ))  deallocate (spack(i)%rk    )
  !    if  (associated(spack(i)%w     ))  deallocate (spack(i)%w     )
  !    if  (associated(spack(i)%sc_p  ))  deallocate (spack(i)%sc_p  )
  !
  !    if  (associated(spack(i)%sc_p_new))deallocate (spack(i)%sc_p_new )
  !    if  (associated(spack(i)%DLr     ))deallocate (spack(i)%DLr      )
  !    !if  (associated(spack(i)%DLk1    ))deallocate (spack(i)%DLk1     )
  !    !if  (associated(spack(i)%DLk2    ))deallocate (spack(i)%DLk2     )
  !    !if (associated(spack(i)%DLb1    ))deallocate (spack(i)%DLb1     )
  !    !if (associated(spack(i)%DLb2    ))deallocate (spack(i)%DLb2     )
  !					
  !    !1D variables
  !    if  (associated(spack(i)%temp   ))deallocate (spack(i)%temp   )
  !    if  (associated(spack(i)%press  ))deallocate (spack(i)%press  )
  !    if  (associated(spack(i)%cosz   ))deallocate (spack(i)%cosz   )
  !    if  (associated(spack(i)%vapp   ))deallocate (spack(i)%vapp   )
  !    if  (associated(spack(i)%volmol ))deallocate (spack(i)%volmol )
  !    if  (associated(spack(i)%xlw    ))deallocate (spack(i)%xlw    )
  !   
  !   do ii=1,maxblock_size
  !    if  (associated(spack_2d(ii,i)%DLmat)) deallocate (spack_2d(ii,i)%DLmat)
  !    if  (associated(spack_2d(ii,i)% DLb1)) deallocate (spack_2d(ii,i)% DLb1)
  !    if  (associated(spack_2d(ii,i)% DLb2)) deallocate (spack_2d(ii,i)% DLb2)
  !   enddo
  !
  !  enddo
  !
  ! !Deallocating  spaces to copy structure
  !  if(allocated(spack)) deallocate(spack)
  !
  !  if(allocated(spack_2d)) deallocate(spack_2d)
  !
  !end subroutine dealloc_spack

END MODULE mem_spack
