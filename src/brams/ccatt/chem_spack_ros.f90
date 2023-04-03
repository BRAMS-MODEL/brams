MODULE mod_chem_spack_ros

  USE mem_chem1, ONLY: &
       chem1_vars          ! Type

  USE mem_spack, ONLY: &
       spack_type,     &   ! Type
       spack_type_2d       ! Type

  !solver sparse
  USE solve_sparse, ONLY: &
       Solve_linear         ! Subroutine, Out: spack_2d(ijk,inob)%DLb1

  USE mod_chem_spack_jacdchemdc, ONLY: &
       jacdchemdc           ! Subroutine
  
  USE mod_chem_spack_kinetic, ONLY: &
       kinetic              ! Subroutine

  USE mod_chem_spack_fexchem, ONLY: &
       fexchem              ! Subroutine

  USE solve_sparse, ONLY: &
       Prepare,           & ! Subroutine
       Eliminate_0,       & ! Subroutine
       reserve              ! Subroutine


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: chem_ros, &
            get_always_zeros, &
            get_number_nonzeros

CONTAINS

  !========================================================================================
  SUBROUTINE chem_ros(m1,m2,m3,dtlt,pp,pi0,theta,rv,dn0,cosz,rcp,nspecies, &
                      nr,nr_photo, weight,PhotojMethod,jphoto,maxnspecies, &
                      nspecies_chem_transported,nspecies_chem_no_transported, &
                      transp_chem_index,no_transp_chem_index,chem1_g,nob,maxblock_size, &
                      block_end,indexk,indexi,indexj,kij_index, &
                      spack,spack_2d,get_non_zeros,N_DYN_CHEM,split_method, &
                      cp,cpor,p00,att)
  !========================================================================================
 

    IMPLICIT NONE
   
    INTEGER , INTENT(IN) :: m1
    INTEGER , INTENT(IN) :: m2
    INTEGER , INTENT(IN) :: m3

    ! mem_grid
    REAL             , INTENT(IN) :: dtlt

    ! mem_basic
    REAL             , INTENT(IN) :: pp(m1,m2,m3)
    REAL             , INTENT(IN) :: pi0(m1,m2,m3)
    REAL             , INTENT(IN) :: theta(m1,m2,m3)
    REAL             , INTENT(IN) :: rv(m1,m2,m3)
    REAL             , INTENT(IN) :: dn0(m1,m2,m3)

    ! mem_radiate
    REAL             , INTENT(IN) :: cosz(m2,m3)

    ! mem_micro
    REAL             , INTENT(IN) :: rcp(m1,m2,m3)

    ! chem1_list
    INTEGER          , INTENT(IN) :: nspecies
    INTEGER          , INTENT(IN) :: nr_photo
    INTEGER          , INTENT(IN) :: nr
    REAL             , INTENT(IN) :: weight(nspecies)
    CHARACTER(LEN=10), INTENT(IN) :: PhotojMethod

    ! FastJX
    DOUBLE PRECISION, INTENT(IN) :: jphoto(nr_photo,m1,m2,m3)

    ! mem_chem1
    INTEGER          , INTENT(IN) :: maxnspecies
    INTEGER          , INTENT(IN) :: nspecies_chem_transported
    INTEGER          , INTENT(IN) :: nspecies_chem_no_transported
    INTEGER          , INTENT(IN) :: transp_chem_index(maxnspecies)
    INTEGER          , INTENT(IN) :: no_transp_chem_index(maxnspecies)
    CHARACTER(LEN=20), INTENT(IN) :: split_method 
    INTEGER          , INTENT(IN) :: N_DYN_CHEM

    ! mem_chem1
    TYPE (chem1_vars), INTENT(INOUT) :: chem1_g(nspecies)

    ! spack_utils
    INTEGER         , INTENT(IN)    :: nob
    INTEGER         , INTENT(IN)    :: maxblock_size
    INTEGER         , INTENT(IN)    :: block_end(nob)
    INTEGER         , INTENT(IN)    :: indexk(maxblock_size,nob)
    INTEGER         , INTENT(IN)    :: indexi(maxblock_size,nob)
    INTEGER         , INTENT(IN)    :: indexj(maxblock_size,nob)
    INTEGER         , INTENT(IN)    :: kij_index(maxblock_size,nob)


    ! mem_spack
    TYPE(spack_type)   , INTENT(INOUT) :: spack(1)
    TYPE(spack_type_2d), INTENT(INOUT) :: spack_2d(maxblock_size,1)

    ! save variable
    LOGICAL , INTENT(INOUT) :: get_non_zeros

    ! rconstants
    REAL, INTENT(IN) :: cp
    REAL, INTENT(IN) :: cpor
    REAL, INTENT(IN) :: p00

    ! uv_atten
    REAL, INTENT(IN) :: att(m2,m3)

  DOUBLE PRECISION,parameter     :: pmar=28.96D0  
  DOUBLE PRECISION,parameter     :: Threshold1= 0.0D0		 
  DOUBLE PRECISION,parameter     :: Threshold2=-1.D-1		  
  
  DOUBLE PRECISION,parameter     :: Igamma = 1.0D0 + 1.0D0/1.4142135381698608 !(1.+ 1./SQRT(2.0)  )
  integer :: kij
  DOUBLE PRECISION dble_dtlt_i,fxc,dble_dtlt, Igamma_dtstep
  integer :: i,ijk,n,j,k,ispc,Ji,Jj,ii
  
  !integer,parameter ::  NOB_MEM=1   ! if 1 : alloc just one block (uses less memory)
  !                                  ! if 0 : alloc all blocks     (uses more memory)

  !- default number of allocatable blocks (re-use scratch arrays spack()% )
  integer,parameter :: inob=1 
  
  !logical           :: get_non_zeros = .false.
  
  integer,parameter :: number_cycling=1
  integer itime,k_,i_,j_,kij_
  
  !- controls memory allocation between differents grids

!  if(.not. spack_alloc) then 
!  !- Determining number of blocks and masks, allocating spack type
!   call allocIndex(block_size,N_DYN_CHEM) 
!   call alloc_spack(CHEMISTRY)
!  !- 
!  end if

  IF(.NOT. get_non_zeros) then
  !- Get always zero elements of jacobian matrix to save 
  !- computation time at solver:
   ! non-opt
   !get_non_zeros = .true.
   !call Prepare(nspecies)
   
   !  opt
   CALL get_always_zeros(nr_photo,nr,nspecies,get_non_zeros,spack(1)%rk(1,1:nr), &
   			 spack(1)%jphoto(1,1:nr_photo),p00)

  ENDIF
  
  dble_dtlt  = dble(    dtlt*N_DYN_CHEM/float(number_cycling))
  dble_dtlt_i= 1.0D0/dble_dtlt 
  
 DO i=1,nob
   
    !- copying structure from input to internal

    DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
    
      k_=indexk(ijk,i) !index_g(ngrid)%indexk(ijk,i) ! k brams index
      i_=indexi(ijk,i) !index_g(ngrid)%indexi(ijk,i) ! i brams index
      j_=indexj(ijk,i) !index_g(ngrid)%indexj(ijk,i) ! j brams index
    
      !- Exner function divided by Cp
       spack(inob)%press(ijk) = ( pp (k_,i_,j_) +   pi0(k_,i_,j_) )/cp
      
      !- Air temperature (K)			      			      
       spack(inob)%temp(ijk)= theta(k_,i_,j_) * spack(inob)%press(ijk)
    	 
      !- transform from Exner function to pressure
       spack(inob)%press(ijk)= spack(inob)%press(ijk)**cpor*p00      

      !- Water vapor 
       spack(inob)%vapp(ijk) =  rv(k_,i_,j_)!&
                             !* dn0(k_,i_,j_)*6.02D20/18.D0

      !- volmol eh usado para converter de ppbm para molecule/cm^3 (fator fxc incluido)
      ! spack(inob)%volmol(ijk)=(spack(inob)%press(ijk))/(8.314e0*spack(inob)%temp(ijk))
      ! fxc = (6.02e23*1e-15*spack(inob)%volmol(ijk))*pmar
      !
       spack(inob)%volmol(ijk)=(6.02D23*1D-15*pmar)*(spack(inob)%press(ijk))/(8.314D0*spack(inob)%temp(ijk))
      
      !- inverse of volmol to get back to ppbm at the end of routine
       spack(inob)%volmol_I(ijk)=1.0d0 / spack(inob)%volmol(ijk)
               
      !- get liquid water content
       spack(inob)%xlw(ijk) =  rcp(k_,i_,j_) * dn0(k_,i_,j_)*1.D-3
             
      !- convert from brams chem (ppbm) arrays to spack (molec/cm3)       
      !- no transported species section
      DO ispc=1,nspecies_chem_no_transported
        
        !- map the species to NO transported ones
        n=no_transp_chem_index(ispc)
        
        !- initialize no-transported species (don't need to convert, because these
        !- species are already saved using molecule/cm^3 units)
        spack(inob)%sc_p(ijk,n) = chem1_g(n)%sc_p(k_,i_,j_) 
      END DO
       
    END DO
    	   
    !- convert from brams chem (ppbm) arrays to spack (molec/cm3)
    !- transported species section
    IF(split_method == 'PARALLEL' .and. N_DYN_CHEM > 1) then 
      DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
         
	 kij_ =kij_index(ijk,i) !index_g(ngrid)%kij_index(ijk,i)!- kij brams tendency index
         k_=indexk(ijk,i) !index_g(ngrid)%indexk(ijk,i) ! k brams index
         i_=indexi(ijk,i) !index_g(ngrid)%indexi(ijk,i) ! i brams index
         j_=indexj(ijk,i) !index_g(ngrid)%indexj(ijk,i) ! j brams index
   
         DO ispc=1,nspecies_chem_transported  
          
           !- map the species to transported ones
           n=transp_chem_index(ispc)
           
           !- conversion from ppbm to molecule/cm^3           
	   !- get back the concentration at the begin of the chemistry timestep integration:
	   
           spack(inob)%sc_p(ijk,n) =  (  chem1_g(n)%sc_p(k_,i_,j_) -                 &  ! updated mixing ratio 
	                                 chem1_g(n)%sc_t_dyn(kij_)*N_DYN_CHEM*dtlt  )&  ! accumulated tendency
        			       * spack(inob)%volmol(ijk)/weight(n)
           !- for testing 		  
           !spack(inob)%sc_p(ijk,n)   = max(0.,spack(inob)%sc_p(ijk,n))
           !spack(inob)%sc_p_4(ijk,n) =(chem1_g(n)%sc_p(k_,i_,j_) - beta*chem1_g(n)%sc_t_dyn(kij_)*N_DYN_CHEM*dtlt)& !dtlt eh o dinamico
           !				* spack(inob)%volmol(ijk)/weight(n)
         END DO 	 
      END DO	   
    ELSE
       
      DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
         
	 kij_ =kij_index(ijk,i) !index_g(ngrid)%kij_index(ijk,i)!- kij brams tendency index
         k_=indexk(ijk,i) !index_g(ngrid)%indexk(ijk,i) ! k brams index
         i_=indexi(ijk,i) !index_g(ngrid)%indexi(ijk,i) ! i brams index
         j_=indexj(ijk,i) !index_g(ngrid)%indexj(ijk,i) ! j brams index

         DO ispc=1,nspecies_chem_transported  
          
           !- map the species to transported ones
           n=transp_chem_index(ispc)
           
           !- conversion from ppbm to molecule/cm^3
           spack(inob)%sc_p(ijk,n) = chem1_g(n)%sc_p(k_,i_,j_) * spack(inob)%volmol(ijk)/weight(n)
           spack(inob)%sc_p(ijk,n) = max(0.D0,spack(inob)%sc_p(ijk,n))

         END DO 	 
      END DO	   
    ENDIF

    !- Photolysis section
    IF(trim(PhotojMethod) == 'FAST-JX' .or. trim(PhotojMethod) == 'FAST-TUV' ) THEN
     
     DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
        
	k_=indexk(ijk,i) !index_g(ngrid)%indexk(ijk,i) ! k brams index
        i_=indexi(ijk,i) !index_g(ngrid)%indexi(ijk,i) ! i brams index
        j_=indexj(ijk,i) !index_g(ngrid)%indexj(ijk,i) ! j brams index
        
	DO n=1,nr_photo
          spack(inob)%jphoto(ijk,n)= jphoto(n,k_,i_,j_)!fast_JX_g(ngrid)%jphoto(n,k_,i_,j_)
        END DO
     
      ENDDO
     
     ELSEIF(trim(PhotojMethod) == 'LUT') then
      
      DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
        
        i_=indexi(ijk,i) !index_g(ngrid)%indexi(ijk,i) ! i brams index
        j_=indexj(ijk,i) !index_g(ngrid)%indexj(ijk,i) ! j brams index
        
	!- UV  attenuation un function AOT
        spack(inob)%att(ijk)= att(i_,j_) 
        
	!- get zenital angle (for LUT PhotojMethod)
        spack(inob)%cosz(ijk)=  cosz(i_,j_)

      ENDDO     
    ENDIF	   

 
!!!!   END DO


  !- Call kinetic
    
   
!!!!       DO i=1,nob
    	CALL kinetic(nr_photo,spack(inob)%Jphoto  &
    			     ,spack(inob)%rk      &
  			     ,spack(inob)%temp    &
    			     ,spack(inob)%vapp    &
    			     ,spack(inob)%Press   &
    			     ,spack(inob)%cosz    &
    			     ,spack(inob)%att     &
			     ,1,block_end(i)      & !index_g(ngrid)%block_end(i)
                             ,maxblock_size,nr)
!!!!       END DO


!- ROSCHEM-------------------------------------------------------------------------------
!- kml Inicio do equivalente a roschem

!     1 - First step
       
!!!!  DO i=1,nob
  do itime=1,number_cycling



!-   Compute the Jacobian (DLRDC).
     CALL jacdchemdc (spack(inob)%sc_p  &
	             ,spack(inob)%rk	& 
      	             ,spack(inob)%DLdrdc& ! Jacobian matrix
		     ,nspecies,1,block_end(i) & !index_g(ngrid)%block_end(i)
                     ,maxblock_size,nr)
    

!-   Compute matrix (1-Igamma*dt*DLRDC) 
     
     Igamma_dtstep = Igamma * dble_dtlt
     
     DO Jj=1,nspecies
   	DO Ji=1,nspecies
     	  DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)

   	   spack_2d(ijk,inob)%DLmat(Ji,Jj) = -  Igamma_dtstep * spack(inob)%DLdrdc(ijk,Ji,Jj)

   	  ENDDO
     	ENDDO
     ENDDO
     DO Jj=1,nspecies
     	DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
     	     spack_2d(ijk,inob)%DLmat(Jj,Jj) = 1.0D0 + spack_2d(ijk,inob)%DLmat(Jj,Jj)
   	ENDDO
     ENDDO
   
!-   Compute chemical production terms at initial time 
     CALL fexchem (spack(inob)%sc_p &
        	  ,spack(inob)%rk   &
        	  ,spack(inob)%DLr  & !production term  
        	  ,nspecies,1,block_end(i) & !index_g(ngrid)%block_end(i)
                  ,maxblock_size,nr)
 
     DO Ji=1,nspecies
       DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
	   spack_2d(ijk,inob)%DLb1(Ji) = spack(inob)%DLr(ijk,Ji)
       ENDDO
     ENDDO

!-   Compute DLb1 by Solving (1-Igamma*dt*DLRDC) DLb1=DLR 

     !-for solver  LU =============================
     !DO ijk=1,index_g(ngrid)%block_end(i)
     ! CALL solve(spack_2d(ijk,inob)%DLmat,spack_2d(ijk,inob)%DLb1,nspecies,keep=.true.)
     !ENDDO

     !- solver sparse
      DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
        CALL Solve_linear(nspecies,spack_2d(ijk,inob)%DLmat&
	                  	  ,spack_2d(ijk,inob)%DLb1 &
			   	  ,spack_2d(ijk,inob)%DLb1 )!&
			   	  !,i)  	        
      ENDDO
!
!---------------------------------------------------------------------------------------------------------------------
!    2- Second step
!    Compute first-order approximation (sc_p_new)

     DO Ji=1,nspecies
         DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)

    	    spack(inob)%sc_p_new(ijk,Ji) = spack(inob)%sc_p(ijk,Ji) + dble_dtlt * spack_2d(ijk,inob)%DLb1(Ji)
    	     	  
     	    IF (spack(inob)%sc_p_new(ijk,Ji) .LT. Threshold1) THEN
     	        
		spack(inob)%sc_p_new(ijk,Ji) = Threshold1
  
     	        !spack(inob)%DLk1(ijk,Ji)   = (spack(inob)%sc_p_new(ijk,Ji) - spack(inob)%sc_p(ijk,Ji)) / dtlt ! original
      	        spack_2d(ijk,inob)%DLb1(Ji) = (spack(inob)%sc_p_new(ijk,Ji) - spack(inob)%sc_p(ijk,Ji)) * dble_dtlt_i 
                
		!print*,'cliping1=',spc_name(Ji),spack(inob)%sc_p_new(ijk,Ji),spack_2d(ijk,inob)%DLb1(Ji)
     	  
	    !ELSE !original

   	      !spack(inob)%DLk1(ijk,Ji) = spack_2d(ijk,inob)%DLb1(Ji)!original

     	    ENDIF

        !if(Ji==XO2.and.ijk==index_g(ngrid)%block_end(i).and.i==1)&
	! print*,'BIS=',maxval(spack(inob)%sc_p_new(1:index_g(ngrid)%block_end(i),ji)),i
       ENDDO
     ENDDO


!-   Compute chemical production terms (DLr) at final time with the first-order approximation
                     
     CALL fexchem (spack(inob)%sc_p_new  &
     		  ,spack(inob)%rk        &
	  	  ,spack(inob)%DLr       &    
     		  ,nspecies,1,block_end(i) & !index_g(ngrid)%block_end(i)
                  ,maxblock_size,nr)
  
!-    Compute DLB2    
      DO Ji=1,nspecies       
        DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
          !spack_2d(ijk,inob)%DLb2(Ji) = spack(inob)%DLr(ijk,Ji) - 2.0 * spack(inob)%DLk1(ijk,Ji) !original
           spack_2d(ijk,inob)%DLb2(Ji) = spack(inob)%DLr(ijk,Ji) - 2.0D0 * spack_2d(ijk,inob)%DLb1(Ji)  
       ENDDO
      ENDDO
	 
!-   Compute DLb2 by solving (1-Igamma*dt*DLRDC) DLB2=DLK1      
     !-for solver  LU	  =============================
     !DO ijk=1,index_g(ngrid)%block_end(i)
     !  CALL solve(spack_2d(ijk,inob)%DLmat,spack_2d(ijk,inob)%DLb2,nspecies,keep=.true.)
     !ENDDO
     !- for solver sparse =============================
      DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
        CALL Solve_linear(nspecies,spack_2d(ijk,inob)%DLmat&
	                  	  ,spack_2d(ijk,inob)%DLb2 &
			   	  ,spack_2d(ijk,inob)%DLb2 )!&
			   	  !,i) 
      ENDDO			  
!-----------------------------------------------------------------------
!    3- Compute concentrations at final time (optimized scheme):
!    sc_p + (3*DLK1 +DLK2) * dt/2 - original
!    sc_p + (3*DLB1 +DLB2) * dt/2 - opt

     DO Ji=1,nspecies
     
       DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
         ! spack(inob)%DLk2(ijk,Ji) = spack_2d(ijk,inob)%DLb2(Ji)  ! original
	                            ! spack(inob)%DLk1(ijk,Ji) +  &! original
                                    ! spack(inob)%DLk2(ijk,Ji)	   ! original
      
          spack(inob)%sc_p(ijk,Ji)= spack(inob)%sc_p(ijk,Ji) + dble_dtlt * &
        			  ( 1.5D0 * spack_2d(ijk,inob)%DLb1(Ji)  + &
        			    0.5D0 * spack_2d(ijk,inob)%DLb2(Ji)  )  
 
          
          spack(inob)%sc_p(ijk,Ji) = max(spack(inob)%sc_p(ijk,Ji),Threshold1)
          
	  
	  
	 ! if(Ji==XO2.and.ijk==index_g(ngrid)%block_end(i).and.i==1) print*,'FIM=',maxval(spack(inob)%sc_p(1:index_g(ngrid)%block_end(i),ji)),i

	  
         ! IF (spack(inob)%sc_p(ijk,Ji) .LT. Threshold1) THEN
         !		
         ! !IF(spack(inob)%sc_p(ijk,Ji).LT.Threshold2) THEN
         ! !  Write (*,*) 'CLIP CHEMISTRY',spack(inob)%sc_p(ijk,Ji), &
         ! !   Ji,indexi(ijk,n),indexj(ijk,n),indexk(ijk,n) !ijk=ijk e n=nob
         ! !  
         ! !ENDIF
         ! 
         ! spack(inob)%sc_p(ijk,Ji) = Threshold1
         ! ENDIF
       ENDDO
         
     ENDDO


  ENDDO ! time-spliting 

!!!!!!!!!!!!!!!!  ENDDO


!!!!!!!!!!!!!!!!  DO i=1,nob
!500 continue 

  !-----------------------------------------------------------------------
  !- Restoring species loss/prod from internal to brams strucuture

  !- transported species section

    if(split_method == 'PARALLEL' .and.  N_DYN_CHEM > 1) then 

       DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)

         kij_ =kij_index(ijk,i) !index_g(ngrid)%kij_index(ijk,i)
         k_   =indexk(ijk,i) !index_g(ngrid)%indexk(ijk,i)   
         i_   =indexi(ijk,i) !index_g(ngrid)%indexi(ijk,i)   
         j_   =indexj(ijk,i) !index_g(ngrid)%indexj(ijk,i)   

         DO ispc=1,nspecies_chem_transported 
	  
            !- map the species to transported ones
            n=transp_chem_index(ispc)
            
            !- include the chemical tendency at total tendency (convert to unit: ppbm/s)
            chem1_g(n)%sc_p(k_,i_,j_ ) =  chem1_g(n)%sc_t_dyn(kij_)*N_DYN_CHEM*dtlt +& 
        				       spack(inob)%sc_p(ijk,n)*weight(n)*spack(inob)%volmol_I(ijk)	
            
            chem1_g(n)%sc_p(k_,i_,j_ ) = max(0.,chem1_g(n)%sc_p(k_,i_,j_ ))

        END DO
      ENDDO
    
    ELSEIF(split_method == 'PARALLEL' .and.  N_DYN_CHEM == 1) then
        
       dble_dtlt_i=1.0d0 / dble( dtlt) 
       DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)

         kij_ =kij_index(ijk,i) !index_g(ngrid)%kij_index(ijk,i)
         k_   =indexk(ijk,i) !index_g(ngrid)%indexk(ijk,i)   
         i_   =indexi(ijk,i) !index_g(ngrid)%indexi(ijk,i)   
         j_   =indexj(ijk,i) !index_g(ngrid)%indexj(ijk,i)   

         DO ispc=1,nspecies_chem_transported  
       
            !- map the species to transported ones
            n=transp_chem_index(ispc)
            
            !- include the chemical tendency at total tendency (convert to unit: ppbm/s)
            chem1_g(n)%sc_t(kij_) = chem1_g(n)%sc_t(kij_)			    +  &! previous tendency
       			       
      			     ( spack(inob)%sc_p(ijk,n)*weight(n)*spack(inob)%volmol_I(ijk)  -  &! new mixing ratio
       			       chem1_g(n)%sc_p(k_,i_,j_ ) 		  )	       &! old mixing ratio
       			     * dble_dtlt_i							! inverse of timestep
         END DO
       END DO

     ELSE
       
       DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)

         kij_ =kij_index(ijk,i) !index_g(ngrid)%kij_index(ijk,i)
         k_   =indexk(ijk,i) !index_g(ngrid)%indexk(ijk,i)   
         i_   =indexi(ijk,i) !index_g(ngrid)%indexi(ijk,i)   
         j_   =indexj(ijk,i) !index_g(ngrid)%indexj(ijk,i)   
         
	 DO ispc=1,nspecies_chem_transported  
       
            !- map the species to transported ones
            n=transp_chem_index(ispc)

  	    chem1_g(n)%sc_p(k_,i_,j_ ) = spack(inob)%sc_p(ijk,n)*weight(n)*spack(inob)%volmol_I(ijk)   
	    chem1_g(n)%sc_p(k_,i_,j_ ) = max(0.,chem1_g(n)%sc_p(k_,i_,j_ ))  
        END DO
       END DO

     ENDIF


  !- no transported species section
     DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)

         kij_ =kij_index(ijk,i) !index_g(ngrid)%kij_index(ijk,i)
         k_   =indexk(ijk,i) !index_g(ngrid)%indexk(ijk,i)   
         i_   =indexi(ijk,i) !index_g(ngrid)%indexi(ijk,i)   
         j_   =indexj(ijk,i) !index_g(ngrid)%indexj(ijk,i)   
         DO ispc=1,nspecies_chem_no_transported       
        
	  !- map the species to no transported ones
          n=no_transp_chem_index(ispc)
        
	  !- save no-transported species (keep current unit : molec/cm3)
          chem1_g(n)%sc_p(k_,i_,j_ )  = max(0., real ( spack(inob)%sc_p(ijk,n) ))   
         END DO

     END DO

 END DO ! enddo loop over all blocks
END SUBROUTINE chem_ros
!--------------------------------------------------------------------------  
  SUBROUTINE get_always_zeros(jppj,nr,nspecies,get_non_zeros,rk,jphoto,p00)

    INTEGER          , INTENT(IN)    :: jppj
    INTEGER          , INTENT(IN)    :: nr
    INTEGER          , INTENT(IN)    :: nspecies
    LOGICAL          , INTENT(INOUT) :: get_non_zeros
    DOUBLE PRECISION , INTENT(INOUT) :: rk(nr)
    DOUBLE PRECISION , INTENT(INOUT) :: jphoto(jppj)
    REAL             , INTENT(IN)    :: p00

    DOUBLE PRECISION ,DIMENSION(nspecies,nspecies) :: def_non_zeros 
    DOUBLE PRECISION ,DIMENSION(nspecies) :: sc_p
    DOUBLE PRECISION  :: xlw,vapp(1),cosz(1),temp(1),press(1),att(1)
    INTEGER :: i

    get_non_zeros = .TRUE.

    CALL Prepare(nspecies)

    jphoto(:) = 2.3333331D0
    xlw	  = 1.D0
    vapp(:) = 1.D15
    cosz(:) = 1.D0
    att(:)  = 1.0D0
    temp(:) = 273.15D0
    press(:)= DBLE(p00)
    sc_p    = 1.D15 ! dummy concentration to get the maximum number
                    ! of possible non zero elements

    call kinetic(jppj,jphoto	 &
     		     ,rk(1:nr)   &
     		     ,temp	 &
     		     ,vapp	 &
     		     ,press	 &
     		     ,cosz	 &
     		     ,att,1,1,1,nr  )

    def_non_zeros = 0.D0

    CALL jacdchemdc(sc_p,rk(1:nr),def_non_zeros, & ! Jacobian matrix
                    nspecies,1,1,1,nr )

    DO i=1,nspecies 
       def_non_zeros(i,i)=1.D0+ def_non_zeros(i,i)
    ENDDO

    CALL Eliminate_0(def_non_zeros,nspecies) 
    CALL reserve()

  END SUBROUTINE get_always_zeros

!--------------------------------------------------------------------------  
  SUBROUTINE get_number_nonzeros(jppj,nr,nspecies,rk,jphoto,p00,nonzeros)

    INTEGER          , INTENT(IN)    :: jppj
    INTEGER          , INTENT(IN)    :: nr
    INTEGER          , INTENT(IN)    :: nspecies
    DOUBLE PRECISION , INTENT(INOUT) :: rk(nr)
    DOUBLE PRECISION , INTENT(INOUT) :: jphoto(jppj)
    REAL             , INTENT(IN)    :: p00
    INTEGER          , INTENT(INOUT) :: nonzeros

    DOUBLE PRECISION ,DIMENSION(nspecies,nspecies) :: def_non_zeros 
    DOUBLE PRECISION ,DIMENSION(nspecies) :: sc_p
    DOUBLE PRECISION  :: xlw,vapp(1),cosz(1),temp(1),press(1),att(1)
    INTEGER :: i,Ji,Jj

    !get_non_zeros = .TRUE.

    !CALL Prepare(nspecies)

    jphoto(:) = 2.3333331D0
    xlw	  = 1.D0
    vapp(:) = 1.D15
    cosz(:) = 1.D0
    att(:)  = 1.0D0
    temp(:) = 273.15D0
    press(:)= DBLE(p00)
    sc_p    = 1.D15 ! dummy concentration to get the maximum number
                    ! of possible non zero elements

    call kinetic(jppj,jphoto	 &
     		     ,rk(1:nr)   &
     		     ,temp	 &
     		     ,vapp	 &
     		     ,press	 &
     		     ,cosz	 &
     		     ,att,1,1,1,nr  )

    def_non_zeros = 0.D0

    CALL jacdchemdc(sc_p,rk(1:nr),def_non_zeros, & ! Jacobian matrix
                    nspecies,1,1,1,nr )

    nonzeros = 0
    DO Jj=1,nspecies
       def_non_zeros(Jj,Jj)=1.D0+ def_non_zeros(Jj,Jj)
       DO Ji=1,nspecies
          if (def_non_zeros(Ji,Jj) .ne. 0.D0) then
             nonzeros = nonzeros + 1
          endif
       ENDDO
    ENDDO

!    DO i=1,nspecies 
!       def_non_zeros(i,i)=1.D0+ def_non_zeros(i,i)
!    ENDDO

    !CALL Eliminate_0(def_non_zeros,nspecies) 
    !CALL reserve()

  END SUBROUTINE get_number_nonzeros


END MODULE mod_chem_spack_ros

  !-----------------------------------------------------------------------
  ! (DMK) Problema: se a subrotina get_always_zeros for inserida dentro 
  !       do modulo mod_chem_spack_ros, a saida de chemistry_driver deixa 
  !       de bater binario. O problema ocorre apenas na versao standalone.
  !       A versao limpa, acoplada esta livre deste problema. Acredito que
  !       esteja faltando algum dado de entrada.
  !-----------------------------------------------------------------------

!!$  SUBROUTINE get_always_zeros(jppj,nr,nspecies,get_non_zeros,rk,jphoto,p00)
!!$
!!$    USE mod_chem_spack_jacdchemdc, ONLY: &
!!$         jacdchemdc   ! Subroutine
!!$  
!!$    USE mod_chem_spack_kinetic, ONLY: &
!!$         kinetic      ! Subroutine
!!$
!!$    USE solve_sparse, ONLY: &
!!$         Prepare,           & ! Subroutine
!!$         Eliminate_0,       & ! Subroutine
!!$         reserve              ! Subroutine
!!$
!!$    IMPLICIT NONE
!!$
!!$    INTEGER          , INTENT(IN)    :: jppj
!!$    INTEGER          , INTENT(IN)    :: nr
!!$    INTEGER          , INTENT(IN)    :: nspecies
!!$    LOGICAL          , INTENT(INOUT) :: get_non_zeros
!!$    DOUBLE PRECISION , INTENT(INOUT) :: rk(nr)
!!$    DOUBLE PRECISION , INTENT(INOUT) :: jphoto(jppj)
!!$    REAL             , INTENT(IN)    :: p00
!!$
!!$    DOUBLE PRECISION ,DIMENSION(nspecies,nspecies) :: def_non_zeros 
!!$    DOUBLE PRECISION ,DIMENSION(nspecies) :: sc_p
!!$    DOUBLE PRECISION  :: xlw,vapp(1),cosz(1),temp(1),press(1)
!!$    INTEGER :: i
!!$
!!$    get_non_zeros = .TRUE.
!!$
!!$    CALL Prepare(nspecies)
!!$
!!$    jphoto(:) = 1.D0
!!$    xlw	  = 1.D0
!!$    vapp(:) = 1.D15
!!$    cosz(:) = 1.D0
!!$    temp(:) = 273.15D0
!!$    press(:)= DBLE(p00)
!!$    sc_p    = 1.D15 ! dummy concentration to get the maximum number
!!$    ! of possible non zero elements
!!$
!!$    CALL kinetic(jppj,jphoto,rk(1:nr),temp,vapp,press,1,1,1,nr)
!!$
!!$    def_non_zeros = 0.D0
!!$
!!$    CALL jacdchemdc(sc_p,rk(1:nr),def_non_zeros, & ! Jacobian matrix
!!$                    nspecies,1,1,1,nr )
!!$
!!$    DO i=1,nspecies 
!!$       def_non_zeros(i,i)=1.D0+ def_non_zeros(i,i)
!!$    ENDDO
!!$
!!$    CALL Eliminate_0(def_non_zeros,nspecies) 
!!$    CALL reserve()
!!$
!!$  END SUBROUTINE get_always_zeros
