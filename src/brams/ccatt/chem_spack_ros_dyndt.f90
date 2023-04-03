MODULE mod_chem_spack_ros_dyndt


  USE mem_chem1, ONLY: &
       chem1_vars        ! Type

  USE mem_spack, ONLY: &
       spack_type,     & ! Type
       spack_type_2d     ! Type

  !solver sparse
  USE solve_sparse, ONLY: &
       Solve_linear      ! Subroutine, OUT: spack_2d(ijk,inob)%DLk1

  USE mod_chem_spack_ros, ONLY: &
       get_always_zeros  ! Subroutine

  USE mod_chem_spack_jacdchemdc, ONLY: &
       jacdchemdc        ! Subroutine

  USE mod_chem_spack_kinetic, ONLY: &
       kinetic           ! Subroutine

  USE mod_chem_spack_fexchem, ONLY: &
       fexchem           ! Subroutine


  IMPLICIT NONE


  PRIVATE


  PUBLIC chem_ros_dyndt


CONTAINS
  !========================================================================================
  SUBROUTINE chem_ros_dyndt(m1,m2,m3,dtlt,pp,pi0,theta,rv,dn0,cosz,rcp, &
                            nspecies,nr,nr_photo,weight,PhotojMethod,maxnspecies, &
                            nspecies_chem_transported,nspecies_chem_no_transported, &
                            transp_chem_index,no_transp_chem_index,chem1_g,nob,maxblock_size, &
                            block_end,indexk,indexi,indexj,kij_index,last_accepted_dt, &
                            spack,spack_2d,atol,rtol,get_non_zeros,&
			    N_DYN_CHEM,split_method,cp,cpor,p00,jphoto,att)
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
    INTEGER         ,   INTENT(IN)    :: nob
    INTEGER         ,   INTENT(IN)    :: maxblock_size
    INTEGER         ,   INTENT(IN)    :: block_end(nob)
    INTEGER         ,   INTENT(IN)    :: indexk(maxblock_size,nob)
    INTEGER         ,   INTENT(IN)    :: indexi(maxblock_size,nob)
    INTEGER         ,   INTENT(IN)    :: indexj(maxblock_size,nob)
    INTEGER         ,   INTENT(IN)    :: kij_index(maxblock_size,nob)
    DOUBLE PRECISION,   INTENT(INOUT) :: last_accepted_dt(nob) 

    ! mem_spack
    TYPE(spack_type)   , INTENT(INOUT) :: spack(1)
    TYPE(spack_type_2d), INTENT(INOUT) :: spack_2d(maxblock_size,1)
    DOUBLE PRECISION   , INTENT(IN)    :: atol(nspecies)
    DOUBLE PRECISION   , INTENT(IN)    :: rtol(nspecies)

    ! save variable
    LOGICAL , INTENT(INOUT) :: get_non_zeros

    ! rconstants
    REAL, INTENT(IN) :: cp
    REAL, INTENT(IN) :: cpor
    REAL, INTENT(IN) :: p00

    ! Fast_JX
    DOUBLE PRECISION, INTENT(IN) :: jphoto(nr_photo,m1,m2,m3)

    ! uv_atten
    REAL, INTENT(IN) :: att(m2,m3)

  
!  INTEGER,PARAMETER :: block_size=2*32
!  INTEGER,PARAMETER :: block_size=4
!  INTEGER,PARAMETER :: block_size=1

  DOUBLE PRECISION,parameter     :: pmar=28.96D0  
  DOUBLE PRECISION,parameter     :: Threshold1= 0.0D0		 
  DOUBLE PRECISION,parameter     :: Threshold2=-1.D-1		  
  
  DOUBLE PRECISION,parameter     :: Igamma = 1.0D0 + 1.0D0/1.4142135381698608 !(1.+ 1./SQRT(2.0)  )
  DOUBLE PRECISION,parameter     :: Igamma_i = 1.0D0/Igamma
  integer :: kij
  DOUBLE PRECISION dble_dtlt_i,fxc, Igamma_dtstep,dt_chem,dt_chem_i
  integer :: i,ijk,n,j,k,ispc,Ji,Jj, k_,i_,j_,kij_
  
  !integer,parameter ::  NOB_MEM=1   ! if 1 : alloc just one block (uses less memory)
  !                                  ! if 0 : alloc all blocks     (uses more memory)

  !- default number of allocatable blocks (re-use scratch arrays spack()% )
  integer,parameter :: inob=1 
  
  !logical           :: get_non_zeros = .false.

  !- parameters for dynamic timestep control
  integer all_accepted
  double precision :: dt_min,dt_max,dt_actual,dt_new
  double precision :: time_f,time_c
  double precision, parameter :: FacMin =0.2d0 ! lower bound on step decrease factor (default=0.2)
  double precision, parameter :: FacMax =6.0D0 ! upper bound on step increase factor (default=6)
  double precision, parameter :: FacRej =0.1D0 ! step decrease factor after multiple rejections
  double precision, parameter :: FacSafe=0.9D0 ! step by which the new step is slightly smaller
                                               ! than the predicted value  (default=0.9)
  double precision, parameter :: uround  = 1.d-15,   elo = 2.0D0, Roundoff = 1.d-8
  double precision Fac, Tol,err1,max_err1
  
  
!  if(.not. spack_alloc) then 
!  !- Determining number of blocks and masks, allocating spack type
!   call allocIndex(block_size,N_DYN_CHEM) 
!   call alloc_spack(CHEMISTRY)
!  !- 
!  end if

  IF(.NOT. get_non_zeros) then
  !- Get always zero elements of Jacobian matrix to save 
  !- computational time at solver:
  ! non-opt
  !get_non_zeros = .true.
  !call Prepare(nspecies)
   
   !  opt
   call get_always_zeros(nr_photo,nr,nspecies,get_non_zeros,spack(1)%rk(1,1:nr), &
                             spack(1)%jphoto(1,1:nr_photo),p00)
  ENDIF
  
  DO i=1,nob !- loop over all blocks
    
    !- copying structure from input to internal

    DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
    
      k_=indexk(ijk,i) !index_g(ngrid)%indexk(ijk,i) ! k brams index
      i_=indexi(ijk,i) !index_g(ngrid)%indexi(ijk,i) ! i brams index
      j_=indexj(ijk,i) !index_g(ngrid)%indexj(ijk,i) ! j brams index
    
      !- Exner function divided by Cp
       spack(inob)%press(ijk) = ( pp (k_,i_,j_) + pi0(k_,i_,j_) )/cp
      
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
          spack(inob)%jphoto(ijk,n)= jphoto(n,k_,i_,j_)!fast_JX_g%jphoto(n,k_,i_,j_)
        END DO
     
      ENDDO
     
     ELSEIF(trim(PhotojMethod) == 'LUT') then
      
      DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
        
        i_=indexi(ijk,i) !index_g(ngrid)%indexi(ijk,i) ! i brams index
        j_=indexj(ijk,i) !index_g(ngrid)%indexj(ijk,i) ! j brams index
        
	!- UV  attenuation un function AOT
        spack(inob)%att(ijk)= att(i_,j_) !uv_atten_g(ngrid)%att(i_,j_) 
        
	!- get zenital angle (for LUT PhotojMethod)
        spack(inob)%cosz(ijk)=  cosz(i_,j_)

      ENDDO     
    ENDIF	   

 
!- compute kinetical and photochemical reactions (array: spack(inob)%rk)
    	CALL kinetic(nr_photo,spack(inob)%Jphoto  &
    			     ,spack(inob)%rk      &
  			     ,spack(inob)%temp    &
    			     ,spack(inob)%vapp    &
    			     ,spack(inob)%Press   &
    			     ,spack(inob)%cosz    &
    			     ,spack(inob)%att     &
			     ,1,block_end(i)      & !index_g(ngrid)%block_end(i)
                             ,maxblock_size,nr)


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     !tmp
!
!  !- srf: to use these arrays, be sure NA_EXTRAD3D RAMSIN is at least 8
!    DO ijk=1,index_g(ngrid)%block_end(i)
!	k_=index_g(ngrid)%indexk(ijk,i) ! k brams index
!        i_=index_g(ngrid)%indexi(ijk,i) ! i brams index
!        j_=index_g(ngrid)%indexj(ijk,i) ! j brams index
!        
!        extra3d(7,ngrid)%d3(k_,i_,j_)=spack(inob)%rk(ijk,1)
!        extra3d(8,ngrid)%d3(k_,i_,j_)=spack(inob)%rk(ijk,3)+spack(inob)%rk(ijk,2)
!
!    ENDDO
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     !tmp




!- ROSCHEM-------------------------------------------------------------------------------
!- kml Inicio do equivalente a roschem
!- srf- including dyn timestep

  dt_chem   = last_accepted_dt(i) !index_g(ngrid)%last_accepted_dt(i) ! dble(dtlt)
  dt_min   = max(1.0D0,1.D-2*dble(dtlt*N_DYN_CHEM))
  dt_max   = dble(dtlt*N_DYN_CHEM) 
  dt_new   = 0.0d0

  time_c   = 0.0d0
  time_f   = dble(dtlt*N_DYN_CHEM) 
  
  run_until_integr_ends:   DO WHILE (time_c + Roundoff < time_f)

!-   Compute the Jacobian (DLRDC).
     CALL jacdchemdc (spack(inob)%sc_p  &
	             ,spack(inob)%rk	& 
      	             ,spack(inob)%DLdrdc& ! Jacobian matrix
		     ,nspecies,1,block_end(i) & !index_g(ngrid)%block_end(i)
                     ,maxblock_size,nr)

!-   Compute chemical net production terms (actually, P-L) at initial time (array spack(inob)%DLr)
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


!-   compute matrix (1/(Igamma*dt)  - Jacobian)  where Jacobian = DLdrdc       
!-   fill DLMAT with non-diagonal (fixed in the timestep) Jacobian
     DO Jj=1,nspecies
   	DO Ji=1,nspecies
     	  DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)

   	   spack_2d(ijk,inob)%DLmat(Ji,Jj) = - spack(inob)%DLdrdc(ijk,Ji,Jj)

   	  ENDDO
     	ENDDO
     ENDDO
     
     
     UntilAccepted: DO

!- fill DLMAT with diagonal Jacobian (which changes in time, because of dt_chem)
     Igamma_dtstep = 1.0D0/(Igamma * dt_chem)
    
     DO Jj=1,nspecies
     	DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
     	     spack_2d(ijk,inob)%DLmat(Jj,Jj) = Igamma_dtstep - spack(inob)%DLdrdc(ijk,Jj,Jj)
   	ENDDO
     ENDDO
     
     

!---------------------------------------------------------------------------------------------------------------------
!    1- First step
!-   Compute DLk1 by Solving (1-Igamma*dt*DLRDC) DLk1=DLR 

     !-for solver  LU =============================
     !DO ijk=1,index_g(ngrid)%block_end(i)
     ! CALL solve(spack_2d(ijk,inob)%DLmat,spack_2d(ijk,inob)%DLb1,nspecies,keep=.true.)
     !ENDDO

     !- solver sparse
      DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
        CALL Solve_linear(nspecies,spack_2d(ijk,inob)%DLmat& !matrix  I - gama dt Jac
	                  	  ,spack_2d(ijk,inob)%DLb1 & !matrix independent (P-L)
			   	  ,spack_2d(ijk,inob)%DLk1 ) !matrix solution
			   	  !,i)  	        
      ENDDO
!
!---------------------------------------------------------------------------------------------------------------------
!    2- Second step
!    Compute first-order approximation (sc_p_new)

     dt_chem_i = 1.0d0/dt_chem
     DO Ji=1,nspecies
         DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)

    	    spack(inob)%sc_p_new(ijk,Ji) = spack(inob)%sc_p(ijk,Ji) + Igamma_i * spack_2d(ijk,inob)%DLk1(Ji) 
    	     	  
    	   ! IF (spack(inob)%sc_p_new(ijk,Ji) .LT. Threshold1) THEN
	   !    spack(inob)%sc_p_new(ijk,Ji) = Threshold1
	   !    spack_2d(ijk,inob)%DLk1(Ji) = (spack(inob)%sc_p_new(ijk,Ji) - spack(inob)%sc_p(ijk,Ji)) * dt_chem_i 
     	   ! ENDIF
       ENDDO
     ENDDO


!-   Compute chemical net production terms (DLr= P-L) at final time with the first-order approximation
     CALL fexchem (spack(inob)%sc_p_new  &
     		  ,spack(inob)%rk        &
	  	  ,spack(inob)%DLr       &    
     		  ,nspecies,1,block_end(i) & !index_g(ngrid)%block_end(i)
                  ,maxblock_size,nr)
  
!-    Compute DLB2    
      DO Ji=1,nspecies       
        DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
           spack_2d(ijk,inob)%DLb2(Ji) = spack(inob)%DLr(ijk,Ji) - &
	                                 2.0D0*Igamma_i* dt_chem_i* spack_2d(ijk,inob)%DLk1(Ji)  
       ENDDO
      ENDDO
	 
!-   Compute DLk2 by solving (1-Igamma*dt*DLRDC) DLk2 = DLr- 2 DLK1      
     
     !-for solver  LU	  =============================
     !DO ijk=1,index_g(ngrid)%block_end(i)
     !  CALL solve(spack_2d(ijk,inob)%DLmat,spack_2d(ijk,inob)%DLb2,nspecies,keep=.true.)
     !ENDDO
     !- for solver sparse ============================= 
      DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
 
 	CALL Solve_linear(nspecies,spack_2d(ijk,inob)%DLmat&
	                  	  ,spack_2d(ijk,inob)%DLb2 &
			   	  ,spack_2d(ijk,inob)%DLk2 )
			   	  !,i) 
      ENDDO			  
!-----------------------------------------------------------------------
!    3- Compute concentrations at final time (optimized scheme):
!    sc_p(t+dt) = sc_p (t)  + (3*DLK1 +DLK2) * 1/( 2 gama)
     DO Ji=1,nspecies
     
       DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
       
        !- the solution 
           spack(inob)%sc_p_new(ijk,Ji)= spack(inob)%sc_p(ijk,Ji) + Igamma_i * &
        			       ( 1.5D0 * spack_2d(ijk,inob)%DLk1(Ji) + &
        			         0.5D0 * spack_2d(ijk,inob)%DLk2(Ji) )  
          
        !- keep non-negative values for concentrations
	!   spack(inob)%sc_p_new(ijk,Ji) = max(spack(inob)%sc_p_new(ijk,Ji),Threshold1)
        !
        !- update DLK2, in case sc_p_new is changed by the statement above 	  
	!   spack_2d(ijk,inob)%DLk2(Ji) = 2.0D0*( (spack(inob)%sc_p_new(ijk,Ji)-spack(inob)%sc_p(ijk,Ji)) & 
	!                                         * dt_chem_i - 1.5D0 * spack_2d(ijk,inob)%DLk1(Ji) )         
	 
       ENDDO
         
     ENDDO

!-  Compute the error estimation : spack(inob)%err  
     DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)

       spack(inob)%err(ijk) = 0.0d0
       DO Ji=1,nspecies
         
	 Tol = ATol(Ji) + RTol(Ji)*DMAX1( DABS(spack(inob)%sc_p(ijk,Ji)),DABS(spack(inob)%sc_p_new(ijk,Ji)) )
	 
	 err1= 0.5D0*Igamma_i * (  spack_2d(ijk,inob)%DLk1(Ji)    + &
        			   spack_2d(ijk,inob)%DLk2(Ji)    )  
         spack(inob)%err(ijk) = spack(inob)%err(ijk) + (err1/Tol)**2.0D0
	
       ENDDO

       spack(inob)%err(ijk) = DMAX1( uround, DSQRT( spack(inob)%err(ijk)/nspecies ) )
     ENDDO

     all_accepted = 1
     DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
        if( spack(inob)%err(ijk) - Roundoff > 1.0D0) then
	   all_accepted = 0 ; exit
	endif
     ENDDO

     !- find the maximum error occurred
     max_err1 = maxval( spack(inob)%err(1:block_end(i))) !index_g(ngrid)%block_end(i) ) )

     !- use it to determine the new time step for all block elements
     !- new step size is bounded by FacMin <= Hnew/H <= FacMax
     Fac  = MIN(FacMax,MAX(FacMin,FacSafe/max_err1**(1.0D0/elo)))
	
     !- possible new timestep
     dt_new = dt_chem * Fac
     dt_new = max(dt_min,min(dt_max,dt_new))
     
     !- to reset the timestep resizing in function of the error estimation, use the statements below:
     ! all_accepted = 1; dt_new=dt_max 

     if( all_accepted == 0 .and. dt_new > dt_min) then  ! current solution is not accepted
 
       !- resize the timestep and try again
       dt_chem = dt_new
          
     else    ! current solution is     accepted

      !- go ahead, updating spack(inob)%sc_p with the solution (spack(inob)%sc_p_new)
      !- next time
      time_c = time_c + dt_chem
      
      !- next timestep (dt_new but limited by the time_f-time_c, the rest of time integration interval)
      dt_chem = min(dt_new, time_f-time_c)

      !- save the accepted timestep for the next integration interval
      if(time_c < time_f) last_accepted_dt(i) = dt_new ! index_g(ngrid)%last_accepted_dt(i) =    dt_new      
      
      !- pointer (does not work yet)
      ! spack(inob)%sc_p=>spack(inob)%sc_p_new     ! POINTER
      !- copy
        DO Ji=1,nspecies
           DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
               spack(inob)%sc_p(ijk,Ji) = spack(inob)%sc_p_new(ijk,Ji)
           ENDDO
        ENDDO

	EXIT UntilAccepted
	   
     endif

   END DO UntilAccepted
   

  ENDDO run_until_integr_ends ! time-spliting 


  !-----------------------------------------------------------------------
  !- Restoring species tendencies from internal to brams strucuture

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

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!DO n=1,nspecies
!print*,'max1=',spc_name(n),maxval(chem1_g(n)%sc_p(:,:,:)),maxloc(chem1_g(n)%sc_p(:,:,:))
!print*,'max2=',n,spc_name(n),maxval(spack(:)%sc_p(:,n)),maxloc(spack(:)%sc_p(:,n))
!call flush(6)
!END DO 	
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX			     

 END SUBROUTINE chem_ros_dyndt
  
!---------------------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------------------

END MODULE mod_chem_spack_ros_dyndt
