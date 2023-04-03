!> @brief This module adapts the data from/to CCATT-BRAMS to MAtrix
!! it uses blocks of athmosferic points to prepare the data. 
!! @author Luiz Flavio
!! @date Oct/2011
MODULE DriverMatrix 

#ifdef MATRIX

      USE memMatrix, only: matrixVar,      &
                           aerfrq,&
                           nColumn,&
                           NumberOfColunms,&
                           nmass_spcs,&
                           iPos,&
                           jPos,&
                           nAeroBox,&
                           nGases,&
                           nemis_spcs,&
                           ndiag_aero         , &
                           allocateMatrix     , &
                           allocateMatrixSetup, &
                           deAllocateMatrix   , &
	                   mw_so4,mw_h2so4,mw_no3,mw_hno3,mw_nh4,mw_nh3, &
			   matrixVarRK,MINCONC
			   
                           
      USE mem_aer1   , ONLY: aer1_g,       &
                             aer1_inorg_g, &
                             aer2_g,       &
                             aer2mp_g,     &
			     aer_timestep
                             
      USE mem_chem1  , ONLY: CHEMISTRY,chem1_g,chem1_src_g

      USE chem1_list , ONLY: NH3,SULF,HNO3
      USE aer1_list  , ONLY: mode_alloc, &
                             nmodes,     &
                             ninorg,     &
                             nspecies,   &
                             aer12matrix, &
                             aer1_inorg2matrix, &
                             aer22matrix
                             
         
      PUBLIC ::  MatrixDriver
      
      LOGICAL :: Test_firstTime=.true.
      LOGICAL, PARAMETER :: test=.false.
      INTEGER, PARAMETER :: test_type=1 !1 - normal , 2= fixed
      LOGICAL :: first_call=.true.
      
      LOGICAL :: testcase=.false.

CONTAINS

   !> @brief To prepare and adapt data from brams to Matrix
   !! @author Luiz Flavio
   !! @date Oct/2011
   !! @todo verify sc_t of aerosols
   !! @todo Verify the wupdraft because I use w wind
   !! @todo aqso4rate is fixed using the box value, must be changed
   !! @todo emis_map is fixed to zero. Must be changed to real values
   !!
   subroutine MatrixDriver(ia,iz,ja,jz,m1,m2,m3)
      
      USE mem_basic  , ONLY: basic_g
      use mem_grid   , ONLY: dtlt,ngrid,time,nzg,npatch
      USE rconstants , ONLY: cp, cpor, p00, pi180, cpi ! INTENT(IN)
      USE ModParticle, ONLY: copyDensRadius
      USE mem_radiate, ONLY: radfrq
      USE mem_leaf   , ONLY: leaf_g
      USE mem_turb   , ONLY: turb_g
      USE mem_micro  , ONLY: micro_g
      USE micphys    , ONLY: level,mcphys_type
      USE mem_grid   , ONLY: grid_g, &
                             zt
      USE Aero_setup, ONLY: &
                            Setup_Config, & !Create and setup species in Matrix
                            Setup_Species_Maps, &
                            Setup_Aero_Mass_Map, & 
                            Setup_Coag_Tensors, &
                            Setup_Dp0, &
                            Setup_Emis !, &  
                            !Setup_Kci   
                            
      USE Aero_coag,   ONLY: &
                            setup_kij_diameters,       & !subroutine
                            setup_kij_tables,          & !subroutine
                            get_kbarnij                  !subroutine
      
      USE Aero_subs,   ONLY: &
                            no_microphysics            !intent()                            

      USE Aero_npf,    ONLY: &
                            setup_npfmass              !subroutine
      use node_mod, only: mynum,  &
                          mchnum, &
                          master_num

      USE Aero_Isrpia, only : InitIsrpia

      implicit none

      !ngases     = 3      ! number of gas-phase species
      !nmass_spcs = 5      ! total number of mass species
      !gas_h2so4  = 1      !\
      !gas_hno3   = 2      !-indices in the gas array
      !gas_nh3    = 3      !/
 
      INTEGER, INTENT(IN) :: ia !< i initial point
      INTEGER, INTENT(IN) :: iz !< i final point
      INTEGER, INTENT(IN) :: ja !< j initial point
      INTEGER, INTENT(IN) :: jz !< j final point
      INTEGER, INTENT(IN) :: m1 !< z size               
      INTEGER, INTENT(IN) :: m2 !< i size               
      INTEGER, INTENT(IN) :: m3 !< j size               
       
      CHARACTER(LEN=20):: aer_split_method ="PARALLEL"
     !CHARACTER(LEN=20):: aer_split_method ="SYMMETRIC"
      INTEGER, PARAMETER :: int_method=1 ! 1,2 or 3
      LOGICAL :: do_Setup_Kci_dynj=.true.
      INTEGER :: i,j,k,noc,point,na,ng,gasOfChem,md,sp,nm,nsp,nmd,npat,ning
      double precision :: picpi,tstep,pblht,dn0i
      DOUBLE PRECision,DIMENSION(m1,ia:iz,ja:jz,nmodes) :: ac  
      DOUBLE PRECision,DIMENSION(m1,ia:iz,ja:jz,nmodes) :: dens_mode_dry 
      REAL,PARAMETER :: tkethrsh=0.001       !   tke threshold for PBL height in m2/s2
                                             !   tkmin    = 5.e-4   minimum TKE in RAMS 
      REAL,PARAMETER :: rcmin=1.e-4          !   min liq water = 0.1 g/kg
      INTEGER :: toci,tocj,fcoi,fcoj
      REAL,    EXTERNAL :: rs        ! Function

      INTEGER, PARAMETER :: isplit =  1 !- use 1 for sequential splitting
                                        !- use 2 for symmetric  splitting
      INTEGER :: N_DYN_AER,kij,ij
      REAL :: DTI

      IF(first_call) THEN
         call InitIsrpia
         CALL allocateMatrixSetup()
         CALL Setup_Config()
         call Setup_Species_Maps
         call Setup_Aero_Mass_Map 
         if ( .not. no_microphysics ) then
            call Setup_Coag_Tensors
            call Setup_Dp0     
            call Setup_Kij_Diameters
            call Setup_Kij_Tables    
            call Setup_Emis  
            !call Setup_Kci   
            call Setup_Npfmass
         END IF
         first_call=.false.
      END IF

     
      aerfrq=aer_timestep
      !--adding control to call matrix aerosol only at N_DYN_AER*dtlong
      !--with the aerosol updating taking place in the middle (0.5*N_DYN_AER*dtlong)
      !--for nested grids, this setting needs further checking.
      N_DYN_AER =max(1,nint(aer_timestep/dtlt)) 
      
      IF (MOD(time-(N_DYN_AER/ISPLIT-1)*dtlt , N_DYN_AER*dtlt) .ne. 0. ) RETURN
!     IF (MOD(time                           , N_DYN_AER*dtlt) .ne. 0. ) RETURN
      !print*,"calling OLD matrix aer model time=",time,dtlt;call flush(6)
   
   
      !IF (.not. (mod(time + .001,aerfrq) .lt. dtlt .or. time .lt. 0.001)) RETURN
     
     
      CALL allocateMatrix(ia,iz,ja,jz,m1,m2,m3,nzg) !Allocating the matrix with right size
      CALL MakeCrossMap(ia,iz,ja,jz,m1)             !Doing point to point cross reference
      
      !Pressure, air temperature and Moisture
      DO i=ia,iz
         do j=ja,jz
            noc=nColumn(i,j)
	    matrixVar(noc)%ustar(1)=0.
            matrixVar(noc)%tstar(1)=0.
            do npat=1,npatch
	       matrixVar(noc)%ustar(1)= matrixVar(noc)%ustar(1) + &
	                     leaf_g(ngrid)%ustar(i,j,npat) *leaf_g(ngrid)%patch_area(i,j,npat)  
            
	       matrixVar(noc)%tstar(1)= matrixVar(noc)%tstar(1) + &
	                     leaf_g(ngrid)%tstar(i,j,npat) *leaf_g(ngrid)%patch_area(i,j,npat) 
            enddo
	    !print*,"surf=",i,j,matrixVar(noc)%ustar(1),matrixVar(noc)%tstar(1)
	    
	    DO k = 1,m1
               !calculating and moving to right position
               picpi = (basic_g(ngrid)%pi0(k,i,j) + basic_g(ngrid)%pp(k,i,j)) * cpi
               !-Pressure (hPa)
               matrixVar(noc)%pres(k)= dble(p00 * picpi ** cpor)
               !-Temperature (K)
               matrixVar(noc)%tk(k) = dble(basic_g(ngrid)%theta(k,i,j) * picpi)
	       
	       !Humidity (0-1)
               matrixVar(noc)%rh(k)=dble(basic_g(ngrid)%rv(k,i,j)&
	                        /rs( real(matrixVar(noc)%pres(k),4),real(matrixVar(noc)%tk(k),4)))
	       matrixVar(noc)%rh(k)=min(1.0d0,max(0.05D0,matrixVar(noc)%rh(k)))
               
	       !Vertical Wind Component (m/s)
               matrixVar(noc)%wupdraft(k)=dble( max(0.01,basic_g(ngrid)%wp(k,i,j)) )
               
	       matrixVar(noc)%dn0(k)=dble( basic_g(ngrid)%dn0(k,i,j) ) 
               matrixVar(noc)%pi0(k)=dble( basic_g(ngrid)%pi0(k,i,j) )
               matrixVar(noc)%pp (k)=dble( basic_g(ngrid)%pp (k,i,j) )
	       !print*,"Column K=",k,matrixVar(noc)%pres(k),matrixVar(noc)%tk(k),matrixVar(noc)%rh(k)
            END DO

            matrixVar(noc)%Zi = zt(2)*grid_g(ngrid)%rtgt(i,j) 
            !- convective layer
            if(turb_g(ngrid)%sflux_t(i,j) >= 1.e-8) then 
              pblht=zt(2)*grid_g(ngrid)%rtgt(i,j) 
              do k=2,m1-1
                 pblht=zt(k)*grid_g(ngrid)%rtgt(i,j) 
                 if( micro_g(ngrid)%rcp(k,i,j) .gt. rcmin     ) EXIT ! dry convective layer
                 if( turb_g(ngrid)%tkep(k,i,j) .le. tkethrsh  ) EXIT 
              enddo
              matrixVar(noc)%Zi=pblht
            endif
         END DO
      END DO
!===== tmp
      !print*,"matrix 1" ;  call flush(6)
!      print*,"stars:",maxval(matrixVar(1:numberOfColunms)%tstar),maxval(matrixVar(:)%ustar)
!      print*,"press:",maxval(matrixVar(:)%pres(:)),minval(matrixVar(:)%pres(:))
!      print*,"tk:",maxval(matrixVar(:)%tk(:)),minval(matrixVar(:)%tk(:))
!      print*,"rh:",maxval(matrixVar(:)%rh(:)),minval(matrixVar(:)%rh(:))
!      call flush(6)
!====

      !Must be changed
      DO noc=1,NumberOfColunms
         matrixVar(noc)%aqso4rate(:)=0.5D-4 !Fixed
         matrixVar(noc)%emis_mas(:,:)=0.0 !No emission of any aerossols
         matrixVar(noc)%aero=0.0
      END DO
      !print*,"matrix 2" ;  call flush(6)
      
      !Gas input
      DO i=ia,iz
         do j=ja,jz
            noc=nColumn(i,j)
            DO k=1,m1
	       !- sulfuric acid express in terms of SO4 concentration
               matrixVar(noc)%gas(k,1)=0.975*matrixVar(noc)%dn0(k)*chem1_g(SULF,ngrid)%sc_p(k,i,j)&
	                               *mw_so4/mw_h2so4
	                               
	       matrixVar(noc)%gas(k,2)=matrixVar(noc)%dn0(k)*chem1_g(HNO3,ngrid)%sc_p(k,i,j)&
	                               *mw_no3/mw_hno3
				       	                                
               matrixVar(noc)%gas(k,3)=matrixVar(noc)%dn0(k)*chem1_g(NH3 ,ngrid)%sc_p(k,i,j)&
	                               *mw_nh4/mw_nh3
	               
            END DO
         END DO
      END DO
      !CALL fillAero(ia,iz,ja,jz,m1) !For the first time copy aerosol data to Matrix

      !Fill the matrixVar type with the data with aerosol ammount from model
      DO i=ia,iz
         do j=ja,jz
            noc=nColumn(i,j)
            DO k = 1,m1
               !Fill inorganic
               DO ning=1,ninorg
                  matrixVar(noc)%aero(k,aer1_inorg2matrix(ning))=matrixVar(noc)%dn0(k)*aer1_inorg_g(ning,ngrid)%sc_p(k,i,j)!+100.
               END DO
               !Fill number
               DO nmd=1,nmodes
                  IF(aer22matrix(nmd)==0) CYCLE
                  matrixVar(noc)%aero(k,aer22matrix(nmd))=matrixVar(noc)%dn0(k)*aer2_g(nmd,ngrid)%sc_p(k,i,j)!+100.
               END DO
               !Fill mass
               DO nmd=1,nmodes
                  DO nsp=1,nspecies
                     IF(aer12matrix(nmd,nsp)==0) CYCLE
                     matrixVar(noc)%aero(k,aer12matrix(nmd,nsp))=matrixVar(noc)%dn0(k)*aer1_g(nmd,nsp,ngrid)%sc_p(k,i,j)!+100.
                  END DO
               END DO
            END DO
         END DO
      END DO
      
      !- setting the aer timestep. 
      tstep=aerfrq               

      !Calling the matrix adapter (will be changed in future to blocking and vectorizing)


      IF(testCase) THEN
         CALL Matrix_test(matrixVar(1),m1)
      ELSE
        DO noc=1,NumberOfColunms
          if	(int_method == 1) then
             do_Setup_Kci_dynj=.true.
	     CALL Matrix(matrixVar(noc),tstep,m1,noc,1,.false.,do_Setup_Kci_dynj)
          
          elseif(int_method == 2) then ! RK2
             
	     !-1st step
	     do_Setup_Kci_dynj=.true.
             matrixVarRK(1)%aero(:,:)= matrixVar(noc)%aero(:,:) ! Y0
             matrixVarRK(1)%gas (:,:)= matrixVar(noc)%gas (:,:) 
             
	     CALL Matrix(matrixVar(noc),0.5D+00*tstep,m1,noc,1,.false.,do_Setup_Kci_dynj)
             
	     matrixVarRK(2)%aero(:,:)= matrixVar(noc)%aero(:,:) ! Y*=Y0+0.5dt F(Y0)
             matrixVarRK(2)%gas (:,:)= matrixVar(noc)%gas (:,:) 

             !-2nd step and final solution
	     do_Setup_Kci_dynj=.false.
             
	     ! - in: matrixVar(noc)=Y* , out: matrixVar(noc)=Y**=Y* + dt F(Y*)
	     CALL Matrix(matrixVar(noc),tstep,m1,noc,1,.false.,do_Setup_Kci_dynj)
             
	     ! Y(t+dt) = Y0 + dt F(Y*), F(Y*)=(Y**-Y*)/dt
	     matrixVar(noc)%aero(:,:)= max(matrixVarRK(1) %aero(:,:)+ &
	                                   matrixVar(noc) %aero(:,:)- &
					   matrixVarRK(2) %aero(:,:), MINCONC)
             
	     matrixVar(noc)%gas (:,:)= max(matrixVarRK(1) %gas(:,:)+ &
	                                   matrixVar(noc) %gas(:,:)- &
					   matrixVarRK(2) %gas(:,:), MINCONC)
            
          
          elseif(int_method == 3) then ! RK3
            !srf needs to be completed
	    stop "int_method must be 1 or 2"
          endif
        END DO
      END IF
    
     
      !LFR  IF(test) THEN
      !LFR  !CALL compareTest(testData,matrixVar(1),1,m1)
      !LFR  END IF
      IF(aer_split_method /= 'PARALLEL')THEN
       !Fill the aerosol with the data from matrix
       DO i=ia,iz
          do j=ja,jz
            noc=nColumn(i,j)
            DO k = 1,m1
	       dn0i=1.0d0/matrixVar(noc)%dn0(k)
	    
               !Fill inorganic
               DO ning=1,ninorg
                  aer1_inorg_g(ning,ngrid)%sc_p(k,i,j)=matrixVar(noc)%aero(k,aer1_inorg2matrix(ning))*dn0i
               END DO
               !Fill number
               DO nmd=1,nmodes
                  IF(aer22matrix(nmd)==0) CYCLE
                  aer2_g(nmd,ngrid)%sc_p(k,i,j)=matrixVar(noc)%aero(k,aer22matrix(nmd))*dn0i
               END DO
               !Fill mass
               DO nmd=1,nmodes
                  DO nsp=1,nspecies
                     IF(aer12matrix(nmd,nsp)==0) CYCLE
                     aer1_g(nmd,nsp,ngrid)%sc_p(k,i,j)=matrixVar(noc)%aero(k,aer12matrix(nmd,nsp))*dn0i
                  END DO
               END DO
            END DO
         END DO
       END DO
      ELSEIF(aer_split_method == 'PARALLEL'  .and. N_DYN_AER== 1) THEN
       DTI=1./DTLT
       DO i=ia,iz
          DO j=ja,jz
            noc=nColumn(i,j)
	    ij =m1*(i-1) + (m1*(m2-1)+m1)*(j-1)   
            
	    DO k = 1,m1
	       dn0i=1.0d0/matrixVar(noc)%dn0(k)
	       kij= k+ ij   

               !Fill inorganic
               DO ning=1,ninorg     
                  aer1_inorg_g(ning,ngrid)%sc_t(kij)=aer1_inorg_g(ning,ngrid)%sc_t(kij)                 + &
		                                    (matrixVar(noc)%aero(k,aer1_inorg2matrix(ning))*dn0i  &
		                                     -aer1_inorg_g(ning,ngrid)%sc_p(k,i,j)                )*DTI
               END DO
               !Fill number
               DO nmd=1,nmodes
                  IF(aer22matrix(nmd)==0) CYCLE
                  aer2_g(nmd,ngrid)%sc_t(kij)=aer2_g(nmd,ngrid)%sc_t(kij)                  + &
		                             (matrixVar(noc)%aero(k,aer22matrix(nmd))*dn0i   &
		                              -aer2_g(nmd,ngrid)%sc_p(k,i,j)                 )*DTI
               END DO
               !Fill mass
               DO nmd=1,nmodes
                  DO nsp=1,nspecies
                     IF(aer12matrix(nmd,nsp)==0) CYCLE
                     aer1_g(nmd,nsp,ngrid)%sc_t(kij)=aer1_g(nmd,nsp,ngrid)%sc_t(kij)                + &
		                                    (matrixVar(noc)%aero(k,aer12matrix(nmd,nsp))*dn0i &
		                                     -aer1_g(nmd,nsp,ngrid)%sc_p(k,i,j)               )*DTI
                  END DO
               END DO
      
            END DO
         END DO
       END DO
      ELSE
      
       STOP "UNKNOW INTEGRATION METHOD FOR MATRIX"
       
      ENDIF 
     
      !-kml/srf - for microphysics activation
      IF(mcphys_type == 3) then
         DO i=ia,iz
           DO j=ja,jz
            noc=nColumn(i,j)
            DO k = 1,m1
	       aer2mp_g(1,ngrid)%kappa_eff (k,i,j)=matrixVar(noc)%aer2mp_eff(k,1)
	       aer2mp_g(1,ngrid)%diam_eff  (k,i,j)=matrixVar(noc)%aer2mp_eff(k,2)
	       aer2mp_g(1,ngrid)%numb_water(k,i,j)=matrixVar(noc)%aer2mp_eff(k,3)
	       aer2mp_g(1,ngrid)%numb_ice  (k,i,j)=matrixVar(noc)%aer2mp_eff(k,4)
            END DO
         END DO
       END DO
      ENDIF
      
!-lixo----------------------------------------------------------------------------------      
     go to 4444
     
     !Fill inorganic
     DO ning=1,ninorg
     	print*,"INOR max min=",maxval(aer1_inorg_g(ning,ngrid)%sc_p),minval(aer1_inorg_g(ning,ngrid)%sc_p);call flush(6)
     END DO
     !Fill number
     DO nmd=1,nmodes
     	IF(aer22matrix(nmd)==0) CYCLE
     	print*,"NUMB max min=",maxval(aer2_g(nmd,ngrid)%sc_p),minval(aer2_g(nmd,ngrid)%sc_p);call flush(6)
     END DO
     !Fill mass
     DO nmd=1,nmodes
     	DO nsp=1,nspecies
     	   IF(aer12matrix(nmd,nsp)==0) CYCLE
     	   print*,"MASS max min=",maxval(aer1_g(nmd,nsp,ngrid)%sc_p),minval(aer1_g(nmd,nsp,ngrid)%sc_p);call flush(6)
     	END DO
     END DO
4444 continue
      !print*,"matrix end", time ; call flush(6)
!     stop "matrix orig"
!-lixo----------------------------------------------------------------------------------      
      
      
      IF (mchnum==master_num) print*,"matrix end", time ; call flush(6)
      
      !CALL DumpAero(ia,iz,ja,jz,m1)
      !***************************************************************!
      !TESTES TESTES TESTES TESTES TESTES TESTES TESTES TESTES        !
      !***************************************************************!
                                !
      CALL deAllocateMatrix !Deallocating the matrix before return    !
      RETURN                                                          !
                                                                      !
      !***************************************************************!
      
   
   END subroutine MatrixDriver
   
   !> @brief create a xreference map between ia,iz,ja,jz,m1 and pointOfBlock  \n
   !! @author Luiz Flavio
   !! @date Oct/2011
   subroutine makeCrossMap(ia,iz,ja,jz,m1)
      implicit none
      
      INTEGER, INTENT(IN) :: ia,iz,ja,jz,m1
      INTEGER :: i,j,k,nb,p,column
      column=0
      do i=ia,iz
         do j=ja,jz
               column=column+1
               iPos(column)=i
               jPos(column)=j
               nColumn(i,j)=column
         END DO
      END DO

   END subroutine makeCrossMap
   
   SUBROUTINE Matrix_test(var,m1)
      USE memMatrix,   only: &
                            matrixType
                            
      TYPE(matrixType),INTENT(INOUT) :: var
      INTEGER, INTENT(IN) :: m1
      
      CHARACTER(len=6) :: label
      INTEGER :: i,it,irep,n,ii,iii,k
      REAL(8) :: laero(naerobox)
      REAL(8) :: lgas(ngases)
      REAL(8) :: ltk                    ! absolute temperature [K]
      REAL(8) :: lrh                    ! relative humidity    [0-1]
      REAL(8) :: lpres                  ! ambient pressure     [Pa]  
      REAL(8) :: ltstep                 ! model time step [s]
      REAL(8) :: ltsubstep              ! increment of model time step [s]
      REAL(8) :: lemis_mass(nemis_spcs) ! mass emission rates [ug/m^3/s]
      REAL(8) :: laqso4rate             ! incloud SO4 production rate [ug/m^3/s]
      REAL(8) :: lwupdraft = 0.5d+00    ! cloud updraft velocity [m/s]
      REAL(8) :: ldiag(ndiag_aero,naerobox)        ! budget or tendency diagnostics [ug/m^3/s] or [#/m^3/s]
      REAL(8) :: ldiag_accum(ndiag_aero,naerobox)  ! budget or tendency diagnostics [ug/m^3/s] or [#/m^3/s]
      integer :: icount
      
      WRITE (*,fmt='(A)') 'Testando. Abrindo arquivo ./tables/matrix/testmatrix.in'
      OPEN(88,file='./tables/matrix/testmatrix.in',status='old')
      OPEN(89,file='./tables/matrix/testcompar.out',status='replace')
       WRITE (*,fmt='(A)') '-------------------------------------------------------'
      WRITE (*,fmt='(A3,1X,4(A2,1X))') 'cnt',' i','it','ir','nr'
      icount = 0
      DO
         icount=icount +1
         READ (88,fmt='(A,4(I2.2,1X))',END=100) label,i,it,irep,n
         WRITE (*,fmt='(i3.3,1X,4(I2.2,1X))') icount,i,it,irep,n
         READ (88,fmt='(51(E13.6,1X))',END=100) (laero(ii),ii=1,naerobox)
         READ (88,fmt='(3(E13.6,1X))',END=100) (lgas(ii),ii=1,ngases)
         READ (88,fmt='(10(E13.6,1X))',END=100) (lemis_mass(ii),ii=1,nemis_spcs)
         READ (88,fmt='(6(E13.6,1X))',END=100) ltsubstep,ltk,lrh,lpres,laqso4rate,lwupdraft
         DO iii=1,ndiag_aero
            READ (88,fmt='((E13.6,1X))',END=100) (ldiag(iii,ii),ii=1,naerobox)
         END DO
      
         DO k=1,m1
            Var%aero(k,:)     = laero(:) !(m1,naerobox))
            Var%gas(k,:)      = lgas(:)  !(m1,ngases))
            Var%emis_mas(k,:) = lemis_mass(:) !(m1,nemis_spcs))
            Var%tk(k)         = ltk !(m1))
            Var%rh(k)         = lrh !(m1))
            Var%pres(k)       = lpres !(m1))
            Var%aqso4rate(k)  = laqso4rate!(m1))
            Var%wupdraft(k)   = lwupdraft!(m1))
            Var%diag(k,:,:)     = ldiag(:,:) !(m1,ndiag_aero,naerobox))
         END DO

         CALL Matrix(Var,ltsubstep,m1,1,1,.true.)

         DO k=1,m1
            laero(:)      = Var%aero(1,:)     
            lgas(:)       = Var%gas(1,:)      
            lemis_mass(:) = Var%emis_mas(1,:) 
            ltk           = Var%tk(1)         
            lrh           = Var%rh(1)         
            lpres         = Var%pres(1)       
            laqso4rate    = Var%aqso4rate(1)  
            lwupdraft     = Var%wupdraft(1)   
            ldiag(:,:)    = Var%diag(1,:,:)   
         END DO


         write (89,fmt='(A,4(I2.2,1X))') label,i,it,irep,n
         write (89,fmt='(51(E13.6,1X))') (laero(ii),ii=1,naerobox)
         write (89,fmt='(3(E13.6,1X))') (lgas(ii),ii=1,ngases)
         write (89,fmt='(10(E13.6,1X))') (lemis_mass(ii),ii=1,nemis_spcs)
         DO iii=1,ndiag_aero
            WRITE (89,fmt='((E13.6,1X))') (ldiag(iii,ii),ii=1,naerobox)
         END DO
         
      END DO 
100   WRITE (*,fmt='(A)') '-----------------------------------------------------------------------'
      WRITE (*,fmt='(A)') '!!! Encerrado. Gerado arquivo de saida ./tables/matrix/testcompar.out !'
      WRITE (*,FMT='(A)') 'Para ser comparado com ./tables/matrix/testmatrix.out (original Matrix)'
      WRITE (*,FMT='(A)') 'Use: diff ./tables/matrix/testcompar.out ./tables/matrix/testmatrix.out'
      WRITE (*,fmt='(A)') '-----------------------------------------------------------------------'
      STOP
      
   END SUBROUTINE Matrix_test

#else 

  PUBLIC ::  MatrixDriver


      LOGICAL :: Test_firstTime=.true.
      LOGICAL, PARAMETER :: test=.false.
      INTEGER, PARAMETER :: test_type=1 !1 - normal , 2= fixed
      LOGICAL :: first_call=.true.
      
      LOGICAL :: testcase=.false.


  CONTAINS

   !> @brief To prepare and adapt data from brams to Matrix
   !! @author Luiz Flavio
   !! @date Oct/2011
   !! @todo verify sc_t of aerosols
   !! @todo Verify the wupdraft because I use w wind
   !! @todo aqso4rate is fixed using the box value, must be changed
   !! @todo emis_map is fixed to zero. Must be changed to real values
   !!
   subroutine MatrixDriver(ia,iz,ja,jz,m1,m2,m3)

      write (*,fmt='(A)') 'Matrix only works with AER=MATRIX and CHEM=RELACS_MX'
      write (*,fmt='(A)') 'Please run config and compile the code from begin'
      STOP

   end subroutine MatrixDriver

#endif
END MODULE DriverMatrix  

