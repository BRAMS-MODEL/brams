MODULE mod_chem_spack_qssa

  USE mem_chem1, ONLY: &
       chem1_vars ! Type

  USE mod_chem_spack_kinetic, ONLY: &
       kinetic    ! Subroutine
  
  USE mod_chem_spack_fexloss, ONLY: &
       fexloss    ! Subroutine

  USE mod_chem_spack_fexprod, ONLY: &
       fexprod    ! Subroutine

  USE mod_chem_spack_dratedc, ONLY: &
       dratedc    ! Subroutine

  USE mod_chem_spack_rates, ONLY: &
       rates      ! Subroutine


  IMPLICIT NONE

  PRIVATE


  PUBLIC :: chem_qssa, & ! Subroutine
            solveur,   & ! Subroutine
            cvmgm_dp,  & ! Fuction
            cvmgz_dp,  & ! Function
            cvmgp_dp     ! Function

CONTAINS


  !========================================================================================
  SUBROUTINE chem_qssa(m1,m2,m3,itchim,nrmax,mynum,dtlt,pp,pi0,theta,rv, &
                       cosz,jphoto,nspecies,nr_photo,weight,PhotoJMethod,chem1_g, &
                       nob,maxblock_size, &
                       block_end,indexk,indexi,indexj,kij_index, &
                       cp,cpor,p00,att)
  !========================================================================================

    ! mem_grid
    REAL             , INTENT(IN)    :: dtlt

    ! node_mod
    INTEGER          , INTENT(IN)    :: mynum
    INTEGER          , INTENT(IN)    :: m1
    INTEGER          , INTENT(IN)    :: m2
    INTEGER          , INTENT(IN)    :: m3

    ! mem_basic
    REAL             , INTENT(IN)    :: pp(m1,m2,m3)
    REAL             , INTENT(IN)    :: pi0(m1,m2,m3)
    REAL             , INTENT(IN)    :: theta(m1,m2,m3)
    REAL             , INTENT(IN)    :: rv(m1,m2,m3)

    ! chemistry_driver() local variable
    INTEGER          , INTENT(IN)    :: itchim

    ! mem_radiate
    REAL             , INTENT(IN)    :: cosz(m2,m3)

    ! chem1_list
    INTEGER          , INTENT(IN)    :: nspecies
    INTEGER          , INTENT(IN)    :: nr_photo
    INTEGER          , INTENT(IN)    :: nrmax
    REAL             , INTENT(IN)    :: weight(nspecies)
    CHARACTER(LEN=10), INTENT(IN)    :: PhotojMethod

    ! FastJX
    ! Fast_JX
    DOUBLE PRECISION, INTENT(IN) :: jphoto(nr_photo,m1,m2,m3)

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

    ! rconstants
    REAL, INTENT(IN) :: cp
    REAL, INTENT(IN) :: cpor
    REAL, INTENT(IN) :: p00

    ! uv_atten
    REAL, INTENT(IN) :: att(m2,m3)


    !  INTEGER,PARAMETER :: block_size=45
    !  INTEGER,PARAMETER :: block_size=4
    !  INTEGER,PARAMETER :: block_size=1
    DOUBLE PRECISION,PARAMETER     :: pmar=28.96  
    ! !LPCE
    !srf
    INTEGER :: i,ijk,n
    DOUBLE PRECISION dtlt_i,fxc
    INTEGER k_,i_,j_,kij_


    TYPE spack_type_
       !3d Real
       DOUBLE PRECISION,POINTER,DIMENSION(:,:,:) :: dw
       !2d Real
       DOUBLE PRECISION,POINTER,DIMENSION(:,:)   :: jphoto 
       DOUBLE PRECISION,POINTER,DIMENSION(:,:)   :: loss	
       DOUBLE PRECISION,POINTER,DIMENSION(:,:)   :: prod	
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
       !DOUBLE PRECISION,POINTER,DIMENSION(:)     :: xlw 
    END TYPE spack_type_

    TYPE(spack_type_) :: spack

!!$  !Determining number of blocks and masks
!!$  CALL AllocIndex(block_size,mynum,mmxp,mmyp,mmzp,mia,miz,mja,mjz,mibcon,mi0,mj0,&
!!$       ngrids,dtlongn)

    dtlt_i=1./dtlt


    !- Allocating Spaces to copy structure
    !ALLOCATE(spack(nob))
    !- 3D Variables
    ALLOCATE(spack%dw    (1:maxblock_size,nrmax,nspecies))
    !- 2D Variables
    ALLOCATE(spack%jphoto(1:maxblock_size,nr_photo))
    ALLOCATE(spack%sc_p  (1:maxblock_size,nspecies))
    ALLOCATE(spack%loss  (1:maxblock_size,nspecies))
    ALLOCATE(spack%prod  (1:maxblock_size,nspecies))
    ALLOCATE(spack%rk    (1:maxblock_size,nrmax))
    ALLOCATE(spack%w     (1:maxblock_size,nrmax))
    !- 1D variables
    ALLOCATE(spack%temp  (1:maxblock_size))
    ALLOCATE(spack%press (1:maxblock_size))
    ALLOCATE(spack%cosz  (1:maxblock_size))
    ALLOCATE(spack%att   (1:maxblock_size))
    ALLOCATE(spack%vapp  (1:maxblock_size))
    ALLOCATE(spack%volmol(1:maxblock_size))
    !ALLOCATE(spack%xlw   (1:maxblock_size))


    DO i=1,nob

       !- Copying structure from input to internal
       DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)

          k_=indexk(ijk,i) !index_g(ngrid)%indexk(ijk,i) ! k brams index
          i_=indexi(ijk,i) !index_g(ngrid)%indexi(ijk,i) ! i brams index
          j_=indexj(ijk,i) !index_g(ngrid)%indexj(ijk,i) ! j brams index

          !- Exner function divided by Cp
          spack%press(ijk)=(  pp (k_,i_,j_) &
               + pi0(k_,i_,j_) )/cp
          !- Air temperature (K)			      			      
          spack%temp(ijk)= theta(k_,i_,j_)*spack%press(ijk)

          !- transform from Exner function to pressure
          spack%press(ijk)= spack%press(ijk)**cpor*p00

          !- Water vapor 
          spack%vapp(ijk)=  rv(k_,i_,j_)!&
          !* basic_g(ngrid)%dn0(k_,i_,j_)*6.02e20/18.


          !spack%volmol(ijk)=(spack%press(ijk)*1.e-3)/(8.314e0*spack%temp(ijk))
          spack%volmol(ijk)=(spack%press(ijk))/(8.314e0*spack%temp(ijk))
          !LPCE
          DO n=1,nspecies

             !LPCE conversion from ppbm to molecule/cm^3
             !spack%sc_p(ijk,n)=chem1_g(n)%sc_p(k_,i_,j_) &
             !                    *(6.02e23*1e-12*spack%volmol(ijk))
             spack%sc_p(ijk,n)=chem1_g(n)%sc_p(k_,i_,j_) &
                  *(6.02e23*1e-15*spack%volmol(ijk))*pmar/weight(n)
          END DO
          !spack%xlw(ijk) =  micro_g(ngrid)%rcp(k_,i_,j_) &
          !                   * basic_g(ngrid)%dn0(k_,i_,j_)*1e-3
          spack%cosz(ijk)=cosz(i_,j_)
       END DO

       !- FAST-JX section
       IF(TRIM(PhotojMethod) == 'FAST-JX') THEN
          DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)

             k_=indexk(ijk,i) !index_g(ngrid)%indexk(ijk,i) ! k brams index
             i_=indexi(ijk,i) !index_g(ngrid)%indexi(ijk,i) ! i brams index
             j_=indexj(ijk,i) !index_g(ngrid)%indexj(ijk,i) ! j brams index

             DO n=1,nr_photo
                spack%jphoto(ijk,n)= jphoto(n,k_,i_,j_) !fjphoto%jphoto(n,k_,i_,j_)
             END DO

          ENDDO
     
       ELSEIF(trim(PhotojMethod) == 'LUT') then
      
          DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
        
            i_=indexi(ijk,i) !index_g(ngrid)%indexi(ijk,i) ! i brams index
            j_=indexj(ijk,i) !index_g(ngrid)%indexj(ijk,i) ! j brams index
        
	   !- UV  attenuation un function AOT
           spack%att(ijk)= att(i_,j_)
        
	   !- get zenital angle (for LUT PhotojMethod)
            spack%cosz(ijk)=  cosz(i_,j_)

          ENDDO     
       ENDIF

       !- call kinetic
       ! IN: nr_photo,spack%Jphoto,spack%temp,spack%vapp,spack%Press,spack%cosz
       !     1.0D0,1,index_g(ngrid)%block_end(i),maxblock_size,nrmax
       ! OUT: spack%rk
       CALL kinetic(nr_photo,spack%Jphoto   &
                           ,spack%rk    &
                           ,spack%temp  &
			   ,spack%vapp   &
                           ,spack%Press &
		           ,spack%cosz  &
    			   ,spack%att     &
			   ,1,block_end(i) & !index_g(ngrid)%block_end(i)
                           ,maxblock_size,nrmax)
!!$       CALL kinetic(nr_photo,spack%Jphoto,spack%rk,spack%temp,spack%vapp, &
!!$                    spack%Press,1,index_g(ngrid)%block_end(i),maxblock_size, &
!!$                    nrmax)

       !- call prod, loss and rates
       ! IN: spack%rk,spack%sc_p,nspecies,1,index_g(ngrid)%block_end(i),maxblock_size,nrmax
       ! OUT: spack%w
       CALL rates(spack%rk,spack%sc_p,spack%w,nspecies,1,block_end(i), & !index_g(ngrid)%block_end(i), &
                  maxblock_size,nrmax)

       ! IN: spack%rk,spack%sc_p,nspecies,1,index_g(ngrid)%block_end(i),maxblock_size,nrmax
       ! OUT: spack%dw
       CALL dratedc(spack%rk,spack%sc_p,spack%dw,nspecies,1, &
                    block_end(i), & !index_g(ngrid)%block_end(i),
                    maxblock_size,nrmax)

       ! IN: spack%w,nspecies,1,index_g(ngrid)%block_end(i),maxblock_size,nrmax
       ! OUT: spack%prod
       CALL fexprod(spack%w,spack%prod,nspecies,1,block_end(i), & !index_g(ngrid)%block_end(i), &
                    maxblock_size,nrmax)

       ! IN: spack%dw,,nspecies,1,index_g(ngrid)%block_end(i),maxblock_size,nrmax
       ! OUT: spack%loss
       CALL fexloss(spack%dw,spack%loss,nspecies,1,block_end(i), & !index_g(ngrid)%block_end(i), &
                    maxblock_size,nrmax)

       !- solving Jacobian Matrix
       ! IN: spack%prod,spack%loss,nspecies,dtlt,maxblock_size,1,index_g(ngrid)%block_end(i),itchim,mynum 
       ! INOUT: spack%sc_p
       CALL solveur(spack%prod,spack%loss,spack%sc_p,nspecies,dtlt,maxblock_size,1, &
                    block_end(i), & !index_g(ngrid)%block_end(i),
                    itchim,mynum)

       !- restoring species loss/prod from internal to ccatt-brams structure
       DO ijk=1,block_end(i) !index_g(ngrid)%block_end(i)
          !- tendency index (mude para calcular somente uma vez, isto e', crie um kij(ijk,i)  )
          k_   =indexk(ijk,i) !index_g(ngrid)%indexk(ijk,i)   !- k   brams index
          i_   =indexi(ijk,i) !index_g(ngrid)%indexi(ijk,i)   !- i   brams index
          j_   =indexj(ijk,i) !index_g(ngrid)%indexj(ijk,i)   !- j   brams index
          kij_   =kij_index(ijk,i) !index_g(ngrid)%kij_index(ijk,i)!- kij brams tendency index

          DO n=1,nspecies
             !LPCE conversion from molecule/cm^3 to ppbm
             !chem1_g(n)%sc_p(k_,i_,j_)= spack%sc_p(ijk,n)&
             !                  /(6.02e23*1e-12*spack%volmol(ijk))
             fxc = 1./(6.02e23*1e-15*spack%volmol(ijk)*(pmar/weight(n)))
             chem1_g(n)%sc_t(kij_) = chem1_g(n)%sc_t(kij_)   +  &    ! previous tendency
                  (  spack%sc_p(ijk,n) * fxc  -  &    ! new mixing ratio
                  chem1_g(n)%sc_p(k_,i_,j_))&! old mixing ratio
                  *  dtlt_i                              ! inverse of timestep

          END DO
       END DO


    END DO

    !Deallocating local structure for later use
    !3D Variables
    DEALLOCATE(spack%dw)
    !2D Variables
    DEALLOCATE(spack%jphoto)
    DEALLOCATE(spack%sc_p)
    DEALLOCATE(spack%loss)
    DEALLOCATE(spack%prod)
    DEALLOCATE(spack%rk)
    DEALLOCATE(spack%w)
    DEALLOCATE(spack%cosz)
    DEALLOCATE(spack%att)
    !1D variables
    DEALLOCATE(spack%temp)
    DEALLOCATE(spack%press)
    DEALLOCATE(spack%vapp)
    DEALLOCATE(spack%volmol)
    !DEALLOCATE(spack%xlw)
    !DEALLOCATE(spack)

!!$  !Deallocating blocks data
!!$  CALL AllocIndex(block_size,mynum,mmxp,mmyp,mmzp,mia,miz,mja,mjz,mibcon,mi0,mj0,&
!!$         ngrids,dtlongn)

  END SUBROUTINE chem_qssa

  !========================================================================================
  SUBROUTINE solveur(prod,loss,specie,nspecies,dtlt,maxblock_size,ijkbegin,ijkend,itchim,mynum)
  !========================================================================================

    DOUBLE PRECISION, INTENT(IN)    :: prod(maxblock_size,nspecies)
    DOUBLE PRECISION, INTENT(IN)    :: loss(maxblock_size,nspecies)
    DOUBLE PRECISION, INTENT(INOUT) :: specie(maxblock_size,nspecies)
    INTEGER,          INTENT(IN)    :: nspecies
    REAL,             INTENT(IN)    :: dtlt
    INTEGER,          INTENT(IN)    :: maxblock_size
    INTEGER,          INTENT(IN)    :: ijkbegin
    INTEGER,          INTENT(IN)    :: ijkend
    INTEGER,          INTENT(IN)    :: itchim
    INTEGER,          INTENT(IN)    :: mynum

    DOUBLE PRECISION :: delt
    INTEGER :: ijk,itt,i
    !       boucle solveur QSSA
    !       !!!!!!!!!!!!!!!!!!!!!

    delt=DBLE(dtlt/itchim)

    DO itt=1,itchim
       CALL qssa(nspecies,specie,prod,loss,delt,maxblock_size,ijkbegin,ijkend)

    END DO

    DO i=1,nspecies
       DO ijk=ijkbegin,ijkend
          IF (specie(ijk,I)<0.d0) THEN
             specie(ijk,i)=0.d0
          END IF

          !     IF ((specie(ijk,i).le.-1.d-10).or.(specie(ijk,i).ge.1.d+20)) THEN
          IF (                               specie(ijk,i).GE.1.d+20 )  THEN
             PRINT*,'Wrong value for specie !!',ijk,i,specie(ijk,i),prod(ijk,i),loss(ijk,i),delt,dtlt
             PRINT*,'mynum=',mynum
             STOP
          END IF
       END DO
    END DO

  END SUBROUTINE solveur

  !========================================================================================
  SUBROUTINE qssa(nspecies,specie,prod,loss,delt,maxblock_size,ijkbegin,ijkend)
  !========================================================================================

    INTEGER                                              , INTENT(IN)    :: maxblock_size
    INTEGER                                              , INTENT(IN)    :: ijkbegin
    INTEGER                                              , INTENT(IN)    :: ijkend
    INTEGER                                              , INTENT(IN)    :: nspecies
    DOUBLE PRECISION                                     , INTENT(IN)    :: delt
    DOUBLE PRECISION , DIMENSION(maxblock_size,nspecies) , INTENT(IN)    :: loss
    DOUBLE PRECISION , DIMENSION(maxblock_size,nspecies) , INTENT(IN)    :: prod
    DOUBLE PRECISION , DIMENSION(maxblock_size,nspecies) , INTENT(INOUT) :: specie

    !- local var
    DOUBLE PRECISION, DIMENSION(ijkbegin:ijkend,nspecies) :: cg1,cg2
    INTEGER ic,ijk
    DOUBLE PRECISION ybid,pr,pr1,pr2
!    DOUBLE PRECISION, EXTERNAL :: cvmgm_dp,cvmgp_dp,cvmgz_dp
    !c **** resolution du systeme par methode de Chang 1987 fonction
    !c **** ctemps caracteristique de destruction et pas de ctemps
    !c ****
    ! les fonctions cvmgm et cvmgp sont dans lib/rfvec.f90

    DO ic=1,nspecies
       DO ijk=ijkbegin,ijkend
          ybid=cvmgm_dp(loss(ijk,ic),1.d0,-loss(ijk,ic))

          cg1(ijk,ic)=specie(ijk,ic)+delt*(prod(ijk,ic)-loss(ijk,ic)*specie(ijk,ic))

          pr=prod(ijk,ic)/ybid

          pr1=pr+(specie(ijk,ic)-pr)*dexp(-ybid*delt)

          pr2=ybid*delt-10.d0

          cg2(ijk,ic)=cvmgp_dp(pr,pr1,pr2)

          specie(ijk,ic)=cvmgp_dp(cg1(ijk,ic),cg2(ijk,ic),-(100.d0*loss(ijk,ic)*delt-1.0D0 ))

          specie(ijk,ic)=cvmgz_dp(cg1(ijk,ic),cg2(ijk,ic),loss(ijk,ic))

       END DO
    END DO

  END SUBROUTINE QSSA




  ! +------------------------------------------------------------------+
  !
  !       Return VCT1 if VCT3 <= 0., else VCT2.
  !
  !========================================================================================
  FUNCTION cvmgm_dp(vct1,vct2,vct3)
    !========================================================================================

    DOUBLE PRECISION , INTENT(IN)  :: vct1
    DOUBLE PRECISION , INTENT(IN)  :: vct2
    DOUBLE PRECISION , INTENT(IN)  :: vct3

    DOUBLE PRECISION :: cvmgm_dp

    IF(vct3.LT.0.D0)THEN
       cvmgm_dp=vct1
    ELSE
       cvmgm_dp=vct2
    ENDIF

  END FUNCTION cvmgm_dp
  ! +------------------------------------------------------------------+

  ! +------------------------------------------------------------------+
  !
  !       Return VCT1 if VCT3 => 0., else VCT2.
  !
  FUNCTION cvmgp_dp(vct1,vct2,vct3)

    DOUBLE PRECISION , INTENT(IN) :: vct1
    DOUBLE PRECISION , INTENT(IN) :: vct2
    DOUBLE PRECISION , INTENT(IN) :: vct3

    DOUBLE PRECISION :: cvmgp_dp

    IF(vct3.GE.0.D0)THEN
       cvmgp_dp=vct1
    ELSE
       cvmgp_dp=vct2
    ENDIF

  END FUNCTION cvmgp_dp

  ! +------------------------------------------------------------------+
  !
  !       Return VCT1 if VCT3 = 0., else VCT2.
  !
  !========================================================================================
  FUNCTION cvmgz_dp(VCT1,VCT2,VCT3)
    !========================================================================================
    DOUBLE PRECISION , INTENT(IN) :: vct1
    DOUBLE PRECISION , INTENT(IN) :: vct2
    DOUBLE PRECISION , INTENT(IN) :: vct3

    DOUBLE PRECISION :: cvmgz_dp

    IF(vct3.EQ.0.D0)THEN
       cvmgz_dp=vct1
    ELSE
       cvmgz_dp=vct2
    END IF
  END FUNCTION cvmgz_dp
  ! +------------------------------------------------------------------+
  !

END MODULE mod_chem_spack_qssa
