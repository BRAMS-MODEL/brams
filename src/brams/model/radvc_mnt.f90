!----------------------------------------------------------------------!
! Optional advection scheme for CCATT-BRAMS/BRAMS models version 4.2+  !
! Based on Walcek, 2000 (JGR) and Walcek and Aleksic, 1998 (ATENV).    !
! The scheme is highly conservative, monotonic and keeps mass mixing   !
! ratio positive definite. 					       !
! Implemented by Saulo Freitas (saulo.freitas@cptec.inpe.br) @ Jun/2009!
! MPI/Paralelized by L. Flavio/J. Panneta                              !
!----------------------------------------------------------------------!

MODULE monotonic_adv

  USE node_mod, only:           &
          	      ibcon,    &  !INTENT(IN)
          	      mynum,    &  !INTENT(IN)
          	      nodei0,   &  !INTENT(IN)
          	      nodej0,   &  !INTENT(IN)
          	      nodemyp,  &  !INTENT(IN)
          	      nodemxp,  &  !INTENT(IN)
          	      nodemzp      !INTENT(IN)
  USE mem_grid, ONLY:        &
  	    	     dtlt,   & !INTENT(IN)
		     time,   &
  	    	     ngrids, & !INTENT(IN)
  	    	     ngrid,  & !INTENT(IN)
	    	     dzt,    & !INTENT(IN)
	    	     dztn,   & !INTENT(IN)
            	     grid_g, & !INTENT(IN)
	    	     grid_g, & !INTENT(IN)
	    	     naddsc, & !INTENT(IN)
	             hw4   , & !INTENT(IN)
		     if_adap,& !INTENT(IN)
		     nzpmax, & !INTENT(IN)
		     dyncore_flag  !INTENT(IN)

  USE mem_basic, ONLY: basic_g  !INTENT(IN)

  USE micphys    ,  ONLY: level !INTENT(IN)

  USE rconstants ,  ONLY: cp,p00,cv,rgas,cpi   !INTENT(IN)

  use mem_aer1, only: &
       aerosol,    &       !INTENT(IN)
       num_scalar_aer_1st !INTENT(IN)

  use mem_chem1, only: &
       nspecies_transported !INTENT(IN)

  use module_dry_dep, only: &
       dd_sedim,            &
       naer_transported,    &
       sedim_type

  use mem_scratch, only  : scratch  ! only scr1, inout

  use var_tables, only : scalar_tab & ! (var_p = IN, var_t = INOUT)
                        ,num_scalar   ! (IN)

  use advMessageMod, ONLY: SendMessageI,RecvMessageI,SendMessageJ,RecvMessageJ, &
                           newM2,newM3,newIa,newIz,newJa,newJz,nRecvI,nRecvJ,nSendI,nSendJ, &
			   totalrecvi,totalsendi,totalrecvj,totalsendj

  USE ParLib, ONLY: parf_send_noblock_real, parf_get_noblock_real, &
                    parf_wait_any_nostatus, &
                    parf_wait_all_nostatus

  USE ModNamelistFile, ONLY: NamelistFile

  use ccatt_start, only: &
       ccatt               ! (IN)

  IMPLICIT NONE

  INTEGER:: advmnt

  INTEGER , PARAMETER :: ON=1,OFF=0

  INTEGER :: GhostZoneLength

  INTEGER , PARAMETER :: use_true_density  = 1 ! 0= OFF, 1=ON

  !- for theoretical experiments
  INTEGER , PARAMETER :: theor_wind = 0        ! 0= OFF, 1=ON

  INTEGER           :: mnt_adv_jnitialized=0

  real, parameter :: c1 = cv/rgas, c2 = p00/rgas !c2 = p00*(cpi**c1)/rgas

  TYPE advmnt_vars
     REAL,POINTER,DIMENSION  (:,:,:)  :: u3d
     REAL,POINTER,DIMENSION  (:,:,:)  :: v3d
     REAL,POINTER,DIMENSION  (:,:,:)  :: w3d

     REAL,POINTER,DIMENSION  (:,:,:)  :: vc3d_in
     REAL,POINTER,DIMENSION  (:,:,:)  :: vc3d_out
     REAL,POINTER,DIMENSION  (:,:,:)  :: vc3d_x
     REAL,POINTER,DIMENSION  (:,:,:)  :: vc3d_y

     REAL,POINTER,DIMENSION  (:,:,:)  :: dd0_3d
     REAL,POINTER,DIMENSION  (:,:,:)  :: dd0_3du
     REAL,POINTER,DIMENSION  (:,:,:)  :: dd0_3dv
     REAL,POINTER,DIMENSION  (:,:,:)  :: dd0_3dw

     REAL,POINTER,DIMENSION  (:,:,:)  :: den0_3d
     REAL,POINTER,DIMENSION  (:,:,:)  :: den1_3d
     REAL,POINTER,DIMENSION  (:,:,:)  :: den2_3d
     REAL,POINTER,DIMENSION  (:,:,:)  :: den3_3d

     REAL,POINTER,DIMENSION  (:,:,:)  :: l_dxtW
     REAL,POINTER,DIMENSION  (:,:,:)  :: l_dytW

     REAL,POINTER,DIMENSION  (:,:)    :: dxtW
     REAL,POINTER,DIMENSION  (:,:)    :: dytW
     REAL,POINTER,DIMENSION  (:)      :: dztW

  END TYPE advmnt_vars

  TYPE(advmnt_vars), ALLOCATABLE,DIMENSION(:)   :: advmnt_g

  PUBLIC :: advmnt_driver  ! Subroutine

  INTEGER :: nSend_i,nSend_j,nRecv_i,nRecv_j
  INTEGER,PARAMETER :: bigdump=1
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: buffcomm
  INTEGER :: nRec_i,nSnd_i,nRec_j,nSnd_j
  INTEGER :: bufSendTotalLength_i,bufSendTotalLength_j
  INTEGER :: bufREcvTotalLength_i,bufRecvTotalLength_j
  real, allocatable :: bufRecv(:)
  real, allocatable :: bufSend(:)

CONTAINS
  !----------------------------------------------------

  SUBROUTINE advmnt_driver(varn,m1 ,m2 ,m3 ,ia,iz,ja,jz,izu,jzv,mynum)
    !use Extras            , only: extra3d,extra2d,na_EXTRA3D

    IMPLICIT NONE
    INTEGER , INTENT(IN) :: m1
    INTEGER , INTENT(IN) :: m2
    INTEGER , INTENT(IN) :: m3
    INTEGER , INTENT(IN) :: ia
    INTEGER , INTENT(IN) :: iz
    INTEGER , INTENT(IN) :: ja
    INTEGER , INTENT(IN) :: jz
    INTEGER , INTENT(IN) :: izu
    INTEGER , INTENT(IN) :: jzv
    INTEGER , INTENT(IN) :: mynum
    character(len=*),intent(IN) :: varn

    !--- local vars
    integer n,ng,mxyzp,i,j,procfile,ibegin,iend,jbegin,jend,i_scl
    real, pointer :: scalarp, scalart
    integer :: sori,sorj,sosi,sosj,current_aer_ispc,current_ndt_z
    integer, dimension(naer_transported) :: ndt_z
    LOGICAL  :: IsThisScalarAer =.false.

    IF(mnt_adv_jnitialized == OFF) THEN
       if(mynum == 0) stop 'ADV MNT called with mynum = 0, try np = 2'
       CALL initialize_advmnt(ngrids,nodemzp(mynum,:), &
                              nodemxp(mynum,:),nodemyp(mynum,:))

       sori=maxval(totalrecvi)
       sorj=maxval(totalrecvj)
       sosi=maxval(totalsendi)
       sosj=maxval(totalsendj)
       allocate(bufRecv(max(sori,sorj)))
       allocate(bufSend(max(sosi,sosj)))

       !WRITE(*,FMT='(A,I3.3,A)') 'Allocated arrays for processor ',mynum,'.'; CALL flush(6)
       do ng=1,ngrids
         iBegin=newIa(ng)-1
         iEnd=newIz(ng)+1
         jBegin=newJa(ng)-1
         jEnd=newJz(ng)+1
         CALL initialize_grid_spacings(ng,nodemzp(mynum,ng), &
	               nodemxp(mynum,ng),nodemyp(mynum,ng) &
                      ,grid_g(ng)%dxt	    &
                      ,grid_g(ng)%dyt	    &
 		      ,grid_g(ng)%fmapt     &
		      ,grid_g(ng)%rtgt      &
!
    		      ,advmnt_g(ng)%dxtW(iBegin:iEnd,jBegin:jEnd) &
    		      ,advmnt_g(ng)%dytW(iBegin:iEnd,jBegin:jEnd) &
		      ,advmnt_g(ng)%dztW )
       enddo

       if(use_true_density == OFF) then
          do ng=1,ngrids
            iBegin=newIa(ng)-1
            iEnd=newIz(ng)+1
            jBegin=newJa(ng)-1
            jEnd=newJz(ng)+1
	    CALL initialize_densities(nodemzp(mynum,ng),&
	               nodemxp(mynum,ng),nodemyp(mynum,ng) &
                      , basic_g(ng)%dn0     &
		      , basic_g(ng)%dn0u    &
                      , basic_g(ng)%dn0v    &
                      ,advmnt_g(ng)%dd0_3d (1:m1,iBegin:iEnd,jBegin:jEnd)  &
    		      ,advmnt_g(ng)%dd0_3du(1:m1,iBegin:iEnd,jBegin:jEnd) &
    		      ,advmnt_g(ng)%dd0_3dv(1:m1,iBegin:iEnd,jBegin:jEnd) &
		      ,advmnt_g(ng)%dd0_3dw(1:m1,iBegin:iEnd,jBegin:jEnd) )
	   enddo
       endif

       mnt_adv_jnitialized= ON
    END IF

    mxyzp=m1*m2*m3

    !- This scheme is not applied to advect  U, V, and W
    if (varn .eq. 'V' .or. varn .eq. 'ALL') then
      stop 'not using mnt to advect u,v,w'
    endif
    if (if_adap /= 0) then
      stop 'MNT advection not ready for shaved eta'
    endif
    !
    ndt_z =0 ! integer initialization
    !
    !- Advect  scalars
    iBegin=newIa(ngrid)-1
    iEnd  =newIz(ngrid)+1
    jBegin=newJa(ngrid)-1
    jEnd  =newJz(ngrid)+1

    !- get actual air densities, if using them instead of basic state fields
    if(use_true_density == ON) then
     call get_true_densities(m1,m2,m3,level &
                      , basic_g(ngrid)%rtp     &
		      , basic_g(ngrid)%rv      &
                      , basic_g(ngrid)%pp      &
		      , basic_g(ngrid)%pi0     &
                      , basic_g(ngrid)%theta   &
                      ,advmnt_g(ngrid)%dd0_3d (1:m1,iBegin:iEnd,jBegin:jEnd)  &
    		      ,advmnt_g(ngrid)%dd0_3du(1:m1,iBegin:iEnd,jBegin:jEnd) &
    		      ,advmnt_g(ngrid)%dd0_3dv(1:m1,iBegin:iEnd,jBegin:jEnd) &
		      ,advmnt_g(ngrid)%dd0_3dw(1:m1,iBegin:iEnd,jBegin:jEnd) )

    endif

    !- prepare wind velocities including map factors
    call prepare_winds(dtlt,m1,m2,m3,ia,iz,ja,jz     &
                      ,basic_g(ngrid)%uc  &
		      ,basic_g(ngrid)%up  &
                      ,basic_g(ngrid)%vc  &
                      ,basic_g(ngrid)%vp  &
		      ,basic_g(ngrid)%wc  &
                      ,basic_g(ngrid)%wp  &
!
 		      ,grid_g(ngrid)%fmapui &
		      ,grid_g(ngrid)%fmapvi &
		      ,grid_g(ngrid)%rtgt   &
                      ,grid_g(ngrid)%rtgu   &
		      ,grid_g(ngrid)%rtgv   &
		      ,grid_g(ngrid)%f13t   &
                      ,grid_g(ngrid)%f23t   &
!
                      ,advmnt_g(ngrid)%u3d(1:m1,iBegin:iEnd,jBegin:jEnd)  &
    		      ,advmnt_g(ngrid)%v3d(1:m1,iBegin:iEnd,jBegin:jEnd)  &
    		      ,advmnt_g(ngrid)%w3d(1:m1,iBegin:iEnd,jBegin:jEnd)  &
		      ,ndt_z                )


    if(theor_wind == on) then

        call prepare_theor_winds(dtlt,m1,m2,m3,ia,iz,ja,jz,time     &
                      ,advmnt_g(ngrid)%u3d(1:m1,iBegin:iEnd,jBegin:jEnd)  &
    		      ,advmnt_g(ngrid)%v3d(1:m1,iBegin:iEnd,jBegin:jEnd)  &
    		      ,advmnt_g(ngrid)%w3d(1:m1,iBegin:iEnd,jBegin:jEnd)  &
                      ,grid_g(ngrid)%dxt    &
                      ,grid_g(ngrid)%dyt    &
                      ,advmnt_g(ngrid)%dd0_3d (1:m1,iBegin:iEnd,jBegin:jEnd)  &
    		      ,advmnt_g(ngrid)%dd0_3du(1:m1,iBegin:iEnd,jBegin:jEnd) &
    		      ,advmnt_g(ngrid)%dd0_3dv(1:m1,iBegin:iEnd,jBegin:jEnd) &
		      ,advmnt_g(ngrid)%dd0_3dw(1:m1,iBegin:iEnd,jBegin:jEnd) )
    endif
    !- prepare Walcek's air densities
    call get_Walceks_densities(dtlt,m1,m2,m3 &
                      ,advmnt_g(ngrid)%u3d(1:m1,iBegin:iEnd,jBegin:jEnd)  &
    		      ,advmnt_g(ngrid)%v3d(1:m1,iBegin:iEnd,jBegin:jEnd)  &
    		      ,advmnt_g(ngrid)%w3d(1:m1,iBegin:iEnd,jBegin:jEnd)  &
                      ,advmnt_g(ngrid)%dd0_3d (1:m1,iBegin:iEnd,jBegin:jEnd)  &
    		      ,advmnt_g(ngrid)%dd0_3du(1:m1,iBegin:iEnd,jBegin:jEnd) &
    		      ,advmnt_g(ngrid)%dd0_3dv(1:m1,iBegin:iEnd,jBegin:jEnd) &
		      ,advmnt_g(ngrid)%dd0_3dw(1:m1,iBegin:iEnd,jBegin:jEnd) &
                      ,advmnt_g(ngrid)%den0_3d(1:m1,iBegin:iEnd,jBegin:jEnd) &
    		      ,advmnt_g(ngrid)%den1_3d(1:m1,iBegin:iEnd,jBegin:jEnd) &
    		      ,advmnt_g(ngrid)%den2_3d(1:m1,iBegin:iEnd,jBegin:jEnd) &
		      ,advmnt_g(ngrid)%den3_3d(1:m1,iBegin:iEnd,jBegin:jEnd) &
    		      ,advmnt_g(ngrid)%dxtW(iBegin:iEnd,jBegin:jEnd) &
    		      ,advmnt_g(ngrid)%dytW(iBegin:iEnd,jBegin:jEnd) &
		      ,advmnt_g(ngrid)%dztW &
                      ,grid_g(ngrid)%dxt    &
                      ,grid_g(ngrid)%dyt    )


    CALL InitialFieldsUpdate(ngrids,m1,m2,m3,newm2(ngrid),newm3(ngrid),ngrid,mynum, &
	   nrecvI(ngrid),RecvMessageI(ngrid)%proc,RecvMessageI(ngrid)%tag, &
	   RecvMessageI(ngrid)%ia,RecvMessageI(ngrid)%iz,RecvMessageI(ngrid)%ja,RecvMessageI(ngrid)%jz, &
	   RecvMessageI(ngrid)%start,RecvMessageI(ngrid)%mSize,TotalRecvI(ngrid), &
	   nSendI(ngrid),sendMessageI(ngrid)%proc,sendMessageI(ngrid)%tag, &
	   sendMessageI(ngrid)%ia,sendMessageI(ngrid)%iz,sendMessageI(ngrid)%ja,sendMessageI(ngrid)%jz, &
	   sendMessageI(ngrid)%start,sendMessageI(ngrid)%mSize,TotalSendI(ngrid))

    CALL InitialFieldsUpdate(ngrids,m1,m2,m3,newm2(ngrid),newm3(ngrid),ngrid,mynum, &
	  nrecvJ(ngrid),RecvMessageJ(ngrid)%proc,RecvMessageJ(ngrid)%tag, &
	  RecvMessageJ(ngrid)%ia,RecvMessageJ(ngrid)%iz,RecvMessageJ(ngrid)%ja,RecvMessageJ(ngrid)%jz, &
	  RecvMessageJ(ngrid)%start,RecvMessageJ(ngrid)%mSize,TotalRecvJ(ngrid), &
	  nSendJ(ngrid),sendMessageJ(ngrid)%proc,sendMessageJ(ngrid)%tag, &
	  sendMessageJ(ngrid)%ia,sendMessageJ(ngrid)%iz,sendMessageJ(ngrid)%ja,sendMessageJ(ngrid)%jz, &
	  sendMessageJ(ngrid)%start,sendMessageJ(ngrid)%mSize,TotalSendJ(ngrid))


     !- ready to do advection, loop over all scalars
     if(advmnt == 1) then
        i_scl=1                                            !- all scalars
     elseif(advmnt == 2) then
        i_scl=num_scalar(ngrid) - NSPECIES_TRANSPORTED +1  !- only chemical + aer species
     elseif(advmnt == 3) then
        i_scl=2                                            !- all scalars, but not theta_il
     endif

     !srf- do n=1,num_scalar(ngrid)     ! original
     do n=i_scl,num_scalar(ngrid)

        !- if RK or ABM3 scheme, THP/THC are not transported here
        if (dyncore_flag == 2 .or. dyncore_flag == 3) then
          if (scalar_tab(n,ngrid)%name == 'THC' .or. &
              scalar_tab(n,ngrid)%name == 'THP') cycle
        endif

       !srf - somente para gases e aerossois
       !     do n=num_scalar(ngrid) - NSPECIES_TRANSPORTED +1,num_scalar(ngrid)
       !      if (scalar_tab(n,ngrid)%name /= 'COP' .and. scalar_tab(n,ngrid)%name /= 'CH4P') cycle
       !          scalar_tab(n,ngrid)%name /= 'O3P'  ) cycle

     !- Aerosol sedimentation
       IsThisScalarAer  = .false.
       current_aer_ispc = 0
       current_ndt_z    = 1
       if(ccatt == 1 .and. aerosol > 0 .and. n >= num_scalar_aer_1st) then
              !srf-  We are going to include sedimentation of aerosols at
              !      vertical advection tendency. It is supposed that scalars
              !      with  N >= num_scalar_aer_1st are _all_ aerosols .
              !
	      IsThisScalarAer=.true.
              current_aer_ispc = n - num_scalar_aer_1st + 1
	      current_ndt_z    = ndt_z (current_aer_ispc)

       endif

      scalarp => scalar_tab(n,ngrid)%var_p
      scalart => scalar_tab(n,ngrid)%var_t

!for intel - descomente
      call atob(mxyzp,scalarp,advmnt_g(ngrid)%vc3d_in(1:m1,iBegin:iEnd,jBegin:jEnd))

      call advect_mnt(ngrid,m1,m2,m3,ia,iz,ja,jz,dtlt,mynum,n, &
                      current_aer_ispc,current_ndt_z,IsThisScalarAer)

      call advtndc(m1,m2,m3,ia,iz,ja,jz        &
                  ,scalarp  ,advmnt_g(ngrid)%vc3d_out(1:m1,iBegin:iEnd,jBegin:jEnd)  &
                  ,scalart  ,dtlt,mynum        )
     enddo

  END SUBROUTINE advmnt_driver
 !--------------------------------------------------------------------------------

  SUBROUTINE initialize_advmnt(ngrids, mmzp,mmxp,mmyp)
   implicit none
   INTEGER , INTENT(IN) :: ngrids
   INTEGER , INTENT(IN) :: mmxp(ngrids)
   INTEGER , INTENT(IN) :: mmyp(ngrids)
   INTEGER , INTENT(IN) :: mmzp(ngrids)

   INTEGER :: ng

   !maxz=maxval(mmzp(1:ngrids))
   !maxx=maxval(mmxp(1:ngrids))
   !maxy=maxval(mmyp(1:ngrids))

   IF(ALLOCATED(advmnt_g)) THEN
      PRINT *,'Error in initialize_advmnt, sub: radvc_mnt: advmnt_g already allocated!'
      STOP 1000
   END IF

   ALLOCATE (advmnt_g(ngrids))
   do ng=1,ngrids

      ALLOCATE(advmnt_g(ng)%u3d    (mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%u3d=0.
      ALLOCATE(advmnt_g(ng)%v3d    (mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%v3d=0.
      ALLOCATE(advmnt_g(ng)%w3d    (mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%w3d=0.

      ALLOCATE(advmnt_g(ng)%dd0_3d (mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%dd0_3d =0.
      ALLOCATE(advmnt_g(ng)%dd0_3du(mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%dd0_3du=0.
      ALLOCATE(advmnt_g(ng)%dd0_3dv(mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%dd0_3dv=0.
      ALLOCATE(advmnt_g(ng)%dd0_3dw(mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%dd0_3dw=0.

      ALLOCATE(advmnt_g(ng)%den0_3d(mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%den0_3d=0.
      ALLOCATE(advmnt_g(ng)%den1_3d(mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%den1_3d=0.
      ALLOCATE(advmnt_g(ng)%den2_3d(mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%den2_3d=0.
      ALLOCATE(advmnt_g(ng)%den3_3d(mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%den3_3d=0.


      ALLOCATE(advmnt_g(ng)%l_dxtW(mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%l_dxtW=0.
      ALLOCATE(advmnt_g(ng)%l_dytW(mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%l_dytW=0.


      ALLOCATE(advmnt_g(ng)%dxtW   (         newM2(ng),newM3(ng))); advmnt_g(ng)%dxtW=0.
      ALLOCATE(advmnt_g(ng)%dytW   (         newM2(ng),newM3(ng))); advmnt_g(ng)%dytW=0.
      ALLOCATE(advmnt_g(ng)%dztW   (mmzp(ng)                    )); advmnt_g(ng)%dztW=0.


      ALLOCATE(advmnt_g(ng)%vc3d_in (mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%vc3d_in =0.
      ALLOCATE(advmnt_g(ng)%vc3d_out(mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%vc3d_out=0.
      !ALLOCATE(advmnt_g(ng)%vc3d_x (mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%vc3d_x  =0.
      !ALLOCATE(advmnt_g(ng)%vc3d_y (mmzp(ng),newM2(ng),newM3(ng))); advmnt_g(ng)%vc3d_y  =0.

   enddo

  END SUBROUTINE initialize_advmnt
 !--------------------------------------------------------------------------------

 SUBROUTINE Deinitialize_advmnt(ngrids, mmzp,mmxp,mmyp)
   implicit none
   INTEGER , INTENT(IN) :: ngrids
   INTEGER , INTENT(IN) :: mmxp(ngrids)
   INTEGER , INTENT(IN) :: mmyp(ngrids)
   INTEGER , INTENT(IN) :: mmzp(ngrids)

   INTEGER :: ng

   !maxz=maxval(mmzp(1:ngrids))
   !maxx=maxval(mmxp(1:ngrids))
   !maxy=maxval(mmyp(1:ngrids))

   ! IF(ALLOCATED(advmnt_g)) THEN
   !    PRINT *,'Error in initialize_advmnt, sub: radvc_mnt: advmnt_g already allocated!'
   !    RETURN
   ! END IF

   do ng=1,ngrids

      deallocate(advmnt_g(ng)%u3d)
      deallocate(advmnt_g(ng)%v3d)
      deallocate(advmnt_g(ng)%w3d)

      deallocate(advmnt_g(ng)%dd0_3d )
      deallocate(advmnt_g(ng)%dd0_3du)
      deallocate(advmnt_g(ng)%dd0_3dv)
      deallocate(advmnt_g(ng)%dd0_3dw)

      deallocate(advmnt_g(ng)%den0_3d)
      deallocate(advmnt_g(ng)%den1_3d)
      deallocate(advmnt_g(ng)%den2_3d)
      deallocate(advmnt_g(ng)%den3_3d)

      deallocate(advmnt_g(ng)%dxtW)
      deallocate(advmnt_g(ng)%dytW)
      deallocate(advmnt_g(ng)%dztW)

      deallocate(advmnt_g(ng)%l_dxtW)
      deallocate(advmnt_g(ng)%l_dytW)

      deallocate(advmnt_g(ng)%vc3d_in )
      deallocate(advmnt_g(ng)%vc3d_out)
      !deallocate(advmnt_g(ng)%vc3d_x)
      !deallocate(advmnt_g(ng)%vc3d_y)

   enddo
   deallocate (advmnt_g)

  END SUBROUTINE Deinitialize_advmnt

 !----------------------------------------------------

  SUBROUTINE initialize_densities(m1,m2,m3,dn0,dn0u,dn0v &
			    ,dd0_3d ,dd0_3du,dd0_3dv,dd0_3dw )
   IMPLICIT NONE
   INTEGER , INTENT(IN) :: m1,m2,m3
   REAL,DIMENSION(m1,m2,m3),intent(IN) ::dn0,dn0u,dn0v
   REAL,DIMENSION(m1,m2,m3),intent(OUT)::dd0_3d,dd0_3du,dd0_3dv,dd0_3dw
   ! local var
   integer i,j,k

   dd0_3d (:,:,:)=  dn0 (:,:,:)
   dd0_3du(:,:,:)=  dn0u(:,:,:)
   dd0_3dv(:,:,:)=  dn0v(:,:,:)
   do j = 1,m3
      do i = 1,m2
         do k = 1,m1-1
            dd0_3dw(k,i,j) = 0.5*(dn0(k,i,j) +dn0(k+1,i,j))
         enddo
	 dd0_3dw(m1,i,j)=dd0_3dw(m1-1,i,j)
   enddo;enddo


  END SUBROUTINE initialize_densities
 !----------------------------------------------------

  SUBROUTINE initialize_grid_spacings(ng,m1,m2,m3,dxt,dyt,fmapt,rtgt &
			    ,dxtW,dytW,dztW )
   IMPLICIT NONE
   INTEGER , INTENT(IN) :: ng,m1,m2,m3
   REAL,DIMENSION(m2,m3),intent(IN) :: dxt,dyt,fmapt,rtgt

   REAL,DIMENSION(m2,m3),intent(OUT):: dxtW,dytW
   REAL,DIMENSION(m1),intent(OUT):: dztW
   ! local var
   integer i,j,k
   real rtgti

   do j = 1,m3
      do i = 1,m2
   	rtgti = 1. / rtgt(i,j)

	!- at init/rams_grid.f90:
        !     dxt(i,j)=fmapt(i,j)/(xmn(i,ngrid)-xmn(i-1,ngrid))
        !     dyt(i,j)=fmapt(i,j)/(ymn(j,ngrid)-ymn(j-1,ngrid))

   	 dxtW(i,j) = 1./(dxt(i,j) * fmapt(i,j) * rtgti)
   	 dytW(i,j) = 1./(dyt(i,j) * fmapt(i,j) * rtgti)
   enddo;enddo
   do k = 1,m1
    !- at init/gridset.f90:
    !  dztn(k,ifm) = 1. / (zmn(k,ifm) - zmn(k-1,ifm))
    ! Por que o Jacobiano nao depende de Z, o dztw depende somente
    ! de z.
    !dztW(k,i,j) = 1./ ( dzt(k) * rtgti * fmapt(i,j)**2 )
     dztW(k)	 = 1./ ( dztn(k,ng) ) !

   enddo
  END SUBROUTINE initialize_grid_spacings
 !----------------------------------------------------

  SUBROUTINE get_true_densities(m1,m2,m3,level,rtp,rv,pp,pi0,theta &
			       ,dd0_3d ,dd0_3du,dd0_3dv,dd0_3dw )
   IMPLICIT NONE
   INTEGER , INTENT(IN) :: m1,m2,m3,level
   REAL,DIMENSION(m1,m2,m3),intent(IN) ::rtp,rv,pp,pi0,theta
   REAL,DIMENSION(:,:,:),intent(OUT):: dd0_3d
   REAL,DIMENSION(m1,m2,m3),intent(OUT)  :: dd0_3du,dd0_3dv,dd0_3dw
   ! local var
   integer i,j,k
   real c3

   c3 = c2 * (cpi**c1)

   !- true air density at points "T"

   if( level == 0 ) then
     dd0_3d(:,:,:) = (c3/theta(:,:,:))*(pi0(:,:,:)+pp(:,:,:))**c1
   else
   do j = 1,m3
      do i = 1,m2
         do k = 1,m1
            dd0_3d(k,i,j) = (c3/theta(k,i,j))* (1. + rtp(k,i,j))/ &
	            (1. + 1.61*rv(k,i,j))*(pi0(k,i,j)+pp(k,i,j))**c1
         end do
      end do
   end do
   endif

   !- true air density at points "U", "V" and "W":

   call fill_dn0uv(m1,m2,m3,dd0_3d,dd0_3du,dd0_3dv)

   do j = 1,m3
      do i = 1,m2
         do k = 1,m1-1
            dd0_3dw(k,i,j) = 0.5*(dd0_3d(k,i,j) + dd0_3d(k+1,i,j))
         enddo
	 dd0_3dw(m1,i,j)=dd0_3dw(m1-1,i,j)
   enddo;enddo

 END SUBROUTINE get_true_densities

 !----------------------------------------------------

 SUBROUTINE prepare_winds(dtlt,m1,m2,m3,ia,iz,ja,jz &
                            ,uc,up,vc,vp,wc,wp &
  		            ,fmapui &
		            ,fmapvi &
                            ,rtgt   &
                            ,rtgu   &
		            ,rtgv   &
		            ,f13t   &
                            ,f23t   &
			    ,u3d,v3d,w3d &
			    ,ndt_z                )

  implicit none
  INTEGER , INTENT(IN) :: m1,m2,m3,ia,iz,ja,jz
  REAL    , INTENT(IN) :: dtlt
  REAL,DIMENSION(m1,m2,m3),intent(IN) :: uc,up,vc,vp,wc,wp
  REAL,DIMENSION(m2,m3)   ,intent(IN) :: rtgt,rtgu,rtgv,fmapui,fmapvi,f13t,f23t

  REAL,DIMENSION(m1,m2,m3),intent(OUT)::u3d,v3d,w3d

  !- aerosol sedimentation
  INTEGER, DIMENSION(naer_transported) , INTENT(INOUT) :: ndt_z

  !- local var
  !real   dtlto2
  integer jm,jp,im,ip , ispc
  integer i,j,k
  real :: cx1,cx2,rtgti,dum(m1)

  ! dtlto2 = .5


  ! u3d, u3d, and w3d are input as the velocity components (averaged
  ! between past and current time levels) times dtlt.
   do j=1,m3
     do i = 1,m2
      do k = 1,m1

          w3d(k,i,j) = ( wc(k,i,j) + wp(k,i,j) )*0.5
          u3d(k,i,j) = ( uc(k,i,j) + up(k,i,j) )*0.5
          v3d(k,i,j) = ( vc(k,i,j) + vp(k,i,j) )*0.5

   enddo;enddo;enddo
  ! after this point w3d is the cartesian vertical velocity


  !return ! for pure cartesian coordinates

  ! here w3d is the cartesian vertical velocity

  ! Add contribution to w3d from horiz winds crossing sloping sigma surfaces,
  ! and include 1/rtgt factor in w3d
  do j = 1,m3
     jm = max(1,j-1)
     jp = min(m3,j+1)
     do i = 1,m2
  	im = max(1,i-1)
  	ip = min(m2,i+1)
        rtgti = 1. / rtgt(i,j)

	do k = 1,m1-1
	    w3d(k,i,j) = ( (u3d(k,i,j) + u3d(k+1,i,j) + u3d(k,im,j) + u3d(k+1,im,j) ) * f13t(i,j)  &
	               +   (v3d(k,i,j) + v3d(k+1,i,j) + v3d(k,i,jm) + v3d(k+1,i,jm) ) * f23t(i,j)  &
		         ) * hw4(k)  &
	               + w3d(k,i,j) * rtgti

  	enddo
     enddo
  enddo
  !- after this point w3d is the sigma_z velocity

  !- including map factors on U,V:
  do j = 1,m3
     do i = 1,m2
  	cx1 = fmapui(i,j) * rtgu(i,j)
  	cx2 = fmapvi(i,j) * rtgv(i,j)
  	do k = 1,m1-1
  	   u3d(k,i,j) = u3d(k,i,j) * cx1
  	   v3d(k,i,j) = v3d(k,i,j) * cx2

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	   !u3d(k,i,j) =  10.
	   !v3d(k,i,j) =  10.
           !w3d(k,i,j) =  0.
	   !wc(k,i,j)=w3d(k,i,j)
	   !uc(k,i,j)=u3d(k,i,j)
	   !vc(k,i,j)=v3d(k,i,j)
	   !wp(k,i,j)=w3d(k,i,j)
	   !up(k,i,j)=u3d(k,i,j)
	   !vp(k,i,j)=v3d(k,i,j)

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  	enddo
     enddo
  enddo
  !-----------------------------------------
  !- control for aerosol sedimentation
  if(AEROSOL > 0 .and. naer_transported > 0) then
    ! very crude estimation of CFL violation and fix for the number of sub-timesteps
    ! for large particles
    do ispc=1,naer_transported
      ndt_z(ispc)=ceiling(maxval(abs(dd_sedim(ispc,ngrid)%v_sed_part))*dtlt*maxval(dzt(1:m1)))
      !print*,'aer specie ndtz=',ispc,ndt_z(ispc),maxval(abs(dd_sedim(ispc,ngrid)%v_sed_part)),&
      !minval(abs(dd_sedim(ispc,ngrid)%v_sed_part))
    enddo
  endif
  !- end of aerosol sedimentation
  !-----------------------------------------


  END SUBROUTINE prepare_winds
 !----------------------------------------------------

  SUBROUTINE get_Walceks_densities(dt,m1,m2,m3,u3d,v3d,w3d &
			    ,dd0_3d ,dd0_3du,dd0_3dv,dd0_3dw &
			    ,den0_3d,den1_3d,den2_3d,den3_3d &
			    ,dxtW,dytW,dztW,dxt,dyt)

   IMPLICIT NONE
   !-in
   INTEGER , INTENT(IN) :: m1,m2,m3
   REAL    , INTENT(IN) :: dt
   REAL, DIMENSION(m1)      , INTENT(IN) :: dztW
   REAL, DIMENSION(m2,m3)   , INTENT(IN) :: dxtW,dytW

   REAL, DIMENSION(m2,m3)   , INTENT(IN) :: dxt,dyt

   REAL, DIMENSION(m1,m2,m3), INTENT(INOUT) :: u3d,v3d,w3d     &
			                   ,dd0_3du &
					   ,dd0_3dv,dd0_3dw
   REAL, DIMENSION(m1,m2,m3), INTENT(INOUT) :: dd0_3d
   !-out
   REAL, DIMENSION(m1,m2,m3), INTENT(OUT):: den0_3d,den1_3d &
                                           ,den2_3d,den3_3d

   ! local var
   integer i,j,k

! for homogeneous u,v,w ----
!dd0_3d=1. ; dd0_3dw=1. ; dd0_3du=1. ; dd0_3dv=1.
!----------------------------
    DO  j=m3,2,-1
      DO  i=2,m2
        DO k = 2,m1

            den0_3d(k,i,j)=dd0_3d(k,i,j)

            den1_3d(k,i,j)=den0_3d(k,i,j)- dt/dxtW(i,j)*(dd0_3du(k,i,j)*u3d(k,i,j)-dd0_3du(k,i-1,j)*u3d(k,i-1,j))

            den2_3d(k,i,j)=den1_3d(k,i,j)- dt/dytW(i,j)*(dd0_3dv(k,i,j)*v3d(k,i,j)-dd0_3dv(k,i,j-1)*v3d(k,i,j-1))

            den3_3d(k,i,j)=den2_3d(k,i,j)- dt/dztW(k)  *(dd0_3dw(k,i,j)*w3d(k,i,j)-dd0_3dw(k-1,i,j)*w3d(k-1,i,j))

       END DO
      END DO
    END DO
    !srf- BC for den3_3d
     den3_3d(:,1,:)=den3_3d(:,2,:)
     den3_3d(:,:,1)=den3_3d(:,:,2)


 END SUBROUTINE get_Walceks_densities
 !----------------------------------------------------
 SUBROUTINE advect_mnt(ngrid,m1,m2,m3,ia,iz,ja,jz,dt,mynum,n,&
                       current_aer_ispc,current_ndt_z,IsThisScalarAer)

  IMPLICIT NONE
  INTEGER , INTENT(IN) :: m1,ngrid
  INTEGER , INTENT(IN) :: m2
  INTEGER , INTENT(IN) :: m3
  INTEGER , INTENT(IN) :: ia
  INTEGER , INTENT(IN) :: iz
  INTEGER , INTENT(IN) :: ja
  INTEGER , INTENT(IN) :: jz,n
  INTEGER , INTENT(IN) :: mynum
  REAL    , INTENT(IN) :: dt
  INTEGER , INTENT(IN) :: current_ndt_z,current_aer_ispc
  LOGICAL , INTENT(IN) :: IsThisScalarAer
  !- local var
  !REAL,DIMENSION(m1)               :: dxx
  !REAL,DIMENSION(m2,m3)            :: dxy
  real masscon,initialmass,vol
  integer nrec,itz
  integer ibegin,iend,jbegin,jend
  !- type of sedimentation scheme (0= Walcek, 1=upwind)
  integer , parameter :: iupwind = 0

    iBegin= newIa(ngrid)-1
    iEnd  = newIz(ngrid)+1
    jBegin= newJa(ngrid)-1
    jEnd  = newJz(ngrid)+1

     !--- do X-advection
     call UpdateBorders(m1, newm2(ngrid), newm3(ngrid),advmnt_g(ngrid)%vc3d_in, &
	  nrecvI(ngrid), RecvMessageI(ngrid)%proc, RecvMessageI(ngrid)%tag, &
	  RecvMessageI(ngrid)%ia, RecvMessageI(ngrid)%iz, RecvMessageI(ngrid)%ja, RecvMessageI(ngrid)%jz, &
	  RecvMessageI(ngrid)%start, RecvMessageI(ngrid)%mSize, TotalRecvI(ngrid), &
	  nSendI(ngrid), sendMessageI(ngrid)%proc, sendMessageI(ngrid)%tag, &
	  sendMessageI(ngrid)%ia, sendMessageI(ngrid)%iz, sendMessageI(ngrid)%ja, sendMessageI(ngrid)%jz, &
	  sendMessageI(ngrid)%start, sendMessageI(ngrid)%mSize, TotalSendI(ngrid))

     call Advec3d_X(m1,newM2(ngrid),newM3(ngrid),2,newM2(ngrid)-1,2,newM3(ngrid)-1 &
                  ,advmnt_g(ngrid)%vc3d_in                             &
                  ,advmnt_g(ngrid)%u3d,advmnt_g(ngrid)%den0_3d         &
                  ,advmnt_g(ngrid)%den1_3d,dt,advmnt_g(ngrid)%dxtW     &
		  ,advmnt_g(ngrid)%dd0_3du                             &
		  ,advmnt_g(ngrid)%vc3d_out    ,mynum )

     !--- do Y-advection

     call UpdateBorders(m1, newm2(ngrid), newm3(ngrid),advmnt_g(ngrid)%vc3d_out, &
              nrecvJ(ngrid), RecvMessageJ(ngrid)%proc, RecvMessageJ(ngrid)%tag, &
              RecvMessageJ(ngrid)%ia, RecvMessageJ(ngrid)%iz, RecvMessageJ(ngrid)%ja, RecvMessageJ(ngrid)%jz, &
              RecvMessageJ(ngrid)%start, RecvMessageJ(ngrid)%mSize, TotalRecvJ(ngrid), &
              nSendJ(ngrid), sendMessageJ(ngrid)%proc, sendMessageJ(ngrid)%tag, &
              sendMessageJ(ngrid)%ia, sendMessageJ(ngrid)%iz, sendMessageJ(ngrid)%ja, sendMessageJ(ngrid)%jz, &
              sendMessageJ(ngrid)%start, sendMessageJ(ngrid)%mSize, TotalSendJ(ngrid))

     call Advec3d_Y(m1,newM2(ngrid),newM3(ngrid),2,newM2(ngrid)-1,2,newM3(ngrid)-1 &
                  ,advmnt_g(ngrid)%vc3d_out                                        &
		  ,advmnt_g(ngrid)%v3d,advmnt_g(ngrid)%den1_3d                     &
                  ,advmnt_g(ngrid)%den2_3d,dt,advmnt_g(ngrid)%dytW                 &
		  ,advmnt_g(ngrid)%dd0_3dv                                         &
		  ,advmnt_g(ngrid)%vc3d_in  ,mynum )

     !--- do k-advection
     call Advec3d_Z(m1,newM2(ngrid),newM3(ngrid),ibegin,iend          ,jbegin,jend &
	                   ,advmnt_g(ngrid)%vc3d_in                                   &
	                   ,advmnt_g(ngrid)%w3d,advmnt_g(ngrid)%den2_3d               &
			   ,advmnt_g(ngrid)%den3_3d,dt,advmnt_g(ngrid)%dztW           &
			   ,advmnt_g(ngrid)%dd0_3dw                                   &
			   ,advmnt_g(ngrid)%vc3d_out  ,mynum )


   !- aerosol section to include sedimentation
   !- the sedimentation process is done using pure cartesian coordinates
   !- so, all sedimentation velocities are treat as cartesian vertical velocities
   !- which are positive downwards.
   IF(AEROSOL > 0 .and. IsThisScalarAer) then

    !-srf introducing a time-splitting for aerosol sedimentation

    IF(iupwind == 0 ) then
      ! - Walcek method
      ! this routine works _only_ for mass concentration or density (kg/m3)
      ! converting mixing ratio (kg/kg) to density (kg/m3)
      advmnt_g(ngrid)%vc3d_in(:,:,:)=advmnt_g(ngrid)%vc3d_out(:,:,:) * advmnt_g(ngrid)%den0_3d(:,:,:)

      !- do time splitting for aerosols with large fall velocities
      do itz=1,current_ndt_z
        call Advec3d_Z_sedim(m1,m2,m3,ia,iz,ja,jz                        &
	           ,advmnt_g(ngrid)%vc3d_in(:,iBegin:iEnd,jBegin:jEnd)	 &
		   ,dd_sedim(current_aer_ispc,ngrid)%v_sed_part          & !fall velocity
		   ,dt/float(current_ndt_z)                              & !subtimestep
		   ,dzt(1:m1),grid_g(ngrid)%rtgt	                 &
		   ,advmnt_g(ngrid)%vc3d_out(:,iBegin:iEnd,jBegin:jEnd)  &
		   ,mynum )

        ! copy output to input array for the next sup-timestep
        if(itz < current_ndt_z) advmnt_g(ngrid)%vc3d_in(:,:,:)=advmnt_g(ngrid)%vc3d_out(:,:,:)

       enddo
       ! converting back mass concentration to mixing ratio
       advmnt_g(ngrid)%vc3d_out(:,:,:)=advmnt_g(ngrid)%vc3d_out(:,:,:)/advmnt_g(ngrid)%den0_3d(:,:,:)

    ELSEIF(iupwind == 1 ) then
       ! - upwind method
       !- do time splitting for aerosols with large fall velocities
       do itz=1,current_ndt_z
         call Advec3d_Z_sedim_upw(m1,m2,m3,ia,iz,ja,jz                          &
		   ,dd_sedim(current_aer_ispc,ngrid)%v_sed_part          & !fall velocity
		   ,dt/float(current_ndt_z)                              & !subtimestep
		   ,dzt(1:m1),grid_g(ngrid)%rtgt	                                 &
		   ,advmnt_g(ngrid)%vc3d_out(:,iBegin:iEnd,jBegin:jEnd)  &
		   ,mynum )

       enddo
     ENDIF

    ENDIF
  !- end of aerosol section

 END SUBROUTINE advect_mnt
 !--------------------------------------------------------------------

  SUBROUTINE prepare_theor_winds(dtlt,m1,m2,m3,ia,iz,ja,jz,time &
                            ,u3d,v3d,w3d&
			    ,dxt,dyt &
			    ,dd0_3d ,dd0_3du,dd0_3dv,dd0_3dw )


   IMPLICIT NONE
   INTEGER , INTENT(IN) :: m1,m2,m3,ia,iz,ja,jz
   REAL, INTENT(IN) :: dtlt,time
   REAL, DIMENSION(m2,m3)   , INTENT(IN) :: dxt,dyt

   REAL,DIMENSION(m1,m2,m3),intent(OUT)::u3d,v3d,w3d     &
                                        ,dd0_3d ,dd0_3du &
					,dd0_3dv,dd0_3dw

   !- local var
   REAL   :: dtlto2
   INTEGER :: i,j,k
   REAL  :: ai0s  =  25.0
   REAL  :: aj0s  =  50.0
   REAL  :: umx   =   80.0
   REAL,PARAMETER :: pii   =   3.141592653589793
   REAL    :: umax  =   0.0
   REAL    :: anrev,curnt,rx,xa,ilop,iwndty,nrec,ya
   REAL    :: periodo  =   6.*3600.

   dtlto2 =  10.!*dtlt

     !WRITE(6,*) ' Wind fields?  0-rotating; or  1-divergent winds'
     !iwndty = 0  ! 0-rotating
     !iwndty = 1  ! 1-divergent winds
     iwndty = 2  ! 1-divergent winds

     IF(iwndty==1) ai0s= 50.5
     ilop= ai0s-21.  ! needed for printouts
     ! Define wind fields (rotation or divergent) and initial mixing ratios
     !  Cone at (25,50) for rotating winds; Cone at (50,50) divergent winds
     DO k = 1,m1
   	DO  j=m3,1,-1
   	   DO  i=1,m2

   	     dd0_3d (k,i,j)=1.
	     dd0_3du(k,i,j)=1.
             dd0_3dv(k,i,j)=1.
	     dd0_3dw(k,i,j)=1.

   	     u3d(k,i,j)= -2.*umx*(REAL(j)-REAL(110)/2.-.5)/REAL(110)
   	     v3d(k,i,j)=  2.*umx*(REAL(i)-REAL(100)/2.-.5)/REAL(100)
 	     w3d(k,i,j)= 0.

   	    IF(iwndty==1) THEN
   	       xa=pii/25.
   	       IF(J>0) u3d(k,i,j)=umx*SIN(xa*REAL(i))*SIN(xa*(REAL(j)))
   	       IF(I>0) v3d(k,i,j)=umx*COS(xa*(REAL(i)-.5))*COS(XA*(REAL(j)+.5))
   	    END IF


   	    IF(iwndty==2) THEN
   	       xa=pii/100. ! m3=m2
   	       IF(J>0) u3d(k,i,j)=  umx* (SIN(xa*REAL(i)))**2 *SIN(2*xa*(REAL(j)))*cos(pii*time/periodo)
   	       IF(I>0) v3d(k,i,j)=- umx* (SIN(xa*REAL(j)))**2 *SIN(2*xa*(REAL(i)))*cos(pii*time/periodo)
   	    END IF


   	    IF(iwndty==3) THEN
   	        xa=pii/100. ! m3=m2
   	        ya=50.
	        IF(J>0) u3d(k,i,j)= -   umx* (SIN(    xa*REAL(i)))**2 * SIN(2.*xa*(REAL(j)-ya)) *cos(pii*time/periodo)
   	        IF(I>0) v3d(k,i,j)= 0.5*umx* (SIN(2.* xa*REAL(i)))    * COS(   xa*(REAL(j)-ya)) *cos(pii*time/periodo)
   	    END IF

   	    umax= MAX(ABS(u3d(k,i,j)),ABS(v3d(k,i,j)),umax)
   	    rx= SQRT((REAL(i)-ai0s)**2.+(REAL(j)-aj0s)**2.)

   	   END DO !i
   	END DO !j
     END DO !k

  END SUBROUTINE prepare_theor_winds
 !--------------------------------------------------------------------

  subroutine UpdateBorders(m1, m2, m3, field, &
	nRecv, procRecv_ext, tagRecv, iaRecv, izRecv, jaRecv, jzRecv, &
	bufRecvStart, bufRecvLength, bufRecvTotalLength, &
	nSend, procSend_ext, tagSend, iaSend, izSend, jaSend, jzSend, &
	bufSendStart, bufSendLength, bufSendTotalLength)
     integer, intent(in) :: m1
     integer, intent(in) :: m2
     integer, intent(in) :: m3
     real,    intent(inout) :: field(m1,m2,m3)
     integer, intent(in) :: nRecv
     integer, intent(in) :: procRecv_ext(nRecv)
     integer, intent(in) :: tagRecv(nRecv)
     integer, intent(in) :: iaRecv(nRecv)
     integer, intent(in) :: izRecv(nRecv)
     integer, intent(in) :: jaRecv(nRecv)
     integer, intent(in) :: jzRecv(nRecv)
     integer, intent(in) :: bufRecvStart(nRecv)
     integer, intent(in) :: bufRecvLength(nRecv)
     integer, intent(in) :: bufRecvTotalLength
     integer, intent(in) :: nSend
     integer, intent(in) :: procSend_ext(nSend)
     integer, intent(in) :: tagSend(nSend)
     integer, intent(in) :: iaSend(nSend)
     integer, intent(in) :: izSend(nSend)
     integer, intent(in) :: jaSend(nSend)
     integer, intent(in) :: jzSend(nSend)
     integer, intent(in) :: bufSendStart(nSend)
     integer, intent(in) :: bufSendLength(nSend)
     integer, intent(in) :: bufSendTotalLength

     integer :: procRecv(nRecv)
     integer :: procSend(nSend)

     integer :: i, j, iCnt, i1, i2, cnt, iRecv, iSend, ierr, iRecS
     integer :: reqRecv(nRecv)
     integer :: reqSend(nSend), recNum

     !print *, 'Inside update: ',procRecv_ext
     !print *, 'Inside update: ',procSend_ext; call flush(6)
     procRecv=procRecv_ext-1
     procSend=procSend_ext-1
     do iRecv = 1, nRecv
        call parf_get_noblock_real(bufRecv(bufRecvStart(iRecv)), bufRecvLength(iRecv),procRecv(iRecv) , &
	                           tagRecv(iRecv), reqRecv(iRecv))
	!call MPI_Irecv(bufRecv(bufRecvStart(iRecv)), bufRecvLength(iRecv), &
	!     MPI_REAL, &
	!     procRecv(iRecv), tagRecv(iRecv), MPI_COMM_WORLD, &
	!     reqRecv(iRecv), ierr)
     end do

     do iSend = 1, nSend
	i1 = bufSendStart(iSend)
	iCnt = bufSendStart(iSend)
	i2 = bufSendLength(iSend)
	do j = jaSend(iSend), jzSend(iSend)
	   do i = iaSend(iSend), izSend(iSend)
	      bufSend(iCnt:iCnt+m1-1) = field(1:m1,i,j)
	      iCnt = iCnt+m1
	   end do
	end do
	call parf_send_noblock_real(bufSend(i1), i2,procSend(iSend) , tagSend(iSend), reqSend(iSend))
	!call MPI_Isend(bufSend(i1), i2, MPI_REAL, &
	!     procSend(iSend), tagSend(iSend), MPI_COMM_WORLD, &
	!     reqSend(iSend), ierr)
     end do

!     do cnt = 1, nRecv
     do iRecV = 1, nRecv
        !WRITE (*,FMT='("Dump: ",2(I3.3,1X),Z9,1X,2(I3.3,1X))') mynum,iRecV,reqRecv(iRecV),iaRecv(iRecv),izRecv(iRecv); CALL flush(6)
 	call parf_wait_any_nostatus(nRecv,reqRecv,recNum)
!        call parf_wait_any_nostatus(nRecv,reqRecv,iRecv)
	!call MPI_Waitany(nRecv, reqRecv, iRecv, status, ierr)
	i1 = bufRecvStart(recNum)
	iCnt = bufRecvStart(recNum)
	i2 = bufRecvLength(recNum)
	do j = jaRecv(recNum), jzRecv(recNum)
	   do i = iaRecv(recNum), izRecv(recNum)
	      field(1:m1,i,j) = bufRecv(iCnt:iCnt+m1-1)
	      iCnt = iCnt+m1
	   end do
	end do
     end do
     call parf_wait_all_nostatus(nSend,reqSend)

   end subroutine UpdateBorders

 !--------------------------------------------------------------------
  SUBROUTINE InitialFieldsUpdate(ngrids,m1,m2,m3,Nm2,Nm3,ng,mynum, &
	nRec, procRecv, tagRecv, iaRecv, izRecv, jaRecv, jzRecv, &
	bufRecvStart, bufRecvLength, bufRecvTotalLength, &
	nSnd, procSend, tagSend, iaSend, izSend, jaSend, jzSend, &
	bufSendStart, bufSendLength, bufSendTotalLength)

     implicit none
     integer, intent(in) :: ngrids
     integer, intent(in) :: m1
     integer, intent(in) :: m2,nm2
     integer, intent(in) :: m3,nm3,ng,mynum
     integer, intent(in) :: nRec
     integer, intent(in) :: procRecv(nRec)
     integer, intent(in) :: tagRecv(nRec)
     integer, intent(in) :: iaRecv(nRec)
     integer, intent(in) :: izRecv(nRec)
     integer, intent(in) :: jaRecv(nRec)
     integer, intent(in) :: jzRecv(nRec)
     integer, intent(in) :: bufRecvStart(nRec)
     integer, intent(in) :: bufRecvLength(nRec)
     integer, intent(in) :: bufRecvTotalLength
     integer, intent(in) :: nSnd
     integer, intent(in) :: procSend(nSnd)
     integer, intent(in) :: tagSend(nSnd)
     integer, intent(in) :: iaSend(nSnd)
     integer, intent(in) :: izSend(nSnd)
     integer, intent(in) :: jaSend(nSnd)
     integer, intent(in) :: jzSend(nSnd)
     integer, intent(in) :: bufSendStart(nSnd)
     integer, intent(in) :: bufSendLength(nSnd)
     integer, intent(in) :: bufSendTotalLength

     INTEGER :: i,j,k

     Integer, save :: iupdate_dxy=0

    IF(bufSendTotalLength==0 .or. bufRecvTotalLength==0) RETURN

    !if(iupdate_dxy==0) then
     DO i=1,nm2
     	DO j=1,nm3
     	   DO k=1,1!m1
     	      advmnt_g(ng)%l_dxtW(k,i,j)=advmnt_g(ng)%dxtW(i,j)
     	      advmnt_g(ng)%l_dytW(k,i,j)=advmnt_g(ng)%dytW(i,j)
     	   END DO
     	END DO
     END DO
    !endif
    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%u3d, &
     nRec, procRecv, tagRecv, &
     iaRecv, izRecv, jaRecv, jzRecv, &
     bufRecvStart, bufRecvLength, bufRecvTotalLength, &
     nSnd, procSend, tagSend, &
     iaSend, izSend, jaSend, jzSend, &
     bufSendStart, bufSendLength, bufSendTotalLength)
!    WRITE (*,FMT='("Em init. Sou o ",I2.2,", na pos. ", I2.2)') myNum,1; CALL flush(6)

    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%v3d, &
     nRec, procRecv, tagRecv, &
     iaRecv, izRecv, jaRecv, jzRecv, &
     bufRecvStart, bufRecvLength, bufRecvTotalLength, &
     nSnd, procSend, tagSend, &
     iaSend, izSend, jaSend, jzSend, &
     bufSendStart, bufSendLength, bufSendTotalLength)
!    WRITE (*,FMT='("Em init. Sou o ",I2.2,", na pos. ", I2.2)') myNum,2; CALL flush(6)

!    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%w3d, &
!     nRec, procRecv, tagRecv, &
!     iaRecv, izRecv, jaRecv, jzRecv, &
!     bufRecvStart, bufRecvLength, bufRecvTotalLength, &
!     nSnd, procSend, tagSend, &
!     iaSend, izSend, jaSend, jzSend, &
!     bufSendStart, bufSendLength, bufSendTotalLength)

!    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%vc3d_x, &
!     nRec, procRecv, tagRecv, &
!     iaRecv, izRecv, jaRecv, jzRecv, &
!     bufRecvStart, bufRecvLength, bufRecvTotalLength, &
!     nSnd, procSend, tagSend, &
!     iaSend, izSend, jaSend, jzSend, &
!     bufSendStart, bufSendLength, bufSendTotalLength)
!    WRITE (*,FMT='("Em init. Sou o ",I2.2,", na pos. ", I2.2)') myNum,6; CALL flush(6)

!    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%vc3d_y, &
!     nRec, procRecv, tagRecv, &
!     iaRecv, izRecv, jaRecv, jzRecv, &
!     bufRecvStart, bufRecvLength, bufRecvTotalLength, &
!     nSnd, procSend, tagSend, &
!     iaSend, izSend, jaSend, jzSend, &
!     bufSendStart, bufSendLength, bufSendTotalLength)
!    WRITE (*,FMT='("Em init. Sou o ",I2.2,", na pos. ", I2.2)') myNum,7; CALL flush(6)

    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%dd0_3d, &
     nRec, procRecv, tagRecv, &
     iaRecv, izRecv, jaRecv, jzRecv, &
     bufRecvStart, bufRecvLength, bufRecvTotalLength, &
     nSnd, procSend, tagSend, &
     iaSend, izSend, jaSend, jzSend, &
     bufSendStart, bufSendLength, bufSendTotalLength)
!    WRITE (*,FMT='("Em init. Sou o ",I2.2,", na pos. ", I2.2)') myNum,8; CALL flush(6)

    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%dd0_3du, &
     nRec, procRecv, tagRecv, &
     iaRecv, izRecv, jaRecv, jzRecv, &
     bufRecvStart, bufRecvLength, bufRecvTotalLength, &
     nSnd, procSend, tagSend, &
     iaSend, izSend, jaSend, jzSend, &
     bufSendStart, bufSendLength, bufSendTotalLength)
!    WRITE (*,FMT='("Em init. Sou o ",I2.2,", na pos. ", I2.2)') myNum,9; CALL flush(6)

    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%dd0_3dv, &
     nRec, procRecv, tagRecv, &
     iaRecv, izRecv, jaRecv, jzRecv, &
     bufRecvStart, bufRecvLength, bufRecvTotalLength, &
     nSnd, procSend, tagSend, &
     iaSend, izSend, jaSend, jzSend, &
     bufSendStart, bufSendLength, bufSendTotalLength)
!    WRITE (*,FMT='("Em init. Sou o ",I2.2,", na pos. ", I2.2)') myNum,10; CALL flush(6)

    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%dd0_3dw, &
     nRec, procRecv, tagRecv, &
     iaRecv, izRecv, jaRecv, jzRecv, &
     bufRecvStart, bufRecvLength, bufRecvTotalLength, &
     nSnd, procSend, tagSend, &
     iaSend, izSend, jaSend, jzSend, &
     bufSendStart, bufSendLength, bufSendTotalLength)
!    WRITE (*,FMT='("Em init. Sou o ",I2.2,", na pos. ", I2.2)') myNum,11; CALL flush(6)

    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%den0_3d, &
     nRec, procRecv, tagRecv, &
     iaRecv, izRecv, jaRecv, jzRecv, &
     bufRecvStart, bufRecvLength, bufRecvTotalLength, &
     nSnd, procSend, tagSend, &
     iaSend, izSend, jaSend, jzSend, &
     bufSendStart, bufSendLength, bufSendTotalLength)
!    WRITE (*,FMT='("Em init. Sou o ",I2.2,", na pos. ", I2.2)') myNum,12; CALL flush(6)

    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%den1_3d, &
     nRec, procRecv, tagRecv, &
     iaRecv, izRecv, jaRecv, jzRecv, &
     bufRecvStart, bufRecvLength, bufRecvTotalLength, &
     nSnd, procSend, tagSend, &
     iaSend, izSend, jaSend, jzSend, &
     bufSendStart, bufSendLength, bufSendTotalLength)
!    WRITE (*,FMT='("Em init. Sou o ",I2.2,", na pos. ", I2.2)') myNum,13; CALL flush(6)

    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%den2_3d, &
     nRec, procRecv, tagRecv, &
     iaRecv, izRecv, jaRecv, jzRecv, &
     bufRecvStart, bufRecvLength, bufRecvTotalLength, &
     nSnd, procSend, tagSend, &
     iaSend, izSend, jaSend, jzSend, &
     bufSendStart, bufSendLength, bufSendTotalLength)
!    WRITE (*,FMT='("Em init. Sou o ",I2.2,", na pos. ", I2.2)') myNum,14; CALL flush(6)

    call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%den3_3d, &
     nRec, procRecv, tagRecv, &
     iaRecv, izRecv, jaRecv, jzRecv, &
     bufRecvStart, bufRecvLength, bufRecvTotalLength, &
     nSnd, procSend, tagSend, &
     iaSend, izSend, jaSend, jzSend, &
     bufSendStart, bufSendLength, bufSendTotalLength)
!    WRITE (*,FMT='("Em init. Sou o ",I2.2,", na pos. ", I2.2)') myNum,15; CALL flush(6)

    !if(iupdate_dxy==0) then
      call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%l_dxtW, &
       nRec, procRecv, tagRecv, &
       iaRecv, izRecv, jaRecv, jzRecv, &
       bufRecvStart, bufRecvLength, bufRecvTotalLength, &
       nSnd, procSend, tagSend, &
       iaSend, izSend, jaSend, jzSend, &
       bufSendStart, bufSendLength, bufSendTotalLength)
!      WRITE (*,FMT='("Em init. Sou o ",I2.2,", na pos. ", I2.2)') myNum,16; CALL flush(6)

      call UpdateBorders(m1, nm2, nm3,advmnt_g(ng)%l_dytW, &
       nRec, procRecv, tagRecv, &
       iaRecv, izRecv, jaRecv, jzRecv, &
       bufRecvStart, bufRecvLength, bufRecvTotalLength, &
       nSnd, procSend, tagSend, &
       iaSend, izSend, jaSend, jzSend, &
       bufSendStart, bufSendLength, bufSendTotalLength)
!      WRITE (*,FMT='("Em init. Sou o ",I2.2,", na pos. ", I2.2)') myNum,17; CALL flush(6)

       DO i=1,nm2
          DO j=1,nm3
        	advmnt_g(ng)%dxtW(i,j)=advmnt_g(ng)%l_dxtW(1,i,j)
        	advmnt_g(ng)%dytW(i,j)=advmnt_g(ng)%l_dytW(1,i,j)
          END DO
       END DO
       !if(ng==ngrids) iupdate_dxy=1
   !endif

  END SUBROUTINE InitialFieldsUpdate

 !--------------------------------------------------------------------
  SUBROUTINE Advec3d_X(m1,m2,m3, ia,iz,ja,jz,q0,u,den0,den1,dt,dxx,dd0,qn,mynum)
  !-------------------------
  ! This subroutine calculates change in mixing ratio (Q0) during time
  !  step DT due to advection along a grid IDIM in length. Mixing ratios
  !  from host code (C) are loaded into Q0 array, which is updated to QN.
  !  Velocities (U) and fluxes (FLUX) are specified at cell FACES, having
  !  dimensions 0:IDIM. U, Q0, QN, DXX and FLUX indices defined here:
  !  Densities at beg, end time (DEN0, DEN1) defined in HOST CODE
  !
  ! I grid->   |  1  |  2  |  I-1  |   I  |..   ..|  IDIM  | <- host grid
  ! U-array-> u(0)  u(1)  u(2)   u(i-1)  u(i)           u(IDIM)
  ! C-array->  | C(1)| C(2)| C(I-1)| C(I) |..   ..| C(IDIM)| mixing ratio
  ! DXX-arry-> | Dx1 | Dx2 | DxI-1 | DxI  |..   ..| DxIDIM |
  ! Density->  | Dd1 | Dd2 | DdI-1 | DdI  |..   ..| DdIDIM |
  !                 Q0 defined along 0 - IDIM+1 cells:
  !    |       | QN  | QN  |  QN   |  QN  |       |   QN   |        |
  !    |   Q0--|-Q0--|-Q0--|--Q0 --|--Q0--|..   ..|-- Q0 --|--Q0    |
  !    |    0  | 1   |  2  | I-1   |  I   |       |  IDIM  | IDIM+1 |
  !   lower BC |             <---   Q0 grid   --->         | upper BC
  !           Boundary conditions are stored in Q0 cells 0 & IDIM+1
  !
  !  Input to this subroutine, provided in common /sub/, and the calling
  !  arguments to this subroutine:
  !     IDIM - #of grid cells being updated
  !     Q0(0:IDIM+1)- Initial mixing ratio along 1-D array, with two
  !                 additional boundary value mixing ratios padded into the
  !                 0th and IDIM+1 cell locations
  !     U(0:IDIM)- velocities BETWEEN grid cells (at the "higher-I" edges of
  !                each grid cell in the array, units consistent with DX, DT
  !     DEN0(IDIM)- Initial fluid density, which needs to be updated during
  !                 multi-dimensional calculations, as noted in Calling code
  !     DEN1(IDIM)- Updated fluid density, which needs to be updated during
  !                 multi-dimensional calculations, as noted in calling code
  !     DT-         time step- units consistent with U
  !     DXX(IDIM)-  Grid cell length along advection direction, Units
  !                   consistent with DT and U
  !     DD0(0:IDIM)- Initial fluid density flowing BETWEEN each grid cell
  !                  (remains constant for all dimensions at the initial
  !                  fluid density of the 1st dimension of a 2-3 D calculation
  !               one can use UPSTREAM density here (DD0(I)= RHO0(I) if u>0
  !               or DD0(I)= RHO0(I+1) if u<0) where RHO0 is the initial
  !               fluid density at the beginning of the 1st dimensional
  !               advection step of a 2 or 3 D advection calculation done one
  !               step at a time
  !
  !  Output of this subroutine is an updated mixing ratio array QN(IDIM)
  !

  IMPLICIT none

  INTEGER,INTENT(IN)                        :: m1,m2,m3, ia,iz,ja,jz,mynum
  REAL   ,INTENT(in),DIMENSION(m1,m2,m3) :: q0
  REAL   ,INTENT(in),DIMENSION(m1,m2,m3)    :: u
  REAL   ,INTENT(in),DIMENSION(m1,m2,m3)    :: den0
  REAL   ,INTENT(in),DIMENSION(m1,m2,m3)    :: den1
  REAL   ,INTENT(in)                        :: dt
  REAL   ,INTENT(in),DIMENSION(m2,m3)       :: dxx
  REAL   ,INTENT(in),DIMENSION(m1,m2,m3)    :: dd0
  REAL   ,INTENT(OUT),DIMENSION(m1,m2,m3)   :: qn

  REAL,DIMENSION(m1,m2,m3)    :: flux
  REAL,DIMENSION(m1,m2,m3)    :: vcmax
  REAL,DIMENSION(m1,m2,m3)    :: vcmin
  LOGICAL,DIMENSION(m1,m2,m3) :: imxmn
  REAL,PARAMETER                       :: zr0=0.0
  REAL,PARAMETER                       :: EPS=1.e-6

  INTEGER :: idime
  INTEGER :: i,j,k
  REAL    :: cf,cf1,ck1,ck2,x1,x1n

  integer, parameter :: debug = 0
  integer, parameter :: on = 1

  INTEGER :: ii,ji,ii0,ji0,ie,je,ie0,je0,ipos,iia,iiz,nvar
  integer :: nf

  idime = m2 ! x dir
  qn=q0
  imxmn=.false.
  ! imxmn(idime+1)=.false.

 DO j=ja,jz

    DO k=2,m1-1
       ! Update mixing ratios and limit Fluxes going UP where u>0
       !  First assume upstream flux at edge of domain
       IF(u(k,1,j)>=zr0) flux(k,1,j)= q0(k,1,j)*u(k,1,j)*dt*dd0(k,1,j)
    END DO
  END DO

 DO j=ja,jz

  ! Identify local max and min, specify mixing ratio limits at new time
  DO  i=2,idime-1 ! ia,iz-1 or 1,iz-1

     DO k=2,m1-1

!        imxmn(k,i,j)=q0(k,i,j)>=max(q0(k,i-1,j),q0(k,i+1,j)) .or. & !=true if local
!                     q0(k,i,j)<=min(q0(k,i-1,j),q0(k,i+1,j))        !       extrema
        imxmn(k,i,j)=q0(k,i,j)>=(max(q0(k,i-1,j),q0(k,i+1,j))-eps) .or. & !=true if local
                     q0(k,i,j)<=(min(q0(k,i-1,j),q0(k,i+1,j))+eps)        !       extrema

        ck1= q0(k,i,j)
        ck2= q0(k,i,j)
        if(u(k,i,j  )< zr0) ck1= q0(k,i+1,j)
        if(u(k,i-1,j)>=zr0) ck2= q0(k,i-1,j)
        vcmax(k,i,j)= max( q0(k,i,j), ck1, ck2 )                      ! Eq-7
        vcmin(k,i,j)= min( q0(k,i,j), ck1, ck2 )                      ! Eq-7
        !  VCMAX and VCMIN are the absolute physical limits to the
        !     mixing ratio at t+dt. If these limits are ever violated,
        !     non-monotonic (oscillatory) behavior in solution results
        !
      enddo
   enddo
 enddo

 DO j=ja,jz

  ! Identify local max and min, specify mixing ratio limits at new time
  DO  i=2,idime-1 ! ia,iz-1 or 1,iz-1

     DO k=2,m1-1

        IF(u(k,i,j)<zr0) CYCLE

        IF(u(k,i-1,j)<zr0) THEN
           flux(k,i,j)= q0(k,i,j)*u(k,i,j)*dt*dd0(k,i,j)    !  outflow-only cell
        ELSE                              !      use upstream
           x1= dt*u(k,i,j)/dxx(i,j)               ! Courant number
           x1n= (1.-x1)*(q0(k,i+1,j)-q0(k,i-1,j))/4.
           !
           ! First, estimate mixing ratio in outflowing fluid (Cf)
           cf= q0(k,i,j) + x1n                                       !Eq-4a
           !
           !   Check to see if there is a peak (min) upwind and/or
           !    downwind of cell face
           IF(imxmn(k,i-1,j)) cf= q0(k,i,j) +MAX(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
           IF(imxmn(k,i+1,j)) cf= q0(k,i,j) +       (1.75 -.45*x1)*x1n   !Eq-10a
           !        CF= Q0(k,i,j) + 5.*X1N   ! uncomment this line for "full sharp"
           !
           !   Limit Cf to be between mixing ratio on either side of edge
           !      where flux is being calculated
           cf1= MIN( MAX( cf, MIN(q0(k,i,j),q0(k,i+1,j))  ), MAX(q0(k,i,j),q0(k,i+1,j)) )
           !- for debug purposes only
	   !if(debug==on)        CF1= CF     ! This statement IGNORES monotonic limitations.
           !                       you should uncomment this line and run your
           !                       advection calculation wth constant initial
           !                       mixing ratios everywhere. If you have properly
           !                       implemented this subroutine constant mixing
           !                       ratios should be maintained
           ! DEBUG    DEBUG   DEBUG    DEBUG   DEBUG
           !
	   !
           !   Calculate mixing ratio at new time, but limit to physically
           !    reasonable values
           qn(k,i,j) = MAX(vcmin(k,i,j),MIN(vcmax(k,i,j),          &   !eq-3&8
                      (q0(k,i,j)*den0(k,i,j)-x1*cf1*dd0(k,i,j)+flux(k,i-1,j)/dxx(i,j))/den1(k,i,j) ))

	   !
	   ! for debug purposes only
           !if(debug==on)  QN(k,i,j)=(Q0(k,i,j)*DEN0(k,i,j)-X1*CF1*DD0(k,i,j)+FLUX(k,I-1,j)/DXX(I,J))/DEN1(k,i,j)
           !            This statement IGNORES monotonic limitations.
           !                       you should uncomment this line and run your
           !                       advection calculation wth constant initial
           !                       mixing ratios everywhere. If you have properly
	   !                       implemented this subroutine constant mixing
           !                       ratios should be maintiained
           ! DEBUG    DEBUG   DEBUG    DEBUG   DEBUG
	   !
           !   Re-calculate OUTFLOWING flux before moving on to next cell
           !    Flux = CF1*X1*DD0 but it must be adjusted if a monotonic limit
           !    is encountered.
           flux(k,i,j)= dxx(i,j)*(q0(k,i,j)*den0(k,i,j) - qn(k,i,j)*den1(k,i,j)) + flux(k,i-1,j)
           !
        END IF                                                  !Eq-9a
     END DO

  END DO
END DO
  !
  ! If periodic boundary conditions are assumed, it is necessary
  !   to recalculate the updated mixing ratio at cell 1 if there
  !   is inflow to that cell from the boundary between IDIM and 1
  !   Here these statements are commented out, but should be uncommented
  !   if this subroutine is needed for periodic boundary conditions,
  !   and then one of the calling arguements to the subroutine is IPERIOD
  !   which is set to "1" if you assume period boundary conditions
  !      IF(IPERIOD==1) THEN
  !        IF(U(IDIM-1)>=ZR0.AND.U(IDIM)>=ZR0)
  !     &  QN(1)=(Q0(1)*DEN0(1)-FLUX(1)/DXX(1)+FLUX(IDIM)/DXX(1))/DEN1(1)
  !      END IF
  !
  ! Update mixing ratios and limit Fluxes going DOWN where u<0
  !  The logic of this loop through the grid line is identical
  !  to the "DO 10" Loop above, only you start at the highest I
  !  edge and work backwards to I=1
  !

 DO j=ja,jz

    DO k=2,m1-1
       IF(u(k,idime-1,j)<zr0) flux(k,idime-1,j)= &
    	       q0(k,idime,j)*u(k,idime-1,j)*dt*dd0(k,idime-1,j)
    END DO
  END DO


  DO j=ja,jz

   DO i=idime-1,2,-1 !iz,ia,-1
     DO k=2,m1-1

        IF(u(k,i-1,j)>=zr0) THEN           ! Inflow-only cell

	   IF(u(k,i,j)<zr0) qn(k,i,j)=  MAX(  vcmin(k,i,j),   MIN(   vcmax(k,i,j),&
                  (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(i,j) + flux(k,i-1,j)/dxx(i,j))/den1(k,i,j) ))
       ELSE
              x1=  dt*ABS(u(k,i-1,j))/dxx(i,j)     ! Courant number
              x1n= (1.-x1)*(q0(k,i-1,j)-q0(k,i+1,j))/4.
              cf= q0(k,i,j) + x1n                                       !Eq-4b
              IF(imxmn(k,i+1,j)) cf= q0(k,i,j) +MAX(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
              IF(imxmn(k,i-1,j)) cf= q0(k,i,j) +   (1.75 -.45*x1)*x1n       !Eq-10a
              cf1= MIN( MAX( cf, MIN(q0(k,i,j),q0(k,i-1,j)) ), MAX(q0(k,i,j),q0(k,i-1,j)) )
              !
	      !- for debug purposes only
	      !if(debug==on)  CF1= CF	! This statement IGNORES monotonic limitations.
              !			   you should uncomment this line and run your
              !			   advection calculation wth constant initial
              !			   mixing ratios everywhere. If you have properly
              !			   implemented this subroutine constant mixing
              !			   ratios should be maintained
	      ! -------------------------------------

	      IF(u(k,i,j)>=zr0) cf1= q0(k,i,j)     ! outflow-only cell upstream
              qn(k,i,j)= MAX(  vcmin(k,i,j),  MIN(   vcmax(k,i,j), 	  &   !Eq-3&8
                       (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(i,j)-x1*cf1*dd0(k,i-1,j))/den1(k,i,j) ))
              !print*,'i q= ', i,qn(i),line,q0(i),q0(i-1),q0(i+1)

	      !- for debug purposes only
              !if(debug==on) QN(k,i,j)=(Q0(k,i,j)*DEN0(k,i,j)-FLUX(k,I,j)/DXX(I,J)-X1*CF1*DD0(k,I-1,j))/DEN1(k,i,j)
              !		This statement IGNORES monotonic limitations.
              !		you should uncomment this line and run your
              !		advection calculation wth constant initial
              !		mixing ratios everywhere. If you have properly
              !		implemented this subroutine constant mixing
	      !		ratios should be maintiained
              ! -------------------------------------

	      flux(k,i-1,j)=dxx(i,j)*(qn(k,i,j)*den1(k,i,j) - q0(k,i,j)*den0(k,i,j)) + flux(k,i,j)!Eq-9b
              !
        END IF
     END DO

  END DO

  !END DO
  !
  ! If periodic boundary conditions are assumed, it is necessary
  !   to recalculate the updated mixing ratio at cell IDIM if there
  !   is inflow to that cell from the boundary between IDIM and 1
  !   Here these statements are commented out, but should be uncommented
  !   if this subroutine is needed for periodic boundary conditions,
  !   and then one of the calling arguements to the subroutine is IPERIOD
  !   which is set to "1" if you assume period boundary conditions
  !      IF(IPERIOD==1) THEN
  !      IF(U(1).LT.ZR0.AND.U(IDIM).LT.ZR0)
  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  !     &  QN(IDIM)=(Q0(IDIM)*DEN0(IDIM)-FLUX(0)/DXX(IDIM)+FLUX(IDIM-1)/
  !     &                DXX(IDIM))/DEN1(IDIM)
  !      END IF
  !
ENDDO !- big loop y-z


END SUBROUTINE Advec3d_X

!---------------------------------------------------------------------------

SUBROUTINE Advec3d_Y(m1,m2,m3, ia,iz,ja,jz,q0,u,den0,den1,dt,dxx,dd0,qn,mynum)
!-------------------------
  ! This subroutine calculates change in mixing ratio (Q0) during time
  !  step DT due to advection along a grid IDIM in length. Mixing ratios
  !  from host code (C) are loaded into Q0 array, which is updated to QN.
  !  Velocities (U) and fluxes (FLUX) are specified at cell FACES, having
  !  dimensions 0:IDIM. U, Q0, QN, DXX and FLUX indices defined here:
  !  Densities at beg, end time (DEN0, DEN1) defined in HOST CODE
  !
  ! I grid->   |  1  |  2  |  I-1  |   I  |..   ..|  IDIM  | <- host grid
  ! U-array-> u(0)  u(1)  u(2)   u(i-1)  u(i)           u(IDIM)
  ! C-array->  | C(1)| C(2)| C(I-1)| C(I) |..   ..| C(IDIM)| mixing ratio
  ! DXX-arry-> | Dx1 | Dx2 | DxI-1 | DxI  |..   ..| DxIDIM |
  ! Density->  | Dd1 | Dd2 | DdI-1 | DdI  |..   ..| DdIDIM |
  !                 Q0 defined along 0 - IDIM+1 cells:
  !    |       | QN  | QN  |  QN   |  QN  |       |   QN   |        |
  !    |   Q0--|-Q0--|-Q0--|--Q0 --|--Q0--|..   ..|-- Q0 --|--Q0    |
  !    |    0  | 1   |  2  | I-1   |  I   |       |  IDIM  | IDIM+1 |
  !   lower BC |             <---   Q0 grid   --->         | upper BC
  !           Boundary conditions are stored in Q0 cells 0 & IDIM+1
  !
  !  Input to this subroutine, provided in common /sub/, and the calling
  !  arguments to this subroutine:
  !     IDIM - #of grid cells being updated
  !     Q0(0:IDIM+1)- Initial mixing ratio along 1-D array, with two
  !                 additional boundary value mixing ratios padded into the
  !                 0th and IDIM+1 cell locations
  !     U(0:IDIM)- velocities BETWEEN grid cells (at the "higher-I" edges of
  !                each grid cell in the array, units consistent with DX, DT
  !     DEN0(IDIM)- Initial fluid density, which needs to be updated during
  !                 multi-dimensional calculations, as noted in Calling code
  !     DEN1(IDIM)- Updated fluid density, which needs to be updated during
  !                 multi-dimensional calculations, as noted in calling code
  !     DT-         time step- units consistent with U
  !     DXX(IDIM)-  Grid cell length along advection direction, Units
  !                   consistent with DT and U
  !     DD0(0:IDIM)- Initial fluid density flowing BETWEEN each grid cell
  !                  (remains constant for all dimensions at the initial
  !                  fluid density of the 1st dimension of a 2-3 D calculation
  !               one can use UPSTREAM density here (DD0(I)= RHO0(I) if u>0
  !               or DD0(I)= RHO0(I+1) if u<0) where RHO0 is the initial
  !               fluid density at the beginning of the 1st dimensional
  !               advection step of a 2 or 3 D advection calculation done one
  !               step at a time
  !
  !  Output of this subroutine is an updated mixing ratio array QN(IDIM)
  !
  IMPLICIT none

  INTEGER,INTENT(IN)                   :: m1,m2,m3, ia,iz,ja,jz,mynum
  INTEGER                  :: idime
  REAL   ,INTENT(in),DIMENSION(m1,m2,m3) :: q0
  REAL   ,INTENT(in),DIMENSION(m1,m2,m3) :: u
  REAL   ,INTENT(in),DIMENSION(m1,m2,m3) :: den0
  REAL   ,INTENT(in),DIMENSION(m1,m2,m3) :: den1
  REAL   ,INTENT(in)                     :: dt
  REAL   ,INTENT(in),DIMENSION(m2,m3)    :: dxx
  REAL   ,INTENT(in),DIMENSION(m1,m2,m3) :: dd0
  REAL   ,INTENT(OUT),DIMENSION(m1,m2,m3):: qn

  REAL,DIMENSION(m1,m2,m3)               :: flux
  REAL,DIMENSION(m1,m2,m3)               :: vcmax
  REAL,DIMENSION(m1,m2,m3)               :: vcmin
  LOGICAL,DIMENSION(m1,m2,m3)            :: imxmn
  REAL,PARAMETER                         :: zr0=0.0
  REAL,PARAMETER                         :: eps=1.e-6

  INTEGER :: i,j,k
  REAL    :: cf,cf1,ck1,ck2,x1,x1n

  integer, parameter :: debug = 0
  integer, parameter :: on = 1
  !
  INTEGER :: ii,ji,ii0,ji0,ie,je,ie0,je0,ipos,iia,iiz,nvar

  idime = m3 ! y dir
  qn= q0
  imxmn=.false.
  ! imxmn(idime+1)=.false.

 DO i=ia,iz

    ! Update mixing ratios and limit Fluxes going UP where u>0
    !  First assume upstream flux at edge of domain
    DO k=2,m1-1
       IF(u(k,i,1)>=zr0) flux(k,i,1)= q0(k,i,1)*u(k,i,1)*dt*dd0(k,i,1)
    END DO

 END DO

 !- big loop y-z
 DO i=ia,iz

  ! Identify local max and min, specify mixing ratio limits at new time
  DO  j=2,idime-1 ! ja,jz
     DO k=2,m1-1
!        imxmn(k,i,j)=q0(k,i,j)>=max(q0(k,i,j-1),q0(k,i,j+1)) .or. & !=true if local
!                     q0(k,i,j)<=min(q0(k,i,j-1),q0(k,i,j+1))	    !	    extrema
        imxmn(k,i,j)=q0(k,i,j)>=(max(q0(k,i,j-1),q0(k,i,j+1))-eps) .or. & !=true if local
                     q0(k,i,j)<=(min(q0(k,i,j-1),q0(k,i,j+1))+eps)	    !	    extrema
        ck1= q0(k,i,j)
        ck2= q0(k,i,j)
        if(u(k,i,j  )< zr0) ck1= q0(k,i,j+1)
        if(u(k,i,j-1)>=zr0) ck2= q0(k,i,j-1)
        vcmax(k,i,j)= max( q0(k,i,j), ck1, ck2 )                      ! Eq-7
        vcmin(k,i,j)= min( q0(k,i,j), ck1, ck2 )                      ! Eq-7
        !  VCMAX and VCMIN are the absolute physical limits to the
        !	mixing ratio at t+dt. If these limits are ever violated,
        !	non-monotonic (oscillatory) behavior in solution results

     enddo
  enddo
 enddo


 DO i=ia,iz

  ! Identify local max and min, specify mixing ratio limits at new time
  DO  j=2,idime-1 ! ja,jz
     DO k=2,m1-1

        IF(u(k,i  ,j)<zr0) CYCLE

        IF(u(k,i,j-1)<zr0) THEN
           flux(k,i,j)= q0(k,i,j)*u(k,i,j)*dt*dd0(k,i,j)    !  outflow-only cell

        ELSE                              !      use upstream
           x1= dt*u(k,i,j)/dxx(i,j)               ! Courant number
           x1n= (1.-x1)*(q0(k,i,j+1)-q0(k,i,j-1))/4.
           !
           ! First, estimate mixing ratio in outflowing fluid (Cf)
           cf= q0(k,i,j) + x1n                                       !Eq-4a
           !
           !   Check to see if there is a peak (min) upwind and/or
           !    downwind of cell face
           IF(imxmn(k,i,j-1)) cf= q0(k,i,j) +MAX(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
           IF(imxmn(k,i,j+1)) cf= q0(k,i,j) +       (1.75 -.45*x1)*x1n   !Eq-10a
           !        CF= Q0(k,i,j) + 5.*X1N   ! uncomment this line for "full sharp"
           !
           !   Limit Cf to be between mixing ratio on either side of edge
           !      where flux is being calculated
           cf1= MIN( MAX( cf, MIN(q0(k,i,j),q0(k,i,j+1))  ), MAX(q0(k,i,j),q0(k,i,j+1)) )
           !
	   !- for debug purposes only
	   !if(debug==on)        CF1= CF     ! This statement IGNORES monotonic limitations.
           !                       you should uncomment this line and run your
           !                       advection calculation wth constant initial
           !                       mixing ratios everywhere. If you have properly
           !                       implemented this subroutine constant mixing
           !                       ratios should be maintained
           !--------------------------------------------------
           !
	   !
           !   Calculate mixing ratio at new time, but limit to physically
           !    reasonable values
           qn(k,i,j)= MAX(  vcmin(k,i,j),   MIN(   vcmax(k,i,j),          &   !eq-3&8
                    (q0(k,i,j)*den0(k,i,j)-x1*cf1*dd0(k,i,j)+flux(k,i,j-1)/dxx(i,j))/den1(k,i,j) ))

	   ! for debug purposes only
           !if(debug==on)  QN(k,i,j)=(Q0(k,i,j)*DEN0(k,i,j)-X1*CF1*DD0(k,i,j)+FLUX(k,i,J-1)/DXX(I,J))/DEN1(k,i,j)
           !            This statement IGNORES monotonic limitations.
           !                       you should uncomment this line and run your
           !                       advection calculation wth constant initial
           !                       mixing ratios everywhere. If you have properly
	   !                       implemented this subroutine constant mixing
           !                       ratios should be maintiained
           !--------------------------------------------------
           !
	   !
           !   Re-calculate OUTFLOWING flux before moving on to next cell
           !    Flux = CF1*X1*DD0 but it must be adjusted if a monotonic limit
           !    is encountered.
	   !
           flux(k,i,j)= dxx(i,j)*(q0(k,i,j)*den0(k,i,j) - qn(k,i,j)*den1(k,i,j)) + flux(k,i,j-1) !Eq-9a
           !
        END IF
     END DO
  END DO
  !
  ! If periodic boundary conditions are assumed, it is necessary
  !   to recalculate the updated mixing ratio at cell 1 if there
  !   is inflow to that cell from the boundary between IDIM and 1
  !   Here these statements are commented out, but should be uncommented
  !   if this subroutine is needed for periodic boundary conditions,
  !   and then one of the calling arguements to the subroutine is IPERIOD
  !   which is set to "1" if you assume period boundary conditions
  !      IF(IPERIOD==1) THEN
  !        IF(U(IDIM-1)>=ZR0.AND.U(IDIM)>=ZR0)
  !     &  QN(1)=(Q0(1)*DEN0(1)-FLUX(1)/DXX(1)+FLUX(IDIM)/DXX(1))/DEN1(1)
  !      END IF
  !
  ! Update mixing ratios and limit Fluxes going DOWN where u<0
  !  The logic of this loop through the grid line is identical
  !  to the "DO 10" Loop above, only you start at the highest I
  !  edge and work backwards to I=1
  !
  END DO

 DO i=ia,iz

    DO k=2,m1-1
       IF(u(k,i,idime-1)<zr0) flux(k,i,idime-1)=q0(k,i,idime)*u(k,i,idime-1)*dt*dd0(k,i,idime-1)
    END DO
  END DO

  DO i=ia,iz

  DO j=idime-1,2,-1 !jz,ja,-1
     DO k=2,m1-1
        IF(u(k,i,j-1)>=zr0) THEN           ! Inflow-only cell
	   IF(u(k,i,j)<zr0) qn(k,i,j)=  MAX(  vcmin(k,i,j),   MIN(   vcmax(k,i,j),&
                  (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(i,j) + flux(k,i,j-1)/dxx(i,j))/den1(k,i,j) ))
           ELSE
              x1=  dt*ABS(u(k,i,j-1))/dxx(i,j)     ! Courant number
              x1n= (1.-x1)*(q0(k,i,j-1)-q0(k,i,j+1))/4.
              cf= q0(k,i,j) + x1n                                       !Eq-4b
              IF(imxmn(k,i,j+1)) cf= q0(k,i,j) +MAX(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
              IF(imxmn(k,i,j-1)) cf= q0(k,i,j) +   (1.75 -.45*x1)*x1n       !Eq-10a
              cf1= MIN( MAX( cf, MIN(q0(k,i,j),q0(k,i,j-1)) ), MAX(q0(k,i,j),q0(k,i,j-1)) )
              !
	      !- for debug purposes only
	      !if(debug==on)  CF1= CF	! This statement IGNORES monotonic limitations.
              !			   you should uncomment this line and run your
              !			   advection calculation wth constant initial
              !			   mixing ratios everywhere. If you have properly
              !			   implemented this subroutine constant mixing
              !			   ratios should be maintained
              ! ---------------------------------------

	      IF(u(k,i,j)>=zr0) cf1= q0(k,i,j)     ! outflow-only cell upstream
              qn(k,i,j)= MAX(  vcmin(k,i,j),  MIN(   vcmax(k,i,j), 	  &   !Eq-3&8
                     (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(i,j)-x1*cf1*dd0(k,i,j-1))/den1(k,i,j) ))
              !
	      !- for debug purposes only
              !if(debug==on) QN(k,i,j)=(Q0(k,i,j)*DEN0(k,i,j)-FLUX(k,i,J)/DXX(I,J)-X1*CF1*DD0(k,I,j-1))/DEN1(k,i,j)
              !		This statement IGNORES monotonic limitations.
              !			   you should uncomment this line and run your
              !			   advection calculation wth constant initial
              !			   mixing ratios everywhere. If you have properly
              !			   implemented this subroutine constant mixing
	      !			   ratios should be maintiained
              ! ---------------------------------------
	      !
	      flux(k,i,j-1)=dxx(i,j)*(qn(k,i,j)*den1(k,i,j) - q0(k,i,j)*den0(k,i,j)) + flux(k,i,j)!Eq-9b
        END IF
     END DO
  END DO
  !
  ! If periodic boundary conditions are assumed, it is necessary
  !   to recalculate the updated mixing ratio at cell IDIM if there
  !   is inflow to that cell from the boundary between IDIM and 1
  !   Here these statements are commented out, but should be uncommented
  !   if this subroutine is needed for periodic boundary conditions,
  !   and then one of the calling arguements to the subroutine is IPERIOD
  !   which is set to "1" if you assume period boundary conditions
  !      IF(IPERIOD==1) THEN
  !      IF(U(1).LT.ZR0.AND.U(IDIM).LT.ZR0)
  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  !     &  QN(IDIM)=(Q0(IDIM)*DEN0(IDIM)-FLUX(0)/DXX(IDIM)+FLUX(IDIM-1)/
  !     &                DXX(IDIM))/DEN1(IDIM)
  !      END IF
  !
ENDDO !- big loop x-z


END SUBROUTINE Advec3d_Y

!---------------------------------------------------------------------------
SUBROUTINE Advec3d_Z(m1,m2,m3, ia,iz,ja,jz,q0,u,den0,den1,dt,dxx,dd0,qn,mynum)
!-------------------------
  ! This subroutine calculates change in mixing ratio (Q0) during time
  !  step DT due to advection along a grid IDIM in length. Mixing ratios
  !  from host code (C) are loaded into Q0 array, which is updated to QN.
  !  Velocities (U) and fluxes (FLUX) are specified at cell FACES, having
  !  dimensions 0:IDIM. U, Q0, QN, DXX and FLUX indices defined here:
  !  Densities at beg, end time (DEN0, DEN1) defined in HOST CODE
  !
  ! I grid->   |  1  |  2  |  I-1  |   I  |..   ..|  IDIM  | <- host grid
  ! U-array-> u(0)  u(1)  u(2)   u(i-1)  u(i)           u(IDIM)
  ! C-array->  | C(1)| C(2)| C(I-1)| C(I) |..   ..| C(IDIM)| mixing ratio
  ! DXX-arry-> | Dx1 | Dx2 | DxI-1 | DxI  |..   ..| DxIDIM |
  ! Density->  | Dd1 | Dd2 | DdI-1 | DdI  |..   ..| DdIDIM |
  !                 Q0 defined along 0 - IDIM+1 cells:
  !    |       | QN  | QN  |  QN   |  QN  |       |   QN   |        |
  !    |   Q0--|-Q0--|-Q0--|--Q0 --|--Q0--|..   ..|-- Q0 --|--Q0    |
  !    |    0  | 1   |  2  | I-1   |  I   |       |  IDIM  | IDIM+1 |
  !   lower BC |             <---   Q0 grid   --->         | upper BC
  !           Boundary conditions are stored in Q0 cells 0 & IDIM+1
  !
  !  Input to this subroutine, provided in common /sub/, and the calling
  !  arguments to this subroutine:
  !     IDIM - #of grid cells being updated
  !     Q0(0:IDIM+1)- Initial mixing ratio along 1-D array, with two
  !                 additional boundary value mixing ratios padded into the
  !                 0th and IDIM+1 cell locations
  !     U(0:IDIM)- velocities BETWEEN grid cells (at the "higher-I" edges of
  !                each grid cell in the array, units consistent with DX, DT
  !     DEN0(IDIM)- Initial fluid density, which needs to be updated during
  !                 multi-dimensional calculations, as noted in Calling code
  !     DEN1(IDIM)- Updated fluid density, which needs to be updated during
  !                 multi-dimensional calculations, as noted in calling code
  !     DT-         time step- units consistent with U
  !     DXX(IDIM)-  Grid cell length along advection direction, Units
  !                   consistent with DT and U
  !     DD0(0:IDIM)- Initial fluid density flowing BETWEEN each grid cell
  !                  (remains constant for all dimensions at the initial
  !                  fluid density of the 1st dimension of a 2-3 D calculation
  !               one can use UPSTREAM density here (DD0(I)= RHO0(I) if u>0
  !               or DD0(I)= RHO0(I+1) if u<0) where RHO0 is the initial
  !               fluid density at the beginning of the 1st dimensional
  !               advection step of a 2 or 3 D advection calculation done one
  !               step at a time
  !
  !  Output of this subroutine is an updated mixing ratio array QN(IDIM)
  !
  IMPLICIT none

  INTEGER,INTENT(IN)                   :: m1,m2,m3, ia,iz,ja,jz,mynum
  INTEGER                  :: idime
  REAL   ,INTENT(in),DIMENSION(m1,m2,m3) :: q0
  REAL   ,INTENT(in),DIMENSION(m1,m2,m3) :: u
  REAL   ,INTENT(in),DIMENSION(m1,m2,m3) :: den0
  REAL   ,INTENT(in),DIMENSION(m1,m2,m3) :: den1
  REAL   ,INTENT(in)                     :: dt
  REAL   ,INTENT(in),DIMENSION(m1)       :: dxx
  REAL   ,INTENT(in),DIMENSION(m1,m2,m3) :: dd0
  REAL   ,INTENT(OUT),DIMENSION(m1,m2,m3):: qn

  REAL,DIMENSION(m1,m2,m3)    :: flux
  REAL,DIMENSION(m1,m2,m3)    :: vcmax
  REAL,DIMENSION(m1,m2,m3)    :: vcmin
  LOGICAL,DIMENSION(m1,m2,m3) :: imxmn
  REAL,PARAMETER                         :: zr0=0.0
  REAL,PARAMETER                         :: eps=1.e-6

  INTEGER :: i,j,k
  REAL    :: cf,cf1,ck1,ck2,x1,x1n

  integer, parameter :: debug = 0
  integer, parameter :: on = 1
  !

  idime = m1 ! z dir
  qn = q0
  imxmn=.false.
  ! imxmn(idime+1)=.false.


  !- big loop y-x
  ! Identify local max and min, specify mixing ratio limits at new time
  DO j=ja,jz ; DO i=ia,iz
    DO  k=2,idime-1 !

     imxmn(k,i,j)=q0(k,i,j)>=(max(q0(k-1,i,j),q0(k+1,i,j))-eps) .or. & !=true if local
                  q0(k,i,j)<=(min(q0(k-1,i,j),q0(k+1,i,j))+eps)	    !	    extrema
     ck1= q0(k,i,j)
     ck2= q0(k,i,j)
     if(u(k  ,i,j)< zr0) ck1= q0(k+1,i,j)
     if(u(k-1,i,j)>=zr0) ck2= q0(k-1,i,j)
     vcmax(k,i,j)= max( q0(k,i,j), ck1, ck2 )                      ! Eq-7
     vcmin(k,i,j)= min( q0(k,i,j), ck1, ck2 )                      ! Eq-7
  ENDDO
  !  VCMAX and VCMIN are the absolute physical limits to the
  !     mixing ratio at t+dt. If these limits are ever violated,
  !     non-monotonic (oscillatory) behavior in solution results
  !
 enddo; enddo

  !- big loop y-x
  DO j=ja,jz ; DO i=ia,iz
  ! Update mixing ratios and limit Fluxes going UP where u>0
  !  First assume upstream flux at edge of domain
  IF(u(1,i,j)>=zr0) flux(1,i,j)= q0(1,i,j)*u(1,i,j)*dt*dd0(1,i,j)

  DO k=2,idime-1

  !print*,'part 1'
     IF(u(k,i  ,j)<zr0) CYCLE

     IF(u(k-1,i,j)<zr0) THEN
        flux(k,i,j)= q0(k,i,j)*u(k,i,j)*dt*dd0(k,i,j)    !  outflow-only cell

     ELSE                              !      use upstream
        x1= dt*u(k,i,j)/dxx(k)               ! Courant number
        x1n= (1.-x1)*(q0(k+1,i,j)-q0(k-1,i,j))/4.
        !
        ! First, estimate mixing ratio in outflowing fluid (Cf)
        cf= q0(k,i,j) + x1n                                       !Eq-4a
        !
        !   Check to see if there is a peak (min) upwind and/or
        !    downwind of cell face
        IF(imxmn(k-1,i,j)) cf= q0(k,i,j) +MAX(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
        IF(imxmn(k+1,i,j)) cf= q0(k,i,j) +       (1.75 -.45*x1)*x1n   !Eq-10a
        !        CF= Q0(k,i,j) + 5.*X1N   ! uncomment this line for "full sharp"
        !
        !   Limit Cf to be between mixing ratio on either side of edge
        !      where flux is being calculated
        cf1= MIN( MAX( cf, MIN(q0(k,i,j),q0(k+1,i,j))  ), MAX(q0(k,i,j),q0(k+1,i,j)) )
        !
	!- for debug purposes only
	!if(debug==on)        CF1= CF     ! This statement IGNORES monotonic limitations.
        !                       you should uncomment this line and run your
        !                       advection calculation wth constant initial
        !                       mixing ratios everywhere. If you have properly
        !                       implemented this subroutine constant mixing
        !                       ratios should be maintained
        ! --------------------------------

	!
        !   Calculate mixing ratio at new time, but limit to physically
        !    reasonable values
        qn(k,i,j)= MAX(  vcmin(k,i,j),   MIN(   vcmax(k,i,j),          &   !eq-3&8
                   (q0(k,i,j)*den0(k,i,j)-x1*cf1*dd0(k,i,j)+flux(k-1,i,j)/dxx(k))/den1(k,i,j) ))

	! for debug purposes only
        !if(debug==on)  QN(k,i,j)=(Q0(k,i,j)*DEN0(k,i,j)-X1*CF1*DD0(k,i,j)+FLUX(k-1,i,j)/dxx(k))/DEN1(k,i,j)
        !            This statement IGNORES monotonic limitations.
        !                       you should uncomment this line and run your
        !                       advection calculation wth constant initial
        !                       mixing ratios everywhere. If you have properly
	!                       implemented this subroutine constant mixing
        !                       ratios should be maintiained
        ! -------------------------------------------------------------------
	!
        !   Re-calculate OUTFLOWING flux before moving on to next cell
        !    Flux = CF1*X1*DD0 but it must be adjusted if a monotonic limit
        !    is encountered.
        flux(k,i,j)= dxx(k)*(q0(k,i,j)*den0(k,i,j) - qn(k,i,j)*den1(k,i,j)) + flux(k-1,i,j)
     END IF                                                  !Eq-9a
   END DO
  END DO
END DO
  !
  ! If periodic boundary conditions are assumed, it is necessary
  !   to recalculate the updated mixing ratio at cell 1 if there
  !   is inflow to that cell from the boundary between IDIM and 1
  !   Here these statements are commented out, but should be uncommented
  !   if this subroutine is needed for periodic boundary conditions,
  !   and then one of the calling arguements to the subroutine is IPERIOD
  !   which is set to "1" if you assume period boundary conditions
  !      IF(IPERIOD==1) THEN
  !        IF(U(IDIM-1)>=ZR0.AND.U(IDIM)>=ZR0)
  !     &  QN(1)=(Q0(1)*DEN0(1)-FLUX(1)/DXX(1)+FLUX(IDIM)/DXX(1))/DEN1(1)
  !      END IF
  !
  ! Update mixing ratios and limit Fluxes going DOWN where u<0
  !  The logic of this loop through the grid line is identical
  !  to the "DO 10" Loop above, only you start at the highest I
  !  edge and work backwards to I=1
  !


 !- big loop y-x
  DO j=ja,jz ; DO i=ia,iz

   IF(u(idime-1,i,j)<zr0) flux(idime-1,i,j)=q0(idime,i,j)*u(idime-1,i,j)*dt*dd0(idime-1,i,j)

   DO k=idime-1,2,-1

       IF(u(k-1,i,j)>=zr0) THEN           ! Inflow-only cell

	IF(u(k,i,j)<zr0) qn(k,i,j)=  MAX(  vcmin(k,i,j),   MIN(   vcmax(k,i,j),&
           (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(k) + flux(k-1,i,j)/dxx(k))/den1(k,i,j) ))
        ELSE
           x1=  dt*ABS(u(k-1,i,j))/dxx(k)     ! Courant number
           x1n= (1.-x1)*(q0(k-1,i,j)-q0(k+1,i,j))/4.
           cf= q0(k,i,j) + x1n                                       !Eq-4b
           IF(imxmn(k+1,i,j)) cf= q0(k,i,j) +MAX(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
           IF(imxmn(k-1,i,j)) cf= q0(k,i,j) +   (1.75 -.45*x1)*x1n       !Eq-10a
           cf1= MIN( MAX( cf, MIN(q0(k,i,j),q0(k-1,i,j)) ), MAX(q0(k,i,j),q0(k-1,i,j)) )
           !
	   !- for debug purposes only
	   !if(debug==on)  CF1= CF	! This statement IGNORES monotonic limitations.
           !			   you should uncomment this line and run your
           !			   advection calculation wth constant initial
           !			   mixing ratios everywhere. If you have properly
           !			   implemented this subroutine constant mixing
           !			   ratios should be maintained
           !--------------------------------------------

	   IF(u(k,i,j)>=zr0) cf1= q0(k,i,j)     ! outflow-only cell upstream
           qn(k,i,j) = MAX(  vcmin(k,i,j),  MIN(   vcmax(k,i,j), 	  &   !Eq-3&8
                    (q0(k,i,j)*den0(k,i,j)-flux(k,i,j)/dxx(k)-x1*cf1*dd0(k-1,i,j))/den1(k,i,j) ))

	   !- for debug purposes only
           !if(debug==on) QN(k,i,j)=(Q0(k,i,j)*DEN0(k,i,j)-FLUX(k,i,j)/dxx(k)-X1*CF1*DD0(k-1,I,j))/DEN1(k,i,j)
           !		This statement IGNORES monotonic limitations.
           !			   you should uncomment this line and run your
           !			   advection calculation wth constant initial
           !			   mixing ratios everywhere. If you have properly
           !			   implemented this subroutine constant mixing
	   !			   ratios should be maintiained
	   !-------------------------------------------

	   flux(k-1,i,j)=dxx(k)*(qn(k,i,j)*den1(k,i,j) - q0(k,i,j)*den0(k,i,j)) + flux(k,i,j)!Eq-9b
     END IF
   END DO
  !
  ! If periodic boundary conditions are assumed, it is necessary
  !   to recalculate the updated mixing ratio at cell IDIM if there
  !   is inflow to that cell from the boundary between IDIM and 1
  !   Here these statements are commented out, but should be uncommented
  !   if this subroutine is needed for periodic boundary conditions,
  !   and then one of the calling arguements to the subroutine is IPERIOD
  !   which is set to "1" if you assume period boundary conditions
  !      IF(IPERIOD==1) THEN
  !      IF(U(1).LT.ZR0.AND.U(IDIM).LT.ZR0)
  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  !     &  QN(IDIM)=(Q0(IDIM)*DEN0(IDIM)-FLUX(0)/DXX(IDIM)+FLUX(IDIM-1)/
  !     &                DXX(IDIM))/DEN1(IDIM)
  !      END IF
  !
  ENDDO;ENDDO !- big loop y-x

END SUBROUTINE Advec3d_Z
!-------------------------------------------------------------------------------------
!---------------------------------------------------------------------------

SUBROUTINE Advec3d_Z_sedim(m1,m2,m3, ia,iz,ja,jz,q0,u,dt,dzt,rtgt,qn,mynum)
  !-------------------------
  ! This subroutine calculates change in mixing ratio (Q0) during time
  !  step DT due to advection along a grid IDIM in length. Mixing ratios
  !  from host code (C) are loaded into Q0 array, which is updated to QN.
  !  Velocities (U) and fluxes (FLUX) are specified at cell FACES, having
  !  dimensions 0:IDIM. U, Q0, QN, DXX and FLUX indices defined here:
  !  Densities at beg, end time (DEN0, DEN1) defined in HOST CODE
  !
  ! I grid->   |  1  |  2  |  I-1  |   I  |..   ..|  IDIM  | <- host grid
  ! U-array-> u(0)  u(1)  u(2)   u(i-1)  u(i)           u(IDIM)
  ! C-array->  | C(1)| C(2)| C(I-1)| C(I) |..   ..| C(IDIM)| mixing ratio
  ! DXX-arry-> | Dx1 | Dx2 | DxI-1 | DxI  |..   ..| DxIDIM |
  ! Density->  | Dd1 | Dd2 | DdI-1 | DdI  |..   ..| DdIDIM |
  !                 Q0 defined along 0 - IDIM+1 cells:
  !    |       | QN  | QN  |  QN   |  QN  |       |   QN   |        |
  !    |   Q0--|-Q0--|-Q0--|--Q0 --|--Q0--|..   ..|-- Q0 --|--Q0    |
  !    |    0  | 1   |  2  | I-1   |  I   |       |  IDIM  | IDIM+1 |
  !   lower BC |             <---   Q0 grid   --->         | upper BC
  !           Boundary conditions are stored in Q0 cells 0 & IDIM+1
  !
  !  Input to this subroutine, provided in common /sub/, and the calling
  !  arguments to this subroutine:
  !     IDIM - #of grid cells being updated
  !     Q0(0:IDIM+1)- Initial mixing ratio along 1-D array, with two
  !                 additional boundary value mixing ratios padded into the
  !                 0th and IDIM+1 cell locations
  !     U(0:IDIM)- velocities BETWEEN grid cells (at the "higher-I" edges of
  !                each grid cell in the array, units consistent with DX, DT
  !     DEN0(IDIM)- Initial fluid density, which needs to be updated during
  !                 multi-dimensional calculations, as noted in Calling code
  !     DEN1(IDIM)- Updated fluid density, which needs to be updated during
  !                 multi-dimensional calculations, as noted in calling code
  !     DT-         time step- units consistent with U
  !     DXX(IDIM)-  Grid cell length along advection direction, Units
  !                   consistent with DT and U
  !     DD0(0:IDIM)- Initial fluid density flowing BETWEEN each grid cell
  !                  (remains constant for all dimensions at the initial
  !                  fluid density of the 1st dimension of a 2-3 D calculation
  !               one can use UPSTREAM density here (DD0(I)= RHO0(I) if u>0
  !               or DD0(I)= RHO0(I+1) if u<0) where RHO0 is the initial
  !               fluid density at the beginning of the 1st dimensional
  !               advection step of a 2 or 3 D advection calculation done one
  !               step at a time
  !
  !  Output of this subroutine is an updated mixing ratio array QN(IDIM)
  !
  IMPLICIT none

  INTEGER,INTENT(IN)                   :: m1,m2,m3, ia,iz,ja,jz,mynum
  INTEGER                  :: idime
  REAL   ,INTENT(in),DIMENSION(m1,m2,m3) :: q0
  REAL   ,INTENT(in),DIMENSION(m1,m2,m3) :: u
  REAL   ,INTENT(in)                     :: dt
  REAL   ,INTENT(in) ,DIMENSION(m1)      :: dzt
  REAL   ,INTENT(in) ,DIMENSION(m2,m3)   :: rtgt
  REAL   ,INTENT(OUT),DIMENSION(m1,m2,m3):: qn

  REAL,DIMENSION(m1,m2,m3)    :: flux
  REAL,DIMENSION(m1,m2,m3)    :: vcmax
  REAL,DIMENSION(m1,m2,m3)    :: vcmin
  LOGICAL,DIMENSION(m1,m2,m3) :: imxmn
  REAL,PARAMETER                         :: zr0=0.0
  REAL,PARAMETER                         :: eps=1.e-6

  INTEGER :: i,j,k
  REAL    :: cf,cf1,ck1,ck2,x1,x1n,rtgti

  integer, parameter :: debug = 0
  integer, parameter :: on = 1
  !

  idime = m1 ! z dir
  qn = q0
  imxmn=.false.
  ! imxmn(idime+1)=.false.


  !- big loop y-x
  ! Identify local max and min, specify mixing ratio limits at new time
  DO j=ja,jz ; DO i=ia,iz
    DO  k=2,idime-1 !

     imxmn(k,i,j)=q0(k,i,j)>=(max(q0(k-1,i,j),q0(k+1,i,j))-eps) .or. & !=true if local
                  q0(k,i,j)<=(min(q0(k-1,i,j),q0(k+1,i,j))+eps)	    !	    extrema
     ck1= q0(k,i,j)
     ck2= q0(k,i,j)
     if(-u(k  ,i,j)< zr0) ck1= q0(k+1,i,j)
     if(-u(k-1,i,j)>=zr0) ck2= q0(k-1,i,j)
     if(k==2) ck2= q0(k,i,j) !for sedim only
     vcmax(k,i,j)= max( q0(k,i,j), ck1, ck2 )                      ! Eq-7
     vcmin(k,i,j)= min( q0(k,i,j), ck1, ck2 )                      ! Eq-7
  ENDDO
  !  VCMAX and VCMIN are the absolute physical limits to the
  !     mixing ratio at t+dt. If these limits are ever violated,
  !     non-monotonic (oscillatory) behavior in solution results
  !
 enddo; enddo


 !- big loop y-x
  DO j=ja,jz ; DO i=ia,iz

   rtgti=1./rtgt(i,j)

   flux(idime-1,i,j)=q0(idime,i,j)*(-u(idime-1,i,j))*dt

   DO k=idime-1,2,-1

!srf       x1=  dt*ABS(u(k-1,i,j))/dxx(k)     ! Courant number
           x1=  dt*ABS(u(k-1,i,j))*dzt(k)*rtgti     ! Courant number

	   if(k==2) x1 = 0. ! no flux below sfc terrain,for sedim only

	   x1n= (1.-x1)*(q0(k-1,i,j)-q0(k+1,i,j))/4.
           cf= q0(k,i,j) + x1n                                       !Eq-4b
           IF(imxmn(k+1,i,j)) cf= q0(k,i,j) +MAX(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
           IF(imxmn(k-1,i,j)) cf= q0(k,i,j) +   (1.75 -.45*x1)*x1n       !Eq-10a
           cf1= MIN( MAX( cf, MIN(q0(k,i,j),q0(k-1,i,j)) ), MAX(q0(k,i,j),q0(k-1,i,j)) )
           !
	   !- for debug purposes only
	   !if(debug==on)  CF1= CF	! This statement IGNORES monotonic limitations.
           !			   you should uncomment this line and run your
           !			   advection calculation wth constant initial
           !			   mixing ratios everywhere. If you have properly
           !			   implemented this subroutine constant mixing
           !			   ratios should be maintained
           !--------------------------------------------
           if(k>2) then  !for sedim only
            qn(k,i,j) = MAX(  vcmin(k,i,j),  MIN(   vcmax(k,i,j), 	  &   !Eq-3&8
!srf                 (q0(k,i,j)-flux(k,i,j)/dxx(k)      -x1*cf1) ))
                     (q0(k,i,j)-flux(k,i,j)*dzt(k)*rtgti-x1*cf1) ))
           else
             qn(k,i,j) = (q0(k,i,j)-flux(k,i,j)*dzt(k)*rtgti-x1*cf1)
	   endif
	   !- for debug purposes only
           !if(debug==on) QN(k,i,j)=(Q0(k,i,j)*DEN0(k,i,j)-FLUX(k,i,j)/dxx(k)-X1*CF1*DD0(k-1,I,j))/DEN1(k,i,j)
           !		This statement IGNORES monotonic limitations.
           !			   you should uncomment this line and run your
           !			   advection calculation wth constant initial
           !			   mixing ratios everywhere. If you have properly
           !			   implemented this subroutine constant mixing
	   !			   ratios should be maintiained
	   !-------------------------------------------

!srf	   flux(k-1,i,j)=dxx(k)             *(qn(k,i,j) - q0(k,i,j)) + flux(k,i,j)!Eq-9b
	   flux(k-1,i,j)=(1./(dzt(k)*rtgti))*(qn(k,i,j) - q0(k,i,j)) + flux(k,i,j)!Eq-9b
    END DO
  !
  ! If periodic boundary conditions are assumed, it is necessary
  !   to recalculate the updated mixing ratio at cell IDIM if there
  !   is inflow to that cell from the boundary between IDIM and 1
  !   Here these statements are commented out, but should be uncommented
  !   if this subroutine is needed for periodic boundary conditions,
  !   and then one of the calling arguements to the subroutine is IPERIOD
  !   which is set to "1" if you assume period boundary conditions
  !      IF(IPERIOD==1) THEN
  !      IF(U(1).LT.ZR0.AND.U(IDIM).LT.ZR0)
  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  !     &  QN(IDIM)=(Q0(IDIM)*DEN0(IDIM)-FLUX(0)/DXX(IDIM)+FLUX(IDIM-1)/
  !     &                DXX(IDIM))/DEN1(IDIM)
  !      END IF
  !
  ENDDO;ENDDO !- big loop y-x

END SUBROUTINE Advec3d_Z_sedim
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
SUBROUTINE Advec3d_Z_sedim_upw(m1,m2,m3, ia,iz,ja,jz,u,dt,dzt,rtgt,qn,mynum)
  !
  IMPLICIT none

  INTEGER,INTENT(IN)                   :: m1,m2,m3, ia,iz,ja,jz,mynum
  INTEGER                  :: idime
  REAL   ,INTENT(in),DIMENSION(m1,m2,m3) :: u
  REAL   ,INTENT(in)                     :: dt
  REAL   ,INTENT(in) ,DIMENSION(m1)      :: dzt
  REAL   ,INTENT(in) ,DIMENSION(m2,m3)   :: rtgt
  REAL   ,INTENT(OUT),DIMENSION(m1,m2,m3):: qn

  !REAL,DIMENSION(m1,m2,m3)    :: flux
  !REAL,DIMENSION(m1,m2,m3)    :: vcmax
  !REAL,DIMENSION(m1,m2,m3)    :: vcmin

  INTEGER :: i,j,k
  REAL    :: cf,cf1,ck1,ck2,x1,x1n,rtgti

  idime = m1 ! z dir

  !- big loop y-x
  DO j=ja,jz ; DO i=ia,iz

     rtgti=1./rtgt(i,j)
!srf dxx = dz = rtgti/dzt
!srf qn(idime-1,i,j) = qn(idime-1,i,j) / (1.0 - dt*u(idime-1,i,j)/dxx(idime-1)      )
     qn(idime-1,i,j) = qn(idime-1,i,j) / (1.0 + dt*u(idime-1,i,j)*dzt(idime-1)*rtgti)

     DO k=idime-2,2,-1 !

!srf    qn(k,i,j)= 1.0/(1.0+dt*u(k,i,j)/dxx(k))&
!srf               *( qn(k,i,j)+ dt*u(k,i,j) /dxx(k+1) * qn(k+1,i,j) )
        qn(k,i,j)= 1.0/(1.0 + dt*u(k,i,j)*dzt(k)*rtgti)&
                   *( qn(k,i,j) + dt*u(k+1,i,j)*dzt(k+1)*rtgti * qn(k+1,i,j) )

     !   tc(i,j,l,k) = 1.0/(1.0+dt_settl(k)*vd_cor/delz(i,j,l2))&
     !  	 *(tc(i,j,l,k) + dt_settl(k)*vd_cor /delz(i,j,l2-1) &
     !  	 * tc(i,j,l+1,k))

     ENDDO

  ENDDO;ENDDO !- big loop y-x

!  !- big loop y-x
!  ! Identify local max and min, specify mixing ratio limits at new time
!  DO j=ja,jz ; DO i=ia,iz
!    DO  k=2,idime-1 !
!
!     imxmn(k,i,j)=q0(k,i,j)>=(max(q0(k-1,i,j),q0(k+1,i,j))-eps) .or. & !=true if local
!                  q0(k,i,j)<=(min(q0(k-1,i,j),q0(k+1,i,j))+eps)	    !	    extrema
!     ck1= q0(k,i,j)
!     ck2= q0(k,i,j)
!     if(u(k  ,i,j)< zr0) ck1= q0(k+1,i,j)
!     if(u(k-1,i,j)>=zr0) ck2= q0(k-1,i,j)
!     vcmax(k,i,j)= max( q0(k,i,j), ck1, ck2 )                      ! Eq-7
!     vcmin(k,i,j)= min( q0(k,i,j), ck1, ck2 )                      ! Eq-7
!  ENDDO
!  VCMAX and VCMIN are the absolute physical limits to the
!     mixing ratio at t+dt. If these limits are ever violated,
!     non-monotonic (oscillatory) behavior in solution results
!
! enddo; enddo


!- big loop y-x
!  DO j=ja,jz ; DO i=ia,iz
!
!   flux(idime-1,i,j)=q0(idime,i,j)*u(idime-1,i,j)*dt
!
!   DO k=idime-1,2,-1
!
!           x1=  dt*ABS(u(k-1,i,j))/dxx(k)     ! Courant number
!           x1n= (1.-x1)*(q0(k-1,i,j)-q0(k+1,i,j))/4.
!           cf= q0(k,i,j) + x1n                                       !Eq-4b
!           IF(imxmn(k+1,i,j)) cf= q0(k,i,j) +MAX(1.5,1.2  +.6 *x1)*x1n   !Eq-10b
!           IF(imxmn(k-1,i,j)) cf= q0(k,i,j) +   (1.75 -.45*x1)*x1n       !Eq-10a
!           cf1= MIN( MAX( cf, MIN(q0(k,i,j),q0(k-1,i,j)) ), MAX(q0(k,i,j),q0(k-1,i,j)) )
!           !
!
!           qn(k,i,j) = MAX(  vcmin(k,i,j),  MIN(   vcmax(k,i,j), 	  &   !Eq-3&8
!                    (q0(k,i,j)-flux(k,i,j)/dxx(k)-x1*cf1) ))
!
!
!	   flux(k-1,i,j)=dxx(k)*(qn(k,i,j) - q0(k,i,j)) + flux(k,i,j)!Eq-9b
!    END DO
!  ENDDO;ENDDO !- big loop y-x

END SUBROUTINE Advec3d_Z_sedim_upw
!-------------------------------------------------------------------------------------

   subroutine CheckBorders(m1, m2, m3, field, &
	nRecv, procRecv, tagRecv, iaRecv, izRecv, jaRecv, jzRecv, &
	bufRecvStart, bufRecvLength, bufRecvTotalLength, &
	nSend, procSend, tagSend, iaSend, izSend, jaSend, jzSend, &
	bufSendStart, bufSendLength, bufSendTotalLength,mynum,op,ie)
     integer, intent(in) :: m1,mynum,op,ie
     integer, intent(in) :: m2
     integer, intent(in) :: m3
     real,    intent(inout) :: field(m1,m2,m3)
     integer, intent(in) :: nRecv
     integer, intent(in) :: procRecv(nRecv)
     integer, intent(in) :: tagRecv(nRecv)
     integer, intent(in) :: iaRecv(nRecv)
     integer, intent(in) :: izRecv(nRecv)
     integer, intent(in) :: jaRecv(nRecv)
     integer, intent(in) :: jzRecv(nRecv)
     integer, intent(in) :: bufRecvStart(nRecv)
     integer, intent(in) :: bufRecvLength(nRecv)
     integer, intent(in) :: bufRecvTotalLength
     integer, intent(in) :: nSend
     integer, intent(in) :: procSend(nSend)
     integer, intent(in) :: tagSend(nSend)
     integer, intent(in) :: iaSend(nSend)
     integer, intent(in) :: izSend(nSend)
     integer, intent(in) :: jaSend(nSend)
     integer, intent(in) :: jzSend(nSend)
     integer, intent(in) :: bufSendStart(nSend)
     integer, intent(in) :: bufSendLength(nSend)
     integer, intent(in) :: bufSendTotalLength

     integer :: fout,i
     CHARACTER :: opc

     fout=80+mynum
     opc='Y'
     IF(op==1) opc='X'

     IF(ie==0) THEN
        WRITE (fout,'("Borders updated, direction ",A)') opc
        RETURN
     END IF


     WRITE (fout,'(" Updating borders, direction  ",A)') opc ; CALL flush(fout)
     WRITE (fout,'("nRecv: ",I3.3," nSend: ",I3.3)') nRecv,nSend; CALL flush(fout)
     WRITE (fout,'("TotRecv: ",I6.6," TotSend: ",I6.6)') bufRecvTotalLength,bufSendTotalLength; CALL flush(fout)
     WRITE (fout,'(A)') '---------------------------------- Send ---------------------------'; CALL flush(fout)
     WRITE (fout,'(7(A3,1X),2(A,1X))') 'nSn','prc','tag','ia ','iz ','ja ','jz ','Start','Length'; CALL flush(fout)
     DO i=1,nRecv
        WRITE (fout,'(7(I3.3,1X),2(I6.6,1X))') i,procRecv(i),tagRecv(i),iaRecv(i),izRecv(i),jaRecv(i),jzRecv(i), &
	                  bufRecvStart(i),bufRecvLength(i); CALL flush(fout)
     END DO
     WRITE (fout,'(A)') '------------------------------- Receive  --------------------------'; CALL flush(fout)
     WRITE (fout,'(7(A3,1X),2(A,1X))') 'nRv','prc','tag','ia ','iz ','ja ','jz ','Start','Length'
     DO i=1,nRecv
        WRITE (fout,'(7(I3.3,1X),2(I6.6,1X))')	i,procSend(i),tagSend(i),iaSend(i),izSend(i),jaSend(i),jzSend(i), &
	                  bufSendStart(i),bufSendLength(i); CALL flush(fout)
     END DO

   end subroutine CheckBorders

  subroutine StoreNamelistFileAtradvc_mnt(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    advmnt = oneNamelistFile%advmnt
    GhostZoneLength=oneNamelistFile%GhostZoneLength

  end subroutine StoreNamelistFileAtRadvc_mnt


 !----------------------------------------------------

END MODULE monotonic_adv
