
MODULE module_chemistry_driver

    USE rconstants , ONLY: &
         cp,               & ! (IN)
         cpor,             & ! (IN)
         pi180,            & ! (IN)
         p00

    USE mem_grid, ONLY: &
         time,          & ! (IN)
         iyear1,        & ! (IN)
         imonth1,       & ! (IN)
         idate1,        & ! (IN)
         itime1,        & ! (IN)
         dtlt,          & ! (IN)
         ngrid,         & ! (IN)
         ngrids,        & ! (IN)
         dtlongn,       & ! (IN)
         dzt,           & ! (IN)
         zm,            & ! (IN)
         initial,       & ! (IN)
         zt,            & ! (IN)
         centlat,       & ! (IN)
         centlon,       & ! (IN)
         grid_g,        & ! %glat(IN), %glon(IN), %rtgt(IN), %lpw(IN), %dxt(IN)
         platn,         & ! (IN)
         plonn,         & ! (IN)
         deltaxn,       & ! (IN)
         deltayn,       & ! (IN)
         xt,            & ! (IN)
         yt               ! (IN)

    USE chem1_list, ONLY :          &
         chem_nspecies =>nspecies,  & ! (IN)
         spc_chem_alloc=>spc_alloc, & ! (IN)
         spc_chem_name =>spc_name,  & ! (IN)
         src,                       & ! (IN)
         on,                        & ! (IN)
         off,                       & ! (IN)
         transport,                 & ! (IN)
         nspecies,                  & ! (IN)
         spc_alloc,                 & ! (IN)
         no,                        & ! (IN)
         weight,                    & ! (IN)
         co,                        & ! (IN)
         photojmethod,              & ! (IN)
         nr,                        & ! (IN)
         nr_photo                     ! (IN)

    USE mem_chem1, ONLY:               &
         CHEMISTRY,                    & ! (IN)
         maxnspecies,                  & ! (IN)
         nspecies_chem_transported,    & ! (IN)
         nspecies_chem_no_transported, & ! (IN)
         transp_chem_index,            & ! (IN)
         no_transp_chem_index,         & ! (IN)
         nsrc,                         & ! (IN)
         src_name,                     & ! (IN)
         antro,                        & ! (IN)
         bburn,                        & ! (IN)
         bioge,                        & ! (IN)
         geoge,                        & ! (IN)
         chem1_src_z_dim_g,            & ! (IN)
         chem1_g,                      & ! %sc_t(INOUT), %sc_p(INOUT)
         chem1_src_g,                  & ! %sc_src(INOUT)
         chem1_src_vars,               & ! Type
         chem1_vars,                   & ! Type
         N_DYN_CHEM_N,                 & ! (IN)
         N_DYN_CHEM,                   & ! (INOUT)
         split_method,                 & ! (IN)
         isplit                          ! (IN)


    USE mem_chem1aq, ONLY: &
         CHEMISTRY_AQ        ! (IN)

    USE mem_basic, ONLY: &
         basic_g,        &   ! %up(IN), %vp(IN), %wp(IN), %theta(IN), %pp(IN), %pi0(IN), %dn0(IN?), %rv(IN), %rtp(IN)
         basic_vars          ! Type

    USE mem_radiate, ONLY : &
         radfrq,            & ! (IN)
         lonrad,            & ! (IN) (nao usado)
         radiate_g,         & ! %cosz(IN), %rlongup(IN), %albedt(IN)
         radiate_vars,      & ! Type
         ilwrtyp,           & ! (IN)
         iswrtyp              ! (IN)

    USE FastJX, ONLY :  &
         FastJX_driver, & ! Subroutine
	 Initialize_Fast_JX, &
	 fast_JX_initialized, &
         fast_JX_g        ! %jphoto(IN)

    USE node_mod, ONLY : &
         mynum,          & ! (IN)
!--(DMK-BRAMS-5.0-INI)----------------------------------------------------
	 nodemxp,        &
	 nodemyp,        &
	 nodemzp,        &
         nodeia,         & ! (IN)
         nodeiz,         & ! (IN)
         nodeja,         & ! (IN)
         nodejz,         & ! (IN)
!--(DMK-BRAMS-5.0-OLD)----------------------------------------------------
!         mmxp,           & ! (IN) substituido por nodemxp no BRAMS 5
!         mmyp,           & ! (IN) substituido por nodemyp no BRAMS 5
!         mmzp,           & ! (IN) substituido por nodemzp no BRAMS 5
!         mia,            & ! (IN) substituido por nodeia no BRAMS 5
!         miz,            & ! (IN) substituido por nodeiz no BRAMS 5
!         mja,            & ! (IN) substituido por nodeja no BRAMS 5
!         mjz,            & ! (IN) substituido por nodejz no BRAMS 5
!         mibcon,         & ! (IN) nao usado
!         mi0,            & ! (IN) nao usado
!         mj0,            & ! (IN) nao usado
!         maxgrds,        & ! (IN) substituido por ngrids no BRAMS 5
!--(DMK-BRAMS-5.0-FIM)----------------------------------------------------
         i0,             & ! (IN)
         j0,             & ! (IN)
         ibcon             ! (IN)

    USE mem_micro, ONLY: &
         micro_g,        & ! %rcp(IN)
         micro_vars        ! Type

    USE micphys, ONLY: &
         level           ! (IN)

    USE grid_dims, ONLY: &
         nzpmax,         &  ! (IN)
         nxpmax,         &  ! (IN)
         nypmax             ! (IN)

    USE mem_scratch1_grell, ONLY: &
         ierr4d                     ! (IN)

    USE mem_grell_param, ONLY: &
         mgmxp,                & ! (IN)
         mgmyp,                & ! (IN)
         maxiens,              & ! (IN)
         ngrids_cp               ! (IN)

    USE mem_cuparm, ONLY: &
         nnqparm            ! (IN)

    USE mem_stilt, ONLY: &
         iexev,          & ! (IN)
         stilt_g           ! %dnp(IN)

    USE mem_globrad, ONLY: &
         raddatfn            ! (IN)

    USE carma_fastjx, ONLY: &
         do3,               & ! (IN)
         daer,              & ! (IN)
         na                   ! (IN)
!----SRF 2015 not need anymore
!    USE aer1_list, ONLY :          &
!         aer_nspecies=>nspecies,   & ! (IN)
!         spc_aer_alloc=>spc_alloc, & ! (IN)
!         nmodes,                   & ! (IN)
!         spc_aer_name =>aer_name,  & ! (IN)
!         aer_bburn => bburn,       & ! (IN)
!         aer_sdust => sdust,       & ! (IN)
!         aer_urban => urban,       & ! (IN)
!         aer_bioge => bioge,       & ! (IN)
!         aer_marin => marin          ! (IN)
!----SRF 2015 not need anymore
    USE mem_aer1, ONLY:                   &
         aerosol, &
         aer_nvert_src=>aer1_src_z_dim_g, & ! (IN)
         aer1_vars                          ! Type

    USE spack_utils, ONLY: &
         allocindex,       & ! Subroutine
         index_g,          & ! IN
         nob,              & ! IN
         maxblock_size       ! IN

    USE mem_spack, ONLY: &
         alloc_spack,    & ! Subroutine
         spack_alloc,    & ! IN
         spack,          & ! OUT
         spack_2d,       & ! INOUT
         atol,           & ! IN
         rtol              ! IN

    USE solve_sparse, ONLY: &
         get_non_zeros

    USE mod_chem_spack_ros_dyndt, ONLY: &
         chem_ros_dyndt    ! Subroutine

    USE mod_chem_spack_ros, ONLY: &
         chem_ros          ! Subroutine

    USE mod_chem_spack_qssa, ONLY: &
         chem_qssa         ! Subroutine

    USE mod_chem_trans_gasaq, ONLY: &
         trans_gaq         ! Subroutine

    USE mod_chem_orage, ONLY: &
         eclair_driver     ! Subroutine

    USE mod_chem_spack_rodas3_dyndt, ONLY: &
         chem_rodas3_dyndt ! Subroutine

    USE uv_atten, ONLY: &
         initialize_uv_atten,  & ! Subroutine
         uv_attenuation,       & ! Subroutine
         uv_atten_g,           & ! (IN)
         uv_atten_initialized    ! (INOUT)

    USE mod_chem_trans_liq, ONLY: &
         trans_liq         ! Subroutine

    USE chem1aq_list, ONLY: &
         nspeciesaq,        & ! (IN)
         ind_gas              ! (IN)

    USE mem_chemic, ONLY: &
         chemic_g            ! (IN)

    USE mem_chem1aq, ONLY: &
         chem1aq_g           ! (INOUT)

    USE mem_aerad, ONLY: &
         nwave              ! (IN)

    USE mem_carma, ONLY: &
         carma              ! (IN)

    USE Extras, ONLY: &
         na_extra3d,  &     ! (IN)
         extra3d            ! (INOUT)

    USE FastJX, ONLY:  &
         fast_JX_g         ! %jphoto(IN)

    USE ModtuvDriver, ONLY: tuvDriver! (IN)

    USE mem_rrtm, ONLY: aot_rrtm_lw
    USE parrrtm, only : nbndlw


    IMPLICIT NONE

    PRIVATE


    PUBLIC ::  chemistry_driver ! Subroutine


CONTAINS

  !------------------------------------------------------------------
  SUBROUTINE chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,task,nt)

    INTEGER , INTENT(IN) :: mzp
    INTEGER , INTENT(IN) :: mxp
    INTEGER , INTENT(IN) :: myp
    INTEGER , INTENT(IN) :: ia
    INTEGER , INTENT(IN) :: iz
    INTEGER , INTENT(IN) :: ja
    INTEGER , INTENT(IN) :: jz
    INTEGER , INTENT(IN) :: task
    INTEGER , INTENT(IN) :: nt

    INTEGER, DIMENSION(ngrids) :: itchim
    INTEGER, DIMENSION(ngrids) :: ncycle

    INTEGER,PARAMETER :: block_size_qssa=45
    INTEGER,PARAMETER :: block_size_ros=64 ! use 64 for pgi / 48 for intel compilers

!     INTEGER,PARAMETER :: blockSize=1 !tmp



    !-  call NO production by "eclair"
    IF( task == 2 .AND. CHEMISTRY >= 0) THEN

       IF(ASSOCIATED(chem1_g(NO,ngrid)%sc_t)) THEN

          CALL eclair_driver(mzp,mxp,myp,ia,iz,ja,jz,dtlt,grid_g(ngrid)%rtgt, &
                             grid_g(ngrid)%dxt,grid_g(ngrid)%dyt,dzt,zt,nzpmax,mgmxp,mgmyp, &
                             maxiens,nnqparm(ngrid),basic_g(ngrid)%wp, &
!                             maxiens,ierr4d(:,:,:,ngrid),nnqparm(ngrid),basic_g(ngrid)%wp, &
                             basic_g(ngrid)%rtp,basic_g(ngrid)%rv,basic_g(ngrid)%pp, &
                             basic_g(ngrid)%pi0,basic_g(ngrid)%theta,chem1_g(NO,ngrid)%sc_t, &
                             weight,nspecies,no,cp,cpor)

       ENDIF

       RETURN

    ELSEIF( task == 3 .AND. CHEMISTRY >= 0) THEN

       IF (CHEMISTRY == 0) THEN

          CALL include_tracer_lifetime(mzp,mxp,myp,spc_alloc,nspecies,transport,on,co, &
                                       chem1_g(CO,ngrid)%sc_p(:,:,:), chem1_g(CO,ngrid)%sc_t(:))

       ELSE

          !-  call chemistry production/loss

          !- Determining number of blocks and masks, allocating spack type
          IF(.NOT. spack_alloc) THEN

!--(DMK-BRAMS-5.0-INI)----------------------------------------------------------------------
             CALL allocIndex(block_size_ros,nodemxp(mynum,:),nodemyp(mynum,:),nodemzp(mynum,:),&
	                     nodeia(mynum,:),nodeiz(mynum,:),nodeja(mynum,:),nodejz(mynum,:),&
!--(DMK-BRAMS-5.0-ORI)----------------------------------------------------------------------
!             CALL allocIndex(block_size_ros,mmxp,mmyp,mmzp,mia,miz,mja,mjz,mibcon,&
!--(DMK-BRAMS-5.0-FIM)----------------------------------------------------------------------
                             ngrids,dtlongn,N_DYN_CHEM)
             CALL alloc_spack(CHEMISTRY)
          END IF

          IF(uv_atten_initialized == 0) THEN

!--(DMK-BRAMS-5.0-INI)----------------------------------------------------------------------
              CALL initialize_uv_atten(ngrids,nodemxp(mynum,:),nodemyp(mynum,:))
!--(DMK-BRAMS-5.0-OLD)----------------------------------------------------------------------
!              CALL initialize_uv_atten(ngrids,mmxp,mmyp)
!--(DMK-BRAMS-5.0-FIM)----------------------------------------------------------------------

                uv_atten_initialized=1
          END IF

	  IF(Fast_JX_initialized == 0) THEN
	       CALL initialize_fast_JX(nodemxp(mynum,:),nodemyp(mynum,:),nodemzp(mynum,:),ngrids)
               Fast_JX_initialized = 1
	  END IF

          !- photolysis or attenuation calculations
          IF (MOD(time + .001,radfrq) .LT. dtlt .OR. time .LT. 0.001) THEN
             IF (TRIM(PhotojMethod) == 'FAST-JX') THEN

                CALL FastJX_driver(mzp,mxp,myp,ia,iz,ja,jz,nzpmax,ngrid,ngrids,grid_g(ngrid)%lpw, &
!--(DMK-BRAMS-5.0-INI)----------------------------------------------------------------------
                                   grid_g(ngrid)%rtgt,zm,imonth1,idate1,iyear1,time,dzt, &
				   nodemzp(mynum,:),nodemxp(mynum,:),nodemyp(mynum,:),i0,j0,ngrids,&
				   basic_g(ngrid)%pp,basic_g(ngrid)%pi0, &
!--(DMK-BRAMS-5.0-OLD)----------------------------------------------------------------------
!                                  grid_g(ngrid)%rtgt,zm,imonth1,idate1,iyear1,time,dzt,mmzp,mmxp, &
!           	  	  	   mmyp,i0,j0,maxgrds,basic_g(ngrid)%pp,basic_g(ngrid)%pi0, &
!--(DMK-BRAMS-5.0-FIM)----------------------------------------------------------------------
           		           basic_g(ngrid)%theta,basic_g(ngrid)%dn0,radiate_g(ngrid)%rlongup, &
           		  	   radiate_g(ngrid)%cosz,radiate_g(ngrid)%albedt,raddatfn,do3,daer,na)

             ELSEIF(TRIM(PhotojMethod) == 'LUT') THEN

                IF(ilwrtyp==4 .or. iswrtyp==4) THEN
                  CALL uv_attenuation(mxp,myp,ia,iz,ja,jz,i0,j0, &
                                    platn(ngrid),plonn(ngrid),deltaxn(ngrid),deltayn(ngrid), &
                                    nxpmax,nypmax,xt,yt,nwave,carma(ngrid)%aot, &
                                    ilwrtyp,iswrtyp,uv_atten_g(ngrid)%att,11)
                ELSEIF(ilwrtyp==6 .or. iswrtyp==6) THEN
                  CALL uv_attenuation(mxp,myp,ia,iz,ja,jz,i0,j0, &
                                    platn(ngrid),plonn(ngrid),deltaxn(ngrid),deltayn(ngrid), &
                                    nxpmax,nypmax,xt,yt,nbndlw,real(aot_rrtm_lw), &
                                    ilwrtyp,iswrtyp,uv_atten_g(ngrid)%att,2)
                END IF

	    ELSEIF(TRIM(PhotojMethod) == 'FAST-TUV') THEN

	       call tuvDriver(mzp,mxp,myp,ia,iz,ja,jz)!,nstr,wstart,wstop,nwint,blockSize,listFiles)

             ELSE
	        STOP ' unknonw photolysis calculation method '

             ENDIF
             IF(TRIM(PhotojMethod) == 'LUT') THEN
               call tuvDriver(mzp,mxp,myp,ia,iz,ja,jz)
          ENDIF
          ENDIF

          !-srf: including options for operator splitting
          !- current dynamic splitting
          N_DYN_CHEM=N_DYN_CHEM_N(ngrid)

          IF(split_method == 'PARALLEL' .AND. N_DYN_CHEM > 1) &
               CALL chem_accum(mzp,mxp,myp,ia,iz,ja,jz,1,nspecies_chem_transported,chem1_g(:,ngrid), &
                               transp_chem_index,N_DYN_CHEM)

          IF (MOD(time-(N_DYN_CHEM/ISPLIT-1)*dtlt , N_DYN_CHEM*dtlt) == 0. ) THEN
             !if(mynum==1)print*,'call chem: time=',time,ngrid; call flush(6)

             IF(CHEMISTRY == 1) THEN
                itchim(ngrid)= 75  ! teste 60
                !STOP 'for now use only CHEMISTRY = 4'
                !Determining number of blocks and masks
                !CALL AllocIndex(block_size_qssa,mmxp,mmyp,mmzp,mia,miz,mja,mjz,mibcon, &
                !	     ngrids,dtlongn,N_DYN_CHEM)

                CALL chem_qssa(mzp,mxp,myp,itchim(ngrid),nr,mynum,dtlt,basic_g(ngrid)%pp, &
           	 	       basic_g(ngrid)%pi0, basic_g(ngrid)%theta,basic_g(ngrid)%rv, &
           		       radiate_g(ngrid)%cosz,fast_JX_g(ngrid)%jphoto, &
           		       nspecies,nr_photo,weight,PhotoJMethod,chem1_g(:,ngrid), &
                               nob(ngrid),maxblock_size, &
                               index_g(ngrid)%block_end,index_g(ngrid)%indexk,index_g(ngrid)%indexi, &
                               index_g(ngrid)%indexj,index_g(ngrid)%kij_index, &
                               cp,cpor,p00, &
                               uv_atten_g(ngrid)%att)

                !Deallocating blocks data
                !CALL AllocIndex(block_size_qssa,mmxp,mmyp,mmzp,mia,miz,mja,mjz,mibcon, &
                !	     ngrids,dtlongn,N_DYN_CHEM)


             ELSEIF(CHEMISTRY == 2) THEN

                CALL chem_ros(mzp,mxp,myp,dtlt,basic_g(ngrid)%pp,basic_g(ngrid)%pi0, &
           	 	      basic_g(ngrid)%theta,basic_g(ngrid)%rv,basic_g(ngrid)%dn0, &
           		      radiate_g(ngrid)%cosz,micro_g(ngrid)%rcp,nspecies,nr,nr_photo,weight, &
           		      PhotojMethod,fast_JX_g(ngrid)%jphoto,maxnspecies,nspecies_chem_transported, &
           		      nspecies_chem_no_transported,transp_chem_index,no_transp_chem_index, &
                              chem1_g(:,ngrid),nob(ngrid),maxblock_size, &
                              index_g(ngrid)%block_end,index_g(ngrid)%indexk,index_g(ngrid)%indexi, &
                              index_g(ngrid)%indexj,index_g(ngrid)%kij_index, &
                              spack,spack_2d,get_non_zeros,&
			      N_DYN_CHEM,split_method,cp,cpor,p00, &
                              uv_atten_g(ngrid)%att)

             ELSEIF(CHEMISTRY == 3) THEN

                CALL chem_ros_dyndt(mzp,mxp,myp,dtlt,basic_g(ngrid)%pp,basic_g(ngrid)%pi0, &
                  		    basic_g(ngrid)%theta,basic_g(ngrid)%rv,basic_g(ngrid)%dn0, &
                     	 	    radiate_g(ngrid)%cosz,micro_g(ngrid)%rcp,nspecies,nr,nr_photo,weight, &
                     		    PhotojMethod,maxnspecies,nspecies_chem_transported, &
                     		    nspecies_chem_no_transported,transp_chem_index,no_transp_chem_index, &
                     		    chem1_g(:,ngrid),nob(ngrid),maxblock_size, &
                                    index_g(ngrid)%block_end,index_g(ngrid)%indexk,index_g(ngrid)%indexi, &
                                    index_g(ngrid)%indexj,index_g(ngrid)%kij_index,index_g(ngrid)%last_accepted_dt, &
                                    spack,spack_2d,atol,rtol,get_non_zeros,&
			            N_DYN_CHEM,split_method,cp,cpor,p00,fast_JX_g(ngrid)%jphoto, &
                                    uv_atten_g(ngrid)%att)
             ELSEIF(CHEMISTRY == 4) THEN

                CALL chem_rodas3_dyndt(mzp,mxp,myp,dtlt,basic_g(ngrid)%pp,basic_g(ngrid)%pi0, &
                     		       basic_g(ngrid)%theta,basic_g(ngrid)%rv,basic_g(ngrid)%dn0, &
                     		       radiate_g(ngrid)%cosz,micro_g(ngrid)%rcp,nspecies,nr,nr_photo,weight, &
                     		       PhotojMethod,maxnspecies,nspecies_chem_transported, &
                     		       nspecies_chem_no_transported,transp_chem_index,no_transp_chem_index, &
                     		       chem1_g(:,ngrid),nob(ngrid),maxblock_size, &
                                       index_g(ngrid)%block_end,index_g(ngrid)%indexk,index_g(ngrid)%indexi, &
                                       index_g(ngrid)%indexj,index_g(ngrid)%kij_index,index_g(ngrid)%last_accepted_dt, &
                                       spack,spack_2d,atol,rtol,get_non_zeros,&
				       N_DYN_CHEM,split_method,cp,cpor,p00,na_extra3d,extra3d(:,ngrid), &
                                       fast_JX_g(ngrid)%jphoto,uv_atten_g(ngrid)%att)

             ENDIF

             !- set dyn tend to zero, only for parallel splitting method
             IF(split_method == 'PARALLEL' .AND. N_DYN_CHEM > 1) &
                  CALL chem_accum(mzp,mxp,myp,ia,iz,ja,jz,2,nspecies_chem_transported,chem1_g(:,ngrid), &
                                  transp_chem_index,N_DYN_CHEM)

          ENDIF

       ENDIF
       ! change MP 11/02/08
    ELSEIF(task == 4  .AND. CHEMISTRY  >= 1 .AND. CHEMISTRY_AQ >= 1 ) THEN
       CALL trans_gaq(mzp,mxp,myp,ia,iz,ja,jz)

    ELSEIF(task == 5  .AND. CHEMISTRY  >= 1 .AND. CHEMISTRY_AQ >= 1 ) THEN
       CALL trans_liq(mzp,mxp,myp,ia,iz,ja,jz,chem1_g(:,ngrid),chem1aq_g(:,ngrid), &
                      nspeciesaq,ind_gas,chemic_g(ngrid)%coll,chemic_g(ngrid)%sedimr, &
                      micro_g(ngrid)%rcp,micro_g(ngrid)%rrp)
       ! change MP 11/02/08 - end

    ENDIF

  END SUBROUTINE chemistry_driver

!------------------------------------------------------------------------------------------

  SUBROUTINE include_tracer_lifetime(m1,m2,m3,spc_alloc_chem,nspecies,&
                                     transport,on,co,sc_p,sc_t)


    INTEGER , INTENT(IN)    :: m1
    INTEGER , INTENT(IN)    :: m2
    INTEGER , INTENT(IN)    :: m3
    INTEGER , INTENT(IN)    :: spc_alloc_chem(6,nspecies)
    INTEGER , INTENT(IN)    :: nspecies
    INTEGER , INTENT(IN)    :: transport
    INTEGER , INTENT(IN)    :: on
    INTEGER , INTENT(IN)    :: co
    REAL    , INTENT(IN)    :: sc_p(m1,m2,m3)
    REAL    , INTENT(INOUT) :: sc_t(m1*m2*m3)

    ! local declarations
    REAL,PARAMETER :: vm_CO_i   = 1./(30.*24.*3600.) ! 1/30 days

    ! - only for CO
    IF(spc_alloc_chem(transport,CO) == ON) THEN

       CALL accum_lifetime(m1*m2*m3     &
                           ,sc_t(1)     & ! tend
                           ,sc_p(1,1,1) & ! mix ratio
                           ,vm_CO_i   )   ! inv of lifetime (1/seconds)
    ENDIF

  END SUBROUTINE include_tracer_lifetime

!------------------------------------------------------------------------------------------

  SUBROUTINE accum_lifetime(m1,sc_t,sc_p,lf_inv)

    INTEGER , INTENT(IN)    :: m1
    REAL    , INTENT(INOUT) :: sc_t(m1)
    REAL    , INTENT(IN)    :: sc_p(m1)
    REAL    , INTENT(IN)    :: lf_inv

    INTEGER i

    DO i=1,m1
       sc_t(i) = sc_t(i) - lf_inv*sc_p(i)
    ENDDO
  END SUBROUTINE accum_lifetime

!------------------------------------------------------------------------------------------

  SUBROUTINE chem_accum(m1,m2,m3,ia,iz,ja,jz,ntask,nspecies_chem_transported,chem1_g, &
                        transp_chem_index,N_DYN_CHEM)

     INTEGER         , INTENT(IN) :: ntask
     INTEGER         , INTENT(IN) :: m1
     INTEGER         , INTENT(IN) :: m2
     INTEGER         , INTENT(IN) :: m3
     INTEGER         , INTENT(IN) :: ia
     INTEGER         , INTENT(IN) :: iz
     INTEGER         , INTENT(IN) :: ja
     INTEGER         , INTENT(IN) :: jz

    ! mem_chem1
    INTEGER          , INTENT(IN)    :: nspecies_chem_transported
    TYPE (chem1_vars), INTENT(INOUT) :: chem1_g(nspecies_chem_transported)
    INTEGER          , INTENT(INOUT) :: transp_chem_index(nspecies_chem_transported)
    INTEGER          , INTENT(IN)    :: N_DYN_CHEM


    !- local var
    INTEGER ispc,n,ixyz,ntps,i,j,k
    REAL N_DYN_CHEM_i

    ntps = m1 * m2 * m3

    IF(ntask == 1) THEN

       N_DYN_CHEM_i= REAL(1./N_DYN_CHEM)

       DO ispc=1,nspecies_chem_transported

          !- map the species to transported ones
          n=transp_chem_index(ispc)

          !- calculate the mean dynamic tendency for the entire chemistry timestep
 	  DO ixyz= 1,ntps
            chem1_g(n)%sc_t_dyn(ixyz) = chem1_g(n)%sc_t_dyn(ixyz) + &
                                              N_DYN_CHEM_i * chem1_g(n)%sc_t(ixyz)
            !call accum(ntps, chem1_g(n)%sc_t_dyn(1), chem1_g(n)%sc_t(1) )
         ENDDO

      ENDDO

      !	do j=ja,jz ; do i=ia,iz ;do k=2,m1
      !
      !
      !	!- Air temperature (K)
      !	extra3d(9,ngrid)%d3(k,i,j)=extra3d(9,ngrid)%d3(k,i,j)+N_DYN_CHEM_i * &
      !				   basic_g(ngrid)%theta(k,i,j)  	* &
      !				( basic_g(ngrid)%pp (k,i,j) + basic_g(ngrid)%pi0(k,i,j) )/cp
      !
      !	!- transform from Exner function to pressure
      !	extra3d(10,ngrid)%d3(k,i,j)=extra3d(10,ngrid)%d3(k,i,j)+N_DYN_CHEM_i * &
      !				   p00* ( ( basic_g(ngrid)%pp (k,i,j) + basic_g(ngrid)%pi0(k,i,j) )/cp )**cpor
      !
      !	!- Water vapor
      !	extra3d(11,ngrid)%d3(k,i,j)=extra3d(11,ngrid)%d3(k,i,j)+N_DYN_CHEM_i * &
      !				    basic_g(ngrid)% rv(k,i,j)
      !	enddo;enddo;enddo

   ELSE
      DO ispc=1,nspecies_chem_transported

         n=transp_chem_index(ispc) !- map the species to transported ones

         !- set to zero arrays for  accumulation over the next chem timestep
         chem1_g(n)%sc_t_dyn(1:ntps) = 0.

      ENDDO
      !	do j=ja,jz ; do i=ia,iz ;do k=2,m1
      !	extra3d(9 ,ngrid)%d3(k,i,j)=0.
      !	extra3d(10,ngrid)%d3(k,i,j)=0.
      !	extra3d(11,ngrid)%d3(k,i,j)=0.
      !	enddo;enddo;enddo

   ENDIF
 END  SUBROUTINE chem_accum
 !------------------------------------------------------------------


END MODULE module_chemistry_driver


SUBROUTINE latset_tracer(m1,m2,m3,ia,iz,ja,jz,ibcon,i0,j0,mynum,ap,uc,vc,dxu, &
     dxm,dyv,dym)
  USE mem_grid, ONLY: dtlt,ibnd,jdim,jbnd,nz
  USE mem_scratch, ONLY: vctr17,vctr18
  IMPLICIT NONE

  INTEGER :: m1,m2,m3,ia,iz,ja,jz,ibcon,i,j,k,lbw,lbe,lbs,lbn,i0,j0,mynum
  REAL :: thresh,dtlx,c1,dxr,dyr
  REAL, DIMENSION(m1,m2,m3) :: ap,uc,vc
  REAL, DIMENSION(m2,m3) :: dxu,dxm,dyv,dym
  !character(len=*) :: vnam

  IF (IAND(ibcon,1) .GT. 0) lbw = ia - 1
  IF (IAND(ibcon,2) .GT. 0) lbe = iz + 1
  IF (IAND(ibcon,4) .GT. 0) lbs = ja - 1
  IF (IAND(ibcon,8) .GT. 0) lbn = jz + 1

  thresh = 0.
  dtlx = dtlt


  IF (ibnd .NE. 4) THEN

     ! Western boundary for lsflg = 2  == constant inflow, radiative b.c. outflow
     ! Veja que o campo ap(k,i,j) na borda so e' modificado se houver outflow, isto e'
     ! em regioes de inflow o valor permanece constante e igual ao valor inicial (t=0).
     ! Para introduzir uma especie de "nudging" na condicao de contorno, modifique o
     ! vetor ap(k,lbw,j) para a 'western boundary', por exemplo. O mesmo para as outras
     ! bordas.

     IF (IAND(ibcon,1) .GT. 0) THEN
	DO j = 1,m3
	   dxr = dxu(ia,j) / dxu(lbw,j)
	   c1 = dtlx * dxu(lbw,j)
	   DO k = 1,m1
	      vctr17(k) = -c1 * uc(k,lbw,j)
	      vctr18(k) = ap(k,ia,j) + dxr * (ap(k,ia,j) - ap(k,ia+1,j))
	   ENDDO
	   DO k = 1,m1
	      IF (vctr17(k) .GE. thresh) ap(k,lbw,j) = vctr18(k)
!	      !---srf-tmp
!	      if(i+i0==12 .and. j+j0==1 .and. k==2) then
!	       print*,'1 ap=',ap(k,lbw,j),lbw,vctr17(k)
!	       call flush(6)
!	      endif
!	      !---srf-tmp
	   ENDDO
	ENDDO
     ENDIF

     !     Eastern Boundary for LSFLG =  2

     IF (IAND(ibcon,2) .GT. 0) THEN
	DO j = 1,m3
	   dxr = dxu(iz-1,j) / dxu(iz,j)
	   c1 = dtlx * dxu(iz,j)
	   DO k = 1,m1
	      vctr17(k) = c1 * uc(k,iz,j)
	      vctr18(k) = ap(k,iz,j) + dxr * (ap(k,iz,j) - ap(k,iz-1,j))
	   ENDDO
	   DO k = 1,m1
	      IF (vctr17(k) .GE. thresh)  ap(k,lbe,j) = vctr18(k)
!	      !---srf-tmp
!	      if(i+i0==12 .and. j+j0==1 .and. k==2) then
!	       print*,'2 ap=',ap(k,lbe,j),lbe,vctr17(k)
!	       call flush(6)
!	      endif
!	      !---srf-tmp
	   ENDDO
	ENDDO
     ENDIF
  ENDIF

  IF(jdim.EQ.1.AND.jbnd.NE.4)THEN

     !     Southern boundary for LSFLG  2

     IF (IAND(ibcon,4) .GT. 0) THEN
	DO i = 1,m2
	   dyr = dyv(i,ja) / dyv(i,lbs)
	   c1 = dtlx * dyv(i,lbs)
!srf - fix from rams 60
!             do k = 1,nz
              DO k = 1,m1
	      vctr17(k) = -c1 * vc(k,i,lbs)
	      vctr18(k) = ap(k,i,ja) + dyr * (ap(k,i,ja) - ap(k,i,ja+1))
	   ENDDO
	   DO k = 1,m1
	      IF (vctr17(k) .GE. thresh) ap(k,i,lbs) = vctr18(k)
!	      !---srf-tmp
!	      if(i+i0==12 .and. k==2) then
!	       print*,'4 ap=',ap(k,i,lbs),lbs,vctr17(k),mynum; call flush(6)
!	       write(mynum*10,100) lbs,i,i0,ja,mynum,ap(k,i,lbs),vctr17(k)*100000,ap(k,i,ja), ap(k,i,ja+1)
!	       100 format(1x,'4 ap=',5i4,5F10.2)
!	       call flush(mynum*10)
!	      endif
!	      !---srf-tmp
	   ENDDO
	ENDDO
     ENDIF

     !     Northern Boundary for LSFLG =  2

     IF (IAND(ibcon,8) .GT. 0) THEN
	DO i = 1,m2
	   dyr = dyv(i,jz-1) / dyv(i,jz)
	   c1 = dtlx * dyv(i,jz)
	   DO k = 1,m1
	      vctr17(k) = c1 * vc(k,i,jz)
	      vctr18(k) = ap(k,i,jz) + dyr * (ap(k,i,jz) - ap(k,i,jz-1))
	   ENDDO
	   DO k = 1,m1
	      IF (vctr17(k) .GE. thresh) ap(k,i,lbn) = vctr18(k)
	      !---srf-tmp
	      !if(i+i0==12  .and. k==2) then
	      ! print*,'8 ap=',ap(k,i,lbn),lbn,vctr17(k)
	      ! call flush(6)
	      !endif
	      !---srf-tmp
	   ENDDO
	ENDDO
     ENDIF
  ENDIF

END SUBROUTINE latset_tracer
!-----------------------------------------------------------------------


!------------------------------------------------------------------
  SUBROUTINE initial_condition(ng,m1,m2,m3)
    USE chem1_list, ONLY : NX=>nspecies,init_ajust, CO2!,CO,O3,ch4,no,no2!,form,no3,pan!,SO2
    USE mem_chem1, ONLY: chem1_g,chem1_src_g,CHEM_ASSIM
    USE mem_basic         , ONLY: basic_g
    USE aer1_list
    USE mem_aer1
    IMPLICIT NONE
    INTEGER :: ng,m1,m2,m3,i,j,k,ispc
    REAL,PARAMETER :: fcu =1.e+9 !=> ppbm
    REAL dummy,dum1,dum2,b
    integer ncount, icx1,jcx1,icx2,jcx2

    return
    chem1_g(CO2,ng)%sc_p(:,:,:)=390.* 44.0 / 28.97

!--(DMK-CCATT-INI)-----------------------------------------------------
    return    ! Corrige distorcao no CO e NO quando CHEM_ASSIM==0
!--(DMK-CCATT-OLD)-----------------------------------------------------
!    IF(CHEM_ASSIM > 0 ) RETURN
!--(DMK-CCATT-FIM)-----------------------------------------------------


    !-attention: values must be in ppbv (check if this correct, should
    !                                    not be in ppbm? SRF)
    !PRINT*,'============================= homogeneous initial. ====='
    !PRINT*,'CO= ',chem1_g(CO ,ng)%sc_p(1:m1,72,68)*fcu
    !PRINT*,'========================================================'

    !   aer1_g(accum,bburn,ng)%sc_p(10:15,1:m2,1:m3) = 300.
    !do k=11,m1
    !  chem1_g(CO,ng)%sc_p(k,:,:)=0.
    !enddo
    !chem1_g(O3,ng)%sc_p=chem1_g(CO,ng)%sc_p
    !----------------------------------  1 square   ----------------------------------
    !chem1_g(CO,ng)%sc_p= 0.
    !chem1_g(CO,ng)%sc_p(5:20,25:50,10:30) = 100.*28./28.96
    !chem1_g(NO,ng)%sc_p= 20.*28./28.96
    !chem1_g(NO,ng)%sc_p(5:20,25:50,10:30) = 100.*28./28.96

    !----------------------------------  1 square fim---------------------------
    !icx1=50
    !jcx1=50
    !icx2=50
    !jcx2=50
    !b=1./90.
    ! b=1./150.
    !b=1./450.
    !b=1./280.
    !do j=1,m3; do i=1,m2
    !  !gaussian
    !  !dum1 = (max( 0., exp(-b*( real(i-icx2)**2+ real(j-jcx2)**2 ) ) )*100)*28./28.96
    !  dum2 = (max( 0., exp(-b*( real(i-icx1)**2+ real(j-jcx1)**2 ) ) )*100)*28./28.96
    !  !
    !  !if(dum1>0.01) chem1_g(CO,ng)%sc_p(:,i,j)=dum1
    !  !if(dum2>0.01)
    !  chem1_g(CO,ng)%sc_p(:,i,j)=dum2
!----------------------------
    !enddo;enddo
    !do j=1,m3; do i=1,m2;do k=1,m1
    !
    !  if(chem1_g(CO,ng)%sc_p(k,i,j)<20.*28./28.96)  chem1_g(CO,ng)%sc_p(k,i,j)=20.*28./28.96
    !enddo;enddo;enddo
    !return
    !----------------------------------  squares   ----------------------------------
    ! chem1_g(CO,ng)%sc_p= 20.*28./28.96
    !- dois quadrados
    ! chem1_g(CO,ng)%sc_p(:,40:60,75:95) = 100.*28./28.96
    ! chem1_g(CO,ng)%sc_p(:,40:60,05:25) = 100.*28./28.96
     !- for linear correlations
    ! chem1_g(Ch4,ng)%sc_p(:,:,:)= 12.4 + 1.53*chem1_g(CO,ng)%sc_p(:,:,:)
    !- for quadr correlation
    ! chem1_g(nO,ng)%sc_p(:,:,:)= 12.4 + 0.01* 1.53*chem1_g(CO,ng)%sc_p(:,:,:)**2.
    !- for 4th correlation
    ! chem1_g(O3,ng)%sc_p(:,:,:)= 12.4 + 0.00025* 1.53*chem1_g(CO,ng)%sc_p(:,:,:)**4.
    !----------------------------------  squares   fim---------------------------
    !return
    !
    icx1=50
    jcx1=25
    icx2=50
    jcx2=75
    dummy= 100.*28./28.96
    !go to 333!<<<<<<<<<<<<<<<
    ! go to 555!<<<<<<<<<<<<<<<
    !----------------------------------
    !----------------------------------  gaussianas   ----------------------------------
    !b=1./90.
    ! b=1./150.
    !b=1./450.
    !b=1./280.
    do j=1,m3/2; do i=1,m2
      !gaussian
      !dum1 = (max( 0., exp(-b*( real(i-icx2)**2+ real(j-jcx2)**2 ) ) )*100)*28./28.96
      !dum2 = (max( 0., exp(-b*( real(i-icx1)**2+ real(j-jcx1)**2 ) ) )*100)*28./28.96
      !
      !if(dum1>0.01) chem1_g(CO,ng)%sc_p(:,i,j)=dum1
      !if(dum2>0.01)
      !chem1_g(CO,ng)%sc_p(:,i,j)=dum2
!--- for sum conservation
      !dum1 = (max( 0., exp(-0.5*b*( real(i-icx1-10)**2+ real(j-jcx1-3)**2 ) ) )*280)*28./28.96
      !chem1_g(form,ng)%sc_p(:,i,j)=dum1
!----------------------------

    enddo;enddo
     do j=m3/2,m3; do i=1,m2
      !gaussian
      !dum1 = (max( 0., exp(-b*( real(i-icx2)**2+ real(j-jcx2)**2 ) ) )*100)*28./28.96
   !dum2 = (max( 0., exp(-b*( real(i-icx1)**2+ real(j-jcx1)**2 ) ) )*100)*28./28.96
      !
      !chem1_g(CO,ng)%sc_p(:,i,j)=chem1_g(CO,ng)%sc_p(:,i,j-m3/2)
!--- for sum conservation
      !chem1_g(form,ng)%sc_p(:,i,j)=chem1_g(form,ng)%sc_p(:,i,j-m3/2)
!----------------------------
    enddo;enddo

    ! ass backgrd
    do j=1,m3; do i=1,m2;do k=1,m1

      !if(chem1_g(CO,ng)%sc_p(k,i,j)<20.*28./28.96)  chem1_g(CO,ng)%sc_p(k,i,j)=20.*28./28.96

!--- for sum conservation
      !if(chem1_g(form,ng)%sc_p(k,i,j)<30.*28./28.96) chem1_g(form,ng)%sc_p(k,i,j)=30.*28./28.96
!----------------------------

    enddo;enddo;enddo

    ! - two tracers
    !  chem1_g(CH4,ng)%sc_p(:,:,:)=120.*28./28.96 - chem1_g(CO,ng)%sc_p(:,:,:)
    ! 3  tracers
    !  chem1_g(NO2,ng)%sc_p(:,:,:)=360.*28./28.96 - chem1_g(form,ng)%sc_p(:,:,:)- chem1_g(CO,ng)%sc_p(:,:,:)
    ! 4  tracers
    ! chem1_g(PAN,ng)%sc_p(:,:,:)=40.*28./28.96 + 0.75*chem1_g(form,ng)%sc_p(:,:,:)+ 1.4* chem1_g(CO,ng)%sc_p(:,:,:)
    ! chem1_g(NO3,ng)%sc_p(:,:,:)=chem1_g(PAN,ng)%sc_p(:,:,:) + chem1_g(FORM,ng)%sc_p(:,:,:)+  chem1_g(CO,ng)%sc_p(:,:,:)
   !
    !- for linear correlations
    ! chem1_g(CH4,ng)%sc_p(:,:,:)= 22.4 + 1.53*chem1_g(CO,ng)%sc_p(:,:,:)
    !- for quadr correlation
    ! chem1_g(NO,ng)%sc_p(:,:,:)= 32.4 + 0.01* 1.53*chem1_g(CO,ng)%sc_p(:,:,:)**2.
    !- for 4th correlation
    ! chem1_g(O3,ng)%sc_p(:,:,:)= 12.4 + 0.5*0.000025* 1.53*chem1_g(CO,ng)%sc_p(:,:,:)**4.
    !----------------------------------  gaussianas fim  ----------------------------------

    return 

    !----------------------------------  triangulos      ----------------------------------
    333 continue
    ! dois triangulos
    do j=1,m3; do i=1,m2
      !gaussian
      dum1 = (max( 0., 1. - (1./15.)*sqrt( real(i-icx2)**2+ real(j-jcx2)**2 ))*100)*28./28.96
      dum2 = (max( 0., 1. - (1./15.)*sqrt( real(i-icx1)**2+ real(j-jcx1)**2 ))*100)*28./28.96
      !
      !if(dum1>0.01) chem1_g(CO,ng)%sc_p(:,i,j)=dum1
      !if(dum2>0.01) chem1_g(CO,ng)%sc_p(:,i,j)=dum2
    enddo;enddo
    ! ass backgrd
    do j=1,m3; do i=1,m2;do k=1,m1

      !if(chem1_g(CO,ng)%sc_p(k,i,j)<20.*28./28.96)  chem1_g(CO,ng)%sc_p(k,i,j)=20.*28./28.96
    enddo;enddo;enddo
    !----------------------------------  triangulos fim  ----------------------------------
    !- for linear correlations
     !chem1_g(Ch4,ng)%sc_p(:,:,:)= 12.4 + 1.53*chem1_g(CO,ng)%sc_p(:,:,:)
    !- for quadr correlation
     !chem1_g(nO,ng)%sc_p(:,:,:)= 12.4 + 0.01* 1.53*chem1_g(CO,ng)%sc_p(:,:,:)**2.
    !- for 4th correlation
     !chem1_g(O3,ng)%sc_p(:,:,:)= 12.4 + 0.00025* 1.53*chem1_g(CO,ng)%sc_p(:,:,:)**4.


    return
!    333 continue
    !--------------------------------------------------------------------
    ! dois circulos para slot
    do j=1,m3; do i=1,m2
      !if( sqrt( (real (i-icx1))**2+(real(j-jcx1))**2 ) .le. 15.) chem1_g(CO,ng)%sc_p(:,i,j)=dummy
      !if( sqrt( (real (i-icx2))**2+(real(j-jcx2))**2 ) .le. 15.) chem1_g(CO,ng)%sc_p(:,i,j)=dummy
      !chem1_g(CO,ng)%sc_p(5:20,10:30,10:30) = 100.*28./28.96
    enddo;enddo
    ! dois retangulos para o slot

    dummy= 20.*28./28.96
    do j=jcx1-2,jcx1+2; do i=icx1-2,icx1+15
      ! 1st retangle
      ! chem1_g(CO,ng)%sc_p(:,i,j)=dummy
    enddo;enddo
    do j=jcx2-2,jcx2+2; do i=icx2-15,icx2+2
      ! 1st retangle
      ! chem1_g(CO,ng)%sc_p(:,i,j)=dummy
    enddo;enddo

    !- for linear correlations
     !chem1_g(Ch4,ng)%sc_p(:,:,:)= 12.4 + 1.53*chem1_g(CO,ng)%sc_p(:,:,:)
    !- for quadr correlation
     !chem1_g(nO,ng)%sc_p(:,:,:)= 12.4 + 0.01* 1.53*chem1_g(CO,ng)%sc_p(:,:,:)**2.
    !- for 4th correlation
     !chem1_g(O3,ng)%sc_p(:,:,:)= 12.4 + 0.00025* 1.53*chem1_g(CO,ng)%sc_p(:,:,:)**4.

    !----------------------------------    slot  fim ----------------------

    555 continue






    return

    !DO ispc=O3,O3

    ! DO k=1,m1
    !  dummy=0.;ncount=0
    !
    ! DO j=1,m3; DO i=1,m2
    !   dummy= dummy+  chem1_g(ispc,ng)%sc_p(k,i,j)
    !	ncount=ncount+1
    ! enddo;enddo
     ! chem1_g(ispc,ng)%sc_p(k,:,:)    = dummy/ncount-10.
    !chem1_g(ispc,ng)%sc_p(k,:,:)=minval(chem1_g(ispc,ng)%sc_p(k,:,:))
    ! DO j=1,m3
    !  DO i=1,m2
    !   chem1_g(ispc,ng)%sc_p(k,i,j)=5.+minval(chem1_g(ispc,ng)%sc_p(k,2:m2-1,2:m3-1))
   !   if(i==2.and.j==2) print*,'o3=',o3,k,chem1_g(ispc,ng)%sc_p(k,i,j)
    ! DO j=7,13
    !  DO i=7,12
    !   DO k=3,20
  !       chem1_g(ispc,ng)%sc_p(k,i,j)=100.
    !ENDDO;ENDDO!;ENDDO;ENDDO

    !return

    DO ispc=1,nx
    !  chem1_g(ispc,ng)%sc_p(:,:,:)=init_ajust(ispc)*chem1_g(ispc ,ng)%sc_p(:,:,:)
    ENDDO

  !stop 33
    !RETURN

    !
    !do i=1,m2
    ! do j=1,m3
    ! chem1_g(O3 ,ng)%sc_p(1:13,i,j)=10. + i*j/100
    !  chem1_g(no ,ng)%sc_p(1:20,i,j)=20. + i*j**2/1000
    !  chem1_g(no2 ,ng)%sc_p(1:10,i,j)=1. + (i**2) * j*2/1000
    !enddo
    ! enddo
    !chem1_g(H2O,ng)%sc_p(:,:,:)=basic_g(ng)%RV(:,:,:)*fcu

    !RETURN

   ! chem1_g(CO,ng)%sc_p(:,:,:) = 60.
    !chem1_g(O3,ng)%sc_p(:,:,:) = 50.
    !chem1_g(O3,ng)%sc_p(25:,:,:) = 100.

    !chem1_g(H2O2,ng)%sc_p(1:13,:,:) = 9.

    !aer1_g(accum,bburn,ng)%sc_p(18:20,5:m2-5,5:m3-5) = 300.
    !aer1_g(accum,bburn,ng)%sc_p(1:3,3:m2-3,3:m3-3) = 300.

    return
    !aer1_g(coarse,bburn,ng)%sc_p(1:10,3:m2-3,3:m3-3) = 10.
    !aer1_g(nucle,urban,ng)%sc_p(1:10,3:m2-3,3:m3-3) = 10.
    !aer1_g(accum,urban,ng)%sc_p(1:10,3:m2-3,3:m3-3) = 10.

    !aer1_g(accum,bburn,ng)%sc_p(17:19,3:m2-3,3:m3-3) = 8.
    !aer1_g(coarse,bburn,ng)%sc_p(17:19,3:m2-3,3:m3-3) = 8.
    !aer1_g(nucle,urban,ng)%sc_p(17:19,3:m2-3,3:m3-3) = 8.
    !aer1_g(accum,urban,ng)%sc_p(17:19,3:m2-3,3:m3-3) = 8.


  END SUBROUTINE initial_condition


 subroutine aer_background(ngrid,m1,m2,m3,ia,iz,ja,jz)

    use aer1_list,       nspecies_aer   =>nspecies &
                        ,spc_alloc_aer  =>spc_alloc
    use mem_aer1 , only: aer1_g , AEROSOL
    USE mem_basic  ,  ONLY: basic_g

   IMPLICIT NONE
   INTEGER,INTENT(IN) :: m1,m2,m3,ia,iz,ja,jz,ngrid
   REAL,PARAMETER :: fcui=1.e-9     !de billion g [gas/part] /kg [ar] para kg/kg
   REAL, PARAMETER :: aer_pbl = 10. , aer_ft=5. ! microg/m^3
   integer  m1pbl,i,j,k,ispc,imode
   real aer_pbl_ppb,aer_ft_ppb

   !-- we are not longer using this approach
   return

   IF(AEROSOL == 0) return

!#ifdef SIMPLE
   !-
   IF(AEROSOL == 1) then
    ispc =urban
    imode=coarse
    if(spc_alloc_aer(transport,imode,ispc) /= ON) return

    m1pbl = min(10,m1)

    do j=ja,jz; do i=ia,iz
      ! pbl
      do k=1,m1pbl
    	aer_pbl_ppb=aer_pbl!*1.e-9 & ! from ug/m3 to kg/m3
    		           ! /basic_g(ngrid)%dn0(k,i,j) & ! to kg/kg
    		           ! *1.e9 ! to ppbv

  	!aer1_g(imode,ispc,ngrid)%sc_p(k,i,j)=max(aer_pbl_ppb,&
    	!		   aer1_g(imode,ispc,ngrid)%sc_p(k,i,j))
  	aer1_g(imode,ispc,ngrid)%sc_p(k,i,j)=aer_pbl_ppb

      enddo
      !free-troposphere
      do k=m1pbl+1,15 !m1
      !
       aer_ft_ppb=aer_ft!*1.e-9 & ! from ug/m3 to kg/m3
      		        !/basic_g(ngrid)%dn0(k,i,j) & ! to kg/kg
       		        !*1.e9 ! to ppbv
      !	aer1_g(imode,ispc,ngrid)%sc_p(k,i,j)=max(aer_ft_ppb,&
      !			   aer1_g(imode,ispc,ngrid)%sc_p(k,i,j))
      	aer1_g(imode,ispc,ngrid)%sc_p(k,i,j)=aer_ft_ppb
      !
      enddo
      do k=16, m1!m1
      !
      	aer1_g(imode,ispc,ngrid)%sc_p(k,i,j)=0.
      !
      enddo

    enddo;enddo
  ELSEIF(AEROSOL == 2) then
!#elif MATRIX
!---srf incluir instrucoes para o MATRIX

!#endif
  ENDIF
end subroutine aer_background
