MODULE Aero_Setup
!-------------------------------------------------------------------------------
!@sum     This module contains various aerosol microphysical variables and routines.
!@auth    Susanne Bauer/Doug Wright
!-------------------------------------------------------------------------------
      !USE aer1_list      , ONLY: nmodes
      USE memMatrix 

      IMPLICIT NONE

      logical :: firstime_indices = .true.
      
      CONTAINS


      SUBROUTINE Setup_Config
!-------------------------------------------------------------------------------
!     Routine to initialize variables that depend upon choice of mechanism.
!-------------------------------------------------------------------------------
      implicit none
      integer :: i,j,k,index,idim
      logical :: found
      !integer :: nmodes=16 !LFR teste
      !-----------------------------------------------------------------------
      ! get the number of modes in the selected mechanism.
      !-----------------------------------------------------------------------
      idim = 0        
      if ( mech .eq. 1 ) idim = nm1
      if ( mech .eq. 2 ) idim = nm2
      if ( mech .eq. 3 ) idim = nm3
      if ( mech .eq. 4 ) idim = nm4
      if ( mech .eq. 5 ) idim = nm5
      if ( mech .eq. 6 ) idim = nm6
      if ( mech .eq. 7 ) idim = nm7
      if ( mech .eq. 8 ) idim = nm8
      if ( idim .ne. nmodes ) then
        write(*,*)'error in mechanism number mech: mech = ', mech
        write(*,*)'nmodes=',nmodes,' idim=',idim
        stop
      endif
      !print *, 'NaeroBox: ',naerobox; CALL flush(6)
      allocate( mode_name(idim) )      
      allocate( mode_spcs(nmass_spcs,idim) )
      allocate( aero_spcs(naerobox) )      
      allocate( icond(idim) )      
      allocate( citable(idim,idim) )
      mode_name(:)   = '   '
      mode_spcs(:,:) = 0
      aero_spcs(:)   = '                '     
      icond(:)       = 0
      citable(:,:)   = '   '

      !-------------------------------------------------------------------------
      ! initialize arrays for coagulation interactions, mode names, mode
      ! species, condensation flags, and whether mode bc3 is present in the
      ! selected mechanism.
      !-------------------------------------------------------------------------
      if     ( mech .eq. 1 ) then
        citable(:,:) = citable1(:,:)
        mode_name(:) = mname(modes1(:))
        mode_spcs(:,:) = mspcs(:,modes1(:))
        icond(:) = icond1(:)
        include_bc3 = .true.

      elseif ( mech .eq. 2 ) then
        citable(:,:) = citable2(:,:)
        mode_name(:) = mname(modes2(:))
        mode_spcs(:,:) = mspcs(:,modes2(:))
        icond(:) = icond2(:)
        include_bc3 = .false.
      elseif ( mech .eq. 3 ) then
        citable(:,:) = citable3(:,:)
        mode_name(:) = mname(modes3(:))
        mode_spcs(:,:) = mspcs(:,modes3(:))
        icond(:) = icond3(:)
        include_bc3 = .false.
      elseif ( mech .eq. 4 ) then
        citable(:,:) = citable4(:,:)
        mode_name(:) = mname(modes4(:))
        mode_spcs(:,:) = mspcs(:,modes4(:))
        icond(:) = icond4(:)
        include_bc3 = .false.
      elseif ( mech .eq. 5 ) then
        citable(:,:) = citable5(:,:)
        mode_name(:) = mname(modes5(:))
        mode_spcs(:,:) = mspcs(:,modes5(:))
        icond(:) = icond5(:)
        include_bc3 = .true.  
      elseif ( mech .eq. 6 ) then
        citable(:,:) = citable6(:,:)
        mode_name(:) = mname(modes6(:))
        mode_spcs(:,:) = mspcs(:,modes6(:))
        icond(:) = icond6(:)
        include_bc3 = .false. 
      elseif ( mech .eq. 7 ) then
        citable(:,:) = citable7(:,:)
        mode_name(:) = mname(modes7(:))
        mode_spcs(:,:) = mspcs(:,modes7(:))
        icond(:) = icond7(:)
        include_bc3 = .false. 
      elseif ( mech .eq. 8 ) then
        citable(:,:) = citable8(:,:)
        mode_name(:) = mname(modes8(:))
        mode_spcs(:,:) = mspcs(:,modes8(:))
        icond(:) = icond8(:)
        include_bc3 = .false. 
      endif
      !-------------------------------------------------------------------------
      ! set intermodal_transfer according to setting in aero_param.f.
      ! if mode akk is not present in the current mechanism, set to .false.
      !-------------------------------------------------------------------------
      intermodal_transfer = set_intermodal_transfer
      found = .false.
      do i=1, nmodes
        if ( mode_name(i) .eq. 'AKK' ) found = .true.
      enddo

      if ( .not. found ) then
        if ( write_log ) write(aunit1,'(/A/)') &
          'INTERMODAL TRANSFER (AKK->ACC) TURNED OFF SINCE MODE AKK IS ABSENT.'
        intermodal_transfer = .false.
      ENDIF
      !-------------------------------------------------------------------------
      ! Check that all receptor modes in the CITABLE array are defined modes
      !   in the present mechanism.
      ! Also check that CITABLE is a symmetric matrix.
      !-------------------------------------------------------------------------
      do i=1, nmodes  
        do j=1, nmodes
          if( citable(i,j) .ne. citable(j,i) ) then
            write(*,*) 'citable(i,j) must be a symmetric matrix.'
            write(*,*) 'the citable(i,j) set in aero_config.f is asymmetric for i, j = ', i, j
            stop
          endif
          found = .false.
          do k=1, nmodes
            if( citable(i,j) .eq. mode_name(k) ) found = .true.
          enddo
          if( citable(i,j) .eq. 'off' ) then 
            found = .true.    ! i-j interaction has been turned off
            ! write(36,'(a,2i5,a5)')'i,j,citable(i,j) = ', i,j,citable(i,j)
          endif
          if( .not. found ) then
            WRITE(*,*)'INVALID RECEPTOR MODE NAME: I, J, CITABLE(I,J) = ', &
                                                  i, j, citable(i,j)
            STOP
          ENDIF
        ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      ! Setup the indices to the AERO array.
      !-------------------------------------------------------------------------
      aero_spcs(1) = 'MASS_NITRATE'
      aero_spcs(2) = 'MASS_AMMONIUM'
      aero_spcs(3) = 'MASS_WATER'
      index = 3          ! The first three values of INDEX (1, 2, 3) are already 
                         ! assigned to NO3, NH4, and H2O.

      IF ( WRITE_LOG ) THEN
        WRITE(aunit1,'(/2A/)') 'MODE #   MODE_NAME   CHEM_SPC #   CHEM_SPC_NAME   location in AERO', &
                              '         AERO_SPCS'
        WRITE(aunit1,90) 0,'NO3',0,'ANO3',1,'MASS_NO3        '
        WRITE(aunit1,90) 0,'NH4',0,'ANH4',2,'MASS_NH4        '
        WRITE(aunit1,90) 0,'H2O',0,'AH2O',3,'MASS_H2O        '
      ENDIF

      do i=1, nmodes
        do j=1, nmass_spcs
          if ( mode_spcs(j,i) .gt. 0 ) then  ! this mode contains species j.
            index = index + 1
            call Setup_Indices(i,j,index,0)
            aero_spcs(index) = 'MASS_'//MODE_NAME(I)//'_'//CHEM_SPC_NAME(J)
            IF ( write_log ) then
              write(aunit1,90) i,mode_name(i),j,chem_spc_name(j),index,aero_spcs(index)
              ! WRITE(31,'(8X,A23)') AERO_SPCS(INDEX)//' = OMIT'
            ENDIF
          ENDIF
        ENDDO

        do j=1, npoints
          index = index + 1
          call Setup_Indices(i,j,index,1)     ! set number conc. indices.
          if ( j .eq. 1 )  aero_spcs(index) = 'NUMB_'//mode_name(i)//'_1'   ! first  quadrature point
          if ( j .eq. 2 )  aero_spcs(index) = 'NUMB_'//mode_name(i)//'_2'   ! second quadrature point
          IF ( WRITE_LOG ) THEN
            write(aunit1,90) i,mode_name(i),j,'numb',index,aero_spcs(index)
            ! WRITE(31,'(8X,A23)') AERO_SPCS(INDEX)//' = OMIT'
          ENDIF
        ENDDO
      ENDDO

      if ( index .ne. naerobox ) then
        WRITE(*,*)'INDEX .NE. NAEROBOX', index, naerobox    ! size of the AERO and AERO_SPCS arrays
        STOP
      ENDIF

      !-------------------------------------------------------------------------------------------------------------------
      ! Setup maps etc. needed to convert sea salt mass concentrations to number concentrations
      !   using mean particle masses for the sea salt modes.
      !-------------------------------------------------------------------------------------------------------------------
      CALL Setup_Seasalt_Maps

!-------------------------------------------------------------------------------------------------------------------------
!     Indices of the AERO array.
!-------------------------------------------------------------------------------------------------------------------------
!     WRITE(*,*)       MASS_NO3,   MASS_NH4,   MASS_H2O
!     WRITE(*,*)       NUMB_AKK_1, NUMB_AKK_2, MASS_AKK_SULF,  
!    &                 NUMB_ACC_1, NUMB_ACC_2, MASS_ACC_SULF,
!    &                 NUMB_DD1_1, NUMB_DD1_2, MASS_DD1_SULF,                               MASS_DD1_DUST, 
!    &                 NUMB_DS1_1, NUMB_DS1_2, MASS_DS1_SULF,                               MASS_DS1_DUST, 
!    &                 NUMB_DD2_1, NUMB_DD2_2, MASS_DD2_SULF,                               MASS_DD2_DUST, 
!    &                 NUMB_DS2_1, NUMB_DS2_2, MASS_DS2_SULF,                               MASS_DS2_DUST, 
!    &                 NUMB_SSA_1, NUMB_SSA_2, MASS_SSA_SULF,                                              MASS_SSA_SEAS, 
!    &                 NUMB_SSC_1, NUMB_SSC_2, MASS_SSC_SULF,                                              MASS_SSC_SEAS,
!    &                 NUMB_SSS_1, NUMB_SSS_2, MASS_SSS_SULF,                                              MASS_SSS_SEAS,
!    &                 NUMB_OCC_1, NUMB_OCC_2, MASS_OCC_SULF,                MASS_OCC_OCAR,
!    &                 NUMB_BC1_1, NUMB_BC1_2, MASS_BC1_SULF, MASS_BC1_BCAR,
!    &                 NUMB_BC2_1, NUMB_BC2_2, MASS_BC2_SULF, MASS_BC2_BCAR,
!    &                 NUMB_BC3_1, NUMB_BC3_2, MASS_BC3_SULF, MASS_BC3_BCAR,
!    &                 NUMB_OCS_1, NUMB_OCS_2, MASS_OCS_SULF,                MASS_OCS_OCAR,
!    &                 NUMB_DBC_1, NUMB_DBC_2, MASS_DBC_SULF, MASS_DBC_BCAR,                MASS_DBC_DUST,
!    &                 NUMB_BOC_1, NUMB_BOC_2, MASS_BOC_SULF, MASS_BOC_BCAR, MASS_BOC_OCAR,
!    &                 NUMB_BCS_1, NUMB_BCS_2, MASS_BCS_SULF, MASS_BCS_BCAR,
!    &                 NUMB_MXX_1, NUMB_MXX_2, MASS_MXX_SULF, MASS_MXX_BCAR, MASS_MXX_OCAR, MASS_MXX_DUST, MASS_MXX_SEAS
!-------------------------------------------------------------------------------------------------------------------------
   90 FORMAT(I6,9X,A3,9X,I4,12X,A4,I19,5X,A16)
      RETURN
      END SUBROUTINE SETUP_CONFIG


      SUBROUTINE Setup_Indices(i,j,index,in)
!-------------------------------------------------------------------------------
!     Routine to initialize indices of the AERO array.
!-------------------------------------------------------------------------------
      implicit none
      integer :: i, j, index, in
      integer, parameter :: omit = 0
      if ( firstime_indices ) then
        firstime_indices = .false.
        !-----------------------------------------------------------------------
        ! set all aero indices to default values indicating that the species 
        ! is not defined in the current mechanism.
        !-----------------------------------------------------------------------
        mass_akk_sulf = omit
        numb_akk_1    = omit
        mass_acc_sulf = omit
        numb_acc_1    = omit
        mass_dd1_sulf = omit
        mass_dd1_dust = omit
        numb_dd1_1    = omit
        mass_ds1_sulf = omit
        mass_ds1_dust = omit
        numb_ds1_1    = omit
        mass_dd2_sulf = omit
        mass_dd2_dust = omit
        numb_dd2_1    = omit
        mass_ds2_sulf = omit
        mass_ds2_dust = omit
        numb_ds2_1    = omit
        mass_ssa_sulf = omit
        mass_ssa_seas = omit
        numb_ssa_1    = omit
        mass_ssc_sulf = omit
        mass_ssc_seas = omit
        numb_ssc_1    = omit
        mass_sss_sulf = omit
        mass_sss_seas = omit
        numb_sss_1    = omit
        mass_occ_sulf = omit
        mass_occ_ocar = omit
        numb_occ_1    = omit
        mass_bc1_sulf = omit
        mass_bc1_bcar = omit
        numb_bc1_1    = omit
        mass_bc2_sulf = omit
        mass_bc2_bcar = omit
        numb_bc2_1    = omit
        mass_bc3_sulf = omit
        mass_bc3_bcar = omit
        numb_bc3_1    = omit
        mass_dbc_sulf = omit
        mass_dbc_bcar = omit
        mass_dbc_dust = omit
        numb_dbc_1    = omit
        mass_boc_sulf = omit
        mass_boc_bcar = omit
        mass_boc_ocar = omit
        numb_boc_1    = omit
        mass_bcs_sulf = omit
        mass_bcs_bcar = omit
        numb_bcs_1    = omit
        mass_ocs_sulf = omit
        mass_ocs_ocar = omit
        numb_ocs_1    = omit
        mass_mxx_sulf = omit
        mass_mxx_bcar = omit
        mass_mxx_ocar = omit
        mass_mxx_dust = omit
        mass_mxx_seas = omit
        numb_mxx_1    = omit
      endif
      if ( in .eq. 0 ) then       ! this is a mass concentration.
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'AKK_SULF' ) mass_akk_sulf = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'ACC_SULF' ) mass_acc_sulf = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'DD1_SULF' ) mass_dd1_sulf = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'DD1_DUST' ) mass_dd1_dust = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'DS1_SULF' ) mass_ds1_sulf = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'DS1_DUST' ) mass_ds1_dust = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'DD2_SULF' ) mass_dd2_sulf = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'DD2_DUST' ) mass_dd2_dust = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'DS2_SULF' ) mass_ds2_sulf = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'DS2_DUST' ) mass_ds2_dust = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'SSA_SULF' ) mass_ssa_sulf = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'SSA_SEAS' ) mass_ssa_seas = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'SSC_SULF' ) mass_ssc_sulf = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'SSC_SEAS' ) mass_ssc_seas = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'SSS_SULF' ) mass_sss_sulf = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'SSS_SEAS' ) mass_sss_seas = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'OCC_SULF' ) mass_occ_sulf = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'OCC_OCAR' ) mass_occ_ocar = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'BC1_SULF' ) mass_bc1_sulf = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'BC1_BCAR' ) mass_bc1_bcar = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'BC2_SULF' ) mass_bc2_sulf = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'BC2_BCAR' ) mass_bc2_bcar = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'BC3_SULF' ) mass_bc3_sulf = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'BC3_BCAR' ) mass_bc3_bcar = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'OCS_SULF' ) mass_ocs_sulf = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'OCS_OCAR' ) mass_ocs_ocar = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'DBC_SULF' ) mass_dbc_sulf = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'DBC_BCAR' ) mass_dbc_bcar = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'DBC_DUST' ) mass_dbc_dust = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'BOC_SULF' ) mass_boc_sulf = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'BOC_BCAR' ) mass_boc_bcar = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'BOC_OCAR' ) mass_boc_ocar = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'BCS_SULF' ) mass_bcs_sulf = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'BCS_BCAR' ) mass_bcs_bcar = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'MXX_SULF' ) mass_mxx_sulf = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'MXX_BCAR' ) mass_mxx_bcar = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'MXX_OCAR' ) mass_mxx_ocar = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'MXX_DUST' ) mass_mxx_dust = index
        if ( mode_name(i)//'_'//chem_spc_name(j) .eq. 'MXX_SEAS' ) mass_mxx_seas = index
      elseif ( in .eq. 1 ) then   ! this is a number concentration.
        if ( mode_name(i).eq.'AKK' .AND. j.eq.1 ) numb_akk_1 = index     
        if ( mode_name(i).eq.'AKK' .AND. j.eq.2 ) numb_akk_2 = index     
        if ( mode_name(i).eq.'ACC' .AND. j.eq.1 ) numb_acc_1 = index     
        if ( mode_name(i).eq.'ACC' .AND. j.eq.2 ) numb_acc_2 = index     
        if ( mode_name(i).eq.'DD1' .AND. j.eq.1 ) numb_dd1_1 = index     
        if ( mode_name(i).eq.'DD1' .AND. j.eq.2 ) numb_dd1_2 = index     
        if ( mode_name(i).eq.'DS1' .AND. j.eq.1 ) numb_ds1_1 = index     
        if ( mode_name(i).eq.'DS1' .AND. j.eq.2 ) numb_ds1_2 = index     
        if ( mode_name(i).eq.'DD2' .AND. j.eq.1 ) numb_dd2_1 = index     
        if ( mode_name(i).eq.'DD2' .AND. j.eq.2 ) numb_dd2_2 = index     
        if ( mode_name(i).eq.'DS2' .AND. j.eq.1 ) numb_ds2_1 = index     
        if ( mode_name(i).eq.'DS2' .AND. j.eq.2 ) numb_ds2_2 = index     
        if ( mode_name(i).eq.'SSA' .AND. j.eq.1 ) numb_ssa_1 = index     
        if ( mode_name(i).eq.'SSA' .AND. j.eq.2 ) numb_ssa_2 = index     
        if ( mode_name(i).eq.'SSC' .AND. j.eq.1 ) numb_ssc_1 = index     
        if ( mode_name(i).eq.'SSC' .AND. j.eq.2 ) numb_ssc_2 = index     
        if ( mode_name(i).eq.'SSS' .AND. j.eq.1 ) numb_sss_1 = index     
        if ( mode_name(i).eq.'SSS' .AND. j.eq.2 ) numb_sss_2 = index     
        if ( mode_name(i).eq.'OCC' .AND. j.eq.1 ) numb_occ_1 = index     
        if ( mode_name(i).eq.'OCC' .AND. j.eq.2 ) numb_occ_2 = index     
        if ( mode_name(i).eq.'BC1' .AND. j.eq.1 ) numb_bc1_1 = index     
        if ( mode_name(i).eq.'BC1' .AND. j.eq.2 ) numb_bc1_2 = index     
        if ( mode_name(i).eq.'BC2' .AND. j.eq.1 ) numb_bc2_1 = index     
        if ( mode_name(i).eq.'BC2' .AND. j.eq.2 ) numb_bc2_2 = index     
        if ( mode_name(i).eq.'BC3' .AND. j.eq.1 ) numb_bc3_1 = index     
        if ( mode_name(i).eq.'BC3' .AND. j.eq.2 ) numb_bc3_2 = index     
        if ( mode_name(i).eq.'OCS' .AND. j.eq.1 ) numb_ocs_1 = index     
        if ( mode_name(i).eq.'OCS' .AND. j.eq.2 ) numb_ocs_2 = index     
        if ( mode_name(i).eq.'DBC' .AND. j.eq.1 ) numb_dbc_1 = index     
        if ( mode_name(i).eq.'DBC' .AND. j.eq.2 ) numb_dbc_2 = index     
        if ( mode_name(i).eq.'BOC' .AND. j.eq.1 ) numb_boc_1 = index     
        if ( mode_name(i).eq.'BOC' .AND. j.eq.2 ) numb_boc_2 = index     
        if ( mode_name(i).eq.'BCS' .AND. j.eq.1 ) numb_bcs_1 = index     
        if ( mode_name(i).eq.'BCS' .AND. j.eq.2 ) numb_bcs_2 = index     
        if ( mode_name(i).eq.'MXX' .AND. j.eq.1 ) numb_mxx_1 = index     
        if ( mode_name(i).eq.'MXX' .AND. j.eq.2 ) numb_mxx_2 = index     
      ENDIF
      RETURN
      END SUBROUTINE Setup_Indices


      SUBROUTINE Setup_Species_Maps
!-------------------------------------------------------------------------------
!     Routine to ... assign a mode number to each mode.
!                ... setup mass maps for species SULF, BCAR, OCAR, DUST, SEAS.
!                ... setup the number of mass species for mode I, NM(I).
!                ... setup the mass species names for mode I, NM_SPC_NAME(I).
!-------------------------------------------------------------------------------
      implicit none
      integer :: i,j
      integer :: isulf,ibcar,iocar,idust,iseas,inm
      integer, parameter :: inactive = 0

      !-----------------------------------------------------------------------
      ! get the number of modes containing each species: sulf, bcar, ocar,
      !   dust, seas.
      !-----------------------------------------------------------------------
      isulf = 0
      ibcar = 0
      iocar = 0
      idust = 0
      iseas = 0

      do i=1, nmodes
        do j=1, naerobox
          if ( aero_spcs(j)(1:13).eq.'MASS_'//mode_name(i)//'_SULF' ) isulf=isulf+1
          if ( aero_spcs(j)(1:13).eq.'MASS_'//mode_name(i)//'_BCAR' ) ibcar=ibcar+1
          if ( aero_spcs(j)(1:13).eq.'MASS_'//mode_name(i)//'_OCAR' ) iocar=iocar+1
          if ( aero_spcs(j)(1:13).eq.'MASS_'//mode_name(i)//'_DUST' ) idust=idust+1
          if ( aero_spcs(j)(1:13).eq.'MASS_'//mode_name(i)//'_SEAS' ) iseas=iseas+1
        ENDDO
      ENDDO

      if( write_log ) then
        WRITE(AUNIT1,'(/A,5I4,/)')'Number of modes containing each species: ISULF,IBCAR,IOCAR,IDUST,ISEAS=', &
                                                                            isulf,ibcar,iocar,idust,iseas
      ENDIF
      nmodes_sulf = isulf
      nmodes_bcar = ibcar
      nmodes_ocar = iocar
      nmodes_dust = idust
      nmodes_seas = iseas

      allocate ( sulf_map(nmodes_sulf) ) 
      allocate ( bcar_map(nmodes_bcar) ) 
      allocate ( ocar_map(nmodes_ocar) ) 
      allocate ( dust_map(nmodes_dust) ) 
      allocate ( seas_map(nmodes_seas) ) 
      allocate ( mode_numb_seas(nmodes_seas) ) 
      sulf_map(:) = 0
      bcar_map(:) = 0
      ocar_map(:) = 0
      dust_map(:) = 0
      seas_map(:) = 0
      mode_numb_seas(:) = 0
      !-------------------------------------------------------------------------
      ! Assign a mode number to each mode.
      !-------------------------------------------------------------------------
      mode_numb_akk = inactive
      mode_numb_acc = inactive
      mode_numb_dd1 = inactive
      mode_numb_dd2 = inactive
      mode_numb_ds1 = inactive
      mode_numb_ds2 = inactive
      mode_numb_ssa = inactive
      mode_numb_ssc = inactive
      mode_numb_sss = inactive
      mode_numb_occ = inactive
      mode_numb_bc1 = inactive
      mode_numb_bc2 = inactive
      mode_numb_bc3 = inactive
      mode_numb_dbc = inactive
      mode_numb_boc = inactive
      mode_numb_bcs = inactive
      mode_numb_ocs = inactive
      mode_numb_mxx = inactive
      do i=1, nmodes
        IF ( mode_name(i) .EQ. 'AKK' ) mode_numb_akk = i
        IF ( mode_name(i) .EQ. 'ACC' ) mode_numb_acc = i
        IF ( mode_name(i) .EQ. 'DD1' ) mode_numb_dd1 = i
        IF ( mode_name(i) .EQ. 'DD2' ) mode_numb_dd2 = i
        IF ( mode_name(i) .EQ. 'DS1' ) mode_numb_ds1 = i
        IF ( mode_name(i) .EQ. 'DS2' ) mode_numb_ds2 = i
        IF ( mode_name(i) .EQ. 'SSA' ) mode_numb_ssa = i
        IF ( mode_name(i) .EQ. 'SSC' ) mode_numb_ssc = i
        IF ( mode_name(i) .EQ. 'SSS' ) mode_numb_sss = i
        IF ( mode_name(i) .EQ. 'OCC' ) mode_numb_occ = i
        IF ( mode_name(i) .EQ. 'BC1' ) mode_numb_bc1 = i
        IF ( mode_name(i) .EQ. 'BC2' ) mode_numb_bc2 = i
        IF ( mode_name(i) .EQ. 'BC3' ) mode_numb_bc3 = i
        IF ( mode_name(i) .EQ. 'DBC' ) mode_numb_dbc = i
        IF ( mode_name(i) .EQ. 'BOC' ) mode_numb_boc = i
        IF ( mode_name(i) .EQ. 'BCS' ) mode_numb_bcs = i
        IF ( mode_name(i) .EQ. 'OCS' ) mode_numb_ocs = i
        IF ( mode_name(i) .EQ. 'MXX' ) mode_numb_mxx = i
      ENDDO  

      if ( write_log ) then
        write(aunit1,'(/a/)')        'mode_name( mode_numb_xxx ), mode_numb_xxx'
        if ( mode_numb_akk .gt. 0 ) then
          write(aunit1,'(3x,a3,3x,i4)') mode_name( mode_numb_akk ), mode_numb_akk
        endif
        if ( mode_numb_acc .gt. 0 ) then
          write(aunit1,'(3x,a3,3x,i4)') mode_name( mode_numb_acc ), mode_numb_acc
        endif
        if ( mode_numb_dd1 .gt. 0 ) then
          write(aunit1,'(3x,a3,3x,i4)') mode_name( mode_numb_dd1 ), mode_numb_dd1
        endif
        if ( mode_numb_ds1 .gt. 0 ) then
          write(aunit1,'(3x,a3,3x,i4)') mode_name( mode_numb_ds1 ), mode_numb_ds1
        endif
        if ( mode_numb_dd2 .gt. 0 ) then
          write(aunit1,'(3x,a3,3x,i4)') mode_name( mode_numb_dd2 ), mode_numb_dd2
        endif
        if ( mode_numb_ds2 .gt. 0 ) then
          write(aunit1,'(3x,a3,3x,i4)') mode_name( mode_numb_ds2 ), mode_numb_ds2
        endif
        if ( mode_numb_ssa .gt. 0 ) then
          write(aunit1,'(3x,a3,3x,i4)') mode_name( mode_numb_ssa ), mode_numb_ssa
        endif
        if ( mode_numb_ssc .gt. 0 ) then
          write(aunit1,'(3x,a3,3x,i4)') mode_name( mode_numb_ssc ), mode_numb_ssc
        endif
        if ( mode_numb_sss .gt. 0 ) then
          write(aunit1,'(3x,a3,3x,i4)') mode_name( mode_numb_sss ), mode_numb_sss
        endif
        if ( mode_numb_occ .gt. 0 ) then
          write(aunit1,'(3x,a3,3x,i4)') mode_name( mode_numb_occ ), mode_numb_occ
        endif
        if ( mode_numb_bc1 .gt. 0 ) then
          write(aunit1,'(3x,a3,3x,i4)') mode_name( mode_numb_bc1 ), mode_numb_bc1
        endif
        if ( mode_numb_bc2 .gt. 0 ) then
          write(aunit1,'(3x,a3,3x,i4)') mode_name( mode_numb_bc2 ), mode_numb_bc2
        endif
        if ( mode_numb_bc3 .gt. 0 ) then
          write(aunit1,'(3x,a3,3x,i4)') mode_name( mode_numb_bc3 ), mode_numb_bc3
        endif
        if ( mode_numb_dbc .gt. 0 ) then
          write(aunit1,'(3x,a3,3x,i4)') mode_name( mode_numb_dbc ), mode_numb_dbc
        endif
        if ( mode_numb_boc .gt. 0 ) then
          write(aunit1,'(3x,a3,3x,i4)') mode_name( mode_numb_boc ), mode_numb_boc
        endif
        if ( mode_numb_bcs .gt. 0 ) then
          write(aunit1,'(3x,a3,3x,i4)') mode_name( mode_numb_bcs ), mode_numb_bcs
        endif
        if ( mode_numb_ocs .gt. 0 ) then
          write(aunit1,'(3x,a3,3x,i4)') mode_name( mode_numb_ocs ), mode_numb_ocs
        endif
        if ( mode_numb_mxx .gt. 0 ) then
          write(aunit1,'(3x,a3,3x,i4)') mode_name( mode_numb_mxx ), mode_numb_mxx
        endif
      endif
      !-------------------------------------------------------------------------
      ! Setup NUMB_MAP(I): location of the Ith number  conc. in the AERO array
      ! Setup SULF_MAP(I): location of the Ith sulfate conc. in the AERO array
      ! Setup BCAR_MAP(I): location of the Ith BC      conc. in the AERO array
      ! Setup OCAR_MAP(I): location of the Ith OC      conc. in the AERO array
      ! Setup DUST_MAP(I): location of the Ith dust    conc. in the AERO array
      ! Setup SEAS_MAP(I): location of the Ith sea salt conc. in the AERO array
      ! Setup NM(I):       number of mass concs. defined for mode I      
      ! Setup NM_SPC_NAME(I,J): name of the Jth mass conc. in mode I
      !-------------------------------------------------------------------------
      ibcar = 1
      iocar = 1
      idust = 1
      iseas = 1
      nm(:) = 0
      nm_spc_name(:,:) = '    '     ! len=4 character variable.

      do i=1, nmodes
        inm = 0
        do j=1, naerobox

          if ( aero_spcs(j)(1:8)  .eq. 'NUMB_'//mode_name(i) )          numb_map(i) = j
          if ( aero_spcs(j)(1:13) .eq. 'MASS_'//mode_name(i)//'_SULF' ) sulf_map(i) = j
          if ( aero_spcs(j)(1:13) .eq. 'MASS_'//mode_name(i)//'_BCAR' ) THEN
            bcar_map(ibcar) = j          ! location of this mass conc. in the aero array
            ibcar = ibcar + 1
          endif

          if ( aero_spcs(j)(1:13) .eq. 'MASS_'//mode_name(i)//'_OCAR' ) THEN
            ocar_map(iocar) = j          ! location of this mass conc. in the AERO array
            iocar = iocar + 1
          endif

          if ( aero_spcs(j)(1:13) .eq. 'MASS_'//mode_name(i)//'_DUST' ) THEN
            dust_map(idust) = j          ! location of this mass conc. in the AERO array
            idust = idust + 1
          endif

          if ( aero_spcs(j)(1:13) .eq. 'MASS_'//mode_name(i)//'_SEAS' ) THEN
            seas_map(iseas) = j          ! location of this mass conc. in the AERO array
            mode_numb_seas(iseas) = i    ! mode number for this sea salt-containing mode
            iseas = iseas + 1
          ENDIF

          if ( aero_spcs(j)(1:8) .eq. 'MASS_'//mode_name(i) ) THEN

            inm = inm + 1
            nm_spc_name(i,inm) = aero_spcs(j)(10:13)

          ENDIF

        ENDDO
      
        nm(i) = inm

      ENDDO

      IF ( write_log ) then
        WRIte(aunit1,'(/a11,16i4)') 'NUMB_MAP = ', numb_map(:)
        WRIte(aunit1,'(/a11,16i4)') 'SULF_MAP = ', sulf_map(:)
        WRIte(aunit1,'(/a11,16i4)') 'BCAR_MAP = ', bcar_map(:)
        WRIte(aunit1,'(/a11,16i4)') 'OCAR_MAP = ', ocar_map(:)
        WRIte(aunit1,'(/a11,16i4)') 'DUST_MAP = ', dust_map(:)
        WRIte(aunit1,'(/a11,16i4)') 'SEAS_MAP = ', seas_map(:)
        WRIte(aunit1,'(/a11,16i4)') 'NM(I)    = ', nm(:)
        do i=1, nmodes
          write(aunit1,'(/A40,A5,I4,5A5)') 'MODE_NAME(I), NM(I), NM_SPC_NAME(I,:) = ', mode_name(i), nm(i), nm_spc_name(i,:)
        enddo
      endif
      !----------------------------------------------------------------------------------------------------
      ! Setup EMIS_MODE_MAP. There are presently 10 emitted species.
      !
      ! EMIS_MODE_MAP and EMIS_SPCS_MAP have elements corresponding to the aerosol types (in this order):
      !   AKK(=1), ACC(=2), BCC(=8), OCC(=7), DD1(=3), SSA(=5), SSC(=6), BOC(BC=8), BOC(OC=9), DD2(=10).
      !
      ! EMIS_MODE_MAP(J) is mode number receiving the emissions held in EMIS_MASS(J).
      ! EMIS_SPCS_MAP(J) is the chemical species number (1-5) of the chemical species held in EMIS_MASS(J).
      ! EMIS_SPCS_MAP = (/1,1,2,3,4,5,5,2,3,4/) is set at the top of the module.
      !----------------------------------------------------------------------------------------------------
      emis_mode_map(:) = 0 
      emis_mode_map(1) = 1  ! Aitken mode sulfate always goes in the first mode, whether it is AKK or ACC.

      do i=1, nmodes        ! If no Aitken mode, then the accumulation mode is the first mode.
        if( mode_name(i) .eq. 'ACC' ) emis_mode_map(2) = i
        if( mode_name(i) .eq. 'BC1' ) emis_mode_map(3) = i
        if( mode_name(i) .eq. 'OCC' ) emis_mode_map(4) = i
        if( mode_name(i) .eq. 'DD1' ) THEN
          emis_mode_map(5) = I
          if( mech .ge. 5 .and. mech .le. 8 ) emis_mode_map(10) = i   ! emissions for both dust modes go into mode DD1
        endif
        if( mode_name(i) .eq. 'SSA' ) emis_mode_map(6) = i
        if( mode_name(i) .eq. 'SSC' ) emis_mode_map(7) = i
        if( mode_name(i) .eq. 'SSS' ) emis_mode_map(6) = i  ! emissions for both sea salt modes go into mode SSS
        if( mode_name(i) .eq. 'SSS' ) emis_mode_map(7) = i  ! emissions for both sea salt modes go into mode SSS
        if( mode_name(i) .eq. 'BOC' ) emis_mode_map(8) = i
        if( mode_name(i) .eq. 'BOC' ) emis_mode_map(9) = i
        if( mode_name(i) .eq. 'DD2' ) emis_mode_map(10) = i
      ENDDO

      !-------------------------------------------------------------------------
      ! If the mechanism does not have the mode BOC to receive the 
      ! mixed BC-OC emissions, put these directly into the BC1 and OCC modes.
      !-------------------------------------------------------------------------
      if ( emis_mode_map(8) .eq. 0 ) then   ! the bc in bc-oc emissions
        do i=1, nmodes   
          if( mode_name(i) .eq. 'BC1' ) then
            emis_mode_map(8) = i
            if ( write_log ) write(aunit1,'(/2A/)')'BC of BO-OC put into mode ',mode_name(i)
          endif
        ENDDO
      ENDIF
      if ( emis_mode_map(9) .eq. 0 ) then   ! the oc in bc-oc emissions
        do i=1, nmodes   
          if( mode_name(i) .eq. 'OCC' ) then
            emis_mode_map(9) = i
            if ( write_log ) write(aunit1,'(/2A/)')'OC of BO-OC put into mode ',mode_name(i)
          ENDIF
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE Setup_Species_Maps


      SUBROUTINE Setup_Seasalt_Maps
!---------------------------------------------------------------------------------------------------------
!     Routine to setup maps and variables needed to derive sea salt 
!     number concentrations from sea salt mass concentrations and
!     characteristic mean masses per particle.
!
!     Three arrays are set up:
!
!       SEAS_MODE_MAP(II)      ! mode number of the IIth SS mode
!       SEAS_MODE_MASS_MAP(II) ! location in the AERO array of the SS mass for the IIth SS mode
!       RECIP_SEAS_MPP(II)     ! reciprocal of the mean mass per particle for the IIth SS mode
!
!     There are either two sea salt modes, SSA and SSC, or only one, SSS.
!---------------------------------------------------------------------------------------------------------
      implicit none
      integer :: i,idim,ii,j
      real(8) :: dps
      logical, save :: firstime = .true.
      if ( firstime ) then
        firstime = .false.
        idim = 0
        do i=1, nmodes
          if ( mode_name(i)(1:2) .eq. 'SS' ) idim = idim + 1
        ENDDO
        if ( idim .lt. 1 .or. idim .gt. 2 ) then
          WRITE(*,*)'SUBROUTIINE SETUP_SEASALT_MAPS:'
          WRITE(*,*)'NUMBER OF SEAS SALT MODES IS ZERO OR GREATER THAN TWO - INCORRECT'
          STOP
        ENDIF
        number_of_seasalt_modes = idim
   !     write (*,FMT='(A,1(I3.3,1X))') 'LFR=DBG Setup: ', &
   !     idim; CALL flush(6)
        allocate( seas_mode_map( idim ) )     
        allocate( seas_mode_mass_map( idim ) )
        allocate( recip_seas_mpp( idim ) )    
        seas_mode_map(:) = 0
        seas_mode_mass_map(:) = 0
        recip_seas_mpp(:) = 0.0d+00
        ii = 0
        do i=1, nmodes
          if ( mode_name(i)(1:2) .eq. 'SS' ) THEN  ! This is a sea salt mode, either SSA, SSC, or SSS.
            ii = ii + 1
            seas_mode_map(ii) = i
            do j=1, naerobox
              if( aero_spcs(j)(1:13) .eq. 'MASS_'//MODE_NAME(I)//'_SEAS' ) then ! this is the sea salt mass in the mode.
                seas_mode_mass_map(ii) = j
                ! write(*,*) aero_spcs(j)(1:13)     ! checked: the correct species are identified. 
              ENDIF
            ENDDO      
            !-------------------------------------------------------------------------------------------------
            ! Calculate the diameter of average mass for each sea salt mode for the emissions lognormals.
            !
            ! The emissions lognormals are used since it is the dry sea salt concentration that
            ! is divided by the average dry mass per particle (in subr. matrix) to obtain the current
            ! number concentration from the dry sea salt mass concentration for each sea salt mode.
            ! The dry NaCl mass per particle changes little over time for sea salt, given low coagulation rates.
            ! Accreted water, sulfate, nitrate, and ammonium are irrelevant here. The dry NaCl per emitted
            ! sea salt particle is the appropriate dry mass to divide into the current dry NaCl concentration
            ! to obtain the particle number concentration.
            ! 
            ! Dam [um] = ( diameter moment 3  /  diameter moment 0    )**(1/3)
            !          = ( dg**3 * sg**9                              )**(1/3)
            !          = ( dg**3 * [ exp( 0.5*(log(sigmag))**2 ) ]**9 )**(1/3)
            !          =   dg    * [ exp( 0.5*(log(sigmag))**2 ) ]**3
            !          =   dg    * [ exp( 1.5*(log(sigmag))**2 ) ]
            !-------------------------------------------------------------------------------------------------
            if( mode_name(i).eq.'SSA') dps = 1.0d-06 * dg_ssa_emis * exp( 1.5d+00 * ( log(sg_ssa_emis) )**2 )
            if( mode_name(i).eq.'SSC') dps = 1.0d-06 * dg_ssc_emis * exp( 1.5d+00 * ( log(sg_ssc_emis) )**2 )
            if( mode_name(i).eq.'SSS') dps = 1.0d-06 * dg_sss_emis * exp( 1.5d+00 * ( log(sg_sss_emis) )**2 )
            if( activation_comparison .and. ( mech.eq.4 .or. mech.eq.8 ) ) then   ! activation test. 
              dps = 1.0d-06 * 0.3308445477d+00 * exp( 1.5d+00 * ( log( sg_sss_emis) )**2 )
              WRITE(*,*)'Special value set for DPS and RECIP_SEAS_MPP for mode SSS in aero_setup.f.'
            ENDIF
            !-------------------------------------------------------------------------------------------------
            ! DPS is the diameter of average mass for mode J for the dry emitted sea salt, and when cubed
            ! and multiplied by pi/6 it yields the average dry sea salt particle volume in emissions mode J.
            ! Multiplication by the emitted dry particle sea salt density then yields the average
            ! dry sea salt mass per particle emitted into the mode. 
            !-------------------------------------------------------------------------------------------------
            if( discrete_eval_option ) then
              recip_seas_mpp(ii) = 1.0d+00 / ( 1.0d+12 * densp          * pi6 * dps**3 )
            else
              recip_seas_mpp(ii) = 1.0d+00 / ( 1.0d+12 * emis_dens_seas * pi6 * dps**3 )
            endif
            if( write_log ) then
              WRITE(AUNIT1,'(/A/)')'II,I,MODE_NAME(I),AERO_SPCS(SEAS_MODE_MASS_MAP(II)),RECIP_SEAS_MPP(II)'
              WRITE(AUNIT1,90000)   ii,i,mode_name(i),aero_spcs(seas_mode_mass_map(ii)),recip_seas_mpp(ii)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
90000 FORMAT(2I4,4X,A6,4X,A16,4X,D15.5)
      RETURN
      END SUBROUTINE Setup_Seasalt_Maps


      SUBROUTINE Setup_Emis        
!-------------------------------------------------------------------------------
!     Routine to calculate the reciprocal of the average particle mass [ug]
!     for each mode for conversion of mode mass emission rates [ug/m^3/s] 
!     to mode number emission rates [#/m^3/s].
!
!     The factor 1.0D+12 converts [m^3] to [cm^3] and [g] to [ug].
!     DP0_EMIS(:) is in [m].
!
!     EMIS_MODE_MAP has elements corresponding to the aerosol types (in this order):
!       AKK(=1), ACC(=2), BCC(=8), OCC(=7),
!       DD1(=3), SSA(=5), SSC(=6), BOC(BC=8), BOC(OC=9), DD2(=10).
!     EMIS_MODE_MAP(J) is mode number receiving the emissions held in EMIS_MASS(J).
!-------------------------------------------------------------------------------
      IMPLICIT NONE
      integer :: i, j

      recip_part_mass(:) = 0.0d+00    ! [1/ug]

      if ( write_log ) then
        WRITE(AUNIT1,'(/2A5,A17,A32,A32/)') 'I','J','  DP0_EMIS(J)[um]', &
                     'RECIP_PART_MASS(I)[1/ug]','EMIS_DENS(I)[g/cm^3]'
      ENDIF
      do i=1, nemis_spcs       ! currently, nemis_spcs = 10 emitted species
        j = emis_mode_map(i)   ! j ranges over the number of modes.
        !-----------------------------------------------------------------------
        ! dp0_emis(j) is the diameter of average mass for mode j, and when cubed
        ! and multiplied by pi/6 it yields the average particle volume in mode j.
        ! multiplication by the emitted particle density then yields the average
        ! mass per particle emitted into the mode.
        ! 
        ! emis_dens(:)       ranges over the emission species only. 
        ! recip_part_mass(:) ranges over the emission species only. 
        ! dp0_emis(:)        ranges over all modes in the mechanism.
        !-----------------------------------------------------------------------
        recip_part_mass(i) = 1.0d+00 / ( 1.0d+12 * emis_dens(i) * pi6 * dp0_emis(j)**3 )
        if( write_log) write(aunit1,90000) i, j, dp0_emis(j)*1.0d+06, recip_part_mass(i), emis_dens(i)
      enddo

90000 FORMAT(2I5,F17.6,D32.6,F32.4)
      RETURN
      END SUBROUTINE Setup_Emis


      SUBROUTINE Setup_Dp0         
!-------------------------------------------------------------------------------
!     Routine to calculate the diameter of average mass [m] for each mode.
!-------------------------------------------------------------------------------
      implicit none
      integer :: i

      !-------------------------------------------------------------------------
      ! set lognormal parameters and activating fraction for each mode.
      !-------------------------------------------------------------------------
      dgn0(:)      = 0.0d+00
      sig0(:)      = 0.0d+00
      lnsig0(:)    = 0.0d+00
      dgn0_emis(:) = 0.0d+00
      sig0_emis(:) = 0.0d+00
      kappai(:)    = 0.0d+00
!KML
      kappapk(:)   = 0.0d+00      
      dp0(:)       = 0.0d+00
      dp0_emis(:)  = 0.0d+00
      do i=1, nmodes
        if ( mode_name(i) .eq. 'AKK') dgn0(i) = dg_akk     ! dgn0 for characteristic lognormal
        if ( mode_name(i) .eq. 'ACC') then
          if( numb_akk_1 .eq. 0 ) then
            if ( activation_comparison .and. ( mech .eq. 4 .or. mech .eq. 8 ) ) then  ! no mode akk 
              dgn0(i) = 0.0506736714d+00        ! for droplet activation test. 
              WRITE(*,*)'DGN0 for mode ACC set to special value of activation test.'
            else
              dgn0(i) = sqrt( dg_akk * dg_acc ) ! no mode akk; reduce dg_akk
              WRITE(*,*)'DGN0 for mode ACC reduced in this mechanism.'
            ENDIF
          ELSE
            dgn0(i) = dg_acc
          endif
        endif
        if ( mode_name(i) .eq. 'DD1') THEN
          if ( activation_comparison .and. ( mech .ge. 5 .and. mech .le. 8 ) ) then  ! no mode dd2 or ds2 
            dgn0(i) = 0.4516304494d+00        ! for droplet activation test.
            WRITE(*,*)'DGN0 for mode DD1 set to special value of activation test.'
          else
            dgn0(i) = dg_dd1
          endif
        endif
        If ( mode_name(i) .eq. 'DD2') dgn0(i) = dg_dd2
        If ( mode_name(i) .eq. 'DS1') dgn0(i) = dg_ds1
        If ( mode_name(i) .eq. 'DS2') dgn0(i) = dg_ds2
        If ( mode_name(i) .eq. 'SSA') dgn0(i) = dg_ssa
        If ( mode_name(i) .eq. 'SSC') dgn0(i) = dg_ssc
        If ( mode_name(i) .eq. 'SSS') THEN 
          if ( activation_comparison .and. ( mech .eq. 4 .or. mech .eq. 8 ) ) then  ! no mode ssa or ssc 
            dgn0(i) = 0.3308445477d+00        ! for droplet activation test. 
            WRITE(*,*)'DGN0 for mode SSS set to special value for activation test.'
          else
            dgn0(i) = dg_sss
          endif 
        endif
        IF ( mode_name(i) .eq. 'OCC') dgn0(i) = dg_occ
        IF ( mode_name(i) .eq. 'BC1') dgn0(i) = dg_bc1
        IF ( mode_name(i) .eq. 'BC2') dgn0(i) = dg_bc2
        IF ( mode_name(i) .eq. 'BC3') dgn0(i) = dg_bc3
        IF ( mode_name(i) .eq. 'DBC') dgn0(i) = dg_dbc
        IF ( mode_name(i) .eq. 'BOC') dgn0(i) = dg_boc
        IF ( mode_name(i) .eq. 'BCS') dgn0(i) = dg_bcs
        IF ( mode_name(i) .eq. 'OCS') dgn0(i) = dg_ocs
        IF ( mode_name(i) .eq. 'MXX') dgn0(i) = dg_mxx
        IF ( mode_name(i) .eq. 'AKK') sig0(i) = sg_akk     ! sig0 for characteristic lognormal
        IF ( mode_name(i) .eq. 'ACC') sig0(i) = sg_acc
        IF ( mode_name(i) .eq. 'DD1') sig0(i) = sg_dd1
        IF ( mode_name(i) .eq. 'DD2') sig0(i) = sg_dd2
        IF ( mode_name(i) .eq. 'DS1') sig0(i) = sg_ds1
        IF ( mode_name(i) .eq. 'DS2') sig0(i) = sg_ds2
        IF ( mode_name(i) .eq. 'SSA') sig0(i) = sg_ssa
        IF ( mode_name(i) .eq. 'SSC') sig0(i) = sg_ssc
        IF ( mode_name(i) .eq. 'SSS') sig0(i) = sg_sss
        IF ( mode_name(i) .eq. 'OCC') sig0(i) = sg_occ
        IF ( mode_name(i) .eq. 'BC1') sig0(i) = sg_bc1
        IF ( mode_name(i) .eq. 'BC2') sig0(i) = sg_bc2
        IF ( mode_name(i) .eq. 'BC3') sig0(i) = sg_bc3
        IF ( mode_name(i) .eq. 'DBC') sig0(i) = sg_dbc
        IF ( mode_name(i) .eq. 'BOC') sig0(i) = sg_boc
        IF ( mode_name(i) .eq. 'BCS') sig0(i) = sg_bcs
        IF ( mode_name(i) .eq. 'OCS') sig0(i) = sg_ocs
        IF ( mode_name(i) .eq. 'MXX') sig0(i) = sg_mxx
        IF ( mode_name(i) .eq. 'AKK') dgn0_emis(i) = dg_akk_emis     ! DGN0_EMIS for emissions lognormal
        IF ( mode_name(i) .eq. 'ACC') THEN
          dgn0_emis(i) = dg_acc_emis
          if( numb_akk_1 .eq. 0 ) then 
            dgn0_emis(i) = sqrt( dg_akk_emis * dg_acc_emis ) ! no mode akk; reduce dg_akk_emis
            WRITE(*,*)'DGN0_EMIS for mode ACC reduced in this mechanism.'
          ENDIF
        ENDIF
        IF ( mode_name(i) .eq. 'DD1') dgn0_emis(i) = dg_dd1_emis
        IF ( mode_name(i) .eq. 'DD2') dgn0_emis(i) = dg_dd2_emis
        IF ( mode_name(i) .eq. 'DS1') dgn0_emis(i) = dg_ds1_emis
        IF ( mode_name(i) .eq. 'DS2') dgn0_emis(i) = dg_ds2_emis
        IF ( mode_name(i) .eq. 'SSA') dgn0_emis(i) = dg_ssa_emis
        IF ( mode_name(i) .eq. 'SSC') dgn0_emis(i) = dg_ssc_emis
        IF ( mode_name(i) .eq. 'SSS') THEN 
          if ( activation_comparison .and. ( mech .eq. 4 .or. mech .eq. 8 ) ) then  ! no mode ssa or ssc 
            dgn0_emis(i) = 0.3308445477d+00        ! for droplet activation test. 
            WRITE(*,*)'DGN0_EMIS for mode SSS set to special value for activation test.'
          else
            dgn0_emis(i) = dg_sss_emis
          endif 
        endif
        IF ( mode_name(i) .eq. 'OCC') dgn0_emis(i) = dg_occ_emis
        IF ( mode_name(i) .eq. 'BC1') dgn0_emis(i) = dg_bc1_emis
        IF ( mode_name(i) .eq. 'BC2') dgn0_emis(i) = dg_bc2_emis
        IF ( mode_name(i) .eq. 'BC3') dgn0_emis(i) = dg_bc3_emis
        IF ( mode_name(i) .eq. 'DBC') dgn0_emis(i) = dg_dbc_emis
        IF ( mode_name(i) .eq. 'BOC') dgn0_emis(i) = dg_boc_emis
        IF ( mode_name(i) .eq. 'BCS') dgn0_emis(i) = dg_bcs_emis
        IF ( mode_name(i) .eq. 'OCS') dgn0_emis(i) = dg_ocs_emis
        IF ( mode_name(i) .eq. 'MXX') dgn0_emis(i) = dg_mxx_emis
        IF ( mode_name(i) .eq. 'AKK') sig0_emis(i) = sg_akk_emis     ! sig0_emis for emissions lognormal
        IF ( mode_name(i) .eq. 'ACC') sig0_emis(i) = sg_acc_emis
        IF ( mode_name(i) .eq. 'DD1') sig0_emis(i) = sg_dd1_emis
        IF ( mode_name(i) .eq. 'DD2') sig0_emis(i) = sg_dd2_emis
        IF ( mode_name(i) .eq. 'DS1') sig0_emis(i) = sg_ds1_emis
        IF ( mode_name(i) .eq. 'DS2') sig0_emis(i) = sg_ds2_emis
        IF ( mode_name(i) .eq. 'SSA') sig0_emis(i) = sg_ssa_emis
        IF ( mode_name(i) .eq. 'SSC') sig0_emis(i) = sg_ssc_emis
        IF ( mode_name(i) .eq. 'SSS') sig0_emis(i) = sg_sss_emis
        IF ( mode_name(i) .eq. 'OCC') sig0_emis(i) = sg_occ_emis
        IF ( mode_name(i) .eq. 'BC1') sig0_emis(i) = sg_bc1_emis
        IF ( mode_name(i) .eq. 'BC2') sig0_emis(i) = sg_bc2_emis
        IF ( mode_name(i) .eq. 'BC3') sig0_emis(i) = sg_bc3_emis
        IF ( mode_name(i) .eq. 'DBC') sig0_emis(i) = sg_dbc_emis
        IF ( mode_name(i) .eq. 'BOC') sig0_emis(i) = sg_boc_emis
        IF ( mode_name(i) .eq. 'BCS') sig0_emis(i) = sg_bcs_emis
        IF ( mode_name(i) .eq. 'OCS') sig0_emis(i) = sg_ocs_emis
        IF ( mode_name(i) .eq. 'MXX') sig0_emis(i) = sg_mxx_emis
        IF ( mode_name(i) .eq. 'AKK') kappai(i) = kappai_akk         ! kappai - activating fraction
        IF ( mode_name(i) .eq. 'ACC') kappai(i) = kappai_acc
        IF ( mode_name(i) .eq. 'DD1') kappai(i) = kappai_dd1
        IF ( mode_name(i) .eq. 'DD2') kappai(i) = kappai_dd2
        IF ( mode_name(i) .eq. 'DS1') kappai(i) = kappai_ds1
        IF ( mode_name(i) .eq. 'DS2') kappai(i) = kappai_ds2
        IF ( mode_name(i) .eq. 'SSA') kappai(i) = kappai_ssa
        IF ( mode_name(i) .eq. 'SSC') kappai(i) = kappai_ssc
        IF ( mode_name(i) .eq. 'SSS') kappai(i) = kappai_sss
        IF ( mode_name(i) .eq. 'OCC') kappai(i) = kappai_occ
        IF ( mode_name(i) .eq. 'BC1') kappai(i) = kappai_bc1
        IF ( mode_name(i) .eq. 'BC2') kappai(i) = kappai_bc2
        IF ( mode_name(i) .eq. 'BC3') kappai(i) = kappai_bc3
        IF ( mode_name(i) .eq. 'DBC') kappai(i) = kappai_dbc
        IF ( mode_name(i) .eq. 'BOC') kappai(i) = kappai_boc
        IF ( mode_name(i) .eq. 'BCS') kappai(i) = kappai_bcs
        IF ( mode_name(i) .eq. 'OCS') kappai(i) = kappai_ocs
        IF ( mode_name(i) .eq. 'MXX') kappai(i) = kappai_mxx

!KML
        IF ( mode_name(i) .eq. 'AKK') kappapk(i) = kappapk_SULF         ! kappapk - hygroscopicity parameter
        IF ( mode_name(i) .eq. 'ACC') kappapk(i) = kappapk_SULF
        IF ( mode_name(i) .eq. 'DD1') kappapk(i) = kappapk_DUST
        IF ( mode_name(i) .eq. 'DD2') kappapk(i) = kappapk_DUST         !To be recalculated
        IF ( mode_name(i) .eq. 'DS1') kappapk(i) = kappapk_DUST
        IF ( mode_name(i) .eq. 'DS2') kappapk(i) = kappapk_DUST         !To be recalculated
        IF ( mode_name(i) .eq. 'SSA') kappapk(i) = kappapk_SEAS
        IF ( mode_name(i) .eq. 'SSC') kappapk(i) = kappapk_SEAS
        IF ( mode_name(i) .eq. 'SSS') kappapk(i) = kappapk_SEAS
        IF ( mode_name(i) .eq. 'OCC') kappapk(i) = kappapk_OCAR
        IF ( mode_name(i) .eq. 'BC1') kappapk(i) = kappapk_BCAR
        IF ( mode_name(i) .eq. 'BC2') kappapk(i) = kappapk_BCAR         !To be recalculated
        IF ( mode_name(i) .eq. 'BC3') kappapk(i) = kappapk_BCAR         !To be recalculated
        IF ( mode_name(i) .eq. 'DBC') kappapk(i) = kappapk_BCAR         
        IF ( mode_name(i) .eq. 'BOC') kappapk(i) = kappapk_OCAR         !To be recalculated
        IF ( mode_name(i) .eq. 'BCS') kappapk(i) = kappapk_BCAR         !To be recalculated
        IF ( mode_name(i) .eq. 'OCS') kappapk(i) = kappapk_OCAR         !To be recalculated
        IF ( mode_name(i) .eq. 'MXX') kappapk(i) = kappapk_SULF         !To be recalculated

        IF ( mode_name(i) .eq. 'AKK') denspi(i) = emis_dens_sulf     ! denspi - default density for mode i
        IF ( mode_name(i) .eq. 'ACC') denspi(i) = emis_dens_sulf
        IF ( mode_name(i) .eq. 'DD1') denspi(i) = emis_dens_dust
        IF ( mode_name(i) .eq. 'DD2') denspi(i) = emis_dens_dust
        IF ( mode_name(i) .eq. 'DS1') denspi(i) = emis_dens_dust
        IF ( mode_name(i) .eq. 'DS2') denspi(i) = emis_dens_dust
        IF ( mode_name(i) .eq. 'SSA') denspi(i) = emis_dens_seas
        IF ( mode_name(i) .eq. 'SSC') denspi(i) = emis_dens_seas
        IF ( mode_name(i) .eq. 'SSS') denspi(i) = emis_dens_seas
        IF ( mode_name(i) .eq. 'OCC') denspi(i) = emis_dens_ocar
        IF ( mode_name(i) .eq. 'BC1') denspi(i) = emis_dens_bcar
        IF ( mode_name(i) .eq. 'BC2') denspi(i) = emis_dens_bcar
        IF ( mode_name(i) .eq. 'BC3') denspi(i) = emis_dens_bcar
        IF ( mode_name(i) .eq. 'DBC') denspi(i) = 0.5d+00 * ( emis_dens_dust + emis_dens_bcar )
        IF ( mode_name(i) .eq. 'BOC') denspi(i) = emis_dens_bocc
        IF ( mode_name(i) .eq. 'BCS') denspi(i) = 0.5d+00 * ( emis_dens_bcar + emis_dens_sulf )
        IF ( mode_name(i) .eq. 'OCS') denspi(i) = 0.5d+00 * ( emis_dens_ocar + emis_dens_sulf )
        IF ( mode_name(i) .eq. 'MXX') denspi(i) = 0.2d+00 * ( emis_dens_sulf + emis_dens_dust &
                                                          +   emis_dens_bcar + emis_dens_ocar + emis_dens_seas )
      ENDDO

      !---------------------------------------------------------------------------------------------------------------------
      ! Calculate the diameter of average mass for both the (default)
      ! characteristic lognormal and the emissions lognormal for each mode.
      ! 
      ! Dam [um] = ( diameter moment 3  /  diameter moment 0    )**(1/3)
      !          = ( dg**3 * sg**9                              )**(1/3)
      !          = ( dg**3 * [ exp( 0.5*(log(sigmag))**2 ) ]**9 )**(1/3)
      !          =   dg    * [ exp( 0.5*(log(sigmag))**2 ) ]**3
      !          =   dg    * [ exp( 1.5*(log(sigmag))**2 ) ]
      !
      ! Also store the natural logarithms of the geo. std. deviations. 
      !---------------------------------------------------------------------------------------------------------------------
      if( write_log ) write(aunit1,'(/8A12/)') 'I','  MODE','DGN0 [um]','SIG0 [1]','DP0 [um]', &
                                               'DGN0_E [um]','SIG0_E [1]','DP0_E [um]'
      do i=1, nweights
        dp0(i)      = 1.0d-06 * dgn0(i)      * exp( 1.5d+00 * ( log(sig0(i))      )**2 )  ! convert from [um] to [m]
        dp0_emis(i) = 1.0d-06 * dgn0_emis(i) * exp( 1.5d+00 * ( log(sig0_emis(i)) )**2 )  ! convert from [um] to [m]
        conv_dpam_to_dgn(i) = 1.0d+00 / exp( 1.5d+00 * ( log(sig0_emis(i)) )**2 ) 
        lnsig0(i) = log( sig0(i) ) 
        if( write_log ) write(aunit1,90000)i,mode_name(i),dgn0(i),sig0(i),dp0(i)*1.0d+06, &
                                           dgn0_emis(i),sig0_emis(i),dp0_emis(i)*1.0d+06
      enddo

      !---------------------------------------------------------------------------------------------------------------------
      ! Set densities and their reciprocals for each chemical component of any mode. 
      !---------------------------------------------------------------------------------------------------------------------
      dens_comp(1) = rho_nh42so4     ! [g/cm^3] sulfate 
      dens_comp(2) = emis_dens_bcar  ! [g/cm^3] bc
      dens_comp(3) = emis_dens_ocar  ! [g/cm^3] oc
      dens_comp(4) = emis_dens_dust  ! [g/cm^3] dust 
      dens_comp(5) = emis_dens_seas  ! [g/cm^3] sea salt
      dens_comp(6) = rho_nh42so4     ! [g/cm^3] nitrate
      dens_comp(7) = rho_nh42so4     ! [g/cm^3] ammonium
      dens_comp(8) = rho_h2o         ! [g/cm^3] water 
      do i=1, 8
        recip_dens_comp(i) = 1.0d+00 / dens_comp(i)   ! [cm^3/g] sulfate 
        ! write(*,'(i4,2f10.4)') i, dens_comp(i), recip_dens_comp(i)
      enddo

      !---------------------------------------------------------------------------------------------------------------------
      ! if doing comparison with the discrete pdf model, set all mode particle densities to the same value.
      !---------------------------------------------------------------------------------------------------------------------
      if( discrete_eval_option ) then
        if( write_log ) write(aunit1,'(/a,f7.3/)') 'setting particle densities for all modes to (g/cm^3) ', densp
        denspi(:) = densp
        dens_comp(:) = densp 
        recip_dens_comp(:) = 1.0d+00 / densp 
      endif

90000 format(i12,a12,6f12.6)
      return
      END SUBROUTINE Setup_Dp0


      subroutine Setup_Coag_Tensors
!-------------------------------------------------------------------------------------------------------------------
!     routine to define the g_ikl,q, the d_ikl, and the d_ij.
!
!     all elements giklq(i,k,l,q), dikl(i,k,l), and dij(i,j) were checked through printouts available below. 
!
!     nm_spc_name(i,:) contains the names of the mass species defined for mode i.
!-------------------------------------------------------------------------------------------------------------------
      implicit none
      integer :: i,j,k,l,q,qq, ntot, n
      logical, parameter :: write_tensors = .false.

      if ( write_tensors ) then
        write(aunit1,'(/a/)') 'CITABLE'
        do i=1, nmodes
          write(aunit1,90000) citable(1:nmodes,i)
        enddo
        write(aunit1,'(a)') '  '
      endif

      giklq(:,:,:,:) = 0
      dikl(:,:,:) = 0
      dij(:,:) = 0

      !-------------------------------------------------------------------------------------------------------------
      ! The tensors g_ikl,q and d_ikl are symmetric in K and L.
      !
      ! GIKLQ is unity if coagulation of modes K and L produce mass of species Q
      !       in mode I, and zero otherwise.
      !
      ! DIKL is unity if coagulation of modes K and L produce particles 
      !      in mode I, and zero otherwise. 
      !      Neither mode K nor mode L can be mode I for a nonzero DIKL:
      !      all three modes I, K, L must be different modes. 
      !-------------------------------------------------------------------------------------------------------------
      do i=1, nmodes
      do k=1, nmodes
      do l=k+1, nmodes                              ! mode l is the same as mode k.
        if ( citable(k,l) .eq. mode_name(i) ) then  ! modes k and l produce mode i
          ! write(36,*)'mode_name(i) = ', mode_name(i)
          if ( i .ne. k  .and. i .ne. l ) then      ! omit intramodal coagulation
            dikl(i,k,l) = 1
            dikl(i,l,k) = 1
          endif
          do q=1, nm(i)                             ! loop over all mass species in mode i
            do qq=1, nmass_spcs                     ! loop over all principal mass species 
              !-----------------------------------------------------------------------------------------------------
              ! compare the name of mass species q in mode i with that of mass species qq in mode k (or l).
              ! the inner loop is over all principal mass species since all species must be checked for 
              !   mode k (or l) for a potential match with species q in mode i.
              !-----------------------------------------------------------------------------------------------------
              if( nm_spc_name(k,qq) .eq. nm_spc_name(i,q) ) then   ! mode k contains q
                if( i .ne. k ) giklq(i,k,l,q) = 1   ! i and k must be different modes
                if( i .ne. k ) giklq(i,l,k,q) = 1   ! i and k must be different modes
              endif
              if( nm_spc_name(l,qq) .eq. nm_spc_name(i,q) ) then   ! mode l contains q
                if( i .ne. l ) giklq(i,k,l,q) = 1   ! i and l must be different modes
                if( i .ne. l ) giklq(i,l,k,q) = 1   ! i and l must be different modes
              endif
            enddo
          enddo
        endif
      enddo
      enddo
      enddo

      if (allocated(dikl_control)) then
        deallocate(dikl_control)
      end if
      allocate(dikl_control(count(dikl /= 0)))

      call initializediklcontrol(dikl_control, dikl)

      if (allocated(giklq_control)) then
        do i = 1, nweights
          deallocate(giklq_control(i)%k)
          deallocate(giklq_control(i)%l)
          deallocate(giklq_control(i)%qq)
        end do
      else
        allocate(giklq_control(nweights))
      end if

      call InitializegIklqControl(giklq_control, giklq)

      !-------------------------------------------------------------------------------------------------------------
      ! the tensor d_ij is not symmetric in i,j.
      !
      ! dij(i,j) is unity if coagulation of mode i with mode j results
      !   in the removal of particles from mode i, and zero otherwise. 
      !-------------------------------------------------------------------------------------------------------------
      do i=1, nmodes
      do j=1, nmodes
        do k=1, nmodes                               ! find the product mode of the i-j coagulation.
          if( i .eq. j ) cycle                       ! omit intramodal interactions: --> i .ne. j .
          if( citable(i,j) .eq. mode_name(k) ) then  ! i-particles and j-particles are lost; k-particles are formed.
            if( i .ne. k ) dij(i,j) = 1              ! the k-particles are not i-particles (but may be j-particles),
          endif                                      !   so i-particles are lost by this i-j interaction.
        enddo
      enddo
      enddo
      xdij = dij

      if( .not. write_tensors ) return

      !-------------------------------------------------------------------------
      ! write the g_ikl,q.
      !-------------------------------------------------------------------------
      do i=1, nmodes
        write(aunit1,'(/2A)') 'g_iklq for MODE ', MODE_NAME(I)
        DO Q=1, NM(I)
          WRITE(AUNIT1,'(/A,I3,3X,3A5/)') 'Q, NM_SPC_NAME(I,Q), MODE', &
                              q, nm_spc_name(i,q), '-->', mode_name(i)
          if ( sum( giklq(i,1:nmodes,1:nmodes,q) ) .eq. 0 ) cycle
          write(aunit1,'(5x,16a5)') mode_name(1:nmodes)
          do k=1, nmodes
            write(aunit1,'(a5,16i5)') mode_name(k),giklq(i,k,1:nmodes,q)
          enddo
        enddo
      enddo

      !-------------------------------------------------------------------------
      ! Write the d_ikl.
      !-------------------------------------------------------------------------
      do i=1, nmodes
        write(aunit1,'(/2A)') 'd_ikl for MODE ', mode_name(i)
        if ( sum( dikl(i,1:nmodes,1:nmodes) ) .eq. 0 ) cycle
        write(aunit1,'(5x,16a5)') mode_name(1:nmodes)
        do k=1, nmodes
          write(aunit1,'(a5,16i5)') mode_name(k),dikl(i,k,1:nmodes)
        enddo
      enddo

      !-------------------------------------------------------------------------
      ! write the d_ij.
      !-------------------------------------------------------------------------
      write(aunit1,'(/2a)') 'd_ij'
      write(aunit1,'(5x,16a5)') mode_name(1:nmodes)
      do i=1, nmodes
        write(aunit1,'(a5,16i5)') mode_name(i),dij(i,1:nmodes)
      enddo

90000 format(14a4)
      return

      contains

      subroutine InitializeDIklControl(control, mask)
      type (dikl_type) :: control(:)
      integer, intent(in) :: mask(:,:,:)
      integer :: i, k, l, n

      n = 0
      do k = 1, nweights
        do l = k+1, nweights
          do i = 1, nweights
            if (mask(i,k,l) /= 0) then
              n = n + 1
              control(n)%i = i
              control(n)%k = k
              control(n)%l = l
            end if          
          end do
        end do
      end do
      NDIKL = n
      end subroutine initializeDiklControl

      subroutine InitializeGiklqControl(control, mask)
      type (GIKLQ_type) :: control(:)
      integer, intent(in) :: mask(:,:,:,:)

      integer :: i, q, k, l, n, nTotal

      do i = 1, nweights
        ! 1) count contributing cases for mode i
        n = 0
        do q = 1, nm(i)
          do k = 1, nmodes
            do l = k+1, nmodes
              if (mask(i,k,l,q) /= 0) then
                if (i /= l) then
                  n = n + 1
                end if
                if (i /= k) then
                  n = n + 1
                end if
              end if
            end do
          end do
        end do
        nTotal = n
        ! 2) allocate nTotal entries
        control(i)%n = nTotal
        allocate(control(i)%k(nTotal))
        allocate(control(i)%l(nTotal))
        allocate(control(i)%qq(nTotal))

        ! 3) repeat sweep, but now assign k,l,qq
        n = 0
        do q = 1, nm(i)
          do k = 1, nmodes
            do l = k+1, nmodes
              if (mask(i,k,l,q) /= 0) then
                if (i /= l) then
                  n = n + 1
                  control(i)%k(n) = k
                  control(i)%l(n) = l
                  qq = prod_index(i,q)
                  control(i)%qq(n) = qq
                end if
                if (i /= k) then
                  n = n + 1
                  control(i)%k(n) = l
                  control(i)%l(n) = k
                  qq = prod_index(i,q)
                  control(i)%qq(n) = qq
                end if
              end if
            end do
          end do
        end do
        control(i)%n = n
      end do

      end subroutine initializeGiklqControl


      END SUBROUTINE Setup_Coag_Tensors


      SUBROUTINE Setup_Aero_Mass_Map
!-------------------------------------------------------------------------------
!     Defines a map giving the AERO locations of the NM(I) masses for each mode.
!
!     MASS_MAP(I,Q) is the location in AERO(:) of the Qth mass in mode I.
!
!     PROD_INDEX(I,Q) is the location in array PIQ(I,Q) of chemical species
!       CHEM_SPC_NAME(Q) for mode (quadrature weight) I.
!-------------------------------------------------------------------------------
      implicit none
      integer :: i,q,j

      mass_map(:,:) = 0
      prod_index(:,:) = 0

      if( write_log ) write(aunit1,'(/a/)')'I,J,Q,MODE_NAME(I),AERO_SPCS(J),PROD_INDEX(I,Q),MASS_MAP(I,Q)'

      do i=1, nweights                                       ! loop over modes (quadrature points)
        q = 1
        do j=1, naerobox                                     ! loop over aero species
          if ( aero_spcs(j)(6:8) .ne. mode_name(i) ) cycle   ! this aero species is not for mode i.
          if ( aero_spcs(j)(1:4) .eq. 'MASS' ) then          ! this is a mass species for mode i. 
            if ( aero_spcs(j)(10:13) .eq. nm_spc_name(i,q) ) mass_map(i,q) = j  ! location of this species in aero
            if ( aero_spcs(j)(10:13) .eq. 'SULF' ) then      ! this mass species is sulfate.
              prod_index(i,q) = prod_index_sulf
            endif
            if ( aero_spcs(j)(10:13) .eq. 'BCAR' ) then      ! this mass species is bc.
              prod_index(i,q) = prod_index_bcar
            endif
            if ( aero_spcs(j)(10:13) .eq. 'OCAR' ) then      ! this mass species is oc.
              prod_index(i,q) = prod_index_ocar
            endif
            if ( aero_spcs(j)(10:13) .eq. 'DUST' ) then      ! this mass species is dust.
              prod_index(i,q) = prod_index_dust
            endif
            if ( aero_spcs(j)(10:13) .eq. 'SEAS' ) then      ! this mass species is sea salt.
              prod_index(i,q) = prod_index_seas
            endif
            if( write_log ) write(aunit1,90000)i,j,q,mode_name(i),aero_spcs(j),prod_index(i,q),mass_map(i,q)
            q = q + 1
            if (q .gt. nm(i) ) goto 10
          endif
        enddo
10      continue
      enddo

90000 format(3i4,a8,4x,a16,5x,2i4)
      return
      end subroutine Setup_Aero_Mass_Map

      SUBROUTINE Setup_Kci_dyn(pres,tk,m1)         
!-----------------------------------------------------------------------------------------------------------------------
!     Routine to calculate the coefficients that multiply the number
!     concentrations, or the number concentrations times the particle diameters,
!     to obtain the condensational sink for each mode or quadrature point.
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: m1
      real(8), INTENT(IN) :: tk(m1),pres(m1)
      
      integer :: i, l     ! indices
      real    :: sigma    ! see subr. atmosphere below.
      real    :: delta    ! see subr. atmosphere below.
      real    :: theta    ! see subr. atmosphere below.
      real(8) :: p        ! ambient pressure [pa]
      real(8) :: t        ! ambient temperature [k]
      real(8) :: d        ! molecular diffusivity of h2so4 in air [m^2/s]
      real(8) :: c        ! mean molecular speed of h2so4 [m/s]
      real(8) :: la       ! mean free path in air [m]
                          ! 6.6328d-08 is the sea level value given in table i.2.8
                          ! on p.10 of u.s. standard atmosphere 1962
      real(8) :: lh       ! mean free path of h2so4 in air [m]
      real(8) :: beta     ! transition regime correction to the condensational flux [1]
      real(8) :: kn       ! knudsen number for h2so4 in air [1]
      real(8) :: thetai   ! monodispersity correction factor [1] (okuyama et al. 1988)
      real(8) :: scale_dp ! scale factor for defined table of ambient particle diameters [1]

      real(8), parameter :: alpha    = 0.86d+00      ! mass accommodation coefficient [1] (hansen, 2005)
      real(8), parameter :: d_stdatm = 1.2250d+00    ! sea-level std. density [kg/m^3]
      real(8), parameter :: p_stdatm = 101325.0d+00  ! sea-level std. pressure [pa]
      real(8), parameter :: t_stdatm = 288.15d+00    ! sea-level std. temperature [k]
      real(8), parameter :: p0       = 101325.0d+00  ! reference pressure    [pa] for d
      real(8), parameter :: t0       = 273.16d+00    ! reference temperature [k]  for d
      real(8), parameter :: d0       = 9.36d-06      ! diffusivity of h2so4 in air [m^2/s]
                                                     ! calculated using eqn 11-4.4 of reid 
                                                     ! et al.(1987) at 273.16 k and 101325 pa

      if( write_log ) write(aunit1,90002)'I','L','ZHEIGHT','P','T','D','C','LA','LH', &
                                        'DP0','KN','THETAI','BETA'
      kci_coef_dp          (:,:) = 0.0d+00   ! mass accommodation coefficient is arbitrary
      kci_coef_dp_aeq1     (:,:) = 0.0d+00   ! mass accommodation coefficient is unity
      kci_dp_condtable     (:,:) = 0.0d+00   ! mass accommodation coefficient is arbitrary
      kci_dp_condtable_aeq1(:,:) = 0.0d+00   ! mass accommodation coefficient is unity
      !----------------------------------------------------------------------------------------------------------------
      ! setup table of ambient particle diameters.
      !----------------------------------------------------------------------------------------------------------------
      scale_dp = ( dp_condtable_max / dp_condtable_min )**(1.0d+00/real(n_dp_condtable-1))
      xln_scale_dp = log( scale_dp )
      do i=1, n_dp_condtable
        dp_condtable(i) = dp_condtable_min * scale_dp**(i-1)              ! [m]
      enddo

      do l=1, m1
        !call atmosphere( real( zheight(l) ), sigma, delta, theta )
!        p = delta * p_stdatm                                           ! [pa}
!        t = theta * t_stdatm                                           ! [k]
        p = pres(l)                                                    ! [pa}
        t = tk(l)                                                      ! [k]
        d = d0 * ( p0 / p ) * ( t / t0 )**1.75                         ! [m^2/s]
        !PRINT *, 'LFR-l,p,t,d:',l,p,t,d
        diffcoef_m2s(l) = d                                            ! [m^2/s]
        c = sqrt( 8.0d+00 * rgas_si * t / ( pi * mw_h2so4*1.0d-03 ) )  ! [m/s]
        la = 6.6332d-08 * ( p_stdatm / p ) * ( t / t_stdatm )          ! [m]
        lh = 3.0d+00 * d / c                                           ! [m]
        do i=1, nweights
          thetai = exp( - ( log(sig0(i)) )**2 )                        ! [1]
          theta_poly(i) = thetai                                       ! [1] polydispersity adjustment factor
          !------------------------------------------------------------------------------------------------------------
          ! for the condensation sink for general use, the mean free path is
          ! that for the condensing vapor (h2so4), and the mass accommodation coefficient is adjustable. 
          !------------------------------------------------------------------------------------------------------------
          kn = 2.0d+00 * lh / dp0(i)                                   ! lh and dp0 in [m]
          beta = ( 1.0d+00 + kn )  &                                   ! [1]
               / ( 1.0d+00 + 0.377d+00*kn + 1.33d+00*kn*(1.0d+00 + kn)/alpha ) 
          kci_coef_dp(i,l)      = 2.0d+00 * pi * thetai * d * beta * dp0(i) ! [m^3/s] may be updated in subr. matrix
          !------------------------------------------------------------------------------------------------------------
          ! for the condensation sink for use in kerminen and kulmala (2002), the mean free path is
          ! that for air, and the mass accommodation coefficient is set to unity. 
          !------------------------------------------------------------------------------------------------------------
          kn = 2.0d+00 * la / dp0(i)                                   ! la and dp0 in [m]
          beta = ( 1.0d+00 + kn )  &                                   ! [1]
               / ( 1.0d+00 + 0.377d+00*kn + 1.33d+00*kn*(1.0d+00 + kn)/1.0d+00 ) 
          kci_coef_dp_aeq1(i,l) = 2.0d+00 * pi * thetai * d * beta * dp0(i) ! [m^3/s] may be updated in subr. matrix
          ! if( write_log ) write(aunit1,90000)i,l,zheight(l),p,t,d,c,la,lh,dp0(i)*1.0d+06,kn,thetai,beta
        enddo
        do i=1, n_dp_condtable
          !------------------------------------------------------------------------------------------------------------
          ! for the condensation sink for general use, the mean free path is
          ! that for the condensing vapor (h2so4), and the mass accommodation coefficient is adjustable. 
          ! the thetai factor is included later in aero_matrix.f.
          !------------------------------------------------------------------------------------------------------------
          kn = 2.0d+00 * lh / dp_condtable(i)                          ! lh and dp in [m]
          beta = ( 1.0d+00 + kn )   &                                  ! [1]
               / ( 1.0d+00 + 0.377d+00*kn + 1.33d+00*kn*(1.0d+00 + kn)/alpha ) 
          kci_dp_condtable(i,l)      = 2.0d+00 * pi * d * beta * dp_condtable(i)  ! [m^3/s]
          !------------------------------------------------------------------------------------------------------------
          ! for the condensation sink for use in kerminen and kulmala (2002), the mean free path is
          ! that for air, and the mass accommodation coefficient is set to unity. 
          ! the thetai factor is included later in aero_matrix.f.
          !------------------------------------------------------------------------------------------------------------
          kn = 2.0d+00 * la / dp_condtable(i)                          ! la and dp in [m]
          beta = ( 1.0d+00 + kn )  &                                   ! [1]
               / ( 1.0d+00 + 0.377d+00*kn + 1.33d+00*kn*(1.0d+00 + kn)/1.0d+00 ) 
          kci_dp_condtable_aeq1(i,l) = 2.0d+00 * pi * d * beta * dp_condtable(i)  ! [m^3/s]
          ! write(aunit1,90000) i, l, zheight(l), p, t, d, c, la, lh, dp_condtable(i)*1.0d+06, kn, thetai, beta
        enddo
      enddo
      if( write_log ) then
        write(aunit1,'(/A/)')'  I  L     KCI_COEF_DP[m^3/s] KCI_COEF_DP_AEQ1[m^3/s]'    
        do l=1, nlays
          do i=1, nweights
            write(aunit1,90001) i, l, kci_coef_dp(i,l), kci_coef_dp_aeq1(i,l)
          enddo
        enddo
      endif
      do i=1, nmodes
        if ( icond(i) .eq. 0 ) then 
          kci_coef_dp     (i,:) = 0.0d+00
          kci_coef_dp_aeq1(i,:) = 0.0d+00
          theta_poly(i)         = 0.0d+00
          if( write_log ) write(aunit1,'(/2A/)') 'Condensational growth turned off for mode ', mode_name(i)
        ENDIF
      ENDDO

90000 FORMAT( 2I3,F8.4,F9.1,F7.2,D10.3,F6.1,2D10.3,4F8.4)
90002 FORMAT(/2A3,A8,  A9,  A7,  A10,  A6,  2A10,  4A8 /)
90001 FORMAT(2I3,4D20.4)
      RETURN
      END SUBROUTINE SETUP_KCI_dyn 

SUBROUTINE setup_kci         
!-----------------------------------------------------------------------------------------------------------------------
!     Routine to calculate the coefficients that multiply the number
!     concentrations, or the number concentrations times the particle diameters,
!     to obtain the condensational sink for each mode or quadrature point.
!-----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i, l     ! indices
      REAL    :: sigma    ! See subr. ATMOSPHERE below.
      REAL    :: delta    ! See subr. ATMOSPHERE below.
      REAL    :: theta    ! See subr. ATMOSPHERE below.
      REAL(8) :: p        ! ambient pressure [Pa]
      REAL(8) :: t        ! ambient temperature [K]
      REAL(8) :: d        ! molecular diffusivity of H2SO4 in air [m^2/s]
      REAL(8) :: c        ! mean molecular speed of H2SO4 [m/s]
      REAL(8) :: la       ! mean free path in air [m]
                          ! 6.6328D-08 is the sea level value given in Table I.2.8
                          ! on p.10 of U.S. Standard Atmosphere 1962
      REAL(8) :: lh       ! mean free path of H2SO4 in air [m]
      REAL(8) :: beta     ! transition regime correction to the condensational flux [1]
      REAL(8) :: kn       ! Knudsen number for H2SO4 in air [1]
      REAL(8) :: thetai   ! monodispersity correction factor [1] (Okuyama et al. 1988)
      REAL(8) :: scale_dp ! scale factor for defined table of ambient particle diameters [1]

      REAL(8), PARAMETER :: alpha    = 0.86D+00      ! mass accommodation coefficient [1] (Hansen, 2005)
      REAL(8), PARAMETER :: d_stdatm = 1.2250D+00    ! sea-level std. density [kg/m^3]
      REAL(8), PARAMETER :: p_stdatm = 101325.0D+00  ! sea-level std. pressure [Pa]
      REAL(8), PARAMETER :: t_stdatm = 288.15D+00    ! sea-level std. temperature [K]
      REAL(8), PARAMETER :: p0       = 101325.0D+00  ! reference pressure    [Pa] for D
      REAL(8), PARAMETER :: t0       = 273.16D+00    ! reference temperature [K]  for D
      REAL(8), PARAMETER :: d0       = 9.36D-06      ! diffusivity of H2SO4 in air [m^2/s]
                                                     ! calculated using Eqn 11-4.4 of Reid 
                                                     ! et al.(1987) at 273.16 K and 101325 Pa

      IF( WRITE_LOG ) WRITE(AUNIT1,90002)'I','L','ZHEIGHT','P','T','D','C','LA','LH', &
                                        'DP0','KN','THETAI','BETA'
      kci_coef_dp          (:,:) = 0.0d+00   ! mass accommodation coefficient is arbitrary
      kci_coef_dp_aeq1     (:,:) = 0.0d+00   ! mass accommodation coefficient is unity
      kci_dp_condtable     (:,:) = 0.0d+00   ! mass accommodation coefficient is arbitrary
      kci_dp_condtable_aeq1(:,:) = 0.0D+00   ! mass accommodation coefficient is unity
      !----------------------------------------------------------------------------------------------------------------
      ! Setup table of ambient particle diameters.
      !----------------------------------------------------------------------------------------------------------------
      scale_dp = ( dp_condtable_max / dp_condtable_min )**(1.0d+00/real(n_dp_condtable-1))
      xln_scale_dp = log( scale_dp )
      do i=1, n_dp_condtable
        dp_condtable(i) = dp_condtable_min * scale_dp**(i-1)              ! [m]
      enddo

      do l=1, nlays
        call atmosphere( real( zheight(l) ), sigma, delta, theta )
        p = delta * p_stdatm                                           ! [pa}
        t = theta * t_stdatm                                           ! [k]
        d = d0 * ( p0 / p ) * ( t / t0 )**1.75                         ! [m^2/s]
        diffcoef_m2s(l) = d                                            ! [m^2/s]
        c = sqrt( 8.0d+00 * rgas_si * t / ( pi * mw_h2so4*1.0d-03 ) )  ! [m/s]
        la = 6.6332d-08 * ( p_stdatm / p ) * ( t / t_stdatm )          ! [m]
        lh = 3.0d+00 * d / c                                           ! [m]
        do i=1, nweights
          thetai = exp( - ( log(sig0(i)) )**2 )                        ! [1]
          theta_poly(i) = thetai                                       ! [1] polydispersity adjustment factor
          !------------------------------------------------------------------------------------------------------------
          ! For the condensation sink for general use, the mean free path is
          ! that for the condensing vapor (H2SO4), and the mass accommodation coefficient is adjustable. 
          !------------------------------------------------------------------------------------------------------------
          kn = 2.0d+00 * lh / dp0(i)                                   ! LH and DP0 in [m]
          BETA = ( 1.0D+00 + KN ) &                                      ! [1]
                / ( 1.0d+00 + 0.377d+00*kn + 1.33d+00*kn*(1.0d+00 + kn)/alpha ) 
          kci_coef_dp(i,l)      = 2.0d+00 * pi * thetai * d * beta * dp0(i) ! [m^3/s] may be updated in subr. matrix
          !------------------------------------------------------------------------------------------------------------
          ! For the condensation sink for use in Kerminen and Kulmala (2002), the mean free path is
          ! that for air, and the mass accommodation coefficient is set to unity. 
          !------------------------------------------------------------------------------------------------------------
          kn = 2.0d+00 * la / dp0(i)                                   ! LA and DP0 in [m]
          beta = ( 1.0d+00 + kn )       &                               ! [1]
                  / ( 1.0d+00 + 0.377d+00*kn + 1.33d+00*kn*(1.0d+00 + kn)/1.0D+00 ) 
          kci_coef_dp_aeq1(i,l) = 2.0d+00 * pi * thetai * d * beta * dp0(i) ! [m^3/s] may be updated in subr. matrix
          ! IF( WRITE_LOG ) WRITE(AUNIT1,90000)I,L,ZHEIGHT(L),P,T,D,C,LA,LH,DP0(I)*1.0D+06,KN,THETAI,BETA
        ENDDO
        DO I=1, n_dp_condtable
          !------------------------------------------------------------------------------------------------------------
          ! For the condensation sink for general use, the mean free path is
          ! that for the condensing vapor (H2SO4), and the mass accommodation coefficient is adjustable. 
          ! The THETAI factor is included later in aero_matrix.f.
          !------------------------------------------------------------------------------------------------------------
          kn = 2.0d+00 * lh / dp_condtable(i)                          ! LH and DP in [m]
          beta = ( 1.0d+00 + kn )     &                                 ! [1]
                 / ( 1.0d+00 + 0.377d+00*kn + 1.33d+00*kn*(1.0d+00 + kn)/alpha ) 
          kci_dp_condtable(i,l)      = 2.0d+00 * pi * d * beta * dp_condtable(i)  ! [m^3/s]
          !------------------------------------------------------------------------------------------------------------
          ! For the condensation sink for use in Kerminen and Kulmala (2002), the mean free path is
          ! that for air, and the mass accommodation coefficient is set to unity. 
          ! The THETAI factor is included later in aero_matrix.f.
          !------------------------------------------------------------------------------------------------------------
          kn = 2.0d+00 * la / dp_condtable(i)                          ! LA and DP in [m]
          beta = ( 1.0d+00 + kn )        &                              ! [1]
              / ( 1.0d+00 + 0.377d+00*kn + 1.33d+00*kn*(1.0d+00 + kn)/1.0d+00 ) 
          kci_dp_condtable_aeq1(i,l) = 2.0d+00 * pi * d * beta * dp_condtable(i)  ! [m^3/s]
          ! WRITE(AUNIT1,90000) I, L, ZHEIGHT(L), P, T, D, C, LA, LH, DP_CONDTABLE(I)*1.0D+06, KN, THETAI, BETA
        ENDDO
      ENDDO
      IF( WRITE_LOG ) THEN
        WRITE(AUNIT1,'(/A/)')'  I  L     KCI_COEF_DP[m^3/s] KCI_COEF_DP_AEQ1[m^3/s]'    
        DO L=1, NLAYS
          DO I=1, NWEIGHTS
            WRITE(AUNIT1,90001) I, L, KCI_COEF_DP(I,L), KCI_COEF_DP_AEQ1(I,L)
          ENDDO
        ENDDO
      ENDIF
      do i=1, nmodes
        if ( icond(i) .eq. 0 ) then 
          kci_coef_dp     (i,:) = 0.0d+00
          kci_coef_dp_aeq1(i,:) = 0.0d+00
          theta_poly(i)         = 0.0d+00
          if( write_log ) write(aunit1,'(/2a/)') 'condensational growth turned off for mode ', mode_name(i)
        endif
      enddo

90000 FORMAT( 2I3,F8.4,F9.1,F7.2,D10.3,F6.1,2D10.3,4F8.4)
90002 FORMAT(/2A3,A8,  A9,  A7,  A10,  A6,  2A10,  4A8 /)
90001 FORMAT(2I3,4D20.4)
      
      END SUBROUTINE setup_kci 


      SUBROUTINE atmosphere( alt, sigma, delta, theta )
!----------------------------------------------------------------------------
! PURPOSE - Compute the properties of the 1976 standard atmosphere to 86 km.
! AUTHOR  - Ralph Carmichael, Public Domain Aeronautical Software
! Reformatted for fixed-form Fortran 90 by D. Wright, 1-9-06.                 
! NOTE - If ALT > 86, the values returned will not be correct, but they will
!   not be too far removed from the correct values for density.
!   The reference document does not use the terms pressure and temperature
!   above 86 km.
!----------------------------------------------------------------------------
      IMPLICIT NONE
!============================================================================
!     A R G U M E N T S                                                     |
!============================================================================
      REAL,INTENT(IN)::  alt    ! geometric ALTitude, km.
      REAL,INTENT(OUT):: sigma  ! density/sea-level standard density
      REAL,INTENT(OUT):: delta  ! pressure/sea-level standard pressure
      REAL,INTENT(OUT):: theta  ! temperature/sea-level standard temperature
!============================================================================
!     L O C A L   C O N S T A N T S                                         |
!============================================================================
      REAL,PARAMETER:: rearth = 6369.           ! radius of the Earth (km)
      REAL,PARAMETER:: gmr = 34.163195                 ! hydrostatic constant
      INTEGER,PARAMETER:: ntab=8   ! number of entries in the defining tables
!============================================================================
!     L O C A L   V A R I A B L E S                                         |
!============================================================================
      INTEGER:: i,j,k                                              ! counters
      REAL:: h                                   ! geopotential ALTitude (km)
      REAL:: tgrad, tbase  ! temperature gradient and base temp of this layer
      REAL:: tlocal                                       ! local temperature
      REAL:: deltah                         ! height above base of this layer
!============================================================================
!     L O C A L   A R R A Y S   ( 1 9 7 6   S T D.  A T M O S P H E R E )   |
!============================================================================
      REAL,DIMENSION(ntab),PARAMETER:: htab= (/0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852/)
      REAL,DIMENSION(ntab),PARAMETER:: ttab= (/288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946/)
      REAL,DIMENSION(ntab),PARAMETER:: ptab= (/1.0, 2.233611E-1, 5.403295E-2, 8.5666784E-3, 1.0945601E-3, &
                                               6.6063531E-4, 3.9046834E-5, 3.68501E-6/)
      REAL,DIMENSION(ntab),PARAMETER:: gtab= (/-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0/)

      h=alt*rearth/(alt+rearth)  ! convert geometric to geopotential altitude

      i=1
      j=ntab                                   ! setting up for binary search
      DO
        k=(i+j)/2                                          ! integer division
        IF (h < htab(k)) THEN
          j=k
        ELSE
          i=k
        END IF
        IF (j <= i+1) EXIT
      END DO

      tgrad=gtab(i)                                 ! I will be in 1...NTAB-1
      tbase=ttab(i)
      deltah=h-htab(i)
      tlocal=tbase+tgrad*deltah
      theta=tlocal/ttab(1)                                ! temperature ratio

      IF (tgrad == 0.0) THEN                                 ! pressure ratio
        delta=ptab(i)*exp(-gmr*deltah/tbase)
      ELSE
        delta=ptab(i)*(tbase/tlocal)**(gmr/tgrad)
      END IF

      sigma=delta/theta                                       ! density ratio

      END SUBROUTINE atmosphere 



      end module Aero_Setup
 
