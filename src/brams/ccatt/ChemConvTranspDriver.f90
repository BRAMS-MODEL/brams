MODULE ChemConvTranspDriver


  USE mem_grid, ONLY:           &
       dtlt,                    & ! (IN) REAL,
                                  !      delta t

       ngrid                      ! (IN) INTEGER,
                                  !      current grid
 
  USE Phys_const, ONLY:         &
       g                          ! (IN) REAL,
                                  !      PARAMETER, g=9.80

  USE mem_basic, ONLY:          &
       basic_g                    ! (IN) basic_g(ngrid)%dn0(mzp,mxp,myp)
                                  !      reference state density [kg/m^3]

  USE mem_grell_param, ONLY:    &

       maxiens,                 & ! (IN) INTEGER,

       mgmzp,                   & ! (IN) INTEGER,
                                  !      max(mmzp(i)), i=1,2,..,ngrids_cp
     
       mgmxp,                   & ! (IN) INTEGER,
                                  !      max(mmxp(i)), i=1,2,..,ngrids_cp

       mgmyp                      ! (IN) INTEGER,
                                  !      max(mmyp(i)), i=1,2,..,ngrids_cp
  
  USE mem_scratch1_grell, ONLY: &
       ierr4d,                  & ! (IN) INTEGER (mgmxp,mgmyp,maxiens,ngrids_cp)
                                  !      ierr4d = 0, there is no convection
                                  !      ierr4d = 1, convection

       jmin4d,                  & ! (IN) INTEGER (mgmxp,mgmyp,maxiens,ngrids_cp)
                                  !      ETL

       kdet4d,                  & ! (IN) INTEGER (mgmxp,mgmyp,maxiens,ngrids_cp)
                                  !      detrainemnt level 

       k224d,                   & ! (IN) INTEGER (mgmxp,mgmyp,maxiens,ngrids_cp)
                                  !

       kbcon4d,                 & ! (IN) INTEGER (mgmxp,mgmyp,maxiens,ngrids_cp)
                                  !      cloud base

       ktop4d,                  & ! (IN) INTEGER (mgmxp,mgmyp,maxiens,ngrids_cp)
                                  !      cloud top
                                  
       kpbl4d,                  & ! (IN) INTEGER (mgmxp,mgmyp,maxiens,ngrids_cp)
                                  !      PBL height

       kstabi4d,                & ! (IN) INTEGER (mgmxp,mgmyp,maxiens,ngrids_cp)
                                  !      cloud base

       kstabm4d,                & ! (IN) INTEGER (mgmxp,mgmyp,maxiens,ngrids_cp)
                                  !      cloud top

       xmb4d,                   & ! (IN) REAL (mgmxp,mgmyp,maxiens,ngrids_cp)
                                  !      updraft mass flux

       edt4d,                   & ! (IN) REAL (mgmxp,mgmyp,maxiens,ngrids_cp)
                                  !      MFLX_DOWN/MFLX_UP

       enup5d,                  & ! (IN) REAL (mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp)
                                  !      entrain up   rate

       endn5d,                  & ! (IN) REAL (mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp)
                                  !      entrain down rate

       deup5d,                  & ! (IN) REAL (mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp)
                                  !      detrain up   rate

       dedn5d,                  & ! (IN) REAL (mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp)
                                  !      detrain down rate

       zup5d,                   & ! (IN) REAL (mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp)
                                  !      norm mass flux up

       zdn5d,                   & ! (IN) REAL (mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp)
                                  !      norm mass flux down

       zcup5d,                  & ! (IN) REAL (mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp)
                                  !      z level

       pcup5d,                  & ! (IN) REAL (mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp)
                                  !      p level

       clwup5d,                 & ! (IN) REAL (mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp)

       tup5d                      ! (IN) REAL (mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp)

  USE mem_chem1, ONLY:            &
       chemistry,                 & ! (IN)
       chem1_g                      ! (INOUT)
  
  USE mem_aer1, ONLY:             &
       aer1_g                       ! (INOUT)

  USE mod_chem_conv_transp, ONLY: &
       trans_conv_mflx              ! Subroutine


  IMPLICIT NONE


  PRIVATE

  PUBLIC :: transconv_driver ! Subroutine



CONTAINS



  SUBROUTINE transconv_driver(mzp,mxp,myp,ia,iz,ja,jz,i0,j0,iens,stcum)

    INTEGER , INTENT(IN)    :: mzp
    INTEGER , INTENT(IN)    :: mxp
    INTEGER , INTENT(IN)    :: myp
    INTEGER , INTENT(IN)    :: ia
    INTEGER , INTENT(IN)    :: iz
    INTEGER , INTENT(IN)    :: ja
    INTEGER , INTENT(IN)    :: jz
    INTEGER , INTENT(IN)    :: i0
    INTEGER , INTENT(IN)    :: j0
    INTEGER , INTENT(IN)    :: iens
    REAL    , INTENT(INOUT) :: stcum(mzp,mxp,myp)
    
    
    if(CHEMISTRY < 0) return

    CALL trans_conv_mflx(mzp,mxp,myp,ia,iz,ja,jz,i0,j0,iens,stcum, &
                         g,maxiens,mgmzp,mgmxp,mgmyp,              &
                         dtlt,chemistry,                           &
                         ierr4d    (:,:,:,ngrid),                  &
                         jmin4d    (:,:,:,ngrid),                  &
                         kdet4d    (:,:,:,ngrid),                  &
                         k224d     (:,:,:,ngrid),                  &
                         kbcon4d   (:,:,:,ngrid),                  &
                         ktop4d    (:,:,:,ngrid),                  &
                         kpbl4d    (:,:,:,ngrid),                  &
                         kstabi4d  (:,:,:,ngrid),                  &
                         kstabm4d  (:,:,:,ngrid),                  &
                         xmb4d     (:,:,:,ngrid),                  &
                         edt4d     (:,:,:,ngrid),                  &
                         enup5d  (:,:,:,:,ngrid),                  &
                         endn5d  (:,:,:,:,ngrid),                  &
                         deup5d  (:,:,:,:,ngrid),                  &
                         dedn5d  (:,:,:,:,ngrid),                  &
                         zup5d   (:,:,:,:,ngrid),                  &
                         zdn5d   (:,:,:,:,ngrid),                  &
                         zcup5d  (:,:,:,:,ngrid),                  &
                         pcup5d  (:,:,:,:,ngrid),                  &
                         clwup5d (:,:,:,:,ngrid),                  &
                         tup5d   (:,:,:,:,ngrid),                  &
                         basic_g (ngrid)%dn0,                      &
                         chem1_g   (:,ngrid),                      &
                         aer1_g  (:,:,ngrid))

  END SUBROUTINE transconv_driver



END MODULE ChemConvTranspDriver
