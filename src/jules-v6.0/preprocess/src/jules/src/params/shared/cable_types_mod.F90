MODULE cable_types_mod

USE max_dimensions,             ONLY: npft_max, nnvg_max
USE cable_other_constants_mod,  ONLY: r_2, nsl, nscs, nvcs, nrb

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Defines variable types and variables for CABLE standalone runs.
!   Based on cable_def_types_mod.F90 from the CABLE trunk.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

REAL, POINTER :: land_alb(:,:)  ! Mean land albedo
INTEGER :: mp          ! # total no of patches/tiles

! Indices for special types
INTEGER, PARAMETER ::                                                         &
  perm_ice_veg = 17,    & ! permanent ice vegetation type index
  perm_ice_soil = 9,    & ! permanent ice soil type index
  non_ice_soil = 2        ! non-permanent ice soil type index

INTEGER, POINTER :: isnow_flg3l(:,:)

REAL, POINTER ::                                                              &
    snow_rho1l(:,:),     & ! mean snow density
    snow_age(:,:)

REAL, POINTER :: albsoil(:)

TYPE derived_rad_bands
  REAL, ALLOCATABLE ::                                                        &
   sw_down_dir(:,:),  & ! Surface downward SW direct radiation (W/m2).
   sw_down_dif(:,:),  & ! Surface downward SW diffuse radiation (W/m2).
   sw_down_vis(:,:),  & ! Surface downward VIS radiation (W/m2).
   sw_down_nir(:,:),  & ! Surface downward NIR radiation (W/m2).
   fbeam(:,:,:)         ! Surface downward SW radiation (W/m2).
END TYPE derived_rad_bands

LOGICAL, ALLOCATABLE :: l_tile_pts(:,:)
REAL, ALLOCATABLE    :: tile_index(:,:)


! Soil parameters:
TYPE soil_parameter_type
  INTEGER, POINTER ::                                                         &
    isoilm(:)     ! integer soil type

  REAL, POINTER ::                                                            &
    soilcol(:), & ! keep color for all patches/tiles
    albsoilf(:)   ! soil reflectance

  REAL, POINTER ::                                                            &
    albsoil(:,:)  ! soil reflectance (2nd dim. BP 21Oct2009)
END TYPE soil_parameter_type


! Soil and snow variables:
TYPE soil_snow_type
  INTEGER, POINTER ::                                                         &
         isflag(:)     ! 0 => no snow 1 => snow

  REAL, POINTER ::                                                            &
         osnowd(:),  & ! snow depth from previous time step
         snage(:),   & ! snow age
         snowd(:),   & ! snow depth (liquid water)
         ssdnn(:)      ! average snow density

  REAL, POINTER ::                                                            &
         tgg(:,:),        & ! soil temperature in K
         tggsn(:,:),      & ! snow temperature in K
         albsoilsn(:,:)     ! soil + snow reflectance

  REAL(r_2), POINTER ::                                                       &
         wb(:,:),      & ! volumetric soil moisture (solid+liq)
         wbice(:,:)      ! soil ice
END TYPE soil_snow_type


! Radiation variables:
TYPE radiation_type
  REAL, POINTER ::                                                            &
     extkb(:),   & ! beam radiation extinction coeff
     extkd(:)      ! diffuse radiation extinction coeff (-)

  REAL, POINTER ::                                                            &
     rhocdf(:,:),  & ! canopy diffuse reflectance (-)
     albedo(:,:),  & ! canopy+soil albedo
     reffdf(:,:),  & ! effective conopy diffuse reflectance
     reffbm(:,:),  & ! effective conopy beam reflectance
     extkbm(:,:),  & ! modified k beam(6.20)(for leaf scattering)
     extkdm(:,:),  & ! modified k diffuse(6.20)(for leaf scattering)
     fbeam(:,:),   & ! beam fraction
     cexpkbm(:,:), & ! canopy beam transmittance
     cexpkdm(:,:), & ! canopy diffuse transmittance
     rhocbm(:,:)     ! modified canopy beam reflectance(6.21)
END TYPE radiation_type


   ! Vegetation parameters:
TYPE veg_parameter_type
  INTEGER, POINTER ::                                                         &
     iveg(:)       ! vegetation type

  REAL, POINTER ::                                                            &
     hc(:),      & ! roughness height of canopy (veg - snow)
     vlai(:)       ! leaf area index

  REAL, POINTER ::                                                            &
     refl(:,:),                                                               &
     taul(:,:)
END TYPE veg_parameter_type


TYPE vegin_type
  REAL ::                                                                     &
       canst1(npft_max),                                                      &
       dleaf(npft_max),                                                       &
       length(npft_max),                                                      &
       width(npft_max),                                                       &
       vcmax(npft_max),                                                       &
       ejmax(npft_max),                                                       &
       hc(npft_max),                                                          &
       xfang(npft_max),                                                       &
       rp20(npft_max),                                                        &
       rpcoef(npft_max),                                                      &
       rs20(npft_max),                                                        &
       wai(npft_max),                                                         &
       rootbeta(npft_max),                                                    &
       shelrb(npft_max),                                                      &
       vegcf(npft_max),                                                       &
       frac4(npft_max),                                                       &
       xalbnir(npft_max),                                                     &
       extkn(npft_max),                                                       &
       tminvj(npft_max),                                                      &
       tmaxvj(npft_max),                                                      &
       vbeta(npft_max),                                                       &
       a1gs(npft_max),                                                        &
       d0gs(npft_max),                                                        &
       alpha(npft_max),                                                       &
       convex(npft_max),                                                      &
       cfrd(npft_max),                                                        &
       gswmin(npft_max),                                                      &
       conkc0(npft_max),                                                      &
       conko0(npft_max),                                                      &
       ekc(npft_max),                                                         &
       eko(npft_max),                                                         &
       g0(npft_max),                                                          &
       g1(npft_max),                                                          &
       zr(npft_max),                                                          &
       clitt(npft_max),                                                       &
       froot(nsl,npft_max),                                                   &
       csoil(nscs,npft_max),                                                  &
       ratecs(nscs,npft_max),                                                 &
       cplant(nvcs,npft_max),                                                 &
       ratecp(nvcs,npft_max),                                                 &
       refl(nrb,npft_max),                                                    &
       taul(nrb,npft_max)
END TYPE vegin_type

TYPE(vegin_type),  SAVE  :: vegin


TYPE soilin_type
  REAL ::                                                                     &
       silt(nnvg_max),                                                        &
       clay(nnvg_max),                                                        &
       sand(nnvg_max),                                                        &
       swilt(nnvg_max),                                                       &
       sfc(nnvg_max),                                                         &
       ssat(nnvg_max),                                                        &
       bch(nnvg_max),                                                         &
       hyds(nnvg_max),                                                        &
       sucs(nnvg_max),                                                        &
       rhosoil(nnvg_max),                                                     &
       css(nnvg_max)
END TYPE soilin_type

TYPE(soilin_type), SAVE :: soilin


! Meterological data:
TYPE met_type
  REAL, POINTER ::                                                            &
         tk(:),      & ! surface air temperature (oK)
         coszen(:)     ! cos(zenith angle of sun)

  REAL, POINTER ::                                                            &
         fsd(:,:)  ! downward short-wave radiation (W/m2)
END TYPE met_type


! Canopy/vegetation variables:
TYPE canopy_type
  REAL, POINTER ::                                                            &
    vlaiw(:)     ! lai adj for snow depth for calc of resistances
END TYPE canopy_type


TYPE(soil_parameter_type), SAVE :: soil       ! soil parameters
TYPE(soil_snow_type), SAVE      :: ssnow
TYPE(derived_rad_bands), SAVE   :: rad_bands
TYPE(radiation_type), SAVE      :: rad
TYPE(veg_parameter_type), SAVE  :: veg        ! vegetation parameters
TYPE(met_type), SAVE            :: met
TYPE(canopy_type), SAVE         :: canopy

REAL, ALLOCATABLE, SAVE  :: sw_down_diag(:,:)


!---CABLE runtime switches def in this type
TYPE cable_internal_switches
  LOGICAL :: um_radiation = .FALSE.
END TYPE cable_internal_switches

! instantiate internal switches
TYPE(cable_internal_switches), SAVE :: cable_runtime


! user switches turned on/off by the user thru namelists
! CABLE-2.0 user switches all in single namelist file cable.nml
! clean these up for new namelist(s) format
TYPE cable_user_switches

  CHARACTER(LEN=20) :: soil_struc    = "default" ! 'default' or 'sli'

END TYPE cable_user_switches


! instantiate internal switches
TYPE(cable_user_switches), SAVE :: cable_user

LOGICAL :: calcsoilalbedo = .FALSE.

INTEGER, SAVE :: kwidth_gl

END MODULE cable_types_mod
