! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  ---------------------------------------------------------------------
!  Sructure containing surface exchange diagnostics.
!  This permits easier addition of new surface exchange
!  diagnostics without additional passing of arguments
!  though the boundary layer tree.
!  It also does not require the addition
!  of extra subroutine arguments when adding a new diagnostic.

!  Code Owner: Please refer to ModuleLeaders.txt
!  This file belongs in section: Surface

!- ----------------------------------------------------------------------

MODULE sf_diags_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE
SAVE

TYPE strnewsfdiag

  ! Need to create a flag and a pointer.
  ! In general initialise the flag to .FALSE. and use code elsewhere to reset
  ! this if required.
  ! Ordered by UM stash code

  LOGICAL :: l_z0h_eff_gb    = .FALSE. ! 27 effective rough len heat
  LOGICAL :: l_z0m_gb        = .FALSE. ! 28 vegetative rough length mom
  LOGICAL :: l_ra            = .FALSE. ! 54 aerodynamic resistance
  LOGICAL :: l_wt_ext        = .FALSE. ! 55 cumulative soil transp
  LOGICAL :: su10            = .FALSE. ! 209 u10m
  LOGICAL :: sv10            = .FALSE. ! 210 v10m
  LOGICAL :: sfme            = .FALSE. ! 224 wind mixing energy
  LOGICAL :: slh             = .FALSE. ! 234 latent heat
  LOGICAL :: sq_t1p5         = .FALSE. ! lots of things require 236 and 237
  LOGICAL :: st1p5           = .FALSE. ! 236, 328, 545, 547 1.5m temp
  LOGICAL :: sq1p5           = .FALSE. ! 237, 329, 546, 548 1.5m q
  LOGICAL :: simlt           = .FALSE. ! 257 sice mlt heat flux
  LOGICAL :: smlt            = .FALSE. ! 258 snowmelt heat flux
  LOGICAL :: l_gpp           = .FALSE. ! 261 gpp gb mean
  LOGICAL :: l_gpp_pft       = .FALSE. ! 289 gpp on pfts
  LOGICAL :: l_rib_surft     = .FALSE. ! 294 bulk ri on tiles
  LOGICAL :: l_fsmc_pft      = .FALSE. ! 313 soil moisture availability
  LOGICAL :: l_g_leaf_pft    = .FALSE. ! 325 leaf turnover
  LOGICAL :: l_rib_ssi       = .FALSE. ! 339 rib_ssi
  LOGICAL :: l_t10m          = .FALSE. ! 344 t at 10m over sea/sea-ice
  LOGICAL :: l_q10m          = .FALSE. ! 345 q at 10m over sea/sea-ice
  LOGICAL :: suv10m_n        = .FALSE. ! package of neutral winds
  LOGICAL :: l_u10m_n        = .FALSE. ! 368 u10m neutral
  LOGICAL :: l_v10m_n        = .FALSE. ! 369 v10m neutral
  LOGICAL :: l_mu10m_n       = .FALSE. ! 370 x pseudostress
  LOGICAL :: l_mv10m_n       = .FALSE. ! 371 y pseudostress
  LOGICAL :: l_radnet_sea    = .FALSE. ! 380 sea net radiation
  LOGICAL :: l_lw_surft      = .FALSE. ! 383, 384 up and down LW
  LOGICAL :: l_vshr_land     = .FALSE. ! 389 wind shear land
  LOGICAL :: l_vshr_ssi      = .FALSE. ! 390 wind shear sea sea ice
  LOGICAL :: l_ftemp         = .FALSE. ! 485 temp modifier resp
  LOGICAL :: l_fsth          = .FALSE. ! 486 soil moisture modifier resp
  LOGICAL :: l_fprf          = .FALSE. ! 489 veg modifier resp
  LOGICAL :: l_lw_up_sice_weighted_cat = .FALSE. ! 530 sea ice LW up
  LOGICAL :: l_lw_up_sice_weighted     = .FALSE. ! 531 sea ice LW up
  LOGICAL :: l_ftl_ice_sm    = .FALSE. ! 533 sea ice heat flux
  LOGICAL :: l_tstar_sice_weighted_cat = .FALSE. ! 534 sice temp
  LOGICAL :: l_tstar_sice_weighted     = .FALSE. ! 535 sice temp
  LOGICAL :: l_ice_present_cat         = .FALSE. ! 536 ice present
  LOGICAL :: l_ice_present             = .FALSE. ! 537 ice present
  LOGICAL :: l_cd_ssi        = .FALSE. ! 538 cd sea and sea ice
  LOGICAL :: l_et_stom       = .FALSE. ! 539 stom transp
  LOGICAL :: l_et_stom_surft = .FALSE. ! 540 stom transp
  LOGICAL :: l_ch_ssi        = .FALSE. ! 541 ch sea and sea ice
  LOGICAL :: l_lh_land       = .FALSE. ! 543 latent heat land
  LOGICAL :: l_lh_ssi        = .FALSE. ! 544 latent heat sea and sea ice
  LOGICAL :: l_snice         = .FALSE. ! 578-583 snow diags
  LOGICAL :: l_resp_p        = .FALSE. ! 663 plant respiration
  LOGICAL :: l_npp_pft       = .FALSE. ! 691 npp on pfts
  LOGICAL :: l_resp_p_pft    = .FALSE. ! 692 plant respiration
  LOGICAL :: l_tau_1         = .FALSE. ! 695 Grid-box mean momentum flux
  LOGICAL :: l_tau_surft     = .FALSE. ! 696 Momentum flux in surface tiles

  REAL(KIND=real_jlslsm), ALLOCATABLE :: z0h_eff_gb(:,:)
  !                    27      effective roughness length for heat
  REAL(KIND=real_jlslsm), ALLOCATABLE :: z0m_gb(:,:)
  !                    28      vegetative roughness length for momentum
  REAL(KIND=real_jlslsm), ALLOCATABLE :: ra(:)
  !                    54      aerodynamic resistance
  REAL(KIND=real_jlslsm), ALLOCATABLE :: wt_ext(:,:)
  !                    55      Culmulative transpiration from soil
  REAL(KIND=real_jlslsm), ALLOCATABLE :: u10m(:,:)
  !                    209      x-cpt of wind at 10m
  REAL(KIND=real_jlslsm), ALLOCATABLE :: v10m(:,:)
  !                    210      y-cpt of wind at 10m
  REAL(KIND=real_jlslsm), ALLOCATABLE :: fme(:,:)
  !                    224      Wind mixing energy
  REAL(KIND=real_jlslsm), ALLOCATABLE :: latent_heat(:,:)
  !                    234      Latent heat flux
  REAL(KIND=real_jlslsm), ALLOCATABLE :: t1p5m(:,:)
  !                    236      Temperature at 1.5m
  REAL(KIND=real_jlslsm), ALLOCATABLE :: q1p5m(:,:)
  !                    237      Specific humidity at 1.5m
  REAL(KIND=real_jlslsm), ALLOCATABLE :: sice_mlt_htf(:,:,:)
  !                    257      Sea ice surface melt heat flux
  REAL(KIND=real_jlslsm), ALLOCATABLE :: snomlt_surf_htf(:,:)
  !                    258      Surface snowmelt heat flux
  REAL(KIND=real_jlslsm), ALLOCATABLE :: gpp(:)
  !                    261      Gross prim prod grid box mean (kg C/m2/s)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: gpp_pft(:,:)
  !                    289      Gross prim prod on plant func types (kg C/m2/s)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: rib_surft(:,:)
  !                    294      bulk richardson number on land tiles
  REAL(KIND=real_jlslsm), ALLOCATABLE :: t1p5m_surft(:,:)
  !                    328      Temperature at 1.5m over tiles
  REAL(KIND=real_jlslsm), ALLOCATABLE :: q1p5m_surft(:,:)
  !                    329      Specific humidity at 1.5m over tiles
  REAL(KIND=real_jlslsm), ALLOCATABLE :: fsmc_pft(:,:)
  !                    313      Soil moisture availability on pfts
  REAL(KIND=real_jlslsm), ALLOCATABLE :: g_leaf_pft(:,:)
  !                    325      Leaf turnover rate on pfts
  REAL(KIND=real_jlslsm), ALLOCATABLE :: rib_ssi(:,:)
  !                    339      bulk richardson number sea and sea ice
  REAL(KIND=real_jlslsm), ALLOCATABLE :: chr10m(:,:)
  !                    CH at 10m over sea/sea-ice needed for 344/345
  REAL(KIND=real_jlslsm), ALLOCATABLE :: t10m(:,:)
  !                    344      Temperature at 10m over sea/sea-ice
  REAL(KIND=real_jlslsm), ALLOCATABLE :: q10m(:,:)
  !                    345      Specific humidity at 10m over sea/sea-ice
  REAL(KIND=real_jlslsm), ALLOCATABLE :: cdr10m_n(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: cdr10m_n_u(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: cdr10m_n_v(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: cd10m_n(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: cd10m_n_u(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: cd10m_n_v(:,:)
  ! coefficients for 10m neutral diagnostics
  REAL(KIND=real_jlslsm), ALLOCATABLE :: u10m_n(:, :)
  !                    368      x-cpt of neutral wind (at 10m)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: v10m_n(:, :)
  !                    369      y-cpt of neutral wind (at 10m)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: mu10m_n(:, :)
  !                    370      x-cpt of pseudostress
  REAL(KIND=real_jlslsm), ALLOCATABLE :: mv10m_n(:, :)
  !                    371      y-cpt of pseudostress
  REAL(KIND=real_jlslsm), ALLOCATABLE :: radnet_sea(:,:)
  !                    380      sea surface net radiation (W/m2)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: lw_up_surft(:,:)
  !                    383      Upwelling LW radiation on tiles
  REAL(KIND=real_jlslsm), ALLOCATABLE :: lw_down_surft(:,:)
  !                    384      Downwelling LW radiation on tiles
  REAL(KIND=real_jlslsm), ALLOCATABLE :: vshr_land(:,:) ! 389 wind shear land
  REAL(KIND=real_jlslsm), ALLOCATABLE :: vshr_ssi(:,:)
                                      ! 390 wind shear sea and sea ice
  REAL(KIND=real_jlslsm), ALLOCATABLE :: ftemp(:,:)
  !                    485      Temperature rate modifier of soil respiration
  REAL(KIND=real_jlslsm), ALLOCATABLE :: fsth(:,:)
  !                    486      Soil moisture rate modifier of soil respiration
  REAL(KIND=real_jlslsm), ALLOCATABLE :: fprf(:)
  !                    489      Soil respiration rate modifier due to vegetation cover
  REAL(KIND=real_jlslsm), ALLOCATABLE :: lw_up_sice_weighted_cat(:,:,:)
  !                    530      category ice area weighted upward LW flux
  !                             over sea ice after boundary layer
  !                             calculation
  REAL(KIND=real_jlslsm), ALLOCATABLE :: lw_up_sice_weighted(:,:)
  !                    531      category ice area weighted upward LW flux
  !                             over sea ice after boundary layer
  !                             calculation
  REAL(KIND=real_jlslsm), ALLOCATABLE :: ftl_ice_sm(:,:)
  !                    533      aggregate ice area weighted sensible
  !                             heat flux over sea ice
  REAL(KIND=real_jlslsm), ALLOCATABLE :: tstar_sice_weighted_cat(:,:,:)
  !                    534      category ice area weighted sea ice surface
  !                             skin temperature
  REAL(KIND=real_jlslsm), ALLOCATABLE :: tstar_sice_weighted(:,:)
  !                    535      category ice area weighted sea ice surface
  !                             skin temperature
  REAL(KIND=real_jlslsm), ALLOCATABLE :: ice_present_cat(:,:,:)
  !                    536      category sea ice time fraction
  REAL(KIND=real_jlslsm), ALLOCATABLE :: ice_present(:,:)
  !                    537      category sea ice time fraction
  REAL(KIND=real_jlslsm), ALLOCATABLE :: cd_ssi(:,:)
  !                    538       Sea and sea ice drag coefficient (momentum)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: et_stom_ij(:,:)
  !                    539      Transpiration through stom (kg/m2/s)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: et_stom_surft(:,:)
  !                    540      Transpiration through stom (kg/m2/s)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: resfs_stom(:,:)
  !                    Combined stomatal and aerodynamic
  !                    resistance factor for fraction 1-FRACA.
  REAL(KIND=real_jlslsm), ALLOCATABLE :: ch_ssi(:,:)
  !                    541       Sea and sea ice drag coefficient (heat)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: lh_land(:)
  !                    543       Land latent heat flux
  REAL(KIND=real_jlslsm), ALLOCATABLE :: lh_ssi(:,:)
  !                    544       Sea and sea ice mean latent heat
  REAL(KIND=real_jlslsm), ALLOCATABLE :: t1p5m_ssi(:,:)
  !                    545      Temperature at 1.5m over sea/sea-ice
  REAL(KIND=real_jlslsm), ALLOCATABLE :: q1p5m_ssi(:,:)
  !                    546      Specific humidity at 1.5m over sea/sea-ice
  REAL(KIND=real_jlslsm), ALLOCATABLE :: snice_smb_surft(:,:)
  !                    578      Tiled snow mass rate of change
  !                             (kg/m2/s)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: snice_m_surft(:,:)
  !                    579      Total internal melt rate of snowpack
  !                             (kg/m2/s)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: snice_freez_surft(:,:)
  !                    580      Total internal refreezing rate in
  !                             snowpack (kg/m2/s)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: snice_runoff_surft(:,:)
  !                    581      Net rate of liquid leaving snowpack
  !                             (kg/m2/s)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: snice_sicerate_surft(:,:)
  !                    582      Rate of change of solid mass in
  !                             snowpack (kg/m2/s)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: snice_sliqrate_surft(:,:)
  !                    583      Rate of change of liquid mass
  !                             snowpack (kg/m2/s)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: resp_p(:)
  !                    663      Plant respiration (kg/m2/s)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: npp_pft(:,:)
  !                    691      Net primary productivity on PFTs (kg C/m2/s)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: resp_p_pft(:,:)
  !                    692      Plant respiration on PFTs (kg/m2/s)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: tau_1(:,:)
  !                    695      Gridbox mean surface momentum flux on P-grid
  !                             (N/m2)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: tau_surft(:,:)
  !                    695      Tiled surface momentum flux on P-grid
  !                             (N/m2)

END TYPE strnewsfdiag

TYPE (Strnewsfdiag) :: sf_diag

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SF_DIAGS_MOD'
! ----------------------------------------------------------------------

END MODULE sf_diags_mod
