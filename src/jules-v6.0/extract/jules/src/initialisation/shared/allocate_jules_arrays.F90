MODULE allocate_jules_arrays_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ALLOCATE_JULES_ARRAYS_MOD'

CONTAINS
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine ALLOCATE_JULES_ARRAYS
!
! Description: Routine that allocates memory to the JULES arrays
! This assume that the values in the jules_surface_types module have been set
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.
!
!   Code Owner: Please refer to ModuleLeaders.txt and UM file CodeOwners.txt
!   This file belongs in section: Technical

SUBROUTINE allocate_jules_arrays(crop_vars_data,psparms_data,top_pdm_data,    &
                                 fire_vars_data,ainfo_data, trif_vars_data,   &
                                 soil_ecosse_vars_data, aero_data,            &
                                 urban_param_data,progs_data,trifctl_data,    &
                                 coastal_data,jules_vars_data)

!Variables- switches
USE jules_vegetation_mod,     ONLY: l_crop
USE jules_vegetation_mod,     ONLY: l_triffid, l_phenol, l_nitrogen,          &
                                    l_use_pft_psi
USE jules_irrig_mod,          ONLY: l_irrig_dmd, irr_crop, irr_crop_doell
USE jules_surface_mod,        ONLY: l_urban2t
USE switches_urban,           ONLY: l_moruses
USE jules_surface_mod,        ONLY: l_flake_model
USE jules_soil_mod,           ONLY: l_bedrock, ns_deep
USE jules_radiation_mod,      ONLY: l_albedo_obs
USE jules_soil_biogeochem_mod,ONLY: soil_model_ecosse, soil_bgc_model,        &
                                    l_layeredc, dim_ch4layer
USE jules_deposition_mod,     ONLY: dry_dep_model, l_deposition,              &
                                    ndry_dep_species

!Variables- dimensions
USE jules_surface_types_mod,  ONLY: ncpft,nnpft
USE jules_snow_mod,           ONLY: nsmax, cansnowtile
USE jules_surface_types_mod,  ONLY: npft, nnvg, ntype
USE theta_field_sizes,        ONLY: t_i_length, t_j_length,                   &
                                    u_i_length,u_j_length,                    &
                                    v_i_length,v_j_length
USE atm_fields_bounds_mod,    ONLY: pdims_s, pdims                                    
USE ancil_info,               ONLY: nsoilt, dim_cs1, dim_cslayer,             &
                                    dim_soil_n_pool,  land_pts,               &
                                    nsurft, nmasst, rad_nband
USE dust_parameters_mod,      ONLY: ndiv
USE nlsizes_namelist_mod,     ONLY: bl_levels
USE jules_sea_seaice_mod,     ONLY: nice, nice_use
USE jules_soil_mod,           ONLY: sm_levels

!Subroutines
USE bvoc_vars,                ONLY: bvocvars_alloc
USE crop_vars_mod,            ONLY: crop_vars_alloc
USE cropparm,                 ONLY: cropparm_alloc
USE c_z0h_z0m,                ONLY: c_z0h_z0m_alloc
USE fire_vars_mod,            ONLY: fire_vars_alloc
USE fluxes,                   ONLY: fluxes_alloc
USE jules_vars_mod,           ONLY: jules_vars_alloc
USE jules_irrig_mod,          ONLY: irrig_vars_alloc
USE metstats_mod,             ONLY: metstats_allocate
USE nvegparm,                 ONLY: nvegparm_alloc
USE ozone_vars,               ONLY: ozone_vars_alloc
USE pftparm,                  ONLY: pftparm_alloc
USE prognostics,              ONLY: prognostics_alloc
USE p_s_parms,                ONLY: psparms_alloc
USE top_pdm,                  ONLY: top_pdm_alloc
USE trif,                     ONLY: trif_alloc
USE trifctl,                  ONLY: trifctl_alloc
USE trif_vars_mod,            ONLY: trif_vars_alloc
USE urban_param_mod,          ONLY: urban_param_alloc
USE lake_mod,                 ONLY: lake_alloc
USE ancil_info,               ONLY: ancil_info_alloc
USE veg3_parm_mod,            ONLY: veg3_parm_allocate
USE veg3_field_mod,           ONLY: veg3_field_allocate
USE jules_rivers_mod,         ONLY: jules_rivers_alloc

#if !defined(UM_JULES)
USE gridmean_fluxes,          ONLY: gridmean_fluxes_alloc
USE forcing,                  ONLY: forcing_alloc
USE jules_deposition_mod,     ONLY: jules_deposition_alloc
USE deposition_species_mod,   ONLY: deposition_species_alloc
USE jules_water_resources_mod, ONLY: water_resources_alloc
#endif
!Moved from ifdef
USE soil_ecosse_vars_mod,     ONLY: soil_ecosse_vars_alloc
USE aero,                     ONLY: aero_alloc
USE coastal,                  ONLY: coastal_alloc

!TYPE definitions
USE crop_vars_mod, ONLY: crop_vars_data_type
USE p_s_parms,     ONLY: psparms_data_type
USE top_pdm,       ONLY: top_pdm_data_type
USE fire_vars_mod, ONLY: fire_vars_data_type
USE ancil_info,    ONLY: ainfo_data_type
USE trif_vars_mod, ONLY: trif_vars_data_type
USE soil_ecosse_vars_mod, ONLY: soil_ecosse_vars_data_type
USE aero,          ONLY: aero_data_type
USE urban_param_mod, ONLY:urban_param_data_type
USE prognostics,   ONLY: progs_data_type
USE trifctl,       ONLY: trifctl_data_type
USE coastal,       ONLY: coastal_data_type
USE jules_vars_mod, ONLY: jules_vars_data_type

USE parkind1,                 ONLY: jprb, jpim
USE yomhook,                  ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments

!TYPES containing field data (IN OUT)
TYPE(crop_vars_data_type), INTENT(IN OUT) :: crop_vars_data
TYPE(psparms_data_type), INTENT(IN OUT) :: psparms_data
TYPE(top_pdm_data_type), INTENT(IN OUT) :: top_pdm_data
TYPE(fire_vars_data_type), INTENT(IN OUT) :: fire_vars_data
TYPE(ainfo_data_type), INTENT(IN OUT) :: ainfo_data
TYPE(trif_vars_data_type), INTENT(IN OUT) :: trif_vars_data
TYPE(soil_ecosse_vars_data_type), INTENT(IN OUT) :: soil_ecosse_vars_data
TYPE(aero_data_type), INTENT(IN OUT) :: aero_data
TYPE(urban_param_data_type), INTENT(IN OUT) :: urban_param_data
TYPE(progs_data_type), INTENT(IN OUT) :: progs_data
TYPE(trifctl_data_type), INTENT(IN OUT) :: trifctl_data
TYPE(coastal_data_type), INTENT(IN OUT) :: coastal_data
TYPE(jules_vars_data_type), INTENT(IN OUT) :: jules_vars_data

!Local variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_JULES_ARRAYS'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!Any dimension sizes should be set before we get here. Some special cases for
!UM mode can be found in surf_couple_allocate.

CALL prognostics_alloc(land_pts, t_i_length, t_j_length,                      &
                      nsurft, npft, nsoilt, sm_levels, ns_deep, nsmax,        &
                      dim_cslayer, dim_cs1, dim_ch4layer,                     &
                      nice, nice_use, soil_bgc_model, soil_model_ecosse,      &
                      l_layeredc, l_triffid, l_phenol, l_bedrock,             &
                      progs_data)

CALL fluxes_alloc(land_pts, t_i_length, t_j_length,                           &
                  nsurft, npft, nsoilt, sm_levels,                            &
                  nice, nice_use)

CALL bvocvars_alloc(land_pts,npft)

CALL crop_vars_alloc(land_pts, t_i_length, t_j_length,                        &
                     nsurft, ncpft,nsoilt, sm_levels, l_crop, irr_crop,       &
                     irr_crop_doell, crop_vars_data)

CALL irrig_vars_alloc(npft, l_irrig_dmd)

CALL cropparm_alloc(ncpft,l_crop)

CALL fire_vars_alloc(land_pts,npft, fire_vars_data)

CALL c_z0h_z0m_alloc(ntype)

CALL jules_vars_alloc(land_pts,ntype,nsurft,rad_nband,nsoilt,sm_levels,       &
                t_i_length, t_j_length, npft, bl_levels, pdims_s, pdims,      &
                l_albedo_obs, cansnowtile, l_deposition,                      &
                jules_vars_data)

CALL metstats_allocate(land_pts)

CALL nvegparm_alloc(nnvg)

CALL ozone_vars_alloc(land_pts,npft)

CALL pftparm_alloc(npft)

CALL psparms_alloc(land_pts,t_i_length,t_j_length,                            &
                   nsoilt,sm_levels,dim_cslayer,nsurft,npft,                  &
                   soil_bgc_model,soil_model_ecosse,l_use_pft_psi,            &
                   psparms_data)

CALL top_pdm_alloc(land_pts,nsoilt, top_pdm_data)

CALL trif_alloc(npft,                                                         &
                l_triffid, l_phenol)

CALL trifctl_alloc(land_pts,                                                  &
                   npft,dim_cslayer,dim_cs1,nsoilt,trifctl_data)

CALL trif_vars_alloc(land_pts,                                                &
                     npft,dim_cslayer,nsoilt,dim_cs1,                         &
                     l_triffid, l_phenol, trif_vars_data)

CALL veg3_parm_allocate(land_pts,nsurft,nnpft,npft)
CALL veg3_field_allocate(land_pts,nsurft,nnpft,npft,nmasst)

CALL urban_param_alloc(land_pts,                                              &
                       l_urban2t, l_moruses, urban_param_data)

CALL lake_alloc(land_pts,                                                     &
                l_flake_model)

CALL ancil_info_alloc(land_pts,t_i_length,t_j_length,                         &
                      nice,nsoilt,ntype,                                      &
                      ainfo_data)

CALL jules_rivers_alloc(land_pts)

#if !defined(UM_JULES)
CALL gridmean_fluxes_alloc(t_i_length,t_j_length)

CALL forcing_alloc(t_i_length,t_j_length)

CALL jules_deposition_alloc(land_pts)

CALL deposition_species_alloc(ntype,ndry_dep_species,                         &
                              l_deposition, dry_dep_model)

CALL water_resources_alloc(land_pts)
#endif

!Moved from ifdef
CALL soil_ecosse_vars_alloc(land_pts,                                         &
                            nsoilt,dim_cslayer,dim_soil_n_pool,sm_levels,     &
                            soil_bgc_model,soil_model_ecosse,                 &
                            soil_ecosse_vars_data)

CALL aero_alloc(land_pts,t_i_length,t_j_length,                               &
                nsurft,ndiv, aero_data)

CALL coastal_alloc(land_pts,t_i_length,t_j_length,                            &
                   u_i_length,u_j_length,                                     &
                   v_i_length,v_j_length,                                     &
                   nice_use,nice,coastal_data)


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE allocate_jules_arrays

END MODULE allocate_jules_arrays_mod
