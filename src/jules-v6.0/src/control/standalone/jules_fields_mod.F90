! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

MODULE jules_fields_mod

! Module containing instances of the data TYPES for JULES standalone

!TYPE definitions
USE crop_vars_mod, ONLY: crop_vars_type, crop_vars_data_type
USE p_s_parms,     ONLY: psparms_type, psparms_data_type
USE top_pdm,       ONLY: top_pdm_type, top_pdm_data_type
USE fire_vars_mod, ONLY: fire_vars_type, fire_vars_data_type
USE ancil_info,    ONLY: ainfo_type, ainfo_data_type
USE trif_vars_mod, ONLY: trif_vars_type, trif_vars_data_type
USE soil_ecosse_vars_mod, ONLY: soil_ecosse_vars_type, soil_ecosse_vars_data_type
USE aero,          ONLY: aero_type, aero_data_type
USE urban_param_mod, ONLY: urban_param_type, urban_param_data_type
USE prognostics, ONLY: progs_type, progs_data_type
USE trifctl,       ONLY: trifctl_type, trifctl_data_type
USE coastal,       ONLY: coastal_data_type, coastal_type
USE jules_vars_mod, ONLY: jules_vars_data_type, jules_vars_type

!Declare instances of the TYPES required to hold the data
!TYPES to hold the data
TYPE(crop_vars_data_type), TARGET :: crop_vars_data
TYPE(psparms_data_type), TARGET :: psparms_data
TYPE(top_pdm_data_type), TARGET :: top_pdm_data
TYPE(fire_vars_data_type), TARGET :: fire_vars_data
TYPE(ainfo_data_type), TARGET :: ainfo_data
TYPE(trif_vars_data_type), TARGET :: trif_vars_data
TYPE(soil_ecosse_vars_data_type), TARGET :: soil_ecosse_vars_data
TYPE(aero_data_type), TARGET :: aero_data
TYPE(urban_param_data_type), TARGET :: urban_param_data
TYPE(progs_data_type), TARGET :: progs_data
TYPE(trifctl_data_type), TARGET :: trifctl_data
TYPE(coastal_data_type), TARGET :: coastal_data
TYPE(jules_vars_data_type), TARGET :: jules_vars_data

!TYPES we pass around. These happen to be pointers to the data types above
!but this should be transparent
TYPE(crop_vars_type) :: crop_vars
TYPE(psparms_type) :: psparms
TYPE(top_pdm_type) :: toppdm
TYPE(fire_vars_type) :: fire_vars
TYPE(ainfo_type) ::  ainfo
TYPE(trif_vars_type) :: trif_vars
TYPE(soil_ecosse_vars_type) :: soilecosse
TYPE(aero_type) :: aerotype
TYPE(urban_param_type) :: urban_param
TYPE(progs_type) :: progs
TYPE(trifctl_type) :: trifctltype
TYPE(coastal_type) :: coast
TYPE(jules_vars_type) :: jules_vars

END MODULE jules_fields_mod
