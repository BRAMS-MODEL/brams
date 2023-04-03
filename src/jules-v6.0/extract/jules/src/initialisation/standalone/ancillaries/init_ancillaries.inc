#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_ancillaries(nml_dir, crop_vars, ainfo, trif_vars, urban_param,&
                            trifctltype)

USE io_constants, ONLY: namelist_unit

USE string_utils_mod, ONLY: to_string

USE aero, ONLY: co2_mmr_mod => co2_mmr

USE dump_mod, ONLY: ancil_dump_read

USE init_flake_ancils_mod, ONLY: init_flake_ancils

USE errormessagelength_mod, ONLY: errormessagelength

USE jules_surface_mod, ONLY: l_urban2t

USE switches_urban, ONLY: l_moruses

USE um_types, ONLY: real_jlslsm

USE jules_model_environment_mod, ONLY: lsm_id, cable

USE allocate_cable_arrays_mod, ONLY: allocate_cable_progs

!TYPE definitions
USE crop_vars_mod, ONLY: crop_vars_type
USE ancil_info,    ONLY: ainfo_type
USE trif_vars_mod, ONLY: trif_vars_type
USE urban_param_mod, ONLY: urban_param_type
USE trifctl,   ONLY: trifctl_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in information about the model ancillaries
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Arguments
CHARACTER(LEN=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                         ! namelists

!TYPES containing field data (IN OUT)
TYPE(crop_vars_type), INTENT(IN OUT) :: crop_vars
TYPE(ainfo_type), INTENT(IN OUT) :: ainfo
TYPE(trif_vars_type), INTENT(IN OUT) :: trif_vars
TYPE(urban_param_type), INTENT(IN OUT) :: urban_param
TYPE(trifctl_type), INTENT(IN OUT) :: trifctltype

CHARACTER(LEN=errormessagelength) :: iomessage

! Work variables
INTEGER :: error  ! Error indicator

LOGICAL :: read_from_dump
REAL(KIND=real_jlslsm)    :: co2_mmr

CHARACTER(LEN=*), PARAMETER :: RoutineName='init_ancillaries'

NAMELIST  / jules_co2 / read_from_dump, co2_mmr

!-----------------------------------------------------------------------------
! Open the ancillaries namelist file
OPEN(namelist_unit, FILE=(TRIM(nml_dir) // '/' // 'ancillaries.nml'),         &
               STATUS='old', POSITION='rewind', ACTION='read', IOSTAT = error,&
               IOMSG = iomessage)
IF ( error /= 0 ) THEN
  CALL log_fatal(RoutineName,                                                 &
                 "Error opening namelist file ancillaries.nml " //            &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")
END IF

! Defer to specialist routines to process each namelist
CALL init_frac()
CALL init_vegetation_props()
CALL init_soil_props()
CALL init_top()
CALL init_pdm()
CALL init_agric(trif_vars%frac_past_gb,trifctltype)
CALL init_crop_props()
CALL init_irrig_props()
CALL init_rivers_props()
CALL init_water_resources_props()
CALL init_flake_ancils()
IF ( l_urban2t .OR. l_moruses ) CALL init_urban_props(ainfo,urban_param)
IF (lsm_id == cable) THEN
  CALL allocate_cable_progs()
  CALL init_cable_progs()
END IF

! Read the JULES_CO2 namelist
! This is so simple it doesn't really need its own subroutine
read_from_dump = .FALSE.

!Copy across the default value set up in the module.
co2_mmr = co2_mmr_mod

CALL log_info(RoutineName, "Reading JULES_CO2 namelist...")
READ(namelist_unit, NML = jules_co2, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 ) THEN
  CALL log_fatal(RoutineName,                                                 &
                 "Error reading namelist JULES_CO2 " //                       &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")
END IF

ancil_dump_read%co2 = read_from_dump

IF ( .NOT. ancil_dump_read%co2) THEN
  co2_mmr_mod = co2_mmr
ELSE
  ! We read from the dump file.
  CALL log_info(RoutineName,                                                  &
                "co2_mmr will be read from the dump file.  " //               &
                "Namelist value ignored")
END IF

! Close the namelist file
CLOSE(namelist_unit, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 ) THEN
  CALL log_fatal(RoutineName,                                                 &
                 "Error closing namelist file ancillaries.nml " //            &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")
END IF

RETURN

END SUBROUTINE init_ancillaries
#endif
