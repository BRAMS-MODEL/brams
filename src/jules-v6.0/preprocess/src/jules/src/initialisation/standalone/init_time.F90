! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_time(nml_dir)

USE io_constants, ONLY: namelist_unit

USE string_utils_mod, ONLY: to_string

USE datetime_mod, ONLY: datetime_str_len, secs_in_day,                        &
                         datetime,                                            &
! Also import the comparison operators on the datetime type
                           operator( == ), operator( /= ), operator( < ),     &
                           operator( > ), operator( <= ), operator( >= ),     &
                           datetime_from_string, datetime_to_string,          &
                           datetime_advance, datetime_subtract,               &
                           l_360, l_leap

USE model_interface_mod, ONLY: identifier_len
USE timestep_mod, ONLY: timestep_real => timestep
USE model_time_mod, ONLY:                                                     &
  run_min_time, run_max_time, current_time, timestep, timestep_len,           &
  print_step,                                                                 &
  is_spinup, max_spinup_cycles, spinup_cycle, terminate_on_spinup_fail,       &
! Alias some variables to different names
    main_run_start_dt => main_run_start, main_run_end_dt => main_run_end,     &
    spinup_start_dt => spinup_start, spinup_end_dt => spinup_end,             &
    start_of_year, end_of_year, timesteps_in_day

USE spinup_mod, ONLY: max_spinup_vars, nvars, spinup_vars

USE trifctl, ONLY: asteps_since_triffid

USE switches, ONLY: l_360_sw => l_360

USE logging_mod, ONLY: log_info, log_warn, log_fatal

USE errormessagelength_mod, ONLY: errormessagelength

USE um_types, ONLY: real_jlslsm

USE mem_brams_jules, ONLY: timestepB,main_run_endB,main_run_startB  !DSM

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in information about the run times and spinup and checks for
!   consistency
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

! Work variables
TYPE(datetime) :: dt  ! Placeholder for a datetime in a calculation

INTEGER :: error  ! Error indicator

INTEGER :: i  ! Loop counter

!-----------------------------------------------------------------------------
! Definition of jules_time namelist - this brings together some variables
! defined in disparate modules, and also reads some local variables that
! are transformed in this routine
!-----------------------------------------------------------------------------
CHARACTER(LEN=datetime_str_len) :: main_run_start, main_run_end
                                   ! Character strings for times of main run
NAMELIST  / jules_time/ l_360, l_leap, timestep_len, main_run_start,          &
                      main_run_end, print_step

!-----------------------------------------------------------------------------
! Definition of jules_spinup namelist - this reads some variables in spinup_mod
! and some local variables that are transformed in this routine
!-----------------------------------------------------------------------------
CHARACTER(LEN=datetime_str_len) :: spinup_start, spinup_end
                              ! Character strings for times of spinup cycles
CHARACTER(LEN=identifier_len) :: var(max_spinup_vars)
                              ! Identifiers of the variables to use for
                              ! spinup
LOGICAL :: use_percent(max_spinup_vars)
                              ! T - use a percentage of the previous value
                              !     as a tolerance
                              ! F - use an absolute tolerance
REAL(KIND=real_jlslsm) :: tolerance(max_spinup_vars)
                              ! The tolerance to use
CHARACTER(LEN=errormessagelength) :: iomessage
NAMELIST  / jules_spinup/ max_spinup_cycles, spinup_start, spinup_end,        &
                        terminate_on_spinup_fail, nvars, var, use_percent,    &
                        tolerance


!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
use_percent(:) = .FALSE.  ! Default is to use absolute tolerance for all
                          ! spinup variables

!-----------------------------------------------------------------------------
! Read namelists
!-----------------------------------------------------------------------------
OPEN(namelist_unit, FILE=(TRIM(nml_dir) // '/' // 'timesteps.nml'),           &
               STATUS='old', POSITION='rewind', ACTION='read', IOSTAT = error,&
               IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_time",                                                 &
                 "Error opening namelist file timesteps.nml " //              &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

! First read the time namelist
CALL log_info("init_time", "Reading JULES_TIME namelist...")
READ(namelist_unit, NML = jules_time, IOSTAT = error, IOMSG = iomessage)

timestep_len=timestepB          !DSM 
main_run_end=main_run_endB      !DSM 
main_run_start=main_run_startB  !DSM 

IF ( error /= 0 )                                                             &
  CALL log_fatal("init_time",                                                 &
                 "Error reading namelist JULES_TIME " //                      &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

! Then read the spinup namelist
CALL log_info("init_time", "Reading JULES_SPINUP namelist...")
READ(namelist_unit, NML = jules_spinup, IOSTAT = error, IOMSG = iomessage)

!{DSM
IF (max_spinup_cycles/=0) THEN
   print*
   print*, 'max_spinup_cycles deve ser igual a 0 (zero) nesta versao acoplada!!!'
   STOP
ENDIF
!DSM}

IF ( error /= 0 )                                                             &
  CALL log_fatal("init_time",                                                 &
                 "Error reading namelist JULES_SPINUP " //                    &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

CLOSE(namelist_unit, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_time",                                                 &
                 "Error closing namelist file timesteps.nml " //              &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")


!-----------------------------------------------------------------------------
! Set values derived from namelists and verify for consistency
!-----------------------------------------------------------------------------
! Set l_360 in switches to the same value as in datetime_mod
l_360_sw = l_360
IF ( l_360 )                                                                  &
  CALL log_info("init_time", "Running with 360 day year")

IF ( .NOT. l_leap )                                                           &
  CALL log_info("init_time", "No leap years")

! For simplicity, we enforce that timestep is a divisor of one day
IF ( MOD(secs_in_day, timestep_len) /= 0 )                                    &
  CALL log_fatal("init_time",                                                 &
                 "A day must contain a whole number of timesteps")
CALL log_info("init_time",                                                    &
              "Timestep is " // TRIM(to_string(timestep_len)) // " seconds")

! Convert character strings for main run times into datetimes
main_run_start_dt = datetime_from_string(main_run_start)
main_run_end_dt = datetime_from_string(main_run_end)

! Check that the main run times are in chronological order
IF ( main_run_end_dt <= main_run_start_dt )                                   &
  CALL log_fatal("init_time",                                                 &
                 "Start time for main run must be before end time")

CALL log_info("init_time",                                                    &
              "Main run start - " // datetime_to_string(main_run_start_dt))
CALL log_info("init_time",                                                    &
              "Main run end - " // datetime_to_string(main_run_end_dt))

! Assume that there is no spinup for now
! I.e. the main run is the entire run, and the time for the first timestep
! is the start of the run
run_min_time = main_run_start_dt
run_max_time = main_run_end_dt
current_time = main_run_start_dt

! Check for spinup
is_spinup = max_spinup_cycles > 0

IF ( is_spinup .AND. nvars < 1 ) THEN
  ! If no spinup variables were specified, that is the same as specifying no
  ! spinup
  CALL log_warn("init_time",                                                  &
                "Spinup requested (max_spinup_cycles > 0) but no " //         &
                "variables specified - spinup will not be performed")

  is_spinup         = .FALSE.
  max_spinup_cycles = 0
END IF

IF ( is_spinup ) THEN
  CALL log_info("init_time",                                                  &
                "Spinup requested - maximum number of spinup cycles is " //   &
                TRIM(to_string(max_spinup_cycles)))

  ! Convert character strings for spinup times into datetimes
  spinup_start_dt = datetime_from_string(spinup_start)
  spinup_end_dt = datetime_from_string(spinup_end)

  ! Check that the spinup times are in chronological order
  IF ( spinup_end_dt <= spinup_start_dt )                                     &
    CALL log_fatal("init_time",                                               &
                   "Start time for spinup must be before end time")

  CALL log_info("init_time",                                                  &
                "Spinup start - " // datetime_to_string(spinup_start_dt))
  CALL log_info("init_time",                                                  &
                "Spinup end - " // datetime_to_string(spinup_end_dt))

  ! Adjust the run times as required
  IF ( spinup_start_dt < main_run_start_dt ) run_min_time = spinup_start_dt
  IF ( spinup_end_dt > main_run_end_dt ) run_max_time = spinup_end_dt
  ! The time for the first timestep is now the start of spinup
  current_time = spinup_start_dt

  ! We start in the first cycle of spinup
  spinup_cycle = 1

  ! Set up the spinup variables
  IF ( nvars > max_spinup_vars )                                              &
    CALL log_fatal("init_time", "Too many spinup variables specified")

  DO i = 1,nvars
    spinup_vars(i)%identifier  = var(i)
    spinup_vars(i)%use_percent = use_percent(i)
    spinup_vars(i)%tolerance   = tolerance(i)
  END DO
ELSE
  CALL log_info("init_time", "No spinup requested")
END IF

! We are (obviously) starting in the first timestep
timestep = 1
! Initialise the number of timesteps since triffid was called
asteps_since_triffid = 1

! Set the variables that indicate if we are starting or ending a year
! If going back one timestep length causes the year to change, then we are in
! the first timestep of a year
dt = datetime_subtract(current_time, timestep_len)
start_of_year = dt%year /= current_time%year
! If advancing the current time causes the year to change, then it is the last
! timestep in a year
dt = datetime_advance(current_time, timestep_len)
end_of_year = dt%year /= current_time%year

timesteps_in_day = secs_in_day / timestep_len

! Initialise real_timestep for use in science routines
timestep_real = REAL(timestep_len)

RETURN

END SUBROUTINE init_time
