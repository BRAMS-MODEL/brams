!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237] 
!******************************COPYRIGHT**************************************
SUBROUTINE imogen_check(                                                      &
  c_emissions,include_co2,include_non_co2,land_feed_co2, land_feed_ch4,       &
  ocean_feed,anlg,anom, wgen                                                  &
)

USE logging_mod, ONLY: log_fatal, log_warn

USE jules_print_mgr, ONLY:                                                    &
  jules_message,                                                              &
  jules_print
IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   This code is designed to make a quick cross-check to ensure that
!   the flags in IMOGEN are set properly for the currently allowed configu
!   The currently allowed configurations should fit with the available doc
!
! Code Owner: Please refer to ModuleLeaders.txt
!             This file belongs in IMOGEN
! Written by: C. Huntingford (4th March 2004) then updated by Edward Comyn-Platt
! and Eleanor Burke to include the checks in imogen_confirmed_run.F90 and remove
! any duplication.
! 
! Code Description:
!   Language: Fortran 90.
!   
!-----------------------------------------------------------------------------

LOGICAL, INTENT(IN) ::                                                        &
  c_emissions,                                                                &
         !IN True: CO2 concentration is calculated from emissions read from
         !         file_scen_emits
         !   False: CO2 concentration read from file_scen_co2_ppmv
  include_co2,                                                                &
         !IN Are adjustments to CO2 values allowed?
  include_non_co2,                                                            &
         !IN Are adjustments to non-CO2 values allowed?
  land_feed_co2,                                                              &
         !IN Are land CO2 feedbacks allowed on atmospheric C
  land_feed_ch4,                                                              &
         !IN Are land CH4 feedbacks allowed on atmospheric C
  ocean_feed,                                                                 &
         !IN Are ocean feedbacks allowed on atmospheric C
  anlg,                                                                       &
         !IN True if the analogue model is used
  anom,                                                                       &
         !IN True if the anomalies are used
  wgen   !IN True if weather generator is used

!-----------------------------------------------------------------------------
! Local variables.
!-----------------------------------------------------------------------------
LOGICAL ::                                                                    &
  check_flag,                                                                 &
         !WORK Checks that configuration is Valid
  confirm_flag
         !WORK Confirms that configuration has been used before


check_flag   = .FALSE.
confirm_flag = .FALSE.

! initialise imogen_check message with line of stars
WRITE(jules_message,*) '***********************************************'
CALL jules_print('imogen_check',jules_message)

! Now make checks on structures which will hopefully all be available at
! To make this code easy to read, simple "if"
! statements comments are used for each instance and case printed.

! User prescribed anomalies
IF (anom .AND. ( .NOT. anlg)) THEN
  WRITE(jules_message,*) 'Run is using user prescribed anomalies'
  CALL jules_print('imogen_check',jules_message)
  check_flag = .TRUE.
END IF

IF ( .NOT. anom) THEN    
  check_flag   = .TRUE.
  WRITE(jules_message,*) 'Run is without meteorology anomalies'
  CALL jules_print('imogen_check',jules_message)
  WRITE(jules_message,*)                                                      &
    'If spin-up, please ensure that flags are compatible with transient run'
  CALL jules_print('imogen_check',jules_message)
    
  IF ( .NOT. anlg) THEN
    WRITE(jules_message,*) 'Hydrology of 20th Century simulation'
    CALL jules_print('imogen_check',jules_message)
    
    IF (c_emissions) THEN
      WRITE(jules_message,*) 'Allows user to provide Anthro. CO2 emissions'
      CALL jules_print('imogen_check',jules_message)
    ELSE ! not c_emissions
      IF (include_co2) THEN
        WRITE(jules_message,*) 'Allows user to provide a file of CO2 concs.'
        CALL jules_print('imogen_check',jules_message)
      ELSE ! not include_co2 and not c_emissions
        WRITE(jules_message,*) 'Simulation with fixed CO2 concentration'
        CALL jules_print('imogen_check',jules_message)
      END IF
    END IF

    ! Confirmed configurations (i.e. people have used this configuration before)
    IF (( .NOT. c_emissions) .AND. ( .NOT. land_feed_ch4) .AND.               &
        ( .NOT. include_non_co2) .AND. ( .NOT. land_feed_co2) .AND.           &
        ( .NOT. ocean_feed) .AND. ( .NOT. wgen))                              &
        THEN !include_co2 is either true or false
      confirm_flag = .TRUE. 
    END IF

  END IF  ! end check for not anlg model
END IF

! Currently coded analogue model possibilities
IF (anom .AND. anlg) THEN
  WRITE(jules_message,*) 'Analogue model simulations'
  CALL jules_print('imogen_check',jules_message)
  
  IF (( .NOT. include_co2) .AND. include_non_co2) THEN
    WRITE(jules_message,*) 'Run is for non-CO2 gases only'
    CALL jules_print('imogen_check',jules_message)
    check_flag = .TRUE.
  END IF
  
  IF (include_co2) THEN
    ! All configurations with co2 are technically valid
    check_flag = .TRUE.

    IF ( .NOT. c_emissions) THEN
      ! Write out CO2 concentrations statement  
      WRITE(jules_message,*) 'Run is for prescribed CO2'
      CALL jules_print('imogen_check',jules_message)
      WRITE(jules_message,*) 'Allows user to provides a file of CO2 concs.'
      CALL jules_print('imogen_check',jules_message)
      
      ! Check for confirmed configurations:
      IF (( .NOT. land_feed_co2) .AND. ( .NOT. ocean_feed) .AND.              &
          ( .NOT. land_feed_ch4) .AND. ( .NOT. wgen)) THEN
        confirm_flag = .TRUE.
      END IF
      
    ELSE ! if c_emissions true
      
      ! Write out CO2 emissions statement  
      WRITE(jules_message,*) 'Carbon cycle driven by prescribed emissions'
      CALL jules_print('imogen_check',jules_message)
      WRITE(jules_message,*) 'Allows user to provides a file of CO2 emissions.'
      CALL jules_print('imogen_check',jules_message)
      
      ! Write which CO2 feedbacks are turned on 
      IF (land_feed_co2) THEN
        WRITE(jules_message,*) 'There are land CO2 feedbacks'
        CALL jules_print('imogen_check',jules_message)
      ELSE
        WRITE(jules_message,*) 'There are NO land CO2 feedbacks'
        CALL jules_print('imogen_check',jules_message)
      END IF
        
      IF (ocean_feed) THEN
        WRITE(jules_message,*) 'There are ocean CO2 feedbacks'
        CALL jules_print('imogen_check',jules_message)
      ELSE
        WRITE(jules_message,*) 'There are NO ocean CO2 feedbacks'
        CALL jules_print('imogen_check',jules_message)
      END IF
        
      ! Check for confirmed configurations:
      IF ( land_feed_co2 .AND. ocean_feed .AND. ( .NOT. land_feed_ch4) .AND.  &
           ( .NOT. wgen) ) THEN
        confirm_flag = .TRUE.
      END IF

    END IF ! end c_emissions true
  END IF ! End include_co2 true
  
  ! Checks for non-co2 components
  IF (include_non_co2) THEN
    WRITE(jules_message,*) 'Run is with non-CO2 radiative forcing'
    CALL jules_print('imogen_check',jules_message)
    WRITE(jules_message,*) 'User provides file of Non-CO2 radiative forcing'
    CALL jules_print('imogen_check',jules_message)

    ! Write if CH4 feedbacks are turned on
    IF (land_feed_ch4) THEN
      WRITE(jules_message,*) 'There are land CH4 feedbacks'
      CALL jules_print('imogen_check',jules_message)
    ELSE
      WRITE(jules_message,*) 'There are NO land CH4 feedbacks'
      CALL jules_print('imogen_check',jules_message)
    END IF
    
    ! Check for confirmed configurations
    IF ( include_co2 .AND. ocean_feed .AND. ( land_feed_ch4) .AND.            &
        ( .NOT. wgen) ) THEN
      confirm_flag = .TRUE.
    END IF
  ELSE
    ! Failsafe for if CH4 feedbacks are on without non-CO2 Radiative Forcing
    IF (land_feed_ch4) THEN
      WRITE(jules_message,*) 'Land CH4 feedbacks require non_CO2 radiative'// &
                             ' forcing to use as baseline'
      CALL jules_print('imogen_check',jules_message)
      check_flag = .FALSE.
    END IF
  END IF

END IF ! end if analogue true.

! Write if weather generator is turned on
IF ( .NOT. wgen) THEN
  WRITE(jules_message,*) 'Weather generator is switched off'
  CALL jules_print('imogen_check',jules_message)
ELSE
  WRITE(jules_message,*) 'Weather generator is switched on'
  CALL jules_print('imogen_check',jules_message)
END IF

IF ( ( .NOT. anom) .AND. anlg ) THEN
  WRITE(jules_message,*) ' The imogen analogue model requires anomalies.' //  &
        ' This combination should only be used for spin-up.'//                &
        ' Please ensure that flags are compatible with transient run.'
  CALL log_warn('imogen_check',jules_message)
END IF

IF ( .NOT. check_flag)                                                        &
  CALL log_fatal("IMOGEN_CHECK",                                              &
                 'Combination not yet allowed')

IF ( .NOT. confirm_flag)                                                      &
  CALL log_warn("IMOGEN_CHECK",                                               &
    'This combination of flags has not yet been tested, proceed with caution')

! Finilise imogen_check message with line of stars
WRITE(jules_message,*) '***********************************************'
CALL jules_print('imogen_check',jules_message)

RETURN

END SUBROUTINE imogen_check
