!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT**************************************

MODULE overbank_inundation_mod

!-----------------------------------------------------------------------------
! Description:
!   Contains river overbank inundation variables and switches
!
!  Code Owner: Please refer to ModuleLeaders.txt
!  This file belongs in section: Hydrology
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE missing_data_mod, ONLY: rmdi

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Declare variables that can be read from the jules_overbank namelist.
!-----------------------------------------------------------------------------
! Overbank inundation parameters
!-----------------------------------------------------------------------------

REAL(KIND=real_jlslsm) ::                                                     &
   riv_a = rmdi,                                                              &
        ! Parameter in Leopold & Maddock (1953: eqn1)
   riv_b = rmdi,                                                              &
        ! Parameter in Leopold & Maddock (1953: eqn1)
   riv_c = rmdi,                                                              &
        ! Parameter in Leopold & Maddock (1953: eqn2)
   riv_f = rmdi,                                                              &
        ! Parameter in Leopold & Maddock (1953: eqn2)
   coef_b = rmdi,                                                             &
        ! Parameter to calculate the bankflow discharge power-law
        ! relationship from "Flood Modeling, Prediction and Mitigation"
        ! by Şen 2018.
   exp_c = rmdi,                                                              &
        ! Parameter to calculate the bankflow discharge power-law
        ! relationship from "Flood Modeling, Prediction and Mitigation"
        ! by Şen 2018.
   ent_ratio = rmdi
        ! Entrenchment ratio (Rosgen 1994). Ratio of the width of
        ! flood-prone area to surface width of bankfull channel.
        ! The flood-prone area width is measured at the elevation that
        ! corresponds to twice the maximum depth of the bankfull channel.

LOGICAL ::                                                                    &
   l_riv_overbank = .FALSE.,                                                  &
                            ! Logical to control overbank inundation
   l_riv_hypsometry = .TRUE.,                                                 &
                            ! TRUE indicates to calculate inundated area
                            ! using hypsometry, which requires additional
                            ! ancillaries logn_mean and logn_stdev
                            ! FALSE indicates to use a simpler width scaling
                            ! method (generally only to be used for testing
                            ! when those ancillaries are not available).
   use_rosgen = .FALSE.
                            ! Modify floodplain width using the Rosgen
                            ! entrenchment ratio (n.b. only for local use)

!-----------------------------------------------------------------------------
! Declare variables that are not in the namelist.
!-----------------------------------------------------------------------------
! Array variables defined on land points updated in overbank inundation
!-----------------------------------------------------------------------------

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  frac_fplain_lp(:)
           ! Overbank inundation area as a fraction of gridcell area
           ! (on land grid)

!-----------------------------------------------------------------------------
! Array variables defined on full 2D rivers grid (as read in from ancillary)
!-----------------------------------------------------------------------------

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  logn_mean(:,:),                                                             &
           ! ln(mean(elevation - elev_min)) for each gridcell (on land grid)
  logn_stdev(:,:)
           ! ln(SD(elevation - elev_min)) for each gridcell (on land grid)

!-----------------------------------------------------------------------------
! Other array variables
!-----------------------------------------------------------------------------

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  logn_mean_rp(:),                                                            &
           ! ln(mean(elevation - elev_min)) for each gridcell (on rivers grid)
  logn_stdev_rp(:),                                                           &
           ! ln(SD(elevation - elev_min)) for each gridcell (on rivers grid)
  qbf(:),                                                                     &
           ! Bankfull discharge rate (m3/s)
  dbf(:),                                                                     &
           ! Channel depth when at bankfull discharge rate (m)
  wbf(:),                                                                     &
           ! Channel width when at bankfull discharge rate (m)
  frac_fplain_rp(:)
           ! Overbank inundation area as a fraction of gridcell area
           ! (on rivers grid)

!-----------------------------------------------------------------------------
! Single namelist definition for UM and standalone
!-----------------------------------------------------------------------------
NAMELIST  / jules_overbank/                                                   &
  l_riv_overbank, l_riv_hypsometry, use_rosgen,                               &
  riv_a, riv_b, riv_c, riv_f, coef_b, exp_c, ent_ratio

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='OVERBANK_INUNDATION_MOD'

CONTAINS

!#############################################################################

SUBROUTINE check_jules_overbank()

USE ereport_mod, ONLY: ereport

USE jules_rivers_mod, ONLY: l_rivers

!-----------------------------------------------------------------------------
! Description:
!   Checks JULES_OVERBANK switches for consistency
!
! Current Code Owner: Toby Marthews (CEH)
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECK_JULES_OVERBANK'

INTEGER :: errcode

! If overbank is not selected, there's nothing more to do here, so leave.
IF ( .NOT. l_riv_overbank ) RETURN

! Overbank code can only be used if rivers are modelled.
IF ( ( .NOT. l_rivers) .AND. l_riv_overbank) THEN
  errcode = 101
  CALL ereport(RoutineName, errcode,                                          &
               'l_riv_overbank=T requires l_rivers=T ')
END IF

!-----------------------------------------------------------------------------
! If l_riv_hypsometry = T, ensure that use_rosgen is F - this means we can
! test use_rosgen alone, without also needing to check that l_riv_hypsometry
! = F.
!-----------------------------------------------------------------------------
IF ( l_riv_hypsometry ) use_rosgen = .FALSE.

!-----------------------------------------------------------------------------
! Check that the required parameters were found.
! We just check that a value is given, not that it is reasonable!
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Check parameters that are required if l_riv_hypsometry = T or
! use_rosgen = T.
!----------------------------------------------------------------------------- 
IF ( l_riv_hypsometry .OR. use_rosgen ) THEN

  IF ( ABS( riv_c - rmdi ) < EPSILON(1.0) ) THEN
    errcode = 101
    CALL ereport(RoutineName, errcode,'No value found for riv_c.')
  END IF

  IF ( ABS( riv_f - rmdi ) < EPSILON(1.0) ) THEN
    errcode = 101
    CALL ereport(RoutineName, errcode,'No value found for riv_f.')
  END IF

END IF

!-----------------------------------------------------------------------------
! Check values that are only required if l_riv_hypsometry = F.
!-----------------------------------------------------------------------------
IF ( .NOT. l_riv_hypsometry ) THEN

  IF ( ABS( riv_a - rmdi ) < EPSILON(1.0) ) THEN
    errcode = 101
    CALL ereport(RoutineName, errcode,'No value found for riv_a.')
  END IF

  IF ( ABS( riv_b - rmdi ) < EPSILON(1.0) ) THEN
    errcode = 101
    CALL ereport(RoutineName, errcode,'No value found for riv_b.')
  END IF

  !---------------------------------------------------------------------------
  ! Check values that are only required if use_rosgen = T.
  !---------------------------------------------------------------------------
  IF ( use_rosgen ) THEN

    IF ( ABS( coef_b - rmdi ) < EPSILON(1.0) ) THEN
      errcode = 101
      CALL ereport(RoutineName, errcode,'No value found for coef_b.')
    END IF

    IF ( ABS( ent_ratio - rmdi ) < EPSILON(1.0) ) THEN
      errcode = 101
      CALL ereport(RoutineName, errcode,'No value found for ent_ratio.')
    END IF

    IF ( ABS( exp_c - rmdi ) < EPSILON(1.0) ) THEN
      errcode = 101
      CALL ereport(RoutineName, errcode,'No value found for exp_c.')
    END IF

  END IF  !  use_rosgen
  
END IF  !  .NOT. l_riv_hypsometry

END SUBROUTINE check_jules_overbank

!#############################################################################

SUBROUTINE print_nlist_jules_overbank()

USE jules_print_mgr, ONLY: jules_print

IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_JULES_OVERBANK'

CALL jules_print('overbank_inundation_mod',                                   &
   'Contents of namelist jules_overbank')

WRITE(lineBuffer,*)' l_riv_overbank = ',l_riv_overbank
CALL jules_print('jules_overbank',lineBuffer)
WRITE(lineBuffer,*)' l_riv_hypsometry = ',l_riv_hypsometry
CALL jules_print('jules_overbank',lineBuffer)
WRITE(lineBuffer,*)' use_rosgen = ',use_rosgen
CALL jules_print('jules_overbank',lineBuffer)
WRITE(lineBuffer,*)' riv_a = ',riv_a
CALL jules_print('jules_overbank',lineBuffer)
WRITE(lineBuffer,*)' riv_b = ',riv_b
CALL jules_print('jules_overbank',lineBuffer)
WRITE(lineBuffer,*)' riv_c = ',riv_c
CALL jules_print('jules_overbank',lineBuffer)
WRITE(lineBuffer,*)' riv_f = ',riv_f
CALL jules_print('jules_overbank',lineBuffer)
WRITE(lineBuffer,*)' coef_b = ',coef_b
CALL jules_print('jules_overbank',lineBuffer)
WRITE(lineBuffer,*)' exp_c = ',exp_c
CALL jules_print('jules_overbank',lineBuffer)
WRITE(lineBuffer,*)' ent_ratio = ',ent_ratio
CALL jules_print('jules_overbank',lineBuffer)

CALL jules_print('overbank_inundation_mod',                                   &
   '- - - - - - end of namelist - - - - - -')

END SUBROUTINE print_nlist_jules_overbank

!#############################################################################


END MODULE overbank_inundation_mod

!-----------------------------------------------------------------------------
