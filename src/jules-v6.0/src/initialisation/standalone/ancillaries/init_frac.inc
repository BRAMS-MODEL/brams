#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_frac()

USE data_cube_mod, ONLY: data_cube, cube_create

USE mem_brams_jules, ONLY: fracB,nxB,nyB !JWS

USE jules_surface_types_mod,       ONLY : ntype

USE io_constants, ONLY: max_file_name_len, max_sdf_name_len, namelist_unit

USE string_utils_mod, ONLY: to_string

USE jules_vegetation_mod, ONLY: l_veg_compete
  
USE dump_mod, ONLY: ancil_dump_read

USE errormessagelength_mod, ONLY: errormessagelength

USE model_interface_mod, ONLY: get_var_id, populate_var

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in the tile fractions and checks for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables

TYPE(data_cube) :: DATA


INTEGER :: error,i,j,k,l  ! Error indicator

!-----------------------------------------------------------------------------
! Definition of the jules_frac namelist
!-----------------------------------------------------------------------------
LOGICAL :: read_from_dump
CHARACTER(LEN=max_file_name_len) :: FILE
CHARACTER(LEN=max_sdf_name_len) :: frac_name
CHARACTER(LEN=errormessagelength) :: iomessage
NAMELIST  / jules_frac/ read_from_dump, FILE, frac_name


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
read_from_dump = .FALSE.
FILE=''
frac_name      = ''

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
CALL log_info("init_frac", "Reading JULES_FRAC namelist...")

!First, we read the namelist
READ(namelist_unit, NML = jules_frac, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_frac",                                                 &
                "Error reading namelist JULES_FRAC " //                      &
                "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                TRIM(iomessage) // ")")

ancil_dump_read%frac = read_from_dump

IF ( .NOT. ancil_dump_read%frac) THEN !we read from the ancil file
  !-------------------------------------------------------------------------
  ! Set frac using namelist values
  !-------------------------------------------------------------------------
  ! Frac is prognostic if competing veg is on, so must be read as an initial
  ! condition
  IF ( l_veg_compete ) THEN
    CALL log_info("init_frac",                                                &
                  "Competing vegetation is enabled - frac will be read " //   &
                  "as an initial condition")
    RETURN
  END IF

  !Check that file name was provided
IF ( LEN_TRIM(FILE) == 0 )                                                  &
  CALL log_fatal("init_frac", "No file name provided")

!JWS CALL log_info("init_frac",                                                  &
!JWS              "Reading tile fractions from file " // TRIM(FILE))


open(unit=66,file='jules.log',position='append',status='old',action='write')
write(unit=66,fmt='(A)') 'Reading tile fractions from BRAMS'
close(unit=66)

!JWS CALL fill_variables_from_file(FILE, (/ 'frac' /), (/ frac_name /))

DATA = cube_create((/nxB, nyB, ntype/))

   DATA%SHAPE(1)=nxB
   DATA%SHAPE(2)=nyB
   DATA%SHAPE(3)=ntype
   l=0
   do k=1,ntype
      do j=1,nyB
         do i=1,nxB
            l=l+1
            DATA%values(l)=fracB(i,j,k)
         enddo
      enddo
   enddo

CALL populate_var(get_var_id('frac'), DATA)

ELSE !We read from the dump file
  CALL log_info("init_frac",                                                  &
                "frac will be read from the dump file.  " //                  &
                "Namelist values ignored")

END IF !.NOT. ancil_dump_read%frac

RETURN

END SUBROUTINE init_frac
#endif
