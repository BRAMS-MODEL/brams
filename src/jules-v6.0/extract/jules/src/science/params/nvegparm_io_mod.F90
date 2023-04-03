! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module contains variables used for reading in nvegparm data
! and initialisations

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.


MODULE nvegparm_io

USE max_dimensions, ONLY:                                                     &
  nnvg_max
USE missing_data_mod, ONLY: rmdi
USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------
! Set up variables to use in IO (a fixed size version of each array
! in nvegparm that we want to initialise).
!-----------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                     &
  albsnc_nvg_io(nnvg_max) = rmdi,                                             &
  albsnf_nvgu_io(nnvg_max) = rmdi,                                            &
  albsnf_nvg_io(nnvg_max) = rmdi,                                             &
  albsnf_nvgl_io(nnvg_max) = rmdi,                                            &
  catch_nvg_io(nnvg_max) = rmdi,                                              &
  gs_nvg_io(nnvg_max) = rmdi,                                                 &
  infil_nvg_io(nnvg_max) = rmdi,                                              &
  z0_nvg_io(nnvg_max) = rmdi,                                                 &
  ch_nvg_io(nnvg_max) = rmdi,                                                 &
  vf_nvg_io(nnvg_max) = rmdi,                                                 &
  emis_nvg_io(nnvg_max) = rmdi,                                               &
  z0hm_nvg_io(nnvg_max) = rmdi,                                               &
  z0hm_classic_nvg_io(nnvg_max) = rmdi


!-----------------------------------------------------------------------
! Set up a namelist for reading and writing these arrays
!-----------------------------------------------------------------------
NAMELIST  / jules_nvegparm/                                                   &
                          albsnc_nvg_io,albsnf_nvgu_io,                       &
                          albsnf_nvg_io, albsnf_nvgl_io,                      &
                          catch_nvg_io,gs_nvg_io,infil_nvg_io,                &
                          z0_nvg_io,ch_nvg_io,vf_nvg_io,                      &
                          emis_nvg_io,z0hm_nvg_io,                            &
                          z0hm_classic_nvg_io

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NVEGPARM_IO'

CONTAINS

SUBROUTINE print_nlist_jules_nvegparm()
USE jules_print_mgr, ONLY: jules_print
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer

CALL jules_print('nvegparm_io',                                               &
    'Contents of namelist jules_nvegparm')

WRITE(lineBuffer,*)' albsnc_nvg_io = ',albsnc_nvg_io
CALL jules_print('nvegparm_io',lineBuffer)
WRITE(lineBuffer,*)' albsnf_nvg_io = ',albsnf_nvg_io
CALL jules_print('nvegparm_io',lineBuffer)
WRITE(lineBuffer,*)' catch_nvg_io = ',catch_nvg_io
CALL jules_print('nvegparm_io',lineBuffer)
WRITE(lineBuffer,*)' gs_nvg_io = ',gs_nvg_io
CALL jules_print('nvegparm_io',lineBuffer)
WRITE(lineBuffer,*)' infil_nvg_io = ',infil_nvg_io
CALL jules_print('nvegparm_io',lineBuffer)
WRITE(lineBuffer,*)' z0_nvg_io = ',z0_nvg_io
CALL jules_print('nvegparm_io',lineBuffer)
WRITE(lineBuffer,*)' ch_nvg_io = ',ch_nvg_io
CALL jules_print('nvegparm_io',lineBuffer)
WRITE(lineBuffer,*)' vf_nvg_io = ',vf_nvg_io
CALL jules_print('nvegparm_io',lineBuffer)
WRITE(lineBuffer,*)' emis_nvg_io = ',emis_nvg_io
CALL jules_print('nvegparm_io',lineBuffer)
WRITE(lineBuffer,*)' z0hm_nvg_io = ',z0hm_nvg_io
CALL jules_print('nvegparm_io',lineBuffer)

CALL jules_print('nvegparm_io',                                               &
    '- - - - - - end of namelist - - - - - -')

END SUBROUTINE print_nlist_jules_nvegparm

#if defined(UM_JULES)
SUBROUTINE read_nml_jules_nvegparm (unitnumber)

! Description:
!  read the JULES_NVEGPARM namelist

USE setup_namelist,   ONLY: setup_nml_type

USE check_iostat_mod, ONLY: check_iostat
USE UM_parcore,       ONLY: mype
USE errormessagelength_mod, ONLY: errormessagelength

USE parkind1,         ONLY: jprb, jpim
USE yomhook,          ONLY: lhook, dr_hook

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: unitnumber
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_JULES_NVEGPARM'
INTEGER(KIND=jpim), PARAMETER          :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER          :: zhook_out = 1

CHARACTER(LEN=errormessagelength) :: iomessage

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 1
INTEGER, PARAMETER :: n_real = 13 * nnvg_max

TYPE my_namelist
  SEQUENCE
  REAL(KIND=real_jlslsm) :: albsnc_nvg_io(nnvg_max)
  REAL(KIND=real_jlslsm) :: albsnf_nvgu_io(nnvg_max)
  REAL(KIND=real_jlslsm) :: albsnf_nvg_io(nnvg_max)
  REAL(KIND=real_jlslsm) :: albsnf_nvgl_io(nnvg_max)
  REAL(KIND=real_jlslsm) :: catch_nvg_io(nnvg_max)
  REAL(KIND=real_jlslsm) :: gs_nvg_io(nnvg_max)
  REAL(KIND=real_jlslsm) :: infil_nvg_io(nnvg_max)
  REAL(KIND=real_jlslsm) :: z0_nvg_io(nnvg_max)
  REAL(KIND=real_jlslsm) :: ch_nvg_io(nnvg_max)
  REAL(KIND=real_jlslsm) :: vf_nvg_io(nnvg_max)
  REAL(KIND=real_jlslsm) :: emis_nvg_io(nnvg_max)
  REAL(KIND=real_jlslsm) :: z0hm_nvg_io(nnvg_max)
  REAL(KIND=real_jlslsm) :: z0hm_classic_nvg_io(nnvg_max)
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_real_in = n_real)

IF (mype == 0) THEN

  READ (UNIT = unitnumber, NML = jules_nvegparm, IOSTAT = errorstatus,        &
        IOMSG = iomessage)
  CALL check_iostat(errorstatus, "namelist jules_nvegparm", iomessage)

  my_nml % albsnc_nvg_io   = albsnc_nvg_io
  my_nml % albsnf_nvgu_io  = albsnf_nvgu_io
  my_nml % albsnf_nvg_io   = albsnf_nvg_io
  my_nml % albsnf_nvgl_io  = albsnf_nvgl_io
  my_nml % catch_nvg_io    = catch_nvg_io
  my_nml % gs_nvg_io       = gs_nvg_io
  my_nml % infil_nvg_io    = infil_nvg_io
  my_nml % z0_nvg_io       = z0_nvg_io
  my_nml % ch_nvg_io       = ch_nvg_io
  my_nml % vf_nvg_io       = vf_nvg_io
  my_nml % emis_nvg_io     = emis_nvg_io
  my_nml % z0hm_nvg_io     = z0hm_nvg_io
  my_nml % z0hm_classic_nvg_io = z0hm_classic_nvg_io
END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  albsnc_nvg_io   = my_nml % albsnc_nvg_io
  albsnf_nvgu_io  = my_nml % albsnf_nvgu_io
  albsnf_nvg_io   = my_nml % albsnf_nvg_io
  albsnf_nvgl_io  = my_nml % albsnf_nvgl_io
  catch_nvg_io    = my_nml % catch_nvg_io
  gs_nvg_io       = my_nml % gs_nvg_io
  infil_nvg_io    = my_nml % infil_nvg_io
  z0_nvg_io       = my_nml % z0_nvg_io
  ch_nvg_io       = my_nml % ch_nvg_io
  vf_nvg_io       = my_nml % vf_nvg_io
  emis_nvg_io     = my_nml % emis_nvg_io
  z0hm_nvg_io     = my_nml % z0hm_nvg_io
  z0hm_classic_nvg_io = my_nml % z0hm_classic_nvg_io
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_nml_jules_nvegparm
#endif

END MODULE nvegparm_io
