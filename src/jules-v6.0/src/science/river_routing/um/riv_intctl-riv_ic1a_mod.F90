#if defined(UM_JULES)
!Huw Lewis (MO), Jan 2015
!DEPRECATED CODE
!This code was transferred from the UM repository at UM vn9.2 / JULES vn 4.1.
!Future developments will supercede these subroutines, and as such they
!should be considered deprecated. They will be retained in the codebase to
!maintain backward compatibility with functionality prior to
!UM vn10.0 / JULES vn 4.2, until such time as they become redundant.
!
!This code is not to be used with soil tiling
!It is possible, with significant effort to make the required modifications,
!but this routine is marked as deprecated.
!
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: River Routing
MODULE riv_intctl_mod_1A

USE riv_rout_mod_1A, ONLY: riv_rout_1A

USE umPrintMgr, ONLY:                                                         &
    umPrint,                                                                  &
    umMessage

USE errormessagelength_mod, ONLY: errormessagelength

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RIV_INTCTL_MOD_1A'

CONTAINS

SUBROUTINE riv_intctl_1a(                                                     &
xpa, xua, xva, ypa, yua, yva,                                                 &
g_p_field, g_r_field, n_proc,                                                 &
gather_pe_trip,land_points,land_index,                                        &
invert_atmos, row_length, rows,                                               &
global_row_length, global_rows,                                               &
river_row_length, river_rows,                                                 &
global_river_row_length, global_river_rows,                                   &
flandg, riv_step, riv_vel, riv_mcoef,                                         &
trivdir, trivseq, twatstor, riverout_rgrid, a_boxareas,                       &
delta_phi, first_entry,                                                       &
r_area, slope, flowobs1, r_inext, r_jnext,                                    &
r_land,  substore, surfstore, flowin, bflowin,                                &
! IN/OUT accumulated runoff
tot_surf_runoff_gb, tot_sub_runoff_gb,                                        &
! OUT
box_outflow, box_inflow, riverout_atmos,                                      &
! Add new arguments for inland basin outflow
! OUT INLAND BASINS
inlandout_atmos,inlandout_riv,                                                &
! Required for soil moisture correction for water conservation
dsm_levels,acc_lake_evap,smvcst,smvcwt,smcl,sthu)

! Purpose:
! New Control routine for River routing for Global Model.
! Parallel River routing
!

! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!-----------------------------------------------------------------

USE water_constants_mod, ONLY: rho_water
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE global_2d_sums_mod, ONLY: global_2d_sums
USE ereport_mod, ONLY: ereport
USE UM_ParVars,   ONLY: lasize, glsize, gc_all_proc_group
USE UM_ParParams, ONLY: halo_type_no_halo
USE Field_Types,  ONLY: fld_type_p
USE all_gather_field_mod, ONLY: all_gather_field
USE jules_soil_mod, ONLY: dzsoil
IMPLICIT NONE



INTEGER, INTENT(IN) ::                                                        &
 row_length                                                                   &
                           ! NO. OF COLUMNS IN ATMOSPHERE
, rows                                                                        &
                           ! NO. OF ROWS IN ATMOSPHERE
, global_row_length                                                           &
                           ! number of points on a row
, global_rows                                                                 &
                           ! NUMBER OF global rows
, land_points                                                                 &
                           ! number of landpoints
, river_row_length                                                            &
                           ! no. of columns in river grid
, river_rows                                                                  &
                           ! no. of rows in river grid
, global_river_row_length                                                     &
                           ! global river row length
, global_river_rows                                                           &
                           ! GLOBAL river rows
, gather_pe_trip                                                              &
                           ! pe River routing to be run on
, n_proc                                                                      &
                           ! Total number of processors
, dsm_levels               ! no. of deep soil moisture levels

INTEGER, INTENT(IN) ::                                                        &
 g_p_field                                                                    &
                            ! IN size of global ATMOS field
, g_r_field                                                                   &
                            ! IN Size of global river field
, land_index (land_points)  ! IN index of land to global points

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
 tot_surf_runoff_gb(land_points)                                              &
                             !IN Surf RUNOFF on land pts(KG/M2/S)
,tot_sub_runoff_gb(land_points)                                               &
                             ! IN Subsurf.RUNOFF (KG/M2/S)

! Water conservation
! Remove global lake evaporation from soil moisture store
     , smvcwt(land_points)                                                    &
!                            ! IN Volumetric wilting point

     , smvcst(land_points)                                                    &
!                            ! IN Volumetric saturation point
     ,xua(0:global_row_length)                                                &
                                  ! IN Atmosphere UV longitude coords
     ,yua(global_rows)                                                        &
                                  ! IN Atmosphere latitude coords
     ,xpa(global_row_length+1)                                                &
                                  ! IN Atmosphere longitude coords
     ,ypa(global_rows)                                                        &
                                  ! IN Atmosphere latitude coords
     ,xva(global_row_length+1)                                                &
                                  ! IN Atmosphere longitude coords
     ,yva(0:global_rows)                                                      &
                                  ! IN Atmosphere latitude coords
     ,a_boxareas(row_length, rows)                                            &
                                  !IN ATMOS gridbox areas
     ,flandg(row_length, rows)
                                  ! IN Land fraction on global field.

     ! ********* NOT USED *********!
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
r_area(global_row_length, global_rows)                                        &
                            ! IN accumalated areas file
,slope(global_row_length, global_rows)                                        &
                            ! IN slopes
,r_inext(global_row_length, global_rows)                                      &
                            ! IN X-cordinate of downstream grid pt
,r_jnext(global_row_length, global_rows)                                      &
                            ! IN Y-cordinate of downstream grid pt
,r_land(global_row_length, global_rows)                                       &
                         ! IN land/river depends on value of a_thresh
,flowobs1(global_row_length, global_rows)                                     &
                         ! IN initialisation for flows
,substore(global_row_length, global_rows)                                     &
                         ! IN routing subsurface store
,surfstore(global_row_length, global_rows)                                    &
                         ! IN routing surface store
,flowin(global_row_length, global_rows)                                       &
                         ! IN surface lateral inflow
,bflowin (global_row_length, global_rows)                                     &
                         ! IN subsurface lateral inflow
,delta_phi
                         ! IN RCM gridsize (radians)

! ********* NOT USED *********!

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
trivdir(river_row_length, river_rows)                                         &
                                         ! IN river direction
,trivseq(river_row_length, river_rows)                                        &
                                         ! IN river sequence
,riv_vel                                                                      &
                            ! IN river velocity
,riv_mcoef                                                                    &
                            ! IN meandering coefficient
,riv_step
                            ! IN river timestep (secs)

REAL(KIND=real_jlslsm), INTENT(INOUT) :: twatstor(river_row_length, river_rows)
! water store(Kg)

LOGICAL, INTENT(IN) ::                                                        &
 invert_atmos               ! IN True if ATMOS fields are S->N
!                                 ! for regridding runoff from ATMOS.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
riverout_atmos(row_length,rows)                                               &
     ! river flow out from each  gridbox(KG/m2/S)
,riverout_rgrid(river_row_length, river_rows)                                 &
           ! river flow out from TRIP grid to ocean (Kg/s)

,box_outflow(river_row_length, river_rows)                                    &
     ! gridbox outflow river grid (Kg/s)

,box_inflow(river_row_length, river_rows)                                     &
     ! gridbox runoff river grid(Kg/s)

! Declare new variables for inland basin outflow

     ,inlandout_riv(river_row_length,river_rows)                              &
           !TRIP OUTFLOW FROM INLAND BASINS ON TRIP GRID Kg/s

     ,inlandout_atmos(row_length,rows)
           ! TRIP OUTFLOW FROM  INLAND BASINS ON atmos GRID Kg/m2/s


! Logical to detect first entry to routing code

LOGICAL, INTENT(IN) ::  first_entry

! local variables

INTEGER ::                                                                    &
 i,j,l

INTEGER ::                                                                    &
 info                                                                         &
                          ! Return code from MPP
, icode                   ! Return code :0 Normal Exit :>0 Error

CHARACTER(LEN=errormessagelength) ::                                          &
 CMessage         ! Error message if return code >0
CHARACTER(LEN=*) :: RoutineName
PARAMETER ( RoutineName='RIV_INTCTL_1A')

LOGICAL ::                                                                    &
 invert_trip                                                                  &
                         ! TRUE WHEN ROW INVERSION IS REQ
, regrid                                                                      &
                         ! TRUE if TRIP grid different to ATMOS
, cyclic_trip                                                                 &
                         ! TRUE WHEN THE TRIP MODEL HAS CYCLIC
, global_trip            ! TRUE WHEN TRIP GRID SURFACE IS SPHER
PARAMETER(invert_trip = .FALSE.,cyclic_trip = .TRUE.,                         &
global_trip = .TRUE.,regrid = .TRUE.)

REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  sthu(land_points)                                                           &
                        ! Frozen soil moisture content of
!                            !    bottom layer as a frac of saturation
     , smcl(land_points)
                             ! Total soil moisture contents
!                            !      of bottom layer (kg/m2).

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
 acc_lake_evap(row_length,rows)
!                                    ! Accumulated lake evap over
!                                    ! river routing timestep (Kg/m2)

REAL(KIND=real_jlslsm) ::                                                     &
 smvcst_g(row_length,rows)                                                    &
!                                    ! SMVCST on global grid
     ,smvcwt_g(row_length,rows)                                               &
!                                    ! SMVCWT on global grid
     ,smcl_dsm_g(row_length,rows)                                             &
!                                    ! SMCL on DSM_LEVELS LAYER
!                                    ! on global grid
     ,sthu_dsm_g(row_length,rows)                                             &
!                                    ! STHU on DSM_LEVELS LAYER
!                                    ! on global grid
     ,surf_runoffin(row_length,rows)                                          &
                                     ! TOTAL RATE OF RUNOFF (KG/M2/S)
     ,sub_runoffin(row_length,rows)
!                                    ! TOTAL RATE OF RUNOFF (KG/M2/S)

REAL(KIND=real_jlslsm) ::                                                     &
 acc_lake_evap_avg
 ! Global daily accumulated total lake evaporation

REAL(KIND=real_jlslsm), DIMENSION(:,:), SAVE, ALLOCATABLE :: land_frac
! land mask for global regridding weights calculation

INTEGER :: my_comm, ierror
REAL(KIND=real_jlslsm) :: acc_lake_evap_avg_g = 0, tot_landarea_g = 0
REAL(KIND=real_jlslsm) :: sum_r(row_length, rows, 2), sum_all(2)

LOGICAL, SAVE :: first = .TRUE.
! use this rather than first entry because of continuation runs
! and the need to recalculate concerns and gather land fractions
! on first fun

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! gather the TRIP variables to PE0 and call the TRIP river routing
! 1. Gather coupling fields from distributed processors onto a single
!    processor for input to river routing routines.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,                &
                        zhook_handle)
info = 0

sum_r = 0.0
sum_all = 0.0

! *******************************************************************
DO j = 1,rows
  DO i = 1,row_length
    surf_runoffin(i,j) = 0.0
    sub_runoffin(i,j) = 0.0
  END DO
END DO
DO j = 1,rows
  DO i = 1,row_length
    smvcst_g(i,j) = 0.0
    smvcwt_g(i,j) = 0.0
    smcl_dsm_g(i,j) = 0.0
    sthu_dsm_g(i,j) = 0.0
  END DO
END DO

! Copy land points output back to full fields array.
DO l = 1, land_points
  j=(land_index(l) - 1) / row_length + 1
  i = land_index(l) - (j-1) * row_length
  surf_runoffin(i,j) = tot_surf_runoff_gb(l)
  sub_runoffin(i,j)  = tot_sub_runoff_gb(l)

  !Note for soil tiling- these varaibles cannot be arithmetically
  !averaged, making conversion to be compatible with nsoilt > 1
  !rather difficult
  smvcst_g(i,j)      = smvcst(l)
  smvcwt_g(i,j)      = smvcwt(l)
  smcl_dsm_g(i,j)    = smcl(l)
  sthu_dsm_g(i,j)    = sthu(l)
END DO

IF (first) THEN
  ALLOCATE(land_frac(global_row_length, global_rows))

  CALL all_gather_field(flandg, land_frac,                                    &
       lasize(1,fld_type_p,halo_type_no_halo),                                &
       lasize(2,fld_type_p,halo_type_no_halo),                                &
       glsize(1,fld_type_p), glsize(2,fld_type_p),                            &
       fld_type_p,halo_type_no_halo,                                          &
       gc_all_proc_group,info,cmessage)
  IF (info /= 0) THEN      ! Check return code
    cmessage='ATMPHB2 : ERROR in gather of flandg'
    icode = 5

    CALL Ereport(RoutineName,icode,Cmessage)
  END IF

END IF

! first get number of river grid points on this proc

! Calculate the water required to be removed from the soil to
! take account of lake evaporation

! Calculate the global total lake evaporation accumulated over a day:

acc_lake_evap_avg = 0.0

DO j = 1, rows
  DO i = 1, row_length

    IF (sthu_dsm_g(i,j) * smvcst_g(i,j) >  smvcwt_g(i,j))                     &
          sum_r(i,j,1) =  a_boxareas(i,j) * flandg(i,j)

    sum_r(i,j,2)= acc_lake_evap(i,j) * a_boxareas(i,j) * flandg(i,j)

  END DO
END DO

! calculate global landarea and acc_lake_evap_avg sums
CALL global_2d_sums(sum_r, row_length, rows, 0, 0, 2, sum_all)

tot_landarea_g = sum_all(1)
acc_lake_evap_avg_g = sum_all(2)

acc_lake_evap_avg = acc_lake_evap_avg_g / tot_landarea_g

! Apply water correction:

DO j = 1, rows
  DO i = 1, row_length
    IF (sthu_dsm_g(i,j) * smvcst_g(i,j) >  smvcwt_g(i,j)) THEN

      smcl_dsm_g(i,j) = smcl_dsm_g(i,j) - acc_lake_evap_avg
      sthu_dsm_g(i,j) = sthu_dsm_g(i,j) - acc_lake_evap_avg                   &
            / (smvcst_g(i,j) * rho_water * dzsoil(dsm_levels))

      IF (sthu_dsm_g(i,j) <  0.0) THEN
        CALL umPrint('*warning layer 4 unfr soil moisture frac<0' ,           &
            src='riv_intctl-riv_ic1a')
        WRITE(umMessage,'(2(F8.2))')STHU_DSM_g(i,j),STHU_DSM_g(i,j)           &
             +acc_lake_evap_avg
        CALL umPrint(umMessage,src='riv_intctl-riv_ic1a')
      END IF

      IF (smcl_dsm_g(i,j) <  0.0) THEN
        CALL umPrint('*warning layer 4 unfrozen soil moisture frac<0' ,       &
            src='riv_intctl-riv_ic1a')
        WRITE(umMessage,'(2(F8.2))') smcl_DSM_g(i,j), smcl_DSM_g(i,j)         &
             +acc_lake_evap_avg
        CALL umPrint(umMessage,src='riv_intctl-riv_ic1a')
      END IF
    END IF
  END DO
END DO



CALL riv_rout_1a(global_row_length, global_rows, surf_runoffin,               &
      sub_runoffin, row_length, rows, invert_atmos, xua, yva,                 &
      river_row_length, river_rows, flandg, regrid, riv_vel,                  &
      riv_mcoef, riv_step, cyclic_trip, a_boxareas, global_trip,              &
      invert_trip, twatstor, riverout_rgrid, trivdir, trivseq,                &
      box_outflow, box_inflow, riverout_atmos, inlandout_riv,                 &
      inlandout_atmos, global_river_row_length, global_river_rows,            &
      land_frac, first, cmessage)

WHERE (flandg == 1.0) riverout_atmos = 0.0

! Copy full fields array back to land points:
DO l = 1, land_points
  j=(land_index(l) - 1) / row_length + 1
  i = land_index(l) - (j-1) * row_length
  smcl(l) = smcl_dsm_g(i,j)
  sthu(l) = sthu_dsm_g(i,j)
END DO

first = .FALSE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,               &
                        zhook_handle)
RETURN
END SUBROUTINE riv_intctl_1a
END MODULE riv_intctl_mod_1A
#endif
