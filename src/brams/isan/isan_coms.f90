!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


Module isan_coms

  use ModNamelistFile, only: namelistFile

  use grid_dims, only : maxfiles !INTENT(IN)

  !---------------------------------------------------------------------------
  !    Configuration COMMON blocks for RAMS isentropic data analysis package.
  !---------------------------------------------------------------------------
  !MAXPR         Maximum number of vertical levels that can be used in 
  !                 the pressure data.
  !MAXISN        Maximum number of vertical levels that can be used in 
  !                 the isentropic analysis.
  !MAXX          Maximum number of grid points (in the x or east-west direction) 
  !                 in any of the RAMS or pressure grids.
  !MAXY          Maximum number of grid points (in the y or north-south direction) 
  !                 in any of the RAMS or pressure grids.
  !MAXTIMES      Maximum number of data analysis times that can be processed 
  !                 in a single run.
  !MAXAGRDS      Maximum number of RAMS grids that can have varfiles generated.
  !MAXSIGZ       Maximum number of vertical levels that can be used 
  !                 in the  _z analysis.
  !MAXLEV	       Maximum number of levels in an input rawinsonde.
  !MAXSNAME      Maximum number of input observations
  !MAXISFILES    Maximum number of input data times
  !------------------------------------------------------------------------------------

  integer, parameter :: maxpr=100 ,maxisn=100 ,maxx=1000 ,maxy=1000  &
       ,maxtimes=100 ,maxagrds=10    ,maxsigz=100  &
       ,maxlev=500   ,maxsname=10000 ,maxisfiles=maxfiles !1000
  !---------------------------------------------------------------------------
  integer :: ioflgisz,ioflgvar,natime,iszstage,ivrstage,iyear,imonth,idate  &
       ,ihour,isan_inc,i1st_flg,iupa_flg,isfc_flg
  !---------------------------------------------------------------------------

  include "files.h"

  character(len=f_name_length)  :: innpr,inrawi,insrfce

  character(len=256) :: varpfx, iapr, iarawi, iasrfce 
  ! Modif. by ALF

  character(len=8)   :: pdata,guess1st

  !---------------------------------------------------------------------------
  !     Input pressure file header
  !---------------------------------------------------------------------------
  integer :: marker,isversion,iyy,imm,idd,ihh,itinc,inproj,ivertcoord
  real    :: xnelat,xnelon,cntlat,cntlon,secondlat

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Input pressure data memory

  real, allocatable, dimension(:,:,:) :: p_u,p_v,p_t,p_z,p_r,p_ur,p_vr
  real, allocatable, dimension(:,:)   :: p_lat,p_lon
  real, allocatable, dimension(:,:)   :: p_slp,p_sfp,p_sft,p_snow,p_sst
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Polar-stereo/pressure grid memory

  real, allocatable, dimension(:,:,:) :: pp_u,pp_v,pp_t,pp_z,pp_r
  real, allocatable, dimension(:,:)   :: pp_sglob
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Polar-stereo/isentropic grid memory

  real, allocatable, dimension(:,:,:) :: pi_u,pi_v,pi_p,pi_s,pi_r  &
       ,pi_scra,pi_scrb
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Polar-stereo/sigma-z grid memory
  !                         :: 
  real, allocatable, dimension(:,:,:) :: ps_u,ps_v,ps_p,ps_t,ps_r  &
       ,ps_scra,ps_scrb
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Polar-stereo/surface grid memory

  real, allocatable, dimension(:,:) :: rs_u,rs_v,rs_p,rs_t,rs_r,rs_s  &
       ,rs_top,rs_qual  &
       ,rs_slp,rs_sfp,rs_sft,rs_snow,rs_sst
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Data type to replace A array memory use in ISAN. 

  type isan_grids
     real, pointer, dimension(:,:,:) :: rr_u,rr_v,rr_t,rr_p,rr_r  
     real, pointer, dimension(:,:,:) :: rr_ug,rr_vg,rr_tg,rr_pg,rr_rg  
     real, pointer, dimension(:,:,:) :: rr_pi0,rr_th0,rr_dn0,rr_dn0u,rr_dn0v  
     real, pointer, dimension(:,:)   :: rr_slp,rr_sfp,rr_sft,rr_snow,rr_sst  
  end type isan_grids
  real, allocatable, dimension(:)    :: rr_scr1,rr_scr2,rr_vt2da

  type (isan_grids)                  :: is_grids(maxagrds)

  real, dimension(maxsigz,maxagrds)  :: piref,thref,dnref,rtref

  integer                            :: maxix,maxiy,maxiz

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Input observation data memory
  !
  real, allocatable, dimension(:,:)  :: up_uz,up_vz,up_ur,up_vr,up_zz  &
       ,up_p,up_t,up_z,up_r
  real, allocatable, dimension(:)    :: up_lat,up_lon,up_top
  real, allocatable, dimension(:,:)  :: up_topg
  integer, allocatable, dimension(:) :: up_lp, up_lz
  character(len=8), allocatable, dimension(:) :: up_chstid

  real, allocatable, dimension(:)    :: sf_u,sf_v,sf_p,sf_t,sf_s,sf_r
  real, allocatable, dimension(:)    :: sf_ur,sf_vr
  real, allocatable, dimension(:)    :: sf_lat,sf_lon,sf_top,sf_scra
  character(len=8), allocatable, dimension(:) :: sf_chstid
  character(len=14), allocatable, dimension(:) :: sf_date
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Upper air-isentropic/sigma-z memory
  !
  real, allocatable, dimension(:,:)  :: upi_u,upi_v,upi_p,upi_s,upi_r
  real, allocatable, dimension(:,:)  :: ups_u,ups_v,ups_p,ups_t,ups_r
  !---------------------------------------------------------------------------

  integer                          :: npdates
  integer, dimension(maxisfiles,5) :: iproc_flag
  !---------------------------------------------------------------------------
  character(len=f_name_length), dimension(maxisfiles,6) :: iproc_names
  !---------------------------------------------------------------------------
  character(len=14), dimension(maxisfiles) :: iproc_dates
  !---------------------------------------------------------------------------
  integer                   :: nprx,npry,nprz,idatelin,iglobew,iglobs,iglobn
  integer, dimension(maxpr) :: levpr
  real                      :: xswlon,xswlat,gdatdx,gdatdy
  real, dimension(maxpr)    :: pnpr
  !---------------------------------------------------------------------------
  integer                   :: nsta,nssfc,notsta  &
       ,maxsta,maxsfc,iobswin
  real                      :: stasep
  real, dimension(maxagrds) :: wvlnth,respon,swvlnth
  !---------------------------------------------------------------------------
  character(len=8), dimension(50)       :: notid
  !---------------------------------------------------------------------------
  integer                    :: nisx,nisy,nisn,interp,igridfl,nigrids,nsigz  &
       ,nfeedvar
  integer, dimension(maxisn) :: levth
  real                       :: gobsep,gobrad,topsigz,hybbot,hybtop,sfcinf  &
       ,sigzwt
  real, dimension(maxsigz)   :: sigz
  real, dimension(maxagrds)  :: gridwt
  !---------------------------------------------------------------------------

  integer  :: icFileType
  character(len=256) :: icPrefix
  !# File prefix
  character(len=32)  :: wind_u_varname   
  !# Name of U wind variable inside grib2 descriptions (m/s)
  character(len=32)  :: wind_v_varname   
  !# Name of V wind variable inside grib2 descriptions (m/s)
  character(len=32)  :: temperature_varname
  !# Name of temperature variable inside grib2 descriptions (k)
  character(len=32)  :: geo_varname    
  !# Name of geopotential variable inside grib2 descriptions (m)
  character(len=32)  :: ur_varname     
  !# Name of humidity variable inside grib2 descriptions (%)
  real  :: initial_latitude
  !# Initial latitude (-90 to 90)
  real  :: final_latitude 
  !# Final latitude (-90 to 90)
  real  :: initial_longitude
  !# Initial longitude (0 to 360)
  real  :: final_longitude
  !# Final longitude (0 to 360)
  integer  :: z_max_level
  !# Maximum number of levels to be processed
  real,dimension(5) :: dlimit
  !# Inferior valid values for each variable (U,V,T,Geo,RU)
  real,dimension(5) :: ulimit
  !# Superior valid values for each variable (U,V,T,Geo,RU)
  real,dimension(5) :: scale_factor
  !# factor to multiplie each field
  integer :: ccGradsWrite
  !# Is to write grads with IC
  character(len=256) :: icGradsPrefix
  !# Grads with IC

  !# Flag to indicate if grads must be write
  logical, allocatable :: mask(:,:)
  character (len=30), allocatable :: slevs(:)
  real,allocatable :: pressLevs(:)
  integer :: nxGrib,nyGrib
  integer :: nprz_grib2
  real :: levpr_grib2(maxpr)

  type nct 
    real, allocatable  :: valuesRead(:,:,:) !i,j,k
  end type nct
  type (nct) :: ncVar(5)

contains

  subroutine StoreNamelistFileAtIsan_coms(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    gobrad = oneNamelistFile%gobrad
    gobsep = oneNamelistFile%gobsep
    gridwt = oneNamelistFile%gridwt
    guess1st = oneNamelistFile%guess1st
    hybbot = oneNamelistFile%hybbot
    hybtop = oneNamelistFile%hybtop
    i1st_flg = oneNamelistFile%i1st_flg
    iapr = oneNamelistFile%iapr
    iarawi = oneNamelistFile%iarawi
    iasrfce = oneNamelistFile%iasrfce
    igridfl = oneNamelistFile%igridfl
    iobswin = oneNamelistFile%iobswin
    ioflgisz = oneNamelistFile%ioflgisz
    ioflgvar = oneNamelistFile%ioflgvar
    isan_inc = oneNamelistFile%isan_inc
    isfc_flg = oneNamelistFile%isfc_flg
    iszstage = oneNamelistFile%iszstage
    iupa_flg = oneNamelistFile%iupa_flg
    ivrstage = oneNamelistFile%ivrstage
    levth = oneNamelistFile%levth
    maxsfc = oneNamelistFile%maxsfc
    maxsta = oneNamelistFile%maxsta
    nfeedvar = oneNamelistFile%nfeedvar
    nigrids = oneNamelistFile%nigrids
    nisn = oneNamelistFile%nisn
    notid = oneNamelistFile%notid
    notsta = oneNamelistFile%notsta
    respon = oneNamelistFile%respon
    sfcinf = oneNamelistFile%sfcinf
    sigzwt = oneNamelistFile%sigzwt
    stasep = oneNamelistFile%stasep
    swvlnth = oneNamelistFile%swvlnth
    topsigz = oneNamelistFile%topsigz
    varpfx = oneNamelistFile%varpfx
    wvlnth = oneNamelistFile%wvlnth
    !
    icFileType = oneNamelistFile%icFileType
    icPrefix = oneNamelistFile%icPrefix
    wind_u_varname = oneNamelistFile%wind_u_varname   
    wind_v_varname = oneNamelistFile%wind_v_varname   
    temperature_varname = oneNamelistFile%temperature_varname
    geo_varname = oneNamelistFile%geo_varname    
    ur_varname = oneNamelistFile%ur_varname     
    initial_latitude = oneNamelistFile%initial_latitude
    final_latitude = oneNamelistFile%final_latitude 
    initial_longitude = oneNamelistFile%initial_longitude
    final_longitude = oneNamelistFile%final_longitude
    z_max_level = oneNamelistFile%z_max_level
    dlimit = oneNamelistFile%dlimit
    ulimit = oneNamelistFile%ulimit
    scale_factor = oneNamelistFile%scale_factor
    CCGradsWrite = oneNamelistFile%CCGradsWrite
    icGradsPrefix = oneNamelistFile%icGradsPrefix

  end subroutine StoreNamelistFileAtIsan_coms
End module isan_coms
