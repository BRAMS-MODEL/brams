!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_grid

  use ModNamelistFile, only: namelistFile

  use grid_dims

  implicit none

  include "i8.h"

  type grid_vars

     ! Variables to be dimensioned by (nxp,nyp)

     real, pointer :: topt(:,:)
     real, pointer :: topu(:,:)
     real, pointer :: topv(:,:)
     real, pointer :: topm(:,:)
     real, pointer :: topma(:,:)
     real, pointer :: topta(:,:)
     real, pointer :: rtgt(:,:)
     real, pointer :: rtgu(:,:)
     real, pointer :: rtgv(:,:)
     real, pointer :: rtgm(:,:)
     real, pointer :: f13t(:,:)
     real, pointer :: f13u(:,:)
     real, pointer :: f13v(:,:)
     real, pointer :: f13m(:,:)
     real, pointer :: f23t(:,:)
     real, pointer :: f23u(:,:)
     real, pointer :: f23v(:,:)
     real, pointer :: f23m(:,:)
     real, pointer :: dxt(:,:)
     real, pointer :: dxu(:,:)
     real, pointer :: dxv(:,:)
     real, pointer :: dxm(:,:)
     real, pointer :: dyt(:,:)
     real, pointer :: dyu(:,:)
     real, pointer :: dyv(:,:)
     real, pointer :: dym(:,:)
     real, pointer :: fmapt(:,:)
     real, pointer :: fmapu(:,:)
     real, pointer :: fmapv(:,:)
     real, pointer :: fmapm(:,:)
     real, pointer :: fmapti(:,:)
     real, pointer :: fmapui(:,:)
     real, pointer :: fmapvi(:,:)
     real, pointer :: fmapmi(:,:)
     real, pointer :: glat(:,:)
     real, pointer :: glon(:,:)
     real, pointer :: topzo(:,:)

     !  Variables for the ADAP coordinate

     real, pointer :: aru(:,:,:)
     real, pointer :: arv(:,:,:)
     real, pointer :: arw(:,:,:)
     real, pointer :: volu(:,:,:)
     real, pointer :: volv(:,:,:)
     real, pointer :: volw(:,:,:)
     real, pointer :: volt(:,:,:)
     real, pointer :: lpu(:,:)
     real, pointer :: lpv(:,:)
     real, pointer :: lpw(:,:)
  end type grid_vars


  type (grid_vars), allocatable :: grid_g(:)
  type (grid_vars), allocatable :: gridm_g(:)

  ! data on entire grid (not domain decomposed)
  ! topography on entire grid (not domain decomposed)

  type GlobalGridData
     real, pointer :: global_topta(:,:)
     real, pointer :: global_glat(:,:)      ! set by GridSetup
     real, pointer :: global_glon(:,:)      ! set by GridSetup
  end type GlobalGridData

  type(GlobalGridData), allocatable, target :: oneGlobalGridData(:)


  character(len=64) :: expnme            ! experiment name; from RAMSIN
  integer :: ngrids                      ! how many grids; from RAMSIN
  integer :: ngridsh
  integer :: nxtnest(maxgrds)            ! next coarser grid number (0 if grid is not nested); from RAMSIN

  real, allocatable :: dtlongn(:)        ! delta t long

  integer, target :: nnxp(maxgrds)       ! global grid cells at x direction; from RAMSIN

  integer, allocatable :: nnx(:)         ! nnxp - 1; set by gridinit
  integer, allocatable :: nnx1(:)        ! nnxp - 2; set by gridinit
  integer, allocatable :: nnx2(:)        ! nnxp - 3; set by gridinit

  integer :: nstratx(maxgrds)            ! nest ratio for next coarser grid; from RAMSIN

  real, allocatable :: deltaxn(:)        ! delta x; set by gridset(1)

  integer, target :: nnyp(maxgrds)       ! grid cells at y direction; from RAMSIN

  integer, allocatable :: nny(:)         ! nnyp - 1; set by gridinit
  integer, allocatable :: nny1(:)        ! nnyp - 2; set by gridinit
  integer, allocatable :: nny2(:)        ! nnyp - 3; set by gridinit

  integer :: nstraty(maxgrds)            ! nest ratio for next coarser grid; from RAMSIN

  real, allocatable :: deltayn(:)        ! delta y; set by gridset(1)

  integer, target :: nnzp(maxgrds)       ! grid points z direction; from RAMSIN
  integer, allocatable :: nnz(:)         ! nnzp - 1; set by gridinit
  integer, allocatable :: nnz1(:)        ! nnzp - 2; set by gridinit
  real, allocatable    :: deltazn(:)     ! delta z; set by gridset(1)

  integer, allocatable :: nnxyp(:)       ! nnxp*nnyp (grid points at each vertical); set by gridinit
  integer, allocatable :: nnxyzp(:)      ! nnxp*nnyp*nnzp (grid points at the air); set by gridinit
  integer, allocatable :: nnxysp(:)      ! nnxp*nnyp*(nzg+nzs+3)*npatch (grid points beneath ground); set by gridinit

  real, allocatable :: platn(:)          ! pole latitude (degrees); set by gridset(1)
  real, allocatable :: plonn(:)          ! pole longitude (degrees); set by gridset(1)

  real :: centlat(maxgrds)               ! grid center latitude (degrees); from RAMSIN
  real :: centlon(maxgrds)               ! grid center longitude (degrees); from RAMSIN

  ! global grid cells coordinates and indices; indexed by global grid, not by a domain decomposed grid

  real, allocatable, target :: xtn(:,:)          ! x coordinate of cell center on polar stereographic projection; set by gridset
  real, allocatable :: xmn(:,:)          ! x coordinate of higher cell boundary on polar stereographic projection; set by gridset

  integer :: ninest(maxgrds)             ! index on next coarser grid where this grid starts (lower southwest corner); from RAMSIN or set by gridset(1)
  !                                      ! ds to interpolate x direction from coarser to finner grids

  integer, allocatable :: ipm(:,:)       ! next coarser grid cell index (icoarser) that contains this finer grid cell; set by gridset
  real, allocatable    :: ei1(:,:)          ! for icoarser-1 on 3 points interpolation; set by cofnest
  real, allocatable    :: ei2(:,:)          ! for icoarser   on 3 points interpolation; set by cofnest
  real, allocatable    :: ei3(:,:)          ! for icoarser+1 on 3 points interpolation; set by cofnest
  real, allocatable    :: ei4(:,:)          ! for icoarser-2 on 4 points interpolation; set by cofnest
  real, allocatable    :: ei5(:,:)          ! for icoarser-1 on 4 points interpolation; set by cofnest
  real, allocatable    :: ei6(:,:)          ! for icoarser   on 4 points interpolation; set by cofnest
  real, allocatable    :: ei7(:,:)          ! for icoarser+1 on 4 points interpolation; set by cofnest

  real, allocatable, target :: ytn(:,:)          ! y coordinate of cell center on polar stereographic projection; set by gridset
  real, allocatable :: ymn(:,:)          ! y coordinate of higher cell boundary on polar stereographic projection; set by gridset

  integer :: njnest(maxgrds)             ! index on next coarser grid where this grid starts (lower southwest corner)
  !                                      ! ds to interpolate y direction from coarser to finner grids; from RAMSIN

  integer, allocatable :: jpm(:,:)       ! next coarser grid cell index (jcoarser) that contains this finer grid cell; set by gridset
  real, allocatable    :: ej1(:,:)       ! for jcoarser-1 on 3 points interpolation; set by cofnest
  real, allocatable    :: ej2(:,:)       ! for jcoarser   on 3 points interpolation; set by cofnest
  real, allocatable    :: ej3(:,:)       ! for jcoarser+1 on 3 points interpolation; set by cofnest
  real, allocatable    :: ej4(:,:)       ! for jcoarser-2 on 4 points interpolation; set by cofnest
  real, allocatable    :: ej5(:,:)       ! for jcoarser-1 on 4 points interpolation; set by cofnest
  real, allocatable    :: ej6(:,:)       ! for jcoarser   on 4 points interpolation; set by cofnest
  real, allocatable    :: ej7(:,:)       ! for jcoarser+1 on 4 points interpolation; set by cofnest

  real, allocatable, target :: ztn(:,:)    ! z coordinate of interval center; set by gridset(1)
  real, allocatable, target :: zmn(:,:)            ! z coordinate of grid point; set by gridset(1)

  integer :: nknest(maxgrds)             ! index on next coarser grid where this grid starts (lower level)
  !                                      ! ds to interpolate z direction (kcoarser) from coarser to finner grids; from RAMSIN

  integer, allocatable :: kpm(:,:)       ! next coarser grid cell index that contains this finer grid cell; set by gridset
  real, allocatable    :: ek1(:,:)       ! for kcoarser-1 on 3 points interpolation; set by cofnest
  real, allocatable    :: ek2(:,:)       ! for kcoarser   on 3 points interpolation; set by cofnest
  real, allocatable    :: ek3(:,:)       ! for kcoarser+1 on 3 points interpolation; set by cofnest
  real, allocatable    :: ek4(:,:)       ! for kcoarser-2 on 4 points interpolation; set by cofnest
  real, allocatable    :: ek5(:,:)       ! for kcoarser-1 on 4 points interpolation; set by cofnest
  real, allocatable    :: ek6(:,:)       ! for kcoarser   on 4 points interpolation; set by cofnest
  real, allocatable    :: ek7(:,:)       ! for kcoarser+1 on 4 points interpolation; set by cofnest

  real, allocatable :: htn(:,:)
  real, allocatable :: ht2n(:,:)
  real, allocatable :: ht4n(:,:)
  real, allocatable :: hwn(:,:)
  real, allocatable :: hw2n(:,:)
  real, allocatable :: hw4n(:,:)
  real, allocatable, target :: dztn(:,:) !  set by gridset(1)
  real, allocatable, target :: dzmn(:,:) !  set by gridset(1)
  real, allocatable :: dzt2n(:,:)
  real, allocatable :: dzm2n(:,:)

  integer :: nxp
  integer :: nx
  integer :: nx1
  integer :: nx2
  integer :: nyp
  integer :: ny
  integer :: ny1
  integer :: ny2
  integer :: nzp
  integer :: nzpp
  integer :: nz
  integer :: nz1
  integer(kind=i8) :: nxyzp

  integer(kind=i8) :: nxyp

  integer(kind=i8) :: nxysp
  integer :: nscl
  integer :: nsttop
  integer :: nstbot
  integer :: ndtrat

  integer :: jdim                        ! all horizontal grids are 1D (jdim=0) or 2D (jdim=1); set by gridinit
  real :: deltax ! from RAMSIN
  real :: deltay ! from RAMSIN
  real :: deltaz ! from RAMSIN

  real, allocatable :: ht(:)
  real, allocatable :: ht2(:)
  real, allocatable :: ht4(:)
  real, allocatable :: hw(:)
  real, allocatable :: hw2(:)
  real, allocatable :: hw4(:)
  real, allocatable :: zt(:)
  real, allocatable :: zm(:)
  real, allocatable :: dzt(:)
  real, allocatable :: dzm(:)
  real, allocatable :: dzt2(:)
  real, allocatable :: dzm2(:)
  real, allocatable :: xt(:)
  real, allocatable :: xm(:)
  real, allocatable :: yt(:)
  real, allocatable :: ym(:)

  integer :: ngrid                       ! current grid;
  integer :: nzg                         ! soil layers; from RAMSIN
  integer :: nzs                         ! snow layers; from RAMSIN
  integer :: npatch                      ! surface patches per grid cell; from RAMSIN
  integer :: if_adap ! from RAMSIN
  integer :: itopo
  integer :: ihtran ! from RAMSIN
  integer :: ngridc
  integer :: ngrido
  integer :: iscr1
  integer :: iscr2
  integer :: memsize
  integer :: iounit
  integer :: maxpro
  integer :: memscr
  integer :: memind
  integer :: iogrid
  integer :: maxpts
  integer :: maxnzp
  integer :: maxnxp
  integer :: maxnyp
  integer :: i2dvar
  real :: time
  real :: ztop  ! set by gridset(1)
  real :: dzrat ! from RAMSIN
  real :: dzmax ! from RAMSIN
  integer :: fixLevels ! From RAMSIN
  real :: eps
  integer :: impl
  integer :: ideltat ! from RAMSIN
  integer :: iyear1 ! from RAMSIN
  integer :: imonth1 ! from RAMSIN
  integer :: idate1 ! from RAMSIN
  integer :: ihour1
  integer :: itime1 ! from RAMSIN
  integer :: nacoust ! from RAMSIN
  integer :: initial ! from RAMSIN
  integer :: iflag

  integer, allocatable :: nnacoust(:)

  real, allocatable :: dimove(:)
  real, allocatable :: djmove(:)

  real :: gridu(maxgrds) ! from RAMSIN
  real :: gridv(maxgrds) ! from RAMSIN
  real :: zz(nzpmax) ! from RAMSIN
  real :: dtlong ! from RAMSIN
  real :: sspct
  real :: polelat ! from RAMSIN
  real :: polelon ! from RAMSIN

  real, allocatable :: cflxy(:)
  real, allocatable :: cflz(:)
  !MB: real, allocatable :: cfl_max_sum(:)

  character(len=16) :: runtype ! from RAMSIN
  character(len=1)  :: timeunit ! from RAMSIN
  integer :: isstp
  integer :: istp
  real    :: timmax ! from RAMSIN
  real    :: dts
  real    :: dtlt
  real    :: dtlv
  integer :: nestz1 ! from RAMSIN
  integer :: nestz2 ! from RAMSIN
  integer :: nndtrat(maxgrds)            ! delta t ratio (coarser/nested), indexed by nested; from RAMSIN

  integer, allocatable :: ngbegun(:)

  integer :: nnsttop(maxgrds) ! from RAMSIN
  integer :: nnstbot(maxgrds) ! from RAMSIN
  integer :: nstratz1(nzpmax) ! from RAMSIN
  integer :: nstratz2(nzpmax) ! from RAMSIN


  integer, parameter :: maxsched=2000    ! maximum number of nested timesteps for a dtlong time advance
  integer, parameter :: maxschent=5      ! number of events to be recorded at each nested timestep
  integer :: nsubs                       ! actual number of nested timesteps for a dtlong time advance
  integer :: isched(maxsched,maxschent)  ! nested timestep events (see modsched for detailed description)


  integer :: nrzflg

  integer, allocatable :: nrz(:,:)         ! set by gridset(1)
  real, allocatable    :: fbcf(:,:,:)         ! set by cofnest

  integer :: iadvl
  integer :: iadvf
  integer :: lsflg ! from RAMSIN
  integer :: ibnd ! from RAMSIN
  integer :: jbnd ! from RAMSIN
  integer :: icorflg ! from RAMSIN
  integer :: dyncore_flag          ! =0 or 1: leapfrog, =2: Runge-Kutta dyn. core
  integer :: pd_or_mnt_constraint
  integer :: order_h
  integer :: order_v

  integer :: vveldamp ! from RAMSIN
  integer :: nfpt ! from RAMSIN
  integer :: naddsc ! from RAMSIN
  integer :: iversion
  real :: distim  ! from RAMSIN
  real :: cphas ! From RAMSIN

  !**(JP)** this should disapear!!!!
  ! Global simulation parameters

  integer :: nhemgrd2                    ! second hemispheric grid (0 if not global simulation); set by gridset(1)
  integer :: nhemt
  integer :: nhemu
  integer :: nhemv


  ! ALF
  ! Flags to set the thermo call on the horizontal boundaries

  logical, allocatable :: f_thermo_e(:)
  logical, allocatable :: f_thermo_w(:)
  logical, allocatable :: f_thermo_n(:)
  logical, allocatable :: f_thermo_s(:)

  type akmintype
     real, pointer :: akmin2d(:,:)
  end type akmintype

  type(akmintype), allocatable :: akminvar(:)

 !RMF
 !digital filter
 real :: begtime

contains

  ! *********************************************************************

  SUBROUTINE createMemGrid(ngrids, nnxp, nnyp, nnzp)
    IMPLICIT NONE
    ! Arguments:
    INTEGER, INTENT(IN) :: ngrids
    INTEGER, INTENT(IN) :: nnxp(maxgrds) ! From RAMSIN
    INTEGER, INTENT(IN) :: nnyp(maxgrds) ! From RAMSIN
    INTEGER, INTENT(IN) :: nnzp(maxgrds) ! From RAMSIN
    ! Local variables:
    INTEGER :: ierr, maxx, maxy, maxz !, maxxyz !ng

    ALLOCATE(dtlongn(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating dtlongn (createMemGrid)")
    ! Initiating dtlongn
    dtlongn = 0

    ALLOCATE(nnx(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating nnx (createMemGrid)")
    ALLOCATE(nnx1(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating nnx1 (createMemGrid)")
    ALLOCATE(nnx2(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating nnx2 (createMemGrid)")
    ALLOCATE(deltaxn(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating deltaxn (createMemGrid)")

    ALLOCATE(nny(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating nny (createMemGrid)")
    ALLOCATE(nny1(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating nny1 (createMemGrid)")
    ALLOCATE(nny2(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating nny2 (createMemGrid)")
    ALLOCATE(deltayn(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating deltayn (createMemGrid)")

    ALLOCATE(nnz(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating nnz (createMemGrid)")
    ALLOCATE(nnz1(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating nnz1 (createMemGrid)")
    ALLOCATE(deltazn(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating deltazn (createMemGrid)")

    ALLOCATE(nnxyp(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating nnxyp (createMemGrid)")
    ALLOCATE(nnxyzp(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating nnxyzp (createMemGrid)")
    ALLOCATE(nnxysp(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating nnxysp (createMemGrid)")

    ALLOCATE(platn(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating platn (createMemGrid)")
    ALLOCATE(plonn(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating plonn (createMemGrid)")

    maxx = maxval(nnxp(1:ngrids))
    ALLOCATE(xtn(maxx,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating xtn (createMemGrid)")
    ALLOCATE(xmn(maxx,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating xmn (createMemGrid)")
    ALLOCATE(ipm(maxx,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ipm (createMemGrid)")
    ALLOCATE(ei1(maxx,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ei1 (createMemGrid)")
    ALLOCATE(ei2(maxx,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ei2 (createMemGrid)")
    ALLOCATE(ei3(maxx,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ei3 (createMemGrid)")
    ALLOCATE(ei4(maxx,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ei4 (createMemGrid)")
    ALLOCATE(ei5(maxx,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ei5 (createMemGrid)")
    ALLOCATE(ei6(maxx,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ei6 (createMemGrid)")
    ALLOCATE(ei7(maxx,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ei7 (createMemGrid)")

    maxy = maxval(nnyp(1:ngrids))
    ALLOCATE(ytn(maxy,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ytn (createMemGrid)")
    ALLOCATE(ymn(maxy,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ymn (createMemGrid)")
    ALLOCATE(jpm(maxy,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating jpm (createMemGrid)")
    ALLOCATE(ej1(maxy,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ej1 (createMemGrid)")
    ALLOCATE(ej2(maxy,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ej2 (createMemGrid)")
    ALLOCATE(ej3(maxy,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ej3 (createMemGrid)")
    ALLOCATE(ej4(maxy,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ej4 (createMemGrid)")
    ALLOCATE(ej5(maxy,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ej5 (createMemGrid)")
    ALLOCATE(ej6(maxy,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ej6 (createMemGrid)")
    ALLOCATE(ej7(maxy,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ej7 (createMemGrid)")

    maxz = maxval(nnzp(1:ngrids))
    ALLOCATE(ztn(maxz,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ztn (createMemGrid)")
    ALLOCATE(zmn(maxz,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating zmn (createMemGrid)")
    ALLOCATE(kpm(maxz,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating kpm (createMemGrid)")
    ALLOCATE(ek1(maxz,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ek1 (createMemGrid)")
    ALLOCATE(ek2(maxz,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ek2 (createMemGrid)")
    ALLOCATE(ek3(maxz,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ek3 (createMemGrid)")
    ALLOCATE(ek4(maxz,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ek4 (createMemGrid)")
    ALLOCATE(ek5(maxz,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ek5 (createMemGrid)")
    ALLOCATE(ek6(maxz,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ek6 (createMemGrid)")
    ALLOCATE(ek7(maxz,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ek7 (createMemGrid)")

    ALLOCATE(htn(maxz,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating htn (createMemGrid)")
    ALLOCATE(ht2n(maxz,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ht2n (createMemGrid)")
    ALLOCATE(ht4n(maxz,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ht4n (createMemGrid)")
    ALLOCATE(hwn(maxz,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating hwn (createMemGrid)")
    ALLOCATE(hw2n(maxz,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating hw2n (createMemGrid)")
    ALLOCATE(hw4n(maxz,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating hw4n (createMemGrid)")
    ALLOCATE(dztn(maxz,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating dztn (createMemGrid)")
    ALLOCATE(dzmn(maxz,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating dzmn (createMemGrid)")
    ALLOCATE(dzt2n(maxz,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating dzt2n (createMemGrid)")
    ALLOCATE(dzm2n(maxz,ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating dzm2n (createMemGrid)")

    ALLOCATE(ht(maxz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ht (createMemGrid)")
    ALLOCATE(ht2(maxz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ht2 (createMemGrid)")
    ALLOCATE(ht4(maxz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ht4 (createMemGrid)")
    ALLOCATE(hw(maxz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating hw (createMemGrid)")
    ALLOCATE(hw2(maxz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating hw2 (createMemGrid)")
    ALLOCATE(hw4(maxz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating hw4 (createMemGrid)")
    ALLOCATE(zt(maxz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating zt (createMemGrid)")
    ALLOCATE(zm(maxz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating zm (createMemGrid)")
    ALLOCATE(dzt(maxz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating dzt (createMemGrid)")
    ALLOCATE(dzm(maxz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating dzm (createMemGrid)")
    ALLOCATE(dzt2(maxz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating dzt2 (createMemGrid)")
    ALLOCATE(dzm2(maxz), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating dzm2 (createMemGrid)")
    ALLOCATE(xt(maxx), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating xt (createMemGrid)")
    ALLOCATE(xm(maxx), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating xm (createMemGrid)")
    ALLOCATE(yt(maxy), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating yt (createMemGrid)")
    ALLOCATE(ym(maxy), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ym (createMemGrid)")

    ALLOCATE(nnacoust(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating nnacoust (createMemGrid)")

    ! Initializing nnacoust for binary reprodutibility pourpouses
    nnacoust = 0

    ALLOCATE(dimove(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating dimove (createMemGrid)")
    ALLOCATE(djmove(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating djmove (createMemGrid)")

    ! Initiating dimove and djmove ! Not using it
    ! Given initial values just for binary reprodutibility pourpouses
    dimove = 0
    djmove = 0

    ALLOCATE(cflxy(ngrids), STAT=ierr)

    !--(DMK-LFR NEC-SX6)----------------------------------------------
    cflxy = 0.
    !--(DMK-LFR NEC-SX6)----------------------------------------------

    IF (ierr/=0) CALL fatal_error("ERROR allocating cflxy (createMemGrid)")

    ALLOCATE(cflz(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating cflz (createMemGrid)")

    !MB: print*, "allocate cfl_max_sum"
    !ALLOCATE(cfl_max_sum(ngrids), STAT=ierr)
    !IF (ierr/=0) CALL fatal_error("ERROR allocating cfl_max_sum (createMemGrid)")

    !--(DMK-LFR NEC-SX6)----------------------------------------------
    cflz = 0.
    !MB: cfl_max_sum(:) = 0.0
    !--(DMK-LFR NEC-SX6)----------------------------------------------

    ALLOCATE(ngbegun(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating ngbegun (createMemGrid)")

    ALLOCATE(nrz(nzpmax,ngrids), STAT=ierr) !maxz
    IF (ierr/=0) CALL fatal_error("ERROR allocating nrz (createMemGrid)")
    ALLOCATE(fbcf(maxz,ngrids,4), STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR allocating fbcf (createMemGrid)")

    ALLOCATE(f_thermo_e(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error (&
         "ERROR allocating f_thermo_e (createMemGrid)")
    ALLOCATE(f_thermo_w(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error (&
         "ERROR allocating f_thermo_w (createMemGrid)")
    ALLOCATE(f_thermo_n(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error (&
         "ERROR allocating f_thermo_n (createMemGrid)")
    ALLOCATE(f_thermo_s(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error (&
         "ERROR allocating f_thermo_s (createMemGrid)")

    ALLOCATE(akminvar(ngrids), STAT=ierr)
    IF (ierr/=0) CALL fatal_error (&
         "ERROR allocating akminvar (createMemGrid)")

  END SUBROUTINE createMemGrid





  SUBROUTINE allocAkmin2d(ngrids, nodemxp, nodemyp)
    IMPLICIT NONE
    ! Arguments:
    INTEGER, INTENT(IN) :: ngrids
    INTEGER, INTENT(IN), TARGET ::  nodemxp(ngrids), nodemyp(ngrids)
    ! Local variables:
    INTEGER :: ifm, ierr

    DO ifm=1,ngrids
       ALLOCATE(akminvar(ifm)%akmin2d(nodemxp(ifm), nodemyp(ifm)), STAT=ierr)
       IF (ierr/=0) CALL fatal_error (&
            "ERROR allocating akmin2d (allocAkmin2d)")
    ENDDO

  END SUBROUTINE allocAkmin2d

  ! *********************************************************************

  SUBROUTINE destroyMemGrid()
    IMPLICIT NONE
    ! Local variables:
    INTEGER :: ierr

    DEALLOCATE(dtlongn, STAT=ierr)
    IF (ierr/=0) CALL fatal_error (&
         "ERROR deallocating dtlongn (destroyMemGrid)")

    DEALLOCATE(nnx, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating nnx (destroyMemGrid)")
    DEALLOCATE(nnx1, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating nnx1 (destroyMemGrid)")
    DEALLOCATE(nnx2, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating nnx2 (destroyMemGrid)")
    DEALLOCATE(deltaxn, STAT=ierr)
    IF (ierr/=0) CALL fatal_error (&
         "ERROR deallocating deltaxn (destroyMemGrid)")

    DEALLOCATE(nny, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating nny (destroyMemGrid)")
    DEALLOCATE(nny1, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating nny1 (destroyMemGrid)")
    DEALLOCATE(nny2, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating nny2 (destroyMemGrid)")
    DEALLOCATE(deltayn, STAT=ierr)
    IF (ierr/=0) CALL fatal_error (&
         "ERROR deallocating deltayn (destroyMemGrid)")

    DEALLOCATE(nnz, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating nnz (destroyMemGrid)")
    DEALLOCATE(nnz1, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating nnz(destroyMemGrid)")
    DEALLOCATE(deltazn, STAT=ierr)
    IF (ierr/=0) CALL fatal_error (&
         "ERROR deallocating deltazn (destroyMemGrid)")

    DEALLOCATE(nnxyp, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating nnxyp (destroyMemGrid)")
    DEALLOCATE(nnxyzp, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating nnxyzp (destroyMemGrid)")
    DEALLOCATE(nnxysp, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating nnxysp (destroyMemGrid)")

    DEALLOCATE(platn, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating platn (destroyMemGrid)")
    DEALLOCATE(plonn, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating plonn (destroyMemGrid)")

    DEALLOCATE(xtn, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating xtn (destroyMemGrid)")
    DEALLOCATE(xmn, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating xmn (destroyMemGrid)")
    DEALLOCATE(ipm, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ipm (destroyMemGrid)")
    DEALLOCATE(ei1, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ei1 (destroyMemGrid)")
    DEALLOCATE(ei2, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ei2 (destroyMemGrid)")
    DEALLOCATE(ei3, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ei3 (destroyMemGrid)")
    DEALLOCATE(ei4, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ei4 (destroyMemGrid)")
    DEALLOCATE(ei5, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ei5 (destroyMemGrid)")
    DEALLOCATE(ei6, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ei6 (destroyMemGrid)")
    DEALLOCATE(ei7, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ei7 (destroyMemGrid)")

    DEALLOCATE(ytn, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ytn (destroyMemGrid)")
    DEALLOCATE(ymn, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ymn (destroyMemGrid)")
    DEALLOCATE(jpm, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating jpm (destroyMemGrid)")
    DEALLOCATE(ej1, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ej1 (destroyMemGrid)")
    DEALLOCATE(ej2, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ej2 (destroyMemGrid)")
    DEALLOCATE(ej3, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ej3 (destroyMemGrid)")
    DEALLOCATE(ej4, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ej4 (destroyMemGrid)")
    DEALLOCATE(ej5, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ej5 (destroyMemGrid)")
    DEALLOCATE(ej6, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ej6 (destroyMemGrid)")
    DEALLOCATE(ej7, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ej7 (destroyMemGrid)")

    DEALLOCATE(ztn, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ztn (destroyMemGrid)")
    DEALLOCATE(zmn, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating zmn (destroyMemGrid)")
    DEALLOCATE(kpm, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating kpm (destroyMemGrid)")
    DEALLOCATE(ek1, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ek1 (destroyMemGrid)")
    DEALLOCATE(ek2, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ek2 (destroyMemGrid)")
    DEALLOCATE(ek3, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ek3 (destroyMemGrid)")
    DEALLOCATE(ek4, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ek4 (destroyMemGrid)")
    DEALLOCATE(ek5, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ek5 (destroyMemGrid)")
    DEALLOCATE(ek6, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ek6 (destroyMemGrid)")
    DEALLOCATE(ek7, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ek7 (destroyMemGrid)")

    DEALLOCATE(htn, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating htn (destroyMemGrid)")
    DEALLOCATE(ht2n, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ht2n (destroyMemGrid)")
    DEALLOCATE(ht4n, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ht4n (destroyMemGrid)")
    DEALLOCATE(hwn, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating hwn (destroyMemGrid)")
    DEALLOCATE(hw2n, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating hw2n (destroyMemGrid)")
    DEALLOCATE(hw4n, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating hw4n (destroyMemGrid)")
    DEALLOCATE(dztn, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating dztn (destroyMemGrid)")
    DEALLOCATE(dzmn, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating dzmn (destroyMemGrid)")
    DEALLOCATE(dzt2n, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating dzt2n (destroyMemGrid)")
    DEALLOCATE(dzm2n, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating dzm2n (destroyMemGrid)")

    DEALLOCATE(ht, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ht (destroyMemGrid)")
    DEALLOCATE(ht2, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ht2 (destroyMemGrid)")
    DEALLOCATE(ht4, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ht4 (destroyMemGrid)")
    DEALLOCATE(hw, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating hw (destroyMemGrid)")
    DEALLOCATE(hw2, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating hw2 (destroyMemGrid)")
    DEALLOCATE(hw4, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating hw4 (destroyMemGrid)")
    DEALLOCATE(zt, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating zt (destroyMemGrid)")
    DEALLOCATE(zm, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating zm (destroyMemGrid)")
    DEALLOCATE(dzt, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating dzt (destroyMemGrid)")
    DEALLOCATE(dzm, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating dzm (destroyMemGrid)")
    DEALLOCATE(dzt2, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating dzt2 (destroyMemGrid)")
    DEALLOCATE(dzm2, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating dzm2 (destroyMemGrid)")
    DEALLOCATE(xt, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating xt (destroyMemGrid)")
    DEALLOCATE(xm, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating xm (destroyMemGrid)")
    DEALLOCATE(yt, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating yt (destroyMemGrid)")
    DEALLOCATE(ym, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating ym (destroyMemGrid)")

    DEALLOCATE(nnacoust, STAT=ierr)
    IF (ierr/=0) CALL fatal_error (&
         "ERROR deallocating nnacoust (destroyMemGrid)")
    DEALLOCATE(dimove, STAT=ierr)
    IF (ierr/=0) CALL fatal_error (&
         "ERROR deallocating dimove (destroyMemGrid)")
    DEALLOCATE(djmove, STAT=ierr)
    IF (ierr/=0) CALL fatal_error (&
         "ERROR deallocating djmove (destroyMemGrid)")

    DEALLOCATE(cflxy, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating cflxy (destroyMemGrid)")

    DEALLOCATE(cflz, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating cflz (destroyMemGrid)")

    !MB: print*, "deallocate cfl_max_sum"
    !DEALLOCATE(cfl_max_sum, STAT=ierr)
    !IF (ierr/=0) CALL fatal_error("ERROR deallocating cfl_max_sum (destroyMemGrid)")

    DEALLOCATE(ngbegun, STAT=ierr)
    IF (ierr/=0) CALL fatal_error (&
         "ERROR deallocating ngbegun (destroyMemGrid)")

    DEALLOCATE(nrz, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating nrz (destroyMemGrid)")
    DEALLOCATE(fbcf, STAT=ierr)
    IF (ierr/=0) CALL fatal_error("ERROR deallocating fbcf (destroyMemGrid)")

    DEALLOCATE(f_thermo_e, STAT=ierr)
    IF (ierr/=0) CALL fatal_error (&
         "ERROR deallocating f_thermo_e (destroyMemGrid)")
    DEALLOCATE(f_thermo_w, STAT=ierr)
    IF (ierr/=0) CALL fatal_error (&
         "ERROR deallocating f_thermo_w (destroyMemGrid)")
    DEALLOCATE(f_thermo_n, STAT=ierr)
    IF (ierr/=0) CALL fatal_error(&
         "ERROR deallocating f_thermo_n (destroyMemGrid)")
    DEALLOCATE(f_thermo_s, STAT=ierr)
    IF (ierr/=0) CALL fatal_error(&
         "ERROR deallocating f_thermo_s (destroyMemGrid)")

  END SUBROUTINE destroyMemGrid


  ! *********************************************************************

  subroutine alloc_grid(grid,n1,n2,n3,ng,if_adap)
    implicit none
    type (grid_vars) :: grid
    integer, intent(in) :: n1,n2,n3,ng,if_adap

    ! Allocate arrays based on options (if necessary)

    allocate (grid%topt(n2,n3))
    allocate (grid%topu(n2,n3))
    allocate (grid%topv(n2,n3))
    allocate (grid%topm(n2,n3))
    allocate (grid%topma(n2,n3))
    allocate (grid%topta(n2,n3))
    allocate (grid%rtgt(n2,n3))
    allocate (grid%rtgu(n2,n3))
    allocate (grid%rtgv(n2,n3))
    allocate (grid%rtgm(n2,n3))
    allocate (grid%f13t(n2,n3))
    allocate (grid%f13u(n2,n3))
    allocate (grid%f13v(n2,n3))
    allocate (grid%f13m(n2,n3))
    allocate (grid%f23t(n2,n3))
    allocate (grid%f23u(n2,n3))
    allocate (grid%f23v(n2,n3))
    allocate (grid%f23m(n2,n3))
    allocate (grid%dxt(n2,n3))
    allocate (grid%dxu(n2,n3))
    allocate (grid%dxv(n2,n3))
    allocate (grid%dxm(n2,n3))
    allocate (grid%dyt(n2,n3))
    allocate (grid%dyu(n2,n3))
    allocate (grid%dyv(n2,n3))
    allocate (grid%dym(n2,n3))
    allocate (grid%fmapt(n2,n3))
    allocate (grid%fmapu(n2,n3))
    allocate (grid%fmapv(n2,n3))
    allocate (grid%fmapm(n2,n3))
    allocate (grid%fmapti(n2,n3))
    allocate (grid%fmapui(n2,n3))
    allocate (grid%fmapvi(n2,n3))
    allocate (grid%fmapmi(n2,n3))
    allocate (grid%glat(n2,n3))
    allocate (grid%glon(n2,n3))
    allocate (grid%topzo(n2,n3))
!TO modify to Xeon-Phi Intel
    if (if_adap == 1) then
       allocate (grid%aru(n1,n2,n3))
       allocate (grid%arv(n1,n2,n3))
       allocate (grid%arw(n1,n2,n3))
       allocate (grid%volu(n1,n2,n3))
       allocate (grid%volv(n1,n2,n3))
       allocate (grid%volw(n1,n2,n3))
       allocate (grid%volt(n1,n2,n3))
    endif
    allocate (grid%lpu(n2,n3))
    allocate (grid%lpv(n2,n3))
    allocate (grid%lpw(n2,n3))

    !--(DMK-LFR NEC-SX6)----------------------------------------------
    grid%topt = 0.
    grid%topu = 0.
    grid%topv = 0.
    grid%topm = 0.
    grid%topma = 0.
    grid%topta = 0.
    grid%rtgt = 0.
    grid%rtgu = 0.
    grid%rtgv = 0.
    grid%rtgm = 0.
    grid%f13t = 0.
    grid%f13u = 0.
    grid%f13v = 0.
    grid%f13m = 0.
    grid%f23t = 0.
    grid%f23u = 0.
    grid%f23v = 0.
    grid%f23m = 0.
    grid%dxt = 0.
    grid%dxu = 0.
    grid%dxv = 0.
    grid%dxm = 0.
    grid%dyt = 0.
    grid%dyu = 0.
    grid%dyv = 0.
    grid%dym = 0.
    grid%fmapt = 0.
    grid%fmapu = 0.
    grid%fmapv = 0.
    grid%fmapm = 0.
    grid%fmapti = 0.
    grid%fmapui = 0.
    grid%fmapvi = 0.
    grid%fmapmi = 0.
    grid%glat = 0.
    grid%glon = 0.
    grid%topzo = 0.
!TO modify to Xeon-Phi Intel
    if (if_adap == 1) then
       grid%aru = 0.
       grid%arv = 0.
       grid%arw = 0.
       grid%volu = 0.
       grid%volv = 0.
       grid%volw = 0.
       grid%volt = 0.
    endif
    grid%lpu = 0
    grid%lpv = 0
    grid%lpw =0
    !--(DMK-LFR NEC-SX6)----------------------------------------------

  end subroutine alloc_grid





  subroutine nullify_grid(grid)
    implicit none
    type (grid_vars) :: grid
    nullify (grid%topt)  ;  nullify (grid%topu)  ;  nullify (grid%topv)
    nullify (grid%topm)  ;  nullify (grid%topma) ;  nullify (grid%topta)
    nullify (grid%rtgt)  ;  nullify (grid%rtgu)
    nullify (grid%rtgv)  ;  nullify (grid%rtgm)  ;  nullify (grid%f13t)
    nullify (grid%f13u)  ;  nullify (grid%f13v)  ;  nullify (grid%f13m)
    nullify (grid%f23t)  ;  nullify (grid%f23u)  ;  nullify (grid%f23v)
    nullify (grid%f23m)  ;  nullify (grid%dxt)   ;  nullify (grid%dxu)
    nullify (grid%dxv)   ;  nullify (grid%dxm)   ;  nullify (grid%dyt)
    nullify (grid%dyu)   ;  nullify (grid%dyv)   ;  nullify (grid%dym)
    nullify (grid%fmapt) ;  nullify (grid%fmapu) ;  nullify (grid%fmapv)
    nullify (grid%fmapm) ;  nullify (grid%fmapti);  nullify (grid%fmapui)
    nullify (grid%fmapvi);  nullify (grid%fmapmi);  nullify (grid%glat)
    nullify (grid%glon)  ;  nullify (grid%topzo)
    nullify (grid%aru)   ;  nullify (grid%arv)   ;  nullify (grid%arw)
    nullify (grid%volu)  ;  nullify (grid%volv)  ;  nullify (grid%volw)
    nullify (grid%volt)
    nullify (grid%lpu)   ;  nullify (grid%lpv)   ;  nullify (grid%lpw)
  end subroutine nullify_grid




  subroutine dealloc_grid(grid)
    implicit none
    type (grid_vars) :: grid
    if (associated(grid%topt)  )    deallocate (grid%topt)
    if (associated(grid%topu)  )    deallocate (grid%topu)
    if (associated(grid%topv)  )    deallocate (grid%topv)
    if (associated(grid%topm)  )    deallocate (grid%topm)
    if (associated(grid%topma) )    deallocate (grid%topma)
    if (associated(grid%topta) )    deallocate (grid%topta)
    if (associated(grid%rtgt)  )    deallocate (grid%rtgt)
    if (associated(grid%rtgu)  )    deallocate (grid%rtgu)
    if (associated(grid%rtgv)  )    deallocate (grid%rtgv)
    if (associated(grid%rtgm)  )    deallocate (grid%rtgm)
    if (associated(grid%f13t)  )    deallocate (grid%f13t)
    if (associated(grid%f13u)  )    deallocate (grid%f13u)
    if (associated(grid%f13v)  )    deallocate (grid%f13v)
    if (associated(grid%f13m)  )    deallocate (grid%f13m)
    if (associated(grid%f23t)  )    deallocate (grid%f23t)
    if (associated(grid%f23u)  )    deallocate (grid%f23u)
    if (associated(grid%f23v)  )    deallocate (grid%f23v)
    if (associated(grid%f23m)  )    deallocate (grid%f23m)
    if (associated(grid%dxt)   )    deallocate (grid%dxt)
    if (associated(grid%dxu)   )    deallocate (grid%dxu)
    if (associated(grid%dxv)   )    deallocate (grid%dxv)
    if (associated(grid%dxm)   )    deallocate (grid%dxm)
    if (associated(grid%dyt)   )    deallocate (grid%dyt)
    if (associated(grid%dyu)   )    deallocate (grid%dyu)
    if (associated(grid%dyv)   )    deallocate (grid%dyv)
    if (associated(grid%dym)   )    deallocate (grid%dym)
    if (associated(grid%fmapt) )    deallocate (grid%fmapt)
    if (associated(grid%fmapu) )    deallocate (grid%fmapu)
    if (associated(grid%fmapv) )    deallocate (grid%fmapv)
    if (associated(grid%fmapm) )    deallocate (grid%fmapm)
    if (associated(grid%fmapti))    deallocate (grid%fmapti)
    if (associated(grid%fmapui))    deallocate (grid%fmapui)
    if (associated(grid%fmapvi))    deallocate (grid%fmapvi)
    if (associated(grid%fmapmi))    deallocate (grid%fmapmi)
    if (associated(grid%glat)  )    deallocate (grid%glat)
    if (associated(grid%glon)  )    deallocate (grid%glon)
    if (associated(grid%topzo) )    deallocate (grid%topzo)
    if (associated(grid%aru)   )    deallocate (grid%aru)
    if (associated(grid%arv)   )    deallocate (grid%arv)
    if (associated(grid%arw)   )    deallocate (grid%arw)
    if (associated(grid%volu)  )    deallocate (grid%volu)
    if (associated(grid%volv)  )    deallocate (grid%volv)
    if (associated(grid%volw)  )    deallocate (grid%volw)
    if (associated(grid%volt)  )    deallocate (grid%volt)
    if (associated(grid%lpu)   )    deallocate (grid%lpu)
    if (associated(grid%lpv)   )    deallocate (grid%lpv)
    if (associated(grid%lpw)   )    deallocate (grid%lpw)
  end subroutine dealloc_grid




  subroutine filltab_grid(grid,gridm,imean,n1,n2,n3,ng)

    use var_tables

    implicit none
    type (grid_vars) :: grid,gridm
    integer, intent(in) :: imean,n1,n2,n3,ng
    integer(kind=i8) :: npts
    real, pointer :: var,varm

    ! Fill pointers to arrays into variable tables

    npts=n2*n3
    if (associated(grid%topt)) &
         call InsertVTab (grid%topt,gridm%topt,ng,npts,imean,  &
         'TOPT :2:hist:anal:mpti')
    if (associated(grid%topu)) &
         call InsertVTab (grid%topu,gridm%topu,ng, npts, imean,  &
         'TOPU :2:mpti')
    if (associated(grid%topv)) &
         call InsertVTab (grid%topv,gridm%topv,ng, npts, imean,  &
         'TOPV :2:mpti')
    if (associated(grid%topm)) &
         call InsertVTab (grid%topm,gridm%topm,ng, npts, imean,  &
         'TOPM :2:mpti')
    if (associated(grid%topma)) &
         call InsertVTab (grid%topma,gridm%topma,ng, npts, imean,  &
         'TOPMA :2:hist:anal:mpti')
    if (associated(grid%topta)) &
         call InsertVTab (grid%topta,gridm%topta,ng, npts, imean,  &
         'TOPTA :2:hist:anal:mpti')
    if (associated(grid%rtgt)) &
         call InsertVTab (grid%rtgt,gridm%rtgt,ng, npts, imean,  &
         'RTGT :2:mpti')
    if (associated(grid%rtgu)) &
         call InsertVTab (grid%rtgu,gridm%rtgu,ng, npts, imean,  &
         'RTGU :2:mpti')
    if (associated(grid%rtgv)) &
         call InsertVTab (grid%rtgv,gridm%rtgv,ng, npts, imean,  &
         'RTGV :2:mpti')
    if (associated(grid%rtgm)) &
         call InsertVTab (grid%rtgm,gridm%rtgm,ng, npts, imean,  &
         'RTGM :2:mpti')
    if (associated(grid%f13t)) &
         call InsertVTab (grid%f13t,gridm%f13t,ng, npts, imean,  &
         'F13T :2:mpti')
    if (associated(grid%f13u)) &
         call InsertVTab (grid%f13u,gridm%f13u,ng, npts, imean,  &
         'F13U :2:mpti')
    if (associated(grid%f13v)) &
         call InsertVTab (grid%f13v,gridm%f13v,ng, npts, imean,  &
         'F13V :2:mpti')
    if (associated(grid%f13m)) &
         call InsertVTab (grid%f13m,gridm%f13m,ng, npts, imean,  &
         'F13M :2:mpti')
    if (associated(grid%f23t)) &
         call InsertVTab (grid%f23t,gridm%f23t,ng, npts, imean,  &
         'F23T :2:mpti')
    if (associated(grid%f23u)) &
         call InsertVTab (grid%f23u,gridm%f23u,ng, npts, imean,  &
         'F23U :2:mpti')
    if (associated(grid%f23v)) &
         call InsertVTab (grid%f23v,gridm%f23v,ng, npts, imean,  &
         'F23V :2:mpti')
    if (associated(grid%f23m)) &
         call InsertVTab (grid%f23m,gridm%f23m,ng, npts, imean,  &
         'F23M :2:mpti')
    if (associated(grid%dxt)) &
         call InsertVTab (grid%dxt,gridm%dxt,ng, npts, imean,  &
         'DXT :2:mpti')
    if (associated(grid%dxu)) &
         call InsertVTab (grid%dxu,gridm%dxu,ng, npts, imean,  &
         'DXU :2:mpti')
    if (associated(grid%dxv)) &
         call InsertVTab (grid%dxv,gridm%dxv,ng, npts, imean,  &
         'DXV :2:mpti')
    if (associated(grid%dxm)) &
         call InsertVTab (grid%dxm,gridm%dxm,ng, npts, imean,  &
         'DXM :2:mpti')
    if (associated(grid%dyt)) &
         call InsertVTab (grid%dyt,gridm%dyt,ng, npts, imean,  &
         'DYT :2:mpti')
    if (associated(grid%dyu)) &
         call InsertVTab (grid%dyu,gridm%dyu,ng, npts, imean,  &
         'DYU :2:mpti')
    if (associated(grid%dyv)) &
         call InsertVTab (grid%dyv,gridm%dyv,ng, npts, imean,  &
         'DYV :2:mpti')
    if (associated(grid%dym)) &
         call InsertVTab (grid%dym,gridm%dym,ng, npts, imean,  &
         'DYM :2:mpti')
    if (associated(grid%fmapt)) &
         call InsertVTab (grid%fmapt,gridm%fmapt,ng, npts, imean,  &
         'FMAPT :2:mpti')
    if (associated(grid%fmapu)) &
         call InsertVTab (grid%fmapu,gridm%fmapu,ng, npts, imean,  &
         'FMAPU :2:mpti')
    if (associated(grid%fmapv)) &
         call InsertVTab (grid%fmapv,gridm%fmapv,ng, npts, imean,  &
         'FMAPV :2:mpti')
    if (associated(grid%fmapm)) &
         call InsertVTab (grid%fmapm,gridm%fmapm,ng, npts, imean,  &
         'FMAPM :2:mpti')
    if (associated(grid%fmapti)) &
         call InsertVTab (grid%fmapti,gridm%fmapti,ng, npts, imean,  &
         'FMAPTI :2:mpti')
    if (associated(grid%fmapui)) &
         call InsertVTab (grid%fmapui,gridm%fmapui,ng, npts, imean,  &
         'FMAPUI :2:mpti')
    if (associated(grid%fmapvi)) &
         call InsertVTab (grid%fmapvi,gridm%fmapvi,ng, npts, imean,  &
         'FMAPVI :2:mpti')
    if (associated(grid%fmapmi)) &
         call InsertVTab (grid%fmapmi,gridm%fmapmi,ng, npts, imean,  &
         'FMAPMI :2:mpti')
    if (associated(grid%glat)) &
         call InsertVTab (grid%glat,gridm%glat,ng, npts, imean,  &
         'GLAT :2:mpti:anal')
    if (associated(grid%glon)) &
         call InsertVTab (grid%glon,gridm%glon,ng, npts, imean,  &
         'GLON :2:mpti:anal')
    if (associated(grid%topzo)) &
         call InsertVTab (grid%topzo,gridm%topzo,ng, npts, imean,  &
         'TOPZO :2:mpti')

    npts=n2*n3
    if (associated(grid%lpu)) &
         call InsertVTab (grid%lpu,gridm%lpu,ng,npts,imean,  &
         'LPU :2:mpti')
    if (associated(grid%lpv)) &
         call InsertVTab (grid%lpv,gridm%lpv,ng,npts,imean,  &
         'LPV :2:mpti')
    if (associated(grid%lpw)) &
         call InsertVTab (grid%lpw,gridm%lpw,ng,npts,imean,  &
         'LPW :2:mpti')

    npts=n1*n2*n3
    if (associated(grid%aru)) &
         call InsertVTab (grid%aru,gridm%aru,ng,npts,imean,  &
         'ARU :3:mpti')
    if (associated(grid%arv)) &
         call InsertVTab (grid%arv,gridm%arv,ng,npts,imean,  &
         'ARV :3:mpti')
    if (associated(grid%arw)) &
         call InsertVTab (grid%arw,gridm%arw,ng,npts,imean,  &
         'ARW :3:mpti')

    if (associated(grid%volu)) &
         call InsertVTab (grid%volu,gridm%volu,ng,npts,imean,  &
         'VOLU :3:mpti')
    if (associated(grid%volv)) &
         call InsertVTab (grid%volv,gridm%volv,ng,npts,imean,  &
         'VOLV :3:mpti')
    if (associated(grid%volw)) &
         call InsertVTab (grid%volw,gridm%volw,ng,npts,imean,  &
         'VOLW :3:mpti')
    if (associated(grid%volt)) &
         call InsertVTab (grid%volt,gridm%volt,ng,npts,imean,  &
         'VOLT :3:anal:mpti')

  end subroutine filltab_grid




!!$  subroutine filltabAkmin2d(ngrids, nodemxp, nodemyp)
!!$
!!$    implicit none
!!$    include "i8.h"
!!$    integer, intent(IN) :: ngrids, nodemxp(ngrids), nodemyp(ngrids)
!!$    integer(kind=i8) :: npts
!!$    integer          :: ifm, imean
!!$
!!$    imean = 0 ! does not change with time
!!$
!!$    do ifm=1,ngrids
!!$       npts = nodemxp(ifm)*nodemyp(ifm)
!!$       if (associated(akminvar(ifm)%akmin2d)) then
!!$          call vtables2 (akminvar(ifm)%akmin2d(1,1), &
!!$               akminvar(ifm)%akmin2d(1,1), ifm, &
!!$               npts, imean, 'AKMIN2D :2:anal')
!!$       endif
!!$    enddo
!!$
!!$  end subroutine filltabAkmin2d




  subroutine alloc_GlobalGridData(global, n2, n3)
    implicit none
    type(GlobalGridData), intent(inout) :: global
    integer, intent(in) :: n2
    integer, intent(in) :: n3

    allocate(global%global_topta(n2,n3))
    allocate(global%global_glat(n2,n3))
    allocate(global%global_glon(n2,n3))
  end subroutine alloc_GlobalGridData






  subroutine nullify_GlobalGridData(global)
    implicit none
    type(GlobalGridData), intent(inout) :: global

    nullify(global%global_topta)
    nullify(global%global_glat)
    nullify(global%global_glon)
  end subroutine nullify_GlobalGridData






  subroutine dealloc_GlobalGridData(global)
    implicit none
    type(GlobalGridData), intent(inout) :: global

    if (associated(global%global_topta)) deallocate(global%global_topta)
    if (associated(global%global_glat)) deallocate(global%global_glat)
    if (associated(global%global_glon)) deallocate(global%global_glon)
  end subroutine dealloc_GlobalGridData






  subroutine dump_mem_grid()

    ! dump_mem_grid: dumps mem_grid sizes and mappings for all grids

    implicit none
    integer :: ind, indGrid, j
    character(len=10) :: c0, c1
    character(len=*), parameter :: h="**(dump_mem_grid)**"

    do indGrid = 1, ngrids
       if (nxtnest(indGrid) == 0) then
          write(c0,"(i10)") indGrid
          write(*,"(a)") h//" Dumping outer grid number "//trim(adjustl(c0))
          write(c0,"(i10)") nnxp(indGrid)-2
          write(c1,"(f10.1)") deltaxn(indGrid)
          write(*,"(a)") h//" x axis has "//trim(adjustl(c0))//&
               " inner cells and 2 boundary cells of length "//trim(adjustl(c1))
          write(*,"(a)") h//" x cells: (index, higher cell boundary coordinate)"
          do ind = 1, nnxp(indGrid), 5
             do j = ind, min(ind+4,nnxp(indGrid))
                write(*,'(" (",i3,",",f10.1,");")',ADVANCE="NO") &
                     j,xmn(j,indGrid)
             end do
             write(*,"(1x)")
          end do
          write(c0,"(i10)") nnyp(indGrid)-2
          write(c1,"(f10.1)") deltayn(indGrid)
          write(*,"(a)") h//" y axis has "//trim(adjustl(c0))//&
               " inner cells and 2 boundary cells of length "//trim(adjustl(c1))
          write(*,"(a)") h//" y cells: (index, higher cell boundary coordinate)"
          do ind = 1, nnyp(indGrid), 5
             do j = ind, min(ind+4,nnyp(indGrid))
                write(*,'(" (",i3,",",f10.1,");")',ADVANCE="NO") &
                     j,ymn(j,indGrid)
             end do
             write(*,"(1x)")
          end do
       else
          write(*,"(a)") h
          write(c0,"(i10)") indGrid
          write(c1,"(i10)") nxtnest(indGrid)
          write(*,"(a)") h//" Dumping inner grid number "//trim(adjustl(c0))//&
               ", nested into grid number "//trim(adjustl(c1))
          write(c0,"(i10)") nnxp(indGrid)-2
          write(c1,"(f10.1)") deltaxn(indGrid)
          write(*,"(a)") h//" x axis has "//trim(adjustl(c0))//&
               " inner cells and 2 boundary cells of length "//trim(adjustl(c1))
          write(*,"(a)") h//" x cells: (index, higher cell boundary coordinate, coarser grid cell index)"
          do ind = 1, nnxp(indGrid), 5
             do j = ind, min(ind+4,nnxp(indGrid))
                write(*,'(" (",i3,",",f10.1,",",i3,");")',ADVANCE="NO") &
                     j,xmn(j,indGrid),ipm(j,indGrid)
             end do
             write(*,"(1x)")
          end do
          write(c0,"(i10)") nnyp(indGrid)-2
          write(c1,"(f10.1)") deltayn(indGrid)
          write(*,"(a)") h//" y axis has "//trim(adjustl(c0))//&
               " inner cells and 2 boundary cells of length "//trim(adjustl(c1))
          write(*,"(a)") h//" y cells: (index, higher cell boundary coordinate, coarser grid cell index)"
          do ind = 1, nnyp(indGrid), 5
             do j = ind, min(ind+4,nnyp(indGrid))
                write(*,'(" (",i3,",",f10.1,",",i3,");")',ADVANCE="NO") &
                     j,ymn(j,indGrid),jpm(j,indGrid)
             end do
             write(*,"(1x)")
          end do
       end if
    end do
  end subroutine dump_mem_grid



  ! GlobalSizes: size of full fields as a function of idim_type



  subroutine GlobalSizes(grid, nmachs_u, nwave, globalSize)
    implicit none
    integer, intent(in ) :: grid
    integer, intent(in ) :: nmachs_u
    integer, intent(in ) :: nwave
    integer, intent(out) :: globalSize(2:7)


    globalSize(2:7) =  (/ &
         nnxp(grid)*nnyp(grid),            nnzp(grid)*nnxp(grid)*nnyp(grid), &
         nzg*nnxp(grid)*nnyp(grid)*npatch, nzs*nnxp(grid)*nnyp(grid)*npatch, &
         nnxp(grid)*nnyp(grid)*npatch,     nnxp(grid)*nnyp(grid)*nwave         /)
  end subroutine GlobalSizes




  ! ExtractLocalFromGlobal: extract local domain of grid_vars from global domain




  subroutine ExtractLocalFromGlobal (global, nzp, nxp, nyp, i0, j0, if_adap, local)
    implicit none
    type(grid_vars), intent(in   ) :: global
    integer,         intent(in   ) :: nzp
    integer,         intent(in   ) :: nxp
    integer,         intent(in   ) :: nyp
    integer,         intent(in   ) :: i0
    integer,         intent(in   ) :: j0
    integer,         intent(in   ) :: if_adap
    type(grid_vars), intent(inout) :: local

    integer :: i, ig
    integer :: j, jg
    integer :: k

    do j= 1, nyp
       jg = j+j0
       do i = 1, nxp
          ig = i+i0
          local%topt(i,j)=global%topt(ig,jg)
          local%topu(i,j)=global%topu(ig,jg)
          local%topv(i,j)=global%topv(ig,jg)
          local%topm(i,j)=global%topm(ig,jg)
          local%topma(i,j)=global%topma(ig,jg)
          local%topta(i,j)=global%topta(ig,jg)
          local%rtgt(i,j)=global%rtgt(ig,jg)
          local%rtgu(i,j)=global%rtgu(ig,jg)
          local%rtgv(i,j)=global%rtgv(ig,jg)
          local%rtgm(i,j)=global%rtgm(ig,jg)
          local%f13t(i,j)=global%f13t(ig,jg)
          local%f13u(i,j)=global%f13u(ig,jg)
          local%f13v(i,j)=global%f13v(ig,jg)
          local%f13m(i,j)=global%f13m(ig,jg)
          local%f23t(i,j)=global%f23t(ig,jg)
          local%f23u(i,j)=global%f23u(ig,jg)
          local%f23v(i,j)=global%f23v(ig,jg)
          local%f23m(i,j)=global%f23m(ig,jg)
          local%dxt(i,j)=global%dxt(ig,jg)
          local%dxu(i,j)=global%dxu(ig,jg)
          local%dxv(i,j)=global%dxv(ig,jg)
          local%dxm(i,j)=global%dxm(ig,jg)
          local%dyt(i,j)=global%dyt(ig,jg)
          local%dyu(i,j)=global%dyu(ig,jg)
          local%dyv(i,j)=global%dyv(ig,jg)
          local%dym(i,j)=global%dym(ig,jg)
          local%fmapt(i,j)=global%fmapt(ig,jg)
          local%fmapu(i,j)=global%fmapu(ig,jg)
          local%fmapv(i,j)=global%fmapv(ig,jg)
          local%fmapm(i,j)=global%fmapm(ig,jg)
          local%fmapti(i,j)=global%fmapti(ig,jg)
          local%fmapui(i,j)=global%fmapui(ig,jg)
          local%fmapvi(i,j)=global%fmapvi(ig,jg)
          local%fmapmi(i,j)=global%fmapmi(ig,jg)
          local%glat(i,j)=global%glat(ig,jg)
          local%glon(i,j)=global%glon(ig,jg)
!!$          local%topzo(i,j)=global%topzo(ig,jg)
          local%lpu(i,j)=global%lpu(ig,jg)
          local%lpv(i,j)=global%lpv(ig,jg)
          local%lpw(i,j)=global%lpw(ig,jg)
       end do
    end do

    if (associated(local%aru) .and. associated(global%aru)) then
       do j = 1, nyp
          jg=j+j0
          do i = 1, nxp
             ig = i+i0
             do k = 1, nzp
                local%aru(k,i,j)=global%aru(k,ig,jg)
             end do
          end do
       end do
    end if

    if (associated(local%arv) .and. associated(global%arv)) then
       do j = 1, nyp
          jg=j+j0
          do i = 1, nxp
             ig = i+i0
             do k = 1, nzp
                local%arv(k,i,j)=global%arv(k,ig,jg)
             end do
          end do
       end do
    end if

    if (associated(local%arw) .and. associated(global%arw)) then
       do j = 1, nyp
          jg=j+j0
          do i = 1, nxp
             ig = i+i0
             do k = 1, nzp
                local%arw(k,i,j)=global%arw(k,ig,jg)
             end do
          end do
       end do
    end if

    if (associated(local%volu) .and. associated(global%volu)) then
       do j = 1, nyp
          jg=j+j0
          do i = 1, nxp
             ig = i+i0
             do k = 1, nzp
                local%volu(k,i,j)=global%volu(k,ig,jg)
             end do
          end do
       end do
    end if

    if (associated(local%volv) .and. associated(global%volv)) then
       do j = 1, nyp
          jg=j+j0
          do i = 1, nxp
             ig = i+i0
             do k = 1, nzp
                local%volv(k,i,j)=global%volv(k,ig,jg)
             end do
          end do
       end do
    end if

    if (associated(local%volw) .and. associated(global%volw)) then
       do j = 1, nyp
          jg=j+j0
          do i = 1, nxp
             ig = i+i0
             do k = 1, nzp
                local%volw(k,i,j)=global%volw(k,ig,jg)
             end do
          end do
       end do
    end if

    if (associated(local%volt) .and. associated(global%volt)) then
       do j = 1, nyp
          jg=j+j0
          do i = 1, nxp
             ig = i+i0
             do k = 1, nzp
                local%volt(k,i,j)=global%volt(k,ig,jg)
             end do
          end do
       end do
    end if
  end subroutine ExtractLocalFromGlobal




  subroutine get_akmin2d(ngr, n2, n3, akmin2d, &
       mynum, nodei0, nodej0,akmin,topo)

    implicit none
    ! Arguments:
    integer, intent(IN) :: n2, n3, ngr
    real, intent(OUT)   :: akmin2d(n2,n3)
    real, intent(IN)   :: topo(n2,n3)
    integer, intent(in) :: mynum
    integer, intent(in) :: nodei0(:,:)
    integer, intent(in) :: nodej0(:,:)
    real, intent(in) :: akmin
    ! Local Variables:
    real    :: rlat(n2,n3), rlon(n2,n3)
    integer :: i, j

    !srf-define diferentes AKMINs para melhorar estabilidade
    !srf-sobre os Andes
    akmin2d = abs(akmin)

    if(akmin == -1.0) then 
       akmin2d = 1.0
       return
    endif


    !----
    !-calculate lat, lon of each grid box T-points
    do j=1,n3
       do i=1,n2
          call xy_ll(rlat(i,j), rlon(i,j), platn(ngr), plonn(ngr), &
               xmn(i+nodei0(mynum,ngr),ngr), ymn(j+nodej0(mynum,ngr),ngr))
       enddo
    enddo

  !srf- versao 4/11/2013 para filtrar chuva espúria nos andes.


     do j=1,n3
       do i=1,n2
!- regiao 1
          if (rlat(i,j)<-15. .and. rlon(i,j)<-60.) then
              if(topo(i,j) > 500. )then
                akmin2d(i,j) = 2.
             endif
          endif
!- regiao 2
          if (rlat(i,j)<-10. .and. rlat(i,j) >= -15.) then
	    if(rlon(i,j)<-62.) then
              if(topo(i,j) > 500.) then
                akmin2d(i,j) = 2.
              endif
            endif
          endif
!- regiao 3
          if (rlat(i,j)> -10.) then
	    if(rlon(i,j)<-69.) then
              if(topo(i,j) > 500.) then
                akmin2d(i,j) = 2.
              endif
            endif
          endif
        enddo
    enddo

    return

!- versao anterior

    do j=1,n3
       do i=1,n2

          if (rlat(i,j)<-15.) then
             if (rlon(i,j)<-60.) then
                !       TOPTWVL_2d(i,j)=15.
!!$              extra2d(5,ngr)%d2(i,j) = 3.
                akmin2d(i,j) = 2.
             endif
          endif

          if (rlat(i,j)>=-15. .and. rlat(i,j)<-9.) then
             if (rlon(i,j)<-62.) then
                !       TOPTWVL_2d(i,j)=15.
!!$              extra2d(5,ngr)%d2(i,j) = 3.
                akmin2d(i,j) = 2.
             endif
          endif

          if (rlat(i,j)>=-9. .and. rlat(i,j)<-1.) then
             if (rlon(i,j)<-70.) then
                !      TOPTWVL_2d(i,j)=15.
!!$              extra2d(5,ngr)%d2(i,j) = 3.
                akmin2d(i,j) = 2.
             endif
          endif

          if (rlat(i,j)>=-1.) then
             if (rlon(i,j)<=-57 .and. rlon(i,j)>-67.) then
                !      TOPTWVL_2d(i,j)=11.
!!$              extra2d(5,ngr)%d2(i,j) = 1.5
                akmin2d(i,j) = 1.5
             endif
          endif

          if (rlat(i,j)>=-1.) then
             if (rlon(i,j)<=-67.) then
                !       TOPTWVL_2d(i,j)=15.
!!$              extra2d(5,ngr)%d2(i,j) = 3.
                akmin2d(i,j) = 2.
             endif
          endif

          if (rlat(i,j)<=-10.) then
             if (rlon(i,j)>=40.) then
                !      TOPTWVL_2d(i,j)=11.
!!$              extra2d(5,ngr)%d2(i,j) = 1.5
                akmin2d(i,j) = 1.5
             endif
          endif
       enddo
    enddo

!!$  print *, "DEBUG-ALF:get_akmin2d:sum(akmin2d),media=", &
!!$       sum(akmin2d), sum(akmin2d)/(n2*n3)
!!$  call flush(6)

  end subroutine get_akmin2d



  subroutine StoreNamelistFileAtMem_grid(oneNamelistFile)
    type(namelistFile), pointer :: oneNamelistFile
    centlat = oneNamelistFile%centlat
    centlon = oneNamelistFile%centlon
    cphas = oneNamelistFile%cphas
    deltax = oneNamelistFile%deltax
    deltay = oneNamelistFile%deltay
    deltaz = oneNamelistFile%deltaz
    distim = oneNamelistFile%distim
    dtlong = oneNamelistFile%dtlong
    dzmax = oneNamelistFile%dzmax
    dzrat = oneNamelistFile%dzrat
    expnme = oneNamelistFile%expnme
    gridu = oneNamelistFile%gridu
    gridv = oneNamelistFile%gridv
    ibnd = oneNamelistFile%ibnd
    icorflg = oneNamelistFile%icorflg
    dyncore_flag = oneNamelistFile%dyncore_flag
    pd_or_mnt_constraint =  oneNamelistFile%pd_or_mnt_constraint
    order_h = oneNamelistFile%order_h
    order_v = oneNamelistFile%order_v
    vveldamp = oneNamelistFile%vveldamp
    idate1 = oneNamelistFile%idate1
    ideltat = oneNamelistFile%ideltat
    if_adap = oneNamelistFile%if_adap
    ihtran = oneNamelistFile%ihtran
    imonth1 = oneNamelistFile%imonth1
    initial = oneNamelistFile%initial
    itime1 = oneNamelistFile%itime1
    iyear1 = oneNamelistFile%iyear1
    jbnd = oneNamelistFile%jbnd
    lsflg = oneNamelistFile%lsflg
    nacoust = oneNamelistFile%nacoust
    naddsc = oneNamelistFile%naddsc
    nestz1 = oneNamelistFile%nestz1
    nestz2 = oneNamelistFile%nestz2
    nfpt = oneNamelistFile%nfpt
    ngrids = oneNamelistFile%ngrids
    ninest = oneNamelistFile%ninest
    njnest = oneNamelistFile%njnest
    nknest = oneNamelistFile%nknest
    nndtrat = oneNamelistFile%nndtrat
    nnstbot = oneNamelistFile%nnstbot
    nnsttop = oneNamelistFile%nnsttop
    nnxp = oneNamelistFile%nnxp
    nnyp = oneNamelistFile%nnyp
    nnzp = oneNamelistFile%nnzp
    npatch = oneNamelistFile%npatch
    nstratx = oneNamelistFile%nstratx
    nstraty = oneNamelistFile%nstraty
    nstratz1 = oneNamelistFile%nstratz1
    nstratz2 = oneNamelistFile%nstratz2
    nxtnest = oneNamelistFile%nxtnest
    nzg = oneNamelistFile%nzg
    nzs = oneNamelistFile%nzs
    polelat = oneNamelistFile%polelat
    polelon = oneNamelistFile%polelon
    runtype = oneNamelistFile%runtype
    timeunit = oneNamelistFile%timeunit
    timmax = oneNamelistFile%timmax
    zz = oneNamelistFile%zz
  end subroutine StoreNamelistFileAtMem_grid
end module mem_grid
