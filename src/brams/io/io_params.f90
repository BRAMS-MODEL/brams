!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module io_params

  use ModNamelistFile, only: namelistFile

  use grid_dims, only: &
       maxgrds,        &
       maxlite,        &
       maxsstfiles,    &
       maxndvifiles

!--(DMK-CCATT-INI)------------------------------------------------------------------
 ! use mem_chem1, only: &
 !      nsrc,           &
 !      maxsrcfiles,    &
 !      max_ntimes_src
!--(DMK-CCATT-INI)------------------------------------------------------------------

  implicit none
  private
  public :: StoreNamelistFileAtIo_Params

 !io_params have ioutput var 
 !mem_chem1 needs ioutput var in order to output recycle vars on vfm.
 !io_params needs mem_chem1 parameters.
 
  integer, parameter :: nsrc=4  !number_sources - cyclic reference error.
  INTEGER, PARAMETER :: maxsrcfiles   = 1500 !cyclic reference error.
    integer, parameter :: max_ntimes_src = 2 !cyclic reference error.


  include "files.h"

  character(len=32), public  :: lite_vars(maxlite) ! from RAMSIN

  character(len=80), public  :: afilin !**(JP)** unused

  character(len=f_name_length), public :: hfilout ! from RAMSIN
  character(len=f_name_length), public :: afilout ! from RAMSIN
  character(len=f_name_length), public :: pastfn ! from RAMSIN
  character(len=f_name_length), public :: hfilin ! from RAMSIN

  character(len=20), public  :: xlite ! from RAMSIN
  character(len=20), public  :: ylite ! from RAMSIN
  character(len=20), public  :: zlite ! from RAMSIN

  integer, public :: ipastin ! from RAMSIN
  integer, public :: ioutput ! from RAMSIN
  integer, public :: ipos    ! from RAMSIN
  integer, public :: iinput
  integer, public :: iopunt
  integer, public :: kwrite ! from RAMSIN
  integer, public :: ihistdel ! from RAMSIN
  integer, public :: iclobber ! from RAMSIN
  integer, public :: nlite_vars ! from RAMSIN

  real, public :: frqhis ! from RAMSIN
  real, public :: frqanl ! from RAMSIN
  real, public :: timstr ! from RAMSIN
  real, public :: avgtim ! from RAMSIN
  real, public :: frqlite ! from RAMSIN
  real, public :: frqmean ! from RAMSIN
  real, public :: frqboth ! from RAMSIN

  integer, public :: itoptflg(maxgrds) ! from RAMSIN
  integer, public :: isstflg(maxgrds) ! from RAMSIN
  integer, public :: ivegtflg(maxgrds) ! from RAMSIN
  integer, public :: isoilflg(maxgrds) ! from RAMSIN
  integer, public :: ndviflg(maxgrds) ! from RAMSIN
  integer, public :: nofilflg(maxgrds) ! from RAMSIN
  integer, public :: itopsflg(maxgrds) ! from RAMSIN
  integer, public :: iz0flg(maxgrds) ! from RAMSIN
  integer, public :: ifusflg(maxgrds) ! from RAMSIN

  real, public  :: z0fact ! from RAMSIN
  integer, public :: ntopsmth
  integer, public :: izflat
  real, public  :: z0max(maxgrds) ! from RAMSIN
  real, public  :: toptenh(maxgrds) ! from RAMSIN
  real, public  :: toptwvl(maxgrds) ! from RAMSIN

  character(len=f_name_length), public :: itoptfn(maxgrds) ! from RAMSIN
  character(len=f_name_length), public :: isstfn(maxgrds) ! from RAMSIN
  character(len=f_name_length), public :: ivegtfn(maxgrds) ! from RAMSIN
  character(len=f_name_length), public :: isoilfn(maxgrds) ! from RAMSIN
  character(len=f_name_length), public :: ndvifn(maxgrds) ! from RAMSIN
  character(len=f_name_length), public :: ifusfn(maxgrds) ! from RAMSIN


  character(len=f_name_length), public :: sfcfiles ! from RAMSIN
  character(len=f_name_length), public :: topfiles ! from RAMSIN
  character(len=f_name_length), public :: fusfiles ! from RAMSIN
!!$  character(len=14), public  :: lastdate_sst(maxgrds)

  ! sst files data structure

  character(len=f_name_length), public :: sstfpfx ! from RAMSIN (sst file prefix)
  integer, public            :: iupdsst ! from RAMSIN (0=no update during integration; 1=update)
  integer, public            :: isstcycdata  ! cyclic (on time) files; set by SstFileInv
  integer, public            :: isstcyclic   ! cyclic files and data will be time updated; set by SstFileInv

!!$  integer, public            :: nsstfiles(maxgrds)                ! number of files for each grid; set by SstFileInv
  integer, allocatable, public :: nsstfiles(:)                ! number of files for each grid; set by SstFileInv

!!$  character(len=256), public :: fnames_sst(maxsstfiles,maxgrds)   ! file names for each grid; set by SstFileInv
  character(len=f_name_length), allocatable, public :: fnames_sst(:,:)   ! file names for each grid; set by SstFileInv

!!$  character(len=14), public  :: itotdate_sst(maxsstfiles,maxgrds) ! file dates for each grid, ordered by increasing date; set by SstFileInv
  character(len=14), allocatable, public  :: itotdate_sst(:,:) ! file dates for each grid, ordered by increasing date; set by SstFileInv

!!$  integer, public            :: isstflp(maxgrds)  ! index of last file prior to (or at) starting date; set by SstFileInv
!!$  integer, public            :: isstflf(maxgrds)  ! index of first file later than starting date; set by SstFileInv
!!$  real, public               :: ssttime1(maxgrds)
!!$  real, public               :: ssttime2(maxgrds)
  integer, allocatable, public            :: isstflp(:)  ! index of last file prior to (or at) starting date; set by SstFileInv
  integer, allocatable, public            :: isstflf(:)  ! index of first file later than starting date; set by SstFileInv
  real, allocatable, public               :: ssttime1(:)
  real, allocatable, public               :: ssttime2(:)

  ! ndvi files data structure

  character(len=f_name_length), public :: ndvifpfx ! from RAMSIN (ndvi file prefix)
  integer, public            :: iupdndvi ! from RAMSIN (0=no update during integration; 1=update)
  integer, public            :: indvicyclic ! cyclic (on time) files; set by NdviFileInv
  integer, public            :: indvicycdata ! cyclic files and data will be time updated; set by NdviFileInv

!!$  integer, public            :: nndvifiles(maxgrds)                 ! number of files for each grid; set by NdviFileInv
!!$  character(len=256), public :: fnames_ndvi(maxndvifiles,maxgrds)   ! file names for each grid; set by NdviFileInv
!!$  character(len=14), public  :: itotdate_ndvi(maxndvifiles,maxgrds) ! file dates for each grid, ordered by increasing date; set by NdviFileInv
!!$  integer, public            :: indviflp(maxgrds)  ! index of last file prior to (or at) starting date; set by NdviFileInv
!!$  integer, public            :: indviflf(maxgrds)  ! index of first file later than starting date; set by NdviFileInv
!!$  real, public               :: ndvitime1(maxgrds)
!!$  real, public               :: ndvitime2(maxgrds)
  integer, allocatable, public            :: nndvifiles(:)                 ! number of files for each grid; set by NdviFileInv
  character(len=f_name_length), allocatable, public :: fnames_ndvi(:,:)   ! file names for each grid; set by NdviFileInv
  character(len=14), allocatable, public  :: itotdate_ndvi(:,:) ! file dates for each grid, ordered by increasing date; set by NdviFileInv
  integer, allocatable, public            :: indviflp(:)  ! index of last file prior to (or at) starting date; set by NdviFileInv
  integer, allocatable, public            :: indviflf(:)  ! index of first file later than starting date; set by NdviFileInv
  real, allocatable, public               :: ndvitime1(:)
  real, allocatable, public               :: ndvitime2(:)

!!$  character(len=8), public   :: plfmt(50)
!!$  character(len=8), public   :: pltit(50)
  character(len=16), public  :: iplfld(50) ! from RAMSIN
  integer, public            :: nplt ! from RAMSIN
  integer, public            :: initfld ! from RAMSIN
  integer, public            :: ixsctn(50) ! from RAMSIN
!!$  integer, public            :: iplvect(50)
  integer, public            :: isbval(50) ! from RAMSIN
!!$  integer, public            :: iaa(50)
!!$  integer, public            :: iab(50)
!!$  integer, public            :: joa(50)
!!$  integer, public            :: job(50)
!!$  integer, public            :: naavg(50)
!!$  integer, public            :: noavg(50)

  real, public             :: frqprt ! from RAMSIN
!!$  real, public             :: plconlo(50)
!!$  real, public             :: plconhi(50)
!!$  real, public             :: plconin(50)

!--(DMK-CCATT-INI)------------------------------------------------------------------
  ! emission source files data structure

  real,                            public :: srctime1=0.
  real,                            public :: srctime2=0.
  character(len=256), allocatable, public :: fnames_src(:,:)     
  character(len=14),  allocatable, public :: itotdate_src(:,:)
  real,               allocatable, public :: src_times(:,:)      
  integer,            allocatable, public :: actual_time_index(:,:)
  integer,            allocatable, public :: nsrcfiles(:)
  integer,            allocatable, public :: next_srcfile(:)
!--(DMK-CCATT-FIM)------------------------------------------------------------------

  ! flag to print wall cpu times for each slave node

  integer, public                :: prtcputime ! from RAMSIN

  ! Subroutines:
  public :: createIoData
  public :: destroyIoData
contains

  subroutine createIoData(ngrids)
    implicit none
    ! Arguments:
    integer, intent(in) :: ngrids
    ! Local Variables:
    integer :: ierr

    allocate(nsstfiles(ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating nsstfiles (createIoData)")
    allocate(fnames_sst(maxsstfiles,ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating fnames_sst (createIoData)")
    allocate(itotdate_sst(maxsstfiles,ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR allocating itotdate_sst (createIoData)")
    allocate(isstflp(ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating isstflp (createIoData)")
    allocate(isstflf(ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating isstflf (createIoData)")
    allocate(ssttime1(ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating ssttime1 (createIoData)")
    allocate(ssttime2(ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating ssttime2 (createIoData)")

    allocate(nndvifiles(ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating nndvifiles (createIoData)")
    allocate(fnames_ndvi(maxndvifiles,ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR allocating fnames_ndvi (createIoData)")
    allocate(itotdate_ndvi(maxndvifiles,ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR allocating itotdate_ndvi (createIoData)")
    allocate(indviflp(ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating indviflp (createIoData)")
    allocate(indviflf(ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating indviflf (createIoData)")
    allocate(ndvitime1(ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating ndvitime1 (createIoData)")
    allocate(ndvitime2(ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating ndvitime2 (createIoData)")

!--(DMK-CCATT-INI)------------------------------------------------------------------
    allocate(nsrcfiles(ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating nsrcfiles (createIoData)")
    allocate(fnames_src(maxsrcfiles,ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating fnames_src (createIoData)")
    allocate(itotdate_src(maxsrcfiles,ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR allocating itotdate_src (createIoData)")
    allocate(next_srcfile(ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating next_srcfile (createIoData)")
    allocate(src_times(maxsrcfiles,ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating src_times (createIoData)")
    allocate(actual_time_index(max_ntimes_src,nsrc), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating actual_time_index (createIoData)")
!--(DMK-CCATT-FIM)------------------------------------------------------------------

  end subroutine createIoData

  ! **************************************************************************

  subroutine destroyIoData()
    implicit none
    ! Local Variables:
    integer :: ierr

    deallocate(nsstfiles, STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR deallocating nsstfiles (destroyIoData)")
    deallocate(fnames_sst, STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR deallocating fnames_sst (destroyIoData)")
    deallocate(itotdate_sst, STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR deallocating itotdate_sst (destroyIoData)")
    deallocate(isstflp, STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR deallocating isstflp (destroyIoData)")
    deallocate(isstflf, STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR deallocating isstflp (destroyIoData)")
    deallocate(ssttime1, STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR deallocating ssttime1 (destroyIoData)")
    deallocate(ssttime2, STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR deallocating ssttime2 (destroyIoData)")

    deallocate(nndvifiles, STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR deallocating nndvifiles (destroyIoData)")
    deallocate(fnames_ndvi, STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR deallocating fnames_ndvi (destroyIoData)")
    deallocate(itotdate_ndvi, STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR deallocating itotdate_ndvi (destroyIoData)")
    deallocate(indviflp, STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR deallocating indviflp (destroyIoData)")
    deallocate(indviflf, STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR deallocating indviflf (destroyIoData)")
    deallocate(ndvitime1, STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR deallocating ndvitime1 (destroyIoData)")
    deallocate(ndvitime2, STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR deallocating ndvitime2 (destroyIoData)")

!--(DMK-CCATT-INI)------------------------------------------------------------------
   deallocate(nsrcfiles, STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR deallocating nsrcfiles (destroyIoData)")
   deallocate(fnames_src, STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR deallocating fnames_src (destroyIoData)")
   deallocate(itotdate_src, STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR deallocating itotdate_src (destroyIoData)")
   deallocate(next_srcfile, STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR deallocating next_srcfile (destroyIoData)")
   deallocate(src_times, STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR deallocating src_times (destroyIoData)")
   deallocate(actual_time_index, STAT=ierr)
    if (ierr/=0) call fatal_error&
         ("ERROR deallocating actual_time_index (destroyIoData)")
!--(DMK-CCATT-FIM)------------------------------------------------------------------

  end subroutine destroyIoData








  subroutine StoreNamelistFileAtIo_Params(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    frqboth = oneNamelistFile%frqboth
    afilout = oneNamelistFile%afilout
    avgtim = oneNamelistFile%avgtim
    frqanl = oneNamelistFile%frqanl
    frqhis = oneNamelistFile%frqhis
    frqlite = oneNamelistFile%frqlite
    frqmean = oneNamelistFile%frqmean
    frqprt = oneNamelistFile%frqprt
    hfilin = oneNamelistFile%hfilin
    hfilout = oneNamelistFile%hfilout
    iclobber = oneNamelistFile%iclobber
    ihistdel = oneNamelistFile%ihistdel
    initfld = oneNamelistFile%initfld
    prtcputime = oneNamelistFile%prtcputime
    ioutput = oneNamelistFile%ioutput
    ipos = oneNamelistFile%ipos
    ipastin = oneNamelistFile%ipastin
    iplfld = oneNamelistFile%iplfld
    isbval = oneNamelistFile%isbval
    isoilflg = oneNamelistFile%isoilflg
    isoilfn = oneNamelistFile%isoilfn
    isstflg = oneNamelistFile%isstflg
    isstfn = oneNamelistFile%isstfn
    itopsflg = oneNamelistFile%itopsflg
    itoptflg = oneNamelistFile%itoptflg
    itoptfn = oneNamelistFile%itoptfn
    iupdndvi = oneNamelistFile%iupdndvi
    iupdsst = oneNamelistFile%iupdsst
    ivegtflg = oneNamelistFile%ivegtflg
    ivegtfn = oneNamelistFile%ivegtfn
    ixsctn = oneNamelistFile%ixsctn
    iz0flg = oneNamelistFile%iz0flg
    kwrite = oneNamelistFile%kwrite
    lite_vars = oneNamelistFile%lite_vars
    ndviflg = oneNamelistFile%ndviflg
    ndvifn = oneNamelistFile%ndvifn
    ndvifpfx = oneNamelistFile%ndvifpfx
    nlite_vars = oneNamelistFile%nlite_vars
    nofilflg = oneNamelistFile%nofilflg
    nplt = oneNamelistFile%nplt
    pastfn = oneNamelistFile%pastfn
    sfcfiles = oneNamelistFile%sfcfiles
    sstfpfx = oneNamelistFile%sstfpfx
    timstr = oneNamelistFile%timstr
    topfiles = oneNamelistFile%topfiles
    toptenh = oneNamelistFile%toptenh
    toptwvl = oneNamelistFile%toptwvl
    xlite = oneNamelistFile%xlite
    ylite = oneNamelistFile%ylite
    z0fact = oneNamelistFile%z0fact
    z0max = oneNamelistFile%z0max
    zlite = oneNamelistFile%zlite
    ifusflg = oneNamelistFile%ifusflg
    ifusfn = oneNamelistFile%ifusfn
    fusfiles = oneNamelistFile%fusfiles
  end subroutine StoreNamelistFileAtIo_Params
end module io_params
