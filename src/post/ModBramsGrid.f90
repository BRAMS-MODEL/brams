module ModBramsGrid

  use ModParallelEnvironment, only: &
       MsgDump

  use node_mod, only: &
       mynum, &
       nodemxp, &
       nodemyp, &
       nodemzp, &
       nmachs, &
       mchnum, &
       mynum, &
       master_num, &
       ixb, &
       ixe, &
       iyb, &
       iye, &
       nodei0, &
       nodej0, &
       nodeibcon

  use mem_grid, only: &
       nnxp, &
       nnyp, &
       nnzp, &
       time, &
       nzg, &
       nzs, &
       npatch, &
       xtn, &
       deltaxn, &
       ytn, &
       deltayn, &
       ztn, &
       dzmn, &
       dztn, &
       zmn, &
       polelat, &
       polelon, &
       oneGlobalGridData, &
       grid_g

  use mem_aerad, only: &
       nwave

  use ref_sounding, only: &
       pi01dn, &
       th01dn

  use ModNamelistFile, only: namelistFile

  use ModPostUtils, only: UpperCase, DumpFloating, &
       DumpIntegerPairs, DumpRealPairs

  implicit none

  private
  public :: BramsGrid
  public :: CreateBramsGrid
  public :: DestroyBramsGrid
  public :: DumpBramsGrid


  type BramsGrid

     ! domain decomposition and grid independent data

     real :: polelat   ! from RAMSIN
     real :: polelon   ! from RAMSIN

     ! current grid number

     integer :: currGrid     ! current grid

     ! grid dependent data but domain decomposition independent data

     ! size of full lat-lon grid as well as polar stereographic projection

     integer :: nnxp      ! nnxp(currGrid), from RAMSIN
     integer :: nnyp      ! nnyp(currGrid), from RAMSIN
     integer :: nnzp      ! nnzp(currGrid), from RAMSIN
     integer :: nzg       ! nzg from mem_grid
     integer :: nzs       ! nzs from mem_grid
     integer :: npatch    ! npatch from mem_grid
     integer :: nwave     ! nwave from mem_aerad

     ! coordinates (meters) of grid points at the tangent plane

     real, pointer :: xtn(:)    ! xtn(1:nnxp,currGrid) from mem_grid
     real          :: deltax    ! deltaxn(currGrid) from mem_grid
     real, pointer :: ytn(:)    ! ytn(1:nnyp,currGrid) from mem_grid
     real          :: deltay    ! deltayn(currGrid) from mem_grid

     ! latitude and longitude of full grid points at earth's surface

     real, pointer :: glat(:,:) ! oneGlobalGridData(currGrid)%global_glat(nnxp,nnyp) from mem_grid
     real, pointer :: glon(:,:) ! oneGlobalGridData(currGrid)%global_glon(nnxp,nnyp) from mem_grid

     ! fields indexed by vertical (domain decomposition invariant data)

     real, pointer :: pi01dn(:) ! pi01dn(1:nnzp,currGrid) from ref_sounding
     real, pointer :: th01dn(:) ! th01dn(1:nnzp,currGrid) from ref_sounding
     real, pointer :: ztn(:)    ! ztn(1:nnzp,currGrid) from mem_grid (dim nnzp)
     real, pointer :: dzmn(:)   ! dzmn(1:nnzp,currGrid) from mem_grid (dim nnzp)
     real, pointer :: dztn(:)   ! dztn(1:nnzp,currGrid) from mem_grid (dim nnzp)
     real, pointer :: zmn(:)    ! zmn(1:nnzp,currGrid) from mem_grid (dim nnzp)
     real          :: ztop      ! zmn(nnzp-1,currGrid) from mem_grid

     ! grid and domain decomposition dependent data
     ! size of the domain decomposed grids (both lat-lon and polar stereographic)
     ! at this process

     integer :: mxp       ! nodemxp(mynum, currGrid) from node_mod
     integer :: myp       ! nodemyp(mynum, currGrid) from node_mod
     integer :: mzp       ! nodemzp(mynum, currGrid) from node_mod

     ! number of MPI and BRAMS processes

     integer :: nmachs    ! nmachs from node_mod

     ! this process number (MPI and BRAMS)

     integer :: mchnum    ! mchnum from node_mod
     integer :: mynum     ! mynum from node_mod

     ! master number (MPI)

     integer :: master_num ! master_num from node_mod

     ! domain decomposition (indexed by BRAMS process number)

     integer, pointer :: ixb(:) ! ixb(:,currGrid) from node_mod
     integer, pointer :: ixe(:) ! ixe(:,currGrid) from node_mod
     integer, pointer :: iyb(:) ! iyb(:,currGrid) from node_mod
     integer, pointer :: iye(:) ! iye(:,currGrid) from node_mod
     integer, pointer :: nodei0(:) ! nodei0(:,currGrid) from node_mod
     integer, pointer :: nodej0(:) ! nodej0(:,currGrid) from node_mod
     integer, pointer :: nodeibcon(:) ! nodeibcon from node_mod
  end type BramsGrid

  logical, parameter :: dumpLocal=.false.
contains





  subroutine CreateBramsGrid(oneBramsGrid, currGrid)
    type(BramsGrid), pointer :: oneBramsGrid
    integer,        intent(in)    :: currGrid

    integer :: ierr
    character(len=*), parameter :: h="**(CreateBramsGrid)**"

    allocate(oneBramsGrid, stat=ierr)
    if (ierr /= 0) then
       call fatal_error(h//" allocating oneBramsGrid")
    end if

    oneBramsGrid%polelat = polelat
    oneBramsGrid%polelon = polelon

    oneBramsGrid%currGrid = currGrid

    oneBramsGrid%nnxp = nnxp(currGrid)
    oneBramsGrid%nnyp = nnyp(currGrid)
    oneBramsGrid%nnzp = nnzp(currGrid)
    oneBramsGrid%nzg = nzg
    oneBramsGrid%nzs = nzs
    oneBramsGrid%npatch = npatch
    oneBramsGrid%nwave = nwave

    oneBramsGrid%xtn => xtn(1:oneBramsGrid%nnxp,currGrid)
    oneBramsGrid%deltax = deltaxn(currGrid)

    oneBramsGrid%ytn => ytn(1:oneBramsGrid%nnyp,currGrid)
    oneBramsGrid%deltay = deltayn(currGrid)

    oneBramsGrid%glat => oneGlobalGridData(currGrid)%global_glat

    oneBramsGrid%glon => oneGlobalGridData(currGrid)%global_glon

    oneBramsGrid%pi01dn => pi01dn(1:oneBramsGrid%nnzp,currGrid)

    oneBramsGrid%th01dn => th01dn(1:oneBramsGrid%nnzp,currGrid)

    oneBramsGrid%ztn => ztn(1:oneBramsGrid%nnzp,currGrid)

    oneBramsGrid%dzmn => dzmn(1:oneBramsGrid%nnzp,currGrid)
    oneBramsGrid%dztn => dztn(1:oneBramsGrid%nnzp,currGrid)

    oneBramsGrid%zmn => zmn(1:oneBramsGrid%nnzp,currGrid)

    oneBramsGrid%ztop = zmn(oneBramsGrid%nnzp-1, currGrid)

    oneBramsGrid%mxp = nodemxp(mynum, currGrid)
    oneBramsGrid%myp = nodemyp(mynum, currGrid)
    oneBramsGrid%mzp = nodemzp(mynum, currGrid)


    oneBramsGrid%nmachs = nmachs
    oneBramsGrid%mchnum = mchnum
    oneBramsGrid%mynum = mynum
    oneBramsGrid%master_num = master_num

    oneBramsGrid%ixb => ixb(:,currGrid)
    oneBramsGrid%ixe => ixe(:,currGrid)
    oneBramsGrid%iyb => iyb(:,currGrid)
    oneBramsGrid%iye => iye(:,currGrid)
    oneBramsGrid%nodei0 => nodei0(:,currGrid)
    oneBramsGrid%nodej0 => nodej0(:,currGrid)
    oneBramsGrid%nodeibcon => nodeibcon(:,currGrid)

    if (dumpLocal) then
       call DumpBramsGrid(oneBramsGrid)
    end if
  end subroutine CreateBramsGrid





  subroutine DestroyBramsGrid(oneBramsGrid)
    type(BramsGrid), pointer :: oneBramsGrid

    integer :: ierr
    character(len=*), parameter :: h="**(DestroyBramsGrid)**"

    ! return if null oneBramsGrid

    if (.not. associated(oneBramsGrid)) then
       return
    end if

    ! deallocate each allocatable; nullify each pointer

    oneBramsGrid%xtn => null()
    oneBramsGrid%ytn => null()
    oneBramsGrid%glat => null()
    oneBramsGrid%glon => null()
    oneBramsGrid%pi01dn => null()
    oneBramsGrid%th01dn => null()
    oneBramsGrid%ztn => null()
    oneBramsGrid%dzmn => null()
    oneBramsGrid%dztn => null()
    oneBramsGrid%zmn => null()
    oneBramsGrid%ixb => null()
    oneBramsGrid%ixe => null()
    oneBramsGrid%iyb => null()
    oneBramsGrid%iye => null()
    oneBramsGrid%nodei0 => null()
    oneBramsGrid%nodej0 => null()
    oneBramsGrid%nodeibcon => null()

    deallocate(oneBramsGrid, stat=ierr)
    if (ierr /= 0) then
       call fatal_error(h//" deallocating oneBramsGrid")
    end if
    oneBramsGrid => null()
  end subroutine DestroyBramsGrid


  

  subroutine DumpBramsGrid(oneBramsGrid)
    type(BramsGrid), pointer :: oneBramsGrid

    integer :: i
    integer :: iy
    integer :: istart
    character(len=8) :: c0, c1, c2
    character(len=16) :: d0, d1
    character(len=80) :: line
    character(len=*), parameter :: h="**(DumpBramsGrid)**"

    call MsgDump (h//" domain decomposition and grid independent data:")
    write(d0,"(e16.7)") oneBramsGrid%polelat
    write(d1,"(e16.7)") oneBramsGrid%polelon
    call MsgDump (h//" polelat="//trim(adjustl(d0))//&
         "; polelon="//trim(adjustl(d1)))

    write(c0,"(i8)") oneBramsGrid%currGrid
    call MsgDump (h//" domain decomposition independent data of grid number "//trim(adjustl(c0)))

    write(c0,"(i8)") oneBramsGrid%nnxp
    write(c1,"(i8)") oneBramsGrid%nnyp
    call MsgDump (h//" grid with number of (longitudes,latitudes) = ("//&
         trim(adjustl(c0))//","//trim(adjustl(c1))//") also (x,y)")
    write(c0,"(i8)") oneBramsGrid%nnzp
    write(c1,"(i8)") oneBramsGrid%nzg
    write(c2,"(i8)") oneBramsGrid%nzs
    call MsgDump (h// " number of verticals layers: atmosphere = "//&
         trim(adjustl(c0))//"; soil = "//trim(adjustl(c1))//"; snow = "//&
         trim(adjustl(c2)))
    write(c0,"(i8)") oneBramsGrid%npatch
    write(c1,"(i8)") oneBramsGrid%nwave
    call MsgDump (h//" number of patches = "//&
         trim(adjustl(c0))//"; waves = "//trim(adjustl(c1)))

    call MsgDump (h//" x (longitudes) axis coordinates (meters) at tangent plane")
    call DumpFloating(h, oneBramsGrid%xtn)

    call MsgDump (h//" y (latitudes) axis coordinates (meters) at tangent plane")
    call DumpFloating(h, oneBramsGrid%ytn)

    call MsgDump (h//" data indexed by verticals: pi01dn")
    call DumpFloating(h, oneBramsGrid%pi01dn)

    call MsgDump (h//" data indexed by verticals: th01dn")
    call DumpFloating(h, oneBramsGrid%th01dn)

    call MsgDump (h//" data indexed by verticals: ztn")
    call DumpFloating(h, oneBramsGrid%ztn)

    call MsgDump (h//" data indexed by verticals: dzmn")
    call DumpFloating(h, oneBramsGrid%dzmn)

    call MsgDump (h//" data indexed by verticals: dztn")
    call DumpFloating(h, oneBramsGrid%dztn)

    call MsgDump (h//" data indexed by verticals: zmn")
    call DumpFloating(h, oneBramsGrid%zmn)

    write(d0,"(e16.7)") oneBramsGrid%ztop
    write(*,"(a)") h//" vertical top layer is "//trim(adjustl(d0))

    call MsgDump (h//" (lon,lat) of full "//&
         "earth surface grid points")
    do iy = 1, oneBramsGrid%nnyp
       write(c0,"(i8)") iy
       call DumpRealPairs(h//" lat="//trim(adjustl(c0)), oneBramsGrid%glon(:,iy), oneBramsGrid%glat(:,iy))
    end do

    write(c0,"(i8)") oneBramsGrid%nmachs
    write(c1,"(i8)") oneBramsGrid%mchnum
    write(c2,"(i8)") oneBramsGrid%master_num
    call MsgDump (h//" this is "//trim(adjustl(c1))//" MPI process out of "//&
         trim(adjustl(c0))//"; "//trim(adjustl(c2))//" is the MPI master process")
    write(c0,"(i8)") mynum
    call MsgDump (h//" this is RAMS process number "//trim(adjustl(c0)))

    write(c0,"(i8)") oneBramsGrid%mxp
    write(c1,"(i8)") oneBramsGrid%myp
    write(c2,"(i8)") oneBramsGrid%mzp
    call MsgDump (h//" size of domain decomposed grid at this process is ("//&
         trim(adjustl(c0))//","//trim(adjustl(c1))//","//trim(adjustl(c2))//")")


    call MsgDump (h//" domain decomposition offsets are")
    call DumpIntegerPairs(h, oneBramsGrid%nodei0, oneBramsGrid%nodej0)
  end subroutine DumpBramsGrid
end module ModBramsGrid
