module ModPostGrid

  use ModParallelEnvironment, only: &
       MsgDump

  use ModNamelistFile, only: namelistFile

  use ModBramsGrid, only: BramsGrid

  use ModPostUtils, only: UpperCase, DumpFixed, &
       DumpFloating, DumpIntegerPairs, DumpInteger, &
       ptransvar, ctransvar

  use ModOutputUtils, only: GetVarFromMemToOutput

  use mem_grid, only: time, timmax, dtlongn

  use modTimeLineFRN, only: &
   isTimeTograds

  use ParLib, only: parf_GatherPostSfc
  use ParLib, only: parf_barrier

  use io_params, only : & ! 
    IPOS
#ifdef cdf
  use ModPostOneFieldNetcdf, only: &
    writeNetCdf2D, &
    writeNetCdf3D
#endif
  use ModPostTypes

  implicit none

  private

  public :: undef
  !public :: PostGrid
  public :: CreatePostGrid
  public :: DestroyPostGrid
  public :: DumpPostGrid

  public :: OpenGradsBinaryFile
  public :: CloseGradsBinaryFile

  public :: OpenGradsControlFile
  public :: FillGradsControlFile
  public :: CloseGradsControlFile

  public :: UpdateVerticals

  public :: OutputGradsField
  interface OutputGradsField
     module procedure OutputGradsField_2D
     module procedure OutputGradsField_3D
     module procedure OutputGradsField_4D
  end interface

  logical, parameter :: dumpLocal=.false.

contains



  ! CreatePostGrid: allocates and fill onePostGrid, except
  !                 grads output file names and info


  subroutine CreatePostGrid(oneNamelistFile, oneBramsGrid, &
       onePostGrid, currGrid)
    type(NamelistFile), pointer :: oneNamelistFile
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid
    integer, intent(in) :: currGrid

    integer :: ierr
    character :: cgrid
    integer :: iFrq
    character(len=3) :: cUnits
    character(len=8) :: c0
    character(len=*), parameter :: h="**(CreatePostGrid)**"

    ! verify input parameters

    if (.not. associated(oneNamelistFile)) then
       call fatal_error(h//" invoked with null oneNamelistFile")
    else if (.not. associated(oneBramsGrid)) then
       call fatal_error(h//" invoked with null oneBramsGrid")
    else if (associated(onePostGrid)) then
       call fatal_error(h//" invoked with already associated onePostGrid")
    end if

    ! allocate onePostGrid area

    allocate(onePostGrid, stat=ierr)
    if (ierr /= 0) then
       call fatal_error(h//" allocate onePostGrid fails")
    end if

    ! title

    onePostGrid%title = oneNamelistFile%expnme

    ! grid projection

    onePostGrid%project = trim(UpperCase(oneNamelistFile%proj)) == "YES"

    ! fieldID list:
    ! build it

    onePostGrid%buildingFieldIDList = .true.

    ! null fields

    onePostGrid%head => null()
    onePostGrid%tail => null()

    onePostGrid%binFileName = ""  ! indicates first use
    onePostGrid%unitBinFile = INT_UNDEF   ! indicates first use
    onePostGrid%lastRec = 0       ! indicates first use

    onePostGrid%ctlFileName = ""  ! indicates first use
    onePostGrid%unitCtlFile = INT_UNDEF   ! indicates first use

    ! grads file options:

    ! one grads file per output iff anl2Gra = "one",
    ! case insensitive

    onePostGrid%oneToOne = &
         trim(UpperCase(oneNamelistFile%anl2gra)) == "ONE"

    ! by default, build a template file whenever one grads
    ! file per output

    onePostGrid%template = .false. ! temporary

    ! chstep for time scale on post control file
    ! frqanl should be in seconds; check if it
    ! represents minutes or hours

    iFrq = oneNamelistFile%frqanl
    cUnits = "sec"
    if (mod(iFrq, 3600) == 0) then
       iFrq = iFrq/3600
       cUnits = "hr "
    else if (mod(iFrq, 60) == 0) then
       iFrq = iFrq/60
       cUnits = "mn "
    end if
    write(onePostGrid%chstep,"(i12)") iFrq
    onePostGrid%chstep = trim(onePostGrid%chstep)//cUnits

    ! fields related to post lat-lon grid:
    ! firstLon, nLon, delLon, lon, &
    ! firstLat, nLat, delLat, lat

    call PostGlobalLatLonGrid(oneNamelistFile, oneBramsGrid, &
         onePostGrid, currGrid)

    ! region of BRAMS grid that will be used
    ! on post grid:
    ! xStart, yStart indices if no projection;
    ! xMap, yMap, weight if projection

    call BramsXYGridToPost(oneNamelistFile, oneBramsGrid, &
         onePostGrid, currGrid)

    ! maps to pack local chunks of one surface of the BRAMS grid used
    ! in one surface of the post grid: packXLocal, packYLocal, localSizeThisProc;
    ! information on sizes of local chunks at every process
    ! and their displacement to gather chunks: localSize, disp;
    ! information to unpack gathered field at master_proc:
    ! unpackMap

    call PackUnpackMaps(oneBramsGrid, onePostGrid)

    ! create vertical components:
    ! in any case: nVert, vertScaleCode, vertScaleValues
    ! when sigma-z is selected: zMap
    ! when pressure levels is selected: topo, pi
    ! when height is selected: topo

    call CreatePostVerticals(oneNamelistFile, oneBramsGrid, &
         onePostGrid)

    if (dumpLocal) then
       call DumpPostGrid(onePostGrid)
    end if
  end subroutine CreatePostGrid




  subroutine OpenGradsBinaryFile(oneNamelistFile, onePostGrid, oneBramsGrid, currGrid)
    type(NamelistFile), pointer :: oneNamelistFile
    type(PostGrid), pointer :: onePostGrid
    type(BramsGrid), pointer :: oneBramsGrid
    integer, intent(in) :: currGrid

    integer, external :: AvailableFileUnit
    character :: cgrid
    integer :: recSize
    integer :: ios
    real :: outputRecord(onePostGrid%nLon*onePostGrid%nLat)
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(OpenGradsBinaryFile)**"

    ! verify input parameters

    if (.not. associated(oneNamelistFile)) then
       call fatal_error(h//" invoked with null oneNamelistFile")
    else if (.not. associated(onePostGrid)) then
       call fatal_error(h//" invoked with null onePostGrid")
    end if

    ! grid in character

    write (cgrid, "(i1)") currGrid

    ! case one file per date and grid

    if (onePostGrid%oneToOne) then

       ! binary file name

       if(oneBramsGrid%mchnum==oneBramsGrid%master_num) then
          CALL makefnam(onePostGrid%binFileName, &
               trim(oneNamelistFile%gprefix), time, &
               oneNamelistFile%iyear1, oneNamelistFile%imonth1,  &
               oneNamelistFile%idate1, oneNamelistFile%itime1*100, &
               'A', 'g'//cgrid, 'gra')

          onePostGrid%unitBinFile = AvailableFileUnit()
          
          ! open binary file
          
          inquire(iolength=recSize) outputRecord
          if (dumpLocal) then
             write(c0,"(i8)") onePostGrid%unitBinFile
             write(c1,"(i8)") recSize
             call MsgDump (h//" at unit "//trim(adjustl(c0))//&
                  " will open file "//trim(onePostGrid%binFileName)//&
                  " with record size "//trim(adjustl(c1)))
          end if

          open(onePostGrid%unitBinFile, &
               file=trim(onePostGrid%binFileName), &
               form="unformatted", access="direct", status="replace", &
               action="write", recl=recSize, iostat=ios)
          if (ios /= 0) then
             call fatal_error(h//" fails opening file "//trim(onePostGrid%binFileName))
          end if

          ! reset record
          
          onePostGrid%lastRec=0
       end if

       if (dumpLocal) then
          if (onePostGrid%unitBinFile == INT_UNDEF) then
             call MsgDump(h//" grads binary file is not opened at this proc")
          else
             write(c0,"(i8)") onePostGrid%unitBinFile
             call MsgDump (h//" opened grads binary file "//&
                  trim(onePostGrid%binFileName)//" at unit "//&
                  trim(adjustl(c0)))
          end if
       end if

    else

       call fatal_error(h//" not ready for a single grads file")

    end if
  end subroutine OpenGradsBinaryFile





  subroutine CloseGradsBinaryFile(onePostGrid, oneBramsGrid)
    type(PostGrid), pointer :: onePostGrid
    type(BramsGrid), pointer :: oneBramsGrid
    character(len=*), parameter :: h="**(CloseGradsBinaryFile)**"

    ! verify input parameters

    if (.not. associated(onePostGrid)) then
       call fatal_error(h//" invoked with null onePostGrid")
    end if

    if(oneBramsGrid%mchnum==oneBramsGrid%master_num) then
       close(onePostGrid%unitBinFile)
    end if
  end subroutine CloseGradsBinaryFile




  subroutine OpenGradsControlFile(oneNamelistFile, onePostGrid, oneBramsGrid, currGrid)
    type(NamelistFile), pointer :: oneNamelistFile
    type(PostGrid), pointer :: onePostGrid
    type(BramsGrid), pointer :: oneBramsGrid
    integer, intent(in) :: currGrid

    character :: cgrid
    integer, external :: AvailableFileUnit
    character(len=8) :: c0
    character(len=*), parameter :: h="**(OpenGradsControlFile)**"

    ! verify input parameters

    if (.not. associated(oneNamelistFile)) then
       call fatal_error(h//" invoked with null oneNamelistFile")
    else if (.not. associated(onePostGrid)) then
       call fatal_error(h//" invoked with null onePostGrid")
    end if

    ! grid in character

    write (cgrid, "(i1)") currGrid

    ! case one file per date and grid

    if (onePostGrid%oneToOne) then

       ! grads control file name

       CALL makefnam(onePostGrid%ctlFileName, trim(oneNamelistFile%gprefix), time, &
            oneNamelistFile%iyear1, oneNamelistFile%imonth1,  &
            oneNamelistFile%idate1, oneNamelistFile%itime1*100, &
            'A', 'g'//cgrid, 'ctl')

       ! grads control file unit

       if(oneBramsGrid%mchnum==oneBramsGrid%master_num) then
          onePostGrid%unitCtlFile = AvailableFileUnit()

          ! open file
          
          open(onePostGrid%unitCtlFile, &
               file=trim(onePostGrid%ctlFileName), &
               status="replace", action="write")
          
          if (dumpLocal) then
             write(c0,"(i8)") onePostGrid%unitCtlFile
             call MsgDump (h//" opened grads control file "//&
                  trim(onePostGrid%binFileName)//" at unit "//&
                  trim(adjustl(c0)))
          end if
       end if
    else

       call fatal_error(h//" not ready for a single grads file")

    end if
  end subroutine OpenGradsControlFile

  subroutine FillGradsControlFile(onePostGrid, oneBramsGrid)
    type(PostGrid), pointer :: onePostGrid
    type(BramsGrid), pointer :: oneBramsGrid

    integer :: lenFName
    integer :: slashIndex
    integer :: izStart

    ! bin fname format is 
    ! <prefix>-A-<4digitYear>-<2digitMonth>-<2digitDay>-<6digitHHMMSS>-g<1digitgrid>.gra
    ! or <prefix>-A-YYYY-MM-DD-HHMMSS-gX.gra
    ! these constants are the number of characters of each field from the string length

    integer, parameter :: afterLastHHMMSS=7
    integer, parameter :: afterFirstHHMMSS=12
    integer, parameter :: afterLastDay=14
    integer, parameter :: afterFirstDay=15
    integer, parameter :: afterLastMonth=17
    integer, parameter :: afterFirstMonth=18
    integer, parameter :: afterLastYear=20
    integer, parameter :: afterFirstYear=23
    character(len=4) :: cYear
    integer :: iMonth
    character(len=2) :: cMonth
    character(len=2) :: cDay
    character(len=6) :: cHHMMSS
    character(len=15) :: chdate
    character(len=3), parameter :: cmo(12) = (/&
         'jan','feb','mar','apr','may','jun', & 
         'jul','aug','sep','oct','nov','dec' /)

    character(len=*), parameter :: h="**(FillGradsControlFile)**"

    ! verify input parameters

    if (.not. associated(onePostGrid)) then
       call fatal_error(h//" invoked with null onePostGrid")
    end if

    ! case one file per date and grid

    if (onePostGrid%oneToOne) then

       ! write binary file name at current directory
       ! extract path from binary file name

       if(oneBramsGrid%mchnum==oneBramsGrid%master_num) then
          lenFName = len_trim(onePostGrid%binFileName)
          slashIndex = index(onePostGrid%binFileName, "/", .true.)
          write(onePostGrid%unitCtlFile, "(a)") "dset ^"//&
               onePostGrid%binFileName(slashIndex+1:lenFName)

          ! write undef and title 

          write(onePostGrid%unitCtlFile, "(a,e15.7)") "undef ",undef
          write(onePostGrid%unitCtlFile, "(a)") "title "//&
               trim(onePostGrid%title)
          
          ! write axis scales: x
          
          write(onePostGrid%unitCtlFile, "(a,i4,a,2f15.7)") &
               "xdef ", onePostGrid%nLon, &
               " linear ", onePostGrid%firstLon, onePostGrid%delLon
          
          ! write axis scales: y
          
          write(onePostGrid%unitCtlFile, "(a,i4,a,2f15.7)") &
               "ydef ", onePostGrid%nLat, &
               " linear ", onePostGrid%firstLat, onePostGrid%delLat
          
          ! write axis scales: z
          
          write(onePostGrid%unitCtlFile, "(a,i4,a)", advance="no") &
               "zdef ", onePostGrid%nVert," levels "
          do izStart = 1, onePostGrid%nVert, 15
             write(onePostGrid%unitCtlFile, "(15f10.1)") &
                  onePostGrid%vertScaleValues(izStart:min(izStart+14,onePostGrid%nVert))
          end do
          
          ! write axis scales: time
          ! a single time output
          ! extract date and time from bin file name
          ! delta 
          
          lenFName = len_trim(onePostGrid%binFileName)
          cYear = onePostGrid%binFileName(lenFName-afterFirstYear:lenFName-afterLastYear)
          read(onePostGrid%binFileName(lenFName-afterFirstMonth:lenFName-afterLastMonth), "(i2)") iMonth
          cDay = onePostGrid%binFileName(lenFName-afterFirstDay:lenFName-afterLastDay)
          cHHMMSS = onePostGrid%binFileName(lenFName-afterFirstHHMMSS:lenFName-afterLastHHMMSS)
          chdate = cHHMMSS(1:2)//":"//cHHMMSS(3:4)//"z"//cDay//cmo(iMonth)//cYear
          if(trim(onePostGrid%chstep)=="mn" .or. trim(onePostGrid%chstep)=="hr" &
            .or.  trim(onePostGrid%chstep)=="dy" .or. trim(onePostGrid%chstep)=="mo") then
               write(onePostGrid%unitCtlFile, "(a4,i4,a8,2a15)") &
                  "tdef ",1," linear ",chdate, onePostGrid%chstep
          else
               write(onePostGrid%unitCtlFile, "(a4,i4,a8,2a15)") &
                  "tdef ",1," linear ",chdate, "1mn"
          endif
          ! post variables
          
          call DumpFieldIDList(onePostGrid%unitCtlFile, onePostGrid)
       end if

    else

       call fatal_error(h//" not ready for a single grads file")

    end if

    ! field ID list is done

    onePostGrid%buildingFieldIDList = .false.

  end subroutine FillGradsControlFile




  subroutine DumpFieldIDList(unit, onePostGrid)
    integer, intent(in) :: unit    ! Fortran i/o unit
    type(postGrid), pointer :: onePostGrid

    integer :: cnt
    type(fieldID), pointer :: curr
    character(len=*), parameter :: mask=" 99    - RAMS : "
    character(len=*), parameter :: h="**(DumpFieldIDList)**"

    if (.not. associated(onePostGrid)) then
       call fatal_error(h//" null onePostGrid")
    end if

    ! count and write number of entries

    if (.not. associated(onePostGrid%head)) then
       call MsgDump (h//" empty FieldIDList")
    end if
    cnt = 0
    curr => onePostGrid%head
    do
       if (.not. associated(curr)) then
          exit
       else
          cnt = cnt + 1
          curr => curr%next
       end if
    end do

    write(unit,"(a5,i4)") 'vars ', cnt

    ! dump each list entry

    curr => onePostGrid%head
    do
       if (.not. associated(curr)) then
          exit
       else
          write(unit, "(a10,i4,a16,a40,a1,a8,a1)") &
               curr%fieldName, curr%nLevCode, mask, &
               curr%fieldDescription, "[",curr%fieldUnits,"]"
          curr => curr%next
       end if
    end do
    write(unit,"(a)") "endvars"
  end subroutine DumpFieldIDList




  subroutine CloseGradsControlFile(onePostGrid, oneBramsGrid)
    type(PostGrid), pointer :: onePostGrid
    type(BramsGrid), pointer :: oneBramsGrid
    character(len=*), parameter :: h="**(CloseGradsControlFile)**"

    ! verify input parameters

    if (.not. associated(onePostGrid)) then
       call fatal_error(h//" invoked with null onePostGrid")
    end if

    if(oneBramsGrid%mchnum==oneBramsGrid%master_num) then
       close(onePostGrid%unitCtlFile)
    end if
  end subroutine CloseGradsControlFile





  subroutine PostGlobalLatLonGrid(oneNamelistFile, oneBramsGrid, &
       onePostGrid, currGrid)
    type(NamelistFile), pointer :: oneNamelistFile
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid
    integer :: currGrid

    integer :: i
    integer :: ierr
    real :: deltaSum
    real :: latMin, latMax, latMinBrams, latMaxBrams
    real :: lonMin, lonMax, lonMinBrams, lonMaxBrams
    character(len=8) :: c0, c1, c2, c3, c4, c5, c6, c7
    character(len=*), parameter :: h="**(PostGlobalLatLonGrid)**"


    ! BRAMS grid is regular at the tangent plane, but it
    ! is not regular at earth's surface (lat-lon) due to the projection.
    ! post grid has to be regular in lat-lon.
    ! this procedure generates the regular post lat-lon grid,
    ! using BRAMS grid as much as possible, subjected to the user
    ! selected projection and grid portion.

    ! Step1: Compute post grid spacing.
    ! post grid will have spacing among points as similar to BRAMS grid 
    ! as possible, to generate a mapping as close to 1 to 1 as possible.

    ! compute average delta longitude over all longitudes:
    ! sum of average delta longitude over all latitudes
    ! divided by number of latitudes

    deltaSum = sum&
         (&
         (oneBramsGrid%glon(oneBramsGrid%nnxp,:) - oneBramsGrid%glon(1,:))&
         /(oneBramsGrid%nnxp-1)&
         )
    onePostGrid%delLon = deltaSum / oneBramsGrid%nnyp

    ! same procedure for latitudes:
    ! sum of average delta latitude over all longitudes
    ! divided by number of longitudes

    deltaSum = sum&
         (&
         (oneBramsGrid%glat(:,oneBramsGrid%nnyp) - oneBramsGrid%glat(:,1))&
         /(oneBramsGrid%nnyp-1)&
         )
    onePostGrid%delLat = deltaSum / oneBramsGrid%nnxp

    ! Step 2: Compute origin and number of points, restricted to
    ! user selected area and projection  specified at namelist file

    if (onePostGrid%project) then

       ! if projection:
       !   find BRAMS grid extremes
       !   intersect with user selected area
       !   define origin and number of points keeping original delta

       ! envelope BRAMS grid area at earth's surface

       lonMinBrams = minval(oneBramsGrid%glon)
       lonMaxBrams = maxval(oneBramsGrid%glon)
       latMinBrams = minval(oneBramsGrid%glat)
       latMaxBrams = maxval(oneBramsGrid%glat)

       ! intersection of BRAMS envelope and user defined area

       lonMin = max(lonMinBrams, oneNamelistFile%loni(currGrid))
       lonMax = min(lonMaxBrams, oneNamelistFile%lonf(currGrid))
       latMin = max(latMinBrams, oneNamelistFile%lati(currGrid))
       latMax = min(latMaxBrams, oneNamelistFile%latf(currGrid))

       ! post grid origin

       onePostGrid%firstLon = lonMin
       onePostGrid%firstLat = latMin

       ! post grid number of points (#intervals + 1)

       onePostGrid%nLon = 1 + &
            ceiling((lonMax-lonMin)/onePostGrid%delLon)
       onePostGrid%nLat = 1 + &
            ceiling((latMax-latMin)/onePostGrid%delLat)

    else

       ! if no projection:
       !   find BRAMS grid origin as average first latitude and average first longitude
       !   compute last latitude and last longitude keeping delta and # points
       !   intersect with user selected area
       !   define origin and number of points keeping original delta

       lonMinBrams = sum(oneBramsGrid%glon(1,:))/oneBramsGrid%nnyp
       latMinBrams = sum(oneBramsGrid%glat(:,1))/oneBramsGrid%nnxp

       lonMaxBrams = lonMinBrams + onePostGrid%delLon*real(oneBramsGrid%nnxp-1)
       latMaxBrams = latMinBrams + onePostGrid%delLat*real(oneBramsGrid%nnyp-1)

       lonMin = max(lonMinBrams, oneNamelistFile%loni(currGrid))
       lonMax = min(lonMaxBrams, oneNamelistFile%lonf(currGrid))
       latMin = max(latMinBrams, oneNamelistFile%lati(currGrid))
       latMax = min(latMaxBrams, oneNamelistFile%latf(currGrid))

       onePostGrid%firstLon = lonMin
       onePostGrid%firstLat = latMin

       if ( (lonMin == lonMinBrams) .and. &
            (lonMax == lonMaxBrams) .and. &
            (latMin == latMinBrams) .and. &
            (latMax == latMaxBrams)         ) then

          onePostGrid%nLon = oneBramsGrid%nnxp
          onePostGrid%nLat = oneBramsGrid%nnyp
       else
          onePostGrid%nLon = 1 + &
               ceiling((lonMax-lonMin)/onePostGrid%delLon)
          onePostGrid%nLat = 1 + &
               ceiling((latMax-latMin)/onePostGrid%delLat)
       end if
    end if

    ! Step 3: Compute longitude and latitude points
    ! given first point, delta and number of points

    allocate(onePostGrid%lon(onePostGrid%nLon), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") onePostGrid%nLon
       call fatal_error(h//" allocate onePostGrid%lon("//&
            trim(adjustl(c0))//") fails")
    end if
    call BuildAxis(onePostGrid%nLon, onePostGrid%firstLon, &
         onePostGrid%delLon, onePostGrid%lon)

    allocate(onePostGrid%lat(onePostGrid%nLat), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") onePostGrid%nLat
       call fatal_error(h//" allocate onePostGrid%lat("//&
            trim(adjustl(c0))//") fails")
    end if
    call BuildAxis(onePostGrid%nLat, onePostGrid%firstLat, &
         onePostGrid%delLat, onePostGrid%lat)


    if (dumpLocal) then
       write(c0,"(i8)") onePostGrid%nLon
       write(c1,"(i8)") onePostGrid%nLat
       write(c2,"(f8.3)") onePostGrid%firstLon
       write(c3,"(f8.3)") onePostGrid%lon(onePostGrid%nLon)
       write(c4,"(f8.3)") onePostGrid%delLon
       write(c5,"(f8.3)") onePostGrid%firstLat
       write(c6,"(f8.3)") onePostGrid%lat(onePostGrid%nLat)
       write(c7,"(f8.3)") onePostGrid%delLat
       call MsgDump (h//" lat-lon Post Grid has size"//&
            " ("//trim(adjustl(c0))//","//trim(adjustl(c1))//")"//&
            " with lon,lat = "//&
            " ("//trim(adjustl(c2))//":"//trim(adjustl(c3))//":"//trim(adjustl(c4))//", "//&
            trim(adjustl(c5))//":"//trim(adjustl(c6))//":"//trim(adjustl(c7))//")")
    end if
  end subroutine PostGlobalLatLonGrid





  subroutine BuildAxis(nVal, firstVal, delVal, val)
    integer, intent(in) :: nVal
    real, intent(in) :: firstVal
    real, intent(in) :: delVal
    real, intent(out) :: val(:)

    integer :: i

    do i = 1, nVal
       val(i) = firstVal + real(i-1)*delVal
    end do
  end subroutine BuildAxis




  subroutine BramsXYGridToPost(oneNamelistFile, oneBramsGrid, &
       onePostGrid, currGrid)
    type(NamelistFile), pointer :: oneNamelistFile
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid
    integer :: currGrid

    integer :: xPost, yPost
    integer :: xBrams, yBrams
    integer :: xEnd, yEnd
    integer :: ierr
    real :: x, y
    real :: distXLow, distYLow
    character(len=8) :: c0, c1, c2, c3
    character(len=16) :: d0, d1, d2
    character(len=*), parameter :: h="**(BramsXYGridToPost)**"


    if (onePostGrid%project) then
       
       ! first case: on projection, post grid points are
       ! interpolated from BRAMS grid vicinity points.
       ! Find mapping from BRAMS to POST and interpolation weights
       
       allocate(onePostGrid%xMap &
            (onePostGrid%nLon,onePostGrid%nLat), &
            stat = ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") onePostGrid%nLon
          write(c1,"(i8)") onePostGrid%nLat
          call fatal_error(h//" fail allocating xMap("//&
               trim(adjustl(c0))//","//trim(adjustl(c1))//")")
       end if
       
       allocate(onePostGrid%yMap &
            (onePostGrid%nLon,onePostGrid%nLat), &
            stat = ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") onePostGrid%nLon
          write(c1,"(i8)") onePostGrid%nLat
          call fatal_error(h//" fail allocating yMap("//&
               trim(adjustl(c0))//","//trim(adjustl(c1))//")")
       end if
       
       allocate(onePostGrid%weight &
            (onePostGrid%nLon,onePostGrid%nLat,4), &
            stat = ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") onePostGrid%nLon
          write(c1,"(i8)") onePostGrid%nLat
          call fatal_error(h//" fail allocating weight("//&
               trim(adjustl(c0))//","//trim(adjustl(c1))//",4)")
       end if
       
       ! set default to undef, map at (1,1)
       
       onePostGrid%xMap = 1
       onePostGrid%yMap = 1
       onePostGrid%weight = undef
       
       ! for every post grid point, find corresponding BRAMS grid cell
       ! and compute weight, according to user selected mean_type (at namelist):

       ! if user selects mean_type="vmp", take the closes brams grid point
       ! if user selects mean_type="bav", take brams points average
       
       select case (UpperCase(oneNamelistFile%mean_type))

       case ("VMP")

          do yPost = 1, onePostGrid%nLat
             do xPost = 1, onePostGrid%nLon
                
                ! project lat-lon post grid point onto BRAMS grid
                
                call ll_xy(onePostGrid%lat(yPost), onePostGrid%lon(xPost), &
                     oneBramsGrid%polelat, oneBramsGrid%polelon, &
                     x, y)
                
                xBrams = floor((x-oneBramsGrid%xtn(1))/oneBramsGrid%deltax) + 1
                yBrams = floor((y-oneBramsGrid%ytn(1))/oneBramsGrid%deltay) + 1
                
                if ( 1 <= xBrams .and. xBrams < oneBramsGrid%nnxp .and. &
                     1 <= yBrams .and. yBrams < oneBramsGrid%nnyp ) then
                   onePostGrid%xMap(xPost,yPost) = xBrams
                   onePostGrid%yMap(xPost,yPost) = yBrams
                   
                   distXLow = (x - oneBramsGrid%xtn(xBrams))/oneBramsGrid%deltax
                   distYLow = (y - oneBramsGrid%ytn(yBrams))/oneBramsGrid%deltay
                   onePostGrid%weight(xPost,yPost,1) = (1.0-distXLow)*(1.0-distYLow)  !(i  , j  )
                   onePostGrid%weight(xPost,yPost,2) = (    distXLow)*(1.0-distYLow)  !(i+1, j  )
                   onePostGrid%weight(xPost,yPost,3) = (1.0-distXLow)*(    distYLow)  !(i  , j+1)
                   onePostGrid%weight(xPost,yPost,4) = (    distXLow)*(    distYLow)  !(i+1, j+1)

                   ! minimize rounding errors: 
                   ! the sum of all weights should be 1.0, but it is usualy not,
                   ! due to rounding. Try to minimize rounding by replacing one
                   ! of the weights by 1-sum(remaining weights)

                   onePostGrid%weight(xPost,yPost,4) = 1.0 - &
                        sum(onePostGrid%weight(xPost,yPost,1:3))

                end if
             end do
          end do

       case ("BAV")

          do yPost = 1, onePostGrid%nLat
             do xPost = 1, onePostGrid%nLon
                
                ! project lat-lon post grid point onto BRAMS grid
                
                call ll_xy(onePostGrid%lat(yPost), onePostGrid%lon(xPost), &
                     oneBramsGrid%polelat, oneBramsGrid%polelon, &
                     x, y)
                
                xBrams = floor((x-oneBramsGrid%xtn(1))/oneBramsGrid%deltax) + 1
                yBrams = floor((y-oneBramsGrid%ytn(1))/oneBramsGrid%deltay) + 1
                
                if ( 1 <= xBrams .and. xBrams < oneBramsGrid%nnxp .and. &
                     1 <= yBrams .and. yBrams < oneBramsGrid%nnyp ) then
                   onePostGrid%xMap(xPost,yPost) = xBrams
                   onePostGrid%yMap(xPost,yPost) = yBrams
                   onePostGrid%weight(xPost,yPost,:) = 0.0
                   
                   distXLow = (x - oneBramsGrid%xtn(xBrams))/oneBramsGrid%deltax
                   distYLow = (y - oneBramsGrid%ytn(yBrams))/oneBramsGrid%deltay

                   if (distXLow > 0.5) then
                      if (distYLow > 0.5) then
                         onePostGrid%weight(xPost,yPost,4) = 1.0  !(i+1, j+1)
                      else                     
                         onePostGrid%weight(xPost,yPost,2) = 1.0  !(i+1, j  )
                      end if
                   else
                      if (distYLow > 0.5) then
                         onePostGrid%weight(xPost,yPost,3) = 1.0  !(i  , j+1)
                      else                     
                         onePostGrid%weight(xPost,yPost,1) = 1.0  !(i  , j  )
                      end if
                   end if
                end if
             end do
          end do

       case default

          call fatal_error(h//" unknown post grid projection option "//&
               trim(UpperCase(oneNamelistFile%mean_type)))
       end select


       if (dumpLocal) then
          write(c0,"(f8.2)") (100.0*count(onePostGrid%weight(:,:,1)==undef))/&
               real(onePostGrid%nLon*onePostGrid%nLat)
          call MsgDump (h//" Before Eliminating Undefined Borders, "//&
               "Post Grid has "//trim(adjustl(c0))//"% of undefs")
          call MsgDump (h//" BRAMS global indices of Post Grib corners are:", .true.)
          write(c2,"(i8)") onePostGrid%xMap(1,1)
          write(c3,"(i8)") onePostGrid%yMap(1,1)
          call MsgDump ("("//trim(adjustl(c2))//","//trim(adjustl(c3))//"); ", .true.)
          write(c2,"(i8)") onePostGrid%xMap(1,onePostGrid%nLat)
          write(c3,"(i8)") onePostGrid%yMap(1,onePostGrid%nLat)
          call MsgDump ("("//trim(adjustl(c2))//","//trim(adjustl(c3))//"); ", .true.)
          write(c2,"(i8)") onePostGrid%xMap(onePostGrid%nLon,1)
          write(c3,"(i8)") onePostGrid%yMap(onePostGrid%nLon,1)
          call MsgDump ("("//trim(adjustl(c2))//","//trim(adjustl(c3))//"); ", .true.)
          write(c2,"(i8)") onePostGrid%xMap(onePostGrid%nLon,onePostGrid%nLat)
          write(c3,"(i8)") onePostGrid%yMap(onePostGrid%nLon,onePostGrid%nLat)
          call MsgDump ("("//trim(adjustl(c2))//","//trim(adjustl(c3))//")")
       end if

       ! Post borders could be fully undef
       ! the next call eliminates borders that are fully undef
       
       call EliminateUndefBorders(onePostGrid, oneBramsGrid)
       
       if (dumpLocal) then
          write(c0,"(f8.2)") (100.0*count(onePostGrid%weight(:,:,1)==undef))/&
               real(onePostGrid%nLon*onePostGrid%nLat)
          call MsgDump (h//" After Eliminating Undefined Borders, "//&
               "Post Grid has "//trim(adjustl(c0))//"% of undefs")
          call MsgDump (h//" BRAMS global indices of Post Grib corners are:", .true.)
          write(c2,"(i8)") onePostGrid%xMap(1,1)
          write(c3,"(i8)") onePostGrid%yMap(1,1)
          call MsgDump ("("//trim(adjustl(c2))//","//trim(adjustl(c3))//"); ", .true.)
          write(c2,"(i8)") onePostGrid%xMap(1,onePostGrid%nLat)
          write(c3,"(i8)") onePostGrid%yMap(1,onePostGrid%nLat)
          call MsgDump ("("//trim(adjustl(c2))//","//trim(adjustl(c3))//"); ", .true.)
          write(c2,"(i8)") onePostGrid%xMap(onePostGrid%nLon,1)
          write(c3,"(i8)") onePostGrid%yMap(onePostGrid%nLon,1)
          call MsgDump ("("//trim(adjustl(c2))//","//trim(adjustl(c3))//"); ", .true.)
          write(c2,"(i8)") onePostGrid%xMap(onePostGrid%nLon,onePostGrid%nLat)
          write(c3,"(i8)") onePostGrid%yMap(onePostGrid%nLon,onePostGrid%nLat)
          call MsgDump ("("//trim(adjustl(c2))//","//trim(adjustl(c3))//")")
       end if

    else if (&
         onePostGrid%nLon == oneBramsGrid%nnxp .and. &
         onePostGrid%nLat == oneBramsGrid%nnyp) then
       
       ! second case: no projection and post grid identical to brams grid
       
       onePostGrid%xStart = 1
       onePostGrid%yStart = 1
       if (dumpLocal) then
          write(c0,"(i8)") onePostGrid%xStart
          write(c1,"(i8)") onePostGrid%yStart
          write(c2,"(i8)") onePostGrid%nLon
          write(c3,"(i8)") onePostGrid%nLat
          call MsgDump (h//" Post Grid indices are BRAMS indices starting at ("//&
               trim(adjustl(c0))//","//trim(adjustl(c1))//") with size ("//&
               trim(adjustl(c2))//","//trim(adjustl(c3))//")")
       end if
       
    else
       
       ! third case: no projection and post grid is a subgrid of brams grid
       ! find indices of first point, which are easier
       ! to find on the tangent plane due to regularity.
       
       ! project first lat-lon post grid point onto BRAMS grid

       call ll_xy(onePostGrid%lat(1), onePostGrid%lon(1), &
            oneBramsGrid%polelat, oneBramsGrid%polelon, &
            x, y)
       
       ! find BRAMS x index corresponding to first and last post x index
       
       onePostGrid%xStart = floor((x-oneBramsGrid%xtn(1))/oneBramsGrid%deltax) + 1
       if ( onePostGrid%xStart < 1 .or. &
            onePostGrid%xStart > oneBramsGrid%nnxp) then
          write(c0,"(f8.3)") onePostGrid%lon(1)
          write(d0,"(e16.7)") x
          write(d1,"(e16.7)") oneBramsGrid%xtn(1)
          write(d2,"(e16.7)") oneBramsGrid%xtn(oneBramsGrid%nnxp)
          call fatal_error(h//" logical error: longitude "//&
               trim(adjustl(c0))//" maps into distance "//&
               trim(adjustl(d0))//" lyes outside BRAMS grid range "//&
               trim(adjustl(d1))//":"//trim(adjustl(d2)))
       end if
       
       xEnd = onePostGrid%xStart + onePostGrid%nLon - 1 
       if (xEnd > oneBramsGrid%nnxp) then
          write(c0,"(i8)") onePostGrid%xStart
          write(c1,"(i8)") onePostGrid%nLon
          write(c2,"(i8)") oneBramsGrid%nnxp
          call fatal_error(h//" logical error: post x axis starting at BRAMS x ="//&
               trim(adjustl(c0))//" and with "//&
               trim(adjustl(c1))//" points exceeds BRAMS number of points ="//&
               trim(adjustl(c2)))
       end if
       
       ! find BRAMS y index corresponding to first and last post y index
       
       onePostGrid%yStart = floor((x-oneBramsGrid%ytn(1))/oneBramsGrid%deltay) + 1
       if ( onePostGrid%yStart < 1 .or. &
            onePostGrid%yStart > oneBramsGrid%nnyp) then
          write(c0,"(f8.3)") onePostGrid%lat(1)
          write(d0,"(e16.7)") y
          write(d1,"(e16.7)") oneBramsGrid%ytn(1)
          write(d2,"(e16.7)") oneBramsGrid%ytn(oneBramsGrid%nnyp)
          call fatal_error(h//" logical error: latitude "//&
               trim(adjustl(c0))//" maps into distance "//&
               trim(adjustl(d0))//" lyes outside BRAMS grid range "//&
               trim(adjustl(d1))//":"//trim(adjustl(d2)))
       end if
       
       yEnd = onePostGrid%yStart + onePostGrid%nLat - 1 
       if (yEnd > oneBramsGrid%nnyp) then
          write(c0,"(i8)") onePostGrid%yStart
          write(c1,"(i8)") onePostGrid%nLat
          write(c2,"(i8)") oneBramsGrid%nnyp
          call fatal_error(h//" logical error: post y axis starting at BRAMS y ="//&
               trim(adjustl(c0))//" and with "//&
               trim(adjustl(c1))//" points exceeds BRAMS number of points ="//&
               trim(adjustl(c2)))
       end if
       
       if (dumpLocal) then
          write(c0,"(i8)") onePostGrid%xStart
          write(c1,"(i8)") onePostGrid%yStart
          write(c2,"(i8)") onePostGrid%nLon
          write(c3,"(i8)") onePostGrid%nLat
          call MsgDump (h//" Post Grid indices are BRAMS indices starting at ("//&
               trim(adjustl(c0))//","//trim(adjustl(c1))//") with size ("//&
               trim(adjustl(c2))//","//trim(adjustl(c3))//")")
       end if
    end if

  end subroutine BramsXYGridToPost




  subroutine EliminateUndefBorders(onePostGrid, oneBramsGrid)
    type(PostGrid), pointer :: onePostGrid
    type(BramsGrid), pointer :: oneBramsGrid

    integer :: firstX, lastX, firstY, lastY
    integer :: nextFirstX, nextLastX, nextFirstY, nextLastY
    integer :: countUndef
    integer :: xSaved(size(onePostGrid%xMap,1),size(onePostGrid%xMap,2))
    integer :: ySaved(size(onePostGrid%yMap,1),size(onePostGrid%yMap,2))
    real    :: wSaved(size(onePostGrid%weight,1),size(onePostGrid%weight,2),4)

    integer :: xPost, yPost
    integer :: xBrams, yBrams
    integer :: ierr
    real :: x, y
    character(len=8) :: c0, c1, c2, c3
    character(len=*), parameter :: h="**(EliminateUndefBorders)**"

    ! eliminate grid borders mostly undef

    firstX=1; lastX=onePostGrid%nLon
    firstY=1; lastY=onePostGrid%nLat

    do

       ! borders (firstX,:), (lastX,:)

       countUndef = count(onePostGrid%weight(firstX      ,firstY:lastY,1)==undef)
       if (2*countUndef > lastY-firstY+1) then
          nextFirstX=firstX+1
       else 
          nextFirstX=firstX
       end if

       countUndef = count(onePostGrid%weight(       lastX,firstY:lastY,1)==undef)
       if (2*countUndef > lastY-firstY+1) then
          nextlastX=lastX-1
       else 
          nextlastX=lastX
       end if

       ! borders (:,firstY), (:,lastY)

       countUndef = count(onePostGrid%weight(firstX:lastX,firstY      ,1)==undef)
       if (2*countUndef > lastX-firstX+1) then
          nextFirstY=firstY+1
       else 
          nextFirstY=firstY
       end if

       countUndef = count(onePostGrid%weight(firstX:lastX,       lastY,1)==undef)
       if (2*countUndef > lastX-firstX+1) then
          nextlastY=lastY-1
       else 
          nextlastY=lastY
       end if

       if ( nextFirstX == firstX .and. &
            nextLastX  == lastX  .and. &
            nextFirstY == firstY .and. &
            nextLastY  == lastY          ) then
          exit
       else
          firstX = nextFirstX
          lastX  = nextLastX
          firstY = nextFirstY
          lastY  = nextLastY
          if (dumpLocal) then
             write(c0,"(i8)") firstX
             write(c1,"(i8)") lastX
             write(c2,"(i8)") firstY
             write(c3,"(i8)") lastY
             call MsgDump (h//" Domain reduced to ("//&
                  trim(adjustl(c0))//":"//trim(adjustl(c1))//", "//&
                  trim(adjustl(c2))//":"//trim(adjustl(c3))//")")
          end if
       end if
    end do


    if ( (firstX /= 1) .or. (lastX /= onePostGrid%nLon) .or. &
         (firstY /= 1) .or. (lastY /= onePostGrid%nLat) ) then

       if (dumpLocal) then
          call MsgDump (h//" axis and maps will be rebuild")
       end if

       ! save current xMap, yMap, weight
 
       xSaved = onePostGrid%xMap
       ySaved = onePostGrid%yMap
       wSaved = onePostGrid%weight
       
       ! deallocate current xMap, yMap, weight

       deallocate(onePostGrid%xMap)
       deallocate(onePostGrid%yMap)
       deallocate(onePostGrid%weight)
       
       ! update nLon, firstLon, lon

       onePostGrid%firstLon = onePostGrid%lon(firstX)
       onePostGrid%nLon = lastX - firstX + 1
       deallocate(onePostGrid%lon)
       allocate(onePostGrid%lon(onePostGrid%nLon))
       call BuildAxis(onePostGrid%nLon, onePostGrid%firstLon, &
            onePostGrid%delLon, onePostGrid%lon)

       ! update nLat, firstLat, lat

       onePostGrid%firstLat = onePostGrid%lat(firstY)
       onePostGrid%nLat = lastY - firstY + 1
       deallocate(onePostGrid%lat)
       allocate(onePostGrid%lat(onePostGrid%nLat))
       call BuildAxis(onePostGrid%nLat, onePostGrid%firstLat, &
            onePostGrid%delLat, onePostGrid%lat)

       ! update xMap, yMap, weight

       allocate(onePostGrid%xMap(onePostGrid%nLon,onePostGrid%nLat))
       allocate(onePostGrid%yMap(onePostGrid%nLon,onePostGrid%nLat))
       allocate(onePostGrid%weight(onePostGrid%nLon,onePostGrid%nLat,4))
       
       do yPost = 1, onePostGrid%nLat
          do xPost = 1, onePostGrid%nLon
             onePostGrid%xMap(xPost,yPost) = xSaved(xPost+firstX-1,yPost+firstY-1)
             onePostGrid%yMap(xPost,yPost) = ySaved(xPost+firstX-1,yPost+firstY-1)
             onePostGrid%weight(xPost,yPost,:) = wSaved(xPost+firstX-1,yPost+firstY-1,:)
             if (dumpLocal) then
                write(*,"(3(a,'(',i2,',',i2,')'))") &
                     h//"grid point ",xPost, yPost,&
                     " was brams point ",xSaved(xPost,yPost),ySaved(xPost,yPost),&
                     " and now is point ",onePostGrid%xMap(xPost,yPost),onePostGrid%yMap(xPost,yPost)
             end if
          end do
       end do
    end if
  end subroutine EliminateUndefBorders




  subroutine PackUnpackMaps(oneBramsGrid, onePostGrid)
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    integer :: proc
    integer :: ix
    integer :: iy
    integer :: ierr
    integer :: xPost, yPost, xBrams, yBrams
    integer :: ownerProc
    integer :: unpackIndexCnt
    integer :: currCnt
    integer :: owner(oneBramsGrid%nnxp, oneBramsGrid%nnyp)
    character(len=8) :: c0, c1, c2, c3, c4, c5, c6
    character(len=*), parameter :: h="**(PackUnpackMaps)**"

    ! which BRAMS process owns each BRAMS grid point

    owner=0
    do proc = 1, oneBramsGrid%nmachs
       do iy = oneBramsGrid%iyb(proc), oneBramsGrid%iye(proc)
          do ix = oneBramsGrid%ixb(proc), oneBramsGrid%ixe(proc)
             owner(ix,iy) = proc
          end do
       end do
       ! west boundary (ix=1)
       if (btest(oneBramsGrid%nodeibcon(proc),0)) then
          do iy = oneBramsGrid%iyb(proc), oneBramsGrid%iye(proc)
             owner(1,iy) = proc
          end do
          ! corner(1,1)
          if (btest(oneBramsGrid%nodeibcon(proc),2)) then
             owner(1,1) = proc
          end if
          ! corner(1,n)
          if (btest(oneBramsGrid%nodeibcon(proc),3)) then
             owner(1,oneBramsGrid%nnyp) = proc
          end if
       end if
       ! east boundary (ix=nnxp)
       if (btest(oneBramsGrid%nodeibcon(proc),1)) then
          do iy = oneBramsGrid%iyb(proc), oneBramsGrid%iye(proc)
             owner(oneBramsGrid%nnxp,iy) = proc
          end do
          ! corner(n,1)
          if (btest(oneBramsGrid%nodeibcon(proc),2)) then
             owner(oneBramsGrid%nnxp,1) = proc
          end if
          ! corner(n,n)
          if (btest(oneBramsGrid%nodeibcon(proc),3)) then
             owner(oneBramsGrid%nnxp,oneBramsGrid%nnyp) = proc
          end if
       end if
       ! south boundary (iy=1)
       if (btest(oneBramsGrid%nodeibcon(proc),2)) then
          do ix = oneBramsGrid%ixb(proc), oneBramsGrid%ixe(proc)
             owner(ix,1) = proc
          end do
       end if
       ! north boundary (iy=nnyp)
       if (btest(oneBramsGrid%nodeibcon(proc),3)) then
          do ix = oneBramsGrid%ixb(proc), oneBramsGrid%ixe(proc)
             owner(ix,oneBramsGrid%nnyp) = proc
          end do
       end if
    end do

    ! check corrrectness of owner

    do iy = 1, oneBramsGrid%nnyp
       do ix = 1, oneBramsGrid%nnxp
          if (owner(ix,iy) == 0) then
             write(c0,"(i8)") ix
             write(c1,"(i8)") iy
             call fatal_error(h//" no BRAMS process owns BRAMS grid point ("//&
                  trim(adjustl(c0))//","//trim(adjustl(c1)))
          end if
       end do
    end do

    ! post grid is spread over BRAMS processes
    ! find post grid size at each process

    allocate(onePostGrid%localSize(oneBramsGrid%nmachs), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") oneBramsGrid%nmachs
       call fatal_error(h//" allocate localSize("//trim(adjustl(c0))//&
            ") fails")
    end if
    onePostGrid%localSize(:) = 0

    if (onePostGrid%project) then

       ! case projection: requires interpolation of BRAMS grid points.
       ! Take only the lower left corner of the interpolation rectangle.
       ! This is the point that will store the interpolated value.
       ! It is garanteed that the entire rectangle lies in the same process.
       ! It may use ghost zone borders that may be outdated.

       do yPost = 1, onePostGrid%nLat
          do xPost = 1, onePostGrid%nLon
             xBrams = onePostGrid%xMap(xPost,yPost)
             yBrams = onePostGrid%yMap(xPost,yPost)
             onePostGrid%localSize(owner(xBrams,yBrams)) = &
                  onePostGrid%localSize(owner(xBrams,yBrams)) + 1
          end do
       end do

    else

       ! case no projection: post grid is a sub-grid of BRAMS grid
       ! starting at BRAMS (xStart,yStart) with size (nLon,nLat).

       do yPost = 1, onePostGrid%nLat
          yBrams = yPost + onePostGrid%yStart - 1
          do xPost = 1, onePostGrid%nLon
             xBrams = xPost + onePostGrid%xStart - 1
             onePostGrid%localSize(owner(xBrams,yBrams)) = &
                  onePostGrid%localSize(owner(xBrams,yBrams)) + 1
          end do
       end do
    end if

    ! save local size of this proc

    onePostGrid%localSizeThisProc=onePostGrid%localSize(oneBramsGrid%mynum)

    ! post grid chunks at every BRAMS process will be gathered
    ! at master_proc for output on a single array;
    ! find the displacements of each chunk first entry
    ! at gathered array

    allocate(onePostGrid%disp(oneBramsGrid%nmachs), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") oneBramsGrid%nmachs
       call fatal_error(h//" allocate disp("//trim(adjustl(c0))//&
            ") fails")
    end if
    onePostGrid%disp(1) = 0
    do proc = 2, oneBramsGrid%nmachs
       onePostGrid%disp(proc) = &
            onePostGrid%disp(proc-1) + &
            onePostGrid%localSize(proc-1)
    end do
    if (dumpLocal) then
       call MsgDump (h//" pairs (localSize, disp) over all processes:")
       call DumpIntegerPairs(h, onePostGrid%localSize, onePostGrid%disp)
    end if

    ! mappings of BRAMS chunk at this process used on post grid onto
    ! a 1D array to be gathered by master_proc;
    ! there are localSizeThisProc entries with local indices (packXLocal, packYLocal);
    ! the packing procedure eliminates unused ghost zones, taking
    ! only BRAMS grid points that are owned by this process.
    ! If the user selects the projection option, interpolation has to
    ! occur at each process before sending the chunk.


    allocate(onePostGrid%packXLocal(onePostGrid%localSizeThisProc), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") onePostGrid%localSizeThisProc
       call fatal_error(h//" allocate packXLocal("//trim(adjustl(c0))//&
            ") fails")
    end if
    allocate(onePostGrid%packYLocal(onePostGrid%localSizeThisProc), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") onePostGrid%localSizeThisProc
       call fatal_error(h//" allocate packYLocal("//trim(adjustl(c0))//&
            ") fails")
    end if
    allocate(onePostGrid%packWeight(onePostGrid%localSizeThisProc,4), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") onePostGrid%localSizeThisProc
       call fatal_error(h//" allocate packWeight("//trim(adjustl(c0))//&
            ",4) fails")
    end if

    allocate(onePostGrid%unpackMap&
         (onePostGrid%nLat*onePostGrid%nLon), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") onePostGrid%nLat*onePostGrid%nLon
       call fatal_error(h//" allocate unpackMap("//trim(adjustl(c0))//&
            ") fails")
    end if


    ! see previous comment on localSize computation

    onePostGrid%localSize(:) = 0
    unpackIndexCnt=0
    do yPost = 1, onePostGrid%nLat
       do xPost = 1, onePostGrid%nLon

          if (onePostGrid%project) then
             xBrams = onePostGrid%xMap(xPost,yPost)
             yBrams = onePostGrid%yMap(xPost,yPost)
          else
             yBrams = yPost + onePostGrid%yStart - 1
             xBrams = xPost + onePostGrid%xStart - 1
          end if

          ownerProc = owner(xBrams,yBrams)

          onePostGrid%localSize(ownerProc) = &
               onePostGrid%localSize(ownerProc) + 1

          if (ownerProc == oneBramsGrid%mynum) then
             currCnt = onePostGrid%localSize(ownerProc) 
             onePostGrid%packXLocal(currCnt) = &
                  xBrams-oneBramsGrid%nodei0(ownerProc)
             onePostGrid%packYLocal(currCnt) = &
                  yBrams-oneBramsGrid%nodej0(ownerProc)

             if (dumpLocal) then
                write(c0,"(i8)") xPost
                write(c1,"(i8)") yPost
                write(c2,"(i8)") xBrams
                write(c3,"(i8)") yBrams
                write(c4,"(i8)") xBrams-oneBramsGrid%nodei0(ownerProc)
                write(c5,"(i8)") yBrams-oneBramsGrid%nodej0(ownerProc)
                write(c6,"(i8)") currCnt
                call MsgDump (h//&
                     " PostXYGrid ("//trim(adjustl(c0))//","//trim(adjustl(c1))//")"//&
                     " globalBrams ("//trim(adjustl(c2))//","//trim(adjustl(c3))//")"//&
                     " localBrams ("//trim(adjustl(c4))//","//trim(adjustl(c5))//")"//&
                     " packed at pos "//trim(adjustl(c6)))
             end if

             if (onePostGrid%project) then
                onePostGrid%packWeight(currCnt,1:4) = &
                     onePostGrid%weight(xPost,yPost,1:4)
             end if
          end if
          unpackIndexCnt=unpackIndexCnt+1
          onePostGrid%unpackMap(unpackIndexCnt) = &
               onePostGrid%localSize(ownerProc) + &
               onePostGrid%disp(ownerProc)
       end do
    end do
    if (dumpLocal) then
       if(oneBramsGrid%mchnum==oneBramsGrid%master_num) then
          do ix = 1, onePostGrid%nLat*onePostGrid%nLon
             write(c0,"(i8)") ix
             write(c1,"(i8)") onePostGrid%unpackMap(ix)
             call MsgDump (h//" unpackMap("//trim(adjustl(c0))//")="//trim(adjustl(c1)))
          end do
       end if
    end if

    if (dumpLocal) then
       call MsgDump (h//" mapping (packXLocal, packYLocal) at this process:")
       call DumpIntegerPairs(h, &
            onePostGrid%packXLocal(1:onePostGrid%localSizeThisProc), &
            onePostGrid%packYLocal(1:onePostGrid%localSizeThisProc))
    end if

  end subroutine PackUnpackMaps






  subroutine CreatePostVerticals(oneNamelistFile, oneBramsGrid, &
       onePostGrid)
    type(NamelistFile), pointer :: oneNamelistFile
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    integer :: iz
    integer :: ierr
    character(len=8) :: c0, c1, c2
    character(len=*), parameter :: h="**(CreatePostVerticals)**"

    ! namelist options processed by this routine:
    ! zlevmax: # vertical levels to output
    ! ipresslev: encoded vertical axis unit
    !         0: BRAMS original sigma-z
    !         1: pressure
    !         2: height
    !         3: selected BRAMS levels (undocumented)
    ! case ipresslev== 1 or 2:
    !      inplevs: # vertical levels to output (same as zlevmax?)
    !      iplev(1:inplevs): vertical pressure or height values 
    ! case ipresslev== 3:
    !      inplevs: # vertical levels to output (same as zlevmax?)
    !      iplev(1:inplevs): BRAMS vertical levels *indices* to output


    onePostGrid%vertScaleCode = oneNamelistFile%ipresslev

    select case (onePostGrid%vertScaleCode)

    case (0, 3) ! original BRAMS verticals

       if (onePostGrid%vertScaleCode == 0) then

          ! case full BRAMS field (vertScaleCode == 0), 
          ! take namelist option zlevmax and limit nVert to nnzp - 1

          onePostGrid%nVert = max(1, &
               min(oneNamelistFile%zlevmax(oneBramsGrid%currGrid), oneBramsGrid%nnzp-1))

       else

          ! case user selected BRAMS verticals (vertScaleCode == 3), 
          ! take namelist option ipresslev and limit nVert to nnzp

          onePostGrid%nVert = max(1, &
               min(oneNamelistFile%inplevs, oneBramsGrid%nnzp))

       end if

       ! create fixed zMap according to user's selection

       allocate(onePostGrid%zMap(onePostGrid%nVert), stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") onePostGrid%nVert
          call fatal_error(h//" allocating zMap("//trim(adjustl(c0))//&
               ") fails")
       end if
       if (onePostGrid%vertScaleCode == 0) then

          ! case full BRAMS field (vertScaleCode == 0), 
          ! the first level has to be skipped on output

          do iz = 1, onePostGrid%nVert
             onePostGrid%zMap(iz) = iz+1
          end do
       else 

          ! case user selected BRAMS verticals (vertScaleCode == 3), 
          ! user selects BRAMS verticals *indices*;
          ! keep the selection up to nnzp levels

          do iz = 1, onePostGrid%nVert
!             onePostGrid%zMap(iz) = max(1, &
!                  min(oneNamelistFile%iplevs(iz), oneBramsGrid%nnzp))
             onePostGrid%zMap(iz) = iz
          end do
       end if

       allocate(onePostGrid%vertScaleValues(onePostGrid%nVert), stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") onePostGrid%nVert
          call fatal_error(h//" allocating vertScaleValues("//&
               trim(adjustl(c0))//") fails")
       end if

       do iz = 1, onePostGrid%nVert
          onePostGrid%vertScaleValues(iz) = oneBramsGrid%ztn(onePostGrid%zMap(iz))
       end do

       if (dumpLocal) then
          write(c0,"(i8)") onePostGrid%nVert
          call MsgDump (h//" vertical scale (sigma-z) has "//&
               trim(adjustl(c0))//" verticals; they are: ")
          call DumpFloating(h, onePostGrid%vertScaleValues)
       end if

    case (1, 2) ! pressure or height levels

       ! limit nVert to nnzp

       onePostGrid%nVert = max(1, &
            min(oneNamelistFile%inplevs, oneBramsGrid%nnzp))

       ! accept user provided vertical levels

       allocate(onePostGrid%vertScaleValues(onePostGrid%nVert), stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") onePostGrid%nVert
          call fatal_error(h//" allocating vertScaleValues("//&
               trim(adjustl(c0))//") fails")
       end if
       do iz = 1, onePostGrid%nVert
          onePostGrid%vertScaleValues(iz) = real(oneNamelistFile%iplevs(iz))
       end do

       ! fetch topo

       allocate(onePostGrid%topo(&
            oneBramsGrid%mxp,oneBramsGrid%myp), stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") oneBramsGrid%mxp
          write(c1,"(i8)") oneBramsGrid%myp
          call fatal_error(h//" allocating topo("//&
               trim(adjustl(c0))//","//&
               trim(adjustl(c1))//") fails")
       end if
       call GetVarFromMemToOutput('TOPT', oneBramsGrid%currGrid, onePostGrid%topo)

       ! reserve area for pi

       if (onePostGrid%vertScaleCode == 1) then

          ! case pressure levels, reserve area for pi, but do not take it,
          ! since pi varies with time; pi is fetch by routine "UpdateVerticals",
          ! that should be invoked prior to any call to OutputPostField

	allocate(onePostGrid%pi(&
               oneBramsGrid%mxp,oneBramsGrid%myp,oneBramsGrid%mzp), stat=ierr)


        !  allocate(onePostGrid%pi(&
        !       oneBramsGrid%mxp,oneBramsGrid%myp,onePostGrid%nVert), stat=ierr)
          
	  if (ierr /= 0) then
             write(c0,"(i8)") oneBramsGrid%mxp
             write(c1,"(i8)") oneBramsGrid%myp
             write(c2,"(i8)") onePostGrid%nVert
             call fatal_error(h//" allocating pi("//&
                  trim(adjustl(c0))//","//&
                  trim(adjustl(c1))//","//&
                  trim(adjustl(c2))//") fails")
          end if
          
       end if

       if (dumpLocal) then
          write(c0,"(i8)") onePostGrid%nVert
          if (onePostGrid%vertScaleCode == 1) then
             call MsgDump (h//" vertical scale (pressure) has "//&
                  trim(adjustl(c0))//" verticals; they are: ")
          else
             call MsgDump (h//" vertical scale (height) has "//&
                  trim(adjustl(c0))//" verticals; they are: ")
          end if
          call DumpFloating(h, onePostGrid%vertScaleValues)
       end if

    case default
       write(c0,"(i8)") onePostGrid%vertScaleCode
       call fatal_error(h//" unknown value of vertScaleCode="//trim(adjustl(c0)))
    end select
  end subroutine CreatePostVerticals






  subroutine UpdateVerticals(oneBramsGrid, onePostGrid)
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    character(len=8) :: c0
    character(len=*), parameter :: h="**(UpdateVerticals)**"

    ! if post uses pressure levels as vertical units, get current pi, 
    ! required to translate BRAMS sigma-z into pressure.
    ! Translation is performed by procedure ptransvar, invoked by 
    ! procedure OutputGradsField on 3D fields.

    if (onePostGrid%vertScaleCode == 1) then
       call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, onePostGrid%pi)
    end if
  end subroutine UpdateVerticals






  subroutine DestroyPostGrid(onePostGrid)
    type(PostGrid), pointer :: onePostGrid

    integer :: ierr
    type(fieldID), pointer :: curr, next
    character(len=*), parameter :: h="**(DestroyPostGrid)**"

    ! return if null onePostGrid

    if (.not. associated(onePostGrid)) then
       return
    end if

    ! destroy field ID list

    curr => onePostGrid%head
    do 
       if (associated(curr)) then
          next => curr%next
          deallocate(curr, stat=ierr)
          if (ierr /= 0) then
             call fatal_error(h//" deallocating field ID list")
          end if
          curr => next
       else
          exit
       end if
    end do
    onePostGrid%head => null()
    onePostGrid%tail => null()

    deallocate(onePostGrid%lon, stat=ierr)
    if (ierr /= 0) then
       call fatal_error(h//" deallocating lon")
    end if

    deallocate(onePostGrid%lat, stat=ierr)
    if (ierr /= 0) then
       call fatal_error(h//" deallocating lat")
    end if

    if (allocated(onePostGrid%xMap)) then
       deallocate(onePostGrid%xMap, stat=ierr)
       if (ierr /= 0) then
          call fatal_error(h//" deallocating xMap")
       end if
    end if

    if (allocated(onePostGrid%yMap)) then
       deallocate(onePostGrid%yMap, stat=ierr)
       if (ierr /= 0) then
          call fatal_error(h//" deallocating yMap")
       end if
    end if

    if (allocated(onePostGrid%weight)) then
       deallocate(onePostGrid%weight, stat=ierr)
       if (ierr /= 0) then
          call fatal_error(h//" deallocating weight")
       end if
    end if

    deallocate(onePostGrid%localSize, stat=ierr)
    if (ierr /= 0) then
       call fatal_error(h//" deallocating localSize")
    end if

    deallocate(onePostGrid%disp, stat=ierr)
    if (ierr /= 0) then
       call fatal_error(h//" deallocating disp")
    end if

    deallocate(onePostGrid%packXLocal, stat=ierr)
    if (ierr /= 0) then
       call fatal_error(h//" deallocating packXLocal")
    end if

    deallocate(onePostGrid%packYLocal, stat=ierr)
    if (ierr /= 0) then
       call fatal_error(h//" deallocating packYLocal")
    end if

    deallocate(onePostGrid%packWeight, stat=ierr)
    if (ierr /= 0) then
       call fatal_error(h//" deallocating packWeight")
    end if

    deallocate(onePostGrid%unpackMap, stat=ierr)
    if (ierr /= 0) then
       call fatal_error(h//" deallocating unpackMap")
    end if

    deallocate(onePostGrid%vertScaleValues, stat=ierr)
    if (ierr /= 0) then
       call fatal_error(h//" deallocating vertScaleValues")
    end if

    if (allocated(onePostGrid%zMap)) then
       deallocate(onePostGrid%zMap, stat=ierr)
       if (ierr /= 0) then
          call fatal_error(h//" deallocating zMap")
       end if
    end if

    if (allocated(onePostGrid%pi)) then
       deallocate(onePostGrid%pi, stat=ierr)
       if (ierr /= 0) then
          call fatal_error(h//" deallocating pi")
       end if
    end if

    if (allocated(onePostGrid%topo)) then
       deallocate(onePostGrid%topo, stat=ierr)
       if (ierr /= 0) then
          call fatal_error(h//" deallocating topo")
       end if
    end if

    deallocate(onePostGrid)
    onePostGrid => null()
  end subroutine DestroyPostGrid






  subroutine DumpPostGrid(onePostGrid)
    type(PostGrid), pointer :: onePostGrid

    integer :: iy
    character(len=8) :: c0, c1
    character(len=16) :: d0, d1
    character(len=*), parameter :: h="**(DumpPostGrid)**"


    if (onePostGrid%unitBinFile == INT_UNDEF) then
       call MsgDump(h//" grads binary file is not opened at this proc")
    else
       write(c0,"(i8)") onePostGrid%unitBinFile
       write(c1,"(i8)") onePostGrid%lastRec
       call MsgDump (h//&
            " grads binary file "//trim(onePostGrid%binFileName)//&
            " at unit "//trim(adjustl(c0))//&
            "; current record size "//trim(adjustl(c1)))
    end if

    if (onePostGrid%unitCtlFile == INT_UNDEF) then
       call MsgDump (h//&
            " grads ascii control file is not opened at this proc")
    else
       write(c0,"(i8)") onePostGrid%unitCtlFile
       call MsgDump (h//&
            " grads ascii control file "//trim(onePostGrid%binFileName)//&
            " at unit "//trim(adjustl(c0)))
    end if

    write(c0,"(i8)") onePostGrid%nLon
    write(d0,"(e16.7)") onePostGrid%firstLon
    write(d1,"(e16.7)") onePostGrid%delLon
    call MsgDump (h//" has "//&
         trim(adjustl(c0))//" longitudes; "//&
         " first longitude is "//trim(adjustl(d0))//&
         " delta longitude is "//trim(adjustl(d1)))
    call MsgDump (h//" all longitudes: ")
    call DumpFixed(h, onePostGrid%lon)

    write(c0,"(i8)") onePostGrid%nLat
    write(d0,"(e16.7)") onePostGrid%firstLat
    write(d1,"(e16.7)") onePostGrid%delLat
    call MsgDump (h//" has "//&
         trim(adjustl(c0))//" latitudes; "//&
         " first latitude is "//trim(adjustl(d0))//&
         " delta latitude is "//trim(adjustl(d1)))
    call MsgDump (h//" all latitudes: ")
    call DumpFixed(h, onePostGrid%lat)


    if (onePostGrid%project) then
       call MsgDump (h//" BRAMS grid is projected onto post grid")
       call MsgDump (h//" post grid is composed from the following "//&
            "pairs of BRAMS indices:")
       do iy = 1, onePostGrid%nLat
          write(c0,"(i8)") iy
          call DumpIntegerPairs(h//" post lat="//trim(adjustl(c0)), &
               onePostGrid%xMap(:,iy), onePostGrid%yMap(:,iy))
       end do
    else
       call MsgDump (h//" BRAMS and post grid have similar spacing")
       write(c0,"(i8)") onePostGrid%xStart
       write(c1,"(i8)") onePostGrid%yStart
       call MsgDump (h//&
            " post grid starts at BRAMS point ("//&
            trim(adjustl(c0))//","//trim(adjustl(c1))//")")
    end if


    call MsgDump (h//" post grid size at each process is")
    call DumpInteger(h, onePostGrid%localSize)

    call MsgDump (h//" post grid displacement at each process is")
    call DumpInteger(h, onePostGrid%disp)

    call MsgDump (h//" post grid packing map at this process is")
    call DumpIntegerPairs(h, &
         onePostGrid%packXLocal(1:onePostGrid%localSizeThisProc), &
         onePostGrid%packYLocal(1:onePostGrid%localSizeThisProc))

    call MsgDump (h//" post grid unpack map is")
    call DumpInteger(h, onePostGrid%unpackMap)

    write(c0,"(i8)") 1
    write(c1,"(i8)") onePostGrid%nVert
    call MsgDump (h//" post z axis build from  BRAMS z indices "//&
         trim(adjustl(c0))//" to "//trim(adjustl(c1)))
  end subroutine DumpPostGrid





  subroutine saveCurrentFieldID(onePostGrid, vertCode)
    type(PostGrid), pointer :: onePostGrid
    integer, intent(in) :: vertCode

    integer :: ierr
    type(fieldID), pointer :: oneFieldID
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(saveCurrentFieldID)**"

    ! return if list is done, to prevent duplicating
    ! list entries on

    if (.not. onePostGrid%buildingFieldIDList) then
       return
    end if

    allocate(oneFieldID, stat=ierr)
    if (ierr /= 0) then
       call fatal_error(h//" allocate oneFieldID fails")
    end if
    oneFieldID%ivar_type = onePostGrid%ivar_type
    oneFieldID%nLevCode  = vertCode
    oneFieldID%fieldName = onePostGrid%fieldName
    oneFieldID%fieldDescription = onePostGrid%fieldDescription
    oneFieldID%fieldUnits = onePostGrid%fieldUnits

    ! this code fragment imposes that tail and head
    ! will be either both not associated or 
    ! both associated.
    ! observe that tail and head are both initialized to null()
    ! and that tail is associated at the end of this procedure
    ! consequently, whenever tail is not associated, so is head
    ! and both will be associated at the end of this procedure.

    if (associated(onePostGrid%tail)) then
       onePostGrid%tail%next => oneFieldID
    else
       onePostGrid%head => oneFieldID
    end if
    onePostGrid%tail => oneFieldID

    if (dumpLocal) then
       write(c0,"(i8)") oneFieldID%ivar_type
       write(c1,"(i8)") oneFieldID%nLevCode
       call MsgDump (h//" node ivar_type="//trim(adjustl(c0))//&
            ", named "//trim(oneFieldID%fieldName)//&
            ", grads vert code "//trim(adjustl(c1))//&
            ", description "//trim(oneFieldID%fieldDescription)//&
            ", units "//trim(oneFieldID%fieldUnits))
    end if
  end subroutine saveCurrentFieldID


  subroutine OutputGradsField_2D (oneBramsGrid, onePostGrid, OutputField)
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid
    real, intent(in) :: OutputField(:,:)

    integer :: i
    integer :: xBrams, yBrams
    real :: localChunk(onePostGrid%localSizeThisProc)
    real :: gathered(onePostGrid%nLon*onePostGrid%nLat)
    real :: OutputArray(onePostGrid%nLon*onePostGrid%nLat)
    character(len=8) :: c0, c1, c2
character(len=16)::d0
integer :: ix, iy
    character(len=*), parameter :: h="**(OutputGradsField_2D)**"

    if (.not. associated(onePostGrid)) then
       call fatal_error(h//" invoked with null onePostGrid")
    end if

   if (dumpLocal) then
       write(*,"(a)") h//" starts for field "//trim(onePostGrid%fieldName)
    end if

    ! outputs selected part of the computed field

    select case (onePostGrid%ivar_type)
    case (2)

       ! save corresponding fieldID entry on PostGrid list

       call saveCurrentFieldID(onePostGrid, 0)

       ! surface fields dimensioned (mxp, myp);

       ! pack field entries at this process:

       call BuildLocalChunk(onePostGrid, OutputField, localChunk)

       ! master_proc gathers packed fields from all process

       if (dumpLocal) then
          write(c0,"(i8)") onePostGrid%localSizeThisProc
          call MsgDump (h//" will gather local chunk of size "//&
               trim(adjustl(c0)))
       end if
       call parf_GatherPostSfc(localChunk, &
            onePostGrid%localSize, onePostGrid%disp, &
            gathered, oneBramsGrid%master_num)

       ! master_proc dumps gathered field

       if (oneBramsGrid%mchnum == oneBramsGrid%master_num) then
          onePostGrid%lastRec = onePostGrid%lastRec + 1
          i = 0
          do iy = 1, onePostGrid%nLat
             do ix = 1, onePostGrid%nLon
                i = i+1
                OutputArray(i)=gathered(onePostGrid%unpackMap(i))
                if (dumpLocal) then
                   write(c0,"(i8)") ix
                   write(c1,"(i8)") iy
                   write(c2,"(i8)") onePostGrid%unpackMap(i)
                   write(d0,"(e16.8)") OutputArray(i)
                   call MsgDump (h//" field("//trim(adjustl(c0))//","//trim(adjustl(c1))//")="//&
                        "gathered("//trim(adjustl(c2))//")="//trim(adjustl(d0)))
                end if
             end do
          end do

          if(IPOS==2 .or. isTimeTograds()) then
            write (onePostGrid%unitBinFile, rec=onePostGrid%lastRec) OutputArray
#ifdef cdf
          elseif(IPOS==3) then
            call writeNetCdf2D(trim(onePostGrid%fieldName),onePostGrid%nLon,onePostGrid%nLat,OutputArray)
#endif
          endif

          if (dumpLocal) then
             write(c0,"(i8)") size(gathered,1)
             write(c1,"(i8)") onePostGrid%lastRec
             call MsgDump (h//" master_proc wrote field of size "//&
                  trim(adjustl(c0))//" at record "//trim(adjustl(c1))//&
                  " of grads binary file")
          end if
       end if

    case default
       write(c0,"(i8)") onePostGrid%ivar_type
       call fatal_error(h//" not prepared to handle ivar_type="//trim(adjustl(c0)))
    end select
  end subroutine OutputGradsField_2D

  subroutine OutputGradsField_3D (oneBramsGrid, onePostGrid, OutputField)
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid
    real, intent(inout) :: OutputField(:,:,:)

    integer :: k
    integer :: i
    integer :: zSize
    integer :: xBrams, yBrams
    integer :: lenName
    integer :: lenDescription
    real :: localChunk(onePostGrid%localSizeThisProc)
    real :: gathered(onePostGrid%nLon*onePostGrid%nLat)
    character(len=8) :: c0, c1, c2
    character(len=*), parameter :: h="**(OutputGradsField_3D)**"

    if (.not. associated(onePostGrid)) then
       call fatal_error(h//" invoked with null onePostGrid")
    end if
    if (dumpLocal) then
       write(*,"(a)") h//" starts for field "//trim(onePostGrid%fieldName)
    end if

    ! outputs selected part of the computed field
    !
    ! if ivar_type == 2, dumps surface of 3D field;
    ! if ivar_type == 3, dumps selected verticals
    ! if ivar_type == 7, dumps all npatches
    

    zSize=size(OutputField,3)

    select case (onePostGrid%ivar_type)
    case (2)

       ! save corresponding fieldID entry on PostGrid list

       call saveCurrentFieldID(onePostGrid, 0)

       ! surface fields dimensioned as 3D fields (mxp, myp, mzp);
       ! dumps just the surface

       ! pack surface field entries at this process
       ! (assumes surface is vertical level 2)

       call BuildLocalChunk(onePostGrid, OutputField(:,:,2), localChunk)

       ! master_proc gathers packed fields from all process

       if (dumpLocal) then
          write(c0,"(i8)") onePostGrid%localSizeThisProc
          call MsgDump (h//" will gather surface local chunk of size "//&
               trim(adjustl(c0)))
       end if
       call parf_GatherPostSfc(localChunk, &
            onePostGrid%localSize, onePostGrid%disp, &
            gathered, oneBramsGrid%master_num)

       ! master_proc dumps gathered field

       if (oneBramsGrid%mchnum == oneBramsGrid%master_num) then
          onePostGrid%lastRec = onePostGrid%lastRec + 1
            if(IPOS==2 .or. isTimeTograds()) then
              write (onePostGrid%unitBinFile, rec=onePostGrid%lastRec) &
                  gathered(onePostGrid%unpackMap(:))
#ifdef cdf
            elseif(IPOS==3) then
               call writeNetCdf2D(trim(onePostGrid%fieldName),onePostGrid%nLon &
                ,onePostGrid%nLat,gathered(onePostGrid%unpackMap(:)))
#endif
            endif          
            if (dumpLocal) then
             write(c0,"(i8)") size(gathered,1)
             write(c1,"(i8)") onePostGrid%lastRec
             call MsgDump (h//" master_proc wrote surface of field with size "//&
                  trim(adjustl(c0))//" at record "//trim(adjustl(c1))//&
                  " of grads binary file")
          end if
       end if

    case (3)

       ! save corresponding fieldID entry on PostGrid list

       call saveCurrentFieldID(onePostGrid, onePostGrid%nVert)

       ! full 3D fields (mxp, myp, mzp);
       ! change vertical scale from Brams to Post
       ! dumps desired verticals, one vertical at a time

       call ChangeVerticalScale(OutputField, onePostGrid, oneBramsGrid)

       do k = 1, onePostGrid%nVert

          ! pack vertical level field entries at this process


          call BuildLocalChunk(onePostGrid, OutputField(:,:,k), localChunk)

          ! master_proc gathers packed fields from all process

          if (dumpLocal) then
             write(c0,"(i8)") onePostGrid%localSizeThisProc
             write(c2,"(i8)") k
             call MsgDump (h//" will gather local chunk of size "//&
                  trim(adjustl(c0))//" of vertical "//trim(adjustl(c2)))
          end if
          call parf_GatherPostSfc(localChunk, &
               onePostGrid%localSize, onePostGrid%disp, &
               gathered, oneBramsGrid%master_num)

          ! master_proc dumps gathered field

          if (oneBramsGrid%mchnum == oneBramsGrid%master_num) then
             onePostGrid%lastRec = onePostGrid%lastRec + 1

            if(IPOS==2 .or. isTimeTograds()) then
              write (onePostGrid%unitBinFile, rec=onePostGrid%lastRec) &
                  gathered(onePostGrid%unpackMap(:))
#ifdef cdf
            elseif(IPOS==3) then
              call writeNetCdf3D(trim(onePostGrid%fieldName),onePostGrid%nLon &
                ,onePostGrid%nLat,k,gathered(onePostGrid%unpackMap(:)))
#endif
            endif

             if (dumpLocal) then
                write(c0,"(i8)") size(gathered,1)
                write(c1,"(i8)") onePostGrid%lastRec
                write(c2,"(i8)") k
                call MsgDump (h//" master_proc wrote field of size "//&
                     trim(adjustl(c0))//" of vertical "//trim(adjustl(c2))//&
                     " at record "//trim(adjustl(c1))//" of grads binary file")
             end if
          end if
       end do


    case (7)

       ! fields (mxp, myp, npatch);
       ! dumps all patches, updating each fieldID entry

       ! saves original current field name and description sizes

       lenName  = len_trim(onePostGrid%fieldName)
       lenDescription = len_trim(onePostGrid%fieldDescription)

       do k = 1, size(OutputField,3)

          ! insert patch number on current field name and description

          write(c0,"(i8)") k
          onePostGrid%fieldName  = onePostGrid%fieldName(1:lenName)//trim(adjustl(c0))
          onePostGrid%fieldDescription = onePostGrid%fieldDescription(1:lenDescription)//&
               ": patch # "//trim(adjustl(c0))

          ! save corresponding fieldID entry on PostGrid list

          call saveCurrentFieldID(onePostGrid, 0)

          ! pack vertical level field entries at this process

          call BuildLocalChunk(onePostGrid, OutputField(:,:,k), localChunk)

          ! master_proc gathers packed fields from all process

          if (dumpLocal) then
             write(c0,"(i8)") onePostGrid%localSizeThisProc
             write(c2,"(i8)") k
             call MsgDump (h//" will gather local chunk of size "//&
                  trim(adjustl(c0))//" of patch "//trim(adjustl(c2)))
          end if
          call parf_GatherPostSfc(localChunk, &
               onePostGrid%localSize, onePostGrid%disp, &
               gathered, oneBramsGrid%master_num)

          ! master_proc dumps gathered field

          if (oneBramsGrid%mchnum == oneBramsGrid%master_num) then
             onePostGrid%lastRec = onePostGrid%lastRec + 1

            if(IPOS==2) then
                write (onePostGrid%unitBinFile, rec=onePostGrid%lastRec) &
                  gathered(onePostGrid%unpackMap(:))
#ifdef cdf
            elseif(IPOS==3) then
              call writeNetCdf2D(trim(onePostGrid%fieldName),onePostGrid%nLon &
                ,onePostGrid%nLat,gathered(onePostGrid%unpackMap(:)))
#endif
            endif

             if (dumpLocal) then
                write(c0,"(i8)") size(gathered,1)
                write(c1,"(i8)") onePostGrid%lastRec
                write(c2,"(i8)") k
                call MsgDump (h//" master_proc wrote field of size "//&
                     trim(adjustl(c0))//" of patch "//trim(adjustl(c2))//&
                     " at record "//trim(adjustl(c1))//" of grads binary file")
             end if
          end if
       end do

    case default
       write(c0,"(i8)") onePostGrid%ivar_type
       call fatal_error(h//" not prepared to handle ivar_type="//trim(adjustl(c0)))
    end select
  end subroutine OutputGradsField_3D

  subroutine OutputGradsField_4D (oneBramsGrid, onePostGrid, OutputField)
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid
    real, intent(in) :: OutputField(:,:,:,:)

    integer :: k
    integer :: l
    integer :: i
    integer :: lenName
    integer :: lenDescription
    integer :: vSize
    real :: localChunk(onePostGrid%localSizeThisProc)
    real :: gathered(onePostGrid%nLon*onePostGrid%nLat)
    character(len=8) :: c0, c1, c2, c3
    character(len=*), parameter :: h="**(OutputGradsField_4D)**"

    if (.not. associated(onePostGrid)) then
       call fatal_error(h//" onePostGrid not associated")
    end if
    if (dumpLocal) then
       write(*,"(a)") h//" starts for field "//trim(onePostGrid%fieldName)
    end if
    vSize=size(OutputField,3)
    ! outputs selected part of the computed field:
    ! all verticals of each patch

    select case (onePostGrid%ivar_type)
    case (8)

       ! fields (mxp, myp, nzg, npatch);
       ! dumps all verticals of each patch, updating each fieldID entry

       ! saves original current field name and description sizes

       lenName  = len_trim(onePostGrid%fieldName)
       lenDescription = len_trim(onePostGrid%fieldDescription)


       do l = 1, size(OutputField,4)

          ! insert patch number on current field name and description

          write(c0,"(i8)") l
          onePostGrid%fieldName  = onePostGrid%fieldName(1:lenName)//trim(adjustl(c0))
          onePostGrid%fieldDescription = onePostGrid%fieldDescription(1:lenDescription)//&
               ": patch # "//trim(adjustl(c0))

          ! save corresponding fieldID entry on PostGrid list

          call saveCurrentFieldID(onePostGrid, oneBramsGrid%nzg)

          do k = 1, size(OutputField,3)

             ! pack vertical level field entries at this process

             call BuildLocalChunk(onePostGrid, OutputField(:,:,k,l), localChunk)

             ! master_proc gathers packed fields from all process

             if (dumpLocal) then
                write(c0,"(i8)") onePostGrid%localSizeThisProc
                write(c2,"(i8)") k
                write(c3,"(i8)") l
                call MsgDump (h//" will gather local chunk of size "//&
                     trim(adjustl(c0))//" of vertical "//trim(adjustl(c2))//" and patch "//&
                     trim(adjustl(c3)))
             end if
             call parf_GatherPostSfc(localChunk, &
                  onePostGrid%localSize, onePostGrid%disp, &
                  gathered, oneBramsGrid%master_num)

             ! master_proc dumps gathered field

             if (oneBramsGrid%mchnum == oneBramsGrid%master_num) then
                onePostGrid%lastRec = onePostGrid%lastRec + 1

              if(IPOS==2 .or. isTimeTograds()) then
                write (onePostGrid%unitBinFile, rec=onePostGrid%lastRec) &
                     gathered(onePostGrid%unpackMap(:))
#ifdef cdf
              elseif(ipos==3) then
                call writeNetCdf3D(trim(onePostGrid%fieldName),onePostGrid%nLon &
                ,onePostGrid%nLat,k,gathered(onePostGrid%unpackMap(:)))
#endif
              endif


                if (dumpLocal) then
                   write(c0,"(i8)") size(gathered,1)
                   write(c1,"(i8)") onePostGrid%lastRec
                   write(c2,"(i8)") k
                   write(c3,"(i8)") l
                   call MsgDump (h//" master_proc wrote field of size "//&
                        trim(adjustl(c0))//" of vertical "//trim(adjustl(c2))//&
                        " and patch "//trim(adjustl(c3))//&
                        " at record "//trim(adjustl(c1))//" of grads binary file")
                end if
             end if

          end do
       end do

    case default
       write(c0,"(i8)") onePostGrid%ivar_type
       call fatal_error(h//" not prepared to handle ivar_type="//trim(adjustl(c0)))
    end select
  end subroutine OutputGradsField_4D




  subroutine ChangeVerticalScale(OutputField, onePostGrid, oneBramsGrid)
    real, intent(inout) :: OutputField(:,:,:)
    type(PostGrid), pointer :: onePostGrid
    type(BramsGrid), pointer :: oneBramsGrid

    integer :: iz
    character(len=8) :: c0
    character(len=*), parameter :: h="**(ChangeVerticalScale)**"

    select case (onePostGrid%vertScaleCode)

    case (0, 3) ! sigma-z

       do iz = 1, onePostGrid%nVert
          OutputField(:,:,iz) = OutputField(:,:,onePostGrid%zMap(iz))
       end do

    case (1) ! pressure levels

       call ptransvar(OutputField, onePostGrid%vertScaleValues, &
            onePostGrid%pi, oneBramsGrid%ztn, onePostGrid%topo, undef)

    case (2) ! height

       call ctransvar(OutputField, onePostGrid%topo, &
            onePostGrid%vertScaleValues, oneBramsGrid%ztn, &
            oneBramsGrid%ztop, undef)

    case default
       write(c0,"(i8)") onePostGrid%vertScaleCode
       call fatal_error(h//" unknown value of vertScaleCode="//trim(adjustl(c0)))
    end select
  end subroutine ChangeVerticalScale




  subroutine BuildLocalChunk(onePostGrid, OutputField, localChunk)
    type(postGrid), pointer :: onePostGrid
    real, intent(in) :: outputField(:,:)
    real, intent(out) :: localChunk(:)

    integer :: i
    integer :: xBrams
    integer :: yBrams
    character(len=8) :: c0, c1, c2
    character(len=16) :: d0, d1, d2, d3, d4
    character(len=4) :: e1, e2, e3, e4
    character(len=*), parameter :: h="**(BuildLocalChunk)**"
    logical :: localChunkHasValue

    if (onePostGrid%project) then

       ! case projection, do local interpolation and pack the
       ! interpolated array

       do i= 1, onePostGrid%localSizeThisProc
       
       
             ! output linear combination of indices
             ! local indices (at this proc)

             xBrams = onePostGrid%packXLocal(i)
             yBrams = onePostGrid%packYLocal(i)

          if ( any(onePostGrid%packWeight(i,:) == undef) .or. &
	       OutputField(xBrams, yBrams) == undef      .or. &
	       OutputField(xBrams+1, yBrams) == undef    .or. &
	       OutputField(xBrams, yBrams+1) == undef    .or. &
	       OutputField(xBrams+1, yBrams+1) == undef              ) then

             ! output undef

             localChunk(i) = undef
             if (dumpLocal) then
                write(c0,"(i8)") i
                write(d0,"(e16.8)") localChunk(i)
                call MsgDump (h//" localChunk("//trim(adjustl(c0))//&
                     ")="//trim(adjustl(d0)))
             end if

             if (dumpLocal) then
                write(c0,"(i8)") i
                write(d0,"(e16.8)") localChunk(i)
                call MsgDump (h//" localChunk("//trim(adjustl(c0))//&
                     ")="//trim(adjustl(d0)))
             end if

          else


             ! interpolate and pack

   	     localChunk(i) = &
	     onePostGrid%packWeight(i,1) * OutputField(xBrams  , yBrams  ) + &
	     onePostGrid%packWeight(i,2) * OutputField(xBrams+1, yBrams  ) + &
	     onePostGrid%packWeight(i,3) * OutputField(xBrams  , yBrams+1) + &
	     onePostGrid%packWeight(i,4) * OutputField(xBrams+1, yBrams+1)

             if (dumpLocal) then
                write(c0,"(i8)") i
                write(c1,"(i8)") xBrams
                write(c2,"(i8)") yBrams
                write(d0,"(e16.8)") localChunk(i)
                call MsgDump (h//" localChunk("//trim(adjustl(c0))//&
                     ") uses Brams grid indices of the rectangle"//&
                     " with lower left corner ("//trim(adjustl(c1))//&
                     ","//trim(adjustl(c2))//") and has value "//trim(adjustl(d0)))
             end if
             if (dumpLocal) then
                write(c0,"(i8)") i
                write(c1,"(i8)") xBrams
                write(c2,"(i8)") yBrams
                call MsgDump (h//" localChunk("//trim(adjustl(c0))//&
                     ") uses Brams grid indices of the rectangle"//&
                     " with lower left corner ("//trim(adjustl(c1))//&
                     ","//trim(adjustl(c2))//")")
                write(d0,"(f16.8)") localChunk(i)
                write(e1,"(f4.2)") onePostGrid%packWeight(i,1)
                write(e2,"(f4.2)") onePostGrid%packWeight(i,2)
                write(e3,"(f4.2)") onePostGrid%packWeight(i,3) 
                write(e4,"(f4.2)") onePostGrid%packWeight(i,4)
                write(d1,"(f16.8)") OutputField(xBrams  , yBrams  )
                write(d2,"(f16.8)") OutputField(xBrams+1, yBrams  )
                write(d3,"(f16.8)") OutputField(xBrams  , yBrams+1)
                write(d4,"(f16.8)") OutputField(xBrams+1, yBrams+1)
                call MsgDump (h//" localChunk("//trim(adjustl(c0))//")="//&
                     trim(adjustl(d0))//"="//&
                     trim(adjustl(e1))//"*"//trim(adjustl(d1))//" + "//&
                     trim(adjustl(e2))//"*"//trim(adjustl(d2))//" + "//&
                     trim(adjustl(e3))//"*"//trim(adjustl(d3))//" + "//&
                     trim(adjustl(e4))//"*"//trim(adjustl(d4)))
             end if
          end if
       end do

    else
	!print*, 'sem projecao'
       ! case no projection, just pack the array

       do i= 1, onePostGrid%localSizeThisProc
          xBrams = onePostGrid%packXLocal(i)
          yBrams = onePostGrid%packYLocal(i)
          localChunk(i) = OutputField(xBrams,yBrams)
       end do
    end if
  end subroutine BuildLocalChunk

end module ModPostGrid
