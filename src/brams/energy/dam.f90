module  dam
  !# modulo para barragens (hidroeletricas)
  !#
  !# @note
  !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
  !#
  !# **Brief**: modulo para barragens (hidroeletricas)
  !#
  !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
  !#
  !# **Author**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
  !#
  !# **Date**: Abr/2019
  !# @endnote
  !#
  !# @changes
  !#
  !# +
  !# @endchanges
  !# @bug
  !#
  !#@endbug
  !#
  !#@todo
  !#  &#9744; <br/>
  !# @endtodo
  !#
  !# @warning
  !# Now is under CC-GPL License, please see
  !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
  !# @endwarning
  !#
  !#---
  use ModNamelistFile, only: namelistFile

  implicit none

  integer :: damModule
  !# Variable to control module run
  real :: frqPrecip
  !# Frequency (seconds) to output dams precipitation
  character(len=256) :: damOutPrefix
  !# Name of outputfile prefix
  character(len=256) :: damTableName
  !# Name of dam table file
  character(len=256) :: damTablePath
  !# Name of dam directory

  logical :: firstTime

  include "constants.f90"
  character(len=*), parameter :: header="**(dam)**"

  !Parameters(constants)

  !module variables
  integer :: nDams
  !# number of dams
  character(len=256) :: damFilesPath
  !# path of damFiles

  type dm
    character(len=256) :: damFileName
    !# Name of dam file
    character(len=256) :: damName
    !# dam name
    integer :: numberOfPoints
    real(kind=kind_rb), pointer :: latsBorder(:)
    real(kind=kind_rb), pointer :: lonsBorder(:)
    logical, pointer :: inside(:,:)
    real :: precipitation
  end type dm
  type(dm), allocatable :: dams(:)
  !# dam information

  contains

  real(kind=kind_rb) function getDistance(LatA,LngA,LatB,LngB)
  !# Get distance between 2 points in lat,lon [m]
  !#
  !# @note
  !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
  !#
  !# **Brief**: using a equation to determine distance between 2 points on
  !# earth surface
  !#
  !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
  !#
  !# **Author**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
  !#
  !# **Date**: 28Mar2019
  !# @endnote
  !#
  !# @changes
  !#
  !# +
  !# @endchanges
  !# @bug
  !#
  !#@endbug
  !#
  !#@todo
  !#  &#9744; <br/>
  !# @endtodo
  !#
  !# @warning
  !# Now is under CC-GPL License, please see
  !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
  !# @endwarning
  !#
  !#---
  use dump, only: &
    dumpMessage

  include "constants.f90"
  character(len=*), parameter :: header="**(getDistance)**"

  !Input/output variables
  real(kind=kind_rb),intent(in) :: latA
  !# Latitude of 1st POINT
  real(kind=kind_rb),intent(in) :: lngA
  !# Longitude of 1st POINT
  real(kind=kind_rb),intent(in) :: latB
   !# Latitude of 2nd POINT
  real(kind=kind_rb),intent(in) :: lngB
  !# Longitude of 2nd POINT

  !Local Variables
  real(kind=kind_rb) :: a,b,c,d,e
  real(kind=kind_rb) :: cosa,cosb,cosd,cose,sinc,sinb,sind
  real(kind=kind_rb) :: part,acospart

  !Code
  a=c_pi*(90.-LatB)/180.
  b=(90.-LatA)*c_pi180
  c=(90.-LatB)*c_pi180
  d=(90.-LatA)*c_pi180
  e=( LngA - LngB)*c_pi180
  cosa=cos(a)
  cosb=cos(b)
  cosd=cos(d)
  sinc=sin(c)
  sind=sin(d)
  sinb=sin(b)
  cose=cos(e)
  part=cosa*cosb+sinc*sind*cose
  !if(part<-1.-1E-10 .or. part>1.0+1E-10) then
  !  write (*,fmt='("Info: ",5(E30.24,1X))') part,LatA,LatB,LngA,LngB
  !  iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion &
  !      ,c_fatal,'Something wrong with Lat and Lon data. Please, check it')
  !endif
  if(part<-1.0) part=-1.0
  if(part>1.0) part=1.0
  acospart=acos(part)
  getDistance = c_erad*acospart

end function getDistance

logical function isInside(latP,lonP,numberOfPoints,latsBorder,lonsBorder)
  use dump, only: &
    dumpMessage
  !# verify if a point (lat,lon) is inside a border points (.true. or .false.)
  !#
  !# @note
  !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
  !#
  !# **Brief**: verify if a point (lat,lon) is inside a border points
  !#
  !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
  !#
  !# **Author**: Luiz Flavio Rodrigues **&#9993;**<<luiz.rodrigues@inpe.br>>
  !#
  !# **Date**: 28Mar2019
  !# @endnote
  !#
  !# @changes
  !#
  !# +
  !# @endchanges
  !# @bug
  !#
  !#@endbug
  !#
  !#@todo
  !#  &#9744; <br/>
  !# @endtodo
  !#
  !# @warning
  !# Now is under CC-GPL License, please see
  !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
  !# @endwarning
  !#
  !#---
  include "constants.f90"
  character(len=*), parameter :: header="**(isInside)**"

  !Parameters(constants)


  !Input/output variables
  integer, intent(in) :: numberOfPoints
  !# Total of frontiers points on border
  real(kind=kind_rb),intent(in) :: latP
  !# latitude of point to be tested
  real(kind=kind_rb),intent(in) :: lonP
  !# longitude of point to be tested
  real(kind=kind_rb),intent(in) :: latsBorder(numberOfPoints)
  !# array with lats of border
  real(kind=kind_rb),intent(in) :: lonsBorder(numberOfPoints)
  !# array with lons of border

  !Local variables
  real(kind=kind_rb) :: degrees
  !# sum degrees over all points
  real(kind=kind_rb) :: aSide,bSide,cSide
  !# sides of triangle from borders
  real(kind=kind_rb) :: auxT
  real(kind=kind_rb) :: angRad
  real(kind=kind_rb) :: angGra
  real(kind=kind_rb) :: ta_x,ta_y,tb_x,tb_y
  real(kind=kind_rb) :: cross
  logical :: clockwise

  integer :: bc
  !# border count
  integer :: bcp1
  !# border count plus one

  !Code
  degrees = 0.
  isInside=.false.

  do bc=1,numberOfPoints-1
    bcp1=bc+1
    !calculate distance of points
    aSide = getDistance(latsBorder(bc),lonsBorder(bc),latsBorder(bcp1) &
            ,lonsBorder(bcp1))
    if(aside<tinyReal) cycle
    bSide = getDistance(latP,lonP,latsBorder(bc),lonsBorder(bc))
    cSide = getDistance(latP,lonP,latsBorder(bcp1),lonsBorder(bcp1))
    if(bside<tinyReal .or. cside<tinyReal) cycle
    !Determining the angle between 2 points and reference point
    auxT=(bSide*bSide+cSide*cSide-aSide*aSide)/(2.0*bSide*cSide)
    angRad=0.0
    if(auxT<1.0) angRad=acos(auxT)
    angGra=c_i_pi180*angRad
    !write (*,fmt='(I5.5,1X,7(E15.5,1X))') bc,aSide,bSide,cSide,auxT,angRad &
    !                                      ,angGra,degrees
    !calculate direction of vector
    ta_x = latsBorder(bc) - latP
    ta_y = lonsBorder(bc) - lonP
    tb_x = latsBorder(bcp1) - latP
    tb_y = lonsBorder(bcp1) - lonP
    cross = tb_y * ta_x - tb_x * ta_y
    !use to test if cw or ccw
    clockwise=.false.
    if(cross<0.0) clockwise=.true.
    if(angGra/=0.) then
      if(clockwise) then
        degrees = degrees + angGra
      else
        degrees = degrees - angGra
      endif
    endif
    !print *,'Deg: ', degrees
  enddo
  if(abs(abs(nint(degrees)) - 360) <= 3) isInside=.true.

end function isInside

subroutine readDamTable(filePath,fileName)
    !# read a file with dams files information
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: read a file with a list of files with dams information
    !#           The file must have the following data:
    !#  nDams !1st line with number of dams files
    !#  folder !2nd line with folder that file are
    !#  file1 !3rd line and (nDams) lines with fileName and dam name
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 03Apr2019
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !#
    !#@endbug
    !#
    !#@todo
    !#  &#9744; <br/>
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !#
    !# @endwarning
    !#
    !#---

  use dump, only: &
    dumpMessage

  implicit none

  include "constants.f90"
  character(len=*), parameter :: header="**(readDamTable)**"

  !Parameters(constants)
  integer, parameter :: funit=33

  !Input/output variables
  character(len=*),intent(in) :: fileName
  !# Name of dam table file
  character(len=*),intent(in) :: filePath
  !# dam file directory

  !Local variables
  logical :: fileExist
  integer :: nd

  !Code
  inquire(file=trim(filePath)//'/'//trim(fileName), exist=fileExist )
  if ( .not. fileExist ) iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion &
      ,c_fatal,'file '//trim(filePath)//'/'//trim(fileName) &
      //'dams table not found. Please, check it!')

  open(unit=funit,file=trim(filePath)//'/'//trim(fileName) &
       ,action='READ',status='OLD',form='FORMATTED')

  write(*,fmt='(A)') '*** Reading the list of dams from '//trim(filePath)//'/'//trim(fileName)
  read(funit,*) nDams
  write(*,fmt='(A,I3.3)') '*** #ofDams=',nDams
  read(funit,fmt='(A)') damFilesPath
  write(*,fmt='(A)') '*** Path of dams files is '//trim(damFilesPath)

  if(.not. allocated(dams)) allocate(dams(nDams))

  do nd=1,nDams
    read(funit,*) dams(nd)%damFileName,dams(nd)%damName
    write(*,fmt='(A,I3.3,A)') '*** Filename for ',nd,' is '//trim(dams(nd)%damFileName) &
         //' nickname: '//trim(dams(nd)%damName)
  enddo

end subroutine readDamTable

subroutine allocDams(nx,ny)
  use dump, only: &
    dumpMessage

  implicit none

  include "constants.f90"
  character(len=*), parameter :: header="**(allocDams)**"

  !Parameters(constants)

  !Input/output variables
  integer, intent(in) :: nx
  integer, intent(in) :: ny

  !Local variables
  integer :: nd

  !Code
  if(.not. allocated(dams)) iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion &
      ,c_fatal,'dams table not read yet. Please, read it first!')

  do nd=1,nDams
    allocate(dams(nd)%inside(nx,ny))
  enddo

end subroutine allocDams


subroutine readBln()
    !# Read the contour file of dams
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: read each contour file of dams named in dam table
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 03Apr2019
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !#
    !#@endbug
    !#
    !#@todo
    !#  &#9744; <br/>
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !#
    !# @endwarning
    !#
    !#---

  use dump, only: &
    dumpMessage

  include "constants.f90"
  character(len=*), parameter :: header="**(readBln)**"

  !Parameters(constants)
  integer, parameter :: funit=33


  !Input/output variables

  !# name of file with contours

  !Local variables
  logical :: fileExist
  integer :: ierr
  integer :: icop
  integer :: dummy
  character(len=64) :: lineHeader
  character(len=2) :: cnumSIze
  real(kind=kind_rb) :: lon,lat
  integer :: nd

  !Code
  do nd=1,nDams
    !write(*,fmt='(A)') '*** Reading '//trim(damFilesPath)//'/'//trim(dams(nd)%damFileName)
    inquire(file=trim(damFilesPath)//'/'//trim(dams(nd)%damFileName), exist=fileExist )
    if (.not. fileExist) iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion &
    ,c_fatal,'file not found. Please, check it!')

    open(unit=funit,file=trim(damFilesPath)//'/'//trim(dams(nd)%damFileName) &
         ,action='READ',status='OLD',form='FORMATTED')

    !Reading the header line
    read(funit,*) dams(nd)%numberOfPoints,dummy

    allocate(dams(nd)%latsBorder(dams(nd)%numberOfPoints))
    allocate(dams(nd)%lonsBorder(dams(nd)%numberOfPoints))

    do icop=1,dams(nd)%numberOfPoints
      read(funit,*) dams(nd)%lonsBorder(icop),dams(nd)%latsBorder(icop)
      !write (*,fmt='(I5,1X,2(F18.14,1X))') icop,lonsBorder(icop),latsBorder(icop)
    enddo

    close(funit)
  enddo

end subroutine  readBln

subroutine fillValidLatLons(nx,ny,glat,glon)
    !# Fill the inside (x,y) array with true if lat,lon is inside a dam
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: If a pojnt lat,lon from model is inside a dam area the inside is true.
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 03Apr2019
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !#
    !#@endbug
    !#
    !#@todo
    !#  &#9744; <br/>
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !#
    !# @endwarning
    !#
    !#---

  use dump, only: &
    dumpMessage

  implicit none

  include "constants.f90"
  character(len=*), parameter :: header="**(fillValidLatLons)**"

  !Parameters(constants)

  !Input/output variables
  integer, intent(in) :: nx
  !# points in lon
  integer, intent(in) :: ny
  !# points in lat
  real, intent(in) :: glat(nx,ny)
  !# lat coordinates of each point
  real, intent(in) :: glon(nx,ny)
  !# lon coordinates of each point

  !Local variables
  integer :: nd
  integer :: i
  integer :: j
  real(kind=kind_rb) :: lon(nx,ny)
  real(kind=kind_rb) :: lat(nx,ny)

  !print *,'glon=',glon
  lon=glon
  lat=glat

  !Code
  do nd=1,nDams
    do i=1,nx
      do j=1,ny
        dams(nd)%inside(i,j)=isInside(lat(i,j),lon(i,j),dams(nd)%numberOfPoints &
                 ,dams(nd)%latsBorder,dams(nd)%lonsBorder)
      enddo
    enddo
  enddo

end subroutine fillValidLatLons


subroutine initDams(nx,ny,glat,glon,mchnum,master_num)
    !# Initializa dam module
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: read damTable, allocate all types, read each dam file and fill all
    !#  data. Send all information for processors.
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 03Apr2019
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !#
    !#@endbug
    !#
    !#@todo
    !#  &#9744; <br/>
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !#
    !# @endwarning
    !#
    !#---

  use dump, only: &
    dumpMessage

  implicit none

  include "constants.f90"
  include "mpif.h"

  character(len=*), parameter :: header="**(initDams)**"

  !Parameters(constants)

  !Input/output variables
  integer, intent(in) :: nx
  !# points in lon direction
  integer, intent(in) :: ny
  !#points in lat direction
  integer, intent(in) :: mchnum
  !# this machine number
  integer, intent(in) :: master_num
  !# the master number
  real, intent(in) :: glat(nx,ny)
  !# latitudes of each x,y point
  real, intent(in) :: glon(nx,ny)
  !# longitude for each x,y point

  !Local variables
  integer :: ierr
  integer :: nd
  integer :: np

  ! character(len=256) :: fn
  ! integer :: i,j
  ! write(fn,fmt='("File",I4.4,".out")') mchnum
  ! open(33,file=fn,status='replace')
  ! do i=1,nx
  !   do j=1,ny
  !     write(33,fmt='(2(I3.3,1X),2(F8.2,1X))') i,j,glat(i,j),glon(i,j)
  !   enddo
  ! enddo
  ! close(33)


  !Code
  firstTime=.true.
  if(mchnum==master_num) call readDamTable('./tables/dam/','damTable.dat')
  !
  CALL MPI_BCAST(nDams,1, MPI_INTEGER, master_num, MPI_COMM_WORLD, ierr)
  !
  if(mchnum/=master_num) allocate(dams(nDams))
  call allocDams(nx,ny)
  do nd=1,nDams
   CALL MPI_BCAST(dams(nd)%damName,256, MPI_CHARACTER, master_num, MPI_COMM_WORLD, ierr)
  enddo
  if(mchnum==master_num) call readBln()
  do nd=1,nDams
    CALL MPI_BCAST(dams(nd)%numberOfPoints,1, MPI_INTEGER, master_num, MPI_COMM_WORLD, ierr)
    if(mchnum/=master_num) then
      allocate(dams(nd)%latsBorder(dams(nd)%numberOfPoints))
      allocate(dams(nd)%lonsBorder(dams(nd)%numberOfPoints))
    endif
    do np=1,dams(nd)%numberOfPoints
      CALL MPI_BCAST(dams(nd)%lonsBorder(np),1, MPI_DOUBLE_PRECISION, master_num, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(dams(nd)%latsBorder(np),1, MPI_DOUBLE_PRECISION, master_num, MPI_COMM_WORLD, ierr)
    enddo
    dams(nd)%precipitation=0.0
  enddo
  if(mchnum==master_num) write(*,fmt='(A)') "*** Filling points inside dam's area"
  call fillValidLatLons(nx,ny,glat,glon)

  !For test: Uncomment bellow only for 1 processor and Paraiba do Sul river.
  !call gradsWrite()

end subroutine initDams

subroutine gradsWrite()
    !# write a grads file and store information
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: reads a grads file and store information in this module
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 01Apr2019
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !#
    !#@endbug
    !#
    !#@todo
    !#  &#9744; Changes to work with another type of grads files<br/>
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !# See the message in todo
    !# @endwarning
    !#
    !#---

  use dump, only: &
    dumpMessage

  include "constants.f90"
  character(len=*), parameter :: header="**(gradsWrite)**"

  !Parameters(constants)
  integer, parameter :: funit=33
  !# unit to manipulate file

  !Local variables
  character(len=255) :: dset
  !# first line of grads header with the name of binary file
  character(len=255) :: graName
  !# Name of binary from dset read
  character(len=255) :: title
  !# title of experiment
  character(len=256) :: dummy1
  !# just a dummy
  character(len=32) :: name
  !#
  integer :: lev
  !# number of levels for each var

  logical :: fileExist
  real(kind=kind_rb) :: undef
  real(kind=kind_rb) :: lonI, lonStep
  real(kind=kind_rb) :: latI, latStep
  integer :: nv
  integer :: time
  integer :: i,j,k,irec,r,recordLen
  real,allocatable :: vvar(:,:)
  real(kind=kind_rb) :: lon,lat

  !Code
  open(unit=funit,file='damTest.ctl' &
       ,action='WRITE',status='replace',form='FORMATTED')

  !writing the name of grads file
  write(funit,*) 'dset ^damTest.gra'
  !writing others infos to ctl
  write(funit,*) 'undef -0.9990000E+34'
  write(funit,*) 'title Rio Paraiba'
  write(funit,*) 'xdef   49 linear     -47.5242233      0.0978927'
  write(funit,*) 'ydef   29 linear     -24.4148064      0.0899690'
  write(funit,*) 'zdef 23 levels 1000.0     975.0     950.0     925.0     900.0     875.0     850.0 '
  write(funit,*) '    825.0     800.0     750.0     700.0     650.0     600.0     550.0     500.0'
  write(funit,*) '   450.0     400.0     350.0     300.0     250.0     200.0     150.0     100.0 '
  write(funit,*) 'tdef 1 linear 00:00z01jan2018     1rh'
  write(funit,*) 'vars 1'
  write(funit,*) 'INSIDE   1 99 Is Inside'
  write(funit,*) 'endvars'

  close(funit)

  recordLen=4*49*29
  allocate(vvar(49,29))
  vvar=0.0
  do i=1,49
    do j=1,29
      if(dams(1)%inside(i,j)) vvar(i,j)=1.0
    enddo
  enddo

  open(unit=funit,file='damTest.gra',&
      action='WRITE',status='REPLACE',form='UNFORMATTED',access='DIRECT', &
      recl=recordLen)

  !# writing grads binary and fill variables
  irec=1
  do time=1,1
    do nv=1,1
      do k=1,1
        write (funit,rec=irec) vvar
        irec=irec+1
      enddo
    enddo
  enddo

  close(funit)

end subroutine  gradsWrite

subroutine StoreNamelistFileAtDams(oneNamelistFile)
  type(namelistFile), pointer :: oneNamelistFile

  damModule = oneNamelistFile%damModule
  frqPrecip = oneNamelistFile%frqPrecip
  damOutPrefix = oneNamelistFile%damOutPrefix

end subroutine StoreNamelistFileAtDams

subroutine acumPrecipInDam(nx,ny,ia,iz,ja,jz,mcphys_type,aconpr,accpr &
                          ,accpp,accps,accpa,accpg,accph)
  !# Sum the total precipitation in each point inside dams
  !#
  !# @note
  !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
  !#
  !# **Brief**: Sum the cumulus precipitation and microphysics precipitation from model
  !# in each dam
  !#
  !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
  !#
  !# **Author**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
  !#
  !# **Date**: 03Apr2019
  !# @endnote
  !#
  !# @changes
  !#
  !# +
  !# @endchanges
  !# @bug
  !#
  !#@endbug
  !#
  !#@todo
  !#  &#9744; <br/>
  !# @endtodo
  !#
  !# @warning
  !# Now is under CC-GPL License, please see
  !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
  !#
  !# @endwarning
  !#
  !#---

  use dump, only: &
    dumpMessage

  implicit none

  include "constants.f90"
  character(len=*), parameter :: header="**(acumPrecipInDam)**"

  !Parameters(constants)

  !Input/output variables
  integer, intent(in) :: nx
  !# points in lon direction
  integer, intent(in) :: ny
  !#points in lat direction
  integer, intent(in) :: ia
  !# first valid point in x
  integer, intent(in) :: iz
  !# last valid point in x
  integer, intent(in) :: ja
  !# first valid point in y
  integer, intent(in) :: jz
  !# last valid point in y
  integer,intent(in) :: mcphys_type
  !# microphysics
  real, intent(in) :: aconpr(nx,ny)
  !# total of precipitation from cuparm
  real, intent(in) :: accpr(nx,ny)
  !# total of precipitation from micphys
  real, intent(in) :: accpp(nx,ny)
  !# part of precipitation frm miphys <=1
  real, intent(in) :: accps(nx,ny)
  !# part of precipitation frm miphys <=1
  real, intent(in) :: accpa(nx,ny)
  !# part of precipitation frm miphys <=1
  real, intent(in) :: accpg(nx,ny)
  !# part of precipitation frm miphys <=1
  real, intent(in) :: accph(nx,ny)
  !# part of precipitation frm miphys <=1

  !local
  integer :: nd
  integer :: i
  integer :: j

  !Code
  do nd=1,nDams
    do i=ia,iz
      do j=ja,jz
        if(dams(nd)%inside(i,j)) then
          dams(nd)%precipitation=dams(nd)%precipitation+aconpr(i,j)+accpr(i,j)
          if(mcphys_type .le. 1) then
            dams(nd)%precipitation=dams(nd)%precipitation &
            +accpp(i,j)+accps(i,j)+accpa(i,j)+accpg(i,j)+accph(i,j)
          endif
        endif
      enddo
    enddo
  enddo

end subroutine acumPrecipInDam

subroutine outputDamPrecip(time,dtlongn,timmax,mchnum,master_num)
    !# write the file with the results of precipitation
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: write the file with the results of precipitation in files.
    !#  Each file has a name of dam plus date and hour initial
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 03Apr2019
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !#
    !#@endbug
    !#
    !#@todo
    !#  &#9744; <br/>
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !#
    !# @endwarning
    !#
    !#---

  use dump, only: &
    dumpMessage

  use mem_grid, only: &
     idate1, &
     imonth1, &
     iyear1, &
     ihour1, &
     itime1

  use ModDateUtils, only: &
       date_add_to

  implicit none

  include "constants.f90"
  include "mpif.h"

  character(len=*), parameter :: header="**(outputDamPrecip)**"

  !Parameters(constants)
  integer, parameter :: funit=33

  !Input/output variables
  real, intent(in) :: Time
  !# current model time
  real, intent(in) :: dtlongn
  !# time quanta
  real, intent(in) :: timmax
  !# time of end model
  integer, intent(in) :: mchnum
  !# this machine number
  integer, intent(in) :: master_num
  !# the master number

  !Local variables
  integer :: nd
  character(len=3) :: cnd
  character(len=256) :: fmt1
  character(len=256) :: fmt2
  integer :: ierr
  real, allocatable :: totalOfPrecipitation(:)
  real :: localPrecip
  real :: masterPrecip
  character(len=256) :: outFileName
  integer :: oyr,omn,ody,otm

  !Code

  ! Do only if is time
  if(.not. (mod(time,frqPrecip)<dtlongn .or. time>=timmax)) return
  allocate(totalOfPrecipitation(nDams))
  totalOfPrecipitation=0.0

  write(cnd,fmt='(I3.3)') nDams-1
  if(nDams>1) then
    fmt1='("date,hour,",'//cnd//'(A,","),A)'
    fmt2='(I4,"/",I2.2,"/",I2.2,",",I6.6,",",'//cnd//'(F15.1,","),F15.1)'
  else
    fmt1='("date,hour",",",A)'
    fmt2='(I4,"/",I2.2,"/",I2.2,",",I6.6,",",F15.1)'
  endif
  !Integrate the precipitation
  do nd=1,nDams
    localPrecip=dams(nd)%precipitation
    !reduce the sum in master
    call MPI_Reduce(localPrecip,masterPrecip,1 &
    ,MPI_REAL,MPI_SUM,master_num,MPI_COMM_WORLD,ierr)
    !Put the total in array
    totalOfPrecipitation(nd)=masterPrecip
  enddo

  if(mchnum==master_num) then
    ! filename is the nickname of dam and initial date and time
    CALL makefnam(outFileName,trim(damOutPrefix), 0, &
        iyear1, imonth1,idate1,itime1*100,'precip', '$', 'csv')
    !DEtermine the advance in time
    call date_add_to(iyear1,imonth1,idate1,itime1*100,time,'s',oyr,omn,ody,otm)
    !at first time create the file and write date, hour and ammount of precip.
    if(firstTime) then
      ! Open file in csv format
      open(unit=funit,file=trim(outFileName) &
        ,action='WRITE',status='REPLACE',form='FORMATTED')
      !write a header in first line
      write(funit,fmt=fmt1) (trim(dams(nd)%damName),nd=1,nDams)
      !write the fields with precipitation for each dam
      write(funit,fmt=fmt2) oyr,omn,ody,otm &
         ,(totalOfPrecipitation(nd),nd=1,nDams)
      close(funit)
      !in other times just append the file for each dam
    else
      open(unit=funit,file=trim(outFileName) &
          ,action='WRITE',status='OLD',form='FORMATTED',POSITION='APPEND')
      write(funit,fmt=fmt2) oyr,omn,ody,otm &
         ,(totalOfPrecipitation(nd),nd=1,nDams)
      close(funit)
    endif
  endif

  firstTime=.false.

end subroutine outputDamPrecip




end module dam
