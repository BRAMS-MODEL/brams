!=============================================================================================
module memoryMod
    !# Memory of all common variables
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: Memory of all common variablesplete
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 26 August 2020 (Wednesday)
    !# @endnote
    !#
    !# @changes
    !# &#9744; <br/>
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
    
    !Use area
    use dump

    implicit none

    character(len=*),parameter :: procedureName='**memory**' !Name of this procedure
    !
    !Local Parameters

    !Global variables
    type tim 
        integer :: year
        !# year
        integer :: month 
        !# month
        integer :: day 
        !# day
        integer :: hour
        !# hour
    end type tim
    type ctim 
        integer :: month
        !# month
        integer :: hour
        !# hour of day
    end type ctim
    type cm 
        integer :: nx
        !#
        integer :: ny
        !#
        integer :: nz
        !#
        integer :: nt
        !#
        integer :: nv
        !#
        integer :: nm 
        !#
        real,allocatable :: lat(:)
        !#
        real,allocatable :: lon(:)
        !#
        real :: dLat
        !#
        real :: dLon
        !#
        real,allocatable :: hoursFrom1900(:)
        !# 
        real, allocatable :: values(:,:,:,:,:,:) !nm,nt,nv,nx,ny,nz
        !#
        real,allocatable :: levels(:)
        !#
        character(len=16),allocatable :: varName(:)
        !#
        character(len=256), allocatable :: longName(:)
        !#
        character(len=32), allocatable :: units(:)
        !#
    end type cm
    type(cm) :: cams

    integer :: init_year
    !# Namelist initial year
    integer :: init_month
    !# Namelist initial month
    integer :: init_day
    !# Namelist initial day
    integer :: init_hour
    !# Namelist initial hour
    integer :: final_year
    !# Namelist final year
    integer :: final_month
    !# Namelist final Month
    integer :: final_day
    !# Namelist Final day
    integer :: final_hour
    !# Namelist final hour
    integer :: step
    !# Namelist step among input in hours
    integer :: atmos_type
    !# Namelist type of atmosphere input (0-dprep, 1-grib2 GFS, 2-NetCDF GEOS)
    character(len=256) :: atmos_prefix
    !# Namelist prefix of atmosphere input
    character(len=256) :: atmos_sufix
    !# Namelist suffix of atmosphere input
    character(len=256) :: atmos_idir
    !# Namelist atmosphere folder
    integer :: levels
    !# Namelist number of levels to be used
    real :: initial_latitude
    !# Namelist initial latitude to be used
    real :: final_latitude
    !# Namelist final latittude to be used
    real :: initial_longitude
    !# Namelist initial latittude to be used
    real :: final_longitude
    !# Namelist final longitude to be used
    integer :: chem_type
    !# Namelist type of chemistry data
    character(len=256) :: chem_idir
    !# Namelist chem file folder
    integer :: out_type
    !# Namelist type of output (0 - text, 1-Grads, 2 - vfm, 3 - NetCDF)
    character(len=256) :: out_prefix
    !# Namelist output file prefix
    character(len=256) :: out_sufix
    !# Namelist output file suffix
    character(len=256) :: out_dir
    !# Namelist output folder
    character(len=256) :: chem1_prefix
    !# Prefix of 1st climatology cams file
    character(len=256) :: chem1_sufix 
    !# sufix of 1st climatology cams file

    character(len=256), allocatable :: grib2FilesNames(:)
    character(len=256), allocatable :: era5FilesNames(:)
    character(len=256), allocatable :: grib2InvFilesNames(:)

    integer :: atmosNz
    !# Number of atmospheric levels
    integer :: atmosNx
    !# Number of atmospheric longitude points
    integer :: atmosNy
    !# Number of atmospheric lattitude points
    real :: atmosLon(2)
    !# First and last longitude coordinates
    real :: atmosLat(2)
    !# First and last latittude coordinates
    real :: atmosXs
    !# Delta longitude
    real :: atmosYs
    !# Delta latitude
    integer, parameter :: atmosNv=5
    !# Number of atmospheric variables
    real, allocatable :: atmosLevels(:) 
    !# atmospheric Leves values  
    character(len=32) :: atmosVarNames(atmosNv)
    !# Name ofatmospheric variables
    real, allocatable :: atmosValues(:,:,:,:,:) 
    !# Value of atmospheric variables (ntime,nvar,nx,ny,nz)
    type(tim), allocatable :: atmosDate(:)
    integer, allocatable :: atmosHoursCount(:)

    integer :: monthCount
    !# Total of months present on data
    character(len=256), allocatable :: cams1FilesNames(:)


    type(ctim),allocatable :: chemDate(:,:)
    !#

    character(len=32),allocatable :: spcName(:)
    !# Chem specie name from equivalence file
    character(len=32),allocatable :: spcCamsName(:,:)
    !# Chem species parts for each chem species
    real, allocatable :: factor(:,:)
    !# Product factor for each species parts
    integer, allocatable :: whichCams(:,:)
    !# Number of file that contains the chem part

    contains

end module memoryMod 