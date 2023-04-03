!=============================================================================================
program pre
    !# Program to perform preprocessing for BRAMS Model
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: This program get data from GRIB2 data from GFS Model and data from 
    !#  CAMS chemical climatology, adjust for the same grid defined by a prepBrams.ini file,
    !#  and concatenate the data in output defined file in grads (bin) or netcdf format.  
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
    
    use utilsMod, only: &
      bramsHeader, &
      fileExist, &
      initAll, &
      stepsBetweenDates, &
      monthsBetweenDates, &
      to_upper, &
      getUnit, &
      releaseUnit

    use filesMod, only: &
      readNamelist, &
      createGrib2FilesNames, &
      createEra5FilesNames, &
      readAtmosGrib2, &
      createCamsFilesNames, &
      readCamsFile, &
      readAtmosEra5

    use memoryMod, only: &
      init_year, &
      init_month, &
      init_day, &
      init_hour, &
      final_year, &
      final_month, &
      final_day, &
      final_hour, &
      step, &
      chemDate, &
      atmosDate, &
      atmos_Type, atmosNx,atmosNy,atmosNz

    use chemMod, only: &
      createChemEquivalence, &
      fillChemArrays

    use engineMod, only: &
      concatenate


    implicit none

    include "constants.f90"
    character(len=*),parameter :: sourceName='pre.f90' !Name of source code
    character(len=*),parameter :: procedureName='**pre**' !Name of this procedure
    !
    !Local Parameters
    character(len=*), parameter :: namelistFile='pre.nml'
    !# Namelist Filename
    character(len=*), parameter :: revision='0.1'
    !# Program revision
    integer, parameter :: gfs=1
    !# GFS type id
    integer, parameter :: era5=2
    !# era5 type id
    

    !Local variables
    integer :: stepsBetDates 
    !# Total of steps between initial and final dates
    logical :: monthBetDates(12)
    integer :: its,i,j,k,iTime
    integer, allocatable :: validDates(:)
    !# time counter
    character(len=15) :: warningFileName
    character(len=8) :: dateNow
    character(len=10) :: timeNow
    character(len=5) :: zoneNow
    integer :: valuesNow(8)
    integer :: lunit


    !Code

    !plot the program header
    iErrNumber=bramsHeader(revision)

    !Initialize utils module
    iErrNumber=initAll()
    
    !Check if namelist file is in current folder
    !then, if exist, read the namelist
    if(.not. fileExist(namelistFile)) then 
      iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
              ,c_fatal,'File '//namelistFile &
              //' not found. Please, verify and solve it!')
    else
      iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Reading '//namelistFile//'!')
      call readNamelist(namelistFile)
    endif

    !Compute how many steps are between initial and final date
    !This step number will be use to compose grib2 files names
    stepsBetDates=stepsBetweenDates(init_year,init_month,init_day,init_hour, &
                          final_year,final_month,final_day,final_hour, &
                          step,'h')

    !Allocating validDates to be processed 
    allocate(validDates(stepsBetDates))

    !Compute the months presents between initial and final date and put it in logical array
    !This months will be used to get the correct cams climatology file name
    monthBetDates=monthsBetweenDates(init_year,init_month,init_day,init_hour, &
                          final_year,final_month,final_day,final_hour, &
                          step,'h')

    select case  (atmos_type)
    case (gfs)
      iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Creating all grib2 files Names!')
      !Create all grib2 files Names and store in memory module
      call createGrib2FilesNames(stepsBetDates)
    case (era5)
      iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Creating all era5 files Names!')
      !Create all era5 files Names and store in memory module
      call createEra5FilesNames(stepsBetDates)
    case default
      iErrNumber=dumpMessage(c_tty,c_yes,'','',c_fatal,'atmos_type not recognized, see namelist!')
    end select

    iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Creating all cams files Names!')
    !Create all cams files names
    call createCamsFilesNames(monthBetDates)

    !Make equivalence for chem variables
    !Variables from chem1_list are readed and checked if is present
    !After it a cross reference between species parts and CAMS species are performed
    !To be used in concatenation
    call createChemEquivalence(.true.)
    
    !Read cams files using specific species for each file
    !The file 1 is the main choice. If specie not there
    !2nd file will be used.
    call readCamsFile() !Cams created by SRF

    select case  (atmos_type)
    case (gfs)
      !Read all grib2 files and cut it for given lats and lons
      !Set in namelist
      do its=1,stepsBetDates
        call readAtmosGrib2(its,stepsBetDates)
      enddo
    case (era5)
      !Read all era5 files and cut it for given lats and lons
      !Set in namelist
      do its=1,stepsBetDates
        call readAtmosEra5(its,stepsBetDates)
      enddo     
    end select

    !Write inventory of dates from Grib2 files
    iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Grib2 dates inventory:')
    write (*,fmt='("+--+----+--+--+--+")') 
    do i=1,stepsBetDates
      write (*,fmt='("|",I2.2,"|",I4.4,"|",I2.2,"|",I2.2,"|",I2.2,"|")') i,atmosDate(i)%year,atmosDate(i)%month &
        ,atmosDate(i)%day,atmosDate(i)%hour
      write (*,fmt='("+--+----+--+--+--+")')
    enddo

    !Concatenates atmos data and chem data, put it in same grid
    !And writes output
    do iTime=1,stepsBetDates
      call concatenate(iTime)
      write(warningFileName,fmt='(I4.4,I2.2,I2.2,I2.2,".done")')  atmosDate(ITime)%year,atmosDate(Itime)%month &
        ,atmosDate(iTime)%day,atmosDate(iTime)%hour
      lunit=getUnit()
      call date_and_time( dateNow, timeNow, zoneNow, valuesNow )
      open(unit=lunit, file=warningFileName)
      write(lunit,fmt='(A8,":",A10,":",A5)') dateNow, timeNow, zoneNow
      ierrNumber=releaseUnit(lunit)
    enddo

end program pre 