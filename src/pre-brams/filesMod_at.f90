!=============================================================================================
module filesMod
    !# Module to manipulate namelist file
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: manipulate namelist file
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

    include "constants.f90"
    character(len=*),parameter :: sourceName='filesMod.f90' !Name of source code
    character(len=*),parameter :: procedureName='**namelist**' !Name of this procedure
    !
    !Local Parameters

    !Local variables
    integer :: chemjIni,chemjFin,chemnpj
    integer :: chemiIni,chemiFin,chemnpi

    contains

    !=============================================================================================
    subroutine readNamelist(namelistFile)
        !# read namelist file from input
        !#
        !# @note
        !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
        !#
        !# **Brief**: Read the namelist file, pre.nml and return data
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
        use memoryMod, only: &
                  init_year ,  init_month,  init_day,  init_hour,  &
                  final_year, final_month, final_day, final_hour,  &
                  step,                                            &
                  atmos_type, atmos_prefix,atmos_sufix,atmos_idir, &
                  levels, initial_latitude,final_latitude,         &
                  initial_longitude,final_longitude,               &
                  chem_type, chem_interp,chem_ppbmconv, chem_idir, &
                  chem1_prefix,chem1_sufix,chem2_prefix,chem2_sufix, & 
                  chem3_prefix,chem3_sufix, &
                  out_type, out_prefix, out_sufix, out_dir
        use utilsMod, only: &
                  getUnit, &
                  releaseUnit
    
        implicit none
    
        include "constants.f90"
        character(len=*),parameter :: procedureName='**readNamelist**' !Name of this procedure
        !
        !Local Parameters
    
        !Input/Output variables
        character(len=*), intent(in) :: namelistFile
    
        !Local variables
        integer :: lunit
    
        !Code
        NAMELIST/args_input/init_year ,  init_month,  init_day,  init_hour,  &
                            final_year, final_month, final_day, final_hour,  &
                            step,                                            &
                            atmos_type, atmos_prefix,atmos_sufix,atmos_idir, &
                            levels, initial_latitude,final_latitude,         &
                            initial_longitude,final_longitude,               &
                            chem_type, chem_interp,chem_ppbmconv, chem_idir, &
                            chem1_prefix,chem1_sufix,&
                            chem2_prefix,chem2_sufix, & 
                            chem3_prefix,chem3_sufix, & 
                            out_type, out_prefix, out_sufix, out_dir

         lunit=getUnit()
         open(unit=lunit, file=trim(namelistFile), status='old')
         read(lunit, args_input)
         ierrNumber=releaseUnit(lunit)


    end subroutine readNamelist 

    !=============================================================================================
    subroutine createGrib2FilesNames(numberOfSteps)
        !# Create a set of files Names for Grib2
        !#
        !# @note
        !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
        !#
        !# **Brief**: Creates a set of filenames for GribFiles accordingly the initial date and
        !# others informations from Namelist
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
        use memoryMod, only: &
            atmos_prefix, &
            atmos_sufix, &
            atmos_idir, &
            step, &
            grib2FilesNames, &
            grib2InvFilesNames, &
            atmosDate, &
            atmosHoursCount
    
        implicit none
    
        include "constants.f90"
        character(len=*),parameter :: procedureName='**createGrib2FilesNames**' !Name of this procedure
        !
        !Local Parameters
        character(len=69) :: line='+--+----------------------------------------------------------------+'
    
        !Input/Output variables
        integer, intent(in) :: numberOfSteps
    
        !Local variables
        integer :: i

        if(allocated(grib2FilesNames)) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
              ,c_fatal,' grib2FilesNames array already allocated with size=',size(grib2FilesNames),"I2.2")

        allocate(grib2FilesNames(numberOfSteps))

        if(allocated(grib2InvFilesNames)) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
              ,c_fatal,' grib2InvFilesNames array already allocated with size=',size(grib2InvFilesNames),"I2.2")

        allocate(grib2InvFilesNames(numberOfSteps))
        allocate(atmosHoursCount(numberOfSteps))

        !Code

        do i=1,numberOfSteps
            write(grib2FilesNames(i),fmt='(A,I3.3,A)') trim(atmos_prefix),(i-1)*step,trim(atmos_sufix)
            write(grib2InvFilesNames(i),fmt='(A,I3.3,A)') trim(atmos_prefix),(i-1)*step,trim(atmos_sufix)//'.inv'
            atmosHoursCount(i)=(i-1)*step
            !print *,trim(grib2FilesNames(i))
            !print *,trim(grib2InvFilesNames(i))
        enddo

        iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Grib2 file names inventory:')
        write(*,fmt='(A)') line
        do i=1,numberOfSteps
            write(*,fmt='("|",I2.2,"|",A64,"|")') i,trim(grib2FilesNames(i))
            write(*,fmt='(A)') line
        enddo

        if(allocated(atmosDate)) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
              ,c_fatal,' atmosDate array already allocated with size=',size(atmosDate),"I2.2")
        allocate(atmosDate(numberOfSteps))
          
    
    end subroutine createGrib2FilesNames 

    !=============================================================================================
    subroutine readAtmosGrib2(iStep,stepsBetDates)
        !# Open a grib2 file and read the contents
        !#
        !# @note
        !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
        !#
        !# **Brief**: Open GFS Grib2 file, read the contents
        !#
        !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
        !#
        !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
        !#
        !# **Date**: 31 July 2020 (Friday)
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

        use wgrib2api

        use memoryMod, only: &
            atmos_idir, &
            grib2FilesNames, &
            grib2InvFilesNames, &
            initial_latitude, &
            final_latitude, &
            initial_longitude, &
            final_longitude, &
            atmosNz, &
            atmosNx, &
            atmosNy, &
            atmosLon, &
            atmosLat, &
            atmosXs, &
            atmosYs, &
            atmosNv, &
            atmosLevels, &
            atmosVarNames,&
            atmosValues  , &
            atmosDate, &
            atmosHoursCount, &
            levels

        use utilsMod, only: &
                  getUnit, &
                  releaseUnit, &
                  fileExist, &
                  date_add_to_dble


        implicit none
    
        include "constants.f90"
        character(len=*),parameter :: procedureName='**openatmosGrib2**' 
        !# Name of this procedure
        !
        !Local Parameters
        character(len=*),parameter :: wind_u_varname='UGRD'
        !# Name of wind u var in grib2 (1)
        character(len=*),parameter :: wind_v_varname='VGRD'
        !# Name of wind v var in grib2 (2)
        character(len=*),parameter :: temperature_varname='TMP'
        !# Name of temperature var in grib2 (3)
        character(len=*),parameter :: geo_varname='HGT'
        !# Name of geopotential var in grib2 (4)
        character(len=*),parameter :: ur_varname='RH'
        !# Name of relative humidity var in grib2 (5)
        real,dimension(5),parameter :: scale_factor = (/1.,1.,1.,1.,0.01/)
        !# Factor to multiplie each var (1,2,3,4,5)
        real,dimension(5), parameter :: limi=(/-100.,-100.,180., -500.,  0./)
        !# inferior limits for each variable
        real,dimension(5), parameter :: lims=(/ 100., 100.,333.,80000.,100./)
        !# inferior limits for each variable
    
        integer, parameter :: maxpr=100
    
        !Input/Output variables
        integer, intent(in)         :: iStep
        !# timestep
        integer, intent(in)         :: stepsBetDates
        !# total of timesteps
    
        !Local variables
        integer :: nprz
        !# Total of vertical levels
        character(len=32) :: varName(5)
        !# Name os GFS variables
        character(len=200) :: metadata
        !# Metadata from grib2 information
        character(len=30), allocatable :: slevs(:)
        !# Levels description from grib2 file
        character(len=30), allocatable :: slevsInv(:)
        !# Levels description INverted
        real,allocatable :: levpr_grib2(:)
        !# levels values from grib2 file
        character (len=300) :: grid_info
        !# grid information of grib2 file
        character (len=99) :: invline
        !# Var information from grib2 file
        real, allocatable :: var(:,:)
        !# Variable value for each (x,y) point
        real, allocatable :: lat(:,:)
        !# lattitude value for each (x,y) point
        real, allocatable :: lon(:,:)
        !# longitude value for each (x,y) point
        integer :: nxGrib
        !# Size of grib data in lon direction
        integer :: nyGrib
        !# Size of grib data in lat direction
        real :: dx
        !# deltaX longitude increment
        real :: dy
        !# deltaY longitude increment 
        integer :: iyear,imonth,iday,ihour
        !# year, month, day and hour read from grib2 file
        real :: latIni,latFin
        !# Initial and final valid latittudes
        real :: lonIni,lonFin
        !# Initial and final valid longitudes
        integer :: jIni,jFin
        !# first and last latittude valid position
        integer :: iIni,Ifin
        !# first and last longitude valid position
        integer :: npi,npj
        !# Total of valid points i,j
        character(len=3) :: clvi
        integer :: iyy,imm,idd,ihh
    
        integer :: i,j,k,lv,nvar,iCount,jCount,ii,jj,lvi
        integer :: lunit

        character(len=256) :: fName,fNameI,fNameB
        character(len=3) :: ncl

        fName=trim(atmos_idir)//trim(grib2FilesNames(iStep))
        fNameI=trim(grib2InvFilesNames(iStep))
        fNameB=trim(grib2FilesNames(iStep))//'.blow'
    
        ! Read the header of input pressure file.
        varName(1)=':'//temperature_varname//':'
    
        !Code

        if(.not. fileExist(trim(fName))) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
              ,c_fatal,'File '//trim(fName) &
              //' not found. Please, verify and solve it!')

        iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Opening/reading Grib2 file: '//trim(fName))
        
        !Making inventory using a wgrib2 function grb2_mk_inv   
        iErrNumber = grb2_mk_inv(fName,fNameI)
        !print *,iErrNumber
        !if (iErrNumber.ne.1) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
        !           ,c_fatal,' Error making grib2 inventory, check it!')
    
        ! get number of first var levels using grb2_inq from wgrib2 lib
        nprz = grb2_inq(fName,fNameI,trim(varname(1)),' mb:')
        allocate(slevs(nprz),levpr_grib2(nprz),slevsInv(nprz))
        !
        !print *,'nprz=',nprz
        !
        do i = 1,nprz
          ! get pressure leves and levels labels using grb2_inq from wgrib2 lib
          iErrNumber=grb2_inq(fName,fNameI,trim(varname(1)),' mb:' &
                   ,sequential=i-1,desc=metadata)
          if (iErrNumber.ne.1) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
                   ,c_fatal,' Error getting pressure levels from grib2 file, check it!')
    
          j = index(metadata,trim(varname(1))) + len(trim(varname(1)))
          k = index(metadata," mb:") + len(" mb:")-1
          !Transform text values to real for each level
          read(metadata(j:),*) levpr_grib2(i)
          !Fill array with level name+suffix
          slevs(i) = metadata(j-1:k)
        enddo
        !
        !print *,index(metadata," mb:"),index(metadata,'hour')
        !print *,'hour: ',metadata(index(metadata," mb:")+4:index(metadata,'hour')-1)
        !print *,'slevs=',slevs
        !print *,'levpr=',levpr_grib2
        !
        !Getting size and geo information in GRIB2 file. Using the first variable (U)
        iErrNumber = grb2_inq(fName,fNameI,trim(varName(1)) &
                   , data2=var &
                   , lat=lat, lon=lon, grid_desc=grid_info, desc=invline)
    
        !Closing grib2 file to save memory
        iErrNumber = grb2_free_file(fName)
    
        !print *,grid_info
        !print *,invline
    
        !print *,var
    
        nxGrib=size(var,1)
        nyGrib=size(var,2)
    
        !print *,nxGrib,nyGrib
    
        ! Increments lat and lon
        dy=(lat(1,nyGrib)-lat(1,1))/(nyGrib-1)
        dx=(lon(nxGrib,1)-lon(1,1))/(nxGrib-1)
    
        !print *,dx,dy
    
        do j=1,nyGrib
            if(lat(1,j)>initial_latitude) then
                latIni=lat(1,j)
                jIni=j
                exit 
            endif
        enddo
        do j=jINi,nyGrib
            if(lat(1,j)<final_latitude) then
                latFin=lat(1,j)
                jFin=j 
            endif
        enddo
        do i=1,nxGrib
            if(lon(i,1)>initial_longitude) then
                lonIni=lon(i,1)
                iIni=i
                exit 
            endif
        enddo
        do i=iINi,nxGrib
            if(lon(i,1)<final_longitude) then
                lonFin=lon(i,1)
                iFin=i 
            endif
        enddo
        npi=iFin-iIni+1
        npj=jFin-jIni+1
        !print *,jIni,jFin,latIni,latFin,npj
        !print *,iIni,iFin,lonIni,lonFin,npi
    
        !convert date string got from grib2 file to integer
        read(invline( 3: 6),*,iostat=iErrNumber)  iyear
        read(invline( 7: 8),*,iostat=iErrNumber)  imonth
        read(invline( 9:10),*,iostat=iErrNumber)  iday
        read(invline(11:12),*,iostat=iErrNumber)  ihour

        call date_add_to_dble(iyear,imonth,iday,ihour,dble(atmosHoursCount(iStep)),'h' &
                       ,iyy,imm,idd,ihh)
        atmosDate(iStep)%year=iyy
        atmosDate(iStep)%month=imm
        atmosDate(iStep)%day=idd
        atmosDate(iStep)%hour=ihh/10000
    
        !print *,iyear,imonth,iday,ihour
    
        !Adjustin dpatmos type vars
        atmosNz=levels
        atmosNx=npi
        atmosNy=npj                  
        atmosLon(1)=lonIni-360.0
        atmosLat(1)=latIni
        atmosXs=dx
        atmosYs=dy
        if(.not. allocated(atmosLevels)) allocate(atmosLevels(atmosNz))
        do lv=1,nprz
          lvi=nprz-lv+1
          if(lvi>atmosNz) cycle
          atmosLevels(lvi)=levpr_grib2(lv)
          slevsInv(lvi)=slevs(lv)
        enddo
        
        if(iStep==1) then
            atmosLat(2) = atmosLat(1) + (atmosYs*atmosNy)
            atmosLon(2) = atmosLon(1) + (atmosXs*atmosNx)
            iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Grib2 sizes and levels inventory: ')
            write(*,fmt='("Ntimes: ",I3.3," Nvars: ",I3.3)') 1,5
            write(*,fmt='("NLons : ",I3.3," NLats: ",I3.3," NLevs : ",I3.3)') atmosNx,atmosNy,atmosNz
            write(*,fmt='("LatI  : ",F8.2," LatF : ",F8.2," DeltaY: ",F6.3)') atmosLat(1),atmosLat(2),atmosYs
            write(*,fmt='("LonI  : ",F8.2," LonF : ",F8.2," DeltaX: ",F6.3)') atmosLon(1),atmosLon(2),atmosXs
            write(ncl,fmt='(I3.3)') atmosNz
            write(*,fmt='("Levels: ",'//ncl//'(F6.1,1X))') atmosLevels             
        endif        
        
        atmosVarNames = (/wind_u_varname,wind_v_varname,temperature_varname &
                             ,geo_varname,ur_varname/)
             
        if(.not. allocated(atmosValues)) allocate(atmosValues(stepsBetDates,atmosNv,atmosNx, atmosNy, atmosNz))

        lunit=getUnit()
        open(unit=lunit,file=fNameB,action='write',status='replace')
    
        do lv=1,levels!nprz
          !lvi=lv !nprz-lv+1
          !if(lvi>atmosnZ) cycle
          write(clvi,fmt='(I3.3)') lv
          !Compose name of var to search inside grib2
          varName(1)=':'//trim(wind_u_varname)//trim(slevsInv(lv))    
          varName(2)=':'//trim(wind_v_varname)//trim(slevsInv(lv)) 
          varName(3)=':'//trim(temperature_varname)//trim(slevsInv(lv)) 
          varName(4)=':'//trim(geo_varname)//trim(slevsInv(lv))   
          varName(5)=':'//trim(ur_varname)//trim(slevsInv(lv)) 
          write(*,fmt='(I1,5(A,1X))') lv,trim(varName(:))
          do nvar=1,5            
            !Getting values for variable nvar 
            iErrNumber = grb2_inq(fName,fNameI,trim(varName(nVar)) &
                   ,data2=var, lat=lat, lon=lon)
            ii=0
            jj=0
            !Fill values only inside valid area
            do i=iIni,iFin
              ii=ii+1
              jj=0
              do j=jIni,jFin
                jj=jj+1
                atmosValues(iStep,nVar,ii,jj,lv)=var(i,j)
                !Verify var limits
                if(var(i,j)<limi(nvar)) &
                   iErrNumber=dumpMessage(lunit,c_yes,sourceName,procedureName &
                   ,c_warning,''//' - lev='//clvi//' - Estouro de limite inf. da var. ' &
                   //trim(atmosVarNames(nvar)),(/var(i,j),limi(nvar)/),'E13.4')
                if(var(i,j)>lims(nvar)) &
                   iErrNumber=dumpMessage(lunit,c_yes,sourceName,procedureName &
                   ,c_warning,''//' - lev='//clvi//' - Estouro de limite sup. da var. ' &
                   //trim(atmosVarNames(nvar)),(/var(i,j),lims(nvar)/),'E13.4')
              enddo
            enddo
          enddo
        end do
    
        !Closing grib2 file to save memory
        !print *, 'Closing ',trim(file)
        iErrNumber = grb2_free_file(fName)
        iErrNumber = releaseUnit(lunit)
    
    end subroutine readAtmosGrib2 

    !=============================================================================================
    subroutine createCamsFilesNames(monthBetDates)
        !# Creates the set of Cams chem climatology names
        !#
        !# @note
        !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
        !#
        !# **Brief**: Creates the set of Cams chem climatology names and store in memoryMod
        !#
        !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
        !#
        !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
        !#
        !# **Date**: 27 August 2020 (Thursday)
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

        use memoryMod, only: &
            monthCount, &
            cams1FilesNames, &
            cams2FilesNames, &
            cams3FilesNames, &    
            chem1_prefix, &
            chem1_sufix, & 
            chem2_prefix, &
            chem2_sufix, & 
            chem3_prefix, &
            chem3_sufix, &
            chemDate

        use utilsMod, only: &
            to_upper
    
        implicit none
    
        include "constants.f90"
        character(len=*),parameter :: procedureName='**createCamsFilesNames**' !Name of this procedure
        character(len=71) :: line='+--+-+----------------------------------------------------------------+'
        !
        !Local Parameters
    
        !Input/Output variables
        logical, intent(in) :: monthBetDates(12)
        !# Array with true for each month present
    
        !Local variables
        integer :: i,iCount
        character(len=3) :: monthCamsName(12)

        !Code
        monthCount=0
        do i=1,12
            if(monthBetDates(i)) then
                monthCount=monthCount+1
                monthCamsName(monthCount)=to_upper(month_Name(i))
            endif
        enddo
        if(allocated(cams1FilesNames)) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
              ,c_fatal,' cams1FilesNames array already allocated with size=',size(cams1FilesNames),"I2.2")
        allocate(cams1FilesNames(monthCount))
        if(allocated(cams2FilesNames)) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
              ,c_fatal,' cams2FilesNames array already allocated with size=',size(cams2FilesNames),"I2.2")
        allocate(cams2FilesNames(monthCount))
        if(allocated(cams3FilesNames)) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
              ,c_fatal,' cams3FilesNames array already allocated with size=',size(cams3FilesNames),"I2.2")
        allocate(cams3FilesNames(monthCount))
        do i=1,monthCount
            write(cams1FilesNames(i),fmt='(A,A,A)') trim(chem1_prefix),monthCamsName(i),trim(chem1_sufix)
            !print *,trim(cams1FilesNames(i))
            write(cams2FilesNames(i),fmt='(A,A,A)') trim(chem2_prefix),monthCamsName(i),trim(chem2_sufix)
            !print *,trim(cams2FilesNames(i))
            write(cams3FilesNames(i),fmt='(A,A,A)') trim(chem3_prefix),monthCamsName(i),trim(chem3_sufix)
            !print *,trim(cams2FilesNames(i))
        enddo

        iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Cams file names inventory:')
        write(*,fmt='(A)') line
        do i=1,monthCount
            write(*,fmt='("|",I2.2,"|",I1,"|",A64,"|")') i,1,trim(cams1FilesNames(i))
            write(*,fmt='("|",A2,"|",I1,"|",A64,"|")') "",2,trim(cams2FilesNames(i))
            write(*,fmt='(A)') line
        enddo
    
    end subroutine createCamsFilesNames 

    !=============================================================================================
    subroutine readCamsFile()
        !# Read and store data from cams files
        !#
        !# @note
        !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
        !#
        !# **Brief**: REad and store data from camsFiles
        !#
        !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
        !#
        !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
        !#
        !# **Date**: 27 August 2020 (Thursday)
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

        use netcdf

        use memoryMod, only: &
            chem1Nz, &
            chem1Nx, &
            chem1Ny, &
            chem1Nt, &
            chem1Lon, &
            chem1Lat, &
            chem1Xs, &
            chem1Ys, &
            chem1Nv, &
            chem1Levels, &
            chem1VarNames, &
            chem1Values, &
            monthCount, &
            cams1FilesNames, &
            chem_idir, &
            chem1Offset, &
            chem1Factor, &
            initial_latitude, &
            final_latitude, &
            initial_longitude, &
            final_longitude, &
            chemDate, &
            chem1LongName, &
            chem1Units, &
            spcCamsName, &
            whichCams, &
            chemllat, &
            chemllon

        use utilsMod, only: &
            fileExist, &
            getUnit, &
            releaseUnit, &
            date_add_to_dble, &
            outRealSize

        use chemMod, only: &
            nSpc
    
        implicit none
    
        include "constants.f90"
        include "netcdf.inc"
        character(len=*),parameter :: procedureName='**readCamsFile**' !Name of this procedure
        !
        !Local Parameters
    
        !Input/Output variables
    
        !Local variables
        integer :: nFile 
        integer :: ncid
        character(len=256) :: fName
        character(len=32), allocatable :: varName(:)
        character(len=32) :: name,atName,cValueAt
        integer, allocatable :: varDim(:),nat(:),lenDim(:)
        integer :: levelVarn,lonVarN,latVarN,timeVarN
        integer :: nVars,nDims,an,nnt,ii,jj,nvar
        real, allocatable :: llat(:),llon(:)
        integer, allocatable :: HoursFrom1900(:)
        integer :: nx,ny,nz,nt,ki,k
        integer :: j
        integer :: i
        real :: latIni,latFin
        real :: lonIni,LonFin
        real, allocatable :: Values(:,:,:,:,:)
        double precision :: xValue 
        integer :: iLen,xtype,lunit
        integer :: iyy,imm,idd,ihh,nlev,levI
        real, allocatable :: localLevels(:)
        character(len=3) :: cnz
        integer :: recordLen, iRec
        character(len=256) :: localFileName,cValue
        character(len=15) :: tDef
        integer, allocatable :: varNum(:)
        integer :: varCount,ji,iCount,jCount


        !Code
        do nfile=1,monthCount
            fName=trim(chem_idir)//trim(cams1FilesNames(nFile))//'.nc'
            if(.not. fileExist(trim(fName))) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
              ,c_fatal,'File '//trim(fName) &
              //' not found. Please, verify and solve it!')
            iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Opening/reading cams file1: '//trim(fName))
            
            !Open the NetCDF file
            iErrNumber = nf90_open(path = trim(fName), mode = nf90_nowrite, ncid = ncid)
            ! get info about netCDF file 
            iErrNumber=nf90_inquire(ncid, ndims, nvars)
            !print *, ndims,nvars
            if(nFile==1) then
                allocate(varName(nvars),varDim(nvars),nat(nvars),lenDim(nvars))
            endif

            do i=1,nvars
                iErrNumber = nf90_Inquire_Dimension(ncid, i, name, lenDim(i))
                select case (trim(name))
                    case ('level') 
                        levelVarn=i
                        nz=lenDim(i)  
                    case ('longitude')
                        lonVarN=i
                        nx=lenDIm(i)
                    case ('latitude')
                        latVarN=i
                        ny=lenDIm(i)
                    case ('time')
                        timeVarN=i
                        nt=lenDIm(i)
                end select
                !print *,i,name,lenDIm(i)
            enddo
            if(nFile==1) then
                allocate(chemllat(ny),chemllon(nx))
                allocate(llat(ny))
                allocate(values(nt,nVars,nx,ny,nz))
                if(allocated(chemDate)) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
                    ,c_fatal,' chemDate array already allocated with size=',size(chemDate),"I2.2")
                allocate(chemDate(monthCount,nt))
            endif
            iErrNumber=nf90_get_var(ncid, latVarN, llat)
            iErrNumber=nf90_get_var(ncid, lonVarN, chemllon)  

            chem1Nx=nx
            chem1Ny=ny
            chem1Nz=nz
            chem1Nt=nt

            do j=1,ny
                ji=ny-j+1
                chemllat(ji)=llat(j)
            enddo
    
            iCount=0
            print *,initial_longitude,final_longitude
            do i=1,nx
                if(chemllon(i)>=initial_longitude .and. chemllon(i)<=final_longitude) then
                    if(icount==0) lonIni=chemllon(i)
                    iCount=iCount+1
                    lonFin=chemllon(i)
                endif
            enddo
            print *,initial_latitude,final_latitude
            jCount=0
            do i=1,ny
                if(chemllat(i)>=initial_latitude .and. chemllat(i)<=final_latitude) then
                    if(jCount==0) latIni=chemllat(i)
                    jCount=jCount+1
                    latFin=chemlLat(i)
                endif
            enddo
            
            print *,iCount-1,jCount-1
            print *,lonIni,latINi
            print *,lonFin,latFin
            print *,(lonFin-LonIni)/(iCount-1),(latFin-LatINi)/(jCount-1)
            stop 12121

            if(nFile==1) then
                !Allocate vars to read 
                allocate(localLevels(nz))
                allocate(HoursFrom1900(chem1Nt))
            endif

            !Read levels, dates, latitudes and longitudes
            iErrNumber=nf90_get_var(ncid, levelVarN, localLevels)
            iErrNumber=nf90_get_var(ncid, timeVarN, HoursFrom1900)


            if(nfile==1) then
                allocate(varNum(nvars))
                allocate(chem1Levels(chem1Nz))
            endif

            chem1Levels=localLevels
         
            !Put in memory vars
            chem1Lat(1)=chemllat(ny)
            chem1Lon(1)=chemllon(1)!-360.0
            chem1Lat(2)=chemllat(1)
            chem1Lon(2)=chemllon(nx)!-360.0
            chem1Xs=chemllon(2)-chemllon(1)
            chem1Ys=chemllat(2)-chemllat(1)


            do i=1,nt
                !write (*,fmt='(A,I2,A,I16)') "Hours for ",i," is ",HoursFrom1900(i)
                call date_add_to_dble(1900,1,1,0,dble(HoursFrom1900(i)),'h' &
                       ,iyy,imm,idd,ihh)
                chemDate(nFile,i)%month=imm
                !print *,'LFR: ',i,chemDate(nFile,i)%month,HoursFrom1900(i)
                chemDate(nFile,i)%hour=ihh/10000
                !print *,'LFR: ',i,chemDate(nFile,i)%hour,HoursFrom1900(i)
            enddo

            !Get var Names
            varCount=0
            do i=1,nvars
                iErrNumber=nf90_Inquire_Variable(ncid, i, name=varName(i) &
                ,ndims=varDim(i), nAtts=nat(i))
                !print *,i,chem1VarNames(i),varDim(i),nat(i)
                !Percorre as variaveis do chem1_list
                do j=1,nSpc
                    !Percorre até 5 membros que compoe a especie
                    do k=1,5
                        if(varName(i)==spcCamsName(j,k)) then !se eh a especie 
                            whichCams(j,k)=1 ! Marca que foi encontrada no arquivo 1
                            varCount=varCount+1 !ACrescenta uma variavel
                            varNum(varCount)=i !Guarda o numero dessa variavel
                            !print *,i,varName(i),varCount
                        endif
                    enddo
                enddo
            enddo
            chem1Nv=varCount
            if(nFile==1) then
                !Allocate vars to read 
                allocate(chem1VarNames(varCount))
                allocate(chem1Offset(varCount))
                allocate(chem1Factor(varCount))
                chem1Factor=1.0
                chem1Offset=0.0
                allocate(chem1Values(monthCount,chem1Nt,varCount,chem1Nx,chem1Ny,chem1Nz))
                allocate(chem1LongName(varCount),chem1Units(varCount))
            endif

            do i=1,varCount
                chem1VarNames(i)=varName(varNum(i))
                do an=1,nat(varNum(i))
                    iErrNumber=nf90_inq_attname(ncid, varNum(i), an, atName)
                    iErrNumber=nf90_Inquire_Attribute(ncid, varNum(i), trim(atName), xtype=xtype, len=iLen)
                    if(trim(atName)=='scale_factor') then
                        iErrNumber = nf90_get_att(ncid, varNum(i), atName, xValue)
                        chem1Factor(i)=xvalue
                    elseif(trim(atName)=='add_offset') then
                        iErrNumber = nf90_get_att(ncid, varNum(i), atName, xValue)
                        chem1Offset(i)=xvalue
                    elseif(trim(atName)=='long_name') then
                        iErrNumber = nf90_get_att(ncid, varNum(i), atName, cValue)
                        chem1LongName(i)=cValue
                    elseif(trim(atName)=='units') then
                        iErrNumber = nf90_get_att(ncid, varNum(i), atName, cValue)
                        chem1Units(i)=cValue
                    endif
                enddo
            end do
            iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Cams sizes and levels inventory for file #1: ')
            write(*,fmt='("Ntimes: ",I3.3," Nvars: ",I3.3)') chem1Nt,chem1Nv
            write(*,fmt='("NLons : ",I3.3," NLats: ",I3.3," NLevs : ",I3.3)') chem1Nx,chem1Ny,chem1Nz
            write(*,fmt='("LatI  : ",F8.2," LatF : ",F8.2," DeltaY: ",F6.3)') chem1Lat(1),chem1Lat(2),chem1Ys
            write(*,fmt='("LonI  : ",F8.2," LonF : ",F8.2," DeltaX: ",F6.3)') chem1Lon(1),chem1Lon(2),chem1Xs
            write(cnz,fmt='(I3.3)') chem1Nz
            write(*,fmt='("Levels: ",'//cnz//'(F6.1,1X))') chem1Levels
            write(*,fmt='("Dates : ",8(I7.7,1X))') HoursFrom1900

            !Reading the data.
            do nnt=1,chem1Nt
                iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Reading Cams #1 variables for t: ',nnt,'I2.2')
                do i=1,varCount
                    iErrNumber = nf90_get_var(ncid, varNum(i), Values(nnt,varNum(i),:,:,:) &
                       ,start = (/ 1, 1, 1,nnt /))
                    write(*,fmt='(I2.2,1X,A16,1X,A40,1X,A16,1X,2(E18.4,1X))') i,chem1VarNames(i),chem1LongName(i) &
                        ,chem1Units(i),maxval(Values(nnt,varNum(i),:,:,:)),minval(Values(nnt,varNum(i),:,:,:))
                enddo
            enddo

            do nnt=1,chem1Nt
                iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Copiando Cams #1 variables for t: ',nnt,'I2.2')
                do nvar=1,chem1Nv
                    do j=1,chem1Ny
                        ji=chem1Ny-j+1 !Invertendo a latitude
                        do k=1,chem1Nz
                            ki=chem1Nz-k+1 !Invertendo os niveis
                            chem1Values(nfile,nnt,nVar,:,ji,ki)=Values(nnt,varNum(nVar),:,j,k) &
                                *chem1Factor(nVar)+chem1Offset(nVar)
                        enddo
                    enddo
                    !write(*,fmt='(I2.2,1X,A16,1X,A40,1X,A16,1X,4(E18.4,1X))') nvar,chem1VarNames(nvar),chem1LongName(nvar) &
                    !    ,chem1Units(nvar),maxval(chem1Values(nfile,nnt,nVar,:,:,:)),minval(chem1Values(nfile,nnt,nVar,:,:,:)), &
                    !    chem1Factor(nVar),chem1Offset(nVar)
                enddo
            enddo
    
            iErrNumber = nf90_close(ncid)

            write(tDef,fmt='(I2.2,":00z",I2.2,A3,I4.4)')  0,1 &
                                                          ,month_Name(chemDate(nFile,1)%month),2019
            write(localFileName,fmt='("Leitura1-",I2.2,"-",I2.2)') nFile,1
            lunit=getUnit()
            open(unit=lunit,file=trim(localFileName)//'.ctl',action='write',status='replace')
            !writing the name of grads file
            write(lunit,*) 'dset ^'//trim(localFileName)//'.gra'
            !writing others infos to ctl
            write(lunit,*) 'undef -0.9990000E+34'
            write(lunit,*) 'title '//procedureName
            write(lunit,*) 'xdef ',chem1Nx,' linear ',chem1Lon(1),chem1Xs
            write(lunit,*) 'ydef ',chem1Ny,' linear ',chem1Lat(1),chem1Ys
            write(lunit,*) 'zdef ',chem1Nz,'levels',(chem1Levels(k),k=chem1Nz,1,-1)
            write(lunit,*) 'tdef 1 linear '//tDef//' 1hr'
            write(lunit,*) 'vars ',chem1Nv
            do i=1,chem1Nv
                write(lunit,*) trim(chem1VarNames(i)),chem1Nz,'99 ',trim(chem1LongName(i))//' '//trim(chem1Units(i))
            enddo
            write(lunit,*) 'endvars'
            ierrNumber=releaseUnit(lunit)
    
            recordLen=outRealSize()*chem1Nx*chem1Ny

            lunit=getUnit()
            open(unit=lunit,file=trim(localFileName)//'.gra',action='WRITE',status='REPLACE' &
                 ,form='UNFORMATTED',access='DIRECT',recl=recordLen)
    
            irec=1
            do nnt=1,chem1Nt
                do nvar=1,chem1Nv
                    do k=1,chem1Nz
                        write(lunit,rec=irec) chem1Values(nfile,nnt,nvar,:,:,k)
                        irec=irec+1
                    enddo
                enddo
            enddo

        enddo

    end subroutine readCamsFile 

    !=============================================================================================
    subroutine readCamsFile2()
        !# Read and store data from cams files
        !#
        !# @note
        !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
        !#
        !# **Brief**: REad and store data from camsFiles
        !#
        !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
        !#
        !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
        !#
        !# **Date**: 27 August 2020 (Thursday)
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

        use netcdf

        use memoryMod, only: &
            chem1Nz, &
            chem1Nx, &
            chem1Ny, &
            chem1Nt, &
            chem1Lon, &
            chem1Lat, &
            chem1Xs, &
            chem1Ys, &
            chem1Nv, &
            chem1Levels, &
            chem2VarNames, &
            chem2Values, &
            monthCount, &
            cams2FilesNames, &
            chem_idir, &
            chem2Offset, &
            chem2Factor, &
            initial_latitude, &
            final_latitude, &
            initial_longitude, &
            final_longitude, &
            chemDate, &
            chem2LongName, &
            chem2units, &
            spcCamsName, &
            whichCams, &
            chem2Nz, &
            chem2Nx, &
            chem2Ny, &
            chem2Nt, &
            chem2Lon, &
            chem2Lat, &
            chem2Xs, &
            chem2Ys, &
            chem2Nv, &
            chem2Levels


        use utilsMod, only: &
            fileExist, &
            getUnit, &
            releaseUnit, &
            date_add_to_dble, &
            outRealSize

        use chemMod, only: &
            nSpc
    
        implicit none
    
        include "constants.f90"
        include "netcdf.inc"
        character(len=*),parameter :: procedureName='**readCamsFile2**' !Name of this procedure
        !
        !Local Parameters
    
        !Input/Output variables
    
        !Local variables
        integer :: nFile 
        integer :: ncid
        character(len=256) :: fName
        character(len=32), allocatable :: varName(:)
        character(len=32) :: name,atName,cValueAt
        integer, allocatable :: varDim(:),nat(:),lenDim(:)
        integer :: levelVarn,lonVarN,latVarN,timeVarN
        integer :: nVars,nDims,an,nnt,ii,jj,nvar
        real, allocatable :: llat(:),llon(:)
        integer, allocatable :: HoursFrom1900(:)
        integer :: nx,ny,nz,nt,ki,k
        integer :: j
        integer :: i
        real :: latIni,latFin
        real :: lonIni,LonFin
        real, allocatable :: Values(:,:,:,:,:)
        double precision :: xValue 
        integer :: iLen,xtype,lunit
        integer :: iyy,imm,idd,ihh,nlev,levI
        real, allocatable :: localLevels(:)
        character(len=3) :: cnz
        integer :: recordLen, iRec
        character(len=256) :: localFileName,cValue
        character(len=15) :: tDef
        integer, allocatable :: varNum(:)
        integer :: varCount,ji


        !Code
        do nfile=1,monthCount
            fName=trim(chem_idir)//trim(cams2FilesNames(nFile))//'.nc'
            if(.not. fileExist(trim(fName))) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
              ,c_fatal,'File '//trim(fName) &
              //' not found. Please, verify and solve it!')
            iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Opening/reading cams file1: '//trim(fName))
            
            !Open the NetCDF file
            iErrNumber = nf90_open(path = trim(fName), mode = nf90_nowrite, ncid = ncid)
            ! get info about netCDF file 
            iErrNumber=nf90_inquire(ncid, ndims, nvars)
            !print *, ndims,nvars
            if(nFile==1) then
                allocate(varName(nvars),varDim(nvars),nat(nvars),lenDim(nvars))
            endif

            do i=1,nvars
                iErrNumber = nf90_Inquire_Dimension(ncid, i, name, lenDim(i))
                select case (trim(name))
                    case ('level') 
                        levelVarn=i
                        nz=lenDim(i)  
                    case ('longitude')
                        lonVarN=i
                        nx=lenDIm(i)
                    case ('latitude')
                        latVarN=i
                        ny=lenDIm(i)
                    case ('time')
                        timeVarN=i
                        nt=lenDIm(i)
                end select
                !print *,i,name,lenDIm(i)
            enddo
            if(nFile==1) then
            !    allocate(llat(ny),llon(nx))
                allocate(values(nt,nVars,nx,ny,nz))
            !    if(allocated(chemDate)) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
            !        ,c_fatal,' chemDate array already allocated with size=',size(chemDate),"I2.2")
            !    allocate(chemDate(monthCount,nt))
            endif
            !iErrNumber=nf90_get_var(ncid, latVarN, llat)
            !iErrNumber=nf90_get_var(ncid, lonVarN, llon)  

            chem2Nx=chem1nx
            chem2Ny=chem1ny
            chem2Nz=chem1nz
            chem2Nt=chem1nt
    
            if(nFile==1) then
                !Allocate vars to read 
                allocate(localLevels(nz))
                allocate(HoursFrom1900(chem1Nt))
            endif

            !Read levels, dates, latitudes and longitudes
            !iErrNumber=nf90_get_var(ncid, levelVarN, localLevels)
            !iErrNumber=nf90_get_var(ncid, timeVarN, HoursFrom1900)


            if(nfile==1) then
                allocate(varNum(nvars))
            !    allocate(chem1Levels(chem1Nz))
            endif

            chem2Levels=chem1Levels
         
            !Put in memory vars
            chem2Lat(1)=chem1Lat(1)
            chem2Lon(1)=chem1Lon(1)
            chem2Lat(2)=chem1Lat(2)
            chem2Lon(2)=chem1Lon(2)
            chem2Xs    =chem1Xs    
            chem2Ys    =chem1Ys    


            !do i=1,nt
            !    !write (*,fmt='(A,I2,A,I16)') "Hours for ",i," is ",HoursFrom1900(i)
            !    call date_add_to_dble(1900,1,1,0,dble(HoursFrom1900(i)),'h' &
            !           ,iyy,imm,idd,ihh)
            !    chemDate(nFile,i)%month=imm
            !    !print *,'LFR: ',i,chemDate(nFile,i)%month,HoursFrom1900(i)
            !    chemDate(nFile,i)%hour=ihh/10000
            !    !print *,'LFR: ',i,chemDate(nFile,i)%hour,HoursFrom1900(i)
            !enddo

            !Get var Names
            varCount=0
            do i=1,nvars
                iErrNumber=nf90_Inquire_Variable(ncid, i, name=varName(i) &
                ,ndims=varDim(i), nAtts=nat(i))
                !print *,i,chem1VarNames(i),varDim(i),nat(i)
                !Percorre as variaveis do chem1_list
                do j=1,nSpc
                    !Percorre até 5 membros que compoe a especie
                    do k=1,5
                        if(varName(i)==spcCamsName(j,k) .and. whichCams(j,k)==0) then !se eh a especie não marcada no arq 1 
                            whichCams(j,k)=2 ! Marca que foi encontrada no arquivo 2
                            varCount=varCount+1 !ACrescenta uma variavel
                            varNum(varCount)=i !Guarda o numero dessa variavel
                            !print *,i,varName(i),varCount
                        endif
                    enddo
                enddo
            enddo
            chem2Nv=varCount
            if(nFile==1) then
                !Allocate vars to read 
                allocate(chem2VarNames(varCount))
                allocate(chem2Offset(varCount))
                allocate(chem2Factor(varCount))
                chem2Factor=1.0
                chem2Offset=0.0
                allocate(chem2Values(monthCount,chem2Nt,varCount,chem2Nx,chem2Ny,chem2Nz))
                allocate(chem2LongName(varCount),chem2Units(varCount))
            endif

            do i=1,varCount
                chem2VarNames(i)=varName(varNum(i))
                do an=1,nat(varNum(i))
                    iErrNumber=nf90_inq_attname(ncid, varNum(i), an, atName)
                    iErrNumber=nf90_Inquire_Attribute(ncid, varNum(i), trim(atName), xtype=xtype, len=iLen)
                    if(trim(atName)=='scale_factor') then
                        iErrNumber = nf90_get_att(ncid, varNum(i), atName, xValue)
                        chem2Factor(i)=xvalue
                    elseif(trim(atName)=='add_offset') then
                        iErrNumber = nf90_get_att(ncid, varNum(i), atName, xValue)
                        chem2Offset(i)=xvalue
                    elseif(trim(atName)=='long_name') then
                        iErrNumber = nf90_get_att(ncid, varNum(i), atName, cValue)
                        chem2LongName(i)=cValue
                    elseif(trim(atName)=='units') then
                        iErrNumber = nf90_get_att(ncid, varNum(i), atName, cValue)
                        chem2Units(i)=cValue
                    endif
                enddo
            end do
            iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Cams sizes and levels inventory for file #1: ')
            write(*,fmt='("Ntimes: ",I3.3," Nvars: ",I3.3)') chem2Nt,chem2Nv
            write(*,fmt='("NLons : ",I3.3," NLats: ",I3.3," NLevs : ",I3.3)') chem2Nx,chem2Ny,chem2Nz
            write(*,fmt='("LatI  : ",F8.2," LatF : ",F8.2," DeltaY: ",F6.3)') chem2Lat(1),chem2Lat(2),chem2Ys
            write(*,fmt='("LonI  : ",F8.2," LonF : ",F8.2," DeltaX: ",F6.3)') chem2Lon(1),chem2Lon(2),chem2Xs
            write(cnz,fmt='(I3.3)') chem2Nz
            write(*,fmt='("Levels: ",'//cnz//'(F6.1,1X))') chem2Levels
            !write(*,fmt='("Dates : ",8(I7.7,1X))') HoursFrom1900

            !Reading the data.
            do nnt=1,chem2Nt
                iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Reading Cams #2 variables for t: ',nnt,'I2.2')
                do i=1,varCount
                    iErrNumber = nf90_get_var(ncid, varNum(i), Values(nnt,varNum(i),:,:,:) &
                       ,start = (/ 1, 1, 1,nnt /))
                    write(*,fmt='(I2.2,1X,A16,1X,A40,1X,A16,1X,2(E18.4,1X))') i,chem2VarNames(i),chem2LongName(i) &
                        ,chem2Units(i),maxval(Values(nnt,varNum(i),:,:,:)),minval(Values(nnt,varNum(i),:,:,:))
                enddo
            enddo

            do nnt=1,chem2Nt
                do nvar=1,chem2Nv
                    do j=1,chem2Ny
                        ji=chem2Ny-j+1 !Invertendo a latitude
                        do k=1,chem2Nz
                            ki=chem2Nz-k+1 !Invertendo os niveis
                            chem2Values(nfile,nnt,nVar,:,ji,ki)=Values(nnt,varNum(nVar),:,j,k) &
                                *chem2Factor(nVar)+chem2Offset(nVar)
                        enddo
                    enddo
                    !write(*,fmt='(I2.2,1X,A16,1X,A40,1X,A16,1X,4(E18.4,1X))') nvar,chem1VarNames(nvar),chem1LongName(nvar) &
                    !    ,chem1Units(nvar),maxval(chem1Values(nfile,nnt,nVar,:,:,:)),minval(chem1Values(nfile,nnt,nVar,:,:,:)), &
                    !    chem1Factor(nVar),chem1Offset(nVar)
                enddo
            enddo
    
            iErrNumber = nf90_close(ncid)

            write(tDef,fmt='(I2.2,":00z",I2.2,A3,I4.4)')  0,1 &
                                                          ,month_Name(chemDate(nFile,1)%month),2019
            write(localFileName,fmt='("Leitura2-",I2.2,"-",I2.2)') nFile,1
            lunit=getUnit()
            open(unit=lunit,file=trim(localFileName)//'.ctl',action='write',status='replace')
            !writing the name of grads file
            write(lunit,*) 'dset ^'//trim(localFileName)//'.gra'
            !writing others infos to ctl
            write(lunit,*) 'undef -0.9990000E+34'
            write(lunit,*) 'title '//procedureName
            write(lunit,*) 'xdef ',chem2Nx,' linear ',chem2Lon(1),chem2Xs
            write(lunit,*) 'ydef ',chem2Ny,' linear ',chem2Lat(1),chem2Ys
            write(lunit,*) 'zdef ',chem2Nz,'levels',(chem2Levels(k),k=chem2Nz,1,-1)
            write(lunit,*) 'tdef 1 linear '//tDef//' 1hr'
            write(lunit,*) 'vars ',chem2Nv
            do i=1,chem2Nv
                write(lunit,*) trim(chem2VarNames(i)),chem2Nz,'99 ',trim(chem2LongName(i))//' '//trim(chem2Units(i))
            enddo
            write(lunit,*) 'endvars'
            ierrNumber=releaseUnit(lunit)
    
            recordLen=outRealSize()*chem2Nx*chem2Ny

            lunit=getUnit()
            open(unit=lunit,file=trim(localFileName)//'.gra',action='WRITE',status='REPLACE' &
                 ,form='UNFORMATTED',access='DIRECT',recl=recordLen)
    
            irec=1
            do nnt=1,chem2Nt
                do nvar=1,chem2Nv
                    do k=1,chem2Nz
                        write(lunit,rec=irec) chem2Values(nfile,nnt,nvar,:,:,k)
                        irec=irec+1
                    enddo
                enddo
            enddo

        enddo

    end subroutine readCamsFile2


    !=============================================================================================
    subroutine writeGradsCtlFile(iTime,nSpc,chemValues,chemVarNames)
        !# Write the output in grads and ctl files
        !#
        !# @note
        !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
        !#
        !# **Brief**: write output in grads and ctl file with meteorological and chem data
        !#
        !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
        !#
        !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
        !#
        !# **Date**: 27 August 2020 (Thursday)
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

        use memoryMod, only: &
            out_type, &
            out_prefix, &
            out_sufix, &
            out_dir, &
            atmosNx, &
            atmosNy, &
            atmosNz, &
            atmosLat, &
            atmosLon, &
            atmosXs, &
            atmosYs, &
            atmosValues, &
            atmosVarNames, &
            atmosDate, &
            atmosNv, &
            atmosLevels  

        use utilsMod, only: &
            getUnit , &
            releaseUnit, &
            outRealSize
    
        implicit none
    
        include "constants.f90"
        character(len=*),parameter :: procedureName='**writeGradsCtlFile**' !Name of this procedure
        !
        !Local Parameters
    
        !Input/Output variables
        integer :: iTime
        !#
        integer, intent(in) :: nSpc
        !#
        real, intent(in) :: chemvalues(nSpc,atmosNx,atmosNy,atmosNz)
        !#
        character(len=*) :: chemVarNames(nSpc)
        !#

        !Local variables
        character(len=256) :: ctlFileName
        character(len=256) :: binFileName
        character(len=15) :: tDef
        integer :: TotalVars,nVar,lunit,recordLen
        integer :: iRec,nz
    
        !Code
        !SUm the amount of vars
        TotalVars=atmosNv+nSpc
        !Creating output names
        write(ctlFileName,fmt='(A,I4.4,I2.2,I2.2,I2.2,A)') trim(out_prefix),atmosDate(ITime)%year,atmosDate(Itime)%month &
                                    ,atmosDate(iTime)%day,atmosDate(iTime)%hour,trim(out_sufix)//'.ctl'
        write(binFileName,fmt='(A,I4.4,I2.2,I2.2,I2.2,A)') trim(out_prefix),atmosDate(ITime)%year,atmosDate(Itime)%month &
                                    ,atmosDate(iTime)%day,atmosDate(iTime)%hour,trim(out_sufix)//'.gra'

        write(tDef,fmt='(I2.2,":00z",I2.2,A3,I4.4)')  atmosDate(iTime)%hour,atmosDate(iTime)%day &
                                                      ,month_Name(atmosDate(Itime)%month),atmosDate(ITime)%year
        !Writing ctl file
        lunit=getUnit()
        open(unit=lunit, file=trim(out_dir)//trim(ctlFileName), action='write', status='replace')

        !writing the name of grads file
        write(lunit,*) 'dset ^'//trim(binFileName)
        !writing others infos to ctl
        write(lunit,*) 'undef -0.9990000E+34'
        write(lunit,*) 'title PRE BRAMS data'
        write(lunit,*) 'xdef ',atmosNx,' linear ',atmosLon(1),atmosXs
        write(lunit,*) 'ydef ',atmosNy,' linear ',atmosLat(1),atmosYs
        write(lunit,*) 'zdef ',atmosNz,'levels',atmosLevels(1:atmosNz)
        write(lunit,*) 'tdef 1 linear '//tDef//' 1mo'
        write(lunit,*) 'vars ',TotalVars
        do nvar=1,atmosNv
            write(lunit,*) atmosVarNames(nvar),atmosNz,'99 ',atmosVarNames(nvar)
        enddo
        do nvar=1,nSpc
            write(lunit,*) chemVarNames(nvar),atmosNz,'99 ',chemVarNames(nvar)
        enddo
        write(lunit,*) 'endvars'
        ierrNumber=releaseUnit(lunit)

        iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'writing grads/ctl file for '//tDef)
        recordLen=outRealSize()*atmosNx*atmosNy
        lunit=getUnit()
        open(unit=lunit,file=trim(out_dir)//trim(binFileName),action='WRITE',status='REPLACE' &
             ,form='UNFORMATTED',access='DIRECT',recl=recordLen)

        irec=1
        do nvar=1,atmosNv
            do nz=1,atmosNz
                write(lunit,rec=irec) atmosValues(iTime,nvar,:,:,nz)
                irec=irec+1
            enddo
        enddo
        do nvar=1,nSpc
            do nz=1,atmosNz
                write(lunit,rec=irec) chemvalues(nvar,:,:,nz)
            enddo
        enddo

    end subroutine writeGradsCtlFile 

end module filesMod 