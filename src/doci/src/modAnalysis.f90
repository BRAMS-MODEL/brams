module modAnalysis
    use dump, only: dumpMessage
    implicit none

    character(len=*),parameter :: srcFile='modAnalysis.f90'

contains

subroutine analysisNcep(innpr,id,outFolder)

    use wgrib2api
    use ModDateUtils, only: date_add_to
    use modUtils, only: str2int

    include "constants.h"
    character(len=*), intent(in) :: innpr,outFolder
    integer, intent(in) :: id

    integer :: nprz_grib2
    character(len=*), parameter :: wind_u_varname='UGRD'   
    character(len=*), parameter :: wind_v_varname='VGRD' 
    character(len=*), parameter :: temperature_varname='TMP'
    character(len=*), parameter :: geo_varname='HGT'     
    character(len=*), parameter :: ur_varname='RH' 
    character,parameter,dimension(5) :: nameV=(/"U","V","T","Z","H"/) 
    character(len=10) :: varName(5)
    real, allocatable :: var(:,:),lat(:,:),lon(:,:),dados(:,:,:,:)
    character (len=200) :: metadata
    character (len=300) :: grid_info
    character (len=99) :: invline
    character (len=30), allocatable :: slevs(:)
    character (len=5), allocatable :: clevs(:)
    integer :: nx
    integer :: ny

    integer :: i,j,k,nvar
    integer :: iyyl,imml,iddl,ihhl
    integer :: stat,irec,recordLen,thisHour
    integer :: iyy,imm,idd,ihh
    character(len=12) :: foutname
    character(len=15) :: gradsData
    real :: dx,dy

    logical :: there

    inquire(file=trim(innpr),exist=there)
    if(.not. there) then
        write (*,fmt='(A)') 'Arquivo '//trim(innpr)//' nao encontrado. Verifique!'
        stop 'Error!'
    endif
    
    write(*,fmt='(I2.2,1X,A)') id," analisando "//trim(innpr)
    varName(1)=':TMP:'

    !Making inventory using a wgrib2 function grb2_mk_inv
    iErrNumber = grb2_mk_inv(innpr(1:len_trim(innpr)),innpr(1:len_trim(innpr))//'.inv')
    
!   get number of first var levels using grb2_inq from wgrib2 lib
    nprz_grib2 = grb2_inq(innpr(1:len_trim(innpr)),innpr(1:len_trim(innpr))//'.inv' &
                          ,trim(varname(1)),' mb:')

    allocate(slevs(nprz_grib2))
    allocate(clevs(nprz_grib2))
    print *,'Nz= ',nprz_grib2

    do i = 1,nprz_grib2
      ! get pressure leves and levels labels using grb2_inq from wgrib2 lib
      iErrNumber=grb2_inq(trim(innpr),trim(innpr)//'.inv',trim(varname(1)),' mb:' &
                 ,sequential=i-1,desc=metadata)
      if (iErrNumber.ne.1) stop 2

      j = index(metadata,trim(varname(1))) + len(trim(varname(1)))
      k = index(metadata," mb:") + len(" mb:")-1

      slevs(i) = metadata(j-1:k)
      clevs(i) = metadata(j:k-3)

    enddo
    
    !Getting information in GRIB2 file. Using the first variable (U)
    iErrNumber = grb2_inq(trim(innpr),trim(innpr)//'.inv',trim(varName(1)) &
                 , data2=var &
                 , lat=lat, lon=lon, grid_desc=grid_info, desc=invline)

    !Nlats e Nlons
    nx=size(var,1)
    ny=size(var,2)
    ! Increments lat and lon
    dy=(lat(1,ny)-lat(1,1))/(ny-1)
    dx=(lon(nx,1)-lon(1,1))/(nx-1)

    allocate(dados(5,nx,ny,nprz_grib2-1))

    !print *,trim(invline)
    !print *,'------------'
    !print *,trim(grid_info)

    !convert date string got from grib2 file to integer
    call str2int(invline(3:6)  ,iyyl,stat)
    call str2int(invline(7:8)  ,imml,stat)
    call str2int(invline(9:10) ,iddl,stat)
    call str2int(invline(11:12),ihhl,stat)
    call str2int(innpr(len_trim(innpr)-2:len_trim(innpr)),thisHour,stat)
    
    call date_add_to(iyyl,imml,iddl,ihhl,real(thisHour*3600),'s' &
        ,iyy,imm,idd,ihh)
    ihh=ihh/10000

    write(foutname,fmt='("IC",I4.4,I2.2,I2.2,I2.2)') iyy,imm,idd,ihh

    open(unit=33,file=trim(outFolder)//foutname//'.dat',&
        action='WRITE',status='REPLACE',form='FORMATTED')

    write(*,fmt='(I2.2,A)') id,' - Escrevendo '//trim(innpr)//'.dat'

    do nvar=1,5
        write (33,fmt='(A20,1X,A8,1X,A8)') 'Var&Level','Min','Max'
        do i=2,nprz_grib2

            !Compose name of var to search inside grib2
            varName(1)=':'//trim(wind_u_varname)     //trim(slevs(i))    
            varName(2)=':'//trim(wind_v_varname)     //trim(slevs(i)) 
            varName(3)=':'//trim(temperature_varname)//trim(slevs(i)) 
            varName(4)=':'//trim(geo_varname)        //trim(slevs(i))   
            varName(5)=':'//trim(ur_varname)         //trim(slevs(i)) 

            iErrNumber = grb2_inq(trim(innpr),trim(innpr)//'.inv',trim(varName(nVar)) &
               ,data2=var, lat=lat, lon=lon)
            
            if(nvar==5) then
              dados(nvar,:,:,i-1)=var/100.0
            else
              dados(nvar,:,:,i-1)=var
            endif

            write (33,fmt='(A20,1X,F8.2,1X,F8.2)') trim(varName(nVar)), &
                  minval(var),maxval(var)
           
            select case (nvar)
                case (1,2)
                    if(minval(dados(nvar,:,:,:))<(-100) .or. maxval(var)>100) &
                    write(33,fmt='(A)') trim(varName(nVar)) &
                    //' - *** Existem valores acima ou abaixo do limite! (-100m/s,+100m/s) ***'
                case (3)
                    if(minval(var)<180 .or. maxval(var)>333) &
                    write(33,fmt='(A)') trim(varName(nVar)) &
                    //' - *** Existem valores acima ou abaixo do limite! (+180K,+333K) ***'
                case (4)
                    if(minval(var)<(-500) .or. maxval(var)>80000) &
                    write(33,fmt='(A)') trim(varName(nVar))&
                    //' - *** Existem valores acima ou abaixo do limite! (-2000gpm,+80000gpm) ***'
                case (5)
                    if(minval(var)<0 .or. maxval(var)>100) &
                    write(33,fmt='(A)') trim(varName(nVar))&
                    //' - *** Existem valores acima ou abaixo do limite! (0%,100%) ***'
            end select

        enddo
 
    enddo

    close(33)

    !Closing grib2 file to save memory
    iErrNumber = grb2_free_file(trim(innpr))

    write(foutname,fmt='("IC",I4.4,I2.2,I2.2,I2.2)') iyy,imm,idd,ihh
    write(gradsData,fmt='(I2.2,":00z",I2.2,A3,I4.4)') ihh,idd,month_name(imm),iyy

    write(*,fmt='(I2.2," - ",5(I4.4,1X),A)') id,iyyl,imml,iddl,ihhl,thisHour,' Escrevendo '//foutname//'...'
    recordLen=4*nx*ny
    open(unit=33,file=trim(outFolder)//foutname//'.gra',&
        action='WRITE',status='REPLACE',form='UNFORMATTED',access='DIRECT', &
        recl=recordLen)
    
    !# writing grads binary and fill variables
    irec=1
    do nvar=1,5
      do k=nprz_grib2,2,-1
!      do k=1,nprz_grib2-1
          write (33,rec=irec) dados(nvar,1:nx,1:ny,k-1)
          irec=irec+1
      enddo
    enddo
    
    close(33) 
    
    open(unit=33,file=trim(outFolder)//foutname//'.ctl' &
         ,action='WRITE',status='replace',form='FORMATTED')
    
    !writing the name of grads file
    write(33,*) 'dset ^'//foutname//'.gra'
    !writing others infos to ctl
    write(33,*) 'undef -0.9990000E+34'
    write(33,*) 'title GRIB2 NCEP GFS'
    write(33,*) 'xdef ',nx,' linear ',lon(1,1),dx
    write(33,*) 'ydef ',ny,' linear ',lat(1,1),dy
    write(33,*) 'zdef ',nprz_grib2-1,'levels ',(trim(clevs(k))//' ',k=nprz_grib2,2,-1)
    write(33,*) 'tdef 1 linear ',gradsData,' 1mo'
    write(33,*) 'vars 5'
    write(33,*) 'U ',nprz_grib2-1,'99 U wind'
    write(33,*) 'V ',nprz_grib2-1,'99 V Wind'
    write(33,*) 'T ',nprz_grib2-1,'99 Temp'
    write(33,*) 'Z ',nprz_grib2-1,'99 Geopot'
    write(33,*) 'R ',nprz_grib2-1,'99 Rel Umid'
    write(33,*) 'endvars'

    close(33)

end subroutine analysisNcep





subroutine analysisNasa(innpr,id,outFolder)
    use netcdf
    use modUtils, only: monthName

    include "netcdf.inc"   
    include "constants.h"
    character(len=*),parameter :: header='**(analysisNasa)**'
    integer, parameter :: validTime=1
    character(len=*), intent(in) :: innpr,outFolder
    integer, intent(in) :: id
    integer :: ncid,ndims,nvars
    integer, allocatable :: nat(:),varDim(:),lenDim(:)
    character(len=32),allocatable :: varName(:)
    integer :: nprz_grib2
    character(len=*), parameter :: wind_u_varname='U'   
    character(len=*), parameter :: wind_v_varname='V' 
    character(len=*), parameter :: temperature_varname='T'
    character(len=*), parameter :: geo_varname='H'     
    character(len=*), parameter :: ur_varname='RH' 
    character(len=32) :: name,atName
    character(len=33) :: timeLong
    integer :: nx,ny,nz,nt,iFound
    integer :: levelVarN,lonVarN,latVarN,timeVarN
    integer :: uVarN,vVarN,tVarN,qVarN,zVarN
    real, allocatable :: PressLevs(:),llat(:),llon(:)
    real :: scaleFactor(32),addOfset(32)
    integer :: beginDate,beginTime,xtype
    integer :: an,iLen,nvar,k,recordLen
    real :: missingValue
    real :: xvalue,dx,dy
    real, allocatable :: dados(:,:,:,:)
    character(len=256) :: fOutName
    character(len=15) :: gradsData

    integer :: i,irec

    iErrNumber = nf90_open(path = trim(innpr), mode = nf90_nowrite, ncid = ncid)
    if (iErrNumber /= nf90_noerr) iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion &
              ,c_fatal,trim(innpr)//' not found. Please, check it!')

    ! get info about netCDF file 
    iErrNumber=nf90_inquire(ncid, ndims, nvars)
    !print *,'File open. Total of variables are ',nvars
    allocate(varName(nvars),varDim(nvars),nat(nvars),lenDim(nvars))

    !Get dimensions in levels, lons, lats and times
    nz=0
    nx=0
    ny=0
    do i=1,nvars
      iErrNumber = nf90_Inquire_Dimension(ncid, i, name, lenDim(i))

      select case (trim(name))
      case ('level') 
        nz=lenDim(i)
        levelVarN=i
      case ('lev')
        nz=lenDim(i)
        levelVarN=i    
      case ('longitude')
        nx=lenDIm(i)
        lonVarN=i
      case ('lon')
        nx=lenDIm(i)
        lonVarN=i
      case ('latitude')
        ny=lenDIm(i)
        latVarN=i
      case ('lat')
        ny=lenDIm(i)
        latVarN=i
      case ('time')
        nt=lenDIm(i)
        timeVarN=i
      end select
    enddo

    !If one of variables not found is error
    if(nz==0 .or. nx==0 .or. ny==0) then
      print *,'Vars expected: level, longitude, latitude and time.'
      iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion &
                ,c_fatal,'Var expected for dimension not found in netCDF file. Please check it!')
    endif

    write(*,fmt='("number of lats points: ",I4)') ny
    write(*,fmt='("number of lons points: ",I4)') nx
    write(*,fmt='("number of levs       : ",I4)') nz
    write(*,fmt='("number of times      : ",I4)') nt 

    !Allocate all vars 
    allocate(PressLevs(nz))
    allocate(llat(ny),llon(nx))
    allocate(dados(5,nx,ny,nz))

    !do i=1,5
    !  allocate(ncVar(i)%valuesRead(nxGrib,nyGrib,nprz_grib2))
    !enddo

    !Get pressure levels
    iErrNumber=nf90_get_var(ncid, levelVarN, PressLevs)
    !do i=1,nz
    !    print *,'PressLevs: ', i,int(pressLevs(i))
    !enddo

    !Get names of vars in netCDF file and total of atributes
    do i=1,nvars
      iErrNumber=nf90_Inquire_Variable(ncid, i, name=varName(i) &
                 ,ndims=varDim(i), nAtts=nat(i))
      !print *,i,trim(varName(i)),varDim(i),nat(i)
    enddo

    !checking if required varName is present
    iFound=0
    do i=1,nvars
      if(trim(varName(i))==trim(WIND_U_VARNAME)) then
        iFound=iFound+1
        uVarN=i 
      elseif(trim(varName(i))==trim(WIND_V_VARNAME)) then
        iFound=iFound+1
        vVarN=i
      elseif(trim(varName(i))==trim(TEMPERATURE_VARNAME)) then
        iFound=iFound+1
        tVarN=i 
      elseif(trim(varName(i))==trim(GEO_VARNAME)) then
        iFound=iFound+1
        zVarN=i 
      elseif(trim(varName(i))==trim(UR_VARNAME)) then
        iFound=iFound+1
        qVarN=i 
      endif
    enddo
    if(iFound/=5) then !Not found all needs
      print *,'Available names of var in netCDF file:'
      do i=1,nvars
        print *,trim(varName(i))
      enddo
      print *,'Expected name os vars from RAMSIN:'
      print *,trim(WIND_U_VARNAME)
      print *,trim(WIND_V_VARNAME)
      print *,trim(TEMPERATURE_VARNAME)
      print *,trim(GEO_VARNAME)
      print *,trim(UR_VARNAME)
      iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion &
                ,c_fatal,'Incompatible var names. See the list above and RAMSIN vars!')
    endif


    missingValue=1.0E+20
    !Getting variables atributes
    do i=1,nvars
      do an=1,nat(i)
        iErrNumber=nf90_inq_attname(ncid, i, an, atName)
        iErrNumber=nf90_Inquire_Attribute(ncid, i, trim(atName), xtype=xtype, len=iLen)
        if(trim(atName)=='scale_factor') then
          iErrNumber = nf90_get_att(ncid, i, atName, xValue)
          scaleFactor(i)=xvalue
        elseif(trim(atName)=='add_offset') then
          iErrNumber = nf90_get_att(ncid, i, atName, xValue)
          addOfset(i)=xvalue
        elseif(trim(atName)=='begin_date') then
          iErrNumber = nf90_get_att(ncid, i, atName, xValue)
          beginDate=xValue
        elseif(trim(atName)=='begin_time') then
          iErrNumber = nf90_get_att(ncid, i, atName, xValue)
          beginTime=xValue
        elseif(trim(atName)=='units' .and. i==timeVarN) then
          iErrNumber = nf90_get_att(ncid, i, atName, timeLong)
        elseif(trim(atName)=='missing_value') then
          iErrNumber = nf90_get_att(ncid, i, atName, missingValue)
        endif
      enddo
    end do

    !get lon and lat info
    iErrNumber = nf90_get_var(ncid, latVarN, llat)
    iErrNumber = nf90_get_var(ncid, lonVarN, llon)
    dy=llat(2)-llat(1)
    dx=llon(2)-llon(1)

    !print *,'dy=',dy
    !print *,'dx=',dx

    !Read all upper air values for u,v,t,z and r 
    do i=1,5
        select case (i)
        case(1)
          iErrNumber = nf90_get_var(ncid, uVarN, dados(i,:,:,:) &
                       ,start = (/ 1, 1, 1,validTime /))
          call adjustMissingValues(i,nx,ny,nz,missingValue &
                                     ,pressLevs,dados(i,:,:,:))
          dados(i,:,:,:)=dados(i,:,:,:)*scaleFactor(uVarN) &
          +addOfset(uVarN)
        case(2)
          iErrNumber = nf90_get_var(ncid, vVarN, dados(i,:,:,:) &
                       ,start = (/ 1, 1, 1,validTime /))
          call adjustMissingValues(i,nx,ny,nz,missingValue &
                               ,pressLevs,dados(i,:,:,:))
          dados(i,:,:,:)=dados(i,:,:,:)*scaleFactor(vVarN) &
          +addOfset(vVarN)
        case(3)
          iErrNumber = nf90_get_var(ncid, tVarN, dados(i,:,:,:) &
                       ,start = (/ 1, 1, 1,validTime /))
          call adjustMissingValues(i,nx,ny,nz,missingValue &
                                     ,pressLevs,dados(i,:,:,:))
          dados(i,:,:,:)=dados(i,:,:,:)*scaleFactor(tVarN) &
          +addOfset(tVarN)
        case(4)
          iErrNumber = nf90_get_var(ncid, zVarN, dados(i,:,:,:) &
                       ,start = (/ 1, 1, 1,validTime /))
          call adjustMissingValues(i,nx,ny,nz,missingValue &
                               ,pressLevs,dados(i,:,:,:))
          dados(i,:,:,:)=dados(i,:,:,:)*scaleFactor(zVarN) &
          +addOfset(zVarN)
        case(5)
          iErrNumber = nf90_get_var(ncid, qVarN, dados(i,:,:,:) &
                       ,start = (/ 1, 1, 1,validTime /))
          call adjustMissingValues(i,nx,ny,nz,missingValue &
                               ,pressLevs,dados(i,:,:,:))
          dados(i,:,:,:)=dados(i,:,:,:)*scaleFactor(qVarN) &
          +addOfset(qVarN)
        end select
    enddo

    foutName=innpr(len_trim(innpr)-20:len_trim(innpr)-13)
    fOutName=trim(fOutName)//innpr(len_trim(innpr)-11:len_trim(innpr)-10)
    fOutName="IC"//trim(fOutName)

    open(unit=33,file=trim(outFolder)//trim(foutname)//'.dat',&
         action='WRITE',status='REPLACE',form='FORMATTED')

    write(*,fmt='(I2.2,A)') id,' - Escrevendo '&
          //trim(outFolder)//trim(foutname)//'.dat'

    do nvar=1,5
        write (33,fmt='(A20,1X,A8,1X,A8)') 'Var&Level','Min','Max'
        do i=1,nz
            write(name,fmt='(":",A,1X,I4,":")') trim(varName(nVar)) &
                  ,int(pressLevs(i))
            select case (nvar)
                case (1,2)
                    if(minval(dados(nvar,:,:,:))<(-100) .or. maxval(dados(nvar,:,:,:))>100) &
                    write(33,fmt='(A)') trim(name) &
                    //' - *** Existem valores acima ou abaixo do limite! (-100m/s,+100m/s) ***'
                case (3)
                    if(minval(dados(nvar,:,:,:))<180 .or. maxval(dados(nvar,:,:,:))>333) &
                    write(33,fmt='(A)') trim(name) &
                    //' - *** Existem valores acima ou abaixo do limite! (+180K,+333K) ***'
                case (4)
                    if(minval(dados(nvar,:,:,:))<(-500) .or. maxval(dados(nvar,:,:,:))>80000) &
                    write(33,fmt='(A)') trim(name)&
                    //' - *** Existem valores acima ou abaixo do limite! (-2000gpm,+80000gpm) ***'
                case (5)
                    if(minval(dados(nvar,:,:,:))<0 .or. maxval(dados(nvar,:,:,:))>100) &
                    write(33,fmt='(A)') trim(name)&
                    //' - *** Existem valores acima ou abaixo do limite! (0%,100%) ***'
            end select
        enddo
    enddo
    close(33)

    recordLen=4*nx*ny
    open(unit=33,file=trim(outFolder)//trim(foutname)//'.gra',&
        action='WRITE',status='REPLACE',form='UNFORMATTED',access='DIRECT', &
        recl=recordLen)
    
    !# writing grads binary and fill variables
    irec=1
    do nvar=1,5
      do k=1,nz
          write (33,rec=irec) dados(nvar,1:nx,1:ny,k)
          irec=irec+1
      enddo
    enddo
    
    close(33) 

    gradsData=innpr(len_trim(innpr)-11:len_trim(innpr)-10)//":00z"
    gradsData=trim(gradsData)//innpr(len_trim(innpr)-14:len_trim(innpr)-13)
    gradsData=trim(gradsData)//monthName(innpr(len_trim(innpr)-16:len_trim(innpr)-14))
    gradsData=trim(gradsData)//innpr(len_trim(innpr)-20:len_trim(innpr)-15)

    open(unit=33,file=trim(outFolder)//trim(foutname)//'.ctl' &
         ,action='WRITE',status='replace',form='FORMATTED')
    
    !writing the name of grads file
    write(33,*) 'dset ^'//trim(foutname)//'.gra'
    !writing others infos to ctl
    write(33,*) 'undef -0.9990000E+34'
    write(33,*) 'title GRIB2 NASA GEOS'
    write(33,*) 'xdef ',nx,' linear ',llon(1),dx
    write(33,*) 'ydef ',ny,' linear ',llat(1),dy
    write(33,*) 'zdef ',nz,'levels ',(int(pressLevs(k)),k=1,nz)
    write(33,*) 'tdef 1 linear ',gradsData,' 1mo'
    write(33,*) 'vars 5'
    write(33,*) 'U ',nz,'99 U wind'
    write(33,*) 'V ',nz,'99 V Wind'
    write(33,*) 'T ',nz,'99 Temp'
    write(33,*) 'Z ',nz,'99 Geopot'
    write(33,*) 'R ',nz,'99 Rel Umid'
    write(33,*) 'endvars'
    close(33)
    
end subroutine analysisNasa

subroutine adjustMissingValues(nVar,nxGrib,nyGrib,nprz_grib2,missingValue &
                                 ,levpr_grib2,valuesRead)
    !# Subroutine to check and fix variable from input
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: This subroutine looks for missing values inside each input variable 
    !# u,v,t,z and r and fix the problem. 
    !#
    !# If variable is T (Temperature) the subroutine does the replacement using an equation:
    !#
    !# \[$T_z=\frac{{P_{z+1}-P_z}}{\rho.C_p}+T_{z+1}$\]
    !#
    !# where:
    !# 
    !# \[$\rho=\frac{P_z}{R.T_{z+1}}$\]
    !#
    !# otherwise it just copy the value of variable of level above for the current 
    !# Z position
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author**: Luiz Flavio Rodrigues **&#9993;**<mailto:luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 2020-02-19
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !# No active bugs reported now
    !# @endbug
    !#
    !# @todo
    !#  &#9744; <br/>
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !# @endwarning
    !#
    !#--- ----------------------------------------------------------------------------------------
    !
    include "constants.h"
    integer, parameter :: p_tVar=3

    !# Type of input data
    integer, intent(in) :: nVar
    !# Number of incoming variable
    integer, intent(in) :: nxGrib
    !# Size of longitudes
    integer, intent(in) :: nyGrib
    !# Size of latitudes
    integer, intent(in) :: nprz_grib2
    !# Size of levels
    real   , intent(in) :: levpr_grib2(nprz_grib2)
    !# Value of level pressure
    real   , intent(in) :: missingValue
    !# Value for data missing
    real   , intent(inout) :: valuesRead(nxGrib,nyGrib,nprz_grib2)
    !# Array with variable

    integer :: i,j,k
    !# Just counters
    real :: rho
    !# $\rho=\frac{P_z}{R.T_{z+1}}$
    real :: rhoCpInv
    !# $rhoCpInv=\frac{1.}{\rho.C_p}$
    real :: diffPress
    !# $diffPress=P_{k}-P_{k+1}$

    do k=nprz_grib2-1,1,-1 !Cycle through all levels from highest to surface
      do j=1,nyGrib !Cycle through all lats and lons
        do i=1,nxGrib
          if(valuesRead(i,j,k)==missingValue) then !If var value is missing     
            if(nvar==p_tVar) then                   !If var is Temperature
              !Equation
              rho=levpr_grib2(k)/(c_rgas*valuesRead(i,j,k+1))
              rhoCpInv=1./(rho*c_cp)
              diffPress=levpr_grib2(k)-levpr_grib2(k+1)
              valuesRead(i,j,k)=valuesRead(i,j,k+1)+diffPress*rhoCpInv
            else                                    !If var is U,V,Z or R
              valuesRead(i,j,k)=valuesRead(i,j,k+1) !Just copy the value from level above
            endif
          endif
        enddo
      enddo
    enddo

  end subroutine adjustMissingValues

end module modAnalysis