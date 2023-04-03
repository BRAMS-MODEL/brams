module smnasa
	!para pegar o arquivo de SM da NASA instale o OpenGrads e use o comando:
	!
	!./lats4d.sh -i https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/fcast/tavg1_2d_lnd_Nx/tavg1_2d_lnd_Nx.20200914_00 
	! -time 00:30z14sep2020 00:30z14sep2020 -ftype sdf -o GEOS.SM.20200914_0030 -format netcdf4 -gzip 2 -vars gwetroot 
	! gwetprof prectot -lon -85 -30 -lat -60 20 -v
	!
	!substituindo a data para a data de interesse

	contains

	subroutine sm_nasa(smFile,outFolder,year,month,day,hour)
		!# _ 
		!#
		!# @note
		!# 
		!#
		!# **Brief**: _ 
		!#
		!# **Documentation**: <LF>
		!#
		!# **Author**: Luiz Flavio **&#9993;**<mailto:lufla@enterprise>
		!#
		!# **Date**: 2020-04-29
		!# @endnote
		!#
		!# @changes
		!#**Changelogs:**
		!# +
		!# @endchanges
		!# @bug
		!# **Open bugs:**
		!#
		!# @endbug
		!#
		!# @todo
		!# **Todo list:**n	!#
		!# @endtodo
		!#
		!# @warning
		!# Now is under CC-GPL License, please see
		!# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
		!# @endwarning
		!#
		!#--- ----------------------------------------------------------------------------------------
		use dump
		use netcdf
		
		implicit none
	
		character(len=*),parameter :: procedureName="sm_nasa"
		character(len=*),parameter :: srcName="modAnalysis.f90"
		include "constants.f90"
	
		character(len=*), intent(in) :: smFile
		character(len=*), intent(in) :: outFolder
		integer,intent(in) :: year,month,day,hour
		
		integer, parameter :: validTime=1
		!Local:
		integer :: ncid
		integer :: ndims
		integer :: nvars
		integer :: beginDate,beginTime,xtype
		integer, allocatable :: nat(:),lenDim(:),varDim(:)
		character(len=32),allocatable :: varName(:)
		character(len=32) :: name,atName
		real(kind=kind_rb), allocatable :: missingV(:)
		real(kind=kind_rb) :: missingValue
		integer :: nx,ny
		integer :: levelVarN,latVarN,lonVarN,timeVarN
		real, allocatable :: llat(:),llon(:),dados(:,:,:)
		integer :: i,iFound,grN,gpN,pcN,an,ilen,irec,recordLen
		integer :: k,nvar
		real :: dx,dy
		character(len=256) :: fOutName
		character(len=15) :: gradsData
		character(len=64) :: cFmt
		integer :: nLevels
		real, allocatable :: us(:,:,:),levels(:)
		!#_
			
		!Code area:
		iErrNumber = nf90_open(path = trim(smFile), mode = nf90_nowrite, ncid = ncid)
		if (iErrNumber /= nf90_noerr) iErrNumber=dumpMessage(c_tty,c_yes,procedureName,srcName &
				,c_fatal,trim(smFile)//' not found. Please, check it!')
		! get info about netCDF file 
		iErrNumber=nf90_inquire(ncid, ndims, nvars)
		allocate(varName(nvars),varDim(nvars),nat(nvars),lenDim(nvars),missingV(nVars))
		
			!Get dimensions in lons, lats 
		nx=0
		ny=0
		do i=1,nvars
		iErrNumber = nf90_Inquire_Dimension(ncid, i, name, lenDim(i))
	
		select case (trim(name))   
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
		end select
		enddo
	
		!If one of variables not found is error
		if(nx==0 .or. ny==0) then
		print *,'Vars expected: longitude, latitude.'
		iErrNumber=dumpMessage(c_tty,c_yes,procedureName,srcName &
					,c_fatal,'Var expected for dimension not found in netCDF file. Please check it!')
		endif
	
		!Allocate all vars 
		allocate(llat(ny),llon(nx))
		allocate(dados(3,nx,ny))
		
	!Forcando os niveis
		nLevels=8
		allocate(us(nx,ny,nLevels))
		allocate(levels(nLevels))
	
	!Fazendo os niveis iguais aos que o BRAMS usava
		levels=(/-3.25,-2.125,-1.375,-0.75,-0.375,-0.19,-0.09,-0.025/)
	
		do i=1,nvars
		iErrNumber=nf90_Inquire_Variable(ncid, i, name=varName(i) &
					,ndims=varDim(i), nAtts=nat(i))
		!print *,i,trim(varName(i)),varDim(i),nat(i)
		enddo
		
		!checking if required varName is present
		iFound=0
		do i=1,nvars
		if(trim(varName(i))=='gwetroot') then
			iFound=iFound+1
			grN=i 
		elseif(trim(varName(i))=='gwetprof') then
			iFound=iFound+1
			gpN=i
		elseif(trim(varName(i))=='prectot') then
			iFound=iFound+1
			pcN=i 
		endif
		enddo
		!print *,grN,gpN,pcN
		
		!Getting variables atributes
		do i=1,nvars
		do an=1,nat(i)
			iErrNumber=nf90_inq_attname(ncid, i, an, atName)
			iErrNumber=nf90_Inquire_Attribute(ncid, i, trim(atName), xtype=xtype, len=iLen)
			!print *,an,trim(atName)
			if(trim(atName)=='missing_value') then
				iErrNumber = nf90_get_att(ncid, i, atName, missingValue)
				missingV(i)=missingValue
			endif
		enddo
		enddo
		
		!get lon and lat info
		iErrNumber = nf90_get_var(ncid, latVarN, llat)
		iErrNumber = nf90_get_var(ncid, lonVarN, llon)
		dy=llat(2)-llat(1)
		dx=llon(2)-llon(1)
	
		do i=1,3
			select case (i)
			case(1)
			iErrNumber = nf90_get_var(ncid, grN, dados(i,:,:) &
						,start = (/ 1, 1, 1,validTime /))
			!call adjustMissingValues(i,nx,ny,nz,missingValue &
			!                           ,pressLevs,dados(i,:,:,:))
			!dados(i,:,:,:)=dados(i,:,:,:)*scaleFactor(uVarN) &
			!+addOfset(uVarN)
			case(2)
			iErrNumber = nf90_get_var(ncid, gpN, dados(i,:,:) &
						,start = (/ 1, 1, 1,validTime /))
			!call adjustMissingValues(i,nx,ny,nz,missingValue &
			!                     ,pressLevs,dados(i,:,:,:))
			!dados(i,:,:,:)=dados(i,:,:,:)*scaleFactor(vVarN) &
			!+addOfset(vVarN)
			case(3)
			iErrNumber = nf90_get_var(ncid, pcN, dados(i,:,:) &
						,start = (/ 1, 1, 1,validTime /))
			!call adjustMissingValues(i,nx,ny,nz,missingValue &
			!                           ,pressLevs,dados(i,:,:,:))
			!dados(i,:,:,:)=dados(i,:,:,:)*scaleFactor(tVarN) &
			!+addOfset(tVarN)
			end select
		enddo
		
		write(fOutName,fmt='("SM.GEOS.",I4.4,3(I2.2))') year,month,day,hour 
		
		print *,'Escrevendo '//trim(outFolder)//trim(foutname)//'.gra/.ctl'
		
		recordLen=4*nx*ny
		open(unit=33,file=trim(outFolder)//trim(foutname)//'.gra',&
			action='WRITE',status='REPLACE',form='UNFORMATTED',access='DIRECT', &
			recl=recordLen)
		
		!# writing grads binary and fill variables
		irec=1
		
		do k=1,nLevels
			if(levels(k)>=-0.10) then ! SE o nivel for maior ou igual a 10 cm
				us(:,:,k)=dados(grn,1:nx,1:ny)
			else
				us(:,:,k)=dados(gpn,1:nx,1:ny)
			endif
		enddo
		
		do k=1,nLevels
			write (33,rec=irec) us(1:nx,1:ny,k)
			irec=irec+1
		enddo
		write (33,rec=irec) dados(pcN,1:nx,1:ny)
	
		
		close(33) 

		write(gradsData,fmt='(I2.2,":00z",I2.2,A3,I4.4)') hour,day,month_Name(month),year

		open(unit=33,file=trim(outFolder)//trim(foutname)//'.ctl' &
			,action='WRITE',status='replace',form='FORMATTED')
		
		write(cFmt,fmt='(A,I2.2,A)') '(A,I2.2,1X,A,1X',nLevels,'(F8.3,1X))'
		
		!writing the name of grads file
		write(33,*) 'dset ^'//trim(foutname)//'.gra'
		!writing others infos to ctl
		write(33,*) 'undef ',1.0e+15
		write(33,*) 'title GRIB2 NASA GEOS'
		write(33,*) 'xdef ',nx,' linear ',llon(1),dx
		write(33,*) 'ydef ',ny,' linear ',llat(1),dy
		write(33,fmt=cFmt) 'zdef ',nlevels,'levels ',levels(:)
		write(33,*) 'tdef 1 linear ',gradsData,' 1mo'
		write(33,*) 'vars 2'
		write(33,*) 'us ',nlevels,'99 umid do solo  [mm^3/mm^3]'
		write(33,*) 'rr ',0,'99 precipitacao[mm]'
		write(33,*) 'endvars'
		close(33)
		
	
	end subroutine sm_nasa
	
end module smnasa

	
