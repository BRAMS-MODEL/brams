Module ModPostgridNetCDF
#ifdef cdf
  use ModNamelistFile, only: namelistFile

  use ModBramsGrid, only: BramsGrid

  use ModPostUtils, only: UpperCase

  use ModPostTypes

  use mem_grid, only: time

  use dump

  !For netCDF
  include "constants.f90"
  character(len=256) :: netCDFFileName
  logical :: netCDFFirstTime
  integer :: ncid,LatDimID,LonDimID,LevDimID,timDimId,SurDimID,SoiDimID
  character(len=256), allocatable :: netCdfFieldName(:)
  integer, allocatable :: netCdfNLevCode(:)
  integer :: VarDimId(1000)
  character(len=256) :: netCdfFieldDescription(1000)
  character(len=256) :: netCdfFieldUnits(1000)
  real,allocatable :: hoursFrom1900(:)

  public :: FillNetcdfVarControlFile
  public :: netCDFFileName
  public :: netCDFFirstTime
  public :: OpenNetCDFBinaryFile

contains

subroutine OpenNetCDFBinaryFile(oneNamelistFile, onePostGrid, oneBramsGrid, igrid) 
    use netCDF 
    use dump
    type(NamelistFile), pointer :: oneNamelistFile
    type(PostGrid), pointer :: onePostGrid
    type(BramsGrid), pointer :: oneBramsGrid  
    integer, intent(in) :: igrid

    integer :: ierr

    character(len=1) :: c0

    write(c0,fmt='(I1)') igrid

    !if(.not. netCDFFirstTime) return
    if (oneBramsGrid%mchnum /= oneBramsGrid%master_num) return

    CALL makefnam(netCDFFileName, trim(oneNamelistFile%gprefix), time, &
         oneNamelistFile%iyear1, oneNamelistFile%imonth1,  &
         oneNamelistFile%idate1, oneNamelistFile%itime1*100, &
          'A', 'g'//c0, 'nc ')

    iErrNumber = nf90_create(path = trim(netCDFFileName), cmode = nf90_write, ncid = ncid)
    !if (iErrNumber /= nf90_noerr) ierr=dumpMessage(c_tty,c_yes,'','',c_fatal,'Error in create file '//trim(netCDFFileName)//' NetCDF err:',iErrNumber)

end subroutine OpenNetCDFBinaryFile

subroutine FillNetcdfVarControlFile(oneNamelistFile, onePostGrid, oneBramsGrid)
    use mem_grid, only: oneGlobalGridData, nnxp, nnyp, &
                        iyear1,imonth1,idate1,ihour1, &
                        timmax, timeunit,npatch
    use io_params, only: frqanl
    use netcdf
    use dump, only: dumpMessage
    use ModDateUtils

    type(NamelistFile), pointer :: oneNamelistFile
    type(PostGrid), pointer :: onePostGrid
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostVarType) :: one_post_variable
    character(len=*), parameter :: h='**(FillNetcdfVarControlFile)**'

    include "constants.f90"
    include "ranks.h"

    integer :: i,ndims,nvars,cnt
    character(len=256) :: name,varname(30),atName(300),atValue(300)
    character(len = 16) :: varNameUpper,vName

    integer, dimension(300) :: lenDim,nat
    integer :: vdim,levs,an,ntimes,id,ivp
    real :: slayer(oneBramsGrid%nzg)
    real(kind=r8) :: seconds
    real :: lon(nnxp(1)-1),lat(nnyp(1)-1)
    character(len=8) :: cVar

    !if(.not. netCDFFirstTime) return
    if (oneBramsGrid%mchnum /= oneBramsGrid%master_num) return


    iErrNumber = nf90_def_dim(ncid, "Longitude", nnxp(1)-1, LonDimID)
    if (iErrNumber /= nf90_noerr) print *,'err'
    iErrNumber = nf90_def_dim(ncid, "Latitude", nnyp(1)-1, LatDimID)
    if (iErrNumber /= nf90_noerr) print *,'err'
    iErrNumber = nf90_def_dim(ncid, "Level",onePostGrid%nVert, LevDimID)
    if (iErrNumber /= nf90_noerr) print *,'err'
    ntimes=1!timmax/frqanl
    iErrNumber = nf90_def_dim(ncid, "Time", ntimes, TimDimID)
    if (iErrNumber /= nf90_noerr) print *,'err'
    iErrNumber = nf90_def_dim(ncid, "soilLevel", oneBramsGrid%nzg, SoiDimID)
    if (iErrNumber /= nf90_noerr) print *,'err'

    ntimes=1
    if (.not. allocated(hoursFrom1900)) allocate(hoursFrom1900(nTimes))

    ! Define the variable Longitude
    iErrNumber = nf90_def_var(ncid, "Longitude", nf90_float, &
    (/LonDimId/), LonDimId)

    ! Define the variable Latitude
    iErrNumber = nf90_def_var(ncid, "Latitude", nf90_float, &
    (/LatDimId/), LatDimId)

    ! Define the variable Level
    iErrNumber = nf90_def_var(ncid, "Level", nf90_float, &
    (/LevDimId/), LevDimId)

    ! Define the variable Level
    iErrNumber = nf90_def_var(ncid, "Time", nf90_float, &
    (/TimDimId/), TimDimId)

    ! Define the variable Level
    iErrNumber = nf90_def_var(ncid, "SoilLevel", nf90_float, &
    (/SoiDimId/), SoiDimId)

    iErrNumber = nf90_redef(ncid)

    iErrNumber = nf90_put_att(ncid, LonDimId, "units", "Degrees_east")
    iErrNumber = nf90_put_att(ncid, LonDimId, "long name", "longitude")
    iErrNumber = nf90_put_att(ncid, LatDimId, "units", "Degrees_north")
    iErrNumber = nf90_put_att(ncid, LatDimId, "long name", "latitude")
    iErrNumber = nf90_put_att(ncid, LevDimId, "units", "mBar")
    iErrNumber = nf90_put_att(ncid, LevDimId, "long name", "pressure level")
    iErrNumber = nf90_put_att(ncid, TimDimId, "units", "hours since 1900-01-01 00:00:00")
    iErrNumber = nf90_put_att(ncid, TimDimId, "long name", "time")   
    iErrNumber = nf90_put_att(ncid, TimDimId, "calendar", "gregorian")
    iErrNumber = nf90_put_att(ncid, soiDimId, "units", "m")
    iErrNumber = nf90_put_att(ncid, soiDimId, "long name", "soil level")
    
    cnt=0
    do ivp = 1, oneNamelistFile%nvp
      vName=oneNamelistFile%vp(ivp)
      varNameUpper = trim(UpperCase(vName))
      one_post_variable = getPostVarible(varNameUpper)
      if (len(trim(one_post_variable%fieldName)) .eq. 0) then
         write(logUnit, "(a)") "**(OnePostField)** Post field " // varName // " does not exists in list of variables."
         write(logUnit, "(a)") "Model will continue ..."
      else
         select case (one_post_variable%ivar_type)
         case (2) !Surface var
            cnt=cnt+1
            netCdfFieldDescription(cnt)=one_post_variable%fieldDescription
            netCdfFieldUnits(cnt)=one_post_variable%fieldUnits      
            iErrNumber = nf90_def_var(ncid,name=varNameUpper, xtype=nf90_float, &
                  dimids=(/LonDimId,LatDimId,TimDimId/),varid=VarDimId(cnt))
         case (3) !Athmos var
            cnt=cnt+1
            netCdfFieldDescription(cnt)=one_post_variable%fieldDescription
            netCdfFieldUnits(cnt)=one_post_variable%fieldUnits
            iErrNumber = nf90_def_var(ncid,name=varNameUpper, xtype=nf90_float, &
                  dimids=(/LonDimId,LatDimId,LevDimId,TimDimId/),varid=VarDimId(cnt)) 
         case (7) !Surface Soil var
            do i=1,npatch
              cnt=cnt+1
              write(cvar,fmt='(I8)') i 
              vName=trim(varNameUpper)//trim(adjustl(cvar))
              netCdfFieldDescription(cnt)=trim(one_post_variable%fieldDescription)//': patch #'//trim(cvar)
              netCdfFieldUnits(cnt)=one_post_variable%fieldUnits
              iErrNumber = nf90_def_var(ncid,name=vName, xtype=nf90_float, &
              dimids=(/LonDimId,LatDimId,TimDimId/),varid=VarDimId(cnt))
            enddo
         case (8) !
            do i=1,npatch
              cnt=cnt+1
              write(cvar,fmt='(I8)') i 
              vName=trim(varNameUpper)//trim(adjustl(cvar))
              netCdfFieldDescription(cnt)=trim(one_post_variable%fieldDescription)//': patch #'//trim(cvar)
              netCdfFieldUnits(cnt)=one_post_variable%fieldUnits
              iErrNumber = nf90_def_var(ncid,name=vName, xtype=nf90_float, &
                  dimids=(/LonDimId,LatDimId,SoiDimId,TimDimId/),varid=VarDimId(cnt))
            enddo
         end select
         do i = 1, size(all_post_variables)
            if(varNameUpper .eq. all_post_variables(i)%fieldName) then
              all_post_variables(i)%netcdfId=VarDimId(cnt)
               exit
            end if
         end do
      end if
    enddo
 
    do i=1,cnt
    !   print *,'vardim(i)=',VarDimId(i)
       iErrNumber = nf90_put_att(ncid, VarDimId(i), "scale_factor", (/1.0/))
    !   if (iErrNumber /= nf90_noerr) print *,'err 1'
       iErrNumber = nf90_put_att(ncid, VarDimId(i), "add_offset", (/0.0/))
    !   if (iErrNumber /= nf90_noerr) print *,'err 2'
       iErrNumber = nf90_put_att(ncid, VarDimId(i), "_FillValue", (/-32767./))
    !   if (iErrNumber /= nf90_noerr) print *,'err 3'
       iErrNumber = nf90_put_att(ncid, VarDimId(i), "missing_value",(/-32767./))
    !   if (iErrNumber /= nf90_noerr) print *,'err 4'
       iErrNumber = nf90_put_att(ncid, VarDimId(i), "units" &
                 ,trim(netCdfFieldUnits(i)))
       if (iErrNumber /= nf90_noerr) print *,'err 5'
       iErrNumber = nf90_put_att(ncid, VarDimId(i), "long_name" &
                 , trim(netCdfFieldDescription(i)))
       if (iErrNumber /= nf90_noerr) print *,'err 6'
    enddo

    iErrNumber = nf90_enddef(ncid)

    !Fiiling lons, lats and pressure levels 
!print *,'LFR->Fiiling lons, lats and pressure levels '
    do i=1,nnyp(1)-1
      lat(i)=oneGlobalGridData(1)%global_glat(1,nnyp(1)-i+1)
    enddo
!    print *,lat
    iErrNumber = nf90_put_var(ncid, LatDimId, lat)

    !Making glon from 0 to 360 east direction
    do i=1,nnxp(1)-1
      if(oneGlobalGridData(1)%global_glon(i,1)<0) &
        lon(i)=360+oneGlobalGridData(1)%global_glon(i,1)
    enddo
!    print *,lon
    iErrNumber = nf90_put_var(ncid, LonDimId, lon)
!    print *,onePostGrid%vertScaleValues
    iErrNumber = nf90_put_var(ncid, LevDimId, onePostGrid%vertScaleValues)

    !Filling Soil layers - Just putting a number to represent it
    do i=1,oneBramsGrid%nzg
      slayer(i)=real(i)
    enddo
    iErrNumber = nf90_put_var(ncid, SoiDimId, slayer)
    
    !Filling Times
    !First convert the initial date to seconds
    call date_abs_secs2(iyear1,imonth1,idate1,ihour1,seconds)
    !Now put increments of frqanl in hours inside array
    do i=1,ntimes
      hoursFrom1900(i)=(((frqanl*(i-1))+seconds)/3600.)-24
    enddo
    !Finally put the array in TimDimId
    iErrNumber = nf90_put_var(ncid, TimDimId, hoursFrom1900)
    !

    netCDFFirstTime=.false.

  end subroutine FillNetcdfVarControlFile

  subroutine netCdfPostField2D(fieldName,nLon,nLat,OutputArray)

    use netcdf
    use dump, only: dumpMessage
    use mem_grid, only: time
    use io_params, only: frqanl

    include "constants.f90"
    integer, intent(in) :: nlon,nlat
    character(len=*), intent(in) :: fieldName
    real, intent(in) :: OutputArray(nlon, nlat)
    integer :: id,nt
    character(len=256) :: varName(500),name
    real :: varArray(nlon, nlat)
    type(PostVarType) :: one_post_variable

    one_post_variable = getPostVarible(fieldName)
    call invertLats(OutputArray,varArray,nlon,nlat)

    iErrNumber = nf90_put_var(ncid,one_post_variable%netcdfId, varArray,start = (/1, 1, 1/) )

  end subroutine netCdfPostField2D

  subroutine netCdfPostField3D(fieldName,nLon,nLat,ilev,OutputArray)

    use netcdf
    use dump, only: dumpMessage
    use mem_grid, only: time
    use io_params, only: frqanl

    include "constants.f90"
    integer, intent(in) :: nlon,nlat,iLev
    character(len=*), intent(in) :: fieldName
    real, intent(in) :: OutputArray(nlon,nlat)
    integer :: id,nt
    character(len=256) :: varName(500),name
    real :: varArray(nlon,nlat)
       type(PostVarType) :: one_post_variable

    one_post_variable = getPostVarible(fieldName)

    call invertLats(OutputArray,varArray,nlon,nlat)

    iErrNumber = nf90_put_var(ncid,one_post_variable%netcdfId, varArray,start=(/1,1,iLev,1/))!, start = (/ 4, 3, 2 /) )

  end subroutine netCdfPostField3D

  subroutine invertLats(iArray,oArray,nLon,nLat)

    integer,intent(in) :: nLon,nLat
    real,intent(in) :: iArray(nLon,nLat)
    real,intent(out) :: oArray(nLon,nLat)

    integer :: i,j

    do i=1,nLon
      do j=1,nLat
        oArray(i,j)=iArray(i,nLat-j+1)
      enddo
    enddo

  end subroutine invertLats
#endif
end Module ModPostgridNetCDF

