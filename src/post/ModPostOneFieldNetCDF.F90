Module ModPostOneFieldNetcdf
#ifdef cdf
  use ModPostGridNetCDF, only: &
      netCdfPostField2D, &
      netCdfPostField3D

  contains

  subroutine writeNetCdf2D(fieldName,nLon,nLat,OutputArray)
  
      integer, intent(in) :: nlon,nlat
      character(len=*), intent(in) :: fieldName
      real, intent(in) :: OutputArray(nlon, nlat)
      
      call netCdfPostField2D(fieldName,nLon,nLat,OutputArray)
  
  end subroutine writeNetCdf2D
  
  subroutine writeNetCdf3D(fieldName,nLon,nLat,iLev,OutputArray)
    integer, intent(in) :: nlon,nlat
    character(len=*), intent(in) :: fieldName
    real, intent(in) :: OutputArray(nlon, nlat)
  
    call netCdfPostField3D(fieldName,nLon,nLat,iLev,OutputArray)
  
  end subroutine writeNetCdf3D
#endif
end Module ModPostOneFieldNetcdf