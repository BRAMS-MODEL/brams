program geraSM
	use smNasa
	
	include "constants.f90"
	character(len=*),parameter :: cmd1='lats4d.sh -i https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/fcast/tavg1_2d_lnd_Nx/tavg1_2d_lnd_Nx.'
	character(len=*),parameter :: cmd2=' -format netcdf4 -gzip 2 -vars gwetroot gwetprof prectot -lon -85 -30 -lat -60 20 -v '
	
	character(len=4) :: arg(5)
	integer :: year,month,day,hour
	character(len=256) :: fileName
	character(len=512) :: command
	character(len=11) :: datac
	character(len=15) :: dataL
	logical :: download
	character(len=32) :: cArg
	
	
	integer :: i

    i=1
    do
		call get_command_argument(i, cArg)
		if (len_trim(cArg) == 0) exit
		!write (*,*) i,trim(cArg)
		arg(i)=trim(cArg)
		i = i+1
		!
		if(trim(cArg)=='-h' .or. trim(cArg)=='--help') then
			call help()
			return
			!
		endif
	end do
	!print *,i
	if(i<5 .or. i>6) then
		call help()
		return
	endif
    
    read(arg(1),fmt='(I4.4)') year
    read(arg(2),fmt='(I4.4)') month
    read(arg(3),fmt='(I4.4)') day
    read(arg(4),fmt='(I4.4)') hour
	
	if(trim(arg(5))=='y' .or. trim(arg(5))=='Y' .or. trim(arg(5))=='s' .or. trim(arg(5))=='S') then
		download=.true.
	else
		download=.false.
	endif

	
    datac=trim(arg(1))//trim(arg(2))//trim(arg(3))//'_00'
    write(dataL,fmt='(I2.2,":30z",I2.2,A3,I4.4)') hour,day,month_name(month),year
    
    write(command,fmt='(A," -time ",A," -ftype sdf -o GEOS.SM.",A,"30 ",A,A)') cmd1//datac,dataL//" "//dataL,datac,cmd2
    
	write(*,*) 'Program to create Grads SM file from Geos5 Soil Moisture'
	write(*,*) ''
	
	if(.not. download) then
		write(*,*) 'Use the command bellow to get Soil Moisture from Nasa Geos5:'
		write(*,*) ''
		write(*,*) command
		write(*,*) ''
	else
		call system(command)
	endif

	fileName='GEOS.SM.'//datac//'30.nc4'
	
	call sm_nasa(trim(fileName),'./',year,month,day,hour)


end program geraSM	

subroutine help()

	write(*,*) ''
	write(*,*) ' geraSM_1.0 YYYY MM DD HH D, where: '
	write(*,*) ''
	write(*,*) '  YYYY - Year  (with 4 digits)'
	write(*,*) '    MM - Month (with 2 digits)'
	write(*,*) '    DD - Day   (with 2 digits)'
	write(*,*) '    HH - Hour  (with 2 digits - Use 00)'
	write(*,*) '     D - Y if You want to download before extract'
	write(*,*) ''
	write(*,*) 'Example:'
	write(*,*) ''
	write(*,*) 'geraSM_1.0 2020 09 14 00 Y'
	write(*,*) ''

end subroutine help
