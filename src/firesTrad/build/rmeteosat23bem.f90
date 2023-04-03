PROGRAM meteosat
	USE HDF5
        USE HDF5UTILS

	IMPLICIT NONE

! USAGE:pgf90 -tp amd64 -mcmodel=medium -Mlarge_arrays -o rmeteosat hdf5utils.f90
! rmeteosat.f90 -L/scratchin/grupos/catt-brams/shared/libs/pgi/hdf5-1.8.8/lib
!-lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -I/scratchin/grupos/catt-brams/shared/libs/pgi/hdf5-1.8.8/include
!-L/scratchin/grupos/catt-brams/shared/libs/pgi/jpeg-6b/lib -ljpeg
!-L/scratchin/grupos/catt-brams/shared/libs/pgi/zlib-1.2.4/lib -lz


! Variables declaration

        CHARACTER(len=250)  :: input_meteosat
	CHARACTER(len=250)  :: output
        CHARACTER(len=250)  :: lulcmap

! Input for HDF5 read


  INTEGER :: hdferr, stat
  LOGICAL         :: avail
  INTEGER(HID_T)  :: file, space, dset, attr ! Handles
  INTEGER, DIMENSION(1:1) :: cd_values


  INTEGER(SIZE_T) :: nelmts
  INTEGER :: flags, filter_info
  REAL, DIMENSION(:), POINTER :: frp
  REAL, DIMENSION(:), POINTER :: lat
  REAL, DIMENSION(:), POINTER :: lon
  REAL, DIMENSION(:), POINTER :: hora
  REAL, DIMENSION(:), POINTER :: btmir
  REAL, DIMENSION(:), POINTER :: bttir
  REAL, DIMENSION(:), POINTER :: bwbtd
  REAL, DIMENSION(:), POINTER :: bwbtmir
  REAL, DIMENSION(:), POINTER :: bwnpix
  REAL, DIMENSION(:), POINTER :: bwsz
  REAL, DIMENSION(:), POINTER :: conf
  REAL, DIMENSION(:), POINTER :: patm
  REAL, DIMENSION(:), POINTER :: siz
  REAL, DIMENSION(:), POINTER :: vza
  REAL :: max
  INTEGER :: i, ny, lim
  INTEGER(SIZE_T), PARAMETER :: MaxChrLen = 80
  CHARACTER(LEN=MaxChrLen) :: name
  CHARACTER(len=4)  :: chora

! Auxiliary variables

        REAL :: lona,lata,ya,xa


! Variables to cross fire location with land use and land cover
! Uso do solo MCD12Q1 com resulu‡Æo variada para cada dado de entrada
! MODIS E VIIRS 2 km e GOES 8 km. 0,072

	INTEGER, PARAMETER ::  n_cols= 867
	INTEGER, PARAMETER ::  n_rows= 1008


! Spatioal Resolution of LULC map
!! DO NOT FORGET TO CHANGE IT IF MODIFICATIONS IN LULC PRODUCT WAS MADE
   REAL  :: resx
   REAL  :: resy
! Position in cartesian matrix
   INTEGER  :: posix
   INTEGER  :: posiy

! Test to lower data
!	INTEGER, PARAMETER ::  n_cols= 7200
!	INTEGER, PARAMETER ::  n_rows= 3600

! Defined lulc variable
        REAL, ALLOCATABLE, DIMENSION(:,:) :: lulc
        REAL :: igbp
        REAL :: fsize_m

! Variables to group fire locations in a specific grid
! declaration values in lon and lat to globe and resolution of 3 km
! Sould be changed acccording to resolution
! 360 / resolution of data in degrees  - EX: 3 / 111.12
! 180 / resolution of data in degrees

	INTEGER, PARAMETER ::  a_cols= 2010
	INTEGER, PARAMETER ::  a_rows= 1005

        INTEGER :: cont(a_cols,a_rows)
	REAL    :: ac_frp(a_cols,a_rows)
        REAL    :: ac_area(a_cols,a_rows)
        REAL    :: ac_fsize(a_cols,a_rows)
        REAL    :: ac_hour(a_cols,a_rows)

        REAL    :: ac_diurnal(a_cols,a_rows)
        REAL    :: idiurnal



        REAL    :: delta


        INTEGER :: x
        INTEGER :: y

        INTEGER :: ihour
        CHARACTER(len=4)    :: spc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dirunal information assimilation
! maps for South America with resolutions according to Input data
! Media e desvio padrao
! Proximo passo usar Skewness e Kurtosis

   CHARACTER(len=250)  :: tave
   CHARACTER(len=250)  :: sdave
   INTEGER, PARAMETER ::  dm_cols= 866
   INTEGER, PARAMETER ::  dm_rows= 1007

   REAL  :: resxmd
   REAL  :: resymd
! Position in cartesian matrix
   INTEGER  :: posixmd
   INTEGER  :: posiymd


! Defined dirunal variable

   REAL, ALLOCATABLE, DIMENSION(:,:) :: tave_diurnal
   REAL, ALLOCATABLE, DIMENSION(:,:) :: sdave_diurnal
   REAL :: tavep
   REAL :: sdavep

! TIME VAR
  CHARACTER(len=2)    :: chour
  CHARACTER(len=2)    :: cminute
  CHARACTER(len=4)    :: yearc
  CHARACTER(len=2)    :: monc
  CHARACTER(len=2)    :: dayc
  REAL(4)             :: dhour,dhh,dyy,dmm,dday
  REAL(4)             :: juldaydec
  INTEGER             :: dmon
  CHARACTER(len=13)   :: date
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Variables calling by command line in download script
! Should be <modis input file > <output file - in YYYYDDD> <LULC map>

  CALL GETARG(1,input_meteosat)
  CALL GETARG(2,output)
! MAP file in binary and float point
  CALL GETARG(3,lulcmap)

  CALL GETARG(4,date)

! MAP file in binary and float point
  CALL GETARG(5,tave)

! MAP file in binary and float point
  CALL GETARG(6,sdave)

!Error if number of described files are lower then need to process
   IF (input_meteosat .EQ. ' ' .OR. output .EQ. ' ' .OR. lulcmap .EQ. ' ' .OR. date .EQ. ' ' .OR. tave .EQ. ' ' .OR. sdave .EQ. ' ' ) THEN
     WRITE(*,*) 'use: program <input_goes> <output> <lulc map> <hour> <diurnal M> <diurnal DP>'
     STOP
   ENDIF
        
 ! EXTRACT HOUR

 ! ARRUMAR A PARTIR DAQUI
 ! FORMAT YYYYMMDDHHmm


   yearc   = date(1:4)
   monc    = date(5:6)
   dayc    = date(7:8)
   chour   = date(9:10)
   cminute = date(11:12)

   read(yearc(1:4),*)   dyy
   read(monc(1:2),*)    dmon
   read(dayc(1:2),*)    dday
   read(chour(1:2),*)   dhh
   read(cminute(1:2),*) dmm
   
! DECIMAL HOUR HH.MMM
   dhour = dhh+(dmm/60.)

              IF ((dyy .EQ. 1996) .OR. (dyy .EQ. 2000) .OR. (dyy .EQ. 2004) .OR. &
                 (dyy .EQ. 2008) .OR. (dyy .EQ. 2012) .OR. (dyy .EQ. 2016) .OR. &
                 (dyy .EQ. 2020) .OR. (dyy .EQ. 2024)) THEN

                  SELECT CASE (dmon)
                  CASE ( 1 )
                    juldaydec = dday
                  CASE ( 2 )
                    juldaydec = 31.+dday
                  CASE ( 3 )
                    juldaydec = 60.+dday
                  CASE ( 4 )
                    juldaydec = 91.+dday
                  CASE ( 5 )
                    juldaydec = 121.+dday
                  CASE ( 6 )
                    juldaydec = 152.+dday
                  CASE ( 7 )
                    juldaydec = 182.+dday
                  CASE ( 8 )
                    juldaydec = 213.+dday
                  CASE ( 9 )
                    juldaydec = 244.+dday
                  CASE ( 10 )
                    juldaydec = 274.+dday
                  CASE ( 11 )
                    juldaydec = 305.+dday
                  CASE ( 12 )
                    juldaydec = 335.+dday
                  END SELECT

              ELSE

                  SELECT CASE (dmon)
                  CASE ( 1 )
                    juldaydec = dday
                  CASE ( 2 )
                    juldaydec = 31.+dday
                  CASE ( 3 )
                    juldaydec = 59.+dday
                  CASE ( 4 )
                    juldaydec = 90.+dday
                  CASE ( 5 )
                    juldaydec = 120.+dday
                  CASE ( 6 )
                    juldaydec = 151.+dday
                  CASE ( 7 )
                    juldaydec = 181.+dday
                  CASE ( 8 )
                    juldaydec = 212.+dday
                  CASE ( 9 )
                    juldaydec = 243.+dday
                  CASE ( 10 )
                    juldaydec = 273.+dday
                  CASE ( 11 )
                    juldaydec = 304.+dday
                  CASE ( 12 )
                    juldaydec = 334.+dday
                  END SELECT

              ENDIF

              juldaydec = juldaydec*100+dhour

! Open METEOSAT FILE

 CALL h5open_f(hdferr)
  WRITE(*,'("Abrindo arquivo ",A)') trim(input_meteosat)
	CALL h5fopen_f(trim(input_meteosat), H5F_ACC_RDWR_F, file, hdferr)
  CALL checkError(hdferr, 'Erro ao abrir arquivo')

  CALL hdf5_read(file, "FRP", frp)
  CALL hdf5_read(file, "LATITUDE", lat)
  CALL hdf5_read(file, "LONGITUDE", lon)
  CALL hdf5_read(file, "ACQTIME", hora)
  CALL hdf5_read(file, "BT_MIR", btmir)
  CALL hdf5_read(file, "BT_TIR", bttir)
  CALL hdf5_read(file, "BW_BTD", bwbtd)
  CALL hdf5_read(file, "BW_BT_MIR", bwbtmir)
  CALL hdf5_read(file, "BW_NUMPIX", bwnpix)
  CALL hdf5_read(file, "BW_SIZE", bwsz)
  CALL hdf5_read(file, "FIRE_CONFIDENCE", conf)
  CALL hdf5_read(file, "PIXEL_ATM_TRANS", patm)
  CALL hdf5_read(file, "PIXEL_SIZE", siz)
  CALL hdf5_read(file, "PIXEL_VZA", vza)

  CALL h5fclose_f(file , hdferr)
  CALL h5close_f(hdferr)

  max = MAXVAL(frp)
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "FRP", max
  max = MAXVAL(lat)
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "lat", max
  max = MAXVAL(lon)
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "lon", max
  max = MAXVAL(hora)
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "hora", max
  max = MAXVAL(btmir)
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "BT_MIR", max
  max = MAXVAL(bttir)
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "BT_TIR", max
  max = MAXVAL(bwbtd)
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "BW_BTD", max
  max = MAXVAL(bwbtmir)
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "BW_BT_MIR", max
  max = MAXVAL(bwnpix)
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "BW_NUMPIX", max
  max = MAXVAL(bwsz)
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "BW_SIZE", max
  max = MAXVAL(conf)
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "FIRE_CONFIDENCE", max
  max = MAXVAL(patm)
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "PIXEL_ATM_TRANS", max
  max = MAXVAL(siz)
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "PIXEL_SIZE", max
  max = MAXVAL(vza)
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "PIXEL_VZA", max

! Output file Sequence order 2  --

	OPEN(22,file=trim(output), status='replace', IOSTAT=stat, ACTION='write')

! Input MAP of Land Use and Land cover
        OPEN(23, file=lulcmap, FORM='UNFORMATTED', ACCESS="STREAM", IOSTAT=stat)

! Allocate MAP file
	ALLOCATE(lulc(n_cols,n_rows))
!        print *,input_viirs,output,lulcmap

! Read the MAP with land use and land cover IGBP

        READ (23) lulc

        print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        print *, '!!!!!  Read of LULC MAP DONE  !!!!!'
        print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

!        print *, ' Processing GOES data', input_goes

! LULC MAP Resolution
!        resx = 0.003871
!        resy = 0.003871

        resx = 0.072
        resy = 0.072

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Input MAP of DIURNAL TIME
   OPEN(24, file=tave, FORM='UNFORMATTED', ACCESS="STREAM", IOSTAT=stat)

! Allocate AVERAGE TIME OF BURNING  file
   ALLOCATE(tave_diurnal(dm_cols,dm_rows))

! Read the MAP with average time according to Landuse
! SEE xxxx for more informations

   READ (24) tave_diurnal

   print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   print *, '!!!!! Read of AVERAGE TIME DONE !!!'
   print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

! LULC MAP Resolution

   resxmd = 0.072
   resymd = 0.072
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Input MAP of DIURNAL TIME  STANDARD DEVIATION
   OPEN(25, file=sdave, FORM='UNFORMATTED', ACCESS="STREAM", IOSTAT=stat)

! Allocate AVERAGE TIME OF BURNING  file
   ALLOCATE(sdave_diurnal(dm_cols,dm_rows))

! Read the MAP with average time according to Landuse
! SEE xxxx for more informations

   READ (25) sdave_diurnal

   print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   print *, '!!!!!  Read of SDEV TIME DONE  !!!!'
   print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'


! DELTA resolution
! COnvert km to degrees
! In near future this variable will be direct assessed in prep_chem.inp

  delta = 20./111.12
  ac_frp(:,:) = 0.
  ac_area(:,:) = 0.
  cont(:,:) = 0
  ac_fsize(:,:) = 0.
  ac_hour(:,:) = 0
  ac_diurnal(:,:) = 0.

! MODIS hdfs have the same dimensions to all variables so just one need.

  lim = SIZE(lat)
  IF(lim .GT. 0) THEN
    DO i=1, lim      
!     	IF (confidence(i) .GE. 40) THEN
!            CONTINUE
!        ELSE

	    IF (frp(i) .GT. 1.) THEN
	       IF (conf(i) .GE. 0.85) THEN
                   IF ((lon(i) .GT. -90.0) .AND. (lon(i) .LT. -30.0)) THEN
                     IF ((lat(i) .GT. -60.0) .AND. (lat(i) .LT. 15.0)) THEN

                        ! GET the position in matrix to determine the LULC

                         WRITE(22,700) lon(i), lat(i),siz(i)


                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
     ENDIF


 CLOSE(22)
 CLOSE(23)
 CLOSE(24)
 CLOSE(25)

 ac_frp(:,:) = 0.
 ac_area(:,:) = 0.
 cont(:,:) = 0
 ac_fsize(:,:) = 0.
 ac_hour(:,:) = 0
 ac_diurnal(:,:) = 0.



        print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        print *, '!!!!!         FINISHED        !!!!!'
        print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

400 FORMAT (1x, 2a)
500 FORMAT (1x, 3a)
600 FORMAT (8x, 1a, 1i7)

700 FORMAT (3f20.8)


800 FORMAT (3f15.6,A10,1f15.6)
!999 CONTINUE

DEALLOCATE(lulc)
DEALLOCATE(tave_diurnal)
DEALLOCATE(sdave_diurnal)

  DEALLOCATE(lat, lon, frp, hora, conf          &
        , siz, btmir, bttir, bwbtd, bwbtmir     &
        , bwnpix, bwsz, patm, vza)

  WRITE(*,'("METEOSAT data processed.")')

END PROGRAM meteosat

SUBROUTINE checkError(hdferr, msg)

  IMPLICIT NONE

  INTEGER :: hdferr
  CHARACTER(*) :: msg

  IF(hdferr .LT. 0) THEN
    WRITE(*,'("ERROR: ",A)') trim(msg)
  ENDIF
END SUBROUTINE checkError


  