PROGRAM rviirs
IMPLICIT NONE

! The preprocessor for the VIIRS FRP data, originally from G.Pereira, 2017
! USAGE: pgf90 -tp amd64 -mcmodel=medium -Mlarge_arrays -o procviirs rviirs.f90

! Variables declaration

  CHARACTER(len=250)  :: input_viirs
  CHARACTER(len=250)  :: output
  CHARACTER(len=250)  :: lulcmap
  CHARACTER(len=80)   :: cab(1)
  CHARACTER(len=255)  :: line
  INTEGER             :: hour
  INTEGER             :: minute

  INTEGER :: stat
  INTEGER :: yy
  INTEGER :: mm
  INTEGER :: dd
  INTEGER :: mask
  INTEGER :: confi
  INTEGER :: posy
  INTEGER :: posx
  INTEGER :: posxx
  INTEGER :: i

  REAL :: lon_vi
  REAL :: lat_vi
  REAL :: t13
  REAL :: frp_v
  REAL :: parea

 ! Variables to cross fire location with land use and land cover
! Uso do solo MCD12Q1 com resulu‡Æo variada para cada dado de entrada
! MODIS E VIIRS 2 km e GOES 8 km. 0,072

  INTEGER, PARAMETER ::  n_cols= 3120
  INTEGER, PARAMETER ::  n_rows= 3626


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
        REAL :: fsize_v

! Variables to group fire locations in a specific grid
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

   REAL :: delta
   INTEGER :: ihour
   REAL :: lon
   REAL :: lat
   INTEGER :: x
   INTEGER :: y
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dirunal information assimilation
! maps for South America with resolutions according to Input data
! Media e desvio padrao
! Proximo passo usar Skewness e Kurtosis

   CHARACTER(len=250)  :: tave
   CHARACTER(len=250)  :: sdave
   INTEGER, PARAMETER ::  dm_cols= 3119
   INTEGER, PARAMETER ::  dm_rows= 3625

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
  CHARACTER(len=13)   :: date
  CHARACTER(len=2)    :: chour
  CHARACTER(len=2)    :: cminute
  CHARACTER(len=4)    :: yearc
  CHARACTER(len=2)    :: monc
  CHARACTER(len=2)    :: dayc
  REAL(4)             :: dhour,dhh,dyy,dmm,dday
  REAL(4)             :: juldaydec
  INTEGER             :: dmon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Call the file data

! VIIRS data
	CALL GETARG(1,input_viirs)

!Output ascii file -- after process
        CALL GETARG(2,output)

! MAP file in binary and float point
        CALL GETARG(3,lulcmap)
        
! GET HOUR AND AVOID TO EXTRACT FROM HEADER
  CALL GETARG(4,date)

! MAP file in binary and float point
   CALL GETARG(5,tave)

! MAP file in binary and float point
   CALL GETARG(6,sdave)

!Error if number of described files are lower then need to process
   IF (input_viirs .EQ. ' ' .OR. output .EQ. ' ' .OR. lulcmap .EQ. ' ' .OR. date .EQ. ' ' .OR. tave .EQ. ' ' .OR. sdave .EQ. ' ' ) THEN
     WRITE(*,*) 'use: program <input_goes> <output> <lulc map> <hour> <diurnal M> <diurnal DP>'
     STOP
   ENDIF

! EXTRACT HOUR
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

! Open Files  --- start in 21

! VIIRS data
	OPEN(21,file=trim(input_viirs),status='old', IOSTAT=stat, ACTION='read')
        	READ(21,'(a80)') (cab(i),i=1,1)
! Output file
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

  resx = 0.02
  resy = 0.02

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

   resxmd = 0.02
   resymd = 0.02
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

! Start to read VIIRS file
! In format
! year,month,day,hh,mm,lon,lat,mask,confidence,bright_t13,frp,line,sample
! Data provided by NESDIS via: ftp://ftp.star.nesdis.noaa.gov/pub/smcd/tsidulko/HRRR/
! A sample VIIRS file: 2017, 10, 09, 18, 01,  -76.404167,   -5.720620,   8,  56,  314.097443,    8.996529,     2,  2563

      DO
  	READ(21, *, iostat=stat, end=999) yy,mm,dd,hour,minute,lon_vi,lat_vi,mask,confi,t13,frp_v,posy,posx
        IF (stat /= 0) EXIT

        IF ((lon_vi .GT. -100.0) .AND. (lon_vi .LT. -30.0)) THEN
          IF ((lat_vi .GT. -60.0) .AND. (lat_vi .LT. 15.0)) THEN
             IF (frp_v .GT. 0.0 ) THEN

                WRITE(22,700) lon_vi, lat_vi
         ENDIF
       ENDIF
      ENDIF

     ENDDO

999 CONTINUE

 CLOSE(21)
 CLOSE(22)
 CLOSE(23)
 CLOSE(24)
 CLOSE(25)


        print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        print *, '!!!!!         FINISHED        !!!!!'
        print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

400 FORMAT (1x, 2a)
500 FORMAT (1x, 3a)
600 FORMAT (8x, 1a, 1i7)

700 FORMAT (2f20.8)
!700	FORMAT (5i6,2f15.6,2i6,2f15.6,2i4)

800 FORMAT (3f15.6,A10,1f15.6)
!999 CONTINUE

DEALLOCATE(lulc)
END PROGRAM rviirs
