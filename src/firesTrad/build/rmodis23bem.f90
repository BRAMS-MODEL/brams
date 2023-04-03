PROGRAM modis
        USE HDF4UTILS
	IMPLICIT NONE

! USAGE: pgf90 -tp amd64 -mcmodel=medium -Mlarge_arrays -o pmodis hdf4utils.f90 rmodis.f90
!-L/HDF-4.2.5/lib -lmfhdf -ldf -L/jpeg-6b/lib -ljpeg -L/zlib-1.2.4/lib -lz
! Variables declaration
! ./modis.exe MOD14.A2015239.0155.006.2015305094007.hdf  20173241745mo.out MCD12_2013_0.072F.bin 201711201745 2013_M.binv2.bin 2013_DP.binv2.bin
	INTEGER	:: sfstart, sfend

! INPUT and OUTPUT variables declaration
	CHARACTER(len=250)  :: input_modis
	CHARACTER(len=250)  :: output
        CHARACTER(len=250)  :: lulcmap
        CHARACTER(len=13)   :: date
        INTEGER :: stat
! Variables to extract MODIS HDF from HDF4 utilities

	REAL, DIMENSION(:), POINTER :: lat
	REAL, DIMENSION(:), POINTER :: lon
	REAL, DIMENSION(:), POINTER :: frp_m
	INTEGER*2, DIMENSION(:), POINTER :: line
	INTEGER*2, DIMENSION(:), POINTER :: sample
	REAL, DIMENSION(:), POINTER :: r2
	REAL, DIMENSION(:), POINTER :: t21
	REAL, DIMENSION(:), POINTER :: t31
	REAL, DIMENSION(:), POINTER :: meant21
	REAL, DIMENSION(:), POINTER :: meant31
	INTEGER*1, DIMENSION(:), POINTER :: cloud
	INTEGER*1, DIMENSION(:), POINTER :: water
	INTEGER*1, DIMENSION(:), POINTER :: window
	INTEGER*2, DIMENSION(:), POINTER :: numvalid
	INTEGER*1, DIMENSION(:), POINTER :: confidence

! Parameter configuration

	INTEGER, PARAMETER :: DFACC_READ = 1

! Auxiliary variables

	INTEGER		:: sd_id, i, ierr, lim
	INTEGER*1	:: uimax
	INTEGER*2	:: imax
	REAL		:: maximum, Re, r, s, N, TETA, Z, YY, XX, W, DeltaS, DeltaA, AREA, frpc
        REAL		:: lona,lata,ya,xa

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
        REAL :: fsize_m


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
  CHARACTER(len=2)    :: chour
  CHARACTER(len=2)    :: cminute
  CHARACTER(len=4)    :: yearc
  CHARACTER(len=2)    :: monc
  CHARACTER(len=2)    :: dayc
  REAL(4)             :: dhour,dhh,dyy,dmm,dday
  REAL(4)             :: juldaydec
  INTEGER             :: dmon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Variables calling by command line in download script
! Should be <modis input file > <output file - in YYYYDDD> <LULC map>

  CALL GETARG(1,input_modis)

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
   IF (input_modis .EQ. ' ' .OR. output .EQ. ' ' .OR. lulcmap .EQ. ' ' .OR. date .EQ. ' ' .OR. tave .EQ. ' ' .OR. sdave .EQ. ' ' ) THEN
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

! Open MODIS14 FILE

        WRITE(*,'("Opening file ",A)') trim(input_modis)
	sd_id = sfstart(input_modis, DFACC_READ)
	CALL checkError(sd_id, "Error in opening input MODIS file")

! Opening the MODIS HDF4 layers
  CALL hdf4_read(sd_id, "FP_latitude", lat)
  CALL hdf4_read(sd_id, "FP_longitude", lon)
  CALL hdf4_read(sd_id, "FP_power", frp_m)
  CALL hdf4_read(sd_id, "FP_line", line)
  CALL hdf4_read(sd_id, "FP_sample", sample)
  CALL hdf4_read(sd_id, "FP_R2", r2)
  CALL hdf4_read(sd_id, "FP_T21", t21)
  CALL hdf4_read(sd_id, "FP_T31", t31)
  CALL hdf4_read(sd_id, "FP_MeanT21", meant21)
  CALL hdf4_read(sd_id, "FP_MeanT31", meant31)
  CALL hdf4_read(sd_id, "FP_AdjCloud", cloud)
  CALL hdf4_read(sd_id, "FP_AdjWater", water)
  CALL hdf4_read(sd_id, "FP_WinSize", window)
  CALL hdf4_read(sd_id, "FP_NumValid", numvalid)
  CALL hdf4_read(sd_id, "FP_confidence", confidence)
	
! Final File acess
	ierr = sfend(sd_id)
	CALL checkError(ierr, "Error in close MODIS file")
	
! Data analysis

  maximum = MAXVAL(lat)
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "FP_latitude", maximum
  maximum = MAXVAL(lon)  
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "FP_longitude", maximum
  maximum = MAXVAL(frp_m)
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "FP_power", maximum
  imax = MAXVAL(line)
!  WRITE(*,'("Maximum value in ",A," is: ",I6)') "FP_line", imax
  imax = MAXVAL(sample)
!  WRITE(*,'("Maximum value in ",A," is: ",I6)') "FP_sample", imax
  maximum = MAXVAL(r2)  
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "FP_R2", maximum
  maximum = MAXVAL(t21)  
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "FP_T21", maximum
  maximum = MAXVAL(t31)  
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "FP_T31", maximum
  maximum = MAXVAL(meant21)  
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "FP_MeanT21", maximum
  maximum = MAXVAL(meant31)  
!  WRITE(*,'("Maximum value in ",A," is: ",f10.5)') "FP_MeanT31", maximum
  uimax = MAXVAL(cloud)
!  WRITE(*,'("Maximum value in ",A," is: ",I6)') "FP_AdjCloud", uimax
  uimax = MAXVAL(water)
!  WRITE(*,'("Maximum value in ",A," is: ",I6)') "FP_AdjWater", uimax
  uimax = MAXVAL(window)
!  WRITE(*,'("Maximum value in ",A," is: ",I6)') "FP_WinSize", uimax
	imax = MAXVAL(numvalid)
!  WRITE(*,'("Maximum value in ",A," is: ",I6)') "FP_NumValid", imax
	uimax = MAXVAL(confidence)
!  WRITE(*,'("Maximum value in ",A," is: ",I6)') "FP_confidence", uimax

! Output file Sequence order 2  --

  OPEN(22,file=trim(output), status='replace', IOSTAT=ierr, ACTION='write')

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

! MODIS hdfs have the same dimensions to all variables so just one need.
  lim = SIZE(lat)
  IF(lim .GT. 0) THEN
    DO i=1, lim      
!     	IF (confidence(i) .GE. 40) THEN
!            CONTINUE
!        ELSE
	    	! Calculo da Â rea do pixel para o MODIS

	      Re = 6378.137
	      r = 7083.137
	      s = 0.0014184
	      N = 1354.0
	      TETA = ((-0.5 * N * s) + (0.5 * s) + ((sample(i) - 1) * s))
	       
	      Z = 0.810842
	      YY = (SIN(TETA) * SIN(TETA))
				XX = Z - YY
	      W = (SQRT (XX))

	      DeltaS = ((Re * s) * ((COS(TETA)/W)-1))
	      DeltaA = ((r * s) * (COS(TETA) - W))

	      AREA = (DeltaS * DeltaA)

	    IF (frp_m(i) .GT. 1.) THEN
               IF ((lon(i) .GT. -90.0) .AND. (lon(i) .LT. -30.0)) THEN
                  IF ((lat(i) .GT. -60.0) .AND. (lat(i) .LT. 15.0)) THEN

                     WRITE(22,700) lon(i), lat(i)

             ENDIF
          ENDIF
        ENDIF
      ENDDO
  ENDIF






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

800 FORMAT (3f15.6,A10,1f15.6)
!999 CONTINUE

DEALLOCATE(lulc)
DEALLOCATE(tave_diurnal)
DEALLOCATE(sdave_diurnal)

DEALLOCATE(lat, lon, frp_m, line, sample, r2, t21, t31, meant21, meant31 &
           ,cloud, water, window, numvalid, confidence)

  WRITE(*,'("MODIS data processed.")')

END PROGRAM modis

