PROGRAM goes
IMPLICIT NONE

! USAGE: pgf90 -tp amd64 -mcmodel=medium -Mlarge_arrays -o rgoes rgoes.f90
! ./goes.exe f20173241745.v65.g13.filt 20173241745go.out MCD12_2013_0.072F.bin 201711201745 2013_M.binv2.bin 2013_DP.binv2.bin
! GOES Variables declaration

	CHARACTER(len=250)  :: input_goes
	CHARACTER(len=250)  :: output
        CHARACTER(len=250)  :: lulcmap
        CHARACTER(len=80)   :: cab(46)
        CHARACTER(len=255)  :: line
        CHARACTER(len=13)    :: date

        INTEGER		    :: eco
        INTEGER		    :: fflag
        INTEGER		    :: conf
        INTEGER		    :: stat
        INTEGER             :: i
        INTEGER             :: nline
        INTEGER             :: nelement

        REAL 		    :: lat_g
        REAL 		    :: lon_g
        REAL 		    :: satzen
        REAL 		    :: solzen
        REAL 		    :: relazi
        REAL 		    :: psize
        REAL 		    :: t4
        REAL 		    :: t11
        REAL 		    :: nbt4
        REAL 		    :: nbt11
        REAL 		    :: fsize
        REAL 		    :: ftemp
        REAL 		    :: frp
        REAL 		    :: frpf
        REAL 		    :: flag
        REAL, PARAMETER     :: c1_planck = 3.74151e8
        REAL, PARAMETER     :: c2_planck = 1.43879e4
        REAL, PARAMETER     :: characteristics_goes = 3.08e-9
        REAL, PARAMETER     :: sigma = 5.67e-8
        REAL, PARAMETER     :: bgv = 185.57
        REAL                :: integralBegin
        REAL                :: integralEnd
        REAL                :: res
        REAL                :: integral
        REAL                :: bma
        REAL                :: bme
        REAL                :: deltat
        LOGICAL             :: existe
  
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
   REAL :: fsize_g

! Variables to group fire locations in a specific grid
! declaration values in lon and lat to globe and resolution of 20 km
! Sould be changed acccording to resolution
! 360 / resolution of data in degrees
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
   REAL    :: lon
   REAL    :: lat
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Call the file data

! VIIRS data
	CALL GETARG(1,input_goes)

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
   IF (input_goes .EQ. ' ' .OR. output .EQ. ' ' .OR. lulcmap .EQ. ' ' .OR. date .EQ. ' ' .OR. tave .EQ. ' ' .OR. sdave .EQ. ' ' ) THEN
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

              print *,juldaydec


! Open Files  --- start in 21

! GOES data
! OPEN GOES DATA
!LOOK in Header for the time
!Convert Character to integer
        INQUIRE(file=trim(input_goes), exist=existe)

        print*, trim(input_goes),existe
        if (.not.existe)then
           print*,'File not found ', trim(input_goes),existe
           stop
        end if
       
        OPEN(21,file=trim(input_goes),form="formatted", status='old', IOSTAT=stat, ERR=100, ACTION='read')

! Output file
        OPEN(22,file=trim(output), status='replace', IOSTAT=stat, ERR=120, ACTION='write')

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

  READ(21,'(a80)') (cab(i),i=1,46)

  WRITE(22,*) 'CABECALHO 1 - APENAS PARA ENTRADA NO PREP-CHEM-TRAD'
  WRITE(22,*) 'GOES-16 WF_ABBA (vs 6.5.008) Experimental Filtered Fire Product'
  WRITE(22,*) 'CABECALHO 3 - APENAS PARA ENTRADA NO PREP-CHEM-TRAD'
  WRITE(22,*) 'CABECALHO 4 - APENAS PARA ENTRADA NO PREP-CHEM-TRAD'
  WRITE(22,*) 'CABECALHO 5 - APENAS PARA ENTRADA NO PREP-CHEM-TRAD'
   
  DO

      READ(21, *, end=999) lat_g, lon_g, fflag, frp, fsize, ftemp, nline, nelement, psize, t4, t11, nbt4, nbt11, solzen, satzen, relazi, eco

   IF ((lon_g .GT. -90.0) .AND. (lon_g .LT. -30.0)) THEN
    IF ((lat_g .GT. -60.0) .AND. (lat_g .LT. 15.0)) THEN
      IF (fflag .EQ. 10) THEN
        IF (fsize .GT. 0) THEN
          psize = psize/1e6
          fsize = fsize/1e6
          fflag = 1
          WRITE(22,700) lon_g, lat_g, satzen, psize,t4,t11,fsize,ftemp,frp,eco,fflag
        ENDIF

700 FORMAT (9f20.4,2I7)
              

      ENDIF
     ENDIF
   ENDIF
  ENDDO



!				write(22,800) lon_g, lat_g, frpf, hour, psize

999 CONTINUE




 CLOSE(21)
 CLOSE(22)
 CLOSE(23)
 CLOSE(24)
 CLOSE(25)



        print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        print *, '!!!!!         FINISHED        !!!!!'
        print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'


  stop 0

100 WRITE (*, 500) 'GOES file end encontered ', trim(input_goes), &
       ' for reading'
  WRITE (*, 600) 'iostat = ', stat
  stop 2

110 WRITE (*, 400) 'Error: while reading from ', trim(input_goes)
  WRITE (*, 600) 'iostat = ', stat
  stop 3

120 WRITE (*, 500) 'Error: while opening ', trim(output), &
       ' for writing'
  WRITE (*, 600) 'iostat = ', stat
  stop 5

        print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        print *, '!!!!!         FINISHED        !!!!!'
        print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

400 FORMAT (1x, 2a)
500 FORMAT (1x, 3a)
600 FORMAT (8x, 1a, 1i7)


800 FORMAT (3f15.6,A10,1f15.6)
900 FORMAT ( 2f10.2, f12.2, f14.4, 2f9.1, f12.4, 2f12.1, i9, i6 )
DEALLOCATE(lulc)
DEALLOCATE(tave_diurnal)
DEALLOCATE(sdave_diurnal)
END PROGRAM goes
