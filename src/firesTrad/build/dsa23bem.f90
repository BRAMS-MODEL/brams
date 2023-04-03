PROGRAM rdsa
  IMPLICIT NONE

  CHARACTER(len=250)  :: input_dsa
  CHARACTER(len=550)  :: output
  CHARACTER(len=15)   :: dsa_date
  REAL                :: latitude,longitude,dsa_frp
  CHARACTER(len=11)    :: dsa_sat
  INTEGER :: stat



  CALL GETARG(1,input_dsa)

  CALL GETARG(2,output)

   IF (input_dsa .EQ. ' ' .OR. output .EQ. ' ' ) THEN
     WRITE(*,*) 'MISSING ARGUMENTS TO PROCESS DSA FILES ::: STOPPED'
     WRITE(*,*) 'NEED: INPUT / OUTPUT'
   STOP
   ENDIF



            OPEN(21,file=trim(input_dsa),status='old', IOSTAT=stat, ACTION='read')
            ! Output file
            OPEN(22,file=trim(output), status='replace', ACTION='write')

            DO
               READ(21, *, iostat=stat, end=999)latitude,longitude,dsa_date,dsa_sat,dsa_frp

                          IF ( (dsa_sat .EQ. "NOAA-18D") .OR. &
                               (dsa_sat .EQ. "NOAA-19D") .OR. &
                               (dsa_sat .EQ. "NOAA-18") .OR. &
                               (dsa_sat .EQ. "NOAA-19") .OR. &
                               (dsa_sat .EQ. "NOAA-15D") .OR. &
                               (dsa_sat .EQ. "NOAA-15") .OR. &
                               (dsa_sat .EQ. "NPP_375") .OR. &
                               (dsa_sat .EQ. "METOP-B") .or. &
                               (dsa_sat .EQ. "TERRA_M-M") .or. &
                               (dsa_sat .EQ. "TERRA_M-T") .or. &
                               (dsa_sat .EQ. "AQUA_M-M") .or. &
                               (dsa_sat .EQ. "AQUA_M-T") &
                               ) THEN

                                WRITE(22,700) latitude, longitude

                          ENDIF

           ENDDO
999 CONTINUE


          CLOSE(21)
          CLOSE(22)

    print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print *, '!!!!!         FINISHED        !!!!!'
    print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

700 FORMAT (2f20.8)

END PROGRAM rdsa
