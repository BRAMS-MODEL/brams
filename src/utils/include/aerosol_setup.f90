INTEGER,PARAMETER :: local_nspecies=6
    INTEGER,PARAMETER :: local_nmodes=10

    ! spaction(specie,[1=source,2=drydep,3=wetdep,4=fdda, 5=offline emission, 6=transport])
    ! attention : for aerosols,  mode_alloc(ispc) = spc_alloc(transport,imode,ispc)
    INTEGER,PARAMETER,DIMENSION(6,local_nmodes,local_nspecies) :: local_spc_alloc=RESHAPE((/ &
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer sdust bin 1
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer sdust bin 2
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer sdust bin 3
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer sdust bin 4
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer sdust bin 5
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer sdust bin 6
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer sdust bin 7
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer sdust bin 8
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer sdust bin 9
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer sdust bin 10

  !
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer bburn bin 1  - PM0
      1 , 1 , 1 , 0 , 0 , 1 ,  & ! aer bburn bin 2  - PM2.5
      1 , 1 , 1 , 0 , 0 , 1 ,  & ! aer bburn bin 3  - PM10
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer bburn bin 4
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer bburn bin 5
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer bburn bin 6
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer bburn bin 7
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer bburn bin 8
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer bburn bin 9
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer bburn bin 10
  !
  !
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer urban bin 1
      1 , 1 , 1 , 0 , 0 , 1 ,  & ! aer urban bin 2 ! urban PM2.5
      1 , 1 , 1 , 0 , 0 , 1 ,  & ! aer urban bin 3 ! continental
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer urban bin 4
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer urban bin 5
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer urban bin 6
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer urban bin 7
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer urban bin 8
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer urban bin 9
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer urban bin 10
  !
  !
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer bioge bin 1
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer bioge bin 2
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer bioge bin 3
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer bioge bin 4
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer bioge bin 5
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer bioge bin 6
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer bioge bin 7
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer bioge bin 8
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer bioge bin 9
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer bioge bin 10
  !
      1 , 1 , 1 , 0 , 1 , 1 ,  & ! aer marin bin 1
      1 , 1 , 1 , 0 , 1 , 1 ,  & ! aer marin bin 2
      1 , 1 , 1 , 0 , 1 , 1 ,  & ! aer marin bin 3
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer marin bin 4
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer marin bin 5
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer marin bin 6
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer marin bin 7
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer marin bin 8
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer marin bin 9
      0 , 0 , 0 , 0 , 0 , 0 ,  & ! aer marin bin 10
  !
      1 , 1 , 0 , 0 , 0 , 1 ,  & ! aer v_ash bin 1
      1 , 1 , 0 , 0 , 0 , 1 ,  & ! aer v_ash bin 2
      1 , 1 , 0 , 0 , 0 , 1 ,  & ! aer v_ash bin 3
      1 , 1 , 0 , 0 , 0 , 1 ,  & ! aer v_ash bin 4
      1 , 1 , 0 , 0 , 0 , 1 ,  & ! aer v_ash bin 5
      1 , 1 , 0 , 0 , 0 , 1 ,  & ! aer v_ash bin 6
      1 , 1 , 0 , 0 , 0 , 1 ,  & ! aer v_ash bin 7
      1 , 1 , 0 , 0 , 0 , 1 ,  & ! aer v_ash bin 8
      1 , 1 , 0 , 0 , 0 , 1 ,  & ! aer v_ash bin 9
      1 , 1 , 0 , 0 , 0 , 1    & ! aer v_ash bin 10
  !
  !
      /),(/6,local_nmodes,local_nspecies/))

    ! effective particle radius (meter)
    REAL,PARAMETER,DIMENSION(local_nmodes,local_nspecies) :: local_part_radius=RESHAPE((/ &
  !------------------------------------------------------------------------------------------------------------------
  !-bin size  1          2        3        4          5       6       7       8     9    10
  !------------------------------------------------------------------------------------------------------------------
       	  1.95e-7 , 1.95e-7 , 1.95e-7 , 999., 999., 999., 999., 999., 999., 999.,   & ! sdust
       	  1.95e-7 , 1.95e-7 , 1.00e-5 , 999., 999., 999., 999., 999., 999., 999.,   & ! bburn (0, pm25, pm10) meters
       	  1.95e-7 , 1.95e-7 , 1.95e-7 , 999., 999., 999., 999., 999., 999., 999.,   & ! urban
       	  1.95e-7 , 1.95e-7 , 1.95e-7 , 999., 999., 999., 999., 999., 999., 999.,   & ! bioge
       	  8.25e-8 , 2.82e-7 , 1.61e-6 , 999., 999., 999., 999., 999., 999., 999.,   & ! marin
       	  0.98e-6,  2.93e-6,  5.89e-6, 11.72e-6, 23.44e-6, 46.88e-6, 93.75e-6, 0.1875e-3, 0.375e-3, 0.750e-3   & ! v_ash
  !---------------------------------------------------------- ---------------------------- ----------------------------
      /),(/local_nmodes,local_nspecies/))

    ! particle density kg/m^3
    !srf-dez2013: changing density for ash from 2500 to 900 kg/m3
    REAL,PARAMETER,DIMENSION(local_nmodes,local_nspecies) :: local_part_dens=RESHAPE((/ &
  !------------------------------------------------------------------------------------------------------------------
  !-bin size 1   2   3  4   5   6   7   8   9  10
  !------------------------------------------------------------------------------------------------------------------
      2.65e+3 , 2.65e+3 , 2.65e+3 , 999., 999., 999., 999., 999., 999., 999.,   & ! sdust
      1.35e+3 , 1.35e+3 , 1.35e+3 , 999., 999., 999., 999., 999., 999., 999.,   & ! bburn (0, pm25, pm10) kg/m^3
      1.35e+3 , 1.35e+3 , 1.35e+3 , 999., 999., 999., 999., 999., 999., 999.,   & ! urban
      1.35e+3 , 1.35e+3 , 1.35e+3 , 999., 999., 999., 999., 999., 999., 999.,   & ! bioge
      2.17e+3 , 2.17e+3 , 2.17e+3 , 999., 999., 999., 999., 999., 999., 999.,   & ! marin
      0.90e+3 , 0.90e+3 , 0.90e+3 , 0.90e+3, 0.90e+3, 0.90e+3, 0.90e+3, 0.90e+3, 0.90e+3, 0.90e+3    & ! v_ash
  !    2.50e+3 , 2.50e+3 , 2.50e+3 , 2.50e+3, 2.50e+3, 2.50e+3, 2.50e+3, 2.50e+3, 2.50e+3, 2.50e+3    & ! v_ash
  !------------------------------------------------------------------------------------------------------------------
      /),(/local_nmodes,local_nspecies/))
