!-----------------------------------------------------------------------------*
!= This is Fast-TUV FTUV4.2
!=  with 11 bins starting from 205nm
!=
!=  In general the error between TUV and FTUV is within 5%, and in
!=     some cases, it is above 5%
!=  There are some limitations to use this FTUV
!=   (1) J for stratospheric species (28-57) is not applied
!=   (2) CH3CHO have only one channel CH3CHO + hv --> CH3 + HCO
!=   (3) O2, N2O, and HO2 have large erros (due to absorption in short W)
!=
!=    Tropospheric Ultraviolet-Visible (TUV) radiation model		   =*
!=    Version 4.2							   =*
!=    May 2003  							   =*
!-----------------------------------------------------------------------------*
!= Developed by Sasha Madronich with important contributions from:	   =*
!= Chris Fischer, Siri Flocke, Julia Lee-Taylor, Bernhard Meyer,	   =*
!= Irina Petropavlovskikh,  Xuexi Tie, and Jun Zen.			   =*
!= Special thanks to Knut Stamnes and co-workers for the development of the  =*
!= Discrete Ordinates code, and to Warren Wiscombe and co-workers for the    =*
!= development of the solar zenith angle subroutine. Citations for the many  =*
!= data bases (e.g. extraterrestrial irradiances, molecular spectra) may be  =*
!= found in the data files headers and/or in the subroutines that read them. =*
!=	      To contact the author, write to:  			   =*
!= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
!= send email to:  sasha@ucar.edu  or tuv@acd.ucar.edu  		   =*
!-----------------------------------------------------------------------------*
!= This program is free software; you can redistribute it and/or modify      =*
!= it under the terms of the GNU General Public License as published by the  =*
!= Free Software Foundation;  either version 2 of the license, or (at your   =*
!= option) any later version.						   =*
!= The TUV package is distributed in the hope that it will be useful, but    =*
!= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
!= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
!= License for more details.						   =*
!= To obtain a copy of the GNU General Public License, write to:	   =*
!= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
!-----------------------------------------------------------------------------*
!= Copyright (C) 1994,95,96,97,98,99,2000,01,02,03  University Corporation   =*
!= for Atmospheric Research						   =*
!-----------------------------------------------------------------------------*
!

MODULE ModTuv

  IMPLICIT NONE

  PRIVATE

  !Public Variables
  PUBLIC :: kz, kw, kt, ks, kj, sw, ns !LFR - doing sw and ns visible outside ModTuv
  PUBLIC  :: initialized

  !Public Subroutines
  PUBLIC :: InitTuv
  PUBLIC :: Tuv
  PUBLIC :: slabel
  PUBLIC :: jlabel
  PUBLIC :: nj
  PUBLIC :: tfiles
  PUBLIC :: files
  PUBLIC :: f
  PUBLIC :: xRef
  PUBLIC :: reverse
  PUBLIC :: wl,wu,wc,nw
  public :: wbioStart,wbioEnd
  PUBLIC :: nrad!,narad
  PUBLIC :: get_zlimit

  LOGICAL,PARAMETER :: isMpi=.true.

  ! BROADLY USED PARAMETERS:
  !_________________________________________________
  ! i/o file unit numbers
  INTEGER, PARAMETER :: kout=53
  INTEGER, PARAMETER :: kin=12
  !_________________________________________________
  ! altitude, wavelength, time (or solar zenith angle) grids
  INTEGER, PARAMETER :: kz=125 ! altitude
  INTEGER, PARAMETER :: kw=20 ! wavelength  --- test Luiz Flavio
  !INTEGER, PARAMETER :: kw=650 ! wavelength
  INTEGER, PARAMETER :: kt=100 ! time/sza
  !_________________________________________________
  ! number of weighting functions
  INTEGER, PARAMETER :: ks=60 !  wavelength dependent
  INTEGER, PARAMETER :: kj=80 !  wavelength and altitude dependent
  REAL,PARAMETER :: deltax = 1.0e-4 ! delta for adding points at
                                    ! beginning or end of data grids
  REAL,PARAMETER :: radius=6.371E+3 ! radius of the earth:
  REAL,PARAMETER :: largest=1.0e+36 ! largest number of the machine:
  REAL,PARAMETER :: pzero = +10./largest ! small numbers (positive and negative)
  REAL,PARAMETER :: nzero = -10./largest
  REAL,PARAMETER :: precis = 1.e-7  ! machine precision
  ! More physical constants:
  !_________________________________________________________________
  ! Na = 6.022142E23  mol-1	= Avogadro constant
  ! kb = 1.38065E-23  J K-1	= Boltzmann constant
  ! R  = 8.31447      J mol-1 K-1 = molar gas constant
  ! h  = 6.626068E-34 J s	= Planck constant
  ! c  = 2.99792458E8 m s-1	= speed of light in vacuum
  ! G  = 6.673E-11    m3 kg-1 s-2 = Netwonian constant of gravitation
  ! sb = 5.67040E-8   W m-2 K-4   = Stefan-Boltzmann constant
  !_________________________________________________________________
  ! (1) From NIST Reference on Constants, Units, and Uncertainty
  ! http://physics.nist.gov/cuu/index.html Oct. 2001.
  ! (2) These constants are not assigned to variable names;  in other
  ! words this is not Fortran code, but only a text table for quick
  ! reference.  To use, you must declare a variable name/type and
  ! assign the value to that variable. Or assign as parameter (see
  ! example for pi above).

  INTEGER, PARAMETER :: totReact=115
  INTEGER, PARAMETER :: nmug = 10
  INTEGER, PARAMETER :: maxstr = 100
  INTEGER, PARAMETER :: maxtrm = 100
  INTEGER, PARAMETER :: maxsqt = 1000
  REAL,DIMENSION(17), PARAMETER :: xslod=(/ &
	  6.2180730E-21, 5.8473627E-22, 5.6996334E-22,  &
	  4.5627094E-22, 1.7668250E-22, 1.1178808E-22,  &
	  1.2040544E-22, 4.0994668E-23, 1.8450616E-23,  &
	  1.5639540E-23, 8.7961075E-24, 7.6475608E-24,  &
	  7.6260556E-24, 7.5565696E-24, 7.6334338E-24,  &
	  7.4371992E-24, 7.3642966E-24/)
  !mz = 5 zenith angle 0,20,40,60,80
  !ms = 73 species
  !mp = 5 pol coeff  0,1,2,3,4
  INTEGER,PARAMETER :: mz=5
  INTEGER,PARAMETER :: ms=73
  INTEGER,PARAMETER :: mp=5
  REAL, PARAMETER :: lbar = 206.214
  CHARACTER(LEN=50),PARAMETER :: jlabel(kj) = &! Photolysis coefficients labels
             (/ &
             'O2 -> O + O                                       ', & !1
             'O3 -> O2 + O(1D)                                  ', & !2 R2 CB07
             'O3 -> O2 + O(3P)                                  ', & !3 R3 CB07
             'NO2 -> NO + O(3P)                                 ', & !4 R1 CB07
             'NO3 -> NO + O2                                    ', & !5 R7 CB07
             'NO3 -> NO2 + O(3P)                                ', & !6 R8 CB07
             'N2O5 -> NO3 + NO + O(3P)                          ', & !7
             'N2O5 -> NO3 + NO2                                 ', & !8
             'N2O -> N2 + O(1D)                                 ', & !9
             'HO2 -> OH + O                                     ', & !10
             'H2O2 -> 2 OH                                      ', & !11 R9 CB07
             'HNO2 -> OH + NO                                   ', & !12 R4 CB07
             'HNO3 -> OH + NO2                                  ', & !13 R5 CB07
             'HNO4 -> HO2 + NO2                                 ', & !14
             'CH2O -> H + HCO                                   ', & !15
             'CH2O -> H2 + CO                                   ', & !16 R10 CB07
             'CH3CHO -> CH3 + HCO                               ', & !17
             'CH3CHO -> CH4 + CO                                ', & !18
             'CH3CHO -> CH3CO + H                               ', & !19
             'C2H5CHO -> C2H5 + HCO                             ', &
             'CHOCHO -> HCO + HCO                               ', &
             'CHOCHO -> CH2O + CO                               ', &
             'CH3COCHO -> CH3CO + HCO                           ', &
             'CH3COCH3 -> CH3CO + CH3                           ', &
             'CH3OOH -> CH3O + OH                               ', &
             'CH3ONO2 -> CH3O + NO2                             ', &
             'CH3CO(OONO2) -> Products                          ', &
             'ClOO -> Products                                  ', &
             'ClONO2 -> Cl + NO3                                ', &
             'ClONO2 -> ClO + NO2                               ', &
             'CH3Cl -> Products                                 ', &
             'CCl2O -> Products                                 ', &
             'CCl4 -> Products                                  ', &
             'CClFO -> Products                                 ', &
             'CF2O -> Products                                  ', &
             'CF2ClCFCl2 (CFC-113) -> Products                  ', &
             'CF2ClCF2Cl (CFC-114) -> Products                  ', &
             'CF3CF2Cl (CFC-115) -> Products                    ', &
             'CCl3F (CFC-11) -> Products                        ', &
             'CCl2F2 (CFC-12) -> Products                       ', &
             'CH3CCl3 -> Products                               ', &
             'CF3CHCl2 (HCFC-123) -> Products                   ', &
             'CF3CHFCl (HCFC-124) -> Products                   ', &
             'CH3CFCl2 (HCFC-141b) -> Products                  ', &
             'CH3CF2Cl (HCFC-142b) -> Products                  ', &
             'CF3CF2CHCl2 (HCFC-225ca) -> Products              ', &
             'CF2ClCF2CHFCl (HCFC-225cb) -> Products            ', &
             'CHClF2 (HCFC-22) -> Products                      ', &
             'BrONO2 -> BrO + NO2                               ', &
             'BrONO2 -> Br + NO3                                ', &
             'CH3Br -> Products                                 ', &
             'CHBr3 -> Products                                 ', &
             'CF3Br (Halon-1301) -> Products                    ', &
             'CF2BrCF2Br (Halon-2402) -> Products               ', &
             'CF2Br2 (Halon-1202) -> Products                   ', &
             'CF2BrCl (Halon-1211) -> Products                  ', &
             'Cl2 -> Cl + Cl                                    ', &
             'CH2(OH)CHO -> Products                            ', &
             'CH3COCOCH3 -> Products                            ', &
             'CH3COCHCH2 -> Products                            ', &
             'CH2C(CH3)CHO -> Products                          ', &
             'CH3COCO(OH) -> Products                           ', &
             'CH3CH2ONO2 -> CH3CH2O + NO2                       ', &
             'CH3CHONO2CH3 -> CH3CHOCH3 + NO2                   ', &
             'CH2(OH)CH2(ONO2) -> CH2(OH)CH2(O.) + NO2          ', &
             'CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2              ', &
             'C(CH3)3(ONO2) -> C(CH3)3(O.) + NO2                ', &
             'ClOOCl -> Cl + ClOO                               ', &
             'CH2(OH)COCH3 -> CH3CO + CH2(OH)                   ', &
             'CH2(OH)COCH3 -> CH2(OH)CO + CH3                   ', &
             'HOBr -> OH + Br                                   ', &
             'BrO -> Br + O                                     ', &
             'Br2 -> Br + Br                                    ', &
             '                                                  ', &
             '                                                  ', &
             '                                                  ', &
             '                                                  ', &
             '                                                  ', &
             '                                                  ', &
             '                                                  '  &
              /)

  !_______________________________________________________________________
  ! select desired extra-terrestrial solar irradiance, using msun:
  !  1 =   extsol.flx:  De Luisi, JGR 80, 345-354, 1975
  !                     280-400 nm, 1 nm steps.
  !  2 =   lowsun3.flx:  Lowtran (John Bahr, priv. comm.)
  !                      173.974-500000 nm, ca. 0.1 nm steps in UV-B
  !  3 =   modtran1.flx:  Modtran (Gail Anderson, priv. comm.)
  !                       200.55-949.40, 0.05 nm steps
  !  4 =   nicolarv.flx:  wvl<300 nm from Nicolet, Plan. Sp. Sci., 29,  951-974, 1981.
  !                       wvl>300 nm supplied by Thekaekera, Arvesen Applied Optics 8,
  !                       11, 2215-2232, 1969 (also see Thekaekera, Applied Optics, 13,
  !                       3, 518, 1974) but with corrections recommended by:
  !                       Nicolet, Plan. Sp. Sci., 37, 1249-1289, 1989.
  !                       270.0-299.0 nm in 0.5 nm steps
  !                       299.6-340.0 nm in ca. 0.4 nm steps
  !                       340.0-380.0 nm in ca. 0.2 nm steps
  !                       380.0-470.0 nm in ca. 0.1 nm steps
  !  5 =  solstice.flx:  From:   MX%"ROTTMAN@virgo.hao.ucar.edu" 12-OCT-1994 13:03:01.62
  !                      Original data gave Wavelength in vacuum
  !                      (Converted to wavelength in air using Pendorf, 1967, J. Opt. Soc. Am.)
  !                      279.5 to 420 nm, 0.24 nm spectral resolution, approx 0.07 nm steps
  !  6 =  suntoms.flx: (from TOMS CD-ROM).  280-340 nm, 0.05 nm steps.
  !  7 =  neckel.flx:  H.Neckel and D.Labs, "The Solar Radiation Between 3300 and 12500 A",
  !                    Solar Physics v.90, pp.205-258 (1984).
  !                    1 nm between 330.5 and 529.5 nm
  !                    2 nm between 631.0 and 709.0 nm
  !                    5 nm between 872.5 and 1247.4 nm
  !                    Units: must convert to W m-2 nm-1 from photons cm-2 s-1 nm-1
  !  8 =  atlas3.flx:  ATLAS3-SUSIM 13 Nov 94 high resolution (0.15 nm FWHM)
  !                    available by ftp from susim.nrl.navy.mil
  !                    atlas3_1994_317_a.dat, downloaded 30 Sept 98.
  !                    150-407.95 nm, in 0.05 nm steps
  !                    (old version from Dianne Prinz through Jim Slusser)
  !                    orig wavelengths in vac, correct here to air.
  !  9 =  solstice.flx:  solstice 1991-1996, average
  !                    119.5-420.5 nm in 1 nm steps

  ! 10 =  susim_hi.flx:  SUSIM SL2 high resolution
  !                      120.5-400.0 in 0.05 nm intervals (0.15 nm resolution)
  ! 11 =  wmo85.flx: from WMO 1995 Ozone Atmospheric Ozone (report no. 16)
  !                  on variable-size bins.  Original values are per bin, not
  !                  per nm.
  ! 12 = combine susim_hi.flx for .lt. 350 nm, neckel.flx for .gt. 350 nm.

  ! 13 = combine
  !     for wl(iw) .lt. 150.01                                susim_hi.flx
  !     for wl(iw) .ge. 150.01 and wl(iw) .le. 400            atlas3.flx
  !     for wl(iw) .gt. 400                                   Neckel & Labs
  INTEGER,PARAMETER :: msun = 13
  INTEGER,PARAMETER :: nTotalFiles=145
  ! quantum yield recommendation:
  !    kjpl87:  JPL recommendation 1987                - JPL 87, 90, 92 do not "tail"
  !    kjpl92:  JPL recommendations 1990/92 (identical) - still with no "tail"
  !    kjpl97:  JPL recommendation 1997, includes tail, similar to Shetter et al.
  !    kmich :  Michelsen et al., 1994
  !    kshet :  Shetter et al., 1996
  !    kjpl00:  JPL 2000
  !    kmats:  Matsumi et al., 2002
  INTEGER,PARAMETER :: kmich = 1
  INTEGER,PARAMETER :: kjpl87 = 2
  INTEGER,PARAMETER :: kjpl92 = 3
  INTEGER,PARAMETER :: kshet = 4
  INTEGER,PARAMETER :: kjpl97 = 5
  INTEGER,PARAMETER :: kjpl00 = 6
  INTEGER,PARAMETER :: kmats = 7
  INTEGER, DIMENSION(20), PARAMETER :: &
           mOption=(/1,7,1,2,6,1,3,1,1,1,4,1,2,5,1,4,2,2,2,0/)
   INTEGER,PARAMETER,DIMENSION(5) :: indSza=(/1,2,3,4,5/)
   REAL,PARAMETER,DIMENSION(5) :: angle=(/0.0,20.0,40.0,60.0,80.0/)
   REAL,PARAMETER,DIMENSION(5) :: fatSum=(/-0.2,-0.2,0.0,0.0,0.0/)
   REAL,PARAMETER,DIMENSION(5) :: ca0= &
                       (/ 4.52372, 4.52372, 4.99378, 0.969867, 1.07801/)
   REAL,PARAMETER,DIMENSION(5) :: ca1= &
                       (/-5.94317,-5.94317,-7.92752,-0.841035,-2.39580/)
   REAL,PARAMETER,DIMENSION(5) :: ca2= &
                       (/ 2.63156, 2.63156, 3.94715, 0.878835, 2.32632/)
   REAL,PARAMETER,DIMENSION(5) :: cb0= &
                       (/ 2.43360, 2.43360, 3.98265, 3.49843,  3.06312/)
   REAL,PARAMETER,DIMENSION(5) :: cb1= &
                       (/-3.61363,-3.61363,-6.90516,-5.98839, -5.26281/)
   REAL,PARAMETER,DIMENSION(5) :: cb2= &
                      (/ 2.19018, 2.19018, 3.93602, 3.50262,  3.20980/)

  !LFR - TUV problem with ozone column
  ! In original code the column of Ozone was obtained from O3 climatology
  ! over USA.
  ! Thefore for the new column (dynamic) from the CCATT model we have
  ! some negative values for Photolisys rate. The cause is the ammount
  ! of ozone at each level. The adjO3 is an attenuation factor for each
  ! reaction
  REAL, PARAMETER :: adjO3(kj)=(/ &
      1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0, & !00
      1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0, & !01
      1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0, & !02
      1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0, & !03
      1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0, & !04
      1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0, & !05
      1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0, & !06
      1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0  & !07
                               /)
  !
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: nrad
!  INTEGER :: narad
  REAL :: s226(kw), s263(kw), s298(kw)
  REAL :: o3xs(kz,kw)
  REAL :: s218(kw), s228(kw), s243(kw), s295(kw)
  REAL :: c0(kw), c1(kw), c2(kw)
  REAL :: yg(kw,totReact), yg1(kw,totReact), yg1n(kw,totReact), yg2(kw,totReact)
  REAL :: yg3(kw,totReact), yg4(kw,totReact), yg5(kw,totReact)
  REAL :: coeff(4,3,totReact)
  REAL :: c( maxtrm )
  REAL :: hugeVar, powmax, powmin, tinyVar
  INTEGER :: maxmsg, nummsg
  LOGICAL :: msglim
  INTEGER :: maxmsg2, nummsg2

  REAL :: tbar(totReact)

  !-------------------------------------------------------
  ! XS_COEFFS containing the Chebyshev
  ! polynomial coeffs necessary to calculate O2 effective
  ! cross-sections
  !-------------------------------------------------------
  REAL*8 :: ac(20,17)
  REAL*8 :: bc(20,17) ! Chebyshev polynomial coeffs
  REAL*4 :: wave_num(17)

  LOGICAL :: initialized
  LOGICAL :: isread
  LOGICAL :: firstcall
  LOGICAL :: call1
  LOGICAL :: pass1
  LOGICAL :: pass2
  LOGICAL :: pass3
  LOGICAL :: pass4
  LOGICAL :: pass5
  LOGICAL :: pass6
  LOGICAL :: doinit
  LOGICAL :: doinit2


  INTEGER :: ila
  INTEGER :: isrb

  REAL :: dither

  REAL :: sqt( maxsqt )

  REAL ::  pi
  REAL ::  twopi
  REAL ::  rpd

  DOUBLE PRECISION :: tol
  DOUBLE PRECISION :: dmach(4)
  REAL :: rmach(4)
  REAL :: gmu( nmug ), gwt( nmug ), ylmg( 0:maxstr, nmug )


  REAL :: coef(mz,ms,mp)  ! pol coefficent for FTUV

  ! Wavelength grid:
  INTEGER :: nw,ns,nj,iw
  REAL :: wl(kw), wc(kw), wu(kw)
  REAL :: f(kw) !extra terrestrial solar flux

  ! O2 absorption cross section
  REAL :: o2xs(kz,kw), o2xs1(kw)

  ! SO2 absorption cross section
  REAL :: so2xs(kw)

  ! NO2 absorption cross section
  REAL :: no2xs(kw)


  REAL :: mm_o3xs(kw)

  integer :: wbioStart,wbioEnd
  REAL :: sw(ks,kw)
  CHARACTER(LEN=50) :: slabel(ks) ! Spectral weighting functions labels

  !REAL              :: sj(kj,kz,kw)
  INTEGER           :: xRef(kj)
  LOGICAL           :: doReaction(kj)
  INTEGER :: ii,jj !Temporary
    ! Parameters for shifting wavelengths air <--> vacuum
  INTEGER           :: mrefr
  LOGICAL           :: lrefr
  REAL              :: airout

  !Files
  TYPE tFiles
     CHARACTER(LEN=200) :: fileName
  END TYPE tFiles
  TYPE(tFiles),ALLOCATABLE,DIMENSION(:) :: files


  REAL,ALLOCATABLE,DIMENSION(:) :: x , y
  REAL,ALLOCATABLE,DIMENSION(:) :: x1, y1
  REAL,ALLOCATABLE,DIMENSION(:) :: x2, y2
  REAL,ALLOCATABLE,DIMENSION(:) :: x3, y3
  REAL,ALLOCATABLE,DIMENSION(:) :: x4, y4
  REAL,ALLOCATABLE,DIMENSION(:) :: x5, y5

  INTEGER :: ierr

   INTEGER :: noPr
   INTEGER :: zLimit
!LFR   REAL :: rgasog,deltap

 CONTAINS


   SUBROUTINE InitTuv(wstart,wstop,nwint,filesName,myNum,chemical_mechanism)

    INCLUDE 'mpif.h'

    INTEGER, INTENT(IN)  :: nwint
    INTEGER, INTENT(IN) :: myNum
    REAL, INTENT(IN)   :: wstart
    REAL, INTENT(IN)   :: wstop
    CHARACTER(LEN=*),INTENT(IN) :: filesName
    CHARACTER(LEN=200) :: filesHome
    CHARACTER(LEN=2) :: fname
    CHARACTER(LEN=*),INTENT(IN) :: chemical_mechanism

    INTEGER :: i,nFiles,nOfFile,iang,js,jp

   IF(myNum==0) THEN !Initialization just for myNum=0

   IF(initialized) THEN
      PRINT *,'ERROR: TUV already initialized!'
      PRINT *,'### Please, check your code ###'
      CALL flush(6)
      STOP
   END IF
   filesHome=''
   IF(LEN(TRIM(filesName))==0) THEN
      !If is to use de deafault files name (test case)
      ALLOCATE (files(nTotalFiles))
      !Default values of files names
      files(  1)%fileName='input/POL.out'
      files(  2)%fileName='datae1/sun/extsol.flx'
      files(  3)%fileName='datae1/sun/lowsun3.flx'
      files(  4)%fileName='datae1/sun/modtran1.flx'
      files(  5)%fileName='datae1/sun/nicolarv.flx'
      files(  6)%fileName='datae2/sun/solstice.flx'
      files(  7)%fileName='datae2/sun/suntoms.flx'
      files(  8)%fileName='datae1/sun/neckel.flx'
      files(  9)%fileName='datae1/sun/atlas3_1994_317_a.dat'
      files( 10)%fileName=files( 6)%fileName
      files( 11)%fileName='datae1/sun/susim_hi.flx'
      files( 12)%fileName='datae1/sun/wmo85.flx'
      files( 13)%fileName=files( 8)%fileName
      files( 14)%fileName=files( 9)%fileName
      files( 15)%fileName=files( 8)%fileName
      files( 16)%fileName='datae1/o2/O2_brasseur.abs'
      files( 17)%fileName='datae1/o2/O2_yoshino.abs'
      files( 18)%fileName='datae1/so2/SO2xs.all'
      files( 19)%fileName='datae1/no2/NO2_ncar_00.abs'
      files( 20)%fileName='datas1/rbm.501'
      files( 21)%fileName='datas1/dna.setlow.new'
      files( 22)%fileName='datas1/SCUP-h'
      files( 23)%fileName='datas1/ery.anders'
      files( 24)%fileName='datas1/acgih.1992'
      files( 25)%fileName='datas1/phaeo.bio'
      files( 26)%fileName='datas1/proro.bio'
      files( 27)%fileName='datas1/cataract_oriowo'
      files( 28)%fileName='dataj1/yld/O3.param_jpl97.yld'
      files( 29)%fileName='dataj1/yld/O3.param.yld'
      files( 30)%fileName='dataj1/yld/O3_shetter.yld'
      files( 31)%fileName='datae1/no2/NO2_jpl94.abs'
      files( 32)%fileName='datae1/no2/NO2_Har.abs'
      files( 33)%fileName='dataj1/yld/NO2_calvert.yld'
      files( 34)%fileName='dataj1/abs/NO3_gj78.abs'
      files( 35)%fileName='dataj1/abs/NO3_jpl94.abs'
      files( 36)%fileName='dataj1/abs/N2O5_jpl97.abs'
      files( 37)%fileName='dataj1/abs/HNO2_jpl92.abs'
      files( 38)%fileName='dataj1/abs/HNO3_burk.abs'
      files( 39)%fileName='dataj1/abs/HNO4_jpl92.abs'
      files( 40)%fileName='dataj1/abs/H2O2_jpl94.abs'
      files( 41)%fileName='dataj1/abs/CHBr3.abs'
      files( 42)%fileName='dataj1/abs/CHBr3.jpl97'
      files( 43)%fileName='dataj1/ch2o/CH2O_nbs.abs'
      files( 44)%fileName='dataj1/CH2O_iupac1.abs'
      files( 45)%fileName='dataj1/ch2o/CH2O_can_hr.abs'
      files( 46)%fileName='dataj1/ch2o/CH2O_can_lr.abs'
      files( 47)%fileName='dataj1/ch2o/CH2O_rog.abs'
      files( 48)%fileName='dataj1/ch2o/CH2O_ncar.abs'
      files( 49)%fileName='dataj1/ch2o/CH2O_i_mad.yld'
      files( 50)%fileName='dataj1/ch2o/CH2O_ii_mad.yld'
      files( 51)%fileName='dataj1/ch2o/CH2O_iupac.yld'
      files( 52)%fileName='dataj1/ch2o/CH2O_jpl97.dat'
      files( 53)%fileName='dataj1/ch3cho/CH3CHO_iup.abs'
      files( 54)%fileName='dataj1/ch3cho/d021_cp.abs'
      files( 55)%fileName='dataj1/ch3cho/CH3CHO_mar.abs'
      files( 56)%fileName='dataj2/kfa/ch3cho.005'
      files( 57)%fileName='dataj1/ch3cho/CH3CHO_iup.yld'
      files( 58)%fileName='dataj1/ch3cho/d021_i.yld'
      files( 59)%fileName='dataj1/ch3cho/d021_ii.yld'
      files( 60)%fileName='dataj1/ch3cho/d021_iii.yld'
      files( 61)%fileName='dataj1/ch3cho/CH3CHO_press.yld'
      files( 62)%fileName='dataj1/c2h5cho/C2H5CHO_iup.abs'
      files( 63)%fileName='dataj2/kfa/c2h5cho.001'
      files(64)%fileName='dataj1/c2h5cho/C2H5CHO_iup.yld'
      files( 65)%fileName='dataj1/chocho/CHOCHO_iup.abs'
      files( 66)%fileName='dataj2/kfa/chocho.001'
      files( 67)%fileName='dataj1/chocho/glyoxal_orl.abs'
      files( 68)%fileName='dataj1/chocho/glyoxal_horowitz.abs'
      files( 69)%fileName='dataj1/ch3cocho/CH3COCHO_iup1.abs'
      files( 70)%fileName='dataj1/ch3cocho/CH3COCHO_iup2.abs'
      files( 71)%fileName='dataj1/ch3cocho/CH3COCHO_ncar.abs'
      files( 72)%fileName='dataj2/kfa/ch3cocho.001'
      files( 73)%fileName='dataj2/kfa/ch3cocho.002'
      files( 74)%fileName='dataj2/kfa/ch3cocho.003'
      files( 75)%fileName='dataj2/kfa/ch3cocho.004'
      files( 76)%fileName='dataj1/ch3cococh3/biacetyl_plum.abs'
      files( 77)%fileName='dataj1/chocho/glyoxal_orl.abs'
      files( 78)%fileName='dataj1/ch3cocho/CH3COCHO_km.yld'
      files( 79)%fileName='dataj1/ch3coch3/CH3COCH3_cp.abs'
      files( 80)%fileName='dataj1/ch3coch3/CH3COCH3_iup.abs'
      files( 81)%fileName='dataj1/ch3coch3/CH3COCH3_noaa.abs'
      files( 82)%fileName='dataj1/ch3coch3/CH3COCH3_iup.yld'
      files( 83)%fileName='dataj1/ch3ooh/CH3OOH_jpl94.abs'
      files( 84)%fileName='dataj1/ch3ooh/CH3OOH_iup.abs'
      files( 85)%fileName='dataj1/ch3ooh/CH3OOH_ct.abs'
      files( 86)%fileName='dataj1/ch3ooh/CH3OOH_ma.abs'
      files( 87)%fileName='dataj1/rono2/CH3ONO2_cp.abs'
      files( 88)%fileName='dataj1/rono2/CH3ONO2_tal.abs'
      files( 89)%fileName='dataj1/rono2/CH3ONO2_iup1.abs'
      files( 90)%fileName='dataj1/rono2/CH3ONO2_iup2.abs'
      files( 91)%fileName='dataj1/rono2/CH3ONO2_tay.abs'
      files( 92)%fileName='dataj1/rono2/CH3ONO2_rat.abs'
      files( 93)%fileName='dataj1/rono2/CH3ONO2_lib.abs'
      files( 94)%fileName='dataj1/rono2/PAN_talukdar.abs'
      files( 95)%fileName='dataj1/abs/CCl2O_jpl94.abs'
      files( 96)%fileName='dataj1/abs/CCl4_jpl94.abs'
      files( 97)%fileName='dataj1/abs/CClFO_jpl94.abs'
      files( 98)%fileName='dataj1/abs/CF2O_jpl94.abs'
      files( 99)%fileName='dataj1/abs/CFC-113_jpl94.abs'
      files(100)%fileName='dataj1/abs/CFC-114_jpl94.abs'
      files(101)%fileName='dataj1/abs/CFC-115_jpl94.abs'
      files(102)%fileName='dataj1/abs/CFC-11_jpl94.abs'
      files(103)%fileName='dataj1/abs/CFC-12_jpl94.abs'
      files(104)%fileName='dataj1/abs/CH3Br_jpl94.abs'
      files(105)%fileName='dataj1/abs/CH3CCl3_jpl94.abs'
      files(106)%fileName='dataj1/abs/CH3Cl_jpl94.abs'
      files(107)%fileName='dataj1/abs/ClOO_jpl94.abs'
      files(108)%fileName='dataj1/abs/HCFCs_orl.abs'
      files(109)%fileName='dataj1/abs/HCFCs_orl.abs'
      files(110)%fileName='dataj1/abs/HCFC-141b_jpl94.abs'
      files(111)%fileName='dataj1/abs/HCFCs_orl.abs'
      files(112)%fileName='dataj1/abs/HCFC-225ca_jpl94.abs'
      files(113)%fileName='dataj1/abs/HCFC-225cb_jpl94.abs'
      files(114)%fileName='dataj1/abs/HCFC-22_jpl94.abs'
      files(115)%fileName='dataj1/abs/HO2_jpl94.abs'
      files(116)%fileName='dataj1/abs/Halon-1202_jpl97.abs'
      files(117)%fileName='dataj1/abs/Halon-1211_jpl97.abs'
      files(118)%fileName='dataj1/abs/Halon-1301_jpl97.abs'
      files(119)%fileName='dataj1/abs/Halon-2402_jpl97.abs'
      files(120)%fileName='dataj1/abs/ClONO2_jpl97.abs'
      files(121)%fileName='dataj1/abs/BrONO2_jpl03.abs'
      files(122)%fileName='dataj1/abs/CL2_fpp.abs'
      files(123)%fileName='dataj1/ch2ohcho/glycolaldehyde.abs'
      files(124)%fileName='dataj1/ch3cococh3/biacetyl_plum.abs'
      files(125)%fileName='dataj1/ch3cococh3/biacetyl_horowitz.abs'
      files(126)%fileName='dataj1/abs/methylvinylketone.abs'
      files(127)%fileName='dataj1/abs/methacrolein.abs'
      files(128)%fileName='dataj1/ch3cocooh/pyruvic_horowitz.abs'
      files(129)%fileName='dataj1/rono2/RONO2_talukdar.abs'
      files(130)%fileName='dataj1/rono2/RONO2_talukdar.abs'
      files(131)%fileName='dataj1/abs/CLOOCL_jpl02.abs'
      files(132)%fileName='dataj1/abs/Hydroxyacetone.abs'
      files(133)%fileName='dataj1/abs/BrO.jpl03'
      files(134)%fileName='dataj1/abs/Br2.abs'
      files(135)%fileName='output/out_11'
      files(136)%fileName='datae1/grids/combined.grid'
      files(137)%fileName='datae1/grids/fast_tuv.grid'
      files(138)%fileName='datae1/o2/effxstex.txt'
      files(139)%fileName='datae1/wmo85'
      files(140)%fileName='datae1/o3/O3.molina.abs'
      files(141)%fileName='datae1/wmo85'
      files(142)%fileName='datae1/o3/o3absqs.dat'
      files(143)%fileName='datae1/wmo85'
      files(144)%fileName='datae1/o3/O3_bass.abs'
!srf  files(145)%filename='input/xreference.dat'
      files(145)%filename='input/'//trim(chemical_mechanism)//'xreference.dat'
    ELSE
       PRINT *, 'TUV: Reading list of input files'

       !begin RMF
       ! filesname already contains a relative directory link to
       ! ./tables/tuvData/listFiles.dat
       ! this was changed to ./tables/tuvData/ and filesHome var
       ! can be defined locally.
       !end RMF

       !Read all files name from filesName file.
       OPEN(kin, FILE=trim(filesName//'listFiles.dat'))
       READ(kin,*) !Header line
       READ(kin,FMT='(I3.3)') nFiles
       PRINT *,'TUV: Number of files: ',nFiles
       READ(kin,*) !Header line
       READ(kin,FMT='(A)') filesHome

       filesHome = trim(filesName)

       PRINT *,'TUV: Home directory of files: ',trim(filesHome)

       READ(kin,*) !Header line
       ALLOCATE (files(nFiles))
       DO i=1,nFiles
          READ(kin,FMT='(I3,1X,A50)') nOfFile,files(i)%fileName
	  PRINT *,i,trim(files(i)%fileName);CALL flush(6)
       END DO
       !-srf - special treatment for xreference.dat
       !-      since it depends on the chemical mechanism
       DO i=1,nFiles
        if(files(i)%fileName == 'input/xreference.dat') then
          files(i)%fileName = 'input/'//trim(chemical_mechanism)//'xreference.dat'
	  print*,'xreference file=',files(i)%fileName
        endif
       END DO
       !-srf- end
    END IF

    !PRINT *,'TUV: puting prefix in filename - ',trim(filesHome) ;CALL flush(6)
    !Putting the prefix of files
    DO i=1,nFiles
       files(i)%fileName=trim(filesHome)//trim(files(i)%fileName)
    END DO

    pi   = 2.*ASIN( 1.0 )
    twopi	= 2.*pi
    rpd   = pi/180.0

    ! wavelengths (creates wavelength grid: lower, center, upper of each bin)
    ! NOTE:  Wavelengths are in vacuum.  To use wavelengths in air, see
    ! Section 3 below, where you must set lrefr= .TRUE.
    !PRINT *,'TUV: creates wavelength grid';CALL flush(6)
    CALL gridw(wstart, wstop, nwint)

    !PRINT *,'TUV: reading pol coef.';CALL flush(6)
    CALL readpol(coef)


    !PRINT *,'TUV: reading all files';CALL flush(6)
    CALL ReadAll(nw,wl)


    !**** Correction for air-vacuum wavelength shift:
    ! The TUV code assumes that all working wavelengths are strictly IN-VACUUM. This is assumed for ALL
    ! spectral data including extraterrestrial fluxes, ozone (and other) absorption cross sections,
    ! and various weighting functons (action spectra, photolysis cross sections, instrument spectral
    ! response functions).

    !  Occasionally, users may want their results to be given for wavelengths measured IN-AIR.
    ! The shift between IN-VACUUM and IN-AIR wavelengths depends on the index of refraction
    ! of air, which in turn depends on the local density of air, which in turn depends on
    ! altitude, temperature, etc.
    !  Here, we provide users with the option to use a wavelength grid IN-AIR, at the air density
    ! corresponding to the output altitude, airout = airden(izout), by setting the logical variable
    ! lrefr = .TRUE.  (default is lrefr = .FALSE.).  The wavelengths specified in gridw.f will be assumed
    ! to be IN-AIR, and will be shifted here to IN-VACUUM values to carry out the calculatons.
    ! The actual radiative transfer calculations will be done strictly with IN-VACUUM values.
    ! If this shift is applied (i.e., if lrefr = .TRUE.), the wavelength grid will be shifted back to air
    ! values just before the output is written.
    !  Note:  if this option is used (lref = .TRUE.), the wavelength values will be correct ONLY at the
    ! selected altitude, iz = iout.  The wavelength shift will be INCORRECT at all other altitudes.
    !  Note:  This option cannot be changed interactively in the input table.  It must be changed here.

    ! ___ SECTION 4: READ SPECTRAL DATA ____________________________
    ! read (and grid) extra terrestrial flux data:

    !PRINT *,'TUV: reading (and grid) extra terrestrial flux data';CALL flush(6)
    CALL rdetfl(nw,wl)

    ! read cross section data for
    !  O2 (will overwrite at Lyman-alpha and SRB wavelengths
    !	       see subroutine la_srb.f)
    !  O3 (temperature-dependent)
    !  SO2
    !  NO2

    !PRINT *,'TUV: reading cross section data for O2, O3, SO2 and NO2';CALL flush(6)
    CALL rdo2xs(nw,wl, o2xs1)
    CALL rdso2xs(nw,wl, so2xs)
    CALL rdno2xs(nw,wl, no2xs)

    !***** Spectral weighting functions
    ! (Some of these depend on temperature T and pressure P, and therefore
    !  on altitude z.  Therefore they are computed only after the T and P profiles
    !  are set above with subroutines settmp and setair.)
    ! Photo-physical   set in swphys.f (transmission functions)
    ! Photo-biological set in swbiol.f (action spectra)
    ! Photo-chemical   set in swchem.f (cross sections x quantum yields)
    ! Physical and biological weigthing functions are assumed to depend
    ! only on wavelength.
    ! Chemical weighting functions (product of cross-section x quantum yield)
    ! for many photolysis reactions are known to depend on temperature
    ! and/or pressure, and therefore are functions of wavelength and altitude.
    ! Output:
    ! from pphys & pbiol:  s(ks,kw) - for each weighting function slabel(ks)
    ! from pchem:  sj(kj,kz,kw) - for each reaction jlabel(kj)
    ! For pchem, need to know temperature and pressure profiles.


    !PRINT *,'TUV:reading physical spectral weighting functions';CALL flush(6)
    CALL swphys(nw,wl,wc, ns,sw,slabel)
    !PRINT *,'TUV:reading biological spectral weighting functions';CALL flush(6)
    CALL swbiol(nw,wl,wc, ns,sw,slabel)

    initialized=.true.

    isread = .false.
    firstcall= .true.
    call1=.true.
    pass1=.true.
    pass2=.true.
    pass3=.true.
    pass4=.true.
    pass5=.true.
    pass5=.true.
    doinit=.true.
    doinit2=.true.

    tol  = 10.*d1mach( 4 )
    maxmsg=100
    nummsg=0
    msglim=.false.
    maxmsg2=50
    nummsg2=0

    !------------------------------------------
    !      Loads Chebyshev polynomial Coeff.
    !------------------------------------------
    CALL init_xs

    dither = 10.*r1mach( 4 )
    !** Must dither more on Cray (14-digit prec)
    IF( dither < 1.e-10 ) dither = 10.*dither

    END IF !(just for processor myNum=0)

    IF (isMpi) THEN
        !PRINT *,'Sending broadcast from initialization to all processors'
        CALL flush(6)
        !MPI broadcast
        !CALL MPI_BCAST(sj,kj*kz*kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(xRef,kj, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(doReaction,kj, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(mrefr,1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(lrefr,1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(airout,1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(coef,mz*ms*mp, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(wl,kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(wc,kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(wu,kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(f,kw,  MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(o2xs,kz*kw,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(o2xs1,kw,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(so2xs,kw,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(no2xs,kw,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(mm_o3xs,kw,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(sqt,maxsqt,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(pi,1 ,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(twopi,1,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(rpd,1,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(tol,1,MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(dmach,4,MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(rmach,4,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(gmu,nmug,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(gwt,nmug,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(ylmg,(maxstr+1)*nmug,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(s226,kw,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(s263,kw,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(s298,kw,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(o3xs,kz*kw,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(s218,kw,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(s228,kw,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(s243,kw,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(s295,kw,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(c0,kw,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(c1,kw,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(c2,kw,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(yg,kw*totReact,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(yg1,kw*totReact,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(yg1n,kw*totReact,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(yg2,kw*totReact,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(yg3,kw*totReact,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(yg4,kw*totReact,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(yg5,kw*totReact,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(coeff,4*3*totReact,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(c,maxtrm,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(hugeVar,1,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(powmax,1,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(powmin,1,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(tinyVar,1,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(maxmsg,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(nummsg,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(msglim,1,MPI_LOGICAL,0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(maxmsg2,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(nummsg2,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(tbar,totReact,MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(ac,20*17,MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(bc,20*17,MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(nw,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(ns,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(nj,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(iw,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call mpi_bcast(sw,ks*kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(initialized,1,MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
     END IF
    PRINT *,'-------------------------------------------------';CALL flush(6)

  END SUBROUTINE InitTuv

  ! Radiative transfer scheme:
  !   nstr = number of streams
  !	   If nstr < 2, will use 2-stream Delta Eddington
  !	   If nstr > 1, will use nstr-stream discrete ordinates
  ! Location (geographic):
  !   lat = LATITUDE (degrees, North = positive)
  !   lon = LONGITUDE (degrees, East = positive)
  !   esfact = 1. (Earth-sun distance = 1.000 AU)
  ! Vertical grid:
  !   zstart = surface elevation above sea level, km
  !   zstop = top of the atmosphere (exospheric), km
  !   nz = number of vertical levels, equally spaced
  !	 (nz will increase by +1 if zout does not match altitude grid)
  ! Wavlength grid:
  !   wstart = starting wavelength, nm
  !   wstop  = final wavelength, nm
  !   nwint = number of wavelength intervals, equally spaced
  !	    if nwint < 0, the standard atmospheric wavelength grid, not
  !	    equally spaced, from 120 to 735 nm, will be used. In this
  !	    case, wstart and wstop values are ignored.
  ! Surface condition:
  !   alsurf = surface albedo, wavelength independent
  ! Column amounts of absorbers (in Dobson Units, from surface to space):
  !	   Vertical profile for O3 from USSA76.  For SO2 and NO2, vertical
  !	   concentration profile is 2.69e10 molec cm-3 between 0 and
  !	   1 km above sea level, very small residual (10/largest) above 1 km.
  !   so2col = sulfur dioxide (SO2)
  !   no2col = nitrogen dioxide (NO2)
  ! Cloud, assumed horizontally uniform, total coverage, single scattering
  !	  albedo = 0.9999, asymmetry factor = 0.85, indep. of wavelength,
  !	  and also uniform vertically between zbase and ztop:
  !   taucld = vertical optical depth, independent of wavelength
  !   zbase = altitude of base, km above sea level
  !   ztop = altitude of top, km above sea level
  ! Aerosols, assumed vertical provile typical of continental regions from
  !	  Elterman (1968):
  !   tauaer = aerosol vertical optical depth at 550 nm, from surface to space.
  !	    If negative, will default to Elterman's values (ca. 0.235
  !	    at 550 nm).
  !   ssaaer = single scattering albedo of aerosols, wavelength-independent.
  !   alpha = Angstrom coefficient = exponent for wavelength dependence of
  !	    tauaer, so that  tauaer1/tauaer2  = (w2/w1)**alpha.
  ! Directional components of radiation, weighting factors:
  !   dirsun = direct sun
  !   difdn = down-welling diffuse
  !   difup = up-welling diffuse
  !	 e.g. use:
  !	 dirsun = difdn = 1.0, difup = 0 for total down-welling irradiance
  !	 dirsun = difdn = difup = 1.0 for actinic flux from all directions
  !	 dirsun = difdn = 1.0, difup = -1 for net irradiance

  SUBROUTINE Tuv(mynum,nbl, &
               nstr,nz,maxZ,zLevel,nwint,sza,albedo,so2col,no2col, &
	       dtcld,     & !
               omcld,     & !
               gcld ,     & !
               dtaer,     & !
               omaer,     & !
               gaer ,     & !
	       alpha,dirsun,difdn, &

	       difup,tlev,tlay,airden,cair,co3,tco3,esfact,is2print, &
	       rate,valj,sirrad,saflux)

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: myNum,maxZ
    INTEGER,INTENT(IN) :: nbl !sub grid domain
    INTEGER,INTENT(IN) :: nstr      ! number of radiation streams
    INTEGER,INTENT(IN) :: nz(nbl) ! Altitude grid
    INTEGER,INTENT(IN) :: nwint     ! Wavelength grid: n
    REAL   ,INTENT(IN) :: zLevel(nbl,kz)     !vector of altitude levels (in km)
    REAL   ,INTENT(IN) :: sza(nbl)	    ! Solar zenith angle and azimuth
    REAL   ,INTENT(IN) :: albedo(nbl,kw)    ! surface albedo
    REAL   ,INTENT(IN) :: so2col(nbl)    ! Total columns of SO2 (Dobson Units)
    REAL   ,INTENT(IN) :: no2col(nbl)    ! Total columns of NO2 (Dobson Units)
    REAL   ,INTENT(IN) :: dtcld(nbl,kz,kw)
    REAL   ,INTENT(IN) :: omcld(nbl,kz,kw)
    REAL   ,INTENT(IN) :: gcld(nbl,kz,kw)
    REAL   ,INTENT(IN) :: dtaer(nbl,kz,kw)
    REAL   ,INTENT(IN) :: omaer(nbl,kz,kw)
    REAL   ,INTENT(IN) :: gaer(nbl,kz,kw)
    REAL   ,INTENT(IN) :: alpha     ! Angstrom alpha
    REAL   ,INTENT(IN) :: dirsun    ! direct sun
    REAL   ,INTENT(IN) :: difdn     ! down-welling diffuse
    REAL   ,INTENT(IN) :: difup     ! up-welling diffuse
    REAL   ,INTENT(IN) :: tlev(nbl,kz)  ! temperature (K) at each specified altitude level
    REAL   ,INTENT(IN) :: tlay(nbl,kz)  ! temperature (K) at each specified altitude layer
    REAL   ,INTENT(IN) :: airden(nbl,kz)! air density (molec/cc) at each specified altitude
    REAL   ,INTENT(IN) :: cair(nbl,kz)  ! number of air molecules per cm^2 at each altitude
    REAL   ,INTENT(IN) :: co3(nbl,kz)   ! Total of molec O3 cm-2
    REAL   ,INTENT(IN) :: tco3(nbl,kz)  ! Total column o3 - molec cm-2
    REAL   ,INTENT(IN) :: esfact    ! sun-earth distance UA

    LOGICAL,INTENT(IN) :: is2print	      !To print a debug on each column

    REAL,INTENT(OUT) :: rate(nbl,ks,kz)	      ! Weighted irradiances (dose rates) W m-2
    REAL,INTENT(OUT) :: valj(nbl,kj,kz)	      ! Photolysis coefficients (j-values)
    REAL,INTENT(OUT) :: sirrad(nbl,kz,kw)	      ! Spectral irradiance, [W m-2 nm-1]
    REAL,INTENT(OUT) :: saflux(nbl,kz,kw)	      ! Spectral actinic flux, quanta s-1 nm-1 cm-2

    ! Altitude grid
    INTEGER :: iz, izout !LFR>izout not defined
    ! Extra terrestrial solar flux
    REAL    :: etf(kw)
    INTEGER :: is, iob
    REAL    :: drdw
    INTEGER :: ij,iaux,jaux
    REAL    :: djdw
    INTEGER :: iw
    INTEGER :: iAng
    ! Other user-defined variables here:
    INTEGER :: nn1,nn2  	   !XUEXI
    REAL    :: adjcoe1(kj,kz),adjcoe2(kj,kz) !adjcoe(kj,kz),

    REAL   ,DIMENSION(nbl,kz)       :: edir
    REAL   ,DIMENSION(nbl,kz)       :: edn
    REAL   ,DIMENSION(nbl,kz)       :: eup
    REAL   ,DIMENSION(nbl,kz)       :: fdir
    REAL   ,DIMENSION(nbl,kz)       :: fdn
    REAL   ,DIMENSION(nbl,kz)       :: fup
    REAL   ,DIMENSION(nbl,kz)       :: scol
    REAL   ,DIMENSION(nbl,kz)       :: vcol
    INTEGER,DIMENSION(nbl,0:kz)     :: nid
    REAL   ,DIMENSION(nbl,kz,kw)    :: dtO2
    REAL   ,DIMENSION(nbl,kz,kw)    :: dtO3
    REAL   ,DIMENSION(nbl,kz,kw)    :: o3xs
    REAL   ,DIMENSION(nbl,kj,kz)    :: adjcoe
    REAL   ,DIMENSION(nbl,kz,kw)    :: o2xs
    REAL   ,DIMENSION(nbl,kz,kw)    :: dtRl
    REAL   ,DIMENSION(nbl,kz,kw)    :: dtSo2
    REAL   ,DIMENSION(nbl,kz,kw)    :: dtNo2
    REAL   ,DIMENSION(nbl,0:kz,kz)  :: dsdh
    REAL   ,DIMENSION(nbl,kj,kz,kw) :: sj

    REAL :: valja


    IF (.not. initialized) THEN
       PRINT *, '+------------------------------------------------------+'
       PRINT *, '|            ERROR: Tuv not initialized!               |'
       PRINT *, '|                 Check your driver                    |'
       PRINT *, '|  Please, before call Tuv you must call the InitTuv   |'
       PRINT *, '| Try using: CALL initTuv(wstart,wstop,nwint,nz,izout) |'
       PRINT *, '+------------------------------------------------------+'
       CALL flush(6)
       STOP
    END IF

    noPr=myNum
    IF (nbl == 0) THEN
        zLimit=0
    ELSE
        ! BUG #4002
        ! MAX z must be at maximum nz(nbl) - check
        ! e.g. - tLev(iob, i) could be 0 at maxZ => check sub r01,
        ! line qy1d(iob) = fo3qy2(wc(iw),tLev(iob,i))
        zLimit=min(maxZ, nz(nbl))
    END IF

    ! correction for earth-sun distance
    !default is 1.0 UA - must be set for each season
    DO iw = 1, nw - 1
         etf(iw) = f(iw) * esfact
    END DO

    rate=0.0
    valJ=0.0
    sirRad=0.0
    saFlux=0.0
    adjcoe=1.0

    ! Calculating sj values of the reaction
    CALL swchem(nw,wl,nz,nj,nbl,tLev,airDen,sj)
    !WRITE(76,FMT='(I3.3,A)') Mynum,'Swchem filled';CALL flush(76)
    !Ozone molecular absorption cross section
    CALL rdo3xs(nw,wl,nz,nbl,tLay,o3xs)
    !WRITE(76,FMT='(I3.3,A)') Mynum,'Ozone absortion OK';CALL flush(76)
    ! Rayleigh optical depth increments:
    CALL odrl(nz, nw, wl,nbl,cAir,dtRl)
    ! O2 vertical profile and O2 absorption optical depths
    ! For now, O2 densitiy assumed as 20.95% of air density, can be change
    ! in subroutine.
    ! Optical depths in Lyman-alpha and SRB will be over-written
    ! in subroutine la_srb.f
    DO iz = 1, zLimit
       DO iw =1, nw - 1
    	  DO iob=1,nbl
     	     IF(iz>nz(iob)) CYCLE
    	     dtO2(iob,iz,iw) = 0.2095 * cAir(iob,iz) * o2xs1(iw)
    	  END DO
       END DO
    END DO
    ! Ozone optical depths
    DO iob=1,nbl
       DO   iw = 1, nw-1
    	  DO   iz = 1, nz(iob) - 1
    	    dtO3(iob,iz,iw) = co3(iob,iz) * o3xs(iob,iz,iw)
    	  END DO
       END DO
    END DO
    ! SO2 vertical profile and optical depths
    CALL setso2(nz,nw,so2xs,nbl,zLevel,so2col,cAir,dtSo2)
    ! NO2 vertical profile and optical depths
    CALL setno2(nz,nw,no2xs,nbl,zLevel,no2col,cAir,dtNo2)
    ! slant path lengths for spherical geometry
    CALL sphers(nz,nbl,zLevel,sza,nid,dsdh)
    ! Calculate vertical and slant air columns
    CALL airmas(nz,nbl,cAir,scol,vcol,nid,dsdh)

    DO iob=1,nbl
       !Obtain index for each angle
       iAng=int(sza(iob)/20.0+1)
       !Adjust for each angle
      IF(iang==1) THEN
          CALL setz(nz(iob),nj,coef,adjcoe1,iAng,iob, &
	        tLev(iob,1),tcO3(iob,:))
      	  DO ij=1,kj
     	     DO iz=1,nz(iob)
     	       adjCoe(iob,ij,iz)= adjcoe1(ij,iz)
    	     END DO
    	  END DO
       ELSE IF(iang<5) THEN
          CALL setz(nz(iob),nj,coef,adjcoe1,iAng,iob, &
     		 tLev(iob,1),tcO3(iob,:))
          CALL setz(nz(iob),nj,coef,adjcoe2,iAng+1,iob, &
     		 tLev(iob,1),tcO3(iob,:))
     	  DO ij=1,kj
     	     DO iz=1,nz(iob)
     		adjCoe(iob,ij,iz)= adjcoe1(ij,iz) + &
     		(adjcoe2(ij,iz) - adjcoe1(ij,iz))  &
     		* (sza(iob)- angle(iAng))/(20.0)
    	    END DO
    	  END DO
       ELSEIF(iang==5) THEN
          CALL setz(nz(iob),nj,coef,adjcoe1,iAng,iob, &
     		 tLev(iob,1),tcO3(iob,:))
          IF(abs(sza(iob))<90) THEN
      	     DO ij=1,kj
     	        DO iz=1,nz(iob)
     		   adjCoe(iob,ij,iz)= adjcoe1(ij,iz)
    	        END DO
    	     END DO
          END IF
       END IF
       IF(iAng > 5) THEN
      	  DO ij=1,kj
     	     DO iz=1,nz(iob)
     		adjCoe(iob,ij,iz)= 1.0
    	    END DO
    	  END DO
       END IF
!Debug beg
!      WRITE(45,FMT='(4(A3,1X),5(A20,1X))') 'iob','ij','iz','iAn','sza(iob)', &
!                                          'angle(iang)', 'adjcoe','adjcoe1','adjcoe2'
!      DO ij=1,kj
!     	  DO iz=1,nz(iob)
!	     IF(adjCoe(iob,ij,iz)<0) THEN
!	         WRITE(45,FMT='(4(i3.3,1X),5(F20.10,1X))') &
!		         iob,ij,iz,iAng,sza(iob),angle(iang),adjCoe(iob,ij,iz), &
!	     		 adjcoe1(ij,iz),adjcoe2(ij,iz)
!             END IF
!          END DO
!      END DO
!Debug end
    END DO
    ! Recalculate effective O2 optical depth and cross sections for Lyman-alpha
    ! and Schumann-Runge bands, must know zenith angle
    ! Then assign O2 cross section to sj(1,*,*)
    CALL la_srb(nz,nw,wl,o2xs1,nbl,tLev,dtO2,o2xs,scol,vcol)
    !=  Update the weighting function (cross section x quantum yield) for O2
    !=  photolysis.  The strong spectral variations in the O2 cross sections are
    !=  parameterized into a few bands for Lyman-alpha (121.4-121.9 nm, one band)
    !=  and Schumann-Runge (174.4-205.8, 17 bands) regions. The parameterizations
    !=  depend on the overhead O2 column, and therefore on altitude and solar
    !=  zenith angle, so they need to be updated at each time/zenith step.
    DO iw = 1, nw-1
       DO iob=1,nbl
    	  DO iz = 1, nz(iob)
    	    sj(iob,1,iz,iw) = o2xs(iob,iz,iw)
    	  END DO
       END DO
    END DO

    !Main wavelength loop:
    DO   iw = 1, nw-1
       !* monochromatic radiative transfer. Outputs are:
       !  normalized irradiances edir(iz), edn(iz), eup(iz)
       !  normalized actinic fluxes  fdir(iz), fdn(zi), fup(iz)
       !  where
       !  dir = direct beam, dn = down-welling diffuse,
       !   up = up-welling diffuse
       CALL rtlink(nstr,nz,iw,nbl,sza,albedo, &
     		   dtCld,omCld,gCld,dtAer,omAer,gAer,dtO2,dtO3, &
     		   eDir,eDn,eUp,fDir,fDn,fUp,dtRl,dtSo2,dtNo2, &
     		   nid,dsdh)
       ! Spectral irradiance, W m-2 nm-1, down-welling:
       DO iob=1,nbl
     	  DO iz = 1, nz(iob)
     	     sirRad(iob,iz,iw) = etf(iw) *  &
     			   (dirsun* eDir(iob,iz)+ &
     			   difdn* eDn(iob,iz) + difup*eUp(iob,iz))
     	  END DO
       END DO
       ! Spectral actinic flux, quanta s-1 nm-1 cm-2, all directions:
       !     units conversion:  1.e-4 * (wc*1e-9) / (hc = 6.62E-34 * 2.998E8)
       DO iob=1,nbl
     	  DO iz = 1, nz(iob)
    	     saFlux(iob,iz,iw) = &
     		etf(iw)* 5.039E11 * wc(iw) *  &
     		(dirsun* fdir(iob,iz) + &
     		difdn*fdn(iob,iz) + difup*fup(iob,iz))
     	  END DO
       END DO
       !** Accumulate weighted integrals over wavelength, at all altitudes:
       DO iob=1,nbl
     	  DO  iz = 1, nz(iob)
     	     ! Weighted irradiances (dose rates) W m-2
     	     DO is = 1, ns
     		drdw = sirRad(iob,iz,iw) * sw(is,iw)
     		rate(iob,is,iz) = rate(iob,is,iz) + drdw * (wu(iw)-wl(iw))
     	     END DO
     	  END DO
       END DO
       DO ij = 1, kj
    	  IF(doReaction(ij)) THEN
    	     DO iob=1,nbl
     		DO  iz = 1,nz(iob)
     		  ! Photolysis rate coefficients (J-values) s-1
		   valja=valJ(iob,ij,iz)
     		   djdw = saFlux(iob,iz,iw) * sj(iob,ij,iz,iw) * &
     			  adjCoe(iob,ij,iz)  ! XUEXI
     		   valJ(iob,ij,iz) = valJ(iob,ij,iz) + &
     			  djdw * (wu(iw) - wl(iw))
!LFR> 	      IF(valJ(iob,ij,iz)<0) THEN
!LFR> 		 WRITE (45,FMT='("Valj<0: ",7(I5.5,1X),19(E14.6,1X))') &
!LFR> 			 iob,nbl,ij,kj,iz,nz(iob),iw,valJ(iob,ij,iz),saFlux(iob,iz,iw), &
!LFR> 			 etf(iw),wc(iw),(dirsun* fdir(iob,iz) + &
!LFR> 				 difdn*fdn(iob,iz) + difup*fup(iob,iz)), &
!LFR> 			 dirsun* fdir(iob,iz),difdn*fdn(iob,iz),difup*fup(iob,iz), &
!LFR> 			 dirsun,fdir(iob,iz),difdn,fdn(iob,iz),difup,fup(iob,iz), &
!LFR> 			 sj(iob,ij,iz,iw),adjCoe(iob,ij,iz),(wu(iw) - wl(iw)),sza(iob), &
!LFR> 			 int(sza(iob)/20.0+1)
!LFR> 	      END IF
!LFR-DEBUG
!IF(valJ(iob,ij,iz)<0.0) THEN
!      WRITE (88,FMT='(5(I4.4,1X),7(E18.10,1X))')     &
!                                   mynum,iob,ij,iz,iw,valJ(iob,ij,iz), &
!                                   saFlux(iob,iz,iw),sj(iob,ij,iz,iw), &
!				   adjCoe(iob,ij,iz),valja,djdw,(wu(iw) - wl(iw))
!				   CALL flush(88)
!      WRITE (88,FMT='("Saflux: ",9(E18.10,1X))') saFlux(iob,iz,iw), &
!     		                  etf(iw),wc(iw),dirsun,fdir(iob,iz),difdn, &
!				  fdn(iob,iz),difup,fup(iob,iz)
!
!
!END IF
!LFR-END-DEBUG

     		END DO
     	      END DO
    	  END IF
       END DO
    END DO
    !end wavelength loop

    !* reset wavelength scale if needed:
    IF(lrefr) THEN
       DO iob=1,nbl
     	  WRITE(*,*) 'applying vacuum to air wavelength shift' , airout
     	  mrefr = -mrefr
     	  CALL wshift(mrefr, nw, wl, airout)
     	  CALL wshift(mrefr, nwint, wc, airout)
     	  CALL wshift(mrefr, nwint, wu, airout)
       END DO
    END IF

901       FORMAT('zenith =   ',f10.1)
902       FORMAT(i10,a30)
903       FORMAT(i10,4E13.3)


  END SUBROUTINE tuv


  SUBROUTINE calcoe(nz,ij,c,xzin,adjin,adjcoe)
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  ADJCOE - REAL, coross section adjust coefficients (in and out)           =*
  !=  c(5,kj)- polynomal coef                                                   =*
  !=  tt     - nomarlized temperature
  !-----------------------------------------------------------------------------*
  !=  EDIT HISTORY:                                                            =*
  !=  11/2000 XUEXI                                                            =*
  !-----------------------------------------------------------------------------*
  != This program is free software;  you can redistribute it and/or modify     =*
  != it under the terms of the GNU General Public License as published by the  =*
  != Free Software Foundation;  either version 2 of the license, or (at your   =*
  != option) any later version.                                                =*
  != The TUV package is distributed in the hope that it will be useful, but    =*
  != WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
  != LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
  != License for more details.                                                 =*
  != To obtain a copy of the GNU General Public License, write to:             =*
  != Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
  !-----------------------------------------------------------------------------*
  != To contact the authors, please mail to:                                   =*
  != Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
  != send email to:  sasha@ucar.edu                                            =*
  !-----------------------------------------------------------------------------*
  != Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nz
  INTEGER, INTENT(IN) :: ij
  REAL, INTENT(IN)    :: c(5,kj)
  REAL, INTENT(IN)    :: xzin(kz)
  REAL, INTENT(IN)    :: adjin
  REAL, INTENT(OUT)   :: adjcoe(kj,kz)
  REAL :: x2,x3,x4

  INTEGER :: k
  REAL :: xz(kz)

  DO k=1,nz
    xz(k)=xzin(k)*adjO3(ij)
    x2=xz(k)*xz(k)
    x3=x2*xz(k)
    x4=x3*xz(k)
    adjcoe(ij,k) =adjin *(1.0 + ( c(1,ij)  &
        +  c(2,ij)*xz(k) +  c(3,ij)*x2  &
        +  c(4,ij)*x3 +  c(5,ij)*x4)*0.01)
!DEBUG BEG
!   IF(adjcoe(ij,k)<0) THEN
!      WRITE(45,FMT= '("calcoe: ",2(I3.3,1X),8(F15.10,1X))') &
!            ij,k,adjcoe(ij,k),adjin,xz(k),c(1,ij),c(2,ij),c(3,ij),c(4,ij),c(5,ij)
!   END IF
!!DEBUG END
  END DO

  END SUBROUTINE calcoe

  REAL FUNCTION fery(w)
    IMPLICIT NONE

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Calculate the action spectrum value for erythema at a given wavelength   =*
  !=  according to: McKinlay, A.F and B.L.Diffey, A reference action spectrum  =*
  !=  for ultraviolet induced erythema in human skin, CIE Journal, vol 6,      =*
  !=  pp 17-22, 1987.                                                          =*
  !=  Value at 300 nm = 0.6486                                                 =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  W - REAL, wavelength (nm)                                             (I)=*
  !-----------------------------------------------------------------------------*

  REAL, INTENT(IN)                     :: w

  IF (w < 250.) THEN
    fery = 1.
  ! outside the ery spectrum range
  ELSE IF ((w >= 250.) .AND. (w < 298)) THEN
    fery = 1.
  ELSE IF ((w >= 298.) .AND. (w < 328.)) THEN
    fery = 10.**( 0.094*(298.-w) )
  ELSE IF ((w >= 328.) .AND. (w < 400.)) THEN
    fery = 10.**( 0.015*(139.-w) )
  ELSE
    fery = 1.e-36
  ! outside the ery spectrum range
  END IF

  END FUNCTION fery

  !=============================================================================*

  REAL FUNCTION fo3qy(w,t)
  IMPLICIT NONE

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  ! function to calculate the quantum yield O3 + hv -> O(1D) + O2,             =*
  ! according to JPL 2000 recommendation:                                      =*
  !-----------------------------------------------------------------------------*


  REAL, INTENT(IN)  :: w
  REAL, INTENT(IN)  :: t
  REAL :: kt
  REAL,PARAMETER :: a(3)=(/0.887, 2.35, 57.0/)
  REAL,PARAMETER :: w0(3)=(/302.0, 311.1, 313.9/)
  REAL,PARAMETER :: nu(3)=(/0.0, 820.0, 1190.0/)
  REAL,PARAMETER :: om(3)=(/7.9, 2.2, 7.4/)

  fo3qy = 0.
  kt = 0.695 * t

  IF(w <= 300.) THEN
    fo3qy = 0.95
  ELSE IF(w > 300. .AND. w <= 330.) THEN
    fo3qy = 0.06 +  &
        a(1)                           *EXP(-((w-w0(1))/om(1))**4)+  &
        a(2)*(t/300.)**4*EXP(-nu(2)/kt)*EXP(-((w-w0(2))/om(2))**2)+  &
        a(3)            *EXP(-nu(3)/kt)*EXP(-((w-w0(3))/om(3))**2)
  ELSE IF(w > 330. .AND. w <= 345.) THEN
    fo3qy = 0.06
  ELSE IF(w > 345.) THEN
    fo3qy = 0.
  END IF

  END FUNCTION fo3qy

  REAL FUNCTION fo3qy2(w,t)
  IMPLICIT NONE

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  ! function to calculate the quantum yield O3 + hv -> O(1D) + O2,             =*
  ! according to:
  ! Matsumi, Y., F. J. Comes, G. Hancock, A. Hofzumanhays, A. J. Hynes,
  ! M. Kawasaki, and A. R. Ravishankara, QUantum yields for production of O(1D)
  ! in the ultraviolet photolysis of ozone:  Recommendation based on evaluation
  ! of laboratory data, J. Geophys. Res., 107, 10.1029/2001JD000510, 2002.
  !-----------------------------------------------------------------------------*


  REAL, INTENT(IN)  :: w
  REAL, INTENT(IN)  :: t
  REAL :: kt

  REAL,PARAMETER :: a(3)=(/ 0.8036, 8.9061, 0.1192/)
  REAL,PARAMETER :: x(3)=(/ 304.225, 314.957, 310.737/)
  REAL,PARAMETER :: om(3)=(/ 5.576, 6.601, 2.187/)
  REAL :: q1, q2

  fo3qy2 = 0.0
  kt = 0.695 * t
  q1 = 1.0
  q2 = EXP(-825.518/kt)

  IF(w <= 305.) THEN
    fo3qy2 = 0.90
  ELSE IF(w > 305. .AND. w <= 328.) THEN

    fo3qy2 = 0.0765 +  &
        a(1)*             (q1/(q1+q2))*EXP(-((x(1)-w)/om(1))**4)+  &
        a(2)*(t/300.)**2 *(q2/(q1+q2))*EXP(-((x(2)-w)/om(2))**2)+  &
        a(3)*(t/300.)**1.5            *EXP(-((x(3)-w)/om(3))**2)

  ELSE IF(w > 328. .AND. w <= 340.) THEN
    fo3qy2 = 0.08
  ELSE IF(w > 340.) THEN
    fo3qy2 = 0.
  END IF

  END FUNCTION fo3qy2

  !=============================================================================*

  REAL FUNCTION fsum(n,x)
  IMPLICIT NONE

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Compute the sum of the first N elements of a floating point vector.      =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  N  - INTEGER, number of elements to sum                               (I)=*
  !=  X  - REAL, vector whose components are to be summed                   (I)=*
  !-----------------------------------------------------------------------------*
  INTEGER, INTENT(IN) :: n
  REAL, INTENT(IN)    :: x(n)

  ! local:
  INTEGER :: i

  fsum = 0.
  DO   i = 1, n
    fsum=fsum+x(i)
  END DO

  END FUNCTION fsum

  !=============================================================================*

  REAL FUNCTION futr(w)
  IMPLICIT NONE

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Calculate the action spectrum value for skin cancer of albino hairless   =*
  !=  mice at a given wavelength according to:  deGRuijl, F.R., H.J.C.M.Steren-=*
  !=  borg, P.D.Forbes, R.E.Davies, C.Colse, G.Kelfkens, H.vanWeelden,         =*
  !=  and J.C.van der Leun, Wavelength dependence of skin cancer induction by  =*
  !=  ultraviolet irradiation of albino hairless mice, Cancer Research, vol 53,=*
  !=  pp. 53-60, 1993                                                          =*
  !=  (Action spectrum for carcinomas)                                         =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  W  - REAL, wavelength (nm)                                            (I)=*
  !-----------------------------------------------------------------------------*


  REAL, INTENT(IN)                     :: w

  ! local:
  REAL :: a1, a2, a3, a4, a5, x1, x2, x3, x4, x5,  &
      t1, t2, t3, t4, t5, b1, b2, b3, b4, b5,  &
      p

  a1 = -10.91
  a2 = - 0.86
  a3 = - 8.60
  a4 = - 9.36
  a5 = -13.15

  x1 = 270.
  x2 = 302.
  x3 = 334.
  x4 = 367.
  x5 = 400.

  t1 = (w-x2)*(w-x3)*(w-x4)*(w-x5)
  t2 = (w-x1)*(w-x3)*(w-x4)*(w-x5)
  t3 = (w-x1)*(w-x2)*(w-x4)*(w-x5)
  t4 = (w-x1)*(w-x2)*(w-x3)*(w-x5)
  t5 = (w-x1)*(w-x2)*(w-x3)*(w-x4)

  b1 = (x1-x2)*(x1-x3)*(x1-x4)*(x1-x5)
  b2 = (x2-x1)*(x2-x3)*(x2-x4)*(x2-x5)

  b3 = (x3-x1)*(x3-x2)*(x3-x4)*(x3-x5)
  b4 = (x4-x1)*(x4-x2)*(x4-x3)*(x4-x5)
  b5 = (x5-x1)*(x5-x2)*(x5-x3)*(x5-x4)

  p = a1*t1/b1 + a2*t2/b2 + a3*t3/b3 + a4*t4/b4 + a5*t5/b5

  futr  = EXP(p)

  END FUNCTION futr


  SUBROUTINE gridw(wstart, wstop, nwint)

  IMPLICIT NONE

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Create the wavelength grid for all interpolations and radiative transfer =*
  !=  calculations.  Grid may be irregularly spaced.  Wavelengths are in nm.   =*
  !=  No gaps are allowed within the wavelength grid.                          =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW  - INTEGER, number of wavelength grid _points_                     (O)=*
  !=  WL  - REAL, vector carrying the lower limit of each wavel. interval   (O)=*
  !=  WC  - REAL, vector carrying the center wavel of each wavel. interval  (O)=*
  !=              (wc(i) = 0.5*(wl(i)+wu(i), i = 1..NW-1)                      =*
  !=  WU  - REAL, vector carrying the upper limit of each wavel. interval   (O)=*
  !=
  !=  MOPT- INTEGER OPTION for wave-length IF 3 good for JO2                (O)=*
  !-----------------------------------------------------------------------------*


  REAL, INTENT(IN)     :: wstart
  REAL, INTENT(IN)     :: wstop
  INTEGER, INTENT(IN)  :: nwint

  INTEGER :: mopt
  REAL :: wincr
  INTEGER :: iw

  CHARACTER (LEN=200) :: fi
  CHARACTER (LEN=20) :: wlabel

  REAL :: dum

  LOGICAL :: ok
  !_______________________________________________________________________

  !*** chose wavelength grid

  ! some pre-set options
  !     mopt = 1    equal spacing
  !     mopt = 2    grid defined in data table
  !     mopt = 3    user-defined
  !     mopt = 4    fast-TUV, troposheric wavelengths only
  mopt = 1
  IF(nwint == -156) mopt = 2
  IF(nwint <= -1 .AND. nwint >= -20)  mopt = 4     ! fast-J XUEXI

  SELECT CASE (mopt)
     CASE(1)
     	wlabel = 'equal spacing'
     	nw = nwint + 1
     	wincr = (wstop - wstart) / FLOAT (nwint)
     	DO   iw = 1, nw-1
     	  wl(iw) = wstart + wincr*FLOAT(iw-1)
     	  wu(iw) = wl(iw) + wincr
     	  wc(iw) = ( wl(iw) + wu(iw) )/2.
     	END DO
     	wl(nw) = wu(nw-1)
     CASE (2)
        ! Input from table.  In this example:
        ! Wavelength grid will be read from a file.
        ! First line of table is:  nw = number of wavelengths (no. of intervals + 1)
        ! Then, nw wavelengths are read in, and assigned to wl(iw)
        ! Finally, wu(iw) and wc(iw) are computed from wl(iw)

        !      wlabel = 'isaksen.grid'
        !wlabel = 'combined.grid'

        fi = trim(files(136)%fileName)
        OPEN(UNIT=kin,FILE=fi,STATUS='old')
        READ(kin,*) nw
        DO iw = 1, nw
           READ(kin,*) wl(iw)
        END DO
        CLOSE(kin)
        DO iw = 1, nw-1
           wu(iw) = wl(iw+1)
           wc(iw) = 0.5*(wl(iw) + wu(iw))
        END DO
     CASE (3)
        ! user-defined grid.  In this example, a single calculation is used to
        ! obtain results for two 1 nm wide intervals centered at 310 and 400 nm:
        ! interval 1 : 1 nm wide, centered at 310 nm
        ! interval 3 : 2 nm wide, centered at 400 nm
        ! (inteval 2 : 310.5 - 399.5 nm, required to connect intervals 1 & 3)

        nw = 4
        wl(1) = 309.5
        wl(2) = 310.5
        wl(3) = 399.5
        wl(4) = 400.5
        DO iw = 1, nw-1
           wu(iw) = wl(iw+1)
           wc(iw) = 0.5*(wl(iw) + wu(iw))
        END DO
     CASE (4)
        wlabel = 'fast-TUV tropospheric grid'

        fi = trim(files(137)%fileName)
        OPEN(UNIT=kin,FILE=fi,STATUS='old')
        DO iw = 1, 4
          READ(kin,*)
        END DO
        !PRINT *, 'LFR->',fi
        ! skip wavelength shorter than 205 nm
        DO iw = 1, 6
           READ(kin,*)
        END DO
        nw = ABS(nwint) + 1
        DO iw = 1, nw-1
           READ(kin,*) dum, wl(iw),dum,dum
        END DO
        wl(nw) = dum
        DO iw = 1, nw-1
           wu(iw) = wl(iw+1)
           wc(iw) = 0.5*(wl(iw) + wu(iw))
        END DO
        CLOSE(kin)
  END SELECT
  CALL gridck(kw,nw,wl,ok)

  IF (.NOT. ok) THEN
    WRITE(*,*)'STOP in GRIDW:  The w-grid does not make sense'
    STOP
  END IF

  !_______________________________________________________________________


  END SUBROUTINE gridw

SUBROUTINE gridck(k,n,x,ok)
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Check a grid X for various improperties.  The values in X have to comply =*
!=  with the following rules:                                                =*
!=  1) Number of actual points cannot exceed declared length of X            =*
!=  2) Number of actual points has to be greater than or equal to 2          =*
!=  3) X-values must be non-negative                                         =*
!=  4) X-values must be unique                                               =*
!=  5) X-values must be in ascending order                                   =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  K  - INTEGER, length of X as declared in the calling program          (I)=*
!=  N  - INTEGER, number of actual points in X                            (I)=*
!=  X  - REAL, vector (grid) to be checked                                (I)=*
!=  OK - LOGICAL, .TRUE. -> X agrees with rules 1)-5)                     (O)=*
!=                .FALSE.-> X violates at least one of 1)-5)                 =*
!-----------------------------------------------------------------------------*
IMPLICIT NONE

INTEGER, INTENT(IN)    :: k
INTEGER, INTENT(IN)    :: n
REAL, INTENT(IN)       :: x(k)
LOGICAL, INTENT(OUT)   :: ok
! local:
INTEGER :: i
INTEGER,PARAMETER :: kout=6
!_______________________________________________________________________

ok = .true.

! check if dimension meaningful and within bounds

IF (n > k) THEN
  ok = .false.
  WRITE(kout,100)
  RETURN
END IF
100 FORMAT('Number of data exceeds dimension')

IF (n < 2) THEN
  ok = .false.
  WRITE(kout,101)
  RETURN
END IF
101 FORMAT('Too few data, number of data points must be >= 2')

! disallow negative grid values

IF(x(1) < 0.) THEN
  ok = .false.
  WRITE(kout,105)
  RETURN
END IF
105 FORMAT('Grid cannot start below zero')

! check sorting

DO   i = 2, n
  IF( x(i) <= x(i-1)) THEN
    ok = .false.
    WRITE(kout,110)
    RETURN
  END IF
END DO
110 FORMAT('Grid is not sorted or contains multiple values')
!_______________________________________________________________________

END SUBROUTINE gridck

  SUBROUTINE la_srb(nz,nw,wl,o2xs1,nbl,tLev,dtO2,o2xs,scol,vcol)
    IMPLICIT NONE

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Compute equivalent optical depths for O2 absorption, and O2 effective    =*
  !=  absorption cross sections, parameterized in the Lyman-alpha and SR bands =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
  !=            grid                                                           =*
  !=  Z       - REAL, specified altitude working grid (km)                  (I)=*
  !=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
  !=            wavelength grid                                                =*
  !=  WL      - REAL, vector of lxower limits of wavelength intervals in    (I)=*
  !=            working wavelength grid                                        =*
  !=  CZ      - REAL, number of air molecules per cm^2 at each specified    (I)=*
  !=            altitude layer                                                 =*
  !=  ZEN     - REAL, solar zenith angle                                    (I)=*
  !=                                                                           =*
  !=  O2XS1   - REAL, O2 cross section from rdo2xs                          (I)=*
  !=                                                                           =*
  !=  DTO2    - REAL, optical depth due to O2 absorption at each specified  (O)=*
  !=            vertical layer at each specified wavelength                    =*
  !=  O2XS    - REAL, molecular absorption cross section in SR bands at     (O)=*
  !=            each specified altitude and wavelength.  Includes Herzberg     =*
  !=            continuum.                                                     =*
  !-----------------------------------------------------------------------------*
  INTEGER, INTENT(IN) :: nbl
  INTEGER, INTENT(IN) :: nz(nbl)
  INTEGER, INTENT(IN) :: nw
  REAL, INTENT(IN)    :: wl(kw)
  REAL, INTENT(IN)    :: o2xs1(kw)
  REAL, INTENT(IN)    :: tLev(nbl,kz)
  REAL, INTENT(IN)    :: scol(nbl,kz)
  REAL, INTENT(IN)    :: vcol(nbl,kz)
  REAL, INTENT(INOUT) :: dtO2(nbl,kz,kw)
  REAL, INTENT(OUT)   :: o2xs(nbl,kz,kw)

  INTEGER :: iz, iw, iob
  REAL,ALLOCATABLE,DIMENSION(:,:) :: o2col!(kz)
  REAL,ALLOCATABLE,DIMENSION(:,:) :: secchi!(kz)

  ! Lyman-alpha variables
  ! O2 optical depth and equivalent cross section in the Lyman-alpha region

  INTEGER, PARAMETER :: nla=1
  INTEGER, PARAMETER :: kla = 2
  REAL,DIMENSION(kla),PARAMETER :: wlla=(/ 121.4, 121.9/) ! Wavelengths for Lyman alpha and SRB parameterizations
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: dto2la!(kz, kla-1)
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: o2xsla!(kz, kla-1)

  ! grid on which Koppers' parameterization is defined
  ! O2 optical depth and equivalent cross section on Koppers' grid

  INTEGER,PARAMETER :: nsrb=17
  INTEGER, PARAMETER :: ksrb = 18
  REAL,DIMENSION(ksrb),PARAMETER :: wlsrb=(/ &
        174.4, 177.0, 178.6, 180.2, 181.8, 183.5, 185.2, 186.9,  &
        188.7, 190.5, 192.3, 194.2, 196.1, 198.0, 200.0, 202.0,  &
        204.1, 205.8/)
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: dto2k!(kz, ksrb-1)
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: o2xsk!(kz, ksrb-1)

  INTEGER :: i

  !----------------------------------------------------------------------
  ! initalize O2 cross sections
  !----------------------------------------------------------------------
  DO iz = 1, zLimit
    DO iw =1, nw - 1
      DO iob=1,nbl
         IF(iz>nz(iob)) CYCLE
         o2xs(iob,iz,iw) = o2xs1(iw)
      END DO
    END DO
  END DO

  IF(wl(1) > wlsrb(nsrb)) RETURN

  !----------------------------------------------------------------------
  ! On first call, check that the user wavelength grid, WL(IW), is compatible
  ! with the wavelengths for the parameterizations of the Lyman-alpha and SRB.
  ! Also compute and save corresponding grid indices (ILA, ISRB)
  !----------------------------------------------------------------------
  IF (call1) THEN
    !* locate Lyman-alpha wavelengths on grid
    ila = 0
    DO iw = 1, nw
      IF(ABS(wl(iw) - wlla(1)) < 10.*precis) THEN
        ila = iw
        EXIT
      END IF
    END DO
    ! check
    IF(ila == 0) STOP ' Lyman alpha grid mis-match - 1'
    DO i = 2, nla + 1
      IF(ABS(wl(ila + i - 1) - wlla(i)) > 10.*precis) THEN
        WRITE(*,*) 'Lyman alpha grid mis-match - 2'
        STOP
      END IF
    END DO
    !* locate Schumann-Runge wavelengths on grid
    isrb = 0
    DO iw = 1, nw
      IF(ABS(wl(iw) - wlsrb(1)) < 10.*precis) THEN
        isrb = iw
        EXIT
      END IF
    END DO
    ! check
    IF(isrb == 0) STOP ' SRB grid mis-match - 1'
    DO i = 2, nsrb + 1
      IF(ABS(wl(isrb + i - 1) - wlsrb(i)) > 10.* precis) THEN
        WRITE(*,*) ' SRB grid mismatch - w'
        STOP
      END IF
    END DO
    call1 = .false.
  END IF

  !Local allocations for vetorization
  ALLOCATE(o2col(nbl,iz))
  ALLOCATE(secchi(nbl,iz))
  ALLOCATE(dto2la(nbl,kz,kla-1))
  ALLOCATE(o2xsla(nbl,kz,kla-1))
  ALLOCATE(dto2k(nbl,kz, ksrb-1))
  ALLOCATE(o2xsk(nbl,kz, ksrb-1))
  !----------------------------------------------------------------------
  ! Slant O2 column and x-sections.
  !----------------------------------------------------------------------
  DO iob=1,nbl
     DO iz = 1, nz(iob)
        o2col(iob,iz) = 0.2095 * scol(iob,iz)
     END DO
  END DO

  !----------------------------------------------------------------------
  ! Effective secant of solar zenith angle.
  ! Use 2.0 if no direct sun (value for isotropic radiation)
  ! For nz, use value at nz-1
  !----------------------------------------------------------------------
  DO iob=1,nbl
     DO i = 1,nz(iob) - 1
        secchi(iob,i) = scol(iob,i)/ vcol(iob,i)
        IF(scol(iob,i) > largest*0.1) secchi(iob,i) = 2.
     END DO
  END DO
  DO iob=1,nbl
     secchi(iob,nz(iob)) = &
               secchi(iob,nz(iob)-1)
  END DO
  !---------------------------------------------------------------------
  ! Lyman-Alpha parameterization, output values of O2 optical depth
  ! and O2 effective (equivalent) cross section
  !----------------------------------------------------------------------
  CALL lymana(nz,o2col,secchi,dto2la,o2xsla,nbl)
  DO iw = ila, ila + nla - 1
    DO iob=1,nbl
       DO iz = 1, nz(iob)
           dtO2(iob,iz,iw) = dto2la(iob,iz, iw - ila + 1)
           o2xs(iob,iz,iw) = o2xsla(iob,iz, iw - ila + 1)
       END DO
    END DO
  END DO

  !------------------------------------------------------------------------------
  ! Koppers' parameterization of the SR bands, output values of O2
  ! optical depth and O2 equivalent cross section
  !------------------------------------------------------------------------------

  CALL schum(nz,o2col,secchi,dto2k,o2xsk,nbl,tLev)
  DO iw = isrb, isrb + nsrb - 1
    DO iob=1,nbl
       DO iz = 1, nz(iob)
          dto2(iob,iz,iw) = dto2k(iob,iz, iw - isrb + 1)
          o2xs(iob,iz,iw) = o2xsk(iob,iz, iw - isrb + 1)
       END DO
    END DO
  END DO

  END SUBROUTINE la_srb

  !=============================================================================*

  SUBROUTINE lymana(nz,o2col,secchi,dto2la,o2xsla,nbl)
  IMPLICIT NONE

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Calculate the effective absorption cross section of O2 in the Lyman-Alpha=*
  !=  bands and an effective O2 optical depth at all altitudes.  Parameterized =*
  !=  after:  Chabrillat, S., and G. Kockarts, Simple parameterization of the  =*
  !=  absorption of the solar Lyman-Alpha line, Geophysical Research Letters,  =*
  !=  Vol.24, No.21, pp 2659-2662, 1997.                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
  !=            grid                                                           =*
  !=  O2COL   - REAL, slant overhead O2 column (molec/cc) at each specified (I)=*
  !=            altitude                                                       =*
  !=  DTO2LA  - REAL, optical depth due to O2 absorption at each specified  (O)=*
  !=            vertical layer                                                 =*
  !=  O2XSLA  - REAL, molecular absorption cross section in LA bands        (O)=*
  !-----------------------------------------------------------------------------*
  INTEGER, INTENT(IN) :: nbl
  INTEGER, INTENT(IN) :: nz(nbl)
  REAL, INTENT(IN)    :: o2col(nbl,kz)
  REAL, INTENT(IN)    :: secchi(nbl,kz)
  REAL, INTENT(OUT)   :: dto2la(nbl,kz,*)
  REAL, INTENT(OUT)   :: o2xsla(nbl,kz,*)

   DOUBLE PRECISION :: rm(nbl,kz), ro2(nbl,kz)

   DOUBLE PRECISION, DIMENSION(3), PARAMETER :: bbb=(/6.8431D-01,  2.29841D-01,  8.65412D-02/)
   DOUBLE PRECISION, DIMENSION(3), PARAMETER :: ccc=(/8.22114D-21, 1.77556D-20,  8.22112D-21/)
   DOUBLE PRECISION, DIMENSION(3), PARAMETER :: ddd=(/6.0073D-21,  4.28569D-21,  1.28059D-20/)
   DOUBLE PRECISION, DIMENSION(3), PARAMETER :: eee=(/8.21666D-21, 1.63296D-20,  4.85121D-17/)

  INTEGER :: iz, i, iob
  REAL :: xsmin

  !------------------------------------------------------------------------------*
  !sm:  set minimum cross section
  xsmin = 1.e-20

  DO iob=1,nbl
     DO iz = 1, nz(iob)
        rm(nbl,iz) = 0.0d+00
        ro2(nbl,iz) = 0.0d+00
     END DO
  END DO
  ! calculate reduction factors at every altitude
  DO iob=1,nbl
     DO iz = 1, nz(iob)
       DO i = 1, 3
          rm(nbl,iz) = rm(nbl,iz) + bbb(i) * DEXP(-ccc(i) * &
	               DBLE(o2col(nbl,iz)))
          ro2(nbl,iz) = ro2(nbl,iz) + ddd(i) * DEXP(-eee(i) * &
	               DBLE(o2col(nbl,iz)))
       END DO
    END DO
  END DO

  ! calculate effective O2 optical depths and effective O2 cross sections
  DO iob=1,nbl
     DO iz = 1, nz(iob)-1
        IF (rm(nbl,iz) > 1.0D-100) THEN
           IF (ro2(nbl,iz) > 1.0d-100) THEN
               o2xsla(nbl,iz,1) = ro2(nbl,iz)/rm(nbl,iz)
           ELSE
               o2xsla(nbl,iz,1) = xsmin
           END IF
           IF (rm(nbl,iz+1) > 0.) THEN
              dto2la(nbl,iz,1) = &
	      LOG(rm(nbl,iz+1))/secchi(nbl,iz+1)- &
              LOG(rm(nbl,iz))/secchi(nbl,iz)
           ELSE
              dto2la(nbl,iz,1) = 1000.
           END IF
        ELSE
           dto2la(nbl,iz,1) = 1000.
           o2xsla(nbl,iz,1) = xsmin
        END IF
     END DO
  END DO

  DO iob=1,nbl
    ! do top layer separately
    dto2la(nbl,nz(iob),1) = 0.
  END DO

  DO iob=1,nbl
     IF(rm(nbl,nz(iob)) > 1.0d-100) THEN
        o2xsla(nbl,nz(iob),1) = &
	                   ro2(nbl,nz(iob))/ &
	                   rm(nbl,nz(iob))
     ELSE
        o2xsla(nbl,nz(iob),1) = xsmin
     END IF
  END DO

  END SUBROUTINE lymana

  SUBROUTINE schum(nz,o2col,secchi,dto2,o2xsk,nobl,tLev)
  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Calculate the equivalent absorption cross section of O2 in the SR bands. =*
  !=  The algorithm is based on parameterization of G.A. Koppers, and          =*
  !=  D.P. Murtagh [ref. Ann.Geophys., 14 68-79, 1996]                         =*
  !=  Final values do include effects from the Herzberg continuum.             =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
  !=            grid                                                           =*
  !=  O2COL   - REAL, slant overhead O2 column (molec/cc) at each specified (I)=*
  !=            altitude                                                       =*
  !=  TLEV    - tmeperature at each level                                   (I)=*
  !=  SECCHI  - ratio of slant to vertical o2 columns                       (I)=*
  !=  DTO2    - REAL, optical depth due to O2 absorption at each specified  (O)=*
  !=            vertical layer at each specified wavelength                    =*
  !=  O2XSK  - REAL, molecular absorption cross section in SR bands at     (O)=*
  !=            each specified wavelength.  Includes Herzberg continuum        =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nobl
  INTEGER, INTENT(IN) :: nz(nobl)
  REAL, INTENT(IN)    :: o2col(nobl,kz)
  REAL, INTENT(IN)    :: secchi(nobl,kz)
  REAL, INTENT(IN)    :: tLev(nobl,kz)

  REAL, INTENT(OUT)   :: dto2(nobl,kz,17)
  REAL, INTENT(OUT)   :: o2xsk(nobl,kz,17)

  REAL :: o2col1(nobl,kz)
  INTEGER :: i, k , iob
  INTEGER :: ktop(nobl), ktop1(nobl), kbot(nobl)
  REAL :: xs(17), x

  !------------------------------------------
  !sm  Initialize cross sections to values
  !sm  at large optical depth
  !------------------------------------------
  DO i = 1, 17
     DO iob=1,nobl
       DO k = 1, nz(iob)
          o2xsk(iob,k,i) = xslod(i)
       END DO
    END DO
  END DO

  !------------------------------------------
  !     Calculate cross sections
  !sm:  Set smallest O2col = exp(38.) molec cm-2
  !sm     to stay in range of parameterization
  !sm     given by Koppers et al. at top of atm.
  !------------------------------------------
  ktop = 121
  kbot = 0

  DO iob=1,nobl
     DO k=1,nz(iob)    !! loop for alt
        o2col1(iob,k) = MAX(o2col(iob,k),EXP(38.))
        x  = ALOG(o2col1(iob,k))
       IF (x < 38.0) THEN
          ktop1(iob) = k-1
          ktop(iob)  = MIN(ktop1(iob),ktop(iob))
       ELSE IF (x > 56.0) THEN
         kbot(iob) = k
       ELSE
          CALL effxs( x, tLev(iob,k), xs )
          DO i=1,17
            o2xsk(iob,k,i) = xs(i)
          END DO
       END IF
     END DO
  END DO                    !! finish loop for alt

  !------------------------------------------
  !  fill in cross section where X is out of range
  !  by repeating edge table values
  !------------------------------------------
  !sm do not allow kbot = nz to avoid division by zero in
  !   no light case.
  DO iob=1,nobl
    IF(kbot(iob) == nz(iob)) kbot(iob) = nz(iob) - 1
  END DO
  DO iob=1,nobl
     DO k=1,kbot(iob)
        DO i=1,17
           o2xsk(iob,k,i) = o2xsk(iob,kbot(iob)+1,i)
        END DO
     END DO
  END DO
  DO iob=1,nobl
     DO k=ktop(iob)+1,nz(iob)
        DO i=1,17
           o2xsk(iob,k,i) = o2xsk(iob,ktop(iob),i)
        END DO
     END DO
  END DO
  !------------------------------------------
  !  Calculate incremental optical depths
  !------------------------------------------
  DO i=1,17                   ! loop over wavelength
     DO iob=1,nobl
        DO k=1,nz(iob)-1            ! loop for alt
          !... calculate an optical depth weighted by density
          !sm:  put in mean value estimate, if in shade
          IF (ABS(1. - o2col1(iob,k+1)/o2col1(iob,k)) <= 2.*precis) THEN
             dto2(iob,k,i) = o2xsk(iob,k+1,i)*o2col1(iob,k+1)/(nz(iob)-1)
          ELSE
             dto2(iob,k,i) = ABS((o2xsk(iob,k+1,i)*o2col1(iob,k+1) - &
	                  o2xsk(iob,k,i)*o2col1(iob,k))/  &
                         (1.0 + ALOG(o2xsk(iob,k+1,i)/o2xsk(iob,k,i)) / &
			 ALOG(o2col1(iob,k+1)/o2col1(iob,k))))
             !... change to vertical optical depth
             dto2(iob,k,i) = 2. * dto2(iob,k,i)/(secchi(iob,k)+secchi(iob,k+1))
          END IF
       END DO
    END DO
  END DO

  DO i=1,17                   ! loop over wavelength
     DO iob=1,nobl
       dto2(iob,nz(iob),i) = 0.0       ! set optical depth to zero at top
     END DO
  END DO

  END SUBROUTINE schum

  !=============================================================================*
  SUBROUTINE effxs( x, t, xs )


  !     Subroutine for evaluating the effective cross section
  !     of O2 in the Schumann-Runge bands using parameterization
  !     of G.A. Koppers, and D.P. Murtagh [ref. Ann.Geophys., 14
  !     68-79, 1996]

  !     method:
  !     ln(xs) = A(X)[T-220]+B(X)
  !     X = log of slant column of O2
  !     A,B calculated from Chebyshev polynomial coeffs
  !     AC and BC using NR routine chebev.  Assume interval
  !     is 38<ln(NO2)<56.

  !     Revision History:

  !     drm 2/97  initial coding

  !-------------------------------------------------------------
  IMPLICIT NONE

  REAL*4, INTENT(IN OUT) :: x
  REAL*4, INTENT(IN)     :: t
  REAL*4, INTENT(OUT)    :: xs(17)

  REAL*4 a(17), b(17)
  INTEGER :: i

  CALL calc_params( x, a, b )

  DO i = 1,17
    xs(i) = EXP( a(i)*( t - 220.) + b(i) )
  END DO

  END SUBROUTINE effxs

  !=============================================================================*
  SUBROUTINE calc_params( x, a, b )

  !-------------------------------------------------------------

  !       calculates coefficients (A,B), used in calculating the
  ! effective cross section, for 17 wavelength intervals
  !       as a function of log O2 column density (X)
  !       Wavelength intervals are defined in WMO1985

  !-------------------------------------------------------------
  IMPLICIT NONE

  REAL*4, INTENT(IN OUT) :: x
  REAL*4, INTENT(OUT)    :: a(17)
  REAL*4, INTENT(OUT)    :: b(17)

  INTEGER :: i

  !       call Chebyshev Evaluation routine to calc A and B from
  ! set of 20 coeficients for each wavelength
  DO i=1,17
    a(i) = chebev(38.0 , 56.0, ac(1,i), 20, x)
    b(i) = chebev(38.0 , 56.0, bc(1,i), 20, x)
  END DO

  END SUBROUTINE calc_params

  !=============================================================================*

  SUBROUTINE init_xs()
  IMPLICIT NONE

  !       locals
  INTEGER*4 in_lun ! file unit number
  INTEGER*4 i, j

  in_lun = 11

  OPEN (UNIT=in_lun, FILE=trim(files(138)%fileName),FORM='FORMATTED')

  READ( in_lun, 901 )
  DO i = 1,20
    READ( in_lun, 903 ) ( ac(i,j), j=1,17 )
  END DO
  READ( in_lun, 901 )
  DO i = 1,20
    READ( in_lun, 903 ) ( bc(i,j), j=1,17 )
  END DO

  901    FORMAT( / )
  903    FORMAT( 17(e23.14,1X))

  998 CLOSE (in_lun)

  DO i=1,17
    wave_num(18-i) = 48250. + (500.*i)
  END DO

  END SUBROUTINE init_xs

  !=============================================================================*
  REAL*4  FUNCTION chebev(a,b,c,m,x)
  !     Chebyshev evaluation algorithm
  !     See Numerical recipes p193
  !-------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: m
  REAL*4, INTENT(IN)  :: a
  REAL*4, INTENT(IN)  :: b
  REAL*8, INTENT(IN)  :: c(m)
  REAL*4, INTENT(IN)  :: x


  INTEGER :: j
  REAL :: d,dd,sv,y,y2

  IF ((x-a)*(x-b) > 0.) THEN
    WRITE(6,*) 'X NOT IN RANGE IN CHEBEV', x
    chebev = 0.0
    RETURN
  END IF

  d=0.
  dd=0.
  y=(2.*x-a-b)/(b-a)
  y2=2.*y
  DO  j=m,2,-1
    sv=d
    d=y2*d-dd+c(j)
    dd=sv
  END DO
  chebev=y*d-dd+0.5*c(1)

  END FUNCTION chebev

  SUBROUTINE inter1(ng,xg,yg, n,x,y)
  IMPLICIT NONE

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Map input data given on single, discrete points, onto a discrete target  =*
  !=  grid.                                                                    =*
  !=  The original input data are given on single, discrete points of an       =*
  !=  arbitrary grid and are being linearly interpolated onto a specified      =*
  !=  discrete target grid.  A typical example would be the re-gridding of a   =*
  !=  given data set for the vertical temperature profile to match the speci-  =*
  !=  fied altitude grid.                                                      =*
  !=  Some caution should be used near the end points of the grids.  If the    =*
  !=  input data set does not span the range of the target grid, the remaining =*
  !=  points will be set to zero, as extrapolation is not permitted.           =*
  !=  If the input data does not encompass the target grid, use ADDPNT to      =*
  !=  expand the input array.                                                  =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NG  - INTEGER, number of points in the target grid                    (I)=*
  !=  XG  - REAL, target grid (e.g. altitude grid)                          (I)=*
  !=  YG  - REAL, y-data re-gridded onto XG                                 (O)=*
  !=  N   - INTEGER, number of points in the input data set                 (I)=*
  !=  X   - REAL, grid on which input data are defined                      (I)=*
  !=  Y   - REAL, input y-data                                              (I)=*
  !-----------------------------------------------------------------------------*

  INTEGER, INTENT(IN)                      :: ng
  INTEGER, INTENT(IN)                      :: n
  REAL, INTENT(IN)                         :: xg(ng)
  REAL, INTENT(IN)                         :: x(n)
  REAL, INTENT(IN)                         :: y(n)
  REAL, INTENT(OUT)                        :: yg(ng)

  ! local:
  REAL :: slope
  INTEGER :: jsave, i, j
  !_______________________________________________________________________

  jsave = 1
  DO   i = 1, ng
    yg(i) = 0.
    j = jsave
    10    CONTINUE
    IF ((x(j) > xg(i)) .OR. (xg(i) >= x(j+1))) THEN
      j = j+1
      IF (j <= n-1) GO TO 10
  !        ---- end of loop 10 ----
    ELSE
      slope = (y(j+1)-y(j)) / (x(j+1)-x(j))
      yg(i) = y(j) + slope * (xg(i) - x(j))
      jsave = j
    END IF
  END DO
  !_______________________________________________________________________

  END SUBROUTINE inter1

  !=============================================================================*

  SUBROUTINE inter2(ng,xg,yg,n,x,y,ierr)
  IMPLICIT NONE

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Map input data given on single, discrete points onto a set of target     =*
  !=  bins.                                                                    =*
  !=  The original input data are given on single, discrete points of an       =*
  !=  arbitrary grid and are being linearly interpolated onto a specified set  =*
  !=  of target bins.  In general, this is the case for most of the weighting  =*
  !=  functions (action spectra, molecular cross section, and quantum yield    =*
  !=  data), which have to be matched onto the specified wavelength intervals. =*
  !=  The average value in each target bin is found by averaging the trapezoi- =*
  !=  dal area underneath the input data curve (constructed by linearly connec-=*
  !=  ting the discrete input values).                                         =*
  !=  Some caution should be used near the endpoints of the grids.  If the     =*
  !=  input data set does not span the range of the target grid, an error      =*
  !=  message is printed and the execution is stopped, as extrapolation of the =*
  !=  data is not permitted.                                                   =*
  !=  If the input data does not encompass the target grid, use ADDPNT to      =*
  !=  expand the input array.                                                  =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NG  - INTEGER, number of bins + 1 in the target grid                  (I)=*
  !=  XG  - REAL, target grid (e.g., wavelength grid);  bin i is defined    (I)=*
  !=        as [XG(i),XG(i+1)] (i = 1..NG-1)                                   =*
  !=  YG  - REAL, y-data re-gridded onto XG, YG(i) specifies the value for  (O)=*
  !=        bin i (i = 1..NG-1)                                                =*
  !=  N   - INTEGER, number of points in input grid                         (I)=*
  !=  X   - REAL, grid on which input data are defined                      (I)=*
  !=  Y   - REAL, input y-data                                              (I)=*
  !-----------------------------------------------------------------------------*

  INTEGER, INTENT(IN)                      :: ng
  INTEGER, INTENT(IN)                      :: n
  REAL, INTENT(IN)                         :: xg(ng)
  REAL, INTENT(OUT)                        :: yg(ng)
  REAL, INTENT(IN)                         :: x(n)
  REAL, INTENT(IN)                         :: y(n)
  INTEGER, INTENT(OUT)                     :: ierr

  ! local:
  REAL :: area, xgl, xgu
  REAL :: darea, slope
  REAL :: a1, a2, b1, b2
  INTEGER :: ngintv
  INTEGER :: i, k, jstart

  !_______________________________________________________________________

  ierr = 0

  !  test for correct ordering of data, by increasing value of x
  DO   i = 2, n
    IF (x(i) <= x(i-1)) THEN
      ierr = 1
      WRITE(*,*)'data not sorted'
      RETURN
    END IF
  END DO

  DO i = 2, ng
    IF (xg(i) <= xg(i-1)) THEN
      ierr = 2
      WRITE(0,*) '>>> ERROR (inter2) <<<  xg-grid not sorted!'
      RETURN
    END IF
  END DO

  ! check for xg-values outside the x-range
  IF ( (x(1) > xg(1)) .OR. (x(n) < xg(ng)) ) THEN
    WRITE(0,*) '>>> ERROR (inter2) <<<  Data do not span '// 'grid.  '
    WRITE(0,*) '                        Use ADDPNT to '//  &
        'expand data and re-run.'
    STOP
  END IF

  !  find the integral of each grid interval and use this to
  !  calculate the average y value for the interval
  !  xgl and xgu are the lower and upper limits of the grid interval
  jstart = 1
  ngintv = ng - 1
  DO   i = 1,ngintv

  ! initalize:
    area = 0.0
    xgl = xg(i)
    xgu = xg(i+1)

  !  discard data before the first grid interval and after the
  !  last grid interval
  !  for internal grid intervals, start calculating area by interpolating
  !  between the last point which lies in the previous interval and the
  !  first point inside the current interval

    k = jstart
    IF (k <= n-1) THEN

  !  if both points are before the first grid, go to the next point
      30         CONTINUE
      IF (x(k+1) <= xgl) THEN
        jstart = k - 1
        k = k+1
        IF (k <= n-1) GO TO 30
      END IF


  !  if the last point is beyond the end of the grid, complete and go to the next
  !  grid
      40         CONTINUE
      IF ((k <= n-1) .AND. (x(k) < xgu)) THEN

        jstart = k-1

  ! compute x-coordinates of increment

        a1 = MAX(x(k),xgl)
        a2 = MIN(x(k+1),xgu)

  !  if points coincide, contribution is zero

        IF (x(k+1) == x(k)) THEN
          darea = 0.e0
        ELSE
          slope = (y(k+1) - y(k))/(x(k+1) - x(k))
          b1 = y(k) + slope*(a1 - x(k))
          b2 = y(k) + slope*(a2 - x(k))
          darea = (a2 - a1)*(b2 + b1)/2.
        END IF


  !  find the area under the trapezoid from a1 to a2

        area = area + darea

  ! go to next point

        k = k+1
        GO TO 40

      END IF

    END IF

  !  calculate the average y after summing the areas in the interval
    yg(i) = area/(xgu - xgl)

  END DO
  !_______________________________________________________________________

  END SUBROUTINE inter2

  !=============================================================================*

  SUBROUTINE inter3(ng,xg,yg, n,x,y, foldin)
  IMPLICIT NONE

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Map input data given on a set of bins onto a different set of target     =*
  !=  bins.                                                                    =*
  !=  The input data are given on a set of bins (representing the integral     =*
  !=  of the input quantity over the range of each bin) and are being matched  =*
  !=  onto another set of bins (target grid).  A typical example would be an   =*
  !=  input data set spcifying the extra-terrestrial flux on wavelength inter- =*
  !=  vals, that has to be matched onto the working wavelength grid.           =*
  !=  The resulting area in a given bin of the target grid is calculated by    =*
  !=  simply adding all fractional areas of the input data that cover that     =*
  !=  particular target bin.                                                   =*
  !=  Some caution should be used near the endpoints of the grids.  If the     =*
  !=  input data do not span the full range of the target grid, the area in    =*
  !=  the "missing" bins will be assumed to be zero.  If the input data extend =*
  !=  beyond the upper limit of the target grid, the user has the option to    =*
  !=  integrate the "overhang" data and fold the remaining area back into the  =*
  !=  last target bin.  Using this option is recommended when re-gridding      =*
  !=  vertical profiles that directly affect the total optical depth of the    =*
  !=  model atmosphere.                                                        =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NG     - INTEGER, number of bins + 1 in the target grid               (I)=*
  !=  XG     - REAL, target grid (e.g. working wavelength grid);  bin i     (I)=*
  !=           is defined as [XG(i),XG(i+1)] (i = 1..NG-1)                     =*
  !=  YG     - REAL, y-data re-gridded onto XG;  YG(i) specifies the        (O)=*
  !=           y-value for bin i (i = 1..NG-1)                                 =*
  !=  N      - INTEGER, number of bins + 1 in the input grid                (I)=*
  !=  X      - REAL, input grid (e.g. data wavelength grid);  bin i is      (I)=*
  !=           defined as [X(i),X(i+1)] (i = 1..N-1)                           =*
  !=  Y      - REAL, input y-data on grid X;  Y(i) specifies the            (I)=*
  !=           y-value for bin i (i = 1..N-1)                                  =*
  !=  FoldIn - Switch for folding option of "overhang" data                 (I)=*
  !=           FoldIn = 0 -> No folding of "overhang" data                     =*
  !=           FoldIn = 1 -> Integerate "overhang" data and fold back into     =*
  !=                         last target bin                                   =*
  !-----------------------------------------------------------------------------*
  INTEGER, INTENT(IN)  :: ng
  INTEGER, INTENT(IN)  :: n
  INTEGER, INTENT(IN)  :: foldin
  REAL, INTENT(IN)     :: xg(ng)
  REAL, INTENT(IN)     :: x(n)
  REAL, INTENT(IN)     :: y(n)
  REAL, INTENT(OUT)    :: yg(ng)

  ! local:
  REAL :: a1, a2, sum
  REAL :: tail
  INTEGER :: jstart, i, j, k
  !_______________________________________________________________________

  ! check whether flag given is legal
  IF ((foldin /= 0) .AND. (foldin /= 1)) THEN
    WRITE(0,*) '>>> ERROR (inter3) <<<  Value for FOLDIN invalid. '
    WRITE(0,*) '                        Must be 0 or 1'
    STOP
  END IF

  ! do interpolation

  jstart = 1

  DO   i = 1, ng - 1

    yg(i) = 0.
    sum = 0.
    j = jstart

    IF (j <= n-1) THEN

      20      CONTINUE

      IF (x(j+1) < xg(i)) THEN
        jstart = j
        j = j+1
        IF (j <= n-1) GO TO 20
      END IF

      25      CONTINUE

      IF ((x(j) <= xg(i+1)) .AND. (j <= n-1)) THEN

        a1 = AMAX1(x(j),xg(i))
        a2 = AMIN1(x(j+1),xg(i+1))

        sum = sum + y(j) * (a2-a1)/(x(j+1)-x(j))
        j = j+1
        GO TO 25

      END IF

      yg(i) = sum

    END IF

  END DO


  ! if wanted, integrate data "overhang" and fold back into last bin

  IF (foldin == 1) THEN

    j = j-1
    a1 = xg(ng)     ! upper limit of last interpolated bin
    a2 = x(j+1)     ! upper limit of last input bin considered

  !        do folding only if grids don't match up and there is more input
    IF ((a2 > a1) .OR. (j+1 < n)) THEN
      tail = y(j) * (a2-a1)/(x(j+1)-x(j))
      DO k = j+1, n-1
        tail = tail + y(k) * (x(k+1)-x(k))
      END DO
      yg(ng-1) = yg(ng-1) + tail
    END IF

  END IF
  !_______________________________________________________________________

  END SUBROUTINE inter3

  !=============================================================================*

  SUBROUTINE inter4(ng,xg,yg, n,x,y, foldin)
  IMPLICIT NONE

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Map input data given on a set of bins onto a different set of target     =*
  !=  bins.                                                                    =*
  !=  The input data are given on a set of bins (representing the integral     =*
  !=  of the input quantity over the range of each bin) and are being matched  =*
  !=  onto another set of bins (target grid).  A typical example would be an   =*
  !=  input data set spcifying the extra-terrestrial flux on wavelength inter- =*
  !=  vals, that has to be matched onto the working wavelength grid.           =*
  !=  The resulting area in a given bin of the target grid is calculated by    =*
  !=  simply adding all fractional areas of the input data that cover that     =*
  !=  particular target bin.                                                   =*
  !=  Some caution should be used near the endpoints of the grids.  If the     =*
  !=  input data do not span the full range of the target grid, the area in    =*
  !=  the "missing" bins will be assumed to be zero.  If the input data extend =*
  !=  beyond the upper limit of the target grid, the user has the option to    =*
  !=  integrate the "overhang" data and fold the remaining area back into the  =*
  !=  last target bin.  Using this option is recommended when re-gridding      =*
  !=  vertical profiles that directly affect the total optical depth of the    =*
  !=  model atmosphere.                                                        =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NG     - INTEGER, number of bins + 1 in the target grid               (I)=*
  !=  XG     - REAL, target grid (e.g. working wavelength grid);  bin i     (I)=*
  !=           is defined as [XG(i),XG(i+1)] (i = 1..NG-1)                     =*
  !=  YG     - REAL, y-data re-gridded onto XG;  YG(i) specifies the        (O)=*
  !=           y-value for bin i (i = 1..NG-1)                                 =*
  !=  N      - INTEGER, number of bins + 1 in the input grid                (I)=*
  !=  X      - REAL, input grid (e.g. data wavelength grid);  bin i is      (I)=*
  !=           defined as [X(i),X(i+1)] (i = 1..N-1)                           =*
  !=  Y      - REAL, input y-data on grid X;  Y(i) specifies the            (I)=*
  !=           y-value for bin i (i = 1..N-1)                                  =*
  !=  FoldIn - Switch for folding option of "overhang" data                 (I)=*
  !=           FoldIn = 0 -> No folding of "overhang" data                     =*
  !=           FoldIn = 1 -> Integerate "overhang" data and fold back into     =*
  !=                         last target bin                                   =*
  !-----------------------------------------------------------------------------*
  INTEGER, INTENT(IN) :: ng
  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(IN) :: foldin
  REAL, INTENT(IN)    :: xg(ng)
  REAL, INTENT(OUT)   :: yg(ng)
  REAL, INTENT(IN)    :: x(n)
  REAL, INTENT(IN)    :: y(n)

  ! local:
  REAL :: a1, a2, sum
  REAL :: tail
  INTEGER :: jstart, i, j, k
  !_______________________________________________________________________

  ! check whether flag given is legal
  IF ((foldin /= 0) .AND. (foldin /= 1)) THEN
    WRITE(0,*) '>>> ERROR (inter3) <<<  Value for FOLDIN invalid. '
    WRITE(0,*) '                        Must be 0 or 1'
    STOP
  END IF

  ! do interpolation

  jstart = 1

  DO   i = 1, ng - 1

    yg(i) = 0.
    sum = 0.
    j = jstart

    IF (j <= n-1) THEN

      20      CONTINUE

      IF (x(j+1) < xg(i)) THEN
        jstart = j
        j = j+1
        IF (j <= n-1) GO TO 20
      END IF

      25      CONTINUE

      IF ((x(j) <= xg(i+1)) .AND. (j <= n-1)) THEN

        a1 = AMAX1(x(j),xg(i))
        a2 = AMIN1(x(j+1),xg(i+1))

        sum = sum + y(j) * (a2-a1)

        j = j+1
        GO TO 25

      END IF

      yg(i) = sum /(xg(i+1)-xg(i))

    END IF

  END DO


  ! if wanted, integrate data "overhang" and fold back into last bin

  IF (foldin == 1) THEN

    j = j-1
    a1 = xg(ng)     ! upper limit of last interpolated bin
    a2 = x(j+1)     ! upper limit of last input bin considered

  !        do folding only if grids don't match up and there is more input
    IF ((a2 > a1) .OR. (j+1 < n)) THEN
      tail = y(j) * (a2-a1)/(x(j+1)-x(j))
      DO k = j+1, n-1
        tail = tail + y(k) * (x(k+1)-x(k))
      END DO
      yg(ng-1) = yg(ng-1) + tail
    END IF

  END IF
  !_______________________________________________________________________

  END SUBROUTINE inter4

  !=============================================================================*

  SUBROUTINE addpnt ( x, y, ld, n, xnew, ynew )
  IMPLICIT NONE

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Add a point <xnew,ynew> to a set of data pairs <x,y>.  x must be in      =*
  !=  ascending order                                                          =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  X    - REAL vector of length LD, x-coordinates                       (IO)=*
  !=  Y    - REAL vector of length LD, y-values                            (IO)=*
  !=  LD   - INTEGER, dimension of X, Y exactly as declared in the calling  (I)=*
  !=         program                                                           =*
  !=  N    - INTEGER, number of elements in X, Y.  On entry, it must be:   (IO)=*
  !=         N < LD.  On exit, N is incremented by 1.                          =*
  !=  XNEW - REAL, x-coordinate at which point is to be added               (I)=*
  !=  YNEW - REAL, y-value of point to be added                             (I)=*
  !-----------------------------------------------------------------------------*
  INTEGER, INTENT(IN)     :: ld
  INTEGER, INTENT(IN OUT) :: n
  REAL, INTENT(IN OUT)    :: x(ld)
  REAL, INTENT(OUT)       :: y(ld)
  REAL, INTENT(IN)        :: xnew
  REAL, INTENT(IN)        :: ynew

  ! local variables

  INTEGER :: insert
  INTEGER :: i

  !-----------------------------------------------------------------------

  ! check n<ld to make sure x will hold another point

  IF (n >= ld) THEN
    WRITE(0,*) '>>> ERROR (ADDPNT) <<<  Cannot expand array '
    WRITE(0,*) '                        All elements used.'
    STOP
  END IF

  insert = 1
  i = 2

  ! check, whether x is already sorted.
  ! also, use this loop to find the point at which xnew needs to be inserted
  ! into vector x, if x is sorted.

  10   CONTINUE
  IF (i < n) THEN
    IF (x(i) < x(i-1)) THEN
      WRITE(0,*) '>>> ERROR (ADDPNT) <<<  x-data must be '//  &
          'in ascending order!'
      STOP
    ELSE
      IF (xnew > x(i)) insert = i + 1
    END IF
    i = i+1
    GO TO 10
  END IF

  ! if <xnew,ynew> needs to be appended at the end, just do so,
  ! otherwise, insert <xnew,ynew> at position INSERT

  IF ( xnew > x(n) ) THEN

    x(n+1) = xnew
    y(n+1) = ynew

  ELSE

  ! shift all existing points one index up

    DO i = n, insert, -1
      x(i+1) = x(i)
      y(i+1) = y(i)
    END DO

  ! insert new point

    x(insert) = xnew
    y(insert) = ynew

  END IF

  ! increase total number of elements in x, y

  n = n+1

  END SUBROUTINE addpnt

  !=============================================================================*

  SUBROUTINE zero1(x,m)
  IMPLICIT NONE

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Initialize all elements of a floating point vector with zero.            =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  X  - REAL, vector to be initialized                                   (O)=*
  !=  M  - INTEGER, number of elements in X                                 (I)=*
  !-----------------------------------------------------------------------------*
  INTEGER, INTENT(IN)                      :: m
  REAL, INTENT(OUT)                        :: x(m)
  INTEGER :: i

  DO  i = 1, m
    x(i) = 0.
  END DO

  END SUBROUTINE zero1

  !=============================================================================*

  SUBROUTINE zero2(x,m,n)
  IMPLICIT NONE

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Initialize all elements of a 2D floating point array with zero.          =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  X  - REAL, array to be initialized                                    (O)=*
  !=  M  - INTEGER, number of elements along the first dimension of X,      (I)=*
  !=       exactly as specified in the calling program                         =*
  !=  N  - INTEGER, number of elements along the second dimension of X,     (I)=*
  !=       exactly as specified in the calling program                         =*
  !-----------------------------------------------------------------------------*

  INTEGER, INTENT(IN)                      :: m
  INTEGER, INTENT(IN)                      :: n
  REAL, INTENT(OUT)                        :: x(m,n)

  ! m,n : dimensions of x, exactly as specified in the calling program

  INTEGER :: i, j


  DO  j = 1, n
    DO  i = 1, m
      x(i,j) = 0.0
    END DO
  END DO

  END SUBROUTINE zero2

  SUBROUTINE odrl(nz,nw,wl,nbl,cAir,dtRl)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Compute Rayleigh optical depths as a function of altitude and wavelength =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
  !=            grid                                                           =*
  !=  Z       - REAL, specified altitude working grid (km)                  (I)=*
  !=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
  !=            wavelength grid                                                =*
  !=  WL      - REAL, vector of lower limits of wavelength intervals in     (I)=*
  !=            working wavelength grid                                        =*
  !=  C       - REAL, number of air molecules per cm^2 at each specified    (O)=*
  !=            altitude layer                                                 =*
  !=  DTRL    - REAL, Rayleigh optical depth at each specified altitude     (O)=*
  !=            and each specified wavelength                                  =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nbl
  INTEGER, INTENT(IN) :: nz(nbl)
  INTEGER, INTENT(IN) :: nw
  REAL, INTENT(IN)    :: wl(kw)
  REAL, INTENT(IN)    :: cAir(nbl,kz)
  REAL, INTENT(OUT)   :: dtRl(nbl,kz,kw)
  REAL,PARAMETER      :: oneover=1.0/1.0e3

  REAL :: srayl, wc, wmicrn, xx
  INTEGER :: iz, iw, iob

  !_______________________________________________________________________

  ! compute Rayleigh cross sections and depths:

  DO   iw = 1, nw - 1
    wc = (wl(iw) + wl(iw+1))/2.

    ! Rayleigh scattering cross section from WMO 1985 (originally from
    ! Nicolet, M., On the molecular scattering in the terrestrial atmosphere:
    ! An empirical formula for its calculation in the homoshpere, Planet.
    ! Space Sci., 32, 1467-1468, 1984.
    xx = 4.04
    wmicrn =  wc*oneover
    IF( wmicrn <= 0.55) xx = 3.6772 + 0.389*wmicrn + 0.09426/wmicrn
    srayl = 4.02E-28/(wmicrn)**xx
    ! alternate (older) expression from
    ! Frohlich and Shaw, Appl.Opt. v.11, p.1773 (1980).
    !     xx = 3.916 + 0.074*wmicrn + 0.050/wmicrn
    !     srayl(iw) = 3.90e-28/(wmicrn)**xx
    DO   iz = 1, zLimit
       DO iob=1,nbl
          IF(iz>(nz(iob)-1)) CYCLE
          dtRl(iob,iz,iw) = cAir(iob,iz)*srayl
       END DO
    END DO

  END DO
  !_______________________________________________________________________

  END SUBROUTINE odrl

  ! subroutines used for calculation of quantum yields for
  ! various photoreactions:
  !     qyacet - q.y. for acetone, based on Blitz et al. (2004)

  !*******************************************************************************

  SUBROUTINE qyacet(w, t, m, fco, fac)
  ! Compute acetone quantum yields according to the parameterization of:
  ! Blitz, M. A., D. E. Heard, M. J. Pilling, S. R. Arnold, and M. P. Chipperfield
  !       (2004), Pressure and temperature-dependent quantum yields for the
  !       photodissociation of acetone between 279 and 327.5 nm, Geophys.
  !       Res. Lett., 31, L06111, doi:10.1029/2003GL018793.
  ! input:
  ! w = wavelength, nm
  ! T = temperature, K
  ! m = air number density, molec. cm-3
  ! output
  ! fco = quantum yield for product CO
  ! fac = quantum yield for product CH3CO (acetyl radical)
  IMPLICIT NONE

  REAL, INTENT(IN)  :: w
  REAL, INTENT(IN)  :: t
  REAL, INTENT(IN)  :: m
  REAL, INTENT(OUT) :: fco
  REAL, INTENT(OUT) :: fac

  REAL :: a0, a1, a2, a3, a4
  REAL :: b0, b1, b2, b3, b4
  REAL :: c3
  REAL :: ca0, ca1, ca2, ca3, ca4

  !** set out-of-range values:
  ! use low pressure limits for shorter wavelengths
  ! set to zero beyound 327.5

  IF(w < 279.) THEN
    fco = 0.05
    fac = 0.95
    RETURN
  END IF

  IF(w > 327.5 ) THEN
    fco = 0.
    fac = 0.
    RETURN
  END IF

  !** CO (carbon monoxide) quantum yields:

  a0 = 0.350 * (t/295.)**(-1.28)
  b0 = 0.068 * (t/295.)**(-2.65)
  ca0 = EXP(b0*(w - 248.)) * a0 / (1. - a0)

  fco = 1. / (1 + ca0)

  !** CH3CO (acetyl radical) quantum yields:

  IF(w >= 279. .AND. w < 302.) THEN

    a1 = 1.600E-19 * (t/295.)**(-2.38)
    b1 = 0.55E-3   * (t/295.)**(-3.19)
    ca1 = a1 * EXP(-b1*((1.e7/w) - 33113.))

    fac = (1. - fco) / (1 + ca1 * m)

  END IF

  IF(w >= 302. .AND. w < 327.5) THEN

    a2 = 1.62E-17 * (t/295.)**(-10.03)
    b2 = 1.79E-3  * (t/295.)**(-1.364)
    ca2 = a2 * EXP(-b2*((1.e7/w) - 30488.))


    a3 = 26.29   * (t/295.)**(-6.59)
    b3 = 5.72E-7 * (t/295.)**(-2.93)
    c3 = 30006   * (t/295.)**(-0.064)
    ca3 = a3 * EXP(-b3*((1.e7/w) - c3)**2)


    a4 = 1.67E-15 * (t/295.)**(-7.25)
    b4 = 2.08E-3  * (t/295.)**(-1.16)
    ca4 = a4 * EXP(-b4*((1.e7/w) - 30488.))

    fac = (1. - fco) * (1. + ca3 + ca4 * m) /  &
        ((1. + ca3 + ca2 * m)*(1. + ca4 * m))

  END IF

  END SUBROUTINE qyacet


  SUBROUTINE rdetfl(nw,wl)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Read and re-grid extra-terrestrial flux data.                            =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  F      - REAL, spectral irradiance at the top of the atmosphere at    (O)=*
  !=           each specified wavelength                                       =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE


  INTEGER, INTENT(IN) :: nw    ! input: (wavelength grid)
  REAL, INTENT(IN)    :: wl(kw)

  INTEGER, PARAMETER :: kdata=20000
  INTEGER :: iw

  ! work arrays for input data files:

  CHARACTER (LEN=200) :: fil
  REAL :: x1(kdata)
  REAL :: y1(kdata)
  INTEGER :: nhead, n, i, ierr
  REAL :: dum

  ! data gridded onto wl(kw) grid:

  REAL :: yg1(kw)
  REAL :: yg2(kw)
  REAL :: yg3(kw)


  REAL, PARAMETER :: hc = 6.62E-34 * 2.998E8

  ! simple files are read and interpolated here in-line. Reading of
  ! more complex files may be done with longer code in a read#.f subroutine.

  SELECT CASE (msun)
    CASE(1)
       fil = trim(files(2)%fileName)
       !WRITE(kout,*) fil
       OPEN(UNIT=kin,FILE=trim(fil),STATUS='old')
       nhead = 3
       n =121
       DO i = 1, nhead
         READ(kin,*)
       END DO
       DO i = 1, n
         READ(kin,*) x1(i), y1(i)
       END DO
       CLOSE (kin)
       CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
       CALL addpnt(x1,y1,kdata,n,	   0.,0.)
       CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
       CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
       CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
       IF (ierr /= 0) THEN
         WRITE(*,*) ierr, fil
         STOP
       END IF
       DO iw = 1, nw-1
         f(iw) = yg1(iw)
       END DO
    CASE(2)
       fil = trim(files(3)%fileName)
       !WRITE(kout,*) fil
       OPEN(UNIT=kin,FILE=trim(fil),STATUS='old')
       nhead = 3
       n = 4327
       DO i = 1, nhead
   	 READ(kin,*)
       END DO
       DO i = 1, n
   	 READ(kin,*) x1(i), y1(i)
       END DO
       CLOSE (kin)
       CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
       CALL addpnt(x1,y1,kdata,n,	   0.,0.)
       CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
       CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
       CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
       IF (ierr /= 0) THEN
   	 WRITE(*,*) ierr, fil
   	 STOP
       END IF
       DO iw = 1, nw-1
   	 f(iw) = yg1(iw)
       END DO
    CASE (3)
       fil = trim(files(4)%fileName)
       !WRITE(kout,*) fil
       OPEN(UNIT=kin,FILE=trim(fil),STATUS='old')
       nhead = 6
       n = 14980
       DO i = 1, nhead
   	 READ(kin,*)
       END DO
       DO i = 1, n
   	 READ(kin,*) x1(i), y1(i)
       END DO
       CLOSE (kin)
       CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
       CALL addpnt(x1,y1,kdata,n,	   0.,0.)
       CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
       CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
       CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
       IF (ierr /= 0) THEN
   	 WRITE(*,*) ierr, fil
   	 STOP
       END IF
       DO iw = 1, nw-1
   	 f(iw) = yg1(iw)
       END DO
    CASE (4)
       fil = trim(files(5)%fileName)
       !WRITE(kout,*) fil
       OPEN(UNIT=kin,FILE=trim(fil),STATUS='old')
       nhead = 8
       n = 1260
       DO i = 1, nhead
   	 READ(kin,*)
       END DO
       DO i = 1, n
   	 READ(kin,*) x1(i), y1(i)
       END DO
       CLOSE (kin)
       CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
       CALL addpnt(x1,y1,kdata,n,	   0.,0.)
       CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
       CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
       CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
       IF (ierr /= 0) THEN
   	 WRITE(*,*) ierr, fil
   	 STOP
       END IF
       DO iw = 1, nw-1
   	 f(iw) = yg1(iw)
       END DO
    CASE (5)
       ! unofficial - do not use
       fil = trim(files(6)%fileName)
       !WRITE(kout,*) fil
       OPEN(UNIT=kin,FILE=trim(fil),STATUS='old')
       nhead = 11
       n = 2047
       DO i = 1, nhead
     	 READ(kin,*)
       END DO
       DO i = 1, n
     	 READ(kin,*) x1(i), y1(i)
       END DO
       CLOSE (kin)
       CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
       CALL addpnt(x1,y1,kdata,n,	   0.,0.)
       CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
       CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
       CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
       IF (ierr /= 0) THEN
     	 WRITE(*,*) ierr, fil
     	 STOP
       END IF
       DO iw = 1, nw-1
     	 f(iw) = yg1(iw)
       END DO
     CASE (6)
       ! unofficial - do not use
       fil = trim(files(7)%fileName)
       !WRITE(kout,*) fil
       OPEN(UNIT=kin,FILE=trim(fil),STATUS='old')
       nhead = 3
       n = 1200
       DO i = 1, nhead
   	 READ(kin,*)
       END DO
       DO i = 1, n
   	 READ(kin,*) x1(i), y1(i)
   	 y1(i) = y1(i)* 1.e-3
       END DO
       CLOSE (kin)
       CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
       CALL addpnt(x1,y1,kdata,n,	   0.,0.)
       CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
       CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
       CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
       IF (ierr /= 0) THEN
   	 WRITE(*,*) ierr, fil
   	 STOP
       END IF
       DO iw = 1, nw-1
          f(iw) = yg1(iw)
       END DO
    CASE (7)
       fil = trim(files(8)%fileName)
       !WRITE(kout,*) fil
       OPEN(UNIT=kin,FILE=trim(fil),STATUS='old')
       nhead = 11
       n = 496
       DO i = 1, nhead
   	 READ(kin,*)
       END DO
       DO i = 1, n
   	 READ(kin,*) dum, y1(i)
   	 IF (dum < 630.0) x1(i) = dum - 0.5
   	 IF (dum > 630.0 .AND. dum < 870.0) x1(i) = dum - 1.0
   	 IF (dum > 870.0) x1(i) = dum - 2.5
   	 y1(i) = y1(i) * 1.e4 * hc / (dum * 1.e-9)
       END DO
       CLOSE (kin)
       x1(n+1) = x1(n) + 2.5
       DO i = 1, n
   	 y1(i) = y1(i) * (x1(i+1)-x1(i))
       END DO
       CALL inter3(nw,wl,yg2,n+1,x1,y1,0)
       DO iw = 1, nw-1
   	 yg1(iw) = yg1(iw) / (wl(iw+1)-wl(iw))
       END DO
       DO iw = 1, nw-1
         f(iw) = yg1(iw)
       END DO
    CASE (8)
       nhead = 5
       fil = trim(files(9)%fileName)
       !WRITE(kout,*) fil
       OPEN(UNIT=kin,FILE=trim(fil),STATUS='old')
       nhead = 13
       n = 5160
       DO i = 1, nhead
   	 READ(kin,*)
       END DO
       DO i = 1, n
   	 READ(kin,*) x1(i), y1(i)
   	 y1(i) = y1(i) * 1.e-3
       END DO
       CLOSE (kin)
       CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
       CALL addpnt(x1,y1,kdata,n,	   0.,0.)
       CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
       CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
       CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
       IF (ierr /= 0) THEN
   	 WRITE(*,*) ierr, fil
   	 STOP
       END IF
       DO iw = 1, nw-1
          f(iw) = yg1(iw)
       END DO
    CASE (9)
       fil = trim(files(10)%filename)
       !WRITE(kout,*) fil
       OPEN(UNIT=kin,FILE=trim(fil),STATUS='old')
       nhead = 2
       n = 302
       DO i = 1, nhead
   	 READ(kin,*)
       END DO
       DO i = 1, n
   	 READ(kin,*) x1(i), y1(i)
       END DO
       CLOSE (kin)
       CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
       CALL addpnt(x1,y1,kdata,n,	   0.,0.)
       CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
       CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
       CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
       IF (ierr /= 0) THEN
   	 WRITE(*,*) ierr, fil
   	 STOP
       END IF
       DO iw = 1, nw-1
   	 f(iw) = yg1(iw)
       END DO
    CASE (10)
       !!WRITE(kout,*) 'datae1/sun/susim_hi.flx'
       CALL read1(nw,wl,yg1)
       DO iw = 1, nw-1
         f(iw) = yg1(iw)
       END DO
    CASE (11)
       !!WRITE(kout,*) 'datae1/sun/wmo85.flx'
       CALL read2(nw,wl,yg1)
       DO iw = 1, nw-1
         f(iw) = yg1(iw)
       END DO
    CASE (12)
       !!WRITE(kout,*) 'datae1/sun/susim_hi.flx'
       CALL read1(nw,wl,yg1)
       fil = trim(files(13)%fileName)
       !WRITE(kout,*) fil
       OPEN(UNIT=kin,FILE=trim(fil),STATUS='old')
       nhead = 11
       n = 496
       DO i = 1, nhead
   	 READ(kin,*)
       END DO
       DO i = 1, n
   	 READ(kin,*) dum, y1(i)
   	 IF (dum < 630.0) x1(i) = dum - 0.5
   	 IF (dum > 630.0 .AND. dum < 870.0) x1(i) = dum - 1.0
   	 IF (dum > 870.0) x1(i) = dum - 2.5
   	 y1(i) = y1(i) * 1.e4 * hc / (dum * 1.e-9)
       END DO
       CLOSE (kin)
       x1(n+1) = x1(n) + 2.5
       DO i = 1, n
   	 y1(i) = y1(i) * (x1(i+1)-x1(i))
       END DO
       CALL inter3(nw,wl,yg2,n+1,x1,y1,0)
       DO iw = 1, nw-1
         yg2(iw) = yg2(iw) / (wl(iw+1)-wl(iw))
       END DO
       DO iw = 1, nw-1
         IF (wl(iw) > 350.) THEN
            f(iw) = yg2(iw)
         ELSE
            f(iw) = yg1(iw)
         END IF
       END DO
    CASE (13)
       !!WRITE(kout,*) 'datae1/sun/susim_hi.flx'
       CALL read1(nw,wl,yg1)
       nhead = 5
       fil = trim(files(14)%fileName)
       !WRITE(kout,*) fil
       OPEN(UNIT=kin,FILE=trim(fil),STATUS='old')
       n = 5160
       DO i = 1, nhead
         READ(kin,*)
       END DO
       DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.e-3
       END DO
       CLOSE (kin)
       CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
       CALL addpnt(x1,y1,kdata,n,	   0.,0.)
       CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
       CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
       CALL inter2(nw,wl,yg2,n,x1,y1,ierr)
       IF (ierr /= 0) THEN
         WRITE(*,*) ierr, fil
         STOP
       END IF
       fil = trim(files(15)%fileName)
       !WRITE(kout,*) fil
       OPEN(UNIT=kin,FILE=trim(fil),STATUS='old')
       nhead = 11
       n = 496
       DO i = 1, nhead
   	 READ(kin,*)
       END DO
       DO i = 1, n
   	 READ(kin,*) dum, y1(i)
   	 IF (dum < 630.0) x1(i) = dum - 0.5
   	 IF (dum > 630.0 .AND. dum < 870.0) x1(i) = dum - 1.0
   	 IF (dum > 870.0) x1(i) = dum - 2.5
   	 y1(i) = y1(i) * 1.e4 * hc / (dum * 1.e-9)
       END DO
       CLOSE (kin)

       x1(n+1) = x1(n) + 2.5
       CALL inter4(nw,wl,yg3,n+1,x1,y1,0)

       DO iw = 1, nw-1
   	 IF (wl(iw) < 150.01) THEN
   	   f(iw) = yg1(iw)
   	 ELSE IF ((wl(iw) >= 150.01) .AND. wl(iw) <= 400.) THEN
           f(iw) = yg2(iw)
         ELSE IF (wl(iw) > 400.) THEN
           f(iw) = yg3(iw)
         END IF
       END DO
    END SELECT

  END SUBROUTINE rdetfl


  SUBROUTINE read1(nw,wl,f)
  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Read extra-terrestrial flux data.  Re-grid data to match specified       =*
  !=  working wavelength grid.                                                 =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  F      - REAL, spectral irradiance at the top of the atmosphere at    (O)=*
  !=           each specified wavelength                                       =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE


  INTEGER, INTENT(IN) :: nw    ! input: (wavelength grid)
  REAL, INTENT(IN)    :: wl(kw)
  REAL, INTENT(OUT)   :: f(kw) ! output: (extra terrestrial solar flux)

  ! local:

  REAL :: lambda_hi(10000),irrad_hi(10000)
  REAL :: lambda
  INTEGER :: ierr
  INTEGER :: i, j, n
  CHARACTER (LEN=200) :: fil

  !_______________________________________________________________________

  !****** SUSIM irradiance
  !_______________________________________________________________________
  ! VanHoosier, M. E., J.-D. F. Bartoe, G. E. Brueckner, and
  ! D. K. Prinz, Absolute solar spectral irradiance 120 nm -
  ! 400 nm (Results from the Solar Ultraviolet Spectral Irradiance
  ! Monitor - SUSIM- Experiment on board Spacelab 2),
  ! Astro. Lett. and Communications, 1988, vol. 27, pp. 163-168.
  !     SUSIM SL2 high resolution (0.15nm) Solar Irridance data.
  !     Irradiance values are given in milliwatts/m^2/nanomenters
  !     and are listed at 0.05nm intervals.  The wavelength given is
  !     the center wavelength of the 0.15nm triangular bandpass.
  !     Normalized to 1 astronomical unit.
  !  DATA for wavelengths > 350 nm are unreliable
  ! (Van Hoosier, personal communication, 1994).
  !_______________________________________________________________________

  !* high resolution

  fil = trim(files( 11)%fileName)
  OPEN(UNIT=kin,FILE=trim(fil),STATUS='old')
  DO   i = 1, 7
    READ(kin,*)
  END DO
  DO   i = 1, 559
    READ(kin,*)lambda,(irrad_hi(10*(i-1)+j), j=1, 10)
  END DO
  CLOSE (kin)

  ! compute wavelengths, convert from mW to W

  n = 559*10
  DO   i = 1, n
    lambda_hi(i)=120.5 + FLOAT(i-1)*.05
    irrad_hi(i) = irrad_hi(i)  /  1000.
  END DO
  !_______________________________________________________________________

  CALL addpnt(lambda_hi,irrad_hi,10000,n, lambda_hi(1)*(1.-deltax),0.)
  CALL addpnt(lambda_hi,irrad_hi,10000,n,                 0.,0.)
  CALL addpnt(lambda_hi,irrad_hi,10000,n, lambda_hi(n)*(1.+deltax),0.)
  CALL addpnt(lambda_hi,irrad_hi,10000,n,              1.e38,0.)
  CALL inter2(nw,wl,f,n,lambda_hi,irrad_hi,ierr)
  IF (ierr /= 0) THEN
    WRITE(*,*) ierr, fil
    STOP
  END IF

  END SUBROUTINE read1

  !=============================================================================*

  SUBROUTINE read2(nw,wl,f)
  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Read extra-terrestrial flux data.  Re-grid data to match specified       =*
  !=  working wavelength grid.                                                 =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  F      - REAL, spectral irradiance at the top of the atmosphere at    (O)=*
  !=           each specified wavelength                                       =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nw     ! input: (wavelength grid)
  REAL, INTENT(IN)    :: wl(kw)
  REAL, INTENT(OUT)   :: f(kw)  ! output: (extra terrestrial solar flux)

  REAL :: yg(kw)
  INTEGER :: iw

  ! local:

  REAL :: x1(1000), y1(1000)
  REAL :: x2(1000)
  REAL :: x3(1000)
  INTEGER :: i, n
  REAL :: dum
  INTEGER :: idum

  !_______________________________________________________________________

  !********WMO 85 irradiance

  OPEN(UNIT=kin,FILE=trim(files(12)%fileName),STATUS='old')
  DO   i = 1, 3
    READ(kin,*)
  END DO
  n = 158
  DO   i = 1, n
    READ(kin,*) idum, x1(i),x2(i),y1(i), dum, dum, dum
    x3(i) = 0.5 * (x1(i) + x2(i))

  ! average value needs to be calculated only if inter2 is
  ! used to interpolate onto wavelength grid (see below)
  !        y1(i) =  y1(i) / (x2(i) - x1(i))

  END DO
  CLOSE (kin)

  x1(n+1) = x2(n)

  ! inter2: INPUT : average value in each bin
  !         OUTPUT: average value in each bin
  ! inter3: INPUT : total area in each bin
  !         OUTPUT: total area in each bin

  CALL inter3(nw,wl,yg, n+1,x1,y1,0)
  !      CALL inter2(nw,wl,yg,n,x3,y1,ierr)

  DO    iw = 1, nw-1
  ! from quanta s-1 cm-2 bin-1 to  watts m-2 nm-1
  ! 1.e4 * ([hc =] 6.62E-34 * 2.998E8)/(wc*1e-9)

  ! the scaling by bin width needs to be done only if
  ! inter3 is used for interpolation

    yg(iw) = yg(iw) / (wl(iw+1)-wl(iw))
    f(iw) = yg(iw) * 1.e4 * (6.62E-34 * 2.998E8) /  &
        ( 0.5 * (wl(iw+1)+wl(iw)) * 1.e-9)

  END DO

  END SUBROUTINE read2

  SUBROUTINE rdno2xs(nw,wl,no2xs)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Read NO2 molecular absorption cross section.  Re-grid data to match      =*
  !=  specified wavelength working grid.                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  NO2XS  - REAL, molecular absoprtion cross section (cm^2) of NO2 at    (O)=*
  !=           each specified wavelength                                       =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE


  INTEGER, INTENT(IN) :: nw
  REAL, INTENT(IN)    :: wl(kw)
  REAL, INTENT(OUT)   :: no2xs(kw)

  INTEGER, PARAMETER :: kdata=1000

  ! input: (altitude working grid)

  ! local:
  REAL :: x1(kdata)
  REAL :: y1(kdata)
  REAL :: yg(kw)
  REAL :: dum
  INTEGER :: ierr
  INTEGER :: i, l, n, idum
  CHARACTER (LEN=200) :: fil
  !_______________________________________________________________________

  !************ absorption cross sections:
  !     measurements by:
  ! Davidson, J. A., C. A. Cantrell, A. H. McDaniel, R. E. Shetter,
  ! S. Madronich, and J. G. Calvert, Visible-ultraviolet absorption
  ! cross sections for NO2 as a function of temperature, J. Geophys.
  ! Res., 93, 7105-7112, 1988.
  !  Values at 273K from 263.8 to 648.8 nm in approximately 0.5 nm intervals

  fil = trim(files(19)%fileName)

  OPEN(UNIT=kin,FILE=trim(fil),STATUS='old')
  n = 750
  DO i = 1, n
    READ(kin,*) x1(i), y1(i), dum, dum, idum
  END DO
  CLOSE(kin)

  CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
  CALL addpnt(x1,y1,kdata,n,          0.,0.)
  CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
  CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
  CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  IF (ierr /= 0) THEN
    WRITE(*,*) ierr, fil
    STOP
  END IF

  DO   l = 1, nw-1
    no2xs(l) = yg(l)
  END DO

  END SUBROUTINE rdno2xs

  SUBROUTINE rdo2xs(nw,wl,o2xs1)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Compute equivalent O2 cross section, except                              =*
  !=  the SR bands and the Lyman-alpha line.                                   =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:
  !=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
  !=            wavelength grid                                                =*
  !=  WL      - REAL, vector of lower limits of wavelength intervals in     (I)=*
  !=            working wavelength grid
  !=            vertical layer at each specified wavelength                    =*
  !=  O2XS1   - REAL, O2 molecular absorption cross section                    =*
  !=
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nw
  REAL, INTENT(IN)    :: wl(kw)
  REAL, INTENT(OUT)   :: o2xs1(kw) ! Output O2 xsect, temporary, will be
                                   ! over-written in Lyman-alpha and
                                   ! Schumann-Runge wavelength bands.

  INTEGER :: i, n
  INTEGER, PARAMETER :: kdata = 200
  REAL :: x1(kdata), y1(kdata)
  REAL :: x, y
  INTEGER :: ierr

  !-----------------------------------------------------

  ! Read O2 absorption cross section data:
  !  116.65 to 203.05 nm = from Brasseur and Solomon 1986
  !  205 to 240 nm = Yoshino et al. 1988

  ! Note that subroutine la_srb.f will over-write values in the spectral regions
  !   corresponding to:
  ! - Lyman-alpha (LA: 121.4-121.9 nm, Chabrillat and Kockarts parameterization)
  ! - Schumann-Runge bands (SRB: 174.4-205.8 nm, Koppers parameteriaztion)

  n = 0

  OPEN(UNIT=kin,FILE=trim(files(16)%fileName))
  DO i = 1, 7
    READ(kin,*)
  END DO
  DO i = 1, 78
    READ(kin,*) x, y
    IF (x <= 204.) THEN
      n = n + 1
      x1(n) = x
      y1(n) = y
    END IF
  END DO
  CLOSE(kin)

  OPEN(UNIT=kin,FILE=trim(files(17)%fileName),STATUS='old')
  DO i = 1, 8
    READ(kin,*)
  END DO
  DO i = 1, 36
    n = n + 1
    READ(kin,*) x, y
    y1(n) = y*1.e-24
    x1(n) = x
  END DO
  CLOSE (kin)

  ! Add termination points and interpolate onto the
  !  user grid (set in subroutine gridw):

  CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
  CALL addpnt(x1,y1,kdata,n,0.               ,y1(1))
  CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
  CALL addpnt(x1,y1,kdata,n,              1.e+38,0.)
  CALL inter2(nw,wl,o2xs1, n,x1,y1, ierr)
  IF (ierr /= 0) THEN
    WRITE(*,*) ierr, 'O2 -> O + O'
    STOP
  END IF

  END SUBROUTINE rdo2xs

  SUBROUTINE rdo3xs(nw,wl,nz,nbl,tLay,o3xs)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Read ozone molecular absorption cross section.  Re-grid data to match    =*
  !=  specified wavelength working grid.                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  O3XS   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (O)=*
  !=           each specified wavelength (WMO value at 273)                    =*
  !=  S226   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (O)=*
  !=           each specified wavelength (value from Molina and Molina at 226K)=*
  !=  S263   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (O)=*
  !=           each specified wavelength (value from Molina and Molina at 263K)=*
  !=  S298   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (O)=*
  !=           each specified wavelength (value from Molina and Molina at 298K)=*
  !=  opt    - if opt=0 read files
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nw
  INTEGER, INTENT(IN) :: nbl
  INTEGER, INTENT(IN) :: nz(nbl)
  REAL, INTENT(IN)    :: wl(kw)
  REAL, INTENT(IN)    :: tLay(nbl,kz)
  REAL, INTENT(OUT)   :: o3xs(nbl,kz,kw)

  INTEGER :: iob
  ! output:
  ! ozone absorption cross section at three different
  ! temperatures: 226, 263, 298 Kelvin.  Can interpolate
  ! to different temperatures. Units are cm2 molecule-1

  SELECT CASE (mOption(1))
     CASE (1)
        DO iob=1,nbl
           CALL o3xs_mm(nw,wl,nz(iob),o3xs(iob,:,:),tLay(iob,:))
	END DO
     CASE (2)
        DO iob=1,nbl
           CALL o3xs_mal(nw,wl,nz(iob),o3xs(iob,:,:),tLay(iob,:))
	END DO
     CASE (3)
        DO iob=1,nbl
	   CALL o3xs_bass(nw,wl,nz(iob),o3xs(iob,:,:),tLay(iob,:))
	END DO
  END SELECT

  END SUBROUTINE rdo3xs

  !=============================================================================*

  SUBROUTINE o3xs_mm(nw,wl,nz,xs,temp)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Interpolate the O3 cross section                                         =*
  !=  Combined data from WMO 85 Ozone Assessment (use 273K value from          =*
  !=  175.439-847.5 nm) and:                                                   =*
  !=  For Hartley and Huggins bands, use temperature-dependent values from     =*
  !=  Molina, L. T., and M. J. Molina, Absolute absorption cross sections      =*
  !=  of ozone in the 185- to 350-nm wavelength range, J. Geophys. Res.,       =*
  !=  vol. 91, 14501-14508, 1986.                                              =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  XS     - REAL, cross section (cm^2) for O3                            (O)=*
  !=           at each defined wavelength and each defined altitude level      =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE


  INTEGER, INTENT(IN)  :: nw
  INTEGER, INTENT(IN)  :: nz
  REAL, INTENT(IN)     :: wl(kw)
  REAL, INTENT(IN)     :: temp(kz)
  REAL, INTENT(OUT)    :: xs(kz,kw)

  INTEGER :: iw
  INTEGER :: iz
  REAL,PARAMETER :: div1=1.0/(263.-226.)
  REAL,PARAMETER :: div2=1.0/(298.-263.)

  DO  iw = 1, nw-1
    DO  iz = 1, nz
      xs(iz,iw) = mm_o3xs(iw)
      IF ( wl(iw) > 240.5  .AND. wl(iw+1) < 350. ) THEN
        IF (temp(iz) < 263.) THEN
          xs(iz,iw) = s226(iw) + (s263(iw)-s226(iw))* &
         	   (temp(iz)-226.)*div1
        ELSE
          xs(iz,iw) = s263(iw) + (s298(iw)-s263(iw))* &
	           (temp(iz)-263.)*div2
        END IF
      END IF
    END DO
  END DO

  END SUBROUTINE o3xs_mm

  !=============================================================================*

  SUBROUTINE o3xs_mal(nw,wl,nz, xs,temp)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Read and interpolate the O3 cross section                                =*
  !=  Combined data from WMO 85 Ozone Assessment (use 273K value from          =*
  !=  175.439-847.5 nm) and:                                                   =*
  !=  For Hartley and Huggins bands, use temperature-dependent values from     =*
  !=  Malicet et al., J. Atmos. Chem.  v.21, pp.263-273, 1995.                 =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  XS     - REAL, cross section (cm^2) for O3                            (O)=*
  !=           at each defined wavelength and each defined altitude level      =*
  !-----------------------------------------------------------------------------*

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nw
  INTEGER, INTENT(IN) :: nz
  REAL, INTENT(IN)    :: wl(kw)
  REAL, INTENT(IN)     :: temp(kz)
  REAL, INTENT(OUT)   :: xs(kz,kw)

  INTEGER :: iw
  INTEGER :: iz
  INTEGER, PARAMETER :: kdata = 16000

  ! assign:
  DO  iw = 1, nw-1
    DO  iz = 1, nz

      xs(iz,iw) = mm_o3xs(iw)
      IF ( wl(iw) > 195.  .AND. wl(iw+1) < 345. ) THEN
        IF(temp(iz) >= 243.) THEN
          xs(iz,iw) = s243(iw) + (s295(iw)-s243(iw))* &
	     (temp(iz)-243.)/(295.-243.)
        END IF
        IF(temp(iz) < 254. .AND. &
	   temp(iz) >= 228.) THEN
          xs(iz,iw) = s228(iw) + (s243(iw)-s228(iw))* &
	    (temp(iz)-228.)/(243.-228.)
        END IF
        IF(temp(iz) < 228.) THEN
          xs(iz,iw) = s218(iw) + (s228(iw)-s218(iw))* &
	    (temp(iz)-218.)/(228.-218.)
        END IF
      END IF

    END DO

  END DO

  END SUBROUTINE o3xs_mal

  !=============================================================================*

  SUBROUTINE o3xs_bass(nw,wl,nz, xs,temp)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Read and interpolate the O3 cross section                                =*
  !=  Combined data from WMO 85 Ozone Assessment (use 273K value from          =*
  !=  175.439-847.5 nm) and:                                                   =*
  !=  For Hartley and Huggins bands, use temperature-dependent values from     =*
  !=  Bass
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  XS     - REAL, cross section (cm^2) for O3                            (O)=*
  !=           at each defined wavelength and each defined altitude level      =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nz
  REAL, INTENT(IN)                         :: wl(kw)
  REAL, INTENT(IN)     :: temp(kz)
  REAL, INTENT(OUT)                        :: xs(kz,kw)

  INTEGER :: iw
  INTEGER :: iz
  INTEGER, PARAMETER :: kdata = 2000
  REAL :: tc

  DO  iw = 1, nw-1
    DO  iz = 1, nz

        tc = temp(iz) - 273.15

        xs(iz,iw) = mm_o3xs(iw)
        IF ( wl(iw) > 245. .AND. wl(iw+1) < 341. ) THEN

    	xs(iz,iw) = 1.e-20 * (c0(iw) + c1(iw)*tc + c2(iw)*tc*tc)

        END IF

      END DO

  END DO

  END SUBROUTINE o3xs_bass

  !=============================================================================*

  SUBROUTINE rdso2xs(nw,wl,so2xs)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Read SO2 molecular absorption cross section.  Re-grid data to match      =*
  !=  specified wavelength working grid.                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  SO2XS  - REAL, molecular absoprtion cross section (cm^2) of SO2 at    (O)=*
  !=           each specified wavelength                                       =*
  !-----------------------------------------------------------------------------*
  !=  EDIT HISTORY:                                                            =*
  !=  02/97  Changed offset for grid-end interpolation to relative number      =*
  !=         (x * (1 +- deltax)                                                =*
  !-----------------------------------------------------------------------------*
  != This program is free software;  you can redistribute it and/or modify     =*
  != it under the terms of the GNU General Public License as published by the  =*
  != Free Software Foundation;  either version 2 of the license, or (at your   =*
  != option) any later version.                                                =*
  != The TUV package is distributed in the hope that it will be useful, but    =*
  != WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
  != LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
  != License for more details.                                                 =*
  != To obtain a copy of the GNU General Public License, write to:             =*
  != Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
  !-----------------------------------------------------------------------------*
  != To contact the authors, please mail to:                                   =*
  != Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
  != send email to:  sasha@ucar.edu                                            =*
  !-----------------------------------------------------------------------------*
  != Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nw
  REAL, INTENT(IN)    :: wl(kw)
  REAL, INTENT(OUT)   :: so2xs(kw)


  INTEGER, PARAMETER :: kdata=1000

  ! local:
  REAL :: x1(kdata)
  REAL :: y1(kdata)
  REAL :: yg(kw)
  INTEGER :: ierr
  INTEGER :: i, l, n!, idum
  CHARACTER (LEN=200) :: fil
  !_______________________________________________________________________

  !************ absorption cross sections:
  ! SO2 absorption cross sections from J. Quant. Spectrosc. Radiat. Transfer
  ! 37, 165-182, 1987, T. J. McGee and J. Burris Jr.
  ! Angstrom vs. cm2/molecule, value at 221 K

  fil = 'data/McGee87'
  OPEN(UNIT=kin,FILE=trim(files(18)%fileName),STATUS='old')
  DO   i = 1,3
    READ(kin,*)
  END DO
  !      n = 681
  n = 704
  DO   i = 1, n
    READ(kin,*) x1(i), y1(i)
    x1(i) = x1(i)/10.
  END DO
  CLOSE (kin)

  CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
  CALL addpnt(x1,y1,kdata,n,          0.,0.)
  CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
  CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
  CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  IF (ierr /= 0) THEN
    WRITE(*,*) ierr, fil
    STOP
  END IF

  DO   l = 1, nw-1
    so2xs(l) = yg(l)
  END DO

  !_______________________________________________________________________

  END SUBROUTINE rdso2xs


  SUBROUTINE rtlink(nstr,nz,iw,nbl,sza,albedo,dtCld,omCld,gCld, &
                    dtAer,omAer,gAer,dtO2,dtO3,eDir,eDn,eUp, &
		    fDir,fDn,fUp,dtRl,dtSo2,dtNo2,nid,dsdh)
     IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nbl
  INTEGER, INTENT(IN)  :: nstr
  INTEGER, INTENT(IN)  :: nz(nbl)
  INTEGER, INTENT(IN)  :: nid(nbl,0:kz)
  INTEGER, INTENT(IN)  :: iw
  REAL   , INTENT(IN)  :: sza(nbl)
  REAL   , INTENT(IN)  :: albedo(nbl,kw)
  REAL   , INTENT(IN)  :: dtCld(nbl,kz,kw)
  REAL   , INTENT(IN)  :: omCld(nbl,kz,kw)
  REAL   , INTENT(IN)  :: gCld(nbl,kz,kw)
  REAL   , INTENT(IN)  :: dtAer(nbl,kz,kw)
  REAL   , INTENT(IN)  :: omAer(nbl,kz,kw)
  REAL   , INTENT(IN)  :: gAer(nbl,kz,kw)
  REAL   , INTENT(IN)  :: dtO2(nbl,kz,kw)
  REAL   , INTENT(IN)  :: dtO3(nbl,kz,kw)
  REAL   , INTENT(IN)  :: dtRl(nbl,kz,kw)
  REAL   , INTENT(IN)  :: dtSo2(nbl,kz,kw)
  REAL   , INTENT(IN)  :: dtNo2(nbl,kz,kw)
  REAL   , INTENT(IN)  :: dsdh(nbl,0:kz,kz)
  REAL   , INTENT(OUT) :: eDir(nbl,kz)
  REAL   , INTENT(OUT) :: eDn(nbl,kz)
  REAL   , INTENT(OUT) :: eUp(nbl,kz)
  REAL   , INTENT(OUT) :: fDir(nbl,kz)
  REAL   , INTENT(OUT) :: fDn(nbl,kz)
  REAL   , INTENT(OUT) :: fUp(nbl,kz)

  REAL :: dt(nbl,kz), om(nbl,kz), g(nbl,kz)
  REAL,DIMENSION(nbl) :: dtabs, dtsct, dscld, dacld, dsaer, daaer
  INTEGER :: i, iii
  ! specific two ps2str
  REAL,DIMENSION(nbl,kz) :: ediri, edni, eupi
  REAL,DIMENSION(nbl,kz) :: fdiri, fdni, fupi
  LOGICAL,PARAMETER :: delta=.true.
  !  specific to psndo:
  REAL :: pmcld, pmray, pmaer
  REAL :: om1
  INTEGER :: istr, iob
  INTEGER, PARAMETER :: maxcly=151
  INTEGER, PARAMETER :: maxulv=151
  INTEGER, PARAMETER :: maxumu=32
  INTEGER, PARAMETER :: maxcmu=32
  INTEGER, PARAMETER :: maxphi=3
  INTEGER :: nlyr, numu
  REAL :: dtauc(nbl,maxcly ), pmom(nbl,0:maxcmu,maxcly),  &
      ssalb(nbl, maxcly ), umu( maxumu ), cwt( maxumu )
  REAL :: umu0(nbl)
  REAL :: rfldir(nbl,maxulv), rfldn(nbl,maxulv), flup(nbl,maxulv),  &
      u0u(maxumu, maxulv),  &
      uavgso(nbl,maxulv), uavgup(nbl,maxulv), uavgdn(nbl,maxulv),  &
      sindir(maxulv), sinup(maxulv), sindn(maxulv)

  DO iob=1,nbl
     DO  i = 1, nz(iob)
        fdir(iob,i)= 0.0
        fup(iob,i) = 0.0
        fdn(iob,i) = 0.0
        edir(iob,i)= 0.0
        eup(iob,i) = 0.0
        edn(iob,i) = 0.0
     END DO
  END DO

  DO iob=1,nbl
     umu0(iob) = COS(sza(iob)*rpd)
  END DO

  DO iob=1,nbl
     DO  i = 1, nz(iob) - 1
        dscld(iob) = dtCld(iob,i,iw)* omCld(iob,i,iw)
        dacld(iob) = dtCld(iob,i,iw)* (1.- omCld(iob,i,iw))
        dsaer(iob) = dtAer(iob,i,iw)* omAer(iob,i,iw)
        daaer(iob) = dtAer(iob,i,iw)* (1.- omAer(iob,i,iw))
        dtsct(iob) = dtRl(iob,i,iw) + dscld(iob) + dsaer(iob)
        dtabs(iob) = dtso2(iob,i,iw) + &
                dtO2(iob,i,iw) + dtO3(iob,i,iw) + &
                dtno2(iob,i,iw) + dacld(iob) + daaer(iob)
        dtabs(iob) = AMAX1(dtabs(iob),1./largest)
        dtsct(iob) = AMAX1(dtsct(iob),1./largest)
        ! invert z-coordinate:
        iii = nz(iob) - i
        dt(iob,iii) = dtsct(iob) + dtabs(iob)
        om(iob,iii) = dtsct(iob)/(dtsct(iob) + dtabs(iob))
        IF(dtsct(iob) == 1./largest) om(iob,iii) = 1./largest
        g(iob,iii) = (gCld(iob,i,iw)*dscld(iob) +&
              gAer(iob,i,iw)*dsaer(iob))/dtsct(iob)

        IF(nstr < 2) CYCLE

        ! DISORD parameters
        om1 = AMIN1(om(iob,iii),1.-precis)
        ssalb(iob,iii) = AMAX1(om1,precis)
        dtauc(iob,iii) = AMAX1(dt(iob,ii),precis)

        !  phase function - assume Henyey-Greenstein for cloud and aerosol
        !  and Rayleigh for molecular scattering
        pmom(iob,0,iii) = 1.0
        DO  istr = 1, nstr
           pmcld = gCld(iob,i,iw)**(istr)
           pmaer = gAer(iob,i,iw)**(istr)
          IF(istr == 2) THEN
             pmray = 0.1
          ELSE
             pmray = 0.0
          END IF
          pmom(iob,istr,iii) = (pmcld*dscld(iob) + pmaer*dsaer(iob) + &
                   pmray*dtRl(iob,i,iw)) / dtsct(iob)
        END DO
     END DO
  END DO

  ! call rt routine:
  IF( nstr < 2 ) THEN
     DO iob=1,nbl
        CALL ps2str(nz(iob),sza(iob),&
           albedo(iob,iw),dt(iob,:),om(iob,:),g(iob,:), &
	   dsdh(iob,:,:), nid(iob,:),&
	   delta,fdiri(iob,:),fupi(iob,:),fdni(iob,:),ediri(iob,:), &
	   eupi(iob,:), edni(iob,:))
     END DO
  ELSE
     DO iob=1,nbl
        nlyr = nz(iob) - 1
        CALL  psndo(dsdh(iob,:,:), &
           nid(iob,:),&
	   nlyr, dtauc(iob,:), ssalb(iob,:), pmom(iob,:,:),  &
           albedo(iob,iw), nstr, numu, umu, cwt, umu0(iob),  &
           maxcly, maxulv, maxumu, maxcmu, maxphi, rfldir(iob,:), &
	   rfldn(iob,:), flup(iob,:), u0u,  &
           uavgso(iob,:), uavgup(iob,:), uavgdn(iob,:), sindir, sinup, sindn)
     END DO
  END IF

  ! output (invert z-coordinate)

    IF( nstr < 2 ) THEN
       DO iob=1,nbl
          DO  i = 1, nz(iob)
             iii = nz(iob) - i + 1
             fdir(iob,i) = fdiri(iob,iii)
             fup(iob,i)  = fupi(iob,iii)
             fdn(iob,i)  = fdni(iob,iii)
             edir(iob,i) = ediri(iob,iii)
             eup(iob,i)  = eupi(iob,iii)
             edn(iob,i)  = edni(iob,iii)
          END DO
       END DO
    ELSE
       DO iob=1,nbl
          DO  i = 1, nz(iob)
             iii = nz(iob) - i + 1
             edir(iob,i) = rfldir(iob,iii)
             edn(iob,i)  = rfldn(iob,iii)
             eup(iob,i)  = flup(iob,iii)
             fdir(iob,i) = 4.* pi * uavgso(iob,iii)
             fdn(iob,i)  = 4.* pi * uavgdn(iob,iii)
             fup(iob,i)  = 4.* pi * uavgup(iob,iii)
          END DO
       END DO
    END IF


  END SUBROUTINE rtlink

  !=============================================================================*

  SUBROUTINE ps2str(nlevel,zen,rsfc,tauu,omu,gu, dsdh, nid, delta,  &
      fdr, fup, fdn, edr, eup, edn)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Solve two-stream equations for multiple layers.  The subroutine is based =*
  !=  on equations from:  Toon et al., J.Geophys.Res., v94 (D13), Nov 20, 1989.=*
  !=  It contains 9 two-stream methods to choose from.  A pseudo-spherical     =*
  !=  correction has also been added.                                          =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NLEVEL  - INTEGER, number of specified altitude levels in the working (I)=*
  !=            grid                                                           =*
  !=  ZEN     - REAL, solar zenith angle (degrees)                          (I)=*
  !=  RSFC    - REAL, surface albedo at current wavelength                  (I)=*
  !=  TAUU    - REAL, unscaled optical depth of each layer                  (I)=*
  !=  OMU     - REAL, unscaled single scattering albedo of each layer       (I)=*
  !=  GU      - REAL, unscaled asymmetry parameter of each layer            (I)=*
  !=  DSDH    - REAL, slant path of direct beam through each layer crossed  (I)=*
  !=            when travelling from the top of the atmosphere to layer i;     =*
  !=            DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                            =*
  !=  NID     - INTEGER, number of layers crossed by the direct beam when   (I)=*
  !=            travelling from the top of the atmosphere to layer i;          =*
  !=            NID(i), i = 0..NZ-1                                            =*
  !=  DELTA   - LOGICAL, switch to use delta-scaling                        (I)=*
  !=            .TRUE. -> apply delta-scaling                                  =*
  !=            .FALSE.-> do not apply delta-scaling                           =*
  !=  FDR     - REAL, contribution of the direct component to the total     (O)=*
  !=            actinic flux at each altitude level                            =*
  !=  FUP     - REAL, contribution of the diffuse upwelling component to    (O)=*
  !=            the total actinic flux at each altitude level                  =*
  !=  FDN     - REAL, contribution of the diffuse downwelling component to  (O)=*
  !=            the total actinic flux at each altitude level                  =*
  !=  EDR     - REAL, contribution of the direct component to the total     (O)=*
  !=            spectral irradiance at each altitude level                     =*
  !=  EUP     - REAL, contribution of the diffuse upwelling component to    (O)=*
  !=            the total spectral irradiance at each altitude level           =*
  !=  EDN     - REAL, contribution of the diffuse downwelling component to  (O)=*
  !=            the total spectral irradiance at each altitude level           =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nlevel
  INTEGER, INTENT(IN) :: nid(0:kz)
  REAL, INTENT(IN)    :: zen
  REAL, INTENT(IN)    :: rsfc
  REAL, INTENT(IN)    :: tauu(kz)
  REAL, INTENT(IN)    :: omu(kz)
  REAL, INTENT(IN)    :: gu(kz)
  REAL, INTENT(IN)    :: dsdh(0:kz,kz)
  LOGICAL, INTENT(IN) :: delta
  REAL, INTENT(OUT)   :: fdr(kz)
  REAL, INTENT(OUT)   :: fup(kz)
  REAL, INTENT(OUT)   :: fdn(kz)
  REAL, INTENT(OUT)   :: edr(kz)
  REAL, INTENT(OUT)   :: eup(kz)
  REAL, INTENT(OUT)   :: edn(kz)


  INTEGER, PARAMETER :: nrows=2*kz

  REAL :: tausla(0:kz), tauc(0:kz)
  REAL :: mu2(0:kz), mu, sum

  ! internal coefficients and matrix
  REAL :: lam(kz),taun(kz),bgam(kz)
  REAL :: e1(kz),e2(kz),e3(kz),e4(kz)
  REAL :: cup(kz),cdn(kz),cuptn(kz),cdntn(kz)
  REAL :: mu1(kz)
  INTEGER :: row
  REAL :: a(nrows),b(nrows),d(nrows),e(nrows),y(nrows)

  REAL :: pifs, fdn0
  REAL :: gi(kz), omi(kz), tempg
  REAL :: f, g, om
  REAL :: gam1, gam2, gam3, gam4

  ! For calculations of Associated Legendre Polynomials for GAMA1,2,3,4
  ! in delta-function, modified quadrature, hemispheric constant,
  ! Hybrid modified Eddington-delta function metods, p633,Table1.
  ! W.E.Meador and W.R.Weaver, GAS,1980,v37,p.630
  ! W.J.Wiscombe and G.W. Grams, GAS,1976,v33,p2440,
  ! uncomment the following two lines and the appropriate statements further
  ! down.
  !     REAL YLM0, YLM2, YLM4, YLM6, YLM8, YLM10, YLM12, YLMS, BETA0,
  !    >     BETA1, BETAn, amu1, subd

  REAL :: expon, expon0, expon1, divisr, temp, up, dn
  REAL :: ssfc
  INTEGER :: nlayer, mrows, lev

  INTEGER :: i, j

  ! Some additional program constants:


  REAL, PARAMETER :: eps = 1.e-3
  !_______________________________________________________________________

  ! MU = cosine of solar zenith angle
  ! RSFC = surface albedo
  ! TAUU =  unscaled optical depth of each layer
  ! OMU  =  unscaled single scattering albedo
  ! GU   =  unscaled asymmetry factor
  ! KLEV = max dimension of number of layers in atmosphere
  ! NLAYER = number of layers in the atmosphere
  ! NLEVEL = nlayer + 1 = number of levels

  ! initial conditions:  pi*solar flux = 1;  diffuse incidence = 0

  pifs = 1.
  fdn0 = 0.

  nlayer = nlevel - 1

  mu = COS(zen*rpd)

  !************* compute coefficients for each layer:
  ! GAM1 - GAM4 = 2-stream coefficients, different for different approximations
  ! EXPON0 = calculation of e when TAU is zero
  ! EXPON1 = calculation of e when TAU is TAUN
  ! CUP and CDN = calculation when TAU is zero
  ! CUPTN and CDNTN = calc. when TAU is TAUN
  ! DIVISR = prevents division by zero

  DO j = 0, kz
    tauc(j) = 0.
    tausla(j) = 0.
    mu2(j) = 1./SQRT(largest)
  END DO

  IF( .NOT. delta ) THEN
    DO i = 1, nlayer
      gi(i) = gu(i)
      omi(i) = omu(i)
      taun(i) = tauu(i)
    END DO
  ELSE

  ! delta-scaling. Have to be done for delta-Eddington approximation,
  ! delta discrete ordinate, Practical Improved Flux Method, delta function,
  ! and Hybrid modified Eddington-delta function methods approximations

    DO i = 1, nlayer
      f = gu(i)*gu(i)
      gi(i) = (gu(i) - f)/(1 - f)
      omi(i) = (1 - f)*omu(i)/(1 - omu(i)*f)
      taun(i) = (1 - omu(i)*f)*tauu(i)
    END DO
  END IF

  ! calculate slant optical depth at the top of the atmosphere when zen>90.
  ! in this case, higher altitude of the top layer is recommended which can
  ! be easily changed in gridz.f.
  IF(zen > 90.0) THEN
    IF(nid(0) < 0) THEN
      tausla(0) = largest
    ELSE
      sum = 0.0
      DO j = 1, nid(0)
        sum = sum + 2.*taun(j)*dsdh(0,j)
      END DO
      tausla(0) = sum
    END IF
  END IF


  DO   i = 1, nlayer
    g = gi(i)
    om = omi(i)
    tauc(i) = tauc(i-1) + taun(i)
    ! stay away from 1 by precision.  For g, also stay away from -1
    tempg = AMIN1(ABS(g),1. - precis)
    g = SIGN(tempg,g)
    om = AMIN1(om,1.-precis)

    ! calculate slant optical depth
    IF(nid(i) < 0) THEN
      tausla(i) = largest
    ELSE
      sum = 0.0
      DO j = 1, MIN(nid(i),i)
        sum = sum + taun(j)*dsdh(i,j)
      END DO
      DO j = MIN(nid(i),i)+1,nid(i)
        sum = sum + 2.*taun(j)*dsdh(i,j)
      END DO
      tausla(i) = sum
      IF(tausla(i) == tausla(i-1)) THEN
        mu2(i) = SQRT(largest)
      ELSE
        mu2(i) = (tauc(i)-tauc(i-1))/(tausla(i)-tausla(i-1))
        mu2(i) = SIGN( AMAX1(ABS(mu2(i)),1./SQRT(largest)), mu2(i) )
      END IF
    END IF

    !** the following gamma equations are from pg 16,289, Table 1
    !** save mu1 for each approx. for use in converting irradiance to actinic flux
    ! Eddington approximation(Joseph et al., 1976, JAS, 33, 2452):
    gam1 =  (7. - om*(4. + 3.*g))/4.
    gam2 = -(1. - om*(4. - 3.*g))/4.
    gam3 = (2. - 3.*g*mu)/4.
    gam4 = 1. - gam3
    mu1(i) = 0.5

    ! quadrature (Liou, 1973, JAS, 30, 1303-1326; 1974, JAS, 31, 1473-1475):
    !	       gam1 = 1.7320508*(2. - om*(1. + g))/2.
    !	       gam2 = 1.7320508*om*(1. - g)/2.
    !	       gam3 = (1. - 1.7320508*g*mu)/2.
    !	       gam4 = 1. - gam3
    !	       mu1(i) = 1./sqrt(3.)

    ! hemispheric mean (Toon et al., 1089, JGR, 94, 16287):

    !	       gam1 = 2. - om*(1. + g)
    !	       gam2 = om*(1. - g)
    !	       gam3 = (2. - g*mu)/4.
    !	       gam4 = 1. - gam3
    !	       mu1(i) = 0.5

    ! PIFM  (Zdunkovski et al.,1980, Conrib.Atmos.Phys., 53, 147-166):
    !	      GAM1 = 0.25*(8. - OM*(5. + 3.*G))
    !	      GAM2 = 0.75*OM*(1.-G)
    !	      GAM3 = 0.25*(2.-3.*G*MU)
    !	      GAM4 = 1. - GAM3
    !	      mu1(i) = 0.5

    ! delta discrete ordinates  (Schaller, 1979, Contrib.Atmos.Phys, 52, 17-26):
    !	      GAM1 = 0.5*1.7320508*(2. - OM*(1. + G))
    !	      GAM2 = 0.5*1.7320508*OM*(1.-G)
    !	      GAM3 = 0.5*(1.-1.7320508*G*MU)
    !	      GAM4 = 1. - GAM3
    !	      mu1(i) = 1./sqrt(3.)

    ! Calculations of Associated Legendre Polynomials for GAMA1,2,3,4
    ! in delta-function, modified quadrature, hemispheric constant,
    ! Hybrid modified Eddington-delta function metods, p633,Table1.
    ! W.E.Meador and W.R.Weaver, GAS,1980,v37,p.630
    ! W.J.Wiscombe and G.W. Grams, GAS,1976,v33,p2440
    !	   YLM0 = 2.
    !	   YLM2 = -3.*G*MU
    !	   YLM4 = 0.875*G**3*MU*(5.*MU**2-3.)
    !	   YLM6=-0.171875*G**5*MU*(15.-70.*MU**2+63.*MU**4)
    !	  YLM8=+0.073242*G**7*MU*(-35.+315.*MU**2-693.*MU**4
    !	 *+429.*MU**6)
    !	  YLM10=-0.008118*G**9*MU*(315.-4620.*MU**2+18018.*MU**4
    !	 *-25740.*MU**6+12155.*MU**8)
    !	  YLM12=0.003685*G**11*MU*(-693.+15015.*MU**2-90090.*MU**4
    !	 *+218790.*MU**6-230945.*MU**8+88179.*MU**10)
    !	   YLMS=YLM0+YLM2+YLM4+YLM6+YLM8+YLM10+YLM12
    !	   YLMS=0.25*YLMS
    !	   BETA0 = YLMS

    !	      amu1=1./1.7320508
    !	   YLM0 = 2.
    !	   YLM2 = -3.*G*amu1
    !	   YLM4 = 0.875*G**3*amu1*(5.*amu1**2-3.)
    !	   YLM6=-0.171875*G**5*amu1*(15.-70.*amu1**2+63.*amu1**4)
    !	  YLM8=+0.073242*G**7*amu1*(-35.+315.*amu1**2-693.*amu1**4
    !	 *+429.*amu1**6)
    !	  YLM10=-0.008118*G**9*amu1*(315.-4620.*amu1**2+18018.*amu1**4
    !	 *-25740.*amu1**6+12155.*amu1**8)
    !	  YLM12=0.003685*G**11*amu1*(-693.+15015.*amu1**2-90090.*amu1**4
    !	 *+218790.*amu1**6-230945.*amu1**8+88179.*amu1**10)
    !	   YLMS=YLM0+YLM2+YLM4+YLM6+YLM8+YLM10+YLM12
    !	   YLMS=0.25*YLMS
    !	   BETA1 = YLMS

    !	      BETAn = 0.25*(2. - 1.5*G-0.21875*G**3-0.085938*G**5
    !	 *-0.045776*G**7)


    ! Hybrid modified Eddington-delta function(Meador and Weaver,1980,JAS,37,630):
    !	      subd=4.*(1.-G*G*(1.-MU))
    !	      GAM1 = (7.-3.*G*G-OM*(4.+3.*G)+OM*G*G*(4.*BETA0+3.*G))/subd
    !	      GAM2 =-(1.-G*G-OM*(4.-3.*G)-OM*G*G*(4.*BETA0+3.*G-4.))/subd
    !	      GAM3 = BETA0
    !	      GAM4 = 1. - GAM3
    !	      mu1(i) = (1. - g*g*(1.- mu) )/(2. - g*g)

    !****
    ! delta function  (Meador, and Weaver, 1980, JAS, 37, 630):
    !	      GAM1 = (1. - OM*(1. - beta0))/MU
    !	      GAM2 = OM*BETA0/MU
    !	      GAM3 = BETA0
    !	      GAM4 = 1. - GAM3
    !	      mu1(i) = mu
    !****
    ! modified quadrature (Meador, and Weaver, 1980, JAS, 37, 630):
    !	      GAM1 = 1.7320508*(1. - OM*(1. - beta1))
    !	      GAM2 = 1.7320508*OM*beta1
    !	      GAM3 = BETA0
    !	      GAM4 = 1. - GAM3
    !	      mu1(i) = 1./sqrt(3.)

    ! hemispheric constant (Toon et al., 1989, JGR, 94, 16287):
    !	      GAM1 = 2.*(1. - OM*(1. - betan))
    !	      GAM2 = 2.*OM*BETAn
    !	      GAM3 = BETA0
    !	      GAM4 = 1. - GAM3
    !	      mu1(i) = 0.5

    !****

    ! lambda = pg 16,290 equation 21
    ! big gamma = pg 16,290 equation 22
    ! if gam2 = 0., then bgam = 0.

    lam(i) = SQRT(gam1*gam1 - gam2*gam2)

    IF( gam2 /= 0.) THEN
      bgam(i) = (gam1 - lam(i))/gam2
    ELSE
      bgam(i) = 0.
    END IF

    expon = EXP(-lam(i)*taun(i))

    ! e1 - e4 = pg 16,292 equation 44
    e1(i) = 1. + bgam(i)*expon
    e2(i) = 1. - bgam(i)*expon
    e3(i) = bgam(i) + expon
    e4(i) = bgam(i) - expon

    ! the following sets up for the C equations 23, and 24
    ! found on page 16,290
    ! prevent division by zero (if LAMBDA=1/MU, shift 1/MU^2 by EPS = 1.E-3
    ! which is approx equiv to shifting MU by 0.5*EPS* (MU)**3
    expon0 = EXP(-tausla(i-1))
    expon1 = EXP(-tausla(i))

    divisr = lam(i)*lam(i) - 1./(mu2(i)*mu2(i))
    temp = AMAX1(eps,ABS(divisr))
    divisr = SIGN(temp,divisr)

    up = om*pifs*((gam1 - 1./mu2(i))*gam3 + gam4*gam2)/divisr
    dn = om*pifs*((gam1 + 1./mu2(i))*gam4 + gam2*gam3)/divisr

    ! cup and cdn are when tau is equal to zero
    ! cuptn and cdntn are when tau is equal to taun
    cup(i) = up*expon0
    cdn(i) = dn*expon0
    cuptn(i) = up*expon1
    cdntn(i) = dn*expon1

  END DO

  !**************** set up matrix ******
  ! ssfc = pg 16,292 equation 37  where pi Fs is one (unity).

  ssfc = rsfc*mu*EXP(-tausla(nlayer))*pifs

  ! MROWS = the number of rows in the matrix

  mrows = 2*nlayer

  ! the following are from pg 16,292  equations 39 - 43.
  ! set up first row of matrix:

  i = 1
  a(1) = 0.
  b(1) = e1(i)
  d(1) = -e2(i)
  e(1) = fdn0 - cdn(i)

  row=1

  ! set up odd rows 3 thru (MROWS - 1):

  i = 0
  DO   row = 3, mrows - 1, 2
    i = i + 1
    a(row) = e2(i)*e3(i) - e4(i)*e1(i)
    b(row) = e1(i)*e1(i + 1) - e3(i)*e3(i + 1)
    d(row) = e3(i)*e4(i + 1) - e1(i)*e2(i + 1)
    e(row) = e3(i)*(cup(i + 1) - cuptn(i)) + e1(i)*(cdntn(i) - cdn(i + 1))
  END DO

  ! set up even rows 2 thru (MROWS - 2):

  i = 0
  DO   row = 2, mrows - 2, 2
    i = i + 1
    a(row) = e2(i + 1)*e1(i) - e3(i)*e4(i + 1)
    b(row) = e2(i)*e2(i + 1) - e4(i)*e4(i + 1)
    d(row) = e1(i + 1)*e4(i + 1) - e2(i + 1)*e3(i + 1)
    e(row) = (cup(i + 1) - cuptn(i))*e2(i + 1) -  &
        (cdn(i + 1) - cdntn(i))*e4(i + 1)
  END DO

  ! set up last row of matrix at MROWS:

  row = mrows
  i = nlayer

  a(row) = e1(i) - rsfc*e3(i)
  b(row) = e2(i) - rsfc*e4(i)
  d(row) = 0.
  e(row) = ssfc - cuptn(i) + rsfc*cdntn(i)

  ! solve tri-diagonal matrix:
  CALL tridag(a, b, d, e, y, mrows)

  !*** unfold solution of matrix, compute output fluxes:
  row = 1
  lev = 1
  j = 1

  ! the following equations are from pg 16,291  equations 31 & 32
  fdr(lev) = EXP( -tausla(0) )
  edr(lev) = mu * fdr(lev)
  edn(lev) = fdn0
  eup(lev) =  y(row)*e3(j) - y(row + 1)*e4(j) + cup(j)
  fdn(lev) = edn(lev)/mu1(lev)
  fup(lev) = eup(lev)/mu1(lev)

  DO   lev = 2, nlayer + 1
    fdr(lev) = EXP(-tausla(lev-1))
    edr(lev) =  mu *fdr(lev)
    edn(lev) =  y(row)*e3(j) + y(row + 1)*e4(j) + cdntn(j)
    eup(lev) =  y(row)*e1(j) + y(row + 1)*e2(j) + cuptn(j)
    !LFR
    IF(fdr(lev)<0) fdr(lev)=0
    IF(edr(lev)<0) edr(lev)=0
    IF(edn(lev)<0) edn(lev)=0
    IF(eup(lev)<0) eup(lev)=0
    !LFR
    fdn(lev) = edn(lev)/mu1(j)
    fup(lev) = eup(lev)/mu1(j)
    row = row + 2
    j = j + 1
  END DO
  !Debug-test LFR
  fup(nlayer+1)=0.0
  fdn(nlayer+1)=0.0
  !End of debug
  !_______________________________________________________________________

  END SUBROUTINE ps2str

  !=============================================================================*

  SUBROUTINE tridag(a,b,c,r,u,n)
  !_______________________________________________________________________
  ! solves tridiagonal system.  From Numerical Recipies, p. 40
  !_______________________________________________________________________
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL, INTENT(IN)    :: a(n)
  REAL, INTENT(IN)    :: b(n)
  REAL, INTENT(IN)    :: c(n)
  REAL, INTENT(IN)    :: r(n)
  REAL, INTENT(OUT)   :: u(n)

  INTEGER :: j

  REAL :: bet
  REAL,DIMENSION(2*kz) :: gam
  !_______________________________________________________________________

  IF (b(1) == 0.) STOP 1001
  bet   = b(1)
  u(1) = r(1)/bet
  DO   j = 2, n
    gam(j) = c(j - 1)/bet
    bet = b(j) - a(j)*gam(j)
    IF (bet == 0.) STOP 2002
    u(j) = (r(j) - a(j)*u(j - 1))/bet
  END DO
  DO   j = n - 1, 1, -1
    u(j) = u(j) - gam(j + 1)*u(j + 1)
  END DO
  !_______________________________________________________________________

  END SUBROUTINE tridag

  SUBROUTINE psndo( dsdh, nid, nlyr, dtauc, ssalb, pmom,  &
      albedo, nstr, numu, umu, cwt, umu0,  &
      maxcly, maxulv, maxumu, maxcmu, maxphi, rfldir, rfldn, flup, u0u,  &
      uavgso, uavgup, uavgdn, sindir, sinup, sindn )
     IMPLICIT NONE
  ! Improved handling of numerical instabilities. Bernhard Mayer on 5/3/99.
  !  disort seems to produce unstable results for certain combinations
  !  of single scattering albedo and phase function. A temporary fix has been
  !  introduced to avoid this problem: The original instability check in
  !  UPBEAM fails on certain compiler/machine combinations (e.g., gcc/LINUX,
  !  or xlf/IBM RS6000). This check has therefore been replaced by a new one.
  !  Whenever UPBEAM reports an instability, the single scattering albedo
  !  of the respective layer is changed by a small amount, and the
  !  calculation is repeated until numerically stable conditions are reached
  !  (all the necessary changes are confined to the new subroutine SOLVEC
  !  and the slighly changed subroutine UPBEAM). To check for potential
  !  instabilities, the variable 'RCOND' returned by SGECO is compared to
  !  a machine-dependent constant, 'MINRCOND'. The value of this constant
  !  determines (a) if really all instabilities are caught; and (b) the
  !  amount by which the single scattering albedo has to be changed. The
  !  value of 'MINRCOND' is therefore a compromise between numerical
  !  stability on the one hand and uncertainties introduced by changing
  !  the atmospheric conditions and increased computational time on the
  !  other hand (an increase of MINRCOND will lead to the detection of
  !  more potential numerical instabilities, and thus to an increase in
  !  computational time; by changing the atmospheric conditions, that is,
  !  the single scattering albedo, the result might however be changed
  !  unfavourably, if the change is too large). From a limited number
  !  of experiments we found that 'MINRCOND = 5000. * R1MACH(4)' seems
  !  to be a good choice if high accuracy is required (more tests are
  !  definitely neccessary!). If an instability is encountered, a message
  !  is printed telling about neccessary changes to the single scattering
  !  albedo. This message may be switched off by setting 'DEBUG = .FALSE.'
  !  in subroutine SOLVEC.


  ! modified to calculate sine-weighted intensities. Bernhard Mayer on 2/12/99.
  ! modified to handle some numerical instabilities. Chris Fischer on 1/22/99.
  ! modified by adding pseudo-spherical correction. Jun Zeng on 3/11/97.
  ! dsdh: slant path of direct beam through each layer crossed
  !       when travelling from the top of the atmosphere to layer i;
  !       dsdh(i,j), i = 0..nlyr, j = 1..nlyr;
  ! nid:  number of layers crossed by the direct beam when
  !       travelling from the top of the atmosphere to layer i;
  !       NID(i), i = 0..nlyr.
  ! uavgso, uvagup, and uvagdn are direct, downward diffuse, and upward
  ! diffuse actinic flux (mean intensity).
  ! u0u is the azimuthally averaged intensity, check DISORT.doc for details.
  ! *******************************************************************
  !       Plane-parallel discrete ordinates radiative transfer program
  !                      V E R S I O N    1.1
  !             ( see DISORT.DOC for complete documentation )
  ! *******************************************************************


  ! +------------------------------------------------------------------+
  !  Calling Tree (omitting calls to ERRMSG):
  !  (routines in parentheses are not in this file)

  !  DISORT-+-(R1MACH)
  !         +-ZEROIT
  !         +-CHEKIN-+-(WRTBAD)
  !         |        +-(WRTDIM)
  !         |        +-DREF
  !         +-ZEROAL
  !         +-SETDIS-+-QGAUSN (1)-+-(D1MACH)
  !         +-PRTINP
  !         +-LEPOLY see 2
  !         +-SURFAC-+-QGAUSN see 1
  !         |        +-LEPOLY see 2
  !         |        +-ZEROIT
  !         +-SOLEIG see 3
  !         +-UPBEAM-+-(SGECO)
  !         |        +-(SGESL)
  !         +-TERPEV
  !         +-TERPSO
  !         +-SETMTX see 4
  !         +-SOLVE0-+-ZEROIT
  !         |        +-(SGBCO)
  !         |        +-(SGBSL)
  !         +-FLUXES--ZEROIT
  !         +-PRAVIN
  !         +-RATIO--(R1MACH)
  !         +-PRTINT

  ! *** Intrinsic Functions used in DISORT package which take
  !     non-negligible amount of time:

  !    EXP :  Called by- ALBTRN, ALTRIN, CMPINT, FLUXES, SETDIS,
  !                      SETMTX, SPALTR, USRINT, PLKAVG

  !    SQRT : Called by- ASYMTX, LEPOLY, SOLEIG

  ! +-------------------------------------------------------------------+

  !  Index conventions (for all DO-loops and all variable descriptions):

  !     IU     :  for user polar angles

  !  IQ,JQ,KQ  :  for computational polar angles ('quadrature angles')

  !   IQ/2     :  for half the computational polar angles (just the ones
  !               in either 0-90 degrees, or 90-180 degrees)

  !     J      :  for user azimuthal angles

  !     K,L    :  for Legendre expansion coefficients or, alternatively,
  !               subscripts of associated Legendre polynomials

  !     LU     :  for user levels

  !     LC     :  for computational layers (each having a different
  !               single-scatter albedo and/or phase function)

  !    LEV     :  for computational levels

  !    MAZIM   :  for azimuthal components in Fourier cosine expansion
  !               of intensity and phase function

  ! +------------------------------------------------------------------+

  !               I N T E R N A L    V A R I A B L E S

  !   AMB(IQ/2,IQ/2)    First matrix factor in reduced eigenvalue problem
  !                     of Eqs. SS(12), STWJ(8E)  (used only in SOLEIG)

  !   APB(IQ/2,IQ/2)    Second matrix factor in reduced eigenvalue problem
  !                     of Eqs. SS(12), STWJ(8E)  (used only in SOLEIG)

  !   ARRAY(IQ,IQ)      Scratch matrix for SOLEIG, UPBEAM and UPISOT
  !                     (see each subroutine for definition)

  !   B()               Right-hand side vector of Eq. SC(5) going into
  !                     SOLVE0,1;  returns as solution vector
  !                     vector  L, the constants of integration

  !   BDR(IQ/2,0:IQ/2)  Bottom-boundary bidirectional reflectivity for a
  !                     given azimuthal component.  First index always
  !                     refers to a computational angle.  Second index:
  !                     if zero, refers to incident beam angle UMU0;
  !                     if non-zero, refers to a computational angle.

  !   BEM(IQ/2)         Bottom-boundary directional emissivity at compu-
  !                     tational angles.

  !   BPLANK            Intensity emitted from bottom boundary

  !   CBAND()           Matrix of left-hand side of the linear system
  !                     Eq. SC(5), scaled by Eq. SC(12);  in banded
  !                     form required by LINPACK solution routines

  !   CC(IQ,IQ)         C-sub-IJ in Eq. SS(5)

  !   CMU(IQ)           Computational polar angles (Gaussian)

  !   CWT(IQ)           Quadrature weights corresponding to CMU

  !   DELM0             Kronecker delta, delta-sub-M0, where M = MAZIM
  !                     is the number of the Fourier component in the
  !                     azimuth cosine expansion

  !   DITHER            Small quantity subtracted from single-scattering
  !                     albedos of unity, in order to avoid using special
  !                     case formulas;  prevents an eigenvalue of exactly
  !                     zero from occurring, which would cause an
  !                     immediate overflow

  !   DTAUCP(LC)        Computational-layer optical depths (delta-M-scaled
  !                     if DELTAM = TRUE, otherwise equal to DTAUC)

  !   EMU(IU)           Bottom-boundary directional emissivity at user
  !                     angles.

  !   EVAL(IQ)          Temporary storage for eigenvalues of Eq. SS(12)

  !   EVECC(IQ,IQ)      Complete eigenvectors of SS(7) on return from
  !                     SOLEIG; stored permanently in  GC

  !   EXPBEA(LC)        Transmission of direct beam in delta-M optical
  !                     depth coordinates

  !   FLYR(LC)          Truncated fraction in delta-M method

  !   GL(K,LC)          Phase function Legendre polynomial expansion
  !                     coefficients, calculated from PMOM by
  !                     including single-scattering albedo, factor
  !                     2K+1, and (if DELTAM=TRUE) the delta-M
  !                     scaling

  !   GC(IQ,IQ,LC)      Eigenvectors at polar quadrature angles,
  !                     g  in Eq. SC(1)

  !   GU(IU,IQ,LC)      Eigenvectors interpolated to user polar angles
  !                     ( g  in Eqs. SC(3) and S1(8-9), i.e.
  !                       G without the L factor )

  !   HLPR()            Legendre coefficients of bottom bidirectional
  !                     reflectivity (after inclusion of 2K+1 factor)

  !   IPVT(LC*IQ)       Integer vector of pivot indices for LINPACK
  !                     routines

  !   KK(IQ,LC)         Eigenvalues of coeff. matrix in Eq. SS(7)

  !   KCONV             Counter in azimuth convergence test

  !   LAYRU(LU)         Computational layer in which user output level
  !                     UTAU(LU) is located

  !   LL(IQ,LC)         Constants of integration L in Eq. SC(1),
  !                     obtained by solving scaled version of Eq. SC(5)

  !   LYRCUT            TRUE, radiation is assumed zero below layer
  !                     NCUT because of almost complete absorption

  !   NAZ               Number of azimuthal components considered

  !   NCUT              Computational layer number in which absorption
  !                     optical depth first exceeds ABSCUT

  !   OPRIM(LC)         Single scattering albedo after delta-M scaling

  !   PASS1             TRUE on first entry, FALSE thereafter

  !   PKAG(0:LC)        Integrated Planck function for internal emission

  !   PSI(IQ)           Sum just after square bracket in  Eq. SD(9)

  !   RMU(IU,0:IQ)      Bottom-boundary bidirectional reflectivity for a
  !                     given azimuthal component.  First index always
  !                     refers to a user angle.  Second index:
  !                     if zero, refers to incident beam angle UMU0;
  !                     if non-zero, refers to a computational angle.

  !   TAUC(0:LC)        Cumulative optical depth (un-delta-M-scaled)

  !   TAUCPR(0:LC)      Cumulative optical depth (delta-M-scaled if
  !                     DELTAM = TRUE, otherwise equal to TAUC)

  !   TPLANK            Intensity emitted from top boundary

  !   UUM(IU,LU)        Expansion coefficients when the intensity
  !                     (u-super-M) is expanded in Fourier cosine series
  !                     in azimuth angle

  !   U0C(IQ,LU)        Azimuthally-averaged intensity

  !   UTAUPR(LU)        Optical depths of user output levels in delta-M
  !                     coordinates;  equal to  UTAU(LU) if no delta-M

  !   WK()              scratch array

  !   XR0(LC)           X-sub-zero in expansion of thermal source func-
  !                     tion preceding Eq. SS(14) (has no mu-dependence)

  !   XR1(LC)           X-sub-one in expansion of thermal source func-
  !                     tion;  see  Eqs. SS(14-16)

  !   YLM0(L)           Normalized associated Legendre polynomial
  !                     of subscript L at the beam angle (not saved
  !                     as function of superscipt M)

  !   YLMC(L,IQ)        Normalized associated Legendre polynomial
  !                     of subscript L at the computational angles
  !                     (not saved as function of superscipt M)

  !   YLMU(L,IU)        Normalized associated Legendre polynomial
  !                     of subscript L at the user angles
  !                     (not saved as function of superscipt M)

  !   Z()               scratch array used in  SOLVE0,1  to solve a
  !                     linear system for the constants of integration

  !   Z0(IQ)            Solution vectors Z-sub-zero of Eq. SS(16)

  !   Z0U(IU,LC)        Z-sub-zero in Eq. SS(16) interpolated to user
  !                     angles from an equation derived from SS(16)

  !   Z1(IQ)            Solution vectors Z-sub-one  of Eq. SS(16)

  !   Z1U(IU,LC)        Z-sub-one in Eq. SS(16) interpolated to user
  !                     angles from an equation derived from SS(16)

  !   ZBEAM(IU,LC)      Particular solution for beam source

  !   ZJ(IQ)            Right-hand side vector  X-sub-zero in
  !                     Eq. SS(19), also the solution vector
  !                     Z-sub-zero after solving that system

  !   ZZ(IQ,LC)         Permanent storage for the beam source vectors ZJ

  !   ZPLK0(IQ,LC)      Permanent storage for the thermal source
  !                     vectors  Z0  obtained by solving  Eq. SS(16)

  !   ZPLK1(IQ,LC)      Permanent storage for the thermal source
  !                     vectors  Z1  obtained by solving  Eq. SS(16)

  ! +-------------------------------------------------------------------+

  !  LOCAL SYMBOLIC DIMENSIONS (have big effect on storage requirements):

  !       MXCLY  = Max no. of computational layers
  !       MXULV  = Max no. of output levels
  !       MXCMU  = Max no. of computation polar angles
  !       MXUMU  = Max no. of output polar angles
  !       MXPHI  = Max no. of output azimuthal angles

  ! +-------------------------------------------------------------------+

  INTEGER, INTENT(IN)   :: maxcly
  INTEGER, INTENT(IN)   :: maxulv
  INTEGER, INTENT(IN)   :: maxumu
  INTEGER, INTENT(IN)   :: maxcmu
  INTEGER, INTENT(IN)   :: maxphi
  INTEGER, INTENT(IN)   :: nstr
  INTEGER, INTENT(INOUT)   :: numu
  INTEGER, INTENT(IN)   :: nlyr
  INTEGER, INTENT(IN)   :: nid(0:kz)
  REAL, INTENT(IN)      :: dsdh(0:kz,kz)
  REAL, INTENT(IN)      :: dtauc( maxcly )
  REAL, INTENT(INOUT)      :: pmom( 0:maxcmu, maxcly )
  REAL, INTENT(IN)      :: albedo
  REAL, INTENT(INOUT)      :: umu( maxumu )
  REAL, INTENT(INOUT)      :: cwt( maxcmu )
  REAL, INTENT(IN)      :: umu0
  REAL, INTENT(INOUT)   :: ssalb( maxcly )
  REAL, INTENT(OUT)     :: rfldir( maxulv )
  REAL, INTENT(OUT)     :: rfldn( maxulv )
  REAL, INTENT(OUT)     :: flup( maxulv )
  REAL, INTENT(OUT)     :: u0u( maxumu, maxulv )
  REAL, INTENT(OUT)     :: uavgso( maxulv )
  REAL, INTENT(OUT)     :: uavgup( maxulv )
  REAL, INTENT(OUT)     :: uavgdn( maxulv )
  REAL, INTENT(OUT)     :: sindir( maxulv )
  REAL, INTENT(OUT)     :: sinup( maxulv )
  REAL, INTENT(OUT)     :: sindn( maxulv )

  INTEGER, PARAMETER :: mxcly = 151
  INTEGER, PARAMETER :: mxulv = 151
  INTEGER, PARAMETER :: mxcmu = 32
  INTEGER, PARAMETER :: mxumu = 32
  INTEGER, PARAMETER :: mxphi = 3
  INTEGER, PARAMETER :: mi = mxcmu / 2
  INTEGER, PARAMETER :: mi9m2 = 9*mi - 2
  INTEGER, PARAMETER :: nnlyri = mxcmu*mxcly
  !     ..
  !     .. Scalar Arguments ..

  INTEGER :: ntau
  REAL :: btemp, temis, ttemp,  wvnmhi, wvnmlo

  !     sherical geometry


  REAL :: tausla(0:kz), tauslau(0:kz), mu2(0:kz)
  !     ..
  !     .. Array Arguments ..

  REAL :: albmed( maxumu ), dfdt( maxulv ),  hl( 0:maxcmu ), phi( maxphi ),  &
        temper( 0:maxcly ),  &
      trnmed( maxumu ), uavg( maxulv ),  utau( maxulv ),  &
      uu( maxumu, maxulv, maxphi )

  !     ..
  !     .. Local Scalars ..

  LOGICAL :: lyrcut
  INTEGER :: iq, iu, j, kconv, l, lc, lu, mazim, naz, ncol, ncos, ncut, nn
  REAL :: azerr, azterm, bplank, cosphi, delm0,  &
      sgn, tplank
  !     ..
  !     .. Local Arrays ..
  REAL :: angcos(maxumu)
  INTEGER :: ipvt( nnlyri ), layru( mxulv )

  REAL :: amb( mi, mi ), apb( mi, mi ), array( mxcmu, mxcmu ),  &
      b( nnlyri ), bdr( mi, 0:mi ), bem( mi ),  &
      cband( mi9m2, nnlyri ), cc( mxcmu, mxcmu ), cmu( mxcmu ), dtaucp( mxcly ),  &
      emu( mxumu ), eval( mi ), evecc( mxcmu, mxcmu ),  &
      expbea( 0:mxcly ), fldir( mxulv ), fldn( mxulv ),  &
      flyr( mxcly ), gc( mxcmu, mxcmu, mxcly ),  &
      gl( 0:mxcmu, mxcly ), gu( mxumu, mxcmu, mxcly ),  &
      hlpr( 0:mxcmu ), kk( mxcmu, mxcly ), ll( mxcmu, mxcly ),  &
      oprim( mxcly ), phirad( mxphi ), pkag( 0:mxcly ),  &
      psi( mxcmu ), rmu( mxumu, 0:mi ), tauc( 0:mxcly ),  &
      taucpr( 0:mxcly ), u0c( mxcmu, mxulv ), utaupr( mxulv ),  &
      uum( mxumu, mxulv ), wk( mxcmu ), xr0( mxcly ),  &
      xr1( mxcly ), ylm0( 0:mxcmu ), ylmc( 0:mxcmu, mxcmu ),  &
      ylmu( 0:mxcmu, mxumu ), z( nnlyri ), z0( mxcmu ),  &
      z0u( mxumu, mxcly ), z1( mxcmu ), z1u( mxumu, mxcly ),  &
      zbeam( mxumu, mxcly ), zj( mxcmu ),  &
      zplk0( mxcmu, mxcly ), zplk1( mxcmu, mxcly ), zz( mxcmu, mxcly )

  !gy added glsave and dgl to allow adjustable dimensioning in SOLVEC
  REAL :: glsave( 0:mxcmu ), dgl( 0:mxcmu )

  DOUBLE PRECISION :: aad( mi, mi ), evald( mi ), eveccd( mi, mi ), wkd( mxcmu )
  !     ..
  !     .. External Functions ..

!  REAL :: plkavg!,  ratio
!LFR>   EXTERNAL  plkavg,  ratio
  !     ..
  !     .. External Subroutines ..

!LFR>   EXTERNAL  chekin, fluxes, lepoly, pravin, prtinp,  &
!LFR>       prtint, setdis, setmtx, soleig, solve0, surfac, upbeam, zeroal, zeroit
  !     ..
  !     .. Intrinsic Functions ..

  INTRINSIC ABS, ASIN, COS, LEN, MAX

  ! Discrete ordinate constants:
  ! For pseudo-spherical DISORT, PLANK, USRTAU and USRANG must be .FALSE.;
  ! ONLYFL must be .TRUE.; FBEAM = 1.; FISOT = 0.; IBCND = 0
  LOGICAL, PARAMETER :: lamber=.true.
  LOGICAL, PARAMETER :: usrtau=.false.
  LOGICAL, PARAMETER :: plank=.false.
  LOGICAL, PARAMETER :: usrang=.false.
  LOGICAL, PARAMETER :: onlyfl=.true.
  LOGICAL, PARAMETER :: deltam=.true. ! delat-M scaling option
  LOGICAL, DIMENSION(7),PARAMETER :: prnt=(/.false.,.false.,.false.,.false., &
                                            .false.,.false.,.false./)
  REAL,PARAMETER :: accur=0.0001
  CHARACTER(LEN=127),PARAMETER :: header=REPEAT(' ',127)
  INTEGER,PARAMETER :: nphi=0
  INTEGER,PARAMETER :: ibcnd=0
  REAL,PARAMETER :: fbeam=1.0
  REAL,PARAMETER :: fisot=0.0
  REAL,PARAMETER :: phi0=0.0

  10 CONTINUE

  !** Calculate cumulative optical depth
  !   and dither single-scatter albedo
  !   to improve numerical behavior of
  !   eigenvalue/vector computation
  CALL zeroit( tauc, mxcly + 1 )

  DO  lc = 1, nlyr

    IF( ssalb( lc ) == 1.0 ) ssalb( lc ) = 1.0 - dither
    tauc( lc ) = tauc( lc - 1 ) + dtauc( lc )

  END DO
  !                                ** Check input dimensions and variables

  CALL chekin( nlyr, dtauc, ssalb, pmom, temper, wvnmlo, wvnmhi,  &
      usrtau, ntau, utau, nstr, usrang, numu, umu, nphi,  &
      phi, ibcnd, fbeam, umu0, phi0, fisot, lamber, albedo,  &
      hl, btemp, ttemp, temis, plank, onlyfl, accur, tauc,  &
      maxcly, maxulv, maxumu, maxcmu, maxphi, mxcly, mxulv, mxumu, mxcmu, mxphi )

  !                                 ** Zero internal and output arrays

  CALL  zeroal( mxcly, expbea(1), flyr, oprim, taucpr(1), xr0, xr1,  &
      mxcmu, cmu, cwt, psi, wk, z0, z1, zj, mxcmu+1, hlpr, ylm0,  &
      mxcmu**2, array, cc, evecc, (mxcmu+1)*mxcly, gl,  &
      (mxcmu+1)*mxcmu, ylmc, (mxcmu+1)*mxumu, ylmu,  &
      mxcmu*mxcly, kk, ll, zz, zplk0, zplk1, mxcmu**2*mxcly, gc,  &
      mxulv, layru, utaupr, mxumu*mxcmu*mxcly, gu,  &
      mxumu*mxcly, z0u, z1u, zbeam, mi, eval,  &
      mi**2, amb, apb, nnlyri, ipvt, z,  &
      maxulv, rfldir, rfldn, flup, uavg, dfdt, maxumu, albmed, trnmed,  &
      maxumu*maxulv, u0u, maxumu*maxulv*maxphi, uu )

  !                                 ** Perform various setup operations

  CALL setdis( dsdh, nid, tausla, tauslau, mu2,  &
      cmu, cwt, deltam, dtauc, dtaucp, expbea, flyr,  &
      gl, hl, hlpr, ibcnd, lamber, layru, lyrcut, maxumu,  &
      maxcmu, mxcmu, ncut, nlyr, ntau, nn, nstr, plank,  &
      numu, onlyfl, oprim, pmom, ssalb, tauc, taucpr, utau,  &
      utaupr, umu, umu0, usrtau, usrang )

  !                                 ** Print input information
  IF ( prnt(1) ) CALL prtinp( nlyr, dtauc, dtaucp, ssalb, pmom, temper,  &
      wvnmlo, wvnmhi, ntau, utau, nstr, numu, umu,  &
      nphi, phi, ibcnd, fbeam, umu0, phi0, fisot,  &
      lamber, albedo, hl, btemp, ttemp, temis,  &
      deltam, plank, onlyfl, accur, flyr, lyrcut,  &
      oprim, tauc, taucpr, maxcmu, prnt(7) )

  !                              ** Handle special case for getting albedo
  !                                 and transmissivity of medium for many
  !                                 beam angles at once
  !                                   ** Calculate Planck functions

  bplank = 0.0
  tplank = 0.0
  CALL zeroit( pkag, mxcly + 1 )

  ! ========  BEGIN LOOP TO SUM AZIMUTHAL COMPONENTS OF INTENSITY  =======
  !           (EQ STWJ 5)

  kconv  = 0
  naz  = nstr - 1
  !                                    ** Azimuth-independent case

  IF( fbeam == 0.0 .OR. ( 1.- umu0 ) < 1.e-5 .OR. onlyfl .OR.  &
      ( numu == 1 .AND. ( 1.- umu(1) ) < 1.e-5 ) ) naz = 0

  DO  mazim = 0, naz

    IF( mazim == 0 ) delm0  = 1.0
    IF( mazim > 0 ) delm0  = 0.0

  !                             ** Get normalized associated Legendre
  !                                polynomials for
  !                                (a) incident beam angle cosine
  !                                (b) computational and user polar angle
  !                                    cosines
    IF( fbeam > 0.0 ) THEN

      ncos   = 1
      angcos = -umu0
      !Saulo: precisamos verificar isso aqui abaixo........
      CALL lepoly( ncos, mazim, mxcmu, nstr-1, angcos, ylm0 )

    END IF


    IF( .NOT.onlyfl .AND. usrang )  &
        CALL lepoly( numu, mazim, mxcmu, nstr-1, umu, ylmu )

    CALL lepoly( nn, mazim, mxcmu, nstr-1, cmu, ylmc )

  !                       ** Get normalized associated Legendre polys.
  !                          with negative arguments from those with
  !                          positive arguments; Dave/Armstrong Eq. (15)
    sgn  = - 1.0

    DO  l = mazim, nstr - 1

      sgn  = - sgn

      DO  iq = nn + 1, nstr
        ylmc( l, iq ) = sgn*ylmc( l, iq - nn )
      END DO

    END DO
  !                                 ** Specify users bottom reflectivity
  !                                    and emissivity properties
    IF ( .NOT.lyrcut ) CALL  surfac( albedo, delm0, fbeam, hlpr, lamber,  &
        mi, mazim, mxcmu, mxumu, nn, numu, nstr, onlyfl,  &
        umu, usrang, ylm0, ylmc, ylmu, bdr, emu, bem, rmu )


  ! ===================  BEGIN LOOP ON COMPUTATIONAL LAYERS  =============

    DO  lc = 1, ncut

      CALL solvec( amb, apb, array, cmu, cwt, gl( 0,lc ), mi,  &
          mazim, mxcmu, nn, nstr, ylm0, ylmc, cc,  &
          evecc, eval, kk( 1,lc ), gc( 1,1,lc ), aad, eveccd,  &
          evald, wk, wkd, delm0, fbeam, ipvt, pi,  &
          zj, zz(1,lc), oprim(lc), lc, dither, mu2(lc), glsave, dgl)
  !gy added glsave and dgl to call to allow adjustable dimensioning

    END DO


  ! ===================  END LOOP ON COMPUTATIONAL LAYERS  ===============


  !                      ** Set coefficient matrix of equations combining
  !                         boundary and layer interface conditions

    CALL setmtx( bdr, cband, cmu, cwt, delm0, dtaucp, gc, kk,  &
        lamber, lyrcut, mi, mi9m2, mxcmu, ncol, ncut, nnlyri, nn, nstr, taucpr, wk )

  !                      ** Solve for constants of integration in homo-
  !                         geneous solution (general boundary conditions)

    CALL solve0( b, bdr, bem, bplank, cband, cmu, cwt, expbea,  &
        fbeam, fisot, ipvt, lamber, ll, lyrcut, mazim, mi,  &
        mi9m2, mxcmu, ncol, ncut, nn, nstr, nnlyri, pi,  &
        tplank, taucpr, umu0, z, zz, zplk0, zplk1 )

  !                                  ** Compute upward and downward fluxes

    IF ( mazim == 0 ) CALL fluxes( tausla, tauslau,  &
        cmu, cwt, fbeam, gc, kk, layru, ll, lyrcut,  &
        maxulv, mxcmu, mxulv, ncut, nn, nstr, ntau,  &
        pi, prnt, ssalb, taucpr, umu0, utau, utaupr,  &
        xr0, xr1, zz, zplk0, zplk1, dfdt, flup,  &
        fldn, fldir, rfldir, rfldn, uavg, u0c, uavgso, uavgup, uavgdn,  &
        sindir, sinup, sindn)

    IF( onlyfl ) THEN

      IF( maxumu >= nstr ) THEN
  !                                     ** Save azimuthal-avg intensities
  !                                        at quadrature angles
        DO  lu = 1, ntau

          DO  iq = 1, nstr
            u0u( iq, lu ) = u0c( iq, lu )
          END DO

        END DO

      END IF

      GO TO  170

    END IF


    CALL zeroit( uum, mxumu*mxulv )

    IF( mazim == 0 ) THEN
  !                               ** Save azimuthally averaged intensities

      DO  lu = 1, ntau

        DO  iu = 1, numu
          u0u( iu, lu ) = uum( iu, lu )

          DO  j = 1, nphi
            uu( iu, lu, j ) = uum( iu, lu )
          END DO

        END DO

      END DO
  !                              ** Print azimuthally averaged intensities
  !                                 at user angles

      IF( prnt( 4 ) ) CALL pravin( umu, numu, maxumu, utau, ntau, u0u )
      IF( naz > 0 ) THEN

        CALL zeroit( phirad, mxphi )
        DO  j = 1, nphi
          phirad( j ) = rpd*( phi( j ) - phi0 )
        END DO

      END IF


    ELSE
  !                                ** Increment intensity by current
  !                                   azimuthal component (Fourier
  !                                   cosine series);  Eq SD(2)
      azerr  = 0.0

      DO  j = 1, nphi

        cosphi = COS( mazim*phirad( j ) )

        DO  lu = 1, ntau

          DO  iu = 1, numu
            azterm = uum( iu, lu )*cosphi
            uu( iu, lu, j ) = uu( iu, lu, j ) + azterm
            azerr = MAX( azerr, ratio( ABS(azterm), ABS(uu(iu,lu,j)) ) )
          END DO

        END DO

      END DO

      IF( azerr <= accur ) kconv  = kconv + 1

      IF( kconv >= 2 ) GO TO  170

    END IF

  END DO

  ! ===================  END LOOP ON AZIMUTHAL COMPONENTS  ===============


  !                                          ** Print intensities
  170 CONTINUE
  IF( prnt( 5 ) .AND. .NOT.onlyfl ) CALL prtint( uu, utau, ntau,  &
      umu, numu, phi, nphi, maxulv, maxumu )

  END SUBROUTINE psndo

  SUBROUTINE asymtx( aa, evec, eval, m, ia, ievec, ier, wkd, aad, evecd, evald )

  !    =======  D O U B L E    P R E C I S I O N    V E R S I O N  ======

  !       Solves eigenfunction problem for real asymmetric matrix
  !       for which it is known a priori that the eigenvalues are real.

  !       This is an adaptation of a subroutine EIGRF in the IMSL
  !       library to use real instead of complex arithmetic, accounting
  !       for the known fact that the eigenvalues and eigenvectors in
  !       the discrete ordinate solution are real.  Other changes include
  !       putting all the called subroutines in-line, deleting the
  !       performance index calculation, updating many DO-loops
  !       to Fortran77, and in calculating the machine precision
  !       TOL instead of specifying it in a data statement.

  !       EIGRF is based primarily on EISPACK routines.  The matrix is
  !       first balanced using the Parlett-Reinsch algorithm.  Then
  !       the Martin-Wilkinson algorithm is applied.

  !       References:
  !          Dongarra, J. and C. Moler, EISPACK -- A Package for Solving
  !             Matrix Eigenvalue Problems, in Cowell, ed., 1984:
  !             Sources and Development of Mathematical Software,
  !             Prentice-Hall, Englewood Cliffs, NJ
  !         Parlett and Reinsch, 1969: Balancing a Matrix for Calculation
  !             of Eigenvalues and Eigenvectors, Num. Math. 13, 293-304
  !         Wilkinson, J., 1965: The Algebraic Eigenvalue Problem,
  !             Clarendon Press, Oxford

  !   I N P U T    V A R I A B L E S:

  !       AA    :  input asymmetric matrix, destroyed after solved
  !        M    :  order of  AA
  !       IA    :  first dimension of  AA
  !    IEVEC    :  first dimension of  EVEC

  !   O U T P U T    V A R I A B L E S:

  !       EVEC  :  (unnormalized) eigenvectors of  AA
  !                   ( column J corresponds to EVAL(J) )

  !       EVAL  :  (unordered) eigenvalues of AA ( dimension at least M )

  !       IER   :  if .NE. 0, signals that EVAL(IER) failed to converge;
  !                   in that case eigenvalues IER+1,IER+2,...,M  are
  !                   correct but eigenvalues 1,...,IER are set to zero.

  !   S C R A T C H   V A R I A B L E S:

  !       WKD   :  work area ( dimension at least 2*M )
  !       AAD   :  double precision stand-in for AA
  !       EVECD :  double precision stand-in for EVEC
  !       EVALD :  double precision stand-in for EVAL

  !   Called by- SOLEIG
  !   Calls- D1MACH, ERRMSG
  ! +-------------------------------------------------------------------+
  IMPLICIT NONE

  INTEGER, INTENT(IN)   :: m
  INTEGER, INTENT(IN)   :: ia
  INTEGER, INTENT(IN)   :: ievec
  REAL, INTENT(IN)      :: aa( ia, m )
  REAL, INTENT(OUT)     :: evec( ievec, m )
  REAL, INTENT(OUT)     :: eval( m )
  INTEGER, INTENT(OUT)  :: ier

   DOUBLE PRECISION, INTENT(OUT)  :: wkd( * )
   DOUBLE PRECISION, INTENT(OUT)  :: aad( ia, m )
   DOUBLE PRECISION, INTENT(OUT)  :: evecd( ia, m )
   DOUBLE PRECISION, INTENT(OUT)  :: evald( m )

  LOGICAL :: noconv, notlas
  INTEGER :: i, iii, in, j, k, ka, kkk, l, lb, lll, n, n1, n2
   DOUBLE PRECISION :: col, discri, f, g, h,  &
      p, q, r, repl, rnorm, row, s, scale, sgn, t,  &
      uu, vv, w, x, y, z

   INTRINSIC ABS, DBLE, MIN, SIGN, SQRT

   DOUBLE PRECISION, PARAMETER :: c1=0.4375D0
   DOUBLE PRECISION, PARAMETER :: c2=0.5D0
   DOUBLE PRECISION, PARAMETER :: c3=0.75D0
   DOUBLE PRECISION, PARAMETER :: c4=0.95D0
   DOUBLE PRECISION, PARAMETER :: c5=16.0d0
   DOUBLE PRECISION, PARAMETER :: c6=256.0d0
   DOUBLE PRECISION, PARAMETER :: zero=0.0d0
   DOUBLE PRECISION, PARAMETER :: one=1.0d0


  ier  = 0

  IF( m < 1 .OR. ia < m .OR. ievec < m )  &
      CALL errmsg( 'ASYMTX--bad input variable(s)', .true. )

  !                           ** Handle 1x1 and 2x2 special cases

  IF( m == 1 ) THEN

    eval( 1 )    = aa( 1, 1 )
    evec( 1, 1 ) = 1.0
    RETURN

  ELSE IF( m == 2 ) THEN

    discri = ( aa( 1,1 ) - aa( 2,2 ) )**2 + 4.*aa( 1, 2 )*aa( 2, 1 )

    IF( discri < 0.0 )  &
        CALL errmsg( 'ASYMTX--complex evals in 2x2 case',.true. )

    sgn  = 1.0

    IF( aa( 1,1 ) < aa( 2,2 ) ) sgn  = - 1.0

    eval( 1 ) = 0.5*( aa( 1,1 ) + aa( 2,2 ) + sgn*SQRT( discri ) )
    eval( 2 ) = 0.5*( aa( 1,1 ) + aa( 2,2 ) - sgn*SQRT( discri ) )
    evec( 1, 1 ) = 1.0
    evec( 2, 2 ) = 1.0

    IF( aa( 1,1 ) == aa( 2,2 ) .AND.  &
          ( aa( 2,1 ) == 0.0 .OR. aa( 1,2 ) == 0.0 ) ) THEN

      rnorm  = ABS( aa( 1,1 ) ) + ABS( aa( 1,2 ) ) +  &
          ABS( aa( 2,1 ) ) + ABS( aa( 2,2 ) )
      w  = tol*rnorm
      evec( 2, 1 ) =   aa( 2, 1 ) / w
      evec( 1, 2 ) = - aa( 1, 2 ) / w

    ELSE

      evec( 2, 1 ) = aa( 2, 1 ) / ( eval( 1 ) - aa( 2,2 ) )
      evec( 1, 2 ) = aa( 1, 2 ) / ( eval( 2 ) - aa( 1,1 ) )

    END IF

    RETURN

  END IF
  !                               ** Put s.p. matrix into d.p. matrix
  DO  j = 1, m

    DO  k = 1, m
      aad( j, k ) = DBLE( aa( j,k ) )
    END DO

  END DO

  !                                ** Initialize output variables
  ier  = 0

  DO  i = 1, m
    evald( i ) = zero

    DO  j = 1, m
      evecd( i, j ) = zero
    END DO

    evecd( i, i ) = one
  END DO

  !                  ** Balance the input matrix and reduce its norm by
  !                     diagonal similarity transformation stored in WK;
  !                     then search for rows isolating an eigenvalue
  !                     and push them down
  rnorm  = zero
  l  = 1
  k  = m

  50 CONTINUE
  kkk  = k

  DO  j = kkk, 1, -1

    row  = zero

    DO  i = 1, k

      IF( i /= j ) row  = row + ABS( aad( j,i ) )

    END DO

    IF( row == zero ) THEN

      wkd( k ) = j

      IF( j /= k ) THEN

        DO  i = 1, k
          repl        = aad( i, j )
          aad( i, j ) = aad( i, k )
          aad( i, k ) = repl
        END DO

        DO  i = l, m
          repl        = aad( j, i )
          aad( j, i ) = aad( k, i )
          aad( k, i ) = repl
        END DO

      END IF

      k  = k - 1
      GO TO  50

    END IF

  END DO
  !                                ** Search for columns isolating an
  !                                   eigenvalue and push them left
  100 CONTINUE
  lll  = l

  DO  j = lll, k

    col  = zero

    DO  i = l, k

      IF( i /= j ) col  = col + ABS( aad( i,j ) )

    END DO

    IF( col == zero ) THEN

      wkd( l ) = j

      IF( j /= l ) THEN

        DO  i = 1, k
          repl        = aad( i, j )
          aad( i, j ) = aad( i, l )
          aad( i, l ) = repl
        END DO

        DO  i = l, m
          repl        = aad( j, i )
          aad( j, i ) = aad( l, i )
          aad( l, i ) = repl
        END DO

      END IF

      l  = l + 1
      GO TO  100

    END IF

  END DO

  !                           ** Balance the submatrix in rows L through K
  DO  i = l, k
    wkd( i ) = one
  END DO

  160 CONTINUE
  noconv = .false.

  DO  i = l, k

    col  = zero
    row  = zero

    DO  j = l, k

      IF( j /= i ) THEN

        col  = col + ABS( aad( j,i ) )
        row  = row + ABS( aad( i,j ) )

      END IF

    END DO

    f  = one
    g  = row / c5
    h  = col + row

    180    CONTINUE
    IF( col < g ) THEN

      f    = f*c5
      col  = col*c6
      GO TO  180

    END IF

    g  = row*c5

    190    CONTINUE
    IF( col >= g ) THEN

      f    = f / c5
      col  = col / c6
      GO TO  190

    END IF
  !                                                ** Now balance
    IF( ( col + row ) / f < c4*h ) THEN

      wkd( i ) = wkd( i )*f
      noconv = .true.

      DO  j = l, m
        aad( i, j ) = aad( i, j ) / f
      END DO

      DO  j = 1, k
        aad( j, i ) = aad( j, i )*f
      END DO

    END IF

  END DO


  IF( noconv ) GO TO  160
  !                                   ** Is A already in Hessenberg form?
  IF( k-1 < l+1 ) GO TO  370

  !                                   ** Transfer A to a Hessenberg form
  DO  n = l + 1, k - 1

    h  = zero
    wkd( n + m ) = zero
    scale  = zero
  !                                                 ** Scale column
    DO  i = n, k
      scale  = scale + ABS( aad( i,n - 1 ) )
    END DO

    IF( scale /= zero ) THEN

      DO  i = k, n, -1
        wkd( i + m ) = aad( i, n - 1 ) / scale
        h  = h + wkd( i + m )**2
      END DO

      g    = - SIGN( SQRT( h ), wkd( n + m ) )
      h    = h - wkd( n + m )*g
      wkd( n + m ) = wkd( n + m ) - g
  !                                            ** Form (I-(U*UT)/H)*A
      DO  j = n, m

        f  = zero

        DO  i = k, n, -1
          f  = f + wkd( i + m )*aad( i, j )
        END DO

        DO  i = n, k
          aad( i, j ) = aad( i, j ) - wkd( i + m )*f / h
        END DO

      END DO
  !                                    ** Form (I-(U*UT)/H)*A*(I-(U*UT)/H)
      DO  i = 1, k

        f  = zero

        DO  j = k, n, -1
          f  = f + wkd( j + m )*aad( i, j )
        END DO

        DO  j = n, k
          aad( i, j ) = aad( i, j ) - wkd( j + m )*f / h
        END DO

      END DO

      wkd( n + m ) = scale*wkd( n + m )
      aad( n, n - 1 ) = scale*g

    END IF

  END DO


  DO  n = k - 2, l, -1

    n1   = n + 1
    n2   = n + 2
    f  = aad( n + 1, n )

    IF( f /= zero ) THEN

      f  = f*wkd( n + 1 + m )

      DO  i = n + 2, k
        wkd( i + m ) = aad( i, n )
      END DO

      IF( n + 1 <= k ) THEN

        DO  j = 1, m

          g  = zero

          DO  i = n + 1, k
            g  = g + wkd( i + m )*evecd( i, j )
          END DO

          g  = g / f

          DO  i = n + 1, k
            evecd( i, j ) = evecd( i, j ) + g*wkd( i + m )
          END DO

        END DO

      END IF

    END IF

  END DO


  370 CONTINUE

  n  = 1

  DO  i = 1, m

    DO  j = n, m
      rnorm  = rnorm + ABS( aad( i,j ) )
    END DO

    n  = i

    IF( i < l .OR. i > k ) evald( i ) = aad( i, i )

  END DO

  n  = k
  t  = zero

  !                                      ** Search for next eigenvalues
  400 CONTINUE
  IF( n < l ) GO TO  550

  in  = 0
  n1  = n - 1
  n2  = n - 2
  !                          ** Look for single small sub-diagonal element
  410 CONTINUE

  DO  i = l, n
    lb  = n + l - i

    IF( lb == l ) GO TO  430

    s  = ABS( aad( lb - 1,lb - 1 ) ) + ABS( aad( lb,lb ) )

    IF( s == zero ) s  = rnorm

    IF( ABS( aad( lb, lb-1 ) ) <= tol*s ) GO TO  430

  END DO


  430 CONTINUE
  x  = aad( n, n )

  IF( lb == n ) THEN
  !                                        ** One eigenvalue found
    aad( n, n ) = x + t
    evald( n ) = aad( n, n )
    n  = n1
    GO TO  400

  END IF

  ! next line has been included to avoid run time error caused by xlf

  IF ( ( n1 <= 0 ).OR.( n <= 0 ) ) THEN
    WRITE(0,*) 'Subscript out of bounds in ASYMTX'
    STOP 9999
  END IF

  y  = aad( n1, n1 )
  w  = aad( n, n1 )*aad( n1, n )

  IF( lb == n1 ) THEN
  !                                        ** Two eigenvalues found
    p  = ( y - x )*c2
    q  = p**2 + w
    z  = SQRT( ABS( q ) )
    aad( n, n ) = x + t
    x  = aad( n, n )
    aad( n1, n1 ) = y + t
  !                                        ** Real pair
    z  = p + SIGN( z, p )
    evald( n1 ) = x + z
    evald( n ) = evald( n1 )

    IF( z /= zero ) evald( n ) = x - w / z

    x  = aad( n, n1 )
  !                                  ** Employ scale factor in case
  !                                     X and Z are very small
    r  = SQRT( x*x + z*z )
    p  = x / r
    q  = z / r
  !                                             ** Row modification
    DO  j = n1, m
      z  = aad( n1, j )
      aad( n1, j ) = q*z + p*aad( n, j )
      aad( n, j ) = q*aad( n, j ) - p*z
    END DO
  !                                             ** Column modification
    DO  i = 1, n
      z  = aad( i, n1 )
      aad( i, n1 ) = q*z + p*aad( i, n )
      aad( i, n ) = q*aad( i, n ) - p*z
    END DO
  !                                          ** Accumulate transformations
    DO  i = l, k
      z  = evecd( i, n1 )
      evecd( i, n1 ) = q*z + p*evecd( i, n )
      evecd( i, n ) = q*evecd( i, n ) - p*z
    END DO

    n  = n2
    GO TO  400

  END IF


  IF( in == 30 ) THEN

  !                    ** No convergence after 30 iterations; set error
  !                       indicator to the index of the current eigenvalue
    ier  = n
    GO TO  700

  END IF
  !                                                  ** Form shift
  IF( in == 10 .OR. in == 20 ) THEN

    t  = t + x

    DO  i = l, n
      aad( i, i ) = aad( i, i ) - x
    END DO

    s  = ABS( aad( n,n1 ) ) + ABS( aad( n1,n2 ) )
    x  = c3*s
    y  = x
    w  = -c1*s**2

  END IF


  in  = in + 1

  !                ** Look for two consecutive small sub-diagonal elements

  ! inhibit vectorization by CF77, as this will cause a run time error

  !DIR$ NEXTSCALAR
  DO  j = lb, n2
    i  = n2 + lb - j
    z  = aad( i, i )
    r  = x - z
    s  = y - z
    p  = ( r*s - w ) / aad( i + 1, i ) + aad( i, i + 1 )
    q  = aad( i + 1, i + 1 ) - z - r - s
    r  = aad( i + 2, i + 1 )
    s  = ABS( p ) + ABS( q ) + ABS( r )
    p  = p / s
    q  = q / s
    r  = r / s

    IF( i == lb ) GO TO  490

    uu   = ABS( aad( i, i-1 ) )*( ABS( q ) + ABS( r ) )
    vv   = ABS( p ) * ( ABS( aad( i-1, i-1 ) ) + ABS( z ) +  &
        ABS( aad( i+1, i+1 ) ) )

    IF( uu <= tol*vv ) GO TO  490

  END DO

  490 CONTINUE
  aad( i+2, i ) = zero

  !                      ** fpp vectorization of this loop triggers
  !                         array bounds errors, so inhibit
  !FPP$ NOVECTOR L
  DO  j = i + 3, n
    aad( j, j - 2 ) = zero
    aad( j, j - 3 ) = zero
  END DO

  !             ** Double QR step involving rows K to N and columns M to N

  DO  ka = i, n1

    notlas = ka /= n1

    IF( ka == i ) THEN

      s  = SIGN( SQRT( p*p + q*q + r*r ), p )

      IF( lb /= i ) aad( ka, ka - 1 ) = -aad( ka, ka - 1 )

    ELSE

      p  = aad( ka, ka - 1 )
      q  = aad( ka + 1, ka - 1 )
      r  = zero

      IF( notlas ) r  = aad( ka + 2, ka - 1 )

      x  = ABS( p ) + ABS( q ) + ABS( r )

      IF( x == zero ) CYCLE

      p  = p / x
      q  = q / x
      r  = r / x
      s  = SIGN( SQRT( p*p + q*q + r*r ), p )
      aad( ka, ka - 1 ) = -s*x

    END IF

    p  = p + s
    x  = p / s
    y  = q / s
    z  = r / s
    q  = q / p
    r  = r / p
  !                                              ** Row modification
    DO  j = ka, m

      p  = aad( ka, j ) + q*aad( ka + 1, j )

      IF( notlas ) THEN

        p  = p + r*aad( ka + 2, j )
        aad( ka + 2, j ) = aad( ka + 2, j ) - p*z

      END IF

      aad( ka + 1, j ) = aad( ka + 1, j ) - p*y
      aad( ka, j ) = aad( ka, j ) - p*x
    END DO
  !                                                 ** Column modification
    DO  iii = 1, MIN( n, ka + 3 )

      p  = x*aad( iii, ka ) + y*aad( iii, ka + 1 )

      IF( notlas ) THEN

        p  = p + z*aad( iii, ka + 2 )
        aad( iii, ka + 2 ) = aad( iii, ka + 2 ) - p*r

      END IF

      aad( iii, ka + 1 ) = aad( iii, ka + 1 ) - p*q
      aad( iii, ka ) = aad( iii, ka ) - p
    END DO
  !                                          ** Accumulate transformations
    DO  iii = l, k

      p  = x*evecd( iii, ka ) + y*evecd( iii, ka + 1 )

      IF( notlas ) THEN

        p  = p + z*evecd( iii, ka + 2 )
        evecd( iii, ka + 2 ) = evecd( iii, ka + 2 ) - p*r

      END IF

      evecd( iii, ka + 1 ) = evecd( iii, ka + 1 ) - p*q
      evecd( iii, ka ) = evecd( iii, ka ) - p
    END DO

  END DO

  GO TO  410
  !                     ** All evals found, now backsubstitute real vector
  550 CONTINUE

  IF( rnorm /= zero ) THEN

    DO  n = m, 1, -1
      n2   = n
      aad( n, n ) = one

      DO  i = n - 1, 1, -1
        w  = aad( i, i ) - evald( n )

        IF( w == zero ) w  = tol*rnorm

        r  = aad( i, n )

        DO  j = n2, n - 1
          r  = r + aad( i, j )*aad( j, n )
        END DO

        aad( i, n ) = -r / w
        n2   = i
      END DO

    END DO
  !                      ** End backsubstitution vectors of isolated evals
    DO  i = 1, m

      IF( i < l .OR. i > k ) THEN

        DO  j = i, m
          evecd( i, j ) = aad( i, j )
        END DO

      END IF

    END DO
  !                                   ** Multiply by transformation matrix
    IF( k /= 0 ) THEN

      DO  j = m, l, -1

        DO  i = l, k
          z  = zero

          DO  n = l, MIN( j, k )
            z  = z + evecd( i, n )*aad( n, j )
          END DO

          evecd( i, j ) = z
        END DO

      END DO

    END IF

  END IF


  DO  i = l, k

    DO  j = 1, m
      evecd( i, j ) = evecd( i, j )*wkd( i )
    END DO
  END DO

  !                           ** Interchange rows if permutations occurred
  DO  i = l-1, 1, -1

    j  = wkd( i )

    IF( i /= j ) THEN

      DO  n = 1, m
        repl   = evecd( i, n )
        evecd( i, n ) = evecd( j, n )
        evecd( j, n ) = repl
      END DO

    END IF

  END DO


  DO  i = k + 1, m

    j  = wkd( i )

    IF( i /= j ) THEN

      DO  n = 1, m
        repl   = evecd( i, n )
        evecd( i, n ) = evecd( j, n )
        evecd( j, n ) = repl
      END DO

    END IF

  END DO

  !                         ** Put results into output arrays
  700 CONTINUE

  DO  j = 1, m

    eval( j ) = evald( j )

    DO  k = 1, m
      evec( j, k ) = evecd( j, k )
    END DO

  END DO

  END SUBROUTINE asymtx

  SUBROUTINE chekin( nlyr, dtauc, ssalb, pmom, temper, wvnmlo,  &
      wvnmhi, usrtau, ntau, utau, nstr, usrang, numu,  &
      umu, nphi, phi, ibcnd, fbeam, umu0, phi0,  &
      fisot, lamber, albedo, hl, btemp, ttemp, temis,  &
      plank, onlyfl, accur, tauc, maxcly, maxulv,  &
      maxumu, maxcmu, maxphi, mxcly, mxulv, mxumu, mxcmu, mxphi )

  !           Checks the input dimensions and variables

  !   Calls- WRTBAD, WRTDIM, DREF, ERRMSG
  !   Called by- DISORT
  ! --------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: maxcly
  INTEGER, INTENT(IN)  :: maxulv
  INTEGER, INTENT(IN)  :: maxumu
  INTEGER, INTENT(IN)  :: maxcmu
  INTEGER, INTENT(IN)  :: maxphi
  INTEGER, INTENT(IN)  :: mxcly
  INTEGER, INTENT(IN)  :: mxulv
  INTEGER, INTENT(IN)  :: mxumu
  INTEGER, INTENT(IN)  :: mxcmu
  INTEGER, INTENT(IN)  :: mxphi
  INTEGER, INTENT(IN)  :: nlyr
  INTEGER, INTENT(IN)  :: ntau
  INTEGER, INTENT(IN)  :: nstr
  INTEGER, INTENT(IN)  :: numu
  INTEGER, INTENT(IN)  :: nphi
  INTEGER, INTENT(IN)  :: ibcnd
  REAL, INTENT(IN)     :: dtauc( maxcly )
  REAL, INTENT(IN)     :: ssalb( maxcly )
  REAL, INTENT(IN)     :: pmom( 0:maxcmu, maxcly )
  REAL, INTENT(IN)     :: temper( 0:maxcly )
  REAL, INTENT(IN)     :: wvnmlo
  REAL, INTENT(IN)     :: wvnmhi
  REAL, INTENT(IN)     :: umu( maxumu )
  REAL, INTENT(IN)     :: phi( maxphi )
  REAL, INTENT(IN)     :: fbeam
  REAL, INTENT(IN)     :: umu0
  REAL, INTENT(IN)     :: phi0
  REAL, INTENT(IN)     :: fisot
  REAL, INTENT(IN)     :: albedo
  REAL, INTENT(IN)     :: hl( 0:maxcmu )
  REAL, INTENT(IN)     :: btemp
  REAL, INTENT(IN)     :: ttemp
  REAL, INTENT(IN)     :: temis
  REAL, INTENT(IN)     :: accur
  REAL, INTENT(IN)     :: tauc( 0:mxcly )
  REAL, INTENT(INOUT)  :: utau( maxulv )
  LOGICAL, INTENT(IN)  :: usrtau
  LOGICAL, INTENT(IN)  :: usrang
  LOGICAL, INTENT(IN)  :: lamber
  LOGICAL, INTENT(IN)  :: plank
  LOGICAL, INTENT(IN)  :: onlyfl

  LOGICAL :: inperr
  INTEGER :: irmu, iu, j, k, lc, lu
  REAL :: flxalb, rmu

  INTRINSIC ABS, MOD

  inperr = .false.

  IF( nlyr < 1 ) inperr = wrtbad( 'NLYR' )

  IF( nlyr > maxcly ) inperr = wrtbad( 'MAXCLY' )

  DO  lc = 1, nlyr

    IF( dtauc( lc ) < 0.0 ) inperr = wrtbad( 'DTAUC' )

    IF( ssalb( lc ) < 0.0 .OR. ssalb( lc ) > 1.0 ) inperr = wrtbad( 'SSALB' )

    IF( plank .AND. ibcnd /= 1 ) THEN

      IF( lc == 1 .AND. temper( 0 ) < 0.0 ) inperr = wrtbad( 'TEMPER' )

      IF( temper( lc ) < 0.0 ) inperr = wrtbad( 'TEMPER' )

    END IF

    DO  k = 0, nstr

      IF( pmom( k,lc ) < -1.0 .OR. pmom( k,lc ) > 1.0 )  &
          inperr = wrtbad( 'PMOM' )

    END DO

  END DO


  IF( ibcnd == 1 ) THEN

    IF( maxulv < 2 ) inperr = wrtbad( 'MAXULV' )

  ELSE IF( usrtau ) THEN

    IF( ntau < 1 ) inperr = wrtbad( 'NTAU' )

    IF( maxulv < ntau ) inperr = wrtbad( 'MAXULV' )

    DO  lu = 1, ntau

      IF( ABS( utau( lu ) - tauc( nlyr ) ) <= 1.e-4 ) utau( lu ) = tauc( nlyr )

      IF( utau( lu ) < 0.0 .OR. utau( lu ) > tauc( nlyr ) )  &
          inperr = wrtbad( 'UTAU' )

    END DO

  ELSE

    IF( maxulv < nlyr + 1 ) inperr = wrtbad( 'MAXULV' )

  END IF


  IF( nstr < 2 .OR. MOD( nstr,2 ) /= 0 ) inperr = wrtbad( 'NSTR' )

  !     IF( NSTR.EQ.2 )
  !    &    CALL ERRMSG( 'CHEKIN--2 streams not recommended;'//
  !    &                 ' use specialized 2-stream code instead',.False.)

  IF( nstr > maxcmu ) inperr = wrtbad( 'MAXCMU' )

  IF( usrang ) THEN

    IF( numu < 0 ) inperr = wrtbad( 'NUMU' )

    IF( .NOT.onlyfl .AND. numu == 0 ) inperr = wrtbad( 'NUMU' )

    IF( numu > maxumu ) inperr = wrtbad( 'MAXUMU' )

    IF( ibcnd == 1 .AND. 2*numu > maxumu ) inperr = wrtbad( 'MAXUMU' )

    DO  iu = 1, numu

      IF( umu( iu ) < -1.0 .OR. umu( iu ) > 1.0 .OR.  &
          umu( iu ) == 0.0 ) inperr = wrtbad( 'UMU' )

      IF( ibcnd == 1 .AND. umu( iu ) < 0.0 ) inperr = wrtbad( 'UMU' )

      IF( iu > 1 ) THEN

        IF( umu( iu ) < umu( iu-1 ) ) inperr = wrtbad( 'UMU' )

      END IF

    END DO

  ELSE

    IF( maxumu < nstr ) inperr = wrtbad( 'MAXUMU' )

  END IF


  IF( .NOT.onlyfl .AND. ibcnd /= 1 ) THEN

    IF( nphi <= 0 ) inperr = wrtbad( 'NPHI' )

    IF( nphi > maxphi ) inperr = wrtbad( 'MAXPHI' )

    DO  j = 1, nphi

      IF( phi( j ) < 0.0 .OR. phi( j ) > 360.0 ) inperr = wrtbad( 'PHI' )

    END DO

  END IF


  IF( ibcnd < 0 .OR. ibcnd > 1 ) inperr = wrtbad( 'ibcnd' )

  IF( ibcnd == 0 ) THEN

    IF( fbeam < 0.0 ) inperr = wrtbad( 'FBEAM' )

    IF( fbeam > 0.0 .AND. ABS(umu0) > 1.0 ) inperr = wrtbad( 'UMU0' )

    IF( fbeam > 0.0 .AND. ( phi0 < 0.0 .OR.phi0 > 360.0 ) )  &
        inperr = wrtbad( 'PHI0' )

    IF( fisot < 0.0 ) inperr = wrtbad( 'FISOT' )

    IF( lamber ) THEN

      IF( albedo < 0.0 .OR. albedo > 1.0 ) inperr = wrtbad( 'ALBEDO' )

    ELSE
  !                    ** Make sure flux albedo at dense mesh of incident
  !                       angles does not assume unphysical values

      DO  irmu = 0, 100
        rmu  = irmu*0.01
        flxalb = dref( rmu, hl, nstr )

        IF( flxalb < 0.0 .OR. flxalb > 1.0 ) inperr = wrtbad( 'HL' )

      END DO

    END IF


  ELSE IF( ibcnd == 1 ) THEN

    IF( albedo < 0.0 .OR. albedo > 1.0 ) inperr = wrtbad( 'ALBEDO' )

  END IF


  IF( plank .AND. ibcnd /= 1 ) THEN

    IF( wvnmlo < 0.0 .OR. wvnmhi <= wvnmlo ) inperr = wrtbad( 'WVNMLO,HI' )

    IF( temis < 0.0 .OR. temis > 1.0 ) inperr = wrtbad( 'temis' )

    IF( btemp < 0.0 ) inperr = wrtbad( 'BTEMP' )

    IF( ttemp < 0.0 ) inperr = wrtbad( 'TTEMP' )

  END IF


  IF( accur < 0.0 .OR. accur > 1.e-2 ) inperr = wrtbad( 'accur' )

  IF( mxcly < nlyr ) inperr = wrtdim( 'MXCLY', nlyr )

  IF( ibcnd /= 1 ) THEN

    IF( usrtau .AND. mxulv < ntau ) inperr = wrtdim( 'MXULV',ntau )

    IF( .NOT.usrtau .AND. mxulv < nlyr + 1 )  &
        inperr = wrtdim( 'MXULV', nlyr + 1 )

  ELSE

    IF( mxulv < 2 ) inperr = wrtdim( 'MXULV', 2 )

  END IF

  IF( mxcmu < nstr ) inperr = wrtdim( 'MXCMU', nstr )

  IF( usrang .AND. mxumu < numu ) inperr = wrtdim( 'MXUMU', numu )

  IF( usrang .AND. ibcnd == 1 .AND.mxumu < 2*numu )  &
      inperr = wrtdim( 'MXUMU', numu )

  IF( .NOT.usrang .AND. mxumu < nstr ) inperr = wrtdim( 'MXUMU', nstr )

  IF( .NOT.onlyfl .AND. ibcnd /= 1 .AND. mxphi < nphi )  &
      inperr = wrtdim( 'MXPHI', nphi )

  IF( inperr ) CALL errmsg( 'DISORT--input and/or dimension errors',.true.)

  IF( plank ) THEN

    DO  lc = 1, nlyr

      IF( ABS( temper( lc ) - temper( lc-1 ) ) > 20.0 )  &
          CALL errmsg('CHEKIN--vertical temperature step may'  &
          // ' be too large for good accuracy', .false.)
    END DO

  END IF

  END SUBROUTINE chekin

  SUBROUTINE fluxes( tausla, tauslau,  &
      cmu, cwt, fbeam, gc, kk, layru, ll, lyrcut,  &
      maxulv, mxcmu, mxulv, ncut, nn, nstr, ntau, pi,  &
      prnt, ssalb, taucpr, umu0, utau, utaupr, xr0,  &
      xr1, zz, zplk0, zplk1, dfdt, flup, fldn, fldir, rfldir, rfldn, uavg, u0c,  &
      uavgso, uavgup, uavgdn, sindir, sinup, sindn)

  !       Calculates the radiative fluxes, mean intensity, and flux
  !       derivative with respect to optical depth from the m=0 intensity
  !       components (the azimuthally-averaged intensity)

  !    I N P U T     V A R I A B L E S:

  !       CMU      :  Abscissae for Gauss quadrature over angle cosine
  !       CWT      :  Weights for Gauss quadrature over angle cosine
  !       GC       :  Eigenvectors at polar quadrature angles, SC(1)
  !       KK       :  Eigenvalues of coeff. matrix in Eq. SS(7)
  !       LAYRU    :  Layer number of user level UTAU
  !       LL       :  Constants of integration in Eq. SC(1), obtained
  !                     by solving scaled version of Eq. SC(5);
  !                     exponential term of Eq. SC(12) not included
  !       LYRCUT   :  Logical flag for truncation of comput. layer
  !       NN       :  Order of double-Gauss quadrature (NSTR/2)
  !       NCUT     :  Number of computational layer where absorption
  !                     optical depth exceeds ABSCUT
  !       TAUCPR   :  Cumulative optical depth (delta-M-scaled)
  !       UTAUPR   :  Optical depths of user output levels in delta-M
  !                     coordinates;  equal to UTAU if no delta-M
  !       XR0      :  Expansion of thermal source function in Eq. SS(14)
  !       XR1      :  Expansion of thermal source function Eqs. SS(16)
  !       ZZ       :  Beam source vectors in Eq. SS(19)
  !       ZPLK0    :  Thermal source vectors Z0, by solving Eq. SS(16)
  !       ZPLK1    :  Thermal source vectors Z1, by solving Eq. SS(16)
  !       (remainder are DISORT input variables)


  !                   O U T P U T     V A R I A B L E S:

  !       U0C      :  Azimuthally averaged intensities
  !                   ( at polar quadrature angles )
  !       (RFLDIR, RFLDN, FLUP, DFDT, UAVG are DISORT output variables)


  !                   I N T E R N A L       V A R I A B L E S:

  !       DIRINT   :  Direct intensity attenuated
  !       FDNTOT   :  Total downward flux (direct + diffuse)
  !       FLDIR    :  Direct-beam flux (delta-M scaled)
  !       FLDN     :  Diffuse down-flux (delta-M scaled)
  !       FNET     :  Net flux (total-down - diffuse-up)
  !       FACT     :  EXP( - UTAUPR / UMU0 )
  !       PLSORC   :  Planck source function (thermal)
  !       ZINT     :  Intensity of m = 0 case, in Eq. SC(1)

  !   Called by- DISORT
  !   Calls- ZEROIT
  ! +-------------------------------------------------------------------+

  INTEGER, INTENT(IN)  :: ncut
  INTEGER, INTENT(IN)  :: nn
  INTEGER, INTENT(IN)  :: nstr
  INTEGER, INTENT(IN)  :: ntau
  INTEGER, INTENT(IN)  :: maxulv
  INTEGER, INTENT(IN)  :: mxcmu
  INTEGER, INTENT(IN)  :: mxulv
  INTEGER, INTENT(IN)  :: layru( mxulv )
  REAL, INTENT(IN)     :: tausla(0:*)
  REAL, INTENT(IN)     :: tauslau(0:*)
  REAL, INTENT(IN)     :: cmu( mxcmu )
  REAL, INTENT(IN)     :: cwt( mxcmu )
  REAL, INTENT(IN)     :: fbeam
  REAL, INTENT(IN)     :: gc( mxcmu, mxcmu, * )
  REAL, INTENT(IN)     :: kk( mxcmu, * )
  REAL, INTENT(IN)     :: ll( mxcmu, * )
  REAL, INTENT(IN)     :: pi
  REAL, INTENT(IN)     :: ssalb( * )
  REAL, INTENT(IN)     :: taucpr( 0:* )
  REAL, INTENT(IN)     :: umu0
  REAL, INTENT(IN)     :: utau( maxulv )
  REAL, INTENT(IN)     :: utaupr( mxulv )
  REAL, INTENT(IN)     :: xr0( * )
  REAL, INTENT(IN)     :: xr1( * )
  REAL, INTENT(IN)     :: zz( mxcmu, * )
  REAL, INTENT(IN)     :: zplk0( mxcmu, * )
  REAL, INTENT(IN)     :: zplk1( mxcmu, * )
  LOGICAL, INTENT(IN)  :: lyrcut
  LOGICAL, INTENT(IN)  :: prnt( * )
  REAL, INTENT(INOUT)  :: uavg( maxulv )
  REAL, INTENT(INOUT)  :: flup( maxulv )
  REAL, INTENT(OUT)    :: dfdt( maxulv )
  REAL, INTENT(OUT)    :: fldn( mxulv )
  REAL, INTENT(OUT)    :: fldir( mxulv )
  REAL, INTENT(OUT)    :: rfldir( maxulv )
  REAL, INTENT(OUT)    :: rfldn( maxulv )
  REAL, INTENT(OUT)    :: u0c( mxcmu, mxulv )
  REAL, INTENT(OUT)    :: uavgso(*)
  REAL, INTENT(OUT)    :: uavgup(*)
  REAL, INTENT(OUT)    :: uavgdn(*)
  REAL, INTENT(OUT)    :: sindir(*)
  REAL, INTENT(OUT)    :: sinup(*)
  REAL, INTENT(OUT)    :: sindn(*)

  INTEGER :: iq, jq, lu, lyu
  REAL :: ang1, ang2, dirint, fact, fdntot, fnet, plsorc, zint
  !     ..

  INTRINSIC ACOS, EXP

  IF( prnt( 2 ) ) WRITE( *, 9000 )
  !                                          ** Zero DISORT output arrays
  CALL zeroit( u0c, mxulv*mxcmu )
  CALL zeroit( fldir, mxulv )
  CALL zeroit( fldn, mxulv )
  CALL  zeroit( uavgso,   maxulv )
  CALL  zeroit( uavgup,   maxulv )
  CALL  zeroit( uavgdn,   maxulv )
  CALL  zeroit( sindir,   maxulv )
  CALL  zeroit( sinup,    maxulv )
  CALL  zeroit( sindn,    maxulv )

  !                                        ** Loop over user levels
  DO  lu = 1, ntau

    lyu  = layru( lu )

    IF( lyrcut .AND. lyu > ncut ) THEN
  !                                                ** No radiation reaches
  !                                                ** this level
      fdntot = 0.0
      fnet   = 0.0
      plsorc = 0.0
      GO TO  70

    END IF

    IF( fbeam > 0.0 ) THEN

      fact  = EXP( - tausla(lu-1) )
      dirint       = fbeam*fact
      fldir( lu )  = umu0*( fbeam*fact )
      rfldir( lu ) = umu0*fbeam * EXP( -tauslau(lu-1) )
      sindir( lu ) = SQRT(1.-umu0*umu0)*fbeam * EXP( -tauslau(lu-1) )

    ELSE

      dirint       = 0.0
      fldir( lu )  = 0.0
      rfldir( lu ) = 0.0
      sindir( lu ) = 0.0

    END IF


    DO  iq = 1, nn

      zint   = 0.0

      DO  jq = 1, nn
        zint   = zint + gc( iq, jq, lyu )*ll( jq, lyu )*  &
            EXP( -kk( jq,lyu )*( utaupr( lu ) - taucpr( lyu ) ) )
      END DO

      DO  jq = nn + 1, nstr
        zint   = zint + gc( iq, jq, lyu )*ll( jq, lyu )*  &
            EXP( -kk( jq,lyu )*( utaupr( lu ) - taucpr( lyu - 1 ) ) )
      END DO

      u0c( iq, lu ) = zint

      IF( fbeam > 0.0 ) u0c( iq, lu ) = zint + zz( iq, lyu )*fact

      u0c( iq, lu ) = u0c( iq, lu ) + zplk0( iq, lyu ) +  &
          zplk1( iq, lyu )*utaupr( lu )
      uavg( lu ) = uavg( lu ) + cwt( nn + 1 - iq )*u0c( iq, lu )
      uavgdn(lu) = uavgdn(lu) + cwt(nn+1-iq) * u0c( iq,lu )
      sindn(lu)  = sindn(lu)  + cwt(nn+1-iq) *  &
          SQRT(1.-cmu(nn+1-iq)*cmu(nn+1-iq))* u0c( iq, lu )
      fldn( lu ) = fldn( lu ) + cwt( nn + 1 - iq )*  &
          cmu( nn + 1 - iq )*u0c( iq, lu )
    END DO


    DO  iq = nn + 1, nstr

      zint   = 0.0

      DO  jq = 1, nn
        zint   = zint + gc( iq, jq, lyu )*ll( jq, lyu )*  &
            EXP( -kk( jq,lyu )*( utaupr( lu ) - taucpr( lyu ) ) )
      END DO

      DO  jq = nn + 1, nstr
        zint   = zint + gc( iq, jq, lyu )*ll( jq, lyu )*  &
            EXP( -kk( jq,lyu )*( utaupr( lu ) - taucpr( lyu - 1 ) ) )
      END DO

      u0c( iq, lu ) = zint

      IF( fbeam > 0.0 ) u0c( iq, lu ) = zint + zz( iq, lyu )*fact

      u0c( iq, lu ) = u0c( iq, lu ) + zplk0( iq, lyu ) +  &
          zplk1( iq, lyu )*utaupr( lu )
      uavg( lu ) = uavg( lu ) + cwt( iq - nn )*u0c( iq, lu )
      uavgup(lu) = uavgup(lu) + cwt(iq-nn) * u0c( iq,lu )
      sinup (lu) = sinup(lu)  + cwt(iq-nn) * SQRT(1.-cmu(iq-nn)*cmu(iq-nn))*  &
          u0c( iq, lu )
      flup( lu ) = flup( lu ) + cwt( iq - nn )*cmu( iq - nn )* u0c( iq, lu )
    END DO


    flup( lu )  = 2.*pi*flup( lu )
    fldn( lu )  = 2.*pi*fldn( lu )
    fdntot      = fldn( lu ) + fldir( lu )
    fnet        = fdntot - flup( lu )
    rfldn( lu ) = fdntot - rfldir( lu )
    uavg( lu )  = ( 2.*pi*uavg( lu ) + dirint ) / ( 4.*pi )
    uavgso( lu ) = dirint / (4.*pi)
    uavgup( lu ) = (2.0 * pi * uavgup(lu) )/ (4.*pi)
    uavgdn( lu)  = (2.0 * pi * uavgdn(lu) )/ (4.*pi)
    sindn ( lu ) = 2.*pi*sindn ( lu )
    sinup ( lu ) = 2.*pi*sinup ( lu )

    plsorc      = xr0( lyu ) + xr1( lyu )*utaupr( lu )
    dfdt( lu )  = ( 1.- ssalb( lyu ) ) * 4.*pi * ( uavg( lu ) - plsorc )

    70    CONTINUE
    IF( prnt( 2 ) ) WRITE( *, FMT = 9010 ) utau( lu ), lyu,  &
        rfldir( lu ), rfldn( lu ), fdntot, flup( lu ), fnet,  &
        uavg( lu ), plsorc, dfdt( lu )

  END DO


  IF( prnt( 3 ) ) THEN

    WRITE( *, FMT = 9020 )

    DO  lu = 1, ntau

      WRITE( *, FMT = 9030 ) utau( lu )

      DO  iq = 1, nn
        ang1   = 180./ pi* ACOS( cmu( 2*nn - iq + 1 ) )
        ang2   = 180./ pi* ACOS( cmu( iq ) )
        WRITE( *, 9040 ) ang1, cmu(2*nn-iq+1), u0c(iq,lu),  &
            ang2, cmu(iq),        u0c(iq+nn,lu)
      END DO

    END DO

  END IF


  9000 FORMAT( //, 21X,  &
      '<----------------------- FLUXES ----------------------->', /,  &
      '   Optical  Compu    Downward    Downward    Downward     ',  &
      ' Upward                    Mean      Planck   d(Net Flux)', /,  &
      '     Depth  Layer      Direct     Diffuse       Total     ',  &
      'Diffuse         Net   Intensity      Source   / d(Op Dep)', / )
  9010 FORMAT( f10.4, i7, 1P, 7E12.3, e14.3 )
  9020 FORMAT( / , / , ' ******** AZIMUTHALLY AVERAGED INTENSITIES',  &
      ' ( at polar quadrature angles ) *******' )
  9030 FORMAT( /, ' Optical depth =', f10.4, //,  &
      '     Angle (deg)   cos(Angle)     Intensity',  &
      '     Angle (deg)   cos(Angle)     Intensity' )
  9040 FORMAT( 2( 0P,f16.4,f13.5,1P,e14.3 ) )

  END SUBROUTINE fluxes

  SUBROUTINE lepoly( nmu, m, maxmu, twonm1, mu, ylm )

  !       Computes the normalized associated Legendre polynomial,
  !       defined in terms of the associated Legendre polynomial
  !       Plm = P-sub-l-super-m as

  !             Ylm(MU) = sqrt( (l-m)!/(l+m)! ) * Plm(MU)

  !       for fixed order m and all degrees from l = m to TWONM1.
  !       When m.GT.0, assumes that Y-sub(m-1)-super(m-1) is available
  !       from a prior call to the routine.

  !       REFERENCE: Dave, J.V. and B.H. Armstrong, Computations of
  !                  High-Order Associated Legendre Polynomials,
  !                  J. Quant. Spectrosc. Radiat. Transfer 10,
  !                  557-562, 1970.  (hereafter D/A)

  !       METHOD: Varying degree recurrence relationship.

  !       NOTE 1: The D/A formulas are transformed by
  !               setting  M = n-1; L = k-1.
  !       NOTE 2: Assumes that routine is called first with  M = 0,
  !               then with  M = 1, etc. up to  M = TWONM1.
  !       NOTE 3: Loops are written in such a way as to vectorize.

  !  I N P U T     V A R I A B L E S:

  !       NMU    :  Number of arguments of YLM
  !       M      :  Order of YLM
  !       MAXMU  :  First dimension of YLM
  !       TWONM1 :  Max degree of YLM
  !       MU(i)  :  Arguments of YLM (i = 1 to NMU)

  !       If M.GT.0, YLM(M-1,i) for i = 1 to NMU is assumed to exist
  !       from a prior call.

  !  O U T P U T     V A R I A B L E:

  !       YLM(l,i) :  l = M to TWONM1, normalized associated Legendre
  !                   polynomials evaluated at argument MU(i)

  !   Called by- DISORT, ALBTRN, SURFAC
  !   Calls- ERRMSG
  ! +-------------------------------------------------------------------+
  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: nmu
  INTEGER, INTENT(IN)    :: m
  INTEGER, INTENT(IN)    :: maxmu
  INTEGER, INTENT(IN)    :: twonm1
  REAL, INTENT(IN)       :: mu( * )
  REAL, INTENT(OUT)      :: ylm( 0:maxmu, * )

  INTEGER :: i, l, ns
  REAL :: tmp1, tmp2

  IF( pass1 ) THEN

    pass1  = .false.

    DO  ns = 1, maxsqt
      sqt( ns ) = SQRT( FLOAT( ns ) )
    END DO

  END IF

  IF( 2*twonm1 > maxsqt )  &
      CALL errmsg('LEPOLY--need to increase param MAXSQT',.true.)


  IF( m == 0 ) THEN
  !                             ** Upward recurrence for ordinary
  !                                Legendre polynomials
    DO  i = 1, nmu
      ylm( 0, i ) = 1.0
      ylm( 1, i ) = mu( i )
    END DO


    DO  l = 2, twonm1

      DO  i = 1, nmu
        ylm( l, i ) = ( ( 2*l - 1 )*mu( i )*ylm( l - 1, i ) -  &
            ( l - 1 )*ylm( l - 2, i ) ) / l
      END DO

    END DO


  ELSE

    DO  i = 1, nmu
  !                               ** Y-sub-m-super-m; derived from
  !                               ** D/A Eqs. (11,12)

      ylm( m, i ) = - sqt( 2*m - 1 ) / sqt( 2*m )*  &
          SQRT( 1.- mu(i)**2 )*ylm( m - 1, i )

  !                              ** Y-sub-(m+1)-super-m; derived from
  !                              ** D/A Eqs.(13,14) using Eqs.(11,12)

      ylm( m + 1, i ) = sqt( 2*m + 1 )*mu( i )*ylm( m, i )

    END DO

  !                                   ** Upward recurrence; D/A EQ.(10)
    DO  l = m + 2, twonm1

      tmp1  = sqt( l - m )*sqt( l + m )
      tmp2  = sqt( l - m - 1 )*sqt( l + m - 1 )

      DO  i = 1, nmu
        ylm( l, i ) = ( ( 2*l - 1 )*mu( i )*ylm( l-1, i ) -  &
            tmp2*ylm( l-2, i ) ) / tmp1
      END DO

    END DO

  END IF


  END SUBROUTINE lepoly

  SUBROUTINE pravin( umu, numu, maxumu, utau, ntau, u0u )

  !        Print azimuthally averaged intensities at user angles

  !   Called by- DISORT

  !     LENFMT   Max number of polar angle cosines UMU that can be
  !                printed on one line, as set in FORMAT statement
  ! --------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: numu
  INTEGER, INTENT(IN)  :: maxumu
  INTEGER, INTENT(IN)  :: ntau
  REAL, INTENT(IN)     :: umu( numu )
  REAL, INTENT(IN)     :: utau( ntau )
  REAL, INTENT(IN)     :: u0u( maxumu, ntau )

  INTEGER :: iu, iumax, iumin, lenfmt, lu, np, npass

  INTRINSIC MIN

  IF( numu < 1 )  RETURN

  WRITE( *, '(//,A)' ) ' *******  AZIMUTHALLY AVERAGED INTENSITIES ' //  &
      '(at user polar angles)  ********'

  lenfmt = 8
  npass  = 1 + (numu-1) / lenfmt

  WRITE( *,'(/,A,/,A)') '   Optical   Polar Angle Cosines', '     Depth'

  DO  np = 1, npass

    iumin  = 1 + lenfmt * ( np - 1 )
    iumax  = MIN( lenfmt*np, numu )
    WRITE( *,'(/,10X,8F14.5)') ( umu(iu), iu = iumin, iumax )

    DO  lu = 1, ntau
      WRITE( *, '(0P,F10.4,1P,8E14.4)' ) utau( lu ),  &
          ( u0u( iu,lu ), iu = iumin, iumax )
    END DO

  END DO


  END SUBROUTINE pravin

  SUBROUTINE prtinp( nlyr, dtauc, dtaucp, ssalb, pmom, temper,  &
      wvnmlo, wvnmhi, ntau, utau, nstr, numu, umu,  &
      nphi, phi, ibcnd, fbeam, umu0, phi0, fisot,  &
      lamber, albedo, hl, btemp, ttemp, temis,  &
      deltam, plank, onlyfl, accur, flyr, lyrcut,  &
      oprim, tauc, taucpr, maxcmu, prtmom )

  !        Print values of input variables

  !   Called by- DISORT
  ! --------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nlyr
  INTEGER, INTENT(IN)  :: nstr
  INTEGER, INTENT(IN)  :: numu
  INTEGER, INTENT(IN)  :: ntau
  INTEGER, INTENT(IN)  :: nphi
  INTEGER, INTENT(IN)  :: ibcnd
  INTEGER, INTENT(IN)  :: maxcmu
  REAL, INTENT(IN)     :: dtauc( * )
  REAL, INTENT(IN)     :: dtaucp( * )
  REAL, INTENT(IN)     :: ssalb( * )
  REAL, INTENT(IN)     :: pmom( 0:maxcmu, * )
  REAL, INTENT(IN)     :: temper( 0:* )
  REAL, INTENT(IN)     :: wvnmlo
  REAL, INTENT(IN)     :: wvnmhi
  REAL, INTENT(IN)     :: utau( * )
  REAL, INTENT(IN)     :: umu( * )
  REAL, INTENT(IN)     :: phi( * )
  REAL, INTENT(IN)     :: fbeam
  REAL, INTENT(IN)     :: umu0
  REAL, INTENT(IN)     :: phi0
  REAL, INTENT(IN)     :: fisot
  REAL, INTENT(IN)     :: albedo
  REAL, INTENT(IN)     :: hl( 0:maxcmu )
  REAL, INTENT(IN)     :: btemp
  REAL, INTENT(IN)     :: ttemp
  REAL, INTENT(IN)     :: temis
  REAL, INTENT(IN)     :: accur
  REAL, INTENT(IN)     :: flyr( * )
  REAL, INTENT(IN)     :: oprim( * )
  REAL, INTENT(IN)     :: tauc( 0:* )
  REAL, INTENT(IN)     :: taucpr( 0:* )
  LOGICAL, INTENT(IN)  :: prtmom
  LOGICAL, INTENT(IN)  :: lamber
  LOGICAL, INTENT(IN)  :: deltam
  LOGICAL, INTENT(IN)  :: plank
  LOGICAL, INTENT(IN)  :: onlyfl
  LOGICAL, INTENT(IN)  :: lyrcut


  INTEGER :: iu, j, k, lc, lu
  REAL :: yessct
  !     ..


  WRITE( *, '(/,A,I4,A,I4)' ) ' No. streams =', nstr,  &
      '     No. computational layers =', nlyr

  IF( ibcnd /= 1 ) WRITE( *, '(I4,A,10F10.4,/,(26X,10F10.4))' )  &
      ntau,' User optical depths :', ( utau(lu), lu = 1, ntau )

  IF( .NOT.onlyfl ) WRITE( *, '(I4,A,10F9.5,/,(31X,10F9.5))' )  &
      numu,' User polar angle cosines :',( umu(iu), iu = 1, numu )

  IF( .NOT.onlyfl .AND. ibcnd /= 1 )  &
      WRITE( *, '(I4,A,10F9.2,/,(28X,10F9.2))' )  &
      nphi,' User azimuthal angles :',( phi(j), j = 1, nphi )

  IF( .NOT.plank .OR. ibcnd == 1 ) WRITE( *, '(A)' ) ' No thermal emission'


  WRITE( *, '(A,I2)' ) ' Boundary condition flag: IBCND =', ibcnd

  IF( ibcnd == 0 ) THEN

    WRITE( *, '(A,1P,E11.3,A,0P,F8.5,A,F7.2,/,A,1P,E11.3)' )  &
        '    Incident beam with intensity =', fbeam,  &
        ' and polar angle cosine = ', umu0, '  and azimuth angle =', phi0,  &
        '    plus isotropic incident intensity =', fisot

    IF( lamber ) WRITE( *, '(A,0P,F8.4)' )  &
        '    Bottom albedo (Lambertian) =', albedo

    IF( .NOT.lamber ) WRITE( *, '(A,/,(10X,10F9.5))' )  &
        '    Legendre coeffs of bottom bidirectional reflectivity :',  &
        ( hl( k ), k = 0, nstr )

    IF( plank ) WRITE( *, '(A,2F14.4,/,A,F10.2,A,F10.2,A,F8.4)' )  &
        '    Thermal emission in wavenumber interval :', wvnmlo, wvnmhi,  &
        '    Bottom temperature =', btemp, '    Top temperature =', ttemp,  &
        '    Top emissivity =',temis

  ELSE IF( ibcnd == 1 ) THEN

    WRITE(*,'(A)') '    Isotropic illumination from top and bottom'
    WRITE( *, '(A,0P,F8.4)' ) '    Bottom albedo (Lambertian) =', albedo
  END IF


  IF( deltam ) WRITE( *, '(A)' ) ' Uses delta-M method'
  IF( .NOT.deltam ) WRITE( *, '(A)' ) ' Does not use delta-M method'


  IF( ibcnd == 1 ) THEN

    WRITE( *, '(A)' ) ' Calculate albedo and transmissivity of'//  &
        ' medium vs. incident beam angle'

  ELSE IF( onlyfl ) THEN

    WRITE( *, '(A)' ) ' Calculate fluxes and azim-averaged intensities only'

  ELSE

    WRITE( *, '(A)' ) ' Calculate fluxes and intensities'

  END IF


  WRITE( *, '(A,1P,E11.2)' )  &
      ' Relative convergence criterion for azimuth series =', accur

  IF( lyrcut ) WRITE( *, '(A)' )  &
      ' Sets radiation = 0 below absorption optical depth 10'


  !                                        ** Print layer variables
  IF( plank ) WRITE( *, FMT = 9180 )
  IF( .NOT.plank ) WRITE( *, FMT = 9190 )

  yessct = 0.0

  DO  lc = 1, nlyr

    yessct = yessct + ssalb( lc )

    IF( plank ) WRITE(*,'(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5,F9.4,F14.3)')  &
        lc, dtauc( lc ), tauc( lc ), ssalb( lc ), flyr( lc ),  &
        dtaucp( lc ), taucpr( lc ), oprim( lc ), pmom(1,lc), temper( lc-1 )

    IF( .NOT.plank ) WRITE(*,'(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5,F9.4)')  &
        lc, dtauc( lc ), tauc( lc ), ssalb( lc ), flyr( lc ),  &
        dtaucp( lc ), taucpr( lc ), oprim( lc ), pmom( 1,lc )
  END DO

  IF( plank ) WRITE( *, '(85X,F14.3)' ) temper( nlyr )


  IF( prtmom .AND. yessct > 0.0 ) THEN

    WRITE( *, '(/,A)' ) ' Layer   Phase Function Moments'

    DO  lc = 1, nlyr

      IF( ssalb( lc ) > 0.0 ) WRITE( *, '(I6,10F11.6,/,(6X,10F11.6))' )  &
          lc, ( pmom( k, lc ), k = 0, nstr )
    END DO

  END IF

  !                ** (Read every other line in these formats)

  9180 FORMAT( /, 37X, '<------------- Delta-M --------------->', /,  &
      '                   Total    Single                           ',  &
      'Total    Single', /, '       Optical   Optical   Scatter   Truncated   ',  &
      'Optical   Optical   Scatter    Asymm', /,  &
      '         Depth     Depth    Albedo    Fraction     ',  &
      'Depth     Depth    Albedo   Factor   Temperature' )
  9190 FORMAT( /, 37X, '<------------- Delta-M --------------->', /,  &
      '                   Total    Single                           ',  &
      'Total    Single', /, '       Optical   Optical   Scatter   Truncated   ',  &
      'Optical   Optical   Scatter    Asymm', /,  &
      '         Depth     Depth    Albedo    Fraction     ',  &
      'Depth     Depth    Albedo   Factor' )

  END SUBROUTINE prtinp

  SUBROUTINE prtint( uu, utau, ntau, umu, numu, phi, nphi, maxulv, maxumu )

  !         Prints the intensity at user polar and azimuthal angles

  !     All arguments are DISORT input or output variables

  !   Called by- DISORT

  !     LENFMT   Max number of azimuth angles PHI that can be printed
  !                on one line, as set in FORMAT statement
  ! +-------------------------------------------------------------------+
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nphi
  INTEGER, INTENT(IN) :: maxulv
  INTEGER, INTENT(IN) :: maxumu
  INTEGER, INTENT(IN) :: ntau
  INTEGER, INTENT(IN) :: numu
  REAL, INTENT(IN)    :: uu( maxumu, maxulv, * )
  REAL, INTENT(IN)    :: utau( * )
  REAL, INTENT(IN)    :: umu( * )
  REAL, INTENT(IN)    :: phi( * )


  INTEGER :: iu, j, jmax, jmin, lenfmt, lu, np, npass

  INTRINSIC MIN
  !     ..


  IF( nphi < 1 )  RETURN

  WRITE( *, '(//,A)' ) ' *********  I N T E N S I T I E S  *********'

  lenfmt = 10
  npass  = 1 + (nphi-1) / lenfmt

  WRITE( *, '(/,A,/,A,/,A)' )  &
      '             Polar   Azimuth angles (degrees)', '   Optical   Angle',  &
      '    Depth   Cosine'

  DO  lu = 1, ntau

    DO  np = 1, npass

      jmin   = 1 + lenfmt * ( np - 1 )
      jmax   = MIN( lenfmt*np, nphi )

      WRITE( *, '(/,18X,10F11.2)' ) ( phi(j), j = jmin, jmax )

      IF( np == 1 ) WRITE( *, '(F10.4,F8.4,1P,10E11.3)' )  &
          utau(lu), umu(1), (uu(1, lu, j), j = jmin, jmax)
      IF( np > 1 ) WRITE( *, '(10X,F8.4,1P,10E11.3)' )  &
          umu(1), (uu(1, lu, j), j = jmin, jmax)

      DO  iu = 2, numu
        WRITE( *, '(10X,F8.4,1P,10E11.3)' )  &
            umu( iu ), ( uu( iu, lu, j ), j = jmin, jmax )
      END DO

    END DO

  END DO


  END SUBROUTINE prtint

  SUBROUTINE qgausn( m, gmu, gwt )

  !       Compute weights and abscissae for ordinary Gaussian quadrature
  !       on the interval (0,1);  that is, such that

  !           sum(i=1 to M) ( GWT(i) f(GMU(i)) )

  !       is a good approximation to

  !           integral(0 to 1) ( f(x) dx )

  !   INPUT :    M       order of quadrature rule

  !   OUTPUT :  GMU(I)   array of abscissae (I = 1 TO M)
  !             GWT(I)   array of weights (I = 1 TO M)

  !   REFERENCE:  Davis, P.J. and P. Rabinowitz, Methods of Numerical
  !                   Integration, Academic Press, New York, pp. 87, 1975

  !   METHOD:  Compute the abscissae as roots of the Legendre
  !            polynomial P-sub-M using a cubically convergent
  !            refinement of Newton's method.  Compute the
  !            weights from EQ. 2.7.3.8 of Davis/Rabinowitz.  Note
  !            that Newton's method can very easily diverge; only a
  !            very good initial guess can guarantee convergence.
  !            The initial guess used here has never led to divergence
  !            even for M up to 1000.

  !   ACCURACY:  relative error no better than TOL or computer
  !              precision (machine epsilon), whichever is larger

  !   INTERNAL VARIABLES:

  !    ITER      : number of Newton Method iterations
  !    MAXIT     : maximum allowed iterations of Newton Method
  !    PM2,PM1,P : 3 successive Legendre polynomials
  !    PPR       : derivative of Legendre polynomial
  !    P2PRI     : 2nd derivative of Legendre polynomial
  !    TOL       : convergence criterion for Legendre poly root iteration
  !    X,XI      : successive iterates in cubically-convergent version
  !                of Newtons Method (seeking roots of Legendre poly.)

  !   Called by- SETDIS, SURFAC
  !   Calls- D1MACH, ERRMSG
  ! +-------------------------------------------------------------------+
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: m
  REAL, INTENT(OUT)   :: gmu( m )
  REAL, INTENT(OUT)   :: gwt( m )

  INTEGER :: iter, k, lim, nn, np1
  REAL :: cona, t
  DOUBLE PRECISION :: en, nnp1, p, p2pri, pm1, pm2, ppr, prod,  &
      tmp, x, xi

  INTRINSIC ABS, ASIN, COS, FLOAT, MOD, TAN

   INTEGER, PARAMETER :: maxit=1000
   DOUBLE PRECISION, PARAMETER :: one=1.0d0
   DOUBLE PRECISION, PARAMETER :: two=2.0d0


  IF( m < 1 ) CALL errmsg( 'QGAUSN--Bad value of M',.true.)

  IF( m == 1 ) THEN

    gmu( 1 ) = 0.5
    gwt( 1 ) = 1.0
    RETURN

  END IF

  en   = m
  np1  = m + 1
  nnp1 = m*np1
  cona = FLOAT( m - 1 ) / ( 8*m**3 )

  lim  = m / 2

  DO  k = 1, lim
  !                                        ** Initial guess for k-th root
  !                                           of Legendre polynomial, from
  !                                           Davis/Rabinowitz (2.7.3.3a)
    t  = ( 4*k - 1 )*pi / ( 4*m + 2 )
    x  = COS( t + cona / TAN( t ) )
    iter = 0
  !                                        ** Upward recurrence for
  !                                           Legendre polynomials
    10    CONTINUE
    iter   = iter + 1
    pm2    = one
    pm1    = x

    DO  nn = 2, m
      p    = ( ( 2*nn - 1 )*x*pm1 - ( nn - 1 )*pm2 ) / nn
      pm2  = pm1
      pm1  = p
    END DO
  !                                              ** Newton Method
    tmp    = one / ( one - x**2 )
    ppr    = en*( pm2 - x*p )*tmp
    p2pri  = ( two*x*ppr - nnp1*p )*tmp
    xi     = x - ( p / ppr )*( one + ( p / ppr )*p2pri / ( two*ppr ) )

  !                                              ** Check for convergence
    IF( ABS( xi - x ) > tol ) THEN

      IF( iter > maxit ) CALL errmsg( 'QGAUSN--max iteration count',.true.)

      x  = xi
      GO TO  10

    END IF
  !                             ** Iteration finished--calculate weights,
  !                                abscissae for (-1,1)
    gmu( k ) = -x
    gwt( k ) = two / ( tmp*( en*pm2 )**2 )
    gmu( np1 - k ) = -gmu( k )
    gwt( np1 - k ) = gwt( k )
  END DO
  !                                    ** Set middle abscissa and weight
  !                                       for rules of odd order
  IF( MOD( m,2 ) /= 0 ) THEN

    gmu( lim + 1 ) = 0.0
    prod   = one

    DO  k = 3, m, 2
      prod   = prod * k / ( k - 1 )
    END DO

    gwt( lim + 1 ) = two / prod**2
  END IF

  !                                        ** Convert from (-1,1) to (0,1)
  DO  k = 1, m
    gmu( k ) = 0.5*gmu( k ) + 0.5
    gwt( k ) = 0.5*gwt( k )
  END DO


  END SUBROUTINE qgausn

  SUBROUTINE setdis( dsdh, nid, tausla, tauslau, mu2,  &
      cmu, cwt, deltam, dtauc, dtaucp, expbea,  &
      flyr, gl, hl, hlpr, ibcnd, lamber, layru,  &
      lyrcut, maxumu, maxcmu, mxcmu, ncut, nlyr,  &
      ntau, nn, nstr, plank, numu, onlyfl, oprim,  &
      pmom, ssalb, tauc, taucpr, utau, utaupr, umu, umu0, usrtau, usrang )

  !          Perform miscellaneous setting-up operations

  !       INPUT :  all are DISORT input variables (see DOC file)

  !       OUTPUT:  NTAU,UTAU   if USRTAU = FALSE
  !                NUMU,UMU    if USRANG = FALSE
  !                CMU,CWT     computational polar angles and
  !                               corresponding quadrature weights
  !                EXPBEA      transmission of direct beam
  !                FLYR        truncated fraction in delta-M method
  !                GL          phase function Legendre coefficients multi-
  !                              plied by (2L+1) and single-scatter albedo
  !                HLPR        Legendre moments of surface bidirectional
  !                              reflectivity, times 2K+1
  !                LAYRU       Computational layer in which UTAU falls
  !                LYRCUT      flag as to whether radiation will be zeroed
  !                              below layer NCUT
  !                NCUT        computational layer where absorption
  !                              optical depth first exceeds  ABSCUT
  !                NN          NSTR / 2
  !                OPRIM       delta-M-scaled single-scatter albedo
  !                TAUCPR      delta-M-scaled optical depth
  !                UTAUPR      delta-M-scaled version of  UTAU

  !   Called by- DISORT
  !   Calls- QGAUSN, ERRMSG
  ! ----------------------------------------------------------------------
    IMPLICIT NONE

  INTEGER, INTENT(IN)  :: maxumu
  INTEGER, INTENT(IN)  :: maxcmu
  INTEGER, INTENT(IN)  :: mxcmu
  INTEGER, INTENT(IN)  :: nlyr
  INTEGER, INTENT(IN)  :: nstr
  INTEGER, INTENT(IN)  :: ibcnd
  INTEGER, INTENT(IN)  :: nid(0:kz)
  REAL, INTENT(IN)     :: dsdh(0:kz,kz)
  REAL, INTENT(IN)     :: dtauc( * )
  REAL, INTENT(IN)     :: hl( 0:maxcmu )
  REAL, INTENT(IN)     :: ssalb( * )
  REAL, INTENT(IN)     :: tauc( 0:* )
  REAL, INTENT(IN)     :: umu0
  LOGICAL, INTENT(IN)  :: deltam
  LOGICAL, INTENT(IN)  :: lamber
  LOGICAL, INTENT(IN)  :: plank
  LOGICAL, INTENT(IN)  :: onlyfl
  LOGICAL, INTENT(IN)  :: usrtau
  LOGICAL, INTENT(IN)  :: usrang
  INTEGER, INTENT(OUT) :: layru( * )
  INTEGER, INTENT(OUT) :: numu
  INTEGER, INTENT(OUT) :: ncut
  INTEGER, INTENT(OUT) :: ntau
  INTEGER, INTENT(OUT) :: nn
  REAL, INTENT(OUT)    :: tausla(0:kz)
  REAL, INTENT(OUT)    :: tauslau(0:kz)
  REAL, INTENT(OUT)    :: mu2(0:kz)
  REAL, INTENT(OUT)    :: cmu( mxcmu )
  REAL, INTENT(OUT)    :: cwt( mxcmu )
  REAL, INTENT(OUT)    :: dtaucp( * )
  REAL, INTENT(OUT)    :: expbea( 0:* )
  REAL, INTENT(OUT)    :: flyr( * )
  REAL, INTENT(OUT)    :: gl( 0:mxcmu, * )
  REAL, INTENT(OUT)    :: hlpr( 0:mxcmu )
  REAL, INTENT(OUT)    :: oprim( * )
  REAL, INTENT(OUT)    :: pmom( 0:maxcmu, * )
  REAL, INTENT(OUT)    :: taucpr( 0:* )
  REAL, INTENT(OUT)    :: utau( * )
  REAL, INTENT(OUT)    :: utaupr( * )
  REAL, INTENT(OUT)    :: umu( maxumu )
  LOGICAL, INTENT(OUT) :: lyrcut


  REAL :: sum, sumu

  INTEGER :: iq, iu, k, lc, lu, i
  REAL :: abstau, f
  REAL,PARAMETER :: abscut=10.0

  INTRINSIC ABS, EXP


  IF( .NOT.usrtau ) THEN
  !                              ** Set output levels at computational
  !                                 layer boundaries
    ntau  = nlyr + 1

    DO  lc = 0, ntau - 1
      utau( lc + 1 ) = tauc( lc )
    END DO

  END IF
  !                        ** Apply delta-M scaling and move description
  !                           of computational layers to local variables
  expbea( 0 ) = 1.0
  taucpr( 0 ) = 0.0
  abstau      = 0.0
  DO i = 0, kz
    tausla( i ) = 0.0
    tauslau( i ) = 0.0
    mu2(i) = 1./largest
  END DO

  DO  lc = 1, nlyr

    pmom( 0, lc ) = 1.0

    IF( abstau < abscut ) ncut  = lc

    abstau = abstau + ( 1.- ssalb( lc ) )*dtauc( lc )

    IF( .NOT.deltam ) THEN

      oprim( lc )  = ssalb( lc )
      dtaucp( lc ) = dtauc( lc )
      taucpr( lc ) = tauc( lc )

      DO  k = 0, nstr - 1
        gl( k, lc ) = ( 2*k + 1 )*oprim( lc )*pmom( k, lc )
      END DO

      f  = 0.0


    ELSE
  !                                    ** Do delta-M transformation

      f  = pmom( nstr, lc )
      oprim(lc) = ssalb(lc) * ( 1.- f ) / ( 1.- f * ssalb(lc) )
      dtaucp( lc ) = ( 1.- f*ssalb( lc ) )*dtauc( lc )
      taucpr( lc ) = taucpr( lc-1 ) + dtaucp( lc )

      DO  k = 0, nstr - 1
        gl( k, lc ) = ( 2*k + 1 ) * oprim( lc ) *  &
            ( pmom( k,lc ) - f ) / ( 1.- f )
      END DO

    END IF

    flyr( lc )   = f
    expbea( lc ) = 0.0

  END DO

  ! calculate slant optical depth

  IF(umu0 < 0.0) THEN
    IF(nid(0) < 0) THEN
      tausla(0) = largest
      tauslau(0) = largest
    ELSE
      sum = 0.0
      sumu = 0.0
      DO lc = 1, nid(0)
        sum = sum + 2.*dtaucp(lc)*dsdh(0,lc)
        sumu = sumu + 2.*dtauc(lc)*dsdh(0,lc)
      END DO
      tausla(0) = sum
      tauslau(0) = sumu
    END IF
  END IF

  expbea( 0 ) = EXP( -tausla( 0 ) )


  DO   lc = 1, nlyr
    IF(nid(lc) < 0) THEN
      tausla(lc) = largest
      tauslau(lc) = largest
    ELSE
      sum = 0.0
      sumu = 0.0
      DO lu = 1, MIN(nid(lc),lc)
        sum = sum + dtaucp(lu)*dsdh(lc,lu)
        sumu = sumu + dtauc(lu)*dsdh(lc,lu)
      END DO
      DO lu = MIN(nid(lc),lc)+1,nid(lc)
        sum = sum + 2.*dtaucp(lu)*dsdh(lc,lu)
        sumu = sumu + 2.*dtauc(lu)*dsdh(lc,lu)
      END DO
      tausla(lc) = sum
      tauslau(lc) = sumu
      IF(tausla(lc) == tausla(lc-1)) THEN
        mu2(lc) = largest
      ELSE
        mu2(lc) = (taucpr(lc)-taucpr(lc-1)) /(tausla(lc)-tausla(lc-1))
        mu2(lc) = SIGN( AMAX1(ABS(mu2(lc)),1./largest), mu2(lc) )
      END IF
    END IF
    expbea(lc) = EXP( -tausla( lc ) )
  END DO

  !                      ** If no thermal emission, cut off medium below
  !                         absorption optical depth = ABSCUT ( note that
  !                         delta-M transformation leaves absorption
  !                         optical depth invariant ).  Not worth the
  !                         trouble for one-layer problems, though.
  lyrcut = .false.

  IF( abstau >= abscut .AND. .NOT.plank .AND. ibcnd /= 1 .AND.  &
      nlyr > 1 ) lyrcut = .true.

  IF( .NOT.lyrcut ) ncut   = nlyr

  !                             ** Set arrays defining location of user
  !                             ** output levels within delta-M-scaled
  !                             ** computational mesh
  DO  lu = 1, ntau

    DO  lc = 1, nlyr

      IF( utau( lu ) >= tauc( lc - 1 ) .AND.  &
          utau( lu ) <= tauc( lc ) ) GO TO  60

    END DO
    lc   = nlyr

    60    CONTINUE
    utaupr( lu ) = utau( lu )
    IF( deltam ) utaupr( lu ) = taucpr( lc - 1 ) +  &
        ( 1.- ssalb( lc )*flyr( lc ) )* ( utau( lu ) - tauc( lc-1 ) )
    layru( lu ) = lc

  END DO
  !                      ** Calculate computational polar angle cosines
  !                         and associated quadrature weights for Gaussian
  !                         quadrature on the interval (0,1) (upward)
  nn   = nstr / 2

  CALL qgausn( nn, cmu, cwt )
  !                                  ** Downward (neg) angles and weights
  DO  iq = 1, nn
    cmu( iq + nn ) = - cmu( iq )
    cwt( iq + nn ) = cwt( iq )
  END DO


  !     IF( FBEAM.GT.0.0 ) THEN
  !                               ** Compare beam angle to comput. angles
  DO  iq = 1, nn

  !                      ** Dither mu2 if it is close to one of the
  !                         quadrature angles.

    DO  lc = 1, nlyr
      IF (  ABS(mu2(lc)) < 1.e5 ) THEN
        IF( ABS( 1. - ABS(mu2(lc))/cmu( iq ) ) < 0.05 ) mu2(lc) = mu2(lc)*0.999
      END IF
    END DO

  END DO

  !     END IF

  IF( .NOT.usrang .OR. ( onlyfl .AND. maxumu >= nstr ) ) THEN

  !                                   ** Set output polar angles to
  !                                      computational polar angles
    numu   = nstr

    DO  iu = 1, nn
      umu( iu ) = - cmu( nn + 1 - iu )
    END DO

    DO  iu = nn + 1, nstr
      umu( iu ) = cmu( iu - nn )
    END DO

  END IF


  IF( usrang .AND. ibcnd == 1 ) THEN

  !                               ** Shift positive user angle cosines to
  !                                  upper locations and put negatives
  !                                  in lower locations
    DO  iu = 1, numu
      umu( iu + numu ) = umu( iu )
    END DO

    DO  iu = 1, numu
      umu( iu ) = -umu( 2*numu + 1 - iu )
    END DO

    numu   = 2*numu

  END IF


  IF( .NOT.lyrcut .AND. .NOT.lamber ) THEN

    DO  k = 0, nstr
      hlpr( k ) = ( 2*k + 1 )*hl( k )
    END DO

  END IF


  END SUBROUTINE setdis

  SUBROUTINE setmtx( bdr, cband, cmu, cwt, delm0, dtaucp, gc, kk,  &
      lamber, lyrcut, mi, mi9m2, mxcmu, ncol, ncut, nnlyri, nn, nstr, taucpr, wk )

  !        Calculate coefficient matrix for the set of equations
  !        obtained from the boundary conditions and the continuity-
  !        of-intensity-at-layer-interface equations;  store in the
  !        special banded-matrix format required by LINPACK routines

  !     I N P U T      V A R I A B L E S:

  !       BDR      :  Surface bidirectional reflectivity
  !       CMU      :  Abscissae for Gauss quadrature over angle cosine
  !       CWT      :  Weights for Gauss quadrature over angle cosine
  !       DELM0    :  Kronecker delta, delta-sub-m0
  !       GC       :  Eigenvectors at polar quadrature angles, SC(1)
  !       KK       :  Eigenvalues of coeff. matrix in Eq. SS(7)
  !       LYRCUT   :  Logical flag for truncation of comput. layer
  !       NN       :  Number of streams in a hemisphere (NSTR/2)
  !       NCUT     :  Total number of computational layers considered
  !       TAUCPR   :  Cumulative optical depth (delta-M-scaled)
  !       (remainder are DISORT input variables)

  !   O U T P U T     V A R I A B L E S:

  !       CBAND    :  Left-hand side matrix of linear system Eq. SC(5),
  !                      scaled by Eq. SC(12); in banded form required
  !                      by LINPACK solution routines
  !       NCOL     :  Counts of columns in CBAND

  !   I N T E R N A L    V A R I A B L E S:

  !       IROW     :  Points to row in CBAND
  !       JCOL     :  Points to position in layer block
  !       LDA      :  Row dimension of CBAND
  !       NCD      :  Number of diagonals below or above main diagonal
  !       NSHIFT   :  For positioning number of rows in band storage
  !       WK       :  Temporary storage for EXP evaluations

  !   Called by- DISORT, ALBTRN
  !   Calls- ZEROIT
  ! +--------------------------------------------------------------------+
  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: mi
  INTEGER, INTENT(IN)    :: mi9m2
  INTEGER, INTENT(IN)    :: mxcmu
  INTEGER, INTENT(IN)    :: ncut
  INTEGER, INTENT(IN)    :: nnlyri
  INTEGER, INTENT(IN)    :: nn
  INTEGER, INTENT(IN)    :: nstr
  REAL, INTENT(IN)       :: bdr( mi, 0:mi )
  REAL, INTENT(IN)       :: cmu( mxcmu )
  REAL, INTENT(IN)       :: cwt( mxcmu )
  REAL, INTENT(IN)       :: delm0
  REAL, INTENT(IN)       :: dtaucp( * )
  REAL, INTENT(IN)       :: gc( mxcmu, mxcmu, * )
  REAL, INTENT(IN)       :: kk( mxcmu, * )
  REAL, INTENT(IN)       :: taucpr( 0:* )
  LOGICAL, INTENT(IN)    :: lamber
  LOGICAL, INTENT(IN)    :: lyrcut
  INTEGER, INTENT(OUT)   :: ncol
  REAL, INTENT(OUT)      :: wk( mxcmu )
  REAL, INTENT(OUT)      :: cband( mi9m2, nnlyri )


  INTEGER :: iq, irow, jcol, jq, k, lc, lda, ncd, nncol, nshift
  REAL :: expa, sum

  INTRINSIC EXP

  CALL zeroit( cband, mi9m2*nnlyri )

  ncd    = 3*nn - 1
  lda    = 3*ncd + 1
  nshift = lda - 2*nstr + 1
  ncol   = 0
  !                         ** Use continuity conditions of Eq. STWJ(17)
  !                            to form coefficient matrix in STWJ(20);
  !                            employ scaling transformation STWJ(22)
  DO  lc = 1, ncut

    DO  iq = 1, nn
      wk( iq ) = EXP( kk( iq,lc )*dtaucp( lc ) )
    END DO

    jcol  = 0

    DO  iq = 1, nn

      ncol  = ncol + 1
      irow  = nshift - jcol

      DO  jq = 1, nstr
        cband( irow + nstr, ncol ) =   gc( jq, iq, lc )
        cband( irow, ncol )        = - gc( jq, iq, lc )*wk( iq )
        irow  = irow + 1
      END DO

      jcol  = jcol + 1

    END DO


    DO  iq = nn + 1, nstr

      ncol  = ncol + 1
      irow  = nshift - jcol

      DO  jq = 1, nstr
        cband( irow + nstr, ncol ) =   gc( jq, iq, lc )* wk( nstr + 1 - iq )
        cband( irow, ncol )        = - gc( jq, iq, lc )
        irow  = irow + 1
      END DO

      jcol  = jcol + 1

    END DO

  END DO
  !                  ** Use top boundary condition of STWJ(20a) for
  !                     first layer

  jcol  = 0

  DO  iq = 1, nn

    expa  = EXP( kk( iq,1 )*taucpr( 1 ) )
    irow  = nshift - jcol + nn

    DO  jq = nn, 1, -1
      cband( irow, jcol + 1 ) = gc( jq, iq, 1 )*expa
      irow  = irow + 1
    END DO

    jcol  = jcol + 1

  END DO


  DO  iq = nn + 1, nstr

    irow  = nshift - jcol + nn

    DO  jq = nn, 1, -1
      cband( irow, jcol + 1 ) = gc( jq, iq, 1 )
      irow  = irow + 1
    END DO

    jcol  = jcol + 1

  END DO
  !                           ** Use bottom boundary condition of
  !                              STWJ(20c) for last layer

  nncol = ncol - nstr
  jcol  = 0

  DO  iq = 1, nn

    nncol  = nncol + 1
    irow   = nshift - jcol + nstr

    DO  jq = nn + 1, nstr

      IF( lyrcut .OR. ( lamber .AND. delm0 == 0 ) ) THEN

  !                          ** No azimuthal-dependent intensity if Lam-
  !                             bert surface; no intensity component if
  !                             truncated bottom layer

        cband( irow, nncol ) = gc( jq, iq, ncut )

      ELSE

        sum  = 0.0

        DO  k = 1, nn
          sum  = sum + cwt( k )*cmu( k )*bdr( jq - nn, k )*  &
              gc( nn + 1 - k, iq, ncut )
        END DO

        cband( irow, nncol ) = gc( jq, iq, ncut ) - ( 1.+ delm0 )*sum
      END IF

      irow  = irow + 1

    END DO

    jcol  = jcol + 1

  END DO


  DO  iq = nn + 1, nstr

    nncol  = nncol + 1
    irow   = nshift - jcol + nstr
    expa   = wk( nstr + 1 - iq )

    DO  jq = nn + 1, nstr

      IF( lyrcut .OR. ( lamber .AND. delm0 == 0 ) ) THEN

        cband( irow, nncol ) = gc( jq, iq, ncut )*expa

      ELSE

        sum  = 0.0

        DO  k = 1, nn
          sum  = sum + cwt( k )*cmu( k )*bdr( jq - nn, k )*  &
              gc( nn + 1 - k, iq, ncut )
        END DO

        cband( irow, nncol ) = ( gc( jq,iq,ncut ) - ( 1.+ delm0 )*sum )*expa
      END IF

      irow  = irow + 1

    END DO

    jcol  = jcol + 1

  END DO

  END SUBROUTINE setmtx


  SUBROUTINE soleig( amb, apb, array, cmu, cwt, gl, mi, mazim,  &
      mxcmu, nn, nstr, ylmc, cc, evecc, eval, kk, gc, aad, eveccd, evald, wkd )

  !         Solves eigenvalue/vector problem necessary to construct
  !         homogeneous part of discrete ordinate solution; STWJ(8b)
  !         ** NOTE ** Eigenvalue problem is degenerate when single
  !                    scattering albedo = 1;  present way of doing it
  !                    seems numerically more stable than alternative
  !                    methods that we tried

  !   I N P U T     V A R I A B L E S:

  !       GL     :  Delta-M scaled Legendre coefficients of phase function
  !                    (including factors 2l+1 and single-scatter albedo)
  !       CMU    :  Computational polar angle cosines
  !       CWT    :  Weights for quadrature over polar angle cosine
  !       MAZIM  :  Order of azimuthal component
  !       NN     :  Half the total number of streams
  !       YLMC   :  Normalized associated Legendre polynomial
  !                    at the quadrature angles CMU
  !       (remainder are DISORT input variables)

  !   O U T P U T    V A R I A B L E S:

  !       CC     :  C-sub-ij in Eq. SS(5); needed in SS(15&18)
  !       EVAL   :  NN eigenvalues of Eq. SS(12) on return from ASYMTX
  !                    but then square roots taken
  !       EVECC  :  NN eigenvectors  (G+) - (G-)  on return
  !                    from ASYMTX ( column j corresponds to EVAL(j) )
  !                    but then  (G+) + (G-)  is calculated from SS(10),
  !                    G+  and  G-  are separated, and  G+  is stacked on
  !                    top of  G-  to form NSTR eigenvectors of SS(7)
  !       GC     :  Permanent storage for all NSTR eigenvectors, but
  !                    in an order corresponding to KK
  !       KK     :  Permanent storage for all NSTR eigenvalues of SS(7),
  !                    but re-ordered with negative values first ( square
  !                    roots of EVAL taken and negatives added )

  !   I N T E R N A L   V A R I A B L E S:

  !       AMB,APB :  Matrices (alpha-beta), (alpha+beta) in reduced
  !                    eigenvalue problem
  !       ARRAY   :  Complete coefficient matrix of reduced eigenvalue
  !                    problem: (alfa+beta)*(alfa-beta)
  !       GPPLGM  :  (G+) + (G-) (cf. Eqs. SS(10-11))
  !       GPMIGM  :  (G+) - (G-) (cf. Eqs. SS(10-11))
  !       WKD     :  Scratch array required by ASYMTX

  !   Called by- DISORT, ALBTRN
  !   Calls- ASYMTX, ERRMSG
  ! +-------------------------------------------------------------------+
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: mi
  INTEGER, INTENT(IN) :: mazim
  INTEGER, INTENT(IN) :: mxcmu
  INTEGER, INTENT(IN) :: nn
  INTEGER, INTENT(IN) :: nstr
  REAL, INTENT(IN)    :: cmu( mxcmu )
  REAL, INTENT(IN)    :: cwt( mxcmu )
  REAL, INTENT(IN)    :: gl( 0:mxcmu )
  REAL, INTENT(IN)    :: ylmc( 0:mxcmu, mxcmu )
  REAL, INTENT(OUT)   :: amb( mi, mi )
  REAL, INTENT(OUT)   :: apb( mi, mi )
  REAL, INTENT(OUT)   :: array( mi, * )
  REAL, INTENT(OUT)   :: cc( mxcmu, mxcmu )
  REAL, INTENT(OUT)   :: evecc( mxcmu, mxcmu )
  REAL, INTENT(OUT)   :: eval( mi )
  REAL, INTENT(OUT)   :: kk( mxcmu )
  REAL, INTENT(OUT)   :: gc( mxcmu, mxcmu )
  DOUBLE PRECISION, INTENT(OUT) :: aad( mi, mi )
  DOUBLE PRECISION, INTENT(OUT) :: eveccd( mi, mi )
  DOUBLE PRECISION, INTENT(OUT) :: evald( mi )
  DOUBLE PRECISION, INTENT(OUT) :: wkd( mxcmu )


  INTEGER :: ier, iq, jq, kq, l
  REAL :: alpha, beta, gpmigm, gpplgm, sum

  INTRINSIC ABS, SQRT

  !                             ** Calculate quantities in Eqs. SS(5-6)
  DO  iq = 1, nn

    DO  jq = 1, nstr

      sum  = 0.0
      DO  l = mazim, nstr - 1
        sum  = sum + gl( l )*ylmc( l, iq )*ylmc( l, jq )
      END DO

      cc( iq, jq ) = 0.5*sum*cwt( jq )

    END DO

    DO  jq = 1, nn
  !                             ** Fill remainder of array using symmetry
  !                                relations  C(-mui,muj) = C(mui,-muj)
  !                                and        C(-mui,-muj) = C(mui,muj)

      cc( iq + nn, jq ) = cc( iq, jq + nn )
      cc( iq + nn, jq + nn ) = cc( iq, jq )

  !                                       ** Get factors of coeff. matrix
  !                                          of reduced eigenvalue problem

      alpha  = cc( iq, jq ) / cmu( iq )
      beta   = cc( iq, jq + nn ) / cmu( iq )
      amb( iq, jq ) = alpha - beta
      apb( iq, jq ) = alpha + beta

    END DO

    amb( iq, iq ) = amb( iq, iq ) - 1.0 / cmu( iq )
    apb( iq, iq ) = apb( iq, iq ) - 1.0 / cmu( iq )

  END DO
  !                      ** Finish calculation of coefficient matrix of
  !                         reduced eigenvalue problem:  get matrix
  !                         product (alfa+beta)*(alfa-beta); SS(12)
  DO  iq = 1, nn
    DO  jq = 1, nn
      sum  = 0.
      DO  kq = 1, nn
        sum  = sum + apb( iq, kq )*amb( kq, jq )
      END DO
      array( iq, jq ) = sum
    END DO
  END DO
  !                      ** Find (real) eigenvalues and eigenvectors

  CALL asymtx( array, evecc, eval, nn, mi, mxcmu, ier, wkd, aad, eveccd, evald )

  IF( ier > 0 ) THEN

    WRITE( *, FMT = '(//,A,I4,A)' ) ' ASYMTX--eigenvalue no. ',  &
        ier, '  didnt converge.  Lower-numbered eigenvalues wrong.'

    CALL errmsg( 'ASYMTX--convergence problems',.true.)

  END IF

  !DIR$ IVDEP
  DO  iq = 1, nn
    eval( iq )    = SQRT( ABS( eval( iq ) ) )
    kk( iq + nn ) = eval( iq )
  !                                      ** Add negative eigenvalue
    kk( nn + 1 - iq ) = -eval( iq )
  END DO

  !                          ** Find eigenvectors (G+) + (G-) from SS(10)
  !                             and store temporarily in APB array
  DO  jq = 1, nn

    DO  iq = 1, nn

      sum  = 0.
      DO  kq = 1, nn
        sum  = sum + amb( iq, kq )*evecc( kq, jq )
      END DO

      apb( iq, jq ) = sum / eval( jq )

    END DO

  END DO


  DO  jq = 1, nn
    DO  iq = 1, nn

      gpplgm = apb( iq, jq )
      gpmigm = evecc( iq, jq )
  !                                ** Recover eigenvectors G+,G- from
  !                                   their sum and difference; stack them
  !                                   to get eigenvectors of full system
  !                                   SS(7) (JQ = eigenvector number)

      evecc( iq,      jq ) = 0.5*( gpplgm + gpmigm )
      evecc( iq + nn, jq ) = 0.5*( gpplgm - gpmigm )

  !                                ** Eigenvectors corresponding to
  !                                   negative eigenvalues (corresp. to
  !                                   reversing sign of 'k' in SS(10) )
      gpplgm = - gpplgm
      evecc(iq,   jq+nn) = 0.5 * ( gpplgm + gpmigm )
      evecc(iq+nn,jq+nn) = 0.5 * ( gpplgm - gpmigm )
      gc( iq+nn,   jq+nn )   = evecc( iq,    jq )
      gc( nn+1-iq, jq+nn )   = evecc( iq+nn, jq )
      gc( iq+nn,   nn+1-jq ) = evecc( iq,    jq+nn )
      gc( nn+1-iq, nn+1-jq ) = evecc( iq+nn, jq+nn )

    END DO

  END DO


  END SUBROUTINE soleig

  SUBROUTINE solve0( b, bdr, bem, bplank, cband, cmu, cwt, expbea,  &
      fbeam, fisot, ipvt, lamber, ll, lyrcut, mazim,  &
      mi, mi9m2, mxcmu, ncol, ncut, nn, nstr, nnlyri,  &
      pi, tplank, taucpr, umu0, z, zz, zplk0, zplk1 )

  !        Construct right-hand side vector B for general boundary
  !        conditions STWJ(17) and solve system of equations obtained
  !        from the boundary conditions and the continuity-of-
  !        intensity-at-layer-interface equations.
  !        Thermal emission contributes only in azimuthal independence.

  !     I N P U T      V A R I A B L E S:

  !       BDR      :  Surface bidirectional reflectivity
  !       BEM      :  Surface bidirectional emissivity
  !       BPLANK   :  Bottom boundary thermal emission
  !       CBAND    :  Left-hand side matrix of linear system Eq. SC(5),
  !                   scaled by Eq. SC(12); in banded form required
  !                   by LINPACK solution routines
  !       CMU      :  Abscissae for Gauss quadrature over angle cosine
  !       CWT      :  Weights for Gauss quadrature over angle cosine
  !       EXPBEA   :  Transmission of incident beam, EXP(-TAUCPR/UMU0)
  !       LYRCUT   :  Logical flag for truncation of comput. layer
  !       MAZIM    :  Order of azimuthal component
  !       ncol     :  Counts of columns in CBAND
  !       NN       :  Order of double-Gauss quadrature (NSTR/2)
  !       NCUT     :  Total number of computational layers considered
  !       TPLANK   :  Top boundary thermal emission
  !       TAUCPR   :  Cumulative optical depth (delta-M-scaled)
  !       ZZ       :  Beam source vectors in Eq. SS(19)
  !       ZPLK0    :  Thermal source vectors Z0, by solving Eq. SS(16)
  !       ZPLK1    :  Thermal source vectors Z1, by solving Eq. SS(16)
  !       (remainder are DISORT input variables)

  !   O U T P U T     V A R I A B L E S:

  !       B        :  Right-hand side vector of Eq. SC(5) going into
  !                   SGBSL; returns as solution vector of Eq. SC(12),
  !                   constants of integration without exponential term

  !      LL        :  Permanent storage for B, but re-ordered

  !   I N T E R N A L    V A R I A B L E S:

  !       IPVT     :  Integer vector of pivot indices
  !       IT       :  Pointer for position in  B
  !       NCD      :  Number of diagonals below or above main diagonal
  !       RCOND    :  Indicator of singularity for CBAND
  !       Z        :  Scratch array required by SGBCO

  !   Called by- DISORT
  !   Calls- ZEROIT, SGBCO, ERRMSG, SGBSL
  ! +-------------------------------------------------------------------+
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: mazim
  INTEGER, INTENT(IN) :: mi
  INTEGER, INTENT(IN) :: mi9m2
  INTEGER, INTENT(IN) :: mxcmu
  INTEGER, INTENT(IN) :: ncut
  INTEGER, INTENT(IN) :: nn
  INTEGER, INTENT(IN) :: nstr
  INTEGER, INTENT(IN) :: nnlyri
  INTEGER, INTENT(IN) :: ncol
  REAL, INTENT(IN)    :: bdr( mi, 0:mi )
  REAL, INTENT(IN)    :: bem( mi )
  REAL, INTENT(IN)    :: bplank
  REAL, INTENT(IN)    :: cmu( mxcmu )
  REAL, INTENT(IN)    :: cwt( mxcmu )
  REAL, INTENT(IN)    :: expbea( 0:* )
  REAL, INTENT(IN)    :: fbeam
  REAL, INTENT(IN)    :: fisot
  REAL, INTENT(IN)    :: pi
  REAL, INTENT(IN)    :: tplank
  REAL, INTENT(IN)    :: taucpr( 0:* )
  REAL, INTENT(IN)    :: umu0
  REAL, INTENT(IN)    :: zz( mxcmu, * )
  REAL, INTENT(IN)    :: zplk0( mxcmu, * )
  REAL, INTENT(IN)    :: zplk1( mxcmu, * )
  LOGICAL, INTENT(IN) :: lamber
  LOGICAL, INTENT(IN) :: lyrcut
  REAL, INTENT(INOUT) :: cband( mi9m2, nnlyri )
  INTEGER, INTENT(OUT):: ipvt( * )
  REAL, INTENT(OUT)   :: ll( mxcmu, * )
  REAL, INTENT(OUT)   :: b( nnlyri )
  REAL, INTENT(OUT)   :: z( nnlyri )

  INTEGER :: ipnt, iq, it, jq, lc, ncd
  REAL :: rcond, sum
  !     ..

  CALL zeroit( b, nnlyri )
  !                              ** Construct B,  STWJ(20a,c) for
  !                                 parallel beam + bottom reflection +
  !                                 thermal emission at top and/or bottom

  IF( mazim > 0 .AND. fbeam > 0.0 ) THEN

  !                                         ** Azimuth-dependent case
  !                                            (never called if FBEAM = 0)
    IF( lyrcut .OR. lamber ) THEN

  !               ** No azimuthal-dependent intensity for Lambert surface;
  !                  no intensity component for truncated bottom layer

      DO  iq = 1, nn
  !                                                  ** Top boundary
        b( iq ) = - zz( nn + 1 - iq, 1 )*expbea( 0 )
  !                                                  ** Bottom boundary

        b( ncol - nn + iq ) = -zz( iq + nn, ncut )*expbea( ncut )

      END DO


    ELSE

      DO  iq = 1, nn

        b( iq ) = - zz( nn + 1 - iq, 1 )*expbea( 0 )

        sum  = 0.
        DO  jq = 1, nn
          sum  = sum + cwt( jq )*cmu( jq )*bdr( iq, jq )*  &
              zz( nn + 1 - jq, ncut )*expbea( ncut )
        END DO

        b( ncol - nn + iq ) = sum
        IF( fbeam > 0.0 ) b( ncol - nn + iq ) = sum +  &
            ( bdr( iq,0 )*umu0*fbeam / pi - zz( iq + nn,ncut ) )* expbea( ncut )

      END DO

    END IF
  !                             ** Continuity condition for layer
  !                                interfaces of Eq. STWJ(20b)
    it   = nn

    DO  lc = 1, ncut - 1

      DO  iq = 1, nstr
        it   = it + 1
        b( it ) = ( zz( iq, lc+1 ) - zz( iq, lc ) )*expbea( lc )
      END DO

    END DO


  ELSE
  !                                   ** Azimuth-independent case

    IF( fbeam == 0.0 ) THEN

      DO  iq = 1, nn
  !                                      ** Top boundary

        b( iq ) = -zplk0( nn + 1 - iq, 1 ) + fisot + tplank

      END DO


      IF( lyrcut ) THEN
  !                               ** No intensity component for truncated
  !                                  bottom layer
        DO  iq = 1, nn
  !                                      ** Bottom boundary

          b( ncol - nn + iq ) = - zplk0( iq + nn, ncut ) -  &
              zplk1( iq + nn, ncut )* taucpr( ncut )
        END DO


      ELSE

        DO  iq = 1, nn

          sum  = 0.
          DO  jq = 1, nn
            sum  = sum + cwt( jq )*cmu( jq )*bdr( iq, jq )*  &
                ( zplk0( nn + 1 - jq,ncut ) +  &
                zplk1( nn + 1 - jq,ncut )*taucpr( ncut ) )
          END DO

          b( ncol - nn + iq ) = 2.*sum + bem( iq )*bplank -  &
              zplk0( iq + nn, ncut ) - zplk1( iq + nn, ncut )*  &
              taucpr( ncut )
        END DO

      END IF
  !                             ** Continuity condition for layer
  !                                interfaces, STWJ(20b)
      it   = nn
      DO  lc = 1, ncut - 1

        DO  iq = 1, nstr
          it   = it + 1
          b( it ) =   zplk0( iq, lc + 1 ) - zplk0( iq, lc ) +  &
              ( zplk1( iq, lc + 1 ) - zplk1( iq, lc ) )* taucpr( lc )
        END DO

      END DO


    ELSE

      DO  iq = 1, nn
        b( iq ) = - zz( nn + 1 - iq, 1 )*expbea( 0 ) -  &
            zplk0( nn + 1 - iq, 1 ) + fisot + tplank
      END DO

      IF( lyrcut ) THEN

        DO  iq = 1, nn
          b(ncol-nn+iq) = - zz(iq+nn, ncut) * expbea(ncut)  &
              - zplk0(iq+nn, ncut) - zplk1(iq+nn, ncut) * taucpr(ncut)
        END DO


      ELSE

        DO  iq = 1, nn

          sum  = 0.
          DO  jq = 1, nn
            sum = sum + cwt(jq) * cmu(jq) * bdr(iq,jq)  &
                * ( zz(nn+1-jq, ncut) * expbea(ncut) + zplk0(nn+1-jq, ncut)  &
                + zplk1(nn+1-jq, ncut) * taucpr(ncut))
          END DO

          b(ncol-nn+iq) = 2.*sum + ( bdr(iq,0) * umu0*fbeam/pi  &
              - zz(iq+nn, ncut) ) * expbea(ncut) + bem(iq) * bplank  &
              - zplk0(iq+nn, ncut) - zplk1(iq+nn, ncut) * taucpr(ncut)
        END DO

      END IF


      it   = nn

      DO  lc = 1, ncut - 1

        DO  iq = 1, nstr

          it   = it + 1
          b(it) = ( zz(iq,lc+1) - zz(iq,lc) ) * expbea(lc)  &
              + zplk0(iq,lc+1) - zplk0(iq,lc) +  &
              ( zplk1(iq,lc+1) - zplk1(iq,lc) ) * taucpr(lc)
        END DO

      END DO

    END IF

  END IF
  !                     ** Find L-U (lower/upper triangular) decomposition
  !                        of band matrix CBAND and test if it is nearly
  !                        singular (note: CBAND is destroyed)
  !                        (CBAND is in LINPACK packed format)
  rcond  = 0.0
  ncd    = 3*nn - 1

  CALL sgbco( cband, mi9m2, ncol, ncd, ncd, ipvt, rcond, z )

  IF( 1.0 + rcond == 1.0 )  &
      CALL errmsg('SOLVE0--SGBCO says matrix near singular',.false.)

  !                   ** Solve linear system with coeff matrix CBAND
  !                      and R.H. side(s) B after CBAND has been L-U
  !                      decomposed.  Solution is returned in B.

  CALL sgbsl( cband, mi9m2, ncol, ncd, ncd, ipvt, b, 0 )

  !                   ** Zero CBAND (it may contain 'foreign'
  !                      elements upon returning from LINPACK);
  !                      necessary to prevent errors

  CALL zeroit( cband, mi9m2*nnlyri )

  DO  lc = 1, ncut

    ipnt  = lc*nstr - nn

    DO  iq = 1, nn
      ll( nn + 1 - iq, lc ) = b( ipnt + 1 - iq )
      ll( iq + nn,     lc ) = b( iq + ipnt )
    END DO

  END DO

  END SUBROUTINE solve0

  SUBROUTINE surfac( albedo, delm0, fbeam, hlpr, lamber, mi, mazim,  &
      mxcmu, mxumu, nn, numu, nstr, onlyfl, umu,  &
      usrang, ylm0, ylmc, ylmu, bdr, emu, bem, rmu )

  !       Specifies user's surface bidirectional properties, STWJ(21)

  !   I N P U T     V A R I A B L E S:

  !       DELM0  :  Kronecker delta, delta-sub-m0
  !       HLPR   :  Legendre moments of surface bidirectional reflectivity
  !                    (with 2K+1 factor included)
  !       MAZIM  :  Order of azimuthal component
  !       NN     :  Order of double-Gauss quadrature (NSTR/2)
  !       YLM0   :  Normalized associated Legendre polynomial
  !                 at the beam angle
  !       YLMC   :  Normalized associated Legendre polynomials
  !                 at the quadrature angles
  !       YLMU   :  Normalized associated Legendre polynomials
  !                 at the user angles
  !       (remainder are DISORT input variables)

  !    O U T P U T     V A R I A B L E S:

  !       BDR :  Surface bidirectional reflectivity (computational angles)
  !       RMU :  Surface bidirectional reflectivity (user angles)
  !       BEM :  Surface directional emissivity (computational angles)
  !       EMU :  Surface directional emissivity (user angles)

  !    I N T E R N A L     V A R I A B L E S:

  !       DREF      Directional reflectivity
  !       NMUG   :  Number of angle cosine quadrature points on (0,1) for
  !                   integrating bidirectional reflectivity to get
  !                   directional emissivity (it is necessary to use a
  !                   quadrature set distinct from the computational
  !                   angles, because the computational angles may not be
  !                   dense enough--NSTR may be too small--to give an
  !                   accurate approximation for the integration).
  !       GMU    :  The NMUG angle cosine quadrature points on (0,1)
  !       GWT    :  The NMUG angle cosine quadrature weights on (0,1)
  !       YLMG   :  Normalized associated Legendre polynomials
  !                   at the NMUG quadrature angles

  !   Called by- DISORT
  !   Calls- QGAUSN, LEPOLY, ZEROIT, ERRMSG
  ! +-------------------------------------------------------------------+
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: mi
  INTEGER, INTENT(IN) :: mazim
  INTEGER, INTENT(IN) :: mxcmu
  INTEGER, INTENT(IN) :: mxumu
  INTEGER, INTENT(IN) :: nn
  INTEGER, INTENT(IN) :: numu
  INTEGER, INTENT(IN) :: nstr
  REAL, INTENT(IN)    :: albedo
  REAL, INTENT(IN)    :: delm0
  REAL, INTENT(IN)    :: fbeam
  REAL, INTENT(IN)    :: hlpr( 0:mxcmu )
  REAL, INTENT(IN)    :: ylm0( 0:mxcmu )
  REAL, INTENT(IN)    :: ylmc( 0:mxcmu, mxcmu )
  REAL, INTENT(IN)    :: ylmu( 0:mxcmu, mxumu )
  REAL, INTENT(IN)    :: umu( * )
  LOGICAL, INTENT(IN) :: lamber
  LOGICAL, INTENT(IN) :: onlyfl
  LOGICAL, INTENT(IN) :: usrang
  REAL, INTENT(OUT)   :: bdr( mi, 0:mi )
  REAL, INTENT(OUT)   :: emu( mxumu )
  REAL, INTENT(OUT)   :: bem( mi )
  REAL, INTENT(OUT)   :: rmu( mxumu, 0:mi )

  INTEGER :: iq, iu, jg, jq, k
  REAL :: dref, sgn, sum

  IF( pass4 ) THEN
    CALL qgausn( nmug, gmu, gwt )
    CALL lepoly( nmug, 0, maxstr, maxstr, gmu, ylmg )
    ! ** Convert Legendre polys. to negative GMU
    sgn  = - 1.0
    DO  k = 0, maxstr
      sgn  = - sgn
      DO  jg = 1, nmug
        ylmg( k, jg ) = sgn*ylmg( k, jg )
      END DO
    END DO
    pass4  = .false.
  END IF

  CALL zeroit( bdr, mi*( mi + 1 ) )
  CALL zeroit( bem, mi )

  IF( lamber .AND. mazim == 0 ) THEN
    DO  iq = 1, nn
      bem( iq ) = 1.- albedo
      DO  jq = 0, nn
        bdr( iq, jq ) = albedo
      END DO
    END DO
  ELSE IF( .NOT.lamber ) THEN
    DO  iq = 1, nn
      DO  jq = 1, nn
        sum  = 0.0
        DO  k = mazim, nstr - 1
          sum  = sum + hlpr( k )*ylmc( k, iq )* ylmc( k, jq + nn )
        END DO
        bdr( iq, jq ) = ( 2.- delm0 )*sum
      END DO
      IF( fbeam > 0.0 ) THEN
        sum  = 0.0
        DO  k = mazim, nstr - 1
          sum  = sum + hlpr( k )*ylmc( k, iq )*ylm0( k )
        END DO
        bdr( iq, 0 ) = ( 2.- delm0 )*sum
      END IF
    END DO


    IF( mazim == 0 ) THEN

      IF( nstr > maxstr )  &
          CALL errmsg('SURFAC--parameter MAXSTR too small',.true.)

  !                              ** Integrate bidirectional reflectivity
  !                                 at reflection polar angles CMU and
  !                                 incident angles GMU to get
  !                                 directional emissivity at
  !                                 computational angles CMU.
      DO  iq = 1, nn

        dref  = 0.0

        DO  jg = 1, nmug

          sum  = 0.0
          DO  k = 0, nstr - 1
            sum  = sum + hlpr( k )*ylmc( k, iq )* ylmg( k, jg )
          END DO

          dref  = dref + 2.*gwt( jg )*gmu( jg )*sum

        END DO

        bem( iq ) = 1.- dref

      END DO

    END IF

  END IF
  !                                       ** Compute surface bidirectional
  !                                          properties at user angles

  IF( .NOT.onlyfl .AND. usrang ) THEN

    CALL zeroit( emu, mxumu )
    CALL zeroit( rmu, mxumu*( mi + 1 ) )

    DO  iu = 1, numu

      IF( umu( iu ) > 0.0 ) THEN

        IF( lamber .AND. mazim == 0 ) THEN

          DO  iq = 0, nn
            rmu( iu, iq ) = albedo
          END DO

          emu( iu ) = 1.- albedo


        ELSE IF( .NOT.lamber ) THEN

          DO  iq = 1, nn

            sum  = 0.0
            DO  k = mazim, nstr - 1
              sum  = sum + hlpr( k )*ylmu( k, iu )* ylmc( k, iq + nn )
            END DO

            rmu( iu, iq ) = ( 2.- delm0 )*sum

          END DO


          IF( fbeam > 0.0 ) THEN

            sum  = 0.0
            DO  k = mazim, nstr - 1
              sum  = sum + hlpr( k )*ylmu( k, iu )*ylm0( k )
            END DO

            rmu( iu, 0 ) = ( 2.- delm0 )*sum

          END IF


          IF( mazim == 0 ) THEN

  !                               ** Integrate bidirectional reflectivity
  !                                  at reflection angles UMU and
  !                                  incident angles GMU to get
  !                                  directional emissivity at
  !                                  user angles UMU.
            dref  = 0.0

            DO  jg = 1, nmug

              sum  = 0.0
              DO  k = 0, nstr - 1
                sum  = sum + hlpr( k )*ylmu( k, iu )* ylmg( k, jg )
              END DO

              dref  = dref + 2.*gwt( jg )*gmu( jg )*sum

            END DO

            emu( iu ) = 1.- dref

          END IF

        END IF

      END IF

    END DO

  END IF

  END SUBROUTINE surfac


  !bm  SOLVEC calls SOLEIG and UPBEAM; if UPBEAM reports a potenially
  !bm  unstable solution, the calculation is repeated with a slightly
  !bm  changed single scattering albedo; this process is iterates
  !bm  until a stable solution is found; as stable solutions may be
  !bm  reached either by increasing or by decreasing the single
  !bm  scattering albedo, both directions are explored ('upward' and
  !bm  'downward' iteration); the solution which required the smaller
  !bm  change in the single scattering albedo is finally returned
  !bm  by SOLVEC.

  SUBROUTINE solvec( amb, apb, array, cmu, cwt, gl, mi,  &
      mazim, mxcmu, nn, nstr, ylm0, ylmc, cc,  &
      evecc, eval, kk, gc, aad, eveccd, evald,  &
      wk, wkd, delm0, fbeam, ipvt, pi, zj, zz,  &
      oprim, lc, dither, mu2, glsave, dgl)
    IMPLICIT NONE

  !gy added glsave and dgl to call to allow adjustable dimensioning

  INTEGER, INTENT(IN)  :: mi
  INTEGER, INTENT(IN)  :: mazim
  INTEGER, INTENT(IN)  :: mxcmu
  INTEGER, INTENT(IN)  :: nn
  INTEGER, INTENT(IN)  :: nstr
  INTEGER, INTENT(IN)  :: lc
  REAL, INTENT(IN)     :: cmu( mxcmu )
  REAL, INTENT(IN)     :: cwt( mxcmu )
  REAL, INTENT(IN)     :: ylm0( 0:mxcmu )
  REAL, INTENT(IN)     :: ylmc( 0:mxcmu, mxcmu )
  REAL, INTENT(IN)     :: oprim
  REAL, INTENT(IN)     :: dither
  REAL, INTENT(IN)     :: mu2
  REAL, INTENT(IN)     :: delm0
  REAL, INTENT(IN)     :: fbeam
  REAL, INTENT(IN)     :: pi
  REAL, INTENT(INOUT)  :: gl( 0:mxcmu )
  INTEGER, INTENT(OUT) :: ipvt( * )
  REAL, INTENT(OUT)    :: evecc( mxcmu, mxcmu )
  REAL, INTENT(OUT)    :: eval( mi )
  REAL, INTENT(OUT)    :: kk( mxcmu )
  REAL, INTENT(OUT)    :: gc( mxcmu, mxcmu )
  REAL, INTENT(OUT)    :: amb( mi, mi )
  REAL, INTENT(OUT)    :: apb( mi, mi )
  REAL, INTENT(OUT)    :: array( mi, * )
  REAL, INTENT(OUT)    :: zj( mxcmu )
  REAL, INTENT(OUT)    :: zz( mxcmu )
  REAL, INTENT(OUT)    :: glsave( 0:mxcmu )
  REAL, INTENT(OUT)    :: dgl( 0:mxcmu )
  REAL, INTENT(OUT)    :: wk( mxcmu )
  REAL, INTENT(OUT)    :: cc( mxcmu, mxcmu )
  DOUBLE PRECISION, INTENT(OUT):: aad( mi, mi )
  DOUBLE PRECISION, INTENT(OUT):: eveccd( mi, mi )
  DOUBLE PRECISION, INTENT(OUT):: evald( mi )
  DOUBLE PRECISION, INTENT(OUT):: wkd( mxcmu )

  !bm   Variabupbeamles for instability fix

  INTEGER :: uagain, dagain
  REAL :: minrcond, add, uadd, dadd, ssa, dssa, factor


  LOGICAL :: done, noup, nodn, debug, instab

  INTEGER :: k

  !bm   reset parameters

  done = .false.
  noup = .false.
  nodn = .false.


  !bm   flag for printing debugging output
  !      DEBUG  = .TRUE.
  debug  = .false.

  !bm   instability parameter; the solution is considered
  !bm   unstable, if the RCOND reported by SGECO is smaller
  !bm   than MINRCOND
  minrcond = 5000. * r1mach(4)

  !bm   if an instability is detected, the single scattering albedo
  !bm   is iterated downwards in steps of DADD and upwards in steps
  !bm   of UADD; in practice, MINRCOND and -MINRCOND should
  !bm   be reasonable choices for these parameters
  dadd    = -minrcond
  uadd    = minrcond

  uagain = 0
  dagain = 0
  add   = dadd


  !bm   save array GL( ) because it will be
  !bm   changed if an iteration should be neccessary
  DO k = mazim, nstr - 1
    glsave( k ) =  gl( k )
  END DO

  ssa = oprim


  !bm   in case of an instability reported by UPBEAM (INSTAB)
  !bm   the single scattering albedo will be changed by a small
  !bm   amount (ADD); this is indicated by DAGAIN or UAGAIN
  !bm   being larger than 0; a change in the single scattering
  !bm   albedo is equivalent to scaling the array GL( )

  666  IF ( dagain > 0 .OR. uagain > 0)  THEN
    factor = (ssa + add) / ssa
    DO k = mazim, nstr - 1
      gl( k ) =  gl( k ) * factor
    END DO

    ssa = ssa + add

  !bm   if the single scattering albedo is now smaller than 0
  !bm   the downward iteration is stopped and upward iteration
  !bm   is forced instead

    IF( ssa < dither) THEN
      nodn = .true.
      dagain = -1
      GO TO 778
    END IF

  !bm   if the single scattering albedo is now larger than its maximum
  !bm   allowed value (1.0 - DITHER), the upward iteration is
  !bm   stopped and downward iteration is forced instead

    IF( ssa > 1.0 - dither) THEN
      noup = .true.
      uagain = -1
      GO TO 888
    END IF
  END IF


  !     ** Solve eigenfunction problem in Eq. STWJ(8B);
  !        return eigenvalues and eigenvectors

  777     CALL soleig( amb, apb, array, cmu, cwt, gl, mi,  &
      mazim, mxcmu, nn, nstr, ylmc, cc, evecc, eval, kk, gc, aad, eveccd, evald,  &
      wkd )

  !     ** Calculate particular solutions of
  !        q.SS(18) for incident beam source

  IF ( fbeam > 0.0 ) THEN
    CALL  upbeam( mu2, array, cc, cmu, delm0, fbeam, gl,  &
        ipvt, mazim, mxcmu, nn, nstr, pi, wk,  &
        ylm0, ylmc, zj, zz, minrcond, instab)
  END IF

  !     ** Calculate particular solutions of
  !        Eq. SS(15) for thermal emission source
  !        (not available in psndo.f)

  !bm   finished if the result is stable on the first try
  IF ( (.NOT. instab) .AND. (uagain == 0) .AND. (dagain == 0)) THEN
    GO TO 999
  END IF

  !bm   downward iteration
  IF( instab .AND. uagain == 0 )  THEN
    dagain = dagain + 1
    GO TO 666
  END IF

  !bm   upward iteration
  IF( instab .AND. uagain > 0 )  THEN
    uagain = uagain + 1
    GO TO 666
  END IF


  !bm   ( DAGAIN .NE. 0 ) at this place means that the downward
  !bm   iteration is finished

  778  IF (dagain /= 0 .AND. uagain == 0) THEN

  !bm   save downward iteration data for later use and
  !bm   restore original input data
    DO k = mazim, nstr - 1
      dgl( k ) =  gl( k )
      gl( k ) =  glsave( k )
    END DO

    dssa = ssa
    ssa = oprim

  !bm   start upward iteration
    add = uadd
    uagain = uagain + 1
    GO TO 666
  END IF

  !bm   both iterations finished
  888  IF (done) THEN
    GO TO 998
  END IF


  !bm  if neither upward nor downward iteration converged, the
  !bm  original conditions are restored and SOLEIG/UPBEAM
  !bm  is called for the last time

  IF (noup .AND. nodn) THEN

    DO k = mazim, nstr - 1
      gl( k ) =  glsave( k )
    END DO

    ssa = oprim

    IF (debug) THEN
      WRITE (*,*) '! *** Neither upward nor downward iteration'
      WRITE (*,*) '! *** converged; using original result.'
    END IF

    done = .true.
    GO TO 777
  END IF

  !bm  if upward iteration did not converge, the stable downward conditions
  !bm  are restored and SOLEIG/UPBEAM is called for the last time
  IF (noup) THEN
    DO k = mazim, nstr - 1
      gl( k ) =  dgl( k )
    END DO

    ssa = dssa

    IF (debug) THEN
      WRITE (*,*) '! *** The upward iteration did not converge.'
      WRITE (*,*) '! *** Had to iterate ', dagain,  &
          ' times in layer LC =', lc,';'
      WRITE (*,*) '! *** changed SSA from ', oprim, ' to ', ssa,','
      WRITE (*,*) '! *** by a factor of ', ssa/oprim
    END IF

    done = .true.
    GO TO 777
  END IF

  !bm  if downward iteration did not converge, we are done
  !bm  (the result of the upward iteration will be used)
  IF (nodn) THEN
    IF (debug) THEN
      WRITE (*,*) '! *** The downward iteration did not converge.'
      WRITE (*,*) '! *** Had to iterate ', uagain,  &
          ' times in layer LC =', lc,';'
      WRITE (*,*) '! *** changed SSA from ', oprim, ' to ', ssa,','
      WRITE (*,*) '! *** by a factor of ', ssa/oprim
    END IF

    done = .true.
    GO TO 998
  END IF


  !bm   if both iterations converged, and if the upward iteration
  !bm   required more steps than the downward iteration, the stable
  !bm   downward conditions are restored and SOLEIG/UPBEAM is
  !bm   called for the last time

  IF (uagain > dagain) THEN
    DO k = mazim, nstr - 1
      gl( k ) =  dgl( k )
    END DO

    ssa = dssa

    IF (debug) THEN
      WRITE (*,*) '! *** Both iterations converged;', ' using downward.'
      WRITE (*,*) '! *** Had to iterate ', dagain,  &
          ' times in layer LC =', lc,';'
      WRITE (*,*) '! *** changed SSA from ', oprim, ' to ', ssa,','
      WRITE (*,*) '! *** by a factor of ', ssa/oprim
    END IF

    done = .true.
    GO TO 777
  ELSE

    IF (debug) THEN
      WRITE (*,*) '! *** Both iterations converged;', ' using upward.'
      WRITE (*,*) '! *** Had to iterate ', uagain,  &
          ' times in layer LC =', lc,';'
      WRITE (*,*) '! *** changed SSA from ', oprim, ' to ', ssa,','
      WRITE (*,*) '! *** by a factor of ', ssa/oprim
    END IF

    done = .true.
    GO TO 998
  END IF

  !bm   finally restore original input data
  998  DO k = mazim, nstr - 1
    gl( k ) =  glsave( k )
  END DO

  999  CONTINUE
  END SUBROUTINE solvec



  SUBROUTINE upbeam( mu2, array, cc, cmu, delm0, fbeam, gl, ipvt, mazim,  &
      mxcmu, nn, nstr, pi, wk, ylm0, ylmc, zj, zz, minrcond, instab )

  !         Finds the incident-beam particular solution of SS(18)

  !   I N P U T    V A R I A B L E S:

  !       CC     :  C-sub-ij in Eq. SS(5)
  !       CMU    :  Abscissae for Gauss quadrature over angle cosine
  !       DELM0  :  Kronecker delta, delta-sub-m0
  !       GL     :  Delta-M scaled Legendre coefficients of phase function
  !                    (including factors 2L+1 and single-scatter albedo)
  !       MAZIM  :  Order of azimuthal component
  !       YLM0   :  Normalized associated Legendre polynomial
  !                    at the beam angle
  !       YLMC   :  Normalized associated Legendre polynomial
  !                    at the quadrature angles
  !       (remainder are DISORT input variables)

  !   O U T P U T    V A R I A B L E S:

  !       ZJ     :  Right-hand side vector X-sub-zero in SS(19); also the
  !                 solution vector Z-sub-zero after solving that system

  !       ZZ     :REAL, INTENT(IN)                         :: pi  Permanent storage for ZJ, but re-ordered

  !   I N T E R N A L    V A R I A B L E S:

  !       ARRAY  :  Coefficient matrix in left-hand side of Eq. SS(19)
  !       IPVT   :  Integer vector of pivot indices required by LINPACK
  !       WK     :  Scratch array required by LINPACK

  !   Called by- DISORT
  !   Calls- SGECO, ERRMSG, SGESL
  ! +-------------------------------------------------------------------+
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: mazim
  INTEGER, INTENT(IN)  :: mxcmu
  INTEGER, INTENT(IN)  :: nn
  INTEGER, INTENT(IN)  :: nstr
  REAL, INTENT(IN)     :: mu2
  REAL, INTENT(IN)     :: cc( mxcmu, mxcmu )
  REAL, INTENT(IN)     :: cmu( mxcmu )
  REAL, INTENT(IN)     :: delm0
  REAL, INTENT(IN)     :: fbeam
  REAL, INTENT(IN)     :: gl( 0:mxcmu )
  REAL, INTENT(IN)     :: pi
  REAL, INTENT(IN)     :: ylm0( 0:mxcmu )
  REAL, INTENT(IN)     :: ylmc( 0:mxcmu, * )
  REAL, INTENT(IN)     :: minrcond
  REAL, INTENT(INOUT)    :: array( mxcmu, mxcmu )
  INTEGER, INTENT(OUT) :: ipvt( * )
  REAL, INTENT(OUT)    :: wk( mxcmu )
  REAL, INTENT(OUT)    :: zj( mxcmu )
  REAL, INTENT(OUT)    :: zz( mxcmu )
  LOGICAL, INTENT(OUT) :: instab


  INTEGER :: iq, job, jq, k
  REAL :: rcond, sum


  DO  iq = 1, nstr

    DO  jq = 1, nstr
      array( iq, jq ) = -cc( iq, jq )
    END DO

    array( iq, iq ) = 1.+ cmu( iq ) / mu2 + array( iq, iq )

    sum  = 0.
    DO  k = mazim, nstr - 1
      sum  = sum + gl( k )*ylmc( k, iq )*ylm0( k )
    END DO

    zj( iq ) = ( 2.- delm0 )*fbeam*sum / ( 4.*pi )
  END DO

  !                  ** Find L-U (lower/upper triangular) decomposition
  !                     of ARRAY and see if it is nearly singular
  !                     (NOTE:  ARRAY is destroyed)
  rcond  = 0.0

  CALL sgeco( array, mxcmu, nstr, ipvt, rcond, wk )

  !bm      IF( 1.0 + RCOND.EQ.1.0 )
  !bm     &    CALL ERRMSG('UPBEAM--SGECO says matrix near singular',.FALSE.)
  !bm
  !bm   replaced original check of RCOND by the following:

  instab = .false.
  IF( ABS(rcond) < minrcond )  THEN
    instab = .true.
    RETURN
  END IF

  !                ** Solve linear system with coeff matrix ARRAY
  !                   (assumed already L-U decomposed) and R.H. side(s)
  !                   ZJ;  return solution(s) in ZJ
  job  = 0

  CALL sgesl( array, mxcmu, nstr, ipvt, zj, job )

  !DIR$ IVDEP
  DO  iq = 1, nn
    zz( iq + nn )     = zj( iq )
    zz( nn + 1 - iq ) = zj( iq + nn )
  END DO

  END SUBROUTINE upbeam


  SUBROUTINE zeroal( nd1, expbea, flyr, oprim, taucpr, xr0, xr1,  &
      nd2, cmu, cwt, psi, wk, z0, z1, zj, nd3, hlpr, ylm0,  &
      nd4, array, cc, evecc, nd5, gl,  &
      nd6, ylmc, nd7, ylmu,  &
      nd8, kk, ll, zz, zplk0, zplk1, nd9, gc,  &
      nd10, layru, utaupr, nd11, gu,  &
      nd12, z0u, z1u, zbeam, nd13, eval,  &
      nd14, amb, apb, nd15, ipvt, z,  &
      nd16, rfldir, rfldn, flup, uavg, dfdt, nd17, albmed, trnmed,  &
      nd18, u0u, nd19, uu )

  !         ZERO ARRAYS; NDn is dimension of all arrays following
  !         it in the argument list

  !   Called by- DISORT
  ! --------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nd1
  INTEGER, INTENT(IN) :: nd2
  INTEGER, INTENT(IN) :: nd3
  INTEGER, INTENT(IN) :: nd4
  INTEGER, INTENT(IN) :: nd5
  INTEGER, INTENT(IN) :: nd6
  INTEGER, INTENT(IN) :: nd7
  INTEGER, INTENT(IN) :: nd8
  INTEGER, INTENT(IN) :: nd9
  INTEGER, INTENT(IN) :: nd10
  INTEGER, INTENT(IN) :: nd11
  INTEGER, INTENT(IN) :: nd12
  INTEGER, INTENT(IN) :: nd13
  INTEGER, INTENT(IN) :: nd14
  INTEGER, INTENT(IN) :: nd15
  INTEGER, INTENT(IN) :: nd16
  INTEGER, INTENT(IN) :: nd17
  INTEGER, INTENT(IN) :: nd18
  INTEGER, INTENT(IN) :: nd19
  INTEGER, INTENT(OUT):: layru( * )
  INTEGER, INTENT(OUT):: ipvt( * )
  REAL, INTENT(OUT)   :: expbea( * )
  REAL, INTENT(OUT)   :: flyr( * )
  REAL, INTENT(OUT)   :: oprim( * )
  REAL, INTENT(OUT)   :: taucpr( * )
  REAL, INTENT(OUT)   :: xr0( * )
  REAL, INTENT(OUT)   :: xr1( * )
  REAL, INTENT(OUT)   :: cmu( * )
  REAL, INTENT(OUT)   :: cwt( * )
  REAL, INTENT(OUT)   :: psi( * )
  REAL, INTENT(OUT)   :: wk( * )
  REAL, INTENT(OUT)   :: z0( * )
  REAL, INTENT(OUT)   :: z1( * )
  REAL, INTENT(OUT)   :: zj( * )
  REAL, INTENT(OUT)   :: hlpr( * )
  REAL, INTENT(OUT)   :: ylm0( * )
  REAL, INTENT(OUT)   :: array( * )
  REAL, INTENT(OUT)   :: cc( * )
  REAL, INTENT(OUT)   :: evecc( * )
  REAL, INTENT(OUT)   :: gl( * )
  REAL, INTENT(OUT)   :: ylmc( * )
  REAL, INTENT(OUT)   :: ylmu( * )
  REAL, INTENT(OUT)   :: kk( * )
  REAL, INTENT(OUT)   :: ll( * )
  REAL, INTENT(OUT)   :: zz( * )
  REAL, INTENT(OUT)   :: zplk0( * )
  REAL, INTENT(OUT)   :: zplk1( * )
  REAL, INTENT(OUT)   :: gc( * )
  REAL, INTENT(OUT)   :: utaupr( * )
  REAL, INTENT(OUT)   :: gu( * )
  REAL, INTENT(OUT)   :: z0u( * )
  REAL, INTENT(OUT)   :: z1u( * )
  REAL, INTENT(OUT)   :: zbeam( * )
  REAL, INTENT(OUT)   :: eval( * )
  REAL, INTENT(OUT)   :: amb( * )
  REAL, INTENT(OUT)   :: apb( * )
  REAL, INTENT(OUT)   :: z( * )
  REAL, INTENT(OUT)   :: rfldir( * )
  REAL, INTENT(OUT)   :: rfldn( * )
  REAL, INTENT(OUT)   :: flup( * )
  REAL, INTENT(OUT)   :: uavg( * )
  REAL, INTENT(OUT)   :: dfdt( * )
  REAL, INTENT(OUT)   :: albmed( * )
  REAL, INTENT(OUT)   :: trnmed( * )
  REAL, INTENT(OUT)   :: u0u( * )
  REAL, INTENT(OUT)   :: uu( * )

  INTEGER :: n

  DO  n = 1, nd1
    expbea( n ) = 0.0
    flyr( n )   = 0.0
    oprim( n )  = 0.0
    taucpr( n ) = 0.0
    xr0( n )    = 0.0
    xr1( n )    = 0.0
  END DO

  DO  n = 1, nd2
    cmu( n ) = 0.0
    cwt( n ) = 0.0
    psi( n ) = 0.0
    wk( n )  = 0.0
    z0( n )  = 0.0
    z1( n )  = 0.0
    zj( n )  = 0.0
  END DO

  DO  n = 1, nd3
    hlpr( n ) = 0.0
    ylm0( n ) = 0.0
  END DO

  DO  n = 1, nd4
    array( n ) = 0.0
    cc( n )    = 0.0
    evecc( n ) = 0.0
  END DO

  DO  n = 1, nd5
    gl( n ) = 0.0
  END DO

  DO  n = 1, nd6
    ylmc( n ) = 0.0
  END DO

  DO  n = 1, nd7
    ylmu( n ) = 0.0
  END DO

  DO  n = 1, nd8
    kk( n )    = 0.0
    ll( n )    = 0.0
    zz( n )    = 0.0
    zplk0( n ) = 0.0
    zplk1( n ) = 0.0
  END DO

  DO  n = 1, nd9
    gc( n ) = 0.0
  END DO

  DO  n = 1, nd10
    layru( n )  = 0
    utaupr( n ) = 0.0
  END DO

  DO  n = 1, nd11
    gu( n ) = 0.0
  END DO

  DO  n = 1, nd12
    z0u( n )   = 0.0
    z1u( n )   = 0.0
    zbeam( n ) = 0.0
  END DO

  DO  n = 1, nd13
    eval( n ) = 0.0
  END DO

  DO  n = 1, nd14
    amb( n ) = 0.0
    apb( n ) = 0.0
  END DO

  DO  n = 1, nd15
    ipvt( n ) = 0
    z( n )    = 0.0
  END DO

  DO  n = 1, nd16
    rfldir( n ) = 0.
    rfldn( n )  = 0.
    flup( n )   = 0.
    uavg( n )   = 0.
    dfdt( n )   = 0.
  END DO

  DO  n = 1, nd17
    albmed( n ) = 0.
    trnmed( n ) = 0.
  END DO

  DO  n = 1, nd18
    u0u( n ) = 0.
  END DO

  DO  n = 1, nd19
    uu( n ) = 0.
  END DO

  END SUBROUTINE zeroal

  SUBROUTINE zeroit( a, length )

  !         Zeros a real array A having LENGTH elements
  ! --------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: length
  REAL, INTENT(OUT)    :: a( length )

  INTEGER :: l
  !     ..

  DO  l = 1, length
    a( l ) = 0.0
  END DO

  END SUBROUTINE zeroit

   REAL FUNCTION dref( mu, hl, nstr )

  !        Exact flux albedo for given angle of incidence, given
  !        a bidirectional reflectivity characterized by its
  !        Legendre coefficients ( NOTE** these will only agree
  !        with bottom-boundary albedos calculated by DISORT in
  !        the limit as number of streams go to infinity, because
  !        DISORT evaluates the integral 'CL' only approximately,
  !        by quadrature, while this routine calculates it exactly.)

  !  INPUT :   MU     Cosine of incidence angle
  !            HL     Legendre coefficients of bidirectional reflectivity
  !          NSTR     Number of elements of HL to consider

  !  INTERNAL VARIABLES (P-sub-L is the L-th Legendre polynomial) :

  !       CL      Integral from 0 to 1 of  MU * P-sub-L(MU)
  !                   (vanishes for  L = 3, 5, 7, ... )
  !       PL      P-sub-L
  !       PLM1    P-sub-(L-1)
  !       PLM2    P-sub-(L-2)

  !   Called by- CHEKIN
  !   Calls- ERRMSG
  ! +-------------------------------------------------------------------+
  IMPLICIT NONE
  INTEGER, INTENT(IN)                      :: nstr
  REAL, INTENT(IN)                         :: mu
  REAL, INTENT(IN)                         :: hl( 0:nstr )

  INTEGER :: l
  REAL :: cl, pl, plm1, plm2

  IF( pass5 ) THEN
    pass5  = .false.
    cl     = 0.125
    c( 2 ) = 10.*cl
    DO  l = 4, maxtrm, 2
      cl     = - cl*( l - 3 ) / ( l + 2 )
      c( l ) = 2.*( 2*l + 1 )*cl
    END DO
  END IF

  IF( nstr < 2 .OR. ABS(mu) > 1.0 )  &
      CALL errmsg( 'DREF--input argument error(s)',.true. )

  IF( nstr > maxtrm ) CALL errmsg( 'DREF--parameter MAXTRM too small',.true. )

  dref  = hl( 0 ) - 2.*hl( 1 )*mu
  plm2  = 1.0
  plm1  = - mu

  DO  l = 2, nstr - 1
  !                                ** Legendre polynomial recurrence

    pl = ( ( 2*l - 1 )*( -mu )*plm1 - ( l-1 )*plm2 ) / l

    IF( MOD( l,2 ) == 0 ) dref   = dref + c( l )*hl( l )*pl

    plm2  = plm1
    plm1  = pl

  END DO

  IF( dref < 0.0 .OR. dref > 1.0 )  &
      CALL errmsg( 'DREF--albedo value not in (0,1)',.false. )

  END FUNCTION dref

   REAL FUNCTION ratio( a, b )

  !        Calculate ratio  A/B  with over- and under-flow protection
  !        (thanks to Prof. Jeff Dozier for some suggestions here).
  !        Since this routine takes two logs, it is no speed demon,
  !        but it is invaluable for comparing results from two runs
  !        of a program under development.

  !        NOTE:  In Fortran90, built-in functions TINY and HUGE
  !               can replace the R1MACH calls.
  ! ---------------------------------------------------------------
  IMPLICIT NONE
  REAL, INTENT(IN)                     :: a
  REAL, INTENT(IN)                     :: b

  REAL :: absa, absb, powa, powb

  INTRINSIC ABS, LOG10, SIGN
  !     ..

  IF( pass6 ) THEN

    tinyVar   = r1mach( 1 )
    hugeVar   = r1mach( 2 )
    powmax = LOG10( hugeVar )
    powmin = LOG10( tinyVar )
    pass6  = .false.

  END IF


  IF( a == 0.0 ) THEN

    IF( b == 0.0 ) THEN

      ratio  = 1.0

    ELSE

      ratio  = 0.0

    END IF


  ELSE IF( b == 0.0 ) THEN

    ratio  = SIGN( hugeVar, a )

  ELSE

    absa   = ABS( a )
    absb   = ABS( b )
    powa   = LOG10( absa )
    powb   = LOG10( absb )

    IF( absa < tinyVar .AND. absb < tinyVar ) THEN

      ratio  = 1.0

    ELSE IF( powa - powb >= powmax ) THEN

      ratio  = hugeVar

    ELSE IF( powa - powb <= powmin ) THEN

      ratio  = tinyVar

    ELSE

      ratio  = absa / absb

    END IF
  !                      ** DONT use old trick of determining sign
  !                      ** from A*B because A*B may (over/under)flow

    IF( ( a > 0.0 .AND. b < 0.0 ) .OR.  &
        ( a < 0.0 .AND. b > 0.0 ) ) ratio = -ratio

  END IF

  END FUNCTION ratio

  SUBROUTINE  errmsg( messag, fatal )

  !        Print out a warning or error message;  abort if error
  !        after making symbolic dump (machine-specific)
  IMPLICIT NONE

  CHARACTER (LEN=*), INTENT(IN)   :: messag
  LOGICAL, INTENT(IN)             :: fatal
  !LOGICAL :: cray

  IF ( fatal )  THEN
    WRITE ( *, '(//,2A,//)' )  ' ******* ERROR >>>>>>  ', messag
    STOP
  END IF

  nummsg = nummsg + 1
  IF( msglim )  RETURN

  IF ( nummsg <= maxmsg )  THEN
    WRITE ( *, '(/,2A,/)' )  ' ******* WARNING >>>>>>  ', messag
  ELSE
    WRITE ( *,99 )
    msglim = .true.
  END IF

  RETURN

  99 FORMAT( //,' >>>>>>  TOO MANY WARNING MESSAGES --  ',  &
      'They will no longer be printed  <<<<<<<', // )
  END SUBROUTINE  errmsg


   LOGICAL FUNCTION  wrtbad ( varnam )
  !          Write names of erroneous variables and return 'TRUE'

  !      INPUT :   VarNam = Name of erroneous variable to be written
  !                         ( CHARACTER, any length )
  IMPLICIT NONE

  CHARACTER (LEN=*), INTENT(IN)        :: varnam

  wrtbad = .true.
  nummsg2 = nummsg2 + 1
  WRITE ( *, '(3A)' )  ' ****  Input variable  ', varnam, '  in error  ****'
  IF ( nummsg2 == maxmsg2 )  &
      CALL  errmsg ( 'Too many input errors.  Aborting...', .true. )

  END FUNCTION  wrtbad

   LOGICAL FUNCTION  wrtdim ( dimnam, minval )

  !          Write name of too-small symbolic dimension and
  !          the value it should be increased to;  return 'TRUE'

  !      INPUT :  DimNam = Name of symbolic dimension which is too small
  !                        ( CHARACTER, any length )
  !               Minval = Value to which that dimension should be
  !                        increased (at least)
  IMPLICIT NONE

  CHARACTER (LEN=*), INTENT(IN)        :: dimnam
  INTEGER, INTENT(IN)                  :: minval

  WRITE ( *, '(3A,I7)' )  ' ****  Symbolic dimension  ', dimnam,  &
      '  should be increased to at least ', minval

  wrtdim = .true.
  END FUNCTION  wrtdim

   LOGICAL FUNCTION  tstbad( varnam, relerr )

  !       Write name (VarNam) of variable failing self-test and its
  !       percent error from the correct value;  return  'FALSE'.
  IMPLICIT NONE

  CHARACTER (LEN=*), INTENT(IN)        :: varnam
  REAL, INTENT(IN)                     :: relerr

  tstbad = .false.
  WRITE( *, '(/,3A,1P,E11.2,A)' )  &
      ' Output variable ', varnam,' differed by ', 100.*relerr,  &
      ' per cent from correct value.  Self-test failed.'

  END FUNCTION  tstbad

  SUBROUTINE  sgbco( abd, lda, n, ml, mu, ipvt, rcond, z )

  !         FACTORS A REAL BAND MATRIX BY GAUSSIAN ELIMINATION
  !         AND ESTIMATES THE CONDITION OF THE MATRIX.

  !         REVISION DATE:  8/1/82
  !         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

  !     IF  RCOND  IS NOT NEEDED, SGBFA IS SLIGHTLY FASTER.
  !     TO SOLVE  A*X = B , FOLLOW SBGCO BY SGBSL.

  !     INPUT:

  !        ABD     REAL(LDA, N)
  !                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS
  !                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND
  !                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
  !                ML+1 THROUGH 2*ML+MU+1 OF  ABD .
  !                SEE THE COMMENTS BELOW FOR DETAILS.

  !        LDA     INTEGER
  !                THE LEADING DIMENSION OF THE ARRAY  ABD .
  !                LDA MUST BE .GE. 2*ML + MU + 1 .

  !        N       INTEGER
  !                THE ORDER OF THE ORIGINAL MATRIX.

  !        ML      INTEGER
  !                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
  !                0 .LE. ML .LT. N .

  !        MU      INTEGER
  !                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
  !                0 .LE. MU .LT. N .
  !                MORE EFFICIENT IF  ML .LE. MU .

  !     ON RETURN

  !        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
  !                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
  !                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
  !                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
  !                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.

  !        IPVT    INTEGER(N)
  !                AN INTEGER VECTOR OF PIVOT INDICES.

  !        RCOND   REAL
  !                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
  !                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
  !                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
  !                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
  !                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
  !                           1.0 + RCOND .EQ. 1.0
  !                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
  !                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
  !                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
  !                UNDERFLOWS.

  !        Z       REAL(N)
  !                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
  !                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
  !                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
  !                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .

  !     BAND STORAGE

  !           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT
  !           WILL SET UP THE INPUT.

  !                   ML = (BAND WIDTH BELOW THE DIAGONAL)
  !                   MU = (BAND WIDTH ABOVE THE DIAGONAL)
  !                   M = ML + MU + 1
  !                   DO 20 J = 1, N
  !                      I1 = MAX0(1, J-MU)
  !                      I2 = MIN0(N, J+ML)
  !                      DO 10 I = I1, I2
  !                         K = I - J + M
  !                         ABD(K,J) = A(I,J)
  !                10    CONTINUE
  !                20 CONTINUE

  !           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .
  !           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR
  !           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.
  !           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .
  !           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE
  !           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.

  !     EXAMPLE:  IF THE ORIGINAL MATRIX IS

  !           11 12 13  0  0  0
  !           21 22 23 24  0  0
  !            0 32 33 34 35  0
  !            0  0 43 44 45 46
  !            0  0  0 54 55 56
  !            0  0  0  0 65 66

  !      THEN  N = 6, ML = 1, MU = 2, LDA .GE. 5  AND ABD SHOULD CONTAIN

  !            *  *  *  +  +  +  , * = NOT USED
  !            *  * 13 24 35 46  , + = USED FOR PIVOTING
  !            * 12 23 34 45 56
  !           11 22 33 44 55 66
  !           21 32 43 54 65  *


  !     ROUTINES CALLED:  FROM LINPACK: SGBFA
  !                       FROM BLAS:    SAXPY, SDOT, SSCAL, SASUM
  !                       FROM FORTRAN: ABS, AMAX1, MAX0, MIN0, SIGN

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lda
  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(IN) :: ml
  INTEGER, INTENT(IN) :: mu
  REAL, INTENT(INOUT)    :: abd(lda,*)
  INTEGER, INTENT(OUT) :: ipvt(*)
  REAL, INTENT(OUT)   :: rcond
  REAL, INTENT(OUT)   :: z(*)

  REAL :: ek, t, wk, wkm
  REAL :: anorm, s, sm, ynorm
  INTEGER :: is, info, j, ju, k, kb, kp1, l, la, lm, lz, m, mm


  !                       ** COMPUTE 1-NORM OF A
  anorm = 0.0E0
  l = ml + 1
  is = l + mu
  DO  j = 1, n
    anorm = AMAX1(anorm, sasum(l,abd(is,j), 1))
    IF (is > ml + 1) is = is - 1
    IF (j <= mu) l = l + 1
    IF (j >= n - ml) l = l - 1
  END DO
  !                                               ** FACTOR
  CALL sgbfa(abd, lda, n, ml, mu, ipvt, info)

  !     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
  !     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
  !     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
  !     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
  !     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
  !     OVERFLOW.

  !                     ** SOLVE TRANS(U)*W = E
  ek = 1.0E0
  DO  j = 1, n
    z(j) = 0.0E0
  END DO

  m = ml + mu + 1
  ju = 0
  DO  k = 1, n
    IF (z(k) /= 0.0E0) ek = SIGN(ek, -z(k))
    IF (ABS(ek-z(k)) > ABS(abd(m,k))) THEN
      s = ABS(abd(m,k))/ABS(ek-z(k))
      CALL sscal(n, s, z, 1)
      ek = s*ek
    END IF
    wk = ek - z(k)
    wkm = -ek - z(k)
    s = ABS(wk)
    sm = ABS(wkm)
    IF (abd(m,k) /= 0.0E0) THEN
      wk  = wk /abd(m,k)
      wkm = wkm/abd(m,k)
    ELSE
      wk  = 1.0E0
      wkm = 1.0E0
    END IF
    kp1 = k + 1
    ju = MIN0(MAX0(ju, mu+ipvt(k)), n)
    mm = m
    IF (kp1 <= ju) THEN
      DO  j = kp1, ju
        mm = mm - 1
        sm = sm + ABS(z(j)+wkm*abd(mm,j))
        z(j) = z(j) + wk*abd(mm,j)
        s = s + ABS(z(j))
      END DO
      IF (s < sm) THEN
        t = wkm - wk
        wk = wkm
        mm = m
        DO  j = kp1, ju
          mm = mm - 1
          z(j) = z(j) + t*abd(mm,j)
        END DO
      END IF
    END IF
    z(k) = wk
  END DO

  s = 1.0E0 / sasum(n, z, 1)
  CALL sscal(n, s, z, 1)

  !                         ** SOLVE TRANS(L)*Y = W
  DO  kb = 1, n
    k = n + 1 - kb
    lm = MIN0(ml, n-k)
    IF (k < n) z(k) = z(k) + sdot(lm, abd(m+1,k), 1, z(k+1), 1)
    IF (ABS(z(k)) > 1.0E0) THEN
      s = 1.0E0 / ABS(z(k))
      CALL sscal(n, s, z, 1)
    END IF
    l = ipvt(k)
    t = z(l)
    z(l) = z(k)
    z(k) = t
  END DO

  s = 1.0E0 / sasum(n, z, 1)
  CALL sscal(n, s, z, 1)

  ynorm = 1.0E0
  !                         ** SOLVE L*V = Y
  DO  k = 1, n
    l = ipvt(k)
    t = z(l)
    z(l) = z(k)
    z(k) = t
    lm = MIN0(ml, n-k)
    IF (k < n) CALL saxpy(lm, t, abd(m+1,k), 1, z(k+1), 1)
    IF (ABS(z(k)) > 1.0E0) THEN
      s = 1.0E0 / ABS(z(k))
      CALL sscal(n, s, z, 1)
      ynorm = s*ynorm
    END IF
  END DO

  s = 1.0E0/sasum(n, z, 1)
  CALL sscal(n, s, z, 1)
  ynorm = s*ynorm
  !                           ** SOLVE  U*Z = W
  DO  kb = 1, n
    k = n + 1 - kb
    IF (ABS(z(k)) > ABS(abd(m,k))) THEN
      s = ABS(abd(m,k)) / ABS(z(k))
      CALL sscal(n, s, z, 1)
      ynorm = s*ynorm
    END IF
    IF (abd(m,k) /= 0.0E0) z(k) = z(k)/abd(m,k)
    IF (abd(m,k) == 0.0E0) z(k) = 1.0E0
    lm = MIN0(k, m) - 1
    la = m - lm
    lz = k - lm
    t = -z(k)
    CALL saxpy(lm, t, abd(la,k), 1, z(lz), 1)
  END DO
  !                              ** MAKE ZNORM = 1.0
  s = 1.0E0 / sasum(n, z, 1)
  CALL sscal(n, s, z, 1)
  ynorm = s*ynorm

  IF (anorm /= 0.0E0) rcond = ynorm/anorm
  IF (anorm == 0.0E0) rcond = 0.0E0

  END SUBROUTINE  sgbco

  SUBROUTINE  sgbfa( abd, lda, n, ml, mu, ipvt, info )

  !         FACTORS A REAL BAND MATRIX BY ELIMINATION.

  !         REVISION DATE:  8/1/82
  !         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

  !     SGBFA IS USUALLY CALLED BY SBGCO, BUT IT CAN BE CALLED
  !     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.

  !     INPUT:  SAME AS 'SGBCO'

  !     ON RETURN:

  !        ABD,IPVT    SAME AS 'SGBCO'

  !        INFO    INTEGER
  !                = 0  NORMAL VALUE.
  !                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
  !                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
  !                     INDICATE THAT SGBSL WILL DIVIDE BY ZERO IF
  !                     CALLED.  USE  RCOND  IN SBGCO FOR A RELIABLE
  !                     INDICATION OF SINGULARITY.

  !     (SEE 'SGBCO' FOR DESCRIPTION OF BAND STORAGE MODE)

  !     ROUTINES CALLED:  FROM BLAS:    SAXPY, SSCAL, ISAMAX
  !                       FROM FORTRAN: MAX0, MIN0
  IMPLICIT NONE

  INTEGER, INTENT(IN)   :: n
  INTEGER, INTENT(IN)   :: ml
  INTEGER, INTENT(IN)   :: mu
  INTEGER, INTENT(IN)   :: lda
  REAL, INTENT(INOUT)     :: abd(lda,*)
  INTEGER, INTENT(OUT)  :: ipvt(*)
  INTEGER, INTENT(OUT)  :: info

  REAL :: t
  INTEGER :: i,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1


  m = ml + mu + 1
  info = 0
  !                        ** ZERO INITIAL FILL-IN COLUMNS
  j0 = mu + 2
  j1 = MIN0(n, m) - 1
  DO  jz = j0, j1
    i0 = m + 1 - jz
    DO  i = i0, ml
      abd(i,jz) = 0.0E0
    END DO
  END DO
  jz = j1
  ju = 0

  !                       ** GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
  nm1 = n - 1
  DO  k = 1, nm1
    kp1 = k + 1
  !                                  ** ZERO NEXT FILL-IN COLUMN
    jz = jz + 1
    IF (jz <= n) THEN
      DO  i = 1, ml
        abd(i,jz) = 0.0E0
      END DO
    END IF
  !                                  ** FIND L = PIVOT INDEX
    lm = MIN0(ml, n-k)
    l = isamax(lm+1, abd(m,k), 1) + m - 1
    ipvt(k) = l + k - m

    IF (abd(l,k) == 0.0E0) THEN
  !                                      ** ZERO PIVOT IMPLIES THIS COLUMN
  !                                      ** ALREADY TRIANGULARIZED
      info = k
    ELSE
  !                                ** INTERCHANGE IF NECESSARY
      IF (l /= m) THEN
        t = abd(l,k)
        abd(l,k) = abd(m,k)
        abd(m,k) = t
      END IF
  !                                   ** COMPUTE MULTIPLIERS
      t = -1.0E0 / abd(m,k)
      CALL sscal(lm, t, abd(m+1,k), 1)

  !                               ** ROW ELIMINATION WITH COLUMN INDEXING

      ju = MIN0(MAX0(ju, mu+ipvt(k)), n)
      mm = m
      DO  j = kp1, ju
        l = l - 1
        mm = mm - 1
        t = abd(l,j)
        IF (l /= mm) THEN
          abd(l,j) = abd(mm,j)
          abd(mm,j) = t
        END IF
        CALL saxpy(lm, t, abd(m+1,k), 1, abd(mm+1,j), 1)
      END DO

    END IF

  END DO

  ipvt(n) = n
  IF (abd(m,n) == 0.0E0) info = n

  END SUBROUTINE  sgbfa

  SUBROUTINE  sgbsl( abd, lda, n, ml, mu, ipvt, b, job )

  !         SOLVES THE REAL BAND SYSTEM
  !            A * X = B  OR  TRANSPOSE(A) * X = B
  !         USING THE FACTORS COMPUTED BY SBGCO OR SGBFA.

  !         REVISION DATE:  8/1/82
  !         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

  !     INPUT:

  !        ABD     REAL(LDA, N)
  !                THE OUTPUT FROM SBGCO OR SGBFA.

  !        LDA     INTEGER
  !                THE LEADING DIMENSION OF THE ARRAY  ABD .

  !        N       INTEGER
  !                THE ORDER OF THE ORIGINAL MATRIX.

  !        ML      INTEGER
  !                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.

  !        MU      INTEGER
  !                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.

  !        IPVT    INTEGER(N)
  !                THE PIVOT VECTOR FROM SBGCO OR SGBFA.

  !        B       REAL(N)
  !                THE RIGHT HAND SIDE VECTOR.

  !        JOB     INTEGER
  !                = 0         TO SOLVE  A*X = B ,
  !                = NONZERO   TO SOLVE  TRANS(A)*X = B , WHERE
  !                            TRANS(A)  IS THE TRANSPOSE.

  !     ON RETURN

  !        B       THE SOLUTION VECTOR  X .

  !     ERROR CONDITION

  !        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
  !        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
  !        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
  !        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
  !        CALLED CORRECTLY AND IF SBGCO HAS SET RCOND .GT. 0.0
  !        OR SGBFA HAS SET INFO .EQ. 0 .

  !     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
  !     WITH  P  COLUMNS
  !           CALL SGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
  !           IF (RCOND IS TOO SMALL) GO TO ...
  !           DO 10 J = 1, P
  !              CALL SGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
  !        10 CONTINUE

  !     ROUTINES CALLED:  FROM BLAS:    SAXPY, SDOT
  !                       FROM FORTRAN: MIN0
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: lda
  INTEGER, INTENT(IN)  :: n
  INTEGER, INTENT(IN)  :: ml
  INTEGER, INTENT(IN)  :: mu
  INTEGER, INTENT(IN)  :: job
  INTEGER, INTENT(IN)  :: ipvt(*)
  REAL, INTENT(IN)     :: abd(lda,*)
  REAL, INTENT(IN OUT) :: b(*)

  REAL :: t
  INTEGER :: k,kb,l,la,lb,lm,m,nm1


  m = mu + ml + 1
  nm1 = n - 1
  IF (job == 0) THEN
  !                               ** JOB = 0 , SOLVE  A * X = B
  !                               ** FIRST SOLVE L*Y = B
    IF (ml /= 0) THEN
      DO  k = 1, nm1
        lm = MIN0(ml, n-k)
        l = ipvt(k)
        t = b(l)
        IF (l /= k) THEN
          b(l) = b(k)
          b(k) = t
        END IF
        CALL saxpy( lm, t, abd(m+1,k), 1, b(k+1), 1 )
      END DO
    END IF
  !                           ** NOW SOLVE  U*X = Y
    DO  kb = 1, n
      k = n + 1 - kb
      b(k) = b(k) / abd(m,k)
      lm = MIN0(k, m) - 1
      la = m - lm
      lb = k - lm
      t = -b(k)
      CALL saxpy(lm, t, abd(la,k), 1, b(lb), 1)
    END DO

  ELSE
  !                          ** JOB = NONZERO, SOLVE  TRANS(A) * X = B
  !                                  ** FIRST SOLVE  TRANS(U)*Y = B
    DO  k = 1, n
      lm = MIN0(k, m) - 1
      la = m - lm
      lb = k - lm
      t = sdot(lm, abd(la,k), 1, b(lb), 1)
      b(k) = (b(k) - t)/abd(m,k)
    END DO
  !                                  ** NOW SOLVE TRANS(L)*X = Y
    IF (ml /= 0) THEN
      DO  kb = 1, nm1
        k = n - kb
        lm = MIN0(ml, n-k)
        b(k) = b(k) + sdot(lm, abd(m+1,k), 1, b(k+1), 1)
        l = ipvt(k)
        IF (l /= k) THEN
          t = b(l)
          b(l) = b(k)
          b(k) = t
        END IF
      END DO
    END IF

  END IF

  END SUBROUTINE  sgbsl

  SUBROUTINE  sgeco( a, lda, n,ipvt, rcond, z )

  !         FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION
  !         AND ESTIMATES THE CONDITION OF THE MATRIX.

  !         REVISION DATE:  8/1/82
  !         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

  !         IF  RCOND  IS NOT NEEDED, SGEFA IS SLIGHTLY FASTER.
  !         TO SOLVE  A*X = B , FOLLOW SGECO BY SGESL.

  !     ON ENTRY

  !        A       REAL(LDA, N)
  !                THE MATRIX TO BE FACTORED.

  !        LDA     INTEGER
  !                THE LEADING DIMENSION OF THE ARRAY  A .

  !        N       INTEGER
  !                THE ORDER OF THE MATRIX  A .

  !     ON RETURN

  !        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
  !                WHICH WERE USED TO OBTAIN IT.
  !                THE FACTORIZATION CAN BE WRITTEN  A = L*U , WHERE
  !                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
  !                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.

  !        IPVT    INTEGER(N)
  !                AN INTEGER VECTOR OF PIVOT INDICES.

  !        RCOND   REAL
  !                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
  !                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
  !                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
  !                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
  !                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
  !                           1.0 + RCOND .EQ. 1.0
  !                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
  !                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
  !                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
  !                UNDERFLOWS.

  !        Z       REAL(N)
  !                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
  !                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
  !                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
  !                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .

  !     ROUTINES CALLED:  FROM LINPACK: SGEFA
  !                       FROM BLAS:    SAXPY, SDOT, SSCAL, SASUM
  !                       FROM FORTRAN: ABS, AMAX1, SIGN
  IMPLICIT NONE

  INTEGER, INTENT(IN)   :: lda
  INTEGER, INTENT(IN)   :: n
  REAL, INTENT(INOUT)   :: a(lda,*)
  INTEGER, INTENT(OUT)  :: ipvt(*)
  REAL, INTENT(OUT)     :: rcond
  REAL, INTENT(OUT)     :: z(*)

  REAL :: ek,t,wk,wkm
  REAL :: anorm,s,sm,ynorm
  INTEGER :: info,j,k,kb,kp1,l

  !                        ** COMPUTE 1-NORM OF A
  anorm = 0.0E0
  DO  j = 1, n
    anorm = AMAX1( anorm, sasum(n,a(1,j),1) )
  END DO
  !                                      ** FACTOR
  CALL sgefa(a,lda,n,ipvt,info)

  !     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
  !     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
  !     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
  !     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
  !     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
  !     OVERFLOW.

  !                        ** SOLVE TRANS(U)*W = E
  ek = 1.0E0
  DO  j = 1, n
    z(j) = 0.0E0
  END DO

  DO  k = 1, n
    IF (z(k) /= 0.0E0) ek = SIGN(ek, -z(k))
    IF (ABS(ek-z(k)) > ABS(a(k,k))) THEN
      s = ABS(a(k,k)) / ABS(ek-z(k))
      CALL sscal(n, s, z, 1)
      ek = s*ek
    END IF
    wk = ek - z(k)
    wkm = -ek - z(k)
    s = ABS(wk)
    sm = ABS(wkm)
    IF (a(k,k) /= 0.0E0) THEN
      wk  = wk  / a(k,k)
      wkm = wkm / a(k,k)
    ELSE
      wk  = 1.0E0
      wkm = 1.0E0
    END IF
    kp1 = k + 1
    IF (kp1 <= n) THEN
      DO  j = kp1, n
        sm = sm + ABS(z(j)+wkm*a(k,j))
        z(j) = z(j) + wk*a(k,j)
        s = s + ABS(z(j))
      END DO
      IF (s < sm) THEN
        t = wkm - wk
        wk = wkm
        DO  j = kp1, n
          z(j) = z(j) + t*a(k,j)
        END DO
      END IF
    END IF
    z(k) = wk
  END DO

  s = 1.0E0 / sasum(n, z, 1)
  CALL sscal(n, s, z, 1)
  !                                ** SOLVE TRANS(L)*Y = W
  DO  kb = 1, n
    k = n + 1 - kb
    IF (k < n) z(k) = z(k) + sdot(n-k, a(k+1,k), 1, z(k+1), 1)
    IF (ABS(z(k)) > 1.0E0) THEN
      s = 1.0E0/ABS(z(k))
      CALL sscal(n, s, z, 1)
    END IF
    l = ipvt(k)
    t = z(l)
    z(l) = z(k)
    z(k) = t
  END DO

  s = 1.0E0 / sasum(n, z, 1)
  CALL sscal(n, s, z, 1)
  !                                 ** SOLVE L*V = Y
  ynorm = 1.0E0
  DO  k = 1, n
    l = ipvt(k)
    t = z(l)
    z(l) = z(k)
    z(k) = t
    IF (k < n) CALL saxpy(n-k, t, a(k+1,k), 1, z(k+1), 1)
    IF (ABS(z(k)) > 1.0E0) THEN
      s = 1.0E0/ABS(z(k))
      CALL sscal(n, s, z, 1)
      ynorm = s*ynorm
    END IF
  END DO

  s = 1.0E0 / sasum(n, z, 1)
  CALL sscal(n, s, z, 1)
  !                                  ** SOLVE  U*Z = V
  ynorm = s*ynorm
  DO  kb = 1, n
    k = n + 1 - kb
    IF (ABS(z(k)) > ABS(a(k,k))) THEN
      s = ABS(a(k,k))/ABS(z(k))
      CALL sscal(n, s, z, 1)
      ynorm = s*ynorm
    END IF
    IF (a(k,k) /= 0.0E0) z(k) = z(k)/a(k,k)
    IF (a(k,k) == 0.0E0) z(k) = 1.0E0
    t = -z(k)
    CALL saxpy(k-1, t, a(1,k), 1, z(1), 1)
  END DO
  !                                   ** MAKE ZNORM = 1.0
  s = 1.0E0 / sasum(n, z, 1)
  CALL sscal(n, s, z, 1)
  ynorm = s*ynorm

  IF (anorm /= 0.0E0) rcond = ynorm/anorm
  IF (anorm == 0.0E0) rcond = 0.0E0

  END SUBROUTINE  sgeco

  SUBROUTINE  sgefa( a, lda, n, ipvt, info )

  !         FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.

  !         REVISION DATE:  8/1/82
  !         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

  !     SGEFA IS USUALLY CALLED BY SGECO, BUT IT CAN BE CALLED
  !     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
  !     (TIME FOR SGECO) = (1 + 9/N)*(TIME FOR SGEFA) .

  !     INPUT:  SAME AS 'SGECO'

  !     ON RETURN:

  !        A,IPVT  SAME AS 'SGECO'

  !        INFO    INTEGER
  !                = 0  NORMAL VALUE.
  !                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
  !                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
  !                     INDICATE THAT SGESL OR SGEDI WILL DIVIDE BY ZERO
  !                     IF CALLED.  USE  RCOND  IN SGECO FOR A RELIABLE
  !                     INDICATION OF SINGULARITY.

  !     ROUTINES CALLED:  FROM BLAS:    SAXPY, SSCAL, ISAMAX
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: lda
  INTEGER, INTENT(IN)  :: n
  REAL, INTENT(IN OUT) :: a(lda,*)
  INTEGER, INTENT(OUT) :: ipvt(*)
  INTEGER, INTENT(OUT) :: info

  REAL :: t
  INTEGER :: j,k,kp1,l,nm1


  !                      ** GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
  info = 0
  nm1 = n - 1
  DO  k = 1, nm1
    kp1 = k + 1
  !                                            ** FIND L = PIVOT INDEX
    l = isamax( n-k+1, a(k,k), 1) + k-1
    ipvt(k) = l

    IF (a(l,k) == 0.0E0) THEN
  !                                     ** ZERO PIVOT IMPLIES THIS COLUMN
  !                                     ** ALREADY TRIANGULARIZED
      info = k
    ELSE
  !                                     ** INTERCHANGE IF NECESSARY
      IF (l /= k) THEN
        t = a(l,k)
        a(l,k) = a(k,k)
        a(k,k) = t
      END IF
  !                                     ** COMPUTE MULTIPLIERS
      t = -1.0E0 / a(k,k)
      CALL sscal( n-k, t, a(k+1,k), 1 )

  !                              ** ROW ELIMINATION WITH COLUMN INDEXING
      DO  j = kp1, n
        t = a(l,j)
        IF (l /= k) THEN
          a(l,j) = a(k,j)
          a(k,j) = t
        END IF
        CALL saxpy( n-k, t, a(k+1,k), 1, a(k+1,j), 1 )
      END DO

    END IF

  END DO

  ipvt(n) = n
  IF (a(n,n) == 0.0E0) info = n

  END SUBROUTINE  sgefa

  SUBROUTINE  sgesl( a, lda, n,ipvt, b, job )

  !         SOLVES THE REAL SYSTEM
  !            A * X = B  OR  TRANS(A) * X = B
  !         USING THE FACTORS COMPUTED BY SGECO OR SGEFA.

  !         REVISION DATE:  8/1/82
  !         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

  !     ON ENTRY

  !        A       REAL(LDA, N)
  !                THE OUTPUT FROM SGECO OR SGEFA.

  !        LDA     INTEGER
  !                THE LEADING DIMENSION OF THE ARRAY  A .

  !        N       INTEGER
  !                THE ORDER OF THE MATRIX  A .

  !        IPVT    INTEGER(N)
  !                THE PIVOT VECTOR FROM SGECO OR SGEFA.

  !        B       REAL(N)
  !                THE RIGHT HAND SIDE VECTOR.

  !        JOB     INTEGER
  !                = 0         TO SOLVE  A*X = B ,
  !                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
  !                            TRANS(A)  IS THE TRANSPOSE.

  !     ON RETURN

  !        B       THE SOLUTION VECTOR  X .

  !     ERROR CONDITION

  !        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
  !        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
  !        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
  !        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
  !        CALLED CORRECTLY AND IF SGECO HAS SET RCOND .GT. 0.0
  !        OR SGEFA HAS SET INFO .EQ. 0 .

  !     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
  !     WITH  P  COLUMNS
  !           CALL SGECO(A,LDA,N,IPVT,RCOND,Z)
  !           IF (RCOND IS TOO SMALL) GO TO ...
  !           DO 10 J = 1, P
  !              CALL SGESL(A,LDA,N,IPVT,C(1,J),0)
  !        10 CONTINUE


  !     ROUTINES CALLED:  FROM BLAS:    SAXPY, SDOT
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: lda
  INTEGER, INTENT(IN)  :: n
  INTEGER, INTENT(IN)  :: job
  INTEGER, INTENT(IN)  :: ipvt(*)
  REAL, INTENT(IN OUT) :: b(*)
  REAL, INTENT(IN)     :: a(lda,*)

  REAL :: t
  INTEGER :: k,kb,l,nm1


  nm1 = n - 1
  IF (job == 0) THEN
  !                                 ** JOB = 0 , SOLVE  A * X = B
  !                                     ** FIRST SOLVE  L*Y = B
    DO  k = 1, nm1
      l = ipvt(k)
      t = b(l)
      IF (l /= k) THEN
        b(l) = b(k)
        b(k) = t
      END IF
      CALL saxpy( n-k, t, a(k+1,k), 1, b(k+1), 1 )
    END DO
  !                                    ** NOW SOLVE  U*X = Y
    DO  kb = 1, n
      k = n + 1 - kb
      b(k) = b(k) / a(k,k)
      t = -b(k)
      CALL saxpy( k-1, t, a(1,k), 1, b(1), 1 )
    END DO

  ELSE
  !                         ** JOB = NONZERO, SOLVE  TRANS(A) * X = B
  !                                    ** FIRST SOLVE  TRANS(U)*Y = B
    DO  k = 1, n
      t = sdot( k-1, a(1,k), 1, b(1), 1 )
      b(k) = (b(k) - t) / a(k,k)
    END DO
  !                                    ** NOW SOLVE  TRANS(L)*X = Y
    DO  kb = 1, nm1
      k = n - kb
      b(k) = b(k) + sdot( n-k, a(k+1,k), 1, b(k+1), 1 )
      l = ipvt(k)
      IF (l /= k) THEN
        t = b(l)
        b(l) = b(k)
        b(k) = t
      END IF
    END DO

  END IF

  END SUBROUTINE  sgesl

   REAL FUNCTION  sasum( n, sx, incx )

  !  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR TO BE SUMMED
  !            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
  !          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'

  ! --OUTPUT-- SASUM   SUM FROM 0 TO N-1 OF  ABS(SX(1+I*INCX))
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(IN) :: incx
  REAL, INTENT(IN)    :: sx(*)

  INTEGER :: i,m

  sasum = 0.0
  IF( n <= 0 )  RETURN
  IF( incx /= 1 ) THEN
  !                                          ** NON-UNIT INCREMENTS
    DO  i = 1, 1+(n-1)*incx, incx
      sasum = sasum + ABS(sx(i))
    END DO
  ELSE
  !                                          ** UNIT INCREMENTS
    m = MOD(n,6)
    IF( m /= 0 ) THEN
  !                             ** CLEAN-UP LOOP SO REMAINING VECTOR
  !                             ** LENGTH IS A MULTIPLE OF 6.
      DO   i = 1, m
        sasum = sasum + ABS(sx(i))
      END DO
    END IF
  !                              ** UNROLL LOOP FOR SPEED
    DO   i = m+1, n, 6
      sasum = sasum + ABS(sx(i))   + ABS(sx(i+1)) + ABS(sx(i+2))  &
          + ABS(sx(i+3)) + ABS(sx(i+4)) + ABS(sx(i+5))
    END DO
  END IF

  END FUNCTION  sasum

  SUBROUTINE     saxpy( n, sa, sx, incx, sy, incy )

  !          Y = A*X + Y  (X, Y = VECTORS, A = SCALAR)

  !  --INPUT--
  !        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
  !       SA  SINGLE PRECISION SCALAR MULTIPLIER 'A'
  !       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
  !     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
  !       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
  !     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'

  ! --OUTPUT--
  !       SY   FOR I = 0 TO N-1, OVERWRITE  SY(LY+I*INCY) WITH
  !                 SA*SX(LX+I*INCX) + SY(LY+I*INCY),
  !            WHERE LX = 1          IF INCX .GE. 0,
  !                     = (-INCX)*N  IF INCX .LT. 0
  !            AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: n
  INTEGER, INTENT(IN)  :: incx
  INTEGER, INTENT(IN)  :: incy
  REAL, INTENT(IN)     :: sa
  REAL, INTENT(IN)     :: sx(*)
  REAL, INTENT(INOUT)    :: sy(*)

  INTEGER :: i,m,ix,iy

  IF( n <= 0 .OR. sa == 0.0 ) RETURN

  IF ( incx == incy .AND. incx > 1 )  THEN

    DO   i = 1, 1+(n-1)*incx, incx
      sy(i) = sy(i) + sa * sx(i)
    END DO

  ELSE IF ( incx == incy .AND. incx == 1 )  THEN

  !                                        ** EQUAL, UNIT INCREMENTS
    m = MOD(n,4)
    IF( m /= 0 ) THEN
  !                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
  !                            ** IS A MULTIPLE OF 4.
      DO   i = 1, m
        sy(i) = sy(i) + sa * sx(i)
      END DO
    END IF
  !                              ** UNROLL LOOP FOR SPEED
    DO   i = m+1, n, 4
      sy(i)   = sy(i)   + sa * sx(i)
      sy(i+1) = sy(i+1) + sa * sx(i+1)
      sy(i+2) = sy(i+2) + sa * sx(i+2)
      sy(i+3) = sy(i+3) + sa * sx(i+3)
    END DO

  ELSE
  !               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
    ix = 1
    iy = 1
    IF( incx < 0 )  ix = 1 + (n-1)*(-incx)
    IF( incy < 0 )  iy = 1 + (n-1)*(-incy)
    DO   i = 1, n
      sy(iy) = sy(iy) + sa*sx(ix)
      ix = ix + incx
      iy = iy + incy
    END DO

  END IF

  END SUBROUTINE     saxpy

   REAL FUNCTION  sdot( n, sx, incx, sy, incy )

  !          S.P. DOT PRODUCT OF VECTORS  'X'  AND  'Y'

  !  --INPUT--
  !        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
  !       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
  !     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
  !       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
  !     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'

  ! --OUTPUT--
  !     SDOT   SUM FOR I = 0 TO N-1 OF  SX(LX+I*INCX) * SY(LY+I*INCY),
  !            WHERE  LX = 1          IF INCX .GE. 0,
  !                      = (-INCX)*N  IF INCX .LT. 0,
  !            AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(IN) :: incx
  INTEGER, INTENT(IN) :: incy
  REAL, INTENT(IN)    :: sx(*)
  REAL, INTENT(IN)    :: sy(*)

  INTEGER :: i,m,ix,iy

  sdot = 0.0
  IF( n <= 0 )  RETURN

  IF ( incx == incy .AND. incx > 1 )  THEN

    DO   i = 1, 1+(n-1)*incx, incx
      sdot = sdot + sx(i) * sy(i)
    END DO

  ELSE IF ( incx == incy .AND. incx == 1 )  THEN

  !                                        ** EQUAL, UNIT INCREMENTS
    m = MOD(n,5)
    IF( m /= 0 ) THEN
  !                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
  !                            ** IS A MULTIPLE OF 4.
      DO   i = 1, m
        sdot = sdot + sx(i) * sy(i)
      END DO
    END IF
  !                              ** UNROLL LOOP FOR SPEED
    DO   i = m+1, n, 5
      sdot = sdot + sx(i)*sy(i)     + sx(i+1)*sy(i+1)  &
          + sx(i+2)*sy(i+2) + sx(i+3)*sy(i+3) + sx(i+4)*sy(i+4)
    END DO

  ELSE
  !               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
    ix = 1
    iy = 1
    IF( incx < 0 )  ix = 1 + (n-1)*(-incx)
    IF( incy < 0 )  iy = 1 + (n-1)*(-incy)
    DO   i = 1, n
      sdot = sdot + sx(ix) * sy(iy)
      ix = ix + incx
      iy = iy + incy
    END DO

  END IF

  END FUNCTION  sdot

  SUBROUTINE     sscal( n, sa, sx, incx )

  !         CALCULATE  X = A*X  (X = VECTOR, A = SCALAR)

  !  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR
  !            SA  SINGLE PRECISION SCALE FACTOR
  !            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
  !          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'

  ! --OUTPUT-- SX  REPLACE  SX(1+I*INCX)  WITH  SA * SX(1+I*INCX)
  !                FOR I = 0 TO N-1
  IMPLICIT NONE

  INTEGER, INTENT(IN)   :: n
  INTEGER, INTENT(IN)   :: incx
  REAL, INTENT(IN)      :: sa
  REAL, INTENT(INOUT)   :: sx(*)

  INTEGER :: i,m

  IF( n <= 0 ) RETURN

  IF( incx /= 1 ) THEN

    DO   i = 1, 1+(n-1)*incx, incx
      sx(i) = sa * sx(i)
    END DO

  ELSE

    m = MOD(n,5)
    IF( m /= 0 ) THEN
  !                           ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
  !                           ** IS A MULTIPLE OF 5.
      DO   i = 1, m
        sx(i) = sa * sx(i)
      END DO
    END IF
  !                             ** UNROLL LOOP FOR SPEED
    DO   i = m+1, n, 5
      sx(i)   = sa * sx(i)
      sx(i+1) = sa * sx(i+1)
      sx(i+2) = sa * sx(i+2)
      sx(i+3) = sa * sx(i+3)
      sx(i+4) = sa * sx(i+4)
    END DO

  END IF

  END SUBROUTINE     sscal

  SUBROUTINE     sswap( n, sx, incx, sy, incy )

  !          INTERCHANGE S.P VECTORS  X  AND  Y

  !  --INPUT--
  !        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
  !       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
  !     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
  !       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
  !     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'

  ! --OUTPUT--
  !       SX  INPUT VECTOR SY (UNCHANGED IF N .LE. 0)
  !       SY  INPUT VECTOR SX (UNCHANGED IF N .LE. 0)

  !     FOR I = 0 TO N-1, INTERCHANGE  SX(LX+I*INCX) AND SY(LY+I*INCY),
  !     WHERE LX = 1          IF INCX .GE. 0,
  !              = (-INCX)*N  IF INCX .LT. 0
  !     AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: n
  INTEGER, INTENT(IN)  :: incx
  INTEGER, INTENT(IN)  :: incy
  REAL, INTENT(IN OUT) :: sx(*)
  REAL, INTENT(IN OUT) :: sy(*)

  REAL :: stemp1, stemp2, stemp3
  INTEGER :: ix,iy,i,m

  IF( n <= 0 ) RETURN

  IF ( incx == incy .AND. incx > 1 )  THEN

    DO   i = 1, 1+(n-1)*incx, incx
      stemp1 = sx(i)
      sx(i) = sy(i)
      sy(i) = stemp1
    END DO

  ELSE IF ( incx == incy .AND. incx == 1 )  THEN

  !                                        ** EQUAL, UNIT INCREMENTS
    m = MOD(n,3)
    IF( m /= 0 ) THEN
  !                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
  !                            ** IS A MULTIPLE OF 3.
      DO   i = 1, m
        stemp1 = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp1
      END DO
    END IF
  !                              ** UNROLL LOOP FOR SPEED
    DO   i = m+1, n, 3
      stemp1  = sx(i)
      stemp2  = sx(i+1)
      stemp3  = sx(i+2)
      sx(i)   = sy(i)
      sx(i+1) = sy(i+1)
      sx(i+2) = sy(i+2)
      sy(i)   = stemp1
      sy(i+1) = stemp2
      sy(i+2) = stemp3
    END DO

  ELSE
  !               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
    ix = 1
    iy = 1
    IF( incx < 0 )  ix = 1 + (n-1)*(-incx)
    IF( incy < 0 )  iy = 1 + (n-1)*(-incy)
    DO   i = 1, n
      stemp1 = sx(ix)
      sx(ix) = sy(iy)
      sy(iy) = stemp1
      ix = ix + incx
      iy = iy + incy
    END DO

  END IF

  END SUBROUTINE     sswap

   INTEGER FUNCTION  isamax( n, sx, incx )

  !  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR OF INTEREST
  !            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
  !          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'

  ! --OUTPUT-- ISAMAX   FIRST I, I = 1 TO N, TO MAXIMIZE
  !                         ABS(SX(1+(I-1)*INCX))
  IMPLICIT NONE

  INTEGER, INTENT(IN)                      :: n
  INTEGER, INTENT(IN)                      :: incx
  REAL, INTENT(IN OUT)                     :: sx(*)

  REAL :: smax, xmag
  INTEGER :: iii, i

  IF( n <= 0 ) THEN
    isamax = 0
  ELSE IF( n == 1 ) THEN
    isamax = 1
  ELSE
    smax = 0.0
    iii = 1
    DO   i = 1, 1+(n-1)*incx, incx
      xmag = ABS(sx(i))
      IF( smax < xmag ) THEN
        smax = xmag
        isamax = iii
      END IF
      iii = iii + 1
    END DO
  END IF

  END FUNCTION  isamax

  DOUBLE PRECISION FUNCTION d1mach(i)

  !-----------------------------------------------------------------------------*
  != PURPOSE:                                                                  =*
  != D1MACH calculates various machine constants in single precision.          =*
  !-----------------------------------------------------------------------------*
  != PARAMETERS:                                                               =*
  !=   I       -  INTEGER, identifies the machine constant (0<I<5)         (I) =*
  !=   D1MACH  -  REAL, machine constant in single precision               (O) =*
  !=      I=1     - the smallest non-vanishing normalized floating-point       =*
  !=                power of the radix, i.e., D1MACH=FLOAT(IBETA)**MINEXP      =*
  !=      I=2     - the largest finite floating-point number.  In              =*
  !=                particular D1MACH=(1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP        =*
  !=                Note - on some machines D1MACH will be only the            =*
  !=                second, or perhaps third, largest number, being            =*
  !=                too small by 1 or 2 units in the last digit of             =*
  !=                the significand.                                           =*
  !=      I=3     - A small positive floating-point number such that           =*
  !=                1.0-D1MACH .NE. 1.0. In particular, if IBETA = 2           =*
  !=                or  IRND = 0, D1MACH = FLOAT(IBETA)**NEGEPS.               =*
  !=                Otherwise,  D1MACH = (IBETA**NEGEPS)/2.  Because           =*
  !=                NEGEPS is bounded below by -(IT+3), D1MACH may not         =*
  !=                be the smallest number that can alter 1.0 by               =*
  !=                subtraction.                                               =*
  !=      I=4     - the smallest positive floating-point number such           =*
  !=                that  1.0+D1MACH .NE. 1.0. In particular, if either        =*
  !=                IBETA = 2  or  IRND = 0, D1MACH=FLOAT(IBETA)**MACHEP.      =*
  !=                Otherwise, D1MACH=(FLOAT(IBETA)**MACHEP)/2                 =*
  !=  (see routine T665D for more information on different constants)          =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER, INTENT(IN)                      :: i

  select case(i)
  case(1)
  	d1mach = tiny(1.0d0)
  case(2)
  	d1mach = huge(1.0d0)
  case(3)
	d1mach = -epsilon(1.0d0)
  case(4)
 	d1mach = epsilon(1.0d0)
  case default
  	WRITE(0,*) '>>> ERROR (D1MACH) <<<  invalid argument'
    	STOP
  end select


  END FUNCTION d1mach


   REAL FUNCTION r1mach(i)

  !-----------------------------------------------------------------------------*
  != PURPOSE:                                                                  =*
  != R1MACH calculates various machine constants in single precision.          =*
  !-----------------------------------------------------------------------------*
  != PARAMETERS:                                                               =*
  !=   I       -  INTEGER, identifies the machine constant (0<I<5)         (I) =*
  !=   R1MACH  -  REAL, machine constant in single precision               (O) =*
  !=      I=1     - the smallest non-vanishing normalized floating-point       =*
  !=                power of the radix, i.e., R1MACH=FLOAT(IBETA)**MINEXP      =*
  !=      I=2     - the largest finite floating-point number.  In              =*
  !=                particular R1MACH=(1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP        =*
  !=                Note - on some machines R1MACH will be only the            =*
  !=                second, or perhaps third, largest number, being            =*
  !=                too small by 1 or 2 units in the last digit of             =*
  !=                the significand.                                           =*
  !=      I=3     - A small positive floating-point number such that           =*
  !=                1.0-R1MACH .NE. 1.0. In particular, if IBETA = 2           =*
  !=                or  IRND = 0, R1MACH = FLOAT(IBETA)**NEGEPS.               =*
  !=                Otherwise,  R1MACH = (IBETA**NEGEPS)/2.  Because           =*
  !=                NEGEPS is bounded below by -(IT+3), R1MACH may not         =*
  !=                be the smallest number that can alter 1.0 by               =*
  !=                subtraction.                                               =*
  !=      I=4     - the smallest positive floating-point number such           =*
  !=                that  1.0+R1MACH .NE. 1.0. In particular, if either        =*
  !=                IBETA = 2  or  IRND = 0, R1MACH=FLOAT(IBETA)**MACHEP.      =*
  !=                Otherwise, R1MACH=(FLOAT(IBETA)**MACHEP)/2                 =*
  !=  (see routine T665R for more information on different constants)          =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER, INTENT(IN)                      :: i

  select case(i)
  case(1)
  	r1mach = tiny(1.0)
  case(2)
  	r1mach = huge(1.0)
  case(3)
	r1mach = -epsilon(1.0)
  case(4)
 	r1mach = epsilon(1.0)
  case default
  	WRITE(0,*) '>>> ERROR (D1MACH) <<<  invalid argument'
    	STOP
  end select

  END FUNCTION r1mach

  SUBROUTINE swchem(nw,wl,nz,nj,nbl,tLev,airDen,sj)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Load various "weighting functions" (products of cross section and        =*
  !=  quantum yield at each altitude and each wavelength).  The altitude       =*
  !=  dependence is necessary to ensure the consideration of pressure and      =*
  !=  temperature dependence of the cross sections or quantum yields.          =*
  !=  The actual reading, evaluation and interpolation is done in separate     =*
  !=  subroutines for ease of management and manipulation.  Please refer to    =*
  !=  the inline documentation of the specific subroutines for detail          =*
  !=  information.                                                             =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section * quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !=  OPT    - If opt=1 read, otherwise just calc
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER, INTENT(IN)   :: nw
  INTEGER, INTENT(IN)   :: nbl
  INTEGER, INTENT(IN)   :: nz(nbl)
  REAL, INTENT(IN)      :: wl(kw)
  REAL, INTENT(IN)      :: tLev(nbl,kz)
  REAL, INTENT(IN)      :: airDen(nbl,kz)
  REAL, INTENT(INOUT)   :: sj(nbl,kj,kz,kw)
  INTEGER, INTENT(OUT)  :: nj

  REAL :: wc(kw)
  INTEGER :: iw

  INTEGER :: j

  ! complete wavelength grid
  DO   iw = 1, nw - 1
    wc(iw) = (wl(iw) + wl(iw+1))/2.
  END DO

  !____________________________________________________________________________

  ! O2 + hv -> O + O
  ! reserve first position.  Cross section parameterization in Schumman-Runge and
  ! Lyman-alpha regions are zenith-angle dependent, will be written in
  ! subroutine seto2.f.

    j = 1
    !jlabel(j) = 'O2 -> O + O'

    ! O3 + hv ->  (both channels)
    IF(doReaction(2) .or. doReaction(3)) CALL r01(nw,wl,wc,nz,j,nbl,sj,tLev) !1,2

    ! NO2 + hv -> NO + O(3P)
    IF(doReaction(4)) CALL r02(nw,j,nbl,sj,tLev) !3

    ! NO3 + hv ->  (both channels)
    IF(doReaction(5) .or. doReaction(6)) CALL r03(nw,wc,j,nbl,sj,tLev) !4,5

    ! N2O5 + hv -> (both channels)
    IF(doReaction(7) .or. doReaction(8)) CALL r04(nw,wc,j,nbl,sj,tLev) !6,7

    ! N2O + hv -> N2 + O(1D)
    IF(doReaction(9)) CALL r44(nw,wc,j,nbl,sj,tLev) !8

    ! HO2 + hv -> OH + O
    IF(doReaction(10)) CALL r39(nw,wc,j,nbl,sj,tLev) !9

    ! H2O2 + hv -> 2 OH
    IF(doReaction(11)) CALL r08(nw,wl,wc,j,nbl,sj,tLev) !10

    ! HNO2 + hv -> OH + NO
    IF(doReaction(12)) CALL r05(nw,j,nbl,sj,tLev) !11

    ! HNO3 + hv -> OH + NO2
    IF(doReaction(13)) CALL r06(nw,j,nbl,sj,tLev) !12

    ! HNO4 + hv -> HO2 + NO2
    IF(doReaction(14)) CALL r07(nw,j,nbl,sj,tLev) !13

    ! CH2O + hv -> (both channels)
    IF(doReaction(15) .or. doReaction(16)) &
                             CALL r10(nw,wl,wc,j,nbl,sj,tLev,airDen) !14,15

    ! CH3CHO + hv -> (all three channels)
    IF(doReaction(17) .or. doReaction(18) .or. doReaction(19)) &
                              CALL r11(nw,j,nbl,sj,tLev,airDen) !16,17,18

    ! C2H5CHO + hv -> C2H5 + HCO
    IF(doReaction(20)) CALL r12(nw,j,nbl,sj,tLev,airDen) !19

    ! CHOCHO + hv -> Products
    IF(doReaction(21) .or. doReaction(22)) CALL r13(nw,wc,j,nbl,sj,tLev) !20,21

    ! CH3COCHO + hv -> Products
    IF(doReaction(23)) CALL r14(nw,wc,j,nbl,sj,tLev,airDen) !22

    ! CH3COCH3 + hv -> Products
    IF(doReaction(24)) CALL r15(nw,wc,j,nbl,sj,tLev,airDen) !23

    ! CH3OOH + hv -> CH3O + OH
    IF(doReaction(25)) CALL r16(nw,j,nbl,sj,tLev) !24

    ! CH3ONO2 + hv -> CH3O + NO2
    IF(doReaction(26)) CALL r17(nw,j,nbl,sj,tLev) !25

    ! PAN + hv -> Products
    IF(doReaction(27)) CALL r18(nw,j,nbl,sj,tLev) !26

    ! ClOO + hv -> Products
    IF(doReaction(28)) CALL r31(nw,j,nbl,sj,tLev) !27

    ! ClONO2 + hv -> Products
    IF(doReaction(29) .or. doReaction(30)) CALL r45(nw,wc,j,nbl,sj,tLev) !28,29

    ! CH3Cl + hv -> Products
    IF(doReaction(31)) CALL r30(nw,j,nbl,sj,tLev) !30

    ! CCl2O + hv -> Products
    IF(doReaction(32)) CALL r19(nw,j,nbl,sj,tLev) !31

    ! CCl4 + hv -> Products
    IF(doReaction(33)) CALL r20(nw,j,nbl,sj,tLev) !32

    ! CClFO + hv -> Products
    IF(doReaction(34)) CALL r21(nw,j,nbl,sj,tLev) !33

    ! CCF2O + hv -> Products
    IF(doReaction(35)) CALL r22(nw,j,nbl,sj,tLev) !34

    ! CF2ClCFCl2 (CFC-113) + hv -> Products
    IF(doReaction(36)) CALL r23(nw,j,nbl,sj,tLev) !35

    ! CF2ClCF2Cl (CFC-114) + hv -> Products
    IF(doReaction(37)) CALL r24(nw,j,nbl,sj,tLev) !36

    ! CF3CF2Cl (CFC-115) + hv -> Products
    IF(doReaction(38)) CALL r25(nw,j,nbl,sj,tLev) !37

    ! CCl3F (CFC-111) + hv -> Products
    IF(doReaction(39)) CALL r26(nw,wc,j,nbl,sj,tLev) !38

    ! CCl2F2 (CFC-112) + hv -> Products
    IF(doReaction(40)) CALL r27(nw,j,nbl,sj,tLev) !39

    ! CH3CCl3 + hv -> Products
    IF(doReaction(41)) CALL r29(nw,j,nbl,sj,tLev) !40

    ! CF3CHCl2 (HCFC-123) + hv -> Products
    IF(doReaction(42)) CALL r32(nw,wc,j,nbl,sj,tLev) !41

    ! CF3CHFCl (HCFC-124) + hv -> Products
    IF(doReaction(43)) CALL r33(nw,wc,j,nbl,sj,tLev) !42

    ! CH3CFCl2 (HCFC-141b) + hv -> Products
    IF(doReaction(44)) CALL r34(nw,j,nbl,sj,tLev) !43

    ! CH3CF2Cl (HCFC-142b) + hv -> Products
    IF(doReaction(45)) CALL r35(nw,wc,j,nbl,sj,tLev) !44

    ! CF3CF2CHCl2 (HCFC-225ca) + hv -> Products
    IF(doReaction(46)) CALL r36(nw,j,nbl,sj,tLev) !45

    ! CF2ClCF2CHFCl (HCFC-225cb) + hv -> Products
    IF(doReaction(47)) CALL r37(nw,j,nbl,sj,tLev) !46

    ! CHClF2 (HCFC-22) + hv -> Products
    IF(doReaction(48)) CALL r38(nw,j,nbl,sj,tLev) !47

    ! BrONO2 + hv -> Products
    IF(doReaction(49) .or. doReaction(50)) CALL r46(nw,j,nbl,sj,tLev) !48,49

    ! CH3Br + hv -> Products
    IF(doReaction(51)) CALL r28(nw,j,nbl,sj,tLev) !50

    ! CHBr3 + hv -> Products
    IF(doReaction(52)) CALL r09(nw,wl,wc,j,nbl,sj,tLev) !51

    ! CF3Br (Halon-1301) + hv -> Products
    IF(doReaction(53)) CALL r42(nw,j,nbl,sj,tLev) !52

    ! CF2BrCF2Br (Halon-2402) + hv -> Products
    IF(doReaction(54)) CALL r43(nw,j,nbl,sj,tLev) !53

    ! CF2Br2 (Halon-1202) + hv -> Products
    IF(doReaction(55)) CALL r40(nw,j,nbl,sj,tLev) !54

    ! CF2BrCl (Halon-1211) + hv -> Products
    IF(doReaction(56)) CALL r41(nw,j,nbl,sj,tLev) !55

    ! CL2 + hc -> CL + CL
    IF(doReaction(57)) CALL r47(nw,j,nbl,sj,tLev) !56

    ! CH2(OH)CH + hv -> Products
    IF(doReaction(58)) CALL r101(nw,j,nbl,sj,tLev) !57

    ! CH3COCOCH3 + hv -> Products
    IF(doReaction(59)) CALL r102(nw,j,nbl,sj,tLev) !58

    ! CH3COCHCH2 + hv -> Products
    IF(doReaction(60)) CALL r103(nw,wc,j,nbl,sj,tLev,airDen) !59

    ! CH2C(CH3)CHO + hv -> Products
    IF(doReaction(61)) CALL r104(nw,j,nbl,sj,tLev) !60

    ! CH3COCO(OH) + hv -> Products
    IF(doReaction(62)) CALL r105(nw,j,nbl,sj,tLev) !61

    ! CH3CH2ONO2 -> CH3CH2O + NO2
    IF(doReaction(63)) CALL r106(nw,j,nbl,sj,tLev) !62

    ! CH3CHONO2CH3 -> CH3CHOCH3 + NO2
    IF(doReaction(64)) CALL r107(nw,j,nbl,sj,tLev) !63

    ! CH2(OH)CH2(ONO2) -> CH2(OH)CH2(O.) + NO2
    IF(doReaction(65)) CALL r108(nw,wc,j,nbl,sj,tLev) !64

    ! CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2
    IF(doReaction(66)) CALL r109(nw,wc,j,nbl,sj,tLev) !65

    ! C(CH3)3(ONO2) -> C(CH3)3(O.) + NO2
    IF(doReaction(67)) CALL r110(nw,wc,j,nbl,sj,tLev) !66

    ! ClOOCl -> Cl + ClOO
    IF(doReaction(68)) CALL r111(nw,j,nbl,sj,tLev) !67

    ! CH2(OH)COCH3 -> CH3CO + CH2(OH)
    ! CH2(OH)COCH3 -> CH2(OH)CO + CH3
    IF(doReaction(69) .or. doReaction(70)) CALL r112(nw,j,nbl,sj,tLev) !68,69

    ! HOBr -> OH + Br'
    IF(doReaction(71)) CALL r113(nw,wc,j,nbl,sj,tLev) !70

    ! BrO -> Br + O'
    IF(doReaction(72)) CALL r114(nw,j,nbl,sj,tLev) !71

    ! Br2 -> Br + Br'
    IF(doReaction(73)) CALL r115(nw,j,nbl,sj,tLev) !72

    !***************************************************************


    IF (j > kj) STOP '1002'

    nj=j

  END SUBROUTINE swchem

  SUBROUTINE r01(nw,wl,wc,nz,j,nbl,sj,tLev)
  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide the product of (cross section) x (quantum yield) for the two     =*
  !=  O3 photolysis reactions:                                                 =*
  !=             (a) O3 + hv -> O2 + O(1D)                                     =*
  !=             (b) O3 + hv -> O2 + O(3P)                                     =*
  !=  Cross section:  Combined data from WMO 85 Ozone Assessment (use 273K     =*
  !=                  value from 175.439-847.5 nm) and data from Molina and    =*
  !=                  Molina (use in Hartley and Huggins bans (240.5-350 nm)   =*
  !=  Quantum yield:  Choice between                                           =*
  !=                   (1) data from Michelsen et al, 1994                     =*
  !=                   (2) JPL 87 recommendation                               =*
  !=                   (3) JPL 90/92 recommendation (no "tail")                =*
  !=                   (4) data from Shetter et al., 1996                      =*
  !=                   (5) JPL 97 recommendation                               =*
  !=                   (6) JPL 00 recommendation                               =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !=  OPT    - If opt=1 read, otherwise just calc
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=01 !1,2

  INTEGER, INTENT(IN)             :: nw
  INTEGER, INTENT(IN)             :: nbl
  INTEGER, INTENT(IN)             :: nz(nbl)
  REAL, INTENT(IN)                :: wl(kw)
  REAL, INTENT(IN)                :: wc(kw)
  REAL, INTENT(IN)                :: tLev(nbl,kz)
  INTEGER, INTENT(INOUT)          :: j
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob
  INTEGER, PARAMETER :: kdata = 500
  REAL :: xs(nbl,kz,kw)
  REAL :: qy1d(nbl)
  REAL :: qy3p
  REAL :: tau, tau2, tau3
  REAL :: a, b, c
  REAL :: a0, a1, a2, a3, a4, a5, a6 !, a7
  REAL :: xl, xl0
  INTEGER :: i, iw

  j = j + 1
  !jlabel(j) = 'O3 -> O2 + O(1D)'
  j = j + 1
  !jlabel(j) = 'O3 -> O2 + O(3P)'

  SELECT CASE (mOption(1))
     CASE (1)
        DO iob=1,nbl
           CALL o3xs_mm(nw,wl,nz(iob),xs(iob,:,:),&
	             tLev(iob,:))
        END DO
     CASE (2)
        DO iob=1,nbl
           CALL o3xs_mal(nw,wl,nz(iob),xs(iob,:,:),&
	             tLev(iob,:))
        END DO
     CASE (3)
        DO iob=1,nbl
           CALL o3xs_bass(nw,wl,nz(iob),xs(iob,:,:),&
	             tLev(iob,:))
        END DO
  END SELECT

  ! compute cross sections and yields at different wavelengths, altitudes:

  ! quantum yields
  DO  iw = 1, nw-1
    DO  i = 1, zLimit
      ! coefficients from jpl 87:
      IF (mOption(2) == kjpl87) THEN
        DO iob=1,nbl
           tau = tLev(iob,i) - 230.
           tau2 = tau*tau
           tau3 = tau2*tau
           xl = wc(iw)
           xl0 = 308.2 + 4.4871E-2*tau + 6.938E-5*tau2 - 2.5452E-6*tau3
           a = 0.9*(0.369 + 2.85E-4*tau + 1.28E-5*tau2 + 2.57E-8*tau3)
           b     = -0.575 + 5.59E-3*tau - 1.439E-5*tau2 - 3.27E-8*tau3
           c = 0.9*(0.518 + 9.87E-4*tau - 3.94E-5*tau2 + 3.91E-7*tau3)
           qy1d(iob) = a*ATAN(b*(xl-xl0)) + c
           qy1d(iob) = AMAX1(0.,qy1d(iob))
           qy1d(iob) = AMIN1(0.9,qy1d(iob))
        END DO
      END IF

      ! from jpl90, jpl92:
      ! (caution: error in JPL92 for first term of a3)
      IF (mOption(2) == kjpl92) THEN
        DO iob=1,nbl
           tau = 298. - tLev(iob,i)
           tau2 = tau*tau
           xl0 = wc(iw) - 305.
           a0 =   .94932   - 1.7039E-4*tau + 1.4072E-6*tau2
           a1 = -2.4052E-2 + 1.0479E-3*tau - 1.0655E-5*tau2
           a2 =  1.8771E-2 - 3.6401E-4*tau - 1.8587E-5*tau2
           a3 = -1.4540E-2 - 4.7787E-5*tau + 8.1277E-6*tau2
           a4 =  2.3287E-3 + 1.9891E-5*tau - 1.1801E-6*tau2
           a5 = -1.4471E-4 - 1.7188E-6*tau + 7.2661E-8*tau2
           a6 =  3.1830E-6 + 4.6209E-8*tau - 1.6266E-9*tau2
           qy1d(iob) = a0 + a1*xl0 + a2*(xl0)**2 + a3*(xl0)**3 +  &
  	      a4*(xl0)**4 + a5*(xl0)**5 + a6*(xl0)**6
           IF (wc(iw) < 305.) qy1d(iob) = 0.95
           IF (wc(iw) > 320.) qy1d(iob) = 0.
           IF (qy1d(iob) < 0.02) qy1d(iob) = 0.
        END DO
      END IF

      ! from JPL'97
      IF (mOption(2) == kjpl97) THEN
        IF (wc(iw) < 271.) THEN
           DO iob=1,nbl
  	      qy1d(iob) = 0.87
           END DO
        ELSE IF (wc(iw) >= 271. .AND. wc(iw) < 290.) THEN
           DO iob=1,nbl
	      qy1d(iob) = 0.87 + (wc(iw)-271.)*(.95-.87)/(290.-271.)
           END DO
	ELSE IF (wc(iw) >= 290. .AND. wc(iw) < 305.) THEN
           DO iob=1,nbl
  	      qy1d(iob) = 0.95
           END DO
        ELSE IF (wc(iw) >= 305. .AND. wc(iw) <= 325.) THEN
           DO iob=1,nbl
	      IF(i>nz(iob)) EXIT
  	      qy1d(iob) = yg1(iw,nReact) * EXP ( -yg2(iw,nReact) / &
	              tLev(iob,i) )
	   END DO
        ELSE
           DO iob=1,nbl
  	      qy1d(iob) = 0.
           END DO
        END IF
      END IF

      ! from Michelsen, H. A., R.J. Salawitch, P. O. Wennber, and J. G. Anderson
      ! Geophys. Res. Lett., 21, 2227-2230, 1994.
      IF (mOption(2) == kmich) THEN
        IF (wc(iw) < 271.) THEN
           DO iob=1,nbl
  	      qy1d(iob) = 0.87
	   END DO
        ELSE IF (wc(iw) >= 271. .AND. wc(iw) < 305.) THEN
           DO iob=1,nbl
  	       qy1d(iob) = 1.98 - 301./wc(iw)
	   END DO
        ELSE IF (wc(iw) >= 305. .AND. wc(iw) <= 325.) THEN
           DO iob=1,nbl
	      IF(i>nz(iob)) EXIT
  	      qy1d(iob) = yg1(iw,nReact) * EXP (-yg2(iw,nReact) / &
	          (0.6951*tLev(iob,i)))
           END DO
        ELSE
           DO iob=1,nbl
  	      qy1d(iob) = 0.
	   END DO
        END IF
      END IF

      ! Shetter et al.:
      ! phi = A * exp(-B/T), A and B are based on meas. at 298 and 230 K
      ! do linear interpolation between phi(298) and phi(230) for wavelengths > 321
      ! as phi(230)=0. for those wavelengths, so there are no A and B factors
      IF (mOption(2) == kshet) THEN
        IF (wl(iw+1) <= 321.) THEN
           DO iob=1,nbl
  	       qy1d(iob) = yg1(iw,nReact) * EXP(-1. * yg2(iw,nReact)/ &
	           tLev(iob,i))
           END DO
        ELSE
           DO iob=1,nbl
  	       qy1d(iob) = (yg3(iw,nReact) - yg4(iw,nReact))/(298.-230.) * &
	           (tLev(iob,i)-230.) + yg4(iw,nReact)
	   END DO
        END IF
      END IF

      ! JPL 2000:
      IF (mOption(2) == kjpl00) THEN
         DO iob=1,nbl
            qy1d(iob) = fo3qy(wc(iw),tLev(iob,i))
	 END DO
      END IF

      ! Matsumi et al.
      IF (mOption(2) == kmats) THEN
         DO iob=1,nbl
            qy1d(iob) = fo3qy2(wc(iw),tLev(iob,i))
         END DO
      END IF

      ! compute product
      DO iob=1,nbl
         sj(iob,2,i,iw) = qy1d(iob)*xs(iob,i,iw) !1
         qy3p = 1.0 - qy1d(iob)
         sj(iob,3,i,iw) = qy3p*xs(iob,i,iw) !2
      END DO
    END DO
  END DO

  END SUBROUTINE r01

  !=============================================================================*
  SUBROUTINE r02(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide the product (cross section) x (quantum yield) for NO2            =*
  !=  photolysis:                                                              =*
  !=         NO2 + hv -> NO + O(3P)                                            =*
  !=  Cross section from JPL94 (can also have Davidson et al.)                 =*
  !=  Quantum yield from Gardiner, Sperry, and Calvert                         =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !=  OPT    - If opt=1 read, otherwise just calc
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=02 !3

  INTEGER, INTENT(IN)             :: nw
  INTEGER, INTENT(IN)             :: nbl
  INTEGER, INTENT(INOUT)          :: j
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob
  REAL, INTENT(IN) :: tLev(nbl,kz)

  INTEGER, PARAMETER :: kdata=200
  INTEGER :: i, iw
  REAL :: r2_no2xs(nbl,kz,kw)

  !*************** NO2 photodissociation
  j = j + 1
  !jlabel(j) = 'NO2 -> NO + O(3P)'

  SELECT CASE(mOption(3))
     CASE (1)
       DO iw = 1, nw-1
         DO i = 1, zLimit
            DO iob=1,nbl
               r2_no2xs(nbl,i,iw) = yg1(iw,nReact) + &
	          yg2(iw,nReact)*(tLev(iob,i)-273.15)
            END DO
	 END DO
       END DO
     CASE (2)
       DO iw = 1, nw-1
         DO i = 1, zLimit
            DO iob=1,nbl
               r2_no2xs(nbl,i,iw) = yg1(iw,nReact)
	    END DO
         END DO
       END DO
  END SELECT

  ! combine
  DO iw = 1, nw - 1
    DO i = 1, zLimit
      DO iob=1,nbl
         sj(iob,4,i,iw) = r2_no2xs(nbl,i,iw)* &
	       yg1n(iw,nReact)	 !3
      END DO
    END DO
  END DO

  END SUBROUTINE r02

  !=============================================================================*

  SUBROUTINE r03(nw,wc,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide the product (absorptioon cross section) x (quantum yield) for    =*
  !=  both channels of NO3 photolysis:                                         =*
  !=          (a) NO3 + hv -> NO2 + O(3P)                                      =*
  !=          (b) NO3 + hv -> NO + O2                                          =*
  !=  Cross section combined from Graham and Johnston (<600 nm) and JPL 94     =*
  !=  Quantum yield from Madronich (1988)                                      =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !=  OPT    - If opt=1 read, otherwise just calc
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=03 !4,5

  INTEGER, INTENT(IN)             :: nw
  INTEGER, INTENT(IN)             :: nbl
  REAL, INTENT(IN)                :: wc(kw)
  INTEGER, INTENT(INOUT)          :: j
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob
  REAL, INTENT(IN) :: tLev(nbl,kz)

  INTEGER, PARAMETER :: kdata=350
  REAL :: qyVal(nbl)
  INTEGER :: i, iw

  !***************      jlabel(j) = 'NO3 -> NO2 + O(3P)'
  !***************      jlabel(j) = 'NO3 -> NO + O2'
   ! quantum yield:
   ! from Madronich (1988) see CEC NO3 book.
      j = j + 1
   !jlabel(j) = 'NO3 -> NO + O2'

   ! for   NO3 ->NO+O2
   DO iw = 1, nw - 1
     IF (wc(iw) < 584.) THEN
        DO iob=1,nbl
           qyVal(iob) = 0.
        END DO
     ELSE IF (wc(iw) >= 640.) THEN
        DO iob=1,nbl
           qyVal(iob) = 0.
        END DO
     ELSE IF (wc(iw) >= 595.) THEN
        DO iob=1,nbl
           qyVal(iob) = 0.35*(1.-(wc(iw)-595.)/45.)
        END DO
     ELSE
        DO iob=1,nbl
           qyVal(iob) = 0.35*(wc(iw)-584.)/11.
	END DO
     END IF
     DO i = 1, zLimit
        DO iob=1,nbl
           sj(iob,5,i,iw) = yg(iw,nReact)*qyVal(iob) !4
	END DO
     END DO
         !PRINT *,'LFR->',sj(iob,j,1,iw),yg(iw,nReact),qyVal(iob)
   END DO

   ! for  NO3 ->NO2+O
      j = j + 1

   !jlabel(j) = 'NO3 -> NO2 + O(3P)'
   ! for  NO3 ->NO2+O
   DO iw = 1, nw - 1
     IF (wc(iw) < 584.) THEN
        DO iob=1,nbl
           qyVal(iob) = 1.
	END DO
     ELSE IF (wc(iw) > 640.) THEN
        DO iob=1,nbl
           qyVal(iob) = 0.
	END DO
     ELSE IF (wc(iw) > 595.) THEN
        DO iob=1,nbl
           qyVal(iob) = 0.65*(1-(wc(iw)-595.)/45.)
	END DO
     ELSE
        DO iob=1,nbl
          qyVal(iob) = 1.-0.35*(wc(iw)-584.)/11.
	END DO
     END IF
     DO i = 1, zLimit
       DO iob=1,nbl
          sj(iob,6,i,iw) = yg(iw,nReact)*qyVal(iob) !5
       END DO
     END DO
   END DO


  END SUBROUTINE r03

  !=============================================================================*

  SUBROUTINE r04(nw,wc,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product of (cross section) x (quantum yiels) for N2O5 photolysis =*
  !=  reactions:                                                               =*
  !=       (a) N2O5 + hv -> NO3 + NO + O(3P)                                   =*
  !=       (b) N2O5 + hv -> NO3 + NO2                                          =*
  !=  Cross section from JPL97: use tabulated values up to 280 nm, use expon.  =*
  !=                            expression for >285nm, linearly interpolate    =*
  !=                            between s(280) and s(285,T) in between         =*
  !=  Quantum yield: Analysis of data in JPL94 (->dataj1/YLD/N2O5.qy)          =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !=  OPT    - If opt=1 read, otherwise just calc
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=04 !6,7

  INTEGER, INTENT(IN)             :: nw
  INTEGER, INTENT(IN)             :: nbl
  REAL, INTENT(IN)                :: wc(kw)
  INTEGER, INTENT(INOUT)          :: j
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob
  REAL, INTENT(IN) :: tLev(nbl,kz)

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  REAL :: xs, xs270, xs280, xst290, xst300
  REAL :: dum1, dum2
  REAL :: t(nbl)
  INTEGER :: i, iw

  !*************** N2O5 photodissociation
  j = j + 1
  !jlabel(j) = 'N2O5 -> NO3 + NO + O(3P)'
  j = j + 1
  !jlabel(j) = 'N2O5 -> NO3 + NO2'
  ! quantum yield : see dataj1/YLD/N2O5.qy for explanation
  ! correct for T-dependence of cross section

  DO iw = 1, nw - 1
    qy = MIN( 1., 3.832441 - 0.012809638 * wc(iw) )
    qy = MAX( 0., qy )
    DO i = 1, zLimit
      ! temperature dependence only valid for 225 - 300 K.
      DO iob=1,nbl
         t(nbl) = MAX(225.,MIN(tLev(iob,i),300.))
      END DO
      ! evaluation of exponential
      IF (wc(iw) >= 285. .AND. wc(iw) <= 380.) THEN
        DO iob=1,nbl
           sj(iob,7,i,iw) = qy * 1.e-20* &
	        EXP( 2.735 + (4728.5-17.127*wc(iw)) / t(nbl) ) !6
           sj(iob,8,i,iw) = (1.-qy) * 1.e-20* &
	        EXP( 2.735 + (4728.5-17.127*wc(iw)) / t(nbl) ) !7
	END DO
      ! between 280 and 285 nm:  Extrapolate from both sides, then average.
      ELSE IF (wc(iw) >= 280. .AND. wc(iw) < 285.) THEN
        DO iob=1,nbl
           xst290 = 1.e-20* EXP( 2.735 + (4728.5-17.127*290.) / t(nbl) )
           xst300 = 1.e-20* EXP( 2.735 + (4728.5-17.127*300.) / t(nbl) )
           dum1 = xs270 + (wc(iw) - 270.)*(xs280 - xs270)/10.
           dum2 = xst290 + (wc(iw) - 290.)*(xst300 - xst290)/10.
           xs = 0.5*(dum1 + dum2)
           sj(iob,7,i,iw) = qy * xs !6
           sj(iob,8,i,iw) = (1.-qy) * xs !7
        END DO
        ! for less than 280 nm, use tabulated values
      ELSE IF (wc(iw) < 280.) THEN
        DO iob=1,nbl
           sj(iob,7,i,iw) = qy * yg(iw,nReact) !6
           sj(iob,8,i,iw) = (1.-qy) * yg(iw,nReact) !7
        END DO
        ! beyond 380 nm, set to zero
      ELSE
        DO iob=1,nbl
           sj(iob,7,i,iw) = 0. !6
           sj(iob,8,i,iw) = 0. !7
        END DO
      END IF
    END DO
  END DO
  END SUBROUTINE r04

  !=============================================================================*

  SUBROUTINE r05(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide the product (cross section) x (quantum yield) for HNO2 photolysis=*
  !=     HNO2 + hv -> NO + OH                                                  =*
  !=  Cross section:  from JPL97                                               =*
  !=  Quantum yield:  assumed to be unity                                      =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  !=  EDIT HISTORY:                                                            =*
  !=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=05 !11

  INTEGER, INTENT(IN)             :: nw
  INTEGER, INTENT(IN)             :: nbl
  INTEGER, INTENT(INOUT)          :: j
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob
  REAL, INTENT(IN) :: tLev(nbl,kz)

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  INTEGER :: i, iw

  !*************** HNO2 photodissociation
  ! cross section from JPL92
  ! (from Bongartz et al., identical to JPL94, JPL97 recommendation)
  j = j + 1
  !jlabel(j) = 'HNO2 -> OH + NO'
  ! quantum yield = 1
  qy = 1.
  DO iw = 1, nw - 1
    DO i = 1, zLimit
       DO iob=1,nbl
          sj(iob,12,i,iw) = yg(iw,nReact)*qy !11
       END DO
    END DO
  END DO

  END SUBROUTINE r05

  !=============================================================================*

  SUBROUTINE r06(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product of (cross section) x (quantum yield) for HNO3 photolysis =*
  !=        HNO3 + hv -> OH + NO2                                              =*
  !=  Cross section: Burkholder et al., 1993                                   =*
  !=  Quantum yield: Assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=06 !12

  INTEGER, INTENT(IN)             :: nw
  INTEGER, INTENT(IN)             :: nbl
  INTEGER, INTENT(INOUT)          :: j
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob
  REAL, INTENT(IN) :: tLev(nbl,kz)

  INTEGER, PARAMETER :: kdata=100
  INTEGER :: i, iw

  !*************** HNO3 photodissociation
  j = j + 1
  !jlabel(j) = 'HNO3 -> OH + NO2'
  ! quantum yield = 1
  ! correct for temperature dependence

  DO iw = 1, nw - 1
    DO i = 1, zLimit
       DO iob=1,nbl
          sj(iob,13,i,iw) = yg1(iw,nReact) * &
	  1.e-20 * EXP( yg2(iw,nReact)/1.e3*&
	  (tLev(iob,i)-298.) ) !12
       END DO
    END DO
  END DO

  END SUBROUTINE r06

  !=============================================================================*

  SUBROUTINE r07(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product of (cross section) x (quantum yield) for HNO4 photolysis =*
  !=       HNO4 + hv -> HO2 + NO2                                              =*
  !=  Cross section:  from JPL97                                               =*
  !=  Quantum yield:  Assumed to be unity                                      =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=07 !13

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob
  REAL, INTENT(IN) :: tLev(nbl,kz)

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  INTEGER :: i, iw

  !*************** HNO4 photodissociation
  ! cross section from JPL85 (identical to JPL92 and JPL94 and JPL97)

  j = j + 1
  !jlabel(j) = 'HNO4 -> HO2 + NO2'

  ! quantum yield = 1

  qy = 1.
  DO iw = 1, nw - 1
    DO i = 1, zLimit
      DO iob=1,nbl
         sj(iob,14,i,iw) = yg(iw,nReact)*qy !13
      END DO
    END DO
  END DO

  END SUBROUTINE r07

  !=============================================================================*

  SUBROUTINE r08(nw,wl,wc,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product of (cross section) x (quantum yield) for H2O2 photolysis =*
  !=         H2O2 + hv -> 2 OH                                                 =*
  !=  Cross section:  From JPL97, tabulated values @ 298K for <260nm, T-depend.=*
  !=                  parameterization for 260-350nm                           =*
  !=  Quantum yield:  Assumed to be unity                                      =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=08 !10

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  REAL, INTENT(IN)                     :: wl(kw)
  REAL, INTENT(IN)                         :: wc(kw)
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob
  REAL, INTENT(IN) :: tLev(nbl,kz)

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  REAL :: a0, a1, a2, a3, a4, a5, a6, a7
  REAL :: b0, b1, b2, b3, b4
  REAL :: xs
  REAL :: t
  INTEGER :: i, iw
  REAL :: lambda
  REAL :: suma, sumb, chi

  !*************** H2O2 photodissociation

  ! cross section from Lin et al. 1978


  j = j + 1
  !jlabel(j) = 'H2O2 -> 2 OH'

  a0 = 6.4761E+04
  a1 = -9.2170972E+02
  a2 = 4.535649
  a3 = -4.4589016E-03
  a4 = -4.035101E-05
  a5 = 1.6878206E-07
  a6 = -2.652014E-10
  a7 = 1.5534675E-13

  b0 = 6.8123E+03
  b1 = -5.1351E+01
  b2 = 1.1522E-01
  b3 = -3.0493E-05
  b4 = -1.0924E-07

  ! quantum yield = 1

  qy = 1.

  DO iw = 1, nw - 1

  ! Parameterization (JPL94)
  ! Range 260-350 nm; 200-400 K

    IF ((wl(iw) >= 260.) .AND. (wl(iw) < 350.)) THEN

      lambda = wc(iw)
      suma = ((((((a7*lambda + a6)*lambda + a5)*lambda +  &
          a4)*lambda +a3)*lambda + a2)*lambda + a1)*lambda + a0
      sumb = (((b4*lambda + b3)*lambda + b2)*lambda + b1)*lambda + b0

      DO i = 1, zLimit
         DO iob=1,nbl
            t = MIN(MAX(tLev(iob,i),200.),400.)
            chi = 1./(1.+EXP(-1265./t))
            xs = (chi * suma + (1.-chi)*sumb)*1E-21
            sj(iob,11,i,iw) = xs*qy !10
         END DO
      END DO
    ELSE
      DO i = 1, zLimit
         DO iob=1,nbl
            sj(iob,11,i,iw) = yg(iw,nReact)*qy !10
         END DO
      END DO
    END IF

  END DO

  END SUBROUTINE r08

  !=============================================================================*

  SUBROUTINE r09(nw,wl,wc,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product of (cross section) x (quantum yield) for CHBr3 photolysis=*
  !=          CHBr3 + hv -> Products                                           =*
  !=  Cross section: Choice of data from Atlas (?Talukdar???) or JPL97         =*
  !=  Quantum yield: Assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=09 !51

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  REAL, INTENT(IN)                         :: wl(kw)
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob
  REAL, INTENT(OUT)                        :: wc(kw)
  REAL, INTENT(IN) :: tLev(nbl,kz)

  INTEGER, PARAMETER :: kdata=200
  REAL :: t
  REAL :: qy
  INTEGER :: iw
  INTEGER :: iz

  DO   iw = 1, nw - 1
    wc(iw) = (wl(iw) + wl(iw+1))/2.
  END DO

  !*************** CHBr3 photodissociation
  j = j + 1
  !jlabel(j) = 'CHBr3 -> Products'

  ! option:
  SELECT CASE(mOption(4))
     CASE (1)
        ! quantum yield = 1
        qy = 1.
        DO iw = 1, nw - 1
           DO iz = 1, zLimit
              DO iob=1,nbl
                 t = tLev(iob,iz)
                 IF (t >= 296.) THEN
                    yg(iw,nReact) = yg1(iw,nReact)
                 ELSE IF(t >= 286.) THEN
                    yg(iw,nReact) = yg1(iw,nReact) + (t-286.)* &
		                 (yg2(iw,nReact)-yg1(iw,nReact))/10.
                 ELSE IF(t >= 276.) THEN
                    yg(iw,nReact) = yg2(iw,nReact) + (t-276.)* &
		                 (yg3(iw,nReact)-yg2(iw,nReact))/10.
                 ELSE IF(t >= 266.) THEN
                    yg(iw,nReact) = yg3(iw,nReact) + (t-266.)* &
		                 (yg4(iw,nReact)-yg3(iw,nReact))/10.
                 ELSE IF(t >= 256.) THEN
                    yg(iw,nReact) = yg4(iw,nReact) + (t-256.)* &
		                 (yg5(iw,nReact)-yg4(iw,nReact))/10.
                 ELSE IF(t < 256.) THEN
                    yg(iw,nReact) = yg5(iw,nReact)
                 END IF
                    sj(iob,52,iz,iw) = yg(iw,nReact)*qy !51
             END DO
	   END DO
        END DO
     ! jpl97, with temperature dependence formula,
     !w = 290 nm to 340 nm,
     !T = 210K to 300 K
     !sigma, cm2 = exp((0.06183-0.000241*w)*(273.-T)-(2.376+0.14757*w))
     CASE (2)
       ! quantum yield = 1
       qy = 1.
       DO iw = 1, nw - 1
          DO iz = 1, zLimit
             DO iob=1,nbl
	        t = tLev(iob,iz)
                yg(iw,nReact) = yg1(iw,nReact)
                IF (wc(iw) > 290. .AND. wc(iw) < 340. .AND.  &
	                                t > 210 .AND. t < 300) THEN
                   yg(iw,nReact) = EXP((0.06183-0.000241*wc(iw))*(273.-t)-  &
                                   (2.376+0.14757*wc(iw)))
                END IF
                sj(iob,52,iz,iw) = yg(iw,nReact)*qy !51
             END DO
	  END DO
       END DO
    END SELECT

  END SUBROUTINE r09

  !=============================================================================*

  SUBROUTINE r10(nw,wl,wc,j,nbl,sj,tLev,AirDen)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product of (cross section) x (quantum yield) for CH2O photolysis =*
  !=        (a) CH2O + hv -> H + HCO                                           =*
  !=        (b) CH2O + hv -> H2 + CO                                           =*
  !=  Cross section: Choice between                                            =*
  !=                 1) Bass et al., 1980 (resolution: 0.025 nm)               =*
  !=                 2) Moortgat and Schneider (resolution: 1 nm)              =*
  !=                 3) Cantrell et al. (orig res.) for > 301 nm,              =*
  !=                    IUPAC 92, 97 elsewhere                                 =*
  !=                 4) Cantrell et al. (2.5 nm res.) for > 301 nm,            =*
  !=                    IUPAC 92, 97 elsewhere                                 =*
  !=                 5) Rogers et al., 1990                                    =*
  !=                 6) new NCAR recommendation, based on averages of          =*
  !=                    Cantrell et al., Moortgat and Schneider, and Rogers    =*
  !=                    et al.                                                 =*
  !=  Quantum yield: Choice between                                            =*
  !=                 1) Evaluation by Madronich 1991 (unpublished)             =*
  !=                 2) IUPAC 89, 92, 97                                       =*
  !=                 3) Madronich, based on 1), updated 1998.                  =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=10

  INTEGER, INTENT(IN)             :: nw
  INTEGER, INTENT(IN)             :: nbl
  REAL, INTENT(IN)                :: wl(kw)
  INTEGER, INTENT(INOUT)          :: j
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob
  REAL, INTENT(OUT)               :: wc(kw)
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(IN) :: airDen(nbl,kz)

  INTEGER, PARAMETER :: kdata=16000
  INTEGER :: iw
  REAL :: x(kdata)
  REAL :: phi1, phi2, phi20, ak300, akt
  REAL :: qy1, qy2(nbl)
  REAL :: sig(nbl), slope
  REAL :: t(nbl)
  INTEGER :: i

  DO   iw = 1, nw - 1
    wc(iw) = (wl(iw) + wl(iw+1))/2.
  END DO

  !***************************************************************
  !*************** CH2O photodissociatation

  j = j+1
  !jlabel(j) = 'CH2O -> H + HCO'

  j = j+1
  !jlabel(j) = 'CH2O -> H2 + CO'

  ! working grid arrays:
  !     yg1(:,nReact) = cross section at a specific temperature
  !     yg2(:,nReact), yg3(:,nReact) = cross sections at different temp or slope, for calculating
  !                temperature depedence
  !     yg4(:,nReact) = quantum yield data for radical channel
  !     yg5(:,nReact) = quantum yield data for molecular channel


  ! combine
  ! y1 = xsect
  ! y2 = xsect(223), Cantrell et al.
  ! y3 = xsect(293), Cantrell et al.
  ! y4 = qy for radical channel
  ! y5 = qy for molecular channel
  ! pressure and temperature dependent for w > 330.
  DO iw = 1, nw - 1
      IF (mOption(5) == 6) THEN
         DO iob=1,nbl
            sig(iob) = yg2(iw,nReact)
	 END DO
      ELSE
         DO iob=1,nbl
            sig(iob) = yg1(iw,nReact)
         END DO
      END IF
      DO i = 1, zLimit
        ! correct cross section for temperature dependence for > 301. nm
        IF (wl(iw) >= 301.) THEN
          DO iob=1,nbl
             t(iob) = MAX(223.15, MIN(tLev(iob,i), 293.15))
  	  END DO
	  IF (mOption(5) == 3 .OR. mOption(5) == 6) THEN
             DO iob=1,nbl
  	        sig(iob) = yg2(iw,nReact) + yg3(iw,nReact) * (t(iob) - 273.15)
  	     END DO
	  ELSE IF (mOption(5) == 4) THEN
             DO iob=1,nbl
  	        slope = (yg3(iw,nReact) - yg2(iw,nReact)) / (293. - 223.)
  	        sig(iob) = yg2(iw,nReact) + slope * (t(iob) - 223.)
             END DO
  	  END IF
        END IF
        DO iob=1,nbl
           sig(iob) = MAX(sig(iob), 0.)
        END DO
        ! quantum yields:
        ! temperature and pressure dependence beyond 330 nm
        qy1 = yg4(iw,nReact)
        IF ( (wc(iw) >= 330.) .AND. (yg5(iw,nReact) > 0.) ) THEN
           DO iob=1,nbl
  	      phi1 = yg4(iw,nReact)
  	      phi2 = yg5(iw,nReact)
  	      phi20 = 1. - phi1
  	      ak300=((1./phi2)-(1./phi20))/2.54E+19
  	      akt=ak300*(1.+61.69*(1.-tLev(iob,i)/300.)*  &
	                   (wc(iw)/329.-1.))
  	      qy2(iob) = 1. / ( (1./phi20) + airDen(iob,i)*akt)
           END DO
        ELSE
           DO iob=1,nbl
  	      qy2(iob) = yg5(iw,nReact)
           END DO
        END IF
        DO iob=1,nbl
           qy2(iob) = MAX(0.,qy2(iob))
           qy2(iob) = MIN(1.,qy2(iob))
           sj(iob,15,i,iw) = sig(iob) * qy1 !14
           sj(iob,16,i,iw) = sig(iob) * qy2(iob) !15
       END DO
    END DO
  END DO

  END SUBROUTINE r10

  !=============================================================================*

  SUBROUTINE r11(nw,j,nbl,sj,tLev,airDen)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for CH3CHO photolysis: =*
  !=      (a)  CH3CHO + hv -> CH3 + HCO                                        =*
  !=      (b)  CH3CHO + hv -> CH4 + CO                                         =*
  !=      (c)  CH3CHO + hv -> CH3CO + H                                        =*
  !=  Cross section:  Choice between                                           =*
  !=                   (1) IUPAC 97 data, from Martinez et al.                 =*
  !=                   (2) Calvert and Pitts                                   =*
  !=                   (3) Martinez et al., Table 1 scanned from paper         =*
  !=                   (4) KFA tabulations                                     =*
  !=  Quantum yields: Choice between                                           =*
  !=                   (1) IUPAC 97, pressure correction using Horowith and    =*
  !=                                 Calvert, 1982                             =*
  !=                   (2) NCAR data file, from Moortgat, 1986                 =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=11 !16,17,18

  INTEGER, INTENT(IN)    :: nw
  INTEGER, INTENT(IN)    :: nbl
  INTEGER, INTENT(INOUT) :: j
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(IN) :: airDen(nbl,kz)

  INTEGER, PARAMETER :: kdata=150
  INTEGER :: i
  REAL :: qy1, qy2, qy3
  REAL :: sig
  INTEGER :: iw

  !CH3CHO photolysis
  j = j+1
  !jlabel(j) = 'CH3CHO -> CH3 + HCO'
  j = j+1
  !jlabel(j) = 'CH3CHO -> CH4 + CO'
  j = j+1
  !jlabel(j) = 'CH3CHO -> CH3CO + H'

  ! combine:
  DO i = 1, zLimit
     DO iw = 1, nw - 1
        sig = yg(iw,nReact)
        ! quantum yields:
        qy2 = yg2(iw,nReact)
        qy3 = yg3(iw,nReact)
        DO iob=1,nbl
           qy1 = yg1(iw,nReact)
           ! pressure correction for channel 1, CH3 + CHO
           ! based on Horowitz and Calvert 1982.
           qy1 = qy1 * (1. + yg4(iw,nReact))/(1. + yg4(iw,nReact)* &
	             airDen(iob,i)/2.465E19)
           qy1 = MIN(1., qy1)
           qy1 = MAX(0., qy1)

           sj(iob,17,i,iw) = sig * qy1 !16


	   sj(iob,18,i,iw) = sig * qy2 !17
!LFR> !lfr debug
!LFR> 	   WRITE (88,FMT='(3(I5.5,1X),6(E18.8,1X))') &
!LFR> 	          iob,i,iw,sj(iob,16,i,iw), &
!LFR> 		  sj(iob,17,i,iw), &
!LFR> 		  yg(iw,nReact),yg1(iw,nReact),yg2(iw,nReact),yg3(iw,nReact)
!LFR> !lfr debug
           sj(iob,19,i,iw) = sig * qy3   !18
       END DO
    END DO
  END DO

  END SUBROUTINE r11

  !=============================================================================*

  SUBROUTINE r12(nw,j,nbl,sj,tLev,airDen)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:  							       =*
  !=  Provide the product (cross section) x (quantum yield) for C2H5CHO	     =*
  !=  photolysis: 							     =*
  !=	   C2H5CHO + hv -> C2H5 + HCO					     =*
  !=									     =*
  !=  Cross section:  Choice between					     =*
  !=		     (1) IUPAC 97 data, from Martinez et al.		     =*
  !=		     (2) Calvert and Pitts, as tabulated by KFA 	     =*
  !=  Quantum yield:  IUPAC 97 recommendation				     =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS: 							     =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working	  (I)=*
  !=	     wavelength grid						     =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in	  (I)=*
  !=	     working wavelength grid					     =*
  !=  WC     - REAL, vector of center points of wavelength intervals in	  (I)=*
  !=	     working wavelength grid					     =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level	  (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level	  (I)=*
  !=  J	   - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each	  (O)=*
  !=	     photolysis reaction defined, at each defined wavelength and     =*
  !=	     at each defined altitude level				     =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=	     defined							     =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=12 !19

  INTEGER, INTENT(IN)    :: nw
  INTEGER, INTENT(IN)    :: nbl
  INTEGER, INTENT(INOUT) :: j
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(IN) :: airDen(nbl,kz)

  INTEGER, PARAMETER :: kdata=150
  INTEGER :: i, n
  INTEGER :: n1
  REAL :: x1(kdata)
  REAL :: y1(kdata)
  REAL :: qy1(nbl)
  REAL :: sig(nbl)
  INTEGER :: ierr
  INTEGER :: iw

  !************************ C2H5CHO photolysis
  ! 1:  C2H5 + HCO

  j = j+1
  !jlabel(j) = 'C2H5CHO -> C2H5 + HCO'

  ! combine:
  DO iw = 1, nw - 1
    DO i = 1, zLimit
       DO iob=1,nbl
          sig(iob) = yg(iw,nReact)
       END DO
       ! quantum yields:
       ! use Ster-Volmer pressure dependence:
       IF (yg1(iw,nReact) < pzero) THEN
          DO iob=1,nbl
             qy1(iob) = 0.
	  END DO
       ELSE
          DO iob=1,nbl
             qy1(iob) = 1./(1. + (1./yg1(iw,nReact) - 1.)* &
	         airDen(iob,i)/2.45E19)
          END DO
       END IF
       DO iob=1,nbl
          qy1(iob) = MIN(qy1(iob),1.)
          sj(iob,20,i,iw) = sig(iob) * qy1(iob) !19
       END DO
    END DO
  END DO

  END SUBROUTINE r12

  !=============================================================================*

  SUBROUTINE r13(nw,wc,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide the product (cross section) x (quantum yield) for CHOCHO         =*
  !=  photolysis:                                                              =*
  !=              CHOCHO + hv -> Products                                      =*
  !=                                                                           =*
  !=  Cross section: Choice between                                            =*
  !=                  (1) Plum et al., as tabulated by IUPAC 97                =*
  !=                  (2) Plum et al., as tabulated by KFA.                    =*
  !=                  (3) Orlando et al.                                       =*
  !=                  (4) Horowitz et al., 2001                                =*
  !=  Quantum yield: IUPAC 97 recommendation                                   =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=13 !20,21


  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  REAL, INTENT(IN)                     :: wc(kw)
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob
  REAL, INTENT(IN) :: tLev(nbl,kz)

  INTEGER, PARAMETER :: kdata=500
  INTEGER :: i
  REAL :: qyii, qyiii
  REAL :: sig
  INTEGER :: iw

  !************************ CHOCHO photolysis
  ! see review by Madronich, Chapter VII in "The Mechansims of
  !  Atmospheric Oxidation of the Alkanes, Calvert et al, Oxford U.
  !  Press, 2000.
  ! Four possible channels:
  !     I     H2 + 2 CO
  !     II    2 HCO
  !     III   HCHO + CO
  !     IV    HCO + H + CO

  !  Based on that review, the following quantum yield assignments are made:

  !     qy_I = 0
  !     qy_II = 0.63 for radiation between 280 and 380 nm
  !     qy_III = 0.2  for radiation between 280 and 380 nm
  !     qy_IV = 0
  ! The yields for channels II and III were determined by Bauerle et al. (personal
  ! communication from G. Moortgat, still unpublished as of Dec 2000).
  ! Bauerle et al. used broad-band irradiation 280-380 nm.
  ! According to Zhu et al., the energetic threshold (for II) is 417 nm.  Therefore,
  ! here the quantum yields were set to zero for wc > 417.  Furthermore, the
  ! qys of Bauerle et al. were reduced to give the same J values when using full solar
  ! spectrum.  The reduction factor was calculated by comparing the J-values (for
  ! high sun) using the 380 and 417 cut offs.  The reduction factor is 7.1

  j = j+1
  !jlabel(j) = 'CHOCHO -> HCO + HCO'
  j = j + 1
  !jlabel(j) = 'CHOCHO -> CH2O + CO'

  ! combine:
  DO iw = 1, nw - 1
    sig = yg(iw,nReact)
    ! quantum yields:
    ! Use values from Bauerle, but corrected to cutoff at 417 rather than 380.
    ! this correction is a reduction by 7.1.
    ! so that qyII = 0.63/7.1  and qyIII = 0.2/7.1
    IF(wc(iw) < 417. ) THEN
      qyii = 0.089
      qyiii = 0.028
    ELSE
      qyii = 0.
      qyiii = 0.
    END IF
    DO i = 1, zLimit
       DO iob=1,nbl
          sj(iob,21,i,iw) = sig * qyii !20
          sj(iob,22,i, iw) = sig * qyiii !21
       END DO
    END DO
  END DO

  END SUBROUTINE r13

  !=============================================================================*

  SUBROUTINE r14(nw,wc,j,nbl,sj,tLev,airDen)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide the product (cross section) x (quantum yield) for CH3COCHO       =*
  !=  photolysis:                                                              =*
  !=           CH3COCHO + hv -> CH3CO + HCO                                    =*
  !=                                                                           =*
  !=  Cross section: Choice between                                            =*
  !=                  (1) from Meller et al., 1991, as tabulated by IUPAC 97   =*
  !=                         5 nm resolution (table 1) for < 402 nm            =*
  !=                         2 nm resolution (table 2) for > 402 nm            =*
  !=                  (2) average at 1 nm of Staffelbach et al., 1995, and     =*
  !=                      Meller et al., 1991                                  =*
  !=                  (3) Plum et al., 1983, as tabulated by KFA              =*
  !=                  (4) Meller et al., 1991 (0.033 nm res.), as tab. by KFA  =*
  !=                  (5) Meller et al., 1991 (1.0 nm res.), as tab. by KFA    =*
  !=                  (6) Staffelbach et al., 1995, as tabulated by KFA        =*
  !=  Quantum yield: Choice between                                            =*
  !=                  (1) Plum et al., fixed at 0.107                          =*
  !=                  (2) Plum et al., divided by 2, fixed at 0.0535           =*
  !=                  (3) Staffelbach et al., 0.45 for < 300 nm, 0 for > 430 nm=*
  !=                      linear interp. in between                            =*
  !=                  (4) Koch and Moortgat, prv. comm., 1997                  =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=14

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  REAL, INTENT(IN)                         :: wc(kw)
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(IN) :: airDen(nbl,kz)

  INTEGER, PARAMETER :: kdata=500
  INTEGER :: i
  REAL :: qy(nbl)
  REAL :: sig
  INTEGER :: iw
  REAL :: phi0, kq

  !************************ CH3COCHO photolysis
  ! 1:  CH3COCHO

  j = j+1
  !jlabel(j) = 'CH3COCHO -> CH3CO + HCO'

  ! combine:
  DO iw = 1, nw - 1
    sig = yg(iw,nReact)
    DO i = 1, zLimit
      ! quantum yields:
      IF (mOption(14) == 1) THEN
        qy = 0.107
      ELSE IF(mOption(14) == 2) THEN
        qy = 0.107/2.
      ELSE IF(mOption(14) == 3) THEN
        IF(wc(iw) <= 300.) THEN
          qy = 0.45
        ELSE IF (wc(iw) >= 430.) THEN
          qy = 0.
        ELSE
          qy = 0.45 + (0-0.45)*(wc(iw)-300.)/(430.-300.)
        END IF
      ELSE IF(mOption(14) == 4) THEN
        IF (yg1(iw,nReact) > 0.) THEN
           DO iob=1,nbl
             qy(iob) = yg2(iw,nReact)/( 1. + &
	           (airDen(iob,i)/2.465E19) * &
	           ( (yg2(iw,nReact)/yg1(iw,nReact)) - 1.))
           END DO
	ELSE
          qy = 0.
        END IF
      ELSE IF(mOption(14) == 5) THEN
        ! zero pressure yield:
        ! 1.0 for wc < 380 nm
        ! 0.0 for wc > 440 nm
        ! linear in between:
        phi0 = 1. - (wc(iw) - 380.)/60.
        phi0 = MIN(phi0,1.)
        phi0 = MAX(phi0,0.)
        ! Pressure correction: quenching coefficient, torr-1
        ! in air, Koch and Moortgat:
        kq = 1.36E8 * EXP(-8793/wc(iw))
        ! in N2, Chen et al:
        !               kq = 1.93e4 * EXP(-5639/wc(iw))
        IF(phi0 > 0.) THEN
          IF (wc(iw) >= 380. .AND. wc(iw) <= 440.) THEN
            DO iob=1,nbl
	       qy(iob) = phi0 / (phi0 + kq * airDen(iob,i) * &
	              760./2.456E19)
            END DO
	  ELSE
            qy = phi0
          END IF
        ELSE
          qy = 0.
        END IF
      END IF
      DO iob=1,nbl
         sj(iob,23,i,iw) = sig * qy(iob) !22
      END DO
    END DO
  END DO

  END SUBROUTINE r14

  !=============================================================================*

  SUBROUTINE r15(nw,wc,j,nbl,sj,tLev,airDen)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for CH3COCH3 photolysis=*
  !=          CH3COCH3 + hv -> Products                                        =*
  !=                                                                           =*
  !=  Cross section:  Choice between                                           =*
  !=                   (1) Calvert and Pitts                                   =*
  !=                   (2) Martinez et al., 1991, alson in IUPAC 97            =*
  !=                   (3) NOAA, 1998, unpublished as of 01/98                 =*
  !=  Quantum yield:  Choice between                                           =*
  !=                   (1) Gardiner et al, 1984                                =*
  !=                   (2) IUPAC 97                                            =*
  !=                   (3) McKeen et al., 1997                                 =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=15 !23


  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  REAL, INTENT(IN)                         :: wc(kw)
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob
  REAL, INTENT(IN) :: airDen(nbl,kz)

  INTEGER, PARAMETER :: kdata=150
  INTEGER :: i
  REAL :: qy(nbl)
  REAL :: sig(nbl)
  INTEGER :: iw
  REAL :: a, b, t, m, w
  REAL :: fco, fac
  REAL,PARAMETER :: deltat=298.-235.

  !*************** CH3COCH3 photodissociation
  j = j + 1
  !jlabel(j) = 'CH3COCH3 -> CH3CO + CH3'


  DO iw = 1, nw - 1
    DO i = 1, zLimit
      sig = yg(iw,nReact)
      IF(mOption(15) == 3) THEN
        DO iob=1,nbl
           t = 298. - tLev(iob,i)
           t = MIN(t, deltaT)
           t = MAX(t, 0.)
           sig(iob) = yg(iw,nReact)*(1. + yg2(iw,nReact)*t + &
	              yg3(iw,nReact)*t*t)
        END DO
      END IF
      IF (mOption(16) == 1) THEN
        DO iob=1,nbl
           qy(iob) = 0.0766 + 0.09415*EXP(-airDen(iob,i)/&
	             3.222E18)
        END DO
      ELSE IF (mOption(16) == 2) THEN
        qy = yg1(iw,nReact)
      ELSE IF (mOption(16) == 3) THEN
        IF (wc(iw) <= 292.) THEN
          qy = 1.
        ELSE IF (wc(iw) >= 292.  .AND. wc(iw) < 308. ) THEN
          a = -15.696 + 0.05707*wc(iw)
          b = EXP(-88.81+0.15161*wc(iw))
          DO iob=1,nbl
             qy(iob) = 1./(a + b*airDen(iob,i))
          END DO
        ELSE IF (wc(iw) >= 308.  .AND. wc(iw) < 337. ) THEN
          a = -130.2 + 0.42884*wc(iw)
          b = EXP(-55.947+0.044913*wc(iw))
          DO iob=1,nbl
             qy(iob) = 1./(a + b*airDen(iob,i))
          END DO
        ELSE IF (wc(iw) >= 337.) THEN
          qy = 0.
        END IF
        qy = MAX(0., qy)
        qy = MIN(1., qy)
      ELSE IF (mOption(16) == 4) THEN
        w = wc(iw)
        DO iob=1,nbl
	   t = tLev(iob,i)
           m = airDen(iob,i)
           CALL qyacet(w, t, m, fco, fac)
           qy(iob) = MAX(0., fac)
           qy(iob) = MIN(1., fac)
        END DO
      END IF
      DO iob=1,nbl
         sj(iob,24,i,iw) = sig(iob)*qy(iob) !23
      END DO
    END DO
  END DO

  END SUBROUTINE r15


  !=============================================================================*

  SUBROUTINE r16(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for CH3OOH photolysis: =*
  !=         CH3OOH + hv -> CH3O + OH                                          =*
  !=                                                                           =*
  !=  Cross section: Choice between                                            =*
  !=                  (1) JPL 97 recommendation (based on Vaghjiana and        =*
  !=                      Ravishankara, 1989), 10 nm resolution                =*
  !=                  (2) IUPAC 97 (from Vaghjiana and Ravishankara, 1989),    =*
  !=                      5 nm resolution                                      =*
  !=                  (3) Cox and Tyndall, 1978; only for wavelengths < 280 nm =*
  !=                  (4) Molina and Arguello, 1979;  might be 40% too high    =*
  !=  Quantum yield: Assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=16 !24

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  INTEGER :: i
  REAL :: qy
  INTEGER :: iw

  !*************** CH3OOH photodissociation
  j = j + 1
  !jlabel(j) = 'CH3OOH -> CH3O + OH'

  ! quantum yield = 1
  qy = 1.
  DO iw = 1, nw - 1

    DO i = 1, zLimit
      DO iob=1,nbl
         sj(iob,25,i,iw) = yg(iw,nReact)*qy !24
      END DO
    END DO
  END DO

  END SUBROUTINE r16

  !=============================================================================*

  SUBROUTINE r17(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for CH3ONO2            =*
  !=  photolysis:                                                              =*
  !=          CH3ONO2 + hv -> CH3O + NO2                                       =*
  !=                                                                           =*
  !=  Cross section: Choice between                                            =*
  !=                  (1) Calvert and Pitts, 1966                              =*
  !=                  (2) Talukdar, Burkholder, Hunter, Gilles, Roberts,       =*
  !=                      Ravishankara, 1997                                   =*
  !=                  (3) IUPAC 97, table of values for 198K                   =*
  !=                  (4) IUPAC 97, temperature-dependent equation             =*
  !=                  (5) Taylor et al, 1980                                   =*
  !=                  (6) fit from Roberts and Fajer, 1989                     =*
  !=                  (7) Rattigan et al., 1992                                =*
  !=                  (8) Libuda and Zabel, 1995                               =*
  !=  Quantum yield: Assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=17 !25

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata = 2000
  INTEGER :: i
  INTEGER :: iw
  REAL :: qy
  REAL :: sig(nbl)

  !*************** CH3ONO2 photodissociation
  j = j + 1
  !jlabel(j) = 'CH3ONO2 -> CH3O + NO2'

  ! quantum yield = 1
  qy = 1.

  DO iw = 1, nw - 1
    sig = yg(iw,nReact)
    DO i = 1, zLimit
      IF(mOption(18) == 2) THEN
        DO iob=1,nbl
           sig(iob) = yg(iw,nReact) * EXP (yg1(iw,nReact) * &
	         (tLev(iob,i)-298.))
        END DO
      ELSE IF (mOption(18) == 4) THEN
        DO iob=1,nbl
           sig(iob) = yg(iw,nReact)*10.**(yg1(iw,nReact)* &
	         tLev(iob,i))
        END DO
      END IF
      DO iob=1,nbl
         sj(iob,26,i,iw) = qy * sig(iob) !25
      END DO
    END DO
  END DO

  END SUBROUTINE r17

  !=============================================================================*

  SUBROUTINE r18(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for PAN photolysis:    =*
  !=       PAN + hv -> Products                                                =*
  !=                                                                           =*
  !=  Cross section: from Talukdar et al., 1995                                =*
  !=  Quantum yield: Assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=18 !26

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100

  INTEGER :: iw
  INTEGER :: i
  REAL :: qy
  REAL :: sig

  !*************** PAN photodissociation
  j = j+1
  !jlabel(j) = 'CH3CO(OONO2) -> Products'

  ! quantum yield
  ! yet unknown, but assumed to be 1.0 (Talukdar et al., 1995)
  qy = 1.0
  DO iw = 1, nw-1
    DO i = 1, zLimit
      DO iob=1,nbl
         sig = yg(iw,nReact) * EXP(yg2(iw,nReact)*&
	       (tLev(iob,i)-298.))
         sj(iob,27,i,iw) = qy * sig !26
      END DO
    END DO
  END DO

  END SUBROUTINE r18

  !=============================================================================*

  SUBROUTINE r19(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for CCl2O photolysis:  =*
  !=        CCl2O + hv -> Products                                             =*
  !=                                                                           =*
  !=  Cross section: JPL 94 recommendation                                     =*
  !=  Quantum yield: Unity (Calvert and Pitts)                                 =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=19 !31

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  INTEGER :: iw
  INTEGER :: iz

  !************ CCl2O photodissociation
  j = j+1
  !jlabel(j) = 'CCl2O -> Products'

  !** quantum yield unity (Calvert and Pitts)
  qy = 1.
  DO iw = 1, nw-1
    DO iz = 1, zLimit
      DO iob=1,nbl
         sj(iob,32,iz,iw) = qy * yg(iw,nReact) !31
      END DO
    END DO
  END DO

  END SUBROUTINE r19

  !=============================================================================*

  SUBROUTINE r20(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for CCl4 photolysis:   =*
  !=      CCl4 + hv -> Products                                                =*
  !=  Cross section: from JPL 97 recommendation                                =*
  !=  Quantum yield: assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=20

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  INTEGER :: iw!, n, idum
  INTEGER :: iz

  !************ CCl4 photodissociation
  j = j+1
  !jlabel(j) = 'CCl4 -> Products'

  !** quantum yield assumed to be unity
  qy = 1.
  DO iw = 1, nw-1
    DO iz = 1, zLimit
      DO iob=1,nbl
         sj(iob,33,iz,iw) = qy * yg(iw,nReact) !32
      END DO
    END DO
  END DO

  END SUBROUTINE r20

  !=============================================================================*

  SUBROUTINE r21(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for CClFO photolysis:  =*
  !=         CClFO + hv -> Products                                            =*
  !=  Cross section: from JPL 97                                               =*
  !=  Quantum yield: assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=21 !33

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  INTEGER :: iw!, n, idum
  INTEGER :: iz

  !************ CClFO photodissociation
  j = j+1
  !jlabel(j) = 'CClFO -> Products'

  !** quantum yield unity
  qy = 1.
  DO iw = 1, nw-1
    DO iz = 1, zLimit
      DO iob=1,nbl
         sj(iob,34,iz,iw) = qy * yg(iw,nReact) !33
      END DO
    END DO
  END DO

  END SUBROUTINE r21

  SUBROUTINE r22(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for CF2O photolysis:   =*
  !=        CF2O + hv -> Products                                              =*
  !=  Cross section:  from JPL 97 recommendation                               =*
  !=  Quantum yield:  unity (Nolle et al.)                                     =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=22 !34

  INTEGER, INTENT(IN)             :: nw
  INTEGER, INTENT(IN)             :: nbl
  INTEGER, INTENT(INOUT)   	  :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  INTEGER :: iw!, n, idum
  INTEGER :: iz

  !************ CF2O photodissociation
  j = j+1
  !jlabel(j) = 'CF2O -> Products'

  !** quantum yield unity (Nolle et al.)
  qy = 1.
  DO iw = 1, nw-1
    DO iz = 1, zLimit
      DO iob=1,nbl
         sj(iob,35,iz,iw) = qy * yg(iw,nReact) !34
      END DO
    END DO
  END DO

  END SUBROUTINE r22

  !=============================================================================*

  SUBROUTINE r23(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for CFC-113 photolysis:=*
  !=          CF2ClCFCl2 + hv -> Products                                      =*
  !=  Cross section:  from JPL 97 recommendation, linear interp. between       =*
  !=                  values at 210 and 295K                                   =*
  !=  Quantum yield:  assumed to be unity                                      =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=23 !35

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  REAL :: t
  INTEGER :: iw!, n, idum
  INTEGER :: iz
  REAL :: slope(nbl)

  !************ CF2ClCFCl2 (CFC-113) photodissociation
  j = j+1
  !jlabel(j) = 'CF2ClCFCl2 (CFC-113) -> Products'

  !** quantum yield assumed to be unity
  qy = 1.

  DO iz = 1, zLimit
    DO iob=1,nbl
       t = MAX(210.,MIN(tLev(iob,iz),295.))
       slope(iob) = (t-210.)/(295.-210.)
    END DO
    DO iw = 1, nw-1
      DO iob=1,nbl
         sj(iob,36,iz,iw) = qy * (yg2(iw,nReact) + &
                      slope(iob)*(yg1(iw,nReact)-yg2(iw,nReact))) !35
      END DO
    END DO
  END DO

  END SUBROUTINE r23

  !=============================================================================*

  SUBROUTINE r24(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for CFC-144 photolysis:=*
  !=              CF2ClCF2Cl + hv -> Products                                  =*
  !=  Cross section: from JPL 97 recommendation, linear interp. between values =*
  !=                 at 210 and 295K                                           =*
  !=  Quantum yield: assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=24 !36

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  REAL :: t
  INTEGER :: iw!, n, idum
  INTEGER :: iz
  REAL :: slope(nbl)

  !************ CF2ClCF2Cl (CFC-114) photodissociation
  j = j+1
  !jlabel(j) = 'CF2ClCF2Cl (CFC-114) -> Products'

  !*** cross sections from JPL97 recommendation (identical to 94 recommendation)

  !** quantum yield assumed to be unity
  qy = 1.

  DO iz = 1, zLimit
    DO iob=1,nbl
       t = MAX(210.,MIN(tLev(iob,iz),295.))
       slope(iob) = (t-210.)/(295.-210.)
    END DO
    DO iw = 1, nw-1
      DO iob=1,nbl
         sj(iob,37,iz,iw) = qy * (yg2(iw,nReact) + &
	         slope(iob)*(yg1(iw,nReact)-yg2(iw,nReact))) !36
      END DO
    END DO
  END DO

  END SUBROUTINE r24

  !=============================================================================*

  SUBROUTINE r25(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for CFC-115 photolysis =*
  !=             CF3CF2Cl + hv -> Products                                     =*
  !=  Cross section: from JPL 97 recommendation                                =*
  !=  Quantum yield: assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=25 !37

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  INTEGER :: iw!, n, idum
  INTEGER :: iz

  !*************************************************************
  !************ CF3CF2Cl (CFC-115) photodissociation

  j = j+1
  !jlabel(j) = 'CF3CF2Cl (CFC-115) -> Products'

  !*** cross sections from JPL97 recommendation (identical to 94 recommendation)

  !*** quantum yield assumed to be unity
  qy = 1.

  DO iw = 1, nw-1
    DO iz = 1, zLimit
      DO iob=1,nbl
         sj(iob,38,iz,iw) = qy * yg(iw,nReact) !37
      END DO
    END DO
  END DO

  END SUBROUTINE r25

  !=============================================================================*

  SUBROUTINE r26(nw,wc,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for CFC-111 photolysis =*
  !=          CCl3F + hv -> Products                                           =*
  !=  Cross section: from JPL 97 recommendation                                =*
  !=  Quantum yield: assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=26 !38

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  REAL, INTENT(IN)                     :: wc(kw)
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  REAL :: t(nbl)
  INTEGER :: iw!, n, idum
  INTEGER :: iz

  !************ CCl3F (CFC-11) photodissociation
  j = j+1
  !jlabel(j) = 'CCl3F (CFC-11) -> Products'

  !*** cross sections from JPL97 recommendation (identical to 94 recommendation)


  !*** quantum yield assumed to be unity
  qy = 1.

  DO iz = 1, zLimit
    DO iob=1,nbl
       t(iob) = 1E-04 * (tLev(iob,iz)-298.)
    END DO
    DO iw = 1, nw-1
      DO iob=1,nbl
         sj(iob,39,iz,iw) = qy * yg(iw,nReact) * &
	            EXP((wc(iw)-184.9) * t(iob)) !38
      END DO
    END DO
  END DO

  END SUBROUTINE r26

  !=============================================================================*

  SUBROUTINE r27(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for CFC-112 photolysis:=*
  !=         CCl2F2 + hv -> Products                                           =*
  !=  Cross section: from JPL 97 recommendation                                =*
  !=  Quantum yield: assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=27 !39

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  REAL :: t(nbl)
  INTEGER :: iw!, n, idum
  INTEGER :: iz

  !************ CCl2F2 (CFC-12) photodissociation
  j = j+1
  !jlabel(j) = 'CCl2F2 (CFC-12) -> Products'

  !*** cross sections from JPL97 recommendation (identical to 94 recommendation)
  !*** quantum yield assumed to be unity
  qy = 1.

  DO iz = 1, zLimit
    DO iob=1,nbl
       t(iob) = 1E-04 * (tLev(iob,iz)-298.)
    END DO
    DO iw = 1, nw-1
      DO iob=1,nbl
         sj(iob,40,iz,iw) = qy * yg(iw,nReact) * &
	        EXP((wc(iw)-184.9) * t(iob)) !39
      END DO
    END DO
  END DO

  END SUBROUTINE r27

  !=============================================================================*

  SUBROUTINE r28(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for CH3Br photolysis:  =*
  !=         CH3Br + hv -> Products                                            =*
  !=  Cross section: from JPL 97 recommendation                                =*
  !=  Quantum yield: assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=28 !50

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  INTEGER :: iw!, n, idum
  INTEGER :: iz

  !************ CH3Br photodissociation
  ! data from JPL97 (identical to 94 recommendation)
  j = j+1
  !jlabel(j) = 'CH3Br -> Products'
  !*** quantum yield assumed to be unity
  qy = 1.

  DO iw = 1, nw-1
    DO iz = 1, zLimit
      DO iob=1,nbl
         sj(iob,51,iz,iw) = qy * yg(iw,nReact) !50
      END DO
    END DO
  END DO

  END SUBROUTINE r28

  !=============================================================================*

  SUBROUTINE r29(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for CH3CCl3 photolysis =*
  !=           CH3CCl3 + hv -> Products                                        =*
  !=  Cross section: from JPL 97 recommendation, piecewise linear interp.      =*
  !=                 of data at 210, 250, and 295K                             =*
  !=  Quantum yield: assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=29

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  REAL :: t
  INTEGER :: iw!, n, idum
  INTEGER :: iz
  REAL :: slope

  !************ CH3CCl3 photodissociation
  j = j+1
  !jlabel(j) = 'CH3CCl3 -> Products'

  !*** cross sections from JPL97 recommendation (identical to 94 recommendation)


  !*** quantum yield assumed to be unity
  qy = 1.

  DO iz = 1, zLimit
    DO iob=1,nbl
       t = MIN(295.,MAX(tLev(iob,iz),210.))
       IF (t <= 250.) THEN
          slope = (t-210.)/(250.-210.)
          DO iw = 1, nw-1
             sj(iob,41,iz,iw) = qy * (yg3(iw,nReact) + &
	              slope*(yg2(iw,nReact)-yg3(iw,nReact))) !40
          END DO
       ELSE
          slope = (t-250.)/(295.-250.)
          DO iw = 1, nw-1
             sj(iob,41,iz,iw) = qy * (yg2(iw,nReact) + &
	              slope*(yg1(iw,nReact)-yg2(iw,nReact))) !40
          END DO
       END IF
    END DO
  END DO

  END SUBROUTINE r29

  !=============================================================================*

  SUBROUTINE r30(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for CH3Cl photolysis:  =*
  !=            CH3Cl + hv -> Products                                         =*
  !=  Cross section: from JPL 97 recommendation, piecewise linear interp.      =*
  !=                 from values at 255, 279, and 296K                         =*
  !=  Quantum yield: assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=30 !30

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  REAL :: t
  INTEGER :: iw!, n, idum
  INTEGER :: iz
  REAL :: slope

  !************ CH3Cl photodissociation
  j = j+1
  !jlabel(j) = 'CH3Cl -> Products'

  !*** cross sections from JPL97 recommendation (identical to 94 recommendation)

  !*** quantum yield assumed to be unity
  qy = 1.

  DO iz = 1, zLimit
    DO iob=1,nbl
       t = MAX(255.,MIN(tLev(iob,iz),296.))
       IF (t <= 279.) THEN
         slope = (t-255.)/(279.-255.)
         DO iw = 1, nw-1
           sj(iob,31,iz,iw) = qy * (yg3(iw,nReact)+ &
	                    slope*(yg2(iw,nReact)-yg3(iw,nReact))) !30
         END DO
       ELSE
         slope = (t-279.)/(296.-279.)
         DO iw = 1, nw-1
            sj(iob,31,iz,iw) = qy * (yg2(iw,nReact)+slope* &
	                   (yg1(iw,nReact)-yg2(iw,nReact))) !30
         END DO
      END IF
    END DO
  END DO

  END SUBROUTINE r30

  !=============================================================================*

  SUBROUTINE r31(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for ClOO photolysis:   =*
  !=          ClOO + hv -> Products                                            =*
  !=  Cross section: from JPL 97 recommendation                                =*
  !=  Quantum yield: assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=31 !27

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  INTEGER :: iw!, n, idum
  INTEGER :: iz

  !*************************************************************
  !************ ClOO photodissociation
  j = j+1
  !jlabel(j) = 'ClOO -> Products'

  !*** cross sections from JPL97 recommendation (identical to 94 recommendation)
  qy = 1.
  DO iw = 1, nw-1
    DO iz = 1, zLimit
      DO iob=1,nbl
         sj(iob,28,iz,iw) = qy * yg(iw,nReact) !27
      END DO
    END DO
  END DO

  END SUBROUTINE r31

  !=============================================================================*

  SUBROUTINE r32(nw,wc,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for HCFC-123 photolysis=*
  !=       CF3CHCl2 + hv -> Products                                           =*
  !=  Cross section: from Orlando et al., 1991                                 =*
  !=  Quantum yield: assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=32

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  REAL, INTENT(IN)                         :: wc(kw)
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  REAL :: qy
  REAL :: t(nbl)
  INTEGER :: i, iw!, idum
  INTEGER :: iz
  REAL :: lambda, sum_(nbl)

  !************ CF3CHCl2 (HCFC-123) photodissociation
  j = j+1
  !jlabel(j) = 'CF3CHCl2 (HCFC-123) -> Products'

  !*** cross sections from JPL94 recommendation

  !     OPEN(kin,FILE='dataj1/abs/HCFC-123_jpl94.abs',STATUS='OLD')
  !     READ(kin,*) idum, n
  !     DO i = 1, idum-2
  !       READ(kin,*)
  !     ENDDO
  !     DO i = 1, n
  !       READ(kin,*) x1(i), y1(i)
  !       y1(i) = y1(i) * 1E-20
  !     ENDDO
  !     CLOSE(kin)

  !     CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
  !     CALL addpnt(x1,y1,kdata,n,          0.,0.)
  !     CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
  !     CALL addpnt(x1,y1,kdata,n,        1E38,0.)

  !     CALL inter2(nw,wl,yg,n,x1,y1,ierr)

  !     IF (ierr .NE. 0) THEN
  !        WRITE(*,*) ierr, jlabel(j)
  !        STOP
  !     ENDIF

  !*** quantum yield assumed to be unity
  !     qy = 1.

  !     DO iw = 1, nw-1
  !       DO iz = 1, zLimit
  !         sq(j,iz,iw) = qy * yg(iw)
  !       ENDDO
  !     ENDDO


  !*** cross section from Orlando et al., 1991
  !*** quantum yield assumed to be unity

  qy = 1.

  DO iw = 1, nw-1
    lambda = wc(iw)
    ! use parameterization only up to 220 nm, as the error bars associated with
    ! the measurements beyond 220 nm are very large (Orlando, priv.comm.)
    IF (lambda >= 190. .AND. lambda <= 220.) THEN
      DO iz = 1, zLimit
        DO iob=1,nbl
           t(iob) = MIN(295.,MAX(tLev(iob,iz),203.))- &
	            tbar(nReact)
           sum_(iob) = 0.
	END DO
        DO i = 1, 4
          DO iob=1,nbl
             sum_(iob) = (coeff(i,1,nReact)+t(iob)*(coeff(i,2,nReact)+&
	                 t(iob)*coeff(i,3,nReact))) *  &
                         (lambda-lbar)**(i-1) + sum_(iob)
          END DO
	END DO
	DO iob=1,nbl
           sj(iob,42,iz,iw) = qy * EXP(sum_(iob)) !41
	END DO
      END DO
    ELSE
      DO iz = 1, zLimit
        DO iob=1,nbl
           sj(iob,42,iz,iw) = 0. !41
	END DO
      END DO
    END IF
  END DO

  END SUBROUTINE r32

  !=============================================================================*

  SUBROUTINE r33(nw,wc,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for HCFC-124 photolysis=*
  !=        CF3CHFCl + hv -> Products                                          =*
  !=  Cross section: from Orlando et al., 1991                                 =*
  !=  Quantum yield: assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=33 !42

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  REAL, INTENT(IN)                         :: wc(kw)
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  REAL :: qy
  REAL :: t(nbl)
  INTEGER :: i, iw!, n, idum
  INTEGER :: iz!, k
  REAL :: lambda, sum_(nbl)

  !*************************************************************
  !************ CF3CHFCl (HCFC-124) photodissociation

  j = j+1
  !jlabel(j) = 'CF3CHFCl (HCFC-124) -> Products'
  !*** cross section from Orlando et al., 1991
  !*** quantum yield assumed to be unity

  qy = 1.

  DO iw = 1, nw-1
    lambda = wc(iw)
    IF (lambda >= 190. .AND. lambda <= 230.) THEN
      DO iz = 1, zLimit
        DO iob=1,nbl
!          t = MIN(295.,MAX(tlev(i),203.))-tbar
           t(iob) = MIN(295.,MAX(tLev(iob,iz),203.))- &
	            tbar(nReact)
           sum_(iob) = 0.0
	END DO
        DO i = 1, 4
	  DO iob=1,nbl
             sum_(iob) = (coeff(i,1,nReact)+t(iob)*(coeff(i,2,nReact)+t(iob)* &
	                  coeff(i,3,nReact))) *  &
                          (lambda-lbar)**(i-1) + sum_(iob)
	  END DO
        END DO
	DO iob=1,nbl
           sj(iob,43,iz,iw) = qy * EXP(sum_(iob)) !42
	END DO
      END DO
    ELSE
      DO iz = 1, zLimit
        DO iob=1,nbl
           sj(iob,43,iz,iw) = 0. !42
        END DO
      END DO
    END IF
  END DO

  END SUBROUTINE r33

  !=============================================================================*

  SUBROUTINE r34(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for HCFC-141b          =*
  !=  photolysis:                                                              =*
  !=         CH3CFCl2 + hv -> Products                                         =*
  !=  Cross section: from JPL97 recommendation                                 =*
  !=  Quantum yield: assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=34 !43

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  INTEGER :: iw!, n, idum
  INTEGER :: iz

  !*************************************************************
  !************ CH3CFCl2 (HCFC-141b) photodissociation

  j = j+1
  !jlabel(j) = 'CH3CFCl2 (HCFC-141b) -> Products'

  !*** cross sections from JPL97 recommendation (identical to 94 recommendation)

  !*** quantum yield assumed to be unity
  qy = 1.
  DO iw = 1, nw-1
    DO iz = 1, zLimit
      DO iob=1,nbl
         sj(iob,44,iz,iw) = qy * yg(iw,nReact) !43
      END DO
    END DO
  END DO

  END SUBROUTINE r34

  !=============================================================================*

  SUBROUTINE r35(nw,wc,j,nbl,sj,tLev)
  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for HCFC-142b          =*
  !=  photolysis:                                                              =*
  !=          CH3CF2Cl + hv -> Products                                        =*
  !=  Cross section: from Orlando et al., 1991                                 =*
  !=  Quantum yield: assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=35 !44

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  REAL, INTENT(IN)                         :: wc(kw)
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  REAL :: qy
  REAL :: t(nbl)
  INTEGER :: i, iw!, n, idum
  INTEGER :: iz!, k
  REAL :: lambda, sum_(nbl)

  !*************************************************************
  !************ CH3CF2Cl (HCFC-142b) photodissociation

  j = j+1
  !jlabel(j) = 'CH3CF2Cl (HCFC-142b) -> Products'

  !*** cross sections from JPL94 recommendation

  !     OPEN(kin,FILE='dataj1/abs/HCFC-142b_jpl94.abs',STATUS='OLD')
  !     READ(kin,*) idum, n
  !     DO i = 1, idum-2
  !       READ(kin,*)
  !     ENDDO
  !     DO i = 1, n
  !       READ(kin,*) x1(i), y1(i)
  !       y1(i) = y1(i) * 1E-20
  !     ENDDO
  !     CLOSE(kin)

  !     CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
  !     CALL addpnt(x1,y1,kdata,n,          0.,0.)
  !     CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
  !     CALL addpnt(x1,y1,kdata,n,        1E38,0.)

  !     CALL inter2(nw,wl,yg,n,x1,y1,ierr)

  !     IF (ierr .NE. 0) THEN
  !       WRITE(*,*) ierr, jlabel(j)
  !       STOP
  !     ENDIF

  !*** quantum yield assumed to be unity
  !     qy = 1.

  !     DO iw = 1, nw-1
  !       DO iz = 1, zLimit
  !         sq(j,iz,iw) = qy * yg(iw)
  !       ENDDO
  !     ENDDO

  !*** cross section from Orlando et al., 1991

  !*** quantum yield assumed to be unity

  qy = 1.
  DO iw = 1, nw-1
    lambda = wc(iw)
    IF (lambda >= 190. .AND. lambda <= 230.) THEN
      DO iz = 1, zLimit
        DO iob=1,nbl
           t(iob) = MIN(295.,MAX(tLev(iob,iz),203.))- &
	            tbar(nReact)
           sum_(iob) = 0.
        END DO
	DO i = 1, 4
	  DO iob=1,nbl
             sum_(iob) = (coeff(i,1,nReact)+t(iob)*(coeff(i,2,nReact)+t(iob)* &
	                  coeff(i,3,nReact))) *(lambda-lbar)**(i-1) + sum_(iob)
          END DO
        END DO
        ! offeset exponent by 40 (exp(-40.) = 4.248e-18) to prevent exp. underflow errors
        ! on some machines.
        !             sq(j,iz,iw) = qy * EXP(sum)
        DO iob=1,nbl
	   sj(iob,45,iz,iw) = qy * 4.248E-18 * EXP(sum_(iob) + &
	                                       40.0) !44
        END DO
      END DO
    ELSE
      DO iz = 1, zLimit
        DO iob=1,nbl
           sj(iob,45,iz,iw) = 0.0 !44
        END DO
      END DO
    END IF
  END DO

  END SUBROUTINE r35

  !=============================================================================*

  SUBROUTINE r36(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for HCFC-225ca         =*
  !=  photolysis:                                                              =*
  !=           CF3CF2CHCl2 + hv -> Products                                    =*
  !=  Cross section: from JPL 97 recommendation                                =*
  !=  Quantum yield: assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=36 !45

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  INTEGER :: iw!, n!, idum
  INTEGER :: iz

  !*************************************************************
  !************ CF3CF2CHCl2 (HCFC-225ca) photodissociation

  j = j+1
  !jlabel(j) = 'CF3CF2CHCl2 (HCFC-225ca) -> Products'

  !*** cross sections from JPL97 recommendation (identical to 94 recommendation)


  !*** quantum yield assumed to be unity
  qy = 1.

  DO iw = 1, nw-1
    DO iz = 1, zLimit
      DO iob=1,nbl
         sj(iob,46,iz,iw) = qy * yg(iw,nReact) !45
      END DO
    END DO
  END DO

  END SUBROUTINE r36

  !=============================================================================*

  SUBROUTINE r37(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for HCFC-225cb         =*
  !=  photolysis:                                                              =*
  !=          CF2ClCF2CHFCl + hv -> Products                                   =*
  !=  Cross section: from JPL 97 recommendation                                =*
  !=  Quantum yield: assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=37 !46

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  INTEGER :: iw!, n, idum
  INTEGER :: iz

  !*************************************************************
  !************ CF2ClCF2CHFCl (HCFC-225cb) photodissociation

  j = j+1
  !jlabel(j) = 'CF2ClCF2CHFCl (HCFC-225cb) -> Products'

  !*** cross sections from JPL97 recommendation (identical to 94 recommendation)

  !*** quantum yield assumed to be unity
  qy = 1.

  DO iw = 1, nw-1
    DO iz = 1, zLimit
      DO iob=1,nbl
         sj(iob,47,iz,iw) = qy * yg(iw,nReact) !46
      END DO
    END DO
  END DO

  END SUBROUTINE r37

  !=============================================================================*

  SUBROUTINE r38(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for HCFC-22 photolysis =*
  !=          CHClF2 + hv -> Products                                          =*
  !=  Cross section: from JPL 97 recommendation, piecewise linear interp.      =*
  !=                 from values at 210, 230, 250, 279, and 295 K              =*
  !=  Quantum yield: assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=38 !47

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  REAL :: t
  INTEGER :: iw!, n, idum
  !INTEGER :: ierr
  INTEGER :: iz
  REAL :: slope

  !*************************************************************
  !************ CHClF2 (HCFC-22) photodissociation

  j = j+1
  !jlabel(j) = 'CHClF2 (HCFC-22) -> Products'

  !*** cross sections from JPL97 recommendation (identical to 94 recommendation)

  !*** quantum yield assumed to be unity
  qy = 1.

  DO iz = 1, zLimit
    DO iob=1,nbl
       t = MIN(295.,MAX(tLev(iob,iz),210.))
       IF (t <= 230.) THEN
          slope = (t-210.)/(230.-210.)
          DO iw = 1, nw-1
             sj(iob,48,iz,iw) = qy * (yg5(iw,nReact)+slope* &
	                    (yg4(iw,nReact)-yg5(iw,nReact))) !47
          END DO
       ELSE IF (t <= 250.) THEN
          slope = (t-230.)/(250.-230.)
          DO iw = 1, nw-1
             sj(iob,48,iz,iw) = qy * (yg4(iw,nReact)+slope*&
	                    (yg3(iw,nReact)-yg4(iw,nReact))) !47
          END DO
       ELSE IF (t <= 270.) THEN
          slope = (t-250.)/(270.-250.)
          DO iw = 1, nw-1
             sj(iob,48,iz,iw) = qy * (yg3(iw,nReact)+slope*&
	                    (yg2(iw,nReact)-yg3(iw,nReact))) !47
          END DO
       ELSE
          slope = (t-270.)/(295.-270.)
          DO iw = 1, nw-1
             sj(iob,48,iz,iw) = qy * (yg2(iw,nReact)+slope*&
	                    (yg1(iw,nReact)-yg2(iw,nReact))) !47
          END DO
       END IF
    END DO
  END DO

  END SUBROUTINE r38

  !=============================================================================*

  SUBROUTINE r39(nw,wc,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for HO2 photolysis:    =*
  !=          HO2 + hv -> OH + O                                               =*
  !=  Cross section: from JPL 97 recommendation                                =*
  !=  Quantum yield: assumed shape based on work by Lee, 1982; normalized      =*
  !=                 to unity at 248 nm                                        =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=39 !9

  INTEGER, INTENT(IN)      :: nw
  INTEGER, INTENT(IN)      :: nbl
  REAL, INTENT(IN)         :: wc(kw)
  INTEGER, INTENT(INOUT)   :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy(nbl)
  INTEGER :: iw!, n, idum
  INTEGER :: iz

  !*************************************************************
  !************ HO2 photodissociation
  j = j+1
  !jlabel(j) = 'HO2 -> OH + O'
  !*** cross sections from JPL97 recommendation (identical to 94 recommendation)

  !*** quantum yield:  absolute quantum yield has not been reported yet, but
  !***                 Lee measured a quantum yield for O(1D) production at 248
  !***                 nm that was 15 time larger than at 193 nm
  !*** here:  a quantum yield of unity is assumed at 248 nm and beyond, for
  !***        shorter wavelengths a linear decrease with lambda is assumed

  DO iw = 1, nw-1
    IF (wc(iw) >= 248.) THEN
       DO iob=1,nbl
          qy = 1.
       END DO
    ELSE
       DO iob=1,nbl
          qy = 1./15. + (wc(iw)-193.)*(14./15.)/(248.-193.)
          qy = MAX(qy,0.)
       END DO
    END IF
    DO iz = 1, zLimit
       DO iob=1,nbl
          sj(iob,10,iz,iw) = qy(iob) * yg(iw,nReact) !9
       END DO
    END DO
  END DO

  END SUBROUTINE r39

  !=============================================================================*

  SUBROUTINE r40(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) Halon-1202 photolysis: =*
  !=         CF2Br2 + hv -> Products                                           =*
  !=  Cross section: from JPL 97 recommendation                                =*
  !=  Quantum yield: unity (Molina and Molina)                                 =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=40 !54

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  INTEGER :: iw!, n, idum
  !INTEGER :: ierr
  INTEGER :: iz

  !*************************************************************
  !************ CF2Br2 (Halon-1202) photodissociation

  j = j+1
  !jlabel(j) = 'CF2Br2 (Halon-1202) -> Products'

  !*** cross sections from JPL97 recommendation (identical to 94 recommendation)

  !*** quantum yield unity (Molina and Molina)
  qy = 1.

  DO iw = 1, nw-1
    DO iz = 1, zLimit
      DO iob=1,nbl
         sj(iob,55,iz,iw) = qy * yg(iw,nReact) !54
      END DO
    END DO
  END DO

  END SUBROUTINE r40

  !=============================================================================*

  SUBROUTINE r41(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for Halon-1211         =*
  !=  photolysis:                                                              =*
  !=           CF2ClBr + hv -> Products                                        =*
  !=  Cross section: from JPL 97 recommendation                                =*
  !=  Quantum yield: assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=41 !55

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  INTEGER :: iw!, n, idum
  !INTEGER :: ierr
  INTEGER :: iz

  !*************************************************************
  !************ CF2BrCl (Halon-1211) photodissociation

  j = j+1
  !jlabel(j) = 'CF2BrCl (Halon-1211) -> Products'

  !*** cross sections from JPL97 recommendation (identical to 94 recommendation)
    !*** quantum yield assumed to be unity
  qy = 1.

  DO iw = 1, nw-1
    DO iz = 1, zLimit
      DO iob=1,nbl
         sj(iob,56,iz,iw) = qy * yg(iw,nReact) !55
      END DO
    END DO
  END DO

  END SUBROUTINE r41

  !=============================================================================*

  SUBROUTINE r42(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for Halon-1301         =*
  !=  photolysis:                                                              =*
  !=         CF3Br + hv -> Products                                            =*
  !=  Cross section: from JPL 97 recommendation                                =*
  !=  Quantum yield: assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=42 !52

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  INTEGER :: iw!, n, idum
  !INTEGER :: ierr
  INTEGER :: iz

  !*************************************************************
  !************ CF3Br (Halon-1301) photodissociation

  j = j+1
  !jlabel(j) = 'CF3Br (Halon-1301) -> Products'

  !*** cross sections from JPL97 recommendation (identical to 94 recommendation)

  !*** quantum yield assumed to be unity
  qy = 1.

  DO iw = 1, nw-1
    DO iz = 1, zLimit
      DO iob=1,nbl
         sj(iob,53,iz,iw) = qy * yg(iw,nReact) !52
      END DO
    END DO
  END DO

  END SUBROUTINE r42

  !=============================================================================*

  SUBROUTINE r43(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for Halon-2402         =*
  !=  photolysis:                                                              =*
  !=           CF2BrCF2Br + hv -> Products                                     =*
  !=  Cross section: from JPL 97 recommendation                                =*
  !=  Quantum yield: assumed to be unity                                       =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=43  !53

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy
  INTEGER :: iw!, n, idum
  !INTEGER :: ierr
  INTEGER :: iz

  !*************************************************************
  !************ CF2BrCF2Br (Halon-2402) photodissociation

  j = j+1
  !jlabel(j) = 'CF2BrCF2Br (Halon-2402) -> Products'

  !*** cross sections from JPL97 recommendation (identical to 94 recommendation)

  !*** quantum yield assumed to be unity
  qy = 1.

  DO iw = 1, nw-1
    DO iz = 1, zLimit
      DO iob=1,nbl
         sj(iob,54,iz,iw) = qy * yg(iw,nReact) !53
      END DO
    END DO
  END DO

  END SUBROUTINE r43

  !=============================================================================*

  SUBROUTINE r44(nw,wc,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for N2O photolysis:    =*
  !=              N2O + hv -> N2 + O(1D)                                       =*
  !=  Cross section: from JPL 97 recommendation                                =*
  !=  Quantum yield: assumed to be unity, based on Greenblatt and Ravishankara =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=44 !8



  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  REAL, INTENT(IN)                         :: wc(kw)
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob


  REAL :: qy
  REAL :: a, b!, c
  REAL :: a0, a1, a2, a3, a4
  REAL :: b0, b1, b2, b3
  REAL :: t
  INTEGER :: iw, iz
  REAL :: lambda

  !*************************************************************
  !************ N2O photodissociation
     j = j+1
     !jlabel(j) = 'N2O -> N2 + O(1D)'

  !*** cross sections according to JPL97 recommendation (identical to 94 rec.)
  !*** see file dataj1/abs/N2O_jpl94.abs for detail

  a0 = 68.21023
  a1 = -4.071805
  a2 = 4.301146E-02
  a3 = -1.777846E-04
  a4 = 2.520672E-07

  b0 = 123.4014
  b1 = -2.116255
  b2 = 1.111572E-02
  b3 = -1.881058E-05

  !*** quantum yield of N(4s) and NO(2Pi) is less than 1% (Greenblatt and
  !*** Ravishankara), so quantum yield of O(1D) is assumed to be unity
  qy = 1.

  DO iw = 1, nw-1
    lambda = wc(iw)
    IF (lambda >= 173. .AND. lambda <= 240.) THEN
      DO iz = 1, zLimit
         DO iob=1,nbl
            t = MAX(194.,MIN(tLev(iob,iz),320.))
            a = (((a4*lambda+a3)*lambda+a2)*lambda+a1)*lambda+a0
            b = (((b3*lambda+b2)*lambda+b1)*lambda+b0)
            b = (t-300.)*EXP(b)
            sj(iob,9,iz,iw) = qy * EXP(a+b) !8
         END DO
      END DO
    ELSE
      DO iz = 1, zLimit
         DO iob=1,nbl
            sj(iob,9,iz,iw) = 0. !8
         END DO
      END DO
    END IF
  END DO

  END SUBROUTINE r44

  !=============================================================================*

  SUBROUTINE r45(nw,wc,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for ClONO2 photolysis: =*
  !=        ClONO2 + hv -> Products                                            =*
  !=                                                                           =*
  !=  Cross section: JPL 97 recommendation                                     =*
  !=  Quantum yield: JPL 97 recommendation                                     =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=45 !28,29

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  REAL, INTENT(IN)                         :: wc(kw)
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=150

  REAL :: qy1, qy2
  REAL :: xs
  INTEGER :: iw!, n!, idum
  INTEGER :: iz

  !************ ClONO2 photodissociation

  j = j+1
  !jlabel(j) = 'ClONO2 -> Cl + NO3'

  !** cross sections from JPL97 recommendation

  DO iw = 1, nw-1
    !** quantum yields (from jpl97)
    IF( wc(iw) < 308.) THEN
      qy1 = 0.6
    ELSE IF( (wc(iw) >= 308) .AND. (wc(iw) <= 364.) ) THEN
      qy1 = 7.143E-3 * wc(iw) - 1.6
    ELSE IF( wc(iw) > 364. ) THEN
      qy1 = 1.0
    END IF
    qy2 = 1. - qy1
    ! compute T-dependent cross section
    DO iz = 1, zLimit
      DO iob=1,nbl
         xs = yg1(iw,nReact)*( 1. + yg2(iw,nReact)* &
	      (tLev(iob,iz)-296) +  &
              yg3(iw,nReact)*(tLev(iob,iz)-296)* &
	      (tLev(iob,iz)-296))
         sj(iob,29,iz,iw) = qy1 * xs !28
         sj(iob,30,iz,iw) = qy2 * xs !29
      END DO
    END DO
  END DO

  j = j+1
  !jlabel(j) = 'ClONO2 -> ClO + NO2'

  END SUBROUTINE r45

  !=============================================================================*

  SUBROUTINE r46(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for BrONO2 photolysis: =*
  !=        BrONO2 + hv -> Products                                            =*
  !=                                                                           =*
  !=  Cross section: JPL 03 recommendation                                     =*
  !=  Quantum yield: JPL 03 recommendation                                     =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=46 !48,49

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=100
  REAL :: qy1, qy2
  INTEGER :: iw!, n, idum
  !INTEGER :: ierr
  INTEGER :: iz

  !************ BrONO2 photodissociation

  j = j+1
  !jlabel(j) = 'BrONO2 -> BrO + NO2'
  j = j+1
  !jlabel(j) = 'BrONO2 -> Br + NO3'

    !** cross sections from JPL03 recommendation

  !** quantum yields (from jpl97)

  qy1 = 0.71
  qy2 = 0.29
  DO iw = 1, nw-1
    DO iz = 1, zLimit
      DO iob=1,nbl
         sj(iob,49,iz,iw) = qy1 * yg1(iw,nReact) !48
         sj(iob,50,iz,iw) = qy2 * yg1(iw,nReact) !49
      END DO
    END DO
  END DO

  END SUBROUTINE r46

  !=============================================================================*

  SUBROUTINE r47(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide product (cross section) x (quantum yield) for Cl2 photolysis:    =*
  !=        Cl2 + hv -> 2 Cl                                                   =*
  !=                                                                           =*
  !=  Cross section: JPL 97 recommendation                                     =*
  !=  Quantum yield: 1     (Calvert and Pitts, 1966)                           =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=47 !56

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=150
  REAL :: qy
  INTEGER :: iz, iw
  !INTEGER :: ierr


  !************ CL2 photodissociation

  j = j+1
  !jlabel(j) = 'Cl2 -> Cl + Cl'

  !** cross sections from JPL97 recommendation (as tab by Finlayson-Pitts
  ! and Pitts, 1999.

  !** quantum yield = 1 (Calvert and Pitts, 1966)

  qy = 1.
  DO iw = 1, nw-1
    DO iz = 1, zLimit
      DO iob=1,nbl
         sj(iob,57,iz,iw) = qy * yg(iw,nReact) !56
      END DO
    END DO
  END DO

  END SUBROUTINE r47

  !=============================================================================*

  SUBROUTINE r101(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide the product (cross section) x (quantum yield) for CH2(OH)CHO     =*
  !=  (glycolaldehye, hydroxy acetaldehyde) photolysis:                        =*
  !=           CH2(OH)CHO + hv -> Products                                     =*
  !=                                                                           =*
  !=  Cross section from                                                       =*
  != The Atmospheric Chemistry of Glycolaldehyde, C. Bacher, G. S. Tyndall     =*
  != and J. J. Orlando, J. Atmos. Chem., 39 (2001) 171-189.                    =*
  !=                                                                           =*
  !=  Quantum yield about 50%                                                  =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=101 !57

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=300
  INTEGER :: i!, n
  REAL :: qy
  INTEGER :: iw

  !************************ CH2(OH)CHO photolysis
  ! 1:  CH2(OH)CHO

  j = j+1
  !jlabel(j) = 'CH2(OH)CHO -> Products'

  ! combine:
  qy = 0.5

  DO iw = 1, nw - 1
    DO i = 1, zLimit
      DO iob=1,nbl
         sj(iob,58,i,iw) = yg(iw,nReact) * qy !57
      END DO
    END DO
  END DO

  END SUBROUTINE r101

  !=============================================================================*

  SUBROUTINE r102(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide the product (cross section) x (quantum yield) for CH3COCOCH3     =*
  !=  (biacetyl) photolysis:                                                   =*
  !=           CH3COCOCH3 + hv -> Products                                     =*
  !=                                                                           =*
  !=  Cross section from either                                                =*
  != 1.  Plum et al., Environ. Sci. Technol., Vol. 17, No. 8, 1983, p.480      =*
  != 2.  Horowitz et al., J. Photochem Photobio A, 146, 19-27, 2001.           =*
  !=                                                                           =*
  !=  Quantum yield =0.158                                                     =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=102 !58

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=300
  INTEGER :: i!, n
  REAL :: qy
  INTEGER :: iw

  !************************ CH3COCOCH3 photolysis
  ! 1:  CH3COCOCH3
  j = j+1
  !jlabel(j) = 'CH3COCOCH3 -> Products'

  ! quantum yield from Plum et al.

  qy = 0.158

  DO iw = 1, nw - 1
    DO i = 1, zLimit
      DO iob=1,nbl
         sj(iob,59,i,iw) = yg(iw,nReact) * qy !58
      END DO
    END DO
  END DO

  END SUBROUTINE r102

  !=============================================================================*

  SUBROUTINE r103(nw,wc,j,nbl,sj,tLev,airDen)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide the product (cross section) x (quantum yield) for CH3COCHCH2     =*
  !=  Methyl vinyl ketone photolysis:                                          =*
  !=           CH3COCHCH2 + hv -> Products                                     =*
  !=                                                                           =*
  !=  Cross section from                                                       =*
  != W. Schneider and G. K. Moorgat, priv. comm, MPI Mainz 1989 as reported by =*
  != Roeth, E.-P., R. Ruhnke, G. Moortgat, R. Meller, and W. Schneider,        =*
  != UV/VIS-Absorption Cross Sections and QUantum Yields for Use in            =*
  != Photochemistry and Atmospheric Modeling, Part 2: Organic Substances,      =*
  != Forschungszentrum Julich, Report Jul-3341, 1997.                          =*
  !=                                                                           =*
  !=  Quantum yield assumed unity                                              =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=103 !59

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  REAL, INTENT(IN)                     :: wc(kw)
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob
  REAL, INTENT(IN) :: airDen(nbl,kz)

  INTEGER, PARAMETER :: kdata=20000
  INTEGER :: i!, n
  REAL :: qy
  INTEGER :: iw

  !************************ CH3COCHCH2 photolysis

  j = j+1
  !jlabel(j) = 'CH3COCHCH2 -> Products'

  ! quantum yield from
  ! Gierczak, T., J. B. Burkholder, R. K. Talukdar, A. Mellouki, S. B. Barone,
  ! and A. R. Ravishankara, Atmospheric fate of methyl vinyl ketone and methacrolein,
  ! J. Photochem. Photobiol A: Chemistry, 110 1-10, 1997.
  ! depends on pressure and wavelength, set upper limit to 1.0

  DO iw = 1, nw - 1
    DO i = 1, zLimit
      DO iob=1,nbl
         qy = EXP(-0.055*(wc(iw)-308.)) / (5.5 + 9.2E-19* &
	     airDen(iob,i))
         qy = MIN(qy, 1.)
         sj(iob,60,i,iw) = yg(iw,nReact) * qy !59
      END DO
    END DO
  END DO

  END SUBROUTINE r103

  !=============================================================================*

  SUBROUTINE r104(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide the product (cross section) x (quantum yield) for CH2C(CH3)CHO   =*
  !=  methacrolein        photolysis:                                          =*
  !=           CH2C(CH3)CHO + hv -> Products                                   =*
  !=                                                                           =*
  !=  Cross section from                                                       =*
  != R. Meller, priv. comm, MPI Mainz 1990 as reported by =*
  != Roeth, E.-P., R. Ruhnke, G. Moortgat, R. Meller, and W. Schneider,        =*
  != UV/VIS-Absorption Cross Sections and QUantum Yields for Use in            =*
  != Photochemistry and Atmospheric Modeling, Part 2: Organic Substances,      =*
  != Forschungszentrum Julich, Report Jul-3341, 1997.                          =*
  !=                                                                           =*
  !=  Quantum yield assumed unity                                              =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=104 !60

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=20000
  INTEGER :: i!, n
  REAL :: qy
  INTEGER :: iw

  !************************ CH2C(CH3)CHO photolysis

  j = j+1
  !jlabel(j) = 'CH2C(CH3)CHO -> Products'

  ! quantum yield from
  ! Gierczak, T., J. B. Burkholder, R. K. Talukdar, A. Mellouki, S. B. Barone,
  ! and A. R. Ravishankara, Atmospheric fate of methyl vinyl ketone and methacrolein,
  ! J. Photochem. Photobiol A: Chemistry, 110 1-10, 1997.
  !   Upper limit, quantum yield < 0.01

  qy = 0.01

  DO iw = 1, nw - 1
    DO i = 1, zLimit
      DO iob=1,nbl
         sj(iob,61,i,iw) = yg(iw,nReact) * qy !60
      END DO
    END DO
  END DO
  !PRINT *,'React 104. sj=',sj(:,61,:,:)
  END SUBROUTINE r104

  !=============================================================================*

  SUBROUTINE r105(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide the product (cross section) x (quantum yield) for CH3COCO(OH)    =*
  !=  pyruvic acid        photolysis:                                          =*
  !=           CH3COCO(OH) + hv -> Products                                    =*
  !=                                                                           =*
  !=  Cross section from                                                       =*
  != Horowitz, A., R. Meller, and G. K. Moortgat, The UV-VIS absorption cross  =*
  != section of the a-dicarbonyl compounds: pyruvic acid, biacetyl, and        =*
  != glyoxal. J. Photochem. Photobiol. A:Chemistry, v.146, pp.19-27, 2001.     =*
  !=                                                                           =*
  !=  Quantum yield assumed unity                                              =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=105 !61

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=20000
  INTEGER :: i!, n
  REAL :: qy
  INTEGER :: iw

  !************************ CH3COCO(OH) photolysis

  j = j+1
  !jlabel(j) = 'CH3COCO(OH) -> Products'

  ! quantum yield  = 1

  qy = 1.

  DO iw = 1, nw - 1
    DO i = 1, zLimit
      DO iob=1,nbl
         sj(iob,62,i,iw) = yg(iw,nReact) * qy !61
      END DO
    END DO
  END DO

  END SUBROUTINE r105

  !=============================================================================*

  SUBROUTINE r106(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide the product (cross section) x (quantum yield) for CH3CH2ONO2     =*
  !=  ethyl nitrate       photolysis:                                          =*
  !=           CH3CH2ONO2 + hv -> CH3CH2O + NO2                                =*
  !=                                                                           =*
  != Absorption cross sections of several organic from                         =*
  != Talukdar, R. K., J. B. Burkholder, M. Hunter, M. K. Gilles,               =*
  != J. M Roberts, and A. R. Ravishankara, Atmospheric fate of several         =*
  != alkyl nitrates, J. Chem. Soc., Faraday Trans., 93(16) 2797-2805, 1997.    =*
  !=                                                                           =*
  !=  Quantum yield assumed unity                                              =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=106 !62

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=200
  INTEGER :: i!, n1, n2
  REAL :: qy, sig
  INTEGER :: iw

  !************************ CH3CH2ONO2 photolysis
  j = j+1
  !jlabel(j) = 'CH3CH2ONO2 -> CH3CH2O + NO2'

  ! quantum yield  = 1

  qy = 1.

  DO iw = 1, nw - 1
    DO i = 1, zLimit
      DO iob=1,nbl
         sig = yg1(iw,nReact)*EXP(yg2(iw,nReact)* &
	      (tLev(iob,i)- 298.))
         sj(iob,63,i,iw) = sig * qy !62
      END DO
    END DO
  END DO

  END SUBROUTINE r106

  !=============================================================================*

  SUBROUTINE r107(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide the product (cross section) x (quantum yield) for CH3CHONO2CH3   =*
  !=  isopropyl nitrate   photolysis:                                          =*
  !=           CH3CHONO2CH3 + hv -> CH3CHOCH3 + NO2                            =*
  !=                                                                           =*
  != Absorption cross sections of several organic from                         =*
  != Talukdar, R. K., J. B. Burkholder, M. Hunter, M. K. Gilles,               =*
  != J. M Roberts, and A. R. Ravishankara, Atmospheric fate of several         =*
  != alkyl nitrates, J. Chem. Soc., Faraday Trans., 93(16) 2797-2805, 1997.    =*
  !=                                                                           =*
  !=  Quantum yield assumed unity                                              =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=107 !63

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=200
  INTEGER :: i!, n1, n2
  REAL :: qy, sig
  !INTEGER :: ierr
  INTEGER :: iw

  !************************ CH3CHONO2CH3 photolysis
  j = j+1
  !jlabel(j) = 'CH3CHONO2CH3 -> CH3CHOCH3 + NO2'

  ! quantum yield  = 1

  qy = 1.

  DO iw = 1, nw - 1
    DO i = 1, zLimit
      DO iob=1,nbl
         sig = yg1(iw,nReact)*EXP(yg2(iw,nReact)* &
	      (tLev(iob,i)-298.))
         sj(iob,64,i,iw) = sig * qy !63
      END DO
    END DO
  END DO

  END SUBROUTINE r107

  !=============================================================================*

  SUBROUTINE r108(nw,wc,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide the product (cross section) x (quantum yield) for                =*
  !=   nitroxy ethanol CH2(OH)CH2(ONO2) + hv -> CH2(OH)CH2(O.) + NO2           =*
  !=                                                                           =*
  !=  Cross section from Roberts, J. R. and R. W. Fajer, UV absorption cross   =*
  !=    sections of organic nitrates of potential atmospheric importance and   =*
  !=    estimation of atmospheric lifetimes, Env. Sci. Tech., 23, 945-951,     =*
  !=    1989.
  !=                                                                           =*
  !=  Quantum yield assumed unity                                              =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=108 !64

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  REAL, INTENT(IN)                         :: wc(kw)
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  REAL :: qy, sig
  INTEGER :: iw, i
  REAL :: a, b, c

  !************************ CH2(OH)CH2(ONO2) photolysis

  j = j+1
  !jlabel(j) = 'CH2(OH)CH2(ONO2) -> CH2(OH)CH2(O.) + NO2'
  ! coefficients from Roberts and Fajer 1989, over 270-306 nm
  a = -2.359E-3
  b = 1.2478
  c = -210.4

  ! quantum yield  = 1

  qy = 1.

  DO iw = 1, nw - 1
    IF (wc(iw) >= 270. .AND. wc(iw) <= 306.) THEN
      sig = EXP(a*wc(iw)*wc(iw) + b*wc(iw) + c)
    ELSE
      sig = 0.
    END IF
    DO i = 1, zLimit
      DO iob=1,nbl
         sj(iob,65,i,iw) = sig * qy !64
      END DO
    END DO
  END DO

  END SUBROUTINE r108

  !=============================================================================*

  SUBROUTINE r109(nw,wc,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide the product (cross section) x (quantum yield) for                =*
  !=   nitroxy acetone CH3COCH2(ONO2) + hv -> CH3COCH2(O.) + NO2               =*
  !=                                                                           =*
  !=  Cross section from Roberts, J. R. and R. W. Fajer, UV absorption cross   =*
  !=    sections of organic nitrates of potential atmospheric importance and   =*
  !=    estimation of atmospheric lifetimes, Env. Sci. Tech., 23, 945-951,     =*
  !=    1989.
  !=                                                                           =*
  !=  Quantum yield assumed unity                                              =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=109 !65

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  REAL, INTENT(IN)                         :: wc(kw)
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  REAL :: qy, sig
  INTEGER :: iw, i
  REAL :: a, b, c

  !************************ CH3COCH2(ONO2) photolysis

  j = j+1
  !jlabel(j) = 'CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2'

  ! coefficients from Roberts and Fajer 1989, over 284-335 nm
  a = -1.365E-3
  b = 0.7834
  c = -156.8

  ! quantum yield  = 1
  qy = 1.
  DO iw = 1, nw - 1
    IF (wc(iw) >= 284. .AND. wc(iw) <= 335.) THEN
      sig = EXP(a*wc(iw)*wc(iw) + b*wc(iw) + c)
    ELSE
      sig = 0.
    END IF
    DO i = 1, zLimit
      DO iob=1,nbl
         sj(iob,66,i,iw) = sig * qy !65
      END DO
    END DO
  END DO

  END SUBROUTINE r109

  !=============================================================================*

  SUBROUTINE r110(nw,wc,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide the product (cross section) x (quantum yield) for                =*
  !=  t-butyl nitrate C(CH3)3(ONO2) + hv -> C(CH3)(O.) + NO2                   =*
  !=                                                                           =*
  !=  Cross section from Roberts, J. R. and R. W. Fajer, UV absorption cross   =*
  !=    sections of organic nitrates of potential atmospheric importance and   =*
  !=    estimation of atmospheric lifetimes, Env. Sci. Tech., 23, 945-951,     =*
  !=    1989.
  !=                                                                           =*
  !=  Quantum yield assumed unity                                              =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=110 !66

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  REAL, INTENT(IN)                         :: wc(kw)
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  REAL :: qy, sig
  INTEGER :: iw, i
  REAL :: a, b, c

  !************************ C(CH3)3(ONO2) photolysis
  j = j+1
  !jlabel(j) = 'C(CH3)3(ONO2) -> C(CH3)3(O.) + NO2'

  ! coefficients from Roberts and Fajer 1989, over 270-330 nm
  a = -0.993E-3
  b = 0.5307
  c = -115.5

  ! quantum yield  = 1
  qy = 1.
  DO iw = 1, nw - 1
    IF (wc(iw) >= 270. .AND. wc(iw) <= 330.) THEN
      sig = EXP(a*wc(iw)*wc(iw) + b*wc(iw) + c)
    ELSE
      sig = 0.
    END IF
    DO i = 1, zLimit
      DO iob=1,nbl
         sj(iob,67,i,iw) = sig * qy !66
      END DO
    END DO
  END DO

  END SUBROUTINE r110

  !=============================================================================*

  SUBROUTINE r111(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide the product (cross section) x (quantum yield) for ClOOCl         =*
  !=  ClO dimer           photolysis:                                          =*
  !=           ClOOCl + hv -> Cl + ClOO                                        =*
  !=                                                                           =*
  !=  Cross section from  JPL2002                                              =*
  !=                                                                           =*
  !=  Quantum yield assumed unity                                              =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=111 !67

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=20000
  INTEGER :: i!, n
  REAL :: qy
  INTEGER :: iw

  !************************ ClOOCl photolysis
  ! from JPL-2002

  j = j+1
  !jlabel(j) = 'ClOOCl -> Cl + ClOO'
  ! quantum yield  = 1
  qy = 1.

  DO iw = 1, nw - 1
    DO i = 1, zLimit
      DO iob=1,nbl
         sj(iob,68,i,iw) = yg(iw,nReact) * qy !67
      END DO
    END DO
  END DO

  END SUBROUTINE r111

  !=============================================================================*

  SUBROUTINE r112(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide the product (cross section) x (quantum yield) for hydroxyacetone =*
  !=  CH2(OH)COCH3        photolysis:                                          =*
  !=           CH2(OH)COCH3  -> CH3CO + CH2OH
  !=                         -> CH2(OH)CO + CH3                                =*
  !=                                                                           =*
  !=  Cross section from Orlando et al. (1999)                                 =*
  !=                                                                           =*
  !=  Quantum yield assumed 0.325 for each channel (J. Orlando, priv.comm.2003)=*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=112 !68,69

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=20000
  INTEGER :: i!, n
  REAL :: qy
  INTEGER :: iw

  !************************ CH2(OH)COCH3 photolysis
  ! from Orlando et al. 1999

  j = j+1
  !jlabel(j) = 'CH2(OH)COCH3 -> CH3CO + CH2(OH)'
  j = j+1
  !jlabel(j) = 'CH2(OH)COCH3 -> CH2(OH)CO + CH3'

  ! Total quantum yield  = 0.65, equal for each of the two channels
  qy = 0.325
  DO iw = 1, nw - 1
    DO i = 1, zLimit
      DO iob=1,nbl
         sj(iob,69,i,iw) = yg(iw,nReact) * qy !68
         sj(iob,70,i,iw) = yg(iw,nReact) * qy   !69
      END DO
    END DO
  END DO

  END SUBROUTINE r112

  !=============================================================================*

  SUBROUTINE r113(nw,wc,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide the product (cross section) x (quantum yield) for HOBr           =*
  !=  HOBr -> OH + Br                                                          =*
  !=  Cross section from JPL 2003                                              =*
  !=  Quantum yield assumed unity as in JPL2003                                =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=113 !70

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  REAL, INTENT(IN)                         :: wc(kw)
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER :: i
  REAL :: qy, sig
  INTEGER :: iw

  !************************ HOBr photolysis
  ! from JPL2003

  j = j+1
  !jlabel(j) = 'HOBr -> OH + Br'

  qy = 1.
  DO iw = 1, nw - 1
    sig = 24.77 * EXP( -109.80*(LOG(284.01/wc(iw)))**2 ) +  &
        12.22 * EXP(  -93.63*(LOG(350.57/wc(iw)))**2 ) +  &
        2.283 * EXP(- 242.40*(LOG(457.38/wc(iw)))**2 )
    sig = sig * 1.e-20
    IF(wc(iw) < 250. .OR. wc(iw) > 550.) sig = 0.
    DO i = 1, zLimit
      DO iob=1,nbl
         sj(iob,71,i,iw) = sig * qy !70
      END DO
    END DO
  END DO

  END SUBROUTINE r113

  !=============================================================================*

  SUBROUTINE r114(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide the product (cross section) x (quantum yield) for BrO            =*
  !=  BrO -> Br + O                                                            =*
  !=  Cross section from JPL 2003                                              =*
  !=  Quantum yield assumed unity as in JPL2003                                =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=114 !71

  INTEGER, INTENT(IN)             :: nw
  INTEGER, INTENT(IN)             :: nbl
  INTEGER, INTENT(INOUT)          :: j
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob
  REAL, INTENT(IN) :: tLev(nbl,kz)

  INTEGER :: i
  INTEGER :: iw
  REAL :: qy!, dum!, yg(kw)

  !************************ HOBr photolysis
  ! from JPL2003

  j = j+1
  !jlabel(j) = 'BrO -> Br + O'
  qy = 1.
  DO iw = 1, nw - 1
    DO i = 1, zLimit
      DO iob=1,nbl
         sj(iob,72,i,iw) = yg(iw,nReact) * qy !71
      END DO
    END DO
  END DO

  END SUBROUTINE r114

  !=============================================================================*

  SUBROUTINE r115(nw,j,nbl,sj,tLev)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Provide the product (cross section) x (quantum yield) for BrO            =*
  !=  Br2 -> Br + Br                                                           =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
  !=           working wavelength grid                                         =*
  !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
  !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
  !=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
  !=           photolysis reaction defined, at each defined wavelength and     =*
  !=           at each defined altitude level                                  =*
  !=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER,PARAMETER :: nReact=115 !72

  INTEGER, INTENT(IN)                      :: nw
  INTEGER, INTENT(IN)                      :: nbl
  INTEGER, INTENT(INOUT)                     :: j
  REAL, INTENT(IN) :: tLev(nbl,kz)
  REAL, INTENT(INOUT)             :: sj(nbl,kj,kz,kw)

  INTEGER :: iob

  INTEGER, PARAMETER :: kdata=50
  INTEGER :: i
  INTEGER :: iw!, ierr
  REAL :: qy !, yg(kw)

  !************************ Br2 photolysis

  j = j + 1
  !jlabel(j) = 'Br2 -> Br + Br'

  ! Absorption cross section from:
  ! Seery, D.J. and D. Britton, The continuous absorption spectra of chlorine,
  ! bromine, bromine chloride, iodine chloride, and iodine bromide, J. Phys.
  ! Chem. 68, p. 2263 (1964).

  qy = 1.
  DO iw = 1, nw - 1
    DO i = 1, zLimit
      DO iob=1,nbl
         sj(iob,73,i,iw) = yg(iw,nReact) * qy !72
      END DO
    END DO
  END DO

  END SUBROUTINE r115

  SUBROUTINE atrim(a,aa,n)
  IMPLICIT NONE

  ! Trim blanks from character string
  ! Input: a
  ! Output:  aa, n
  ! Internal: i

  CHARACTER (LEN=6), INTENT(IN)            :: a
  CHARACTER (LEN=6), INTENT(OUT)           :: aa
  INTEGER, INTENT(OUT)                     :: n

  INTEGER :: i

  aa = ' '
  n = 0
  DO i = 1, 6
    IF(a(i:i) /= ' ') THEN
      n = n + 1
      aa(n:n) = a(i:i)
    END IF
  END DO

  END SUBROUTINE atrim


  SUBROUTINE setz(nz,nj,coef,adjcoe,iang,iob,tLev,tcO3)
  !-----------------------------------------------------------------------------*
  !=  ZENITH = 0
  !-----------------------------------------------------------------------------*
  !=  INPUT:                                                                   =*
  !=
  !=  nz = height   = 121
  !=  nj = species  = 73
  !=  cz = overhead total ozone
  !=  tlev=Temperature
  !=  coef=POL coef(MZ,MS,MP)
  !=
  !=  OUTPUT:
  !=
  !=  ADJCOE - REAL, coross section adjust coefficients                        =*
  !=                                                                           =*
  !-----------------------------------------------------------------------------*
  !=  EDIT HISTORY:                                                            =*
  !=  08/2005 XUEXI                                                            =*
  !-----------------------------------------------------------------------------*
  != This program is free software;  you can redistribute it and/or modify     =*
  != it under the terms of the GNU General Public License as published by the  =*
  != Free Software Foundation;  either version 2 of the license, or (at your   =*
  != option) any later version.                                                =*
  != The TUV package is distributed in the hope that it will be useful, but    =*
  != WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
  != LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
  != License for more details.                                                 =*
  != To obtain a copy of the GNU General Public License, write to:             =*
  != Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
  !-----------------------------------------------------------------------------*
  != To contact the authors, please mail to:                                   =*
  != Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
  != send email to:  sasha@ucar.edu                                            =*
  !-----------------------------------------------------------------------------*
  != Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER, PARAMETER :: mz=5
  INTEGER, PARAMETER :: ms=73
  INTEGER, PARAMETER :: mp=5

  INTEGER, INTENT(IN) :: iang
  INTEGER, INTENT(IN) :: nz
  INTEGER, INTENT(IN) :: nj
  INTEGER, INTENT(IN) :: iob
  REAL, INTENT(IN)    :: coef(mz,ms,mp)
  REAL, INTENT(IN)    :: tLev
  REAL, INTENT(IN)    :: tco3(kz)

  REAL, INTENT(OUT)   :: adjcoe(kj,kz)


  INTEGER :: k,js,jp

  REAL :: xz(kz), c(5,kj)

  REAL :: c0,c1,c2
  REAL :: adjin
  INTEGER :: ij
  REAL :: tt

  !===============================================
  !  SET-UP
  !===============================================


  DO k=1,nz
    DO ij=1,nj
      adjcoe(ij,k)=1.00
    END DO
  END DO

  IF(iang>5) RETURN

  DO k=1,nz
    xz(k)=tcO3(k)*1.e-18
  END DO

  tt=tLev/281.0

  DO js=1,ms
    DO jp=1,mp
!      c(jp,js) = coef(1,js,jp)
      c(jp,js) = coef(iang,js,jp)
    END DO
  END DO

  ! CAL ADJ COEF

  !All species except tropospheric
  DO ij=1,27
    adjin=1.00
    CALL calcoe(nz,ij,c,xz,adjin,adjcoe)
  END DO
  DO ij=58,73
    adjin=1.00
    CALL calcoe(nz,ij,c,xz,adjin,adjcoe)
  END DO
  !(2) O3 -> O2 + O(1D)
  !----------------------------------------------------------------------
  !      Temperature Modification
  !      T0.9 (1.3) T0.95(1.25)  T1.0(1.2)  T1.15(1.18)  T1.1(1.16)
  !----------------------------------------------------------------------
  c0=ca0(iang)
  c1=ca1(iang)
  c2=ca2(iang)
  adjin = c0+c1*tt+c2*tt*tt+fatSum(iAng)
  CALL calcoe(nz,2,c,xz,adjin,adjcoe)
  ! 11 H2O2 -> 2 OH
  !----------------------------------------------------------------------
  !      Temperature Modification
  !      T0.9(0.95)  T0.95(0.975)    T1.0(1.0)  T1.15(0.105)   T1.1(1.10)
  !----------------------------------------------------------------------
  c0= cb0(iAng)
  c1= cb1(iAng)
  c2= cb2(iAng)
  adjin = c0+c1*tt+c2*tt*tt
  CALL calcoe(nz,11,c,xz,adjin,adjcoe)

  END SUBROUTINE setz

  SUBROUTINE setno2(nz,nw,no2xs,nbl,zLevel,no2col,cAir,dtNo2)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Set up an altitude profile of NO2 molecules, and corresponding absorption=*
  !=  optical depths.  Subroutine includes a shape-conserving scaling method   =*
  !=  that allows scaling of the entire profile to a given overhead NO2        =*
  !=  column amount.                                                           =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NO2NEW - REAL, overhead NO2 column amount (molec/cm^2) to which       (I)=*
  !=           profile should be scaled.  If NO2NEW < 0, no scaling is done    =*
  !=  NZ     - INTEGER, number of specified altitude levels in the working  (I)=*
  !=           grid                                                            =*
  !=  Z      - REAL, specified altitude working grid (km)                   (I)=*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  NO2XS  - REAL, molecular absoprtion cross section (cm^2) of O2 at     (I)=*
  !=           each specified wavelength                                       =*
  !=  TLAY   - REAL, temperature (K) at each specified altitude layer       (I)=*
  !=  DTNO2  - REAL, optical depth due to NO2 absorption at each            (O)=*
  !=           specified altitude at each specified wavelength                 =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nbl
  INTEGER, INTENT(IN) :: nz(nbl)
  INTEGER, INTENT(IN) :: nw
  REAL, INTENT(IN)    :: no2xs(kw)
  REAL, INTENT(IN)    :: zLevel(nbl,kz)
  REAL, INTENT(IN)    :: cAir(nbl,kz)
  REAL, INTENT(IN)    :: no2col(nbl)
  REAL, INTENT(OUT)   :: dtNo2(nbl,kz,kw)

  INTEGER, PARAMETER :: kdata=51
  REAL :: cz(nbl,kz)
  ! nitrogen dioxide profile data:
  REAL :: zd(kdata), no2(kdata)
  REAL :: cd(kdata)
  REAL :: hscale
  REAL :: colold
  REAL :: scale_(nbl)
  REAL :: sno2

  INTEGER :: i, l, nd, iob

  ! Example:  set to 1 ppb in lowest 1 km, set to zero above that.
  ! - do by specifying concentration at 3 altitudes.

  nd = 3
  zd(1) = 0.
  no2(1) = 1. * 2.69E10
  zd(2) = 1.
  no2(2) = 1. * 2.69E10
  zd(3) = zd(2)* 1.000001
  no2(3) = 10./largest
  ! compute column increments (alternatively, can specify these directly)
  DO   i = 1, nd - 1
    cd(i) = (no2(i+1)+no2(i)) * 1.e5 * (zd(i+1)-zd(i)) / 2.
  END DO
  ! Include exponential tail integral from top level to infinity.
  ! fold tail integral into top layer
  ! specify scale height near top of data (use ozone value)
  hscale = 4.50E5
  cd(nd-1) = cd(nd-1) + hscale * no2(nd)
  !********** end data input.

  ! Compute column increments and total column on standard z-grid.
  DO iob=1,nbl
     CALL inter3(nz(iob),zLevel(iob,:),cz(iob,:),nd,zd,cd, 1)
  END DO
  !*** Scaling of vertical profile by ratio of new to old column:
  ! If old column is near zero (less than 1 molec cm-2),
  ! use constant mixing ratio profile (nominal 1 ppt before scaling)
  ! to avoid numerical problems when scaling.

  DO iob=1,nbl
     IF(fsum(nz(iob)-1,cz(iob,:)) < 1.) THEN
        DO i = 1, nz(iob)-1
           cz(iob,i) = 1.e-12 * cAir(iob,i)
        END DO
     END IF
  END DO
  DO iob=1,nbl
     colold = fsum(nz(iob)-1, cz(iob,:))
     scale_(iob) =  2.687E16 * no2col(iob) / colold
  END DO
  DO i = 1, zLimit
     DO iob=1,nbl
        IF(i>nz(iob)-1) CYCLE
        cz(iob,i) = cz(iob,i) * scale_(iob)
     END DO
  END DO

  !***********************************
  ! calculate optical depth for each layer.  Output: dtno2(kz,kw)
  DO   l = 1, nw-1
    sno2 = no2xs(l)
    DO   i = 1, zLimit
       DO iob=1,nbl
          IF(i>nz(iob)-1) CYCLE
          dtNo2(iob,i,l) = cz(iob,i)*sno2
       END DO
    END DO
  END DO
  !_______________________________________________________________________

  END SUBROUTINE setno2
  !=============================================================================*

  SUBROUTINE setso2(nz,nw,so2xs,nbl,zLevel,so2col,cAir,dtSo2)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Set up an altitude profile of SO2 molecules, and corresponding absorption=*
  !=  optical depths.  Subroutine includes a shape-conserving scaling method   =*
  !=  that allows scaling of the entire profile to a given overhead SO2        =*
  !=  column amount.                                                           =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  SO2NEW - REAL, overhead SO2 column amount (molec/cm^2) to which       (I)=*
  !=           profile should be scaled.  If SO2NEW < 0, no scaling is done    =*
  !=  NZ     - INTEGER, number of specified altitude levels in the working  (I)=*
  !=           grid                                                            =*
  !=  Z      - REAL, specified altitude working grid (km)                   (I)=*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  SO2XS  - REAL, molecular absoprtion cross section (cm^2) of O2 at     (I)=*
  !=           each specified wavelength                                       =*
  !=  TLAY   - REAL, temperature (K) at each specified altitude layer       (I)=*
  !=  DTSO2  - REAL, optical depth due to SO2 absorption at each            (O)=*
  !=           specified altitude at each specified wavelength                 =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nw
  INTEGER, INTENT(IN) :: nbl
  INTEGER, INTENT(IN) :: nz(nbl)
  REAL, INTENT(IN)    :: so2xs(kw)
  REAL, INTENT(IN)    :: zLevel(nbl,kz)
  REAL, INTENT(IN)    :: cAir(nbl,kz)
  REAL, INTENT(IN)    :: so2col(nbl)
  REAL, INTENT(OUT)   :: dtSo2(nbl,kz,kw)

  INTEGER, PARAMETER :: kdata=51
  REAL :: cz(nbl,kz)
  ! sulfur dioxide profile data:
  REAL :: zd(kdata), so2(kdata)
  REAL :: cd(kdata)
  REAL :: hscale
  REAL :: colold
  REAL :: scale_(nbl)
  REAL :: sso2

  INTEGER :: i, l, nd,iob

  ! Example:  set to 1 ppb in lowest 1 km, set to zero above that.
  ! - do by specifying concentration at 3 altitudes.

  nd = 3
  zd(1) = 0.
  so2(1) = 1. * 2.69E10
  zd(2) = 1.
  so2(2) = 1. * 2.69E10
  zd(3) = zd(2)* 1.000001
  so2(3) = 10./largest
  ! compute column increments (alternatively, can specify these directly)
  DO   i = 1, nd - 1
    cd(i) = (so2(i+1)+so2(i)) * 1.e5 * (zd(i+1)-zd(i)) / 2.
  END DO
  ! Include exponential tail integral from top level to infinity.
  ! fold tail integral into top layer
  ! specify scale height near top of data (use ozone value)
  hscale = 4.50E5
  cd(nd-1) = cd(nd-1) + hscale * so2(nd)
  !********** end data input.

  ! Compute column increments on standard z-grid.
  DO iob=1,nbl
     CALL inter3(nz(iob),zLevel(iob,:),cz(iob,:), nd,zd,cd, 1)
  END DO
  !*** Scaling of vertical profile by ratio of new to old column:
  ! If old column is near zero (less than 1 molec cm-2),
  ! use constant mixing ratio profile (nominal 1 ppt before scaling)
  ! to avoid numerical problems when scaling.

  DO iob=1,nbl
     IF(fsum(nz(iob)-1,cz(iob,:)) < 1.) THEN
        DO i = 1, nz(iob)-1
           cz(iob,i) = 1.e-12 * cAir(iob,i)
        END DO
     END IF
  END DO
  DO iob=1,nbl
     colold = fsum(nz(iob)-1,cz(iob,:))
     scale_(iob) =  2.687E16 * so2col(iob) / colold
  END DO
  DO i = 1, zLimit
     DO iob=1,nbl
        IF(i>nz(iob)-1) CYCLE
        cz(iob,i) = cz(iob,i) * scale_(iob)
     END DO
  END DO

  !***********************************
  ! calculate sulfur optical depth for each layer, with optional temperature
  ! correction.  Output, dtso2(kz,kw)

  DO   l = 1, nw-1
    sso2 = so2xs(l)
    DO   i = 1, zLimit
       ! Leaving this part in in case i want to interpolate between
       ! the 221K and 298K data.
       !            IF ( wl(l) .GT. 240.5  .AND. wl(l+1) .LT. 350. ) THEN
       !               IF (tlay(i) .LT. 263.) THEN
       !                  sso2 = s221(l) + (s263(l)-s226(l)) / (263.-226.) *
       !     $                 (tlay(i)-226.)
       !               ELSE
       !                  sso2 = s263(l) + (s298(l)-s263(l)) / (298.-263.) *
       !     $              (tlay(i)-263.)
       !               ENDIF
       !            ENDIF
       DO iob=1,nbl
          IF(i>nz(iob)-1) CYCLE
          dtSo2(iob,i,l) = cz(iob,i)*sso2
       END DO
    END DO
  END DO
  !_______________________________________________________________________

  END SUBROUTINE setso2

  SUBROUTINE sphers(nz,nbl,zLevel,sza,nid,dsdh)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Calculate slant path over vertical depth ds/dh in spherical geometry.    =*
  !=  Calculation is based on:  A.Dahlback, and K.Stamnes, A new spheric model =*
  !=  for computing the radiation field available for photolysis and heating   =*
  !=  at twilight, Planet.Space Sci., v39, n5, pp. 671-683, 1991 (Appendix B)  =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
  !=            grid                                                           =*
  !=  Z       - REAL, specified altitude working grid (km)                  (I)=*
  !=  ZEN     - REAL, solar zenith angle (degrees)                          (I)=*
  !=  DSDH    - REAL, slant path of direct beam through each layer crossed  (O)=*
  !=            when travelling from the top of the atmosphere to layer i;     =*
  !=            DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                            =*
  !=  NID     - INTEGER, number of layers crossed by the direct beam when   (O)=*
  !=            travelling from the top of the atmosphere to layer i;          =*
  !=            NID(i), i = 0..NZ-1                                            =*
  !-----------------------------------------------------------------------------*
  !=  EDIT HISTORY:                                                            =*
  !=  double precision fix for shallow layers - Julia Lee-Taylor Dec 2000      =*
  !-----------------------------------------------------------------------------*
  != This program is free software;  you can redistribute it and/or modify     =*
  != it under the terms of the GNU General Public License as published by the  =*
  != Free Software Foundation;  either version 2 of the license, or (at your   =*
  != option) any later version.                                                =*
  != The TUV package is distributed in the hope that it will be useful, but    =*
  != WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
  != LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
  != License for more details.                                                 =*
  != To obtain a copy of the GNU General Public License, write to:             =*
  != Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
  !-----------------------------------------------------------------------------*
  != To contact the authors, please mail to:                                   =*
  != Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
  != send email to:  sasha@ucar.edu                                            =*
  !-----------------------------------------------------------------------------*
  != Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nbl
  INTEGER, INTENT(IN) :: nz(nbl)
  REAL   , INTENT(IN) :: zLevel(nbl,kz)
  REAL   , INTENT(IN) :: sza(nbl)
  INTEGER,INTENT(OUT) :: nid(nbl,0:kz)
  REAL   ,INTENT(OUT) :: dsdh(nbl,0:kz,kz)

  REAL :: re(nbl)
  REAL :: ze(nbl,kz)
  DOUBLE PRECISION :: zenrad(nbl)
  DOUBLE PRECISION :: rpsinz, rj, rjp1, dsj, dhj, ga, gb, sm
  INTEGER :: i, j, k, iob
  INTEGER :: id

  INTEGER :: nlayer
  REAL :: zd(nbl,0:kz-1)

  DO iob=1,nbl
     zenrad(iob) = sza(iob)*rpd
  END DO
  ! number of layers:
  !nlayer = nz(iob) - 1

  ! include the elevation above sea level to the radius of the earth:
  DO iob=1,nbl
     re(iob) = radius + zLevel(iob,1)
  END DO
  ! correspondingly z changed to the elevation above earth surface:
  DO k = 1, zLimit
     DO iob=1,nbl
        IF(k>nz(iob)) CYCLE
        ze(iob,k) = zLevel(iob,k) - zLevel(iob,1)
     END DO
  END DO

  ! inverse coordinate of z
  DO iob=1,nbl
     zd(iob,0) = ze(iob,nz(iob))
     DO k = 1, nz(iob)-1
        zd(iob,k) = ze(iob,nz(iob) - k)
     END DO
  END DO

  ! initialize dsdh(i,j), nid(i)
  DO i = 0, kz
     DO iob=1,nbl
        nid(iob,i) = 0
     END DO
     DO j = 1, kz
        DO iob=1,nbl
           dsdh(iob,i,j) = 0.
        END DO
     END DO
  END DO

  ! calculate ds/dh of every layer
  DO  i = 0, zLimit
     DO iob=1,nbl
        IF(i>nz(iob)-1) CYCLE
        rpsinz = (re(iob) + zd(iob,i)) * SIN(zenrad(iob))
        nid(iob,i) = -1
        IF ( .not. (sza(iob) > 90.0) .AND. &
                    (rpsinz < re(iob)) ) THEN
        !ELSE
           ! Find index of layer in which the screening height lies
           id = i
           IF( sza(iob) > 90.0 ) THEN
              DO  j = 1, nz(iob) - 1
                 IF( (rpsinz < ( zd(iob,j-1) + re(iob) ) ) .AND.  &
                     (rpsinz >= ( zd(iob,j) + re(iob) )) ) id = j
              END DO
           END IF
           DO  j = 1, id
              sm = 1.0
              IF(j == id .AND. id == i .AND. &
	           sza(iob) > 90.0) sm = -1.0
              rj = re(iob) + zd(iob,j-1)
              rjp1 = re(iob) + zd(iob,j)
              dhj = zd(iob,j-1) - zd(iob,j)
	      IF(dhj==0.0) CYCLE !LFR for teste
              ga = rj*rj - rpsinz*rpsinz
              gb = rjp1*rjp1 - rpsinz*rpsinz
              IF (ga < 0.0) ga = 0.0
              IF (gb < 0.0) gb = 0.0
              IF(id > i .AND. j == id) THEN
                  dsj = SQRT( ga )
              ELSE
                  dsj = SQRT( ga ) - sm*SQRT( gb )
              END IF
              dsdh(iob,i,j) = dsj / dhj
           END DO
           nid(iob,i) = id
         END IF
      END DO
  END DO

  END SUBROUTINE sphers

  !=============================================================================*
  SUBROUTINE airmas(nz,nbl,cAir,scol,vcol,nid,dsdh)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Calculate vertical and slant air columns, in spherical geometry, as a    =*
  !=  function of altitude.                                                    =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
  !=            grid                                                           =*
  !=  DSDH    - REAL, slant path of direct beam through each layer crossed  (O)=*
  !=            when travelling from the top of the atmosphere to layer i;     =*
  !=            DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                            =*
  !=  NID     - INTEGER, number of layers crossed by the direct beam when   (O)=*
  !=            travelling from the top of the atmosphere to layer i;          =*
  !=            NID(i), i = 0..NZ-1                                            =*
  !=  VCOL    - REAL, output, vertical air column, molec cm-2, above level iz  =*
  !=  SCOL    - REAL, output, slant air column in direction of sun, above iz   =*
  !=            also in molec cm-2                                             =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nbl
  INTEGER, INTENT(IN) :: nz(nbl)
  INTEGER,INTENT(IN)  :: nid(nbl,0:kz)
  REAL, INTENT(IN)    :: cAir(nbl,kz)
  REAL, INTENT(IN)    :: dsdh(nbl,0:kz,kz)
  REAL, INTENT(OUT)   :: scol(nbl,kz)
  REAL, INTENT(OUT)   :: vcol(nbl,kz)

  INTEGER :: id, j, iob
  REAL :: sum(nbl), vsum(nbl)

  ! calculate vertical and slant column from each level:
  ! work downward
  vsum = 0.
  DO iob=1,nbl
     DO id = 0, nz(iob) - 1
        vsum(iob) = vsum(iob) + &
	         cAir(iob,nz(iob)-id)
        vcol(iob,nz(iob)-id) = vsum(iob)
        sum(iob) = 0.
        IF(nid(iob,id) < 0) THEN
          sum = largest
        ELSE
           ! single pass layers:
           DO j = 1, MIN(nid(iob,id), id)
               sum(iob) = sum(iob) + &
	          cAir(iob,nz(iob)-j)* dsdh(iob,id,j)
           END DO
           ! double pass layers:
           DO j = MIN(nid(iob,id),id)+1,nid(iob,id)
              sum(iob) = sum(iob) + 2.0 * &
	           cAir(iob,nz(iob)-j)* dsdh(iob,id,j)
           END DO
        END IF
        scol(iob,nz(iob)-id) = sum(iob)
     END DO
  END DO

  END SUBROUTINE airmas


  SUBROUTINE swbiol(nw,wl,wc,j,s,label)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Create or read various weighting functions, e.g. biological action       =*
  !=  spectra, UV index, etc.                                                  =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of central wavelength of wavelength intervals    I)=*
  !=           in working wavelength grid                                      =*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  S      - REAL, value of each defined weighting function at each       (O)=*
  !=           defined wavelength                                              =*
  !=  LABEL  - CHARACTER*50, string identifier for each weighting function  (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER, INTENT(IN)               :: nw
  REAL, INTENT(IN)                  :: wl(kw)
  REAL, INTENT(IN)                  :: wc(kw)
  INTEGER, INTENT(INOUT)            :: j
  REAL, INTENT(INOUT)               :: s(ks,kw)
  CHARACTER (LEN=50), INTENT(INOUT) :: label(ks)

  INTEGER, PARAMETER :: kdata=1000

  ! internal:
  REAL :: x1(kdata)
  REAL :: y1(kdata)
  REAL :: yg(kw)

!LFR>   REAL :: fery, futr
  !EXTERNAL fery, futr
  INTEGER :: i, iw, n, iob

  INTEGER :: ierr

  INTEGER :: idum
  REAL :: dum1, dum2
  REAL :: em, a, b, c
  !REAL :: sum

  REAL :: a0, a1, a2, a3

  !_______________________________________________________________________

  !******** Photosynthetic Active Radiation (400 < PAR < 700 nm)
  ! conversion to micro moles m-2 s-1:
  !  s = s * (1e6/6.022142E23)(w/1e9)/(6.626068E-34*2.99792458E8)

  j = j + 1

  wbioStart=j
  label(j) = 'PAR, 400-700 nm, umol m-2 s-1'
  print *, j, label(j),wbioStart; call flush(6)
  DO iw = 1, nw-1
    IF (wc(iw) > 400. .AND. wc(iw) < 700.) THEN
      s(j,iw) = 8.36E-3 * wc(iw)
    ELSE
      s(j,iw) = 0.
    END IF
  END DO

  !********* unity raf constant slope:

  j = j + 1
  label(j) = 'Exponential decay, 14 nm/10'
  DO iw = 1, nw-1
    s(j,iw) = 10.**(-(wc(iw) -300.)/14.)
  END DO

  !*********** DNA damage action spectrum
  ! from: Setlow, R. B., The wavelengths in sunlight effective in
  !       producing skin cancer: a theoretical analysis, Proceedings
  !       of the National Academy of Science, 71, 3363 -3366, 1974.
  ! normalize to unity at 300 nm
  ! Data read from original hand-drawn plot by Setlow
  ! received from R. Setlow in May 1995
  ! data is per quantum (confirmed with R. Setlow in May 1995).
  ! Therefore must put on energy basis if irradiance is is energy
  ! (rather than quanta) units.

  j = j + 1
  label(j) = 'DNA damage, in vitro (Setlow, 1974)'
  OPEN(UNIT=kin,FILE=trim(files(21)%fileName),STATUS='old')
  DO i = 1, 11
    READ(kin,*)
  END DO
  n = 55
  DO i = 1, n
    READ(kin,*) x1(i), y1(i)
    y1(i) = y1(i) / 2.4E-02  *  x1(i)/300.
  END DO
  CLOSE (kin)

  CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
  CALL addpnt(x1,y1,kdata,n,          0.,y1(1))
  CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
  CALL addpnt(x1,y1,kdata,n,      1.e+38,   0.)
  CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  IF (ierr /= 0) THEN
    WRITE(*,*) ierr, label(j)
    STOP
  END IF

  DO iw = 1, nw-1
    s(j,iw) = yg(iw)
  END DO

  !******** skin cancer in mice,  Utrecht/Phildelphia study
  !from de Gruijl, F. R., H. J. C. M. Sterenborg, P. D. Forbes,
  !     R. E. Davies, C. Cole, G. Kelfkens, H. van Weelden, H. Slaper,
  !     and J. C. van der Leun, Wavelength dependence of skin cancer
  !     induction by ultraviolet irradiation of albino hairless mice,
  !     Cancer Res., 53, 53-60, 1993.
  ! normalize at 300 nm.

  j = j + 1
  label(j) = 'SCUP-mice (de Gruijl et al., 1993)'
  DO iw = 1, nw-1
    s(j,iw) =  futr(wc(iw)) / futr(300.)
  END DO

  !********** Utrecht/Philadelphia mice spectrum corrected for humans skin.
  ! From de Gruijl, F.R. and J. C. van der Leun, Estimate of the wavelength
  ! dependency of ultraviolet carcinogenesis and its relevance to the
  ! risk assessment of a stratospheric ozone depletion, Health Phys., 4,
  ! 317-323, 1994.

  j = j + 1
  label(j) = 'SCUP-human (de Gruijl and van der Leun, 1994)'
  OPEN(UNIT=kin,FILE=trim(files(22)%fileName),STATUS='old')
  n = 28
  DO i = 1, n
    READ(kin,*) x1(i), y1(i)
  END DO

  CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
  CALL addpnt(x1,y1,kdata,n,          0.,y1(1))
  CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
  CALL addpnt(x1,y1,kdata,n,      1.e+38,   0.)
  CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  IF (ierr /= 0) THEN
    WRITE(*,*) ierr, label(j)
    STOP
  END IF

  DO iw = 1, nw-1
    s(j,iw) = yg(iw)
  END DO
  CLOSE (kin)

  !**************** CIE standard human erythema action spectrum
  !from:
  ! McKinlay, A. F., and B. L. Diffey, A reference action spectrum for
  ! ultraviolet induced erythema in human skin, in Human Exposure to
  ! Ultraviolet Radiation: Risks and Regulations, W. R. Passchler
  ! and B. F. M. Bosnajokovic, (eds.), Elsevier, Amsterdam, 1987.

  j = j + 1
  label(j) = 'CIE human erythema (McKinlay and Diffey, 1987)'
  DO iw = 1, nw-1
    s(j,iw) = fery(wc(iw))
  END DO

  !**************** UV index (Canadian - WMO/WHO)
  ! from:
  ! Report of the WMO Meeting of experts on UV-B measurements, data quality
  ! and standardization of UV indices, World Meteorological Organization
  ! (WMO), report No. 95, Geneva, 1994.
  ! based on the CIE erythema weighting, multiplied by 40.

  j = j + 1
  label(j) = 'UV index (WMO, 1994)'
  DO iw = 1, nw-1
    s(j,iw) = 40. * fery(wc(iw))
  END DO

  !************ Human erythema - Anders et al.
  ! from:
  ! Anders, A., H.-J. Altheide, M. Knalmann, and H. Tronnier,
  ! Action spectrum for erythema in humands investigated with dye lasers,
  ! Photochem. and Photobiol., 61, 200-203, 1995.
  ! for skin types II and III, Units are J m-2.

  j = j + 1
  label(j) = 'Erythema, humans (Anders et al., 1995)'
  OPEN(UNIT=kin,FILE=trim(files(23)%fileName),STATUS='old')
  DO i = 1, 5
    READ(kin,*)
  END DO
  n = 28
  DO i = 1, n
    READ(kin,*) x1(i), y1(i)
    y1(i) = 1./y1(i)
  END DO

  CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
  CALL addpnt(x1,y1,kdata,n,          0.,y1(1))
  CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
  CALL addpnt(x1,y1,kdata,n,      1.e+38,   0.)
  CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  IF (ierr /= 0) THEN
    WRITE(*,*) ierr, label(j)
    STOP
  END IF

  DO iw = 1, nw-1
    s(j,iw) = yg(iw)
  END DO
  CLOSE (kin)

  !******** 1991-92 ACGIH threshold limit values
  ! from
  ! ACGIH, 1991-1992 Threshold Limit Values, American Conference
  !  of Governmental and Industrial Hygienists, 1992.

  j = j + 1
  label(j) = 'Occupational TLV (ACGIH, 1992)'
  OPEN(UNIT=kin,FILE=trim(files(24)%fileName),STATUS='old')
  n = 56
  DO i = 1, n
    READ(kin,*) x1(i), y1(i)
    y1(i) = y1(i)
  END DO

  CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
  CALL addpnt(x1,y1,kdata,n,          0.,y1(1))
  CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
  CALL addpnt(x1,y1,kdata,n,      1.e+38,   0.)
  CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  IF (ierr /= 0) THEN
    WRITE(*,*) ierr, label(j)
    STOP
  END IF

  DO iw = 1, nw-1
    s(j,iw) = yg(iw)
  END DO
  CLOSE (kin)

  !******** phytoplankton, Boucher et al. (1994)
  ! from Boucher, N., Prezelin, B.B., Evens, T., Jovine, R., Kroon, B., Moline, M.A.,
  ! and Schofield, O., Icecolors '93: Biological weighting function for the ultraviolet
  !  inhibition  of carbon fixation in a natural antarctic phytoplankton community,
  ! Antarctic Journal, Review 1994, pp. 272-275, 1994.
  ! In original paper, value of b and m (em below are given as positive.  Correct values
  ! are negative. Also, limit to positive values.

  j = j + 1
  label(j) = 'Phytoplankton (Boucher et al., 1994)'
  a = 112.5
  b = -6.223E-01
  c = 7.670E-04
  em = -3.17E-06
  DO iw = 1, nw-1
    IF (wc(iw) > 290. .AND. wc(iw) < 400.) THEN
      s(j,iw) = em + EXP(a+b*wc(iw)+c*wc(iw)*wc(iw))
    ELSE
      s(j,iw) = 0.
    END IF
    s(j,iw) = MAX(s(j,iw),0.)
  END DO

  !******** phytoplankton, Cullen et al.
  ! Cullen, J.J., Neale, P.J., and Lesser, M.P., Biological weighting function for the
  !  inhibition of phytoplankton photosynthesis by ultraviolet radiation, Science, 25,
  !  646-649, 1992.
  ! phaeo

  j = j + 1
  label(j) = 'Phytoplankton, phaeo (Cullen et al., 1992)'
  OPEN(UNIT=kin,FILE=trim(files(25)%fileName),STATUS='old')
  n = 106
  DO i = 1, n
    READ(kin,*) idum, dum1, dum2, y1(i)
    x1(i) = (dum1+dum2)/2.
  END DO

  CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
  CALL addpnt(x1,y1,kdata,n,          0.,y1(1))
  CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
  CALL addpnt(x1,y1,kdata,n,      1.e+38,   0.)
  CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  IF (ierr /= 0) THEN
    WRITE(*,*) ierr, label(j)
    STOP
  END IF

  DO iw = 1, nw-1
    s(j,iw) = yg(iw)
  END DO
  CLOSE(kin)

  ! proro

  j = j + 1
  label(j) = 'Phytoplankton, proro (Cullen et al., 1992)'
  OPEN(UNIT=kin,FILE=trim(files(26)%fileName),STATUS='old')
  n = 100
  DO i = 1, n
    READ(kin,*) idum, dum1, dum2, y1(i)
    x1(i) = (dum1+dum2)/2.
  END DO

  CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
  CALL addpnt(x1,y1,kdata,n,          0.,y1(1))
  CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
  CALL addpnt(x1,y1,kdata,n,      1.e+38,   0.)
  CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  IF (ierr /= 0) THEN
    WRITE(*,*) ierr, label(j)
    STOP
  END IF

  DO iw = 1, nw-1
    s(j,iw) = yg(iw)
  END DO
  CLOSE (kin)

  !*** Damage to lens of pig eyes, from
  ! Oriowo, M. et al. (2001). Action spectrum for in vitro
  ! UV-induced cataract using whole lenses. Invest. Ophthalmol. & Vis. Sci. 42,
  ! 2596-2602.  For pig eyes. Last two columns computed by L.O.Bjorn.

  j = j + 1
  label(j) = 'Cataract, pig (Oriowo et al., 2001)'
  OPEN(UNIT=kin,FILE=trim(files(27)%fileName),STATUS='old')
  DO i = 1, 7
    READ(kin,*)
  END DO
  n = 18
  DO i = 1, n
    READ(kin,*) x1(i), dum1, dum1, y1(i)
  END DO

  ! extrapolation to 400 nm (has very little effect on raf):
  !      do i = 1, 30
  !         n = n + 1
  !         x1(n) = x1(n-1) + 1.
  !         y1(n) = 10**(5.7666 - 0.0254*x1(n))
  !      enddo

  CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
  CALL addpnt(x1,y1,kdata,n,          0.,y1(1))
  CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
  CALL addpnt(x1,y1,kdata,n,      1.e+38,   0.)
  CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  IF (ierr /= 0) THEN
    WRITE(*,*) ierr, label(j)
    STOP
  END IF

  DO iw = 1, nw-1
    s(j,iw) = yg(iw)
  END DO
  CLOSE(kin)

  !***** Plant damage - Caldwell 1971
  !  Caldwell, M. M., Solar ultraviolet radiation and the growth and
  ! development of higher plants, Photophysiology 6:131-177, 1971.

  j = j + 1
  label(j) = 'Plant damage (Caldwell, 1971)'

  ! Fit to Caldwell (1971) data by
  ! Green, A. E. S., T. Sawada, and E. P. Shettle, The middle
  ! ultraviolet reaching the ground, Photochem. Photobiol., 19,
  ! 251-259, 1974.

  DO iw = 1, nw-1
    s(j,iw) = 2.628*(1. - (wc(iw)/313.3)**2)* EXP(-(wc(iw)-300.)/31.08)
    IF( s(j,iw) < 0. .OR. wc(iw) > 313.) THEN
      s(j,iw) = 0.
    END IF
  END DO

  ! Alternative fit to Caldwell (1971) by
  ! Micheletti, M. I. and R. D. Piacentini, Photochem. Photobiol.,
  ! 76, pp.?, 2002.

  a0 = 570.25
  a1 = -4.70144
  a2 = 0.01274
  a3 = -1.13118E-5
  DO iw = 1, nw-1
    s(j,iw) = a0 + a1*wc(iw) + a2*wc(iw)**2  + a3*wc(iw)**3
    IF( s(j,iw) < 0. .OR. wc(iw) > 313.) THEN
      s(j,iw) = 0.
    END IF
  END DO

  !***** Plant damage - Flint & Caldwell 2003
  !  Flint, S. D. and M. M. Caldwell, A biological spectral weigthing
  !  function for ozone depletion research with higher plants, Physiologia
  !  Plantorum, in press, 2003.
  !  Data available to 366 nm

  j = j + 1
  label(j) = 'Plant damage (Flint & Caldwell, 2003)'

  DO iw = 1, nw-1
    s(j,iw) = EXP( 4.688272*EXP( -EXP(0.1703411*(wc(iw)-307.867)/1.15))+  &
        ((390-wc(iw))/121.7557-4.183832) )

  ! put on per joule (rather than per quantum) basis:

    s(j,iw) = s(j,iw) * wc(iw)/300.

    IF( s(j,iw) < 0. .OR. wc(iw) > 366.) THEN
      s(j,iw) = 0.
    END IF

  END DO
  wbioEnd=j
  !***************************************************************
  !***************************************************************

  !_______________________________________________________________________

  IF (j > ks) STOP '1001'
  !_______________________________________________________________________

  END SUBROUTINE swbiol

  SUBROUTINE swphys(nw,wl,wc,j,s,label)

  !-----------------------------------------------------------------------------*
  !=  PURPOSE:                                                                 =*
  !=  Create or read various spectral weighting functions, physically-based    =*
  !=  e.g. UV-B, UV-A, visible ranges, instrument responses, etc.              =*
  !-----------------------------------------------------------------------------*
  !=  PARAMETERS:                                                              =*
  !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
  !=           wavelength grid                                                 =*
  !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
  !=           working wavelength grid                                         =*
  !=  WC     - REAL, vector of central wavelength of wavelength intervals    I)=*
  !=           in working wavelength grid                                      =*
  !=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
  !=  S      - REAL, value of each defined weighting function at each       (O)=*
  !=           defined wavelength                                              =*
  !=  LABEL  - CHARACTER*40, string identifier for each weighting function  (O)=*
  !=           defined                                                         =*
  !-----------------------------------------------------------------------------*
  IMPLICIT NONE

  INTEGER, INTENT(IN)             :: nw
  REAL, INTENT(IN)                :: wl(kw)
  REAL, INTENT(IN)                :: wc(kw)
  INTEGER, INTENT(OUT)            :: j
  REAL, INTENT(OUT)               :: s(ks,kw)
  CHARACTER (LEN=50), INTENT(OUT) :: label(ks)

  INTEGER, PARAMETER :: kdata=1000

  ! internal:
  REAL :: x1(kdata)
  REAL :: y1(kdata)
  REAL :: yg(kw)

  INTEGER :: i, iw, n

  INTEGER :: ierr

  !INTEGER :: idum
  !REAL :: dum1, dum2
  !REAL :: em, a, b, c
  REAL :: sum

  !_______________________________________________________________________

  j = 0

  !******** UV-B (280-315 nm)

  j = j + 1
  label(j) = 'UV-B, 280-315 nm'
  DO iw = 1, nw-1
    IF (wc(iw) > 280. .AND. wc(iw) < 315.) THEN
      s(j,iw) = 1.
    ELSE
      s(j,iw) = 0.
    END IF
  END DO

  !******** UV-B* (280-320 nm)

  j = j + 1
  label(j) = 'UV-B*, 280-320 nm'
  DO iw = 1, nw-1
    IF (wc(iw) > 280. .AND. wc(iw) < 320.) THEN
      s(j,iw) = 1.
    ELSE
      s(j,iw) = 0.
    END IF
  END DO

  !******** UV-A (315-400 nm)

  j = j + 1
  label(j) = 'UV-A, 315-400 nm'
  DO iw = 1, nw-1
    IF (wc(iw) > 315. .AND. wc(iw) < 400.) THEN
      s(j,iw) = 1.
    ELSE
      s(j,iw) = 0.
    END IF
  END DO

  !******** visible+ (> 400 nm)

  j = j + 1
  label(j) = 'vis+, > 400 nm'
  DO iw = 1, nw-1
    IF (wc(iw) > 400.) THEN
      s(j,iw) = 1.
    ELSE
      s(j,iw) = 0.
    END IF
  END DO

  !*********  Gaussian transmission functions

  j = j + 1
  label(j) = 'Gaussian, 305 nm, 10 nm FWHM'
  sum = 0.
  DO iw = 1, nw-1
    !srf -avoiding floating-point exception for single precision
    if( ( LOG(2.) * ((wc(iw)-305.)/(5.))**2)  < 80.) then
     s(j,iw) = EXP(- ( LOG(2.) * ((wc(iw)-305.)/(5.))**2) )
    else
      s(j,iw) = 0.0
    endif

    sum = sum + s(j,iw)
  END DO
  DO iw = 1, nw-1
    s(j,iw) = s(j,iw)/sum
  END DO

  j = j + 1
  label(j) = 'Gaussian, 320 nm, 10 nm FWHM'
  sum = 0.
  DO iw = 1, nw-1
    if( ( LOG(2.) * ((wc(iw)-320.)/(5.))**2) < 80.) then
      s(j,iw) = EXP(- ( LOG(2.) * ((wc(iw)-320.)/(5.))**2) )
    else
      s(j,iw) = 0.0
    endif

    sum = sum + s(j,iw)
  END DO
  DO iw = 1, nw-1
    s(j,iw) = s(j,iw)/sum
  END DO

  j = j + 1
  label(j) = 'Gaussian, 340 nm, 10 nm FWHM'
  sum = 0.
  DO iw = 1, nw-1
    if( (LOG(2.) * ((wc(iw)-340.)/(5.))**2)< 80.) then
      s(j,iw) = EXP(- ( LOG(2.) * ((wc(iw)-340.)/(5.))**2) )
    else
      s(j,iw) = 0.0
    endif
    sum = sum + s(j,iw)
  END DO
  DO iw = 1, nw-1
    s(j,iw) = s(j,iw)/sum
  END DO

  j = j + 1
  label(j) = 'Gaussian, 380 nm, 10 nm FWHM'
  sum = 0.
  DO iw = 1, nw-1
    if(( LOG(2.) * ((wc(iw)-380.)/(5.))**2)< 80.) then
      s(j,iw) = EXP(- ( LOG(2.) * ((wc(iw)-380.)/(5.))**2) )
    else
      s(j,iw) = 0.0
    endif
    sum = sum + s(j,iw)
  END DO
  DO iw = 1, nw-1
    s(j,iw) = s(j,iw)/sum
  END DO

  !********* RB Meter, model 501
  !  private communication, M. Morys (Solar Light Co.), 1994.
  ! From: morys@omni.voicenet.com (Marian Morys)
  ! Received: from acd.ucar.edu by sasha.acd.ucar.edu (AIX 3.2/UCB 5.64/4.03)
  !          id AA17274; Wed, 21 Sep 1994 11:35:44 -0600

  j = j + 1
  label(j) = 'RB Meter, model 501'
  OPEN(UNIT=kin,FILE=trim(files(20)%fileName),STATUS='old')
  n = 57
  DO i = 1, n
    READ(kin,*) x1(i), y1(i)
  END DO

  CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
  CALL addpnt(x1,y1,kdata,n,          0.,y1(1))
  CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
  CALL addpnt(x1,y1,kdata,n,      1.e+38,   0.)
  CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  IF (ierr /= 0) THEN
    WRITE(*,*) ierr, label(j)
    STOP
  END IF

  DO iw = 1, nw-1
    s(j,iw) = yg(iw)
  END DO
  CLOSE (kin)

  !***************************************************************
  !***************************************************************

  !_______________________________________________________________________

  IF (j > ks) STOP '1001'
  !_______________________________________________________________________

  END SUBROUTINE swphys

  !_______________________________________________________________________
  REAL FUNCTION refrac(w,airden)
  ! input vacuum wavelength, nm and air density, molec cm-3
  ! output refractive index for standard air
  ! (dry air at 15 deg. C, 101.325 kPa, 0.03% CO2)


  REAL, INTENT(IN) :: w
  REAL, INTENT(IN) :: airden


  ! internal

  REAL :: sig,  dum

  ! from CRC Handbook, originally from Edlen, B., Metrologia, 2, 71, 1966.
  ! valid from 200 nm to 2000 nm
  ! beyond this range, use constant value

  sig = 1.e3/w

  IF (w < 200.) sig = 1.e3/200.
  IF (w > 2000.) sig = 1.e3/2000.

  dum = 8342.13 + 2406030./(130. - sig*sig) + 15997./(38.9 - sig*sig)

  ! adjust to local air density

  dum = dum * airden/(2.69E19 * 273.15/288.15)

  ! index of refraction:

  refrac = 1. + 1.e-8 * dum

  END FUNCTION refrac

  !_______________________________________________________________________
  SUBROUTINE wshift(mrefr, n, w, airden)

  ! Shift wavelength scale between air and vacuum.
  ! if mrefr = 1, shift input waveelengths in air to vacuum.
  ! if mrefr = -1, shift input wavelengths from vacuum to air
  ! if any other number, don't shift

  INTEGER, INTENT(IN)   :: n
  INTEGER, INTENT(IN)   :: mrefr
  REAL, INTENT(IN)      :: airden
  REAL, INTENT(INOUT)   :: w(n)

  ! internal
  INTEGER :: i

  IF(mrefr == 1) THEN
    DO i = 1, n
      w(i) = w(i) / refrac(w(i),airden)
    END DO
  ELSE IF(mrefr == -1) THEN
    DO i = 1, n
      w(i) = w(i) * refrac(w(i),airden)
    END DO
  END IF

  END SUBROUTINE wshift

  SUBROUTINE readpol(coef)
    !-----------------------------------------------------------------------------*
    !=  PARAMETERS:  (XUEXI)                                                     =*
    !=  coef(mz,ms,mp) - REAL, pol coefficent for FTUV                           =*
    !=  mz = 5 zenith anagle 0,20,40,60,80				       =*
    !=  ms = 73 species
    !=  mp = 5 pol coeff  0,1,2,3,4
    IMPLICIT NONE

    REAL, INTENT(OUT) :: coef(mz,ms,mp)

    INTEGER :: iz,is,ip
    OPEN(81, FILE=trim(files(  1)%fileName))

    DO iz=1,mz
      READ(81,*)
      DO is=1,ms
        READ(81,*)
        DO ip=1,mp
          READ(81,*) coef(iz,is,ip)
        END DO
      END DO
    END DO

  END SUBROUTINE readpol
  SUBROUTINE ChangeKData(kData)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: kData

      IF(ALLOCATED(x1)) DEALLOCATE(x1)
      IF(ALLOCATED(y1)) DEALLOCATE(y1)
      IF(ALLOCATED(x2)) DEALLOCATE(x2)
      IF(ALLOCATED(y2)) DEALLOCATE(y2)
      IF(ALLOCATED(x3)) DEALLOCATE(x3)
      IF(ALLOCATED(y3)) DEALLOCATE(y3)
      IF(ALLOCATED(x4)) DEALLOCATE(x4)
      IF(ALLOCATED(y4)) DEALLOCATE(y4)
      IF(ALLOCATED(x5)) DEALLOCATE(x5)
      IF(ALLOCATED(y5)) DEALLOCATE(y5)
      IF(ALLOCATED(x)) DEALLOCATE(x)
      IF(ALLOCATED(y)) DEALLOCATE(y)

      ALLOCATE(x1(kData))
      ALLOCATE(y1(kData))
      ALLOCATE(x2(kData))
      ALLOCATE(y2(kData))
      ALLOCATE(x3(kData))
      ALLOCATE(y3(kData))
      ALLOCATE(x4(kData))
      ALLOCATE(y4(kData))
      ALLOCATE(x5(kData))
      ALLOCATE(y5(kData))
      ALLOCATE(x (kData))
      ALLOCATE(y (kData))

  END SUBROUTINE ChangeKData

  SUBROUTINE ReadAll(nw,wl)
    IMPLICIT NONE

    !-----------------------------------------------------------------------------*
    !=  PURPOSE:								 =*
    !=  Read the O3 cross section				 =*
    !=  Combined data from WMO 85 Ozone Assessment (use 273K value from 	 =*
    !=  175.439-847.5 nm) and:  						 =*
    !=  For Hartley and Huggins bands, use temperature-dependent values from	 =*
    !=  Molina, L. T., and M. J. Molina, Absolute absorption cross sections	 =*
    !=  of ozone in the 185- to 350-nm wavelength range, J. Geophys. Res.,	 =*
    !=  vol. 91, 14501-14508, 1986.						 =*
    !-----------------------------------------------------------------------------*
    !=  PARAMETERS:								 =*
    !=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
    !=  	 wavelength grid						 =*
    !=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
    !=  	 working wavelength grid					 =*
    !=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
    !=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
    !=  XS     - REAL, cross section (cm^2) for O3			      (O)=*
    !=  	 at each defined wavelength and each defined altitude level	 =*
    !-----------------------------------------------------------------------------*


    INTEGER, INTENT(IN)  :: nw
    !INTEGER, INTENT(IN)  :: nz
    REAL, INTENT(IN)	 :: wl(kw)

    INTEGER :: iw
    !INTEGER :: iz


    INTEGER :: kdata
    INTEGER :: nReact

    INTEGER :: i, n, k, idum
    REAL    :: a1, a2, dum

    INTEGER :: ierr,irev
    INTEGER :: n1, n2, n3, n4, n5
    REAL    :: yglocal(kw)
    !REAL :: qyVal
    INTEGER :: irow, icol
    REAL :: xs270
    REAL :: xs280
    CHARACTER (LEN=120) :: inline

    !Reading xreference file
    OPEN(UNIT=kin,FILE=trim(files(145)%fileName),STATUS='old')
    !PRINT *,'Reading crossreference file for reactions'
    DO i=1,3
       READ(kin,*)
    END DO
    DO i=1,kj
       READ(kin,FMT='(A57,I3.3,1X,L1)') inline,xRef(i),doReaction(i)
    END DO
    CLOSE(UNIT=kin)


    kData = 250
    CALL ChangeKData(kData)
    SELECT CASE (mOption(1))
      CASE(1)
  	!----------------------------------------------------------
  	! cross sections from WMO 1985 Ozone Assessment
  	! from 175.439 to 847.500 nm
  	! use value at 273 K
  	!PRINT *,'Reading cross sections from WMO 1985 Ozone Assessment'
	!PRINT *,'LFR->'//trim(files(139)%fileName)//'!'
  	OPEN(UNIT=kin,FILE=trim(files(139)%fileName),STATUS='old')
  	DO i = 1, 3
  	  READ(kin,*)
  	END DO
  	n = 158
  	DO i = 1, n
  	  READ(kin,*) idum, a1, a2, dum, dum, dum, dum, y1(i)
  	  x1(i) = (a1+a2)/2.
  	END DO
  	CLOSE (kin)

  	CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
  	CALL addpnt(x1,y1,kdata,n,		 0.,0.)
  	CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
  	CALL addpnt(x1,y1,kdata,n,	     1.e+38,0.)
  	CALL inter2(nw,wl,yglocal,n,x1,y1,ierr)
  	IF (ierr /= 0) THEN
  	  WRITE(*,*) ierr, 'O3 cross section - WMO'
  	  STOP
  	END IF

  	DO iw = 1, nw-1
  	  mm_o3xs(iw) = yglocal(iw)
  	END DO

  	! For Hartley and Huggins bands, use temperature-dependent values from
  	! Molina, L. T., and M. J. Molina, Absolute absorption cross sections
  	! of ozone in the 185- to 350-nm wavelength range,
  	! J. Geophys. Res., vol. 91, 14501-14508, 1986.

  	OPEN(UNIT=kin,FILE=trim(files(140)%fileName),STATUS='old')
  	DO i = 1, 5
  	  READ(kin,*)
  	END DO
  	n1 = 220
  	n2 = 220
  	n3 = 220
  	DO i = 1, n1
  	  READ(kin,*) x1(i), y1(i), y2(i), y3(i)
  	  x2(i) = x1(i)
  	  x3(i) = x1(i)
  	END DO
  	CLOSE (kin)

  	CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
  	CALL addpnt(x1,y1,kdata,n1,		  0.,0.)
  	CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
  	CALL addpnt(x1,y1,kdata,n1,	       1.e+38,0.)
  	CALL inter2(nw,wl,yglocal,n1,x1,y1,ierr)
  	IF (ierr /= 0) THEN
  	  WRITE(*,*) ierr, 'O3 xsect - 226K Molina'
  	  STOP
  	END IF
  	DO iw = 1, nw-1
  	  s226(iw) = yglocal(iw)*1.e-20
  	END DO

  	CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
  	CALL addpnt(x2,y2,kdata,n2,		  0.,0.)
  	CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
  	CALL addpnt(x2,y2,kdata,n2,	       1.e+38,0.)
  	CALL inter2(nw,wl,yglocal,n2,x2,y2,ierr)
  	IF (ierr /= 0) THEN
  	  WRITE(*,*) ierr, 'O3 xsect - 263K Molina'
  	  STOP
  	END IF
  	DO iw = 1, nw-1
  	  s263(iw) = yglocal(iw)*1.e-20
  	END DO

  	CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
  	CALL addpnt(x3,y3,kdata,n3,		  0.,0.)
  	CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
  	CALL addpnt(x3,y3,kdata,n3,	       1.e+38,0.)
  	CALL inter2(nw,wl,yglocal,n3,x3,y3,ierr)
  	IF (ierr /= 0) THEN
  	  WRITE(*,*) ierr, 'O3 xsect - 298K Molina'
  	  STOP
  	END IF
  	DO iw = 1, nw-1
  	  s298(iw) = yglocal(iw)*1.e-20
  	END DO
      CASE (2)
  	!----------------------------------------------------------
  	! cross sections from WMO 1985 Ozone Assessment
  	! from 175.439 to 847.500 nm
  	! use value at 273 K
  	PRINT *, 'Reading cross sections from WMO 1985 Ozone Assessment'
  	OPEN(UNIT=kin,FILE=trim(files(141)%fileName),STATUS='old')
  	DO i = 1, 3
  	  READ(kin,*)
  	END DO
  	n = 158
  	DO i = 1, n
  	  READ(kin,*) idum, a1, a2, dum, dum, dum, dum, y1(i)
  	  x1(i) = (a1+a2)/2.
  	END DO
  	CLOSE (kin)

  	CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
  	CALL addpnt(x1,y1,kdata,n,	       0.,0.)
  	CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
  	CALL addpnt(x1,y1,kdata,n,	   1.e+38,0.)
  	CALL inter2(nw,wl,yglocal,n,x1,y1,ierr)
  	IF (ierr /= 0) THEN
  	  WRITE(*,*) ierr, 'O3 cross section - WMO'
  	  STOP
  	END IF

  	DO iw = 1, nw-1
  	  mm_o3xs(iw) = yglocal(iw)
  	END DO

  	!=  For Hartley and Huggins bands, use temperature-dependent values from     =*
  	!=  Malicet et al., J. Atmos. Chem.  v.21, pp.263-273, 1995.		   =*

  	OPEN(UNIT=kin,FILE=trim(files(142)%fileName),STATUS='old')
  	DO i = 1, 1
  	  READ(kin,*)
  	END DO
  	n1 = 15001
  	n2 = 15001
  	n3 = 15001
  	n4 = 15001

  	DO i = 1, n1
  	  READ(kin,*) x1(i), y1(i), y2(i), y3(i), y4(i)
  	  x2(i) = x1(i)
  	  x3(i) = x1(i)
  	  x4(i) = x1(i)
  	END DO
  	CLOSE (kin)

  	CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
  	CALL addpnt(x1,y1,kdata,n1,		0.,0.)
  	CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
  	CALL addpnt(x1,y1,kdata,n1,	     1.e+38,0.)
  	CALL inter2(nw,wl,yglocal,n1,x1,y1,ierr)
  	IF (ierr /= 0) THEN
  	  WRITE(*,*) ierr, 'O3 xsect - 295K Malicet'
  	  STOP
  	END IF
  	DO iw = 1, nw-1
  	  s295(iw) = yglocal(iw)
  	END DO

  	CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
  	CALL addpnt(x2,y2,kdata,n2,		0.,0.)
  	CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
  	CALL addpnt(x2,y2,kdata,n2,	     1.e+38,0.)
  	CALL inter2(nw,wl,yglocal,n2,x2,y2,ierr)
  	IF (ierr /= 0) THEN
  	  WRITE(*,*) ierr, 'O3 xsect - 243K Malicet'
  	  STOP
  	END IF
  	DO iw = 1, nw-1
  	  s243(iw) = yglocal(iw)
  	END DO

  	CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
  	CALL addpnt(x3,y3,kdata,n3,		0.,0.)
  	CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
  	CALL addpnt(x3,y3,kdata,n3,	     1.e+38,0.)
  	CALL inter2(nw,wl,yglocal,n3,x3,y3,ierr)
  	IF (ierr /= 0) THEN
  	  WRITE(*,*) ierr, 'O3 xsect - 228K Malicet'
  	  STOP
  	END IF
  	DO iw = 1, nw-1
  	  s228(iw) = yglocal(iw)
  	END DO

  	CALL addpnt(x4,y4,kdata,n4,x4(1)*(1.-deltax),0.)
  	CALL addpnt(x4,y4,kdata,n4,		0.,0.)
  	CALL addpnt(x4,y4,kdata,n4,x4(n4)*(1.+deltax),0.)
  	CALL addpnt(x4,y4,kdata,n4,	     1.e+38,0.)
  	CALL inter2(nw,wl,yglocal,n4,x4,y4,ierr)
  	IF (ierr /= 0) THEN
  	  WRITE(*,*) ierr, 'O3 xsect - 218K Malicet'
  	  STOP
  	END IF
  	DO iw = 1, nw-1
  	  s218(iw) = yglocal(iw)
  	END DO

      CASE (3)
    	!----------------------------------------------------------
    	! cross sections from WMO 1985 Ozone Assessment
    	! from 175.439 to 847.500 nm
    	! use value at 273 K

    	OPEN(UNIT=kin,FILE=trim(files(143)%fileName),STATUS='old')
    	DO i = 1, 3
    	  READ(kin,*)
    	END DO
    	n = 158
    	DO i = 1, n
    	  READ(kin,*) idum, a1, a2, dum, dum, dum, dum, y1(i)
    	  x1(i) = (a1+a2)/2.
    	END DO
    	CLOSE (kin)

    	CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
    	CALL addpnt(x1,y1,kdata,n,	       0.,0.)
    	CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
    	CALL addpnt(x1,y1,kdata,n,	   1.e+38,0.)
    	CALL inter2(nw,wl,yglocal,n,x1,y1,ierr)
    	IF (ierr /= 0) THEN
    	  WRITE(*,*) ierr, 'O3 cross section - WMO'
    	  STOP
    	END IF

    	DO iw = 1, nw-1
    	  mm_o3xs(iw) = yglocal(iw)
    	END DO

    	! For Hartley and Huggins bands, use temperature-dependent values from
    	!  Bass et al.

    	OPEN(UNIT=kin,FILE=trim(files(144)%fileName),STATUS='old')
    	DO i = 1, 8
    	  READ(kin,*)
    	END DO
    	n1 = 1915
    	n2 = 1915
    	n3 = 1915
    	DO i = 1, n1
    	  READ(kin,*) x1(i), y1(i), y2(i), y3(i)
    	  x2(i) = x1(i)
    	  x3(i) = x1(i)
    	END DO
    	CLOSE (kin)

    	CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
    	CALL addpnt(x1,y1,kdata,n1,		0.,0.)
    	CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
    	CALL addpnt(x1,y1,kdata,n1,	     1.e+38,0.)
    	CALL inter2(nw,wl,yglocal,n1,x1,y1,ierr)
    	IF (ierr /= 0) THEN
    	  WRITE(*,*) ierr, 'O3 xsect - c0 Bass'
    	  STOP
    	END IF
    	DO iw = 1, nw-1
    	  c0(iw) = yglocal(iw)
    	END DO

    	CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
    	CALL addpnt(x2,y2,kdata,n2,		0.,0.)
    	CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
    	CALL addpnt(x2,y2,kdata,n2,	     1.e+38,0.)
    	CALL inter2(nw,wl,yglocal,n2,x2,y2,ierr)
    	IF (ierr /= 0) THEN
    	  WRITE(*,*) ierr, 'O3 xsect - c1 Bass'
    	  STOP
    	END IF
    	DO iw = 1, nw-1
    	  c1(iw) = yglocal(iw)
    	END DO

    	CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
    	CALL addpnt(x3,y3,kdata,n3,		0.,0.)
    	CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
    	CALL addpnt(x3,y3,kdata,n3,	     1.e+38,0.)
    	CALL inter2(nw,wl,yglocal,n3,x3,y3,ierr)
    	IF (ierr /= 0) THEN
    	  WRITE(*,*) ierr, 'O3 xsect - c2 Bass'
    	  STOP
    	END IF
        DO iw = 1, nw-1
          c2(iw) = yglocal(iw)
        END DO

      END SELECT

      !===============================================================
      !R01
      !===============================================================
      ! read parameters from JPL'97
      nReact=01
      kData = 500
      CALL ChangeKData(kData)
      SELECT CASE(mOption(2))
         CASE(kjpl97)
            OPEN(UNIT=kin,FILE=trim(files(28)%fileName),STATUS='old')
            READ(kin,*)
            READ(kin,*)
            READ(kin,*)
            n1 = 21
            n2 = n1
            DO i = 1, n1
              READ(kin,*) x1(i), y1(i), y2(i)
              x2(i) = x1(i)
            END DO
            CLOSE(kin)

            CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),y1(1))
            CALL addpnt(x1,y1,kdata,n1, 	    0.,y1(1))
            CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),y1(n1))
            CALL addpnt(x1,y1,kdata,n1, 	 1.e+38,y1(n1))
            CALL inter2(nw,wl,yg1(:,nReact),n1,x1,y1,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr,'Reading R01' ! jlabel(j)
              STOP
            END IF

            CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
      	    CALL addpnt(x2,y2,kdata,n2, 	    0.,y2(1))
            CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),y2(n2))
            CALL addpnt(x2,y2,kdata,n2, 	 1.e+38,y2(n2))
            CALL inter2(nw,wl,yg2(:,nReact),n2,x2,y2,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr,'Reading R01'!, jlabel(j)
              STOP
            END IF
    !       read parameters from Michelsen, H. A., R.J. Salawitch, P. O. Wennber,
    !       and J. G. Anderson, Geophys. Res. Lett., 21, 2227-2230, 1994.
        CASE(kmich)
            OPEN(UNIT=kin,FILE=trim(files(29)%fileName),STATUS='old')
            READ(kin,*)
            READ(kin,*)
            READ(kin,*)
            n1 = 21
            n2 = n1
            DO i = 1, n1
              READ(kin,*) x1(i), y1(i), y2(i)
              x2(i) = x1(i)
            END DO
            CLOSE(kin)

            CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),y1(1))
            CALL addpnt(x1,y1,kdata,n1, 	    0.,y1(1))
            CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),y1(n1))
            CALL addpnt(x1,y1,kdata,n1, 	 1.e+38,y1(n1))
            CALL inter2(nw,wl,yg1(:,nReact),n1,x1,y1,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr,'Reading R01'!, jlabel(j)
              STOP
            END IF

            CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
            CALL addpnt(x2,y2,kdata,n2, 	    0.,y2(1))
            CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),y2(n2))
            CALL addpnt(x2,y2,kdata,n2, 	 1.e+38,y2(n2))
            CALL inter2(nw,wl,yg2(:,nReact),n2,x2,y2,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr,'Reading R01'!, jlabel(j)
              STOP
            END IF
        ! quantum yield data from
        ! Shetter et al, J.Geophys.Res., v 101 (D9), pg. 14,631-14,641, June 20, 1996
	CASE(kshet)
    	   OPEN(UNIT=kin,FILE=trim(files(30)%fileName),STATUS='OLD')
    	   READ(kin,*) idum, n
    	   DO i = 1, idum-2
    	      READ(kin,*)
    	   END DO
    	   n = n-2
    	   DO i = 1, n
    	       READ(kin,*) x1(i),y3(i),y4(i),y1(i),y2(i)
    	      x2(i) = x1(i)
    	      x3(i) = x1(i)
    	      x4(i) = x1(i)
    	   END DO
    	   DO i = n+1, n+2
    	     READ(kin,*) x3(i),y3(i),y4(i)
    	     x4(i) = x3(i)
    	   END DO
    	   CLOSE(kin)

    	   n1 = n
    	   n2 = n
    	   n3 = n+2
    	   n4 = n+2

    	   ! coefficients for exponential fit:

    	   CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax), y1(1))
    	   CALL addpnt(x1,y1,kdata,n1,		   0., y1(1))
    	   CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
    	   CALL addpnt(x1,y1,kdata,n1,		 1E38,0.)

    	   CALL inter2(nw,wl,yg1(:,nReact), n1,x1,y1, ierr)
    	   IF (ierr /= 0) THEN
    	     WRITE(*,*) ierr,'Reading R01'!, jlabel(j)
    	     STOP
    	   END IF

    	   CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
    	   CALL addpnt(x2,y2,kdata,n2,		  0.,y2(1))
    	   CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
    	   CALL addpnt(x2,y2,kdata,n2,		 1E38,0.)

    	   CALL inter2(nw,wl,yg2(:,nReact), n2,x2,y2, ierr)
    	   IF (ierr /= 0) THEN
    	     WRITE(*,*) ierr,'Reading R01'!, jlabel(j)
    	     STOP
    	   END IF

    	   ! phi data at 298 and 230 K

    	   CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),y3(1))
    	   CALL addpnt(x3,y3,kdata,n3,		  0.,y3(1))
    	   CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
    	   CALL addpnt(x3,y3,kdata,n3,		 1E38,0.)

    	   CALL inter2(nw,wl,yg3(:,nReact), n3,x3,y3, ierr)
    	   IF (ierr /= 0) THEN
    	     WRITE(*,*) ierr,'Reading R01'!,jlabel(j)
    	     STOP
    	   END IF

    	   CALL addpnt(x4,y4,kdata,n4,x4(1)*(1.-deltax),y4(1))
    	   CALL addpnt(x4,y4,kdata,n4,		  0.,y4(1))
    	   CALL addpnt(x4,y4,kdata,n4,x4(n4)*(1.+deltax),0.)
    	   CALL addpnt(x4,y4,kdata,n4,		 1E38,0.)

    	   CALL inter2(nw,wl,yg4(:,nReact), n4,x4,y4, ierr)
    	   IF (ierr /= 0) THEN
    	     WRITE(*,*) ierr,'Reading R01'!,jlabel(j)
    	     STOP
    	   END IF
      END SELECT

      !================================================================
      !R02
      !================================================================
      nReact=02
      kData = 200
      CALL ChangeKData(kData)
      ! cross section
      !------------NEED TO CHANGE kdata = 1000 FOR DAVIDSON ET AL. DATA---------
      ! measurements by:
      ! Davidson, J. A., C. A. Cantrell, A. H. McDaniel, R. E. Shetter,
      ! S. Madronich, and J. G. Calvert, Visible-ultraviolet absorption
      ! cross sections for NO2 as a function of temperature, J. Geophys.
      ! Res., 93, 7105-7112, 1988.
      !     from 263.8 to 648.8 nm in approximately 0.5 nm intervals
      !     OPEN(UNIT=kin,FILE='DATAE1/NO2/NO2_ncar_00.abs',STATUS='old')
      !     n = 750
      !     DO i = 1, n
      !        READ(kin,*) x1(i), y1(i), dum, dum, idum
      !     ENDDO
      !     CLOSE(kin)

      !     CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      !     CALL addpnt(x1,y1,kdata,n,               0.,0.)
      !     CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      !     CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      !     CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      !     IF (ierr .NE. 0) THEN
      !        WRITE(*,*) ierr, jlabel(j)
      !        STOP
      !     ENDIF

      ! cross section data from JPL 94 recommendation
      ! JPL 97 recommendation is identical
      SELECT CASE (mOption(3))
         CASE (1)
            OPEN(UNIT=kin,FILE=trim(files(31)%fileName),STATUS='old')
            READ(kin,*) idum, n
            DO i = 1, idum-2
              READ(kin,*)
            END DO
            ! read in wavelength bins, cross section at T0 and temperature correction
            ! coefficient a;  see input file for details.
            ! data need to be scaled to total area per bin so that they can be used with
            ! inter3

            DO i = 1, n
              READ(kin,*) x1(i), x3(i), y1(i), dum, y2(i)
              y1(i) = (x3(i)-x1(i)) * y1(i)*1.e-20
              y2(i) = (x3(i)-x1(i)) * y2(i)*1.e-22
              x2(i) = x1(i)
            END DO
            CLOSE(kin)

            x1(n+1) = x3(n)
            x2(n+1) = x3(n)
            n = n+1
            n1 = n

            CALL inter3(nw,wl,yg1(:,nReact),n,x1,y1,0)
            CALL inter3(nw,wl,yg2(:,nReact),n1,x2,y2,0)

            ! yg1(:,nReact), yg2(:,nReact) are per nm, so rescale by bin widths

            DO iw = 1, nw-1
              yg1(iw,nReact) = yg1(iw,nReact)/(wl(iw+1)-wl(iw))
              yg2(iw,nReact) = yg2(iw,nReact)/(wl(iw+1)-wl(iw))
            END DO

        CASE (2)
            OPEN(UNIT=kin,FILE=trim(files(32)%fileName),STATUS='old')
            DO i = 1, 9
              READ(kin,*)
            END DO
            n = 135
            DO i = 1, n
              READ(kin,*) idum, y1(i)
              x1(i) = FLOAT(idum)
            END DO

            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
            CALL addpnt(x1,y1,kdata,n,  	     0.,y1(1))
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
            CALL addpnt(x1,y1,kdata,n,  	 1.e+38,   0.)
            CALL inter2(nw,wl,yg1(:,nReact),n,x1,y1,ierr)
      END SELECT

      ! quantum yield
      ! from Gardiner, Sperry, and Calvert
      OPEN(UNIT=kin,FILE=trim(files(33)%fileName),STATUS='old')
      DO i = 1, 8
 	 READ(kin,*)
      END DO
      n = 66
      DO i = 1, n
 	 READ(kin,*) x1(i),y1(i)
      END DO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,	       0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
      CALL addpnt(x1,y1,kdata,n,	   1.e+38,   0.)
      CALL inter2(nw,wl,yg1n(:,nReact),n,x1,y1,ierr)
      IF (ierr /= 0) THEN
 	 WRITE(*,*) ierr,'Reading R02'!, jlabel(j)
 	 STOP
      END IF

      !================================================================
      !R03
      !================================================================
      nReact=03
      kData = 350
      CALL ChangeKData(kData)

      ! cross section
      !     measurements of Graham and Johnston 1978
      OPEN(UNIT=kin,FILE=trim(files(34)%fileName),STATUS='old')
      DO i = 1, 9
        READ(kin,*)
      END DO
      n = 305
      DO irow = 1, 30
        READ(kin,*) ( y1(10*(irow-1) + icol), icol =  1, 10 )
      END DO
      READ(kin,*) ( y1(300 + icol), icol = 1, 5 )
      CLOSE (kin)
      DO i = 1, n
        y1(i) =  y1(i) * 1.e-19
        x1(i) = 400. + 1.*FLOAT(i-1)
      END DO

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	       0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	   1.e+38,0.)
      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
      IF (ierr /= 0) THEN
        WRITE(*,*) ierr,'Reading R03'!, jlabel(j)
        STOP
      END IF

      !     cross section from JPL94:
      OPEN(UNIT=kin,FILE=trim(files(35)%fileName),STATUS='old')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      END DO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i)*1E-20
      END DO
      CLOSE (kin)
      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	       0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	   1.e+38,0.)
      CALL inter2(nw,wl,yg1(:,nReact),n,x1,y1,ierr)
      IF (ierr /= 0) THEN
        WRITE(*,*) ierr,'Reading R03'!, jlabel(j)
        STOP
      END IF

      ! use JPL94 for wavelengths longer than 600 nm
      DO iw = 1, nw-1
         IF(wl(iw) > 600.) yg(iw,nReact) = yg1(iw,nReact)
      END DO

      !================================================================
      !R04
      !================================================================
      nReact=04
      kData = 100
      CALL ChangeKData(kData)
      ! cross section from jpl97, table up to 280 nm

      OPEN(UNIT=kin,FILE=trim(files(36)%fileName),STATUS='old')
      READ(kin,*) idum, n
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      DO i = 1,  n
  	READ(kin,*) x1(i), y1(i)
  	y1(i) = y1(i) * 1.e-20
      END DO
      xs270 = y1(n-2)
      xs280 = y1(n)

      CLOSE(kin)

      CALL addpnt(x1,y1,kdata, n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata, n,		0.,0.)
      CALL addpnt(x1,y1,kdata, n,x1(n)*(1.+deltax),y1(n))
      CALL addpnt(x1,y1,kdata, n,	     1.e36,y1(n))

      CALL inter2(nw,wl,yg(:,nReact), n,x1,y1, ierr)
      IF (ierr /= 0) THEN
  	WRITE(0,*) ierr,'Reading R04'!,jlabel(j)
  	STOP
      END IF

      !================================================================
      !R05
      !================================================================
      nReact=05
      kData = 100
      CALL ChangeKData(kData)
      OPEN(UNIT=kin,FILE=trim(files(37)%fileName),STATUS='old')
      DO i = 1, 13
  	READ(kin,*)
      END DO
      n = 91
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i)
  	y1(i) = y1(i) * 1.e-20
      END DO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	       0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	   1.e+38,0.)
      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr,'Reading R05'!, jlabel(j)
  	STOP
      END IF

      !================================================================
      !R06
      !================================================================
      nReact=06
      kData = 100
      CALL ChangeKData(kData)
      !* cross section from JPL85

      !      OPEN(UNIT=kin,FILE='dataj1/abs/HNO3.abs',STATUS='old')
      !      DO i = 1, 9
      ! 	READ(kin,*)
      !      ENDDO
      !      n = 29
      !      DO i = 1, n
      ! 	READ(kin,*) x1(i), y1(i)
      ! 	y1(i) = y1(i) * 1.E-20
      !      ENDDO
      !      CLOSE (kin)

      !      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      !      CALL addpnt(x1,y1,kdata,n, 	      0.,0.)
      !      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      !      CALL addpnt(x1,y1,kdata,n, 	  1.e+38,0.)
      !      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      !      IF (ierr .NE. 0) THEN
      ! 	WRITE(*,*) ierr, jlabel(j)
      ! 	STOP
      !      ENDIF

      !* quantum yield = 1

      !      qy = 1.
      !      DO iw = 1, nw - 1
      ! 	DO i = 1, nz
      ! 	   sq(j,i,iw) = yg(iw)*qy
      ! 	ENDDO
      !      ENDDO


      ! HNO3 cross section parameters from Burkholder et al. 1993

      OPEN(UNIT=kin,FILE=trim(files(38)%fileName),STATUS='old')
      DO i = 1, 6
        READ(kin,*)
      END DO
      n1 =  83
      n2 = n1
      DO i = 1, n1
        READ(kin,*) y1(i), y2(i)
        x1(i) = 184. + i*2.
        x2(i) = x1(i)
      END DO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,		0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,	     1.e+38,0.)
      CALL inter2(nw,wl,yg1(:,nReact),n1,x1,y1,ierr)
      IF (ierr /= 0) THEN
        WRITE(*,*) ierr,'Reading R06'!, jlabel(j)
        STOP
      END IF

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
      CALL addpnt(x2,y2,kdata,n2,		0.,y2(1))
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),y2(n2))
      CALL addpnt(x2,y2,kdata,n2,	     1.e+38,y2(n2))
      CALL inter2(nw,wl,yg2(:,nReact),n2,x2,y2,ierr)
      IF (ierr /= 0) THEN
        WRITE(*,*) ierr,'Reading R06'!, jlabel(j)
        STOP
      END IF

      !================================================================
      !R07
      !================================================================
      nReact=07
      kData = 100
      CALL ChangeKData(kData)
      OPEN(UNIT=kin,FILE=trim(files(39)%fileName),STATUS='old')
      DO i = 1, 4
 	READ(kin,*)
      END DO
      n = 31
      DO i = 1, n
 	READ(kin,*) x1(i), y1(i)
 	y1(i) = y1(i) * 1.e-20
      END DO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	       0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	   1.e+38,0.)
      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
      IF (ierr /= 0) THEN
 	WRITE(*,*) ierr,'Reading R07'!, jlabel(j)
 	STOP
      END IF

      !================================================================
      !R08
      !================================================================
      nReact=08
      kData = 100
      CALL ChangeKData(kData)
      !     OPEN(UNIT=kin,FILE='dataj1/abs/H2O2_lin.abs',STATUS='old')
      !     DO i = 1, 7
      !        READ(kin,*)
      !     ENDDO
      !     n = 32
      !     DO i = 1, n
      !        READ(kin,*) x1(i), y1(i)
      !        y1(i) = y1(i) * 1.E-20
      !     ENDDO
      !     CLOSE (kin)

      !      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      !      CALL addpnt(x1,y1,kdata,n, 	      0.,0.)
      !      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      !      CALL addpnt(x1,y1,kdata,n, 	  1.e+38,0.)
      !      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      !      IF (ierr .NE. 0) THEN
      ! 	WRITE(*,*) ierr, jlabel(j)
      ! 	STOP
      !      ENDIF

      ! cross section from JPL94 (identical to JPL97)
      ! tabulated data up to 260 nm
      OPEN(UNIT=kin,FILE=trim(files(40)%filename),STATUS='old')
      READ(kin,*) idum,n
      DO i = 1, idum-2
        READ(kin,*)
      END DO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1.e-20
      END DO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	       0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	   1.e+38,0.)
      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
      IF (ierr /= 0) THEN
        WRITE(*,*) ierr,'Reading R08'!, jlabel(j)
        STOP
      END IF

      !================================================================
      !R09
      !================================================================
      nReact=09
      kData = 100
      CALL ChangeKData(kData)
      SELECT CASE (mOption(4))
         CASE (1)
            OPEN(UNIT=kin,FILE=trim(files(41)%fileName),STATUS='old')
            DO i = 1, 5
              READ(kin,*)
            END DO

            n5 = 25
            n4 = 27
            n3 = 29
            n2 = 31
            n1 = 39
            DO i = 1, n5
              READ(kin,*) x1(i), y1(i), y2(i), y3(i), y4(i), y5(i)
            END DO
            DO i = n5 + 1, n4
      	      READ(kin,*) x1(i), y1(i), y2(i), y3(i), y4(i)
            END DO
            DO i = n4 + 1, n3
              READ(kin,*) x1(i), y1(i), y2(i), y3(i)
            END DO
            DO i = n3 + 1, n2
               READ(kin,*) x1(i), y1(i), y2(i)
            END DO
            DO i = n2 + 1, n1
              READ(kin,*) x1(i), y1(i)
            END DO
      	    CLOSE (kin)

            DO i = 1, n1
              y1(i) = y1(i) * 1.e-23
            END DO
            DO i = 1, n2
              x2(i) = x1(i)
              y2(i) = y2(i) * 1.e-23
            END DO
            DO i = 1, n3
      	      x3(i) = x1(i)
              y3(i) = y3(i) * 1.e-23
            END DO
            DO i = 1, n4
              x4(i) = x1(i)
              y4(i) = y4(i) * 1.e-23
            END DO
            DO i = 1, n5
              x5(i) = x1(i)
              y5(i) = y5(i) * 1.e-23
      	    END DO

            CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),y1(1))
            CALL addpnt(x1,y1,kdata,n1, 	    0.,y1(1))
            CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n1, 	1.e+38,0.)
            CALL inter2(nw,wl,yg1(:,nReact),n1,x1,y1,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr,'Reading R09'!, jlabel(j)
              STOP
      	    END IF

            CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
            CALL addpnt(x2,y2,kdata,n2, 	    0.,y2(1))
            CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
            CALL addpnt(x2,y2,kdata,n2, 	1.e+38,0.)
            CALL inter2(nw,wl,yg2(:,nReact),n2,x2,y2,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr,'Reading R09'!, jlabel(j)
              STOP
      	    END IF

            CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),y3(1))
            CALL addpnt(x3,y3,kdata,n3, 	    0.,y3(1))
            CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
            CALL addpnt(x3,y3,kdata,n3, 	1.e+38,0.)
            CALL inter2(nw,wl,yg3(:,nReact),n3,x3,y3,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr,'Reading R09'!, jlabel(j)

      	    END IF

            CALL addpnt(x4,y4,kdata,n4,x4(1)*(1.-deltax),y4(1))
            CALL addpnt(x4,y4,kdata,n4, 	    0.,y4(1))
            CALL addpnt(x4,y4,kdata,n4,x4(n4)*(1.+deltax),0.)
            CALL addpnt(x4,y4,kdata,n4, 	1.e+38,0.)
            CALL inter2(nw,wl,yg4(:,nReact),n4,x4,y4,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr,'Reading R09'!, jlabel(j)
              STOP
      	    END IF

            CALL addpnt(x5,y5,kdata,n5,x5(1)*(1.-deltax),y5(1))
            CALL addpnt(x5,y5,kdata,n5, 	    0.,y5(1))
            CALL addpnt(x5,y5,kdata,n5,x5(n5)*(1.+deltax),0.)
            CALL addpnt(x5,y5,kdata,n5, 	1.e+38,0.)
   	    CALL inter2(nw,wl,yg5(:,nReact),n5,x5,y5,ierr)
   	    IF (ierr /= 0) THEN
   	      WRITE(*,*) ierr,'Reading R09'!, jlabel(j)
   	      STOP
   	    END IF
         ! jpl97, with temperature dependence formula,
   	 !w = 290 nm to 340 nm,
   	 !T = 210K to 300 K
         !sigma, cm2 = exp((0.06183-0.000241*w)*(273.-T)-(2.376+0.14757*w))
         CASE (2)
            OPEN(UNIT=kin,FILE=trim(files(42)%fileName),STATUS='old')
            DO i = 1, 6
               READ(kin,*)
            END DO
            n1 = 87
            DO i = 1, n1
               READ(kin,*) x1(i), y1(i)
               y1(i) = y1(i) * 1.e-20
            END DO
            CLOSE(kin)

            CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),y1(1))
            CALL addpnt(x1,y1,kdata,n1, 	      0.,y1(1))
            CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n1, 	  1.e+38,0.)
            CALL inter2(nw,wl,yg1(:,nReact),n1,x1,y1,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr,'Reading R09'!, jlabel(j)
              STOP
            END IF
      END SELECT

      !================================================================
      !R10
      !================================================================
      nReact=10
      kData = 16000
      CALL ChangeKData(kData)
      SELECT CASE(mOption(5))
         CASE(1)
            ! read NBS/Bass data
            OPEN(UNIT=kin,FILE=trim(files(43)%fileName),STATUS='old')
            n = 4032
            DO i = 1, n
              READ(kin,*) x(i), y(i)
            END DO
            CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
            CALL addpnt(x,y,kdata,n,		   0.,0.)
            CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
            CALL addpnt(x,y,kdata,n,	       1.e+38,0.)

            CALL inter2(nw,wl,yg1(:,nReact),n,x,y,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading 10'
              STOP
            END IF
         CASE(2:4)
            OPEN(UNIT=kin,FILE=trim(files(44)%fileName),STATUS='old')
            DO i = 1, 4
              READ(kin,*)
            END DO
            n = 121
            DO i = 1, n
              READ(kin,*) x(i), y(i)
              y(i) = y(i) * 1.e-20
            END DO
            CLOSE(kin)
            CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
      	    CALL addpnt(x,y,kdata,n,		 0.,0.)
            CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
            CALL addpnt(x,y,kdata,n,	     1.e+38,0.)
            CALL inter2(nw,wl,yg1(:,nReact),n,x,y,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading 10'
              STOP
            END IF
       END SELECT
       SELECT CASE(mOption(5))
	 CASE (3)
    	  ! data are on wavenumber grid (cm-1), so convert to wavelength in nm:
    	  ! grid was on increasing wavenumbers, so need to reverse to get increasing
    	  ! wavelengths
    	  ! cross section assumed to be zero for wavelengths longer than 360 nm
    	  ! if y1 < 0, then make = 0 (some negative cross sections, actually 273 K intercepts
    	  ! are in the original data,  Here, make equal to zero)
    	      OPEN(kin,FILE=trim(files(45)%fileName),STATUS='old')
    	      READ(kin,*) idum, n
    	      DO i = 1, idum-2
    	      READ(kin,*)
    	      END DO
    	      DO i = 1, n
    	         READ(kin,*) x1(i), y1(i), y2(i)
    	         x1(i) = 1./x1(i) * 1E7
    	         IF (x1(i) > 360.) THEN
    	  	    y1(i) = 0.
    	            y2(i) = 0.
    	         END IF
              END DO
              CLOSE(kin)

              DO i = 1, n/2
    	         irev = n+1-i
    	         dum = x1(i)
    	         x1(i) = x1(irev)
    	         x1(irev) = dum
    	         dum = y1(i)
    	         y1(i) = y1(irev)
    	         y1(irev) = dum
    	         dum = y2(i)
    	         y2(i) = y2(irev)
    	         y2(irev) = dum
              END DO
              DO i = 1, n
    	        x2(i) = x1(i)
    	        y1(i) = MAX(y1(i),0.)
              END DO
              n1 = n
              n2 = n

              CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
              CALL addpnt(x1,y1,kdata,n1,	      0.,0.)
              CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
              CALL addpnt(x1,y1,kdata,n1,	    1E38,0.)
              CALL inter2(nw,wl,yg2(:,nReact),n1,x1,y1,ierr)
              IF (ierr /= 0) THEN
    	         WRITE(*,*) ierr,'Reading R10'
    	         STOP
              END IF

              CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
              CALL addpnt(x2,y2,kdata,n2,		0.,0.)
              CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
              CALL addpnt(x2,y2,kdata,n2,	       1E38,0.)
              CALL inter2(nw,wl,yg3(:,nReact),n2,x2,y2,ierr)
              IF (ierr /= 0) THEN
    	         WRITE(*,*) ierr,'Reading R10'
    	         STOP
              END IF
           CASE (4)
              OPEN(UNIT=kin,FILE=trim(files(46)%fileName), STATUS='old')
              DO i = 1, 4
    	         READ(kin,*)
              END DO
              n = 23
              DO i = 1, n
    	         READ(kin,*) x2(i), y2(i), y3(i), dum, dum
    	         x3(i) = x2(i)
              END DO
              CLOSE(kin)
              n2 = n
              n3 = n

              CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
              CALL addpnt(x2,y2,kdata,n2,	      0.,0.)
              CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
              CALL addpnt(x2,y2,kdata,n2,	    1E38,0.)
              CALL inter2(nw,wl,yg2(:,nReact),n2,x2,y2,ierr)
              IF (ierr /= 0) THEN
    	      WRITE(*,*) ierr,'Reading R10'
    	      STOP
              END IF

              CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
              CALL addpnt(x3,y3,kdata,n3,	      0.,0.)
              CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
              CALL addpnt(x3,y3,kdata,n3,	     1E38,0.)
              CALL inter2(nw,wl,yg3(:,nReact),n3,x3,y3,ierr)
              IF (ierr /= 0) THEN
    	         WRITE(*,*) ierr,'Reading R10'
    	         STOP
              END IF
           CASE (5)
              ! read Rodgers data
              OPEN(UNIT=kin,FILE=trim(files(47)%fileName) ,STATUS='old')
              DO i = 1, 10
                  READ(kin,*)
              END DO
              n = 261
              DO i = 1, n
                  READ(kin,*) x(i), y(i), dum
                  y(i) = y(i) * 1.e-20
              END DO
              CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
              CALL addpnt(x,y,kdata,n,		   0.,0.)
              CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
              CALL addpnt(x,y,kdata,n,		 1.e+38,0.)
              CALL inter2(nw,wl,yg1(:,nReact),n,x,y,ierr)
              IF (ierr /= 0) THEN
                WRITE(*,*) ierr, 'Reading 10'
                STOP
              END IF
           CASE (6)
              OPEN(UNIT=kin,FILE=trim(files(48)%fileName),STATUS='old')
              DO i = 1, 3
                READ(kin,*)
              END DO
              n = 126
              DO i = 1, n
                READ(kin,*) x2(i), y2(i), y3(i)
                x3(i) = x2(i)
              END DO
              CLOSE(kin)
              n2 = n
              n3 = n

              CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
              CALL addpnt(x2,y2,kdata,n2, 	      0.,0.)
              CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
              CALL addpnt(x2,y2,kdata,n2, 	    1E38,0.)
              CALL inter2(nw,wl,yg2(:,nReact),n2,x2,y2,ierr)
              IF (ierr /= 0) THEN
              WRITE(*,*) ierr,'Reading R10'
                 STOP
              END IF

              CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
              CALL addpnt(x3,y3,kdata,n3, 	      0.,0.)
              CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
              CALL addpnt(x3,y3,kdata,n3, 	     1E38,0.)
              CALL inter2(nw,wl,yg3(:,nReact),n3,x3,y3,ierr)
              IF (ierr /= 0) THEN
                 WRITE(*,*) ierr,'Reading R10'
                 STOP
              END IF
        END SELECT
        ! quantum yield
        SELECT CASE (mOption(6))
           CASE (1)
              OPEN(UNIT=kin,FILE=trim(files(49)%fileName),STATUS='old')
              DO i = 1, 11
                 READ(kin,*)
              END DO
              n = 20
              DO i = 1, n
                 READ(kin,*) x(i), y(i)
              END DO
              CLOSE(kin)
              CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),y(1))
              CALL addpnt(x,y,kdata,n,		   0.,y(1))
              CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
              CALL addpnt(x,y,kdata,n,	       1.e+38,0.)
              CALL inter2(nw,wl,yg4(:,nReact),n,x,y,ierr)
              IF (ierr /= 0) THEN
                WRITE(*,*) ierr, 'Reading 10'
                STOP
              END IF
              OPEN(UNIT=kin,FILE=trim(files(50)%fileName),STATUS='old')
              DO i = 1, 9
                READ(kin,*)
              END DO
              n = 33
              DO i = 1, n
                READ(kin,*) x(i), y(i)
              END DO
              CLOSE(kin)
              CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),y(1))
              CALL addpnt(x,y,kdata,n,		   0.,y(1))
              CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
              CALL addpnt(x,y,kdata,n,		 1.e+38,0.)
              CALL inter2(nw,wl,yg5(:,nReact),n,x,y,ierr)
              IF (ierr /= 0) THEN
                 WRITE(*,*) ierr,'Reading R10'
                 STOP
              END IF
           CASE (2)
              OPEN(UNIT=kin,FILE=trim(files(51)%fileName),STATUS='old')
  	      DO i = 1, 7
  		 READ(kin,*)
              END DO
  	      n = 13
  	      DO i = 1, n
  		READ(kin,*) x1(i), y1(i), y2(i)
  		x2(i) = x1(i)
  	      END DO
  	      CLOSE(kin)
  	      n1 = n
  	      n2 = n

  	      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),y1(1))
  	      CALL addpnt(x1,y1,kdata,n1,	      0.,y1(1))
  	      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
  	      CALL addpnt(x1,y1,kdata,n1,	  1.e+38,0.)
  	      CALL inter2(nw,wl,yg4(:,nReact),n1,x1,y1,ierr)
  	      IF (ierr /= 0) THEN
  		WRITE(*,*) ierr,'Reading R10'
  		STOP
  	      END IF

  	      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
  	      CALL addpnt(x2,y2,kdata,n2,	      0.,y2(1))
  	      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
  	      CALL addpnt(x2,y2,kdata,n2,	  1.e+38,0.)
  	      CALL inter2(nw,wl,yg5(:,nReact),n2,x2,y2,ierr)
  	      IF (ierr /= 0) THEN
  		WRITE(*,*) ierr,'Reading R10'
  		STOP
  	      END IF
	   CASE (3)
              OPEN(UNIT=kin,FILE=trim(files(52)%fileName),STATUS='old')
              DO i = 1, 4
                READ(kin,*)
              END DO
              n = 23
              DO i = 1, n
                READ(kin,*) x1(i), dum, dum, dum, dum, y1(i), y2(i)
                x2(i) = x1(i)
              END DO
              CLOSE(kin)
              n1 = n
              n2 = n

              CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),y1(1))
              CALL addpnt(x1,y1,kdata,n1, 	      0.,y1(1))
              CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
              CALL addpnt(x1,y1,kdata,n1, 	  1.e+38,0.)
              CALL inter2(nw,wl,yg4(:,nReact),n1,x1,y1,ierr)
              IF (ierr /= 0) THEN
                WRITE(*,*) ierr,'Reading R10'
                STOP
              END IF

              CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
              CALL addpnt(x2,y2,kdata,n2, 	      0.,y2(1))
              CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
              CALL addpnt(x2,y2,kdata,n2, 	  1.e+38,0.)
              CALL inter2(nw,wl,yg5(:,nReact),n2,x2,y2,ierr)
              IF (ierr /= 0) THEN
                WRITE(*,*) ierr,'Reading R10'
                STOP
              END IF
        END SELECT

      !================================================================
      !R11
      !================================================================
      nReact=11
      kData = 150
      CALL ChangeKData(kData)
      SELECT CASE(mOption(7))
         CASE (1)
              OPEN(UNIT=kin,FILE=trim(files(53)%fileName),STATUS='old')
              DO i = 1, 4
                READ(kin,*)
              END DO
              n = 106
              DO i = 1, n
                READ(kin,*) x1(i), y1(i)
                y1(i) = y1(i) * 1.e-20
              END DO
              CLOSE(kin)

              CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
              CALL addpnt(x1,y1,kdata,n,  	     0.,0.)
              CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
              CALL addpnt(x1,y1,kdata,n,  	 1.e+38,0.)
              CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
              IF (ierr /= 0) THEN
                WRITE(*,*) ierr, 'Reading R11'
                STOP
              END IF
         CASE (2)
              ! cross section from Calvert and  Pitts
              OPEN(UNIT=kin,FILE=trim(files(54)%fileName),STATUS='old')
              DO i = 1, 14
                READ(kin,*)
              END DO
              n = 54
              DO i = 1, n
                READ(kin,*) x1(i), y1(i)
                x1(i) = x1(i)/10.
                y1(i) = y1(i) * 3.82E-21
              END DO
              CLOSE (kin)

              CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
              CALL addpnt(x1,y1,kdata,n,  	     0.,0.)
              CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
              CALL addpnt(x1,y1,kdata,n,  	 1.e+38,0.)
              CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
              IF (ierr /= 0) THEN
                WRITE(*,*) ierr, 'Reading R11'
                STOP
              END IF
          CASE (3)
              OPEN(UNIT=kin,FILE=trim(files(55)%fileName),STATUS='old')
              DO i = 1, 3
                READ(kin,*)
              END DO
              n = 106
              DO i = 1, n
                READ(kin,*) x1(i), y1(i)
              END DO
              CLOSE (kin)

              CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
              CALL addpnt(x1,y1,kdata,n,  	     0.,0.)
              CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
              CALL addpnt(x1,y1,kdata,n,  	 1.e+38,0.)
              CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
              IF (ierr /= 0) THEN
                WRITE(*,*) ierr, 'Reading R11'
                STOP
              END IF
          CASE (4)
    	      !cross section from KFA tables
    	      !ch3cho.001 - Calvert and Pitts 1966
    	      !ch3cho.002 - Meyrahn thesis 1984
    	      !ch3cho.003 - Schneider and Moortgat, priv comm. MPI Mainz 1989, 0.012 nm resol.
    	      !ch3cho.004 - Schneider and Moortgat, priv comm. MPI Mainz 1989, 0.08  nm resol.
    	      !ch3cho.005 - IUPAC'92
    	      !ch3cho.006 - Libuda, thesis Wuppertal 1992

    	      ! OPEN(UNIT=kin,FILE='dataj2/kfa/ch3cho.001',STATUS='old')
    	      ! n = 217
    	      ! OPEN(UNIT=kin,FILE='dataj2/kfa/ch3cho.002',STATUS='old')
    	      ! n = 63
    	      ! OPEN(UNIT=kin,FILE='dataj2/kfa/ch3cho.003',STATUS='old')
    	      ! n = 13738
    	      ! OPEN(UNIT=kin,FILE='dataj2/kfa/ch3cho.004',STATUS='old')
    	      ! n = 2053
    	      OPEN(UNIT=kin,FILE=trim(files(56)%fileName),STATUS='old')
    	      n = 18
    	      ! OPEN(UNIT=kin,FILE='dataj2/kfa/ch3cho.006',STATUS='old')
    	      ! n = 1705

    	      DO i = 1, n
    	        READ(kin,*) x1(i), y1(i)
    	      END DO
    	      CLOSE (kin)

    	      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
    	      CALL addpnt(x1,y1,kdata,n,  	     0.,0.)
    	      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
    	      CALL addpnt(x1,y1,kdata,n,  	 1.e+38,0.)
    	      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
    	      IF (ierr /= 0) THEN
                 WRITE(*,*) ierr, 'Reading R11'
                 STOP
              END IF
        END SELECT
        ! quantum yields
        SELECT CASE (mOption(8))
           CASE (1)
              OPEN(UNIT=kin,FILE=trim(files(57)%fileName),STATUS='old')
              DO i = 1, 4
                READ(kin,*)
              END DO
              n = 12
              DO i = 1, n
                READ(kin,*) x1(i), y2(i), y1(i)
                x2(i) = x1(i)
              END DO
              CLOSE(kin)
              n1 = n
              n2 = n

              CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
              CALL addpnt(x1,y1,kdata,n1, 	      0.,0.)
              CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
              CALL addpnt(x1,y1,kdata,n1, 	  1.e+38,0.)
              CALL inter2(nw,wl,yg1(:,nReact),n1,x1,y1,ierr)
              IF (ierr /= 0) THEN
                WRITE(*,*) ierr, 'Reading R11'
                STOP
              END IF

              CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
              CALL addpnt(x2,y2,kdata,n2, 	      0.,0.)
              CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
              CALL addpnt(x2,y2,kdata,n2, 	  1.e+38,0.)
              CALL inter2(nw,wl,yg2(:,nReact),n2,x2,y2,ierr)
              IF (ierr /= 0) THEN
                WRITE(*,*) ierr, 'Reading R11'
                STOP
              END IF

              DO iw = 1, nw-1
                 yg3(:,nReact) = 0.
              END DO
         CASE (2)
              OPEN(UNIT=kin,FILE=trim(files(58)%fileName),STATUS='old')
              DO i = 1, 18
                READ(kin,*)
              END DO
              n = 10
              DO i = 1, n
                READ(kin,*) x1(i), y1(i)
              END DO
              CLOSE (kin)

              CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
              CALL addpnt(x1,y1,kdata,n,  	     0.,y1(1))
              CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
              CALL addpnt(x1,y1,kdata,n,  	 1.e+38,0.)
              CALL inter2(nw,wl,yg1(:,nReact),n,x1,y1,ierr)
              IF (ierr /= 0) THEN
                WRITE(*,*) ierr, 'Reading R11'
                STOP
              END IF

              OPEN(UNIT=kin,FILE=trim(files(59)%fileName),STATUS='old')
              DO i = 1, 10
                READ(kin,*)
              END DO
              n = 9
              DO i = 1, n
                READ(kin,*) x1(i), y1(i)
              END DO
              CLOSE (kin)

              CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
              CALL addpnt(x1,y1,kdata,n,  	     0.,y1(1))
              CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
              CALL addpnt(x1,y1,kdata,n,  	 1.e+38,0.)
              CALL inter2(nw,wl,yg2(:,nReact),n,x1,y1,ierr)
              IF (ierr /= 0) THEN
                WRITE(*,*) ierr, 'Reading R11'
                STOP
              END IF

              OPEN(UNIT=kin,FILE=trim(files(60)%fileName),STATUS='old')
              DO i = 1, 10
                READ(kin,*)
              END DO
              n = 9
              DO i = 1, n
                READ(kin,*) x1(i), y1(i)
              END DO
              CLOSE (kin)

              CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
              CALL addpnt(x1,y1,kdata,n,  	     0.,y1(1))
              CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
              CALL addpnt(x1,y1,kdata,n,  	 1.e+38,0.)
              CALL inter2(nw,wl,yg3(:,nReact),n,x1,y1,ierr)
              IF (ierr /= 0) THEN
                WRITE(*,*) ierr, 'Reading R11'
                STOP
              END IF
       END SELECT
       !pressure-dependence parameters
       OPEN(UNIT=kin,FILE=trim(files(61)%fileName), STATUS='old')
       DO i = 1, 4
         READ(kin,*)
       END DO
       n = 5
       DO i = 1, n
         READ(kin,*) x1(i), dum, dum, y1(i)
       END DO
       CLOSE (kin)

       CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
       CALL addpnt(x1,y1,kdata,n,	      0.,0.)
       CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
       CALL addpnt(x1,y1,kdata,n,	  1.e+38,0.)
       CALL inter2(nw,wl,yg4(:,nReact),n,x1,y1,ierr)
       IF (ierr /= 0) THEN
         WRITE(*,*) ierr, 'Reading R11'
         STOP
       END IF

      !================================================================
      !R12
      !================================================================
      nReact=12
      kData = 150
      CALL ChangeKData(kData)

      SELECT CASE( mOption(9))
         CASE (1)
            OPEN(UNIT=kin,FILE=trim(files( 62)%fileName), STATUS='old')
            DO i = 1, 4
              READ(kin,*)
            END DO
            n = 106
            DO i = 1, n
              READ(kin,*) x1(i), y1(i)
              y1(i) = y1(i) * 1.e-20
            END DO
            CLOSE(kin)

            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,  	     0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,  	 1.e+38,0.)
            CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading 12'
              STOP
            END IF
         CASE (2)
    	    ! cross section from KFA tables
    	    ! c2h5cho.001 - Calvert and Pitts 1966

    	    OPEN(UNIT=kin,FILE=trim(files(63)%fileName),STATUS='old')
    	    n = 83

    	    DO i = 1, n
    	      READ(kin,*) x1(i), y1(i)
    	    END DO
    	    CLOSE (kin)

    	    CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
    	    CALL addpnt(x1,y1,kdata,n,  	   0.,0.)
    	    CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
    	    CALL addpnt(x1,y1,kdata,n,         1.e+38,0.)
    	    CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
    	    IF (ierr /= 0) THEN
    	      WRITE(*,*) ierr, 'Reading 12'
    	      STOP
    	    END IF
      END SELECT
      ! quantum yields
      SELECT CASE (mOption(10))
         CASE (1)
	    !PRINT *,'LFR->Abrindo "'//trim(files(64)%fileName)//'"'
            OPEN(UNIT=kin,FILE=trim(files(64)%fileName), STATUS='old')
	    !PRINT *,'LFR->Abriu "'//trim(files(64)%fileName)//'"'
            DO i = 1, 4
              READ(kin,*)
            END DO
            n = 5
            DO i = 1, n
              READ(kin,*) x1(i), y1(i)
            END DO
            CLOSE(kin)
            n1 = n

            CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n1, 	    0.,0.)
            CALL addpnt(x1,y1,kdata,n1,340.,0.)
            CALL addpnt(x1,y1,kdata,n1, 	1.e+38,0.)
              CALL inter2(nw,wl,yg1(:,nReact),n1,x1,y1,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading 12'
              STOP
            END IF
         CASE (2)
            STOP
      END SELECT

      !================================================================
      !R13
      !================================================================
      nReact=13
      kData = 500
      CALL ChangeKData(kData)

      SELECT CASE (mOption(11))
         CASE (1)
            OPEN(UNIT=kin,FILE=trim(files(65)%fileName), STATUS='old')
            DO i = 1, 4
              READ(kin,*)
            END DO
            n = 110
            DO i = 1, n
              READ(kin,*) x1(i), y1(i)
              y1(i) = y1(i) * 1.e-20
            END DO
            CLOSE(kin)


            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,  	   0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,         1.e+38,0.)
            CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading R13'
              STOP
            END IF
         CASE (2)
            ! cross section from KFA tables
            ! chocho.001 - Plum et al. 1983
            OPEN(UNIT=kin,FILE=trim(files(66)%fileName),STATUS='old')
            n = 219

            DO i = 1, n
              READ(kin,*) x1(i), y1(i)
            END DO
            CLOSE (kin)

            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,  	   0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,         1.e+38,0.)
            CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading R13'
              STOP
            END IF
         CASE (3)
            ! cross section from Orlando et la.
            ! Orlando, J. J.; G. S. Tyndall, 2001:  The atmospheric chemistry of the
            ! HC(O)CO radical. Int. J. Chem. Kinet., 33, 149-156.
            OPEN(UNIT=kin, FILE=trim(files(67)%fileName),STATUS='old')

            DO i = 1, 6
              READ(kin,*)
            END DO
            n = 481
            DO i = 1, n
              READ(kin,*) x1(i), y1(i)
            END DO
            CLOSE (kin)

            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,  	   0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,         1.e+38,0.)
            CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading R13'
              STOP
            END IF
         CASE (4)
            OPEN(UNIT=kin, FILE=trim(files(68)%fileName),STATUS='old')

            DO i = 1, 8
              READ(kin,*)
            END DO
            n = 270
            DO i = 1, n
              READ(kin,*) x1(i), y1(i)
              y1(i) = y1(i) * 1.e-20
            END DO
            CLOSE (kin)

            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,  	     0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,  	 1.e+38,0.)
            CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading R13'
              STOP
            END IF
      END SELECT

      !================================================================
      !R14
      !================================================================
      nReact=14
      kData = 500
      CALL ChangeKData(kData)

      SELECT CASE (mOption(13))
         CASE (1)
            OPEN(UNIT=kin,FILE=trim(files(69)%fileName), STATUS='old')
            DO i = 1, 4
              READ(kin,*)
            END DO
            n = 38
            DO i = 1, n
              READ(kin,*) x1(i), y1(i)
              y1(i) = y1(i) * 1.e-20
            END DO
            CLOSE(kin)

            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,  	   0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,         1.e+38,0.)
            CALL inter2(nw,wl,yg1(:,nReact),n,x1,y1,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading 14'
              STOP
            END IF

            OPEN(UNIT=kin,FILE=trim(files(70)%fileName), STATUS='old')
            DO i = 1, 4
              READ(kin,*)
            END DO
            n = 75
            DO i = 1, n
              READ(kin,*) x1(i), y1(i)
              y1(i) = y1(i) * 1.e-20
            END DO
            CLOSE(kin)

            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,  	   0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,         1.e+38,0.)
            CALL inter2(nw,wl,yg2(:,nReact),n,x1,y1,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading 14'
              STOP
            END IF

            DO iw = 1, nw-1
              IF(wc(iw) < 402.) THEN
              yg(iw,nReact) = yg1(iw,nReact)
              ELSE
              yg(iw,nReact) = yg2(iw,nReact)
              END IF
            END DO
         CASE (2)
            OPEN(UNIT=kin,FILE=trim(files(71)%fileName), STATUS='old')
            n = 271
            DO i = 1, n
              READ(kin,*) x1(i), y1(i)
            END DO
            CLOSE(kin)

            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,  	   0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,         1.e+38,0.)
            CALL inter2(nw,wl,yg(:,14),n,x1,y1,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading 14'
              STOP
            END IF
         CASE (3:6)
            !       cross section from KFA tables
            !       ch3cocho.001 - Plum et al. 1983
            !       ch3cocho.002 - Meller et al. 1991, 0.033 nm resolution
            !       ch3cocho.003 - Meller et al. 1991, 1.0   nm resolution
            !       ch3cocho.004 - Staffelbach et al. 1995
            SELECT CASE (mOption(13))
	       CASE (3)
                  OPEN(UNIT=kin,FILE=trim(files(72)%fileName),STATUS='old')
                  n = 136
               CASE (4)
                  OPEN(UNIT=kin,FILE=trim(files(73)%fileName),STATUS='old')
                  n = 8251
	       CASE (5)
                  OPEN(UNIT=kin,FILE=trim(files(74)%fileName),STATUS='old')
                  n = 275
	       CASE (6)
                  OPEN(UNIT=kin,FILE=trim(files(75)%fileName),STATUS='old')
                  n = 162
	    END SELECT

            DO i = 1, n
              READ(kin,*) x1(i), y1(i)
            END DO
            CLOSE (kin)

            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,  	   0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,         1.e+38,0.)
            CALL inter2(nw,wl,yg(:,14),n,x1,y1,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading 14'
              STOP
            END IF
         CASE (7)
            OPEN(UNIT=kin,FILE=trim(files(76)%fileName), STATUS='old')
            DO i = 1, 7
              READ(kin,*)
            END DO
            n = 55
            DO i = 1, n
              READ(kin,*) x(i), y(i)
              y(i) = y(i) * 1.e-20
            END DO
            CLOSE(kin)

            CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
            CALL addpnt(x,y,kdata,n,		 0.,0.)
            CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
            CALL addpnt(x,y,kdata,n,	     1.e+38,0.)
            CALL inter2(nw,wl,yg1(:,nReact),n,x,y,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading 14'
              STOP
            END IF


            OPEN(UNIT=kin, FILE=trim(files(77)%fileName),STATUS='old')
            DO i = 1, 6
              READ(kin,*)
            END DO
            n = 481
            DO i = 1, n
              READ(kin,*) x1(i), y1(i)
            END DO
            CLOSE (kin)

            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,  	   0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,         1.e+38,0.)
            CALL inter2(nw,wl,yg2(:,nReact),n,x1,y1,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading 14'
              STOP
            END IF

            DO iw = 1, nw-1
              yg(iw,nReact) = 0.5*(yg1(iw,nReact) + yg2(iw,nReact))
            END DO
      END SELECT
      ! quantum yields
      SELECT CASE (mOption(14))
         CASE (4)
            OPEN(UNIT=kin,FILE=trim(files(78)%fileName), STATUS='old')
            DO i = 1, 5
              READ(kin,*)
            END DO
            n = 5
            DO i = 1, n
              READ(kin,*) x1(i), y1(i), y2(i)
              x2(i) = x1(i)
            END DO
            CLOSE (kin)
            n1 = n
            n2 = n

            CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),1.)
            CALL addpnt(x1,y1,kdata,n1, 	    0.,1.)
            CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n1, 	1.e+38,0.)
            CALL inter2(nw,wl,yg1(:,nReact),n1,x1,y1,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading 14'
              STOP
            END IF

            CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),1.)
            CALL addpnt(x2,y2,kdata,n2, 	    0.,1.)
            CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
            CALL addpnt(x2,y2,kdata,n2, 	1.e+38,0.)
            CALL inter2(nw,wl,yg2(:,nReact),n2,x2,y2,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading 14'
              STOP
            END IF
      END SELECT

      !================================================================
      !R15
      !================================================================
      nReact=15
      kData = 150
      CALL ChangeKData(kData)

      SELECT CASE (mOption(15))
         CASE (1)
     	   OPEN(UNIT=kin,FILE=trim(files(79)%fileName), STATUS='old')
     	   DO i = 1, 6
     	     READ(kin,*)
     	   END DO
     	   n = 35
     	   DO i = 1, n
     	     READ(kin,*) x1(i), y1(i)
     	     y1(i) = y1(i) * 3.82E-21
     	   END DO
     	   CLOSE (kin)

     	   CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
     	   CALL addpnt(x1,y1,kdata,n,		  0.,0.)
     	   CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
     	   CALL addpnt(x1,y1,kdata,n,	      1.e+38,0.)
     	   CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
     	   IF (ierr /= 0) THEN
     	     WRITE(*,*) ierr, 'Reading R15'
     	     STOP
     	   END IF
     	CASE (2)
     	   OPEN(UNIT=kin,FILE=trim(files(80)%fileName), STATUS='old')
     	   DO i = 1, 4
     	     READ(kin,*)
     	   END DO
     	   n = 96
     	   DO i = 1, n
     	     READ(kin,*) x1(i), y1(i)
     	     y1(i) = y1(i) * 1.e-20
     	   END DO
     	   CLOSE (kin)

     	   CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
     	   CALL addpnt(x1,y1,kdata,n,		  0.,0.)
     	   CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
     	   CALL addpnt(x1,y1,kdata,n,	      1.e+38,0.)
     	   CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
     	   IF (ierr /= 0) THEN
     	     WRITE(*,*) ierr, 'Reading R15'
     	     STOP
     	   END IF
     	CASE (3)
     	   OPEN(UNIT=kin,FILE=trim(files(81)%fileName), STATUS='old')
     	   DO i = 1, 12
     	     READ(kin,*)
     	   END DO
     	   n = 135
     	   DO i = 1, n
     	     READ(kin,*) x1(i), y1(i), y2(i), y3(i)
     	     x2(i) = x1(i)
     	     x3(i) = x1(i)
     	   END DO
     	   CLOSE (kin)
     	   n1 = n
     	   n2 = n
     	   n3 = n

     	   CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
     	   CALL addpnt(x1,y1,kdata,n1,  	   0.,0.)
     	   CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
     	   CALL addpnt(x1,y1,kdata,n1,         1.e+38,0.)
     	   CALL inter2(nw,wl,yg(:,nReact),n1,x1,y1,ierr)
     	   IF (ierr /= 0) THEN
     	     WRITE(*,*) ierr, 'Reading R15'
     	     STOP
     	   END IF

     	   CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
     	   CALL addpnt(x2,y2,kdata,n2,  	   0.,0.)
     	   CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
     	   CALL addpnt(x2,y2,kdata,n2,         1.e+38,0.)
     	   CALL inter2(nw,wl,yg2(:,nReact),n2,x2,y2,ierr)
     	   IF (ierr /= 0) THEN
     	     WRITE(*,*) ierr, 'Reading R15'
     	     STOP
     	   END IF


     	   CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
     	   CALL addpnt(x3,y3,kdata,n3,  	   0.,0.)
     	   CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
     	   CALL addpnt(x3,y3,kdata,n3,         1.e+38,0.)
     	   CALL inter2(nw,wl,yg3(:,nReact),n3,x3,y3,ierr)
     	   IF (ierr /= 0) THEN
     	     WRITE(*,*) ierr, 'Reading R15'
     	     STOP
     	   END IF
     END SELECT
     SELECT CASE(mOption(16))
        CASE (2)
     	   OPEN(UNIT=kin,FILE=trim(files(82)%fileName), STATUS='old')
     	   DO i = 1, 4
     	     READ(kin,*)
     	   END DO
     	   n = 9
     	   DO i = 1, n
     	     READ(kin,*) x1(i), y1(i)
     	   END DO
     	   CLOSE (kin)

     	   CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
     	   CALL addpnt(x1,y1,kdata,n,		  0.,0.)
     	   CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
     	   CALL addpnt(x1,y1,kdata,n,	      1.e+38,0.)
     	   CALL inter2(nw,wl,yg1(:,nReact),n,x1,y1,ierr)
     	   IF (ierr /= 0) THEN
     	     WRITE(*,*) ierr, 'Reading R15'
     	     STOP
     	   END IF
     END SELECT

      !================================================================
      !R16
      !================================================================
      nReact=16
      kData = 100
      CALL ChangeKData(kData)

      SELECT CASE (mOption(17))
         CASE (1)
            !	    OPEN(UNIT=kin,FILE='dataj1/CH3OOH/CH3OOH_jpl85.abs',
            !	$	 STATUS='old')
            !	    OPEN(UNIT=kin,FILE='dataj1/CH3OOH/CH3OOH_jpl92.abs',
            !	$	 STATUS='old')
            !	    OPEN(UNIT=kin,FILE='dataj1/CH3OOH/CH3OOH_jpl94.abs',
            !	$	 STATUS='old')
            OPEN(UNIT=kin,FILE=trim(files(83)%fileName), STATUS='old')
            READ(kin,*) idum, n
            DO i = 1, idum-2
              READ(kin,*)
            END DO
            DO i = 1, n
              READ(kin,*) x1(i), y1(i)
              y1(i) = y1(i) * 1.e-20
            END DO
            CLOSE (kin)

            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,  	   0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,         1.e+38,0.)
            CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading 16'
              STOP
            END IF
         CASE (2)
            OPEN(UNIT=kin,FILE=trim(files(84)%fileName), STATUS='old')
            DO i = 1, 4
              READ(kin,*)
            END DO
            n = 32
            DO i = 1, n
              READ(kin,*) x1(i), y1(i)
              y1(i) = y1(i) * 1.e-20
            END DO
            CLOSE (kin)

            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,  	   0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,         1.e+38,0.)
            CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading 16'
              STOP
            END IF
         CASE (3)
            OPEN(UNIT=kin,FILE=trim(files(85)%fileName), STATUS='old')
            DO i = 1, 4
              READ(kin,*)
            END DO
            n = 12
            DO i = 1, n
              READ(kin,*) x1(i), y1(i)
            END DO
            CLOSE (kin)

            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,  	   0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,         1.e+38,0.)
            CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading 16'
              STOP
            END IF
         CASE (4)
            OPEN(UNIT=kin,FILE=trim(files(86)%fileName), STATUS='old')
            DO i = 1, 4
              READ(kin,*)
            END DO
            n = 15
            DO i = 1, n
              READ(kin,*) x1(i), y1(i)
            END DO
            CLOSE (kin)

            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,  	   0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,         1.e+38,0.)
            CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading 16'
              STOP
            END IF
      END SELECT

      !================================================================
      !R17
      !================================================================
      nReact=17
      kData = 2000
      CALL ChangeKData(kData)

      SELECT CASE (mOption(18))
         CASE (1)
     	   OPEN(UNIT=kin,FILE=trim(files(87)%fileName),STATUS='old')
     	   DO i = 1, 3
     	     READ(kin,*)
     	   END DO
     	   n = 15
     	   DO i = 1, n
     	     READ(kin,*) x1(i), y1(i)
     	   END DO
     	   CLOSE (kin)

     	   n1 = n
     	   CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
     	   CALL addpnt(x1,y1,kdata,n1,  	   0.,0.)
     	   CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
     	   CALL addpnt(x1,y1,kdata,n1,         1.e+38,0.)
     	   CALL inter2(nw,wl,yg(:,nReact),n1,x1,y1,ierr)
     	   IF (ierr /= 0) THEN
     	     WRITE(*,*) ierr, 'Reading R17'
     	     STOP
           END IF
        CASE (2)
           !	  sigma(T,lambda) = sigma(298,lambda) * exp(B * (T-298))
     	   OPEN(UNIT=kin,FILE=trim(files(88)%fileName),STATUS='old')
     	   DO i = 1, 4
     	     READ(kin,*)
     	   END DO
     	   n = 55
     	   DO i = 1, n
     	     READ(kin,*) x1(i), y1(i), y2(i)
     	     x2(i) = x1(i)
     	     y1(i) = y1(i) * 1.e-20
     	   END DO
     	   CLOSE (kin)

     	   n1 = n
     	   n2 = n
     	   CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
     	   CALL addpnt(x1,y1,kdata,n1,  	   0.,0.)
     	   CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
     	   CALL addpnt(x1,y1,kdata,n1,         1.e+38,0.)
     	   CALL inter2(nw,wl,yg(:,nReact),n1,x1,y1,ierr)
     	   IF (ierr /= 0) THEN
     	     WRITE(*,*) ierr, 'Reading R17'
     	     STOP
     	   END IF

     	   CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
     	   CALL addpnt(x2,y2,kdata,n2,  	   0.,y2(1))
     	   CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),y2(n2))
     	   CALL addpnt(x2,y2,kdata,n2,  	1.e+38,y2(n2))
     	   CALL inter2(nw,wl,yg1(:,nReact),n2,x2,y2,ierr)
     	   IF (ierr /= 0) THEN
     	     WRITE(*,*) ierr, 'Reading R17'
     	     STOP
     	   END IF
     	CASE (3)
     	   OPEN(UNIT=kin,FILE=trim(files(89)%fileName), STATUS='old')
     	   DO i = 1, 4
     	     READ(kin,*)
     	   END DO
     	   n = 13
     	   DO i = 1, n
     	     READ(kin,*) x1(i), y1(i)
     	     y1(i) = y1(i)*1E-20
     	   END DO
     	   CLOSE (kin)

     	   n1 = n
     	   CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
     	   CALL addpnt(x1,y1,kdata,n1,  	   0.,0.)
     	   CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
     	   CALL addpnt(x1,y1,kdata,n1,         1.e+38,0.)
     	   CALL inter2(nw,wl,yg(:,nReact),n1,x1,y1,ierr)
     	   IF (ierr /= 0) THEN
     	     WRITE(*,*) ierr, 'Reading R17'
     	     STOP
     	   END IF
        CASE (4)
           !	  sigma(T,lambda) = sigma(298,lambda) * 10**(B * T)
     	   OPEN(UNIT=kin,FILE=trim(files(90)%fileName), STATUS='old')
     	   DO i = 1, 4
     	     READ(kin,*)
     	   END DO
     	   n = 7
     	   DO i = 1, n
     	     READ(kin,*) x1(i), y1(i), y2(i)
     	     x2(i) = x1(i)
     	     y1(i) = y1(i) * 1.e-21
     	     y2(i) = y2(i) * 1.e-3
     	   END DO
     	   CLOSE (kin)

     	   n1 = n
     	   n2 = n
     	   CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),-36.)
     	   CALL addpnt(x1,y1,kdata,n1,  	   0.,-36.)
     	   CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),-36.)
     	   CALL addpnt(x1,y1,kdata,n1,         1.e+38,-36.)
     	   CALL inter2(nw,wl,yg(:,nReact),n1,x1,y1,ierr)
     	   IF (ierr /= 0) THEN
     	     WRITE(*,*) ierr, 'Reading R17'
     	     STOP
     	   END IF

     	   CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
     	   CALL addpnt(x2,y2,kdata,n2,  	   0.,y2(1))
     	   CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),y2(n2))
     	   CALL addpnt(x2,y2,kdata,n2,  	1.e+38,y2(n2))
     	   CALL inter2(nw,wl,yg1(:,nReact),n2,x2,y2,ierr)
     	   IF (ierr /= 0) THEN
     	     WRITE(*,*) ierr, 'Reading R17'
     	     STOP
     	   END IF
     	CASE (5)
     	   OPEN(UNIT=kin,FILE=trim(files(91)%fileName), STATUS='old')
     	   DO i = 1, 4
     	     READ(kin,*)
     	   END DO
     	   n = 13
     	   DO i = 1, n
     	     READ(kin,*) x1(i), y1(i)
     	   END DO
     	   CLOSE (kin)

     	   n1 = n
     	   CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
     	   CALL addpnt(x1,y1,kdata,n1,  	   0.,0.)
     	   CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
     	   CALL addpnt(x1,y1,kdata,n1,         1.e+38,0.)
     	   CALL inter2(nw,wl,yg(:,nReact),n1,x1,y1,ierr)
     	   IF (ierr /= 0) THEN
     	     WRITE(*,*) ierr, 'Reading R17'
     	     STOP
     	   END IF
        CASE (6)
     	   DO iw = 1, nw-1
     	     IF(wc(iw) > 284.) THEN
     	       yg(iw,nReact) = EXP(-1.044E-3*wc(iw)*wc(iw) + 0.5309*wc(iw) - 112.4)
     	     ELSE
     	       yg(iw,nReact) = 0.
     	     END IF
     	   END DO
     	CASE (7)
     	   OPEN(UNIT=kin,FILE=trim(files(92)%fileName), STATUS='old')
     	   DO i = 1, 4
     	     READ(kin,*)
     	   END DO
     	   n = 24
     	   DO i = 1, n
     	     READ(kin,*) x1(i), y1(i)
     	   END DO
     	   CLOSE (kin)

     	   n1 = n
     	   CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
     	   CALL addpnt(x1,y1,kdata,n1,  	   0.,0.)
     	   CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
     	   CALL addpnt(x1,y1,kdata,n1,         1.e+38,0.)
     	   CALL inter2(nw,wl,yg(:,nReact),n1,x1,y1,ierr)
     	   IF (ierr /= 0) THEN
     	     WRITE(*,*) ierr, 'Reading R17'
     	     STOP
     	   END IF
     	CASE (8)
     	   OPEN(UNIT=kin,FILE=trim(files(93)%fileName), STATUS='old')
     	   DO i = 1, 4
     	     READ(kin,*)
     	   END DO
     	   n = 1638
     	   DO i = 1, n
     	     READ(kin,*) x1(i), y1(i)
     	   END DO
     	   CLOSE (kin)

     	   n1 = n
     	   CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
     	   CALL addpnt(x1,y1,kdata,n1,  	   0.,0.)
     	   CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
     	   CALL addpnt(x1,y1,kdata,n1,         1.e+38,0.)
     	   CALL inter2(nw,wl,yg(:,nReact),n1,x1,y1,ierr)
     	   IF (ierr /= 0) THEN
     	     WRITE(*,*) ierr, 'Reading R17'
     	     STOP
     	   END IF
      END SELECT

      !================================================================
      !R18
      !================================================================
      nReact=18
      kData = 100
      CALL ChangeKData(kData)

      ! cross section from Senum et al., 1984, J.Phys.Chem. 88/7, 1269-1270

      !     OPEN(UNIT=kin,FILE='dataj1/RONO2/PAN_senum.abs',STATUS='OLD')
      !     DO i = 1, 14
      !        READ(kin,*)
      !     ENDDO
      !     n = 21
      !     DO i = 1, n
      !        READ(kin,*) x1(i), y1(i)
      !        y1(i) = y1(i) * 1.E-20
      !     ENDDO
      !     CLOSE(kin)

      !      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      !      CALL addpnt(x1,y1,kdata,n, 	      0.,0.)
      !      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      !      CALL addpnt(x1,y1,kdata,n, 	  1.e+38,0.)
      !      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      !      IF (ierr .NE. 0) THEN
      ! 	WRITE(*,*) ierr, 'Reading 18'
      ! 	STOP
      !      ENDIF

      ! cross section from
      !      Talukdar et al., 1995, J.Geophys.Res. 100/D7, 14163-14174
      OPEN(UNIT=kin,FILE=trim(files(94)%fileName),STATUS='OLD')
      DO i = 1, 14
        READ(kin,*)
      END DO
      n = 78
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i)
        y1(i) = y1(i) * 1.e-20
        y2(i) = y2(i) * 1E-3
        x2(i) = x1(i)
      END DO
      n2 = n
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	     0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	 1.e+38,0.)
      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
      IF (ierr /= 0) THEN
        WRITE(*,*) ierr, 'Reading 18'
        STOP
      END IF

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,	 0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,    1.e+38,0.)
      CALL inter2(nw,wl,yg2(:,nReact),n2,x2,y2,ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) ierr, 'Reading 18'
         STOP
      END IF

      !================================================================
      !R19
      !================================================================
      nReact=19
      kData = 100
      CALL ChangeKData(kData)

      !** cross sections from JPL94 recommendation
      OPEN(kin,FILE=trim(files(95)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      END DO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      END DO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	      0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	    1E38,0.)

      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
      IF (ierr /= 0) THEN
        WRITE(*,*) ierr, 'Reading R19'
        STOP
      END IF

      !================================================================
      !R20
      !================================================================
      nReact=20
      kData = 100
      CALL ChangeKData(kData)

      !** cross sections from JPL97 recommendation (identical to 94 data)

      OPEN(kin,FILE=trim(files(96)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i)
  	y1(i) = y1(i) * 1E-20
      END DO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1E38,0.)

      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R20'
  	STOP
      END IF

      !================================================================
      !R21
      !================================================================
      nReact=21
      kData = 100
      CALL ChangeKData(kData)

      !** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE=trim(files(97)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i)
  	y1(i) = y1(i) * 1E-20
      END DO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1E38,0.)

      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading 21'
  	STOP
      END IF

      !================================================================
      !R22
      !================================================================
      nReact=22
      kData = 100
      CALL ChangeKData(kData)

      !*** cross sections from JPL97 recommendation (identical to 94 recommendation)
      OPEN(kin,FILE=trim(files(98)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i)
  	y1(i) = y1(i) * 1E-20
      END DO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1E38,0.)

      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading 22'
  	STOP
      END IF

      !================================================================
      !R23
      !================================================================
      nReact=23
      kData = 100
      CALL ChangeKData(kData)

     !** cross sections from JPL97 recommendation (identical to 94 recommendation)
      OPEN(kin,FILE=trim(files(99)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i), y2(i)
  	y1(i) = y1(i) * 1E-20
  	y2(i) = y2(i) * 1E-20
  	x2(i) = x1(i)
      END DO
      CLOSE(kin)

      n1 = n
      n2 = n

      !* sigma @ 295 K

      CALL addpnt(x1,y1,kdata,n1, x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,	  0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,	1E38,0.)

      CALL inter2(nw,wl,yg1(:,nReact),n1,x1,y1,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R23'
  	STOP
      END IF

      ! sigma @ 210 K

      CALL addpnt(x2,y2,kdata,n2, x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,	  0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,	1E38,0.)

      CALL inter2(nw,wl,yg2(:,nReact),n2,x2,y2,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R23'
  	STOP
      END IF

      !================================================================
      !R24
      !================================================================
      nReact=24
      kData = 100
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(100)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i), y2(i)
  	y1(i) = y1(i) * 1E-20
  	y2(i) = y2(i) * 1E-20
  	x2(i) = x1(i)
      END DO
      CLOSE(kin)

      n1 = n
      n2 = n

      !* sigma @ 295 K

      CALL addpnt(x1,y1,kdata,n1, x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,	  0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,	1E38,0.)

      CALL inter2(nw,wl,yg1(:,nReact),n1,x1,y1,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R24'
  	STOP
      END IF

      ! sigma @ 210 K

      CALL addpnt(x2,y2,kdata,n2, x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,	  0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,	1E38,0.)

      CALL inter2(nw,wl,yg2(:,nReact),n2,x2,y2,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R24'
  	STOP
      END IF

      !================================================================
      !R25
      !================================================================
      nReact=25
      kData = 100
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(101)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i)
  	y1(i) = y1(i) * 1E-20
      END DO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1E38,0.)

      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading 25'
  	STOP
      END IF

      !================================================================
      !R26
      !================================================================
      nReact=26
      kData = 100
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(102)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i)
  	y1(i) = y1(i) * 1E-20
      END DO
      CLOSE(kin)

      !* sigma @ 298 K

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1E38,0.)

      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R26'
  	STOP
      END IF

      !================================================================
      !R27
      !================================================================
      nReact=27
      kData = 100
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(103)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i)
  	y1(i) = y1(i) * 1E-20
      END DO
      CLOSE(kin)

      !* sigma @ 298 K

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1E38,0.)

      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R27'
  	STOP
      END IF

      !================================================================
      !R28
      !================================================================
      nReact=28
      kData = 100
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(104)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i)
  	y1(i) = y1(i) * 1E-20
      END DO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1E38,0.)

      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'REading R28'
  	STOP
      END IF

      !================================================================
      !R29
      !================================================================
      nReact=29
      kData = 100
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(105)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i), y2(i), y3(i)
  	y1(i) = y1(i) * 1E-20
  	y2(i) = y2(i) * 1E-20
  	y3(i) = y3(i) * 1E-20
  	x2(i) = x1(i)
  	x3(i) = x1(i)
      END DO
      CLOSE(kin)

      n1 = n
      n2 = n
      n3 = n

      !* sigma @ 295 K

      CALL addpnt(x1,y1,kdata,n1, x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,	  0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,	1E38,0.)

      CALL inter2(nw,wl,yg1(:,nReact),n1,x1,y1,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R29'
  	STOP
      END IF

      !* sigma @ 250 K

      CALL addpnt(x2,y2,kdata,n2, x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,	  0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,	1E38,0.)

      CALL inter2(nw,wl,yg2(:,nReact),n2,x2,y2,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R29'
  	STOP
      END IF

      !* sigma @ 210 K

      CALL addpnt(x3,y3,kdata,n3, x3(1)*(1.-deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,	  0.,0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,	1E38,0.)

      CALL inter2(nw,wl,yg3(:,nReact),n3,x3,y3,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R29'
  	STOP
      END IF

      !================================================================
      !R30
      !================================================================
      nReact=30
      kData = 100
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(106)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i), y2(i), y3(i)
  	y1(i) = y1(i) * 1E-20
  	y2(i) = y2(i) * 1E-20
  	y3(i) = y3(i) * 1E-20
  	x2(i) = x1(i)
  	x3(i) = x1(i)
      END DO
      CLOSE(kin)

      n1 = n
      n2 = n
      n3 = n

      !* sigma @ 296 K

      CALL addpnt(x1,y1,kdata,n1, x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,	  0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,	1E38,0.)

      CALL inter2(nw,wl,yg1(:,nReact),n1,x1,y1,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R30'
  	STOP
      END IF

      !* sigma @ 279 K

      CALL addpnt(x2,y2,kdata,n2, x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,	  0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,	1E38,0.)

      CALL inter2(nw,wl,yg2(:,nReact),n2,x2,y2,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R30'
  	STOP
      END IF

      !* sigma @ 255 K

      CALL addpnt(x3,y3,kdata,n3, x3(1)*(1.-deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,	  0.,0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,	1E38,0.)

      CALL inter2(nw,wl,yg3(:,nReact),n3,x3,y3,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R30'
  	STOP
      END IF

      !================================================================
      !R31
      !================================================================
      nReact=31
      kData = 100
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(107)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i)
  	y1(i) = y1(i) * 1E-20
      END DO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1E38,0.)

      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R31'
  	STOP
      END IF

      !================================================================
      !R32
      !================================================================
      nReact=32
      kData = 100
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(108)%fileName),STATUS='OLD')
      READ(kin,*) idum
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      READ(kin,100) inline
      READ(inline(6:),*) tbar(nReact),i,(coeff(i,k,nReact),k=1,3)
      READ(kin,*)	  i,(coeff(i,k,nReact),k=1,3)
      READ(kin,*)	  i,(coeff(i,k,nReact),k=1,3)
      READ(kin,*)	  i,(coeff(i,k,nReact),k=1,3)
      CLOSE(kin)

      !================================================================
      !R33
      !================================================================
      nReact=33
      kData = 100
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(109)%fileName),STATUS='OLD')
      READ(kin,*) idum
      idum = idum+5
      DO i = 1, idum-2
        READ(kin,*)
      END DO
      READ(kin,100) inline
      READ(inline(6:),*) tbar(nReact),i,(coeff(i,k,nReact),k=1,3)
      READ(kin,*)           i,(coeff(i,k,nReact),k=1,3)
      READ(kin,*)           i,(coeff(i,k,nReact),k=1,3)
      READ(kin,*)           i,(coeff(i,k,nReact),k=1,3)
      CLOSE(kin)

      !================================================================
      !R34
      !================================================================
      nReact=34
      kData = 100
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(110)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i)
  	y1(i) = y1(i) * 1E-20
      END DO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1E38,0.)

      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'REading R34'
  	STOP
      END IF

      !================================================================
      !R35
      !================================================================
      nReact=35
      kData = 100
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(111)%fileName),STATUS='OLD')
      READ(kin,*) idum
      idum = idum+10
      DO i = 1, idum-2
        READ(kin,*)
      END DO
      READ(kin,101) inline
      READ(inline(6:),*) tbar(nReact),i,(coeff(i,k,nReact),k=1,3)
      READ(kin,*)           i,(coeff(i,k,nReact),k=1,3)
      READ(kin,*)           i,(coeff(i,k,nReact),k=1,3)
      READ(kin,*)           i,(coeff(i,k,nReact),k=1,3)
      CLOSE(kin)

      !================================================================
      !R36
      !================================================================
      nReact=36
      kData = 100
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(112)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
 	READ(kin,*)
      END DO
      DO i = 1, n
 	READ(kin,*) x1(i), y1(i)
 	y1(i) = y1(i) * 1E-20
      END DO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1E38,0.)

      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)

      IF (ierr /= 0) THEN
 	WRITE(*,*) ierr, 'Reading R36'
 	STOP
      END IF

      !================================================================
      !R37
      !================================================================
      nReact=37
      kData = 100
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(113)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i)
  	y1(i) = y1(i) * 1E-20
      END DO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1E38,0.)

      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R37'
  	STOP
      END IF

      !================================================================
      !R38
      !================================================================
      nReact=38
      kData = 100
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(114)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i), y2(i), y3(i), y4(i), y5(i)
  	y1(i) = y1(i) * 1E-20
  	y2(i) = y2(i) * 1E-20
  	y3(i) = y3(i) * 1E-20
  	y4(i) = y4(i) * 1E-20
  	y5(i) = y5(i) * 1E-20
  	x2(i) = x1(i)
  	x3(i) = x1(i)
  	x4(i) = x1(i)
  	x5(i) = x1(i)
      END DO
      CLOSE(kin)

      n1 = n
      n2 = n
      n3 = n
      n4 = n
      n5 = n

      !* sigma @ 295 K

      CALL addpnt(x1,y1,kdata,n1, x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,	  0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,	1E38,0.)

      CALL inter2(nw,wl,yg1(:,nReact),n1,x1,y1,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R38'
  	STOP
      END IF

      !* sigma @ 270 K

      CALL addpnt(x2,y2,kdata,n2, x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,	  0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,	1E38,0.)

      CALL inter2(nw,wl,yg2(:,nReact),n2,x2,y2,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R38'
  	STOP
      END IF

      !* sigma @ 250 K

      CALL addpnt(x3,y3,kdata,n3, x3(1)*(1.-deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,	  0.,0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,	1E38,0.)

      CALL inter2(nw,wl,yg3(:,nReact),n3,x3,y3,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R38'
  	STOP
      END IF

      !* sigma @ 230 K

      CALL addpnt(x4,y4,kdata,n4, x4(1)*(1.-deltax),0.)
      CALL addpnt(x4,y4,kdata,n4,	  0.,0.)
      CALL addpnt(x4,y4,kdata,n4,x4(n4)*(1.+deltax),0.)
      CALL addpnt(x4,y4,kdata,n4,	1E38,0.)

      CALL inter2(nw,wl,yg4(:,nReact),n4,x4,y4,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R38'
  	STOP
      END IF

      !* sigma @ 210 K

      CALL addpnt(x5,y5,kdata,n5, x5(1)*(1.-deltax),0.)
      CALL addpnt(x5,y5,kdata,n5,	  0.,0.)
      CALL addpnt(x5,y5,kdata,n5,x5(n5)*(1.+deltax),0.)
      CALL addpnt(x5,y5,kdata,n5,	1E38,0.)

      CALL inter2(nw,wl,yg5(:,nReact),n5,x5,y5,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R38'
  	STOP
      END IF

      !================================================================
      !R39
      !================================================================
      nReact=39
      kData = 100
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(115)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i)
  	y1(i) = y1(i) * 1E-20
      END DO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1E38,0.)

      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R39'
  	STOP
      END IF

      !================================================================
      !R40
      !================================================================
      nReact=40
      kData = 100
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(116)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i)
  	y1(i) = y1(i) * 1E-20
      END DO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1E38,0.)

      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R40'
  	STOP
      END IF

      !================================================================
      !R41
      !================================================================
      nReact=41
      kData = 100
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(117)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i)
  	y1(i) = y1(i) * 1E-20
      END DO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1E38,0.)

      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading 41'
  	STOP
      END IF

      !================================================================
      !R42
      !================================================================
      nReact=42
      kData = 100
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(118)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i)
  	y1(i) = y1(i) * 1E-20
      END DO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1E38,0.)

      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R42'
  	STOP
      END IF

      !================================================================
      !R43
      !================================================================
      nReact=43
      kData = 100
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(119)%fileName),STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
  	READ(kin,*)
      END DO
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i)
  	y1(i) = y1(i) * 1E-20
      END DO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)+deltax,0.)
      CALL addpnt(x1,y1,kdata,n,      1E38,0.)

      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)

      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading 43'
  	STOP
      END IF

      !================================================================
      !R44
      !================================================================
      ! Fixed, not read.

      !================================================================
      !R45
      !================================================================
      nReact=45
      kData = 150
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(120)%fileName),STATUS='OLD')
      n = 119
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i), y2(i), y3(i)
  	y1(i) = y1(i) * 1E-20
  	x2(i) = x1(i)
  	x3(i) = x1(i)
      END DO
      CLOSE(kin)

      n1 = n
      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,	 0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,      1E38,0.)
      CALL inter2(nw,wl,yg1(:,nReact),n1,x1,y1,ierr)
      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading 45'
  	STOP
      END IF

      n2 = n
      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,	 0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,      1E38,0.)
      CALL inter2(nw,wl,yg2(:,nReact),n2,x2,y2,ierr)
      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading 45'
  	STOP
      END IF

      n3 = n
      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,	 0.,0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,      1E38,0.)
      CALL inter2(nw,wl,yg3(:,nReact),n3,x3,y3,ierr)
      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading 45'
  	STOP
      END IF

      !================================================================
      !R46
      !================================================================
      nReact=46
      kData = 100
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(121)%fileName),STATUS='OLD')
      DO i = 1, 13
  	READ(kin,*)
      END DO
      n = 61
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i)
  	y1(i) = y1(i) * 1E-20
      END DO
      CLOSE(kin)

      n1 = n
      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,	   0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,	 1E38,0.)
      CALL inter2(nw,wl,yg1(:,nReact),n1,x1,y1,ierr)
      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading 46'
  	STOP
      END IF

      !================================================================
      !R47
      !================================================================
      nReact=47
      kData = 150
      CALL ChangeKData(kData)

      OPEN(kin,FILE=trim(files(122)%fileName),STATUS='OLD')
      DO i = 1, 5
  	READ(kin,*)
      END DO
      n = 22
      DO i = 1, n
  	READ(kin,*) x1(i), y1(i)
  	y1(i) = y1(i) * 1E-20
      END DO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,	0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1E38,0.)
      CALL inter2(nw,wl,yg(:,nReact),n,x1,y1,ierr)
      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R47'
  	STOP
      END IF

      !================================================================
      !R48 until R100
      !================================================================
      !Does not exist!

      !================================================================
      !R101
      !================================================================
      nReact=101
      kData = 300
      CALL ChangeKData(kData)

      OPEN(UNIT=kin,FILE=trim(files(123)%fileName), STATUS='old')
      DO i = 1, 15
  	READ(kin,*)
      END DO
      n = 131
      DO i = 1, n
  	READ(kin,*) x(i), y(i)
      END DO
      CLOSE(kin)

      CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
      CALL addpnt(x,y,kdata,n,  	   0.,0.)
      CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
      CALL addpnt(x,y,kdata,n,         1.e+38,0.)
      CALL inter2(nw,wl,yg(:,nReact),n,x,y,ierr)
      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R101'
  	STOP
      END IF

      !================================================================
      !R102
      !================================================================
      nReact=102
      kData = 300
      CALL ChangeKData(kData)

      SELECT CASE (mOption(19))
         CASE (1)
            OPEN(UNIT=kin,FILE=trim(files(124)%fileName), STATUS='old')
            DO i = 1, 7
              READ(kin,*)
            END DO
            n = 55
            DO i = 1, n
              READ(kin,*) x(i), y(i)
              y(i) = y(i) * 1.e-20
            END DO
            CLOSE(kin)

            CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
            CALL addpnt(x,y,kdata,n,		   0.,0.)
            CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
            CALL addpnt(x,y,kdata,n,	       1.e+38,0.)
            CALL inter2(nw,wl,yg(:,nReact),n,x,y,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading R102'
              STOP
            END IF
         CASE (2)
            OPEN(UNIT=kin,FILE=trim(files(125)%fileName), STATUS='old')
            DO i = 1, 8
              READ(kin,*)
            END DO
            n = 287
            DO i = 1, n
              READ(kin,*) x(i), y(i)
              y(i) = y(i) * 1.e-20
            END DO
            CLOSE(kin)

            CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
            CALL addpnt(x,y,kdata,n,		 0.,0.)
            CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
            CALL addpnt(x,y,kdata,n,	     1.e+38,0.)
            CALL inter2(nw,wl,yg(:,nReact),n,x,y,ierr)
            IF (ierr /= 0) THEN
              WRITE(*,*) ierr, 'Reading R102'
              STOP
            END IF
      END SELECT

      !================================================================
      !R103
      !================================================================
      nReact=103
      kData = 20000
      CALL ChangeKData(kData)

      OPEN(UNIT=kin,FILE=trim(files(126)%fileName), STATUS='old')
      DO i = 1, 9
  	READ(kin,*)
      END DO
      n = 19682
      DO i = 1, n
  	READ(kin,*) x(i), y(i)
      END DO
      CLOSE(kin)

      CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
      CALL addpnt(x,y,kdata,n,  	   0.,0.)
      CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
      CALL addpnt(x,y,kdata,n,         1.e+38,0.)
      CALL inter2(nw,wl,yg(:,nReact),n,x,y,ierr)
      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R103'
  	STOP
      END IF

      !================================================================
      !R104
      !================================================================
      nReact=104
      kData = 20000
      CALL ChangeKData(kData)

      OPEN(UNIT=kin,FILE=trim(files(127)%fileName), STATUS='old')
      DO i = 1, 10
  	READ(kin,*)
      END DO
      n = 15213
      DO i = 1, n
  	READ(kin,*) x(i), y(i)
      END DO
      CLOSE(kin)

      CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
      CALL addpnt(x,y,kdata,n,  	   0.,0.)
      CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
      CALL addpnt(x,y,kdata,n,         1.e+38,0.)
      CALL inter2(nw,wl,yg(:,nReact),n,x,y,ierr)
      !PRINT *,'React104 yg=',yg(:,nReact); CALL flush(6)
      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R104'
  	STOP
      END IF

      !================================================================
      !R105
      !================================================================
      nReact=105
      kData = 20000
      CALL ChangeKData(kData)

      OPEN(UNIT=kin,FILE=trim(files(128)%fileName), STATUS='old')
      DO i = 1, 8
  	READ(kin,*)
      END DO
      n = 148
      DO i = 1, n
  	READ(kin,*) x(i), y(i)
  	y(i) = y(i) * 1.e-20
      END DO
      CLOSE(kin)

      CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
      CALL addpnt(x,y,kdata,n,  	   0.,0.)
      CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
      CALL addpnt(x,y,kdata,n,         1.e+38,0.)
      CALL inter2(nw,wl,yg(:,nReact),n,x,y,ierr)
      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R105'
  	STOP
      END IF

      !================================================================
      !R106
      !================================================================
      nReact=106
      kData = 200
      CALL ChangeKData(kData)

      OPEN(UNIT=kin,FILE=trim(files(129)%fileName), STATUS='old')
      DO i = 1, 10
  	READ(kin,*)
      END DO
      n1 = 0
      n2 = 0
      DO i = 1, 63
  	READ(kin,*) x1(i), dum, dum, y1(i), y2(i), dum, dum
  	IF (y1(i) > 0.) n1 = n1 + 1
  	IF (y2(i) > 0.) n2 = n2 + 1
  	x2(i) = x1(i)
  	y1(i) = y1(i) * 1.e-20
  	y2(i) = y2(i) * 1.e-3
      END DO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,	      0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,	  1.e+38,0.)
      CALL inter2(nw,wl,yg1(:,nReact),n1,x1,y1,ierr)
      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R106'
  	STOP
      END IF

      CALL addpnt(x2,y2,kdata,n2,	      0.,y2(1))
      CALL addpnt(x2,y2,kdata,n2,	  1.e+38,y2(n2))
      CALL inter2(nw,wl,yg2(:,nReact),n2,x2,y2,ierr)
      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R106'
  	STOP
      END IF

      !================================================================
      !R107
      !================================================================
      nReact=107
      kData = 200
      CALL ChangeKData(kData)

      OPEN(UNIT=kin,FILE=trim(files(130)%fileName), STATUS='old')
      DO i = 1, 10
  	READ(kin,*)
      END DO
      n1 = 0
      n2 = 0
      DO i = 1, 63
  	READ(kin,*) x1(i), dum, dum, dum, dum, y1(i), y2(i)
  	IF (y1(i) > 0.) n1 = n1 + 1
  	IF (y2(i) > 0.) n2 = n2 + 1
  	x2(i) = x1(i)
  	y1(i) = y1(i) * 1.e-20
  	y2(i) = y2(i) * 1.e-3
      END DO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,	      0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,	  1.e+38,0.)
      CALL inter2(nw,wl,yg1(:,nReact),n1,x1,y1,ierr)
      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R107'
  	STOP
      END IF

      CALL addpnt(x2,y2,kdata,n2,	      0.,y2(1))
      CALL addpnt(x2,y2,kdata,n2,	  1.e+38,y2(n2))
      CALL inter2(nw,wl,yg2(:,nReact),n2,x2,y2,ierr)
      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R107'
  	STOP
      END IF

      !================================================================
      !R108 until R110
      !================================================================
      !Do not read - fixed

      !================================================================
      !R111
      !================================================================
      nReact=111
      kData = 20000
      CALL ChangeKData(kData)

      OPEN(UNIT=kin,FILE=trim(files(131)%fileName), STATUS='old')
      DO i = 1, 25
  	READ(kin,*)
      END DO
      n = 131
      DO i = 1, n
  	READ(kin,*) x(i), y(i)
  	y(i) = y(i) * 1.e-20
      END DO
      CLOSE(kin)

      CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
      CALL addpnt(x,y,kdata,n,  	   0.,0.)
      CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
      CALL addpnt(x,y,kdata,n,         1.e+38,0.)
      CALL inter2(nw,wl,yg(:,nReact),n,x,y,ierr)
      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R111'
  	STOP
      END IF

      !================================================================
      !R112
      !================================================================
      nReact=112
      kData = 20000
      CALL ChangeKData(kData)

      OPEN(UNIT=kin,FILE=trim(files(132)%fileName), STATUS='old')
      DO i = 1, 8
  	READ(kin,*)
      END DO
      n = 101
      DO i = 1, n
  	READ(kin,*) x(i), y(i)
      END DO
      CLOSE(kin)

      CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
      CALL addpnt(x,y,kdata,n,  	   0.,0.)
      CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
      CALL addpnt(x,y,kdata,n,         1.e+38,0.)
      CALL inter2(nw,wl,yg(:,nReact),n,x,y,ierr)
      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading 112'
  	STOP
      END IF

      !================================================================
      !R113
      !================================================================
      !Do not read - fixed

      !================================================================
      !R114
      !================================================================
      nReact=114
      kData = 20000
      CALL ChangeKData(kData)

      OPEN(UNIT=kin,FILE=trim(files(133)%fileName), STATUS='old')
      DO i = 1, 14
  	READ(kin,*)
      END DO
      n = 15
      DO i = 1, n
  	READ(kin,*) x(i), dum, y(i)
  	y(i) = y(i) * 1.e-20
      END DO
      n = n + 1
      x(n) = dum
      CLOSE(kin)

      ! use bin-to-bin interpolation

      CALL inter4(nw,wl,yg(:,nReact),n,x,y,1)

      !================================================================
      !R115
      !================================================================
      nReact=115
      kData = 50
      CALL ChangeKData(kData)

      OPEN(UNIT=kin,FILE=trim(files(134)%fileName), STATUS='old')

      DO i = 1, 6
  	READ(kin,*)
      END DO
      n = 29
      DO i = 1, n
  	READ(kin,*) x(i),  y(i)
      END DO
      CLOSE(kin)

      CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
      CALL addpnt(x,y,kdata,n,  	   0.,0.)
      CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
      CALL addpnt(x,y,kdata,n,         1.e+38,0.)
      CALL inter2(nw,wl,yg(:,nReact),n,x,y,ierr)
      IF (ierr /= 0) THEN
  	WRITE(*,*) ierr, 'Reading R115'
  	STOP
      END IF

100   FORMAT(a120)
101   FORMAT(a80)


  END SUBROUTINE ReadAll


    INTEGER FUNCTION reverse(pos,pmax)
       INTEGER,INTENT(IN) :: pos,pmax
       reverse=pmax-pos+1
    END FUNCTION reverse

    SUBROUTINE alert(message)
      CHARACTER(LEN=*),INTENT(IN) :: message

      WRITE(noPr,*) message
      CALL flush(noPr)

    END SUBROUTINE alert

    INTEGER FUNCTION get_zlimit()
       get_zlimit=zLimit
    END FUNCTION


END MODULE ModTuv
