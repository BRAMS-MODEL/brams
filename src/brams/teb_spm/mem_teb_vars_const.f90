!############################# Change Log ##################################
! 5.0.2
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Brazilian Regional Atmospheric Modeling System - BRAMS
!###########################################################################

MODULE teb_vars_const

  use ModNamelistFile, only: namelistFile

  USE grid_dims, ONLY: &
       maxsteb,        & ! INTENT(IN)
       maxubtp           ! INTENT(IN)

  IMPLICIT NONE

  !**From Edmilson Freitas, IAG/USP
  INTEGER         :: iteb !Flag for urban parameterization using TEB (1=on, 0=off), from RAMSIN
  INTEGER         :: nteb !number of roof, road and wall layers used in TEB  - EDF, from RAMSIN
  REAL            :: d_road(maxsteb) ! from RAMSIN
  REAL            :: tc_road(maxsteb) ! from RAMSIN
  REAL            :: hc_road(maxsteb) ! from RAMSIN
  REAL            :: d_wall(maxsteb) ! from RAMSIN
  REAL            :: tc_wall(maxsteb) ! from RAMSIN
  REAL            :: hc_wall(maxsteb) ! from RAMSIN
  REAL            :: d_roof(maxsteb) ! from RAMSIN
  REAL            :: tc_roof(maxsteb) ! from RAMSIN
  REAL            :: hc_roof(maxsteb) ! from RAMSIN
  REAL            :: tminbld   !Minimum internal building temperature, from RAMSIN
  REAL            :: rushh1    !Morning Rush Hour, from RAMSIN
  REAL            :: rushh2    !Afternoon/Evening Rush Hour, from RAMSIN
  REAL            :: daylight  !Daylight saving time (horario de verao), from RAMSIN
  INTEGER         :: nurbtype    !Number of urban types, from RAMSIN
  INTEGER         :: ileafcod(maxubtp)   !Leaf class code to identify each urban type, from RAMSIN
  REAL            :: z0_town(maxubtp) !  from RAMSIN
  REAL            :: bld(maxubtp) !  from RAMSIN
  REAL            :: bld_height(maxubtp) !  from RAMSIN
  REAL            :: bld_hl_ratio(maxubtp) !  from RAMSIN
  REAL            :: aroof(maxubtp) !  from RAMSIN
  REAL            :: eroof(maxubtp) !  from RAMSIN
  REAL            :: aroad(maxubtp) !  from RAMSIN
  REAL            :: eroad(maxubtp) !  from RAMSIN
  REAL            :: awall(maxubtp) !  from RAMSIN
  REAL            :: ewall(maxubtp) !  from RAMSIN
  REAL            :: htraf(maxubtp) !  from RAMSIN
  REAL            :: hindu(maxubtp) !  from RAMSIN
  REAL            :: pletraf(maxubtp) !  from RAMSIN
  REAL            :: pleindu(maxubtp) !  from RAMSIN
  REAL, PARAMETER :: xpi= 3.1415
  REAL, PARAMETER :: xkarman= 0.4
  REAL, PARAMETER :: xstefan = 5.6697E-08
  REAL, PARAMETER :: xday= 86400.
  REAL, PARAMETER :: xg=9.80665
  REAL, PARAMETER :: xrd=287.0597
  REAL, PARAMETER :: xrv=461.5250
  REAL, PARAMETER :: xcpd=1004.70895
  REAL, PARAMETER :: xlvtt=2.5008E+6
  REAL, PARAMETER :: xtt=273.16


contains

  subroutine StoreNamelistFileAtTeb_vars_const(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    rushh1 = oneNamelistFile%rushh1
    rushh2 = oneNamelistFile%rushh2
    daylight = oneNamelistFile%daylight
    iteb = oneNamelistFile%iteb
    tminbld = oneNamelistFile%tminbld
    nteb = oneNamelistFile%nteb
    hc_roof = oneNamelistFile%hc_roof
    tc_roof = oneNamelistFile%tc_roof
    d_roof = oneNamelistFile%d_roof
    hc_road = oneNamelistFile%hc_road
    d_road = oneNamelistFile%d_road
    tc_road = oneNamelistFile%tc_road
    d_wall = oneNamelistFile%d_wall
    tc_wall = oneNamelistFile%tc_wall
    hc_wall = oneNamelistFile%hc_wall
    nurbtype = oneNamelistFile%nurbtype
    ileafcod = oneNamelistFile%ileafcod
    z0_town = oneNamelistFile%z0_town
    bld = oneNamelistFile%bld
    bld_height = oneNamelistFile%bld_height
    bld_hl_ratio = oneNamelistFile%bld_hl_ratio
    aroof = oneNamelistFile%aroof
    eroof = oneNamelistFile%eroof
    aroad = oneNamelistFile%aroad
    eroad = oneNamelistFile%eroad
    awall = oneNamelistFile%awall
    ewall = oneNamelistFile%ewall
    htraf = oneNamelistFile%htraf
    hindu = oneNamelistFile%hindu
    pletraf = oneNamelistFile%pletraf
    pleindu = oneNamelistFile%pleindu
  end subroutine StoreNamelistFileAtTeb_vars_const
END MODULE teb_vars_const
