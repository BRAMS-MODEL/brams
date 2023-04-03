!############################# Change Log ##################################
! 5.0.2
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Brazilian Regional Atmospheric Modeling System - BRAMS
!###########################################################################

module mem_emiss

  use ModNamelistFile, only: namelistFile

  include "files.h"

  !**From Edmilson Freitas, IAG/USP
  integer             :: isource     ! Flag for using emission module,  (1=on, 0=0ff) - EDF, from RAMSIN
  real                :: eindno  !Industrial emissions for NO, from RAMSIN
  real                :: eindno2 !Industrial emissions for NO2, from RAMSIN
  real                :: eindpm !Industrial emissions for PM25, from RAMSIN
  real                :: eindco !Industrial emissions for CO, from RAMSIN
  real                :: eindso2 !Industrial emissions for SO2, from RAMSIN
  real                :: eindvoc !Industrial emissions for VOC's, from RAMSIN
  real                :: eveino  !Vehicular emissions for NO, from RAMSIN
  real                :: eveino2 !Vehicular emissions for NO2, from RAMSIN
  real                :: eveipm  !Vehicular emissions for PM25, from RAMSIN
  real                :: eveico  !Vehicular emissions for CO, from RAMSIN
  real                :: eveiso2 !Vehicular emissions for SO2, from RAMSIN
  real                :: eveivoc !Vehicular emissions for VOC's, from RAMSIN
  character(len=3)    :: weekdayin   !Initial weekday of simulation to findout rates of anthropogenic and emission sources, from RAMSIN
  real                :: efsat !emission fractions for saturdays, from RAMSIN
  real                :: efsun !emission fractions for sundays, from RAMSIN
  integer             :: ichemi    !for photochemical module activation - EDF, from RAMSIN
  integer             :: ichemi_in   !flag for reading a previous run as initial values, from RAMSIN
  character (len=f_name_length) :: chemdata_in !path for initial values reading, from RAMSIN

contains
  subroutine StoreNamelistFileAtMem_emiss(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    ichemi = oneNamelistFile%ichemi
    ichemi_in = oneNamelistFile%ichemi_in
    chemdata_in = oneNamelistFile%chemdata_in
    isource = oneNamelistFile%isource
    weekdayin = oneNamelistFile%weekdayin
    efsat = oneNamelistFile%efsat
    efsun = oneNamelistFile%efsun
    eindno = oneNamelistFile%eindno
    eindno2 = oneNamelistFile%eindno2
    eindpm = oneNamelistFile%eindpm
    eindco = oneNamelistFile%eindco
    eindso2 = oneNamelistFile%eindso2
    eindvoc = oneNamelistFile%eindvoc
    eveino = oneNamelistFile%eveino
    eveino2 = oneNamelistFile%eveino2
    eveipm = oneNamelistFile%eveipm
    eveico = oneNamelistFile%eveico
    eveiso2 = oneNamelistFile%eveiso2
    eveivoc = oneNamelistFile%eveivoc
  end subroutine StoreNamelistFileAtMem_emiss
end module mem_emiss
