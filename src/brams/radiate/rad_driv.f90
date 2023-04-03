
MODULE RADIATION

CONTAINS

  SUBROUTINE radiate(mzp, mxp, myp, ia, iz, ja, jz, mynum)

    USE mem_radiate, only: ilwrtyp, iswrtyp
    USE rrtm_driv  , only: rrtm_driver
    USE carma_driv , only: carma_driver

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: mzp, mxp, myp, ia, iz, ja, jz, mynum

    if &
       ((ilwrtyp + iswrtyp)==0) return ! teste

    if &
      ( ilwrtyp==6 .and. iswrtyp==6) &
      then

       call rrtm_driver(mzp, mxp, myp, ia, iz, ja, jz, mynum)

    elseif( ilwrtyp==4 .and. iswrtyp==4) then

       call carma_driver(mzp,mxp,myp,ia,iz,ja,jz,mynum) !teste 2

    else
       stop "unknown radiation scheme"
    endif

  END SUBROUTINE radiate

END MODULE RADIATION
