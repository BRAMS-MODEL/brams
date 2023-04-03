! Module necessary to Soil Moisture Init.

MODULE memSoilMoisture

  IMPLICIT NONE

  include "files.h"

  CHARACTER(len=1)   :: soil_moist       ! from RAMSIN
  CHARACTER(len=1)   :: soil_moist_fail  ! from RAMSIN
  CHARACTER(len=f_name_length) :: usdata_in        ! from RAMSIN
  CHARACTER(len=f_name_length) :: usmodel_in       ! from RAMSIN

END MODULE memSoilMoisture
