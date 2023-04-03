MODULE chem1aq_list
  IMPLICIT NONE
  
  
  INTEGER,PARAMETER :: maxnspeciesaq = 01
  INTEGER,PARAMETER :: nspeciesaq    = 01
  CHARACTER(LEN=8) :: spcaq_name(maxnspeciesaq)
  INTEGER,PARAMETER,DIMENSION(nspeciesaq) :: ind_gas=(/&
    001    & ! O3aq - 001
    /)

END MODULE chem1aq_list
