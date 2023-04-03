module domain_decomp
  use ModNamelistFile, only: namelistFile

  implicit none

  include "files.h"
  
  character(len=f_name_length) :: domain_fname ! from RAMSIN

contains
  subroutine StoreNamelistFileAtDomain_decomp(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    domain_fname = oneNamelistFile%domain_fname
  end subroutine StoreNamelistFileAtDomain_decomp

end module domain_decomp
