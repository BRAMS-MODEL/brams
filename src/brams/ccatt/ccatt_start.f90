! Active or not the CCATT module
module ccatt_start

  use ModNamelistFile, only: namelistFile

  implicit none

  integer :: ccatt ! from CHEM_RAMSIN

contains
  subroutine StoreNamelistFileAtCCatt_start(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    ccatt = oneNamelistFile%ccatt
  end subroutine StoreNamelistFileAtCCatt_start

end module ccatt_start
