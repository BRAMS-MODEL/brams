! Flag to Active the TEB, Emission and/or SPM schemes
module teb_spm_start

  use ModNamelistFile, only: namelistFile

  implicit none

  integer :: teb_spm ! from RAMSIN

contains
  subroutine StoreNamelistFileAtTeb_spm_start(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    teb_spm = oneNamelistFile%teb_spm
  end subroutine StoreNamelistFileAtTeb_spm_start

end module teb_spm_start

