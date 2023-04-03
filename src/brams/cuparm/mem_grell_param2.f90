! Module necessary to Grell Cumulus param.
! Variables used to define scratch array dimensions

module mem_grell_param
  use ModNamelistFile, only: namelistFile

  !Memory for Grell's Cumulus scheme
  integer :: ngrids_cp, mgmxp, mgmyp, mgmzp
  integer, parameter ::               &
       maxiens    = 3,  & !Cloud spectral size
       maxens     = 3,  & ! 3  ensemble one on cap_max
       maxens2    = 3,  & ! 3  ensemble two on precip efficiency
       maxens_sh  = 3,  & ! 3  ensemble one on mbdt
       maxens2_sh = 1,  & ! 1  ensemble two on precip efficiency
       maxens3_sh = 10, & !10 ensemble three done in cup_forcing_ens16
       maxens3    = 16, & !16 ensemble three done in cup_forcing_ens16

       maxens_g3d  = 1, & ! 1  ensemble one on cap_max for G3d
       maxens2_g3d = 1, & ! 1  ensemble two on precip efficiency for G3d
       maxens3_g3d = 16   !16 ensemble three done in cup_forcing_ens16 for G3d

  integer :: ensdim,ensdim_g3d                     !Ensemble dimension

  integer :: icoic        ! Closure choice for deep
  integer :: icoic_sh     ! Closure choice for shallow

  integer :: Flag_Grell = 0  ! = 0 Grell Arrays not allocated

  !     mgm*p not determined
  ! = 1 Grell Arrays not allocated
  !     mgm*p DETERMINED
  ! = 2 Grell Arrays ALLOCATED
  !     mgm*p DETERMINED

  character (len=2) :: CLOSURE_TYPE  ! For new G.Grell Parameterization, From RAMSIN

contains

  subroutine define_memory(mmxp, mmyp, mmzp, ngrids, nnqparm,nnschu)

    implicit none
    integer, dimension (*) ::   &
        mmxp,                  &  ! Number of points in X direction
        mmyp,                  &  ! Number of points in Y direction
        mmzp,                  &  ! Number of points in Z direction
        nnqparm,	       &
	nnschu                 ! Flags for cumulus parameterization
    ! indexed by number of grids
    ! The above integers data are passed by arguments and can be the amount
    ! of points in a node or the total points in a grid
    integer, intent(in) :: ngrids  ! Number of grids (nested)

    ! Local Variables
    integer :: i

    !Ensemble dimension
!srf - out 2004
    !ensdim=maxiens*maxens*maxens2*maxens3
     ensdim     =1*maxens    *maxens2    *maxens3
     ensdim_g3d =1*maxens_g3d*maxens2_g3d*maxens3_g3d

    ! Allocate arrays based on options (if necessary)
    !Memory for Grell's Cumulus Scheme

    mgmxp     = 0
    mgmyp     = 0
    mgmzp     = 0
    ngrids_cp = 0

    do i=1, ngrids
!srf-g3d
       if (nnqparm(i) >= 2 .or. nnschu(i) >= 2)  then

          mgmxp = max(mgmxp,mmxp(i))
          mgmyp = max(mgmyp,mmyp(i))
          mgmzp = max(mgmzp,mmzp(i))
          ngrids_cp = ngrids_cp + 1
       endif
    enddo

    Flag_Grell = 1  ! Seting the Flag

  end subroutine define_memory

  subroutine StoreNamelistFileAtMem_grell_param(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    closure_type = oneNamelistFile%closure_type
  end subroutine StoreNamelistFileAtMem_grell_param
end module mem_grell_param
