!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_turb

  use ModNamelistFile, only: namelistFile

  use grid_dims, only:  maxgrds ! INTENT(IN)
  use mem_stilt, only:  imassflx! INTENT(IN)

  type turb_vars

     ! Variables to be dimensioned by (nzp,nxp,nyp)
     real, pointer, dimension(:,:,:) :: &
          tkep, epsp, hkm, vkm, vkh, cdrag

     ! Variables to be dimensioned by (nxp,nyp)
     real, pointer, dimension(:,:) :: &
          sflux_u, sflux_v, sflux_w, sflux_t, sflux_r

    ![MLO - For Nakanishi/Niino
    !srf     integer, pointer, dimension(:,:) :: kpbl
     real, pointer, dimension(:,:) :: kpbl
    !MLO]

  end type turb_vars

  type (turb_vars), allocatable :: turb_g(:), turbm_g(:)

  integer :: if_urban_canopy ! from RAMSIN
  integer :: ihorgrad        ! from RAMSIN
  integer :: idiffk(maxgrds) ! from RAMSIN
  real    :: zkhkm(maxgrds)  ! from RAMSIN
  real    :: xkhkm(maxgrds)  ! from RAMSIN
  real    :: csz(maxgrds)    ! from RAMSIN
  real    :: csx(maxgrds)    ! from RAMSIN
  real    :: akmin(maxgrds)  ! from RAMSIN

  real    :: brunt
  real    :: rmax
  real    :: rmin

contains

  subroutine alloc_turb(turb, n1, n2, n3, ng)

    implicit none
    ! Arguments:
    type (turb_vars), intent(INOUT) :: turb
    integer, intent(IN) :: n1, n2, n3, ng
    ! Local Varaibles:
    integer :: ierr
    character(len=*), parameter :: h="**(alloc_turb)**"

    ! Allocate arrays based on options (if necessary)
    ierr = 0
    if(idiffk(ng) == 1 .or. idiffk(ng) == 4 .or.  &
       idiffk(ng) == 5 .or. idiffk(ng) == 6 .or.  &
       idiffk(ng) == 7                           ) &
        allocate (turb%tkep(n1,n2,n3), STAT=ierr)
         if (ierr/=0) call fatal_error(h//"Allocating turb%tkep")
!TO modify to Xeon-Phi
!    if (.not. associated(turb%epsp)) allocate (turb%epsp(n1,n2,n3))
    if(idiffk(ng) == 5 .or. idiffk(ng) == 6)  & 
         allocate (turb%epsp(n1,n2,n3), STAT=ierr)
         if (ierr/=0) call fatal_error(h//"Allocating turb%epsp")
    
!srf if(idiffk(ng) == 7 .and. IMASSFLX == 1) & 
    if(idiffk(ng) == 7 ) & 
         allocate (turb%kpbl(n2,n3), STAT=ierr)
         if (ierr/=0) call fatal_error(h//"Allocating turb%kpbl")
    
    allocate (turb%hkm(n1,n2,n3), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating turb%hkm")
    allocate (turb%vkm(n1,n2,n3), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating turb%vkm")
    allocate (turb%vkh(n1,n2,n3), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating turb%vkh")

    if(if_urban_canopy > 0) &  
         allocate (turb%cdrag(n1,n2,n3), STAT=ierr)
         if (ierr/=0) call fatal_error(h//"Allocating turb%cdrag")
    
    allocate (turb%sflux_u(n2,n3), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating turb%sflux_u")
    allocate (turb%sflux_v(n2,n3), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating turb%sflux_v")
    allocate (turb%sflux_w(n2,n3), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating turb%sflux_w")
    allocate (turb%sflux_t(n2,n3), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating turb%sflux_t")
    allocate (turb%sflux_r(n2,n3), STAT=ierr)
    if (ierr/=0) call fatal_error(h//"Allocating turb%sflux_r")

  end subroutine alloc_turb

  ! ***********************************************************************

  subroutine nullify_turb(turb)

    implicit none
    type (turb_vars), intent(INOUT) :: turb

    nullify (turb%tkep)
    nullify (turb%epsp)
    nullify (turb%hkm)
    nullify (turb%vkm)
    nullify (turb%vkh)
    nullify (turb%cdrag)
    nullify (turb%sflux_u)
    nullify (turb%sflux_v)
    nullify (turb%sflux_w)
    nullify (turb%sflux_t)
    nullify (turb%sflux_r)
    nullify (turb%kpbl   )
  end subroutine nullify_turb

  ! ***********************************************************************

  subroutine dealloc_turb(turb)

    implicit none
    type (turb_vars), intent(INOUT) :: turb

    if (associated(turb%tkep))    deallocate (turb%tkep)
    if (associated(turb%epsp))    deallocate (turb%epsp)
    if (associated(turb%hkm))     deallocate (turb%hkm)
    if (associated(turb%vkm))     deallocate (turb%vkm)
    if (associated(turb%vkh))     deallocate (turb%vkh)
    if (associated(turb%cdrag))   deallocate (turb%cdrag)
    if (associated(turb%sflux_u)) deallocate (turb%sflux_u)
    if (associated(turb%sflux_v)) deallocate (turb%sflux_v)
    if (associated(turb%sflux_w)) deallocate (turb%sflux_w)
    if (associated(turb%sflux_t)) deallocate (turb%sflux_t)
    if (associated(turb%sflux_r)) deallocate (turb%sflux_r)
    if (associated(turb%kpbl    )) deallocate (turb%kpbl  )

  end subroutine dealloc_turb

  ! ***********************************************************************

  subroutine filltab_turb(turb, turbm, imean, n1, n2, n3, ng)
    use var_tables, only: InsertVTab
    implicit none
    include "i8.h"
    ! Arguments:
    type (turb_vars), intent(IN) :: turb, turbm
    integer, intent(IN) :: imean,n1,n2,n3,ng
    ! Local Variables:
    integer(kind=i8) :: npts

    ! Fill pointers to arrays into variable tables

    npts = n1*n2*n3

    if (associated(turb%tkep))  &
         call InsertVTab(turb%tkep, turbm%tkep,  &
         ng, npts, imean, 'TKEP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(turb%epsp))  &
         call InsertVTab(turb%epsp, turbm%epsp,  &
         ng, npts, imean, 'EPSP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(turb%hkm))  &
         call InsertVTab(turb%hkm, turbm%hkm,  &
         ng, npts, imean, 'HKM :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(turb%vkm))  &
         call InsertVTab(turb%vkm, turbm%vkm,  &
         ng, npts, imean, 'VKM :3:hist:mpti:mpt3:mpt1')
    if (associated(turb%vkh))  &
         call InsertVTab(turb%vkh, turbm%vkh,  &
         ng, npts, imean, 'VKH :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(turb%cdrag))  &
         call InsertVTab(turb%cdrag, turbm%cdrag,  &
         ng, npts, imean, 'CDRAG :3:hist:anal:mpti')

    npts = n2*n3
    if (associated(turb%sflux_u))  &
         call InsertVTab(turb%sflux_u, turbm%sflux_u,  &
         ng, npts, imean, 'SFLUX_U :2:anal:mpt3:mpt1')
    if (associated(turb%sflux_v))  &
         call InsertVTab(turb%sflux_v, turbm%sflux_v,  &
         ng, npts, imean, 'SFLUX_V :2:anal:mpt3:mpt1')
    if (associated(turb%sflux_w))  &
         call InsertVTab (turb%sflux_w, turbm%sflux_w,  &
         ng, npts, imean, 'SFLUX_W :2:anal:mpt3')
    if (associated(turb%sflux_t))  &
         call InsertVTab (turb%sflux_t, turbm%sflux_t,  &
         ng, npts, imean, 'SFLUX_T :2:anal:mpt3')
    if (associated(turb%sflux_r))  &
         call InsertVTab (turb%sflux_r, turbm%sflux_r,  &
         ng, npts, imean, 'SFLUX_R :2:anal:mpt3')
    if (associated(turb%kpbl))  &
         call InsertVTab (turb%kpbl, turbm%kpbl,  &
         ng, npts, imean, 'KPBL :2:anal:mpt3')


  end subroutine filltab_turb

  subroutine StoreNamelistFileAtMem_turb(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    akmin = oneNamelistFile%akmin
    csx = oneNamelistFile%csx
    csz = oneNamelistFile%csz
    idiffk = oneNamelistFile%idiffk
    if_urban_canopy = oneNamelistFile%if_urban_canopy
    ihorgrad = oneNamelistFile%ihorgrad
    xkhkm = oneNamelistFile%xkhkm
    zkhkm = oneNamelistFile%zkhkm
  end subroutine StoreNamelistFileAtMem_turb
end module mem_turb
