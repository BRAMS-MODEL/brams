!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_cuparm

  use ModNamelistFile, only: namelistFile

  use grid_dims, only: &
       maxgrds,        & ! INTENT(IN)
       maxfiles          ! INTENT(IN)

  type cuparm_vars

     ! Variables to be dimensioned by (nzp,nxp,nyp)
     real, pointer :: thsrc(:,:,:)
     real, pointer :: rtsrc(:,:,:)
     real, pointer :: thsrcf(:,:,:)
     real, pointer :: rtsrcf(:,:,:)
     real, pointer :: thsrcp(:,:,:)
     real, pointer :: rtsrcp(:,:,:)
	 real, pointer :: clsrc(:,:,:) !srf-cloud/ice source term

     ! Variables to be dimensioned by (nxp,nyp)
     real, pointer :: aconpr(:,:)
     real, pointer :: conprr(:,:)
     real, pointer :: conprrp(:,:)
     real, pointer :: conprrf(:,:)

  end type cuparm_vars

  type (cuparm_vars), allocatable :: cuparm_g(:)
  type (cuparm_vars), allocatable :: cuparmm_g(:)
  type (cuparm_vars), allocatable :: cuparm_g_sh(:)
  type (cuparm_vars), allocatable :: cuparmm_g_sh(:)

  include "files.h"

  integer, parameter :: maxcufiles = maxfiles
  integer, parameter :: maxcugrids = 10

  integer :: if_cuinv             ! from RAMSIN
  real :: tcu_beg                 ! from RAMSIN
  real :: tcu_end                 ! from RAMSIN
  real :: cu_til                  ! from RAMSIN
  real :: cu_tel                  ! from RAMSIN
  real :: tnudcu                  ! from RAMSIN
  real :: wt_cu_grid(maxcugrids)  ! from RAMSIN
  character(len=128) :: cu_prefix ! from RAMSIN

  character(len=f_name_length) :: fnames_cu(maxcufiles)
  character(len=14)  :: itotdate_cu(maxcufiles)
  real :: cu_times(maxcufiles)

  integer :: ncufiles
  integer :: ncufl
  real :: cutime1
  real :: cutime2

  integer :: nnqparm(maxgrds) ! from RAMSIN
  real :: wcldbs              ! from RAMSIN
  real :: confrq              ! from RAMSIN

  public :: hasAconpr

contains

  logical function hasAconpr(ng)
    implicit none
    integer, intent(IN) :: ng

    hasAconpr = nnqparm(ng)>= 1 .or. if_cuinv == 1

  end function hasAconpr

  subroutine alloc_cuparm(cuparm, n1, n2, n3, ng)

    implicit none
    ! Arguments:
    type (cuparm_vars), intent(INOUT) :: cuparm
    integer, intent(IN)               :: n1, n2, n3, ng
    ! Local Variables:
    character(len=*), parameter :: h="**(alloc_cuparm)**"
    integer :: ierr

    ! Allocate arrays based on options (if necessary)

    if( hasAconpr(ng) )  then
    	if(nnqparm(ng) < 3) then !srf: for conv 3 and 4, special arrays
                                !srf:  are allocated instead of these ones

       		allocate (cuparm%thsrc(n1,n2,n3), STAT=ierr)
       		if (ierr/=0) call fatal_error(h//"Allocating cuparm%thsrc")
			cuparm%thsrc=0.0
       		allocate (cuparm%rtsrc(n1,n2,n3), STAT=ierr)
       		if (ierr/=0) call fatal_error(h//"Allocating cuparm%rtsrc")
	   		cuparm%rtsrc=0.0
       		allocate (cuparm%clsrc(n1,n2,n3), STAT=ierr)
       		if (ierr/=0) call fatal_error(h//"Allocating cuparm%clsrc")
			cuparm%clsrc=0.0 !srf-cloud/ice source term

       endif
       allocate (cuparm%aconpr(n2,n3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating cuparm%aconpr")
       cuparm%aconpr = 0.
       allocate (cuparm%conprr(n2,n3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating cuparm%conprr")
       cuparm%conprr = 0.  
   
       if (if_cuinv == 1) then
          allocate (cuparm%thsrcp(n1,n2,n3), STAT=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating cuparm%thsrcp")
			cuparm%thsrcp=0.0
          allocate (cuparm%rtsrcp(n1,n2,n3), STAT=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating cuparm%rtsrcp")
			cuparm%rtsrcp=0.0
          allocate (cuparm%thsrcf(n1,n2,n3), STAT=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating cuparm%thsrcf")
			cuparm%thsrcf=0.0
          allocate (cuparm%rtsrcf(n1,n2,n3), STAT=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating cuparm%rtsrcf")
			cuparm%rtsrcf=0.0
          allocate (cuparm%conprrp(n2,n3), STAT=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating cuparm%conprp")
			cuparm%conprrp=0.0
          allocate (cuparm%conprrf(n2,n3), STAT=ierr)
          if (ierr/=0) call fatal_error(h//"Allocating cuparm%conprrf")
			cuparm%conprrf=0.0
       endif
    endif

  end subroutine alloc_cuparm

  ! ----------------------------------------------------------------------
  subroutine alloc_cuparm_sh(cuparm,n1,n2,n3,ng)

    implicit none
    type (cuparm_vars) :: cuparm
    integer, intent(in) :: n1,n2,n3,ng
    character(len=*), parameter :: h="**(alloc_cuparm)**"
    integer :: ierr
    ! Allocate arrays for shallow cum feedback

        allocate (cuparm%thsrc(n1,n2,n3), STAT=ierr)
       	if (ierr/=0) call fatal_error(h//"Allocating cuparm%thsrc")
		cuparm%thsrc=0.0
       	allocate (cuparm%rtsrc(n1,n2,n3), STAT=ierr)
       	if (ierr/=0) call fatal_error(h//"Allocating cuparm%rtsrc")
	   	cuparm%rtsrc=0.0

  end subroutine alloc_cuparm_sh

  ! ----------------------------------------------------------------------

  subroutine nullify_cuparm(cuparm)

    implicit none
    ! Arguments:
    type (cuparm_vars), intent(INOUT) :: cuparm

    if (associated(cuparm%thsrc))    nullify (cuparm%thsrc)
    if (associated(cuparm%rtsrc))    nullify (cuparm%rtsrc)
    if (associated(cuparm%clsrc))    nullify (cuparm%clsrc)
    if (associated(cuparm%thsrcp))    nullify (cuparm%thsrcp)
    if (associated(cuparm%rtsrcp))    nullify (cuparm%rtsrcp)
    if (associated(cuparm%thsrcf))    nullify (cuparm%thsrcf)
    if (associated(cuparm%rtsrcf))    nullify (cuparm%rtsrcf)
    if (associated(cuparm%aconpr))   nullify (cuparm%aconpr)
    if (associated(cuparm%conprr))   nullify (cuparm%conprr)
    if (associated(cuparm%conprrp))   nullify (cuparm%conprrp)
    if (associated(cuparm%conprrf))   nullify (cuparm%conprrf)

  end subroutine nullify_cuparm

  ! ----------------------------------------------------------------------

  subroutine dealloc_cuparm(cuparm)

    implicit none
    ! Arguments:
    type (cuparm_vars), intent(INOUT) :: cuparm

    if (associated(cuparm%thsrc))    deallocate (cuparm%thsrc)
    if (associated(cuparm%rtsrc))    deallocate (cuparm%rtsrc)
    if (associated(cuparm%clsrc))    deallocate (cuparm%clsrc)
    if (associated(cuparm%thsrcp))    deallocate (cuparm%thsrcp)
    if (associated(cuparm%rtsrcp))    deallocate (cuparm%rtsrcp)
    if (associated(cuparm%thsrcf))    deallocate (cuparm%thsrcf)
    if (associated(cuparm%rtsrcf))    deallocate (cuparm%rtsrcf)
    if (associated(cuparm%aconpr))   deallocate (cuparm%aconpr)
    if (associated(cuparm%conprr))   deallocate (cuparm%conprr)
    if (associated(cuparm%conprrp))   deallocate (cuparm%conprrp)
    if (associated(cuparm%conprrf))   deallocate (cuparm%conprrf)

  end subroutine dealloc_cuparm

  ! ----------------------------------------------------------------------

  subroutine filltab_cuparm_sh(cuparm, cuparmm, imean, n1, n2, n3, ng)
    use var_tables
    implicit none
    include "i8.h"
    ! Arguments:
    type (cuparm_vars), intent(IN) :: cuparm, cuparmm
    integer, intent(IN)            :: imean, n1, n2, n3, ng
    ! Local Variables:
    integer(kind=i8) :: npts

    ! Fill pointers to arrays into variable tables

    npts = n1*n2*n3

    if (associated(cuparm%thsrc))  &
         call InsertVTab(cuparm%thsrc, cuparmm%thsrc,  &
         ng, npts, imean, 'THSRC_SH :3:hist:anal:mpti:mpt3')
    if (associated(cuparm%rtsrc))  &
         call InsertVTab(cuparm%rtsrc, cuparmm%rtsrc,  &
         ng, npts, imean, 'RTSRC_SH :3:hist:anal:mpti:mpt3')

  end subroutine filltab_cuparm_sh

  ! ----------------------------------------------------------------------

  subroutine filltab_cuparm(cuparm, cuparmm, imean, n1, n2, n3, ng)
    use var_tables

    implicit none
    include "i8.h"
    ! Arguments:
    type (cuparm_vars), intent(IN) :: cuparm, cuparmm
    integer, intent(IN)            :: imean, n1, n2, n3, ng
    ! Local Variables:
    integer(kind=i8) :: npts

    ! Fill pointers to arrays into variable tables

    npts = n1*n2*n3

    if (associated(cuparm%thsrc))  &
         call InsertVTab(cuparm%thsrc, cuparmm%thsrc,  &
         ng, npts, imean, 'THSRC :3:hist:anal:mpti:mpt3')
    if (associated(cuparm%rtsrc))  &
         call InsertVTab(cuparm%rtsrc, cuparmm%rtsrc,  &
         ng, npts, imean, 'RTSRC :3:hist:anal:mpti:mpt3')
    if (associated(cuparm%clsrc))  &
         call InsertVTab(cuparm%clsrc,cuparmm%clsrc  &
         ,ng, npts, imean, 'CLSRC :3:hist:anal:mpti:mpt3')
    if (associated(cuparm%thsrcp))  &
         call InsertVTab(cuparm%thsrcp, cuparmm%thsrcp,  &
         ng, npts, imean, 'THSRCP :3:mpti:')
    if (associated(cuparm%rtsrcp))  &
         call InsertVTab(cuparm%rtsrcp, cuparmm%rtsrcp,  &
         ng, npts, imean, 'RTSRCP :3:mpti:')
    if (associated(cuparm%thsrcf))  &
         call InsertVTab(cuparm%thsrcf, cuparmm%thsrcf,  &
         ng, npts, imean, 'THSRCF :3:mpti:')
    if (associated(cuparm%rtsrcf))  &
         call InsertVTab(cuparm%rtsrcf, cuparmm%rtsrcf,  &
         ng, npts, imean, 'RTSRCF :3:mpti:')

    npts = n2*n3
    if (associated(cuparm%aconpr))  &
         call InsertVTab(cuparm%aconpr, cuparmm%aconpr,  &
         ng, npts, imean, 'ACONPR :2:hist:anal:mpti:mpt3')
    if (associated(cuparm%conprr))  &
         call InsertVTab(cuparm%conprr, cuparmm%conprr,  &
         ng, npts, imean, 'CONPRR :2:hist:anal:mpt3')
    if (associated(cuparm%conprrp))  &
         call InsertVTab(cuparm%conprrp, cuparmm%conprrp,  &
         ng, npts, imean, 'CONPRRP :2:mpti')
    if (associated(cuparm%conprrf))  &
         call InsertVTab(cuparm%conprrf, cuparmm%conprrf,  &
         ng, npts, imean, 'CONPRRF :2:mpti')

  end subroutine filltab_cuparm

  subroutine StoreNamelistFileAtMem_cuparm(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    confrq = oneNamelistFile%confrq
    cu_prefix = oneNamelistFile%cu_prefix
    cu_tel = oneNamelistFile%cu_tel
    cu_til = oneNamelistFile%cu_til
    if_cuinv = oneNamelistFile%if_cuinv
    nnqparm = oneNamelistFile%nnqparm
    tcu_beg = oneNamelistFile%tcu_beg
    tcu_end = oneNamelistFile%tcu_end
    tnudcu = oneNamelistFile%tnudcu
    wcldbs = oneNamelistFile%wcldbs
    wt_cu_grid = oneNamelistFile%wt_cu_grid
  end subroutine StoreNamelistFileAtMem_cuparm

end module mem_cuparm
