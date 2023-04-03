!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_opt

  type opt_vars

     integer, pointer :: ind1_x_a(:,:,:), ind1_x_b(:,:,:)
     integer, pointer :: ind2_x_a(:,:,:), ind2_x_b(:,:,:)
     real   , pointer :: weight_x_a(:,:,:), weight_x_b(:,:,:)
     integer, pointer :: ind1_y_a(:,:,:), ind1_y_b(:,:,:)
     integer, pointer :: ind2_y_a(:,:,:), ind2_y_b(:,:,:)
     real   , pointer :: weight_y_a(:,:,:), weight_y_b(:,:,:)

  end type opt_vars

  type (opt_vars) :: opt

  integer :: jstep, istep

contains

  subroutine alloc_opt_scratch(proc_type,ngrs,nnzp,nnxp,nnyp,lnxp,lnyp)

    implicit none
    include "i8.h"
    ! Arguments:

    integer, intent(in) :: proc_type, ngrs
    integer, intent(in), dimension (ngrs) :: nnzp, nnxp, nnyp
    integer, intent(in) :: lnxp, lnyp

    ! Local Variables:

    integer :: maxx, maxy, maxz, ng !, npts

!!$    integer(kind=i8) :: npts

    !print*, 'enter alloc_opt_scratch'

    if(proc_type==1) then  ! Master in Parallel runs does not need this
       maxx = 1
       maxy = 1
       maxz = 1
    else
       maxx = lnxp
       maxy = lnyp
       maxz = 0
       do ng=1,ngrs
          ! Looking for the grid with minimum number of points in x and y
          ! and comparing with the values lnxp and lnyp
          ! Avoiding the use of a number of points bigger then the local grid
          maxx = min(maxx,nnxp(ng))
          maxy = min(maxy,nnyp(ng))
          ! Looking for the grid with maximum number of points in z
          maxz = max(maxz,nnzp(ng))
       enddo
    endif

    !print *, "alloc_opt_scratch: proc_type,maxx,maxy,maxz=", &
    !     proc_type,maxx,maxy,maxz

    istep = maxx
    jstep = maxy

    allocate (opt%ind1_x_a(maxz,maxx,maxy))
    allocate (opt%ind1_x_b(maxz,maxx,maxy))
    allocate (opt%ind2_x_a(maxz,maxx,maxy))
    allocate (opt%ind2_x_b(maxz,maxx,maxy))
    allocate (opt%weight_x_a(maxz,maxx,maxy))
    allocate (opt%weight_x_b(maxz,maxx,maxy))
    allocate (opt%ind1_y_a(maxz,maxx,maxy))
    allocate (opt%ind1_y_b(maxz,maxx,maxy))
    allocate (opt%ind2_y_a(maxz,maxx,maxy))
    allocate (opt%ind2_y_b(maxz,maxx,maxy))
    allocate (opt%weight_y_a(maxz,maxx,maxy))
    allocate (opt%weight_y_b(maxz,maxx,maxy))

    ! ALF - Putting zero in all variables
!!$    npts = maxz*maxx*maxy
!!$    call azero(npts, opt%ind1_x_a(1,1,1))
!!$    call azero(npts, opt%ind1_x_b(1,1,1))
!!$    call azero(npts, opt%ind2_x_a(1,1,1))
!!$    call azero(npts, opt%ind2_x_b(1,1,1))
!!$    call azero(npts, opt%weight_x_a(1,1,1))
!!$    call azero(npts, opt%weight_x_b(1,1,1))
!!$    call azero(npts, opt%ind1_y_a(1,1,1))
!!$    call azero(npts, opt%ind1_y_b(1,1,1))
!!$    call azero(npts, opt%ind2_y_a(1,1,1))
!!$    call azero(npts, opt%ind2_y_b(1,1,1))
!!$    call azero(npts, opt%weight_y_a(1,1,1))
!!$    call azero(npts, opt%weight_y_b(1,1,1))
    opt%ind1_x_a   = 0.
    opt%ind1_x_b   = 0.
    opt%ind2_x_a   = 0.
    opt%ind2_x_b   = 0.
    opt%weight_x_a = 0.
    opt%weight_x_b = 0.
    opt%ind1_y_a   = 0.
    opt%ind1_y_b   = 0.
    opt%ind2_y_a   = 0.
    opt%ind2_y_b   = 0.
    opt%weight_y_a = 0.
    opt%weight_y_b = 0.

  end subroutine alloc_opt_scratch

  !---------------------------------------------------------------

  subroutine nullify_opt_scratch()

    implicit none

    ! Deallocate all scratch arrays

    if (associated(opt%ind1_x_a  ))  nullify (opt%ind1_x_a  )
    if (associated(opt%ind1_x_b  ))  nullify (opt%ind1_x_b  )
    if (associated(opt%ind2_x_a  ))  nullify (opt%ind2_x_a  )
    if (associated(opt%ind2_x_b  ))  nullify (opt%ind2_x_b  )
    if (associated(opt%weight_x_a))  nullify (opt%weight_x_a)
    if (associated(opt%weight_x_b))  nullify (opt%weight_x_b)
    if (associated(opt%ind1_y_a  ))  nullify (opt%ind1_y_a  )
    if (associated(opt%ind1_y_b  ))  nullify (opt%ind1_y_b  )
    if (associated(opt%ind2_y_a  ))  nullify (opt%ind2_y_a  )
    if (associated(opt%ind2_y_b  ))  nullify (opt%ind2_y_b  )
    if (associated(opt%weight_y_a))  nullify (opt%weight_y_a)
    if (associated(opt%weight_y_b))  nullify (opt%weight_y_b)


    return
  end subroutine nullify_opt_scratch
  !---------------------------------------------------------------

  subroutine dealloc_opt_scratch()

    implicit none

    ! Deallocate all scratch arrays

    if (associated(opt%ind1_x_a  ))  deallocate (opt%ind1_x_a  )
    if (associated(opt%ind1_x_b  ))  deallocate (opt%ind1_x_b  )
    if (associated(opt%ind2_x_a  ))  deallocate (opt%ind2_x_a  )
    if (associated(opt%ind2_x_b  ))  deallocate (opt%ind2_x_b  )
    if (associated(opt%weight_x_a))  deallocate (opt%weight_x_a)
    if (associated(opt%weight_x_b))  deallocate (opt%weight_x_b)
    if (associated(opt%ind1_y_a  ))  deallocate (opt%ind1_y_a  )
    if (associated(opt%ind1_y_b  ))  deallocate (opt%ind1_y_b  )
    if (associated(opt%ind2_y_a  ))  deallocate (opt%ind2_y_a  )
    if (associated(opt%ind2_y_b  ))  deallocate (opt%ind2_y_b  )
    if (associated(opt%weight_y_a))  deallocate (opt%weight_y_a)
    if (associated(opt%weight_y_b))  deallocate (opt%weight_y_b)

    return
  end subroutine dealloc_opt_scratch

end module mem_opt
