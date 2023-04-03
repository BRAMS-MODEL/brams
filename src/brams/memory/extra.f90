!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
module Extras

  use ModNamelistFile, only: namelistFile

  include "i8.h"

  ! Used in CCATT

  integer :: na_extra2d, na_extra3d

  type ext2d
     real, pointer :: d2(:,:)
  end type ext2d

  type ext3d
     real, pointer :: d3(:,:,:)
  end type ext3d

  type(ext2d), allocatable :: extra2d(:,:), extra2dm(:,:)
  ! extrad3d(indice,ngrid)
  type(ext3d), allocatable :: extra3d(:,:), extra3dm(:,:)

contains

  subroutine alloc_extra2d(scal,m1,m2,na2d,ngrid)
    use dump, only: &
      dumpMessage

    implicit none
    include "constants.f90"
    ! Arguments:
    type (ext2d), intent(INOUT) :: scal(:,:)
    integer, intent(IN) :: m1,m2 !Dimension of arrays
    integer, intent(IN) :: na2d ! number of 2d extras arrays without ngrid
    integer, intent(IN) :: ngrid
    ! Local Variables:
    character(len=*), parameter :: h="**(alloc_extra2d)**"
    integer :: ierr
    integer :: j

    do j=1,na2d
       allocate(scal(j,ngrid)%d2(m1,m2), STAT=ierr)
       if (ierr/=0) &! call fatal_error(h//"Allocating scal%d2")
       iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
         "error allocating scal%d2")
    end do

  end subroutine alloc_extra2d

  !---------------------------------------------------------------

  subroutine alloc_extra3d(scal,m1,m2,m3,na3d,ngrid)
    use dump, only: &
      dumpMessage

    implicit none
    include "constants.f90"
    ! Arguments:
    type (ext3d), intent(INOUT) :: scal(:,:)
    integer, intent(IN) :: m1,m2,m3 !Dimension of arrays
    integer, intent(IN) :: na3d ! number of 2d extras arrays without ngrid
    integer, intent(IN) :: ngrid
    ! Local Variables:
    character(len=*), parameter :: h="**(alloc_extra3d)**"
    integer :: ierr
    integer :: j

    do j=1,na3d
       allocate(scal(j,ngrid)%d3(m1,m2,m3), STAT=ierr)
       if (ierr/=0) &! call fatal_error(h//"Allocating scal%d3")
       iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
         "error allocating scal%d3")
    end do

  end subroutine alloc_extra3d

  !---------------------------------------------------------------

  subroutine dealloc_extra2d(scal,na2d,ngrid)

    implicit none
    ! Arguments:
    type (ext2d), intent(INOUT) :: scal(:,:)
    integer, intent(in) :: ngrid, na2d
    ! Local Variables:
    integer :: nsc, nd

    do nd=1,na2d
       do nsc=1,ngrid
          if (associated(scal(nd,nsc)%d2))   deallocate (scal(nd,nsc)%d2)
       enddo
    end do

  end subroutine dealloc_extra2d

  !---------------------------------------------------------------

  subroutine dealloc_extra3d(scal,na3d,ngrid)

    implicit none
    ! Arguments:
    type (ext3d), intent(INOUT) :: scal(:,:)
    integer, intent(in) :: ngrid, na3d
    ! Local Variables:
    integer :: nsc, nd

    do nd=1,na3d
       do nsc=1,ngrid
          if (associated(scal(nd,nsc)%d3))   deallocate (scal(nd,nsc)%d3)
       enddo
    end do

  end subroutine dealloc_extra3d

  !---------------------------------------------------------------

  subroutine nullify_extra2d(scal,na2d,ngrid)

    implicit none
    ! Arguments:
    type (ext2d), intent(INOUT) :: scal(:,:)
    integer, intent(in) :: na2d, ngrid
    ! Local Variables:
    integer :: nsc, nd

    do nd=1,na2d
       do nsc=1,ngrid
          nullify (scal(nd,ngrid)%d2)
       enddo
    enddo

  end subroutine nullify_extra2d

  !---------------------------------------------------------------

  subroutine nullify_extra3d(scal,na3d,ngrid)

    implicit none
    ! Arguments:
    type (ext3d), intent(INOUT) :: scal(:,:)
    integer, intent(in) :: na3d, ngrid
    ! Local Variables:
    integer :: nsc, nd

    do nd=1,na3d
       do nsc=1,ngrid
          nullify (scal(nd,ngrid)%d3)
       enddo
    enddo

  end subroutine nullify_extra3d

  !---------------------------------------------------------------

  subroutine filltab_extra2d(scal2, scalm2, imean, n1, n2, ng, na)
    use var_tables

    implicit none
    include "i8.h"
    ! Arguments:
    type (ext2d), intent(IN) :: scal2, scalm2
    integer, intent(IN) :: imean, n1, n2, ng, na
    ! Local Variables:
    integer(kind=i8)  :: npts
    character (len=7) :: sname

    ! Fill pointers to arrays into variable tables

    if (associated(scal2%d2)) then
       npts = n1*n2
       write(sname, '(a2,i3.3)') 'd2', na
       call InsertVTab(scal2%d2, scalm2%d2,  &
            ng, npts, imean,  &
            trim(sname)//' :2:hist:anal:mpti:mpt3') ! Default - Column oriented Proc.

    endif

  end subroutine filltab_extra2d

  !---------------------------------------------------------------
  subroutine filltab_extra3d(scal3, scalm3, imean, n1, n2, n3, ng, na)
    use var_tables

    implicit none
    include "i8.h"
    ! Arguments:
    type (ext3d), intent(IN) :: scal3, scalm3
    integer, intent(IN) :: imean, n1, n2, n3, ng, na
    ! Local Variables:
    integer(kind=i8)  :: npts
    character (len=7) :: sname

    ! Fill pointers to arrays into variable tables

    if (associated(scal3%d3)) then
       npts = n1*n2*n3
       write(sname, '(a2,i3.3)') 'd3', na
       call InsertVTab(scal3%d3, scalm3%d3,  &
            ng, npts, imean,  &
            trim(sname)//' :3:hist:anal:mpti:mpt3') ! Default - Column oriented Proc.

    endif

  end subroutine filltab_extra3d

  !-----------------------------------------------------------------
  subroutine zero_extra3d(scal, na3d, ngrid)

    implicit none
    ! Arguments:
    type (ext3d), intent(INOUT) :: scal(:,:)
    integer, intent(IN) :: na3d, ngrid
    ! Local Variables:
    integer :: nsc

    do nsc=1,na3d
       scal(nsc,ngrid)%d3(:,:,:) = 0.
    enddo

  end subroutine zero_extra3d

  !-----------------------------------------------------------------
  subroutine zero_extra2d(scal, na2d, ngrid)

    implicit none
    ! Arguments:
    type (ext2d), intent(INOUT) :: scal(:,:)
    integer, intent(IN) :: na2d, ngrid
    ! Local Variables:
    integer :: nsc

    do nsc=1,na2d
       scal(nsc,ngrid)%d2(:,:) = 0.
    enddo

  end subroutine zero_extra2d


  subroutine StoreNamelistFileAtExtras(oneNamelistFile)
    type(namelistFile), pointer :: oneNamelistFile
    na_extra2d = oneNamelistFile%na_extra2d
    na_extra3d = oneNamelistFile%na_extra3d
  end subroutine StoreNamelistFileAtExtras

end module Extras
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------
