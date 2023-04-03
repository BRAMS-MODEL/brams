! Module necessary to Shallow Cumulus param.

module mem_shcu

  type shcu_vars

     ! Variables to be dimensioned by (nzp,nxp,nyp)
     real, pointer, dimension(:,:,:) :: &
          THSRCSH, RTSRCSH

     ! Variables to be dimensioned by (nxp,nyp)
     real, pointer, dimension(:,:) :: &
          SHMF

  end type shcu_vars

  type (shcu_vars), allocatable :: shcu_g(:), shcum_g(:)

contains

  subroutine alloc_shcu(shcu,n1,n2,n3,ng)

    implicit none
    type (shcu_vars) :: shcu
    integer, intent(in) :: n1,n2,n3,ng

    ! INCLUDE 'rcommons.h'

    ! Allocate arrays based on options (if necessary)

!    IF( nnqparm(ng)>= 1 )  THEN
       allocate (shcu%THSRCSH(n1,n2,n3))
       allocate (shcu%RTSRCSH(n1,n2,n3))
       allocate (shcu%SHMF(n2,n3))
       
!--(DMK-LFR NEC-SX6)----------------------------------------------
       shcu%thsrcsh = 0.
       shcu%rtsrcsh = 0.
       shcu%shmf = 0.
!--(DMK-LFR NEC-SX6)----------------------------------------------       
       
!    ENDIF

    return
  end subroutine alloc_shcu


  subroutine nullify_shcu(shcu)

    implicit none
    type (shcu_vars) :: shcu


    if (associated(shcu%THSRCSH))  nullify (shcu%THSRCSH)
    if (associated(shcu%RTSRCSH))  nullify (shcu%RTSRCSH)
    if (associated(shcu%SHMF))     nullify (shcu%SHMF)

    return
  end subroutine nullify_shcu

  subroutine dealloc_shcu(shcu)

    implicit none
    type (shcu_vars) :: shcu


    if (associated(shcu%THSRCSH))  deallocate (shcu%THSRCSH)
    if (associated(shcu%RTSRCSH))  deallocate (shcu%RTSRCSH)
    if (associated(shcu%SHMF))     deallocate (shcu%SHMF)

    return
  end subroutine dealloc_shcu


  subroutine filltab_shcu(shcu,shcum,imean,n1,n2,n3,ng)
    use var_tables, only: InsertVTab
    implicit none
    include "i8.h"
    type (shcu_vars) :: shcu, shcum
    integer, intent(in) :: imean,n1,n2,n3,ng
    integer(kind=i8) :: npts
    real, pointer :: var,varm

    ! Fill pointers to arrays into variable tables

    npts=n1*n2*n3

    if (associated(shcu%THSRCSH))  &
         call InsertVTab (shcu%THSRCSH,shcum%THSRCSH &
         ,ng, npts, imean,  &
         'THSRCSH :3:hist:anal:mpti:mpt3')
    if (associated(shcu%RTSRCSH))  &
         call InsertVTab (shcu%RTSRCSH,shcum%RTSRCSH &
         ,ng, npts, imean,  &
         'RTSRCSH :3:hist:anal:mpti:mpt3')

    npts=n2*n3
    if (associated(shcu%SHMF))  &
         call InsertVTab (shcu%SHMF,shcum%SHMF  &
         ,ng, npts, imean,  &
         'SHMF :2:hist:anal:mpti:mpt3')

    return
  end subroutine filltab_shcu

end module mem_shcu
