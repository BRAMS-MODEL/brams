!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


Module mem_basic

!!$  use MPI_IO_ENGINE, only: &
!!$       WRITE_BRAMS_DATA

  use ModNamelistFile, only:   namelistfile
 
  use mem_grid, only: dyncore_flag
 
  implicit none

  Type basic_vars
   
     ! Variables to be dimensioned by (nzp,nxp,nyp)
     real, pointer, dimension(:,:,:) :: &
          up,uc,vp,vc,wp,wc,pp,pc  &
          ,rv,theta,thp,thc,rtp &
          ,pi0,th0,dn0,dn0u,dn0v
                         
     ! Variables to be dimensioned by (nxp,nyp)
     real, pointer, dimension(:,:) :: &
          fcoru,fcorv,cputime

  End Type basic_vars
   
  type (basic_vars), allocatable :: basic_g(:), basicm_g(:)

  
Contains

  subroutine alloc_basic(basic, n1, n2, n3, ng)
    
    implicit none

    type (basic_vars), intent(inout) :: basic
    integer, intent(in) :: n1, n2, n3, ng
    
    integer :: ierr
    
    ! Allocate arrays based on options (if necessary)
    
    allocate (basic%up(n1,n2,n3), STAT=ierr);    IF (ierr/=0) CALL fatal_error("ERROR allocating up (alloc_basic)")
    allocate (basic%uc(n1,n2,n3), STAT=ierr);    IF (ierr/=0) CALL fatal_error("ERROR allocating uc (alloc_basic)")
    allocate (basic%vp(n1,n2,n3), STAT=ierr);    IF (ierr/=0) CALL fatal_error("ERROR allocating vp (alloc_basic)")
    allocate (basic%vc(n1,n2,n3), STAT=ierr);    IF (ierr/=0) CALL fatal_error("ERROR allocating vc (alloc_basic)")
    allocate (basic%wp(n1,n2,n3), STAT=ierr);    IF (ierr/=0) CALL fatal_error("ERROR allocating wp (alloc_basic)")
    allocate (basic%wc(n1,n2,n3), STAT=ierr);    IF (ierr/=0) CALL fatal_error("ERROR allocating wc (alloc_basic)")
    allocate (basic%pp(n1,n2,n3), STAT=ierr);    IF (ierr/=0) CALL fatal_error("ERROR allocating pp (alloc_basic)")
    allocate (basic%pc(n1,n2,n3), STAT=ierr);    IF (ierr/=0) CALL fatal_error("ERROR allocating pc (alloc_basic)")
    allocate (basic%thp(n1,n2,n3), STAT=ierr);   IF (ierr/=0) CALL fatal_error("ERROR allocating thp (alloc_basic)")
    allocate (basic%rtp(n1,n2,n3), STAT=ierr);   IF (ierr/=0) CALL fatal_error("ERROR allocating rtp (alloc_basic)")
    allocate (basic%rv(n1,n2,n3),  STAT=ierr);   IF (ierr/=0) CALL fatal_error("ERROR allocating rv (alloc_basic)")
    allocate (basic%theta(n1,n2,n3), STAT=ierr); IF (ierr/=0) CALL fatal_error("ERROR allocating theta (alloc_basic)")
    if ( dyncore_flag == 2 .or. dyncore_flag == 3 ) then
      allocate (basic%thc(n1,n2,n3), STAT=ierr); IF (ierr/=0) CALL fatal_error("ERROR allocating thc (alloc_basic)")
    end if
    allocate (basic%pi0(n1,n2,n3), STAT=ierr);   IF (ierr/=0) CALL fatal_error("ERROR allocating pi0 (alloc_basic)")
    allocate (basic%th0(n1,n2,n3), STAT=ierr);   IF (ierr/=0) CALL fatal_error("ERROR allocating th0 (alloc_basic)")
    allocate (basic%dn0(n1,n2,n3), STAT=ierr);   IF (ierr/=0) CALL fatal_error("ERROR allocating dn0 (alloc_basic)")
    allocate (basic%dn0u(n1,n2,n3), STAT=ierr);  IF (ierr/=0) CALL fatal_error("ERROR allocating dn0u (alloc_basic)")
    allocate (basic%dn0v(n1,n2,n3), STAT=ierr);  IF (ierr/=0) CALL fatal_error("ERROR allocating dn0v (alloc_basic)")
    
    allocate (basic%fcoru(n2,n3), STAT=ierr);    IF (ierr/=0) CALL fatal_error("ERROR allocating fcoru (alloc_basic)")
    allocate (basic%fcorv(n2,n3), STAT=ierr);    IF (ierr/=0) CALL fatal_error("ERROR allocating fcorv (alloc_basic)")
    allocate (basic%cputime(n2,n3), STAT=ierr);  IF (ierr/=0) CALL fatal_error("ERROR allocating cputime (alloc_basic)")

    ! Putting zeros on cputime array
    basic%cputime = 0.

!--(DMK-LFR NEC-SX6)----------------------------------------------
    basic%up    = 0.
    basic%uc    = 0.
    basic%vp    = 0.
    basic%vc    = 0.
    basic%wp    = 0.
    basic%wc    = 0.
    basic%pp    = 0.
    basic%pc    = 0.
    basic%thp   = 0.
    basic%rtp   = 0.
    basic%rv    = 0.
    basic%theta = 0.
    if ( dyncore_flag == 2 .or. dyncore_flag == 3 ) then
      basic%thc = 0.
    end if
    basic%pi0   = 0.
    basic%th0   = 0.
    basic%dn0   = 0.
    basic%dn0u  = 0.
    basic%dn0v  = 0.
    basic%fcoru = 0.
    basic%fcorv = 0.    
!--(DMK-LFR NEC-SX6)----------------------------------------------
   
  end subroutine alloc_basic
  
  
  subroutine nullify_basic(basic)
    
    implicit none
    type (basic_vars), intent(inout) :: basic
    
    if (associated(basic%up   ))    nullify (basic%up   )
    if (associated(basic%uc   ))    nullify (basic%uc   )
    if (associated(basic%vp   ))    nullify (basic%vp   )
    if (associated(basic%vc   ))    nullify (basic%vc   )
    if (associated(basic%wp   ))    nullify (basic%wp   )
    if (associated(basic%wc   ))    nullify (basic%wc   )
    if (associated(basic%pp   ))    nullify (basic%pp   )
    if (associated(basic%pc   ))    nullify (basic%pc   )
    if (associated(basic%thp  ))    nullify (basic%thp  )
    if (associated(basic%rtp  ))    nullify (basic%rtp  )
    if (associated(basic%rv   ))    nullify (basic%rv   )
    if (associated(basic%theta))    nullify (basic%theta)
    !if ( dyncore_flag == 2 .or. dyncore_flag == 3 ) then
      if (associated(basic%thc))    nullify (basic%thc)
    !end if
    if (associated(basic%pi0  ))    nullify (basic%pi0  )
    if (associated(basic%th0  ))    nullify (basic%th0  )
    if (associated(basic%dn0  ))    nullify (basic%dn0  )
    if (associated(basic%dn0u ))    nullify (basic%dn0u )
    if (associated(basic%dn0v ))    nullify (basic%dn0v )
    if (associated(basic%fcoru ))   nullify (basic%fcoru )
    if (associated(basic%fcorv ))   nullify (basic%fcorv )
    if (associated(basic%cputime )) nullify (basic%cputime )
    
  end subroutine nullify_basic
  

  subroutine dealloc_basic(basic)
    
    implicit none
    type (basic_vars), intent(inout) :: basic
    
    if (associated(basic%up   ))    deallocate (basic%up   )
    if (associated(basic%uc   ))    deallocate (basic%uc   )
    if (associated(basic%vp   ))    deallocate (basic%vp   )
    if (associated(basic%vc   ))    deallocate (basic%vc   )
    if (associated(basic%wp   ))    deallocate (basic%wp   )
    if (associated(basic%wc   ))    deallocate (basic%wc   )
    if (associated(basic%pp   ))    deallocate (basic%pp   )
    if (associated(basic%pc   ))    deallocate (basic%pc   )
    if (associated(basic%thp  ))    deallocate (basic%thp  )
    if (associated(basic%rtp  ))    deallocate (basic%rtp  )
    if (associated(basic%rv   ))    deallocate (basic%rv   )
    if (associated(basic%theta))    deallocate (basic%theta)
    !if ( dyncore_flag == 2 ) then
       if (associated(basic%thc))    deallocate (basic%thc)
    !end if	
    if (associated(basic%pi0  ))    deallocate (basic%pi0  )
    if (associated(basic%th0  ))    deallocate (basic%th0  )
    if (associated(basic%dn0  ))    deallocate (basic%dn0  )
    if (associated(basic%dn0u ))    deallocate (basic%dn0u )
    if (associated(basic%dn0v ))    deallocate (basic%dn0v )
    if (associated(basic%fcoru ))   deallocate (basic%fcoru )
    if (associated(basic%fcorv ))   deallocate (basic%fcorv )
    if (associated(basic%cputime )) deallocate (basic%cputime )
    
  end subroutine dealloc_basic
  
  
  subroutine filltab_basic(basic, basicm, imean, n1, n2, n3, ng)

    use var_tables, only: InsertVTab
    use mem_stilt, only: iexev
    implicit none
    include "i8.h"
    type (basic_vars), intent(in) :: basic, basicm
    integer, intent(in)           :: imean, n1, n2, n3, ng
    integer(kind=i8) :: npts
    real, pointer :: var, varm
    
    ! Fill pointers to arrays into variable tables
    
    npts = n1*n2*n3
    
    if (associated(basic%up))  &
         call InsertVTab (basic%up,basicm%up  &
         ,ng, npts, imean,  &
         'UP :3:hist:anal:mpti:mpt3:mpt2')
    if (associated(basic%vp))  &
         call InsertVTab (basic%vp,basicm%vp  &
         ,ng, npts, imean,  &
         'VP :3:hist:anal:mpti:mpt3:mpt2')
    if (associated(basic%wp))  &
         call InsertVTab (basic%wp,basicm%wp  &
         ,ng, npts, imean,  &
         'WP :3:hist:anal:mpti:mpt3:mpt2')
    if (associated(basic%pp))  &
         call InsertVTab (basic%pp,basicm%pp  &
         ,ng, npts, imean,  &
         'PP :3:hist:anal:mpti:mpt3:mpt2')
    if (associated(basic%uc))  &
         call InsertVTab (basic%uc,basicm%uc  &
         ,ng, npts, imean,  &
         'UC :3:hist:mpti:mpt3:mpt2')
    if (associated(basic%vc))  &
         call InsertVTab (basic%vc,basicm%vc  &
         ,ng, npts, imean,  &
         'VC :3:hist:mpti:mpt3:mpt2')
    if (associated(basic%wc))  &
         call InsertVTab (basic%wc,basicm%wc  &
         ,ng, npts, imean,  &
         'WC :3:hist:mpti:mpt3:mpt2')
    if (associated(basic%pc))  &
         call InsertVTab (basic%pc,basicm%pc  &
         ,ng, npts, imean,  &
         'PC :3:hist:mpti:mpt3:mpt2')
    
    if (associated(basic%thp)) &
         call InsertVTab (basic%thp,basicm%thp  &
         ,ng, npts, imean,  &
         'THP :3:hist:mpti:mpt3:mpt1')
    if (associated(basic%rtp)) &
         call InsertVTab (basic%rtp,basicm%rtp  &
         ,ng, npts, imean,  &
         'RTP :3:hist:mpti:mpt3:mpt1')
    
    if(iexev == 2) then
       if (associated(basic%theta)) &
         call InsertVTab (basic%theta,basicm%theta  &
         ,ng, npts, imean,  &
         'THETA :3:hist:anal:mpti:mpt3:mpt1')
    else
       if (associated(basic%theta)) &
         call InsertVTab (basic%theta,basicm%theta  &
         ,ng, npts, imean,  &
         'THETA :3:hist:anal:mpti:mpt3')
    endif
    
    !if ( dyncore_flag == 2 ) then
    if (associated(basic%thc)) &
         call InsertVTab (basic%thc,basicm%thc  &
         ,ng, npts, imean,  &
         'THC :3:hist:mpti:mpt3:mpt1')
    !endif
    
    if(iexev == 2) then
       if (associated(basic%rv)) &
         call InsertVTab (basic%rv,basicm%rv  &
         ,ng, npts, imean,  &
         'RV :3:hist:anal:mpti:mpt3:mpt1')
    else
       if (associated(basic%rv)) &
         call InsertVTab (basic%rv,basicm%rv  &
         ,ng, npts, imean,  &
         'RV :3:hist:anal:mpti:mpt3')
    endif
       
    if (associated(basic%pi0)) &
         call InsertVTab (basic%pi0,basicm%pi0  &
         ,ng, npts, imean,  &
         'PI0 :3:mpti')
    if (associated(basic%th0)) &
         call InsertVTab (basic%th0,basicm%th0  &
         ,ng, npts, imean,  &
         'TH0 :3:mpti')
    if (associated(basic%dn0)) &
         call InsertVTab (basic%dn0,basicm%dn0  &
         ,ng, npts, imean,  &
         'DN0 :3:mpti')
    if (associated(basic%dn0u)) &
         call InsertVTab (basic%dn0u,basicm%dn0u  &
         ,ng, npts, imean,  &
         'DN0U :3:mpti')
    if (associated(basic%dn0v)) &
         call InsertVTab (basic%dn0v,basicm%dn0v  &
         ,ng, npts, imean,  &
         'DN0V :3:mpti')
    
    npts = n2*n3

    if (associated(basic%fcoru)) &
         call InsertVTab (basic%fcoru,basicm%fcoru  &
         ,ng, npts, imean,  &
         'FCORU :2:mpti')      
    if (associated(basic%fcorv)) &
         call InsertVTab (basic%fcorv,basicm%fcorv  &
         ,ng, npts, imean,  &
         'FCORV :2:mpti')
    if (associated(basic%cputime)) &
         call InsertVTab (basic%cputime,basicm%cputime  &
         ,ng, npts, imean,  &
         'CPUTIME :2:anal:mpti:mpt3')
    
  end subroutine filltab_basic


!!$  subroutine save_inv_basic(basic, fileTypeWrite, fileTypeWriteGlobal)
!!$
!!$    implicit none
!!$    type (basic_vars), intent(in) :: basic
!!$    integer, intent(in)           :: fileTypeWrite, fileTypeWriteGlobal
!!$
!!$    call WRITE_BRAMS_DATA(basic%thp, fileTypeWrite, fileTypeWriteGlobal)
!!$    call WRITE_BRAMS_DATA(basic%rtp, fileTypeWrite, fileTypeWriteGlobal)
!!$    call WRITE_BRAMS_DATA(basic%rv, fileTypeWrite, fileTypeWriteGlobal)
!!$    call WRITE_BRAMS_DATA(basic%theta, fileTypeWrite, fileTypeWriteGlobal)
!!$    if ( dyncore_flag == 2 ) then
!!$      call WRITE_BRAMS_DATA(basic%thc, fileTypeWrite, fileTypeWriteGlobal)
!!$    end if
!!$    call WRITE_BRAMS_DATA(basic%pi0, fileTypeWrite, fileTypeWriteGlobal)
!!$    call WRITE_BRAMS_DATA(basic%th0, fileTypeWrite, fileTypeWriteGlobal)
!!$    call WRITE_BRAMS_DATA(basic%dn0, fileTypeWrite, fileTypeWriteGlobal)
!!$    call WRITE_BRAMS_DATA(basic%dn0u, fileTypeWrite, fileTypeWriteGlobal)
!!$    call WRITE_BRAMS_DATA(basic%dn0v, fileTypeWrite, fileTypeWriteGlobal)    
!!$    call WRITE_BRAMS_DATA(basic%fcoru, fileTypeWrite, fileTypeWriteGlobal)
!!$    call WRITE_BRAMS_DATA(basic%fcorv, fileTypeWrite, fileTypeWriteGlobal)
!!$
!!$  end subroutine save_inv_basic
!!$
!!$
!!$  subroutine save_his_basic(basic, fileTypeWrite, fileTypeWriteGlobal)
!!$
!!$    implicit none
!!$
!!$    type (basic_vars), intent(in) :: basic
!!$    integer, intent(in)           :: fileTypeWrite, fileTypeWriteGlobal
!!$
!!$    call WRITE_BRAMS_DATA(basic%up, fileTypeWrite, fileTypeWriteGlobal)
!!$    call WRITE_BRAMS_DATA(basic%uc, fileTypeWrite, fileTypeWriteGlobal)
!!$    call WRITE_BRAMS_DATA(basic%vp, fileTypeWrite, fileTypeWriteGlobal)
!!$    call WRITE_BRAMS_DATA(basic%vc, fileTypeWrite, fileTypeWriteGlobal)
!!$    call WRITE_BRAMS_DATA(basic%wp, fileTypeWrite, fileTypeWriteGlobal)
!!$    call WRITE_BRAMS_DATA(basic%wc, fileTypeWrite, fileTypeWriteGlobal)
!!$    call WRITE_BRAMS_DATA(basic%pp, fileTypeWrite, fileTypeWriteGlobal)
!!$    call WRITE_BRAMS_DATA(basic%pc, fileTypeWrite, fileTypeWriteGlobal)
!!$
!!$  end subroutine save_his_basic
!!$
!!$
!!$  subroutine save_all_basic(basic, fileTypeWrite, fileTypeWriteGlobal)
!!$
!!$    implicit none
!!$
!!$    type (basic_vars), intent(in) :: basic
!!$    integer, intent(in)           :: fileTypeWrite, fileTypeWriteGlobal
!!$
!!$    call save_inv_basic(basic, fileTypeWrite, fileTypeWriteGlobal)
!!$    call save_his_basic(basic, fileTypeWrite, fileTypeWriteGlobal)
!!$
!!$  end subroutine save_all_basic
  
End Module mem_basic
