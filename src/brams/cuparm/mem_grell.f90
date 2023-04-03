! Module necessary to GRELL param.

module mem_grell

  type grell_vars

     ! Variables to be dimensioned by (m2,m3)
     real, pointer, dimension(:,:) :: &
          UPMF,                       &
          DNMF,                       &
          XIACT_C,                    &
          XIACT_P,                    &
          XIERR,                      &
          XKDT,                       &
          XKTOP,                      &
          XKBCON,                     &
          XJMIN,                      &
	  XK22					
	  
     ! Variables to be dimensioned by (m1,m2,m3)
     !real, pointer, dimension(:,:,:) :: &
          !lsfth, &
	  !lsfrt

  end type grell_vars

  type (grell_vars), allocatable :: grell_g   (:), grellm_g   (:)
  type (grell_vars), allocatable :: grell_g_sh(:), grellm_g_sh(:)


  type cuforc_vars
     real, pointer, dimension(:,:,:) :: &
          lsfth, &
	  lsfrt
  end type cuforc_vars

  type (cuforc_vars), allocatable :: cuforc_g   (:), cuforcm_g   (:)
  type (cuforc_vars), allocatable :: cuforc_sh_g(:), cuforcm_sh_g(:)



contains

  subroutine alloc_grell(grell, m1, m2, m3, ng)

    use mem_cuparm, only : nnqparm  ! INTENT(IN)

    implicit none
    ! Arguments
    type (grell_vars), intent(INOUT) :: grell
    integer, intent(IN)              :: m1, m2, m3, ng
    ! Local Variables:
    integer :: ierr
    character(len=*), parameter :: h="**(alloc_grell)**"

    ! Allocate arrays based on options (if necessary)
    if (abs(nnqparm(ng))==2)  then
       allocate (grell%UPMF     (m2, m3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating grell%UPMF")
       allocate (grell%DNMF     (m2, m3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating grell%DNMF")
       allocate (grell%XIACT_C  (m2, m3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating grell%XIACT_C")
       allocate (grell%XIACT_P  (m2, m3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating grell%XIACT_P")
       allocate (grell%XIERR    (m2, m3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating grell%XIERR")
       allocate (grell%XKDT     (m2, m3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating grell%XKDT")
       allocate (grell%XKTOP    (m2, m3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating grell%XKTOP")
       allocate (grell%XKBCON   (m2, m3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating grell%XKBCON")
       allocate (grell%XJMIN    (m2, m3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating grell%XJMIN")
       allocate (grell%XK22     (m2, m3), STAT=ierr)
       if (ierr/=0) call fatal_error(h//"Allocating grell%XK22")
  

!--(DMK-LFR NEC-SX6)----------------------------------------------
       grell%upmf = 0.
       grell%dnmf = 0.
       grell%xiact_c = 0.
       grell%xiact_p = 0.
       grell%xierr = 0.
       grell%xkdt = 0.
       grell%xktop = 0.
       grell%xkbcon = 0.
       grell%xjmin = 0.
       grell%xk22 = 0.
!--(DMK-LFR NEC-SX6)----------------------------------------------

    endif

  end subroutine alloc_grell
!---------------------------------------------------------------
!---------------------------------------------------------------
  subroutine alloc_grell_sh(grell, m1, m2, m3, ng)

    use mem_cuparm, only : nnqparm  ! INTENT(IN)
    use shcu_vars_const, only : nnshcu   ! INTENT(IN)

    implicit none
    type (grell_vars) :: grell
    integer, intent(in) :: m1, m2, m3, ng

    ! Allocate arrays based on options (if necessary)

       allocate (grell%UPMF     (m2, m3));grell%UPMF    =0.0
       allocate (grell%DNMF     (m2, m3));grell%DNMF    =0.0
       allocate (grell%XIACT_C  (m2, m3));grell%XIACT_C =0.0
       allocate (grell%XIACT_P  (m2, m3));grell%XIACT_P =0.0
       allocate (grell%XIERR    (m2, m3));grell%XIERR   =0.0
       allocate (grell%XKDT     (m2, m3));grell%XKDT    =0.0
       allocate (grell%XKTOP    (m2, m3));grell%XKTOP   =0.0
       allocate (grell%XKBCON   (m2, m3));grell%XKBCON  =0.0
       allocate (grell%XJMIN    (m2, m3));grell%XJMIN   =0.0
       allocate (grell%XK22     (m2, m3));grell%XK22	=0.0

    return
  end subroutine alloc_grell_sh
!---------------------------------------------------------------
!---------------------------------------------------------------
 subroutine alloc_cu_forcings(cuforc, m1, m2, m3, ng)

!    use mem_cuparm, only : nnqparm 
!    use shcu_vars_const, only : nnshcu   

    implicit none
    type (cuforc_vars) :: cuforc
    integer, intent(in) :: m1, m2, m3, ng
    
    allocate (cuforc%lsfth(m1, m2, m3));cuforc%lsfth=0.0
    allocate (cuforc%lsfrt(m1, m2, m3));cuforc%lsfrt=0.0

  end subroutine alloc_cu_forcings
!---------------------------------------------------------------
  subroutine nullify_cuforc(cuforc)

    implicit none
    type (cuforc_vars) :: cuforc

    if (associated(cuforc%lsfth))   nullify (cuforc%lsfth)
    if (associated(cuforc%lsfrt))   nullify (cuforc%lsfrt)

  end subroutine nullify_cuforc
!---------------------------------------------------------------
!---------------------------------------------------------------

  subroutine filltab_cuforc_sh(cuforc, cuforcm,imean, m1, m2, m3, ng)

    use var_tables, only: InsertVTab
    implicit none
    include "i8.h"
    
    type (cuforc_vars) :: cuforc, cuforcm
    integer, intent(in) :: imean, m1,  m2, m3, ng
    integer(kind=i8) :: npts
     npts=m1*m2*m3

     if (associated(cuforc%lsfth))  &
         call InsertVTab(cuforc%lsfth,cuforcm%lsfth,ng, npts, imean, 'LSFTH_SH :3:hist:anal:mpti:mpt3')

     if (associated(cuforc%lsfrt))  &
         call InsertVTab(cuforc%lsfrt,cuforcm%lsfrt,ng, npts, imean, 'LSFRT_SH :3:hist:anal:mpti:mpt3')
	 
  end subroutine filltab_cuforc_sh

!---------------------------------------------------------------
!---------------------------------------------------------------

  subroutine filltab_cuforc(cuforc, cuforcm,imean, m1, m2, m3, ng)

    use var_tables, only: InsertVTab
    implicit none
    include "i8.h"

   
    type (cuforc_vars) :: cuforc, cuforcm
    integer, intent(in) :: imean, m1,  m2, m3, ng
    integer(kind=i8) :: npts
     npts=m1*m2*m3

     if (associated(cuforc%lsfth))  &
         call InsertVTab (cuforc%lsfth,cuforcm%lsfth &
         ,ng, npts, imean,  &
         'LSFTH :3:hist:anal:mpti:mpt3')

     if (associated(cuforc%lsfrt))  &
         call InsertVTab (cuforc%lsfrt,cuforcm%lsfrt &
         ,ng, npts, imean,  &
         'LSFRT :3:hist:anal:mpti:mpt3')
  end subroutine filltab_cuforc

!---------------------------------------------------------------
!---------------------------------------------------------------

  subroutine nullify_grell(grell)

    implicit none
    type (grell_vars), intent(INOUT) :: grell

    nullify (grell%UPMF)
    nullify (grell%DNMF)
    nullify (grell%XIACT_C)
    nullify (grell%XIACT_P)
    nullify (grell%XIERR)
    nullify (grell%XKDT)
    nullify (grell%XKTOP)
    nullify (grell%XKBCON)
    nullify (grell%XJMIN)
    nullify (grell%XK22)


  end subroutine nullify_grell

  ! *************************************************************************

  subroutine dealloc_grell(grell)

    implicit none
    type (grell_vars), intent(INOUT) :: grell

    if (associated(grell%UPMF))    deallocate (grell%UPMF)
    if (associated(grell%DNMF))    deallocate (grell%DNMF)
    if (associated(grell%XIACT_C)) deallocate (grell%XIACT_C)
    if (associated(grell%XIACT_P)) deallocate (grell%XIACT_P)
    if (associated(grell%XIERR))   deallocate (grell%XIERR)
    if (associated(grell%XKDT))    deallocate (grell%XKDT)
    if (associated(grell%XKTOP))   deallocate (grell%XKTOP)
    if (associated(grell%XKBCON))  deallocate (grell%XKBCON)
    if (associated(grell%XJMIN))   deallocate (grell%XJMIN)
    if (associated(grell%XK22))    deallocate (grell%XK22)

  end subroutine dealloc_grell

  ! *************************************************************************

  subroutine filltab_grell(grell, grellm, imean, m1, m2, m3, ng)

     use mem_cuparm, only : nnqparm  
	 use var_tables, only: InsertVTab
    implicit none
    include "i8.h"
    ! Arguments:
    type (grell_vars), intent(IN) :: grell, grellm
    integer, intent(IN)           :: imean, m1,  m2, m3, ng
    ! Local Variables:
    integer(kind=i8) :: npts

    ! Fill pointers to arrays into variable tables

    npts = m2*m3

    if( nnqparm(ng) == 2 )  then
    
    if (associated(grell%UPMF))  &
         call InsertVTab(grell%UPMF, grellm%UPMF, &
         ng, npts, imean, 'UPMF :2:hist:anal:mpti:mpt3')

    if (associated(grell%DNMF))  &
         call InsertVTab(grell%DNMF, grellm%DNMF, &
         ng, npts, imean, 'DNMF :2:hist:anal:mpti:mpt3')

    if (associated(grell%XIACT_C))  &
         call InsertVTab(grell%XIACT_C, grellm%XIACT_C, &
         ng, npts, imean, 'XIACT_C :2:hist:anal:mpti:mpt3')

    if (associated(grell%XIACT_P))  &
         call InsertVTab(grell%XIACT_P, grellm%XIACT_P, &
         ng, npts, imean, 'XIACT_P :2:hist:anal:mpti:mpt3')

    if (associated(grell%XIERR))  &
         call InsertVTab(grell%XIERR, grellm%XIERR, &
         ng, npts, imean, 'XIERR :2:hist:anal:mpti:mpt3')

    if (associated(grell%XKDT))  &
         call InsertVTab(grell%XKDT, grellm%XKDT, &
         ng, npts, imean, 'XKDT :2:hist:anal:mpti:mpt3')

    if (associated(grell%XKTOP))  &
         call InsertVTab(grell%XKTOP, grellm%XKTOP, &
         ng, npts, imean, 'XKTOP :2:hist:anal:mpti:mpt3')

    if (associated(grell%XKBCON))  &
         call InsertVTab(grell%XKBCON, grellm%XKBCON, &
         ng, npts, imean, 'XKBCON :2:hist:anal:mpti:mpt3')

    if (associated(grell%XJMIN))  &
         call InsertVTab(grell%XJMIN, grellm%XJMIN, &
         ng, npts, imean, 'XJMIN :2:hist:anal:mpti:mpt3')

    if (associated(grell%XK22))  &
         call InsertVTab(grell%XK22, grellm%XK22, &
         ng, npts, imean, 'XK22 :2:hist:anal:mpti:mpt3')
    endif

  end subroutine filltab_grell

!---------------------------------------------------------------
!---------------------------------------------------------------
  subroutine filltab_grell_sh(grell_sh, grellm_sh, imean, m1, m2, m3, ng)

   	use var_tables, only: InsertVTab
    USE shcu_vars_const, ONLY: NNSHCU 

    implicit none
    include "i8.h"
    type (grell_vars) :: grell_sh, grellm_sh
    integer, intent(in) :: imean, m1,  m2, m3, ng
    integer(kind=i8) :: npts

    ! Fill pointers to arrays into variable tables

    if(NNSHCU(ng) == 1 .or. NNSHCU(ng) ==2) then
     npts=m2*m3
   
    if (associated(grell_sh%UPMF))  &
         call InsertVTab (grell_sh%UPMF,grellm_sh%UPMF &
         ,ng, npts, imean,  &
         'UPMFSH :2:hist:anal:mpti:mpt3')

    if (associated(grell_sh%XIERR))  &
         call InsertVTab (grell_sh%XIERR,grellm_sh%XIERR &
         ,ng, npts, imean,  &
         'XIERRSH :2:hist:anal:mpti:mpt3')

    if (associated(grell_sh%XKTOP))  &
         call InsertVTab (grell_sh%XKTOP,grellm_sh%XKTOP &
         ,ng, npts, imean,  &
         'XKTOPSH :2:hist:anal:mpti:mpt3')

    if (associated(grell_sh%XKBCON))  &
         call InsertVTab (grell_sh%XKBCON,grellm_sh%XKBCON &
         ,ng, npts, imean,  &
         'XKBCONSH :2:hist:mpti:mpt3')

    if (associated(grell_sh%XK22))  &
         call InsertVTab(grell_sh%XK22,grellm_sh%XK22 &
         ,ng, npts, imean,  &
         'XK22SH :2:hist:mpti:mpt3')

    endif
  end subroutine filltab_grell_sh

end module mem_grell
