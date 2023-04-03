!############################# Change Log ##################################
! 3.1.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003, 2004 - All Rights Reserved
!  Brazilian Regional Atmospheric Modeling System - BRAMS
!###########################################################################

module mem_stilt
 
 use grid_dims
 use io_params, only: frqanl
 
!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  use ModNamelistFile, only: namelistFile

  include "i8.h"
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

 type stilt_vars
   real, pointer, dimension(:,:,:) :: &
                            thvlast,lnthvadv, lnthetav, lnthvtend, afxu, afxv, afxw &
                           ,ltscale, sigw, ltscaleb, sigwb &
                           ,tkepb, afxub, afxvb, afxwb &
                           ,cfxup1, cfxdn1, dfxup1, efxup1 &
                           ,dfxdn1, efxdn1, cfxup2, dfxup2 &
                           ,efxup2                         & 
!-srf : for the true air density
                           ,dnp

   real, pointer, dimension(:,:) :: pblhgt,lmo
 end type

 type (stilt_vars), allocatable :: stilt_g(:), stiltm_g(:)
 integer                        :: iexev,imassflx   
 
 real                          :: frqmassave
 !----- These variables control the time when the averages should be reset. -------------!
 real   , dimension(maxgrds)   :: etime_adve
 real   , dimension(maxgrds)   :: etime_turb

contains

 subroutine alloc_stilt(idiffk,stilt,n1,n2,n3,ng)
   implicit none
   type (stilt_vars) :: stilt
   integer, intent(in) :: n1, n2, n3, ng,idiffk
 
! Allocate arrays based on options (if necessary)
   if (iexev == 2) then
     allocate (stilt%thvlast(n1,n2,n3))
     allocate (stilt%lnthvadv(n1,n2,n3))
     allocate (stilt%lnthetav(n1,n2,n3))
     allocate (stilt%lnthvtend(n1,n2,n3))
!-srf : for the true air density
     allocate (stilt%dnp(n1,n2,n3))
   endif
   
   if (idiffk == 7) then
     allocate (stilt%ltscale(n1,n2,n3))
     allocate (stilt%sigw(n1,n2,n3))
     allocate (stilt%pblhgt(n2,n3))
!
     allocate (stilt%lmo(n2,n3))
   end if

   if (imassflx == 1) then
    if (idiffk == 7) then
      allocate (stilt%ltscale(n1,n2,n3))
      allocate (stilt%sigw(n1,n2,n3))
      allocate (stilt%pblhgt(n2,n3))
     end if
     allocate (stilt%afxu(n1,n2,n3))
     allocate (stilt%afxv(n1,n2,n3))
     allocate (stilt%afxw(n1,n2,n3))
     allocate (stilt%afxub(n1,n2,n3))
     allocate (stilt%afxvb(n1,n2,n3))
     allocate (stilt%afxwb(n1,n2,n3))
     allocate (stilt%cfxup1(n1,n2,n3))
     allocate (stilt%cfxdn1(n1,n2,n3))
     allocate (stilt%dfxup1(n1,n2,n3))
     allocate (stilt%efxup1(n1,n2,n3))
     allocate (stilt%dfxdn1(n1,n2,n3))
     allocate (stilt%efxdn1(n1,n2,n3))
     allocate (stilt%cfxup2(n1,n2,n3))
     allocate (stilt%dfxup2(n1,n2,n3))
     allocate (stilt%efxup2(n1,n2,n3))
     if (idiffk /= 2 .and. idiffk /= 3) then
       allocate (stilt%tkepb(n1,n2,n3))
     end if
     if (idiffk == 7) then
       allocate (stilt%ltscaleb(n1,n2,n3))
       allocate (stilt%sigwb(n1,n2,n3))
     end if
   end if

  !srf- initialization
  frqmassave =frqanl
  etime_adve =0.0
  etime_turb =0.0

 end subroutine alloc_stilt
 
 
 subroutine nullify_stilt(stilt)
   implicit none
   type (stilt_vars)   :: stilt

   if (associated(stilt%thvlast ))  nullify (stilt%thvlast )
   if (associated(stilt%lnthvadv ))  nullify (stilt%lnthvadv )
   if (associated(stilt%lnthetav ))  nullify (stilt%lnthetav )
   if (associated(stilt%lnthvtend ))  nullify (stilt%lnthvtend )
!-srf : for the true air density
   if (associated(stilt%dnp     ))  nullify (stilt%dnp     )

   if (associated(stilt%ltscale ))  nullify (stilt%ltscale )
   if (associated(stilt%sigw    ))  nullify (stilt%sigw    )
   if (associated(stilt%pblhgt  ))  nullify (stilt%pblhgt  )
   if (associated(stilt%lmo     ))  nullify (stilt%lmo     )
   if (associated(stilt%afxu    ))  nullify (stilt%afxu    )
   if (associated(stilt%afxv    ))  nullify (stilt%afxv    )
   if (associated(stilt%afxw    ))  nullify (stilt%afxw    )
   if (associated(stilt%ltscaleb))  nullify (stilt%ltscaleb)
   if (associated(stilt%sigwb   ))  nullify (stilt%sigwb   )
   if (associated(stilt%tkepb   ))  nullify (stilt%tkepb   ) 
   if (associated(stilt%afxub   ))  nullify (stilt%afxub   )
   if (associated(stilt%afxvb   ))  nullify (stilt%afxvb   )
   if (associated(stilt%afxwb   ))  nullify (stilt%afxwb   )
   if (associated(stilt%cfxup1  ))  nullify (stilt%cfxup1  )
   if (associated(stilt%cfxdn1  ))  nullify (stilt%cfxdn1  )
   if (associated(stilt%dfxup1  ))  nullify (stilt%dfxup1  )
   if (associated(stilt%efxup1  ))  nullify (stilt%efxup1  )
   if (associated(stilt%dfxdn1  ))  nullify (stilt%dfxdn1  )
   if (associated(stilt%efxdn1  ))  nullify (stilt%efxdn1  )
   if (associated(stilt%cfxup2  ))  nullify (stilt%cfxup2  )
   if (associated(stilt%dfxup2  ))  nullify (stilt%dfxup2  )
   if (associated(stilt%efxup2  ))  nullify (stilt%efxup2  )
   return
 end subroutine

 subroutine dealloc_stilt(stilt)
   implicit none
   type (stilt_vars)   :: stilt

   if (associated(stilt%thvlast ))  deallocate (stilt%thvlast )
   if (associated(stilt%lnthvadv ))  deallocate (stilt%lnthvadv )
   if (associated(stilt%lnthetav ))  deallocate (stilt%lnthetav ) 
   if (associated(stilt%lnthvtend ))  deallocate (stilt%lnthvtend )
!-srf : for the true air density
   if (associated(stilt%dnp     ))  deallocate (stilt%dnp     )

   if (associated(stilt%ltscale ))  deallocate (stilt%ltscale )
   if (associated(stilt%sigw    ))  deallocate (stilt%sigw    )
   if (associated(stilt%pblhgt  ))  deallocate (stilt%pblhgt  )
   if (associated(stilt%lmo     ))  deallocate (stilt%lmo     )
   if (associated(stilt%afxu    ))  deallocate (stilt%afxu    )
   if (associated(stilt%afxv    ))  deallocate (stilt%afxv    )
   if (associated(stilt%afxw    ))  deallocate (stilt%afxw    )
   if (associated(stilt%ltscaleb))  deallocate (stilt%ltscaleb)
   if (associated(stilt%sigwb   ))  deallocate (stilt%sigwb   )
   if (associated(stilt%tkepb   ))  deallocate (stilt%tkepb   )
   if (associated(stilt%afxub   ))  deallocate (stilt%afxub   )
   if (associated(stilt%afxvb   ))  deallocate (stilt%afxvb   )
   if (associated(stilt%afxwb   ))  deallocate (stilt%afxwb   )
   if (associated(stilt%cfxup1  ))  deallocate (stilt%cfxup1  )
   if (associated(stilt%cfxdn1  ))  deallocate (stilt%cfxdn1  )
   if (associated(stilt%dfxup1  ))  deallocate (stilt%dfxup1  )
   if (associated(stilt%efxup1  ))  deallocate (stilt%efxup1  )
   if (associated(stilt%dfxdn1  ))  deallocate (stilt%dfxdn1  )
   if (associated(stilt%efxdn1  ))  deallocate (stilt%efxdn1  )
   if (associated(stilt%cfxup2  ))  deallocate (stilt%cfxup2  )
   if (associated(stilt%dfxup2  ))  deallocate (stilt%dfxup2  )
   if (associated(stilt%efxup2  ))  deallocate (stilt%efxup2  )
   return
 end subroutine 

 subroutine filltab_stilt(stilt,stiltm,imean,n1,n2,n3,ng)
 
!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
    use var_tables, only: InsertVTab
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------
   
   implicit none
   type (stilt_vars)   :: stilt,stiltm
   integer             :: imean, n1, n2, n3, ng
   real, pointer       :: var, varm

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
    integer(kind=i8) :: npts
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!   integer :: npts
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   npts=n1*n2*n3

   if (associated(stilt%thvlast)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%thvlast,stiltm%thvlast &
                 ,ng, npts, imean, &
                 'THVLAST :3:hist:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%thvlast(1,1,1),stiltm%thvlast(1,1,1) &
!                 ,ng, npts, imean, &
!                 'THVLAST :3:hist:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------
		 
   if (associated(stilt%lnthvadv)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%lnthvadv,stiltm%lnthvadv &
                 ,ng, npts, imean, &
                 'lnthvadv :3:hist:mpti:mpt3:mpt1')		 
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%lnthvadv(1,1,1),stiltm%lnthvadv(1,1,1) &
!                 ,ng, npts, imean, &
!                 'lnthvadv :3:hist:mpti:mpt3:mpt1')		 
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------
		 
   if (associated(stilt%lnthetav)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%lnthetav,stiltm%lnthetav &
                 ,ng, npts, imean, &
                 'lnthetav :3:hist:mpti:mpt3:mpt1')	
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%lnthetav(1,1,1),stiltm%lnthetav(1,1,1) &
!                 ,ng, npts, imean, &
!                 'lnthetav :3:hist:mpti:mpt3:mpt1')	
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%lnthvtend)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%lnthvtend,stiltm%lnthvtend &
                 ,ng, npts, imean, &
                 'lnthvtend :3:hist:mpti:mpt3:mpt1')		 
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%lnthvtend(1,1,1),stiltm%lnthvtend(1,1,1) &
!                 ,ng, npts, imean, &
!                 'lnthvtend :3:hist:mpti:mpt3:mpt1')		 
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------
		  
!-srf : for the true air density
   if (associated(stilt%dnp)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%dnp,stiltm%dnp &
                 ,ng, npts, imean, &
!srf             'DNP :3:hist:anal:mpti:mpt3:mpt1')
                 'DNP :3:hist:anal:mpti:mpt3')

!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%dnp(1,1,1),stiltm%dnp(1,1,1) &
!                 ,ng, npts, imean, &
!                 'DNP :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%ltscale)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%ltscale,stiltm%ltscale &
                 ,ng, npts, imean, &
                 'TL :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%ltscale(1,1,1),stiltm%ltscale(1,1,1) &
!                 ,ng, npts, imean, &
!                 'TL :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%sigw)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%sigw,stiltm%sigw &
                 ,ng, npts, imean, &
                 'SIGW :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%sigw(1,1,1),stiltm%sigw(1,1,1) &
!                 ,ng, npts, imean, &
!                 'SIGW :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%afxu)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%afxu,stiltm%afxu &
                 ,ng, npts, imean, &
                 'AFXU :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%afxu(1,1,1),stiltm%afxu(1,1,1) &
!                 ,ng, npts, imean, &
!                 'AFXU :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%afxv)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%afxv,stiltm%afxv &
                 ,ng, npts, imean, &
                 'AFXV :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%afxv(1,1,1),stiltm%afxv(1,1,1) &
!                 ,ng, npts, imean, &
!                 'AFXV :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%afxw)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%afxw,stiltm%afxw &
                 ,ng, npts, imean, &
                 'AFXW :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%afxw(1,1,1),stiltm%afxw(1,1,1) &
!                 ,ng, npts, imean, &
!                 'AFXW :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%ltscaleb)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%ltscaleb,stiltm%ltscaleb &
                 ,ng, npts, imean, &
                 'TLB :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%ltscaleb(1,1,1),stiltm%ltscaleb(1,1,1) &
!                 ,ng, npts, imean, &
!                 'TLB :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%sigwb)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%sigwb,stiltm%sigwb &
                 ,ng, npts, imean, &
                 'SIGWB :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%sigwb(1,1,1),stiltm%sigwb(1,1,1) &
!                 ,ng, npts, imean, &
!                 'SIGWB :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%tkepb)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%tkepb,stiltm%tkepb &
                 ,ng, npts, imean, &
                 'TKEPB :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%tkepb(1,1,1),stiltm%tkepb(1,1,1) &
!                 ,ng, npts, imean, &
!                 'TKEPB :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%afxub)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%afxub,stiltm%afxub &
                 ,ng, npts, imean, &
                 'AFXUB :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%afxub(1,1,1),stiltm%afxub(1,1,1) &
!                 ,ng, npts, imean, &
!                 'AFXUB :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%afxvb)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%afxvb,stiltm%afxvb &
                 ,ng, npts, imean, &
                 'AFXVB :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%afxvb(1,1,1),stiltm%afxvb(1,1,1) &
!                 ,ng, npts, imean, &
!                 'AFXVB :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%afxwb)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%afxwb,stiltm%afxwb &
                 ,ng, npts, imean, &
                 'AFXWB :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%afxwb(1,1,1),stiltm%afxwb(1,1,1) &
!                 ,ng, npts, imean, &
!                 'AFXWB :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%cfxup1)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%cfxup1,stiltm%cfxup1 &
                 ,ng, npts, imean, &
                 'CFXUP1 :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%cfxup1(1,1,1),stiltm%cfxup1(1,1,1) &
!                 ,ng, npts, imean, &
!                 'CFXUP1 :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%cfxdn1)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%cfxdn1,stiltm%cfxdn1 &
                 ,ng, npts, imean, &
                 'CFXDN1 :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%cfxdn1(1,1,1),stiltm%cfxdn1(1,1,1) &
!                 ,ng, npts, imean, &
!                 'CFXDN1 :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%dfxup1)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%dfxup1,stiltm%dfxup1 &
                 ,ng, npts, imean, &
                 'DFXUP1 :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%dfxup1(1,1,1),stiltm%dfxup1(1,1,1) &
!                 ,ng, npts, imean, &
!                 'DFXUP1 :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%efxup1)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%efxup1,stiltm%efxup1 &
                 ,ng, npts, imean, &
                 'EFXUP1 :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%efxup1(1,1,1),stiltm%efxup1(1,1,1) &
!                 ,ng, npts, imean, &
!                 'EFXUP1 :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%dfxdn1)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%dfxdn1,stiltm%dfxdn1 &
                 ,ng, npts, imean, &
                 'DFXDN1 :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%dfxdn1(1,1,1),stiltm%dfxdn1(1,1,1) &
!                 ,ng, npts, imean, &
!                 'DFXDN1 :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%efxdn1)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%efxdn1,stiltm%efxdn1 &
                 ,ng, npts, imean, &
                 'EFXDN1 :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%efxdn1(1,1,1),stiltm%efxdn1(1,1,1) &
!                 ,ng, npts, imean, &
!                 'EFXDN1 :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%cfxup2)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%cfxup2,stiltm%cfxup2 &
                 ,ng, npts, imean, &
                 'CFXUP2 :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%cfxup2(1,1,1),stiltm%cfxup2(1,1,1) &
!                 ,ng, npts, imean, &
!                 'CFXUP2 :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%dfxup2)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%dfxup2,stiltm%dfxup2 &
                 ,ng, npts, imean, &
                 'DFXUP2 :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%dfxup2(1,1,1),stiltm%dfxup2(1,1,1) &
!                 ,ng, npts, imean, &
!                 'DFXUP2 :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%efxup2)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%efxup2,stiltm%efxup2 &
                 ,ng, npts, imean, &
                 'EFXUP2 :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%efxup2(1,1,1),stiltm%efxup2(1,1,1) &
!                 ,ng, npts, imean, &
!                 'EFXUP2 :3:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   npts=n2*n3
   if (associated(stilt%pblhgt)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%pblhgt,stiltm%pblhgt &
                 ,ng, npts, imean, &
                 'PBLHGT :2:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%pblhgt(1,1),stiltm%pblhgt(1,1) &
!                 ,ng, npts, imean, &
!                 'PBLHGT :2:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   if (associated(stilt%lmo)) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
      call InsertVTab (stilt%lmo,stiltm%lmo &
                 ,ng, npts, imean, &
                 'LMO    :2:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!      call vtables2 (stilt%lmo(1,1),stiltm%lmo(1,1) &
!                 ,ng, npts, imean, &
!                 'LMO    :2:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

   return
 end subroutine filltab_stilt
 
 !=======================================================================================!
 !=======================================================================================!
 !     This function finds the virtual temperature based on the temperature and mixing   !
 ! ratio. Two notes: 
 ! 1. It will use the condensation in case the total mixing ratio is provided.           !
 ! 2. This can be used for virtual potential temperature, just give potential tempera-   !
 !    ture instead of temperature. 
 !---------------------------------------------------------------------------------------!
real function virtt(temp,rvap,rtot)
  
  use rconstants, only: ep
  
  implicit none
  
  !----- Arguments--------------------------------------------------------------------!
  real, intent(in)           :: temp     ! Temperature         [    K]
  real, intent(in)           :: rvap     ! Vapour mixing ratio [kg/kg]
  real, intent(in)           :: rtot     ! Total mixing ratio  [kg/kg]
 !------------------------------------------------------------------------------------!

  virtt = temp * (1. + rvap / ep) / (1. + rtot)

  return

end function virtt
!=======================================================================================!
!=======================================================================================!
  subroutine zero_average_mass_adve(stilt)
      implicit none
      type (stilt_vars)   :: stilt

      if (associated(stilt%afxub   ))  stilt%afxub   = 0.0
      if (associated(stilt%afxvb   ))  stilt%afxvb   = 0.0
      if (associated(stilt%afxwb   ))  stilt%afxwb   = 0.0

      return
   end subroutine zero_average_mass_adve
!=======================================================================================!  
!==========================================================================================!
!    This subroutine flushes the average variables to zero whenever necessary. This is     !
! done at the output subroutine, right after the analysis (lite and full) are saved.       !
!------------------------------------------------------------------------------------------!
   subroutine zero_average_mass_turb(stilt)
      implicit none
      type (stilt_vars)   :: stilt

      if (associated(stilt%ltscaleb))  stilt%ltscaleb= 0.0
      if (associated(stilt%sigwb   ))  stilt%sigwb   = 0.0
      if (associated(stilt%tkepb   ))  stilt%tkepb   = 0.0

      return
   end subroutine zero_average_mass_turb

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  subroutine StoreNamelistFileAtMem_stilt(oneNamelistFile)
  implicit none
    type(namelistFile), pointer :: oneNamelistFile
    iexev = oneNamelistFile%iexev
    imassflx = oneNamelistFile%imassflx
  end subroutine StoreNamelistFileAtMem_stilt
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

end module mem_stilt
