!############################# Change Log ##################################
! 5.0.0
!
! Demerval S. Moreira - 27/Ago/2012
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


Module mem_jules

  use ModNamelistFile, only: namelistFile

  use grid_dims, only: maxgrds ! INTENT(IN)

  Type jules_vars

     ! Variables to be dimensioned by (nxp,nyp)

     !real, pointer :: u_s(:,:,:)    ! DSM 11/04/2014
     real, pointer :: gpp(:,:)      
     real, pointer :: resp_s(:,:)   
     real, pointer :: resp_p(:,:)   
     real, pointer :: npp(:,:)      
     real, pointer :: u10mj(:,:)
     real, pointer :: v10mj(:,:)
     real, pointer :: u10mj1hr(:,:,:)
     real, pointer :: v10mj1hr(:,:,:)
     real, pointer :: fracj(:,:,:)
     real, pointer :: t2mj(:,:)
!LFR
     real, pointer :: t2mj_max(:,:)
     real, pointer :: t2mj_min(:,:)
!LFR

     real, pointer :: rv2mj(:,:)
     real, pointer :: csj(:,:)
     real, pointer :: anthrop_heatj(:,:,:)
     real, pointer :: radnet_tilej(:,:,:)
     real, pointer :: ftl_tilej(:,:,:)
     real, pointer :: le_tilej(:,:,:)
     real, pointer :: htf_tilej(:,:,:)
     real, pointer :: snowdepthj(:,:,:)
     real, pointer :: ht_fluxj(:,:)
     real, pointer :: temp_surfj(:,:)

  End Type jules_vars

  type (jules_vars), allocatable :: jules_g(:), julesm_g(:)
  logical :: firstDayHour

Contains

  subroutine alloc_jules(jules,nz,nx,ny,nzg,nzs,np,ng)

    implicit none
    type (jules_vars) :: jules
    integer, intent(in) :: nz,nx,ny,nzg,nzs,np,ng

    ! Allocate arrays based on options (if necessary)

    !allocate (jules%u_s       (nx,ny,np))  ! DSM 11/04/2014
    allocate (jules%gpp       (nx,ny))    
    allocate (jules%resp_s    (nx,ny))   
    allocate (jules%resp_p    (nx,ny))   
    allocate (jules%npp       (nx,ny))  
    allocate (jules%u10mj     (nx,ny))
    allocate (jules%v10mj     (nx,ny))
    allocate (jules%u10mj1hr  (nz,nx,ny))
    allocate (jules%v10mj1hr  (nz,nx,ny))
    allocate (jules%fracj     (nz,nx,ny))
    allocate (jules%t2mj      (nx,ny))
    allocate (jules%t2mj_max  (nx,ny))
    allocate (jules%t2mj_min  (nx,ny))
    allocate (jules%rv2mj     (nx,ny))
    allocate (jules%csj       (nx,ny))
    allocate (jules%anthrop_heatj  (nx,ny,np))
    allocate (jules%radnet_tilej   (nx,ny,np))
    allocate (jules%ftl_tilej      (nx,ny,np))
    allocate (jules%le_tilej       (nx,ny,np))
    allocate (jules%htf_tilej      (nx,ny,np))
    allocate (jules%snowdepthj     (nx,ny,np))
    allocate (jules%ht_fluxj       (nx,ny))
    allocate (jules%temp_surfj     (nx,ny))

  end subroutine alloc_jules

  !************************************************************************

  subroutine nullify_jules(jules)

    implicit none
    type (jules_vars) :: jules

    !if(associated(jules%u_s))        nullify (jules%u_s)    ! DSM 10/04/2014
    if(associated(jules%gpp))        nullify (jules%gpp)   
    if(associated(jules%resp_s))     nullify (jules%resp_s) 
    if(associated(jules%resp_p))     nullify (jules%resp_p) 
    if(associated(jules%npp))        nullify (jules%npp)   
    if(associated(jules%u10mj))      nullify (jules%u10mj)
    if(associated(jules%v10mj))      nullify (jules%v10mj)
    if(associated(jules%u10mj1hr))   nullify (jules%u10mj1hr)
    if(associated(jules%v10mj1hr))   nullify (jules%v10mj1hr)
    if(associated(jules%fracj))      nullify (jules%fracj)
    if(associated(jules%t2mj))       nullify (jules%t2mj)
    if(associated(jules%t2mj_max))   nullify (jules%t2mj_max)
    if(associated(jules%t2mj_min))   nullify (jules%t2mj_min)
    if(associated(jules%rv2mj))      nullify (jules%rv2mj)
    if(associated(jules%csj))        nullify (jules%csj)
    if(associated(jules%anthrop_heatj))   nullify (jules%anthrop_heatj)
    if(associated(jules%radnet_tilej ))   nullify (jules%radnet_tilej )
    if(associated(jules%ftl_tilej    ))   nullify (jules%ftl_tilej    )
    if(associated(jules%le_tilej     ))   nullify (jules%le_tilej     )
    if(associated(jules%htf_tilej    ))   nullify (jules%htf_tilej    )
    if(associated(jules%snowdepthj   ))   nullify (jules%snowdepthj   )
    if(associated(jules%ht_fluxj     ))   nullify (jules%ht_fluxj     )
    if(associated(jules%temp_surfj   ))   nullify (jules%temp_surfj   )

  end subroutine nullify_jules

  ! ********************************************************************

  subroutine dealloc_jules(jules)

    implicit none
    type (jules_vars) :: jules

    !if(associated(jules%u_s))        deallocate (jules%u_s)   ! DSM 10/04/2014
    if(associated(jules%gpp))        deallocate (jules%gpp)    
    if(associated(jules%npp))        deallocate (jules%npp)   
    if(associated(jules%resp_s))     deallocate (jules%resp_s) 
    if(associated(jules%resp_p))     deallocate (jules%resp_p) 
    if(associated(jules%u10mj))      deallocate (jules%u10mj)
    if(associated(jules%v10mj))      deallocate (jules%v10mj)
    if(associated(jules%u10mj1hr))   deallocate (jules%u10mj1hr)
    if(associated(jules%v10mj1hr))   deallocate (jules%v10mj1hr)
    if(associated(jules%fracj))      deallocate (jules%fracj)
    if(associated(jules%t2mj))       deallocate (jules%t2mj)
    if(associated(jules%t2mj_max))       deallocate (jules%t2mj_max)    
    if(associated(jules%t2mj_min))       deallocate (jules%t2mj_min)
    if(associated(jules%rv2mj))      deallocate (jules%rv2mj)
    if(associated(jules%csj))        deallocate (jules%csj)
    if(associated(jules%anthrop_heatj))   deallocate (jules%anthrop_heatj)
    if(associated(jules%radnet_tilej ))   deallocate (jules%radnet_tilej )
    if(associated(jules%ftl_tilej    ))   deallocate (jules%ftl_tilej    )
    if(associated(jules%le_tilej     ))   deallocate (jules%le_tilej     )
    if(associated(jules%htf_tilej    ))   deallocate (jules%htf_tilej    )
    if(associated(jules%snowdepthj   ))   deallocate (jules%snowdepthj   )
    if(associated(jules%ht_fluxj     ))   deallocate (jules%ht_fluxj     )
    if(associated(jules%temp_surfj   ))   deallocate (jules%temp_surfj   )

  end subroutine dealloc_jules

  ! ********************************************************************

  subroutine filltab_jules(jules,julesm,imean,nz,nx,ny,nzg,nzs,np,ng)
    ! ALF
    use io_params, only: ipastin ! INTENT(IN)
    use var_tables, only: InsertVTab

    implicit none
    include "i8.h"
    type (jules_vars) :: jules,julesm
    integer, intent(in) :: imean,nz,nx,ny,nzg,nzs,np,ng
    integer(kind=i8) :: npts
    real, pointer :: var,varm
    ! ALF
    character(len=8) :: str_recycle

    ! ALF
    str_recycle = ''
    if (ipastin == 1) then
       str_recycle = ':recycle'
    endif

    ! Fill pointers to arrays into variable tables

    !npts=nx*ny*np
    !call InsertVTab (jules%u_s,julesm%u_s  &            
    !     ,ng, npts, imean, 'U_S :6:anal:mpti:mpt3')    

    npts=nx*ny
    call InsertVTab (jules%gpp,julesm%gpp  &           
         ,ng, npts, imean, 'GPP :2:anal:mpti:mpt3')  
	 
    call InsertVTab (jules%resp_s,julesm%resp_s  &    
         ,ng, npts, imean, 'RESP_S :2:anal:mpti:mpt3') 
	 
    call InsertVTab (jules%resp_p,julesm%resp_p  &     
         ,ng, npts, imean, 'RESP_P :2:anal:mpti:mpt3') 
	 
    call InsertVTab (jules%npp,julesm%npp  &           
         ,ng, npts, imean, 'NPP :2:anal:mpti:mpt3')    

    call InsertVTab (jules%u10mj,julesm%u10mj  &
         ,ng, npts, imean, 'U10MJ :2:anal:mpti:mpt3')

    call InsertVTab (jules%v10mj,julesm%v10mj  &
         ,ng, npts, imean, 'V10MJ :2:anal:mpti:mpt3')

    call InsertVTab (jules%t2mj,julesm%t2mj  &
         ,ng, npts, imean, 'T2MJ :2:anal:mpti:mpt3')

    call InsertVTab (jules%t2mj_max,julesm%t2mj_max  &
         ,ng, npts, imean, 'T2MJ_MAX :2:anal:mpti:mpt3')

    call InsertVTab (jules%t2mj_min,julesm%t2mj_min  &
         ,ng, npts, imean, 'T2MJ_MIN :2:anal:mpti:mpt3')

    call InsertVTab (jules%rv2mj,julesm%rv2mj  &
         ,ng, npts, imean, 'RV2MJ :2:anal:mpti:mpt3')

    call InsertVTab (jules%csj,julesm%csj  &
         ,ng, npts, imean, 'CSJ :2:anal:mpti:mpt3')

    call InsertVTab (jules%ht_fluxj,julesm%ht_fluxj  &
         ,ng, npts, imean, 'ht_fluxj :2:anal:mpti:mpt3')

    call InsertVTab (jules%temp_surfj,julesm%temp_surfj  &
         ,ng, npts, imean, 'temp_surfj :2:anal:mpti:mpt3')

    npts=nz*nx*ny
    call InsertVTab (jules%u10mj1hr,julesm%u10mj1hr  &
         ,ng, npts, imean, 'U10MJ1hr :3:hist:anal:mpti:mpt3:mpt2')

    call InsertVTab (jules%v10mj1hr,julesm%v10mj1hr  &
         ,ng, npts, imean, 'V10MJ1hr :3:hist:anal:mpti:mpt3:mpt2')

    call InsertVTab (jules%fracj,julesm%fracj  &
         ,ng, npts, imean, 'fracj :3:hist:anal:mpti:mpt3:mpt2')

    npts=nx*ny*np
    call InsertVTab (jules%anthrop_heatj,julesm%anthrop_heatj  &
         ,ng, npts, imean,  &
         'anthrop_heatj :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (jules%radnet_tilej,julesm%radnet_tilej  &
         ,ng, npts, imean,  &
         'radnet_tilej :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (jules%ftl_tilej,julesm%ftl_tilej  &
         ,ng, npts, imean,  &
         'ftl_tilej :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (jules%le_tilej,julesm%le_tilej  &
         ,ng, npts, imean,  &
         'le_tilej :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (jules%htf_tilej,julesm%htf_tilej  &
         ,ng, npts, imean,  &
         'htf_tilej :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (jules%snowdepthj,julesm%snowdepthj  &
         ,ng, npts, imean,  &
         'snowdepthj :6:hist:anal:mpti:mpt3'//trim(str_recycle))

  end subroutine filltab_jules

End Module mem_jules
