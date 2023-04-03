!###########################################################################
! CCATT- B - Regional Atmospheric Modeling System - RAMS
!###########################################################################

module mem_volc_chem1

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  use ModNamelistFile, only: namelistFile

  include "i8.h"
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

  type volc_mean_vars   
     real, pointer, dimension(:,:)  :: plum_heigth    
     real, pointer, dimension(:,:)  :: vent_elev     
     real, pointer, dimension(:,:)  :: duration     
     real, pointer, dimension(:,:)  :: begin_time     
  end type volc_mean_vars  

  type (volc_mean_vars)    , allocatable :: volc_mean_g(:), volc_meanm_g(:)
   
  integer:: volcanoes 

contains
  !---------------------------------------------------------------

  subroutine alloc_volc_chem1(volc_mean,n1,n2,n3)

    !use chem1_list, only : spc_alloc,spc_name,src,on
    implicit none

    integer,intent(in) :: n1,n2,n3
    type (volc_mean_vars) :: volc_mean
        
    allocate (volc_mean%plum_heigth(n2,n3));  volc_mean%plum_heigth(:,:)=0.
    allocate (volc_mean%vent_elev   (n2,n3)); volc_mean%vent_elev  (:,:)=0.
    allocate (volc_mean%duration   (n2,n3)); volc_mean%duration   (:,:)=0.
    !allocate (volc_mean%begin_time (n2,n3)); volc_mean%begin_time (:,:)=0.

  end subroutine alloc_volc_chem1

  !---------------------------------------------------------------
  !subroutine dealloc_volc_chem1(volc_mean)
  !
  ! implicit none
  !
  ! 
  !end subroutine dealloc_volc_chem1
  !---------------------------------------------------------------

  subroutine nullify_volc_chem1(volc_mean)

    implicit none
    type (volc_mean_vars) :: volc_mean

       if (associated(volc_mean%plum_heigth)) nullify(volc_mean%plum_heigth)
       if (associated(volc_mean%vent_elev  )) nullify(volc_mean%vent_elev   )
       if (associated(volc_mean%duration )) nullify(volc_mean%duration )
       !if (associated(volc_mean%begin_time )) nullify(volc_mean%begin_time )
     
  end subroutine nullify_volc_chem1

  !---------------------------------------------------------------

  subroutine filltab_volc_chem1(volc_mean,volc_meanm&
                          ,imean,n1,n2,n3,ng)

    !use chem1_list, only: spc_alloc,spc_name,src,on 
    !use mem_chem1, only: chem1_g

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
    use var_tables, only: InsertVTab
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

    implicit none

    integer, intent(in) :: imean,n1,n2,n3,ng
    type (volc_mean_vars) ::volc_mean,volc_meanm

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
    integer(kind=i8) :: npts
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!   integer npts
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

    ! Fill pointers to arrays into variable tables
    ! 2d var
    npts=n2*n3

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
    call InsertVTab(volc_mean%plum_heigth, volc_meanm%plum_heigth, &
   		    ng, npts, imean,                               &
                    'plum_heigth'//'_volc'//                       &
                    ' :2:hist:anal:mpti:mpt3:mpt1')
    call InsertVTab(volc_mean%vent_elev, volc_meanm%vent_elev, &
   		    ng, npts, imean, 'vent_elev'//'_volc'//    &
                    ' :2:hist:anal:mpti:mpt3:mpt1')
    call InsertVTab(volc_mean%vent_elev, volc_meanm%vent_elev, &
   		    ng, npts, imean, 'duration'//'_volc'//    &
                    ' :2:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!    call vtables2 (volc_mean%plum_heigth(1,1), volc_meanm%plum_heigth(1,1)&
!   		   ,ng, npts, imean, 'plum_heigth'//'_volc'//&
!        	   ' :2:hist:anal:mpti:mpt3:mpt1')
!    call vtables2 (volc_mean%vent_elev(1,1), volc_meanm%vent_elev(1,1)&
!   		   ,ng, npts, imean, 'vent_elev'//'_volc'//&
!        	   ' :2:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

  end subroutine filltab_volc_chem1
  !---------------------------------------------------------------


!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  subroutine StoreNamelistFileAtMem_volcChem1(oneNamelistFile)
    type(namelistFile), pointer :: oneNamelistFile
    volcanoes = oneNamelistFile%volcanoes
  end subroutine StoreNamelistFileAtMem_volcChem1
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------


end module mem_volc_chem1
