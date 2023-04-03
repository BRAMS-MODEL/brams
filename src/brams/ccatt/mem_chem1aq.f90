!###########################################################################
!  B - Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_chem1aq
  use grid_dims, only : maxgrds

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  use ModNamelistFile, only: namelistFile

  include "i8.h"
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

  type chem1aq_vars   
!--- All families
     real, pointer, dimension(:,:,:)  :: sc_pr,sc_pc ! r:rain, c:cloud
     real, pointer, dimension(:    )  :: sc_tr,sc_tc
!-----------
  end type chem1aq_vars
  type (chem1aq_vars)    , allocatable :: chem1aq_g(:,:) , chem1maq_g(:,:)
  
  integer :: CHEMISTRY_AQ

contains
  !---------------------------------------------------------------

  subroutine alloc_chem1aq(chem1aq,n1,n2,n3,nspeciesaq)

    use chem1aq_list, only : spcaq_name
    
    implicit none

    integer,intent(in) :: n1,n2,n3,nspeciesaq
    integer :: ispcaq
    
    type (chem1aq_vars)  ,dimension(     nspeciesaq) :: chem1aq
    
    !print*,'----------------------------------------------------------------'
    !print*,' memory allocation for aqueous chemical species:'
    
    do ispcaq=1,nspeciesaq
     !print*,'spc=',spcaq_name(ispcaq),'size=',n1,n2,n3
     
     !- allocate memory for the past time tracer mixing ratio
     allocate (chem1aq(ispcaq)%sc_pr  (n1,n2,n3));chem1aq(ispcaq)%sc_pr = 0.
     allocate (chem1aq(ispcaq)%sc_pc  (n1,n2,n3));chem1aq(ispcaq)%sc_pc = 0.

    enddo
    return
  end subroutine alloc_chem1aq

  !--------------------------------------------------------------------------

  subroutine dealloc_chem1aq(chem1aq,nspeciesaq)

   implicit none

   integer,intent(in) :: nspeciesaq
   type (chem1aq_vars)  ,dimension( nspeciesaq) :: chem1aq
   integer :: ispcaq
  
    !  Deallocate arrays
    do ispcaq=1,nspeciesaq

     if (associated(chem1aq(ispcaq)%sc_pr )) deallocate(chem1aq(ispcaq)%sc_pr  )
     if (associated(chem1aq(ispcaq)%sc_pc )) deallocate(chem1aq(ispcaq)%sc_pc  )
     !if (associated(chem1(ispc)%sc_s ))    deallocate(chem1(ispc)%sc_s  )
    enddo

    return
  end subroutine dealloc_chem1aq

  !---------------------------------------------------------------
  !---------------------------------------------------------------

  subroutine nullify_chem1aq(chem1aq,nspeciesaq)

    implicit none

    integer,intent(in) :: nspeciesaq
    type (chem1aq_vars),dimension(nspeciesaq) :: chem1aq
    integer :: ispcaq

    do ispcaq=1,nspeciesaq

     if (associated(chem1aq(ispcaq)%sc_pr )) nullify (chem1aq(ispcaq)%sc_pr )
     if (associated(chem1aq(ispcaq)%sc_pc )) nullify (chem1aq(ispcaq)%sc_pc )
     !if (associated(chem1(ispc)%sc_s ))    nullify (chem1(ispc)%sc_s  )
    enddo
    return
  end subroutine nullify_chem1aq

  !---------------------------------------------------------------

  subroutine filltab_chem1aq(chem1aq,chem1maq,imean,n1,n2,n3,nspeciesaq,ng)

!    use var_tables
    use chem1aq_list, only: spcaq_name

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
    use var_tables, only: InsertVTab
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

    implicit none

    integer, intent(in) :: imean,n1,n2,n3,nspeciesaq,ng
    type (chem1aq_vars)  ,dimension(   nspeciesaq) :: chem1aq,chem1maq

    integer :: ispcaq  

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
    integer(kind=i8) :: npts
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------
    
    character(len=8) :: str_recycle
   
    str_recycle = ''

    !- Fill pointers to arrays into variable tables
    do ispcaq=1,nspeciesaq

     if (associated(chem1aq(ispcaq)%sc_pr)) then
!--- tracer mixing ratio (dimension 3d)
       npts = n1 * n2 * n3

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
        call InsertVTab(chem1aq(ispcaq)%sc_pr, chem1maq(ispcaq)%sc_pr,  &
                        ng, npts, imean,                                &
                        trim(spcaq_name(ispcaq)) //'PR :3:hist:anal:mpti:mpt3:mpt1'//trim(str_recycle))
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!        call vtables2 (chem1aq(ispcaq)%sc_pr(1,1,1), chem1maq(ispcaq)%sc_pr(1,1,1)  &
!         ,ng, npts, imean, trim(spcaq_name(ispcaq)) //'PR :3:hist:anal:mpti:mpt3:mpt1'//trim(str_recycle))
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

     end if
!     
     if (associated(chem1aq(ispcaq)%sc_pc)) then
!--- tracer mixing ratio (dimension 3d)
       npts = n1 * n2 * n3

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
        call InsertVTab(chem1aq(ispcaq)%sc_pc, chem1maq(ispcaq)%sc_pc,  &
                        ng, npts, imean,                                &
                        trim(spcaq_name(ispcaq)) //'PC :3:hist:anal:mpti:mpt3:mpt1'//trim(str_recycle))
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!        call vtables2 (chem1aq(ispcaq)%sc_pc(1,1,1), chem1maq(ispcaq)%sc_pc(1,1,1)  &
!         ,ng, npts, imean, trim(spcaq_name(ispcaq)) //'PC :3:hist:anal:mpti:mpt3:mpt1'//trim(str_recycle))
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

     end if

    enddo
  end subroutine filltab_chem1aq


  !---------------------------------------------------------------
subroutine alloc_tend_chem1aq(nmzp,nmxp,nmyp,ngrs,nspeciesaq,proc_type)

   implicit none
   integer,intent(in)                   :: ngrs,proc_type,nspeciesaq
   integer,intent(in), dimension (ngrs) :: nmzp,nmxp,nmyp
   integer :: ng,ntpts,ispcaq

!         Find the maximum number of grid points needed for any grid.

   if(proc_type==1) then
      ntpts=1
   else
      ntpts=0
      do ng=1,ngrs
         ntpts=max( nmxp(ng)*nmyp(ng)*nmzp(ng),ntpts )
      enddo
   endif

    do ispcaq=1,nspeciesaq

     
      if (associated(chem1aq_g(ispcaq,1)%sc_pr)) allocate (chem1aq_g(ispcaq,1)%sc_tr(ntpts)) 
      if (associated(chem1aq_g(ispcaq,1)%sc_pc)) allocate (chem1aq_g(ispcaq,1)%sc_tc(ntpts))
      do ng=2,ngrs
         chem1aq_g(ispcaq,ng)%sc_tr => chem1aq_g(ispcaq,1)%sc_tr
         chem1aq_g(ispcaq,ng)%sc_tc => chem1aq_g(ispcaq,1)%sc_tc
      enddo
      
    enddo

  end subroutine alloc_tend_chem1aq
  !---------------------------------------------------------------

  subroutine nullify_tend_chem1aq(nspeciesaq)

    implicit none
    integer,intent(in) :: nspeciesaq
    integer ::ispcaq

    do ispcaq=1,nspeciesaq
      if (associated(chem1aq_g(ispcaq,1)%sc_tr)) nullify (chem1aq_g(ispcaq,1)%sc_tr)
      if (associated(chem1aq_g(ispcaq,1)%sc_tc)) nullify (chem1aq_g(ispcaq,1)%sc_tc)
    enddo


  end subroutine nullify_tend_chem1aq

  !---------------------------------------------------------------
  subroutine dealloc_tend_chem1aq(nspeciesaq)
    implicit none
    integer,intent(in) :: nspeciesaq
    integer ::ispcaq

    do ispcaq=1,nspeciesaq
     if (associated(chem1aq_g(ispcaq,1)%sc_tr)) deallocate (chem1aq_g(ispcaq,1)%sc_tr)
     if (associated(chem1aq_g(ispcaq,1)%sc_tc)) deallocate (chem1aq_g(ispcaq,1)%sc_tc)
    enddo

  end  subroutine dealloc_tend_chem1aq
   
  !---------------------------------------------------------------

  subroutine filltab_tend_chem1aq(nspeciesaq,ng)
    use chem1aq_list, only:spcaq_name
    use mem_chem1, only: nspecies_transported ! this is first calculated at chemistry 
                                              ! "filltab_tend_chem1" routine
    implicit none

    integer,intent(in) :: nspeciesaq,ng
    integer ::ispcaq
    integer :: elements

    
    do ispcaq=1,nspeciesaq

! Fill pointers to scalar arrays into scalar tables 


      if ( associated(chem1aq_g(ispcaq,ng)%sc_tr)) then
      	call vtables_scalar (chem1aq_g(ispcaq,ng)%sc_pr(1,1,1),&
        chem1aq_g(ispcaq,ng)%sc_tr(1),ng,trim(spcaq_name(ispcaq))//'PR')
        
	elements = size(chem1aq_g(ispcaq,ng)%sc_tr)
        
        call vtables_scalar_new (chem1aq_g(ispcaq,ng)%sc_pr(1,1,1),&
        chem1aq_g(ispcaq,ng)%sc_tr(1),ng, trim(spcaq_name(ispcaq))//'PR',elements)

	!- total number of transported species (CHEM + CHEM_AQ)
	nspecies_transported = nspecies_transported + 1 		


      endif
!      
      if ( associated(chem1aq_g(ispcaq,ng)%sc_tc)) then
    	call vtables_scalar (chem1aq_g(ispcaq,ng)%sc_pc(1,1,1),&
       chem1aq_g(ispcaq,ng)%sc_tc(1),ng,trim(spcaq_name(ispcaq))//'PC')
       
	elements = size(chem1aq_g(ispcaq,ng)%sc_tc)
        
        call vtables_scalar_new (chem1aq_g(ispcaq,ng)%sc_pc(1,1,1),&
        chem1aq_g(ispcaq,ng)%sc_tc(1),ng, trim(spcaq_name(ispcaq))//'PC',elements)

	!- total number of transported species (CHEM + CHEM_AQ)
	nspecies_transported = nspecies_transported + 1 		


      endif
!      
    enddo



  end subroutine filltab_tend_chem1aq
  

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  subroutine StoreNamelistFileAtMem_chem1aq(oneNamelistFile)
    type(namelistFile), pointer :: oneNamelistFile
    chemistry_aq  = oneNamelistFile%chemistry_aq
  end subroutine StoreNamelistFileAtMem_chem1aq
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------


end module mem_chem1aq

!--------------------------------------------------------------------------

