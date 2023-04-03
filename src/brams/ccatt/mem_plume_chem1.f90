!###########################################################################
!  B - Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_plume_chem1

  USE ModNamelistFile, only: namelistFile

  include "i8.h"

  type plume_vars   
     real, pointer, dimension(:,:)  :: fct     
  end type plume_vars  
  type (plume_vars)    , allocatable :: plume_g (:,:,:), plumem_g (:,:,:)

  type plume_mean_vars   
     real, pointer, dimension(:,:)  :: fire_size     
     real, pointer, dimension(:,:)  :: flam_frac     
  end type plume_mean_vars  

  type (plume_mean_vars)    , allocatable :: plume_mean_g(:,:), plume_meanm_g(:,:)


  integer, parameter :: nveg_agreg      = 4
  integer, parameter :: tropical_forest = 1
  integer, parameter :: boreal_forest   = 2
  integer, parameter :: savannah        = 3
  integer, parameter :: grassland       = 4 ! must be equal to nveg_agreg
  character(len=20), parameter :: veg_name(nveg_agreg) = (/ &
                               'Tropical-Forest', &
                               'Boreal-Forest  ', &
                               'Savanna        ', &
                               'Grassland      ' /)
  character(len=20), parameter :: spc_suf(nveg_agreg) = (/ &
                               'agtf' , &  ! trop forest
                               'agef' , &  ! extratrop forest
                               'agsv' , &  ! savanna
                               'aggr'   /) ! grassland
  REAL                         :: prfrq 	    
  INTEGER                      :: plumerise 


 !-- this section is for the FRP methodology

  type plume_fre_vars   
     real, pointer, dimension(:,:)  :: pvar     
  end type plume_fre_vars  
  type (plume_fre_vars)    , allocatable :: plume_fre_g (:,:), plumem_fre_g (:,:)
 
  integer, parameter ::      &
           iflam_frac  =1  &
          ,imean_frp   =2  &
          ,istd_frp    =3  &
          ,imean_size  =4  &
          ,istd_size   =5  

  character(len=10), parameter :: fre_var_name(5) = (/ &
                               'flam_frac ' , &  ! 
                               'mean_frp  ' , &  ! 
                               'std_frp   ' , &  ! 
                               'mean_size ' , &  ! 
                               'std_size  '   /) ! 




contains
  !---------------------------------------------------------------

  subroutine alloc_plume_chem1(plume,plume_mean,plume_fre,n1,n2,n3,nspecies )
                              
    use chem1_list, only : spc_alloc,spc_name,src,on
    implicit none

    integer,intent(in) :: n1,n2,n3,nspecies
    integer :: ispc,iv
    
    type (plume_vars)      ,dimension(nveg_agreg,nspecies) :: plume
    type (plume_mean_vars) ,dimension(nveg_agreg         ) :: plume_mean

    type (plume_fre_vars)  ,dimension(5                  ) :: plume_fre

    integer:: imean_plume
    imean_plume = 1 !change this at emis_flam_smold routine also
         
    !print*,'----------------------------------------------------------------'
    !print*,' memory allocation for plumerise sources:'
    IF(PLUMERISE == 1) THEN
      if(imean_plume /= on) then
    
    	  do ispc=1,nspecies
    	   !print*,'spc=',spc_name(ispc),'size=',n1,n2,n3

    	   if(spc_alloc(src,ispc) == on) then 

    	      do iv=1,nveg_agreg
    		  allocate (plume(iv,ispc)%fct(n2,n3))
    		  plume(iv,ispc)%fct(:,:) = 0.
    	      enddo
    	   endif
          enddo
       else 
 
          do iv=1,nveg_agreg
            allocate (plume_mean(iv)%flam_frac(n2,n3))
            plume_mean(iv)%flam_frac(:,:)=0.
          enddo
      endif


      do iv=1,nveg_agreg
        allocate (plume_mean(iv)%fire_size(n2,n3))
        plume_mean(iv)%fire_size(:,:)=0.
      enddo
    ENDIF
    !- for FRP methodology
    IF(PLUMERISE == 2)THEN
      do iv=1,5
         allocate (plume_fre(iv)%pvar(n2,n3))
         plume_fre(iv)%pvar(:,:)=0.
      enddo
    
    ENDIF

    return
  end subroutine alloc_plume_chem1

  !---------------------------------------------------------------

  subroutine dealloc_plume_chem1(plume,plume_mean,plume_fre,nspecies)

   implicit none

   integer,intent(in) :: nspecies
   type (plume_vars)    ,dimension(nveg_agreg,nspecies) :: plume
   type (plume_mean_vars) ,dimension(nveg_agreg         ) :: plume_mean
   type (plume_fre_vars)  ,dimension(5                  ) :: plume_fre
   integer :: ispc,iv
  
    !  Deallocate arrays
    do ispc=1,nspecies
     do iv=1,nveg_agreg
       if (associated(plume(iv,ispc)%fct)) deallocate(plume(iv,ispc)%fct)
     enddo
    enddo
    do iv=1,nveg_agreg
       if (associated(plume_mean(iv)%fire_size)) deallocate(plume_mean(iv)%fire_size)
       if (associated(plume_mean(iv)%flam_frac)) deallocate(plume_mean(iv)%flam_frac)
    enddo
     
    do iv=1,5
       if (associated(plume_fre(iv)%pvar)) deallocate(plume_fre(iv)%pvar)
    enddo
    
     
 end subroutine dealloc_plume_chem1

  !---------------------------------------------------------------

  subroutine nullify_plume_chem1(plume,plume_mean,plume_fre,nspecies)

    implicit none

    integer,intent(in) :: nspecies
    type (plume_vars)    ,dimension(nveg_agreg,nspecies) :: plume
    type (plume_mean_vars) ,dimension(nveg_agreg         ) :: plume_mean
    type (plume_fre_vars)  ,dimension(5                  ) :: plume_fre
    
    integer :: ispc,iv
 
    do ispc=1,nspecies

      do iv=1,nveg_agreg
        if (associated(plume(iv,ispc)%fct)) nullify(plume(iv,ispc)%fct)
      enddo
    
    enddo
    do iv=1,nveg_agreg
       if (associated(plume_mean(iv)%fire_size)) nullify(plume_mean(iv)%fire_size)
       if (associated(plume_mean(iv)%flam_frac)) nullify(plume_mean(iv)%flam_frac)
    enddo
 
    !- this is for FRP method
    do iv=1,5
       if (associated(plume_fre(iv)%pvar)) nullify(plume_fre(iv)%pvar)
    enddo

    return
  end subroutine nullify_plume_chem1

  !---------------------------------------------------------------

  subroutine filltab_plume_chem1(plume,plumem,plume_mean,plume_meanm,    &
                                 plume_fre,plumem_fre,imean,n1,n2,n3,nspecies,ng)

    use chem1_list, only: spc_alloc,spc_name,src,on 
    use mem_chem1, only: chem1_g

    use var_tables, only: InsertVTab

    implicit none

    integer, intent(in) :: imean,n1,n2,n3,nspecies,ng
    type (plume_vars)     ,dimension(nveg_agreg,nspecies) :: plume,plumem
    type (plume_mean_vars),dimension(nveg_agreg) ::plume_mean,plume_meanm
    type (plume_fre_vars) ,dimension(5)          ::plume_fre,plumem_fre

    integer :: ispc,iv  

    integer(kind=i8) :: npts

    integer:: imean_plume
    imean_plume = 1 !change this at emis_flam_smold routine also
         

    ! Fill pointers to arrays into variable tables
    ! 2d var
    npts=n2*n3

    IF(PLUMERISE == 1) THEN
      if(imean_plume /= on) then

    	do ispc=1,nspecies
    	  if (associated(chem1_g(ispc,ng)%sc_p)) then  ! check this latter

          !------- sources 
    	   
    	    if(spc_alloc(src,ispc) == on) then
    		do iv=1,nveg_agreg

    		    call InsertVTab(plume(iv,ispc)%fct,       &
                                    plumem(iv,ispc)%fct,      &
    				    ng, npts, imean,          &
                                    trim(spc_name(ispc))//'_'//trim(spc_suf(iv))//&
                                    ' :2:hist:anal:mpti:mpt3:mpt1')

    		enddo
    	    endif
    	   
    	  endif

    	enddo
      !---  flam frac mean
      else
        do iv=1,nveg_agreg

           call InsertVTab(plume_mean(iv)%flam_frac,  &
             	           plume_meanm(iv)%flam_frac,  &
                           ng, npts, imean, 'flam_frac'&
                           //'_'//trim(spc_suf(iv))//&
                           ' :2:hist:anal:mpti:mpt3:mpt1')

        enddo
      endif

      !--- fire size   
      do iv=1,nveg_agreg

       call InsertVTab(plume_mean(iv)%fire_size,  &
                       plume_meanm(iv)%fire_size, &
                       ng, npts, imean, 'firesize'&
                       //'_'//trim(spc_suf(iv))//&
                       ' :2:hist:anal:mpti:mpt3:mpt1')
   
      enddo
    ENDIF
    IF(PLUMERISE == 2)THEN
       do iv=1,5

        call InsertVTab(plume_fre(iv)%pvar,plumem_fre(iv)%pvar, &
                        ng, npts, imean, fre_var_name(iv)//&
                        ' :2:hist:anal:mpti:mpt3:mpt1')
   
       enddo
    ENDIF
      
  end subroutine filltab_plume_chem1

  !---------------------------------------------------------------

  subroutine StoreNamelistFileAtMem_plumeChem1(oneNamelistFile)
    type(namelistFile), pointer :: oneNamelistFile
    plumerise  = oneNamelistFile%plumerise
    prfrq = oneNamelistFile%prfrq
  end subroutine StoreNamelistFileAtMem_plumeChem1


  !---------------------------------------------------------------
end module mem_plume_chem1
