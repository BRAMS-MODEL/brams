!###########################################################################
! CCATT- B - Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_chem1

  use grid_dims, only : maxgrds
  use chem1_list, only : maxnspecies

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  use ModNamelistFile, only: namelistFile

  include "i8.h"
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

  TYPE chem1_vars   
     !--- All family
     REAL, POINTER, DIMENSION(:,:,:)  :: sc_p,sc_pp,sc_pf
     REAL, POINTER, DIMENSION(:,:  )  :: sc_dd,sc_wd
     REAL, POINTER, DIMENSION(:    )  :: sc_t,sc_t_dyn
     !-----------
  END TYPE chem1_vars
  type (chem1_vars)    , allocatable :: chem1_g(:,:), chem1m_g(:,:)

  INTEGER, PARAMETER :: maxsrcfiles   = 1500 

  integer, parameter :: nsrc=4  !number_sources
  
  type chem1_src_vars   
     real, pointer, dimension(:,:,:)  :: sc_src
  end type chem1_src_vars

!!$  type (chem1_src_vars), allocatable :: chem1_src_g(:,:,:),chem1m_src_g(:,:,:)
  type (chem1_src_vars), allocatable :: chem1_src_g(:,:,:,:),chem1m_src_g(:,:,:,:)

  !- dimension of sources arrays (=1 for 2dim, =m1 for 3dim)
  integer :: chem1_src_z_dim_g(nsrc,maxgrds)

  character(LEN=20),dimension(nsrc),parameter :: src_name= &
! '12345678901234567890'
(/                      &
  'antro               '&    
, 'bburn               '&
, 'bioge               '&
, 'geoge               '/)

  integer, parameter     :: &
   antro   = 01, & ! anthropogenic sources
   bburn   = 02, & ! biomass burning sources 
   bioge   = 03, & ! biogenic sources 
   geoge   = 04    ! geogenic/volc sources ! must be equal to "nsrc"
  
  !- use of the prescribed diurnal cycle or linear interpolation for the instantaneous
  !- emission rate
  !- for biomass burning, the diur_cycle must be always 1 (the 2nd element of diur_cycle array)
  !- 1=on, 0=off (=> will use linterp.)
  integer, dimension(nsrc) :: diur_cycle !diur_cycle(1)== antro; diur_cycle(2) == bburn
                                         !diur_cycle(3)== bioge; diur_cycle(4) == geoge
						    
  integer, parameter :: max_ntimes_src = 2  !- number maximum of src files for linterp.
  integer, dimension(nsrc) :: ntimes_src    !- actual number used

  integer :: RECYCLE_TRACERS, NSPECIES_TRANSPORTED,NSPECIES_CHEM_TRANSPORTED &
            ,NSPECIES_CHEM_NO_TRANSPORTED
  integer, dimension(maxnspecies) :: TRANSP_CHEM_INDEX,NO_TRANSP_CHEM_INDEX
  integer :: CHEM_ASSIM & ! determine if 4DDA will or  not be used
            ,CHEMISTRY  & ! define if the chemistry will run and in this case
	        	  !  the type of solver 
            ,ISPLIT	  ! to control when chemistry will be called
	    
  character(LEN=20) :: SPLIT_METHOD  ! determine the type of splliting method
  real CHEM_TIMESTEP                 ! timestep of chemistry integration
  integer :: N_DYN_CHEM_N(maxgrds)&  ! number of dynamic timesteps per chem. timestep
            ,N_DYN_CHEM              ! current N_DYN_CHEM_N

contains
  !---------------------------------------------------------------
  subroutine  define_n_dyn_chem(ngrids,dtlong,nndtrat,mynum)
  !  use mem_grid , only : nndtrat
    implicit none
    integer, intent(in) :: ngrids,nndtrat(ngrids),mynum
    real, intent(in) :: dtlong
    integer ng
    

    !- get the number of dynamics cycles inside each chemistry cycle:
    !- observe that 'dtlong' (timestep of coarser grid) is used, instead
    !- of dtlongn(ngrid) , the grid dependent long timestep.
    
    !- for the coarser grid
    N_DYN_CHEM_N(1) =max(1,nint(CHEM_TIMESTEP/DTLONG)) 

    !- for nested grids
    do ng =2 , ngrids
       N_DYN_CHEM_N(ng)=min(nndtrat(ng),  N_DYN_CHEM_N(1))  
    end do
    
    !- define split control to setup when chem. will be called   
    isplit = 1
    if(trim(adjustl(SPLIT_METHOD))=='SYMMETRIC') isplit = 2

    if((mynum == 0 .or. mynum ==1) .and. chemistry > 0) then
     print*,'----------------------------------------------------------------'
     print*,' -- > chemistry splitting: N_DYN_CHEM=',N_DYN_CHEM_N(1:NGRIDS)
     print*,'----------------------------------------------------------------'
    endif
    
  end subroutine  define_n_dyn_chem
  !---------------------------------------------------------------

  subroutine alloc_chem1(chem1,chem1_src,nvert_src,n1,n2,n3,nspecies,ng,volcanoes)

    use chem1_list, only : spc_alloc,spc_name, src, ddp, wdp, fdda, on ,off,&
                            transport
    
    implicit none

    integer,intent(in) :: n1,n2,n3,nspecies,ng,volcanoes
    integer :: ispc,isrc,itime
    
    type (chem1_vars)    ,dimension(     nspecies) :: chem1
    type (chem1_src_vars),dimension(max_ntimes_src,nsrc,nspecies) :: chem1_src
    integer,dimension(nsrc)    :: nvert_src
    
    !-determine the time dimension for src arrays
    do isrc=1,nsrc
      ! if   : ntimes_src(isrc) = 1, will use the prescribed diurnal cycle
      ! else : the array will be allocated with 2 dimensions to
      !        read 2 files and make linear interpolation
      ntimes_src(isrc)= 2 - diur_cycle(isrc)

      if(ntimes_src(isrc) .gt. max_ntimes_src  .or. &
         ntimes_src(isrc) .lt. 1  ) stop 'ntimes_src > max_ntimes_src or < 1'
    
      if(ntimes_src(bburn) > 1) stop "bburn must run with diurnal cycle ON"
      if(ntimes_src(geoge) > 1) stop "geoge must run with diurnal cycle ON"
    enddo
    
    !- current dynamic splliting
    N_DYN_CHEM=N_DYN_CHEM_N(ng)
    
    !- memory allocation for chemical species     
    do ispc=1,nspecies
     !print*,'spc=',spc_name(ispc),'size=',n1,n2,n3
     
     !- allocate memory for the past time tracer mixing ratio
     allocate (chem1(ispc)%sc_p  (n1,n2,n3)) ;chem1(ispc)%sc_p = 0.
     
     !- allocate memory for the tendency tracer mixing ratio if parallel spliting 
     !- operator will be used
     if(spc_alloc(transport,ispc)==on.and.trim(adjustl(SPLIT_METHOD))=='PARALLEL' &
          .and. N_DYN_CHEM > 1) then 
         allocate (chem1(ispc)%sc_t_dyn  (n1*n2*n3)) ;chem1(ispc)%sc_t_dyn = 0.
     else
         allocate (chem1(ispc)%sc_t_dyn  (1)       ) ;chem1(ispc)%sc_t_dyn = 0.
    endif
     
     !- allocate memory for the dry deposition of the tracer mixing ratio
     if(spc_alloc(ddp,ispc) == on) then 
           allocate (chem1(ispc)%sc_dd    (n2,n3)) ;chem1(ispc)%sc_dd= 0.
     endif

     !- allocate memory for the wet deposition of the tracer mixing ratio
     if(spc_alloc(wdp,ispc) == on) then
           allocate (chem1(ispc)%sc_wd    (n2,n3)) ;chem1(ispc)%sc_wd= 0.
     endif
     
     !- allocate memory for the past/future time tracer mixing ratio 
     !- (4D data assimilation)
     if(chem_assim == on) then
        if(spc_alloc(fdda,ispc) == on) then 
           allocate (chem1(ispc)%sc_pp (n1,n2,n3),chem1(ispc)%sc_pf (n1,n2,n3))
           chem1(ispc)%sc_pp =0.; chem1(ispc)%sc_pf = 0.
        endif
     endif

     !- allocate memory for the sources (3d=bburn and 2d= antro/bioge)
     if(spc_alloc(src,ispc) == on) then 
        do isrc=1,nsrc
     
	    do itime=1, ntimes_src(isrc) 

             !-for geoge, only alloc src if volcanoes is ON.
	     if(isrc==geoge .and. volcanoes == OFF) cycle  

             allocate (chem1_src(itime,isrc,ispc)%sc_src(nvert_src(isrc),n2,n3))
	     chem1_src(itime,isrc,ispc)%sc_src = 0.
	    
            enddo
        enddo
     endif

    enddo
    return
  end subroutine alloc_chem1

  !--------------------------------------------------------------------------

  subroutine dealloc_chem1(chem1,chem1_src,nspecies)

   implicit none

   integer,intent(in) :: nspecies
   type (chem1_vars)    ,dimension(     nspecies) :: chem1
   type (chem1_src_vars),dimension(max_ntimes_src,nsrc,nspecies) :: chem1_src
   integer :: ispc,isrc,itime
  
    !  Deallocate arrays
    do ispc=1,nspecies

     if (associated(chem1(ispc)%sc_p ))    deallocate(chem1(ispc)%sc_p  )
     if (associated(chem1(ispc)%sc_t_dyn ))    deallocate(chem1(ispc)%sc_t_dyn  )
     if (associated(chem1(ispc)%sc_dd))    deallocate(chem1(ispc)%sc_dd )
     if (associated(chem1(ispc)%sc_wd))    deallocate(chem1(ispc)%sc_wd )
     if (associated(chem1(ispc)%sc_pp))    deallocate(chem1(ispc)%sc_pp )
     if (associated(chem1(ispc)%sc_pf))    deallocate(chem1(ispc)%sc_pf )
     do isrc=1,nsrc
       do itime=1,max_ntimes_src
        if (associated(chem1_src(itime,isrc,ispc)%sc_src)) deallocate(chem1_src(itime,isrc,ispc)%sc_src)
      enddo
     enddo
    enddo

    return
  end subroutine dealloc_chem1

  !---------------------------------------------------------------
  !---------------------------------------------------------------

  subroutine nullify_chem1(chem1,chem1_src,nspecies)

    implicit none

    integer,intent(in) :: nspecies
    type (chem1_vars),dimension(nspecies) :: chem1
    type (chem1_src_vars),dimension(max_ntimes_src,nsrc,nspecies) :: chem1_src
    integer :: ispc,isrc,itime

    do ispc=1,nspecies

     if (associated(chem1(ispc)%sc_p ))    nullify (chem1(ispc)%sc_p  )
     if (associated(chem1(ispc)%sc_t_dyn ))    nullify (chem1(ispc)%sc_t_dyn  )
     if (associated(chem1(ispc)%sc_dd))    nullify (chem1(ispc)%sc_dd )
     if (associated(chem1(ispc)%sc_wd))    nullify (chem1(ispc)%sc_wd )
     if (associated(chem1(ispc)%sc_pp))    nullify (chem1(ispc)%sc_pp )
     if (associated(chem1(ispc)%sc_pf))    nullify (chem1(ispc)%sc_pf )
     do isrc=1,nsrc
       do itime=1,max_ntimes_src
         if (associated(chem1_src(itime,isrc,ispc)%sc_src)) nullify(chem1_src(itime,isrc,ispc)%sc_src)
       enddo
     enddo
    
    enddo

    return
  end subroutine nullify_chem1

  !---------------------------------------------------------------

 subroutine filltab_chem1(chem1,chem1m,chem1_src,chem1m_src&
                          ,imean,nvert_src,n1,n2,n3,nspecies,ng,volcanoes)

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
    use var_tables, only: InsertVTab
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

    use chem1_list, only: spc_alloc,spc_name, src, ddp, wdp, fdda, on ,off, transport

        use io_params, only : ioutput         ! INTENT(IN)

    implicit none

    integer, intent(in) :: imean,n1,n2,n3,nspecies,ng,volcanoes
    type (chem1_vars)    ,dimension(     nspecies) :: chem1,chem1m
    type (chem1_src_vars),dimension(max_ntimes_src,nsrc,nspecies) :: chem1_src,chem1m_src
    integer,dimension(nsrc)    :: nvert_src

    integer :: ispc,isrc,itime  

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
    integer(kind=i8) :: npts
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!   integer :: npts
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------


    character(len=8) :: str_recycle
    character(len=25) :: str_output
    character(len=1) :: str_src_dim,str_src_num
   
    str_recycle = ''; str_src_dim = '';str_src_num= ''
    if (RECYCLE_TRACERS == 1 .or. RECYCLE_TRACERS == 2 .or. ioutput ==5) then
       str_recycle = ':recycle'
    endif

    !- Fill pointers to arrays into variable tables
    do ispc=1,nspecies

     if (associated(chem1(ispc)%sc_p)) then !--- tracer mixing ratio (dimension 3d)
       npts = n1 * n2 * n3
       if(spc_alloc(transport,ispc) == on) then 

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
        call InsertVTab(chem1(ispc)%sc_p,chem1m(ispc)%sc_p,  &
                        ng, npts, imean,                     &
                        trim(spc_name(ispc))//'P :3:hist:anal:mpti:mpt3:mpt1'//trim(str_recycle))
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!        call vtables2 (chem1(ispc)%sc_p(1,1,1),chem1m(ispc)%sc_p(1,1,1)  &
!             ,ng, npts, imean, trim(spc_name(ispc))//'P :3:hist:anal:mpti:mpt3:mpt1'//trim(str_recycle))
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

       endif

!--- sources (3 and 2 dimension)
      !- old way
      ! npts=n1*n2*n3
      ! if(spc_alloc(1,ispc) == 1) &
      ! call vtables2 (chem1(ispc)%sc_s(1,1,1),chem1m(ispc)%sc_s(1,1,1)  &
      !      ,ng, npts, imean, trim(spc_name(ispc))//'S :3:hist:anal:mpti:mpt3:mpt1')
      !- new way
       if(spc_alloc(src,ispc) == on) then
            do isrc=1,nsrc
	
	      !- only alloc geoge emissions if volcanoes is ON.
	      if(isrc==geoge .and. volcanoes == off) cycle

	      npts= nvert_src(isrc)  * n2 * n3 * ntimes_src(isrc)
	      if(nvert_src(isrc) ==  1) str_src_dim = '2'  ! for 2d sources
	      if(nvert_src(isrc) == n1) str_src_dim = '3'  ! for 3d sources

	      do itime=1,ntimes_src(isrc)
	   
	        if(itime > 1) write(str_src_num,'(i1)')itime
!	        print *,'LFR: ',trim(spc_name(ispc))//'_'//trim(src_name(isrc))//  &
!			        '_SRC'//trim(str_src_num)//' :'//trim(str_src_dim) &
!			        //':hist:anal:mpti:mpt3:mpt1'
!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
	        call InsertVTab(chem1_src(itime,isrc,ispc)%sc_src,                 &
		                chem1m_src(itime,isrc,ispc)%sc_src,                &
                                ng, npts, imean,                                   &
                                trim(spc_name(ispc))//'_'//trim(src_name(isrc))//  &
			        '_SRC'//trim(str_src_num)//' :'//trim(str_src_dim) &
			        //':hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!	        call vtables2 ( chem1_src(itime,isrc,ispc)%sc_src(1,1,1)  &
!		              ,chem1m_src(itime,isrc,ispc)%sc_src(1,1,1)  &
!                              ,ng, npts, imean, trim(spc_name(ispc))      &
!			      //'_'//trim(src_name(isrc))//               &
!			      '_SRC'//trim(str_src_num)//' :'//trim(str_src_dim)&
!			      //':hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------
            
	      enddo
	    !print*,'src alloc=', trim(spc_name(ispc))//'_'//trim(src_name(isrc)),' npts=',npts
            enddo
       endif
       
!--- dry and wet deposition (dimension 2d)
       npts = n2 * n3
       if(spc_alloc(ddp,ispc) == on) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
       call InsertVTab(chem1(ispc)%sc_dd,chem1m(ispc)%sc_dd,   &
                       ng, npts, imean,                        &
                       trim(spc_name(ispc))//'DD :2:hist:anal:mpti:mpt3')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!       call vtables2 (chem1(ispc)%sc_dd(1,1),chem1m(ispc)%sc_dd (1,1)  &
!            ,ng, npts, imean, trim(spc_name(ispc))//'DD :2:hist:anal:mpti:mpt3')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

       npts = n2 * n3
       if(spc_alloc(wdp,ispc) == on) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
       call InsertVTab(chem1(ispc)%sc_wd,chem1m(ispc)%sc_wd,  &
                       ng, npts, imean,                       &
                       trim(spc_name(ispc))//'WD :2:hist:anal:mpti:mpt3')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!       call vtables2 (chem1(ispc)%sc_wd(1,1),chem1m(ispc)%sc_wd(1,1)  &
!           ,ng, npts, imean, trim(spc_name(ispc))//'WD :2:hist:anal:mpti:mpt3')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------


!---  data assimilation (dimension 3d)
       if(chem_assim == on) then
          npts = n1 * n2 * n3
          if(spc_alloc(fdda,ispc) == on) then

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
             call InsertVTab(chem1(ispc)%sc_pp,chem1m(ispc)%sc_pp,  &
                             ng, npts, imean,                       &
                             trim(spc_name(ispc))//'PP :3:mpti')
             call InsertVTab(chem1(ispc)%sc_pf,chem1m(ispc)%sc_pf,  &
                             ng, npts, imean,                       &
                             trim(spc_name(ispc))//'PF :3:mpti')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!             call vtables2 (chem1(ispc)%sc_pp(1,1,1),chem1m(ispc)%sc_pp(1,1,1)  &
!              ,ng, npts, imean, trim(spc_name(ispc))//'PP :3:mpti')
!             call vtables2 (chem1(ispc)%sc_pf(1,1,1),chem1m(ispc)%sc_pf(1,1,1)  &
!              ,ng, npts, imean, trim(spc_name(ispc))//'PF :3:mpti')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

          endif
        endif

     endif

    enddo
  end subroutine filltab_chem1

  !---------------------------------------------------------------

  subroutine alloc_tend_chem1(nmzp,nmxp,nmyp,ngrs,nspecies,proc_type)
   use chem1_list, only:spc_name,spc_alloc,transport,on
   implicit none
   integer,intent(in)                   :: ngrs,proc_type,nspecies
   integer,intent(in), dimension (ngrs) :: nmzp,nmxp,nmyp
   integer :: ng,ntpts,ispc

!         Find the maximum number of grid points needed for any grid.

   if(proc_type==1) then
      ntpts=1
   else
      ntpts=0
      do ng=1,ngrs
         ntpts=max( nmxp(ng)*nmyp(ng)*nmzp(ng),ntpts )
      enddo
   endif

   do ispc=1,nspecies

      !-srf: only tendencies arrays for transported species are allocated (save memory)
      if (spc_alloc(transport,ispc) == on) then  
          if (associated(chem1_g(ispc,1)%sc_p)) then
	      allocate  (chem1_g(ispc,1)%sc_t(ntpts))
	      chem1_g(ispc,1)%sc_t(1:ntpts)=0.
	  endif
	  !print*,'transp=',ispc,spc_name(ispc)
      else
      !- for non-transported species, arrays are allocated with one-dimension
          if (associated(chem1_g(ispc,1)%sc_p)) then 
	      allocate  (chem1_g(ispc,1)%sc_t(1    ))
	      chem1_g(ispc,1)%sc_t(1)= 0.
	  endif
	  !print*,'NO transp=',ispc,spc_name(ispc)
      endif
      
      do ng=2,ngrs
         chem1_g(ispc,ng)%sc_t => chem1_g(ispc,1)%sc_t
      enddo
      
    enddo

  end subroutine alloc_tend_chem1
  !---------------------------------------------------------------

  subroutine nullify_tend_chem1(nspecies)

    implicit none
    integer,intent(in) :: nspecies
    integer ::ispc

    do ispc=1,nspecies
      if (associated(chem1_g(ispc,1)%sc_t)) nullify (chem1_g(ispc,1)%sc_t)
    enddo


  end subroutine nullify_tend_chem1

  !---------------------------------------------------------------
  subroutine dealloc_tend_chem1(nspecies)
    implicit none
    integer,intent(in) :: nspecies
    integer ::ispc

    do ispc=1,nspecies
     if (associated(chem1_g(ispc,1)%sc_t)) deallocate (chem1_g(ispc,1)%sc_t)
    enddo

  end  subroutine dealloc_tend_chem1
   
  !---------------------------------------------------------------

  subroutine filltab_tend_chem1(nspecies,ng)
   use chem1_list, only:spc_name,spc_alloc,transport,on
   implicit none

    integer,intent(in) :: nspecies,ng
    integer ::ispc
    integer :: elements

    
    nspecies_chem_transported = 0
    transp_chem_index   (:)   = 0
    
    nspecies_chem_no_transported = 0
    no_transp_chem_index(:)      = 0
    
    do ispc=1,nspecies

    !- Fill pointers to scalar arrays into scalar tables 
    !-srf - only for the "transported" species


      if (spc_alloc(transport,ispc) == on .and. associated(chem1_g(ispc,ng)%sc_t)) then
      	call vtables_scalar (chem1_g(ispc,ng)%sc_p(1,1,1),chem1_g(ispc,ng)%sc_t(1),&
                           ng,trim(spc_name(ispc))//'P')
	elements = size(chem1_g(ispc,ng)%sc_t)
        call vtables_scalar_new (chem1_g(ispc,ng)%sc_p(1,1,1),chem1_g(ispc,ng)%sc_t(1),&
                           ng,trim(spc_name(ispc))//'P',elements)
	
	!- number of chem transported species
	nspecies_chem_transported = nspecies_chem_transported + 1	
        
	!- mapping between ispc and transported chem species
	transp_chem_index(nspecies_chem_transported) = ispc 	   

      else
        
	nspecies_chem_no_transported = nspecies_chem_no_transported + 1
        no_transp_chem_index(nspecies_chem_no_transported) = ispc
      
      
      endif		   

    enddo

    !- save the total number of chemical species to be transported in the
    !- total number of transported scalar (chemical + aerosols).
    !- this will be changed in mem_aer1.f90 routine to include aerosols.
    nspecies_transported = nspecies_chem_transported

  end subroutine filltab_tend_chem1
  

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  subroutine StoreNamelistFileAtMem_chem1(oneNamelistFile)
    type(namelistFile), pointer :: oneNamelistFile
    chemistry = oneNamelistFile%chemistry
    split_method = oneNamelistFile%split_method
    chem_timestep = oneNamelistFile%chem_timestep
    chem_assim = oneNamelistFile%chem_assim
    recycle_tracers = oneNamelistFile%recycle_tracers
    diur_cycle = oneNamelistFile%diur_cycle
  end subroutine StoreNamelistFileAtMem_chem1
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------


end module mem_chem1

!--------------------------------------------------------------------------

subroutine define_chem1_src_zdim(chem1_src_z_dim,n1)

  use mem_chem1, only: nsrc,bburn,antro,bioge,geoge
  implicit none
  integer,intent(in) :: n1
  
  integer,dimension(nsrc)    :: chem1_src_z_dim

   !- determination of the dimension of Z-dir of source field array
   chem1_src_z_dim(antro) = n1	     ! 2d 
   chem1_src_z_dim(bburn) = n1       ! 3d
   chem1_src_z_dim(bioge) = n1	     ! 2d
   chem1_src_z_dim(geoge) = n1       ! 3d for volcanoes
  return
end subroutine define_chem1_src_zdim

!--------------------------------------------------------------------------

