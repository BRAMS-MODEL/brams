!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
!:DOC%BEGIN
!:DOC%TITLE Modulo Mem_Radiate

Module mem_radiate

  use ModNamelistFile, only: namelistFile

   Type radiate_vars
   
      ! Variables to be dimensioned by (nzp,nxp,nyp)
   real, pointer, dimension(:,:,:) :: &
                          fthrd, &
!--(DMK-CCATT-INI)-----------------------------------------------------------
                          cldamnt,cluamnt,fthrd_sw,fthrd_lw  !NER
!--(DMK-CCATT-FIM)-----------------------------------------------------------
			  
                          
      ! Variables to be dimensioned by (nxp,nyp)
   real, pointer, dimension(:,:) :: &
                          rshort,rlong,rlongup,albedt,cosz, &
!--(DMK-CCATT-INI)-----------------------------------------------------------
                          rshortdif,sw_up_toa,lw_up_toa !NER
!--(DMK-CCATT-FIM)-----------------------------------------------------------
			  
    real, pointer, dimension(:,:,:) :: &
                          cloud_fraction
   End Type
   
   type (radiate_vars), allocatable :: radiate_g(:), radiatem_g(:)
   
   integer :: lonrad ! from RAMSIN
   integer :: ilwrtyp ! from RAMSIN
   integer :: iswrtyp ! from RAMSIN
   real    :: radfrq ! from RAMSIN
   real    :: radtun ! from RAMSIN
   integer :: ncall_i !Indica primeira chamada
   real    :: prsnz,prsnzp !Calculadas na primeira chamada
     
Contains

   subroutine alloc_radiate(radiate,n1,n2,n3,ng)

   implicit none
   type (radiate_vars) :: radiate
   integer, intent(in) :: n1,n2,n3,ng

! Allocate arrays based on options (if necessary)
      
      if( ilwrtyp+iswrtyp > 0)  then
                         allocate (radiate%fthrd(n1,n2,n3))
                         allocate (radiate%cloud_fraction(n1,n2,n3))
                         allocate (radiate%rshort(n2,n3))
                         allocate (radiate%rlong(n2,n3))
                         allocate (radiate%rlongup(n2,n3))
                         allocate (radiate%albedt(n2,n3))
                         allocate (radiate%cosz(n2,n3))


!-20/10/2015: srf - not being used, actually
!
!--(DMK-CCATT-INI)-----------------------------------------------------------
!NER
!			 allocate (radiate%sw_up_toa(n2,n3))
!			 allocate (radiate%lw_up_toa(n2,n3))
!			 allocate (radiate%rshortdif(n2,n3))
!			 allocate (radiate%cldamnt(n1,n2,n3))
!			 allocate (radiate%cluamnt(n1,n2,n3))  
!			 allocate (radiate%fthrd_sw(n1,n2,n3)) 
!			 allocate (radiate%fthrd_lw(n1,n2,n3)) 
!--(DMK-CCATT-FIM)-----------------------------------------------------------

!--(DMK-LFR NEC-SX6)----------------------------------------------
                         radiate%fthrd   = 0.
			 radiate%rshort  = 0.
			 radiate%rlong   = 0.
			 radiate%rlongup = 0.
			 radiate%albedt  = 0.
			 radiate%cosz    = 0.
!-20/10/2015: srf - not being used, actually
!--(DMK-CCATT-INI)-----------------------------------------------------------
!			 !NER
!			 radiate%cldamnt  = 0.
!			 radiate%cluamnt  = 0.
!			 radiate%fthrd_sw  = 0.
!			 radiate%fthrd_sw  = 0.
!			 radiate%sw_up_toa  = 0.
!			 radiate%lw_up_toa  = 0.
!			 radiate%rshortdif  = 0. 
!--(DMK-CCATT-FIM)-----------------------------------------------------------			 
!--(DMK-LFR NEC-SX6)----------------------------------------------
			 
      endif
                         
   return
   end subroutine


   subroutine nullify_radiate(radiate)

   implicit none
   type (radiate_vars) :: radiate
   

   if (associated(radiate%fthrd))    nullify (radiate%fthrd)
   if (associated(radiate%cloud_fraction)) nullify (radiate%cloud_fraction)
   if (associated(radiate%rshort))   nullify (radiate%rshort)
   if (associated(radiate%rlong))    nullify (radiate%rlong)
   if (associated(radiate%rlongup))  nullify (radiate%rlongup)
   if (associated(radiate%albedt))   nullify (radiate%albedt)
   if (associated(radiate%cosz))     nullify (radiate%cosz)

!--(DMK-CCATT-INI)-----------------------------------------------------------
   !NER
   if (associated(radiate%cldamnt))  nullify (radiate%cldamnt)
   if (associated(radiate%cluamnt))  nullify (radiate%cluamnt)
   if (associated(radiate%rshortdif))nullify (radiate%rshortdif)
   if (associated(radiate%sw_up_toa))nullify (radiate%sw_up_toa)
   if (associated(radiate%lw_up_toa))nullify (radiate%lw_up_toa)
   if (associated(radiate%fthrd_sw)) nullify (radiate%fthrd_sw) 
   if (associated(radiate%fthrd_lw)) nullify (radiate%fthrd_lw)
!--(DMK-CCATT-FIM)-----------------------------------------------------------

   return
   end subroutine

   subroutine dealloc_radiate(radiate)

   implicit none
   type (radiate_vars) :: radiate
   

   if (associated(radiate%fthrd))    deallocate (radiate%fthrd)
   if (associated(radiate%cloud_fraction)) deallocate (radiate%cloud_fraction)
   if (associated(radiate%rshort))   deallocate (radiate%rshort)
   if (associated(radiate%rlong))    deallocate (radiate%rlong)
   if (associated(radiate%rlongup))  deallocate (radiate%rlongup)
   if (associated(radiate%albedt))   deallocate (radiate%albedt)
   if (associated(radiate%cosz))     deallocate (radiate%cosz)

!--(DMK-CCATT-INI)-----------------------------------------------------------
   !NER
   if (associated(radiate%cldamnt))  deallocate (radiate%cldamnt)
   if (associated(radiate%cluamnt))  deallocate (radiate%cluamnt) 
   if (associated(radiate%rshortdif))deallocate (radiate%rshortdif)
   if (associated(radiate%sw_up_toa))deallocate (radiate%sw_up_toa)
   if (associated(radiate%lw_up_toa))deallocate (radiate%lw_up_toa)
   if (associated(radiate%fthrd_sw)) deallocate (radiate%fthrd_sw)
   if (associated(radiate%fthrd_lw)) deallocate (radiate%fthrd_lw)
!--(DMK-CCATT-FIM)-----------------------------------------------------------

   return
   end subroutine


subroutine filltab_radiate(radiate,radiatem,imean,n1,n2,n3,ng)
  use var_tables, only: InsertVTab
   implicit none
   include "i8.h"
   type (radiate_vars) :: radiate,radiatem
   integer, intent(in) :: imean,n1,n2,n3,ng
   integer(kind=i8) :: npts
   real, pointer :: var,varm

! Fill pointers to arrays into variable tables

   npts=n1*n2*n3

   if (associated(radiate%cloud_fraction))  &
      call InsertVTab (radiate%cloud_fraction,radiatem%cloud_fraction  &
                 ,ng, npts, imean,  &
                 'CLOUD_FRACTION :3:anal:mpti:mpt3')

   if (associated(radiate%fthrd))  &
      call InsertVTab (radiate%fthrd,radiatem%fthrd  &
                 ,ng, npts, imean,  &
                 'FTHRD :3:hist:anal:mpti:mpt3')

!--(DMK-CCATT-INI)-----------------------------------------------------------
   if (associated(radiate%cldamnt))  &                     !NER
      call InsertVTab (radiate%cldamnt,radiatem%cldamnt  &
                 ,ng, npts, imean,  &
                 'CLDAMNT :3:hist:anal:mpti:mpt3')
		 
   if (associated(radiate%fthrd_sw))  &                     !NER
      call InsertVTab (radiate%fthrd_sw,radiatem%fthrd_sw  &
                 ,ng, npts, imean,  &
                 'FTHRD_SW :3:hist:anal:mpti:mpt3')
   
   if (associated(radiate%fthrd_lw))  &                     !NER
      call InsertVTab (radiate%fthrd_lw,radiatem%fthrd_lw  &
                 ,ng, npts, imean,  &
                 'FTHRD_LW :3:hist:anal:mpti:mpt3')
  
  if (associated(radiate%cluamnt))  &                     !NER
      call InsertVTab (radiate%cluamnt,radiatem%cluamnt  &
                 ,ng, npts, imean,  &
                 'CLUAMNT :3:hist:anal:mpti:mpt3')			 
!--(DMK-CCATT-FIM)-----------------------------------------------------------

   npts=n2*n3
   if (associated(radiate%rshort))  &
      call InsertVTab (radiate%rshort,radiatem%rshort  &
                 ,ng, npts, imean,  &
                 'RSHORT :2:hist:anal:mpti:mpt3')
   if (associated(radiate%rlong))  &
      call InsertVTab (radiate%rlong,radiatem%rlong  &
                 ,ng, npts, imean,  &
                 'RLONG :2:hist:anal:mpti:mpt3')
   if (associated(radiate%rlongup))  &
      call InsertVTab (radiate%rlongup,radiatem%rlongup  &
                 ,ng, npts, imean,  &
                 'RLONGUP :2:hist:anal:mpti:mpt3')
   if (associated(radiate%albedt))  &
      call InsertVTab (radiate%albedt,radiatem%albedt  &
                 ,ng, npts, imean,  &
                 'ALBEDT :2:hist:anal:mpti:mpt3')
   if (associated(radiate%cosz))  &
      call InsertVTab (radiate%cosz,radiatem%cosz  &
                 ,ng, npts, imean,  &
                 'COSZ :2:hist:anal:mpt3')

!--(DMK-CCATT-INI)-----------------------------------------------------------
  if (associated(radiate%rshortdif))  &
      call InsertVTab (radiate%rshortdif,radiatem%rshortdif  &
                 ,ng, npts, imean,  &
                 'RSHORTDIF :2:hist:anal:mpt3')	
		 	 
  if (associated(radiate%sw_up_toa))  &
      call InsertVTab (radiate%sw_up_toa,radiatem%sw_up_toa  &
                 ,ng, npts, imean,  &
                 'SW_UP_TOA :2:hist:anal:mpt3')	
		 
  if (associated(radiate%lw_up_toa))  &
      call InsertVTab (radiate%lw_up_toa,radiatem%lw_up_toa  &
                 ,ng, npts, imean,  &
                 'LW_UP_TOA :2:hist:anal:mpt3')	
!--(DMK-CCATT-FIM)-----------------------------------------------------------

   return
   end subroutine


  subroutine StoreNamelistFileAtMem_radiate(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    ilwrtyp = oneNamelistFile%ilwrtyp
    iswrtyp = oneNamelistFile%iswrtyp
    lonrad = oneNamelistFile%lonrad
    radfrq = oneNamelistFile%radfrq
    radtun = oneNamelistFile%radtun
  end subroutine StoreNamelistFileAtMem_radiate
End Module mem_radiate
