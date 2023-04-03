!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_mksfc


  type sfcfile_vars

     !(nxp,nyp,nzg,npatch)
     real, pointer :: soil_text(:,:,:,:)

     !(nxp,nyp,npatch)
     real, pointer :: patch_area(:,:,:)
     real, pointer :: leaf_class(:,:,:)
     real, pointer :: veg_ndvif(:,:,:)
     
     !(nxp,nyp)
     real, pointer :: topt(:,:)
     real, pointer :: seatf(:,:)
     real, pointer :: topzo(:,:)
     
     ! TEB_SPM
     real, pointer :: fuso(:,:)
     
  end type sfcfile_vars
   

  type (sfcfile_vars), allocatable :: sfcfile_p(:)
  
  !(nxpmax,nypmax)
  real, allocatable :: mksfc_scr1(:,:)
  real, allocatable :: mksfc_scr2(:,:)
  real, allocatable :: mksfc_vt2da(:,:)
  real, allocatable :: mksfc_vt2db(:,:)
  
  !(max(nxp*nyp)*nzg)
  real, allocatable :: scrx(:)
  
  !(np,np,nxp,nyp)
  real, allocatable :: glatp(:,:,:,:)
  real, allocatable :: glonp(:,:,:,:)
  real, allocatable :: datp(:,:,:,:)

  !(np,np,nxp,nyp)
  integer, allocatable :: datq_patch(:,:,:,:)
  
  !(np*np*nxp*nyp)
  integer, allocatable :: ptable(:)
  
  !(iblksizo,iblksizo)
  real, allocatable :: dato(:,:)
  
  !(iblksizo,iblksizo)
  character(len=1), allocatable :: cdato(:,:)
  
  !(ifile_max,jfile_max)
  integer, allocatable :: nump(:,:)
  integer, allocatable :: numpind(:,:)
  integer, allocatable :: numpind1(:,:)
  integer, allocatable :: numpind2(:,:)

  
  integer :: npq
  
  integer, parameter :: maxsfcgrids=10
  
  include "files.h"

  ! SST file creation variables
  integer, parameter :: maxsstdata=100
  integer :: iyearvs(maxsstdata,maxsfcgrids)
  integer :: imonthvs(maxsstdata,maxsfcgrids)
  integer :: idatevs(maxsstdata,maxsfcgrids)
  integer :: ihourvs(maxsstdata,maxsfcgrids)

  integer         :: nvsstf(maxsfcgrids)
  character(len=f_name_length)     :: vsstfil(maxsstdata,maxsfcgrids)
  
  ! NDVI file creation variables
  integer, parameter :: maxndvidata=100
  integer :: iyearvn(maxndvidata,maxsfcgrids)
  integer :: imonthvn(maxndvidata,maxsfcgrids)
  integer :: idatevn(maxndvidata,maxsfcgrids)
  integer :: ihourvn(maxndvidata,maxsfcgrids)

  integer         :: nvndvif(maxsfcgrids)
  character(len=f_name_length)     :: vndvifil(maxndvidata,maxsfcgrids)
  
contains







  subroutine alloc_sfcfile(sfcfile,nx,ny,nzg,npat)
    use teb_spm_start, only: TEB_SPM ! INTENT(IN)

    implicit none

    type (sfcfile_vars), intent(inout) :: sfcfile
    integer, intent(in) :: nx,ny,nzg,npat
    integer :: ierr
    character(len=8) :: c0, c1, c2, c3, c4
    character(len=*), parameter :: h="**(alloc_sfcfile)**"
     
    allocate (sfcfile%soil_text(nzg,nx,ny,npat), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       write(c1,"(i8)") nzg
       write(c2,"(i8)") nx
       write(c3,"(i8)") ny
       write(c4,"(i8)") npat
       call fatal_error(h//" allocate soil_text("//&
            trim(adjustl(c1))//","//trim(adjustl(c2))//","//&
            trim(adjustl(c3))//","//trim(adjustl(c4))//&
            ") fails with stat="//trim(adjustl(c0)))
    end if
    
    allocate (sfcfile%patch_area(nx,ny,npat), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       write(c2,"(i8)") nx
       write(c3,"(i8)") ny
       write(c4,"(i8)") npat
       call fatal_error(h//" allocate patch_area("//&
            trim(adjustl(c2))//","//&
            trim(adjustl(c3))//","//trim(adjustl(c4))//&
            ") fails with stat="//trim(adjustl(c0)))
    end if

    allocate (sfcfile%leaf_class(nx,ny,npat), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       write(c2,"(i8)") nx
       write(c3,"(i8)") ny
       write(c4,"(i8)") npat
       call fatal_error(h//" allocate leaf_class("//&
            trim(adjustl(c2))//","//&
            trim(adjustl(c3))//","//trim(adjustl(c4))//&
            ") fails with stat="//trim(adjustl(c0)))
    end if

    allocate (sfcfile%veg_ndvif(nx,ny,npat), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       write(c2,"(i8)") nx
       write(c3,"(i8)") ny
       write(c4,"(i8)") npat
       call fatal_error(h//" allocate veg_ndvif("//&
            trim(adjustl(c2))//","//&
            trim(adjustl(c3))//","//trim(adjustl(c4))//&
            ") fails with stat="//trim(adjustl(c0)))
    end if
    
    allocate (sfcfile%topt(nx,ny), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       write(c2,"(i8)") nx
       write(c3,"(i8)") ny
       call fatal_error(h//" allocate topt("//trim(adjustl(c2))//","//&
            trim(adjustl(c3))//") fails with stat="//trim(adjustl(c0)))
    end if
    
    allocate (sfcfile%seatf(nx,ny), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       write(c2,"(i8)") nx
       write(c3,"(i8)") ny
       call fatal_error(h//" allocate seatf("//trim(adjustl(c2))//","//&
            trim(adjustl(c3))//") fails with stat="//trim(adjustl(c0)))
    end if
    
    allocate (sfcfile%topzo(nx,ny), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       write(c2,"(i8)") nx
       write(c3,"(i8)") ny
       call fatal_error(h//" allocate topzo("//trim(adjustl(c2))//","//&
            trim(adjustl(c3))//") fails with stat="//trim(adjustl(c0)))
    end if
    
    if (TEB_SPM==1) then
       allocate (sfcfile%fuso(nx,ny), stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") ierr
          write(c2,"(i8)") nx
          write(c3,"(i8)") ny
          call fatal_error(h//" allocate fuso("//trim(adjustl(c2))//","//&
               trim(adjustl(c3))//") fails with stat="//trim(adjustl(c0)))
       end if
    end if
  end subroutine alloc_sfcfile
  
  ! ******************************************************************
  
  subroutine dealloc_sfcfile(sfcfile)
    
    use teb_spm_start, only: TEB_SPM ! INTENT(IN)
    
    implicit none
    
    type (sfcfile_vars), intent(inout) :: sfcfile
    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(dealloc_sfcfile)**"
    
    deallocate (sfcfile%soil_text, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" deallocate soil_text fails with stat="//trim(adjustl(c0)))
    end if
    
    deallocate (sfcfile%patch_area, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" deallocate patch_area fails with stat="//trim(adjustl(c0)))
    end if

    deallocate (sfcfile%leaf_class, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" deallocate leaf_class fails with stat="//trim(adjustl(c0)))
    end if

    deallocate (sfcfile%veg_ndvif, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" deallocate veg_ndvif fails with stat="//trim(adjustl(c0)))
    end if
    
    deallocate (sfcfile%topt, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" deallocate topt fails with stat="//trim(adjustl(c0)))
    end if

    deallocate (sfcfile%seatf, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" deallocate seatf fails with stat="//trim(adjustl(c0)))
    end if

    deallocate (sfcfile%topzo, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" deallocate topzo fails with stat="//trim(adjustl(c0)))
    end if

    if (TEB_SPM==1) then
       deallocate (sfcfile%fuso, stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") ierr
          call fatal_error(h//" deallocate fuso fails with stat="//trim(adjustl(c0)))
       end if
    endif
  end subroutine dealloc_sfcfile
end module mem_mksfc
