!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


Module mem_leaf

  use ModNamelistFile, only: namelistFile

  use grid_dims

  Type leaf_vars

     ! Variables to be dimensioned by (nxp,nyp,nzg,npatch)

     real, pointer :: soil_water(:,:,:,:)
     real, pointer :: soil_energy(:,:,:,:)
     real, pointer :: soil_text(:,:,:,:)

     ! Variables to be dimensioned by (nxp,nyp,nzs,npatch)

     real, pointer :: sfcwater_mass(:,:,:,:)
     real, pointer :: sfcwater_energy(:,:,:,:)
     real, pointer :: sfcwater_depth(:,:,:,:)

     ! Variables to be dimensioned by (nxp,nyp,npatch)

     real, pointer :: ustar(:,:,:)
     real, pointer :: tstar(:,:,:)
     real, pointer :: rstar(:,:,:)
     real, pointer :: veg_fracarea(:,:,:)
     real, pointer :: veg_lai(:,:,:)
     real, pointer :: veg_rough(:,:,:)
     real, pointer :: veg_height(:,:,:)
     real, pointer :: veg_albedo(:,:,:)
     real, pointer :: veg_tai(:,:,:)
     real, pointer :: patch_area(:,:,:)
     real, pointer :: patch_rough(:,:,:)
     real, pointer :: patch_wetind(:,:,:)
     real, pointer :: leaf_class(:,:,:)
     real, pointer :: soil_rough(:,:,:)
     real, pointer :: sfcwater_nlev(:,:,:)
     real, pointer :: stom_resist(:,:,:)
     real, pointer :: ground_rsat(:,:,:)
     real, pointer :: ground_rvap(:,:,:)
     real, pointer :: veg_water(:,:,:)
     real, pointer :: veg_temp(:,:,:)
     real, pointer :: can_rvap(:,:,:)
     real, pointer :: can_temp(:,:,:)
     real, pointer :: veg_ndvip(:,:,:)
     real, pointer :: veg_ndvic(:,:,:)
     real, pointer :: veg_ndvif(:,:,:)


     ! TEB_SPM
     real, pointer :: G_URBAN(:,:,:)

     real, pointer :: R_aer(:,:,:)   !kml drydep

     ! Variables to be dimensioned by (nxp,nyp)

     real, pointer :: snow_mass(:,:)
     real, pointer :: snow_depth(:,:)
     real, pointer :: seatp(:,:)
     real, pointer :: seatf(:,:)
  End Type leaf_vars

  type (leaf_vars), allocatable :: leaf_g(:), leafm_g(:)

  !----------------------------------------------------------------------------
  integer :: nslcon ! from RAMSIN
  integer :: nvgcon ! from RAMSIN
  integer :: nvegpat ! from RAMSIN
  integer :: isfcl ! from RAMSIN
  integer :: isfcl_ocean ! from RAMSIN
  real    :: zrough ! from RAMSIN
  real    :: pctlcon ! from RAMSIN
  real    :: ubmin
  real    :: albedo ! from RAMSIN
  real    :: drtcon ! from RAMSIN
  real    :: dthcon ! from RAMSIN
  real    :: seatmp ! from RAMSIN
  real    :: stgoff(nzgmax) ! from RAMSIN
  real    :: slmstr(nzgmax) ! from RAMSIN
  real    :: slz(nzgmax) ! from RAMSIN

Contains

  subroutine alloc_leaf(leaf,nz,nx,ny,nzg,nzs,np,ng)

    ! TEB_SPM
    use teb_spm_start, only: TEB_SPM ! INTENT(IN)

    implicit none
    type (leaf_vars) :: leaf
    integer, intent(in) :: nz,nx,ny,nzg,nzs,np,ng

    ! Allocate arrays based on options (if necessary)

    allocate (leaf%soil_water     (nzg,nx,ny,np));leaf%soil_water =0.0
    allocate (leaf%soil_energy    (nzg,nx,ny,np));leaf%soil_energy=0.0
    allocate (leaf%soil_text      (nzg,nx,ny,np));leaf%soil_text  =0.0

    allocate (leaf%sfcwater_mass  (nzs,nx,ny,np));leaf%sfcwater_mass  =0.0
    allocate (leaf%sfcwater_energy(nzs,nx,ny,np));leaf%sfcwater_energy=0.0
    allocate (leaf%sfcwater_depth (nzs,nx,ny,np));leaf%sfcwater_depth =0.0

    allocate (leaf%ustar        (nx,ny,np));leaf%ustar=0.0
    allocate (leaf%tstar        (nx,ny,np));leaf%tstar=0.0
    allocate (leaf%rstar        (nx,ny,np));leaf%rstar=0.0

    allocate (leaf%veg_fracarea (nx,ny,np));leaf%veg_fracarea=0.0
    allocate (leaf%veg_lai      (nx,ny,np));leaf%veg_lai     =0.0
    allocate (leaf%veg_rough    (nx,ny,np));leaf%veg_rough   =0.0
    allocate (leaf%veg_height   (nx,ny,np));leaf%veg_height  =0.0
    allocate (leaf%veg_albedo   (nx,ny,np));leaf%veg_albedo  =0.0
    allocate (leaf%veg_tai      (nx,ny,np));leaf%veg_tai     =0.0

    allocate (leaf%patch_area   (nx,ny,np));leaf%patch_area  =0.0
    allocate (leaf%patch_rough  (nx,ny,np));leaf%patch_rough =0.0
    allocate (leaf%patch_wetind (nx,ny,np));leaf%patch_wetind=0.0
    allocate (leaf%leaf_class   (nx,ny,np));leaf%leaf_class  =0.0

    ! TEB_SPM
    if (TEB_SPM==1) then
       allocate (leaf%G_URBAN   (nx,ny,np));leaf%G_URBAN=0.0
    endif

    allocate (leaf%soil_rough   (nx,ny,np));leaf%soil_rough   =0.0
    allocate (leaf%sfcwater_nlev(nx,ny,np));leaf%sfcwater_nlev=0.0
    allocate (leaf%stom_resist  (nx,ny,np));leaf%stom_resist  =0.0

    allocate (leaf%ground_rsat  (nx,ny,np));leaf%ground_rsat=0.0
    allocate (leaf%ground_rvap  (nx,ny,np));leaf%ground_rvap=0.0
;
    allocate (leaf%veg_water    (nx,ny,np));leaf%veg_water  =0.0
    allocate (leaf%veg_temp     (nx,ny,np));leaf%veg_temp   =0.0
;                                          
    allocate (leaf%can_rvap     (nx,ny,np));leaf%can_rvap   =0.0
    allocate (leaf%can_temp     (nx,ny,np));leaf%can_temp   =0.0
;                                          
    allocate (leaf%veg_ndvip    (nx,ny,np));leaf%veg_ndvip  =0.0
    allocate (leaf%veg_ndvic    (nx,ny,np));leaf%veg_ndvic  =0.0
    allocate (leaf%veg_ndvif    (nx,ny,np));leaf%veg_ndvif  =0.0

    allocate (leaf%R_aer        (nx,ny,np));leaf%R_aer =0.0  !kml drydep

    allocate (leaf%snow_mass    (nx,ny));leaf%snow_mass =0.0
    allocate (leaf%snow_depth   (nx,ny));leaf%snow_depth=0.0
    allocate (leaf%seatp        (nx,ny));leaf%seatp     =0.0
    allocate (leaf%seatf        (nx,ny));leaf%seatf     =0.0

  end subroutine alloc_leaf

  !************************************************************************

  subroutine nullify_leaf(leaf)

    ! TEB_SPM
    use teb_spm_start, only: TEB_SPM ! INTENT(IN)

    implicit none
    type (leaf_vars) :: leaf

    if(associated(leaf%soil_water))      nullify (leaf%soil_water)
    if(associated(leaf%soil_energy))     nullify (leaf%soil_energy)
    if(associated(leaf%soil_text))       nullify (leaf%soil_text)

    if(associated(leaf%sfcwater_mass))   nullify (leaf%sfcwater_mass)
    if(associated(leaf%sfcwater_energy)) nullify (leaf%sfcwater_energy)
    if(associated(leaf%sfcwater_depth))  nullify (leaf%sfcwater_depth)

    if(associated(leaf%ustar))           nullify (leaf%ustar)
    if(associated(leaf%tstar))           nullify (leaf%tstar)
    if(associated(leaf%rstar))           nullify (leaf%rstar)

    if(associated(leaf%veg_fracarea))    nullify (leaf%veg_fracarea)
    if(associated(leaf%veg_lai))         nullify (leaf%veg_lai)
    if(associated(leaf%veg_rough))       nullify (leaf%veg_rough)
    if(associated(leaf%veg_height))      nullify (leaf%veg_height)
    if(associated(leaf%veg_albedo))      nullify (leaf%veg_albedo)
    if(associated(leaf%veg_tai))         nullify (leaf%veg_tai)

    if(associated(leaf%patch_area))      nullify (leaf%patch_area)
    if(associated(leaf%patch_rough))     nullify (leaf%patch_rough)
    if(associated(leaf%patch_wetind))    nullify (leaf%patch_wetind)
    if(associated(leaf%leaf_class))      nullify (leaf%leaf_class)

    ! TEB_SPM
    if (TEB_SPM==1) then
       if(associated(leaf%G_URBAN))      nullify (leaf%G_URBAN)
    endif

    if(associated(leaf%soil_rough))      nullify (leaf%soil_rough)
    if(associated(leaf%sfcwater_nlev))   nullify (leaf%sfcwater_nlev)
    if(associated(leaf%stom_resist))     nullify (leaf%stom_resist)

    if(associated(leaf%ground_rsat))     nullify (leaf%ground_rsat)
    if(associated(leaf%ground_rvap))     nullify (leaf%ground_rvap)

    if(associated(leaf%veg_water))       nullify (leaf%veg_water)
    if(associated(leaf%veg_temp))        nullify (leaf%veg_temp)

    if(associated(leaf%can_rvap))        nullify (leaf%can_rvap)
    if(associated(leaf%can_temp))        nullify (leaf%can_temp)

    if(associated(leaf%veg_ndvip))       nullify (leaf%veg_ndvip)
    if(associated(leaf%veg_ndvic))       nullify (leaf%veg_ndvic)
    if(associated(leaf%veg_ndvif))       nullify (leaf%veg_ndvif)

    if(associated(leaf%R_aer))           nullify (leaf%R_aer)   !kml drydep

    if(associated(leaf%snow_mass))       nullify (leaf%snow_mass)
    if(associated(leaf%snow_depth))      nullify (leaf%snow_depth)
    if(associated(leaf%seatp))           nullify (leaf%seatp)
    if(associated(leaf%seatf))           nullify (leaf%seatf)

  end subroutine nullify_leaf

  ! ********************************************************************

  subroutine dealloc_leaf(leaf)

    ! TEB_SPM
    use teb_spm_start, only: TEB_SPM ! INTENT(IN)

    implicit none
    type (leaf_vars) :: leaf

    if(associated(leaf%soil_water))      deallocate (leaf%soil_water)
    if(associated(leaf%soil_energy))     deallocate (leaf%soil_energy)
    if(associated(leaf%soil_text))       deallocate (leaf%soil_text)

    if(associated(leaf%sfcwater_mass))   deallocate (leaf%sfcwater_mass)
    if(associated(leaf%sfcwater_energy)) deallocate (leaf%sfcwater_energy)
    if(associated(leaf%sfcwater_depth))  deallocate (leaf%sfcwater_depth)

    if(associated(leaf%ustar))           deallocate (leaf%ustar)
    if(associated(leaf%tstar))           deallocate (leaf%tstar)
    if(associated(leaf%rstar))           deallocate (leaf%rstar)

    if(associated(leaf%veg_fracarea))    deallocate (leaf%veg_fracarea)
    if(associated(leaf%veg_lai))         deallocate (leaf%veg_lai)
    if(associated(leaf%veg_rough))       deallocate (leaf%veg_rough)
    if(associated(leaf%veg_height))      deallocate (leaf%veg_height)
    if(associated(leaf%veg_albedo))      deallocate (leaf%veg_albedo)
    if(associated(leaf%veg_tai))         deallocate (leaf%veg_tai)

    if(associated(leaf%patch_area))      deallocate (leaf%patch_area)
    if(associated(leaf%patch_rough))     deallocate (leaf%patch_rough)
    if(associated(leaf%patch_wetind))    deallocate (leaf%patch_wetind)
    if(associated(leaf%leaf_class))      deallocate (leaf%leaf_class)

    ! TEB_SPM
    if (TEB_SPM==1) then
       if(associated(leaf%G_URBAN))      deallocate (leaf%G_URBAN)
    endif

    if(associated(leaf%soil_rough))      deallocate (leaf%soil_rough)
    if(associated(leaf%sfcwater_nlev))   deallocate (leaf%sfcwater_nlev)
    if(associated(leaf%stom_resist))     deallocate (leaf%stom_resist)

    if(associated(leaf%ground_rsat))     deallocate (leaf%ground_rsat)
    if(associated(leaf%ground_rvap))     deallocate (leaf%ground_rvap)

    if(associated(leaf%veg_water))       deallocate (leaf%veg_water)
    if(associated(leaf%veg_temp))        deallocate (leaf%veg_temp)

    if(associated(leaf%can_rvap))        deallocate (leaf%can_rvap)
    if(associated(leaf%can_temp))        deallocate (leaf%can_temp)

    if(associated(leaf%veg_ndvip))       deallocate (leaf%veg_ndvip)
    if(associated(leaf%veg_ndvic))       deallocate (leaf%veg_ndvic)
    if(associated(leaf%veg_ndvif))       deallocate (leaf%veg_ndvif)

    if(associated(leaf%R_aer))           deallocate (leaf%R_aer)   !kml drydep

    if(associated(leaf%snow_mass))       deallocate (leaf%snow_mass)
    if(associated(leaf%snow_depth))      deallocate (leaf%snow_depth)
    if(associated(leaf%seatp))           deallocate (leaf%seatp)
    if(associated(leaf%seatf))           deallocate (leaf%seatf)

  end subroutine dealloc_leaf

  ! ********************************************************************

  subroutine filltab_leaf(leaf,leafm,imean,nz,nx,ny,nzg,nzs,np,ng)
    ! TEB_SPM
    use teb_spm_start, only: TEB_SPM ! INTENT(IN)
    ! ALF
    use io_params, only: ipastin ! INTENT(IN)
    use var_tables, only: InsertVTab

    implicit none
    include "i8.h"
    type (leaf_vars) :: leaf,leafm
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

    npts=nzg*nx*ny*np
    call InsertVTab (leaf%soil_water,leafm%soil_water  &
         ,ng, npts, imean,  &
         'SOIL_WATER :4:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (leaf%soil_energy,leafm%soil_energy  &
         ,ng, npts, imean,  &
         'SOIL_ENERGY :4:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (leaf%soil_text,leafm%soil_text  &
         ,ng, npts, imean,  &
         'SOIL_TEXT :4:hist:anal:mpti:mpt3'//trim(str_recycle))

    npts=nzs*nx*ny*np
    call InsertVTab (leaf%sfcwater_mass,leafm%sfcwater_mass  &
         ,ng, npts, imean,  &
         'SFCWATER_MASS :5:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (leaf%sfcwater_energy, leafm%sfcwater_energy &
         ,ng, npts, imean,  &
         'SFCWATER_ENERGY :5:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (leaf%sfcwater_depth,leafm%sfcwater_depth  &
         ,ng, npts, imean,  &
         'SFCWATER_DEPTH :5:hist:anal:mpti:mpt3'//trim(str_recycle))

    npts=nx*ny*np
    call InsertVTab (leaf%ustar,leafm%ustar  &
         ,ng, npts, imean,  &
         'USTAR :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (leaf%tstar,leafm%tstar  &
         ,ng, npts, imean,  &
         'TSTAR :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (leaf%rstar,leafm%rstar  &
         ,ng, npts, imean,  &
         'RSTAR :6:hist:anal:mpti:mpt3'//trim(str_recycle))

    call InsertVTab (leaf%veg_fracarea,leafm%veg_fracarea  &
         ,ng, npts, imean,  &
         'VEG_FRACAREA :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (leaf%veg_lai,leafm%veg_lai  &
         ,ng, npts, imean,  &
         'VEG_LAI :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (leaf%veg_rough,leafm%veg_rough  &
         ,ng, npts, imean,  &
         'VEG_ROUGH :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (leaf%veg_height,leafm%veg_height  &
         ,ng, npts, imean,  &
         'VEG_HEIGHT :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (leaf%veg_albedo,leafm%veg_albedo  &
         ,ng, npts, imean,  &
         'VEG_ALBEDO :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (leaf%veg_tai,leafm%veg_tai  &
         ,ng, npts, imean,  &
         'VEG_TAI :6:hist:anal:mpti:mpt3'//trim(str_recycle))

    call InsertVTab (leaf%patch_area,leafm%patch_area  &
         ,ng, npts, imean,  &
         'PATCH_AREA :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (leaf%patch_rough,leafm%patch_rough  &
         ,ng, npts, imean,  &
         'PATCH_ROUGH :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (leaf%patch_wetind,leafm%patch_wetind  &
         ,ng, npts, imean,  &
         'PATCH_WETIND :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (leaf%leaf_class,leafm%leaf_class  &
         ,ng, npts, imean,  &
         'LEAF_CLASS :6:hist:anal:mpti:mpt3'//trim(str_recycle))

    ! TEB_SPM
    if (TEB_SPM==1) then
       call InsertVTab (leaf%G_URBAN,leafm%G_URBAN  &
            ,ng, npts, imean,  &
            'G_URBAN :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    endif

    call InsertVTab (leaf%soil_rough,leafm%soil_rough  &
         ,ng, npts, imean,  &
         'SOIL_ROUGH :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (leaf%sfcwater_nlev,leafm%sfcwater_nlev  &
         ,ng, npts, imean,  &
         'SFCWATER_NLEV :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (leaf%stom_resist,leafm%stom_resist  &
         ,ng, npts, imean,  &
         'STOM_RESIST :6:hist:anal:mpti:mpt3'//trim(str_recycle))

    call InsertVTab (leaf%ground_rsat,leafm%ground_rsat  &
         ,ng, npts, imean,  &
         'GROUND_RSAT :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (leaf%ground_rvap,leafm%ground_rvap  &
         ,ng, npts, imean,  &
         'GROUND_RVAP :6:hist:anal:mpti:mpt3'//trim(str_recycle))

    call InsertVTab (leaf%veg_water,leafm%veg_water  &
         ,ng, npts, imean,  &
         'VEG_WATER :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (leaf%veg_temp,leafm%veg_temp  &
         ,ng, npts, imean,  &
         'VEG_TEMP :6:hist:anal:mpti:mpt3'//trim(str_recycle))

    call InsertVTab (leaf%can_rvap,leafm%can_rvap  &
         ,ng, npts, imean,  &
         'CAN_RVAP :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (leaf%can_temp,leafm%can_temp  &
         ,ng, npts, imean,  &
         'CAN_TEMP :6:hist:anal:mpti:mpt3'//trim(str_recycle))

    call InsertVTab (leaf%veg_ndvip,leafm%veg_ndvip  &
         ,ng, npts, imean,  &
         'VEG_NDVIP :6:hist:mpti')
    call InsertVTab (leaf%veg_ndvic,leafm%veg_ndvic  &
         ,ng, npts, imean,  &
         'VEG_NDVIC :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    call InsertVTab (leaf%veg_ndvif,leafm%veg_ndvif  &
         ,ng, npts, imean,  &
         'VEG_NDVIF :6:hist:mpti')

    call InsertVTab (leaf%R_aer,leafm%R_aer  &      !kml drydep
         ,ng, npts, imean,  &                          !kml drydep
         'R_aer :6:hist:mpti')                         !kml drydep

    npts=nx*ny
    call InsertVTab (leaf%snow_mass,leafm%snow_mass  &
         ,ng, npts, imean,  &
         'SNOW_MASS :2:mpti')
    call InsertVTab (leaf%snow_depth,leafm%snow_depth  &
         ,ng, npts, imean,  &
         'SNOW_DEPTH :2:mpti')
    call InsertVTab (leaf%seatp,leafm%seatp  &
         ,ng, npts, imean,  &
         'SEATP :2:mpti')
    call InsertVTab (leaf%seatf,leafm%seatf  &
         ,ng, npts, imean,  &
         'SEATF :2:mpti')

  end subroutine filltab_leaf

  subroutine StoreNamelistFileAtMem_leaf(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    albedo = oneNamelistFile%albedo
    drtcon = oneNamelistFile%drtcon
    dthcon = oneNamelistFile%dthcon
    isfcl = oneNamelistFile%isfcl
    isfcl_ocean = oneNamelistFile%isfcl_ocean
    nslcon = oneNamelistFile%nslcon
    nvegpat = oneNamelistFile%nvegpat
    nvgcon = oneNamelistFile%nvgcon
    pctlcon = oneNamelistFile%pctlcon
    seatmp = oneNamelistFile%seatmp
    slmstr = oneNamelistFile%slmstr
    slz = oneNamelistFile%slz
    stgoff = oneNamelistFile%stgoff
    zrough = oneNamelistFile%zrough
    !-- if the ocean model is undefined, use the one 
    !-- in JULES or LEAF.
    if(isfcl_ocean == -999 ) then
       isfcl_ocean = isfcl
    endif
  end subroutine StoreNamelistFileAtMem_leaf
End Module mem_leaf
