!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


Module chem_isan_coms
  use    isan_coms , only : maxagrds
  use    chem1_list, only : total_nspecies=>nspecies, spc_alloc, fdda
  use    aer1_list , only : total_speciesAer=>nspecies,nmodes,spc_alloc_aer=>spc_alloc

  !---------------------------------------------------------------------------
  !    Configuration COMMON blocks for RAMS isentropic data analysis package.
  !---------------------------------------------------------------------------
  !MAXPR         Maximum number of vertical levels that can be used in 
  !                 the pressure data.
  !MAXISN        Maximum number of vertical levels that can be used in 
  !                 the isentropic analysis.
  !MAXX          Maximum number of grid points (in the x or east-west direction) 
  !                 in any of the RAMS or pressure grids.
  !MAXY          Maximum number of grid points (in the y or north-south direction) 
  !                 in any of the RAMS or pressure grids.
  !MAXTIMES      Maximum number of data analysis times that can be processed 
  !                 in a single run.
  !MAXAGRDS      Maximum number of RAMS grids that can have varfiles generated.
  !MAXSIGZ       Maximum number of vertical levels that can be used 
  !                 in the  _z analysis.
  !MAXLEV	       Maximum number of levels in an input rawinsonde.
  !MAXSNAME      Maximum number of input observations
  !MAXISFILES    Maximum number of input data times
  !------------------------------------------------------------------------------------



  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Input pressure data memory
  integer :: nspecies
  integer :: nspecies_aer_in
!srf-chem  
  real, allocatable, dimension(:,:,:,:)   :: p_sc

  real, allocatable, dimension(:,:,:,:)   :: p_aer_sc


  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Polar-stereo/pressure grid memory

!srf-chem  
  real, allocatable, dimension(:,:,:,:)   :: pp_sc
  real, allocatable, dimension(:,:,:,:)   :: pp_aer_sc
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Polar-stereo/isentropic grid memory

!srf-chem  
  real, allocatable, dimension(:,:,:,:)   :: pi_sc

 !srf-chem  
  real, allocatable, dimension(:,:,:,:)   :: pi_aer_sc 
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Polar-stereo/sigma-z grid memory
  !                         :: 
!srf-chem  
  real, allocatable, dimension(:,:,:,:)   :: ps_sc

!lfr-chem  
  real, allocatable, dimension(:,:,:,:)   :: ps_aer_sc
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Polar-stereo/surface grid memory
   ! ...
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Data type to replace A array memory use in ISAN. 

  type chem_isan_grids
!srf-chem  
     real, pointer, dimension(:,:,:,:) :: rr_sc,rr_scg,rr_sc0  

  end type chem_isan_grids

  type (chem_isan_grids)               :: chem_is_grids(maxagrds)
  type (chem_isan_grids)               :: aer_is_grids(maxagrds)


  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Input observation data memory
  !
  ! ...
End module chem_isan_coms
