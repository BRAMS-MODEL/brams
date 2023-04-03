MODULE ChemDryDepDriver


  USE rconstants, ONLY: &
       cpi,             &
       cpor,            &
       p00,             &
       g,               &
       vonk
       
  USE mem_grid, ONLY: &
       grid_g,        &
       jdim,          &
       dzt,           &
       zt,            &
       nzpmax,        &
       npatch,        &
       dtlt,          &
       imonth1,       &
       idate1,        &
       iyear1,        &
       ngrid

  USE micphys, ONLY: &
       level

  USE mem_cuparm, ONLY: &
       cuparm_g,        &
       nnqparm

  USE mem_basic, ONLY: &
       basic_g

  USE mem_turb, ONLY: &
       turb_g

  USE mem_leaf, ONLY: &
       leaf_g

  USE mem_micro, ONLY: &
       micro_g

  USE mem_radiate, ONLY: &
       radiate_g

  USE mem_chem1, ONLY: &
       chem1_g,        &
       chemistry
  
  USE mem_aer1, ONLY: &
       aerosol, &
       aer1_g

  USE module_dry_dep, ONLY: &
       dd_sedim,            &
       dry_dep                 ! Subroutine



  IMPLICIT NONE

  PRIVATE


  PUBLIC :: drydep_driver



CONTAINS


  !========================================================================
  SUBROUTINE drydep_driver(m1,m2,m3,ia,iz,ja,jz)

    INTEGER,              INTENT(IN)    :: m1
    INTEGER,              INTENT(IN)    :: m2
    INTEGER,              INTENT(IN)    :: m3
    INTEGER,              INTENT(IN)    :: ia
    INTEGER,              INTENT(IN)    :: iz
    INTEGER,              INTENT(IN)    :: ja
    INTEGER,              INTENT(IN)    :: jz


    RETURN
  END SUBROUTINE drydep_driver
  !========================================================================


END MODULE ChemDryDepDriver
