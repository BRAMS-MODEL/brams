!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
! 
!---  M. Pirre 12/08/08 ------------------------------

module mem_chemic

   type chemic_vars
   
      ! Variables to be dimensioned by (nzp,nxp,nyp)
   real, pointer, dimension(:,:,:) :: coll,sedimr 

                          
   end type               
                          
   type (chemic_vars), allocatable :: chemic_g(:)
                          
contains                  
                          
  subroutine alloc_chemic(chemic,n1,n2,n3)

    use micphys

    implicit none          
    type (chemic_vars) :: chemic
    integer, intent(in) :: n1,n2,n3

    ! Allocate arrays based on options (if necessary)

    if (level >= 3) then
       if(irain >= 1)  then
          allocate (chemic%coll(n1,n2,n3)); chemic%coll=0.
          allocate (chemic%sedimr(n1,n2,n3)); chemic%sedimr=0.
       endif
    endif

    return
  end subroutine alloc_chemic


  subroutine nullify_chemic(chemic)

   implicit none
   type (chemic_vars) :: chemic 
   
   if (associated(chemic%coll))       nullify (chemic%coll)
   if (associated(chemic%sedimr))     nullify (chemic%sedimr)

   return
  end subroutine nullify_chemic

  subroutine dealloc_chemic(chemic)

   implicit none
   type (chemic_vars) :: chemic 
   
   if (associated(chemic%coll))     deallocate (chemic%coll)
   if (associated(chemic%sedimr))     deallocate (chemic%sedimr)

   return
  end subroutine dealloc_chemic

end module mem_chemic
