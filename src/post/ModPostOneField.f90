module ModPostOneField

   use ModBramsGrid, only : BramsGrid
   use ModPostTypes, only : PostGrid
   use ModPostOneField2d, only : Brams2Post_2d
   use ModPostOneField3d, only : Brams2Post_3d
   use ModPostOneField7d, only : Brams2Post_7d
   use ModPostOneField8d, only : Brams2Post_8d

   use ModPostOneFieldUtils, only : initialize_all_post_variables
   use ModPostOneFieldUtils, only : finalize_all_post_variables
   use ModPostOneFieldUtils, only : add_post_variable
   !use ModPostOneFieldUtils, only : getPostVarible
   !use ModPostOneFieldUtils, only : PostVarType
   use ModPostTypes, only: PostVarType
   use ModPostTypes, only : getPostVarible
   use ModPostUtils, only : UpperCase
   use ModNamelistFile, only: namelistFile

   implicit none

   private

   public :: PostOneField
   public :: initialize_post_variables
   public :: finalize_post_variables


contains

   ! TODO - Read the variables from a file to exchange with Brabu
   subroutine initialize_post_variables(oneNamelistFile)
      type(NamelistFile), pointer :: oneNamelistFile
      type(PostVarType) :: postVar
      integer :: varfileUnit, err

      call initialize_all_post_variables()
      varfileUnit = 135
      open(varfileUnit, file = oneNamelistFile%csvFile, status = "old", action = "read")
      do
         read(varfileUnit, *, iostat = err) postVar%ivar_type, postVar%fieldName, postVar%fieldDescription, postVar%fieldUnits
         if (err /= 0) then
            exit
         end if
         call add_post_variable(postVar%ivar_type, postVar%fieldName, postVar%fieldDescription, postVar%fieldUnits)
      enddo
      close(varfileUnit)

   end subroutine initialize_post_variables


   subroutine finalize_post_variables()
      call finalize_all_post_variables()
   end subroutine


subroutine PostOneField(varName, oneBramsGrid, onePostGrid, oneNamelistFile)
      use dump
      use node_mod, only: &
       mchnum,        &
       master_num

      include "constants.f90"
      character(len = *), intent(in) :: varName
      type(NamelistFile), pointer :: oneNamelistFile
      type(BramsGrid), pointer :: oneBramsGrid
      type(PostGrid), pointer :: onePostGrid
      type(PostVarType) :: one_post_variable
      character(len = 16) :: varNameUpper
      integer :: err

      varNameUpper = trim(UpperCase(varName))
      one_post_variable = getPostVarible(varNameUpper)
      if (len(trim(one_post_variable%fieldName)) .eq. 0) then
         if (mchnum==master_num) &
            err=dumpMessage(c_tty,c_yes,'Post','',c_warning,'Post field "' // varName &
               // '" does not exists in list of variables. Model will continue!')
      else
         select case (one_post_variable%ivar_type)
         case (2)
            call Brams2Post_2d(one_post_variable, oneBramsGrid, onePostGrid)
         case (3)
            call Brams2Post_3d(one_post_variable, oneBramsGrid, onePostGrid, oneNamelistFile)
         case (7)
            call Brams2Post_7d(one_post_variable, oneBramsGrid, onePostGrid)
         case (8)
            call Brams2Post_8d(one_post_variable, oneBramsGrid, onePostGrid)
         end select
      end if

   end subroutine PostOneField

end module ModPostOneField
