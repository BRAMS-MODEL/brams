module ModPostOneFieldUtils

   use ModOutputUtils, only : GetVarFromMemToOutput

   use ModPostGrid, only : OutputGradsField
   use ModPostTypes, only : PostGrid
   use ModBramsGrid, only : BramsGrid
   use ModPostTypes, only: &
      PostVarType, &
      all_post_variables

   implicit none

   private

   !type(PostVarType), allocatable :: all_post_variables(:)

   public :: PrepareAndOutputGradsField
   interface PrepareAndOutputGradsField
      module procedure PrepareAndOutputGradsField_2D
      module procedure PrepareAndOutputGradsField_3D
      module procedure PrepareAndOutputGradsField_4D
   end interface

   public :: initialize_all_post_variables
   public :: finalize_all_post_variables
   public :: add_post_variable
   !public :: getPostVarible
   public :: PrepareGradsField

contains


   subroutine PrepareGradsField(one_post_variable, onePostGrid)
      type(PostVarType), intent(inout) :: one_post_variable
      type(PostGrid), intent(inout), pointer :: onePostGrid

      onePostGrid%ivar_type = one_post_variable%ivar_type
      onePostGrid%fieldName = one_post_variable%fieldName
      onePostGrid%fieldDescription = one_post_variable%fieldDescription
      onePostGrid%fieldUnits = one_post_variable%fieldUnits
   end subroutine PrepareGradsField


   subroutine PrepareAndOutputGradsField_2d(one_post_variable, oneBramsGrid, onePostGrid, OutputField)
      type(PostVarType) :: one_post_variable
      type(BramsGrid), pointer :: oneBramsGrid
      type(PostGrid), pointer :: onePostGrid
      real, intent(in) :: OutputField(:, :)

      call PrepareGradsField(one_post_variable, onePostGrid)
      call OutputGradsField(oneBramsGrid, onePostGrid, OutputField)
   end subroutine PrepareAndOutputGradsField_2d


   subroutine PrepareAndOutputGradsField_3d(one_post_variable, oneBramsGrid, onePostGrid, OutputField)
      type(PostVarType) :: one_post_variable
      type(BramsGrid), pointer :: oneBramsGrid
      type(PostGrid), pointer :: onePostGrid
      real, intent(inout) :: OutputField(:, :, :)

      call PrepareGradsField(one_post_variable, onePostGrid)
      call OutputGradsField(oneBramsGrid, onePostGrid, OutputField)
   end subroutine PrepareAndOutputGradsField_3d


   subroutine PrepareAndOutputGradsField_4d(one_post_variable, oneBramsGrid, onePostGrid, OutputField)
      type(PostVarType) :: one_post_variable
      type(BramsGrid), pointer :: oneBramsGrid
      type(PostGrid), pointer :: onePostGrid
      real, intent(in) :: OutputField(:, :, :, :)

      call PrepareGradsField(one_post_variable, onePostGrid)
      call OutputGradsField(oneBramsGrid, onePostGrid, OutputField)
   end subroutine PrepareAndOutputGradsField_4d


   subroutine initialize_all_post_variables()
      if (.not. allocated(all_post_variables)) allocate(all_post_variables(0))
   end subroutine


   subroutine finalize_all_post_variables()
      deallocate(all_post_variables)
   end subroutine


   subroutine add_post_variable(varType, varName, varDescription, varUnit)
      integer, intent(in) :: varType        ! Dimension
      character(len = *), intent(in) :: varName        ! post field name
      character(len = *), intent(in) :: varDescription ! post field description
      character(len = *), intent(in) :: varUnit       ! post field units
      type(PostVarType) :: new_post_variable
      type(PostVarType), allocatable :: all_post_variables_copy(:)
      integer :: all_size

      new_post_variable%ivar_type = varType
      new_post_variable%fieldName = varName
      new_post_variable%fieldDescription = varDescription
      new_post_variable%fieldUnits = varUnit
      all_size = size(all_post_variables) + 1
      if(all_size .gt. 1) then
         allocate(all_post_variables_copy(all_size - 1))
         all_post_variables_copy = all_post_variables
         deallocate(all_post_variables)
         allocate(all_post_variables(all_size))
         all_post_variables(1 : all_size - 1) = all_post_variables_copy(1 : all_size - 1)
         deallocate(all_post_variables_copy)
      else
         deallocate(all_post_variables)
         allocate(all_post_variables(all_size))
      end if
      all_post_variables(all_size) = new_post_variable

   end subroutine add_post_variable


end module ModPostOneFieldUtils
