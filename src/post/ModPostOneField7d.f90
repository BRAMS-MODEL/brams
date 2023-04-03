module ModPostOneField7d

   use ModOutputUtils, only : GetVarFromMemToOutput
   use ModBramsGrid, only : BramsGrid
   use ModPostTypes, only : PostGrid
   use ModPostOneFieldUtils, only : PrepareAndOutputGradsField
   !use ModPostOneFieldUtils, only : PostVarType
   use ModPostTypes, only: PostVarType

   use ModPostUtils, only : rams_comp_tempc
   use ModPostUtils, only : rams_comp_vegclass

   implicit none

   private

   public :: Brams2Post_7d

contains


   subroutine Brams2Post_7d (one_post_variable, oneBramsGrid, onePostGrid)
      type(PostVarType) :: one_post_variable
      type(BramsGrid), pointer :: oneBramsGrid
      type(PostGrid), pointer :: onePostGrid

      real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)

      select case (one_post_variable%fieldName)
      case ('TVEG')
         call GetVarFromMemToOutput ('VEG_TEMP', oneBramsGrid%currGrid, OutputField)
         call rams_comp_tempc (OutputField)
      case ('VTYPE')
         call GetVarFromMemToOutput ('LEAF_CLASS', oneBramsGrid%currGrid, OutputField)
         call rams_comp_vegclass(OutputField)
      case ('LAI')
         call GetVarFromMemToOutput ('VEG_LAI', oneBramsGrid%currGrid, OutputField)
      case ('H_ANTRHJ')
         call GetVarFromMemToOutput ('anthrop_heatj', oneBramsGrid%currGrid, OutputField)
      case ('RADNETJ')
         call GetVarFromMemToOutput ('radnet_tilej', oneBramsGrid%currGrid, OutputField)
      case ('H_TILEJ')
         call GetVarFromMemToOutput ('ftl_tilej', oneBramsGrid%currGrid, OutputField)
      case ('LE_TILEJ')
         call GetVarFromMemToOutput ('le_tilej', oneBramsGrid%currGrid, OutputField)
      case ('HTF_TILEJ')
         call GetVarFromMemToOutput ('htf_tilej', oneBramsGrid%currGrid, OutputField)
      case ('NDVI')
         call GetVarFromMemToOutput ('VEG_NDVIC', oneBramsGrid%currGrid, OutputField)
      case default
         write(*, "(a)") "**(OnePostField)** Post field 7d " // one_post_variable%fieldName // " not implemented!"
         stop
      end select

      ! most of cases run this. Some just returns before
      call PrepareAndOutputGradsField(one_post_variable, oneBramsGrid, onePostGrid, OutputField)

   end subroutine Brams2Post_7d


end module ModPostOneField7d
