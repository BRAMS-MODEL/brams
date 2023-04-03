module ModPostProcess

   use ModParallelEnvironment, only : &
         MsgDump

   use ModPostOneField, only : PostOneField
   use ModPostOneField, only : initialize_post_variables
   use ModPostOneField, only : finalize_post_variables

   use ModNamelistFile, only : namelistFile

   use ModBramsGrid, only : &
         BramsGrid, CreateBramsGrid, DestroyBramsGrid

   use ModPostGrid, only : &
         CreatePostGrid, DestroyPostGrid, &
         OpenGradsBinaryFile, CloseGradsBinaryFile, &
         OpenGradsControlFile, CloseGradsControlFile, &
         FillGradsControlFile, UpdateVerticals

   use ModPostTypes, only: &
      PostGrid
#ifdef cdf
   use ModPOstGridNetCDF, only: &
         FillNetcdfVarControlFile, &
         ncid, OpenNetCDFBinaryFile

   use netcdf
#endif
   use ModGridTree, only : &
         GridTree, &
         GridTreeRoot, &
         NextOnGridTree

   use ModGrid, only : &
         Grid

   use ModMessageSet, only : &
         PostRecvSendMsgs, &
         WaitRecvMsgs

     use io_params, only : & ! 
       ipos
       
   use modTimeLineFRN, only: &
      isTimeTograds

   use mem_grid, only: &
      timmax, &
      dtlongn, &
      time



   !use CUPARM_GRELL3, only: g3d_g

   implicit none
   include "tsNames.h"
   private
   public :: AllPostTypes
   public :: CreatePostProcess
   public :: PostProcess
   public :: DestroyPostProcess

   type PostTypePair
      type(BramsGrid), pointer :: bg
      type(PostGrid), pointer :: pg
   end type PostTypePair

   type AllPostTypes
      type(PostTypePair), allocatable :: allGrids(:)
   end type AllPostTypes

   logical, parameter :: dumpLocal = .false.

contains


   subroutine CreatePostProcess(oneNamelistFile, oneAllPostTypes)
      type(namelistFile), pointer :: oneNamelistFile
      type(AllPostTypes), pointer :: oneAllPostTypes

      integer :: igrid
      integer :: ierr
      character(len = 8) :: c0
      character(len = *), parameter :: h = "**(CreatePostProcess)**"

      if (.not. associated(oneNamelistFile)) then
         call fatal_error(h // " invoked with null oneNamelistFile")
      else if (associated(oneAllPostTypes)) then
         call fatal_error(h // " invoked with already associated oneAllPostTypes")
      end if

      allocate(oneAllPostTypes, stat = ierr)
      if (ierr /= 0) then
         call fatal_error(h // " allocate oneAllPostTypes fails")
      end if

      allocate(oneAllPostTypes%allGrids(oneNamelistFile%ngrids), stat = ierr)
      if (ierr /= 0) then
         call fatal_error(h // " allocate oneAllPostTypes%allGrids fails")
      end if

      do igrid = 1, oneNamelistFile%ngrids
         if (dumpLocal) then
            write(c0, "(i8)") igrid
            call MsgDump(h // " will create types for grid " // &
                  trim(adjustl(c0)))
         end if

         oneAllPostTypes%allGrids(igrid)%bg => null()
         call CreateBramsGrid(oneAllPostTypes%allGrids(igrid)%bg, igrid)

         oneAllPostTypes%allGrids(igrid)%pg => null()
         call CreatePostGrid(oneNamelistFile, &
               oneAllPostTypes%allGrids(igrid)%bg, &
               oneAllPostTypes%allGrids(igrid)%pg, &
               igrid)
      end do

   end subroutine CreatePostProcess


   subroutine PostProcess(AllGrids, oneNamelistFile, oneAllPostTypes)
      type(GridTree), pointer :: AllGrids
      type(namelistFile), pointer :: oneNamelistFile
      type(AllPostTypes), pointer :: oneAllPostTypes

      include "constants.f90"

      integer :: igrid, ivp, ierr
      character(len = 8) :: c0, c1
      type(GridTree), pointer :: OneGridTreeNode => null ()
      type(Grid), pointer :: OneGrid => null()
      character(len = *), parameter :: h = "**(PostProcess)**"
      if (.not. associated(oneNamelistFile)) then
         call fatal_error(h // " invoked with null oneNamelistFile")
      else if (.not. associated(oneAllPostTypes)) then
         call fatal_error(h // " invoked with null oneAllPostTypes")
      else if (.not. allocated(oneAllPostTypes%allGrids)) then
         call fatal_error(h // " invoked with null oneAllPostTypes%allGrids")
      else if (size(oneAllPostTypes%allGrids) /= oneNamelistFile%ngrids) then
         write(c0, "(i8)") size(oneAllPostTypes%allGrids)
         write(c1, "(i8)") oneNamelistFile%ngrids
         call fatal_error(h // " number of grids at oneAllPostTypes (" // &
               trim(adjustl(c0)) // ") differs from required (" // trim(adjustl(c1)) // ")")
      else

         do igrid = 1, oneNamelistFile%ngrids
            if (.not. associated(oneAllPostTypes%allGrids(igrid)%bg)) then
               write(c0, "(i8)") igrid
               call fatal_error(h // " Brams grid type for grid " // &
                     trim(adjustl(c0)) // " not created")
            else if (.not. associated(oneAllPostTypes%allGrids(igrid)%pg)) then
               write(c0, "(i8)") igrid
               call fatal_error(h // " Post grid type for grid " // &
                     trim(adjustl(c0)) // " not created")
            end if
         end do
      end if
      ! for each grid

      OneGridTreeNode => GridTreeRoot(AllGrids)
      do while (associated(OneGridTreeNode))

         OneGrid => OneGridTreeNode%curr

         ! update Ghost Zone of all vartables variables part 1:
         ! post receives and send messages
         call PostRecvSendMsgs(&
               OneGrid%AllGhostZoneSend, &
               OneGrid%AllGhostZoneRecv)
         igrid = OneGrid%Id

         ! open grads files
         if(IPOS==2 .or. isTimeTograds()) then
            call OpenGradsBinaryFile(oneNamelistFile, &
               oneAllPostTypes%allGrids(igrid)%pg, &
               oneAllPostTypes%allGrids(igrid)%bg, igrid)
            call OpenGradsControlFile(oneNamelistFile, &
               oneAllPostTypes%allGrids(igrid)%pg, &
               oneAllPostTypes%allGrids(igrid)%bg, igrid)
         endif
#ifdef cdf
         if(IPOS==3) then
            call OpenNetCDFBinaryFile(oneNamelistFile, &
               oneAllPostTypes%allGrids(igrid)%pg, &
               oneAllPostTypes%allGrids(igrid)%bg, igrid)
         endif
#endif
         ! master_process dumps

         if (dumpLocal) then
            call MsgDump(h // " Creating Post Processed file " // &
                  trim(oneAllPostTypes%allGrids(igrid)%pg%binFileName))
         end if

         ! update Ghost Zone of all vartables variables part 2:
         ! receive messages and wait on sends

         call WaitRecvMsgs(&
               OneGrid%AllGhostZoneSend, &
               OneGrid%AllGhostZoneRecv)

         ! update verticals according to namelist options

         call UpdateVerticals(&
               oneAllPostTypes%allGrids(igrid)%bg, &
               oneAllPostTypes%allGrids(igrid)%pg)
         ! post process each desired field and
         ! write resulting field to grads binary file
         call initialize_post_variables(oneNamelistFile)
#ifdef cdf
         if(IPOS==3) then
            call FillNetcdfVarControlFile(oneNamelistFile, &
               oneAllPostTypes%allGrids(igrid)%pg, &
               oneAllPostTypes%allGrids(igrid)%bg)
         endif
#endif
         do ivp = 1, oneNamelistFile%nvp
            if (dumpLocal) then
               call MsgDump (h // " variable " // &
                     trim(oneNamelistFile%vp(ivp)))
            end if
            call PostOneField(trim(oneNamelistFile%vp(ivp)), &
                  oneAllPostTypes%allGrids(igrid)%bg, &
                  oneAllPostTypes%allGrids(igrid)%pg,oneNamelistFile)
         end do
         call finalize_post_variables()
         ! control file contents

         if(IPOS==2 .or. isTimeTograds()) then
            call FillGradsControlFile(&
               oneAllPostTypes%allGrids(igrid)%pg, &
               oneAllPostTypes%allGrids(igrid)%bg)

            ! close grads files
            call CloseGradsBinaryFile(oneAllPostTypes%allGrids(igrid)%pg, &
               oneAllPostTypes%allGrids(igrid)%bg)

            call CloseGradsControlFile(oneAllPostTypes%allGrids(igrid)%pg, &
               oneAllPostTypes%allGrids(igrid)%bg)
         endif
#ifdef cdf 
         if(IPOS==3) ierr=nf90_close(ncid)
#endif         
         OneGridTreeNode => NextOnGridTree(OneGridTreeNode)
         

      end do

   end subroutine PostProcess


   subroutine DestroyPostProcess(oneNamelistFile, oneAllPostTypes)
      type(namelistFile), pointer :: oneNamelistFile
      type(AllPostTypes), pointer :: oneAllPostTypes

      integer :: igrid
      integer :: ierr
      character(len = 8) :: c0, c1
      character(len = *), parameter :: h = "**(DestroyPostProcess)**"

      if (.not. associated(oneNamelistFile)) then
         call fatal_error(h // " invoked with null oneNamelistFile")
      else if (.not. associated(oneAllPostTypes)) then
         return
      else if (.not. allocated(oneAllPostTypes%allGrids)) then
         return
      else if (size(oneAllPostTypes%allGrids) /= oneNamelistFile%ngrids) then
         write(c0, "(i8)") size(oneAllPostTypes%allGrids)
         write(c1, "(i8)") oneNamelistFile%ngrids
         call fatal_error(h // " number of grids at oneAllPostTypes (" // &
               trim(adjustl(c0)) // ") differs from required (" // trim(adjustl(c1)) // ")")
      else
         do igrid = 1, oneNamelistFile%ngrids
            if (.not. associated(oneAllPostTypes%allGrids(igrid)%bg)) then
               write(c0, "(i8)") igrid
               call fatal_error(h // " Brams grid type for grid " // &
                     trim(adjustl(c0)) // " not created")
            end if
            call DestroyBramsGrid(oneAllPostTypes%allGrids(igrid)%bg)

            if (.not. associated(oneAllPostTypes%allGrids(igrid)%pg)) then
               write(c0, "(i8)") igrid
               call fatal_error(h // " Post grid type for grid " // &
                     trim(adjustl(c0)) // " not created")
            end if
            call DestroyPostGrid(oneAllPostTypes%allGrids(igrid)%pg)

         end do
      end if

      if (dumpLocal) then
         call MsgDump(h // " executing ")
      end if

      deallocate(oneAllPostTypes%allGrids, stat = ierr)
      if (ierr /= 0) then
         call fatal_error(h // " deallocate oneAllPostTypes%allGrids fails")
      end if

      deallocate(oneAllPostTypes, stat = ierr)
      if (ierr /= 0) then
         call fatal_error(h // " deallocate oneAllPostTypes fails")
      end if

      oneAllPostTypes => null()
   end subroutine DestroyPostProcess
end module ModPostProcess
