!############################# Change Log ##################################
!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


!=============================================================================================
module modSchedUtils
   !# Some common vars for modesched
   !#
   !# @note
   !#
   !# **Documentation/Version Control**: <documents>
   !#
   !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
   !#
   !# **Date**: 14 October 2021 (Thursday)
   !#
   !# **Full Description**: Some common vars for modesched
   !#
   !# @endnote
   !#
   !# @warning
   !# ![](https://licensebuttons.net/l/by-sa/3.0/88x31.png"")
   !# Now is under CC Attribution-ShareAlike 4.0 International, please see:
   !# &copy; <https://creativecommons.org/licenses/by-sa/4.0/>
   !# @endwarning
   !#
    
   !Use area
   use dump !Dump contains a lot of functions for debugs and formated printouts

   implicit none

   include "constants.f90"
   character(len=*),parameter :: sourceName='modsched.f90' !Name of this source code
   character(len=*),parameter :: procedureName='**modSchedUtils**' !Name of this procedure
   !

   real :: savedDtlongn
   logical :: dtlong_saved=.false.
   logical :: firstTime=.true.
   real :: last_isan_inc=0.0
   real :: last_tnudcent=0.0
   real :: last_tnudlat =0.0 
   real :: last_tnudtop =0.0 
   real :: last_frqhis  =0.0
   real :: last_frqanl  =0.0
   real :: last_radfrq  =0.0
   real :: last_confrq  =0.0
   real :: last_shcufrq =0.0
   real :: local_isan_inc
   real :: local_tnudcent
   real :: local_tnudlat 
   real :: local_tnudtop 
   real :: local_frqhis  
   real :: local_frqanl  
   real :: local_radfrq  
   real :: local_confrq  
   real :: local_shcufrq 

   real :: futureTime=0.0

   private

   public nextMainPoint

   contains

   !=============================================================================================
   function ismultiple(value) result(ismult)
      !# Verified is is multiple or divisor of hour
      !#
      !# @note
      !#
      !# **Documentation/Version Control**: <documents>
      !#
      !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 15 October 2021 (Friday)
      !#
      !# **Full Description**: Verified is is multiple or divisor of hour
      !#
      !# @endnote
      !#
      !# @warning
      !# ![](https://licensebuttons.net/l/by-sa/3.0/88x31.png"")
      !# Now is under CC Attribution-ShareAlike 4.0 International, please see:
      !# &copy; <https://creativecommons.org/licenses/by-sa/4.0/>
      !# @endwarning
      !#
       
      !Use area
      use dump !Dump contains a lot of functions for debugs and formated printouts
   
      implicit none
   
      include "constants.f90"
      character(len=*),parameter :: sourceName='modsched.f90' !Name of this source code
      character(len=*),parameter :: procedureName='**ismultiple**' !Name of this procedure
      !
      !Local Parameters
   
      !Input/Output variables
      real, intent(in) :: value
      logical :: ismult
      real :: divid

      !Code
      ismult=.false.
      if(value>3600.0) then 
         divid=mod(value,3600.)
         if(divid/=0.) ismult=.true.
      else
         divid=mod(3600.0,value)
         if(divid/=0.) ismult=.true.
      endif 

   
   end function ismultiple 

   !=============================================================================================
function nextMainPoint(ngrid) result(dt)
   !# return the time distance for the next main step point
   !#
   !# @note
   !#
   !# **Documentation/Version Control**: <documents>
   !#
   !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
   !#
   !# **Date**: 15 October 2021 (Friday)
   !#
   !# **Full Description**: return the time distance for the next main step point
   !#      The points are some special function to execute as nuding, output, radiation, etc.
   !#
   !# @endnote
   !#
   !# @warning
   !# ![](https://licensebuttons.net/l/by-sa/3.0/88x31.png"")
   !# Now is under CC Attribution-ShareAlike 4.0 International, please see:
   !# &copy; <https://creativecommons.org/licenses/by-sa/4.0/>
   !# @endwarning
   !#
    
   !Use area
   use dump !Dump contains a lot of functions for debugs and formated printouts
   use isan_coms, only: &
       isan_inc
   use mem_varinit, only: &
       tnudcent, tnudtop, tnudlat
   use io_params, only : & ! 
       frqanl,          & !INTENT(IN)
       frqlite,         & !INTENT(IN)
       frqmean,         & !INTENT(IN)
       frqhis             !INTENT(IN)
   use mem_radiate, only: &
      radfrq              ! INTENT(IN)
   use mem_cuparm, only: &
      confrq
   use shcu_vars_const, only: &
      shcufrq   
  use mem_grid, only: &
       dtlongn, &
       time

   implicit none

   include "constants.f90"
   character(len=*),parameter :: sourceName='modsched.f90' !Name of this source code
   character(len=*),parameter :: procedureName='**nextMainPoint**' !Name of this procedure
   !
   !Local Parameters

   !Input/Output variables
   integer, intent(in) :: ngrid

   real :: dt

   !Local variables
   real :: nextTime(9)
   real :: mxval
   real :: aux
   real :: timeDistance

   !Code
   !Na primeira passada deve-se eliminar os valores zerados

   if(firsttime) then 
     aux=(isan_inc/100)*3600.0
     mxval=maxval((/ &
       aux        &
      ,tnudcent   &
      ,tnudlat    &
      ,tnudtop    &
      ,frqhis     &
      ,frqanl     &
      ,radfrq     &
      ,confrq     &
      ,shcufrq    &
      /))

     local_isan_inc=aux
     local_tnudcent=tnudcent
     local_tnudlat =tnudlat 
     local_tnudtop =tnudtop 
     local_frqhis  =frqhis  
     local_frqanl  =frqanl  
     local_radfrq  =radfrq  
     local_confrq  =confrq  
     local_shcufrq =shcufrq 

     if(aux     ==0.) local_isan_inc=mxval
     if(tnudcent==0.) local_tnudcent=mxval
     if(tnudlat ==0.) local_tnudlat =mxval
     if(tnudtop ==0.) local_tnudtop =mxval
     if(frqhis  ==0.) local_frqhis  =mxval
     if(frqanl  ==0.) local_frqanl  =mxval
     if(radfrq  ==0.) local_radfrq  =mxval
     if(confrq  ==0.) local_confrq  =mxval
     if(shcufrq ==0.) local_shcufrq =mxval
     firstTime=.false.

     if(ismultiple(local_isan_inc)) iErrNumber=dumpMessage(c_tty,c_yes,'','',c_warning,'isan_inc is not multiple of 1 hour, ',local_isan_inc,"F12.1")
     if(ismultiple(local_tnudcent)) iErrNumber=dumpMessage(c_tty,c_yes,'','',c_warning,'tnudcent is not multiple of 1 hour, ',local_tnudcent,"F12.1")
     if(ismultiple(local_tnudlat )) iErrNumber=dumpMessage(c_tty,c_yes,'','',c_warning,'tnudlat  is not multiple of 1 hour, ',local_tnudlat ,"F12.1")
     if(ismultiple(local_tnudtop )) iErrNumber=dumpMessage(c_tty,c_yes,'','',c_warning,'tnudtop  is not multiple of 1 hour, ',local_tnudtop ,"F12.1")
     if(ismultiple(local_frqhis  )) iErrNumber=dumpMessage(c_tty,c_yes,'','',c_warning,'frqhis   is not multiple of 1 hour, ',local_frqhis  ,"F12.1")
     if(ismultiple(local_frqanl  )) iErrNumber=dumpMessage(c_tty,c_yes,'','',c_warning,'frqanl   is not multiple of 1 hour, ',local_frqanl  ,"F12.1")
     if(ismultiple(local_radfrq  )) iErrNumber=dumpMessage(c_tty,c_yes,'','',c_warning,'radfrq   is not multiple of 1 hour, ',local_radfrq  ,"F12.1")
     if(ismultiple(local_confrq  )) iErrNumber=dumpMessage(c_tty,c_yes,'','',c_warning,'confrq   is not multiple of 1 hour, ',local_confrq  ,"F12.1")
     if(ismultiple(local_shcufrq )) iErrNumber=dumpMessage(c_tty,c_yes,'','',c_warning,'shcufrq  is not multiple of 1 hour, ',local_shcufrq ,"F12.1")

   endif 

   if(mod(time,local_isan_inc)<dtlongn(ngrid)) last_isan_inc=time
   if(mod(time,local_tnudcent)<dtlongn(ngrid)) last_tnudcent=time
   if(mod(time,local_tnudlat )<dtlongn(ngrid)) last_tnudlat =time
   if(mod(time,local_tnudtop )<dtlongn(ngrid)) last_tnudtop =time
   if(mod(time,local_frqhis  )<dtlongn(ngrid)) last_frqhis  =time
   if(mod(time,local_frqanl  )<dtlongn(ngrid)) last_frqanl  =time
   if(mod(time,local_radfrq  )<dtlongn(ngrid)) last_radfrq  =time
   if(mod(time,local_confrq  )<dtlongn(ngrid)) last_confrq  =time
   if(mod(time,local_shcufrq )<dtlongn(ngrid)) last_shcufrq =time

   nextTime=(/ &
       last_isan_inc+local_isan_inc-futureTime &
      ,last_tnudcent+local_tnudcent-futureTime &
      ,last_tnudlat +local_tnudlat -futureTime &
      ,last_tnudtop +local_tnudtop -futureTime &
      ,last_frqhis  +local_frqhis  -futureTime &
      ,last_frqanl  +local_frqanl  -futureTime &
      ,last_radfrq  +local_radfrq  -futureTime &
      ,last_confrq  +local_confrq  -futureTime &
      ,last_shcufrq +local_shcufrq -futureTime &
      /)



   timeDistance=minval(nextTime)
   if(timeDistance<dtlongn(ngrid)) then 
      if(.not. dtlong_saved) then
         savedDtlongn=dtlongn(ngrid)
         dtlong_saved=.true.
      endif
      dt=timeDistance
   else
      if(dtlong_saved) then 
         dt=savedDtlongn
         dtlong_saved=.false.
      else
         dt=dtlongn(ngrid)
      endif 
   endif 

   futureTime=time+dt

end function nextMainPoint 

end module modSchedUtils 


subroutine modsched(isched, maxsched, maxschent, ngrids, nxtnest, nndtrat, nsubs)

  ! modsched purpose:
  !  (1) schedules nested time step iterations of multiple grids for synchronous
  !      time advance;
  !  (2) schedules data exchange among grids to keep fields consistent among grids
  !
  ! modsched assumptions:
  ! Grids are organized in a tree with:
  !  (1) exactly one node for each grid and
  !  (2) all grids at any nonempty subtree are fully nested at the root grid
  !
  ! Grid nesting geography is restricted by:
  !  (a) There is a single outermost grid that contains all nested grids;
  !  (b) Sibling grids are disjoint (including boundary cells);
  !  (c) If grid b is nested into grid a, then a fixed integer number
  !      of gird b cells partitions a cell of grid a
  !  (d) If grid b is nested into grid a, then the inner part
  !      of grid b (grid b without boundary cells) fully partitions a
  !      set of inner cells of grid a.
  !
  ! Grid nested timestep is restricted by:
  !  (a) A nested grid timestep is an integer fraction of the parent's grid timestep.
  !
  ! modsched reasoning:
  !   A single timestep of any grid requires multiple timesteps of nested grids
  !   (since nested grids have a smaller timestep than coarser grids). Once a grid
  !   runs a timestep, all nested grids should run a set of timesteps to synchronous
  !   advance.
  !   To maintain fields consistent among grids, nested grids should receive
  !   boundary conditions from the parent grid (since the parent -- coarser -- grid is
  !   advanced in time) before running its set of own timesteps. After synchronizing
  !   in time with the parent grid, it should feedback the parent with all detailed
  !   fields to be integrated back into the parent grid (feedback fields).
  !
  ! modsched algorithm, giving the tree and the timestep ratio of each grid
  ! to the parent's grid:
  !   Noting to do on an empty tree;
  !   Visit the tree root as many times as the timestep ratio (for time synchronization
  !   of the root). At each visit, schedule a timestep, send boundary conditions
  !   to all direct sons (to keep fields consistent) and recursivelly run the algorithm
  !   at all sub-trees rooted at sons. Finally, if the root has a parent grid, send feedback
  !   fields to keep parent's field consistent.
  !
  ! modsched data structure:
  !   nsubs:       is the total number of nested timesteps required to advance the full
  !                tree a single timestep;
  !   isched:      what to do at nested timestep i (first index, 1 <= i <= nsubs);
  !   isched(i,1): grid to run this nested timestep
  !   isched(i,2): number of grids to send boundary conditions
  !   isched(i,3): nested timestep counter for this nested timestep grid
  !                during a single timestep of the parent's grid
  !   isched(i,4): number of grids to send feedback fields
  !   isched(i,5): nested timestep counter for this nested timestep grid
  !                during a single timestep of grid 1 (that is, among all
  !                nested timesteps of a single outermost grid timestep)
  !
  ! modsched software assumptions: (not verified by modshched:)
  !   ngrids > 0                       (at least grid 1 exists)
  !   nxtnest(1)=0                     (grid 1 is the single outermost grid)
  !   1 <= nxtnest(2:ngrids) <= ngrids (tree is self contained and is not a forest)
  !   nndtrat(1)=1                     (grid 1 is timestep ratio reference)
  !   nndtrat(2:ngrids) >= 2           (nested grids should have smaller dt);
  !   maxsched >= ngrids               (space allocation)
  !   maxschent >= 5                   (space allocation)

  implicit none
  integer, intent(in)  :: maxsched                   ! maximum number of grids in the tree of nested grids
  integer, intent(in)  :: maxschent                  ! number of features per nested timestep stored at isched
  integer, intent(in)  :: ngrids                     ! number of grids in the tree of nested grids
  integer, intent(in)  :: nxtnest(ngrids)            ! parent grid: tree representation
  integer, intent(in)  :: nndtrat(ngrids)            ! delta t ratio: parent grid/this grid
  integer, intent(out) :: isched(maxsched,maxschent) ! scheduling table (see above)
  integer, intent(out) :: nsubs                      ! total number of nested timesteps for a dtlong

  integer :: nTS(ngrids)      ! number of nested timesteps at each grid

  ! initialize counters

  isched(:,:) = 0
  nsubs = 0
  nTS(:) = 0

  ! visit tree in pre-order to schedule timesteps
  ! (assumes node 1 is the root of a single tree)

  call operTS(1, isched, maxsched, maxschent, ngrids, nxtnest, nndtrat, nsubs, nTS)
end subroutine modsched

! *************************************************************************

recursive subroutine operTS(root, isched, maxsched, maxschent, ngrids, nxtnest, nndtrat, nsubs, nTS)

  ! operTS:
  !   visit the tree of grids with root "root" in preorder;
  !   requires (1 <= root <= ngrids);

  implicit none
  integer, intent(in)    :: root                       ! grid id of the sub-tree root
  integer, intent(in)    :: maxsched                   ! maximum number of grids in the tree of nested grids
  integer, intent(in)    :: maxschent                  ! number of features per timestep stored at isched
  integer, intent(in)    :: ngrids                     ! number of grids in the tree of nested grids
  integer, intent(in)    :: nxtnest(ngrids)            ! parent grid
  integer, intent(in)    :: nndtrat(ngrids)            ! delta t ratio: parent grid/this grid
  integer, intent(inout) :: isched(maxsched,maxschent) ! scheduling table
  integer, intent(inout) :: nsubs                      ! total number of timesteps for a dtlong
  integer, intent(inout) :: nTS(ngrids)                ! timesteps of each grid

  integer :: iter  ! nested timestep count for the current grid
  integer :: node  ! auxiliar grid count
  character(len=10) :: c0, c1
  character(len=*), parameter :: h="**(operTS)**"

  ! should be a non-empty tree

  if (root <= 0 .or. root > ngrids) then
     write(c0,"(i10)") root
     write(c1,"(i10)") ngrids
     write(*,"(a)") h//" INTERNAL ERROR: invoked with root ="//&
          &trim(adjustl(c0))//" outside bounds [1:"//&
          &trim(adjustl(c1))//"]"
     stop
  end if

  ! run all nested timesteps to synchronize

  do iter = 1, nndtrat(root)

     ! check bounds
     ! schedule execution of this grid nested timestep

     nsubs = nsubs + 1
     if (nsubs > maxsched) then
        write(c0,"(i10)") maxsched
        write(*,"(a)") h//" INTERNAL ERROR: increment exausted maxsched ="//&
             &trim(adjustl(c0))
        stop
     end if
     isched(nsubs,1) = root

     ! nested timestep counters relative to next coarser grid and to most coarser grid

     isched(nsubs,3) = iter
     nTS(root) = nTS(root)+1
     isched(nsubs,5) = nTS(root)

     ! # grids to receive boundary conditions
     ! (these are the direct sons)

     isched(nsubs,2) = count(nxtnest(1:ngrids) == root)

     ! schedule execution of trees rooted at each son

     do node = 1, ngrids
        if (nxtnest(node) == root) then
           call operTS(node, isched, maxsched, maxschent, ngrids, nxtnest, nndtrat, nsubs, nTS)
        end if
     end do
  end do

  ! at the adequate nested timestep iteration (after all finer grids are synchronized)
  ! increase count of sending feedback fields to the father, if any.

  if (nxtnest(root) /= 0) then
     isched(nsubs,4) = isched(nsubs,4) + 1
  end if
end subroutine operTS

! *************************************************************************

subroutine dump_modsched(isched, maxsched, maxschent, ngrids, nxtnest, nndtrat, nsubs)
  implicit none
  integer, intent(in) :: maxsched                   ! maximum number of grids in the tree of nested grids
  integer, intent(in) :: maxschent                  ! number of features per timestep stored at isched
  integer, intent(in) :: ngrids                     ! number of grids in the tree of nested grids
  integer, intent(in) :: nxtnest(ngrids)            ! parent grid
  integer, intent(in) :: nndtrat(ngrids)            ! delta t ratio: parent grid/this grid
  integer, intent(in) :: isched(maxsched,maxschent) ! scheduling table (see below)
  integer, intent(in) :: nsubs

  integer :: is                   ! nested timestep count
  integer :: grid                 ! grid at this nested timestep
  integer :: ifg                  ! some finner grid
  integer :: icg                  ! some coarser grid
  integer :: ng                   ! auxiliar grid count
  integer :: bcgrid(ngrids)       ! grids to receive bc
  integer :: dtratouter(ngrids)   ! # grid timesteps to advance dtlong
  character(len=10) :: c0, c1, c2, c3, c4, c5, c6
  character(len=*), parameter :: h="**(dump_modsched)**"

  ! grid timesteps to advance dtlong

  dtratouter = 1
  do grid = 1, ngrids
     icg = grid
     do
        dtratouter(grid) = nndtrat(icg)*dtratouter(grid)
        icg = nxtnest(icg)
        if (icg == 0) exit
     end do
  end do

  ! banner

  open(unit=22,file='brams.log',position='append',action='write')
  write(unit=22,fmt="(a)")
  write(c0,"(i10)") nsubs
  write(unit=22,fmt="(a)") " === Timestep Schedule requires "//trim(adjustl(c0))//&
       &" nested timesteps: ===="


  ! for each timestep:

  do is = 1, nsubs
     grid = isched(is,1)

     ! dump grid to execute timestep and timestep counters

     write(c0,"(i10)") dtratouter(grid)
     write(c1,"(i10)") isched(is,1)
     write(c2,"(i10)") isched(is,3)
     write(c3,"(i10)") isched(is,5)
     write(c4,"(i10)") nndtrat(grid)
     write(c5,"(i10)") nxtnest(grid)
     write(c6,"(i10)") is
     write(unit=22,fmt="(a)", advance="NO") &
          &" nts "//trim(adjustl(c6))//&
          &": Grid "//trim(adjustl(c1))
     if (dtratouter(grid) == 1) then
        write(unit=22,fmt="(a)", advance="NO") " advances dtlong"
     else if (dtratouter(grid) == isched(is,5)) then
        write(unit=22,fmt="(a)", advance="NO") " reaches dtlong"
     else if (isched(is,3) == nndtrat(grid)) then
        write(unit=22,fmt="(a)", advance="NO") " reaches "//&
             &trim(adjustl(c3))//"/"//trim(adjustl(c0))//" of dtlong"//&
             &" (one grid "//trim(adjustl(c5))//" dt)"
     else
        write(unit=22,fmt="(a)", advance="NO") " reaches "//&
             &trim(adjustl(c3))//"/"//trim(adjustl(c0))//" of dtlong"//&
             &" ("//trim(adjustl(c2))//"/"//trim(adjustl(c4))//&
             &" of grid "//trim(adjustl(c5))//" dt)"
     end if

     ! find all direct sons

     ng=0
     do ifg = 1, ngrids
        if (nxtnest(ifg) == grid) then
           ng = ng + 1
           bcgrid(ng) = ifg
        end if
     end do

     ! dump grid(s) to receive boundary conditions

     if (isched(is,2) /= ng) then
        write(c0,"(i10)") is
        write(c1,"(i10)") isched(is,2)
        write(c2,"(i10)") ng
        write(*, "(/,a)") h//"**ERROR**: isched("//trim(adjustl(c0))//&
             &") = "//trim(adjustl(c1))//" and ng = "//trim(adjustl(c2))//&
             &" disagree"
        stop
     else if (isched(is,2) == 1) then
        write(c2,"(i10)") bcgrid(ng)
        write(unit=22,fmt="(a)", advance="NO") "; changes boundary of grid "//trim(adjustl(c2))
     else if (isched(is,2) > 1) then
        write(unit=22,fmt="(a)", advance="NO") "; changes boundary of grids "
        do ifg = 1, ng - 1
           write(c2,"(i10)") bcgrid(ifg)
           write(unit=22,fmt="(a)", advance="NO") trim(adjustl(c2))//", "
        end do
        write(c2,"(i10)") bcgrid(ng)
        write(unit=22,fmt="(a)", advance="NO") trim(adjustl(c2))
     end if

     ! dump grid(s) to receive feedback fields

     if (isched(is,4) /= 0) then
        ifg = isched(is,1)
        icg = nxtnest(ifg)
        write(c4,"(i10)") icg
        write(unit=22,fmt="(a)", advance="NO") "; feedbacks grid "//trim(adjustl(c4))
        do ng = 2, isched(is,4)
           ifg = icg
           icg = nxtnest(ifg)
           write(c4,"(i10)") icg
           write(unit=22,fmt="(a)", advance="NO") " that feedbacks grid "//trim(adjustl(c4))
        end do
     end if
     write(unit=22,fmt="(a)")
  end do
  write(unit=22,fmt="(a)")
  close(unit=22)
end subroutine dump_modsched

!     *****************************************************************

subroutine cfl(n1,n2,n3,i0,j0)

  use mem_basic, only: &
       basic_g           ! intent(in)

  use mem_grid, only: &
       ngrid,         & ! intent(in)
       grid_g           ! intent(in)

  implicit none

  integer, intent(in) :: n1
  integer, intent(in) :: n2
  integer, intent(in) :: n3
  integer, intent(in) :: i0
  integer, intent(in) :: j0

  call cfll(n1,n2,n3,i0,j0  &
       ,basic_g(ngrid)%up   (1,1,1)  ,basic_g(ngrid)%vp   (1,1,1)  &
       ,basic_g(ngrid)%wp   (1,1,1)  ,grid_g(ngrid)%rtgt    (1,1)  &
       ,grid_g(ngrid)%f13t    (1,1)  ,grid_g(ngrid)%f23t    (1,1)  &
       ,grid_g(ngrid)%dxt     (1,1)  ,grid_g(ngrid)%dyt     (1,1)  )
  return
end subroutine cfl

! **********************************************************************

subroutine cfll(n1,n2,n3,i0,j0,up,vp,wp,rtgt,f13t,f23t,dxt,dyt)
  use dump, only: &
    dumpMessage

  use mem_grid, only: &
       nnxp,          & ! intent(in)
       nnyp,          & ! intent(in)
       cflxy,         & ! intent(out)
       cflz,          & ! intent(out)
       ! cfl_max_sum,   &
       jdim,          & ! intent(in)
       ngrids,        & ! intent(in)
       ngrid,         & ! intent(in)
       nxtnest,       & ! intent(in)
       ipm,           & ! intent(in)
       jpm,           & ! intent(in)
       dtlt,          & ! intent(in)
       ht,            & ! intent(in)
       dzt,           & ! intent(in)
       dyncore_flag,  & ! intent(in)
       dtlongn, &
       akminvar, &
       time , &
       grid_g         ! intent(inout)

  use mem_basic, only: basic_g

  use node_mod, only: &
       mynum, &        ! INTENT(IN)
       nMachs

  use mem_turb, only: &
       akmin        ! INTENT(in)

  implicit none
  include "constants.f90"
  character(len=*),parameter :: h='**(cfll)**'
  integer, intent(in) :: n1
  integer, intent(in) :: n2
  integer, intent(in) :: n3
  integer, intent(in) :: i0
  integer, intent(in) :: j0
  real,    intent(in) :: up(n1,n2,n3)
  real,    intent(in) :: vp(n1,n2,n3)
  real,    intent(in) :: wp(n1,n2,n3)
  real,    intent(in) :: rtgt(n2,n3)
  real,    intent(in) :: f13t(n2,n3)
  real,    intent(in) :: f23t(n2,n3)
  real,    intent(in) :: dxt(n2,n3)
  real,    intent(in) :: dyt(n2,n3)


  integer :: i,j,k,ifm,icm,innest
  real :: c1x,c1y,c1z,cflnumh,cflnumv
  real :: cfl_sum
  real :: cfl_max_sum   !MB: should be treated analogous to cflxy and cflz (as ngrid-wide field)
  real :: vctr1(n1)
  real :: vctr2(n1)
  real :: vctr3(n1)
  character(len=5) :: cproc

  !     This routine returns flags the model to bring itself down when the CFL
  !     linear stability criteria on advection is exceeded.
  !     (Actually check on 90% of CFL)

  cflnumh = .90
  cflnumv = .90
  cflxy(ngrid) = 0.
  cflz(ngrid) = 0.
  cfl_max_sum = 0.0


  ! Let's try a new thing... if we have a grid point that is on a
  !   coarse grid, but it is under a nested grid, we will ignore it
  !   under the assumption that the fine grid values will overwrite
  !   it, hence not allowing it to go numerically unstable.

  jloop: do j = 1+jdim,n3-jdim
     iloop: do i = 2,n2-1

        ! See if this is under a fine grid horizontally... ignore vertical for now
        innest=0
        if (ngrids > ngrid) then
           do ifm=ngrid+1,ngrids
              icm=nxtnest(ifm)
              if(icm == ngrid .and. &
                   i+i0 >= ipm(1,ifm) .and. i+i0 <= ipm(nnxp(ifm),ifm) .and. &
                   j+j0 >= jpm(1,ifm) .and. j+j0 <= jpm(nnyp(ifm),ifm) ) then
                 innest=1
                 exit
              endif
           enddo
        endif

        if(innest == 1) then
           cycle iloop
        endif

        kloop: do k = 2,n1-1

           vctr1(k) = .5*(up(k,i,j)+up(k,i-1,j   ))*dtlt*dxt(i,j)
           vctr2(k) = .5*(vp(k,i,j)+vp(k,i,j-jdim))*dtlt*dyt(i,j)
           vctr3(k) = ((wp(k,i,j)+wp(k-1,i,j   ))  &
                      +(up(k,i,j)+up(k,i-1,j   ))*f13t(i,j)*ht(k)*rtgt(i,j)  &
                      +(vp(k,i,j)+vp(k,i,j-jdim))*f23t(i,j)*ht(k)*rtgt(i,j)  &
                      )*.5*dtlt*dzt(k)
        enddo kloop


        if (AKMIN(ngrid) == -1.0) then
            c1x = maxval(abs(vctr1(2:n1-1)))
           c1y = maxval(abs(vctr2(2:n1-1)))
           c1z = maxval(abs(vctr3(2:n1-1)))
           cfl_sum = c1x + c1y + c1z
            akminvar(ngrid)%akmin2d(i,j)=1.0+max(0.,atan(cfl_sum-0.1)**3.)/0.310
           if(cfl_sum>1.6) print*,"cfl_sum=",cfl_sum,max(0.,atan(cfl_sum-0.1)**3.)/0.310,maxval(vctr3(2:n1-1)),mynum;call flush(6)
        endif

        do k = 2,n1-1
           c1x = abs(vctr1(k))
           c1y = abs(vctr2(k))
           c1z = abs(vctr3(k))

           if (c1x .gt. cflxy(ngrid)) cflxy(ngrid) = c1x
           if (c1y .gt. cflxy(ngrid)) cflxy(ngrid) = c1y
           if (c1z .gt. cflz(ngrid)) cflz(ngrid) = c1z

           cfl_sum = c1x + c1y + c1z
            if ( cfl_sum > cfl_max_sum ) then
               cfl_max_sum = cfl_sum
              if ( cfl_max_sum > 1.5 * 1.62 ) then
               write(cproc,fmt='(I5.5)') mynum
               call dump_cfl(cproc,n1,n2,n3,jdim,i0,j0,i,j,k,c1x,c1y,c1z,cfl_max_sum &
                   ,dxt(i,j),dtlt,up(k,i,j),up(k,i-1,j),dyt(i,j),vp(k,i,j),vp(k,i,j-jdim) &
                   , f13t(i,j),f23t(i,j),ht(k),rtgt(i,j) &
                   ,wp(k,i,j),wp(k-1,i,j),dzt(k))
               endif
            end if

        end do
     enddo iloop
  enddo jloop


!Se for tempo de integracao os processadores enviam seu cfl_max pro mestre (proc1)
!Este verifica de algum valor ultrapassa o limite de 75% do máximo permitido e, 
!caso ultrapasse, ajusta o dtlong para valores menores para evitar instabilizacao
if (time>0.0) then 

   call commCFL(cfl_max_sum,ngrid)

endif

  !MB:
  if ( (dyncore_flag==0) .or. (dyncore_flag==1) ) then
    if ( cfl_max_sum > 0.9 * 1.0 ) then
      write(cproc,fmt='(I5.5)') mynum
      call dump_cfl(cproc,n1,n2,n3,jdim,i0,j0,i,j,k,c1x,c1y,c1z,cfl_max_sum &
          ,dxt(i,j),dtlt,up(k,i,j),up(k,i-1,j),dyt(i,j),vp(k,i,j),vp(k,i,j-jdim) &
          , f13t(i,j),f23t(i,j),ht(k),rtgt(i,j) &
          ,wp(k,i,j),wp(k-1,i,j),dzt(k))
      open(unit=66,file='brams.log',status='replace',action='write')
      write(unit=66,fmt='(A)') 'WARNING in cfl1: cfl_max_sum has exceeded 90% of the max. allowed value!'
      write(unit=66,fmt='(A,I8)') "WARNING in cfl1 - mynum:",mynum
      close(unit=66)
    end if
  else if ( dyncore_flag==2 ) then
    if ( cfl_max_sum > 0.9 * 1.62 ) then
      ! criterion for RK3 and upwind 3
      write(cproc,fmt='(I5.5)') mynum
      call dump_cfl(cproc,n1,n2,n3,jdim,i0,j0,i,j,k,c1x,c1y,c1z,cfl_max_sum &
          ,dxt(i,j),dtlt,up(k,i,j),up(k,i-1,j),dyt(i,j),vp(k,i,j),vp(k,i,j-jdim) &
          , f13t(i,j),f23t(i,j),ht(k),rtgt(i,j) &
          ,wp(k,i,j),wp(k-1,i,j),dzt(k)) 
      open(unit=66,file='brams.log',status='replace',action='write')
      write(unit=66,fmt='(A)') 'WARNING in cfl1 RK3: cfl_max_sum has exceeded 90% of the max. allowed value!'
      write(unit=66,fmt='(A,I8)') "WARNING in cfl1 RK3 - mynum:",mynum
      close(unit=66)     
    end if
    if ( cfl_max_sum > 1.5 * 1.62 ) then
      !- criterion for RK3 and upwind 3 (for upwind 5 -> the maximum allowed is 1.4)
      !call fatal_error("SEVERE PROBLEM in cfl1 RK3: cfl_max_sum has exceeded 150% of the max. allowed value!")
      iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
        "SEVERE PROBLEM in cfl1 RK3: cfl_max_sum has exceeded "// &
        " 150% of the max. allowed value!")
    end if
  else if ( dyncore_flag==3 ) then
    if ( cfl_max_sum > 0.9 * 1.0 ) then
      ! criterion for ABM3 and upwind 3
      write(cproc,fmt='(I5.5)') mynum
      call dump_cfl(cproc,n1,n2,n3,jdim,i0,j0,i,j,k,c1x,c1y,c1z,cfl_max_sum &
          ,dxt(i,j),dtlt,up(k,i,j),up(k,i-1,j),dyt(i,j),vp(k,i,j),vp(k,i,j-jdim) &
          , f13t(i,j),f23t(i,j),ht(k),rtgt(i,j) &
          ,wp(k,i,j),wp(k-1,i,j),dzt(k))      
      print*, "WARNING in cfl1 ABM3: cfl_max_sum has exceeded 90% of the max. allowed value!"
    end if
    if ( cfl_max_sum > 1.5 * 1.0 ) then
      !- criterion for ABM3 and upwind 3 (for upwind 5 -> the maximum allowed is ??)
      !call fatal_error("SEVERE PROBLEM in cfl1 ABM3: cfl_max_sum has exceeded 150% of the max. allowed value!")
      write(cproc,fmt='(I5.5)') mynum
      call dump_cfl(cproc,n1,n2,n3,jdim,i0,j0,i,j,k,c1x,c1y,c1z,cfl_max_sum &
          ,dxt(i,j),dtlt,up(k,i,j),up(k,i-1,j),dyt(i,j),vp(k,i,j),vp(k,i,j-jdim) &
          , f13t(i,j),f23t(i,j),ht(k),rtgt(i,j) &
          ,wp(k,i,j),wp(k-1,i,j),dzt(k))      
      iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
        "SEVERE PROBLEM in cfl1 ABM3: cfl_max_sum has exceeded"//&
        " 150% of the max. allowed value!")
    end if
 else
    stop "ERROR in cfl1: false value for dyncore_flag!"
  end if


  return
end subroutine cfll


!=============================================================================================
subroutine commCFL(cfl_max_sum,ngrid)
   !# COmunicate de CFL_max 
   !#
   !# @note
   !#
   !# **Documentation/Version Control**: <documents>
   !#
   !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
   !#
   !# **Date**: 14 October 2021 (Thursday)
   !#
   !# **Full Description**: COmunicate de CFL_max and adjust dtlogn down or up for
   !#                       limits of <15% of CFLMax and >60% CFLMax compute the new
   !#                       courant number
   !#
   !# @endnote
   !#
   !# @warning
   !# ![](https://licensebuttons.net/l/by-sa/3.0/88x31.png"")
   !# Now is under CC Attribution-ShareAlike 4.0 International, please see:
   !# &copy; <https://creativecommons.org/licenses/by-sa/4.0/>
   !# @endwarning
   !#
    
   !Use area
   use dump !Dump contains a lot of functions for debugs and formated printouts
   use ParLib  , only: & !SUbroutines for parallel comunications
    parf_send_noblock_real, &
    parf_get_noblock_real , &
    parf_wait_any_nostatus, &
    parf_wait_all_nostatus

   use node_mod, only: &
       mynum, &        ! INTENT(IN)
       nMachs, &
       nodemxp, &
       nodemyp

   use mem_grid, only: &
       dtlt,          & ! intent(in)
       dyncore_flag,  & ! intent(in)
       dtlongn, &
       ideltat, &
       grid_g, &
       nnxp,           & ! INTENT(IN)
       nnyp,           & ! INTENT(IN)
       nnzp              ! INTENT(IN)

   use ReadBcst, only: Broadcast

   use ref_sounding, only : &
         th01dn,             & ! INTENT(IN)
         pi01dn                ! INTENT(IN)

   use dtset, only: &
       sscourn, &
       maxCflPercent
   
   use modSchedUtils, only: &
       nextMainPoint

   implicit none

   include "constants.f90"
   character(len=*),parameter :: sourceName='modsched.f90' !Name of this source code
   character(len=*),parameter :: procedureName='**commCFL**' !Name of this procedure
   !
   !Local Parameters

   !Input/Output variables
   integer, intent(in) :: ngrid
   real, intent(in) :: cfl_max_sum

   !Local variables
   real :: mcfl(1),cfl_max(nMachs),cfl_max_rec(nMachs)
   integer :: reqsend,recNum,iRecV,reqRecv(nMachs)
   integer :: k,nn2,nn3
   real :: maxCfl
   real :: stepUp,stepDown,tmax,ssmax,dxtmax,ssodx
   real :: aux(200)

   !Code
   if(mynum==1) then
       do iRecv = 1,nmachs-1
          call parf_get_noblock_real(cfl_max_rec(iRecv),1,iRecv,iRecv+10000,reqRecv(iRecv))
          !print *, 'get ',iRecv,iRecv+10000; call flush(6)
       enddo 
   else
       !print *, 'Sou o ',mynum,' enviando ',cfl_max_sum; call flush()
       call parf_send_noblock_real((/cfl_max_sum/),1,0,mynum-1+10000,reqsend)
       !print *,mynum,mynum-1+10000; call flush(6)
   endif 
   if(mynum==1) then
       cfl_max(1)=cfl_max_sum
       do iRecV = 1,nmachs-1
         !receive_acou the message recNum from request
         call parf_wait_any_nostatus(nmachs-1,reqRecv,recNum)
         !print *,' Sou o mestre e recebi de ',recNum+1,' com valor ',cfl_max_rec(recNum)
         cfl_max(recNum+1)=cfl_max_rec(recNum)
       enddo
       maxCfl=maxval(cfl_max)
       
       !Determina valores para subir ou descer o dtlong
       if(dtlongn(ngrid)>99) then
         stepup=5.0
         stepDown=10.0
       endif 
       if(dtlongn(ngrid)<100 .and. dtlongn(ngrid)>10) then 
         stepup=1.0
         stepDown=5.0
      endif
      if(dtlongn(ngrid)<=10) then
         stepup=0.0
         stepDown=1.0
      endif

       !Determina o CFL máximo
       maxCflPercent=maxCfl/(0.9 * 1.62)


       !Conforme o quão próximo do estouro está o CFL ajusta o timestep
       if(ideltat>0) then
         if(maxCflPercent>1.00) dtlongn(ngrid)=dtlongn(ngrid)-(StepDown*3) !call adjustDtlongDown(dtlt,nGrid)
         if(maxCflPercent>0.75) dtlongn(ngrid)=dtlongn(ngrid)-(StepDown*2) !call adjustDtlongDown(dtlt,nGrid)
         if(maxCflPercent>0.60) dtlongn(ngrid)=dtlongn(ngrid)-(StepDown) !call adjustDtlongDown(dtlt,nGrid)
       endif
       if(ideltat>1) then
         if(maxCflPercent<0.15) dtlongn(ngrid)=dtlongn(ngrid)+stepUp !call adjustDtlongUp(dtlt,nGrid)
       endif
        
        !Determina o novo courant
        nn2 = nodemxp(mynum,nGrid)
        nn3 = nodemyp(mynum,nGrid)
        do k = 1,nnzp(nGrid)
          aux(k) = th01dn(k,1) * pi01dn(k,1) / c_cp
        enddo
        tmax = maxval(aux(1:nnzp(ngrid)))
        ssmax = sqrt(c_cp / c_cv * c_rgas * tmax)

        dxtmax = max(grid_g(nGrid)%dxt(1,1), &
                  grid_g(nGrid)%dxt(nn2,1),       &
                  grid_g(nGrid)%dxt(nn2,nn3),     &
                  grid_g(nGrid)%dxt(1,nn3)        )
        ssodx = ssmax * dxtmax
        sscourn(nGrid) = 2. * ssodx * dtlongn(nGrid)

        dtlongn(ngrid)=nextMainPoint(ngrid)

   endif 

   call Broadcast(dtlongn(ngrid),0,'dtlong')

end subroutine commCFL 


! ***************************************************************************

subroutine dump_dtset(nndtflg)
  use mem_grid, only: &
       ngrids,        &
       nndtrat,       &
       nnacoust,      &
       dtlongn

  implicit none
  include "constants.f90"
  integer, intent(in) :: nndtflg
  integer :: ifm

  open(unit=22,file='brams.log',position='append',action='write')
  write(unit=22,fmt="(a)")
  if (nndtflg == 0) then
     write(unit=22,fmt="(a)") " === Grids Delta T ===="
  else
     write(unit=22,fmt="(a)") " === Changed Grids Delta T ===="
  end if
  write(unit=22,fmt="(a)") " Grid, delta t, fraction of next coarser grid dt, acoustic steps per delta t"
  do ifm = 1, ngrids
     write(unit=22,fmt="(i5,f9.2,15x,i4,24x,i4)") ifm, dtlongn(ifm), nndtrat(ifm), nnacoust(ifm)
  end do
  close(unit=22)
end subroutine dump_dtset

subroutine dump_cfl(cproc,n1,n2,n3,jdim,i0,j0,i,j,k,c1x,c1y,c1z,cfl_max_sum &
                   ,dxt,dtlt,up,upm1,dyt,vp,vpm1,f13t,f23t,ht,rtgt &
                   ,wp,wpm1,dzt)

   use mem_grid, only: grid_g,nnxp,nnyp,nnzp
   use mem_basic, only: basic_g


   integer, intent(in) :: n1,n2,n3,jdim,i0,j0,i,j,k
   real, intent(in)    :: cfl_max_sum,c1x,c1y,c1z,wp,wpm1,dzt
   real, intent(in)    :: dxt,dtlt,up,upm1,dyt,vp,vpm1, f13t,f23t,ht,rtgt
   character(len=*), intent(in) :: cproc

   open(unit=22,file='blow.'//cproc//'.log',position='append',action='write')

   write(22,fmt='("n1,n2,n3,nnxp,nnyp,nnzp=",6(I6.6,1X))') n1,n2,n3,nnzp(1),nnxp(1),nnyp(1)
   write(22,fmt='("Position of error point at local  grid: i,j,k:",3(I5.5,1X))') i,j,k
   write(22,fmt='("Position of error point at global grid: i,j  :",2(I5.5,1X))') i+i0,j+j0
   write(22,fmt='("Limits Lat: ",3(F8.2,1X))') grid_g(1)%glat(1,1),grid_g(1)%glat(n2,n3),grid_g(1)%glat(1,1)-grid_g(1)%glat(1,2)
   write(22,fmt='("Limits Lon: ",3(F8.2,1X))') grid_g(1)%glon(1,1),grid_g(1)%glon(n2,n3),grid_g(1)%glon(1,1)-grid_g(1)%glon(2,1)
   write(22,fmt='("Lat,lon of error point: ",2(F8.2,1X))') grid_g(1)%glat(i,j),grid_g(1)%glon(i,j)
   write(22,*) 'Diagnostics: ============================================================================='
   write(22,fmt='(3(A,E20.5,1X))') 'U: ',c1x,'V: ',c1y,'W: ',c1z
   write(22,fmt='("cfl_max_sum: ",E20.5,", Max allowed is: ",E20.5,", Percent: ",F8.3,"%")') cfl_max_sum,1.5 * 1.62,(cfl_max_sum/(1.5 * 1.62))*100.
   write(22,fmt='("U: dxt(i,j),dtlt,up(k,i,j),up(k,i-1,j): ",4(E20.5,1X))') dxt,dtlt,up,upm1
   write(22,fmt='("V: dyt(i,j),vp(k,i,j),vp(k,i,j-jdim)",4(E20.5,1X))') dyt,dtlt,vp,vpm1
   write(22,fmt='("W: f13t(i,j),f23t(i,j),ht(k),rtgt(i,j),wp(k,i,j),wp(k-1,i,j): ",4(E20.5,1X))') f13t,f23t,ht,rtgt
   write(22,fmt='("W: wp(k,i,j),wp(k-1,i,j): ",2(E20.5,1X))') wp,wpm1
   write(22,fmt='("dzt(k): ",E20.5)') dzt 
   write(22,fmt='("W U_component: ",E20.5)') (up+upm1)*f13t*ht*rtgt
   write(22,fmt='("W V_component: ",E20.5)') (vp+vpm1)*f23t*ht*rtgt
   write(22,fmt='("W W_component: ",E20.5)') wp+wpm1
   write(22,fmt='("W Product_factor: ",E20.5)') .5*dtlt*dzt
   write(22,fmt='("PI0,Theta,thp,th0: ",4(E20.5,1X))') basic_g(1)%pi0(k,i,j),basic_g(1)%theta(k,i,j),basic_g(1)%thp(k,i,j),basic_g(1)%th0(k,i,j)
   write(22,*) "###########################################################################################################"
   close(unit=22)             
end subroutine dump_cfl

