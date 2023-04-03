!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine fcorio(n2,n3,fcoru,fcorv,glat)
!> @brief: This routine calculates the Coriolis parameter.
!! @author:  unknow
!! @date:  17/Nov/2015
!! @details: modified by Luiz Flavio in order to prevent pointer in calls  
!! @version:  5.2
!! @param: n2,n3,fcoru,fcorv,glat
!! 
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
!!
  use mem_grid, only: jdim
  use rconstants

  implicit none
  
  integer,intent(in) :: n2,n3
  real,intent(in)    :: glat(n2,n3)
  real,intent(inout) :: fcoru(n2,n3)
  real,intent(inout) :: fcorv(n2,n3)

  real :: omega2
  integer :: i,j
  
  omega2 = 2. * omega

  do j = 1,max(1,n3-1)
     do i = 1,n2-1
        fcoru(i,j) = omega2*sin((glat(i,j)+glat(i+1,j))  &
             *.5*pi180)
        fcorv(i,j) = omega2*sin((glat(i,j)+glat(i,j+jdim))  &
             *.5*pi180)
     enddo
  enddo

  !--(DMK)--------------------------------------
  call CorrectBorderFcorx(1,n2,n3,fcoru)
  call CorrectBorderFcory(1,n2,n3,fcoru)

  if (jdim==1) then
     call CorrectBorderFcorx(1,n2,n3,fcorv)
     call CorrectBorderFcory(1,n2,n3,fcorv)
  end if
  !--(DMK)--------------------------------------
end subroutine fcorio

!     ******************************************************************

!subroutine corlos(mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv, ut_ptr, vt_ptr)
subroutine corlos(mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv, ut, vt)
!> @brief: This routine is the coriolis driver.  Its purpose is to compute
!!coriolis accelerations for u and v and add them into
!!the accumulated tendency arrays of ut_ptr and vt_ptr.
!! @author:  unknow
!! @date:  17/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
!! 
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
!!
use mem_tend, only: &
      tend
use mem_grid, only: &
      icorflg
use mem_scratch, only: &
    scratch
implicit none

integer,intent(in) :: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv

!real, pointer :: ut_ptr(:), vt_ptr(:)
integer :: i,j,k,n
real, intent(inout) :: ut(mzp,mxp,myp)
real, intent(inout) :: vt(mzp,mxp,myp)
real :: vt3da(mzp,mxp,myp)

if(icorflg.eq.0) return

n=0
do j=1,myp
  do i=1,mxp
    do k=1,mzp
      n=n+1
      !ut(k,i,j)=tend%ut(n)  !MB: (original code)
      !vt(k,i,j)=tend%vt(n)
      !ut(k,i,j) = ut_ptr(n)
      !vt(k,i,j) = vt_ptr(n)
      vt3da(k,i,j)=scratch%scr1(n)   !MB: really necessary ???
    end do
  end do
end do

call corlsu(mzp,mxp,myp,i0,j0,ia,izu,ja,jz,ut,vt3da)

call corlsv(mzp,mxp,myp,i0,j0,ia,iz,ja,jzv,vt,vt3da)

n=0
do j=1,myp
  do i=1,mxp
    do k=1,mzp
      n=n+1
      !tend%ut(n)=ut(k,i,j)   !MB: (original code)
      !tend%vt(n)=vt(k,i,j)
      !ut_ptr(n) = ut(k,i,j) 
      !vt_ptr(n) = vt(k,i,j)
      scratch%scr1(n)=vt3da(k,i,j)  !MB: really necessary ???
    end do
  end do
end do

end

!     ******************************************************************

subroutine corlsu(m1,m2,m3,i0,j0,ia,iz,ja,jz,ut,vt3da)
!> @brief: This routine is the coriolis tendencies U direction
!! @author:  unknow
!! @date:  17/Nov/2015
!! @version:  5.2
!! @param: m1,m2,m3,i0,j0,ia,iz,ja,jz
!! 
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
!!
use mem_grid
use rconstants
use mem_scratch
use ref_sounding
use mem_basic, only: &
     basic_g   

implicit none

integer,intent(in) :: m1,m2,m3,i0,j0,ia,iz,ja,jz
real,intent(inout) :: ut(m1,m2,m3),vt3da(m1,m2,m3)

integer :: i,j,k
real :: c1

do j=ja,jz
   do i=ia,iz
      do k=2,m1-1
         vt3da(k,i,j)=(basic_g(ngrid)%vc(k,i,j)+basic_g(ngrid)%vc(k,i,j-jdim)  &
              +basic_g(ngrid)%vc(k,i+1,j)+basic_g(ngrid)%vc(k,i+1,j-jdim))*.25
      enddo
   enddo
enddo

c1=1./(erad*erad*2.)
if(ihtran.eq.0) c1=0.
do j=ja,jz
   do i=ia,iz
      do k=2,m1-1
         ut(k,i,j)=ut(k,i,j)-vt3da(k,i,j)*(-basic_g(ngrid)%fcoru(i,j)  &
                  +c1*(vt3da(k,i,j)*xm(i+i0)-basic_g(ngrid)%uc(k,i,j)*yt(j+j0)))
      enddo
   enddo
enddo
!
if (initial /= 1) return

if (if_adap == 0 .and. itopo == 1) then

   do j = ja,jz
      do i = ia,iz
         do k = 1,m1
            vctr2(k) = zt(k) * grid_g(ngrid)%rtgu(i,j) + grid_g(ngrid)%topu(i,j)
         enddo
         call htint(nzp,v01dn(1,ngrid),zt,nz,vctr5,vctr2)
         do k = 2,m1-1
            ut(k,i,j) = ut(k,i,j) - basic_g(ngrid)%fcoru(i,j) * vctr5(k)
         enddo
      enddo
   enddo

else

   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-1
            ut(k,i,j) = ut(k,i,j) - basic_g(ngrid)%fcoru(i,j) * v01dn(k,ngrid)
         enddo
      enddo
   enddo

endif

return
end

! **************************************************************

subroutine corlsv(m1,m2,m3,i0,j0,ia,iz,ja,jz,vt,vt3da)
!> @brief: This routine is the coriolis tendencies V direction
!! @author:  unknow
!! @date:  17/Nov/2015
!! @version:  5.2
!! @param: m1,m2,m3,i0,j0,ia,iz,ja,jz
!! 
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
!!
use mem_grid
use rconstants
use mem_scratch
use ref_sounding
use mem_basic, only: &
     basic_g   
implicit none

integer,intent(in) :: m1,m2,m3,i0,j0,ia,iz,ja,jz
real,intent(inout) :: vt(m1,m2,m3),vt3da(m1,m2,m3)

integer :: i,j,k
real :: c1

do j = ja,jz
   do i = ia,iz
      do k = 2,m1-1
         vt3da(k,i,j) = (basic_g(ngrid)%uc(k,i,j) + basic_g(ngrid)%uc(k,i-1,j)  &
            + basic_g(ngrid)%uc(k,i,j+jdim) + basic_g(ngrid)%uc(k,i-1,j+jdim)) * .25
      enddo
   enddo
enddo

c1 = 1. / (erad * erad * 2.)
if (ihtran .eq. 0) c1 = 0.
do j = ja,jz
   do i = ia,iz
      do k = 2,m1-1
         vt(k,i,j) = vt(k,i,j) - vt3da(k,i,j) * (basic_g(ngrid)%fcorv(i,j)  &
            - c1 * (basic_g(ngrid)%vc(k,i,j) * xt(i+i0) - vt3da(k,i,j) * ym(j+j0)))
      enddo
   enddo
enddo

if (initial /= 1) return

if (if_adap == 0 .and. itopo == 1) then
   do j = ja,jz
      do i = ia,iz
         do k = 1,m1
            vctr2(k) = zt(k) * grid_g(ngrid)%rtgv(i,j) + grid_g(ngrid)%topv(i,j)
         enddo
         call htint(nzp,u01dn(1,ngrid),zt,nz,vctr5,vctr2)
         do k = 2,m1-1
            vt(k,i,j) = vt(k,i,j) + basic_g(ngrid)%fcorv(i,j) * vctr5(k)
         enddo
      enddo
   enddo

else

   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-1
            vt(k,i,j) = vt(k,i,j) + basic_g(ngrid)%fcorv(i,j) * u01dn(k,ngrid)
         enddo
      enddo
   enddo

endif

return
end




!--(DMK)-----------------------------------------------------------------------

subroutine CorrectBorderFcorx(grid, n2, n3, fcor)
!> @brief: CorrectBorderFcorx
!! @author:  unknow
!! @date:  17/Nov/2015
!! @version:  5.2
!! @param: grid, n2, n3, fcor
!! 
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
!!
  use ParLib, only: &
       parf_get_noblock, &
       parf_pack, &
       parf_send_block, &
       parf_wait_nostatus, &
       parf_unpack

  use mem_grid, only: &
       nnyp
  use node_mod, only: &
       mynum, nmachs, &
       !--(DMK)-------------------------------
       nxbeg, nxend, nybeg, nyend, &
       ixb, ixe, iyb, iye, &
       nodeibcon, nodei0, nodej0, &
       !--(DMK)-------------------------------
       f_ndmd_size,       &
       machs
  use ModBuffering, only: &
       FieldSection2Buffer, &
       Buffer2FieldSection

  implicit none
  integer, intent(in)  :: grid
  integer, intent(in)  :: n2
  integer, intent(in)  :: n3
  real,    intent(inout) :: fcor(n2,n3)

  include "i8.h"
  integer(i8) :: cnt_i8
  integer :: StartRecvJ, EndRecvJ
  integer :: StartRecvI, EndRecvI
  integer :: StartSendJ, EndSendJ
  integer :: StartSendI, EndSendI
  integer :: StartOwnJ, EndOwnJ
  integer :: StartOwnI, EndOwnI
  integer :: sendProc, iSend, nSend
  integer :: recvProc, iRecv, nRecv
  integer :: maxSize
  integer :: cnt, cnt2
  integer :: ierr
  integer, parameter :: tag=33
  integer :: nodeRecv(nmachs)
  integer :: sizeRecv(nmachs)
  integer :: firstIRecv(nmachs)
  integer :: lastIRecv(nmachs)
  integer :: firstJRecv(nmachs)
  integer :: lastJRecv(nmachs)
  integer :: reqRecv(nmachs)
  integer :: nodeSend(nmachs)
  integer :: sizeSend(nmachs)
  integer :: firstISend(nmachs)
  integer :: lastISend(nmachs)
  integer :: firstJSend(nmachs)
  integer :: lastJSend(nmachs)
  real, allocatable :: buffRecv(:,:)
  real, allocatable :: buffUnpack(:)
  real, allocatable :: buffSend(:)
  real, allocatable :: buffToPack(:)

  character(len=8) :: c0, c1, c2, c3, c4, c5, c6, c7, c8, c9
  character(len=*), parameter :: h="**(CorrectBorderFcorx)**"
  logical, parameter :: dumpLocal=.false.

  ! send and receive counts

  nSend=0; nRecv=0

  ! This process fcor border is wrong if
  ! high x bound is not a full domain border.

  ! Find region to be corrected, to be received from other processes;
  ! mark as [StartRecvI:EndRecvI,StartRecvJ:EndRecvJ] (global indices).

  ! Find which processes will send data to correct fcor;
  ! Store at receive info data structure (local indices)

  if (.not. btest(nodeibcon(mynum,grid),1)) then

     ! portion of global domain to be received

     StartRecvI = n2 + nodei0(mynum,grid)
     EndRecvI   = StartRecvI
     StartRecvJ =  1 + nodej0(mynum,grid)
     EndRecvJ   = n3 + nodej0(mynum,grid)
     if (dumpLocal) then
        write(c0,"(i8)") mynum
        write(c1,"(i8)") StartRecvI
        write(c2,"(i8)") EndRecvI
        write(c3,"(i8)") StartRecvJ
        write(c4,"(i8)") EndRecvJ
        write(*,"(a)") h//" proc "//trim(adjustl(c0))//&
             " has to correct fcor["//&
             trim(adjustl(c1))//":"//trim(adjustl(c2))//","//&
             trim(adjustl(c3))//":"//trim(adjustl(c4))//"]"
     end if

     ! processes that will send data to correct fcor at this process 

     do sendProc = 1, nmachs
        if (sendProc == mynum) cycle
        StartSendI = max(ixb(sendProc, grid), StartRecvI)
        EndSendI   = min(ixe(sendProc, grid), EndRecvI)
        StartSendJ = max(iyb(sendProc, grid), StartRecvJ)
        EndSendJ   = min(iye(sendProc, grid), EndRecvJ)
        if (StartSendI > EndSendI .or. StartSendJ > EndSendJ) cycle

        if (StartSendJ == 2) then
           StartSendJ = 1
        end if
        if (EndSendJ == nnyp(grid)-1) then
           EndSendJ = nnyp(grid)
        end if

        ! build receives from sendProc to this proc, including
        ! positions to store data (local indices) at this proc

        nRecv = nRecv + 1
        nodeRecv(nRecv) = machs(sendProc)
        sizeRecv(nRecv) = (EndSendI-StartSendI+1)*(EndSendJ-StartSendJ+1)
        firstIRecv(nRecv) = StartSendI - nodei0(mynum,grid)
        lastIRecv(nRecv)  = EndSendI   - nodei0(mynum,grid)
        firstJRecv(nRecv) = StartSendJ - nodej0(mynum,grid)
        lastJRecv(nRecv)  = EndSendJ   - nodej0(mynum,grid)

        if (dumpLocal) then
           write(c0,"(i8)") sendProc
           write(c1,"(i8)") StartSendI
           write(c2,"(i8)") EndSendI
           write(c3,"(i8)") StartSendJ
           write(c4,"(i8)") EndSendJ
           write(c5,"(i8)") mynum
           write(c6,"(i8)") firstIRecv(nRecv)
           write(c7,"(i8)") lastIRecv(nRecv)
           write(c8,"(i8)") firstJRecv(nRecv)
           write(c9,"(i8)") lastJRecv(nRecv)
           write(*,"(a)") h//" proc "//trim(adjustl(c5))//&
                " will issue recv from proc "//trim(adjustl(c0))//" fcor["//&
                trim(adjustl(c1))//":"//trim(adjustl(c2))//","//&
                trim(adjustl(c3))//":"//trim(adjustl(c4))//"] to store at ["//&
                trim(adjustl(c6))//":"//trim(adjustl(c7))//","//&
                trim(adjustl(c8))//":"//trim(adjustl(c9))//"]"
        end if
     end do
  end if

  ! This process may send the portion of the global domain that it owns,
  ! including global domain boundaries, to correct fcor at other processes. 
  ! Owned points are (global indices) [StartOwnI:EndOwnI,StartOwnJ:EndOwnJ].
  ! It owns internal points + ghost zone points that are full domain boundaries.
  ! Ties on full domain boundaries (arise on more than one process) are given
  ! to processes where inner area "extends" to tied point.

  if (btest(nodeibcon(mynum,grid),0)) then
     StartOwnI = 1 + nodei0(mynum,grid)
  else
     StartOwnI = ixb(mynum,grid)
  end if
  if (btest(nodeibcon(mynum,grid),1)) then
     EndOwnI = n2 + nodei0(mynum,grid)
  else
     EndOwnI = ixe(mynum,grid)
  end if
  if (btest(nodeibcon(mynum,grid),2)) then
     StartOwnJ = 1 + nodej0(mynum,grid)
  else
     StartOwnJ = iyb(mynum,grid)
  end if
  if (btest(nodeibcon(mynum,grid),3)) then
     EndOwnJ = n3 + nodej0(mynum,grid)
  else
     EndOwnJ = iye(mynum,grid)
  end if
  
  ! Find processes that will receive fcor points from this process.
  ! These are processes where some of the fcor points that are wrong
  ! belong to this process (and so, are correct).

  ! A process fcor border has to be corrected iff its
  ! high x bound is not a full domain border (global indices)

  do recvProc = 1, nmachs
     if (recvProc == mynum .or. btest(nodeibcon(recvProc,grid),1)) cycle
     StartRecvI = max(nxend(recvProc, grid), StartOwnI)
     EndRecvI   = min(nxend(recvProc, grid), EndOwnI)
     StartRecvJ = max(nybeg(recvProc, grid), StartOwnJ)
     EndRecvJ   = min(nyend(recvProc, grid), EndOwnJ)
     if (StartRecvI > EndRecvI .or. StartRecvJ > EndRecvJ) cycle

     ! build sends from this process to recvProc, and store
     ! at send data structure (local indices)

     nSend = nSend + 1
     nodeSend(nSend) = machs(recvProc)
     sizeSend(nSend) = (EndRecvI-StartRecvI+1)*(EndRecvJ-StartRecvJ+1)
     firstISend(nSend) = StartRecvI - nodei0(mynum,grid)
     lastISend(nSend)  = EndRecvI   - nodei0(mynum,grid)
     firstJSend(nSend) = StartRecvJ - nodej0(mynum,grid)
     lastJSend(nSend)  = EndRecvJ   - nodej0(mynum,grid)
     
     if (dumpLocal) then
        write(c0,"(i8)") recvProc
        write(c1,"(i8)") StartRecvI
        write(c2,"(i8)") EndRecvI
        write(c3,"(i8)") StartRecvJ
        write(c4,"(i8)") EndRecvJ
        write(c5,"(i8)") mynum
        write(c6,"(i8)") firstISend(nSend)
        write(c7,"(i8)") lastISend(nSend)
        write(c8,"(i8)") firstJSend(nSend)
        write(c9,"(i8)") lastJSend(nSend)
        write(*,"(a)") h//" proc "//trim(adjustl(c5))//&
             " will send to proc "//trim(adjustl(c0))//" fcor["//&
             trim(adjustl(c1))//":"//trim(adjustl(c2))//","//&
             trim(adjustl(c3))//":"//trim(adjustl(c4))//"] retrieved from ["//&
             trim(adjustl(c6))//":"//trim(adjustl(c7))//","//&
             trim(adjustl(c8))//":"//trim(adjustl(c9))//"]"
     end if
  end do

  ! allocate receive and unpack buffers

  if (nRecv > 0) then
     maxSize = maxval(sizeRecv(1:nRecv))*f_ndmd_size
     allocate(buffRecv(maxSize,nRecv), stat=ierr)
     if (ierr /= 0) then
        write(c0,"(i8)") maxSize
        write(c1,"(i8)") nRecv
        write(c2,"(i8)") ierr
        call fatal_error(h//" allocate buffRecv("//&
             trim(adjustl(c0))//","//trim(adjustl(c1))//&
             ") fails with stat="//trim(adjustl(c2)))
     end if
     allocate(buffUnpack(maxSize), stat=ierr)
     if (ierr /= 0) then
        write(c0,"(i8)") maxSize
        write(c2,"(i8)") ierr
        call fatal_error(h//" allocate buffUnpack("//&
             trim(adjustl(c0))//") fails with stat="//trim(adjustl(c2)))
     end if
  end if

  ! allocate send and pack buffers

  if (nSend > 0) then
     maxSize = maxval(sizeSend(1:nSend))*f_ndmd_size
     allocate(buffSend(maxSize), stat=ierr)
     if (ierr /= 0) then
        write(c0,"(i8)") maxSize
        write(c2,"(i8)") ierr
        call fatal_error(h//" allocate buffSend("//&
             trim(adjustl(c0))//") fails with stat="//&
             trim(adjustl(c2)))
     end if
     allocate(buffToPack(maxSize), stat=ierr)
     if (ierr /= 0) then
        write(c0,"(i8)") maxSize
        write(c2,"(i8)") ierr
        call fatal_error(h//" allocate buffToPack("//&
             trim(adjustl(c0))//") fails with stat="//&
             trim(adjustl(c2)))
     end if
  end if

  ! post non-blocking receives to avoid deadlock

  do iRecv = 1, nRecv
     call parf_get_noblock(buffRecv(1,iRecv), &
          int(sizeRecv(iRecv),i8)*int(f_ndmd_size,i8), &
          nodeRecv(iRecv), &
          tag+grid, &
          reqRecv(iRecv))
  end do

  ! issue blocking sends (receives are already posted by target process)

  do iSend = 1, nSend

     ! gather info to send

     cnt=0
     call FieldSection2Buffer(fcor, 2, &
          firstISend(iSend), lastISend(iSend), &
          firstJSend(iSend), lastJSend(iSend), &
          buffToPack, cnt)

     ! pack gathered info

     cnt_i8 = 0_i8
     call parf_pack(buffToPack, int(sizeSend(iSend),i8), &
          buffSend, int(sizeSend(iSend),i8)*int(f_ndmd_size,i8), cnt_i8)

     ! send

     call parf_send_block(buffSend, int(sizeSend(iSend),i8)*int(f_ndmd_size,i8), &
          nodeSend(iSend), tag+grid)
  end do

  ! wait for all receives

  do iRecv = 1, nRecv
     call parf_wait_nostatus(reqRecv(iRecv))
  end do

  ! unpack received info 

  do iRecv = 1, nRecv

     ! unpack buffer

     cnt_i8=0
     call parf_unpack(buffRecv(:,iRecv), &
          (int(sizeRecv(iRecv),i8)*int(f_ndmd_size,i8)), &
          cnt_i8, buffUnpack, int(sizeRecv(iRecv),i8))

     ! store at desired positions of fcor

     cnt=0
     call Buffer2FieldSection(fcor, 2, &
          firstIRecv(iRecv), lastIRecv(iRecv), &
          firstJRecv(iRecv), lastJRecv(iRecv), &
          buffUnpack, cnt)
  end do

  ! deallocate receive buffers

  if (nRecv > 0) then
     deallocate(buffRecv, stat=ierr)
     if (ierr /= 0) then
        write(c2,"(i8)") ierr
        call fatal_error(h//" deallocate buffRecv "//&
             "fails with stat="//trim(adjustl(c2)))
     end if
     deallocate(buffUnpack, stat=ierr)
     if (ierr /= 0) then
        write(c2,"(i8)") ierr
        call fatal_error(h//" deallocate buffUnpack "//&
             "fails with stat="//trim(adjustl(c2)))
     end if
  end if

  ! deallocate send buffer

  if (nSend > 0) then
     deallocate(buffSend, stat=ierr)
     if (ierr /= 0) then
        write(c2,"(i8)") ierr
        call fatal_error(h//" deallocate buffSend "//&
             "fails with stat="//trim(adjustl(c2)))
     end if
     deallocate(buffToPack, stat=ierr)
     if (ierr /= 0) then
        write(c2,"(i8)") ierr
        call fatal_error(h//" deallocate buffToPack "//&
             "fails with stat="//trim(adjustl(c2)))
     end if
  end if
end subroutine CorrectBorderFcorx

!--(DMK)-----------------------------------------------------------------------




!--(DMK)-----------------------------------------------------------------------

subroutine CorrectBorderFcory(grid, n2, n3, fcor)
!> @brief: CorrectBorderFcory
!! @author:  unknow
!! @date:  17/Nov/2015
!! @version:  5.2
!! @param: grid, n2, n3, fcor
!! 
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
!!
  use ParLib, only: &
       parf_get_noblock, &
       parf_pack, &
       parf_send_block, &
       parf_wait_nostatus, &
       parf_unpack
  use mem_grid, only: &
       nnxp
  use node_mod, only: &
       mynum, nmachs, &
       !--(DMK)-------------------------------
       nxbeg, nxend, nybeg, nyend, &
       ixb, ixe, iyb, iye, &
       nodeibcon, nodei0, nodej0, &
       !--(DMK)-------------------------------
       f_ndmd_size,       &
       machs
  use ModBuffering, only: &
       FieldSection2Buffer, &
       Buffer2FieldSection

  implicit none
  integer, intent(in)  :: grid
  integer, intent(in)  :: n2
  integer, intent(in)  :: n3
  real,    intent(inout) :: fcor(n2,n3)

  include "i8.h"
  integer(i8) :: cnt_i8
  integer :: StartRecvJ, EndRecvJ
  integer :: StartRecvI, EndRecvI
  integer :: StartSendJ, EndSendJ
  integer :: StartSendI, EndSendI
  integer :: StartOwnJ, EndOwnJ
  integer :: StartOwnI, EndOwnI
  integer :: sendProc, iSend, nSend
  integer :: recvProc, iRecv, nRecv
  integer :: maxSize
  integer :: cnt, cnt2
  integer :: ierr
  integer, parameter :: tag=33
  integer :: nodeRecv(nmachs)
  integer :: sizeRecv(nmachs)
  integer :: firstIRecv(nmachs)
  integer :: lastIRecv(nmachs)
  integer :: firstJRecv(nmachs)
  integer :: lastJRecv(nmachs)
  integer :: reqRecv(nmachs)
  integer :: nodeSend(nmachs)
  integer :: sizeSend(nmachs)
  integer :: firstISend(nmachs)
  integer :: lastISend(nmachs)
  integer :: firstJSend(nmachs)
  integer :: lastJSend(nmachs)
  real, allocatable :: buffRecv(:,:)
  real, allocatable :: buffUnpack(:)
  real, allocatable :: buffSend(:)
  real, allocatable :: buffToPack(:)

  character(len=8) :: c0, c1, c2, c3, c4, c5, c6, c7, c8, c9
  character(len=*), parameter :: h="**(CorrectBorderFcory)**"
  logical, parameter :: dumpLocal=.false.

  ! send and receive counts

  nSend=0; nRecv=0

  ! This process fcor border is wrong if
  ! high y bound is not a full domain border.

  ! Find region to be corrected, to be received from other processes;
  ! mark as [StartRecvI:EndRecvI,StartRecvJ:EndRecvJ] (global indices).

  ! Find which processes will send data to correct fcor;
  ! Store at receive info data structure (local indices)

  if (.not. btest(nodeibcon(mynum,grid),3)) then

     ! portion of global domain to be received

     StartRecvI = 1  + nodei0(mynum,grid)
     EndRecvI   = n2 + nodei0(mynum,grid)
     StartRecvJ = n3 + nodej0(mynum,grid)
     EndRecvJ   = n3 + nodej0(mynum,grid)
     if (dumpLocal) then
        write(c0,"(i8)") mynum
        write(c1,"(i8)") StartRecvI
        write(c2,"(i8)") EndRecvI
        write(c3,"(i8)") StartRecvJ
        write(c4,"(i8)") EndRecvJ
        write(*,"(a)") h//" proc "//trim(adjustl(c0))//&
             " has to correct fcor["//&
             trim(adjustl(c1))//":"//trim(adjustl(c2))//","//&
             trim(adjustl(c3))//":"//trim(adjustl(c4))//"]"
     end if

     ! processes that will send data to correct fcor at this process 

     do sendProc = 1, nmachs
        if (sendProc == mynum) cycle
        StartSendI = max(ixb(sendProc, grid), StartRecvI)
        EndSendI   = min(ixe(sendProc, grid), EndRecvI)
        StartSendJ = max(iyb(sendProc, grid), StartRecvJ)
        EndSendJ   = min(iye(sendProc, grid), EndRecvJ)
        if (StartSendI > EndSendI .or. StartSendJ > EndSendJ) cycle

        if (StartSendI == 2) then
           StartSendI = 1
        end if
        if (EndSendI == nnxp(grid)-1) then
           EndSendI = nnxp(grid)
        end if

        ! build receives from sendProc to this proc, including
        ! positions to store data (local indices) at this proc

        nRecv = nRecv + 1
        nodeRecv(nRecv) = machs(sendProc)
        sizeRecv(nRecv) = (EndSendI-StartSendI+1)*(EndSendJ-StartSendJ+1)
        firstIRecv(nRecv) = StartSendI - nodei0(mynum,grid)
        lastIRecv(nRecv)  = EndSendI   - nodei0(mynum,grid)
        firstJRecv(nRecv) = StartSendJ - nodej0(mynum,grid)
        lastJRecv(nRecv)  = EndSendJ   - nodej0(mynum,grid)

        if (dumpLocal) then
           write(c0,"(i8)") sendProc
           write(c1,"(i8)") StartSendI
           write(c2,"(i8)") EndSendI
           write(c3,"(i8)") StartSendJ
           write(c4,"(i8)") EndSendJ
           write(c5,"(i8)") mynum
           write(c6,"(i8)") firstIRecv(nRecv)
           write(c7,"(i8)") lastIRecv(nRecv)
           write(c8,"(i8)") firstJRecv(nRecv)
           write(c9,"(i8)") lastJRecv(nRecv)
           write(*,"(a)") h//" proc "//trim(adjustl(c5))//&
                " will issue recv from proc "//trim(adjustl(c0))//" fcor["//&
                trim(adjustl(c1))//":"//trim(adjustl(c2))//","//&
                trim(adjustl(c3))//":"//trim(adjustl(c4))//"] to store at ["//&
                trim(adjustl(c6))//":"//trim(adjustl(c7))//","//&
                trim(adjustl(c8))//":"//trim(adjustl(c9))//"]"
        end if
     end do
  end if

  ! This process may send the portion of the global domain that it owns,
  ! including global domain boundaries, to correct fcor at other processes. 
  ! Owned points are (global indices) [StartOwnI:EndOwnI,StartOwnJ:EndOwnJ].
  ! It owns internal points + ghost zone points that are full domain boundaries.
  ! Ties on full domain boundaries (arise on more than one process) are given
  ! to processes where inner area "extends" to tied point.

  if (btest(nodeibcon(mynum,grid),0)) then
     StartOwnI = 1 + nodei0(mynum,grid)
  else
     StartOwnI = ixb(mynum,grid)
  end if
  if (btest(nodeibcon(mynum,grid),1)) then
     EndOwnI = n2 + nodei0(mynum,grid)
  else
     EndOwnI = ixe(mynum,grid)
  end if
  if (btest(nodeibcon(mynum,grid),2)) then
     StartOwnJ = 1 + nodej0(mynum,grid)
  else
     StartOwnJ = iyb(mynum,grid)
  end if
  if (btest(nodeibcon(mynum,grid),3)) then
     EndOwnJ = n3 + nodej0(mynum,grid)
  else
     EndOwnJ = iye(mynum,grid)
  end if
  
  ! Find processes that will receive fcor points from this process.
  ! These are processes where some of the fcor points that are wrong
  ! belong to this process (and so, are correct).

  ! A process fcor border has to be corrected iff its
  ! high x bound is not a full domain border (global indices)

  do recvProc = 1, nmachs
     if (recvProc == mynum .or. btest(nodeibcon(recvProc,grid),3)) cycle
     StartRecvI = max(nxbeg(recvProc, grid), StartOwnI)
     EndRecvI   = min(nxend(recvProc, grid), EndOwnI)
     StartRecvJ = max(nyend(recvProc, grid), StartOwnJ)
     EndRecvJ   = min(nyend(recvProc, grid), EndOwnJ)
     if (StartRecvI > EndRecvI .or. StartRecvJ > EndRecvJ) cycle

     ! build sends from this process to recvProc, and store
     ! at send data structure (local indices)

     nSend = nSend + 1
     nodeSend(nSend) = machs(recvProc)
     sizeSend(nSend) = (EndRecvI-StartRecvI+1)*(EndRecvJ-StartRecvJ+1)
     firstISend(nSend) = StartRecvI - nodei0(mynum,grid)
     lastISend(nSend)  = EndRecvI   - nodei0(mynum,grid)
     firstJSend(nSend) = StartRecvJ - nodej0(mynum,grid)
     lastJSend(nSend)  = EndRecvJ   - nodej0(mynum,grid)
     
     if (dumpLocal) then
        write(c0,"(i8)") recvProc
        write(c1,"(i8)") StartRecvI
        write(c2,"(i8)") EndRecvI
        write(c3,"(i8)") StartRecvJ
        write(c4,"(i8)") EndRecvJ
        write(c5,"(i8)") mynum
        write(c6,"(i8)") firstISend(nSend)
        write(c7,"(i8)") lastISend(nSend)
        write(c8,"(i8)") firstJSend(nSend)
        write(c9,"(i8)") lastJSend(nSend)
        write(*,"(a)") h//" proc "//trim(adjustl(c5))//&
             " will send to proc "//trim(adjustl(c0))//" fcor["//&
             trim(adjustl(c1))//":"//trim(adjustl(c2))//","//&
             trim(adjustl(c3))//":"//trim(adjustl(c4))//"] retrieved from ["//&
             trim(adjustl(c6))//":"//trim(adjustl(c7))//","//&
             trim(adjustl(c8))//":"//trim(adjustl(c9))//"]"
     end if
  end do

  ! allocate receive and unpack buffers

  if (nRecv > 0) then
     maxSize = maxval(sizeRecv(1:nRecv))*f_ndmd_size
     allocate(buffRecv(maxSize,nRecv), stat=ierr)
     if (ierr /= 0) then
        write(c0,"(i8)") maxSize
        write(c1,"(i8)") nRecv
        write(c2,"(i8)") ierr
        call fatal_error(h//" allocate buffRecv("//&
             trim(adjustl(c0))//","//trim(adjustl(c1))//&
             ") fails with stat="//trim(adjustl(c2)))
     end if
     allocate(buffUnpack(maxSize), stat=ierr)
     if (ierr /= 0) then
        write(c0,"(i8)") maxSize
        write(c2,"(i8)") ierr
        call fatal_error(h//" allocate buffUnpack("//&
             trim(adjustl(c0))//") fails with stat="//trim(adjustl(c2)))
     end if
  end if

  ! allocate send and pack buffers

  if (nSend > 0) then
     maxSize = maxval(sizeSend(1:nSend))*f_ndmd_size
     allocate(buffSend(maxSize), stat=ierr)
     if (ierr /= 0) then
        write(c0,"(i8)") maxSize
        write(c2,"(i8)") ierr
        call fatal_error(h//" allocate buffSend("//&
             trim(adjustl(c0))//") fails with stat="//&
             trim(adjustl(c2)))
     end if
     allocate(buffToPack(maxSize), stat=ierr)
     if (ierr /= 0) then
        write(c0,"(i8)") maxSize
        write(c2,"(i8)") ierr
        call fatal_error(h//" allocate buffToPack("//&
             trim(adjustl(c0))//") fails with stat="//&
             trim(adjustl(c2)))
     end if
  end if

  ! post non-blocking receives to avoid deadlock

  do iRecv = 1, nRecv
     call parf_get_noblock(buffRecv(1,iRecv), &
          int(sizeRecv(iRecv),i8)*int(f_ndmd_size,i8), &
          nodeRecv(iRecv), &
          tag+grid, &
          reqRecv(iRecv))
  end do

  ! issue blocking sends (receives are already posted by target process)

  do iSend = 1, nSend

     ! gather info to send

     cnt=0
     call FieldSection2Buffer(fcor, 2, &
          firstISend(iSend), lastISend(iSend), &
          firstJSend(iSend), lastJSend(iSend), &
          buffToPack, cnt)

     ! pack gathered info

     cnt_i8 = 0
     call parf_pack(buffToPack, int(sizeSend(iSend),i8), &
          buffSend, int(sizeSend(iSend),i8)*int(f_ndmd_size,i8), cnt_i8)

     ! send

     call parf_send_block(buffSend, int(sizeSend(iSend),i8)*int(f_ndmd_size,i8), &
          nodeSend(iSend), tag+grid)
  end do

  ! wait for all receives

  do iRecv = 1, nRecv
     call parf_wait_nostatus(reqRecv(iRecv))
  end do

  ! unpack received info 

  do iRecv = 1, nRecv

     ! unpack buffer

     cnt_i8=0
     call parf_unpack(buffRecv(:,iRecv), &
          (int(sizeRecv(iRecv),i8)*int(f_ndmd_size,i8)), &
          cnt_i8, buffUnpack, int(sizeRecv(iRecv),i8))

     ! store at desired positions of fcor

     cnt=0
     call Buffer2FieldSection(fcor, 2, &
          firstIRecv(iRecv), lastIRecv(iRecv), &
          firstJRecv(iRecv), lastJRecv(iRecv), &
          buffUnpack, cnt)
  end do

  ! deallocate receive buffers

  if (nRecv > 0) then
     deallocate(buffRecv, stat=ierr)
     if (ierr /= 0) then
        write(c2,"(i8)") ierr
        call fatal_error(h//" deallocate buffRecv "//&
             "fails with stat="//trim(adjustl(c2)))
     end if
     deallocate(buffUnpack, stat=ierr)
     if (ierr /= 0) then
        write(c2,"(i8)") ierr
        call fatal_error(h//" deallocate buffUnpack "//&
             "fails with stat="//trim(adjustl(c2)))
     end if
  end if

  ! deallocate send buffer

  if (nSend > 0) then
     deallocate(buffSend, stat=ierr)
     if (ierr /= 0) then
        write(c2,"(i8)") ierr
        call fatal_error(h//" deallocate buffSend "//&
             "fails with stat="//trim(adjustl(c2)))
     end if
     deallocate(buffToPack, stat=ierr)
     if (ierr /= 0) then
        write(c2,"(i8)") ierr
        call fatal_error(h//" deallocate buffToPack "//&
             "fails with stat="//trim(adjustl(c2)))
     end if
  end if
end subroutine CorrectBorderFcory

!--(DMK)-----------------------------------------------------------------------
