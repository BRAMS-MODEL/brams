!
! Copyright (C) 1991-2004  ; All Rights Reserved ; Colorado State University
! Colorado State University Research Foundation ; ATMET, LLC
!
! This file is free software; you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later version.
!
! This software is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with this
! program; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!===============================================================================
!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################

subroutine node_sendnbc(ifm, icm)

  use ParLib, only: &
       parf_get_noblock, &
       parf_pack, &
!!$       parf_pack_int, &
!!$       parf_pack_real, &
       parf_send_noblock

  use mem_grid, only: &
       NGRIDS,        & ! INTENT(IN)
       NNZP             ! INTENT(IN)
  use node_mod, only: &
       NMACHS,        & ! INTENT(IN)
       IRECV_REQ,     & ! INTENT(IN)
       IGET_PATHS,    & ! INTENT(IN)
       NODE_BUFFS,    & ! INTENT(OUT)
       F_NDMD_SIZE,   & ! INTENT(IN)
       MACHS,         & ! INTENT(IN)
       IPATHS,        & ! INTENT(IN)
       MYNUM,         & ! INTENT(IN)
       ISEND_REQ,     & ! INTENT(INOUT)
       nodemxp,       & ! INTENT(IN)
       nodemyp,       & ! INTENT(IN)
       nodei0,        & ! INTENT(IN)
       nodej0,        & ! INTENT(IN)
       mynum
  use var_tables, only: &
       NUM_SCALAR,      & ! INTENT(IN)
       SCALAR_TAB         ! INTENT(IN)
  use mem_basic, only:  &
       BASIC_G            ! INTENT(IN)

  implicit none
  include "i8.h"
  ! Arguments:
  integer, intent(in) :: ifm, icm
  ! Local Variables:
  integer(i8) :: ipos, nsize
  integer :: nm,i1,i2,j1,j2,k1,k2,ng,itype,mtp,iptr,nv
!!$  real, allocatable, save :: buffnest(:)
!!$  integer, save :: ncall=0, membuff,membuff_extra,nvar
  real, allocatable :: buffnest(:)
  integer :: membuff, nvar

  itype=5

  !______________________
  !
  !   First, before we send anything, let's post the receives.

  do nm=1,nmachs
     irecv_req(nm)=0
     if (iget_paths(itype,ifm,nm).ne.0) then
        call parf_get_noblock(node_buffs(nm)%lbc_recv_buff, &
             int(node_buffs(nm)%nrecv*f_ndmd_size,i8), machs(nm), 5000+icm, &
             irecv_req(nm))
     endif
  enddo

  ! Send coarse grid points necessary for fine grid boundary interpolation
  !   to fine grid nodes. Note that even though coarse grid points are sent,
  !   ipaths is referenced by the fine grid, since all nests only have one
  !   parent, not vice versa.


  ! Compute size of buffer needed and allocate if necessary
!!$  if(ncall == 0) then
!!$     ncall=1
!!$     membuff_extra=nvar*2+100
     membuff=0
     do ng=1,ngrids
        do nm=1,nmachs
           if(ipaths(1,itype,ng,nm)/=0) then
              i1=ipaths(1,itype,ng,nm)
              i2=ipaths(2,itype,ng,nm)
              j1=ipaths(3,itype,ng,nm)
              j2=ipaths(4,itype,ng,nm)
              k1=1
              k2=nnzp(ng)
              nvar=4 + num_scalar(ng)
              mtp=(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
              membuff=max(membuff, mtp*nvar)
           endif
        enddo
     enddo
     membuff=membuff + nvar*2 + 100
     allocate (buffnest(membuff))
!!$     print*,'sending nesting condition: Allocate buffer for:',mynum &
!!$              ,membuff,nvar
!!$  endif


  do nm=1,nmachs

     isend_req(nm)=0

     if(ipaths(1,itype,ifm,nm).ne.0) then

        i1=ipaths(1,itype,ifm,nm)
        i2=ipaths(2,itype,ifm,nm)
        j1=ipaths(3,itype,ifm,nm)
        j2=ipaths(4,itype,ifm,nm)
        k1=1
        k2=nnzp(icm)

        mtp=(i2-i1+1)*(j2-j1+1)*(k2-k1+1)

  ! Put variables into buffer. All need coarse grid density weighting first.

        iptr=0
        call mknest_buff(1,basic_g(icm)%uc(1,1,1),buffnest(1+iptr)  &
            ,basic_g(icm)%dn0(1,1,1),nnzp(icm),nodemxp(mynum,icm),nodemyp(mynum,icm)  &
            ,nodei0(mynum,icm),nodej0(mynum,icm),i1,i2,j1,j2,k1,k2,mynum,nm,nv)
        iptr=iptr+mtp
        call mknest_buff(2,basic_g(icm)%vc(1,1,1),buffnest(1+iptr)  &
            ,basic_g(icm)%dn0(1,1,1),nnzp(icm),nodemxp(mynum,icm),nodemyp(mynum,icm)  &
            ,nodei0(mynum,icm),nodej0(mynum,icm),i1,i2,j1,j2,k1,k2,mynum,nm,nv)
        iptr=iptr+mtp
        call mknest_buff(3,basic_g(icm)%wc(1,1,1),buffnest(1+iptr)  &
            ,basic_g(icm)%dn0(1,1,1),nnzp(icm),nodemxp(mynum,icm),nodemyp(mynum,icm)  &
            ,nodei0(mynum,icm),nodej0(mynum,icm),i1,i2,j1,j2,k1,k2,mynum,nm,nv)
        iptr=iptr+mtp
        call mknest_buff(4,basic_g(icm)%pc(1,1,1),buffnest(1+iptr)  &
            ,basic_g(icm)%dn0(1,1,1),nnzp(icm),nodemxp(mynum,icm),nodemyp(mynum,icm)  &
            ,nodei0(mynum,icm),nodej0(mynum,icm),i1,i2,j1,j2,k1,k2,mynum,nm,nv)
        iptr=iptr+mtp

        do nv=1,num_scalar(ifm)
           call mknest_buff(5,scalar_tab(nv,icm)%var_p,buffnest(1+iptr)  &
               ,basic_g(icm)%dn0(1,1,1),nnzp(icm),nodemxp(mynum,icm),nodemyp(mynum,icm)  &
               ,nodei0(mynum,icm),nodej0(mynum,icm),i1,i2,j1,j2,k1,k2,mynum,nm,nv)
           iptr=iptr+mtp
        enddo


        ipos = 0
        nsize = node_buffs(nm)%nsend*f_ndmd_size
        call parf_pack(i1, node_buffs(nm)%lbc_send_buff, nsize, ipos)
        call parf_pack(i2, node_buffs(nm)%lbc_send_buff, nsize, ipos)
        call parf_pack(j1, node_buffs(nm)%lbc_send_buff, nsize, ipos)
        call parf_pack(j2, node_buffs(nm)%lbc_send_buff, nsize, ipos)
        call parf_pack(k1, node_buffs(nm)%lbc_send_buff, nsize, ipos)
        call parf_pack(k2, node_buffs(nm)%lbc_send_buff, nsize, ipos)
        call parf_pack(mynum, node_buffs(nm)%lbc_send_buff, nsize, ipos)
        call parf_pack(nvar, node_buffs(nm)%lbc_send_buff, nsize, ipos)
        call parf_pack(iptr, node_buffs(nm)%lbc_send_buff, nsize, ipos)

        call parf_pack(buffnest, int(iptr,i8), &
             node_buffs(nm)%lbc_send_buff, &
             nsize, ipos)

        call parf_send_noblock(node_buffs(nm)%lbc_send_buff, &
             ipos, ipaths(5,itype,ifm,nm), 5000+icm, isend_req(nm))
     endif

  enddo

  if (allocated(buffnest)) deallocate(buffnest)

end subroutine node_sendnbc
!
!     ****************************************************************
!
subroutine mknest_buff(ivarn, ac, acs, den, m1, m2, m3, i0, j0,  &
     i1, i2, j1, j2, k1, k2, mynum, nm, nv)

  implicit none
  ! Arguments:
  integer, intent(in) :: ivarn, m1, m2, m3, i0, j0, i1, i2, j1, j2, k1, k2, &
       mynum, nm, nv
  real, intent(in)    :: ac(m1,m2,m3), den(m1,m2,m3)
  real, intent(inout) :: acs(0:k2-k1,0:i2-i1,0:j2-j1)
  ! Local Variables:
  integer :: i, j, k

  !     ivarn = variable types 1- u
  !                            2- v
  !                            3- w
  !                            4- p
  !                            5- scalar

  if(ivarn.eq.5) then
     do j=j1,j2
        do i=i1,i2
           do k=k1,k2
              acs(k-k1,i-i1,j-j1)=ac(k,i-i0,j-j0)*den(k,i-i0,j-j0)
           enddo
        enddo
     enddo
  elseif(ivarn.eq.1) then
     do j=j1,j2
        do i=i1,i2
           do k=k1,k2
              acs(k-k1,i-i1,j-j1)=ac(k,i-i0,j-j0)  &
                   *((den(k,i-i0,j-j0)+den(k,i+1-i0,j-j0))*.5)
           enddo
        enddo
     enddo
  elseif(ivarn.eq.2) then
     do j=j1,j2
        do i=i1,i2
           do k=k1,k2
              acs(k-k1,i-i1,j-j1)=ac(k,i-i0,j-j0)  &
                   *((den(k,i-i0,j-j0)+den(k,i-i0,j+1-j0))*.5)
           enddo
        enddo
     enddo
  elseif(ivarn.eq.3) then
     do j=j1,j2
        do i=i1,i2
           do k=k1,k2-1
              acs(k-k1,i-i1,j-j1)=ac(k,i-i0,j-j0)  &
                   *((den(k,i-i0,j-j0)+den(k+1,i-i0,j-j0))*.5)
           enddo
        enddo
     enddo
  elseif(ivarn.eq.4) then
     do j=j1,j2
        do i=i1,i2
           do k=k1,k2
              acs(k-k1,i-i1,j-j1)=ac(k,i-i0,j-j0)
           enddo
        enddo
     enddo
  endif

end subroutine mknest_buff
!
!     ****************************************************************
!
subroutine node_getnbc(ifm, icm)

  use dump, only: &
    dumpMessage

  use ParLib, only: &
       parf_wait_nostatus, &
       parf_unpack
!!$       parf_unpack_int, &
!!$       parf_unpack_real

  use mem_grid, only: &
       MAXNZP,        & ! INTENT(IN)
       MAXNXP,        & ! INTENT(IN)
       MAXNYP,        & ! INTENT(IN)
       NNZP             ! INTENT(IN)
  use node_mod, only: &
       NMACHS,        & ! INTENT(IN)
       IPATHS,        & ! INTENT(IN)
       ISEND_REQ,     & ! INTENT(INOUT)
       IGET_PATHS,    & ! INTENT(IN)
       IRECV_REQ,     & ! INTENT(INOUT)
       NBUFF_NEST,    & ! INTENT(IN)
       NODE_BUFFS,    & ! INTENT(OUT)
       F_NDMD_SIZE,   & ! INTENT(IN)
       MYNUM,         & ! INTENT(IN)
       nodemxp,       & ! INTENT(IN)
       nodemyp,       & ! INTENT(IN)
       nodei0,        & ! INTENT(IN)
       nodej0,        & ! INTENT(IN)
       nodeibcon,     & ! INTENT(IN)
       mynum
  use grid_dims, only: &
       maxmach             ! INTENT(IN)
  use mem_basic, only: &
       BASIC_G             ! INTENT(IN)
  use mem_nestb, only: &
       NBOUNDS           ! INTENT(IN)

  implicit none
  include "i8.h"
  include "constants.f90"
  ! Arguments:
  integer, intent(in) :: ifm, icm
  ! Local Variables:
  integer(i8) :: ipos, nsize
  integer :: IntArr(1)
  integer, dimension(maxmach) :: i1c,i2c,j1c,j2c,k1c,k2c,iptv,iptc
  integer :: itype,nm,iptr,machf,nvar,nwords,nv,nxc,nyc,nzc,mtp
!!$  real, allocatable, save :: buffnest(:)
!!$  integer, save :: ncall=0, nbuff_save=0
  real, allocatable :: buffnest(:)

  integer :: ierr
  real, allocatable :: scr1(:), scr2(:)
  character(len=*), parameter :: h="**(node_getnbc)**"
  itype      = 5

  !_____________________________________________________________________
  !
  !  First, let's make sure our sends are all finished and de-allocated
  do nm=1,nmachs
     if(ipaths(1,itype,ifm,nm).ne.0)then
        call parf_wait_nostatus(isend_req(nm))
     endif
  enddo
  !_____________________________________________________________________
  !
  !  Now, let's wait on our receives

  do nm=1,nmachs
     if(iget_paths(itype,ifm,nm).ne.0)then
        call parf_wait_nostatus(irecv_req(nm))
     endif
  enddo
  !_____________________________________________________________________
  !


  ! Compute size of buffer needed and allocate if necessary

!!$  if(nbuff_nest > nbuff_save) then
!!$     nbuff_save = nbuff_nest
  allocate (buffnest(nbuff_nest))
!!$  endif

  !     From the fine grid nodes, get the coarse grid buffers,
  !      interpolate the boundaries, and put them in the "b" array.

  iptr=0
  do nm=1,nmachs

     if(iget_paths(itype,ifm,nm).ne.0) then
        ipos = 0
        nsize=node_buffs(nm)%nrecv*f_ndmd_size
        call parf_unpack(node_buffs(nm)%lbc_recv_buff, &
             nsize, ipos, i1c(nm:), 1_i8)
        call parf_unpack(node_buffs(nm)%lbc_recv_buff, &
             nsize, ipos, i2c(nm:), 1_i8)
        call parf_unpack(node_buffs(nm)%lbc_recv_buff, &
             nsize, ipos, j1c(nm:), 1_i8)
        call parf_unpack(node_buffs(nm)%lbc_recv_buff, &
             nsize, ipos, j2c(nm:), 1_i8)
        call parf_unpack(node_buffs(nm)%lbc_recv_buff, &
             nsize, ipos, k1c(nm:), 1_i8)
        call parf_unpack(node_buffs(nm)%lbc_recv_buff, &
             nsize, ipos, k2c(nm:), 1_i8)
        call parf_unpack(node_buffs(nm)%lbc_recv_buff, &
             nsize, ipos, IntArr, 1_i8); machf=IntArr(1)
        call parf_unpack(node_buffs(nm)%lbc_recv_buff, &
             nsize, ipos, IntArr, 1_i8); nvar=IntArr(1)
        call parf_unpack(node_buffs(nm)%lbc_recv_buff, &
             nsize, ipos, IntArr, 1_i8); nwords=IntArr(1)

        call parf_unpack(node_buffs(nm)%lbc_recv_buff, &
             nsize, ipos, buffnest(1+iptr:), int(nwords,i8))

        iptc(nm)=1+iptr
        iptv(nm)=0

        iptr=iptr+nwords
     endif
  enddo

  ! scratch area

  allocate(scr1(maxnzp*maxnxp*maxnyp), stat=ierr)
  if (ierr /= 0) &!call fatal_error(h//" allocating scr1")
    iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
    "ERROR allocating scr1")

  allocate(scr2(maxnzp*maxnxp*maxnyp), stat=ierr)
  if (ierr /= 0) &!call fatal_error(h//" allocating scr2")
    iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
    "ERROR allocating scr2")

  !   We have all the coarse grid info. Start looping through each variable.

  do nv=1,nvar

  !            First, construct coarse grid variable in scr1.

     scr1 = 0.
     scr2 = 0.
     do nm=1,nmachs
        if(iget_paths(itype,ifm,nm).ne.0) then
           nzc=k2c(nm)-k1c(nm)+1
           nxc=i2c(nm)-i1c(nm)+1
           nyc=j2c(nm)-j1c(nm)+1
           mtp=nzc*nxc*nyc
           call unmkbuff(scr1(1),buffnest(iptc(nm)+iptv(nm))  &
                ,maxnzp,maxnxp,maxnyp,nzc,nxc,nyc  &
                ,i1c(nm),i2c(nm),j1c(nm),j2c(nm),k1c(nm),k2c(nm),mynum)

           iptv(nm)=iptv(nm)+mtp
        endif
     enddo

  !            Do the actual interpolation and put stuff into the "b" array

     if(nv.eq.1) then
        call par_bintp(scr1(1),scr2(1)  &
             ,basic_g(ifm)%dn0(1,1,1)  &
             ,maxnzp,maxnxp,maxnyp,nnzp(ifm)  &
             ,nnzp(ifm),nodemxp(mynum,ifm),nodemyp(mynum,ifm)  &
             ,ifm,1,nodei0(mynum,ifm),nodej0(mynum,ifm),nodeibcon(mynum,ifm)  &
             ,nbounds(ifm)%bux(1,1,1),nbounds(ifm)%buy(1,1,1)  &
             ,nbounds(ifm)%buz(1,1,1),mynum)
     elseif(nv.eq.2) then
        call par_bintp(scr1(1),scr2(1)  &
             ,basic_g(ifm)%dn0(1,1,1)  &
             ,maxnzp,maxnxp,maxnyp,nnzp(ifm)  &
             ,nnzp(ifm),nodemxp(mynum,ifm),nodemyp(mynum,ifm)  &
             ,ifm,2,nodei0(mynum,ifm),nodej0(mynum,ifm),nodeibcon(mynum,ifm)  &
             ,nbounds(ifm)%bvx(1,1,1),nbounds(ifm)%bvy(1,1,1)  &
             ,nbounds(ifm)%bvz(1,1,1),mynum)
     elseif(nv.eq.3) then
        call par_bintp(scr1(1),scr2(1)  &
             ,basic_g(ifm)%dn0(1,1,1)  &
             ,maxnzp,maxnxp,maxnyp,nnzp(ifm)  &
             ,nnzp(ifm),nodemxp(mynum,ifm),nodemyp(mynum,ifm)  &
             ,ifm,3,nodei0(mynum,ifm),nodej0(mynum,ifm),nodeibcon(mynum,ifm)  &
             ,nbounds(ifm)%bwx(1,1,1),nbounds(ifm)%bwy(1,1,1)  &
             ,nbounds(ifm)%bwz(1,1,1),mynum)
     elseif(nv.eq.4) then
        call par_bintp(scr1(1),scr2(1)  &
             ,basic_g(ifm)%dn0(1,1,1)  &
             ,maxnzp,maxnxp,maxnyp,nnzp(ifm)  &
             ,nnzp(ifm),nodemxp(mynum,ifm),nodemyp(mynum,ifm)  &
             ,ifm,4,nodei0(mynum,ifm),nodej0(mynum,ifm),nodeibcon(mynum,ifm)  &
             ,nbounds(ifm)%bpx(1,1,1),nbounds(ifm)%bpy(1,1,1)  &
             ,nbounds(ifm)%bpz(1,1,1),mynum)
     else
        call par_bintp(scr1(1),scr2(1)  &
             ,basic_g(ifm)%dn0(1,1,1)  &
             ,maxnzp,maxnxp,maxnyp,nnzp(ifm)  &
             ,nnzp(ifm),nodemxp(mynum,ifm),nodemyp(mynum,ifm)  &
             ,ifm,5,nodei0(mynum,ifm),nodej0(mynum,ifm),nodeibcon(mynum,ifm)  &
             ,nbounds(ifm)%bsx(1,1,1,nv-4),nbounds(ifm)%bsy(1,1,1,nv-4)  &
             ,nbounds(ifm)%bsz(1,1,1,nv-4),mynum)
     endif

  enddo

  deallocate(buffnest)
  deallocate(scr1)
  deallocate(scr2)
end subroutine node_getnbc
!
!     ****************************************************************
!
subroutine unmkbuff(ac,buff,max1,max2,max3,m1,m2,m3,  &
     i1,i2,j1,j2,k1,k2,mynum)
  implicit none
  ! Arguments:
  integer, intent(in) :: max1, max2, max3, m1, m2, m3, i1, i2, j1, j2, k1, k2,&
       mynum
  real, intent(in)    :: buff(0:m1-1,0:m2-1,0:m3-1)
  real, intent(inout) :: ac(max1,max2,max3)
  ! Local Variables:
  integer :: i, j, k

  do j=j1,j2
     do i=i1,i2
        do k=k1,k2
           ac(k,i,j)=buff(k-k1,i-i1,j-j1)
        enddo
     enddo
  enddo

end subroutine unmkbuff
!!$!
!!$!     ****************************************************************
!!$!
!!$subroutine prtlev(ac,m1,m2,m3,lev,my)
!!$  implicit none
!!$  ! Arguments:
!!$  integer, intent(in) :: m1, m2, m3, lev, my
!!$  real, intent(in)    :: ac(m1,m2,m3)
!!$  ! Local variables:
!!$  integer :: j, i
!!$
!!$  print*,'======',my,'====================',m1,m2,m3
!!$  do j=m3,10,-1
!!$     print '(62f6.0)',(ac(lev,i,j),i=10,m2)
!!$  enddo
!!$end subroutine prtlev
