!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine nstbdriv()

  use mem_tend, only: &
       tend

  use var_tables, only: &
       num_scalar,      &
       scalar_tab

  use mem_basic, only: &
       basic_g

  use mem_nestb, only: &
       nbounds

  use node_mod, only: &
       mxp,           &
       myp,           &
       mzp,           &
       ia,            &
       iz,            &
       ja,            &
       jz,            &
       ibcon,mynum

  use mem_grid, only: &
       dtlv,          &
       dtlt,          &
       nndtrat,       &
       isstp,         &
       ngrid,         &
       nstbot,        &
       nsttop,        &
       jdim

  implicit none


  integer :: n,k
  real :: tymeinvv,tymeinvs

  tymeinvv = 1.0 / (dtlv * 0.5 * float(nndtrat(ngrid)+2-isstp))
  tymeinvs = 1.0 / (dtlt * float(nndtrat(ngrid)+1-isstp))
  
  write (*,fmt='("wt-nest-in: ",I2.2,1X,I1.1,1X,4(E9.3,1X),"/",4(E9.3,1X))') mynum,1,(tend%wt(k),k=1,4),(tend%wt(k),k=1,4)
  call flush(6)
  
  call nstbtnd(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
       ,basic_g(ngrid)%up(1,1,1),tend%ut(1)  &
       ,nbounds(ngrid)%bux(1,1,1),nbounds(ngrid)%buy(1,1,1)  &
       ,nbounds(ngrid)%buz(1,1,1)  &
       ,'u',tymeinvv,nstbot,nsttop,jdim)

  call nstbtnd(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
       ,basic_g(ngrid)%vp(1,1,1),tend%vt(1)  &
       ,nbounds(ngrid)%bvx(1,1,1),nbounds(ngrid)%bvy(1,1,1)  &
       ,nbounds(ngrid)%bvz(1,1,1)  &
       ,'v',tymeinvv,nstbot,nsttop,jdim)

  call nstbtnd(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
       ,basic_g(ngrid)%wp(1,1,1),tend%wt(1)  &
       ,nbounds(ngrid)%bwx(1,1,1),nbounds(ngrid)%bwy(1,1,1)  &
       ,nbounds(ngrid)%bwz(1,1,1)  &
       ,'w',tymeinvv,nstbot,nsttop,jdim)

  call nstbtnd(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
       ,basic_g(ngrid)%pp(1,1,1),tend%pt(1)  &
       ,nbounds(ngrid)%bpx(1,1,1),nbounds(ngrid)%bpy(1,1,1)  &
       ,nbounds(ngrid)%bpz(1,1,1)  &
       ,'p',tymeinvv,nstbot,nsttop,jdim)

  do n = 1,num_scalar(ngrid)
     call nstbtnd(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
          ,scalar_tab(n,ngrid)%var_p,scalar_tab(n,ngrid)%var_t  &
          ,nbounds(ngrid)%bsx(1,1,1,n),nbounds(ngrid)%bsy(1,1,1,n)  &
          ,nbounds(ngrid)%bsz(1,1,1,n)  &
          ,'t',tymeinvs,nstbot,nsttop,jdim)
  enddo
    write (*,fmt='("wt-nest-ou: ",I2.2,1X,I1.1,1X,4(E9.3,1X),"/",4(E9.3,1X))') mynum,1,(tend%wt(k),k=1,4),(tend%wt(k),k=1,4)
  call flush(6)
end subroutine nstbdriv

!******************************************************************************
!******************************************************************************

!from initlz:
!      call prgintrp(nnzp(icm),nnxp(icm),nnyp(icm),nnzp(icm),nnxp(icm)
!     +   ,nnyp(icm),0,0,ifm,1,b)
!from model:
!      call prgintrp(nnzp(icm),nnxp(icm),nnyp(icm),nnzp(icm),mmxp(icm)
!     +   ,nodemyp(mynum,icm),nodei0(mynum,icm),nodej0(mynum,icm),ifm,0,b)
!from parallel model:
!      call prgintrp(nnzp(icm),nnxp(icm),nnyp(icm),nnzp(icm),mmxp(icm)
!     +   ,nodemyp(mynum,icm),nodei0(mynum,icm),nodej0(mynum,icm),ifm,-1,b)

subroutine prgintrp(n1c,n2c,n3c,m1c,m2c,m3c,i0c,j0c,ifm,initflg,mynum)

  use var_tables, only: &
       num_scalar,      &
       scalar_tab

  use mem_scratch, only: &
       scratch

  use mem_basic, only: &
       basic_g

  use mem_grid, only: &
       nnxp,          &
       nnyp,          &
       nnzp,          &
       maxnxp,        &
       maxnyp,        &
       maxnzp,        &
       grid_g,        &
       nxtnest,       &
       nnstbot,       &
       nnsttop,       &
       icorflg,       &
       jdim

  use mem_nestb, only: &
       nbounds

  implicit none

  integer, intent(in) :: n1c       ! unused!!!
  integer, intent(in) :: n2c
  integer, intent(in) :: n3c
  integer, intent(in) :: m1c
  integer, intent(in) :: m2c
  integer, intent(in) :: m3c
  integer, intent(in) :: i0c       ! unused!!!
  integer, intent(in) :: j0c       ! unused!!!
  integer, intent(in) :: initflg
  integer, intent(in) :: ifm
  integer, intent(in) :: mynum     ! unused!!!

  integer :: icm
  integer :: nf
  integer :: nc
  integer :: ierr
  real, allocatable :: scr1(:), scr2(:)
  character(len=10) :: c0, c1
  character(len=*), parameter :: h="**(prgintrp)**"

  ! scratch area

  allocate(scr1(maxnzp*maxnxp*maxnyp), stat=ierr)
  if (ierr /= 0) call fatal_error(h//" allocating scr1")
  allocate(scr2(maxnzp*maxnxp*maxnyp), stat=ierr)
  if (ierr /= 0) call fatal_error(h//" allocating scr2")

  !     Temporarily fill VT2DA with interpolated topography from coarser grid

  icm = nxtnest(ifm)
  if (icm .eq. 0) return
  if (initflg .eq. 1) then
     call fillscr(1,maxnxp,maxnyp,1,n2c,n3c,1,1  &
          ,scr1(1),grid_g(icm)%topt(1,1))
     call eintp(scr1(1),scr2(1)  &
          ,1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
     call fillvar(1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1  &
          ,scr2(1),scratch%vt2da(1))
  endif

  !     Interpolate atmospheric variables

  call fmint3(m1c,m2c,m3c,nnzp(ifm),nnxp(ifm),nnyp(ifm),maxnzp,maxnxp,maxnyp  &
       ,ifm,icm,nnstbot(ifm),nnsttop(ifm),jdim,initflg,1,1,'u'  &
       ,basic_g(icm)%uc(1,1,1),basic_g(ifm)%uc(1,1,1)  &
       ,basic_g(icm)%dn0u(1,1,1),basic_g(ifm)%dn0u(1,1,1)  &
       ,scr1(1),scr2(1)  &
       ,grid_g(ifm)%topt(1,1),scratch%vt2da(1)  &
       ,nbounds(ifm)%bux(1,1,1),nbounds(ifm)%buy(1,1,1)  &
       ,nbounds(ifm)%buz(1,1,1))

  if (jdim .eq. 1 .or. icorflg .eq. 1)  &
       call fmint3(m1c,m2c,m3c,nnzp(ifm),nnxp(ifm),nnyp(ifm),maxnzp,maxnxp,maxnyp  &
       ,ifm,icm,nnstbot(ifm),nnsttop(ifm),jdim,initflg,1,1,'v'  &
       ,basic_g(icm)%vc(1,1,1),basic_g(ifm)%vc(1,1,1)  &
       ,basic_g(icm)%dn0v(1,1,1),basic_g(ifm)%dn0v(1,1,1)  &
       ,scr1(1),scr2(1)  &
       ,grid_g(ifm)%topt(1,1),scratch%vt2da(1)  &
       ,nbounds(ifm)%bvx(1,1,1),nbounds(ifm)%bvy(1,1,1)  &
       ,nbounds(ifm)%bvz(1,1,1))

  call fmint3(m1c,m2c,m3c,nnzp(ifm),nnxp(ifm),nnyp(ifm),maxnzp,maxnxp,maxnyp  &
       ,ifm,icm,nnstbot(ifm),nnsttop(ifm),jdim,initflg,1,1,'w'  &
       ,basic_g(icm)%wc(1,1,1),basic_g(ifm)%wc(1,1,1)  &
       ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
       ,scr1(1),scr2(1)  &
       ,grid_g(ifm)%topt(1,1),scratch%vt2da(1)  &
       ,nbounds(ifm)%bwx(1,1,1),nbounds(ifm)%bwy(1,1,1)  &
       ,nbounds(ifm)%bwz(1,1,1))

  call fmint3(m1c,m2c,m3c,nnzp(ifm),nnxp(ifm),nnyp(ifm),maxnzp,maxnxp,maxnyp  &
       ,ifm,icm,nnstbot(ifm),nnsttop(ifm),jdim,initflg,0,1,'t'  &
       ,basic_g(icm)%pc(1,1,1),basic_g(ifm)%pc(1,1,1)  &
       ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
       ,scr1(1),scr2(1)  &
       ,grid_g(ifm)%topt(1,1),scratch%vt2da(1)  &
       ,nbounds(ifm)%bpx(1,1,1),nbounds(ifm)%bpy(1,1,1)  &
       ,nbounds(ifm)%bpz(1,1,1))

  do nf = 1,num_scalar(ifm)
     do nc = 1,num_scalar(icm)
        if (scalar_tab(nf,ifm)%name == scalar_tab(nc,icm)%name) then
           call fmint3(m1c,m2c,m3c,nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
                ,maxnzp,maxnxp,maxnyp,ifm,icm  &
                ,nnstbot(ifm),nnsttop(ifm),jdim,initflg,1,1,'t'  &
                ,scalar_tab(nc,icm)%var_p &
                ,scalar_tab(nf,ifm)%var_p  &
                ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
                ,scr1(1),scr2(1)  &
                ,grid_g(ifm)%topt(1,1),scratch%vt2da(1)  &
                ,nbounds(ifm)%bsx(1,1,1,nf),nbounds(ifm)%bsy(1,1,1,nf)  &
                ,nbounds(ifm)%bsz(1,1,1,nf))
        endif
     enddo
  enddo

  deallocate(scr1)
  deallocate(scr2)
end subroutine prgintrp

!******************************************************************************

!!$subroutine nstfeed(ifm,icm)
!!$
!!$  use var_tables, only: &
!!$       num_scalar,      &
!!$       scalar_tab
!!$
!!$  use mem_basic, only: &
!!$       basic_g
!!$
!!$  use mem_scratch, only: &
!!$       scratch
!!$
!!$  use node_mod, only: &
!!$       nodeia,        &
!!$       nodeiz,        &
!!$       nodeja,        &
!!$       nodejz,        &
!!$       mynum,         &
!!$       ibcon
!!$
!!$  use mem_grid, only: &
!!$       nnzp,          &
!!$       nnxp,          &
!!$       nnyp,          &
!!$       nnstbot,       &
!!$       nnsttop
!!$
!!$  implicit none
!!$
!!$  integer :: ifm,icm
!!$
!!$  integer :: nf,nc
!!$  real, pointer :: scalarf, scalarc
!!$  character(len=10) :: c0, c1
!!$  character(len=*), parameter :: h="**(nstfeed)**"
!!$
!!$  !     feed back the finer mesh to the coarser mesh.
!!$
!!$  if (icm .eq. 0) return
!!$
!!$  call fdback(basic_g(icm)%uc   (1,1,1)  ,basic_g(ifm)%uc   (1,1,1)  &
!!$       ,basic_g(icm)%dn0u (1,1,1)  ,basic_g(ifm)%dn0u (1,1,1)  &
!!$       ,nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm),nnyp(ifm)    &
!!$       ,ifm,'u',scratch%scr1(1))
!!$  if (nnstbot(icm) .eq. 1) then
!!$     call botset(nnzp(icm),nnxp(icm),nnyp(icm)  &
!!$          ,nodeia(mynum,icm),nodeiz(mynum,icm),nodeja(mynum,icm),nodejz(mynum,icm),ibcon  &
!!$          ,basic_g(icm)%uc(1,1,1),'U')
!!$  end if
!!$  if (nnsttop(icm) .eq. 1) then
!!$     call topset(nnzp(icm),nnxp(icm),nnyp(icm)  &
!!$          ,nodeia(mynum,icm),nodeiz(mynum,icm),nodeja(mynum,icm),nodejz(mynum,icm),ibcon  &
!!$          ,basic_g(icm)%uc(1,1,1),basic_g(icm)%uc(1,1,1),'U')
!!$  end if
!!$
!!$  call fdback(basic_g(icm)%vc   (1,1,1)  ,basic_g(ifm)%vc   (1,1,1)  &
!!$       ,basic_g(icm)%dn0v (1,1,1)  ,basic_g(ifm)%dn0v (1,1,1)  &
!!$       ,nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm),nnyp(ifm)    &
!!$       ,ifm,'v',scratch%scr1(1))
!!$  if (nnstbot(icm) .eq. 1) then
!!$     call botset(nnzp(icm),nnxp(icm),nnyp(icm)  &
!!$          ,nodeia(mynum,icm),nodeiz(mynum,icm),nodeja(mynum,icm),nodejz(mynum,icm),ibcon  &
!!$          ,basic_g(icm)%vc(1,1,1),'V')
!!$  end if
!!$  if (nnsttop(icm) .eq. 1) then
!!$     call topset(nnzp(icm),nnxp(icm),nnyp(icm)  &
!!$          ,nodeia(mynum,icm),nodeiz(mynum,icm),nodeja(mynum,icm),nodejz(mynum,icm),ibcon  &
!!$          ,basic_g(icm)%vc(1,1,1),basic_g(icm)%vc(1,1,1),'V')
!!$  end if
!!$
!!$  call fdback(basic_g(icm)%wc  (1,1,1)  ,basic_g(ifm)%wc  (1,1,1)  &
!!$       ,basic_g(icm)%dn0 (1,1,1)  ,basic_g(ifm)%dn0 (1,1,1)  &
!!$       ,nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
!!$       ,ifm,'w',scratch%scr1(1))
!!$
!!$  call fdback(basic_g(icm)%pc  (1,1,1)  ,basic_g(ifm)%pc   (1,1,1)  &
!!$       ,basic_g(icm)%dn0 (1,1,1)  ,basic_g(ifm)%dn0  (1,1,1)  &
!!$       ,nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm),nnyp(ifm)   &
!!$       ,ifm,'p',scratch%scr1(1))
!!$  if (nnstbot(icm) .eq. 1) then
!!$     call botset(nnzp(icm),nnxp(icm),nnyp(icm)  &
!!$          ,nodeia(mynum,icm),nodeiz(mynum,icm),nodeja(mynum,icm),nodejz(mynum,icm),ibcon  &
!!$          ,basic_g(icm)%pc(1,1,1),'P')
!!$  end if
!!$  if (nnsttop(icm) .eq. 1) then
!!$     call topset(nnzp(icm),nnxp(icm),nnyp(icm)  &
!!$          ,nodeia(mynum,icm),nodeiz(mynum,icm),nodeja(mynum,icm),nodejz(mynum,icm),ibcon  &
!!$          ,basic_g(icm)%pc(1,1,1),basic_g(icm)%pc(1,1,1),'P')
!!$  end if
!!$
!!$  do nf = 1,num_scalar(ifm)
!!$     do nc = 1,num_scalar(icm)
!!$        if (scalar_tab(nf,ifm)%name == scalar_tab(nc,icm)%name) then
!!$           scalarc => scalar_tab(nc,icm)%var_p
!!$           scalarf => scalar_tab(nf,ifm)%var_p
!!$
!!$           call fdback(scalarc,scalarf  &
!!$                ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
!!$                ,nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm)  &
!!$                ,nnyp(ifm),ifm,'t',scratch%scr1(1))
!!$           if (nnstbot(icm) .eq. 1) then
!!$              call botset(nnzp(icm),nnxp(icm),nnyp(icm)  &
!!$                   ,nodeia(mynum,icm),nodeiz(mynum,icm),nodeja(mynum,icm),nodejz(mynum,icm),ibcon,scalarc,'T')
!!$           end if
!!$           if (nnsttop(icm) .eq. 1) then
!!$              call topset(nnzp(icm),nnxp(icm),nnyp(icm)  &
!!$                   ,nodeia(mynum,icm),nodeiz(mynum,icm),nodeja(mynum,icm),nodejz(mynum,icm),ibcon,scalarc,scalarc,'T')
!!$           end if
!!$        endif
!!$     enddo
!!$  enddo
!!$end subroutine nstfeed

