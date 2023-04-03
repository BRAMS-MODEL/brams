!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine decomp_node(init)

  use grid_dims, only: &
       maxmach,        &
       maxgrds

  use mem_grid, only: &
       runtype,       &
       ngrids,        &
       nnxp,          &
       nnyp,          &
       nnzp,          &
       nnxyp

  use mem_basic, only: &
       basic_g

  use node_mod, only: &
       nmachs,        &
       mchnum,        &
       master_num,    &
       dumpUnit,      &
       ixb,           &
       ixe,           &
       iyb,           &
       iye

  implicit none
  ! Arguments:
  integer, intent(in) :: init
  ! Local Variables:
  integer :: ngr
  logical :: failed
  integer :: largestDomain
  real, allocatable :: fakeCPU(:)
  real, allocatable :: LoadWeight(:)
  integer :: nm
  logical, parameter :: dumpLocal=.false.
!!$
!!$  ! allocate local scratch area
!!$
!!$  largestDomain = maxval(nnxyp(1:ngrids))
!!$  allocate(LoadWeight(largestDomain))
!!$  allocate(fakeCPU(largestDomain))
!!$
!!$  ! try to read domain decomposition data from input file
!!$  call decomp_input_par(maxmach, maxgrds, nmachs, ngrids, nnxp, nnyp, &
!!$       ixb, ixe, iyb, iye, failed)
!!$
!!$  ! escape if domain decomposition specified by input file;
!!$  ! else, Decompose all grids into subdomains
!!$
!!$  if (failed) then
!!$
!!$     ! Obtain estimates of the fraction of computational time (LoadWeight) required
!!$     ! for each column in the region of the domain, and
!!$     ! Decompose the grid taking into account the work numbers.
!!$
!!$     if (init == 1) then
!!$
!!$        ! first domain decomposition occurs prior to allocate memory
!!$        ! for basic_g and scratch; due to that, use a fake cputime
!!$
!!$        fakeCPU = 1.0
!!$
!!$        do ngr = 1,ngrids
!!$           call est_time_par(nnxp(ngr),nnyp(ngr),LoadWeight(1), fakeCPU(1))
!!$           call decomp_par(nnxp(ngr),nnyp(ngr),nmachs,LoadWeight(1),  &
!!$                ixb(1,ngr),ixe(1,ngr),iyb(1,ngr),iye(1,ngr))
!!$        end do
!!$
!!$        if (runtype(1:7)=='MAKESFC' .or. runtype(1:9)=='MAKEVFILE') then
!!$           do nm = 1, nmachs
!!$              ixb(nm,1:ngrids) = 2
!!$              ixe(nm,1:ngrids) = nnxp(1:ngrids) - 1
!!$              iyb(nm,1:ngrids) = 2
!!$              iye(nm,1:ngrids) = nnyp(1:ngrids) - 1
!!$           enddo
!!$        endif
!!$
!!$     else
!!$
!!$        ! use collected cpu execution time for non-initialization calls
!!$        !**(JP)** but cpu time is not collected at a point basis! I do not
!!$        !         see how computing load may be attributed to a domain point.
!!$
!!$        do ngr = 1,ngrids
!!$
!!$           call est_time_par(nnxp(ngr),nnyp(ngr),LoadWeight(1),  &
!!$                basic_g(ngr)%cputime(1,1))
!!$
!!$           call decomp_par(nnxp(ngr),nnyp(ngr),nmachs,LoadWeight(1),  &
!!$                ixb(1,ngr),ixe(1,ngr),iyb(1,ngr),iye(1,ngr))
!!$        end do
!!$     end if
!!$  end if

  ! Compute various bounds for the subdomains, storing results at module node_mod

  call decomp_bounds_par(maxgrds, ngrids, nnxp, nnyp, nnzp, 1, 1)

  ! deallocate scratch areas

!!$  deallocate(LoadWeight)
!!$  deallocate(fakeCPU)
!!$
!!$  if (dumpLocal) then
!!$     if (mchnum==master_num) then
!!$        call domain_decomposition_dump(dumpUnit)
!!$     end if
!!$  end if

end subroutine decomp_node


subroutine NodePathsBuffAlloc()

  use grid_dims, only: &
       maxmach,        &
       maxgrds,        &
       nxpmax,         &
       nypmax

  use mem_grid, only: &
       ngrids,        &
       nnxp,          &
       nnyp,          &
       nnzp,          &
       nxtnest,       &
       ipm,           &
       jpm,           &
       ibnd,          &
       jbnd

  use var_tables, only: &
       num_var,         &
       vtab_r,          &
       num_scalar,      &
       scalar_tab

  use node_mod, only: &
       mynum,         &
       nmachs,        &
       machs,         &
       ipaths,        &
       iget_paths,    &
       nxbeg, nxend, nybeg, nyend, &
       ixb, ixe, iyb, iye, &
       nodeibcon, &
       nodebounds, &
       nbuff_nest, &
       node_buffs
  !--(DMK)-------------------------

  implicit none

  integer :: idn
  integer :: isn
  integer :: ng
  integer :: mp_nzp
  integer :: numbuff
  integer :: icm
  integer :: ifm
  integer :: nestvar
  integer :: nc
  integer :: nf
  integer :: npvar2
  integer :: npvar3
  integer :: nv
  integer :: num_lbc_buff
  integer :: num_nest_buff
  integer :: num_feed_buff
  integer :: itype
  integer :: i1
  integer :: i2
  integer :: j1
  integer :: j2
  integer :: ixy
  integer :: ixyz
  integer :: memf


  ! Determine node sending paths and numbers of receiving nodes

  call node_paths_par(maxgrds, ngrids, nxpmax, nypmax, nnxp, nnyp,            &
       nxtnest, maxmach, nmachs, machs, nodeibcon, nxbeg, nxend, nybeg, nyend,&
                                ! For reproducibility - Saulo Barros
       ipm, jpm, ixb, ixe, iyb, iye, nodebounds, &
       2, ibnd, jbnd)



  !     Compute  send and receive buffer sizes. These will be maximum of
  !       long timestep, turbulence, nest boundaries, and nest feedback.
  !       Small timestep will use same buffers as they are always smaller.

  node_buffs(:)%nsend = 0
  node_buffs(:)%nrecv = 0


  mp_nzp=0
  do ng=1,ngrids
     mp_nzp=max(mp_nzp,nnzp(ng))
  enddo


  do ng=1,ngrids

     !          Find number of nested variables to be communicated.
     icm=ng
     ifm=ng
     if(ng /= 1) icm=nxtnest(ifm)
     nestvar=4
     do nf=1,num_scalar(ifm)
        do nc=1,num_scalar(icm)
           if(scalar_tab(nf,ifm)%name==scalar_tab(nc,icm)%name)  &
                nestvar=nestvar+1
        enddo
     enddo

     !  Find number of lbc variables to be communicated.
     npvar3=0 ; npvar2=0
     do nv = 1,num_var(ng)
        if(vtab_r(nv,ng)%impt1 == 1 ) then
           if (vtab_r(nv,ng)%idim_type==2) npvar2=npvar2+1
           if (vtab_r(nv,ng)%idim_type==3) npvar3=npvar3+1
        endif
     enddo

     isn = mynum
     do idn=1,nmachs
        num_lbc_buff=0
        num_nest_buff=0
        num_feed_buff=0

        itype=1
        i1=ipaths(1,itype,ng,idn)
        i2=ipaths(2,itype,ng,idn)
        j1=ipaths(3,itype,ng,idn)
        j2=ipaths(4,itype,ng,idn)

        if(i1/=0) then
           ixy          = (i2-i1+1)*(j2-j1+1)
           ixyz         = (i2-i1+1)*(j2-j1+1)*(mp_nzp)
           num_lbc_buff = ixyz*npvar3+ixy*npvar2 + 2*(npvar3+npvar2+100)
        endif

        itype=5
        i1=ipaths(1,itype,ng,idn)
        i2=ipaths(2,itype,ng,idn)
        j1=ipaths(3,itype,ng,idn)
        j2=ipaths(4,itype,ng,idn)

        if(i1/=0) then
           ixyz          = (i2-i1+1)*(j2-j1+1)*(mp_nzp)
           num_nest_buff = ixyz*nestvar+2*(nestvar+100)
        endif

        !  itype=6
        itype=7                  !For reproducibility - Saulo Barros
        i1=ipaths(1,itype,ng,idn)
        i2=ipaths(2,itype,ng,idn)
        j1=ipaths(3,itype,ng,idn)
        j2=ipaths(4,itype,ng,idn)

        if(i1/=0) then
           ixyz          = (i2-i1+1)*(j2-j1+1)*(mp_nzp)
           num_feed_buff = ixyz*nestvar+2*(nestvar+100)
        endif

        node_buffs(idn)%nsend = &
             max(node_buffs(idn)%nsend, &
             num_lbc_buff, num_nest_buff,     &
             num_feed_buff)

        node_buffs(idn)%nrecv = &
             max(node_buffs(idn)%nrecv, &
             num_lbc_buff, num_nest_buff,     &
             num_feed_buff)

     enddo !idn

  enddo !ng

  !       Check nest boundary receive buffer size

  itype=5
  do idn=1,nmachs


     if (mynum==idn) nbuff_nest = 0

     do ng=1,ngrids
        numbuff=0
        icm=nxtnest(ng)

        i1=ipaths(1,itype,ng,idn)
        i2=ipaths(2,itype,ng,idn)
        j1=ipaths(3,itype,ng,idn)
        j2=ipaths(4,itype,ng,idn)

        memf    = (i2-i1+1)*(j2-j1+1)*(mp_nzp)*nestvar
        numbuff = numbuff+memf

        if (mynum==idn) nbuff_nest = max(nbuff_nest, numbuff)
     enddo !ng

  enddo !idn

end subroutine NodePathsBuffAlloc
!
!     ****************************************************************
!
subroutine decomp_bounds_par(maxgrds, ngrids, nnxp, nnyp, nnzp, nbndx, nbndy)

  use node_mod, only: &
       nmachs,        &
       nodeia,        &
       nodeiz,        &
       nodeja,        &
       nodejz,        &
       nodei0,        &
       nodej0,        &
  !--(DMK)---------------------
       nodeibcon,     &
       nxbeg,         &
       nxend,         &
       nybeg,         &
       nyend,         &
       ixb,           &
       ixe,           &
       iyb,           &
       iye,           &
       npxy,          &
                                !--(DMK)---------------------
       nodemxp,       &
       nodemyp,       &
       nodemzp

  implicit none
  ! Arguments:
  integer, intent(in) :: maxgrds
  integer, intent(in) :: ngrids
  integer, intent(in) :: nnxp(maxgrds)
  integer, intent(in) :: nnyp(maxgrds)
  integer, intent(in) :: nnzp(maxgrds)
  integer, intent(in) :: nbndx
  integer, intent(in) :: nbndy
  ! Local variables:
  integer :: ng, nm, nx, ny

  !        Compute various subdomain boundary numbers for the nodes
  !             nxbeg,nybeg,nxend,nyend - portions of full domain that nodes will have
  !                                     - includes overlap region
  !             nodei0,nodej0  - subdomain offsets relative to full domain
  !             nodeia,nodeiz,nodeja,nodejz - subdomain "compute" points,
  !                         or normal thermodynamic tendency points (2-nx, 2-ny for
  !                          non-parallel run
  !             nodeibcon - flag denoting if real boundary is on subdomain
  !                       bit 1=west, bit 2=east, bit 3=south, bit 4=north

  do ng=1,ngrids

     nx=nnxp(ng)
     ny=nnyp(ng)

     do nm=1,nmachs

        nodeibcon(nm,ng)=0

        if(ixb(nm,ng).eq.2) then
           nxbeg(nm,ng)=1
           nodei0(nm,ng)=0
           nodeia(nm,ng)=2
           nodeibcon(nm,ng)=nodeibcon(nm,ng)+1
        else
           nxbeg(nm,ng)=ixb(nm,ng)-nbndx
           nodei0(nm,ng)=nxbeg(nm,ng)-1
           nodeia(nm,ng)=nbndx+1
        endif

        if(ixe(nm,ng).eq.nx-1) then
           nxend(nm,ng)=nx
           nodeibcon(nm,ng)=nodeibcon(nm,ng)+2
           nodeiz(nm,ng)=(ixe(nm,ng)-ixb(nm,ng))+nodeia(nm,ng)
        else
           nxend(nm,ng)=ixe(nm,ng)+nbndx
           nodeiz(nm,ng)=(ixe(nm,ng)-ixb(nm,ng))+nodeia(nm,ng)
        endif

        if(iyb(nm,ng).eq.2) then
           nybeg(nm,ng)=1
           nodej0(nm,ng)=0
           nodeja(nm,ng)=2
           nodeibcon(nm,ng)=nodeibcon(nm,ng)+4
        else
           nybeg(nm,ng)=iyb(nm,ng)-nbndy
           nodej0(nm,ng)=nybeg(nm,ng)-1
           nodeja(nm,ng)=nbndy+1
        endif

        if(iye(nm,ng).eq.ny-1) then
           nyend(nm,ng)=ny
           nodeibcon(nm,ng)=nodeibcon(nm,ng)+8
           nodejz(nm,ng)=(iye(nm,ng)-iyb(nm,ng))+nodeja(nm,ng)
        else
           nyend(nm,ng)=iye(nm,ng)+nbndy
           nodejz(nm,ng)=(iye(nm,ng)-iyb(nm,ng))+nodeja(nm,ng)
        endif
        !print *,'LFR-DBG: ',nm,ng,nxbeg(nm,ng),nxend(nm,ng),nybeg(nm,ng),nyend(nm,ng)
        !call flush(6)
     enddo

     do nm=1,nmachs
        npxy(nm,ng)=(nxend(nm,ng)-nxbeg(nm,ng)+1)  &
             *(nyend(nm,ng)-nybeg(nm,ng)+1)
     enddo
  enddo

  nodemxp(:,:) = nxend(:,:) - nxbeg(:,:) + 1
  nodemyp(:,:) = nyend(:,:) - nybeg(:,:) + 1
  !do nm=1,nmachs
  !  print *,'LFR-DBG inside decomp: ',nodemxp(nm,1); call flush(6)
  !end do
  nodemzp=0
  do ng = 1, ngrids
     nodemzp(:,ng) = nnzp(ng)
  end do

end subroutine decomp_bounds_par

!============================================================================


subroutine decomp_par(nxp,nyp,nodes,work,ixb,ixe,iyb,iye)

  implicit none
  ! Arguments:
  integer, intent(in)  :: nxp
  integer, intent(in)  :: nyp
  integer, intent(in)  :: nodes
  real,    intent(in)  :: work(nxp,nyp)
  integer, intent(out) :: ixb(nodes)
  integer, intent(out) :: ixe(nodes)
  integer, intent(out) :: iyb(nodes)
  integer, intent(out) :: iye(nodes)
  ! Local variables:
  real :: workrow(nyp)
  real :: workcol(nxp)
  real :: workload(nodes)
  real :: workblock(nodes)
  integer :: jrows(nodes)
  integer :: jrow(nodes)
  integer :: nblocks(nodes)
  real :: relspeed(nodes)

  integer :: inode,i,j,islab,jnodes,nslabs,min_blocks,nbigslabs,iblock &
       ,jnode,knode
  real :: anodes,aslabs,totspeed,workdom,workaccum,worksofar &
       ,slabspeed,workslab

  ! default relspeed = 1.0 for nodes of uniform speed.
  relspeed = 1.0

  ! This routine decomposes grid domains of size (nnxp,nnyp) into a number,
  ! specified by nodes, of rectangular subdomains.  The convention is followed
  ! that any internal boundaries (between subdomains) that are parallel to
  ! the x-axis run continuously across the full domain, while boundaries
  ! parallel to the y-axis may or may not run the full distance across the
  ! domain.  For convenience, regions of the domain bounded by adjacent
  ! east-west internal boundaries are termed "slabs", while smaller divisions
  ! within each slab are termed "blocks".  Each block is required to have
  ! a minimum dimension of 6 by 6 grid cells.  If this cannot be satisfied
  ! with the given input parameters, the subroutine stops.


  ! Estimate the number of slabs to be used (aslabs), and compute a final
  ! nearest integer value (nslabs) which is limited to allowable values.
  ! Zero out array for accumulating number of columns for each node.

  anodes = float(nodes)
  aslabs = sqrt(anodes * float(nyp) / float(nxp))
  nslabs = min(nodes,max(1,nint(aslabs)))

  !          print*, 'nslabs',nslabs

  totspeed = 0.
  do inode = 1,nodes
     ixe(inode) = 0
     totspeed = totspeed + relspeed(inode)
  enddo

  !          print*, 'totspeed',totspeed

  ! Compute total work load over each row and over entire domain.

  workdom = 0.
  do j = 1,nyp
     workrow(j) = 0.
     do i = 1,nxp
        workrow(j) = workrow(j) + work(i,j)
     enddo
     workdom = workdom + workrow(j)

     !          print*, 'j,workdom,workrow(j)',j,workdom,workrow(j)

  enddo
  workrow(2) = workrow(2) + workrow(1)
  workrow(nyp-1) = workrow(nyp-1) + workrow(nyp)

  ! Determine number of blocks and the average workload for each slab.

  min_blocks = nodes / nslabs
  nbigslabs = nodes - min_blocks * nslabs
  inode = 0
  do islab = 1,nslabs
     workload(islab) = 0.
     nblocks(islab) = min_blocks
     if (islab .le. nbigslabs) nblocks(islab) = min_blocks + 1
     do iblock = 1,nblocks(islab)
        inode = inode + 1
        workload(islab) = workload(islab)  &
             + workdom * relspeed(inode) / totspeed

        !           print*, 'islab,iblock,workload(islab),workdom,inode'
        !           print*,  islab,iblock,workload(islab),workdom,inode

     enddo
  enddo

  ! Assign all j-rows to their respective slabs in a way that balances the work
  ! load among slabs according to their respective numbers of nodes (blocks).
  ! The array jrows counts the number of rows in each slab, and the array
  ! jrow is the index of the southernmost row in each slab.

  do islab = 1,nslabs
     jrows(islab) = 0
  enddo

  workaccum = 0.
  worksofar = 0.
  islab = 0

  do j = 2,nyp-1
     workaccum = workaccum + workrow(j)
     if (workaccum - .5 * workrow(j) .gt. worksofar .and.  &
          islab .lt. nslabs) then
        islab = islab + 1
        jrow(islab) = j
        worksofar = worksofar + workload(islab)
     endif
     jrows(islab) = jrows(islab) + 1
  enddo

  inode = 0
  jnode = 0
  knode = 0
  do islab = 1,nslabs

     ! Compute the total work load for each slab and for each i-column in the
     ! slab.

     slabspeed = 0.
     workslab = 0.
     do i = 1,nxp
        workcol(i) = 0.
        do j = jrow(islab),jrow(islab)+jrows(islab)-1
           workcol(i) = workcol(i) + work(i,j)
        enddo
        workslab = workslab + workcol(i)
     enddo
     workcol(2) = workcol(2) + workcol(1)
     workcol(nxp-1) = workcol(nxp-1) + workcol(nxp)

     ! Determine average workload for each block.

     do iblock = 1,nblocks(islab)
        jnode = jnode + 1
        slabspeed = slabspeed + relspeed(jnode)

        !           print*, 'r1:iblock,jnode,slabspeed,relspeed(jnode)'
        !           print*,     iblock,jnode,slabspeed,relspeed(jnode)

     enddo
     do iblock = 1,nblocks(islab)
        knode = knode + 1
        workblock(iblock) = workslab  &
             * relspeed(knode) / slabspeed

        !       print*, 'islab,iblock,workblock,workslab,relspeed,slabspeed'
        !       print*, islab,iblock,workblock(iblock),workslab,relspeed(knode)
        !     +       ,slabspeed
        !       print*, 'knode',knode

     enddo

     ! Assign the i-columns of each slab to their respective blocks in a way that
     ! balances the work load among the blocks.  The array ncols counts the number
     ! of i-columns on each node, and the array ncol is the index of the
     ! westernmost i-column on each node.

     workaccum = 0.
     worksofar = 0.

     iblock = 0
     do i = 2,nxp-1
        workaccum = workaccum + workcol(i)

        !        print*, 'islab',islab
        !        print*, 'i,workaccum,workcol(i),worksofar,iblock,nblocks'
        !        print*, i,workaccum,workcol(i),worksofar,iblock,nblocks(islab)

        if (workaccum - .5 * workcol(i) .gt. worksofar .and.  &
             iblock .lt. nblocks(islab)) then
           iblock = iblock + 1

           !ccccccc defining node variables here ccccccccccccccccccccccccccccccccccc
           inode = inode + 1
           iyb(inode) = jrow(islab)
           ixb(inode) = i
           iye(inode) = iyb(inode) + jrows(islab) - 1

           !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

           worksofar = worksofar + workblock(iblock)
        endif
        ixe(inode) = ixe(inode) + 1
     enddo
  enddo

  !ccccccc defining node variable here ccccccccccccccccccccccccccccccccccc
  do jnode = 1,nodes
     ixe(jnode) = ixb(jnode) + ixe(jnode) - 1

     !           print*, 'jnode,ixb,ixe',jnode,ixb(jnode),ixe(jnode)
     !     +        ,(ixe(jnode)-ixb(jnode)+1),(iye(jnode)-iyb(jnode)+1)
     !     +        ,(ixe(jnode)-ixb(jnode)+1)*(iye(jnode)-iyb(jnode)+1)

  enddo
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! Check to make sure that each subdomain has at least 2 interior
  ! rows and columns.

  do jnode = 1,nodes
     if (iye(jnode) - iyb(jnode) .lt. 1 .or.  &
          ixe(jnode) - ixb(jnode) .lt. 1) then
        print*, 'grid:',nxp,nyp,'  subdomain too small on node ',jnode
        print*, '(ixb,ixe,iyb,iye) = '  &
             ,ixb(jnode),ixe(jnode),iyb(jnode),iye(jnode)
        !stop 'small_nodes'
        call fatal_error('small_nodes')
     endif
  enddo

  return
end subroutine decomp_par

!******************************************************************************


subroutine est_time_par(nxp,nyp,work,cput)
  implicit none
  ! Arguments:
  integer, intent(in)  :: nxp
  integer, intent(in)  :: nyp
  real,    intent(out) :: work(nxp,nyp)
  real,    intent(in)  :: cput(nxp,nyp)
  ! Local variables:
  integer :: i,j
  real :: bfact

  ! Sample routine to fill work elements with values proportional to the time
  ! required to perform model operations.

  do j = 2,nyp-1
     do i = 2,nxp-1
        work(i,j) = cput(i,j)
     enddo
  enddo

  ! Fill real boundaries with .2 of interior points

  bfact=.2
  do j = 1,nyp
     work(1,j) = bfact * work(2,j)
     work(nxp,j) = bfact * work(nxp-1,j)
  enddo

  do i = 1,nxp
     work(i,1) = bfact * work(i,2)
     work(i,nyp) = bfact * work(i,nyp-1)
  enddo
end subroutine est_time_par


!******************************************************************************
subroutine node_paths_par(maxgrds,ngrids,nxpmax,nypmax,nnxp,nnyp  &
     ,nxtnest,maxmach,nmachs,machs,ibcflg,nxbeg,nxend,nybeg,nyend  &
     ,ipm,jpm,ixb,ixe,iyb,iye,nodebounds                   &  !For Reproducibility - Saulo Barros
     ,iwid,ibnd,jbnd)

  use node_mod, only: &
       dumpUnit,      &
       DumpNodePaths, &
       mynum,         &
       ipaths,        &
       iget_paths

  implicit none
  ! Arguments:
  integer, intent(in) :: maxgrds,ngrids,nxpmax,nypmax,maxmach,nmachs
  integer, intent(in) :: nnxp(maxgrds),nnyp(maxgrds), nxtnest(maxgrds), &
       machs(maxmach), ibcflg(nmachs,ngrids), nxbeg(nmachs,ngrids), &
       nxend(nmachs,ngrids), nybeg(nmachs,ngrids), nyend(nmachs,ngrids), &
       ipm(nxpmax,ngrids), jpm(nypmax,ngrids), ixb(nmachs,ngrids),  &
       ixe(nmachs,ngrids), iyb(nmachs,ngrids), iye(nmachs,ngrids)
  integer, intent(in) :: iwid, ibnd, jbnd
  integer, intent(inout) :: nodebounds(maxgrds,8)
  ! Local variables:
  integer :: is0t(3),is0u(3),is0v(3),js0t(3),js0u(3),js0v(3)
  integer :: ngr,isend_type,idn,isn,i,j,info,nnn &
       ,indt,indu,indv,nxp,nyp,id,jd,nijst,nijsu,nijsv &
       ,iselft,iselfu,iselfv,mijs,is,js
  integer :: ig, ih !For Reproducibility - Saulo Barros
  character(len=8) :: c0, c1, c2, c3, c4, c5, c6
  character(len=*), parameter :: h="**(node_paths_par)**"
  logical, parameter :: dumpLocal=.false.

  ! Zero out ipaths array

  ipaths(:,:,:,:) = 0; iget_paths(:,:,:) = 0

  ! for all grids:
  !  take the next coarser grid (parent grid)

  do ngr=1,ngrids
     nnn = nxtnest(ngr)

     ! For reproducibility - Saulo Barros
     nodebounds(ngr,1) = nnxp(ngr)
     nodebounds(ngr,2) = 1
     nodebounds(ngr,3) = nnyp(ngr)
     nodebounds(ngr,4) = 1
     nodebounds(ngr,5) = nnxp(ngr)
     nodebounds(ngr,6) = 1
     nodebounds(ngr,7) = nnyp(ngr)
     nodebounds(ngr,8) = 1

     ! for all message passing destination nodes

     do idn=1,nmachs

        if (dumpLocal) then
           write(c0,"(i8)") idn
           write(c2,"(i8)") nxbeg(idn,ngr)
           write(c3,"(i8)") nxend(idn,ngr)
           write(c4,"(i8)") nybeg(idn,ngr)
           write(c5,"(i8)") nyend(idn,ngr)
           write(dumpUnit,"(a)") h//" destination node "//trim(adjustl(c0))//&
                " has domain ["//trim(adjustl(c2))//":"//trim(adjustl(c3))//&
                ","//trim(adjustl(c4))//":"//trim(adjustl(c5))//"]"
        end if
        ! for all message passing source nodes

        do isn = 1,nmachs

           ! Compute sub-domains overlap among source and destination nodes.
           !
           ! Skips computation if source and destination are the same node,
           ! since there is full overlap but no message passing.
           !
           ! Skips computation if sub-domains of send and destination nodes do not overlap.

           if (isn==idn) go to 6

           if (dumpLocal) then
              write(c0,"(i8)") isn
              write(c2,"(i8)") ixb(isn,ngr)
              write(c3,"(i8)") ixe(isn,ngr)
              write(c4,"(i8)") iyb(isn,ngr)
              write(c5,"(i8)") iye(isn,ngr)
              write(dumpUnit,"(a)") h//" source node "//trim(adjustl(c0))//&
                   " owns domain ["//trim(adjustl(c2))//":"//trim(adjustl(c3))//&
                   ","//trim(adjustl(c4))//":"//trim(adjustl(c5))//"] + borders"
           end if

           if ( &
                nxbeg(idn,ngr) > ixe(isn,ngr) .or.  &
                nybeg(idn,ngr) > iye(isn,ngr) .or.  &
                nxend(idn,ngr) < ixb(isn,ngr) .or.  &
                nyend(idn,ngr) < iyb(isn,ngr)          ) go to 6

           ! There is intersection among source and destination node sub-domains.
           ! Destination will receive message from source at the current grid

           !**(JP)** change semantics of iget_paths
           !         it only indicates whenever there will be communication
           !           iget_paths_master(1,ngr,isn,idn) = machs(isn)
           if (idn==mynum) then
              iget_paths(1,ngr,isn) = 1
              if (dumpLocal) then
                 write(c0,"(i8)") idn
                 write(c1,"(i8)") isn
                 write(dumpUnit,"(a)") h//" node "//trim(adjustl(c0))//&
                   " will receive wp, pp on small ts "//&
                   "or ghost zone on large ts from node "//trim(adjustl(c1))
              end if
           end if
           !**(JP)** end of modification

           ! Sub-domain region overlap among source and destination processes
           ! for the current grid.

           if (isn==mynum) then
              ipaths(1,1,ngr,idn)=max(ixb(isn,ngr),nxbeg(idn,ngr))
              ipaths(2,1,ngr,idn)=min(ixe(isn,ngr),nxend(idn,ngr))
              ipaths(3,1,ngr,idn)=max(iyb(isn,ngr),nybeg(idn,ngr))
              ipaths(4,1,ngr,idn)=min(iye(isn,ngr),nyend(idn,ngr))
              ipaths(5,1,ngr,idn)=machs(idn)
           endif

           ! Expand ipaths to include [coarse] grid external boundary points.

           if (ipaths(1,1,ngr,idn)==2) then
              ipaths(1,1,ngr,idn) = 1
           endif

           if (ipaths(2,1,ngr,idn)==nnxp(ngr)-1) then
              ipaths(2,1,ngr,idn) = nnxp(ngr)
           endif

           if (ipaths(3,1,ngr,idn)==2) then
              ipaths(3,1,ngr,idn) = 1
           endif

           if (ipaths(4,1,ngr,idn)==nnyp(ngr)-1) then
              ipaths(4,1,ngr,idn) = nnyp(ngr)
           endif

           if (ipaths(1,1,ngr,idn) /= 0 .and. dumpLocal) then
              write(c0,"(i8)") isn
              write(c1,"(i8)") idn
              write(c2,"(i8)") ipaths(1,1,ngr,idn)
              write(c3,"(i8)") ipaths(2,1,ngr,idn)
              write(c4,"(i8)") ipaths(3,1,ngr,idn)
              write(c5,"(i8)") ipaths(4,1,ngr,idn)
              write(dumpUnit,"(a)") h//" node "//trim(adjustl(c0))//&
                   " sends to node "//trim(adjustl(c1))//&
                   " ["//trim(adjustl(c2))//":"//trim(adjustl(c3))//&
                   ","//trim(adjustl(c4))//":"//trim(adjustl(c5))//"]"//&
                   " wp, pp on small ts or ghost zone on large ts"
           end if

           ! Small timestep overlap regions for u

           if (ixb(idn,ngr)-1  ==  ixe(isn,ngr) .and.  &
                iye(idn,ngr)   >=  iyb(isn,ngr) .and.  &
                iyb(idn,ngr)   <=  iye(isn,ngr)) then

              !**(JP)** change semantics of iget_paths
              !         it only indicates whenever there will be communication
              !              iget_paths_master(2,ngr,isn,idn) = machs(isn)
              if (idn==mynum) then
                 iget_paths(2,ngr,isn) = 1
                 if (dumpLocal) then
                    write(c0,"(i8)") idn
                    write(c1,"(i8)") isn
                    write(dumpUnit,"(a)") h//" node "//trim(adjustl(c0))//&
                         " will receive up on small ts from node "//trim(adjustl(c1))
                 end if
              end if
              !**(JP)** end of modification

              if (isn==mynum) then
                 ipaths(1,2,ngr,idn)=ixe(isn,ngr)
                 ipaths(2,2,ngr,idn)=ixe(isn,ngr)
                 ipaths(3,2,ngr,idn)=max(iyb(isn,ngr),iyb(idn,ngr))
                 ipaths(4,2,ngr,idn)=min(iye(isn,ngr),iye(idn,ngr))
                 ipaths(5,2,ngr,idn)=machs(idn)

                 if (dumpLocal) then
                    write(c0,"(i8)") isn
                    write(c1,"(i8)") idn
                    write(c2,"(i8)") ipaths(1,2,ngr,idn)
                    write(c3,"(i8)") ipaths(2,2,ngr,idn)
                    write(c4,"(i8)") ipaths(3,2,ngr,idn)
                    write(c5,"(i8)") ipaths(4,2,ngr,idn)
                    write(dumpUnit,"(a)") h//" node "//&
                         trim(adjustl(c0))//&
                         " sends to node "//trim(adjustl(c1))//&
                         " ["//trim(adjustl(c2))//":"//trim(adjustl(c3))//&
                         ","//trim(adjustl(c4))//":"//trim(adjustl(c5))//"]"//&
                         " up on small ts "
                 end if
              end if
           end if

           ! Small timestep overlap regions for v

           if (iyb(idn,ngr)-1  ==  iye(isn,ngr) .and.  &
                ixe(idn,ngr)   >=  ixb(isn,ngr) .and.  &
                ixb(idn,ngr)   <=  ixe(isn,ngr)) then

              !**(JP)** change semantics of iget_paths
              !         it only indicates whenever there will be communication
              !              iget_paths_master(3,ngr,isn,idn) = machs(isn)
              if (idn==mynum) then
                 iget_paths(3,ngr,isn) = 1
                 if (dumpLocal) then
                    write(c0,"(i8)") idn
                    write(c1,"(i8)") isn
                    write(dumpUnit,"(a)") h//" node "//trim(adjustl(c0))//&
                         " will receive vp on small ts from node "//trim(adjustl(c1))
                 end if
              end if
              !**(JP)** end of modification

              if (isn==mynum) then
                 ipaths(1,3,ngr,idn)=max(ixb(isn,ngr),ixb(idn,ngr))
                 ipaths(2,3,ngr,idn)=min(ixe(isn,ngr),ixe(idn,ngr))
                 ipaths(3,3,ngr,idn)=iye(isn,ngr)
                 ipaths(4,3,ngr,idn)=iye(isn,ngr)
                 ipaths(5,3,ngr,idn)=machs(idn)

                 if (dumpLocal) then
                    write(c0,"(i8)") isn
                    write(c1,"(i8)") idn
                    write(c2,"(i8)") ipaths(1,3,ngr,idn)
                    write(c3,"(i8)") ipaths(2,3,ngr,idn)
                    write(c4,"(i8)") ipaths(3,3,ngr,idn)
                    write(c5,"(i8)") ipaths(4,3,ngr,idn)
                    write(dumpUnit,"(a)") h//" node "//&
                         trim(adjustl(c0))//&
                         " sends to node "//trim(adjustl(c1))//&
                         " ["//trim(adjustl(c2))//":"//trim(adjustl(c3))//&
                         ","//trim(adjustl(c4))//":"//trim(adjustl(c5))//"]"//&
                         " vp on small ts"
                 end if
              endif

           endif

           ! Small timestep overlap regions for pi'

           if (ixe(idn,ngr)+1  ==  ixb(isn,ngr) .and.  &
                iye(idn,ngr)   >=  iyb(isn,ngr) .and.  &
                iyb(idn,ngr)   <=  iye(isn,ngr)) then

              !**(JP)** change semantics of iget_paths
              !         it only indicates whenever there will be communication
              !              iget_paths_master(4,ngr,isn,idn) = machs(isn)
              if (idn==mynum) then
                 iget_paths(4,ngr,isn) = 1
                 if (dumpLocal) then
                    write(c0,"(i8)") idn
                    write(c1,"(i8)") isn
                    write(dumpUnit,"(a)") h//" node "//trim(adjustl(c0))//&
                         " will receive pp on small ts from node "//trim(adjustl(c1))
                 end if
              end if
              !**(JP)** end of modification

              if (isn==mynum) then
                 ipaths(1,4,ngr,idn)=ixb(isn,ngr)
                 ipaths(2,4,ngr,idn)=ixb(isn,ngr)
                 ipaths(3,4,ngr,idn)=max(iyb(isn,ngr),iyb(idn,ngr))
                 ipaths(4,4,ngr,idn)=min(iye(isn,ngr),iye(idn,ngr))
                 ipaths(5,4,ngr,idn)=machs(idn)

                 if (dumpLocal) then
                    write(c0,"(i8)") isn
                    write(c1,"(i8)") idn
                    write(c2,"(i8)") ipaths(1,4,ngr,idn)
                    write(c3,"(i8)") ipaths(2,4,ngr,idn)
                    write(c4,"(i8)") ipaths(3,4,ngr,idn)
                    write(c5,"(i8)") ipaths(4,4,ngr,idn)
                    write(dumpUnit,"(a)") h//" node "//&
                         trim(adjustl(c0))//&
                         " sends to node "//trim(adjustl(c1))//&
                         " ["//trim(adjustl(c2))//":"//trim(adjustl(c3))//&
                         ","//trim(adjustl(c4))//":"//trim(adjustl(c5))//"]"//&
                         " pp on small ts "
                 end if
              endif

           elseif (iye(idn,ngr)+1  ==  iyb(isn,ngr) .and.  &
                ixe(idn,ngr)       >=  ixb(isn,ngr) .and.  &
                ixb(idn,ngr)       <=  ixe(isn,ngr)) then

              !**(JP)** change semantics of iget_paths
              !         it only indicates whenever there will be communication
              !              iget_paths_master(4,ngr,isn,idn) = machs(isn)
              if (idn==mynum) then
                 iget_paths(4,ngr,isn) = 1
                 if (dumpLocal) then
                    write(c0,"(i8)") idn
                    write(c1,"(i8)") isn
                    write(dumpUnit,"(a)") h//" node "//trim(adjustl(c0))//&
                         " will receive pp on small ts from node "//trim(adjustl(c1))
                 end if
              end if
              !**(JP)** end of modification

              if (isn==mynum) then
                 ipaths(1,4,ngr,idn)=max(ixb(isn,ngr),ixb(idn,ngr))
                 ipaths(2,4,ngr,idn)=min(ixe(isn,ngr),ixe(idn,ngr))
                 ipaths(3,4,ngr,idn)=iyb(isn,ngr)
                 ipaths(4,4,ngr,idn)=iyb(isn,ngr)
                 ipaths(5,4,ngr,idn)=machs(idn)
                 if (dumpLocal) then
                    write(c0,"(i8)") isn
                    write(c1,"(i8)") idn
                    write(c2,"(i8)") ipaths(1,4,ngr,idn)
                    write(c3,"(i8)") ipaths(2,4,ngr,idn)
                    write(c4,"(i8)") ipaths(3,4,ngr,idn)
                    write(c5,"(i8)") ipaths(4,4,ngr,idn)
                    write(dumpUnit,"(a)") h//" node "//&
                         trim(adjustl(c0))//&
                         " sends to node "//trim(adjustl(c1))//&
                         " ["//trim(adjustl(c2))//":"//trim(adjustl(c3))//&
                         ","//trim(adjustl(c4))//":"//trim(adjustl(c5))//"]"//&
                         " pp on small ts "
                 end if
              end if
           end if

6          continue

           ! Coarse grid to fine grid interpolation communication

           if (nnn==0) go to 26

           if(ipm(ixb(idn,ngr)-1,ngr)-2    >  ixe(isn,nnn) .or.  &
                jpm(iyb(idn,ngr)-1,ngr)-2  >  iye(isn,nnn) .or.  &
                ipm(ixe(idn,ngr)+1,ngr)+1  <  ixb(isn,nnn) .or.  &
                jpm(iye(idn,ngr)+1,ngr)+1  <  iyb(isn,nnn)) go to 16

           !**(JP)** change semantics of iget_paths
           !         it only indicates whenever there will be communication
           !           iget_paths_master(5,ngr,isn,idn) = machs(isn)
           if (idn==mynum) then
              iget_paths(5,ngr,isn) = 1
              if (dumpLocal) then
                 write(c0,"(i8)") idn
                 write(c1,"(i8)") isn
                 write(dumpUnit,"(a)") h//" node "//trim(adjustl(c0))//&
                      " will receive coarse to fine grid data from node "//trim(adjustl(c1))
              end if
           end if
           !**(JP)** end of modification

           if (isn==mynum) then
              ipaths(1,5,ngr,idn) = max(ixb(isn,nnn),ipm(ixb(idn,ngr)-1,ngr)-2)
              ipaths(2,5,ngr,idn) = min(ixe(isn,nnn),ipm(ixe(idn,ngr)+1,ngr)+1)
              ipaths(3,5,ngr,idn) = max(iyb(isn,nnn),jpm(iyb(idn,ngr)-1,ngr)-2)
              ipaths(4,5,ngr,idn) = min(iye(isn,nnn),jpm(iye(idn,ngr)+1,ngr)+1)
              ipaths(5,5,ngr,idn) = machs(idn)
              if (dumpLocal) then
                 write(c0,"(i8)") isn
                 write(c1,"(i8)") idn
                 write(c2,"(i8)") ipaths(1,5,ngr,idn)
                 write(c3,"(i8)") ipaths(2,5,ngr,idn)
                 write(c4,"(i8)") ipaths(3,5,ngr,idn)
                 write(c5,"(i8)") ipaths(4,5,ngr,idn)
                 write(dumpUnit,"(a)") h//" node "//&
                      trim(adjustl(c0))//&
                      " sends to node "//trim(adjustl(c1))//&
                      " ["//trim(adjustl(c2))//":"//trim(adjustl(c3))//&
                      ","//trim(adjustl(c4))//":"//trim(adjustl(c5))//"]"//&
                      " data to interpolate from coarse to fine grid"
              end if
           end if


16         continue

           ! Fine grid to coarse grid averaging communication

           !c rtimh - micphys: in the following lines: change ixb, ixe, iyb, and iye
           !c to nxbeg, nxend, nybeg, and nyend on coarse grid destination nodes
           !c in order to keep internal CG boundaries up to date, since they have
           !c already been updated by their CG neighboring nodes??

           if(nxbeg(idn,nnn)    >  ipm(ixe(isn,ngr),ngr).or.  &
                nybeg(idn,nnn)  >  jpm(iye(isn,ngr),ngr).or.  &
                nxend(idn,nnn)  <  ipm(ixb(isn,ngr),ngr).or.  &
                nyend(idn,nnn)  <  jpm(iyb(isn,ngr),ngr)) go to 26

           !**(JP)** change semantics of iget_paths
           !         it only indicates whenever there will be communication
           !           iget_paths_master(6,ngr,isn,idn) = machs(isn)
           if (idn==mynum) then
              iget_paths(6,ngr,isn) = 1
              if (dumpLocal) then
                 write(c0,"(i8)") idn
                 write(c1,"(i8)") isn
                 write(dumpUnit,"(a)") h//" node "//trim(adjustl(c0))//&
                      " will receive fine to coarse grid data from node "//trim(adjustl(c1))
              end if
           end if

           !**(JP)** end of modification

           if (isn==mynum) then
              ipaths(1,6,ngr,idn) = max(ipm(ixb(isn,ngr),ngr),nxbeg(idn,nnn))
              ipaths(2,6,ngr,idn) = min(ipm(ixe(isn,ngr),ngr),nxend(idn,nnn))
              ipaths(3,6,ngr,idn) = max(jpm(iyb(isn,ngr),ngr),nybeg(idn,nnn))
              ipaths(4,6,ngr,idn) = min(jpm(iye(isn,ngr),ngr),nyend(idn,nnn))
              ipaths(5,6,ngr,idn) = machs(idn)

              if (dumpLocal) then
                 write(c0,"(i8)") isn
                 write(c1,"(i8)") idn
                 write(c2,"(i8)") ipaths(1,6,ngr,idn)
                 write(c3,"(i8)") ipaths(2,6,ngr,idn)
                 write(c4,"(i8)") ipaths(3,6,ngr,idn)
                 write(c5,"(i8)") ipaths(4,6,ngr,idn)
                 write(dumpUnit,"(a)") h//" node "//&
                      trim(adjustl(c0))//&
                      " sends to node "//trim(adjustl(c1))//&
                      " ["//trim(adjustl(c2))//":"//trim(adjustl(c3))//&
                      ","//trim(adjustl(c4))//":"//trim(adjustl(c5))//"]"//&
                      " data to interpolate from fine to coarse grid"
              end if
           end if

           ! For reproducibility - Saulo Barros
           if (ngr>=2) then
              if (idn==mynum) then
                 nodebounds(ngr-1,5) =   &
                      min(nodebounds(ngr-1,5),ipaths(1,6,ngr,idn))
                 nodebounds(ngr-1,6) =   &
                      max(nodebounds(ngr-1,6),ipaths(2,6,ngr,idn))
                 nodebounds(ngr-1,7) =   &
                      min(nodebounds(ngr-1,7),ipaths(3,6,ngr,idn))
                 nodebounds(ngr-1,8) =   &
                      max(nodebounds(ngr-1,8),ipaths(4,6,ngr,idn))
              endif

           endif
           ! --------------------------------------

           ! A second index value of 7 of the ipaths array is used to determine
           ! the loop limits in fdbackp for averaging the fm over the overlap
           ! between the cm node and fm node, rather than always over the full
           ! fm node.  It is not used for actually sending stuff.  The
           ! ipaths(*,6,*,*,*) part of the array is still used for sending the
           ! block of averaged cm points from the fm node to the cm node.

           do i = ixb(isn,ngr),ixe(isn,ngr)
              if (ipm(i,ngr)==ipaths(1,6,ngr,idn)) then
                 ipaths(1,7,ngr,idn) = i
                 exit
              endif
           enddo

           do i = ixe(isn,ngr),ixb(isn,ngr),-1
              if (ipm(i,ngr)==ipaths(2,6,ngr,idn)) then
                 ipaths(2,7,ngr,idn) = i
                 exit
              endif
           enddo

           do j = iyb(isn,ngr),iye(isn,ngr)
              if (jpm(j,ngr)==ipaths(3,6,ngr,idn)) then
                 ipaths(3,7,ngr,idn) = j
                 exit
              endif
           enddo

           do j = iye(isn,ngr),iyb(isn,ngr),-1
              if (jpm(j,ngr)==ipaths(4,6,ngr,idn)) then
                 ipaths(4,7,ngr,idn) = j
                 exit
              endif
           enddo

26         continue


           ! For reproducibility - Saulo Barros
           if (ipaths(1,7,ngr,idn)/=0) then

              if (idn==mynum) then
                 nodebounds(ngr,1) =   &
                      min(nodebounds(ngr,1),ipaths(1,7,ngr,idn))
                 nodebounds(ngr,2) =   &
                      max(nodebounds(ngr,2),ipaths(2,7,ngr,idn))
                 nodebounds(ngr,3) =   &
                      min(nodebounds(ngr,3),ipaths(3,7,ngr,idn))
                 nodebounds(ngr,4) =   &
                      max(nodebounds(ngr,4),ipaths(4,7,ngr,idn))
                 if (dumpLocal) then
                    write(c0,"(i8)") isn
                    write(c1,"(i8)") idn
                    write(c2,"(i8)") ipaths(1,7,ngr,idn)
                    write(c3,"(i8)") ipaths(2,7,ngr,idn)
                    write(c4,"(i8)") ipaths(3,7,ngr,idn)
                    write(c5,"(i8)") ipaths(4,7,ngr,idn)
                    write(dumpUnit,"(a)") h//" on messages of type 7, node "//&
                         trim(adjustl(c0))//&
                         " sends to node "//trim(adjustl(c1))//&
                         " ["//trim(adjustl(c2))//":"//trim(adjustl(c3))//&
                         ","//trim(adjustl(c4))//":"//trim(adjustl(c5))//"]"
              end if
              endif

           endif
           ! ----------------------------------------

        enddo

     enddo

     if (dumpLocal) then
        call DumpNodePaths(ngr)
     end if
  end do
end subroutine node_paths_par

!*******************************************************************************

subroutine decomp_input_par (maxmach, maxgrds, nmachs, ngr, nnxp, nnyp, ixb, ixe, iyb, iye, failed)

  ! par_decomp_input:
  !  reads domain decomposition from input file

  use domain_decomp, only: &
       domain_fname         ! intent(in)

  implicit none
  ! Arguments:
  integer, intent(in)  :: maxmach
  integer, intent(in)  :: maxgrds
  integer, intent(in)  :: nmachs
  integer, intent(in)  :: ngr
  integer, intent(in)  :: nnxp(maxgrds)
  integer, intent(in)  :: nnyp(maxgrds)
  integer, intent(out) :: ixb(nmachs,ngr)
  integer, intent(out) :: ixe(nmachs,ngr)
  integer, intent(out) :: iyb(nmachs,ngr)
  integer, intent(out) :: iye(nmachs,ngr)
  logical, intent(out) :: failed
  ! Local variables:
  integer, parameter :: un=48
  logical :: ex
  integer :: err
  integer :: grid
  integer :: node
  integer :: xbeg
  integer :: xend
  integer :: ybeg
  integer :: yend
  integer :: occur
  integer :: igr
  integer :: imachs
  integer :: iy
  integer :: ix
  integer, parameter :: UndefIndex=-234
  !character(len=*), parameter :: fname="Domains"
  character(len=*), parameter :: h="**(par_decomp_input)**"
  character(len=10) :: c0, c1, c2
  logical, parameter :: dumpLocal=.FALSE.

  ! fail is the default
  failed = .true.

  ! Checking if domain_fname is valid
  if (len_trim(domain_fname)==0) then
     if (dumpLocal) then
        write(*,"(a)") h//" Explicit domain decomposition file invalid!"
        write(*,"(a)") h//" Original domain decomposition mode activated"
     end if
     return
  endif

  if (dumpLocal) then
     write(*,"(a)") h//" fetching file "//domain_fname
  end if

  ! Checking file existence
  inquire(file=domain_fname(1:len_trim(domain_fname)), exist=ex)
  if (.not. ex) then
     write(*,"(a)") h//" Explicit domain decomposition file "//domain_fname// &
          " not found"
     write(*,"(a)") h//" Original domain decomposition mode activated"
     return
  end if

  ! initialization
  ixb(1:nmachs, 1:ngr) = UndefIndex
  ixe(1:nmachs, 1:ngr) = UndefIndex
  iyb(1:nmachs, 1:ngr) = UndefIndex
  iye(1:nmachs, 1:ngr) = UndefIndex

  ! store file contents at desired places
  open(un, file=domain_fname(1:len_trim(domain_fname)), status="old", action="read")
  do

     ! read until eof
     read(un, *, iostat=err) grid, node, xbeg, xend, ybeg, yend
     if (err /= 0) then
        exit
     end if

     ! Checking data consistency
     if (grid < 1 .or. grid > ngr) then
        write(c0,"(i10)") grid
        write(c1,"(i10)") ngr
        write(*,"(a)") h//" unexpected grid number on input file "//domain_fname
        write(*,"(a)") h//" grid range [1:"//trim(adjustl(c1))//&
             &"]; unexpected grid number="//trim(adjustl(c0))
        write(*,"(a)") h//" Explicit domain decomposition file invalid!"
        write(*,"(a)") h//" Original domain decomposition mode activated"
        return
     end if
     if (node < 1 .or. node > nmachs) then
        write(c0,"(i10)") node
        write(c1,"(i10)") nmachs
        write(*,"(a)") h//" unexpected node number on input file "//domain_fname
        write(*,"(a)") h//" node range [1:"//trim(adjustl(c1))//&
             &"]; unexpected node number="//trim(adjustl(c0))
        write(*,"(a)") h//" Explicit domain decomposition file invalid!"
        write(*,"(a)") h//" Original domain decomposition mode activated"
        return
     end if

     ! store

     ixb(node,grid) = xbeg
     ixe(node,grid) = xend
     iyb(node,grid) = ybeg
     iye(node,grid) = yend
  end do
  close(un)

  ! all data stored?
  ! it suffices to verify ixb (since ixe,iyb,iye are at the same read)
  do igr = 1, ngr
     do imachs = 1, nmachs
        if (ixb(imachs,igr) == UndefIndex) then
           write(c0,"(i10)") imachs
           write(c1,"(i10)") igr
           write(*,"(a)") h//" missing data for machine "&
                &//trim(adjustl(c0))//" and grid "//trim(adjustl(c1))
           write(*,"(a)") h//" on input file "//domain_fname
           write(*,"(a)") h//" Explicit domain decomposition file invalid!"
           write(*,"(a)") h//" Original domain decomposition mode activated"
           return
        end if
     end do
  end do

  ! input data partitions the domain?
  do igr = 1, ngr
     do iy = 2, nnyp(igr)-1
        do ix = 2, nnxp(igr)-1
           occur = count(&
                (ixb(1:nmachs,igr) <= ix) .and. &
                (ixe(1:nmachs,igr) >= ix) .and. &
                (iyb(1:nmachs,igr) <= iy) .and. &
                (iye(1:nmachs,igr) >= iy)          )
           if (occur /= 1) then
              write(*,"(a)") h//" domain decomposition is not a partition on file "//domain_fname
              write(*,"(a,3i10)") h//" point igr, ix, iy=",igr,ix,iy
              write(*,"(a,i10)") h//" number of occurences=",occur
              write(*,"(a)") h//" Explicit domain decomposition file invalid!"
              write(*,"(a)") h//" Original domain decomposition mode activated"
              return
           end if
        end do
     end do
  end do

  ! accepted domain decomposition by external file
  failed = .false.
  if (dumpLocal) then
     write(*,"(a)") h//" domain decomposition imposed by file "//domain_fname
  end if
end subroutine decomp_input_par




! domain_decomposition_dump: dumps domain decomposition at
!                            file "fname" (if present) or at
!                            stdout (if absent)





subroutine domain_decomposition_dump(fUnit)

  use mem_grid, only: &
       ngrids, &
       dyncore_flag, &
       order_h

  use node_mod, only: &
  !--(DMK)-------------
       ixb,           &
       ixe,           &
       iyb,           &
       iye,           &
  !--(DMK)-------------
       nmachs

  use dump, only: &
    dumpMessage

  implicit none
  integer, intent(in) :: fUnit

  include "constants.f90"
  integer :: ngr
  integer :: jnode
  integer :: ncols
  logical, parameter :: dumpLocal=.false.
  character(len=*), parameter :: h="**(domain_decomposition_dump)**"
  character(len=*), parameter :: header="**(domain_decomposition_dump)**"
  character(len=*), parameter :: version="para_init.f90"
  integer :: err
  integer :: colsInx
  integer :: colsInY
  integer :: minColValid
  character(len=8) :: c0

  ! dumps at selected unit

  open(unit=22,file='brams.log',position='append',action='write')
  write(unit=22,fmt="(a30,i8.8,a12,i8.8,a11)") ' === Domain decomposition for ',&
       ngrids,' grids with ',nmachs,' cores ==='

  write(unit=22,fmt="(a)") '      grid    node   x-beg   x-end   y-beg   y-end      cols'
  do ngr = 1, ngrids
     do jnode = 1,nmachs
        ncols= (1+ixe(jnode,ngr)-ixb(jnode,ngr))  &
             *(1+iye(jnode,ngr)-iyb(jnode,ngr))
        write(unit=22,fmt="('  ',7i8)") ngr,jnode,ixb(jnode,ngr),ixe(jnode,ngr)  &
             ,iyb(jnode,ngr),iye(jnode,ngr),ncols
     enddo
     write(unit=22,fmt="(a)")
  end do
  close(unit=22)

  !Checking number of columns for Runge-Kupta
  if(dyncore_flag==2) then
    minColValid=1+2*order_h
    do ngr = 1, ngrids
      do jnode = 1,nmachs
        colsInX=1+ixe(jnode,ngr)-ixb(jnode,ngr)
        colsInY=1+iye(jnode,ngr)-iyb(jnode,ngr)
        if(colsInx<minColValid) then
          write(*,fmt='(A)') 'ERR: Grid set incompatible for Runge Kupta'
          write(*,fmt='(A,I5.5,A,I5,A,I5)') 'Processor: ',jnode,', Cols:',colsInX,' for a minimum of ',minColValid
          err=dumpMessage(c_tty,c_yes,header,' !!! ',c_fatal&
          ,'Number of columns in X less than allowed. Please, increase NNXP in RAMSIN '//&
          'or decrease the amount of processors')
        endif
        if(colsInY<minColValid) then
          write(*,fmt='(A)') 'ERR: Grid set incompatible for Runge Kupta'
          write(*,fmt='(A,I5.5,A,I5,A,I5)') 'Processor: ',jnode,', Cols:',colsInY,' for a minimum of ',minColValid
          err=dumpMessage(c_tty,c_yes,header,' !!! ',c_fatal&
             ,'Number of columns in Y less than allowed. Please, increase NNYP in RAMSIN '//&
             'or decrease the amount of processors')
        endif
      enddo
    enddo
  endif

end subroutine domain_decomposition_dump




subroutine domain_decomposition_summary

  use ISO_FORTRAN_ENV

  use mem_grid, only: &
       nnxp, &
       nnyp, &
       ngrids

  use node_mod, only: &
       ixb,           &
       ixe,           &
       iyb,           &
       iye,           &
       nmachs

  implicit none
  integer :: ngr
  integer :: jnode
  integer :: ncols(nmachs)
  integer, allocatable :: bins(:)
  integer :: low, high, ibin
  logical, parameter :: dumpLocal=.false.
  character(len=*), parameter :: h="**(domain_decomposition_summary)**"
  character(len=8) :: c0, c1, c2, c3

  ! dumps at selected unit

  write(OUTPUT_UNIT,"(//,61('-'))")
  write(OUTPUT_UNIT,"(a)") &
       "-----------  DOMAIN DECOMPOSITION HISTOGRAM ----------------"
  do ngr = 1, ngrids

     ! grid points per rank

     jnode=1
     ncols(jnode)= (1+ixe(jnode,ngr)-ixb(jnode,ngr))  &
          *(1+iye(jnode,ngr)-iyb(jnode,ngr))
     low=ncols(1); high=ncols(1)
     do jnode = 2,nmachs
        ncols(jnode)= (1+ixe(jnode,ngr)-ixb(jnode,ngr))  &
             *(1+iye(jnode,ngr)-iyb(jnode,ngr))
        if (ncols(jnode) < low) then
           low=ncols(jnode)
        else if (ncols(jnode) > high) then
           high=ncols(jnode)
        end if
     enddo
     write(OUTPUT_UNIT,"(61('-'))")
     write(c0,"(i8)") ngr
     write(c1,"(i8)") nmachs
     write(c2,"(i8)") nnxp(ngr)
     write(c3,"(i8)") nnyp(ngr)
     write(OUTPUT_UNIT,"(a)") "Grid "//trim(adjustl(c0))//&
          " with "//trim(adjustl(c2))//" x "//trim(adjustl(c3))//&
          " surface points decomposed on "//trim(adjustl(c1))//" MPI ranks:"

     ! build histogram

     allocate(bins(low:high))
     bins = 0
     do jnode=1, nmachs
        bins(ncols(jnode)) = bins(ncols(jnode)) + 1
     end do

     ! dump histogram

     do ibin = low, high
        if (bins(ibin) == 0) cycle
        write(c0,"(i8)") bins(ibin)
        write(c1,"(i8)") ibin
        write(c2,"(i8)") (100*bins(ibin))/nmachs
        write(OUTPUT_UNIT,"(2x,a)") trim(adjustl(c0))//" ranks with "//&
             trim(adjustl(c1))//" surface points; "//trim(adjustl(c2))//"% of ranks"
     end do
     deallocate(bins)
     write(OUTPUT_UNIT,"(61('-'),//)")
  end do
end subroutine domain_decomposition_summary




integer function GetAvailableUnit()
  implicit none

  integer, parameter :: unitLow=10
  integer, parameter :: unitHigh=99
  integer :: iunit
  logical :: op
  character(len=*), parameter :: h="**(GetAvailableUnit)**"

  do iunit = unitLow, unitHigh
     inquire (unit=iunit, opened=op)
     if (.not. op) exit
  end do
  if (iunit > unitHigh) then
     call fatal_error(h//" Fortran i/o units exausted")
  else
     GetAvailableUnit = iunit
  end if
end function GetAvailableUnit
