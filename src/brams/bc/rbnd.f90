!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine latbnd()

use mem_tend
use mem_basic
use mem_grid
use node_mod

implicit none

!     This routine drives the computation of the radiative lateral
!     boundary condition for normal velocity on the coarsest grid
!     and the recomputation of all boundary tendencies on nested grids
!     after the first nested grid timestep.

if (nxtnest(ngrid) .eq. 0) then

!         Radiative and/or mesoscale compensation region lateral
!            boundary conditions.

   if (ibnd .le. 3 .or. jbnd .le. 3) then

      call latnormv(mzp,mxp,myp,ia,iz,ja,jz,ibcon                  &
         ,grid_g(ngrid)%lpu  (1,1)    ,grid_g(ngrid)%lpv(1,1)      &
         ,basic_g(ngrid)%up  (1,1,1)  ,basic_g(ngrid)%uc  (1,1,1)  &
         ,tend%ut            (1)      ,basic_g(ngrid)%vp  (1,1,1)  &
         ,basic_g(ngrid)%vc  (1,1,1)  ,tend%vt            (1)      &
         ,grid_g(ngrid)%dxt  (1,1)    ,grid_g(ngrid)%dyt  (1,1)    )

   endif

endif
return
end

!     *****************************************************************

subroutine latnormv(m1,m2,m3,ia,iz,ja,jz,ibcon,lpu_R,lpv_R  &
   ,up,uc,ut,vp,vc,vt,dxt,dyt)

use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ibcon,i,j,k
real, dimension(m2,m3) :: lpu_R,lpv_R

real :: dxl,dxr,cphx,cphy
real, dimension(m1,m2,m3) :: up,uc,ut,vp,vc,vt
real, dimension(m2,m3) :: dxt,dyt
integer, dimension(m2,m3) :: lpu,lpv

!     This routine ultimately updates tendencies at lateral boundaries
!     after first diagnosing appropriate phase speeds.
!
!     IBND and JBND are flags for the radiative type in the X and Y
!     direction. Their meaning is:
!
!        IBND=1......Klemp-Wilhelmson (1978) type; phase speed given
!                    by CPHAS
!        IBND=2......Klemp-Lilly (1980) type; doppler shifted phase
!                    speed constant with height and diagnosed from
!                    average of Orlanski speeds, i.e. function of
!                    only (X,Y)
!        IBND=3......Orlanski(1974) type; Phase speeds diagnosed
!                    from local conditions and function of (x,y,z)

!     Calculate the diagnostic phase
!     speed. The Orlanski(1976) leapfrog method is to use three time
!     levels of information, namely the T-2, T-1, and T level to
!     evaluate the phase speed given by - du/dt / du/dx = u + C.
!     If this is to be an Orlanski or Klemp-Lilly type boundary in the x
!     direction then this following diagnostic procedure is necessary.

!     If this is the first call to a routine, initialize the phase
!       speed arrays if necessary.

lpu=int(lpu_R);lpv=int(lpv_R)

if (ibcon.eq.0) return

!     first compute "X" boundaries.

if (iand(ibcon,1) .ne. 0) then
   do j = 1,m3
      dxl = 1. / (dtlv * dxt(2,j))
      do k = lpu(1,j),m1

         cphx = min(0.,max(-dxl,(up(k,1,j)-cphas)))
         ut(k,1,j) = ut(k,1,j) - cphx * dxt(2,j)  &
            * (up(k,2,j) + ut(k,2,j) * dtlv - up(k,1,j))

      enddo
   enddo
endif

if (iand(ibcon,2) .ne. 0) then
   do j = 1,m3
      dxr = 1. / (dtlv * dxt(m2-1,j))
      do k = lpu(m2-1,j),m1

         cphx = max(0.,min(dxr,(up(k,m2-1,j)+cphas)))
         ut(k,m2-1,j) = ut(k,m2-1,j) - cphx * dxt(m2-1,j)  &
            * (up(k,m2-1,j) - (up(k,m2-2,j) + ut(k,m2-2,j) * dtlv))

      enddo
   enddo
endif

!     South and north boundaries.

if (jdim .eq. 1) then

   if (iand(ibcon,4) .ne. 0) then
      do i = 1,m2
         dxl = 1. / (dtlv * dyt(i,2))
         do k = lpv(i,1),m1
            cphy = min(0.,max(-dxl,(vp(k,i,1)-cphas)))
            vt(k,i,1) = vt(k,i,1) - cphy * dyt(i,2)  &
               * (vp(k,i,2) + vt(k,i,2) * dtlv - vp(k,i,1))
         enddo
      enddo
   endif

   if (iand(ibcon,8) .ne. 0) then
      do i = 1,m2
         dxr = 1. / (dtlv * dyt(i,m3-1))
         do k = lpv(i,m3-1),m1
            cphy = max(0.,min(dxr,(vp(k,i,m3-1)+cphas)))
            vt(k,i,m3-1) = vt(k,i,m3-1) - cphy * dyt(i,m3-1)  &
               * (vp(k,i,m3-1) - (vp(k,i,m3-2) + vt(k,i,m3-2) * dtlv))
         enddo
      enddo
   endif

endif
return
end

!*****************************************************************************

subroutine vpsets()

use mem_basic
use mem_grid
use node_mod

implicit none

if (nxtnest(ngrid) .eq. 0) then
   call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'U'              &
      ,basic_g(ngrid)%up (1,1,1)  ,basic_g(ngrid)%up (1,1,1)  &
      ,basic_g(ngrid)%vp (1,1,1)  ,grid_g(ngrid)%dxu (1,1)    &
      ,grid_g(ngrid)%dxm (1,1)    ,grid_g(ngrid)%dyv (1,1)    &
      ,grid_g(ngrid)%dym (1,1)    )
   call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'V'              &
      ,basic_g(ngrid)%vp (1,1,1)  ,basic_g(ngrid)%up (1,1,1)  &
      ,basic_g(ngrid)%vp (1,1,1)  ,grid_g(ngrid)%dxu (1,1)    &
      ,grid_g(ngrid)%dxm (1,1)    ,grid_g(ngrid)%dyv (1,1)    &
      ,grid_g(ngrid)%dym (1,1)    )
   call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'W'              &
      ,basic_g(ngrid)%wp (1,1,1)  ,basic_g(ngrid)%up (1,1,1)  &
      ,basic_g(ngrid)%vp (1,1,1)  ,grid_g(ngrid)%dxu (1,1)    &
      ,grid_g(ngrid)%dxm (1,1)    ,grid_g(ngrid)%dyv (1,1)    &
      ,grid_g(ngrid)%dym (1,1)    )
   call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'P'              &
      ,basic_g(ngrid)%pp (1,1,1)  ,basic_g(ngrid)%up (1,1,1)  &
      ,basic_g(ngrid)%vp (1,1,1)  ,grid_g(ngrid)%dxu (1,1)    &
      ,grid_g(ngrid)%dxm (1,1)    ,grid_g(ngrid)%dyv (1,1)    &
      ,grid_g(ngrid)%dym (1,1)    )
endif

if (nsttop .eq. 1) then
   call topset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
      ,basic_g(ngrid)%up(1,1,1)  ,basic_g(ngrid)%up(1,1,1),'U')
   call topset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
      ,basic_g(ngrid)%vp(1,1,1)  ,basic_g(ngrid)%vp(1,1,1),'V')
endif

if (nstbot .eq. 1) then
   if (if_adap == 0) then
      call botset(mzp,mxp,myp,ia,iz,ja,jz,ibcon   &
         ,basic_g(ngrid)%up(1,1,1),'U')
      call botset(mzp,mxp,myp,ia,iz,ja,jz,ibcon   &
         ,basic_g(ngrid)%vp(1,1,1),'V')
   else
      call botset_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon,grid_g(ngrid)%lpu(1,1)  &
         ,basic_g(ngrid)%up(1,1,1),'U')
      call botset_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon,grid_g(ngrid)%lpv(1,1)  &
         ,basic_g(ngrid)%vp(1,1,1),'V')
   endif
endif

call topset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
      ,basic_g(ngrid)%pp(1,1,1)  ,basic_g(ngrid)%pp(1,1,1),'P')

if (if_adap == 0) then
   call botset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
      ,basic_g(ngrid)%pp(1,1,1),'P')
else
   call botset_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon,grid_g(ngrid)%lpw(1,1)  &
      ,basic_g(ngrid)%pp(1,1,1),'P')
endif

call dumset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
      ,basic_g(ngrid)%wp(1,1,1),'W')

call dumset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
      ,basic_g(ngrid)%up(1,1,1),'U')
call dumset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
      ,basic_g(ngrid)%vp(1,1,1),'V')

if (nxtnest(ngrid) .eq. 0) then
   call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'U'  &
      ,basic_g(ngrid)%uc (1,1,1)  ,basic_g(ngrid)%up (1,1,1)  &
      ,basic_g(ngrid)%vp (1,1,1)  ,grid_g(ngrid)%dxu (1,1)    &
      ,grid_g(ngrid)%dxm (1,1)    ,grid_g(ngrid)%dyv (1,1)    &
      ,grid_g(ngrid)%dym (1,1)    )
   call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'V'  &
      ,basic_g(ngrid)%vc (1,1,1)  ,basic_g(ngrid)%up (1,1,1)  &
      ,basic_g(ngrid)%vp (1,1,1)  ,grid_g(ngrid)%dxu (1,1)    &
      ,grid_g(ngrid)%dxm (1,1)    ,grid_g(ngrid)%dyv (1,1)    &
      ,grid_g(ngrid)%dym (1,1)    )
   call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon ,'W' &
      ,basic_g(ngrid)%wc (1,1,1)  ,basic_g(ngrid)%up (1,1,1)  &
      ,basic_g(ngrid)%vp (1,1,1)  ,grid_g(ngrid)%dxu (1,1)    &
      ,grid_g(ngrid)%dxm (1,1)    ,grid_g(ngrid)%dyv (1,1)    &
      ,grid_g(ngrid)%dym (1,1)    )
   call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'P'  &
      ,basic_g(ngrid)%pc (1,1,1)  ,basic_g(ngrid)%up (1,1,1)  &
      ,basic_g(ngrid)%vp (1,1,1)  ,grid_g(ngrid)%dxu (1,1)    &
      ,grid_g(ngrid)%dxm (1,1)    ,grid_g(ngrid)%dyv (1,1)    &
      ,grid_g(ngrid)%dym (1,1)    )
endif

if (nsttop .eq. 1) then
   call topset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
      ,basic_g(ngrid)%uc(1,1,1)  ,basic_g(ngrid)%uc(1,1,1),'U')
   call topset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
      ,basic_g(ngrid)%vc(1,1,1)  ,basic_g(ngrid)%vc(1,1,1),'V')
endif

if (nstbot .eq. 1) then
   if (if_adap == 0) then
      call botset(mzp,mxp,myp,ia,iz,ja,jz,ibcon   &
         ,basic_g(ngrid)%uc(1,1,1),'U')
      call botset(mzp,mxp,myp,ia,iz,ja,jz,ibcon   &
         ,basic_g(ngrid)%vc(1,1,1),'V')
   else
      call botset_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon,grid_g(ngrid)%lpu(1,1)  &
         ,basic_g(ngrid)%uc(1,1,1),'U')
      call botset_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon,grid_g(ngrid)%lpv(1,1)  &
         ,basic_g(ngrid)%vc(1,1,1),'V')
   endif
endif


call topset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
      ,basic_g(ngrid)%pc(1,1,1)  &
      ,basic_g(ngrid)%pc(1,1,1),'P')

if (if_adap == 0) then
   call botset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
      ,basic_g(ngrid)%pc(1,1,1),'P')
else
   call botset_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon,grid_g(ngrid)%lpw(1,1)  &
      ,basic_g(ngrid)%pc(1,1,1),'P')
endif

call dumset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
      ,basic_g(ngrid)%wc(1,1,1),'W')
call dumset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
      ,basic_g(ngrid)%uc(1,1,1),'U')
call dumset(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
      ,basic_g(ngrid)%vc(1,1,1),'V')

return
end

!*****************************************************************************

subroutine trsets()

use var_tables
use mem_basic
use mem_grid
use mem_turb
use node_mod

use ccatt_start, only: ccatt
use mem_chem1, only: CHEMISTRY,NSPECIES_TRANSPORTED

implicit none

integer :: n, &

!--(DMK-CCATT-INI)-----------------------------------------------------
     mxyzp
!--(DMK-CCATT-FIM)-----------------------------------------------------

real, pointer :: scalarp, scalart

!     Apply lateral, top, and bottom boundary conditions.

do n = 1,num_scalar(ngrid)
   scalarp => scalar_tab(n,ngrid)%var_p
   scalart => scalar_tab(n,ngrid)%var_t
 
   if (nxtnest(ngrid) .eq. 0) then
      call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'TR'   &
           ,scalarp                    ,basic_g(ngrid)%up (1,1,1)  &
           ,basic_g(ngrid)%vp (1,1,1)  ,grid_g(ngrid)%dxu (1,1)    &
           ,grid_g(ngrid)%dxm (1,1)    ,grid_g(ngrid)%dyv (1,1)    &
           ,grid_g(ngrid)%dym (1,1)    )
   endif
   if (nsttop .eq. 1)  then 
      
     if(n .le. (num_scalar(ngrid) - ( NADDSC + NSPECIES_TRANSPORTED) ))  then
           call topset (mzp,mxp,myp,ia,iz,ja,jz,ibcon,scalarp,scalarp,'T')
     else
           call topset2(mzp,mxp,myp,ia,iz,ja,jz,ibcon,scalarp,scalarp,'T')
     endif
   endif
   if (nstbot .eq. 1)  then
      if (if_adap == 0) then
         call botset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,scalarp,'T')
      else
         call botset_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
            ,grid_g(ngrid)%lpw(1,1),scalarp,'T')
      endif
   endif
enddo

!       Make sure all positive definite quantities remain such.

call tkeinit(mzp,mxp,myp)

call negadj1(mzp,mxp,myp)

!--(DMK-CCATT-INI)-----------------------------------------------------
!-srf for chem - aerosol quantities
if( ccatt == 1 .and. CHEMISTRY >= 0) then
   mxyzp = mzp*mxp*myp
   do n = 1,num_scalar(ngrid)
     
     if(n .le. (num_scalar(ngrid) - ( NADDSC + NSPECIES_TRANSPORTED) )) cycle 

      scalarp => scalar_tab(n,ngrid)%var_p
      
      call keep_tracers_nonneg(mxyzp,scalarp)

   enddo
endif
!--(DMK-CCATT-FIM)-----------------------------------------------------

return
end

!--(DMK-CCATT-INI)-----------------------------------------------------
!******************************************************************************
subroutine keep_tracers_nonneg(mxyzp,scp)
implicit none
integer :: mxyzp,i
real, dimension(mxyzp) :: scp

do i = 1,mxyzp
          scp(i) = max(0.,scp(i))
enddo
end subroutine keep_tracers_nonneg
!--(DMK-CCATT-FIM)-----------------------------------------------------

!******************************************************************************

subroutine latset(m1,m2,m3,ia,iz,ja,jz,ibcon,vnam,ap,uc,vc,dxu,dxm,dyv,dym)

use mem_grid
use mem_scratch

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ibcon,i,j,k,lbw,lbe,lbs,lbn
real :: thresh,dtlx,c1,dxr,dyr
real, dimension(m1,m2,m3) :: ap,uc,vc
real, dimension(m2,m3) :: dxu,dxm,dyv,dym
character(len=*) :: vnam

if (iand(ibcon,1) .gt. 0) lbw = ia - 1
if (iand(ibcon,2) .gt. 0) lbe = iz + 1
if (iand(ibcon,4) .gt. 0) lbs = ja - 1
if (iand(ibcon,8) .gt. 0) lbn = jz + 1

!!$print *, "DEBUG-ALF:latset:m1,m2,m3,lbw,lbe,lbs,lbn=", &
!!$     m1,m2,m3,lbw,lbe,lbs,lbn
!!$call flush(8)

thresh = 0.
if (vnam .eq. 'U' .or. vnam .eq. 'V' .or. vnam .eq. 'W' .or. vnam .eq.'P') then
   dtlx = dtlv
else
   dtlx = dtlt
endif

if (ibnd .ne. 4 .and. vnam .ne. 'U' .and. lsflg .ne. 3) then

!     Western and Eastern boundaries for zero gradient option

   if (lsflg .eq. 0) then
      if (iand(ibcon,1) .gt. 0) then
         do j = 1,m3
            do k = 1,m1
               ap(k,lbw,j) = ap(k,ia,j)
            enddo
         enddo
      endif
      if (iand(ibcon,2) .gt. 0) then
         do j = 1,m3
            do k = 1,m1
               ap(k,lbe,j) = ap(k,iz,j)
            enddo
         enddo
      endif
   else

!     Western boundary for lsflg = 1 or 2

      if (iand(ibcon,1) .gt. 0) then
         do j = 1,m3-1 !m3  !ALF
            if (vnam .eq. 'V') then
               dxr = dxm(ia,j) / dxm(lbw,j)
               c1 = .5 * dtlx * dxm(lbw,j)
               do k = 1,m1
                  vctr17(k) = -c1 * (uc(k,lbw,j) + uc(k,lbw,j+jdim))
               enddo
            elseif (vnam .eq. 'W') then
               dxr = dxu(ia,j) / dxu(lbw,j)
               c1 = .5 * dtlx * dxu(lbw,j)
               do k = 1,m1-1 !m1  !ALF
                  vctr17(k) = -c1 * (uc(k,lbw,j) + uc(k+1,lbw,j))
               enddo
            else
               dxr = dxu(ia,j) / dxu(lbw,j)
               c1 = dtlx * dxu(lbw,j)
               do k = 1,m1
                  vctr17(k) = -c1 * uc(k,lbw,j)
               enddo
            endif
            do k = 1,m1
               vctr18(k) = ap(k,ia,j) + dxr * (ap(k,ia,j) - ap(k,ia+1,j))
            enddo
            do k = 1,m1
               if (vctr17(k) .ge. thresh) then
                  ap(k,lbw,j) = vctr18(k)
               elseif (lsflg .eq. 1) then
                  ap(k,lbw,j) = ap(k,ia,j)
               endif
            enddo
         enddo
      endif

!     Eastern Boundary for LSFLG = 1 or 2

      if (iand(ibcon,2) .gt. 0) then
         do j = 1,m3-1 !m3  !ALF
            if (vnam .eq. 'V') then
               dxr = dxm(iz-1,j) / dxm(iz,j)
               c1 = .5 * dtlx * dxm(iz,j)
               do k = 1,m1
                  vctr17(k) = c1 * (uc(k,iz,j) + uc(k,iz,j+jdim))
               enddo
            elseif (vnam .eq. 'W') then
               dxr = dxu(iz-1,j) / dxu(iz,j)
               c1 = .5 * dtlx * dxu(iz,j)
               do k = 1,m1-1 !m1  !ALF
                  vctr17(k) = c1 * (uc(k,iz,j) + uc(k+1,iz,j))
               enddo
            else
               dxr = dxu(iz-1,j) / dxu(iz,j)
               c1 = dtlx * dxu(iz,j)
               do k = 1,m1
                  vctr17(k) = c1 * uc(k,iz,j)
               enddo
            endif
            do k = 1,m1
               vctr18(k) = ap(k,iz,j) + dxr * (ap(k,iz,j) - ap(k,iz-1,j))
            enddo
            do k = 1,m1
               if (vctr17(k) .ge. thresh) then
                  ap(k,lbe,j) = vctr18(k)
               elseif (lsflg .eq. 1) then
                  ap(k,lbe,j) = ap(k,iz,j)
               endif
            enddo
         enddo
      endif
   endif
endif

if(jdim.eq.1.and.jbnd.ne.4.and.vnam.ne.'V'.and.lsflg.ne.3)then

!     Southern and Northern boundaries for zero gradient option

  if (lsflg .eq. 0) then
     if (iand(ibcon,4) .gt. 0) then
        do i = 1,m2
           do k = 1,m1
              ap(k,i,lbs) = ap(k,i,ja)
           enddo
        enddo
     endif
     if (iand(ibcon,8) .gt. 0) then
        do i = 1,m2
           do k = 1,m1
              ap(k,i,lbn) = ap(k,i,jz)
           enddo
        enddo
     endif
  else

!     Southern boundary for LSFLG = 1 or 2

     if (iand(ibcon,4) .gt. 0) then
        do i = 1,m2-1 !m2 !ALF
           if (vnam .eq. 'U') then
              dyr = dym(i,ja) / dym(i,lbs)
              c1 = .5 * dtlx * dym(i,lbs)
              do k = 1,m1
                 vctr17(k) = -c1 * (vc(k,i,lbs) + vc(k,i+1,lbs))
              enddo
           elseif (vnam .eq. 'W') then
              dyr = dyv(i,ja) / dyv(i,lbs)
              c1 = .5 * dtlx * dyv(i,lbs)
              do k = 1,m1-1 !m1  !ALF
                 vctr17(k) = -c1 * (vc(k,i,lbs) + vc(k+1,i,lbs))
              enddo
           else
              dyr = dyv(i,ja) / dyv(i,lbs)
              c1 = dtlx * dyv(i,lbs)

!--(DMK-CCATT-INI)-----------------------------------------------------
              !srf - fix from rams 60
              do k = 1,m1
!--(DMK-CCATT-OLD)-----------------------------------------------------
!              do k = 1,nz
!--(DMK-CCATT-FIM)-----------------------------------------------------	      

                 vctr17(k) = -c1 * vc(k,i,lbs)
              enddo
           endif
           do k = 1,m1
              vctr18(k) = ap(k,i,ja) + dyr * (ap(k,i,ja) - ap(k,i,ja+1))
           enddo
           do k = 1,m1
              if (vctr17(k) .ge. thresh) then
                 ap(k,i,lbs) = vctr18(k)
              elseif (lsflg .eq. 1) then
                 ap(k,i,lbs) = ap(k,i,ja)
              endif
           enddo
        enddo
     endif

!     Northern Boundary for LSFLG = 1 or 2

     if (iand(ibcon,8) .gt. 0) then
        do i = 1,m2-1 !m2 !ALF
           if (vnam .eq. 'U') then
              dyr = dym(i,jz-1) / dym(i,jz)
              c1 = .5 * dtlx * dym(i,jz)
              do k = 1,m1
                 vctr17(k) = c1 * (vc(k,i,jz) + vc(k,i+1,jz))
              enddo
           elseif (vnam .eq. 'W') then
              dyr = dyv(i,jz-1) / dyv(i,jz)
              c1 = .5 * dtlx * dyv(i,jz)
              do k = 1,m1-1 !m1  !ALF
                 vctr17(k) = c1 * (vc(k,i,jz) + vc(k+1,i,jz))
              enddo
           else
              dyr = dyv(i,jz-1) / dyv(i,jz)
              c1 = dtlx * dyv(i,jz)
              do k = 1,m1
                 vctr17(k) = c1 * vc(k,i,jz)
              enddo
           endif
           do k = 1,m1
              vctr18(k) = ap(k,i,jz) + dyr * (ap(k,i,jz) - ap(k,i,jz-1))
           enddo
           do k = 1,m1
              if (vctr17(k) .ge. thresh) then
                 ap(k,i,lbn) = vctr18(k)
              elseif (lsflg .eq. 1) then
                 ap(k,i,lbn) = ap(k,i,jz)
              endif
           enddo
        enddo
     endif
  endif
endif

return
end subroutine latset

!     ******************************************************************

subroutine topset(m1,m2,m3,ia,iz,ja,jz,ibcon,ap,fa,vnam)

use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ibcon,i,j

real :: dzmr,dztr
real, dimension(m1,m2,m3) :: ap,fa
character(len=*) :: vnam

dzmr = dzm(m1-2) / dzm(m1-1)
dztr = dzt(m1-2) / dzt(m1-1)

!     Computation of all prognostic variables (other than W) at
!       level NZP by extrapolation from below

if (vnam .eq. 'U' .or. vnam .eq. 'V' .or. vnam .eq. 'P') then
   do j = 1,m3
      do i = 1,m2
         ap(m1,i,j) = ap(m1-1,i,j) + dzmr * (ap(m1-1,i,j) - ap(m1-2,i,j))
      enddo
   enddo
endif
if (vnam .eq. 'T') then
   do j = 1,m3
      do i = 1,m2
         ap(m1,i,j) = max(0.,ap(m1-1,i,j)+dzmr*(ap(m1-1,i,j)-ap(m1-2,i,j)))
      enddo
   enddo
endif

return
end
!     ******************************************************************

subroutine topset2(m1,m2,m3,ia,iz,ja,jz,ibcon,ap,fa,vnam)

use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ibcon,i,j

real :: dzmr,dztr
real, dimension(m1,m2,m3) :: ap,fa
character(len=*) :: vnam

!dzmr = dzm(m1-2) / dzm(m1-1)
!dztr = dzt(m1-2) / dzt(m1-1)

!     Computation of all prognostic variables (other than W) at
!       level NZP by extrapolation from below

!if (vnam .eq. 'U' .or. vnam .eq. 'V' .or. vnam .eq. 'P') then
!  do j = 1,m3
!      do i = 1,m2
!         ap(m1,i,j) = ap(m1-1,i,j) + dzmr * (ap(m1-1,i,j) - ap(m1-2,i,j))
!      enddo
!   enddo
!endif
!if (vnam .eq. 'T') then
   do j = 1,m3
      do i = 1,m2
         ap(m1,i,j) = ap(m1-1,i,j)
      enddo
   enddo
!endif

return
end

!     ******************************************************************

subroutine botset(m1,m2,m3,ia,iz,ja,jz,ibcon,aa,vnam)

use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ibcon,i,j
real :: dzmr
real, dimension(m1,m2,m3) :: aa
character(len=*) :: vnam

if (vnam .eq. 'P') then
   dzmr = dzm(2) / dzm(1)
   do i = 1,m2
      do j = 1,m3
         aa(1,i,j) = aa(2,i,j) + (aa(2,i,j) - aa(3,i,j)) * dzmr
      enddo
   enddo
else
   do i = 1,m2
      do j = 1,m3
         aa(1,i,j) = aa(2,i,j)
      enddo
   enddo
endif

return
end

!     ******************************************************************

subroutine dumset(m1,m2,m3,ia,iz,ja,jz,ibcon,aa,vnam)

use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ibcon,i,j,k
real, dimension(m1,m2,m3) :: aa
character(len=*) :: vnam

if (vnam .eq. 'U' .and. iand(ibcon,2) .gt. 0) then
   do j = 1,m3
      do k = 1,m1
         aa(k,m2,j) = aa(k,m2-1,j)
      enddo
   enddo
elseif (vnam .eq. 'V' .and. iand(ibcon,8) .gt. 0) then
   do i = 1,m2
      do k = 1,m1
         aa(k,i,m3) = aa(k,i,m3-jdim)
      enddo
   enddo
elseif (vnam .eq. 'W') then
   do j = 1,m3
      do i = 1,m2
         aa(m1,i,j) = aa(m1-1,i,j)
      enddo
   enddo
endif

return
end

!     *****************************************************************

subroutine rayft()

use mem_tend, only: tend
use mem_scratch, only: scratch
use mem_basic, only: basic_g
use mem_grid, only: nfpt, distim, ngrid, if_adap, grid_g, nnzp
use node_mod
use micphys

implicit none
include "i8.h"
integer :: ii,i,j,k

integer(kind=i8) :: mxyzp, ind

!     This routine is the rayleigh friction driver for the
!     theta friction and is called from the long timestep.

if (nfpt .eq. 0 .or. distim .eq. 0.) return

mxyzp = mxp * myp * mzp

!     First load past virtual theta into temporary.

if (level .ge. 1) then
   ind = 0
   do j = 1,nodemyp(mynum,ngrid)
      do i = 1,nodemxp(mynum,ngrid)
         do k = 1,nnzp(ngrid)
            ind = ind + 1
            scratch%vt3da(ind) = basic_g(ngrid)%theta(k,i,j)  &
               * (1. + .61 * basic_g(ngrid)%rv(k,i,j))
         enddo
      enddo
   enddo
else
  call atob_long(mxyzp, basic_g(ngrid)%theta(1,1,1), scratch%vt3da(1))
endif

!     Now get rayleigh friction tendency

if (if_adap == 0) then

   call rayf(4,mzp,mxp,myp,ia,iz,ja,jz,ibcon                  &
      ,scratch%vt3da      (1)    ,basic_g(ngrid)%th0 (1,1,1)  &
      ,tend%tht           (1)    ,grid_g(ngrid)%rtgt (1,1)    &
      ,grid_g(ngrid)%topt (1,1)                               )

else

   call rayf_adap(4,mzp,mxp,myp,ia,iz,ja,jz,ibcon     &
      ,grid_g(ngrid)%lpw  (1,1)   ,scratch%vt3da (1)  &
      ,basic_g(ngrid)%th0 (1,1,1) ,tend%tht      (1)  )

endif

return
end

!******************************************************************************

subroutine rayf(ifrom,m1,m2,m3,ia,iz,ja,jz,ibcon,var,th0,tht,rtgx,topx)

use mem_grid
use mem_scratch
use ref_sounding

implicit none

integer :: ifrom,m1,m2,m3,ia,iz,ja,jz,ibcon
real, dimension(m1,m2,m3) :: var,th0,tht
real, dimension(m2,m3) :: rtgx,topx

real :: zmkf,c1,c2
integer :: kf,i,j,k

!     This routine calculates rayleigh friction terms velocity and theta_il

if (nfpt .eq. 0 .or. distim .le. 0) return
kf = nnz(1) - nfpt
zmkf = zmn(kf,1)
c1 = 1. / (distim * (ztop - zmkf))
c2 = dts * c1
goto(100,200,300,400) ifrom
100   continue

!     u friction

do j = ja,jz
   do i = ia,iz
      do k = 1,nzp
         vctr2(k) = zt(k) * rtgx(i,j) + topx(i,j)
      enddo
      call htint(nzp,u01dn(1,ngrid),zt,nzp,vctr5,vctr2)
      do k = nz,2,-1
         if (vctr2(k) .le. zmkf) go to 10

         var(k,i,j) = var(k,i,j) + c2 * (vctr2(k) - zmkf)  &
            * (vctr5(k) - var(k,i,j))

      enddo
10         continue
   enddo
enddo
return
200   continue

!     V friction

if (jdim .eq. 0 .and. icorflg .eq. 0) return
do j = ja,jz
   do i = ia,iz
      do k = 1,nzp
         vctr2(k) = zt(k) * rtgx(i,j) + topx(i,j)
      enddo
      call htint(nzp,v01dn(1,ngrid),zt,nzp,vctr5,vctr2)
      do k = nz,2,-1
         if (vctr2(k) .le. zmkf) go to 20
         var(k,i,j) = var(k,i,j) + c2 * (vctr2(k) - zmkf)  &
            * (vctr5(k) - var(k,i,j))
      enddo
20         continue
   enddo
enddo
return
300   continue

!     W friction

do j = ja,jz
   do i = ia,iz
      do k = nz,2,-1
         vctr2(k) = zt(k) * rtgx(i,j) + topx(i,j)
         if (vctr2(k) .le. zmkf) go to 30
         var(k,i,j) = var(k,i,j) - c2 * (vctr2(k) - zmkf) * var(k,i,j)
      enddo
30         continue
   enddo
enddo
return
400   continue

!     THETA FRICTION

do j = ja,jz
   do i = ia,iz
      do k = nz,2,-1
         vctr2(k) = zt(k) * rtgx(i,j) + topx(i,j)
         if (vctr2(k) .le. zmkf) go to 40
         tht(k,i,j) = tht(k,i,j) + c1 * (vctr2(k) - zmkf)  &
              * (th0(k,i,j) - var(k,i,j))
      enddo
40         continue
   enddo
enddo
return
end

