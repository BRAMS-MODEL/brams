module ModPostUtils
   use ModParallelEnvironment, only : &
         MsgDump

   implicit none

   private

   public :: seaprs_0
   public :: rams_comp_tempc
   public :: rams_comp_tempk
   public :: rams_comp_dewk
   public :: rams_comp_dewk_2m
   public :: get_ZItheta
   public :: rams_comp_thetv
   public :: rams_get_surface
   public :: rams_comp_1minus
   public :: rams_comp_slpmm5
   public :: rams_comp_rh
   public :: relative_humidity_2m
   public :: rams_comp_press
   public :: rams_comp_5050
   public :: rams_comp_avgu
   public :: rams_comp_avgv
   public :: rams_comp_avgw
   public :: rams_comp_bigpatch
   public :: rams_comp_bowen
   public :: rams_comp_copysst
   public :: rams_comp_sfc_press
   public :: cape_cine
   public :: htint
   public :: rams_comp_dn0
   public :: rams_comp_z
   public :: xy_ll
   public :: ll_xy
   public :: uvtoueve
   public :: rams_comp_rotate
   public :: rams_comp_speed
   public :: rams_reduced_temp
   public :: rams_fill_sst
   public :: get_leaf_soil
   public :: rams_comp_pbl
   public :: rs
   public :: td
   public :: integrando
   public :: potencialeq
   public :: potencial
   public :: presdoncl
   public :: razaodemistura
   public :: rparcela
   public :: tparcela
   public :: tempvirtual
   public :: calccine
   public :: calccape
   public :: UpperCase
   public :: DumpFloating
   public :: DumpFixed
   public :: DumpInteger
   public :: DumpRealPairs
   public :: DumpIntegerPairs
   public :: ptransvar
   public :: ctransvar
   public :: calc_omeg
   public :: comp_vertint
   public :: comp_vertint_press
   public :: comp_slp_metar
   public :: rams_reduced_rv
   public :: calc_v10m
   public :: calc_u10m
   public :: rams_comp_dir
   public :: rams_reduced_wind
   public :: rams_comp_vegclass
   public :: calc_poda_index
   public :: copy_x_to_y
   public :: undef
   public :: checkUsingJules

   include "files.h"

   interface rams_comp_tempc
      module procedure rams_comp_tempc_2d, rams_comp_tempc_3d
   end interface

   interface rams_comp_tempk
      module procedure rams_comp_tempk_2d, rams_comp_tempk_3d
   end interface

   logical, parameter :: dumpLocal = .false.

   real, parameter :: undef = -9.99e+33
contains


   subroutine seaprs_0(t, pp, ter, sfp, ts, imx, jmx, kx, slp)
      integer, intent(in) :: imx
      integer, intent(in) :: jmx
      integer, intent(in) :: kx
      real, intent(in) :: t(:, :, :)
      real, intent(in) :: pp(:, :, :)
      real, intent(in) :: ter(:, :)
      real, intent(in) :: sfp(:, :)
      real, intent(out) :: ts(:, :)
      real, intent(out) :: slp(:, :)

      integer :: i
      integer :: j
      integer :: k
      integer :: kupto
      integer :: klo
      integer :: khi

      real, parameter :: r = 287.04
      real, parameter :: g = 9.8
      real, parameter :: gamma = 6.5e-3
      real, parameter :: tc = 273.16 + 17.5
      real, parameter :: pconst = 100.
      real :: xterm
      real :: xk
      real :: xkhold
      real :: plo
      real :: phi
      real :: tlo
      real :: thi
      real :: tl
      real :: tbar
      real :: hl
      real :: t0hold
      real :: ps(imx, jmx)
      real :: pl(imx, jmx)
      real :: t0(imx, jmx)
      real :: xklev(imx, jmx)
      !
      logical l1, l2, l3, l4
      !
      !
      !
      !
      !     ... sea level pressure
      !
      xterm = gamma * r / g
      !
      !     ... compute pressure at pconst mb above surface (pl)
      !
      kupto = kx / 2
      99  continue
      do j = 1, jmx
         do i = 1, imx
            pl(i, j) = sfp(i, j) - pconst
            xklev(i, j) = 0.
         end do
      end do
      !
      !     ... find 2 levels on sigma surfaces surrounding pl at each i,j
      !
      do j = 1, jmx
         do i = 1, imx
            do k = kx - 1, kupto, -1
               xk = float(k)
               xkhold = xklev(i, j)
               xklev(i, j) = merge(xk, xkhold, &
                     (((pp(i, j, k)).lt.pl(i, j)) .and.  &
                           ((pp(i, j, k + 1)).ge.pl(i, j))))
            end do
            if(xklev(i, j).lt.1.) then
               print *, 'error finding pressure level ', pconst, ' mb ', &
                     'above the surface'
               print *, 'last k level =', kupto
               if(kupto.ne.1) then
                  print *, 'trying again with kupto=1'
                  kupto = 1
                  goto 99
               else
                  print *, 'i,j=', i, j
                  print *, 'pl=', pl(i, j)
                  print *, 'psfc=', sfp(i, j)
                  stop
               end if
            end if
         end do
      end do
      !
      !     ... get temperature at pl (tl), extrapolate t at surface (ts)
      !         and t at sea level (t0) with 6.5 k/km lapse rate
      !
      do j = 1, jmx
         do i = 1, imx
            klo = nint(xklev(i, j)) + 1
            khi = nint(xklev(i, j))
            plo = pp(i, j, klo)
            phi = pp(i, j, khi)
            tlo = t(i, j, klo)
            thi = t(i, j, khi)
            tl = thi - (thi - tlo) * alog(pl(i, j) / phi) / alog(plo / phi)
            ts(i, j) = tl * (sfp(i, j) / pl(i, j))**xterm
            tbar = (ts(i, j) + tl) * 0.5
            hl = ter(i, j) - r / g * alog(pl(i, j) / sfp(i, j)) * tbar
            t0(i, j) = tl + gamma * hl
         end do
      end do
      !
      !     ... correct sea level temperature if too hot
      !
      do j = 1, jmx
         do i = 1, imx
            l1 = t0(i, j).lt.tc
            l2 = ts(i, j).le.tc
            l3 = .not.l1
            t0hold = t0(i, j)
            t0(i, j) = merge(t0hold, &
                  merge(tc, tc - 0.005 * (ts(i, j) - tc)**2, l2.and.l3), l1.and.l2)
         end do
      end do
      !
      !     ... compute sea level pressure
      !
      do j = 1, jmx
         do i = 1, imx
            slp(i, j) = sfp(i, j) * exp(2. * g * ter(i, j) / (r * (ts(i, j) + t0(i, j))))
         end do
      end do
   end subroutine seaprs_0


   real function rs(p, t)
      real, intent(in) :: p
      real, intent(in) :: t
      real :: es

      if (t>100.) then
         es = 610.78 * exp (17.269 * (t - 273.16) / (t - 35.86))
      else
         es = 0.
      endif
      rs = .622 * es / (p - es)
   end function rs


   real function td(p, rs)
      real, intent(in) :: p
      real, intent(in) :: rs
      real :: rr
      real :: es
      real :: esln
      rr = rs + 1e-8
      es = p * rr / (.622 + rr)
      esln = log(es)
      td = (35.86 * esln - 4947.2325) / (esln - 23.6837)
   end function td


   subroutine rams_comp_tempk_3d(a, b)
      real, intent(INOUT) :: a(:, :, :)
      real, intent(IN) :: b(:, :, :)
      integer :: i, j, k

      include "post_rconstants.h"

      do k = 1, size(a, 3)
         do j = 1, size(a, 2)
            do i = 1, size(a, 1)
               a(i, j, k) = a(i, j, k) * b(i, j, k) / cp
            enddo
         enddo
      enddo
   end subroutine rams_comp_tempk_3d


   subroutine rams_comp_tempk_2d(a, b)
      real, intent(INOUT) :: a(:, :)
      real, intent(IN) :: b(:, :)
      integer :: i, j

      include "post_rconstants.h"

      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            a(i, j) = a(i, j) * b(i, j) / cp
         enddo
      enddo
   end subroutine rams_comp_tempk_2d


   subroutine rams_comp_tempc_3d(a)
      real, intent(INOUT) :: a(:, :, :)
      integer :: i, j, k

      do k = 1, size(a, 3)
         do j = 1, size(a, 2)
            do i = 1, size(a, 1)
               a(i, j, k) = a(i, j, k) - 273.16
            enddo
         enddo
      enddo
   end subroutine rams_comp_tempc_3d


   subroutine rams_comp_tempc_2d(a)
      real, intent(INOUT) :: a(:, :)
      integer :: i, j

      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            a(i, j) = a(i, j) - 273.16
         enddo
      enddo
   end subroutine rams_comp_tempc_2d


   subroutine rams_comp_dewk(a, b, c)
      real, intent(INOUT) :: a(:, :, :)
      real, intent(IN) :: b(:, :, :)
      real, intent(IN) :: c(:, :, :)

      include "post_rconstants.h"

      integer :: i, j, k
      real :: xpress, xtemp, xwatsat

      do k = 1, size(a, 3)
         do j = 1, size(a, 2)
            do i = 1, size(a, 1)
               xpress = (b(i, j, k) / cp)**cpor * p00
               xtemp = c(i, j, k) * b(i, j, k) / cp
               xwatsat = rs(xpress, xtemp)
               a(i, j, k) = td(xpress, min (a(i, j, k), xwatsat))
            enddo
         enddo
      enddo
   end subroutine rams_comp_dewk


   subroutine rams_comp_dewk_2m(a, b, c)
      real, intent(INOUT) :: a(:, :)
      real, intent(IN) :: b(:, :)
      real, intent(IN) :: c(:, :)

      include "post_rconstants.h"

      integer :: i, j, k
      real :: xpress, xtemp, xwatsat

      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            xpress = (b(i, j) / cp)**cpor * p00
            xtemp = c(i, j)
            xwatsat = rs(xpress, xtemp)
            a(i, j) = td(xpress, min(a(i, j) / 1000., xwatsat))
         enddo
      enddo
   end subroutine rams_comp_dewk_2m

   subroutine get_ZItheta(a, c, e, ztn)
      real, intent(OUT) :: a(:, :)
      real, intent(IN) :: c(:, :, :)
      real, intent(IN) :: e(:, :, :)
      real, intent(IN) :: ztn(:) !ztn(:,currGrid) from mem_grid
      real :: dtheta, x
      integer :: kzi, i, j, k

      a = 0.

      do i = 1, size(c, 1)
         do j = 1, size(c, 2)
            do k = 3, size(c, 3) - 5
               dtheta = c(i, j, k) - c(i, j, k - 1)
               x = 0.6
               if (dtheta>x .or. abs(e(i, j, k))>1.e-5) exit
            enddo
            kzi = k
            a(i, j) = 0.5 * (ztn(kzi) + ztn(kzi - 1))

            if(a(i, j)<0. .or. kzi<=2) a(i, j) = 0.
         enddo
      enddo

   end subroutine get_ZItheta


   subroutine rams_comp_thetv(a, b)
      real, intent(inout) :: a(:, :, :)
      real, intent(in) :: b(:, :, :)
      integer :: i, j, k

      do k = 1, size(a, 3)
         do j = 1, size(a, 2)
            do i = 1, size(a, 1)
               a(i, j, k) = a(i, j, k) * (1. + .61 * b(i, j, k))
            enddo
         enddo
      enddo
   end subroutine rams_comp_thetv


   ! rams_get_surface extracts the surface from a 3D field,
   !                  storing result in a 2D field.
   !                  Assumes that surface is the second level
   !                  of the 3D field

   subroutine rams_get_surface(a, b)
      real, intent(inout) :: a(:, :)
      real, intent(in) :: b(:, :, :)
      integer :: j, k

      do k = 1, size(a, 2)
         do j = 1, size(a, 1)
            a(j, k) = b(j, k, 2)
         enddo
      enddo
   end subroutine rams_get_surface

   subroutine rams_comp_1minus(a, b)
      real, intent(inout) :: a(:, :)
      real, intent(in) :: b(:, :, :)
      integer :: i, j

      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            a(i, j) = 1 - b(i, j, 1)
         enddo
      enddo
   end subroutine rams_comp_1minus


   subroutine rams_comp_slpmm5(theta, pp, z, slp)
      real, intent(in) :: theta(:, :, :)
      real, intent(in) :: pp(:, :, :)
      real, intent(in) :: z(:, :)
      real, intent(out) :: slp(:, :)

      integer :: n1, n2, n3
      integer :: i, j, k, kk
      real :: cp
      real :: rgas
      real :: cpor
      real :: p00
      real :: sfp(size(theta, 1), size(theta, 2))
      real :: ts(size(theta, 1), size(theta, 2))
      real :: t_mm5(size(theta, 1), size(theta, 2), size(theta, 3) - 1)
      real :: p_mm5(size(theta, 1), size(theta, 2), size(theta, 3) - 1)

      n1 = size(theta, 1)
      n2 = size(theta, 2)
      n3 = size(theta, 3)

      cp = 1004
      rgas = 287
      cpor = cp / rgas
      p00 = 1.e5

      do j = 1, n2
         do i = 1, n1
            !! calculate surface pressure
            sfp(i, j) = (0.5 * (pp(i, j, 1) + pp(i, j, 2)) / cp)**cpor * p00 * .01
            !! calculate surface temp
            ts(i, j) = (0.5 / cp) * (theta(i, j, 1) * pp(i, j, 1) + &
                  theta(i, j, 2) * pp(i, j, 2))
         enddo
      enddo

      do k = 2, n3
         kk = n3 - k + 1
         do j = 1, n2
            !! flip array upside down for input to GRAPH SUBROUTINE
            do i = 1, n1
               t_mm5(i, j, kk) = theta(i, j, k) * pp(i, j, k) / cp
               p_mm5(i, j, kk) = (pp(i, j, k) / cp)**cpor * p00 * .01
            enddo
         enddo
      enddo

      call seaprs_0(t_mm5, p_mm5, z, sfp, ts, n1, n2, n3 - 1, slp)
   end subroutine rams_comp_slpmm5


   subroutine rams_comp_rh(a, b, c)
      real, intent(INOUT) :: a(:, :, :)
      real, intent(IN) :: b(:, :, :)
      real, intent(IN) :: c(:, :, :)
      integer :: i, j, k
      real :: xpress, xtemp

      include "post_rconstants.h"

      do k = 1, size(a, 3)
         do j = 1, size(a, 2)
            do i = 1, size(a, 1)
               xtemp = c(i, j, k) * b(i, j, k) / cp
               xpress = (b(i, j, k) / cp)**cpor * p00
               a(i, j, k) = 100. * min (1.                                      &
                     , max (0., a(i, j, k) / rs(xpress, xtemp)))
            enddo
         enddo
      enddo

   end subroutine rams_comp_rh


   subroutine satured_vapor_pressure(svp, temp)
      implicit none
      ! !INPUT PARAMETERS:
      !
      ! Air temperature [Celsius]
      real, intent(in) :: temp(:, :)

      ! !OUTPUT PARAMETERS:
      !
      ! saturated vapor pressure [Pa]
      real, intent(inout) :: svp(:, :)

      ! !REVISION HISTORY:
      !  19 Feb 2013 - J. G. de Mattos - Initial Version
      !
      ! !SEE ALSO:
      !     Bolton, D., The computation of equivalent potential temperature,
      !                  Monthly Weather Review, 108, 1046-1053, 1980..

      !character(len=1024),parameter :: myname_=trim(myname)//'::es'
      integer :: i, j

      do j = 1, size(temp, 2)
         do i = 1, size(temp, 1)
            svp(i, j) = 611.2 * exp((17.67 * temp(i, j)) / (temp(i, j) + 243.5))
         end do
      end do

   end subroutine satured_vapor_pressure


   subroutine relative_humidity_2m(relative_humidity, td2mj, t2mj)
      ! INPUT PARAMETERS:
      ! td2mj - Dewpoint temp at 2m - from JULES [Celsius]
      real, intent(in) :: td2mj(:, :)
      ! t2mj - Temperature at 2m - from JULES [Celsius]
      real, intent(in) :: t2mj(:, :)
      ! OUTPUT
      ! relative_humidity [percent]
      real, intent(inout) :: relative_humidity(:, :)

      integer :: i, j, k

      real, allocatable :: svp1(:, :)
      real, allocatable :: svp2(:, :)

      ! all variables must have same dimensions, td2mj, t2mj, svp*, relative_humidity
      allocate(svp1(size(td2mj, 1), size(td2mj, 2)))
      allocate(svp2(size(td2mj, 1), size(td2mj, 2)))

      call satured_vapor_pressure(svp1, td2mj)
      call satured_vapor_pressure(svp2, t2mj)

      do j = 1, size(svp1, 2)
         do i = 1, size(svp1, 1)
            relative_humidity(i, j) = (svp1(i, j) / svp2(i, j)) * 100
         end do
      end do

   end subroutine relative_humidity_2m


   subroutine rams_comp_sfc_press(a, c)
      real, intent(INOUT) :: a(:, :)
      real, intent(IN) :: c(:, :, :)
      integer :: i, j

      include "post_rconstants.h"

      do j = 1, size(c, 2)
         do i = 1, size(c, 1)
            a(i, j) = (0.5 * (c(i, j, 1) + c(i, j, 2)) / cp)**cpor * p00 * .01
         enddo
      enddo

   end subroutine rams_comp_sfc_press

   subroutine rams_comp_press(a)
      real, intent(inout) :: a(:, :, :)
      integer :: i, j, k

      include "post_rconstants.h"

      do k = 1, size(a, 3)
         do j = 1, size(a, 2)
            do i = 1, size(a, 1)
               a(i, j, k) = (a(i, j, k) / cp)**cpor * p00 * .01
            enddo
         enddo
      enddo

   end subroutine rams_comp_press


   subroutine rams_comp_5050(a, d)
      real, intent(inout) :: a(:, :)
      real, intent(in) :: d(:, :, :)
      integer :: i, j

      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            a(i, j) = .5 * (a(i, j) + d(i, j, 2))
         enddo
      enddo

   end subroutine rams_comp_5050


   subroutine rams_comp_avgu(a)
      real, intent(inout) :: a(:, :, :)
      integer :: i, j, k

      do k = 1, size(a, 3)
         do j = 1, size(a, 2)
            do i = size(a, 1), 2, -1
               a(i, j, k) = 0.5 * (a(i, j, k) + a(i - 1, j, k))
            enddo
         enddo
      enddo

   end subroutine rams_comp_avgu


   subroutine rams_comp_avgv(a)
      real, intent(inout) :: a(:, :, :)
      integer :: i, j, k

      do k = 1, size(a, 3)
         do j = size(a, 2), 2, -1
            do i = 1, size(a, 1)
               a(i, j, k) = 0.5 * (a(i, j, k) + a(i, j - 1, k))
            enddo
         enddo
      enddo

   end subroutine rams_comp_avgv


   subroutine rams_comp_avgw(a)
      real, intent(INOUT) :: a(:, :, :)
      integer :: i, j, k

      do k = size(a, 3), 2, -1
         do j = 1, size(a, 2)
            do i = 1, size(a, 1)
               a(i, j, k) = 0.5 * (a(i, j, k) + a(i, j, k - 1))
            enddo
         enddo
      enddo

   end subroutine rams_comp_avgw


   subroutine rams_comp_bigpatch(a, f, b)
      real, intent(in) :: a(:, :, :, :)
      real, intent(inout) :: f(:, :, :)
      real, intent(out) :: b(:, :, :)
      integer :: i, j, k

      ! Extract LSP value from largest patch
      do k = 1, size(a, 3)
         do j = 1, size(a, 2)
            do i = 1, size(a, 1)
               if (f(i, j, 2) >= f(i, j, 1)) then
                  b(i, j, k) = a(i, j, k, 2)
               else
                  b(i, j, k) = a(i, j, k, 1)
               endif
            enddo
         enddo
      enddo

      ! Copy b into f, which was passed in as a(1).  n3 may exceed n4 but this
      ! should be ok.
      do k = 1, min(size(f, 3), size(b, 3))
         do j = 1, size(f, 2)
            do i = 1, size(f, 1)
               f(i, j, k) = b(i, j, k)
            enddo
         enddo
      enddo

   end subroutine rams_comp_bigpatch


   subroutine rams_comp_bowen(a, b)
      real, intent(inout) :: a(:, :, :)
      real, intent(in) :: b(:, :, :)
      integer :: i, j, k, n1, n2, n3

      n3 = size(a, 3)
      n2 = size(a, 2)
      n1 = size(a, 1)
      do k = 1, n3
         do j = 1, n2
            do i = 1, n1
               a(i, j, k) = a(i, j, k) / max (1.e-12, b(i, j, k)) * 1004. / 2.5e6
            enddo
         enddo
      enddo

   end subroutine rams_comp_bowen


   subroutine rams_comp_copysst(a)
      real, intent(INOUT) :: a(:, :, :)
      integer :: i, j, k, n1, n2, n3

      n3 = size(a, 3)
      n2 = size(a, 2)
      n1 = size(a, 1)

      do k = 1, n3
         do j = 1, n2
            do i = 1, n1
               a(i, j, k) = a(i, j, n3) - 273.15
            enddo
         enddo
      enddo

   end subroutine rams_comp_copysst

   real function integrando(tvamb, tvpar, indef)
      real, intent(in) :: tvamb
      real, intent(in) :: tvpar
      real, intent(in) :: indef

      include "post_rconfig.h"

      if ((tvamb == indef) .or. (tvpar == indef)) then
         integrando = indef
      else
         integrando = -ra * (tvpar - tvamb)
      endif

   end function integrando


   real function potencialeq(pres0, temp0, rmis0, indef)
      real, intent(in) :: pres0
      real, intent(in) :: temp0
      real, intent(in) :: rmis0
      real, intent(in) :: indef
      real :: pvap, tncl, tpot

      include "post_rconfig.h"

      if ((pres0 == indef) .or. (temp0 == indef) .or. (rmis0 == indef)) then
         potencialeq = indef
      else
         pvap = pres0 * rmis0 / (epsi + rmis0)
         tncl = 2840. / (3.5 * log (temp0) - log (pvap) - 4.805) + 55.
         tpot = temp0 * ((1000. / pres0)**(0.2854 * (1 - 0.28 * rmis0)))
         potencialeq = tpot * exp ((3.376 / tncl - 0.00254) * &
               1000. * rmis0 * (1 + 0.81 * rmis0))
      endif

   end function potencialeq


   real function potencial(pres0, temp0, rmis0, indef)
      real, intent(in) :: pres0
      real, intent(in) :: temp0
      real, intent(in) :: rmis0
      real, intent(in) :: indef

      if ((pres0 == indef) .or. (temp0 == indef)) then
         potencial = indef
      elseif (rmis0 == indef) then
         potencial = temp0 * ((1000. / pres0)**0.2854)
      else
         potencial = temp0 * ((1000. / pres0)**(0.2854 * (1 - 0.28 * rmis0)))
      endif

   end function potencial


   real function presdoncl(pres0, temp0, urel0, indef)
      real, intent(in) :: pres0
      real, intent(in) :: temp0
      real, intent(in) :: urel0
      real, intent(in) :: indef
      real :: tempk, tpot, tncl

      if ((pres0 == indef) .or. (temp0 == indef) .or. (urel0 == indef)) then
         presdoncl = indef
      else
         tempk = temp0 + 273.16
         tpot = tempk * ((1000. / pres0)**0.286)
         tncl = 1 / (1 / (tempk - 55.) - log (urel0 / 100.) / 2840.) + 55.
         presdoncl = 1000. * ((tncl / tpot)**3.4965035)
      endif

   end function presdoncl


   real function razaodemistura(pres0, temp0, urel0, indef)
      real, intent(in) :: pres0
      real, intent(in) :: temp0
      real, intent(in) :: urel0
      real, intent(in) :: indef
      real :: pvap

      include "post_rconfig.h"

      if ((pres0 == indef) .or. (temp0 == indef) .or. (urel0 == indef)) then
         razaodemistura = indef
      else
         pvap = 0.01 * urel0 * 6.112 * exp (17.67 * temp0 / (temp0 + 243.5))
         razaodemistura = epsi * pvap / (pres0 - pvap)
      endif

   end function razaodemistura


   real function rparcela(pres0, pncl0, rmis0, tpar, indef)
      real, intent(in) :: pres0
      real, intent(in) :: pncl0
      real, intent(in) :: rmis0
      real, intent(in) :: tpar
      real, intent(in) :: indef

      if ((pres0 == indef) .or. (pncl0 == 0)) then
         rparcela = indef

      elseif (pres0 > pncl0) then
         rparcela = rmis0
      else
         rparcela = razaodemistura(pres0, tpar - 273.16, 0. * tpar + 100., indef)
      endif

   end function rparcela


   real function tparcela(pres0, pncl0, tpot0, tpeq0, rmis0, erro0, indef)
      real, intent(in) :: pres0
      real, intent(in) :: pncl0
      real, intent(in) :: tpot0
      real, intent(in) :: tpeq0
      real, intent(in) :: rmis0
      real, intent(in) :: erro0
      real, intent(in) :: indef
      real :: erro, tparnovo, tparm0, tparm1, esatm0, esatm1, rsatm0
      real :: rsatm1, tpotm0, tpotm1, tpeqm0, tpeqm1

      include "post_rconfig.h"

      tparnovo = 273.16
      erro = 2. * erro0

      if ((pres0 == indef) .or. (pncl0 == indef)) then
         tparnovo = indef

      elseif (pncl0 == 0.) then !Pressão DO NCL inatingível
         tparnovo = indef

      elseif (pres0 > pncl0) then !Iterage com temperatura potencial

         do while (erro > erro0)
            tparm0 = tparnovo
            tparm1 = tparm0 + 1.
            rsatm0 = 1000. * rmis0 !Só por via das dúvidas
            tpotm0 = tparm0 * ((1000. / pres0)**(0.2854 * (1 - 0.28e-3 * rsatm0)))
            tpotm1 = tparm1 * ((1000. / pres0)**(0.2854 * (1 - 0.28e-3 * rsatm0)))
            tparnovo = tparm0 + (tpot0 - tpotm0) / (tpotm1 - tpotm0)
            erro = 200. * (tparnovo - tparm0) / (tparnovo + tparm0) !O erro est� em %
         enddo

      else !Iterage com temperatura potencial Equivalente

         do while (erro > erro0)
            tparm0 = tparnovo
            tparm1 = tparm0 + 1.
            esatm0 = 6.112 * exp (17.67 * (tparm0 - 273.16) / (tparm0 - 29.66))
            esatm1 = 6.112 * exp (17.67 * (tparm1 - 273.16) / (tparm1 - 29.66))
            rsatm0 = 1000. * epsi * esatm0 / (pres0 - esatm0)
            rsatm1 = 1000. * epsi * esatm1 / (pres0 - esatm1)
            tpotm0 = tparm0 * ((1000. / pres0)**(0.2854 * (1 - 0.28e-3 * rsatm0)))
            tpotm1 = tparm1 * ((1000. / pres0)**(0.2854 * (1 - 0.28e-3 * rsatm1)))
            tpeqm0 = tpotm0 * exp ((3.376 / tparm0 - 0.00254) * rsatm0 * (1 + 0.81e-3 * rsatm0))
            tpeqm1 = tpotm1 * exp ((3.376 / tparm1 - 0.00254) * rsatm1 * (1 + 0.81e-3 * rsatm1))
            tparnovo = tparm0 + (tpeq0 - tpeqm0) / (tpeqm1 - tpeqm0)
            erro = abs (200. * (tparnovo - tparm0) / (tparnovo + tparm0)) !O erro est� em %
         enddo

      endif
      tparcela = tparnovo

   end function tparcela


   real function tempvirtual(pres0, temp0, rmis0, indef)
      real, intent(in) :: pres0
      real, intent(in) :: temp0
      real, intent(in) :: rmis0
      real, intent(in) :: indef
      real :: umes, pvap

      if ((pres0 == indef) .or. (temp0 == indef) .or. (rmis0 == indef)) then
         tempvirtual = indef
      else
         umes = rmis0 / (rmis0 + 1)
         tempvirtual = temp0 * (1 + 0.61 * umes)
      endif

   end function tempvirtual


   real function calccine(num, i0, pres, temp, urel, erro0, indef)
      integer, intent(IN) :: num
      integer, intent(IN) :: i0
      real, intent(IN) :: pres(:)
      real, intent(IN) :: temp(:)
      real, intent(IN) :: urel(:)
      real, intent(IN) :: erro0
      real, intent(IN) :: indef
      real :: pncl0, rmis0, tpot0, tpeq0, tamb, ramb, tvamb
      real :: pres1, pres2, inte1, inte2
      real :: tpar1, tpar2, rpar1, rpar2, tvpar
      real :: cine
      integer :: i
      logical :: fim

      cine = 0.
      fim = .false.
      if ((pres(i0) == indef) .or. (temp(i0) == indef) .or. (urel(i0) == indef)) then
         calccine = indef
         return
      endif

      !* Tomo os primeiros valores, para depois jogá-los aos valores velhos...
      pncl0 = presdoncl(pres(i0), temp(i0), urel(i0), indef)
      rmis0 = razaodemistura(pres(i0), temp(i0), urel(i0), indef)
      tpot0 = potencial(pres(i0), temp(i0) + 273.16, rmis0, indef)
      tpeq0 = potencialeq(pres(i0), temp(i0) + 273.16, rmis0, indef)
      pres2 = pres(i0)
      tamb = temp(i0) + 273.16
      ramb = razaodemistura(pres2, tamb - 273.16, urel(i0), indef)
      tvamb = tempvirtual(pres2, tamb, ramb, indef)
      tpar2 = temp(i0) + 273.16 !Começa com mesma temperatura DO ambiente
      rpar2 = rmis0           !Começa com mesmo rmis DO ambiente
      tvpar = tempvirtual(pres2, tpar2, rpar2, indef)
      inte2 = 0.
      i = i0 + 1

      do while ((.not. fim) .and. (i <= num))
         if (pres(i) /= indef .and.temp(i) /= indef .and.urel(i) /= indef)then

            ! Passo os valores de algumas Variáveis para o valor "velho"
            pres1 = pres2
            tpar1 = tpar2
            rpar1 = rpar2
            inte1 = inte2

            ! Recalculo estas Variáveis e calculo a contribuição para o CINE
            pres2 = pres(i)
            tamb = temp(i) + 273.16
            ramb = razaodemistura(pres2, tamb - 273.16, urel(i), indef)
            tvamb = tempvirtual(pres2, tamb, ramb, indef)
            tpar2 = tparcela(pres2, pncl0, tpot0, tpeq0, rmis0, erro0, indef)
            rpar2 = rparcela(pres2, pncl0, rmis0, tpar2, indef)
            tvpar = tempvirtual(pres2, tpar2, rpar2, indef)
            inte2 = integrando(tvamb, tvpar, indef)

            if (inte2 < 0.) then
               fim = .true.
            else
               cine = cine - 0.5 * (inte1 + inte2) * log (pres2 / pres1)
            endif
         endif
         i = i + 1
      enddo

      !   Caso tenha acabado até aqui, indefini-lo-ei, pois na REALidade ele
      ! vale infinito.....
      if (.not. fim) then
         cine = indef
      end if
      calccine = cine

   end function calccine


   real function calccape(num, i0, pres, temp, urel, erro0, indef)
      integer, intent(IN) :: num
      integer, intent(IN) :: i0
      real, intent(IN) :: pres(:)
      real, intent(IN) :: temp(:)
      real, intent(IN) :: urel(:)
      real, intent(IN) :: erro0
      real, intent(IN) :: indef
      real :: pncl0, rmis0, tpot0, tpeq0, tamb, ramb, tvamb               !Ambiente
      real :: pres1, pres2, inte1, inte2                               !Comuns
      real :: tpar1, tpar2, rpar1, rpar2, tvpar                         !Parcela
      real :: cape
      integer :: i
      logical :: fim, embaixo

      cape = 0.
      fim = .false.
      embaixo = .true.

      if ((pres(i0) == indef) .or. (temp(i0) == indef) .or. (urel(i0) == indef))then
         calccape = indef
         return
      endif

      !* Tomo os primeiros valores, para depois jogá-los aos valores velhos...
      pncl0 = presdoncl(pres(i0), temp(i0), urel(i0), indef)
      rmis0 = razaodemistura(pres(i0), temp(i0), urel(i0), indef)
      tpot0 = potencial(pres(i0), temp(i0) + 273.16, rmis0, indef)
      tpeq0 = potencialeq(pres(i0), temp(i0) + 273.16, rmis0, indef)
      pres2 = pres(i0)
      tamb = temp(i0) + 273.16
      ramb = razaodemistura(pres2, tamb - 273.16, urel(i0), indef)
      tvamb = tempvirtual(pres2, tamb, ramb, indef)
      tpar2 = temp(i0) + 273.16 !Começa com mesma temperatura DO ambiente
      rpar2 = rmis0           !Começa com mesmo rmis DO ambiente
      tvpar = tempvirtual(pres2, tpar2, rpar2, indef)
      inte2 = 0.
      i = i0 + 1

      do while ((.not. fim) .and. (i <= num))
         if ((pres(i) /= indef) .and. (temp(i) /= indef) .and. (urel(i) /= indef))then
            !* Passo os valores de algumas Variáveis para o valor "velho"
            pres1 = pres2
            tpar1 = tpar2
            rpar1 = rpar2
            inte1 = inte2

            !* Recalculo estas Variáveis e calculo a contribuição para o CINE
            pres2 = pres(i)
            tamb = temp(i) + 273.16
            ramb = razaodemistura(pres2, tamb - 273.16, urel(i), indef)
            tvamb = tempvirtual(pres2, tamb, ramb, indef)
            tpar2 = tparcela(pres2, pncl0, tpot0, tpeq0, rmis0, erro0, indef)
            rpar2 = rparcela(pres2, pncl0, rmis0, tpar2, indef)
            tvpar = tempvirtual(pres2, tpar2, rpar2, indef)
            inte2 = integrando(tvamb, tvpar, indef)

            if ((.not. embaixo) .and. (inte2 > 0.)) then
               fim = .true.
            elseif (inte2 < 0) then
               embaixo = .false.
               cape = cape + 0.5 * (inte1 + inte2) * log (pres2 / pres1)
            endif
         endif
         i = i + 1
      enddo

      !*   Caso tenha acabaDO até aqui, indefini-lo-ei, pois na REALidade ele
      !* vale infinito.....
      if (.not. fim) then
         cape = 0.
      end if
      calccape = cape

   end function calccape


   subroutine cape_cine(press, tempk, ur, dummy, nome)
      real, intent(in) :: press(:, :, :)
      real, intent(in) :: tempk(:, :, :)
      real, intent(in) :: ur(:, :, :)
      real, intent(out) :: dummy(:, :)
      character(len = *), intent(in) :: nome
      real :: cine(size(press, 1), size(press, 2))
      real :: cape(size(press, 1), size(press, 2))
      real :: press2(size(press, 3))
      real :: temp2(size(press, 3))
      real :: ur2(size(press, 3))
      character :: nomeIN*100, nomeOUT*100
      integer :: t, z, x, y
      real :: erro0, indef
      integer :: i, j, k, nx, ny, nz

      indef = -9.99e33
      nx = size(press, 1)
      ny = size(press, 2)
      nz = size(press, 3)
      erro0 = 5.e-5

      do i = 1, nx
         do j = 1, ny
            do z = 1, nz - 1
               press2(z) = PRESS(i, j, z + 1)
               temp2 (z) = TEMPK(i, j, z + 1) - 273.16
               ur2   (z) = UR(i, j, z + 1)
            enddo

            if (nome  ==  'cape') then
               cape(i, j) = calccape(nz - 1, 1, press2, temp2, ur2, erro0, indef)
               dummy(i, j) = cape(i, j)
            endif

            if (nome  ==  'cine') then
               cine(i, j) = calccine(nz - 1, 1, press2, temp2, ur2, erro0, indef)
               dummy(i, j) = cine(i, j)
            endif
         enddo
      enddo

   end subroutine cape_cine


   ! htint: vertical interpolation;
   !        given vector values vin at heights htin,
   !        interpolate vector values vectb at heights htout,
   !        with vin and htin of size nin and
   !        with vctrb and htout of size nout.
   !        htin and htout are assumed non-decreasing, with
   !        sizes >= 2

   subroutine htint(nin, vin, htin, nout, vout, htout)
      integer, intent(in) :: nin
      real, intent(in) :: vin(:)   ! (nin)
      real, intent(in) :: htin(:)  ! (nin)
      integer, intent(in) :: nout
      real, intent(out) :: vout(:)  ! (nout)
      real, intent(in) :: htout(:) ! (nout)
      integer :: in
      integer :: out
      real :: wt
      character(len = 8) :: c0, c1, c2
      character(len = 16) :: d0, d1, d2
      character(len = *), parameter :: h = "**(htint)**"

      in = 1
      do out = 1, nout
         if (htout(out) < htin(1)) then
            ! htout is bellow range of htin
            wt = (htout(out) - htin(1)) / (htin(2) - htin(1))
            vout(out) = vin(1) + (vin(2) - vin(1)) * wt
            if (dumpLocal) then
               write(c0, "(i8)") out
               write(d0, "(e15.7)") htout(out)
               write(c1, "(i8)") 1
               write(d1, "(e15.7)") htin(1)
               call MsgDump (h // &
                     " htout(" // trim(adjustl(c0)) // ") (" // trim(adjustl(d0)) // ")" // &
                     " < " // &
                     "htin(" // trim(adjustl(c1)) // ") (" // trim(adjustl(d1)) // ")")
            end if
         else if (htout(out) > htin(nin)) then
            ! htout is above range of htin
            wt = (htout(out) - htin(nin)) / (htin(nin - 1) - htin(nin))
            vout(out) = vin(nin) + (vin(nin - 1) - vin(nin)) * wt
            if (dumpLocal) then
               write(c0, "(i8)") out
               write(d0, "(e15.7)") htout(out)
               write(c1, "(i8)") nin
               write(d1, "(e15.7)") htin(nin)
               call MsgDump (h // &
                     " htout(" // trim(adjustl(c0)) // ") (" // trim(adjustl(d0)) // ")" // &
                     " > " // &
                     "htin(" // trim(adjustl(c1)) // ") (" // trim(adjustl(d1)) // ")")
            end if
         else
            ! htout is within range of htin
            ! search starts at last htin interval that matched
            ! previous htout, or at 1 in the first search
            do in = in, nin - 1
               if (htout(out)>=htin(in) .and. htout(out) <=htin(in + 1)) then
                  ! htout is at the current interval of htin
                  wt = (htout(out) - htin(in)) / (htin(in + 1) - htin(in))
                  vout(out) = vin(in) + (vin(in + 1) - vin(in)) * wt
                  if (dumpLocal) then
                     write(c0, "(i8)") out
                     write(d0, "(e15.7)") htout(out)
                     write(c1, "(i8)") in
                     write(d1, "(e15.7)") htin(in)
                     write(c2, "(i8)") in + 1
                     write(d2, "(e15.7)") htin(in + 1)
                     call MsgDump (h // &
                           " htin(" // trim(adjustl(c1)) // ") (" // trim(adjustl(d1)) // ") <= " // &
                           "htout(" // trim(adjustl(c0)) // ") (" // trim(adjustl(d0)) // ") < " // &
                           "htin(" // trim(adjustl(c2)) // ") (" // trim(adjustl(d2)) // ")")
                  end if
                  exit ! in loop
               end if
            end do
         end if
      end do
   end subroutine htint


   subroutine rams_comp_dn0(a, b, c, topt, pi01dn, th01dn, ztn, ztop, dzmn)
      real, intent(OUT) :: a(:, :, :)
      real, intent(OUT) :: b(:, :, :) ! scratch com dim. errada; (n3+1)
      real, intent(OUT) :: c(:, :, :)
      real, intent(IN) :: topt(:, :)
      real, intent(IN) :: pi01dn(:) !pi01dn(:,currGrid) from ref_sounding
      real, intent(IN) :: th01dn(:) !th01dn(:,currGrid) from ref_sounding
      real, intent(IN) :: ztn(:) !ztn(:,currGrid) from mem_grid
      real, intent(IN) :: ztop !zmn(nnzp(1)-1,currGrid) from mem_grid
      real, intent(IN) :: dzmn(:) !dzmn(:,currGrid) from mem_grid
      integer :: i, j, k, n1, n2, n3
      real :: c1, c2, c3
      real :: vctr2(size(a, 3))
      real :: vctr11(size(a, 3) + 1)
      real :: vctr12(size(a, 3) + 1)

      include "post_rconstants.h"

      n1 = size(a, 1)
      n2 = size(a, 2)
      n3 = size(a, 3)
      do j = 1, n2
         do i = 1, n1
            do k = 1, n3
               vctr2(k) = ztn(k) * (1. - topt(i, j) / ztop) + topt(i, j)
            enddo
            call htint(n3, pi01dn, ztn, n3, vctr11, vctr2)
            call htint(n3, th01dn, ztn, n3, vctr12, vctr2)

            do k = 1, n3
               b(i, j, k) = vctr12(k)
            enddo
            a(i, j, n3) = vctr11(n3)
            c1 = g * 2. * (1. - topt(i, j) / ztop)
            c2 = (1 - cpor)
            c3 = cp**c2
            do k = n3 - 1, 1, -1
               a(i, j, k) = a(i, j, k + 1) + c1 / ((b(i, j, k) + b(i, j, k + 1)) * dzmn(k))
            enddo
            do k = 1, n3
               c(i, j, k) = (c3 * p00) / (rgas * b(i, j, k) * a(i, j, k)**c2)
            enddo
         enddo
      enddo

   end subroutine rams_comp_dn0


   subroutine rams_comp_z(a, c, ztn, ztop)
      real, intent(OUT) :: a(:, :, :)
      real, intent(IN) :: c(:, :)
      real, intent(IN) :: ztn(:) ! ztn(:,currGrid) from mem_grid
      real, intent(IN) :: ztop ! zmn(nnzp(1)-1,currGrid) from mem_grid
      integer :: i, j, k

      do k = 1, size(a, 3)
         do j = 1, size(a, 2)
            do i = 1, size(a, 1)
               a(i, j, k) = c(i, j) + ztn(k) * (1. - c(i, j) / ztop)
            enddo
         enddo
      enddo

   end subroutine rams_comp_z


   subroutine xy_ll(qlat, qlon, polelat, polelon, x, y)
      real, intent(OUT) :: qlat
      real, intent(OUT) :: qlon
      real, intent(IN) :: polelat
      real, intent(IN) :: polelon
      real, intent(IN) :: x
      real, intent(IN) :: y
      real :: sinplat, sinplon
      real :: cosplat, cosplon
      real :: x3p, y3p, z3p
      real :: x3q, y3q, z3q, r3q
      real :: xq, yq, zq
      real :: d, t
      real :: alpha

      include "post_rconfig.h"

      ! Evaluate sine and cosine of latitude and longitude of pole point p.
      sinplat = sin (polelat * pi180)
      cosplat = cos (polelat * pi180)
      sinplon = sin (polelon * pi180)
      cosplon = cos (polelon * pi180)

      ! Compute (x3,y3,z3) coordinates of the pole point where the origin is the
      ! center of the earth, the z axis is the north pole, the x axis is the
      ! equator and prime meridian, and the y axis is the equator and 90 E.
      x3p = erad * cosplat * cosplon
      y3p = erad * cosplat * sinplon
      z3p = erad * sinplat

      ! Compute distance d from given point R on the polar stereographic plane
      ! to the pole point P:
      d = sqrt (x ** 2 + y ** 2)

      ! Compute angle QCP where C is the center of the Earth.  This is twice
      ! angle QAP where A is the antipodal point.  Angle QAP is the same as
      ! angle RAP:

      alpha = 2. * atan2 (d, erad2)

      ! Compute zq, the height of Q relative to the polar stereographic plane:
      zq = erad * (cos (alpha) - 1.)

      ! Compute the parameter t which is the the distance ratio AQ:AR
      t = (erad2 + zq) / erad2

      ! Compute xq and yq, the x and y coordinates of Q in polar stereographic space:
      xq = t * x
      yq = t * y

      ! Transform location of Q from (x,y,z) coordinates to (x3,y3,z3):
      x3q = x3p - xq * sinplon - yq * cosplon * sinplat  &
            + zq * cosplat * cosplon
      y3q = y3p + xq * cosplon - yq * sinplon * sinplat  &
            + zq * cosplat * sinplon
      z3q = z3p + yq * cosplat + zq * sinplat

      ! Compute the latitude and longitude of Q:
      qlon = atan2 (y3q, x3q) / pi180
      r3q = sqrt (x3q ** 2 + y3q ** 2)

      qlat = atan2 (z3q, r3q) / pi180

   end subroutine xy_ll


   subroutine ll_xy(qlat, qlon, polelat, polelon, x, y)
      real, intent(in) :: qlat
      real, intent(in) :: qlon
      real, intent(in) :: polelat
      real, intent(in) :: polelon
      real, intent(out) :: x
      real, intent(out) :: y
      real :: sinplat, sinplon, sinqlat, sinqlon
      real :: cosplat, cosplon, cosqlat, cosqlon
      real :: x3p, y3p, z3p
      real :: x3q, y3q, z3q
      real :: xq, yq, zq
      real :: t

      include "post_rconfig.h"
      !-----------------------------------------------------------

      ! Evaluate sine and cosine of latitude and longitude of pole point p and
      ! input point q.
      sinplat = sin (polelat * pi180)
      cosplat = cos (polelat * pi180)
      sinplon = sin (polelon * pi180)
      cosplon = cos (polelon * pi180)

      sinqlat = sin (qlat * pi180)
      cosqlat = cos (qlat * pi180)
      sinqlon = sin (qlon * pi180)
      cosqlon = cos (qlon * pi180)

      ! Compute (x3,y3,z3) coordinates where the origin is the center of the earth,
      ! the z axis is the north pole, the x axis is the equator and prime
      ! meridian, and the y axis is the equator and 90 E.

      ! For the pole point, these are:
      x3p = erad * cosplat * cosplon
      y3p = erad * cosplat * sinplon
      z3p = erad * sinplat

      ! For the given lat,lon point, these are:
      z3q = erad * sinqlat
      x3q = erad * cosqlat * cosqlon
      y3q = erad * cosqlat * sinqlon

      ! Transform q point from (x3,y3,z3) coordinates in the above system to
      ! polar stereographic coordinates (x,y,z):
      xq = - sinplon * (x3q - x3p) + cosplon * (y3q - y3p)
      yq = cosplat * (z3q - z3p)                                      &
            - sinplat * (cosplon * (x3q - x3p) + sinplon * (y3q - y3p))
      zq = sinplat * (z3q - z3p)                                      &
            + cosplat * (cosplon * (x3q - x3p) + sinplon * (y3q - y3p))

      ! Parametric equation for line from antipodal point at (0,0,-2 erad) to
      ! point q has the following parameter (t) value on the polar stereographic
      ! plane:
      t = erad2 / (erad2 + zq)

      ! This gives the following x and y coordinates for the projection of point q
      ! onto the polar stereographic plane:
      x = xq * t
      y = yq * t

   end subroutine ll_xy


   subroutine uvtoueve(u, v, ue, ve, qlat, qlon, polelat, polelon)
      real, intent(IN) :: u
      real, intent(IN) :: v
      real, intent(OUT) :: ue
      real, intent(OUT) :: ve
      real, intent(IN) :: qlat
      real, intent(IN) :: qlon
      real, intent(IN) :: polelat
      real, intent(IN) :: polelon
      real :: angle
      real :: x0, y0
      real :: x1, y1

      call ll_xy(qlat, qlon, polelat, polelon, x0, y0)
      call ll_xy(qlat, qlon + .1, polelat, polelon, x1, y1)

      angle = -atan2 (y1 - y0, x1 - x0)
      ue = u * cos (angle) - v * sin (angle)
      ve = u * sin (angle) + v * cos (angle)

   end subroutine uvtoueve


   subroutine rams_comp_rotate(a, b, xtn, ytn, polelat, polelon)
      real, intent(inout) :: a(:, :, :)
      real, intent(inout) :: b(:, :, :)
      real, intent(in) :: xtn(:) ! xtn(:,currGrid) from mod_grid
      real, intent(in) :: ytn(:) ! ytn(:,currGrid) from mod_grid
      real, intent(in) :: polelat ! platn(currGrid) from mem_grid
      real, intent(in) :: polelon ! plonn(currGrid) from mem_grid
      real :: qlat, qlon
      real :: u, v
      integer :: i, j, k

      do k = 1, size(a, 3)
         do j = 1, size(a, 2)
            do i = 1, size(a, 1)
               call xy_ll(qlat, qlon, polelat, polelon      &
                     , xtn(i), ytn(j))
               u = a(i, j, k)
               v = b(i, j, k)
               call uvtoueve(u, v, a(i, j, k), b(i, j, k)               &
                     , qlat, qlon, polelat, polelon)
            enddo
         enddo
      enddo

   end subroutine rams_comp_rotate

   subroutine rams_comp_speed(a, b)
      real, intent(INOUT) :: a(:, :, :)
      real, intent(IN) :: b(:, :, :)
      integer :: i, j, k

      do k = 1, size(a, 3)
         do j = 1, size(a, 2)
            do i = 1, size(a, 1)
               a(i, j, k) = sqrt (a(i, j, k)**2 + b(i, j, k)**2)
            enddo
         enddo
      enddo

   end subroutine rams_comp_speed


   subroutine rams_reduced_temp (tempnew, speed, ustar  &
         , tstar, znew, zold, zrough, patfrac  &
         , cantemp, theta, pi, topo, ztop)
      real, intent(OUT) :: tempnew(:, :)
      real, intent(IN) :: speed(:, :, :)
      real, intent(IN) :: ustar(:, :, :)
      real, intent(IN) :: tstar(:, :, :)
      real, intent(IN) :: znew
      real, intent(IN) :: zold
      real, intent(IN) :: zrough(:, :, :)
      real, intent(IN) :: patfrac(:, :, :)
      real, intent(IN) :: cantemp(:, :, :)
      real, intent(IN) :: theta(:, :, :)
      real, intent(IN) :: pi(:, :, :)
      real, intent(IN) :: topo(:, :)
      real, intent(IN) :: ztop
      integer :: i, j, np, n1, n2, n3, n4
      real :: richno, rtgt, zagl, rtemp, rtempw
      real :: z0, a2, spd, cantheta, sfcpi, fh

      include "post_rconstants.h"

      n1 = size(speed, 1)
      n2 = size(speed, 2)
      n3 = size(speed, 3)
      n4 = size(ustar, 3)

      !SRF: consistent with Louis 1981
      fh = 1.
      do j = 1, n2
         do i = 1, n1
            rtgt = 1. - topo(i, j) / ztop
            zagl = zold * rtgt
            sfcpi = .5 * (pi(i, j, 1) + pi(i, j, 2))
            rtempw = 0.
            do np = 1, n4
               z0 = zrough(i, j, np)
               if (np  ==  1) then
                  z0 = .001
               end if
               spd = max (speed(i, j, 2), .25)
               cantheta = cantemp(i, j, np) * cp / sfcpi
               !SRF
               richno = g * zagl * (theta(i, j, 2) - cantheta)                   &
                     / (.5 * (theta(i, j, 2) + cantheta) * spd**2)
               a2 = (vonk / log (znew / z0)) ** 2
               if (richno>0.) then
                  rtemp = cantheta                                         &
                        !SRF
                        + (ustar(i, j, np) * tstar(i, j, np) * fh) / (a2 * spd)    &
                              * (1. + 15. * richno * sqrt (1 + 5 * richno))
                  rtemp = min (max (rtemp, cantheta), theta(i, j, 2))
               else
                  rtemp = cantheta                                                    &
                        !SRF
                        + ((ustar(i, j, np) * tstar(i, j, np) * fh) / (a2 * spd))   &
                              / (1. - 15. * richno / (1. + 75. * a2                     &
                                    * sqrt (-znew * richno / z0)))
                  rtemp = max (min (rtemp, cantheta), theta(i, j, 2))
               endif
               rtempw = rtempw + rtemp * patfrac(i, j, np)
            enddo
            tempnew(i, j) = rtempw ! temperatura potencial
         enddo
      enddo
   end subroutine rams_reduced_temp

   subroutine rams_fill_sst(kp, a, c)
      integer, intent(in) :: kp
      real, intent(out) :: a(:, :)
      real, intent(in) :: c(:, :, :, :)
      integer :: i, j

      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            a(i, j) = (c(i, j, kp, 1) - 334000) / 4186
         enddo
      enddo

   end subroutine rams_fill_sst

   subroutine get_leaf_soil(a, a2)
      real, intent(IN) :: a(:, :, :, :)
      real, intent(OUT) :: a2(:, :, :, :)
      integer :: i, j, k, ip, kip

      do ip = 1, size(a2, 4)
         do k = 1, size(a2, 3)
            do i = 1, size(a2, 1)
               do j = 1, size(a2, 2)
                  a2(i, j, k, ip) = a(i, j, k, ip)
               enddo
            enddo
         enddo
      enddo

   end subroutine get_leaf_soil

   subroutine rams_comp_pbl(a, c, ztn, ztop)
      real, intent(inout) :: a(:, :, :)
      real, intent(in) :: c(:, :)
      real, intent(in) :: ztn(:)! ztn(:,currGrid) from mod_grid
      real, intent(in) :: ztop
      integer :: i, j, k, k2, n1, n2, n3
      real :: tkethrsh
      real :: pblht

      tkethrsh = 0.01   !tke threshold for PBL height in m2/s2

      n1 = size(a, 1)
      n2 = size(a, 2)
      n3 = size(a, 3)
      do j = 1, n2
         do i = 1, n1
            pblht = 0.
            do k = 2, n3
               pblht = ztn(k - 1) * (1. - c(i, j) / ztop)
               if(a(i, j, k) <= tkethrsh) then
                  do k2 = 1, n3
                     a(i, j, k2) = pblht
                  end do
               end if
            end do
         end do
      end do
   end subroutine rams_comp_pbl


   function UpperCase(strIn) result(strOut)

      ! converts lower case letters in "strIn" into upper case letters at strOut
      ! all remaining symbols in "strIn" are copied to strOut

      character(len = *), intent(in) :: strIn
      character(len = len(strIn)) :: strOut

      integer :: i
      integer :: charInInt
      integer, parameter :: firstLowCase = iachar("a")
      integer, parameter :: lastLowCase = iachar("z")
      integer, parameter :: firstUpperCase = iachar("A")
      character(len = *), parameter :: h = "**(UpperCase)**"

      strOut = ""
      do i = 1, len_trim(strIn)
         charInInt = iachar(strIn(i : i))
         if (charInInt>=firstLowCase .and. charInInt<=lastLowCase) then
            strOut(i : i) = achar(charInInt - firstLowCase + firstUpperCase)
         else
            strOut(i : i) = strIn(i : i)
         end if
      end do
   end function UpperCase


   subroutine DumpInteger(h, array)
      character(len = *), intent(in) :: h
      integer, intent(in) :: array(:)

      integer :: i
      integer :: lenArray
      integer :: lenLine
      character(len = 8) :: d0
      integer, parameter :: lineLength = 100

      lenLine = 0
      do i = 1, size(array, 1)
         write(d0, "(i8)") array(i)
         if (lenLine == 0) then
            call MsgDump (h // d0, .true.)
            lenLine = lenLine + len(h) + len(d0)
         else if (lenLine + len(d0) > lineLength) then
            call MsgDump (" ")
            call MsgDump (h // d0, .true.)
            lenLine = len(h) + len(d0)
         else
            call MsgDump (d0, .true.)
            lenLine = lenLine + len(d0)
         end if
      end do
      if (lenLine /= 0) then
         call MsgDump (" ")
      end if
   end subroutine DumpInteger


   subroutine DumpFloating(h, array)
      character(len = *), intent(in) :: h
      real, intent(in) :: array(:)

      integer :: i
      integer :: lenArray
      integer :: lenLine
      character(len = 16) :: d0
      integer, parameter :: lineLength = 100

      lenLine = 0
      do i = 1, size(array, 1)
         write(d0, "(e16.7)") array(i)
         if (lenLine == 0) then
            call MsgDump (h // d0, .true.)
            lenLine = lenLine + len(h) + len(d0)
         else if (lenLine + len(d0) > lineLength) then
            call MsgDump (" ")
            call MsgDump (h // d0, .true.)
            lenLine = len(h) + len(d0)
         else
            call MsgDump (d0, .true.)
            lenLine = lenLine + len(d0)
         end if
      end do
      if (lenLine /= 0) then
         call MsgDump (" ")
      end if
   end subroutine DumpFloating


   subroutine DumpFixed(h, array)
      character(len = *), intent(in) :: h
      real, intent(in) :: array(:)

      integer :: i
      integer :: lenArray
      integer :: lenLine
      character(len = 8) :: d0
      integer, parameter :: lineLength = 100

      lenLine = 0
      do i = 1, size(array, 1)
         write(d0, "(f8.2)") array(i)
         if (lenLine == 0) then
            call MsgDump (h // d0, .true.)
            lenLine = lenLine + len(h) + len(d0)
         else if (lenLine + len(d0) > lineLength) then
            call MsgDump (" ")
            call MsgDump (h // d0, .true.)
            lenLine = len(h) + len(d0)
         else
            call MsgDump (d0, .true.)
            lenLine = lenLine + len(d0)
         end if
      end do
      if (lenLine /= 0) then
         call MsgDump (" ")
      end if
   end subroutine DumpFixed
  
   subroutine DumpRealPairs(h, x, y)
     use dump, only: &
      dumpMessage

      include "constants.f90"
      character(len=*), parameter :: header="**(DumpRealPairs)**"

      character(len = *), intent(in) :: h
      real, intent(in) :: x(:)
      real, intent(in) :: y(:)

      integer :: i
      integer :: lenLine
      character(len = 18) :: d0
      integer, parameter :: lineLength = 100

      if (size(x, 1) /= size(y, 1)) then
        iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion,c_fatal, &
          "invoked with unequal array sizes")
         !call fatal_error("**(DumpRealPairs)** invoked with unequal array sizes")
      end if

      lenLine = 0
      do i = 1, size(x, 1)
         write(d0, "(' (',f7.2,',',f7.2,')')") x(i), y(i)
         if (lenLine == 0) then
            call MsgDump (h // d0, .true.)
            lenLine = lenLine + len(h) + len(d0)
         else if (lenLine + len(d0) > lineLength) then
            call MsgDump (" ")
            call MsgDump (h // d0, .true.)
            lenLine = len(h) + len(d0)
         else
            call MsgDump (d0, .true.)
            lenLine = lenLine + len(d0)
         end if
      end do
      if (lenLine /= 0) then
         call MsgDump (" ")
      end if
   end subroutine DumpRealPairs


   subroutine DumpIntegerPairs(h, x, y)
      character(len = *), intent(in) :: h
      integer, intent(in) :: x(:)
      integer, intent(in) :: y(:)

      integer :: i
      integer :: lenLine
      character(len = 12) :: d0
      integer, parameter :: lineLength = 100

      if (size(x, 1) /= size(y, 1)) then
         call fatal_error("**(DumpIntegerPairs)** invoked with unequal array sizes")
      end if

      lenLine = 0
      do i = 1, size(x, 1)
         write(d0, "(' (',i4,',',i4,')')") x(i), y(i)
         if (lenLine == 0) then
            call MsgDump (h // d0, .true.)
            lenLine = lenLine + len(h) + len(d0)
         else if (lenLine + len(d0) > lineLength) then
            call MsgDump (" ")
            call MsgDump (h // d0, .true.)
            lenLine = len(h) + len(d0)
         else
            call MsgDump (d0, .true.)
            lenLine = lenLine + len(d0)
         end if
      end do
      if (lenLine /= 0) then
         call MsgDump (" ")
      end if
   end subroutine DumpIntegerPairs


   subroutine ptransvar(OutputField, vertScaleValues, pi, ztn, topo, undef)
      include "post_rconfig.h"

      real, intent(inout) :: OutputField(:, :, :)
      real, intent(in) :: vertScaleValues(:)
      real, intent(in) :: pi(:, :, :)
      real, intent(in) :: ztn(:)
      real, intent(in) :: topo(:, :)
      real, intent(in) :: undef
      real :: vctra(size(OutputField, 3))
      real :: eleva(size(OutputField, 3))
      real :: vctrb(size(OutputField, 3))
      real :: elevb(size(OutputField, 3))
      integer :: i, j, k, kk
      integer :: nz
      integer :: nVert
      real :: var_OK   !DSM

      nz = size(OutputField, 3)
      nVert = size(vertScaleValues, 1)

      do i = 1, nVert
         elevb(nVert - i + 1) = 1004. * (vertScaleValues(i) / 1000.)**.286
      enddo

      var_OK = 0.0  !DSM
      do j = 1, size(OutputField, 2)
         do i = 1, size(OutputField, 1)
            do k = 1, nz
               kk = nz - k + 1
               vctra(kk) = OutputField(i, j, k)
               eleva(kk) = pi(i, j, k)
            enddo

            call htint(nz, vctra, eleva, nVert, vctrb, elevb)

            do k = 1, nVert
               if (elevb(k)<eleva(1) .or. elevb(k)>eleva(nz)) then
                  !DSM  OutputField(i,j,nVert-k+1) = undef

                  !DSM{
                  !Em lugares onde a topografia eh muito elevada o modelo fornece valores muito fora
                  !do esperado (p. ex. tempc=-189C em 1000mb), isso faz com que o campo plotado fique
                  !parecendo incorreto. Assim, estes valores serao substituidos por
                  !algum valor que esta fora da topografia (valor dentro do esperado).
                  if (abs(vctrb(k) - var_OK) > abs(var_OK * 0.3)) then
                     OutputField(i, j, nVert - k + 1) = var_OK
                  else
                     OutputField(i, j, nVert - k + 1) = vctrb(k)
                  end if

                  !DSM}

               else
                  OutputField(i, j, nVert - k + 1) = vctrb(k)
                  var_OK = vctrb(k)  !DSM
               end if
            enddo
         enddo
      enddo

   end subroutine ptransvar


   subroutine ctransvar(OutputField, topo, vertScaleValues, ztn, ztop, undef)
      include "post_rconfig.h"

      real, intent(inout) :: OutputField(:, :, :)
      real, intent(in) :: topo(:, :)
      real, intent(in) :: vertScaleValues(:)
      real, intent(in) :: ztn(:)
      real, intent(in) :: ztop
      real, intent(in) :: undef

      real :: vctra(size(OutputField, 3))
      real :: eleva(size(OutputField, 3))
      real :: vctrb(size(OutputField, 3))
      real :: elevb(size(OutputField, 3))
      integer :: i, j, k
      integer :: n3
      integer :: nVert
      character(len = 8) :: c0, c1, c2
      character(len = *), parameter :: h = "**(ctransvar)**"

      n3 = size(OutputField, 3)
      nVert = size(vertScaleValues, 1)

      do k = 1, nVert
         elevb(k) = vertScaleValues(k)
      enddo

      do j = 1, size(OutputField, 2)
         do i = 1, size(OutputField, 1)
            if (dumpLocal) then
               write(c0, "(i8)") i
               write(c1, "(i8)") j
               write(c2, "(i8)") n3
               call MsgDump (h // " at entrance, k, OutputField(" // &
                     trim(adjustl(c0)) // "," // trim(adjustl(c1)) // ",1:" // &
                     trim(adjustl(c2)) // ") is")
               call DumpFloating(h, OutputField(i, j, :))
            end if
            do k = 1, n3
               vctra(k) = OutputField(i, j, k)
               eleva(k) = topo(i, j) + ztn(k) * (1. - topo(i, j) / ztop)
            enddo
            call htint(n3, vctra, eleva, nVert, vctrb, elevb)
            do k = 1, nVert
               if (elevb(k)<eleva(1) .or. elevb(k)>eleva(n3)) then
                  OutputField(i, j, k) = undef
               else
                  OutputField(i, j, k) = vctrb(k)
               endif
            enddo
            if (dumpLocal) then
               write(c0, "(i8)") i
               write(c1, "(i8)") j
               write(c2, "(i8)") n3
               call MsgDump (h // " at exit, k, OutputField(" // &
                     trim(adjustl(c0)) // "," // trim(adjustl(c1)) // ",1:" // &
                     trim(adjustl(c2)) // ") is")
               call DumpFloating(h, OutputField(i, j, 1 : nVert))
            end if
         enddo
      enddo
   end subroutine ctransvar


   subroutine calc_omeg(b, a, c)
      real, intent(OUT) :: b(:, :, :)
      real, intent(IN) :: a(:, :, :)
      real, intent(IN) :: c(:, :, :)
      integer :: i, j, k

      do k = 1, size(a, 3)
         do j = 1, size(a, 2)
            do i = 1, size(a, 1)
               b(i, j, k) = -a(i, j, k) * 9.8 * c(i, j, k)
            enddo
         enddo
      enddo
   end subroutine calc_omeg


   subroutine comp_vertint (b, a, topt, ztop, zmn)
      real, intent(OUT) :: b(:, :, :)
      real, intent(IN) :: a(:, :, :)
      real, intent(IN) :: zmn(:)
      real, intent(IN) :: topt(:, :)
      real, intent(IN) :: ztop
      real :: rtgt
      integer :: i, j, k, altmax

      altmax = 100000
      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            rtgt = 1. - topt(i, j) / ztop
            b(i, j, 1) = 0.
            do k = 2, size(a, 3) - 1
               if(zmn(k) > altmax) cycle
               b(i, j, 2) = b(i, j, 2) + a(i, j, k) * (zmn(k) - zmn(k - 1)) * rtgt
            enddo
         enddo
      enddo
   end subroutine comp_vertint

   !---RB
   subroutine comp_vertint_press (b, a, c, topt, ztop, zmn)
      real, intent(OUT) :: b(:, :)
      real, intent(IN) :: a(:, :, :)
      real, intent(IN) :: c(:, :, :)
      real, intent(IN) :: zmn(:)
      real, intent(IN) :: topt(:, :)
      real, intent(IN) :: ztop
      real :: rtgt
      integer :: i, j, k, altmax, press_min

      press_min = 300
      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            rtgt = 1. - topt(i, j) / ztop
            b(i, j) = 0.
            do k = 2, size(a, 3) - 1
               if(c(i, j, k) >= press_min) &
                     b(i, j) = b(i, j) + a(i, j, k) * (zmn(k) - zmn(k - 1)) * rtgt
            enddo
         enddo
      enddo
   end subroutine comp_vertint_press

   !---RB END

   subroutine comp_slp_metar(c, b, topt, a, ztn)
      real, intent(IN) :: a(:, :, :)
      real, intent(IN) :: b(:, :, :)
      real, intent(OUT) :: c(:, :)
      real, intent(IN) :: topt(:, :)
      real, intent(IN) :: ztn(:)
      real :: g = 9.8, R = 287.04, altura, Tm
      integer :: i, j

      do j = 1, size(a, 2)
         do i = 1, size(a, 1)

            altura = ztn(2) + topt(i, j)
            Tm = 288.15 - 0.0065 * altura / 2
            c(i, j) = b(i, j, 2) / exp(-(g * altura) / (R * Tm))

         enddo
      enddo
   end subroutine comp_slp_metar


   subroutine rams_reduced_rv (rvnew, speed, ustar  &
         , rstar, znew, zold, zrough, patfrac  &
         , canrvap, rv, pi, topo, ztop, cantemp, theta)

      integer :: i, j, np
      real, intent(OUT) :: rvnew(:, :)
      real, intent(IN) :: speed(:, :, :)
      real, intent(IN) :: ustar(:, :, :)
      real, intent(IN) :: rstar(:, :, :)
      real, intent(IN) :: znew
      real, intent(IN) :: zold
      real, intent(IN) :: zrough(:, :, :)
      real, intent(IN) :: patfrac(:, :, :)
      real, intent(IN) :: canrvap(:, :, :)
      real, intent(IN) :: rv(:, :, :)
      real, intent(IN) :: pi(:, :, :)
      real, intent(IN) :: topo(:, :)
      real, intent(IN) :: ztop
      real, intent(IN) :: cantemp(:, :, :)
      real, intent(IN) :: theta(:, :, :)
      real :: richno, rtgt, zagl, rrv, rrvw
      real :: z0, a2, spd, cantheta, canrv, sfcpi, fh

      include "post_rconstants.h"

      fh = 1.
      do j = 1, size(speed, 2)
         do i = 1, size(speed, 1)
            rtgt = 1. - topo(i, j) / ztop
            zagl = zold * rtgt
            sfcpi = .5 * (pi(i, j, 1) + pi(i, j, 2))
            rrvw = 0.

            do np = 1, size(ustar, 3)
               z0 = zrough(i, j, np)
               if (np  ==  1) then
                  z0 = .001
               end if
               spd = max(speed(i, j, 2), .25)
               canrv = canrvap(i, j, np)
               cantheta = cantemp(i, j, np) * cp / sfcpi
               richno = g * zagl * (theta(i, j, 2) - cantheta) &
                     / (.5 * (theta(i, j, 2) + cantheta) * spd**2)
               a2 = (vonk / log (znew / z0)) ** 2

               if (richno>0.) then
                  rrv = canrv &
                        + (ustar(i, j, np) * rstar(i, j, np) * fh) / (a2 * spd) &
                              * (1. + 15. * richno * sqrt (1 + 5 * richno))
                  rrv = min (max (rrv, canrv), rv(i, j, 2))
               else
                  rrv = canrv &
                        + ((ustar(i, j, np) * rstar(i, j, np) * fh) / (a2 * spd))   &
                              / (1. - 15. * richno / (1. + 75. * a2 &
                                    * sqrt (-znew * richno / z0)))
                  rrv = max (min (rrv, canrv), rv(i, j, 2))
               endif
               rrvw = rrvw + rrv * patfrac(i, j, np)

            enddo
            rvnew(i, j) = rrvw

         enddo
      enddo
   end subroutine rams_reduced_rv


   subroutine rams_reduced_wind (velnew, speed, ustar, znew, zold, zrough, patfrac  &
         , cantemp, theta, pi, topo, ztop)

      real, intent(OUT) :: velnew(:, :)
      real, intent(IN) :: speed(:, :, :)
      real, intent(IN) :: ustar(:, :, :)
      real, intent(IN) :: znew
      real, intent(IN) :: zold
      real, intent(IN) :: zrough(:, :, :)
      real, intent(IN) :: patfrac(:, :, :)
      real, intent(IN) :: cantemp(:, :, :)
      real, intent(IN) :: theta(:, :, :)
      real, intent(IN) :: pi(:, :, :)
      real, intent(IN) :: topo(:, :)
      real, intent(IN) :: ztop
      integer :: i, j, np
      real :: richno, rtgt, zagl, rwind, rwindw
      real :: z0, a2, spd, cantheta, sfcpi

      include "post_rconstants.h"

      do j = 1, size(speed, 2)
         do i = 1, size(speed, 1)
            rtgt = 1. - topo(i, j) / ztop
            zagl = zold * rtgt
            sfcpi = .5 * (pi(i, j, 1) + pi(i, j, 2))
            rwindw = 0.
            do np = 1, size(ustar, 3)
               z0 = zrough(i, j, np)
               if (np  ==  1) then
                  z0 = .001
               end if
               spd = max (speed(i, j, 2), .25)
               cantheta = cantemp(i, j, np) * cp / sfcpi
               richno = g * zagl * (theta(i, j, 2) - cantheta) / (theta(i, j, 2) * spd**2)
               a2 = (vonk / log (znew / z0)) ** 2
               if (richno>0.) then
                  rwind = sqrt(ustar(i, j, np)**2 / a2 * (1. + 10. * richno / sqrt(1 + 5 * richno)))
               else
                  rwind = sqrt(ustar(i, j, np)**2 / a2 / (1. - 10. * richno / (1. + 75. * a2 * sqrt(-znew * richno / z0))))
               endif
               rwind = max(min(rwind, speed(i, j, 2)), 0.)
               rwindw = rwindw + rwind * patfrac(i, j, np)

            enddo
            velnew(i, j) = rwindw
         enddo
      enddo
   end subroutine rams_reduced_wind


   subroutine rams_comp_dir(a, b, xtn, ytn, polelat, polelon)
      real, intent(INOUT) :: a(:, :, :)
      real, intent(INOUT) :: b(:, :, :)
      real, intent(IN) :: xtn(:)
      real, intent(IN) :: ytn(:)
      real, intent(IN) :: polelat
      real, intent(IN) :: polelon
      real :: qlat, qlon
      real :: u, v, ff
      integer :: i, j, k

      do k = 1, size(a, 3)
         do j = 1, size(a, 2)
            do i = 1, size(a, 1)
               call xy_ll(qlat, qlon, polelat, polelon, xtn(i), ytn(j))
               u = a(i, j, k)
               v = b(i, j, k)

               call uvtoueve(u, v, a(i, j, k), b(i, j, k), qlat, qlon, polelat, polelon)

               call winddf(a(i, j, k), ff, a(i, j, k), b(i, j, k))

            enddo
         enddo
      enddo
   end subroutine rams_comp_dir


   subroutine calc_u10m(a, c)
      real, intent(INOUT) :: a(:, :)
      real, intent(INOUT) :: c(:, :, :)
      integer :: i, j

      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            a(i, j) = a(i, j) * cos((270.0 - c(i, j, 2)) * 3.1415 / 180.0)
         enddo
      enddo
   end subroutine calc_u10m


   subroutine calc_v10m(a, c)
      real, intent(INOUT) :: a(:, :)
      real, intent(INOUT) :: c(:, :, :)
      integer :: i, j

      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            a(i, j) = a(i, j) * sin((270.0 - c(i, j, 2)) * 3.1415 / 180.0)
         enddo
      enddo
   end subroutine calc_v10m

   subroutine rams_comp_vegclass(a)
      real, intent(INOUT) :: a(:, :, :)
      integer :: i

      do i = 1, size(a, 1)
         a(i, 1, 1) = float(int(a(i, 1, 1) + .1))
      enddo

   end subroutine rams_comp_vegclass

   subroutine calc_poda_index(n1, n2, n3, c, f, g, a)
      implicit none
      integer :: n1, n2, n3, i, j
      real :: a(n1, n2, n3), c(n1, n2, n3), f(n1, n2, n3), g(n1, n2, n3)

      real x, e
      ! g= dew point T in K
      ! f= sl pressure in hPa
      ! c=  temp in K

      !print*,"max.minf",maxval(f),minval(f)
      !print*,"max.ming",maxval(g),minval(g)
      !print*,"max.minc",maxval(c),minval(c)

      do j = 2, n2 - 1
         do i = 2, n1 - 1
            g(i, j, 1) = g(i, j, 1) - 273.15
            x = 235. * log(f(i, j, 1) / 1.333224) - 1936.4 - g(i, j, 1) * 25.34 + &
                  g(i, j, 1) * log(f(i, j, 1) / 1.333224)

            e = exp(- x / (g(i, j, 1) + 235.) + log(f(i, j, 1) / 1.333224) - log(622.0))

            a(i, j, 1) = 80.51 * f(i, j, 1) / c(i, j, 1) * (1. - e / f(i, j, 1))
            a(i, j, 1) = max(0., a(i, j, 1))
         enddo
      enddo
      a(1, :, 1) = a(2, :, 1)
      a(n1, :, 1) = a(n1 - 1, :, 1)
      a(:, 1, 1) = a(:, 2, 1)
      a(:, n2, 1) = a(:, n2 - 1, 1)

      return

   end subroutine calc_poda_index

   subroutine copy_x_to_y(n1, n2, n3, c, a)
      integer :: n1, n2, n3, i, j, k
      real, intent(IN) :: c(n1, n2, n3)
      real, intent(INOUT) :: a(n1, n2, n3)

      do k = 1, n3
         do j = 1, n2
            do i = 1, n1
               a(i, j, k) = c(i, j, k)
            enddo
         enddo
      enddo
   end subroutine copy_x_to_y

   subroutine checkUsingJules(varName)
      use mem_leaf, only: &
         ISFCL

      character(len=*), intent(in) :: varName
      if (.not. ISFCL == 5) then
         !TODO call fatal_error("Variable requires JULES (ISFCL=5): ") - gera erro não identificado
         print*, 'Variable requires JULES (ISFCL=5): ', varName
         stop
      end if
   end subroutine checkUsingJules

end module ModPostUtils
