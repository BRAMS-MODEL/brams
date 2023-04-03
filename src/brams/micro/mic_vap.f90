!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine thrmstr(m1,k1,k2,lpw,pp,thp,theta,pi0,rtp,rv,i,j)

use rconstants
use micphys

implicit none

integer :: m1,i,j,k,lcat,lpw
real :: fracliq,tcoal,tairstr
integer, dimension(10) :: k1,k2
real, dimension(m1) :: pp,thp,theta,pi0,rtp,rv

do k = lpw,m1
   pitot(k) = pi0(k) + pp(k)
   press(k) = p00 * (pitot(k) * cpi) ** cpor
   tair(k) = theta(k) * pitot(k) * cpi
enddo

do k = 1,k1(10)-1
   theta(k) = thp(k)
   rv(k) = rtp(k)
enddo

do k = k2(10)+1,m1
   theta(k) = thp(k)
   rv(k) = rtp(k)
enddo

do k = k1(10),k2(10)
   til(k) = thp(k) * pitot(k) * cpi
   rliq(k) = 0.
   rice(k) = 0.
enddo

do lcat = 1,2
   do k = k1(lcat),k2(lcat)
      rliq(k) = rliq(k) + rx(k,lcat)
   enddo
enddo

do lcat = 3,5
   do k = k1(lcat),k2(lcat)
      rice(k) = rice(k) + rx(k,lcat)
   enddo
enddo

do lcat = 6,7
   do k = k1(lcat),k2(lcat)
      call qtc(qx(k,lcat),tcoal,fracliq)
      rliq(k) = rliq(k) + rx(k,lcat) * fracliq
      rice(k) = rice(k) + rx(k,lcat) * (1. - fracliq)
   enddo
enddo

do k = k1(10),k2(10)
   qhydm(k) = alvl * rliq(k) + alvi * rice(k)
   rvstr(k) = rtp(k) - rliq(k) - rice(k)
   sa(k,1) = til(k) * qhydm(k) / (1.e-12 + rliq(k) + rice(k))
enddo

do k = k1(10),k2(10)
   if (tair(k) .gt. 253.) then
      tairstr = 0.5 * (til(k)  &
         + sqrt(til(k) * (til(k) + cpi4 * qhydm(k))))
      sa(k,1) = sa(k,1) * cpi / (2. * tairstr - til(k))
   else
      tairstr = til(k) * (1. + qhydm(k) * cp253i)
      sa(k,1) = sa(k,1) * cp253i
  endif
  tairstrc(k) = tairstr - 273.16
  !LFR
  IF(tairstr<100.0) THEN
     PRINT *,'tairstr: ',tairstr,k,i,j
     print *,'til  qhydm ='  ,til(k) , qhydm(k) 
     stop 'mic_vap.f90 routine thrmstr'
  END IF
enddo

return
end

!******************************************************************************

subroutine diffprep(m1,lcat,k1,k2,rv,dn0,i,j,mynum)

use rconstants
use micphys

implicit none

integer :: m1,lcat,k1,k2,i,j,k,mynum,if1,if4,if6,if8,lhcat
real :: fre,scdei
real, dimension(m1) :: rv,dn0

if (lcat .le. 2) then
   if1 = 1
   if4 = 4
   if6 = 6
   if8 = 8
else
   if1 = 2
   if4 = 5
   if6 = 7
   if8 = 9
endif

do k = k1,k2
   lhcat = jhcat(k,lcat)

   if (rx(k,lcat) .lt. 1.e-9) go to 229

   fre = frefac1(lhcat) * emb(k,lcat) ** pwmasi(lhcat)  &
      + rdynvsci(k) * frefac2(lhcat) * emb(k,lcat) ** cdp1(lhcat)

   sb(k,lcat) = cx(k,lcat) * dn0(k) * fre * pi4dt
   su(k,lcat) = vapdif(k) * sb(k,lcat)
   sd(k,lcat) = sh(k,lcat) * rx(k,lcat)
   se(k,lcat) = su(k,lcat) * sa(k,if6) + sb(k,lcat) * thrmcon(k)
   sf(k,lcat) = su(k,lcat) * sl(if1) - sb(k,lcat) * sa(k,2)
   sg(k,lcat) = su(k,lcat) * sa(k,if8) + sb(k,lcat) * sa(k,3)  &
              + sj(lcat) * qr(k,lcat)
!     + lambda_j [Joules/kg_air added by radiative heating this timestep]
   scdei = 1. / (sc(if1) * sd(k,lcat) + se(k,lcat))
   ss(k,lcat) = sf(k,lcat) * scdei
   sw(k,lcat) = (sg(k,lcat) - sk(if1) * sd(k,lcat)) * scdei
   ttest(k,lcat) = ss(k,lcat) * rv(k) + sw(k,lcat)

229    continue

enddo

if (lcat .ge. 3 .and. lcat .le. 5) then
   do k = k1,k2
      if (rx(k,lcat) .lt. 1.e-9) go to 228
      if (ttest(k,lcat) .ge. 0.) then
         sm(k,lcat) = 0.
         sh(k,lcat) = 1.
         sd(k,lcat) = sh(k,lcat) * rx(k,lcat)
         scdei = 1. / (sc(if1) * sd(k,lcat) + se(k,lcat))
         ss(k,lcat) = sf(k,lcat) * scdei
         sw(k,lcat) = (sg(k,lcat) - sk(if1) * sd(k,lcat)) * scdei
      else
         sm(k,lcat) = 1.
      endif
228        continue
   enddo
endif

if (lcat .ge. 6) then
   do k = k1,k2
      if (rx(k,lcat) .lt. 1.e-9) go to 227
      if (ttest(k,lcat) .ge. 0.) then
         sm(k,lcat) = 0.
      else
         sm(k,lcat) = 1.
      endif
227        continue
   enddo
endif

do k = k1,k2
   if (rx(k,lcat) .lt. 1.e-9) go to 226
   sy(k,lcat) = rvsrefp(k,if1) * sm(k,lcat) * sw(k,lcat) - sa(k,if4)
   sz(k,lcat) = 1. - rvsrefp(k,if1) * ss(k,lcat) * sm(k,lcat)
   sumuy(k) = sumuy(k) + su(k,lcat) * sy(k,lcat)
   sumuz(k) = sumuz(k) + su(k,lcat) * sz(k,lcat)

226      continue
enddo

return
end

!******************************************************************************

subroutine vapdiff (m1,kf1,kf2,rv,i,j,mynum)

use micphys

implicit none

integer :: m1,kf1,kf2,i,j,k,mynum
real, dimension(m1) :: rv

do k = kf1,kf2
   rv(k) = (rvstr(k) + sumuy(k)) / (1.0 + sumuz(k))
enddo

return
end

!******************************************************************************

subroutine vapflux(m1,lcat,i,j,mynum,k1,k2,dn0,rv)

use micphys

implicit none

integer :: m1,lcat,i,j,k,mynum,k1,k2,if1,if4
real :: rxx
real, dimension(m1) :: dn0,rv

rxx=0.0 !LFR

if (lcat .le. 2) then
   if1 = 1
   if4 = 4
else
   if1 = 2
   if4 = 5
endif

do k = k1,k2

   if (rx(k,lcat) .lt. 1.e-9) go to 229

   tx(k,lcat) = (ss(k,lcat) * rv(k) + sw(k,lcat)) * sm(k,lcat)
   vap(k,lcat) = su(k,lcat) * (rv(k) + sa(k,if4) - rvsrefp(k,if1) * tx(k,lcat))

   if (vap(k,lcat) .gt. -rx(k,lcat)) then

      rxx = rx(k,lcat) + vap(k,lcat)

      if (sm(k,lcat) .gt. .5) then
         qx(k,lcat) = sc(if1) * tx(k,lcat) + sk(if1)
         qr(k,lcat) = qx(k,lcat) * rxx
      else
         qx(k,lcat) = (rv(k) * sf(k,lcat) + sg(k,lcat)  &
                    - tx(k,lcat) * se(k,lcat)) / sd(k,lcat)
         qx(k,lcat) = min(350000.,max(-100000.,qx(k,lcat)))
         qr(k,lcat) = qx(k,lcat) * rxx
      endif

   endif

!bob Now also do the following section if pristine ice totally melts:
! evaporate it too.

   if ((lcat .eq. 3 .and. qx(k,lcat) .gt. 330000.) .or.  &
      vap(k,lcat) .le. -rx(k,lcat)) then

      sumuy(k) = sumuy(k) - su(k,lcat) * sy(k,lcat)
      sumuz(k) = sumuz(k) - su(k,lcat) * sz(k,lcat)
      sumvr(k) = sumvr(k) + rx(k,lcat)
      rv(k) = (rvstr(k) + sumuy(k) + sumvr(k)) / (1.0 + sumuz(k))

      vap(k,lcat) = - rx(k,lcat)
      tx(k,lcat) = 0.
      rx(k,lcat) = 0.
      qx(k,lcat) = 0.
      qr(k,lcat) = 0.
   else
      rx(k,lcat) = rxx
   endif

229     continue

enddo
return
end

!******************************************************************************

subroutine psxfer(m1,k1,k2,dn0,i,j)

use micphys

implicit none

integer :: m1,k1,k2,i,j,k,lhcat,it
real :: embx,dn,xlim,dvap,dqr,dnum
real, dimension(m1) :: dn0

do k = k1,k2

   if (vap(k,3) .gt. 0. .or. vap(k,4) .lt. 0.) then

      if (vap(k,3) .gt. 0.) then
         lhcat = jhcat(k,3)
         embx = max(1.e-9,rx(k,3)) / max(1.e-3,cx(k,3))
         dn = dnfac(lhcat) * embx ** pwmasi(lhcat)
         it = nint(dn * 1.e6)

!srf
         it=min(it, 5000)   
         xlim = gam(it,3) * dps2 * (dps / dn) ** (gnu(3) - 1.)  &
            / (gamn1(3) * pwmas(lhcat) * dn ** 2)

         dvap = min(rx(k,3),  &
                    vap(k,3) * (xlim + gam(it,1) / gamn1(3)))
         dqr = dvap * qx(k,3)
         dnum = dvap * min(dpsmi(lhcat),1./embx)
      else
         lhcat = jhcat(k,4)
         embx = max(1.e-9,rx(k,4)) / max(1.e-3,cx(k,4))
         dn = dnfac(lhcat) * embx ** pwmasi(lhcat)
         it = nint(dn * 1.e6)
!srf
         it=min(it, 5000)   
         xlim = gam(it,3) * dps2 * (dps / dn) ** (gnu(4) - 1.)  &
            / (gamn1(4) * pwmas(lhcat) * dn ** 2)

         dvap = max(-rx(k,4),vap(k,4) * xlim)
         dqr = dvap * qx(k,4)
         dnum = dvap * max(dpsmi(lhcat),1./embx)
      endif

      rx(k,3) = rx(k,3) - dvap
      cx(k,3) = cx(k,3) - dnum
      qr(k,3) = qr(k,3) - dqr
      rx(k,4) = rx(k,4) + dvap
      cx(k,4) = cx(k,4) + dnum
      qr(k,4) = qr(k,4) + dqr

   endif
enddo
return
end

!******************************************************************************

subroutine newtemp(m1,kf1,kf2,rv,theta,i,j)

use rconstants
use micphys

implicit none

real rslf,rsif

integer :: m1,kf1,kf2,i,j,k
real, dimension(m1) :: rv,theta

do k = kf1,kf2
   tairc(k) = tairstrc(k) + sa(k,1) * (rvstr(k) - rv(k))
   tair(k)  = tairc(k) + 273.16
   theta(k) = tair(k) * cp / pitot(k)

   rvlsair(k) = rslf(press(k),tair(k))
   rvisair(k) = rsif (press(k),tair(k))
   !LFR
   IF(theta(k)<100.0) THEN
      PRINT *,theta(k)
    END IF
enddo

return
end
