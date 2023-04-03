!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

!ALF
subroutine azero_l(nlong, along)
  implicit none
  include "i8.h"
  integer(kind=i8), intent(in) :: nlong
  real, intent(inout) :: along(nlong)

  along = 0.

end subroutine azero_l

!****************************************

!!$subroutine azerov(n1)
!!$
!!$implicit none
!!$integer :: n,n1
!!$real :: a1(n1),a2(n1),a3(n1),a4(n1),a5(n1)
!!$
!!$entry azero(n1,a1)
!!$   do n=1,n1
!!$      a1(n)=0.
!!$   enddo
!!$return
!!$
!!$entry azero2(n1,a1,a2)
!!$   do n=1,n1
!!$      a1(n)=0.
!!$      a2(n)=0.
!!$   enddo
!!$return
!!$entry azero3(n1,a1,a2,a3)
!!$   do n=1,n1
!!$      a1(n)=0.
!!$      a2(n)=0.
!!$      a3(n)=0.
!!$   enddo
!!$return
!!$entry azero4(n1,a1,a2,a3,a4)
!!$   do n=1,n1
!!$      a1(n)=0.
!!$      a2(n)=0.
!!$      a3(n)=0.
!!$      a4(n)=0.
!!$   enddo
!!$return
!!$entry azero5(n1,a1,a2,a3,a4,a5)
!!$   do n=1,n1
!!$      a1(n)=0.
!!$      a2(n)=0.
!!$      a3(n)=0.
!!$      a4(n)=0.
!!$      a5(n)=0.
!!$   enddo
!!$return
!!$end subroutine azerov

!!$subroutine ae1t0(n1,a,b,c)
!!$implicit none
!!$integer :: n1
!!$real, dimension(n1) :: a,b
!!$real :: c
!!$
!!$integer :: n
!!$
!!$do n=1,n1
!!$   a(n)=b(n)*c
!!$enddo
!!$
!!$return
!!$end
subroutine ae1t0_l(n1, a, b, c)
  implicit none
  include "i8.h"
  integer(kind=i8), intent(IN) :: n1
  real, intent(INOUT) :: a(n1)
  real, intent(IN)    :: b(n1)
  real, intent(IN)    :: c

  a(:) = b(:)*c

end subroutine ae1t0_l

!!$subroutine ae1p1(npts,a,b,c)
!!$implicit none
!!$integer :: npts
!!$real :: a(npts),b(npts),c(npts)
!!$integer :: i
!!$do i=1,npts
!!$  a(i)=b(i)+c(i)
!!$enddo
!!$return
!!$end
subroutine ae1p1_l(npts,a,b,c)
  implicit none
  ! Arguments:
  include "i8.h"
  integer(kind=i8), intent(IN) :: npts
  real, intent(INOUT)          :: a(npts)
  real, intent(IN)             :: b(npts), c(npts)

  a(:) = b(:) + c(:)
end subroutine ae1p1_l
!
!     ******************************************************************
!
!!$subroutine ae1m1(npts,a,b,c)
!!$implicit none
!!$integer :: npts
!!$real :: a(npts),b(npts),c(npts)
!!$integer :: i
!!$do i=1,npts
!!$  a(i)=b(i)-c(i)
!!$enddo
!!$return
!!$end
!!$!
subroutine ae1m1_l(npts, a, b, c)
  implicit none
  include "i8.h"
  integer(kind=i8), intent(IN) :: npts
  real, intent(OUT)            :: a(npts)
  real, intent(IN)             :: b(npts)
  real, intent(IN)             :: c(npts)
  a(:) = b(:) - c(:)
end subroutine ae1m1_l
!     ******************************************************************
!!$!
!!$subroutine aen1(npts,a,b)
!!$implicit none
!!$integer :: npts
!!$real :: a(npts),b(npts)
!!$integer :: i
!!$do i=1,npts
!!$  a(i)=-b(i)
!!$enddo
!!$return
!!$end
!!$!
!     ******************************************************************
!!$!
!!$subroutine ae1t1(npts,a,b,c)
!!$implicit none
!!$integer :: npts
!!$real :: a(npts),b(npts),c(npts)
!!$integer :: i
!!$do i=1,npts
!!$  a(i)=b(i)*c(i)
!!$enddo
!!$return
!!$end
!!$!
!     ******************************************************************
!
!!$subroutine ae1(npts,a,b)
!!$implicit none
!!$integer :: npts
!!$real :: a(npts),b(npts)
!!$integer :: i
!!$do i=1,npts
!!$  a(i)=b(i)
!!$enddo
!!$return
!!$end
!
subroutine ae1_l(npts, a, b)
  implicit none
  include "i8.h"
  integer(kind=i8), intent(IN) :: npts
  real, intent(OUT)            :: a(npts)
  real, intent(IN)             :: b(npts)
  a = b
end subroutine ae1_l
!     ******************************************************************
!!$!
!!$subroutine ae1tn1(npts,a,b,c)
!!$implicit none
!!$integer :: npts
!!$real :: a(npts),b(npts),c(npts)
!!$integer :: i
!!$do i=1,npts
!!$  a(i)=-b(i)*c(i)
!!$enddo
!!$return
!!$end
!!$!
!!$!     ******************************************************************
!!$!
!!$subroutine ae1t0p1(npts,a,b,c,d)
!!$implicit none
!!$integer :: npts
!!$real :: a(npts),b(npts),c,d(npts)
!!$integer :: i
!!$do i=1,npts
!!$  a(i)=b(i)*c+d(i)
!!$enddo
!!$return
!!$end
!!$!
!!$!     ******************************************************************
!!$!
!!$subroutine ae3(n1,n2,n3,i1,i2,j1,j2,k1,k2,a,b)
!!$implicit none
!!$integer :: n1,n2,n3,i1,i2,j1,j2,k1,k2
!!$real :: a(n1,n2,n3),b(n1,n2,n3)
!!$integer :: i,j,k
!!$do j=j1,j2
!!$  do i=i1,i2
!!$    do k=k1,k2
!!$      a(k,i,j)=b(k,i,j)
!!$    enddo
!!$  enddo
!!$enddo
!!$return
!!$end
!!$!
!!$!     ******************************************************************
!!$!
!!$subroutine ae3p3(n1,n2,n3,i1,i2,j1,j2,k1,k2,a,b,c)
!!$implicit none
!!$integer :: n1,n2,n3,i1,i2,j1,j2,k1,k2
!!$real :: a(n1,n2,n3),b(n1,n2,n3),c(n1,n2,n3)
!!$integer :: i,j,k
!!$do j=j1,j2
!!$  do i=i1,i2
!!$    do k=k1,k2
!!$      a(k,i,j)=b(k,i,j)+c(k,i,j)
!!$    enddo
!!$  enddo
!!$enddo
!!$return
!!$end
!!$!
!!$!     ******************************************************************
!!$!
!!$subroutine ae3m3(n1,n2,n3,i1,i2,j1,j2,k1,k2,a,b,c)
!!$implicit none
!!$integer :: n1,n2,n3,i1,i2,j1,j2,k1,k2
!!$real :: a(n1,n2,n3),b(n1,n2,n3),c(n1,n2,n3)
!!$integer :: i,j,k
!!$do j=j1,j2
!!$  do i=i1,i2
!!$    do k=k1,k2
!!$      a(k,i,j)=b(k,i,j)-c(k,i,j)
!!$    enddo
!!$  enddo
!!$enddo
!!$return
!!$end
!!$!
!!$!     ******************************************************************
!!$!
!!$subroutine ae3t3(n1,n2,n3,i1,i2,j1,j2,k1,k2,a,b,c)
!!$implicit none
!!$integer :: n1,n2,n3,i1,i2,j1,j2,k1,k2
!!$real :: a(n1,n2,n3),b(n1,n2,n3),c(n1,n2,n3)
!!$integer :: i,j,k
!!$do j=j1,j2
!!$  do i=i1,i2
!!$    do k=k1,k2
!!$      a(k,i,j)=b(k,i,j)*c(k,i,j)
!!$    enddo
!!$  enddo
!!$enddo
!!$return
!!$end
!!$!
!!$!     ******************************************************************
!!$!
!!$subroutine ae1p1p1(npts,a,b,c,f)
!!$implicit none
!!$integer :: npts
!!$real :: a(npts),b(npts),c(npts),f(npts)
!!$integer :: i
!!$do i=1,npts
!!$  a(i)=b(i)+c(i)+f(i)
!!$enddo
!!$return
!!$end
!!$!
!!$!     ******************************************************************
!!$!
!!$subroutine ae1t1p1(npts, a, b, c, f)
!!$  implicit none
!!$  include "i8.h"
!!$  ! Arguments:
!!$  integer(kind=i8), intent(IN) :: npts
!!$  real, intent(INOUT) :: a(npts)
!!$  real, intent(IN)    :: b(npts), c(npts), f(npts)
!!$  ! Local variables:
!!$  integer :: i
!!$  do i=1,npts
!!$     a(i) = b(i)*c(i) + f(i)
!!$  enddo
!!$end subroutine ae1t1p1
!!$!
!!$!     ******************************************************************
!!$!
!!$subroutine ae2(n2,n3,i1,i2,j1,j2,a,b)
!!$implicit none
!!$integer :: n2,n3,i1,i2,j1,j2
!!$real :: a(n2,n3),b(n2,n3)
!!$integer :: i,j
!!$do j=j1,j2
!!$  do i=i1,i2
!!$    a(i,j)=b(i,j)
!!$  enddo
!!$enddo
!!$return
!!$end
!!$!
!!$!     ******************************************************************
!!$!
!!$subroutine ae3t3p3(n1,n2,n3,i1,i2,j1,j2,k1,k2,a,b,c,f)
!!$implicit none
!!$integer :: n1,n2,n3,i1,i2,j1,j2,k1,k2
!!$real :: a(n1,n2,n3),b(n1,n2,n3),c(n1,n2,n3),f(n1,n2,n3)
!!$integer :: i,j,k
!!$do j=j1,j2
!!$  do i=i1,i2
!!$    do k=k1,k2
!!$      a(k,i,j)=b(k,i,j)*c(k,i,j)+f(k,i,j)
!!$    enddo
!!$  enddo
!!$enddo
!!$return
!!$end
!!$!
!!$!     ******************************************************************
!!$!
!!$subroutine ae3t0p3(n1,n2,n3,i1,i2,j1,j2,k1,k2,a,b,c,f)
!!$implicit none
!!$integer :: n1,n2,n3,i1,i2,j1,j2,k1,k2
!!$real :: a(n1,n2,n3),b(n1,n2,n3),c,f(n1,n2,n3)
!!$integer :: i,j,k
!!$do j=j1,j2
!!$  do i=i1,i2
!!$    do k=k1,k2
!!$      a(k,i,j)=b(k,i,j)*c+f(k,i,j)
!!$    enddo
!!$  enddo
!!$enddo
!!$return
!!$end
!!$!
!!$!     ******************************************************************
!!$!
!!$subroutine aen3t0p3(n1,n2,n3,i1,i2,j1,j2,k1,k2,a,b,c,f)
!!$implicit none
!!$integer :: n1,n2,n3,i1,i2,j1,j2,k1,k2
!!$real :: a(n1,n2,n3),b(n1,n2,n3),c,f(n1,n2,n3)
!!$integer :: i,j,k
!!$do j=j1,j2
!!$  do i=i1,i2
!!$    do k=k1,k2
!!$      a(k,i,j)=-b(k,i,j)*c+f(k,i,j)
!!$    enddo
!!$  enddo
!!$enddo
!!$return
!!$end
!!$!
!!$!     ******************************************************************
!!$!
!!$subroutine ae3m3d0(n1,n2,n3,i1,i2,j1,j2,k1,k2,a,b,c,f)
!!$implicit none
!!$integer :: n1,n2,n3,i1,i2,j1,j2,k1,k2
!!$real :: a(n1,n2,n3),b(n1,n2,n3),c(n1,n2,n3),f
!!$integer :: i,j,k
!!$
!!$do j=j1,j2
!!$  do i=i1,i2
!!$    do k=k1,k2
!!$      a(k,i,j)=(b(k,i,j)-c(k,i,j))/f
!!$    enddo
!!$  enddo
!!$enddo
!!$return
!!$end
!!$!
!!$!     ******************************************************************
!!$!
!!$subroutine a3e2(n1,n2,n3,i1,i2,j1,j2,k,a,b)
!!$implicit none
!!$integer :: n1,n2,n3,i1,i2,j1,j2,k
!!$real :: a(n1,n2,n3),b(n2,n3)
!!$integer :: i,j
!!$do j=j1,j2
!!$  do i=i1,i2
!!$    a(k,i,j)=b(i,j)
!!$  enddo
!!$enddo
!!$return
!!$end
!!$!
!!$!     ******************************************************************
!!$!
!!$subroutine a3e0(n1,n2,n3,i1,i2,j1,j2,k,a,b)
!!$implicit none
!!$integer :: n1,n2,n3,i1,i2,j1,j2,k
!!$real :: a(n1,n2,n3),b
!!$integer :: i,j
!!$do j=j1,j2
!!$  do i=i1,i2
!!$    a(k,i,j)=b
!!$  enddo
!!$enddo
!!$return
!!$end
!!$!
!!$!     ******************************************************************
!!$!
!!$subroutine alebl(n1,n2,n3,ka,kb,a,b)
!!$implicit none
!!$integer :: n1,n2,n3,ka,kb
!!$real :: a(n1,n2,n3),b(n1,n2,n3)
!!$integer :: i,j
!!$do j=1,n3
!!$  do i=1,n2
!!$    a(ka,i,j)=b(kb,i,j)
!!$  enddo
!!$enddo
!!$return
!!$end
!!$!
!!$!     ******************************************************************
!!$!
!!$subroutine ae0(npts,a,b)
!!$implicit none
!!$integer :: npts
!!$real :: a(npts),b
!!$integer :: i
!!$do i=1,npts
!!$  a(i)=b
!!$enddo
!!$return
!!$end
!!$!
!!$subroutine adivb(nnn,a,b,c)
!!$implicit none
!!$integer :: nnn,nn
!!$real :: a(nnn),b(nnn),c(nnn)
!!$do 1 nn=1,nnn
!!$c(nn)=a(nn)/b(nn)
!!$1 continue
!!$return
!!$end
!!$subroutine atimb(nnn,a,b,c)
!!$implicit none
!!$integer :: nnn,nn
!!$real :: a(nnn),b(nnn),c(nnn)
!!$do 1 nn=1,nnn
!!$c(nn)=a(nn)*b(nn)
!!$1 continue
!!$return
!!$end
!!$
!!$!     ******************************************************************

SUBROUTINE trid(var, cim1, ci, cip1, rhs, npts)

  !     SOLVES A DIAGONALLY-DOMINANT TRIDIAGONAL MATRIX EQUATION BY
  !     THE STANDARD QUICK METHOD
  !
  !     VAR   - VARIABLE BEING SOLVED FOR
  !     CIM1  - VECTOR OF COEFFICIENTS AT THE I-1 POINT
  !     CI    -   "     "       "       "  "  I     "
  !     CIP1  -   "     "       "       "  "  I+1   "
  !     RHS   -   "     "  THE RIGHT HAND SIDE OF THE EQUATION
  !     NPTS  - NUMBER OF EQUATIONS IN THE MATRIX
  !
  !     WARNING: THE ARRAYS CIM1,CI, AND RHS ARE REUSED IN THIS ROUTINE.

  IMPLICIT NONE
  ! Arguments:
  INTEGER, INTENT(in) :: npts
  REAL, INTENT(inout) :: var(npts)
  REAL, INTENT(in)    :: cim1(npts)
  REAL, INTENT(inout) :: ci(npts)
  REAL, INTENT(inout) :: cip1(npts)
  REAL, INTENT(inout) :: rhs(npts)
  ! Local variables:
  INTEGER :: k

  cip1(1) = cip1(1)/ci(1)
  DO k=2,npts
     ci(k)   = ci(k) - cim1(k)*cip1(k-1)
     cip1(k) = cip1(k)/ci(k)
  ENDDO

  rhs(1) = rhs(1)/ci(1)
  DO k=2,npts
     rhs(k) = (rhs(k) - cim1(k)*rhs(k-1))/ci(k)
  ENDDO

  var(npts) = rhs(npts)
  DO k=npts-1,1,-1
     var(k) = rhs(k)-cip1(k)*var(k+1)
  ENDDO

END SUBROUTINE trid

!!$!     ******************************************************************
!!$
!!$subroutine trid2(var,cim1,ci,cip1,rhs,npts,scr1,scr2,scr3)
!!$
!!$!     SOLVES A DIAGONALLY-DOMINANT TRIDIAGONAL MATRIX EQUATION BY
!!$!     THE STANDARD QUICK METHOD
!!$!
!!$!     VAR   - VARIABLE BEING SOLVED FOR
!!$!     CIM1  - VECTOR OF COEFFICIENTS AT THE I-1 POINT
!!$!     CI    -   "     "       "       "  "  I     "
!!$!     CIP1  -   "     "       "       "  "  I+1   "
!!$!     RHS   -   "     "  THE RIGHT HAND SIDE OF THE EQUATION
!!$!     NPTS  - NUMBER OF EQUATIONS IN THE MATRIX
!!$!     SCR1  - SCRATCH ARRAY AT LEAST NPTS LONG
!!$!     SCR2  - SCRATCH ARRAY "    "    "    "
!!$!     SCR3  - SCRATCH ARRAY "    "    "    "
!!$!
!!$implicit none
!!$integer :: npts
!!$real :: var(npts),cim1(npts),ci(npts),cip1(npts),rhs(npts)
!!$real :: scr1(npts),scr2(npts),scr3(npts)
!!$integer :: k
!!$
!!$scr1(1)=cip1(1)/ci(1)
!!$scr2(1)=ci(1)
!!$do 10 k=2,npts
!!$scr2(k)=ci(k)-cim1(k)*scr1(k-1)
!!$scr1(k)=cip1(k)/scr2(k)
!!$10 continue
!!$
!!$scr3(1)=rhs(1)/scr2(1)
!!$do 20 k=2,npts
!!$scr3(k)=(rhs(k)-cim1(k)*scr3(k-1))/scr2(k)
!!$20 continue
!!$
!!$var(npts)=scr3(npts)
!!$do 30 k=npts-1,1,-1
!!$var(k)=scr3(k)-scr1(k)*var(k+1)
!!$30 continue
!!$
!!$return
!!$end
!!$
!!$!     ******************************************************************

subroutine update(n, a, fa, dt)
  implicit none
  include "i8.h"
  integer :: n, nn
  real :: a(n), fa(n), dt
  integer(kind=i8) :: nlong
  real :: along(nlong), falong(nlong)

  entry update_single(n,a,fa,dt)
!!$do 10 nn=1,n
!!$  a(nn)=a(nn)+fa(nn)*dt
!!$10 continue
  a = a + fa * dt
  return

  entry update_long(nlong, along, falong, dt)
  along = along + falong * dt
  return

end subroutine update

!     ****************************************************************

subroutine accum(nxyz, arr1, arr2)
  implicit none
  include "i8.h"
  ! Arguments:
  integer(kind=i8), intent(IN) :: nxyz
  real, intent(INOUT)          :: arr1(nxyz)
  real, intent(IN)             :: arr2(nxyz)

  arr1 = arr1 + arr2

!!$  ! Local variables:
!!$  integer(kind=i8) :: n
!!$  do n=1,nxyz
!!$     arr1(n) = arr1(n) + arr2(n)
!!$  enddo
end subroutine accum

!     ******************************************************************

subroutine atob(n,a,b)
  implicit none
  include "i8.h"
  integer :: n,i
  real :: a(n),b(n)
  integer(kind=i8) :: nlong
  real :: along(nlong), blong(nlong)

  entry atob_single(n,a,b)
!!$do 100 i=1,n
!!$b(i)=a(i)
!!$100 continue
    b = a
  return

  entry atob_long(nlong, along, blong)
    blong = along
  return

end subroutine atob

!!$!     ******************************************************************
!!$
!!$subroutine acnst(n,a,cnst)
!!$implicit none
!!$integer :: n
!!$real :: a(n),cnst
!!$integer :: nn
!!$do 10 nn=1,n
!!$  a(nn)=cnst
!!$10 continue
!!$return
!!$end
!!$
!!$!     ******************************************************************
!!$
!!$real function valugp(n1,n2,n3,k,i,j,a)
!!$implicit none
!!$integer :: n1,n2,n3,k,i,j
!!$real :: a(n1,n2,n3)
!!$  valugp=a(k,i,j)
!!$return
!!$end
!!$
!!$!     ******************************************************************
!!$
!!$integer function ivalugp(n1,n2,n3,k,i,j,ia)
!!$implicit none
!!$integer :: n1,n2,n3,k,i,j
!!$real :: ia(n1,n2,n3)
!!$  ivalugp=ia(k,i,j)
!!$return
!!$end
!!$
!!$integer function ibindec(str)
!!$implicit none
!!$character(len=*) :: str
!!$integer :: inc,i,ic,l
!!$ibindec=0
!!$inc=1
!!$l=len(str)
!!$do ic=l,1,-1
!!$   if(str(ic:ic).eq.'1') ibindec=ibindec+inc
!!$   inc=inc+inc
!!$enddo
!!$return
!!$end
!!$
!!$!     ******************************************************************
!!$
!!$integer function ibias(y,l,n)
!!$implicit none
!!$integer :: l,n
!!$real :: y(l)
!!$integer :: jd,jpow,jpow1
!!$real :: ymin,ymax
!!$ymin=1.e10
!!$ymax=-1.e10
!!$do jd=1,l
!!$ymin=min(ymin,abs(y(jd)))
!!$ymax=max(ymax,abs(y(jd)))
!!$enddo
!!$jpow=int(log10(ymax+1.e-20)+.999999)
!!$jpow1=int(log10(ymax-ymin+1.e-20)+.999999)
!!$ibias=2-jpow1
!!$if(ibias+jpow.gt.4) ibias=4-jpow
!!$10    continue
!!$if(ibias+jpow.lt.4.and.ibias.lt.0)then
!!$ibias=ibias+1
!!$go to 10
!!$endif
!!$return
!!$end
!!$
!!$!     ******************************************************************
!!$
!!$real function heav(x)
!!$implicit none
!!$real :: x
!!$if(x.gt.0.)then
!!$heav=1.
!!$else
!!$heav=0.
!!$endif
!!$return
!!$end
!!$
!!$!     ******************************************************************
!!$
!!$integer function iprim(m)
!!$implicit none
!!$integer :: m,n
!!$n=m
!!$if(n.le.0)then
!!$  print 1,n
!!$1       format(' n=',i5,' in iprim')
!!$  stop
!!$endif
!!$10    continue
!!$if(mod(n,2).ne.0) go to 20
!!$  n=n/2
!!$  go to 10
!!$20    continue
!!$if(mod(n,3).ne.0) go to 30
!!$  n=n/3
!!$  go to 20
!!$30    continue
!!$if(mod(n,5).ne.0) go to 40
!!$  n=n/5
!!$  go to 30
!!$40    continue
!!$if(n.eq.1.and.mod(m,2).eq.0)then
!!$  iprim=1
!!$else
!!$  iprim =0
!!$endif
!!$return
!!$end
!!$
!!$!     ******************************************************************

subroutine sort3(a, b, c, n)
  implicit none
  integer, intent(IN) :: n
  real, intent(INOUT) :: a(n), b(n), c(n)
  integer :: np1, k, i
!!$  integer, external :: ismin
  real :: at, bt, ct
  np1 = n+1
  do k=1,n !10
!!$     i    = ismin(np1-k,a(k),1) + k-1
     i    = minloc(a(k:np1),1) + k-1
     at   = a(i)
     bt   = b(i)
     ct   = c(i)
     a(i) = a(k)
     b(i) = b(k)
     c(i) = c(k)
     a(k) = at
     b(k) = bt
     c(k) = ct
!!$10   continue
  enddo
end subroutine sort3


!!$!
!!$!-----------------------------------------------------------------------
!!$!        The following functions are the FORTRAN replacements for the
!!$!          CRAY intrinsic vector functions that perform various tasks.
!!$!          Since they are non-existent on
!!$!          other machines, they need to be replaced by calls to that
!!$!          machines's functions or simulated in standard FORTRAN.
!!$!
!!$!       Most or all now exist in f90 and could be replaced in the
!!$!        few places in the code where they are used.
!!$!-----------------------------------------------------------------------
!!$!
!!$!       Return sum of vector.
!!$!
!!$real function ssum(nn,vctr,inc)
!!$implicit none
!!$integer :: nn,inc
!!$real :: vctr(*)
!!$integer :: n,nnn
!!$real :: sum
!!$sum=0.
!!$nnn=nn*inc
!!$do 10 n=1,nnn,inc
!!$  sum=sum+vctr(n)
!!$10 continue
!!$ssum=sum
!!$return
!!$end
!!$! +------------------------------------------------------------------+
!!$!
!!$!       return index of maximum of vector
!!$!
!!$integer function ismax(nn,vctr,inc)
!!$implicit none
!!$integer :: nn,inc
!!$real :: vctr(*)
!!$integer :: ism,nnn
!!$real :: smax
!!$ism=0
!!$smax=-1e10
!!$do 10 nnn=1,nn,inc
!!$  if(vctr(nnn).gt.smax)then
!!$    ism=nnn
!!$    smax=vctr(nnn)
!!$  endif
!!$10 continue
!!$ismax=ism
!!$return
!!$end
!!$! +------------------------------------------------------------------+
!!$!
!!$!       return index of minimum of vector
!!$!
!!$integer function ismin(nn,vctr,inc)
!!$implicit none
!!$integer :: nn,inc
!!$real :: vctr(*)
!!$integer :: ism,nnn
!!$real :: smin
!!$ism=0
!!$smin=1e10
!!$do 10 nnn=1,nn,inc
!!$  if(vctr(nnn).lt.smin)then
!!$    ism=nnn
!!$    smin=vctr(nnn)
!!$  endif
!!$10 continue
!!$ismin=ism
!!$return
!!$end
! +------------------------------------------------------------------+
!
!       return vct1 if vct3 => 0., else vct2.
!
real function cvmgp(vct1, vct2, vct3)
  implicit none
  real, intent(IN) :: vct1, vct2, vct3
  if (vct3>=0) then
     cvmgp = vct1
  else
     cvmgp = vct2
  endif
end function cvmgp
! +------------------------------------------------------------------+
!
!       return vct1 if vct3 < 0., else vct2.
!
real function cvmgm(vct1, vct2, vct3)
  implicit none
  real, intent(IN) :: vct1, vct2, vct3
  if (vct3<0) then
     cvmgm = vct1
  else
     cvmgm = vct2
  endif
end function cvmgm
!!$! +------------------------------------------------------------------+
!!$!
!!$!       return vct1 if vct3 = 0., else vct2.
!!$!
!!$real function cvmgz(vct1, vct2, vct3)
!!$  implicit none
!!$  real, intent(IN) :: vct1, vct2, vct3
!!$  if (vct3==0) then
!!$     cvmgz = vct1
!!$  else
!!$     cvmgz = vct2
!!$  endif
!!$end function cvmgz
!!$! +------------------------------------------------------------------+
!!$!
!!$!       return vct1 if vct3 ne 0., else vct2.
!!$!
!!$real function cvmgn(vct1, vct2, vct3)
!!$  implicit none
!!$  real, intent(IN) :: vct1, vct2, vct3
!!$  if (vct3/=0) then
!!$     cvmgn = vct1
!!$  else
!!$     cvmgn = vct2
!!$  endif
!!$end function cvmgn

subroutine rcopy_i8(dest, source, npts)

  !**(JP)** para copiar um ponteiro escalar da var_tables
  !         como se fosse um array...

  implicit none
  include "i8.h"
  real, intent(out) :: dest(*)
  real, intent(in ) :: source(*)
  integer(kind=i8), intent(in) :: npts

  integer(kind=i8) :: i

  do i = 1, npts
     dest(i) = source(i)
  end do
end subroutine rcopy_i8

subroutine update_long_rk(nlong, dt,rk, a, a0, fa)
  implicit none
  include "i8.h"
  integer(kind=i8) :: nlong
  real ::dt,rk,dtxrk
  real :: a(nlong),a0(nlong), fa(nlong)

  dtxrk = dt *rk
  a = a0 + fa * dtxrk
  return

end subroutine update_long_rk

!     ****************************************************************
subroutine copy_long_rk(nlong, a, a0)
  implicit none
  include "i8.h"
  integer(kind=i8) :: nlong
  real :: a(nlong),a0(nlong)

  a0 = a
  return

end subroutine copy_long_rk
!     ****************************************************************
subroutine get_ab2_tend(nlong, c1,a1,c2,a2, a0)
  implicit none
  include "i8.h"
  integer(kind=i8) :: nlong
  real :: a1(nlong),a2(nlong),a0(nlong)
  real :: c1,c2
  a0 = c1*a1+c2*a2
  return

end subroutine get_ab2_tend
!     ****************************************************************
subroutine get_am3_tend(nlong, c1,a1,c2,a2,c3,a3, a0)
  implicit none
  include "i8.h"
  integer(kind=i8) :: nlong
  real :: a1(nlong),a2(nlong),a3(nlong),a0(nlong)
  real :: c1,c2,c3
  a0 = c1*a1+c2*a2+c3*a3
  return

end subroutine get_am3_tend

real function getDistance(LatA,LngA,LatB,LngB)
  !# Get distance between 2 points in lat,lon [m]
  !#
  !# @note
  !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
  !#
  !# **Brief**: using a equation to determine distance between 2 points on
  !# earth surface
  !#
  !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
  !#
  !# **Author**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
  !#
  !# **Date**: 28Mar2019
  !# @endnote
  !#
  !# @changes
  !#
  !# +
  !# @endchanges
  !# @bug
  !#
  !#@endbug
  !#
  !#@todo
  !#  &#9744; <br/>
  !# @endtodo
  !#
  !# @warning
  !# Now is under CC-GPL License, please see
  !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
  !# @endwarning
  !#
  !#---
  use dump, only: &
    dumpMessage

  include "constants.f90"
  character(len=*), parameter :: header="**(getDistance)**"

  !Input/output variables
  real,intent(in) :: latA
  !# Latitude of 1st POINT
  real,intent(in) :: lngA
  !# Longitude of 1st POINT
  real,intent(in) :: latB
   !# Latitude of 2nd POINT
  real,intent(in) :: lngB
  !# Longitude of 2nd POINT

  !Code
  getDistance = c_erad*acos(cos(c_pi*(90.-LatB)/180.)*cos((90.-LatA)*c_pi180) &
                + sin((90.-LatB)*c_pi180)*sin((90.-LatA)*c_pi180) &
                * cos(( LngA - LngB)*c_pi180))

end function getDistance
