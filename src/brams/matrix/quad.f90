!>----------------------------------------------------------------------------------------------------------------------
!! @brief    perform an n-point quadrature from 2n moments to yield  n abscissas and n weights.
!! @author   susanne bauer/doug wright
!!----------------------------------------------------------------------------------------------------------------------
subroutine Gauss_Matrix(n,x,r,w,zf)
   implicit none

   ! arguments.

   integer, intent(in   ) :: n               ! number of quadrature points
   real(8), intent(in   ) :: x(2*n)          ! moments
   real(8), intent(  out) :: r(n)            ! abscissas
   real(8), intent(  out) :: w(n)            ! weights
   real(8), intent(  out) :: zf              ! =0.0 successful quadrature, =1.0 failed quadrature

   ! local variables.

   integer :: ifailtql
   real(8) :: a(n),b(n),anu(2*n),amu0

   zf = 0.0d+00                              ! successful quadrature
   amu0 = x(1)                               ! normalizing moment
   anu(:) = x(:)/amu0                        ! normalize the moments
   call orthog(n,anu,a,b)
   call gaucof(n,a,b,amu0,r,w,ifailtql)

   
   if(     minval(r(:)) .lt. 0.0d+00 ) then  ! failed quadrature
     zf = 1.0d+00
   elseif( minval(w(:)) .lt. 0.0d+00 ) then  ! failed quadrature
     zf = 1.0d+00
   elseif( ifailtql .gt. 0 ) then            ! failed quadrature
     zf = 1.0d+00
   endif

end subroutine Gauss_matrix

!>----------------------------------------------------------------------------------------------------------------------
!!  @brief  see numerical recipes, w. press et al., 2nd edition.
!!----------------------------------------------------------------------------------------------------------------------
subroutine Orthog(n,anu,a,b)

   implicit none
   integer :: l,n,k,nmax
   parameter (nmax=30)
   real(8) :: a(n),anu(2*n),b(n),sig(2*nmax+1,2*nmax+1)

   do l=3,2*n
     sig(1,l)=0.d+00
   ENDDO
   do l=2,2*n+1
     sig(2,l)=anu(l-1)
   ENDDO
   a(1)=anu(2)/anu(1)
   b(1)=0.d+00
   do k=3,n+1
     do l=k,2*n-k+3
       sig(k,l)=sig(k-1,l+1)-a(k-2)*sig(k-1,l)-b(k-2)*sig(k-2,l)
       if(sig(k,k).le.0.d+00) sig(k,k) = 1.d-20
     END DO
     a(k-1)=sig(k,k+1)/sig(k,k)-sig(k-1,k)/sig(k-1,k-1)
     b(k-1)=sig(k,k)/sig(k-1,k-1)
   END DO

end subroutine Orthog


!>----------------------------------------------------------------------------------------------------------------------
!!  @brief  see numerical recipes, w. press et al., 2nd edition.
!!----------------------------------------------------------------------------------------------------------------------
subroutine GauCof(n,a,b,amu0,x,w,ifailtql)

   implicit none
   integer :: i,j,n,nmax, ifailtql
   parameter (nmax=30)
   real(8) :: a(n),b(n),w(n),x(n),z(nmax,nmax),amu0
   ifailtql = 0
   do i=1,n
     if(i.ne.1)b(i)=sqrt(b(i))
     do j=1,n
       if(i.eq.j)then
         z(i,j)=1.d+00
       else
         z(i,j)=0.d+00
       endif
     ENDDO
   ENDDO
   call tqli(a,b,n,nmax,z,ifailtql)
   if(ifailtql.gt.0) return
!-------------------------------------------------------------------------------------------------------------------
!  ordering of the abscissas is usually not needed.
!-------------------------------------------------------------------------------------------------------------------
!  call eigsrt(a,z,n,nmax)
!-------------------------------------------------------------------------------------------------------------------
   do i=1,n
     x(i)=a(i)
     w(i)=amu0*z(1,i)**2
     !--------------------------------------------------------------------------------------------------------------
     ! avoid zero weights.
     !--------------------------------------------------------------------------------------------------------------
     ! if(w(i).eq.0.d+00) w(i) = 1.d-30
     !--------------------------------------------------------------------------------------------------------------
   END DO
   return
   end subroutine GauCof


   subroutine Tqli(d,e,n,np,z,ifailtql)
   implicit none
!-------------------------------------------------------------------------------------------------------------------
!  see numerical recipes, w. press et al., 2nd edition.
!-------------------------------------------------------------------------------------------------------------------
   integer :: i,n,np,m,l,k,iter,ifailtql
   real(8) :: d(np),e(np),z(np,np),dd,g,r,s,c,p,f,b,pythag
   do i=2,n
     e(i-1)=e(i)
   ENDDO
   e(n)=0.d+00
   do l=1,n
     iter=0
1    do  m=l,n-1
       dd=abs(d(m))+abs(d(m+1))
       if (abs(e(m))+dd.eq.dd) goto 2
     ENDDO
     m=n
2    if(m.ne.l)then
       if(iter.eq.300) then
         ifailtql = 1
         return
       endif
       iter=iter+1
       g=(d(l+1)-d(l))/(2.d+00*e(l))
       r=pythag(g,1.0d0)
       g=d(m)-d(l)+e(l)/(g+sign(r,g))
       s=1.d+00
       c=1.d+00
       p=0.d+00
       do i=m-1,l,-1
         f=s*e(i)
         b=c*e(i)
         r=pythag(f,g)
         e(i+1)=r
         if(r.eq.0.d0)then
           d(i+1)=d(i+1)-p
           e(m)=0.d+00
           goto 1
         endif
         s=f/r
         c=g/r
         g=d(i+1)-p
         r=(d(i)-g)*s+2.d0*c*b
         p=s*r
         d(i+1)=g+p
         g=c*r-b
!  omit lines from here ...
         do k=1,n
           f=z(k,i+1)
           z(k,i+1)=s*z(k,i)+c*f
           z(k,i)=c*z(k,i)-s*f
         ENDDO
!  ... to here when finding only eigenvalues.
       ENDDO
       d(l)=d(l)-p
       e(l)=g
       e(m)=0.d+00
       goto 1
     endif
   ENDDO

end subroutine Tqli


!>----------------------------------------------------------------------------------------------------------------------
!! @brief   see numerical recipes, w. press et al., 2nd edition.
!!----------------------------------------------------------------------------------------------------------------------
 double precision function pythag(a,b)

   implicit none
   real(8) :: absa,absb, a,b,phytag
   absa=abs(a)
   absb=abs(b)
   if(absa.gt.absb)then
     pythag=absa*sqrt(1.d+00+(absb/absa)**2)
   else
     if(absb.eq.0.d+00)then
       pythag=0.d+00
     else
       pythag=absb*sqrt(1.d+00+(absa/absb)**2)
     endif
   endif

end function pythag


!>----------------------------------------------------------------------------------------------------------------------
!!  @brief   see numerical recipes, w. press et al., 2nd edition.
!!----------------------------------------------------------------------------------------------------------------------
subroutine Eigsrt(d,v,n,np)

   implicit none
   integer :: n,np,k,j,i
   real(8) :: d(np),v(np,np),p
   do i=1,n-1
     k=i
     p=d(i)
     do j=i+1,n
       if(d(j).ge.p)then
         k=j
         p=d(j)
       endif
     ENDDO
     if(k.ne.i)then
       d(k)=d(i)
       d(i)=p
       do j=1,n
         p=v(j,i)
         v(j,i)=v(j,k)
         v(j,k)=p
        ENDDO
     endif
   ENDDO

end subroutine Eigsrt


!>----------------------------------------------------------------------------------------------------------------------
!! @brief    dlw: 091206: computes the moments from the abscissas and weights.
!!                  this routine is independent of the units used.
!!----------------------------------------------------------------------------------------------------------------------
subroutine GaussInv(n,x,w,u)
   implicit none

   ! arguments.

   integer :: n          ! number of quadrature points
   real(8) :: u(2*n)     ! moments
   real(8) :: x(n)       ! abscissas
   real(8) :: w(n)       ! weights
 
   ! local variables.

   integer :: i 

   do i=1, 2*n
     u(i) = sum( w(:) * ( x(:)**(i-1) ) )
   enddo

end subroutine GaussInv


!>----------------------------------------------------------------------------------------------------------------------
!!  @brief  dlw, 091306: check quadrature routines. 
!!----------------------------------------------------------------------------------------------------------------------
subroutine Test_Quad
   implicit none
   integer, parameter :: n = 6               ! number of quadrature points
   real(8)            :: x(2*n)              ! moments
   real(8)            :: r(n)                ! abscissas
   real(8)            :: w(n)                ! weights
   real(8)            :: zf                  ! =0.0 successful quadrature, =1.0 failed quadrature
   integer            :: i,k
   real(8)            :: n0,dg,sigmag,sg 

   zf = 0.0d+00
   n0 = 1.0d+03
   dg = 0.1d+00
   sigmag = 1.6d+00
   sg = exp( 0.5d+00 * ( log(sigmag) )**2 )
   do i=1, 2*n
     k = i-1
     x(i) = n0 * dg**k * sg**(k*k)
   enddo
   write(*,90) x(:)
   call Gauss_matrix(n,x,r,w,zf)
   write(*,90) r(:),w(:)
   call GaussInv(n,r,w,x)
   write(*,90) x(:)

90 format(50d18.10)

end subroutine Test_Quad


