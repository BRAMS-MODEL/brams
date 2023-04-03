!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

! +-------------------------------------------------------------------+
! _                     SYSTEM DEPENDENT ROUTINES                     _
! _                                                                   _
! _  This module contains short utility routines that are not         _
! _  of the FORTRAN 77 standard and may differ from system to system. _
! _  These include bit manipulation, I/O, JCL calls, and vector       _
! _  functions.                                                       _
! +-------------------------------------------------------------------+

!----------------------------------------------------------------------------

integer function iran_recsize()
  implicit none
      
  !     Returns a factor for the length of a random access record 
  !     length unit.
  !     For example, IBM specifies bytes while SGI uses words, so
  !     specify the open statement in WORDS, then multiply by this 
  !     returned value.

#if defined(IBM) || defined(SUN) || defined(HP)
  iran_recsize=4
#endif

#if defined(PC_LINUX1) || defined(PC_NT1) || defined(SGI) || defined(NEC_SX)
  iran_recsize=4
#endif

#if defined(ALPHA) || defined(CRAY) 
  iran_recsize=1
#endif

  return
end function iran_recsize

!----------------------------------------------------------------------------

subroutine timing(icall,t1)
  !     Routine returns CPU time.  Called with ICALL=1 at beginning
  !     of timestep, ICALL=2 at end of timestep.
  implicit none
  !Arguments:
  integer, intent(in) :: icall
  real, intent(out)   :: t1
  ! Local variables:

  call cpu_time(t1)

end subroutine timing

!***************************************************************************

subroutine dcw_swap16 (a,n)
  implicit none

  ! reverse order of two bytes in integer*2
  integer :: n
  integer(kind=2) ::  a(n)

#if defined(SGI)

  integer(kind=2) :: itemp
  character(len=1) :: jtemp(2),ktemp
  equivalence  (itemp,jtemp(1))

  integer :: i

  do i=1,n
     itemp    = a(i)
     ktemp    = jtemp(1)
     jtemp(1) = jtemp(2)
     jtemp(2) = ktemp
     a(i)     = itemp
  enddo

#endif

  return
end subroutine dcw_swap16

!***************************************************************************

subroutine dcw_swap32 (a,n)
  implicit none

  ! reverse order of bytes in integer*4 word, or real*4
  integer :: n
  integer(kind=4) :: a(n)

#if defined(SGI)

  integer(kind=4) :: itemp

  character(len=1) :: jtemp(4), ktemp
  equivalence (jtemp(1),itemp)

  integer :: i

  do i=1,n
     itemp    = a(i)
     ktemp    = jtemp(4)
     jtemp(4) = jtemp(1)
     jtemp(1) = ktemp
     ktemp    = jtemp(3)
     jtemp(3) = jtemp(2)
     jtemp(2) = ktemp
     a(i)     = itemp
  enddo

#endif

  return
end subroutine dcw_swap32

!***************************************************************************

subroutine dcw_swap64 (a,n)
  implicit none

  ! reverse order of eight bytes in real*8 word
  integer :: n
  real(kind=8) ::  a(n)

#if defined(SGI)

  real(kind=8) :: itemp

  character(len=1) :: jtemp(8), ktemp
  equivalence (jtemp(1),itemp)

  integer :: i

  do i = 1,n
     itemp    = a(i)
     ktemp    = jtemp(8)
     jtemp(8) = jtemp(1)
     jtemp(1) = ktemp
     ktemp    = jtemp(7)
     jtemp(7) = jtemp(2)
     jtemp(2) = ktemp
     ktemp    = jtemp(6)
     jtemp(6) = jtemp(3)
     jtemp(3) = ktemp
     ktemp    = jtemp(5)
     jtemp(5) = jtemp(4)
     jtemp(4) = ktemp
     a(i)     = itemp
  enddo

#endif

  return
end subroutine dcw_swap64
