!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine error_mess(msg)

  ! error_mess: dumps a message on stdout

  implicit none
  character(len=*), intent(in) :: msg
  write(*,"(a)") msg
end subroutine error_mess
