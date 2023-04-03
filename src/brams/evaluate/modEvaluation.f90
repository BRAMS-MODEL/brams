module evaluation
!# Module for statistics of model
!#
!# @note
!# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
!#
!# **Brief**: This module evaluate 5 variables by compare they with CI/CC
!#
!# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
!#
!# **Author**: Luiz Flavio Rodrigues **&#9993;**<mailto:luiz.rodrigues@inpe.br>
!#
!# **Date**: 2019-11-22
!# @endnote
!#
!# @changes
!#
!# +
!# @endchanges
!# @bug
!# No active bugs reported now
!# @endbug
!#
!# @todo
!#  &#9744; <br/>
!# @endtodo
!#
!# @warning
!# Now is under CC-GPL License, please see
!# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
!# @endwarning
!#
!#--- ----------------------------------------------------------------------------------------
!  
use ParLib  , only: & !SUbroutines for parallel comunications
    parf_send_noblock_real, &
    parf_get_noblock_real , &
    parf_wait_any_nostatus, &
    parf_wait_all_nostatus

use node_mod, only:  &
       mynum, &
       mchnum, &
       master_num, &
       nmachs, &
       mxp,            &
       myp, &
       nzp

use mem_basic: &
       nnxp, &
       nnyp, &
       nnzp

implicit none

private

public allocStatistic,comunicateStatistic,statistic

type st_vars
  real,pointer,dimension(:) :: rmse !nLeves
  !# Root mean square error
  real,pointer,dimension(:) :: bias
  !# bias error
end type st_vars
type(st_vars),dimension(:,:) :: statistic !nTimes,nvars
!# variable with statistics
integer,parameter :: nVars=5
character(len=2),dimension(nVars),parameter=(/"UP","VP","TH","RP","PP"/)

contains

subroutine allocStatistic(nTimes)
  !# Allocate statistic variable
  !#
  !# @note
  !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
  !#
  !# **Brief**: Allocate statistic Variables RMSE e BIAS. RMSE and Bias is calculated
  !# over all points x and y and one for each layer.
  !#
  !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
  !#
  !# **Author**: Luiz Flavio Rodrigues **&#9993;**<mailto:luiz.rodrigues@inpe.br>
  !#
  !# **Date**: 2019-11-29
  !# @endnote
  !#
  !# @changes
  !#
  !# +
  !# @endchanges
  !# @bug
  !# No active bugs reported now
  !# @endbug
  !#
  !# @todo
  !#  &#9744; <br/>
  !# @endtodo
  !#
  !# @warning
  !# Now is under CC-GPL License, please see
  !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
  !# @endwarning
  !#
  !#--- ----------------------------------------------------------------------------------------
  !  
  implicit none

	integer, intent(in) :: nTimes
	!# total of timesteps
  integer :: nt,nv,ierr

	allocate (statistic(nTimes,nVars), STAT=ierr)
  if (ierr/=0) CALL fatal_error("ERROR allocating up (statistic(nTimes,nVars))")
	
  do nt=1,nTimes
    do nv=1,nVars
	    allocate(statistic(nt,nv)%rmse(nnzp), STAT=ierr)
      if (ierr/=0) CALL fatal_error("ERROR allocating up (rmse)")
      statistic(nt,nv)%rmse=0.0

      allocate(statistic(nt,nv)%bias(nnzp), STAT=ierr)
      if (ierr/=0) CALL fatal_error("ERROR allocating up (bias)")
      statistic(nt,nv)%bias=0.0

    enddo
  enddo

end subroutine allocStatistic

subroutine comunicateStatistic(soma,somaQ)
  !# Send SUM local information for master 
  !#
  !# @note
  !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
  !#
  !# **Brief**: This subroutine send sum total from all processor to master and
  !# calculates and store RMSE and BIAS in master processor
  !#
  !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
  !#
  !# **Author**: Luiz Flavio Rodrigues **&#9993;**<mailto:luiz.rodrigues@inpe.br>
  !#
  !# **Date**: 2019-11-22
  !# @endnote
  !#
  !# @changes
  !#
  !# +
  !# @endchanges
  !# @bug
  !# No active bugs reported now
  !# @endbug
  !#
  !# @todo
  !#  &#9744; <br/>
  !# @endtodo
  !#
  !# @warning
  !# Now is under CC-GPL License, please see
  !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
  !# @endwarning
  !#
  !#--- ----------------------------------------------------------------------------------------
  !   
  real, intent(in) :: soma(nnzp,nVars)
  !# Sum of error for each var
  real, intent(in) :: somaQ(nnzp,nVars)
  !# Sum of quadratic error for each var
  type dm 
    real,pointer,allocatable :: dados(:)
  end type dm
  type(dm),allocatable :: dataMess(:)

  !for example $a^b = \int{\phi^2(x) dx}$

  integer :: sizeMess,iRecv,iProc,iCount,tag,recNum
  real,allocatable :: dataMess
  integer, allocatable :: reqRecv,reqSend

  sizeMess=nVars*nnzp*2

  allocate(dataMess(nnzp))
  do iProc=1,nMachs
    allocate(dataMess(iProc)%dados(sizeMess))
  enddo
   
  iCount=0

  if(mynum/=master_num) then
    do nv=1,nVars
      do z=1,nnzp
        iCount=iCount+1
        sendMess(iCount)=soma(z,nv)
        iCount=iCount+1
        sendMess(iCount)=somaQ(z,nv)
      enddo
    enddo
  endif

  if(mynum==master_num) then
    do iRecv = 1, nMachs      
      tag=1000+iRecv   
      call parf_get_noblock_real(dataMess(iRecv)%dados,sizeMess,iRecv,tag,reqRecv(iRecv))
    end do
  endif
  !
  if(mynum==master_num) then
    do iRecV = 1, nMachs
      call parf_wait_any_nostatus(iRecv,reqRecv,recNum)
      iCount=0
      do nv=1,nVars
        do z=1,nnzp
          iCount=iCount+1
          soma(z,nv)=dataMess(recNum)+soma(z,nv)
          iCount=iCount+1
          somaQ(z,nv)=dataMess(recNum)+somaQ(z,nv)
        enddo
      enddo
    enddo
  else
    tag=1000+mynum
    call parf_send_noblock_real(sendMess,sizeMess,0,tag,reqSend)
  endif

  call computeStatistic(soma,somaQ)

end subroutine comunicateStatistic

subroutine computeStatistic(soma,somaQ)
  real, intent(in) :: soma(nnzp,nVars)
  !# Sum of error for each var
  real, intent(in) :: somaQ(nnzp,nVars)
  !# Sum of quadratic error for each var

  integer :: totColumns,k,v

  totColumns=nnxp*nnyp

  do k=1,nnzp
    do v=1,5
      statistic(t,v)%rmse(k)=sqrt(somaQ(k,v)/totColumns)
      statistic(t,v)%bias(k)=soma(k,v)/totColumns
    enddo
  enddo

end subroutine computeStatistic

subroutine writeStatistics()


end subroutine writeStatistics

end module evaluation
