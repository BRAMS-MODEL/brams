module ModEvaluation
!# Module for statistics of model
!#
!# @note
!# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
!#
!# **Brief**: This module evaluate 5 variables by compare they with CI/CC. The results are some
!# files with RMSE and BIAS for UP (U Wind), VP (VWind), TH (Theta), PP (Pressure) and RP.
!# The quadratic sum of diferences and sum of diferences each level are computed in 
!# in in nud_analisys.f90. For each timestep the values are:
!# 
!# \[${RMSE}_{k_var_{var}}={\sqrt{\sum_{k=1}^{nzp} {(V_{p_{k,i,j}} - V_{o_{k,i,j}})}^{2}}\over{nzp}}$\]
!#
!# \[${BIAS}_{k_{var}}={{\sum_{k=1}^{nzp} {(V_{p_{k,i,j}} - V_{o_{k,i,j}})}}\over{nzp}}$\]
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
       nmachs

use mem_grid, only: &
       nnxp, &
       nnyp, &
       nnzp, &
       time, &
       expnme, &
       ztn

use ModNamelistFile, only: namelistFile

implicit none

private

public allocStatistic,comunicateStatistic,statistic,timeCount,nTimes, &
       StoreNamelistFileAtEvaluate,RMSE_average,evaluate

integer :: evaluate
!# if is to run (1) or not (.not. 1) the evaluation. 
character(len=256) :: evaluatePrefix
!# File prefix for output file
integer, parameter :: fileNum=25
character(len=*), parameter :: outputFormat="csv"
!# Define the output format that may be gnp (gnuplot) or csv 

type st_vars
  real,pointer,dimension(:) :: rmse !nLeves
  !# Root mean square error
  real,pointer,dimension(:) :: bias
  !# bias error
  real :: average_rmse
  !# RMSE average in all colunmns
end type st_vars
type(st_vars),allocatable,dimension(:,:) :: statistic !nTimes,nvars
!# variable with statistics
real, allocatable :: RMSE_average(:)
!# RMSE average in all levels and vars

integer,parameter :: nVars=5
!# total of vars
character(len=2),dimension(nVars),parameter :: VarName=(/"UP","VP","TH","RP","PP"/)
!# name of each var
integer :: timeCount
!# sequential (from 1 to end) count of timesteps
integer :: nTimes
!# total number of timesteps
real,allocatable :: timeOfError(:)
!# Time in seconds for each timeCount

contains

subroutine allocStatistic(nTimes)
  !# Allocate statistic variable
  !#
  !# @note
  !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
  !#
  !# **Brief**: Allocate statistic Variables RMSE, BIAS, etc. RMSE and Bias is calculated
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

  if(evaluate/=1) return

  !write(*,fmt='(A,2(I5.5,1X))') 'Allocating statistic with ',nTimes,nVars
	allocate (statistic(nTimes,nVars), STAT=ierr)
  if (ierr/=0) CALL fatal_error("ERROR allocating up (statistic(nTimes,nVars))")
	
  !write(*,fmt='(A,I3.3)') 'Allocating RMSE and BIAS with ',nnzp(1)
  do nt=1,nTimes
    do nv=1,nVars
	    allocate(statistic(nt,nv)%rmse(nnzp(1)), STAT=ierr)
      if (ierr/=0) CALL fatal_error("ERROR allocating up (rmse)")
      statistic(nt,nv)%rmse=0.0

      allocate(statistic(nt,nv)%bias(nnzp(1)), STAT=ierr)
      if (ierr/=0) CALL fatal_error("ERROR allocating up (bias)")
      statistic(nt,nv)%bias=0.0

    enddo
  enddo
  allocate(timeOfError(nTimes))
  allocate(RMSE_average(nTimes))

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
  real, intent(inout) :: soma(nnzp(1),nVars)
  !# Sum of error for each var
  real, intent(inout) :: somaQ(nnzp(1),nVars)
  !# Sum of quadratic error for each var
  type dm 
    real,pointer,dimension(:) :: dados
    !# Data to be communicated
  end type dm
  type(dm),allocatable :: dataMess(:)
  !# Comunicated data from each processor


  integer :: sizeMess,iRecv,iProc,iCount,tag,recNum,nv,z
  integer, allocatable :: reqRecv(:)
  integer :: reqSend
  real, allocatable :: sendMess(:)

  if(evaluate/=1) return !No evaluate are set in RAMSIN

  sizeMess=nVars*nnzp(1)*2 !Define size of message for all vars, levels and RMSE & BIAS

  allocate(sendMess(nnzp(1)*nVars*2))
  allocate(reqRecv(nMachs))

  allocate(dataMess(nnzp(1)))
  do iProc=1,nMachs
    allocate(dataMess(iProc)%dados(sizeMess))
  enddo

  iCount=0
  !Fill and pack the message to be send to master with Soma and SomaQ values
  if(mchnum/=master_num) then
    do nv=1,nVars
      do z=1,nnzp(1)
        iCount=iCount+1
        sendMess(iCount)=soma(z,nv)
      enddo
    enddo
    do nv=1,nVars
      do z=1,nnzp(1)
        iCount=iCount+1
        sendMess(iCount)=somaQ(z,nv)
      enddo
    enddo
  endif

  !If master ask to each processor for local soma and somaQ
  if(mchnum==master_num) then
    do iRecv = 1, nMachs-1      
      tag=1000+iRecv   
      !write(40+mchnum,fmt='(2(A,I5.5,1X))') 'get from ',irecV, 'with tag=',tag; call flush(40+mchnum)
      call parf_get_noblock_real(dataMess(iRecv)%dados,sizeMess,iRecv,tag,reqRecv(iRecv))
    end do
  endif
  !
  !If master wait for messages sent by processors one-by-one in order of comming
  if(mchnum==master_num) then
    do iRecV = 1, nMachs-1
      call parf_wait_any_nostatus(iRecv,reqRecv,recNum)
      !For each recNum message unpack it and add in local (master) soma and somaQ
      iCount=0
      !
      do nv=1,nVars
        do z=1,nnzp(1)
          iCount=iCount+1
          soma(z,nv)=dataMess(recNum)%dados(iCount)+soma(z,nv)
        enddo
      enddo
      do nv=1,nVars
        do z=1,nnzp(1)
          iCount=iCount+1
          somaQ(z,nv)=dataMess(recNum)%dados(iCount)+somaQ(z,nv)
        enddo
      enddo
    enddo

    ! All the soma and somaQ from slaves are summed. Compute now RMSE and BIAS
    call computeStatistic(soma,somaQ)

  else
    !If slave send message to master
    tag=1000+mchnum
    !write(40+mchnum,fmt='(3(A,I5.5,1X))') 'I am ', mchnum,'send to ',master_num, 'with tag=',tag; call flush(6)
    call parf_send_noblock_real(sendMess,sizeMess,master_num,tag,reqSend)
  endif

end subroutine comunicateStatistic

subroutine computeStatistic(soma,somaQ)
  !# Computes the RMSE and Bias and make the average
  !#
  !# @note
  !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
  !#
  !# **Brief**: This subroutine compute RMSE and BIAS using the equation
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
  real, intent(in) :: soma(nnzp(1),nVars)
  !# Sum of error for each var
  real, intent(in) :: somaQ(nnzp(1),nVars)
  !# Sum of quadratic error for each var

  integer :: totColumns,k,v
  real :: totalRmse

  totColumns=nnxp(1)*nnyp(1)
  !Increment the timeCount
  timeCount=timeCount+1

  do v=1,nVars
    totalRmse=0.0
    !Compute RMSE and BIAS for each k
    do k=1,nnzp(1)
      statistic(timeCount,v)%rmse(k)=sqrt(somaQ(k,v)/totColumns)
      statistic(timeCount,v)%bias(k)=soma(k,v)/totColumns
      !Sum the RMSE in all levels for each var
      totalRmse=totalRmse+statistic(timeCount,v)%rmse(k)
    enddo
    !Compute the average for each var
    statistic(timeCount,v)%average_rmse=totalRmse/nnzp(1)
  enddo
  !Put the time of timestep in timeOfError for this timeCount
  timeOfError(timeCount)=time
  totalRmse=0.0
  !Sum average for each var
  do v=1,nVars
    totalRmse=totalRmse+statistic(timeCount,v)%average_rmse
  enddo
  !Make the all average to put in message of each timestep
  RMSE_average(timeCount)=totalRmse/5.0
  
  !Verify if is time to write data
  if(timeCount==nTimes) call writeStatistics()

end subroutine computeStatistic

subroutine writeStatistics()
  !# write statiscs files
  !#
  !# @note
  !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
  !#
  !# **Brief**: This subroutine write RMSE and BIAS in specific files
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
  integer :: nt,nv,nz,z,nLev
  character(len=256) :: fileName,file2Open
  character(len=3) :: cLev
  character(len=3) :: cPos

  write(cLev,fmt='(I3.3)') nnzp(1)

  if(outputFormat=="csv") then
    do nv=1,nVars
      !Open the CSV file if CSV is set in module header
      fileName=trim(evaluatePrefix)//trim(expnme)//"-"//varName(nv)//".RMSE.csv"
      !Open the output file
      call rams_f_open(fileNum,fileName,'FORMATTED','REPLACE','WRITE',1)
      !Write a header in first line
      write(fileNum,fmt='("time",'//cLev//'(";",F8.2))') (ztn(z,1), z=1, nnzp(1))
      !Write RMSE in files with time and each level RMSE value
      do nt=1,nTimes
        write(filenum,fmt='(F12.1,1X,'//cLev//'(";",E9.2))') timeOfError(nt),statistic(nt,nv)%rmse
      enddo
      close(fileNum)
      !Make the same thing with BIAS
      fileName=trim(evaluatePrefix)//trim(expnme)//"-"//varName(nv)//".BIAS.csv"
      !Open the output file
      call rams_f_open(fileNum,fileName,'FORMATTED','REPLACE','WRITE',1)
      write(fileNum,fmt='("time",'//cLev//'(";",F8.2))') (ztn(z,1), z=1, nnzp(1))
      do nt=1,nTimes
        write(filenum,fmt='(F12.1,1X,'//cLev//'(";",E9.2))') timeOfError(nt),statistic(nt,nv)%bias
      enddo
      close(fileNum)
    enddo

    do nv=1,nVars
      fileName=trim(evaluatePrefix)//trim(expnme)//"-"//varName(nv)//".RMSE.py"
      call rams_f_open(fileNum,fileName,'FORMATTED','REPLACE','WRITE',1)

      file2Open=trim(evaluatePrefix)//trim(expnme)//"-"//varName(nv)//".RMSE.csv"
      write(fileNum,fmt='(A)') '#!/usr/bin/python'
      write(fileNum,fmt='(A)') 'import csv'
      write(fileNum,fmt='(A)') 'import matplotlib.pyplot as plt'
      write(fileNum,fmt='(A)') 'import numpy as np'
      write(fileNum,fmt='(A)') 'from scipy.ndimage.filters import gaussian_filter'
      write(fileNum,fmt='(A,A,A)') 'arquivo = open("',trim(file2Open),'")'
      write(fileNum,fmt='(A)') 'linhas = csv.reader(arquivo,delimiter=";")'
      write(fileNum,fmt='(A)') 'cabec=True'
      write(fileNum,fmt='(A)') 'valTime =[]'
      write(fileNum,fmt='(A)') 'valCampo=[]'
      write(fileNum,fmt='(A)') 'lev=[]'
      write(fileNum,fmt='(A)') 'for linha in linhas:'
      write(fileNum,fmt='(A)') '  if cabec:'
      write(fileNum,fmt='(A)') '    xlabel=linha[0]'
      write(fileNum,fmt='(A)') '    ylabel=linha[1:]'
      write(fileNum,fmt='(A)') '    for i in range(len(ylabel)):'
      write(fileNum,fmt='(A)') '      lev.append(float(ylabel[i]))'
      write(fileNum,fmt='(A)') '      valCampo.append([])'
      write(fileNum,fmt='(A)') '    cabec=False'
      write(fileNum,fmt='(A)') '    totalFields=len(ylabel)'
      write(fileNum,fmt='(A)') '    continue'
      write(fileNum,fmt='(A)') '  valTime.append(float(linha[0]))'
      write(fileNum,fmt='(A)') '  for i in range(totalFields):'
      write(fileNum,fmt='(A)') '    valCampo[i].append(float(linha[i+1]))'
      write(fileNum,fmt='(A)') ''
      write(fileNum,fmt='(A)') 'cmap = plt.get_cmap("rainbow")'
      write(fileNum,fmt='(A)') 'colors = [cmap(i) for i in np.linspace(0, 1, totalFields)]'
      write(fileNum,fmt='(A)') '#for i in range(totalFields):'
      write(fileNum,fmt='(A)') '#  plt.plot(valTime,valCampo[i][:],color=colors[i])'
      write(fileNum,fmt='(A)') ''
      write(fileNum,fmt='(A)') 'valCampo=gaussian_filter(valCampo, 1.5)' 
      write(fileNum,fmt='(A)') 'plt.contourf(valTime,lev,valCampo)'
      write(fileNum,fmt='(A)') 'plt.colorbar()'
      write(fileNum,fmt='(A)') ''
      write(fileNum,fmt='(A,A,A,I3.3,A)') 'plt.title("BRAMS - Evaluate Module - ',varName(nv),' - RMSE") '
      write(fileNum,fmt='(A)') 'plt.ylabel("Levels [m]")'
      write(fileNum,fmt='(A)') 'plt.xlabel("Time [Seconds]")'
      write(fileNum,fmt='(A)') 'plt.minorticks_on()'
      write(fileNum,fmt='(A)') 'plt.xticks(rotation=45)'
      if(timeOfError(ntimes)<=3600.0) then
        write(fileNum,fmt='(A,F18.1,A)') 'plt.xticks(np.arange(0, ',timeOfError(ntimes),', 600.0))'
      else
        write(fileNum,fmt='(A,F18.1,A)') 'plt.xticks(np.arange(0, ',timeOfError(ntimes),', 3600.0))'
      endif
      write(fileNum,fmt='(A)') 'plt.grid(b=True, which="major", color="#666666", linestyle="-")'
      write(fileNum,fmt='(A)') 'plt.grid(b=True, which="minor", color="#666666", linestyle=":")'
      write(fileNum,fmt='(A)') '#plt.legend(ylabel,loc=7,ncol=2,bbox_to_anchor=(1.1, 0.5))'
      write(fileNum,fmt='(A)') 'plt.show()'
      close(fileNum)
    enddo

    do nv=1,nVars
      fileName=trim(evaluatePrefix)//trim(expnme)//"-"//varName(nv)//".BIAS.py"
      call rams_f_open(fileNum,fileName,'FORMATTED','REPLACE','WRITE',1)

      file2Open=trim(evaluatePrefix)//trim(expnme)//"-"//varName(nv)//".BIAS.csv"
      write(fileNum,fmt='(A)') '#!/usr/bin/python'
      write(fileNum,fmt='(A)') 'import csv'
      write(fileNum,fmt='(A)') 'import matplotlib.pyplot as plt'
      write(fileNum,fmt='(A)') 'import numpy as np'
      write(fileNum,fmt='(A)') 'from scipy.ndimage.filters import gaussian_filter'
      write(fileNum,fmt='(A,A,A)') 'arquivo = open("',trim(file2Open),'")'
      write(fileNum,fmt='(A)') 'linhas = csv.reader(arquivo,delimiter=";")'
      write(fileNum,fmt='(A)') 'cabec=True'
      write(fileNum,fmt='(A)') 'valTime =[]'
      write(fileNum,fmt='(A)') 'valCampo=[]'
      write(fileNum,fmt='(A)') 'lev=[]'
      write(fileNum,fmt='(A)') 'for linha in linhas:'
      write(fileNum,fmt='(A)') '  if cabec:'
      write(fileNum,fmt='(A)') '    xlabel=linha[0]'
      write(fileNum,fmt='(A)') '    ylabel=linha[1:]'
      write(fileNum,fmt='(A)') '    for i in range(len(ylabel)):'
      write(fileNum,fmt='(A)') '      lev.append(float(ylabel[i]))'
      write(fileNum,fmt='(A)') '      valCampo.append([])'
      write(fileNum,fmt='(A)') '    cabec=False'
      write(fileNum,fmt='(A)') '    totalFields=len(ylabel)'
      write(fileNum,fmt='(A)') '    continue'
      write(fileNum,fmt='(A)') '  valTime.append(float(linha[0]))'
      write(fileNum,fmt='(A)') '  for i in range(totalFields):'
      write(fileNum,fmt='(A)') '    valCampo[i].append(float(linha[i+1]))'
      write(fileNum,fmt='(A)') ''
      write(fileNum,fmt='(A)') 'cmap = plt.get_cmap("rainbow")'
      write(fileNum,fmt='(A)') 'colors = [cmap(i) for i in np.linspace(0, 1, totalFields)]'
      write(fileNum,fmt='(A)') '#for i in range(totalFields):'
      write(fileNum,fmt='(A)') '#  plt.plot(valTime,valCampo[i][:],color=colors[i])'
      write(fileNum,fmt='(A)') ''
      write(fileNum,fmt='(A)') 'valCampo=gaussian_filter(valCampo, 1.5)' 
      write(fileNum,fmt='(A)') 'plt.contourf(valTime,lev,valCampo)'
      write(fileNum,fmt='(A)') 'plt.colorbar()'
      write(fileNum,fmt='(A)') ''
      write(fileNum,fmt='(A,A,A,I3.3,A)') 'plt.title("BRAMS - Evaluate Module - ',varName(nv),' - BIAS") '
      write(fileNum,fmt='(A)') 'plt.ylabel("Levels [m]")'
      write(fileNum,fmt='(A)') 'plt.xlabel("Time [Seconds]")'
      write(fileNum,fmt='(A)') 'plt.minorticks_on()'
      write(fileNum,fmt='(A)') 'plt.xticks(rotation=45)'
      if(timeOfError(ntimes)<=3600.0) then
        write(fileNum,fmt='(A,F18.1,A)') 'plt.xticks(np.arange(0, ',timeOfError(ntimes),', 600.0))'
      else
        write(fileNum,fmt='(A,F18.1,A)') 'plt.xticks(np.arange(0, ',timeOfError(ntimes),', 3600.0))'
      endif
      write(fileNum,fmt='(A)') 'plt.grid(b=True, which="major", color="#666666", linestyle="-")'
      write(fileNum,fmt='(A)') 'plt.grid(b=True, which="minor", color="#666666", linestyle=":")'
      write(fileNum,fmt='(A)') '#plt.legend(ylabel,loc=7,ncol=2,bbox_to_anchor=(1.1, 0.5))'
      write(fileNum,fmt='(A)') 'plt.show()'
      close(fileNum)
    enddo


  elseif(outputFormat=="gnp") then

    do nv=1,nVars
      !Make the same thing as CSV but using gnuplot file.
      fileName=trim(evaluatePrefix)//"-"//varName(nv)//".RMSE.gnuplot"
      !Open the output file
      call rams_f_open(fileNum,fileName,'FORMATTED','REPLACE','WRITE',1)
      write(fileNum,fmt='("#time",'//cLev//'(" L",I2.2))') ((z), z=1, nnzp(1))
      do nt=1,nTimes
        write(filenum,fmt='(F12.1,1X,'//cLev//'(" ",E9.2))') timeOfError(nt),statistic(nt,nv)%rmse
      enddo
      close(fileNum)
  
      fileName=trim(evaluatePrefix)//"-"//varName(nv)//".BIAS.gnuplot"
      !Open the output file
      call rams_f_open(fileNum,fileName,'FORMATTED','REPLACE','WRITE',1)
      write(fileNum,fmt='("#time",'//cLev//'(" L",I2.2))') ((z), z=1, nnzp(1))
      do nt=1,nTimes
        write(filenum,fmt='(F12.1,1X,'//cLev//'(" ",E9.2))') timeOfError(nt),statistic(nt,nv)%bias
      enddo
      close(fileNum)
    enddo
    ! If gnuplot write a file with gnuplot commands to generate graphs from data
    ! OBS: need gnuplot instalation to run the script.
    do nv=1,nVars
      filename=trim(evaluatePrefix)//"-"//varName(nv)//".RMSE.gnp"
      call rams_f_open(fileNum,fileName,'FORMATTED','REPLACE','WRITE',1)
      do nLev=1,nnzp(1)
        write(cLev,fmt='(I3.3)') nLev
        write(cPos,fmt='(I3.3)') nLev+1
        if(nLev==1) then
          write(fileNum,fmt='(A,I3,A)') 'plot "'//trim(evaluatePrefix)//"-"//varName(nv)// &
                                  '.RMSE.gnuplot" using 1:',nLev+1,' t"'//cLev//'" with lines'
        else
          write(fileNum,fmt='(A,I3,A)') 'rep  "'//trim(evaluatePrefix)//"-"//varName(nv)// &
                                  '.RMSE.gnuplot" using 1:',nLev+1,' t"'//cLev//'" with lines'
        endif
      enddo
      write(fileNum,fmt='(A)') 'set xlabel "Time [s]"'
      write(fileNum,fmt='(A)') 'set ylabel "RMSE"'
      write(fileNum,fmt='(A,I3.3,A)') 'set title "'//varName(nv)//' \n Levels: 1 to ',nnzp(1),'"'
      write(fileNum,fmt='(A)') 'set mxtics 21600'
      write(fileNum,fmt='(A)') 'set isosamples 100,100'
      write(fileNum,fmt='(A)') 'set output "'//trim(evaluatePrefix)//"-"//varName(nv)//'.RMSE.ps"'
      write(fileNum,fmt='(A)') 'set terminal postscript color landscape'
      write(fileNum,fmt='(A)') 'set noborder'
      write(fileNum,fmt='(A)') 'set grid'
      write(fileNum,fmt='(A)') 'replot'
      close(fileNum)
    enddo


    do nv=1,nVars
      filename=trim(evaluatePrefix)//"-"//varName(nv)//".BIAS.gnp"
      call rams_f_open(fileNum,fileName,'FORMATTED','REPLACE','WRITE',1)
      do nLev=1,nnzp(1)
        write(cLev,fmt='(I3.3)') nLev
        write(cPos,fmt='(I3.3)') nLev+1
        if(nLev==1) then
          write(fileNum,fmt='(A,I3,A)') 'plot "'//trim(evaluatePrefix)//"-"//varName(nv)// &
                                  '.BIAS.gnuplot" using 1:',nLev+1,' t"'//cLev//'" with lines'
        else
          write(fileNum,fmt='(A,I3,A)') 'rep  "'//trim(evaluatePrefix)//"-"//varName(nv)// &
                                  '.BIAS.gnuplot" using 1:',nLev+1,' t"'//cLev//'" with lines'
        endif
      enddo
      write(fileNum,fmt='(A)') 'set xlabel "Time [s]"'
      write(fileNum,fmt='(A)') 'set ylabel "BIAS"'
      write(fileNum,fmt='(A,I3.3,A)') 'set title "'//varName(nv)//' \n Levels: 1 to ',nnzp(1),'"'
      write(fileNum,fmt='(A)') 'set mxtics 21600'
      write(fileNum,fmt='(A)') 'set isosamples 100,100'
      write(fileNum,fmt='(A)') 'set output "'//trim(evaluatePrefix)//"-"//varName(nv)//'.BIAS.ps"'
      write(fileNum,fmt='(A)') 'set terminal postscript color landscape'
      write(fileNum,fmt='(A)') 'set noborder'
      write(fileNum,fmt='(A)') 'set grid'
      write(fileNum,fmt='(A)') 'replot'
      close(fileNum)
    enddo


  endif

end subroutine writeStatistics


subroutine StoreNamelistFileAtEvaluate(oneNamelistFile)
  type(namelistFile), pointer :: oneNamelistFile

  evaluatePrefix = oneNamelistFile%evaluatePrefix
  evaluate = oneNamelistFile%evaluate

end subroutine StoreNamelistFileAtEvaluate

end module ModEvaluation
