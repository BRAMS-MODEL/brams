program main
  !# BRAMS Model Main program
  !#
  !# @note
  !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
  !#
  !# **Brief**:BRAMS Main Model - Based on RAMS (Regional Atmospheric Modeling System) & under permission of Atmet
  !# &copy; 1990, 1995, 1999, 2000, 2003 - All Rights Reserved
  !#
  !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
  !#
  !# **Author**: Drs. Tremback, Walko **&#9993;**<mailto:tremback@atmet.com>
  !#
  !# **Date**: 1990
  !# @endnote
  !#
  !# @changes
  !#
  !# + (1998Feb - By Dias, P.L. - USP) Created the BRAMS 1.0 (Brasilian version) of RAMS
  !# + (2002Feb - By Panetta. J. - CPTEC) BRAMS 2.0 - Development of Model by CPTEC/INPE
  !# + (2003Jan - By Panetta. J. - CPTEC) Versions 3.1 and 3.2 (Based on RAMS 5.0.4)
  !# + (2005Apr - CPTEC/IAG/UFCG/UFRJ/FURG) The BRAMSnet created to continuos Development of model
  !# + (2006Jul - By Freitas, S.R. & all - CPTEC) Version 4.0 - Included the Chemical Model (CCATT)
  !# + (2008Jan - By Panetta, J. & all - CPTEC) Version 4.2 - Changes/adapt for enhanced performance
  !# + (2010May - BRAMS' Team - CPTEC) Version 5.0 - Master node elimination - High scalability
  !# + (2011Jun - BRAMS'Team - CPTEC) Version 5.1 - Code modernization
  !# + (2015Mar - BRAMS'Team - CPTEC) Version 5.2 - New surface model JULES, Grell Freitas scheme and a lot of features
  !# + (2017Sep - BRAMS'Team - CPTEC) Version 5.3 - No use of HDF inputs, RRTMG radiation and a lot of features
  !#
  !# @endchanges
  !# @bug
  !# No active bugs reported now
  !# @endbug
  !#
  !# @todo
  !# 1.  &#9744; <br/>
  !# @endtodo
  !#
  !# @warning
  !# Now is under CC-GPL License, please see
  !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
  !# @endwarning
  !#
  !#--- ----------------------------------------------------------------------------------------
  use ModMemory, only: &
       CreateMemory,   & ! subroutine
       DumpMemory,     & ! subroutine
       DestroyMemory     ! subroutine

  use ModTimeStamp, only: &
       CreateTimeStamp,    & ! subroutine
       DestroyTimeStamp      ! subroutine

  use ParLib, only: &
       parf_init_mpi, & ! subroutine
       parf_barrier, &  ! subroutine
       parf_exit_mpi    ! subroutine

  use ModOneProc, only: &
       OneProc

  use dump, only: &
        dumpMessage ! Subroutine

  ! Main:
  !   starts MPI
  !   initializes memory use and execution time instrumentation
  !   dispatches processes
  !   invoking master/slave processes or full model process
  !   dumps and destroys instrumentation
  !   finishes MPI

  implicit none

  include "constants.f90"
  ! process rank and size (local variables)
  integer :: nmachs_in
  integer :: mchnum_in
  integer :: master_num_in
  character(len=*), parameter :: h="**(main)**"
  character(len=*), parameter :: header="**(main)**"
  character(len=*), parameter :: version="5.4"
  character(len=8) :: c0, c1

  ! execution time instrumentation

  include "tsNames.h"
  
  integer :: ierr

  ! enroll MPI; get number of processes, process Id, master Id

  call parf_init_mpi(mchnum_in, nmachs_in, master_num_in)

  ! initialize memory instrumentation

  if (mchnum_in==master_num_in) then
     call CreateMemory(mchnum_in, nmachs_in, master_num_in)
  end if

  ! initialize execution time instrumentation

  call CreateTimeStamp(nmachs_in, mchnum_in, names)

  ! dispatch processes

  call OneProc(nmachs_in, mchnum_in, master_num_in)

  ! finishes execution

  call parf_barrier(0)

  call DestroyTimeStamp()

  if (mchnum_in==master_num_in) then
     call DumpMemory("Fim")
     call DestroyMemory()
  end if

  if (mchnum_in==master_num_in) then!write(*,"(a)") " ****** BRAMS execution ends ******"
#ifdef color
    iErrNumber=dumpMessage(c_tty,c_yes,header,version,c_notice,'BRAMS execution normal ends!')
    iErrNumber=dumpMessage(c_tty,c_yes,header,version,c_notice &
         ,'for more information about submission(non fatal errors, notices, warnings), please, see the file brams.log and jules.log !')
#else
    iErrNumber=dumpMessage(c_tty,c_no,header,version,c_notice,'BRAMS execution normal ends!')
    iErrNumber=dumpMessage(c_tty,c_no,header,version,c_notice &
     ,'for more information about submission(non fatal errors, notices, warnings), please, see the file brams.log and jules.log  !')
#endif
  endif
  call parf_exit_mpi()
end program main
