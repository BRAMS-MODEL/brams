module dump
  !# Auxiliary module to dump messages.
  !#
  !# @note
  !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
  !#
  !# **Brief**: The module can print 4 types of messages:
  !# 0) No Error message - just a print,
  !# 1) Notice,
  !# 2) warning or
  !# 3) Fatal.
  !# The interfaces selects the correct function and print the message using
  !# color or no. May be used simple messages only with text, text and
  !# integers, text and real, text and logical or combination of text and arrays
  !# 1D or 2D. The message may be printed inside a file or in ordinary TTY.
  !# See the examples of functions invokes:
  !#
  !# err=dumpMessage(c_tty,c_yes,header,version,c_notice,'Just a notice!')
  !# err=dumpMessage(c_tty,c_yes,header,version,c_warning,'Just a warning!')
  !# err=dumpMessage(c_tty,c_yes,header,version,c_warning,'Warning with an integer:',j,'I2')
  !# err=dumpMessage(c_tty,c_yes,header,version,c_notice,'Notice with a real:',x,'E18.4')
  !# err=dumpMessage(c_tty,c_yes,header,version,c_notice,'Do You love colors?',l1(2))
  !# err=dumpMessage(c_tty,c_yes,header,version,c_warning,'Warning with an array of reals: ',y,'F4.1')
  !# err=dumpMessage(c_tty,c_yes,header,version,c_notice,'Notice of 3 integers a,b,c: ',(/a,b,c/),'I1')
  !# err=dumpMessage(c_tty,c_yes,header,version,c_warning,'Array of Logicals: ',l1)
  !# err=dumpMessage(c_tty,c_yes,header,version,c_notice,'Notice with 2 d array:',a2d,'F3.1')
  !# err=dumpMessage(c_tty,c_no,header,version,c_notice,'Notice with 2 d array:',i2d,'I1')
  !# err=dumpMessage(c_tty,c_yes,header,version,c_warning,'Error in logical array:',l2d)
  !# err=dumpMessage(c_tty,c_yes,header,version,c_fatal,'Fatal - System crash!')
  !# err=dumpMessage(66,c_yes,c_empty,c_empty,c_noError,'Error in logical array:',l2d)
  !#
  !# c_tty, c_yes, c_no, c_notice, c_warning and c_fatal are parameters from contants.f90
  !#
  !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
  !#
  !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
  !#
  !# **Date**: 2015Sep
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

  interface dumpMessage
    module procedure dumpSingle
    module procedure dumpZeroInteger
    module procedure dumpZeroReal
    module procedure dump1DInteger
    module procedure dump1DReal
    module procedure dumpZeroLogical
    module procedure dump1DLogical
    module procedure dump2DReal
    module procedure dump2DInteger
    module procedure dump2DLogical
  end interface

  character(len=*), parameter :: noticeColor=achar(27)//'[97m'&
                                 //'Notice.! '//achar(27)//'[0m'
  !# Use to put the notice word in bold color
  character(len=*), parameter :: warningColor=achar(27)//'[97m'&
                                // 'Warning! '//achar(27)//'[0m'
  !# Use to put the warning word in bold color
  character(len=*), parameter :: fatalColor=achar(27)//'[97m'&
                                // 'Fatal..! '//achar(27)//'[0m'
  !# Use to put the fatal word in bold color

  private

  public :: dumpMessage, emptyLine, logical2Int,openLogFile

contains

  integer function openLogFile()
    include "constants.h"

    openLogFile=0
    open(unit=logUnit,file='bramsLog.out',ERR=100)
    return

100 openLogFile=-1

  end function openLogFile

  integer function logical2Int(logicalValue)
    implicit none

    logical,intent(in) :: logicalValue
    logical2Int=0
    if(logicalValue) logical2Int=1

  end function logical2Int

  integer function emptyLine(tty)
    implicit none

    include "constants.h"
    integer, intent(in) :: tty
    !# terminal/file to print

    write(tty,fmt='(A)') c_empty

    emptyLine=0

  end function emptyLine

  integer function dumpSingle(tty,color,header,version,dumpType,message)
    implicit none

    include "constants.h"

    integer, intent(in) :: tty
    !# terminal/file to print
    logical, intent(in) :: color
    !# if is to color message
    character(len=*), intent(in) :: header
    !# Name of procedure that invokes
    character(len=*), intent(in) :: version
    !# version of the procedure
    integer, intent(in) :: dumpType
    !# Type of dump message: notice, warning or fatal
    character(len=*), intent(in) :: message
    !# message text to print

    select case(dumpType)
    case(c_noError)
      write(tty,fmt='(A)') header//' '//version//' '//message
      dumpSingle=0
    case(c_notice)
      if(color) then
        write(tty,fmt='(A)') noticeColor//header//' ' //version//' '//achar(27)&
        //'[32m'//message//achar(27)//'[0m'
      else
        write(tty,fmt='(A)') 'Notice! '//header//' '//version//' '//message
      endif
      dumpSingle=1
    case(c_warning)
      if(color) then
        write(tty,fmt='(A)') warningColor//header//' '//version//' ' &
        //achar(27)//'[93m'//message//achar(27)//'[0m'
      else
        write(tty,fmt='(A)') 'Warning!! '//header//' '//version//' '//message
      endif
      dumpSingle=2
    case(c_fatal)
      if(color) then
        write(tty,fmt='(A)') fatalColor//header//' '//version//' '&
        //achar(27)//'[31m'//message//achar(27)//'[0m'
        write(tty,fmt='(A)') c_blink//'Fatal Error! See the message above!' &
                           //c_noColor
      else
        write(tty,fmt='(A)') 'Fatal!!! from '//header//' '//version//' '&
        //message
        write(tty,fmt='(A)') '!!!!!!! Fatal Error! See the message above !!!!!!!!'
      endif

      print *,c_empty
      stop
    end select

  end function dumpSingle

  integer function dumpZeroInteger(tty,color,header,version,dumpType,message &
                   ,value,cFormat)
    implicit none

    include "constants.h"

    integer, intent(in) :: tty
    !# terminal/file to print
    logical, intent(in) :: color
    !# if is to color message
    character(len=*), intent(in) :: header
    !# Name of procedure that invokes
    character(len=*), intent(in) :: version
    !# version of the procedure
    integer, intent(in) :: dumpType
    !# Type of dump message: notice, warning or fatal
    character(len=*), intent(in) :: message
    !# message text to print
    integer,intent(in) :: value
    !# Integer value to print
    character(len=*), intent(in) :: cFormat
    !# Format of number output

    select case(dumpType)
    case(c_noError)
      write(tty,fmt='(A,'//cFormat//')') header//' '//version&
              //' '//message//' ',value
      dumpZeroInteger=0
    case(c_notice)
      if(color) then
        write(tty,fmt='(A,'//cFormat//',A)') noticeColor//header//' '//version &
             //' '//achar(27)//'[32m'//message//achar(27)//'[0m'//' ' &
             //achar(27)//'[95m',value,achar(27)//'[0m'
      else
        write(tty,fmt='(A,'//cFormat//')') 'Notice! '//header//' '//version&
              //' '//message//' ',value
      endif
      dumpZeroInteger=1
    case(c_warning)
      if(color) then
        write(tty,fmt='(A,'//cFormat//',A)') warningColor//header//' '//version&
             //' '//achar(27)//'[93m'//message//achar(27)//'[0m'//' ' &
             //achar(27)//'[95m',value,achar(27)//'[0m'
      else
        write(tty,fmt='(A,'//cFormat//')') 'Notice! '//header//' '//version&
              //' '//message//' ',value
      endif
      dumpZeroInteger=2
    case(c_fatal)
      if(color) then
        write(tty,fmt='(A,'//cFormat//',A)') fatalColor//header//' '//version &
             //' '//achar(27)//'[31m'//message//achar(27)//'[0m'//' ' &
             //achar(27)//'[95m',value,achar(27)//'[0m'
      else
        write(tty,fmt='(A,'//cFormat//')') 'Fatal!!! from '//header//'|' &
             //version//' '//message//' ',value
      endif
      write(tty,fmt='(A)') c_blink//'Fatal Error! See the message above!'&
                         //c_noColor
      print *,c_empty
      stop
    end select

  end function dumpZeroInteger

  integer function dumpZeroReal(tty,color,header,version,dumpType,message &
                   ,value,cFormat)
    implicit none

    include "constants.h"

    integer, intent(in) :: tty
    !# terminal/file to print
    logical, intent(in) :: color
    !# if is to color message
    character(len=*), intent(in) :: header
    !# Name of procedure that invokes
    character(len=*), intent(in) :: version
    !# version of the procedure
    integer, intent(in) :: dumpType
    !# Type of dump message: notice, warning or fatal
    character(len=*), intent(in) :: message
    !# message text to print
    real,intent(in) :: value
    !# real value to print
    character(len=*), intent(in) :: cFormat
    !# Format of number output

    select case(dumpType)
    case(c_noError)
      write(tty,fmt='(A,'//cFormat//')') header//' '//version&
              //' '//message//' ',value
      dumpZeroReal=0
    case(c_notice)
      if(color) then
        write(tty,fmt='(A,'//cFormat//',A)') noticeColor//header//' '//version &
             //' '//achar(27)//'[32m'//message//achar(27)//'[0m'//' ' &
             //achar(27)//'[35m',value,achar(27)//'[0m'
      else
        write(tty,fmt='(A,'//cFormat//')') 'Notice! '//header//' '//version&
              //' '//message//' ',value
      endif
      dumpZeroReal=1
    case(c_warning)
      if(color) then
        write(tty,fmt='(A,'//cFormat//',A)') warningColor//header//' '//version&
             //' '//achar(27)//'[93m'//message//achar(27)//'[0m'//' ' &
             //achar(27)//'[35m',value,achar(27)//'[0m'
      else
        write(tty,fmt='(A,'//cFormat//')') 'Notice! '//header//' '//version&
              //' '//message//' ',value
      endif
      dumpZeroReal=2
    case(c_fatal)
      if(color) then
        write(tty,fmt='(A,'//cFormat//',A)') fatalColor//header//' '//version &
             //' '//achar(27)//'[31m'//message//achar(27)//'[0m'//' ' &
             //achar(27)//'[35m',value,achar(27)//'[0m'
      else
        write(tty,fmt='(A,'//cFormat//')') 'Fatal!!! from '//header//'|' &
             //version//' '//message//' ',value
      endif
      write(tty,fmt='(A)') c_blink//'Fatal Error! See the message above!'&
                         //c_noColor
      print *,c_empty
      stop
    end select

  end function dumpZeroReal

  integer function dump1DInteger(tty,color,header,version,dumpType,message &
                   ,value,cFormat)
    implicit none

    include "constants.h"

    integer, intent(in) :: tty
    !# terminal/file to print
    logical, intent(in) :: color
    !# if is to color message
    character(len=*), intent(in) :: header
    !# Name of procedure that invokes
    character(len=*), intent(in) :: version
    !# version of the procedure
    integer, intent(in) :: dumpType
    !# Type of dump message: notice, warning or fatal
    character(len=*), intent(in) :: message
    !# message text to print
    integer,intent(in) :: value(:)
    !# Integer value to print
    character(len=*), intent(in) :: cFormat
    !# Format of number output

    integer,allocatable :: value2d(:,:)
    allocate(value2d(size(value,1),1))

    value2d(:,1)=value(:)

    dump1DInteger=dump2dInteger(tty,color,header,version,dumpType,message &
                  ,value2d,cFormat)

  end function dump1DInteger

  integer function dump1DReal(tty,color,header,version,dumpType,message,value &
                   ,cFormat)
    implicit none

    include "constants.h"

    integer, intent(in) :: tty
    !# terminal/file to print
    logical, intent(in) :: color
    !# if is to color message
    character(len=*), intent(in) :: header
    !# Name of procedure that invokes
    character(len=*), intent(in) :: version
    !# version of the procedure
    integer, intent(in) :: dumpType
    !# Type of dump message: notice, warning or fatal
    character(len=*), intent(in) :: message
    !# message text to print
    real,intent(in) :: value(:)
    !# Integer value to print
    character(len=*), intent(in) :: cFormat
    !# Format of number output

    real,allocatable :: value2d(:,:)
    allocate(value2d(size(value,1),1))

    value2d(:,1)=value(:)

    dump1DReal=dump2dReal(tty,color,header,version,dumpType,message,value2d &
              ,cFormat)

  end function dump1DReal

  integer function dump1DLogical(tty,color,header,version,dumpType,message &
                   ,value)
    implicit none

    include "constants.h"

    integer, intent(in) :: tty
    !# terminal/file to print
    logical, intent(in) :: color
    !# if is to color message
    character(len=*), intent(in) :: header
    !# Name of procedure that invokes
    character(len=*), intent(in) :: version
    !# version of the procedure
    integer, intent(in) :: dumpType
    !# Type of dump message: notice, warning or fatal
    character(len=*), intent(in) :: message
    !# message text to print
    logical,intent(in) :: value(:)
    !# Integer value to print

    character(len=*), parameter :: cFormat='L1'
    !# Format of number output
    logical,allocatable :: value2d(:,:)
    !# var with one more dimension
    integer :: err

    allocate(value2d(size(value,1),1))

    value2d(:,1)=value(:)

    dump1DLogical=dump2DLogical(tty,color,header,version,dumpType,message &
                  ,value2d)

  end function dump1DLogical

  integer function dumpZeroLogical(tty,color,header,version,dumpType,message &
                   ,value)
    implicit none

    include "constants.h"

    integer, intent(in) :: tty
    !# terminal/file to print
    logical, intent(in) :: color
    !# if is to color message
    character(len=*), intent(in) :: header
    !# Name of procedure that invokes
    character(len=*), intent(in) :: version
    !# version of the procedure
    integer, intent(in) :: dumpType
    !# Type of dump message: notice, warning or fatal
    character(len=*), intent(in) :: message
    !# message text to print
    logical,intent(in) :: value
    !# real value to print

    character(len=*), parameter :: cFormat='L1'
    !# Format of number output

    select case(dumpType)
    case(c_noError)
      write(tty,fmt='(A,'//cFormat//')') header//' '//version&
              //' '//message//' ',value
      dumpZeroLogical=0
    case(c_notice)
      if(color) then
        write(tty,fmt='(A,'//cFormat//',A)') noticeColor//header//' '//version &
             //' '//achar(27)//'[32m'//message//achar(27)//'[0m'//' ' &
             //achar(27)//'[34m',value,achar(27)//'[0m'
      else
        write(tty,fmt='(A,'//cFormat//')') 'Notice! '//header//' '//version&
              //' '//message//' ',value
      endif
      dumpZeroLogical=1
    case(c_warning)
      if(color) then
        write(tty,fmt='(A,'//cFormat//',A)') warningColor//header//' '//version&
             //' '//achar(27)//'[93m'//message//achar(27)//'[0m'//' ' &
             //achar(27)//'[34m',value,achar(27)//'[0m'
      else
        write(tty,fmt='(A,'//cFormat//')') 'Warning! '//header//' '//version&
              //' '//message//' ',value
      endif
      dumpZeroLogical=2
    case(c_fatal)
      if(color) then
        write(tty,fmt='(A,'//cFormat//',A)') fatalColor//header//' '//version &
             //' '//achar(27)//'[31m'//message//achar(27)//'[0m'//' ' &
             //achar(27)//'[34m',value,achar(27)//'[0m'
      else
        write(tty,fmt='(A,'//cFormat//')') 'Fatal!!! from '//header//'|' &
             //version//' '//message//' ',value
      endif
      write(tty,fmt='(A)') c_blink//'Fatal Error! See the message above!'&
                         //c_noColor
      print *,c_empty
      stop
    end select

  end function dumpZeroLogical

  integer function dump2DReal(tty,color,header,version,dumpType,message,value &
                   ,cFormat)
    implicit none

    include "constants.h"

    integer, intent(in) :: tty
    !# terminal/file to print
    logical, intent(in) :: color
    !# if is to color message
    character(len=*), intent(in) :: header
    !# Name of procedure that invokes
    character(len=*), intent(in) :: version
    !# version of the procedure
    integer, intent(in) :: dumpType
    !# Type of dump message: notice, warning or fatal
    character(len=*), intent(in) :: message
    !# message text to print
    real,intent(in) :: value(:,:)
    !# real value to print
    character(len=*), intent(in) :: cFormat
    !# Format of number output


    character(len=4) :: c_sizeArray
    integer :: iCount

    write(c_sizeArray,fmt='(I4.4)') size(value,1)

    select case(dumpType)
    case(c_noError)
      write(tty,fmt='(A)') header//' '//version//' '//message
      dump2DReal=0
    case(c_notice)
      if(color) then
        write(tty,fmt='(A)') noticeColor//header//' ' //version//' '//achar(27)&
        //'[32m'//message//achar(27)//'[0m'
      else
        write(tty,fmt='(A)') 'Notice! '//header//' '//version//' '//message
      endif

      dump2DReal=1
    case(c_warning)
      if(color) then
        write(tty,fmt='(A)') warningColor//header//' '//version//' ' &
             //achar(27)//'[93m'//message//achar(27)//'[0m'
      else
        write(tty,fmt='(A)') 'Warning! '//header//' '//version//' '//message
      endif

      dump2DReal=2
    case(c_fatal)
      if(color) then
        write(tty,fmt='(A)') fatalColor//header//' ' //version//' '//achar(27) &
             //'[31m'//message//achar(27)//'[0m'
      else
        write(tty,fmt='(A)') 'Fatal!!! '//header//' '//version//' '//message
      endif

    end select
    !
    do iCount=1,size(value,2)
      write(tty,fmt='(A,'//c_sizeArray//'('//cFormat//',1X),A)') achar(27) &
           //'[35m',value(:,iCount),achar(27)//'[0m'
    enddo

    !
    if(dumpType==c_fatal) then
      write(tty,fmt='(A)') c_blink//'Fatal Error! See the message above!' &
                         //c_noColor
      print *,c_empty
      stop
    endif

  end function dump2DReal

  integer function dump2DInteger(tty,color,header,version,dumpType,message&
                   ,value,cFormat)
    implicit none

    include "constants.h"

    integer, intent(in) :: tty
    !# terminal/file to print
    logical, intent(in) :: color
    !# if is to color message
    character(len=*), intent(in) :: header
    !# Name of procedure that invokes
    character(len=*), intent(in) :: version
    !# version of the procedure
    integer, intent(in) :: dumpType
    !# Type of dump message: notice, warning or fatal
    character(len=*), intent(in) :: message
    !# message text to print
    integer,intent(in) :: value(:,:)
    !# real value to print
    character(len=*), intent(in) :: cFormat
    !# Format of number output


    character(len=4) :: c_sizeArray
    integer :: iCount

    write(c_sizeArray,fmt='(I4.4)') size(value,1)

    select case(dumpType)
    case(c_noError)
      write(tty,fmt='(A)') header//' '//version//' '//message
      dump2DInteger=0
    case(c_notice)
      if(color) then
        write(tty,fmt='(A)') noticeColor//header//' ' //version//' '//achar(27)&
        //'[32m'//message//achar(27)//'[0m'
      else
        write(tty,fmt='(A)') 'Notice! '//header//' '//version//' '//message
      endif

      dump2DInteger=1
    case(c_warning)
      if(color) then
        write(tty,fmt='(A)') warningColor//header//' '//version//' ' &
             //achar(27)//'[93m'//message//achar(27)//'[0m'
      else
        write(tty,fmt='(A)') 'Warning! '//header//' '//version//' '//message
      endif

      dump2DInteger=2
    case(c_fatal)
      if(color) then
        write(tty,fmt='(A)') fatalColor//header//' ' //version//' '//achar(27) &
             //'[31m'//message//achar(27)//'[0m'
      else
        write(tty,fmt='(A)') 'Fatal!!! '//header//' '//version//' '//message
      endif

    end select
    !
    do iCount=1,size(value,2)
      write(tty,fmt='(A,'//c_sizeArray//'('//cFormat//',1X),A)') achar(27) &
           //'[95m',value(:,iCount),achar(27)//'[0m'
    enddo
    !
    if(dumpType==c_fatal) then
      write(tty,fmt='(A)') c_blink//'Fatal Error! See the message above!' &
                         //c_noColor
      print *,c_empty
      stop
    endif

  end function dump2DInteger


  integer function dump2DLogical(tty,color,header,version,dumpType,message &
                  ,value)
    implicit none

    include "constants.h"

    integer, intent(in) :: tty
    !# terminal/file to print
    logical, intent(in) :: color
    !# if is to color message
    character(len=*), intent(in) :: header
    !# Name of procedure that invokes
    character(len=*), intent(in) :: version
    !# version of the procedure
    integer, intent(in) :: dumpType
    !# Type of dump message: notice, warning or fatal
    character(len=*), intent(in) :: message
    !# message text to print
    logical,intent(in) :: value(:,:)
    !# real value to print

    character(len=*), parameter :: cFormat='L1'
    !# Format of number output


    character(len=4) :: c_sizeArray
    integer :: iCount

    write(c_sizeArray,fmt='(I4.4)') size(value,1)

    select case(dumpType)
    case(c_noError)
      write(tty,fmt='(A)') header//' '//version//' '//message
      dump2DLogical=0
    case(c_notice)
      if(color) then
        write(tty,fmt='(A)') noticeColor//header//' ' //version//' '//achar(27)&
        //'[32m'//message//achar(27)//'[0m'
      else
        write(tty,fmt='(A)') 'Notice! '//header//' '//version//' '//message
      endif

      dump2DLogical=1
    case(c_warning)
      if(color) then
        write(tty,fmt='(A)') warningColor//header//' '//version//' ' &
             //achar(27)//'[93m'//message//achar(27)//'[0m'
      else
        write(tty,fmt='(A)') 'Warning! '//header//' '//version//' '//message
      endif

      dump2DLogical=2
    case(c_fatal)
      if(color) then
        write(tty,fmt='(A)') fatalColor//header//' ' //version//' '//achar(27) &
             //'[31m'//message//achar(27)//'[0m'
      else
        write(tty,fmt='(A)') 'Fatal!!! '//header//' '//version//' '//message
      endif

    end select
    !
    do iCount=1,size(value,2)
      write(tty,fmt='(A,'//c_sizeArray//'('//cFormat//',1X),A)') achar(27) &
           //'[34m',value(:,iCount),achar(27)//'[0m'
    enddo

    !
    if(dumpType==c_fatal) then
      write(tty,fmt='(A)') c_blink//'Fatal Error! See the message above!' &
                         //c_noColor
      print *,c_empty
      stop
    endif

  end function dump2DLogical

end module dump
