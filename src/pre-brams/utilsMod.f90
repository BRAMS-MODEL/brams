!=============================================================================================
module utilsMod
    !# Module with a sort of utilities
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: utilities are
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 26 August 2020 (Wednesday)
    !# @endnote
    !#
    !# @changes
    !# &#9744; <br/>
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
    
    !Use area
    use dump

    implicit none

    include "constants.f90"
    character(len=*),parameter :: sourceName='utils.f90' !Name of source code
    character(len=*),parameter :: procedureName='**utils**' !Name of this procedure
    !
    !Local Parameters
    integer, parameter :: i8 = selected_int_kind(14)   !Kind for 64-bits Integer Numbers
    integer, parameter :: r8 = selected_real_kind(15)  !Kind for 64-bits Real Numbers
    !Local variables

    logical :: fileUnits(20:99)
    logical :: firstTime

    Contains

    !=============================================================================================
    integer function initAll()
        !# Initialize utils
        !#
        !# @note
        !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
        !#
        !# **Brief**: initialize utils
        !#
        !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
        !#
        !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
        !#
        !# **Date**: 26 August 2020 (Wednesday)
        !# @endnote
        !#
        !# @changes
        !# &#9744; <br/>
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
        
        !Use area
        use dump
    
        implicit none
    
        include "constants.f90"
        character(len=*),parameter :: sourceName='utils.f90' !Name of source code
        character(len=*),parameter :: procedureName='**initAll**' !Name of this procedure
        !
        !Local Parameters
    
        !Input/Output variables
    
        !Local variables
    
        !Code
        fileUnits=.false.
        firstTime=.true.

        initAll=0
    
    end function initAll 

    !=============================================================================================
    integer function getUnit()
        !# Get a free unit to use
        !#
        !# @note
        !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
        !#
        !# **Brief**: Get a free unit to open file. 
        !#
        !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
        !#
        !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
        !#
        !# **Date**: 26 August 2020 (Wednesday)
        !# @endnote
        !#
        !# @changes
        !# &#9744; <br/>
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
        
        !Use area
        use dump
    
        implicit none
    
        include "constants.f90"
        character(len=*),parameter :: sourceName='utils.f90' !Name of source code
        character(len=*),parameter :: procedureName='**getUnit**' !Name of this procedure
        !
        !Local Parameters
    
        !Input/Output variables
    
        !Local variables
        integer :: i
    
        !Code
        do i=20,99
            if(.not. fileUnits(i)) then
               getUnit=i
               fileUnits(i)=.true.
               exit
            endif
         enddo

    
    end function getUnit 

    !=============================================================================================
    integer function releaseUnit(unitNum)
        !# Release a unit for file
        !#
        !# @note
        !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
        !#
        !# **Brief**: release a unit for file and close the file unitNum
        !#
        !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
        !#
        !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
        !#
        !# **Date**: 26 August 2020 (Wednesday)
        !# @endnote
        !#
        !# @changes
        !# &#9744; <br/>
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
        
        !Use area
        use dump
    
        implicit none
    
        include "constants.f90"
        character(len=*),parameter :: sourceName='utils.f90' !Name of source code
        character(len=*),parameter :: procedureName='**releaseUnit**' !Name of this procedure
        !
        !Local Parameters
    
        !Input/Output variables
        integer, intent(in) :: unitNum
    
        !Local variables
    
        !Code
        if(.not. fileUnits(unitNum)) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
              ,c_fatal,'Unit not used. Please, verify and solve it!',unitNum,"I2.2")
        close(unitNum)
        fileUnits(unitNum)=.false.
    
    end function releaseUnit 


    !=============================================================================================
    integer function bramsHeader(rev) 
        !# write a header in screen
        !#
        !# @note
        !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
        !#
        !# **Brief**: write a header in screen
        !#
        !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
        !#
        !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
        !#
        !# **Date**: 26 August 2020 (Wednesday)
        !# @endnote
        !#
        !# @changes
        !# &#9744; <br/>
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
        
        !Use area
        use dump
    
        implicit none
    
        include "constants.f90"
        character(len=*),parameter :: sourceName='utils.f90' !Name of source code
        character(len=*),parameter :: procedureName='**bramsHeader**' !Name of this procedure
        !
        !Input/Output variables
        character(len=*),intent(in) :: rev
    
        !Code
         write (*,fmt='(A)') ''
         write (*,fmt='(A)') '            ######  ######     #    #     #  #####'
         write (*,fmt='(A)') '            #     # #     #   # #   ##   ## #     #'
         write (*,fmt='(A)') '            #     # #     #  #   #  # # # # #'
         write (*,fmt='(A)') '            ######  ######  #     # #  #  #  #####'
         write (*,fmt='(A)') '            #     # #   #   ####### #     #       #'
         write (*,fmt='(A)') '            #     # #    #  #     # #     # #     #'
         write (*,fmt='(A)') '            ######  #     # #     # #     #  #####'
         write (*,fmt='(A)') '------------------------------------------------------------------'
         write (*,fmt='(A)') 'Brazilian developments on the Regional Atmospheric Modeling System'
         write (*,fmt='(A)') '                   PREP CHEM APP - Rev. '//rev
         write (*,fmt='(A)') '------------------------------------------------------------------' 
         write (*,fmt='(A)') ''

         bramsHeader=0
    
    end function bramsHeader 

    !=============================================================================================
    function to_upper(strIn) result(strOut)
        !# Convert case to upper case
        !#
        !# @note
        !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
        !#
        !# **Brief**: Convert case to upper
        !#
        !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
        !#
        !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
        !#
        !# **Date**: 26 August 2020 (Wednesday)
        !# @endnote
        !#
        !# @changes
        !# &#9744; <br/>
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
        
        !Use area
        use dump
    
        implicit none
    
        include "constants.f90"
        character(len=*),parameter :: sourceName='utils.f90' !Name of source code
        character(len=*),parameter :: procedureName='**to_upper**' !Name of this procedure
        !
        !Local Parameters
    
        !Input/Output variables
        character(*), intent(in) :: strIn
        !# String to be converted
        character(len=len(strIn)) :: strOut
        !# Return string converted
     
        !Local variables
        integer :: i
        !Code
        do i = 1, len(strIn)
            select case(strIn(i:i))
            case("a":"z")
               strOut(i:i) = achar(iachar(strIn(i:i))-32)
            end select
         end do
    
    end function to_upper 

    !=============================================================================================
    function to_lower(strIn) result(strOut)
        !# Convert case to lower case
        !#
        !# @note
        !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
        !#
        !# **Brief**: Convert case to lower case
        !#
        !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
        !#
        !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
        !#
        !# **Date**: 26 August 2020 (Wednesday)
        !# @endnote
        !#
        !# @changes
        !# &#9744; <br/>
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
        
        !Use area
        use dump
    
        implicit none
    
        include "constants.f90"
        character(len=*),parameter :: sourceName='utils.f90' !Name of source code
        character(len=*),parameter :: procedureName='**to_lower**' !Name of this procedure
        !
        !Local Parameters
    
        !Input/Output variables
        character(*), intent(in) :: strIn
        !# String to be converted
        character(len=len(strIn)) :: strOut
        !# Return string converted
     
        !Local variables
        integer :: i
        !Code
        do i = 1, len(strIn)
            select case(strIn(i:i))
            case("A":"Z")
               strOut(i:i) = achar(iachar(strIn(i:i))+32)
            end select
         end do
    
    end function to_lower

   !=============================================================================================
   integer function julday (imonth,iday,iyear)
      !# returns which day of the year is the input date
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: returns which day of the year is the input date
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 26 August 2020 (Wednesday)
      !# @endnote
      !#
      !# @changes
      !# &#9744; <br/>
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
        
      ! 

      integer, intent(in) :: imonth
      integer, intent(in) :: iday
      integer, intent(in) :: iyear
   
      julday= iday  &
           + min(1,max(0,imonth-1))*31  &
           + min(1,max(0,imonth-2))*(28+(1-min(1,mod(iyear,4))))  &
           + min(1,max(0,imonth-3))*31  &
           + min(1,max(0,imonth-4))*30  &
           + min(1,max(0,imonth-5))*31  &
           + min(1,max(0,imonth-6))*30  &
           + min(1,max(0,imonth-7))*31  &
           + min(1,max(0,imonth-8))*31  &
           + min(1,max(0,imonth-9))*30  &
           + min(1,max(0,imonth-10))*31  &
           + min(1,max(0,imonth-11))*30  &
           + min(1,max(0,imonth-12))*31
   
   end function julday

   !=============================================================================================
   subroutine date_make_big (inyear,inmonth,indate,inhour,outdate)
      !# convert integers year, month, date and hour into YYYYMMDDHHHHHH
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: convert integers representing year, month, date and hour into
      !#  a character string YYYYMMDDHHHHHH
      !# input hour is an integer with 6 digits in the form HHMMSS(hour, minute, second)
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 26 August 2020 (Wednesday)
      !# @endnote
      !#
      !# @changes
      !# &#9744; <br/>
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

      integer,           intent(in ) :: inyear
      integer,           intent(in ) :: inmonth
      integer,           intent(in ) :: indate
      integer,           intent(in ) :: inhour
      character(len=14), intent(out) :: outdate
   
      write(outdate(1:4), "(i4.4)") inyear
      write(outdate(5:6), "(i2.2)") inmonth
      write(outdate(7:8), "(i2.2)") indate
      write(outdate(9:14),"(i6.6)") inhour
   
   end subroutine date_make_big

   !=============================================================================================
   subroutine date_unmake_big (inyear,inmonth,indate,inhour,outdate)
      !# Convert a character string YYYYMMDDHHHHHH into integers
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: convert a character string YYYYMMDDHHHHHH into integers
      !# representing year, month, date and hour.
      !# output hour is an integer with 6 digits in the form HHMMSS
      !# (hour, minute, second)
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 26 August 2020 (Wednesday)
      !# @endnote
      !#
      !# @changes
      !# &#9744; <br/>
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

      integer,           intent(out) :: inyear
      integer,           intent(out) :: inmonth
      integer,           intent(out) :: indate
      integer,           intent(out) :: inhour
      character(len=14), intent(in ) :: outdate
   
      read(outdate(1:4), "(i4)") inyear
      read(outdate(5:6), "(i2)") inmonth
      read(outdate(7:8), "(i2)") indate
      read(outdate(9:14),"(i6)") inhour
   
   end subroutine date_unmake_big

   !=============================================================================================
   subroutine date_abs_secs(indate1, seconds)
      !# compute number of seconds past 1 January 1900 12:00 am to string
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: compute number of seconds past 1 January 1900 12:00 am
      !# from an input string in the form YYYYMMDDHHHHHH
      !# returns a real(kind=r8)
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 26 August 2020 (Wednesday)
      !# @endnote
      !#
      !# @changes
      !# &#9744; <br/>
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

      use dump, only: &
         dumpMessage

      implicit none
   
      include "constants.f90"
      ! 
   
      character(len=14), intent(in ) :: indate1
      real(kind=r8),     intent(out) :: seconds
   
      real(kind=r8) :: s1, s2, s3, s4
      integer :: year1, month1, date1, hour1, iy, ndays
      character(len=8) :: c0
      character(len=*), parameter :: h="**(date_abs_secs)**"
   
      call date_unmake_big(year1, month1, date1, hour1, indate1)
   
      if (year1 < 1900) then
         write(c0,"(i8)") year1
         !call fatal_error(h//" input year should be <= 1970; it was "&
         !      &//trim(adjustl(c0)))
         iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
                   " input year should be <= 1970; it was "&
                        &//trim(adjustl(c0)))
      end if
   
      iy      = year1 - 1900
      ndays   = iy*365 + max(0,(iy-1)/4) + julday(month1,date1,iy)
      s1      = dble(ndays)*86400.
      s2      = dble(hour1/10000)*3600.
      s3      = dble(mod(hour1,10000)/100)*60.
      s4      = dble(mod(hour1,100))
      seconds = s1 + s2 + s3 + s4
  
   end subroutine date_abs_secs

   !=============================================================================================
   subroutine date_abs_secs2 (year1,month1,date1,hour1,seconds)
      !# compute number of seconds past 1 January 1900 12:00 am to integer
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: compute number of seconds past 1 January 1900 12:00 am
      !# from integers representing year, month, date and hour
      !# returns a real(kind=r8)
      !# input hour is an integer with 6 digits in the form HHMMSS
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 26 August 2020 (Wednesday)
      !# @endnote
      !#
      !# @changes
      !# &#9744; <br/>
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

      use dump, only: &
         dumpMessage

      implicit none
   
      include "constants.f90"
      ! 
   
      integer,       intent(in ) :: year1
      integer,       intent(in ) :: month1
      integer,       intent(in ) :: date1
      integer,       intent(in ) :: hour1
      real(kind=r8), intent(out) :: seconds
   
      real(kind=r8) :: s1,s2,s3,s4
      integer :: iy,ndays
      character(len=8) :: c0
      character(len=*), parameter :: h="**(date_abs_secs2)**"
   
      if (year1 < 1900) then
         write(c0,"(i8)") year1
         !call fatal_error(h//" input year should be <= 1970; it was "&
        !      &//trim(adjustl(c0)))
         iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
                        " input year should be <= 1970; it was "&
                             &//trim(adjustl(c0)))
      end if
   
      iy = year1 - 1900
      ndays = iy * 365 + max(0,(iy-1)/4) + julday(month1,date1,iy)
      s1= dble(ndays) *86400.
      s2= dble(hour1/10000)*3600.
      s3= dble(mod(hour1,10000)/100)*60.
      s4= dble(mod(hour1,100))
      seconds= s1+s2+s3+s4
   
   end subroutine date_abs_secs2

   !=============================================================================================
   subroutine date_secs_ymdt (seconds,iyear1,imonth1,idate1,ihour1)
      !# compute real time given number of seconds past 1 January 1900 12:00 am
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: compute real time given number of seconds past 1 January 1900 12:00 am
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 26 August 2020 (Wednesday)
      !# @endnote
      !#
      !# @changes
      !# &#9744; <br/>
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
      real(kind=r8), intent(in ) :: seconds
      integer,       intent(out) :: iyear1
      integer,       intent(out) :: imonth1
      integer,       intent(out) :: idate1
      integer,       intent(out) :: ihour1

      integer, parameter :: mondays(12)=&
           (/31,28,31,30,31,30,31,31,30,31,30,31/)
      integer :: ny,nyr,ileap,nm,nd,ihr,imn,isc
      real(kind=r8) :: s1
   
      ! Get what year it is
   
      s1=seconds
      do ny=0,10000
         ileap=0
         if(mod(1900+ny,4) == 0) ileap=1
         s1=s1-(365.+ileap)*86400.
         if(s1 < 0.) then
            nyr=ny
            s1=s1+(365.+ileap)*86400.
            exit
         endif
      enddo
      iyear1=1900+nyr
   
      ! s1 is now number of secs into the year
      !   Get month
   
      do nm=1,12
         ileap=0
         if(mod(1900+ny,4) == 0 .and. nm == 2) ileap=1
         s1=s1-(mondays(nm)+ileap)*86400.
         if(s1 < 0.) then
            s1=s1+(mondays(nm)+ileap)*86400.
            exit
         endif
      enddo
      imonth1=nm
   
      ! s1 is now number of secs into the month
      !   Get date and time
   
      idate1=int(s1/86400.)
      s1=s1-idate1*86400.
      idate1=idate1+1 ! Since date starts at 1
   
      ihr=int(s1/3600.)
      s1=s1-ihr*3600.
      imn=int(s1/60.)
      s1=s1-imn*60.
      isc=s1
      ihour1=ihr*10000+imn*100+isc

   end subroutine date_secs_ymdt

   !=============================================================================================
   subroutine date_add_to_big (cindate,tinc,tunits,coutdate)
      !# adds/subtracts a time increment to a date and output new date
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: adds/subtracts a time increment to a date and output new date
      !#   uses hhmmss for hours, 4 digit year
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 26 August 2020 (Wednesday)
      !# @endnote
      !#
      !# @changes
      !# &#9744; <br/>
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
      character(len=14), intent(in ) :: cindate
      real,              intent(in ) :: tinc
      character(len=1),  intent(in ) :: tunits
      character(len=14), intent(out) :: coutdate
   
      real(kind=8) :: ttinc,secs
      integer :: inyear,inmonth,indate,inhour  &
           ,outyear,outmonth,outdate,outhour
   
      ! convert input increment to seconds
   
      select case(tunits)
      case("d","D")
         ttinc = tinc*86400.0
      case("h","H")
         ttinc = tinc*3600.0
      case("m","M")
         ttinc = tinc*60.0
      case default
         ttinc = tinc
      end select
   
      ! convert input time to seconds
   
      call date_unmake_big(inyear,inmonth,indate,inhour,cindate)
   
      call date_abs_secs2(inyear,inmonth,indate,inhour,secs)
   
      secs=secs+ttinc
   
      call date_secs_ymdt(secs,outyear,outmonth,outdate,outhour)
      call date_make_big(outyear,outmonth,outdate,outhour,coutdate)

   end subroutine date_add_to_big

   !=============================================================================================
   subroutine date_add_to (inyear,inmonth,indate,inhour,  &
       tinc,tunits,outyear,outmonth,outdate,outhour)
      !# convert input increment to seconds
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: convert input increment to seconds
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 26 August 2020 (Wednesday)
      !# @endnote
      !#
      !# @changes
      !# &#9744; <br/>
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
      integer,          intent(in ) :: inyear
      integer,          intent(in ) :: inmonth
      integer,          intent(in ) :: indate
      integer,          intent(in ) :: inhour
      real,             intent(in ) :: tinc
      character(len=1), intent(in ) :: tunits
      integer,          intent(out) :: outyear
      integer,          intent(out) :: outmonth
      integer,          intent(out) :: outdate
      integer,          intent(out) :: outhour
   
      real(kind=8) :: ttinc,secs
   
      ! convert input increment to seconds
   
      select case(tunits)
      case("d","D")
         ttinc = tinc*86400.0
      case("h","H")
         ttinc = tinc*3600.0
      case("m","M")
         ttinc = tinc*60.0
      case default
         ttinc = tinc
      end select
   
      ! convert input time to seconds
   
      call date_abs_secs2(inyear,inmonth,indate,inhour,secs)
   
      ! add increment
   
      secs=secs+ttinc
   
      ! convert seconds into date
   
      call date_secs_ymdt(secs,outyear,outmonth,outdate,outhour)
  
   end subroutine date_add_to

   !=============================================================================================
   subroutine date_add_to_dble (inyear,inmonth,indate,inhour,  &
       tinc,tunits,outyear,outmonth,outdate,outhour)
      !# convert input increment to seconds
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: convert input increment to seconds
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 26 August 2020 (Wednesday)
      !# @endnote
      !#
      !# @changes
      !# &#9744; <br/>
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
      integer,          intent(in ) :: inyear
      integer,          intent(in ) :: inmonth
      integer,          intent(in ) :: indate
      integer,          intent(in ) :: inhour
      real(kind=8),             intent(in ) :: tinc
      character(len=1), intent(in ) :: tunits
      integer,          intent(out) :: outyear
      integer,          intent(out) :: outmonth
      integer,          intent(out) :: outdate
      integer,          intent(out) :: outhour
   
      real(kind=8) :: ttinc,secs
   
      ! 
   
      select case(tunits)
      case("d","D")
         ttinc = tinc*86400.0
      case("h","H")
         ttinc = tinc*3600.0
      case("m","M")
         ttinc = tinc*60.0
      case default
         ttinc = tinc
      end select
   
      ! convert input time to seconds
   
      call date_abs_secs2(inyear,inmonth,indate,inhour,secs)
   
      ! add increment
   
      secs=secs+ttinc
   
      ! convert seconds into date
   
      call date_secs_ymdt(secs,outyear,outmonth,outdate,outhour)
   
   end subroutine date_add_to_dble

   !=============================================================================================
   logical function fileExist(fileName)
       !# Check if fileName exist
       !#
       !# @note
       !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
       !#
       !# **Brief**: Check if fileName exist
       !#
       !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
       !#
       !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
       !#
       !# **Date**: 26 August 2020 (Wednesday)
       !# @endnote
       !#
       !# @changes
       !# &#9744; <br/>
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
       
       !Use area
       use dump
   
       implicit none
   
       include "constants.f90"
       character(len=*),parameter :: sourceName='utils.f90' !Name of source code
       character(len=*),parameter :: procedureName='**fileExist**' !Name of this procedure
       !
       !Local Parameters
   
       !Input/Output variables
       character(len=*), intent(in) :: fileName
   
       !Local variables
   
       !Code

       inquire(file=fileName(1:len_trim(fileName)),exist=fileExist)

   
   end function fileExist 

   !=============================================================================================
   integer function stepsBetweenDates(iY,iM,iD,iH,fY,fM,fD,fH,step,cTime)
       !# Calc the number of steps between two dates
       !#
       !# @note
       !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
       !#
       !# **Brief**: Cal the numer of steps between two dates in intial and final
       !#
       !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
       !#
       !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
       !#
       !# **Date**: 26 August 2020 (Wednesday)
       !# @endnote
       !#
       !# @changes
       !# &#9744; <br/>
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
       
       !Use area
       use dump
   
       implicit none
   
       include "constants.f90"
       character(len=*),parameter :: sourceName='utils.f90' !Name of source code
       character(len=*),parameter :: procedureName='**stepsBetweenDates**' !Name of this procedure
       !
       !Local Parameters
   
       !Input/Output variables
       integer, intent(in) :: iY
       !# initial year
       integer, intent(in) :: iM
       !# initial month
       integer, intent(in) :: iD
       !# initial day
       integer, intent(in) :: iH
       !# initial hour
       integer, intent(in) :: fY
       !# final year
       integer, intent(in) :: fM
       !# final month
       integer, intent(in) :: fD
       !# final day
       integer, intent(in) :: fH
       !# final hour
       integer, intent(in) :: step
       !# step to be calulated
       character, intent(in) :: cTime
       !# unit of time (h,m,s)
       
       !Local variables
       integer :: iyy,imm,idd,ihh
   
       if(mod(fH,step)/=0) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
              ,c_fatal,' time difference isnt divided by step: ',step,"I2.2")

       stepsBetweenDates=0
       !Code
       do while(.true.)
         stepsBetweenDates=stepsBetweenDates+1
         call date_add_to_dble(iY,iM,iD,iH,dble(step)*(stepsBetweenDates-1),cTime &
                       ,iyy,imm,idd,ihh)

         if(iyy==fY .and. imm==fM .and. idd==fD .and. ihh/10000==fH) exit
       enddo 
   
   end function stepsBetweenDates 

   ! !=============================================================================================
   ! function validateDates(iY,iM,iD,iH,fY,fM,fD,fH,step,cTime,stepsBetDates) result(validDates)
   !     !# Validate dates inside given dates and step in hours
   !     !#
   !     !# @note
   !     !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
   !     !#
   !     !# **Brief**: validate date array inside two given dates and step hour
   !     !#
   !     !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
   !     !#
   !     !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
   !     !#
   !     !# **Date**: 28 August 2020 (Friday)
   !     !# @endnote
   !     !#
   !     !# @changes
   !     !# &#9744; <br/>
   !     !# @endchanges
   !     !# @bug
   !     !#
   !     !#@endbug
   !     !#
   !     !#@todo
   !     !#  &#9744; <br/>
   !     !# @endtodo
   !     !#
   !     !# @warning
   !     !# Now is under CC-GPL License, please see
   !     !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
   !     !# @endwarning
   !     !#
       
   !     !Use area
   !     use dump
   
   !     implicit none
   
   !     include "constants.f90"
   !     character(len=*),parameter :: procedureName='**validateDates**' !Name of this procedure
   !     !
   !     !Local Parameters
   
   !     !Input/Output variables
   !     integer, intent(in) :: iY
   !     !# initial year
   !     integer, intent(in) :: iM
   !     !# initial month
   !     integer, intent(in) :: iD
   !     !# initial day
   !     integer, intent(in) :: iH
   !     !# initial hour
   !     integer, intent(in) :: fY
   !     !# final year
   !     integer, intent(in) :: fM
   !     !# final month
   !     integer, intent(in) :: fD
   !     !# final day
   !     integer, intent(in) :: fH
   !     !# final hour
   !     integer, intent(in) :: step
   !     !# step to be calulated
   !     character, intent(in) :: cTime
   !     !# unit of time (h,m,s)  

   !     integer, intent(in) :: stepsBetDates
   !     !# Number of steps between dates 
   
   !     !Local variables
   !     integer :: validDates(stepsBetDates)
   
   !     !Code
   !     do while(.true.)
   !       sbd=sbd+1
   !       call date_add_to_dble(iY,iM,iD,iH,dble(step)*(sbd-1),cTime &
   !                     ,iyy,imm,idd,ihh)

   !       if(iyy==fY .and. imm==fM .and. idd==fD .and. ihh/10000==fH) exit
   !       lmonth(imm)=.true.

   !     enddo 
   
   ! end function validateDates 

   !=============================================================================================
   function monthsBetweenDates(iY,iM,iD,iH,fY,fM,fD,fH,step,cTime) result(lmonth)
       !# Fill an array with months presents between two dates
       !#
       !# @note
       !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
       !#
       !# **Brief**: Fill an array with months presents between two dates
       !#
       !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
       !#
       !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
       !#
       !# **Date**: 27 August 2020 (Thursday)
       !# @endnote
       !#
       !# @changes
       !# &#9744; <br/>
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
       
       !Use area
       use dump
   
       implicit none
   
       include "constants.f90"
       character(len=*),parameter :: procedureName='**monthsBetweenDates**' !Name of this procedure
       !
       !Local Parameters
   
       !Input/Output variables
       integer, intent(in) :: iY
       !# initial year
       integer, intent(in) :: iM
       !# initial month
       integer, intent(in) :: iD
       !# initial day
       integer, intent(in) :: iH
       !# initial hour
       integer, intent(in) :: fY
       !# final year
       integer, intent(in) :: fM
       !# final month
       integer, intent(in) :: fD
       !# final day
       integer, intent(in) :: fH
       !# final hour
       integer, intent(in) :: step
       !# step to be calulated
       character, intent(in) :: cTime
       !# unit of time (h,m,s)  

       !Local variables
       integer :: iyy,imm,idd,ihh
       integer :: sbd
       logical :: lMonth(12)
   
       !Code
       if(mod(fH,step)/=0) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
              ,c_fatal,' time difference isnt divided by step: ',step,"I2.2")
       lMonth=.false.

       sbd=0
       do while(.true.)
         sbd=sbd+1
         call date_add_to_dble(iY,iM,iD,iH,dble(step)*(sbd-1),cTime &
                       ,iyy,imm,idd,ihh)

         if(iyy==fY .and. imm==fM .and. idd==fD .and. ihh/10000==fH) exit
         lmonth(imm)=.true.

       enddo 

   end function monthsBetweenDates 

   !=============================================================================================
   subroutine interpolationBilinear(iValues, iNX, iNY, iNZ, iLat, iLon, iXS, iYS, oValues, oNX, oNY, oNZ, oLat, oLon, oXS, oYS)
       !# do a interpolatino bilinear of two data
       !#
       !# @note
       !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
       !#
       !# **Brief**: do a interpolation bilinear of two data
       !#
       !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
       !#
       !# **Author(s)**: Rafael Mello  **&#9993;**<rafael.mello@inpe.br>
       !#
       !# **Date**: 28 August 2020 (Friday)
       !# @endnote
       !#
       !# @changes
       !# &#9744; Modified by Luiz Flávio Rodrigues &#9993;**<luiz.rodrigues@inpe.br><br/>
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
       
       !Use area
       use dump
   
       implicit none
   
       include "constants.f90"
       character(len=*),parameter :: procedureName='**interpolationBilinear**' !Name of this procedure
       !
       !Local Parameters
       integer, parameter :: localundef = -999

       !Input/Output variables
       integer, intent(in) :: inx
       !#
       integer, intent(in) :: iny
       !#
       integer, intent(in) :: inz
       !#
       integer, intent(in) :: onx
       !#
       integer, intent(in) :: ony
       !#
       integer, intent(in) :: onz
       !#
       real, intent(in) :: ilat(2)
       !#
       real, intent(in) :: ilon(2)
       !#
       real, intent(in) :: olat(2)
       !#
       real, intent(in) :: olon(2)
       !#
       real, intent(in) :: ixs
       !#
       real, intent(in) :: iys
       !#
       real, intent(in) :: oxs
       !#
       real, intent(in) :: oys
       !#
       real, intent(in) :: ivalues(inx,iny,inz)
       !#
       real, intent(out) :: ovalues(onx,ony,onz)
    
      !Local variables
      integer :: i,j,z,px, py,idxlon(1)
      integer :: idxlat(1),signx,signy
      real    :: lon,lat,tmplat(2),tmplon(2)
      real    :: x0, x1, x2, y0, y1, y2, f11, f12, f21, f22, fr1, fr2, p
      real, allocatable :: latrepro(:),lonrepro(:),latorig(:),lonorig(:)
      integer, allocatable :: indexreprox(:,:),indexreproy(:,:)
   
   
       !Code
       if(.not. allocated(indexreprox))then  
          allocate(indexreprox(2,onx), indexreproy(2,ony))
          allocate(latrepro(ony), lonrepro(onx))
          allocate(latorig(iny), lonorig(inx))
          
          indexreprox = localundef
          indexreproy = localundef
          
          !creates a original data lon/lat vector
          do i = 1, inx
             if(ilon(1) .lt. ilon(2))then
                lon = ilon(1) + (ixs*(i-1))
             else
                lon = ilon(1) - (ixs*(i-1))                  
             end if
             lon = lon - (min(1,int(lon/181))*360)
             lonorig(i) = lon
          end do
          do j = 1, iny
             if(ilat(1) .lt. ilat(2))then
                lat = ilat(1) + (iys*(j-1))
             else
                lat = ilat(1) - (iys*(j-1))                  
             end if
             lat = lat - (min(1,int(lat/91))*180)
             latorig(j) = lat
          end do            
          
          !find max/min bounds in orig lat/lon vector
          !trying to avoid problems in reprojection when the dimension ends in values greater than 180 degrees.
          !0.......180.....360 => 181.....359.0......180
          idxlat = minloc(latorig)
          tmplat(1) = minval(latorig)
          tmplat(2) = maxval(latorig)  
          idxlon = minloc(lonorig)
          tmplon(1) = minval(lonorig)
          tmplon(2) = maxval(lonorig)            
          signx = 1
          signy = 1
          
          !if lon(1)/lat(1) are not the lowest/left bound then change the direction
          if(ilon(1) .gt. ilon(2)) signx = -1               
          if(ilat(1) .gt. ilat(2)) signy = -1
          
          !creates a matrix that defines the relation between original and reprojected data.             
          do i = 1, onx             
             lon = olon(1) + (oxs*(i-1))
             if(lon .ge. tmplon(1) .and. lon .le. tmplon(2))then
                px = ((lon - tmplon(1))/ixs)
             else
                px = -1
             end if
             lonrepro(i) = lon
             if(px .ge. 0)then
                px = idxlon(1) + (int(px) * signx)
                if(px .gt. inx) px = mod(px, inx)
             end if      
             if(px .gt. 0 .and. px .le. inx)                indexreprox(1,i) = px
             if((px+signx) .gt. 0 .and. (px+signx) .le. inx)indexreprox(2,i) = px + signx
          end do
          do j = 1, ony            
             lat = olat(1) + (oys*(j-1))
             if(lat .ge. tmplat(1) .and. lat .le. tmplat(2))then
                py = ((lat - tmplat(1))/iys)
             else
                py = -1
             end if
             latrepro(j) = lat
             if(py .ge. 0)then
                py = idxlat(1) + (int(py) * signy)
                if(py .gt. iny) py = mod(py, iny)
             end if
             if(py .gt. 0 .and. py .le. iny)                indexreproy(1,j) = py      
             if((py+signy) .gt. 0 .and. (py+signy) .le. iny)indexreproy(2,j) = py + signy
          end do            
       end if         
         
       do z = 1, onz
         do j = 1, ony
            do i = 1, onx               
               if(indexreprox(1,i) .ne. localundef .and. indexreproy(1,j) .ne. localundef)then
                  if(indexreprox(2,i) .ne. localundef .and. indexreproy(2,j) .ne. localundef)then
                     x0 = lonrepro(i)
                     y0 = latrepro(j)
                     x1 = lonorig(indexreprox(1,i))
                     y1 = latorig(indexreproy(1,j))
                     x2 = lonorig(indexreprox(2,i))
                     y2 = latorig(indexreproy(2,j))
                     f11 = ivalues(indexreprox(1,i), indexreproy(1,j), z)
                     f12 = ivalues(indexreprox(1,i), indexreproy(2,j), z)
                     f21 = ivalues(indexreprox(2,i), indexreproy(1,j), z)
                     f22 = ivalues(indexreprox(2,i), indexreproy(2,j), z)
                     if(f11 .ne. localundef .and. f12 .ne. localundef .and. &
                        f21 .ne. localundef .and. f22 .ne. localundef)then
                        fr1 = ( (((x2-x0)/(x2-x1))*f11) + (((x0-x1)/(x2-x1))*f21) )
                        fr2 = ( (((x2-x0)/(x2-x1))*f12) + (((x0-x1)/(x2-x1))*f22) )
                        p = ( ((y2-y0)/(y2-y1)*fr1) + ((y0-y1)/(y2-y1)*fr2) ) 
                        !nan check
                        if(p .ne. p) then
                           ovalues(i,j,z) = localundef
                        else
                           ovalues(i,j,z) = p
                        end if                           
                     end if
                  else
                     !no neighboors available
                     ovalues(i,j,z) = ivalues(indexreprox(1,i), indexreproy(1,j), z)
                  end if
               end if
            end do
         end do
      end do

      if(allocated(indexreprox))then
         deallocate(indexreprox, indexreproy)
         deallocate(latrepro,lonrepro)
         deallocate(latorig, lonorig)
      end if
   
   end subroutine interpolationBilinear 
  
   !=============================================================================================
   integer function outRealSize()
       !# get the output byte size accordingly the machine
       !#
       !# @note
       !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
       !#
       !# **Brief**: Return the output real size (1,4,8,etc) accordingly the machine
       !#
       !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
       !#
       !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
       !#
       !# **Date**: 28 July 2020 (Tuesday)
       !# @endnote
       !#
       !# @changes
       !# &#9744; <br/>
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
       
       !Use area
       use dump
   
       implicit none
   
       include "constants.f90"
       character(len=*),parameter :: sourceName='generic.f90' !Name of source code
       character(len=*),parameter :: procedureName='**getOutputByteSize**' !Name of this procedure
       !
       !Local Parameters
   
       !Input/Output variables
   
       !Local variables
       integer :: output_byte_size
       !# output_byte_size
       real :: bytes_in_float
       !# bytes_in_float
   
       !Code
       inquire(iolength=output_byte_size) bytes_in_float
   
       outRealSize=output_byte_size
   
   end function outRealSize 

   !=============================================================================================
   subroutine sortz(qtLevels, levels, pLevel, order, outIndex)
       !# A simple sort that returns the 'order' indexes nearest of 'plevel'
       !#
       !# @note
       !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
       !#
       !# **Brief**: simple sort that returns the 'order' indexes nearest of 'plevel' 
       !#
       !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
       !#
       !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
       !#
       !# **Date**: 01 September 2020 (Tuesday)
       !# @endnote
       !#
       !# @changes
       !# &#9744; <br/>
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
       
       !Use area
       use dump
   
       implicit none
   
       include "constants.f90"
       character(len=*),parameter :: procedureName='**sortz**' !Name of this procedure
       !
       !Local Parameters
   
       !Input/Output variables
       integer, intent(in) :: qtLevels
       !#
       real, intent(in)  :: levels(qtlevels)
       !#
       real, intent(in) :: pLevel
       !#
       integer, intent(in) :: order
       !#
       integer, intent(out) :: outIndex(order)
   
       !Local variables
       integer :: i,j,minidx
       real                      :: currmin
       real, dimension(qtLevels) :: ztemp
   
       !Code
       ztemp = levels
       do j = 1, order
          currMin = abs(ztemp(1) - pLevel)
          minIdx  = 1
          do i = 2, qtLevels
             if( abs(ztemp(i) - pLevel) .le. currMin)then
                currMin = abs(ztemp(i) - pLevel)
                minIdx  = i
             end if
          end do            
          outIndex(j) = minIdx
          ztemp(minIdx) = -9E5
       end do   
   
   end subroutine sortz 

end module utilsMod 