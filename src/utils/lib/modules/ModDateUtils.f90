module ModDateUtils
  implicit none

  ! provides gregorian calendar utilities

  ! gregorian dates are expressed in two forms:
  ! (1) a single character string YYYYMMDDHHMMSS
  ! (2) four integers year, month, day, hour(as HHMMSS)

  ! utilities for dates expressed as (1):
  ! date_add_to_big increments input date by given amount
  ! date_abs_secs converts input date into seconds since
  !               January 1st, 1900

  ! utilities for dates expressed as (2):
  ! date_add_to increments input date by given amount
  ! date_abs_secs2 converts input date into seconds since
  !                January 1st, 1900

  ! conversion utilities:
  ! date_make_big converts form (2) into form (1)
  ! date_unmake_big converts form (1) into form (2)
  ! julday returns which day of the year is the input date

  private
  public :: date_add_to
  public :: date_add_to_dble
  public :: date_add_to_big
  public :: date_make_big
  public :: date_unmake_big
  public :: julday
  public :: date_abs_secs
  public :: date_abs_secs2
  public :: date_secs_ymdt
  include "ranks.h"


contains





  integer function julday (imonth,iday,iyear)

    ! returns which day of the year is the input date

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





  subroutine date_make_big (inyear,inmonth,indate,inhour,outdate)

    ! convert integers representing year, month, date and hour into
    ! a character string YYYYMMDDHHHHHH
    ! input hour is an integer with 6 digits in the form HHMMSS
    ! (hour, minute, second)

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





  subroutine date_unmake_big (inyear,inmonth,indate,inhour,outdate)

    ! convert a character string YYYYMMDDHHHHHH into integers
    ! representing year, month, date and hour.
    ! output hour is an integer with 6 digits in the form HHMMSS
    ! (hour, minute, second)

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





  subroutine date_abs_secs(indate1, seconds)
    use dump, only: &
      dumpMessage

    implicit none

    include "constants.f90"
    ! compute number of seconds past 1 January 1900 12:00 am
    ! from an input string in the form YYYYMMDDHHHHHH
    ! returns a real(kind=r8)

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





  subroutine date_abs_secs2 (year1,month1,date1,hour1,seconds)
    use dump, only: &
      dumpMessage

    implicit none

    include "constants.f90"
    ! compute number of seconds past 1 January 1900 12:00 am
    ! from integers representing year, month, date and hour
    ! returns a real(kind=r8)
    ! input hour is an integer with 6 digits in the form HHMMSS

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





  subroutine date_secs_ymdt (seconds,iyear1,imonth1,idate1,ihour1)
    real(kind=r8), intent(in ) :: seconds
    integer,       intent(out) :: iyear1
    integer,       intent(out) :: imonth1
    integer,       intent(out) :: idate1
    integer,       intent(out) :: ihour1

    ! compute real time given number of seconds past 1 January 1900 12:00 am


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





  subroutine date_add_to_big (cindate,tinc,tunits,coutdate)
    character(len=14), intent(in ) :: cindate
    real,              intent(in ) :: tinc
    character(len=1),  intent(in ) :: tunits
    character(len=14), intent(out) :: coutdate

    ! adds/subtracts a time increment to a date and output new date
    ! -> uses hhmmss for hours, 4 digit year

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

    return
  end subroutine date_add_to_big






  subroutine date_add_to (inyear,inmonth,indate,inhour,  &
       tinc,tunits,outyear,outmonth,outdate,outhour)

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

    subroutine date_add_to_dble (inyear,inmonth,indate,inhour,  &
       tinc,tunits,outyear,outmonth,outdate,outhour)

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
  end subroutine date_add_to_dble

end module ModDateUtils
