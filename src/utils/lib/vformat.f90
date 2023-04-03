!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine vfinit()

  implicit none
  character(len=1) :: vc(0:63) !vc, vcscr(0:63)
  common/vform/vc !(0:63)
  integer :: n
!!$  data vcscr/&
!!$       '0','1','2','3','4','5','6','7','8','9',  &
!!$       'A','B','C','D','E','F','G','H','I','J',  &
!!$       'K','L','M','N','O','P','Q','R','S','T',  &
!!$       'U','V','W','X','Y','Z','a','b','c','d',  &
!!$       'e','f','g','h','i','j','k','l','m','n',  &
!!$       'o','p','q','r','s','t','u','v','w','x',  &
!!$       'y','z','{','|'/
  character(len=1), parameter :: vcscr(0:63) = (/&
       '0','1','2','3','4','5','6','7','8','9',  &
       'A','B','C','D','E','F','G','H','I','J',  &
       'K','L','M','N','O','P','Q','R','S','T',  &
       'U','V','W','X','Y','Z','a','b','c','d',  &
       'e','f','g','h','i','j','k','l','m','n',  &
       'o','p','q','r','s','t','u','v','w','x',  &
       'y','z','{','|'/)

  do n=0,63
     vc(n) = vcscr(n)
  enddo

end subroutine vfinit

!---------------------------------------------------------

!!$subroutine vctran(iun1,iun2,type)
!!$implicit none
!!$integer :: iun1,iun2
!!$character(len=78) :: line,type*(*)
!!$integer :: n,nlines,nbits,nl
!!$
!!$read(iun1,'(a78)') line
!!$write(iun2,'(a78)') line
!!$
!!$if(type.eq.'C') then
!!$   read(line,'(i8)') n
!!$   nlines=n/78+1
!!$elseif(type.eq.'I'.or.type.eq.'F') then
!!$   read(line,'(2i8)')n,nbits
!!$   nlines=(n*nbits/6-1)/78+1
!!$endif
!!$
!!$do 10 nl=1,nlines
!!$   read(iun1,'(a78)') line
!!$   write(iun2,'(a78)') line
!!$10   continue
!!$
!!$return
!!$end

!--------------------------------------------------------

!!$subroutine vcorec(iunit,a,n)
!!$implicit none
!!$integer :: iunit,n
!!$character(len=1) :: a(*)
!!$integer :: nc,ne,i
!!$
!!$
!!$write(iunit,11)n
!!$11   format(i8)
!!$
!!$do 10 nc=1,n,78
!!$   ne=min(n,nc+78)
!!$   write(iunit,'(78a1)') (a(i),i=nc,ne)
!!$10   continue
!!$
!!$return
!!$end

!--------------------------------------------------------

!!$subroutine vcirec(iunit,a,n)
!!$implicit none
!!$integer :: iunit,n
!!$character(len=1) :: a(*)
!!$integer :: nn,nc,ne,i
!!$
!!$read(iunit,11)nn
!!$11   format(i8)
!!$
!!$if(nn.ne.n) then
!!$   print*,' Character count mismatch on vcirec record '
!!$   print*,' Characters on record - ',nn
!!$   print*,' Characters expected  - ',n
!!$   stop 'vcirec'
!!$endif
!!$
!!$do 10 nc=1,n,78
!!$   ne=min(n,nc+78)
!!$   read(iunit,'(78a1)') (a(i),i=nc,ne)
!!$10   continue
!!$
!!$return
!!$end

!--------------------------------------------------------

subroutine vforec(iunit, a, n, nbits, scr, type)
  implicit none
  ! Arguments:
  integer, intent(IN) :: iunit, n, nbits
  real, intent(IN)    :: a(*)
  real, intent(OUT)   :: scr(*)
  character(len=*), intent(IN) :: type
  ! Local Variables:
  character(len=1) :: vc
  common/vform/vc(0:63)
  real(kind=8) :: bias, fact
  real         :: sbias, sfact, amin,amax

  if (vc(0)/='0') call vfinit()

  !     log scaling assumes range of +/- 10e10

  call cscale(a, n, amin, amax)
  bias  = dble(-amin + 1.e-20)
  call cfact(bias, amax, nbits, fact)
  sbias = sngl(bias)
  sfact = sngl(fact)

  write(iunit, 10) n, nbits, bias, fact
!srf increasing size of format for parameter 'n'
!10 format(2i8,2e20.10)
10 format(i12,i8,2e20.10)

  call vwrt(iunit, a, n, bias, fact, nbits, scr, type)

end subroutine vforec

!--------------------------------------------------------

subroutine vwrt(iunit, a, n, bias, fact, nbits, scr, type)
  implicit none
  ! Arguments:
  integer, intent(IN) :: iunit, n, nbits
  real, intent(IN)    :: a(n)
  real, intent(OUT)   :: scr(n)
  character(len=*), intent(IN) :: type
  real(kind=8), intent(IN) ::  bias, fact
  ! Local variables:
  character(len=1) :: vc
  common/vform/vc(0:63)
  character(len=80) :: line
  character(len=05) :: form
!!$  character :: line*80,form*5
  integer :: i, nvalline, nchs, ic, ii, isval, iii, iscr
  real    :: scfct

  if (type=='LIN') then
     do i=1,n !10
        scr(i) = sngl((dble(a(i)) + bias)*fact)
        !10      continue
     enddo
  elseif (type=='LOG') then
     scfct = 2.**(nbits - 1)
     do i=1,n !!11
        scr(i) = (sign(1., a(i))*(log10(max(1.e-10, abs(a(i)))) + 10.)/ &
             20. + 1.)*scfct
        !11      continue
     enddo
  endif

  nvalline = (78*6)/nbits
  nchs     = nbits/6
  do i=1,n,nvalline !20
     ic = 0
     do ii=i,i+nvalline-1 !30
        if (ii>n) goto 31
        isval = int(scr(ii))
        do iii=1,nchs !40
           iscr        = iand(ishft(isval, -6*(nchs-iii)), 63)
           ic          = ic + 1
           line(ic:ic) = vc(iscr)
           !40         continue
        enddo
        !30         continue
     enddo
31   continue
     write (form, 100)   ic
     write (iunit, form) line(1:ic)
     !20   continue
  enddo

100 format('(a',i2,')')

end subroutine vwrt

!--------------------------------------------------------

subroutine vfirec(iunit, a, n, type)
  use dump, only: &
    dumpMessage
    
  implicit none
  include "constants.f90"
  character(len=*),parameter :: header="***(vfirec)***"
  ! Arguments:
  integer, intent(IN)          :: iunit, n
  real, intent(INOUT)          :: a(n)
  character(len=*), intent(IN) :: type
  ! Local variables:
  character(len=1) :: vc
  common/vform/vc(0:63)
!!$  character :: line*80, cs*1
  character(len=80) :: line
  character(len=01) :: cs
  integer :: ich0, ich9, ichcz, ichca, ichla, ichlz
  integer :: i, nvalline, nchs, ic, ii, isval, iii, ics, nn, nbits, nc
  real    :: bias, fact, facti, scfct

  if (vc(0)/='0') call vfinit()

  ich0  = ichar('0')
  ich9  = ichar('9')
  ichcz = ichar('Z')
  ichlz = ichar('z')
  ichca = ichar('A')
  ichla = ichar('a')

  !read (UNIT=iunit, FMT='(2i8)',     ADVANCE='NO') nn, nbits
  !read (UNIT=iunit, FMT='(2e20.10)', ADVANCE='NO') bias, fact
  read (iunit, *) nn, nbits, bias, fact

!!$  print *, "DEBUG-ALF:nn, nbits, bias, fact=", nn, nbits, bias, fact
!!$  call flush(6)

  !!10   format(2i8,2e20.10)
  if (nn/=n) then
     print*,' Words on record - ',nn
     print*,' Words expected  - ',n
     iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion,c_fatal, &
       " Word count mismatch on vfirec record")
  endif

  nvalline = (78*6)/nbits
  nchs     = nbits/6
  do i=1,n,nvalline !20
     read (iunit, '(a78)') line
     ic = 0
     do ii=i,i+nvalline-1 !30
        isval = 0
        if (ii>n) goto 20
        do iii=1,nchs !40
           ic  = ic + 1
           cs  = line(ic:ic)
           ics = ichar(cs)
           if (ics<=ich9) then
              nc = ics - ich0
           elseif (ics<=ichcz) then
              nc = ics - ichca + 10
           else
              nc = ics -ichla  + 36
           endif
           isval = ior(ishft(nc, 6*(nchs-iii)), isval)
           !!40         continue
        enddo
        a(ii) = isval
        !!30      continue
     enddo
20   continue
  enddo
  IF(fact==0.0) THEN
                                !0        1         2         3       0        1         2         3
                                !1234567890123456789012345678901234   1234567890123456789012345678901234
    write(*,*) 'Vfirec: fact is zero. Making facti 1.0E38.','Please, check input files         '
    facti=0.0
  ELSE
    facti = 1./fact
  END IF
  if (type=='LIN') then
     do i=1,n !48
        a(i) = a(i)*facti-bias
        !!48      continue
     enddo
  elseif (type=='LOG') then
     scfct = 2.**(nbits-1)
     do i=1,n !55
        a(i) = sign(1., (a(i)-scfct))*(10.**(abs(20.*(a(i)/scfct-1.)) - 10.))
        !!55      continue
     enddo
  endif

end subroutine vfirec

!--------------------------------------------------------

subroutine cscale(a, n, amin, amax)
  implicit none
  ! Arguments:
  integer, intent(IN) :: n
  real, intent(IN)    :: a(n)
  real, intent(OUT)   :: amin, amax
  ! Local Variables:
  integer :: nn

!!$  amin =  1.e30
!!$  amax = -1.e30
!!$  do nn=1,n !10
!!$     amin = min(amin, a(nn))
!!$     amax = max(amax, a(nn))
!!$     !!10   continue
!!$  enddo

  amin = minval(a(1:n))
  amax = maxval(a(1:n))

  amin = min( 1.e30, amin)
  amax = max(-1.e30, amax)

end subroutine cscale

!---------------------------------------------------------

subroutine cfact(bias, amax, nbits, fact)
  implicit none
  ! Arguments:
  !!double precision, intent(IN)  :: bias
  real(kind=8), intent(IN)      :: bias
  real, intent(IN)              :: amax
  integer, intent(IN)           :: nbits
  !!double precision, intent(OUT) :: fact
  real(kind=8), intent(OUT)     :: fact
  ! Local Variables:
  !double precision :: bignum, tnum
  real(kind=8) :: bignum, tnum

  bignum = dble(2**nbits - 1)
  tnum   = bias + dble(amax)
!-srf - avoid div by zero
  tnum = max(tiny(bignum),tnum)
!-srf
  fact   = bignum/(tnum + 1.d-20)

end subroutine cfact

!--------------------------------------------------------

!!$subroutine viorec(iunit,ia,n,nbits,scr)
!!$implicit none
!!$integer :: iunit,n,nbits
!!$integer :: ia(*),scr(*)
!!$integer :: iamin,iamax
!!$real :: bias,fact
!!$
!!$character(len=1) :: vc
!!$common/vform/vc(0:63)
!!$
!!$if(vc(0).ne.'0') call vfinit
!!$
!!$call cscalei(ia,n,iamin,iamax)
!!$bias=-iamin
!!$fact=1.
!!$if((iamax+bias).gt.(2**nbits-1)) then
!!$  print*,'!! Warning from viorec !! - truncation will occur !!'
!!$  print*,'   Maximum- ',iamax,'  bias- ',bias,'  nbits-',nbits
!!$endif
!!$
!!$write(iunit,10)n,nbits,bias,fact
!!$10   format(2i8,2e20.10)
!!$call vwrti(iunit,ia,n,bias,fact,nbits,scr)
!!$
!!$return
!!$end

!--------------------------------------------------------

!!$subroutine vwrti(iunit,ia,n,bias,fact,nbits,iscr)
!!$implicit none
!!$integer :: iunit,n,nbits
!!$integer :: ia(n),iscr(n)
!!$real :: bias,fact
!!$integer :: i,nvalline,nchs,ic,ii,isval,iii,iiscr
!!$
!!$character(len=1) :: vc
!!$common/vform/vc(0:63)
!!$character line*80,form*5
!!$
!!$do 10 i=1,n
!!$   iscr(i)=(ia(i)+bias)*fact+.001
!!$10   continue
!!$
!!$nvalline=(78*6)/nbits
!!$nchs=nbits/6
!!$do 20 i=1,n,nvalline
!!$   ic=0
!!$   do 30 ii=i,i+nvalline-1
!!$      if(ii.gt.n) goto 31
!!$      isval=iscr(ii)
!!$      do 40 iii=1,nchs
!!$         iiscr=iand(ishft(isval,-6*(nchs-iii)),63)
!!$         ic=ic+1
!!$         line(ic:ic)=vc(iiscr)
!!$40         continue
!!$30      continue
!!$31      continue
!!$   write(form,100) ic
!!$   write(iunit,form) line(1:ic)
!!$20   continue
!!$
!!$100  format('(a',i2,')')
!!$
!!$return
!!$end

!--------------------------------------------------------

subroutine viirec(iunit,ia,n)
implicit none
integer :: n,iunit
integer :: ia(n)

character(len=1) :: vc
common/vform/vc(0:63)
character line*80, cs*1
integer :: ich0,ich9,ichcz,ichca,ichla,ichlz

integer :: nn,nbits,nvalline,nchs,i,ic,ii,isval,iii,ics,nc
real :: bias,fact,facti

ich0=ichar('0')
ich9=ichar('9')
ichcz=ichar('Z')
ichlz=ichar('z')
ichca=ichar('A')
ichla=ichar('a')

if(vc(0).ne.'0') call vfinit

read(iunit,10,end=12) nn,nbits,bias,fact
10   format(2i8,2e20.10)
goto 15
12   continue
n=-1
return
15   continue
if(nn.ne.n) then
   print*,' Word count mismatch on viirec record '
   print*,' Words on record - ',nn
   print*,' Words expected  - ',n
   stop 'viirec'
endif

nvalline=(78*6)/nbits
nchs=nbits/6
do 20 i=1,n,nvalline
   read(iunit,'(a78)') line
   ic=0
   do 30 ii=i,i+nvalline-1
      isval=0
      if(ii.gt.n) goto 20
      do 40 iii=1,nchs
         ic=ic+1
         cs=line(ic:ic)
         ics=ichar(cs)
         if(ics.le.ich9) then
            nc=ics-ich0
         elseif(ics.le.ichcz) then
            nc=ics-ichca+10
         else
            nc=ics-ichla+36
         endif
         isval=ior(ishft(nc,6*(nchs-iii)),isval)
40         continue
      ia(ii)=isval
30      continue
20   continue

facti=1./fact
do 48 i=1,n
   ia(i)=ia(i)*facti-bias
48   continue

return
end

!--------------------------------------------------------

!!$subroutine cscalei(ia,n,iamin,iamax)
!!$implicit none
!!$integer :: n,iamin,iamax
!!$integer :: ia(n)
!!$
!!$iamin= minval(ia(1:n))
!!$iamax= maxval(ia(1:n))
!!$
!!$return
!!$end
