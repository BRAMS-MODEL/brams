!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

!!$integer function lastchar(str)
!!$implicit none
!!$character(len=*) :: str
!!$integer :: n,ln
!!$! returns last non-blank character position from a string
!!$
!!$ln=len(str)
!!$do n=ln,1,-1
!!$   if(str(n:n).ne.' ') then
!!$      lastchar=n
!!$      return
!!$   endif
!!$enddo
!!$lastchar=0
!!$
!!$return
!!$end
!!$
!!$!***************************************************************************
!!$
!!$integer function ifirstchar(str)
!!$implicit none
!!$character(len=*) :: str
!!$integer :: n,ln
!!$
!!$! returns first non-blank character position from a string
!!$
!!$ln=len(str)
!!$do n=1,ln
!!$   if(str(n:n).ne.' ') then
!!$      ifirstchar=n
!!$      return
!!$   endif
!!$enddo
!!$ifirstchar=1
!!$
!!$return
!!$end
!!$
!***************************************************************************

subroutine deblank(str1, str2, nch)
  implicit none
  character(len=*), intent(IN)  :: str1
  character(len=*), intent(OUT) :: str2
  integer, intent(OUT)          :: nch
  integer :: n, ln

  ! strips blanks from a string and returns number of chars
  
  str2 = ' '
  ln   = len(str1)
  nch  = 0
  do n=1,ln
     if (str1(n:n)/=' ') then
        nch = nch + 1
        str2(nch:nch) = str1(n:n)
     endif
  enddo
end subroutine deblank

!!$!***************************************************************************
!!$
!!$subroutine detab(str1,str2,nch)
!!$implicit none
!!$character(len=*) :: str1,str2
!!$integer :: n,ln,nch
!!$character(len=1) ::  tab
!!$
!!$tab=achar( 9)
!!$
!!$! strips tabs from a string and returns number of chars
!!$
!!$str2=' '
!!$ln=len_trim(str1)
!!$nch=0
!!$do n=1,ln
!!$   if(str1(n:n).ne.tab) then
!!$      !print*,'no tab:',str1(n:n)
!!$      nch=nch+1
!!$      str2(nch:nch)=str1(n:n)
!!$   else
!!$      print*,'found one:',str1
!!$      str2(nch+1:nch+6)='      '
!!$      nch=nch+6
!!$   endif
!!$enddo
!!$
!!$return
!!$end
!!$
!!$!***************************************************************************
!!$
!!$integer function lastslash(str)
!!$implicit none
!!$character(len=*) :: str
!!$integer :: n,ln
!!$
!!$! returns last slash character position from a string
!!$
!!$ln=len(str)
!!$do n=ln,1,-1
!!$   if(str(n:n).eq.'/') then
!!$      lastslash=n
!!$      return
!!$   endif
!!$enddo
!!$lastslash=0
!!$
!!$return
!!$end
!!$
!***************************************************************************

subroutine char_strip_var(line, var, line2)

  implicit none
  ! Arguments
  character(len=*), intent(IN)  :: line
  character(len=*), intent(OUT) :: var, line2
  ! Local variables:
  integer :: nn, ncl, nb

  ! removes instances of a substring from a string

  ncl = len(line)
  do nn=1,ncl
     if (line(nn:nn)/=' ') then
        nb  = index(line(nn:),' ')
        var = line(nn:nn+nb-1)
        exit !goto 25
     endif
  enddo
  !!25 continue
  line2 = line(nn+nb-1:)

end subroutine char_strip_var

!!$!***************************************************************************
!!$
!!$subroutine findln(text,ltext,order)
!!$implicit none
!!$character(len=*) :: text
!!$integer :: ltext,order
!!$integer :: i
!!$
!!$! find first non-blank character if order=0, last non-blank if order=1
!!$
!!$if(order.eq.1) then
!!$   do i=len(text),1,-1
!!$      if(text(i:i).ne.' ') then
!!$         ltext=i
!!$         goto 10
!!$      endif
!!$   enddo
!!$   10 continue
!!$else
!!$   do i=1,len(text)
!!$      if(text(i:i).ne.' ') then
!!$         ltext=i
!!$         goto 20
!!$      endif
!!$   enddo
!!$   20 continue
!!$endif
!!$
!!$return
!!$end
!!$
!!$!***************************************************************************

subroutine parse(str, tokens, ntok)
  implicit none
  ! Arguments
  integer, intent(OUT)            :: ntok
  character(len=*), intent(IN)    :: str
  character(len=*), intent(INOUT) ::tokens(*)
  ! Local Variables:
  character(len=1) :: sep
  integer, parameter :: ntokmax = 100
  integer :: n, nc, npt, nch, ntbeg, ntend

  ! this routine "parses" character string str into different pieces
  ! or tokens by looking for  possible token separators (toks
  ! str contains nch characters.  the number of tokens identified is nto
  ! the character string tokens are stored in tokens.
  
  sep  = ' '
  ntok = 0
  npt  = 1
  nch  = len_trim(str)
  nc   = 1
  do ntok=1,ntokmax
     do n=nc,nch
        if (str(n:n)/=sep) then
           ntbeg = n
           exit
!!$           goto 21
        endif
     enddo
!!$21   continue
     do n=ntbeg,nch
        if (str(n:n)==sep) then
           ntend = n-1
           exit
!!$           goto 22
        endif
        if (n==nch) then
           ntend = n
           exit
!!$           goto 22
        endif
     enddo
!!$22   continue
     tokens(ntok) = str(ntbeg:ntend)
     nc           = ntend + 1
     if (nc>=nch) exit !goto 25
  enddo

!!$25 continue

  !do nc=1,nch
  !   if(str(nc:nc).eq.sep.or.nc.eq.nch)then
  !      if(nc-npt.ge.1)then
  !         ntok=ntok+1
  !         tokens(ntok)=str(npt:nc-1)
  !         if(nc.eq.nch.and.str(nc:nc).ne.sep)then
  !            tokens(ntok)=str(npt:nc)
  !            go to 10
  !         endif
  !      endif
  !      ntok=ntok+1
  !      tokens(ntok)=str(nc:nc)
  !      npt=nc+1
  !      go to 10
  !   endif
  !   10 continue
  !enddo
  
end subroutine parse

!!$!***************************************************************************
!!$
!!$subroutine tokenize(str1,tokens,ntok,toksep,nsep)
!!$implicit none
!!$integer :: nsep,ntok
!!$character(len=*) :: str1,tokens(*)
!!$character(len=1) :: toksep(nsep)
!!$
!!$character(len=256) :: str
!!$integer :: npt,nch,nc,ns
!!$
!!$! this routine "parses" character string str into different pieces
!!$! or tokens by looking for  possible token separators (toks
!!$! str contains nch characters.  the number of tokens identified is nto
!!$! the character string tokens are stored in tokens.
!!$
!!$ntok=0
!!$npt=1
!!$call deblank(str1,str,nch)
!!$do nc=1,nch
!!$   do ns=1,nsep
!!$      if(str(nc:nc).eq.toksep(ns).or.nc.eq.nch) then
!!$         if(nc-npt.ge.1)then
!!$            ntok=ntok+1
!!$            tokens(ntok)=str(npt:nc-1)
!!$            if(nc.eq.nch.and.str(nc:nc).ne.toksep(ns)) then
!!$               tokens(ntok)=str(npt:nc)
!!$               goto 10
!!$            endif
!!$         endif
!!$         ntok=ntok+1
!!$         tokens(ntok)=str(nc:nc)
!!$         npt=nc+1
!!$         goto 10
!!$      endif
!!$   enddo
!!$10      continue
!!$enddo
!!$return
!!$end

!***************************************************************************

subroutine tokenize1(str1, tokens, ntok, toksep)
  implicit none
  ! Arguments:
  character(len=*), intent(IN)  :: str1
  character(len=1), intent(IN)  :: toksep
  integer, intent(OUT)          :: ntok
  character(len=*), intent(OUT) :: tokens(10) !,tokens(*)
!!$  integer :: ntok
!!$  character(len=*) :: str1, tokens(10) !,tokens(*)
!!$  character(len=1) :: toksep
  ! Local Variables:
  character(len=256) :: str
  integer            :: nch, ist, npt, nc

  ! this routine "parses" character string str into different pieces
  ! or tokens by looking for  possible token separators (toks
  ! str contains nch characters.  the number of tokens identified is nto
  ! the character string tokens are stored in tokens.
  
  call deblank(str1,str,nch)
  
  ist  = 1
  if (str(1:1)==toksep) ist = 2
  npt  = ist
  ntok = 0
  do nc=ist,nch
     if (str(nc:nc)==toksep .or. nc==nch) then
        if (nc-npt>=1) then
           ntok         = ntok+1
           tokens(ntok) = str(npt:nc-1)
           if(nc==nch .and. str(nc:nc)/=toksep) then
              tokens(ntok)=str(npt:nc)
              exit !goto 10
           endif
           npt = nc+1
        endif
     endif
  enddo
!10 continue

end subroutine tokenize1

!!$!***************************************************************************
!!$
!!$subroutine tokfind(toks,ntok,str,iff)
!!$implicit none
!!$integer :: ntok,iff
!!$character(len=*) :: toks(*),str
!!$
!!$integer :: n
!!$
!!$! looks for a number of tokens (substrings) within a string
!!$
!!$do n=1,ntok
!!$   !print*,'tokfind-',n,toks(n)(1:lastchar(toks(n)))  &
!!$   !      ,'=====',str(1:lastchar(str))
!!$   if(trim(str) == trim(toks(n))) then
!!$      iff=1
!!$      return
!!$   endif
!!$enddo
!!$iff=0
!!$
!!$return
!!$end
!!$
!!$!***************************************************************************
!!$
!!$subroutine rams_intsort(ni,nums,cstr)
!!$implicit none
!!$integer :: nums(*),ni
!!$character(len=*) :: cstr(*)
!!$
!!$character(len=200) :: cscr
!!$integer :: n,mini,nm,nmm,nscr
!!$
!!$! sort an array of character strings by an associated integer field
!!$
!!$do n=1,ni
!!$   mini=1000000
!!$   do nm=n,ni
!!$      if(nums(nm).lt.mini) then
!!$         nmm=nm
!!$         mini=nums(nm)
!!$      endif
!!$   enddo
!!$   nscr=nums(n)
!!$   nums(n)=nums(nmm)
!!$   nums(nmm)=nscr
!!$   cscr=cstr(n)
!!$   cstr(n)=cstr(nmm)
!!$   cstr(nmm)=cscr
!!$enddo
!!$
!!$return
!!$end
!!$
!!$!***************************************************************************
!!$
!!$subroutine rams_fltsort(ni,xnums,cstr)
!!$implicit none
!!$integer :: ni
!!$real :: xnums(*)
!!$character(len=*) :: cstr(*)
!!$
!!$character(len=200) :: cscr
!!$integer :: n,nm,nmm
!!$real :: xmini,xnscr
!!$
!!$! sort an array of character strings by an associated float field
!!$
!!$do n=1,ni
!!$   xmini=1.e30
!!$   do nm=n,ni
!!$      if(xnums(nm).lt.xmini) then
!!$         nmm=nm
!!$         xmini=xnums(nm)
!!$      endif
!!$   enddo
!!$   xnscr=xnums(n)
!!$   xnums(n)=xnums(nmm)
!!$   xnums(nmm)=xnscr
!!$   cscr=cstr(n)
!!$   cstr(n)=cstr(nmm)
!!$   cstr(nmm)=cscr
!!$enddo
!!$
!!$return
!!$end
