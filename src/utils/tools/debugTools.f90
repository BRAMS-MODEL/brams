module  debugData
    integer :: &
    DbgOn, &
    ivalue, &
    jvalue , &
    kvalue
    real :: endTimeDebug
    real :: beginDumpTime
    character(len=256) :: dirDump
    integer, parameter :: fof=80
    character(len=4),parameter :: lev(40)=(/'l 1=','l 2=','l 3=','l 4=','l 5=','l 6=','l 7=','l 8=','l 9=','l10=', &
                                            'l11=','l12=','l13=','l14=','l15=','l16=','l17=','l18=','l19=','l20=', &
                                            'l21=','l22=','l23=','l24=','l25=','l26=','l27=','l28=','l29=','l30=', &
                                            'l31=','l32=','l33=','l34=','l35=','l36=','l37=','l38=','l39=','l40='/)

end module debugData

subroutine readDebug()
  use debugData, only: &
    DbgOn, &
    ivalue, &
    jvalue , &
    kvalue, &
    endTimeDebug, &
    beginDumpTime, &
    dirdump,fof

  open(unit=55,file='./debug.data')
  read(55,fmt='(I1.1)') DbgOn
  read(55,fmt='(I2.2)') ivalue
  read(55,fmt='(I2.2)') jvalue
  read(55,fmt='(I2.2)') kvalue
  read(55,fmt='(F7.1)') endTimeDebug
  read(55,fmt='(F7.1)') beginDumpTime
  read(55,fmt='(A256)') dirdump
  close(55)

  print *,'Degug information:'
  print *,DbgOn, &
  ivalue, &
  jvalue , &
  kvalue, endTimeDebug, &
  beginDumpTime,trim(dirDump)

  open(fof+1,file=trim(dirdump)//'output-1P.out')
  open(fof+2,file=trim(dirdump)//'output-2P.out')

end subroutine readDebug

subroutine finaliza(time)
  use node_mod, only: nmachs,ixb,iyb,mynum
  use debugData, only: &
    DbgOn, &
    ivalue, &
    jvalue , &
    kvalue, &
    endTimeDebug,fof
  include 'mpif.h'
  real,intent(in) :: time
  integer :: ierr
  !call flush()
  if(time==endtimeDebug .and. DbgOn==1) then
    close(fof)
    close(fof+1)
    close(fof+2)
    print *,'+=============================================================+'
    print *,'####### #     # ######          ####### #######  #####  #######'
    print *,'#       ##    # #     #            #    #       #     #    #'
    print *,'#       # #   # #     #            #    #       #          #'
    print *,'#####   #  #  # #     #            #    #####    #####     #'
    print *,'#       #   # # #     #            #    #             #    #'
    print *,'#       #    ## #     #            #    #       #     #    #'
    print *,'####### #     # ######             #    #######  #####     #'
    print *,'+=============================================================+'
    print *,'           for a total of ',nmachs,' cores'
    call flush(6)
    stop 80808080
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_FINALIZE()
  end if
  return
end subroutine finaliza


subroutine dumpVar3d(var,varn,label,dir,beg,maxsize,ic,all,p1,m1)
  use debugData, only: &
    DbgOn, &
    ivalue, &
    jvalue , &
    kvalue
  use node_mod, only :  nmachs,  & ! intent(in)
       ia,                       & ! intent(in)
       iz,                       & ! intent(in)
       ja,                       & ! intent(in)
       jz,                       & ! intent(in)
       mynum,iyb,ixb,mxp,myp,mzp

  integer, intent(in) :: maxsize,ic,all,p1,m1,beg
  character, intent(in) :: dir !'x' or 'y'
  character(len=*), intent(in) :: varn,label
  real,intent(in) :: var(mzp,mxp,myp)

  character(len=30),dimension(:), allocatable :: cpos
  character(len=2) cmxp
  integer :: i,j,k,ix,ms,ini,fim,fp
  ms=maxsize
  if(DbgOn==0) return
  i=ivalue-ixb(mynum,1)+2
  j=jvalue-iyb(mynum,1)+2
  fp=1
  if(beg/=0) fp=beg

  if((dir=='x' .or. dir=='X') .and. i>=2 .and. j>=2 .and. i<=mxp .and. j<=myp) then
    if(maxsize==0) ms=mxp
    write(cmxp,fmt='(I2.2)') ms
    allocate(cpos(fp:fp+ms))
    if(beg/=0) then
      ini=beg-ixb(mynum,1)+2
      fim=mxp
      ms=fim-ini
    endif


    write (80+nmachs,fmt='(2(A20,1X),2(A5,I2.2),1(A5,I5.5),A5,F30.15)') trim(label),trim(varn), &
                                   ',pi=',ivalue,', pj=',jvalue,',cnt=',ic,',var=',var(kvalue,i,j)
    if(all==1) then
      write (80+nmachs,fmt='(A2,1X'//cmxp//'(I30,1X))') 'X:',(ix,ix=fp,fp+ms)
      write (80+nmachs,fmt='(A2,1X'//cmxp//'(F30.15,1X))') 'X:',var(kvalue,ini:fim,j)
    endif
    if(p1==1 .and. m1==0) write (80+nmachs,fmt='(A2,1X,2(F30.15,1X))') '=>',var(kvalue,i,j),var(kvalue,i+1,j)
    if(p1==0 .and. m1==1) write (80+nmachs,fmt='(A2,1X,2(F30.15,1X))') '=>',var(kvalue,i-1,j),var(kvalue,i,j)
    if(p1==1 .and. m1==1) write (80+nmachs,fmt='(A2,1X,3(F30.15,1X))') '=>',var(kvalue,i-1,j),var(kvalue,i,j),var(kvalue,i+1,j)
  elseif(i>=2 .and. j>=2 .and. i<=mxp .and. j<=myp) then
    if(maxsize==0) ms=myp
    write(cmxp,fmt='(I2.2)') ms
    allocate(cpos(fp:fp+ms))
    if(beg/=0) then
      ini=beg-iyb(mynum,1)+2
      fim=myp
      ms=fim-ini
    endif
    write (80+nmachs,fmt='(A,1X,A,1X,8(A2,1X))') trim(label),trim(varn),'mx','my','mz','ivalue','jvalue','i','j','pn'
    write (80+nmachs,fmt='(A,1X,A,1X,8(I2.2,1X))') trim(label),trim(varn),mxp,myp,mzp,ivalue,jvalue,i,j,mynum
    if(all==1) then
      write (80+nmachs,fmt='(A2,1X'//cmxp//'(I30,1X))') 'Y:',(ix,ix=fp,fp+ms)
      write (80+nmachs,fmt='(A2,1X'//cmxp//'(F30.15,1X))') 'Y:',var(kvalue,i,ini:fim)
    endif
    if(p1==0 .and. m1==0) write (80+nmachs,fmt='(A2,1X,1(F30.15,1X))') '=>',var(kvalue,i,j)
    if(p1==1 .and. m1==0) write (80+nmachs,fmt='(A2,1X,2(F30.15,1X))') '=>',var(kvalue,i,j),var(kvalue,i,j+1)
    if(p1==0 .and. m1==1) write (80+nmachs,fmt='(A2,1X,2(F30.15,1X))') '=>',var(kvalue,i,j-1),var(kvalue,i,j)
    if(p1==1 .and. m1==1) write (80+nmachs,fmt='(A2,1X,3(F30.15,1X))') '=>',var(kvalue,i,j-1),var(kvalue,i,j),var(kvalue,i,j+1)
  endif


end subroutine dumpVar3d

subroutine dumpVar2d(var,varn,label,dir,beg,maxsize,ic,all,p1,m1)
  use debugData, only: &
    DbgOn, &
    ivalue, &
    jvalue , &
    kvalue
  use node_mod, only :  nmachs,  & ! intent(in)
       ia,                       & ! intent(in)
       iz,                       & ! intent(in)
       ja,                       & ! intent(in)
       jz,                       & ! intent(in)
       mynum,iyb,ixb,mxp,myp,mzp

  integer, intent(in) :: maxsize,ic,all,p1,m1,beg
  character, intent(in) :: dir !'x' or 'y'
  character(len=*), intent(in) :: varn,label
  real,intent(in) :: var(mxp,myp)

  character(len=30),dimension(:), allocatable :: cpos
  character(len=2) cmxp
  integer :: i,j,k,ix,ms,ini,fim,fp
  ms=maxsize
  if(DbgOn==0) return
  i=ivalue-ixb(mynum,1)+2
  j=jvalue-iyb(mynum,1)+2
  fp=1
  if(beg/=0) fp=beg

  if((dir=='x' .or. dir=='X') .and. i>=2 .and. j>=2 .and. i<=mxp .and. j<=myp) then
    if(maxsize==0) ms=mxp
    write(cmxp,fmt='(I2.2)') ms
    allocate(cpos(fp:fp+ms))
    if(beg/=0) then
      ini=beg-ixb(mynum,1)+2
      fim=mxp
      ms=fim-ini
    endif
    write (80+nmachs,fmt='(2(A20,1X),2(A5,I2.2),1(A5,I5.5),A5,F30.15)') trim(label),trim(varn), &
                                   ',pi=',ivalue,', pj=',jvalue,',cnt=',ic,',var=',var(i,j)
    if(all==1) then
      write (80+nmachs,fmt='(A2,1X'//cmxp//'(I30,1X))') 'X:',(ix,ix=fp,fp+ms)
      write (80+nmachs,fmt='(A2,1X'//cmxp//'(F30.15,1X))') 'X:',var(ini:fim,j)
    endif
    if(p1==1 .and. m1==0) write (80+nmachs,fmt='(A2,1X,2(F30.15,1X))') '=>',var(i,j),var(i+1,j)
    if(p1==0 .and. m1==1) write (80+nmachs,fmt='(A2,1X,2(F30.15,1X))') '=>',var(i-1,j),var(i,j)
    if(p1==1 .and. m1==1) write (80+nmachs,fmt='(A2,1X,3(F30.15,1X))') '=>',var(i-1,j),var(i,j),var(i+1,j)
  elseif(i>=2 .and. j>=2 .and. i<=mxp .and. j<=myp) then
    if(maxsize==0) ms=myp
    write(cmxp,fmt='(I2.2)') ms
    allocate(cpos(fp:fp+ms))
    if(beg/=0) then
      ini=beg-iyb(mynum,1)+2
      fim=myp
      ms=fim-ini
    endif
    write (80+nmachs,fmt='(A,1X,A,1X,7(A2,1X))') trim(label),trim(varn),'mx','my','ivalue','jvalue','i','j','pn'
    write (80+nmachs,fmt='(A,1X,A,1X,7(I2.2,1X))') trim(label),trim(varn),mxp,myp,ivalue,jvalue,i,j,mynum
    if(all==1) then
      write (80+nmachs,fmt='(A2,1X'//cmxp//'(I30,1X))') 'Y:',(ix,ix=fp,fp+ms)
      write (80+nmachs,fmt='(A2,1X'//cmxp//'(F30.15,1X))') 'Y:',var(i,ini:fim)
    endif
    if(p1==0 .and. m1==0) write (80+nmachs,fmt='(A2,1X,1(F30.15,1X))') '=>',var(i,j)
    if(p1==1 .and. m1==0) write (80+nmachs,fmt='(A2,1X,2(F30.15,1X))') '=>',var(i,j),var(i,j+1)
    if(p1==0 .and. m1==1) write (80+nmachs,fmt='(A2,1X,2(F30.15,1X))') '=>',var(i,j-1),var(i,j)
    if(p1==1 .and. m1==1) write (80+nmachs,fmt='(A2,1X,3(F30.15,1X))') '=>',var(i,j-1),var(i,j),var(i,j+1)
  endif


end subroutine dumpVar2d


subroutine dumpVar1d(var,varn,label,ip,jp)
  use debugData, only: &
    DbgOn, &
    ivalue, &
    jvalue , &
    kvalue
  use node_mod, only :  nmachs,  & ! intent(in)
       ia,                       & ! intent(in)
       iz,                       & ! intent(in)
       ja,                       & ! intent(in)
       jz,                       & ! intent(in)
       mynum,iyb,ixb,mxp,myp,mzp

  integer, intent(in) :: ip,jp
  character(len=*), intent(in) :: varn,label
  real,intent(in) :: var(mzp)
  if(DbgOn==0) return

  i=ivalue-ixb(mynum,1)+2
  j=jvalue-iyb(mynum,1)+2
  fp=1
  if(beg/=0) fp=beg

  if(ip==i .and. jp==j .and. i>=ia .and. j>=ja .and. i<=iz .and. j<=jz) then

    write (80+nmachs,fmt='(2(A20,1X),5(A5,I2.2),A5,F30.15)') trim(label),trim(varn), &
                                   ',pi=',ivalue,', pj=',jvalue,', pk=',kvalue,',ip=',ip,', jp=',jp,',var=',var(kvalue)
  endif


end subroutine dumpVar1d

subroutine dumpVar0d(var,varn,label,ip,jp)
  use debugData, only: &
    DbgOn, &
    ivalue, &
    jvalue , &
    kvalue
  use node_mod, only :  nmachs,  & ! intent(in)
       ia,                       & ! intent(in)
       iz,                       & ! intent(in)
       ja,                       & ! intent(in)
       jz,                       & ! intent(in)
       mynum,iyb,ixb

  integer, intent(in) :: ip,jp
  character(len=*), intent(in) :: varn,label
  real,intent(in) :: var
  if(DbgOn==0) return

  i=ivalue-ixb(mynum,1)+2
  j=jvalue-iyb(mynum,1)+2
  fp=1
  if(beg/=0) fp=beg

  if(ip==i .and. jp==j .and. i>=ia .and. j>=ja .and. i<=iz .and. j<=jz) then

    write (80+nmachs,fmt='(2(A10,1X),2(A5,I2.2),A5,F30.15)') trim(label),trim(varn), &
                                   ',pi=',ivalue,', pj=',jvalue,',var=',var
  endif


end subroutine dumpVar0d


subroutine dumpBar(cnt)
  use debugData, only: &
    DbgOn, &
    ivalue, &
    jvalue , &
    kvalue,beginDumpTime
  use node_mod, only :  nmachs,  & ! intent(in)
       ia,                       & ! intent(in)
       iz,                       & ! intent(in)
       ja,                       & ! intent(in)
       jz,                       & ! intent(in)
       mynum,iyb,ixb,mxp,myp,mzp
  use mem_grid, only: &
                          ngrids,     & ! INTENT(IN)
                          ngrid,      & ! INTENT(IN)
                          npatch,     & ! INTENT(IN)
                          time

  integer, intent(in) :: cnt
  integer :: i,j

  if(DbgOn==0) return
  if(beginDumpTime>time) return

  i=ivalue-ixb(mynum,1)+2
  j=jvalue-iyb(mynum,1)+2

  if(i>=2 .and. j>=2 .and. i<=mxp .and. j<=myp) &
    write (80+nmachs,fmt='(A,F6.1,A,I3.3,A)') "==========",time,"=======",cnt,"========"

end subroutine dumpBar

character(len=3) function n2char(p)
  integer, intent(in) :: p

  write(n2char,fmt='(I3.3)') p

end function n2char

integer function ifl2g(i)
  use node_mod, only :  nmachs,  & ! intent(in)
       ia,                       & ! intent(in)
       iz,                       & ! intent(in)
       ja,                       & ! intent(in)
       jz,                       & ! intent(in)
       mynum,iyb,ixb,mxp,myp,mzp
  integer, intent(in) :: i

  ifl2g=i-ixb(mynum,1)+2
end function ifl2g

integer function jfl2g(i)
  use node_mod, only :  nmachs,  & ! intent(in)
       ia,                       & ! intent(in)
       iz,                       & ! intent(in)
       ja,                       & ! intent(in)
       jz,                       & ! intent(in)
       mynum,iyb,ixb,mxp,myp,mzp
  integer, intent(in) :: i

  jfl2g=i-iyb(mynum,1)+2
end function jfl2g

logical function ispdb(i,j)
  use debugData, only: &
    DbgOn, &
    ivalue, &
    jvalue , &
    kvalue
  use node_mod, only :  nmachs,  & ! intent(in)
         ia,                       & ! intent(in)
         iz,                       & ! intent(in)
         ja,                       & ! intent(in)
         jz,                       & ! intent(in)
         mynum,iyb,ixb,mxp,myp,mzp

  integer, intent(in) :: i,j
  integer :: il,jl

  il=i+ixb(mynum,1)-2
  jl=j+iyb(mynum,1)-2

  if(il==ivalue .and. jl==jvalue) then
    ispdb=.true.
  else
    ispdb=.false.
  end if

end function ispdb

real function cpntarr3d(m1,m2,m3,var,k,i,j)
  integer,intent(in) :: m1,m2,m3,i,j,k
  real, intent(in) :: var(m1,m2,m3)

  cpntarr3d=var(k,i,j)

end function cpntarr3d

subroutine dumpInK3d(var,labelVar,labelPos,ipos,ai,aj)
  use debugData, only: &
    DbgOn, &
    ivalue, &
    jvalue , &
    kvalue,beginDumpTime, lev
  use node_mod, only :  nmachs,  & ! intent(in)
         ia,                       & ! intent(in)
         iz,                       & ! intent(in)
         ja,                       & ! intent(in)
         jz,                       & ! intent(in)
         mynum,iyb,ixb,mxp,myp,mzp
  use mem_grid, only: &
                     ngrids,     & ! INTENT(IN)
                     ngrid,      & ! INTENT(IN)
                     npatch,     & ! INTENT(IN)
                     time

  real, intent(in) :: var(mzp,mxp,myp)
  integer, intent(in) :: ipos,ai,aj
  character(len=*),intent(in) :: labelVar,labelPos
  logical,external :: ispdb
  character(len=2) :: cmzp,cipos
  character(len=8) :: ctim
  character(len=12) :: ci
  integer :: i,j,il,jl,k

  if(beginDumpTime>time) return
  write(cmzp,fmt='(I2.2)') mzp
  write(cipos,fmt='(I2.2)') ipos
  write(ctim,fmt='(F6.1)') time


  do i=ia,iz
    do j=ja,jz
      il=i+ixb(mynum,1)-2
      jl=j+iyb(mynum,1)-2
      if(ispdb(i,j)) then
        write (ci,fmt='(A,I2.2,A,I2.2,A)') '(:,',il+ai,',',jl+aj,')'
        !print *,'('//ctim//') '//trim(labelPos)//' '//trim(labelVar)//ci//' '//cipos//': '; call flush(6)

        write(80+nmachs,fmt='(A,'//cmzp//'(E20.10,1X))') &
        '('//ctim//') '//trim(labelPos)//' '//trim(labelVar)//ci//' '//cipos//': ',var   (:,i+ai,j+aj)
      endif
    enddo
  enddo

end subroutine dumpInK3d

subroutine dumpInK2d(var,labelVar,labelPos,ipos,ai,aj)
  use debugData, only: &
    DbgOn, &
    ivalue, &
    jvalue , &
    kvalue, beginDumpTime
  use node_mod, only :  nmachs,  & ! intent(in)
         ia,                       & ! intent(in)
         iz,                       & ! intent(in)
         ja,                       & ! intent(in)
         jz,                       & ! intent(in)
         mynum,iyb,ixb,mxp,myp,mzp
  use mem_grid, only: &
                     ngrids,     & ! INTENT(IN)
                     ngrid,      & ! INTENT(IN)
                     npatch,     & ! INTENT(IN)
                     time

  real, intent(in) :: var(mxp,myp)
  integer, intent(in) :: ipos,ai,aj
  character(len=*),intent(in) :: labelVar,labelPos
  logical,external :: ispdb
  character(len=2) :: cmzp,cipos
  character(len=8) :: ctim
  character(len=12) :: ci
  integer :: i,j,il,jl

  if(beginDumpTime>time) return
  write(cmzp,fmt='(I2.2)') 1
  write(cipos,fmt='(I2.2)') ipos
  write(ctim,fmt='(F6.1)') time

  do i=ia,iz
    do j=ja,jz
      il=i+ixb(mynum,1)-2
      jl=j+iyb(mynum,1)-2
      write (ci,fmt='(A,I2.2,A,I2.2,A)') '(:,',il+ai,',',jl+aj,')'
      if(ispdb(i,j)) write(80+nmachs,fmt='(A,'//cmzp//'(E20.10,1X))') &
        '('//ctim//') '//trim(labelPos)//' '//trim(labelVar)//ci//' '//cipos//': ',var   (i+ai,j+aj)
    enddo
  enddo

end subroutine dumpInK2d

subroutine dumpInK1d(var,labelVar,labelPos,ipos)
  use debugData, only: &
    DbgOn, &
    ivalue, &
    jvalue , &
    kvalue, beginDumpTime
  use node_mod, only :  nmachs,  & ! intent(in)
         ia,                       & ! intent(in)
         iz,                       & ! intent(in)
         ja,                       & ! intent(in)
         jz,                       & ! intent(in)
         mynum,iyb,ixb,mxp,myp,mzp
  use mem_grid, only: &
              ngrids,     & ! INTENT(IN)
              ngrid,      & ! INTENT(IN)
              npatch,     & ! INTENT(IN)
              time

  real, intent(in) :: var(mzp)
  integer, intent(in) :: ipos
  character(len=*),intent(in) :: labelVar,labelPos
  logical,external :: ispdb
  character(len=2) :: cmzp,cipos
  character(len=8) :: ctim
  integer :: i,j


  if(beginDumpTime>time) return
  write(cmzp,fmt='(I2.2)') mzp
  write(cipos,fmt='(I2.2)') ipos
  write(ctim,fmt='(F6.1)') time

  do i=ia,iz
    do j=ja,jz
      if(ispdb(i,j)) write(80+nmachs,fmt='(A,'//cmzp//'(E20.10,1X))') &
        '('//ctim//') '//trim(labelPos)//' '//trim(labelVar)//' '//cipos//': ',var   (:)
    enddo
  enddo

end subroutine dumpInK1d

subroutine dumpInK0d(var,labelVar,labelPos,ipos)
  use debugData, only: &
    DbgOn, &
    ivalue, &
    jvalue , &
    kvalue, beginDumpTime
  use node_mod, only :  nmachs,  & ! intent(in)
         ia,                       & ! intent(in)
         iz,                       & ! intent(in)
         ja,                       & ! intent(in)
         jz,                       & ! intent(in)
         mynum,iyb,ixb,mxp,myp,mzp
  use mem_grid, only: &
              ngrids,     & ! INTENT(IN)
              ngrid,      & ! INTENT(IN)
              npatch,     & ! INTENT(IN)
              time

  real, intent(in) :: var
  integer, intent(in) :: ipos
  character(len=*),intent(in) :: labelVar,labelPos
  logical,external :: ispdb
  character(len=2) :: cmzp,cipos
  character(len=8) :: ctim
  integer :: i,j


  if(beginDumpTime>time) return
  write(cmzp,fmt='(I2.2)') 1
  write(cipos,fmt='(I2.2)') ipos
  write(ctim,fmt='(F6.1)') time

  do i=ia,iz
    do j=ja,jz
      if(ispdb(i,j)) write(80+nmachs,fmt='(A,'//cmzp//'(E20.10,1X))') &
        '('//ctim//') '//trim(labelPos)//' '//trim(labelVar)//' '//cipos//': ',var
    enddo
  enddo

end subroutine dumpInK0d



subroutine dumpInK1dforce(var,labelVar,labelPos,ipos,i,j)
  use debugData, only: &
    DbgOn, &
    ivalue, &
    jvalue , &
    kvalue, beginDumpTime
  use node_mod, only :  nmachs,  & ! intent(in)
         ia,                       & ! intent(in)
         iz,                       & ! intent(in)
         ja,                       & ! intent(in)
         jz,                       & ! intent(in)
         mynum,iyb,ixb,mxp,myp,mzp
  use mem_grid, only: &
              ngrids,     & ! INTENT(IN)
              ngrid,      & ! INTENT(IN)
              npatch,     & ! INTENT(IN)
              time

  real, intent(in) :: var(mzp)
  integer, intent(in) :: i,j
  character(len=*),intent(in) :: labelVar,labelPos
  logical,external :: ispdb
  character(len=2) :: cmzp,cipos
  character(len=8) :: ctim

  if(beginDumpTime>time) return
  write(cmzp,fmt='(I2.2)') mzp
  write(cipos,fmt='(I2.2)') ipos
  write(ctim,fmt='(F6.1)') time

  if(ispdb(i,j)) write(80+nmachs,fmt='(A,'//cmzp//'(E20.10,1X))') &
        '('//ctim//') '//trim(labelPos)//' '//trim(labelVar)//' '//cipos//': ',var   (:)


end subroutine dumpInK1dforce

subroutine dumpInK0dforce(var,labelVar,labelPos,ipos,i,j)
  use debugData, only: &
    DbgOn, &
    ivalue, &
    jvalue , &
    kvalue, beginDumpTime
  use node_mod, only :  nmachs,  & ! intent(in)
         ia,                       & ! intent(in)
         iz,                       & ! intent(in)
         ja,                       & ! intent(in)
         jz,                       & ! intent(in)
         mynum,iyb,ixb,mxp,myp,mzp
  use mem_grid, only: &
              ngrids,     & ! INTENT(IN)
              ngrid,      & ! INTENT(IN)
              npatch,     & ! INTENT(IN)
              time

  real, intent(in) :: var
  integer, intent(in) :: i,j,ipos
  character(len=*),intent(in) :: labelVar,labelPos
  logical,external :: ispdb
  character(len=2) :: cmzp,cipos
  character(len=8) :: ctim

  if(beginDumpTime>time) return
  write(cmzp,fmt='(I2.2)') 1
  write(cipos,fmt='(I2.2)') ipos
  write(ctim,fmt='(F6.1)') time

  if(ispdb(i,j)) write(80+nmachs,fmt='(A,'//cmzp//'(E20.10,1X))') &
        '('//ctim//') '//trim(labelPos)//' '//trim(labelVar)//' '//cipos//': ',var

  call flush(80+nmachs)
end subroutine dumpInK0dforce

subroutine dumpInK0dLforce(var,labelVar,labelPos,ipos,i,j)
  use debugData, only: &
    DbgOn, &
    ivalue, &
    jvalue , &
    kvalue, beginDumpTime
  use node_mod, only :  nmachs,  & ! intent(in)
         ia,                       & ! intent(in)
         iz,                       & ! intent(in)
         ja,                       & ! intent(in)
         jz,                       & ! intent(in)
         mynum,iyb,ixb,mxp,myp,mzp
  use mem_grid, only: &
              ngrids,     & ! INTENT(IN)
              ngrid,      & ! INTENT(IN)
              npatch,     & ! INTENT(IN)
              time

  logical, intent(in) :: var
  integer, intent(in) :: i,j,ipos
  character(len=*),intent(in) :: labelVar,labelPos
  logical,external :: ispdb
  character(len=2) :: cmzp,cipos
  character(len=8) :: ctim

  if(beginDumpTime>time) return
  write(cmzp,fmt='(I2.2)') 1
  write(cipos,fmt='(I2.2)') ipos
  write(ctim,fmt='(F6.1)') time

  if(ispdb(i,j)) write(80+nmachs,fmt='(A,'//cmzp//'(L1,1X))') &
        '('//ctim//') '//trim(labelPos)//' '//trim(labelVar)//' '//cipos//': ',var

  call flush(80+nmachs)
end subroutine dumpInK0dLforce

subroutine dumpInK0dIforce(var,labelVar,labelPos,ipos,i,j)
  use debugData, only: &
    DbgOn, &
    ivalue, &
    jvalue , &
    kvalue, beginDumpTime
  use node_mod, only :  nmachs,  & ! intent(in)
         ia,                       & ! intent(in)
         iz,                       & ! intent(in)
         ja,                       & ! intent(in)
         jz,                       & ! intent(in)
         mynum,iyb,ixb,mxp,myp,mzp
  use mem_grid, only: &
              ngrids,     & ! INTENT(IN)
              ngrid,      & ! INTENT(IN)
              npatch,     & ! INTENT(IN)
              time

  integer, intent(in) :: var
  integer, intent(in) :: i,j,ipos
  character(len=*),intent(in) :: labelVar,labelPos
  logical,external :: ispdb
  character(len=2) :: cmzp,cipos
  character(len=8) :: ctim

  if(beginDumpTime>time) return
  write(cmzp,fmt='(I2.2)') 1
  write(cipos,fmt='(I2.2)') ipos
  write(ctim,fmt='(F6.1)') time

  if(ispdb(i,j)) write(80+nmachs,fmt='(A,'//cmzp//'(I10.10,1X))') &
        '('//ctim//') '//trim(labelPos)//' '//trim(labelVar)//' '//cipos//': ',var

  call flush(80+nmachs)
end subroutine dumpInK0dIforce


subroutine dumpFullVar3d(var,varn,label,pos)
  use debugData, only: &
    DbgOn, &
    ivalue, &
    jvalue , &
    kvalue, dirDump,fof,beginDumpTime
  use node_mod, only :  nmachs,  & ! intent(in)
         ia,                       & ! intent(in)
         iz,                       & ! intent(in)
         ja,                       & ! intent(in)
         jz,                       & ! intent(in)
         mynum,iyb,ixb,mxp,myp,mzp
  use mem_grid, only: &
              ngrids,     & ! INTENT(IN)
              ngrid,      & ! INTENT(IN)
              npatch,     & ! INTENT(IN)
              time
  character(len=*),intent(in) :: varn,label
  real, intent(in) :: var(mzp,mxp,myp)
  integer, intent(in) :: pos
  integer :: i,j,k,il,jl
  character(len=2) :: cmzp
  character(len=2) :: cipos
  character(len=8) :: ctim
  character :: cnum,cmyn

  if(beginDumpTime>time) return

  write(cmzp,fmt='(I2.2)') mzp
  write(cipos,fmt='(I2.2)') pos
  write(ctim,fmt='(F6.1)') time
  write(cnum,fmt='(I1)') nmachs
  write(cmyn,fmt='(I1)') mynum

  open(unit=90,file=trim(adjustl(dirdump))//trim(adjustl(label))//'-'// &
        trim(adjustl(varn))//'-'//cipos//'-'//trim(adjustl(ctim))//'_'//cnum//cmyn)

    do i=ia,iz
      do j=ja,jz
        il=i+ixb(mynum,1)-2
        jl=j+iyb(mynum,1)-2
        write (90,fmt='(2(I2.2,1X),'//cmzp//'(E20.10,1X))') il,jl,var(:,i,j)
      enddo
    enddo

  close(90)

end subroutine dumpFullVar3d

subroutine dumpFullVar2d(var,varn,label,pos)
  use debugData, only: &
    DbgOn, &
    ivalue, &
    jvalue , &
    kvalue, dirDump,fof,beginDumpTime
  use node_mod, only :  nmachs,  & ! intent(in)
         ia,                       & ! intent(in)
         iz,                       & ! intent(in)
         ja,                       & ! intent(in)
         jz,                       & ! intent(in)
         mynum,iyb,ixb,mxp,myp,mzp
  use mem_grid, only: &
              ngrids,     & ! INTENT(IN)
              ngrid,      & ! INTENT(IN)
              npatch,     & ! INTENT(IN)
              time
  character(len=*),intent(in) :: varn,label
  real, intent(in) :: var(mxp,myp)
  integer, intent(in) :: pos
  integer :: i,j,k,il,jl
  character(len=2) :: cmzp
  character(len=2) :: cipos
  character(len=8) :: ctim
  character :: cnum,cmyn

  if(beginDumpTime>time) return

  write(cmzp,fmt='(I2.2)') 1
  write(cipos,fmt='(I2.2)') pos
  write(ctim,fmt='(F6.1)') time
  write(cnum,fmt='(I1)') nmachs
  write(cmyn,fmt='(I1)') mynum

  open(unit=90,file=trim(adjustl(dirdump))//trim(adjustl(label))//'-'// &
        trim(adjustl(varn))//'-'//cipos//'-'//trim(adjustl(ctim))//'_'//cnum//cmyn)

    do i=ia,iz
      do j=ja,jz
        il=i+ixb(mynum,1)-2
        jl=j+iyb(mynum,1)-2
        write (90,fmt='(2(I2.2,1X),'//cmzp//'(E20.10,1X))') il,jl,var(i,j)
      enddo
    enddo

  close(90)

end subroutine dumpFullVar2d
