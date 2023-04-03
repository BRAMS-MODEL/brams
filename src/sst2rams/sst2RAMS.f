
c Este Programa eh disparado pelo script get_sst.sh e escreve os arquivos de sst
c no formato do RAMS.

c Demerval S. Moreira - 11/Nov/2002

      Program escreve_sst
c  read OIv2 SST file and land mask and print selected values.

c  sst   - sea surface temperature array (deg C)
c  err   - normalized error variance
c  ice   - ice concentration array (%)  (0-100,   >100 = land or coast)
c  iyrst - year of start date of analysis
c  imst  - month of start date of analysis
c  idst  - day of start date of analysis
c  iyrnd - year of end date of analysis
c  imnd  - month of end date of analysis
c  idnd  - day of end date of analysis
c  ndays - number of days in analysis (start date thru enddate)
c  index - analysis version for reference
c  xlon  - longitude of center of grid square
c  xlat  - latitude of center of grid square
c  tagls - land/sea tag array (0=land, 1=water)

c  NOTES:  
c   - land values for sst do not necessarily coincide with land values 
c       from ice analysis
c   - recl definition for the direct access open depends on the compiler
c     (number of 4-byte words OR number of bytes).
c     Choose the appropiate parameter statement for krecl
c  

      
      implicit none
      integer krecl
      parameter(krecl=360*180*4)      ! number of bytes
      real*4 sst(360,180),sst2(360,180),err(360,180),xlon,xlat
      integer*4 iyrst,imst,idst,iyrnd,imnd,idnd,ndays,index,i,j
      integer a(11)

      INTEGER IDM,JDM
      PARAMETER (IDM=180,JDM=360)
      real SCR(32761) !SCR(40000)
      character FPREF1*9,arq*80,dia*8,mes*2    !mesc*3,
      
      call getarg(1,arq)
      call getarg(2,dia)
      
      if (dia(1:1).eq."") then
        print*,"Entre com o arq. de sst (ex: /p1/ramsop/oisst.20020911)"
        read(*,'(A)') arq
        print*,"Entre com a data (ex: 20020911)"
        read(*,'(A)') dia
      endif
      !arq='oisst.20050601'
      !dia='20050601'
      

      open(11,file=arq,STATUS='old'
     &       ,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=64811*4)
      
      read(11,rec=1) a,sst2
      call SWAP32(a,11)
      call SWAP32(int(sst2),360*180)
      
      do i=1,180
        do j=1,180  ! 360
          sst(j,i)=sst2(j+180,i)+273.15                  
          sst(j+180,i)=sst2(j,i)+273.15                  
          write(*,'(f8.2,1x,$)') sst(i,j)
        enddo
      enddo 

      iyrst=a(2)
      imst=a(3)
      idst=a(4)
      iyrnd=a(5)
      imnd=a(6)
      idnd=a(7)
      ndays=a(8)
      index=a(9)
      print*,iyrst,imst,idst,iyrnd,imnd,idnd,ndays,index,i,j
           
      if (imst.eq.imnd.or.idnd.lt.4) then
        write(mes,'(i2.2)') imst
      else
        write(mes,'(i2.2)') imnd
      endif   
          
c      if (mes.eq.01) mesc="JAN"     
c      if (mes.eq.02) mesc="FEB"
c      if (mes.eq.03) mesc="MAR"
c      if (mes.eq.04) mesc="APR"
c      if (mes.eq.05) mesc="MAY"
c      if (mes.eq.06) mesc="JUN"
c      if (mes.eq.07) mesc="JUL"
c      if (mes.eq.08) mesc="AUG"
c      if (mes.eq.09) mesc="SEP"
c      if (mes.eq.10) mesc="OCT"
c      if (mes.eq.11) mesc="NOV"
c      if (mes.eq.12) mesc="DEC"
      
      if (mes.ne.dia(5:6)) STOP "mes nao confere, saindo..."
      
      FPREF1="W"//dia
      print*
      print'(2a)', "gerando os arquivos: ",FPREF1
      !write(*,'(360(f6.3,1x))') ((sst(i,j),i=1,360),j=1,180)

      !print*,sst(:,:)

      CALL MKTOPO(JDM,IDM,sst,-180,-90,1.,180,FPREF1,SCR,-0.5,-0.5)
      
      END

C-------------------------------------------------------------------------
C         
      SUBROUTINE MKTOPO(NX,NY,TOPO,IWLON,ISLAT,TRES,IBLSIZE,FPREF,SCR
     +    ,offlat,offlon)
      CHARACTER*9 FPREF
      DIMENSION TOPO(NX,NY),SCR(*)
      CHARACTER*80 TITLE1,TITLE2,TITLE3
C
      IBLDIM=INT(FLOAT(IBLSIZE)/TRES+.001)+1
C
      NSQX=NX/(IBLDIM-1) !+1
      NSQY=NY/(IBLDIM-1) !+1
cccccccccccc
            WRITE(TITLE1,'(A15)') FPREF//'HEADER'
c            TITLE3=FPREF(1:NBL)//TITLE1(1:10)
c            OPEN(29,STATUS='UNKNOWN',FILE=TITLE3,FORM='FORMATTED')
c            PRINT*, 'Making file-',TITLE3
            OPEN(29,STATUS='UNKNOWN',FILE=TITLE1,FORM='FORMATTED')
            PRINT*, 'Making file-',TITLE1
            WRITE(29,14)IBLSIZE,IBLDIM,ISLAT,IWLON,offlat,offlon
 14         FORMAT(4I5,2f10.6)
            CLOSE(29)
cccccccccccc
      DO 1 NSX=1,NSQX
         DO 2 NSY=1,NSQY
            I=(NSX-1)*(IBLDIM-1)+1
            J=(NSY-1)*(IBLDIM-1)+1
            LAT=ISLAT+(NSY-1)*IBLSIZE
            LON=IWLON+(NSX-1)*IBLSIZE

            LATT=ABS(LAT)/10
            LATO=ABS(LAT)-LATT*10
            LONH=ABS(LON)/100
            LONT=(ABS(LON)-LONH*100)/10
            LONO=ABS(LON)-LONH*100-LONT*10

            IF(LAT.GE.0)THEN
               WRITE(TITLE1,'(2I1,A1)')LATT,LATO,'N'
            ELSE
               WRITE(TITLE1,'(2I1,A1)')LATT,LATO,'S'
            ENDIF
            IF(LON.GE.0)THEN
               WRITE(TITLE2,'(3I1,A1)')LONH,LONT,LONO,'E'
            ELSE
               WRITE(TITLE2,'(3I1,A1)')LONH,LONT,LONO,'W'
            ENDIF

            NBL=INDEX(FPREF,' ')-1
            IF(NBL.EQ.-1) NBL=LEN(FPREF)
            TITLE3=FPREF(1:NBL)//TITLE1(1:3)//TITLE2(1:4)
            PRINT*, NSX,NSY,i,j,LAT,LON
            PRINT*, 'Making file-',TITLE3
            
            NPT=1
            DO 10 JJ=J,(IBLDIM-1)+J
               DO 11 II=I,(IBLDIM-1)+I
                  SCR(NPT)=TOPO(II,JJ)
                  !write(*,'(f8.2,1x,$)') SCR(NPT)
!                  if(SCR(NPT)  > 271.35) print*,'sst=',SCR(NPT)-271.35
                  NPT=NPT+1
 11            CONTINUE
 10         CONTINUE

            OPEN(29,STATUS='UNKNOWN',FILE=TITLE3,FORM='FORMATTED')
            CALL VFOREC(29,SCR,IBLDIM*IBLDIM,12,SCR,'LIN')
            CLOSE(29)
 2       CONTINUE
 1    CONTINUE
C
      END
C
      subroutine vfinit
      character*1 vc, vcscr(0:63)
      common/vform/vc(0:63)
      data vcscr/'0','1','2','3','4','5','6','7','8','9'
     +          ,'A','B','C','D','E','F','G','H','I','J'
     +          ,'K','L','M','N','O','P','Q','R','S','T'
     +          ,'U','V','W','X','Y','Z','a','b','c','d'
     +          ,'e','f','g','h','i','j','k','l','m','n'
     +          ,'o','p','q','r','s','t','u','v','w','x'
     +          ,'y','z','{','|'/

      do 10 n=0,63
         vc(n)=vcscr(n)
  10  continue

      return
      end
c---------------------------------------------------------
      subroutine vctran(iun1,iun2,type)
      character*78 line,type*(*)

      read(iun1,'(a78)') line
      write(iun2,'(a78)') line
      
      if(type.eq.'C') then
         read(line,'(i8)') n
         nlines=n/78+1
      elseif(type.eq.'I'.or.type.eq.'F') then
         read(line,'(2i8)')n,nbits
         nlines=(n*nbits/6-1)/78+1
      endif
      
      do 10 nl=1,nlines
         read(iun1,'(a78)') line
         write(iun2,'(a78)') line
 10   continue
      
      return
      end
c--------------------------------------------------------
      subroutine vcorec(iunit,a,n)
      character*1 a(1)

      write(iunit,11)n
 11   format(i8)

      do 10 nc=1,n,78
         ne=min(n,nc+78)
         write(iunit,'(78a1)') (a(i),i=nc,ne)
 10   continue
      
      return
      end
c--------------------------------------------------------
      subroutine vcirec(iunit,a,n)
      character*1 a(1)

      read(iunit,11)nn
 11   format(i8)
      
      if(nn.ne.n) then
         print*,' Character count mismatch on vcirec record '
         print*,' Characters on record - ',nn
         print*,' Characters expected  - ',n
         stop 'vcirec'
      endif

      do 10 nc=1,n,78
         ne=min(n,nc+78)
         read(iunit,'(78a1)') (a(i),i=nc,ne)
 10   continue
      
      return
      end
c--------------------------------------------------------
      subroutine vforec(iunit,a,n,nbits,scr,type)
      double precision  bias, fact
      dimension a(*),scr(*)
      character*(*) type
      character*1 vc
      common/vform/vc(0:63)
c
      if(vc(0).ne.'0') call vfinit
c
c        log scaling assumes range of +/- 10e10
c
      call cscale(a,n,amin,amax)
      bias=dble(-amin+1.e-20)
      call cfact(bias,amax,nbits,fact)
      sbias=sngl(bias)
      sfact=sngl(fact)
c      print*,' amin,amax = ',amin,amax
c      print*,' bias,fact = ',bias,fact
      
      write(iunit,10)n,nbits,bias,fact
 10   format(2i8,2e20.10)
      call vwrt(iunit,a,n,bias,fact,nbits,scr,type)

      return
      end
c--------------------------------------------------------
      subroutine vwrt(iunit,a,n,bias,fact,nbits,scr,type)
      double precision  bias, fact
      character*(*) type
      character*1 vc
      common/vform/vc(0:63)
      dimension a(n),scr(n)
      character line*80,form*5

      if(type.eq.'LIN') then
         do 10 i=1,n
            scr(i)=sngl( ( dble(a(i))+bias ) *fact )
 10      continue
      elseif(type.eq.'LOG') then
         scfct=2.**(nbits-1)
         do 11 i=1,n
            scr(i)=(sign(1.,a(i))*(log10(max(1.e-10,abs(a(i))))+10.)
     +           /20.+1.)*scfct
 11      continue
      endif

      nvalline=(78*6)/nbits
      nchs=nbits/6
      do 20 i=1,n,nvalline
         ic=0
         do 30 ii=i,i+nvalline-1
            if(ii.gt.n) go to 31
            isval=int(scr(ii))
            do 40 iii=1,nchs
               iscr=intand(intrshft(isval,6*(nchs-iii)),63)
               ic=ic+1
               line(ic:ic)=vc(iscr)
 40         continue
 30      continue
 31      continue
         write(form,100) ic
         write(iunit,form) line(1:ic)
 20   continue                 

 100  format('(a',i2,')')

      return
      end
c--------------------------------------------------------
      subroutine vfirec(iunit,a,n,type)
      character*1 vc
      character*(*) type
      common/vform/vc(0:63)
      character line*80, cs*1
      dimension a(*)

      if(vc(0).ne.'0') call vfinit

      ich0=ichar('0')
      ich9=ichar('9')
      ichcz=ichar('Z')
      ichlz=ichar('z')
      ichca=ichar('A')
      ichla=ichar('a')
      
      read(iunit,10)nn,nbits,bias,fact
 10   format(2i8,2e20.10)
      if(nn.ne.n) then
         print*,' Word count mismatch on vfirec record '
         print*,' Words on record - ',nn
         print*,' Words expected  - ',n
         stop 'vfirec'
      endif

      nvalline=(78*6)/nbits
      nchs=nbits/6
      do 20 i=1,n,nvalline
         read(iunit,'(a78)') line
         ic=0
         do 30 ii=i,i+nvalline-1
            isval=0
            if(ii.gt.n) go to 20
            do 40 iii=1,nchs
               ic=ic+1
               cs=line(ic:ic)
               ics=ichar(cs)
               if(ics.le.ich9)then
                  nc=ics-ich0
               elseif(ics.le.ichcz) then
                  nc=ics-ichca+10
               else
                  nc=ics-ichla+36
               endif
               isval=intor(intlshft(nc,6*(nchs-iii)),isval)
 40         continue
            a(ii)=isval
 30      continue
 20   continue

      facti=1./fact
      if(type.eq.'LIN') then
         do 48 i=1,n
            a(i)=a(i)*facti-bias
 48      continue
      elseif(type.eq.'LOG') then
         scfct=2.**(nbits-1)
         do 55 i=1,n
            a(i)=sign(1.,a(i)-scfct)
     +           *(10.**(abs(20.*(a(i)/scfct-1.))-10.))
 55      continue
      endif

      return
      end
c--------------------------------------------------------
      subroutine cscale(a,n,amin,amax)
      dimension a(n)

      amin=1.e30
      amax=-1.e30
      do 10 nn=1,n
         amin=min(amin,a(nn))
         amax=max(amax,a(nn))
 10   continue

      return
      end
c---------------------------------------------------------      
      subroutine cfact(bias,amax,nbits,fact)
      double precision  bias, bignum, fact, tnum

      bignum=dble(2**nbits-1)
      tnum=bias+dble(amax)
      fact=bignum/(tnum+1.d-20)

      return
      end
      
      subroutine viorec(iunit,ia,n,nbits,scr)
      dimension ia(*),scr(*)
      character*1 vc
      common/vform/vc(0:63)
      
      if(vc(0).ne.'0') call vfinit


      call cscalei(ia,n,iamin,iamax)
      bias=-iamin
      fact=1.
      if((iamax+bias).gt.(2**nbits-1)) then
        print*,'!! Warning from viorec !! - truncation will occur !!'
        print*,'   Maximum- ',iamax,'  bias- ',bias,'  nbits-',nbits
      endif

      write(iunit,10)n,nbits,bias,fact
 10   format(2i8,2e20.10)
      call vwrti(iunit,ia,n,bias,fact,nbits,scr)

      return
      end
c--------------------------------------------------------
      subroutine vwrti(iunit,ia,n,bias,fact,nbits,scr)
      character*1 vc
      common/vform/vc(0:63)
      dimension ia(n),scr(n),iscr(n)
      character line*80,form*5

      do 10 i=1,n
         iscr(i)=(ia(i)+bias)*fact+.001
	 scr(i)=iscr(i)
 10   continue

      nvalline=(78*6)/nbits
      nchs=nbits/6
      do 20 i=1,n,nvalline
         ic=0
         do 30 ii=i,i+nvalline-1
            if(ii.gt.n) go to 31
            isval=iscr(ii)
            do 40 iii=1,nchs
               iiscr=intand(intrshft(isval,6*(nchs-iii)),63)
               ic=ic+1
               line(ic:ic)=vc(iiscr)
 40         continue
 30      continue
 31      continue
         write(form,100) ic
         write(iunit,form) line(1:ic)
 20   continue

 100  format('(a',i2,')')

      return
      end
c--------------------------------------------------------
      subroutine viirec(iunit,ia,n)
      character*1 vc
      common/vform/vc(0:63)
      character line*80, cs*1
      dimension ia(*)
      data ich0/48/,ich9/57/,ichcz/90/,ichca/65/,ichla/97/
      ich0=ichar('0')
      ich9=ichar('9')
      ichcz=ichar('Z')
      ichlz=ichar('z')
      ichca=ichar('A')
      ichla=ichar('a')
      
      if(vc(0).ne.'0') call vfinit

      read(iunit,10,end=12)nn,nbits,bias,fact
 10   format(2i8,2e20.10)
      go to 15
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
            if(ii.gt.n) go to 20
            do 40 iii=1,nchs
               ic=ic+1
               cs=line(ic:ic)
               ics=ichar(cs)
               if(ics.le.ich9)then
                  nc=ics-ich0
               elseif(ics.le.ichcz) then
                  nc=ics-ichca+10
               else
                  nc=ics-ichla+36
               endif
               isval=intor(intlshft(nc,6*(nchs-iii)),isval)
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
c--------------------------------------------------------
      subroutine cscalei(ia,n,iamin,iamax)
      dimension ia(n)

      iamin= 1000000000
      iamax=-1000000000
      do 10 nn=1,n
         iamin=min(iamin,ia(nn))
         iamax=max(iamax,ia(nn))
 10   continue

      return
      end
      FUNCTION INTRSHFT(IWORD,NSHFT)
      INTRSHFT=ISHFT(IWORD,-NSHFT)
      RETURN
      END
      FUNCTION INTLSHFT(IWORD,NSHFT)
      INTLSHFT=ISHFT(IWORD,NSHFT)
      RETURN
      END
      FUNCTION INTAND(IWORD1,IWORD2)
      INTAND=IAND(IWORD1,IWORD2)
      RETURN
      END
      FUNCTION INTOR(IWORD1,IWORD2)
      INTOR=IOR(IWORD1,IWORD2)
      RETURN
      END

       
c============================================================================
       SUBROUTINE SWAP32(A,N)
c============================================================================
C
C      REVERSE ORDER OF BYTES IN INTEGER*4 WORD, or REAL*4
C
       INTEGER*4   A(N)
C
       CHARACTER*1 JTEMP(4)
       CHARACTER*1 KTEMP
C
       EQUIVALENCE (JTEMP(1),ITEMP)
C
       SAVE
C
       DO 10 I = 1,N
         ITEMP    = A(I)
         KTEMP    = JTEMP(4)
         JTEMP(4) = JTEMP(1)
         JTEMP(1) = KTEMP
         KTEMP    = JTEMP(3)
         JTEMP(3) = JTEMP(2)
         JTEMP(2) = KTEMP
         A(I)     = ITEMP
 10    CONTINUE
       RETURN
       END
c

