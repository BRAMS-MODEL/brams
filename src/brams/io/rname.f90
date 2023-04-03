!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine NAMEOUT()

  use dump

  use mem_all

  ! For Shallow Cumulus
  use shcu_vars_const, only : nnshcu,  & ! INTENT(IN)
       shcufrq                           ! INTENT(IN)

  ! For Grell Paramet.
  use mem_grell_param, only : CLOSURE_TYPE ! INTENT(IN)

  ! For Soil Moisture Init.
  use memSoilMoisture, only : SOIL_MOIST,  & ! INTENT(IN)
       SOIL_MOIST_FAIL,                      & ! INTENT(IN)
       USDATA_IN,                            & ! INTENT(IN)
       USMODEL_IN                              ! INTENT(IN)

  ! CATT
  use catt_start, only: CATT           ! intent(in)
  use emission_source_map, only: FIREMAPFN, & ! intent(in)
       plumerise                              ! intent(in)

  use mem_globrad, only: RADDATFN !, &   !  intent(in)
!!$                         RDATFNUM    !  intent(in)

  use plume_utils, only: prfrq !INTENT(IN)

  ! Explicit domain decomposition
  use domain_decomp, only: domain_fname

  ! TEB_SPM
  use teb_spm_start, only: TEB_SPM ! intent(in)

!--(DMK-CCATT-INI)-----------------------------------------------------
  USE monotonic_adv, ONLY: advmnt
!--(DMK-CCATT-FIM)-----------------------------------------------------

  implicit none

  include "contants.f90"
  integer :: ng,np,k,m

  ! This routine prints out a listing of the values of all variables
  ! in the NAMELISTS

  write(*,fmt='(A)') '----------------------------NAMELIST VARIABLES-------'  &
       ,'------------------------'
  
  iErrNumber=dumpMessage(c_tty,c_yes,'RAMSIN','Namelist',c_notice,'NNXP',nnxp(ng),"I4")
  iErrNumber=dumpMessage(c_tty,c_yes,'RAMSIN','Namelist',c_notice,'NNYP',nnyp(ng),"I4")
  iErrNumber=dumpMessage(c_tty,c_yes,'RAMSIN','Namelist',c_notice,'NNZP',nnzp(ng),"I4")
  iErrNumber=dumpMessage(c_tty,c_yes,'RAMSIN','Namelist',c_notice,'NXTNEST',nxtnest(ng),"I16")
  iErrNumber=dumpMessage(c_tty,c_yes,'RAMSIN','Namelist',c_notice,'NSTRATX',nstratx(ng),"I4")


!   'NNXP=',I4
! 'NNYP=',I4,
! 'NNZP=',I4  
! 'NXTNEST=',I4,
! 'NSTRATX=',I4,999(A1,/,I13,4I16))



!   write(6,101)(' ',NNXP(NG),NNYP(NG),NNZP(NG)  &
!        ,NXTNEST(NG),NSTRATX(NG),NG=1,NGRIDS)

  if (len_trim(domain_fname)>0) then
     write(6, FMT='("Domain File Name =",A80)') domain_fname
  else
     write(6, '("No defined Domain File Name for explicit grid partition.")')
  endif


  write(6,102)(' ',NSTRATY(NG),IDIFFK(NG),NNDTRAT(NG)  &
       ,NNQPARM(NG),NINEST(NG),NG=1,NGRIDS)

  print *, "CLOSURE_TYPE (for Grell Param.): ", CLOSURE_TYPE

  write(6,999) (' ', NNSHCU(NG), NG=1, NGRIDS)  ! For Shallow Cumulus Param.

  write(6,103)(' ',NJNEST(NG),NKNEST(NG),NNSTTOP(NG)  &
       ,NNSTBOT(NG),ITOPTFLG(NG),NG=1,NGRIDS)
  write(6,104)(' ',ISSTFLG(NG),IVEGTFLG(NG),ISOILFLG(NG)  &
       ,NOFILFLG(NG),NDVIFLG(NG),NG=1,NGRIDS)

101 format(A1,'   NNXP=',I4,'       NNYP=',I4,'       NNZP=',I4  &
       ,'    NXTNEST=',I4,'    NSTRATX=',I4,999(A1,/,I13,4I16))
102 format(A1,'NSTRATY=',I4,'     IDIFFK=',I4,'    NNDTRAT=',I4  &
       ,'    NNQPARM=',I4,'     NINEST=',I4,999(A1,/,I13,4I16))

999 format(A1,' NNSHCU=',I4,999(A1,/,I13))  ! For Shallow Cumulus Param.

103 format(A1,' NJNEST=',I4,'     NKNEST=',I4,'    NNSTTOP=',I4  &
       ,'    NNSTBOT=',I4,'   ITOPTFLG=',I4,999(A1,/,I13,4I16))
104 format(A1,'ISSTFLG=',I4,'   IVEGTFLG=',I4,'   ISOILFLG=',I4  &
       ,'   NOFILFLG=',I4,'    NDVIFLG=',I4,999(A1,/,I13,4I16))

  print*, ' '
  write(6,201)NGRIDS,NESTZ1,NESTZ2,INITIAL,IOUTPUT,NUDLAT,if_adap
  write(6,202)INITFLD,IHTRAN,NACOUST,NTOPSMTH,KWRITE
  write(6,203)IUPDSST,IZFLAT,IMPL,ICORFLG,NSLCON,IBND
  write(6,204)JBND,LSFLG,NFPT,IDELTAT,ISWRTYP,ILWRTYP
  write(6,205)LONRAD,IMONTH1,IDATE1,IYEAR1,ITIME1,ISFCL
  write(6,206)NVGCON,NPLT,IPSFLG,ITSFLG,IRTSFLG
  write(6,207)IUSFLG,MKCOLTAB,NZG,NZS,IUPDNDVI
  write(6,208)NPATCH,NVEGPAT,LEVEL,ICLOUD
  write(6,209)IRAIN,IPRIS,ISNOW,IAGGR,IGRAUP
  write(6,210)IHAIL,NADDSC,IADVL,IADVF

201 format('  NGRIDS=',I4,'     NESTZ1=',I4,'     NESTZ2=',I4  &
       ,'    INITIAL=',I4,'    IOUTPUT=',I4,'     NUDLAT=',I4,'     IF_ADAP=',I4)
202 format(' INITFLD=',I4,'     IHTRAN=',I4  &
       ,'    NACOUST=',I4,'   NTOPSMTH=',I4,'     KWRITE=',I4)
203 format(' IUPDSST=',I4,'     IZFLAT=',I4,'       IMPL=',I4  &
       ,'    ICORFLG=',I4,'     NSLCON=',I4,'       IBND=',I4)
204 format('    JBND=',I4,'      LSFLG=',I4,'       NFPT=',I4  &
       ,'    IDELTAT=',I4,'    ISWRTYP=',I4,'    ILWRTYP=',I4)
205 format('  LONRAD=',I4,'    IMONTH1=',I4,'     IDATE1=',I4  &
       ,'     IYEAR1=',I4,'     ITIME1=',I4,'      ISFCL=',I4)
206 format('  NVGCON=',I4,'       NPLT=',I4,'     IPSFLG=',I4  &
       ,'     ITSFLG=',I4,'    IRTSFLG=',I4)
207 format('  IUSFLG=',I4,'   MKCOLTAB=',I4,'        NZG=',I4  &
       ,'        NZS=',I4,'   IUPDNDVI=',I4)
208 format('  NPATCH=',I4,'    NVEGPAT=',I4,'      LEVEL=',I4  &
       ,'     ICLOUD=',I4)
209 format('   IRAIN=',I4,'      IPRIS=',I4,'      ISNOW=',I4  &
       ,'      IAGGR=',I4,'     IGRAUP=',I4)
210 format('   IHAIL=',I4,'     NADDSC=',I4,'      IADVL=',I4  &
       ,'      IADVF=',I4)

  print*, ' '
  write(6,301)(' ',TOPTENH(NG),TOPTWVL(NG),CENTLAT(NG),NG=1,NGRIDS)
  write(6,302)(' ',CENTLON(NG),CSX(NG),CSZ(NG),NG=1,NGRIDS)
  write(6,303)(' ',XKHKM(NG),ZKHKM(NG),AKMIN(NG),NG=1,NGRIDS)
  write(6,304)(' ',GRIDU(NG),GRIDV(NG),NG=1,NGRIDS)

301 format(A1,'TOPTENH=',E12.5,'        TOPTWVL=',E12.5  &
       ,'        CENTLAT=',E12.5,999(A1,/,E21.5,2E28.5))
302 format(A1,'CENTLON=',E12.5,'            CSX=',E12.5  &
       ,'            CSZ=',E12.5,999(A1,/,E21.5,2E28.5))
303 format(A1,'  XKHKM=',E12.5,'          ZKHKM=',E12.5  &
       ,'          AKMIN=',E12.5,999(A1,/,E21.5,2E28.5))
304 format(A1,'  GRIDU=',E12.5,'          GRIDV=',E12.5  &
       ,999(A1,/,E21.5,E28.5))

  print*, ' '

  write(6,401)TIMMAX,TIMSTR,FRQHIS
  write(6,402)FRQANL,FRQPRT,TNUDLAT
  write(6,403)TNUDTOP,TNUDCENT,ZNUDTOP
  write(6,404)DTLONG,DELTAX,DELTAY
  write(6,405)POLELAT,POLELON,DELTAZ
  write(6,406)DZRAT,DZMAX,SSPCT
  write(6,407)CPHAS,DISTIM,RADFRQ
  write(6,408)CONFRQ,WCLDBS

  write(6,409)SHCUFRQ  ! For Shallow Cumulus

  write(6,410)PCTLCON,ZROUGH,ALBEDO
  write(6,411)SEATMP,DTHCON,DRTCON
  write(6,412)CPARM,RPARM,PPARM
  write(6,413)SPARM,APARM,GPARM
  write(6,414)HPARM

401 format('   TIMMAX=',E12.5,'         TIMSTR=',E12.5  &
       ,'         FRQHIS=',E12.5)
402 format('   FRQANL=',E12.5,'         FRQPRT=',E12.5  &
       ,'        TNUDLAT=',E12.5)
403 format('  TNUDTOP=',E12.5,'       TNUDCENT=',E12.5  &
       ,'        ZNUDTOP=',E12.5)
404 format('   DTLONG=',E12.5,'         DELTAX=',E12.5  &
       ,'         DELTAY=',E12.5)
405 format('  POLELAT=',E12.5,'        POLELON=',E12.5  &
       ,'         DELTAZ=',E12.5)
406 format('    DZRAT=',E12.5,'          DZMAX=',E12.5  &
       ,'          SSPCT=',E12.5)
407 format('    CPHAS=',E12.5,'         DISTIM=',E12.5  &
       ,'         RADFRQ=',E12.5)
408 format('   CONFRQ=',E12.5,'         WCLDBS=',E12.5)

409 format('   SHCUFRQ=',E12.5)  ! For Shallow Cumulus

410 format('  PCTLCON=',E12.5,'         ZROUGH=',E12.5  &
       ,'         ALBEDO=',E12.5)
411 format('   SEATMP=',E12.5,'         DTHCON=',E12.5  &
       ,'         DRTCON=',E12.5)
412 format('    CPARM=',E12.5,'          RPARM=',E12.5  &
       ,'          PPARM=',E12.5)
413 format('    SPARM=',E12.5,'          APARM=',E12.5  &
       ,'          GPARM=',E12.5)
414 format('    HPARM=',E12.5)

!!$  print*, ' '
!!$  write(6,501)(' ',IPLFLD(NP),PLFMT(NP),IXSCTN(NP)  &
!!$       ,ISBVAL(NP),NP=1,NPLT)
!!$
!!$501 format(A1,' IPLFLD=',A10,'   PLFMT=',A10  &
!!$       ,'   IXSCTN=',I2,'   ISBVAL=',I4  &
!!$       ,999(A1,/,9X,A10,9X,A10,I12,I14))

  print*, ' '
  write(6,601)(' ',ITOPTFN(M),M=1,NGRIDS)
  write(6,602)(' ',ISSTFN(M),M=1,NGRIDS)
  write(6,603)(' ',IVEGTFN(M),M=1,NGRIDS)
  write(6,604)(' ',ISOILFN(M),M=1,NGRIDS)
  write(6,605)(' ',NDVIFN(M),M=1,NGRIDS)
  write(6,606) ' ',trim(COLTABFN)
  write(6,607) ' ',trim(VARFPFX)
  write(6,608) ' ',trim(SSTFPFX)
  write(6,609) ' ',trim(NDVIFPFX)

601 format(A1,'  ITOPTFN=',A40,999(A1,/,11X,A40))
602 format(A1,'   ISSTFN=',A40,999(A1,/,11X,A40))
603 format(A1,'  IVEGTFN=',A40,999(A1,/,11X,A40))
604 format(A1,'  ISOILFN=',A40,999(A1,/,11X,A40))
605 format(A1,'   NDVIFN=',A40,999(A1,/,11X,A40))
606 format(A1,' COLTABFN=',A)
607 format(A1,'  VARFPFX=',A)
608 format(A1,'  SSTFPFX=',A)
609 format(A1,' NDVIFPFX=',A)

  print*, ' '
  write(6,701)EXPNME
  print*, ' '
  write(6,702)HFILIN
  write(6,703)HFILOUT
  write(6,704)AFILOUT
  print*, ' '
  write(6,705)RUNTYPE,TIMEUNIT

701 format('  EXPNME=',A40)
702 format('  HFILIN=',A40)
703 format(' HFILOUT=',A40)
704 format(' AFILOUT=',A40)
705 format(' RUNTYPE=',A10,'      TIMEUNIT=',A3)

  ! For Soil Moisture Init.
  write(6,*) 'SOIL_MOIST = ', SOIL_MOIST
  write(6,*) 'SOIL_MOIST_FAIL = ', SOIL_MOIST_FAIL
  write(6,711) USDATA_IN
  write(6,712) USMODEL_IN
711 format(' USDATA_IN =',A40)
712 format(' USMODEL_IN=',A40)

  print*, ' '
  write(6,901)(' ',SLMSTR(K),SLZ(K),STGOFF(K),K=1,NZG)

901 format(A1,' SLMSTR=',F6.2,'      SLZ=',F7.3,'      STGOFF=',F8.2  &
       ,999(A1,/,F15.2,F17.3,F21.2))

  print*, ' '
  write(6,1001)(ZZ(K),K=1,NNZP(1))

1001 format('ZZ=',8F9.1,/,(F12.1,7F9.1))

  write(6,1101)(NSTRATZ1(K),K=1,NNZP(1))
1101 format(/,'NSTRATZ1=',(t9,23i3) )

  write(6,1201)(NSTRATZ2(K),K=1,NNZP(1))
1201 format(/,'NSTRATZ2=',(t9,23i3) )

  write(6,1301)(gnu(k),k=1,7)
1301 format(/,'GNU!!=',(t9,7f5.2))

  !Extras arrays in CATT
  if (CATT == 1) then
     write(6,FMT='(A,I2)') 'CATT Activated: CATT=', CATT
     write(6,FMT='("FIREMAPFN=",A80)')  trim(FIREMAPFN)
     !!write(6,FMT='("TRACERSFN=",A80)')  trim(TRACERSFN)
     !write(6,FMT='(A,A)')  'FIREMAPFN= ', FIREMAPFN
     !write(6,FMT='(A,A)')  'TRACERSFN= ', TRACERSFN
     !write(6,FMT='(A,I2)') 'NA_EXTRA2D= ', NA_EXTRA2D
     !write(6,FMT='(A,I2)') 'NA_EXTRA3D= ', NA_EXTRA3D
     write(6,FMT='(A,I2)') 'PLUMERISE= ', PLUMERISE
     write(6,FMT='(A,F8.2)') 'PRFRQ= ', PRFRQ
  else
     write(6,FMT='(A,I2)') 'CATT Not activated: CATT=', CATT
  endif

  !New rad carma
  if (ISWRTYP==4 .or. ILWRTYP==4) then
     !write(6,FMT='("RADDATFN  =",A80)')  trim(RADDATFN)
     write(6,FMT='("RADDATFN = ",A80)')  trim(RADDATFN)
     !write(6,FMT='("RDATFNUM=",I3)')  RDATFNUM
  endif

  if (TEB_SPM/=1) then
     write(6,FMT='(A,I2)') 'TEB_SPM Deactivated: TEB_SPM=', TEB_SPM
  else
     write(6,FMT='(A,I2)') 'TEB_SPM Activated: TEB_SPM=', TEB_SPM
  endif
  write(6,"(/,'----------------------------NAMELIST VARIABLES ENDS--',&
       &'------------------------',/)")
end subroutine NAMEOUT
