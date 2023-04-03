module aerClimMod

   include "constants.f90"

  integer, parameter :: nlat = 46
  !# number of latitude points in aerosol file
  integer, parameter :: nlon = 72
  !# number of longitude points in aerosol file
  integer, parameter :: nlev=12
  !# number of leves in aerosol file
  integer, parameter :: no_months=12
  !# number of months in aerosol file
  integer, parameter :: no_src_types=6
  !# number of aerosol sources in aerosol file
  real, parameter :: latni=-90.0
  !# initial latitude in aerosol file
  real, parameter :: latnf=90.0
  !# final latitude in aerosol file
  real, parameter :: lonni=-180.0
  !# initial longitude in aerosol file
  real, Parameter :: lonnf=180.0
  !# final longitude in aerosol file
  real, parameter :: ilatn=4.0
  !# delta latitude in aerosol file
  real, parameter :: ilonn=5.0
  !# delta longitude in aerosol file
  !# number of each aerosol specie:
  integer, parameter :: OCcam   = 1 
  integer, parameter :: SScam   = 2
  integer, parameter :: DUSTcam = 3
  integer, parameter :: BCcam   = 4
  integer, parameter :: SULFcam = 5
  integer, parameter :: ASHEScam= 6
  character(len=14),dimension(no_src_types) :: specieName=(/&
     'organic carbon' &
    ,'sea salt      ' &
    ,'dust          ' &
    ,'black carbon  ' &
    ,'sulfalte      ' &
    ,'volcanic ashes' &
      /)
  real,parameter,dimension(nlev) :: pressCam=(/ &
    959.0,884.0,787.0,635.0,470.0,337.5,247.5,180.0,125.0,80.0,45.0,20.0/)

  type varT
    real(kind=kind_rb),pointer ::  value(:,:,:) !k,x,y
    real(kind=kind_rb),pointer ::  aer(:,:,:)   !k,x,y
  end type varT
  type(varT) :: aerCam(no_src_types,no_months) !aer,time
  !# type with variables information frm grads file

  type gll
    real(kind=kind_rb) :: lat
    real(kind=kind_rb) :: lon
  end type gll
  type(gll),allocatable :: position(:,:)
  !# grads calculated geographic information

contains

  subroutine gradsRead(filePath,fileName,glat,glon)
    !# Read a grads file and store information
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: reads a grads file and store information in this module
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 01Apr2019
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !#
    !#@endbug
    !#
    !#@todo
    !#  &#9744; Changes to work with another type of grads files<br/>
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !# See the message in todo
    !# @endwarning
    !#
    !#---

  use dump, only: &
    dumpMessage     !subroutine
  
  use parlib, only: &
         parf_bcast ! Subroutine

  use node_mod, only: &
     nodei0, & ! intent(in)
     nodej0, & ! intent(in)
     nodemxp, & ! intent(in)
     nodemyp, & ! intent(in)
     nmachs,  & ! intent(in)
     mynum,  &  ! intent(in)
     mchnum, &  ! intent(in)
     master_num, &! intent(in)
     ia,iz,ja,jz

  use mem_grid, only : &
     runtype,  &  ! intent(in)
     iyear1, &    ! intent(in)
     imonth1, &   ! intent(in)
     idate1, &    ! intent(in)
     itime1, &    ! intent(in)
     nnxp, &      ! intent(in)
     nnyp, &      ! intent(in)
     nnzp,    &   ! intent(in)
     GlobalSizes  ! Subroutine

  use mem_basic, only: &
      basic_g

  use readbcst, only: &
     gatherData   ! Subroutine

  include "constants.f90"
  include "i8.h"
  character(len=*), parameter :: header="**(gradsREad)**"

  !Parameters(constants)
  integer, parameter :: funit=33
  !# unit to manipulate file
  integer, parameter :: ifm=1 
  !# working with Only one grid


  !Input/output variables
  character(len=*),intent(in) :: fileName
  !# Name of grads file
  character(len=*),intent(in) :: filePath
  !# Grads file directory
  real,intent(in) :: glat(:,:)
  real,intent(in) :: glon(:,:)

  !Local variables
  !#

  logical :: fileExist
  integer :: nv
  integer :: time
  integer :: i,j,k,irec,r,recordLen,kk
  integer :: i1, i2, ic, j1, j2, jc
  integer :: qi1,qi2,qj1,qj2, ncount
  real,allocatable :: vvar(:,:)
  real,allocatable :: prlat(:,:)
  real, allocatable :: prlon(:,:)
  real :: dlonr,dlatr
  real(kind=kind_rb) :: lon,lat
  character(len=16) :: varn
  real :: globalGlon(nnxp(ifm), nnyp(ifm))
  real :: globalGlat(nnxp(ifm), nnyp(ifm))
  real :: usdum(nlev)
  real :: picpi
  real :: pressBrams(nnzp(ifm))
  character(len=2) :: cnzp

  !Code

  !allocate(aerCam(no_src_types,no_months))
  do i=1,no_src_types
    do j=1,no_months
      allocate(aercam(i,j)%value(nlev,nlon,nlat)) 
      aercam(i,j)%value=0.0
      allocate(aercam(i,j)%aer(nnzp(ifm),nnxp(ifm),nnyp(ifm)))
      aercam(i,j)%aer=0.0
    enddo
  enddo
  allocate(vvar(nlon,nlat))

  if(mchnum==master_num) then
  
    inquire(file=trim(filePath)//'/'//trim(fileName), exist=fileExist )
    if ( .not. fileExist ) iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion &
            ,c_fatal,'file '//trim(fileName)//' not found. Please, check it!')
  
    recordLen=4*nlon*nlat
  
    open(unit=funit,file=trim(filePath)//'/'//trim(fileName),&
        action='READ',status='OLD',form='UNFORMATTED',access='DIRECT', &
        recl=recordLen)
  
    iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion &
            ,c_Notice,'Reading grads binary file for climatologic Aerosols!')
  
    irec=1
    do time=1,no_months
      do nv=1,no_src_types
        do k=1,nlev
          read (funit,rec=irec) vvar
          aercam(nv,time)%value(k,:,:)=vvar
          irec=irec+1
        enddo
      enddo
    enddo
  endif


  !Send aerosol values to all slaves
  do time=1,no_months
    do nv=1,no_src_types
      do k=1,nlev
        vvar=aercam(nv,time)%value(k,:,:)
        call parf_bcast(vvar,int(nlon,i8), int(nlat,i8) &
           ,master_num)
        aercam(nv,time)%value(k,:,:)=vvar
      enddo
    enddo 
  enddo

  allocate(prlat(nlon,nlat), stat=ierrnumber)
  if (ierrnumber/=0) ierrnumber=dumpmessage(c_tty,c_yes,header,modelversion &
          ,c_fatal,"error allocating prlats")
  allocate(prlon(nlon,nlat), stat=ierrnumber)
  if (ierrnumber/=0) iErrNumber=dumpmessage(c_tty,c_yes,header,modelversion &
          ,c_fatal,"error allocating prlon")

  call apiPrlatlon(nlon, nlat, prlat, prlon, ilatn, ilonn, latni, lonni)

  varn = 'GLON'
  call gatherData(2, varn, ifm, nnxp(ifm), nnyp(ifm), &
       nmachs, mchnum, mynum, master_num,             &
       glon, globalGlon)
  varn = 'GLAT'
  call gatherData(2, varn, ifm, nnxp(ifm), nnyp(ifm), &
       nmachs, mchnum, mynum, master_num,             &
       glat, globalGlat)

  !write(88,fmt='("Grades: ",2(A2,1X),8(A8,1X),6(A8,1X),8(A8,1X))') &
  !    'i','j','globGlat','globGlon','latni','latnf','lonni','lonnf','prlat','prlon' &
  !    ,'i1', 'i2', 'ic', 'j1', 'j2', 'jc','qi1','qi2','qj1','qj2','max(jj)','min(jj)', &
  !            'max(ii)','min(ii)'

  ! loop in global domain of model
  do time=1,no_months
    do nv=1,no_src_types
      if(nv==ASHEScam .or. nv==SScam) cycle !Volcanic ashes and sea salt from other sources
      do j=ja,jz
        do i=ia,iz  
          do k=2,nnzp(ifm)
            picpi = (basic_g(ifm)%pi0(k,i,j) + basic_g(ifm)%pp(k,i,j)) * c_cpi
            pressBrams(k-1) = (c_p00 * picpi ** c_cpor)*c_adjust
          enddo

          !Avoiding points out of domain grid                
          if ( globalGlat(i,j)<latni .or. globalGlat(i,j)>latnf & 
            .or. globalGlon(i,j)<lonni .or. globalGlon(i,j)>lonnf) cycle 

          call interpolacao(globalGlon(i,j), globalGlat(i,j), nlon, nlat, &
                  prlat, prlon, i1, i2, ic, j1, j2, jc)

          if (ic>=0 .and. jc>=0) then
            !print*,ic,jc,i,j,ifm
            dlonr = 0.5*(globalGlon(nnxp(ifm)-1,j) - globalGlon(1,j)) &
                    /float(nnxp(ifm)-1)
            dlatr = 0.5*(globalGlat(i,nnyp(ifm)-1) - globalGlat(i,1)) &
                    /float(nnyp(ifm)-1) 
            qi1   = int(dlonr/ilonn+0.5)
            qi2   = int(dlonr/ilonn+0.5)
            qj1   = int(dlatr/ilatn+0.5)
            qj2   = int(dlatr/ilatn+0.5)

            !write(88,fmt='("....... ",2(I2.2,1X),8(F8.2,1X),6(I8,1X),8(I8,1X))') &
            !  i,j,globalGlat(i,j),globalGlon(i,j),latni,latnf,lonni,lonnf,prlat(i,j),prlon(i,j) &
            !  ,i1, i2, ic, j1, j2, jc,qi1,qi2,qj1,qj2,max(1,jc-qj1),min(nlat,jc+qj2), &
            !  max(1,ic-qi1),min(nlon,ic+qi2)

            do k=1,nlev
              ncount = 0
              usdum(k)=0.
              do jj=max(1,jc-qj1),min(nlat,jc+qj2)
                do ii=max(1,ic-qi1),min(nlon,ic+qi2)
                  if (aercam(nv,time)%value(k,ii,jj)>0.0 &
                      .and. aercam(nv,time)%value(k,ii,jj)<1.0e+38) then
                    ncount = ncount + 1
                    usdum(k) = usdum(k) + aercam(nv,time)%value(k,ii,jj)
                  endif
                enddo
              enddo
              usdum(k) = usdum(k)/(float(ncount))
              !write(*,fmt='("Usdum:",5(I2.2,1X),E18.8)') time,nv,i,j,k,usdum(k)
            enddo

            do k=nnzp(ifm),1,-1
              do kk=nlev,1,-1
                if(pressBrams(k)>=pressCam(kk)) aercam(nv,time)%aer(k,i,j)=usdum(kk)
              enddo
            enddo

          endif
        enddo       
      enddo
    enddo
  enddo

end subroutine  gradsREad

end module aerClimMod

