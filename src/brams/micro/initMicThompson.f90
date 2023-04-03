module initMicThompson

   use ModDateUtils, only: &
             date_add_to
   
   use genericFunctions
   
   real,allocatable :: dataValues(:,:,:,:,:)
   character(len=256), allocatable :: varName(:)
   integer :: xdef,ydef,zdef,zbeg,zend,tdef,vars
   real :: xbeg,xstep,ybeg,ystep

   real, allocatable :: qnifa(:,:,:,:) !T,Lev,i,j
   real, allocatable :: qnwfa(:,:,:,:) !T,Lev,i,j
   real, allocatable :: qpress(:,:,:,:) !T,Lev,i,j

   integer :: pastMonth

   logical, parameter :: dumpGrads=.true.
   integer, parameter :: funit=33

   contains

   subroutine readDataFriendly()
        !# Read aerosol friendly information for ice and water
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# *Read aerosol friendly information for ice and water
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 07 May 2019 (Tuesday)
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
      !#  &#9744; <br/>
      !# @endtodo
      !#
      !# @warning
      !# Now is under CC-GPL License, please see
      !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
      !# @endwarning
      !#
      !#---
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

      use dump, only: &
           dumpMessage

      use genericFunctions

      use ReadBcst, only: & !Just for broadcast comm infos
          Broadcast

      use parlib, only: &
           parf_bcast ! Subroutine
      
      implicit none
      
      include "constants.f90"
      include "i8.h"
      
      !Parameters (constants)
      character(len=*), parameter :: header="**(readDataFriendly)**"
      character(len=*), parameter :: fileDir='./tables/micro'
      character(len=*), parameter :: fileCtl=fileDir//'/friend_data.ctl'
      integer, parameter :: fNum=33
      
      
      ! Input/Output variables
      
      !Local variables
      character(len=256) :: lixo
      character(len=256) :: dset
      character(len=256) :: title
      real :: undef
      real, allocatable :: vvar(:,:)
      logical :: fileExist
      integer :: recSize
      integer :: rec,nVar,nTime,nLev
      
      
      !Code
      pastMonth=0
      
      if(mchnum==master_num) then

         iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion &
           ,c_notice,'Reading fliendly ice and water aerosols')
         
         inquire(file=trim(filectl), exist=fileExist )
         if ( .not. fileExist ) iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion &
           ,c_fatal,'file '//trim(filectl)//' not found. Please, check it!')
         
         open (unit=fNum, file=trim(filectl), status='old' &
            , access='sequential', form='formatted', action='read' )
            
         read(fNum, *) lixo,dset  !'DSET ^friend_data.dat'
         read(fNum, *) lixo,title !'TITLE QNWFA_QNIFA_SIGMA_MONTHLY_.dat '
         read(fNum, *) lixo,undef !'UNDEF -9.99E33'
         read(fNum, *) lixo,xdef,lixo,xbeg,xstep!'XDEF 288 LINEAR -180 1.25'
         read(fNum, *) lixo,ydef,lixo,ybeg,ystep !'YDEF 181 LINEAR -90 1.0'
         read(fNum, *) lixo,zdef,lixo,zbeg,zend !'ZDEF 30 LINEAR 1 1'
         read(fNum, *) lixo,tdef !'TDEF 12 LINEAR 00z01Jan2007 1mo'
         read(fNum, *) lixo,vars !'VARS 3'
         
         if(.not. allocated(dataValues)) then
            allocate(dataValues(vars,tdef,zdef,xdef,ydef))
            allocate(varName(vars))
         endif
         
         do nVar=1,vars
            read(fNum, *) varName(nVar)
         enddo
         
         close(fNum)
         
         
         dset=trim(fileDir)//'/'//trim(dset(2:))
         inquire(file=trim(dset), exist=fileExist )
         if ( .not. fileExist ) iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion &
              ,c_fatal,'file '//trim(dset)//' not found. Please, check it!')
         
         recSize=xdef*ydef*outRealSize()
         open (fNum,FILE=trim(dSet),FORM='UNFORMATTED' &
                ,STATUS='OLD',ACCESS='DIRECT',RECL=recSize, action='read' )
         
         rec=1
         do nTime=1,tdef
            do nVar=1,vars
               do nLev=1,zdef
                  read (fNum,REC=rec) dataValues(nvar,nTime,nLev,:,:)
                  rec = rec + 1
               enddo
            enddo
         enddo
         
         close(fNum)
         
      endif

      !Send dimensions to all slaves
      call Broadcast(vars, master_num, "vars")
      call Broadcast(tdef, master_num, "tdef")
      call Broadcast(zdef, master_num, "zdef")
      call Broadcast(xdef, master_num, "xdef")
      call Broadcast(ydef, master_num, "ydef")
      call Broadcast(xstep, master_num,"xstep")
      call Broadcast(ystep, master_num,"ystep")


        call interpolateDataFriendly()

   end subroutine readDataFriendly

   subroutine interpolateDataFriendly()
      !# Interpolate the data from firendly aerosol do model's grid
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: Interpolate the data from firendly aerosol do model's grid
      !# use the formula of inverse quadratic distance to interpolate:
      !# ${{\sum_{i=1}^{n}({1\over{d_i^2}}X_i)}}\over{{\sum_{i=1}^{n}({1\over{d_i^2}})}}$
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 07 May 2019 (Tuesday)
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
      !#  &#9744; <br/>
      !# @endtodo
      !#
      !# @warning
      !# Now is under CC-GPL License, please see
      !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
      !# @endwarning
      !#
      !#---
       use node_mod, only: &
           nodei0, & ! intent(in)
           nodej0, & ! intent(in)
           nodemxp, & ! intent(in)
           nodemyp, & ! intent(in)
           nmachs,  & ! intent(in)
           mynum,  &  ! intent(in)
           mchnum, &  ! intent(in)
           master_num, &! intent(in)
           ia,iz,ja,jz, &
           mzp,mxp,myp

        use dump, only: &
           dumpMessage

        use genericFunctions

        use mem_grid, only : &
           grid_g, &
        zt, &
           runtype,  &  ! intent(in)
           iyear1, &    ! intent(in)
           imonth1, &   ! intent(in)
           idate1, &    ! intent(in)
           itime1, &    ! intent(in)
           timmax, &    ! intent(in)
           timeunit, &  ! intent(in)
           nnxp, &      ! intent(in)
           nnyp, &      ! intent(in)
           nnzp,    &   ! intent(in)
           deltaxn, &
           deltayn, &
           oneGlobalGridData, &
           GlobalSizes  ! Subroutine

        use readbcst, only: &
           gatherData   ! Subroutine

        use mem_basic, only: &
            basic_g
          
       use ParLib, only: &
       parf_bcast

       implicit none
   
       include "constants.f90"
       include "i8.h"
   
       !Parameters (constants)
       character(len=*), parameter :: header="**(interpolateDataFriendly)**"
       integer, parameter :: ifm=1
      !# Grid number
      integer, parameter :: nMaxLatLons=10
      !#
      !# Maximum points to distance 
      integer, parameter :: nifa=1
      integer, parameter :: nwfa=2
      integer, parameter :: npre=3
      !# Ids of data inside arrays
   
       !Local variables
      real(kind=kind_rb) :: mLat
      real(kind=kind_rb) :: mLon
      real :: MaxLats
      real :: MaxLons
      real,allocatable :: Glat(:)
      real,allocatable :: Glon(:)
      real :: pressBrams(nnzp(ifm),nnxp(ifm),nnyp(ifm))
      real :: picpi

      integer :: t,i,j,ii,jj,k,kk,recordLen,iRec,funit
      integer :: m2,m3,icount
      integer :: ial,izl,jal,jzl
      real,allocatable :: qnifaGlobal(:,:,:,:)
      real,allocatable :: qnwfaGlobal(:,:,:,:)
      real,allocatable :: qPressGlobal(:,:,:,:)
      real :: deltaLat,limLat
      real :: deltaLon,limLon
      real :: dens,peso,a1,a2
      real :: deltaS
      real :: dist
      real :: piR_180
      integer :: piLon,pfLon,pilat,pflat
      integer :: nMAxLat,nMaxLon
      logical :: lMonth(12)
      integer :: oyr,omn,ody,otm,nz
      character(len=2) :: ct
      integer :: tim
      integer :: ldimx,ldimy,z

      lMonth=.false.
      lMonth(imonth1)=.true.
      call date_add_to(iyear1,imonth1,idate1,itime1 &
               *100,timmax,'s',oyr,omn,ody,otm)
      
      !Preenche o array de meses validos
      if(oyr==iyear1 .and. omn>imonth1) then
         do i=imonth1,omn
            lMonth(i)=.true.
         enddo
      elseif(oyr==iyear1+1) then
         do i=imonth1,12
            lMonth(i)=.true.
         enddo
         do i=1,omn
            lMonth(i)=.true.
         enddo
      elseif(oyr>iyear1+1) then
         lmonth=.true.
      endif
      
      !Allocating the 3 variables qpress, qnifa and qnwfa with model's levels
      !Using local size
      allocate(qnifa(tdef,zdef,mxp,myp)) 
      allocate(qnwfa(tdef,zdef,mxp,myp))
      allocate(qpress(tdef,zdef,mxp,myp))
      
      !Allocating the 3 variables qpress, qnifa and qnwfa with friendly aerosol levels
      !Using the global size
      allocate(qnifaGlobal(tdef,zdef,nnxp(ifm),nnyp(ifm))) 
      allocate(qnwfaGlobal(tdef,zdef,nnxp(ifm),nnyp(ifm)))
      allocate(qPressGlobal(tdef,zdef,nnxp(ifm),nnyp(ifm)))

      if(mchnum==master_num) then
           iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion &
            ,c_notice,'Interpolate fliendly ice and water aerosols')

         !Allocating the list of lons and lats in friendly aerosol
         allocate(Glon(xdef),Glat(ydef))
   
         !Limiting the space around the points to interpolate
         !The code will compute nMaxLatLons points with deltaLat and
         !deltaLon size in each direction 
         deltaLat=oneGlobalGridData(ifm)%global_glat(1,2) &
                 -oneGlobalGridData(ifm)%global_glat(1,1)
         deltaLon=oneGlobalGridData(ifm)%global_glon(2,1) &
                 -oneGlobalGridData(ifm)%global_glon(1,1)
         deltaS=(deltaLat+deltaLon)/2.0
   
         !Fill the list of lats & lons from friendly aerosol
         !Consider equal distances (delta) in all latitudes&long.
         do i=1,xdef
            Glon(i)=xbeg+(i-1)*xstep
         enddo
         do j=1,ydef
            Glat(j)=ybeg+(j-1)*ystep
         enddo
         
         do i=1,xdef
            if(Glon(i)<oneGlobalGridData(ifm)%global_glon(1,1)) pilon=i
            if(Glon(i)>oneGlobalGridData(ifm)%global_glon(nnxp(ifm),1)) then
               pflon=i
               exit
            endif
         enddo
         
         do i=1,ydef
            if(Glat(i)<oneGlobalGridData(ifm)%global_glat(1,1)) pilat=i
            if(Glat(i)>oneGlobalGridData(ifm)%global_glat(1,nnyp(ifm))) then
               pflat=i
               exit
            endif
         enddo   

         nmaxLat=max(int(nnyp(ifm)/real(pflat-pilat))*2,2)
         nmaxLon=max(int(nnxp(ifm)/real(pflon-pilon))*2,2)
         limLat=nMAxLat*deltaLat
         limLon=nMaxLon*deltaLon

         qnifaGlobal=0.0   
         qnwfaGlobal=0.0 
         qPressGlobal=0.0
         do t=1,tdef
            if(.not. lMonth(t)) cycle
            do k=1,zdef 
                  do  i=1,nnxp(ifm)
                     do j=1,nnyp(ifm)
                        dens=0.0
                        do ii=1,xdef
                            do jj=1,ydef
                              if(   Glon(ii)>oneGlobalGridData(ifm)%global_glon(i,j)-limLon &
                              .and. Glon(ii)<oneGlobalGridData(ifm)%global_glon(i,j)+limLon &
                              .and. Glat(jj)>oneGlobalGridData(ifm)%global_glat(i,j)-limLat &
                              .and. Glat(jj)<oneGlobalGridData(ifm)%global_glat(i,j)+limLat) then

                                peso=1.0/arcDistance(oneGlobalGridData(ifm)%global_glat(i,j) &
                                                    ,oneGlobalGridData(ifm)%global_glat(i,j) &
                                                    ,Glat(jj),Glon(ii))**2
    
                                qnifaGlobal(t,k,i,j)=qnifaGlobal(t,k,i,j)+dataValues(nifa,t,k,ii,jj)*peso
                                qnwfaGlobal(t,k,i,j)=qnwfaGlobal(t,k,i,j)+dataValues(nwfa,t,k,ii,jj)*peso
                                qPressGlobal(t,k,i,j)=qPressGlobal(t,k,i,j)+dataValues(npre,t,k,ii,jj)*peso
    
                                dens=dens+peso 
                              endif
                            end do
                        enddo
                        !if (dens.lt.1e-20) dens=1e-20 
                        qnifaGlobal(t,k,i,j)=qnifaGlobal(t,k,i,j)/dens
                        qnwfaGlobal(t,k,i,j)=qnwfaGlobal(t,k,i,j)/dens
                        qPressGlobal(t,k,i,j)=qPressGlobal(t,k,i,j)/dens
                     enddo
                  enddo
            enddo
         enddo


         if(dumpGrads) then
            do t=1,tdef
               if(.not. lMonth(t)) cycle
         
               write(ct,fmt='(I2.2)') t
            iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion &
               ,c_notice,'Writing fliendly ice and water aerosols 1st, after horiz. Interpolation')
            open(unit=funit,file='qnifa_I'//ct//'.ctl' &
               ,action='WRITE',status='replace',form='FORMATTED')
      
            !writing the name of grads file
            write(funit,*) 'dset ^qnifa_I'//ct//'.gra'
            !writing others infos to ctl
            write(funit,*) 'undef -0.9990000E+34'
            write(funit,*) 'title Friendly Data'
            write(funit,*) 'xdef ',nnxp(ifm),' linear ' &
            ,oneGlobalGridData(ifm)%global_glon(1,1),' ',deltaLon
            write(funit,*) 'ydef ',nnyp(ifm),' linear ' &
            ,oneGlobalGridData(ifm)%global_glat(1,1),' ',deltaLat
            write(funit,*) 'zdef ',zdef,'levels ',(z,z=1,zdef)
            write(funit,*) 'tdef 1 linear 00:00z01jan2018     1mo'
            write(funit,*) 'vars 1'
            write(funit,*) 'QNIFA       ',nnzp(ifm),' 99 Ice'
            write(funit,*) 'endvars'
            close(funit)  
      
            recordLen=outRealSize()*nnxp(ifm)*nnyp(ifm)
      
            open(unit=funit,file='qnifa_I'//ct//'.gra',&
            action='WRITE',status='REPLACE',form='UNFORMATTED',access='DIRECT', &
            recl=recordLen)
      
            !# writing grads binary and fill variables
            irec=1
            do k=1,zdef
               write (funit,rec=irec) qnifaGlobal(t,k,:,:)
               irec=irec+1
            enddo
            close(funit)
      
            write(ct,fmt='(I2.2)') t
            iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion &
               ,c_notice,'Writing fliendly ice and water aerosols 1st, after horiz. Interpolation')
               
            open(unit=funit,file='qnwfa_I'//ct//'.ctl' &
               ,action='WRITE',status='replace',form='FORMATTED')
      
            !writing the name of grads file
            write(funit,*) 'dset ^qnwfa_I'//ct//'.gra'
            !writing others infos to ctl
            write(funit,*) 'undef -0.9990000E+34'
            write(funit,*) 'title Friendly Data'
            write(funit,*) 'xdef ',nnxp(ifm),' linear ' &
            ,oneGlobalGridData(ifm)%global_glon(1,1),' ',deltaLon
            write(funit,*) 'ydef ',nnyp(ifm),' linear ' &
            ,oneGlobalGridData(ifm)%global_glat(1,1),' ',deltaLat
            write(funit,*) 'zdef ',zdef,'levels ',(z,z=1,zdef)
            write(funit,*) 'tdef 1 linear 00:00z01jan2018     1mo'
            write(funit,*) 'vars 1'
            write(funit,*) 'QNWFA       ',nnzp(ifm),' 99 Ice'
            write(funit,*) 'endvars'
            close(funit)  
      
            recordLen=4*nnxp(ifm)*nnyp(ifm)
      
            open(unit=funit,file='qnwfa_I'//ct//'.gra',&
            action='WRITE',status='REPLACE',form='UNFORMATTED',access='DIRECT', &
            recl=recordLen)
      
            !# writing grads binary and fill variables
            irec=1
            do k=1,zdef
               write (funit,rec=irec) qnwfaGlobal(t,k,:,:)
               irec=irec+1
            enddo
            close(funit)
      
            open(unit=funit,file='qPressGlobal_I'//ct//'.ctl' &
               ,action='WRITE',status='replace',form='FORMATTED')
      
            !writing the name of grads file
            write(funit,*) 'dset ^qPressGlobal_I'//ct//'.gra'
            !writing others infos to ctl
            write(funit,*) 'undef -0.9990000E+34'
            write(funit,*) 'title Friendly Data'
            write(funit,*) 'xdef ',nnxp(ifm),' linear ' &
            ,oneGlobalGridData(ifm)%global_glon(1,1),' ',deltaLat
            write(funit,*) 'ydef ',nnyp(ifm),' linear ' &
            ,oneGlobalGridData(ifm)%global_glat(1,1),' ',deltaLon
            write(funit,*) 'zdef ',zdef,'levels ',(z,z=1,zdef)
            write(funit,*) 'tdef 1 linear 00:00z01jan2018     1mo'
            write(funit,*) 'vars 1'
            write(funit,*) 'p_wif       ',nnzp(ifm),' 99 Press [Pa]'
            write(funit,*) 'endvars'
            close(funit)  
      
            recordLen=4*nnxp(ifm)*nnyp(ifm)
      
            open(unit=funit,file='qPressGlobal_I'//ct//'.gra',&
            action='WRITE',status='REPLACE',form='UNFORMATTED',access='DIRECT', &
            recl=recordLen)
      
            !# writing grads binary and fill variables
            irec=1
            do k=1,zdef
               write (funit,rec=irec) qPressGlobal(t,k,:,:)
               irec=irec+1
            enddo
            close(funit)
      
            
            enddo
      
         endif
      endif
      
       !Enviando as 3 variaveis para todos os processadores
      do t=1,tdef
         call parf_bcast(qnifaGlobal(t,:,:,:),int(zdef,i8),int(nnxp(ifm),i8),int(nnyp(ifm),i8), master_num)
         call parf_bcast(qnwfaGlobal(t,:,:,:),int(zdef,i8),int(nnxp(ifm),i8),int(nnyp(ifm),i8), master_num)
         call parf_bcast(qPressGlobal(t,:,:,:),int(zdef,i8),int(nnxp(ifm),i8),int(nnyp(ifm),i8), master_num)

         ial = nodei0(mynum,1)+1
         izl = nodei0(mynum,1)+nodemxp(mynum,1)
         jal = nodej0(mynum,1)+1
         jzl = nodej0(mynum,1)+nodemyp(mynum,1)
         
         !REtirando o pedaco que cabe a cada 1
         do k=1,zdef
            qnifa(t,k,1:mxp,1:myp)=qnifaGlobal(t,k,ial:izl,jal:jzl)
            qnwfa(t,k,1:mxp,1:myp)=qnwfaGlobal(t,k,ial:izl,jal:jzl)
            qPress(t,k,1:mxp,1:myp)=qPressGlobal(t,k,ial:izl,jal:jzl)
         enddo

      enddo 
      
       
   end subroutine interpolateDataFriendly

   subroutine adJustFriendlyForMonth(time)
      !# Adjust aerosol data friendly for month and pressures
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: Interpolate the data from firendly aerosol do model's grid
      !# use the formula of inverse quadratic pressure difs to interpolate:
      !# {{\sum_{i=1}^{n}({1\over{deltaP_i^2}}X_i)}}\over{{\sum_{i=1}^{n}({1\over{deltaP_i^2}})}}
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 07 May 2019 (Tuesday)
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
      !#  &#9744; <br/>
      !# @endtodo
      !#
      !# @warning
      !# Now is under CC-GPL License, please see
      !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
      !# @endwarning
      !#
      use mem_micro, only: &
            micro_g

          use mem_grid, only: &
           idate1, &
           imonth1, &
           iyear1, &
           ihour1, &
           itime1, &
           grid_g

        use node_mod, only: &
           nodei0, & ! intent(in)
           nodej0, & ! intent(in)
           nodemxp, & ! intent(in)
           nodemyp, & ! intent(in)
           nmachs,  & ! intent(in)
           mynum,  &  ! intent(in)
           mchnum, &  ! intent(in)
           master_num, &! intent(in)
           ia,iz,ja,jz, &
           mzp,mxp,myp


          use mem_grid, only : &
           grid_g, &
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

       implicit none
   
       include "constants.f90"
   
       !Parameters (constants)
      integer, parameter :: ifm=1
      

       ! Input/Output variables
       real   ,intent(in)    :: time
       !# current time
   
       !Local variables
      integer :: oyr,omn,ody,otm
      integer :: i,j,k,ldimx,ldimy,kk,z
      real :: dens,dist,peso
      real :: pressBrams(mzp,mxp,myp)
      integer :: irec,recordLen
      real :: deltaLat,deltaLon
      character :: cmn

      !Code
      call date_add_to(iyear1,imonth1,idate1,itime1 &
              *100,time,'s',oyr,omn,ody,otm)
            
      if(omn==pastMonth) return

      deltaLat=grid_g(1)%glat(1,2)-grid_g(1)%glat(1,1)
      deltaLon=grid_g(1)%glon(2,1)-grid_g(1)%glon(1,1)
      
      write (cmn,fmt='(I1)') mynum

      pastMonth=omn
      ldimx=nodemxp(mynum,ifm)
      ldimy=nodemyp(mynum,ifm)
      
      micro_g(ifm)%cccnp=0.0
      micro_g(ifm)%cifnp=0.0
      
      do j=1,myp
         do i=1,mxp
            !Calc the model pressure in each layer
            do k=1,mzp
               dens=0.0
               pressBrams(k,i,j)=((basic_g(ifm)%pp(k,i,j) &
                   +basic_g(ifm)%pi0(k,i,j))/c_cp)**c_cpor*C_p00 !Pa
               do kk=k-4,k+4
                  if(kk<1 .or. kk>zdef) cycle
                  dist=pressBrams(k,i,j)-qPress(omn,kk,i,j)
                  peso=1.0/dist**2
                  micro_g(ifm)%cccnp(k,i,j)=micro_g(ifm)%cccnp(k,i,j)+qnwfa(omn,kk,i,j)*peso
                  micro_g(ifm)%cifnp(k,i,j)=micro_g(ifm)%cifnp(k,i,j)+qnifa(omn,kk,i,j)*peso
                  dens=dens+peso
               enddo
               micro_g(ifm)%cccnp(k,i,j)=micro_g(ifm)%cccnp(k,i,j)/dens
               micro_g(ifm)%cifnp(k,i,j)=micro_g(ifm)%cifnp(k,i,j)/dens
            enddo
         enddo
      enddo

      if(dumpGrads) then
        open(unit=funit,file='microg'//cmn//'.ctl' &
                 ,action='WRITE',status='replace',form='FORMATTED')
      
        !writing the name of grads file
        write(funit,*) 'dset ^microg'//cmn//'.gra'
        !writing others infos to ctl
        write(funit,*) 'undef -0.9990000E+34'
        write(funit,*) 'title Microg Data'
        write(funit,*) 'xdef ',mxp,' linear ' &
        ,grid_g(1)%glon(1,1),' ',deltaLon
        write(funit,*) 'ydef ',myp,' linear ' &
        ,grid_g(1)%glat(1,1),' ',deltaLat
        write(funit,*) 'zdef ',mzp,'levels ',(z,z=1,mzp)
        write(funit,*) 'tdef 1 linear 00:00z01jan2018     1mo'
        write(funit,*) 'vars 2'
        write(funit,*) 'cccnp       ',mzp,' 99 cccnp'
        write(funit,*) 'cifnp       ',mzp,' 99 cifnp'
        write(funit,*) 'endvars'
        close(funit)  

        recordLen=4*mxp*myp

        open(unit=funit,file='microg'//cmn//'.gra',&
        action='WRITE',status='REPLACE',form='UNFORMATTED',access='DIRECT', &
        recl=recordLen)

        !# writing grads binary and fill variables
        irec=1
        do k=1,mzp
           write (funit,rec=irec) micro_g(ifm)%cccnp(k,:,:)
           irec=irec+1
        enddo
        do k=1,mzp
           write (funit,rec=irec) micro_g(ifm)%cifnp(k,:,:)
           irec=irec+1
        enddo      
        close(funit)
      endif

   end subroutine adJustFriendlyForMonth


end module initMicThompson
