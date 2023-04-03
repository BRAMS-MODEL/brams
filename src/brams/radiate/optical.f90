module optical
   !# Module to setup the optical properties of the model
   !#
   !# @note
   !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
   !#
   !# **Brief**: The module was take from BRAMS' model original and modified Carma radiation
   !#  (Community Aerosol and Radiation Model for Atmosphere) and setup the atmosphere optical properties.
   !# ---
   !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
   !# ---
   !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
   !#             Longo, K.M. **&#9993;**<mailto:karla.longo@inpe.br>
   !#             Rosario, N. E. **&#9993;**<niltoncvbr@gmail.com>
   !# ---
   !# **Date**: 2018Aug
   !#
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
   use dump
   
   use node_mod, only:  &
       mynum, &
       mchnum, &
       master_num
       
   implicit none

   private

   !parameters
   double precision,parameter :: one = 1.d0
   !#  define double precision 1.0 to be multiplied by literal constants
   !#  that are embedded within parenthes
   double precision,parameter :: small_pc = 1.d-25
   !# define small particle number concentration [ # / x_units / y_units / z_units ]
   double precision,parameter :: almost_zero = 1.d-7
   !#  define smallest possible number such that one + almost_zero > one
   double precision,parameter :: almost_one = one - almost_zero
   !#  define value slightly < 1.0 to test for total evaporation [ dimensionless ]
   !double precision,parameter :: epsilon=almost_zero
   double precision,parameter :: epsilon=1e-12
   !#  logical unit number for output print file
   
   integer, parameter :: nsitesaot  = 12
   
   logical, parameter :: plotGrads=.false.
   logical, parameter :: justTest=.false.
   logical, parameter :: testHomog=.false.
   
   integer, parameter :: im = selected_int_kind(6)
   !# 4 byte integer
   integer, parameter :: rb = selected_real_kind(12)
   !# 8 byte real
   integer(kind=im), parameter :: nbndsw = 14
   integer(kind=im), parameter :: nbndlw = 16   !KML
   integer(kind=im), parameter :: nwt=nbndsw+nbndlw   
!KML--------
   real :: casee(nwt,nsitesaot)
   real :: casew(nwt,nsitesaot)
   real :: caseg(nwt,nsitesaot)

!KML-------
   
   real :: caser(nsitesaot)
   real :: cased(nsitesaot)
   real :: cases(nsitesaot)
   
   type aotMap_t
      real, pointer, dimension(:,:)    :: aotMap
      integer, pointer, dimension(:,:) :: currSite
      real, pointer, dimension(:,:,:)  :: extCoef
      real, pointer, dimension(:,:)  :: r0
      real, pointer, dimension(:,:)  :: pdens
      real, pointer, dimension(:,:)  :: rsig
      !# geometric standard deviation
   end type aotMap_t
   type(aotmap_t), allocatable, dimension(:) :: opt_aotmap, &
                  opt_aotMapm
   
   type aod_t
      real, pointer, dimension(:,:,:) :: tauaer
      !# Aod
      real, pointer, dimension(:,:,:) :: ssa
      !# single scattering albedo
      real, pointer, dimension(:,:,:) :: asp
      !# asymmetry parameter
   end type aod_t
   type(aod_t), allocatable, dimension(:,:) :: optProp
   !# tauaer in each cell (ngrid,nbands)
   
   logical :: aodFirstTime
   !# initialization control
   character(len=*),parameter :: mapAotFile = './tables/rad_carma/infMapAOT.vfm'
   character(len=*),parameter :: opticalPropertiesFile='./tables/rrtmg/opticalProperties.csv'
   !# File with AOT data
   
   character(len=*),parameter :: author='LFR;KML;NER'
   character(len=*),parameter :: sourceName='optical.f90'
   
   
   public aodDriver,aodFirstTime,setOptMemory,optProp

   
   contains


   subroutine aodDriver(m1,m2,m3,ia,iz,ja,jz,ngrids)
      !# Driver to compute tauaer to RRTGM
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: Compute tauaer for all points. If first Time compute the site information in all
      !# grid according aot map and/or path area. The tauaer is return back to RRTMG 
      !# using the memory type optProp
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#             Longo, K.M. **&#9993;**<mailto:karla.longo@inpe.br>
      !#             Rosario, N. E. **&#9993;**<niltoncvbr@gmail.com>      !#
      !#
      !# **Date**: 2020Jul
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
      !
      use mem_leaf, only: &
                     leaf_g
      use aer1_list, only: &
                     mode_alloc, &
                     nspecies,    &! Number of aerosol species
                     nmodes,      &! Number of aerosol modes
                     spc_alloc,   &! Species table allocation
                     part_radius, &! Radius of particles
                     part_dens, &  ! Density of particles (rho)
                     transport, &
                     accum,&
                     bburn, &
                     coarse, &
                     urban, &
                     on, &
                     spc_alloc_aer  =>spc_alloc
      use ccatt_start, only: &
                     ccatt
      use mem_aer1, only: &
                     aer1_g, & !Aerosol data
                     aerosol
      use mem_grid, only: &
                     grid_g, &
                     dzt, &
                     dtlt, &
                     time, &
                     ngrid
      use mem_basic, only: &
                     basic_g
      use mem_radiate, only:        &
                     radfrq, &
                     iswrtyp
      use dump
      
      integer, parameter :: aer_dir_effect = 1
      real,parameter :: fcui=1.e-9
      !# de billion g [gas/part] /kg [ar] para kg/kg
      integer, parameter:: ntotal_aer = 3 
      character(len=*),parameter :: procedureName='aodDriver'
      
      
      integer, intent(in) :: m1
      !# number of vertical levels
      integer, intent(in) :: m2
      !#
      integer, intent(in) :: m3
      !#
      integer, intent(in) :: ia
      !#
      integer, intent(in) :: iz
      !#
      integer, intent(in) :: ja
      !#
      integer, intent(in) :: jz
      !#
      integer, intent(in) :: ngrids
      !# Current grid
      
      include "constants.f90"
      integer :: n_aer, nb
      integer :: i,j,k,kr
      real :: pmr(m1,m2,m3)
!KML
      real :: tauaer(m1,m2,m3,nwt)
!KML
      real :: totm(m2,m3,m1)
      real :: caer(m2,m3,m1)
      
      !--- radiation calculation is updated only every radfrq seconds
      IF ( .not. (mod(time+.001, radfrq) < dtlt .or. time<0.001) .or. iswrtyp/=6) return
      
      totm=0.0
      tauaer=0.0

      ! Prepare some variables and read aotMap.
      ! !Only in the first call because the site have the same type
      if(aodFirstTime) then 
         !Fill the arrays with optical characteristics
         call setupraddata() 
         !Read aotMap
         call opt_read_aotmap()
         
         !Adjust the site accordongly aot map end veg patch area
         !compute particles R0 and Particles density
         !Compute Extinction Coeficient
         call initCurrSite(ia,iz,jz,jz,m1,m2,m3,ngrid &
                           ,opt_aotMap(1)%aotmap &
                           ,opt_aotMap(1)%currSite &
                           ,leaf_g(1)%patch_area(:,:,1) &
                           ,opt_aotMap(1)%extCoef &
                           ,opt_aotMap(1)%r0 &
                           ,opt_aotMap(1)%pdens &
                           ,opt_aotMap(1)%rsig &
                           ) !Init the site type
         
         aodFirstTime=.false.
      endif
      
      !Get the total mass concentration of aerosol and put it in pmr 
      pmr=0.0
      do n_aer=1,ntotal_aer
         if(aer_dir_effect == 1) then
            if (n_aer==1 ) then
               if(CCATT==1 .and. AEROSOL ==1 .and. spc_alloc_aer(transport,accum,bburn) == ON) then
                  pmr( 1:10,:,:)=aer1_g(accum,bburn,ngrid)%sc_p( 1:10,:,:)*fcui + 2.e-9  ! plus 10ug/m^3
                  pmr(11:20,:,:)=aer1_g(accum,bburn,ngrid)%sc_p(11:20,:,:)*fcui + 4.e-10 ! plus 5ug/m^3
                  pmr(21:m1,:,:)=aer1_g(accum,bburn,ngrid)%sc_p(21:m1,:,:)*fcui
               endif
            elseif (n_aer==2 ) then
               if(CCATT==1 .and. AEROSOL ==1  .and. spc_alloc_aer(transport,accum,urban) == ON) then
                  pmr(:,:,:)=pmr(:,:,:)+aer1_g(accum,urban,ngrid)%sc_p(:,:,:)*fcui
               endif
            elseif (n_aer==3 ) then
               if(CCATT==1 .and. AEROSOL ==1 .and. spc_alloc_aer(transport,coarse,urban) == ON) then
                  pmr(:,:,:)=pmr(:,:,:)+aer1_g(coarse,urban,ngrid)%sc_p(:,:,:)*fcui
               endif
            endif
         endif
      enddo
         
      !Compute total mass particle concentration (kg/m3)
      do k = 1,m1
            do i=ia,iz
               do j=ja,jz
                  totm(i,j,k) = pmr(k,i,j)*basic_g(ngrid)%dn0(k,i,j)
               enddo
            enddo
      enddo
         
      if(justTest) then
            if (mynum>1) iErrNumber &
            =dumpMessage(c_tty,c_yes,sourceName,procedureName,c_fatal, &
            'For test You must run with only one processor!')
            call fillTest(m1,m2,m3,totm,n_aer)
      endif

      ! Computes tauaer 
      call computeTauaer(m1,m2,m3,ia,iz,ja,jz,totm,tauaer,caer,n_aer)
         
      !Copy accum tauaer to a memory type - it will be used in RRTMG
!KML
      do nb=1,nwt
            optProp(1,nb)%tauaer(:,:,:)=tauaer(:,:,:,nb)
      enddo
!KML
  
      !If is to plot grads (see the constant above) plot a grads file
      if (plotGrads) call gradsTauaer(m1,m2,m3,ntotal_aer,totm)

   end subroutine aodDriver
  

   subroutine computeTauaer(m1,m2,m3,ia,iz,ja,jz,totm,tauaer,caer,n_aer)
      !# Compute the tauaer
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: this subroutine transform the mass concentration of
      !# aerosol in number and calculate the tauaer for nbands and level.
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#             Longo, K.M. **&#9993;**<mailto:karla.longo@inpe.br>
      !#             Rosario, N. E. **&#9993;**<niltoncvbr@gmail.com>      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 2020Jul
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
      !
      use mem_grid, only: &
                     dtlt, &
                     time, &
                     dzt, &
                     grid_g, &
                     ngrid
      use mem_radiate, only:        &
                     radfrq
      use mem_leaf, only: &
                     leaf_g
      
      implicit none
      include "constants.f90"
      
      integer,intent(in) :: m1,m2,m3
      !# Size of dimensions of variables for this proc
      integer, intent(in) :: ia,iz,ja,jz
      !# limits os lon and lat directions for this processor
      integer, intent(in) :: n_aer
      
      real,intent(inout) :: totm(m2,m3,m1)
      !# total mass particle concentration (kg/m3)
!KML
      real,intent(out)    :: tauaer(m1,m2,m3,nwt)
!KML
      !#
      real,intent(out) :: caer(m2,m3,m1)
      !# 
      
      integer :: i,j,nw,nl,ng,nr,nt,k,nb,kr,inb
      integer :: icount,kk

      real :: totn(m1,m2,m3)
      real :: arg1,arg2
      real :: sump(m2,m3)
      real :: r0(m2,m3)
      real :: pdens(m2,m3)
!KML
      real :: qex(m2,m3,nwt)

!KML
      real :: rsig(m2,m3)
      real :: xymet
      real :: dztr
      real :: xsecta
!KML
      real :: tau(m1,m2,m3,nwt)
!KML      
      character(len=8) :: ctime
      character :: cn_aer
      
!KML
      !integer, dimension(nbndsw) :: adjust=(/12,11,10,9,8,6,5,4,4,3,2,1,1,12/)
!      integer, dimension(nbndsw) :: adjust=(/1,2,3,4,5,6,7,8,9,10,11,12,13,14/)      
!KML      
      !Pega os valores de qext, r0, rsig e pdens obtidos em initCurrSite
      !e coloca em variaveis locais para melhorar a visualizacao
      do i=ia,iz
         do j=ja,jz
            r0(i,j)=opt_aotMap(1)%r0(i,j) !From caser
            pdens(i,j)=opt_aotMap(1)%pdens(i,j) !From cased
            rsig(i,j)=opt_aotMap(1)%rsig(i,j) !From cases
!KML
            do nb=1,nwt
               qex(i,j,nb)=opt_aotMap(1)%extCoef(i,j,nb) !Form casee
	    enddo
!KML	    
         enddo
      enddo
      !    
      !Calcula o concentracao de numero de particulas por ponto de grade #/m2
      do k = 1,m1
         kr=k+1
         do i=ia,iz
            do j=ja,jz
               totn(k,i,j)=(6.*totm(i,j,k)/(pdens(i,j)*c_pi*r0(i,j)**3)) &
                        *exp((-9./2)*log(rsig(i,j))**2)
            end do
         enddo
      enddo

      do k=1,m1
         do i=ia,iz
            do j=ja,jz
               dztr=dzt(k)/grid_g(ngrid)%rtgt(i,j) !* 1.e-2
               xymet = grid_g(ngrid)%fmapt(i,j)*grid_g(ngrid)%fmapt(i,j)
               caer(i,j,k) = totn(k,i,j)*(1./dztr)/xymet
            end do
         end do
      enddo
  
      !Calcula tauaer onde T=ExtCoef*Pi*R**2*Npart
      do i=ia,iz
         do j=ja,jz
            xsecta = c_pi * r0(i,j)**2.
            do k=1,m1
!KML
               do  nb = 1,nwt
                  tau(k,i,j,nb)=qex(i,j,nb)*xsecta*caer(i,j,k)
                  tauaer(k,i,j,nb)=max(real(epsilon),tau(k,i,j,nb))
               enddo

            end do
         end do
      end do

      
if(justTest) then 
   write(ctime,fmt='(I8.8)') int(time)
   write(cn_aer,fmt='(I1)') n_aer
   open(unit=33,file='dump'//ctime//'-'//cn_aer//'.txt') 
      
   write(33,fmt='(5(A2,1X),10(A13,1X))') 'nb','i','j','k','Site' &
           ,'xymet','dztr','Tot mass','Tot Number','caer','extCoef','r0','pdens','tau','tauaer'    
   icount=0
   do nb=1,nbndsw
      do k=1,m1
         do i=ia,iz
            do j=ja,jz
               if(nb>9 .and. nb<14) then
                  icount=icount+1
                  if(icount==39) then
                     icount=0
                     write(33,fmt='(5(A2,1X),10(A13,1X))') 'nb','i','j','k','Site' &
                     ,'xymet','dztr','Tot mass','Tot Number','caer','extCoef','r0','pdens','tau','tauaer' 
                  endif
                  write(33,fmt='(5(I2.2,1X),10(E13.6,1X))') nb &
                     ,i,j,k,opt_aotMap(1)%currSite(i,j) &
                     ,grid_g(ngrid)%fmapt(i,j)*grid_g(ngrid)%fmapt(i,j),dzt(k)/grid_g(ngrid)%rtgt(i,j) &
                     ,totm(i,j,k),totn(k,i,j),caer(i,j,k),qex(i,j,nb),r0(i,j),pdens(i,j),tau(k,i,j,nb),tauaer(k,i,j,nb)
               end if
            enddo
         enddo
      enddo
   enddo
   do i=ia,iz
      do j=ja,jz
         write(33,*) i,j,sum(totm(:,i,j)),sum(tauaer(:,i,j,10))
      enddo
   enddo
   close(unit=33)
endif    

   end subroutine computeTauaer


   subroutine initCurrSite(ia,iz,ja,jz,m1,m2,m3,ngrid,aotMap,currSite,patch_area &
                           ,qexti,r0,pdens,rsig)
      !# Initialize crossreference pointer to Qext tables
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: this subroutine prepare the pointer to Qext table
      !# Tropospheric Aerosol Optical Thickness from the GOCART Model and Comparisons
      !# with Satellite and Sun Photometer Measurements, Chin et al. 2002. and
      !# Regional representativity of AERONET observation sites during the biomass burning
      !# season in South America determined by correlation studies with MODIS Aerosol Optical Depth.
      !# Hoelzemann, J. J. and al. JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 114, D13301,
      !# doi:10.1029/2008JD010369, 2009
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#             Longo, K.M. **&#9993;**<mailto:karla.longo@inpe.br>
      !#             Rosario, N. E. **&#9993;**<niltoncvbr@gmail.com>      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 2018Aug
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
      !
      character(len=*),parameter :: procedureName='initCurrSite'
      include "constants.f90"
      integer, intent(in) :: ia
      !#
      integer, intent(in) :: iz
      !#
      integer, intent(in) :: ja
      !#
      integer, intent(in) :: jz
      !#
      integer, intent(in) :: m2
      !#
      integer, intent(in) :: m3
      !
      integer, intent(in) :: m1
      !
      integer, intent(in) :: ngrid
      !
      real,intent(in) :: patch_area(m2,m3)
      !# Path of soil characteritics
      real, intent(inout) :: aotMap(m2,m3)
      !# Aot map
      integer,dimension(m2,m3),intent(out) :: currSite
      !# Used to determine which table used for Ext. Coef.
!KML
      real,intent(out) :: qexti(m2,m3,nwt)
!KML
      real, intent(out) :: r0(m2,m3)
      real, intent(out) :: pdens(m2,m3)
      real, intent(out) :: rsig(m2,m3)
      
      integer :: i,j,k,nr,ns
      

      !do i=1,m2
      !   do j=1,m3
      !	 
      ! if(aotMap(i,j)<0 .or. aotMap(i,j)>11 ) write (88,*) i,j,mynum,aotMap(i,j);call flush(88)
      !    enddo
      ! enddo

      where( aotMap < 0. .or. aotMap > 10. ) aotMap = 11.
     
!      do i=1,m2
!         do j=1,m3
!        
!          if(aotMap(i,j)<0 .or. aotMap(i,j)>11 ) write (88,*) i,j,mynum,aotMap(i,j);call flush(88)
!          enddo
!       enddo

      !Loop in all Cols
      do i=1,m2
         do j=1,m3
	 
	         if(aotMap(i,j)<0 .or. aotMap(i,j)>11 ) write (88,*) i,j,mynum,aotMap(i,j);call flush(88)	

         !print*,"aotMap",i,j,mynum,aotMap(i,j);call flush(6)
	 
            if(int(aotMap(i,j)) > 0 .and. int(aotMap(i,j)) < 11) then !If not mapped inside list
               !Is one of the mapped
               currSite(i,j)=aotMap(i,j)
            else
               !If not mapped inside list
               if (patch_area(i,j) .gt. .009) then !Is ocean
                  currSite(i,j)=11
               else !Is continental
                  currSite(i,j)=12
               endif
            endif
!KML
!            do ns=1,nsitesaot
!LFR               r0(i,j)   =caser(currsite(i,j))
               r0(i,j)=2.62533E-7
!Valor acima foi alterado para ficar fixo ateh que se tenha um valor
!adequado lido do CASER!
               pdens(i,j)=cased(currsite(i,j))
               rsig(i,j) =cases(currsite(i,j))
!            enddo

            do nr=1,nwt
               qexti(i,j,nr) = caseE(nr,currSite(i,j))
               do k=1,m1
                  optProp(ngrid,nr)%ssa(k,i,j)=caseW(nr,currSite(i,j))
                  optProp(ngrid,nr)%asp(k,i,j)=caseG(nr,currSite(i,j))
               enddo
            enddo
!KML
            if(rsig(i,j)<1.e-12) then 
	    
	      print*,"rsig",i,j,rsig(i,j),mynum
	      call flush(6)
	      stop 333
	    endif
	    
	    
         enddo
      enddo

   end subroutine initCurrSite


   subroutine setupraddata()
      !# Fill the cases array values
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: Fill the cases array values
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#             Longo, K.M. **&#9993;**<mailto:karla.longo@inpe.br>
      !#             Rosario, N. E. **&#9993;**<niltoncvbr@gmail.com>      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 2018Aug
      !#
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

      use dump
      
      include "constants.f90"
      character(len=*),parameter :: procedureName='setupraddata'
      logical, parameter :: inventory=.true.

      character(len=20) :: text
      character(len=256) :: header
      integer :: iValue,iSites,iBands,nS,nB
      logical :: there
      real, allocatable :: wavenumberIni(:)
      real, allocatable :: wavenumberEnd(:)
      real, allocatable :: wavelenghtIni(:)
      real, allocatable :: wavelenghtEnd(:)
      real, allocatable :: wavelenghtCen(:)
      
      
      inquire(file=trim(opticalPropertiesFile),exist=there)
      if(.not.there) then
         iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName,c_fatal, &
          "File '"//trim(opticalPropertiesFile)//' needed for RTMG radiation not found. ')
      endif
      
      if(inventory) open(unit=34,file='radInventory.dat',status='replace' &
           ,action='WRITE',form='FORMATTED')
      
      open(unit=33,file=trim(opticalPropertiesFile))
      
      !Read Header Information
      read(33,*) header
      if(inventory) write(34,*) header
      
      !Read the number of Aeronet Sites
      read(33,*) text,ivalue
      if(trim(text)/='sites') iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName,c_fatal, &
          "Snd line of opticalPropertiesFile must be number of sites 'sites'")
      iSites=iValue
      if(inventory) write(34,*) 'SItes=',iSites
      
      !Read the number of sw bands
      read(33,*) text,ivalue
      if(trim(text)/='bands') iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName,c_fatal, &
          "Trd line of opticalPropertiesFile must be number of radiation bands 'bands'")
      ibands=iValue 
      if(inventory) write(34,*) 'Bands=',iBands
      
      !Check if #of bands is equal of RRTMG
!KML
!      if(ibands/=nbndsw) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName,c_fatal, &
!          "Number of bands read is different of RRTMG bands, ",nbndsw,'I2.2')

      if(ibands/=nwt) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName,c_fatal, &
          "Number of bands read is different of RRTMG bands, ",nwt,'I2.2')
!KML   
      allocate(wavenumberIni(iBands))
      allocate(wavenumberEnd(iBands))
      allocate(wavelenghtIni(iBands))
      allocate(wavelenghtEnd(iBands))
      allocate(wavelenghtCen(iBands))
      
      !Readind wavenumber Initial
      read(33,*) text,wavenumberIni
      if(inventory) write(34,*) 'wavenumberIni=',wavenumberIni
      !Reading wavenumber Final
      read(33,*) text,wavenumberEnd
      if(inventory) write(34,*) 'wavenumberEnd=',wavenumberEnd
      !Reading wavelenght Initial
      read(33,*) text,wavelenghtIni
      if(inventory) write(34,*) 'wavelenghtIni=',wavelenghtIni
      !Reading wavelenght Final
      read(33,*) text,wavelenghtEnd
      if(inventory) write(34,*) 'wavelenghtEnd=',wavelenghtEnd
      !Compute wavelenght Center
      wavelenghtCen=(wavelenghtIni+wavelenghtEnd)/2.0
      if(inventory) write(34,*) 'wavelenghtCen=',wavelenghtCen
      
      !Skip 2 lines
      read(33,*) text
      read(33,*) text
      
      !Reading casee (Extinction Coef)
      do ns=1,iSites
         read(33,*) text,casee(:,ns)
      enddo
      
      do ns=1,iSites
         if(inventory) write(34,*) 'casee(',ns,')',casee(:,ns)
      enddo
      !Skip 1 lines
      read(33,*) text
      !Reading casew (Single scat. Albedo)
      do ns=1,iSites
         read(33,*) text,casew(:,ns)
      enddo
      do ns=1,iSites
         if(inventory) write(34,*) 'casew(',ns,')',casew(:,ns)
      enddo
      !Skip 1 lines
      read(33,*) text
      !Reading caseG (Assim. Parameter)
      do ns=1,iSites
         read(33,*) text,caseG(:,ns)
      enddo
      do ns=1,iSites
         if(inventory) write(34,*) 'caseg(',ns,')',caseg(:,ns)
      enddo
      !Skip 1 lines
      read(33,*) text
      !Reading caseR (Efec. Radius)
      do ns=1,iSites
         read(33,*) text,caseR(ns)
      enddo  
      if(inventory) write(34,*) 'caser:',caser
      !Skip 1 lines
      read(33,*) text
      !Reading caseS (Sigma)
      do ns=1,iSites
         read(33,*) text,caseS(ns)
      enddo 
      if(inventory) write(34,*) 'cases:',cases
      !Skip 1 lines
      read(33,*) text
      !Reading caseD (Density)
      do ns=1,iSites
         read(33,*) text,caseD(ns)
      enddo 
      if(inventory) write(34,*) 'cased:',cased
      close(unit=33)    
      if(inventory) close(unit=34)
      
   end subroutine setupraddata
  
   subroutine setOptMemory(ngrids,imean,nmzp,nmxp,nmyp)
      !# Allocate and set aotMapvariables in memory
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: Sets the variable aotmap in memory. Must be called in ModeOneproc
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#             Longo, K.M. **&#9993;**<mailto:karla.longo@inpe.br>
      !#             Rosario, N. E. **&#9993;**<niltoncvbr@gmail.com>      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 2018Aug
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
      !
      integer, intent(in) :: ngrids
      !# Total of grids
      integer, intent(in) :: imean
      !# Mean type of variables
      integer, dimension(ngrids) :: nmzp
      !# Levels in each grid
      integer, dimension(ngrids) :: nmxp
      !# Lon pointes in each grid
      integer, dimension(ngrids) :: nmyp
      !#Lat points in each grid
      
      integer :: ierr
      integer :: ng
      integer :: band
      character(len=*),parameter :: h='(**setOptMemory**)'
      
      !print *,'i. nWave=',nWave
!KML      !print *,'Allocating memory for optical, grids=',ngrids,', bands=',nwt; call flush(6)
      allocate(opt_aotMap(ngrids), STAT=ierr)
      if (ierr/=0) call fatal_error(h//"Error allocating aotMap in optical")
      allocate(opt_aotMapm(ngrids), STAT=ierr)
      if (ierr/=0) call fatal_error(h//"Error allocating aotMapm in optical")
!KML
!      allocate(optProp(ngrids,nbndsw), STAT=ierr)
      allocate(optProp(ngrids,nwt), STAT=ierr)
!KML
      if (ierr/=0) call fatal_error(h//"Error llocating aod in optical")
      do ng=1,ngrids
      call opt_nullify_aotMap(opt_aotMap,ng)
      call opt_nullify_aotMap(opt_aotMapm,ng)
      call opt_alloc_aotMap  (opt_aotMap(ng),nmxp(ng), nmyp(ng))
      call opt_alloc_aotMap  (opt_aotMapm(ng),nmxp(ng), nmyp(ng))
      call opt_filltab_aotMap(opt_aotMap(ng), opt_aotMapm(ng), ng, imean,nmxp(ng), nmyp(ng) )
      !
!KML
!      do band=1,nbndsw
      do band=1,nwt
         call opt_nullify_aod(optProp,ng,band)
         call opt_alloc_aod  (optProp(ng,band),nmzp(ng),nmxp(ng), nmyp(ng))
      enddo
      
      enddo
      
      aodFirstTime=.true.

   end subroutine setOptMemory

   subroutine opt_alloc_aod(c_aod,nz,nx,ny)
      !# Allocate  aod variable
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: Allocate and initilize to zero the variable aod
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#             Longo, K.M. **&#9993;**<mailto:karla.longo@inpe.br>
      !#             Rosario, N. E. **&#9993;**<niltoncvbr@gmail.com>      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 2018Aug
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
      !
      type(aod_t) :: c_aod
      !# Variable to be allocated Map of AOT
      integer, intent(in) :: nz
      !# Number of levels
      integer, intent(in) :: nx
      !# Number of lats
      integer, intent(in) :: ny
      !# Number of lons
      allocate(c_aod%tauaer(nz,nx,ny))
      allocate(c_aod%ssa(nz,nx,ny))
      allocate(c_aod%asp(nz,nx,ny))

      c_aod%tauaer(:,:,:) = 0.
      c_aod%ssa(:,:,:) = 0.
      c_aod%asp(:,:,:) = 0.
      
   end subroutine opt_alloc_aod

   subroutine opt_nullify_aod(aod_c, ng, band)
      !# Nullify  aod variable
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: Allocate and initiliza to zero the variable aod
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#             Longo, K.M. **&#9993;**<mailto:karla.longo@inpe.br>
      !#             Rosario, N. E. **&#9993;**<niltoncvbr@gmail.com>      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 2018Aug
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
      !
      implicit none
      
      type(aod_t), intent(inout) :: aod_c(:,:)
      !# Type to be nullified
      integer, intent(in)          :: ng
      !# Grid
      integer, intent(in)          :: band

      !# band number

      if (associated(aod_c(ng,band)%tauaer))  nullify (aod_c(ng,band)%tauaer)
      if (associated(aod_c(ng,band)%asp))     nullify (aod_c(ng,band)%asp)
      if (associated(aod_c(ng,band)%ssa))     nullify (aod_c(ng,band)%ssa)

   end subroutine opt_nullify_aod

   subroutine opt_alloc_aotMap(c_aotMap,nx,ny)
      !# Allocate  aotMap variable
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: Allocate and initiliza to zero the variable aotmap
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#             Longo, K.M. **&#9993;**<mailto:karla.longo@inpe.br>
      !#             Rosario, N. E. **&#9993;**<niltoncvbr@gmail.com>      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 2018Aug
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
      !
      type(aotMap_t) :: c_aotMap
      !# Variable to be allocated Map of AOT
      integer, intent(in) :: nx
      !# Number of lats
      integer, intent(in) :: ny
      !# Number of lons
      
      allocate(c_aotMap%aotMap(nx,ny))
      allocate(c_aotMap%currSite(nx,ny))
!KML
      allocate(c_aotMap%extCoef(nx,ny,nwt))
!KML
      allocate(c_aotMap%r0(nx,ny))
      allocate(c_aotMap%pdens(nx,ny))
      allocate(c_aotMap%rsig(nx,ny))
      
      c_aotMap%aotMap(:,:) = 0.
      c_aotMap%currSite(:,:)=0
      c_aotMap%extCoef(:,:,:)=0.
      c_aotMap%r0=0.
      c_aotMap%pdens=0.
      c_aotMap%rsig=0.
      
   end subroutine opt_alloc_aotMap

  subroutine opt_nullify_aotMap(aot, ng)
      !# Nullify  aotMap variable
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: Allocate and initiliza to zero the variable aotmap
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#             Longo, K.M. **&#9993;**<mailto:karla.longo@inpe.br>
      !#             Rosario, N. E. **&#9993;**<niltoncvbr@gmail.com>      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 2018Aug
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
      !
      implicit none
      
      type(aotMap_t), intent(inout) :: aot(:)
      !# Type to be nullified
      integer, intent(in)          :: ng
      !# Grid
      if (associated(aot(ng)%aotMap))  nullify (aot(ng)%aotMap)
      if (associated(aot(ng)%currSite))  nullify (aot(ng)%currSite)
      if (associated(aot(ng)%extCoef))  nullify (aot(ng)%extCoef)
      if (associated(aot(ng)%r0))  nullify (aot(ng)%r0)
      if (associated(aot(ng)%pdens))  nullify (aot(ng)%pdens)
      if (associated(aot(ng)%rsig))  nullify (aot(ng)%rsig)

   end subroutine opt_nullify_aotMap

   subroutine opt_dealloc_aotMap()
      !# Dellocate  aotMap variable
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: Dellocate and the variable aotmap
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#             Longo, K.M. **&#9993;**<mailto:karla.longo@inpe.br>
      !#             Rosario, N. E. **&#9993;**<niltoncvbr@gmail.com>      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 2018Aug
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
      !
      use mem_grid, only: &
         ngrids,         &!(in)
         nnxp,           &!(in)
         nnyp,           &!(in)
         nnzp             !(in)
      
      
      integer :: idx
      
      if(ALLOCATED(opt_aotMap))then
      do idx = 1, ngrids
         deallocate(opt_aotMap(idx)%aotMap)
         deallocate(opt_aotMap(idx)%currSite)
         deallocate(opt_aotMap(idx)%extCoef)
         deallocate(opt_aotMap(idx)%r0)
         deallocate(opt_aotMap(idx)%pdens)  
         deallocate(opt_aotMap(idx)%rsig)  
      end do
      deallocate(opt_aotMap)
      end if
      
      if(allocated(opt_aotMapm))then
      do idx = 1, ngrids
         deallocate(opt_aotMapm(idx)%aotMap)
         deallocate(opt_aotMapm(idx)%currSite)
         deallocate(opt_aotMapm(idx)%extCoef)
         deallocate(opt_aotMapm(idx)%r0)
         deallocate(opt_aotMapm(idx)%pdens)
         deallocate(opt_aotMapm(idx)%rsig)
      end do
      deallocate(opt_aotMapm)
      end if

   end subroutine opt_dealloc_aotMap

   subroutine opt_filltab_aotMap(imap, imapm, ng, imean, n1, n2)
      !# Fill the  aotMap variable in var tables of model
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: Fill the  aotMap variable in var tables of model
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#             Longo, K.M. **&#9993;**<mailto:karla.longo@inpe.br>
      !#             Rosario, N. E. **&#9993;**<niltoncvbr@gmail.com>      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 2018Aug
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
      !
      use var_tables, only: InsertVTab
      include "i8.h"
      integer, intent(in) :: ng
      !# grid
      integer, intent(in) :: n1
      !# Lon number of points
      integer, intent(in) :: n2
      !# Lat number of points
      integer, intent(in) :: imean
      !# Mean of variable
      type(aotMap_t) :: imap
      !# Var to be inserted
      type(aotMap_t) :: imapm
      !# Var to be inserted (m)
      
      integer(kind=i8)       :: npts
      !# Total of points 1D
      
      if(associated(imap%aotMap))then
         npts = n1*n2
         call InsertVTab(imap%aotMap, imapm%aotMap, ng, npts, imean, &
               'AOTMAP :2:hist:anal:mpti')
      end if

   end subroutine opt_filltab_aotMap

   subroutine opt_read_aotmap()
      !# Read  aot Map from input file
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: Read aot map from file. The file contas the site type according with:
      !# 00 undef
      !# 01 Abracos_Hill.lev15
      !# 02 CUIABA-MIRANDA.lev15
      !# 03 Rio_Branco.lev15
      !# 04 Alta_Floresta.lev15
      !# 05 Campo_Grande_SONDA.lev15
      !# 06 CEILAP-BA.lev15
      !# 07 Cordoba-CETT.lev15
      !# 08 SANTA_CRUZ.lev15
      !# 09 Sao_Paulo.lev15
      !# 10 belterra.
      !#
      !# For the following we set inside model:
      !# 
      !# 11 1 oceanico
      !# 12 continental
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#             Longo, K.M. **&#9993;**<mailto:karla.longo@inpe.br>
      !#             Rosario, N. E. **&#9993;**<niltoncvbr@gmail.com>      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 2018Aug
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
      !
      use mem_grid
      use node_mod, only: &
            nodei0, &
            nodej0, & ! intent(in)
            nodemxp, &
            nodemyp, & ! intent(in)
            nmachs,  & ! intent(in)
            mynum,  &  ! intent(in)
            mchnum, &  ! intent(in)
            master_num, & ! intent(in)
            ia, &
            iz, &
            ja, &
            jz, &
            mxp, &
            myp, &
            mzp, &
            nxbeg, &
            nxend, &
            nybeg, &
            nyend
      use parlib, only: &
            parf_bcast ! subroutine
      use readbcst, only: &
         gatherdata
      
      integer :: i
      integer :: j
      integer :: nlon
      integer :: nlat
      integer :: ii
      integer :: jj
      integer :: effsite
      integer :: nsites
      integer :: ifm
      integer :: i1
      integer :: i2
      integer :: ic
      integer :: j1
      integer :: j2
      integer :: jc
      integer :: qi1
      integer :: qi2
      integer :: qj1
      integer :: qj2
      real    :: latni
      real    :: latnf
      real    :: lonni
      real    :: lonnf
      real    :: latstep
      real    :: lonstep
      real    :: dlonr
      real    :: dlatr
      real    :: undef
      real    :: comm(10)
      integer, allocatable, dimension(:,:) :: infmap
      real, allocatable, dimension(:,:)    :: infmapreal
      real, allocatable, dimension(:)      :: moda
      real, allocatable, dimension(:,:)    :: prlon
      real, allocatable, dimension(:,:)    :: prlat
      real :: globalglon(nnxp(1), nnyp(1)) !apenas para uma grade teste
      real :: globalglat(nnxp(1), nnyp(1)) !apenas para uma grade teste
      real :: globalaot(nnxp(1), nnyp(1))
      character(len=16)  :: varn
      integer, external :: outRealSize
      integer :: recordLen
      
      if (mchnum==master_num) then !Only master do the read
      open(unit=17, file=trim(mapaotfile), status='old')
      !print *,'Reading the header of AOTMap vfm file: ', trim(mapaotfile); call flush(6)
      read(17,*) lonni, latni, lonnf, latnf
      read(17,*) nlon, nlat, lonstep, latstep
      read(17,*) nsites, undef
      !Prepare array to be broadcast for all processors
      comm(1) =lonni
      comm(2) =latni
      comm(3) =lonnf
      comm(4) =latnf
      comm(5) =lonstep
      comm(6) =latstep
      comm(7) =undef
      comm(8) =real(nlon)
      comm(9) =real(nlat)
      comm(10)=real(nsites)
      end if
      !print *, 'broadcasting AOT Map data'; call flush(6)
      call parf_bcast(comm, int(size(comm,1),i8), master_num)
      !Getting information from array
      lonni    =comm(1)
      latni    =comm(2)
      lonnf    =comm(3)
      latnf    =comm(4)
      lonstep  =comm(5)
      latstep  =comm(6)
      undef    =comm(7)
      nlon     =int(comm(8) )
      nlat     =int(comm(9) )
      nsites   =int(comm(10))
      !print *, 'Allocating information'; call flush(6)
      allocate(infmap(nlon, nlat),infmapreal(nlon, nlat))
      if (mchnum==master_num) then !Only master read
         !print *,'Reading the data inside AOTMap file'; call flush(6)
         call viirec(17, infmap, nlon*nlat, 'lin')
         !Make it real
         infmapreal=real(infmap)
         close(17)

         recordLen=outRealSize()*nlon*nlat
         open(unit=33,file='aotMap_check.gra',&
           action='WRITE',status='REPLACE',form='UNFORMATTED',access='DIRECT', &
           recl=recordLen)
         write (33,rec=1) infmapreal(:,:)
         close(unit=33)

         open(unit=33,file='aotMap_check.ctl' &
           ,action='WRITE',status='replace',form='FORMATTED')

         !writing the name of grads file
         write(33,*) 'dset ^aotMap_check.gra'
         !writing others infos to ctl
         write(33,*) 'undef -0.9990000E+34'
         write(33,*) 'title AOTMAP INPUT ',nsites
         write(33,*) 'xdef ',nlon,' linear ',lonni,lonstep
         write(33,*) 'ydef ',nlat,' linear ',latni,latstep
         write(33,*) 'zdef ',1,'levels',1000.
         write(33,*) 'tdef 1 linear 00:00z01jan2018     1mo'
         write(33,*) 'vars ',1
         write(33,*) 'map',1,'99 ','aot map'
         write(33,*) 'endvars'
   
         close(33)


      end if
      !print *, 'broadcasting AOT data'; call flush(6)
      call parf_bcast(infmapreal, int(size(infmap,1),i8), int(size(infmap,2),i8), master_num)
      !Make it int again
      infmap=int(infmapreal)
      !print *,'Processing local data'; call flush(6)
      allocate(prlat(nlon, nlat), prlon(nlon,nlat), moda(nsites))
      call apiprlatlon(nlon, nlat, prlat, prlon, latstep, lonstep, latni, lonni)
      do ifm =1, ngrids
      opt_aotMap(ifm)%aotmap(:,:) = 0.
      call newgrid(ifm)
      varn = 'glon'
      call gatherdata(2, varn, ifm, nnxp(ifm), nnyp(ifm), &
            nmachs, mchnum, mynum, master_num,             &
               grid_g(ifm)%glon, globalglon)
      varn = 'glat'
      call gatherdata(2, varn, ifm, nnxp(ifm), nnyp(ifm), &
            nmachs, mchnum, mynum, master_num,             &
            grid_g(ifm)%glat, globalglat)

      globalaot=0

      do i =1, nnxp(ifm)
         do j=1, nnyp(ifm)
            ! evitando pontos fora do dominio da grade
            if ( globalglat(i,j)<latni .or. globalglat(i,j)>latnf .or. &
               globalglon(i,j)<lonni .or. globalglon(i,j)>lonnf) cycle
            call interpolacao(globalglon(i,j), globalglat(i,j), nlon, nlat, &
                              prlat, prlon, i1, i2, ic, j1, j2, jc)
            if(ic.ge.0 .and. jc .ge. 0)then
               dlonr = 0.5*(globalglon(nodemxp(mynum,ifm),j) - globalglon(1,j))/&
                                                         float(nnxp(ifm)-1)
               dlatr = 0.5*(globalglat(i,nodemyp(mynum,ifm)) - globalglat(i,1))/&
                                                         float(nnyp(ifm)-1)
               qi1=int(dlonr/lonstep+0.5)
               qi2=int(dlonr/lonstep+0.5)
               qj1=int(dlatr/latstep+0.5)
               qj2=int(dlatr/latstep+0.5)
               moda(:) = 0
               do jj =max(1,jc-qj1),min(nlat,jc+qj2)
                  do ii = max(1,ic-qi1),min(nlon,ic+qi2)
                     if( infmap(ii,jj) .ne. undef)then
                        moda(infmap(ii,jj)) = moda(infmap(ii,jj)) + 1
                     end if
                  end do
               end do
               effsite = 1
               do ii = 1, nsites
                  do jj = 1, nsites-1
                     if(moda(ii) .gt. moda(jj))then
                        effsite = ii
                     end if
                  end do
               end do
               globalaot(i,j)= effsite !to remove 0 values
            end if
         end do
      end do
      
      if (mchnum==master_num) then !Only master test

         recordLen=outRealSize()*nnxp(1)*nnyp(1)
         open(unit=33,file='aotMap_global.gra',&
           action='WRITE',status='REPLACE',form='UNFORMATTED',access='DIRECT', &
           recl=recordLen)
         write (33,rec=1) globalaot(:,:)
         close(unit=33)

         open(unit=33,file='aotMap_global.ctl' &
           ,action='WRITE',status='replace',form='FORMATTED')

         !writing the name of grads file
         write(33,*) 'dset ^aotMap_global.gra'
         !writing others infos to ctl
         write(33,*) 'undef -0.9990000E+34'
         write(33,*) 'title AOTMAP INTERNAL ',nsites
         write(33,*) 'xdef ',nnxp(1),' linear ',globalglon(1,1),globalglon(2,1)-globalglon(1,1)
         write(33,*) 'ydef ',nnyp(1),' linear ',globalglat(1,1),globalglon(1,2)-globalglon(1,1)
         write(33,*) 'zdef ',1,'levels',1000.
         write(33,*) 'tdef 1 linear 00:00z01jan2018     1mo'
         write(33,*) 'vars ',1
         write(33,*) 'map',1,'99 ','aot map'
         write(33,*) 'endvars'
   
         close(33)


      end if


      !print *, "Making buffer for grid",ifm; call flush(6)
      call mk_2_buff(globalaot,opt_aotMap(ifm)%aotmap, &
                              nnxp(ifm), nnyp(ifm), &
                              nodemxp(mynum,ifm), nodemyp(mynum,ifm), &
                              nxbeg(mynum,ifm),nxend(mynum,ifm),nybeg(mynum,ifm),nyend(mynum,ifm))
      end do
      !print *, 'End AOT Reading'; call flush(6)
   end subroutine opt_read_aotMap

   subroutine gradsTauaer(m1,m2,m3,n_aer,aer1)
      use mem_grid, only: &
            grid_g, &
            time
      
      integer, intent(in) :: m1,m2,m3,n_aer
      real,intent(in) ::aer1(n_aer,m2,m3,m1)
      real :: somaAer(m2,m3),aot10(m2,m3),aot11(m2,m3),aot12(m2,m3),aot13(m2,m3)
      
      character(len=5) :: cnb
      character(len=6) :: ctime
      integer :: nb,k,recordLen,irec,i,j
      
      write(ctime,fmt='(I6.6)') int(time)
      
      open(unit=33,file='tauaer'//ctime//'.ctl' &
         ,action='WRITE',status='replace',form='FORMATTED')
      
      !writing the name of grads file
      write(33,*) 'dset ^tauaer'//ctime//'.gra'
      !writing others infos to ctl
      write(33,*) 'undef -0.9990000E+34'
      write(33,*) 'title tauaer'
      write(33,*) 'xdef ',m2,' linear ',grid_g(1)%glon(1,1),grid_g(1)%glon(2,1)-grid_g(1)%glon(1,1)
      write(33,*) 'ydef ',m3,' linear ',grid_g(1)%glat(1,1),grid_g(1)%glat(1,2)-grid_g(1)%glat(1,1)
      write(33,*) 'zdef ',m1,'levels ',(k,k=1,m1)
      write(33,*) 'tdef   1 linear 06:00z12jun2020           1hr'
      write(33,*) 'vars ',22

      do nb=1,nbndsw
      write(cnb,fmt='("Tau",I2.2)') nb
      write(33,*) cnb,m1,99,'Tauaer Band '//cnb
      enddo
      write(33,*) 'aer1',m1,99,'concentration accum,bburn (kg/m3)'
      write(33,*) 'aer2',m1,99,'concentration accum,urban (kg/m3)'
      write(33,*) 'csite',m1,99,'AotMap SIte'
      write(33,*) 'aerCol',m1,99,'Sum Aer Conc 1 & 2 at column'
      write(33,*) 'aot10',m1,99,'AOT column 533nm'
      write(33,*) 'aot11',m1,99,'AOT column 393nm'
      write(33,*) 'aot12',m1,99,'AOT column 303nm'
      write(33,*) 'aot13',m1,99,'AOT column 231nm'
      write(33,*) 'endvars'
      
      close(33)
      
      recordLen=4*m2*m3
      open(unit=33,file='tauaer'//ctime//'.gra',&
         action='WRITE',status='REPLACE',form='UNFORMATTED',access='DIRECT', &
         recl=recordLen)
      
      !# writing grads binary and fill variables
      irec=1
!KML      do nb=1,nbndsw
      do nb=1,nwt
         do k=1,m1
            write (33,rec=irec) optProp(1,nb)%tauaer(k,:,:)
            irec=irec+1
         enddo
!KML
	 
      enddo
      do k=1,m1
         write (33,rec=irec) aer1(1,:,:,k)
         irec=irec+1
      enddo
      do k=1,m1
         write (33,rec=irec) aer1(2,:,:,k)
         irec=irec+1
      enddo
      do k=1,m1
         write (33,rec=irec) real(opt_aotMap(1)%currSite(:,:))
         irec=irec+1
      enddo
      aot10=0.0
      aot11=0.0
      aot12=0.0
      aot13=0.0
      somaAer=0.0
      do i=1,m2
         do j=1,m3
            do k=1,m1-1
               somaAer(i,j)=somaAer(i,j)+aer1(1,i,j,k)+aer1(2,i,j,k)
               aot10(i,j)=aot10(i,j)+optProp(1,10)%tauaer(k,i,j)
               aot11(i,j)=aot11(i,j)+optProp(1,11)%tauaer(k,i,j)
               aot12(i,j)=aot12(i,j)+optProp(1,12)%tauaer(k,i,j)
               aot13(i,j)=aot13(i,j)+optProp(1,13)%tauaer(k,i,j)
            enddo
         enddo
      enddo
      do k=1,m1
         write (33,rec=irec) somaAer(:,:)
      enddo
      do k=1,m1
         write (33,rec=irec) aot10(:,:)
      enddo  
      do k=1,m1
         write (33,rec=irec) aot11(:,:)
      enddo 
      do k=1,m1
         write (33,rec=irec) aot12(:,:)
      enddo 
      do k=1,m1
         write (33,rec=irec) aot13(:,:)
      enddo     
      close(33) 
  
   end subroutine gradsTauaer


   subroutine fillTest(m1,m2,m3,totm,n_aer)
      !# Fill the  aer concentration for test of module response
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: If justTest is set the routine fills the totm (mass concentration)
      !# using fixed values:
      !# 1)If testHomog is true the concentration is setup to 1 ug/m3 in first level
      !#   and for next each level this ammount is increasing by 0.01 ug/m3
      !# 2)If testHomog is false the concentration goes from 0 to 1.5 ug/m3
      !#   in steps of 0.15ug/m3 and all levels are the same.
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author**: Rodrigues, L.F. **&#9993;**<mailto:luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 2018Aug
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
      !
   
      real, parameter :: factor=1.0e-9 !(convert ug/m3 to kg/m3)
      integer, intent(in) :: m1,m2,m3,n_aer
      real, intent(inout) :: totm(2,m2,m3,m1)
      
      integer :: i,j,k,icount
      real :: value
      
      totm(n_aer,:,:,:)=0.0
      icount=1
      
      if(testHomog) then
         do k=1,m1
            totm(n_aer,:,:,1:11)=30.0*factor
         enddo
         return
      else
         do i=5,35
            do j=5,35
               totm(n_aer,i,j,1:11)=7.5*factor
            enddo
         enddo
         do i=10,30
            do j=10,30
               totm(n_aer,i,j,1:11)=15.0*factor
            enddo
         enddo  
         do i=15,25
            do j=15,25
               totm(n_aer,i,j,1:11)=30.0*factor
            enddo
         enddo                  
      endif
   end subroutine fillTest

end module optical
