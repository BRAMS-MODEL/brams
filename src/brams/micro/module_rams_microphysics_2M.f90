!
! Copyright (C) 1991-2004  ; All Rights Reserved ; Colorado State University
! Colorado State University Research Foundation ; ATMET, LLC
!
! This file is free software; you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later version.
!
! This software is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with this
! program; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!======================================================================================
!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################


MODULE rams_microphysics_2M

use grid_dims, only: &
               NZPMAX, maxgrds

use micphys, only : &
    mcphys_type ,&
    aparm       ,&
    coltabfn    ,&
    cparm       ,&
    gnu         ,&
    gparm       ,&
    hparm       ,&
    iaggr       ,&
    icloud      ,&
    igraup      ,&
    ihail       ,&
    ipris       ,&
    irain       ,&
    isnow       ,&
    level       ,&
    mkcoltab    ,&
    pparm       ,&
    rparm       ,&
    sparm       ,&
    idriz       ,&
    irime       ,&
    iplaws      ,&
    idust       ,&
    isalt       ,&
    imbudget    ,&
    imbudtot    ,&
    iccnlev     ,&
    dparm,epsil,jnmb,cnparm,gnparm

IMPLICIT NONE

!--------------------------------------------------------------------------
!     The product [(nthz-1)  * dthz ] must equal 25.0.
!     The product [(nrhhz-1) * drhhz] must equal 0.18.
!     The product [(ntc-1)   * dtc  ] must equal 20.0.
!     The product [(ndnc-1)  * ddnc ] must equal 20.e-6.

integer, parameter :: nthz=26,nrhhz=10,ngam=5000,ninc=201   &
                     ,ndns=15,ntc=21,ndnc=11                &
                     ,ndcc=60,ndcd=20,ndccr=15,nrrcr=10     &
                     ,ndrcr=30,nrrr=20,ndrr=20              &
                     ,ncat=8,nhcat=16,npairc=101,npairr=147 &
                     ,nembc=20
real, parameter    :: dtc=1.,ddnc=2.e-6 ,dthz=1.,drhhz=.02
real, parameter    :: budget_scale=1000., budget_scalet=1.
!--------------------------------------------------------------------------
!integer :: mcphys_type,level,icloud,idriz,irain,ipris,isnow,iaggr,igraup,ihail  &
!          ,mkcoltab,irime,iplaws,idust,isalt,imbudget,imbudtot

!integer, dimension(ncat)        :: jnmb
integer, dimension(nhcat,nhcat) :: ipairc,ipairr
integer, dimension(31,100,2)    :: jhabtab
integer, dimension(nzpmax,ncat) :: jhcat,ict1,ict2
integer :: iamsflg,iss1flg,iss2flg,imd1flg,imd2flg
logical, parameter :: lhrtheta = .true.


real :: &! cparm,dparm,rparm,pparm,sparm,aparm,gparm,hparm,epsil         &
         !cnparm,gnparm,
	d1parm,d2parm                             &
       ,rictmin,rictmax,dps,dps2                                &
       ,d1min,d1max,d2min,d2max,d3min,d3minx,d3max,r3min,r3max  &
       ,d1ecr,d2ecr,d3ecr,r3ecr,d3err,r3err                     &
       ,colf,pi4dt,sedtime0,sedtime1                            &
       ,dimin,diminx,dimax,rimin,rimax,dieci,rieci

!real, dimension(ncat) :: emb0,emb1,gnu,parm,emb0log,emb1log,dict
real, dimension(ncat)  :: emb0,emb1    ,parm,emb0log,emb1log,dict
real, dimension(nhcat) :: var_shape,cfmas,pwmas,cfvt,pwvt,dpsmi,cfden,pwden  &
                         ,cfemb0,cfen0,pwemb0,pwen0,vtfac,frefac1,frefac2  &
                         ,cfmasft,dnfac,sipfac,pwmasi,ch1,ch3,cdp1,pwvtmasi
real, dimension(nzpmax) :: tair,tairc,tairstrc,til,rvstr,press,pitot  &
                          ,rliq,rice,qhydm,rvlsair,rvisair,rvs0,thrmcon  &
                          ,vapdif,dynvisc,rdynvsci,denfac,dn0i,colfacr  &
                          ,colfacr2,colfacc,colfacc2,sumuy,sumuz,sumvr  &
                          ,scrmic1,scrmic2,scrmic3,cccnx,gccnx,cifnx &
                          ,cccmx,gccmx,md1nx,md2nx,saltfx,saltjx,saltsx
real, dimension(nzpmax,ncat) :: rx,cx,qr,qx,tx,emb,vterm,vap,ttest,wct1  &
                               ,wct2,sb,sd,se,sf,sg,sh,sm,ss,su,sw,sy,sz &
                               ,rx_lhr,qx_lhr,cnmhx,pcpvx
real, dimension(nzpmax,2)  :: tref,rvsref,rvsrefp
real, dimension(nzpmax,9)  :: sa
real, dimension(nzpmax,12) :: eff

real, dimension(nzpmax,ncat,ncat) :: rxfer,qrxfer,enxfer
real, dimension(nhcat,maxgrds)    :: dispemb0,dispemb1,ch2

real, dimension(nembc,nembc,npairc) :: coltabc
real, dimension(nembc,nembc,npairr) :: coltabr

real, dimension(nrhhz,nthz)       :: frachz
real, dimension(ndnc,ntc,maxgrds) :: fracc
real, dimension(4)                :: gamm,gamn1
real, dimension(ngam,3)           :: gam
real, dimension(ngam,2)           :: gaminc
real, dimension(2,ngam)           :: gamsip13,gamsip24
real, dimension(ninc,nhcat)       :: rmlttab
real, dimension(ninc,nhcat)       :: enmlttab
real, dimension(ninc,ndns)        :: shedtab
real, dimension(2)                :: sc,sk,sl
real, dimension(8)                :: sj,pcprx,accpx

!Lookup tables for droplet autoconversion and liquid collection
real, dimension(ndcc)        :: r1tabcc,c1tabcc,c2tabcc
real, dimension(ndcc)        :: r2tabdd,c2tabdd,c3tabdd
real, dimension(ndcd,ndcd)   :: r1tabcd,c1tabcd,r2tabcd,c2tabcd
real, dimension(ndccr,nrrcr,ndrcr)  :: r1tabcr,c1tabcr,r2tabcr,c2tabcr
real, dimension(nrrr,ndrr)          :: c3tabrr


!Lookup table arrays for Binned Riming Scheme
real, dimension(ndccr,nrrcr,ndrcr,4) :: r1tabci,c1tabci,r2tabci,c2tabci &
                                       ,r1rimer,r2rimer

!******Variables Needed for BUBBLE SIMULATION ************************
integer :: ibubble,ibdxia,ibdxiz,ibdyja,ibdyjz,ibdzk1,ibdzk2
real :: bthp,brtp
!*********************************************************************

!******Variables Needed for CCN nucleation and restore **************
!integer :: iccnlev,ic,rgb
integer :: ic,rgb
real :: rxferratio,ccnmass,ccnnum,rxtemp,fracmass,cxloss
real :: rhosol,rg,ant,rcm,rmlar,rmsma,power,pctccnmas

!Number of bins in lognormal aerosol distribution
integer, parameter :: itbin=100

!Solubility fraction (epsilon)
integer,parameter :: maxeps=7
real, dimension(maxeps) :: epsfrac
data epsfrac / 0.05,0.1,0.2,0.4,0.6,0.8,1.0 /

!Median radii (cm) for CCN
integer, parameter :: maxrg=9
real, dimension(maxrg) :: rg_ccn
data rg_ccn / 0.01e-4,0.02e-4,0.04e-4,0.08e-4 &
             ,0.16e-4,0.32e-4,0.48e-4,0.64e-4 &
             ,0.96e-4 /

real, dimension(itbin,maxrg) :: ccncon,ccnmas
real, dimension(itbin)       :: binrad,smass,rs,con0

!*********************************************************************

!******Variables Needed for Printing Micro Budgets********************
real :: cldnuct,cldnuct2,cld2pris,cld2pris2,prisnuc
real :: raincollc1,raincollc2,c2collc1
real :: sumrim35_1,sumrim45_1
real :: sumrim35_2,sumrim45_2
real :: sumrim35_3,sumrim36_3,sumrim37_3,sumrim45_3,sumrim46_3,sumrim47_3 &
       ,sumrim56_3,sumrim57_3,sumrim67_3
real ::     sumrim141,sumrim142,sumrim143,sumrim151,sumrim152,sumrim153 &
 ,sumrim161,sumrim162,sumrim163,sumrim171,sumrim172,sumrim841,sumrim842 &
 ,sumrim843,sumrim851,sumrim852,sumrim853,sumrim861,sumrim862,sumrim863 &
 ,sumrim871,sumrim872,icesplintC1toPR,icesplintC2toPR,sumrimSPLlow &
 ,sumrimSPLhi
real :: sumrim231,sumrim232,sumrim233,sumrim234,sumrim241,sumrim242 &
       ,sumrim243,sumrim244,sumrim251,sumrim252,sumrim253,sumrim254 &
       ,sumrim261,sumrim262,sumrim263,sumrim264,sumrim271,sumrim272 &
       ,sumrim273,sumrim274
real, dimension(8) :: rgainlcat,rlossjcat,rxferorig
real :: cldevap1,cldevap2,cldevap3,cldevap4,cldevap5,cldevap6,cldevap7 &
     ,cldevap8,cldvap1,cldvap2,cldvap3,cldvap4,cldvap5,cldvap6,cldvap7 &
     ,cldvap8,psrxfer
real :: rmelt46,ricetor46,rmelt56,ricetor56,rmelt6,rmelt7,rmelt3
double precision, dimension(8) :: rxolder,rxnewer
real :: mixingC1,mixingRA,mixingPR,mixingSN,mixingAG,mixingGR,mixingHA &
    ,mixingDR,mixing1C1,mixing1RA,mixing1PR,mixing1SN,mixing1AG        &
    ,mixing1GR,mixing1HA,mixing1DR,numberC1,numberRA,numberPR,numberSN &
    ,numberAG,numberGR,numberHA,numberDR,number1C1,number1RA           &
    ,number1PR,number1SN,number1AG,number1GR,number1HA,number1DR       &
    ,mixingVP,mixingRT,mixingSL,mixingSI                               &
    ,accumuRA,accumuPR,accumuSN,accumuAG,accumuGR,accumuHA,accumuDR    &
    ,mixing1VP,mixing1RT,mixing1SL,mixing1SI


PRIVATE ::                  &
        copyback            &
       ,calc_lhr_vap	    &
       ,calc_lhr_collmelt   &
       ,calc_lhr_cldnuc     &
       ,calc_lhr_icenuc     &
       ,aero_nuc_tab	    &
       ,cldnuc 	            &
       ,icenuc 	            &
       ,contnuc             &
       ,gcf		    &
       ,gser		    &
       ,avint		    &
       ,getict 	            &
       ,auto_accret	    &
       ,auto_accret_ice     &
       ,effxy		    &
       ,cols                &
       ,col3344	            &
       ,col3443	            &
       ,col1		    &
       ,col2		    &
       ,col3		    &
       ,colxfers	    &
       ,each_call	    &
       ,range_check         &
       ,each_column	    &
       ,enemb		    &
       ,x02		    &
       ,c03		    &
       ,pc03		    &
       ,sedim		    &
       ,adj1_2M_rams60      &
       ,ae1mic 	            &
       ,ae1kmic	            &
       ,micinit	            &
       ,haznuc 	            &
       ,homfrzcl	    &
       ,mksedim_tab	    &
       ,tabmelt	            &
       ,mkcoltb	            &
       ,make_autotab	    &
       ,sxy                 &
       ,mcphys_data	    &
       ,data_cs	            &
       ,data_ca	            &
       ,data_cg	            &
       ,data_ch	            &
       ,initg2mode	    &
       ,sumn		    &
       ,initg1mode	    &
       ,tabhab 	            &
       ,thrmstr	            &
       ,vapdiff 	    &
       ,vapflux	            &
       ,newtemp             &
       ,gammp               &
       ,gammq		    &
       ,gammln		    &
       ,xj


PUBLIC :: micro_2M_rams60,negadj1_2M_rams60 &
       ,jnmbinit            &
       ,initqin	            &
       ,initqin2	    &
       ,initqin3	    &
       ,initqin4	    &
       ,initqin5	    &
       ,micro_master


CONTAINS

!###########################################################################

subroutine micro_2M_rams60()

  use mem_basic, only : &
       basic_g            ! INTENT(INOUT)

  use mem_micro, only:  &
       micro_g            ! INTENT(INOUT)

  use mem_grid, only:   &
       ngrids,          & ! INTENT(IN)
       ngrid,           & ! INTENT(IN)
       zm,              & ! INTENT(IN)
       dzt,             & ! INTENT(IN)
       dtlt,            & ! INTENT(IN)
       jdim,            & ! INTENT(IN)
       maxnzp,          & ! INTENT(IN)
       time,            & ! INTENT(IN)
       zt,              & ! INTENT(IN)
       itime1,          & ! INTENT(IN)
       if_adap,         & ! INTENT(IN)
       grid_g,          & ! INTENT(IN)
       nnzp,npatch,imonth1

!  use mem_radiate, only:&
!       radiate_g          ! INTENT(INOUT)

  use node_mod, only :  &
       mzp,             & ! INTENT(IN)
       mxp,             & ! INTENT(IN)
       myp,             & ! INTENT(IN)
       ja,              & ! INTENT(IN)
       jz,              & ! INTENT(IN)
       ia,              & ! INTENT(IN)
       iz,              & ! INTENT(IN)
       mynum              ! INTENT(IN)


use mem_leaf, only : leaf_g
!use rrad3

use rconstants


implicit none

integer :: nembfall,maxkfall,ngr,lhcat,i,j,k
integer, dimension(11)  :: k1,k2,k3
integer, save :: ncall = 0
integer, save, dimension(16)  :: ncall2g
real :: dtlti,albedt,cosz,rlongup,rshort,rlong
data ncall2g/16*0/

type pcp_tab_type
  real, pointer, dimension(:,:,:,:) :: pcpfillc,pcpfillr
  real, pointer, dimension(:,:,:)   :: sfcpcp
  real, pointer, dimension(:,:,:)   :: allpcp
end type

type (pcp_tab_type), save :: pcp_tab(maxgrds)

if (level .ne. 3) return

nembfall = 20
maxkfall = 4

if(ncall == 0) then
   ncall = 1

   do ngr = 1,ngrids
      allocate (pcp_tab(ngr)%pcpfillc(nnzp(ngr),maxkfall,nembfall,nhcat))
      allocate (pcp_tab(ngr)%pcpfillr(nnzp(ngr),maxkfall,nembfall,nhcat))
      allocate (pcp_tab(ngr)%sfcpcp(maxkfall,nembfall,nhcat))
      allocate (pcp_tab(ngr)%allpcp(nnzp(ngr),nembfall,nhcat))
   enddo

   call micinit()
   call make_autotab()
   call haznuc()
   call tabmelt()
   call tabhab()

   do lhcat = 1,nhcat
      ch3(lhcat) = pwvt(lhcat) * pwmasi(lhcat)
      cdp1(lhcat) = pwmasi(lhcat) * (1.5 + .5 * pwvt(lhcat))
      pwvtmasi(lhcat) = pwvt(lhcat) * pwmasi(lhcat)
   enddo
endif

!-srf 2015------------------ not included in BRAMS
!Seigel (9-23-10)
!Run the SEASALT and DUST Source model before the call to Micro
!if(idust==1) then
!  call dust_sources(mzp,mxp,myp,ia,iz,ja,jz                  &
!    ,leaf_g(ngrid)%leaf_class  ,grid_g(ngrid)%rtgm           &
!    ,leaf_g(ngrid)%patch_area                                &
!    ,basic_g(ngrid)%up         ,basic_g(ngrid)%vp            &
!    ,micro_g(ngrid)%cccnp      ,micro_g(ngrid)%gccnp         &
!    ,micro_g(ngrid)%cccmp      ,micro_g(ngrid)%gccmp         &
!    ,micro_g(ngrid)%md1np      ,micro_g(ngrid)%md2np         &
!    ,leaf_g(ngrid)%soil_water  ,leaf_g(ngrid)%soil_text)
!endif
!if(isalt==1) then
!  call salt_sources(mzp,mxp,myp,ia,iz,ja,jz                  &
!    ,leaf_g(ngrid)%leaf_class  ,grid_g(ngrid)%rtgm           &
!    ,leaf_g(ngrid)%patch_area                                &
!    ,basic_g(ngrid)%up         ,basic_g(ngrid)%vp            &
!    ,micro_g(ngrid)%cccnp      ,micro_g(ngrid)%gccnp         &
!    ,micro_g(ngrid)%cccmp      ,micro_g(ngrid)%gccmp         &
!    ,micro_g(ngrid)%salt_filmp ,micro_g(ngrid)%salt_jetp     &
!    ,micro_g(ngrid)%salt_spmp)
!endif
!-srf 2015------------------ not included in BRAMS

if (ncall2g(ngrid) .ne. 5) then
   ncall2g(ngrid) = 5

   call mksedim_tab(mzp,mxp,myp,ngrid,nembfall,maxkfall,zm,dzt  &
      ,pcp_tab(ngrid)%pcpfillc(1:mzp,1,1,1),pcp_tab(ngrid)%pcpfillr(1:mzp,1,1,1)  &
	  ,pcp_tab(ngrid)%sfcpcp(1:maxkfall,1,1),pcp_tab(ngrid)%allpcp(1:mzp,1,1),dtlt)

   do lhcat = 1,nhcat
      ch2(lhcat,ngrid) = float(nembfall-1) &
                     / log10(dispemb1(lhcat,ngrid) / dispemb0(lhcat,ngrid))
   enddo

   call homfrzcl(dtlt,ngrid)
endif

call each_call(mzp,dtlt)
dtlti = 1. / dtlt
ngr = ngrid

do j = ja,jz
   do i = ia,iz

      !Saleeby: add checking here instead of just in Harrington radiation
      do k = 1,mzp
        if( basic_g(ngr)%rv(k,i,j)<0.0 .or.  &
             basic_g(ngr)%rtp(k,i,j)<0.0 .or. &
             (jnmb(1)>0 .and. micro_g(ngr)%rcp(k,i,j)<0.0) .or. &
             (jnmb(2)>0 .and. micro_g(ngr)%rrp(k,i,j)<0.0) .or. &
             (jnmb(3)>0 .and. micro_g(ngr)%rpp(k,i,j)<0.0) .or. &
             (jnmb(4)>0 .and. micro_g(ngr)%rsp(k,i,j)<0.0) .or. &
             (jnmb(5)>0 .and. micro_g(ngr)%rap(k,i,j)<0.0) .or. &
             (jnmb(6)>0 .and. micro_g(ngr)%rgp(k,i,j)<0.0) .or. &
             (jnmb(7)>0 .and. micro_g(ngr)%rhp(k,i,j)<0.0) ) then
             print*,'Negative Condensate MICRO (ngr,k,i,j):',ngr,k,i,j
             print*,'vapor:  ',basic_g(ngr)%rv(k,i,j)
             print*,'rtp:    ',basic_g(ngr)%rtp(k,i,j)
             if(jnmb(1)>0) print*,'cloud:  ',micro_g(ngr)%rcp(k,i,j)
             if(jnmb(2)>0) print*,'rain:   ',micro_g(ngr)%rrp(k,i,j)
             if(jnmb(3)>0) print*,'ice:    ',micro_g(ngr)%rpp(k,i,j)
             if(jnmb(4)>0) print*,'snow:   ',micro_g(ngr)%rsp(k,i,j)
             if(jnmb(5)>0) print*,'aggr:   ',micro_g(ngr)%rap(k,i,j)
             if(jnmb(6)>0) print*,'graup:  ',micro_g(ngr)%rgp(k,i,j)
             if(jnmb(7)>0) print*,'hail:   ',micro_g(ngr)%rhp(k,i,j)
             stop
        endif
      enddo

      call range_check(mzp,k1,k2,k3,i,j,grid_g(ngr)%lpw(i,j),micro_g(ngr))


      call mcphys(&
          mzp&
	 ,k1&
	 ,k2&
	 ,k3&
	 ,i&
	 ,j&
	 ,ngrid&
	 ,jdim&
	 ,maxnzp&
	 ,imonth1&
	 ,npatch     &
         ,nembfall&
	 ,maxkfall&
	 ,mynum&
	 ,dtlt&
	 ,dtlti&
	 ,time&
	 ,zm&
	 ,dzt    &
         ,zt&
	 ,itime1&
	 ,basic_g(ngr)%thp     (1:mzp,i,j)   &
	 ,basic_g(ngr)%theta   (1:mzp,i,j)   &
         ,basic_g(ngr)%pp      (1:mzp,i,j)   &
	 ,basic_g(ngr)%rtp     (1:mzp,i,j)   &
         ,basic_g(ngr)%rv      (1:mzp,i,j)   &
	 ,basic_g(ngr)%wp      (1:mzp,i,j)   &
         ,basic_g(ngr)%dn0     (1:mzp,i,j)   &
	 ,basic_g(ngr)%pi0     (1:mzp,i,j)   &
	 ,grid_g(ngr)%rtgt     (i,j)     &
	 ,grid_g(ngr)%lpw      (i,j)     &
         ,micro_g(ngr)%pcpg    (i,j)     &
	 ,micro_g(ngr)%qpcpg   (i,j)     & !31
         ,micro_g(ngr)%dpcpg   (i,j)     &
	 ,pcp_tab(ngr)%pcpfillc(1:mzp,1,1,1) &
	 ,pcp_tab(ngr)%pcpfillr(1:mzp,1,1,1) &
	 ,pcp_tab(ngr)%sfcpcp  (1:maxkfall,1,1)   &
         ,pcp_tab(ngr)%allpcp  (1:mzp,1,1)   &
	 ,grid_g(ngr)%glat     (i,j)  &
	 ,grid_g(ngr)%topt     (i,j)  &
	 ,if_adap &!40
         ,leaf_g(ngr)%leaf_class(i,j,1:npatch) &
         ,leaf_g(ngr)%ustar(i,j,1:npatch)      &
         ,leaf_g(ngr)%patch_area(i,j,1:npatch) &
         ,leaf_g(ngr)%veg_rough(i,j,1:npatch)  &
         ,leaf_g(ngr)%soil_rough(i,j,1:npatch) &
         ,basic_g(ngr)%up(1:mzp,i,j)               &
         ,basic_g(ngr)%vp(1:mzp,i,j)               &
     !---------------------------------------------------------
      ,micro_g(ngr)%md1np(:,i,j)	    &
      ,micro_g(ngr)%md2np(:,i,j)	    &
      ,micro_g(ngr)%salt_filmp(:,i,j)	    &
      ,micro_g(ngr)%salt_jetp(:,i,j)	    &
      ,micro_g(ngr)%salt_spmp(:,i,j)	    &
      !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES
      ,micro_g(ngr)%nuccldr (:,i,j)      ,micro_g(ngr)%nuccldc (:,i,j)      &
      ,micro_g(ngr)%nucicer (:,i,j)      ,micro_g(ngr)%nucicec (:,i,j)      &
      ,micro_g(ngr)%inuchomr (:,i,j)     ,micro_g(ngr)%inuchomc (:,i,j)     &
      ,micro_g(ngr)%inuccontr (:,i,j)    ,micro_g(ngr)%inuccontc (:,i,j)    &
      ,micro_g(ngr)%inucifnr (:,i,j)     ,micro_g(ngr)%inucifnc (:,i,j)     &
      ,micro_g(ngr)%inuchazr (:,i,j)     ,micro_g(ngr)%inuchazc (:,i,j)     &
      ,micro_g(ngr)%vapliq (:,i,j)       ,micro_g(ngr)%vapice (:,i,j)       &
      ,micro_g(ngr)%vapcld (:,i,j)       ,micro_g(ngr)%vaprain (:,i,j)      &
      ,micro_g(ngr)%vappris (:,i,j)      ,micro_g(ngr)%vapsnow (:,i,j)      &
      ,micro_g(ngr)%vapaggr (:,i,j)      ,micro_g(ngr)%vapgrau (:,i,j)      &
      ,micro_g(ngr)%vaphail (:,i,j)      ,micro_g(ngr)%vapdriz (:,i,j)      &
      ,micro_g(ngr)%meltice (:,i,j)      ,micro_g(ngr)%meltpris (:,i,j)     &
      ,micro_g(ngr)%meltsnow (:,i,j)     ,micro_g(ngr)%meltaggr (:,i,j)     &
      ,micro_g(ngr)%meltgrau (:,i,j)     ,micro_g(ngr)%melthail (:,i,j)     &
      ,micro_g(ngr)%cld2rain (:,i,j)     ,micro_g(ngr)%rimecld (:,i,j)      &
      ,micro_g(ngr)%rimecldsnow (:,i,j)  ,micro_g(ngr)%rimecldaggr (:,i,j)  &
      ,micro_g(ngr)%rimecldgrau (:,i,j)  ,micro_g(ngr)%rimecldhail (:,i,j)  &
      ,micro_g(ngr)%rain2ice (:,i,j)     ,micro_g(ngr)%rain2pr (:,i,j)      &
      ,micro_g(ngr)%rain2sn (:,i,j)      ,micro_g(ngr)%rain2ag (:,i,j)      &
      ,micro_g(ngr)%rain2gr (:,i,j)      ,micro_g(ngr)%rain2ha (:,i,j)      &
      ,micro_g(ngr)%rain2ha_xtra (:,i,j) ,micro_g(ngr)%ice2rain (:,i,j)     &
      ,micro_g(ngr)%aggregate (:,i,j)    ,micro_g(ngr)%aggrselfpris (:,i,j) &
      ,micro_g(ngr)%aggrselfsnow (:,i,j) ,micro_g(ngr)%aggrprissnow (:,i,j) &
      ,micro_g(ngr)%latheatvap (:,i,j)   ,micro_g(ngr)%latheatfrz (:,i,j)   &
      !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES (totals)
      ,micro_g(ngr)%nuccldrt (:,i,j)      ,micro_g(ngr)%nuccldct (:,i,j)      &
      ,micro_g(ngr)%nucicert (:,i,j)      ,micro_g(ngr)%nucicect (:,i,j)      &
      ,micro_g(ngr)%inuchomrt (:,i,j)     ,micro_g(ngr)%inuchomct (:,i,j)     &
      ,micro_g(ngr)%inuccontrt (:,i,j)    ,micro_g(ngr)%inuccontct (:,i,j)    &
      ,micro_g(ngr)%inucifnrt (:,i,j)     ,micro_g(ngr)%inucifnct (:,i,j)     &
      ,micro_g(ngr)%inuchazrt (:,i,j)     ,micro_g(ngr)%inuchazct (:,i,j)     &
      ,micro_g(ngr)%vapliqt (:,i,j)       ,micro_g(ngr)%vapicet (:,i,j)       &
      ,micro_g(ngr)%vapcldt (:,i,j)       ,micro_g(ngr)%vapraint (:,i,j)      &
      ,micro_g(ngr)%vapprist (:,i,j)      ,micro_g(ngr)%vapsnowt (:,i,j)      &
      ,micro_g(ngr)%vapaggrt (:,i,j)      ,micro_g(ngr)%vapgraut (:,i,j)      &
      ,micro_g(ngr)%vaphailt (:,i,j)      ,micro_g(ngr)%vapdrizt (:,i,j)      &
      ,micro_g(ngr)%melticet (:,i,j)      ,micro_g(ngr)%meltprist (:,i,j)     &
      ,micro_g(ngr)%meltsnowt (:,i,j)     ,micro_g(ngr)%meltaggrt (:,i,j)     &
      ,micro_g(ngr)%meltgraut (:,i,j)     ,micro_g(ngr)%melthailt (:,i,j)     &
      ,micro_g(ngr)%cld2raint (:,i,j)     ,micro_g(ngr)%rimecldt (:,i,j)      &
      ,micro_g(ngr)%rimecldsnowt (:,i,j)  ,micro_g(ngr)%rimecldaggrt (:,i,j)  &
      ,micro_g(ngr)%rimecldgraut (:,i,j)  ,micro_g(ngr)%rimecldhailt (:,i,j)  &
      ,micro_g(ngr)%rain2icet (:,i,j)     ,micro_g(ngr)%rain2prt (:,i,j)      &
      ,micro_g(ngr)%rain2snt (:,i,j)      ,micro_g(ngr)%rain2agt (:,i,j)      &
      ,micro_g(ngr)%rain2grt (:,i,j)      ,micro_g(ngr)%rain2hat (:,i,j)      &
      ,micro_g(ngr)%rain2ha_xtrat (:,i,j) ,micro_g(ngr)%ice2raint (:,i,j)     &
      ,micro_g(ngr)%aggregatet (:,i,j)    ,micro_g(ngr)%aggrselfprist (:,i,j) &
      ,micro_g(ngr)%aggrselfsnowt (:,i,j) ,micro_g(ngr)%aggrprissnowt (:,i,j) &
      ,micro_g(ngr)%latheatvapt (:,i,j)   ,micro_g(ngr)%latheatfrzt (:,i,j))

      call copyback(mzp,k1,k2,k3,grid_g(ngr)%lpw(i,j),i,j,micro_g(ngr))

   enddo
enddo

return
end subroutine micro_2M_rams60

!     *****************************************************************

subroutine mcphys(m1,k1,k2,k3,i,j,ngr,jdim,maxnzp ,imonthx,npatch   &
   ,nembfall,maxkfall,mynum,dtlt,dtlti,time,zm,dzt,zt  &
   ,itime1  &!radiate  &
   ,thp,theta,pp,rtp,rv,wp,dn0,pi0  &
   ,rtgt,lpw_R,pcpg,qpcpg,dpcpg  &
   ,pcpfillc,pcpfillr,sfcpcp,allpcp  &
   ,glat,topt,if_adap &
   ,lclass,ustar,parea,vrough,srough &
   ,uup,vvp,md1np,md2np,salt_filmp,salt_jetp,salt_spmp &
   !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES
   ,nuccldr,nuccldc,nucicer,nucicec,inuchomr                   &
   ,inuchomc,inuccontr,inuccontc,inucifnr,inucifnc             &
   ,inuchazr,inuchazc,vapliq,vapice,vapcld                     &
   ,vaprain,vappris,vapsnow,vapaggr,vapgrau                    &
   ,vaphail,vapdriz,meltice,meltpris,meltsnow                  &
   ,meltaggr,meltgrau,melthail,cld2rain,rimecld                &
   ,rimecldsnow,rimecldaggr,rimecldgrau,rimecldhail,rain2ice   &
   ,rain2pr,rain2sn,rain2ag,rain2gr,rain2ha                    &
   ,rain2ha_xtra,ice2rain,aggregate,aggrselfpris,aggrselfsnow  &
   ,aggrprissnow,latheatvap,latheatfrz                         &
   !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES (totals)
   ,nuccldrt,nuccldct,nucicert,nucicect,inuchomrt                   &
   ,inuchomct,inuccontrt,inuccontct,inucifnrt,inucifnct             &
   ,inuchazrt,inuchazct,vapliqt,vapicet,vapcldt                     &
   ,vapraint,vapprist,vapsnowt,vapaggrt,vapgraut                    &
   ,vaphailt,vapdrizt,melticet,meltprist,meltsnowt                  &
   ,meltaggrt,meltgraut,melthailt,cld2raint,rimecldt                &
   ,rimecldsnowt,rimecldaggrt,rimecldgraut,rimecldhailt,rain2icet   &
   ,rain2prt,rain2snt,rain2agt,rain2grt,rain2hat                    &
   ,rain2ha_xtrat,ice2raint,aggregatet,aggrselfprist,aggrselfsnowt  &
   ,aggrprissnowt,latheatvapt,latheatfrzt)

!-srf - not using radiation type 3
!use mem_radiate
!use rconstants
!use rrad3
!!use micphys
!-srf
implicit none

!type (radiate_vars) :: radiate

integer :: i,j,k,lcat,jcat,icv,icx,mc1,mc2,mc3,mc4,m1,lpw,if_adap  &
          ,ngr,jdim,nembfall,maxkfall  &
          ,mynum,maxnzp,mcat  &
          ,k1cnuc,k2cnuc,k1dnuc,k2dnuc,k1pnuc,k2pnuc,lhcat,itime1
real :: lpw_r

real,    dimension(8)   :: dpcp0
integer, dimension(8)   :: mcats,mivap,mix02
integer, dimension(9,4) :: mcat1
integer, dimension(8,2) :: mcat2,mcat22
integer, dimension(4)   :: mcat33
integer, dimension(11)  :: k1,k2,k3

real                       :: dtlt,dtlti,time

real                       :: rtgt,pcpg,qpcpg,dpcpg

real, dimension(m1)        :: zm,dzt,thp,theta,pp,rtp,rv,wp,dn0,pi0
!COMPUTE AND OUTPUT MICRO BUDGET PROCESSES
real, dimension(m1)        :: &
            nuccldr,nuccldc,nucicer,nucicec,inuchomr                  &
           ,inuchomc,inuccontr,inuccontc,inucifnr,inucifnc            &
           ,inuchazr,inuchazc,vapliq,vapice,vapcld                    &
           ,vaprain,vappris,vapsnow,vapaggr,vapgrau                   &
           ,vaphail,vapdriz,meltice,meltpris,meltsnow                 &
           ,meltaggr,meltgrau,melthail,cld2rain,rimecld               &
           ,rimecldsnow,rimecldaggr,rimecldgrau,rimecldhail,rain2ice  &
           ,rain2pr,rain2sn,rain2ag,rain2gr,rain2ha                   &
           ,rain2ha_xtra,ice2rain,aggregate,aggrselfpris,aggrselfsnow &
           ,aggrprissnow,latheatvap,latheatfrz
!COMPUTE AND OUTPUT MICRO BUDGET PROCESSES (totals)
real, dimension(m1)        :: &
            nuccldrt,nuccldct,nucicert,nucicect,inuchomrt                  &
           ,inuchomct,inuccontrt,inuccontct,inucifnrt,inucifnct            &
           ,inuchazrt,inuchazct,vapliqt,vapicet,vapcldt                    &
           ,vapraint,vapprist,vapsnowt,vapaggrt,vapgraut                   &
           ,vaphailt,vapdrizt,melticet,meltprist,meltsnowt                 &
           ,meltaggrt,meltgraut,melthailt,cld2raint,rimecldt               &
           ,rimecldsnowt,rimecldaggrt,rimecldgraut,rimecldhailt,rain2icet  &
           ,rain2prt,rain2snt,rain2agt,rain2grt,rain2hat                   &
           ,rain2ha_xtrat,ice2raint,aggregatet,aggrselfprist,aggrselfsnowt &
           ,aggrprissnowt,latheatvapt,latheatfrzt

! Variables needed for Harrington radiation scheme

real                       :: glat,topt
real, dimension(m1)        :: zt

real, dimension(m1,maxkfall,nembfall,nhcat) :: pcpfillc,pcpfillr
real, dimension(maxkfall,nembfall,nhcat) :: sfcpcp
real, dimension(m1,nembfall,nhcat) :: allpcp

!Saleeby(2011-04-20)
!Variables needed for aerosol sedimentation
integer :: npatch,imonthx,nnzp
real, dimension(npatch) :: lclass,ustar,parea,vrough,srough
real, dimension(m1) :: uup,vvp,md1np,md2np,salt_filmp,salt_jetp,salt_spmp

! (mcats) is for rain, agg, graupel, hail self-collection (cols)
data mcats /0,3,0,0,6,7,10,0/     ! (effxy) number for coll efficiency
! (mcat1) is for ice-ice interactions (col1)
data mcat1 /3,3,3,4,4,4,5,5,6  &  ! 1st variable to interact
           ,5,6,7,5,6,7,6,7,7  &  ! 2nd variable to interact
           ,5,6,7,5,6,7,6,7,7  &  ! transfer variable for (r) & (q)
           ,4,7,8,5,7,8,7,8,8/    ! (effxy) number for coll efficiency
! (mcat2) is for cloud_1 - ice interactions (col2)
data mcat2 /0,0,0,6,6,7,7,0  &
           ,0,0,0,2,2,9,9,0/
! (mcat22) is for cloud_2 - ice interactions (col2)
data mcat22 /0,0,0,6,6,7,7,0  &
            ,0,0,0,11,11,12,12,0/
! (mcat33) is for pristine, snow self-collection (col3344)
data mcat33 /0,0,4,5/  ! (effxy) number for coll efficiency
! Order of vapor flux calculation
data mivap /1,8,3,4,5,2,6,7/
! Order of melting calculation and transfer
! though Cloud and Drizzle not considered
data mix02 /3,1,8,4,5,6,7,2/
! Multiplier in (sedim)
data dpcp0 /.001,.001,.010,.010,.010,.003,.001,.001/
save

!ZERO OUT MICRO BUDGET PROCESSES
if(imbudget>=1 .or. imbudtot>=1) then
 do k = 1,m1
  latheatvap(k) = 0.
  latheatfrz(k) = 0.
 enddo
endif
if(imbudget>=1) then
 do k = 1,m1
  nuccldr(k) = 0.
  nuccldc(k) = 0.
  cld2rain(k) = 0.
  ice2rain(k) = 0.
  nucicer(k) = 0.
  nucicec(k) = 0.
  vapliq(k) = 0.
  vapice(k) = 0.
  meltice(k) = 0.
  rimecld(k) = 0.
  aggregate(k) = 0.
  rain2ice(k) = 0.
 enddo
endif
if(imbudget==2) then
 do k = 1,m1
  inuchomr(k) = 0.
  inuchomc(k) = 0.
  inuccontr(k) = 0.
  inuccontc(k) = 0.
  inucifnr(k) = 0.
  inucifnc(k) = 0.
  inuchazr(k) = 0.
  inuchazc(k) = 0.
  vapcld(k) = 0.
  vaprain(k) = 0.
  vappris(k) = 0.
  vapsnow(k) = 0.
  vapaggr(k) = 0.
  vapgrau(k) = 0.
  vaphail(k) = 0.
  vapdriz(k) = 0.
  meltpris(k) = 0.
  meltsnow(k) = 0.
  meltaggr(k) = 0.
  meltgrau(k) = 0.
  melthail(k) = 0.
  rimecldsnow(k) = 0.
  rimecldaggr(k) = 0.
  rimecldgrau(k) = 0.
  rimecldhail(k) = 0.
  rain2pr(k) = 0.
  rain2sn(k) = 0.
  rain2ag(k) = 0.
  rain2gr(k) = 0.
  rain2ha(k) = 0.
  rain2ha_xtra(k) = 0.
  aggrselfpris(k) = 0.
  aggrselfsnow(k) = 0.
  aggrprissnow(k) = 0.
 enddo
endif

!ZERO OUT MICRO BUDGET PROCESSES (totals)
if(time .lt. .001) then
if(imbudtot>=1) then
 do k = 1,m1
  nuccldrt(k) = 0.
  nuccldct(k) = 0.
  cld2raint(k) = 0.
  ice2raint(k) = 0.
  nucicert(k) = 0.
  nucicect(k) = 0.
  vapliqt(k) = 0.
  vapicet(k) = 0.
  melticet(k) = 0.
  rimecldt(k) = 0.
  aggregatet(k) = 0.
  rain2icet(k) = 0.
  latheatvapt(k) = 0.
  latheatfrzt(k) = 0.
 enddo
endif
if(imbudtot==2) then
 do k = 1,m1
  inuchomrt(k) = 0.
  inuchomct(k) = 0.
  inuccontrt(k) = 0.
  inuccontct(k) = 0.
  inucifnrt(k) = 0.
  inucifnct(k) = 0.
  inuchazrt(k) = 0.
  inuchazct(k) = 0.
  vapcldt(k) = 0.
  vapraint(k) = 0.
  vapprist(k) = 0.
  vapsnowt(k) = 0.
  vapaggrt(k) = 0.
  vapgraut(k) = 0.
  vaphailt(k) = 0.
  vapdrizt(k) = 0.
  meltprist(k) = 0.
  meltsnowt(k) = 0.
  meltaggrt(k) = 0.
  meltgraut(k) = 0.
  melthailt(k) = 0.
  rimecldsnowt(k) = 0.
  rimecldaggrt(k) = 0.
  rimecldgraut(k) = 0.
  rimecldhailt(k) = 0.
  rain2prt(k) = 0.
  rain2snt(k) = 0.
  rain2agt(k) = 0.
  rain2grt(k) = 0.
  rain2hat(k) = 0.
  rain2ha_xtrat(k) = 0.
  aggrselfprist(k) = 0.
  aggrselfsnowt(k) = 0.
  aggrprissnowt(k) = 0.
 enddo
endif
endif

 call thrmstr(m1,k1,k2,lpw,pp(1),thp(1),theta(1),pi0(1),rtp(1),rv(1),i,j)

 call each_column(m1,k1,k2,i,j,lpw,rv(1),dn0(1))

! Diagnose hydrometeor mean mass emb, and if necessary, number concentration.

do lcat = 1,8
   if (jnmb(lcat) .ge. 1) then
      call enemb(m1,k1(lcat),k2(lcat),lcat,dn0(1),i,j)
   endif
enddo

! Evaluate radiative heating rates if using Harrington radiation scheme
!
!if (iswrtyp .eq. 3 .or. ilwrtyp .eq. 3) then
!   if (mod(time + .001,radfrq) .lt. dtlt .or. time .lt. .001) then
!      !Saleeby(2008): Change passing of 7 to 8 if adding drizzle mode
!      ! and modify locations in radcalc3 and radcomp3 to match
!
!      call radcalc3(m1,maxnzp,7,iswrtyp,ilwrtyp,if_adap,lpw  &
!	     ,glat,rtgt,topt  &
!         ,radiate%albedt  (i,j) ,radiate%cosz  (i,j)  &
!         ,radiate%rlongup (i,j) ,radiate%rshort(i,j)  &
!         ,radiate%rlong   (i,j)  &
!        ,zm,zt,rv(1),dn0(1),radiate%fthrd(1,i,j),i,j,time,ngr &
!         !Saleeby(2011): Variables for the radiatively active aerosol
!         ,radiate%bext(1,i,j),radiate%swup(1,i,j),radiate%swdn(1,i,j) &
!         ,radiate%lwup(1,i,j),radiate%lwdn(1,i,j))
!
!   endif
!endif

! Save rx and qx before vapor diffusion...
rx_lhr = rx
qx_lhr = qx

do lcat = 1,8
   if (jnmb(lcat) .ge. 1) then
      call diffprep(m1,lcat,k1(lcat),k2(lcat),rv(1),dn0(1),i,j,mynum)
   endif
enddo

! Pre-compute variables for vapor diffusion
 call vapdiff(m1,k1(11),k2(11),rv(1),i,j,mynum)

!***********************************************************************
! Calculate vapor flux in order of species given by (mivap)
do icv = 1,8
   lcat = mivap(icv)
   if (jnmb(lcat) .ge. 1) then
      call vapflux(m1,lcat,i,j,mynum,k1(lcat),k2(lcat),dn0(1),rv(1) &
               ,vapliq(1),vapice(1),vapcld(1),vaprain(1),vappris(1) &
               ,vapsnow(1),vapaggr(1),vapgrau(1),vaphail(1),vapdriz(1) &
          ,vapliqt(1),vapicet(1),vapcldt(1),vapraint(1),vapprist(1) &
          ,vapsnowt(1),vapaggrt(1),vapgraut(1),vaphailt(1),vapdrizt(1))
   endif
enddo

!**********************************************************************
! Pristine ice to snow transfer
if (jnmb(4) .ge. 1) then
   call psxfer (m1,min(k1(3),k1(4)),max(k2(3),k2(4)),dn0(1),i,j)
endif

do lcat = 1,8
   if (jnmb(lcat) .ge. 1) then
      call enemb(m1,k1(lcat),k2(lcat),lcat,dn0(1),i,j)
      call getict(k1(lcat),k2(lcat),lcat,i,j,mynum)
   endif
enddo

! Update temperature and moisture following vapor growth
 call newtemp(m1,k1(11),k2(11),rv(1),theta(1),i,j)

! Add call to update LHR theta' after vapor diffusion
 call calc_lhr_vap(m1,k1,k2,i,j,latheatvap(1),latheatvapt(1))

!***********************************************************************
! Auto-accretion considered only for rain
if (jnmb(2) .ge. 1) then
   call auto_accret(m1,k1(1),k2(1),k1(8),k2(8),dn0(1),dtlt,i,j &
                   ,cld2rain(1),cld2raint(1))
endif

! Binned riming of cloud species by the ice species
if(irime==1)then
 do jcat = 1,8,7
 do lcat = 4,7
  if (jnmb(jcat) .ge. 1 .and. jnmb(lcat) .ge. 1) then
   call auto_accret_ice(m1,jcat,lcat,k1(jcat),k2(jcat),dn0(1),dtlt,i,j &
     ,rimecld(1),rimecldsnow(1),rimecldaggr(1),rimecldgrau(1),rimecldhail(1) &
     ,rimecldt(1),rimecldsnowt(1),rimecldaggrt(1),rimecldgraut(1) &
     ,rimecldhailt(1))
  endif
 enddo
 enddo
endif

!***********************************************************************
! Calculate collision/coalescence efficiencies before collection subroutines
 call effxy(m1,k1,k2,i,j)

!***********************************************************************
! Self collection of rain, aggregates, graupel, hail:  number change only
do lcat = 2,7
   if (lcat .eq. 3 .or. lcat .eq. 4) go to 29
   mc1 = mcats(lcat)
   if (jnmb(lcat) >= 5) then
      call cols (m1,lcat,mc1,k1(lcat),k2(lcat),i,j)
   endif
29 continue
enddo

! Self collection of pristine ice, snow (transfers to aggregates)
do lcat = 3,4
   mc1 = mcat33(lcat)
   if (jnmb(lcat) .ge. 1 .and. jnmb(5) .ge. 1) then
      call col3344 (m1,lcat,5,mc1,k1(lcat),k2(lcat),i,j,aggregate(1) &
                   ,aggrselfpris(1),aggrselfsnow(1),aggregatet(1) &
                   ,aggrselfprist(1),aggrselfsnowt(1))
   endif
enddo

! Collection between pristine ice and snow
if (jnmb(5) .ge. 1) then
    call col3443 (m1,3,4,5,max(k1(3),k1(4)),min(k2(3),k2(4)),i,j,aggregate(1) &
                 ,aggrprissnow(1),aggregatet(1),aggrprissnowt(1))
endif

! Ice-ice collisions
do icx = 1,9
   mc1 = mcat1(icx,1)
   mc2 = mcat1(icx,2)
   mc3 = mcat1(icx,3)
   mc4 = mcat1(icx,4)
   if (jnmb(mc1) .ge. 1 .and. jnmb(mc3) .ge. 1) then
      call col1 (m1,mc1,mc2,mc3,mc4,max(k1(mc1),k1(mc2))  &
                ,min(k2(mc1),k2(mc2)),i,j)
   endif
enddo

! Ice-cloud collisions
if(irime==0) then
do jcat = 1,8,7
 do lcat = 4,7
   if(jcat==1)then
     mc1=mcat2(lcat,1)
     mc2=mcat2(lcat,2)
   endif
   if(jcat==8)then
     mc1=mcat22(lcat,1)
     mc2=mcat22(lcat,2)
   endif
   if (jnmb(jcat) .ge. 1 .and. jnmb(lcat).ge.1 .and. jnmb(mc1).ge.1) then
      call col2 (m1,jcat,lcat,mc1,mc2 &
    ,max(k1(jcat),k1(lcat)),min(k2(jcat),k2(lcat)),dn0(1),dtlt,i,j &
    ,rimecld(1),rimecldsnow(1),rimecldaggr(1),rimecldgrau(1),rimecldhail(1) &
    ,rimecldt(1),rimecldsnowt(1),rimecldaggrt(1),rimecldgraut(1),rimecldhailt(1))
   endif
 enddo
enddo
endif

! Ice-rain collisions
do lcat = 3,7
   if (jnmb(2) .ge. 1 .and. jnmb(lcat) .ge. 1 .and. jnmb(7) .ge. 1) then
      call col3 (m1,2,lcat,7,max(k1(2),k1(lcat)),min(k2(2),k2(lcat)),i,j &
       ,ice2rain(1),rain2ice(1),rain2pr(1),rain2sn(1)         &
       ,rain2ag(1),rain2gr(1),rain2ha(1),rain2ha_xtra(1)      &
       ,ice2raint(1),rain2icet(1),rain2prt(1),rain2snt(1)     &
       ,rain2agt(1),rain2grt(1),rain2hat(1),rain2ha_xtrat(1))
   endif
enddo

! Save rx and qx before collisions and melting...
! (colxfers is where mixing ratio is changed after collisions)
rx_lhr = rx
qx_lhr = qx

! Make hydrometeor transfers due to collision-coalescence
 call colxfers(m1,k1,k2,i,j,scrmic1,scrmic2)

! Calcs r,q,c for each category considering melting processes
! in the order of pristine,cloud,drizzle,snow,agg,graupel,hail,rain
! though nothing done for cloud or drizzle in (x02) subroutine
do mcat = 1,8
   lcat = mix02(mcat)
   if (jnmb(lcat) .ge. 1) then
      call x02(m1,k1,k2,lcat,dn0(1),i,j,meltice(1),meltpris(1),meltsnow(1) &
              ,meltaggr(1),meltgrau(1),melthail(1),melticet(1),meltprist(1) &
              ,meltsnowt(1),meltaggrt(1),meltgraut(1),melthailt(1))
   endif
enddo

! Add call to update LHR theta' after collisions and melting
 call calc_lhr_collmelt(m1,k1,k2,i,j,latheatfrz(1),latheatfrzt(1))

! Save rx and qx before cloud nucleation...
rx_lhr = rx
qx_lhr = qx

! Calcs cloud mix ratio and concentration for cloud nucleation from CCN
if (jnmb(1) .ge. 1) then
   call cldnuc(m1,k1cnuc,k2cnuc,k1dnuc,k2dnuc,lpw,rv(1),wp(1),i,j,dn0(1) &
        ,nuccldr(1),nuccldc(1),nuccldrt(1),nuccldct(1))
endif

if(jnmb(1) .ge. 1) then
 !Finds bottom and top layer of cloud water
 k1(1) = min(k1(1),k1cnuc)
 k2(1) = max(k2(1),k2cnuc)
 k3(1) = max(k2(1),k3(1))
endif
if(jnmb(8) .ge. 5) then
 !Finds bottom and top layer of drizzle water
 k1(8) = min(k1(8),k1dnuc)
 k2(8) = max(k2(8),k2dnuc)
 k3(8) = max(k2(8),k3(8))
endif

! Add call to update LHR theta' after cloud nucleation...
 call calc_lhr_cldnuc(m1,k1,k2,i,j,latheatvap(1),latheatvapt(1))

! Calcs mass of cloud water
if (jnmb(1) .ge. 1) then
   call c03(m1,k1(1),k2(1),1,dn0(1),i,j)
endif
! Calcs mass of drizzle water
if (jnmb(8) .ge. 5) then
   call c03(m1,k1(8),k2(8),8,dn0(1),i,j)
endif

! Save rx and qx before ice nucleation...
rx_lhr = rx
qx_lhr = qx

if (jnmb(3) .ge. 1) then
   call icenuc(m1,k1(1),k2(1),k1(8),k2(8),k1pnuc,k2pnuc,lpw,ngr,rv(1)  &
   ,dn0(1),dtlt,i,j                                                    &
   ,nucicer(1),nucicec(1),inuchomr(1),inuchomc(1),inuccontr(1)         &
   ,inuccontc(1),inucifnr(1),inucifnc(1),inuchazr(1),inuchazc(1)       &
   ,nucicert(1),nucicect(1),inuchomrt(1),inuchomct(1),inuccontrt(1)    &
   ,inuccontct(1),inucifnrt(1),inucifnct(1),inuchazrt(1),inuchazct(1))
endif

! Finds bottom and top later of pristine ice
k1(3) = min(k1(3),k1pnuc)
k2(3) = max(k2(3),k2pnuc)
k3(3) = max(k2(3),k3(3))

!Update latent heating budgets after ice nucleation
 call calc_lhr_icenuc(m1,k1,k2,i,j,latheatvap(1),latheatvapt(1))

! Calcs mass of pristine, cloud, and drizzle in order and adjusts fall speed
if (jnmb(3) .ge. 1) then
   call pc03(m1,k1(3),k2(3),3,dn0(1),i,j)
endif
if (jnmb(1) .ge. 1) then
   call pc03(m1,k1(1),k2(1),1,dn0(1),i,j)
endif
if (jnmb(8) .ge. 1) then
   call pc03(m1,k1(8),k2(8),8,dn0(1),i,j)
endif

!  Zero out precip arrays.

pcpg  = 0.
qpcpg = 0.
dpcpg = 0.

! tairc is used here to accumulate changes to thp from sedim

do k = lpw,m1
   tairc(k) = 0.
enddo

do lhcat = 2,nhcat
   ch1(lhcat) = dtlt * cfvt(lhcat) / rtgt
enddo

do lcat = 2,8
   if (jnmb(lcat) .ge. 1) then
      call sedim (m1,lcat,ngr,nembfall,maxkfall  &
         ,k1(lcat),k2(lcat),lpw,i,j  &
         ,rtp(1),thp(1),theta(1),dn0(1),dpcp0(lcat)  &
         ,pcpg,qpcpg,dpcpg,dtlti,scrmic1,scrmic2,scrmic3  &
         ,pcpfillc,pcpfillr,sfcpcp,allpcp,dzt,if_adap)
   endif
enddo

do k = lpw,m1
   thp(k) = thp(k) + tairc(k)
enddo


!
!! Dust and sea-salt dry and wet deposition driver
! call deposition_driver(i,j,m1,imonthx,ustar(1),lclass(1) &
!     ,topt,zm,rv(1),pi0(1),pp(1),theta(1),uup(1),vvp(1)          &
!     ,dn0(1),k1(1),parea(1),vrough(1),srough(1)                  &
!     ,md1np(1),md2np(1),salt_filmp(1),salt_jetp(1))
!
return
end subroutine mcphys

!******************************************************************************

subroutine copyback(m1,k1,k2,k3,lpw_R,i,j,micro)

use mem_micro
!srf - !use micphys

implicit none

type (micro_vars) :: micro

integer, dimension(11)  :: k1,k2,k3
real :: lpw_R
integer :: m1,i,j,lpw

integer :: k

lpw=int(lpw_R)

if (jnmb(1) >= 1) then
   call ae1kmic(lpw,k3(1),micro%rcp(1:m1,i,j),rx(1:m1,1))
   if (jnmb(1) >= 5) then
     call ae1kmic(lpw,k3(1),micro%ccp(1:m1,i,j),cx(1:m1,1))
     call ae1kmic(lpw,k3(1),micro%cccnp(1:m1,i,j),cccnx(1:m1))
     call ae1kmic(lpw,k3(1),micro%cccmp(1:m1,i,j),cccmx(1:m1))
     if(iccnlev >= 2) call ae1kmic(lpw,k3(1),micro%cnm1p(1:m1,i,j),cnmhx(1:m1,1))
   endif
endif

if (jnmb(2) >= 1) then
   call ae1kmic(lpw,k2(11),micro%rrp(1:m1,i,j),rx(1:m1,2))
   call ae1kmic(lpw,k2(11),micro%q2(1:m1,i,j),qx(1:m1,2))
   micro%accpr(i,j) = micro%accpr(i,j) + accpx(2)
   micro%pcprr(i,j) = pcprx(2)
   call ae1kmic(lpw,k2(11),micro%pcpvr(1:m1,i,j),pcpvx(1:m1,2))
   if (jnmb(2) >= 5) call ae1kmic(lpw,k2(11),micro%crp(1:m1,i,j),cx(1:m1,2))
   if (iccnlev >= 2) call ae1kmic(lpw,k2(11),micro%cnm2p(1:m1,i,j),cnmhx(1:m1,2))
endif

if (jnmb(3) >= 1) then
   call ae1kmic(lpw,k3(3),micro%rpp(1:m1,i,j),rx(1:m1,3))
   micro%accpp(i,j) = micro%accpp(i,j) + accpx(3)
   micro%pcprp(i,j) = pcprx(3)
   call ae1kmic(lpw,k3(3),micro%pcpvp(1:m1,i,j),pcpvx(1:m1,3))
   if (jnmb(3) >= 5) call ae1kmic(lpw,k3(3),micro%cpp(1:m1,i,j),cx(1:m1,3))
   if (iccnlev >= 2) call ae1kmic(lpw,k3(3),micro%cnm3p(1:m1,i,j),cnmhx(1:m1,3))
endif

if (jnmb(4) >= 1) then
   call ae1kmic(lpw,k2(11),micro%rsp(1:m1,i,j),rx(1:m1,4))
   micro%accps(i,j) = micro%accps(i,j) + accpx(4)
   micro%pcprs(i,j) = pcprx(4)
   call ae1kmic(lpw,k2(11),micro%pcpvs(1:m1,i,j),pcpvx(1:m1,4))
   if (jnmb(4) >= 5) call ae1kmic(lpw,k2(11),micro%csp(1:m1,i,j),cx(1:m1,4))
endif

if (jnmb(5) >= 1) then
   call ae1kmic(lpw,k2(11),micro%rap(1:m1,i,j),rx(1:m1,5))
   micro%accpa(i,j) = micro%accpa(i,j) + accpx(5)
   micro%pcpra(i,j) = pcprx(5)
   call ae1kmic(lpw,k2(11),micro%pcpva(1:m1,i,j),pcpvx(1:m1,5))
   if (jnmb(5) >= 5) call ae1kmic(lpw,k2(11),micro%cap(1:m1,i,j),cx(1:m1,5))
endif

if (jnmb(6) >= 1) then
   call ae1kmic(lpw,k2(11),micro%rgp(1:m1,i,j),rx(1:m1,6))
   call ae1kmic(lpw,k2(11),micro%q6(1:m1,i,j),qx(1:m1,6))
   micro%accpg(i,j) = micro%accpg(i,j) + accpx(6)
   micro%pcprg(i,j) = pcprx(6)
   call ae1kmic(lpw,k2(11),micro%pcpvg(1:m1,i,j),pcpvx(1:m1,6))
   if (jnmb(6) >= 5) call ae1kmic(lpw,k2(11),micro%cgp(1:m1,i,j),cx(1:m1,6))
endif

if (jnmb(7) >= 1) then
   call ae1kmic(lpw,k2(11),micro%rhp(1:m1,i,j),rx(1:m1,7))
   call ae1kmic(lpw,k2(11),micro%q7(1:m1,i,j),qx(1:m1,7))
   micro%accph(i,j) = micro%accph(i,j) + accpx(7)
   micro%pcprh(i,j) = pcprx(7)
   call ae1kmic(lpw,k2(11),micro%pcpvh(1:m1,i,j),pcpvx(1:m1,7))
   if (jnmb(7) >= 5) call ae1kmic(lpw,k2(11),micro%chp(1:m1,i,j),cx(1:m1,7))
endif

if (jnmb(8) >= 1) then
   if(jnmb(8) <= 4) k3(8) = k2(11)
   call ae1kmic(lpw,k3(8),micro%rdp(1:m1,i,j),rx(1:m1,8))
   micro%accpd(i,j) = micro%accpd(i,j) + accpx(8)
   micro%pcprd(i,j) = pcprx(8)
   call ae1kmic(lpw,k3(8),micro%pcpvd(1:m1,i,j),pcpvx(1:m1,8))
   if (jnmb(8) >= 5) then
     call ae1kmic(lpw,k3(8),micro%cdp(1:m1,i,j),cx(1:m1,8))
     call ae1kmic(lpw,k3(8),micro%gccnp(1:m1,i,j),gccnx(1:m1))
     call ae1kmic(lpw,k3(8),micro%gccmp(1:m1,i,j),gccmx(1:m1))
   endif
   if (iccnlev >= 2) call ae1kmic(lpw,k3(8),micro%cnm8p(1:m1,i,j),cnmhx(1:m1,8))
endif

!Set bottom level with first level above ground
if(imbudget>=1 .or. imbudtot>=1) then
  micro%latheatvap(1,i,j) = micro%latheatvap(2,i,j)
  micro%latheatfrz(1,i,j) = micro%latheatfrz(2,i,j)
endif
if(imbudget>=1) then
  micro%nuccldr(1,i,j)   = micro%nuccldr(2,i,j)
  micro%nuccldc(1,i,j)   = micro%nuccldc(2,i,j)
  micro%cld2rain(1,i,j)  = micro%cld2rain(2,i,j)
  micro%ice2rain(1,i,j)  = micro%ice2rain(2,i,j)
  micro%nucicer(1,i,j)   = micro%nucicer(2,i,j)
  micro%nucicec(1,i,j)   = micro%nucicec(2,i,j)
  micro%vapliq(1,i,j)    = micro%vapliq(2,i,j)
  micro%vapice(1,i,j)    = micro%vapice(2,i,j)
  micro%meltice(1,i,j)   = micro%meltice(2,i,j)
  micro%rimecld(1,i,j)   = micro%rimecld(2,i,j)
  micro%rain2ice(1,i,j)  = micro%rain2ice(2,i,j)
  micro%aggregate(1,i,j) = micro%aggregate(2,i,j)
endif
if(imbudget==2) then
  micro%inuchomr(1,i,j)     = micro%inuchomr(2,i,j)
  micro%inuchomc(1,i,j)     = micro%inuchomc(2,i,j)
  micro%inuccontr(1,i,j)    = micro%inuccontr(2,i,j)
  micro%inuccontc(1,i,j)    = micro%inuccontc(2,i,j)
  micro%inucifnr(1,i,j)     = micro%inucifnr(2,i,j)
  micro%inucifnc(1,i,j)     = micro%inucifnc(2,i,j)
  micro%inuchazr(1,i,j)     = micro%inuchazr(2,i,j)
  micro%inuchazc(1,i,j)     = micro%inuchazc(2,i,j)
  micro%vapcld(1,i,j)       = micro%vapcld(2,i,j)
  micro%vaprain(1,i,j)      = micro%vaprain(2,i,j)
  micro%vappris(1,i,j)      = micro%vappris(2,i,j)
  micro%vapsnow(1,i,j)      = micro%vapsnow(2,i,j)
  micro%vapaggr(1,i,j)      = micro%vapaggr(2,i,j)
  micro%vapgrau(1,i,j)      = micro%vapgrau(2,i,j)
  micro%vaphail(1,i,j)      = micro%vaphail(2,i,j)
  micro%vapdriz(1,i,j)      = micro%vapdriz(2,i,j)
  micro%meltpris(1,i,j)     = micro%meltpris(2,i,j)
  micro%meltsnow(1,i,j)     = micro%meltsnow(2,i,j)
  micro%meltaggr(1,i,j)     = micro%meltaggr(2,i,j)
  micro%meltgrau(1,i,j)     = micro%meltgrau(2,i,j)
  micro%melthail(1,i,j)     = micro%melthail(2,i,j)
  micro%rimecldsnow(1,i,j)  = micro%rimecldsnow(2,i,j)
  micro%rimecldaggr(1,i,j)  = micro%rimecldaggr(2,i,j)
  micro%rimecldgrau(1,i,j)  = micro%rimecldgrau(2,i,j)
  micro%rimecldhail(1,i,j)  = micro%rimecldhail(2,i,j)
  micro%rain2pr(1,i,j)      = micro%rain2pr(2,i,j)
  micro%rain2sn(1,i,j)      = micro%rain2sn(2,i,j)
  micro%rain2ag(1,i,j)      = micro%rain2ag(2,i,j)
  micro%rain2gr(1,i,j)      = micro%rain2gr(2,i,j)
  micro%rain2ha(1,i,j)      = micro%rain2ha(2,i,j)
  micro%rain2ha_xtra(1,i,j) = micro%rain2ha_xtra(2,i,j)
  micro%aggrselfpris(1,i,j) = micro%aggrselfpris(2,i,j)
  micro%aggrselfsnow(1,i,j) = micro%aggrselfsnow(2,i,j)
  micro%aggrprissnow(1,i,j) = micro%aggrprissnow(2,i,j)
endif
!COMPUTE AND OUTPUT MICRO BUDGET PROCESSES (totals)
if(imbudtot>=1) then
  micro%nuccldrt(1,i,j)    = micro%nuccldrt(2,i,j)
  micro%nuccldct(1,i,j)    = micro%nuccldct(2,i,j)
  micro%cld2raint(1,i,j)   = micro%cld2raint(2,i,j)
  micro%ice2raint(1,i,j)   = micro%ice2raint(2,i,j)
  micro%nucicert(1,i,j)    = micro%nucicert(2,i,j)
  micro%nucicect(1,i,j)    = micro%nucicect(2,i,j)
  micro%vapliqt(1,i,j)     = micro%vapliqt(2,i,j)
  micro%vapicet(1,i,j)     = micro%vapicet(2,i,j)
  micro%melticet(1,i,j)    = micro%melticet(2,i,j)
  micro%rimecldt(1,i,j)    = micro%rimecldt(2,i,j)
  micro%rain2icet(1,i,j)   = micro%rain2icet(2,i,j)
  micro%aggregatet(1,i,j)  = micro%aggregatet(2,i,j)
  micro%latheatvapt(1,i,j) = micro%latheatvapt(2,i,j)
  micro%latheatfrzt(1,i,j) = micro%latheatfrzt(2,i,j)
endif
if(imbudtot==2) then
  micro%inuchomrt(1,i,j)     = micro%inuchomrt(2,i,j)
  micro%inuchomct(1,i,j)     = micro%inuchomct(2,i,j)
  micro%inuccontrt(1,i,j)    = micro%inuccontrt(2,i,j)
  micro%inuccontct(1,i,j)    = micro%inuccontct(2,i,j)
  micro%inucifnrt(1,i,j)     = micro%inucifnrt(2,i,j)
  micro%inucifnct(1,i,j)     = micro%inucifnct(2,i,j)
  micro%inuchazrt(1,i,j)     = micro%inuchazrt(2,i,j)
  micro%inuchazct(1,i,j)     = micro%inuchazct(2,i,j)
  micro%vapcldt(1,i,j)       = micro%vapcldt(2,i,j)
  micro%vapraint(1,i,j)      = micro%vapraint(2,i,j)
  micro%vapprist(1,i,j)      = micro%vapprist(2,i,j)
  micro%vapsnowt(1,i,j)      = micro%vapsnowt(2,i,j)
  micro%vapaggrt(1,i,j)      = micro%vapaggrt(2,i,j)
  micro%vapgraut(1,i,j)      = micro%vapgraut(2,i,j)
  micro%vaphailt(1,i,j)      = micro%vaphailt(2,i,j)
  micro%vapdrizt(1,i,j)      = micro%vapdrizt(2,i,j)
  micro%meltprist(1,i,j)     = micro%meltprist(2,i,j)
  micro%meltsnowt(1,i,j)     = micro%meltsnowt(2,i,j)
  micro%meltaggrt(1,i,j)     = micro%meltaggrt(2,i,j)
  micro%meltgraut(1,i,j)     = micro%meltgraut(2,i,j)
  micro%melthailt(1,i,j)     = micro%melthailt(2,i,j)
  micro%rimecldsnowt(1,i,j)  = micro%rimecldsnowt(2,i,j)
  micro%rimecldaggrt(1,i,j)  = micro%rimecldaggrt(2,i,j)
  micro%rimecldgraut(1,i,j)  = micro%rimecldgraut(2,i,j)
  micro%rimecldhailt(1,i,j)  = micro%rimecldhailt(2,i,j)
  micro%rain2prt(1,i,j)      = micro%rain2prt(2,i,j)
  micro%rain2snt(1,i,j)      = micro%rain2snt(2,i,j)
  micro%rain2agt(1,i,j)      = micro%rain2agt(2,i,j)
  micro%rain2grt(1,i,j)      = micro%rain2grt(2,i,j)
  micro%rain2hat(1,i,j)      = micro%rain2hat(2,i,j)
  micro%rain2ha_xtrat(1,i,j) = micro%rain2ha_xtrat(2,i,j)
  micro%aggrselfprist(1,i,j) = micro%aggrselfprist(2,i,j)
  micro%aggrselfsnowt(1,i,j) = micro%aggrselfsnowt(2,i,j)
  micro%aggrprissnowt(1,i,j) = micro%aggrprissnowt(2,i,j)
endif

return
end subroutine copyback

!******************************************************************************
! This F90 subroutine computes the change to potential temperature associated
! with all of the RAMS microphysics vapor diffusion processes. It assumes
! that the heat storage and mixing ratios were stored before calls to vapor
! diffusion and accounts for fractional liquid/ice on hail and graupel. It uses
! the "jnmb" vector to determine whether changes to a given species are
! predicted. At this time, it does not differentiate between liquid and
! ice diffusion, nor does it keep track of condensation vs. evaporation.

subroutine calc_lhr_vap(m1,k1,k2,i,j,latheatvap,latheatvapt)

use rconstants
!use micphys

implicit none

integer :: m1,lcat,k,i,j
real, dimension(m1) :: latheatvap,latheatvapt
integer, dimension(11) :: k1,k2
real :: temp,fracliq1,fracliq2

if(imbudget>=1 .or. imbudtot>=1)then

 ! For each category, compute the change due to vapor transfer and apply to theta
 do lcat=1,ncat
  if (jnmb(lcat) .ge. 1) then
    do k=k1(lcat),k2(lcat)

      ! Vapor to liquid
      if (lcat .eq. 1 .or. lcat .eq. 2 .or. lcat .eq. 8) then
        if (lhrtheta) then
            latheatvap(k) = latheatvap(k) + (alvl/pitot(k))*(rx(k,lcat)-rx_lhr(k,lcat))
            latheatvapt(k) = latheatvapt(k) + (alvl/pitot(k))*(rx(k,lcat)-rx_lhr(k,lcat))
          else
            latheatvap(k) = latheatvap(k) +  alvl * (rx(k,lcat)-rx_lhr(k,lcat))
            latheatvapt(k) = latheatvapt(k) +  alvl * (rx(k,lcat)-rx_lhr(k,lcat))
        endif
      endif

      ! Vapor to ice
      if (lcat .eq. 3 .or. lcat .eq. 4 .or. lcat .eq. 5) then
        if (lhrtheta) then
            latheatvap(k) = latheatvap(k) + (alvi/pitot(k))*(rx(k,lcat)-rx_lhr(k,lcat))
            latheatvapt(k) = latheatvapt(k) + (alvi/pitot(k))*(rx(k,lcat)-rx_lhr(k,lcat))
          else
            latheatvap(k) = latheatvap(k) +  alvi * (rx(k,lcat)-rx_lhr(k,lcat))
            latheatvapt(k) = latheatvapt(k) +  alvi * (rx(k,lcat)-rx_lhr(k,lcat))
        endif
      endif

      ! Account for fractional liquid in graupel and hail categories
      if (lcat .eq. 6 .or. lcat .eq. 7) then
        ! Get fraction of liquid from old values
        call qtc(qx_lhr(k,lcat),temp,fracliq1)
        ! Get fraction of liquid from new values
        call qtc(qx(k,lcat),    temp,fracliq2)
        ! Update theta for liquid portion
        if (lhrtheta) then
            latheatvap(k) = latheatvap(k) + (alvl/pitot(k)) &
              * (rx(k,lcat)*fracliq2-rx_lhr(k,lcat)*fracliq1)
            latheatvapt(k) = latheatvapt(k) + (alvl/pitot(k)) &
              * (rx(k,lcat)*fracliq2-rx_lhr(k,lcat)*fracliq1)
        else
            latheatvap(k) = latheatvap(k) +  alvl &
              * (rx(k,lcat)*fracliq2-rx_lhr(k,lcat)*fracliq1)
            latheatvapt(k) = latheatvapt(k) +  alvl &
              * (rx(k,lcat)*fracliq2-rx_lhr(k,lcat)*fracliq1)
        endif
        ! Update theta for ice portion
        if (lhrtheta) then
            latheatvap(k) = latheatvap(k) + (alvi/pitot(k)) &
              * (rx(k,lcat)*(1.-fracliq2)-rx_lhr(k,lcat)*(1.-fracliq1))
            latheatvapt(k) = latheatvapt(k) + (alvi/pitot(k)) &
              * (rx(k,lcat)*(1.-fracliq2)-rx_lhr(k,lcat)*(1.-fracliq1))
        else
            latheatvap(k) = latheatvap(k) +  alvi &
              * (rx(k,lcat)*(1.-fracliq2)-rx_lhr(k,lcat)*(1.-fracliq1))
            latheatvapt(k) = latheatvapt(k) +  alvi &
              * (rx(k,lcat)*(1.-fracliq2)-rx_lhr(k,lcat)*(1.-fracliq1))
        endif
      endif

    enddo
  endif
 enddo

endif

return
end subroutine calc_lhr_vap

!******************************************************************************
! This F90 subroutine computes the change to potential temperature associated
! with all of the RAMS microphysics collision freezing/melting processes. It assumes
! that the heat storage and mixing ratios were stored before calls to colxfers
! and accounts for fractional liquid/ice on hail and graupel. It uses
! the "jnmb" vector to determine whether changes to a given species are
! predicted. At this time, it does not differentiate between freezing and
! melting--it only computes the net effect of all freezing/melting processes
! on the potential temperature.

subroutine calc_lhr_collmelt(m1,k1,k2,i,j,latheatfrz,latheatfrzt)

use rconstants
!use micphys

implicit none

integer :: m1,lcat,k,i,j
integer, dimension(11) :: k1,k2
real :: fracliq1, fracliq2, temp
real, dimension(m1) :: rliq1,rice1,rliq2,rice2,latheatfrz,latheatfrzt

if(imbudget>=1 .or. imbudtot>=1)then

 ! Compute total liquid and ice amounts before and after collisions/melting
 rliq1 = 0.
 rice1 = 0.
 rliq2 = 0.
 rice2 = 0.
 do lcat=1,ncat
  if (jnmb(lcat) .ge. 1) then
    do k=2,m1
      if (lcat.eq.1 .or. lcat.eq.2 .or. lcat.eq.8) then ! liquid only
        rliq1(k) = rliq1(k) + rx_lhr(k,lcat)
        rliq2(k) = rliq2(k) + rx(k,lcat)
      elseif (lcat.eq.3 .or. lcat.eq.4 .or. lcat.eq.5) then ! ice only
        rice1(k) = rice1(k) + rx_lhr(k,lcat)
        rice2(k) = rice2(k) + rx(k,lcat)
      elseif (lcat.eq.6 .or. lcat.eq.7) then !mixed phase
        !Get liquid fractions
        call qtc(qx_lhr(k,lcat),temp,fracliq1)
        call qtc(qx(k,lcat),    temp,fracliq2)
        !Accumulate liquid and ice
        rliq1(k) = rliq1(k) + rx_lhr(k,lcat) * fracliq1
        rliq2(k) = rliq2(k) + rx(k,lcat)     * fracliq2
        rice1(k) = rice1(k) + rx_lhr(k,lcat) * (1.-fracliq1)
        rice2(k) = rice2(k) + rx(k,lcat)     * (1.-fracliq2)
      endif
    enddo
  endif
 enddo

 ! Heating/cooling is due to change in net liquid vs. ice mass only
 do k=2,m1
   if (lhrtheta) then
    latheatfrz(k) = latheatfrz(k) + (alli/pitot(k))*(rice2(k)-rice1(k))
    latheatfrzt(k) = latheatfrzt(k) + (alli/pitot(k))*(rice2(k)-rice1(k))
   else
    latheatfrz(k) = latheatfrz(k) +  alli * (rice2(k)-rice1(k))
    latheatfrzt(k) = latheatfrzt(k) +  alli * (rice2(k)-rice1(k))
   endif
 enddo

endif

return
end subroutine calc_lhr_collmelt

!******************************************************************************
! This F90 subroutine computes the change to potential temperature associated
! with cloud nucleation in RAMS microphysics. It assumes that the mass mixing
! ratios were stored before calls to cldnuc (hail/graupel are not considered
! here, and so there is no need to store the heat storage vectors). It uses
! the "jnmb" vector to determine whether changes to a given species are
! predicted. At this time, it does not differentiate between nucleation and
! evaporation--it only computes the net effect of the nucleation processes
! on the potential temperature.

subroutine calc_lhr_cldnuc(m1,k1,k2,i,j,latheatvap,latheatvapt)

use rconstants
!use micphys

implicit none

integer :: m1,lcat,k,i,j
integer, dimension(11) :: k1,k2
real, dimension(m1) :: latheatvap,latheatvapt

if(imbudget>=1 .or. imbudtot>=1)then

 ! Cloud nucleation only changes the mixing ratio
 ! (and number concentration) of cloud1 and cloud2
 if (jnmb(1) .ge. 1) then
  do k=k1(1),k2(1)
    if (lhrtheta) then
      latheatvap(k) = latheatvap(k) + (alvl/pitot(k))*(rx(k,1)-rx_lhr(k,1))
      latheatvapt(k) = latheatvapt(k) + (alvl/pitot(k))*(rx(k,1)-rx_lhr(k,1))
    else
      latheatvap(k) = latheatvap(k) +  alvl * (rx(k,1)-rx_lhr(k,1))
      latheatvapt(k) = latheatvapt(k) +  alvl * (rx(k,1)-rx_lhr(k,1))
    endif
  enddo
 endif
 if (jnmb(8) .ge. 1) then
  do k=k1(8),k2(8)
    if (lhrtheta) then
      latheatvap(k) = latheatvap(k) + (alvl/pitot(k))*(rx(k,8)-rx_lhr(k,8))
      latheatvapt(k) = latheatvapt(k) + (alvl/pitot(k))*(rx(k,8)-rx_lhr(k,8))
    else
      latheatvap(k) = latheatvap(k) +  alvl * (rx(k,8)-rx_lhr(k,8))
      latheatvapt(k) = latheatvapt(k) +  alvl * (rx(k,8)-rx_lhr(k,8))
    endif
  enddo
 endif

endif

return
end subroutine calc_lhr_cldnuc

!******************************************************************************
! This F90 subroutine computes the change to potential temperature associated
! with ice nucleation in RAMS microphysics. It assumes that the mass mixing
! ratios were stored before calls to icenuc (hail/graupel are not considered
! here, and so there is no need to store the heat storage vectors). It uses
! the "jnmb" vector to determine whether changes to a given species are
! predicted. At this time, it does not differentiate between nucleation and
! sublimation--it only computes the net effect of the nucleation processes
! on the potential temperature.

subroutine calc_lhr_icenuc(m1,k1,k2,i,j,latheatvap,latheatvapt)

use rconstants
!use micphys

implicit none

integer :: m1,lcat,k,i,j
integer, dimension(11) :: k1,k2
real, dimension(m1) :: latheatvap,latheatvapt
real, dimension(nzpmax) :: vapdep

if(imbudget>=1 .or. imbudtot>=1)then

 ! Ice nucleation only changes the mixing ratio
 ! (and number concentration) of pristine ice
 if (jnmb(3) .ge. 1) then
  ! However, some comes from contact nucleation of cloud1 and cloud2,
  ! and some from vapor deposition...
  ! These need to be treated seperately due to different latent heats...

  ! First, changes in cloud1 and cloud2 (contact nucleation)
  if (jnmb(1) .ge. 1) then
    do k=k1(1),k2(1) ! Loss of cloud1
      if (lhrtheta) then
        latheatvap(k) = latheatvap(k) - (alli/pitot(k))*(rx(k,1)-rx_lhr(k,1))
        latheatvapt(k) = latheatvapt(k) - (alli/pitot(k))*(rx(k,1)-rx_lhr(k,1))
      else
        latheatvap(k) = latheatvap(k) - alli * (rx(k,1)-rx_lhr(k,1))
        latheatvapt(k) = latheatvapt(k) - alli * (rx(k,1)-rx_lhr(k,1))
      endif
    enddo
  endif
  if (jnmb(8) .ge. 1) then
    do k=k1(8),k2(8) ! Loss of cloud2
      if (lhrtheta) then
        latheatvap(k) = latheatvap(k) - (alli/pitot(k))*(rx(k,8)-rx_lhr(k,8))
        latheatvapt(k) = latheatvapt(k) - (alli/pitot(k))*(rx(k,8)-rx_lhr(k,8))
      else
        latheatvap(k) = latheatvap(k) - alli * (rx(k,8)-rx_lhr(k,8))
        latheatvapt(k) = latheatvapt(k) - alli * (rx(k,8)-rx_lhr(k,8))
      endif
    enddo
  endif

  ! Then, changes in pristine ice minus changes in
  ! cloud1 and cloud2 (vapor deposition)
  do k=k1(3),k2(3)
    ! First, compute total change in pristine ice (should be positive)
    vapdep(k) = (rx(k,3)-rx_lhr(k,3))
    ! Then, subtract off contributions from cloud1 and cloud2
    if (jnmb(1) .ge. 1) vapdep(k) = vapdep(k) - (rx_lhr(k,1)-rx(k,1))
    if (jnmb(8) .ge. 1) vapdep(k) = vapdep(k) - (rx_lhr(k,8)-rx(k,8))
    ! Finally, compute change to theta due to vapor deposition...
    if (lhrtheta) then
      latheatvap(k) = latheatvap(k) + (alvi/pitot(k)) * vapdep(k)
      latheatvapt(k) = latheatvapt(k) + (alvi/pitot(k)) * vapdep(k)
    else
      latheatvap(k) = latheatvap(k) +  alvi * vapdep(k)
      latheatvapt(k) = latheatvapt(k) +  alvi * vapdep(k)
    endif
  enddo

 endif ! End test for whether pristine ice is predicted

endif

return
end subroutine calc_lhr_icenuc

!******************************************************************************
subroutine aero_nuc_tab(rrv,rrvlsair,eps1,eps2,wtw1,wtw2,wtcon1,wtcon2 &
                    ,jrg1,jrg2,epstab,jw,jconcen,jtemp,rgccn1,tabvalue)

implicit none

integer :: iw,iconc
integer :: epstab,jw,jconcen,jtemp,rgccn1,ss
real :: eps1,eps2,wtw1,wtw2,wtcon1,wtcon2,jrg1,jrg2
real :: tabvalue,nucvalue,nucvalue1,nucvalue2,ssvalue,ssvalue1,ssvalue2
real, dimension(9,9,7,7,7) :: cldnuctab,supersat
real :: supsatpcnt,ssvalue11,ssvalue12,ssvalue21,ssvalue22
real :: mult1,mult2,nucss1,nucss2,nucss,rrv,rrvlsair

!Saleeby (6/3/02)
!The following tables are percent of prescibed CCN to activate.
!Table dimensions are:
! (median radius,updraft,concentration,temperature,solubility).
!Large aerosols or GCCN number are prescibed by the user and these
!automatically nucleate and they have no lookup tables.
!Small aerosols or CCN with median radius from 0.01 - 0.96 microns
!are user prescribed in concentration from 1 - 10,000 /cm3.
!Median radius is in centimeters as is the parcel model.

!cldnuctab(a,b,c,d,e) as produced in parcel model
!a = 1-9 median radius of aerosol dist (.01,.02,.04,.08,.16,.32,.48,.64,.96 microns)
!b = 1-9 vertical velocity (1 - 10,000 cm/s)
!c = 1-7 aerosol concentration (1 - 10,000 /cm3)
!d = 1-7 air temperature (-30 to 30 C)
!e = 1-7 epsilon solubility fraction (0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0)

!Input W (m/s), rg (cm), concentration (/cm3), T (C)

data ((cldnuctab( 1,iw,iconc,1, 1),iw=1,9),iconc=1,7)/    &
0.0167, 0.0467, 0.1256, 0.2824, 0.5471, 0.8354, 0.9846, 0.9970, 0.9970,  &
0.0088, 0.0272, 0.0759, 0.1944, 0.4175, 0.7147, 0.9391, 0.9970, 0.9970,  &
0.0050, 0.0167, 0.0508, 0.1342, 0.2963, 0.5792, 0.8548, 0.9861, 0.9970,  &
0.0031, 0.0099, 0.0328, 0.0883, 0.2177, 0.4496, 0.7548, 0.9483, 0.9970,  &
0.0018, 0.0063, 0.0204, 0.0599, 0.1625, 0.3551, 0.6418, 0.8809, 0.9899,  &
0.0011, 0.0039, 0.0136, 0.0428, 0.1174, 0.2687, 0.5146, 0.7675, 0.9483,  &
0.0006, 0.0024, 0.0088, 0.0299, 0.0819, 0.2059, 0.3551, 0.4983, 0.7007 /
data ((cldnuctab( 1,iw,iconc,2, 1),iw=1,9),iconc=1,7)/    &
0.0110, 0.0328, 0.0883, 0.2059, 0.4335, 0.7418, 0.9525, 0.9970, 0.9970,  &
0.0063, 0.0185, 0.0508, 0.1342, 0.3106, 0.5792, 0.8639, 0.9875, 0.9970,  &
0.0035, 0.0110, 0.0328, 0.0883, 0.2177, 0.4496, 0.7418, 0.9483, 0.9970,  &
0.0018, 0.0063, 0.0204, 0.0552, 0.1527, 0.3400, 0.6108, 0.8726, 0.9887,  &
0.0011, 0.0035, 0.0122, 0.0359, 0.1021, 0.2425, 0.4983, 0.7798, 0.9600,  &
0.0006, 0.0021, 0.0079, 0.0248, 0.0702, 0.1834, 0.3859, 0.6718, 0.8964,  &
0.0003, 0.0012, 0.0050, 0.0167, 0.0508, 0.1342, 0.2963, 0.5308, 0.7418 /
data ((cldnuctab( 1,iw,iconc,3, 1),iw=1,9),iconc=1,7)/    &
0.0088, 0.0248, 0.0649, 0.1625, 0.3551, 0.6418, 0.9036, 0.9942, 0.9970,  &
0.0044, 0.0136, 0.0392, 0.1021, 0.2425, 0.4820, 0.7675, 0.9634, 0.9970,  &
0.0024, 0.0071, 0.0225, 0.0649, 0.1625, 0.3400, 0.6264, 0.8889, 0.9919,  &
0.0012, 0.0039, 0.0136, 0.0392, 0.1021, 0.2425, 0.4983, 0.7798, 0.9600,  &
0.0007, 0.0024, 0.0079, 0.0248, 0.0702, 0.1728, 0.3704, 0.6569, 0.9036,  &
0.0003, 0.0014, 0.0044, 0.0151, 0.0467, 0.1256, 0.2824, 0.5471, 0.8143,  &
0.0002, 0.0008, 0.0031, 0.0099, 0.0299, 0.0883, 0.2177, 0.4335, 0.6864 /
data ((cldnuctab( 1,iw,iconc,4, 1),iw=1,9),iconc=1,7)/    &
0.0071, 0.0204, 0.0552, 0.1342, 0.2963, 0.5632, 0.8453, 0.9846, 0.9970,  &
0.0035, 0.0110, 0.0299, 0.0819, 0.1944, 0.4016, 0.6864, 0.9286, 0.9968,  &
0.0018, 0.0056, 0.0167, 0.0467, 0.1174, 0.2824, 0.5308, 0.8143, 0.9747,  &
0.0009, 0.0031, 0.0099, 0.0272, 0.0759, 0.1834, 0.4016, 0.6864, 0.9168,  &
0.0005, 0.0016, 0.0056, 0.0167, 0.0467, 0.1256, 0.2824, 0.5471, 0.8250,  &
0.0003, 0.0009, 0.0031, 0.0099, 0.0299, 0.0883, 0.2059, 0.4335, 0.7147,  &
0.0001, 0.0005, 0.0018, 0.0063, 0.0204, 0.0599, 0.1527, 0.3251, 0.5951 /
data ((cldnuctab( 1,iw,iconc,5, 1),iw=1,9),iconc=1,7)/    &
0.0056, 0.0167, 0.0467, 0.1174, 0.2554, 0.4983, 0.7917, 0.9695, 0.9970,  &
0.0031, 0.0088, 0.0248, 0.0649, 0.1625, 0.3551, 0.6264, 0.8889, 0.9909,  &
0.0014, 0.0044, 0.0136, 0.0392, 0.1021, 0.2299, 0.4658, 0.7548, 0.9525,  &
0.0007, 0.0024, 0.0071, 0.0225, 0.0599, 0.1527, 0.3251, 0.5951, 0.8639,  &
0.0004, 0.0012, 0.0039, 0.0122, 0.0359, 0.0950, 0.2299, 0.4496, 0.7418,  &
0.0002, 0.0006, 0.0024, 0.0071, 0.0225, 0.0599, 0.1527, 0.3400, 0.6108,  &
0.0001, 0.0003, 0.0012, 0.0044, 0.0136, 0.0392, 0.1095, 0.2554, 0.4983 /
data ((cldnuctab( 1,iw,iconc,6, 1),iw=1,9),iconc=1,7)/    &
0.0050, 0.0151, 0.0428, 0.1021, 0.2299, 0.4658, 0.7548, 0.9564, 0.9970,  &
0.0027, 0.0079, 0.0225, 0.0599, 0.1433, 0.3106, 0.5792, 0.8548, 0.9846,  &
0.0012, 0.0039, 0.0122, 0.0328, 0.0883, 0.2059, 0.4175, 0.7007, 0.9286,  &
0.0006, 0.0021, 0.0063, 0.0185, 0.0508, 0.1256, 0.2824, 0.5308, 0.8143,  &
0.0003, 0.0011, 0.0035, 0.0099, 0.0299, 0.0759, 0.1834, 0.3859, 0.6718,  &
0.0002, 0.0005, 0.0018, 0.0056, 0.0167, 0.0467, 0.1256, 0.2824, 0.5308,  &
0.0001, 0.0003, 0.0009, 0.0031, 0.0099, 0.0299, 0.0819, 0.1944, 0.4016 /
data ((cldnuctab( 1,iw,iconc,7, 1),iw=1,9),iconc=1,7)/    &
0.0050, 0.0136, 0.0392, 0.0950, 0.2177, 0.4335, 0.7284, 0.9438, 0.9970,  &
0.0024, 0.0071, 0.0204, 0.0552, 0.1342, 0.2824, 0.5471, 0.8250, 0.9770,  &
0.0012, 0.0035, 0.0110, 0.0299, 0.0759, 0.1834, 0.3704, 0.6569, 0.9036,  &
0.0005, 0.0018, 0.0056, 0.0167, 0.0428, 0.1095, 0.2554, 0.4820, 0.7675,  &
0.0003, 0.0009, 0.0031, 0.0088, 0.0248, 0.0649, 0.1625, 0.3400, 0.6108,  &
0.0001, 0.0005, 0.0016, 0.0050, 0.0151, 0.0392, 0.1021, 0.2425, 0.4658,  &
0.0001, 0.0002, 0.0008, 0.0027, 0.0079, 0.0248, 0.0649, 0.1625, 0.3400 /
data ((cldnuctab( 1,iw,iconc,1, 2),iw=1,9),iconc=1,7)/    &
0.0248, 0.0649, 0.1527, 0.3400, 0.6108, 0.8726, 0.9887, 0.9970, 0.9970,  &
0.0136, 0.0392, 0.1021, 0.2425, 0.4658, 0.7548, 0.9483, 0.9970, 0.9970,  &
0.0088, 0.0248, 0.0702, 0.1728, 0.3551, 0.6264, 0.8809, 0.9887, 0.9970,  &
0.0050, 0.0167, 0.0467, 0.1174, 0.2687, 0.5146, 0.7917, 0.9600, 0.9970,  &
0.0031, 0.0110, 0.0328, 0.0883, 0.2059, 0.4175, 0.6864, 0.9036, 0.9919,  &
0.0018, 0.0071, 0.0225, 0.0599, 0.1527, 0.3251, 0.5792, 0.7917, 0.9600,  &
0.0012, 0.0044, 0.0151, 0.0428, 0.1174, 0.2425, 0.4016, 0.5308, 0.7284 /
data ((cldnuctab( 1,iw,iconc,2, 2),iw=1,9),iconc=1,7)/    &
0.0167, 0.0467, 0.1174, 0.2554, 0.4983, 0.7798, 0.9600, 0.9970, 0.9970,  &
0.0099, 0.0272, 0.0702, 0.1728, 0.3704, 0.6418, 0.8889, 0.9909, 0.9970,  &
0.0056, 0.0167, 0.0467, 0.1174, 0.2687, 0.5146, 0.7798, 0.9600, 0.9970,  &
0.0031, 0.0099, 0.0299, 0.0819, 0.1944, 0.4016, 0.6718, 0.9036, 0.9919,  &
0.0018, 0.0063, 0.0185, 0.0552, 0.1433, 0.3106, 0.5632, 0.8143, 0.9695,  &
0.0011, 0.0039, 0.0122, 0.0392, 0.1021, 0.2299, 0.4496, 0.7147, 0.9168,  &
0.0006, 0.0024, 0.0079, 0.0248, 0.0702, 0.1728, 0.3551, 0.5792, 0.7798 /
data ((cldnuctab( 1,iw,iconc,3, 2),iw=1,9),iconc=1,7)/    &
0.0122, 0.0359, 0.0883, 0.2059, 0.4016, 0.6864, 0.9229, 0.9959, 0.9970,  &
0.0071, 0.0204, 0.0552, 0.1342, 0.2963, 0.5471, 0.8143, 0.9722, 0.9970,  &
0.0039, 0.0110, 0.0328, 0.0883, 0.2059, 0.4016, 0.6864, 0.9103, 0.9935,  &
0.0021, 0.0071, 0.0204, 0.0552, 0.1342, 0.2963, 0.5471, 0.8143, 0.9695,  &
0.0012, 0.0039, 0.0122, 0.0359, 0.0950, 0.2177, 0.4335, 0.7147, 0.9229,  &
0.0007, 0.0024, 0.0079, 0.0248, 0.0649, 0.1625, 0.3400, 0.6108, 0.8453,  &
0.0004, 0.0014, 0.0050, 0.0167, 0.0467, 0.1174, 0.2687, 0.4983, 0.7284 /
data ((cldnuctab( 1,iw,iconc,4, 2),iw=1,9),iconc=1,7)/    &
0.0099, 0.0272, 0.0702, 0.1728, 0.3551, 0.6264, 0.8726, 0.9875, 0.9970,  &
0.0056, 0.0151, 0.0428, 0.1095, 0.2425, 0.4658, 0.7418, 0.9438, 0.9970,  &
0.0027, 0.0088, 0.0248, 0.0649, 0.1625, 0.3400, 0.5951, 0.8548, 0.9811,  &
0.0016, 0.0050, 0.0151, 0.0428, 0.1021, 0.2299, 0.4658, 0.7284, 0.9340,  &
0.0008, 0.0027, 0.0088, 0.0248, 0.0702, 0.1625, 0.3551, 0.6108, 0.8548,  &
0.0005, 0.0016, 0.0050, 0.0167, 0.0467, 0.1174, 0.2687, 0.4983, 0.7675,  &
0.0003, 0.0009, 0.0031, 0.0110, 0.0299, 0.0819, 0.1944, 0.4016, 0.6569 /
data ((cldnuctab( 1,iw,iconc,5, 2),iw=1,9),iconc=1,7)/    &
0.0088, 0.0248, 0.0599, 0.1433, 0.3106, 0.5632, 0.8354, 0.9770, 0.9970,  &
0.0044, 0.0136, 0.0359, 0.0883, 0.2059, 0.4016, 0.6864, 0.9103, 0.9935,  &
0.0024, 0.0071, 0.0204, 0.0552, 0.1342, 0.2824, 0.5308, 0.7917, 0.9634,  &
0.0012, 0.0039, 0.0110, 0.0328, 0.0819, 0.1944, 0.3859, 0.6569, 0.8889,  &
0.0006, 0.0021, 0.0063, 0.0185, 0.0508, 0.1256, 0.2824, 0.5146, 0.7917,  &
0.0003, 0.0012, 0.0039, 0.0122, 0.0359, 0.0883, 0.2059, 0.4016, 0.6718,  &
0.0002, 0.0006, 0.0024, 0.0071, 0.0225, 0.0599, 0.1433, 0.3106, 0.5632 /
data ((cldnuctab( 1,iw,iconc,6, 2),iw=1,9),iconc=1,7)/    &
0.0079, 0.0225, 0.0552, 0.1342, 0.2824, 0.5308, 0.8032, 0.9666, 0.9970,  &
0.0039, 0.0110, 0.0328, 0.0819, 0.1834, 0.3704, 0.6418, 0.8809, 0.9887,  &
0.0021, 0.0063, 0.0167, 0.0467, 0.1095, 0.2425, 0.4820, 0.7548, 0.9438,  &
0.0011, 0.0035, 0.0099, 0.0272, 0.0702, 0.1625, 0.3400, 0.5951, 0.8548,  &
0.0005, 0.0016, 0.0056, 0.0151, 0.0428, 0.1095, 0.2425, 0.4496, 0.7284,  &
0.0003, 0.0009, 0.0031, 0.0088, 0.0272, 0.0702, 0.1625, 0.3400, 0.5951,  &
0.0002, 0.0005, 0.0016, 0.0056, 0.0167, 0.0467, 0.1174, 0.2554, 0.4820 /
data ((cldnuctab( 1,iw,iconc,7, 2),iw=1,9),iconc=1,7)/    &
0.0071, 0.0204, 0.0508, 0.1256, 0.2687, 0.4983, 0.7798, 0.9564, 0.9970,  &
0.0035, 0.0110, 0.0299, 0.0702, 0.1728, 0.3400, 0.6108, 0.8639, 0.9829,  &
0.0018, 0.0056, 0.0151, 0.0428, 0.1021, 0.2299, 0.4335, 0.7147, 0.9286,  &
0.0009, 0.0027, 0.0088, 0.0248, 0.0599, 0.1433, 0.3106, 0.5471, 0.8143,  &
0.0005, 0.0014, 0.0044, 0.0136, 0.0359, 0.0950, 0.2059, 0.4016, 0.6718,  &
0.0002, 0.0008, 0.0024, 0.0079, 0.0225, 0.0599, 0.1433, 0.2963, 0.5308,  &
0.0001, 0.0004, 0.0014, 0.0044, 0.0136, 0.0359, 0.0950, 0.2177, 0.4175 /
data ((cldnuctab( 1,iw,iconc,1, 3),iw=1,9),iconc=1,7)/    &
0.0359, 0.0950, 0.2059, 0.4175, 0.6718, 0.9036, 0.9919, 0.9970, 0.9970,  &
0.0225, 0.0599, 0.1433, 0.3106, 0.5471, 0.8032, 0.9634, 0.9970, 0.9970,  &
0.0136, 0.0392, 0.1021, 0.2299, 0.4335, 0.7007, 0.9103, 0.9919, 0.9970,  &
0.0088, 0.0272, 0.0702, 0.1728, 0.3400, 0.5951, 0.8354, 0.9722, 0.9970,  &
0.0056, 0.0185, 0.0508, 0.1256, 0.2687, 0.4983, 0.7548, 0.9286, 0.9942,  &
0.0039, 0.0122, 0.0359, 0.0950, 0.2059, 0.4016, 0.6418, 0.8354, 0.9695,  &
0.0024, 0.0079, 0.0248, 0.0649, 0.1527, 0.2963, 0.4496, 0.5308, 0.7548 /
data ((cldnuctab( 1,iw,iconc,2, 3),iw=1,9),iconc=1,7)/    &
0.0248, 0.0649, 0.1527, 0.3251, 0.5792, 0.8250, 0.9722, 0.9970, 0.9970,  &
0.0151, 0.0428, 0.1021, 0.2299, 0.4496, 0.7007, 0.9168, 0.9935, 0.9970,  &
0.0088, 0.0272, 0.0702, 0.1625, 0.3400, 0.5792, 0.8354, 0.9722, 0.9970,  &
0.0056, 0.0167, 0.0467, 0.1174, 0.2554, 0.4820, 0.7284, 0.9286, 0.9942,  &
0.0035, 0.0110, 0.0328, 0.0819, 0.1944, 0.3859, 0.6264, 0.8639, 0.9770,  &
0.0021, 0.0071, 0.0225, 0.0599, 0.1433, 0.3106, 0.5308, 0.7798, 0.9391,  &
0.0012, 0.0044, 0.0151, 0.0428, 0.1095, 0.2425, 0.4335, 0.6418, 0.8143 /
data ((cldnuctab( 1,iw,iconc,3, 3),iw=1,9),iconc=1,7)/    &
0.0185, 0.0508, 0.1256, 0.2687, 0.4820, 0.7548, 0.9438, 0.9970, 0.9970,  &
0.0110, 0.0299, 0.0759, 0.1834, 0.3551, 0.6108, 0.8548, 0.9811, 0.9970,  &
0.0063, 0.0185, 0.0508, 0.1256, 0.2687, 0.4820, 0.7418, 0.9340, 0.9959,  &
0.0035, 0.0110, 0.0328, 0.0819, 0.1944, 0.3859, 0.6264, 0.8639, 0.9791,  &
0.0021, 0.0071, 0.0204, 0.0552, 0.1342, 0.2963, 0.5308, 0.7798, 0.9438,  &
0.0012, 0.0044, 0.0136, 0.0392, 0.1021, 0.2299, 0.4335, 0.6864, 0.8889,  &
0.0008, 0.0027, 0.0088, 0.0272, 0.0702, 0.1728, 0.3400, 0.5792, 0.7798 /
data ((cldnuctab( 1,iw,iconc,4, 3),iw=1,9),iconc=1,7)/    &
0.0151, 0.0428, 0.1021, 0.2177, 0.4335, 0.6864, 0.9103, 0.9919, 0.9970,  &
0.0088, 0.0248, 0.0599, 0.1433, 0.2963, 0.5471, 0.8032, 0.9600, 0.9970,  &
0.0050, 0.0136, 0.0392, 0.0950, 0.2059, 0.4016, 0.6718, 0.8889, 0.9875,  &
0.0027, 0.0079, 0.0225, 0.0599, 0.1433, 0.3106, 0.5471, 0.7917, 0.9525,  &
0.0014, 0.0050, 0.0151, 0.0392, 0.1021, 0.2299, 0.4335, 0.6864, 0.8964,  &
0.0008, 0.0031, 0.0088, 0.0272, 0.0702, 0.1625, 0.3400, 0.5792, 0.8143,  &
0.0005, 0.0018, 0.0056, 0.0185, 0.0508, 0.1256, 0.2687, 0.4820, 0.7147 /
data ((cldnuctab( 1,iw,iconc,5, 3),iw=1,9),iconc=1,7)/    &
0.0136, 0.0359, 0.0883, 0.1944, 0.3859, 0.6418, 0.8809, 0.9846, 0.9970,  &
0.0071, 0.0204, 0.0508, 0.1256, 0.2687, 0.4820, 0.7418, 0.9391, 0.9959,  &
0.0039, 0.0110, 0.0299, 0.0759, 0.1728, 0.3551, 0.5951, 0.8453, 0.9747,  &
0.0021, 0.0063, 0.0185, 0.0508, 0.1174, 0.2554, 0.4658, 0.7284, 0.9229,  &
0.0011, 0.0035, 0.0110, 0.0299, 0.0819, 0.1834, 0.3551, 0.5951, 0.8354,  &
0.0006, 0.0021, 0.0071, 0.0204, 0.0552, 0.1256, 0.2687, 0.4983, 0.7418,  &
0.0003, 0.0012, 0.0044, 0.0122, 0.0359, 0.0950, 0.2059, 0.4016, 0.6418 /
data ((cldnuctab( 1,iw,iconc,6, 3),iw=1,9),iconc=1,7)/    &
0.0122, 0.0328, 0.0759, 0.1728, 0.3551, 0.6108, 0.8548, 0.9770, 0.9970,  &
0.0063, 0.0167, 0.0467, 0.1095, 0.2425, 0.4496, 0.7147, 0.9168, 0.9927,  &
0.0035, 0.0099, 0.0272, 0.0649, 0.1527, 0.3106, 0.5471, 0.8032, 0.9600,  &
0.0016, 0.0056, 0.0151, 0.0428, 0.1021, 0.2177, 0.4175, 0.6718, 0.8889,  &
0.0009, 0.0031, 0.0088, 0.0248, 0.0649, 0.1527, 0.3106, 0.5308, 0.7917,  &
0.0005, 0.0016, 0.0050, 0.0151, 0.0428, 0.1021, 0.2299, 0.4175, 0.6718,  &
0.0003, 0.0009, 0.0031, 0.0099, 0.0272, 0.0702, 0.1625, 0.3251, 0.5632 /
data ((cldnuctab( 1,iw,iconc,7, 3),iw=1,9),iconc=1,7)/    &
0.0110, 0.0299, 0.0759, 0.1625, 0.3400, 0.5792, 0.8354, 0.9722, 0.9970,  &
0.0056, 0.0167, 0.0428, 0.1021, 0.2177, 0.4175, 0.6864, 0.8964, 0.9887,  &
0.0031, 0.0088, 0.0248, 0.0599, 0.1433, 0.2963, 0.5146, 0.7798, 0.9483,  &
0.0016, 0.0044, 0.0136, 0.0359, 0.0883, 0.1944, 0.3859, 0.6264, 0.8639,  &
0.0008, 0.0027, 0.0079, 0.0204, 0.0552, 0.1342, 0.2687, 0.4983, 0.7418,  &
0.0005, 0.0014, 0.0044, 0.0136, 0.0359, 0.0883, 0.1944, 0.3704, 0.6108,  &
0.0002, 0.0008, 0.0027, 0.0079, 0.0225, 0.0599, 0.1342, 0.2824, 0.4983 /
data ((cldnuctab( 1,iw,iconc,1, 4),iw=1,9),iconc=1,7)/    &
0.0552, 0.1342, 0.2824, 0.4983, 0.7548, 0.9340, 0.9954, 0.9970, 0.9970,  &
0.0359, 0.0950, 0.2059, 0.4016, 0.6418, 0.8639, 0.9770, 0.9970, 0.9970,  &
0.0248, 0.0649, 0.1527, 0.3106, 0.5308, 0.7798, 0.9391, 0.9954, 0.9970,  &
0.0167, 0.0428, 0.1095, 0.2425, 0.4335, 0.6864, 0.8889, 0.9829, 0.9970,  &
0.0110, 0.0299, 0.0819, 0.1834, 0.3551, 0.5951, 0.8143, 0.9525, 0.9968,  &
0.0063, 0.0204, 0.0552, 0.1433, 0.2824, 0.4983, 0.7147, 0.8809, 0.9811,  &
0.0035, 0.0122, 0.0359, 0.0950, 0.2059, 0.3400, 0.4820, 0.5471, 0.7675 /
data ((cldnuctab( 1,iw,iconc,2, 4),iw=1,9),iconc=1,7)/    &
0.0392, 0.0950, 0.2177, 0.4175, 0.6569, 0.8809, 0.9829, 0.9970, 0.9970,  &
0.0248, 0.0649, 0.1527, 0.3106, 0.5308, 0.7798, 0.9438, 0.9964, 0.9970,  &
0.0151, 0.0428, 0.1095, 0.2299, 0.4335, 0.6718, 0.8809, 0.9829, 0.9970,  &
0.0099, 0.0299, 0.0759, 0.1728, 0.3400, 0.5792, 0.8032, 0.9525, 0.9968,  &
0.0063, 0.0185, 0.0552, 0.1256, 0.2687, 0.4820, 0.7147, 0.9036, 0.9875,  &
0.0039, 0.0122, 0.0359, 0.0950, 0.2059, 0.4016, 0.6264, 0.8354, 0.9600,  &
0.0024, 0.0079, 0.0248, 0.0649, 0.1625, 0.3106, 0.5146, 0.6864, 0.8639 /
data ((cldnuctab( 1,iw,iconc,3, 4),iw=1,9),iconc=1,7)/    &
0.0299, 0.0759, 0.1728, 0.3400, 0.5792, 0.8250, 0.9634, 0.9970, 0.9970,  &
0.0185, 0.0467, 0.1174, 0.2425, 0.4496, 0.7007, 0.9036, 0.9887, 0.9970,  &
0.0110, 0.0299, 0.0759, 0.1728, 0.3400, 0.5792, 0.8143, 0.9600, 0.9970,  &
0.0063, 0.0185, 0.0508, 0.1256, 0.2687, 0.4820, 0.7147, 0.9103, 0.9875,  &
0.0039, 0.0122, 0.0359, 0.0883, 0.2059, 0.3859, 0.6264, 0.8354, 0.9666,  &
0.0024, 0.0079, 0.0248, 0.0649, 0.1527, 0.3106, 0.5308, 0.7675, 0.9229,  &
0.0014, 0.0050, 0.0167, 0.0467, 0.1095, 0.2425, 0.4335, 0.6569, 0.8453 /
data ((cldnuctab( 1,iw,iconc,4, 4),iw=1,9),iconc=1,7)/    &
0.0248, 0.0599, 0.1433, 0.2963, 0.5146, 0.7675, 0.9391, 0.9954, 0.9970,  &
0.0136, 0.0359, 0.0950, 0.2059, 0.3859, 0.6264, 0.8548, 0.9747, 0.9970,  &
0.0079, 0.0225, 0.0599, 0.1342, 0.2824, 0.4983, 0.7548, 0.9286, 0.9927,  &
0.0044, 0.0136, 0.0392, 0.0950, 0.2059, 0.4016, 0.6418, 0.8548, 0.9722,  &
0.0027, 0.0088, 0.0248, 0.0649, 0.1527, 0.3106, 0.5308, 0.7675, 0.9340,  &
0.0016, 0.0056, 0.0167, 0.0467, 0.1095, 0.2425, 0.4335, 0.6718, 0.8726,  &
0.0009, 0.0035, 0.0110, 0.0299, 0.0819, 0.1834, 0.3551, 0.5792, 0.7917 /
data ((cldnuctab( 1,iw,iconc,5, 4),iw=1,9),iconc=1,7)/    &
0.0204, 0.0508, 0.1256, 0.2554, 0.4658, 0.7284, 0.9168, 0.9919, 0.9970,  &
0.0110, 0.0299, 0.0759, 0.1728, 0.3400, 0.5792, 0.8143, 0.9600, 0.9970,  &
0.0063, 0.0185, 0.0467, 0.1174, 0.2425, 0.4496, 0.6864, 0.8964, 0.9861,  &
0.0035, 0.0110, 0.0299, 0.0759, 0.1728, 0.3400, 0.5632, 0.8032, 0.9525,  &
0.0021, 0.0063, 0.0185, 0.0508, 0.1174, 0.2554, 0.4496, 0.7007, 0.8889,  &
0.0012, 0.0039, 0.0122, 0.0328, 0.0819, 0.1834, 0.3551, 0.5951, 0.8143,  &
0.0007, 0.0024, 0.0079, 0.0225, 0.0599, 0.1433, 0.2824, 0.4983, 0.7284 /
data ((cldnuctab( 1,iw,iconc,6, 4),iw=1,9),iconc=1,7)/    &
0.0185, 0.0467, 0.1095, 0.2425, 0.4335, 0.6864, 0.9036, 0.9875, 0.9970,  &
0.0099, 0.0272, 0.0702, 0.1527, 0.3106, 0.5471, 0.7917, 0.9483, 0.9959,  &
0.0056, 0.0151, 0.0392, 0.0950, 0.2177, 0.4016, 0.6418, 0.8639, 0.9770,  &
0.0031, 0.0088, 0.0248, 0.0649, 0.1433, 0.2963, 0.5146, 0.7548, 0.9286,  &
0.0016, 0.0050, 0.0151, 0.0392, 0.0950, 0.2059, 0.4016, 0.6418, 0.8548,  &
0.0009, 0.0031, 0.0088, 0.0272, 0.0649, 0.1527, 0.3106, 0.5308, 0.7548,  &
0.0005, 0.0018, 0.0056, 0.0167, 0.0467, 0.1095, 0.2299, 0.4335, 0.6569 /
data ((cldnuctab( 1,iw,iconc,7, 4),iw=1,9),iconc=1,7)/    &
0.0167, 0.0428, 0.1021, 0.2177, 0.4175, 0.6718, 0.8889, 0.9829, 0.9970,  &
0.0088, 0.0248, 0.0599, 0.1433, 0.2963, 0.5146, 0.7675, 0.9340, 0.9942,  &
0.0050, 0.0136, 0.0359, 0.0883, 0.1944, 0.3704, 0.6108, 0.8453, 0.9695,  &
0.0027, 0.0079, 0.0204, 0.0552, 0.1256, 0.2687, 0.4820, 0.7147, 0.9103,  &
0.0014, 0.0044, 0.0122, 0.0359, 0.0819, 0.1834, 0.3551, 0.5951, 0.8143,  &
0.0008, 0.0024, 0.0079, 0.0225, 0.0552, 0.1256, 0.2687, 0.4658, 0.7147,  &
0.0005, 0.0014, 0.0044, 0.0136, 0.0359, 0.0883, 0.1944, 0.3704, 0.5951 /
data ((cldnuctab( 1,iw,iconc,1, 5),iw=1,9),iconc=1,7)/    &
0.0702, 0.1625, 0.3251, 0.5632, 0.8032, 0.9525, 0.9968, 0.9970, 0.9970,  &
0.0508, 0.1174, 0.2554, 0.4496, 0.7007, 0.8964, 0.9846, 0.9970, 0.9970,  &
0.0328, 0.0819, 0.1834, 0.3704, 0.5951, 0.8250, 0.9564, 0.9968, 0.9970,  &
0.0225, 0.0599, 0.1433, 0.2824, 0.4983, 0.7418, 0.9168, 0.9875, 0.9970,  &
0.0151, 0.0428, 0.1021, 0.2299, 0.4175, 0.6569, 0.8548, 0.9666, 0.9970,  &
0.0088, 0.0272, 0.0759, 0.1728, 0.3400, 0.5632, 0.7675, 0.9036, 0.9875,  &
0.0044, 0.0151, 0.0467, 0.1174, 0.2299, 0.3704, 0.4983, 0.5308, 0.8143 /
data ((cldnuctab( 1,iw,iconc,2, 5),iw=1,9),iconc=1,7)/    &
0.0508, 0.1256, 0.2554, 0.4658, 0.7147, 0.9103, 0.9887, 0.9970, 0.9970,  &
0.0328, 0.0819, 0.1834, 0.3704, 0.5951, 0.8250, 0.9600, 0.9970, 0.9970,  &
0.0225, 0.0599, 0.1342, 0.2824, 0.4983, 0.7284, 0.9103, 0.9875, 0.9970,  &
0.0136, 0.0392, 0.0950, 0.2177, 0.4016, 0.6418, 0.8453, 0.9666, 0.9970,  &
0.0088, 0.0272, 0.0702, 0.1625, 0.3251, 0.5471, 0.7675, 0.9286, 0.9909,  &
0.0056, 0.0167, 0.0508, 0.1174, 0.2554, 0.4658, 0.6864, 0.8726, 0.9722,  &
0.0031, 0.0099, 0.0299, 0.0819, 0.1944, 0.3551, 0.5632, 0.7418, 0.8889 /
data ((cldnuctab( 1,iw,iconc,3, 5),iw=1,9),iconc=1,7)/    &
0.0392, 0.0950, 0.2059, 0.4016, 0.6418, 0.8639, 0.9747, 0.9970, 0.9970,  &
0.0248, 0.0599, 0.1433, 0.2963, 0.5146, 0.7548, 0.9286, 0.9919, 0.9970,  &
0.0151, 0.0392, 0.1021, 0.2177, 0.4016, 0.6418, 0.8548, 0.9722, 0.9970,  &
0.0099, 0.0272, 0.0702, 0.1625, 0.3106, 0.5308, 0.7675, 0.9340, 0.9919,  &
0.0056, 0.0185, 0.0467, 0.1174, 0.2425, 0.4496, 0.6864, 0.8726, 0.9747,  &
0.0035, 0.0110, 0.0328, 0.0819, 0.1944, 0.3704, 0.5951, 0.8032, 0.9438,  &
0.0021, 0.0071, 0.0204, 0.0599, 0.1433, 0.2963, 0.4983, 0.7147, 0.8726 /
data ((cldnuctab( 1,iw,iconc,4, 5),iw=1,9),iconc=1,7)/    &
0.0299, 0.0759, 0.1728, 0.3400, 0.5792, 0.8143, 0.9564, 0.9970, 0.9970,  &
0.0185, 0.0467, 0.1174, 0.2425, 0.4496, 0.6864, 0.8889, 0.9829, 0.9970,  &
0.0110, 0.0299, 0.0759, 0.1728, 0.3400, 0.5632, 0.8032, 0.9483, 0.9954,  &
0.0063, 0.0185, 0.0508, 0.1174, 0.2554, 0.4496, 0.7007, 0.8889, 0.9811,  &
0.0039, 0.0122, 0.0328, 0.0883, 0.1944, 0.3704, 0.5951, 0.8143, 0.9525,  &
0.0024, 0.0079, 0.0225, 0.0599, 0.1433, 0.2824, 0.4983, 0.7284, 0.9036,  &
0.0014, 0.0044, 0.0151, 0.0428, 0.1021, 0.2299, 0.4175, 0.6418, 0.8354 /
data ((cldnuctab( 1,iw,iconc,5, 5),iw=1,9),iconc=1,7)/    &
0.0272, 0.0649, 0.1527, 0.3106, 0.5308, 0.7798, 0.9391, 0.9942, 0.9970,  &
0.0151, 0.0392, 0.0950, 0.2059, 0.4016, 0.6418, 0.8548, 0.9722, 0.9970,  &
0.0088, 0.0248, 0.0599, 0.1433, 0.2824, 0.4983, 0.7418, 0.9229, 0.9899,  &
0.0050, 0.0151, 0.0392, 0.0950, 0.2059, 0.3859, 0.6264, 0.8453, 0.9666,  &
0.0031, 0.0088, 0.0248, 0.0649, 0.1527, 0.2963, 0.5146, 0.7548, 0.9168,  &
0.0018, 0.0056, 0.0167, 0.0467, 0.1095, 0.2299, 0.4175, 0.6569, 0.8548,  &
0.0011, 0.0035, 0.0110, 0.0299, 0.0759, 0.1728, 0.3400, 0.5632, 0.7798 /
data ((cldnuctab( 1,iw,iconc,6, 5),iw=1,9),iconc=1,7)/    &
0.0225, 0.0599, 0.1342, 0.2824, 0.4983, 0.7418, 0.9286, 0.9919, 0.9970,  &
0.0136, 0.0359, 0.0819, 0.1834, 0.3551, 0.5951, 0.8354, 0.9634, 0.9970,  &
0.0071, 0.0204, 0.0508, 0.1256, 0.2554, 0.4658, 0.7007, 0.8964, 0.9846,  &
0.0039, 0.0122, 0.0328, 0.0819, 0.1834, 0.3400, 0.5792, 0.8032, 0.9483,  &
0.0024, 0.0071, 0.0204, 0.0552, 0.1256, 0.2554, 0.4658, 0.7007, 0.8889,  &
0.0014, 0.0044, 0.0122, 0.0359, 0.0883, 0.1944, 0.3551, 0.5951, 0.8032,  &
0.0008, 0.0027, 0.0079, 0.0225, 0.0599, 0.1433, 0.2824, 0.4820, 0.7147 /
data ((cldnuctab( 1,iw,iconc,7, 5),iw=1,9),iconc=1,7)/    &
0.0204, 0.0552, 0.1256, 0.2687, 0.4820, 0.7284, 0.9168, 0.9887, 0.9970,  &
0.0122, 0.0328, 0.0759, 0.1728, 0.3400, 0.5792, 0.8143, 0.9564, 0.9959,  &
0.0063, 0.0185, 0.0467, 0.1095, 0.2299, 0.4335, 0.6718, 0.8809, 0.9791,  &
0.0035, 0.0099, 0.0272, 0.0702, 0.1625, 0.3106, 0.5308, 0.7675, 0.9340,  &
0.0021, 0.0063, 0.0167, 0.0428, 0.1095, 0.2299, 0.4175, 0.6569, 0.8639,  &
0.0012, 0.0035, 0.0110, 0.0299, 0.0702, 0.1625, 0.3251, 0.5308, 0.7675,  &
0.0006, 0.0021, 0.0063, 0.0185, 0.0508, 0.1174, 0.2425, 0.4335, 0.6718 /
data ((cldnuctab( 1,iw,iconc,1, 6),iw=1,9),iconc=1,7)/    &
0.0883, 0.1944, 0.3704, 0.6108, 0.8354, 0.9634, 0.9970, 0.9970, 0.9970,  &
0.0599, 0.1433, 0.2824, 0.4983, 0.7418, 0.9168, 0.9887, 0.9970, 0.9970,  &
0.0392, 0.1021, 0.2177, 0.4016, 0.6418, 0.8548, 0.9666, 0.9970, 0.9970,  &
0.0272, 0.0702, 0.1625, 0.3251, 0.5471, 0.7798, 0.9340, 0.9909, 0.9970,  &
0.0185, 0.0508, 0.1256, 0.2687, 0.4658, 0.7007, 0.8809, 0.9747, 0.9970,  &
0.0110, 0.0328, 0.0883, 0.2059, 0.3859, 0.5951, 0.7917, 0.9229, 0.9899,  &
0.0050, 0.0167, 0.0508, 0.1342, 0.2425, 0.3704, 0.5146, 0.5632, 0.8032 /
data ((cldnuctab( 1,iw,iconc,2, 6),iw=1,9),iconc=1,7)/    &
0.0649, 0.1433, 0.2963, 0.5146, 0.7548, 0.9286, 0.9919, 0.9970, 0.9970,  &
0.0392, 0.1021, 0.2177, 0.4016, 0.6418, 0.8548, 0.9695, 0.9970, 0.9970,  &
0.0272, 0.0702, 0.1625, 0.3106, 0.5308, 0.7675, 0.9286, 0.9909, 0.9970,  &
0.0185, 0.0467, 0.1174, 0.2425, 0.4496, 0.6864, 0.8726, 0.9747, 0.9970,  &
0.0110, 0.0328, 0.0883, 0.1944, 0.3704, 0.5951, 0.8032, 0.9438, 0.9935,  &
0.0071, 0.0204, 0.0599, 0.1433, 0.2963, 0.4983, 0.7284, 0.8964, 0.9791,  &
0.0035, 0.0110, 0.0359, 0.0950, 0.2177, 0.3859, 0.5951, 0.7548, 0.9103 /
data ((cldnuctab( 1,iw,iconc,3, 6),iw=1,9),iconc=1,7)/    &
0.0467, 0.1095, 0.2425, 0.4335, 0.6864, 0.8889, 0.9811, 0.9970, 0.9970,  &
0.0299, 0.0759, 0.1728, 0.3251, 0.5632, 0.7917, 0.9438, 0.9942, 0.9970,  &
0.0185, 0.0508, 0.1174, 0.2554, 0.4496, 0.6864, 0.8809, 0.9791, 0.9970,  &
0.0122, 0.0328, 0.0819, 0.1834, 0.3551, 0.5792, 0.8032, 0.9483, 0.9942,  &
0.0071, 0.0225, 0.0599, 0.1433, 0.2824, 0.4983, 0.7284, 0.8964, 0.9811,  &
0.0044, 0.0136, 0.0392, 0.1021, 0.2177, 0.4016, 0.6418, 0.8354, 0.9564,  &
0.0024, 0.0079, 0.0248, 0.0702, 0.1625, 0.3251, 0.5471, 0.7548, 0.8964 /
data ((cldnuctab( 1,iw,iconc,4, 6),iw=1,9),iconc=1,7)/    &
0.0359, 0.0883, 0.1944, 0.3859, 0.6264, 0.8453, 0.9666, 0.9970, 0.9970,  &
0.0225, 0.0599, 0.1342, 0.2824, 0.4820, 0.7284, 0.9103, 0.9875, 0.9970,  &
0.0136, 0.0359, 0.0883, 0.2059, 0.3859, 0.6108, 0.8354, 0.9600, 0.9964,  &
0.0079, 0.0248, 0.0599, 0.1433, 0.2963, 0.4983, 0.7418, 0.9103, 0.9861,  &
0.0050, 0.0151, 0.0428, 0.1021, 0.2177, 0.4016, 0.6418, 0.8453, 0.9634,  &
0.0031, 0.0099, 0.0272, 0.0759, 0.1625, 0.3251, 0.5471, 0.7675, 0.9229,  &
0.0018, 0.0056, 0.0185, 0.0508, 0.1256, 0.2554, 0.4658, 0.6864, 0.8548 /
data ((cldnuctab( 1,iw,iconc,5, 6),iw=1,9),iconc=1,7)/    &
0.0328, 0.0759, 0.1728, 0.3400, 0.5792, 0.8143, 0.9525, 0.9959, 0.9970,  &
0.0185, 0.0467, 0.1095, 0.2425, 0.4335, 0.6864, 0.8809, 0.9791, 0.9970,  &
0.0110, 0.0299, 0.0759, 0.1625, 0.3251, 0.5471, 0.7798, 0.9391, 0.9927,  &
0.0063, 0.0185, 0.0467, 0.1174, 0.2425, 0.4335, 0.6718, 0.8726, 0.9747,  &
0.0039, 0.0110, 0.0328, 0.0819, 0.1728, 0.3400, 0.5632, 0.7917, 0.9391,  &
0.0024, 0.0071, 0.0204, 0.0552, 0.1256, 0.2687, 0.4658, 0.7007, 0.8809,  &
0.0014, 0.0044, 0.0136, 0.0359, 0.0950, 0.2059, 0.3859, 0.6108, 0.8143 /
data ((cldnuctab( 1,iw,iconc,6, 6),iw=1,9),iconc=1,7)/    &
0.0272, 0.0702, 0.1527, 0.3106, 0.5471, 0.7798, 0.9438, 0.9942, 0.9970,  &
0.0151, 0.0428, 0.1021, 0.2177, 0.4016, 0.6418, 0.8639, 0.9722, 0.9970,  &
0.0088, 0.0248, 0.0649, 0.1433, 0.2963, 0.4983, 0.7418, 0.9229, 0.9887,  &
0.0050, 0.0151, 0.0392, 0.0950, 0.2059, 0.3859, 0.6264, 0.8354, 0.9600,  &
0.0031, 0.0088, 0.0248, 0.0649, 0.1433, 0.2963, 0.4983, 0.7418, 0.9103,  &
0.0018, 0.0056, 0.0167, 0.0428, 0.1021, 0.2177, 0.4016, 0.6264, 0.8354,  &
0.0009, 0.0035, 0.0099, 0.0272, 0.0702, 0.1625, 0.3251, 0.5308, 0.7548 /
data ((cldnuctab( 1,iw,iconc,7, 6),iw=1,9),iconc=1,7)/    &
0.0248, 0.0649, 0.1433, 0.2963, 0.5146, 0.7675, 0.9340, 0.9919, 0.9970,  &
0.0136, 0.0359, 0.0883, 0.1944, 0.3704, 0.6108, 0.8453, 0.9666, 0.9970,  &
0.0079, 0.0225, 0.0552, 0.1256, 0.2687, 0.4658, 0.7147, 0.9036, 0.9846,  &
0.0044, 0.0122, 0.0328, 0.0819, 0.1834, 0.3551, 0.5792, 0.8032, 0.9483,  &
0.0027, 0.0079, 0.0225, 0.0552, 0.1256, 0.2554, 0.4658, 0.7007, 0.8889,  &
0.0014, 0.0044, 0.0136, 0.0359, 0.0883, 0.1944, 0.3551, 0.5792, 0.8032,  &
0.0008, 0.0027, 0.0079, 0.0225, 0.0599, 0.1342, 0.2824, 0.4820, 0.7147 /
data ((cldnuctab( 1,iw,iconc,1, 7),iw=1,9),iconc=1,7)/    &
0.1021, 0.2177, 0.4016, 0.6418, 0.8548, 0.9695, 0.9970, 0.9970, 0.9970,  &
0.0702, 0.1625, 0.3106, 0.5308, 0.7675, 0.9286, 0.9909, 0.9970, 0.9970,  &
0.0467, 0.1174, 0.2425, 0.4496, 0.6718, 0.8726, 0.9722, 0.9970, 0.9970,  &
0.0328, 0.0819, 0.1944, 0.3704, 0.5951, 0.8032, 0.9438, 0.9927, 0.9970,  &
0.0225, 0.0599, 0.1433, 0.2963, 0.4983, 0.7284, 0.8964, 0.9791, 0.9970,  &
0.0122, 0.0359, 0.1021, 0.2177, 0.4175, 0.6264, 0.8143, 0.9391, 0.9919,  &
0.0056, 0.0185, 0.0552, 0.1433, 0.2554, 0.3859, 0.4820, 0.5792, 0.8250 /
data ((cldnuctab( 1,iw,iconc,2, 7),iw=1,9),iconc=1,7)/    &
0.0702, 0.1625, 0.3251, 0.5471, 0.7798, 0.9391, 0.9935, 0.9970, 0.9970,  &
0.0467, 0.1174, 0.2425, 0.4335, 0.6718, 0.8726, 0.9770, 0.9970, 0.9970,  &
0.0328, 0.0819, 0.1834, 0.3551, 0.5792, 0.8032, 0.9438, 0.9935, 0.9970,  &
0.0204, 0.0552, 0.1342, 0.2824, 0.4820, 0.7147, 0.8964, 0.9791, 0.9970,  &
0.0136, 0.0392, 0.0950, 0.2177, 0.4016, 0.6264, 0.8354, 0.9564, 0.9948,  &
0.0079, 0.0248, 0.0702, 0.1625, 0.3251, 0.5471, 0.7675, 0.9103, 0.9829,  &
0.0039, 0.0122, 0.0392, 0.1095, 0.2425, 0.4175, 0.6264, 0.7917, 0.9229 /
data ((cldnuctab( 1,iw,iconc,3, 7),iw=1,9),iconc=1,7)/    &
0.0552, 0.1256, 0.2687, 0.4658, 0.7147, 0.9036, 0.9846, 0.9970, 0.9970,  &
0.0328, 0.0883, 0.1944, 0.3704, 0.5951, 0.8143, 0.9525, 0.9954, 0.9970,  &
0.0225, 0.0599, 0.1342, 0.2824, 0.4820, 0.7147, 0.9036, 0.9829, 0.9970,  &
0.0136, 0.0392, 0.0950, 0.2059, 0.3859, 0.6264, 0.8354, 0.9564, 0.9954,  &
0.0088, 0.0272, 0.0702, 0.1625, 0.3106, 0.5308, 0.7548, 0.9168, 0.9861,  &
0.0056, 0.0167, 0.0467, 0.1174, 0.2425, 0.4496, 0.6718, 0.8639, 0.9634,  &
0.0027, 0.0088, 0.0272, 0.0759, 0.1834, 0.3551, 0.5792, 0.7675, 0.9103 /
data ((cldnuctab( 1,iw,iconc,4, 7),iw=1,9),iconc=1,7)/    &
0.0428, 0.1021, 0.2177, 0.4175, 0.6569, 0.8639, 0.9747, 0.9970, 0.9970,  &
0.0272, 0.0649, 0.1527, 0.3106, 0.5308, 0.7675, 0.9286, 0.9909, 0.9970,  &
0.0167, 0.0428, 0.1021, 0.2177, 0.4175, 0.6418, 0.8548, 0.9666, 0.9970,  &
0.0099, 0.0272, 0.0702, 0.1625, 0.3251, 0.5471, 0.7675, 0.9229, 0.9887,  &
0.0063, 0.0185, 0.0508, 0.1174, 0.2425, 0.4496, 0.6718, 0.8639, 0.9695,  &
0.0039, 0.0122, 0.0328, 0.0819, 0.1944, 0.3551, 0.5792, 0.8032, 0.9391,  &
0.0021, 0.0063, 0.0204, 0.0552, 0.1342, 0.2824, 0.4983, 0.7147, 0.8809 /
data ((cldnuctab( 1,iw,iconc,5, 7),iw=1,9),iconc=1,7)/    &
0.0359, 0.0883, 0.1944, 0.3704, 0.6108, 0.8354, 0.9634, 0.9970, 0.9970,  &
0.0204, 0.0552, 0.1256, 0.2687, 0.4658, 0.7147, 0.9036, 0.9846, 0.9970,  &
0.0122, 0.0328, 0.0819, 0.1834, 0.3551, 0.5792, 0.8143, 0.9483, 0.9948,  &
0.0079, 0.0225, 0.0552, 0.1342, 0.2687, 0.4658, 0.7007, 0.8889, 0.9791,  &
0.0044, 0.0136, 0.0359, 0.0950, 0.1944, 0.3704, 0.5951, 0.8143, 0.9483,  &
0.0027, 0.0088, 0.0248, 0.0649, 0.1433, 0.2963, 0.4983, 0.7284, 0.9036,  &
0.0016, 0.0050, 0.0151, 0.0428, 0.1021, 0.2299, 0.4175, 0.6418, 0.8354 /
data ((cldnuctab( 1,iw,iconc,6, 7),iw=1,9),iconc=1,7)/    &
0.0328, 0.0759, 0.1728, 0.3400, 0.5792, 0.8143, 0.9525, 0.9954, 0.9970,  &
0.0185, 0.0467, 0.1095, 0.2425, 0.4335, 0.6718, 0.8809, 0.9770, 0.9970,  &
0.0110, 0.0299, 0.0702, 0.1625, 0.3251, 0.5471, 0.7798, 0.9340, 0.9919,  &
0.0063, 0.0185, 0.0467, 0.1095, 0.2299, 0.4175, 0.6569, 0.8639, 0.9695,  &
0.0035, 0.0110, 0.0299, 0.0759, 0.1625, 0.3251, 0.5471, 0.7675, 0.9229,  &
0.0021, 0.0063, 0.0185, 0.0508, 0.1174, 0.2425, 0.4335, 0.6718, 0.8639,  &
0.0012, 0.0039, 0.0122, 0.0328, 0.0819, 0.1834, 0.3551, 0.5792, 0.7917 /
data ((cldnuctab( 1,iw,iconc,7, 7),iw=1,9),iconc=1,7)/    &
0.0299, 0.0702, 0.1625, 0.3251, 0.5471, 0.7917, 0.9438, 0.9942, 0.9970,  &
0.0167, 0.0428, 0.1021, 0.2177, 0.4016, 0.6569, 0.8639, 0.9722, 0.9970,  &
0.0088, 0.0248, 0.0649, 0.1433, 0.2963, 0.5146, 0.7418, 0.9229, 0.9887,  &
0.0056, 0.0151, 0.0392, 0.0950, 0.2059, 0.3859, 0.6108, 0.8354, 0.9600,  &
0.0031, 0.0088, 0.0248, 0.0649, 0.1433, 0.2824, 0.4983, 0.7284, 0.9036,  &
0.0018, 0.0056, 0.0167, 0.0428, 0.1021, 0.2177, 0.3859, 0.6264, 0.8250,  &
0.0011, 0.0031, 0.0099, 0.0272, 0.0702, 0.1527, 0.3106, 0.5146, 0.7418 /
data ((cldnuctab( 2,iw,iconc,1, 1),iw=1,9),iconc=1,7)/    &
0.0643, 0.1516, 0.3234, 0.5773, 0.8441, 0.9767, 0.9994, 0.9999, 0.9999,  &
0.0424, 0.1086, 0.2410, 0.4639, 0.7403, 0.9334, 0.9953, 0.9999, 0.9999,  &
0.0269, 0.0752, 0.1715, 0.3686, 0.6246, 0.8629, 0.9809, 0.9995, 0.9999,  &
0.0182, 0.0503, 0.1246, 0.2808, 0.5290, 0.7784, 0.9520, 0.9967, 0.9999,  &
0.0121, 0.0355, 0.0942, 0.2163, 0.4316, 0.6847, 0.9027, 0.9859, 0.9996,  &
0.0070, 0.0222, 0.0643, 0.1614, 0.3382, 0.5773, 0.8019, 0.9433, 0.9959,  &
0.0039, 0.0134, 0.0424, 0.1086, 0.2285, 0.3686, 0.4964, 0.5290, 0.8130 /
data ((cldnuctab( 2,iw,iconc,2, 1),iw=1,9),iconc=1,7)/    &
0.0424, 0.1086, 0.2539, 0.4801, 0.7533, 0.9433, 0.9967, 0.9999, 0.9999,  &
0.0269, 0.0752, 0.1715, 0.3533, 0.6246, 0.8629, 0.9827, 0.9996, 0.9999,  &
0.0165, 0.0462, 0.1246, 0.2672, 0.4964, 0.7660, 0.9478, 0.9963, 0.9999,  &
0.0108, 0.0324, 0.0875, 0.2045, 0.3998, 0.6701, 0.8879, 0.9859, 0.9996,  &
0.0070, 0.0201, 0.0593, 0.1516, 0.3089, 0.5613, 0.8130, 0.9596, 0.9974,  &
0.0044, 0.0134, 0.0424, 0.1086, 0.2410, 0.4639, 0.7269, 0.9096, 0.9873,  &
0.0023, 0.0087, 0.0269, 0.0752, 0.1821, 0.3686, 0.5773, 0.7784, 0.9222 /
data ((cldnuctab( 2,iw,iconc,3, 1),iw=1,9),iconc=1,7)/    &
0.0324, 0.0812, 0.1931, 0.3841, 0.6551, 0.8955, 0.9897, 0.9998, 0.9999,  &
0.0182, 0.0503, 0.1246, 0.2808, 0.5127, 0.7903, 0.9559, 0.9980, 0.9999,  &
0.0108, 0.0324, 0.0875, 0.2045, 0.3998, 0.6701, 0.8955, 0.9873, 0.9997,  &
0.0070, 0.0201, 0.0593, 0.1422, 0.3089, 0.5452, 0.8019, 0.9596, 0.9977,  &
0.0044, 0.0134, 0.0388, 0.1012, 0.2285, 0.4477, 0.7131, 0.9096, 0.9897,  &
0.0027, 0.0087, 0.0269, 0.0752, 0.1715, 0.3533, 0.6090, 0.8441, 0.9691,  &
0.0016, 0.0055, 0.0165, 0.0503, 0.1246, 0.2808, 0.5127, 0.7403, 0.9096 /
data ((cldnuctab( 2,iw,iconc,4, 1),iw=1,9),iconc=1,7)/    &
0.0245, 0.0643, 0.1614, 0.3234, 0.5933, 0.8441, 0.9789, 0.9994, 0.9999,  &
0.0149, 0.0388, 0.1012, 0.2285, 0.4477, 0.7131, 0.9222, 0.9941, 0.9999,  &
0.0087, 0.0245, 0.0643, 0.1516, 0.3234, 0.5773, 0.8342, 0.9719, 0.9988,  &
0.0049, 0.0149, 0.0424, 0.1086, 0.2410, 0.4477, 0.7131, 0.9222, 0.9926,  &
0.0030, 0.0097, 0.0269, 0.0752, 0.1715, 0.3533, 0.6090, 0.8441, 0.9744,  &
0.0018, 0.0055, 0.0182, 0.0503, 0.1246, 0.2672, 0.4964, 0.7533, 0.9334,  &
0.0011, 0.0034, 0.0108, 0.0324, 0.0875, 0.2045, 0.4156, 0.6551, 0.8716 /
data ((cldnuctab( 2,iw,iconc,5, 1),iw=1,9),iconc=1,7)/    &
0.0222, 0.0546, 0.1332, 0.2947, 0.5290, 0.8019, 0.9630, 0.9985, 0.9999,  &
0.0121, 0.0324, 0.0812, 0.1931, 0.3841, 0.6551, 0.8879, 0.9873, 0.9997,  &
0.0062, 0.0182, 0.0503, 0.1246, 0.2672, 0.4964, 0.7660, 0.9478, 0.9967,  &
0.0039, 0.0108, 0.0324, 0.0812, 0.1931, 0.3841, 0.6400, 0.8716, 0.9827,  &
0.0021, 0.0070, 0.0201, 0.0546, 0.1332, 0.2808, 0.5127, 0.7784, 0.9478,  &
0.0012, 0.0039, 0.0121, 0.0355, 0.0942, 0.2163, 0.4156, 0.6701, 0.8879,  &
0.0007, 0.0023, 0.0078, 0.0245, 0.0643, 0.1516, 0.3234, 0.5613, 0.8019 /
data ((cldnuctab( 2,iw,iconc,6, 1),iw=1,9),iconc=1,7)/    &
0.0182, 0.0503, 0.1246, 0.2672, 0.4964, 0.7660, 0.9520, 0.9974, 0.9999,  &
0.0097, 0.0296, 0.0752, 0.1715, 0.3382, 0.6090, 0.8537, 0.9809, 0.9994,  &
0.0055, 0.0165, 0.0424, 0.1086, 0.2410, 0.4477, 0.7131, 0.9279, 0.9941,  &
0.0030, 0.0087, 0.0269, 0.0696, 0.1614, 0.3234, 0.5773, 0.8238, 0.9691,  &
0.0018, 0.0055, 0.0165, 0.0424, 0.1086, 0.2410, 0.4477, 0.7131, 0.9160,  &
0.0009, 0.0030, 0.0097, 0.0269, 0.0752, 0.1715, 0.3382, 0.5933, 0.8342,  &
0.0006, 0.0018, 0.0062, 0.0182, 0.0503, 0.1246, 0.2672, 0.4801, 0.7403 /
data ((cldnuctab( 2,iw,iconc,7, 1),iw=1,9),iconc=1,7)/    &
0.0165, 0.0462, 0.1086, 0.2410, 0.4639, 0.7403, 0.9385, 0.9959, 0.9999,  &
0.0087, 0.0245, 0.0643, 0.1516, 0.3234, 0.5613, 0.8342, 0.9719, 0.9990,  &
0.0049, 0.0134, 0.0388, 0.0942, 0.2163, 0.4156, 0.6847, 0.9027, 0.9908,  &
0.0027, 0.0078, 0.0222, 0.0593, 0.1422, 0.2947, 0.5290, 0.7903, 0.9559,  &
0.0014, 0.0044, 0.0134, 0.0355, 0.0875, 0.2045, 0.3998, 0.6551, 0.8800,  &
0.0008, 0.0027, 0.0078, 0.0222, 0.0593, 0.1422, 0.2947, 0.5290, 0.7784,  &
0.0004, 0.0016, 0.0049, 0.0149, 0.0388, 0.1012, 0.2163, 0.4156, 0.6701 /
data ((cldnuctab( 2,iw,iconc,1, 2),iw=1,9),iconc=1,7)/    &
0.0942, 0.2045, 0.3998, 0.6551, 0.8800, 0.9844, 0.9996, 0.9999, 0.9999,  &
0.0643, 0.1516, 0.3089, 0.5452, 0.7903, 0.9520, 0.9967, 0.9999, 0.9999,  &
0.0424, 0.1086, 0.2410, 0.4477, 0.6991, 0.9027, 0.9873, 0.9996, 0.9999,  &
0.0296, 0.0752, 0.1821, 0.3533, 0.6090, 0.8342, 0.9662, 0.9977, 0.9999,  &
0.0182, 0.0546, 0.1332, 0.2808, 0.5127, 0.7533, 0.9279, 0.9908, 0.9997,  &
0.0108, 0.0324, 0.0942, 0.2163, 0.4156, 0.6400, 0.8441, 0.9559, 0.9971,  &
0.0049, 0.0182, 0.0546, 0.1332, 0.2539, 0.3841, 0.4801, 0.5773, 0.8537 /
data ((cldnuctab( 2,iw,iconc,2, 2),iw=1,9),iconc=1,7)/    &
0.0643, 0.1516, 0.3234, 0.5613, 0.8019, 0.9596, 0.9977, 0.9999, 0.9999,  &
0.0424, 0.1086, 0.2285, 0.4316, 0.6991, 0.9027, 0.9873, 0.9997, 0.9999,  &
0.0269, 0.0752, 0.1715, 0.3382, 0.5773, 0.8238, 0.9630, 0.9977, 0.9999,  &
0.0182, 0.0503, 0.1246, 0.2672, 0.4801, 0.7269, 0.9160, 0.9897, 0.9997,  &
0.0121, 0.0355, 0.0875, 0.2045, 0.3998, 0.6400, 0.8629, 0.9719, 0.9982,  &
0.0070, 0.0222, 0.0643, 0.1516, 0.3234, 0.5452, 0.7903, 0.9334, 0.9917,  &
0.0034, 0.0121, 0.0388, 0.1012, 0.2285, 0.4156, 0.6246, 0.8130, 0.9385 /
data ((cldnuctab( 2,iw,iconc,3, 2),iw=1,9),iconc=1,7)/    &
0.0503, 0.1164, 0.2539, 0.4639, 0.7269, 0.9279, 0.9934, 0.9999, 0.9999,  &
0.0296, 0.0752, 0.1821, 0.3533, 0.5933, 0.8441, 0.9719, 0.9986, 0.9999,  &
0.0182, 0.0503, 0.1246, 0.2672, 0.4801, 0.7403, 0.9222, 0.9917, 0.9998,  &
0.0121, 0.0355, 0.0875, 0.1931, 0.3841, 0.6246, 0.8537, 0.9719, 0.9985,  &
0.0078, 0.0222, 0.0593, 0.1422, 0.3089, 0.5290, 0.7784, 0.9385, 0.9934,  &
0.0049, 0.0149, 0.0424, 0.1086, 0.2410, 0.4477, 0.6847, 0.8879, 0.9789,  &
0.0023, 0.0078, 0.0245, 0.0696, 0.1715, 0.3533, 0.5773, 0.7903, 0.9334 /
data ((cldnuctab( 2,iw,iconc,4, 2),iw=1,9),iconc=1,7)/    &
0.0388, 0.0942, 0.2045, 0.3998, 0.6551, 0.8879, 0.9859, 0.9996, 0.9999,  &
0.0222, 0.0593, 0.1422, 0.2947, 0.5290, 0.7784, 0.9478, 0.9959, 0.9999,  &
0.0134, 0.0388, 0.0942, 0.2163, 0.3998, 0.6551, 0.8716, 0.9809, 0.9993,  &
0.0087, 0.0245, 0.0643, 0.1516, 0.3089, 0.5452, 0.7784, 0.9433, 0.9953,  &
0.0055, 0.0165, 0.0424, 0.1086, 0.2410, 0.4316, 0.6847, 0.8879, 0.9827,  &
0.0034, 0.0097, 0.0296, 0.0752, 0.1715, 0.3533, 0.5933, 0.8130, 0.9559,  &
0.0018, 0.0062, 0.0182, 0.0503, 0.1246, 0.2808, 0.4964, 0.7269, 0.9027 /
data ((cldnuctab( 2,iw,iconc,5, 2),iw=1,9),iconc=1,7)/    &
0.0324, 0.0812, 0.1821, 0.3533, 0.6090, 0.8537, 0.9767, 0.9991, 0.9999,  &
0.0182, 0.0503, 0.1164, 0.2539, 0.4639, 0.7269, 0.9222, 0.9917, 0.9998,  &
0.0108, 0.0296, 0.0752, 0.1715, 0.3382, 0.5933, 0.8238, 0.9662, 0.9980,  &
0.0062, 0.0182, 0.0503, 0.1164, 0.2539, 0.4639, 0.7131, 0.9096, 0.9886,  &
0.0039, 0.0121, 0.0324, 0.0812, 0.1821, 0.3686, 0.6090, 0.8342, 0.9630,  &
0.0023, 0.0070, 0.0222, 0.0546, 0.1332, 0.2808, 0.4964, 0.7403, 0.9222,  &
0.0014, 0.0044, 0.0134, 0.0388, 0.0942, 0.2163, 0.3998, 0.6551, 0.8537 /
data ((cldnuctab( 2,iw,iconc,6, 2),iw=1,9),iconc=1,7)/    &
0.0269, 0.0696, 0.1614, 0.3234, 0.5773, 0.8238, 0.9662, 0.9982, 0.9999,  &
0.0165, 0.0424, 0.1012, 0.2285, 0.4316, 0.6847, 0.8955, 0.9873, 0.9996,  &
0.0087, 0.0245, 0.0643, 0.1516, 0.3089, 0.5290, 0.7903, 0.9478, 0.9963,  &
0.0055, 0.0149, 0.0388, 0.1012, 0.2163, 0.4156, 0.6551, 0.8800, 0.9809,  &
0.0030, 0.0097, 0.0269, 0.0643, 0.1516, 0.3089, 0.5290, 0.7784, 0.9433,  &
0.0018, 0.0055, 0.0165, 0.0462, 0.1086, 0.2285, 0.4316, 0.6701, 0.8800,  &
0.0011, 0.0034, 0.0097, 0.0296, 0.0752, 0.1715, 0.3382, 0.5773, 0.8019 /
data ((cldnuctab( 2,iw,iconc,7, 2),iw=1,9),iconc=1,7)/    &
0.0269, 0.0643, 0.1516, 0.3089, 0.5452, 0.8019, 0.9596, 0.9974, 0.9999,  &
0.0149, 0.0388, 0.0942, 0.2045, 0.3998, 0.6551, 0.8800, 0.9827, 0.9994,  &
0.0078, 0.0222, 0.0593, 0.1332, 0.2808, 0.4964, 0.7533, 0.9334, 0.9941,  &
0.0044, 0.0134, 0.0355, 0.0875, 0.1931, 0.3686, 0.6090, 0.8441, 0.9719,  &
0.0027, 0.0078, 0.0222, 0.0546, 0.1332, 0.2672, 0.4801, 0.7269, 0.9222,  &
0.0016, 0.0049, 0.0134, 0.0355, 0.0875, 0.1931, 0.3841, 0.6090, 0.8441,  &
0.0009, 0.0030, 0.0087, 0.0245, 0.0593, 0.1422, 0.2947, 0.5127, 0.7533 /
data ((cldnuctab( 2,iw,iconc,1, 3),iw=1,9),iconc=1,7)/    &
0.1332, 0.2808, 0.4964, 0.7403, 0.9222, 0.9908, 0.9997, 0.9999, 0.9999,  &
0.0942, 0.2163, 0.3998, 0.6400, 0.8537, 0.9691, 0.9980, 0.9999, 0.9999,  &
0.0696, 0.1614, 0.3234, 0.5452, 0.7784, 0.9334, 0.9917, 0.9997, 0.9999,  &
0.0462, 0.1164, 0.2539, 0.4477, 0.6991, 0.8879, 0.9789, 0.9986, 0.9999,  &
0.0296, 0.0812, 0.1931, 0.3686, 0.6090, 0.8238, 0.9520, 0.9941, 0.9998,  &
0.0149, 0.0462, 0.1246, 0.2808, 0.4964, 0.7131, 0.8800, 0.9719, 0.9985,  &
0.0062, 0.0222, 0.0696, 0.1614, 0.2947, 0.3998, 0.5127, 0.5127, 0.8537 /
data ((cldnuctab( 2,iw,iconc,2, 3),iw=1,9),iconc=1,7)/    &
0.1012, 0.2163, 0.3998, 0.6551, 0.8629, 0.9767, 0.9988, 0.9999, 0.9999,  &
0.0696, 0.1516, 0.3089, 0.5452, 0.7784, 0.9385, 0.9926, 0.9998, 0.9999,  &
0.0462, 0.1086, 0.2410, 0.4316, 0.6847, 0.8800, 0.9767, 0.9986, 0.9999,  &
0.0296, 0.0812, 0.1821, 0.3533, 0.5773, 0.8130, 0.9478, 0.9941, 0.9998,  &
0.0201, 0.0546, 0.1332, 0.2808, 0.4964, 0.7269, 0.9027, 0.9844, 0.9991,  &
0.0097, 0.0324, 0.0875, 0.2163, 0.3998, 0.6400, 0.8441, 0.9559, 0.9953,  &
0.0044, 0.0149, 0.0503, 0.1332, 0.2947, 0.4801, 0.6991, 0.8441, 0.9596 /
data ((cldnuctab( 2,iw,iconc,3, 3),iw=1,9),iconc=1,7)/    &
0.0752, 0.1715, 0.3382, 0.5613, 0.8019, 0.9520, 0.9963, 0.9999, 0.9999,  &
0.0503, 0.1164, 0.2410, 0.4477, 0.6991, 0.8955, 0.9827, 0.9992, 0.9999,  &
0.0324, 0.0812, 0.1821, 0.3533, 0.5773, 0.8130, 0.9520, 0.9953, 0.9999,  &
0.0201, 0.0546, 0.1332, 0.2808, 0.4801, 0.7269, 0.9027, 0.9844, 0.9992,  &
0.0134, 0.0388, 0.0942, 0.2163, 0.3998, 0.6246, 0.8441, 0.9630, 0.9963,  &
0.0070, 0.0222, 0.0643, 0.1516, 0.3234, 0.5452, 0.7660, 0.9279, 0.9873,  &
0.0034, 0.0108, 0.0355, 0.1012, 0.2285, 0.4316, 0.6701, 0.8537, 0.9520 /
data ((cldnuctab( 2,iw,iconc,4, 3),iw=1,9),iconc=1,7)/    &
0.0593, 0.1332, 0.2808, 0.4964, 0.7533, 0.9279, 0.9917, 0.9998, 0.9999,  &
0.0355, 0.0875, 0.2045, 0.3841, 0.6246, 0.8441, 0.9691, 0.9977, 0.9999,  &
0.0222, 0.0593, 0.1422, 0.2808, 0.4964, 0.7403, 0.9222, 0.9886, 0.9996,  &
0.0149, 0.0388, 0.1012, 0.2163, 0.3998, 0.6400, 0.8537, 0.9662, 0.9974,  &
0.0097, 0.0269, 0.0696, 0.1614, 0.3234, 0.5452, 0.7660, 0.9279, 0.9897,  &
0.0055, 0.0165, 0.0462, 0.1164, 0.2410, 0.4477, 0.6847, 0.8716, 0.9744,  &
0.0027, 0.0087, 0.0245, 0.0752, 0.1821, 0.3533, 0.5933, 0.8019, 0.9334 /
data ((cldnuctab( 2,iw,iconc,5, 3),iw=1,9),iconc=1,7)/    &
0.0503, 0.1164, 0.2410, 0.4477, 0.6991, 0.9027, 0.9859, 0.9995, 0.9999,  &
0.0296, 0.0752, 0.1715, 0.3382, 0.5613, 0.8019, 0.9520, 0.9959, 0.9999,  &
0.0182, 0.0462, 0.1164, 0.2410, 0.4316, 0.6847, 0.8879, 0.9809, 0.9990,  &
0.0108, 0.0296, 0.0752, 0.1715, 0.3382, 0.5613, 0.7903, 0.9433, 0.9941,  &
0.0070, 0.0201, 0.0546, 0.1246, 0.2539, 0.4639, 0.6991, 0.8879, 0.9789,  &
0.0039, 0.0121, 0.0355, 0.0875, 0.1931, 0.3686, 0.5933, 0.8130, 0.9520,  &
0.0023, 0.0062, 0.0201, 0.0546, 0.1422, 0.2947, 0.5127, 0.7403, 0.9027 /
data ((cldnuctab( 2,iw,iconc,6, 3),iw=1,9),iconc=1,7)/    &
0.0424, 0.1012, 0.2163, 0.4156, 0.6701, 0.8800, 0.9809, 0.9991, 0.9999,  &
0.0245, 0.0643, 0.1422, 0.2947, 0.5127, 0.7660, 0.9334, 0.9934, 0.9998,  &
0.0149, 0.0388, 0.0942, 0.2045, 0.3998, 0.6400, 0.8537, 0.9691, 0.9980,  &
0.0087, 0.0245, 0.0643, 0.1422, 0.2947, 0.5127, 0.7403, 0.9222, 0.9897,  &
0.0055, 0.0149, 0.0424, 0.1012, 0.2163, 0.3998, 0.6400, 0.8441, 0.9662,  &
0.0034, 0.0097, 0.0269, 0.0696, 0.1614, 0.3089, 0.5290, 0.7660, 0.9222,  &
0.0016, 0.0055, 0.0165, 0.0424, 0.1086, 0.2410, 0.4316, 0.6701, 0.8629 /
data ((cldnuctab( 2,iw,iconc,7, 3),iw=1,9),iconc=1,7)/    &
0.0388, 0.0942, 0.2045, 0.3998, 0.6400, 0.8629, 0.9767, 0.9986, 0.9999,  &
0.0222, 0.0593, 0.1332, 0.2672, 0.4964, 0.7403, 0.9222, 0.9908, 0.9996,  &
0.0134, 0.0355, 0.0875, 0.1821, 0.3533, 0.5933, 0.8238, 0.9630, 0.9967,  &
0.0078, 0.0222, 0.0546, 0.1246, 0.2672, 0.4639, 0.7131, 0.9027, 0.9844,  &
0.0049, 0.0134, 0.0355, 0.0875, 0.1931, 0.3533, 0.5933, 0.8130, 0.9520,  &
0.0027, 0.0078, 0.0222, 0.0593, 0.1332, 0.2672, 0.4801, 0.7131, 0.8955,  &
0.0014, 0.0044, 0.0134, 0.0355, 0.0942, 0.2045, 0.3841, 0.6090, 0.8238 /
data ((cldnuctab( 2,iw,iconc,1, 4),iw=1,9),iconc=1,7)/    &
0.2045, 0.3841, 0.6090, 0.8238, 0.9559, 0.9953, 0.9998, 0.9999, 0.9999,  &
0.1516, 0.2947, 0.5127, 0.7403, 0.9096, 0.9844, 0.9991, 0.9999, 0.9999,  &
0.1086, 0.2285, 0.4156, 0.6551, 0.8537, 0.9630, 0.9959, 0.9998, 0.9999,  &
0.0696, 0.1715, 0.3382, 0.5613, 0.7784, 0.9279, 0.9886, 0.9993, 0.9999,  &
0.0388, 0.1086, 0.2539, 0.4639, 0.6991, 0.8879, 0.9719, 0.9971, 0.9999,  &
0.0182, 0.0593, 0.1614, 0.3533, 0.5773, 0.7784, 0.9279, 0.9844, 0.9992,  &
0.0070, 0.0245, 0.0812, 0.1821, 0.2947, 0.3686, 0.4316, 0.4316, 0.8955 /
data ((cldnuctab( 2,iw,iconc,2, 4),iw=1,9),iconc=1,7)/    &
0.1422, 0.2947, 0.5127, 0.7403, 0.9160, 0.9873, 0.9994, 0.9999, 0.9999,  &
0.1012, 0.2285, 0.4156, 0.6400, 0.8537, 0.9630, 0.9967, 0.9999, 0.9999,  &
0.0752, 0.1715, 0.3234, 0.5452, 0.7784, 0.9279, 0.9886, 0.9994, 0.9999,  &
0.0503, 0.1246, 0.2539, 0.4639, 0.6847, 0.8716, 0.9719, 0.9971, 0.9999,  &
0.0269, 0.0812, 0.1931, 0.3686, 0.6090, 0.8130, 0.9433, 0.9917, 0.9995,  &
0.0134, 0.0424, 0.1164, 0.2808, 0.4964, 0.7403, 0.9027, 0.9767, 0.9974,  &
0.0049, 0.0182, 0.0593, 0.1614, 0.3382, 0.5452, 0.7403, 0.8955, 0.9744 /
data ((cldnuctab( 2,iw,iconc,3, 4),iw=1,9),iconc=1,7)/    &
0.1086, 0.2410, 0.4316, 0.6701, 0.8716, 0.9744, 0.9982, 0.9999, 0.9999,  &
0.0752, 0.1715, 0.3382, 0.5613, 0.7903, 0.9385, 0.9917, 0.9996, 0.9999,  &
0.0503, 0.1246, 0.2539, 0.4639, 0.6847, 0.8800, 0.9744, 0.9977, 0.9999,  &
0.0355, 0.0875, 0.1931, 0.3686, 0.5933, 0.8130, 0.9433, 0.9917, 0.9996,  &
0.0201, 0.0593, 0.1422, 0.2947, 0.5127, 0.7269, 0.9027, 0.9789, 0.9982,  &
0.0097, 0.0296, 0.0875, 0.2163, 0.4156, 0.6400, 0.8441, 0.9559, 0.9934,  &
0.0039, 0.0134, 0.0424, 0.1246, 0.2947, 0.5127, 0.7403, 0.8955, 0.9744 /
data ((cldnuctab( 2,iw,iconc,4, 4),iw=1,9),iconc=1,7)/    &
0.0875, 0.1931, 0.3686, 0.6090, 0.8342, 0.9596, 0.9963, 0.9999, 0.9999,  &
0.0593, 0.1332, 0.2808, 0.4801, 0.7269, 0.9027, 0.9827, 0.9990, 0.9999,  &
0.0388, 0.0942, 0.2045, 0.3841, 0.6090, 0.8238, 0.9559, 0.9948, 0.9998,  &
0.0245, 0.0643, 0.1516, 0.2947, 0.5127, 0.7403, 0.9096, 0.9827, 0.9988,  &
0.0149, 0.0424, 0.1086, 0.2285, 0.4156, 0.6551, 0.8441, 0.9596, 0.9953,  &
0.0078, 0.0222, 0.0643, 0.1614, 0.3382, 0.5613, 0.7784, 0.9279, 0.9859,  &
0.0030, 0.0097, 0.0324, 0.0942, 0.2285, 0.4477, 0.6847, 0.8629, 0.9630 /
data ((cldnuctab( 2,iw,iconc,5, 4),iw=1,9),iconc=1,7)/    &
0.0752, 0.1614, 0.3234, 0.5613, 0.7903, 0.9433, 0.9934, 0.9997, 0.9999,  &
0.0462, 0.1086, 0.2410, 0.4316, 0.6701, 0.8716, 0.9744, 0.9980, 0.9999,  &
0.0296, 0.0752, 0.1715, 0.3234, 0.5452, 0.7784, 0.9334, 0.9897, 0.9995,  &
0.0201, 0.0503, 0.1164, 0.2410, 0.4477, 0.6701, 0.8716, 0.9691, 0.9971,  &
0.0121, 0.0324, 0.0812, 0.1821, 0.3533, 0.5773, 0.7903, 0.9334, 0.9897,  &
0.0062, 0.0182, 0.0503, 0.1246, 0.2672, 0.4801, 0.6991, 0.8879, 0.9744,  &
0.0027, 0.0078, 0.0245, 0.0752, 0.1821, 0.3686, 0.6090, 0.8130, 0.9433 /
data ((cldnuctab( 2,iw,iconc,6, 4),iw=1,9),iconc=1,7)/    &
0.0643, 0.1516, 0.2947, 0.5127, 0.7660, 0.9279, 0.9908, 0.9996, 0.9999,  &
0.0388, 0.0942, 0.2045, 0.3841, 0.6246, 0.8441, 0.9662, 0.9967, 0.9999,  &
0.0245, 0.0643, 0.1422, 0.2808, 0.4964, 0.7269, 0.9096, 0.9859, 0.9991,  &
0.0165, 0.0424, 0.1012, 0.2045, 0.3841, 0.6246, 0.8342, 0.9559, 0.9948,  &
0.0097, 0.0269, 0.0643, 0.1516, 0.2947, 0.5127, 0.7403, 0.9096, 0.9827,  &
0.0049, 0.0149, 0.0424, 0.1012, 0.2285, 0.4156, 0.6400, 0.8441, 0.9559,  &
0.0021, 0.0070, 0.0201, 0.0593, 0.1516, 0.3234, 0.5452, 0.7660, 0.9160 /
data ((cldnuctab( 2,iw,iconc,7, 4),iw=1,9),iconc=1,7)/    &
0.0593, 0.1332, 0.2808, 0.4964, 0.7403, 0.9160, 0.9886, 0.9994, 0.9999,  &
0.0355, 0.0875, 0.1931, 0.3686, 0.5933, 0.8238, 0.9559, 0.9959, 0.9998,  &
0.0222, 0.0546, 0.1246, 0.2539, 0.4639, 0.6991, 0.8955, 0.9809, 0.9986,  &
0.0134, 0.0355, 0.0875, 0.1821, 0.3533, 0.5773, 0.8019, 0.9433, 0.9926,  &
0.0078, 0.0222, 0.0546, 0.1332, 0.2672, 0.4639, 0.6991, 0.8800, 0.9744,  &
0.0044, 0.0134, 0.0355, 0.0875, 0.1931, 0.3686, 0.5933, 0.8019, 0.9385,  &
0.0018, 0.0062, 0.0182, 0.0503, 0.1246, 0.2672, 0.4801, 0.7131, 0.8879 /
data ((cldnuctab( 2,iw,iconc,1, 5),iw=1,9),iconc=1,7)/    &
0.2410, 0.4477, 0.6701, 0.8716, 0.9691, 0.9971, 0.9999, 0.9999, 0.9999,  &
0.1821, 0.3533, 0.5773, 0.7903, 0.9385, 0.9897, 0.9995, 0.9999, 0.9999,  &
0.1332, 0.2808, 0.4801, 0.7131, 0.8879, 0.9744, 0.9974, 0.9999, 0.9999,  &
0.0875, 0.2045, 0.3998, 0.6246, 0.8342, 0.9520, 0.9926, 0.9996, 0.9999,  &
0.0462, 0.1332, 0.2947, 0.5290, 0.7660, 0.9160, 0.9827, 0.9982, 0.9999,  &
0.0182, 0.0643, 0.1821, 0.3841, 0.6246, 0.8238, 0.9385, 0.9897, 0.9995,  &
0.0070, 0.0269, 0.0812, 0.1931, 0.2947, 0.3533, 0.4316, 0.4477, 0.8716 /
data ((cldnuctab( 2,iw,iconc,2, 5),iw=1,9),iconc=1,7)/    &
0.1821, 0.3533, 0.5773, 0.8019, 0.9433, 0.9917, 0.9996, 0.9999, 0.9999,  &
0.1332, 0.2672, 0.4801, 0.7131, 0.8879, 0.9767, 0.9980, 0.9999, 0.9999,  &
0.0942, 0.2045, 0.3841, 0.6090, 0.8238, 0.9478, 0.9926, 0.9996, 0.9999,  &
0.0643, 0.1516, 0.3089, 0.5290, 0.7533, 0.9096, 0.9809, 0.9982, 0.9999,  &
0.0324, 0.0942, 0.2285, 0.4316, 0.6701, 0.8537, 0.9630, 0.9948, 0.9997,  &
0.0134, 0.0462, 0.1332, 0.3234, 0.5613, 0.7903, 0.9279, 0.9844, 0.9985,  &
0.0055, 0.0201, 0.0643, 0.1821, 0.3686, 0.5773, 0.7784, 0.9027, 0.9827 /
data ((cldnuctab( 2,iw,iconc,3, 5),iw=1,9),iconc=1,7)/    &
0.1422, 0.2808, 0.4964, 0.7269, 0.9096, 0.9844, 0.9990, 0.9999, 0.9999,  &
0.1012, 0.2163, 0.3998, 0.6246, 0.8342, 0.9559, 0.9948, 0.9998, 0.9999,  &
0.0696, 0.1516, 0.3089, 0.5290, 0.7533, 0.9096, 0.9827, 0.9986, 0.9999,  &
0.0462, 0.1086, 0.2410, 0.4316, 0.6551, 0.8537, 0.9630, 0.9953, 0.9997,  &
0.0245, 0.0696, 0.1715, 0.3533, 0.5773, 0.7903, 0.9279, 0.9873, 0.9990,  &
0.0108, 0.0355, 0.1012, 0.2539, 0.4639, 0.7131, 0.8879, 0.9719, 0.9963,  &
0.0044, 0.0149, 0.0503, 0.1422, 0.3234, 0.5613, 0.7784, 0.9222, 0.9827 /
data ((cldnuctab( 2,iw,iconc,4, 5),iw=1,9),iconc=1,7)/    &
0.1164, 0.2410, 0.4316, 0.6701, 0.8716, 0.9744, 0.9977, 0.9999, 0.9999,  &
0.0752, 0.1715, 0.3234, 0.5452, 0.7784, 0.9334, 0.9897, 0.9994, 0.9999,  &
0.0503, 0.1164, 0.2539, 0.4477, 0.6847, 0.8716, 0.9691, 0.9967, 0.9998,  &
0.0324, 0.0812, 0.1821, 0.3533, 0.5773, 0.7903, 0.9334, 0.9886, 0.9993,  &
0.0182, 0.0546, 0.1332, 0.2808, 0.4801, 0.7131, 0.8879, 0.9744, 0.9971,  &
0.0087, 0.0269, 0.0752, 0.1931, 0.3841, 0.6246, 0.8238, 0.9478, 0.9908,  &
0.0034, 0.0108, 0.0355, 0.1086, 0.2672, 0.4964, 0.7269, 0.8955, 0.9744 /
data ((cldnuctab( 2,iw,iconc,5, 5),iw=1,9),iconc=1,7)/    &
0.0942, 0.2045, 0.3841, 0.6246, 0.8342, 0.9630, 0.9963, 0.9998, 0.9999,  &
0.0593, 0.1422, 0.2808, 0.4964, 0.7269, 0.9096, 0.9827, 0.9988, 0.9999,  &
0.0388, 0.0942, 0.2045, 0.3841, 0.6090, 0.8238, 0.9520, 0.9941, 0.9997,  &
0.0269, 0.0643, 0.1516, 0.2947, 0.5127, 0.7403, 0.9027, 0.9809, 0.9982,  &
0.0149, 0.0424, 0.1012, 0.2285, 0.4156, 0.6400, 0.8342, 0.9559, 0.9934,  &
0.0070, 0.0222, 0.0593, 0.1516, 0.3234, 0.5452, 0.7660, 0.9160, 0.9827,  &
0.0027, 0.0087, 0.0296, 0.0812, 0.2163, 0.4316, 0.6701, 0.8629, 0.9596 /
data ((cldnuctab( 2,iw,iconc,6, 5),iw=1,9),iconc=1,7)/    &
0.0812, 0.1821, 0.3533, 0.5773, 0.8130, 0.9520, 0.9948, 0.9997, 0.9999,  &
0.0546, 0.1246, 0.2539, 0.4477, 0.6847, 0.8800, 0.9767, 0.9982, 0.9999,  &
0.0324, 0.0812, 0.1821, 0.3382, 0.5613, 0.7903, 0.9385, 0.9908, 0.9995,  &
0.0222, 0.0546, 0.1246, 0.2539, 0.4477, 0.6847, 0.8716, 0.9719, 0.9971,  &
0.0121, 0.0355, 0.0875, 0.1821, 0.3533, 0.5773, 0.7903, 0.9334, 0.9886,  &
0.0062, 0.0182, 0.0503, 0.1246, 0.2672, 0.4801, 0.6991, 0.8800, 0.9719,  &
0.0023, 0.0078, 0.0245, 0.0696, 0.1715, 0.3686, 0.6090, 0.8130, 0.9433 /
data ((cldnuctab( 2,iw,iconc,7, 5),iw=1,9),iconc=1,7)/    &
0.0752, 0.1715, 0.3234, 0.5613, 0.7903, 0.9433, 0.9934, 0.9997, 0.9999,  &
0.0462, 0.1086, 0.2285, 0.4156, 0.6551, 0.8629, 0.9719, 0.9974, 0.9999,  &
0.0296, 0.0696, 0.1614, 0.3089, 0.5290, 0.7533, 0.9222, 0.9873, 0.9992,  &
0.0182, 0.0462, 0.1086, 0.2285, 0.4156, 0.6400, 0.8441, 0.9630, 0.9953,  &
0.0108, 0.0296, 0.0752, 0.1614, 0.3234, 0.5290, 0.7533, 0.9160, 0.9844,  &
0.0055, 0.0149, 0.0424, 0.1086, 0.2285, 0.4316, 0.6551, 0.8441, 0.9596,  &
0.0021, 0.0070, 0.0201, 0.0593, 0.1516, 0.3234, 0.5452, 0.7660, 0.9160 /
data ((cldnuctab( 2,iw,iconc,1, 6),iw=1,9),iconc=1,7)/    &
0.2808, 0.4964, 0.7131, 0.8955, 0.9767, 0.9980, 0.9999, 0.9999, 0.9999,  &
0.2163, 0.3998, 0.6246, 0.8342, 0.9520, 0.9926, 0.9996, 0.9999, 0.9999,  &
0.1614, 0.3234, 0.5290, 0.7533, 0.9096, 0.9827, 0.9985, 0.9999, 0.9999,  &
0.1012, 0.2285, 0.4477, 0.6701, 0.8629, 0.9630, 0.9948, 0.9997, 0.9999,  &
0.0503, 0.1422, 0.3234, 0.5773, 0.8019, 0.9334, 0.9873, 0.9988, 0.9999,  &
0.0201, 0.0696, 0.1931, 0.3998, 0.6551, 0.8441, 0.9520, 0.9926, 0.9996,  &
0.0078, 0.0269, 0.0875, 0.1931, 0.2808, 0.3382, 0.3382, 0.4639, 0.8879 /
data ((cldnuctab( 2,iw,iconc,2, 6),iw=1,9),iconc=1,7)/    &
0.2163, 0.3998, 0.6246, 0.8342, 0.9559, 0.9948, 0.9997, 0.9999, 0.9999,  &
0.1614, 0.3089, 0.5290, 0.7533, 0.9096, 0.9827, 0.9986, 0.9999, 0.9999,  &
0.1164, 0.2410, 0.4316, 0.6701, 0.8537, 0.9630, 0.9948, 0.9997, 0.9999,  &
0.0696, 0.1715, 0.3533, 0.5773, 0.7903, 0.9279, 0.9873, 0.9988, 0.9999,  &
0.0355, 0.1086, 0.2539, 0.4801, 0.7131, 0.8879, 0.9719, 0.9963, 0.9998,  &
0.0149, 0.0503, 0.1516, 0.3382, 0.5933, 0.8130, 0.9433, 0.9886, 0.9990,  &
0.0055, 0.0201, 0.0696, 0.1931, 0.3841, 0.5933, 0.7784, 0.9222, 0.9873 /
data ((cldnuctab( 2,iw,iconc,3, 6),iw=1,9),iconc=1,7)/    &
0.1614, 0.3234, 0.5452, 0.7784, 0.9279, 0.9886, 0.9993, 0.9999, 0.9999,  &
0.1164, 0.2410, 0.4477, 0.6701, 0.8629, 0.9662, 0.9963, 0.9998, 0.9999,  &
0.0812, 0.1821, 0.3533, 0.5773, 0.7903, 0.9334, 0.9886, 0.9991, 0.9999,  &
0.0546, 0.1332, 0.2672, 0.4801, 0.7131, 0.8800, 0.9719, 0.9967, 0.9998,  &
0.0269, 0.0812, 0.1931, 0.3841, 0.6246, 0.8238, 0.9478, 0.9908, 0.9993,  &
0.0121, 0.0388, 0.1086, 0.2672, 0.5127, 0.7403, 0.9096, 0.9789, 0.9974,  &
0.0044, 0.0149, 0.0503, 0.1516, 0.3533, 0.5933, 0.8130, 0.9334, 0.9859 /
data ((cldnuctab( 2,iw,iconc,4, 6),iw=1,9),iconc=1,7)/    &
0.1332, 0.2672, 0.4801, 0.7131, 0.8955, 0.9809, 0.9985, 0.9999, 0.9999,  &
0.0942, 0.1931, 0.3686, 0.5933, 0.8130, 0.9478, 0.9926, 0.9996, 0.9999,  &
0.0593, 0.1422, 0.2808, 0.4964, 0.7269, 0.8955, 0.9789, 0.9980, 0.9999,  &
0.0388, 0.1012, 0.2163, 0.3998, 0.6246, 0.8342, 0.9520, 0.9926, 0.9995,  &
0.0222, 0.0593, 0.1516, 0.3089, 0.5290, 0.7533, 0.9096, 0.9809, 0.9980,  &
0.0087, 0.0296, 0.0875, 0.2163, 0.4316, 0.6701, 0.8537, 0.9596, 0.9941,  &
0.0034, 0.0121, 0.0388, 0.1164, 0.2808, 0.5290, 0.7660, 0.9160, 0.9809 /
data ((cldnuctab( 2,iw,iconc,5, 6),iw=1,9),iconc=1,7)/    &
0.1086, 0.2285, 0.4316, 0.6701, 0.8716, 0.9719, 0.9974, 0.9999, 0.9999,  &
0.0752, 0.1614, 0.3234, 0.5452, 0.7660, 0.9279, 0.9886, 0.9992, 0.9999,  &
0.0503, 0.1164, 0.2410, 0.4316, 0.6551, 0.8629, 0.9662, 0.9959, 0.9998,  &
0.0324, 0.0812, 0.1715, 0.3382, 0.5613, 0.7784, 0.9222, 0.9859, 0.9990,  &
0.0182, 0.0503, 0.1164, 0.2539, 0.4639, 0.6847, 0.8716, 0.9662, 0.9959,  &
0.0078, 0.0245, 0.0696, 0.1715, 0.3533, 0.5933, 0.8019, 0.9334, 0.9873,  &
0.0030, 0.0097, 0.0296, 0.0875, 0.2285, 0.4639, 0.7131, 0.8879, 0.9691 /
data ((cldnuctab( 2,iw,iconc,6, 6),iw=1,9),iconc=1,7)/    &
0.0942, 0.2045, 0.3998, 0.6246, 0.8441, 0.9630, 0.9963, 0.9998, 0.9999,  &
0.0643, 0.1422, 0.2808, 0.4964, 0.7269, 0.9096, 0.9844, 0.9988, 0.9999,  &
0.0424, 0.0942, 0.2045, 0.3841, 0.6090, 0.8238, 0.9520, 0.9941, 0.9996,  &
0.0269, 0.0643, 0.1422, 0.2947, 0.4964, 0.7269, 0.8955, 0.9789, 0.9980,  &
0.0149, 0.0388, 0.1012, 0.2163, 0.3998, 0.6246, 0.8238, 0.9520, 0.9926,  &
0.0062, 0.0201, 0.0546, 0.1422, 0.2947, 0.5290, 0.7403, 0.9027, 0.9789,  &
0.0023, 0.0078, 0.0245, 0.0752, 0.1931, 0.3998, 0.6400, 0.8441, 0.9559 /
data ((cldnuctab( 2,iw,iconc,7, 6),iw=1,9),iconc=1,7)/    &
0.0875, 0.1931, 0.3686, 0.6090, 0.8238, 0.9596, 0.9953, 0.9998, 0.9999,  &
0.0546, 0.1246, 0.2672, 0.4639, 0.6991, 0.8955, 0.9809, 0.9985, 0.9999,  &
0.0355, 0.0875, 0.1821, 0.3533, 0.5773, 0.8019, 0.9433, 0.9917, 0.9995,  &
0.0222, 0.0546, 0.1246, 0.2672, 0.4639, 0.6847, 0.8716, 0.9719, 0.9971,  &
0.0134, 0.0355, 0.0875, 0.1931, 0.3533, 0.5773, 0.7903, 0.9334, 0.9886,  &
0.0055, 0.0165, 0.0503, 0.1246, 0.2672, 0.4639, 0.6991, 0.8800, 0.9691,  &
0.0021, 0.0070, 0.0222, 0.0643, 0.1614, 0.3533, 0.5933, 0.8019, 0.9385 /
data ((cldnuctab( 2,iw,iconc,1, 7),iw=1,9),iconc=1,7)/    &
0.3089, 0.5290, 0.7533, 0.9096, 0.9827, 0.9986, 0.9999, 0.9999, 0.9999,  &
0.2410, 0.4316, 0.6701, 0.8537, 0.9596, 0.9948, 0.9997, 0.9999, 0.9999,  &
0.1715, 0.3533, 0.5773, 0.7903, 0.9279, 0.9859, 0.9988, 0.9999, 0.9999,  &
0.1086, 0.2539, 0.4801, 0.7131, 0.8800, 0.9719, 0.9963, 0.9998, 0.9999,  &
0.0503, 0.1516, 0.3533, 0.6090, 0.8238, 0.9478, 0.9897, 0.9991, 0.9999,  &
0.0201, 0.0696, 0.2045, 0.4316, 0.6701, 0.8629, 0.9630, 0.9941, 0.9997,  &
0.0078, 0.0269, 0.0875, 0.2045, 0.2808, 0.3234, 0.3382, 0.3533, 0.8955 /
data ((cldnuctab( 2,iw,iconc,2, 7),iw=1,9),iconc=1,7)/    &
0.2410, 0.4316, 0.6701, 0.8629, 0.9662, 0.9959, 0.9998, 0.9999, 0.9999,  &
0.1821, 0.3382, 0.5613, 0.7784, 0.9279, 0.9873, 0.9990, 0.9999, 0.9999,  &
0.1246, 0.2672, 0.4801, 0.6991, 0.8800, 0.9691, 0.9963, 0.9998, 0.9999,  &
0.0812, 0.1931, 0.3841, 0.6090, 0.8130, 0.9433, 0.9897, 0.9992, 0.9999,  &
0.0388, 0.1164, 0.2672, 0.5127, 0.7403, 0.9027, 0.9789, 0.9974, 0.9998,  &
0.0149, 0.0503, 0.1516, 0.3686, 0.6246, 0.8342, 0.9520, 0.9908, 0.9993,  &
0.0055, 0.0201, 0.0696, 0.1931, 0.3998, 0.6090, 0.7903, 0.9334, 0.9908 /
data ((cldnuctab( 2,iw,iconc,3, 7),iw=1,9),iconc=1,7)/    &
0.1821, 0.3533, 0.5773, 0.8019, 0.9433, 0.9917, 0.9995, 0.9999, 0.9999,  &
0.1332, 0.2672, 0.4801, 0.7131, 0.8879, 0.9744, 0.9974, 0.9999, 0.9999,  &
0.0942, 0.2045, 0.3841, 0.6090, 0.8130, 0.9433, 0.9908, 0.9994, 0.9999,  &
0.0593, 0.1422, 0.2947, 0.5127, 0.7403, 0.9027, 0.9789, 0.9974, 0.9998,  &
0.0296, 0.0875, 0.2163, 0.4156, 0.6551, 0.8441, 0.9559, 0.9934, 0.9995,  &
0.0121, 0.0388, 0.1164, 0.2947, 0.5452, 0.7784, 0.9222, 0.9827, 0.9980,  &
0.0044, 0.0149, 0.0546, 0.1614, 0.3686, 0.6246, 0.8238, 0.9433, 0.9897 /
data ((cldnuctab( 2,iw,iconc,4, 7),iw=1,9),iconc=1,7)/    &
0.1516, 0.2947, 0.5127, 0.7533, 0.9160, 0.9859, 0.9990, 0.9999, 0.9999,  &
0.1012, 0.2163, 0.3998, 0.6400, 0.8441, 0.9596, 0.9948, 0.9997, 0.9999,  &
0.0696, 0.1614, 0.3234, 0.5290, 0.7533, 0.9160, 0.9827, 0.9985, 0.9999,  &
0.0462, 0.1086, 0.2410, 0.4316, 0.6701, 0.8537, 0.9596, 0.9948, 0.9996,  &
0.0222, 0.0643, 0.1614, 0.3382, 0.5613, 0.7784, 0.9279, 0.9859, 0.9986,  &
0.0097, 0.0296, 0.0942, 0.2285, 0.4639, 0.6991, 0.8800, 0.9691, 0.9953,  &
0.0034, 0.0121, 0.0388, 0.1246, 0.3089, 0.5613, 0.7903, 0.9334, 0.9844 /
data ((cldnuctab( 2,iw,iconc,5, 7),iw=1,9),iconc=1,7)/    &
0.1246, 0.2539, 0.4639, 0.6991, 0.8879, 0.9789, 0.9982, 0.9999, 0.9999,  &
0.0875, 0.1821, 0.3533, 0.5773, 0.8019, 0.9433, 0.9917, 0.9995, 0.9999,  &
0.0546, 0.1332, 0.2672, 0.4639, 0.6991, 0.8800, 0.9719, 0.9971, 0.9998,  &
0.0355, 0.0875, 0.1931, 0.3686, 0.5933, 0.8019, 0.9385, 0.9897, 0.9992,  &
0.0182, 0.0546, 0.1332, 0.2808, 0.4964, 0.7131, 0.8879, 0.9744, 0.9967,  &
0.0078, 0.0245, 0.0752, 0.1821, 0.3841, 0.6246, 0.8238, 0.9478, 0.9908,  &
0.0030, 0.0097, 0.0324, 0.0942, 0.2410, 0.4801, 0.7403, 0.9027, 0.9767 /
data ((cldnuctab( 2,iw,iconc,6, 7),iw=1,9),iconc=1,7)/    &
0.1086, 0.2285, 0.4316, 0.6701, 0.8716, 0.9719, 0.9974, 0.9999, 0.9999,  &
0.0752, 0.1614, 0.3089, 0.5290, 0.7660, 0.9279, 0.9886, 0.9992, 0.9999,  &
0.0462, 0.1086, 0.2285, 0.4156, 0.6400, 0.8537, 0.9630, 0.9953, 0.9997,  &
0.0296, 0.0752, 0.1614, 0.3234, 0.5290, 0.7533, 0.9160, 0.9844, 0.9986,  &
0.0165, 0.0462, 0.1086, 0.2410, 0.4316, 0.6551, 0.8537, 0.9596, 0.9941,  &
0.0070, 0.0201, 0.0593, 0.1516, 0.3234, 0.5613, 0.7784, 0.9222, 0.9844,  &
0.0027, 0.0078, 0.0269, 0.0752, 0.2045, 0.4156, 0.6701, 0.8629, 0.9630 /
data ((cldnuctab( 2,iw,iconc,7, 7),iw=1,9),iconc=1,7)/    &
0.1012, 0.2163, 0.3998, 0.6400, 0.8537, 0.9662, 0.9967, 0.9998, 0.9999,  &
0.0643, 0.1422, 0.2947, 0.4964, 0.7403, 0.9096, 0.9844, 0.9988, 0.9999,  &
0.0424, 0.1012, 0.2045, 0.3841, 0.6090, 0.8238, 0.9520, 0.9941, 0.9996,  &
0.0269, 0.0643, 0.1422, 0.2808, 0.4964, 0.7269, 0.8955, 0.9789, 0.9980,  &
0.0149, 0.0388, 0.0942, 0.2045, 0.3841, 0.6090, 0.8238, 0.9478, 0.9917,  &
0.0062, 0.0182, 0.0546, 0.1332, 0.2808, 0.4964, 0.7269, 0.8955, 0.9767,  &
0.0021, 0.0070, 0.0222, 0.0643, 0.1715, 0.3686, 0.6246, 0.8342, 0.9478 /
data ((cldnuctab( 3,iw,iconc,1, 1),iw=1,9),iconc=1,7)/    &
0.2107, 0.4002, 0.6519, 0.8697, 0.9795, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.1616, 0.3148, 0.5457, 0.7892, 0.9485, 0.9952, 0.9998, 0.9999, 0.9999,  &
0.1117, 0.2527, 0.4542, 0.7013, 0.8967, 0.9836, 0.9992, 0.9999, 0.9999,  &
0.0739, 0.1852, 0.3652, 0.6173, 0.8383, 0.9614, 0.9968, 0.9999, 0.9999,  &
0.0385, 0.1117, 0.2676, 0.5091, 0.7617, 0.9260, 0.9898, 0.9994, 0.9999,  &
0.0183, 0.0619, 0.1732, 0.3652, 0.6173, 0.8268, 0.9531, 0.9945, 0.9998,  &
0.0070, 0.0254, 0.0806, 0.1852, 0.2986, 0.3481, 0.4181, 0.4361, 0.8882 /
data ((cldnuctab( 3,iw,iconc,2, 1),iw=1,9),iconc=1,7)/    &
0.1507, 0.3148, 0.5457, 0.7892, 0.9531, 0.9958, 0.9999, 0.9999, 0.9999,  &
0.1117, 0.2382, 0.4361, 0.6852, 0.8967, 0.9836, 0.9993, 0.9999, 0.9999,  &
0.0739, 0.1732, 0.3481, 0.5997, 0.8268, 0.9575, 0.9963, 0.9999, 0.9999,  &
0.0515, 0.1302, 0.2676, 0.4908, 0.7324, 0.9193, 0.9884, 0.9994, 0.9999,  &
0.0283, 0.0806, 0.1977, 0.4002, 0.6519, 0.8598, 0.9716, 0.9976, 0.9999,  &
0.0129, 0.0425, 0.1207, 0.2986, 0.5457, 0.7757, 0.9380, 0.9909, 0.9996,  &
0.0054, 0.0183, 0.0619, 0.1732, 0.3481, 0.5638, 0.7757, 0.9122, 0.9898 /
data ((cldnuctab( 3,iw,iconc,3, 1),iw=1,9),iconc=1,7)/    &
0.1207, 0.2527, 0.4542, 0.7170, 0.9122, 0.9898, 0.9997, 0.9999, 0.9999,  &
0.0806, 0.1852, 0.3652, 0.5997, 0.8268, 0.9651, 0.9976, 0.9999, 0.9999,  &
0.0515, 0.1302, 0.2676, 0.4908, 0.7324, 0.9193, 0.9898, 0.9996, 0.9999,  &
0.0348, 0.0952, 0.2107, 0.4002, 0.6347, 0.8493, 0.9685, 0.9979, 0.9999,  &
0.0205, 0.0619, 0.1507, 0.3148, 0.5457, 0.7757, 0.9380, 0.9920, 0.9997,  &
0.0102, 0.0314, 0.0877, 0.2242, 0.4361, 0.6852, 0.8882, 0.9795, 0.9984,  &
0.0041, 0.0129, 0.0425, 0.1302, 0.3148, 0.5457, 0.7757, 0.9322, 0.9898 /
data ((cldnuctab( 3,iw,iconc,4, 1),iw=1,9),iconc=1,7)/    &
0.0952, 0.1977, 0.4002, 0.6519, 0.8697, 0.9795, 0.9991, 0.9999, 0.9999,  &
0.0619, 0.1402, 0.2986, 0.5091, 0.7617, 0.9380, 0.9938, 0.9998, 0.9999,  &
0.0385, 0.0952, 0.2107, 0.4002, 0.6519, 0.8697, 0.9771, 0.9986, 0.9999,  &
0.0254, 0.0677, 0.1616, 0.3148, 0.5457, 0.7892, 0.9380, 0.9938, 0.9997,  &
0.0163, 0.0425, 0.1117, 0.2382, 0.4542, 0.6852, 0.8882, 0.9795, 0.9988,  &
0.0079, 0.0228, 0.0677, 0.1732, 0.3481, 0.5997, 0.8148, 0.9531, 0.9952,  &
0.0031, 0.0102, 0.0348, 0.0952, 0.2382, 0.4725, 0.7170, 0.9047, 0.9795 /
data ((cldnuctab( 3,iw,iconc,5, 1),iw=1,9),iconc=1,7)/    &
0.0739, 0.1732, 0.3481, 0.5818, 0.8268, 0.9685, 0.9982, 0.9999, 0.9999,  &
0.0468, 0.1117, 0.2527, 0.4542, 0.7013, 0.9047, 0.9884, 0.9996, 0.9999,  &
0.0314, 0.0806, 0.1732, 0.3481, 0.5818, 0.8148, 0.9575, 0.9968, 0.9999,  &
0.0205, 0.0515, 0.1207, 0.2527, 0.4725, 0.7170, 0.9047, 0.9854, 0.9993,  &
0.0115, 0.0348, 0.0877, 0.1852, 0.3652, 0.5997, 0.8268, 0.9614, 0.9963,  &
0.0061, 0.0183, 0.0515, 0.1302, 0.2829, 0.5091, 0.7473, 0.9193, 0.9884,  &
0.0027, 0.0079, 0.0254, 0.0739, 0.1977, 0.4002, 0.6519, 0.8598, 0.9651 /
data ((cldnuctab( 3,iw,iconc,6, 1),iw=1,9),iconc=1,7)/    &
0.0677, 0.1507, 0.3148, 0.5457, 0.8022, 0.9575, 0.9968, 0.9999, 0.9999,  &
0.0425, 0.0952, 0.2107, 0.4002, 0.6519, 0.8792, 0.9816, 0.9992, 0.9999,  &
0.0254, 0.0677, 0.1507, 0.2986, 0.5274, 0.7757, 0.9380, 0.9938, 0.9998,  &
0.0163, 0.0425, 0.1032, 0.2242, 0.4181, 0.6519, 0.8697, 0.9745, 0.9986,  &
0.0102, 0.0283, 0.0677, 0.1616, 0.3148, 0.5457, 0.7757, 0.9380, 0.9929,  &
0.0054, 0.0145, 0.0425, 0.1032, 0.2382, 0.4361, 0.6852, 0.8792, 0.9771,  &
0.0020, 0.0070, 0.0205, 0.0619, 0.1616, 0.3312, 0.5638, 0.8022, 0.9434 /
data ((cldnuctab( 3,iw,iconc,7, 1),iw=1,9),iconc=1,7)/    &
0.0619, 0.1402, 0.2829, 0.5091, 0.7757, 0.9434, 0.9958, 0.9999, 0.9999,  &
0.0348, 0.0877, 0.1977, 0.3826, 0.6173, 0.8598, 0.9745, 0.9988, 0.9999,  &
0.0228, 0.0565, 0.1302, 0.2676, 0.4908, 0.7324, 0.9193, 0.9909, 0.9997,  &
0.0129, 0.0348, 0.0877, 0.1977, 0.3652, 0.5997, 0.8383, 0.9651, 0.9976,  &
0.0079, 0.0228, 0.0565, 0.1302, 0.2829, 0.4908, 0.7324, 0.9122, 0.9870,  &
0.0047, 0.0129, 0.0348, 0.0877, 0.1977, 0.3826, 0.6173, 0.8383, 0.9614,  &
0.0020, 0.0061, 0.0183, 0.0515, 0.1302, 0.2829, 0.5091, 0.7473, 0.9193 /
data ((cldnuctab( 3,iw,iconc,1, 2),iw=1,9),iconc=1,7)/    &
0.2829, 0.5091, 0.7473, 0.9193, 0.9870, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.2242, 0.4181, 0.6519, 0.8598, 0.9685, 0.9972, 0.9999, 0.9999, 0.9999,  &
0.1616, 0.3312, 0.5638, 0.7757, 0.9322, 0.9909, 0.9995, 0.9999, 0.9999,  &
0.1032, 0.2382, 0.4542, 0.7013, 0.8882, 0.9771, 0.9982, 0.9999, 0.9999,  &
0.0468, 0.1402, 0.3312, 0.5997, 0.8268, 0.9531, 0.9938, 0.9997, 0.9999,  &
0.0205, 0.0677, 0.1977, 0.4181, 0.6687, 0.8697, 0.9685, 0.9968, 0.9999,  &
0.0070, 0.0283, 0.0877, 0.1977, 0.2829, 0.3312, 0.3312, 0.4542, 0.9047 /
data ((cldnuctab( 3,iw,iconc,2, 2),iw=1,9),iconc=1,7)/    &
0.2242, 0.4181, 0.6519, 0.8598, 0.9716, 0.9979, 0.9999, 0.9999, 0.9999,  &
0.1616, 0.3148, 0.5457, 0.7757, 0.9322, 0.9909, 0.9996, 0.9999, 0.9999,  &
0.1117, 0.2527, 0.4542, 0.6852, 0.8792, 0.9745, 0.9979, 0.9999, 0.9999,  &
0.0739, 0.1732, 0.3652, 0.5997, 0.8148, 0.9485, 0.9929, 0.9997, 0.9999,  &
0.0348, 0.1032, 0.2527, 0.4908, 0.7324, 0.9047, 0.9836, 0.9986, 0.9999,  &
0.0145, 0.0515, 0.1507, 0.3481, 0.6173, 0.8383, 0.9575, 0.9945, 0.9997,  &
0.0054, 0.0205, 0.0677, 0.1977, 0.3826, 0.5997, 0.7892, 0.9380, 0.9938 /
data ((cldnuctab( 3,iw,iconc,3, 2),iw=1,9),iconc=1,7)/    &
0.1732, 0.3312, 0.5638, 0.7892, 0.9434, 0.9945, 0.9998, 0.9999, 0.9999,  &
0.1207, 0.2527, 0.4542, 0.6852, 0.8882, 0.9795, 0.9986, 0.9999, 0.9999,  &
0.0806, 0.1852, 0.3652, 0.5997, 0.8148, 0.9485, 0.9938, 0.9997, 0.9999,  &
0.0515, 0.1302, 0.2829, 0.4908, 0.7324, 0.9047, 0.9816, 0.9988, 0.9999,  &
0.0283, 0.0806, 0.1977, 0.4002, 0.6347, 0.8493, 0.9614, 0.9958, 0.9998,  &
0.0115, 0.0385, 0.1117, 0.2829, 0.5274, 0.7617, 0.9260, 0.9870, 0.9991,  &
0.0041, 0.0145, 0.0515, 0.1507, 0.3481, 0.6173, 0.8268, 0.9485, 0.9929 /
data ((cldnuctab( 3,iw,iconc,4, 2),iw=1,9),iconc=1,7)/    &
0.1302, 0.2829, 0.4908, 0.7324, 0.9193, 0.9884, 0.9995, 0.9999, 0.9999,  &
0.0952, 0.1977, 0.3826, 0.6173, 0.8383, 0.9614, 0.9968, 0.9999, 0.9999,  &
0.0619, 0.1402, 0.2986, 0.5091, 0.7473, 0.9122, 0.9870, 0.9992, 0.9999,  &
0.0385, 0.1032, 0.2242, 0.4181, 0.6519, 0.8493, 0.9651, 0.9963, 0.9998,  &
0.0205, 0.0619, 0.1507, 0.3148, 0.5457, 0.7757, 0.9260, 0.9884, 0.9993,  &
0.0090, 0.0283, 0.0877, 0.2242, 0.4361, 0.6852, 0.8792, 0.9716, 0.9972,  &
0.0036, 0.0115, 0.0385, 0.1207, 0.2986, 0.5457, 0.7892, 0.9380, 0.9884 /
data ((cldnuctab( 3,iw,iconc,5, 2),iw=1,9),iconc=1,7)/    &
0.1117, 0.2382, 0.4361, 0.6852, 0.8882, 0.9816, 0.9990, 0.9999, 0.9999,  &
0.0739, 0.1616, 0.3312, 0.5457, 0.7892, 0.9434, 0.9938, 0.9997, 0.9999,  &
0.0515, 0.1207, 0.2382, 0.4361, 0.6852, 0.8792, 0.9771, 0.9984, 0.9999,  &
0.0314, 0.0806, 0.1732, 0.3481, 0.5638, 0.7892, 0.9380, 0.9920, 0.9996,  &
0.0183, 0.0468, 0.1207, 0.2676, 0.4725, 0.7013, 0.8882, 0.9771, 0.9982,  &
0.0079, 0.0228, 0.0677, 0.1732, 0.3652, 0.5997, 0.8148, 0.9485, 0.9929,  &
0.0027, 0.0090, 0.0314, 0.0877, 0.2382, 0.4725, 0.7170, 0.9047, 0.9795 /
data ((cldnuctab( 3,iw,iconc,6, 2),iw=1,9),iconc=1,7)/    &
0.0952, 0.2107, 0.4002, 0.6347, 0.8598, 0.9745, 0.9984, 0.9999, 0.9999,  &
0.0619, 0.1402, 0.2986, 0.5091, 0.7473, 0.9260, 0.9898, 0.9996, 0.9999,  &
0.0425, 0.0952, 0.2107, 0.4002, 0.6173, 0.8383, 0.9651, 0.9968, 0.9999,  &
0.0254, 0.0677, 0.1507, 0.2986, 0.5091, 0.7473, 0.9122, 0.9870, 0.9992,  &
0.0145, 0.0425, 0.1032, 0.2242, 0.4002, 0.6347, 0.8493, 0.9614, 0.9963,  &
0.0061, 0.0205, 0.0565, 0.1402, 0.3148, 0.5274, 0.7617, 0.9193, 0.9870,  &
0.0024, 0.0079, 0.0254, 0.0739, 0.1977, 0.4002, 0.6519, 0.8598, 0.9651 /
data ((cldnuctab( 3,iw,iconc,7, 2),iw=1,9),iconc=1,7)/    &
0.0877, 0.1977, 0.3652, 0.6173, 0.8383, 0.9685, 0.9979, 0.9999, 0.9999,  &
0.0565, 0.1302, 0.2676, 0.4725, 0.7170, 0.9047, 0.9870, 0.9993, 0.9999,  &
0.0348, 0.0877, 0.1852, 0.3481, 0.5818, 0.8148, 0.9531, 0.9952, 0.9998,  &
0.0228, 0.0565, 0.1302, 0.2676, 0.4725, 0.7013, 0.8882, 0.9795, 0.9988,  &
0.0129, 0.0348, 0.0877, 0.1852, 0.3652, 0.5818, 0.8022, 0.9485, 0.9938,  &
0.0061, 0.0183, 0.0468, 0.1207, 0.2676, 0.4725, 0.7170, 0.8967, 0.9795,  &
0.0024, 0.0070, 0.0205, 0.0619, 0.1616, 0.3481, 0.5997, 0.8268, 0.9485 /
data ((cldnuctab( 3,iw,iconc,1, 3),iw=1,9),iconc=1,7)/    &
0.3826, 0.6173, 0.8268, 0.9531, 0.9938, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.2986, 0.5274, 0.7473, 0.9122, 0.9836, 0.9986, 0.9999, 0.9999, 0.9999,  &
0.2107, 0.4181, 0.6519, 0.8598, 0.9614, 0.9952, 0.9997, 0.9999, 0.9999,  &
0.1207, 0.2986, 0.5457, 0.7892, 0.9322, 0.9884, 0.9991, 0.9999, 0.9999,  &
0.0565, 0.1732, 0.4002, 0.6852, 0.8792, 0.9716, 0.9968, 0.9998, 0.9999,  &
0.0205, 0.0739, 0.2107, 0.4542, 0.7013, 0.8967, 0.9816, 0.9982, 0.9999,  &
0.0079, 0.0283, 0.0952, 0.1977, 0.2382, 0.2527, 0.3148, 0.3652, 0.9193 /
data ((cldnuctab( 3,iw,iconc,2, 3),iw=1,9),iconc=1,7)/    &
0.2986, 0.5091, 0.7473, 0.9122, 0.9854, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.2242, 0.4181, 0.6519, 0.8493, 0.9614, 0.9958, 0.9998, 0.9999, 0.9999,  &
0.1616, 0.3312, 0.5638, 0.7757, 0.9260, 0.9870, 0.9991, 0.9999, 0.9999,  &
0.0952, 0.2382, 0.4542, 0.7013, 0.8792, 0.9716, 0.9968, 0.9998, 0.9999,  &
0.0425, 0.1302, 0.3148, 0.5818, 0.8148, 0.9434, 0.9909, 0.9993, 0.9999,  &
0.0163, 0.0565, 0.1732, 0.4002, 0.6852, 0.8882, 0.9745, 0.9972, 0.9998,  &
0.0061, 0.0228, 0.0739, 0.2107, 0.4181, 0.5997, 0.8148, 0.9575, 0.9958 /
data ((cldnuctab( 3,iw,iconc,3, 3),iw=1,9),iconc=1,7)/    &
0.2382, 0.4361, 0.6687, 0.8697, 0.9716, 0.9972, 0.9999, 0.9999, 0.9999,  &
0.1732, 0.3481, 0.5638, 0.7892, 0.9322, 0.9898, 0.9993, 0.9999, 0.9999,  &
0.1207, 0.2676, 0.4725, 0.7013, 0.8792, 0.9716, 0.9972, 0.9998, 0.9999,  &
0.0739, 0.1852, 0.3652, 0.5997, 0.8148, 0.9434, 0.9909, 0.9994, 0.9999,  &
0.0314, 0.0952, 0.2527, 0.4908, 0.7324, 0.9047, 0.9795, 0.9979, 0.9999,  &
0.0129, 0.0425, 0.1302, 0.3312, 0.6173, 0.8383, 0.9575, 0.9938, 0.9996,  &
0.0047, 0.0163, 0.0565, 0.1732, 0.4002, 0.6519, 0.8598, 0.9685, 0.9963 /
data ((cldnuctab( 3,iw,iconc,4, 3),iw=1,9),iconc=1,7)/    &
0.1977, 0.3652, 0.5997, 0.8148, 0.9531, 0.9945, 0.9997, 0.9999, 0.9999,  &
0.1402, 0.2829, 0.4908, 0.7170, 0.8967, 0.9795, 0.9984, 0.9999, 0.9999,  &
0.0952, 0.2107, 0.3826, 0.6173, 0.8268, 0.9531, 0.9938, 0.9996, 0.9999,  &
0.0565, 0.1402, 0.2986, 0.5091, 0.7473, 0.9047, 0.9816, 0.9984, 0.9999,  &
0.0254, 0.0806, 0.1977, 0.4002, 0.6519, 0.8493, 0.9614, 0.9945, 0.9997,  &
0.0102, 0.0348, 0.1032, 0.2676, 0.5274, 0.7757, 0.9260, 0.9854, 0.9986,  &
0.0036, 0.0129, 0.0425, 0.1302, 0.3312, 0.6173, 0.8383, 0.9575, 0.9938 /
data ((cldnuctab( 3,iw,iconc,5, 3),iw=1,9),iconc=1,7)/    &
0.1616, 0.3148, 0.5457, 0.7757, 0.9322, 0.9909, 0.9996, 0.9999, 0.9999,  &
0.1117, 0.2382, 0.4361, 0.6687, 0.8598, 0.9685, 0.9972, 0.9999, 0.9999,  &
0.0806, 0.1732, 0.3312, 0.5457, 0.7757, 0.9260, 0.9884, 0.9992, 0.9999,  &
0.0468, 0.1117, 0.2527, 0.4542, 0.6687, 0.8697, 0.9685, 0.9963, 0.9998,  &
0.0228, 0.0619, 0.1616, 0.3481, 0.5638, 0.7892, 0.9322, 0.9884, 0.9991,  &
0.0090, 0.0283, 0.0806, 0.2107, 0.4361, 0.7013, 0.8882, 0.9716, 0.9968,  &
0.0031, 0.0102, 0.0348, 0.1032, 0.2676, 0.5457, 0.7892, 0.9380, 0.9884 /
data ((cldnuctab( 3,iw,iconc,6, 3),iw=1,9),iconc=1,7)/    &
0.1402, 0.2829, 0.5091, 0.7473, 0.9193, 0.9884, 0.9993, 0.9999, 0.9999,  &
0.0952, 0.2107, 0.3826, 0.6173, 0.8268, 0.9575, 0.9958, 0.9998, 0.9999,  &
0.0677, 0.1507, 0.2829, 0.4908, 0.7324, 0.9047, 0.9816, 0.9986, 0.9999,  &
0.0385, 0.0952, 0.2107, 0.4002, 0.6173, 0.8268, 0.9531, 0.9938, 0.9997,  &
0.0183, 0.0515, 0.1402, 0.2986, 0.5091, 0.7324, 0.9047, 0.9816, 0.9982,  &
0.0079, 0.0228, 0.0677, 0.1852, 0.3826, 0.6347, 0.8383, 0.9575, 0.9938,  &
0.0036, 0.0090, 0.0283, 0.0877, 0.2242, 0.4725, 0.7324, 0.9122, 0.9816 /
data ((cldnuctab( 3,iw,iconc,7, 3),iw=1,9),iconc=1,7)/    &
0.1302, 0.2676, 0.4725, 0.7170, 0.9047, 0.9854, 0.9991, 0.9999, 0.9999,  &
0.0877, 0.1852, 0.3481, 0.5818, 0.8022, 0.9485, 0.9938, 0.9997, 0.9999,  &
0.0565, 0.1302, 0.2676, 0.4542, 0.6852, 0.8792, 0.9771, 0.9982, 0.9999,  &
0.0348, 0.0877, 0.1852, 0.3481, 0.5818, 0.7892, 0.9380, 0.9909, 0.9994,  &
0.0163, 0.0468, 0.1207, 0.2527, 0.4542, 0.6852, 0.8792, 0.9716, 0.9972,  &
0.0070, 0.0205, 0.0619, 0.1616, 0.3312, 0.5818, 0.8022, 0.9380, 0.9898,  &
0.0036, 0.0079, 0.0254, 0.0739, 0.1977, 0.4181, 0.6852, 0.8792, 0.9716 /
data ((cldnuctab( 3,iw,iconc,1, 4),iw=1,9),iconc=1,7)/    &
0.4908, 0.7170, 0.8967, 0.9771, 0.9976, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.4002, 0.6347, 0.8383, 0.9531, 0.9920, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.2676, 0.5274, 0.7617, 0.9122, 0.9816, 0.9982, 0.9999, 0.9999, 0.9999,  &
0.1402, 0.3652, 0.6519, 0.8598, 0.9614, 0.9945, 0.9996, 0.9999, 0.9999,  &
0.0619, 0.1852, 0.4542, 0.7473, 0.9260, 0.9870, 0.9986, 0.9999, 0.9999,  &
0.0228, 0.0806, 0.2242, 0.4908, 0.7473, 0.9260, 0.9884, 0.9991, 0.9999,  &
0.0115, 0.0348, 0.0952, 0.1852, 0.2242, 0.2382, 0.3148, 0.4002, 0.8792 /
data ((cldnuctab( 3,iw,iconc,2, 4),iw=1,9),iconc=1,7)/    &
0.4002, 0.6347, 0.8383, 0.9531, 0.9938, 0.9996, 0.9999, 0.9999, 0.9999,  &
0.3148, 0.5274, 0.7617, 0.9122, 0.9816, 0.9982, 0.9999, 0.9999, 0.9999,  &
0.2107, 0.4181, 0.6687, 0.8598, 0.9614, 0.9945, 0.9996, 0.9999, 0.9999,  &
0.1117, 0.2829, 0.5457, 0.7892, 0.9322, 0.9870, 0.9986, 0.9999, 0.9999,  &
0.0468, 0.1507, 0.3652, 0.6687, 0.8792, 0.9716, 0.9963, 0.9997, 0.9999,  &
0.0183, 0.0619, 0.1852, 0.4542, 0.7473, 0.9260, 0.9870, 0.9988, 0.9999,  &
0.0090, 0.0254, 0.0806, 0.2107, 0.4002, 0.5997, 0.8383, 0.9575, 0.9979 /
data ((cldnuctab( 3,iw,iconc,3, 4),iw=1,9),iconc=1,7)/    &
0.3312, 0.5457, 0.7757, 0.9260, 0.9870, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.2527, 0.4542, 0.6687, 0.8697, 0.9651, 0.9958, 0.9997, 0.9999, 0.9999,  &
0.1732, 0.3481, 0.5818, 0.7892, 0.9322, 0.9870, 0.9988, 0.9999, 0.9999,  &
0.0877, 0.2242, 0.4542, 0.7013, 0.8882, 0.9716, 0.9963, 0.9997, 0.9999,  &
0.0348, 0.1117, 0.2986, 0.5818, 0.8148, 0.9434, 0.9909, 0.9991, 0.9999,  &
0.0129, 0.0468, 0.1507, 0.3826, 0.6852, 0.8967, 0.9771, 0.9972, 0.9998,  &
0.0090, 0.0183, 0.0619, 0.1852, 0.4181, 0.6852, 0.8967, 0.9836, 0.9984 /
data ((cldnuctab( 3,iw,iconc,4, 4),iw=1,9),iconc=1,7)/    &
0.2676, 0.4725, 0.7013, 0.8882, 0.9771, 0.9979, 0.9999, 0.9999, 0.9999,  &
0.1977, 0.3826, 0.5997, 0.8148, 0.9434, 0.9920, 0.9994, 0.9999, 0.9999,  &
0.1302, 0.2829, 0.4908, 0.7170, 0.8967, 0.9771, 0.9976, 0.9998, 0.9999,  &
0.0739, 0.1852, 0.3826, 0.6173, 0.8268, 0.9485, 0.9920, 0.9993, 0.9999,  &
0.0283, 0.0877, 0.2382, 0.4908, 0.7473, 0.9122, 0.9795, 0.9979, 0.9998,  &
0.0115, 0.0348, 0.1117, 0.2986, 0.5997, 0.8383, 0.9575, 0.9938, 0.9995,  &
0.0090, 0.0145, 0.0468, 0.1402, 0.3652, 0.6687, 0.8882, 0.9771, 0.9976 /
data ((cldnuctab( 3,iw,iconc,5, 4),iw=1,9),iconc=1,7)/    &
0.2382, 0.4181, 0.6519, 0.8598, 0.9685, 0.9963, 0.9998, 0.9999, 0.9999,  &
0.1616, 0.3312, 0.5457, 0.7617, 0.9193, 0.9870, 0.9990, 0.9999, 0.9999,  &
0.1117, 0.2382, 0.4361, 0.6687, 0.8598, 0.9614, 0.9952, 0.9997, 0.9999,  &
0.0619, 0.1507, 0.3312, 0.5457, 0.7757, 0.9260, 0.9854, 0.9986, 0.9999,  &
0.0254, 0.0739, 0.1977, 0.4181, 0.6687, 0.8697, 0.9651, 0.9952, 0.9997,  &
0.0102, 0.0283, 0.0877, 0.2527, 0.5274, 0.7892, 0.9322, 0.9870, 0.9988,  &
0.0102, 0.0115, 0.0348, 0.1117, 0.2986, 0.5997, 0.8493, 0.9651, 0.9952 /
data ((cldnuctab( 3,iw,iconc,6, 4),iw=1,9),iconc=1,7)/    &
0.2107, 0.3826, 0.6173, 0.8268, 0.9575, 0.9952, 0.9997, 0.9999, 0.9999,  &
0.1402, 0.2829, 0.4908, 0.7170, 0.8967, 0.9816, 0.9984, 0.9999, 0.9999,  &
0.0952, 0.2107, 0.3826, 0.6173, 0.8148, 0.9485, 0.9929, 0.9995, 0.9999,  &
0.0515, 0.1302, 0.2829, 0.4908, 0.7170, 0.8967, 0.9771, 0.9976, 0.9998,  &
0.0228, 0.0619, 0.1732, 0.3652, 0.6173, 0.8268, 0.9485, 0.9920, 0.9993,  &
0.0102, 0.0254, 0.0739, 0.2107, 0.4542, 0.7170, 0.9047, 0.9795, 0.9976,  &
0.0102, 0.0115, 0.0314, 0.0952, 0.2527, 0.5274, 0.8022, 0.9485, 0.9920 /
data ((cldnuctab( 3,iw,iconc,7, 4),iw=1,9),iconc=1,7)/    &
0.1852, 0.3652, 0.5818, 0.8148, 0.9485, 0.9938, 0.9996, 0.9999, 0.9999,  &
0.1302, 0.2676, 0.4542, 0.6852, 0.8792, 0.9745, 0.9979, 0.9999, 0.9999,  &
0.0877, 0.1852, 0.3481, 0.5638, 0.7892, 0.9322, 0.9898, 0.9992, 0.9999,  &
0.0468, 0.1207, 0.2527, 0.4542, 0.6852, 0.8697, 0.9685, 0.9963, 0.9998,  &
0.0205, 0.0565, 0.1507, 0.3312, 0.5638, 0.7892, 0.9322, 0.9870, 0.9990,  &
0.0115, 0.0228, 0.0677, 0.1852, 0.4002, 0.6687, 0.8697, 0.9685, 0.9958,  &
0.0115, 0.0129, 0.0254, 0.0806, 0.2107, 0.4725, 0.7617, 0.9322, 0.9870 /
data ((cldnuctab( 3,iw,iconc,1, 5),iw=1,9),iconc=1,7)/    &
0.5638, 0.7757, 0.9260, 0.9854, 0.9986, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4542, 0.6852, 0.8792, 0.9685, 0.9958, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.2986, 0.5818, 0.8022, 0.9380, 0.9884, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.1507, 0.4002, 0.7013, 0.8967, 0.9745, 0.9968, 0.9997, 0.9999, 0.9999,  &
0.0619, 0.1977, 0.4725, 0.7892, 0.9434, 0.9909, 0.9992, 0.9999, 0.9999,  &
0.0254, 0.0877, 0.2382, 0.4908, 0.7617, 0.9380, 0.9920, 0.9995, 0.9999,  &
0.0183, 0.0385, 0.0952, 0.1852, 0.2107, 0.2242, 0.2527, 0.3312, 0.8792 /
data ((cldnuctab( 3,iw,iconc,2, 5),iw=1,9),iconc=1,7)/    &
0.4725, 0.7013, 0.8792, 0.9716, 0.9963, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.3652, 0.5997, 0.8148, 0.9380, 0.9898, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.2382, 0.4725, 0.7170, 0.8967, 0.9745, 0.9968, 0.9998, 0.9999, 0.9999,  &
0.1207, 0.3148, 0.5997, 0.8268, 0.9531, 0.9920, 0.9992, 0.9999, 0.9999,  &
0.0468, 0.1507, 0.4002, 0.7170, 0.9122, 0.9816, 0.9979, 0.9998, 0.9999,  &
0.0183, 0.0619, 0.1977, 0.4725, 0.7757, 0.9485, 0.9920, 0.9993, 0.9999,  &
0.0163, 0.0283, 0.0806, 0.2107, 0.4002, 0.5818, 0.8022, 0.9651, 0.9986 /
data ((cldnuctab( 3,iw,iconc,3, 5),iw=1,9),iconc=1,7)/    &
0.3826, 0.6173, 0.8268, 0.9485, 0.9920, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.2986, 0.5091, 0.7324, 0.9047, 0.9795, 0.9976, 0.9998, 0.9999, 0.9999,  &
0.1977, 0.4002, 0.6347, 0.8383, 0.9531, 0.9929, 0.9994, 0.9999, 0.9999,  &
0.0952, 0.2527, 0.5091, 0.7617, 0.9193, 0.9816, 0.9979, 0.9998, 0.9999,  &
0.0385, 0.1207, 0.3148, 0.6173, 0.8598, 0.9614, 0.9945, 0.9995, 0.9999,  &
0.0163, 0.0468, 0.1507, 0.4002, 0.7170, 0.9193, 0.9854, 0.9984, 0.9999,  &
0.0163, 0.0205, 0.0619, 0.1852, 0.4181, 0.7013, 0.9122, 0.9854, 0.9988 /
data ((cldnuctab( 3,iw,iconc,4, 5),iw=1,9),iconc=1,7)/    &
0.3312, 0.5457, 0.7617, 0.9260, 0.9870, 0.9988, 0.9999, 0.9999, 0.9999,  &
0.2382, 0.4361, 0.6687, 0.8598, 0.9651, 0.9952, 0.9997, 0.9999, 0.9999,  &
0.1616, 0.3312, 0.5638, 0.7757, 0.9260, 0.9854, 0.9986, 0.9999, 0.9999,  &
0.0806, 0.2107, 0.4361, 0.6852, 0.8697, 0.9685, 0.9952, 0.9996, 0.9999,  &
0.0314, 0.0952, 0.2676, 0.5457, 0.7892, 0.9380, 0.9884, 0.9988, 0.9999,  &
0.0163, 0.0385, 0.1207, 0.3312, 0.6347, 0.8792, 0.9716, 0.9963, 0.9997,  &
0.0163, 0.0205, 0.0515, 0.1507, 0.3826, 0.6852, 0.9047, 0.9854, 0.9984 /
data ((cldnuctab( 3,iw,iconc,5, 5),iw=1,9),iconc=1,7)/    &
0.2829, 0.4908, 0.7170, 0.8967, 0.9795, 0.9982, 0.9999, 0.9999, 0.9999,  &
0.1977, 0.3826, 0.5997, 0.8148, 0.9485, 0.9920, 0.9994, 0.9999, 0.9999,  &
0.1302, 0.2829, 0.4908, 0.7170, 0.8967, 0.9771, 0.9972, 0.9998, 0.9999,  &
0.0677, 0.1732, 0.3652, 0.6173, 0.8268, 0.9485, 0.9909, 0.9992, 0.9999,  &
0.0254, 0.0806, 0.2107, 0.4725, 0.7324, 0.9047, 0.9771, 0.9972, 0.9998,  &
0.0183, 0.0314, 0.0952, 0.2676, 0.5638, 0.8268, 0.9531, 0.9920, 0.9993,  &
0.0183, 0.0205, 0.0385, 0.1207, 0.3148, 0.6347, 0.8792, 0.9745, 0.9972 /
data ((cldnuctab( 3,iw,iconc,6, 5),iw=1,9),iconc=1,7)/    &
0.2527, 0.4542, 0.6852, 0.8792, 0.9745, 0.9972, 0.9998, 0.9999, 0.9999,  &
0.1732, 0.3481, 0.5638, 0.7757, 0.9322, 0.9884, 0.9991, 0.9999, 0.9999,  &
0.1117, 0.2527, 0.4361, 0.6687, 0.8598, 0.9651, 0.9958, 0.9997, 0.9999,  &
0.0565, 0.1507, 0.3312, 0.5638, 0.7757, 0.9260, 0.9854, 0.9986, 0.9999,  &
0.0228, 0.0677, 0.1852, 0.4002, 0.6687, 0.8697, 0.9651, 0.9952, 0.9996,  &
0.0183, 0.0254, 0.0806, 0.2242, 0.4908, 0.7757, 0.9322, 0.9870, 0.9986,  &
0.0183, 0.0205, 0.0314, 0.0952, 0.2676, 0.5638, 0.8383, 0.9651, 0.9952 /
data ((cldnuctab( 3,iw,iconc,7, 5),iw=1,9),iconc=1,7)/    &
0.2242, 0.4181, 0.6519, 0.8598, 0.9685, 0.9968, 0.9998, 0.9999, 0.9999,  &
0.1616, 0.3148, 0.5274, 0.7473, 0.9193, 0.9854, 0.9988, 0.9999, 0.9999,  &
0.1032, 0.2242, 0.4002, 0.6347, 0.8383, 0.9575, 0.9945, 0.9996, 0.9999,  &
0.0515, 0.1402, 0.2986, 0.5091, 0.7473, 0.9047, 0.9816, 0.9982, 0.9999,  &
0.0205, 0.0619, 0.1616, 0.3652, 0.6173, 0.8383, 0.9531, 0.9929, 0.9994,  &
0.0205, 0.0228, 0.0677, 0.1977, 0.4361, 0.7170, 0.9047, 0.9795, 0.9976,  &
0.0205, 0.0228, 0.0283, 0.0806, 0.2242, 0.5091, 0.8022, 0.9485, 0.9920 /
data ((cldnuctab( 3,iw,iconc,1, 6),iw=1,9),iconc=1,7)/    &
0.5997, 0.8148, 0.9434, 0.9898, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4908, 0.7324, 0.8967, 0.9771, 0.9972, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.3312, 0.6173, 0.8383, 0.9531, 0.9920, 0.9993, 0.9999, 0.9999, 0.9999,  &
0.1616, 0.4181, 0.7324, 0.9122, 0.9816, 0.9979, 0.9998, 0.9999, 0.9999,  &
0.0619, 0.1977, 0.4908, 0.8022, 0.9575, 0.9938, 0.9995, 0.9999, 0.9999,  &
0.0314, 0.0877, 0.2382, 0.4908, 0.7617, 0.9434, 0.9938, 0.9996, 0.9999,  &
0.0254, 0.0468, 0.1032, 0.1852, 0.2107, 0.2242, 0.2527, 0.3481, 0.8148 /
data ((cldnuctab( 3,iw,iconc,2, 6),iw=1,9),iconc=1,7)/    &
0.5091, 0.7324, 0.9047, 0.9795, 0.9976, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.4002, 0.6347, 0.8383, 0.9531, 0.9929, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.2527, 0.5091, 0.7617, 0.9193, 0.9816, 0.9979, 0.9998, 0.9999, 0.9999,  &
0.1207, 0.3312, 0.6347, 0.8598, 0.9651, 0.9945, 0.9995, 0.9999, 0.9999,  &
0.0515, 0.1616, 0.4181, 0.7473, 0.9260, 0.9870, 0.9986, 0.9999, 0.9999,  &
0.0254, 0.0677, 0.1977, 0.4725, 0.7892, 0.9575, 0.9945, 0.9995, 0.9999,  &
0.0254, 0.0348, 0.0877, 0.2107, 0.3826, 0.5818, 0.7892, 0.9685, 0.9986 /
data ((cldnuctab( 3,iw,iconc,3, 6),iw=1,9),iconc=1,7)/    &
0.4361, 0.6687, 0.8598, 0.9614, 0.9945, 0.9996, 0.9999, 0.9999, 0.9999,  &
0.3312, 0.5638, 0.7757, 0.9260, 0.9854, 0.9984, 0.9999, 0.9999, 0.9999,  &
0.2107, 0.4361, 0.6852, 0.8697, 0.9651, 0.9952, 0.9996, 0.9999, 0.9999,  &
0.1032, 0.2676, 0.5457, 0.7892, 0.9322, 0.9870, 0.9986, 0.9999, 0.9999,  &
0.0385, 0.1207, 0.3312, 0.6519, 0.8792, 0.9716, 0.9963, 0.9997, 0.9999,  &
0.0254, 0.0515, 0.1507, 0.4002, 0.7473, 0.9322, 0.9898, 0.9990, 0.9999,  &
0.0228, 0.0283, 0.0677, 0.1852, 0.4181, 0.7013, 0.9193, 0.9884, 0.9992 /
data ((cldnuctab( 3,iw,iconc,4, 6),iw=1,9),iconc=1,7)/    &
0.3652, 0.5818, 0.8022, 0.9434, 0.9909, 0.9993, 0.9999, 0.9999, 0.9999,  &
0.2676, 0.4725, 0.7170, 0.8882, 0.9745, 0.9968, 0.9998, 0.9999, 0.9999,  &
0.1732, 0.3652, 0.5997, 0.8148, 0.9434, 0.9898, 0.9991, 0.9999, 0.9999,  &
0.0806, 0.2242, 0.4725, 0.7170, 0.8967, 0.9771, 0.9968, 0.9997, 0.9999,  &
0.0314, 0.1032, 0.2676, 0.5638, 0.8268, 0.9485, 0.9920, 0.9992, 0.9999,  &
0.0254, 0.0385, 0.1207, 0.3312, 0.6687, 0.8967, 0.9795, 0.9976, 0.9998,  &
0.0254, 0.0283, 0.0515, 0.1507, 0.3826, 0.6852, 0.9122, 0.9884, 0.9988 /
data ((cldnuctab( 3,iw,iconc,5, 6),iw=1,9),iconc=1,7)/    &
0.3148, 0.5274, 0.7617, 0.9193, 0.9854, 0.9988, 0.9999, 0.9999, 0.9999,  &
0.2382, 0.4181, 0.6519, 0.8493, 0.9614, 0.9945, 0.9996, 0.9999, 0.9999,  &
0.1507, 0.3148, 0.5457, 0.7617, 0.9193, 0.9836, 0.9984, 0.9999, 0.9999,  &
0.0677, 0.1852, 0.4002, 0.6519, 0.8598, 0.9614, 0.9945, 0.9995, 0.9999,  &
0.0283, 0.0806, 0.2242, 0.4908, 0.7617, 0.9260, 0.9836, 0.9984, 0.9999,  &
0.0254, 0.0314, 0.0952, 0.2676, 0.5818, 0.8493, 0.9651, 0.9945, 0.9996,  &
0.0254, 0.0283, 0.0425, 0.1207, 0.3148, 0.6347, 0.8882, 0.9816, 0.9982 /
data ((cldnuctab( 3,iw,iconc,6, 6),iw=1,9),iconc=1,7)/    &
0.2829, 0.4908, 0.7170, 0.9047, 0.9816, 0.9984, 0.9999, 0.9999, 0.9999,  &
0.2107, 0.3826, 0.5997, 0.8148, 0.9485, 0.9920, 0.9994, 0.9999, 0.9999,  &
0.1302, 0.2829, 0.4908, 0.7170, 0.8882, 0.9771, 0.9972, 0.9998, 0.9999,  &
0.0619, 0.1616, 0.3481, 0.5997, 0.8148, 0.9434, 0.9909, 0.9992, 0.9999,  &
0.0283, 0.0739, 0.1977, 0.4361, 0.7013, 0.8882, 0.9745, 0.9968, 0.9997,  &
0.0283, 0.0314, 0.0806, 0.2382, 0.5091, 0.8022, 0.9434, 0.9909, 0.9991,  &
0.0283, 0.0314, 0.0348, 0.1032, 0.2676, 0.5818, 0.8598, 0.9716, 0.9968 /
data ((cldnuctab( 3,iw,iconc,7, 6),iw=1,9),iconc=1,7)/    &
0.2676, 0.4542, 0.6852, 0.8882, 0.9771, 0.9979, 0.9999, 0.9999, 0.9999,  &
0.1852, 0.3481, 0.5638, 0.7892, 0.9380, 0.9898, 0.9992, 0.9999, 0.9999,  &
0.1207, 0.2527, 0.4542, 0.6852, 0.8697, 0.9685, 0.9963, 0.9997, 0.9999,  &
0.0565, 0.1507, 0.3148, 0.5638, 0.7757, 0.9260, 0.9870, 0.9988, 0.9999,  &
0.0283, 0.0619, 0.1732, 0.4002, 0.6519, 0.8598, 0.9651, 0.9952, 0.9996,  &
0.0283, 0.0314, 0.0739, 0.1977, 0.4542, 0.7473, 0.9260, 0.9854, 0.9984,  &
0.0283, 0.0314, 0.0385, 0.0877, 0.2382, 0.5274, 0.8268, 0.9614, 0.9945 /
data ((cldnuctab( 3,iw,iconc,1, 7),iw=1,9),iconc=1,7)/    &
0.6347, 0.8383, 0.9531, 0.9929, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5274, 0.7617, 0.9193, 0.9816, 0.9979, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.3312, 0.6347, 0.8598, 0.9651, 0.9945, 0.9995, 0.9999, 0.9999, 0.9999,  &
0.1616, 0.4181, 0.7473, 0.9260, 0.9870, 0.9986, 0.9999, 0.9999, 0.9999,  &
0.0677, 0.2107, 0.5091, 0.8268, 0.9651, 0.9952, 0.9996, 0.9999, 0.9999,  &
0.0348, 0.0952, 0.2382, 0.4908, 0.7617, 0.9485, 0.9952, 0.9997, 0.9999,  &
0.0314, 0.0515, 0.1032, 0.1732, 0.2107, 0.2242, 0.2676, 0.3652, 0.8148 /
data ((cldnuctab( 3,iw,iconc,2, 7),iw=1,9),iconc=1,7)/    &
0.5457, 0.7757, 0.9193, 0.9836, 0.9984, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4361, 0.6687, 0.8697, 0.9651, 0.9945, 0.9996, 0.9999, 0.9999, 0.9999,  &
0.2676, 0.5457, 0.7892, 0.9322, 0.9870, 0.9986, 0.9999, 0.9999, 0.9999,  &
0.1302, 0.3481, 0.6519, 0.8792, 0.9716, 0.9958, 0.9997, 0.9999, 0.9999,  &
0.0515, 0.1616, 0.4181, 0.7617, 0.9380, 0.9898, 0.9990, 0.9999, 0.9999,  &
0.0348, 0.0677, 0.1977, 0.4908, 0.8022, 0.9614, 0.9958, 0.9996, 0.9999,  &
0.0314, 0.0425, 0.0952, 0.2107, 0.3826, 0.5638, 0.7892, 0.9685, 0.9988 /
data ((cldnuctab( 3,iw,iconc,3, 7),iw=1,9),iconc=1,7)/    &
0.4725, 0.7013, 0.8792, 0.9716, 0.9963, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.3481, 0.5818, 0.8022, 0.9380, 0.9884, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.2242, 0.4542, 0.7170, 0.8882, 0.9745, 0.9963, 0.9997, 0.9999, 0.9999,  &
0.1032, 0.2829, 0.5638, 0.8148, 0.9485, 0.9909, 0.9991, 0.9999, 0.9999,  &
0.0385, 0.1302, 0.3481, 0.6687, 0.8967, 0.9795, 0.9972, 0.9998, 0.9999,  &
0.0314, 0.0515, 0.1616, 0.4181, 0.7617, 0.9434, 0.9920, 0.9992, 0.9999,  &
0.0314, 0.0385, 0.0739, 0.1977, 0.4181, 0.7013, 0.9260, 0.9909, 0.9994 /
data ((cldnuctab( 3,iw,iconc,4, 7),iw=1,9),iconc=1,7)/    &
0.4002, 0.6347, 0.8383, 0.9531, 0.9929, 0.9995, 0.9999, 0.9999, 0.9999,  &
0.2986, 0.5091, 0.7473, 0.9047, 0.9795, 0.9979, 0.9998, 0.9999, 0.9999,  &
0.1852, 0.3826, 0.6347, 0.8383, 0.9531, 0.9929, 0.9994, 0.9999, 0.9999,  &
0.0877, 0.2242, 0.4908, 0.7473, 0.9122, 0.9816, 0.9979, 0.9998, 0.9999,  &
0.0314, 0.1032, 0.2829, 0.5818, 0.8493, 0.9614, 0.9938, 0.9994, 0.9999,  &
0.0314, 0.0425, 0.1207, 0.3481, 0.6852, 0.9122, 0.9836, 0.9982, 0.9998,  &
0.0314, 0.0385, 0.0565, 0.1616, 0.3826, 0.7013, 0.9193, 0.9909, 0.9992 /
data ((cldnuctab( 3,iw,iconc,5, 7),iw=1,9),iconc=1,7)/    &
0.3481, 0.5638, 0.7892, 0.9380, 0.9898, 0.9992, 0.9999, 0.9999, 0.9999,  &
0.2527, 0.4542, 0.6852, 0.8697, 0.9685, 0.9963, 0.9997, 0.9999, 0.9999,  &
0.1616, 0.3312, 0.5638, 0.7892, 0.9322, 0.9870, 0.9988, 0.9999, 0.9999,  &
0.0739, 0.1977, 0.4181, 0.6852, 0.8792, 0.9685, 0.9958, 0.9997, 0.9999,  &
0.0348, 0.0877, 0.2382, 0.5091, 0.7892, 0.9380, 0.9884, 0.9988, 0.9999,  &
0.0348, 0.0385, 0.1032, 0.2829, 0.5997, 0.8697, 0.9716, 0.9963, 0.9997,  &
0.0348, 0.0385, 0.0468, 0.1302, 0.3148, 0.6519, 0.8967, 0.9854, 0.9986 /
data ((cldnuctab( 3,iw,iconc,6, 7),iw=1,9),iconc=1,7)/    &
0.3148, 0.5274, 0.7473, 0.9193, 0.9854, 0.9988, 0.9999, 0.9999, 0.9999,  &
0.2242, 0.4181, 0.6347, 0.8383, 0.9575, 0.9945, 0.9996, 0.9999, 0.9999,  &
0.1402, 0.2986, 0.5274, 0.7473, 0.9122, 0.9816, 0.9982, 0.9999, 0.9999,  &
0.0619, 0.1732, 0.3826, 0.6347, 0.8383, 0.9531, 0.9929, 0.9994, 0.9999,  &
0.0348, 0.0739, 0.1977, 0.4542, 0.7324, 0.9122, 0.9816, 0.9979, 0.9998,  &
0.0348, 0.0385, 0.0877, 0.2382, 0.5274, 0.8148, 0.9531, 0.9929, 0.9994,  &
0.0348, 0.0385, 0.0468, 0.1032, 0.2676, 0.5818, 0.8697, 0.9771, 0.9976 /
data ((cldnuctab( 3,iw,iconc,7, 7),iw=1,9),iconc=1,7)/    &
0.2986, 0.4908, 0.7324, 0.9047, 0.9816, 0.9984, 0.9999, 0.9999, 0.9999,  &
0.2107, 0.3826, 0.5997, 0.8148, 0.9485, 0.9929, 0.9995, 0.9999, 0.9999,  &
0.1302, 0.2676, 0.4908, 0.7170, 0.8882, 0.9745, 0.9972, 0.9998, 0.9999,  &
0.0565, 0.1507, 0.3481, 0.5818, 0.8022, 0.9434, 0.9898, 0.9991, 0.9999,  &
0.0385, 0.0677, 0.1852, 0.4181, 0.6852, 0.8792, 0.9716, 0.9963, 0.9997,  &
0.0385, 0.0425, 0.0739, 0.2107, 0.4725, 0.7757, 0.9380, 0.9884, 0.9990,  &
0.0385, 0.0425, 0.0468, 0.0952, 0.2382, 0.5274, 0.8383, 0.9685, 0.9963 /
data ((cldnuctab( 4,iw,iconc,1, 1),iw=1,9),iconc=1,7)/    &
0.5136, 0.7509, 0.9140, 0.9858, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4046, 0.6560, 0.8519, 0.9660, 0.9965, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.2867, 0.5319, 0.7791, 0.9337, 0.9901, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.1427, 0.3694, 0.6728, 0.8815, 0.9751, 0.9977, 0.9999, 0.9999, 0.9999,  &
0.0578, 0.1882, 0.4587, 0.7652, 0.9447, 0.9922, 0.9996, 0.9999, 0.9999,  &
0.0235, 0.0823, 0.2276, 0.4953, 0.7509, 0.9394, 0.9932, 0.9997, 0.9999,  &
0.0118, 0.0357, 0.0971, 0.1882, 0.2276, 0.2276, 0.2417, 0.4046, 0.8987 /
data ((cldnuctab( 4,iw,iconc,2, 1),iw=1,9),iconc=1,7)/    &
0.4046, 0.6560, 0.8519, 0.9660, 0.9969, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3188, 0.5501, 0.7791, 0.9276, 0.9901, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.2140, 0.4405, 0.6892, 0.8815, 0.9724, 0.9977, 0.9999, 0.9999, 0.9999,  &
0.1139, 0.2867, 0.5683, 0.8053, 0.9447, 0.9922, 0.9996, 0.9999, 0.9999,  &
0.0479, 0.1427, 0.3869, 0.6892, 0.8987, 0.9800, 0.9985, 0.9999, 0.9999,  &
0.0168, 0.0633, 0.1882, 0.4587, 0.7652, 0.9447, 0.9932, 0.9996, 0.9999,  &
0.0093, 0.0261, 0.0823, 0.2140, 0.4046, 0.6040, 0.8519, 0.9693, 0.9991 /
data ((cldnuctab( 4,iw,iconc,3, 1),iw=1,9),iconc=1,7)/    &
0.3353, 0.5683, 0.7925, 0.9394, 0.9922, 0.9996, 0.9999, 0.9999, 0.9999,  &
0.2563, 0.4587, 0.6892, 0.8815, 0.9751, 0.9982, 0.9999, 0.9999, 0.9999,  &
0.1761, 0.3522, 0.5862, 0.8053, 0.9447, 0.9932, 0.9996, 0.9999, 0.9999,  &
0.0895, 0.2276, 0.4770, 0.7208, 0.8987, 0.9800, 0.9985, 0.9999, 0.9999,  &
0.0357, 0.1139, 0.3025, 0.5862, 0.8296, 0.9585, 0.9947, 0.9997, 0.9999,  &
0.0133, 0.0479, 0.1427, 0.3869, 0.7052, 0.9140, 0.9858, 0.9989, 0.9999,  &
0.0093, 0.0188, 0.0578, 0.1761, 0.4225, 0.7052, 0.9140, 0.9901, 0.9994 /
data ((cldnuctab( 4,iw,iconc,4, 1),iw=1,9),iconc=1,7)/    &
0.2713, 0.4953, 0.7208, 0.9066, 0.9858, 0.9991, 0.9999, 0.9999, 0.9999,  &
0.2008, 0.3869, 0.6216, 0.8296, 0.9585, 0.9954, 0.9998, 0.9999, 0.9999,  &
0.1326, 0.2867, 0.5136, 0.7361, 0.9140, 0.9841, 0.9990, 0.9999, 0.9999,  &
0.0692, 0.1882, 0.3869, 0.6390, 0.8519, 0.9624, 0.9959, 0.9998, 0.9999,  &
0.0290, 0.0895, 0.2417, 0.4953, 0.7652, 0.9276, 0.9874, 0.9991, 0.9999,  &
0.0105, 0.0357, 0.1139, 0.3025, 0.6216, 0.8623, 0.9693, 0.9969, 0.9998,  &
0.0093, 0.0133, 0.0479, 0.1427, 0.3694, 0.6728, 0.9066, 0.9858, 0.9989 /
data ((cldnuctab( 4,iw,iconc,5, 1),iw=1,9),iconc=1,7)/    &
0.2417, 0.4225, 0.6728, 0.8721, 0.9751, 0.9985, 0.9999, 0.9999, 0.9999,  &
0.1644, 0.3353, 0.5501, 0.7791, 0.9337, 0.9912, 0.9996, 0.9999, 0.9999,  &
0.1139, 0.2417, 0.4405, 0.6728, 0.8721, 0.9724, 0.9977, 0.9999, 0.9999,  &
0.0633, 0.1533, 0.3353, 0.5683, 0.7925, 0.9394, 0.9912, 0.9994, 0.9999,  &
0.0261, 0.0755, 0.2008, 0.4225, 0.6892, 0.8815, 0.9751, 0.9977, 0.9999,  &
0.0105, 0.0290, 0.0895, 0.2563, 0.5319, 0.7925, 0.9447, 0.9922, 0.9995,  &
0.0105, 0.0118, 0.0357, 0.1139, 0.3025, 0.6040, 0.8721, 0.9751, 0.9977 /
data ((cldnuctab( 4,iw,iconc,6, 1),iw=1,9),iconc=1,7)/    &
0.2140, 0.3869, 0.6216, 0.8410, 0.9660, 0.9977, 0.9999, 0.9999, 0.9999,  &
0.1427, 0.2867, 0.4953, 0.7361, 0.9140, 0.9874, 0.9994, 0.9999, 0.9999,  &
0.0971, 0.2140, 0.3869, 0.6216, 0.8296, 0.9585, 0.9959, 0.9998, 0.9999,  &
0.0527, 0.1326, 0.2867, 0.5136, 0.7361, 0.9066, 0.9841, 0.9990, 0.9999,  &
0.0210, 0.0633, 0.1761, 0.3694, 0.6216, 0.8410, 0.9585, 0.9954, 0.9998,  &
0.0105, 0.0261, 0.0755, 0.2140, 0.4587, 0.7361, 0.9140, 0.9858, 0.9989,  &
0.0105, 0.0118, 0.0290, 0.0895, 0.2563, 0.5501, 0.8177, 0.9585, 0.9954 /
data ((cldnuctab( 4,iw,iconc,7, 1),iw=1,9),iconc=1,7)/    &
0.1882, 0.3694, 0.5862, 0.8177, 0.9585, 0.9965, 0.9999, 0.9999, 0.9999,  &
0.1326, 0.2713, 0.4587, 0.7052, 0.8987, 0.9821, 0.9990, 0.9999, 0.9999,  &
0.0895, 0.1882, 0.3522, 0.5862, 0.8053, 0.9447, 0.9940, 0.9997, 0.9999,  &
0.0479, 0.1230, 0.2563, 0.4587, 0.6892, 0.8815, 0.9777, 0.9982, 0.9999,  &
0.0188, 0.0578, 0.1533, 0.3353, 0.5683, 0.8053, 0.9447, 0.9922, 0.9996,  &
0.0118, 0.0235, 0.0692, 0.1882, 0.4046, 0.6892, 0.8815, 0.9751, 0.9980,  &
0.0118, 0.0133, 0.0261, 0.0823, 0.2140, 0.4953, 0.7791, 0.9394, 0.9922 /
data ((cldnuctab( 4,iw,iconc,1, 2),iw=1,9),iconc=1,7)/    &
0.6216, 0.8296, 0.9496, 0.9932, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4953, 0.7361, 0.9066, 0.9821, 0.9982, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3188, 0.6216, 0.8519, 0.9624, 0.9947, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.1533, 0.4046, 0.7361, 0.9210, 0.9858, 0.9989, 0.9999, 0.9999, 0.9999,  &
0.0633, 0.2008, 0.4953, 0.8177, 0.9624, 0.9959, 0.9998, 0.9999, 0.9999,  &
0.0290, 0.0895, 0.2417, 0.4953, 0.7652, 0.9496, 0.9959, 0.9998, 0.9999,  &
0.0261, 0.0435, 0.0971, 0.1882, 0.2140, 0.2276, 0.2563, 0.3522, 0.8296 /
data ((cldnuctab( 4,iw,iconc,2, 2),iw=1,9),iconc=1,7)/    &
0.5136, 0.7509, 0.9140, 0.9841, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4046, 0.6560, 0.8519, 0.9624, 0.9947, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.2563, 0.5136, 0.7652, 0.9276, 0.9858, 0.9989, 0.9999, 0.9999, 0.9999,  &
0.1230, 0.3353, 0.6390, 0.8721, 0.9693, 0.9965, 0.9998, 0.9999, 0.9999,  &
0.0479, 0.1533, 0.4225, 0.7509, 0.9337, 0.9901, 0.9993, 0.9999, 0.9999,  &
0.0261, 0.0692, 0.2008, 0.4770, 0.8053, 0.9585, 0.9959, 0.9998, 0.9999,  &
0.0235, 0.0357, 0.0895, 0.2140, 0.3869, 0.5683, 0.7925, 0.9724, 0.9993 /
data ((cldnuctab( 4,iw,iconc,3, 2),iw=1,9),iconc=1,7)/    &
0.4405, 0.6728, 0.8623, 0.9660, 0.9965, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.3353, 0.5683, 0.7791, 0.9337, 0.9888, 0.9991, 0.9999, 0.9999, 0.9999,  &
0.2140, 0.4405, 0.6892, 0.8815, 0.9724, 0.9969, 0.9998, 0.9999, 0.9999,  &
0.0971, 0.2713, 0.5501, 0.8053, 0.9394, 0.9901, 0.9993, 0.9999, 0.9999,  &
0.0394, 0.1230, 0.3353, 0.6560, 0.8903, 0.9777, 0.9977, 0.9998, 0.9999,  &
0.0235, 0.0527, 0.1533, 0.4046, 0.7509, 0.9394, 0.9922, 0.9994, 0.9999,  &
0.0235, 0.0290, 0.0692, 0.1882, 0.4225, 0.7052, 0.9276, 0.9912, 0.9996 /
data ((cldnuctab( 4,iw,iconc,4, 2),iw=1,9),iconc=1,7)/    &
0.3694, 0.6040, 0.8177, 0.9496, 0.9932, 0.9996, 0.9999, 0.9999, 0.9999,  &
0.2713, 0.4770, 0.7208, 0.8903, 0.9777, 0.9980, 0.9999, 0.9999, 0.9999,  &
0.1761, 0.3694, 0.6040, 0.8177, 0.9496, 0.9922, 0.9995, 0.9999, 0.9999,  &
0.0823, 0.2276, 0.4770, 0.7208, 0.9066, 0.9800, 0.9980, 0.9999, 0.9999,  &
0.0322, 0.0971, 0.2713, 0.5683, 0.8296, 0.9542, 0.9940, 0.9996, 0.9999,  &
0.0261, 0.0394, 0.1230, 0.3353, 0.6728, 0.9066, 0.9841, 0.9985, 0.9999,  &
0.0261, 0.0290, 0.0527, 0.1533, 0.3869, 0.6892, 0.9210, 0.9912, 0.9994 /
data ((cldnuctab( 4,iw,iconc,5, 2),iw=1,9),iconc=1,7)/    &
0.3188, 0.5319, 0.7652, 0.9276, 0.9888, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.2276, 0.4225, 0.6560, 0.8519, 0.9660, 0.9965, 0.9998, 0.9999, 0.9999,  &
0.1533, 0.3188, 0.5501, 0.7652, 0.9210, 0.9858, 0.9990, 0.9999, 0.9999,  &
0.0692, 0.1882, 0.4046, 0.6560, 0.8623, 0.9660, 0.9959, 0.9998, 0.9999,  &
0.0261, 0.0823, 0.2276, 0.4953, 0.7652, 0.9276, 0.9874, 0.9990, 0.9999,  &
0.0261, 0.0322, 0.0971, 0.2713, 0.5862, 0.8519, 0.9693, 0.9965, 0.9998,  &
0.0261, 0.0290, 0.0435, 0.1230, 0.3188, 0.6390, 0.8987, 0.9841, 0.9989 /
data ((cldnuctab( 4,iw,iconc,6, 2),iw=1,9),iconc=1,7)/    &
0.2867, 0.4953, 0.7208, 0.9066, 0.9841, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.2008, 0.3869, 0.6040, 0.8177, 0.9542, 0.9940, 0.9997, 0.9999, 0.9999,  &
0.1326, 0.2713, 0.4953, 0.7208, 0.8987, 0.9800, 0.9982, 0.9999, 0.9999,  &
0.0633, 0.1644, 0.3522, 0.6040, 0.8177, 0.9496, 0.9922, 0.9995, 0.9999,  &
0.0290, 0.0692, 0.2008, 0.4405, 0.7052, 0.8987, 0.9777, 0.9980, 0.9999,  &
0.0290, 0.0322, 0.0823, 0.2276, 0.5136, 0.8053, 0.9496, 0.9932, 0.9995,  &
0.0261, 0.0290, 0.0357, 0.0971, 0.2713, 0.5862, 0.8623, 0.9751, 0.9977 /
data ((cldnuctab( 4,iw,iconc,7, 2),iw=1,9),iconc=1,7)/    &
0.2713, 0.4587, 0.7052, 0.8903, 0.9800, 0.9985, 0.9999, 0.9999, 0.9999,  &
0.1882, 0.3522, 0.5683, 0.7925, 0.9394, 0.9922, 0.9996, 0.9999, 0.9999,  &
0.1139, 0.2563, 0.4587, 0.6892, 0.8721, 0.9724, 0.9973, 0.9998, 0.9999,  &
0.0578, 0.1427, 0.3188, 0.5501, 0.7791, 0.9337, 0.9888, 0.9993, 0.9999,  &
0.0290, 0.0633, 0.1761, 0.3869, 0.6560, 0.8721, 0.9693, 0.9965, 0.9998,  &
0.0290, 0.0322, 0.0755, 0.2008, 0.4587, 0.7652, 0.9276, 0.9874, 0.9990,  &
0.0290, 0.0322, 0.0394, 0.0895, 0.2276, 0.5136, 0.8296, 0.9624, 0.9959 /
data ((cldnuctab( 4,iw,iconc,1, 3),iw=1,9),iconc=1,7)/    &
0.7208, 0.8903, 0.9751, 0.9969, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5862, 0.8177, 0.9496, 0.9912, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3694, 0.6892, 0.8987, 0.9800, 0.9977, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.1644, 0.4405, 0.7925, 0.9542, 0.9940, 0.9995, 0.9999, 0.9999, 0.9999,  &
0.0755, 0.2140, 0.5136, 0.8519, 0.9777, 0.9982, 0.9999, 0.9999, 0.9999,  &
0.0578, 0.1053, 0.2417, 0.4770, 0.7509, 0.9542, 0.9973, 0.9999, 0.9999,  &
0.0527, 0.0692, 0.1139, 0.1761, 0.2008, 0.2008, 0.2713, 0.3869, 0.8177 /
data ((cldnuctab( 4,iw,iconc,2, 3),iw=1,9),iconc=1,7)/    &
0.6216, 0.8296, 0.9542, 0.9922, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4953, 0.7361, 0.9066, 0.9800, 0.9980, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.3025, 0.6040, 0.8410, 0.9585, 0.9940, 0.9995, 0.9999, 0.9999, 0.9999,  &
0.1326, 0.3694, 0.7052, 0.9140, 0.9841, 0.9985, 0.9999, 0.9999, 0.9999,  &
0.0578, 0.1644, 0.4405, 0.7925, 0.9624, 0.9954, 0.9996, 0.9999, 0.9999,  &
0.0578, 0.0823, 0.2140, 0.4953, 0.8177, 0.9724, 0.9982, 0.9999, 0.9999,  &
0.0527, 0.0633, 0.1053, 0.2140, 0.3694, 0.5319, 0.7652, 0.9447, 0.9995 /
data ((cldnuctab( 4,iw,iconc,3, 3),iw=1,9),iconc=1,7)/    &
0.5501, 0.7652, 0.9210, 0.9858, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4225, 0.6560, 0.8623, 0.9624, 0.9954, 0.9996, 0.9999, 0.9999, 0.9999,  &
0.2417, 0.5136, 0.7791, 0.9276, 0.9858, 0.9987, 0.9999, 0.9999, 0.9999,  &
0.1053, 0.3025, 0.6216, 0.8623, 0.9693, 0.9959, 0.9997, 0.9999, 0.9999,  &
0.0578, 0.1326, 0.3694, 0.7208, 0.9276, 0.9888, 0.9990, 0.9999, 0.9999,  &
0.0578, 0.0692, 0.1644, 0.4225, 0.7791, 0.9624, 0.9965, 0.9998, 0.9999,  &
0.0527, 0.0633, 0.0895, 0.2008, 0.4225, 0.6892, 0.9276, 0.9947, 0.9998 /
data ((cldnuctab( 4,iw,iconc,4, 3),iw=1,9),iconc=1,7)/    &
0.4770, 0.7052, 0.8815, 0.9751, 0.9973, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.3522, 0.5862, 0.8053, 0.9394, 0.9901, 0.9993, 0.9999, 0.9999, 0.9999,  &
0.2008, 0.4405, 0.7052, 0.8903, 0.9751, 0.9969, 0.9998, 0.9999, 0.9999,  &
0.0895, 0.2417, 0.5319, 0.8053, 0.9447, 0.9912, 0.9993, 0.9999, 0.9999,  &
0.0578, 0.1053, 0.3025, 0.6390, 0.8903, 0.9777, 0.9973, 0.9998, 0.9999,  &
0.0578, 0.0692, 0.1326, 0.3522, 0.7052, 0.9394, 0.9922, 0.9994, 0.9999,  &
0.0578, 0.0633, 0.0823, 0.1644, 0.3869, 0.7052, 0.9337, 0.9940, 0.9997 /
data ((cldnuctab( 4,iw,iconc,5, 3),iw=1,9),iconc=1,7)/    &
0.4225, 0.6560, 0.8519, 0.9624, 0.9954, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.3025, 0.5319, 0.7509, 0.9140, 0.9841, 0.9987, 0.9999, 0.9999, 0.9999,  &
0.1761, 0.3869, 0.6390, 0.8519, 0.9585, 0.9947, 0.9996, 0.9999, 0.9999,  &
0.0755, 0.2140, 0.4770, 0.7509, 0.9140, 0.9841, 0.9985, 0.9999, 0.9999,  &
0.0578, 0.0895, 0.2417, 0.5501, 0.8410, 0.9624, 0.9947, 0.9996, 0.9999,  &
0.0578, 0.0692, 0.1053, 0.2867, 0.6216, 0.8987, 0.9841, 0.9985, 0.9999,  &
0.0578, 0.0633, 0.0755, 0.1427, 0.3353, 0.6560, 0.9210, 0.9912, 0.9994 /
data ((cldnuctab( 4,iw,iconc,6, 3),iw=1,9),iconc=1,7)/    &
0.3869, 0.6040, 0.8177, 0.9496, 0.9932, 0.9996, 0.9999, 0.9999, 0.9999,  &
0.2713, 0.4770, 0.7052, 0.8903, 0.9777, 0.9977, 0.9999, 0.9999, 0.9999,  &
0.1533, 0.3522, 0.5862, 0.8053, 0.9447, 0.9912, 0.9994, 0.9999, 0.9999,  &
0.0692, 0.1882, 0.4225, 0.6892, 0.8903, 0.9751, 0.9969, 0.9998, 0.9999,  &
0.0633, 0.0755, 0.2140, 0.4953, 0.7791, 0.9394, 0.9901, 0.9991, 0.9999,  &
0.0633, 0.0692, 0.0895, 0.2417, 0.5501, 0.8623, 0.9724, 0.9969, 0.9998,  &
0.0633, 0.0692, 0.0755, 0.1230, 0.2867, 0.6040, 0.8903, 0.9874, 0.9990 /
data ((cldnuctab( 4,iw,iconc,7, 3),iw=1,9),iconc=1,7)/    &
0.3522, 0.5683, 0.7925, 0.9394, 0.9912, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.2563, 0.4587, 0.6728, 0.8721, 0.9693, 0.9969, 0.9998, 0.9999, 0.9999,  &
0.1427, 0.3188, 0.5501, 0.7791, 0.9276, 0.9874, 0.9990, 0.9999, 0.9999,  &
0.0692, 0.1644, 0.3869, 0.6560, 0.8623, 0.9660, 0.9954, 0.9997, 0.9999,  &
0.0692, 0.0755, 0.1882, 0.4405, 0.7361, 0.9210, 0.9858, 0.9987, 0.9999,  &
0.0692, 0.0692, 0.0823, 0.2140, 0.4953, 0.8177, 0.9585, 0.9947, 0.9996,  &
0.0633, 0.0692, 0.0823, 0.1053, 0.2563, 0.5501, 0.8623, 0.9800, 0.9982 /
data ((cldnuctab( 4,iw,iconc,1, 4),iw=1,9),iconc=1,7)/    &
0.8053, 0.9394, 0.9888, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6560, 0.8815, 0.9751, 0.9969, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3869, 0.7509, 0.9394, 0.9912, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.1761, 0.4587, 0.8296, 0.9751, 0.9973, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.1230, 0.2276, 0.5136, 0.8623, 0.9874, 0.9991, 0.9999, 0.9999, 0.9999,  &
0.1139, 0.1427, 0.2563, 0.4587, 0.7052, 0.9337, 0.9973, 0.9999, 0.9999,  &
0.0971, 0.1139, 0.1533, 0.2008, 0.2140, 0.2140, 0.2867, 0.3694, 0.7509 /
data ((cldnuctab( 4,iw,iconc,2, 4),iw=1,9),iconc=1,7)/    &
0.7208, 0.8987, 0.9777, 0.9973, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5683, 0.8177, 0.9496, 0.9922, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3188, 0.6728, 0.8987, 0.9800, 0.9977, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.1427, 0.3869, 0.7509, 0.9496, 0.9932, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.1139, 0.1882, 0.4587, 0.8296, 0.9777, 0.9980, 0.9998, 0.9999, 0.9999,  &
0.1139, 0.1326, 0.2417, 0.4770, 0.8177, 0.9800, 0.9990, 0.9999, 0.9999,  &
0.1053, 0.1139, 0.1427, 0.2276, 0.3522, 0.4953, 0.7052, 0.9276, 0.9996 /
data ((cldnuctab( 4,iw,iconc,3, 4),iw=1,9),iconc=1,7)/    &
0.6390, 0.8519, 0.9585, 0.9940, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4953, 0.7509, 0.9210, 0.9841, 0.9982, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2713, 0.5862, 0.8410, 0.9624, 0.9947, 0.9995, 0.9999, 0.9999, 0.9999,  &
0.1139, 0.3188, 0.6728, 0.9140, 0.9858, 0.9985, 0.9999, 0.9999, 0.9999,  &
0.1139, 0.1427, 0.3694, 0.7509, 0.9585, 0.9954, 0.9996, 0.9999, 0.9999,  &
0.1139, 0.1326, 0.2008, 0.4405, 0.7925, 0.9751, 0.9982, 0.9999, 0.9999,  &
0.1053, 0.1139, 0.1427, 0.2276, 0.4046, 0.6560, 0.9210, 0.9959, 0.9999 /
data ((cldnuctab( 4,iw,iconc,4, 4),iw=1,9),iconc=1,7)/    &
0.5862, 0.7925, 0.9394, 0.9888, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4225, 0.6892, 0.8815, 0.9724, 0.9965, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.2276, 0.5136, 0.7925, 0.9394, 0.9888, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.1230, 0.2713, 0.5862, 0.8721, 0.9724, 0.9965, 0.9997, 0.9999, 0.9999,  &
0.1230, 0.1326, 0.3188, 0.6728, 0.9276, 0.9888, 0.9990, 0.9999, 0.9999,  &
0.1139, 0.1326, 0.1644, 0.3694, 0.7361, 0.9585, 0.9965, 0.9998, 0.9999,  &
0.1139, 0.1230, 0.1427, 0.2008, 0.3869, 0.6728, 0.9337, 0.9965, 0.9998 /
data ((cldnuctab( 4,iw,iconc,5, 4),iw=1,9),iconc=1,7)/    &
0.5136, 0.7509, 0.9140, 0.9841, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3869, 0.6216, 0.8410, 0.9585, 0.9940, 0.9995, 0.9999, 0.9999, 0.9999,  &
0.2008, 0.4587, 0.7361, 0.9066, 0.9821, 0.9980, 0.9998, 0.9999, 0.9999,  &
0.1230, 0.2276, 0.5136, 0.8177, 0.9542, 0.9932, 0.9994, 0.9999, 0.9999,  &
0.1230, 0.1326, 0.2563, 0.5862, 0.8815, 0.9800, 0.9980, 0.9998, 0.9999,  &
0.1230, 0.1326, 0.1533, 0.3025, 0.6560, 0.9337, 0.9922, 0.9994, 0.9999,  &
0.1139, 0.1230, 0.1427, 0.1882, 0.3522, 0.6390, 0.9276, 0.9954, 0.9998 /
data ((cldnuctab( 4,iw,iconc,6, 4),iw=1,9),iconc=1,7)/    &
0.4770, 0.7052, 0.8903, 0.9777, 0.9977, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.3353, 0.5862, 0.8053, 0.9394, 0.9901, 0.9993, 0.9999, 0.9999, 0.9999,  &
0.1761, 0.4046, 0.6892, 0.8815, 0.9724, 0.9969, 0.9998, 0.9999, 0.9999,  &
0.1230, 0.2008, 0.4770, 0.7652, 0.9337, 0.9888, 0.9990, 0.9999, 0.9999,  &
0.1230, 0.1427, 0.2276, 0.5319, 0.8410, 0.9693, 0.9959, 0.9997, 0.9999,  &
0.1230, 0.1326, 0.1533, 0.2713, 0.5862, 0.8987, 0.9858, 0.9989, 0.9999,  &
0.1230, 0.1326, 0.1427, 0.1761, 0.3025, 0.5862, 0.9066, 0.9922, 0.9996 /
data ((cldnuctab( 4,iw,iconc,7, 4),iw=1,9),iconc=1,7)/    &
0.4587, 0.6728, 0.8721, 0.9724, 0.9969, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.3188, 0.5501, 0.7791, 0.9276, 0.9874, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.1644, 0.3694, 0.6390, 0.8519, 0.9624, 0.9954, 0.9997, 0.9999, 0.9999,  &
0.1326, 0.1761, 0.4225, 0.7361, 0.9140, 0.9841, 0.9985, 0.9999, 0.9999,  &
0.1326, 0.1427, 0.2008, 0.4770, 0.8053, 0.9542, 0.9940, 0.9995, 0.9999,  &
0.1326, 0.1427, 0.1533, 0.2417, 0.5319, 0.8623, 0.9777, 0.9980, 0.9998,  &
0.1326, 0.1326, 0.1427, 0.1761, 0.2867, 0.5501, 0.8721, 0.9888, 0.9993 /
data ((cldnuctab( 4,iw,iconc,1, 5),iw=1,9),iconc=1,7)/    &
0.8410, 0.9585, 0.9940, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6892, 0.9140, 0.9841, 0.9982, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4046, 0.7791, 0.9585, 0.9947, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2008, 0.4770, 0.8519, 0.9821, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.1644, 0.2563, 0.5136, 0.8721, 0.9901, 0.9995, 0.9999, 0.9999, 0.9999,  &
0.1533, 0.1761, 0.2713, 0.4405, 0.6728, 0.9140, 0.9969, 0.9999, 0.9999,  &
0.1326, 0.1533, 0.1882, 0.2417, 0.2140, 0.2140, 0.2713, 0.3869, 0.7509 /
data ((cldnuctab( 4,iw,iconc,2, 5),iw=1,9),iconc=1,7)/    &
0.7791, 0.9276, 0.9858, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6040, 0.8623, 0.9660, 0.9954, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3353, 0.6892, 0.9210, 0.9874, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.1761, 0.3869, 0.7791, 0.9624, 0.9959, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.1644, 0.2140, 0.4587, 0.8410, 0.9841, 0.9989, 0.9999, 0.9999, 0.9999,  &
0.1533, 0.1761, 0.2563, 0.4770, 0.8053, 0.9800, 0.9993, 0.9999, 0.9999,  &
0.1427, 0.1533, 0.1882, 0.2563, 0.3353, 0.4587, 0.6728, 0.9066, 0.9991 /
data ((cldnuctab( 4,iw,iconc,3, 5),iw=1,9),iconc=1,7)/    &
0.7052, 0.8903, 0.9751, 0.9969, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5319, 0.8053, 0.9447, 0.9901, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2713, 0.6216, 0.8815, 0.9751, 0.9969, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.1644, 0.3188, 0.6892, 0.9337, 0.9912, 0.9991, 0.9999, 0.9999, 0.9999,  &
0.1644, 0.1882, 0.3869, 0.7652, 0.9660, 0.9969, 0.9998, 0.9999, 0.9999,  &
0.1644, 0.1761, 0.2276, 0.4405, 0.7925, 0.9800, 0.9990, 0.9999, 0.9999,  &
0.1533, 0.1644, 0.1882, 0.2563, 0.4046, 0.6216, 0.8987, 0.9954, 0.9999 /
data ((cldnuctab( 4,iw,iconc,4, 5),iw=1,9),iconc=1,7)/    &
0.6390, 0.8410, 0.9585, 0.9940, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4587, 0.7361, 0.9140, 0.9821, 0.9980, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.2417, 0.5319, 0.8296, 0.9585, 0.9940, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.1761, 0.2713, 0.6216, 0.8903, 0.9821, 0.9980, 0.9998, 0.9999, 0.9999,  &
0.1644, 0.1882, 0.3188, 0.6892, 0.9394, 0.9932, 0.9994, 0.9999, 0.9999,  &
0.1644, 0.1761, 0.2140, 0.3869, 0.7361, 0.9660, 0.9977, 0.9998, 0.9999,  &
0.1644, 0.1644, 0.1882, 0.2417, 0.3869, 0.6728, 0.9276, 0.9969, 0.9999 /
data ((cldnuctab( 4,iw,iconc,5, 5),iw=1,9),iconc=1,7)/    &
0.5862, 0.8053, 0.9394, 0.9901, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4225, 0.6892, 0.8815, 0.9724, 0.9965, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.2140, 0.4770, 0.7791, 0.9337, 0.9888, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.1761, 0.2276, 0.5501, 0.8519, 0.9693, 0.9959, 0.9997, 0.9999, 0.9999,  &
0.1761, 0.1882, 0.2713, 0.6040, 0.9066, 0.9874, 0.9989, 0.9999, 0.9999,  &
0.1761, 0.1882, 0.2140, 0.3353, 0.6560, 0.9394, 0.9954, 0.9997, 0.9999,  &
0.1644, 0.1761, 0.1882, 0.2276, 0.3694, 0.6390, 0.9210, 0.9965, 0.9998 /
data ((cldnuctab( 4,iw,iconc,6, 5),iw=1,9),iconc=1,7)/    &
0.5501, 0.7652, 0.9210, 0.9858, 0.9989, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3869, 0.6390, 0.8519, 0.9624, 0.9947, 0.9996, 0.9999, 0.9999, 0.9999,  &
0.1882, 0.4405, 0.7361, 0.9140, 0.9821, 0.9982, 0.9999, 0.9999, 0.9999,  &
0.1761, 0.2008, 0.4953, 0.8053, 0.9542, 0.9932, 0.9994, 0.9999, 0.9999,  &
0.1761, 0.1882, 0.2417, 0.5501, 0.8623, 0.9777, 0.9980, 0.9998, 0.9999,  &
0.1761, 0.1882, 0.2140, 0.2867, 0.5862, 0.9066, 0.9912, 0.9994, 0.9999,  &
0.1761, 0.1761, 0.1882, 0.2276, 0.3353, 0.5862, 0.8987, 0.9947, 0.9997 /
data ((cldnuctab( 4,iw,iconc,7, 5),iw=1,9),iconc=1,7)/    &
0.5136, 0.7361, 0.9066, 0.9821, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3522, 0.6040, 0.8177, 0.9496, 0.9932, 0.9995, 0.9999, 0.9999, 0.9999,  &
0.1882, 0.4046, 0.6892, 0.8903, 0.9777, 0.9973, 0.9998, 0.9999, 0.9999,  &
0.1882, 0.2008, 0.4405, 0.7652, 0.9394, 0.9901, 0.9991, 0.9999, 0.9999,  &
0.1882, 0.2008, 0.2276, 0.4953, 0.8296, 0.9693, 0.9965, 0.9998, 0.9999,  &
0.1882, 0.2008, 0.2140, 0.2713, 0.5319, 0.8721, 0.9858, 0.9989, 0.9999,  &
0.1761, 0.1882, 0.2008, 0.2276, 0.3188, 0.5501, 0.8721, 0.9912, 0.9996 /
data ((cldnuctab( 4,iw,iconc,1, 6),iw=1,9),iconc=1,7)/    &
0.8721, 0.9693, 0.9959, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7052, 0.9276, 0.9888, 0.9989, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4046, 0.7925, 0.9660, 0.9965, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2276, 0.4770, 0.8519, 0.9858, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2008, 0.2713, 0.5136, 0.8623, 0.9912, 0.9996, 0.9999, 0.9999, 0.9999,  &
0.1882, 0.2140, 0.2867, 0.4405, 0.6560, 0.8987, 0.9965, 0.9999, 0.9999,  &
0.1644, 0.1761, 0.2140, 0.2867, 0.2417, 0.2417, 0.2713, 0.4046, 0.7052 /
data ((cldnuctab( 4,iw,iconc,2, 6),iw=1,9),iconc=1,7)/    &
0.8053, 0.9447, 0.9901, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6216, 0.8815, 0.9751, 0.9969, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3353, 0.7052, 0.9394, 0.9912, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2140, 0.4046, 0.7925, 0.9693, 0.9973, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.2008, 0.2417, 0.4770, 0.8410, 0.9858, 0.9993, 0.9999, 0.9999, 0.9999,  &
0.1882, 0.2140, 0.2867, 0.4770, 0.7925, 0.9800, 0.9994, 0.9999, 0.9999,  &
0.1761, 0.1882, 0.2140, 0.2867, 0.3522, 0.4405, 0.6560, 0.8903, 0.9989 /
data ((cldnuctab( 4,iw,iconc,3, 6),iw=1,9),iconc=1,7)/    &
0.7361, 0.9066, 0.9821, 0.9980, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5501, 0.8296, 0.9542, 0.9932, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2867, 0.6390, 0.8987, 0.9821, 0.9980, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.2140, 0.3353, 0.7052, 0.9447, 0.9932, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.2140, 0.2417, 0.4046, 0.7791, 0.9724, 0.9980, 0.9998, 0.9999, 0.9999,  &
0.2008, 0.2140, 0.2713, 0.4405, 0.7925, 0.9821, 0.9993, 0.9999, 0.9999,  &
0.1882, 0.2008, 0.2140, 0.2713, 0.4046, 0.6040, 0.8815, 0.9947, 0.9999 /
data ((cldnuctab( 4,iw,iconc,4, 6),iw=1,9),iconc=1,7)/    &
0.6728, 0.8721, 0.9693, 0.9959, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4953, 0.7652, 0.9276, 0.9874, 0.9989, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2417, 0.5501, 0.8519, 0.9660, 0.9959, 0.9996, 0.9999, 0.9999, 0.9999,  &
0.2140, 0.2713, 0.6216, 0.9066, 0.9858, 0.9987, 0.9999, 0.9999, 0.9999,  &
0.2140, 0.2276, 0.3353, 0.6892, 0.9496, 0.9954, 0.9996, 0.9999, 0.9999,  &
0.2008, 0.2140, 0.2563, 0.4046, 0.7361, 0.9693, 0.9985, 0.9999, 0.9999,  &
0.2008, 0.2140, 0.2276, 0.2713, 0.4046, 0.6390, 0.9140, 0.9973, 0.9999 /
data ((cldnuctab( 4,iw,iconc,5, 6),iw=1,9),iconc=1,7)/    &
0.6216, 0.8296, 0.9542, 0.9932, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4405, 0.7208, 0.8987, 0.9800, 0.9980, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.2276, 0.4953, 0.8053, 0.9496, 0.9922, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.2140, 0.2417, 0.5501, 0.8721, 0.9751, 0.9973, 0.9998, 0.9999, 0.9999,  &
0.2140, 0.2417, 0.2867, 0.6216, 0.9210, 0.9901, 0.9993, 0.9999, 0.9999,  &
0.2140, 0.2276, 0.2563, 0.3522, 0.6560, 0.9496, 0.9965, 0.9998, 0.9999,  &
0.2140, 0.2140, 0.2276, 0.2563, 0.3869, 0.6216, 0.9140, 0.9969, 0.9999 /
data ((cldnuctab( 4,iw,iconc,6, 6),iw=1,9),iconc=1,7)/    &
0.5862, 0.8053, 0.9394, 0.9912, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4046, 0.6728, 0.8721, 0.9724, 0.9965, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.2276, 0.4587, 0.7652, 0.9276, 0.9874, 0.9989, 0.9999, 0.9999, 0.9999,  &
0.2276, 0.2417, 0.4953, 0.8296, 0.9624, 0.9954, 0.9997, 0.9999, 0.9999,  &
0.2276, 0.2417, 0.2713, 0.5501, 0.8815, 0.9841, 0.9987, 0.9999, 0.9999,  &
0.2276, 0.2276, 0.2563, 0.3188, 0.6040, 0.9140, 0.9932, 0.9996, 0.9999,  &
0.2140, 0.2276, 0.2417, 0.2713, 0.3522, 0.5862, 0.8903, 0.9954, 0.9998 /
data ((cldnuctab( 4,iw,iconc,7, 6),iw=1,9),iconc=1,7)/    &
0.5501, 0.7791, 0.9276, 0.9888, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3694, 0.6390, 0.8519, 0.9624, 0.9954, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.2417, 0.4225, 0.7208, 0.9140, 0.9841, 0.9985, 0.9999, 0.9999, 0.9999,  &
0.2417, 0.2563, 0.4587, 0.7925, 0.9496, 0.9932, 0.9995, 0.9999, 0.9999,  &
0.2276, 0.2417, 0.2713, 0.4953, 0.8410, 0.9751, 0.9977, 0.9998, 0.9999,  &
0.2276, 0.2417, 0.2563, 0.3025, 0.5501, 0.8815, 0.9888, 0.9993, 0.9999,  &
0.2276, 0.2276, 0.2417, 0.2713, 0.3353, 0.5501, 0.8623, 0.9922, 0.9997 /
data ((cldnuctab( 4,iw,iconc,1, 7),iw=1,9),iconc=1,7)/    &
0.8815, 0.9751, 0.9969, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7208, 0.9394, 0.9912, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4046, 0.8053, 0.9724, 0.9973, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2563, 0.4770, 0.8623, 0.9888, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2417, 0.3025, 0.5136, 0.8519, 0.9922, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.2140, 0.2417, 0.3188, 0.4405, 0.6390, 0.8815, 0.9959, 0.9999, 0.9999,  &
0.1882, 0.2008, 0.2276, 0.3025, 0.2713, 0.2417, 0.2867, 0.3869, 0.7052 /
data ((cldnuctab( 4,iw,iconc,2, 7),iw=1,9),iconc=1,7)/    &
0.8296, 0.9542, 0.9932, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6390, 0.8987, 0.9800, 0.9980, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3353, 0.7208, 0.9447, 0.9932, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2563, 0.4046, 0.7925, 0.9751, 0.9980, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.2417, 0.2867, 0.4770, 0.8410, 0.9888, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.2276, 0.2417, 0.3025, 0.4770, 0.7652, 0.9777, 0.9995, 0.9999, 0.9999,  &
0.2140, 0.2276, 0.2417, 0.3025, 0.3694, 0.4405, 0.6216, 0.8721, 0.9985 /
data ((cldnuctab( 4,iw,iconc,3, 7),iw=1,9),iconc=1,7)/    &
0.7652, 0.9276, 0.9858, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5683, 0.8519, 0.9660, 0.9954, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2867, 0.6390, 0.9140, 0.9858, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2563, 0.3353, 0.7208, 0.9542, 0.9947, 0.9996, 0.9999, 0.9999, 0.9999,  &
0.2417, 0.2713, 0.4046, 0.7791, 0.9777, 0.9985, 0.9999, 0.9999, 0.9999,  &
0.2276, 0.2563, 0.2867, 0.4587, 0.7791, 0.9821, 0.9994, 0.9999, 0.9999,  &
0.2276, 0.2276, 0.2563, 0.3025, 0.4225, 0.6040, 0.8623, 0.9940, 0.9999 /
data ((cldnuctab( 4,iw,iconc,4, 7),iw=1,9),iconc=1,7)/    &
0.7052, 0.8903, 0.9751, 0.9973, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4953, 0.7925, 0.9447, 0.9912, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2563, 0.5683, 0.8623, 0.9724, 0.9969, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.2563, 0.2867, 0.6390, 0.9210, 0.9901, 0.9991, 0.9999, 0.9999, 0.9999,  &
0.2563, 0.2713, 0.3522, 0.6892, 0.9542, 0.9965, 0.9998, 0.9999, 0.9999,  &
0.2417, 0.2563, 0.2867, 0.4046, 0.7208, 0.9724, 0.9989, 0.9999, 0.9999,  &
0.2276, 0.2417, 0.2563, 0.3025, 0.4225, 0.6216, 0.9066, 0.9973, 0.9999 /
data ((cldnuctab( 4,iw,iconc,5, 7),iw=1,9),iconc=1,7)/    &
0.6560, 0.8623, 0.9660, 0.9954, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4587, 0.7509, 0.9210, 0.9841, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2563, 0.5136, 0.8177, 0.9585, 0.9940, 0.9995, 0.9999, 0.9999, 0.9999,  &
0.2563, 0.2867, 0.5683, 0.8815, 0.9821, 0.9982, 0.9998, 0.9999, 0.9999,  &
0.2563, 0.2713, 0.3025, 0.6216, 0.9276, 0.9922, 0.9994, 0.9999, 0.9999,  &
0.2563, 0.2563, 0.2867, 0.3694, 0.6560, 0.9496, 0.9973, 0.9998, 0.9999,  &
0.2417, 0.2563, 0.2713, 0.3025, 0.4046, 0.6216, 0.9140, 0.9969, 0.9999 /
data ((cldnuctab( 4,iw,iconc,6, 7),iw=1,9),iconc=1,7)/    &
0.6040, 0.8296, 0.9542, 0.9932, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4225, 0.7052, 0.8903, 0.9777, 0.9977, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.2713, 0.4587, 0.7791, 0.9394, 0.9912, 0.9993, 0.9999, 0.9999, 0.9999,  &
0.2713, 0.2867, 0.5136, 0.8410, 0.9693, 0.9969, 0.9998, 0.9999, 0.9999,  &
0.2713, 0.2867, 0.3188, 0.5501, 0.8903, 0.9874, 0.9990, 0.9999, 0.9999,  &
0.2563, 0.2713, 0.2867, 0.3522, 0.6040, 0.9210, 0.9947, 0.9997, 0.9999,  &
0.2563, 0.2563, 0.2713, 0.3025, 0.3694, 0.5862, 0.8815, 0.9959, 0.9998 /
data ((cldnuctab( 4,iw,iconc,7, 7),iw=1,9),iconc=1,7)/    &
0.5862, 0.8053, 0.9394, 0.9912, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3869, 0.6728, 0.8721, 0.9724, 0.9969, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.2713, 0.4225, 0.7509, 0.9276, 0.9874, 0.9989, 0.9999, 0.9999, 0.9999,  &
0.2713, 0.2867, 0.4587, 0.8053, 0.9585, 0.9954, 0.9996, 0.9999, 0.9999,  &
0.2713, 0.2867, 0.3188, 0.5136, 0.8519, 0.9800, 0.9982, 0.9999, 0.9999,  &
0.2713, 0.2713, 0.3025, 0.3522, 0.5501, 0.8903, 0.9912, 0.9994, 0.9999,  &
0.2563, 0.2713, 0.2867, 0.3025, 0.3694, 0.5501, 0.8519, 0.9932, 0.9998 /
data ((cldnuctab( 5,iw,iconc,1, 1),iw=1,9),iconc=1,7)/    &
0.8084, 0.9459, 0.9925, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6601, 0.8924, 0.9783, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3912, 0.7544, 0.9459, 0.9934, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.1790, 0.4632, 0.8437, 0.9783, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.1253, 0.2310, 0.5181, 0.8744, 0.9904, 0.9995, 0.9999, 0.9999, 0.9999,  &
0.1161, 0.1452, 0.2599, 0.4632, 0.7246, 0.9407, 0.9980, 0.9999, 0.9999,  &
0.0991, 0.1161, 0.1560, 0.2040, 0.2173, 0.2173, 0.2906, 0.3737, 0.7686 /
data ((cldnuctab( 5,iw,iconc,2, 1),iw=1,9),iconc=1,7)/    &
0.7246, 0.9085, 0.9806, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5727, 0.8325, 0.9553, 0.9942, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3228, 0.6768, 0.9085, 0.9845, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.1452, 0.3912, 0.7686, 0.9553, 0.9949, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.1253, 0.1912, 0.4632, 0.8325, 0.9806, 0.9989, 0.9999, 0.9999, 0.9999,  &
0.1161, 0.1350, 0.2310, 0.4815, 0.8207, 0.9826, 0.9994, 0.9999, 0.9999,  &
0.1073, 0.1161, 0.1452, 0.2310, 0.3564, 0.4998, 0.7091, 0.9351, 0.9995 /
data ((cldnuctab( 5,iw,iconc,3, 1),iw=1,9),iconc=1,7)/    &
0.6432, 0.8545, 0.9633, 0.9955, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4998, 0.7544, 0.9227, 0.9862, 0.9989, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2750, 0.5906, 0.8545, 0.9668, 0.9961, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.1253, 0.3228, 0.6768, 0.9227, 0.9877, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.1161, 0.1452, 0.3737, 0.7544, 0.9633, 0.9966, 0.9998, 0.9999, 0.9999,  &
0.1161, 0.1350, 0.1912, 0.4450, 0.7957, 0.9783, 0.9990, 0.9999, 0.9999,  &
0.1073, 0.1253, 0.1452, 0.2310, 0.4090, 0.6601, 0.9227, 0.9970, 0.9999 /
data ((cldnuctab( 5,iw,iconc,4, 1),iw=1,9),iconc=1,7)/    &
0.5906, 0.8084, 0.9407, 0.9915, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4269, 0.6931, 0.8837, 0.9758, 0.9974, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2310, 0.5181, 0.7957, 0.9407, 0.9915, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.1253, 0.2599, 0.5906, 0.8744, 0.9758, 0.9974, 0.9999, 0.9999, 0.9999,  &
0.1253, 0.1350, 0.3065, 0.6768, 0.9291, 0.9915, 0.9995, 0.9999, 0.9999,  &
0.1161, 0.1350, 0.1672, 0.3737, 0.7397, 0.9633, 0.9974, 0.9999, 0.9999,  &
0.1161, 0.1253, 0.1452, 0.2040, 0.3912, 0.6768, 0.9407, 0.9974, 0.9999 /
data ((cldnuctab( 5,iw,iconc,5, 1),iw=1,9),iconc=1,7)/    &
0.5181, 0.7544, 0.9158, 0.9862, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3737, 0.6432, 0.8437, 0.9595, 0.9955, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.2040, 0.4450, 0.7397, 0.9158, 0.9845, 0.9987, 0.9999, 0.9999, 0.9999,  &
0.1253, 0.2310, 0.5181, 0.8207, 0.9595, 0.9949, 0.9997, 0.9999, 0.9999,  &
0.1253, 0.1350, 0.2599, 0.5906, 0.8924, 0.9826, 0.9985, 0.9999, 0.9999,  &
0.1253, 0.1350, 0.1560, 0.3065, 0.6601, 0.9351, 0.9942, 0.9997, 0.9999,  &
0.1253, 0.1253, 0.1452, 0.1912, 0.3564, 0.6432, 0.9291, 0.9966, 0.9999 /
data ((cldnuctab( 5,iw,iconc,6, 1),iw=1,9),iconc=1,7)/    &
0.4815, 0.7091, 0.8924, 0.9806, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3395, 0.5906, 0.8084, 0.9459, 0.9925, 0.9995, 0.9999, 0.9999, 0.9999,  &
0.1790, 0.4090, 0.6931, 0.8837, 0.9758, 0.9977, 0.9999, 0.9999, 0.9999,  &
0.1350, 0.2040, 0.4632, 0.7686, 0.9407, 0.9904, 0.9994, 0.9999, 0.9999,  &
0.1350, 0.1452, 0.2310, 0.5364, 0.8437, 0.9701, 0.9970, 0.9998, 0.9999,  &
0.1253, 0.1350, 0.1560, 0.2750, 0.5906, 0.9007, 0.9891, 0.9993, 0.9999,  &
0.1253, 0.1350, 0.1452, 0.1790, 0.3065, 0.5906, 0.9085, 0.9942, 0.9997 /
data ((cldnuctab( 5,iw,iconc,7, 1),iw=1,9),iconc=1,7)/    &
0.4632, 0.6931, 0.8744, 0.9758, 0.9977, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3228, 0.5546, 0.7824, 0.9291, 0.9891, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.1672, 0.3737, 0.6432, 0.8545, 0.9668, 0.9966, 0.9998, 0.9999, 0.9999,  &
0.1350, 0.1790, 0.4269, 0.7397, 0.9158, 0.9862, 0.9989, 0.9999, 0.9999,  &
0.1350, 0.1452, 0.2040, 0.4815, 0.8084, 0.9595, 0.9949, 0.9997, 0.9999,  &
0.1350, 0.1452, 0.1560, 0.2452, 0.5364, 0.8647, 0.9806, 0.9985, 0.9999,  &
0.1350, 0.1350, 0.1560, 0.1790, 0.2906, 0.5546, 0.8837, 0.9904, 0.9995 /
data ((cldnuctab( 5,iw,iconc,1, 2),iw=1,9),iconc=1,7)/    &
0.8744, 0.9731, 0.9966, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7091, 0.9291, 0.9891, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4090, 0.7957, 0.9701, 0.9970, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2310, 0.4815, 0.8545, 0.9877, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2040, 0.2750, 0.5181, 0.8647, 0.9925, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.1912, 0.2173, 0.2906, 0.4269, 0.6601, 0.9007, 0.9974, 0.9999, 0.9999,  &
0.1672, 0.1790, 0.2173, 0.2906, 0.2452, 0.2452, 0.2750, 0.4090, 0.7091 /
data ((cldnuctab( 5,iw,iconc,2, 2),iw=1,9),iconc=1,7)/    &
0.8084, 0.9459, 0.9915, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6259, 0.8837, 0.9783, 0.9974, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3395, 0.7091, 0.9407, 0.9925, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2173, 0.4090, 0.7957, 0.9731, 0.9977, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2040, 0.2452, 0.4632, 0.8437, 0.9877, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.1912, 0.2173, 0.2906, 0.4815, 0.7957, 0.9826, 0.9996, 0.9999, 0.9999,  &
0.1790, 0.1912, 0.2173, 0.2906, 0.3564, 0.4450, 0.6601, 0.8924, 0.9992 /
data ((cldnuctab( 5,iw,iconc,3, 2),iw=1,9),iconc=1,7)/    &
0.7397, 0.9158, 0.9826, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5546, 0.8325, 0.9595, 0.9942, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2906, 0.6259, 0.9007, 0.9826, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2173, 0.3395, 0.7091, 0.9459, 0.9942, 0.9996, 0.9999, 0.9999, 0.9999,  &
0.2040, 0.2452, 0.3912, 0.7824, 0.9758, 0.9983, 0.9999, 0.9999, 0.9999,  &
0.2040, 0.2173, 0.2599, 0.4450, 0.7957, 0.9826, 0.9995, 0.9999, 0.9999,  &
0.1912, 0.2040, 0.2173, 0.2750, 0.4090, 0.6084, 0.8837, 0.9955, 0.9999 /
data ((cldnuctab( 5,iw,iconc,4, 2),iw=1,9),iconc=1,7)/    &
0.6768, 0.8744, 0.9701, 0.9966, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4815, 0.7686, 0.9351, 0.9891, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2452, 0.5546, 0.8545, 0.9701, 0.9961, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.2173, 0.2750, 0.6259, 0.9158, 0.9877, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.2173, 0.2310, 0.3395, 0.6931, 0.9508, 0.9961, 0.9997, 0.9999, 0.9999,  &
0.2040, 0.2173, 0.2599, 0.3912, 0.7397, 0.9701, 0.9987, 0.9999, 0.9999,  &
0.2040, 0.2040, 0.2310, 0.2750, 0.4090, 0.6432, 0.9227, 0.9977, 0.9999 /
data ((cldnuctab( 5,iw,iconc,5, 2),iw=1,9),iconc=1,7)/    &
0.6259, 0.8325, 0.9553, 0.9942, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4450, 0.7246, 0.9007, 0.9806, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.2173, 0.4998, 0.8084, 0.9508, 0.9934, 0.9995, 0.9999, 0.9999, 0.9999,  &
0.2173, 0.2452, 0.5546, 0.8744, 0.9783, 0.9977, 0.9999, 0.9999, 0.9999,  &
0.2173, 0.2310, 0.2906, 0.6084, 0.9227, 0.9915, 0.9994, 0.9999, 0.9999,  &
0.2173, 0.2310, 0.2599, 0.3564, 0.6601, 0.9508, 0.9970, 0.9998, 0.9999,  &
0.2040, 0.2173, 0.2310, 0.2599, 0.3737, 0.6259, 0.9158, 0.9974, 0.9999 /
data ((cldnuctab( 5,iw,iconc,6, 2),iw=1,9),iconc=1,7)/    &
0.5906, 0.8084, 0.9407, 0.9915, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4090, 0.6768, 0.8744, 0.9731, 0.9970, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.2310, 0.4450, 0.7544, 0.9291, 0.9891, 0.9992, 0.9999, 0.9999, 0.9999,  &
0.2310, 0.2452, 0.4998, 0.8325, 0.9668, 0.9961, 0.9997, 0.9999, 0.9999,  &
0.2310, 0.2452, 0.2750, 0.5546, 0.8837, 0.9845, 0.9989, 0.9999, 0.9999,  &
0.2173, 0.2310, 0.2599, 0.3228, 0.5906, 0.9158, 0.9942, 0.9997, 0.9999,  &
0.2173, 0.2310, 0.2310, 0.2599, 0.3564, 0.5906, 0.8924, 0.9961, 0.9999 /
data ((cldnuctab( 5,iw,iconc,7, 2),iw=1,9),iconc=1,7)/    &
0.5546, 0.7824, 0.9291, 0.9891, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3737, 0.6432, 0.8545, 0.9633, 0.9961, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.2452, 0.4269, 0.7246, 0.9158, 0.9845, 0.9987, 0.9999, 0.9999, 0.9999,  &
0.2310, 0.2452, 0.4632, 0.7957, 0.9508, 0.9942, 0.9996, 0.9999, 0.9999,  &
0.2310, 0.2452, 0.2750, 0.4998, 0.8437, 0.9758, 0.9980, 0.9999, 0.9999,  &
0.2310, 0.2452, 0.2599, 0.3065, 0.5546, 0.8837, 0.9891, 0.9994, 0.9999,  &
0.2310, 0.2310, 0.2452, 0.2750, 0.3395, 0.5546, 0.8647, 0.9934, 0.9998 /
data ((cldnuctab( 5,iw,iconc,1, 3),iw=1,9),iconc=1,7)/    &
0.9158, 0.9862, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7397, 0.9595, 0.9955, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4269, 0.8084, 0.9806, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3395, 0.4998, 0.8545, 0.9925, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3065, 0.3564, 0.5364, 0.8325, 0.9925, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.2906, 0.3065, 0.3737, 0.4632, 0.5906, 0.8437, 0.9925, 0.9999, 0.9999,  &
0.2750, 0.2750, 0.3065, 0.3737, 0.4090, 0.2599, 0.3065, 0.4090, 0.6768 /
data ((cldnuctab( 5,iw,iconc,2, 3),iw=1,9),iconc=1,7)/    &
0.8647, 0.9731, 0.9966, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6601, 0.9227, 0.9891, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3564, 0.7397, 0.9633, 0.9966, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3395, 0.4269, 0.7957, 0.9826, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3228, 0.3564, 0.4998, 0.8325, 0.9915, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.3065, 0.3228, 0.3737, 0.4998, 0.7544, 0.9701, 0.9997, 0.9999, 0.9999,  &
0.2750, 0.2906, 0.3065, 0.3737, 0.4632, 0.4450, 0.5906, 0.7957, 0.9915 /
data ((cldnuctab( 5,iw,iconc,3, 3),iw=1,9),iconc=1,7)/    &
0.8084, 0.9508, 0.9925, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5906, 0.8837, 0.9783, 0.9977, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3395, 0.6601, 0.9351, 0.9925, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3395, 0.3737, 0.7246, 0.9668, 0.9974, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.3228, 0.3564, 0.4450, 0.7824, 0.9826, 0.9993, 0.9999, 0.9999, 0.9999,  &
0.3228, 0.3228, 0.3564, 0.4815, 0.7544, 0.9806, 0.9997, 0.9999, 0.9999,  &
0.3065, 0.3065, 0.3228, 0.3737, 0.4632, 0.5546, 0.8207, 0.9806, 0.9999 /
data ((cldnuctab( 5,iw,iconc,4, 3),iw=1,9),iconc=1,7)/    &
0.7544, 0.9291, 0.9862, 0.9989, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5364, 0.8325, 0.9633, 0.9955, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3564, 0.5906, 0.8924, 0.9845, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3395, 0.3737, 0.6432, 0.9351, 0.9942, 0.9996, 0.9999, 0.9999, 0.9999,  &
0.3395, 0.3564, 0.4090, 0.6931, 0.9633, 0.9980, 0.9999, 0.9999, 0.9999,  &
0.3228, 0.3395, 0.3737, 0.4632, 0.7246, 0.9731, 0.9993, 0.9999, 0.9999,  &
0.3228, 0.3228, 0.3395, 0.3737, 0.4632, 0.6259, 0.8647, 0.9966, 0.9999 /
data ((cldnuctab( 5,iw,iconc,5, 3),iw=1,9),iconc=1,7)/    &
0.7091, 0.9007, 0.9806, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4815, 0.7957, 0.9459, 0.9915, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3564, 0.5181, 0.8545, 0.9731, 0.9974, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.3564, 0.3737, 0.5727, 0.9007, 0.9891, 0.9992, 0.9999, 0.9999, 0.9999,  &
0.3395, 0.3564, 0.4090, 0.6259, 0.9351, 0.9955, 0.9998, 0.9999, 0.9999,  &
0.3395, 0.3564, 0.3737, 0.4269, 0.6601, 0.9508, 0.9985, 0.9999, 0.9999,  &
0.3228, 0.3395, 0.3564, 0.3737, 0.4450, 0.6259, 0.8837, 0.9970, 0.9999 /
data ((cldnuctab( 5,iw,iconc,6, 3),iw=1,9),iconc=1,7)/    &
0.6768, 0.8744, 0.9731, 0.9970, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4450, 0.7544, 0.9291, 0.9877, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3564, 0.4815, 0.8207, 0.9633, 0.9955, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.3564, 0.3737, 0.5181, 0.8647, 0.9826, 0.9985, 0.9999, 0.9999, 0.9999,  &
0.3564, 0.3737, 0.4090, 0.5727, 0.9007, 0.9925, 0.9995, 0.9999, 0.9999,  &
0.3564, 0.3564, 0.3737, 0.4269, 0.6259, 0.9227, 0.9966, 0.9999, 0.9999,  &
0.3395, 0.3395, 0.3564, 0.3737, 0.4269, 0.6084, 0.8647, 0.9961, 0.9999 /
data ((cldnuctab( 5,iw,iconc,7, 3),iw=1,9),iconc=1,7)/    &
0.6432, 0.8545, 0.9633, 0.9961, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4090, 0.7246, 0.9085, 0.9845, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.3737, 0.4450, 0.7824, 0.9508, 0.9934, 0.9995, 0.9999, 0.9999, 0.9999,  &
0.3737, 0.3912, 0.4815, 0.8325, 0.9731, 0.9977, 0.9998, 0.9999, 0.9999,  &
0.3737, 0.3737, 0.4090, 0.5181, 0.8744, 0.9877, 0.9992, 0.9999, 0.9999,  &
0.3564, 0.3737, 0.3912, 0.4269, 0.5727, 0.8924, 0.9942, 0.9997, 0.9999,  &
0.3564, 0.3564, 0.3737, 0.3912, 0.4269, 0.5727, 0.8325, 0.9934, 0.9999 /
data ((cldnuctab( 5,iw,iconc,1, 4),iw=1,9),iconc=1,7)/    &
0.9459, 0.9942, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7544, 0.9731, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5181, 0.8084, 0.9877, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4632, 0.5546, 0.8437, 0.9942, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4450, 0.4815, 0.5727, 0.7957, 0.9877, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4090, 0.4269, 0.4632, 0.5727, 0.5727, 0.7686, 0.9553, 0.9999, 0.9999,  &
0.3912, 0.3912, 0.4090, 0.4632, 0.5906, 0.3395, 0.3395, 0.4269, 0.6601 /
data ((cldnuctab( 5,iw,iconc,2, 4),iw=1,9),iconc=1,7)/    &
0.9085, 0.9862, 0.9989, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6931, 0.9508, 0.9955, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4998, 0.7397, 0.9758, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4815, 0.5364, 0.7957, 0.9877, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4632, 0.4815, 0.5546, 0.8084, 0.9915, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.4269, 0.4450, 0.4815, 0.5727, 0.6931, 0.9351, 0.9995, 0.9999, 0.9999,  &
0.4090, 0.4269, 0.4450, 0.4815, 0.5727, 0.4998, 0.5546, 0.7397, 0.9731 /
data ((cldnuctab( 5,iw,iconc,3, 4),iw=1,9),iconc=1,7)/    &
0.8647, 0.9758, 0.9974, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6084, 0.9158, 0.9904, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4998, 0.6768, 0.9553, 0.9966, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4815, 0.5181, 0.7246, 0.9758, 0.9989, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4632, 0.4815, 0.5364, 0.7686, 0.9845, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.4450, 0.4632, 0.4815, 0.5546, 0.7397, 0.9668, 0.9998, 0.9999, 0.9999,  &
0.4269, 0.4450, 0.4632, 0.4815, 0.5546, 0.6084, 0.7397, 0.9459, 0.9999 /
data ((cldnuctab( 5,iw,iconc,4, 4),iw=1,9),iconc=1,7)/    &
0.8207, 0.9595, 0.9949, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5546, 0.8837, 0.9806, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4998, 0.6084, 0.9227, 0.9925, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4998, 0.5181, 0.6601, 0.9508, 0.9974, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.4815, 0.4998, 0.5364, 0.7091, 0.9701, 0.9992, 0.9999, 0.9999, 0.9999,  &
0.4632, 0.4815, 0.4998, 0.5546, 0.7246, 0.9633, 0.9996, 0.9999, 0.9999,  &
0.4450, 0.4632, 0.4632, 0.4998, 0.5546, 0.6601, 0.8084, 0.9891, 0.9999 /
data ((cldnuctab( 5,iw,iconc,5, 4),iw=1,9),iconc=1,7)/    &
0.7824, 0.9407, 0.9915, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4998, 0.8437, 0.9701, 0.9970, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4998, 0.5364, 0.8837, 0.9862, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4998, 0.5181, 0.5906, 0.9227, 0.9942, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.4815, 0.4998, 0.5364, 0.6601, 0.9459, 0.9980, 0.9999, 0.9999, 0.9999,  &
0.4815, 0.4815, 0.4998, 0.5364, 0.6931, 0.9407, 0.9990, 0.9999, 0.9999,  &
0.4632, 0.4815, 0.4815, 0.4998, 0.5546, 0.6768, 0.8325, 0.9934, 0.9999 /
data ((cldnuctab( 5,iw,iconc,6, 4),iw=1,9),iconc=1,7)/    &
0.7397, 0.9227, 0.9877, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5181, 0.8084, 0.9595, 0.9949, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5181, 0.5364, 0.8545, 0.9783, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5181, 0.5364, 0.5727, 0.8837, 0.9904, 0.9995, 0.9999, 0.9999, 0.9999,  &
0.4998, 0.5181, 0.5364, 0.6259, 0.9085, 0.9961, 0.9998, 0.9999, 0.9999,  &
0.4998, 0.4998, 0.5181, 0.5546, 0.6601, 0.9158, 0.9980, 0.9999, 0.9999,  &
0.4815, 0.4815, 0.4998, 0.5181, 0.5546, 0.6601, 0.8325, 0.9915, 0.9999 /
data ((cldnuctab( 5,iw,iconc,7, 4),iw=1,9),iconc=1,7)/    &
0.7246, 0.9085, 0.9845, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5364, 0.7824, 0.9459, 0.9934, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5181, 0.5546, 0.8207, 0.9731, 0.9974, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.5181, 0.5364, 0.5727, 0.8545, 0.9862, 0.9992, 0.9999, 0.9999, 0.9999,  &
0.5181, 0.5181, 0.5546, 0.6084, 0.8744, 0.9934, 0.9997, 0.9999, 0.9999,  &
0.4998, 0.5181, 0.5181, 0.5546, 0.6432, 0.8837, 0.9961, 0.9999, 0.9999,  &
0.4998, 0.4998, 0.5181, 0.5181, 0.5546, 0.6432, 0.8084, 0.9891, 0.9999 /
data ((cldnuctab( 5,iw,iconc,1, 5),iw=1,9),iconc=1,7)/    &
0.9595, 0.9966, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7544, 0.9783, 0.9989, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5906, 0.8084, 0.9904, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5546, 0.6084, 0.8325, 0.9942, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4998, 0.5364, 0.6259, 0.7824, 0.9783, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4998, 0.5181, 0.5364, 0.6259, 0.6259, 0.7246, 0.9291, 0.9998, 0.9999,  &
0.4450, 0.4632, 0.4815, 0.5181, 0.6259, 0.5546, 0.3737, 0.4632, 0.6601 /
data ((cldnuctab( 5,iw,iconc,2, 5),iw=1,9),iconc=1,7)/    &
0.9291, 0.9915, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6931, 0.9595, 0.9970, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5906, 0.7397, 0.9783, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5546, 0.6084, 0.7957, 0.9891, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5364, 0.5546, 0.6084, 0.7957, 0.9891, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5181, 0.5181, 0.5364, 0.6084, 0.7246, 0.9158, 0.9990, 0.9999, 0.9999,  &
0.4998, 0.4998, 0.5181, 0.5546, 0.6259, 0.6601, 0.5546, 0.7246, 0.9407 /
data ((cldnuctab( 5,iw,iconc,3, 5),iw=1,9),iconc=1,7)/    &
0.8924, 0.9845, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6259, 0.9351, 0.9934, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5906, 0.6768, 0.9595, 0.9977, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5727, 0.6084, 0.7397, 0.9758, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5546, 0.5727, 0.6084, 0.7824, 0.9826, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.5364, 0.5364, 0.5546, 0.6084, 0.7544, 0.9553, 0.9997, 0.9999, 0.9999,  &
0.5181, 0.5181, 0.5364, 0.5546, 0.6084, 0.7091, 0.7246, 0.9158, 0.9995 /
data ((cldnuctab( 5,iw,iconc,4, 5),iw=1,9),iconc=1,7)/    &
0.8545, 0.9731, 0.9970, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5906, 0.9007, 0.9877, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5906, 0.6259, 0.9351, 0.9955, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5727, 0.6084, 0.6768, 0.9553, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5727, 0.5727, 0.6084, 0.7397, 0.9668, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.5546, 0.5546, 0.5727, 0.6084, 0.7397, 0.9508, 0.9996, 0.9999, 0.9999,  &
0.5364, 0.5364, 0.5546, 0.5727, 0.6084, 0.7246, 0.7824, 0.9758, 0.9999 /
data ((cldnuctab( 5,iw,iconc,5, 5),iw=1,9),iconc=1,7)/    &
0.8084, 0.9595, 0.9949, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5906, 0.8647, 0.9806, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5906, 0.6259, 0.9007, 0.9915, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5906, 0.6084, 0.6601, 0.9227, 0.9966, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.5727, 0.5906, 0.6084, 0.6931, 0.9407, 0.9985, 0.9999, 0.9999, 0.9999,  &
0.5727, 0.5727, 0.5906, 0.6084, 0.7091, 0.9291, 0.9992, 0.9999, 0.9999,  &
0.5546, 0.5546, 0.5727, 0.5906, 0.6084, 0.7091, 0.8207, 0.9891, 0.9999 /
data ((cldnuctab( 5,iw,iconc,6, 5),iw=1,9),iconc=1,7)/    &
0.7824, 0.9459, 0.9925, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6084, 0.8325, 0.9701, 0.9974, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6084, 0.6259, 0.8647, 0.9862, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5906, 0.6084, 0.6601, 0.8924, 0.9934, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.5906, 0.5906, 0.6259, 0.6768, 0.9085, 0.9970, 0.9999, 0.9999, 0.9999,  &
0.5727, 0.5906, 0.5906, 0.6259, 0.6931, 0.9007, 0.9983, 0.9999, 0.9999,  &
0.5727, 0.5727, 0.5727, 0.5906, 0.6259, 0.6931, 0.8325, 0.9877, 0.9999 /
data ((cldnuctab( 5,iw,iconc,7, 5),iw=1,9),iconc=1,7)/    &
0.7544, 0.9351, 0.9904, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6259, 0.7957, 0.9633, 0.9961, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6084, 0.6259, 0.8325, 0.9806, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6084, 0.6259, 0.6601, 0.8647, 0.9891, 0.9995, 0.9999, 0.9999, 0.9999,  &
0.5906, 0.6084, 0.6259, 0.6768, 0.8837, 0.9949, 0.9998, 0.9999, 0.9999,  &
0.5906, 0.5906, 0.6084, 0.6259, 0.6931, 0.8744, 0.9966, 0.9999, 0.9999,  &
0.5906, 0.5906, 0.5906, 0.6084, 0.6259, 0.6931, 0.8207, 0.9826, 0.9999 /
data ((cldnuctab( 5,iw,iconc,1, 6),iw=1,9),iconc=1,7)/    &
0.9633, 0.9974, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7544, 0.9826, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6432, 0.8084, 0.9915, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6084, 0.6601, 0.8325, 0.9942, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5727, 0.5906, 0.6601, 0.7824, 0.9668, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5546, 0.5546, 0.5906, 0.6601, 0.7246, 0.7246, 0.9085, 0.9994, 0.9999,  &
0.4998, 0.4998, 0.5181, 0.5546, 0.6432, 0.8325, 0.3912, 0.4815, 0.6768 /
data ((cldnuctab( 5,iw,iconc,2, 6),iw=1,9),iconc=1,7)/    &
0.9407, 0.9942, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6931, 0.9668, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6432, 0.7544, 0.9806, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6084, 0.6601, 0.8084, 0.9891, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5906, 0.6084, 0.6601, 0.7957, 0.9845, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5727, 0.5906, 0.6084, 0.6432, 0.7686, 0.8837, 0.9980, 0.9999, 0.9999,  &
0.5546, 0.5546, 0.5727, 0.5906, 0.6601, 0.7686, 0.5727, 0.7091, 0.9227 /
data ((cldnuctab( 5,iw,iconc,3, 6),iw=1,9),iconc=1,7)/    &
0.9085, 0.9877, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6601, 0.9407, 0.9955, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6432, 0.6931, 0.9633, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6259, 0.6601, 0.7544, 0.9783, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6084, 0.6259, 0.6601, 0.7824, 0.9806, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.5906, 0.5906, 0.6084, 0.6601, 0.7686, 0.9351, 0.9997, 0.9999, 0.9999,  &
0.5727, 0.5906, 0.5906, 0.6084, 0.6601, 0.7686, 0.7091, 0.8924, 0.9985 /
data ((cldnuctab( 5,iw,iconc,4, 6),iw=1,9),iconc=1,7)/    &
0.8647, 0.9783, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6601, 0.9085, 0.9915, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6432, 0.6768, 0.9351, 0.9966, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6432, 0.6601, 0.7246, 0.9553, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6259, 0.6432, 0.6601, 0.7544, 0.9633, 0.9995, 0.9999, 0.9999, 0.9999,  &
0.6084, 0.6259, 0.6259, 0.6601, 0.7544, 0.9407, 0.9996, 0.9999, 0.9999,  &
0.6084, 0.6084, 0.6084, 0.6259, 0.6601, 0.7544, 0.7824, 0.9595, 0.9999 /
data ((cldnuctab( 5,iw,iconc,5, 6),iw=1,9),iconc=1,7)/    &
0.8325, 0.9701, 0.9966, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6601, 0.8744, 0.9845, 0.9989, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6601, 0.6768, 0.9085, 0.9934, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6432, 0.6601, 0.7091, 0.9291, 0.9974, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6432, 0.6432, 0.6601, 0.7397, 0.9407, 0.9989, 0.9999, 0.9999, 0.9999,  &
0.6259, 0.6259, 0.6432, 0.6601, 0.7397, 0.9227, 0.9992, 0.9999, 0.9999,  &
0.6259, 0.6259, 0.6259, 0.6432, 0.6601, 0.7397, 0.8207, 0.9806, 0.9999 /
data ((cldnuctab( 5,iw,iconc,6, 6),iw=1,9),iconc=1,7)/    &
0.7957, 0.9595, 0.9949, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6768, 0.8437, 0.9783, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6601, 0.6768, 0.8744, 0.9891, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6601, 0.6768, 0.7091, 0.8924, 0.9949, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.6432, 0.6601, 0.6768, 0.7246, 0.9085, 0.9974, 0.9999, 0.9999, 0.9999,  &
0.6432, 0.6432, 0.6601, 0.6768, 0.7246, 0.9007, 0.9983, 0.9999, 0.9999,  &
0.6259, 0.6432, 0.6432, 0.6432, 0.6768, 0.7246, 0.8437, 0.9758, 0.9999 /
data ((cldnuctab( 5,iw,iconc,7, 6),iw=1,9),iconc=1,7)/    &
0.7686, 0.9459, 0.9934, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6768, 0.8084, 0.9701, 0.9974, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6768, 0.6931, 0.8437, 0.9845, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6601, 0.6768, 0.7091, 0.8647, 0.9915, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.6601, 0.6601, 0.6768, 0.7246, 0.8837, 0.9955, 0.9999, 0.9999, 0.9999,  &
0.6432, 0.6601, 0.6601, 0.6768, 0.7246, 0.8744, 0.9966, 0.9999, 0.9999,  &
0.6432, 0.6432, 0.6601, 0.6601, 0.6768, 0.7246, 0.8545, 0.9758, 0.9999 /
data ((cldnuctab( 5,iw,iconc,1, 7),iw=1,9),iconc=1,7)/    &
0.9668, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7686, 0.9826, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6768, 0.8207, 0.9915, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6432, 0.6931, 0.8325, 0.9925, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6084, 0.6259, 0.6768, 0.7957, 0.9633, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5906, 0.6084, 0.6259, 0.6931, 0.7824, 0.7091, 0.8924, 0.9985, 0.9999,  &
0.5181, 0.5364, 0.5364, 0.5727, 0.6432, 0.8207, 0.4090, 0.4815, 0.6601 /
data ((cldnuctab( 5,iw,iconc,2, 7),iw=1,9),iconc=1,7)/    &
0.9459, 0.9955, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7091, 0.9668, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6768, 0.7686, 0.9826, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6601, 0.6931, 0.8084, 0.9891, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6259, 0.6432, 0.6931, 0.8084, 0.9826, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6259, 0.6259, 0.6432, 0.6931, 0.7957, 0.8647, 0.9961, 0.9999, 0.9999,  &
0.5906, 0.5906, 0.6084, 0.6259, 0.6931, 0.8437, 0.5727, 0.7091, 0.9158 /
data ((cldnuctab( 5,iw,iconc,3, 7),iw=1,9),iconc=1,7)/    &
0.9158, 0.9904, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6931, 0.9459, 0.9966, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6931, 0.7246, 0.9633, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6768, 0.6931, 0.7686, 0.9758, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6601, 0.6601, 0.6931, 0.7957, 0.9783, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.6432, 0.6432, 0.6601, 0.6931, 0.7957, 0.9291, 0.9995, 0.9999, 0.9999,  &
0.6259, 0.6259, 0.6259, 0.6601, 0.6931, 0.8084, 0.7091, 0.8744, 0.9966 /
data ((cldnuctab( 5,iw,iconc,4, 7),iw=1,9),iconc=1,7)/    &
0.8744, 0.9826, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6931, 0.9158, 0.9925, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6931, 0.7246, 0.9407, 0.9974, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6768, 0.6931, 0.7544, 0.9553, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6601, 0.6768, 0.6931, 0.7686, 0.9633, 0.9996, 0.9999, 0.9999, 0.9999,  &
0.6601, 0.6601, 0.6768, 0.6931, 0.7824, 0.9291, 0.9995, 0.9999, 0.9999,  &
0.6432, 0.6432, 0.6601, 0.6601, 0.6931, 0.7824, 0.7957, 0.9553, 0.9999 /
data ((cldnuctab( 5,iw,iconc,5, 7),iw=1,9),iconc=1,7)/    &
0.8437, 0.9758, 0.9977, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7091, 0.8837, 0.9877, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6931, 0.7246, 0.9085, 0.9949, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6931, 0.7091, 0.7397, 0.9291, 0.9977, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6768, 0.6931, 0.7091, 0.7544, 0.9351, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.6601, 0.6768, 0.6931, 0.7091, 0.7686, 0.9158, 0.9992, 0.9999, 0.9999,  &
0.6601, 0.6601, 0.6768, 0.6768, 0.7091, 0.7686, 0.8437, 0.9668, 0.9999 /
data ((cldnuctab( 5,iw,iconc,6, 7),iw=1,9),iconc=1,7)/    &
0.8084, 0.9668, 0.9966, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7091, 0.8545, 0.9826, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7091, 0.7246, 0.8744, 0.9915, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6931, 0.7091, 0.7397, 0.8924, 0.9955, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.6931, 0.6931, 0.7091, 0.7544, 0.9085, 0.9977, 0.9999, 0.9999, 0.9999,  &
0.6768, 0.6931, 0.6931, 0.7091, 0.7544, 0.9007, 0.9980, 0.9999, 0.9999,  &
0.6768, 0.6768, 0.6768, 0.6931, 0.7091, 0.7544, 0.8647, 0.9701, 0.9999 /
data ((cldnuctab( 5,iw,iconc,7, 7),iw=1,9),iconc=1,7)/    &
0.7824, 0.9553, 0.9949, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7246, 0.8207, 0.9758, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7091, 0.7397, 0.8437, 0.9862, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7091, 0.7246, 0.7544, 0.8647, 0.9925, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.7091, 0.7091, 0.7246, 0.7544, 0.8837, 0.9961, 0.9999, 0.9999, 0.9999,  &
0.6931, 0.6931, 0.7091, 0.7246, 0.7544, 0.8837, 0.9966, 0.9999, 0.9999,  &
0.6931, 0.6931, 0.6931, 0.7091, 0.7246, 0.7544, 0.8647, 0.9633, 0.9999 /
data ((cldnuctab( 6,iw,iconc,1, 1),iw=1,9),iconc=1,7)/    &
0.9579, 0.9968, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7631, 0.9797, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5836, 0.8160, 0.9910, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5475, 0.6190, 0.8394, 0.9953, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5109, 0.5475, 0.6190, 0.7770, 0.9818, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.4926, 0.5109, 0.5475, 0.6190, 0.6190, 0.7338, 0.9328, 0.9999, 0.9999,  &
0.4560, 0.4560, 0.4743, 0.5292, 0.6190, 0.5292, 0.3843, 0.4743, 0.6868 /
data ((cldnuctab( 6,iw,iconc,2, 1),iw=1,9),iconc=1,7)/    &
0.9266, 0.9921, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6868, 0.9618, 0.9972, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5836, 0.7487, 0.9797, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5475, 0.6014, 0.8035, 0.9899, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5292, 0.5475, 0.6190, 0.8035, 0.9910, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5109, 0.5292, 0.5475, 0.6190, 0.7186, 0.9266, 0.9993, 0.9999, 0.9999,  &
0.4926, 0.5109, 0.5109, 0.5475, 0.6190, 0.6364, 0.5656, 0.7338, 0.9618 /
data ((cldnuctab( 6,iw,iconc,3, 1),iw=1,9),iconc=1,7)/    &
0.8891, 0.9838, 0.9986, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6190, 0.9328, 0.9939, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5836, 0.6703, 0.9618, 0.9979, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5656, 0.6014, 0.7487, 0.9773, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5475, 0.5656, 0.6014, 0.7770, 0.9855, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.5292, 0.5475, 0.5656, 0.6190, 0.7487, 0.9536, 0.9998, 0.9999, 0.9999,  &
0.5109, 0.5292, 0.5292, 0.5656, 0.6014, 0.7029, 0.7338, 0.9266, 0.9996 /
data ((cldnuctab( 6,iw,iconc,4, 1),iw=1,9),iconc=1,7)/    &
0.8503, 0.9719, 0.9972, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5836, 0.8975, 0.9886, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5836, 0.6190, 0.9328, 0.9953, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5656, 0.6014, 0.6703, 0.9579, 0.9984, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5656, 0.5836, 0.6014, 0.7338, 0.9688, 0.9995, 0.9999, 0.9999, 0.9999,  &
0.5475, 0.5656, 0.5656, 0.6190, 0.7338, 0.9536, 0.9997, 0.9999, 0.9999,  &
0.5292, 0.5475, 0.5475, 0.5656, 0.6190, 0.7186, 0.8035, 0.9818, 0.9999 /
data ((cldnuctab( 6,iw,iconc,5, 1),iw=1,9),iconc=1,7)/    &
0.8160, 0.9579, 0.9953, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6014, 0.8608, 0.9797, 0.9984, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5836, 0.6190, 0.8975, 0.9910, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5836, 0.6014, 0.6535, 0.9266, 0.9964, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.5656, 0.5836, 0.6190, 0.7028, 0.9439, 0.9988, 0.9999, 0.9999, 0.9999,  &
0.5656, 0.5656, 0.5836, 0.6190, 0.7186, 0.9328, 0.9993, 0.9999, 0.9999,  &
0.5475, 0.5656, 0.5656, 0.5836, 0.6190, 0.7186, 0.8160, 0.9921, 0.9999 /
data ((cldnuctab( 6,iw,iconc,6, 1),iw=1,9),iconc=1,7)/    &
0.7770, 0.9439, 0.9930, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6014, 0.8279, 0.9719, 0.9972, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6014, 0.6190, 0.8707, 0.9855, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.5836, 0.6014, 0.6535, 0.8975, 0.9939, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.5836, 0.6014, 0.6190, 0.6703, 0.9130, 0.9972, 0.9999, 0.9999, 0.9999,  &
0.5836, 0.5836, 0.6014, 0.6190, 0.7028, 0.9055, 0.9986, 0.9999, 0.9999,  &
0.5656, 0.5656, 0.5836, 0.5836, 0.6190, 0.7028, 0.8279, 0.9910, 0.9999 /
data ((cldnuctab( 6,iw,iconc,7, 1),iw=1,9),iconc=1,7)/    &
0.7487, 0.9328, 0.9899, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6190, 0.8035, 0.9618, 0.9964, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6014, 0.6364, 0.8394, 0.9797, 0.9986, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6014, 0.6190, 0.6535, 0.8608, 0.9899, 0.9996, 0.9999, 0.9999, 0.9999,  &
0.6014, 0.6014, 0.6190, 0.6703, 0.8801, 0.9953, 0.9998, 0.9999, 0.9999,  &
0.5836, 0.6014, 0.6014, 0.6190, 0.6868, 0.8801, 0.9972, 0.9999, 0.9999,  &
0.5836, 0.5836, 0.5836, 0.6014, 0.6190, 0.6868, 0.8279, 0.9871, 0.9999 /
data ((cldnuctab( 6,iw,iconc,1, 2),iw=1,9),iconc=1,7)/    &
0.9719, 0.9984, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7631, 0.9838, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7029, 0.8279, 0.9921, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6703, 0.7028, 0.8394, 0.9921, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6535, 0.6703, 0.7029, 0.8160, 0.9536, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.6190, 0.6364, 0.6535, 0.7186, 0.8160, 0.7186, 0.8801, 0.9979, 0.9999,  &
0.5475, 0.5475, 0.5656, 0.5836, 0.6535, 0.8160, 0.4560, 0.5109, 0.6868 /
data ((cldnuctab( 6,iw,iconc,2, 2),iw=1,9),iconc=1,7)/    &
0.9489, 0.9964, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7338, 0.9688, 0.9988, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7028, 0.7770, 0.9818, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6868, 0.7186, 0.8160, 0.9886, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6535, 0.6703, 0.7186, 0.8160, 0.9797, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6535, 0.6535, 0.6703, 0.7029, 0.8160, 0.8707, 0.9953, 0.9999, 0.9999,  &
0.6190, 0.6190, 0.6190, 0.6535, 0.7028, 0.8394, 0.6014, 0.7186, 0.9130 /
data ((cldnuctab( 6,iw,iconc,3, 2),iw=1,9),iconc=1,7)/    &
0.9200, 0.9921, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7186, 0.9489, 0.9972, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7186, 0.7487, 0.9655, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7029, 0.7186, 0.7905, 0.9773, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6868, 0.6868, 0.7186, 0.8035, 0.9748, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6703, 0.6703, 0.6868, 0.7029, 0.8035, 0.9200, 0.9994, 0.9999, 0.9999,  &
0.6535, 0.6535, 0.6535, 0.6703, 0.7186, 0.8279, 0.7338, 0.8707, 0.9959 /
data ((cldnuctab( 6,iw,iconc,4, 2),iw=1,9),iconc=1,7)/    &
0.8891, 0.9855, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7338, 0.9200, 0.9939, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7186, 0.7487, 0.9386, 0.9976, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7029, 0.7186, 0.7631, 0.9579, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7029, 0.7029, 0.7186, 0.7905, 0.9618, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.6868, 0.6868, 0.7029, 0.7186, 0.7905, 0.9266, 0.9995, 0.9999, 0.9999,  &
0.6703, 0.6703, 0.6868, 0.6868, 0.7186, 0.7905, 0.8160, 0.9536, 0.9999 /
data ((cldnuctab( 6,iw,iconc,5, 2),iw=1,9),iconc=1,7)/    &
0.8503, 0.9773, 0.9982, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7338, 0.8891, 0.9899, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7338, 0.7487, 0.9055, 0.9953, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7186, 0.7338, 0.7631, 0.9266, 0.9979, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7029, 0.7186, 0.7338, 0.7770, 0.9386, 0.9991, 0.9999, 0.9999, 0.9999,  &
0.7029, 0.7029, 0.7029, 0.7338, 0.7770, 0.9130, 0.9992, 0.9999, 0.9999,  &
0.6868, 0.6868, 0.7028, 0.7029, 0.7338, 0.7770, 0.8707, 0.9618, 0.9999 /
data ((cldnuctab( 6,iw,iconc,6, 2),iw=1,9),iconc=1,7)/    &
0.8160, 0.9688, 0.9972, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7338, 0.8503, 0.9838, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7338, 0.7487, 0.8801, 0.9921, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7186, 0.7338, 0.7631, 0.8975, 0.9964, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7186, 0.7186, 0.7338, 0.7770, 0.9130, 0.9982, 0.9999, 0.9999, 0.9999,  &
0.7028, 0.7186, 0.7186, 0.7338, 0.7770, 0.8975, 0.9979, 0.9999, 0.9999,  &
0.7029, 0.7029, 0.7029, 0.7186, 0.7338, 0.7770, 0.8801, 0.9655, 0.9999 /
data ((cldnuctab( 6,iw,iconc,7, 2),iw=1,9),iconc=1,7)/    &
0.7905, 0.9618, 0.9959, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7487, 0.8279, 0.9797, 0.9984, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7487, 0.7631, 0.8503, 0.9886, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7338, 0.7487, 0.7770, 0.8707, 0.9939, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.7338, 0.7338, 0.7487, 0.7770, 0.8891, 0.9964, 0.9999, 0.9999, 0.9999,  &
0.7186, 0.7186, 0.7338, 0.7487, 0.7770, 0.8891, 0.9964, 0.9999, 0.9999,  &
0.7186, 0.7186, 0.7186, 0.7338, 0.7487, 0.7770, 0.8801, 0.9655, 0.9999 /
data ((cldnuctab( 6,iw,iconc,1, 3),iw=1,9),iconc=1,7)/    &
0.9748, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8394, 0.9855, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8035, 0.8608, 0.9899, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7770, 0.7905, 0.8608, 0.9818, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7631, 0.7770, 0.8035, 0.8707, 0.9200, 0.9988, 0.9999, 0.9999, 0.9999,  &
0.7186, 0.7186, 0.7338, 0.7770, 0.8707, 0.7631, 0.8503, 0.9797, 0.9999,  &
0.6190, 0.6190, 0.6364, 0.6535, 0.6868, 0.7905, 0.9655, 0.5656, 0.7029 /
data ((cldnuctab( 6,iw,iconc,2, 3),iw=1,9),iconc=1,7)/    &
0.9579, 0.9982, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8394, 0.9719, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8160, 0.8503, 0.9797, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7905, 0.8160, 0.8608, 0.9818, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7770, 0.7905, 0.8035, 0.8608, 0.9536, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.7631, 0.7631, 0.7770, 0.8035, 0.8707, 0.8707, 0.9719, 0.9999, 0.9999,  &
0.7029, 0.7029, 0.7029, 0.7186, 0.7631, 0.8503, 0.8503, 0.7186, 0.8801 /
data ((cldnuctab( 6,iw,iconc,3, 3),iw=1,9),iconc=1,7)/    &
0.9328, 0.9964, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8394, 0.9489, 0.9986, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8160, 0.8503, 0.9618, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8035, 0.8160, 0.8503, 0.9719, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7905, 0.8035, 0.8160, 0.8608, 0.9579, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.7905, 0.7905, 0.7905, 0.8160, 0.8608, 0.9266, 0.9964, 0.9999, 0.9999,  &
0.7487, 0.7487, 0.7631, 0.7631, 0.7905, 0.8608, 0.9200, 0.8394, 0.9773 /
data ((cldnuctab( 6,iw,iconc,4, 3),iw=1,9),iconc=1,7)/    &
0.9055, 0.9930, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8394, 0.9200, 0.9968, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8279, 0.8503, 0.9386, 0.9986, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8160, 0.8279, 0.8503, 0.9489, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8035, 0.8160, 0.8279, 0.8503, 0.9489, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.8035, 0.8035, 0.8035, 0.8279, 0.8503, 0.9328, 0.9976, 0.9999, 0.9999,  &
0.7770, 0.7770, 0.7905, 0.7905, 0.8160, 0.8608, 0.9328, 0.9130, 0.9979 /
data ((cldnuctab( 6,iw,iconc,5, 3),iw=1,9),iconc=1,7)/    &
0.8707, 0.9886, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8394, 0.8891, 0.9946, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8394, 0.8503, 0.9055, 0.9972, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8279, 0.8394, 0.8503, 0.9266, 0.9986, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8160, 0.8160, 0.8279, 0.8503, 0.9328, 0.9992, 0.9999, 0.9999, 0.9999,  &
0.8160, 0.8160, 0.8160, 0.8279, 0.8503, 0.9328, 0.9972, 0.9999, 0.9999,  &
0.7905, 0.8035, 0.8035, 0.8035, 0.8279, 0.8608, 0.9328, 0.9386, 0.9997 /
data ((cldnuctab( 6,iw,iconc,6, 3),iw=1,9),iconc=1,7)/    &
0.8503, 0.9838, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8503, 0.8608, 0.9910, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8394, 0.8503, 0.8801, 0.9953, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8394, 0.8394, 0.8608, 0.9055, 0.9972, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8279, 0.8279, 0.8394, 0.8608, 0.9200, 0.9982, 0.9999, 0.9999, 0.9999,  &
0.8279, 0.8279, 0.8279, 0.8394, 0.8608, 0.9200, 0.9953, 0.9999, 0.9999,  &
0.8035, 0.8160, 0.8160, 0.8160, 0.8279, 0.8608, 0.9200, 0.9489, 0.9999 /
data ((cldnuctab( 6,iw,iconc,7, 3),iw=1,9),iconc=1,7)/    &
0.8608, 0.9773, 0.9984, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8503, 0.8707, 0.9871, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8503, 0.8608, 0.8801, 0.9921, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8394, 0.8503, 0.8608, 0.8975, 0.9953, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8394, 0.8394, 0.8503, 0.8608, 0.9130, 0.9964, 0.9999, 0.9999, 0.9999,  &
0.8279, 0.8394, 0.8394, 0.8394, 0.8608, 0.9130, 0.9910, 0.9999, 0.9999,  &
0.8160, 0.8160, 0.8279, 0.8279, 0.8394, 0.8608, 0.9130, 0.9618, 0.9998 /
data ((cldnuctab( 6,iw,iconc,1, 4),iw=1,9),iconc=1,7)/    &
0.9748, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9055, 0.9818, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8801, 0.9130, 0.9818, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8707, 0.8801, 0.9055, 0.9719, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8503, 0.8503, 0.8707, 0.9130, 0.9618, 0.9886, 0.9999, 0.9999, 0.9999,  &
0.7905, 0.7905, 0.8035, 0.8279, 0.8891, 0.9818, 0.8394, 0.9489, 0.9999,  &
0.6868, 0.6868, 0.6868, 0.7028, 0.7338, 0.7905, 0.9200, 0.6535, 0.7487 /
data ((cldnuctab( 6,iw,iconc,2, 4),iw=1,9),iconc=1,7)/    &
0.9579, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9055, 0.9688, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8975, 0.9130, 0.9748, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8801, 0.8891, 0.9130, 0.9719, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8707, 0.8707, 0.8891, 0.9130, 0.9655, 0.9992, 0.9999, 0.9999, 0.9999,  &
0.8394, 0.8394, 0.8503, 0.8707, 0.9055, 0.9838, 0.9439, 0.9991, 0.9999,  &
0.7631, 0.7631, 0.7631, 0.7770, 0.8035, 0.8608, 0.9618, 0.7487, 0.8608 /
data ((cldnuctab( 6,iw,iconc,3, 4),iw=1,9),iconc=1,7)/    &
0.9328, 0.9982, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9130, 0.9489, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8975, 0.9130, 0.9618, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8891, 0.8975, 0.9130, 0.9655, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8801, 0.8891, 0.8975, 0.9130, 0.9655, 0.9996, 0.9999, 0.9999, 0.9999,  &
0.8707, 0.8707, 0.8707, 0.8891, 0.9130, 0.9688, 0.9818, 0.9999, 0.9999,  &
0.8160, 0.8160, 0.8279, 0.8279, 0.8503, 0.8891, 0.9719, 0.8503, 0.9489 /
data ((cldnuctab( 6,iw,iconc,4, 4),iw=1,9),iconc=1,7)/    &
0.9200, 0.9964, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9130, 0.9328, 0.9982, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9055, 0.9130, 0.9439, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8975, 0.9055, 0.9130, 0.9536, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8891, 0.8975, 0.8975, 0.9130, 0.9579, 0.9992, 0.9999, 0.9999, 0.9999,  &
0.8801, 0.8801, 0.8891, 0.8975, 0.9130, 0.9579, 0.9899, 0.9999, 0.9999,  &
0.8503, 0.8503, 0.8503, 0.8608, 0.8707, 0.9055, 0.9655, 0.9130, 0.9855 /
data ((cldnuctab( 6,iw,iconc,5, 4),iw=1,9),iconc=1,7)/    &
0.9200, 0.9939, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9130, 0.9266, 0.9964, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9130, 0.9200, 0.9386, 0.9982, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9055, 0.9055, 0.9200, 0.9489, 0.9988, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8975, 0.9055, 0.9055, 0.9200, 0.9489, 0.9984, 0.9999, 0.9999, 0.9999,  &
0.8891, 0.8975, 0.8975, 0.9055, 0.9200, 0.9489, 0.9910, 0.9999, 0.9999,  &
0.8707, 0.8707, 0.8707, 0.8707, 0.8891, 0.9055, 0.9536, 0.9719, 0.9953 /
data ((cldnuctab( 6,iw,iconc,6, 4),iw=1,9),iconc=1,7)/    &
0.9266, 0.9899, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9200, 0.9328, 0.9946, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9130, 0.9200, 0.9386, 0.9964, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9130, 0.9130, 0.9200, 0.9439, 0.9972, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9055, 0.9055, 0.9130, 0.9200, 0.9439, 0.9968, 0.9999, 0.9999, 0.9999,  &
0.8975, 0.9055, 0.9055, 0.9130, 0.9200, 0.9439, 0.9910, 0.9999, 0.9999,  &
0.8801, 0.8801, 0.8801, 0.8891, 0.8891, 0.9130, 0.9489, 0.9871, 0.9982 /
data ((cldnuctab( 6,iw,iconc,7, 4),iw=1,9),iconc=1,7)/    &
0.9266, 0.9871, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9266, 0.9328, 0.9921, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9200, 0.9266, 0.9386, 0.9939, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9130, 0.9200, 0.9266, 0.9439, 0.9953, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9130, 0.9130, 0.9200, 0.9266, 0.9439, 0.9939, 0.9999, 0.9999, 0.9999,  &
0.9055, 0.9055, 0.9130, 0.9130, 0.9200, 0.9439, 0.9910, 0.9999, 0.9999,  &
0.8891, 0.8891, 0.8891, 0.8891, 0.8975, 0.9130, 0.9439, 0.9886, 0.9986 /
data ((cldnuctab( 6,iw,iconc,1, 5),iw=1,9),iconc=1,7)/    &
0.9719, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9328, 0.9797, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9200, 0.9328, 0.9797, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9055, 0.9200, 0.9328, 0.9773, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8801, 0.8891, 0.8975, 0.9328, 0.9871, 0.9797, 0.9999, 0.9999, 0.9999,  &
0.8160, 0.8160, 0.8279, 0.8503, 0.8891, 0.9719, 0.8503, 0.9386, 0.9997,  &
0.7186, 0.7186, 0.7186, 0.7338, 0.7487, 0.8035, 0.9055, 0.7186, 0.7770 /
data ((cldnuctab( 6,iw,iconc,2, 5),iw=1,9),iconc=1,7)/    &
0.9536, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9386, 0.9688, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9200, 0.9386, 0.9748, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9200, 0.9200, 0.9328, 0.9748, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9055, 0.9130, 0.9200, 0.9386, 0.9773, 0.9972, 0.9999, 0.9999, 0.9999,  &
0.8707, 0.8707, 0.8801, 0.8891, 0.9200, 0.9797, 0.9328, 0.9964, 0.9999,  &
0.7905, 0.7905, 0.8035, 0.8035, 0.8279, 0.8608, 0.9489, 0.8035, 0.8608 /
data ((cldnuctab( 6,iw,iconc,3, 5),iw=1,9),iconc=1,7)/    &
0.9489, 0.9986, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9386, 0.9579, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9328, 0.9386, 0.9655, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9266, 0.9266, 0.9386, 0.9688, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9200, 0.9200, 0.9266, 0.9386, 0.9719, 0.9988, 0.9999, 0.9999, 0.9999,  &
0.8975, 0.8975, 0.9055, 0.9130, 0.9328, 0.9773, 0.9748, 0.9999, 0.9999,  &
0.8503, 0.8503, 0.8503, 0.8608, 0.8707, 0.8975, 0.9655, 0.8975, 0.9386 /
data ((cldnuctab( 6,iw,iconc,4, 5),iw=1,9),iconc=1,7)/    &
0.9489, 0.9972, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9439, 0.9536, 0.9986, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9386, 0.9439, 0.9618, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9328, 0.9328, 0.9439, 0.9655, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9266, 0.9266, 0.9328, 0.9386, 0.9655, 0.9984, 0.9999, 0.9999, 0.9999,  &
0.9130, 0.9130, 0.9200, 0.9266, 0.9386, 0.9688, 0.9886, 0.9999, 0.9999,  &
0.8801, 0.8801, 0.8801, 0.8891, 0.8975, 0.9200, 0.9655, 0.9818, 0.9748 /
data ((cldnuctab( 6,iw,iconc,5, 5),iw=1,9),iconc=1,7)/    &
0.9489, 0.9953, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9439, 0.9536, 0.9972, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9386, 0.9439, 0.9579, 0.9982, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9386, 0.9386, 0.9439, 0.9618, 0.9984, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9328, 0.9328, 0.9386, 0.9439, 0.9618, 0.9972, 0.9999, 0.9999, 0.9999,  &
0.9266, 0.9266, 0.9266, 0.9328, 0.9439, 0.9618, 0.9939, 0.9999, 0.9999,  &
0.8975, 0.8975, 0.8975, 0.9055, 0.9055, 0.9266, 0.9618, 0.9939, 0.9899 /
data ((cldnuctab( 6,iw,iconc,6, 5),iw=1,9),iconc=1,7)/    &
0.9536, 0.9921, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9489, 0.9536, 0.9953, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9439, 0.9489, 0.9579, 0.9964, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9386, 0.9439, 0.9489, 0.9618, 0.9968, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9386, 0.9386, 0.9439, 0.9489, 0.9618, 0.9959, 0.9999, 0.9999, 0.9999,  &
0.9328, 0.9328, 0.9328, 0.9386, 0.9439, 0.9618, 0.9939, 0.9999, 0.9999,  &
0.9055, 0.9055, 0.9055, 0.9130, 0.9200, 0.9328, 0.9579, 0.9946, 0.9939 /
data ((cldnuctab( 6,iw,iconc,7, 5),iw=1,9),iconc=1,7)/    &
0.9536, 0.9899, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9489, 0.9579, 0.9930, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9489, 0.9489, 0.9579, 0.9939, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9439, 0.9439, 0.9489, 0.9618, 0.9946, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9439, 0.9439, 0.9439, 0.9489, 0.9618, 0.9930, 0.9999, 0.9999, 0.9999,  &
0.9386, 0.9386, 0.9386, 0.9386, 0.9489, 0.9618, 0.9921, 0.9999, 0.9999,  &
0.9130, 0.9130, 0.9130, 0.9200, 0.9200, 0.9328, 0.9579, 0.9930, 0.9959 /
data ((cldnuctab( 6,iw,iconc,1, 6),iw=1,9),iconc=1,7)/    &
0.9719, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9489, 0.9797, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9386, 0.9489, 0.9797, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9266, 0.9328, 0.9489, 0.9818, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8975, 0.9055, 0.9130, 0.9386, 0.9855, 0.9748, 0.9998, 0.9999, 0.9999,  &
0.8394, 0.8394, 0.8503, 0.8608, 0.8975, 0.9655, 0.8608, 0.9328, 0.9990,  &
0.7487, 0.7487, 0.7487, 0.7631, 0.7770, 0.8160, 0.8975, 0.9838, 0.8035 /
data ((cldnuctab( 6,iw,iconc,2, 6),iw=1,9),iconc=1,7)/    &
0.9655, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9536, 0.9719, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9439, 0.9536, 0.9773, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9386, 0.9439, 0.9489, 0.9773, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9266, 0.9266, 0.9328, 0.9489, 0.9818, 0.9939, 0.9999, 0.9999, 0.9999,  &
0.8891, 0.8891, 0.8975, 0.9055, 0.9328, 0.9797, 0.9386, 0.9921, 0.9999,  &
0.8160, 0.8160, 0.8160, 0.8279, 0.8394, 0.8707, 0.9439, 0.9995, 0.8707 /
data ((cldnuctab( 6,iw,iconc,3, 6),iw=1,9),iconc=1,7)/    &
0.9618, 0.9988, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9536, 0.9688, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9489, 0.9536, 0.9719, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9439, 0.9489, 0.9536, 0.9748, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9386, 0.9386, 0.9439, 0.9536, 0.9748, 0.9976, 0.9999, 0.9999, 0.9999,  &
0.9200, 0.9200, 0.9200, 0.9266, 0.9439, 0.9797, 0.9818, 0.9997, 0.9999,  &
0.8707, 0.8707, 0.8707, 0.8707, 0.8801, 0.9055, 0.9618, 0.9998, 0.9386 /
data ((cldnuctab( 6,iw,iconc,4, 6),iw=1,9),iconc=1,7)/    &
0.9618, 0.9976, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9579, 0.9655, 0.9986, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9536, 0.9579, 0.9688, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9489, 0.9489, 0.9579, 0.9719, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9439, 0.9489, 0.9489, 0.9579, 0.9719, 0.9976, 0.9999, 0.9999, 0.9999,  &
0.9328, 0.9328, 0.9328, 0.9386, 0.9536, 0.9748, 0.9946, 0.9999, 0.9999,  &
0.8975, 0.8975, 0.8975, 0.8975, 0.9055, 0.9266, 0.9655, 0.9995, 0.9719 /
data ((cldnuctab( 6,iw,iconc,5, 6),iw=1,9),iconc=1,7)/    &
0.9655, 0.9959, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9618, 0.9688, 0.9972, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9579, 0.9618, 0.9688, 0.9979, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9536, 0.9536, 0.9579, 0.9688, 0.9982, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9489, 0.9536, 0.9536, 0.9579, 0.9688, 0.9968, 0.9999, 0.9999, 0.9999,  &
0.9439, 0.9439, 0.9439, 0.9489, 0.9536, 0.9719, 0.9959, 0.9999, 0.9999,  &
0.9130, 0.9130, 0.9130, 0.9200, 0.9200, 0.9386, 0.9655, 0.9986, 0.9855 /
data ((cldnuctab( 6,iw,iconc,6, 6),iw=1,9),iconc=1,7)/    &
0.9655, 0.9939, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9618, 0.9688, 0.9953, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9579, 0.9618, 0.9688, 0.9964, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9579, 0.9579, 0.9618, 0.9688, 0.9964, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9536, 0.9536, 0.9579, 0.9618, 0.9688, 0.9953, 0.9999, 0.9999, 0.9999,  &
0.9489, 0.9489, 0.9489, 0.9536, 0.9579, 0.9719, 0.9946, 0.9999, 0.9999,  &
0.9200, 0.9200, 0.9266, 0.9266, 0.9328, 0.9386, 0.9618, 0.9968, 0.9921 /
data ((cldnuctab( 6,iw,iconc,7, 6),iw=1,9),iconc=1,7)/    &
0.9688, 0.9910, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9655, 0.9688, 0.9930, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9618, 0.9655, 0.9719, 0.9939, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9618, 0.9618, 0.9655, 0.9719, 0.9939, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9579, 0.9579, 0.9618, 0.9618, 0.9719, 0.9939, 0.9999, 0.9999, 0.9999,  &
0.9489, 0.9536, 0.9536, 0.9536, 0.9618, 0.9719, 0.9930, 0.9999, 0.9999,  &
0.9266, 0.9266, 0.9266, 0.9328, 0.9328, 0.9439, 0.9618, 0.9946, 0.9959 /
data ((cldnuctab( 6,iw,iconc,1, 7),iw=1,9),iconc=1,7)/    &
0.9748, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9618, 0.9818, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9536, 0.9579, 0.9818, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9386, 0.9489, 0.9579, 0.9838, 0.9988, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9130, 0.9130, 0.9266, 0.9439, 0.9855, 0.9719, 0.9996, 0.9999, 0.9999,  &
0.8503, 0.8503, 0.8608, 0.8707, 0.8975, 0.9618, 0.8975, 0.9328, 0.9972,  &
0.7770, 0.7770, 0.7770, 0.7770, 0.7905, 0.8279, 0.8975, 0.9748, 0.8160 /
data ((cldnuctab( 6,iw,iconc,2, 7),iw=1,9),iconc=1,7)/    &
0.9719, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9618, 0.9773, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9536, 0.9618, 0.9797, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9489, 0.9536, 0.9618, 0.9797, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9386, 0.9386, 0.9439, 0.9579, 0.9838, 0.9939, 0.9999, 0.9999, 0.9999,  &
0.8975, 0.9055, 0.9055, 0.9130, 0.9386, 0.9773, 0.9719, 0.9871, 0.9999,  &
0.8279, 0.8279, 0.8279, 0.8394, 0.8503, 0.8801, 0.9386, 0.9984, 0.8801 /
data ((cldnuctab( 6,iw,iconc,3, 7),iw=1,9),iconc=1,7)/    &
0.9719, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9655, 0.9748, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9579, 0.9655, 0.9773, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9579, 0.9579, 0.9618, 0.9773, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9489, 0.9489, 0.9536, 0.9618, 0.9797, 0.9972, 0.9999, 0.9999, 0.9999,  &
0.9266, 0.9266, 0.9328, 0.9386, 0.9536, 0.9797, 0.9946, 0.9992, 0.9999,  &
0.8801, 0.8801, 0.8801, 0.8801, 0.8891, 0.9130, 0.9579, 0.9995, 0.9386 /
data ((cldnuctab( 6,iw,iconc,4, 7),iw=1,9),iconc=1,7)/    &
0.9719, 0.9979, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9688, 0.9748, 0.9986, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9618, 0.9655, 0.9748, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9618, 0.9618, 0.9655, 0.9773, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9579, 0.9579, 0.9618, 0.9655, 0.9748, 0.9976, 0.9999, 0.9999, 0.9999,  &
0.9439, 0.9439, 0.9439, 0.9489, 0.9579, 0.9797, 0.9976, 0.9999, 0.9999,  &
0.9055, 0.9055, 0.9055, 0.9130, 0.9200, 0.9328, 0.9655, 0.9991, 0.9688 /
data ((cldnuctab( 6,iw,iconc,5, 7),iw=1,9),iconc=1,7)/    &
0.9719, 0.9964, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9688, 0.9748, 0.9972, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9655, 0.9688, 0.9748, 0.9979, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9655, 0.9655, 0.9688, 0.9773, 0.9979, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9618, 0.9618, 0.9655, 0.9688, 0.9773, 0.9968, 0.9999, 0.9999, 0.9999,  &
0.9536, 0.9536, 0.9536, 0.9579, 0.9618, 0.9773, 0.9968, 0.9999, 0.9999,  &
0.9200, 0.9266, 0.9266, 0.9266, 0.9328, 0.9439, 0.9655, 0.9979, 0.9838 /
data ((cldnuctab( 6,iw,iconc,6, 7),iw=1,9),iconc=1,7)/    &
0.9748, 0.9939, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9719, 0.9748, 0.9953, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9688, 0.9719, 0.9773, 0.9959, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9655, 0.9688, 0.9719, 0.9773, 0.9959, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9655, 0.9655, 0.9655, 0.9688, 0.9773, 0.9953, 0.9999, 0.9999, 0.9999,  &
0.9579, 0.9579, 0.9579, 0.9618, 0.9655, 0.9773, 0.9959, 0.9999, 0.9999,  &
0.9328, 0.9328, 0.9328, 0.9328, 0.9386, 0.9489, 0.9655, 0.9959, 0.9921 /
data ((cldnuctab( 6,iw,iconc,7, 7),iw=1,9),iconc=1,7)/    &
0.9748, 0.9910, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9719, 0.9773, 0.9930, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9719, 0.9719, 0.9773, 0.9939, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9688, 0.9688, 0.9719, 0.9773, 0.9939, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9688, 0.9688, 0.9688, 0.9719, 0.9773, 0.9939, 0.9999, 0.9999, 0.9999,  &
0.9618, 0.9618, 0.9618, 0.9618, 0.9688, 0.9773, 0.9939, 0.9999, 0.9999,  &
0.9386, 0.9386, 0.9386, 0.9386, 0.9439, 0.9489, 0.9655, 0.9939, 0.9979 /
data ((cldnuctab( 7,iw,iconc,1, 1),iw=1,9),iconc=1,7)/    &
0.9742, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7741, 0.9852, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7306, 0.8370, 0.9919, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6994, 0.7306, 0.8480, 0.9919, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6833, 0.6994, 0.7306, 0.8370, 0.9479, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.6499, 0.6667, 0.6833, 0.7306, 0.8686, 0.7306, 0.8781, 0.9963, 0.9999,  &
0.5797, 0.5797, 0.5976, 0.6153, 0.6667, 0.8133, 0.4887, 0.5617, 0.6994 /
data ((cldnuctab( 7,iw,iconc,2, 1),iw=1,9),iconc=1,7)/    &
0.9526, 0.9972, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7600, 0.9713, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7455, 0.8007, 0.9834, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7152, 0.7455, 0.8254, 0.9883, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6994, 0.7152, 0.7455, 0.8254, 0.9742, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6833, 0.6994, 0.6994, 0.7455, 0.8370, 0.8586, 0.9919, 0.9999, 0.9999,  &
0.6499, 0.6499, 0.6499, 0.6833, 0.7306, 0.8480, 0.6327, 0.7306, 0.9114 /
data ((cldnuctab( 7,iw,iconc,3, 1),iw=1,9),iconc=1,7)/    &
0.9252, 0.9937, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7600, 0.9479, 0.9979, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7455, 0.7741, 0.9681, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7306, 0.7455, 0.8007, 0.9768, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7152, 0.7306, 0.7455, 0.8133, 0.9713, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.6994, 0.7152, 0.7152, 0.7455, 0.8254, 0.9185, 0.9992, 0.9999, 0.9999,  &
0.6833, 0.6833, 0.6994, 0.6994, 0.7455, 0.8370, 0.7600, 0.8781, 0.9928 /
data ((cldnuctab( 7,iw,iconc,4, 1),iw=1,9),iconc=1,7)/    &
0.8957, 0.9883, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7600, 0.9185, 0.9951, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7455, 0.7741, 0.9428, 0.9981, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7455, 0.7600, 0.7877, 0.9570, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7306, 0.7306, 0.7600, 0.8007, 0.9570, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.7152, 0.7306, 0.7306, 0.7455, 0.8133, 0.9252, 0.9996, 0.9999, 0.9999,  &
0.6994, 0.7152, 0.7152, 0.7306, 0.7455, 0.8133, 0.8480, 0.9428, 0.9998 /
data ((cldnuctab( 7,iw,iconc,5, 1),iw=1,9),iconc=1,7)/    &
0.8586, 0.9814, 0.9986, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7600, 0.8872, 0.9908, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7600, 0.7741, 0.9114, 0.9963, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7455, 0.7600, 0.7877, 0.9252, 0.9984, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7455, 0.7455, 0.7600, 0.8007, 0.9374, 0.9993, 0.9999, 0.9999, 0.9999,  &
0.7306, 0.7306, 0.7455, 0.7600, 0.8007, 0.9185, 0.9991, 0.9999, 0.9999,  &
0.7306, 0.7306, 0.7306, 0.7455, 0.7600, 0.8007, 0.8872, 0.9647, 0.9999 /
data ((cldnuctab( 7,iw,iconc,6, 1),iw=1,9),iconc=1,7)/    &
0.8254, 0.9742, 0.9979, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7741, 0.8586, 0.9868, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7600, 0.7877, 0.8781, 0.9937, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7600, 0.7741, 0.7877, 0.8957, 0.9968, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7455, 0.7600, 0.7741, 0.8007, 0.9114, 0.9984, 0.9999, 0.9999, 0.9999,  &
0.7455, 0.7455, 0.7600, 0.7600, 0.8007, 0.9038, 0.9981, 0.9999, 0.9999,  &
0.7306, 0.7455, 0.7455, 0.7455, 0.7600, 0.8007, 0.8957, 0.9647, 0.9999 /
data ((cldnuctab( 7,iw,iconc,7, 1),iw=1,9),iconc=1,7)/    &
0.8007, 0.9681, 0.9972, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7741, 0.8254, 0.9814, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7741, 0.7877, 0.8480, 0.9896, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7741, 0.7741, 0.8007, 0.8686, 0.9945, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.7600, 0.7600, 0.7741, 0.8007, 0.8957, 0.9968, 0.9999, 0.9999, 0.9999,  &
0.7600, 0.7600, 0.7600, 0.7741, 0.8007, 0.8957, 0.9957, 0.9999, 0.9999,  &
0.7455, 0.7455, 0.7455, 0.7600, 0.7741, 0.8007, 0.8872, 0.9610, 0.9999 /
data ((cldnuctab( 7,iw,iconc,1, 2),iw=1,9),iconc=1,7)/    &
0.9768, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8586, 0.9852, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8370, 0.8686, 0.9883, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8133, 0.8254, 0.8781, 0.9792, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7877, 0.8007, 0.8254, 0.8781, 0.9185, 0.9979, 0.9999, 0.9999, 0.9999,  &
0.7455, 0.7455, 0.7600, 0.8007, 0.8781, 0.8133, 0.8480, 0.9742, 0.9999,  &
0.6499, 0.6499, 0.6499, 0.6667, 0.6994, 0.8007, 0.9526, 0.5976, 0.7306 /
data ((cldnuctab( 7,iw,iconc,2, 2),iw=1,9),iconc=1,7)/    &
0.9570, 0.9986, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8586, 0.9713, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8370, 0.8686, 0.9792, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8133, 0.8370, 0.8686, 0.9792, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8133, 0.8133, 0.8254, 0.8781, 0.9526, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.7877, 0.7877, 0.8007, 0.8254, 0.8872, 0.9114, 0.9647, 0.9999, 0.9999,  &
0.7152, 0.7306, 0.7306, 0.7455, 0.7741, 0.8586, 0.9814, 0.7306, 0.8781 /
data ((cldnuctab( 7,iw,iconc,3, 2),iw=1,9),iconc=1,7)/    &
0.9315, 0.9968, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8586, 0.9526, 0.9988, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8480, 0.8686, 0.9610, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8370, 0.8480, 0.8686, 0.9681, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8254, 0.8254, 0.8370, 0.8686, 0.9570, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.8133, 0.8133, 0.8254, 0.8370, 0.8781, 0.9374, 0.9945, 0.9999, 0.9999,  &
0.7741, 0.7741, 0.7741, 0.7877, 0.8133, 0.8686, 0.9792, 0.8480, 0.9742 /
data ((cldnuctab( 7,iw,iconc,4, 2),iw=1,9),iconc=1,7)/    &
0.9038, 0.9937, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8586, 0.9252, 0.9975, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8480, 0.8686, 0.9374, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8370, 0.8480, 0.8686, 0.9526, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8370, 0.8370, 0.8480, 0.8686, 0.9526, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.8254, 0.8254, 0.8370, 0.8480, 0.8686, 0.9428, 0.9972, 0.9999, 0.9999,  &
0.8007, 0.8007, 0.8007, 0.8133, 0.8370, 0.8686, 0.9526, 0.9038, 0.9968 /
data ((cldnuctab( 7,iw,iconc,5, 2),iw=1,9),iconc=1,7)/    &
0.8781, 0.9896, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8686, 0.8957, 0.9951, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8586, 0.8686, 0.9114, 0.9975, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8480, 0.8586, 0.8686, 0.9315, 0.9988, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8370, 0.8480, 0.8586, 0.8686, 0.9374, 0.9992, 0.9999, 0.9999, 0.9999,  &
0.8370, 0.8370, 0.8480, 0.8480, 0.8686, 0.9374, 0.9963, 0.9999, 0.9999,  &
0.8254, 0.8254, 0.8254, 0.8254, 0.8370, 0.8686, 0.9428, 0.9374, 0.9995 /
data ((cldnuctab( 7,iw,iconc,6, 2),iw=1,9),iconc=1,7)/    &
0.8781, 0.9852, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8686, 0.8872, 0.9919, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8586, 0.8781, 0.9038, 0.9957, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8586, 0.8586, 0.8781, 0.9185, 0.9975, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8480, 0.8586, 0.8586, 0.8781, 0.9252, 0.9979, 0.9999, 0.9999, 0.9999,  &
0.8480, 0.8480, 0.8480, 0.8586, 0.8781, 0.9252, 0.9928, 0.9999, 0.9999,  &
0.8370, 0.8370, 0.8370, 0.8370, 0.8480, 0.8781, 0.9315, 0.9526, 0.9998 /
data ((cldnuctab( 7,iw,iconc,7, 2),iw=1,9),iconc=1,7)/    &
0.8781, 0.9814, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8781, 0.8872, 0.9883, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8686, 0.8781, 0.9038, 0.9928, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8686, 0.8686, 0.8781, 0.9114, 0.9951, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8586, 0.8586, 0.8686, 0.8781, 0.9185, 0.9957, 0.9999, 0.9999, 0.9999,  &
0.8586, 0.8586, 0.8586, 0.8686, 0.8781, 0.9185, 0.9908, 0.9999, 0.9999,  &
0.8370, 0.8370, 0.8480, 0.8480, 0.8586, 0.8781, 0.9185, 0.9713, 0.9998 /
data ((cldnuctab( 7,iw,iconc,1, 3),iw=1,9),iconc=1,7)/    &
0.9742, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9185, 0.9814, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8957, 0.9185, 0.9814, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8872, 0.8957, 0.9185, 0.9742, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8686, 0.8686, 0.8872, 0.9252, 0.9742, 0.9868, 0.9999, 0.9999, 0.9999,  &
0.8007, 0.8007, 0.8133, 0.8370, 0.8872, 0.9792, 0.8370, 0.9479, 0.9999,  &
0.6994, 0.6994, 0.7152, 0.7152, 0.7455, 0.8007, 0.9114, 0.6833, 0.7600 /
data ((cldnuctab( 7,iw,iconc,2, 3),iw=1,9),iconc=1,7)/    &
0.9570, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9252, 0.9681, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9114, 0.9252, 0.9742, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8957, 0.9038, 0.9185, 0.9742, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8872, 0.8957, 0.9038, 0.9252, 0.9713, 0.9986, 0.9999, 0.9999, 0.9999,  &
0.8586, 0.8586, 0.8586, 0.8781, 0.9185, 0.9834, 0.9374, 0.9984, 0.9999,  &
0.7741, 0.7877, 0.7877, 0.7877, 0.8133, 0.8586, 0.9570, 0.7741, 0.8586 /
data ((cldnuctab( 7,iw,iconc,3, 3),iw=1,9),iconc=1,7)/    &
0.9315, 0.9984, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9252, 0.9479, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9114, 0.9252, 0.9647, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9038, 0.9114, 0.9252, 0.9681, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9038, 0.9038, 0.9114, 0.9252, 0.9681, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.8781, 0.8872, 0.8872, 0.8957, 0.9252, 0.9768, 0.9792, 0.9999, 0.9999,  &
0.8370, 0.8370, 0.8370, 0.8480, 0.8586, 0.8957, 0.9681, 0.8586, 0.9428 /
data ((cldnuctab( 7,iw,iconc,4, 3),iw=1,9),iconc=1,7)/    &
0.9315, 0.9968, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9252, 0.9428, 0.9984, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9185, 0.9252, 0.9526, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9114, 0.9185, 0.9252, 0.9610, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9114, 0.9114, 0.9185, 0.9252, 0.9610, 0.9991, 0.9999, 0.9999, 0.9999,  &
0.8957, 0.8957, 0.9038, 0.9114, 0.9252, 0.9647, 0.9868, 0.9999, 0.9999,  &
0.8686, 0.8686, 0.8686, 0.8686, 0.8872, 0.9114, 0.9647, 0.9374, 0.9834 /
data ((cldnuctab( 7,iw,iconc,5, 3),iw=1,9),iconc=1,7)/    &
0.9374, 0.9945, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9315, 0.9428, 0.9968, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9252, 0.9315, 0.9479, 0.9981, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9185, 0.9252, 0.9315, 0.9526, 0.9986, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9185, 0.9185, 0.9185, 0.9315, 0.9526, 0.9981, 0.9999, 0.9999, 0.9999,  &
0.9114, 0.9114, 0.9114, 0.9185, 0.9315, 0.9570, 0.9919, 0.9999, 0.9999,  &
0.8872, 0.8872, 0.8872, 0.8872, 0.8957, 0.9185, 0.9570, 0.9834, 0.9945 /
data ((cldnuctab( 7,iw,iconc,6, 3),iw=1,9),iconc=1,7)/    &
0.9374, 0.9919, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9315, 0.9428, 0.9945, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9315, 0.9315, 0.9479, 0.9963, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9252, 0.9252, 0.9315, 0.9526, 0.9972, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9252, 0.9252, 0.9252, 0.9315, 0.9526, 0.9963, 0.9999, 0.9999, 0.9999,  &
0.9185, 0.9185, 0.9185, 0.9252, 0.9315, 0.9526, 0.9919, 0.9999, 0.9999,  &
0.8957, 0.8957, 0.8957, 0.8957, 0.9038, 0.9185, 0.9526, 0.9908, 0.9968 /
data ((cldnuctab( 7,iw,iconc,7, 3),iw=1,9),iconc=1,7)/    &
0.9374, 0.9883, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9374, 0.9428, 0.9919, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9315, 0.9374, 0.9479, 0.9945, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9315, 0.9315, 0.9374, 0.9526, 0.9951, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9252, 0.9252, 0.9315, 0.9374, 0.9526, 0.9937, 0.9999, 0.9999, 0.9999,  &
0.9185, 0.9185, 0.9252, 0.9252, 0.9374, 0.9526, 0.9919, 0.9999, 0.9999,  &
0.9038, 0.9038, 0.9038, 0.9038, 0.9114, 0.9252, 0.9479, 0.9908, 0.9972 /
data ((cldnuctab( 7,iw,iconc,1, 4),iw=1,9),iconc=1,7)/    &
0.9742, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9570, 0.9814, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9479, 0.9570, 0.9814, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9374, 0.9428, 0.9570, 0.9834, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9114, 0.9114, 0.9252, 0.9428, 0.9852, 0.9713, 0.9997, 0.9999, 0.9999,  &
0.8480, 0.8480, 0.8586, 0.8686, 0.9038, 0.9610, 0.8872, 0.9374, 0.9979,  &
0.7741, 0.7741, 0.7741, 0.7741, 0.7877, 0.8254, 0.8957, 0.9768, 0.8254 /
data ((cldnuctab( 7,iw,iconc,2, 4),iw=1,9),iconc=1,7)/    &
0.9713, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9610, 0.9768, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9526, 0.9570, 0.9792, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9479, 0.9526, 0.9610, 0.9792, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9374, 0.9374, 0.9428, 0.9570, 0.9834, 0.9937, 0.9999, 0.9999, 0.9999,  &
0.8957, 0.9038, 0.9038, 0.9114, 0.9374, 0.9768, 0.9570, 0.9883, 0.9999,  &
0.8254, 0.8254, 0.8370, 0.8370, 0.8480, 0.8781, 0.9374, 0.9988, 0.8872 /
data ((cldnuctab( 7,iw,iconc,3, 4),iw=1,9),iconc=1,7)/    &
0.9713, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9647, 0.9742, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9570, 0.9647, 0.9768, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9526, 0.9570, 0.9610, 0.9768, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9479, 0.9479, 0.9526, 0.9610, 0.9792, 0.9972, 0.9999, 0.9999, 0.9999,  &
0.9252, 0.9252, 0.9315, 0.9374, 0.9526, 0.9792, 0.9928, 0.9994, 0.9999,  &
0.8781, 0.8781, 0.8781, 0.8781, 0.8872, 0.9114, 0.9610, 0.9996, 0.9374 /
data ((cldnuctab( 7,iw,iconc,4, 4),iw=1,9),iconc=1,7)/    &
0.9713, 0.9979, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9647, 0.9742, 0.9986, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9610, 0.9647, 0.9742, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9570, 0.9610, 0.9647, 0.9742, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9526, 0.9570, 0.9570, 0.9647, 0.9742, 0.9975, 0.9999, 0.9999, 0.9999,  &
0.9428, 0.9428, 0.9428, 0.9479, 0.9570, 0.9768, 0.9968, 0.9999, 0.9999,  &
0.9038, 0.9038, 0.9038, 0.9114, 0.9185, 0.9315, 0.9647, 0.9992, 0.9713 /
data ((cldnuctab( 7,iw,iconc,5, 4),iw=1,9),iconc=1,7)/    &
0.9713, 0.9963, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9681, 0.9742, 0.9975, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9647, 0.9681, 0.9742, 0.9979, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9610, 0.9647, 0.9681, 0.9742, 0.9979, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9610, 0.9610, 0.9610, 0.9647, 0.9742, 0.9968, 0.9999, 0.9999, 0.9999,  &
0.9479, 0.9479, 0.9526, 0.9526, 0.9610, 0.9768, 0.9968, 0.9999, 0.9999,  &
0.9185, 0.9185, 0.9252, 0.9252, 0.9315, 0.9428, 0.9647, 0.9981, 0.9852 /
data ((cldnuctab( 7,iw,iconc,6, 4),iw=1,9),iconc=1,7)/    &
0.9713, 0.9937, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9681, 0.9742, 0.9957, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9681, 0.9681, 0.9742, 0.9963, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9647, 0.9647, 0.9681, 0.9742, 0.9963, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9610, 0.9647, 0.9647, 0.9681, 0.9742, 0.9957, 0.9999, 0.9999, 0.9999,  &
0.9526, 0.9570, 0.9570, 0.9570, 0.9647, 0.9742, 0.9957, 0.9999, 0.9999,  &
0.9315, 0.9315, 0.9315, 0.9315, 0.9374, 0.9479, 0.9647, 0.9963, 0.9919 /
data ((cldnuctab( 7,iw,iconc,7, 4),iw=1,9),iconc=1,7)/    &
0.9742, 0.9919, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9713, 0.9742, 0.9928, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9681, 0.9713, 0.9768, 0.9937, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9681, 0.9681, 0.9713, 0.9768, 0.9937, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9647, 0.9647, 0.9681, 0.9713, 0.9768, 0.9937, 0.9999, 0.9999, 0.9999,  &
0.9570, 0.9570, 0.9610, 0.9610, 0.9647, 0.9742, 0.9937, 0.9999, 0.9999,  &
0.9374, 0.9374, 0.9374, 0.9374, 0.9428, 0.9479, 0.9647, 0.9937, 0.9972 /
data ((cldnuctab( 7,iw,iconc,1, 5),iw=1,9),iconc=1,7)/    &
0.9834, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9713, 0.9852, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9681, 0.9713, 0.9852, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9570, 0.9610, 0.9713, 0.9883, 0.9979, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9315, 0.9315, 0.9374, 0.9526, 0.9834, 0.9883, 0.9984, 0.9999, 0.9999,  &
0.8781, 0.8781, 0.8781, 0.8872, 0.9114, 0.9570, 0.9994, 0.9374, 0.9945,  &
0.8133, 0.8133, 0.8133, 0.8133, 0.8254, 0.8480, 0.9038, 0.9681, 0.8586 /
data ((cldnuctab( 7,iw,iconc,2, 5),iw=1,9),iconc=1,7)/    &
0.9814, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9768, 0.9834, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9713, 0.9742, 0.9852, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9647, 0.9681, 0.9742, 0.9852, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9526, 0.9526, 0.9570, 0.9681, 0.9868, 0.9972, 0.9999, 0.9999, 0.9999,  &
0.9185, 0.9185, 0.9185, 0.9252, 0.9428, 0.9742, 0.9998, 0.9852, 0.9999,  &
0.8586, 0.8586, 0.8586, 0.8586, 0.8686, 0.8872, 0.9374, 0.9945, 0.8957 /
data ((cldnuctab( 7,iw,iconc,3, 5),iw=1,9),iconc=1,7)/    &
0.9814, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9768, 0.9834, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9742, 0.9742, 0.9834, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9713, 0.9713, 0.9768, 0.9834, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9647, 0.9647, 0.9681, 0.9742, 0.9852, 0.9984, 0.9999, 0.9999, 0.9999,  &
0.9428, 0.9428, 0.9479, 0.9479, 0.9610, 0.9814, 0.9997, 0.9981, 0.9999,  &
0.8957, 0.8957, 0.8957, 0.8957, 0.9038, 0.9252, 0.9570, 0.9981, 0.9428 /
data ((cldnuctab( 7,iw,iconc,4, 5),iw=1,9),iconc=1,7)/    &
0.9814, 0.9981, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9792, 0.9834, 0.9986, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9768, 0.9792, 0.9834, 0.9988, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9742, 0.9742, 0.9768, 0.9834, 0.9986, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9713, 0.9713, 0.9713, 0.9768, 0.9834, 0.9981, 0.9999, 0.9999, 0.9999,  &
0.9570, 0.9570, 0.9570, 0.9610, 0.9681, 0.9834, 0.9993, 0.9997, 0.9999,  &
0.9252, 0.9252, 0.9252, 0.9252, 0.9315, 0.9428, 0.9647, 0.9979, 0.9681 /
data ((cldnuctab( 7,iw,iconc,5, 5),iw=1,9),iconc=1,7)/    &
0.9834, 0.9968, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9792, 0.9834, 0.9972, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9768, 0.9792, 0.9834, 0.9975, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9768, 0.9768, 0.9792, 0.9834, 0.9975, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9742, 0.9742, 0.9768, 0.9792, 0.9834, 0.9972, 0.9999, 0.9999, 0.9999,  &
0.9647, 0.9647, 0.9647, 0.9681, 0.9713, 0.9834, 0.9984, 0.9999, 0.9999,  &
0.9374, 0.9374, 0.9374, 0.9374, 0.9428, 0.9526, 0.9681, 0.9968, 0.9852 /
data ((cldnuctab( 7,iw,iconc,6, 5),iw=1,9),iconc=1,7)/    &
0.9834, 0.9945, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9814, 0.9852, 0.9951, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9792, 0.9814, 0.9852, 0.9957, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9792, 0.9792, 0.9814, 0.9852, 0.9957, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9768, 0.9768, 0.9768, 0.9792, 0.9852, 0.9957, 0.9999, 0.9999, 0.9999,  &
0.9681, 0.9681, 0.9681, 0.9713, 0.9742, 0.9834, 0.9968, 0.9999, 0.9999,  &
0.9428, 0.9428, 0.9428, 0.9479, 0.9479, 0.9570, 0.9713, 0.9945, 0.9975 /
data ((cldnuctab( 7,iw,iconc,7, 5),iw=1,9),iconc=1,7)/    &
0.9852, 0.9919, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9834, 0.9852, 0.9928, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9814, 0.9834, 0.9852, 0.9937, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9814, 0.9814, 0.9814, 0.9852, 0.9945, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9792, 0.9792, 0.9792, 0.9814, 0.9852, 0.9945, 0.9999, 0.9999, 0.9999,  &
0.9713, 0.9713, 0.9713, 0.9742, 0.9768, 0.9834, 0.9951, 0.9999, 0.9999,  &
0.9479, 0.9479, 0.9479, 0.9526, 0.9526, 0.9570, 0.9713, 0.9928, 0.9998 /
data ((cldnuctab( 7,iw,iconc,1, 6),iw=1,9),iconc=1,7)/    &
0.9883, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9792, 0.9883, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9768, 0.9792, 0.9883, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9647, 0.9681, 0.9768, 0.9896, 0.9986, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9374, 0.9428, 0.9479, 0.9570, 0.9834, 0.9999, 0.9968, 0.9999, 0.9999,  &
0.8872, 0.8872, 0.8957, 0.9038, 0.9185, 0.9570, 0.9979, 0.9428, 0.9919,  &
0.8370, 0.8370, 0.8370, 0.8480, 0.8480, 0.8686, 0.9114, 0.9647, 0.9908 /
data ((cldnuctab( 7,iw,iconc,2, 6),iw=1,9),iconc=1,7)/    &
0.9868, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9834, 0.9883, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9792, 0.9814, 0.9883, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9742, 0.9768, 0.9814, 0.9883, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9610, 0.9610, 0.9647, 0.9742, 0.9883, 0.9996, 0.9999, 0.9999, 0.9999,  &
0.9252, 0.9315, 0.9315, 0.9374, 0.9479, 0.9742, 0.9995, 0.9814, 0.9999,  &
0.8781, 0.8781, 0.8781, 0.8781, 0.8872, 0.9038, 0.9374, 0.9896, 0.9185 /
data ((cldnuctab( 7,iw,iconc,3, 6),iw=1,9),iconc=1,7)/    &
0.9868, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9834, 0.9883, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9814, 0.9834, 0.9883, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9792, 0.9814, 0.9834, 0.9883, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9713, 0.9742, 0.9742, 0.9792, 0.9883, 0.9991, 0.9999, 0.9999, 0.9999,  &
0.9526, 0.9526, 0.9526, 0.9570, 0.9647, 0.9834, 0.9995, 0.9963, 0.9999,  &
0.9114, 0.9114, 0.9114, 0.9114, 0.9185, 0.9315, 0.9570, 0.9957, 0.9479 /
data ((cldnuctab( 7,iw,iconc,4, 6),iw=1,9),iconc=1,7)/    &
0.9883, 0.9981, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9852, 0.9883, 0.9984, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9834, 0.9852, 0.9883, 0.9986, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9814, 0.9834, 0.9852, 0.9883, 0.9984, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9768, 0.9792, 0.9792, 0.9834, 0.9883, 0.9984, 0.9999, 0.9999, 0.9999,  &
0.9647, 0.9647, 0.9647, 0.9681, 0.9742, 0.9852, 0.9991, 0.9991, 0.9999,  &
0.9315, 0.9315, 0.9315, 0.9315, 0.9374, 0.9479, 0.9681, 0.9963, 0.9742 /
data ((cldnuctab( 7,iw,iconc,5, 6),iw=1,9),iconc=1,7)/    &
0.9883, 0.9968, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9868, 0.9883, 0.9972, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9852, 0.9868, 0.9883, 0.9972, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9834, 0.9852, 0.9852, 0.9883, 0.9975, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9814, 0.9814, 0.9834, 0.9852, 0.9883, 0.9975, 0.9999, 0.9999, 0.9999,  &
0.9713, 0.9713, 0.9713, 0.9742, 0.9768, 0.9852, 0.9981, 0.9997, 0.9999,  &
0.9428, 0.9428, 0.9479, 0.9479, 0.9479, 0.9570, 0.9713, 0.9951, 0.9945 /
data ((cldnuctab( 7,iw,iconc,6, 6),iw=1,9),iconc=1,7)/    &
0.9896, 0.9945, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9883, 0.9896, 0.9951, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9868, 0.9868, 0.9896, 0.9957, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9852, 0.9868, 0.9868, 0.9896, 0.9963, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9834, 0.9834, 0.9834, 0.9852, 0.9896, 0.9963, 0.9999, 0.9999, 0.9999,  &
0.9742, 0.9742, 0.9768, 0.9768, 0.9792, 0.9868, 0.9972, 0.9999, 0.9999,  &
0.9526, 0.9526, 0.9526, 0.9526, 0.9570, 0.9610, 0.9713, 0.9928, 0.9999 /
data ((cldnuctab( 7,iw,iconc,7, 6),iw=1,9),iconc=1,7)/    &
0.9896, 0.9928, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9883, 0.9896, 0.9945, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9868, 0.9883, 0.9896, 0.9951, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9868, 0.9868, 0.9883, 0.9896, 0.9957, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9852, 0.9852, 0.9852, 0.9868, 0.9896, 0.9957, 0.9999, 0.9999, 0.9999,  &
0.9768, 0.9768, 0.9792, 0.9792, 0.9814, 0.9868, 0.9957, 0.9999, 0.9999,  &
0.9570, 0.9570, 0.9570, 0.9570, 0.9610, 0.9647, 0.9742, 0.9919, 0.9999 /
data ((cldnuctab( 7,iw,iconc,1, 7),iw=1,9),iconc=1,7)/    &
0.9908, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9852, 0.9908, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9814, 0.9852, 0.9908, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9713, 0.9742, 0.9792, 0.9908, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9479, 0.9479, 0.9526, 0.9610, 0.9834, 0.9998, 0.9951, 0.9999, 0.9999,  &
0.9038, 0.9038, 0.9038, 0.9114, 0.9252, 0.9570, 0.9963, 0.9479, 0.9908,  &
0.8586, 0.8586, 0.8586, 0.8686, 0.8686, 0.8872, 0.9185, 0.9610, 0.9868 /
data ((cldnuctab( 7,iw,iconc,2, 7),iw=1,9),iconc=1,7)/    &
0.9908, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9868, 0.9908, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9852, 0.9868, 0.9896, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9792, 0.9814, 0.9852, 0.9908, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9681, 0.9681, 0.9713, 0.9768, 0.9896, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.9374, 0.9374, 0.9374, 0.9428, 0.9526, 0.9742, 0.9990, 0.9792, 0.9997,  &
0.8872, 0.8872, 0.8872, 0.8957, 0.8957, 0.9114, 0.9428, 0.9868, 0.9374 /
data ((cldnuctab( 7,iw,iconc,3, 7),iw=1,9),iconc=1,7)/    &
0.9908, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9883, 0.9908, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9868, 0.9868, 0.9908, 0.9991, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9834, 0.9852, 0.9868, 0.9908, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9768, 0.9768, 0.9792, 0.9834, 0.9908, 0.9996, 0.9999, 0.9999, 0.9999,  &
0.9570, 0.9570, 0.9610, 0.9610, 0.9681, 0.9834, 0.9992, 0.9951, 0.9999,  &
0.9185, 0.9185, 0.9185, 0.9185, 0.9252, 0.9374, 0.9570, 0.9937, 0.9570 /
data ((cldnuctab( 7,iw,iconc,4, 7),iw=1,9),iconc=1,7)/    &
0.9908, 0.9981, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9896, 0.9908, 0.9984, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9883, 0.9883, 0.9908, 0.9984, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9868, 0.9868, 0.9883, 0.9908, 0.9984, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9834, 0.9834, 0.9834, 0.9868, 0.9908, 0.9988, 0.9999, 0.9999, 0.9999,  &
0.9713, 0.9713, 0.9713, 0.9742, 0.9768, 0.9868, 0.9990, 0.9988, 0.9999,  &
0.9374, 0.9374, 0.9374, 0.9428, 0.9428, 0.9526, 0.9681, 0.9945, 0.9852 /
data ((cldnuctab( 7,iw,iconc,5, 7),iw=1,9),iconc=1,7)/    &
0.9919, 0.9963, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9896, 0.9919, 0.9968, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9896, 0.9896, 0.9919, 0.9975, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9883, 0.9883, 0.9896, 0.9908, 0.9975, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9852, 0.9852, 0.9868, 0.9883, 0.9908, 0.9979, 0.9999, 0.9999, 0.9999,  &
0.9768, 0.9768, 0.9768, 0.9792, 0.9814, 0.9883, 0.9981, 0.9997, 0.9999,  &
0.9526, 0.9526, 0.9526, 0.9526, 0.9526, 0.9610, 0.9713, 0.9937, 0.9999 /
data ((cldnuctab( 7,iw,iconc,6, 7),iw=1,9),iconc=1,7)/    &
0.9919, 0.9945, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9908, 0.9919, 0.9957, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9896, 0.9908, 0.9919, 0.9968, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9896, 0.9896, 0.9908, 0.9919, 0.9968, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9868, 0.9868, 0.9883, 0.9896, 0.9919, 0.9968, 0.9999, 0.9999, 0.9999,  &
0.9792, 0.9792, 0.9792, 0.9814, 0.9834, 0.9883, 0.9972, 0.9999, 0.9999,  &
0.9570, 0.9570, 0.9570, 0.9570, 0.9610, 0.9647, 0.9742, 0.9919, 0.9999 /
data ((cldnuctab( 7,iw,iconc,7, 7),iw=1,9),iconc=1,7)/    &
0.9928, 0.9951, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9919, 0.9928, 0.9957, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9908, 0.9919, 0.9928, 0.9963, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9908, 0.9908, 0.9908, 0.9928, 0.9963, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9883, 0.9883, 0.9883, 0.9896, 0.9919, 0.9963, 0.9999, 0.9999, 0.9999,  &
0.9814, 0.9814, 0.9814, 0.9834, 0.9852, 0.9883, 0.9963, 0.9999, 0.9999,  &
0.9610, 0.9610, 0.9610, 0.9610, 0.9647, 0.9681, 0.9742, 0.9908, 0.9999 /
data ((cldnuctab( 8,iw,iconc,1, 1),iw=1,9),iconc=1,7)/    &
0.9759, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8841, 0.9846, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8550, 0.8841, 0.9878, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8330, 0.8443, 0.8928, 0.9759, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8213, 0.8213, 0.8443, 0.8928, 0.9294, 0.9966, 0.9999, 0.9999, 0.9999,  &
0.7693, 0.7693, 0.7830, 0.8090, 0.8841, 0.9011, 0.8443, 0.9702, 0.9999,  &
0.6776, 0.6776, 0.6776, 0.6939, 0.7253, 0.7963, 0.9462, 0.6609, 0.7551 /
data ((cldnuctab( 8,iw,iconc,2, 1),iw=1,9),iconc=1,7)/    &
0.9596, 0.9989, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8749, 0.9702, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8652, 0.8841, 0.9784, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8443, 0.8550, 0.8841, 0.9784, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8330, 0.8443, 0.8550, 0.8841, 0.9555, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.8090, 0.8090, 0.8213, 0.8443, 0.8928, 0.9410, 0.9596, 0.9999, 0.9999,  &
0.7404, 0.7551, 0.7551, 0.7693, 0.7963, 0.8550, 0.9784, 0.7551, 0.8841 /
data ((cldnuctab( 8,iw,iconc,3, 1),iw=1,9),iconc=1,7)/    &
0.9354, 0.9974, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8841, 0.9510, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8652, 0.8841, 0.9635, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8550, 0.8652, 0.8841, 0.9670, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8443, 0.8550, 0.8652, 0.8841, 0.9596, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.8330, 0.8330, 0.8443, 0.8550, 0.8928, 0.9510, 0.9915, 0.9999, 0.9999,  &
0.7963, 0.7963, 0.7963, 0.8090, 0.8330, 0.8841, 0.9759, 0.8550, 0.9702 /
data ((cldnuctab( 8,iw,iconc,4, 1),iw=1,9),iconc=1,7)/    &
0.9088, 0.9949, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8841, 0.9230, 0.9977, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8749, 0.8841, 0.9410, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8652, 0.8749, 0.8928, 0.9510, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8550, 0.8652, 0.8652, 0.8928, 0.9510, 0.9997, 0.9999, 0.9999, 0.9999,  &
0.8550, 0.8550, 0.8550, 0.8652, 0.8841, 0.9510, 0.9955, 0.9999, 0.9999,  &
0.8213, 0.8213, 0.8330, 0.8330, 0.8550, 0.8841, 0.9670, 0.9088, 0.9949 /
data ((cldnuctab( 8,iw,iconc,5, 1),iw=1,9),iconc=1,7)/    &
0.8928, 0.9915, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8841, 0.9011, 0.9961, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8749, 0.8928, 0.9161, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8749, 0.8749, 0.8928, 0.9354, 0.9989, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8652, 0.8652, 0.8749, 0.8928, 0.9410, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.8652, 0.8652, 0.8652, 0.8749, 0.8928, 0.9410, 0.9949, 0.9999, 0.9999,  &
0.8443, 0.8443, 0.8443, 0.8550, 0.8652, 0.8841, 0.9510, 0.9410, 0.9990 /
data ((cldnuctab( 8,iw,iconc,6, 1),iw=1,9),iconc=1,7)/    &
0.8928, 0.9878, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8928, 0.9011, 0.9934, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8841, 0.8928, 0.9161, 0.9961, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8749, 0.8841, 0.8928, 0.9294, 0.9974, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8749, 0.8749, 0.8841, 0.8928, 0.9354, 0.9977, 0.9999, 0.9999, 0.9999,  &
0.8652, 0.8749, 0.8749, 0.8841, 0.8928, 0.9354, 0.9925, 0.9999, 0.9999,  &
0.8550, 0.8550, 0.8550, 0.8550, 0.8652, 0.8928, 0.9354, 0.9635, 0.9996 /
data ((cldnuctab( 8,iw,iconc,7, 1),iw=1,9),iconc=1,7)/    &
0.9011, 0.9827, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8928, 0.9088, 0.9904, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8928, 0.9011, 0.9161, 0.9934, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8841, 0.8928, 0.9011, 0.9230, 0.9955, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8841, 0.8841, 0.8841, 0.9011, 0.9294, 0.9955, 0.9999, 0.9999, 0.9999,  &
0.8749, 0.8749, 0.8841, 0.8841, 0.8928, 0.9294, 0.9904, 0.9999, 0.9999,  &
0.8652, 0.8652, 0.8652, 0.8652, 0.8749, 0.8928, 0.9294, 0.9807, 0.9996 /
data ((cldnuctab( 8,iw,iconc,1, 2),iw=1,9),iconc=1,7)/    &
0.9732, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9354, 0.9807, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9161, 0.9354, 0.9807, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9088, 0.9161, 0.9294, 0.9759, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.8841, 0.8841, 0.9011, 0.9294, 0.9878, 0.9827, 0.9999, 0.9999, 0.9999,  &
0.8213, 0.8213, 0.8330, 0.8550, 0.8928, 0.9759, 0.8550, 0.9462, 0.9999,  &
0.7253, 0.7253, 0.7404, 0.7404, 0.7693, 0.8090, 0.9088, 0.7253, 0.7830 /
data ((cldnuctab( 8,iw,iconc,2, 2),iw=1,9),iconc=1,7)/    &
0.9555, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9354, 0.9670, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9230, 0.9354, 0.9759, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9161, 0.9230, 0.9294, 0.9759, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9011, 0.9088, 0.9161, 0.9354, 0.9759, 0.9980, 0.9999, 0.9999, 0.9999,  &
0.8652, 0.8749, 0.8749, 0.8928, 0.9230, 0.9807, 0.9410, 0.9974, 0.9999,  &
0.7963, 0.7963, 0.7963, 0.8090, 0.8213, 0.8652, 0.9510, 0.7963, 0.8749 /
data ((cldnuctab( 8,iw,iconc,3, 2),iw=1,9),iconc=1,7)/    &
0.9462, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9354, 0.9555, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9294, 0.9354, 0.9670, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9230, 0.9230, 0.9354, 0.9702, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9161, 0.9161, 0.9230, 0.9354, 0.9702, 0.9992, 0.9999, 0.9999, 0.9999,  &
0.9011, 0.9011, 0.9011, 0.9088, 0.9354, 0.9759, 0.9784, 0.9999, 0.9999,  &
0.8443, 0.8550, 0.8550, 0.8550, 0.8749, 0.9011, 0.9670, 0.8841, 0.9462 /
data ((cldnuctab( 8,iw,iconc,4, 2),iw=1,9),iconc=1,7)/    &
0.9462, 0.9970, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9410, 0.9510, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9354, 0.9410, 0.9596, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9294, 0.9294, 0.9410, 0.9635, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9230, 0.9230, 0.9294, 0.9354, 0.9635, 0.9987, 0.9999, 0.9999, 0.9999,  &
0.9088, 0.9161, 0.9161, 0.9230, 0.9354, 0.9670, 0.9878, 0.9999, 0.9999,  &
0.8749, 0.8841, 0.8841, 0.8841, 0.8928, 0.9161, 0.9670, 0.9670, 0.9784 /
data ((cldnuctab( 8,iw,iconc,5, 2),iw=1,9),iconc=1,7)/    &
0.9462, 0.9949, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9410, 0.9510, 0.9970, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9354, 0.9410, 0.9555, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9294, 0.9354, 0.9410, 0.9596, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9294, 0.9294, 0.9354, 0.9410, 0.9596, 0.9977, 0.9999, 0.9999, 0.9999,  &
0.9230, 0.9230, 0.9230, 0.9294, 0.9410, 0.9596, 0.9934, 0.9999, 0.9999,  &
0.8928, 0.8928, 0.9011, 0.9011, 0.9088, 0.9230, 0.9596, 0.9915, 0.9925 /
data ((cldnuctab( 8,iw,iconc,6, 2),iw=1,9),iconc=1,7)/    &
0.9462, 0.9925, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9462, 0.9510, 0.9949, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9410, 0.9462, 0.9555, 0.9966, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9354, 0.9410, 0.9462, 0.9596, 0.9970, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9354, 0.9354, 0.9410, 0.9462, 0.9596, 0.9961, 0.9999, 0.9999, 0.9999,  &
0.9294, 0.9294, 0.9294, 0.9354, 0.9410, 0.9596, 0.9934, 0.9999, 0.9999,  &
0.9088, 0.9088, 0.9088, 0.9088, 0.9161, 0.9294, 0.9555, 0.9934, 0.9961 /
data ((cldnuctab( 8,iw,iconc,7, 2),iw=1,9),iconc=1,7)/    &
0.9510, 0.9892, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9462, 0.9555, 0.9925, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9462, 0.9462, 0.9555, 0.9942, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9410, 0.9410, 0.9462, 0.9596, 0.9949, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9410, 0.9410, 0.9410, 0.9462, 0.9596, 0.9934, 0.9999, 0.9999, 0.9999,  &
0.9354, 0.9354, 0.9354, 0.9410, 0.9462, 0.9596, 0.9925, 0.9999, 0.9999,  &
0.9161, 0.9161, 0.9161, 0.9161, 0.9230, 0.9294, 0.9555, 0.9925, 0.9966 /
data ((cldnuctab( 8,iw,iconc,1, 3),iw=1,9),iconc=1,7)/    &
0.9784, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9670, 0.9827, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9596, 0.9635, 0.9827, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9462, 0.9510, 0.9635, 0.9863, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9230, 0.9230, 0.9294, 0.9510, 0.9846, 0.9732, 0.9994, 0.9999, 0.9999,  &
0.8652, 0.8652, 0.8652, 0.8841, 0.9088, 0.9596, 0.9998, 0.9354, 0.9966,  &
0.7963, 0.7963, 0.7963, 0.7963, 0.8090, 0.8443, 0.9011, 0.9732, 0.8443 /
data ((cldnuctab( 8,iw,iconc,2, 3),iw=1,9),iconc=1,7)/    &
0.9759, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9670, 0.9807, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9635, 0.9635, 0.9807, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9555, 0.9596, 0.9670, 0.9827, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9410, 0.9462, 0.9510, 0.9635, 0.9846, 0.9934, 0.9999, 0.9999, 0.9999,  &
0.9088, 0.9088, 0.9161, 0.9230, 0.9410, 0.9759, 0.9904, 0.9863, 0.9999,  &
0.8443, 0.8443, 0.8443, 0.8550, 0.8652, 0.8841, 0.9410, 0.9977, 0.8928 /
data ((cldnuctab( 8,iw,iconc,3, 3),iw=1,9),iconc=1,7)/    &
0.9759, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9702, 0.9784, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9635, 0.9702, 0.9807, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9635, 0.9635, 0.9670, 0.9807, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9555, 0.9555, 0.9596, 0.9670, 0.9827, 0.9974, 0.9999, 0.9999, 0.9999,  &
0.9354, 0.9354, 0.9354, 0.9410, 0.9555, 0.9807, 0.9974, 0.9989, 0.9999,  &
0.8841, 0.8841, 0.8841, 0.8928, 0.9011, 0.9161, 0.9596, 0.9993, 0.9410 /
data ((cldnuctab( 8,iw,iconc,4, 3),iw=1,9),iconc=1,7)/    &
0.9759, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9732, 0.9784, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9702, 0.9702, 0.9784, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9670, 0.9670, 0.9702, 0.9784, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9635, 0.9635, 0.9635, 0.9702, 0.9784, 0.9977, 0.9999, 0.9999, 0.9999,  &
0.9510, 0.9510, 0.9510, 0.9555, 0.9635, 0.9807, 0.9980, 0.9999, 0.9999,  &
0.9161, 0.9161, 0.9161, 0.9161, 0.9230, 0.9354, 0.9670, 0.9989, 0.9702 /
data ((cldnuctab( 8,iw,iconc,5, 3),iw=1,9),iconc=1,7)/    &
0.9759, 0.9966, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9732, 0.9784, 0.9974, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9702, 0.9732, 0.9784, 0.9977, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9702, 0.9702, 0.9732, 0.9784, 0.9977, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9670, 0.9670, 0.9702, 0.9732, 0.9784, 0.9970, 0.9999, 0.9999, 0.9999,  &
0.9555, 0.9555, 0.9596, 0.9596, 0.9670, 0.9784, 0.9974, 0.9999, 0.9999,  &
0.9294, 0.9294, 0.9294, 0.9294, 0.9354, 0.9462, 0.9670, 0.9977, 0.9827 /
data ((cldnuctab( 8,iw,iconc,6, 3),iw=1,9),iconc=1,7)/    &
0.9784, 0.9942, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9759, 0.9784, 0.9955, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9732, 0.9759, 0.9807, 0.9961, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9732, 0.9732, 0.9732, 0.9807, 0.9961, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9702, 0.9702, 0.9702, 0.9732, 0.9784, 0.9955, 0.9999, 0.9999, 0.9999,  &
0.9596, 0.9635, 0.9635, 0.9635, 0.9702, 0.9784, 0.9961, 0.9999, 0.9999,  &
0.9354, 0.9354, 0.9410, 0.9410, 0.9410, 0.9510, 0.9670, 0.9955, 0.9915 /
data ((cldnuctab( 8,iw,iconc,7, 3),iw=1,9),iconc=1,7)/    &
0.9784, 0.9915, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9759, 0.9807, 0.9925, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9759, 0.9759, 0.9807, 0.9934, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9732, 0.9759, 0.9759, 0.9807, 0.9942, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9732, 0.9732, 0.9732, 0.9759, 0.9807, 0.9942, 0.9999, 0.9999, 0.9999,  &
0.9635, 0.9635, 0.9670, 0.9670, 0.9702, 0.9784, 0.9942, 0.9999, 0.9999,  &
0.9410, 0.9410, 0.9410, 0.9462, 0.9462, 0.9555, 0.9670, 0.9934, 0.9987 /
data ((cldnuctab( 8,iw,iconc,1, 4),iw=1,9),iconc=1,7)/    &
0.9904, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9846, 0.9904, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9807, 0.9846, 0.9904, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9702, 0.9732, 0.9784, 0.9915, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9462, 0.9462, 0.9510, 0.9635, 0.9827, 0.9998, 0.9955, 0.9999, 0.9999,  &
0.9011, 0.9011, 0.9011, 0.9088, 0.9230, 0.9555, 0.9966, 0.9462, 0.9915,  &
0.8550, 0.8550, 0.8652, 0.8652, 0.8652, 0.8841, 0.9161, 0.9635, 0.9878 /
data ((cldnuctab( 8,iw,iconc,2, 4),iw=1,9),iconc=1,7)/    &
0.9904, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9846, 0.9904, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9827, 0.9846, 0.9892, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9784, 0.9807, 0.9846, 0.9904, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9670, 0.9670, 0.9702, 0.9759, 0.9892, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.9354, 0.9354, 0.9354, 0.9410, 0.9510, 0.9759, 0.9992, 0.9807, 0.9998,  &
0.8928, 0.8928, 0.8928, 0.8928, 0.9011, 0.9088, 0.9410, 0.9878, 0.9294 /
data ((cldnuctab( 8,iw,iconc,3, 4),iw=1,9),iconc=1,7)/    &
0.9904, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9878, 0.9904, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9863, 0.9863, 0.9904, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9827, 0.9846, 0.9863, 0.9904, 0.9989, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9759, 0.9759, 0.9784, 0.9827, 0.9904, 0.9995, 0.9999, 0.9999, 0.9999,  &
0.9555, 0.9555, 0.9596, 0.9596, 0.9670, 0.9827, 0.9993, 0.9955, 0.9999,  &
0.9161, 0.9161, 0.9161, 0.9230, 0.9230, 0.9354, 0.9596, 0.9942, 0.9555 /
data ((cldnuctab( 8,iw,iconc,4, 4),iw=1,9),iconc=1,7)/    &
0.9904, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9892, 0.9904, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9863, 0.9878, 0.9904, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9863, 0.9863, 0.9878, 0.9904, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9807, 0.9827, 0.9827, 0.9863, 0.9904, 0.9987, 0.9999, 0.9999, 0.9999,  &
0.9702, 0.9702, 0.9702, 0.9732, 0.9759, 0.9863, 0.9989, 0.9987, 0.9999,  &
0.9410, 0.9410, 0.9410, 0.9410, 0.9410, 0.9510, 0.9670, 0.9949, 0.9807 /
data ((cldnuctab( 8,iw,iconc,5, 4),iw=1,9),iconc=1,7)/    &
0.9904, 0.9966, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9892, 0.9915, 0.9966, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9878, 0.9892, 0.9915, 0.9974, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9878, 0.9878, 0.9892, 0.9904, 0.9974, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9846, 0.9846, 0.9863, 0.9878, 0.9904, 0.9977, 0.9999, 0.9999, 0.9999,  &
0.9759, 0.9759, 0.9759, 0.9784, 0.9807, 0.9878, 0.9980, 0.9997, 0.9999,  &
0.9510, 0.9510, 0.9510, 0.9510, 0.9555, 0.9596, 0.9732, 0.9942, 0.9997 /
data ((cldnuctab( 8,iw,iconc,6, 4),iw=1,9),iconc=1,7)/    &
0.9915, 0.9942, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9904, 0.9915, 0.9955, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9892, 0.9904, 0.9915, 0.9966, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9892, 0.9892, 0.9892, 0.9915, 0.9966, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9863, 0.9863, 0.9878, 0.9892, 0.9915, 0.9970, 0.9999, 0.9999, 0.9999,  &
0.9784, 0.9784, 0.9784, 0.9807, 0.9827, 0.9878, 0.9970, 0.9999, 0.9999,  &
0.9555, 0.9555, 0.9555, 0.9596, 0.9596, 0.9635, 0.9732, 0.9925, 0.9999 /
data ((cldnuctab( 8,iw,iconc,7, 4),iw=1,9),iconc=1,7)/    &
0.9915, 0.9942, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9915, 0.9925, 0.9955, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9904, 0.9904, 0.9925, 0.9961, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9892, 0.9904, 0.9904, 0.9915, 0.9961, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9878, 0.9878, 0.9878, 0.9892, 0.9915, 0.9961, 0.9999, 0.9999, 0.9999,  &
0.9807, 0.9807, 0.9807, 0.9827, 0.9846, 0.9892, 0.9961, 0.9999, 0.9999,  &
0.9596, 0.9596, 0.9596, 0.9596, 0.9635, 0.9670, 0.9759, 0.9904, 0.9999 /
data ((cldnuctab( 8,iw,iconc,1, 5),iw=1,9),iconc=1,7)/    &
0.9942, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9904, 0.9934, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9878, 0.9904, 0.9942, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9784, 0.9807, 0.9846, 0.9925, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9555, 0.9555, 0.9596, 0.9670, 0.9827, 0.9994, 0.9942, 0.9999, 0.9999,  &
0.9230, 0.9230, 0.9230, 0.9294, 0.9410, 0.9596, 0.9925, 0.9635, 0.9892,  &
0.8928, 0.8928, 0.8928, 0.8928, 0.9011, 0.9088, 0.9354, 0.9670, 0.9846 /
data ((cldnuctab( 8,iw,iconc,2, 5),iw=1,9),iconc=1,7)/    &
0.9942, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9915, 0.9942, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9904, 0.9915, 0.9942, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9863, 0.9878, 0.9892, 0.9942, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9732, 0.9759, 0.9759, 0.9807, 0.9904, 0.9996, 0.9995, 0.9999, 0.9999,  &
0.9462, 0.9462, 0.9510, 0.9510, 0.9596, 0.9759, 0.9974, 0.9827, 0.9990,  &
0.9161, 0.9161, 0.9161, 0.9161, 0.9161, 0.9294, 0.9510, 0.9846, 0.9997 /
data ((cldnuctab( 8,iw,iconc,3, 5),iw=1,9),iconc=1,7)/    &
0.9942, 0.9989, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9925, 0.9942, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9915, 0.9925, 0.9942, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9892, 0.9904, 0.9915, 0.9942, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9827, 0.9846, 0.9846, 0.9878, 0.9925, 0.9995, 0.9999, 0.9999, 0.9999,  &
0.9670, 0.9670, 0.9670, 0.9702, 0.9732, 0.9846, 0.9985, 0.9955, 0.9999,  &
0.9354, 0.9354, 0.9354, 0.9354, 0.9410, 0.9462, 0.9635, 0.9904, 0.9999 /
data ((cldnuctab( 8,iw,iconc,4, 5),iw=1,9),iconc=1,7)/    &
0.9949, 0.9977, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9934, 0.9949, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9925, 0.9934, 0.9949, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9915, 0.9915, 0.9925, 0.9942, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9878, 0.9878, 0.9892, 0.9904, 0.9934, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.9759, 0.9759, 0.9759, 0.9784, 0.9807, 0.9878, 0.9983, 0.9995, 0.9999,  &
0.9510, 0.9510, 0.9510, 0.9510, 0.9555, 0.9596, 0.9702, 0.9925, 0.9999 /
data ((cldnuctab( 8,iw,iconc,5, 5),iw=1,9),iconc=1,7)/    &
0.9949, 0.9970, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9942, 0.9949, 0.9977, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9934, 0.9942, 0.9949, 0.9977, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9925, 0.9934, 0.9934, 0.9949, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9904, 0.9904, 0.9904, 0.9915, 0.9942, 0.9983, 0.9999, 0.9999, 0.9999,  &
0.9807, 0.9807, 0.9827, 0.9827, 0.9846, 0.9904, 0.9977, 0.9999, 0.9999,  &
0.9596, 0.9596, 0.9596, 0.9596, 0.9635, 0.9670, 0.9759, 0.9915, 0.9999 /
data ((cldnuctab( 8,iw,iconc,6, 5),iw=1,9),iconc=1,7)/    &
0.9955, 0.9970, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9949, 0.9955, 0.9974, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9942, 0.9942, 0.9955, 0.9974, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9934, 0.9934, 0.9942, 0.9949, 0.9977, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9915, 0.9915, 0.9915, 0.9925, 0.9949, 0.9977, 0.9999, 0.9999, 0.9999,  &
0.9846, 0.9846, 0.9846, 0.9846, 0.9878, 0.9904, 0.9970, 0.9999, 0.9999,  &
0.9635, 0.9635, 0.9670, 0.9670, 0.9670, 0.9702, 0.9784, 0.9904, 0.9999 /
data ((cldnuctab( 8,iw,iconc,7, 5),iw=1,9),iconc=1,7)/    &
0.9955, 0.9970, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9949, 0.9955, 0.9974, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9949, 0.9949, 0.9955, 0.9974, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9942, 0.9942, 0.9949, 0.9955, 0.9974, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9925, 0.9925, 0.9925, 0.9934, 0.9949, 0.9974, 0.9999, 0.9999, 0.9999,  &
0.9863, 0.9863, 0.9863, 0.9863, 0.9878, 0.9915, 0.9966, 0.9999, 0.9999,  &
0.9670, 0.9670, 0.9670, 0.9702, 0.9702, 0.9732, 0.9784, 0.9904, 0.9999 /
data ((cldnuctab( 8,iw,iconc,1, 6),iw=1,9),iconc=1,7)/    &
0.9955, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9942, 0.9955, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9915, 0.9925, 0.9955, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9827, 0.9846, 0.9878, 0.9934, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9635, 0.9635, 0.9670, 0.9732, 0.9827, 0.9987, 0.9934, 0.9999, 0.9999,  &
0.9354, 0.9354, 0.9354, 0.9410, 0.9462, 0.9635, 0.9904, 0.9999, 0.9904,  &
0.9161, 0.9161, 0.9161, 0.9161, 0.9230, 0.9294, 0.9462, 0.9702, 0.9846 /
data ((cldnuctab( 8,iw,iconc,2, 6),iw=1,9),iconc=1,7)/    &
0.9961, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9949, 0.9955, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9934, 0.9942, 0.9961, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9892, 0.9904, 0.9925, 0.9955, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9784, 0.9784, 0.9807, 0.9846, 0.9904, 0.9994, 0.9992, 0.9999, 0.9999,  &
0.9555, 0.9555, 0.9555, 0.9596, 0.9635, 0.9759, 0.9961, 0.9999, 0.9985,  &
0.9294, 0.9294, 0.9294, 0.9294, 0.9354, 0.9410, 0.9555, 0.9827, 0.9992 /
data ((cldnuctab( 8,iw,iconc,3, 6),iw=1,9),iconc=1,7)/    &
0.9961, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9949, 0.9961, 0.9989, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9942, 0.9949, 0.9961, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9925, 0.9934, 0.9942, 0.9961, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9863, 0.9878, 0.9878, 0.9904, 0.9942, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.9702, 0.9702, 0.9732, 0.9732, 0.9784, 0.9846, 0.9977, 0.9999, 0.9999,  &
0.9462, 0.9462, 0.9462, 0.9462, 0.9510, 0.9555, 0.9670, 0.9892, 0.9999 /
data ((cldnuctab( 8,iw,iconc,4, 6),iw=1,9),iconc=1,7)/    &
0.9966, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9955, 0.9966, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9955, 0.9955, 0.9961, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9942, 0.9942, 0.9949, 0.9961, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9904, 0.9904, 0.9915, 0.9925, 0.9949, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.9807, 0.9807, 0.9807, 0.9827, 0.9846, 0.9892, 0.9980, 0.9999, 0.9999,  &
0.9596, 0.9596, 0.9596, 0.9596, 0.9596, 0.9635, 0.9732, 0.9904, 0.9999 /
data ((cldnuctab( 8,iw,iconc,5, 6),iw=1,9),iconc=1,7)/    &
0.9966, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9961, 0.9966, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9961, 0.9961, 0.9966, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9949, 0.9955, 0.9955, 0.9966, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9925, 0.9925, 0.9934, 0.9942, 0.9955, 0.9987, 0.9999, 0.9999, 0.9999,  &
0.9846, 0.9846, 0.9846, 0.9863, 0.9878, 0.9915, 0.9977, 0.9999, 0.9999,  &
0.9670, 0.9670, 0.9670, 0.9670, 0.9670, 0.9702, 0.9784, 0.9904, 0.9999 /
data ((cldnuctab( 8,iw,iconc,6, 6),iw=1,9),iconc=1,7)/    &
0.9970, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9966, 0.9970, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9961, 0.9966, 0.9970, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9955, 0.9961, 0.9961, 0.9970, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9934, 0.9942, 0.9942, 0.9949, 0.9961, 0.9983, 0.9999, 0.9999, 0.9999,  &
0.9878, 0.9878, 0.9878, 0.9878, 0.9892, 0.9925, 0.9974, 0.9999, 0.9999,  &
0.9702, 0.9702, 0.9702, 0.9702, 0.9702, 0.9732, 0.9784, 0.9904, 0.9999 /
data ((cldnuctab( 8,iw,iconc,7, 6),iw=1,9),iconc=1,7)/    &
0.9974, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9970, 0.9974, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9966, 0.9970, 0.9970, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9961, 0.9961, 0.9966, 0.9970, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9942, 0.9942, 0.9949, 0.9955, 0.9966, 0.9983, 0.9999, 0.9999, 0.9999,  &
0.9892, 0.9892, 0.9892, 0.9892, 0.9904, 0.9925, 0.9970, 0.9999, 0.9999,  &
0.9732, 0.9732, 0.9732, 0.9732, 0.9732, 0.9759, 0.9807, 0.9904, 0.9998 /
data ((cldnuctab( 8,iw,iconc,1, 7),iw=1,9),iconc=1,7)/    &
0.9970, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9955, 0.9966, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9925, 0.9942, 0.9966, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9846, 0.9863, 0.9892, 0.9942, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9670, 0.9702, 0.9702, 0.9759, 0.9846, 0.9980, 0.9942, 0.9999, 0.9999,  &
0.9462, 0.9462, 0.9462, 0.9510, 0.9555, 0.9670, 0.9904, 0.9999, 0.9925,  &
0.9294, 0.9294, 0.9294, 0.9354, 0.9354, 0.9410, 0.9510, 0.9732, 0.9846 /
data ((cldnuctab( 8,iw,iconc,2, 7),iw=1,9),iconc=1,7)/    &
0.9970, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9961, 0.9970, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9949, 0.9955, 0.9970, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9915, 0.9925, 0.9934, 0.9961, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9807, 0.9827, 0.9827, 0.9863, 0.9915, 0.9992, 0.9993, 0.9999, 0.9999,  &
0.9635, 0.9635, 0.9635, 0.9635, 0.9702, 0.9784, 0.9949, 0.9999, 0.9983,  &
0.9410, 0.9410, 0.9410, 0.9410, 0.9462, 0.9510, 0.9635, 0.9827, 0.9985 /
data ((cldnuctab( 8,iw,iconc,3, 7),iw=1,9),iconc=1,7)/    &
0.9974, 0.9989, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9966, 0.9974, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9961, 0.9966, 0.9970, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9942, 0.9949, 0.9955, 0.9970, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9892, 0.9892, 0.9904, 0.9915, 0.9949, 0.9993, 0.9999, 0.9999, 0.9999,  &
0.9759, 0.9759, 0.9759, 0.9759, 0.9807, 0.9863, 0.9970, 0.9999, 0.9999,  &
0.9555, 0.9555, 0.9555, 0.9555, 0.9555, 0.9596, 0.9702, 0.9878, 0.9999 /
data ((cldnuctab( 8,iw,iconc,4, 7),iw=1,9),iconc=1,7)/    &
0.9974, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9970, 0.9974, 0.9989, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9966, 0.9970, 0.9974, 0.9989, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9955, 0.9961, 0.9966, 0.9974, 0.9989, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9925, 0.9925, 0.9934, 0.9942, 0.9961, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.9827, 0.9827, 0.9827, 0.9846, 0.9863, 0.9904, 0.9974, 0.9999, 0.9999,  &
0.9635, 0.9635, 0.9635, 0.9635, 0.9670, 0.9702, 0.9759, 0.9904, 0.9999 /
data ((cldnuctab( 8,iw,iconc,5, 7),iw=1,9),iconc=1,7)/    &
0.9977, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9974, 0.9977, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9970, 0.9974, 0.9977, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9966, 0.9966, 0.9970, 0.9974, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9942, 0.9942, 0.9942, 0.9949, 0.9966, 0.9989, 0.9999, 0.9999, 0.9999,  &
0.9863, 0.9863, 0.9878, 0.9878, 0.9892, 0.9925, 0.9974, 0.9999, 0.9999,  &
0.9702, 0.9702, 0.9702, 0.9702, 0.9702, 0.9732, 0.9807, 0.9904, 0.9999 /
data ((cldnuctab( 8,iw,iconc,6, 7),iw=1,9),iconc=1,7)/    &
0.9980, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9977, 0.9980, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9974, 0.9974, 0.9977, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9970, 0.9970, 0.9974, 0.9977, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9949, 0.9949, 0.9955, 0.9961, 0.9970, 0.9987, 0.9999, 0.9999, 0.9999,  &
0.9892, 0.9892, 0.9892, 0.9892, 0.9904, 0.9934, 0.9974, 0.9999, 0.9999,  &
0.9732, 0.9732, 0.9732, 0.9732, 0.9759, 0.9759, 0.9807, 0.9904, 0.9998 /
data ((cldnuctab( 8,iw,iconc,7, 7),iw=1,9),iconc=1,7)/    &
0.9980, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9977, 0.9980, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9977, 0.9977, 0.9980, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9970, 0.9974, 0.9974, 0.9980, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9955, 0.9955, 0.9961, 0.9961, 0.9970, 0.9987, 0.9999, 0.9999, 0.9999,  &
0.9904, 0.9904, 0.9904, 0.9904, 0.9915, 0.9934, 0.9970, 0.9999, 0.9999,  &
0.9759, 0.9759, 0.9759, 0.9759, 0.9759, 0.9784, 0.9827, 0.9904, 0.9995 /
data ((cldnuctab( 9,iw,iconc,1, 1),iw=1,9),iconc=1,7)/    &
0.9756, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9592, 0.9804, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9505, 0.9592, 0.9825, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9404, 0.9457, 0.9592, 0.9844, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9154, 0.9223, 0.9287, 0.9505, 0.9861, 0.9756, 0.9997, 0.9999, 0.9999,  &
0.8642, 0.8642, 0.8739, 0.8832, 0.9080, 0.9666, 0.9080, 0.9457, 0.9985,  &
0.8077, 0.8077, 0.8077, 0.8077, 0.8200, 0.8539, 0.9080, 0.9804, 0.8539 /
data ((cldnuctab( 9,iw,iconc,2, 1),iw=1,9),iconc=1,7)/    &
0.9729, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9631, 0.9782, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9551, 0.9631, 0.9804, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9505, 0.9551, 0.9592, 0.9804, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9404, 0.9404, 0.9457, 0.9592, 0.9844, 0.9941, 0.9999, 0.9999, 0.9999,  &
0.9080, 0.9080, 0.9080, 0.9223, 0.9404, 0.9782, 0.9592, 0.9914, 0.9999,  &
0.8539, 0.8539, 0.8539, 0.8539, 0.8642, 0.8919, 0.9457, 0.9992, 0.9002 /
data ((cldnuctab( 9,iw,iconc,3, 1),iw=1,9),iconc=1,7)/    &
0.9699, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9631, 0.9756, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9592, 0.9631, 0.9756, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9551, 0.9592, 0.9631, 0.9782, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9505, 0.9505, 0.9551, 0.9631, 0.9782, 0.9977, 0.9999, 0.9999, 0.9999,  &
0.9287, 0.9348, 0.9348, 0.9404, 0.9551, 0.9804, 0.9903, 0.9996, 0.9999,  &
0.8919, 0.8919, 0.8919, 0.8919, 0.9002, 0.9223, 0.9631, 0.9997, 0.9457 /
data ((cldnuctab( 9,iw,iconc,4, 1),iw=1,9),iconc=1,7)/    &
0.9699, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9666, 0.9729, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9631, 0.9666, 0.9756, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9592, 0.9631, 0.9631, 0.9756, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9551, 0.9592, 0.9592, 0.9666, 0.9756, 0.9977, 0.9999, 0.9999, 0.9999,  &
0.9457, 0.9457, 0.9457, 0.9505, 0.9592, 0.9782, 0.9965, 0.9999, 0.9999,  &
0.9154, 0.9154, 0.9154, 0.9154, 0.9223, 0.9348, 0.9666, 0.9994, 0.9756 /
data ((cldnuctab( 9,iw,iconc,5, 1),iw=1,9),iconc=1,7)/    &
0.9729, 0.9960, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9699, 0.9729, 0.9974, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9666, 0.9699, 0.9756, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9631, 0.9631, 0.9666, 0.9756, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9592, 0.9631, 0.9631, 0.9666, 0.9756, 0.9970, 0.9999, 0.9999, 0.9999,  &
0.9505, 0.9505, 0.9551, 0.9551, 0.9631, 0.9756, 0.9965, 0.9999, 0.9999,  &
0.9287, 0.9287, 0.9287, 0.9287, 0.9348, 0.9457, 0.9666, 0.9985, 0.9861 /
data ((cldnuctab( 9,iw,iconc,6, 1),iw=1,9),iconc=1,7)/    &
0.9729, 0.9941, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9699, 0.9756, 0.9955, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9666, 0.9699, 0.9756, 0.9960, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9666, 0.9666, 0.9699, 0.9756, 0.9965, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9631, 0.9631, 0.9666, 0.9699, 0.9756, 0.9955, 0.9999, 0.9999, 0.9999,  &
0.9551, 0.9592, 0.9592, 0.9592, 0.9666, 0.9756, 0.9955, 0.9999, 0.9999,  &
0.9348, 0.9348, 0.9348, 0.9348, 0.9404, 0.9505, 0.9666, 0.9965, 0.9933 /
data ((cldnuctab( 9,iw,iconc,7, 1),iw=1,9),iconc=1,7)/    &
0.9756, 0.9914, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9729, 0.9756, 0.9933, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9699, 0.9729, 0.9756, 0.9941, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9699, 0.9699, 0.9729, 0.9756, 0.9941, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9666, 0.9666, 0.9699, 0.9699, 0.9756, 0.9941, 0.9999, 0.9999, 0.9999,  &
0.9592, 0.9592, 0.9631, 0.9631, 0.9666, 0.9756, 0.9941, 0.9999, 0.9999,  &
0.9404, 0.9404, 0.9404, 0.9404, 0.9457, 0.9505, 0.9666, 0.9941, 0.9965 /
data ((cldnuctab( 9,iw,iconc,1, 2),iw=1,9),iconc=1,7)/    &
0.9876, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9804, 0.9890, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9782, 0.9804, 0.9876, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9666, 0.9699, 0.9782, 0.9903, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9404, 0.9457, 0.9505, 0.9592, 0.9844, 0.9999, 0.9970, 0.9999, 0.9999,  &
0.9002, 0.9002, 0.9002, 0.9080, 0.9223, 0.9592, 0.9983, 0.9505, 0.9924,  &
0.8539, 0.8539, 0.8642, 0.8642, 0.8642, 0.8832, 0.9154, 0.9666, 0.9914 /
data ((cldnuctab( 9,iw,iconc,2, 2),iw=1,9),iconc=1,7)/    &
0.9876, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9825, 0.9876, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9804, 0.9825, 0.9890, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9756, 0.9782, 0.9804, 0.9890, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9631, 0.9631, 0.9666, 0.9756, 0.9890, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.9348, 0.9348, 0.9348, 0.9404, 0.9505, 0.9756, 0.9996, 0.9844, 0.9999,  &
0.8832, 0.8919, 0.8919, 0.8919, 0.9002, 0.9080, 0.9404, 0.9914, 0.9287 /
data ((cldnuctab( 9,iw,iconc,3, 2),iw=1,9),iconc=1,7)/    &
0.9876, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9844, 0.9876, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9825, 0.9825, 0.9890, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9804, 0.9804, 0.9825, 0.9876, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9729, 0.9729, 0.9756, 0.9804, 0.9890, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.9551, 0.9551, 0.9551, 0.9592, 0.9666, 0.9825, 0.9995, 0.9965, 0.9999,  &
0.9154, 0.9154, 0.9154, 0.9154, 0.9223, 0.9348, 0.9592, 0.9965, 0.9551 /
data ((cldnuctab( 9,iw,iconc,4, 2),iw=1,9),iconc=1,7)/    &
0.9876, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9861, 0.9890, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9844, 0.9844, 0.9890, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9825, 0.9825, 0.9844, 0.9890, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9782, 0.9782, 0.9804, 0.9825, 0.9890, 0.9985, 0.9999, 0.9999, 0.9999,  &
0.9666, 0.9666, 0.9666, 0.9699, 0.9756, 0.9861, 0.9992, 0.9994, 0.9999,  &
0.9348, 0.9348, 0.9348, 0.9404, 0.9404, 0.9505, 0.9699, 0.9965, 0.9729 /
data ((cldnuctab( 9,iw,iconc,5, 2),iw=1,9),iconc=1,7)/    &
0.9890, 0.9965, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9861, 0.9890, 0.9970, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9861, 0.9861, 0.9890, 0.9974, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9844, 0.9844, 0.9861, 0.9890, 0.9974, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9825, 0.9825, 0.9825, 0.9844, 0.9890, 0.9977, 0.9999, 0.9999, 0.9999,  &
0.9729, 0.9729, 0.9729, 0.9756, 0.9782, 0.9861, 0.9983, 0.9997, 0.9999,  &
0.9457, 0.9505, 0.9505, 0.9505, 0.9505, 0.9592, 0.9729, 0.9955, 0.9914 /
data ((cldnuctab( 9,iw,iconc,6, 2),iw=1,9),iconc=1,7)/    &
0.9890, 0.9948, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9876, 0.9890, 0.9948, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9861, 0.9876, 0.9890, 0.9960, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9861, 0.9861, 0.9876, 0.9890, 0.9965, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9844, 0.9844, 0.9844, 0.9861, 0.9890, 0.9965, 0.9999, 0.9999, 0.9999,  &
0.9756, 0.9756, 0.9756, 0.9782, 0.9804, 0.9861, 0.9970, 0.9999, 0.9999,  &
0.9551, 0.9551, 0.9551, 0.9551, 0.9592, 0.9631, 0.9729, 0.9933, 0.9998 /
data ((cldnuctab( 9,iw,iconc,7, 2),iw=1,9),iconc=1,7)/    &
0.9903, 0.9924, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9890, 0.9903, 0.9941, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9876, 0.9890, 0.9903, 0.9955, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9876, 0.9876, 0.9876, 0.9903, 0.9955, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9844, 0.9861, 0.9861, 0.9876, 0.9903, 0.9955, 0.9999, 0.9999, 0.9999,  &
0.9782, 0.9782, 0.9782, 0.9804, 0.9825, 0.9876, 0.9960, 0.9999, 0.9999,  &
0.9592, 0.9592, 0.9592, 0.9592, 0.9592, 0.9666, 0.9756, 0.9914, 0.9999 /
data ((cldnuctab( 9,iw,iconc,1, 3),iw=1,9),iconc=1,7)/    &
0.9948, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9924, 0.9941, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9890, 0.9914, 0.9948, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9804, 0.9825, 0.9861, 0.9933, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9592, 0.9631, 0.9631, 0.9699, 0.9844, 0.9992, 0.9933, 0.9999, 0.9999,  &
0.9348, 0.9348, 0.9348, 0.9348, 0.9457, 0.9631, 0.9924, 0.9999, 0.9903,  &
0.9080, 0.9080, 0.9080, 0.9154, 0.9154, 0.9223, 0.9404, 0.9699, 0.9861 /
data ((cldnuctab( 9,iw,iconc,2, 3),iw=1,9),iconc=1,7)/    &
0.9948, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9933, 0.9948, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9914, 0.9924, 0.9948, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9876, 0.9890, 0.9903, 0.9948, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9756, 0.9782, 0.9782, 0.9825, 0.9903, 0.9995, 0.9995, 0.9999, 0.9999,  &
0.9551, 0.9551, 0.9551, 0.9551, 0.9631, 0.9756, 0.9970, 0.9861, 0.9989,  &
0.9223, 0.9223, 0.9287, 0.9287, 0.9287, 0.9348, 0.9551, 0.9844, 0.9996 /
data ((cldnuctab( 9,iw,iconc,3, 3),iw=1,9),iconc=1,7)/    &
0.9955, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9933, 0.9955, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9933, 0.9933, 0.9948, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9914, 0.9914, 0.9924, 0.9948, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9844, 0.9861, 0.9861, 0.9890, 0.9933, 0.9995, 0.9999, 0.9999, 0.9999,  &
0.9699, 0.9699, 0.9699, 0.9729, 0.9756, 0.9844, 0.9983, 0.9960, 0.9999,  &
0.9404, 0.9404, 0.9404, 0.9457, 0.9457, 0.9505, 0.9666, 0.9903, 0.9999 /
data ((cldnuctab( 9,iw,iconc,4, 3),iw=1,9),iconc=1,7)/    &
0.9955, 0.9977, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9948, 0.9955, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9941, 0.9941, 0.9955, 0.9985, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9924, 0.9933, 0.9941, 0.9955, 0.9987, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9890, 0.9890, 0.9903, 0.9914, 0.9941, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.9782, 0.9782, 0.9782, 0.9804, 0.9825, 0.9890, 0.9983, 0.9999, 0.9999,  &
0.9551, 0.9551, 0.9551, 0.9551, 0.9592, 0.9631, 0.9729, 0.9914, 0.9999 /
data ((cldnuctab( 9,iw,iconc,5, 3),iw=1,9),iconc=1,7)/    &
0.9960, 0.9974, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9948, 0.9960, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9948, 0.9948, 0.9960, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9941, 0.9941, 0.9948, 0.9955, 0.9983, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9914, 0.9914, 0.9924, 0.9933, 0.9948, 0.9985, 0.9999, 0.9999, 0.9999,  &
0.9825, 0.9825, 0.9844, 0.9844, 0.9861, 0.9903, 0.9977, 0.9999, 0.9999,  &
0.9631, 0.9631, 0.9631, 0.9631, 0.9666, 0.9699, 0.9756, 0.9914, 0.9999 /
data ((cldnuctab( 9,iw,iconc,6, 3),iw=1,9),iconc=1,7)/    &
0.9960, 0.9974, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9955, 0.9960, 0.9977, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9955, 0.9955, 0.9960, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9948, 0.9948, 0.9955, 0.9960, 0.9980, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9924, 0.9924, 0.9933, 0.9941, 0.9955, 0.9980, 0.9999, 0.9999, 0.9999,  &
0.9861, 0.9861, 0.9861, 0.9861, 0.9876, 0.9914, 0.9974, 0.9999, 0.9999,  &
0.9699, 0.9699, 0.9699, 0.9699, 0.9699, 0.9729, 0.9782, 0.9914, 0.9999 /
data ((cldnuctab( 9,iw,iconc,7, 3),iw=1,9),iconc=1,7)/    &
0.9965, 0.9974, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9960, 0.9965, 0.9977, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9955, 0.9960, 0.9965, 0.9977, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9948, 0.9955, 0.9955, 0.9960, 0.9977, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9933, 0.9933, 0.9933, 0.9941, 0.9955, 0.9980, 0.9999, 0.9999, 0.9999,  &
0.9876, 0.9876, 0.9876, 0.9876, 0.9890, 0.9924, 0.9970, 0.9999, 0.9999,  &
0.9699, 0.9699, 0.9699, 0.9729, 0.9729, 0.9756, 0.9804, 0.9903, 0.9999 /
data ((cldnuctab( 9,iw,iconc,1, 4),iw=1,9),iconc=1,7)/    &
0.9977, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9970, 0.9980, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9948, 0.9955, 0.9974, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9876, 0.9890, 0.9914, 0.9948, 0.9995, 0.9998, 0.9999, 0.9999, 0.9999,  &
0.9756, 0.9756, 0.9756, 0.9804, 0.9861, 0.9974, 0.9999, 0.9997, 0.9999,  &
0.9592, 0.9592, 0.9592, 0.9631, 0.9666, 0.9756, 0.9903, 0.9996, 0.9999,  &
0.9505, 0.9505, 0.9505, 0.9505, 0.9551, 0.9592, 0.9666, 0.9782, 0.9876 /
data ((cldnuctab( 9,iw,iconc,2, 4),iw=1,9),iconc=1,7)/    &
0.9983, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9977, 0.9980, 0.9994, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9965, 0.9970, 0.9980, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9933, 0.9941, 0.9948, 0.9970, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9844, 0.9861, 0.9861, 0.9876, 0.9924, 0.9987, 0.9999, 0.9999, 0.9999,  &
0.9699, 0.9699, 0.9699, 0.9729, 0.9756, 0.9825, 0.9941, 0.9999, 0.9980,  &
0.9592, 0.9592, 0.9592, 0.9592, 0.9592, 0.9631, 0.9699, 0.9861, 0.9977 /
data ((cldnuctab( 9,iw,iconc,3, 4),iw=1,9),iconc=1,7)/    &
0.9983, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9980, 0.9983, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9974, 0.9977, 0.9983, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9960, 0.9960, 0.9965, 0.9977, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9914, 0.9914, 0.9924, 0.9933, 0.9955, 0.9992, 0.9999, 0.9999, 0.9999,  &
0.9804, 0.9804, 0.9804, 0.9804, 0.9825, 0.9876, 0.9965, 0.9999, 0.9998,  &
0.9666, 0.9666, 0.9666, 0.9666, 0.9666, 0.9699, 0.9756, 0.9890, 0.9997 /
data ((cldnuctab( 9,iw,iconc,4, 4),iw=1,9),iconc=1,7)/    &
0.9985, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9983, 0.9985, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9980, 0.9980, 0.9985, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9970, 0.9974, 0.9974, 0.9983, 0.9993, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9941, 0.9941, 0.9948, 0.9955, 0.9970, 0.9992, 0.9999, 0.9999, 0.9999,  &
0.9861, 0.9861, 0.9861, 0.9876, 0.9890, 0.9914, 0.9974, 0.9999, 0.9999,  &
0.9729, 0.9729, 0.9729, 0.9729, 0.9729, 0.9756, 0.9804, 0.9903, 0.9999 /
data ((cldnuctab( 9,iw,iconc,5, 4),iw=1,9),iconc=1,7)/    &
0.9987, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9983, 0.9985, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9983, 0.9983, 0.9985, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9977, 0.9977, 0.9980, 0.9985, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9955, 0.9955, 0.9960, 0.9965, 0.9974, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.9890, 0.9890, 0.9890, 0.9903, 0.9914, 0.9933, 0.9974, 0.9999, 0.9999,  &
0.9756, 0.9756, 0.9756, 0.9756, 0.9782, 0.9782, 0.9825, 0.9914, 0.9998 /
data ((cldnuctab( 9,iw,iconc,6, 4),iw=1,9),iconc=1,7)/    &
0.9987, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9985, 0.9987, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9985, 0.9985, 0.9987, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9980, 0.9980, 0.9983, 0.9985, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9965, 0.9965, 0.9965, 0.9970, 0.9977, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.9914, 0.9914, 0.9914, 0.9914, 0.9924, 0.9941, 0.9974, 0.9999, 0.9999,  &
0.9782, 0.9782, 0.9782, 0.9804, 0.9804, 0.9804, 0.9844, 0.9914, 0.9995 /
data ((cldnuctab( 9,iw,iconc,7, 4),iw=1,9),iconc=1,7)/    &
0.9989, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9987, 0.9989, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9985, 0.9987, 0.9989, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9983, 0.9983, 0.9985, 0.9987, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9970, 0.9970, 0.9970, 0.9974, 0.9980, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.9924, 0.9924, 0.9924, 0.9924, 0.9933, 0.9948, 0.9974, 0.9999, 0.9999,  &
0.9804, 0.9804, 0.9804, 0.9804, 0.9804, 0.9825, 0.9861, 0.9914, 0.9990 /
data ((cldnuctab( 9,iw,iconc,1, 5),iw=1,9),iconc=1,7)/    &
0.9989, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9983, 0.9989, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9965, 0.9970, 0.9983, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9914, 0.9914, 0.9933, 0.9955, 0.9992, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9825, 0.9825, 0.9825, 0.9844, 0.9890, 0.9965, 0.9999, 0.9995, 0.9999,  &
0.9729, 0.9729, 0.9729, 0.9756, 0.9756, 0.9825, 0.9914, 0.9992, 0.9999,  &
0.9699, 0.9699, 0.9699, 0.9699, 0.9699, 0.9729, 0.9756, 0.9844, 0.9903 /
data ((cldnuctab( 9,iw,iconc,2, 5),iw=1,9),iconc=1,7)/    &
0.9990, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9987, 0.9990, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9980, 0.9983, 0.9989, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9955, 0.9955, 0.9965, 0.9977, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9890, 0.9890, 0.9890, 0.9914, 0.9933, 0.9983, 0.9999, 0.9999, 0.9999,  &
0.9782, 0.9804, 0.9804, 0.9804, 0.9825, 0.9861, 0.9941, 0.9999, 0.9999,  &
0.9729, 0.9729, 0.9729, 0.9729, 0.9729, 0.9756, 0.9804, 0.9876, 0.9970 /
data ((cldnuctab( 9,iw,iconc,3, 5),iw=1,9),iconc=1,7)/    &
0.9992, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9989, 0.9990, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9985, 0.9987, 0.9990, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9974, 0.9974, 0.9977, 0.9985, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9933, 0.9933, 0.9941, 0.9948, 0.9965, 0.9990, 0.9999, 0.9999, 0.9999,  &
0.9844, 0.9861, 0.9861, 0.9861, 0.9876, 0.9903, 0.9960, 0.9999, 0.9996,  &
0.9756, 0.9756, 0.9756, 0.9756, 0.9782, 0.9782, 0.9825, 0.9903, 0.9993 /
data ((cldnuctab( 9,iw,iconc,4, 5),iw=1,9),iconc=1,7)/    &
0.9993, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9990, 0.9992, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9989, 0.9990, 0.9992, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9983, 0.9983, 0.9985, 0.9989, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9960, 0.9960, 0.9960, 0.9965, 0.9974, 0.9992, 0.9999, 0.9999, 0.9999,  &
0.9890, 0.9890, 0.9903, 0.9903, 0.9914, 0.9933, 0.9970, 0.9999, 0.9999,  &
0.9804, 0.9804, 0.9804, 0.9804, 0.9804, 0.9825, 0.9844, 0.9914, 0.9995 /
data ((cldnuctab( 9,iw,iconc,5, 5),iw=1,9),iconc=1,7)/    &
0.9993, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9992, 0.9993, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9990, 0.9992, 0.9993, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9985, 0.9987, 0.9987, 0.9990, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9970, 0.9970, 0.9970, 0.9974, 0.9980, 0.9992, 0.9999, 0.9999, 0.9999,  &
0.9924, 0.9924, 0.9924, 0.9924, 0.9933, 0.9948, 0.9974, 0.9999, 0.9999,  &
0.9825, 0.9825, 0.9825, 0.9825, 0.9844, 0.9844, 0.9876, 0.9924, 0.9993 /
data ((cldnuctab( 9,iw,iconc,6, 5),iw=1,9),iconc=1,7)/    &
0.9994, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9993, 0.9994, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9992, 0.9993, 0.9994, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9989, 0.9989, 0.9990, 0.9992, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9974, 0.9977, 0.9977, 0.9980, 0.9983, 0.9993, 0.9999, 0.9999, 0.9999,  &
0.9933, 0.9933, 0.9933, 0.9933, 0.9941, 0.9955, 0.9974, 0.9999, 0.9999,  &
0.9844, 0.9844, 0.9844, 0.9844, 0.9844, 0.9861, 0.9876, 0.9924, 0.9989 /
data ((cldnuctab( 9,iw,iconc,7, 5),iw=1,9),iconc=1,7)/    &
0.9995, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9994, 0.9994, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9993, 0.9993, 0.9994, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9990, 0.9990, 0.9990, 0.9993, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9977, 0.9977, 0.9980, 0.9980, 0.9985, 0.9993, 0.9999, 0.9999, 0.9999,  &
0.9941, 0.9941, 0.9941, 0.9941, 0.9948, 0.9955, 0.9974, 0.9998, 0.9999,  &
0.9861, 0.9861, 0.9861, 0.9861, 0.9861, 0.9876, 0.9890, 0.9924, 0.9983 /
data ((cldnuctab( 9,iw,iconc,1, 6),iw=1,9),iconc=1,7)/    &
0.9993, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9989, 0.9992, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9974, 0.9977, 0.9985, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9933, 0.9933, 0.9941, 0.9960, 0.9990, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9861, 0.9861, 0.9861, 0.9876, 0.9914, 0.9965, 0.9999, 0.9995, 0.9999,  &
0.9804, 0.9804, 0.9804, 0.9804, 0.9825, 0.9861, 0.9924, 0.9989, 0.9999,  &
0.9782, 0.9782, 0.9782, 0.9782, 0.9782, 0.9804, 0.9825, 0.9876, 0.9914 /
data ((cldnuctab( 9,iw,iconc,2, 6),iw=1,9),iconc=1,7)/    &
0.9994, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9992, 0.9994, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9985, 0.9987, 0.9992, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9965, 0.9965, 0.9970, 0.9980, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9914, 0.9914, 0.9914, 0.9924, 0.9948, 0.9980, 0.9999, 0.9999, 0.9999,  &
0.9844, 0.9844, 0.9844, 0.9844, 0.9861, 0.9890, 0.9948, 0.9998, 0.9999,  &
0.9804, 0.9804, 0.9804, 0.9804, 0.9804, 0.9825, 0.9844, 0.9903, 0.9970 /
data ((cldnuctab( 9,iw,iconc,3, 6),iw=1,9),iconc=1,7)/    &
0.9995, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9994, 0.9995, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9990, 0.9992, 0.9994, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9980, 0.9980, 0.9983, 0.9989, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9948, 0.9948, 0.9948, 0.9955, 0.9970, 0.9989, 0.9999, 0.9999, 0.9999,  &
0.9890, 0.9890, 0.9890, 0.9890, 0.9903, 0.9924, 0.9960, 0.9999, 0.9998,  &
0.9825, 0.9825, 0.9825, 0.9825, 0.9825, 0.9844, 0.9861, 0.9924, 0.9989 /
data ((cldnuctab( 9,iw,iconc,4, 6),iw=1,9),iconc=1,7)/    &
0.9995, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9995, 0.9995, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9993, 0.9994, 0.9995, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9987, 0.9987, 0.9989, 0.9992, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9965, 0.9965, 0.9970, 0.9974, 0.9980, 0.9992, 0.9999, 0.9999, 0.9999,  &
0.9914, 0.9914, 0.9914, 0.9924, 0.9924, 0.9941, 0.9970, 0.9999, 0.9999,  &
0.9844, 0.9844, 0.9844, 0.9844, 0.9861, 0.9861, 0.9890, 0.9933, 0.9992 /
data ((cldnuctab( 9,iw,iconc,5, 6),iw=1,9),iconc=1,7)/    &
0.9996, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9995, 0.9995, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9994, 0.9995, 0.9995, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9990, 0.9990, 0.9992, 0.9994, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9977, 0.9977, 0.9977, 0.9980, 0.9985, 0.9993, 0.9999, 0.9999, 0.9999,  &
0.9933, 0.9933, 0.9933, 0.9941, 0.9941, 0.9955, 0.9974, 0.9999, 0.9999,  &
0.9876, 0.9876, 0.9876, 0.9876, 0.9876, 0.9876, 0.9903, 0.9933, 0.9989 /
data ((cldnuctab( 9,iw,iconc,6, 6),iw=1,9),iconc=1,7)/    &
0.9996, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9996, 0.9996, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9995, 0.9995, 0.9996, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9992, 0.9992, 0.9993, 0.9995, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9980, 0.9980, 0.9983, 0.9983, 0.9987, 0.9993, 0.9999, 0.9999, 0.9999,  &
0.9948, 0.9948, 0.9948, 0.9948, 0.9948, 0.9960, 0.9977, 0.9998, 0.9999,  &
0.9876, 0.9876, 0.9876, 0.9876, 0.9890, 0.9890, 0.9903, 0.9933, 0.9985 /
data ((cldnuctab( 9,iw,iconc,7, 6),iw=1,9),iconc=1,7)/    &
0.9996, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9996, 0.9996, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9995, 0.9995, 0.9996, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9993, 0.9993, 0.9994, 0.9995, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9983, 0.9983, 0.9983, 0.9985, 0.9989, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.9948, 0.9948, 0.9955, 0.9955, 0.9955, 0.9960, 0.9977, 0.9996, 0.9999,  &
0.9890, 0.9890, 0.9890, 0.9890, 0.9890, 0.9903, 0.9914, 0.9933, 0.9980 /
data ((cldnuctab( 9,iw,iconc,1, 7),iw=1,9),iconc=1,7)/    &
0.9995, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9992, 0.9994, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9977, 0.9980, 0.9987, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9941, 0.9941, 0.9948, 0.9965, 0.9989, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9890, 0.9890, 0.9890, 0.9903, 0.9924, 0.9965, 0.9999, 0.9999, 0.9999,  &
0.9844, 0.9844, 0.9861, 0.9861, 0.9861, 0.9890, 0.9941, 0.9987, 0.9998,  &
0.9825, 0.9825, 0.9825, 0.9844, 0.9844, 0.9844, 0.9876, 0.9903, 0.9933 /
data ((cldnuctab( 9,iw,iconc,2, 7),iw=1,9),iconc=1,7)/    &
0.9996, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9994, 0.9995, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9989, 0.9990, 0.9993, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9970, 0.9970, 0.9974, 0.9983, 0.9995, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9924, 0.9924, 0.9933, 0.9941, 0.9955, 0.9980, 0.9999, 0.9999, 0.9999,  &
0.9876, 0.9876, 0.9876, 0.9876, 0.9890, 0.9914, 0.9955, 0.9997, 0.9999,  &
0.9825, 0.9825, 0.9844, 0.9844, 0.9844, 0.9861, 0.9876, 0.9924, 0.9974 /
data ((cldnuctab( 9,iw,iconc,3, 7),iw=1,9),iconc=1,7)/    &
0.9996, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9995, 0.9996, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9993, 0.9994, 0.9995, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9983, 0.9985, 0.9987, 0.9990, 0.9996, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9955, 0.9955, 0.9960, 0.9960, 0.9974, 0.9989, 0.9999, 0.9999, 0.9999,  &
0.9903, 0.9903, 0.9903, 0.9914, 0.9914, 0.9933, 0.9965, 0.9998, 0.9999,  &
0.9861, 0.9861, 0.9861, 0.9861, 0.9861, 0.9876, 0.9890, 0.9933, 0.9987 /
data ((cldnuctab( 9,iw,iconc,4, 7),iw=1,9),iconc=1,7)/    &
0.9997, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9996, 0.9996, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9995, 0.9995, 0.9996, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9990, 0.9990, 0.9992, 0.9994, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9974, 0.9974, 0.9974, 0.9977, 0.9983, 0.9992, 0.9999, 0.9999, 0.9999,  &
0.9933, 0.9933, 0.9933, 0.9933, 0.9941, 0.9948, 0.9974, 0.9998, 0.9999,  &
0.9876, 0.9876, 0.9876, 0.9876, 0.9890, 0.9890, 0.9903, 0.9941, 0.9990 /
data ((cldnuctab( 9,iw,iconc,5, 7),iw=1,9),iconc=1,7)/    &
0.9997, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9996, 0.9997, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9996, 0.9996, 0.9996, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9993, 0.9993, 0.9994, 0.9995, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9980, 0.9980, 0.9980, 0.9983, 0.9987, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.9948, 0.9948, 0.9948, 0.9948, 0.9948, 0.9960, 0.9977, 0.9998, 0.9999,  &
0.9890, 0.9890, 0.9890, 0.9903, 0.9903, 0.9903, 0.9914, 0.9941, 0.9987 /
data ((cldnuctab( 9,iw,iconc,6, 7),iw=1,9),iconc=1,7)/    &
0.9997, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9997, 0.9997, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9996, 0.9996, 0.9997, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9994, 0.9994, 0.9995, 0.9995, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9983, 0.9985, 0.9985, 0.9985, 0.9989, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.9955, 0.9955, 0.9955, 0.9955, 0.9960, 0.9965, 0.9977, 0.9997, 0.9999,  &
0.9903, 0.9903, 0.9903, 0.9903, 0.9903, 0.9914, 0.9924, 0.9941, 0.9983 /
data ((cldnuctab( 9,iw,iconc,7, 7),iw=1,9),iconc=1,7)/    &
0.9998, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9997, 0.9997, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9996, 0.9997, 0.9997, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9995, 0.9995, 0.9995, 0.9996, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999,  &
0.9985, 0.9987, 0.9987, 0.9987, 0.9990, 0.9994, 0.9999, 0.9999, 0.9999,  &
0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9970, 0.9980, 0.9996, 0.9999,  &
0.9914, 0.9914, 0.9914, 0.9914, 0.9914, 0.9914, 0.9924, 0.9948, 0.9980 /

data ((supersat( 1,iw,iconc,1, 1),iw=1,9),iconc=1,7)/    &
1.0171, 1.0249, 1.0378, 1.0598, 1.1011, 1.1902, 1.4374, 2.7352, 5.0027,  &
1.0142, 1.0206, 1.0308, 1.0478, 1.0780, 1.1379, 1.2786, 1.7430, 4.0530,  &
1.0120, 1.0172, 1.0255, 1.0392, 1.0625, 1.1060, 1.1982, 1.4418, 2.5803,  &
1.0099, 1.0147, 1.0217, 1.0327, 1.0514, 1.0848, 1.1507, 1.3035, 1.7889,  &
1.0088, 1.0125, 1.0186, 1.0278, 1.0433, 1.0699, 1.1197, 1.2257, 1.5066,  &
1.0077, 1.0109, 1.0161, 1.0238, 1.0370, 1.0586, 1.0980, 1.1763, 1.3613,  &
1.0056, 1.0093, 1.0138, 1.0207, 1.0320, 1.0502, 1.0820, 1.1423, 1.2736 /
data ((supersat( 1,iw,iconc,2, 1),iw=1,9),iconc=1,7)/    &
1.0139, 1.0201, 1.0298, 1.0463, 1.0760, 1.1361, 1.2844, 1.8220, 3.7340,  &
1.0115, 1.0163, 1.0242, 1.0367, 1.0587, 1.0999, 1.1892, 1.4354, 2.6145,  &
1.0097, 1.0135, 1.0198, 1.0300, 1.0468, 1.0772, 1.1373, 1.2793, 1.7414,  &
1.0078, 1.0114, 1.0166, 1.0250, 1.0384, 1.0618, 1.1057, 1.1990, 1.4465,  &
1.0070, 1.0100, 1.0140, 1.0211, 1.0322, 1.0509, 1.0846, 1.1515, 1.3076,  &
1.0061, 1.0082, 1.0124, 1.0181, 1.0275, 1.0429, 1.0697, 1.1203, 1.2287,  &
1.0044, 1.0073, 1.0106, 1.0159, 1.0236, 1.0365, 1.0585, 1.0984, 1.1784 /
data ((supersat( 1,iw,iconc,3, 1),iw=1,9),iconc=1,7)/    &
1.0117, 1.0166, 1.0245, 1.0373, 1.0598, 1.1036, 1.2038, 1.5019, 2.7537,  &
1.0097, 1.0135, 1.0197, 1.0294, 1.0461, 1.0763, 1.1382, 1.2917, 1.8312,  &
1.0080, 1.0111, 1.0162, 1.0239, 1.0366, 1.0588, 1.1011, 1.1937, 1.4490,  &
1.0066, 1.0093, 1.0135, 1.0196, 1.0299, 1.0469, 1.0780, 1.1403, 1.2886,  &
1.0056, 1.0081, 1.0114, 1.0166, 1.0249, 1.0386, 1.0625, 1.1077, 1.2052,  &
1.0050, 1.0066, 1.0096, 1.0142, 1.0211, 1.0324, 1.0515, 1.0862, 1.1557,  &
1.0040, 1.0058, 1.0085, 1.0119, 1.0180, 1.0277, 1.0433, 1.0709, 1.1234 /
data ((supersat( 1,iw,iconc,4, 1),iw=1,9),iconc=1,7)/    &
1.0101, 1.0142, 1.0206, 1.0311, 1.0491, 1.0830, 1.1565, 1.3541, 2.1048,  &
1.0083, 1.0115, 1.0165, 1.0245, 1.0376, 1.0610, 1.1072, 1.2143, 1.5329,  &
1.0067, 1.0093, 1.0135, 1.0197, 1.0297, 1.0468, 1.0784, 1.1442, 1.3096,  &
1.0058, 1.0079, 1.0112, 1.0162, 1.0240, 1.0371, 1.0602, 1.1048, 1.2041,  &
1.0044, 1.0067, 1.0094, 1.0135, 1.0197, 1.0303, 1.0480, 1.0805, 1.1468,  &
1.0042, 1.0056, 1.0078, 1.0115, 1.0168, 1.0252, 1.0394, 1.0644, 1.1122,  &
1.0036, 1.0045, 1.0069, 1.0098, 1.0144, 1.0215, 1.0331, 1.0530, 1.0895 /
data ((supersat( 1,iw,iconc,5, 1),iw=1,9),iconc=1,7)/    &
1.0089, 1.0125, 1.0179, 1.0268, 1.0417, 1.0694, 1.1273, 1.2738, 1.7469,  &
1.0072, 1.0100, 1.0143, 1.0210, 1.0319, 1.0509, 1.0874, 1.1688, 1.3906,  &
1.0059, 1.0082, 1.0116, 1.0167, 1.0249, 1.0388, 1.0637, 1.1140, 1.2333,  &
1.0050, 1.0069, 1.0096, 1.0137, 1.0201, 1.0305, 1.0486, 1.0826, 1.1549,  &
1.0041, 1.0055, 1.0079, 1.0114, 1.0164, 1.0247, 1.0384, 1.0630, 1.1114,  &
1.0031, 1.0049, 1.0067, 1.0096, 1.0137, 1.0204, 1.0313, 1.0501, 1.0851,  &
1.0031, 1.0039, 1.0056, 1.0081, 1.0117, 1.0172, 1.0262, 1.0411, 1.0677 /
data ((supersat( 1,iw,iconc,6, 1),iw=1,9),iconc=1,7)/    &
1.0078, 1.0111, 1.0159, 1.0235, 1.0365, 1.0599, 1.1079, 1.2247, 1.5716,  &
1.0064, 1.0090, 1.0126, 1.0184, 1.0277, 1.0438, 1.0742, 1.1399, 1.3105,  &
1.0054, 1.0073, 1.0101, 1.0146, 1.0216, 1.0332, 1.0538, 1.0944, 1.1877,  &
1.0042, 1.0060, 1.0083, 1.0118, 1.0172, 1.0259, 1.0407, 1.0680, 1.1245,  &
1.0035, 1.0050, 1.0068, 1.0098, 1.0140, 1.0207, 1.0319, 1.0515, 1.0891,  &
1.0031, 1.0042, 1.0059, 1.0082, 1.0116, 1.0170, 1.0258, 1.0406, 1.0675,  &
1.0022, 1.0037, 1.0046, 1.0067, 1.0098, 1.0143, 1.0213, 1.0331, 1.0535 /
data ((supersat( 1,iw,iconc,7, 1),iw=1,9),iconc=1,7)/    &
1.0071, 1.0099, 1.0142, 1.0209, 1.0323, 1.0528, 1.0940, 1.1919, 1.4692,  &
1.0058, 1.0080, 1.0113, 1.0164, 1.0245, 1.0385, 1.0646, 1.1200, 1.2596,  &
1.0047, 1.0065, 1.0091, 1.0130, 1.0190, 1.0291, 1.0467, 1.0810, 1.1577,  &
1.0040, 1.0053, 1.0074, 1.0104, 1.0151, 1.0225, 1.0351, 1.0580, 1.1045,  &
1.0031, 1.0045, 1.0062, 1.0086, 1.0122, 1.0179, 1.0273, 1.0435, 1.0743,  &
1.0028, 1.0036, 1.0051, 1.0070, 1.0101, 1.0146, 1.0219, 1.0341, 1.0558,  &
1.0024, 1.0031, 1.0043, 1.0059, 1.0085, 1.0122, 1.0180, 1.0275, 1.0439 /
data ((supersat( 1,iw,iconc,1, 2),iw=1,9),iconc=1,7)/    &
1.0147, 1.0218, 1.0340, 1.0552, 1.0957, 1.1844, 1.4324, 2.7307, 5.0001,  &
1.0121, 1.0181, 1.0276, 1.0438, 1.0732, 1.1323, 1.2730, 1.7384, 4.0482,  &
1.0101, 1.0151, 1.0227, 1.0356, 1.0582, 1.1009, 1.1926, 1.4365, 2.5746,  &
1.0088, 1.0126, 1.0192, 1.0297, 1.0476, 1.0802, 1.1453, 1.2979, 1.7835,  &
1.0076, 1.0111, 1.0163, 1.0253, 1.0399, 1.0657, 1.1147, 1.2201, 1.5009,  &
1.0055, 1.0092, 1.0142, 1.0218, 1.0338, 1.0548, 1.0932, 1.1708, 1.3554,  &
1.0055, 1.0086, 1.0125, 1.0190, 1.0293, 1.0467, 1.0775, 1.1370, 1.2676 /
data ((supersat( 1,iw,iconc,2, 2),iw=1,9),iconc=1,7)/    &
1.0119, 1.0174, 1.0265, 1.0421, 1.0709, 1.1305, 1.2793, 1.8185, 3.7311,  &
1.0098, 1.0142, 1.0213, 1.0332, 1.0542, 1.0947, 1.1837, 1.4309, 2.6107,  &
1.0079, 1.0118, 1.0175, 1.0270, 1.0430, 1.0725, 1.1321, 1.2742, 1.7372,  &
1.0069, 1.0100, 1.0146, 1.0224, 1.0350, 1.0577, 1.1008, 1.1938, 1.4416,  &
1.0059, 1.0082, 1.0126, 1.0190, 1.0293, 1.0473, 1.0802, 1.1463, 1.3024,  &
1.0041, 1.0073, 1.0106, 1.0163, 1.0251, 1.0396, 1.0657, 1.1154, 1.2234,  &
1.0040, 1.0066, 1.0096, 1.0139, 1.0215, 1.0337, 1.0548, 1.0938, 1.1731 /
data ((supersat( 1,iw,iconc,3, 2),iw=1,9),iconc=1,7)/    &
1.0099, 1.0143, 1.0214, 1.0334, 1.0544, 1.0983, 1.1984, 1.4984, 2.7508,  &
1.0080, 1.0115, 1.0171, 1.0262, 1.0420, 1.0714, 1.1329, 1.2872, 1.8286,  &
1.0066, 1.0095, 1.0139, 1.0212, 1.0332, 1.0546, 1.0962, 1.1886, 1.4451,  &
1.0056, 1.0081, 1.0114, 1.0173, 1.0269, 1.0433, 1.0735, 1.1353, 1.2839,  &
1.0049, 1.0065, 1.0098, 1.0145, 1.0224, 1.0354, 1.0585, 1.1031, 1.2003,  &
1.0036, 1.0057, 1.0085, 1.0123, 1.0190, 1.0296, 1.0479, 1.0819, 1.1509,  &
1.0027, 1.0052, 1.0073, 1.0109, 1.0165, 1.0253, 1.0401, 1.0670, 1.1187 /
data ((supersat( 1,iw,iconc,4, 2),iw=1,9),iconc=1,7)/    &
1.0084, 1.0121, 1.0179, 1.0276, 1.0448, 1.0780, 1.1515, 1.3505, 2.1027,  &
1.0067, 1.0097, 1.0142, 1.0215, 1.0339, 1.0565, 1.1022, 1.2097, 1.5297,  &
1.0057, 1.0080, 1.0115, 1.0172, 1.0266, 1.0429, 1.0738, 1.1393, 1.3055,  &
1.0044, 1.0066, 1.0094, 1.0141, 1.0214, 1.0337, 1.0561, 1.1001, 1.1994,  &
1.0041, 1.0055, 1.0079, 1.0117, 1.0175, 1.0274, 1.0444, 1.0762, 1.1421,  &
1.0034, 1.0045, 1.0068, 1.0097, 1.0148, 1.0227, 1.0363, 1.0605, 1.1077,  &
1.0024, 1.0043, 1.0058, 1.0086, 1.0125, 1.0194, 1.0303, 1.0496, 1.0853 /
data ((supersat( 1,iw,iconc,5, 2),iw=1,9),iconc=1,7)/    &
1.0073, 1.0105, 1.0154, 1.0235, 1.0378, 1.0647, 1.1225, 1.2702, 1.7451,  &
1.0059, 1.0084, 1.0122, 1.0182, 1.0284, 1.0467, 1.0827, 1.1642, 1.3876,  &
1.0050, 1.0069, 1.0097, 1.0144, 1.0221, 1.0351, 1.0594, 1.1092, 1.2293,  &
1.0040, 1.0055, 1.0079, 1.0117, 1.0176, 1.0274, 1.0448, 1.0781, 1.1503,  &
1.0032, 1.0048, 1.0067, 1.0097, 1.0144, 1.0220, 1.0351, 1.0590, 1.1069,  &
1.0030, 1.0037, 1.0056, 1.0080, 1.0121, 1.0182, 1.0284, 1.0466, 1.0809,  &
1.0026, 1.0032, 1.0049, 1.0070, 1.0102, 1.0152, 1.0237, 1.0380, 1.0640 /
data ((supersat( 1,iw,iconc,6, 2),iw=1,9),iconc=1,7)/    &
1.0065, 1.0093, 1.0135, 1.0206, 1.0327, 1.0555, 1.1032, 1.2212, 1.5698,  &
1.0053, 1.0075, 1.0107, 1.0159, 1.0245, 1.0398, 1.0696, 1.1354, 1.3075,  &
1.0042, 1.0060, 1.0085, 1.0125, 1.0189, 1.0298, 1.0497, 1.0899, 1.1836,  &
1.0037, 1.0049, 1.0069, 1.0101, 1.0150, 1.0230, 1.0371, 1.0638, 1.1201,  &
1.0029, 1.0041, 1.0059, 1.0082, 1.0121, 1.0183, 1.0289, 1.0478, 1.0848,  &
1.0021, 1.0035, 1.0045, 1.0066, 1.0101, 1.0149, 1.0231, 1.0374, 1.0636,  &
1.0021, 1.0026, 1.0042, 1.0059, 1.0085, 1.0126, 1.0191, 1.0303, 1.0500 /
data ((supersat( 1,iw,iconc,7, 2),iw=1,9),iconc=1,7)/    &
1.0059, 1.0083, 1.0121, 1.0182, 1.0289, 1.0486, 1.0895, 1.1882, 1.4673,  &
1.0047, 1.0066, 1.0095, 1.0141, 1.0216, 1.0348, 1.0603, 1.1157, 1.2566,  &
1.0039, 1.0053, 1.0076, 1.0111, 1.0166, 1.0259, 1.0428, 1.0766, 1.1538,  &
1.0031, 1.0045, 1.0062, 1.0088, 1.0130, 1.0199, 1.0318, 1.0540, 1.1002,  &
1.0028, 1.0035, 1.0050, 1.0071, 1.0105, 1.0157, 1.0245, 1.0400, 1.0702,  &
1.0022, 1.0031, 1.0043, 1.0060, 1.0086, 1.0127, 1.0195, 1.0310, 1.0522,  &
1.0016, 1.0027, 1.0033, 1.0050, 1.0072, 1.0106, 1.0160, 1.0250, 1.0407 /
data ((supersat( 1,iw,iconc,1, 3),iw=1,9),iconc=1,7)/    &
1.0123, 1.0191, 1.0304, 1.0505, 1.0900, 1.1782, 1.4271, 2.7259, 4.9977,  &
1.0104, 1.0156, 1.0245, 1.0398, 1.0682, 1.1265, 1.2671, 1.7335, 4.0428,  &
1.0088, 1.0128, 1.0202, 1.0323, 1.0539, 1.0956, 1.1866, 1.4309, 2.5689,  &
1.0075, 1.0113, 1.0169, 1.0267, 1.0439, 1.0754, 1.1397, 1.2919, 1.7776,  &
1.0055, 1.0093, 1.0147, 1.0228, 1.0366, 1.0614, 1.1094, 1.2141, 1.4946,  &
1.0054, 1.0086, 1.0126, 1.0195, 1.0309, 1.0510, 1.0883, 1.1649, 1.3489,  &
1.0052, 1.0071, 1.0110, 1.0165, 1.0265, 1.0431, 1.0729, 1.1312, 1.2610 /
data ((supersat( 1,iw,iconc,2, 3),iw=1,9),iconc=1,7)/    &
1.0099, 1.0150, 1.0234, 1.0380, 1.0658, 1.1247, 1.2740, 1.8151, 3.7276,  &
1.0081, 1.0123, 1.0187, 1.0298, 1.0498, 1.0893, 1.1780, 1.4261, 2.6064,  &
1.0069, 1.0101, 1.0154, 1.0241, 1.0393, 1.0678, 1.1266, 1.2687, 1.7325,  &
1.0058, 1.0084, 1.0130, 1.0198, 1.0319, 1.0536, 1.0958, 1.1883, 1.4364,  &
1.0041, 1.0072, 1.0107, 1.0168, 1.0265, 1.0437, 1.0756, 1.1410, 1.2968,  &
1.0040, 1.0065, 1.0096, 1.0141, 1.0226, 1.0365, 1.0615, 1.1104, 1.2177,  &
1.0039, 1.0050, 1.0075, 1.0126, 1.0195, 1.0309, 1.0511, 1.0890, 1.1674 /
data ((supersat( 1,iw,iconc,3, 3),iw=1,9),iconc=1,7)/    &
1.0081, 1.0122, 1.0185, 1.0298, 1.0505, 1.0928, 1.1930, 1.4948, 2.7479,  &
1.0066, 1.0099, 1.0148, 1.0232, 1.0381, 1.0666, 1.1275, 1.2824, 1.8252,  &
1.0055, 1.0081, 1.0120, 1.0186, 1.0298, 1.0503, 1.0911, 1.1834, 1.4409,  &
1.0048, 1.0066, 1.0101, 1.0153, 1.0242, 1.0397, 1.0690, 1.1302, 1.2789,  &
1.0035, 1.0057, 1.0085, 1.0129, 1.0199, 1.0323, 1.0545, 1.0983, 1.1951,  &
1.0028, 1.0051, 1.0073, 1.0111, 1.0170, 1.0270, 1.0445, 1.0774, 1.1458,  &
1.0027, 1.0040, 1.0058, 1.0090, 1.0145, 1.0228, 1.0370, 1.0629, 1.1138 /
data ((supersat( 1,iw,iconc,4, 3),iw=1,9),iconc=1,7)/    &
1.0068, 1.0102, 1.0154, 1.0243, 1.0406, 1.0729, 1.1463, 1.3470, 2.1002,  &
1.0057, 1.0082, 1.0122, 1.0188, 1.0304, 1.0521, 1.0970, 1.2050, 1.5267,  &
1.0043, 1.0066, 1.0097, 1.0149, 1.0235, 1.0391, 1.0691, 1.1343, 1.3014,  &
1.0040, 1.0055, 1.0081, 1.0121, 1.0189, 1.0306, 1.0520, 1.0953, 1.1946,  &
1.0033, 1.0044, 1.0068, 1.0100, 1.0156, 1.0246, 1.0409, 1.0719, 1.1372,  &
1.0022, 1.0041, 1.0058, 1.0087, 1.0131, 1.0205, 1.0333, 1.0566, 1.1031,  &
1.0016, 1.0035, 1.0045, 1.0072, 1.0113, 1.0173, 1.0277, 1.0461, 1.0810 /
data ((supersat( 1,iw,iconc,5, 3),iw=1,9),iconc=1,7)/    &
1.0060, 1.0088, 1.0132, 1.0206, 1.0339, 1.0599, 1.1175, 1.2665, 1.7433,  &
1.0049, 1.0070, 1.0103, 1.0157, 1.0252, 1.0425, 1.0778, 1.1596, 1.3846,  &
1.0039, 1.0055, 1.0082, 1.0124, 1.0194, 1.0316, 1.0550, 1.1044, 1.2252,  &
1.0031, 1.0048, 1.0068, 1.0100, 1.0154, 1.0245, 1.0410, 1.0736, 1.1457,  &
1.0029, 1.0037, 1.0056, 1.0080, 1.0125, 1.0196, 1.0319, 1.0550, 1.1024,  &
1.0024, 1.0032, 1.0049, 1.0071, 1.0104, 1.0162, 1.0258, 1.0432, 1.0766,  &
1.0017, 1.0030, 1.0037, 1.0057, 1.0088, 1.0136, 1.0215, 1.0350, 1.0601 /
data ((supersat( 1,iw,iconc,6, 3),iw=1,9),iconc=1,7)/    &
1.0054, 1.0077, 1.0115, 1.0178, 1.0291, 1.0510, 1.0985, 1.2176, 1.5680,  &
1.0042, 1.0061, 1.0090, 1.0136, 1.0215, 1.0360, 1.0651, 1.1309, 1.3045,  &
1.0036, 1.0049, 1.0072, 1.0107, 1.0164, 1.0265, 1.0456, 1.0853, 1.1796,  &
1.0027, 1.0041, 1.0059, 1.0084, 1.0129, 1.0204, 1.0337, 1.0595, 1.1156,  &
1.0021, 1.0034, 1.0046, 1.0069, 1.0104, 1.0162, 1.0260, 1.0441, 1.0805,  &
1.0020, 1.0024, 1.0041, 1.0060, 1.0086, 1.0131, 1.0208, 1.0342, 1.0597,  &
1.0019, 1.0021, 1.0035, 1.0046, 1.0072, 1.0109, 1.0171, 1.0276, 1.0466 /
data ((supersat( 1,iw,iconc,7, 3),iw=1,9),iconc=1,7)/    &
1.0047, 1.0068, 1.0102, 1.0157, 1.0256, 1.0445, 1.0850, 1.1847, 1.4655,  &
1.0038, 1.0054, 1.0080, 1.0119, 1.0188, 1.0313, 1.0560, 1.1114, 1.2537,  &
1.0031, 1.0045, 1.0064, 1.0094, 1.0143, 1.0229, 1.0390, 1.0722, 1.1499,  &
1.0027, 1.0034, 1.0050, 1.0074, 1.0112, 1.0174, 1.0286, 1.0500, 1.0959,  &
1.0021, 1.0030, 1.0042, 1.0061, 1.0089, 1.0137, 1.0218, 1.0366, 1.0661,  &
1.0013, 1.0026, 1.0032, 1.0051, 1.0072, 1.0110, 1.0173, 1.0282, 1.0485,  &
1.0011, 1.0019, 1.0031, 1.0043, 1.0062, 1.0092, 1.0142, 1.0225, 1.0376 /
data ((supersat( 1,iw,iconc,1, 4),iw=1,9),iconc=1,7)/    &
1.0109, 1.0168, 1.0272, 1.0462, 1.0846, 1.1723, 1.4222, 2.7212, 4.9952,  &
1.0088, 1.0137, 1.0218, 1.0362, 1.0635, 1.1209, 1.2614, 1.7287, 4.0375,  &
1.0076, 1.0116, 1.0181, 1.0292, 1.0499, 1.0905, 1.1808, 1.4254, 2.5629,  &
1.0055, 1.0092, 1.0153, 1.0243, 1.0404, 1.0709, 1.1342, 1.2861, 1.7716,  &
1.0054, 1.0086, 1.0127, 1.0204, 1.0334, 1.0573, 1.1042, 1.2081, 1.4883,  &
1.0051, 1.0072, 1.0113, 1.0176, 1.0283, 1.0472, 1.0835, 1.1590, 1.3423,  &
1.0043, 1.0055, 1.0091, 1.0152, 1.0239, 1.0396, 1.0683, 1.1255, 1.2542 /
data ((supersat( 1,iw,iconc,2, 4),iw=1,9),iconc=1,7)/    &
1.0085, 1.0128, 1.0206, 1.0343, 1.0610, 1.1192, 1.2690, 1.8118, 3.7245,  &
1.0069, 1.0104, 1.0163, 1.0267, 1.0458, 1.0843, 1.1726, 1.4216, 2.6020,  &
1.0058, 1.0089, 1.0134, 1.0216, 1.0359, 1.0634, 1.1214, 1.2636, 1.7284,  &
1.0041, 1.0069, 1.0111, 1.0178, 1.0290, 1.0498, 1.0910, 1.1829, 1.4314,  &
1.0040, 1.0065, 1.0098, 1.0151, 1.0242, 1.0404, 1.0713, 1.1358, 1.2913,  &
1.0038, 1.0055, 1.0077, 1.0130, 1.0202, 1.0335, 1.0576, 1.1054, 1.2120,  &
1.0034, 1.0040, 1.0072, 1.0107, 1.0171, 1.0283, 1.0475, 1.0843, 1.1618 /
data ((supersat( 1,iw,iconc,3, 4),iw=1,9),iconc=1,7)/    &
1.0069, 1.0104, 1.0161, 1.0265, 1.0462, 1.0877, 1.1882, 1.4915, 2.7449,  &
1.0055, 1.0083, 1.0129, 1.0206, 1.0345, 1.0620, 1.1224, 1.2780, 1.8225,  &
1.0048, 1.0068, 1.0105, 1.0164, 1.0269, 1.0465, 1.0864, 1.1785, 1.4369,  &
1.0034, 1.0057, 1.0086, 1.0136, 1.0217, 1.0364, 1.0648, 1.1253, 1.2742,  &
1.0028, 1.0051, 1.0075, 1.0114, 1.0180, 1.0296, 1.0509, 1.0937, 1.1901,  &
1.0027, 1.0039, 1.0058, 1.0095, 1.0151, 1.0246, 1.0412, 1.0733, 1.1409,  &
1.0027, 1.0028, 1.0056, 1.0084, 1.0130, 1.0206, 1.0341, 1.0590, 1.1090 /
data ((supersat( 1,iw,iconc,4, 4),iw=1,9),iconc=1,7)/    &
1.0058, 1.0086, 1.0132, 1.0214, 1.0367, 1.0682, 1.1416, 1.3436, 2.0980,  &
1.0043, 1.0067, 1.0105, 1.0164, 1.0272, 1.0479, 1.0923, 1.2007, 1.5239,  &
1.0039, 1.0056, 1.0085, 1.0131, 1.0210, 1.0357, 1.0648, 1.1296, 1.2976,  &
1.0032, 1.0044, 1.0069, 1.0107, 1.0168, 1.0277, 1.0483, 1.0908, 1.1901,  &
1.0021, 1.0041, 1.0059, 1.0090, 1.0139, 1.0223, 1.0378, 1.0678, 1.1327,  &
1.0016, 1.0034, 1.0045, 1.0072, 1.0117, 1.0186, 1.0306, 1.0531, 1.0987,  &
1.0016, 1.0023, 1.0043, 1.0066, 1.0098, 1.0155, 1.0252, 1.0428, 1.0768 /
data ((supersat( 1,iw,iconc,5, 4),iw=1,9),iconc=1,7)/    &
1.0050, 1.0073, 1.0112, 1.0179, 1.0303, 1.0556, 1.1130, 1.2632, 1.7415,  &
1.0039, 1.0057, 1.0088, 1.0136, 1.0223, 1.0387, 1.0733, 1.1554, 1.3818,  &
1.0031, 1.0048, 1.0071, 1.0107, 1.0171, 1.0285, 1.0510, 1.1000, 1.2215,  &
1.0028, 1.0036, 1.0056, 1.0087, 1.0136, 1.0220, 1.0376, 1.0694, 1.1414,  &
1.0023, 1.0032, 1.0049, 1.0073, 1.0110, 1.0175, 1.0292, 1.0514, 1.0981,  &
1.0015, 1.0029, 1.0037, 1.0057, 1.0093, 1.0145, 1.0235, 1.0400, 1.0727,  &
1.0009, 1.0025, 1.0032, 1.0053, 1.0079, 1.0122, 1.0194, 1.0323, 1.0566 /
data ((supersat( 1,iw,iconc,6, 4),iw=1,9),iconc=1,7)/    &
1.0042, 1.0064, 1.0097, 1.0154, 1.0259, 1.0470, 1.0942, 1.2144, 1.5665,  &
1.0036, 1.0051, 1.0076, 1.0116, 1.0189, 1.0325, 1.0609, 1.1269, 1.3019,  &
1.0027, 1.0041, 1.0060, 1.0091, 1.0143, 1.0237, 1.0419, 1.0811, 1.1760,  &
1.0021, 1.0034, 1.0047, 1.0073, 1.0113, 1.0180, 1.0306, 1.0557, 1.1115,  &
1.0020, 1.0024, 1.0041, 1.0061, 1.0091, 1.0143, 1.0235, 1.0407, 1.0765,  &
1.0018, 1.0021, 1.0034, 1.0048, 1.0076, 1.0117, 1.0187, 1.0314, 1.0561,  &
1.0015, 1.0021, 1.0024, 1.0043, 1.0064, 1.0098, 1.0153, 1.0252, 1.0434 /
data ((supersat( 1,iw,iconc,7, 4),iw=1,9),iconc=1,7)/    &
1.0039, 1.0057, 1.0086, 1.0136, 1.0226, 1.0407, 1.0810, 1.1816, 1.4640,  &
1.0030, 1.0045, 1.0066, 1.0102, 1.0164, 1.0281, 1.0521, 1.1075, 1.2512,  &
1.0026, 1.0035, 1.0052, 1.0080, 1.0124, 1.0203, 1.0356, 1.0683, 1.1464,  &
1.0020, 1.0030, 1.0043, 1.0064, 1.0097, 1.0153, 1.0257, 1.0464, 1.0920,  &
1.0013, 1.0025, 1.0033, 1.0051, 1.0077, 1.0120, 1.0195, 1.0336, 1.0624,  &
1.0011, 1.0018, 1.0030, 1.0043, 1.0065, 1.0098, 1.0155, 1.0257, 1.0452,  &
1.0011, 1.0011, 1.0026, 1.0033, 1.0052, 1.0081, 1.0126, 1.0204, 1.0347 /
data ((supersat( 1,iw,iconc,1, 5),iw=1,9),iconc=1,7)/    &
1.0098, 1.0154, 1.0254, 1.0440, 1.0817, 1.1692, 1.4196, 2.7187, 4.9938,  &
1.0083, 1.0124, 1.0205, 1.0343, 1.0610, 1.1179, 1.2585, 1.7263, 4.0342,  &
1.0066, 1.0107, 1.0168, 1.0278, 1.0477, 1.0878, 1.1778, 1.4225, 2.5599,  &
1.0055, 1.0090, 1.0143, 1.0228, 1.0386, 1.0684, 1.1312, 1.2829, 1.7684,  &
1.0052, 1.0080, 1.0122, 1.0193, 1.0320, 1.0551, 1.1014, 1.2049, 1.4848,  &
1.0047, 1.0061, 1.0102, 1.0162, 1.0266, 1.0454, 1.0808, 1.1558, 1.3386,  &
1.0036, 1.0054, 1.0088, 1.0140, 1.0226, 1.0378, 1.0658, 1.1222, 1.2503 /
data ((supersat( 1,iw,iconc,2, 5),iw=1,9),iconc=1,7)/    &
1.0074, 1.0120, 1.0191, 1.0323, 1.0585, 1.1163, 1.2665, 1.8102, 3.7233,  &
1.0065, 1.0098, 1.0154, 1.0251, 1.0437, 1.0816, 1.1698, 1.4193, 2.6000,  &
1.0050, 1.0077, 1.0127, 1.0202, 1.0341, 1.0611, 1.1186, 1.2609, 1.7262,  &
1.0041, 1.0069, 1.0104, 1.0166, 1.0277, 1.0478, 1.0884, 1.1801, 1.4287,  &
1.0039, 1.0058, 1.0090, 1.0138, 1.0228, 1.0386, 1.0690, 1.1330, 1.2884,  &
1.0036, 1.0041, 1.0074, 1.0121, 1.0193, 1.0319, 1.0555, 1.1028, 1.2090,  &
1.0029, 1.0040, 1.0068, 1.0103, 1.0164, 1.0267, 1.0455, 1.0818, 1.1587 /
data ((supersat( 1,iw,iconc,3, 5),iw=1,9),iconc=1,7)/    &
1.0060, 1.0095, 1.0150, 1.0249, 1.0440, 1.0851, 1.1857, 1.4899, 2.7433,  &
1.0052, 1.0077, 1.0118, 1.0191, 1.0327, 1.0596, 1.1198, 1.2758, 1.8210,  &
1.0041, 1.0058, 1.0097, 1.0154, 1.0254, 1.0444, 1.0839, 1.1760, 1.4348,  &
1.0028, 1.0054, 1.0082, 1.0127, 1.0205, 1.0347, 1.0626, 1.1227, 1.2719,  &
1.0027, 1.0045, 1.0066, 1.0107, 1.0169, 1.0280, 1.0489, 1.0913, 1.1876,  &
1.0027, 1.0031, 1.0057, 1.0088, 1.0142, 1.0232, 1.0395, 1.0710, 1.1382,  &
1.0025, 1.0028, 1.0053, 1.0079, 1.0118, 1.0196, 1.0326, 1.0570, 1.1064 /
data ((supersat( 1,iw,iconc,4, 5),iw=1,9),iconc=1,7)/    &
1.0051, 1.0078, 1.0122, 1.0199, 1.0347, 1.0658, 1.1392, 1.3420, 2.0970,  &
1.0042, 1.0063, 1.0094, 1.0153, 1.0256, 1.0458, 1.0898, 1.1985, 1.5227,  &
1.0036, 1.0048, 1.0077, 1.0120, 1.0197, 1.0338, 1.0625, 1.1272, 1.2957,  &
1.0026, 1.0043, 1.0065, 1.0096, 1.0158, 1.0262, 1.0463, 1.0885, 1.1879,  &
1.0016, 1.0038, 1.0051, 1.0083, 1.0129, 1.0212, 1.0361, 1.0657, 1.1303,  &
1.0016, 1.0028, 1.0044, 1.0070, 1.0109, 1.0173, 1.0291, 1.0511, 1.0964,  &
1.0016, 1.0018, 1.0041, 1.0061, 1.0093, 1.0146, 1.0241, 1.0412, 1.0747 /
data ((supersat( 1,iw,iconc,5, 5),iw=1,9),iconc=1,7)/    &
1.0045, 1.0067, 1.0103, 1.0166, 1.0285, 1.0533, 1.1107, 1.2617, 1.7409,  &
1.0033, 1.0052, 1.0079, 1.0126, 1.0208, 1.0368, 1.0711, 1.1534, 1.3806,  &
1.0030, 1.0043, 1.0064, 1.0098, 1.0159, 1.0269, 1.0489, 1.0977, 1.2196,  &
1.0025, 1.0032, 1.0053, 1.0079, 1.0125, 1.0206, 1.0359, 1.0673, 1.1393,  &
1.0018, 1.0030, 1.0043, 1.0066, 1.0102, 1.0165, 1.0277, 1.0495, 1.0959,  &
1.0011, 1.0027, 1.0032, 1.0056, 1.0084, 1.0136, 1.0223, 1.0383, 1.0706,  &
1.0007, 1.0020, 1.0031, 1.0049, 1.0074, 1.0114, 1.0184, 1.0308, 1.0546 /
data ((supersat( 1,iw,iconc,6, 5),iw=1,9),iconc=1,7)/    &
1.0039, 1.0058, 1.0089, 1.0142, 1.0243, 1.0449, 1.0921, 1.2129, 1.5658,  &
1.0031, 1.0044, 1.0069, 1.0107, 1.0176, 1.0308, 1.0587, 1.1249, 1.3007,  &
1.0022, 1.0038, 1.0055, 1.0082, 1.0133, 1.0223, 1.0400, 1.0790, 1.1742,  &
1.0020, 1.0029, 1.0043, 1.0065, 1.0103, 1.0169, 1.0290, 1.0537, 1.1095,  &
1.0019, 1.0021, 1.0038, 1.0056, 1.0084, 1.0134, 1.0221, 1.0390, 1.0744,  &
1.0016, 1.0021, 1.0029, 1.0044, 1.0068, 1.0108, 1.0177, 1.0300, 1.0543,  &
1.0011, 1.0020, 1.0021, 1.0041, 1.0060, 1.0089, 1.0145, 1.0240, 1.0418 /
data ((supersat( 1,iw,iconc,7, 5),iw=1,9),iconc=1,7)/    &
1.0034, 1.0050, 1.0078, 1.0125, 1.0211, 1.0388, 1.0790, 1.1802, 1.4633,  &
1.0029, 1.0041, 1.0061, 1.0094, 1.0152, 1.0265, 1.0501, 1.1057, 1.2500,  &
1.0023, 1.0031, 1.0048, 1.0072, 1.0114, 1.0190, 1.0339, 1.0663, 1.1447,  &
1.0016, 1.0027, 1.0038, 1.0057, 1.0088, 1.0143, 1.0243, 1.0446, 1.0901,  &
1.0011, 1.0021, 1.0031, 1.0048, 1.0070, 1.0112, 1.0184, 1.0320, 1.0605,  &
1.0011, 1.0013, 1.0028, 1.0037, 1.0058, 1.0090, 1.0145, 1.0244, 1.0436,  &
1.0011, 1.0011, 1.0021, 1.0032, 1.0051, 1.0073, 1.0119, 1.0194, 1.0333 /
data ((supersat( 1,iw,iconc,1, 6),iw=1,9),iconc=1,7)/    &
1.0089, 1.0148, 1.0244, 1.0425, 1.0799, 1.1672, 1.4180, 2.7173, 4.9930,  &
1.0079, 1.0121, 1.0195, 1.0331, 1.0594, 1.1160, 1.2566, 1.7249, 4.0324,  &
1.0059, 1.0099, 1.0159, 1.0267, 1.0464, 1.0860, 1.1759, 1.4206, 2.5576,  &
1.0054, 1.0087, 1.0134, 1.0221, 1.0374, 1.0669, 1.1293, 1.2809, 1.7663,  &
1.0050, 1.0074, 1.0117, 1.0187, 1.0309, 1.0536, 1.0996, 1.2028, 1.4823,  &
1.0043, 1.0055, 1.0093, 1.0157, 1.0258, 1.0440, 1.0791, 1.1536, 1.3361,  &
1.0031, 1.0053, 1.0085, 1.0131, 1.0218, 1.0363, 1.0641, 1.1201, 1.2477 /
data ((supersat( 1,iw,iconc,2, 6),iw=1,9),iconc=1,7)/    &
1.0070, 1.0113, 1.0182, 1.0311, 1.0569, 1.1145, 1.2650, 1.8093, 3.7223,  &
1.0061, 1.0093, 1.0146, 1.0241, 1.0423, 1.0800, 1.1680, 1.4179, 2.5985,  &
1.0043, 1.0073, 1.0120, 1.0193, 1.0330, 1.0595, 1.1169, 1.2592, 1.7248,  &
1.0040, 1.0066, 1.0101, 1.0161, 1.0266, 1.0464, 1.0868, 1.1783, 1.4269,  &
1.0038, 1.0052, 1.0083, 1.0134, 1.0221, 1.0374, 1.0674, 1.1312, 1.2866,  &
1.0033, 1.0040, 1.0072, 1.0113, 1.0186, 1.0309, 1.0541, 1.1010, 1.2070,  &
1.0026, 1.0039, 1.0065, 1.0099, 1.0157, 1.0257, 1.0441, 1.0800, 1.1566 /
data ((supersat( 1,iw,iconc,3, 6),iw=1,9),iconc=1,7)/    &
1.0056, 1.0088, 1.0142, 1.0238, 1.0426, 1.0834, 1.1842, 1.4889, 2.7425,  &
1.0049, 1.0072, 1.0111, 1.0183, 1.0315, 1.0581, 1.1181, 1.2744, 1.8204,  &
1.0036, 1.0057, 1.0090, 1.0146, 1.0244, 1.0432, 1.0823, 1.1744, 1.4335,  &
1.0028, 1.0052, 1.0078, 1.0119, 1.0195, 1.0336, 1.0612, 1.1211, 1.2703,  &
1.0027, 1.0040, 1.0059, 1.0102, 1.0163, 1.0271, 1.0477, 1.0897, 1.1859,  &
1.0026, 1.0028, 1.0056, 1.0085, 1.0138, 1.0224, 1.0384, 1.0696, 1.1365,  &
1.0023, 1.0027, 1.0049, 1.0074, 1.0114, 1.0188, 1.0315, 1.0555, 1.1047 /
data ((supersat( 1,iw,iconc,4, 6),iw=1,9),iconc=1,7)/    &
1.0046, 1.0072, 1.0114, 1.0190, 1.0335, 1.0642, 1.1377, 1.3410, 2.0967,  &
1.0040, 1.0059, 1.0090, 1.0145, 1.0245, 1.0444, 1.0883, 1.1972, 1.5218,  &
1.0032, 1.0044, 1.0070, 1.0116, 1.0188, 1.0327, 1.0611, 1.1257, 1.2945,  &
1.0021, 1.0041, 1.0061, 1.0094, 1.0150, 1.0253, 1.0451, 1.0870, 1.1865,  &
1.0016, 1.0034, 1.0045, 1.0076, 1.0122, 1.0203, 1.0350, 1.0643, 1.1288,  &
1.0016, 1.0024, 1.0043, 1.0067, 1.0103, 1.0168, 1.0282, 1.0499, 1.0949,  &
1.0016, 1.0016, 1.0039, 1.0056, 1.0089, 1.0141, 1.0232, 1.0400, 1.0732 /
data ((supersat( 1,iw,iconc,5, 6),iw=1,9),iconc=1,7)/    &
1.0041, 1.0062, 1.0095, 1.0157, 1.0274, 1.0519, 1.1093, 1.2607, 1.7405,  &
1.0031, 1.0050, 1.0074, 1.0118, 1.0199, 1.0356, 1.0696, 1.1521, 1.3798,  &
1.0028, 1.0038, 1.0058, 1.0094, 1.0152, 1.0259, 1.0477, 1.0963, 1.2185,  &
1.0022, 1.0032, 1.0050, 1.0076, 1.0120, 1.0199, 1.0347, 1.0659, 1.1380,  &
1.0015, 1.0029, 1.0039, 1.0060, 1.0098, 1.0159, 1.0268, 1.0483, 1.0945,  &
1.0009, 1.0024, 1.0032, 1.0054, 1.0080, 1.0129, 1.0214, 1.0373, 1.0693,  &
1.0006, 1.0017, 1.0031, 1.0045, 1.0070, 1.0107, 1.0176, 1.0299, 1.0534 /
data ((supersat( 1,iw,iconc,6, 6),iw=1,9),iconc=1,7)/    &
1.0037, 1.0054, 1.0083, 1.0134, 1.0232, 1.0436, 1.0909, 1.2120, 1.5654,  &
1.0028, 1.0041, 1.0063, 1.0100, 1.0167, 1.0296, 1.0574, 1.1237, 1.2999,  &
1.0021, 1.0035, 1.0051, 1.0079, 1.0126, 1.0213, 1.0388, 1.0776, 1.1731,  &
1.0020, 1.0025, 1.0041, 1.0063, 1.0099, 1.0161, 1.0280, 1.0524, 1.1082,  &
1.0018, 1.0021, 1.0035, 1.0051, 1.0081, 1.0127, 1.0213, 1.0379, 1.0732,  &
1.0014, 1.0020, 1.0025, 1.0043, 1.0065, 1.0104, 1.0169, 1.0290, 1.0531,  &
1.0009, 1.0019, 1.0021, 1.0038, 1.0057, 1.0086, 1.0139, 1.0231, 1.0407 /
data ((supersat( 1,iw,iconc,7, 6),iw=1,9),iconc=1,7)/    &
1.0030, 1.0047, 1.0073, 1.0118, 1.0201, 1.0376, 1.0777, 1.1793, 1.4630,  &
1.0026, 1.0037, 1.0056, 1.0088, 1.0144, 1.0254, 1.0489, 1.1045, 1.2492,  &
1.0020, 1.0030, 1.0045, 1.0067, 1.0108, 1.0181, 1.0328, 1.0650, 1.1437,  &
1.0013, 1.0025, 1.0034, 1.0053, 1.0084, 1.0136, 1.0234, 1.0434, 1.0889,  &
1.0011, 1.0018, 1.0030, 1.0045, 1.0068, 1.0106, 1.0176, 1.0310, 1.0593,  &
1.0011, 1.0011, 1.0026, 1.0033, 1.0053, 1.0087, 1.0139, 1.0235, 1.0425,  &
1.0010, 1.0011, 1.0018, 1.0031, 1.0048, 1.0071, 1.0113, 1.0186, 1.0323 /
data ((supersat( 1,iw,iconc,1, 7),iw=1,9),iconc=1,7)/    &
1.0088, 1.0143, 1.0237, 1.0415, 1.0786, 1.1659, 1.4168, 2.7165, 4.9924,  &
1.0075, 1.0117, 1.0188, 1.0322, 1.0583, 1.1146, 1.2553, 1.7239, 4.0311,  &
1.0055, 1.0092, 1.0156, 1.0258, 1.0453, 1.0848, 1.1745, 1.4193, 2.5563,  &
1.0053, 1.0085, 1.0127, 1.0215, 1.0364, 1.0657, 1.1279, 1.2795, 1.7648,  &
1.0049, 1.0069, 1.0113, 1.0181, 1.0299, 1.0525, 1.0982, 1.2012, 1.4806,  &
1.0039, 1.0054, 1.0090, 1.0153, 1.0251, 1.0429, 1.0778, 1.1520, 1.3343,  &
1.0028, 1.0052, 1.0082, 1.0124, 1.0211, 1.0355, 1.0628, 1.1184, 1.2457 /
data ((supersat( 1,iw,iconc,2, 7),iw=1,9),iconc=1,7)/    &
1.0069, 1.0107, 1.0176, 1.0302, 1.0557, 1.1132, 1.2639, 1.8086, 3.7216,  &
1.0057, 1.0088, 1.0140, 1.0235, 1.0413, 1.0787, 1.1668, 1.4169, 2.5977,  &
1.0041, 1.0072, 1.0115, 1.0188, 1.0322, 1.0584, 1.1156, 1.2581, 1.7238,  &
1.0039, 1.0063, 1.0098, 1.0156, 1.0258, 1.0455, 1.0856, 1.1771, 1.4257,  &
1.0037, 1.0047, 1.0077, 1.0131, 1.0215, 1.0367, 1.0663, 1.1299, 1.2852,  &
1.0031, 1.0040, 1.0070, 1.0106, 1.0180, 1.0301, 1.0530, 1.0997, 1.2055,  &
1.0023, 1.0038, 1.0061, 1.0095, 1.0152, 1.0251, 1.0432, 1.0788, 1.1550 /
data ((supersat( 1,iw,iconc,3, 7),iw=1,9),iconc=1,7)/    &
1.0054, 1.0083, 1.0135, 1.0229, 1.0415, 1.0822, 1.1831, 1.4884, 2.7420,  &
1.0046, 1.0068, 1.0108, 1.0177, 1.0306, 1.0571, 1.1170, 1.2735, 1.8195,  &
1.0031, 1.0056, 1.0086, 1.0140, 1.0237, 1.0422, 1.0812, 1.1733, 1.4328,  &
1.0028, 1.0049, 1.0074, 1.0114, 1.0191, 1.0328, 1.0601, 1.1200, 1.2693,  &
1.0027, 1.0036, 1.0058, 1.0097, 1.0158, 1.0265, 1.0467, 1.0886, 1.1847,  &
1.0025, 1.0028, 1.0054, 1.0083, 1.0133, 1.0218, 1.0375, 1.0685, 1.1352,  &
1.0022, 1.0027, 1.0046, 1.0070, 1.0112, 1.0182, 1.0307, 1.0546, 1.1034 /
data ((supersat( 1,iw,iconc,4, 7),iw=1,9),iconc=1,7)/    &
1.0043, 1.0067, 1.0110, 1.0183, 1.0325, 1.0631, 1.1367, 1.3403, 2.0962,  &
1.0038, 1.0056, 1.0087, 1.0138, 1.0238, 1.0435, 1.0872, 1.1963, 1.5213,  &
1.0029, 1.0043, 1.0069, 1.0111, 1.0182, 1.0319, 1.0601, 1.1247, 1.2937,  &
1.0018, 1.0040, 1.0058, 1.0091, 1.0144, 1.0246, 1.0442, 1.0860, 1.1855,  &
1.0016, 1.0031, 1.0044, 1.0071, 1.0119, 1.0197, 1.0342, 1.0633, 1.1277,  &
1.0016, 1.0020, 1.0042, 1.0065, 1.0098, 1.0163, 1.0275, 1.0490, 1.0938,  &
1.0016, 1.0016, 1.0037, 1.0052, 1.0086, 1.0137, 1.0225, 1.0391, 1.0721 /
data ((supersat( 1,iw,iconc,5, 7),iw=1,9),iconc=1,7)/    &
1.0038, 1.0058, 1.0092, 1.0151, 1.0266, 1.0509, 1.1083, 1.2601, 1.7402,  &
1.0030, 1.0047, 1.0072, 1.0114, 1.0192, 1.0347, 1.0686, 1.1512, 1.3793,  &
1.0027, 1.0035, 1.0055, 1.0090, 1.0146, 1.0252, 1.0467, 1.0953, 1.2178,  &
1.0020, 1.0031, 1.0048, 1.0073, 1.0116, 1.0192, 1.0340, 1.0650, 1.1370,  &
1.0013, 1.0028, 1.0035, 1.0057, 1.0095, 1.0153, 1.0261, 1.0474, 1.0935,  &
1.0007, 1.0022, 1.0032, 1.0052, 1.0078, 1.0125, 1.0209, 1.0365, 1.0683,  &
1.0006, 1.0015, 1.0030, 1.0042, 1.0066, 1.0103, 1.0171, 1.0292, 1.0525 /
data ((supersat( 1,iw,iconc,6, 7),iw=1,9),iconc=1,7)/    &
1.0035, 1.0051, 1.0078, 1.0129, 1.0225, 1.0426, 1.0898, 1.2114, 1.5651,  &
1.0025, 1.0040, 1.0060, 1.0096, 1.0161, 1.0288, 1.0564, 1.1229, 1.2994,  &
1.0021, 1.0032, 1.0048, 1.0076, 1.0121, 1.0206, 1.0380, 1.0767, 1.1724,  &
1.0019, 1.0022, 1.0040, 1.0061, 1.0095, 1.0156, 1.0273, 1.0516, 1.1074,  &
1.0016, 1.0021, 1.0032, 1.0047, 1.0077, 1.0123, 1.0207, 1.0372, 1.0723,  &
1.0012, 1.0020, 1.0022, 1.0042, 1.0064, 1.0101, 1.0164, 1.0283, 1.0522,  &
1.0008, 1.0018, 1.0021, 1.0036, 1.0053, 1.0084, 1.0134, 1.0225, 1.0399 /
data ((supersat( 1,iw,iconc,7, 7),iw=1,9),iconc=1,7)/    &
1.0030, 1.0045, 1.0069, 1.0112, 1.0195, 1.0368, 1.0769, 1.1787, 1.4627,  &
1.0025, 1.0034, 1.0053, 1.0083, 1.0139, 1.0247, 1.0480, 1.1038, 1.2487,  &
1.0018, 1.0029, 1.0043, 1.0065, 1.0103, 1.0175, 1.0320, 1.0642, 1.1430,  &
1.0011, 1.0023, 1.0032, 1.0051, 1.0081, 1.0131, 1.0227, 1.0426, 1.0881,  &
1.0011, 1.0016, 1.0029, 1.0042, 1.0065, 1.0103, 1.0171, 1.0303, 1.0584,  &
1.0011, 1.0011, 1.0024, 1.0032, 1.0052, 1.0084, 1.0135, 1.0229, 1.0417,  &
1.0011, 1.0011, 1.0016, 1.0030, 1.0046, 1.0069, 1.0109, 1.0181, 1.0316 /
data ((supersat( 2,iw,iconc,1, 1),iw=1,9),iconc=1,7)/    &
1.0106, 1.0161, 1.0261, 1.0441, 1.0808, 1.1669, 1.4170, 2.7133, 4.9906,  &
1.0087, 1.0132, 1.0212, 1.0347, 1.0608, 1.1164, 1.2559, 1.7237, 4.0288,  &
1.0074, 1.0114, 1.0175, 1.0283, 1.0479, 1.0870, 1.1758, 1.4200, 2.5549,  &
1.0055, 1.0092, 1.0150, 1.0235, 1.0390, 1.0683, 1.1299, 1.2807, 1.7658,  &
1.0054, 1.0085, 1.0126, 1.0198, 1.0326, 1.0554, 1.1008, 1.2031, 1.4825,  &
1.0050, 1.0070, 1.0111, 1.0169, 1.0274, 1.0459, 1.0807, 1.1545, 1.3364,  &
1.0043, 1.0055, 1.0091, 1.0149, 1.0232, 1.0385, 1.0661, 1.1216, 1.2485 /
data ((supersat( 2,iw,iconc,2, 1),iw=1,9),iconc=1,7)/    &
1.0082, 1.0126, 1.0198, 1.0327, 1.0582, 1.1147, 1.2642, 1.8075, 3.7194,  &
1.0069, 1.0102, 1.0159, 1.0257, 1.0438, 1.0809, 1.1678, 1.4173, 2.5959,  &
1.0057, 1.0085, 1.0132, 1.0209, 1.0345, 1.0609, 1.1173, 1.2588, 1.7239,  &
1.0041, 1.0072, 1.0107, 1.0172, 1.0282, 1.0480, 1.0878, 1.1785, 1.4267,  &
1.0040, 1.0064, 1.0096, 1.0146, 1.0234, 1.0391, 1.0689, 1.1319, 1.2865,  &
1.0038, 1.0049, 1.0075, 1.0128, 1.0199, 1.0324, 1.0557, 1.1022, 1.2075,  &
1.0034, 1.0040, 1.0071, 1.0106, 1.0169, 1.0275, 1.0461, 1.0817, 1.1576 /
data ((supersat( 2,iw,iconc,3, 1),iw=1,9),iconc=1,7)/    &
1.0067, 1.0101, 1.0157, 1.0254, 1.0441, 1.0840, 1.1837, 1.4882, 2.7399,  &
1.0055, 1.0082, 1.0125, 1.0198, 1.0331, 1.0594, 1.1183, 1.2739, 1.8188,  &
1.0047, 1.0066, 1.0103, 1.0160, 1.0261, 1.0446, 1.0832, 1.1742, 1.4332,  &
1.0033, 1.0056, 1.0086, 1.0133, 1.0211, 1.0351, 1.0625, 1.1216, 1.2700,  &
1.0028, 1.0050, 1.0073, 1.0112, 1.0174, 1.0286, 1.0491, 1.0907, 1.1861,  &
1.0027, 1.0038, 1.0059, 1.0092, 1.0146, 1.0239, 1.0400, 1.0710, 1.1372,  &
1.0026, 1.0028, 1.0056, 1.0084, 1.0127, 1.0201, 1.0332, 1.0573, 1.1060 /
data ((supersat( 2,iw,iconc,4, 1),iw=1,9),iconc=1,7)/    &
1.0057, 1.0084, 1.0129, 1.0205, 1.0350, 1.0652, 1.1376, 1.3405, 2.0945,  &
1.0044, 1.0067, 1.0102, 1.0159, 1.0261, 1.0458, 1.0888, 1.1968, 1.5214,  &
1.0039, 1.0055, 1.0083, 1.0127, 1.0203, 1.0342, 1.0622, 1.1259, 1.2942,  &
1.0031, 1.0044, 1.0068, 1.0104, 1.0164, 1.0268, 1.0465, 1.0878, 1.1864,  &
1.0022, 1.0041, 1.0058, 1.0088, 1.0136, 1.0217, 1.0365, 1.0656, 1.1293,  &
1.0016, 1.0034, 1.0045, 1.0072, 1.0115, 1.0180, 1.0296, 1.0513, 1.0959,  &
1.0016, 1.0023, 1.0043, 1.0066, 1.0097, 1.0151, 1.0247, 1.0417, 1.0746 /
data ((supersat( 2,iw,iconc,5, 1),iw=1,9),iconc=1,7)/    &
1.0049, 1.0071, 1.0109, 1.0172, 1.0290, 1.0531, 1.1094, 1.2603, 1.7391,  &
1.0039, 1.0055, 1.0086, 1.0132, 1.0214, 1.0371, 1.0704, 1.1519, 1.3794,  &
1.0031, 1.0047, 1.0070, 1.0104, 1.0165, 1.0274, 1.0489, 1.0968, 1.2182,  &
1.0028, 1.0035, 1.0055, 1.0084, 1.0132, 1.0212, 1.0363, 1.0670, 1.1381,  &
1.0023, 1.0031, 1.0048, 1.0072, 1.0108, 1.0170, 1.0282, 1.0497, 1.0952,  &
1.0016, 1.0029, 1.0036, 1.0057, 1.0091, 1.0142, 1.0229, 1.0387, 1.0705,  &
1.0009, 1.0025, 1.0032, 1.0053, 1.0078, 1.0120, 1.0190, 1.0314, 1.0549 /
data ((supersat( 2,iw,iconc,6, 1),iw=1,9),iconc=1,7)/    &
1.0042, 1.0062, 1.0095, 1.0149, 1.0248, 1.0448, 1.0911, 1.2116, 1.5646,  &
1.0035, 1.0050, 1.0074, 1.0113, 1.0182, 1.0311, 1.0584, 1.1237, 1.2997,  &
1.0027, 1.0041, 1.0059, 1.0089, 1.0138, 1.0228, 1.0402, 1.0783, 1.1730,  &
1.0020, 1.0034, 1.0046, 1.0072, 1.0110, 1.0175, 1.0295, 1.0536, 1.1085,  &
1.0019, 1.0023, 1.0041, 1.0060, 1.0088, 1.0140, 1.0227, 1.0393, 1.0741,  &
1.0018, 1.0021, 1.0034, 1.0047, 1.0074, 1.0114, 1.0183, 1.0305, 1.0544,  &
1.0015, 1.0020, 1.0024, 1.0042, 1.0064, 1.0096, 1.0150, 1.0246, 1.0422 /
data ((supersat( 2,iw,iconc,7, 1),iw=1,9),iconc=1,7)/    &
1.0039, 1.0056, 1.0083, 1.0131, 1.0217, 1.0389, 1.0781, 1.1790, 1.4626,  &
1.0029, 1.0045, 1.0065, 1.0099, 1.0159, 1.0269, 1.0500, 1.1047, 1.2490,  &
1.0026, 1.0035, 1.0050, 1.0078, 1.0120, 1.0196, 1.0342, 1.0658, 1.1436,  &
1.0020, 1.0030, 1.0043, 1.0063, 1.0094, 1.0149, 1.0248, 1.0447, 1.0894,  &
1.0014, 1.0025, 1.0033, 1.0051, 1.0076, 1.0118, 1.0190, 1.0324, 1.0603,  &
1.0009, 1.0018, 1.0030, 1.0043, 1.0063, 1.0095, 1.0151, 1.0249, 1.0438,  &
1.0009, 1.0011, 1.0025, 1.0033, 1.0052, 1.0080, 1.0124, 1.0200, 1.0338 /
data ((supersat( 2,iw,iconc,1, 2),iw=1,9),iconc=1,7)/    &
1.0089, 1.0146, 1.0240, 1.0417, 1.0780, 1.1645, 1.4154, 2.7117, 4.9895,  &
1.0079, 1.0120, 1.0191, 1.0324, 1.0581, 1.1137, 1.2538, 1.7219, 4.0268,  &
1.0057, 1.0097, 1.0158, 1.0261, 1.0454, 1.0843, 1.1733, 1.4179, 2.5518,  &
1.0054, 1.0087, 1.0132, 1.0219, 1.0367, 1.0655, 1.1271, 1.2781, 1.7632,  &
1.0050, 1.0073, 1.0116, 1.0185, 1.0303, 1.0527, 1.0978, 1.2002, 1.4795,  &
1.0043, 1.0055, 1.0092, 1.0156, 1.0255, 1.0432, 1.0777, 1.1513, 1.3330,  &
1.0031, 1.0052, 1.0085, 1.0129, 1.0215, 1.0359, 1.0629, 1.1180, 1.2447 /
data ((supersat( 2,iw,iconc,2, 2),iw=1,9),iconc=1,7)/    &
1.0070, 1.0111, 1.0179, 1.0304, 1.0555, 1.1122, 1.2625, 1.8066, 3.7184,  &
1.0060, 1.0092, 1.0144, 1.0237, 1.0414, 1.0782, 1.1656, 1.4158, 2.5939,  &
1.0042, 1.0073, 1.0119, 1.0191, 1.0324, 1.0583, 1.1149, 1.2569, 1.7221,  &
1.0040, 1.0066, 1.0100, 1.0159, 1.0262, 1.0456, 1.0852, 1.1761, 1.4247,  &
1.0038, 1.0052, 1.0082, 1.0134, 1.0218, 1.0369, 1.0663, 1.1292, 1.2841,  &
1.0033, 1.0040, 1.0072, 1.0111, 1.0184, 1.0305, 1.0531, 1.0994, 1.2046,  &
1.0026, 1.0040, 1.0064, 1.0098, 1.0156, 1.0255, 1.0435, 1.0787, 1.1544 /
data ((supersat( 2,iw,iconc,3, 2),iw=1,9),iconc=1,7)/    &
1.0055, 1.0086, 1.0139, 1.0232, 1.0415, 1.0815, 1.1819, 1.4874, 2.7394,  &
1.0048, 1.0072, 1.0110, 1.0181, 1.0308, 1.0568, 1.1161, 1.2725, 1.8175,  &
1.0035, 1.0057, 1.0088, 1.0144, 1.0240, 1.0423, 1.0807, 1.1723, 1.4319,  &
1.0028, 1.0051, 1.0077, 1.0117, 1.0194, 1.0330, 1.0600, 1.1193, 1.2683,  &
1.0027, 1.0040, 1.0059, 1.0100, 1.0162, 1.0268, 1.0468, 1.0882, 1.1839,  &
1.0026, 1.0027, 1.0056, 1.0085, 1.0136, 1.0222, 1.0378, 1.0684, 1.1347,  &
1.0023, 1.0027, 1.0049, 1.0074, 1.0114, 1.0186, 1.0310, 1.0547, 1.1032 /
data ((supersat( 2,iw,iconc,4, 2),iw=1,9),iconc=1,7)/    &
1.0045, 1.0071, 1.0112, 1.0186, 1.0327, 1.0627, 1.1357, 1.3395, 2.0938,  &
1.0040, 1.0059, 1.0089, 1.0142, 1.0241, 1.0434, 1.0866, 1.1953, 1.5205,  &
1.0032, 1.0044, 1.0070, 1.0114, 1.0185, 1.0320, 1.0598, 1.1239, 1.2929,  &
1.0021, 1.0041, 1.0061, 1.0093, 1.0148, 1.0249, 1.0442, 1.0855, 1.1846,  &
1.0016, 1.0034, 1.0045, 1.0075, 1.0121, 1.0200, 1.0345, 1.0632, 1.1271,  &
1.0016, 1.0024, 1.0043, 1.0067, 1.0101, 1.0166, 1.0278, 1.0491, 1.0935,  &
1.0015, 1.0016, 1.0039, 1.0055, 1.0089, 1.0140, 1.0229, 1.0394, 1.0721 /
data ((supersat( 2,iw,iconc,5, 2),iw=1,9),iconc=1,7)/    &
1.0041, 1.0061, 1.0094, 1.0154, 1.0266, 1.0507, 1.1076, 1.2594, 1.7385,  &
1.0031, 1.0049, 1.0074, 1.0117, 1.0195, 1.0348, 1.0682, 1.1504, 1.3787,  &
1.0028, 1.0038, 1.0058, 1.0093, 1.0149, 1.0254, 1.0467, 1.0947, 1.2170,  &
1.0022, 1.0031, 1.0050, 1.0075, 1.0119, 1.0196, 1.0342, 1.0648, 1.1363,  &
1.0015, 1.0029, 1.0038, 1.0060, 1.0097, 1.0157, 1.0264, 1.0475, 1.0931,  &
1.0009, 1.0024, 1.0032, 1.0053, 1.0080, 1.0127, 1.0212, 1.0367, 1.0682,  &
1.0006, 1.0017, 1.0031, 1.0045, 1.0069, 1.0106, 1.0173, 1.0295, 1.0526 /
data ((supersat( 2,iw,iconc,6, 2),iw=1,9),iconc=1,7)/    &
1.0036, 1.0054, 1.0082, 1.0132, 1.0227, 1.0426, 1.0892, 1.2107, 1.5642,  &
1.0028, 1.0041, 1.0062, 1.0099, 1.0164, 1.0290, 1.0562, 1.1222, 1.2989,  &
1.0021, 1.0035, 1.0051, 1.0078, 1.0124, 1.0209, 1.0380, 1.0763, 1.1718,  &
1.0020, 1.0025, 1.0041, 1.0062, 1.0098, 1.0159, 1.0275, 1.0515, 1.1068,  &
1.0018, 1.0021, 1.0035, 1.0051, 1.0080, 1.0125, 1.0210, 1.0373, 1.0720,  &
1.0014, 1.0020, 1.0025, 1.0043, 1.0065, 1.0104, 1.0167, 1.0286, 1.0522,  &
1.0009, 1.0019, 1.0021, 1.0038, 1.0056, 1.0086, 1.0138, 1.0228, 1.0401 /
data ((supersat( 2,iw,iconc,7, 2),iw=1,9),iconc=1,7)/    &
1.0030, 1.0047, 1.0072, 1.0116, 1.0198, 1.0368, 1.0764, 1.1781, 1.4621,  &
1.0026, 1.0037, 1.0056, 1.0087, 1.0142, 1.0249, 1.0479, 1.1032, 1.2483,  &
1.0020, 1.0030, 1.0045, 1.0067, 1.0107, 1.0178, 1.0321, 1.0639, 1.1424,  &
1.0013, 1.0025, 1.0033, 1.0052, 1.0083, 1.0134, 1.0230, 1.0426, 1.0876,  &
1.0011, 1.0018, 1.0030, 1.0045, 1.0067, 1.0105, 1.0174, 1.0305, 1.0583,  &
1.0011, 1.0011, 1.0026, 1.0033, 1.0053, 1.0086, 1.0138, 1.0232, 1.0418,  &
1.0011, 1.0011, 1.0018, 1.0031, 1.0047, 1.0071, 1.0111, 1.0184, 1.0319 /
data ((supersat( 2,iw,iconc,1, 3),iw=1,9),iconc=1,7)/    &
1.0083, 1.0130, 1.0219, 1.0391, 1.0753, 1.1622, 1.4139, 2.7103, 4.9885,  &
1.0065, 1.0108, 1.0177, 1.0303, 1.0555, 1.1111, 1.2518, 1.7203, 4.0242,  &
1.0054, 1.0089, 1.0146, 1.0244, 1.0431, 1.0816, 1.1708, 1.4158, 2.5489,  &
1.0051, 1.0078, 1.0122, 1.0200, 1.0346, 1.0629, 1.1244, 1.2755, 1.7601,  &
1.0043, 1.0057, 1.0101, 1.0166, 1.0283, 1.0501, 1.0948, 1.1972, 1.4763,  &
1.0032, 1.0053, 1.0086, 1.0141, 1.0234, 1.0407, 1.0746, 1.1479, 1.3294,  &
1.0023, 1.0047, 1.0074, 1.0117, 1.0192, 1.0333, 1.0597, 1.1143, 1.2404 /
data ((supersat( 2,iw,iconc,2, 3),iw=1,9),iconc=1,7)/    &
1.0063, 1.0098, 1.0162, 1.0283, 1.0530, 1.1099, 1.2611, 1.8055, 3.7170,  &
1.0047, 1.0078, 1.0129, 1.0217, 1.0391, 1.0758, 1.1636, 1.4144, 2.5925,  &
1.0040, 1.0068, 1.0104, 1.0176, 1.0303, 1.0559, 1.1125, 1.2550, 1.7205,  &
1.0038, 1.0055, 1.0090, 1.0144, 1.0244, 1.0433, 1.0827, 1.1737, 1.4227,  &
1.0033, 1.0040, 1.0073, 1.0122, 1.0199, 1.0346, 1.0637, 1.1266, 1.2816,  &
1.0025, 1.0039, 1.0066, 1.0102, 1.0166, 1.0284, 1.0506, 1.0965, 1.2017,  &
1.0018, 1.0036, 1.0053, 1.0086, 1.0137, 1.0234, 1.0409, 1.0757, 1.1511 /
data ((supersat( 2,iw,iconc,3, 3),iw=1,9),iconc=1,7)/    &
1.0050, 1.0078, 1.0125, 1.0214, 1.0392, 1.0793, 1.1804, 1.4865, 2.7381,  &
1.0038, 1.0058, 1.0100, 1.0163, 1.0288, 1.0545, 1.1141, 1.2712, 1.8169,  &
1.0027, 1.0053, 1.0081, 1.0132, 1.0222, 1.0401, 1.0785, 1.1704, 1.4306,  &
1.0027, 1.0042, 1.0065, 1.0109, 1.0179, 1.0311, 1.0577, 1.1171, 1.2665,  &
1.0026, 1.0028, 1.0057, 1.0087, 1.0145, 1.0249, 1.0446, 1.0859, 1.1817,  &
1.0023, 1.0027, 1.0050, 1.0078, 1.0122, 1.0204, 1.0356, 1.0659, 1.1321,  &
1.0018, 1.0026, 1.0038, 1.0060, 1.0104, 1.0168, 1.0291, 1.0521, 1.1003 /
data ((supersat( 2,iw,iconc,4, 3),iw=1,9),iconc=1,7)/    &
1.0040, 1.0063, 1.0100, 1.0169, 1.0305, 1.0605, 1.1342, 1.3388, 2.0929,  &
1.0033, 1.0048, 1.0079, 1.0129, 1.0222, 1.0413, 1.0845, 1.1941, 1.5197,  &
1.0023, 1.0042, 1.0065, 1.0102, 1.0170, 1.0301, 1.0577, 1.1221, 1.2917,  &
1.0016, 1.0035, 1.0050, 1.0084, 1.0136, 1.0232, 1.0422, 1.0834, 1.1830,  &
1.0016, 1.0025, 1.0043, 1.0069, 1.0112, 1.0186, 1.0325, 1.0610, 1.1251,  &
1.0016, 1.0016, 1.0039, 1.0059, 1.0093, 1.0151, 1.0260, 1.0469, 1.0912,  &
1.0015, 1.0016, 1.0032, 1.0044, 1.0076, 1.0125, 1.0213, 1.0373, 1.0696 /
data ((supersat( 2,iw,iconc,5, 3),iw=1,9),iconc=1,7)/    &
1.0031, 1.0052, 1.0083, 1.0139, 1.0248, 1.0485, 1.1060, 1.2586, 1.7381,  &
1.0028, 1.0041, 1.0065, 1.0105, 1.0178, 1.0328, 1.0662, 1.1491, 1.3780,  &
1.0023, 1.0031, 1.0052, 1.0081, 1.0135, 1.0236, 1.0446, 1.0930, 1.2159,  &
1.0015, 1.0029, 1.0042, 1.0067, 1.0106, 1.0180, 1.0322, 1.0627, 1.1348,  &
1.0009, 1.0024, 1.0032, 1.0055, 1.0087, 1.0143, 1.0247, 1.0454, 1.0912,  &
1.0006, 1.0017, 1.0031, 1.0047, 1.0073, 1.0118, 1.0196, 1.0348, 1.0661,  &
1.0006, 1.0011, 1.0027, 1.0034, 1.0057, 1.0098, 1.0161, 1.0276, 1.0504 /
data ((supersat( 2,iw,iconc,6, 3),iw=1,9),iconc=1,7)/    &
1.0030, 1.0044, 1.0071, 1.0118, 1.0209, 1.0405, 1.0877, 1.2099, 1.5637,  &
1.0021, 1.0037, 1.0055, 1.0088, 1.0148, 1.0271, 1.0543, 1.1209, 1.2982,  &
1.0020, 1.0027, 1.0042, 1.0068, 1.0112, 1.0193, 1.0361, 1.0745, 1.1707,  &
1.0017, 1.0021, 1.0037, 1.0055, 1.0086, 1.0145, 1.0257, 1.0495, 1.1053,  &
1.0013, 1.0020, 1.0027, 1.0043, 1.0070, 1.0115, 1.0195, 1.0354, 1.0702,  &
1.0009, 1.0018, 1.0021, 1.0039, 1.0059, 1.0093, 1.0154, 1.0269, 1.0503,  &
1.0006, 1.0015, 1.0020, 1.0031, 1.0046, 1.0078, 1.0125, 1.0213, 1.0381 /
data ((supersat( 2,iw,iconc,7, 3),iw=1,9),iconc=1,7)/    &
1.0027, 1.0040, 1.0062, 1.0103, 1.0181, 1.0349, 1.0749, 1.1773, 1.4616,  &
1.0021, 1.0030, 1.0048, 1.0076, 1.0127, 1.0231, 1.0460, 1.1019, 1.2477,  &
1.0013, 1.0026, 1.0037, 1.0059, 1.0095, 1.0163, 1.0303, 1.0622, 1.1414,  &
1.0011, 1.0019, 1.0031, 1.0047, 1.0073, 1.0121, 1.0213, 1.0408, 1.0861,  &
1.0011, 1.0011, 1.0026, 1.0036, 1.0059, 1.0095, 1.0159, 1.0288, 1.0565,  &
1.0011, 1.0011, 1.0019, 1.0031, 1.0049, 1.0076, 1.0125, 1.0216, 1.0400,  &
1.0010, 1.0011, 1.0012, 1.0028, 1.0040, 1.0064, 1.0102, 1.0170, 1.0301 /
data ((supersat( 2,iw,iconc,1, 4),iw=1,9),iconc=1,7)/    &
1.0073, 1.0118, 1.0204, 1.0370, 1.0731, 1.1604, 1.4126, 2.7080, 4.9873,  &
1.0054, 1.0093, 1.0161, 1.0283, 1.0533, 1.1089, 1.2500, 1.7186, 4.0215,  &
1.0051, 1.0082, 1.0131, 1.0226, 1.0409, 1.0793, 1.1686, 1.4138, 2.5462,  &
1.0044, 1.0064, 1.0112, 1.0185, 1.0324, 1.0605, 1.1218, 1.2731, 1.7572,  &
1.0033, 1.0053, 1.0088, 1.0153, 1.0263, 1.0477, 1.0921, 1.1942, 1.4729,  &
1.0023, 1.0048, 1.0078, 1.0122, 1.0215, 1.0382, 1.0716, 1.1445, 1.3255,  &
1.0018, 1.0039, 1.0060, 1.0103, 1.0175, 1.0308, 1.0565, 1.1104, 1.2358 /
data ((supersat( 2,iw,iconc,2, 4),iw=1,9),iconc=1,7)/    &
1.0054, 1.0089, 1.0149, 1.0263, 1.0509, 1.1081, 1.2598, 1.8046, 3.7161,  &
1.0040, 1.0067, 1.0119, 1.0203, 1.0371, 1.0737, 1.1619, 1.4132, 2.5909,  &
1.0038, 1.0061, 1.0097, 1.0161, 1.0286, 1.0538, 1.1104, 1.2533, 1.7189,  &
1.0034, 1.0048, 1.0076, 1.0132, 1.0227, 1.0413, 1.0805, 1.1717, 1.4207,  &
1.0026, 1.0034, 1.0068, 1.0107, 1.0186, 1.0328, 1.0614, 1.1242, 1.2792,  &
1.0018, 1.0033, 1.0056, 1.0092, 1.0153, 1.0265, 1.0482, 1.0937, 1.1988,  &
1.0014, 1.0029, 1.0041, 1.0071, 1.0125, 1.0215, 1.0385, 1.0726, 1.1477 /
data ((supersat( 2,iw,iconc,3, 4),iw=1,9),iconc=1,7)/    &
1.0043, 1.0068, 1.0112, 1.0197, 1.0373, 1.0774, 1.1792, 1.4855, 2.7372,  &
1.0028, 1.0055, 1.0087, 1.0151, 1.0271, 1.0526, 1.1124, 1.2701, 1.8159,  &
1.0027, 1.0046, 1.0073, 1.0119, 1.0207, 1.0382, 1.0765, 1.1689, 1.4294,  &
1.0026, 1.0032, 1.0057, 1.0098, 1.0164, 1.0293, 1.0557, 1.1152, 1.2649,  &
1.0023, 1.0027, 1.0052, 1.0082, 1.0135, 1.0234, 1.0426, 1.0837, 1.1797,  &
1.0018, 1.0026, 1.0042, 1.0067, 1.0111, 1.0189, 1.0337, 1.0636, 1.1297,  &
1.0013, 1.0023, 1.0029, 1.0054, 1.0090, 1.0155, 1.0270, 1.0498, 1.0975 /
data ((supersat( 2,iw,iconc,4, 4),iw=1,9),iconc=1,7)/    &
1.0036, 1.0055, 1.0088, 1.0153, 1.0287, 1.0587, 1.1330, 1.3381, 2.0927,  &
1.0025, 1.0042, 1.0068, 1.0115, 1.0206, 1.0394, 1.0829, 1.1930, 1.5189,  &
1.0016, 1.0037, 1.0057, 1.0091, 1.0157, 1.0284, 1.0559, 1.1206, 1.2908,  &
1.0016, 1.0027, 1.0044, 1.0072, 1.0123, 1.0216, 1.0404, 1.0817, 1.1816,  &
1.0016, 1.0016, 1.0040, 1.0063, 1.0100, 1.0171, 1.0307, 1.0591, 1.1233,  &
1.0015, 1.0016, 1.0033, 1.0048, 1.0084, 1.0140, 1.0243, 1.0449, 1.0891,  &
1.0013, 1.0016, 1.0024, 1.0041, 1.0067, 1.0114, 1.0196, 1.0353, 1.0673 /
data ((supersat( 2,iw,iconc,5, 4),iw=1,9),iconc=1,7)/    &
1.0029, 1.0046, 1.0073, 1.0125, 1.0231, 1.0468, 1.1048, 1.2579, 1.7377,  &
1.0024, 1.0032, 1.0055, 1.0093, 1.0164, 1.0310, 1.0646, 1.1481, 1.3774,  &
1.0017, 1.0030, 1.0046, 1.0073, 1.0123, 1.0221, 1.0429, 1.0916, 1.2151,  &
1.0009, 1.0025, 1.0032, 1.0056, 1.0097, 1.0166, 1.0306, 1.0611, 1.1335,  &
1.0006, 1.0018, 1.0031, 1.0050, 1.0078, 1.0132, 1.0231, 1.0437, 1.0896,  &
1.0006, 1.0012, 1.0028, 1.0038, 1.0064, 1.0106, 1.0183, 1.0331, 1.0642,  &
1.0006, 1.0008, 1.0022, 1.0031, 1.0053, 1.0087, 1.0146, 1.0259, 1.0484 /
data ((supersat( 2,iw,iconc,6, 4),iw=1,9),iconc=1,7)/    &
1.0022, 1.0039, 1.0061, 1.0106, 1.0193, 1.0389, 1.0865, 1.2093, 1.5633,  &
1.0020, 1.0030, 1.0047, 1.0078, 1.0135, 1.0255, 1.0527, 1.1199, 1.2978,  &
1.0018, 1.0021, 1.0039, 1.0060, 1.0100, 1.0179, 1.0345, 1.0732, 1.1699,  &
1.0014, 1.0020, 1.0030, 1.0047, 1.0079, 1.0134, 1.0242, 1.0480, 1.1040,  &
1.0009, 1.0018, 1.0021, 1.0040, 1.0063, 1.0104, 1.0181, 1.0338, 1.0686,  &
1.0007, 1.0015, 1.0020, 1.0033, 1.0052, 1.0084, 1.0142, 1.0254, 1.0486,  &
1.0005, 1.0011, 1.0019, 1.0024, 1.0042, 1.0068, 1.0115, 1.0198, 1.0364 /
data ((supersat( 2,iw,iconc,7, 4),iw=1,9),iconc=1,7)/    &
1.0022, 1.0033, 1.0054, 1.0091, 1.0166, 1.0333, 1.0737, 1.1768, 1.4613,  &
1.0015, 1.0027, 1.0042, 1.0066, 1.0116, 1.0216, 1.0446, 1.1010, 1.2473,  &
1.0011, 1.0020, 1.0031, 1.0050, 1.0084, 1.0150, 1.0288, 1.0609, 1.1406,  &
1.0011, 1.0013, 1.0027, 1.0041, 1.0066, 1.0111, 1.0200, 1.0393, 1.0850,  &
1.0011, 1.0011, 1.0021, 1.0032, 1.0051, 1.0086, 1.0148, 1.0273, 1.0551,  &
1.0010, 1.0011, 1.0013, 1.0029, 1.0044, 1.0069, 1.0115, 1.0203, 1.0384,  &
1.0008, 1.0010, 1.0011, 1.0023, 1.0032, 1.0055, 1.0092, 1.0158, 1.0285 /
data ((supersat( 2,iw,iconc,1, 5),iw=1,9),iconc=1,7)/    &
1.0066, 1.0113, 1.0195, 1.0359, 1.0720, 1.1595, 1.4120, 2.7073, 4.9867,  &
1.0053, 1.0089, 1.0153, 1.0274, 1.0521, 1.1078, 1.2491, 1.7178, 4.0198,  &
1.0049, 1.0077, 1.0123, 1.0217, 1.0398, 1.0780, 1.1675, 1.4128, 2.5441,  &
1.0040, 1.0056, 1.0104, 1.0178, 1.0314, 1.0593, 1.1205, 1.2718, 1.7553,  &
1.0028, 1.0051, 1.0085, 1.0146, 1.0252, 1.0464, 1.0906, 1.1926, 1.4708,  &
1.0020, 1.0044, 1.0071, 1.0117, 1.0205, 1.0369, 1.0699, 1.1426, 1.3231,  &
1.0017, 1.0034, 1.0052, 1.0094, 1.0164, 1.0293, 1.0547, 1.1082, 1.2330 /
data ((supersat( 2,iw,iconc,2, 5),iw=1,9),iconc=1,7)/    &
1.0047, 1.0082, 1.0141, 1.0255, 1.0499, 1.1072, 1.2593, 1.8043, 3.7155,  &
1.0039, 1.0068, 1.0111, 1.0195, 1.0361, 1.0727, 1.1611, 1.4126, 2.5900,  &
1.0037, 1.0055, 1.0092, 1.0155, 1.0275, 1.0528, 1.1094, 1.2525, 1.7181,  &
1.0030, 1.0040, 1.0072, 1.0127, 1.0218, 1.0403, 1.0793, 1.1706, 1.4196,  &
1.0022, 1.0038, 1.0064, 1.0102, 1.0178, 1.0317, 1.0601, 1.1229, 1.2779,  &
1.0016, 1.0034, 1.0050, 1.0086, 1.0145, 1.0253, 1.0469, 1.0922, 1.1972,  &
1.0013, 1.0027, 1.0036, 1.0066, 1.0117, 1.0204, 1.0371, 1.0710, 1.1457 /
data ((supersat( 2,iw,iconc,3, 5),iw=1,9),iconc=1,7)/    &
1.0037, 1.0062, 1.0106, 1.0189, 1.0363, 1.0766, 1.1786, 1.4852, 2.7370,  &
1.0028, 1.0052, 1.0083, 1.0144, 1.0262, 1.0516, 1.1116, 1.2695, 1.8153,  &
1.0027, 1.0041, 1.0067, 1.0112, 1.0200, 1.0373, 1.0756, 1.1682, 1.4288,  &
1.0024, 1.0027, 1.0056, 1.0091, 1.0159, 1.0285, 1.0547, 1.1142, 1.2641,  &
1.0020, 1.0027, 1.0048, 1.0078, 1.0129, 1.0224, 1.0416, 1.0826, 1.1786,  &
1.0015, 1.0025, 1.0036, 1.0059, 1.0106, 1.0182, 1.0326, 1.0624, 1.1284,  &
1.0012, 1.0021, 1.0025, 1.0051, 1.0082, 1.0146, 1.0260, 1.0484, 1.0960 /
data ((supersat( 2,iw,iconc,4, 5),iw=1,9),iconc=1,7)/    &
1.0032, 1.0049, 1.0084, 1.0147, 1.0278, 1.0579, 1.1324, 1.3378, 2.0924,  &
1.0021, 1.0041, 1.0065, 1.0110, 1.0198, 1.0386, 1.0822, 1.1925, 1.5185,  &
1.0016, 1.0034, 1.0051, 1.0087, 1.0150, 1.0276, 1.0550, 1.1200, 1.2903,  &
1.0016, 1.0022, 1.0042, 1.0069, 1.0117, 1.0208, 1.0395, 1.0808, 1.1810,  &
1.0015, 1.0016, 1.0038, 1.0059, 1.0094, 1.0164, 1.0299, 1.0581, 1.1224,  &
1.0014, 1.0016, 1.0029, 1.0043, 1.0078, 1.0134, 1.0235, 1.0439, 1.0880,  &
1.0012, 1.0014, 1.0021, 1.0039, 1.0064, 1.0108, 1.0187, 1.0342, 1.0660 /
data ((supersat( 2,iw,iconc,5, 5),iw=1,9),iconc=1,7)/    &
1.0027, 1.0042, 1.0069, 1.0119, 1.0223, 1.0460, 1.1042, 1.2577, 1.7376,  &
1.0021, 1.0031, 1.0052, 1.0089, 1.0157, 1.0302, 1.0639, 1.1476, 1.3771,  &
1.0013, 1.0028, 1.0042, 1.0070, 1.0117, 1.0213, 1.0421, 1.0909, 1.2147,  &
1.0008, 1.0022, 1.0032, 1.0054, 1.0092, 1.0160, 1.0298, 1.0603, 1.1329,  &
1.0006, 1.0015, 1.0029, 1.0046, 1.0075, 1.0125, 1.0224, 1.0428, 1.0888,  &
1.0006, 1.0010, 1.0025, 1.0033, 1.0059, 1.0100, 1.0176, 1.0322, 1.0633,  &
1.0005, 1.0008, 1.0019, 1.0029, 1.0050, 1.0081, 1.0140, 1.0250, 1.0473 /
data ((supersat( 2,iw,iconc,6, 5),iw=1,9),iconc=1,7)/    &
1.0021, 1.0036, 1.0058, 1.0100, 1.0186, 1.0382, 1.0860, 1.2091, 1.5631,  &
1.0019, 1.0026, 1.0042, 1.0074, 1.0129, 1.0247, 1.0521, 1.1195, 1.2976,  &
1.0016, 1.0021, 1.0036, 1.0057, 1.0096, 1.0172, 1.0337, 1.0726, 1.1695,  &
1.0011, 1.0019, 1.0026, 1.0043, 1.0075, 1.0128, 1.0235, 1.0472, 1.1035,  &
1.0007, 1.0017, 1.0021, 1.0038, 1.0060, 1.0100, 1.0175, 1.0330, 1.0679,  &
1.0005, 1.0013, 1.0020, 1.0029, 1.0047, 1.0080, 1.0136, 1.0246, 1.0477,  &
1.0005, 1.0010, 1.0017, 1.0020, 1.0039, 1.0062, 1.0108, 1.0191, 1.0354 /
data ((supersat( 2,iw,iconc,7, 5),iw=1,9),iconc=1,7)/    &
1.0019, 1.0030, 1.0049, 1.0086, 1.0159, 1.0326, 1.0732, 1.1766, 1.4611,  &
1.0012, 1.0025, 1.0038, 1.0063, 1.0110, 1.0209, 1.0439, 1.1006, 1.2471,  &
1.0011, 1.0017, 1.0030, 1.0048, 1.0080, 1.0143, 1.0281, 1.0603, 1.1403,  &
1.0011, 1.0011, 1.0025, 1.0037, 1.0062, 1.0105, 1.0193, 1.0386, 1.0845,  &
1.0010, 1.0011, 1.0017, 1.0031, 1.0049, 1.0082, 1.0142, 1.0266, 1.0544,  &
1.0009, 1.0011, 1.0011, 1.0026, 1.0040, 1.0065, 1.0110, 1.0196, 1.0377,  &
1.0007, 1.0009, 1.0010, 1.0020, 1.0030, 1.0050, 1.0086, 1.0151, 1.0277 /
data ((supersat( 2,iw,iconc,1, 6),iw=1,9),iconc=1,7)/    &
1.0061, 1.0109, 1.0189, 1.0353, 1.0713, 1.1589, 1.4115, 2.7068, 4.9861,  &
1.0052, 1.0087, 1.0150, 1.0269, 1.0514, 1.1072, 1.2486, 1.7171, 4.0183,  &
1.0047, 1.0073, 1.0121, 1.0212, 1.0391, 1.0773, 1.1668, 1.4120, 2.5431,  &
1.0036, 1.0054, 1.0099, 1.0172, 1.0308, 1.0584, 1.1196, 1.2709, 1.7540,  &
1.0025, 1.0049, 1.0082, 1.0140, 1.0246, 1.0455, 1.0895, 1.1915, 1.4692,  &
1.0019, 1.0041, 1.0066, 1.0112, 1.0197, 1.0359, 1.0688, 1.1413, 1.3214,  &
1.0016, 1.0031, 1.0048, 1.0087, 1.0155, 1.0282, 1.0533, 1.1066, 1.2310 /
data ((supersat( 2,iw,iconc,2, 6),iw=1,9),iconc=1,7)/    &
1.0043, 1.0078, 1.0136, 1.0249, 1.0492, 1.1067, 1.2589, 1.8039, 3.7151,  &
1.0039, 1.0066, 1.0106, 1.0189, 1.0355, 1.0721, 1.1606, 1.4121, 2.5894,  &
1.0035, 1.0051, 1.0088, 1.0151, 1.0270, 1.0521, 1.1088, 1.2519, 1.7174,  &
1.0028, 1.0039, 1.0071, 1.0123, 1.0214, 1.0395, 1.0786, 1.1699, 1.4188,  &
1.0019, 1.0037, 1.0061, 1.0099, 1.0172, 1.0309, 1.0593, 1.1220, 1.2770,  &
1.0014, 1.0031, 1.0046, 1.0081, 1.0138, 1.0246, 1.0460, 1.0912, 1.1960,  &
1.0012, 1.0024, 1.0034, 1.0063, 1.0111, 1.0197, 1.0360, 1.0698, 1.1443 /
data ((supersat( 2,iw,iconc,3, 6),iw=1,9),iconc=1,7)/    &
1.0034, 1.0057, 1.0103, 1.0184, 1.0357, 1.0761, 1.1783, 1.4850, 2.7369,  &
1.0027, 1.0050, 1.0081, 1.0139, 1.0256, 1.0511, 1.1112, 1.2692, 1.8151,  &
1.0026, 1.0037, 1.0064, 1.0109, 1.0194, 1.0367, 1.0751, 1.1678, 1.4283,  &
1.0023, 1.0027, 1.0054, 1.0086, 1.0154, 1.0279, 1.0541, 1.1137, 1.2637,  &
1.0018, 1.0026, 1.0045, 1.0074, 1.0124, 1.0218, 1.0409, 1.0819, 1.1779,  &
1.0013, 1.0024, 1.0033, 1.0055, 1.0102, 1.0176, 1.0319, 1.0616, 1.1276,  &
1.0011, 1.0019, 1.0024, 1.0048, 1.0079, 1.0139, 1.0252, 1.0475, 1.0949 /
data ((supersat( 2,iw,iconc,4, 6),iw=1,9),iconc=1,7)/    &
1.0029, 1.0046, 1.0081, 1.0142, 1.0273, 1.0574, 1.1321, 1.3376, 2.0924,  &
1.0018, 1.0039, 1.0063, 1.0107, 1.0193, 1.0380, 1.0817, 1.1923, 1.5182,  &
1.0016, 1.0031, 1.0047, 1.0084, 1.0146, 1.0270, 1.0545, 1.1196, 1.2900,  &
1.0016, 1.0019, 1.0042, 1.0067, 1.0114, 1.0204, 1.0389, 1.0803, 1.1805,  &
1.0015, 1.0016, 1.0035, 1.0055, 1.0092, 1.0160, 1.0294, 1.0575, 1.1218,  &
1.0013, 1.0015, 1.0026, 1.0042, 1.0074, 1.0129, 1.0229, 1.0432, 1.0873,  &
1.0010, 1.0014, 1.0019, 1.0037, 1.0061, 1.0103, 1.0182, 1.0334, 1.0652 /
data ((supersat( 2,iw,iconc,5, 6),iw=1,9),iconc=1,7)/    &
1.0026, 1.0039, 1.0066, 1.0115, 1.0218, 1.0456, 1.1040, 1.2576, 1.7374,  &
1.0019, 1.0031, 1.0051, 1.0085, 1.0152, 1.0297, 1.0635, 1.1474, 1.3769,  &
1.0011, 1.0027, 1.0039, 1.0067, 1.0113, 1.0208, 1.0416, 1.0905, 1.2144,  &
1.0006, 1.0020, 1.0031, 1.0053, 1.0089, 1.0156, 1.0293, 1.0598, 1.1325,  &
1.0006, 1.0013, 1.0028, 1.0043, 1.0072, 1.0121, 1.0219, 1.0423, 1.0883,  &
1.0006, 1.0009, 1.0023, 1.0031, 1.0055, 1.0097, 1.0170, 1.0316, 1.0626,  &
1.0005, 1.0008, 1.0017, 1.0028, 1.0047, 1.0076, 1.0135, 1.0244, 1.0466 /
data ((supersat( 2,iw,iconc,6, 6),iw=1,9),iconc=1,7)/    &
1.0020, 1.0034, 1.0055, 1.0096, 1.0181, 1.0377, 1.0857, 1.2089, 1.5631,  &
1.0018, 1.0023, 1.0041, 1.0071, 1.0125, 1.0243, 1.0517, 1.1193, 1.2974,  &
1.0014, 1.0021, 1.0034, 1.0055, 1.0092, 1.0167, 1.0333, 1.0722, 1.1693,  &
1.0009, 1.0019, 1.0023, 1.0042, 1.0071, 1.0124, 1.0231, 1.0468, 1.1032,  &
1.0006, 1.0015, 1.0020, 1.0036, 1.0058, 1.0096, 1.0171, 1.0326, 1.0674,  &
1.0005, 1.0011, 1.0019, 1.0027, 1.0043, 1.0077, 1.0132, 1.0241, 1.0472,  &
1.0005, 1.0009, 1.0016, 1.0020, 1.0037, 1.0060, 1.0103, 1.0185, 1.0348 /
data ((supersat( 2,iw,iconc,7, 6),iw=1,9),iconc=1,7)/    &
1.0017, 1.0029, 1.0047, 1.0083, 1.0155, 1.0322, 1.0730, 1.1764, 1.4610,  &
1.0011, 1.0023, 1.0035, 1.0060, 1.0106, 1.0205, 1.0436, 1.1004, 1.2469,  &
1.0011, 1.0015, 1.0029, 1.0046, 1.0078, 1.0140, 1.0277, 1.0600, 1.1401,  &
1.0010, 1.0011, 1.0023, 1.0034, 1.0060, 1.0101, 1.0189, 1.0382, 1.0842,  &
1.0009, 1.0011, 1.0015, 1.0030, 1.0048, 1.0079, 1.0138, 1.0261, 1.0540,  &
1.0008, 1.0010, 1.0011, 1.0025, 1.0037, 1.0063, 1.0105, 1.0191, 1.0372,  &
1.0006, 1.0009, 1.0010, 1.0018, 1.0029, 1.0048, 1.0083, 1.0147, 1.0272 /
data ((supersat( 2,iw,iconc,1, 7),iw=1,9),iconc=1,7)/    &
1.0057, 1.0106, 1.0184, 1.0348, 1.0708, 1.1586, 1.4112, 2.7065, 4.9856,  &
1.0052, 1.0085, 1.0147, 1.0264, 1.0509, 1.1067, 1.2482, 1.7166, 4.0174,  &
1.0045, 1.0070, 1.0118, 1.0208, 1.0386, 1.0767, 1.1662, 1.4114, 2.5420,  &
1.0033, 1.0053, 1.0094, 1.0168, 1.0302, 1.0578, 1.1190, 1.2702, 1.7530,  &
1.0023, 1.0048, 1.0080, 1.0136, 1.0241, 1.0448, 1.0888, 1.1907, 1.4680,  &
1.0018, 1.0038, 1.0062, 1.0108, 1.0191, 1.0351, 1.0679, 1.1402, 1.3200,  &
1.0015, 1.0028, 1.0044, 1.0081, 1.0148, 1.0273, 1.0524, 1.1054, 1.2294 /
data ((supersat( 2,iw,iconc,2, 7),iw=1,9),iconc=1,7)/    &
1.0040, 1.0074, 1.0132, 1.0245, 1.0488, 1.1063, 1.2587, 1.8037, 3.7149,  &
1.0038, 1.0064, 1.0102, 1.0185, 1.0351, 1.0717, 1.1603, 1.4118, 2.5889,  &
1.0034, 1.0047, 1.0085, 1.0147, 1.0266, 1.0516, 1.1084, 1.2516, 1.7169,  &
1.0025, 1.0039, 1.0069, 1.0119, 1.0210, 1.0390, 1.0781, 1.1694, 1.4182,  &
1.0017, 1.0036, 1.0058, 1.0097, 1.0167, 1.0304, 1.0587, 1.1214, 1.2763,  &
1.0014, 1.0029, 1.0042, 1.0076, 1.0133, 1.0241, 1.0453, 1.0905, 1.1952,  &
1.0012, 1.0022, 1.0033, 1.0060, 1.0106, 1.0190, 1.0353, 1.0689, 1.1432 /
data ((supersat( 2,iw,iconc,3, 7),iw=1,9),iconc=1,7)/    &
1.0031, 1.0056, 1.0100, 1.0180, 1.0353, 1.0757, 1.1781, 1.4848, 2.7368,  &
1.0027, 1.0048, 1.0079, 1.0135, 1.0252, 1.0507, 1.1109, 1.2690, 1.8149,  &
1.0025, 1.0035, 1.0060, 1.0107, 1.0190, 1.0363, 1.0747, 1.1674, 1.4280,  &
1.0022, 1.0027, 1.0053, 1.0084, 1.0151, 1.0275, 1.0536, 1.1132, 1.2632,  &
1.0017, 1.0026, 1.0043, 1.0071, 1.0121, 1.0214, 1.0405, 1.0814, 1.1774,  &
1.0012, 1.0023, 1.0030, 1.0054, 1.0098, 1.0171, 1.0314, 1.0610, 1.1269,  &
1.0010, 1.0018, 1.0023, 1.0046, 1.0076, 1.0134, 1.0246, 1.0468, 1.0942 /
data ((supersat( 2,iw,iconc,4, 7),iw=1,9),iconc=1,7)/    &
1.0027, 1.0043, 1.0078, 1.0139, 1.0269, 1.0571, 1.1319, 1.3374, 2.0923,  &
1.0016, 1.0038, 1.0061, 1.0104, 1.0190, 1.0376, 1.0815, 1.1921, 1.5182,  &
1.0016, 1.0029, 1.0044, 1.0082, 1.0142, 1.0266, 1.0541, 1.1193, 1.2898,  &
1.0016, 1.0017, 1.0041, 1.0066, 1.0112, 1.0201, 1.0385, 1.0799, 1.1802,  &
1.0014, 1.0016, 1.0034, 1.0053, 1.0090, 1.0157, 1.0289, 1.0571, 1.1214,  &
1.0012, 1.0015, 1.0024, 1.0041, 1.0070, 1.0125, 1.0225, 1.0427, 1.0867,  &
1.0010, 1.0013, 1.0018, 1.0035, 1.0058, 1.0100, 1.0177, 1.0328, 1.0645 /
data ((supersat( 2,iw,iconc,5, 7),iw=1,9),iconc=1,7)/    &
1.0025, 1.0037, 1.0064, 1.0111, 1.0214, 1.0453, 1.1038, 1.2574, 1.7375,  &
1.0017, 1.0030, 1.0049, 1.0083, 1.0149, 1.0294, 1.0632, 1.1472, 1.3767,  &
1.0010, 1.0025, 1.0036, 1.0064, 1.0111, 1.0205, 1.0413, 1.0903, 1.2142,  &
1.0006, 1.0018, 1.0031, 1.0052, 1.0086, 1.0153, 1.0289, 1.0595, 1.1323,  &
1.0006, 1.0012, 1.0027, 1.0041, 1.0070, 1.0118, 1.0216, 1.0419, 1.0879,  &
1.0005, 1.0009, 1.0022, 1.0030, 1.0054, 1.0095, 1.0166, 1.0311, 1.0622,  &
1.0005, 1.0008, 1.0016, 1.0027, 1.0045, 1.0074, 1.0132, 1.0240, 1.0461 /
data ((supersat( 2,iw,iconc,6, 7),iw=1,9),iconc=1,7)/    &
1.0020, 1.0033, 1.0054, 1.0093, 1.0178, 1.0375, 1.0855, 1.2089, 1.5630,  &
1.0017, 1.0021, 1.0040, 1.0068, 1.0122, 1.0239, 1.0514, 1.1192, 1.2972,  &
1.0013, 1.0020, 1.0032, 1.0053, 1.0090, 1.0164, 1.0330, 1.0720, 1.1692,  &
1.0008, 1.0018, 1.0021, 1.0041, 1.0069, 1.0120, 1.0227, 1.0465, 1.1029,  &
1.0006, 1.0014, 1.0020, 1.0034, 1.0056, 1.0094, 1.0168, 1.0322, 1.0671,  &
1.0005, 1.0011, 1.0018, 1.0025, 1.0042, 1.0075, 1.0129, 1.0237, 1.0468,  &
1.0005, 1.0008, 1.0015, 1.0019, 1.0036, 1.0058, 1.0100, 1.0181, 1.0344 /
data ((supersat( 2,iw,iconc,7, 7),iw=1,9),iconc=1,7)/    &
1.0015, 1.0028, 1.0046, 1.0080, 1.0152, 1.0319, 1.0728, 1.1764, 1.4610,  &
1.0011, 1.0022, 1.0034, 1.0058, 1.0104, 1.0202, 1.0433, 1.1002, 1.2468,  &
1.0011, 1.0013, 1.0028, 1.0045, 1.0075, 1.0137, 1.0274, 1.0598, 1.1400,  &
1.0010, 1.0011, 1.0022, 1.0032, 1.0058, 1.0099, 1.0186, 1.0379, 1.0840,  &
1.0009, 1.0011, 1.0014, 1.0029, 1.0046, 1.0077, 1.0135, 1.0258, 1.0537,  &
1.0007, 1.0010, 1.0011, 1.0023, 1.0035, 1.0061, 1.0103, 1.0188, 1.0368,  &
1.0006, 1.0008, 1.0009, 1.0017, 1.0028, 1.0047, 1.0081, 1.0143, 1.0268 /
data ((supersat( 3,iw,iconc,1, 1),iw=1,9),iconc=1,7)/    &
1.0071, 1.0117, 1.0200, 1.0362, 1.0720, 1.1593, 1.4119, 2.7073, 4.9866,  &
1.0054, 1.0091, 1.0157, 1.0278, 1.0523, 1.1078, 1.2490, 1.7178, 4.0203,  &
1.0051, 1.0081, 1.0128, 1.0222, 1.0402, 1.0782, 1.1675, 1.4129, 2.5447,  &
1.0044, 1.0063, 1.0110, 1.0183, 1.0319, 1.0596, 1.1206, 1.2719, 1.7557,  &
1.0032, 1.0053, 1.0088, 1.0151, 1.0257, 1.0469, 1.0909, 1.1929, 1.4713,  &
1.0023, 1.0048, 1.0077, 1.0121, 1.0212, 1.0376, 1.0705, 1.1430, 1.3236,  &
1.0018, 1.0039, 1.0059, 1.0101, 1.0172, 1.0302, 1.0555, 1.1089, 1.2337 /
data ((supersat( 3,iw,iconc,2, 1),iw=1,9),iconc=1,7)/    &
1.0053, 1.0087, 1.0146, 1.0259, 1.0500, 1.1071, 1.2591, 1.8042, 3.7156,  &
1.0040, 1.0069, 1.0116, 1.0199, 1.0365, 1.0727, 1.1610, 1.4126, 2.5898,  &
1.0038, 1.0060, 1.0096, 1.0159, 1.0280, 1.0530, 1.1095, 1.2525, 1.7182,  &
1.0034, 1.0042, 1.0074, 1.0131, 1.0223, 1.0407, 1.0795, 1.1707, 1.4198,  &
1.0026, 1.0039, 1.0067, 1.0105, 1.0184, 1.0323, 1.0605, 1.1231, 1.2781,  &
1.0018, 1.0036, 1.0056, 1.0091, 1.0151, 1.0260, 1.0475, 1.0926, 1.1975,  &
1.0014, 1.0030, 1.0040, 1.0069, 1.0123, 1.0212, 1.0379, 1.0716, 1.1463 /
data ((supersat( 3,iw,iconc,3, 1),iw=1,9),iconc=1,7)/    &
1.0042, 1.0067, 1.0109, 1.0194, 1.0366, 1.0765, 1.1785, 1.4853, 2.7366,  &
1.0027, 1.0054, 1.0085, 1.0148, 1.0266, 1.0518, 1.1116, 1.2695, 1.8156,  &
1.0027, 1.0045, 1.0072, 1.0117, 1.0204, 1.0376, 1.0757, 1.1682, 1.4289,  &
1.0026, 1.0031, 1.0057, 1.0096, 1.0163, 1.0289, 1.0550, 1.1143, 1.2642,  &
1.0022, 1.0027, 1.0052, 1.0081, 1.0134, 1.0230, 1.0420, 1.0828, 1.1788,  &
1.0018, 1.0026, 1.0041, 1.0065, 1.0110, 1.0187, 1.0332, 1.0628, 1.1287,  &
1.0013, 1.0023, 1.0029, 1.0054, 1.0089, 1.0153, 1.0267, 1.0490, 1.0965 /
data ((supersat( 3,iw,iconc,4, 1),iw=1,9),iconc=1,7)/    &
1.0035, 1.0054, 1.0087, 1.0151, 1.0281, 1.0579, 1.1323, 1.3378, 2.0925,  &
1.0025, 1.0042, 1.0067, 1.0114, 1.0202, 1.0388, 1.0821, 1.1925, 1.5186,  &
1.0016, 1.0037, 1.0056, 1.0091, 1.0154, 1.0279, 1.0552, 1.1199, 1.2903,  &
1.0016, 1.0027, 1.0044, 1.0071, 1.0121, 1.0213, 1.0398, 1.0809, 1.1810,  &
1.0016, 1.0017, 1.0040, 1.0063, 1.0098, 1.0168, 1.0303, 1.0584, 1.1225,  &
1.0015, 1.0016, 1.0033, 1.0048, 1.0083, 1.0138, 1.0239, 1.0444, 1.0883,  &
1.0013, 1.0015, 1.0024, 1.0041, 1.0067, 1.0113, 1.0192, 1.0348, 1.0665 /
data ((supersat( 3,iw,iconc,5, 1),iw=1,9),iconc=1,7)/    &
1.0029, 1.0045, 1.0072, 1.0123, 1.0226, 1.0461, 1.1042, 1.2577, 1.7375,  &
1.0024, 1.0032, 1.0054, 1.0092, 1.0161, 1.0305, 1.0639, 1.1476, 1.3772,  &
1.0016, 1.0030, 1.0046, 1.0073, 1.0121, 1.0217, 1.0423, 1.0909, 1.2147,  &
1.0009, 1.0025, 1.0032, 1.0056, 1.0096, 1.0163, 1.0301, 1.0604, 1.1329,  &
1.0006, 1.0018, 1.0031, 1.0050, 1.0077, 1.0130, 1.0228, 1.0431, 1.0889,  &
1.0006, 1.0012, 1.0027, 1.0038, 1.0064, 1.0104, 1.0181, 1.0326, 1.0636,  &
1.0005, 1.0009, 1.0022, 1.0031, 1.0053, 1.0086, 1.0144, 1.0255, 1.0478 /
data ((supersat( 3,iw,iconc,6, 1),iw=1,9),iconc=1,7)/    &
1.0022, 1.0039, 1.0060, 1.0104, 1.0189, 1.0383, 1.0859, 1.2091, 1.5633,  &
1.0020, 1.0030, 1.0046, 1.0077, 1.0133, 1.0250, 1.0521, 1.1195, 1.2975,  &
1.0017, 1.0021, 1.0038, 1.0060, 1.0099, 1.0176, 1.0340, 1.0726, 1.1695,  &
1.0014, 1.0020, 1.0029, 1.0046, 1.0078, 1.0132, 1.0239, 1.0474, 1.1035,  &
1.0009, 1.0018, 1.0021, 1.0040, 1.0062, 1.0103, 1.0179, 1.0334, 1.0680,  &
1.0006, 1.0015, 1.0020, 1.0033, 1.0051, 1.0083, 1.0141, 1.0250, 1.0480,  &
1.0005, 1.0011, 1.0018, 1.0023, 1.0041, 1.0067, 1.0113, 1.0196, 1.0359 /
data ((supersat( 3,iw,iconc,7, 1),iw=1,9),iconc=1,7)/    &
1.0022, 1.0032, 1.0054, 1.0090, 1.0163, 1.0328, 1.0732, 1.1765, 1.4612,  &
1.0015, 1.0027, 1.0041, 1.0065, 1.0113, 1.0212, 1.0440, 1.1005, 1.2471,  &
1.0010, 1.0020, 1.0031, 1.0050, 1.0083, 1.0147, 1.0283, 1.0604, 1.1403,  &
1.0010, 1.0013, 1.0027, 1.0041, 1.0065, 1.0109, 1.0196, 1.0388, 1.0845,  &
1.0010, 1.0011, 1.0021, 1.0031, 1.0051, 1.0085, 1.0146, 1.0269, 1.0546,  &
1.0009, 1.0011, 1.0013, 1.0029, 1.0043, 1.0068, 1.0114, 1.0200, 1.0380,  &
1.0008, 1.0010, 1.0011, 1.0023, 1.0032, 1.0054, 1.0090, 1.0156, 1.0282 /
data ((supersat( 3,iw,iconc,1, 2),iw=1,9),iconc=1,7)/    &
1.0060, 1.0108, 1.0187, 1.0349, 1.0708, 1.1585, 1.4113, 2.7062, 4.9854,  &
1.0052, 1.0087, 1.0149, 1.0267, 1.0511, 1.1067, 1.2482, 1.7168, 4.0178,  &
1.0046, 1.0073, 1.0120, 1.0211, 1.0388, 1.0768, 1.1663, 1.4116, 2.5423,  &
1.0036, 1.0053, 1.0098, 1.0171, 1.0305, 1.0580, 1.1191, 1.2703, 1.7534,  &
1.0025, 1.0049, 1.0082, 1.0139, 1.0244, 1.0451, 1.0890, 1.1909, 1.4686,  &
1.0019, 1.0041, 1.0066, 1.0112, 1.0196, 1.0356, 1.0683, 1.1406, 1.3204,  &
1.0016, 1.0031, 1.0047, 1.0086, 1.0153, 1.0279, 1.0529, 1.1058, 1.2299 /
data ((supersat( 3,iw,iconc,2, 2),iw=1,9),iconc=1,7)/    &
1.0042, 1.0077, 1.0135, 1.0247, 1.0488, 1.1063, 1.2587, 1.8037, 3.7147,  &
1.0039, 1.0065, 1.0105, 1.0187, 1.0352, 1.0717, 1.1603, 1.4119, 2.5893,  &
1.0035, 1.0050, 1.0088, 1.0150, 1.0268, 1.0517, 1.1084, 1.2516, 1.7172,  &
1.0028, 1.0040, 1.0070, 1.0122, 1.0212, 1.0392, 1.0782, 1.1695, 1.4185,  &
1.0019, 1.0037, 1.0061, 1.0099, 1.0170, 1.0307, 1.0589, 1.1215, 1.2764,  &
1.0014, 1.0031, 1.0045, 1.0080, 1.0137, 1.0245, 1.0457, 1.0907, 1.1954,  &
1.0012, 1.0024, 1.0034, 1.0062, 1.0110, 1.0195, 1.0358, 1.0693, 1.1436 /
data ((supersat( 3,iw,iconc,3, 2),iw=1,9),iconc=1,7)/    &
1.0033, 1.0057, 1.0102, 1.0182, 1.0354, 1.0757, 1.1780, 1.4850, 2.7367,  &
1.0027, 1.0050, 1.0080, 1.0137, 1.0254, 1.0507, 1.1109, 1.2690, 1.8150,  &
1.0026, 1.0037, 1.0062, 1.0109, 1.0193, 1.0365, 1.0747, 1.1674, 1.4282,  &
1.0023, 1.0027, 1.0054, 1.0085, 1.0153, 1.0277, 1.0538, 1.1133, 1.2634,  &
1.0018, 1.0026, 1.0045, 1.0074, 1.0124, 1.0217, 1.0407, 1.0815, 1.1776,  &
1.0013, 1.0024, 1.0032, 1.0055, 1.0101, 1.0174, 1.0317, 1.0612, 1.1271,  &
1.0011, 1.0019, 1.0024, 1.0048, 1.0079, 1.0138, 1.0251, 1.0472, 1.0945 /
data ((supersat( 3,iw,iconc,4, 2),iw=1,9),iconc=1,7)/    &
1.0029, 1.0045, 1.0080, 1.0141, 1.0270, 1.0571, 1.1319, 1.3375, 2.0919,  &
1.0017, 1.0039, 1.0063, 1.0106, 1.0192, 1.0377, 1.0814, 1.1921, 1.5182,  &
1.0016, 1.0031, 1.0047, 1.0084, 1.0144, 1.0268, 1.0542, 1.1193, 1.2898,  &
1.0016, 1.0019, 1.0041, 1.0067, 1.0114, 1.0203, 1.0387, 1.0800, 1.1803,  &
1.0015, 1.0016, 1.0035, 1.0055, 1.0091, 1.0159, 1.0292, 1.0572, 1.1215,  &
1.0013, 1.0015, 1.0026, 1.0042, 1.0073, 1.0128, 1.0228, 1.0430, 1.0869,  &
1.0010, 1.0013, 1.0019, 1.0037, 1.0061, 1.0103, 1.0181, 1.0332, 1.0648 /
data ((supersat( 3,iw,iconc,5, 2),iw=1,9),iconc=1,7)/    &
1.0026, 1.0039, 1.0065, 1.0113, 1.0216, 1.0453, 1.1037, 1.2574, 1.7374,  &
1.0019, 1.0030, 1.0051, 1.0085, 1.0151, 1.0295, 1.0632, 1.1472, 1.3769,  &
1.0011, 1.0027, 1.0039, 1.0066, 1.0113, 1.0207, 1.0413, 1.0903, 1.2143,  &
1.0006, 1.0020, 1.0031, 1.0053, 1.0089, 1.0155, 1.0291, 1.0595, 1.1323,  &
1.0006, 1.0013, 1.0028, 1.0043, 1.0072, 1.0120, 1.0218, 1.0420, 1.0880,  &
1.0006, 1.0009, 1.0023, 1.0031, 1.0055, 1.0097, 1.0169, 1.0314, 1.0624,  &
1.0005, 1.0008, 1.0017, 1.0028, 1.0047, 1.0075, 1.0135, 1.0243, 1.0464 /
data ((supersat( 3,iw,iconc,6, 2),iw=1,9),iconc=1,7)/    &
1.0020, 1.0034, 1.0055, 1.0095, 1.0180, 1.0375, 1.0855, 1.2089, 1.5630,  &
1.0018, 1.0023, 1.0041, 1.0070, 1.0124, 1.0241, 1.0514, 1.1191, 1.2974,  &
1.0014, 1.0020, 1.0034, 1.0054, 1.0092, 1.0166, 1.0331, 1.0720, 1.1692,  &
1.0009, 1.0019, 1.0023, 1.0042, 1.0071, 1.0123, 1.0229, 1.0465, 1.1030,  &
1.0006, 1.0015, 1.0020, 1.0036, 1.0057, 1.0096, 1.0170, 1.0324, 1.0672,  &
1.0005, 1.0011, 1.0019, 1.0026, 1.0043, 1.0077, 1.0131, 1.0239, 1.0469,  &
1.0005, 1.0009, 1.0016, 1.0019, 1.0037, 1.0060, 1.0102, 1.0184, 1.0346 /
data ((supersat( 3,iw,iconc,7, 2),iw=1,9),iconc=1,7)/    &
1.0017, 1.0029, 1.0047, 1.0082, 1.0154, 1.0320, 1.0728, 1.1764, 1.4611,  &
1.0011, 1.0023, 1.0035, 1.0060, 1.0105, 1.0203, 1.0434, 1.1002, 1.2469,  &
1.0011, 1.0015, 1.0029, 1.0046, 1.0077, 1.0139, 1.0275, 1.0598, 1.1400,  &
1.0010, 1.0011, 1.0023, 1.0034, 1.0059, 1.0101, 1.0188, 1.0380, 1.0840,  &
1.0009, 1.0011, 1.0015, 1.0030, 1.0047, 1.0079, 1.0137, 1.0260, 1.0538,  &
1.0008, 1.0010, 1.0011, 1.0025, 1.0037, 1.0063, 1.0105, 1.0190, 1.0370,  &
1.0006, 1.0009, 1.0010, 1.0018, 1.0029, 1.0048, 1.0083, 1.0146, 1.0270 /
data ((supersat( 3,iw,iconc,1, 3),iw=1,9),iconc=1,7)/    &
1.0054, 1.0099, 1.0177, 1.0338, 1.0698, 1.1578, 1.4107, 2.7050, 4.9841,  &
1.0049, 1.0081, 1.0140, 1.0255, 1.0499, 1.1057, 1.2473, 1.7157, 4.0149,  &
1.0040, 1.0062, 1.0113, 1.0199, 1.0374, 1.0755, 1.1651, 1.4102, 2.5398,  &
1.0028, 1.0051, 1.0086, 1.0157, 1.0291, 1.0564, 1.1175, 1.2686, 1.7505,  &
1.0020, 1.0043, 1.0073, 1.0125, 1.0229, 1.0433, 1.0871, 1.1886, 1.4654,  &
1.0016, 1.0033, 1.0054, 1.0099, 1.0177, 1.0335, 1.0659, 1.1377, 1.3166,  &
1.0013, 1.0023, 1.0037, 1.0070, 1.0133, 1.0256, 1.0501, 1.1024, 1.2252 /
data ((supersat( 3,iw,iconc,2, 3),iw=1,9),iconc=1,7)/    &
1.0039, 1.0070, 1.0125, 1.0236, 1.0479, 1.1056, 1.2582, 1.8032, 3.7138,  &
1.0036, 1.0059, 1.0098, 1.0178, 1.0341, 1.0708, 1.1596, 1.4113, 2.5879,  &
1.0030, 1.0040, 1.0078, 1.0139, 1.0257, 1.0506, 1.1074, 1.2507, 1.7159,  &
1.0021, 1.0038, 1.0066, 1.0112, 1.0201, 1.0380, 1.0769, 1.1683, 1.4170,  &
1.0015, 1.0033, 1.0052, 1.0091, 1.0157, 1.0294, 1.0574, 1.1199, 1.2746,  &
1.0012, 1.0025, 1.0036, 1.0067, 1.0124, 1.0230, 1.0438, 1.0887, 1.1930,  &
1.0011, 1.0018, 1.0028, 1.0053, 1.0095, 1.0176, 1.0337, 1.0668, 1.1406 /
data ((supersat( 3,iw,iconc,3, 3),iw=1,9),iconc=1,7)/    &
1.0027, 1.0054, 1.0095, 1.0173, 1.0345, 1.0751, 1.1777, 1.4845, 2.7359,  &
1.0026, 1.0044, 1.0074, 1.0130, 1.0244, 1.0498, 1.1102, 1.2685, 1.8145,  &
1.0024, 1.0029, 1.0056, 1.0102, 1.0183, 1.0354, 1.0738, 1.1667, 1.4275,  &
1.0019, 1.0026, 1.0049, 1.0081, 1.0143, 1.0265, 1.0527, 1.1123, 1.2624,  &
1.0014, 1.0024, 1.0038, 1.0065, 1.0112, 1.0206, 1.0394, 1.0802, 1.1762,  &
1.0011, 1.0020, 1.0026, 1.0051, 1.0090, 1.0160, 1.0303, 1.0596, 1.1254,  &
1.0009, 1.0015, 1.0021, 1.0040, 1.0070, 1.0125, 1.0233, 1.0452, 1.0923 /
data ((supersat( 3,iw,iconc,4, 3),iw=1,9),iconc=1,7)/    &
1.0022, 1.0042, 1.0073, 1.0132, 1.0262, 1.0564, 1.1315, 1.3373, 2.0918,  &
1.0016, 1.0035, 1.0057, 1.0099, 1.0183, 1.0369, 1.0809, 1.1917, 1.5179,  &
1.0016, 1.0024, 1.0043, 1.0076, 1.0135, 1.0259, 1.0534, 1.1187, 1.2894,  &
1.0015, 1.0016, 1.0038, 1.0062, 1.0107, 1.0193, 1.0377, 1.0791, 1.1795,  &
1.0013, 1.0015, 1.0029, 1.0047, 1.0085, 1.0150, 1.0280, 1.0561, 1.1205,  &
1.0010, 1.0014, 1.0021, 1.0039, 1.0065, 1.0116, 1.0216, 1.0416, 1.0856,  &
1.0008, 1.0011, 1.0016, 1.0031, 1.0053, 1.0091, 1.0167, 1.0316, 1.0631 /
data ((supersat( 3,iw,iconc,5, 3),iw=1,9),iconc=1,7)/    &
1.0022, 1.0032, 1.0059, 1.0106, 1.0208, 1.0447, 1.1034, 1.2573, 1.7372,  &
1.0013, 1.0028, 1.0046, 1.0077, 1.0143, 1.0287, 1.0627, 1.1469, 1.3767,  &
1.0007, 1.0023, 1.0032, 1.0059, 1.0106, 1.0198, 1.0406, 1.0898, 1.2139,  &
1.0006, 1.0015, 1.0029, 1.0049, 1.0081, 1.0146, 1.0282, 1.0588, 1.1317,  &
1.0006, 1.0010, 1.0025, 1.0036, 1.0065, 1.0113, 1.0208, 1.0411, 1.0871,  &
1.0005, 1.0008, 1.0019, 1.0029, 1.0051, 1.0089, 1.0159, 1.0302, 1.0612,  &
1.0004, 1.0007, 1.0014, 1.0024, 1.0040, 1.0069, 1.0123, 1.0229, 1.0449 /
data ((supersat( 3,iw,iconc,6, 3),iw=1,9),iconc=1,7)/    &
1.0019, 1.0029, 1.0050, 1.0088, 1.0172, 1.0369, 1.0852, 1.2087, 1.5629,  &
1.0015, 1.0020, 1.0038, 1.0063, 1.0117, 1.0233, 1.0509, 1.1189, 1.2972,  &
1.0011, 1.0019, 1.0028, 1.0048, 1.0085, 1.0159, 1.0323, 1.0715, 1.1689,  &
1.0007, 1.0016, 1.0020, 1.0039, 1.0064, 1.0115, 1.0221, 1.0459, 1.1025,  &
1.0005, 1.0012, 1.0019, 1.0031, 1.0052, 1.0088, 1.0161, 1.0315, 1.0665,  &
1.0005, 1.0009, 1.0016, 1.0021, 1.0040, 1.0069, 1.0122, 1.0229, 1.0459,  &
1.0004, 1.0007, 1.0012, 1.0017, 1.0032, 1.0054, 1.0094, 1.0173, 1.0334 /
data ((supersat( 3,iw,iconc,7, 3),iw=1,9),iconc=1,7)/    &
1.0012, 1.0026, 1.0043, 1.0075, 1.0146, 1.0314, 1.0725, 1.1762, 1.4609,  &
1.0011, 1.0018, 1.0031, 1.0054, 1.0099, 1.0196, 1.0429, 1.0999, 1.2467,  &
1.0010, 1.0011, 1.0026, 1.0042, 1.0071, 1.0132, 1.0268, 1.0594, 1.1397,  &
1.0010, 1.0011, 1.0018, 1.0031, 1.0053, 1.0095, 1.0180, 1.0374, 1.0836,  &
1.0008, 1.0010, 1.0012, 1.0027, 1.0043, 1.0072, 1.0130, 1.0252, 1.0531,  &
1.0006, 1.0009, 1.0010, 1.0020, 1.0031, 1.0056, 1.0098, 1.0181, 1.0361,  &
1.0005, 1.0007, 1.0008, 1.0015, 1.0026, 1.0043, 1.0076, 1.0135, 1.0259 /
data ((supersat( 3,iw,iconc,1, 4),iw=1,9),iconc=1,7)/    &
1.0052, 1.0089, 1.0169, 1.0329, 1.0690, 1.1572, 1.4100, 2.7041, 4.9828,  &
1.0045, 1.0075, 1.0131, 1.0244, 1.0489, 1.1048, 1.2464, 1.7145, 4.0118,  &
1.0033, 1.0053, 1.0104, 1.0188, 1.0362, 1.0744, 1.1639, 1.4086, 2.5364,  &
1.0022, 1.0047, 1.0081, 1.0146, 1.0276, 1.0549, 1.1159, 1.2666, 1.7470,  &
1.0017, 1.0036, 1.0063, 1.0112, 1.0213, 1.0416, 1.0851, 1.1862, 1.4613,  &
1.0013, 1.0025, 1.0043, 1.0085, 1.0162, 1.0315, 1.0634, 1.1346, 1.3119,  &
1.0009, 1.0015, 1.0028, 1.0056, 1.0115, 1.0232, 1.0469, 1.0985, 1.2196 /
data ((supersat( 3,iw,iconc,2, 4),iw=1,9),iconc=1,7)/    &
1.0038, 1.0066, 1.0119, 1.0228, 1.0471, 1.1051, 1.2578, 1.8028, 3.7130,  &
1.0033, 1.0052, 1.0093, 1.0170, 1.0333, 1.0700, 1.1590, 1.4105, 2.5866,  &
1.0025, 1.0039, 1.0070, 1.0129, 1.0247, 1.0496, 1.1065, 1.2498, 1.7145,  &
1.0017, 1.0035, 1.0060, 1.0101, 1.0190, 1.0368, 1.0757, 1.1670, 1.4152,  &
1.0013, 1.0027, 1.0044, 1.0081, 1.0148, 1.0280, 1.0559, 1.1183, 1.2724,  &
1.0011, 1.0020, 1.0031, 1.0059, 1.0113, 1.0214, 1.0421, 1.0866, 1.1903,  &
1.0008, 1.0012, 1.0021, 1.0043, 1.0080, 1.0159, 1.0315, 1.0642, 1.1372 /
data ((supersat( 3,iw,iconc,3, 4),iw=1,9),iconc=1,7)/    &
1.0027, 1.0050, 1.0087, 1.0166, 1.0339, 1.0746, 1.1774, 1.4842, 2.7352,  &
1.0025, 1.0038, 1.0068, 1.0123, 1.0237, 1.0492, 1.1097, 1.2681, 1.8138,  &
1.0021, 1.0027, 1.0053, 1.0095, 1.0175, 1.0346, 1.0731, 1.1661, 1.4267,  &
1.0015, 1.0025, 1.0044, 1.0076, 1.0133, 1.0256, 1.0517, 1.1114, 1.2614,  &
1.0011, 1.0021, 1.0031, 1.0056, 1.0104, 1.0196, 1.0381, 1.0790, 1.1749,  &
1.0009, 1.0016, 1.0022, 1.0045, 1.0079, 1.0150, 1.0288, 1.0580, 1.1235,  &
1.0006, 1.0010, 1.0016, 1.0032, 1.0060, 1.0113, 1.0217, 1.0432, 1.0899 /
data ((supersat( 3,iw,iconc,4, 4),iw=1,9),iconc=1,7)/    &
1.0017, 1.0039, 1.0066, 1.0125, 1.0255, 1.0560, 1.1312, 1.3369, 2.0916,  &
1.0016, 1.0030, 1.0051, 1.0091, 1.0176, 1.0363, 1.0804, 1.1914, 1.5175,  &
1.0015, 1.0019, 1.0041, 1.0069, 1.0129, 1.0252, 1.0527, 1.1182, 1.2889,  &
1.0013, 1.0015, 1.0034, 1.0057, 1.0099, 1.0185, 1.0368, 1.0784, 1.1789,  &
1.0011, 1.0014, 1.0024, 1.0041, 1.0078, 1.0140, 1.0270, 1.0551, 1.1195,  &
1.0008, 1.0011, 1.0017, 1.0034, 1.0059, 1.0108, 1.0203, 1.0403, 1.0842,  &
1.0005, 1.0008, 1.0013, 1.0024, 1.0044, 1.0080, 1.0153, 1.0301, 1.0612 /
data ((supersat( 3,iw,iconc,5, 4),iw=1,9),iconc=1,7)/    &
1.0017, 1.0030, 1.0052, 1.0100, 1.0202, 1.0442, 1.1032, 1.2571, 1.7372,  &
1.0010, 1.0026, 1.0041, 1.0072, 1.0137, 1.0281, 1.0623, 1.1467, 1.3763,  &
1.0006, 1.0018, 1.0030, 1.0053, 1.0099, 1.0191, 1.0400, 1.0894, 1.2136,  &
1.0006, 1.0011, 1.0027, 1.0044, 1.0074, 1.0139, 1.0275, 1.0581, 1.1312,  &
1.0005, 1.0008, 1.0021, 1.0030, 1.0058, 1.0106, 1.0199, 1.0402, 1.0863,  &
1.0004, 1.0007, 1.0015, 1.0026, 1.0046, 1.0081, 1.0150, 1.0292, 1.0601,  &
1.0004, 1.0007, 1.0011, 1.0019, 1.0033, 1.0061, 1.0112, 1.0216, 1.0434 /
data ((supersat( 3,iw,iconc,6, 4),iw=1,9),iconc=1,7)/    &
1.0017, 1.0024, 1.0044, 1.0083, 1.0166, 1.0365, 1.0850, 1.2086, 1.5628,  &
1.0012, 1.0020, 1.0035, 1.0058, 1.0111, 1.0228, 1.0506, 1.1187, 1.2969,  &
1.0008, 1.0017, 1.0023, 1.0042, 1.0078, 1.0152, 1.0318, 1.0712, 1.1687,  &
1.0005, 1.0013, 1.0020, 1.0036, 1.0059, 1.0109, 1.0214, 1.0453, 1.1021,  &
1.0005, 1.0010, 1.0017, 1.0026, 1.0045, 1.0081, 1.0153, 1.0307, 1.0658,  &
1.0004, 1.0008, 1.0013, 1.0018, 1.0036, 1.0062, 1.0114, 1.0220, 1.0450,  &
1.0004, 1.0006, 1.0009, 1.0014, 1.0026, 1.0047, 1.0086, 1.0162, 1.0322 /
data ((supersat( 3,iw,iconc,7, 4),iw=1,9),iconc=1,7)/    &
1.0011, 1.0023, 1.0039, 1.0070, 1.0141, 1.0310, 1.0723, 1.1761, 1.4609,  &
1.0011, 1.0014, 1.0029, 1.0049, 1.0093, 1.0191, 1.0425, 1.0998, 1.2466,  &
1.0010, 1.0011, 1.0022, 1.0037, 1.0065, 1.0126, 1.0263, 1.0590, 1.1395,  &
1.0008, 1.0010, 1.0014, 1.0029, 1.0048, 1.0089, 1.0174, 1.0369, 1.0832,  &
1.0006, 1.0009, 1.0010, 1.0023, 1.0038, 1.0066, 1.0123, 1.0245, 1.0526,  &
1.0005, 1.0007, 1.0009, 1.0017, 1.0028, 1.0049, 1.0091, 1.0173, 1.0353,  &
1.0003, 1.0004, 1.0007, 1.0012, 1.0021, 1.0038, 1.0068, 1.0127, 1.0249 /
data ((supersat( 3,iw,iconc,1, 5),iw=1,9),iconc=1,7)/    &
1.0050, 1.0086, 1.0165, 1.0324, 1.0686, 1.1568, 1.4095, 2.7030, 4.9815,  &
1.0042, 1.0071, 1.0125, 1.0238, 1.0484, 1.1043, 1.2459, 1.7137, 4.0094,  &
1.0029, 1.0051, 1.0099, 1.0181, 1.0356, 1.0737, 1.1632, 1.4076, 2.5341,  &
1.0020, 1.0043, 1.0077, 1.0141, 1.0269, 1.0541, 1.1150, 1.2654, 1.7447,  &
1.0015, 1.0032, 1.0057, 1.0106, 1.0204, 1.0405, 1.0838, 1.1846, 1.4586,  &
1.0011, 1.0021, 1.0038, 1.0077, 1.0152, 1.0302, 1.0619, 1.1326, 1.3088,  &
1.0007, 1.0012, 1.0023, 1.0047, 1.0103, 1.0216, 1.0450, 1.0959, 1.2158 /
data ((supersat( 3,iw,iconc,2, 5),iw=1,9),iconc=1,7)/    &
1.0037, 1.0064, 1.0116, 1.0224, 1.0468, 1.1048, 1.2575, 1.8024, 3.7122,  &
1.0031, 1.0048, 1.0090, 1.0166, 1.0328, 1.0696, 1.1586, 1.4101, 2.5855,  &
1.0022, 1.0037, 1.0068, 1.0126, 1.0242, 1.0491, 1.1060, 1.2493, 1.7137,  &
1.0015, 1.0032, 1.0056, 1.0097, 1.0183, 1.0362, 1.0750, 1.1662, 1.4141,  &
1.0012, 1.0024, 1.0039, 1.0076, 1.0142, 1.0272, 1.0550, 1.1173, 1.2711,  &
1.0009, 1.0016, 1.0027, 1.0054, 1.0106, 1.0204, 1.0410, 1.0853, 1.1886,  &
1.0006, 1.0010, 1.0018, 1.0036, 1.0072, 1.0149, 1.0301, 1.0625, 1.1349 /
data ((supersat( 3,iw,iconc,3, 5),iw=1,9),iconc=1,7)/    &
1.0026, 1.0048, 1.0084, 1.0162, 1.0335, 1.0743, 1.1772, 1.4840, 2.7350,  &
1.0024, 1.0034, 1.0065, 1.0120, 1.0233, 1.0488, 1.1095, 1.2678, 1.8135,  &
1.0019, 1.0026, 1.0051, 1.0091, 1.0171, 1.0342, 1.0727, 1.1657, 1.4261,  &
1.0013, 1.0024, 1.0040, 1.0072, 1.0130, 1.0251, 1.0512, 1.1108, 1.2607,  &
1.0010, 1.0019, 1.0027, 1.0052, 1.0100, 1.0189, 1.0375, 1.0783, 1.1740,  &
1.0008, 1.0013, 1.0019, 1.0041, 1.0073, 1.0143, 1.0279, 1.0571, 1.1224,  &
1.0005, 1.0008, 1.0014, 1.0027, 1.0053, 1.0104, 1.0207, 1.0420, 1.0884 /
data ((supersat( 3,iw,iconc,4, 5),iw=1,9),iconc=1,7)/    &
1.0016, 1.0037, 1.0064, 1.0123, 1.0252, 1.0558, 1.1311, 1.3368, 2.0912,  &
1.0016, 1.0027, 1.0048, 1.0088, 1.0172, 1.0360, 1.0802, 1.1913, 1.5174,  &
1.0015, 1.0016, 1.0039, 1.0066, 1.0126, 1.0248, 1.0524, 1.1179, 1.2885,  &
1.0012, 1.0015, 1.0031, 1.0053, 1.0095, 1.0180, 1.0364, 1.0780, 1.1785,  &
1.0009, 1.0013, 1.0021, 1.0039, 1.0073, 1.0135, 1.0265, 1.0545, 1.1189,  &
1.0006, 1.0009, 1.0015, 1.0031, 1.0055, 1.0103, 1.0197, 1.0396, 1.0833,  &
1.0004, 1.0007, 1.0011, 1.0021, 1.0039, 1.0074, 1.0145, 1.0291, 1.0602 /
data ((supersat( 3,iw,iconc,5, 5),iw=1,9),iconc=1,7)/    &
1.0015, 1.0029, 1.0051, 1.0097, 1.0199, 1.0441, 1.1031, 1.2570, 1.7370,  &
1.0008, 1.0024, 1.0038, 1.0070, 1.0134, 1.0279, 1.0621, 1.1465, 1.3762,  &
1.0006, 1.0016, 1.0030, 1.0052, 1.0095, 1.0188, 1.0397, 1.0892, 1.2134,  &
1.0006, 1.0010, 1.0025, 1.0041, 1.0072, 1.0135, 1.0271, 1.0578, 1.1309,  &
1.0005, 1.0008, 1.0018, 1.0029, 1.0054, 1.0101, 1.0194, 1.0398, 1.0859,  &
1.0004, 1.0007, 1.0013, 1.0023, 1.0042, 1.0076, 1.0145, 1.0286, 1.0595,  &
1.0004, 1.0006, 1.0009, 1.0016, 1.0029, 1.0056, 1.0105, 1.0209, 1.0426 /
data ((supersat( 3,iw,iconc,6, 5),iw=1,9),iconc=1,7)/    &
1.0016, 1.0021, 1.0041, 1.0080, 1.0163, 1.0363, 1.0849, 1.2085, 1.5627,  &
1.0011, 1.0019, 1.0032, 1.0056, 1.0108, 1.0225, 1.0504, 1.1186, 1.2968,  &
1.0007, 1.0016, 1.0020, 1.0041, 1.0076, 1.0149, 1.0316, 1.0710, 1.1685,  &
1.0005, 1.0012, 1.0019, 1.0033, 1.0057, 1.0106, 1.0211, 1.0450, 1.1019,  &
1.0005, 1.0009, 1.0016, 1.0023, 1.0042, 1.0078, 1.0149, 1.0303, 1.0655,  &
1.0004, 1.0007, 1.0011, 1.0016, 1.0033, 1.0057, 1.0110, 1.0215, 1.0445,  &
1.0003, 1.0005, 1.0007, 1.0012, 1.0023, 1.0043, 1.0080, 1.0156, 1.0315 /
data ((supersat( 3,iw,iconc,7, 5),iw=1,9),iconc=1,7)/    &
1.0011, 1.0021, 1.0037, 1.0068, 1.0139, 1.0309, 1.0722, 1.1761, 1.4608,  &
1.0010, 1.0012, 1.0028, 1.0047, 1.0090, 1.0189, 1.0424, 1.0997, 1.2465,  &
1.0009, 1.0011, 1.0020, 1.0034, 1.0063, 1.0123, 1.0261, 1.0589, 1.1394,  &
1.0007, 1.0010, 1.0013, 1.0028, 1.0047, 1.0086, 1.0171, 1.0366, 1.0831,  &
1.0006, 1.0008, 1.0010, 1.0021, 1.0035, 1.0063, 1.0120, 1.0242, 1.0523,  &
1.0004, 1.0006, 1.0008, 1.0015, 1.0026, 1.0046, 1.0087, 1.0169, 1.0348,  &
1.0003, 1.0004, 1.0006, 1.0011, 1.0018, 1.0034, 1.0063, 1.0121, 1.0243 /
data ((supersat( 3,iw,iconc,1, 6),iw=1,9),iconc=1,7)/    &
1.0049, 1.0084, 1.0162, 1.0321, 1.0684, 1.1566, 1.4093, 2.7028, 4.9808,  &
1.0039, 1.0068, 1.0122, 1.0235, 1.0480, 1.1040, 1.2456, 1.7132, 4.0079,  &
1.0027, 1.0050, 1.0095, 1.0178, 1.0352, 1.0733, 1.1627, 1.4069, 2.5323,  &
1.0018, 1.0041, 1.0074, 1.0137, 1.0264, 1.0536, 1.1143, 1.2645, 1.7430,  &
1.0014, 1.0029, 1.0053, 1.0102, 1.0198, 1.0398, 1.0829, 1.1834, 1.4565,  &
1.0009, 1.0017, 1.0034, 1.0071, 1.0144, 1.0293, 1.0607, 1.1311, 1.3064,  &
1.0006, 1.0010, 1.0018, 1.0042, 1.0094, 1.0204, 1.0436, 1.0940, 1.2128 /
data ((supersat( 3,iw,iconc,2, 6),iw=1,9),iconc=1,7)/    &
1.0036, 1.0062, 1.0114, 1.0222, 1.0466, 1.1047, 1.2574, 1.8022, 3.7116,  &
1.0029, 1.0045, 1.0087, 1.0163, 1.0326, 1.0694, 1.1584, 1.4098, 2.5849,  &
1.0020, 1.0037, 1.0066, 1.0123, 1.0239, 1.0488, 1.1057, 1.2489, 1.7130,  &
1.0014, 1.0030, 1.0053, 1.0094, 1.0179, 1.0358, 1.0746, 1.1657, 1.4132,  &
1.0011, 1.0022, 1.0036, 1.0072, 1.0137, 1.0267, 1.0544, 1.1165, 1.2700,  &
1.0008, 1.0014, 1.0024, 1.0050, 1.0101, 1.0197, 1.0402, 1.0843, 1.1872,  &
1.0005, 1.0008, 1.0015, 1.0031, 1.0066, 1.0141, 1.0291, 1.0612, 1.1332 /
data ((supersat( 3,iw,iconc,3, 6),iw=1,9),iconc=1,7)/    &
1.0026, 1.0047, 1.0081, 1.0160, 1.0334, 1.0742, 1.1771, 1.4838, 2.7346,  &
1.0023, 1.0032, 1.0062, 1.0117, 1.0231, 1.0486, 1.1093, 1.2677, 1.8134,  &
1.0017, 1.0026, 1.0050, 1.0088, 1.0169, 1.0339, 1.0725, 1.1655, 1.4258,  &
1.0012, 1.0022, 1.0038, 1.0069, 1.0127, 1.0247, 1.0508, 1.1105, 1.2603,  &
1.0009, 1.0017, 1.0025, 1.0050, 1.0097, 1.0185, 1.0370, 1.0777, 1.1733,  &
1.0007, 1.0011, 1.0018, 1.0038, 1.0070, 1.0138, 1.0273, 1.0564, 1.1215,  &
1.0004, 1.0007, 1.0012, 1.0024, 1.0048, 1.0098, 1.0200, 1.0411, 1.0873 /
data ((supersat( 3,iw,iconc,4, 6),iw=1,9),iconc=1,7)/    &
1.0016, 1.0036, 1.0062, 1.0121, 1.0251, 1.0557, 1.1310, 1.3367, 2.0914,  &
1.0015, 1.0025, 1.0046, 1.0086, 1.0170, 1.0358, 1.0801, 1.1911, 1.5171,  &
1.0014, 1.0016, 1.0038, 1.0065, 1.0124, 1.0246, 1.0522, 1.1178, 1.2883,  &
1.0011, 1.0014, 1.0029, 1.0051, 1.0092, 1.0177, 1.0361, 1.0777, 1.1782,  &
1.0008, 1.0012, 1.0020, 1.0037, 1.0070, 1.0132, 1.0261, 1.0542, 1.1184,  &
1.0005, 1.0008, 1.0014, 1.0028, 1.0052, 1.0099, 1.0193, 1.0391, 1.0827,  &
1.0004, 1.0006, 1.0010, 1.0018, 1.0035, 1.0070, 1.0140, 1.0284, 1.0594 /
data ((supersat( 3,iw,iconc,5, 6),iw=1,9),iconc=1,7)/    &
1.0014, 1.0028, 1.0050, 1.0095, 1.0197, 1.0440, 1.1030, 1.2569, 1.7369,  &
1.0007, 1.0022, 1.0036, 1.0068, 1.0132, 1.0277, 1.0620, 1.1465, 1.3761,  &
1.0006, 1.0015, 1.0029, 1.0051, 1.0093, 1.0186, 1.0396, 1.0890, 1.2133,  &
1.0005, 1.0009, 1.0024, 1.0039, 1.0070, 1.0132, 1.0269, 1.0576, 1.1307,  &
1.0004, 1.0008, 1.0017, 1.0028, 1.0051, 1.0098, 1.0191, 1.0394, 1.0856,  &
1.0004, 1.0007, 1.0012, 1.0021, 1.0040, 1.0073, 1.0141, 1.0281, 1.0590,  &
1.0003, 1.0005, 1.0008, 1.0013, 1.0027, 1.0052, 1.0101, 1.0203, 1.0420 /
data ((supersat( 3,iw,iconc,6, 6),iw=1,9),iconc=1,7)/    &
1.0015, 1.0020, 1.0040, 1.0078, 1.0162, 1.0362, 1.0848, 1.2085, 1.5627,  &
1.0010, 1.0019, 1.0031, 1.0055, 1.0106, 1.0224, 1.0503, 1.1185, 1.2967,  &
1.0006, 1.0015, 1.0020, 1.0040, 1.0075, 1.0148, 1.0314, 1.0709, 1.1684,  &
1.0005, 1.0011, 1.0018, 1.0032, 1.0056, 1.0104, 1.0209, 1.0449, 1.1017,  &
1.0004, 1.0008, 1.0015, 1.0021, 1.0040, 1.0076, 1.0147, 1.0301, 1.0652,  &
1.0004, 1.0006, 1.0010, 1.0015, 1.0031, 1.0055, 1.0107, 1.0212, 1.0441,  &
1.0003, 1.0004, 1.0006, 1.0011, 1.0021, 1.0040, 1.0077, 1.0151, 1.0309 /
data ((supersat( 3,iw,iconc,7, 6),iw=1,9),iconc=1,7)/    &
1.0011, 1.0020, 1.0035, 1.0067, 1.0138, 1.0308, 1.0722, 1.1760, 1.4608,  &
1.0010, 1.0011, 1.0027, 1.0046, 1.0089, 1.0187, 1.0423, 1.0996, 1.2464,  &
1.0009, 1.0010, 1.0019, 1.0033, 1.0062, 1.0121, 1.0260, 1.0588, 1.1394,  &
1.0007, 1.0009, 1.0012, 1.0027, 1.0045, 1.0084, 1.0169, 1.0365, 1.0829,  &
1.0005, 1.0007, 1.0009, 1.0020, 1.0033, 1.0061, 1.0117, 1.0240, 1.0521,  &
1.0003, 1.0005, 1.0007, 1.0014, 1.0024, 1.0044, 1.0084, 1.0166, 1.0345,  &
1.0002, 1.0004, 1.0006, 1.0009, 1.0016, 1.0031, 1.0060, 1.0118, 1.0238 /
data ((supersat( 3,iw,iconc,1, 7),iw=1,9),iconc=1,7)/    &
1.0048, 1.0083, 1.0160, 1.0319, 1.0682, 1.1564, 1.4090, 2.7022, 4.9802,  &
1.0037, 1.0065, 1.0119, 1.0233, 1.0478, 1.1037, 1.2453, 1.7126, 4.0065,  &
1.0025, 1.0048, 1.0092, 1.0175, 1.0349, 1.0729, 1.1623, 1.4063, 2.5309,  &
1.0018, 1.0039, 1.0071, 1.0134, 1.0261, 1.0531, 1.1138, 1.2637, 1.7415,  &
1.0013, 1.0026, 1.0050, 1.0098, 1.0193, 1.0392, 1.0822, 1.1824, 1.4548,  &
1.0008, 1.0015, 1.0031, 1.0066, 1.0138, 1.0285, 1.0598, 1.1298, 1.3043,  &
1.0005, 1.0008, 1.0016, 1.0037, 1.0087, 1.0195, 1.0424, 1.0923, 1.2103 /
data ((supersat( 3,iw,iconc,2, 7),iw=1,9),iconc=1,7)/    &
1.0035, 1.0061, 1.0112, 1.0220, 1.0465, 1.1045, 1.2572, 1.8021, 3.7113,  &
1.0027, 1.0043, 1.0086, 1.0161, 1.0324, 1.0692, 1.1582, 1.4095, 2.5844,  &
1.0019, 1.0036, 1.0065, 1.0121, 1.0236, 1.0485, 1.1054, 1.2485, 1.7124,  &
1.0013, 1.0029, 1.0050, 1.0092, 1.0177, 1.0354, 1.0742, 1.1652, 1.4125,  &
1.0010, 1.0020, 1.0034, 1.0068, 1.0133, 1.0263, 1.0539, 1.1159, 1.2691,  &
1.0007, 1.0012, 1.0022, 1.0047, 1.0096, 1.0192, 1.0395, 1.0835, 1.1861,  &
1.0004, 1.0007, 1.0013, 1.0027, 1.0061, 1.0134, 1.0283, 1.0602, 1.1318 /
data ((supersat( 3,iw,iconc,3, 7),iw=1,9),iconc=1,7)/    &
1.0025, 1.0045, 1.0079, 1.0159, 1.0332, 1.0741, 1.1770, 1.4838, 2.7346,  &
1.0022, 1.0030, 1.0060, 1.0115, 1.0229, 1.0485, 1.1091, 1.2675, 1.8131,  &
1.0016, 1.0025, 1.0049, 1.0086, 1.0167, 1.0337, 1.0723, 1.1653, 1.4255,  &
1.0011, 1.0021, 1.0036, 1.0067, 1.0125, 1.0245, 1.0506, 1.1102, 1.2599,  &
1.0009, 1.0016, 1.0024, 1.0048, 1.0094, 1.0182, 1.0367, 1.0773, 1.1728,  &
1.0006, 1.0010, 1.0017, 1.0035, 1.0067, 1.0134, 1.0269, 1.0558, 1.1208,  &
1.0004, 1.0006, 1.0010, 1.0021, 1.0045, 1.0094, 1.0194, 1.0404, 1.0864 /
data ((supersat( 3,iw,iconc,4, 7),iw=1,9),iconc=1,7)/    &
1.0016, 1.0035, 1.0061, 1.0119, 1.0250, 1.0556, 1.1310, 1.3367, 2.0912,  &
1.0015, 1.0024, 1.0044, 1.0085, 1.0168, 1.0357, 1.0800, 1.1911, 1.5170,  &
1.0013, 1.0015, 1.0037, 1.0064, 1.0122, 1.0244, 1.0521, 1.1176, 1.2881,  &
1.0010, 1.0014, 1.0027, 1.0049, 1.0090, 1.0175, 1.0359, 1.0775, 1.1779,  &
1.0007, 1.0011, 1.0018, 1.0036, 1.0068, 1.0129, 1.0258, 1.0538, 1.1181,  &
1.0004, 1.0007, 1.0013, 1.0026, 1.0050, 1.0096, 1.0189, 1.0387, 1.0822,  &
1.0003, 1.0005, 1.0008, 1.0016, 1.0033, 1.0066, 1.0136, 1.0279, 1.0587 /
data ((supersat( 3,iw,iconc,5, 7),iw=1,9),iconc=1,7)/    &
1.0012, 1.0028, 1.0049, 1.0094, 1.0197, 1.0439, 1.1030, 1.2568, 1.7369,  &
1.0007, 1.0021, 1.0035, 1.0067, 1.0130, 1.0276, 1.0619, 1.1464, 1.3760,  &
1.0006, 1.0014, 1.0028, 1.0050, 1.0092, 1.0185, 1.0394, 1.0889, 1.2131,  &
1.0005, 1.0009, 1.0022, 1.0037, 1.0069, 1.0131, 1.0267, 1.0574, 1.1305,  &
1.0004, 1.0007, 1.0016, 1.0027, 1.0049, 1.0096, 1.0189, 1.0392, 1.0853,  &
1.0004, 1.0006, 1.0011, 1.0020, 1.0038, 1.0070, 1.0138, 1.0278, 1.0586,  &
1.0003, 1.0004, 1.0007, 1.0012, 1.0024, 1.0049, 1.0098, 1.0199, 1.0415 /
data ((supersat( 3,iw,iconc,6, 7),iw=1,9),iconc=1,7)/    &
1.0014, 1.0020, 1.0039, 1.0077, 1.0161, 1.0361, 1.0848, 1.2084, 1.5626,  &
1.0009, 1.0018, 1.0030, 1.0054, 1.0105, 1.0223, 1.0502, 1.1185, 1.2967,  &
1.0006, 1.0014, 1.0020, 1.0039, 1.0074, 1.0146, 1.0313, 1.0708, 1.1684,  &
1.0005, 1.0010, 1.0018, 1.0031, 1.0054, 1.0103, 1.0207, 1.0447, 1.1016,  &
1.0004, 1.0008, 1.0014, 1.0020, 1.0038, 1.0074, 1.0145, 1.0299, 1.0650,  &
1.0004, 1.0006, 1.0009, 1.0014, 1.0029, 1.0053, 1.0104, 1.0209, 1.0438,  &
1.0003, 1.0004, 1.0006, 1.0010, 1.0019, 1.0037, 1.0074, 1.0148, 1.0305 /
data ((supersat( 3,iw,iconc,7, 7),iw=1,9),iconc=1,7)/    &
1.0011, 1.0019, 1.0034, 1.0065, 1.0137, 1.0307, 1.0721, 1.1760, 1.4607,  &
1.0010, 1.0011, 1.0026, 1.0045, 1.0088, 1.0187, 1.0423, 1.0996, 1.2463,  &
1.0008, 1.0010, 1.0018, 1.0031, 1.0061, 1.0120, 1.0259, 1.0588, 1.1393,  &
1.0006, 1.0009, 1.0011, 1.0026, 1.0044, 1.0083, 1.0168, 1.0364, 1.0828,  &
1.0005, 1.0007, 1.0009, 1.0019, 1.0032, 1.0060, 1.0116, 1.0238, 1.0519,  &
1.0003, 1.0004, 1.0007, 1.0013, 1.0023, 1.0042, 1.0082, 1.0164, 1.0343,  &
1.0002, 1.0003, 1.0005, 1.0008, 1.0015, 1.0029, 1.0057, 1.0115, 1.0235 /
data ((supersat( 4,iw,iconc,1, 1),iw=1,9),iconc=1,7)/    &
1.0051, 1.0088, 1.0168, 1.0327, 1.0688, 1.1570, 1.4097, 2.7037, 4.9826,  &
1.0044, 1.0074, 1.0129, 1.0242, 1.0487, 1.1046, 1.2462, 1.7141, 4.0109,  &
1.0032, 1.0052, 1.0103, 1.0186, 1.0361, 1.0741, 1.1636, 1.4081, 2.5354,  &
1.0022, 1.0046, 1.0080, 1.0146, 1.0274, 1.0547, 1.1156, 1.2661, 1.7462,  &
1.0017, 1.0036, 1.0063, 1.0111, 1.0211, 1.0413, 1.0847, 1.1856, 1.4602,  &
1.0013, 1.0025, 1.0043, 1.0084, 1.0160, 1.0312, 1.0630, 1.1340, 1.3106,  &
1.0009, 1.0015, 1.0028, 1.0055, 1.0114, 1.0229, 1.0464, 1.0977, 1.2183 /
data ((supersat( 4,iw,iconc,2, 1),iw=1,9),iconc=1,7)/    &
1.0038, 1.0066, 1.0119, 1.0226, 1.0470, 1.1050, 1.2576, 1.8027, 3.7127,  &
1.0033, 1.0051, 1.0093, 1.0169, 1.0331, 1.0698, 1.1588, 1.4103, 2.5863,  &
1.0025, 1.0039, 1.0070, 1.0129, 1.0246, 1.0494, 1.1063, 1.2495, 1.7142,  &
1.0017, 1.0034, 1.0059, 1.0100, 1.0188, 1.0366, 1.0755, 1.1667, 1.4147,  &
1.0013, 1.0027, 1.0043, 1.0081, 1.0147, 1.0279, 1.0557, 1.1180, 1.2718,  &
1.0011, 1.0019, 1.0030, 1.0059, 1.0112, 1.0212, 1.0418, 1.0862, 1.1898,  &
1.0008, 1.0012, 1.0021, 1.0042, 1.0079, 1.0158, 1.0312, 1.0637, 1.1365 /
data ((supersat( 4,iw,iconc,3, 1),iw=1,9),iconc=1,7)/    &
1.0027, 1.0050, 1.0087, 1.0165, 1.0337, 1.0745, 1.1773, 1.4841, 2.7355,  &
1.0025, 1.0037, 1.0068, 1.0122, 1.0235, 1.0490, 1.1096, 1.2680, 1.8138,  &
1.0021, 1.0027, 1.0053, 1.0095, 1.0174, 1.0344, 1.0730, 1.1660, 1.4264,  &
1.0015, 1.0025, 1.0044, 1.0075, 1.0133, 1.0254, 1.0515, 1.1112, 1.2611,  &
1.0011, 1.0021, 1.0030, 1.0055, 1.0104, 1.0194, 1.0380, 1.0788, 1.1746,  &
1.0009, 1.0015, 1.0022, 1.0045, 1.0078, 1.0149, 1.0286, 1.0578, 1.1232,  &
1.0006, 1.0010, 1.0016, 1.0032, 1.0059, 1.0112, 1.0216, 1.0429, 1.0895 /
data ((supersat( 4,iw,iconc,4, 1),iw=1,9),iconc=1,7)/    &
1.0017, 1.0039, 1.0065, 1.0125, 1.0254, 1.0559, 1.1312, 1.3369, 2.0915,  &
1.0016, 1.0030, 1.0051, 1.0091, 1.0175, 1.0362, 1.0803, 1.1913, 1.5174,  &
1.0015, 1.0018, 1.0041, 1.0068, 1.0129, 1.0251, 1.0526, 1.1181, 1.2887,  &
1.0013, 1.0015, 1.0033, 1.0056, 1.0098, 1.0184, 1.0367, 1.0783, 1.1787,  &
1.0011, 1.0014, 1.0024, 1.0041, 1.0077, 1.0139, 1.0269, 1.0549, 1.1193,  &
1.0008, 1.0011, 1.0017, 1.0034, 1.0059, 1.0107, 1.0202, 1.0401, 1.0840,  &
1.0005, 1.0007, 1.0013, 1.0024, 1.0044, 1.0080, 1.0152, 1.0299, 1.0610 /
data ((supersat( 4,iw,iconc,5, 1),iw=1,9),iconc=1,7)/    &
1.0017, 1.0030, 1.0052, 1.0099, 1.0201, 1.0441, 1.1031, 1.2569, 1.7370,  &
1.0010, 1.0025, 1.0041, 1.0071, 1.0136, 1.0280, 1.0622, 1.1466, 1.3762,  &
1.0006, 1.0018, 1.0030, 1.0053, 1.0098, 1.0190, 1.0399, 1.0893, 1.2135,  &
1.0006, 1.0012, 1.0027, 1.0044, 1.0074, 1.0138, 1.0274, 1.0580, 1.1311,  &
1.0005, 1.0008, 1.0021, 1.0030, 1.0057, 1.0105, 1.0198, 1.0401, 1.0862,  &
1.0004, 1.0007, 1.0015, 1.0025, 1.0045, 1.0080, 1.0149, 1.0290, 1.0599,  &
1.0004, 1.0007, 1.0011, 1.0019, 1.0033, 1.0060, 1.0111, 1.0215, 1.0433 /
data ((supersat( 4,iw,iconc,6, 1),iw=1,9),iconc=1,7)/    &
1.0017, 1.0024, 1.0044, 1.0082, 1.0165, 1.0364, 1.0849, 1.2085, 1.5628,  &
1.0012, 1.0020, 1.0034, 1.0058, 1.0110, 1.0227, 1.0505, 1.1186, 1.2968,  &
1.0008, 1.0017, 1.0023, 1.0042, 1.0078, 1.0152, 1.0317, 1.0711, 1.1686,  &
1.0005, 1.0013, 1.0020, 1.0035, 1.0059, 1.0109, 1.0213, 1.0452, 1.1020,  &
1.0005, 1.0010, 1.0017, 1.0026, 1.0045, 1.0080, 1.0152, 1.0306, 1.0657,  &
1.0004, 1.0008, 1.0013, 1.0018, 1.0036, 1.0061, 1.0114, 1.0219, 1.0449,  &
1.0004, 1.0006, 1.0009, 1.0013, 1.0026, 1.0047, 1.0085, 1.0161, 1.0320 /
data ((supersat( 4,iw,iconc,7, 1),iw=1,9),iconc=1,7)/    &
1.0011, 1.0023, 1.0038, 1.0070, 1.0140, 1.0309, 1.0722, 1.1760, 1.4608,  &
1.0011, 1.0014, 1.0029, 1.0049, 1.0092, 1.0190, 1.0425, 1.0997, 1.2464,  &
1.0010, 1.0011, 1.0022, 1.0037, 1.0065, 1.0125, 1.0262, 1.0590, 1.1395,  &
1.0008, 1.0010, 1.0014, 1.0029, 1.0048, 1.0089, 1.0173, 1.0368, 1.0832,  &
1.0006, 1.0009, 1.0010, 1.0023, 1.0038, 1.0065, 1.0123, 1.0245, 1.0525,  &
1.0005, 1.0007, 1.0008, 1.0017, 1.0028, 1.0049, 1.0091, 1.0173, 1.0352,  &
1.0003, 1.0004, 1.0007, 1.0012, 1.0021, 1.0037, 1.0067, 1.0126, 1.0248 /
data ((supersat( 4,iw,iconc,1, 2),iw=1,9),iconc=1,7)/    &
1.0049, 1.0084, 1.0161, 1.0320, 1.0683, 1.1566, 1.4092, 2.7024, 4.9808,  &
1.0039, 1.0067, 1.0121, 1.0235, 1.0479, 1.1039, 1.2454, 1.7129, 4.0075,  &
1.0027, 1.0050, 1.0095, 1.0177, 1.0352, 1.0732, 1.1626, 1.4066, 2.5319,  &
1.0018, 1.0041, 1.0073, 1.0137, 1.0264, 1.0535, 1.1142, 1.2642, 1.7426,  &
1.0014, 1.0029, 1.0053, 1.0101, 1.0197, 1.0397, 1.0828, 1.1832, 1.4560,  &
1.0009, 1.0017, 1.0034, 1.0070, 1.0144, 1.0292, 1.0605, 1.1308, 1.3055,  &
1.0006, 1.0010, 1.0018, 1.0041, 1.0094, 1.0203, 1.0434, 1.0936, 1.2121 /
data ((supersat( 4,iw,iconc,2, 2),iw=1,9),iconc=1,7)/    &
1.0036, 1.0062, 1.0113, 1.0221, 1.0465, 1.1046, 1.2573, 1.8022, 3.7116,  &
1.0029, 1.0045, 1.0087, 1.0163, 1.0325, 1.0693, 1.1583, 1.4097, 2.5848,  &
1.0020, 1.0036, 1.0066, 1.0123, 1.0238, 1.0487, 1.1056, 1.2487, 1.7128,  &
1.0014, 1.0030, 1.0053, 1.0094, 1.0179, 1.0357, 1.0745, 1.1656, 1.4129,  &
1.0011, 1.0022, 1.0036, 1.0071, 1.0137, 1.0266, 1.0543, 1.1164, 1.2696,  &
1.0008, 1.0014, 1.0024, 1.0050, 1.0100, 1.0197, 1.0400, 1.0841, 1.1869,  &
1.0005, 1.0008, 1.0015, 1.0031, 1.0066, 1.0140, 1.0290, 1.0610, 1.1329 /
data ((supersat( 4,iw,iconc,3, 2),iw=1,9),iconc=1,7)/    &
1.0026, 1.0047, 1.0081, 1.0160, 1.0333, 1.0742, 1.1770, 1.4839, 2.7348,  &
1.0023, 1.0032, 1.0062, 1.0117, 1.0230, 1.0486, 1.1093, 1.2675, 1.8133,  &
1.0017, 1.0026, 1.0050, 1.0088, 1.0168, 1.0338, 1.0724, 1.1655, 1.4256,  &
1.0012, 1.0022, 1.0038, 1.0069, 1.0127, 1.0247, 1.0508, 1.1104, 1.2600,  &
1.0009, 1.0017, 1.0025, 1.0050, 1.0097, 1.0184, 1.0370, 1.0777, 1.1732,  &
1.0007, 1.0011, 1.0018, 1.0038, 1.0070, 1.0138, 1.0273, 1.0563, 1.1214,  &
1.0004, 1.0007, 1.0012, 1.0024, 1.0048, 1.0098, 1.0199, 1.0410, 1.0871 /
data ((supersat( 4,iw,iconc,4, 2),iw=1,9),iconc=1,7)/    &
1.0016, 1.0036, 1.0062, 1.0120, 1.0250, 1.0556, 1.1310, 1.3367, 2.0912,  &
1.0015, 1.0025, 1.0045, 1.0086, 1.0170, 1.0358, 1.0801, 1.1911, 1.5171,  &
1.0014, 1.0016, 1.0038, 1.0065, 1.0123, 1.0245, 1.0522, 1.1177, 1.2882,  &
1.0011, 1.0014, 1.0029, 1.0050, 1.0092, 1.0177, 1.0361, 1.0777, 1.1781,  &
1.0008, 1.0012, 1.0020, 1.0037, 1.0070, 1.0132, 1.0261, 1.0541, 1.1184,  &
1.0005, 1.0008, 1.0014, 1.0028, 1.0052, 1.0099, 1.0192, 1.0390, 1.0826,  &
1.0004, 1.0006, 1.0010, 1.0018, 1.0035, 1.0070, 1.0140, 1.0284, 1.0593 /
data ((supersat( 4,iw,iconc,5, 2),iw=1,9),iconc=1,7)/    &
1.0013, 1.0028, 1.0049, 1.0095, 1.0197, 1.0439, 1.1030, 1.2568, 1.7370,  &
1.0007, 1.0022, 1.0036, 1.0068, 1.0132, 1.0277, 1.0620, 1.1464, 1.3760,  &
1.0006, 1.0014, 1.0029, 1.0050, 1.0093, 1.0186, 1.0395, 1.0890, 1.2132,  &
1.0005, 1.0010, 1.0024, 1.0039, 1.0070, 1.0132, 1.0268, 1.0576, 1.1307,  &
1.0004, 1.0008, 1.0017, 1.0027, 1.0051, 1.0098, 1.0191, 1.0394, 1.0855,  &
1.0004, 1.0007, 1.0012, 1.0021, 1.0040, 1.0072, 1.0140, 1.0281, 1.0589,  &
1.0003, 1.0005, 1.0008, 1.0013, 1.0026, 1.0052, 1.0101, 1.0203, 1.0419 /
data ((supersat( 4,iw,iconc,6, 2),iw=1,9),iconc=1,7)/    &
1.0015, 1.0020, 1.0040, 1.0078, 1.0162, 1.0362, 1.0848, 1.2084, 1.5627,  &
1.0010, 1.0019, 1.0031, 1.0055, 1.0106, 1.0224, 1.0503, 1.1185, 1.2967,  &
1.0006, 1.0015, 1.0020, 1.0040, 1.0075, 1.0147, 1.0314, 1.0709, 1.1684,  &
1.0005, 1.0011, 1.0018, 1.0032, 1.0056, 1.0104, 1.0208, 1.0448, 1.1017,  &
1.0004, 1.0008, 1.0014, 1.0021, 1.0040, 1.0076, 1.0146, 1.0300, 1.0652,  &
1.0004, 1.0006, 1.0010, 1.0015, 1.0031, 1.0055, 1.0107, 1.0211, 1.0441,  &
1.0003, 1.0004, 1.0006, 1.0011, 1.0021, 1.0040, 1.0076, 1.0151, 1.0309 /
data ((supersat( 4,iw,iconc,7, 2),iw=1,9),iconc=1,7)/    &
1.0011, 1.0020, 1.0035, 1.0066, 1.0137, 1.0308, 1.0721, 1.1760, 1.4607,  &
1.0010, 1.0011, 1.0027, 1.0046, 1.0089, 1.0187, 1.0423, 1.0996, 1.2464,  &
1.0009, 1.0010, 1.0019, 1.0033, 1.0061, 1.0121, 1.0259, 1.0588, 1.1393,  &
1.0007, 1.0009, 1.0012, 1.0026, 1.0045, 1.0084, 1.0169, 1.0365, 1.0829,  &
1.0005, 1.0007, 1.0009, 1.0020, 1.0033, 1.0061, 1.0117, 1.0240, 1.0521,  &
1.0003, 1.0005, 1.0007, 1.0014, 1.0024, 1.0044, 1.0084, 1.0166, 1.0345,  &
1.0002, 1.0004, 1.0006, 1.0009, 1.0016, 1.0031, 1.0059, 1.0117, 1.0238 /
data ((supersat( 4,iw,iconc,1, 3),iw=1,9),iconc=1,7)/    &
1.0046, 1.0080, 1.0155, 1.0315, 1.0678, 1.1561, 1.4085, 2.7010, 4.9786,  &
1.0034, 1.0061, 1.0114, 1.0228, 1.0472, 1.1032, 1.2445, 1.7115, 4.0030,  &
1.0022, 1.0046, 1.0086, 1.0170, 1.0342, 1.0722, 1.1613, 1.4048, 2.5274,  &
1.0016, 1.0034, 1.0065, 1.0127, 1.0252, 1.0521, 1.1125, 1.2617, 1.7379,  &
1.0011, 1.0021, 1.0043, 1.0090, 1.0183, 1.0379, 1.0806, 1.1800, 1.4506,  &
1.0006, 1.0011, 1.0024, 1.0056, 1.0125, 1.0268, 1.0576, 1.1268, 1.2991,  &
1.0004, 1.0006, 1.0011, 1.0027, 1.0072, 1.0174, 1.0397, 1.0885, 1.2040 /
data ((supersat( 4,iw,iconc,2, 3),iw=1,9),iconc=1,7)/    &
1.0033, 1.0058, 1.0109, 1.0217, 1.0462, 1.1043, 1.2569, 1.8017, 3.7103,  &
1.0025, 1.0039, 1.0082, 1.0157, 1.0320, 1.0688, 1.1578, 1.4089, 2.5831,  &
1.0016, 1.0034, 1.0062, 1.0117, 1.0230, 1.0480, 1.1048, 1.2476, 1.7111,  &
1.0012, 1.0025, 1.0046, 1.0087, 1.0171, 1.0347, 1.0734, 1.1642, 1.4108,  &
1.0009, 1.0017, 1.0030, 1.0062, 1.0126, 1.0253, 1.0528, 1.1145, 1.2668,  &
1.0005, 1.0009, 1.0019, 1.0040, 1.0087, 1.0181, 1.0381, 1.0817, 1.1833,  &
1.0003, 1.0005, 1.0009, 1.0021, 1.0051, 1.0120, 1.0266, 1.0578, 1.1284 /
data ((supersat( 4,iw,iconc,3, 3),iw=1,9),iconc=1,7)/    &
1.0024, 1.0043, 1.0077, 1.0156, 1.0330, 1.0739, 1.1768, 1.4836, 2.7339,  &
1.0020, 1.0027, 1.0057, 1.0112, 1.0226, 1.0482, 1.1089, 1.2671, 1.8126,  &
1.0014, 1.0024, 1.0046, 1.0082, 1.0162, 1.0333, 1.0719, 1.1648, 1.4248,  &
1.0010, 1.0019, 1.0032, 1.0063, 1.0120, 1.0239, 1.0500, 1.1095, 1.2588,  &
1.0007, 1.0013, 1.0021, 1.0044, 1.0089, 1.0174, 1.0359, 1.0764, 1.1716,  &
1.0004, 1.0007, 1.0014, 1.0030, 1.0061, 1.0126, 1.0259, 1.0546, 1.1192,  &
1.0003, 1.0004, 1.0008, 1.0016, 1.0037, 1.0083, 1.0181, 1.0388, 1.0842 /
data ((supersat( 4,iw,iconc,4, 3),iw=1,9),iconc=1,7)/    &
1.0016, 1.0033, 1.0059, 1.0117, 1.0247, 1.0554, 1.1309, 1.3365, 2.0909,  &
1.0014, 1.0021, 1.0041, 1.0083, 1.0166, 1.0355, 1.0798, 1.1908, 1.5168,  &
1.0012, 1.0015, 1.0035, 1.0061, 1.0118, 1.0241, 1.0518, 1.1173, 1.2877,  &
1.0009, 1.0013, 1.0025, 1.0045, 1.0086, 1.0171, 1.0355, 1.0771, 1.1773,  &
1.0006, 1.0009, 1.0016, 1.0033, 1.0063, 1.0124, 1.0252, 1.0532, 1.1173,  &
1.0004, 1.0006, 1.0011, 1.0022, 1.0045, 1.0089, 1.0181, 1.0377, 1.0811,  &
1.0003, 1.0004, 1.0006, 1.0012, 1.0027, 1.0059, 1.0126, 1.0267, 1.0572 /
data ((supersat( 4,iw,iconc,5, 3),iw=1,9),iconc=1,7)/    &
1.0011, 1.0027, 1.0047, 1.0092, 1.0195, 1.0438, 1.1029, 1.2567, 1.7368,  &
1.0006, 1.0019, 1.0032, 1.0065, 1.0128, 1.0274, 1.0618, 1.1462, 1.3758,  &
1.0005, 1.0012, 1.0027, 1.0047, 1.0089, 1.0182, 1.0392, 1.0887, 1.2128,  &
1.0004, 1.0008, 1.0020, 1.0034, 1.0066, 1.0127, 1.0263, 1.0571, 1.1302,  &
1.0004, 1.0007, 1.0014, 1.0024, 1.0047, 1.0091, 1.0184, 1.0387, 1.0848,  &
1.0003, 1.0005, 1.0009, 1.0017, 1.0033, 1.0064, 1.0131, 1.0271, 1.0578,  &
1.0002, 1.0003, 1.0005, 1.0009, 1.0020, 1.0043, 1.0090, 1.0190, 1.0403 /
data ((supersat( 4,iw,iconc,6, 3),iw=1,9),iconc=1,7)/    &
1.0012, 1.0020, 1.0038, 1.0075, 1.0159, 1.0360, 1.0847, 1.2083, 1.5625,  &
1.0008, 1.0017, 1.0028, 1.0052, 1.0103, 1.0221, 1.0501, 1.1183, 1.2965,  &
1.0005, 1.0013, 1.0019, 1.0038, 1.0071, 1.0144, 1.0311, 1.0707, 1.1681,  &
1.0005, 1.0009, 1.0016, 1.0028, 1.0052, 1.0099, 1.0204, 1.0445, 1.1013,  &
1.0004, 1.0007, 1.0012, 1.0018, 1.0036, 1.0071, 1.0141, 1.0295, 1.0646,  &
1.0003, 1.0004, 1.0007, 1.0012, 1.0026, 1.0049, 1.0099, 1.0203, 1.0432,  &
1.0002, 1.0003, 1.0004, 1.0008, 1.0015, 1.0032, 1.0067, 1.0140, 1.0296 /
data ((supersat( 4,iw,iconc,7, 3),iw=1,9),iconc=1,7)/    &
1.0010, 1.0017, 1.0032, 1.0064, 1.0135, 1.0306, 1.0721, 1.1759, 1.4606,  &
1.0009, 1.0010, 1.0024, 1.0043, 1.0086, 1.0185, 1.0422, 1.0995, 1.2462,  &
1.0007, 1.0010, 1.0016, 1.0029, 1.0059, 1.0118, 1.0257, 1.0586, 1.1392,  &
1.0005, 1.0008, 1.0010, 1.0024, 1.0042, 1.0080, 1.0165, 1.0362, 1.0827,  &
1.0004, 1.0005, 1.0008, 1.0017, 1.0029, 1.0057, 1.0112, 1.0235, 1.0516,  &
1.0002, 1.0004, 1.0006, 1.0011, 1.0020, 1.0039, 1.0078, 1.0159, 1.0338,  &
1.0002, 1.0002, 1.0004, 1.0006, 1.0012, 1.0025, 1.0052, 1.0108, 1.0228 /
data ((supersat( 4,iw,iconc,1, 4),iw=1,9),iconc=1,7)/    &
1.0042, 1.0077, 1.0150, 1.0310, 1.0673, 1.1555, 1.4078, 2.6992, 4.9761,  &
1.0029, 1.0054, 1.0109, 1.0222, 1.0465, 1.1024, 1.2435, 1.7097, 3.9971,  &
1.0019, 1.0041, 1.0078, 1.0161, 1.0333, 1.0710, 1.1599, 1.4025, 2.5217,  &
1.0013, 1.0027, 1.0056, 1.0116, 1.0240, 1.0506, 1.1106, 1.2588, 1.7319,  &
1.0007, 1.0014, 1.0034, 1.0077, 1.0167, 1.0358, 1.0780, 1.1762, 1.4438,  &
1.0004, 1.0007, 1.0015, 1.0041, 1.0105, 1.0243, 1.0542, 1.1219, 1.2908,  &
1.0003, 1.0004, 1.0007, 1.0015, 1.0050, 1.0143, 1.0353, 1.0823, 1.1940 /
data ((supersat( 4,iw,iconc,2, 4),iw=1,9),iconc=1,7)/    &
1.0030, 1.0054, 1.0104, 1.0213, 1.0458, 1.1039, 1.2565, 1.8011, 3.7086,  &
1.0021, 1.0036, 1.0077, 1.0151, 1.0315, 1.0683, 1.1572, 1.4081, 2.5810,  &
1.0014, 1.0030, 1.0057, 1.0111, 1.0224, 1.0472, 1.1040, 1.2465, 1.7090,  &
1.0010, 1.0020, 1.0039, 1.0080, 1.0162, 1.0336, 1.0722, 1.1625, 1.4081,  &
1.0006, 1.0011, 1.0024, 1.0052, 1.0114, 1.0239, 1.0511, 1.1123, 1.2634,  &
1.0003, 1.0006, 1.0012, 1.0029, 1.0073, 1.0164, 1.0359, 1.0788, 1.1789,  &
1.0002, 1.0003, 1.0006, 1.0012, 1.0036, 1.0098, 1.0237, 1.0540, 1.1227 /
data ((supersat( 4,iw,iconc,3, 4),iw=1,9),iconc=1,7)/    &
1.0023, 1.0039, 1.0074, 1.0152, 1.0327, 1.0737, 1.1765, 1.4832, 2.7333,  &
1.0017, 1.0025, 1.0052, 1.0107, 1.0221, 1.0478, 1.1085, 1.2666, 1.8119,  &
1.0012, 1.0022, 1.0042, 1.0076, 1.0157, 1.0327, 1.0713, 1.1641, 1.4238,  &
1.0008, 1.0015, 1.0028, 1.0056, 1.0113, 1.0232, 1.0491, 1.1085, 1.2574,  &
1.0005, 1.0008, 1.0017, 1.0038, 1.0080, 1.0164, 1.0347, 1.0750, 1.1695,  &
1.0003, 1.0005, 1.0009, 1.0022, 1.0050, 1.0113, 1.0243, 1.0526, 1.1165,  &
1.0002, 1.0003, 1.0005, 1.0010, 1.0025, 1.0067, 1.0160, 1.0362, 1.0807 /
data ((supersat( 4,iw,iconc,4, 4),iw=1,9),iconc=1,7)/    &
1.0015, 1.0030, 1.0057, 1.0114, 1.0245, 1.0553, 1.1307, 1.3363, 2.0907,  &
1.0013, 1.0018, 1.0039, 1.0080, 1.0162, 1.0352, 1.0795, 1.1905, 1.5163,  &
1.0010, 1.0014, 1.0031, 1.0058, 1.0114, 1.0236, 1.0513, 1.1169, 1.2871,  &
1.0007, 1.0010, 1.0021, 1.0040, 1.0081, 1.0165, 1.0348, 1.0764, 1.1764,  &
1.0004, 1.0007, 1.0013, 1.0028, 1.0056, 1.0117, 1.0244, 1.0522, 1.1160,  &
1.0003, 1.0004, 1.0007, 1.0016, 1.0037, 1.0079, 1.0170, 1.0363, 1.0794,  &
1.0002, 1.0003, 1.0004, 1.0008, 1.0018, 1.0046, 1.0111, 1.0248, 1.0548 /
data ((supersat( 4,iw,iconc,5, 4),iw=1,9),iconc=1,7)/    &
1.0009, 1.0025, 1.0045, 1.0089, 1.0193, 1.0436, 1.1027, 1.2566, 1.7366,  &
1.0006, 1.0017, 1.0029, 1.0062, 1.0125, 1.0272, 1.0616, 1.1460, 1.3756,  &
1.0005, 1.0010, 1.0024, 1.0044, 1.0086, 1.0178, 1.0389, 1.0884, 1.2124,  &
1.0004, 1.0008, 1.0017, 1.0030, 1.0061, 1.0123, 1.0258, 1.0566, 1.1296,  &
1.0003, 1.0006, 1.0011, 1.0020, 1.0042, 1.0085, 1.0177, 1.0379, 1.0839,  &
1.0002, 1.0003, 1.0006, 1.0012, 1.0027, 1.0056, 1.0122, 1.0260, 1.0566,  &
1.0002, 1.0002, 1.0003, 1.0006, 1.0013, 1.0033, 1.0078, 1.0175, 1.0385 /
data ((supersat( 4,iw,iconc,6, 4),iw=1,9),iconc=1,7)/    &
1.0010, 1.0019, 1.0036, 1.0073, 1.0158, 1.0359, 1.0846, 1.2082, 1.5625,  &
1.0006, 1.0015, 1.0025, 1.0050, 1.0100, 1.0219, 1.0500, 1.1182, 1.2964,  &
1.0005, 1.0011, 1.0018, 1.0035, 1.0068, 1.0141, 1.0309, 1.0704, 1.1679,  &
1.0004, 1.0007, 1.0014, 1.0025, 1.0048, 1.0095, 1.0200, 1.0441, 1.1009,  &
1.0003, 1.0005, 1.0009, 1.0016, 1.0032, 1.0066, 1.0135, 1.0289, 1.0640,  &
1.0002, 1.0003, 1.0005, 1.0010, 1.0020, 1.0043, 1.0091, 1.0194, 1.0422,  &
1.0001, 1.0002, 1.0003, 1.0005, 1.0010, 1.0024, 1.0057, 1.0128, 1.0283 /
data ((supersat( 4,iw,iconc,7, 4),iw=1,9),iconc=1,7)/    &
1.0010, 1.0015, 1.0029, 1.0061, 1.0133, 1.0305, 1.0720, 1.1758, 1.4606,  &
1.0008, 1.0010, 1.0022, 1.0041, 1.0084, 1.0183, 1.0420, 1.0994, 1.2461,  &
1.0006, 1.0009, 1.0014, 1.0028, 1.0056, 1.0115, 1.0255, 1.0585, 1.1390,  &
1.0004, 1.0006, 1.0009, 1.0021, 1.0039, 1.0076, 1.0162, 1.0359, 1.0824,  &
1.0003, 1.0004, 1.0007, 1.0014, 1.0025, 1.0052, 1.0107, 1.0230, 1.0511,  &
1.0002, 1.0003, 1.0004, 1.0008, 1.0016, 1.0033, 1.0071, 1.0152, 1.0330,  &
1.0001, 1.0002, 1.0003, 1.0004, 1.0008, 1.0018, 1.0044, 1.0098, 1.0216 /
data ((supersat( 4,iw,iconc,1, 5),iw=1,9),iconc=1,7)/    &
1.0040, 1.0074, 1.0147, 1.0307, 1.0671, 1.1552, 1.4073, 2.6979, 4.9739,  &
1.0026, 1.0051, 1.0106, 1.0218, 1.0461, 1.1019, 1.2428, 1.7085, 3.9931,  &
1.0017, 1.0037, 1.0074, 1.0156, 1.0327, 1.0703, 1.1589, 1.4009, 2.5176,  &
1.0010, 1.0023, 1.0051, 1.0109, 1.0232, 1.0496, 1.1093, 1.2567, 1.7278,  &
1.0005, 1.0011, 1.0028, 1.0068, 1.0157, 1.0346, 1.0763, 1.1736, 1.4390,  &
1.0003, 1.0005, 1.0011, 1.0032, 1.0092, 1.0226, 1.0520, 1.1186, 1.2851,  &
1.0002, 1.0003, 1.0005, 1.0011, 1.0036, 1.0123, 1.0325, 1.0780, 1.1870 /
data ((supersat( 4,iw,iconc,2, 5),iw=1,9),iconc=1,7)/    &
1.0029, 1.0052, 1.0102, 1.0211, 1.0456, 1.1037, 1.2562, 1.8006, 3.7074,  &
1.0019, 1.0035, 1.0074, 1.0148, 1.0312, 1.0680, 1.1568, 1.4075, 2.5795,  &
1.0013, 1.0027, 1.0054, 1.0107, 1.0220, 1.0468, 1.1034, 1.2456, 1.7075,  &
1.0008, 1.0017, 1.0035, 1.0075, 1.0156, 1.0330, 1.0714, 1.1614, 1.4063,  &
1.0004, 1.0009, 1.0020, 1.0047, 1.0106, 1.0231, 1.0501, 1.1108, 1.2611,  &
1.0003, 1.0004, 1.0009, 1.0024, 1.0064, 1.0153, 1.0344, 1.0768, 1.1759,  &
1.0002, 1.0003, 1.0005, 1.0009, 1.0027, 1.0084, 1.0218, 1.0514, 1.1188 /
data ((supersat( 4,iw,iconc,3, 5),iw=1,9),iconc=1,7)/    &
1.0021, 1.0037, 1.0073, 1.0151, 1.0326, 1.0736, 1.1764, 1.4830, 2.7327,  &
1.0015, 1.0025, 1.0051, 1.0105, 1.0219, 1.0476, 1.1083, 1.2663, 1.8114,  &
1.0010, 1.0020, 1.0039, 1.0074, 1.0154, 1.0324, 1.0709, 1.1636, 1.4231,  &
1.0007, 1.0013, 1.0025, 1.0053, 1.0109, 1.0227, 1.0486, 1.1078, 1.2564,  &
1.0004, 1.0007, 1.0015, 1.0034, 1.0075, 1.0158, 1.0340, 1.0740, 1.1682,  &
1.0002, 1.0004, 1.0007, 1.0017, 1.0044, 1.0104, 1.0233, 1.0513, 1.1147,  &
1.0002, 1.0002, 1.0004, 1.0007, 1.0019, 1.0058, 1.0147, 1.0345, 1.0783 /
data ((supersat( 4,iw,iconc,4, 5),iw=1,9),iconc=1,7)/    &
1.0015, 1.0029, 1.0055, 1.0113, 1.0244, 1.0552, 1.1306, 1.3362, 2.0905,  &
1.0012, 1.0017, 1.0038, 1.0078, 1.0161, 1.0350, 1.0794, 1.1903, 1.5161,  &
1.0009, 1.0013, 1.0029, 1.0056, 1.0111, 1.0234, 1.0511, 1.1166, 1.2867,  &
1.0005, 1.0009, 1.0019, 1.0037, 1.0078, 1.0161, 1.0345, 1.0759, 1.1758,  &
1.0003, 1.0006, 1.0012, 1.0025, 1.0052, 1.0112, 1.0238, 1.0515, 1.1151,  &
1.0002, 1.0003, 1.0006, 1.0013, 1.0032, 1.0073, 1.0162, 1.0354, 1.0782,  &
1.0001, 1.0002, 1.0003, 1.0006, 1.0014, 1.0039, 1.0101, 1.0236, 1.0532 /
data ((supersat( 4,iw,iconc,5, 5),iw=1,9),iconc=1,7)/    &
1.0008, 1.0023, 1.0044, 1.0088, 1.0192, 1.0436, 1.1027, 1.2565, 1.7365,  &
1.0005, 1.0015, 1.0028, 1.0060, 1.0124, 1.0270, 1.0615, 1.1459, 1.3754,  &
1.0005, 1.0010, 1.0023, 1.0043, 1.0084, 1.0176, 1.0387, 1.0882, 1.2122,  &
1.0004, 1.0007, 1.0015, 1.0028, 1.0059, 1.0120, 1.0255, 1.0563, 1.1292,  &
1.0003, 1.0005, 1.0009, 1.0018, 1.0039, 1.0081, 1.0173, 1.0375, 1.0833,  &
1.0002, 1.0003, 1.0005, 1.0010, 1.0023, 1.0052, 1.0116, 1.0253, 1.0557,  &
1.0001, 1.0002, 1.0003, 1.0005, 1.0010, 1.0027, 1.0070, 1.0166, 1.0374 /
data ((supersat( 4,iw,iconc,6, 5),iw=1,9),iconc=1,7)/    &
1.0010, 1.0018, 1.0035, 1.0071, 1.0157, 1.0359, 1.0846, 1.2081, 1.5624,  &
1.0006, 1.0014, 1.0023, 1.0048, 1.0099, 1.0218, 1.0499, 1.1181, 1.2963,  &
1.0005, 1.0010, 1.0017, 1.0034, 1.0066, 1.0139, 1.0307, 1.0703, 1.1677,  &
1.0004, 1.0007, 1.0012, 1.0023, 1.0046, 1.0093, 1.0198, 1.0439, 1.1006,  &
1.0002, 1.0004, 1.0007, 1.0014, 1.0029, 1.0063, 1.0132, 1.0285, 1.0636,  &
1.0002, 1.0002, 1.0004, 1.0008, 1.0017, 1.0039, 1.0087, 1.0189, 1.0416,  &
1.0001, 1.0002, 1.0002, 1.0004, 1.0008, 1.0020, 1.0051, 1.0121, 1.0273 /
data ((supersat( 4,iw,iconc,7, 5),iw=1,9),iconc=1,7)/    &
1.0009, 1.0014, 1.0028, 1.0060, 1.0133, 1.0305, 1.0720, 1.1758, 1.4605,  &
1.0007, 1.0010, 1.0021, 1.0040, 1.0083, 1.0183, 1.0420, 1.0993, 1.2460,  &
1.0005, 1.0008, 1.0013, 1.0027, 1.0054, 1.0114, 1.0254, 1.0583, 1.1388,  &
1.0003, 1.0005, 1.0008, 1.0019, 1.0037, 1.0074, 1.0160, 1.0357, 1.0822,  &
1.0002, 1.0003, 1.0006, 1.0012, 1.0023, 1.0050, 1.0104, 1.0227, 1.0508,  &
1.0001, 1.0002, 1.0003, 1.0006, 1.0014, 1.0030, 1.0067, 1.0147, 1.0325,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0006, 1.0015, 1.0039, 1.0092, 1.0209 /
data ((supersat( 4,iw,iconc,1, 6),iw=1,9),iconc=1,7)/    &
1.0038, 1.0073, 1.0145, 1.0305, 1.0668, 1.1549, 1.4069, 2.6969, 4.9721,  &
1.0025, 1.0049, 1.0104, 1.0215, 1.0458, 1.1015, 1.2422, 1.7075, 3.9899,  &
1.0016, 1.0035, 1.0071, 1.0153, 1.0322, 1.0698, 1.1581, 1.3997, 2.5144,  &
1.0008, 1.0020, 1.0047, 1.0104, 1.0226, 1.0489, 1.1082, 1.2552, 1.7245,  &
1.0004, 1.0009, 1.0023, 1.0062, 1.0149, 1.0336, 1.0749, 1.1715, 1.4353,  &
1.0003, 1.0004, 1.0009, 1.0026, 1.0083, 1.0213, 1.0502, 1.1159, 1.2806,  &
1.0002, 1.0003, 1.0005, 1.0009, 1.0028, 1.0108, 1.0303, 1.0747, 1.1815 /
data ((supersat( 4,iw,iconc,2, 6),iw=1,9),iconc=1,7)/    &
1.0027, 1.0051, 1.0101, 1.0209, 1.0455, 1.1036, 1.2560, 1.8003, 3.7066,  &
1.0018, 1.0034, 1.0072, 1.0146, 1.0309, 1.0677, 1.1565, 1.4071, 2.5783,  &
1.0012, 1.0026, 1.0052, 1.0105, 1.0217, 1.0464, 1.1030, 1.2450, 1.7063,  &
1.0007, 1.0015, 1.0033, 1.0072, 1.0152, 1.0325, 1.0708, 1.1605, 1.4049,  &
1.0004, 1.0007, 1.0018, 1.0043, 1.0101, 1.0224, 1.0492, 1.1097, 1.2592,  &
1.0002, 1.0004, 1.0007, 1.0019, 1.0057, 1.0144, 1.0333, 1.0752, 1.1735,  &
1.0002, 1.0002, 1.0004, 1.0007, 1.0021, 1.0074, 1.0204, 1.0494, 1.1158 /
data ((supersat( 4,iw,iconc,3, 6),iw=1,9),iconc=1,7)/    &
1.0021, 1.0036, 1.0072, 1.0150, 1.0325, 1.0735, 1.1763, 1.4829, 2.7323,  &
1.0014, 1.0024, 1.0050, 1.0103, 1.0218, 1.0474, 1.1081, 1.2661, 1.8110,  &
1.0010, 1.0019, 1.0037, 1.0072, 1.0151, 1.0321, 1.0706, 1.1632, 1.4225,  &
1.0006, 1.0011, 1.0023, 1.0050, 1.0106, 1.0224, 1.0482, 1.1073, 1.2556,  &
1.0003, 1.0006, 1.0013, 1.0030, 1.0071, 1.0154, 1.0334, 1.0733, 1.1672,  &
1.0002, 1.0003, 1.0006, 1.0014, 1.0039, 1.0098, 1.0226, 1.0503, 1.1133,  &
1.0001, 1.0002, 1.0003, 1.0006, 1.0015, 1.0050, 1.0137, 1.0332, 1.0764 /
data ((supersat( 4,iw,iconc,4, 6),iw=1,9),iconc=1,7)/    &
1.0014, 1.0028, 1.0054, 1.0112, 1.0243, 1.0551, 1.1305, 1.3361, 2.0903,  &
1.0012, 1.0016, 1.0037, 1.0077, 1.0159, 1.0349, 1.0792, 1.1901, 1.5159,  &
1.0008, 1.0012, 1.0028, 1.0054, 1.0110, 1.0232, 1.0509, 1.1163, 1.2864,  &
1.0004, 1.0008, 1.0018, 1.0035, 1.0075, 1.0159, 1.0342, 1.0756, 1.1754,  &
1.0003, 1.0005, 1.0010, 1.0022, 1.0049, 1.0108, 1.0234, 1.0510, 1.1145,  &
1.0002, 1.0003, 1.0005, 1.0011, 1.0028, 1.0069, 1.0156, 1.0347, 1.0773,  &
1.0001, 1.0002, 1.0003, 1.0005, 1.0011, 1.0034, 1.0094, 1.0226, 1.0519 /
data ((supersat( 4,iw,iconc,5, 6),iw=1,9),iconc=1,7)/    &
1.0008, 1.0023, 1.0043, 1.0087, 1.0191, 1.0435, 1.1026, 1.2564, 1.7364,  &
1.0005, 1.0014, 1.0028, 1.0059, 1.0123, 1.0270, 1.0614, 1.1458, 1.3753,  &
1.0004, 1.0009, 1.0022, 1.0041, 1.0082, 1.0175, 1.0386, 1.0880, 1.2120,  &
1.0004, 1.0007, 1.0014, 1.0026, 1.0057, 1.0118, 1.0253, 1.0560, 1.1289,  &
1.0002, 1.0004, 1.0008, 1.0016, 1.0036, 1.0078, 1.0170, 1.0371, 1.0829,  &
1.0002, 1.0002, 1.0004, 1.0008, 1.0020, 1.0048, 1.0111, 1.0248, 1.0551,  &
1.0001, 1.0002, 1.0002, 1.0004, 1.0008, 1.0023, 1.0065, 1.0159, 1.0365 /
data ((supersat( 4,iw,iconc,6, 6),iw=1,9),iconc=1,7)/    &
1.0009, 1.0018, 1.0034, 1.0071, 1.0156, 1.0358, 1.0845, 1.2081, 1.5623,  &
1.0006, 1.0014, 1.0022, 1.0047, 1.0098, 1.0218, 1.0498, 1.1180, 1.2962,  &
1.0004, 1.0009, 1.0016, 1.0033, 1.0065, 1.0138, 1.0306, 1.0702, 1.1676,  &
1.0003, 1.0006, 1.0011, 1.0022, 1.0044, 1.0091, 1.0196, 1.0437, 1.1004,  &
1.0002, 1.0003, 1.0006, 1.0013, 1.0027, 1.0060, 1.0129, 1.0282, 1.0632,  &
1.0001, 1.0002, 1.0003, 1.0006, 1.0015, 1.0036, 1.0083, 1.0185, 1.0411,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0006, 1.0017, 1.0047, 1.0115, 1.0267 /
data ((supersat( 4,iw,iconc,7, 6),iw=1,9),iconc=1,7)/    &
1.0009, 1.0013, 1.0028, 1.0060, 1.0132, 1.0305, 1.0719, 1.1757, 1.4605,  &
1.0007, 1.0009, 1.0020, 1.0039, 1.0082, 1.0182, 1.0419, 1.0993, 1.2460,  &
1.0005, 1.0007, 1.0012, 1.0026, 1.0053, 1.0113, 1.0253, 1.0583, 1.1388,  &
1.0003, 1.0004, 1.0008, 1.0018, 1.0036, 1.0073, 1.0158, 1.0355, 1.0820,  &
1.0002, 1.0003, 1.0005, 1.0011, 1.0022, 1.0048, 1.0102, 1.0225, 1.0506,  &
1.0001, 1.0002, 1.0003, 1.0005, 1.0012, 1.0028, 1.0064, 1.0144, 1.0321,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0005, 1.0013, 1.0035, 1.0087, 1.0203 /
data ((supersat( 4,iw,iconc,1, 7),iw=1,9),iconc=1,7)/    &
1.0037, 1.0071, 1.0144, 1.0303, 1.0667, 1.1547, 1.4066, 2.6960, 4.9707,  &
1.0024, 1.0047, 1.0102, 1.0213, 1.0455, 1.1012, 1.2418, 1.7066, 3.9872,  &
1.0015, 1.0033, 1.0069, 1.0150, 1.0319, 1.0693, 1.1575, 1.3987, 2.5116,  &
1.0007, 1.0017, 1.0044, 1.0100, 1.0221, 1.0483, 1.1074, 1.2539, 1.7217,  &
1.0004, 1.0007, 1.0020, 1.0057, 1.0143, 1.0328, 1.0738, 1.1698, 1.4320,  &
1.0002, 1.0004, 1.0007, 1.0022, 1.0075, 1.0203, 1.0488, 1.1136, 1.2767,  &
1.0002, 1.0003, 1.0005, 1.0008, 1.0022, 1.0097, 1.0286, 1.0720, 1.1769 /
data ((supersat( 4,iw,iconc,2, 7),iw=1,9),iconc=1,7)/    &
1.0026, 1.0050, 1.0099, 1.0208, 1.0454, 1.1034, 1.2559, 1.8000, 3.7058,  &
1.0017, 1.0033, 1.0070, 1.0144, 1.0308, 1.0675, 1.1562, 1.4067, 2.5773,  &
1.0011, 1.0024, 1.0050, 1.0102, 1.0214, 1.0461, 1.1026, 1.2445, 1.7054,  &
1.0006, 1.0013, 1.0031, 1.0069, 1.0149, 1.0321, 1.0703, 1.1598, 1.4036,  &
1.0003, 1.0006, 1.0015, 1.0040, 1.0097, 1.0219, 1.0485, 1.1087, 1.2577,  &
1.0002, 1.0003, 1.0006, 1.0016, 1.0052, 1.0137, 1.0324, 1.0739, 1.1715,  &
1.0001, 1.0002, 1.0003, 1.0006, 1.0017, 1.0066, 1.0193, 1.0478, 1.1132 /
data ((supersat( 4,iw,iconc,3, 7),iw=1,9),iconc=1,7)/    &
1.0020, 1.0035, 1.0071, 1.0149, 1.0324, 1.0734, 1.1762, 1.4827, 2.7319,  &
1.0014, 1.0023, 1.0049, 1.0102, 1.0216, 1.0473, 1.1079, 1.2659, 1.8106,  &
1.0009, 1.0017, 1.0036, 1.0071, 1.0149, 1.0319, 1.0704, 1.1630, 1.4220,  &
1.0005, 1.0010, 1.0022, 1.0048, 1.0104, 1.0221, 1.0479, 1.1068, 1.2550,  &
1.0003, 1.0005, 1.0012, 1.0028, 1.0067, 1.0150, 1.0329, 1.0727, 1.1663,  &
1.0002, 1.0003, 1.0005, 1.0012, 1.0036, 1.0093, 1.0219, 1.0495, 1.1121,  &
1.0001, 1.0002, 1.0003, 1.0005, 1.0012, 1.0044, 1.0129, 1.0321, 1.0749 /
data ((supersat( 4,iw,iconc,4, 7),iw=1,9),iconc=1,7)/    &
1.0014, 1.0027, 1.0054, 1.0111, 1.0243, 1.0550, 1.1304, 1.3360, 2.0901,  &
1.0011, 1.0015, 1.0037, 1.0076, 1.0159, 1.0348, 1.0791, 1.1900, 1.5157,  &
1.0007, 1.0011, 1.0027, 1.0053, 1.0108, 1.0231, 1.0507, 1.1161, 1.2861,  &
1.0004, 1.0007, 1.0017, 1.0034, 1.0074, 1.0157, 1.0339, 1.0753, 1.1750,  &
1.0002, 1.0004, 1.0009, 1.0020, 1.0047, 1.0106, 1.0230, 1.0506, 1.1139,  &
1.0002, 1.0002, 1.0004, 1.0009, 1.0025, 1.0065, 1.0152, 1.0341, 1.0765,  &
1.0001, 1.0002, 1.0002, 1.0004, 1.0009, 1.0030, 1.0088, 1.0219, 1.0509 /
data ((supersat( 4,iw,iconc,5, 7),iw=1,9),iconc=1,7)/    &
1.0007, 1.0022, 1.0042, 1.0087, 1.0191, 1.0435, 1.1026, 1.2564, 1.7363,  &
1.0005, 1.0014, 1.0027, 1.0058, 1.0122, 1.0269, 1.0613, 1.1457, 1.3752,  &
1.0004, 1.0009, 1.0021, 1.0040, 1.0081, 1.0174, 1.0384, 1.0879, 1.2119,  &
1.0003, 1.0006, 1.0013, 1.0025, 1.0055, 1.0116, 1.0251, 1.0558, 1.1287,  &
1.0002, 1.0003, 1.0007, 1.0015, 1.0034, 1.0076, 1.0167, 1.0368, 1.0825,  &
1.0001, 1.0002, 1.0003, 1.0007, 1.0018, 1.0046, 1.0108, 1.0244, 1.0545,  &
1.0001, 1.0001, 1.0002, 1.0004, 1.0007, 1.0020, 1.0060, 1.0153, 1.0357 /
data ((supersat( 4,iw,iconc,6, 7),iw=1,9),iconc=1,7)/    &
1.0008, 1.0017, 1.0034, 1.0070, 1.0156, 1.0358, 1.0845, 1.2081, 1.5623,  &
1.0005, 1.0013, 1.0022, 1.0047, 1.0098, 1.0217, 1.0498, 1.1180, 1.2961,  &
1.0004, 1.0009, 1.0015, 1.0032, 1.0064, 1.0137, 1.0305, 1.0701, 1.1675,  &
1.0003, 1.0005, 1.0010, 1.0021, 1.0043, 1.0090, 1.0195, 1.0435, 1.1003,  &
1.0002, 1.0003, 1.0006, 1.0012, 1.0026, 1.0058, 1.0127, 1.0280, 1.0630,  &
1.0001, 1.0002, 1.0003, 1.0006, 1.0013, 1.0034, 1.0080, 1.0181, 1.0407,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0006, 1.0014, 1.0043, 1.0110, 1.0261 /
data ((supersat( 4,iw,iconc,7, 7),iw=1,9),iconc=1,7)/    &
1.0009, 1.0012, 1.0027, 1.0059, 1.0132, 1.0304, 1.0719, 1.1757, 1.4605,  &
1.0006, 1.0009, 1.0020, 1.0039, 1.0081, 1.0182, 1.0419, 1.0992, 1.2459,  &
1.0004, 1.0007, 1.0012, 1.0025, 1.0053, 1.0112, 1.0252, 1.0582, 1.1387,  &
1.0002, 1.0004, 1.0008, 1.0017, 1.0035, 1.0072, 1.0157, 1.0354, 1.0819,  &
1.0002, 1.0003, 1.0005, 1.0010, 1.0021, 1.0046, 1.0100, 1.0223, 1.0504,  &
1.0001, 1.0002, 1.0003, 1.0005, 1.0010, 1.0026, 1.0062, 1.0141, 1.0318,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0005, 1.0011, 1.0032, 1.0083, 1.0199 /
data ((supersat( 5,iw,iconc,1, 1),iw=1,9),iconc=1,7)/    &
1.0042, 1.0076, 1.0150, 1.0309, 1.0673, 1.1554, 1.4077, 2.6990, 4.9755,  &
1.0029, 1.0054, 1.0109, 1.0221, 1.0465, 1.1023, 1.2434, 1.7095, 3.9965,  &
1.0019, 1.0041, 1.0078, 1.0161, 1.0332, 1.0710, 1.1597, 1.4023, 2.5209,  &
1.0013, 1.0027, 1.0056, 1.0115, 1.0239, 1.0505, 1.1104, 1.2585, 1.7312,  &
1.0007, 1.0014, 1.0033, 1.0076, 1.0166, 1.0357, 1.0778, 1.1758, 1.4431,  &
1.0004, 1.0007, 1.0015, 1.0041, 1.0104, 1.0241, 1.0540, 1.1214, 1.2899,  &
1.0003, 1.0004, 1.0007, 1.0015, 1.0049, 1.0142, 1.0350, 1.0817, 1.1928 /
data ((supersat( 5,iw,iconc,2, 1),iw=1,9),iconc=1,7)/    &
1.0030, 1.0054, 1.0104, 1.0212, 1.0458, 1.1039, 1.2565, 1.8010, 3.7084,  &
1.0021, 1.0036, 1.0076, 1.0151, 1.0314, 1.0682, 1.1571, 1.4080, 2.5808,  &
1.0014, 1.0030, 1.0057, 1.0111, 1.0223, 1.0472, 1.1039, 1.2463, 1.7088,  &
1.0010, 1.0020, 1.0039, 1.0080, 1.0162, 1.0336, 1.0721, 1.1624, 1.4079,  &
1.0006, 1.0011, 1.0024, 1.0052, 1.0113, 1.0239, 1.0510, 1.1121, 1.2631,  &
1.0003, 1.0006, 1.0012, 1.0029, 1.0072, 1.0163, 1.0358, 1.0786, 1.1785,  &
1.0002, 1.0003, 1.0006, 1.0012, 1.0035, 1.0097, 1.0236, 1.0537, 1.1221 /
data ((supersat( 5,iw,iconc,3, 1),iw=1,9),iconc=1,7)/    &
1.0023, 1.0039, 1.0074, 1.0152, 1.0327, 1.0737, 1.1765, 1.4832, 2.7332,  &
1.0017, 1.0025, 1.0052, 1.0107, 1.0221, 1.0478, 1.1085, 1.2666, 1.8118,  &
1.0012, 1.0021, 1.0041, 1.0076, 1.0157, 1.0327, 1.0713, 1.1640, 1.4237,  &
1.0008, 1.0015, 1.0027, 1.0056, 1.0113, 1.0231, 1.0491, 1.1084, 1.2573,  &
1.0005, 1.0008, 1.0017, 1.0038, 1.0080, 1.0164, 1.0347, 1.0749, 1.1694,  &
1.0003, 1.0005, 1.0009, 1.0021, 1.0050, 1.0112, 1.0243, 1.0525, 1.1163,  &
1.0002, 1.0003, 1.0005, 1.0010, 1.0025, 1.0067, 1.0160, 1.0361, 1.0804 /
data ((supersat( 5,iw,iconc,4, 1),iw=1,9),iconc=1,7)/    &
1.0015, 1.0030, 1.0057, 1.0114, 1.0245, 1.0552, 1.1307, 1.3363, 2.0906,  &
1.0013, 1.0018, 1.0039, 1.0079, 1.0162, 1.0352, 1.0795, 1.1905, 1.5163,  &
1.0010, 1.0014, 1.0031, 1.0058, 1.0114, 1.0236, 1.0513, 1.1168, 1.2871,  &
1.0007, 1.0010, 1.0021, 1.0040, 1.0081, 1.0165, 1.0348, 1.0763, 1.1764,  &
1.0004, 1.0007, 1.0013, 1.0028, 1.0055, 1.0117, 1.0243, 1.0521, 1.1159,  &
1.0003, 1.0004, 1.0007, 1.0016, 1.0036, 1.0079, 1.0169, 1.0363, 1.0793,  &
1.0002, 1.0003, 1.0004, 1.0008, 1.0018, 1.0046, 1.0110, 1.0248, 1.0547 /
data ((supersat( 5,iw,iconc,5, 1),iw=1,9),iconc=1,7)/    &
1.0009, 1.0025, 1.0045, 1.0089, 1.0192, 1.0436, 1.1027, 1.2566, 1.7365,  &
1.0006, 1.0016, 1.0029, 1.0062, 1.0125, 1.0272, 1.0616, 1.1460, 1.3756,  &
1.0005, 1.0010, 1.0024, 1.0044, 1.0086, 1.0178, 1.0389, 1.0884, 1.2124,  &
1.0004, 1.0008, 1.0017, 1.0030, 1.0061, 1.0122, 1.0258, 1.0566, 1.1296,  &
1.0003, 1.0006, 1.0011, 1.0020, 1.0042, 1.0085, 1.0177, 1.0379, 1.0839,  &
1.0002, 1.0003, 1.0006, 1.0012, 1.0027, 1.0056, 1.0122, 1.0260, 1.0565,  &
1.0002, 1.0002, 1.0003, 1.0006, 1.0013, 1.0033, 1.0078, 1.0175, 1.0385 /
data ((supersat( 5,iw,iconc,6, 1),iw=1,9),iconc=1,7)/    &
1.0010, 1.0019, 1.0036, 1.0072, 1.0157, 1.0359, 1.0846, 1.2082, 1.5624,  &
1.0006, 1.0015, 1.0025, 1.0049, 1.0100, 1.0219, 1.0500, 1.1182, 1.2964,  &
1.0005, 1.0011, 1.0018, 1.0035, 1.0068, 1.0141, 1.0309, 1.0704, 1.1679,  &
1.0004, 1.0007, 1.0014, 1.0025, 1.0048, 1.0095, 1.0200, 1.0441, 1.1009,  &
1.0003, 1.0005, 1.0009, 1.0016, 1.0032, 1.0066, 1.0135, 1.0289, 1.0640,  &
1.0002, 1.0003, 1.0005, 1.0010, 1.0020, 1.0043, 1.0091, 1.0194, 1.0422,  &
1.0001, 1.0002, 1.0003, 1.0005, 1.0010, 1.0024, 1.0057, 1.0128, 1.0282 /
data ((supersat( 5,iw,iconc,7, 1),iw=1,9),iconc=1,7)/    &
1.0010, 1.0015, 1.0029, 1.0061, 1.0133, 1.0305, 1.0720, 1.1758, 1.4606,  &
1.0008, 1.0010, 1.0022, 1.0041, 1.0084, 1.0183, 1.0420, 1.0994, 1.2461,  &
1.0006, 1.0009, 1.0014, 1.0028, 1.0056, 1.0115, 1.0255, 1.0584, 1.1390,  &
1.0004, 1.0006, 1.0009, 1.0021, 1.0039, 1.0076, 1.0162, 1.0359, 1.0823,  &
1.0003, 1.0004, 1.0007, 1.0014, 1.0025, 1.0052, 1.0107, 1.0230, 1.0511,  &
1.0002, 1.0003, 1.0004, 1.0008, 1.0016, 1.0033, 1.0071, 1.0151, 1.0330,  &
1.0001, 1.0002, 1.0003, 1.0004, 1.0008, 1.0018, 1.0044, 1.0098, 1.0217 /
data ((supersat( 5,iw,iconc,1, 2),iw=1,9),iconc=1,7)/    &
1.0038, 1.0072, 1.0145, 1.0305, 1.0668, 1.1549, 1.4069, 2.6969, 4.9722,  &
1.0025, 1.0049, 1.0104, 1.0215, 1.0458, 1.1015, 1.2422, 1.7074, 3.9896,  &
1.0016, 1.0035, 1.0071, 1.0153, 1.0322, 1.0698, 1.1581, 1.3996, 2.5141,  &
1.0008, 1.0020, 1.0047, 1.0104, 1.0225, 1.0488, 1.1082, 1.2550, 1.7242,  &
1.0004, 1.0009, 1.0023, 1.0062, 1.0149, 1.0336, 1.0748, 1.1713, 1.4349,  &
1.0003, 1.0004, 1.0009, 1.0026, 1.0083, 1.0213, 1.0501, 1.1156, 1.2802,  &
1.0002, 1.0003, 1.0005, 1.0009, 1.0028, 1.0108, 1.0302, 1.0745, 1.1810 /
data ((supersat( 5,iw,iconc,2, 2),iw=1,9),iconc=1,7)/    &
1.0027, 1.0051, 1.0100, 1.0209, 1.0455, 1.1036, 1.2560, 1.8003, 3.7065,  &
1.0018, 1.0034, 1.0072, 1.0146, 1.0309, 1.0677, 1.1564, 1.4071, 2.5783,  &
1.0012, 1.0026, 1.0052, 1.0104, 1.0216, 1.0464, 1.1029, 1.2450, 1.7063,  &
1.0007, 1.0015, 1.0033, 1.0072, 1.0152, 1.0325, 1.0707, 1.1605, 1.4047,  &
1.0004, 1.0007, 1.0018, 1.0043, 1.0101, 1.0224, 1.0492, 1.1096, 1.2591,  &
1.0002, 1.0004, 1.0007, 1.0019, 1.0057, 1.0144, 1.0333, 1.0751, 1.1733,  &
1.0002, 1.0002, 1.0004, 1.0007, 1.0021, 1.0074, 1.0204, 1.0494, 1.1156 /
data ((supersat( 5,iw,iconc,3, 2),iw=1,9),iconc=1,7)/    &
1.0021, 1.0036, 1.0072, 1.0150, 1.0325, 1.0735, 1.1763, 1.4828, 2.7323,  &
1.0014, 1.0024, 1.0050, 1.0103, 1.0218, 1.0474, 1.1081, 1.2661, 1.8110,  &
1.0010, 1.0019, 1.0037, 1.0072, 1.0151, 1.0321, 1.0706, 1.1632, 1.4225,  &
1.0006, 1.0011, 1.0023, 1.0050, 1.0106, 1.0224, 1.0482, 1.1073, 1.2556,  &
1.0003, 1.0006, 1.0013, 1.0030, 1.0071, 1.0153, 1.0334, 1.0733, 1.1671,  &
1.0002, 1.0003, 1.0006, 1.0014, 1.0039, 1.0098, 1.0226, 1.0503, 1.1132,  &
1.0001, 1.0002, 1.0003, 1.0006, 1.0015, 1.0050, 1.0137, 1.0331, 1.0763 /
data ((supersat( 5,iw,iconc,4, 2),iw=1,9),iconc=1,7)/    &
1.0014, 1.0028, 1.0054, 1.0112, 1.0243, 1.0551, 1.1305, 1.3361, 2.0902,  &
1.0012, 1.0016, 1.0037, 1.0077, 1.0159, 1.0349, 1.0792, 1.1901, 1.5159,  &
1.0008, 1.0012, 1.0028, 1.0054, 1.0110, 1.0232, 1.0509, 1.1163, 1.2864,  &
1.0004, 1.0008, 1.0018, 1.0035, 1.0075, 1.0159, 1.0342, 1.0756, 1.1754,  &
1.0003, 1.0005, 1.0010, 1.0022, 1.0049, 1.0108, 1.0234, 1.0510, 1.1145,  &
1.0002, 1.0003, 1.0005, 1.0011, 1.0028, 1.0069, 1.0156, 1.0347, 1.0772,  &
1.0001, 1.0002, 1.0003, 1.0005, 1.0011, 1.0034, 1.0094, 1.0226, 1.0519 /
data ((supersat( 5,iw,iconc,5, 2),iw=1,9),iconc=1,7)/    &
1.0008, 1.0023, 1.0043, 1.0087, 1.0191, 1.0435, 1.1026, 1.2564, 1.7364,  &
1.0005, 1.0014, 1.0028, 1.0059, 1.0123, 1.0270, 1.0614, 1.1458, 1.3753,  &
1.0004, 1.0009, 1.0022, 1.0041, 1.0082, 1.0175, 1.0386, 1.0880, 1.2120,  &
1.0004, 1.0007, 1.0014, 1.0026, 1.0057, 1.0118, 1.0253, 1.0560, 1.1289,  &
1.0002, 1.0004, 1.0008, 1.0016, 1.0036, 1.0078, 1.0170, 1.0371, 1.0829,  &
1.0002, 1.0002, 1.0004, 1.0008, 1.0020, 1.0048, 1.0111, 1.0248, 1.0551,  &
1.0001, 1.0002, 1.0002, 1.0004, 1.0008, 1.0023, 1.0065, 1.0159, 1.0365 /
data ((supersat( 5,iw,iconc,6, 2),iw=1,9),iconc=1,7)/    &
1.0009, 1.0018, 1.0034, 1.0071, 1.0156, 1.0358, 1.0845, 1.2081, 1.5623,  &
1.0006, 1.0014, 1.0022, 1.0047, 1.0098, 1.0218, 1.0498, 1.1180, 1.2962,  &
1.0004, 1.0009, 1.0016, 1.0033, 1.0065, 1.0138, 1.0306, 1.0702, 1.1676,  &
1.0003, 1.0006, 1.0011, 1.0022, 1.0044, 1.0091, 1.0196, 1.0437, 1.1004,  &
1.0002, 1.0003, 1.0006, 1.0013, 1.0027, 1.0060, 1.0129, 1.0282, 1.0632,  &
1.0001, 1.0002, 1.0003, 1.0006, 1.0015, 1.0036, 1.0083, 1.0185, 1.0412,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0006, 1.0017, 1.0047, 1.0115, 1.0267 /
data ((supersat( 5,iw,iconc,7, 2),iw=1,9),iconc=1,7)/    &
1.0009, 1.0013, 1.0028, 1.0060, 1.0132, 1.0305, 1.0719, 1.1757, 1.4605,  &
1.0007, 1.0009, 1.0020, 1.0039, 1.0082, 1.0182, 1.0419, 1.0993, 1.2460,  &
1.0005, 1.0007, 1.0012, 1.0026, 1.0053, 1.0113, 1.0253, 1.0583, 1.1388,  &
1.0003, 1.0004, 1.0008, 1.0018, 1.0036, 1.0073, 1.0158, 1.0355, 1.0820,  &
1.0002, 1.0003, 1.0005, 1.0011, 1.0022, 1.0048, 1.0102, 1.0225, 1.0506,  &
1.0001, 1.0002, 1.0003, 1.0005, 1.0012, 1.0028, 1.0064, 1.0144, 1.0322,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0005, 1.0013, 1.0035, 1.0087, 1.0204 /
data ((supersat( 5,iw,iconc,1, 3),iw=1,9),iconc=1,7)/    &
1.0034, 1.0068, 1.0140, 1.0300, 1.0663, 1.1542, 1.4059, 2.6940, 4.9676,  &
1.0021, 1.0044, 1.0098, 1.0208, 1.0450, 1.1005, 1.2407, 1.7047, 3.9805,  &
1.0012, 1.0028, 1.0063, 1.0143, 1.0311, 1.0683, 1.1561, 1.3962, 2.5051,  &
1.0005, 1.0013, 1.0037, 1.0091, 1.0210, 1.0469, 1.1054, 1.2507, 1.7150,  &
1.0003, 1.0005, 1.0014, 1.0046, 1.0129, 1.0309, 1.0711, 1.1657, 1.4244,  &
1.0002, 1.0003, 1.0006, 1.0014, 1.0059, 1.0180, 1.0454, 1.1085, 1.2676,  &
1.0002, 1.0002, 1.0004, 1.0007, 1.0015, 1.0072, 1.0247, 1.0658, 1.1660 /
data ((supersat( 5,iw,iconc,2, 3),iw=1,9),iconc=1,7)/    &
1.0024, 1.0047, 1.0097, 1.0206, 1.0451, 1.1031, 1.2555, 1.7994, 3.7039,  &
1.0015, 1.0031, 1.0067, 1.0141, 1.0304, 1.0671, 1.1557, 1.4058, 2.5750,  &
1.0009, 1.0021, 1.0046, 1.0098, 1.0209, 1.0455, 1.1018, 1.2433, 1.7031,  &
1.0004, 1.0010, 1.0026, 1.0063, 1.0142, 1.0312, 1.0691, 1.1582, 1.4007,  &
1.0002, 1.0004, 1.0011, 1.0033, 1.0087, 1.0207, 1.0470, 1.1064, 1.2540,  &
1.0002, 1.0003, 1.0005, 1.0011, 1.0040, 1.0122, 1.0303, 1.0709, 1.1668,  &
1.0001, 1.0002, 1.0003, 1.0005, 1.0011, 1.0049, 1.0167, 1.0440, 1.1073 /
data ((supersat( 5,iw,iconc,3, 3),iw=1,9),iconc=1,7)/    &
1.0019, 1.0033, 1.0069, 1.0147, 1.0322, 1.0732, 1.1759, 1.4824, 2.7311,  &
1.0012, 1.0022, 1.0047, 1.0100, 1.0214, 1.0470, 1.1076, 1.2654, 1.8098,  &
1.0007, 1.0015, 1.0033, 1.0068, 1.0146, 1.0315, 1.0699, 1.1623, 1.4209,  &
1.0004, 1.0008, 1.0019, 1.0044, 1.0098, 1.0215, 1.0471, 1.1059, 1.2535,  &
1.0002, 1.0004, 1.0008, 1.0023, 1.0061, 1.0141, 1.0319, 1.0713, 1.1642,  &
1.0001, 1.0002, 1.0004, 1.0008, 1.0028, 1.0082, 1.0205, 1.0476, 1.1093,  &
1.0001, 1.0001, 1.0002, 1.0004, 1.0009, 1.0032, 1.0112, 1.0296, 1.0712 /
data ((supersat( 5,iw,iconc,4, 3),iw=1,9),iconc=1,7)/    &
1.0013, 1.0025, 1.0052, 1.0110, 1.0241, 1.0549, 1.1303, 1.3358, 2.0898,  &
1.0010, 1.0014, 1.0035, 1.0074, 1.0157, 1.0346, 1.0789, 1.1897, 1.5153,  &
1.0006, 1.0010, 1.0024, 1.0050, 1.0105, 1.0227, 1.0504, 1.1157, 1.2855,  &
1.0003, 1.0006, 1.0014, 1.0031, 1.0070, 1.0152, 1.0334, 1.0746, 1.1741,  &
1.0002, 1.0003, 1.0006, 1.0017, 1.0042, 1.0099, 1.0223, 1.0497, 1.1127,  &
1.0001, 1.0002, 1.0003, 1.0007, 1.0019, 1.0057, 1.0141, 1.0328, 1.0747,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0007, 1.0021, 1.0074, 1.0201, 1.0485 /
data ((supersat( 5,iw,iconc,5, 3),iw=1,9),iconc=1,7)/    &
1.0007, 1.0021, 1.0041, 1.0085, 1.0190, 1.0434, 1.1025, 1.2563, 1.7362,  &
1.0005, 1.0013, 1.0026, 1.0057, 1.0120, 1.0267, 1.0612, 1.1456, 1.3750,  &
1.0004, 1.0008, 1.0019, 1.0038, 1.0079, 1.0171, 1.0382, 1.0876, 1.2115,  &
1.0003, 1.0005, 1.0011, 1.0023, 1.0052, 1.0112, 1.0247, 1.0554, 1.1281,  &
1.0002, 1.0003, 1.0005, 1.0012, 1.0030, 1.0071, 1.0162, 1.0361, 1.0816,  &
1.0001, 1.0002, 1.0003, 1.0005, 1.0014, 1.0040, 1.0100, 1.0234, 1.0533,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0005, 1.0014, 1.0050, 1.0140, 1.0340 /
data ((supersat( 5,iw,iconc,6, 3),iw=1,9),iconc=1,7)/    &
1.0008, 1.0017, 1.0033, 1.0069, 1.0155, 1.0357, 1.0844, 1.2080, 1.5622,  &
1.0005, 1.0012, 1.0021, 1.0045, 1.0096, 1.0216, 1.0497, 1.1178, 1.2960,  &
1.0004, 1.0008, 1.0014, 1.0030, 1.0062, 1.0135, 1.0304, 1.0699, 1.1672,  &
1.0002, 1.0004, 1.0008, 1.0019, 1.0040, 1.0087, 1.0192, 1.0432, 1.0999,  &
1.0001, 1.0002, 1.0004, 1.0010, 1.0023, 1.0054, 1.0122, 1.0275, 1.0624,  &
1.0001, 1.0001, 1.0002, 1.0004, 1.0010, 1.0029, 1.0074, 1.0174, 1.0398,  &
1.0001, 1.0001, 1.0001, 1.0002, 1.0004, 1.0010, 1.0035, 1.0100, 1.0248 /
data ((supersat( 5,iw,iconc,7, 3),iw=1,9),iconc=1,7)/    &
1.0008, 1.0012, 1.0027, 1.0059, 1.0131, 1.0304, 1.0718, 1.1757, 1.4604,  &
1.0006, 1.0008, 1.0019, 1.0038, 1.0080, 1.0181, 1.0418, 1.0992, 1.2458,  &
1.0003, 1.0006, 1.0011, 1.0024, 1.0051, 1.0111, 1.0251, 1.0580, 1.1385,  &
1.0002, 1.0003, 1.0007, 1.0015, 1.0032, 1.0070, 1.0155, 1.0352, 1.0816,  &
1.0001, 1.0002, 1.0004, 1.0008, 1.0018, 1.0043, 1.0097, 1.0219, 1.0499,  &
1.0001, 1.0001, 1.0002, 1.0004, 1.0008, 1.0022, 1.0057, 1.0135, 1.0311,  &
1.0001, 1.0001, 1.0001, 1.0002, 1.0004, 1.0008, 1.0026, 1.0075, 1.0189 /
data ((supersat( 5,iw,iconc,1, 4),iw=1,9),iconc=1,7)/    &
1.0030, 1.0064, 1.0135, 1.0295, 1.0656, 1.1534, 1.4046, 2.6905, 4.9617,  &
1.0017, 1.0039, 1.0091, 1.0200, 1.0440, 1.0992, 1.2389, 1.7012, 3.9687,  &
1.0008, 1.0021, 1.0054, 1.0132, 1.0297, 1.0665, 1.1535, 1.3919, 2.4937,  &
1.0003, 1.0008, 1.0026, 1.0077, 1.0192, 1.0444, 1.1020, 1.2451, 1.7033,  &
1.0002, 1.0003, 1.0008, 1.0031, 1.0105, 1.0278, 1.0666, 1.1585, 1.4110,  &
1.0002, 1.0002, 1.0004, 1.0008, 1.0035, 1.0142, 1.0397, 1.0995, 1.2518,  &
1.0000, 1.0000, 1.0002, 1.0005, 1.0011, 1.0037, 1.0185, 1.0555, 1.1480 /
data ((supersat( 5,iw,iconc,2, 4),iw=1,9),iconc=1,7)/    &
1.0022, 1.0044, 1.0093, 1.0202, 1.0447, 1.1027, 1.2548, 1.7983, 3.7006,  &
1.0013, 1.0027, 1.0062, 1.0135, 1.0298, 1.0663, 1.1547, 1.4043, 2.5708,  &
1.0006, 1.0015, 1.0039, 1.0090, 1.0200, 1.0444, 1.1003, 1.2412, 1.6990,  &
1.0003, 1.0006, 1.0019, 1.0053, 1.0130, 1.0297, 1.0671, 1.1552, 1.3956,  &
1.0002, 1.0003, 1.0006, 1.0022, 1.0072, 1.0187, 1.0442, 1.1024, 1.2475,  &
1.0001, 1.0002, 1.0003, 1.0007, 1.0024, 1.0096, 1.0268, 1.0658, 1.1585,  &
1.0001, 1.0001, 1.0002, 1.0004, 1.0008, 1.0025, 1.0126, 1.0378, 1.0972 /
data ((supersat( 5,iw,iconc,3, 4),iw=1,9),iconc=1,7)/    &
1.0016, 1.0031, 1.0067, 1.0144, 1.0319, 1.0729, 1.1756, 1.4819, 2.7295,  &
1.0010, 1.0019, 1.0044, 1.0096, 1.0209, 1.0465, 1.1070, 1.2646, 1.8084,  &
1.0005, 1.0011, 1.0028, 1.0062, 1.0139, 1.0308, 1.0690, 1.1611, 1.4190,  &
1.0002, 1.0005, 1.0014, 1.0037, 1.0090, 1.0204, 1.0458, 1.1041, 1.2509,  &
1.0001, 1.0002, 1.0005, 1.0016, 1.0049, 1.0127, 1.0301, 1.0689, 1.1606,  &
1.0001, 1.0002, 1.0003, 1.0005, 1.0017, 1.0065, 1.0181, 1.0444, 1.1045,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0005, 1.0017, 1.0084, 1.0255, 1.0650 /
data ((supersat( 5,iw,iconc,4, 4),iw=1,9),iconc=1,7)/    &
1.0012, 1.0024, 1.0050, 1.0108, 1.0239, 1.0547, 1.1301, 1.3356, 2.0892,  &
1.0008, 1.0013, 1.0032, 1.0071, 1.0153, 1.0343, 1.0785, 1.1893, 1.5146,  &
1.0004, 1.0008, 1.0021, 1.0046, 1.0101, 1.0222, 1.0498, 1.1150, 1.2845,  &
1.0002, 1.0004, 1.0010, 1.0026, 1.0063, 1.0145, 1.0325, 1.0735, 1.1726,  &
1.0001, 1.0002, 1.0004, 1.0011, 1.0034, 1.0089, 1.0210, 1.0480, 1.1104,  &
1.0001, 1.0001, 1.0002, 1.0004, 1.0012, 1.0043, 1.0124, 1.0306, 1.0716,  &
1.0001, 1.0001, 1.0001, 1.0002, 1.0004, 1.0012, 1.0054, 1.0172, 1.0443 /
data ((supersat( 5,iw,iconc,5, 4),iw=1,9),iconc=1,7)/    &
1.0006, 1.0019, 1.0039, 1.0084, 1.0188, 1.0432, 1.1023, 1.2561, 1.7359,  &
1.0005, 1.0011, 1.0024, 1.0054, 1.0118, 1.0265, 1.0609, 1.1452, 1.3746,  &
1.0003, 1.0007, 1.0016, 1.0035, 1.0075, 1.0167, 1.0378, 1.0871, 1.2109,  &
1.0002, 1.0003, 1.0007, 1.0020, 1.0047, 1.0107, 1.0241, 1.0546, 1.1271,  &
1.0001, 1.0002, 1.0003, 1.0008, 1.0024, 1.0063, 1.0152, 1.0350, 1.0801,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0008, 1.0029, 1.0087, 1.0217, 1.0511,  &
1.0001, 1.0001, 1.0001, 1.0002, 1.0003, 1.0008, 1.0035, 1.0118, 1.0310 /
data ((supersat( 5,iw,iconc,6, 4),iw=1,9),iconc=1,7)/    &
1.0007, 1.0015, 1.0032, 1.0068, 1.0154, 1.0356, 1.0843, 1.2079, 1.5620,  &
1.0005, 1.0010, 1.0019, 1.0043, 1.0095, 1.0214, 1.0495, 1.1176, 1.2957,  &
1.0003, 1.0006, 1.0011, 1.0027, 1.0059, 1.0132, 1.0300, 1.0695, 1.1668,  &
1.0002, 1.0003, 1.0006, 1.0015, 1.0036, 1.0082, 1.0187, 1.0426, 1.0992,  &
1.0001, 1.0002, 1.0003, 1.0006, 1.0018, 1.0047, 1.0115, 1.0266, 1.0613,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0006, 1.0021, 1.0063, 1.0161, 1.0382,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0003, 1.0006, 1.0023, 1.0083, 1.0225 /
data ((supersat( 5,iw,iconc,7, 4),iw=1,9),iconc=1,7)/    &
1.0007, 1.0011, 1.0025, 1.0058, 1.0131, 1.0303, 1.0718, 1.1756, 1.4603,  &
1.0004, 1.0007, 1.0017, 1.0036, 1.0079, 1.0179, 1.0416, 1.0990, 1.2456,  &
1.0002, 1.0004, 1.0010, 1.0021, 1.0048, 1.0108, 1.0248, 1.0578, 1.1382,  &
1.0001, 1.0002, 1.0005, 1.0012, 1.0029, 1.0066, 1.0151, 1.0347, 1.0811,  &
1.0001, 1.0001, 1.0002, 1.0005, 1.0014, 1.0037, 1.0090, 1.0212, 1.0491,  &
1.0001, 1.0001, 1.0001, 1.0002, 1.0005, 1.0016, 1.0048, 1.0124, 1.0299,  &
1.0000, 1.0001, 1.0001, 1.0001, 1.0002, 1.0005, 1.0016, 1.0061, 1.0171 /
data ((supersat( 5,iw,iconc,1, 5),iw=1,9),iconc=1,7)/    &
1.0028, 1.0061, 1.0132, 1.0291, 1.0652, 1.1528, 1.4037, 2.6879, 4.9574,  &
1.0015, 1.0036, 1.0086, 1.0195, 1.0434, 1.0984, 1.2377, 1.6988, 3.9603,  &
1.0006, 1.0017, 1.0049, 1.0125, 1.0288, 1.0654, 1.1518, 1.3888, 2.4856,  &
1.0003, 1.0006, 1.0020, 1.0068, 1.0179, 1.0428, 1.0996, 1.2413, 1.6951,  &
1.0002, 1.0003, 1.0006, 1.0022, 1.0090, 1.0257, 1.0635, 1.1536, 1.4016,  &
1.0001, 1.0002, 1.0004, 1.0007, 1.0022, 1.0118, 1.0361, 1.0935, 1.2409,  &
0.9998, 0.9998, 1.0000, 1.0003, 1.0009, 1.0024, 1.0148, 1.0490, 1.1360 /
data ((supersat( 5,iw,iconc,2, 5),iw=1,9),iconc=1,7)/    &
1.0020, 1.0042, 1.0091, 1.0200, 1.0444, 1.1023, 1.2544, 1.7975, 3.6982,  &
1.0011, 1.0025, 1.0059, 1.0132, 1.0294, 1.0658, 1.1540, 1.4032, 2.5679,  &
1.0005, 1.0013, 1.0035, 1.0085, 1.0194, 1.0437, 1.0994, 1.2397, 1.6961,  &
1.0002, 1.0005, 1.0015, 1.0047, 1.0121, 1.0286, 1.0657, 1.1532, 1.3921,  &
1.0001, 1.0002, 1.0005, 1.0016, 1.0062, 1.0173, 1.0424, 1.0997, 1.2429,  &
1.0001, 1.0002, 1.0003, 1.0005, 1.0016, 1.0080, 1.0245, 1.0623, 1.1528,  &
1.0000, 1.0000, 1.0001, 1.0003, 1.0007, 1.0017, 1.0101, 1.0337, 1.0904 /
data ((supersat( 5,iw,iconc,3, 5),iw=1,9),iconc=1,7)/    &
1.0015, 1.0030, 1.0065, 1.0143, 1.0318, 1.0727, 1.1753, 1.4815, 2.7284,  &
1.0009, 1.0017, 1.0042, 1.0093, 1.0207, 1.0462, 1.1066, 1.2641, 1.8074,  &
1.0004, 1.0009, 1.0025, 1.0059, 1.0135, 1.0303, 1.0684, 1.1603, 1.4176,  &
1.0002, 1.0004, 1.0011, 1.0033, 1.0084, 1.0197, 1.0449, 1.1029, 1.2490,  &
1.0001, 1.0002, 1.0004, 1.0011, 1.0042, 1.0118, 1.0288, 1.0672, 1.1580,  &
1.0001, 1.0001, 1.0002, 1.0004, 1.0012, 1.0053, 1.0165, 1.0421, 1.1011,  &
1.0000, 1.0001, 1.0001, 1.0003, 1.0005, 1.0012, 1.0066, 1.0228, 1.0608 /
data ((supersat( 5,iw,iconc,4, 5),iw=1,9),iconc=1,7)/    &
1.0011, 1.0022, 1.0049, 1.0107, 1.0238, 1.0546, 1.1299, 1.3353, 2.0887,  &
1.0006, 1.0012, 1.0031, 1.0069, 1.0151, 1.0341, 1.0783, 1.1890, 1.5141,  &
1.0003, 1.0007, 1.0019, 1.0043, 1.0098, 1.0219, 1.0494, 1.1145, 1.2838,  &
1.0002, 1.0003, 1.0008, 1.0023, 1.0059, 1.0140, 1.0319, 1.0728, 1.1715,  &
1.0001, 1.0002, 1.0003, 1.0008, 1.0029, 1.0082, 1.0201, 1.0469, 1.1089,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0008, 1.0035, 1.0113, 1.0290, 1.0694,  &
1.0001, 1.0001, 1.0001, 1.0002, 1.0004, 1.0009, 1.0041, 1.0153, 1.0416 /
data ((supersat( 5,iw,iconc,5, 5),iw=1,9),iconc=1,7)/    &
1.0006, 1.0018, 1.0038, 1.0083, 1.0187, 1.0431, 1.1022, 1.2560, 1.7357,  &
1.0004, 1.0011, 1.0023, 1.0053, 1.0117, 1.0263, 1.0607, 1.1450, 1.3743,  &
1.0003, 1.0006, 1.0014, 1.0033, 1.0073, 1.0165, 1.0375, 1.0868, 1.2104,  &
1.0001, 1.0003, 1.0006, 1.0017, 1.0044, 1.0103, 1.0236, 1.0541, 1.1264,  &
1.0001, 1.0001, 1.0003, 1.0006, 1.0020, 1.0058, 1.0145, 1.0342, 1.0791,  &
1.0001, 1.0001, 1.0001, 1.0003, 1.0006, 1.0023, 1.0078, 1.0206, 1.0496,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0003, 1.0006, 1.0026, 1.0103, 1.0291 /
data ((supersat( 5,iw,iconc,6, 5),iw=1,9),iconc=1,7)/    &
1.0007, 1.0015, 1.0031, 1.0067, 1.0153, 1.0356, 1.0842, 1.2078, 1.5619,  &
1.0004, 1.0009, 1.0018, 1.0042, 1.0093, 1.0213, 1.0494, 1.1175, 1.2955,  &
1.0002, 1.0005, 1.0010, 1.0025, 1.0057, 1.0130, 1.0298, 1.0693, 1.1665,  &
1.0001, 1.0002, 1.0005, 1.0013, 1.0033, 1.0079, 1.0183, 1.0422, 1.0987,  &
1.0001, 1.0001, 1.0002, 1.0005, 1.0015, 1.0043, 1.0109, 1.0260, 1.0605,  &
1.0001, 1.0001, 1.0001, 1.0002, 1.0005, 1.0016, 1.0056, 1.0152, 1.0371,  &
1.0000, 1.0001, 1.0001, 1.0001, 1.0002, 1.0005, 1.0017, 1.0072, 1.0210 /
data ((supersat( 5,iw,iconc,7, 5),iw=1,9),iconc=1,7)/    &
1.0006, 1.0010, 1.0025, 1.0057, 1.0130, 1.0302, 1.0717, 1.1755, 1.4602,  &
1.0004, 1.0006, 1.0016, 1.0035, 1.0078, 1.0178, 1.0415, 1.0989, 1.2455,  &
1.0002, 1.0004, 1.0009, 1.0020, 1.0047, 1.0107, 1.0247, 1.0576, 1.1380,  &
1.0001, 1.0002, 1.0004, 1.0010, 1.0026, 1.0063, 1.0148, 1.0344, 1.0807,  &
1.0001, 1.0001, 1.0002, 1.0004, 1.0012, 1.0034, 1.0086, 1.0207, 1.0485,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0004, 1.0012, 1.0042, 1.0117, 1.0290,  &
1.0000, 1.0001, 1.0001, 1.0001, 1.0002, 1.0004, 1.0012, 1.0052, 1.0158 /
data ((supersat( 5,iw,iconc,1, 6),iw=1,9),iconc=1,7)/    &
1.0026, 1.0059, 1.0130, 1.0288, 1.0649, 1.1523, 1.4029, 2.6858, 4.9540,  &
1.0013, 1.0034, 1.0083, 1.0191, 1.0429, 1.0977, 1.2367, 1.6969, 3.9536,  &
1.0005, 1.0015, 1.0046, 1.0119, 1.0281, 1.0644, 1.1504, 1.3864, 2.4791,  &
1.0002, 1.0005, 1.0016, 1.0061, 1.0170, 1.0415, 1.0977, 1.2382, 1.6885,  &
1.0002, 1.0003, 1.0005, 1.0016, 1.0079, 1.0241, 1.0611, 1.1497, 1.3942,  &
1.0000, 1.0001, 1.0003, 1.0007, 1.0017, 1.0101, 1.0333, 1.0889, 1.2324,  &
0.9996, 0.9996, 0.9998, 1.0001, 1.0007, 1.0022, 1.0123, 1.0442, 1.1269 /
data ((supersat( 5,iw,iconc,2, 6),iw=1,9),iconc=1,7)/    &
1.0019, 1.0041, 1.0090, 1.0198, 1.0442, 1.1021, 1.2540, 1.7969, 3.6963,  &
1.0010, 1.0023, 1.0057, 1.0129, 1.0291, 1.0655, 1.1535, 1.4024, 2.5655,  &
1.0004, 1.0011, 1.0032, 1.0082, 1.0190, 1.0431, 1.0986, 1.2385, 1.6937,  &
1.0002, 1.0004, 1.0012, 1.0042, 1.0115, 1.0279, 1.0646, 1.1516, 1.3892,  &
1.0001, 1.0002, 1.0004, 1.0012, 1.0054, 1.0163, 1.0409, 1.0975, 1.2393,  &
1.0001, 1.0001, 1.0003, 1.0005, 1.0012, 1.0068, 1.0227, 1.0595, 1.1482,  &
0.9999, 0.9999, 1.0000, 1.0002, 1.0006, 1.0014, 1.0083, 1.0307, 1.0852 /
data ((supersat( 5,iw,iconc,3, 6),iw=1,9),iconc=1,7)/    &
1.0014, 1.0029, 1.0064, 1.0142, 1.0316, 1.0725, 1.1751, 1.4812, 2.7276,  &
1.0008, 1.0016, 1.0040, 1.0091, 1.0205, 1.0460, 1.1063, 1.2636, 1.8065,  &
1.0003, 1.0008, 1.0023, 1.0056, 1.0132, 1.0299, 1.0679, 1.1596, 1.4166,  &
1.0002, 1.0003, 1.0009, 1.0029, 1.0080, 1.0192, 1.0442, 1.1020, 1.2475,  &
1.0001, 1.0002, 1.0003, 1.0009, 1.0037, 1.0110, 1.0279, 1.0659, 1.1560,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0009, 1.0045, 1.0153, 1.0404, 1.0984,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0004, 1.0010, 1.0053, 1.0208, 1.0576 /
data ((supersat( 5,iw,iconc,4, 6),iw=1,9),iconc=1,7)/    &
1.0010, 1.0022, 1.0048, 1.0106, 1.0237, 1.0545, 1.1298, 1.3352, 2.0884,  &
1.0005, 1.0012, 1.0030, 1.0067, 1.0150, 1.0339, 1.0781, 1.1887, 1.5137,  &
1.0003, 1.0006, 1.0017, 1.0041, 1.0095, 1.0216, 1.0491, 1.1141, 1.2832,  &
1.0001, 1.0003, 1.0007, 1.0021, 1.0055, 1.0136, 1.0314, 1.0721, 1.1706,  &
1.0001, 1.0001, 1.0003, 1.0007, 1.0025, 1.0076, 1.0195, 1.0461, 1.1076,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0007, 1.0029, 1.0104, 1.0279, 1.0677,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0003, 1.0007, 1.0033, 1.0139, 1.0394 /
data ((supersat( 5,iw,iconc,5, 6),iw=1,9),iconc=1,7)/    &
1.0006, 1.0018, 1.0037, 1.0082, 1.0187, 1.0431, 1.1021, 1.2559, 1.7355,  &
1.0004, 1.0010, 1.0022, 1.0052, 1.0115, 1.0262, 1.0606, 1.1449, 1.3740,  &
1.0002, 1.0005, 1.0012, 1.0031, 1.0071, 1.0163, 1.0373, 1.0865, 1.2101,  &
1.0001, 1.0002, 1.0005, 1.0015, 1.0041, 1.0100, 1.0233, 1.0537, 1.1258,  &
1.0001, 1.0001, 1.0002, 1.0005, 1.0017, 1.0054, 1.0140, 1.0335, 1.0783,  &
1.0001, 1.0001, 1.0001, 1.0002, 1.0005, 1.0019, 1.0071, 1.0198, 1.0484,  &
1.0000, 1.0001, 1.0001, 1.0001, 1.0003, 1.0005, 1.0020, 1.0093, 1.0275 /
data ((supersat( 5,iw,iconc,6, 6),iw=1,9),iconc=1,7)/    &
1.0006, 1.0014, 1.0030, 1.0067, 1.0153, 1.0355, 1.0842, 1.2077, 1.5618,  &
1.0004, 1.0009, 1.0018, 1.0041, 1.0092, 1.0212, 1.0493, 1.1174, 1.2954,  &
1.0002, 1.0004, 1.0009, 1.0024, 1.0056, 1.0129, 1.0297, 1.0691, 1.1663,  &
1.0001, 1.0002, 1.0004, 1.0011, 1.0031, 1.0077, 1.0180, 1.0419, 1.0983,  &
1.0001, 1.0001, 1.0002, 1.0004, 1.0013, 1.0040, 1.0105, 1.0255, 1.0599,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0004, 1.0013, 1.0051, 1.0145, 1.0362,  &
1.0000, 1.0000, 1.0001, 1.0001, 1.0002, 1.0004, 1.0013, 1.0063, 1.0199 /
data ((supersat( 5,iw,iconc,7, 6),iw=1,9),iconc=1,7)/    &
1.0006, 1.0010, 1.0024, 1.0056, 1.0129, 1.0302, 1.0717, 1.1754, 1.4601,  &
1.0003, 1.0006, 1.0015, 1.0034, 1.0077, 1.0178, 1.0415, 1.0988, 1.2454,  &
1.0002, 1.0003, 1.0008, 1.0019, 1.0046, 1.0106, 1.0245, 1.0574, 1.1378,  &
1.0001, 1.0002, 1.0003, 1.0009, 1.0024, 1.0061, 1.0146, 1.0342, 1.0804,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0010, 1.0031, 1.0083, 1.0203, 1.0480,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0003, 1.0010, 1.0038, 1.0111, 1.0283,  &
1.0000, 1.0000, 1.0001, 1.0001, 1.0002, 1.0004, 1.0010, 1.0045, 1.0149 /
data ((supersat( 5,iw,iconc,1, 7),iw=1,9),iconc=1,7)/    &
1.0025, 1.0057, 1.0128, 1.0286, 1.0646, 1.1519, 1.4023, 2.6841, 4.9511,  &
1.0011, 1.0032, 1.0080, 1.0188, 1.0425, 1.0972, 1.2358, 1.6952, 3.9477,  &
1.0004, 1.0012, 1.0042, 1.0115, 1.0275, 1.0636, 1.1492, 1.3844, 2.4735,  &
1.0002, 1.0004, 1.0013, 1.0055, 1.0163, 1.0405, 1.0962, 1.2356, 1.6829,  &
1.0001, 1.0003, 1.0005, 1.0013, 1.0071, 1.0228, 1.0591, 1.1464, 1.3879,  &
1.0000, 1.0001, 1.0002, 1.0006, 1.0014, 1.0087, 1.0310, 1.0850, 1.2253,  &
0.9994, 0.9995, 0.9996, 0.9999, 1.0005, 1.0021, 1.0103, 1.0405, 1.1196 /
data ((supersat( 5,iw,iconc,2, 7),iw=1,9),iconc=1,7)/    &
1.0018, 1.0040, 1.0088, 1.0196, 1.0441, 1.1019, 1.2537, 1.7964, 3.6947,  &
1.0009, 1.0021, 1.0055, 1.0127, 1.0288, 1.0651, 1.1530, 1.4017, 2.5635,  &
1.0003, 1.0010, 1.0030, 1.0078, 1.0186, 1.0426, 1.0980, 1.2375, 1.6917,  &
1.0002, 1.0003, 1.0010, 1.0038, 1.0110, 1.0272, 1.0637, 1.1502, 1.3867,  &
1.0001, 1.0002, 1.0003, 1.0010, 1.0048, 1.0154, 1.0397, 1.0957, 1.2362,  &
1.0001, 1.0001, 1.0002, 1.0004, 1.0010, 1.0059, 1.0212, 1.0573, 1.1445,  &
0.9998, 0.9999, 0.9999, 1.0001, 1.0005, 1.0013, 1.0070, 1.0283, 1.0809 /
data ((supersat( 5,iw,iconc,3, 7),iw=1,9),iconc=1,7)/    &
1.0013, 1.0028, 1.0063, 1.0140, 1.0315, 1.0724, 1.1750, 1.4809, 2.7268,  &
1.0007, 1.0015, 1.0039, 1.0090, 1.0203, 1.0458, 1.1060, 1.2633, 1.8058,  &
1.0003, 1.0007, 1.0021, 1.0054, 1.0130, 1.0296, 1.0675, 1.1590, 1.4156,  &
1.0001, 1.0003, 1.0007, 1.0027, 1.0076, 1.0187, 1.0437, 1.1012, 1.2462,  &
1.0001, 1.0002, 1.0003, 1.0008, 1.0033, 1.0105, 1.0271, 1.0648, 1.1543,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0008, 1.0039, 1.0144, 1.0390, 1.0962,  &
1.0000, 1.0000, 1.0000, 1.0002, 1.0004, 1.0009, 1.0044, 1.0192, 1.0550 /
data ((supersat( 5,iw,iconc,4, 7),iw=1,9),iconc=1,7)/    &
1.0009, 1.0021, 1.0047, 1.0105, 1.0237, 1.0544, 1.1297, 1.3350, 2.0881,  &
1.0005, 1.0011, 1.0028, 1.0066, 1.0149, 1.0337, 1.0779, 1.1885, 1.5133,  &
1.0002, 1.0006, 1.0016, 1.0040, 1.0093, 1.0214, 1.0488, 1.1137, 1.2827,  &
1.0001, 1.0002, 1.0006, 1.0019, 1.0053, 1.0133, 1.0310, 1.0716, 1.1699,  &
1.0001, 1.0001, 1.0002, 1.0006, 1.0022, 1.0072, 1.0189, 1.0453, 1.1066,  &
1.0001, 1.0001, 1.0001, 1.0002, 1.0006, 1.0025, 1.0097, 1.0269, 1.0663,  &
1.0000, 1.0000, 1.0001, 1.0001, 1.0003, 1.0006, 1.0026, 1.0127, 1.0376 /
data ((supersat( 5,iw,iconc,5, 7),iw=1,9),iconc=1,7)/    &
1.0006, 1.0017, 1.0037, 1.0082, 1.0186, 1.0430, 1.1021, 1.2558, 1.7354,  &
1.0004, 1.0010, 1.0021, 1.0051, 1.0114, 1.0261, 1.0605, 1.1447, 1.3738,  &
1.0002, 1.0004, 1.0011, 1.0030, 1.0070, 1.0161, 1.0371, 1.0863, 1.2098,  &
1.0001, 1.0002, 1.0004, 1.0014, 1.0039, 1.0097, 1.0230, 1.0533, 1.1254,  &
1.0001, 1.0001, 1.0002, 1.0004, 1.0015, 1.0051, 1.0136, 1.0330, 1.0776,  &
1.0001, 1.0001, 1.0001, 1.0002, 1.0004, 1.0016, 1.0066, 1.0190, 1.0475,  &
1.0000, 1.0000, 1.0001, 1.0001, 1.0002, 1.0005, 1.0016, 1.0084, 1.0263 /
data ((supersat( 5,iw,iconc,6, 7),iw=1,9),iconc=1,7)/    &
1.0006, 1.0014, 1.0030, 1.0067, 1.0152, 1.0355, 1.0841, 1.2076, 1.5617,  &
1.0004, 1.0008, 1.0017, 1.0040, 1.0092, 1.0211, 1.0492, 1.1173, 1.2952,  &
1.0002, 1.0004, 1.0009, 1.0023, 1.0055, 1.0128, 1.0295, 1.0689, 1.1661,  &
1.0001, 1.0002, 1.0004, 1.0010, 1.0029, 1.0075, 1.0178, 1.0417, 1.0980,  &
1.0001, 1.0001, 1.0002, 1.0004, 1.0011, 1.0038, 1.0102, 1.0251, 1.0594,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0004, 1.0011, 1.0047, 1.0140, 1.0355,  &
1.0000, 1.0000, 1.0001, 1.0001, 1.0002, 1.0004, 1.0011, 1.0057, 1.0189 /
data ((supersat( 5,iw,iconc,7, 7),iw=1,9),iconc=1,7)/    &
1.0006, 1.0010, 1.0024, 1.0056, 1.0129, 1.0302, 1.0716, 1.1754, 1.4601,  &
1.0003, 1.0006, 1.0015, 1.0033, 1.0076, 1.0177, 1.0414, 1.0987, 1.2453,  &
1.0001, 1.0003, 1.0008, 1.0018, 1.0044, 1.0105, 1.0244, 1.0573, 1.1376,  &
1.0001, 1.0001, 1.0003, 1.0008, 1.0023, 1.0060, 1.0144, 1.0340, 1.0802,  &
1.0001, 1.0001, 1.0001, 1.0003, 1.0008, 1.0029, 1.0080, 1.0200, 1.0476,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0003, 1.0008, 1.0034, 1.0107, 1.0277,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0002, 1.0003, 1.0009, 1.0040, 1.0141 /
data ((supersat( 6,iw,iconc,1, 1),iw=1,9),iconc=1,7)/    &
1.0028, 1.0061, 1.0133, 1.0292, 1.0653, 1.1529, 1.4038, 2.6882, 4.9580,  &
1.0015, 1.0036, 1.0087, 1.0196, 1.0435, 1.0985, 1.2379, 1.6992, 3.9613,  &
1.0006, 1.0018, 1.0050, 1.0126, 1.0290, 1.0655, 1.1520, 1.3892, 2.4866,  &
1.0003, 1.0006, 1.0021, 1.0069, 1.0181, 1.0430, 1.0999, 1.2417, 1.6960,  &
1.0002, 1.0003, 1.0006, 1.0023, 1.0093, 1.0260, 1.0638, 1.1541, 1.4025,  &
1.0001, 1.0002, 1.0004, 1.0007, 1.0024, 1.0121, 1.0364, 1.0940, 1.2418,  &
0.9998, 0.9999, 1.0000, 1.0003, 1.0009, 1.0025, 1.0151, 1.0493, 1.1365 /
data ((supersat( 6,iw,iconc,2, 1),iw=1,9),iconc=1,7)/    &
1.0020, 1.0042, 1.0092, 1.0200, 1.0445, 1.1024, 1.2545, 1.7977, 3.6985,  &
1.0011, 1.0025, 1.0060, 1.0133, 1.0294, 1.0659, 1.1541, 1.4034, 2.5683,  &
1.0005, 1.0013, 1.0036, 1.0086, 1.0195, 1.0438, 1.0995, 1.2399, 1.6965,  &
1.0002, 1.0005, 1.0015, 1.0048, 1.0123, 1.0288, 1.0659, 1.1535, 1.3925,  &
1.0001, 1.0002, 1.0005, 1.0017, 1.0063, 1.0175, 1.0426, 1.1001, 1.2435,  &
1.0001, 1.0002, 1.0003, 1.0005, 1.0017, 1.0082, 1.0248, 1.0627, 1.1534,  &
1.0000, 1.0000, 1.0001, 1.0003, 1.0007, 1.0018, 1.0104, 1.0341, 1.0909 /
data ((supersat( 6,iw,iconc,3, 1),iw=1,9),iconc=1,7)/    &
1.0015, 1.0030, 1.0065, 1.0143, 1.0318, 1.0727, 1.1754, 1.4816, 2.7286,  &
1.0009, 1.0018, 1.0042, 1.0094, 1.0207, 1.0463, 1.1067, 1.2642, 1.8075,  &
1.0004, 1.0009, 1.0026, 1.0059, 1.0136, 1.0304, 1.0685, 1.1604, 1.4179,  &
1.0002, 1.0004, 1.0011, 1.0033, 1.0085, 1.0198, 1.0451, 1.1031, 1.2493,  &
1.0001, 1.0002, 1.0004, 1.0012, 1.0043, 1.0119, 1.0291, 1.0675, 1.1585,  &
1.0001, 1.0001, 1.0002, 1.0004, 1.0012, 1.0055, 1.0168, 1.0424, 1.1016,  &
1.0001, 1.0001, 1.0001, 1.0003, 1.0005, 1.0013, 1.0068, 1.0232, 1.0613 /
data ((supersat( 6,iw,iconc,4, 1),iw=1,9),iconc=1,7)/    &
1.0011, 1.0023, 1.0049, 1.0107, 1.0238, 1.0546, 1.1300, 1.3354, 2.0888,  &
1.0006, 1.0012, 1.0031, 1.0069, 1.0152, 1.0341, 1.0784, 1.1890, 1.5142,  &
1.0003, 1.0007, 1.0019, 1.0044, 1.0098, 1.0220, 1.0495, 1.1146, 1.2839,  &
1.0002, 1.0003, 1.0008, 1.0024, 1.0060, 1.0141, 1.0320, 1.0729, 1.1717,  &
1.0001, 1.0002, 1.0003, 1.0009, 1.0030, 1.0083, 1.0203, 1.0472, 1.1092,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0009, 1.0036, 1.0115, 1.0293, 1.0698,  &
1.0001, 1.0001, 1.0001, 1.0002, 1.0004, 1.0009, 1.0043, 1.0156, 1.0420 /
data ((supersat( 6,iw,iconc,5, 1),iw=1,9),iconc=1,7)/    &
1.0006, 1.0018, 1.0038, 1.0083, 1.0188, 1.0432, 1.1022, 1.2560, 1.7357,  &
1.0004, 1.0011, 1.0023, 1.0053, 1.0117, 1.0264, 1.0608, 1.1451, 1.3743,  &
1.0003, 1.0006, 1.0014, 1.0033, 1.0073, 1.0165, 1.0376, 1.0869, 1.2105,  &
1.0001, 1.0003, 1.0006, 1.0018, 1.0044, 1.0103, 1.0237, 1.0542, 1.1266,  &
1.0001, 1.0001, 1.0003, 1.0006, 1.0021, 1.0059, 1.0147, 1.0343, 1.0794,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0006, 1.0024, 1.0080, 1.0209, 1.0499,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0003, 1.0007, 1.0027, 1.0106, 1.0294 /
data ((supersat( 6,iw,iconc,6, 1),iw=1,9),iconc=1,7)/    &
1.0007, 1.0015, 1.0031, 1.0068, 1.0153, 1.0356, 1.0842, 1.2078, 1.5620,  &
1.0004, 1.0009, 1.0018, 1.0042, 1.0093, 1.0213, 1.0494, 1.1175, 1.2956,  &
1.0002, 1.0005, 1.0010, 1.0026, 1.0058, 1.0131, 1.0299, 1.0694, 1.1666,  &
1.0001, 1.0002, 1.0005, 1.0014, 1.0034, 1.0080, 1.0184, 1.0423, 1.0989,  &
1.0001, 1.0001, 1.0002, 1.0005, 1.0015, 1.0044, 1.0110, 1.0261, 1.0607,  &
1.0001, 1.0001, 1.0001, 1.0002, 1.0005, 1.0017, 1.0057, 1.0154, 1.0374,  &
1.0000, 1.0001, 1.0001, 1.0001, 1.0002, 1.0005, 1.0018, 1.0074, 1.0214 /
data ((supersat( 6,iw,iconc,7, 1),iw=1,9),iconc=1,7)/    &
1.0006, 1.0010, 1.0025, 1.0057, 1.0130, 1.0303, 1.0717, 1.1755, 1.4602,  &
1.0004, 1.0006, 1.0016, 1.0035, 1.0078, 1.0179, 1.0416, 1.0989, 1.2456,  &
1.0002, 1.0004, 1.0009, 1.0020, 1.0047, 1.0107, 1.0247, 1.0576, 1.1381,  &
1.0001, 1.0002, 1.0004, 1.0010, 1.0027, 1.0064, 1.0149, 1.0345, 1.0808,  &
1.0001, 1.0001, 1.0002, 1.0004, 1.0012, 1.0034, 1.0087, 1.0208, 1.0487,  &
1.0001, 1.0001, 1.0001, 1.0002, 1.0004, 1.0013, 1.0043, 1.0119, 1.0292,  &
1.0000, 1.0001, 1.0001, 1.0001, 1.0002, 1.0004, 1.0013, 1.0053, 1.0162 /
data ((supersat( 6,iw,iconc,1, 2),iw=1,9),iconc=1,7)/    &
1.0024, 1.0056, 1.0127, 1.0285, 1.0644, 1.1518, 1.4020, 2.6832, 4.9495,  &
1.0010, 1.0031, 1.0079, 1.0186, 1.0422, 1.0969, 1.2354, 1.6943, 3.9445,  &
1.0004, 1.0011, 1.0041, 1.0112, 1.0272, 1.0632, 1.1486, 1.3832, 2.4703,  &
1.0002, 1.0004, 1.0012, 1.0052, 1.0158, 1.0399, 1.0953, 1.2340, 1.6796,  &
1.0001, 1.0002, 1.0004, 1.0012, 1.0066, 1.0220, 1.0580, 1.1445, 1.3842,  &
0.9999, 1.0000, 1.0002, 1.0006, 1.0013, 1.0080, 1.0297, 1.0827, 1.2210,  &
0.9993, 0.9993, 0.9994, 0.9997, 1.0004, 1.0020, 1.0093, 1.0383, 1.1152 /
data ((supersat( 6,iw,iconc,2, 2),iw=1,9),iconc=1,7)/    &
1.0017, 1.0039, 1.0088, 1.0196, 1.0440, 1.1018, 1.2536, 1.7961, 3.6938,  &
1.0008, 1.0021, 1.0054, 1.0126, 1.0287, 1.0650, 1.1528, 1.4013, 2.5624,  &
1.0003, 1.0009, 1.0028, 1.0077, 1.0184, 1.0424, 1.0976, 1.2370, 1.6907,  &
1.0002, 1.0003, 1.0009, 1.0036, 1.0107, 1.0269, 1.0632, 1.1494, 1.3854,  &
1.0001, 1.0002, 1.0003, 1.0009, 1.0045, 1.0150, 1.0391, 1.0947, 1.2345,  &
1.0000, 1.0001, 1.0002, 1.0004, 1.0009, 1.0054, 1.0204, 1.0560, 1.1423,  &
0.9998, 0.9998, 0.9999, 1.0000, 1.0004, 1.0013, 1.0063, 1.0269, 1.0784 /
data ((supersat( 6,iw,iconc,3, 2),iw=1,9),iconc=1,7)/    &
1.0013, 1.0028, 1.0063, 1.0140, 1.0315, 1.0723, 1.1749, 1.4808, 2.7264,  &
1.0006, 1.0015, 1.0038, 1.0089, 1.0202, 1.0457, 1.1059, 1.2631, 1.8055,  &
1.0003, 1.0007, 1.0020, 1.0053, 1.0128, 1.0295, 1.0673, 1.1588, 1.4152,  &
1.0001, 1.0003, 1.0007, 1.0025, 1.0074, 1.0185, 1.0434, 1.1008, 1.2456,  &
1.0001, 1.0001, 1.0003, 1.0007, 1.0030, 1.0101, 1.0267, 1.0642, 1.1534,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0007, 1.0035, 1.0138, 1.0382, 1.0950,  &
0.9999, 1.0000, 1.0000, 1.0001, 1.0004, 1.0009, 1.0039, 1.0183, 1.0535 /
data ((supersat( 6,iw,iconc,4, 2),iw=1,9),iconc=1,7)/    &
1.0009, 1.0021, 1.0047, 1.0105, 1.0236, 1.0544, 1.1296, 1.3350, 2.0880,  &
1.0004, 1.0011, 1.0028, 1.0066, 1.0148, 1.0337, 1.0779, 1.1884, 1.5132,  &
1.0002, 1.0005, 1.0015, 1.0039, 1.0092, 1.0213, 1.0487, 1.1136, 1.2824,  &
1.0001, 1.0002, 1.0005, 1.0018, 1.0052, 1.0131, 1.0308, 1.0714, 1.1695,  &
1.0001, 1.0001, 1.0002, 1.0005, 1.0020, 1.0070, 1.0186, 1.0450, 1.1060,  &
1.0001, 1.0001, 1.0001, 1.0002, 1.0005, 1.0022, 1.0093, 1.0264, 1.0656,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0003, 1.0006, 1.0023, 1.0121, 1.0367 /
data ((supersat( 6,iw,iconc,5, 2),iw=1,9),iconc=1,7)/    &
1.0006, 1.0017, 1.0037, 1.0082, 1.0186, 1.0430, 1.1020, 1.2557, 1.7353,  &
1.0004, 1.0009, 1.0020, 1.0050, 1.0114, 1.0261, 1.0604, 1.1447, 1.3738,  &
1.0002, 1.0004, 1.0011, 1.0029, 1.0069, 1.0161, 1.0370, 1.0862, 1.2096,  &
1.0001, 1.0002, 1.0004, 1.0013, 1.0038, 1.0096, 1.0229, 1.0532, 1.1252,  &
1.0001, 1.0001, 1.0002, 1.0004, 1.0014, 1.0049, 1.0134, 1.0328, 1.0772,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0004, 1.0015, 1.0063, 1.0187, 1.0470,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0002, 1.0004, 1.0015, 1.0080, 1.0257 /
data ((supersat( 6,iw,iconc,6, 2),iw=1,9),iconc=1,7)/    &
1.0006, 1.0013, 1.0029, 1.0066, 1.0152, 1.0354, 1.0841, 1.2076, 1.5617,  &
1.0003, 1.0008, 1.0017, 1.0040, 1.0091, 1.0211, 1.0491, 1.1172, 1.2952,  &
1.0002, 1.0003, 1.0009, 1.0022, 1.0054, 1.0127, 1.0295, 1.0689, 1.1660,  &
1.0001, 1.0002, 1.0003, 1.0010, 1.0028, 1.0074, 1.0177, 1.0415, 1.0979,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0010, 1.0036, 1.0101, 1.0249, 1.0592,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0003, 1.0010, 1.0044, 1.0137, 1.0352,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0002, 1.0004, 1.0010, 1.0054, 1.0185 /
data ((supersat( 6,iw,iconc,7, 2),iw=1,9),iconc=1,7)/    &
1.0005, 1.0010, 1.0024, 1.0056, 1.0129, 1.0301, 1.0716, 1.1754, 1.4601,  &
1.0003, 1.0005, 1.0015, 1.0033, 1.0076, 1.0177, 1.0414, 1.0987, 1.2453,  &
1.0001, 1.0003, 1.0007, 1.0018, 1.0044, 1.0104, 1.0244, 1.0572, 1.1376,  &
1.0001, 1.0001, 1.0003, 1.0008, 1.0022, 1.0059, 1.0143, 1.0339, 1.0801,  &
1.0001, 1.0001, 1.0001, 1.0003, 1.0008, 1.0028, 1.0079, 1.0198, 1.0475,  &
1.0000, 1.0001, 1.0001, 1.0001, 1.0003, 1.0008, 1.0033, 1.0105, 1.0275,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0001, 1.0003, 1.0008, 1.0037, 1.0138 /
data ((supersat( 6,iw,iconc,1, 3),iw=1,9),iconc=1,7)/    &
1.0019, 1.0050, 1.0120, 1.0276, 1.0634, 1.1503, 1.3997, 2.6765, 4.9384,  &
1.0006, 1.0024, 1.0069, 1.0174, 1.0407, 1.0948, 1.2321, 1.6880, 3.9223,  &
1.0002, 1.0006, 1.0029, 1.0096, 1.0251, 1.0603, 1.1442, 1.3754, 2.4492,  &
1.0002, 1.0003, 1.0007, 1.0034, 1.0131, 1.0361, 1.0894, 1.2242, 1.6586,  &
1.0000, 1.0001, 1.0003, 1.0008, 1.0038, 1.0175, 1.0509, 1.1326, 1.3609,  &
0.9996, 0.9996, 0.9998, 1.0002, 1.0011, 1.0041, 1.0224, 1.0697, 1.1959,  &
0.9983, 0.9984, 0.9985, 0.9988, 0.9994, 1.0011, 1.0052, 1.0272, 1.0918 /
data ((supersat( 6,iw,iconc,2, 3),iw=1,9),iconc=1,7)/    &
1.0014, 1.0035, 1.0083, 1.0190, 1.0433, 1.1010, 1.2524, 1.7940, 3.6874,  &
1.0005, 1.0017, 1.0048, 1.0118, 1.0277, 1.0638, 1.1511, 1.3985, 2.5546,  &
1.0002, 1.0005, 1.0021, 1.0066, 1.0170, 1.0406, 1.0952, 1.2331, 1.6831,  &
1.0001, 1.0002, 1.0005, 1.0024, 1.0089, 1.0245, 1.0598, 1.1443, 1.3761,  &
1.0001, 1.0001, 1.0003, 1.0005, 1.0026, 1.0120, 1.0347, 1.0880, 1.2230,  &
0.9999, 0.9999, 1.0000, 1.0003, 1.0007, 1.0028, 1.0156, 1.0481, 1.1287,  &
0.9993, 0.9993, 0.9994, 0.9995, 0.9999, 1.0009, 1.0030, 1.0193, 1.0643 /
data ((supersat( 6,iw,iconc,3, 3),iw=1,9),iconc=1,7)/    &
1.0010, 1.0025, 1.0059, 1.0136, 1.0310, 1.0719, 1.1743, 1.4798, 2.7235,  &
1.0004, 1.0012, 1.0033, 1.0084, 1.0196, 1.0449, 1.1049, 1.2616, 1.8028,  &
1.0002, 1.0004, 1.0015, 1.0046, 1.0119, 1.0283, 1.0658, 1.1566, 1.4116,  &
1.0001, 1.0002, 1.0004, 1.0017, 1.0061, 1.0169, 1.0412, 1.0977, 1.2407,  &
1.0001, 1.0001, 1.0002, 1.0004, 1.0018, 1.0081, 1.0238, 1.0601, 1.1470,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0005, 1.0018, 1.0105, 1.0331, 1.0869,  &
0.9997, 0.9997, 0.9997, 0.9998, 1.0001, 1.0006, 1.0019, 1.0131, 1.0445 /
data ((supersat( 6,iw,iconc,4, 3),iw=1,9),iconc=1,7)/    &
1.0007, 1.0019, 1.0045, 1.0102, 1.0233, 1.0540, 1.1293, 1.3344, 2.0868,  &
1.0003, 1.0009, 1.0024, 1.0061, 1.0143, 1.0332, 1.0772, 1.1875, 1.5118,  &
1.0001, 1.0003, 1.0011, 1.0032, 1.0086, 1.0205, 1.0477, 1.1123, 1.2804,  &
1.0001, 1.0001, 1.0003, 1.0012, 1.0043, 1.0119, 1.0294, 1.0694, 1.1667,  &
1.0001, 1.0001, 1.0002, 1.0003, 1.0012, 1.0055, 1.0166, 1.0422, 1.1021,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0004, 1.0012, 1.0069, 1.0229, 1.0603,  &
0.9998, 0.9998, 0.9999, 0.9999, 1.0001, 1.0005, 1.0013, 1.0084, 1.0306 /
data ((supersat( 6,iw,iconc,5, 3),iw=1,9),iconc=1,7)/    &
1.0005, 1.0015, 1.0035, 1.0080, 1.0184, 1.0427, 1.1018, 1.2554, 1.7347,  &
1.0003, 1.0007, 1.0018, 1.0047, 1.0111, 1.0257, 1.0600, 1.1441, 1.3730,  &
1.0001, 1.0003, 1.0008, 1.0024, 1.0064, 1.0155, 1.0363, 1.0853, 1.2084,  &
1.0001, 1.0001, 1.0003, 1.0008, 1.0031, 1.0087, 1.0218, 1.0518, 1.1233,  &
1.0000, 1.0001, 1.0001, 1.0003, 1.0008, 1.0038, 1.0119, 1.0308, 1.0746,  &
1.0000, 1.0000, 1.0001, 1.0001, 1.0003, 1.0008, 1.0045, 1.0161, 1.0433,  &
0.9999, 0.9999, 0.9999, 1.0000, 1.0001, 1.0003, 1.0009, 1.0052, 1.0213 /
data ((supersat( 6,iw,iconc,6, 3),iw=1,9),iconc=1,7)/    &
1.0005, 1.0012, 1.0028, 1.0065, 1.0151, 1.0353, 1.0839, 1.2074, 1.5613,  &
1.0002, 1.0005, 1.0015, 1.0037, 1.0089, 1.0208, 1.0488, 1.1169, 1.2947,  &
1.0001, 1.0002, 1.0006, 1.0019, 1.0050, 1.0122, 1.0289, 1.0682, 1.1652,  &
1.0001, 1.0001, 1.0002, 1.0006, 1.0023, 1.0067, 1.0169, 1.0405, 1.0966,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0006, 1.0027, 1.0089, 1.0235, 1.0573,  &
1.0000, 1.0000, 1.0001, 1.0001, 1.0002, 1.0006, 1.0031, 1.0117, 1.0324,  &
0.9999, 0.9999, 0.9999, 1.0000, 1.0001, 1.0003, 1.0007, 1.0033, 1.0151 /
data ((supersat( 6,iw,iconc,7, 3),iw=1,9),iconc=1,7)/    &
1.0004, 1.0009, 1.0023, 1.0055, 1.0127, 1.0300, 1.0714, 1.1752, 1.4598,  &
1.0002, 1.0005, 1.0013, 1.0031, 1.0073, 1.0174, 1.0411, 1.0984, 1.2449,  &
1.0001, 1.0002, 1.0005, 1.0015, 1.0040, 1.0100, 1.0239, 1.0568, 1.1370,  &
1.0001, 1.0001, 1.0002, 1.0005, 1.0018, 1.0053, 1.0136, 1.0331, 1.0791,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0005, 1.0020, 1.0069, 1.0186, 1.0460,  &
1.0000, 1.0000, 1.0001, 1.0001, 1.0002, 1.0005, 1.0022, 1.0088, 1.0253,  &
0.9999, 0.9999, 0.9999, 1.0000, 1.0001, 1.0002, 1.0005, 1.0022, 1.0111 /
data ((supersat( 6,iw,iconc,1, 4),iw=1,9),iconc=1,7)/    &
1.0013, 1.0043, 1.0112, 1.0266, 1.0621, 1.1485, 1.3967, 2.6677, 4.9235,  &
1.0004, 1.0016, 1.0058, 1.0159, 1.0388, 1.0922, 1.2280, 1.6798, 3.8932,  &
1.0002, 1.0004, 1.0017, 1.0077, 1.0226, 1.0567, 1.1386, 1.3654, 2.4219,  &
1.0001, 1.0002, 1.0005, 1.0018, 1.0101, 1.0315, 1.0823, 1.2119, 1.6318,  &
0.9998, 0.9999, 1.0001, 1.0006, 1.0019, 1.0127, 1.0427, 1.1182, 1.3323,  &
0.9989, 0.9990, 0.9991, 0.9995, 1.0005, 1.0027, 1.0152, 1.0555, 1.1671,  &
0.9969, 0.9969, 0.9970, 0.9973, 0.9980, 0.9996, 1.0037, 1.0172, 1.0686 /
data ((supersat( 6,iw,iconc,2, 4),iw=1,9),iconc=1,7)/    &
1.0010, 1.0031, 1.0078, 1.0184, 1.0425, 1.1000, 1.2510, 1.7913, 3.6792,  &
1.0003, 1.0012, 1.0041, 1.0109, 1.0265, 1.0622, 1.1490, 1.3949, 2.5444,  &
1.0001, 1.0003, 1.0013, 1.0053, 1.0154, 1.0385, 1.0921, 1.2283, 1.6733,  &
1.0001, 1.0002, 1.0003, 1.0013, 1.0069, 1.0216, 1.0557, 1.1377, 1.3643,  &
1.0000, 1.0000, 1.0002, 1.0004, 1.0013, 1.0087, 1.0296, 1.0799, 1.2087,  &
0.9995, 0.9996, 0.9997, 0.9999, 1.0005, 1.0017, 1.0106, 1.0392, 1.1127,  &
0.9984, 0.9984, 0.9985, 0.9986, 0.9991, 1.0001, 1.0025, 1.0122, 1.0496 /
data ((supersat( 6,iw,iconc,3, 4),iw=1,9),iconc=1,7)/    &
1.0007, 1.0022, 1.0055, 1.0132, 1.0305, 1.0712, 1.1735, 1.4785, 2.7197,  &
1.0003, 1.0009, 1.0029, 1.0077, 1.0188, 1.0439, 1.1037, 1.2598, 1.7993,  &
1.0001, 1.0003, 1.0009, 1.0037, 1.0108, 1.0269, 1.0640, 1.1539, 1.4070,  &
1.0001, 1.0001, 1.0003, 1.0009, 1.0047, 1.0150, 1.0386, 1.0939, 1.2345,  &
1.0000, 1.0001, 1.0001, 1.0003, 1.0010, 1.0058, 1.0205, 1.0550, 1.1389,  &
0.9998, 0.9998, 0.9999, 1.0001, 1.0004, 1.0011, 1.0070, 1.0273, 1.0772,  &
0.9992, 0.9992, 0.9992, 0.9993, 0.9996, 1.0002, 1.0016, 1.0079, 1.0348 /
data ((supersat( 6,iw,iconc,4, 4),iw=1,9),iconc=1,7)/    &
1.0006, 1.0017, 1.0042, 1.0099, 1.0230, 1.0536, 1.1288, 1.3337, 2.0853,  &
1.0002, 1.0006, 1.0021, 1.0056, 1.0138, 1.0325, 1.0764, 1.1865, 1.5100,  &
1.0001, 1.0002, 1.0007, 1.0027, 1.0077, 1.0195, 1.0465, 1.1106, 1.2779,  &
1.0001, 1.0001, 1.0002, 1.0007, 1.0032, 1.0105, 1.0276, 1.0670, 1.1631,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0007, 1.0038, 1.0142, 1.0389, 1.0971,  &
0.9999, 0.9999, 1.0000, 1.0001, 1.0003, 1.0007, 1.0044, 1.0188, 1.0540,  &
0.9995, 0.9995, 0.9995, 0.9996, 0.9998, 1.0002, 1.0010, 1.0047, 1.0240 /
data ((supersat( 6,iw,iconc,5, 4),iw=1,9),iconc=1,7)/    &
1.0005, 1.0013, 1.0032, 1.0077, 1.0181, 1.0425, 1.1014, 1.2550, 1.7340,  &
1.0002, 1.0005, 1.0016, 1.0043, 1.0106, 1.0252, 1.0594, 1.1434, 1.3719,  &
1.0001, 1.0002, 1.0005, 1.0020, 1.0058, 1.0147, 1.0354, 1.0842, 1.2068,  &
1.0001, 1.0001, 1.0002, 1.0005, 1.0023, 1.0077, 1.0205, 1.0501, 1.1210,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0005, 1.0026, 1.0101, 1.0284, 1.0712,  &
0.9999, 1.0000, 1.0000, 1.0001, 1.0002, 1.0006, 1.0027, 1.0131, 1.0389,  &
0.9996, 0.9997, 0.9997, 0.9997, 0.9998, 1.0001, 1.0007, 1.0028, 1.0164 /
data ((supersat( 6,iw,iconc,6, 4),iw=1,9),iconc=1,7)/    &
1.0004, 1.0010, 1.0026, 1.0063, 1.0148, 1.0350, 1.0837, 1.2071, 1.5609,  &
1.0002, 1.0004, 1.0013, 1.0034, 1.0085, 1.0204, 1.0484, 1.1164, 1.2940,  &
1.0001, 1.0002, 1.0004, 1.0015, 1.0045, 1.0116, 1.0282, 1.0674, 1.1641,  &
1.0000, 1.0001, 1.0002, 1.0004, 1.0017, 1.0058, 1.0159, 1.0393, 1.0949,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0004, 1.0018, 1.0075, 1.0216, 1.0548,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0002, 1.0004, 1.0018, 1.0094, 1.0291,  &
0.9997, 0.9997, 0.9997, 0.9998, 0.9999, 1.0001, 1.0005, 1.0019, 1.0114 /
data ((supersat( 6,iw,iconc,7, 4),iw=1,9),iconc=1,7)/    &
1.0003, 1.0008, 1.0021, 1.0053, 1.0126, 1.0298, 1.0712, 1.1750, 1.4595,  &
1.0001, 1.0003, 1.0010, 1.0028, 1.0071, 1.0171, 1.0408, 1.0980, 1.2444,  &
1.0001, 1.0001, 1.0003, 1.0012, 1.0036, 1.0095, 1.0234, 1.0561, 1.1362,  &
1.0000, 1.0001, 1.0001, 1.0003, 1.0013, 1.0046, 1.0128, 1.0321, 1.0778,  &
1.0000, 1.0000, 1.0001, 1.0001, 1.0003, 1.0013, 1.0057, 1.0171, 1.0441,  &
1.0000, 1.0000, 1.0000, 1.0000, 1.0001, 1.0003, 1.0013, 1.0069, 1.0227,  &
0.9997, 0.9997, 0.9998, 0.9998, 0.9999, 1.0000, 1.0004, 1.0014, 1.0082 /
data ((supersat( 6,iw,iconc,1, 5),iw=1,9),iconc=1,7)/    &
1.0011, 1.0038, 1.0106, 1.0259, 1.0612, 1.1472, 1.3945, 2.6614, 4.9128,  &
1.0003, 1.0012, 1.0051, 1.0150, 1.0375, 1.0903, 1.2251, 1.6740, 3.8722,  &
1.0002, 1.0003, 1.0012, 1.0065, 1.0209, 1.0542, 1.1347, 1.3583, 2.4023,  &
1.0000, 1.0001, 1.0004, 1.0012, 1.0082, 1.0285, 1.0774, 1.2034, 1.6131,  &
0.9995, 0.9996, 0.9999, 1.0004, 1.0017, 1.0098, 1.0375, 1.1087, 1.3130,  &
0.9983, 0.9984, 0.9985, 0.9989, 0.9999, 1.0023, 1.0112, 1.0471, 1.1491,  &
0.9959, 0.9959, 0.9960, 0.9963, 0.9969, 0.9985, 1.0025, 1.0123, 1.0561 /
data ((supersat( 6,iw,iconc,2, 5),iw=1,9),iconc=1,7)/    &
1.0008, 1.0028, 1.0074, 1.0179, 1.0420, 1.0992, 1.2499, 1.7894, 3.6731,  &
1.0002, 1.0009, 1.0036, 1.0103, 1.0257, 1.0612, 1.1474, 1.3923, 2.5370,  &
1.0001, 1.0003, 1.0009, 1.0045, 1.0143, 1.0370, 1.0899, 1.2248, 1.6663,  &
1.0001, 1.0001, 1.0003, 1.0009, 1.0056, 1.0197, 1.0528, 1.1332, 1.3560,  &
0.9998, 0.9999, 1.0001, 1.0004, 1.0010, 1.0067, 1.0264, 1.0744, 1.1990,  &
0.9992, 0.9993, 0.9994, 0.9996, 1.0002, 1.0015, 1.0077, 1.0338, 1.1023,  &
0.9977, 0.9977, 0.9978, 0.9979, 0.9983, 0.9994, 1.0018, 1.0085, 1.0412 /
data ((supersat( 6,iw,iconc,3, 5),iw=1,9),iconc=1,7)/    &
1.0006, 1.0020, 1.0052, 1.0129, 1.0302, 1.0708, 1.1729, 1.4776, 2.7170,  &
1.0002, 1.0007, 1.0026, 1.0073, 1.0182, 1.0432, 1.1028, 1.2584, 1.7968,  &
1.0001, 1.0002, 1.0007, 1.0032, 1.0100, 1.0259, 1.0627, 1.1520, 1.4036,  &
1.0001, 1.0001, 1.0002, 1.0007, 1.0038, 1.0136, 1.0368, 1.0912, 1.2302,  &
1.0000, 1.0000, 1.0001, 1.0003, 1.0007, 1.0044, 1.0183, 1.0516, 1.1334,  &
0.9996, 0.9997, 0.9997, 0.9999, 1.0002, 1.0010, 1.0049, 1.0236, 1.0709,  &
0.9987, 0.9987, 0.9987, 0.9989, 0.9991, 0.9998, 1.0012, 1.0053, 1.0291 /
data ((supersat( 6,iw,iconc,4, 5),iw=1,9),iconc=1,7)/    &
1.0005, 1.0015, 1.0040, 1.0096, 1.0227, 1.0533, 1.1284, 1.3332, 2.0842,  &
1.0002, 1.0005, 1.0019, 1.0053, 1.0134, 1.0320, 1.0759, 1.1857, 1.5087,  &
1.0001, 1.0002, 1.0005, 1.0022, 1.0072, 1.0188, 1.0456, 1.1094, 1.2761,  &
1.0001, 1.0001, 1.0002, 1.0005, 1.0026, 1.0096, 1.0263, 1.0653, 1.1605,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0005, 1.0029, 1.0127, 1.0366, 1.0936,  &
0.9998, 0.9998, 0.9999, 1.0000, 1.0002, 1.0007, 1.0030, 1.0163, 1.0499,  &
0.9992, 0.9992, 0.9992, 0.9993, 0.9995, 0.9999, 1.0008, 1.0032, 1.0200 /
data ((supersat( 6,iw,iconc,5, 5),iw=1,9),iconc=1,7)/    &
1.0004, 1.0012, 1.0031, 1.0076, 1.0179, 1.0422, 1.1012, 1.2546, 1.7335,  &
1.0001, 1.0004, 1.0014, 1.0040, 1.0103, 1.0249, 1.0590, 1.1429, 1.3711,  &
1.0001, 1.0001, 1.0004, 1.0016, 1.0054, 1.0142, 1.0348, 1.0834, 1.2057,  &
1.0000, 1.0001, 1.0001, 1.0004, 1.0018, 1.0070, 1.0196, 1.0489, 1.1193,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0004, 1.0019, 1.0089, 1.0268, 1.0688,  &
0.9999, 0.9999, 0.9999, 1.0000, 1.0001, 1.0005, 1.0019, 1.0112, 1.0360,  &
0.9994, 0.9994, 0.9994, 0.9995, 0.9996, 0.9999, 1.0005, 1.0020, 1.0135 /
data ((supersat( 6,iw,iconc,6, 5),iw=1,9),iconc=1,7)/    &
1.0003, 1.0009, 1.0024, 1.0061, 1.0147, 1.0349, 1.0835, 1.2068, 1.5605,  &
1.0001, 1.0003, 1.0011, 1.0032, 1.0083, 1.0201, 1.0481, 1.1160, 1.2935,  &
1.0001, 1.0001, 1.0003, 1.0012, 1.0042, 1.0112, 1.0278, 1.0668, 1.1633,  &
1.0000, 1.0001, 1.0001, 1.0003, 1.0013, 1.0053, 1.0151, 1.0384, 1.0937,  &
1.0000, 1.0000, 1.0001, 1.0001, 1.0003, 1.0013, 1.0065, 1.0203, 1.0531,  &
0.9999, 0.9999, 0.9999, 1.0000, 1.0001, 1.0004, 1.0013, 1.0079, 1.0269,  &
0.9995, 0.9995, 0.9995, 0.9996, 0.9997, 0.9999, 1.0003, 1.0014, 1.0092 /
data ((supersat( 6,iw,iconc,7, 5),iw=1,9),iconc=1,7)/    &
1.0003, 1.0008, 1.0020, 1.0052, 1.0124, 1.0297, 1.0711, 1.1748, 1.4592,  &
1.0001, 1.0003, 1.0009, 1.0026, 1.0069, 1.0169, 1.0405, 1.0977, 1.2440,  &
1.0001, 1.0001, 1.0003, 1.0010, 1.0033, 1.0092, 1.0230, 1.0556, 1.1356,  &
1.0000, 1.0001, 1.0001, 1.0003, 1.0010, 1.0041, 1.0122, 1.0313, 1.0769,  &
1.0000, 1.0000, 1.0001, 1.0001, 1.0003, 1.0010, 1.0049, 1.0161, 1.0427,  &
0.9999, 0.9999, 0.9999, 1.0000, 1.0001, 1.0003, 1.0010, 1.0057, 1.0209,  &
0.9995, 0.9996, 0.9996, 0.9996, 0.9997, 0.9999, 1.0002, 1.0011, 1.0064 /
data ((supersat( 6,iw,iconc,1, 6),iw=1,9),iconc=1,7)/    &
1.0009, 1.0035, 1.0102, 1.0254, 1.0605, 1.1462, 1.3928, 2.6561, 4.9039,  &
1.0002, 1.0009, 1.0046, 1.0142, 1.0365, 1.0888, 1.2227, 1.6692, 3.8552,  &
1.0001, 1.0003, 1.0009, 1.0057, 1.0195, 1.0522, 1.1316, 1.3527, 2.3866,  &
0.9999, 1.0001, 1.0004, 1.0010, 1.0068, 1.0262, 1.0737, 1.1967, 1.5983,  &
0.9993, 0.9994, 0.9997, 1.0002, 1.0015, 1.0079, 1.0338, 1.1017, 1.2982,  &
0.9978, 0.9978, 0.9980, 0.9984, 0.9994, 1.0018, 1.0087, 1.0414, 1.1361,  &
0.9951, 0.9952, 0.9953, 0.9955, 0.9961, 0.9977, 1.0015, 1.0124, 1.0480 /
data ((supersat( 6,iw,iconc,2, 6),iw=1,9),iconc=1,7)/    &
1.0007, 1.0025, 1.0071, 1.0175, 1.0416, 1.0987, 1.2491, 1.7878, 3.6682,  &
1.0002, 1.0007, 1.0032, 1.0098, 1.0251, 1.0603, 1.1462, 1.3902, 2.5309,  &
1.0001, 1.0002, 1.0007, 1.0039, 1.0135, 1.0358, 1.0882, 1.2220, 1.6606,  &
1.0000, 1.0001, 1.0003, 1.0007, 1.0047, 1.0182, 1.0506, 1.1296, 1.3494,  &
0.9997, 0.9998, 0.9999, 1.0003, 1.0010, 1.0054, 1.0239, 1.0702, 1.1913,  &
0.9989, 0.9990, 0.9991, 0.9993, 0.9999, 1.0013, 1.0059, 1.0300, 1.0947,  &
0.9971, 0.9971, 0.9971, 0.9973, 0.9977, 0.9987, 1.0012, 1.0077, 1.0356 /
data ((supersat( 6,iw,iconc,3, 6),iw=1,9),iconc=1,7)/    &
1.0005, 1.0019, 1.0050, 1.0126, 1.0299, 1.0705, 1.1724, 1.4769, 2.7147,  &
1.0002, 1.0005, 1.0023, 1.0070, 1.0178, 1.0428, 1.1021, 1.2573, 1.7946,  &
1.0001, 1.0002, 1.0005, 1.0028, 1.0094, 1.0252, 1.0616, 1.1504, 1.4009,  &
1.0000, 1.0001, 1.0002, 1.0006, 1.0032, 1.0126, 1.0354, 1.0891, 1.2267,  &
0.9999, 0.9999, 1.0000, 1.0002, 1.0006, 1.0035, 1.0166, 1.0490, 1.1290,  &
0.9995, 0.9995, 0.9996, 0.9997, 1.0001, 1.0009, 1.0037, 1.0210, 1.0662,  &
0.9983, 0.9983, 0.9983, 0.9984, 0.9987, 0.9993, 1.0009, 1.0046, 1.0252 /
data ((supersat( 6,iw,iconc,4, 6),iw=1,9),iconc=1,7)/    &
1.0004, 1.0014, 1.0038, 1.0095, 1.0225, 1.0531, 1.1281, 1.3328, 2.0832,  &
1.0001, 1.0004, 1.0017, 1.0051, 1.0131, 1.0317, 1.0754, 1.1850, 1.5076,  &
1.0001, 1.0001, 1.0004, 1.0019, 1.0067, 1.0183, 1.0449, 1.1085, 1.2746,  &
1.0000, 1.0001, 1.0002, 1.0004, 1.0021, 1.0089, 1.0254, 1.0639, 1.1584,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0005, 1.0022, 1.0115, 1.0348, 1.0908,  &
0.9997, 0.9997, 0.9998, 0.9999, 1.0001, 1.0006, 1.0023, 1.0144, 1.0467,  &
0.9989, 0.9989, 0.9989, 0.9990, 0.9992, 0.9996, 1.0005, 1.0028, 1.0172 /
data ((supersat( 6,iw,iconc,5, 6),iw=1,9),iconc=1,7)/    &
1.0003, 1.0010, 1.0029, 1.0074, 1.0178, 1.0421, 1.1010, 1.2544, 1.7330,  &
1.0001, 1.0003, 1.0013, 1.0038, 1.0101, 1.0246, 1.0587, 1.1425, 1.3705,  &
1.0001, 1.0001, 1.0003, 1.0014, 1.0050, 1.0138, 1.0343, 1.0828, 1.2047,  &
1.0000, 1.0001, 1.0001, 1.0003, 1.0015, 1.0064, 1.0188, 1.0480, 1.1179,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0004, 1.0015, 1.0081, 1.0255, 1.0669,  &
0.9998, 0.9998, 0.9998, 0.9999, 1.0001, 1.0004, 1.0016, 1.0098, 1.0337,  &
0.9992, 0.9992, 0.9992, 0.9993, 0.9994, 0.9997, 1.0003, 1.0018, 1.0115 /
data ((supersat( 6,iw,iconc,6, 6),iw=1,9),iconc=1,7)/    &
1.0003, 1.0008, 1.0023, 1.0060, 1.0146, 1.0347, 1.0833, 1.2067, 1.5602,  &
1.0001, 1.0003, 1.0010, 1.0030, 1.0081, 1.0199, 1.0478, 1.1157, 1.2930,  &
1.0001, 1.0001, 1.0003, 1.0011, 1.0039, 1.0109, 1.0274, 1.0664, 1.1626,  &
1.0000, 1.0001, 1.0001, 1.0003, 1.0011, 1.0048, 1.0146, 1.0376, 1.0927,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0003, 1.0011, 1.0058, 1.0194, 1.0517,  &
0.9998, 0.9998, 0.9999, 0.9999, 1.0001, 1.0003, 1.0011, 1.0068, 1.0252,  &
0.9993, 0.9993, 0.9993, 0.9994, 0.9995, 0.9997, 1.0002, 1.0013, 1.0077 /
data ((supersat( 6,iw,iconc,7, 6),iw=1,9),iconc=1,7)/    &
1.0002, 1.0007, 1.0020, 1.0051, 1.0123, 1.0296, 1.0710, 1.1746, 1.4590,  &
1.0001, 1.0002, 1.0008, 1.0025, 1.0067, 1.0167, 1.0403, 1.0975, 1.2437,  &
1.0000, 1.0001, 1.0002, 1.0008, 1.0031, 1.0089, 1.0226, 1.0553, 1.1351,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0008, 1.0038, 1.0117, 1.0307, 1.0761,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0002, 1.0008, 1.0044, 1.0153, 1.0416,  &
0.9998, 0.9998, 0.9999, 0.9999, 1.0000, 1.0003, 1.0009, 1.0049, 1.0196,  &
0.9994, 0.9994, 0.9994, 0.9994, 0.9995, 0.9997, 1.0001, 1.0009, 1.0052 /
data ((supersat( 6,iw,iconc,1, 7),iw=1,9),iconc=1,7)/    &
1.0008, 1.0033, 1.0098, 1.0249, 1.0598, 1.1453, 1.3913, 2.6517, 4.8965,  &
1.0002, 1.0008, 1.0041, 1.0135, 1.0356, 1.0876, 1.2207, 1.6652, 3.8407,  &
1.0001, 1.0003, 1.0008, 1.0049, 1.0185, 1.0506, 1.1290, 1.3479, 2.3731,  &
0.9998, 1.0000, 1.0003, 1.0010, 1.0058, 1.0244, 1.0706, 1.1911, 1.5859,  &
0.9991, 0.9992, 0.9994, 1.0000, 1.0014, 1.0064, 1.0309, 1.0960, 1.2860,  &
0.9973, 0.9974, 0.9975, 0.9980, 0.9990, 1.0014, 1.0069, 1.0370, 1.1260,  &
0.9946, 0.9946, 0.9947, 0.9949, 0.9955, 0.9970, 1.0007, 1.0111, 1.0422 /
data ((supersat( 6,iw,iconc,2, 7),iw=1,9),iconc=1,7)/    &
1.0006, 1.0023, 1.0069, 1.0173, 1.0412, 1.0982, 1.2484, 1.7865, 3.6641,  &
1.0002, 1.0006, 1.0029, 1.0094, 1.0246, 1.0596, 1.1451, 1.3884, 2.5258,  &
1.0001, 1.0002, 1.0006, 1.0035, 1.0128, 1.0348, 1.0868, 1.2196, 1.6558,  &
1.0000, 1.0001, 1.0003, 1.0006, 1.0040, 1.0171, 1.0488, 1.1266, 1.3438,  &
0.9996, 0.9997, 0.9998, 1.0002, 1.0009, 1.0043, 1.0220, 1.0669, 1.1851,  &
0.9987, 0.9987, 0.9988, 0.9990, 0.9997, 1.0011, 1.0046, 1.0270, 1.0887,  &
0.9965, 0.9966, 0.9966, 0.9968, 0.9972, 0.9982, 1.0007, 1.0071, 1.0314 /
data ((supersat( 6,iw,iconc,3, 7),iw=1,9),iconc=1,7)/    &
1.0005, 1.0017, 1.0049, 1.0124, 1.0297, 1.0702, 1.1720, 1.4762, 2.7127,  &
1.0002, 1.0005, 1.0021, 1.0067, 1.0175, 1.0423, 1.1014, 1.2564, 1.7928,  &
1.0001, 1.0002, 1.0005, 1.0024, 1.0089, 1.0245, 1.0607, 1.1491, 1.3986,  &
1.0000, 1.0001, 1.0002, 1.0005, 1.0027, 1.0118, 1.0342, 1.0874, 1.2237,  &
0.9998, 0.9999, 1.0000, 1.0002, 1.0006, 1.0028, 1.0153, 1.0469, 1.1254,  &
0.9993, 0.9993, 0.9994, 0.9996, 0.9999, 1.0008, 1.0030, 1.0190, 1.0623,  &
0.9979, 0.9979, 0.9979, 0.9980, 0.9983, 0.9990, 1.0005, 1.0043, 1.0223 /
data ((supersat( 6,iw,iconc,4, 7),iw=1,9),iconc=1,7)/    &
1.0004, 1.0013, 1.0037, 1.0093, 1.0223, 1.0529, 1.1279, 1.3324, 2.0825,  &
1.0001, 1.0004, 1.0015, 1.0049, 1.0129, 1.0314, 1.0750, 1.1845, 1.5067,  &
1.0001, 1.0001, 1.0004, 1.0017, 1.0064, 1.0178, 1.0443, 1.1076, 1.2733,  &
1.0000, 1.0001, 1.0001, 1.0004, 1.0018, 1.0083, 1.0246, 1.0628, 1.1566,  &
0.9999, 1.0000, 1.0000, 1.0002, 1.0004, 1.0018, 1.0105, 1.0334, 1.0886,  &
0.9996, 0.9996, 0.9996, 0.9998, 1.0000, 1.0005, 1.0019, 1.0129, 1.0442,  &
0.9986, 0.9986, 0.9987, 0.9987, 0.9989, 0.9993, 1.0003, 1.0026, 1.0151 /
data ((supersat( 6,iw,iconc,5, 7),iw=1,9),iconc=1,7)/    &
1.0003, 1.0009, 1.0028, 1.0073, 1.0176, 1.0419, 1.1008, 1.2541, 1.7326,  &
1.0001, 1.0003, 1.0011, 1.0037, 1.0099, 1.0244, 1.0584, 1.1421, 1.3700,  &
1.0001, 1.0001, 1.0003, 1.0012, 1.0047, 1.0135, 1.0339, 1.0822, 1.2039,  &
1.0000, 1.0001, 1.0001, 1.0003, 1.0013, 1.0060, 1.0182, 1.0471, 1.1167,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0003, 1.0013, 1.0073, 1.0244, 1.0654,  &
0.9997, 0.9997, 0.9998, 0.9998, 1.0000, 1.0004, 1.0013, 1.0087, 1.0319,  &
0.9990, 0.9990, 0.9990, 0.9990, 0.9992, 0.9995, 1.0001, 1.0017, 1.0100 /
data ((supersat( 6,iw,iconc,6, 7),iw=1,9),iconc=1,7)/    &
1.0002, 1.0008, 1.0023, 1.0060, 1.0145, 1.0346, 1.0832, 1.2065, 1.5600,  &
1.0001, 1.0002, 1.0009, 1.0029, 1.0080, 1.0198, 1.0476, 1.1154, 1.2927,  &
1.0000, 1.0001, 1.0002, 1.0009, 1.0037, 1.0106, 1.0270, 1.0659, 1.1621,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0009, 1.0045, 1.0141, 1.0370, 1.0919,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0003, 1.0009, 1.0053, 1.0185, 1.0505,  &
0.9998, 0.9998, 0.9998, 0.9999, 1.0000, 1.0003, 1.0010, 1.0060, 1.0239,  &
0.9991, 0.9991, 0.9991, 0.9992, 0.9993, 0.9995, 1.0000, 1.0011, 1.0065 /
data ((supersat( 6,iw,iconc,7, 7),iw=1,9),iconc=1,7)/    &
1.0002, 1.0007, 1.0019, 1.0050, 1.0123, 1.0295, 1.0709, 1.1745, 1.4588,  &
1.0001, 1.0002, 1.0007, 1.0024, 1.0066, 1.0166, 1.0402, 1.0973, 1.2434,  &
1.0000, 1.0001, 1.0002, 1.0007, 1.0029, 1.0087, 1.0224, 1.0549, 1.1347,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0007, 1.0035, 1.0113, 1.0302, 1.0755,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0002, 1.0007, 1.0039, 1.0146, 1.0407,  &
0.9998, 0.9998, 0.9998, 0.9999, 1.0000, 1.0002, 1.0008, 1.0042, 1.0185,  &
0.9992, 0.9992, 0.9992, 0.9992, 0.9993, 0.9995, 1.0000, 1.0008, 1.0044 /
data ((supersat( 7,iw,iconc,1, 1),iw=1,9),iconc=1,7)/    &
1.0023, 1.0055, 1.0125, 1.0283, 1.0642, 1.1515, 1.4015, 2.6817, 4.9472,  &
1.0009, 1.0029, 1.0077, 1.0183, 1.0419, 1.0964, 1.2347, 1.6930, 3.9396,  &
1.0003, 1.0010, 1.0038, 1.0109, 1.0268, 1.0626, 1.1476, 1.3815, 2.4659,  &
1.0002, 1.0003, 1.0010, 1.0048, 1.0152, 1.0390, 1.0940, 1.2318, 1.6749,  &
1.0001, 1.0002, 1.0004, 1.0010, 1.0059, 1.0210, 1.0563, 1.1417, 1.3788,  &
0.9999, 0.9999, 1.0001, 1.0005, 1.0013, 1.0070, 1.0278, 1.0794, 1.2147,  &
0.9991, 0.9991, 0.9992, 0.9995, 1.0002, 1.0018, 1.0079, 1.0351, 1.1086 /
data ((supersat( 7,iw,iconc,2, 1),iw=1,9),iconc=1,7)/    &
1.0017, 1.0038, 1.0087, 1.0194, 1.0438, 1.1016, 1.2534, 1.7957, 3.6925,  &
1.0007, 1.0020, 1.0053, 1.0124, 1.0285, 1.0647, 1.1525, 1.4007, 2.5608,  &
1.0003, 1.0008, 1.0026, 1.0074, 1.0181, 1.0420, 1.0971, 1.2362, 1.6891,  &
1.0001, 1.0003, 1.0008, 1.0033, 1.0104, 1.0263, 1.0625, 1.1483, 1.3834,  &
1.0001, 1.0002, 1.0003, 1.0008, 1.0041, 1.0143, 1.0381, 1.0932, 1.2319,  &
1.0000, 1.0001, 1.0002, 1.0004, 1.0008, 1.0047, 1.0193, 1.0540, 1.1390,  &
0.9997, 0.9997, 0.9998, 0.9999, 1.0003, 1.0012, 1.0053, 1.0249, 1.0747 /
data ((supersat( 7,iw,iconc,3, 1),iw=1,9),iconc=1,7)/    &
1.0012, 1.0027, 1.0062, 1.0139, 1.0314, 1.0722, 1.1748, 1.4806, 2.7259,  &
1.0005, 1.0014, 1.0037, 1.0088, 1.0201, 1.0455, 1.1057, 1.2628, 1.8050,  &
1.0002, 1.0006, 1.0019, 1.0051, 1.0126, 1.0292, 1.0670, 1.1583, 1.4145,  &
1.0001, 1.0002, 1.0006, 1.0023, 1.0071, 1.0182, 1.0429, 1.1001, 1.2446,  &
1.0001, 1.0001, 1.0002, 1.0006, 1.0027, 1.0097, 1.0261, 1.0633, 1.1520,  &
1.0000, 1.0001, 1.0002, 1.0003, 1.0006, 1.0031, 1.0131, 1.0370, 1.0932,  &
0.9999, 0.9999, 1.0000, 1.0001, 1.0003, 1.0008, 1.0033, 1.0170, 1.0513 /
data ((supersat( 7,iw,iconc,4, 1),iw=1,9),iconc=1,7)/    &
1.0009, 1.0020, 1.0046, 1.0104, 1.0235, 1.0543, 1.1296, 1.3349, 2.0878,  &
1.0004, 1.0011, 1.0027, 1.0065, 1.0147, 1.0336, 1.0777, 1.1882, 1.5129,  &
1.0002, 1.0005, 1.0014, 1.0037, 1.0091, 1.0211, 1.0485, 1.1133, 1.2821,  &
1.0001, 1.0002, 1.0005, 1.0016, 1.0050, 1.0128, 1.0306, 1.0710, 1.1690,  &
1.0001, 1.0001, 1.0002, 1.0005, 1.0018, 1.0067, 1.0182, 1.0444, 1.1052,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0005, 1.0019, 1.0088, 1.0256, 1.0645,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0002, 1.0006, 1.0020, 1.0113, 1.0353 /
data ((supersat( 7,iw,iconc,5, 1),iw=1,9),iconc=1,7)/    &
1.0006, 1.0016, 1.0036, 1.0081, 1.0185, 1.0429, 1.1020, 1.2557, 1.7352,  &
1.0004, 1.0009, 1.0020, 1.0050, 1.0113, 1.0260, 1.0603, 1.1446, 1.3736,  &
1.0002, 1.0004, 1.0010, 1.0028, 1.0068, 1.0159, 1.0369, 1.0861, 1.2094,  &
1.0001, 1.0002, 1.0004, 1.0012, 1.0036, 1.0094, 1.0227, 1.0529, 1.1248,  &
1.0001, 1.0001, 1.0002, 1.0004, 1.0013, 1.0047, 1.0131, 1.0324, 1.0767,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0004, 1.0013, 1.0060, 1.0182, 1.0463,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0002, 1.0004, 1.0013, 1.0074, 1.0247 /
data ((supersat( 7,iw,iconc,6, 1),iw=1,9),iconc=1,7)/    &
1.0006, 1.0013, 1.0029, 1.0066, 1.0152, 1.0354, 1.0841, 1.2076, 1.5616,  &
1.0003, 1.0007, 1.0016, 1.0039, 1.0091, 1.0210, 1.0491, 1.1172, 1.2951,  &
1.0001, 1.0003, 1.0008, 1.0021, 1.0053, 1.0126, 1.0294, 1.0688, 1.1659,  &
1.0001, 1.0001, 1.0003, 1.0009, 1.0027, 1.0072, 1.0176, 1.0414, 1.0977,  &
1.0001, 1.0001, 1.0001, 1.0003, 1.0009, 1.0034, 1.0098, 1.0247, 1.0589,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0003, 1.0009, 1.0042, 1.0133, 1.0346,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0002, 1.0003, 1.0009, 1.0049, 1.0178 /
data ((supersat( 7,iw,iconc,7, 1),iw=1,9),iconc=1,7)/    &
1.0005, 1.0010, 1.0023, 1.0056, 1.0129, 1.0301, 1.0716, 1.1754, 1.4600,  &
1.0002, 1.0005, 1.0014, 1.0032, 1.0075, 1.0176, 1.0413, 1.0986, 1.2452,  &
1.0001, 1.0003, 1.0007, 1.0017, 1.0043, 1.0103, 1.0243, 1.0572, 1.1375,  &
1.0001, 1.0001, 1.0003, 1.0007, 1.0022, 1.0058, 1.0142, 1.0337, 1.0799,  &
1.0000, 1.0001, 1.0001, 1.0003, 1.0007, 1.0026, 1.0077, 1.0196, 1.0473,  &
1.0000, 1.0001, 1.0001, 1.0001, 1.0003, 1.0007, 1.0030, 1.0102, 1.0271,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0001, 1.0003, 1.0007, 1.0034, 1.0133 /
data ((supersat( 7,iw,iconc,1, 2),iw=1,9),iconc=1,7)/    &
1.0018, 1.0049, 1.0118, 1.0274, 1.0631, 1.1499, 1.3991, 2.6747, 4.9352,  &
1.0006, 1.0022, 1.0066, 1.0171, 1.0403, 1.0942, 1.2313, 1.6863, 3.9162,  &
1.0002, 1.0006, 1.0026, 1.0092, 1.0246, 1.0595, 1.1430, 1.3733, 2.4436,  &
1.0001, 1.0003, 1.0006, 1.0030, 1.0125, 1.0351, 1.0879, 1.2216, 1.6529,  &
1.0000, 1.0001, 1.0003, 1.0008, 1.0032, 1.0164, 1.0490, 1.1293, 1.3546,  &
0.9994, 0.9995, 0.9997, 1.0001, 1.0010, 1.0034, 1.0207, 1.0663, 1.1892,  &
0.9980, 0.9981, 0.9982, 0.9984, 0.9991, 1.0008, 1.0049, 1.0246, 1.0859 /
data ((supersat( 7,iw,iconc,2, 2),iw=1,9),iconc=1,7)/    &
1.0013, 1.0034, 1.0082, 1.0189, 1.0432, 1.1008, 1.2522, 1.7935, 3.6858,  &
1.0005, 1.0016, 1.0047, 1.0116, 1.0275, 1.0634, 1.1507, 1.3978, 2.5525,  &
1.0002, 1.0005, 1.0019, 1.0063, 1.0167, 1.0402, 1.0945, 1.2321, 1.6811,  &
1.0001, 1.0002, 1.0005, 1.0021, 1.0085, 1.0238, 1.0589, 1.1429, 1.3737,  &
1.0000, 1.0001, 1.0002, 1.0005, 1.0022, 1.0113, 1.0336, 1.0863, 1.2199,  &
0.9998, 0.9999, 1.0000, 1.0002, 1.0007, 1.0023, 1.0144, 1.0461, 1.1251,  &
0.9991, 0.9991, 0.9992, 0.9994, 0.9998, 1.0007, 1.0030, 1.0176, 1.0607 /
data ((supersat( 7,iw,iconc,3, 2),iw=1,9),iconc=1,7)/    &
1.0010, 1.0024, 1.0058, 1.0135, 1.0309, 1.0717, 1.1741, 1.4796, 2.7229,  &
1.0004, 1.0012, 1.0032, 1.0082, 1.0194, 1.0447, 1.1047, 1.2613, 1.8021,  &
1.0002, 1.0004, 1.0013, 1.0044, 1.0117, 1.0280, 1.0655, 1.1561, 1.4107,  &
1.0001, 1.0002, 1.0004, 1.0015, 1.0058, 1.0165, 1.0407, 1.0970, 1.2395,  &
1.0001, 1.0001, 1.0002, 1.0004, 1.0015, 1.0076, 1.0231, 1.0590, 1.1453,  &
0.9999, 1.0000, 1.0000, 1.0002, 1.0005, 1.0015, 1.0097, 1.0318, 1.0848,  &
0.9996, 0.9996, 0.9996, 0.9997, 1.0000, 1.0006, 1.0018, 1.0118, 1.0422 /
data ((supersat( 7,iw,iconc,4, 2),iw=1,9),iconc=1,7)/    &
1.0006, 1.0018, 1.0044, 1.0101, 1.0232, 1.0540, 1.1292, 1.3343, 2.0865,  &
1.0003, 1.0009, 1.0023, 1.0060, 1.0142, 1.0330, 1.0771, 1.1873, 1.5115,  &
1.0001, 1.0003, 1.0010, 1.0031, 1.0084, 1.0203, 1.0475, 1.1120, 1.2800,  &
1.0001, 1.0001, 1.0003, 1.0010, 1.0041, 1.0116, 1.0290, 1.0690, 1.1660,  &
1.0001, 1.0001, 1.0001, 1.0003, 1.0010, 1.0051, 1.0160, 1.0415, 1.1011,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0004, 1.0011, 1.0063, 1.0220, 1.0590,  &
0.9998, 0.9998, 0.9998, 0.9999, 1.0000, 1.0004, 1.0012, 1.0075, 1.0291 /
data ((supersat( 7,iw,iconc,5, 2),iw=1,9),iconc=1,7)/    &
1.0005, 1.0015, 1.0034, 1.0079, 1.0183, 1.0427, 1.1017, 1.2553, 1.7346,  &
1.0002, 1.0007, 1.0018, 1.0046, 1.0110, 1.0256, 1.0599, 1.1440, 1.3728,  &
1.0001, 1.0002, 1.0007, 1.0023, 1.0062, 1.0153, 1.0361, 1.0851, 1.2082,  &
1.0001, 1.0001, 1.0002, 1.0007, 1.0029, 1.0085, 1.0215, 1.0515, 1.1229,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0007, 1.0035, 1.0115, 1.0303, 1.0739,  &
1.0000, 1.0000, 1.0001, 1.0001, 1.0003, 1.0008, 1.0041, 1.0155, 1.0424,  &
0.9998, 0.9999, 0.9999, 0.9999, 1.0000, 1.0003, 1.0008, 1.0046, 1.0202 /
data ((supersat( 7,iw,iconc,6, 2),iw=1,9),iconc=1,7)/    &
1.0005, 1.0011, 1.0027, 1.0064, 1.0150, 1.0352, 1.0839, 1.2073, 1.5613,  &
1.0002, 1.0005, 1.0015, 1.0037, 1.0088, 1.0207, 1.0487, 1.1168, 1.2946,  &
1.0001, 1.0002, 1.0006, 1.0018, 1.0048, 1.0121, 1.0288, 1.0681, 1.1650,  &
1.0001, 1.0001, 1.0002, 1.0006, 1.0022, 1.0065, 1.0167, 1.0403, 1.0963,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0006, 1.0025, 1.0086, 1.0231, 1.0568,  &
1.0000, 1.0000, 1.0001, 1.0001, 1.0002, 1.0006, 1.0027, 1.0112, 1.0318,  &
0.9999, 0.9999, 0.9999, 0.9999, 1.0000, 1.0002, 1.0006, 1.0029, 1.0143 /
data ((supersat( 7,iw,iconc,7, 2),iw=1,9),iconc=1,7)/    &
1.0003, 1.0009, 1.0022, 1.0054, 1.0127, 1.0300, 1.0714, 1.1752, 1.4598,  &
1.0002, 1.0004, 1.0012, 1.0030, 1.0073, 1.0174, 1.0410, 1.0983, 1.2448,  &
1.0001, 1.0002, 1.0005, 1.0015, 1.0039, 1.0099, 1.0238, 1.0566, 1.1369,  &
1.0001, 1.0001, 1.0002, 1.0005, 1.0017, 1.0051, 1.0135, 1.0329, 1.0789,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0005, 1.0019, 1.0066, 1.0183, 1.0456,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0002, 1.0005, 1.0019, 1.0084, 1.0248,  &
0.9999, 0.9999, 0.9999, 0.9999, 1.0000, 1.0002, 1.0005, 1.0019, 1.0105 /
data ((supersat( 7,iw,iconc,1, 3),iw=1,9),iconc=1,7)/    &
1.0013, 1.0041, 1.0110, 1.0263, 1.0617, 1.1480, 1.3959, 2.6654, 4.9197,  &
1.0003, 1.0014, 1.0055, 1.0155, 1.0383, 1.0914, 1.2269, 1.6776, 3.8853,  &
1.0002, 1.0004, 1.0015, 1.0073, 1.0219, 1.0557, 1.1371, 1.3627, 2.4145,  &
1.0001, 1.0002, 1.0005, 1.0015, 1.0093, 1.0303, 1.0804, 1.2086, 1.6246,  &
0.9997, 0.9998, 1.0000, 1.0005, 1.0016, 1.0115, 1.0406, 1.1145, 1.3248,  &
0.9987, 0.9987, 0.9989, 0.9993, 1.0003, 1.0026, 1.0135, 1.0521, 1.1598,  &
0.9965, 0.9965, 0.9966, 0.9969, 0.9976, 0.9992, 1.0033, 1.0151, 1.0634 /
data ((supersat( 7,iw,iconc,2, 3),iw=1,9),iconc=1,7)/    &
1.0009, 1.0030, 1.0076, 1.0182, 1.0423, 1.0997, 1.2506, 1.7906, 3.6770,  &
1.0003, 1.0011, 1.0039, 1.0106, 1.0262, 1.0618, 1.1484, 1.3939, 2.5416,  &
1.0001, 1.0003, 1.0011, 1.0050, 1.0150, 1.0379, 1.0913, 1.2269, 1.6707,  &
1.0001, 1.0002, 1.0003, 1.0011, 1.0064, 1.0208, 1.0546, 1.1360, 1.3611,  &
0.9999, 1.0000, 1.0001, 1.0004, 1.0012, 1.0079, 1.0283, 1.0777, 1.2050,  &
0.9994, 0.9995, 0.9996, 0.9998, 1.0004, 1.0016, 1.0094, 1.0371, 1.1086,  &
0.9981, 0.9981, 0.9982, 0.9984, 0.9988, 0.9998, 1.0022, 1.0107, 1.0461 /
data ((supersat( 7,iw,iconc,3, 3),iw=1,9),iconc=1,7)/    &
1.0007, 1.0021, 1.0054, 1.0131, 1.0304, 1.0711, 1.1732, 1.4782, 2.7187,  &
1.0002, 1.0008, 1.0028, 1.0076, 1.0186, 1.0437, 1.1033, 1.2593, 1.7983,  &
1.0001, 1.0002, 1.0008, 1.0035, 1.0105, 1.0265, 1.0635, 1.1532, 1.4057,  &
1.0001, 1.0001, 1.0002, 1.0008, 1.0044, 1.0145, 1.0379, 1.0929, 1.2329,  &
1.0000, 1.0000, 1.0001, 1.0003, 1.0008, 1.0053, 1.0196, 1.0537, 1.1368,  &
0.9997, 0.9998, 0.9998, 1.0000, 1.0003, 1.0011, 1.0061, 1.0258, 1.0748,  &
0.9990, 0.9990, 0.9990, 0.9991, 0.9994, 1.0000, 1.0014, 1.0068, 1.0325 /
data ((supersat( 7,iw,iconc,4, 3),iw=1,9),iconc=1,7)/    &
1.0005, 1.0016, 1.0041, 1.0098, 1.0229, 1.0535, 1.1286, 1.3335, 2.0849,  &
1.0002, 1.0006, 1.0020, 1.0055, 1.0137, 1.0323, 1.0762, 1.1862, 1.5096,  &
1.0001, 1.0002, 1.0006, 1.0025, 1.0075, 1.0192, 1.0462, 1.1102, 1.2772,  &
1.0001, 1.0001, 1.0002, 1.0006, 1.0030, 1.0102, 1.0271, 1.0664, 1.1621,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0006, 1.0034, 1.0136, 1.0380, 1.0958,  &
0.9999, 0.9999, 0.9999, 1.0000, 1.0002, 1.0007, 1.0038, 1.0178, 1.0524,  &
0.9994, 0.9994, 0.9994, 0.9995, 0.9997, 1.0001, 1.0009, 1.0040, 1.0224 /
data ((supersat( 7,iw,iconc,5, 3),iw=1,9),iconc=1,7)/    &
1.0004, 1.0012, 1.0032, 1.0077, 1.0180, 1.0424, 1.1014, 1.2549, 1.7338,  &
1.0002, 1.0005, 1.0015, 1.0042, 1.0105, 1.0251, 1.0593, 1.1433, 1.3716,  &
1.0001, 1.0002, 1.0005, 1.0018, 1.0056, 1.0145, 1.0352, 1.0839, 1.2064,  &
1.0000, 1.0001, 1.0002, 1.0005, 1.0021, 1.0074, 1.0201, 1.0497, 1.1204,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0005, 1.0023, 1.0097, 1.0278, 1.0703,  &
0.9999, 0.9999, 1.0000, 1.0000, 1.0002, 1.0005, 1.0024, 1.0124, 1.0378,  &
0.9996, 0.9996, 0.9996, 0.9996, 0.9998, 1.0000, 1.0006, 1.0024, 1.0153 /
data ((supersat( 7,iw,iconc,6, 3),iw=1,9),iconc=1,7)/    &
1.0003, 1.0009, 1.0025, 1.0062, 1.0148, 1.0350, 1.0836, 1.2070, 1.5608,  &
1.0001, 1.0004, 1.0012, 1.0033, 1.0084, 1.0203, 1.0483, 1.1162, 1.2938,  &
1.0001, 1.0001, 1.0004, 1.0014, 1.0044, 1.0115, 1.0281, 1.0672, 1.1638,  &
1.0000, 1.0001, 1.0001, 1.0004, 1.0015, 1.0056, 1.0156, 1.0390, 1.0945,  &
1.0000, 1.0000, 1.0001, 1.0001, 1.0004, 1.0016, 1.0071, 1.0211, 1.0542,  &
0.9999, 0.9999, 1.0000, 1.0000, 1.0001, 1.0004, 1.0016, 1.0088, 1.0283,  &
0.9996, 0.9996, 0.9997, 0.9997, 0.9998, 1.0000, 1.0004, 1.0017, 1.0106 /
data ((supersat( 7,iw,iconc,7, 3),iw=1,9),iconc=1,7)/    &
1.0003, 1.0008, 1.0021, 1.0052, 1.0125, 1.0298, 1.0712, 1.1749, 1.4594,  &
1.0001, 1.0003, 1.0009, 1.0027, 1.0070, 1.0170, 1.0407, 1.0979, 1.2443,  &
1.0001, 1.0001, 1.0003, 1.0011, 1.0035, 1.0094, 1.0232, 1.0560, 1.1360,  &
1.0000, 1.0001, 1.0001, 1.0003, 1.0012, 1.0044, 1.0126, 1.0318, 1.0775,  &
1.0000, 1.0000, 1.0001, 1.0001, 1.0003, 1.0012, 1.0054, 1.0167, 1.0436,  &
0.9999, 0.9999, 1.0000, 1.0000, 1.0001, 1.0003, 1.0012, 1.0065, 1.0220,  &
0.9997, 0.9997, 0.9997, 0.9997, 0.9998, 1.0000, 1.0003, 1.0012, 1.0075 /
data ((supersat( 7,iw,iconc,1, 4),iw=1,9),iconc=1,7)/    &
1.0008, 1.0033, 1.0099, 1.0250, 1.0600, 1.1455, 1.3917, 2.6530, 4.8988,  &
1.0002, 1.0008, 1.0042, 1.0137, 1.0359, 1.0879, 1.2213, 1.6663, 3.8447,  &
1.0001, 1.0003, 1.0008, 1.0051, 1.0188, 1.0510, 1.1297, 1.3492, 2.3768,  &
0.9999, 1.0000, 1.0003, 1.0010, 1.0060, 1.0249, 1.0714, 1.1926, 1.5893,  &
0.9991, 0.9992, 0.9995, 1.0001, 1.0014, 1.0068, 1.0317, 1.0975, 1.2892,  &
0.9974, 0.9975, 0.9977, 0.9981, 0.9991, 1.0015, 1.0073, 1.0382, 1.1286,  &
0.9947, 0.9947, 0.9948, 0.9950, 0.9956, 0.9971, 1.0009, 1.0115, 1.0437 /
data ((supersat( 7,iw,iconc,2, 4),iw=1,9),iconc=1,7)/    &
1.0006, 1.0024, 1.0069, 1.0173, 1.0413, 1.0983, 1.2486, 1.7869, 3.6653,  &
1.0002, 1.0006, 1.0030, 1.0095, 1.0247, 1.0598, 1.1454, 1.3889, 2.5273,  &
1.0001, 1.0002, 1.0006, 1.0036, 1.0130, 1.0351, 1.0872, 1.2203, 1.6572,  &
1.0000, 1.0001, 1.0003, 1.0007, 1.0042, 1.0174, 1.0493, 1.1275, 1.3453,  &
0.9997, 0.9997, 0.9999, 1.0002, 1.0009, 1.0046, 1.0225, 1.0678, 1.1867,  &
0.9987, 0.9988, 0.9989, 0.9991, 0.9997, 1.0011, 1.0049, 1.0278, 1.0903,  &
0.9967, 0.9967, 0.9967, 0.9969, 0.9973, 0.9983, 1.0008, 1.0072, 1.0325 /
data ((supersat( 7,iw,iconc,3, 4),iw=1,9),iconc=1,7)/    &
1.0005, 1.0018, 1.0049, 1.0125, 1.0297, 1.0703, 1.1721, 1.4764, 2.7133,  &
1.0002, 1.0005, 1.0021, 1.0068, 1.0176, 1.0424, 1.1016, 1.2567, 1.7933,  &
1.0001, 1.0002, 1.0005, 1.0025, 1.0091, 1.0247, 1.0610, 1.1495, 1.3993,  &
1.0000, 1.0001, 1.0002, 1.0005, 1.0028, 1.0121, 1.0345, 1.0878, 1.2245,  &
0.9999, 0.9999, 1.0000, 1.0002, 1.0006, 1.0030, 1.0157, 1.0475, 1.1264,  &
0.9993, 0.9994, 0.9994, 0.9996, 1.0000, 1.0008, 1.0031, 1.0195, 1.0634,  &
0.9980, 0.9980, 0.9980, 0.9981, 0.9984, 0.9991, 1.0006, 1.0044, 1.0231 /
data ((supersat( 7,iw,iconc,4, 4),iw=1,9),iconc=1,7)/    &
1.0004, 1.0013, 1.0037, 1.0093, 1.0224, 1.0530, 1.1279, 1.3325, 2.0827,  &
1.0001, 1.0004, 1.0016, 1.0049, 1.0129, 1.0315, 1.0751, 1.1846, 1.5070,  &
1.0001, 1.0001, 1.0004, 1.0018, 1.0065, 1.0179, 1.0445, 1.1079, 1.2737,  &
1.0000, 1.0001, 1.0001, 1.0004, 1.0019, 1.0085, 1.0248, 1.0631, 1.1571,  &
0.9999, 1.0000, 1.0000, 1.0002, 1.0004, 1.0019, 1.0108, 1.0338, 1.0892,  &
0.9996, 0.9996, 0.9997, 0.9998, 1.0000, 1.0005, 1.0020, 1.0133, 1.0449,  &
0.9987, 0.9987, 0.9987, 0.9988, 0.9990, 0.9994, 1.0004, 1.0027, 1.0157 /
data ((supersat( 7,iw,iconc,5, 4),iw=1,9),iconc=1,7)/    &
1.0003, 1.0010, 1.0028, 1.0073, 1.0177, 1.0420, 1.1009, 1.2542, 1.7327,  &
1.0001, 1.0003, 1.0012, 1.0037, 1.0100, 1.0244, 1.0585, 1.1422, 1.3701,  &
1.0001, 1.0001, 1.0003, 1.0013, 1.0048, 1.0136, 1.0340, 1.0824, 1.2042,  &
1.0000, 1.0001, 1.0001, 1.0003, 1.0013, 1.0061, 1.0184, 1.0474, 1.1171,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0003, 1.0013, 1.0075, 1.0247, 1.0658,  &
0.9997, 0.9997, 0.9998, 0.9998, 1.0000, 1.0004, 1.0014, 1.0090, 1.0324,  &
0.9990, 0.9990, 0.9990, 0.9991, 0.9992, 0.9995, 1.0002, 1.0017, 1.0104 /
data ((supersat( 7,iw,iconc,6, 4),iw=1,9),iconc=1,7)/    &
1.0002, 1.0008, 1.0023, 1.0060, 1.0145, 1.0347, 1.0832, 1.2065, 1.5601,  &
1.0001, 1.0002, 1.0009, 1.0029, 1.0080, 1.0198, 1.0477, 1.1155, 1.2928,  &
1.0001, 1.0001, 1.0002, 1.0010, 1.0037, 1.0107, 1.0271, 1.0661, 1.1623,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0010, 1.0045, 1.0142, 1.0372, 1.0921,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0003, 1.0010, 1.0054, 1.0188, 1.0509,  &
0.9998, 0.9998, 0.9998, 0.9999, 1.0000, 1.0003, 1.0010, 1.0062, 1.0243,  &
0.9992, 0.9992, 0.9992, 0.9992, 0.9993, 0.9996, 1.0001, 1.0012, 1.0068 /
data ((supersat( 7,iw,iconc,7, 4),iw=1,9),iconc=1,7)/    &
1.0002, 1.0007, 1.0019, 1.0050, 1.0123, 1.0295, 1.0709, 1.1746, 1.4589,  &
1.0001, 1.0002, 1.0007, 1.0024, 1.0066, 1.0166, 1.0402, 1.0973, 1.2435,  &
1.0000, 1.0001, 1.0002, 1.0007, 1.0030, 1.0087, 1.0225, 1.0550, 1.1348,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0007, 1.0035, 1.0114, 1.0304, 1.0757,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0002, 1.0008, 1.0040, 1.0148, 1.0410,  &
0.9998, 0.9998, 0.9998, 0.9999, 1.0000, 1.0002, 1.0008, 1.0044, 1.0188,  &
0.9992, 0.9992, 0.9993, 0.9993, 0.9994, 0.9996, 1.0000, 1.0009, 1.0046 /
data ((supersat( 7,iw,iconc,1, 5),iw=1,9),iconc=1,7)/    &
1.0006, 1.0028, 1.0092, 1.0241, 1.0588, 1.1438, 1.3887, 2.6440, 4.8835,  &
1.0002, 1.0006, 1.0034, 1.0125, 1.0342, 1.0854, 1.2173, 1.6582, 3.8154,  &
1.0001, 1.0002, 1.0006, 1.0038, 1.0167, 1.0479, 1.1246, 1.3397, 2.3499,  &
0.9997, 0.9998, 1.0002, 1.0009, 1.0042, 1.0215, 1.0656, 1.1818, 1.5648,  &
0.9987, 0.9988, 0.9990, 0.9996, 1.0010, 1.0044, 1.0264, 1.0869, 1.2661,  &
0.9965, 0.9966, 0.9967, 0.9971, 0.9981, 1.0006, 1.0068, 1.0307, 1.1105,  &
0.9937, 0.9937, 0.9938, 0.9940, 0.9945, 0.9959, 0.9994, 1.0090, 1.0341 /
data ((supersat( 7,iw,iconc,2, 5),iw=1,9),iconc=1,7)/    &
1.0005, 1.0020, 1.0065, 1.0168, 1.0406, 1.0974, 1.2471, 1.7841, 3.6567,  &
1.0002, 1.0005, 1.0024, 1.0087, 1.0237, 1.0584, 1.1433, 1.3853, 2.5168,  &
1.0001, 1.0002, 1.0005, 1.0027, 1.0116, 1.0332, 1.0843, 1.2156, 1.6475,  &
0.9999, 1.0000, 1.0002, 1.0006, 1.0029, 1.0151, 1.0457, 1.1216, 1.3342,  &
0.9994, 0.9995, 0.9996, 1.0000, 1.0008, 1.0030, 1.0189, 1.0614, 1.1746,  &
0.9981, 0.9981, 0.9982, 0.9985, 0.9991, 1.0006, 1.0043, 1.0225, 1.0792,  &
0.9957, 0.9957, 0.9957, 0.9959, 0.9963, 0.9973, 0.9997, 1.0059, 1.0255 /
data ((supersat( 7,iw,iconc,3, 5),iw=1,9),iconc=1,7)/    &
1.0004, 1.0015, 1.0046, 1.0121, 1.0292, 1.0697, 1.1713, 1.4750, 2.7093,  &
1.0001, 1.0004, 1.0017, 1.0062, 1.0168, 1.0415, 1.1004, 1.2548, 1.7897,  &
1.0001, 1.0001, 1.0004, 1.0019, 1.0081, 1.0234, 1.0592, 1.1468, 1.3946,  &
1.0000, 1.0000, 1.0002, 1.0004, 1.0020, 1.0105, 1.0323, 1.0843, 1.2186,  &
0.9997, 0.9998, 0.9999, 1.0001, 1.0005, 1.0020, 1.0132, 1.0434, 1.1193,  &
0.9990, 0.9990, 0.9991, 0.9992, 0.9996, 1.0005, 1.0026, 1.0158, 1.0563,  &
0.9972, 0.9972, 0.9972, 0.9973, 0.9976, 0.9983, 0.9999, 1.0037, 1.0180 /
data ((supersat( 7,iw,iconc,4, 5),iw=1,9),iconc=1,7)/    &
1.0003, 1.0011, 1.0034, 1.0090, 1.0220, 1.0526, 1.1274, 1.3318, 2.0810,  &
1.0001, 1.0003, 1.0013, 1.0045, 1.0124, 1.0308, 1.0743, 1.1835, 1.5051,  &
1.0001, 1.0001, 1.0003, 1.0013, 1.0058, 1.0170, 1.0433, 1.1062, 1.2711,  &
1.0000, 1.0000, 1.0001, 1.0003, 1.0014, 1.0073, 1.0232, 1.0608, 1.1536,  &
0.9998, 0.9999, 0.9999, 1.0001, 1.0004, 1.0014, 1.0090, 1.0310, 1.0847,  &
0.9994, 0.9994, 0.9994, 0.9995, 0.9998, 1.0004, 1.0017, 1.0107, 1.0401,  &
0.9981, 0.9981, 0.9981, 0.9982, 0.9984, 0.9988, 0.9998, 1.0022, 1.0120 /
data ((supersat( 7,iw,iconc,5, 5),iw=1,9),iconc=1,7)/    &
1.0002, 1.0008, 1.0026, 1.0071, 1.0174, 1.0417, 1.1005, 1.2537, 1.7319,  &
1.0001, 1.0002, 1.0010, 1.0034, 1.0096, 1.0240, 1.0580, 1.1415, 1.3690,  &
1.0000, 1.0001, 1.0002, 1.0010, 1.0043, 1.0129, 1.0331, 1.0813, 1.2025,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0010, 1.0052, 1.0172, 1.0458, 1.1147,  &
0.9999, 0.9999, 1.0000, 1.0001, 1.0003, 1.0010, 1.0062, 1.0227, 1.0627,  &
0.9995, 0.9996, 0.9996, 0.9997, 0.9998, 1.0002, 1.0011, 1.0070, 1.0290,  &
0.9985, 0.9985, 0.9986, 0.9986, 0.9987, 0.9991, 0.9998, 1.0014, 1.0077 /
data ((supersat( 7,iw,iconc,6, 5),iw=1,9),iconc=1,7)/    &
1.0002, 1.0007, 1.0021, 1.0058, 1.0143, 1.0344, 1.0830, 1.2062, 1.5595,  &
1.0001, 1.0002, 1.0007, 1.0027, 1.0077, 1.0194, 1.0473, 1.1150, 1.2920,  &
1.0000, 1.0001, 1.0002, 1.0007, 1.0033, 1.0101, 1.0265, 1.0652, 1.1611,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0007, 1.0039, 1.0133, 1.0360, 1.0904,  &
0.9999, 0.9999, 1.0000, 1.0001, 1.0002, 1.0008, 1.0044, 1.0172, 1.0486,  &
0.9996, 0.9996, 0.9997, 0.9997, 0.9999, 1.0002, 1.0008, 1.0047, 1.0216,  &
0.9987, 0.9987, 0.9988, 0.9988, 0.9989, 0.9992, 0.9997, 1.0009, 1.0049 /
data ((supersat( 7,iw,iconc,7, 5),iw=1,9),iconc=1,7)/    &
1.0002, 1.0006, 1.0018, 1.0049, 1.0121, 1.0293, 1.0707, 1.1743, 1.4585,  &
1.0001, 1.0002, 1.0006, 1.0022, 1.0064, 1.0163, 1.0398, 1.0969, 1.2429,  &
1.0000, 1.0001, 1.0002, 1.0006, 1.0026, 1.0083, 1.0219, 1.0544, 1.1339,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0006, 1.0030, 1.0106, 1.0294, 1.0743,  &
0.9999, 0.9999, 1.0000, 1.0000, 1.0002, 1.0006, 1.0032, 1.0135, 1.0391,  &
0.9997, 0.9997, 0.9997, 0.9997, 0.9999, 1.0001, 1.0006, 1.0033, 1.0167,  &
0.9989, 0.9989, 0.9989, 0.9989, 0.9990, 0.9992, 0.9996, 1.0006, 1.0034 /
data ((supersat( 7,iw,iconc,1, 6),iw=1,9),iconc=1,7)/    &
1.0005, 1.0024, 1.0087, 1.0234, 1.0578, 1.1423, 1.3862, 2.6366, 4.8711,  &
1.0002, 1.0005, 1.0028, 1.0116, 1.0329, 1.0835, 1.2140, 1.6516, 3.7916,  &
1.0000, 1.0002, 1.0006, 1.0030, 1.0151, 1.0454, 1.1205, 1.3321, 2.3282,  &
0.9995, 0.9997, 1.0000, 1.0008, 1.0031, 1.0190, 1.0612, 1.1735, 1.5457,  &
0.9982, 0.9983, 0.9986, 0.9992, 1.0007, 1.0042, 1.0228, 1.0792, 1.2488,  &
0.9958, 0.9958, 0.9960, 0.9964, 0.9974, 0.9998, 1.0059, 1.0258, 1.0981,  &
0.9930, 0.9930, 0.9931, 0.9933, 0.9938, 0.9950, 0.9982, 1.0072, 1.0353 /
data ((supersat( 7,iw,iconc,2, 6),iw=1,9),iconc=1,7)/    &
1.0004, 1.0018, 1.0061, 1.0163, 1.0400, 1.0966, 1.2459, 1.7819, 3.6497,  &
1.0001, 1.0004, 1.0020, 1.0081, 1.0228, 1.0573, 1.1417, 1.3824, 2.5083,  &
1.0001, 1.0002, 1.0004, 1.0021, 1.0106, 1.0317, 1.0820, 1.2118, 1.6397,  &
0.9998, 0.9999, 1.0001, 1.0006, 1.0022, 1.0135, 1.0430, 1.1170, 1.3254,  &
0.9992, 0.9992, 0.9994, 0.9998, 1.0006, 1.0024, 1.0164, 1.0566, 1.1653,  &
0.9976, 0.9976, 0.9977, 0.9980, 0.9986, 1.0001, 1.0038, 1.0190, 1.0713,  &
0.9949, 0.9949, 0.9950, 0.9951, 0.9955, 0.9964, 0.9988, 1.0048, 1.0210 /
data ((supersat( 7,iw,iconc,3, 6),iw=1,9),iconc=1,7)/    &
1.0003, 1.0013, 1.0044, 1.0118, 1.0289, 1.0692, 1.1706, 1.4739, 2.7061,  &
1.0001, 1.0003, 1.0014, 1.0057, 1.0163, 1.0408, 1.0994, 1.2532, 1.7867,  &
1.0001, 1.0001, 1.0003, 1.0015, 1.0074, 1.0224, 1.0578, 1.1447, 1.3909,  &
0.9999, 1.0000, 1.0001, 1.0004, 1.0015, 1.0094, 1.0305, 1.0816, 1.2139,  &
0.9996, 0.9996, 0.9997, 1.0000, 1.0005, 1.0017, 1.0114, 1.0403, 1.1139,  &
0.9986, 0.9986, 0.9987, 0.9989, 0.9993, 1.0002, 1.0024, 1.0133, 1.0511,  &
0.9965, 0.9965, 0.9966, 0.9967, 0.9969, 0.9976, 0.9992, 1.0031, 1.0147 /
data ((supersat( 7,iw,iconc,4, 6),iw=1,9),iconc=1,7)/    &
1.0002, 1.0010, 1.0032, 1.0088, 1.0218, 1.0523, 1.1270, 1.3312, 2.0797,  &
1.0001, 1.0002, 1.0011, 1.0042, 1.0120, 1.0304, 1.0737, 1.1826, 1.5035,  &
1.0001, 1.0001, 1.0002, 1.0011, 1.0053, 1.0163, 1.0424, 1.1049, 1.2689,  &
1.0000, 1.0000, 1.0001, 1.0003, 1.0011, 1.0065, 1.0220, 1.0590, 1.1507,  &
0.9998, 0.9998, 0.9999, 1.0000, 1.0003, 1.0012, 1.0077, 1.0289, 1.0812,  &
0.9991, 0.9991, 0.9992, 0.9993, 0.9996, 1.0002, 1.0015, 1.0088, 1.0366,  &
0.9976, 0.9976, 0.9976, 0.9977, 0.9978, 0.9983, 0.9994, 1.0018, 1.0096 /
data ((supersat( 7,iw,iconc,5, 6),iw=1,9),iconc=1,7)/    &
1.0002, 1.0008, 1.0025, 1.0069, 1.0172, 1.0415, 1.1003, 1.2533, 1.7312,  &
1.0001, 1.0002, 1.0008, 1.0032, 1.0093, 1.0236, 1.0575, 1.1409, 1.3681,  &
1.0000, 1.0001, 1.0002, 1.0008, 1.0039, 1.0123, 1.0325, 1.0803, 1.2012,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0008, 1.0046, 1.0163, 1.0445, 1.1128,  &
0.9998, 0.9999, 0.9999, 1.0000, 1.0002, 1.0008, 1.0052, 1.0211, 1.0603,  &
0.9994, 0.9994, 0.9994, 0.9995, 0.9997, 1.0001, 1.0010, 1.0057, 1.0265,  &
0.9981, 0.9981, 0.9981, 0.9982, 0.9983, 0.9986, 0.9994, 1.0011, 1.0060 /
data ((supersat( 7,iw,iconc,6, 6),iw=1,9),iconc=1,7)/    &
1.0002, 1.0006, 1.0020, 1.0056, 1.0141, 1.0343, 1.0828, 1.2059, 1.5591,  &
1.0001, 1.0002, 1.0006, 1.0025, 1.0074, 1.0191, 1.0469, 1.1145, 1.2914,  &
1.0000, 1.0001, 1.0002, 1.0006, 1.0030, 1.0097, 1.0259, 1.0645, 1.1601,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0006, 1.0034, 1.0126, 1.0350, 1.0890,  &
0.9999, 0.9999, 0.9999, 1.0000, 1.0002, 1.0006, 1.0036, 1.0160, 1.0468,  &
0.9995, 0.9995, 0.9995, 0.9996, 0.9997, 1.0000, 1.0007, 1.0037, 1.0197,  &
0.9984, 0.9984, 0.9984, 0.9984, 0.9985, 0.9988, 0.9994, 1.0006, 1.0039 /
data ((supersat( 7,iw,iconc,7, 6),iw=1,9),iconc=1,7)/    &
1.0001, 1.0005, 1.0017, 1.0047, 1.0120, 1.0292, 1.0705, 1.1741, 1.4582,  &
1.0001, 1.0001, 1.0005, 1.0020, 1.0061, 1.0160, 1.0396, 1.0966, 1.2424,  &
1.0000, 1.0001, 1.0001, 1.0005, 1.0023, 1.0079, 1.0215, 1.0538, 1.1332,  &
1.0000, 1.0000, 1.0001, 1.0001, 1.0005, 1.0026, 1.0100, 1.0286, 1.0733,  &
0.9999, 0.9999, 0.9999, 1.0000, 1.0002, 1.0005, 1.0026, 1.0125, 1.0377,  &
0.9995, 0.9995, 0.9995, 0.9996, 0.9997, 1.0000, 1.0006, 1.0027, 1.0151,  &
0.9985, 0.9985, 0.9985, 0.9985, 0.9986, 0.9988, 0.9993, 1.0004, 1.0028 /
data ((supersat( 7,iw,iconc,1, 7),iw=1,9),iconc=1,7)/    &
1.0004, 1.0021, 1.0082, 1.0228, 1.0570, 1.1411, 1.3841, 2.6302, 4.8603,  &
1.0002, 1.0004, 1.0023, 1.0108, 1.0318, 1.0818, 1.2113, 1.6460, 3.7712,  &
0.9999, 1.0001, 1.0006, 1.0024, 1.0139, 1.0435, 1.1172, 1.3257, 2.3098,  &
0.9993, 0.9995, 0.9999, 1.0007, 1.0025, 1.0171, 1.0576, 1.1667, 1.5299,  &
0.9978, 0.9979, 0.9982, 0.9988, 1.0003, 1.0039, 1.0201, 1.0733, 1.2350,  &
0.9952, 0.9953, 0.9954, 0.9958, 0.9968, 0.9991, 1.0051, 1.0223, 1.0889,  &
0.9926, 0.9926, 0.9926, 0.9928, 0.9933, 0.9944, 0.9974, 1.0058, 1.0320 /
data ((supersat( 7,iw,iconc,2, 7),iw=1,9),iconc=1,7)/    &
1.0003, 1.0016, 1.0058, 1.0159, 1.0395, 1.0959, 1.2449, 1.7799, 3.6436,  &
1.0001, 1.0003, 1.0017, 1.0076, 1.0221, 1.0563, 1.1402, 1.3799, 2.5009,  &
1.0000, 1.0001, 1.0004, 1.0017, 1.0097, 1.0304, 1.0801, 1.2086, 1.6331,  &
0.9997, 0.9998, 1.0000, 1.0005, 1.0018, 1.0122, 1.0408, 1.1132, 1.3181,  &
0.9989, 0.9990, 0.9992, 0.9996, 1.0004, 1.0025, 1.0145, 1.0529, 1.1578,  &
0.9971, 0.9972, 0.9972, 0.9975, 0.9982, 0.9997, 1.0034, 1.0164, 1.0653,  &
0.9943, 0.9944, 0.9944, 0.9945, 0.9949, 0.9958, 0.9981, 1.0039, 1.0179 /
data ((supersat( 7,iw,iconc,3, 7),iw=1,9),iconc=1,7)/    &
1.0003, 1.0012, 1.0042, 1.0115, 1.0285, 1.0688, 1.1701, 1.4730, 2.7032,  &
1.0001, 1.0003, 1.0012, 1.0054, 1.0158, 1.0402, 1.0986, 1.2519, 1.7841,  &
1.0000, 1.0001, 1.0003, 1.0013, 1.0068, 1.0216, 1.0566, 1.1429, 1.3877,  &
0.9999, 1.0000, 1.0001, 1.0004, 1.0013, 1.0084, 1.0290, 1.0793, 1.2099,  &
0.9995, 0.9995, 0.9996, 0.9998, 1.0004, 1.0016, 1.0100, 1.0378, 1.1094,  &
0.9983, 0.9983, 0.9984, 0.9986, 0.9990, 0.9999, 1.0022, 1.0114, 1.0472,  &
0.9960, 0.9960, 0.9960, 0.9961, 0.9964, 0.9970, 0.9986, 1.0025, 1.0124 /
data ((supersat( 7,iw,iconc,4, 7),iw=1,9),iconc=1,7)/    &
1.0002, 1.0009, 1.0031, 1.0086, 1.0216, 1.0520, 1.1267, 1.3306, 2.0785,  &
1.0001, 1.0002, 1.0009, 1.0039, 1.0117, 1.0299, 1.0731, 1.1818, 1.5022,  &
1.0000, 1.0001, 1.0002, 1.0009, 1.0048, 1.0158, 1.0416, 1.1038, 1.2671,  &
1.0000, 1.0000, 1.0001, 1.0003, 1.0009, 1.0058, 1.0210, 1.0575, 1.1483,  &
0.9997, 0.9997, 0.9998, 0.9999, 1.0003, 1.0010, 1.0067, 1.0272, 1.0784,  &
0.9989, 0.9989, 0.9990, 0.9991, 0.9994, 1.0000, 1.0014, 1.0074, 1.0338,  &
0.9971, 0.9971, 0.9971, 0.9972, 0.9974, 0.9978, 0.9989, 1.0014, 1.0079 /
data ((supersat( 7,iw,iconc,5, 7),iw=1,9),iconc=1,7)/    &
1.0002, 1.0007, 1.0024, 1.0068, 1.0170, 1.0413, 1.1000, 1.2530, 1.7307,  &
1.0001, 1.0002, 1.0007, 1.0030, 1.0090, 1.0233, 1.0571, 1.1404, 1.3672,  &
1.0000, 1.0001, 1.0002, 1.0007, 1.0035, 1.0119, 1.0319, 1.0795, 1.2000,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0007, 1.0041, 1.0156, 1.0434, 1.1112,  &
0.9998, 0.9998, 0.9998, 1.0000, 1.0002, 1.0007, 1.0045, 1.0199, 1.0583,  &
0.9992, 0.9992, 0.9992, 0.9993, 0.9995, 1.0000, 1.0009, 1.0047, 1.0245,  &
0.9977, 0.9977, 0.9977, 0.9978, 0.9979, 0.9982, 0.9990, 1.0007, 1.0050 /
data ((supersat( 7,iw,iconc,6, 7),iw=1,9),iconc=1,7)/    &
1.0001, 1.0005, 1.0019, 1.0055, 1.0140, 1.0341, 1.0826, 1.2057, 1.5587,  &
1.0001, 1.0001, 1.0005, 1.0023, 1.0072, 1.0189, 1.0466, 1.1141, 1.2908,  &
1.0000, 1.0001, 1.0002, 1.0005, 1.0027, 1.0093, 1.0255, 1.0640, 1.1593,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0005, 1.0030, 1.0120, 1.0342, 1.0879,  &
0.9998, 0.9998, 0.9999, 1.0000, 1.0002, 1.0006, 1.0031, 1.0150, 1.0453,  &
0.9993, 0.9993, 0.9993, 0.9994, 0.9996, 0.9999, 1.0006, 1.0032, 1.0182,  &
0.9980, 0.9980, 0.9980, 0.9981, 0.9982, 0.9984, 0.9990, 1.0004, 1.0035 /
data ((supersat( 7,iw,iconc,7, 7),iw=1,9),iconc=1,7)/    &
1.0001, 1.0004, 1.0016, 1.0046, 1.0119, 1.0290, 1.0704, 1.1739, 1.4579,  &
1.0001, 1.0001, 1.0004, 1.0019, 1.0060, 1.0158, 1.0393, 1.0963, 1.2420,  &
1.0000, 1.0001, 1.0001, 1.0004, 1.0021, 1.0076, 1.0211, 1.0533, 1.1325,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0004, 1.0022, 1.0095, 1.0279, 1.0723,  &
0.9998, 0.9998, 0.9999, 1.0000, 1.0001, 1.0005, 1.0023, 1.0117, 1.0365,  &
0.9994, 0.9994, 0.9994, 0.9995, 0.9996, 0.9999, 1.0005, 1.0023, 1.0139,  &
0.9982, 0.9982, 0.9982, 0.9982, 0.9983, 0.9985, 0.9990, 1.0001, 1.0026 /
data ((supersat( 8,iw,iconc,1, 1),iw=1,9),iconc=1,7)/    &
1.0017, 1.0047, 1.0116, 1.0272, 1.0629, 1.1496, 1.3985, 2.6729, 4.9320,  &
1.0005, 1.0020, 1.0064, 1.0168, 1.0399, 1.0937, 1.2304, 1.6847, 3.9100,  &
1.0002, 1.0005, 1.0023, 1.0088, 1.0240, 1.0587, 1.1418, 1.3712, 2.4379,  &
1.0001, 1.0003, 1.0005, 1.0026, 1.0118, 1.0340, 1.0862, 1.2188, 1.6470,  &
0.9999, 1.0000, 1.0002, 1.0007, 1.0027, 1.0152, 1.0470, 1.1258, 1.3479,  &
0.9993, 0.9993, 0.9995, 0.9999, 1.0008, 1.0029, 1.0188, 1.0626, 1.1818,  &
0.9977, 0.9977, 0.9978, 0.9981, 0.9988, 1.0005, 1.0046, 1.0218, 1.0795 /
data ((supersat( 8,iw,iconc,2, 1),iw=1,9),iconc=1,7)/    &
1.0012, 1.0033, 1.0081, 1.0187, 1.0430, 1.1006, 1.2519, 1.7930, 3.6841,  &
1.0004, 1.0015, 1.0045, 1.0114, 1.0272, 1.0631, 1.1503, 1.3971, 2.5505,  &
1.0002, 1.0004, 1.0017, 1.0060, 1.0163, 1.0397, 1.0939, 1.2311, 1.6791,  &
1.0001, 1.0002, 1.0004, 1.0018, 1.0080, 1.0232, 1.0580, 1.1415, 1.3712,  &
1.0000, 1.0001, 1.0002, 1.0005, 1.0019, 1.0105, 1.0324, 1.0844, 1.2167,  &
0.9997, 0.9998, 0.9999, 1.0001, 1.0006, 1.0020, 1.0132, 1.0439, 1.1212,  &
0.9989, 0.9989, 0.9990, 0.9992, 0.9996, 1.0005, 1.0029, 1.0157, 1.0568 /
data ((supersat( 8,iw,iconc,3, 1),iw=1,9),iconc=1,7)/    &
1.0009, 1.0024, 1.0057, 1.0134, 1.0308, 1.0716, 1.1740, 1.4793, 2.7221,  &
1.0003, 1.0011, 1.0031, 1.0081, 1.0192, 1.0445, 1.1044, 1.2609, 1.8015,  &
1.0001, 1.0003, 1.0012, 1.0042, 1.0114, 1.0277, 1.0651, 1.1556, 1.4098,  &
1.0001, 1.0001, 1.0003, 1.0013, 1.0055, 1.0161, 1.0401, 1.0962, 1.2382,  &
1.0000, 1.0001, 1.0002, 1.0004, 1.0013, 1.0071, 1.0224, 1.0579, 1.1436,  &
0.9999, 0.9999, 1.0000, 1.0002, 1.0005, 1.0014, 1.0089, 1.0304, 1.0826,  &
0.9995, 0.9995, 0.9995, 0.9996, 0.9999, 1.0005, 1.0018, 1.0106, 1.0398 /
data ((supersat( 8,iw,iconc,4, 1),iw=1,9),iconc=1,7)/    &
1.0006, 1.0018, 1.0043, 1.0100, 1.0232, 1.0539, 1.1291, 1.3342, 2.0863,  &
1.0003, 1.0008, 1.0023, 1.0059, 1.0141, 1.0329, 1.0769, 1.1872, 1.5112,  &
1.0001, 1.0003, 1.0009, 1.0030, 1.0082, 1.0201, 1.0472, 1.1117, 1.2795,  &
1.0001, 1.0001, 1.0003, 1.0009, 1.0038, 1.0113, 1.0287, 1.0685, 1.1653,  &
1.0000, 1.0001, 1.0001, 1.0003, 1.0009, 1.0047, 1.0156, 1.0408, 1.1001,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0003, 1.0009, 1.0058, 1.0211, 1.0576,  &
0.9997, 0.9997, 0.9997, 0.9998, 1.0000, 1.0003, 1.0011, 1.0067, 1.0276 /
data ((supersat( 8,iw,iconc,5, 1),iw=1,9),iconc=1,7)/    &
1.0005, 1.0014, 1.0034, 1.0079, 1.0183, 1.0426, 1.1017, 1.2553, 1.7345,  &
1.0002, 1.0006, 1.0017, 1.0045, 1.0109, 1.0255, 1.0598, 1.1439, 1.3726,  &
1.0001, 1.0002, 1.0006, 1.0022, 1.0061, 1.0152, 1.0360, 1.0850, 1.2079,  &
1.0001, 1.0001, 1.0002, 1.0007, 1.0027, 1.0083, 1.0213, 1.0512, 1.1225,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0007, 1.0032, 1.0112, 1.0299, 1.0733,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0003, 1.0007, 1.0037, 1.0148, 1.0416,  &
0.9998, 0.9998, 0.9998, 0.9999, 1.0000, 1.0002, 1.0008, 1.0040, 1.0192 /
data ((supersat( 8,iw,iconc,6, 1),iw=1,9),iconc=1,7)/    &
1.0005, 1.0011, 1.0027, 1.0064, 1.0150, 1.0352, 1.0838, 1.2073, 1.5612,  &
1.0002, 1.0005, 1.0014, 1.0036, 1.0087, 1.0206, 1.0487, 1.1167, 1.2945,  &
1.0001, 1.0002, 1.0005, 1.0017, 1.0047, 1.0120, 1.0287, 1.0680, 1.1648,  &
1.0001, 1.0001, 1.0002, 1.0005, 1.0020, 1.0063, 1.0165, 1.0401, 1.0960,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0005, 1.0023, 1.0083, 1.0228, 1.0564,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0002, 1.0005, 1.0024, 1.0108, 1.0312,  &
0.9998, 0.9998, 0.9999, 0.9999, 1.0000, 1.0002, 1.0006, 1.0025, 1.0136 /
data ((supersat( 8,iw,iconc,7, 1),iw=1,9),iconc=1,7)/    &
1.0003, 1.0009, 1.0022, 1.0054, 1.0127, 1.0299, 1.0714, 1.1751, 1.4597,  &
1.0002, 1.0004, 1.0011, 1.0029, 1.0072, 1.0173, 1.0410, 1.0983, 1.2448,  &
1.0001, 1.0002, 1.0004, 1.0014, 1.0039, 1.0098, 1.0237, 1.0565, 1.1368,  &
1.0000, 1.0001, 1.0002, 1.0004, 1.0016, 1.0050, 1.0133, 1.0327, 1.0787,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0004, 1.0017, 1.0064, 1.0181, 1.0453,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0002, 1.0004, 1.0017, 1.0081, 1.0244,  &
0.9998, 0.9999, 0.9999, 0.9999, 1.0000, 1.0001, 1.0005, 1.0017, 1.0099 /
data ((supersat( 8,iw,iconc,1, 2),iw=1,9),iconc=1,7)/    &
1.0012, 1.0039, 1.0108, 1.0261, 1.0614, 1.1475, 1.3951, 2.6630, 4.9154,  &
1.0003, 1.0012, 1.0052, 1.0152, 1.0378, 1.0907, 1.2258, 1.6754, 3.8772,  &
1.0002, 1.0003, 1.0013, 1.0068, 1.0212, 1.0548, 1.1356, 1.3600, 2.4072,  &
1.0000, 1.0002, 1.0004, 1.0013, 1.0086, 1.0291, 1.0785, 1.2053, 1.6174,  &
0.9996, 0.9997, 0.9999, 1.0005, 1.0017, 1.0104, 1.0386, 1.1107, 1.3171,  &
0.9984, 0.9985, 0.9986, 0.9991, 1.0000, 1.0024, 1.0120, 1.0487, 1.1525,  &
0.9961, 0.9961, 0.9962, 0.9965, 0.9971, 0.9987, 1.0028, 1.0132, 1.0583 /
data ((supersat( 8,iw,iconc,2, 2),iw=1,9),iconc=1,7)/    &
1.0008, 1.0028, 1.0075, 1.0180, 1.0421, 1.0994, 1.2502, 1.7899, 3.6747,  &
1.0003, 1.0009, 1.0037, 1.0104, 1.0259, 1.0615, 1.1478, 1.3930, 2.5389,  &
1.0001, 1.0003, 1.0010, 1.0047, 1.0145, 1.0373, 1.0904, 1.2256, 1.6681,  &
1.0001, 1.0001, 1.0003, 1.0010, 1.0059, 1.0201, 1.0535, 1.1343, 1.3580,  &
0.9999, 0.9999, 1.0001, 1.0004, 1.0010, 1.0072, 1.0271, 1.0756, 1.2011,  &
0.9993, 0.9993, 0.9994, 0.9997, 1.0002, 1.0016, 1.0083, 1.0349, 1.1044,  &
0.9978, 0.9978, 0.9979, 0.9981, 0.9985, 0.9995, 1.0019, 1.0092, 1.0428 /
data ((supersat( 8,iw,iconc,3, 2),iw=1,9),iconc=1,7)/    &
1.0006, 1.0021, 1.0053, 1.0129, 1.0303, 1.0709, 1.1730, 1.4779, 2.7177,  &
1.0002, 1.0007, 1.0026, 1.0074, 1.0184, 1.0434, 1.1030, 1.2588, 1.7975,  &
1.0001, 1.0002, 1.0007, 1.0033, 1.0102, 1.0262, 1.0630, 1.1525, 1.4045,  &
1.0001, 1.0001, 1.0002, 1.0007, 1.0040, 1.0140, 1.0372, 1.0919, 1.2313,  &
1.0000, 1.0000, 1.0001, 1.0003, 1.0008, 1.0047, 1.0188, 1.0524, 1.1347,  &
0.9997, 0.9997, 0.9998, 0.9999, 1.0003, 1.0010, 1.0054, 1.0244, 1.0723,  &
0.9988, 0.9988, 0.9988, 0.9989, 0.9992, 0.9998, 1.0013, 1.0058, 1.0303 /
data ((supersat( 8,iw,iconc,4, 2),iw=1,9),iconc=1,7)/    &
1.0005, 1.0015, 1.0040, 1.0097, 1.0228, 1.0534, 1.1285, 1.3334, 2.0845,  &
1.0002, 1.0005, 1.0019, 1.0054, 1.0135, 1.0322, 1.0760, 1.1859, 1.5091,  &
1.0001, 1.0002, 1.0005, 1.0023, 1.0073, 1.0190, 1.0458, 1.1098, 1.2766,  &
1.0001, 1.0001, 1.0002, 1.0005, 1.0027, 1.0098, 1.0267, 1.0657, 1.1612,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0006, 1.0031, 1.0130, 1.0372, 1.0945,  &
0.9998, 0.9998, 0.9999, 1.0000, 1.0002, 1.0007, 1.0033, 1.0169, 1.0509,  &
0.9992, 0.9993, 0.9993, 0.9994, 0.9995, 0.9999, 1.0008, 1.0034, 1.0209 /
data ((supersat( 8,iw,iconc,5, 2),iw=1,9),iconc=1,7)/    &
1.0004, 1.0012, 1.0031, 1.0076, 1.0180, 1.0423, 1.1013, 1.2547, 1.7337,  &
1.0001, 1.0004, 1.0015, 1.0041, 1.0104, 1.0250, 1.0592, 1.1431, 1.3714,  &
1.0001, 1.0001, 1.0004, 1.0017, 1.0054, 1.0143, 1.0350, 1.0837, 1.2061,  &
1.0000, 1.0001, 1.0002, 1.0004, 1.0019, 1.0071, 1.0198, 1.0493, 1.1198,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0004, 1.0020, 1.0092, 1.0272, 1.0695,  &
0.9999, 0.9999, 0.9999, 1.0000, 1.0002, 1.0005, 1.0020, 1.0117, 1.0368,  &
0.9995, 0.9995, 0.9995, 0.9995, 0.9997, 0.9999, 1.0005, 1.0022, 1.0142 /
data ((supersat( 8,iw,iconc,6, 2),iw=1,9),iconc=1,7)/    &
1.0003, 1.0009, 1.0024, 1.0062, 1.0147, 1.0349, 1.0835, 1.2069, 1.5607,  &
1.0001, 1.0003, 1.0011, 1.0032, 1.0083, 1.0202, 1.0482, 1.1161, 1.2937,  &
1.0001, 1.0001, 1.0003, 1.0013, 1.0042, 1.0113, 1.0279, 1.0670, 1.1636,  &
1.0000, 1.0001, 1.0001, 1.0003, 1.0014, 1.0054, 1.0153, 1.0386, 1.0941,  &
1.0000, 1.0000, 1.0001, 1.0001, 1.0003, 1.0014, 1.0068, 1.0207, 1.0536,  &
0.9999, 0.9999, 0.9999, 1.0000, 1.0001, 1.0004, 1.0014, 1.0083, 1.0275,  &
0.9995, 0.9996, 0.9996, 0.9996, 0.9997, 0.9999, 1.0004, 1.0015, 1.0098 /
data ((supersat( 8,iw,iconc,7, 2),iw=1,9),iconc=1,7)/    &
1.0003, 1.0008, 1.0021, 1.0052, 1.0125, 1.0297, 1.0711, 1.1749, 1.4593,  &
1.0001, 1.0003, 1.0009, 1.0027, 1.0069, 1.0169, 1.0406, 1.0978, 1.2442,  &
1.0001, 1.0001, 1.0003, 1.0010, 1.0034, 1.0093, 1.0231, 1.0558, 1.1358,  &
1.0000, 1.0001, 1.0001, 1.0003, 1.0011, 1.0042, 1.0123, 1.0316, 1.0772,  &
1.0000, 1.0000, 1.0001, 1.0001, 1.0003, 1.0010, 1.0051, 1.0164, 1.0431,  &
0.9999, 0.9999, 0.9999, 1.0000, 1.0001, 1.0003, 1.0011, 1.0061, 1.0214,  &
0.9996, 0.9996, 0.9996, 0.9996, 0.9997, 0.9999, 1.0003, 1.0011, 1.0069 /
data ((supersat( 8,iw,iconc,1, 3),iw=1,9),iconc=1,7)/    &
1.0007, 1.0031, 1.0097, 1.0247, 1.0596, 1.1449, 1.3906, 2.6497, 4.8930,  &
1.0002, 1.0007, 1.0039, 1.0132, 1.0353, 1.0870, 1.2198, 1.6633, 3.8339,  &
1.0001, 1.0003, 1.0007, 1.0046, 1.0180, 1.0499, 1.1278, 1.3457, 2.3670,  &
0.9998, 0.9999, 1.0003, 1.0010, 1.0053, 1.0236, 1.0692, 1.1885, 1.5801,  &
0.9990, 0.9991, 0.9993, 0.9999, 1.0013, 1.0058, 1.0296, 1.0933, 1.2803,  &
0.9971, 0.9971, 0.9973, 0.9977, 0.9987, 1.0011, 1.0075, 1.0351, 1.1213,  &
0.9943, 0.9943, 0.9944, 0.9946, 0.9952, 0.9966, 1.0003, 1.0105, 1.0396 /
data ((supersat( 8,iw,iconc,2, 3),iw=1,9),iconc=1,7)/    &
1.0005, 1.0022, 1.0068, 1.0171, 1.0410, 1.0980, 1.2480, 1.7859, 3.6621,  &
1.0002, 1.0005, 1.0027, 1.0092, 1.0243, 1.0593, 1.1447, 1.3876, 2.5235,  &
1.0001, 1.0002, 1.0006, 1.0032, 1.0125, 1.0344, 1.0861, 1.2186, 1.6537,  &
1.0000, 1.0001, 1.0002, 1.0006, 1.0036, 1.0165, 1.0479, 1.1253, 1.3412,  &
0.9996, 0.9996, 0.9998, 1.0001, 1.0009, 1.0039, 1.0211, 1.0653, 1.1821,  &
0.9985, 0.9985, 0.9986, 0.9989, 0.9995, 1.0009, 1.0042, 1.0257, 1.0859,  &
0.9963, 0.9963, 0.9963, 0.9965, 0.9969, 0.9979, 1.0004, 1.0067, 1.0296 /
data ((supersat( 8,iw,iconc,3, 3),iw=1,9),iconc=1,7)/    &
1.0004, 1.0016, 1.0048, 1.0123, 1.0295, 1.0700, 1.1718, 1.4759, 2.7119,  &
1.0001, 1.0004, 1.0019, 1.0065, 1.0173, 1.0421, 1.1012, 1.2560, 1.7921,  &
1.0001, 1.0001, 1.0004, 1.0022, 1.0087, 1.0242, 1.0603, 1.1485, 1.3976,  &
1.0000, 1.0001, 1.0002, 1.0005, 1.0024, 1.0115, 1.0337, 1.0866, 1.2224,  &
0.9998, 0.9999, 0.9999, 1.0002, 1.0006, 1.0025, 1.0147, 1.0459, 1.1238,  &
0.9992, 0.9992, 0.9993, 0.9995, 0.9998, 1.0007, 1.0027, 1.0181, 1.0606,  &
0.9977, 0.9977, 0.9977, 0.9978, 0.9981, 0.9988, 1.0003, 1.0041, 1.0210 /
data ((supersat( 8,iw,iconc,4, 3),iw=1,9),iconc=1,7)/    &
1.0003, 1.0012, 1.0036, 1.0092, 1.0223, 1.0528, 1.1278, 1.3323, 2.0821,  &
1.0001, 1.0003, 1.0014, 1.0048, 1.0128, 1.0312, 1.0748, 1.1843, 1.5063,  &
1.0001, 1.0001, 1.0003, 1.0016, 1.0062, 1.0176, 1.0441, 1.1073, 1.2728,  &
1.0000, 1.0001, 1.0001, 1.0003, 1.0016, 1.0080, 1.0242, 1.0623, 1.1558,  &
0.9999, 0.9999, 1.0000, 1.0001, 1.0004, 1.0017, 1.0101, 1.0327, 1.0876,  &
0.9995, 0.9995, 0.9996, 0.9997, 0.9999, 1.0005, 1.0018, 1.0123, 1.0430,  &
0.9985, 0.9985, 0.9985, 0.9986, 0.9987, 0.9992, 1.0002, 1.0025, 1.0143 /
data ((supersat( 8,iw,iconc,5, 3),iw=1,9),iconc=1,7)/    &
1.0003, 1.0009, 1.0028, 1.0072, 1.0176, 1.0419, 1.1008, 1.2540, 1.7325,  &
1.0001, 1.0003, 1.0011, 1.0036, 1.0098, 1.0243, 1.0583, 1.1420, 1.3697,  &
1.0001, 1.0001, 1.0003, 1.0011, 1.0046, 1.0133, 1.0337, 1.0820, 1.2036,  &
1.0000, 1.0001, 1.0001, 1.0003, 1.0012, 1.0058, 1.0180, 1.0468, 1.1163,  &
0.9999, 1.0000, 1.0000, 1.0001, 1.0003, 1.0012, 1.0070, 1.0240, 1.0647,  &
0.9997, 0.9997, 0.9997, 0.9998, 1.0000, 1.0003, 1.0012, 1.0083, 1.0312,  &
0.9988, 0.9988, 0.9989, 0.9989, 0.9990, 0.9993, 1.0000, 1.0016, 1.0093 /
data ((supersat( 8,iw,iconc,6, 3),iw=1,9),iconc=1,7)/    &
1.0002, 1.0008, 1.0022, 1.0059, 1.0144, 1.0346, 1.0831, 1.2064, 1.5599,  &
1.0001, 1.0002, 1.0008, 1.0029, 1.0079, 1.0197, 1.0475, 1.1153, 1.2926,  &
1.0000, 1.0001, 1.0002, 1.0009, 1.0035, 1.0105, 1.0269, 1.0658, 1.1619,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0009, 1.0043, 1.0139, 1.0368, 1.0916,  &
0.9999, 1.0000, 1.0000, 1.0001, 1.0002, 1.0009, 1.0050, 1.0182, 1.0501,  &
0.9997, 0.9997, 0.9998, 0.9998, 1.0000, 1.0002, 1.0009, 1.0056, 1.0233,  &
0.9990, 0.9990, 0.9990, 0.9991, 0.9992, 0.9994, 0.9999, 1.0011, 1.0060 /
data ((supersat( 8,iw,iconc,7, 3),iw=1,9),iconc=1,7)/    &
1.0002, 1.0006, 1.0019, 1.0050, 1.0122, 1.0294, 1.0708, 1.1745, 1.4588,  &
1.0001, 1.0002, 1.0007, 1.0023, 1.0065, 1.0165, 1.0401, 1.0972, 1.2433,  &
1.0000, 1.0001, 1.0002, 1.0007, 1.0028, 1.0085, 1.0223, 1.0548, 1.1345,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0007, 1.0033, 1.0111, 1.0301, 1.0752,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0002, 1.0007, 1.0037, 1.0143, 1.0403,  &
0.9997, 0.9997, 0.9998, 0.9998, 0.9999, 1.0002, 1.0007, 1.0039, 1.0180,  &
0.9991, 0.9991, 0.9991, 0.9991, 0.9992, 0.9994, 0.9999, 1.0008, 1.0041 /
data ((supersat( 8,iw,iconc,1, 4),iw=1,9),iconc=1,7)/    &
1.0004, 1.0022, 1.0083, 1.0230, 1.0572, 1.1414, 1.3847, 2.6321, 4.8631,  &
1.0002, 1.0004, 1.0024, 1.0110, 1.0321, 1.0823, 1.2121, 1.6476, 3.7770,  &
1.0000, 1.0002, 1.0006, 1.0025, 1.0142, 1.0440, 1.1181, 1.3276, 2.3151,  &
0.9994, 0.9995, 0.9999, 1.0007, 1.0027, 1.0177, 1.0586, 1.1685, 1.5343,  &
0.9979, 0.9980, 0.9983, 0.9989, 1.0004, 1.0040, 1.0208, 1.0749, 1.2388,  &
0.9954, 0.9954, 0.9956, 0.9960, 0.9969, 0.9993, 1.0053, 1.0232, 1.0914,  &
0.9927, 0.9927, 0.9928, 0.9929, 0.9934, 0.9946, 0.9976, 1.0061, 1.0328 /
data ((supersat( 8,iw,iconc,2, 4),iw=1,9),iconc=1,7)/    &
1.0003, 1.0016, 1.0059, 1.0160, 1.0397, 1.0961, 1.2452, 1.7805, 3.6454,  &
1.0001, 1.0003, 1.0018, 1.0077, 1.0223, 1.0566, 1.1406, 1.3806, 2.5030,  &
1.0000, 1.0002, 1.0004, 1.0018, 1.0100, 1.0308, 1.0807, 1.2095, 1.6350,  &
0.9998, 0.9999, 1.0001, 1.0005, 1.0019, 1.0125, 1.0415, 1.1143, 1.3202,  &
0.9990, 0.9991, 0.9992, 0.9996, 1.0005, 1.0025, 1.0150, 1.0539, 1.1599,  &
0.9972, 0.9973, 0.9974, 0.9976, 0.9983, 0.9998, 1.0035, 1.0171, 1.0669,  &
0.9945, 0.9945, 0.9946, 0.9947, 0.9951, 0.9960, 0.9983, 1.0042, 1.0187 /
data ((supersat( 8,iw,iconc,3, 4),iw=1,9),iconc=1,7)/    &
1.0003, 1.0012, 1.0042, 1.0116, 1.0286, 1.0689, 1.1703, 1.4733, 2.7040,  &
1.0001, 1.0003, 1.0013, 1.0055, 1.0160, 1.0404, 1.0988, 1.2523, 1.7849,  &
1.0000, 1.0001, 1.0003, 1.0013, 1.0070, 1.0219, 1.0570, 1.1434, 1.3886,  &
0.9999, 1.0000, 1.0001, 1.0004, 1.0013, 1.0087, 1.0295, 1.0799, 1.2111,  &
0.9995, 0.9995, 0.9996, 0.9999, 1.0004, 1.0016, 1.0104, 1.0385, 1.1107,  &
0.9984, 0.9984, 0.9985, 0.9986, 0.9991, 1.0000, 1.0022, 1.0119, 1.0483,  &
0.9961, 0.9961, 0.9962, 0.9963, 0.9965, 0.9972, 0.9988, 1.0026, 1.0130 /
data ((supersat( 8,iw,iconc,4, 4),iw=1,9),iconc=1,7)/    &
1.0002, 1.0009, 1.0031, 1.0087, 1.0216, 1.0521, 1.1268, 1.3308, 2.0788,  &
1.0001, 1.0002, 1.0009, 1.0040, 1.0118, 1.0301, 1.0733, 1.1821, 1.5026,  &
1.0000, 1.0001, 1.0002, 1.0010, 1.0050, 1.0159, 1.0418, 1.1041, 1.2676,  &
1.0000, 1.0000, 1.0001, 1.0003, 1.0010, 1.0060, 1.0213, 1.0579, 1.1490,  &
0.9997, 0.9997, 0.9998, 1.0000, 1.0003, 1.0011, 1.0070, 1.0276, 1.0792,  &
0.9990, 0.9990, 0.9990, 0.9991, 0.9994, 1.0000, 1.0014, 1.0078, 1.0346,  &
0.9972, 0.9972, 0.9973, 0.9973, 0.9975, 0.9980, 0.9990, 1.0015, 1.0084 /
data ((supersat( 8,iw,iconc,5, 4),iw=1,9),iconc=1,7)/    &
1.0002, 1.0007, 1.0024, 1.0068, 1.0171, 1.0413, 1.1001, 1.2531, 1.7308,  &
1.0001, 1.0002, 1.0007, 1.0030, 1.0091, 1.0234, 1.0572, 1.1405, 1.3675,  &
1.0000, 1.0001, 1.0002, 1.0007, 1.0036, 1.0120, 1.0320, 1.0798, 1.2004,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0007, 1.0042, 1.0158, 1.0437, 1.1117,  &
0.9998, 0.9998, 0.9999, 1.0000, 1.0002, 1.0008, 1.0047, 1.0202, 1.0589,  &
0.9992, 0.9992, 0.9993, 0.9994, 0.9996, 1.0000, 1.0009, 1.0050, 1.0251,  &
0.9978, 0.9978, 0.9978, 0.9979, 0.9980, 0.9983, 0.9991, 1.0008, 1.0052 /
data ((supersat( 8,iw,iconc,6, 4),iw=1,9),iconc=1,7)/    &
1.0002, 1.0005, 1.0020, 1.0055, 1.0140, 1.0341, 1.0826, 1.2058, 1.5588,  &
1.0001, 1.0002, 1.0006, 1.0024, 1.0073, 1.0190, 1.0467, 1.1143, 1.2910,  &
1.0000, 1.0001, 1.0002, 1.0006, 1.0028, 1.0095, 1.0256, 1.0641, 1.1596,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0006, 1.0031, 1.0121, 1.0344, 1.0882,  &
0.9998, 0.9998, 0.9999, 1.0000, 1.0002, 1.0006, 1.0032, 1.0153, 1.0457,  &
0.9994, 0.9994, 0.9994, 0.9995, 0.9996, 1.0000, 1.0007, 1.0033, 1.0186,  &
0.9981, 0.9981, 0.9981, 0.9982, 0.9983, 0.9985, 0.9991, 1.0004, 1.0036 /
data ((supersat( 8,iw,iconc,7, 4),iw=1,9),iconc=1,7)/    &
1.0001, 1.0004, 1.0016, 1.0046, 1.0119, 1.0291, 1.0704, 1.1739, 1.4580,  &
1.0001, 1.0001, 1.0004, 1.0019, 1.0060, 1.0159, 1.0394, 1.0964, 1.2422,  &
1.0000, 1.0001, 1.0001, 1.0004, 1.0022, 1.0077, 1.0212, 1.0535, 1.1327,  &
1.0000, 1.0000, 1.0001, 1.0001, 1.0005, 1.0023, 1.0097, 1.0281, 1.0726,  &
0.9998, 0.9999, 0.9999, 1.0000, 1.0001, 1.0005, 1.0024, 1.0119, 1.0369,  &
0.9994, 0.9994, 0.9994, 0.9995, 0.9996, 0.9999, 1.0005, 1.0024, 1.0143,  &
0.9982, 0.9983, 0.9983, 0.9983, 0.9984, 0.9986, 0.9991, 1.0002, 1.0026 /
data ((supersat( 8,iw,iconc,1, 5),iw=1,9),iconc=1,7)/    &
1.0003, 1.0017, 1.0074, 1.0218, 1.0556, 1.1390, 1.3805, 2.6191, 4.8415,  &
1.0001, 1.0004, 1.0017, 1.0096, 1.0300, 1.0790, 1.2066, 1.6363, 3.7360,  &
0.9998, 1.0000, 1.0005, 1.0018, 1.0119, 1.0402, 1.1116, 1.3150, 2.2781,  &
0.9990, 0.9992, 0.9995, 1.0004, 1.0025, 1.0142, 1.0520, 1.1555, 1.5034,  &
0.9971, 0.9972, 0.9975, 0.9981, 0.9996, 1.0033, 1.0161, 1.0642, 1.2132,  &
0.9943, 0.9944, 0.9945, 0.9949, 0.9958, 0.9980, 1.0037, 1.0175, 1.0756,  &
0.9919, 0.9920, 0.9920, 0.9922, 0.9926, 0.9936, 0.9962, 1.0036, 1.0268 /
data ((supersat( 8,iw,iconc,2, 5),iw=1,9),iconc=1,7)/    &
1.0003, 1.0013, 1.0052, 1.0152, 1.0387, 1.0948, 1.2431, 1.7765, 3.6330,  &
1.0001, 1.0003, 1.0013, 1.0067, 1.0210, 1.0547, 1.1378, 1.3756, 2.4882,  &
1.0000, 1.0001, 1.0004, 1.0013, 1.0084, 1.0284, 1.0769, 1.2032, 1.6217,  &
0.9996, 0.9997, 0.9999, 1.0004, 1.0016, 1.0101, 1.0373, 1.1069, 1.3057,  &
0.9985, 0.9986, 0.9987, 0.9991, 1.0001, 1.0022, 1.0116, 1.0470, 1.1457,  &
0.9963, 0.9963, 0.9964, 0.9967, 0.9973, 0.9989, 1.0026, 1.0128, 1.0564,  &
0.9935, 0.9935, 0.9936, 0.9937, 0.9940, 0.9948, 0.9969, 1.0024, 1.0185 /
data ((supersat( 8,iw,iconc,3, 5),iw=1,9),iconc=1,7)/    &
1.0002, 1.0009, 1.0038, 1.0110, 1.0280, 1.0681, 1.1691, 1.4714, 2.6982,  &
1.0001, 1.0002, 1.0009, 1.0048, 1.0150, 1.0392, 1.0971, 1.2496, 1.7796,  &
1.0000, 1.0001, 1.0003, 1.0010, 1.0059, 1.0203, 1.0546, 1.1398, 1.3821,  &
0.9998, 0.9999, 1.0000, 1.0003, 1.0010, 1.0070, 1.0267, 1.0755, 1.2032,  &
0.9992, 0.9992, 0.9993, 0.9996, 1.0002, 1.0014, 1.0079, 1.0339, 1.1021,  &
0.9977, 0.9977, 0.9978, 0.9980, 0.9984, 0.9994, 1.0017, 1.0087, 1.0411,  &
0.9951, 0.9951, 0.9952, 0.9952, 0.9955, 0.9961, 0.9977, 1.0014, 1.0118 /
data ((supersat( 8,iw,iconc,4, 5),iw=1,9),iconc=1,7)/    &
1.0002, 1.0007, 1.0029, 1.0083, 1.0212, 1.0515, 1.1261, 1.3297, 2.0764,  &
1.0001, 1.0002, 1.0007, 1.0035, 1.0111, 1.0292, 1.0722, 1.1805, 1.4998,  &
1.0000, 1.0001, 1.0002, 1.0007, 1.0042, 1.0148, 1.0403, 1.1018, 1.2639,  &
0.9999, 0.9999, 1.0000, 1.0002, 1.0008, 1.0048, 1.0193, 1.0550, 1.1442,  &
0.9995, 0.9995, 0.9996, 0.9998, 1.0001, 1.0009, 1.0052, 1.0244, 1.0736,  &
0.9985, 0.9985, 0.9985, 0.9987, 0.9989, 0.9996, 1.0011, 1.0055, 1.0296,  &
0.9963, 0.9963, 0.9963, 0.9964, 0.9966, 0.9970, 0.9981, 1.0006, 1.0072 /
data ((supersat( 8,iw,iconc,5, 5),iw=1,9),iconc=1,7)/    &
1.0001, 1.0005, 1.0022, 1.0065, 1.0167, 1.0409, 1.0996, 1.2524, 1.7296,  &
1.0001, 1.0001, 1.0005, 1.0026, 1.0086, 1.0228, 1.0565, 1.1395, 1.3658,  &
1.0000, 1.0001, 1.0002, 1.0005, 1.0030, 1.0112, 1.0309, 1.0782, 1.1980,  &
0.9999, 1.0000, 1.0000, 1.0002, 1.0006, 1.0033, 1.0143, 1.0416, 1.1085,  &
0.9996, 0.9997, 0.9997, 0.9998, 1.0001, 1.0007, 1.0034, 1.0179, 1.0550,  &
0.9988, 0.9989, 0.9989, 0.9990, 0.9992, 0.9997, 1.0007, 1.0036, 1.0214,  &
0.9970, 0.9970, 0.9970, 0.9971, 0.9972, 0.9975, 0.9983, 1.0001, 1.0045 /
data ((supersat( 8,iw,iconc,6, 5),iw=1,9),iconc=1,7)/    &
1.0001, 1.0004, 1.0018, 1.0053, 1.0138, 1.0338, 1.0822, 1.2053, 1.5580,  &
1.0001, 1.0001, 1.0004, 1.0021, 1.0068, 1.0185, 1.0461, 1.1135, 1.2898,  &
1.0000, 1.0001, 1.0001, 1.0004, 1.0023, 1.0087, 1.0247, 1.0629, 1.1579,  &
0.9999, 1.0000, 1.0000, 1.0001, 1.0004, 1.0024, 1.0110, 1.0328, 1.0859,  &
0.9997, 0.9997, 0.9998, 0.9999, 1.0001, 1.0005, 1.0024, 1.0134, 1.0428,  &
0.9990, 0.9990, 0.9991, 0.9991, 0.9993, 0.9997, 1.0005, 1.0025, 1.0158,  &
0.9974, 0.9974, 0.9974, 0.9974, 0.9975, 0.9978, 0.9984, 0.9998, 1.0030 /
data ((supersat( 8,iw,iconc,7, 5),iw=1,9),iconc=1,7)/    &
1.0001, 1.0004, 1.0014, 1.0044, 1.0116, 1.0288, 1.0701, 1.1736, 1.4574,  &
1.0000, 1.0001, 1.0004, 1.0016, 1.0056, 1.0155, 1.0389, 1.0957, 1.2413,  &
1.0000, 1.0000, 1.0001, 1.0004, 1.0018, 1.0071, 1.0205, 1.0525, 1.1314,  &
0.9999, 1.0000, 1.0000, 1.0001, 1.0004, 1.0018, 1.0087, 1.0268, 1.0708,  &
0.9997, 0.9997, 0.9998, 0.9999, 1.0000, 1.0004, 1.0018, 1.0104, 1.0345,  &
0.9991, 0.9991, 0.9991, 0.9992, 0.9993, 0.9996, 1.0003, 1.0020, 1.0120,  &
0.9976, 0.9976, 0.9976, 0.9976, 0.9977, 0.9979, 0.9984, 0.9996, 1.0021 /
data ((supersat( 8,iw,iconc,1, 6),iw=1,9),iconc=1,7)/    &
1.0003, 1.0013, 1.0068, 1.0209, 1.0543, 1.1370, 1.3771, 2.6086, 4.8234,  &
1.0001, 1.0004, 1.0013, 1.0085, 1.0283, 1.0765, 1.2023, 1.6272, 3.7027,  &
0.9997, 0.9999, 1.0004, 1.0014, 1.0103, 1.0373, 1.1065, 1.3051, 2.2485,  &
0.9986, 0.9988, 0.9992, 1.0001, 1.0022, 1.0118, 1.0472, 1.1457, 1.4798,  &
0.9965, 0.9966, 0.9968, 0.9974, 0.9990, 1.0027, 1.0131, 1.0569, 1.1950,  &
0.9936, 0.9937, 0.9938, 0.9941, 0.9950, 0.9971, 1.0025, 1.0185, 1.0654,  &
0.9915, 0.9915, 0.9916, 0.9917, 0.9921, 0.9930, 0.9953, 1.0019, 1.0227 /
data ((supersat( 8,iw,iconc,2, 6),iw=1,9),iconc=1,7)/    &
1.0002, 1.0010, 1.0048, 1.0147, 1.0379, 1.0937, 1.2414, 1.7733, 3.6230,  &
1.0001, 1.0002, 1.0010, 1.0060, 1.0199, 1.0532, 1.1354, 1.3715, 2.4760,  &
0.9999, 1.0000, 1.0003, 1.0011, 1.0072, 1.0266, 1.0740, 1.1981, 1.6110,  &
0.9994, 0.9995, 0.9997, 1.0002, 1.0015, 1.0084, 1.0342, 1.1013, 1.2945,  &
0.9981, 0.9981, 0.9983, 0.9987, 0.9997, 1.0019, 1.0093, 1.0421, 1.1353,  &
0.9956, 0.9956, 0.9957, 0.9960, 0.9966, 0.9981, 1.0018, 1.0120, 1.0493,  &
0.9929, 0.9929, 0.9929, 0.9930, 0.9933, 0.9941, 0.9960, 1.0011, 1.0160 /
data ((supersat( 8,iw,iconc,3, 6),iw=1,9),iconc=1,7)/    &
1.0002, 1.0008, 1.0035, 1.0106, 1.0275, 1.0674, 1.1682, 1.4698, 2.6933,  &
1.0001, 1.0002, 1.0008, 1.0043, 1.0143, 1.0382, 1.0957, 1.2475, 1.7753,  &
1.0000, 1.0001, 1.0002, 1.0008, 1.0051, 1.0190, 1.0528, 1.1369, 1.3769,  &
0.9997, 0.9998, 0.9999, 1.0002, 1.0010, 1.0058, 1.0246, 1.0720, 1.1970,  &
0.9989, 0.9990, 0.9991, 0.9993, 0.9999, 1.0013, 1.0063, 1.0306, 1.0957,  &
0.9971, 0.9972, 0.9972, 0.9974, 0.9978, 0.9988, 1.0012, 1.0075, 1.0362,  &
0.9944, 0.9944, 0.9944, 0.9945, 0.9947, 0.9953, 0.9968, 1.0004, 1.0104 /
data ((supersat( 8,iw,iconc,4, 6),iw=1,9),iconc=1,7)/    &
1.0001, 1.0006, 1.0026, 1.0080, 1.0208, 1.0511, 1.1255, 1.3288, 2.0744,  &
1.0001, 1.0002, 1.0006, 1.0031, 1.0106, 1.0286, 1.0713, 1.1792, 1.4975,  &
1.0000, 1.0001, 1.0002, 1.0006, 1.0036, 1.0139, 1.0390, 1.1000, 1.2609,  &
0.9998, 0.9999, 1.0000, 1.0002, 1.0006, 1.0039, 1.0178, 1.0527, 1.1404,  &
0.9993, 0.9994, 0.9994, 0.9996, 1.0000, 1.0008, 1.0041, 1.0221, 1.0694,  &
0.9980, 0.9980, 0.9981, 0.9982, 0.9985, 0.9992, 1.0008, 1.0046, 1.0261,  &
0.9956, 0.9956, 0.9956, 0.9957, 0.9959, 0.9963, 0.9974, 0.9999, 1.0063 /
data ((supersat( 8,iw,iconc,5, 6),iw=1,9),iconc=1,7)/    &
1.0001, 1.0005, 1.0020, 1.0062, 1.0165, 1.0406, 1.0992, 1.2518, 1.7286,  &
1.0001, 1.0001, 1.0005, 1.0023, 1.0082, 1.0223, 1.0558, 1.1386, 1.3644,  &
1.0000, 1.0001, 1.0001, 1.0005, 1.0026, 1.0105, 1.0300, 1.0769, 1.1960,  &
0.9999, 0.9999, 1.0000, 1.0001, 1.0005, 1.0027, 1.0132, 1.0400, 1.1059,  &
0.9995, 0.9995, 0.9996, 0.9997, 1.0000, 1.0006, 1.0027, 1.0161, 1.0521,  &
0.9985, 0.9985, 0.9985, 0.9986, 0.9988, 0.9993, 1.0004, 1.0030, 1.0188,  &
0.9963, 0.9963, 0.9964, 0.9964, 0.9965, 0.9969, 0.9977, 0.9995, 1.0039 /
data ((supersat( 8,iw,iconc,6, 6),iw=1,9),iconc=1,7)/    &
1.0001, 1.0004, 1.0016, 1.0051, 1.0135, 1.0336, 1.0819, 1.2048, 1.5573,  &
1.0000, 1.0001, 1.0004, 1.0018, 1.0065, 1.0181, 1.0456, 1.1128, 1.2889,  &
1.0000, 1.0000, 1.0001, 1.0004, 1.0019, 1.0082, 1.0240, 1.0620, 1.1565,  &
0.9999, 0.9999, 1.0000, 1.0001, 1.0004, 1.0019, 1.0101, 1.0315, 1.0840,  &
0.9996, 0.9996, 0.9997, 0.9998, 1.0000, 1.0004, 1.0020, 1.0120, 1.0406,  &
0.9987, 0.9987, 0.9987, 0.9988, 0.9990, 0.9994, 1.0002, 1.0022, 1.0138,  &
0.9967, 0.9967, 0.9968, 0.9968, 0.9969, 0.9972, 0.9978, 0.9992, 1.0025 /
data ((supersat( 8,iw,iconc,7, 6),iw=1,9),iconc=1,7)/    &
1.0001, 1.0003, 1.0013, 1.0043, 1.0114, 1.0286, 1.0698, 1.1732, 1.4569,  &
1.0000, 1.0001, 1.0003, 1.0015, 1.0054, 1.0151, 1.0385, 1.0952, 1.2406,  &
1.0000, 1.0000, 1.0001, 1.0003, 1.0015, 1.0066, 1.0198, 1.0517, 1.1303,  &
0.9999, 0.9999, 1.0000, 1.0001, 1.0003, 1.0015, 1.0079, 1.0258, 1.0693,  &
0.9996, 0.9996, 0.9997, 0.9998, 1.0000, 1.0003, 1.0015, 1.0092, 1.0327,  &
0.9988, 0.9988, 0.9988, 0.9989, 0.9990, 0.9994, 1.0001, 1.0017, 1.0103,  &
0.9969, 0.9970, 0.9970, 0.9970, 0.9971, 0.9973, 0.9979, 0.9991, 1.0017 /
data ((supersat( 8,iw,iconc,1, 7),iw=1,9),iconc=1,7)/    &
1.0002, 1.0011, 1.0062, 1.0201, 1.0532, 1.1354, 1.3741, 2.5994, 4.8081,  &
1.0001, 1.0003, 1.0011, 1.0076, 1.0270, 1.0743, 1.1986, 1.6194, 3.6744,  &
0.9996, 0.9998, 1.0003, 1.0015, 1.0090, 1.0350, 1.1024, 1.2968, 2.2234,  &
0.9983, 0.9985, 0.9989, 0.9998, 1.0020, 1.0101, 1.0435, 1.1379, 1.4604,  &
0.9959, 0.9960, 0.9963, 0.9969, 0.9984, 1.0021, 1.0109, 1.0514, 1.1809,  &
0.9931, 0.9932, 0.9933, 0.9936, 0.9944, 0.9963, 1.0015, 1.0166, 1.0582,  &
0.9912, 0.9913, 0.9913, 0.9914, 0.9917, 0.9925, 0.9947, 1.0007, 1.0196 /
data ((supersat( 8,iw,iconc,2, 7),iw=1,9),iconc=1,7)/    &
1.0002, 1.0009, 1.0045, 1.0142, 1.0373, 1.0927, 1.2400, 1.7704, 3.6143,  &
1.0001, 1.0002, 1.0009, 1.0054, 1.0191, 1.0519, 1.1335, 1.3680, 2.4655,  &
0.9999, 1.0000, 1.0003, 1.0009, 1.0063, 1.0251, 1.0716, 1.1939, 1.6019,  &
0.9992, 0.9993, 0.9996, 1.0001, 1.0014, 1.0071, 1.0318, 1.0967, 1.2852,  &
0.9977, 0.9977, 0.9979, 0.9983, 0.9993, 1.0015, 1.0077, 1.0384, 1.1270,  &
0.9950, 0.9951, 0.9951, 0.9954, 0.9960, 0.9975, 1.0011, 1.0110, 1.0441,  &
0.9924, 0.9924, 0.9925, 0.9926, 0.9928, 0.9936, 0.9953, 1.0001, 1.0140 /
data ((supersat( 8,iw,iconc,3, 7),iw=1,9),iconc=1,7)/    &
1.0002, 1.0007, 1.0032, 1.0103, 1.0271, 1.0669, 1.1673, 1.4684, 2.6892,  &
1.0001, 1.0002, 1.0007, 1.0039, 1.0137, 1.0374, 1.0946, 1.2456, 1.7716,  &
1.0000, 1.0000, 1.0002, 1.0007, 1.0044, 1.0180, 1.0513, 1.1345, 1.3724,  &
0.9996, 0.9997, 0.9998, 1.0002, 1.0009, 1.0048, 1.0230, 1.0692, 1.1919,  &
0.9987, 0.9987, 0.9988, 0.9991, 0.9997, 1.0011, 1.0051, 1.0280, 1.0906,  &
0.9966, 0.9967, 0.9967, 0.9969, 0.9973, 0.9983, 1.0007, 1.0070, 1.0325,  &
0.9938, 0.9938, 0.9939, 0.9939, 0.9942, 0.9947, 0.9961, 0.9996, 1.0092 /
data ((supersat( 8,iw,iconc,4, 7),iw=1,9),iconc=1,7)/    &
1.0001, 1.0005, 1.0024, 1.0078, 1.0205, 1.0507, 1.1249, 1.3280, 2.0726,  &
1.0001, 1.0001, 1.0005, 1.0028, 1.0102, 1.0280, 1.0705, 1.1780, 1.4956,  &
1.0000, 1.0000, 1.0002, 1.0005, 1.0031, 1.0132, 1.0380, 1.0984, 1.2584,  &
0.9998, 0.9998, 0.9999, 1.0001, 1.0006, 1.0032, 1.0166, 1.0508, 1.1372,  &
0.9992, 0.9992, 0.9993, 0.9994, 0.9998, 1.0007, 1.0034, 1.0202, 1.0661,  &
0.9976, 0.9977, 0.9977, 0.9978, 0.9981, 0.9988, 1.0004, 1.0043, 1.0235,  &
0.9950, 0.9950, 0.9951, 0.9951, 0.9953, 0.9957, 0.9968, 0.9993, 1.0056 /
data ((supersat( 8,iw,iconc,5, 7),iw=1,9),iconc=1,7)/    &
1.0001, 1.0004, 1.0019, 1.0061, 1.0162, 1.0403, 1.0988, 1.2513, 1.7277,  &
1.0000, 1.0001, 1.0004, 1.0021, 1.0078, 1.0218, 1.0553, 1.1378, 1.3632,  &
1.0000, 1.0000, 1.0001, 1.0004, 1.0022, 1.0099, 1.0292, 1.0759, 1.1944,  &
0.9998, 0.9999, 0.9999, 1.0001, 1.0004, 1.0023, 1.0123, 1.0386, 1.1037,  &
0.9994, 0.9994, 0.9995, 0.9996, 0.9999, 1.0005, 1.0023, 1.0147, 1.0497,  &
0.9982, 0.9982, 0.9982, 0.9983, 0.9985, 0.9990, 1.0002, 1.0028, 1.0169,  &
0.9958, 0.9958, 0.9958, 0.9959, 0.9960, 0.9963, 0.9971, 0.9990, 1.0033 /
data ((supersat( 8,iw,iconc,6, 7),iw=1,9),iconc=1,7)/    &
1.0001, 1.0003, 1.0015, 1.0049, 1.0133, 1.0333, 1.0817, 1.2045, 1.5567,  &
1.0000, 1.0001, 1.0003, 1.0016, 1.0062, 1.0177, 1.0452, 1.1123, 1.2880,  &
1.0000, 1.0000, 1.0001, 1.0003, 1.0017, 1.0077, 1.0234, 1.0612, 1.1553,  &
0.9999, 0.9999, 1.0000, 1.0001, 1.0004, 1.0017, 1.0093, 1.0305, 1.0824,  &
0.9995, 0.9995, 0.9995, 0.9996, 0.9999, 1.0004, 1.0017, 1.0109, 1.0387,  &
0.9984, 0.9984, 0.9984, 0.9985, 0.9987, 0.9991, 1.0000, 1.0020, 1.0123,  &
0.9962, 0.9962, 0.9962, 0.9963, 0.9964, 0.9966, 0.9973, 0.9987, 1.0020 /
data ((supersat( 8,iw,iconc,7, 7),iw=1,9),iconc=1,7)/    &
1.0001, 1.0003, 1.0012, 1.0041, 1.0113, 1.0284, 1.0696, 1.1729, 1.4565,  &
1.0000, 1.0001, 1.0003, 1.0013, 1.0051, 1.0148, 1.0381, 1.0947, 1.2399,  &
1.0000, 1.0000, 1.0001, 1.0003, 1.0013, 1.0062, 1.0193, 1.0511, 1.1294,  &
0.9999, 0.9999, 1.0000, 1.0001, 1.0003, 1.0013, 1.0073, 1.0249, 1.0680,  &
0.9995, 0.9995, 0.9996, 0.9997, 0.9999, 1.0003, 1.0013, 1.0083, 1.0312,  &
0.9985, 0.9985, 0.9986, 0.9986, 0.9988, 0.9991, 0.9999, 1.0015, 1.0091,  &
0.9964, 0.9965, 0.9965, 0.9965, 0.9966, 0.9968, 0.9973, 0.9986, 1.0013 /
data ((supersat( 9,iw,iconc,1, 1),iw=1,9),iconc=1,7)/    &
1.0008, 1.0034, 1.0100, 1.0251, 1.0602, 1.1458, 1.3922, 2.6544, 4.9005,  &
1.0002, 1.0008, 1.0043, 1.0138, 1.0361, 1.0882, 1.2218, 1.6675, 3.8483,  &
1.0001, 1.0003, 1.0008, 1.0052, 1.0189, 1.0513, 1.1302, 1.3503, 2.3806,  &
0.9999, 1.0000, 1.0003, 1.0010, 1.0062, 1.0251, 1.0718, 1.1934, 1.5916,  &
0.9991, 0.9992, 0.9995, 1.0001, 1.0014, 1.0070, 1.0318, 1.0978, 1.2904,  &
0.9975, 0.9975, 0.9977, 0.9981, 0.9991, 1.0015, 1.0075, 1.0382, 1.1286,  &
0.9948, 0.9948, 0.9949, 0.9951, 0.9957, 0.9972, 1.0010, 1.0116, 1.0434 /
data ((supersat( 9,iw,iconc,2, 1),iw=1,9),iconc=1,7)/    &
1.0006, 1.0024, 1.0070, 1.0174, 1.0414, 1.0985, 1.2488, 1.7874, 3.6665,  &
1.0002, 1.0006, 1.0030, 1.0096, 1.0248, 1.0600, 1.1457, 1.3895, 2.5289,  &
1.0001, 1.0002, 1.0006, 1.0036, 1.0131, 1.0353, 1.0875, 1.2209, 1.6586,  &
1.0000, 1.0001, 1.0003, 1.0007, 1.0043, 1.0176, 1.0496, 1.1280, 1.3465,  &
0.9997, 0.9997, 0.9999, 1.0002, 1.0009, 1.0047, 1.0227, 1.0682, 1.1876,  &
0.9987, 0.9988, 0.9989, 0.9991, 0.9997, 1.0011, 1.0051, 1.0280, 1.0906,  &
0.9967, 0.9967, 0.9968, 0.9969, 0.9974, 0.9984, 1.0009, 1.0073, 1.0325 /
data ((supersat( 9,iw,iconc,3, 1),iw=1,9),iconc=1,7)/    &
1.0005, 1.0018, 1.0049, 1.0125, 1.0298, 1.0704, 1.1723, 1.4766, 2.7141,  &
1.0002, 1.0005, 1.0021, 1.0068, 1.0176, 1.0426, 1.1018, 1.2570, 1.7941,  &
1.0001, 1.0002, 1.0005, 1.0025, 1.0092, 1.0248, 1.0612, 1.1499, 1.4001,  &
1.0000, 1.0001, 1.0002, 1.0005, 1.0029, 1.0122, 1.0348, 1.0883, 1.2253,  &
0.9999, 0.9999, 1.0000, 1.0002, 1.0006, 1.0031, 1.0159, 1.0478, 1.1271,  &
0.9993, 0.9994, 0.9994, 0.9996, 1.0000, 1.0008, 1.0032, 1.0198, 1.0638,  &
0.9980, 0.9980, 0.9980, 0.9982, 0.9984, 0.9991, 1.0006, 1.0044, 1.0233 /
data ((supersat( 9,iw,iconc,4, 1),iw=1,9),iconc=1,7)/    &
1.0004, 1.0013, 1.0037, 1.0094, 1.0224, 1.0530, 1.1280, 1.3327, 2.0831,  &
1.0001, 1.0004, 1.0016, 1.0049, 1.0130, 1.0316, 1.0753, 1.1849, 1.5074,  &
1.0001, 1.0001, 1.0004, 1.0018, 1.0065, 1.0181, 1.0447, 1.1082, 1.2742,  &
1.0000, 1.0001, 1.0001, 1.0004, 1.0019, 1.0086, 1.0250, 1.0634, 1.1577,  &
0.9999, 1.0000, 1.0000, 1.0002, 1.0004, 1.0020, 1.0110, 1.0341, 1.0898,  &
0.9996, 0.9996, 0.9997, 0.9998, 1.0000, 1.0005, 1.0021, 1.0136, 1.0453,  &
0.9987, 0.9987, 0.9987, 0.9988, 0.9990, 0.9994, 1.0004, 1.0027, 1.0160 /
data ((supersat( 9,iw,iconc,5, 1),iw=1,9),iconc=1,7)/    &
1.0003, 1.0010, 1.0029, 1.0074, 1.0177, 1.0420, 1.1010, 1.2544, 1.7330,  &
1.0001, 1.0003, 1.0012, 1.0037, 1.0100, 1.0245, 1.0586, 1.1424, 1.3704,  &
1.0001, 1.0001, 1.0003, 1.0013, 1.0049, 1.0137, 1.0342, 1.0826, 1.2046,  &
1.0000, 1.0001, 1.0001, 1.0003, 1.0013, 1.0062, 1.0186, 1.0477, 1.1176,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0003, 1.0013, 1.0077, 1.0251, 1.0664,  &
0.9997, 0.9997, 0.9998, 0.9998, 1.0000, 1.0004, 1.0014, 1.0093, 1.0329,  &
0.9990, 0.9990, 0.9990, 0.9991, 0.9992, 0.9995, 1.0002, 1.0017, 1.0107 /
data ((supersat( 9,iw,iconc,6, 1),iw=1,9),iconc=1,7)/    &
1.0002, 1.0008, 1.0023, 1.0060, 1.0145, 1.0347, 1.0833, 1.2067, 1.5602,  &
1.0001, 1.0002, 1.0009, 1.0030, 1.0080, 1.0199, 1.0478, 1.1157, 1.2931,  &
1.0001, 1.0001, 1.0002, 1.0010, 1.0038, 1.0107, 1.0273, 1.0663, 1.1626,  &
1.0000, 1.0001, 1.0001, 1.0002, 1.0010, 1.0046, 1.0144, 1.0375, 1.0926,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0003, 1.0010, 1.0056, 1.0191, 1.0514,  &
0.9998, 0.9998, 0.9998, 0.9999, 1.0000, 1.0003, 1.0010, 1.0064, 1.0247,  &
0.9992, 0.9992, 0.9992, 0.9992, 0.9993, 0.9996, 1.0001, 1.0012, 1.0072 /
data ((supersat( 9,iw,iconc,7, 1),iw=1,9),iconc=1,7)/    &
1.0002, 1.0007, 1.0019, 1.0050, 1.0123, 1.0295, 1.0710, 1.1747, 1.4590,  &
1.0001, 1.0002, 1.0007, 1.0024, 1.0067, 1.0166, 1.0403, 1.0975, 1.2437,  &
1.0000, 1.0001, 1.0002, 1.0008, 1.0030, 1.0088, 1.0226, 1.0552, 1.1351,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0008, 1.0036, 1.0116, 1.0306, 1.0761,  &
1.0000, 1.0000, 1.0000, 1.0001, 1.0002, 1.0008, 1.0041, 1.0151, 1.0414,  &
0.9998, 0.9998, 0.9998, 0.9999, 1.0000, 1.0002, 1.0008, 1.0046, 1.0192,  &
0.9992, 0.9992, 0.9993, 0.9993, 0.9994, 0.9996, 1.0000, 1.0009, 1.0049 /
data ((supersat( 9,iw,iconc,1, 2),iw=1,9),iconc=1,7)/    &
1.0005, 1.0025, 1.0087, 1.0235, 1.0580, 1.1425, 1.3867, 2.6380, 4.8731,  &
1.0002, 1.0005, 1.0028, 1.0117, 1.0331, 1.0838, 1.2145, 1.6527, 3.7952,  &
1.0000, 1.0002, 1.0006, 1.0030, 1.0153, 1.0457, 1.1211, 1.3332, 2.3318,  &
0.9995, 0.9997, 1.0000, 1.0008, 1.0032, 1.0193, 1.0616, 1.1744, 1.5482,  &
0.9982, 0.9983, 0.9986, 0.9992, 1.0007, 1.0043, 1.0231, 1.0798, 1.2504,  &
0.9958, 0.9959, 0.9960, 0.9964, 0.9974, 0.9998, 1.0060, 1.0262, 1.0989,  &
0.9931, 0.9931, 0.9932, 0.9934, 0.9939, 0.9952, 0.9984, 1.0073, 1.0356 /
data ((supersat( 9,iw,iconc,2, 2),iw=1,9),iconc=1,7)/    &
1.0004, 1.0018, 1.0061, 1.0163, 1.0401, 1.0967, 1.2462, 1.7823, 3.6509,  &
1.0001, 1.0004, 1.0020, 1.0081, 1.0229, 1.0574, 1.1419, 1.3829, 2.5099,  &
1.0001, 1.0002, 1.0004, 1.0022, 1.0107, 1.0319, 1.0823, 1.2124, 1.6410,  &
0.9998, 0.9999, 1.0001, 1.0006, 1.0022, 1.0137, 1.0434, 1.1176, 1.3266,  &
0.9992, 0.9992, 0.9994, 0.9998, 1.0006, 1.0024, 1.0167, 1.0571, 1.1663,  &
0.9976, 0.9976, 0.9977, 0.9980, 0.9986, 1.0002, 1.0038, 1.0194, 1.0719,  &
0.9950, 0.9950, 0.9950, 0.9952, 0.9956, 0.9965, 0.9989, 1.0049, 1.0214 /
data ((supersat( 9,iw,iconc,3, 2),iw=1,9),iconc=1,7)/    &
1.0003, 1.0013, 1.0044, 1.0118, 1.0289, 1.0693, 1.1708, 1.4742, 2.7067,  &
1.0001, 1.0003, 1.0015, 1.0058, 1.0164, 1.0409, 1.0996, 1.2536, 1.7873,  &
1.0001, 1.0001, 1.0003, 1.0015, 1.0075, 1.0226, 1.0580, 1.1451, 1.3916,  &
0.9999, 1.0000, 1.0001, 1.0004, 1.0016, 1.0095, 1.0307, 1.0820, 1.2147,  &
0.9996, 0.9996, 0.9997, 1.0000, 1.0005, 1.0017, 1.0116, 1.0407, 1.1146,  &
0.9986, 0.9987, 0.9987, 0.9989, 0.9993, 1.0002, 1.0024, 1.0136, 1.0517,  &
0.9966, 0.9966, 0.9966, 0.9967, 0.9970, 0.9976, 0.9992, 1.0031, 1.0151 /
data ((supersat( 9,iw,iconc,4, 2),iw=1,9),iconc=1,7)/    &
1.0002, 1.0010, 1.0033, 1.0088, 1.0218, 1.0523, 1.1271, 1.3313, 2.0801,  &
1.0001, 1.0002, 1.0011, 1.0042, 1.0121, 1.0304, 1.0738, 1.1828, 1.5039,  &
1.0001, 1.0001, 1.0002, 1.0011, 1.0053, 1.0164, 1.0426, 1.1052, 1.2694,  &
1.0000, 1.0000, 1.0001, 1.0003, 1.0011, 1.0066, 1.0222, 1.0593, 1.1513,  &
0.9998, 0.9998, 0.9999, 1.0000, 1.0003, 1.0012, 1.0079, 1.0292, 1.0818,  &
0.9991, 0.9991, 0.9992, 0.9993, 0.9996, 1.0002, 1.0015, 1.0091, 1.0371,  &
0.9976, 0.9976, 0.9976, 0.9977, 0.9979, 0.9983, 0.9994, 1.0018, 1.0100 /
data ((supersat( 9,iw,iconc,5, 2),iw=1,9),iconc=1,7)/    &
1.0002, 1.0008, 1.0025, 1.0069, 1.0172, 1.0415, 1.1003, 1.2535, 1.7315,  &
1.0001, 1.0002, 1.0008, 1.0032, 1.0093, 1.0237, 1.0576, 1.1411, 1.3683,  &
1.0000, 1.0001, 1.0002, 1.0008, 1.0039, 1.0124, 1.0326, 1.0806, 1.2016,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0008, 1.0047, 1.0165, 1.0447, 1.1133,  &
0.9998, 0.9999, 0.9999, 1.0000, 1.0002, 1.0008, 1.0053, 1.0214, 1.0608,  &
0.9994, 0.9994, 0.9994, 0.9995, 0.9997, 1.0001, 1.0010, 1.0059, 1.0269,  &
0.9981, 0.9981, 0.9981, 0.9982, 0.9983, 0.9986, 0.9994, 1.0011, 1.0062 /
data ((supersat( 9,iw,iconc,6, 2),iw=1,9),iconc=1,7)/    &
1.0002, 1.0006, 1.0020, 1.0056, 1.0142, 1.0343, 1.0828, 1.2060, 1.5592,  &
1.0001, 1.0002, 1.0006, 1.0025, 1.0074, 1.0192, 1.0470, 1.1147, 1.2916,  &
1.0000, 1.0001, 1.0002, 1.0006, 1.0030, 1.0098, 1.0260, 1.0647, 1.1604,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0006, 1.0034, 1.0127, 1.0352, 1.0894,  &
0.9999, 0.9999, 0.9999, 1.0000, 1.0002, 1.0006, 1.0037, 1.0162, 1.0472,  &
0.9995, 0.9995, 0.9995, 0.9996, 0.9997, 1.0000, 1.0007, 1.0039, 1.0201,  &
0.9984, 0.9984, 0.9984, 0.9984, 0.9985, 0.9988, 0.9994, 1.0006, 1.0040 /
data ((supersat( 9,iw,iconc,7, 2),iw=1,9),iconc=1,7)/    &
1.0001, 1.0005, 1.0017, 1.0047, 1.0120, 1.0292, 1.0706, 1.1742, 1.4583,  &
1.0001, 1.0001, 1.0005, 1.0020, 1.0062, 1.0161, 1.0396, 1.0967, 1.2426,  &
1.0000, 1.0001, 1.0001, 1.0005, 1.0024, 1.0080, 1.0216, 1.0540, 1.1334,  &
1.0000, 1.0000, 1.0001, 1.0001, 1.0005, 1.0026, 1.0101, 1.0288, 1.0736,  &
0.9999, 0.9999, 0.9999, 1.0000, 1.0002, 1.0005, 1.0027, 1.0127, 1.0381,  &
0.9995, 0.9995, 0.9995, 0.9996, 0.9997, 1.0000, 1.0006, 1.0027, 1.0155,  &
0.9985, 0.9985, 0.9985, 0.9985, 0.9986, 0.9988, 0.9993, 1.0004, 1.0029 /
data ((supersat( 9,iw,iconc,1, 3),iw=1,9),iconc=1,7)/    &
1.0003, 1.0015, 1.0072, 1.0215, 1.0552, 1.1384, 1.3794, 2.6159, 4.8358,  &
1.0001, 1.0004, 1.0015, 1.0092, 1.0294, 1.0782, 1.2052, 1.6334, 3.7253,  &
0.9998, 1.0000, 1.0005, 1.0016, 1.0113, 1.0392, 1.1099, 1.3117, 2.2686,  &
0.9989, 0.9990, 0.9994, 1.0003, 1.0024, 1.0134, 1.0503, 1.1521, 1.4954,  &
0.9969, 0.9970, 0.9972, 0.9979, 0.9994, 1.0031, 1.0150, 1.0615, 1.2067,  &
0.9941, 0.9941, 0.9942, 0.9946, 0.9955, 0.9976, 1.0033, 1.0199, 1.0718,  &
0.9918, 0.9918, 0.9919, 0.9920, 0.9924, 0.9934, 0.9959, 1.0030, 1.0253 /
data ((supersat( 9,iw,iconc,2, 3),iw=1,9),iconc=1,7)/    &
1.0002, 1.0012, 1.0051, 1.0151, 1.0384, 1.0944, 1.2426, 1.7755, 3.6300,  &
1.0001, 1.0003, 1.0012, 1.0065, 1.0206, 1.0542, 1.1370, 1.3743, 2.4844,  &
1.0000, 1.0001, 1.0003, 1.0012, 1.0080, 1.0278, 1.0759, 1.2015, 1.6183,  &
0.9995, 0.9996, 0.9998, 1.0003, 1.0015, 1.0095, 1.0362, 1.1050, 1.3020,  &
0.9983, 0.9984, 0.9986, 0.9990, 0.9999, 1.0021, 1.0108, 1.0453, 1.1421,  &
0.9960, 0.9961, 0.9962, 0.9964, 0.9970, 0.9986, 1.0023, 1.0118, 1.0538,  &
0.9933, 0.9933, 0.9933, 0.9934, 0.9938, 0.9946, 0.9966, 1.0019, 1.0176 /
data ((supersat( 9,iw,iconc,3, 3),iw=1,9),iconc=1,7)/    &
1.0002, 1.0009, 1.0037, 1.0109, 1.0278, 1.0679, 1.1688, 1.4709, 2.6968,  &
1.0001, 1.0002, 1.0009, 1.0046, 1.0148, 1.0389, 1.0967, 1.2490, 1.7783,  &
1.0000, 1.0001, 1.0002, 1.0009, 1.0056, 1.0198, 1.0541, 1.1389, 1.3805,  &
0.9998, 0.9998, 1.0000, 1.0003, 1.0010, 1.0065, 1.0260, 1.0743, 1.2012,  &
0.9991, 0.9991, 0.9992, 0.9995, 1.0001, 1.0014, 1.0074, 1.0328, 1.1000,  &
0.9975, 0.9975, 0.9976, 0.9977, 0.9982, 0.9992, 1.0015, 1.0079, 1.0394,  &
0.9948, 0.9948, 0.9949, 0.9950, 0.9952, 0.9958, 0.9973, 1.0010, 1.0113 /
data ((supersat( 9,iw,iconc,4, 3),iw=1,9),iconc=1,7)/    &
1.0002, 1.0007, 1.0028, 1.0082, 1.0211, 1.0514, 1.1259, 1.3295, 2.0758,  &
1.0001, 1.0002, 1.0007, 1.0034, 1.0109, 1.0290, 1.0719, 1.1801, 1.4992,  &
1.0000, 1.0001, 1.0002, 1.0007, 1.0039, 1.0145, 1.0399, 1.1013, 1.2631,  &
0.9999, 0.9999, 1.0000, 1.0002, 1.0007, 1.0044, 1.0188, 1.0542, 1.1431,  &
0.9994, 0.9995, 0.9995, 0.9997, 1.0001, 1.0009, 1.0048, 1.0237, 1.0723,  &
0.9983, 0.9983, 0.9984, 0.9985, 0.9988, 0.9995, 1.0010, 1.0050, 1.0284,  &
0.9960, 0.9960, 0.9961, 0.9961, 0.9963, 0.9967, 0.9978, 1.0004, 1.0069 /
data ((supersat( 9,iw,iconc,5, 3),iw=1,9),iconc=1,7)/    &
1.0001, 1.0005, 1.0021, 1.0064, 1.0166, 1.0408, 1.0995, 1.2522, 1.7293,  &
1.0001, 1.0001, 1.0005, 1.0025, 1.0084, 1.0226, 1.0563, 1.1392, 1.3654,  &
1.0000, 1.0001, 1.0001, 1.0005, 1.0029, 1.0109, 1.0306, 1.0779, 1.1975,  &
0.9999, 0.9999, 1.0000, 1.0002, 1.0005, 1.0031, 1.0139, 1.0411, 1.1077,  &
0.9996, 0.9996, 0.9997, 0.9998, 1.0001, 1.0006, 1.0031, 1.0173, 1.0541,  &
0.9987, 0.9987, 0.9987, 0.9988, 0.9991, 0.9995, 1.0006, 1.0033, 1.0206,  &
0.9967, 0.9967, 0.9968, 0.9968, 0.9969, 0.9973, 0.9981, 0.9999, 1.0043 /
data ((supersat( 9,iw,iconc,6, 3),iw=1,9),iconc=1,7)/    &
1.0001, 1.0004, 1.0017, 1.0052, 1.0137, 1.0337, 1.0822, 1.2051, 1.5579,  &
1.0000, 1.0001, 1.0004, 1.0020, 1.0067, 1.0183, 1.0460, 1.1133, 1.2896,  &
1.0000, 1.0001, 1.0001, 1.0004, 1.0021, 1.0085, 1.0245, 1.0627, 1.1575,  &
0.9999, 1.0000, 1.0000, 1.0001, 1.0004, 1.0022, 1.0107, 1.0324, 1.0853,  &
0.9997, 0.9997, 0.9997, 0.9998, 1.0000, 1.0005, 1.0022, 1.0130, 1.0421,  &
0.9989, 0.9989, 0.9989, 0.9990, 0.9992, 0.9996, 1.0004, 1.0024, 1.0152,  &
0.9971, 0.9971, 0.9971, 0.9972, 0.9973, 0.9975, 0.9982, 0.9996, 1.0028 /
data ((supersat( 9,iw,iconc,7, 3),iw=1,9),iconc=1,7)/    &
1.0001, 1.0003, 1.0014, 1.0044, 1.0116, 1.0287, 1.0700, 1.1735, 1.4573,  &
1.0000, 1.0001, 1.0003, 1.0016, 1.0055, 1.0154, 1.0388, 1.0956, 1.2411,  &
1.0000, 1.0000, 1.0001, 1.0003, 1.0017, 1.0069, 1.0203, 1.0523, 1.1311,  &
0.9999, 1.0000, 1.0000, 1.0001, 1.0003, 1.0017, 1.0084, 1.0265, 1.0704,  &
0.9997, 0.9997, 0.9997, 0.9998, 1.0000, 1.0004, 1.0017, 1.0100, 1.0340,  &
0.9990, 0.9990, 0.9990, 0.9991, 0.9992, 0.9995, 1.0002, 1.0018, 1.0115,  &
0.9973, 0.9973, 0.9973, 0.9974, 0.9974, 0.9977, 0.9982, 0.9994, 1.0020 /
data ((supersat( 9,iw,iconc,1, 4),iw=1,9),iconc=1,7)/    &
1.0002, 1.0009, 1.0055, 1.0190, 1.0517, 1.1330, 1.3699, 2.5864, 4.7866,  &
1.0000, 1.0003, 1.0009, 1.0065, 1.0251, 1.0714, 1.1935, 1.6085, 3.6342,  &
0.9993, 0.9996, 1.0001, 1.0014, 1.0073, 1.0319, 1.0967, 1.2854, 2.1883,  &
0.9978, 0.9980, 0.9984, 0.9994, 1.0016, 1.0080, 1.0388, 1.1276, 1.4343,  &
0.9952, 0.9952, 0.9955, 0.9961, 0.9976, 1.0012, 1.0111, 1.0447, 1.1630,  &
0.9925, 0.9926, 0.9927, 0.9929, 0.9937, 0.9955, 1.0002, 1.0142, 1.0598,  &
0.9910, 0.9910, 0.9910, 0.9911, 0.9914, 0.9921, 0.9940, 0.9993, 1.0160 /
data ((supersat( 9,iw,iconc,2, 4),iw=1,9),iconc=1,7)/    &
1.0002, 1.0007, 1.0039, 1.0135, 1.0363, 1.0914, 1.2379, 1.7664, 3.6019,  &
1.0001, 1.0002, 1.0007, 1.0046, 1.0179, 1.0501, 1.1307, 1.3630, 2.4506,  &
0.9997, 0.9999, 1.0002, 1.0009, 1.0052, 1.0231, 1.0683, 1.1880, 1.5892,  &
0.9989, 0.9990, 0.9993, 0.9999, 1.0012, 1.0056, 1.0286, 1.0906, 1.2725,  &
0.9971, 0.9971, 0.9973, 0.9977, 0.9987, 1.0010, 1.0070, 1.0337, 1.1163,  &
0.9943, 0.9943, 0.9944, 0.9946, 0.9952, 0.9966, 1.0001, 1.0096, 1.0379,  &
0.9919, 0.9919, 0.9920, 0.9921, 0.9923, 0.9930, 0.9946, 0.9988, 1.0115 /
data ((supersat( 9,iw,iconc,3, 4),iw=1,9),iconc=1,7)/    &
1.0001, 1.0005, 1.0029, 1.0098, 1.0264, 1.0660, 1.1662, 1.4664, 2.6833,  &
1.0001, 1.0002, 1.0005, 1.0033, 1.0129, 1.0363, 1.0929, 1.2430, 1.7663,  &
0.9999, 1.0000, 1.0002, 1.0006, 1.0036, 1.0166, 1.0492, 1.1311, 1.3661,  &
0.9994, 0.9995, 0.9997, 1.0000, 1.0008, 1.0038, 1.0208, 1.0653, 1.1848,  &
0.9982, 0.9983, 0.9984, 0.9987, 0.9993, 1.0008, 1.0044, 1.0247, 1.0839,  &
0.9959, 0.9959, 0.9960, 0.9962, 0.9966, 0.9976, 1.0000, 1.0062, 1.0281,  &
0.9931, 0.9932, 0.9932, 0.9933, 0.9935, 0.9940, 0.9953, 0.9986, 1.0076 /
data ((supersat( 9,iw,iconc,4, 4),iw=1,9),iconc=1,7)/    &
1.0001, 1.0004, 1.0022, 1.0074, 1.0201, 1.0501, 1.1242, 1.3269, 2.0700,  &
1.0000, 1.0001, 1.0004, 1.0024, 1.0096, 1.0272, 1.0695, 1.1764, 1.4928,  &
1.0000, 1.0000, 1.0001, 1.0004, 1.0025, 1.0122, 1.0365, 1.0963, 1.2547,  &
0.9997, 0.9997, 0.9998, 1.0001, 1.0006, 1.0026, 1.0150, 1.0482, 1.1328,  &
0.9989, 0.9989, 0.9990, 0.9992, 0.9996, 1.0005, 1.0028, 1.0179, 1.0616,  &
0.9970, 0.9970, 0.9971, 0.9972, 0.9975, 0.9982, 0.9999, 1.0039, 1.0202,  &
0.9943, 0.9943, 0.9943, 0.9944, 0.9945, 0.9949, 0.9960, 0.9984, 1.0045 /
data ((supersat( 9,iw,iconc,5, 4),iw=1,9),iconc=1,7)/    &
1.0001, 1.0003, 1.0017, 1.0058, 1.0159, 1.0399, 1.0983, 1.2506, 1.7264,  &
1.0000, 1.0001, 1.0003, 1.0018, 1.0073, 1.0212, 1.0545, 1.1368, 1.3615,  &
1.0000, 1.0000, 1.0001, 1.0003, 1.0018, 1.0091, 1.0282, 1.0744, 1.1920,  &
0.9998, 0.9998, 0.9999, 1.0000, 1.0004, 1.0018, 1.0111, 1.0367, 1.1008,  &
0.9991, 0.9992, 0.9992, 0.9994, 0.9997, 1.0004, 1.0020, 1.0129, 1.0465,  &
0.9976, 0.9977, 0.9977, 0.9978, 0.9980, 0.9985, 0.9997, 1.0025, 1.0145,  &
0.9950, 0.9951, 0.9951, 0.9951, 0.9952, 0.9956, 0.9963, 0.9982, 1.0025 /
data ((supersat( 9,iw,iconc,6, 4),iw=1,9),iconc=1,7)/    &
1.0001, 1.0003, 1.0013, 1.0047, 1.0130, 1.0330, 1.0813, 1.2039, 1.5559,  &
1.0000, 1.0001, 1.0003, 1.0014, 1.0058, 1.0172, 1.0446, 1.1115, 1.2868,  &
1.0000, 1.0000, 1.0001, 1.0003, 1.0014, 1.0071, 1.0225, 1.0600, 1.1536,  &
0.9998, 0.9998, 0.9999, 1.0000, 1.0003, 1.0014, 1.0084, 1.0290, 1.0802,  &
0.9993, 0.9993, 0.9993, 0.9995, 0.9997, 1.0002, 1.0015, 1.0095, 1.0363,  &
0.9979, 0.9980, 0.9980, 0.9981, 0.9982, 0.9987, 0.9996, 1.0017, 1.0104,  &
0.9955, 0.9955, 0.9955, 0.9955, 0.9956, 0.9959, 0.9965, 0.9980, 1.0013 /
data ((supersat( 9,iw,iconc,7, 4),iw=1,9),iconc=1,7)/    &
1.0001, 1.0002, 1.0011, 1.0039, 1.0110, 1.0281, 1.0693, 1.1725, 1.4558,  &
1.0000, 1.0001, 1.0002, 1.0011, 1.0048, 1.0144, 1.0376, 1.0941, 1.2390,  &
1.0000, 1.0000, 1.0001, 1.0002, 1.0011, 1.0057, 1.0186, 1.0501, 1.1281,  &
0.9998, 0.9998, 0.9999, 1.0000, 1.0003, 1.0011, 1.0065, 1.0237, 1.0662,  &
0.9993, 0.9994, 0.9994, 0.9995, 0.9997, 1.0002, 1.0012, 1.0072, 1.0293,  &
0.9981, 0.9981, 0.9981, 0.9982, 0.9983, 0.9987, 0.9995, 1.0012, 1.0076,  &
0.9957, 0.9957, 0.9957, 0.9958, 0.9959, 0.9961, 0.9966, 0.9979, 1.0006 /
data ((supersat( 9,iw,iconc,1, 5),iw=1,9),iconc=1,7)/    &
1.0002, 1.0007, 1.0044, 1.0174, 1.0492, 1.1292, 1.3631, 2.5647, 4.7504,  &
0.9999, 1.0002, 1.0009, 1.0049, 1.0223, 1.0668, 1.1854, 1.5909, 3.5691,  &
0.9989, 0.9992, 0.9998, 1.0011, 1.0052, 1.0275, 1.0883, 1.2679, 2.1322,  &
0.9970, 0.9971, 0.9975, 0.9986, 1.0009, 1.0069, 1.0323, 1.1128, 1.3956,  &
0.9941, 0.9942, 0.9944, 0.9950, 0.9964, 0.9998, 1.0092, 1.0361, 1.1390,  &
0.9918, 0.9918, 0.9919, 0.9922, 0.9928, 0.9944, 0.9985, 1.0107, 1.0505,  &
0.9907, 0.9907, 0.9907, 0.9908, 0.9910, 0.9916, 0.9931, 0.9975, 1.0115 /
data ((supersat( 9,iw,iconc,2, 5),iw=1,9),iconc=1,7)/    &
1.0002, 1.0005, 1.0032, 1.0124, 1.0348, 1.0893, 1.2345, 1.7597, 3.5813,  &
1.0000, 1.0002, 1.0006, 1.0035, 1.0161, 1.0474, 1.1263, 1.3549, 2.4259,  &
0.9995, 0.9997, 1.0000, 1.0008, 1.0037, 1.0201, 1.0632, 1.1787, 1.5688,  &
0.9984, 0.9985, 0.9988, 0.9994, 1.0008, 1.0044, 1.0241, 1.0816, 1.2531,  &
0.9961, 0.9962, 0.9963, 0.9967, 0.9978, 1.0001, 1.0061, 1.0275, 1.1013,  &
0.9933, 0.9933, 0.9934, 0.9936, 0.9941, 0.9954, 0.9987, 1.0075, 1.0350,  &
0.9913, 0.9914, 0.9914, 0.9915, 0.9917, 0.9922, 0.9936, 0.9972, 1.0082 /
data ((supersat( 9,iw,iconc,3, 5),iw=1,9),iconc=1,7)/    &
1.0001, 1.0004, 1.0024, 1.0091, 1.0255, 1.0647, 1.1643, 1.4631, 2.6733,  &
1.0000, 1.0001, 1.0004, 1.0025, 1.0116, 1.0345, 1.0903, 1.2386, 1.7574,  &
0.9998, 0.9999, 1.0001, 1.0006, 1.0026, 1.0146, 1.0459, 1.1257, 1.3559,  &
0.9991, 0.9992, 0.9994, 0.9998, 1.0006, 1.0028, 1.0176, 1.0595, 1.1738,  &
0.9975, 0.9976, 0.9977, 0.9980, 0.9987, 1.0002, 1.0039, 1.0202, 1.0742,  &
0.9949, 0.9949, 0.9950, 0.9951, 0.9955, 0.9965, 0.9989, 1.0049, 1.0223,  &
0.9923, 0.9923, 0.9924, 0.9924, 0.9926, 0.9931, 0.9942, 0.9971, 1.0052 /
data ((supersat( 9,iw,iconc,4, 5),iw=1,9),iconc=1,7)/    &
1.0001, 1.0003, 1.0018, 1.0068, 1.0194, 1.0493, 1.1230, 1.3250, 2.0656,  &
1.0000, 1.0001, 1.0003, 1.0018, 1.0086, 1.0259, 1.0677, 1.1738, 1.4880,  &
0.9999, 1.0000, 1.0001, 1.0004, 1.0019, 1.0106, 1.0343, 1.0928, 1.2488,  &
0.9995, 0.9995, 0.9996, 0.9999, 1.0004, 1.0020, 1.0127, 1.0442, 1.1259,  &
0.9983, 0.9984, 0.9985, 0.9986, 0.9991, 1.0001, 1.0025, 1.0145, 1.0550,  &
0.9961, 0.9961, 0.9962, 0.9963, 0.9966, 0.9973, 0.9990, 1.0030, 1.0159,  &
0.9933, 0.9933, 0.9933, 0.9934, 0.9935, 0.9939, 0.9948, 0.9971, 1.0028 /
data ((supersat( 9,iw,iconc,5, 5),iw=1,9),iconc=1,7)/    &
1.0001, 1.0003, 1.0014, 1.0053, 1.0153, 1.0392, 1.0974, 1.2493, 1.7241,  &
1.0000, 1.0001, 1.0003, 1.0014, 1.0066, 1.0203, 1.0532, 1.1349, 1.3585,  &
0.9999, 1.0000, 1.0001, 1.0003, 1.0014, 1.0080, 1.0265, 1.0719, 1.1881,  &
0.9996, 0.9996, 0.9997, 0.9999, 1.0003, 1.0014, 1.0093, 1.0338, 1.0959,  &
0.9987, 0.9988, 0.9988, 0.9990, 0.9993, 1.0001, 1.0017, 1.0104, 1.0417,  &
0.9968, 0.9968, 0.9968, 0.9969, 0.9972, 0.9977, 0.9990, 1.0018, 1.0112,  &
0.9940, 0.9940, 0.9940, 0.9941, 0.9942, 0.9945, 0.9953, 0.9970, 1.0012 /
data ((supersat( 9,iw,iconc,6, 5),iw=1,9),iconc=1,7)/    &
1.0001, 1.0002, 1.0011, 1.0043, 1.0126, 1.0325, 1.0806, 1.2030, 1.5544,  &
1.0000, 1.0001, 1.0002, 1.0011, 1.0052, 1.0164, 1.0436, 1.1101, 1.2847,  &
0.9999, 1.0000, 1.0001, 1.0002, 1.0011, 1.0061, 1.0212, 1.0581, 1.1507,  &
0.9997, 0.9997, 0.9998, 0.9999, 1.0002, 1.0011, 1.0069, 1.0267, 1.0765,  &
0.9989, 0.9989, 0.9990, 0.9991, 0.9994, 1.0000, 1.0013, 1.0075, 1.0326,  &
0.9972, 0.9972, 0.9972, 0.9973, 0.9975, 0.9979, 0.9989, 1.0011, 1.0078,  &
0.9944, 0.9945, 0.9945, 0.9945, 0.9946, 0.9949, 0.9955, 0.9969, 1.0002 /
data ((supersat( 9,iw,iconc,7, 5),iw=1,9),iconc=1,7)/    &
1.0001, 1.0002, 1.0009, 1.0036, 1.0106, 1.0276, 1.0687, 1.1718, 1.4547,  &
1.0000, 1.0001, 1.0002, 1.0009, 1.0043, 1.0137, 1.0368, 1.0930, 1.2374,  &
0.9999, 1.0000, 1.0000, 1.0002, 1.0009, 1.0049, 1.0175, 1.0486, 1.1258,  &
0.9997, 0.9997, 0.9998, 0.9999, 1.0002, 1.0009, 1.0053, 1.0218, 1.0633,  &
0.9990, 0.9990, 0.9991, 0.9992, 0.9994, 0.9999, 1.0010, 1.0055, 1.0263,  &
0.9974, 0.9974, 0.9974, 0.9975, 0.9976, 0.9980, 0.9989, 1.0007, 1.0057,  &
0.9947, 0.9947, 0.9947, 0.9947, 0.9948, 0.9950, 0.9956, 0.9968, 0.9996 /
data ((supersat( 9,iw,iconc,1, 6),iw=1,9),iconc=1,7)/    &
1.0002, 1.0006, 1.0036, 1.0161, 1.0473, 1.1261, 1.3575, 2.5468, 4.7211,  &
0.9997, 1.0001, 1.0008, 1.0038, 1.0203, 1.0633, 1.1790, 1.5768, 3.5169,  &
0.9986, 0.9988, 0.9995, 1.0009, 1.0040, 1.0244, 1.0821, 1.2545, 2.0882,  &
0.9963, 0.9964, 0.9969, 0.9979, 1.0002, 1.0062, 1.0280, 1.1024, 1.3674,  &
0.9934, 0.9935, 0.9937, 0.9943, 0.9956, 0.9988, 1.0077, 1.0307, 1.1232,  &
0.9914, 0.9914, 0.9915, 0.9917, 0.9923, 0.9937, 0.9974, 1.0084, 1.0441,  &
0.9900, 0.9900, 0.9905, 0.9906, 0.9908, 0.9913, 0.9926, 0.9965, 1.0087 /
data ((supersat( 9,iw,iconc,2, 6),iw=1,9),iconc=1,7)/    &
1.0001, 1.0004, 1.0027, 1.0116, 1.0337, 1.0876, 1.2318, 1.7542, 3.5643,  &
0.9999, 1.0001, 1.0006, 1.0028, 1.0147, 1.0453, 1.1228, 1.3484, 2.4058,  &
0.9993, 0.9995, 0.9999, 1.0007, 1.0029, 1.0180, 1.0593, 1.1716, 1.5527,  &
0.9979, 0.9980, 0.9983, 0.9990, 1.0004, 1.0040, 1.0210, 1.0751, 1.2385,  &
0.9954, 0.9954, 0.9956, 0.9960, 0.9970, 0.9993, 1.0052, 1.0234, 1.0910,  &
0.9927, 0.9927, 0.9928, 0.9930, 0.9935, 0.9947, 0.9976, 1.0059, 1.0314,  &
0.9910, 0.9910, 0.9911, 0.9911, 0.9913, 0.9918, 0.9930, 0.9962, 1.0060 /
data ((supersat( 9,iw,iconc,3, 6),iw=1,9),iconc=1,7)/    &
1.0001, 1.0003, 1.0020, 1.0085, 1.0247, 1.0637, 1.1627, 1.4603, 2.6649,  &
1.0000, 1.0001, 1.0004, 1.0020, 1.0107, 1.0331, 1.0882, 1.2351, 1.7501,  &
0.9997, 0.9998, 1.0000, 1.0005, 1.0021, 1.0131, 1.0434, 1.1214, 1.3477,  &
0.9988, 0.9989, 0.9991, 0.9995, 1.0004, 1.0026, 1.0153, 1.0552, 1.1653,  &
0.9969, 0.9970, 0.9971, 0.9974, 0.9981, 0.9997, 1.0034, 1.0172, 1.0674,  &
0.9942, 0.9942, 0.9942, 0.9944, 0.9948, 0.9957, 0.9980, 1.0038, 1.0186,  &
0.9918, 0.9918, 0.9919, 0.9919, 0.9921, 0.9925, 0.9936, 0.9962, 1.0036 /
data ((supersat( 9,iw,iconc,4, 6),iw=1,9),iconc=1,7)/    &
1.0001, 1.0003, 1.0015, 1.0064, 1.0188, 1.0485, 1.1220, 1.3234, 2.0619,  &
1.0000, 1.0001, 1.0003, 1.0015, 1.0079, 1.0249, 1.0663, 1.1716, 1.4840,  &
0.9998, 0.9999, 1.0000, 1.0004, 1.0015, 1.0095, 1.0325, 1.0900, 1.2440,  &
0.9993, 0.9993, 0.9994, 0.9997, 1.0003, 1.0018, 1.0110, 1.0412, 1.1205,  &
0.9979, 0.9979, 0.9980, 0.9982, 0.9987, 0.9998, 1.0022, 1.0123, 1.0502,  &
0.9954, 0.9954, 0.9954, 0.9955, 0.9958, 0.9966, 0.9983, 1.0022, 1.0132,  &
0.9927, 0.9927, 0.9927, 0.9928, 0.9929, 0.9933, 0.9941, 0.9962, 1.0016 /
data ((supersat( 9,iw,iconc,5, 6),iw=1,9),iconc=1,7)/    &
1.0001, 1.0002, 1.0011, 1.0050, 1.0149, 1.0387, 1.0967, 1.2483, 1.7222,  &
1.0000, 1.0001, 1.0002, 1.0011, 1.0060, 1.0195, 1.0522, 1.1334, 1.3560,  &
0.9999, 0.9999, 1.0000, 1.0003, 1.0012, 1.0071, 1.0252, 1.0699, 1.1849,  &
0.9995, 0.9995, 0.9996, 0.9998, 1.0002, 1.0013, 1.0080, 1.0316, 1.0922,  &
0.9984, 0.9984, 0.9984, 0.9986, 0.9990, 0.9998, 1.0015, 1.0086, 1.0382,  &
0.9961, 0.9962, 0.9962, 0.9963, 0.9965, 0.9971, 0.9983, 1.0012, 1.0091,  &
0.9933, 0.9934, 0.9934, 0.9934, 0.9935, 0.9938, 0.9945, 0.9962, 1.0002 /
data ((supersat( 9,iw,iconc,6, 6),iw=1,9),iconc=1,7)/    &
1.0001, 1.0002, 1.0009, 1.0040, 1.0122, 1.0320, 1.0801, 1.2022, 1.5531,  &
1.0000, 1.0001, 1.0002, 1.0009, 1.0048, 1.0158, 1.0428, 1.1090, 1.2830,  &
0.9999, 0.9999, 1.0000, 1.0002, 1.0009, 1.0054, 1.0201, 1.0566, 1.1483,  &
0.9995, 0.9996, 0.9996, 0.9998, 1.0002, 1.0010, 1.0059, 1.0250, 1.0737,  &
0.9986, 0.9986, 0.9987, 0.9988, 0.9991, 0.9997, 1.0011, 1.0061, 1.0299,  &
0.9966, 0.9966, 0.9966, 0.9967, 0.9968, 0.9973, 0.9983, 1.0006, 1.0064,  &
0.9938, 0.9938, 0.9938, 0.9938, 0.9939, 0.9941, 0.9947, 0.9961, 0.9993 /
data ((supersat( 9,iw,iconc,7, 6),iw=1,9),iconc=1,7)/    &
1.0000, 1.0002, 1.0007, 1.0033, 1.0103, 1.0273, 1.0683, 1.1712, 1.4537,  &
1.0000, 1.0001, 1.0002, 1.0007, 1.0039, 1.0132, 1.0361, 1.0921, 1.2360,  &
0.9999, 0.9999, 1.0000, 1.0002, 1.0007, 1.0043, 1.0166, 1.0473, 1.1239,  &
0.9996, 0.9996, 0.9997, 0.9998, 1.0001, 1.0008, 1.0045, 1.0203, 1.0611,  &
0.9987, 0.9987, 0.9988, 0.9989, 0.9991, 0.9997, 1.0008, 1.0046, 1.0241,  &
0.9968, 0.9968, 0.9968, 0.9969, 0.9970, 0.9974, 0.9983, 1.0002, 1.0049,  &
0.9940, 0.9940, 0.9940, 0.9940, 0.9941, 0.9943, 0.9948, 0.9960, 0.9987 /
data ((supersat( 9,iw,iconc,1, 7),iw=1,9),iconc=1,7)/    &
1.0002, 1.0006, 1.0030, 1.0151, 1.0458, 1.1236, 1.3528, 2.5314, 4.6960,  &
0.9996, 1.0000, 1.0008, 1.0031, 1.0187, 1.0604, 1.1737, 1.5649, 3.4727,  &
0.9982, 0.9985, 0.9992, 1.0006, 1.0042, 1.0221, 1.0771, 1.2436, 2.0516,  &
0.9957, 0.9959, 0.9963, 0.9973, 0.9997, 1.0056, 1.0248, 1.0945, 1.3454,  &
0.9929, 0.9930, 0.9932, 0.9937, 0.9950, 0.9980, 1.0065, 1.0327, 1.1238,  &
0.9912, 0.9912, 0.9913, 0.9915, 0.9920, 0.9932, 0.9966, 1.0067, 1.0394,  &
0.9900, 0.9900, 0.9900, 0.9905, 0.9907, 0.9911, 0.9923, 0.9958, 1.0068 /
data ((supersat( 9,iw,iconc,2, 7),iw=1,9),iconc=1,7)/    &
1.0001, 1.0004, 1.0023, 1.0109, 1.0327, 1.0861, 1.2294, 1.7494, 3.5497,  &
0.9999, 1.0001, 1.0005, 1.0023, 1.0136, 1.0435, 1.1199, 1.3429, 2.3886,  &
0.9992, 0.9993, 0.9997, 1.0006, 1.0025, 1.0163, 1.0563, 1.1656, 1.5392,  &
0.9975, 0.9976, 0.9979, 0.9986, 1.0001, 1.0037, 1.0187, 1.0700, 1.2268,  &
0.9948, 0.9949, 0.9950, 0.9954, 0.9964, 0.9987, 1.0045, 1.0205, 1.0833,  &
0.9923, 0.9923, 0.9923, 0.9925, 0.9930, 0.9941, 0.9969, 1.0046, 1.0286,  &
0.9900, 0.9900, 0.9909, 0.9909, 0.9911, 0.9915, 0.9926, 0.9956, 1.0044 /
data ((supersat( 9,iw,iconc,3, 7),iw=1,9),iconc=1,7)/    &
1.0001, 1.0003, 1.0017, 1.0080, 1.0240, 1.0628, 1.1613, 1.4580, 2.6577,  &
1.0000, 1.0001, 1.0004, 1.0017, 1.0099, 1.0319, 1.0864, 1.2321, 1.7438,  &
0.9996, 0.9997, 0.9999, 1.0004, 1.0018, 1.0119, 1.0414, 1.1179, 1.3409,  &
0.9986, 0.9986, 0.9988, 0.9993, 1.0002, 1.0025, 1.0136, 1.0518, 1.1585,  &
0.9964, 0.9965, 0.9966, 0.9969, 0.9976, 0.9992, 1.0030, 1.0150, 1.0621,  &
0.9936, 0.9936, 0.9937, 0.9938, 0.9942, 0.9951, 0.9973, 1.0029, 1.0193,  &
0.9915, 0.9915, 0.9915, 0.9916, 0.9917, 0.9921, 0.9931, 0.9955, 1.0024 /
data ((supersat( 9,iw,iconc,4, 7),iw=1,9),iconc=1,7)/    &
1.0001, 1.0002, 1.0013, 1.0060, 1.0183, 1.0479, 1.1211, 1.3220, 2.0586,  &
1.0000, 1.0001, 1.0003, 1.0013, 1.0073, 1.0241, 1.0651, 1.1697, 1.4806,  &
0.9998, 0.9998, 1.0000, 1.0003, 1.0013, 1.0086, 1.0311, 1.0877, 1.2399,  &
0.9991, 0.9991, 0.9993, 0.9996, 1.0002, 1.0017, 1.0098, 1.0389, 1.1161,  &
0.9975, 0.9975, 0.9976, 0.9978, 0.9983, 0.9994, 1.0019, 1.0106, 1.0466,  &
0.9948, 0.9948, 0.9949, 0.9950, 0.9953, 0.9960, 0.9976, 1.0016, 1.0125,  &
0.9923, 0.9923, 0.9923, 0.9923, 0.9925, 0.9928, 0.9936, 0.9956, 1.0006 /
data ((supersat( 9,iw,iconc,5, 7),iw=1,9),iconc=1,7)/    &
1.0001, 1.0002, 1.0010, 1.0047, 1.0145, 1.0382, 1.0961, 1.2474, 1.7205,  &
1.0000, 1.0001, 1.0002, 1.0010, 1.0056, 1.0189, 1.0513, 1.1322, 1.3538,  &
0.9998, 0.9999, 1.0000, 1.0002, 1.0010, 1.0064, 1.0241, 1.0682, 1.1821,  &
0.9993, 0.9994, 0.9994, 0.9997, 1.0001, 1.0012, 1.0070, 1.0298, 1.0891,  &
0.9980, 0.9980, 0.9981, 0.9983, 0.9986, 0.9995, 1.0013, 1.0074, 1.0355,  &
0.9956, 0.9956, 0.9956, 0.9957, 0.9959, 0.9965, 0.9978, 1.0007, 1.0083,  &
0.9929, 0.9929, 0.9929, 0.9929, 0.9930, 0.9933, 0.9940, 0.9955, 0.9994 /
data ((supersat( 9,iw,iconc,6, 7),iw=1,9),iconc=1,7)/    &
1.0001, 1.0002, 1.0008, 1.0038, 1.0119, 1.0317, 1.0796, 1.2016, 1.5519,  &
1.0000, 1.0001, 1.0002, 1.0008, 1.0044, 1.0153, 1.0421, 1.1080, 1.2814,  &
0.9999, 0.9999, 1.0000, 1.0002, 1.0008, 1.0048, 1.0193, 1.0553, 1.1463,  &
0.9994, 0.9995, 0.9995, 0.9997, 1.0001, 1.0009, 1.0051, 1.0236, 1.0714,  &
0.9983, 0.9983, 0.9984, 0.9985, 0.9988, 0.9995, 1.0009, 1.0053, 1.0278,  &
0.9960, 0.9960, 0.9961, 0.9961, 0.9963, 0.9968, 0.9978, 1.0002, 1.0059,  &
0.9933, 0.9932, 0.9933, 0.9933, 0.9934, 0.9936, 0.9942, 0.9955, 0.9986 /
data ((supersat( 9,iw,iconc,7, 7),iw=1,9),iconc=1,7)/    &
1.0000, 1.0001, 1.0006, 1.0031, 1.0101, 1.0270, 1.0679, 1.1706, 1.4528,  &
1.0000, 1.0000, 1.0001, 1.0006, 1.0035, 1.0127, 1.0355, 1.0913, 1.2348,  &
0.9999, 0.9999, 1.0000, 1.0002, 1.0007, 1.0038, 1.0159, 1.0463, 1.1223,  &
0.9995, 0.9995, 0.9996, 0.9997, 1.0000, 1.0007, 1.0039, 1.0192, 1.0592,  &
0.9984, 0.9984, 0.9985, 0.9986, 0.9988, 0.9994, 1.0006, 1.0040, 1.0223,  &
0.9963, 0.9963, 0.9963, 0.9963, 0.9965, 0.9969, 0.9978, 0.9998, 1.0044,  &
0.9935, 0.9935, 0.9935, 0.9935, 0.9936, 0.9938, 0.9943, 0.9954, 0.9981 /


nucvalue1 = wtw1*(wtcon1*(jrg1*cldnuctab(rgccn1  ,jw  ,jconcen  ,jtemp,epstab )   &
                        + jrg2*cldnuctab(rgccn1+1,jw  ,jconcen  ,jtemp,epstab))   &
                + wtcon2*(jrg1*cldnuctab(rgccn1  ,jw  ,jconcen+1,jtemp,epstab )   &
                        + jrg2*cldnuctab(rgccn1+1,jw  ,jconcen+1,jtemp,epstab)))  &
           +wtw2*(wtcon1*(jrg1*cldnuctab(rgccn1  ,jw+1,jconcen  ,jtemp,epstab )   &
                        + jrg2*cldnuctab(rgccn1+1,jw+1,jconcen  ,jtemp,epstab))   &
                + wtcon2*(jrg1*cldnuctab(rgccn1  ,jw+1,jconcen+1,jtemp,epstab )   &
                        + jrg2*cldnuctab(rgccn1+1,jw+1,jconcen+1,jtemp,epstab)))

nucvalue2 = wtw1*(wtcon1*(jrg1*cldnuctab(rgccn1  ,jw  ,jconcen  ,jtemp,epstab+1 )   &
                        + jrg2*cldnuctab(rgccn1+1,jw  ,jconcen  ,jtemp,epstab+1))   &
                + wtcon2*(jrg1*cldnuctab(rgccn1  ,jw  ,jconcen+1,jtemp,epstab+1 )   &
                        + jrg2*cldnuctab(rgccn1+1,jw  ,jconcen+1,jtemp,epstab+1)))  &
           +wtw2*(wtcon1*(jrg1*cldnuctab(rgccn1  ,jw+1,jconcen  ,jtemp,epstab+1 )   &
                        + jrg2*cldnuctab(rgccn1+1,jw+1,jconcen  ,jtemp,epstab+1))   &
                + wtcon2*(jrg1*cldnuctab(rgccn1  ,jw+1,jconcen+1,jtemp,epstab+1 )   &
                        + jrg2*cldnuctab(rgccn1+1,jw+1,jconcen+1,jtemp,epstab+1)))

nucvalue = eps1 * nucvalue1 + eps2 * nucvalue2

supsatpcnt = rrv / (1.0001 * rrvlsair)
nucss=0.0
do ss=1,8
  ssvalue11 = wtcon1*(jrg1*supersat(rgccn1  ,ss,jconcen  ,jtemp,epstab ) &
                    + jrg2*supersat(rgccn1+1,ss,jconcen  ,jtemp,epstab)) &
            + wtcon2*(jrg1*supersat(rgccn1  ,ss,jconcen+1,jtemp,epstab ) &
                    + jrg2*supersat(rgccn1+1,ss,jconcen+1,jtemp,epstab))
  ssvalue12 = wtcon1*(jrg1*supersat(rgccn1  ,ss,jconcen  ,jtemp,epstab+1 ) &
                    + jrg2*supersat(rgccn1+1,ss,jconcen  ,jtemp,epstab+1)) &
            + wtcon2*(jrg1*supersat(rgccn1  ,ss,jconcen+1,jtemp,epstab+1 ) &
                    + jrg2*supersat(rgccn1+1,ss,jconcen+1,jtemp,epstab+1))
  ssvalue1  = eps1 * ssvalue11 + eps2 * ssvalue12

  ssvalue21 = wtcon1*(jrg1*supersat(rgccn1  ,ss+1,jconcen  ,jtemp,epstab ) &
                    + jrg2*supersat(rgccn1+1,ss+1,jconcen  ,jtemp,epstab)) &
            + wtcon2*(jrg1*supersat(rgccn1  ,ss+1,jconcen+1,jtemp,epstab ) &
                    + jrg2*supersat(rgccn1+1,ss+1,jconcen+1,jtemp,epstab))
  ssvalue22 = wtcon1*(jrg1*supersat(rgccn1  ,ss+1,jconcen  ,jtemp,epstab+1 ) &
                    + jrg2*supersat(rgccn1+1,ss+1,jconcen  ,jtemp,epstab+1)) &
            + wtcon2*(jrg1*supersat(rgccn1  ,ss+1,jconcen+1,jtemp,epstab+1 ) &
                    + jrg2*supersat(rgccn1+1,ss+1,jconcen+1,jtemp,epstab+1))
  ssvalue2  = eps1 * ssvalue21 + eps2 * ssvalue22

  if((supsatpcnt>=ssvalue1.and.supsatpcnt<=ssvalue2).or. &
     (supsatpcnt>ssvalue2).and.ss==8) then
     mult1=(ssvalue2-supsatpcnt)/(ssvalue2-ssvalue1)
     mult2=(supsatpcnt-ssvalue1)/(ssvalue2-ssvalue1)

     if(ssvalue1==ssvalue2 .or. ssvalue1>ssvalue2 .or. supsatpcnt>ssvalue2) then
        mult1=0.0
        mult2=1.0
     endif

     nucss1 =mult1*(wtcon1*(jrg1*cldnuctab(rgccn1  ,ss  ,jconcen  ,jtemp,epstab ) &
                          + jrg2*cldnuctab(rgccn1+1,ss  ,jconcen  ,jtemp,epstab))   &
                  + wtcon2*(jrg1*cldnuctab(rgccn1  ,ss  ,jconcen+1,jtemp,epstab )   &
                          + jrg2*cldnuctab(rgccn1+1,ss  ,jconcen+1,jtemp,epstab)))  &
            +mult2*(wtcon1*(jrg1*cldnuctab(rgccn1  ,ss+1,jconcen  ,jtemp,epstab )   &
                          + jrg2*cldnuctab(rgccn1+1,ss+1,jconcen  ,jtemp,epstab))   &
                  + wtcon2*(jrg1*cldnuctab(rgccn1  ,ss+1,jconcen+1,jtemp,epstab )   &
                          + jrg2*cldnuctab(rgccn1+1,ss+1,jconcen+1,jtemp,epstab)))
     nucss2 =mult1*(wtcon1*(jrg1*cldnuctab(rgccn1  ,ss  ,jconcen  ,jtemp,epstab+1 ) &
                          + jrg2*cldnuctab(rgccn1+1,ss  ,jconcen  ,jtemp,epstab+1))   &
                  + wtcon2*(jrg1*cldnuctab(rgccn1  ,ss  ,jconcen+1,jtemp,epstab+1 )   &
                          + jrg2*cldnuctab(rgccn1+1,ss  ,jconcen+1,jtemp,epstab+1)))  &
            +mult2*(wtcon1*(jrg1*cldnuctab(rgccn1  ,ss+1,jconcen  ,jtemp,epstab+1 )   &
                          + jrg2*cldnuctab(rgccn1+1,ss+1,jconcen  ,jtemp,epstab+1))   &
                  + wtcon2*(jrg1*cldnuctab(rgccn1  ,ss+1,jconcen+1,jtemp,epstab+1 )   &
                          + jrg2*cldnuctab(rgccn1+1,ss+1,jconcen+1,jtemp,epstab+1)))
     nucss = eps1 * nucss1 + eps2 * nucss2
     go to 102
  endif
enddo
102  continue
tabvalue=max(nucvalue,nucss)
return
end subroutine aero_nuc_tab

!###########################################################################
subroutine cldnuc(m1,k1cnuc,k2cnuc,k1dnuc,k2dnuc,lpw,rv,wp,i,j,dn0 &
          ,nuccldr,nuccldc,nuccldrt,nuccldct)

!use micphys
use mem_grid, only : ngrid

implicit none

integer :: m1,i,j,k,k1cnuc,k2cnuc,k1dnuc,k2dnuc,lpw &
          ,epstab,jw,jconcen,jtemp,rgccn1,ct,ctc,maxct,epsnum
real :: eps1,eps2,wtw1,wtw2,wtcon1,wtcon2,jrg1,jrg2,rg1,rg2 &
       ,vaprccn,vaprgccn,grat,crat,rnuc,excessrv,rcnew,tab,concen_tab &
       ,cvap,gvap,concen_nuc,tairc_nuc,w_nuc,rjw,rjconcen
real, dimension(m1) :: rv,wp,dn0,nuccldr,nuccldc,nuccldrt,nuccldct

k1cnuc = lpw
k2cnuc = 1
k1dnuc = lpw
k2dnuc = 1
maxct=1
if(jnmb(8)>=5) maxct=2

!*********************************************************
!***** If NOT pronosing number concentration of cloud1****
!*********************************************************
if (jnmb(1) == 1 .or. jnmb(1) == 4) then
   rnuc = parm(1) * emb0(1)
   do k = lpw,m1-1
      excessrv = rv(k) - 1.0001 * rvlsair(k)
      rcnew = 0.
      if (excessrv > 0.) then
         rcnew = min(rnuc,.5*excessrv)
         rx(k,1) = rx(k,1) + rcnew
         rv(k) = rv(k) - rcnew
         k2cnuc = k
         cx(k,1) = min(parm(1),rx(k,1) / emb0(1))
      elseif (k2cnuc == 1) then
         k1cnuc = k + 1
      endif
   enddo

!*************************************************************************
!Saleeby(6/3/02) Prognosing number concentration of cloud 1 and/or 2 *****
!*************************************************************************
elseif (jnmb(1) >= 5) then
 do k = lpw,m1-1
   excessrv = rv(k) - 1.0001 * rvlsair(k)
   if (excessrv > 0.) then

!**************CCN MEDIAN RADIUS AND CONCENTRATIONS*************************
    crat = 1.0
    grat = 0.0
    cvap = 0.0
    gvap = 0.0

    if(iccnlev==0) then
     rg1=cnparm
     if(cccnx(k) < 1.0) cccnx(k) = 0.0
     if(jnmb(8)>=5) then
       if(gccnx(k) < 1.e-7) gccnx(k) = 0.0
       cvap= (4.0 * 3.14159 * rg1**2) * cccnx(k)
       rg2=gnparm
       gvap= (4.0 * 3.14159 * rg2**2) * gccnx(k)
     endif
    endif

    if(iccnlev>=1) then
     if(cccnx(k) <= 1.0 .or. cccmx(k) <= 1.e-16) then
       !Do this in case cccnp is a small number or zero
       rg1 = 0.0
       cccnx(k) = 0.0
       cccmx(k) = 0.0
     else
       !From Saleeby&Cotton 2004 Equation 10 for median radius approx.
       !Use avg solute density (1.967g/cm3) of NH42SO4=1.769,NaCl=2.165
       rg1=(0.02523*cccmx(k)/cccnx(k))**.3333
     endif
     if((gccnx(k) <= 1.e-7 .or. gccmx(k) <= 1.e-16) .and. jnmb(8)>=5) then
       !Do this in case gccnp is a small number or zero
       rg2 = 0.0
       gccnx(k) = 0.0
       gccmx(k) = 0.0
     else
      !From Saleeby&Cotton 2004 Equation 10 for median radius approx.
      !Use avg solute density (1.967g/cm3) of NH42SO4=1.769,NaCl=2.165
      rg2=(0.02523*gccmx(k)/gccnx(k))**.3333
     endif
    endif

!*********************CCN TOLERANCES****************************************
    if(iccnlev>=1 .and. cccnx(k)>0.) then
      if(rg1 < 0.01e-4) rg1 = 0.010e-4
      if(rg1 > 0.96e-4) rg1 = 0.960e-4
      cccmx(k) = (rg1**3)*cccnx(k)/0.02523
      cvap=(4.0 * 3.14159 * rg1**2) * cccnx(k)
    endif
!**********************GCCN TOLERANCES**************************************
    if(iccnlev>=1 .and. gccnx(k)>0. .and. jnmb(8)>=5) then
      if(rg2 <= 0.96e-4) rg2 = 0.96e-4
      if(rg2 >  5.00e-4) rg2 = 5.00e-4
      gccmx(k) = (rg2**3)*gccnx(k)/0.02523
      gvap=(4.0 * 3.14159 * rg2**2) * gccnx(k)
    endif

!******************LOOP OVER GCCN AND CCN IF BOTH ARE PRESENT***************
    do ct = 1,maxct
      ctc = 0
!**************TEMPERATURE DEPENDENCY***************************************
      if(ct==1) then
        tairc_nuc = tairc(k)
        if (tairc_nuc < -30.) then
           tairc_nuc = -30.
        elseif (tairc_nuc > 30.) then
           tairc_nuc = 30.
        endif
        jtemp = nint(.1 * (tairc_nuc + 30.)) + 1
!*******VERTICAL VELOCITY DEPENDENCY****************************************
        w_nuc = wp(k)
        if (w_nuc < .010001) then
           w_nuc = .010001
        elseif (w_nuc > 99.99) then
           w_nuc = 99.99
        endif
        rjw = 2. * log10(100. * w_nuc) + 1.
        jw = int(rjw)
        wtw2 = rjw - float(jw)
        wtw1 = 1. - wtw2
!********** CCN NUMBER, MASS, & MEDIAN RADIUS CONSTRAINTS ******************
        concen_nuc = cccnx(k)
        if(concen_nuc > 10000.) then
          print*,"Too many CCN (#/cm3) (grid,concen):",ngrid,concen_nuc
          print*,"at point (k,i,j):",k,i,j
        elseif(concen_nuc < 0.0) then
          print*,"Negative CCN (#/cm3) (grid,concen):",ngrid,concen_nuc
          print*,"at point (k,i,j):",k,i,j
        endif

        !Convert from #/cm3 to #/mg since lookup tables based on #/mg of air
        !Also convert from #/mg to #/kg so we corresspond to cloud droplet
        !number concentration in #/kg
        concen_nuc = concen_nuc / dn0(k) * 1.e6

        if(concen_nuc > 1.0)then
         rjconcen = max(1., min(7., 2. * log10(1.e-7 * concen_nuc) + 1.))
         jconcen = int(rjconcen)
         wtcon2 = rjconcen - float(jconcen)
         wtcon1 = 1. - wtcon2
        endif
!********** MEDIAN RADIUS DEPENDENCY FOR CCN ********************************
        rg=rg1
        if (rg < 0.01e-4) then
           rg = 0.01e-4
        elseif (rg > 0.96e-4) then
           rg = 0.96e-4
        endif
        do rgb=1,maxrg-1
         if((rg>=rg_ccn(rgb)) .and. (rg<=rg_ccn(rgb+1))) then
           rgccn1=rgb
           jrg2 = (rg-rg_ccn(rgb)) / (rg_ccn(rgb+1)-rg_ccn(rgb))
           jrg1 = 1. - jrg2
         endif
        enddo
!********** EPSILON SOLUBILITY FRACTION FOR CCN *****************************
        !Determine weights for interpolating between epsilon table values
        if (epsil < 0.05) then
           epsil = 0.05
        elseif (epsil > 1.00) then
           epsil = 1.00
        endif
        do epsnum=1,maxeps-1
         if((epsil>=epsfrac(epsnum)) .and. (epsil<=epsfrac(epsnum+1))) then
          epstab=epsnum
          eps2 = (epsil-epsfrac(epsnum)) / (epsfrac(epsnum+1)-epsfrac(epsnum))
          eps1 = 1. - eps2
         endif
        enddo
!******************* DETERMINE LOOKUP TABLE VALUES **************************
        if(concen_nuc > 1.0)then
         call aero_nuc_tab(rv(k),rvlsair(k),eps1,eps2,wtw1,wtw2,wtcon1 &
             ,wtcon2,jrg1,jrg2,epstab,jw,jconcen,jtemp,rgccn1,tab)
        else
         tab = 0.0
        endif !concen_nuc positive

        !Get number to nucleate
        concen_tab = concen_nuc * tab

        if(jnmb(8) >= 5 .and. (gvap.ne.0. .or. cvap.ne.0.)) then
         grat = gvap / (gvap + cvap*tab)
         crat = 1.0 - grat
        endif
!****************************************************************************
      endif !if ct=1

!********** GCCN MEDIAN RADIUS AND NUCLEATION *******************************
      if(ct==2) then
        concen_nuc = gccnx(k)
        if(concen_nuc > 100.) then
          print*,"Too many GCCN (#/cm3) (grid,concen):",ngrid,concen_nuc
          print*,"at point (k,i,j):",k,i,j
        elseif(concen_nuc < 0.0) then
          print*,"Negative GCCN (#/cm3) (grid,concen):",ngrid,concen_nuc
          print*,"at point (k,i,j):",k,i,j
        endif

        !Convert from #/cm3 to #/mg since lookup tables based on #/mg of air
        !Also convert from #/mg to #/kg so we corresspond to cloud droplet
        !number concentration in #/kg
        concen_nuc = concen_nuc / dn0(k) * 1.e6
        concen_tab = concen_nuc
!********** MEDIAN RADIUS DEPENDENCY FOR CCN ********************************
        rg = rg2
        if (rg < 0.96e-4) then
           rg = 0.96e-4
        elseif (rg > 5.0e-4) then
           rg = 5.0e-4
        endif
      endif !if ct=2

!*********FOR NUMBER CONCENTRATION PREDICTION OF CLOUD_1*********************
      !If you're not depleting aerosols then nucleate cloud droplets
      !only if CCN and/or GCCN concentration > existing cloud concen
      if((ct==1 .and. iccnlev==0 .and. concen_tab < cx(k,1)) .or. &
         (ct==2 .and. iccnlev==0 .and. concen_tab < cx(k,8)) .or. &
          concen_tab <= 0.0) ctc=1

      if(ctc==0) then

        if(ct==1) then
          if(iccnlev==0) concen_tab = concen_tab - cx(k,1)
          vaprccn = excessrv*crat
          if(concen_tab > vaprccn / emb0(1)) concen_tab = vaprccn / emb0(1)
          if(concen_tab < vaprccn / emb1(1)) vaprccn = concen_tab * emb1(1)
          cx(k,1) = cx(k,1) + concen_tab
          rx(k,1) = rx(k,1) + vaprccn
          if(imbudget >= 1) then
            nuccldr(k) = nuccldr(k) + vaprccn * budget_scale
            nuccldc(k) = nuccldc(k) + concen_tab * 1.e-6
          endif
          if(imbudtot >= 1) then
            nuccldrt(k) = nuccldrt(k) + vaprccn * budget_scalet
            nuccldct(k) = nuccldct(k) + concen_tab * 1.e-6
          endif
        endif

        if(ct==2) then
          if(iccnlev==0) concen_tab = concen_tab - cx(k,8)
          vaprgccn = excessrv*grat
          if(vaprccn < excessrv*crat) &
            vaprgccn = vaprgccn + excessrv*crat - vaprccn
          if(concen_tab > vaprgccn / emb0(8)) concen_tab = vaprgccn / emb0(8)
          if(concen_tab < vaprgccn / emb1(8)) vaprgccn = concen_tab * emb1(8)
          cx(k,8) = cx(k,8) + concen_tab
          rx(k,8) = rx(k,8) + vaprgccn
          if(imbudget >= 1) then
            nuccldr(k) = nuccldr(k) + vaprgccn * budget_scale
            nuccldc(k) = nuccldc(k) + concen_tab * 1.e-6
          endif
          if(imbudtot >= 1) then
            nuccldrt(k) = nuccldrt(k) + vaprgccn * budget_scalet
            nuccldct(k) = nuccldct(k) + concen_tab * 1.e-6
          endif
        endif

        if(iccnlev>=1) then
          ccnnum  = 0.0 !Variable ccnnum is #/cm3
          ccnmass = 0.0 !Variable ccnmass is g/cm3
          concen_tab = concen_tab / 1.e6 * dn0(k) !Convert #/kg to #/cm3
          concen_nuc = concen_nuc / 1.e6 * dn0(k) !Convert #/kg to #/cm3
          !For determining mass of CCN or GCCN to remove, use "concen_nuc"
          if(ct==1) then
            do ic=1,itbin-1
             ccnnum  = jrg1*ccncon(ic,rgccn1  )*concen_nuc &
                     + jrg2*ccncon(ic,rgccn1+1)*concen_nuc
             ccnmass = jrg1*ccnmas(ic,rgccn1  )*concen_nuc &
                     + jrg2*ccnmas(ic,rgccn1+1)*concen_nuc
             if(ccnnum>concen_tab .or. ccnmass>cccmx(k)) then
               ccnnum  = jrg1*ccncon(ic-1,rgccn1  )*concen_nuc &
                       + jrg2*ccncon(ic-1,rgccn1+1)*concen_nuc
               ccnmass = jrg1*ccnmas(ic-1,rgccn1  )*concen_nuc &
                       + jrg2*ccnmas(ic-1,rgccn1+1)*concen_nuc
               go to 101
             endif
            enddo
          endif

          !For determining mass of CCN or GCCN to remove, use "concen_nuc"
          if(ct==2) then
           ccnmass = gccmx(k) * (concen_tab/concen_nuc)
          endif

101       continue
          !If either is zero, set both to zero
          if(concen_tab==0.0 .or. ccnmass==0.0) then
             concen_tab=0.0
             ccnmass=0.0
          endif
          if(ct==1) then
            cccnx(k) = cccnx(k) - concen_tab
            cccmx(k) = cccmx(k) - ccnmass
            if(iccnlev>=2) cnmhx(k,1) = cnmhx(k,1) + ccnmass
          endif
          if(ct==2) then
            gccnx(k) = gccnx(k) - concen_tab
            gccmx(k) = gccmx(k) - ccnmass
            if(iccnlev>=2) cnmhx(k,8) = cnmhx(k,8) + ccnmass
          endif

        endif !CCNLEV CHECK

      endif !CT TYPE OF NUCLEATION
   enddo !Loop over CT 1 to 2
  endif !endif no excess vapor

  if (rx(k,1) .ge. 1.e-12) k2cnuc = k
  if (k2cnuc .eq. 1 .and. rx(k,1) .lt. 1.e-12) k1cnuc = k + 1
  if(jnmb(8)>=5) then
   if (rx(k,8) .ge. 1.e-12) k2dnuc = k
   if (k2dnuc .eq. 1 .and. rx(k,8) .lt. 1.e-12) k1dnuc = k + 1
  endif

 enddo !loop over all vertical levels

else
 print*, 'icloud not allowed to be 2 or 3'
 print*, 'stopping model '
 stop 'icloud'
endif !if number prediction of cloud droplets

return
end subroutine cldnuc

!******************************************************************************

subroutine icenuc(m1,kc1,kc2,kd1,kd2,k1pnuc,k2pnuc,lpw,ngr,rv,dn0,dtlt,i,j &
    ,nucicer,nucicec,inuchomr,inuchomc,inuccontr,inuccontc,inucifnr,inucifnc &
    ,inuchazr,inuchazc,nucicert,nucicect,inuchomrt,inuchomct,inuccontrt &
    ,inuccontct,inucifnrt,inucifnct,inuchazrt,inuchazct)

use rconstants
!use micphys

implicit none

integer :: m1,kc1,kc2,kd1,kd2,k1pnuc,k2pnuc,lpw,ngr,i,j,k,idnc,itc,irhhz,ithz
real :: dn1,fraccld,ridnc,dtlt,ssi0,wdnc2,tc,ritc,wtc2  &
       ,pbvi,ptvi,pdvi,ptotvi,fracifn,cldnuc,cldnucr,rhhz,haznuc  &
       ,rirhhz,wrhhz2,thz,rithz,wthz2,frachaz,ssi,diagni  &
       ,vapnuc,vapnucr,availvap,cont_nuc,homo_nuc,nucfrac,pcthaze

real, dimension(m1) :: rv,dn0,nucicer,nucicec,inuchomr,inuchomc,inuccontr &
                      ,inuccontc,inucifnr,inucifnc,inuchazr,inuchazc &
                      ,nucicert,nucicect,inuchomrt,inuchomct,inuccontrt &
                      ,inuccontct,inucifnrt,inucifnct,inuchazrt,inuchazct

! Define ssi0 to be maximum supersaturation with respect to ice for
! determining total number of IFN that can nucleate in Meyers' formula
data ssi0/0.40/
save

! implement paul's immersion freezing of rain here.  This would
! replace mike's homogeneous freezing of rain which was in h03.

!************************************************************************
!************* CLOUD DROPLET HOMOGENEOUS ICE NUCLEATION******************
!************************************************************************
do k = kc1,kc2

 !If cloud water exists at a minimum quantity
 if (rx(k,1) .gt. 1.e-10) then

   !define dn locally from emb
   dn1 = dnfac(1) * emb(k,1) ** pwmasi(1)
   fraccld=0.
   if (rx(k,1) .gt. 1.e-10 .and. tairc(k) .le. -30.01) then
      ridnc = max(1.,min(float(ndnc-1),dn1 / ddnc))
      idnc = int(ridnc)
      wdnc2 = ridnc - float(idnc)
      tc = max(-49.99,tairc(k))
      ritc = (tc + 50.00) / dtc + 1.0
      itc = int(ritc)
      wtc2 = ritc - float(itc)
      fraccld = (1.-wdnc2) * (1.-wtc2) * fracc(idnc  ,itc  ,ngr)  &
              +     wdnc2  * (1.-wtc2) * fracc(idnc+1,itc  ,ngr)  &
              + (1.-wdnc2) *     wtc2  * fracc(idnc  ,itc+1,ngr)  &
              +     wdnc2  *     wtc2  * fracc(idnc+1,itc+1,ngr)
      if(fraccld > 0.990) fraccld=1.0
   endif

!  Heterogeneous contact ice nucleation of cloud droplets by diffusio-
!  phoresis, thermophoresis, and Brownian motion (transport of IN)
   call contnuc (rx(k,1),cx(k,1),tx(k,1),vap(k,1),press(k)  &
      ,dynvisc(k),thrmcon(k),tair(k),tairc(k)  &
      ,pbvi,ptvi,pdvi,ptotvi,dn1,dtlt,i,k)

! progIFN: Scale ptotvi returned from contnuc by prognosed IFN fraction
!::later   ptotvi = ptotvi * fracifn
! MIKE ADDED THIS COMMENTED ccinp(k)=ccinp(k)-ptotvi, but
! probably do not want sink of ccinp here.
   !Saleeby(2009): Need separate homogeneous freezing options for
   !1-moment and 2-moment cloud and drizzle droplet treatments
   if(jnmb(1) <  5) cldnuc = max(0.,fraccld * cx(k,1) - cx(k,3))
   if(jnmb(1) >= 5) cldnuc = max(0.,fraccld * cx(k,1))

   cont_nuc = ptotvi * emb(k,1)
   homo_nuc = fraccld * rx(k,1)

   if(cont_nuc + homo_nuc > rx(k,1)) then
     nucfrac=rx(k,1)/(cont_nuc + homo_nuc)
     cont_nuc=cont_nuc*nucfrac
     homo_nuc=homo_nuc*nucfrac
     cldnuc=cldnuc*nucfrac
     ptotvi=ptotvi*nucfrac
   endif

   !Transfering aerosol mass from cloud to pristine ice
   if(iccnlev>=2 .and. rx(k,1).gt.1.e-10 .and. tairc(k).le.-2.) then
    rxferratio = min(1.0,(cont_nuc + homo_nuc) / rx(k,1))
    ccnmass  = cnmhx(k,1) * rxferratio
    cnmhx(k,1) = cnmhx(k,1) - ccnmass
    cnmhx(k,3) = cnmhx(k,3) + ccnmass
   endif

   rx(k,3) = rx(k,3) + min(rx(k,1),cont_nuc + homo_nuc)
   rx(k,1) = rx(k,1) - min(rx(k,1),cont_nuc + homo_nuc)
   cx(k,3) = cx(k,3) + min(cx(k,1),cldnuc + ptotvi)
   cx(k,1) = cx(k,1) - min(cx(k,1),cldnuc + ptotvi)

   if(imbudget >= 1) then
     nucicer(k) = nucicer(k) + (cont_nuc + homo_nuc) * budget_scale
     nucicec(k) = nucicec(k) + cldnuc + ptotvi
   endif
   if(imbudget == 2) then
     inuchomr(k)  = inuchomr(k)  + homo_nuc * budget_scale
     inuccontr(k) = inuccontr(k) + cont_nuc * budget_scale
     inuchomc(k)  = inuchomc(k)  + cldnuc
     inuccontc(k) = inuccontc(k) + ptotvi
   endif

   if(imbudtot >= 1) then
     nucicert(k) = nucicert(k) + (cont_nuc + homo_nuc) * budget_scalet
     nucicect(k) = nucicect(k) + cldnuc + ptotvi
   endif
   if(imbudtot == 2) then
     inuchomrt(k)  = inuchomrt(k)  + homo_nuc * budget_scalet
     inuccontrt(k) = inuccontrt(k) + cont_nuc * budget_scalet
     inuchomct(k)  = inuchomct(k)  + cldnuc
     inuccontct(k) = inuccontct(k) + ptotvi
   endif

 endif
enddo

!************************************************************************
!************* DRIZZLE DROPLET HOMOGENEOUS ICE NUCLEATION****************
!************************************************************************
do k = kd1,kd2

 !If drizzle water exists at a minimum quantity
 if (rx(k,8) .gt. 1.e-10) then

   !define dn locally from emb
   dn1 = dnfac(16) * emb(k,8) ** pwmasi(16)
   fraccld=0.
   if (rx(k,8) .gt. 1.e-10 .and. tairc(k) .le. -30.01) then
      ridnc = max(1.,min(float(ndnc-1),dn1 / ddnc))
      idnc = int(ridnc)
      wdnc2 = ridnc - float(idnc)
      tc = max(-49.99,tairc(k))
      ritc = (tc + 50.00) / dtc + 1.0
      itc = int(ritc)
      wtc2 = ritc - float(itc)
      fraccld = (1.-wdnc2) * (1.-wtc2) * fracc(idnc  ,itc  ,ngr)  &
              +     wdnc2  * (1.-wtc2) * fracc(idnc+1,itc  ,ngr)  &
              + (1.-wdnc2) *     wtc2  * fracc(idnc  ,itc+1,ngr)  &
              +     wdnc2  *     wtc2  * fracc(idnc+1,itc+1,ngr)
      if(fraccld > 0.990) fraccld=1.0
   endif

!  Heterogeneous contact ice nucleation of cloud droplets by diffusio-
!  phoresis, thermophoresis, and Brownian motion (transport of IN)
   call contnuc (rx(k,8),cx(k,8),tx(k,8),vap(k,8),press(k)  &
      ,dynvisc(k),thrmcon(k),tair(k),tairc(k)  &
      ,pbvi,ptvi,pdvi,ptotvi,dn1,dtlt,i,k)

! progIFN: Scale ptotvi returned from contnuc by prognosed IFN fraction
!::later   ptotvi = ptotvi * fracifn
! MIKE ADDED THIS COMMENTED ccinp(k)=ccinp(k)-ptotvi, but
! probably do not want sink of ccinp here.
   !Saleeby(2009): Need separate homogeneous freezing options for
   !1-moment and 2-moment cloud and drizzle droplet treatments
   if(jnmb(8) <  5) cldnuc = max(0.,fraccld * cx(k,8) - cx(k,3))
   if(jnmb(8) >= 5) cldnuc = max(0.,fraccld * cx(k,8))

   cont_nuc = ptotvi * emb(k,8)
   homo_nuc = fraccld * rx(k,8)

   if(cont_nuc + homo_nuc > rx(k,8)) then
     nucfrac=rx(k,8)/(cont_nuc + homo_nuc)
     cont_nuc=cont_nuc*nucfrac
     homo_nuc=homo_nuc*nucfrac
     cldnuc=cldnuc*nucfrac
     ptotvi=ptotvi*nucfrac
   endif

   !Transfering aerosol mass from cloud to pristine ice
   if(iccnlev>=2 .and. rx(k,8).gt.1.e-10 .and. tairc(k).le.-2.) then
    rxferratio = min(1.0,(cont_nuc + homo_nuc) / rx(k,8))
    ccnmass  = cnmhx(k,8) * rxferratio
    cnmhx(k,8) = cnmhx(k,8) - ccnmass
    cnmhx(k,3) = cnmhx(k,3) + ccnmass
   endif

   rx(k,3) = rx(k,3) + min(rx(k,8),cont_nuc + homo_nuc)
   rx(k,8) = rx(k,8) - min(rx(k,8),cont_nuc + homo_nuc)
   cx(k,3) = cx(k,3) + min(cx(k,8),cldnuc + ptotvi)
   cx(k,8) = cx(k,8) - min(cx(k,8),cldnuc + ptotvi)

   if(imbudget >= 1) then
     nucicer(k) = nucicer(k) + (cont_nuc + homo_nuc) * budget_scale
     nucicec(k) = nucicec(k) + cldnuc + ptotvi
   endif
   if(imbudget == 2) then
     inuchomr(k)  = inuchomr(k)  + homo_nuc * budget_scale
     inuccontr(k) = inuccontr(k) + cont_nuc * budget_scale
     inuchomc(k)  = inuchomc(k)  + cldnuc
     inuccontc(k) = inuccontc(k) + ptotvi
   endif

   if(imbudtot >= 1) then
     nucicert(k) = nucicert(k) + (cont_nuc + homo_nuc) * budget_scalet
     nucicect(k) = nucicect(k) + cldnuc + ptotvi
   endif
   if(imbudtot == 2) then
     inuchomrt(k)  = inuchomrt(k)  + homo_nuc * budget_scalet
     inuccontrt(k) = inuccontrt(k) + cont_nuc * budget_scalet
     inuchomct(k)  = inuchomct(k)  + cldnuc
     inuccontct(k) = inuccontct(k) + ptotvi
   endif

 endif
enddo

!************************************************************************
!  Homogeneous nucleation of haze
!************************************************************************
k1pnuc = 2
k2pnuc = 1

do k = lpw,m1-1
   rhhz = rv(k) / rvlsair(k)
   haznuc = 0.
   if (rhhz .gt. 0.82 .and. tairc(k) .le. -35.01) then
      rirhhz = min(0.1799,rhhz-0.82) / drhhz + 1.0
      irhhz = int(rirhhz)
      wrhhz2 = rirhhz - float(irhhz)
      thz = max(-59.99,tairc(k))
      rithz = (thz + 60.00) / dthz + 1.0
      ithz = int(rithz)
      wthz2 = rithz - float(ithz)
      frachaz = (1.-wrhhz2) * (1.-wthz2) * frachz(irhhz  ,ithz  )  &
              +     wrhhz2  * (1.-wthz2) * frachz(irhhz+1,ithz  )  &
              + (1.-wrhhz2) *     wthz2  * frachz(irhhz  ,ithz+1)  &
              +     wrhhz2  *     wthz2  * frachz(irhhz+1,ithz+1)
      frachaz = 1. - exp(-frachaz * dtlt)

      !Saleeby(2009): Haze nuclei can be too plentiful here compared
      ! to reality. For 2-moment cloud droplet prediction I scale the
      ! haze nuclei to the CCN concentration. Need better option here.
      if(jnmb(1)>=5) haznuc = frachaz * (cccnx(k) * 1.e6 / dn0(k))
      if(jnmb(1)< 5) haznuc = frachaz * 300.e6
   endif

!  Heterogeneous nucleation by deposition condensation freezing
!  with deposition nuclei.  In 4.3 and beyond, assume that it gives #/kg.
   ssi = min(ssi0,rv(k) / rvisair(k) - 1.)

   !Saleeby(2009): Meyers formula seems over-aggresive. DeMott formula based
   ! on IFN measurement from SPL show much lower nucleation rate. Need better
   ! IFN nucleation scheme, perhaps from ice parcel model.
   if (ssi .gt. 0. .and. tairc(k) .le. -5.) then
      !Meyers formula
      fracifn = exp(12.96 * (ssi - ssi0))
      !DeMott SPL modification
      !fracifn = (10 ** (-4.077421 * ssi + 0.097562)) * fracifn
   else
      fracifn = 0.
   endif

! Diagnose maximum number of IFN to activate (#/kg)
   if (ipris .ge. 5) then
      diagni = fracifn * cifnx(k)
   endif

! Orig Meyers formula: + diagni = exp(6.269 + 12.96 * ssi)
! Combine nucleation types, and limit amounts
! vapnuc is #/kg_air and vapnucr is kg/kg_air

! BEGIN MIKE'S SECTION FOR LIMITING NUMBER OF CRYSTALS NUCLEATED
! BY NUMBER OF ICE CRYSTALS PRESENT ALREADY

   vapnuc = max(0.,haznuc + diagni - cx(k,3))
   vapnucr = vapnuc * emb0(3)
   if (vapnucr .gt. 0.) then
      availvap = .5 * (rv(k) - rvisair(k))
      if (vapnucr .gt. availvap) then
         vapnucr = min(vapnucr, max(0.,availvap))
      endif
   endif
   vapnuc = vapnucr / emb0(3)

!(Saleeby02-21-07) Note that we only want to subtract of the portion
!of IFN that contribute to ice nucleation. Don't remove the IFN until
!we implement a restorative option. At that point, alter "vapnuc" above.
!   if((haznuc.gt.0.0 .or. diagni.gt.0.0) .and. jnmb(3).ge.5) then
!    diagni = (diagni/(haznuc+diagni))*vapnuc
!    cifnx(k) = cifnx(k) - diagni
!   endif

   rx(k,3) = rx(k,3) + vapnucr
   cx(k,3) = cx(k,3) + vapnuc

   pcthaze = haznuc / max(1.e-30,(haznuc + diagni))
   if(imbudget >= 1) then
     nucicer(k) = nucicer(k) + vapnucr * budget_scale
     nucicec(k) = nucicec(k) + vapnuc
   endif
   if(imbudget == 2) then
     inucifnr(k) = vapnucr * (1.0 - pcthaze) * budget_scale
     inucifnc(k) = vapnuc  * (1.0 - pcthaze)
     inuchazr(k) = vapnucr * pcthaze * budget_scale
     inuchazc(k) = vapnuc  * pcthaze
   endif

   if(imbudtot >= 1) then
     nucicert(k) = nucicert(k) + vapnucr * budget_scalet
     nucicect(k) = nucicect(k) + vapnuc
   endif
   if(imbudtot == 2) then
     inucifnrt(k) = inucifnrt(k) + vapnucr * (1.0 - pcthaze) * budget_scalet
     inucifnct(k) = inucifnct(k) + vapnuc  * (1.0 - pcthaze)
     inuchazrt(k) = inuchazrt(k) + vapnucr * pcthaze * budget_scalet
     inuchazct(k) = inuchazct(k) + vapnuc  * pcthaze
   endif

   if (rx(k,3) .gt. 1.e-12) k2pnuc = k
   if (k2pnuc .eq. 1 .and. rx(k,3) .lt. 1.e-12) k1pnuc = k + 1

enddo

! here mike has the habit diagnosis. option 1 is to use habit
! at cloud top, option 2 is to use new habit at each level.
! need to consider other options.  how about method of formation?
! my question about how much of habit is due to existing ice
! structure, and how much is due to current growth environment
! (temp and supsat). relevant supsat is wrt liquid?

return
end subroutine icenuc

!******************************************************************************

subroutine contnuc (rx,cx,tx,vap,press  &
   ,dynvisc,thrmcon,tair,tairc,pbvi,ptvi,pdvi,ptotvi,dn1,dtlt,i,k)

implicit none

integer :: i,k
real :: rx,cx,tx,vap,press,dynvisc,thrmcon,tair,tairc,pbvi,ptvi,pdvi,ptotvi  &
       ,dn1,dtlt,aka,raros,ana,akn,dfar,f1,f2,ft
data aka,raros/5.39e-3,3.e-7/

!  Heterogeneous contact ice nucleation of cloud droplets by diffusio-
!  phoresis, thermophoresis, and Brownian motion (transport of IN)
!
!  ana   = # IN per kg available for contact freezing (from Meyers et al. 1992
!          where ana was interpreted as # per m^3)
!  akn   = Knudsen number (Walko et al. 1995, Eq. 58)
!          [2.28e-5 = mfp * p00 / 293.15]
!  raros = aerosol radius = 3.e-7 m from Cotton et al. (1986)
!  dfar  = aerosol diffusivity (Pruppacher and Klett Eq. 12-15)
!          [7.32e-25 = Boltzmann constant / (6 pi)]
!  f1    = "function 1" (Walko et al. 1995 Eq. 55) multiplied by delta t
!           but now cld concen in #/kg_air so (pvbi, ptvi, pdvi) all per kg_air
!  f2    = "function 2" (Walko et al. 1995 Eq. 56)
!  ft    = "function ft" (Walko et al. 1995 Eq. 57)
!  pbvi  = Brownian motion nucleation amount this timestep [#/kg_air]
!  ptvi  = Thermophoretic nucleation amount this timestep [#/kg_air]
!  pdvi  = Diffusiophoretic nucleation amount this timestep [#/kg_air],
!          reformulated to use vapor diffusion directly.  Factor of 1.2
!          is (1+sigma_va x_a) from Pruppacher and Klett Eq. 12-102
!          divided by .622, the molecular weight ratio between water and air.

   ptotvi = 0.

   if (tx .le. -2. .and. rx .gt. 1.e-10) then

      ana = exp(4.11 - 0.262 * tx)
      akn = 2.28e-5 * tair / (press * raros)
      dfar = 7.32e-25 * tair * (1.+ akn) / (raros * dynvisc)
      f1 = 6.28318 * dn1 * cx * ana * dtlt
      f2 = thrmcon * (tairc - tx) / press
      ft = 0.4 * (1. + 1.45 * akn + 0.4 * akn * exp(-1. / akn))  &
         * (thrmcon + 2.5 * akn * aka)  &
         / ((1. + 3. * akn)  &
         * (2. * thrmcon + 5. * aka * akn + aka))
      pbvi = f1 * dfar
      ptvi = f1 * f2 * ft
      pdvi = 1.2 * ana * vap
      ptotvi = max(0.,pbvi + ptvi + pdvi)

   endif
   return
end subroutine contnuc
!
!###########################################################################

real function gammp(a,x)
implicit none
real :: a,x,gln,gammcf

if(x .lt. a+1.) then
    call gser(gammp,a,x,gln)
else
    call gcf(gammcf,a,x,gln)
    gammp = 1. - gammcf
endif

return
end function

!     *************************************************************

real function gammq(a,x)
implicit none
real :: a,x,gamser,gln

if(x .lt. a+1.) then
    call gser(gamser,a,x,gln)
    gammq = 1. - gamser
else
    call gcf(gammq,a,x,gln)
endif

return
end function

!     ****************************************************************

subroutine gcf(gammcf,a,x,gln)
implicit none
real :: a,x,gammcf,gln

integer, parameter :: itmax=100
real, parameter :: eps=3.e-7
real,external ::gammln

real :: gold,a0,a1,b0,b1,fac,an,ana,anf,gaccel
integer :: n

gln = gammln(a)
gold = 0.
a0 = 1.
a1 = x
b0 = 0.
b1 = 1.
fac = 1.
do n = 1, itmax
    an = float(n)
    ana = an - a
    a0 = (a1+a0*ana)*fac
    b0 = (b1+b0*ana)*fac
    anf = an*fac
    a1 = x*a0 + anf*a1
    b1 = x*b0 + anf*b1
    if (a1 .ne. 0.) then
        fac = 1./a1
        gaccel = b1*fac
        if(abs((gaccel-gold)/gaccel) .lt. eps) goto 20
        gold = gaccel
    endif
enddo
20 continue
gammcf = exp(-x+a*alog(x)-gln)*gaccel
if((-x+a*log(x)-gln) .gt. -38.) then
  gammcf = exp(-x+a*alog(x)-gln)*gaccel
else
  gammcf = 0.
endif
return
end subroutine gcf

!     ****************************************************************

subroutine gser(gamser,a,x,gln)
implicit none
real :: a,x,gamser,gln

integer, parameter :: itmax=100
real, parameter :: eps=3.e-7
real,external ::gammln

real :: ap,sum,del
integer :: n

gln = gammln(a)
if(x .le. 0.) then
    !if(x .lt. 0.) pause
    gamser = 0.
    return
endif
ap = a
sum = 1./a
del = sum
do n = 1, itmax
    ap = ap + 1.
    del = del*x/ap
    sum =  sum  +  del
    if(abs(del) .lt. abs(sum)*eps) goto 20
enddo
20 continue
if((-x+a*log(x)-gln) .gt. -38.) then
  gamser = sum*exp(-x+a*log(x)-gln)
else
  gamser = 0.
endif
return
end subroutine gser

!     ***************************************************************

real function gammln(xx)
implicit none
real :: xx

real(kind=8) :: cof(6),stp
data cof, stp/76.18009173d0, -86.50532033d0, 24.01409822d0,  &
     -1.231739516d0, .120858003d-2, -.536382d-5, 2.50662827465d0/
real(kind=8), parameter :: half=0.5d0, one=1.0d0, fpf=5.5d0

real :: x,tmp,ser
integer :: j

x=xx-one
tmp=x+fpf
tmp=(x+half)*log(tmp)-tmp
ser=one
do j=1,6
    x=x+one
    ser=ser+cof(j)/x
enddo
gammln=tmp+log(stp*ser)
return
end function

!     ****************************************************************

!476
subroutine avint(x,y,n,xlo,xup,ans)
implicit none
integer :: n
real :: x(n),y(n),xlo,xup,ans

real(kind=8) ::r3,rp5,sum,syl,syl2,syl3,syu,syu2,syu3,x1,x2,x3  &
,x12,x13,x23,term1,term2,term3,a,b,c,ca,cb,cc

integer :: i,inlft,inrt,istart,istop
real :: slope,fl,fr

ans =0.0
if(xlo.lt.xup) goto 3
if(xlo.eq.xup) goto 100
if(xlo.gt.xup) goto 200
3 if(n.lt.2) goto 215
do i=2,n
   if(x(i).le.x(i-1)) goto 210
   if(x(i).gt.xup) goto 6
enddo
6 continue
if(n.ge.3) goto 9

!     special n=2 case

slope = (y(2)-y(1))/(x(2)-x(1))
fl = y(1) + slope*(xlo-x(1))
fr = y(2) + slope*(xup-x(2))
ans = 0.5*(fl+fr)*(xup-xlo)
return
9 continue
if(x(n-2).lt.xlo)  goto 205
if(x(3).gt.xup)    goto 205
i = 1
10 if(x(i).ge.xlo) goto 15
i = i+1
goto 10
15 inlft = i
i = n
20 if(x(i).le.xup) goto 25
i = i-1
goto 20
25 inrt = i
if((inrt-inlft).lt.2) goto 205
istart = inlft
if(inlft.eq.1) istart = 2
istop  = inrt
if(inrt.eq.n)  istop  = n-1

r3 = 3.0d0
rp5= 0.5d0
sum = 0.0
syl = xlo
syl2= syl*syl
syl3= syl2*syl

do i=istart,istop
   x1 = x(i-1)
   x2 = x(i)
   x3 = x(i+1)
   x12 = x1-x2
   x13 = x1-x3
   x23 = x2-x3
   term1 = dble(y(i-1))/(x12*x13)
   term2 =-dble(y(i)) /(x12*x23)
   term3 = dble(y(i+1))/(x13*x23)
   a = term1+term2+term3
   b = -(x2+x3)*term1 - (x1+x3)*term2 - (x1+x2)*term3
   c = x2*x3*term1 + x1*x3*term2 + x1*x2*term3
   if(i.le.istart) goto 30
   if(i.gt.istart) goto 35
30    ca = a
   cb = b
   cc = c
   goto 40
35    ca = 0.5*(a+ca)
   cb = 0.5*(b+cb)
   cc = 0.5*(c+cc)
40    syu = x2
   syu2= syu*syu
   syu3= syu2*syu
   sum = sum + ca*(syu3-syl3)/r3 + cb*rp5*(syu2-syl2)  &
             + cc*(syu-syl)
   ca  = a
   cb  = b
   cc  = c
   syl = syu
   syl2= syu2
   syl3= syu3
enddo
syu = xup
ans = sum + ca*(syu**3-syl3)/r3 + cb*rp5*(syu**2-syl2)  &
          + cc*(syu-syl)
100 return
200 print*, 'Upper limit of integration not greater than lower limit.'
stop 'avint2'
205 print*, 'Less than 3 function values between integration limits.'
stop 'avint3'
210 print*, 'Abscissas not strictly increasing.'
stop 'avint4'
215 print*, 'Less than 2 function values were supplied.'
stop 'avint5'
end subroutine avint
!
!###########################################################################

subroutine getict(k1,k2,lcat,i,j,mynum)

!use micphys

implicit none

integer :: k1,k2,lcat,i,j,k,mynum
real :: rict,rictmm

do k = k1,k2
  if (rx(k,lcat) .ge. 1.e-12) then
   rict = dict(lcat) * (log(emb(k,lcat)) - emb0log(lcat)) + 1.
   rictmm = max(rictmin,min(rictmax,rict))
   ict1(k,lcat) = int(rictmm)
   ict2(k,lcat) = ict1(k,lcat) + 1
   wct2(k,lcat) = rictmm - float(ict1(k,lcat))
   wct1(k,lcat) = 1.0 - wct2(k,lcat)
  endif
enddo
return
end subroutine getict

!******************************************************************************

subroutine auto_accret(m1,k1,k2,k3,k4,dn0,dtlt,i,j,cld2rain,cld2raint)

!use micphys

implicit none

integer :: m1,i,j,k,k1,k2,k3,k4,kbot,ktop &
          ,id1cc,id2dd,id1cd,id2cd,id1cr,id1crn,id2cr,id2crn &
          ,ir3cr,id3cr,ir3rr,id3rr

real :: dtlt,dtlt3,dtlt6,dmb1cgs,dmb2cgs,dmb3cgs,r3cgs,en1cgs,en2cgs &
   ,ad1,ad2,ar3,ad3,bd1,bd2,br3,bd3,d3e  &
   ,bd1cc,bd2dd,bd1cd,bd2cd,bd1cr,bd2cr,br3cr,bd3cr,br3rr,bd3rr &
   ,wd1cc,wd2dd,wd1cd,wd2cd,wd1cr,wd2cr,wr3cr,wr3rr,wd3rr &
   ,tm1cc,tn1cc,tn2cc,tm2dd,tn2dd,tn3dd,tm1cd,tn1cd,tm2cd,tn2cd &
   ,tm1cr,tn1cr,tm2cr,tn2cr,tn3rr &
   ,en1cgs_2,en2cgs_2,um1cc,un1cc,un2cc,um2dd,un2dd,un3dd,um1cd &
   ,un1cd,um2cd,un2cd,um1cr,un1cr,um2cr,un2cr,un3rr,um12,um82 &
   ,um18,un1,un8,cfmasi1,cfmasi2,cfmasi3,pwmasi1,pwmasi2,pwmasi3

real, dimension(m1) :: dn0,cld2rain,cld2raint

dtlt3 = 1.e3 * dtlt
dtlt6 = 1.e6 * dtlt

cfmasi1 = 1. / cfmas(1)
cfmasi2 = 1. / cfmas(16)
cfmasi3 = 1. / cfmas(2)

pwmasi1 = 1. / pwmas(1)
pwmasi2 = 1. / pwmas(16)
pwmasi3 = 1. / pwmas(2)

if(jnmb(8) .ne. 0) then
  kbot=min(k1,k3)
  ktop=max(k2,k4)
else
  kbot=k1
  ktop=k2
endif

do k = kbot,ktop
if(rx(k,1) .ge. 1.e-12 .or. rx(k,8) .ge. 1.e-12) then

! This subroutine works in cgs units, so convert inputs from mks
!mean diameter cloud,cloud2,rain - convert (meters) to (cm)
   dmb1cgs = 100. * (emb(k,1) * cfmasi1) ** pwmasi1
   dmb2cgs = 100. * (emb(k,8) * cfmasi2) ** pwmasi2
   dmb3cgs = 100. * (emb(k,2) * cfmasi3) ** pwmasi3

!mixing ratio rain - convert (kg/kg) to (g/cm3)
   r3cgs  = 1.e-3 * rx(k,2) * dn0(k)
   if(rx(k,2) .lt. 1.e-12) r3cgs = 0.0

!number concentration cloud 1 & 2 - convert (#/kg) to (#/cm3)
   en1cgs = 1.e-6 * cx(k,1) * dn0(k)
   en2cgs = 1.e-6 * cx(k,8) * dn0(k)

!max diameter of cloud 1 & 2
   ad1 = max(d1min,min(d1max,dmb1cgs))
   ad2 = max(d2min,min(d2max,dmb2cgs))

!max mixing ratio rain
   ar3 = max(r3min,min(r3max,r3cgs))

!max diameter of rain
   d3minx = max(d3min,(r3cgs / (.1 * .5236)) ** pwmasi3)
   ad3 = max(d3minx,min(d3max,dmb3cgs))

!log of CLOUD 1 & 2 diameter range ratio
   bd1 = alog10(ad1/d1min)
   bd2 = alog10(ad2/d2min)
!log of RAIN mix ratio range ratio
   br3 = alog10(ar3/r3min)
!log of RAIN diamter range ratio
   bd3 = alog10(ad3/d3minx)
!log of RAIN diamter of ratio of given max to max from dmin and mix ratio
   d3e = alog10(d3max/d3minx)

!Calculate ratios of true diameter and mixratio to range in (mkautotab)
!to determine which bin on curve is closest
   bd1cc   = float(ndcc-1) * (ad1 - d1min) / (d1max - d1min) + 1.
   !print*,"bd1cc=",bd1cc,float(ndcc-1),ad1, d1min,d1max , d1min

   bd2dd   = float(ndcc-1) * (ad2 - d2min) / (d2max - d2min) + 1.
   bd1cd   = float(ndcd-1) * (ad1 - d1min) / (d1max - d1min) + 1.
   bd2cd   = float(ndcd-1) * (ad2 - d2min) / (d2max - d2min) + 1.
   bd1cr   = bd1 / d1ecr + 1.
   bd2cr   = bd2 / d2ecr + 1.
   br3cr   = br3 / r3ecr + 1.
   bd3cr   = bd3 / d3e * float(ndrcr-1) + 1.
   !br3rr   = br3 / r3err + 1.                !rain-rain done in "cols"
   !bd3rr   = bd3 / d3e * float(ndrr-1)  + 1. !rain-rain done in "cols"

!Find closest integer bin
   id1cc   =  int(bd1cc)
   id2dd   =  int(bd2dd)
   id1cd   =  int(bd1cd)
   id2cd   =  int(bd2cd)
   id1cr   =  int(bd1cr)
   id1crn  = nint(bd1cr)
   id2cr   =  int(bd2cr)
   id2crn  = nint(bd2cr)
   ir3cr   =  int(br3cr)
   id3cr   = nint(bd3cr)
   !ir3rr   =  int(br3rr) !rain-rain done in "cols"
   !id3rr   =  int(bd3rr) !rain-rain done in "cols"

!Find distance to closest integer bin
   wd1cc   = bd1cc   - float(id1cc)
   wd2dd   = bd2dd   - float(id2dd)
   wd1cd   = bd1cd   - float(id1cd)
   wd2cd   = bd2cd   - float(id2cd)
   wd1cr   = bd1cr   - float(id1cr)
   wd2cr   = bd2cr   - float(id2cr)
   wr3cr   = br3cr   - float(ir3cr)
   !wr3rr   = br3rr   - float(ir3rr) !rain-rain done in "cols"
   !wd3rr   = bd3rr   - float(id3rr) !rain-rain done in "cols"
!To fix array bounds exceptions
   if(id1cc.gt.ndcc-1) then
     id1cc=ndcc-1
     wd1cc=1.-wd1cc
   endif
   if(id2dd.gt.ndcc-1) then
     id2dd=ndcc-1
     wd2dd=1.-wd2dd
   endif
   if(id1cr.gt.ndccr-1) then
     id1cr=ndccr-1
     wd1cr=1.-wd1cr
   endif
   if(ir3cr.gt.nrrcr-1) then
     ir3cr=nrrcr-1
     wr3cr=1.-wr3cr
   endif
   if(id1cd.gt.ndcd-1)  then
     id1cd=ndcd-1
     wd1cd=1.-wd1cd
   endif
   if(id2cd.gt.ndcd-1)  then
     id2cd=ndcd-1
     wd2cd=1.-wd2cd
   endif
   if(id2cr.gt.ndccr-1) then
     id2cr=ndccr-1
     wd2cr=1.-wd2cr
   endif

!***************************************************************************
!cc-> effect m,n of cloud,drizzle and possibly n of rain
tm1cc =   (1.-wd1cc) * r1tabcc(id1cc) + wd1cc * r1tabcc(id1cc+1)
tn1cc =   (1.-wd1cc) * c1tabcc(id1cc) + wd1cc * c1tabcc(id1cc+1)
tn2cc =   (1.-wd1cc) * c2tabcc(id1cc) + wd1cc * c2tabcc(id1cc+1)
!***************************************************************************
!cr-> effect m,n of cloud on rain
tm1cr = (1.-wd1cr) * ((1.-wr3cr) * r1tabcr(id1cr  ,ir3cr  ,id3cr)  &
      +                   wr3cr  * r1tabcr(id1cr  ,ir3cr+1,id3cr)) &
      +     wd1cr  * ((1.-wr3cr) * r1tabcr(id1cr+1,ir3cr  ,id3cr)  &
      +                   wr3cr  * r1tabcr(id1cr+1,ir3cr+1,id3cr))

tn1cr =               (1.-wr3cr) * c1tabcr(id1crn,ir3cr  ,id3cr)   &
      +                   wr3cr  * c1tabcr(id1crn,ir3cr+1,id3cr)
!***************************************************************************
!Currently rain-rain done in "cols"
!rr-> effect n of rain
!tn3rr = (1.-wd3rr) * ((1.-wr3rr) * c3tabrr(ir3rr  ,id3rr  )   &
!      +                   wr3rr  * c3tabrr(ir3rr+1,id3rr  ))  &
!      +     wd3rr  * ((1.-wr3rr) * c3tabrr(ir3rr  ,id3rr+1)   &
!      +                   wr3rr  * c3tabrr(ir3rr+1,id3rr+1))
!***************************************************************************

if(jnmb(8) .ne. 0) then
!***************************************************************************
!dd-> effect m,n of drizzle and n on rain
tm2dd =  (1.-wd2dd) * r2tabdd(id2dd) + wd2dd * r2tabdd(id2dd+1)
tn2dd =  (1.-wd2dd) * c2tabdd(id2dd) + wd2dd * c2tabdd(id2dd+1)
tn3dd =  (1.-wd2dd) * c3tabdd(id2dd) + wd2dd * c3tabdd(id2dd+1)
!***************************************************************************
!cd-> effect m,n of cloud,drizzle and possibly n on rain
!For cloud mixing ratio portion to rain
tm1cd = (1.-wd2cd) * ((1.-wd1cd) * r1tabcd(id1cd  ,id2cd  )  &
       +                  wd1cd  * r1tabcd(id1cd+1,id2cd  )) &
       +    wd2cd  * ((1.-wd1cd) * r1tabcd(id1cd  ,id2cd+1)  &
       +                  wd1cd  * r1tabcd(id1cd+1,id2cd+1))
!For cloud loss of number
tn1cd = (1.-wd2cd) * ((1.-wd1cd) * c1tabcd(id1cd  ,id2cd  )  &
       +                  wd1cd  * c1tabcd(id1cd+1,id2cd  )) &
       +    wd2cd  * ((1.-wd1cd) * c1tabcd(id1cd  ,id2cd+1)  &
       +                  wd1cd  * c1tabcd(id1cd+1,id2cd+1))
!For drizzle mixing ratio portion to rain
tm2cd = (1.-wd2cd) * ((1.-wd1cd) * r2tabcd(id1cd  ,id2cd  )  &
       +                  wd1cd  * r2tabcd(id1cd+1,id2cd  )) &
       +    wd2cd  * ((1.-wd1cd) * r2tabcd(id1cd  ,id2cd+1)  &
       +                  wd1cd  * r2tabcd(id1cd+1,id2cd+1))
!For drizzle loss of number to rain
tn2cd = (1.-wd2cd) * ((1.-wd1cd) * c2tabcd(id1cd  ,id2cd  )  &
       +                  wd1cd  * c2tabcd(id1cd+1,id2cd  )) &
       +    wd2cd  * ((1.-wd1cd) * c2tabcd(id1cd  ,id2cd+1)  &
       +                  wd1cd  * c2tabcd(id1cd+1,id2cd+1))
!***************************************************************************
!dr-> effect m,n of drizzle
tm2cr =  (1.-wd2cr) * ((1.-wr3cr) * r2tabcr(id2cr  ,ir3cr  ,id3cr)  &
      +                    wr3cr  * r2tabcr(id2cr  ,ir3cr+1,id3cr)) &
      +      wd2cr  * ((1.-wr3cr) * r2tabcr(id2cr+1,ir3cr  ,id3cr)  &
      +                    wr3cr  * r2tabcr(id2cr+1,ir3cr+1,id3cr))

tn2cr =                (1.-wr3cr) * c2tabcr(id2crn,ir3cr  ,id3cr)  &
      +                    wr3cr  * c2tabcr(id2crn,ir3cr+1,id3cr)
!***************************************************************************
else
!Prevent uninitialized
tm2dd=0.; tn2dd=0.; tn3dd=0.; tm1cd=0.; tn1cd=0.
tm2cd=0.; tn2cd=0.; tm2cr=0.; tn2cr=0.
endif

!Put tables values in correct units
   en1cgs_2 = en1cgs ** 2
   en2cgs_2 = en2cgs ** 2

   um1cc = tm1cc * en1cgs_2 * dtlt3
   un1cc = tn1cc * en1cgs_2 * dtlt6
   un2cc = tn2cc * en1cgs_2 * dtlt6
   um1cr = 10. ** tm1cr * en1cgs * dtlt3
   un1cr = 10. ** tn1cr * en1cgs * dtlt6
   !un3rr = 10. ** tn3rr * dtlt6 !rain-rain done in "cols"
   !Saleeby(1-13-06) If rain mixing ratio is zero, do not collect cloud
   if(r3cgs == 0.) then
     um1cr=0.
     un1cr=0.
   endif

  if(jnmb(8) .ne. 0) then
   um2dd = tm2dd * en2cgs_2 * dtlt3
   un2dd = tn2dd * en2cgs_2 * dtlt6
   un3dd = tn3dd * en2cgs_2 * dtlt6
   um1cd = tm1cd * en1cgs * en2cgs * dtlt3
   un1cd = tn1cd * en1cgs * en2cgs * dtlt6
   um2cd = tm2cd * en1cgs * en2cgs * dtlt3
   un2cd = tn2cd * en1cgs * en2cgs * dtlt6
   um2cr = 10. ** tm2cr * en2cgs * dtlt3
   un2cr = 10. ** tn2cr * en2cgs * dtlt6
   !Saleeby(1-13-06) If rain mixing ratio is zero, do not collect drizzle
   if(r3cgs == 0.) then
     um2cr=0.
     un2cr=0.
   endif
  endif

! The above values are amounts in kg/m^3 or #/m^3 converted in the
! present timestep, but must still be corrected for the effect of
! density on fall velocity.  Thus, they must be multiplied by
! (dn0i ** .5) which fall velocity is proportional to.  Also, since
! rxfer and enxfer are in units of kg/kg and #/kg, respectively, the
! above transfer amounts must also be multiplied by dn0i.  Together,
! these factors make (dn0i ** 1.5).

!****************** FOR CLOUD1 - CLOUD2 - RAIN TRANSFERS *********************
if(jnmb(8) .eq. 0) then
   um12 = min(rx(k,1),(um1cc + um1cr) * dn0i(k))
   un1 = min(cx(k,1)*dn0(k),(un1cc + un1cr))
   rxfer(k,1,2)  =  rxfer(k,1,2) + um12
   qrxfer(k,1,2) = qrxfer(k,1,2) + um12 * qx(k,1)
   enxfer(k,1,1) = enxfer(k,1,1) + un1 - un2cc
   enxfer(k,1,2) = enxfer(k,1,2) + un2cc

   if(imbudget >= 1) cld2rain(k)  = um12 * budget_scale
   if(imbudtot >= 1) cld2raint(k) = cld2raint(k) + um12 * budget_scalet
endif

if(jnmb(8) .ne. 0) then
   if(um2cd <= 0.0) &
     um12 = min(rx(k,1),(um1cr + um1cd + um2cd) * dn0i(k))
   if(um2cd > 0.0) &
     um12 = min(rx(k,1),(um1cr + um1cd) * dn0i(k))

   if(um2cd <= 0.0) &
     um82 = min(rx(k,8),(um2cr + um2dd) * dn0i(k))
   if(um2cd > 0.0) &
     um82 = min(rx(k,8),(um2cr + um2dd + um2cd) * dn0i(k))

   if(um2cd <= 0.0) &
     um18 = min(rx(k,1),(um1cc + abs(um2cd)) * dn0i(k))
   if(um2cd > 0.0) &
     um18 = min(rx(k,1),(um1cc) * dn0i(k))

    un1 = min(cx(k,1)*dn0(k),(un1cc + un1cd + un1cr))
    un8 = min(cx(k,8)*dn0(k),(un2dd + un2cd + un2cr))

    rxfer(k,1,2) =  rxfer(k,1,2) + um12
    rxfer(k,8,2) =  rxfer(k,8,2) + um82
    rxfer(k,1,8) =  rxfer(k,1,8) + um18
   qrxfer(k,1,2) = qrxfer(k,1,2) + um12 * qx(k,1)
   qrxfer(k,8,2) = qrxfer(k,8,2) + um82 * qx(k,8)
   qrxfer(k,1,8) = qrxfer(k,1,8) + um18 * qx(k,1)
   enxfer(k,1,1) = enxfer(k,1,1) + un1 - un2cc
   enxfer(k,1,8) = enxfer(k,1,8) + un2cc
   enxfer(k,8,8) = enxfer(k,8,8) + un8 - un3dd - un2cd
   enxfer(k,8,2) = enxfer(k,8,2) + un3dd + un2cd

   if(imbudget >= 1) cld2rain(k)  = (um12 + um82) * budget_scale
   if(imbudtot >= 1) cld2raint(k) = cld2raint(k) + (um12 + um82) * budget_scalet
endif

!USING CLOUD_2 TRANSFER TO TRANSFER CC2 COLLISION NUMBER TO
!RAIN SINCE YOU CAN'T JUST SUM THE LOSS OF NUMBER OF CLOUD1 AND CLOUD2
!CLOUD1 COLLISIONS WITH ITSELF DOESN'T AFFECT RAIN NUMBER.

endif !if cloud mixing ratio greater than min threshold
enddo !loop of vertical levels
return
end subroutine auto_accret

!****************************************************************************
subroutine auto_accret_ice(m1,jcat,lcat,k1,k2,dn0,dtlt,i,j &
       ,rimecld,rimecldsnow,rimecldaggr,rimecldgrau,rimecldhail &
       ,rimecldt,rimecldsnowt,rimecldaggrt,rimecldgraut,rimecldhailt)

!use micphys

implicit none

integer :: m1,i,j,k,k1,k2,id1cr,id1crn,irici,idici &
          ,lcat,jhcaty,mx,it,ccat,jcat

real :: dtlt,dtlt3,dtlt6,coeff1,cfmasi1,pwmasi1,dmb1cgs,dmbicgs &
   ,ricgs,en1cgs,ad1,ari,adi,bd1,bri,bdi,die,bd1cr,brici,bdici  &
   ,wd1cr,wrici,tm1ci,tn1ci,um1ci,un1ci,umcld,uncld,trime,urime,rimer &
   ,qcoal,tcoal,fracliq,area,cn13,cn24,sip,rsip,qrsip,coalliq   &
   ,rfinlz,xtoz,xtoy,ytoz,dcmin,dcmax,dcecr,rcoal

real, dimension(m1) :: dn0,rimecld,rimecldsnow,rimecldaggr,rimecldgrau &
                      ,rimecldhail,rimecldt,rimecldsnowt,rimecldaggrt  &
                      ,rimecldgraut,rimecldhailt

integer, dimension(8) :: mcatc

data mcatc /0,0,0,6,6,7,7,0/

dtlt3 = 1.e3 * dtlt
dtlt6 = 1.e6 * dtlt

if(jcat==1)then
  cfmasi1 = 1. / cfmas(1)
  pwmasi1 = 1. / pwmas(1)
  dcmin = d1min
  dcmax = d1max
  dcecr = d1ecr
endif
if(jcat==8)then
  cfmasi1 = 1. / cfmas(16)
  pwmasi1 = 1. / pwmas(16)
  dcmin = d2min
  dcmax = d2max
  dcecr = d2ecr
endif

do k = k1,k2
if(rx(k,jcat).ge.1.e-12 .and. rx(k,lcat).ge.1.e-12) then
if(  ((lcat==4 .or. lcat==5) .and. emb(k,jcat) .gt. 9.0e-13) .or. &
     ((lcat==6 .or. lcat==7) .and. emb(k,jcat) .gt. 3.4e-14) ) then

! This subroutine works in cgs units, so convert inputs from mks
!mean diameter cloud,cloud2,rain - convert (meters) to (cm)
   dmb1cgs = 100. * (emb(k,jcat) * cfmasi1) ** pwmasi1
   dmbicgs = 100. * (emb(k,lcat) / cfmas(lcat)) ** (1./pwmas(lcat))

!mixing ratio ice species - convert (kg/kg) to (g/cm3)
   ricgs  = 1.e-3 * rx(k,lcat) * dn0(k)

!number concentration cloud 1 & 2 - convert (#/kg) to (#/cm3)
   en1cgs = 1.e-6 * cx(k,jcat) * dn0(k)

!max diameter of cloud 1 & 2
   ad1 = max(dcmin,min(dcmax,dmb1cgs))

!max mixing ratio ice species
   ari = max(rimin,min(rimax,ricgs))

!max diameter of ice species
   diminx=max(dimin,((ricgs/1000./cfmas(lcat))**(1./pwmas(lcat)))*200.)
   adi = max(diminx,min(dimax,dmbicgs))

!log of CLOUD 1 & 2 diameter range ratio
   bd1 = alog10(ad1/dcmin)
!log of ICE mix ratio range ratio
   bri = alog10(ari/rimin)
!log of ICE diamter range ratio
   bdi = alog10(adi/diminx)
!log of ICE diamter of ratio of given max to max from dmin and mix ratio
   die = alog10(dimax/diminx)

!Calculate ratios of true diameter and mixratio to range in (mkautotab)
!to determine which bin on curve is closest
   bd1cr   = bd1 / dcecr + 1.
   brici   = bri / rieci + 1.
   bdici   = bdi / die * float(ndrcr-1)  + 1.

!Find closest integer bin
   id1cr   =  int(bd1cr)
   id1crn  = nint(bd1cr)
   irici   =  int(brici)
   idici   = nint(bdici)

!Find distance to closest integer bin
   wd1cr  = bd1cr  - float(id1cr)
   wrici  = brici  - float(irici)

!To fix array bounds exceptions
   if(id1cr.gt.ndccr-1) then
     id1cr=ndccr-1
     wd1cr=1.-wd1cr
   endif
   if(irici.gt.nrrcr-1) then
     irici=nrrcr-1
     wrici=1.-wrici
   endif
!Prevent uninitialized
tm1ci=0.; tn1ci=0.; trime=0.

!***************************************************************************
!ci-> effect m,n of cloud on ice
if(jcat==1)then
 tm1ci = (1.-wd1cr) * ((1.-wrici) * r1tabci(id1cr  ,irici  ,idici,lcat-3)  &
       +                   wrici  * r1tabci(id1cr  ,irici+1,idici,lcat-3)) &
       +     wd1cr  * ((1.-wrici) * r1tabci(id1cr+1,irici  ,idici,lcat-3)  &
       +                   wrici  * r1tabci(id1cr+1,irici+1,idici,lcat-3))
 tn1ci =               (1.-wrici) * c1tabci(id1crn,irici  ,idici,lcat-3)   &
       +                   wrici  * c1tabci(id1crn,irici+1,idici,lcat-3)
 trime = (1.-wd1cr) * ((1.-wrici) * r1rimer(id1cr  ,irici  ,idici,lcat-3)  &
       +                   wrici  * r1rimer(id1cr  ,irici+1,idici,lcat-3)) &
       +     wd1cr  * ((1.-wrici) * r1rimer(id1cr+1,irici  ,idici,lcat-3)  &
       +                   wrici  * r1rimer(id1cr+1,irici+1,idici,lcat-3))
endif
!***************************************************************************

!***************************************************************************
!c2i-> effect m,n of cloud2 on ice
if(jcat==8)then
 tm1ci = (1.-wd1cr) * ((1.-wrici) * r2tabci(id1cr  ,irici  ,idici,lcat-3)  &
       +                   wrici  * r2tabci(id1cr  ,irici+1,idici,lcat-3)) &
       +     wd1cr  * ((1.-wrici) * r2tabci(id1cr+1,irici  ,idici,lcat-3)  &
       +                   wrici  * r2tabci(id1cr+1,irici+1,idici,lcat-3))
 tn1ci =               (1.-wrici) * c2tabci(id1crn,irici  ,idici,lcat-3)   &
       +                   wrici  * c2tabci(id1crn,irici+1,idici,lcat-3)
 trime = (1.-wd1cr) * ((1.-wrici) * r2rimer(id1cr  ,irici  ,idici,lcat-3)  &
       +                   wrici  * r2rimer(id1cr  ,irici+1,idici,lcat-3)) &
       +     wd1cr  * ((1.-wrici) * r2rimer(id1cr+1,irici  ,idici,lcat-3)  &
       +                   wrici  * r2rimer(id1cr+1,irici+1,idici,lcat-3))
endif
!***************************************************************************

 um1ci = 10. ** tm1ci * en1cgs * dtlt3
 un1ci = 10. ** tn1ci * en1cgs * dtlt6
 urime = 10. ** trime * en1cgs * dtlt3

!For Cloud-ice collisions 4=snow, 5=aggregates, 6=graupel, 7=hail
 umcld = min(rx(k,jcat),um1ci * dn0i(k))
 uncld = min(cx(k,jcat)*dn0(k),un1ci)
 rimer = min(rx(k,lcat),urime * dn0i(k))

!*****************************************************************
! Do secondary ice production and transfers due to collection
!****************************************************************
 if(jcat==1) mx=1 !uses gamsip(1,it) for cloud1 set in mic_init.f90
 if(jcat==8) mx=2 !uses gamsip(2,it) for cloud2 set in mic_init.f90

 !Compute fraction liquid of combined ice and cloud1 droplets
 rcoal = umcld + rimer
 qcoal = (qx(k,jcat) * umcld + qx(k,lcat) * rimer) / (1.e-20 + rcoal)

 call qtc(qcoal,tcoal,fracliq)

 coalliq = rcoal * fracliq

 !Secondary Ice Production Based on Hydrometeor Internal Temperature
 if(tcoal.gt.-8.0.and.tcoal.lt.-3.0 .and. (jnmb(jcat)>=5.or.jnmb(lcat)>=5))then
   jhcaty=jhcat(k,lcat)
   area = cx(k,lcat) * dn0(k) * sipfac(jhcaty) * emb(k,lcat)  &
      ** (2.*pwmasi(jhcaty))
   it = nint(emb(k,jcat) / emb1(jcat) * 5000.)
   cn13 = uncld * gamsip13(mx,it) / (area * dtlt)
   cn24 = min(cx(k,jcat)*dn0(k),uncld) * gamsip24(mx,it)
   sip = 9.1e-10 * cn24 * cn13 ** .93
   if (tcoal .lt. -5.) then
      sip = 0.33333 * (tcoal + 8.) * sip
   else
      sip = -0.5 * (tcoal + 3.) * sip
   endif
   rsip = sip * emb0(3) * dn0i(k)
   qrsip = qcoal * rsip
   rcoal = rcoal - rsip
   enxfer(k,jcat,3) = enxfer(k,jcat,3) + sip
   rxfer(k,jcat,3)  = rxfer(k,jcat,3) + rsip
   qrxfer(k,jcat,3) = qrxfer(k,jcat,3) + qrsip
 endif

 !For Cloud-ice collisions 4=snow, 5=aggregates, 6=graupel, 7=hail

 !Saleeby(5-10-2006) For graupel, use 0.5*coalliq so that not all
 ! liquid is sent to hail. Allow some to remain on graupel. The
 ! routine "x02" will handle any further melting of graupel to hail.
 ! This comes from Meyers (1997) two-moment micro scheme. This method
 ! allows melted graupel to transfer to hail if the amount of melted
 ! graupel is greater than the collected liquid cloud water.
 if(lcat.ne.6) rfinlz = min(rcoal,coalliq)
 if(lcat.eq.6) rfinlz = min(rcoal,0.5*coalliq)

 xtoz = min(umcld,rfinlz)
 xtoy = umcld - xtoz
 ytoz = rfinlz - xtoz
 ccat=mcatc(lcat)

 if(imbudget >= 1) rimecld(k) = rimecld(k) + umcld * budget_scale
 if(imbudget == 2) then
    if(lcat.eq.4) rimecldsnow(k) = rimecldsnow(k) + umcld * budget_scale
    if(lcat.eq.5) rimecldaggr(k) = rimecldaggr(k) + umcld * budget_scale
    if(lcat.eq.6) rimecldgrau(k) = rimecldgrau(k) + umcld * budget_scale
    if(lcat.eq.7) rimecldhail(k) = rimecldhail(k) + umcld * budget_scale
 endif

 if(imbudtot >= 1) rimecldt(k) = rimecldt(k) + umcld * budget_scalet
 if(imbudtot == 2) then
    if(lcat.eq.4) rimecldsnowt(k) = rimecldsnowt(k) + umcld * budget_scalet
    if(lcat.eq.5) rimecldaggrt(k) = rimecldaggrt(k) + umcld * budget_scalet
    if(lcat.eq.6) rimecldgraut(k) = rimecldgraut(k) + umcld * budget_scalet
    if(lcat.eq.7) rimecldhailt(k) = rimecldhailt(k) + umcld * budget_scalet
 endif

 rxfer(k,jcat,ccat)  =  rxfer(k,jcat,ccat) + xtoz
 rxfer(k,jcat,lcat)  =  rxfer(k,jcat,lcat) + xtoy
 if(lcat.ne.ccat) rxfer(k,lcat,ccat) = rxfer(k,lcat,ccat) + ytoz

 qrxfer(k,jcat,ccat) = qrxfer(k,jcat,ccat) + qx(k,jcat) * xtoz
 qrxfer(k,jcat,lcat) = qrxfer(k,jcat,lcat) + qx(k,jcat) * xtoy
 if(lcat.ne.ccat) &
   qrxfer(k,lcat,ccat) = qrxfer(k,lcat,ccat) + qx(k,lcat) * ytoz

 enxfer(k,jcat,jcat) = enxfer(k,jcat,jcat) + min(uncld,cx(k,mx))
 if(lcat.ne.ccat) enxfer(k,lcat,ccat) = enxfer(k,lcat,ccat) &
   + ytoz * min(uncld,cx(k,lcat)) / (1.e-20 + rx(k,lcat))

endif !if cloud mean mass is greater than min threshold
endif !if cloud mixing ratio greater than min threshold
enddo !loop of vertical levels

return
end subroutine auto_accret_ice


!******************************************************************************

subroutine effxy(m1,k1,k2,i,j)

!use micphys

implicit none

integer :: m1,i,j,k,ncall7
integer, dimension(11) :: k1,k2
real :: dmr
data ncall7/0/
save

! 1 = rp,rs,ra,rg,rh
if (ncall7 .eq. 0 .and. jnmb(2) .ge. 1 .and. jnmb(3) .ge. 1) then
   ncall7 = 7
   do k = 2,m1-1
      eff(k,1) = 1.0
   enddo
endif

! 2 = cs,ca
! Rough fit from Pruppacher and Klett Fig. 14-14 p. 496:
! close to curve for 404 microns.  Replace with auto_accret eventually.
if (jnmb(2) .ge. 1 .or. jnmb(3) .ge. 1) then
   do k = k1(1),k2(1)
      if (emb(k,1) .gt. 9.e-13) then
         eff(k,2) = min(1.,30. * (emb(k,1) - 9.e-13) ** .15)
      else
         eff(k,2) = 0.
      endif
   enddo
endif

! 3 = rr
if (jnmb(2) .ge. 1) then
   do k = k1(2),k2(2)
    if (rx(k,2) .ge. 1.e-12) then
      if (emb(k,2) .lt. .113e-6) then
         eff(k,3) = 1.0
      elseif (emb(k,2) .gt. .158e-5) then
         eff(k,3) = -5.0
      else
         eff(k,3) = 2. - exp(.1326e7 * (emb(k,2) - .113e-6))
      endif
    endif
   enddo
endif

! 4 = pp,ps,pa
if (jnmb(5) .ge. 1) then
   do k = k1(3),k2(3)
      if (abs(tx(k,3)+14.) .le. 2.) then
         eff(k,4) = 1.4
      else
         eff(k,4) = min(0.2,10. ** (0.035 * tx(k,3) - 0.7))
      endif
   enddo

! 5 = ss,sa
   do k = k1(4),k2(4)
      if (abs(tx(k,4)+14.) .le. 2.) then
         eff(k,5) = 1.4
      else
         eff(k,5) = min(0.2,10. ** (0.035 * tx(k,4) - 0.7))
      endif
   enddo

! 6 = aa
   do k = k1(5),k2(5)
    if (rx(k,5) .ge. 1.e-12) then
      if (abs(tx(k,5)+14.) .le. 2.) then
         eff(k,6) = 1.4
      elseif (tx(k,5) .ge. -1.) then
         eff(k,6) = 1.
      else
         eff(k,6) = min(0.2,10. ** (0.035 * tx(k,5) - 0.7))
      endif
    endif
   enddo
endif

! 7 = pg,sg,ag,gg,gh
if (jnmb(6) .ge. 1) then
   do k = k1(6),k2(6)
      if (qr(k,6) .gt. 0.) then
         eff(k,7) = 1.0
      else
         eff(k,7) = min(0.2,10. ** (0.035 * tx(k,6) - 0.7))
      endif
   enddo
endif

! 8 = ph,sh,ah,gh
if (jnmb(7) .ge. 1) then
   do k = k1(7),k2(7)
    if (rx(k,7) .ge. 1.e-12) then
      if (qr(k,7) .gt. 0.) then
         eff(k,8) = 1.0
      else
         eff(k,8) = min(0.2,10. ** (0.035 * tx(k,7) - 0.7))
      endif
    endif
   enddo
endif

! 9 = cg,ch
! Rough fit from Pruppacher and Klett Fig. 14-11 p. 485:
! close to curves for 142 and 305 microns.  Replace with auto_accret eventually.
if (jnmb(2) .ge. 1 .or. jnmb(3) .ge. 1) then
   do k = k1(1),k2(1)
      if (emb(k,1) .gt. 3.4e-14) then
         eff(k,9) = min(1.,1426. * (emb(k,1) - 3.4e-14) ** .28)
      else
         eff(k,9) = 0.
      endif
   enddo
endif

! 10 = hh (trial)
if (jnmb(7) .ge. 1) then
   do k = k1(7),k2(7)
      eff(k,10) = max(0.,.1 + .005 * tx(k,7))
   enddo
endif

! 11 = ds,da
! Rough fit from Pruppacher and Klett Fig. 14-14 p. 496:
! close to curve for 404 microns.  Replace with auto_accret eventually.
if (jnmb(2) .ge. 1 .or. jnmb(3) .ge. 1) then
   do k = k1(8),k2(8)
      if (emb(k,8) .gt. 9.e-13) then
         eff(k,11) = min(1.,30. * (emb(k,8) - 9.e-13) ** .15)
      else
         eff(k,11) = 0.
      endif
   enddo
endif

! 12 = dg,dh
! Rough fit from Pruppacher and Klett Fig. 14-11 p. 485:
! close to curves for 142 and 305 microns.  Replace with auto_accret eventually.
if (jnmb(2) .ge. 1 .or. jnmb(3) .ge. 1) then
   do k = k1(8),k2(8)
      if (emb(k,8) .gt. 3.4e-14) then
         eff(k,12) = min(1.,1426. * (emb(k,8) - 3.4e-14) ** .28)
      else
         eff(k,12) = 0.
      endif
   enddo
endif

return
end subroutine effxy

!******************************************************************************

subroutine cols(m1,mx,mc1,k1,k2,i,j)

!use micphys

implicit none

integer :: ipc,m1,mx,mc1,k1,k2,i,j,k
real :: colnum,tabval

do k = k1,k2
if(rx(k,mx) .ge. 1.e-12) then
   ipc = ipairc(jhcat(k,mx),jhcat(k,mx))

   tabval  &
   = wct1(k,mx) ** 2               * coltabc(ict1(k,mx),ict1(k,mx),ipc)  &
   + 2. * wct1(k,mx) * wct2(k,mx) * coltabc(ict1(k,mx),ict2(k,mx),ipc)  &
   + wct2(k,mx) ** 2               * coltabc(ict2(k,mx),ict2(k,mx),ipc)

   colnum = colfacc(k) * eff(k,mc1) * cx(k,mx) ** 2 * 10. ** (-tabval)
   enxfer(k,mx,mx) = enxfer(k,mx,mx) + min(0.5 * cx(k,mx),colnum)
endif
enddo
return
end subroutine cols

!******************************************************************************

subroutine col3344(m1,mx,mz,mc1,k1,k2,i,j,aggregate,aggrselfpris,aggrselfsnow &
                  ,aggregatet,aggrselfprist,aggrselfsnowt)

!use micphys

implicit none

integer :: m1,mx,mz,mc1,k1,k2,i,j,k,ip,ipc
real :: c1,tabvalx,colamt,tabvaln,colnum
real, dimension(m1) :: aggregate,aggrselfpris,aggrselfsnow &
                      ,aggregatet,aggrselfprist,aggrselfsnowt

do k = k1,k2
if(rx(k,mx) .ge. 1.e-12) then
   ip = ipairr(jhcat(k,mx),jhcat(k,mx))
   ipc = ipairc(jhcat(k,mx),jhcat(k,mx))
   c1 = eff(k,mc1) * cx(k,mx) ** 2

   tabvalx  &
    = wct1(k,mx) ** 2               * coltabr(ict1(k,mx),ict1(k,mx),ip)  &
    + 2. * wct1(k,mx) * wct2(k,mx) * coltabr(ict1(k,mx),ict2(k,mx),ip)  &
    + wct2(k,mx) ** 2               * coltabr(ict2(k,mx),ict2(k,mx),ip)

   colamt = min(rx(k,mx),colfacr2(k) * c1 * 10. ** (-tabvalx))
   rxfer(k,mx,mz) = rxfer(k,mx,mz) + colamt
   qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + colamt * qx(k,mx)

   if(imbudget >= 1) aggregate(k) = aggregate(k) + colamt*budget_scale
   if(imbudget == 2) then
    if(mx.eq.3) aggrselfpris(k) = colamt*budget_scale
    if(mx.eq.4) aggrselfsnow(k) = colamt*budget_scale
   endif

   if(imbudtot >= 1) aggregatet(k) = aggregatet(k) + colamt*budget_scalet
   if(imbudtot == 2) then
    if(mx.eq.3) aggrselfprist(k) = aggrselfprist(k) + colamt*budget_scalet
    if(mx.eq.4) aggrselfsnowt(k) = aggrselfsnowt(k) + colamt*budget_scalet
   endif

   if (jnmb(mz) >= 5) then

   tabvaln  &
   = wct1(k,mx) ** 2               * coltabc(ict1(k,mx),ict1(k,mx),ipc)  &
   + 2. * wct1(k,mx) * wct2(k,mx) * coltabc(ict1(k,mx),ict2(k,mx),ipc)  &
   + wct2(k,mx) ** 2               * coltabc(ict2(k,mx),ict2(k,mx),ipc)

   colnum = min(0.5 * cx(k,mx),colfacc2(k) * c1 * 10. ** (-tabvaln))
   enxfer(k,mx,mz) = enxfer(k,mx,mz) + colnum
   enxfer(k,mx,mx) = enxfer(k,mx,mx) + colnum

   endif
endif
enddo
return
end subroutine col3344

!******************************************************************************

subroutine col3443(m1,mx,my,mz,k1,k2,i,j,aggregate,aggrprissnow &
                  ,aggregatet,aggrprissnowt)

!use micphys

implicit none

integer :: m1,mx,my,mz,k1,k2,i,j,k,jhcatx,jhcaty,ipxy,ipyx,ipc
real :: c1,tabvalx,rcx,tabvaly,rcy,tabvaln,colnum
real, dimension(m1) :: aggregate,aggrprissnow,aggregatet,aggrprissnowt

do k = k1,k2
if(rx(k,mx) .ge. 1.e-12 .and. rx(k,my) .ge. 1.e-12) then
   jhcatx = jhcat(k,mx)
   jhcaty = jhcat(k,my)
   ipxy = ipairr(jhcatx,jhcaty)
   ipyx = ipairr(jhcaty,jhcatx)
   ipc  = ipairc(jhcatx,jhcaty)
   c1 = eff(k,4) * cx(k,mx) * cx(k,my)

   tabvalx  &
     = wct1(k,mx) * wct1(k,my) * coltabr (ict1(k,mx),ict1(k,my),ipxy)  &
     + wct2(k,mx) * wct1(k,my) * coltabr (ict2(k,mx),ict1(k,my),ipxy)  &
     + wct1(k,mx) * wct2(k,my) * coltabr (ict1(k,mx),ict2(k,my),ipxy)  &
     + wct2(k,mx) * wct2(k,my) * coltabr (ict2(k,mx),ict2(k,my),ipxy)
   rcx = min(rx(k,mx),c1 * colfacr(k) * 10. ** (-tabvalx))

   tabvaly  &
     = wct1(k,my) * wct1(k,mx) * coltabr (ict1(k,my),ict1(k,mx),ipyx)  &
     + wct2(k,my) * wct1(k,mx) * coltabr (ict2(k,my),ict1(k,mx),ipyx)  &
     + wct1(k,my) * wct2(k,mx) * coltabr (ict1(k,my),ict2(k,mx),ipyx)  &
     + wct2(k,my) * wct2(k,mx) * coltabr (ict2(k,my),ict2(k,mx),ipyx)
   rcy = min(rx(k,my),c1 * colfacr(k) * 10. ** (-tabvaly))

   rxfer(k,mx,mz) = rxfer(k,mx,mz) + rcx
   qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + rcx * qx(k,mx)

   rxfer(k,my,mz) = rxfer(k,my,mz) + rcy
   qrxfer(k,my,mz) = qrxfer(k,my,mz) + rcy * qx(k,my)

   if(imbudget >= 1) aggregate(k) = aggregate(k) + (rcx + rcy)*budget_scale
   if(imbudget == 2) aggrprissnow(k) = (rcx + rcy)*budget_scale

   if(imbudtot >= 1) aggregatet(k) = aggregatet(k) + (rcx + rcy)*budget_scalet
   if(imbudtot == 2) aggrprissnowt(k) = aggrprissnowt(k) + (rcx + rcy)*budget_scalet

   tabvaln  &
       = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc)  &
       + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc)  &
       + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc)  &
       + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)
   colnum = c1 * colfacc(k) * 10. ** (-tabvaln)

   if (cx(k,mx) .gt. cx(k,my)) then
      enxfer(k,my,mz) = enxfer(k,my,mz) + min(cx(k,my),colnum)
      enxfer(k,mx,mx) = enxfer(k,mx,mx) + min(cx(k,mx),colnum)
   else
      enxfer(k,mx,mz) = enxfer(k,mx,mz) + min(cx(k,mx),colnum)
      enxfer(k,my,my) = enxfer(k,my,my) + min(cx(k,my),colnum)
   endif

! also loss for aerosol

endif
enddo
return
end subroutine col3443

!******************************************************************************

subroutine col1(m1,mx,my,mz,mc4,k1,k2,i,j)

!use micphys

implicit none

integer :: m1,mx,my,mz,mc4,k1,k2,i,j,k,ipxy,ipc
real :: c1,tabvalx,rcx,tabvaln,colnum

do k = k1,k2
if(rx(k,mx) .ge. 1.e-12 .and. rx(k,my) .ge. 1.e-12) then
   ipxy = ipairr(jhcat(k,mx),jhcat(k,my))
   ipc  = ipairc(jhcat(k,mx),jhcat(k,my))
   c1 = eff(k,mc4) * cx(k,mx) * cx(k,my)

   tabvalx  &
     = wct1(k,mx) * wct1(k,my) * coltabr (ict1(k,mx),ict1(k,my),ipxy)  &
     + wct2(k,mx) * wct1(k,my) * coltabr (ict2(k,mx),ict1(k,my),ipxy)  &
     + wct1(k,mx) * wct2(k,my) * coltabr (ict1(k,mx),ict2(k,my),ipxy)  &
     + wct2(k,mx) * wct2(k,my) * coltabr (ict2(k,mx),ict2(k,my),ipxy)

   rcx = min(rx(k,mx),c1 * colfacr(k) * 10. ** (-tabvalx))
   rxfer(k,mx,mz) = rxfer(k,mx,mz) + rcx
   qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + rcx * qx(k,mx)

   if (jnmb(mx) >= 5) then
      tabvaln  &
        = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc)  &
        + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc)  &
        + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc)  &
        + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)

      colnum = c1 * colfacc(k) * 10. ** (-tabvaln)
      enxfer(k,mx,mx) = enxfer(k,mx,mx) + min(colnum,cx(k,mx))

! also loss for aerosol

   endif

endif
enddo
return
end subroutine col1

!******************************************************************************

subroutine col2(m1,mx,my,mz,mc2,k1,k2,dn0,dtlt,i,j,rimecld &
               ,rimecldsnow,rimecldaggr,rimecldgrau,rimecldhail &
               ,rimecldt,rimecldsnowt,rimecldaggrt,rimecldgraut,rimecldhailt)

use rconstants
!use micphys

implicit none

integer :: m1,mx,my,mz,mc2,k1,k2,i,j,k,jhcatx,jhcaty,ipxy,ipyx,ipc,it,mxx
real :: c1,c2,tabvalx,rcx,tabvaly,rcy,tabvaln,colnum0,colnum,rcoal  &
       ,qrcx,qrcy,qrcoal,qcoal,fracliq,tcoal,coalliq,coalice,area,cn13,cn24  &
       ,sip,rsip,qrsip,rfinlz,xtoz,dtlt

real, dimension(m1) :: dn0 &
          ,rimecld,rimecldsnow,rimecldaggr,rimecldgrau,rimecldhail &
          ,rimecldt,rimecldsnowt,rimecldaggrt,rimecldgraut,rimecldhailt

real, dimension(16) ::  alpha,beta
!            1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
data alpha /00.,00.,00., 1., 1., 1., 1.,00.,00.,00.,00., 1., 1., 1., 1.,00./
data beta  /00.,00.,00.,1.5,1.1,0.0,0.0,00.,00.,00.,00.,1.2,1.1,1.1,1.3,00./

do k = k1,k2
if(rx(k,mx) .ge. 1.e-12 .and. rx(k,my) .ge. 1.e-12) then
   jhcatx = jhcat(k,mx)
   jhcaty = jhcat(k,my)
   ipxy = ipairr(jhcatx,jhcaty)
   ipyx = ipairr(jhcaty,jhcatx)
   ipc  = ipairc(jhcatx,jhcaty)
   c2 = cx(k,mx) * cx(k,my)
   c1 = eff(k,mc2) * c2

   tabvalx  &
     = wct1(k,mx) * wct1(k,my) * coltabr (ict1(k,mx),ict1(k,my),ipxy)  &
     + wct2(k,mx) * wct1(k,my) * coltabr (ict2(k,mx),ict1(k,my),ipxy)  &
     + wct1(k,mx) * wct2(k,my) * coltabr (ict1(k,mx),ict2(k,my),ipxy)  &
     + wct2(k,mx) * wct2(k,my) * coltabr (ict2(k,mx),ict2(k,my),ipxy)

   rcx = min(rx(k,mx),c1 * colfacr(k) * 10. ** (-tabvalx))

   tabvaly  &
     = wct1(k,my) * wct1(k,mx) * coltabr (ict1(k,my),ict1(k,mx),ipyx)  &
     + wct2(k,my) * wct1(k,mx) * coltabr (ict2(k,my),ict1(k,mx),ipyx)  &
     + wct1(k,my) * wct2(k,mx) * coltabr (ict1(k,my),ict2(k,mx),ipyx)  &
     + wct2(k,my) * wct2(k,mx) * coltabr (ict2(k,my),ict2(k,mx),ipyx)

   rcy = min(rx(k,my),c1 * colfacr(k) * 10. ** (-tabvaly))

   if (jnmb(mx) >= 5 .or. jnmb(my) >= 5) then

      tabvaln  &
       = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc)  &
       + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc)  &
       + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc)  &
       + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)

      colnum0 = c2 * colfacc(k) * 10. ** (-tabvaln)
      colnum = colnum0 * eff(k,mc2)
   endif

   rcoal = rcx + rcy
   qrcx = rcx * qx(k,mx)
   qrcy = rcy * qx(k,my)
   qrcoal = qrcx + qrcy

   !Saleeby (3-15-06) Lowering the threshold to prevent over-freezing
   qcoal = qrcoal / (1.e-20 + rcoal)

   call qtc(qcoal,tcoal,fracliq)

   coalliq = rcoal * fracliq
   coalice = rcoal - coalliq

! Secondary ice production: cn24 is the number fraction of collected cloud
! droplets larger than 24 microns and is obtained from an incomplete gamma
! function table. cn13 is the fraction of collected cloud droplets
! smaller than 13 microns. "area" is cross section area of collecting ice
! per m^3 of atmospheric volume.

! Saleeby(6/3/02): Hallett-Mossop is done for both cloud droplet modes, though
! contribution from the large droplet mode is minimal compared to the small
! droplet mode. Ice splintering is only done if number concentration is
! prognostic for at least one of the two hydromet species involved. This is
! specified above in the calculations for "colnum".

   if (tcoal.gt.-8.0 .and. tcoal.lt.-3.0 .and. (jnmb(mx)>=5.or.jnmb(my)>=5)) then
      if(mx==1) mxx=1 !uses gamsip(1,it) for cloud1 set in mic_init.f90
      if(mx==8) mxx=2 !uses gamsip(2,it) for cloud2 set in mic_init.f90
      area = cx(k,my) * dn0(k) * sipfac(jhcaty) * emb(k,my)  &
         ** (2.*pwmasi(jhcaty))
      it = nint(emb(k,mx) / emb1(mx) * 5000.)
      cn13 = colnum * gamsip13(mxx,it) / (area * dtlt)
      cn24 = min(cx(k,mx)*dn0(k),colnum0) * gamsip24(mxx,it)
      sip = 9.1e-10 * cn24 * cn13 ** .93
      if (tcoal .lt. -5.) then
         sip = 0.33333 * (tcoal + 8.) * sip
      else
         sip = -0.5 * (tcoal + 3.) * sip
      endif
      rsip = sip * emb0(3) * dn0i(k)
      qrsip = qcoal * rsip
      rcoal = rcoal - rsip
      qrcoal = qrcoal - qrsip
      enxfer(k,mx,3) = enxfer(k,mx,3) + sip
      rxfer(k,mx,3) = rxfer(k,mx,3) + rsip
      qrxfer(k,mx,3) = qrxfer(k,mx,3) + qrsip
   endif

! ALWAYS NEED (ALPHA + BETA) .GE. 1 but in the (rare) case that
! fracliq may be a little larger than fracx due to collected
! liquid being above 0C, need (ALPHA + BETA) to be at least 1.1
! or 1.2, or need ALPHA itself to be at least 1.0.

   rfinlz = min(rcoal,  &
      alpha(jhcaty) * coalliq + beta(jhcaty) * rcx)

   xtoz = min(rcx,rfinlz)

   if(imbudget >= 1) rimecld(k) = rimecld(k) + rcx * budget_scale
   if(imbudget == 2) then
    if(my.eq.4) rimecldsnow(k) = rimecldsnow(k) + rcx * budget_scale
    if(my.eq.5) rimecldaggr(k) = rimecldaggr(k) + rcx * budget_scale
    if(my.eq.6) rimecldgrau(k) = rimecldgrau(k) + rcx * budget_scale
    if(my.eq.7) rimecldhail(k) = rimecldhail(k) + rcx * budget_scale
   endif

   if(imbudtot >= 1) rimecldt(k) = rimecldt(k) + rcx * budget_scalet
   if(imbudtot == 2) then
    if(my.eq.4) rimecldsnowt(k) = rimecldsnowt(k) + rcx * budget_scalet
    if(my.eq.5) rimecldaggrt(k) = rimecldaggrt(k) + rcx * budget_scalet
    if(my.eq.6) rimecldgraut(k) = rimecldgraut(k) + rcx * budget_scalet
    if(my.eq.7) rimecldhailt(k) = rimecldhailt(k) + rcx * budget_scalet
   endif

   rxfer(k,mx,mz) = rxfer(k,mx,mz) + xtoz
   rxfer(k,mx,my) = rxfer(k,mx,my) + rcx - xtoz
   if (my .ne. mz) rxfer(k,my,mz) = rxfer(k,my,mz)  &
      + rfinlz - xtoz

   qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + qx(k,mx) * xtoz
   qrxfer(k,mx,my) = qrxfer(k,mx,my) + qx(k,mx) * (rcx - xtoz)
   if (my .ne. mz) qrxfer(k,my,mz) = qrxfer(k,my,mz)  &
      + qx(k,my) * (rfinlz - xtoz)

   enxfer(k,mx,mx) = enxfer(k,mx,mx) + min(colnum,cx(k,mx))
   if (my .ne. mz) enxfer(k,my,mz) = enxfer(k,my,mz)  &
      + (rfinlz - xtoz) * min(colnum,cx(k,my)) / (1.e-20 + rcy)

! BUT NEED TO CHANGE THE ABOVE FOR 177 COLLECTION BECAUSE X = Y
! also include loss of aerosol

endif
enddo
return
end subroutine col2

!******************************************************************************

subroutine col3(m1,mx,my,mz,k1,k2,i,j,ice2rain,rain2ice,rain2pr,rain2sn &
               ,rain2ag,rain2gr,rain2ha,rain2ha_xtra,ice2raint,rain2icet &
               ,rain2prt,rain2snt,rain2agt,rain2grt,rain2hat,rain2ha_xtrat)

!use micphys

implicit none

integer :: m1,mx,my,mz,k1,k2,i,j,k,ipxy,ipyx,ipc,jhcaty
real :: c1,tabvalx,rcx,tabvaly,rcy,tabvaln,colnum,colnumx,colnumy,coalnum  &
       ,rcoal,qrcx,qrcy,qrcoal,qcoal,fracliq,coalliq,coalice,xtoz  &
       ,rfinlz,tcoal,cfinlz
real, dimension(16) :: alpha,beta
real, dimension(m1) :: ice2rain,rain2ice,rain2pr,rain2sn,rain2ag,rain2gr &
                      ,rain2ha,rain2ha_xtra,ice2raint,rain2icet,rain2prt &
                      ,rain2snt,rain2agt,rain2grt,rain2hat,rain2ha_xtrat

!            1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
data alpha /00.,00., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,00./
data beta  /00.,00., 2., 2., 2., 1., 0., 2., 2., 2., 2., 2., 2., 2., 2.,00./

!Note: mx=2 (rain), my=lcat, mz=7 (hail)
do k = k1,k2
if(rx(k,mx) .ge. 1.e-12 .and. rx(k,my) .ge. 1.e-12) then
   jhcaty = jhcat(k,my)
   ipxy = ipairr(jhcat(k,mx),jhcaty)
   ipyx = ipairr(jhcaty,jhcat(k,mx))
   ipc  = ipairc(jhcat(k,mx),jhcaty)
   c1 = eff(k,1) * cx(k,mx) * cx(k,my)

   tabvalx  &
     = wct1(k,mx) * wct1(k,my) * coltabr (ict1(k,mx),ict1(k,my),ipxy)  &
     + wct2(k,mx) * wct1(k,my) * coltabr (ict2(k,mx),ict1(k,my),ipxy)  &
     + wct1(k,mx) * wct2(k,my) * coltabr (ict1(k,mx),ict2(k,my),ipxy)  &
     + wct2(k,mx) * wct2(k,my) * coltabr (ict2(k,mx),ict2(k,my),ipxy)

   rcx = min(rx(k,mx),c1 * colfacr(k) * 10. ** (-tabvalx))

   tabvaly  &
     = wct1(k,my) * wct1(k,mx) * coltabr (ict1(k,my),ict1(k,mx),ipyx)  &
     + wct2(k,my) * wct1(k,mx) * coltabr (ict2(k,my),ict1(k,mx),ipyx)  &
     + wct1(k,my) * wct2(k,mx) * coltabr (ict1(k,my),ict2(k,mx),ipyx)  &
     + wct2(k,my) * wct2(k,mx) * coltabr (ict2(k,my),ict2(k,mx),ipyx)

   rcy = min(rx(k,my),c1 * colfacr(k) * 10. ** (-tabvaly))

   if (jnmb(mx) >= 5) then
      tabvaln  &
       = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc)  &
       + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc)  &
       + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc)  &
       + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)

      colnum = c1 * colfacc(k) * 10. ** (-tabvaln)
      colnumx = min(cx(k,mx),colnum)
      colnumy = min(cx(k,my),colnum)
      coalnum = min(colnumx,colnumy)
   endif

   rcoal = rcx + rcy
   qrcx = rcx * qx(k,mx)
   qrcy = rcy * qx(k,my)
   qrcoal = qrcx + qrcy
   qcoal = qrcoal / (1.e-20 + rcoal)

   call qtc(qcoal,tcoal,fracliq)

   coalliq = rcoal * fracliq
   coalice = rcoal - coalliq

   if (fracliq .ge. .99) then

      rxfer(k,my,mx) = rxfer(k,my,mx) + rcy
      qrxfer(k,my,mx) = qrxfer(k,my,mx) + qrcy
      if (jnmb(mx) >= 5)  &
         enxfer(k,my,my) = enxfer(k,my,my) + colnumy

      if(imbudget >= 1) ice2rain(k) = ice2rain(k) + rcy * budget_scale
      if(imbudtot >= 1) ice2raint(k) = ice2raint(k) + rcy * budget_scalet

   else

      rfinlz = min(rcoal,  &
         alpha(jhcaty) * coalliq + beta(jhcaty) * rcx)

      xtoz = min(rcx,rfinlz)

      if(imbudget >= 1) rain2ice(k) = rain2ice(k) + rcx * budget_scale
      if(imbudget == 2) then
        if(my==3) rain2pr(k) = rcx * budget_scale
        if(my==4) rain2sn(k) = rcx * budget_scale
        if(my==5) rain2ag(k) = rcx * budget_scale
        if(my==6) rain2gr(k) = rcx * budget_scale
        if(my==7) rain2ha(k) = rcx * budget_scale
        if(my .ne. mz) rain2ha_xtra(k) = rain2ha_xtra(k) + xtoz * budget_scale
      endif

      if(imbudtot >= 1) rain2icet(k) = rain2icet(k) + rcx * budget_scalet
      if(imbudtot == 2) then
        if(my==3) rain2prt(k) = rain2prt(k) + rcx * budget_scalet
        if(my==4) rain2snt(k) = rain2snt(k) + rcx * budget_scalet
        if(my==5) rain2agt(k) = rain2agt(k) + rcx * budget_scalet
        if(my==6) rain2grt(k) = rain2grt(k) + rcx * budget_scalet
        if(my==7) rain2hat(k) = rain2hat(k) + rcx * budget_scalet
        if(my .ne. mz) rain2ha_xtrat(k) = rain2ha_xtrat(k) + xtoz * budget_scalet
      endif

      rxfer(k,mx,mz) = rxfer(k,mx,mz) + xtoz
      rxfer(k,mx,my) = rxfer(k,mx,my) + rcx - xtoz
      if (my .ne. mz) rxfer(k,my,mz) = rxfer(k,my,mz)  &
         + rfinlz - xtoz

! NEED TO USE QCOAL TO TRANSFER Q?

      qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + qx(k,mx) * xtoz
      qrxfer(k,mx,my) = qrxfer(k,mx,my) + qx(k,mx) * (rcx - xtoz)
      if (my .ne. mz) qrxfer(k,my,mz) = qrxfer(k,my,mz)  &
         + qx(k,my) * (rfinlz - xtoz)

      if (jnmb(mx) >= 5) then
         if (my .eq. mz) then
            enxfer(k,mx,mx) = enxfer(k,mx,mx) + colnumx
         elseif (colnumy .ge. colnumx) then
            cfinlz = coalnum * rfinlz / (rcoal + 1.e-20)
            enxfer(k,mx,mz) = enxfer(k,mx,mz) + cfinlz
            enxfer(k,mx,mx) = enxfer(k,mx,mx) + colnumx - cfinlz
            enxfer(k,my,my) = enxfer(k,my,my) + colnumy
         else
            cfinlz = coalnum * rfinlz / (rcoal + 1.e-20)
            enxfer(k,my,mz) = enxfer(k,my,mz) + cfinlz
            enxfer(k,mx,mx) = enxfer(k,mx,mx) + colnumx
            enxfer(k,my,my) = enxfer(k,my,my) + colnumy - cfinlz
         endif
      endif

   endif
endif
enddo

! also include loss of aerosol

return
end subroutine col3

!******************************************************************************

subroutine colxfers(m1,k1,k2,i,j,rloss,enloss)

!use micphys

implicit none

integer :: m1,i,j,k,lcat,kd1,kd2,jcat
integer, dimension(11) :: k1,k2
real, dimension(m1) :: rloss,enloss

!  All rxfer values are nonnegative.

do lcat = 1,8
   if (jnmb(lcat) .ge. 1) then
      kd1 = k1(lcat)
      kd2 = k2(lcat)

      do k = kd1,kd2
         rloss(k) = 0.
         enloss(k) = 0.
      enddo

      do jcat = 1,8
! change this to include enxfer of the same categories
         if (jnmb(jcat) .ge. 1) then
            if (lcat .ne. jcat) then
               do k = kd1,kd2
                  rloss(k) = rloss(k) + rxfer(k,lcat,jcat)
               enddo
            endif
            do k = kd1,kd2
               enloss(k) = enloss(k) + enxfer(k,lcat,jcat)
            enddo
         endif
      enddo

      do k = kd1,kd2
         rloss(k) = min(1.,rx(k,lcat) / max(1.e-20,rloss(k)))
         enloss(k) = min(1.,cx(k,lcat) / max(1.e-10,enloss(k)))
      enddo

      do jcat = 1,8
         if (jnmb(jcat) .ge. 1) then
            if (lcat .ne. jcat) then
               do k = kd1,kd2
                  rxfer(k,lcat,jcat) = rxfer(k,lcat,jcat)*rloss(k)
                  qrxfer(k,lcat,jcat)=qrxfer(k,lcat,jcat)*rloss(k)
               enddo
            endif
            do k = kd1,kd2
               enxfer(k,lcat,jcat) = enxfer(k,lcat,jcat)*enloss(k)
            enddo
         endif
      enddo
   endif
enddo

do lcat = 1,8

   if (jnmb(lcat) .ge. 1) then

      kd1 = k1(lcat)
      kd2 = k2(lcat)

      do jcat = 1,8

         !TRANSFER MIXING RATIO AND NUMBER AND CCN MASS BETWEEN CATEGORIES
         if (jnmb(jcat) .ge. 1 .and. lcat .ne. jcat) then
            do k = kd1,kd2
               if(rxfer(k,lcat,jcat)>rx(k,lcat)) then
                 rxfer(k,lcat,jcat)=rx(k,lcat)
               endif
               !Transfer aerosol-in-hydrometeor masses between categories
               if(iccnlev>=2 .and. rx(k,lcat)>0.0) then
                rxferratio = rxfer(k,lcat,jcat)/rx(k,lcat)
                if(lcat==1.or.lcat==2.or.lcat==3.or.lcat==8) then
                  ccnmass  = cnmhx(k,lcat) * rxferratio
                  cnmhx(k,lcat) = cnmhx(k,lcat) - ccnmass
                  if(jcat==2.or.jcat==3.or.jcat==8) &
                    cnmhx(k,jcat) = cnmhx(k,jcat) + ccnmass
                endif
               endif
               rx(k,lcat) = rx(k,lcat) - rxfer(k,lcat,jcat)
               rx(k,jcat) = rx(k,jcat) + rxfer(k,lcat,jcat)
               qr(k,lcat) = qr(k,lcat) - qrxfer(k,lcat,jcat)
               qr(k,jcat) = qr(k,jcat) + qrxfer(k,lcat,jcat)
               cx(k,lcat) = cx(k,lcat) - enxfer(k,lcat,jcat)
               cx(k,jcat) = cx(k,jcat) + enxfer(k,lcat,jcat)
               if(rx(k,lcat)<0.0) rx(k,lcat)=0.0
               if(cx(k,lcat)<0.0) cx(k,lcat)=0.0
            enddo
         endif
      enddo

      !CHANGE TO ".GE. 5" TO INCLUDE OPTION 6 & 7 FOR CLOUD AND PRIS
      if (jnmb(lcat) >= 5) then
         do k = kd1,kd2
            cx(k,lcat) = cx(k,lcat) - enxfer(k,lcat,lcat)
            if(cx(k,lcat)<0.0) cx(k,lcat)=0.0
         enddo
      endif

   endif
enddo
return
end subroutine colxfers


!###########################################################################
!
!  FOR ICNFLG=3, DEBATING WHETHER TO KEEP IT #/M4 OR CHANGE
!  PARM TO #/KG/M.  NEED TO DEFINE AVMIPSA, ETC. FOR ALL CATEGORIES.
!  MAY WANT TO DEFINE C1 TOO AND RENAME IT.
!  IMPORTANT ISSUE: k loop limits for the jnmb == 5 sections
!  need to consider collection efficiencies for different habits?
!  collection efficiency for hail too high.  big hail should not
!  coallesce.

subroutine each_call(m1,dtlt)

use rconstants
!use micphys

implicit none

integer :: m1,lcat,k,lhcat
real :: dtlt
integer, dimension(8) :: lcat0
data lcat0 /1,2,3,4,5,6,7,16/ ! lcat corressponding to lhcat

! Initialize constants for vapor diffusion and, for fixed diameter cases, emb.

colf = .785 * dtlt
pi4dt = pi4 * dtlt
sl(1) = alvl
sl(2) = alvi
sc(1) = 4186.
sc(2) = 2093.
sj(1) = 0
sj(2) = 1
sj(3) = 0
sj(4) = 0
sj(5) = 0
sj(6) = 1
sj(7) = 1
sj(8) = 0
sk(1) = alli
sk(2) = 0.

do lcat = 1,8
   lhcat = lcat0(lcat)
   if (jnmb(lcat) == 2) then
      do k = 2,m1-1
         emb(k,lcat) = cfmas(lhcat) * parm(lcat) ** pwmas(lhcat)
      enddo
   endif
   do k = 2,m1-1
      jhcat(k,lcat) = lhcat
   enddo
enddo

do k = 2,m1-1
   sh(k,1) = 0.
   sh(k,2) = 1.
   sh(k,6) = 1.
   sh(k,7) = 1.
   sh(k,8) = 0.

   sm(k,1) = 1.
   sm(k,2) = 1.
   sm(k,8) = 1.
enddo

return
end subroutine each_call

!******************************************************************************

subroutine range_check(m1,k1,k2,k3,i,j,lpw_R,micro)

use mem_micro
!use micphys

implicit none

type (micro_vars) :: micro

real :: lpw_R
integer :: m1,i,j,k,lcatt,lcat,l,jcat,lpw
integer, dimension(11) :: k1,k2,k3

! zero out microphysics scratch arrays for the present i,j column
lpw=int(lpw_R)

do lcat = 1,ncat
   do k = 2,m1-1
      rx(k,lcat) = 0.
      cx(k,lcat) = 0.
      qr(k,lcat) = 0.
      qx(k,lcat) = 0.
      vap(k,lcat) = 0.
      tx(k,lcat) = 0.
   enddo

   if (jnmb(lcat) >= 3) then
      do k = 2,m1-1
         emb(k,lcat) = 0.
      enddo
   endif

   do jcat = 1,ncat
      do k = 2,m1-1
         rxfer(k,lcat,jcat) = 0.
         qrxfer(k,lcat,jcat) = 0.
         enxfer(k,lcat,jcat) = 0.
      enddo
   enddo
enddo

do l = 1,8
   k1(l) = lpw
   k2(l) = 1
enddo

! fill scratch arrays for dust and salt modes
do k = lpw,m1-1
  if (idust == 1 .or. imd1flg == 1) md1nx(k) = micro%md1np(k,i,j)
  if (idust == 1 .or. imd2flg == 1) md2nx(k) = micro%md2np(k,i,j)
  if (isalt == 1) saltfx(k) = micro%salt_filmp(k,i,j)
  if (isalt == 1) saltjx(k) = micro%salt_jetp(k,i,j)
  if (isalt == 1) saltsx(k) = micro%salt_spmp(k,i,j)
enddo

! fill scratch arrays for cloud water

if (jnmb(1) >= 1) then
   do k = lpw,m1-1
      if (micro%rcp(k,i,j) >= 1.e-12) then
         k2(1) = k
         rx(k,1) = micro%rcp(k,i,j)
         if (jnmb(1) >= 5) cx(k,1) = micro%ccp(k,i,j)
         if (iccnlev >= 2) cnmhx(k,1) = micro%cnm1p(k,i,j)
      else
         if (k2(1) == 1) k1(1) = k + 1
      endif
      if (jnmb(1) >= 5) then
        cccnx(k) = micro%cccnp(k,i,j)
        cccmx(k) = micro%cccmp(k,i,j)
      endif
   enddo
endif

! fill scratch arrays for rain

if (jnmb(2) >= 1) then
   do k = lpw,m1-1
      if (micro%rrp(k,i,j) >= 1.e-12) then
         k2(2) = k
         rx(k,2) = micro%rrp(k,i,j)
         qx(k,2) = micro%q2(k,i,j)
         qr(k,2) = qx(k,2) * rx(k,2)
         if (jnmb(2) >= 5) cx(k,2) = micro%crp(k,i,j)
         if (iccnlev >= 2) cnmhx(k,2) = micro%cnm2p(k,i,j)
      else
         if (k2(2) == 1) k1(2) = k + 1
      endif
   enddo
endif

! fill scratch arrays for pristine ice

if (jnmb(3) >= 1) then
   do k = lpw,m1-1
      if (micro%rpp(k,i,j) >= 1.e-12) then
         k2(3) = k
         rx(k,3) = micro%rpp(k,i,j)
         cx(k,3) = micro%cpp(k,i,j)
         if (iccnlev >= 2) cnmhx(k,3) = micro%cnm3p(k,i,j)
      else
         if (k2(3) == 1) k1(3) = k + 1
      endif
      if (jnmb(3) >= 5) cifnx(k) = micro%cifnp(k,i,j)
   enddo
endif

! fill scratch arrays for snow

if (jnmb(4) >= 1) then
   do k = lpw,m1-1
      if (micro%rsp(k,i,j) >= 1.e-12) then
         k2(4) = k
         rx(k,4) = micro%rsp(k,i,j)
         if (jnmb(4) >= 5) cx(k,4) = micro%csp(k,i,j)
      else
         if (k2(4) == 1) k1(4) = k + 1
      endif
   enddo
endif

! fill scratch arrays for aggregates

if (jnmb(5) >= 1) then
   do k = lpw,m1-1
      if (micro%rap(k,i,j) >= 1.e-12) then
         k2(5) = k
         rx(k,5) = micro%rap(k,i,j)
         if (jnmb(5) >= 5) cx(k,5) = micro%cap(k,i,j)
      else
         if (k2(5) == 1) k1(5) = k + 1
      endif
   enddo
endif

! fill scratch arrays for graupel

if (jnmb(6) >= 1) then
   do k = lpw,m1-1
      if (micro%rgp(k,i,j) >= 1.e-12) then
         k2(6) = k
         rx(k,6) = micro%rgp(k,i,j)
         qx(k,6) = micro%q6(k,i,j)
         qr(k,6) = qx(k,6) * rx(k,6)
         if (jnmb(6) >= 5) cx(k,6) = micro%cgp(k,i,j)
      else
         if (k2(6) == 1) k1(6) = k + 1
      endif
   enddo
endif

! fill scratch arrays for hail

if (jnmb(7) >= 1) then
   do k = lpw,m1-1
      if (micro%rhp(k,i,j) >= 1.e-12) then
         k2(7) = k
         rx(k,7) = micro%rhp(k,i,j)
         qx(k,7) = micro%q7(k,i,j)
         qr(k,7) = qx(k,7) * rx(k,7)
         if (jnmb(7) >= 5) cx(k,7) = micro%chp(k,i,j)
      else
         if (k2(7) == 1) k1(7) = k + 1
      endif
   enddo
endif

! fill scratch arrays for drizzle
if (jnmb(8) >= 1) then
   do k = lpw,m1-1
      if (micro%rdp(k,i,j) >= 1.e-12) then
         k2(8) = k
         rx(k,8) = micro%rdp(k,i,j)
         if (jnmb(8) >= 5) cx(k,8) = micro%cdp(k,i,j)
         if (iccnlev >= 2) cnmhx(k,8) = micro%cnm8p(k,i,j)
      else
         if (k2(8) == 1) k1(8) = k + 1
      endif
      if (jnmb(8) >= 5) then
        gccnx(k) = micro%gccnp(k,i,j)
        gccmx(k) = micro%gccmp(k,i,j)
      endif
   enddo
endif

k3(1) = k2(1)
k3(3) = k2(3)
k3(8) = k2(8)

k1(9) = min(k1(1),k1(2),k1(8))
k2(9) = max(k2(1),k2(2),k2(8))
k1(10) = min(k1(3),k1(4),k1(5),k1(6),k1(7))
k2(10) = max(k2(3),k2(4),k2(5),k2(6),k2(7))
k1(11) = min(k1(9),k1(10))
k2(11) = max(k2(9),k2(10))

return
end subroutine range_check

!******************************************************************************

subroutine each_column(m1,k1,k2,i,j,lpw,rv,dn0)

use rconstants
!use micphys

implicit none

integer :: m1,i,j,k,nt,ns,lpw
integer, dimension(11) :: k1,k2
real :: ck1,ck2,ck3,elsref,elsrefp,dplinv,eisref,eisrefp,dpiinv,relhum
real, dimension(m1) :: rv,dn0
real :: rslf,rsif,eslf,eslpf,esif,esipf

data ck1,ck2,ck3/-4.818544e-3,1.407892e-4,-1.249986e-7/

do k = lpw,m1-1
   rvlsair(k) = rslf (press(k),tair(k))
   rvisair(k) = rsif (press(k),tair(k))
   dn0i(k) = 1. / dn0(k)
   tairc(k)   = tair(k) - 273.15
   tx(k,1) = tairc(k)
   thrmcon(k) = ck1 + (ck2 + ck3 * tair(k)) * tair(k)
   dynvisc(k) = .1718e-4 + .49e-7 * tairc(k)

   ! Diagnose habit of pristine ice and snow

   nt = max(1,min(31,-nint(tairc(k))))
   relhum = min(1.,rv(k) / rvlsair(k))
   ns = max(1,nint(100. * relhum))
   jhcat(k,3) = jhabtab(nt,ns,1)
   jhcat(k,4) = jhabtab(nt,ns,2)

enddo

do k = k1(11),k2(11)
   vapdif(k)     = 2.14 * (tair(k) / 273.15) ** 1.94 / press(k)
   rdynvsci(k) = sqrt(1. / dynvisc(k))
   denfac(k) = sqrt(dn0i(k))

   colfacr(k) = colf * denfac(k) * dn0(k)
   colfacr2(k) = 2. * colfacr(k)

!Saleeby(2010): Loftus: remove density from: colfacc(k) = colfacr(k) * dn0(k)
   colfacc(k) = colfacr(k)

   colfacc2(k) = 2. * colfacc(k)

   tref(k,1)   = tairc(k) - min(25.,700. * (rvlsair(k) - rv(k)))
   sa(k,2) = thrmcon(k) * sa(k,1)
   sa(k,3) = thrmcon(k) * (tairstrc(k) + sa(k,1) * rvstr(k))

   sumuy(k) = 0.
   sumuz(k) = 0.
   sumvr(k) = 0.
enddo

do k = k1(9),k2(9)
   elsref       = eslf(tref(k,1))
   elsrefp      = eslpf(tref(k,1))
   dplinv       = 1. / (press(k) - elsref)
   rvsref (k,1) = .622 * elsref * dplinv
   rvsrefp(k,1) = .622 * elsrefp * dplinv * (1. + elsref * dplinv)

   sa(k,4) = rvsrefp(k,1) * tref(k,1) - rvsref(k,1)
   sa(k,6) = alvl * rvsrefp(k,1)
   sa(k,8) = alvl * sa(k,4)
enddo

do k = k1(10),k2(10)
   tref(k,2)    = min(0.,tref(k,1))
   eisref       = esif (tref(k,2))
   eisrefp      = esipf(tref(k,2))
   dpiinv       = 1. / (press(k) - eisref)
   rvsref (k,2) = .622 * eisref * dpiinv
   rvsrefp(k,2) = .622 * eisrefp * dpiinv * (1. + eisref * dpiinv)
   rvs0(k)      = 379.4 / (press(k) - 610.)

   sa(k,5) = rvsrefp(k,2) * tref(k,2) - rvsref(k,2)
   sa(k,7) = alvi * rvsrefp(k,2)
   sa(k,9) = alvi * sa(k,5)
   sh(k,3) = 0.
   sh(k,4) = 0.
   sh(k,5) = 0.

enddo

return
end subroutine each_column

!******************************************************************************

subroutine enemb(m1,k1,k2,lcat,dn0,i,j)

!use micphys

implicit none

integer :: m1,k1,k2,lcat,i,j,k,lhcat
real :: embi,parmi,embtemp
real, dimension(m1) :: dn0

if (jnmb(lcat) == 2) then
   embi = 1. / emb(2,lcat)
   do k = k1,k2
      cx(k,lcat) = rx(k,lcat) * embi
   enddo
elseif (jnmb(lcat) == 3) then
   do k = k1,k2
      lhcat = jhcat(k,lcat)
      emb(k,lcat) = cfemb0(lhcat) * (dn0(k) * rx(k,lcat)) ** pwemb0(lhcat)
      cx(k,lcat) = cfen0(lhcat) * dn0i(k)  &
         * (dn0(k) * rx(k,lcat)) ** pwen0(lhcat)
   enddo
elseif (jnmb(lcat) == 4) then
   parmi = 1. / parm(lcat)
   do k = k1,k2
      emb(k,lcat) = max(emb0(lcat),min(emb1(lcat),rx(k,lcat) * parmi))
      cx(k,lcat) = rx(k,lcat) / emb(k,lcat)
   enddo
elseif (jnmb(lcat) >= 5) then
   do k = k1,k2
      embtemp=rx(k,lcat)/max(1.e-12,cx(k,lcat))
      emb(k,lcat) = max(emb0(lcat),min(emb1(lcat),rx(k,lcat)  &
         / max(1.e-12,cx(k,lcat))))
      !Saleeby(2011): Use of single precision here allows for enemb to produce
      !artificial number concentration even when mass is in bounds. This is
      !generally small but could accumulate over time. Added the IF statement
      !to stop adjustment if mean size is in bounds.
      if( (embtemp > 0.0 .and. (embtemp < 0.999*emb0(lcat)   .or.  &
                                embtemp > 1.001*emb1(lcat))) .or.  &
          (rx(k,lcat) == 0.0 .and. cx(k,lcat) > 0.0 )        .or.  &
          (embtemp >= emb0(lcat) .and. embtemp <= emb1(lcat) .and. &
               cx(k,lcat)<1.e-12 .and. rx(k,lcat)>0.0))      then
        cx(k,lcat) = rx(k,lcat) / emb(k,lcat)
      endif
   enddo
endif

return
end subroutine enemb

!******************************************************************************

subroutine x02(m1,k1,k2,lcat,dn0,i,j,meltice,meltpris,meltsnow,meltaggr &
              ,meltgrau,melthail,melticet,meltprist,meltsnowt,meltaggrt &
              ,meltgraut,melthailt)

use rconstants
!use micphys

implicit none

integer :: m1,lcat,i,j,k,lhcat,inc,idns
integer, dimension(11) :: k1,k2
real :: rinv,closs,rxinv,rmelt,fracliq,cmelt,tcoal,ricetor6,rshed,rmltshed  &
       ,qrmltshed,shedmass,fracmloss,dn
real, dimension(m1) :: dn0,meltice,meltpris,meltsnow,meltaggr,meltgrau,melthail &
                    ,melticet,meltprist,meltsnowt,meltaggrt,meltgraut,melthailt

k1(lcat) = k1(11)
k2(lcat) = 1
do k = k1(11),k2(11)
   if (rx(k,lcat) >= 1.e-12) k2(lcat) = k
   if (k2(lcat) == 1 .and. rx(k,lcat) < 1.e-12) k1(lcat) = k + 1
enddo

if ((lcat == 2 .or. lcat >= 4) .and. (lcat .ne. 8)) then

   call enemb(m1,k1(lcat),k2(lcat),lcat,dn0,i,j)

   do k = k1(lcat),k2(lcat)

      if (rx(k,lcat) >= 1.e-12) then

      lhcat = jhcat(k,lcat)
      vterm(k,lcat) = -vtfac(lhcat) * emb(k,lcat) ** pwvtmasi(lhcat) &
         * denfac(k)

      endif

   enddo
endif

if (lcat == 2) then

   do k = k1(lcat),k2(lcat)

      if (rx(k,lcat) >= 1.e-12) then

      rxinv = 1. / rx(k,lcat)
      qx(k,lcat) = qr(k,lcat) * rxinv
! limit rain to under 48C and over -80C
      qx(k,lcat) = max(0.,min(1.6*alli,qx(k,lcat)))

      endif

   enddo

elseif (lcat == 3) then
!Allow pristine ice to melt to cloud1 since we assume that smaller
! particles will melt first. Perhaps need a way in the future to treat
! melting of larger pristine ice into cloud2.

   do k = k1(lcat),k2(lcat)

      if (rx(k,lcat) >= 1.e-12) then

      rinv = 1. / rx(k,lcat)
      qx(k,lcat) = qr(k,lcat) * rinv

      call qtc(qx(k,lcat),tcoal,fracliq)

      rmelt = rx(k,lcat) * fracliq
      cmelt = cx(k,lcat) * fracliq

      rx(k,lcat) = rx(k,lcat) - rmelt
      rx(k,1) = rx(k,1) + rmelt
      cx(k,lcat) = cx(k,lcat) - cmelt
      cx(k,1) = cx(k,1) + cmelt

      if(imbudget >= 1) meltice(k)  = meltice(k) + rmelt * budget_scale
      if(imbudget == 2) meltpris(k) = rmelt * budget_scale

      if(imbudtot >= 1) melticet(k)  = melticet(k)  + rmelt * budget_scalet
      if(imbudtot == 2) meltprist(k) = meltprist(k) + rmelt * budget_scalet

      endif

   enddo
!
! meyers - source for cloud aerosol number here?
!
elseif (lcat == 4 .or. lcat == 5) then
!Allow snow and aggregates to melt to graupel.
!Perhaps snow should melt to cloud2 and aggregates to rain??
!Change this??? move to rain instead ??? look at melting decisions in col2

   do k = k1(lcat),k2(lcat)

     if (rx(k,lcat) >= 1.e-12) then

      rinv = 1. / rx(k,lcat)
      qx(k,lcat) = qr(k,lcat) * rinv
      call qtc(qx(k,lcat),tcoal,fracliq)

      if (fracliq > 1.e-6) then
         rmelt = rx(k,lcat) * fracliq
         ricetor6 = min(rx(k,lcat) - rmelt,rmelt)
         rx(k,lcat) = rx(k,lcat) - rmelt - ricetor6
         rx(k,6) = rx(k,6) + rmelt + ricetor6
         qr(k,6) = qr(k,6) + rmelt * alli
         qx(k,lcat) = 0.

! keep the above the same with ricetor6
! meyers - use sa melt table here? yes
         fracmloss = (rmelt + ricetor6) * rinv
         closs = enmlttab(int(200. * fracmloss) + 1,jhcat(k,lcat)) * cx(k,lcat)
         cx(k,lcat) = cx(k,lcat) - closs
         cx(k,6) = cx(k,6) + closs

         if(imbudget >= 1) meltice(k) = meltice(k)+(rmelt+ricetor6)*budget_scale
         if(imbudget == 2) then
           if(lcat==4) meltsnow(k) = (rmelt + ricetor6) * budget_scale
           if(lcat==5) meltaggr(k) = (rmelt + ricetor6) * budget_scale
         endif

         if(imbudtot >= 1) melticet(k) = melticet(k)+(rmelt+ricetor6)*budget_scalet
         if(imbudtot == 2) then
           if(lcat==4) meltsnowt(k) = meltsnowt(k) + (rmelt + ricetor6) * budget_scalet
           if(lcat==5) meltaggrt(k) = meltaggrt(k) + (rmelt + ricetor6) * budget_scalet
         endif

      endif

     endif
   enddo


elseif (lcat == 6) then

   do k = k1(lcat),k2(lcat)

      if (rx(k,lcat) >= 1.e-12) then

      rxinv = 1. / rx(k,lcat)
      qx(k,lcat) = qr(k,lcat) * rxinv
      call qtc(qx(k,lcat),tcoal,fracliq)

      if (fracliq > 0.95) then
         rx(k,2) = rx(k,2) + rx(k,6)
         qr(k,2) = qr(k,2) + rx(k,6) * alli
         cx(k,2) = cx(k,2) + cx(k,6)

         if(imbudget >= 1) meltice(k)  = meltice(k) + rx(k,6) * budget_scale
         if(imbudget == 2) meltgrau(k) = rx(k,6) * budget_scale

         if(imbudtot >= 1) melticet(k)  = melticet(k) + rx(k,6) * budget_scalet
         if(imbudtot == 2) meltgraut(k) = meltgraut(k) + rx(k,6) * budget_scalet

         rx(k,6) = 0.
         qr(k,6) = 0.
         qx(k,6) = 0.
         cx(k,6) = 0.
      endif

      endif

   enddo

elseif (lcat == 7) then

   shedmass = 5.236e-7
   do k = k1(lcat),k2(lcat)

      if (rx(k,lcat) >= 1.e-12) then

      rxinv = 1. / rx(k,lcat)
      qx(k,lcat) = qr(k,lcat) * rxinv
      call qtc(qx(k,lcat),tcoal,fracliq)

      if (fracliq > 0.95) then
         rx(k,2) = rx(k,2) + rx(k,7)
         qr(k,2) = qr(k,2) + rx(k,7) * alli
         cx(k,2) = cx(k,2) + cx(k,7)

         if(imbudget >= 1) meltice(k)  = meltice(k) + rx(k,7) * budget_scale
         if(imbudget == 2) melthail(k) = rx(k,7) * budget_scale

         if(imbudtot >= 1) melticet(k)  = melticet(k) + rx(k,7) * budget_scalet
         if(imbudtot == 2) melthailt(k) = melthailt(k) + rx(k,7) * budget_scalet

         rx(k,7) = 0.
         qr(k,7) = 0.
         qx(k,7) = 0.
         cx(k,7) = 0.

!  take out following IF statement?

      elseif (fracliq > 0.3) then

         lhcat = jhcat(k,lcat)
         inc = nint(200. * fracliq) + 1
         dn = dnfac(lhcat) * emb(k,lcat) ** pwmasi(lhcat)
         idns = max(1,nint(1.e3 * dn * gnu(lcat)))
         rshed = rx(k,lcat) * shedtab(inc,idns)
         rmltshed = rshed
         qrmltshed = rmltshed * alli

         rx(k,2) = rx(k,2) + rmltshed
         qr(k,2) = qr(k,2) + qrmltshed
         cx(k,2) = cx(k,2) + rshed / shedmass

         rx(k,lcat) = rx(k,lcat) - rmltshed
         qr(k,lcat) = qr(k,lcat) - qrmltshed
         qx(k,lcat) = qr(k,lcat) * (1./rx(k,lcat))

         if(imbudget >= 1) meltice(k) = meltice(k) + rmltshed * budget_scale
         if(imbudget == 2) melthail(k) = rmltshed * budget_scale

         if(imbudtot >= 1) melticet(k)  = melticet(k) + rmltshed * budget_scalet
         if(imbudtot == 2) melthailt(k) = melthailt(k) + rmltshed * budget_scalet

      endif

      endif

   enddo

endif
return
end subroutine x02

!******************************************************************************

subroutine c03(m1,k1,k2,lcat,dn0,i,j)

!use micphys

implicit none

integer :: m1,k1,k2,lcat,i,j
real, dimension(m1) :: dn0

if (jnmb(lcat) >= 3) call enemb(m1,k1,k2,lcat,dn0,i,j)

return
end subroutine c03

!******************************************************************************

subroutine pc03(m1,k1,k2,lcat,dn0,i,j)

!use micphys

implicit none

integer :: m1,k1,k2,lcat,i,j,k,lhcat
real, dimension(m1) :: dn0

if (jnmb(lcat) >= 3) call enemb(m1,k1,k2,lcat,dn0,i,j)

do k = k1,k2

   if (rx(k,lcat) >= 1.e-12) then

   lhcat = jhcat(k,lcat)
   vterm(k,lcat) = -vtfac(lhcat) * emb(k,lcat) ** pwvtmasi(lhcat) * denfac(k)

   endif

enddo

return
end subroutine pc03

!******************************************************************************

subroutine sedim(m1,lcat,ngr,nembfall,maxkfall,k1,k2,lpw,i,j  &
   ,rtp,thp,theta,dn0,alphasfc  &
   ,pcpg,qpcpg,dpcpg,dtlti,cnew,rnew,qrnew  &
   ,pcpfillc,pcpfillr,sfcpcp,allpcp,dzt,if_adap)

use rconstants
!use micphys

implicit none

integer :: m1,lcat,ngr,nembfall,maxkfall, if_adap  &
          ,k1,k2,lpw,i,j,k,lhcat,iemb,iemb2,kkf,kk
real :: colddn0,rolddn0,qrolddn0,dispemb,riemb,wt2,psfc,qpcpg,pcpg,dpcpg  &
   ,dtlti,qnew,alphasfc,cnmhddn0
real, dimension(m1) :: rtp,thp,theta,dn0,cnew,rnew,qrnew,dzt,cnmhpnew
real, dimension(m1,maxkfall,nembfall,nhcat) :: pcpfillc,pcpfillr
real, dimension(maxkfall,nembfall,nhcat) :: sfcpcp
real, dimension(m1,nembfall,nhcat) :: allpcp

pcprx(lcat) = 0.
do k = 1,m1
   rnew(k) = 0.
   cnew(k) = 0.
   qrnew(k) = 0.
   cnmhpnew(k) = 0.
   pcpvx(k,lcat) = 0.
enddo

do k = k1,k2
   lhcat = jhcat(k,lcat)

   if (rx(k,lcat) > 1.e-20) then
      colddn0 = cx(k,lcat) * dn0(k)
      rolddn0 = rx(k,lcat) * dn0(k)
      qrolddn0 = qx(k,lcat) * rolddn0
      if(iccnlev>=2.and.(lcat==2.or.lcat==3.or.lcat==8)) cnmhddn0 = cnmhx(k,lcat) * dn0(k)

      !Here determine which set of powerlaws to use: the original
      ! ones in RAMS or the ones from R.Carver adapted from Mitchell 1996.
      !The Mitchell power laws are not based at sea level so we adjust the
      ! density factor based at 0.7 kg/m3 instead of 1.0 kg/m3.
      if(iplaws==0) then
        dispemb = ch1(lhcat)  &
          * (emb(k,lcat)/cfmas(lhcat)) ** ch3(lhcat) * sqrt(dn0i(k))
      else
        dispemb = ch1(lhcat)  &
          * (emb(k,lcat)/cfmas(lhcat)) ** ch3(lhcat) * (0.7*dn0i(k))**.362
      endif

      riemb = 1. + ch2(lhcat,ngr) * log10(dispemb / dispemb0(lhcat,ngr))

      !Limiting iemb to max of nembfall
      iemb = min(nint(riemb),nembfall)

      if (k <= maxkfall) then
         psfc = rolddn0 * sfcpcp(k,iemb,lhcat)
      endif

      do kkf = 1,min(maxkfall,k-1)
         kk = k + 1 - kkf

         cnew(kk)  = cnew(kk)   &
            +  colddn0 * dn0i(kk) * pcpfillc(k,kkf,iemb,lhcat)
         rnew(kk)  = rnew(kk)   &
            +  rolddn0 * dn0i(kk) * pcpfillr(k,kkf,iemb,lhcat)
         qrnew(kk) = qrnew(kk)  &
            + qrolddn0 * dn0i(kk) * pcpfillr(k,kkf,iemb,lhcat)
         if(iccnlev>=2.and.(lcat==2.or.lcat==3.or.lcat==8)) then
           cnmhpnew(kk) = cnmhpnew(kk) &
            + cnmhddn0 * dn0i(kk) * pcpfillr(k,kkf,iemb,lhcat)
         endif
      enddo

      !Surface Precip rate
      if (k <= maxkfall) then
         qpcpg = qpcpg + psfc * qx(k,lcat)
         pcprx(lcat) = pcprx(lcat) + psfc
      endif

      !Precip rate at all levels (mm/s)
      pcpvx(k,lcat)= rolddn0 * allpcp(k,iemb,lhcat) * dtlti

   endif
enddo

if (if_adap == 1) then
   do k = 2,lpw-1
      pcprx(lcat) = pcprx(lcat) + rnew(k) * dn0(k) / dzt(k)
      qpcpg = qpcpg + qrnew(k) * dn0(k) / dzt(k)
      cnew(k) = 0.
      rnew(k) = 0.
      qrnew(k) = 0.
   enddo
endif

pcpg = pcpg + pcprx(lcat)
accpx(lcat) = pcprx(lcat)
dpcpg = dpcpg + pcprx(lcat) * alphasfc
pcprx(lcat) = pcprx(lcat) * dtlti

do k = lpw,k2
   rtp(k) = rtp(k) + rnew(k) - rx(k,lcat)
   qnew = qrnew(k) / max(1.e-20, rnew(k))

!         if (iqflag == 1) then
      tairc(k) = tairc(k) - thp(k) * thp(k)  &
         * (2820. * (rnew(k) - rx(k,lcat))  &
         - cpi * (qrnew(k) - qx(k,lcat) * rx(k,lcat)))  &
         / (max(tair(k), 253.) * theta(k))
!         else
!            tairc(k) = tairc(k) - thp(k) * thp(k) * 2820.
!     +          * (rnew(k) - rx(k,lcat)) / (max(tair(k), 253.) * theta(k))
!         endif

   if(iccnlev>=2.and.(lcat==2.or.lcat==3.or.lcat==8)) cnmhx(k,lcat) = cnmhpnew(k)
   rx(k,lcat) = rnew(k)
   cx(k,lcat) = cnew(k)
   qx(k,lcat) = qnew

   if (rx(k,lcat) < 1.e-12) then
      rx(k,lcat) = 0.
      cx(k,lcat) = 0.
      qx(k,lcat) = 0.
      pcpvx(k,lcat) = 0.
   endif

enddo
return
end subroutine sedim

!******************************************************************************

subroutine negadj1_2M_rams60(m1,m2,m3)

  use mem_basic, only: basic_g ! INTENT(OUT)

  use mem_micro, only: micro_g ! INTENT(OUT)

  use mem_grid,  only: &
       grid_g,         & ! INTENT(IN)
       ngrid             ! INTENT(IN)

 !!use micphys, only:   level ! INTENT(IN)

  use mem_scratch, only: vctr9 ! INTENT(OUT)

implicit none

integer :: m1,m2,m3

if (level == 0) return

call adj1_2M_rams60(m1,m2,m3,grid_g(ngrid)%lpw(1:m1,1),basic_g(ngrid)%rtp(1:m1,1,1)  &
   ,basic_g(ngrid)%thp(1:m1,1,1),micro_g(ngrid),vctr9)

return
end subroutine negadj1_2M_rams60

!******************************************************************************

subroutine adj1_2M_rams60(m1,m2,m3,lpw_R,rtp,thp,micro,vctr9)

use mem_micro
!use micphys

implicit none

integer :: m1,m2,m3
real, dimension(m2,m3) :: lpw_R

type (micro_vars) :: micro

integer :: i,j,k,lcat,ka
real :: frac
real, dimension(m1,m2,m3) :: rtp,thp

real, dimension(m1) :: vctr9
integer, dimension(m2,m3) :: lpw

lpw=int(lpw_R)

if (level .eq. 0) return

do lcat = 1,ncat
   do k = 1,m1
      rx(k,lcat) = 0.
      cx(k,lcat) = 0.
   enddo
enddo

do j = 1,m3
   do i = 1,m2

      !!!ka = lpw(i,j)
      ! Do this for all levels, regardless of ADAP
      ka = 1

      if (jnmb(1) > 0) then
        call ae1kmic(ka,m1,rx(1:m1,1),micro%rcp(1:m1,i,j))
        if(iccnlev >= 2) call ae1kmic(ka,m1,cnmhx(1:m1,1),micro%cnm1p(1:m1,i,j))
      endif
      if (jnmb(2) > 0) then
        call ae1kmic(ka,m1,rx(1:m1,2),micro%rrp(1:m1,i,j))
        if(iccnlev >= 2) call ae1kmic(ka,m1,cnmhx(1:m1,2),micro%cnm2p(1:m1,i,j))
      endif
      if (jnmb(3) > 0) then
        call ae1kmic(ka,m1,rx(1:m1,3),micro%rpp(1:m1,i,j))
        if(iccnlev >= 2) call ae1kmic(ka,m1,cnmhx(1:m1,3),micro%cnm3p(1:m1,i,j))
      endif
      if (jnmb(4) > 0) call ae1kmic(ka,m1,rx(1:m1,4),micro%rsp(1:m1,i,j))
      if (jnmb(5) > 0) call ae1kmic(ka,m1,rx(1:m1,5),micro%rap(1:m1,i,j))
      if (jnmb(6) > 0) call ae1kmic(ka,m1,rx(1:m1,6),micro%rgp(1:m1,i,j))
      if (jnmb(7) > 0) call ae1kmic(ka,m1,rx(1:m1,7),micro%rhp(1:m1,i,j))
      if (jnmb(8) > 0) then
        call ae1kmic(ka,m1,rx(1:m1,8),micro%rdp(1:m1,i,j))
        if(iccnlev >= 2) call ae1kmic(ka,m1,cnmhx(1:m1,8),micro%cnm8p(1:m1,i,j))
      endif

      do lcat = 1,ncat
         do k = ka,m1
            if (rx(k,lcat) < 1.e-12) then
              rx(k,lcat) = 0.
              cnmhx(k,lcat) = 0.
            endif
         enddo
      enddo

      do k = ka,m1
         rtp(k,i,j) = max(0.,rtp(k,i,j))
         vctr9(k) = 1.001 * (rx(k,1)+ rx(k,2) + rx(k,3)  &
            + rx(k,4) + rx(k,5) + rx(k,6) + rx(k,7) + rx(k,8))
      enddo

      do k = ka,m1
         if (vctr9(k) > rtp(k,i,j)) then
            frac = rtp(k,i,j) / (1.e-12 + vctr9(k))
            do lcat = 1,ncat
               rx(k,lcat) = rx(k,lcat) * frac
               cnmhx(k,lcat) = cnmhx(k,lcat) * frac
            enddo
         endif
      enddo

      if (jnmb(1) > 0) then
        call ae1kmic(ka,m1,micro%rcp(1:m1,i,j),rx(1:m1,1))
        if(iccnlev >= 2) call ae1kmic(ka,m1,micro%cnm1p(1:m1,i,j),cnmhx(1:m1,1))
      endif
      if (jnmb(2) > 0) then
        call ae1kmic(ka,m1,micro%rrp(1:m1,i,j),rx(1:m1,2))
        if(iccnlev >= 2) call ae1kmic(ka,m1,micro%cnm2p(1:m1,i,j),cnmhx(1:m1,2))
      endif
      if (jnmb(3) > 0) then
        call ae1kmic(ka,m1,micro%rpp(1:m1,i,j),rx(1:m1,3))
        if(iccnlev >= 2) call ae1kmic(ka,m1,micro%cnm3p(1:m1,i,j),cnmhx(1:m1,3))
      endif
      if (jnmb(4) > 0) call ae1kmic(ka,m1,micro%rsp(1:m1,i,j),rx(1:m1,4))
      if (jnmb(5) > 0) call ae1kmic(ka,m1,micro%rap(1:m1,i,j),rx(1:m1,5))
      if (jnmb(6) > 0) call ae1kmic(ka,m1,micro%rgp(1:m1,i,j),rx(1:m1,6))
      if (jnmb(7) > 0) call ae1kmic(ka,m1,micro%rhp(1:m1,i,j),rx(1:m1,7))
      if (jnmb(8) > 0) then
        call ae1kmic(ka,m1,micro%rdp(1:m1,i,j),rx(1:m1,8))
        if(iccnlev >= 2) call ae1kmic(ka,m1,micro%cnm8p(1:m1,i,j),cnmhx(1:m1,8))
      endif

      if (jnmb(1) >= 5)  &
         call ae1mic(ka,m1,micro%ccp(1:m1,i,j),micro%rcp(1:m1,i,j),rx(1,1))
      if (jnmb(2) >= 5)  &
         call ae1mic(ka,m1,micro%crp(1:m1,i,j),micro%rrp(1:m1,i,j),rx(1,2))
      if (jnmb(3) >= 5)  &
         call ae1mic(ka,m1,micro%cpp(1:m1,i,j),micro%rpp(1:m1,i,j),rx(1,3))
      if (jnmb(4) >= 5)  &
         call ae1mic(ka,m1,micro%csp(1:m1,i,j),micro%rsp(1:m1,i,j),rx(1,4))
      if (jnmb(5) >= 5)  &
         call ae1mic(ka,m1,micro%cap(1:m1,i,j),micro%rap(1:m1,i,j),rx(1,5))
      if (jnmb(6) >= 5)  &
         call ae1mic(ka,m1,micro%cgp(1:m1,i,j),micro%rgp(1:m1,i,j),rx(1,6))
      if (jnmb(7) >= 5)  &
         call ae1mic(ka,m1,micro%chp(1:m1,i,j),micro%rhp(1:m1,i,j),rx(1,7))
      if (jnmb(8) >= 5)  &
         call ae1mic(ka,m1,micro%cdp(1:m1,i,j),micro%rdp(1:m1,i,j),rx(1,8))

      do k = 1,m1
         if(imd1flg==1 .or. idust==1) then
           if(micro%md1np(k,i,j) .lt. 0.) micro%md1np(k,i,j) = 0.0
         endif
         if(imd2flg==1 .or. idust==1) then
           if(micro%md2np(k,i,j) .lt. 0.) micro%md2np(k,i,j) = 0.0
         endif
      enddo

!
!  Think about how thp should change here - should it be due to a change in
!     rtp or to a change in the condensate?
!
!               vctr10(k) = rrp(k,i,j +rpp(k,i,j)  + rsp(k,i,j) + rap(k,i,j)
!     +                   + rgp(k,i,j) + rhp(k,i,j)
!               thp(k,i,j) = thp(k,i,j)
!     +                    * (1. - aklv * (vctr8(k) - rtp(k,i,j))
!c or +                    * (1. - aklv * (vctr10(k) - vctr9(k,i,j))
!     +                    /(max(temp, 253.)))

   enddo
enddo
return
end subroutine adj1_2M_rams60

!---------------------------------------------------------------------------

subroutine ae1mic(ka,m1,c3,r3,r1)

implicit none

integer :: m1,ka
real, dimension(m1) :: c3,r3,r1
integer :: k

do k = ka,m1
   c3(k) = c3(k) * r1(k) / (1.e-12 + r3(k))
   if (c3(k) < 0.) c3(k) = 0.
enddo

return
end subroutine ae1mic

!---------------------------------------------------------------------------

subroutine ae1kmic(ka,kb,cr3,cr1)

implicit none

integer :: ka,kb
real, dimension(kb) :: cr3,cr1
integer :: k

do k = ka,kb
   cr3(k) = cr1(k)
enddo

return
end subroutine ae1kmic
!###########################################################################

subroutine micro_master

use dump, only: &
  dumpMessage
!use micphys
use rconstants

implicit none
include "constants.f90"
character(len=*),parameter :: h='**(micro_master)**'
integer :: lhcat,khcat,lcat,nd1,nd2,nip,ilcat,ilhcat,idum
integer, dimension(8) :: lcat0

real, dimension(7,16) :: dstprms,dstprms1,dstprms2,dstprms3
real, dimension(16,16) :: jpairr,jpairc
character*80 dataline,cname
logical :: there

data lcat0 /1,2,3,4,5,6,7,16/ ! lcat corressponding to lhcat

data dstprms1/ &
!----------------------------------------------------------------------
! shape      cfmas   pwmas      cfvt    pwvt     dmb0      dmb1
!----------------------------------------------------------------------
    .5,      524.,     3.,    3173.,     2.,   2.e-6,   50.e-6,  & !cloud
    .5,      524.,     3.,     149.,     .5,   .1e-3,    5.e-3,  & !rain
  .179,     110.8,   2.91,  5.769e5,   1.88,  15.e-6,  125.e-6,  & !pris col
  .179,  2.739e-3,   1.74,  188.146,   .933,   .1e-3,   10.e-3,  & !snow col
    .5,      .496,    2.4,    3.084,     .2,   .1e-3,   10.e-3,  & !aggreg
    .5,      157.,     3.,     93.3,     .5,   .1e-3,    5.e-3,  & !graup
    .5,      471.,     3.,     161.,     .5,   .8e-3,   10.e-3,  & !hail
  .429,     .8854,    2.5,     316.,   1.01,      00,       00,  & !pris hex
 .3183,   .377e-2,     2.,     316.,   1.01,      00,       00,  & !pris den
 .1803,   1.23e-3,    1.8,  5.769e5,   1.88,      00,       00,  & !pris ndl
    .5,     .1001,  2.256,   3.19e4,   1.66,      00,       00,  & !pris ros
 .0429,     .8854,    2.5,    4.836,    .25,      00,       00,  & !snow hex
 .3183,   .377e-2,     2.,    4.836,    .25,      00,       00,  & !snow den
 .1803,   1.23e-3,    1.8,  188.146,   .933,      00,       00,  & !snow ndl
    .5,     .1001,  2.256,  1348.38,  1.241,      00,       00,  & !snow ros
    .5,      524.,     3.,    3173.,     2.,  65.e-6,  100.e-6/    !drizzle

data dstprms2/ &
!----------------------------------------------------------------------
! shape      cfmas   pwmas      cfvt    pwvt     dmb0      dmb1
!----------------------------------------------------------------------
    .5,      524.,     3.,    3173.,     2.,   2.e-6,   50.e-6,  & !cloud
    .5,      524.,     3.,     144.,   .497,   .1e-3,    5.e-3,  & !rain
  .179,     110.8,   2.91,    1538.,   1.00,  15.e-6,  125.e-6,  & !pris col
  .179,  2.739e-3,   1.74,     27.7,   .484,   .1e-3,   10.e-3,  & !snow col
    .5,      .496,    2.4,     16.1,   .416,   .1e-3,   10.e-3,  & !aggreg
    .5,      157.,     3.,     332.,   .786,   .1e-3,    5.e-3,  & !graup
    .5,      471.,     3.,    152.1,   .497,   .8e-3,   10.e-3,  & !hail
  .429,     .8854,    2.5,   20801.,  1.377,      00,       00,  & !pris hex
 .3183,   .377e-2,     2.,     56.4,   .695,      00,       00,  & !pris den
 .1803,   1.23e-3,    1.8,   1617.9,   .983,      00,       00,  & !pris ndl
    .5,     .1001,  2.256,    6239.,   1.24,      00,       00,  & !pris ros
  .429,     .8854,    2.5,    30.08,   .563,      00,       00,  & !snow hex
 .3183,   .377e-2,     2.,     3.39,   .302,      00,       00,  & !snow den
 .1803,   1.23e-3,    1.8,     44.6,   .522,      00,       00,  & !snow ndl
    .5,     .1001,  2.256,    125.7,   .716,      00,       00,  & !snow ros
    .5,      524.,     3.,   1.26e7,   1.91,  65.e-6,  100.e-6/    !drizzle

data dstprms3/ &
!----------------------------------------------------------------------
! shape      cfmas   pwmas      cfvt    pwvt     dmb0      dmb1
!----------------------------------------------------------------------
    .5,      524.,     3.,    3173.,     2.,   2.e-6,   50.e-6,  & !cloud
    .5,      524.,     3.,     144.,   .497,   .1e-3,    5.e-3,  & !rain
  .179,     110.8,   2.91,    1538.,   1.00,  15.e-6,  125.e-6,  & !pris col
  .179,  2.739e-3,   1.74,  188.146,   .933,   .1e-3,   10.e-3,  & !snow col
    .5,      .496,    2.4,    3.084,     .2,   .1e-3,   10.e-3,  & !aggreg
    .5,      157.,     3.,     332.,   .786,   .1e-3,    5.e-3,  & !graup
    .5,      471.,     3.,    152.1,   .497,   .8e-3,   10.e-3,  & !hail
  .429,     .8854,    2.5,   20801.,  1.377,      00,       00,  & !pris hex
 .3183,   .377e-2,     2.,     56.4,   .695,      00,       00,  & !pris den
 .1803,   1.23e-3,    1.8,   1617.9,   .983,      00,       00,  & !pris ndl
    .5,     .1001,  2.256,    6239.,   1.24,      00,       00,  & !pris ros
  .429,     .8854,    2.5,    30.08,   .563,      00,       00,  & !snow hex
 .3183,   .377e-2,     2.,     3.39,   .302,      00,       00,  & !snow den
 .1803,   1.23e-3,    1.8,     44.6,   .522,      00,       00,  & !snow ndl
    .5,     .1001,  2.256,    125.7,   .716,      00,       00,  & !snow ros
    .5,      524.,     3.,   1.26e7,   1.91,  65.e-6,  100.e-6/    !drizzle

data jpairr/  &
     0,  0,  0,  1,  2,  3,  4,  0,  0,  0,  0,  5,  6,  7,  8,  0,  &
     0,  0,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,  0,  &
     0, 22, 23, 24,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
    25, 26, 27, 28,  0,  0,  0, 29, 30, 31, 32,  0,  0,  0,  0, 33,  &
    34, 35, 36, 37,  0,  0,  0, 38, 39, 40, 41, 42, 43, 44, 45, 46,  &
    47, 48, 49, 50, 51,  0,  0, 52, 53, 54, 55, 56, 57, 58, 59, 60,  &
    61, 62, 63, 64, 65, 66,  0, 67, 68, 69, 70, 71, 72, 73, 74, 75,  &
     0, 76,  0, 77,  0,  0,  0, 78,  0,  0,  0, 79, 80, 81, 82,  0,  &
     0, 83,  0, 84,  0,  0,  0,  0, 85,  0,  0, 86, 87, 88, 89,  0,  &
     0, 90,  0, 91,  0,  0,  0,  0,  0, 92,  0, 93, 94, 95, 96,  0,  &
     0, 97,  0, 98,  0,  0,  0,  0,  0,  0, 99,100,101,102,103,  0,  &
   104,105,106,  0,  0,  0,  0,107,108,109,110,111,  0,  0,  0,112,  &
   113,114,115,  0,  0,  0,  0,116,117,118,119,  0,120,  0,  0,121,  &
   122,123,124,  0,  0,  0,  0,125,126,127,128,  0,  0,129,  0,130,  &
   131,132,133,  0,  0,  0,  0,134,135,136,137,  0,  0,  0,138,139,  &
     0,  0,  0,140,141,142,143,  0,  0,  0,  0,144,145,146,147,  0/

data jpairc/  &
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
     0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
     0,  2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
     4,  5,  6,  7,  0,  0,  0,  8,  9, 10, 11,  0,  0,  0,  0, 12,  &
    13, 14, 15, 16, 17,  0,  0, 18, 19, 20, 21, 22, 23, 24, 25, 26,  &
    27, 28, 29, 30, 31, 32,  0, 33, 34, 35, 36, 37, 38, 39, 40, 41,  &
    42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,  &
     0, 58,  0,  0,  0,  0,  0, 59,  0,  0,  0,  0,  0,  0,  0,  0,  &
     0, 60,  0,  0,  0,  0,  0,  0, 61,  0,  0,  0,  0,  0,  0,  0,  &
     0, 62,  0,  0,  0,  0,  0,  0,  0, 63,  0,  0,  0,  0,  0,  0,  &
     0, 64,  0,  0,  0,  0,  0,  0,  0,  0, 65,  0,  0,  0,  0,  0,  &
    66, 67, 68,  0,  0,  0,  0, 69, 70, 71, 72, 73,  0,  0,  0, 74,  &
    75, 76, 77,  0,  0,  0,  0, 78, 79, 80, 81,  0, 82,  0,  0, 83,  &
    84, 85, 86,  0,  0,  0,  0, 87, 88, 89, 90,  0,  0, 91,  0, 92,  &
    93, 94, 95,  0,  0,  0,  0, 96, 97, 98, 99,  0,  0,  0,100,101,  &
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/

!  Define several parameters from above data list

do lhcat=1,nhcat
   !Using original RAMS 4.3 power laws
   if(iplaws==0) then
     dstprms(1,lhcat) = dstprms1(1,lhcat)
     dstprms(2,lhcat) = dstprms1(2,lhcat)
     dstprms(3,lhcat) = dstprms1(3,lhcat)
     dstprms(4,lhcat) = dstprms1(4,lhcat)
     dstprms(5,lhcat) = dstprms1(5,lhcat)
     dstprms(6,lhcat) = dstprms1(6,lhcat)
     dstprms(7,lhcat) = dstprms1(7,lhcat)
   !Using Carver/Mitchell 1996 power laws
   elseif(iplaws==1) then
     dstprms(1,lhcat) = dstprms2(1,lhcat)
     dstprms(2,lhcat) = dstprms2(2,lhcat)
     dstprms(3,lhcat) = dstprms2(3,lhcat)
     dstprms(4,lhcat) = dstprms2(4,lhcat)
     dstprms(5,lhcat) = dstprms2(5,lhcat)
     dstprms(6,lhcat) = dstprms2(6,lhcat)
     dstprms(7,lhcat) = dstprms2(7,lhcat)
   !Using mix of original and Carver/Mitchell 1996 power laws
   !Faster falling P,G,H but slower S,A
   elseif(iplaws==2) then
     dstprms(1,lhcat) = dstprms3(1,lhcat)
     dstprms(2,lhcat) = dstprms3(2,lhcat)
     dstprms(3,lhcat) = dstprms3(3,lhcat)
     dstprms(4,lhcat) = dstprms3(4,lhcat)
     dstprms(5,lhcat) = dstprms3(5,lhcat)
     dstprms(6,lhcat) = dstprms3(6,lhcat)
     dstprms(7,lhcat) = dstprms3(7,lhcat)
   endif

   var_shape(lhcat) = dstprms(1,lhcat)
   cfmas(lhcat) = dstprms(2,lhcat)
   pwmas(lhcat) = dstprms(3,lhcat)
   cfvt (lhcat) = dstprms(4,lhcat)
   pwvt (lhcat) = dstprms(5,lhcat)

   do khcat=1,nhcat
      ipairc(lhcat,khcat) = jpairc(lhcat,khcat)
      ipairr(lhcat,khcat) = jpairr(lhcat,khcat)
   enddo
enddo

do lcat=1,ncat
   lhcat = lcat0(lcat)
   emb0 (lcat) = cfmas(lhcat) * dstprms(6,lhcat) ** pwmas(lhcat)
   emb1 (lcat) = cfmas(lhcat) * dstprms(7,lhcat) ** pwmas(lhcat)
enddo

if (level .ne. 3) return

if(mkcoltab.lt.0.or.mkcoltab.gt.1)then
   print*, 'mkcoltab set to ',mkcoltab, 'which is out of bounds'
   stop 'mkcoltab'
endif

cname=coltabfn(1:len_trim(coltabfn))

if(mkcoltab.eq.1)then

! Make collection table and write to file

   call mkcoltb
   open(91,file=cname,form='formatted',status='unknown')
   rewind(91)
   write(91,181)
   do lcat = 1,ncat
      write(91,182)lcat,gnu(lcat),emb0(lcat),emb1(lcat)
   enddo
   write(91,180)
   write(91,183)
   do lhcat = 1,nhcat
      write(91,182)lhcat,cfmas(lhcat),pwmas(lhcat)  &
         ,cfvt(lhcat),pwvt(lhcat)
   enddo
   write(91,180)
   do nip=1,npairc
      write(91,186)nip
      write(91,184)(nd2,(coltabc(nd1,nd2,nip)  &
         ,nd1=1,nembc),nd2=1,nembc)
   enddo
   write(91,180)
   do nip=1,npairr
      write(91,187)nip
      write(91,184)(nd2,(coltabr(nd1,nd2,nip)  &
         ,nd1=1,nembc),nd2=1,nembc)
   enddo

else

!  Read collection table
   inquire(file=cname,exist=there)
   if(.NOT.there) &!call fatal_error("col_table not found")
   iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
     "col_table not found")

   open(91,file=cname,form='formatted',status='old')
   read(91,185)dataline
   do ilcat = 1,ncat
      read(91,182)lcat,gnu(lcat),emb0(lcat),emb1(lcat)
   enddo
   read(91,185)dataline
   read(91,185)dataline
   do ilhcat = 1,nhcat
      read(91,182)lhcat,cfmas(lhcat),pwmas(lhcat)  &
         ,cfvt(lhcat),pwvt(lhcat)
   enddo
   read(91,185)dataline
   do nip=1,npairc
      read(91,185)dataline
      read(91,184)(idum,(coltabc(nd1,nd2,nip)  &
         ,nd1=1,nembc),nd2=1,nembc)
   enddo
   read(91,185)dataline
   do nip=1,npairr
      read(91,185)dataline
      read(91,184)(idum,(coltabr(nd1,nd2,nip)  &
         ,nd1=1,nembc),nd2=1,nembc)
   enddo

endif

180  format(' ')
181  format(' lcat    gnu        emb0       emb1    ')
182  format(i4,7e11.4)
183  format(' lhcat  cfmas      pwmas       cfvt       pwvt')
184  format(i3,20f6.2)
185  format(a80)
186  format('ipairc',i4)
1186 format(6x,i4)
187  format('ipairr',i4)

close(91)

return
end subroutine micro_master

!******************************************************************************


subroutine initqin(n1,n2,n3,q2,q6,q7,pi0,pp,theta,dn0)

!use micphys
use rconstants

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: q2,q6,q7,pi0,pp,theta,dn0

! Initialize Q2, Q6, Q7, CCN, IFN.

do j = 1,n3
   do i = 1,n2
      do k = 1,n1
         pitot(k) = pi0(k,i,j) + pp(k,i,j)
         tair(k) = theta(k,i,j) * pitot(k) / cp

         if(irain .ge. 1) q2(k,i,j) = tair(k) - 193.16
         if(igraup .ge. 1) q6(k,i,j) = 0.5 * min(0.,tair(k) - 273.15)
         if(ihail .ge. 1) q7(k,i,j) = 0.5 * min(0.,tair(k) - 273.15)

      enddo
   enddo
enddo
return
end subroutine initqin
!******************************************************************************

subroutine initqin2(n1,n2,n3,cccnp,cccmp,dn0)

!use micphys
use rconstants
use mem_grid

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: cccnp,cccmp,dn0

! Initialize  CCN

do j = 1,n3
 do i = 1,n2
  do k = 1,n1

   !Set up 3D homogeneous field from RAMSIN (#/cm3)
   if(icloud.eq. 5) cccnp(k,i,j) = cparm

   !Set up Vertical profile of CCN (#/cm3)
   if(icloud.eq.6)then
    if(k<=2) cccnp(k,i,j)=cparm
    if(k>2.and.zt(k)<=4000.) cccnp(k,i,j)=max(100.,cparm * (1.-zt(k)/4000.))
    if(zt(k)>4000.) cccnp(k,i,j) = 100.
   endif

   !Set up 3D Heterogeneously varying field of CCN (#/cm3)
   if(icloud.eq.7)then
     cccnp(k,i,j)=cparm
     print*,'You must set up a 3D field of CCN in mic_init.f90'
     stop
   endif

   !Set up 3D Field of CCN mass (g/cm3)
   if(icloud.ge.5) then
     rg = cnparm  !CCN median radius in centimeters
     rhosol=1.967 !Avg of NH42SO4=1.769 & NaCl=2.165
     ant=cccnp(k,i,j)
     if(rg<=0.02e-4) then
       rmsma = 1.0e-20
       rmlar = 1.0e-13
     elseif(rg>0.02e-4 .and. rg<=0.04e-4) then
       rmsma = 1.0e-19
       rmlar = 1.0e-12
     elseif(rg>0.04e-4 .and. rg<=0.08e-4) then
       rmsma = 1.0e-18
       rmlar = 1.0e-11
     elseif(rg>0.08e-4 .and. rg<=0.16e-4) then
       rmsma = 1.0e-17
       rmlar = 1.0e-10
     elseif(rg>0.16e-4 .and. rg<=0.36e-4) then
       rmsma = 5.0e-17
       rmlar = 5.0e-10
     elseif(rg>0.36e-4 .and. rg<=0.64e-4) then
       rmsma = 5.0e-16
       rmlar = 5.0e-09
     elseif(rg>0.64e-4) then
       rmsma = 1.0e-15
       rmlar = 1.0e-08
     endif
     power = alog10(rmsma/rmlar) / float(itbin-1)
     do ic=1,itbin
      smass(ic) = rmlar * 10.**(float(ic-1) * power) ! solute masses (g)
      rs(ic)=(3.*smass(ic)/(12.56637*rhosol) )**(.3333333) !radii (cm)
     enddo
     ccnnum = 0.0
     ccnmass = 0.0
     do ic=1,itbin-1
      rcm=0.5*(rs(ic)+rs(ic+1))
      con0(ic)=ant/(1.473362674*rcm)*exp(-(alog(rcm/rg))**2/0.690986327)
      con0(ic)=con0(ic)*(rs(ic)-rs(ic+1))
      if(ccnnum <= cccnp(k,i,j))ccnmass=ccnmass+smass(ic)*con0(ic)
      ccnnum = ccnnum + con0(ic)
     enddo
     cccmp(k,i,j) = ccnmass
   endif

  enddo
 enddo
enddo
return
end subroutine initqin2

!******************************************************************************

subroutine initqin3(n1,n2,n3,gccnp,gccmp,dn0)

!use micphys
use rconstants
use mem_grid

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: gccnp,gccmp,dn0

! Initialize  Giant-CCN

do j = 1,n3
 do i = 1,n2
  do k = 1,n1

   !Set up 3D homogeneous field from RAMSIN (#/cm3)
   if(idriz .eq. 5) gccnp(k,i,j) = dparm

   !Set up Vertical profile of GCCN (#/cm3)
   if(idriz .eq.6)then
    if(k<=4) gccnp(k,i,j)=dparm
    if(k>4.and.zt(k)<=4000.) gccnp(k,i,j)=max(1.e-5,dparm * (1.-zt(k)/4000.))
    if(zt(k)>4000.) gccnp(k,i,j) = 1.e-5
   endif

   !Set up 3D Heterogeneously varying field of GCCN (#/cm3)
   if(idriz .eq.7)then
     gccnp(k,i,j)=dparm
     print*,'You must set up a 3D field of Giant-CCN in mic_init.f90'
     stop
   endif

   !Set up 3D Field of CCN mass (g/cm3)
   if(idriz.ge.5) then
     rg = gnparm  !GCCN median radius in centimeters
     rhosol=1.967 !Avg of NH42SO4=1.769 & NaCl=2.165
     ant=gccnp(k,i,j)
     rmsma = 1.0e-14
     rmlar = 1.0e-05
     power = alog10(rmsma/rmlar) / float(itbin-1)
     do ic=1,itbin
      smass(ic) = rmlar * 10.**(float(ic-1) * power) ! solute masses (g)
      rs(ic)=(3.*smass(ic)/(12.56637*rhosol) )**(.3333333) !radii (cm)
     enddo
     ccnnum = 0.0
     ccnmass = 0.0
     do ic=1,itbin-1
      rcm=0.5*(rs(ic)+rs(ic+1))
      con0(ic)=ant/(1.473362674*rcm)*exp(-(alog(rcm/rg))**2/0.690986327)
      con0(ic)=con0(ic)*(rs(ic)-rs(ic+1))
      if(ccnnum <= gccnp(k,i,j))ccnmass=ccnmass+smass(ic)*con0(ic)
      ccnnum = ccnnum + con0(ic)
     enddo
     gccmp(k,i,j) = ccnmass
   endif

  enddo
 enddo
enddo
return
end subroutine initqin3

!******************************************************************************

subroutine initqin4(n1,n2,n3,cifnp,dn0)

!use micphys
use rconstants

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: cifnp,dn0

! Initialize  IFN

do j = 1,n3
   do i = 1,n2
      do k = 1,n1
        if(ipris .ge. 5) cifnp(k,i,j) = pparm * dn0(k,i,j) ** 5.4
      enddo
   enddo
enddo
return
end subroutine initqin4

!******************************************************************************
subroutine initqin5(n1,n2,n3,md1np,md2np)

!use micphys
use rconstants
use mem_grid

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: md1np,md2np

print*,'idust,imd2flg,imd1flg',idust,imd2flg,imd1flg

do j = 1,n3
   do i = 1,n2
      do k = 1,n1
        !Set up concentration of Mineral Dust (#/cm3)
        if(idust .eq. 0) then
            !If not using dust source model, initialize background dust
            if(imd1flg .eq. 1) then
             if(k<=2) md1np(k,i,j)=1000.
             if(k>2) md1np(k,i,j)=max(20.,1000.*(1.-zm(k)/4000.))
            endif
            if(imd2flg .eq. 1) then
             if(k<=2) md2np(k,i,j)=100.
             if(k>2) md2np(k,i,j)=max(2.e-2,100.*(1.-zm(k)/4000.))
            endif
        elseif(idust .eq. 1) then
            !Set dust median radii (cm) based on dust model
            !These are used in radiation and dust deposition & scavenging
            d1parm = 0.699e-4 !(weighted in situ mean is 0.699e-4)
            d2parm = 2.95e-4  !(weighted in situ mean is 2.95e-4)
            !If using dust source model, do not initialize with background dust
            md1np(k,i,j) = 0.
            md2np(k,i,j) = 0.
        endif
      enddo
   enddo
enddo
return
end subroutine initqin5

!******************************************************************************

subroutine jnmbinit()

!use micphys

implicit none

if (level /= 3) then

   if (level <= 1) then
      jnmb(1) = 0
   else
      jnmb(1) = 4
   endif

   jnmb(2) = 0
   jnmb(3) = 0
   jnmb(4) = 0
   jnmb(5) = 0
   jnmb(6) = 0
   jnmb(7) = 0
   jnmb(8) = 0

else

   jnmb(1) = icloud
   jnmb(2) = irain
   jnmb(3) = ipris
   jnmb(4) = isnow
   jnmb(5) = iaggr
   jnmb(6) = igraup
   jnmb(7) = ihail
   jnmb(8) = idriz

   if (icloud .eq. 1) jnmb(1) = 4
   if (irain  .eq. 1) jnmb(2) = 2
   if (ipris  .ge. 1) jnmb(3) = 5
   if (isnow  .eq. 1) jnmb(4) = 2
   if (iaggr  .eq. 1) jnmb(5) = 2
   if (igraup .eq. 1) jnmb(6) = 2
   if (ihail  .eq. 1) jnmb(7) = 2
   if (idriz  .eq. 1) jnmb(8) = 4

   if (irain == 5 .or. isnow == 5 .or. iaggr == 5 .or.  &
      igraup == 5 .or. ihail == 5) then

      if (irain  >= 1) jnmb(2) = 5
      if (isnow  >= 1) jnmb(4) = 5
      if (iaggr  >= 1) jnmb(5) = 5
      if (igraup >= 1) jnmb(6) = 5
      if (ihail  >= 1) jnmb(7) = 5

   endif

endif
PRINT*," 2M microphysics"
PRINT*,'JNMB(1)=',jnmb(1),' JNMB(2)=',jnmb(2),' JNMB(3)=',jnmb(3)
PRINT*,'JNMB(4)=',jnmb(4),' JNMB(5)=',jnmb(5),' JNMB(6)=',jnmb(6)
PRINT*,'JNMB(7)=',jnmb(7),' JNMB(8)=',jnmb(8)

end subroutine jnmbinit

!******************************************************************************

subroutine micinit()

!use micphys

implicit none

integer :: lhcat,lcat,ia
integer, dimension(16) :: lcat0
real :: cfmasi,c1,glg,glg1,glg2,glgm,glgc,glgmv,gym,flngi,dpsi,embsip,dnsip
real :: gammln,gammp,gammq

data lcat0 /1,2,3,4,5,6,7,3,3,3,3,4,4,4,4,8/ ! lcat corressponding to lhcat

! Initialize arrays based on microphysics namelist parameters

parm(1) = cparm
parm(2) = rparm
!     parm(3) = pparm   [obsolete]
parm(4) = sparm
parm(5) = aparm
parm(6) = gparm
parm(7) = hparm
parm(8) = dparm

if (icloud .le. 1) parm(1) = .3e9
if (irain  .eq. 1) parm(2) = .1e-2
!     if (ipris  .eq. 1) parm(3) = .1e4     [obsolete]
if (isnow  .eq. 1) parm(4) = .1e-2
if (iaggr  .eq. 1) parm(5) = .1e-2
if (igraup .eq. 1) parm(6) = .1e-2
if (ihail  .eq. 1) parm(7) = .3e-2
if (idriz  .eq. 1) parm(8) = .1e6  ! # per kg ~ m^3
                                   !(mid-range avg from Feingold(99)

dps = 125.e-6
dps2 = dps ** 2
rictmin = 1.0001
rictmax = 0.9999 * float(nembc)

do lhcat = 1,nhcat
   lcat=lcat0(lhcat)

   cfden(lhcat) = cfmas(lhcat) * 6.0 / 3.14159
   pwden(lhcat) = pwmas(lhcat) - 3.
   emb0log(lcat) = log(emb0(lcat))
   emb1log(lcat) = log(emb1(lcat))

! Define coefficients [vtfac, frefac1, frefac2] used for terminal velocity
! and Reynolds number

   cfmasi = 1. / cfmas(lhcat)
   pwmasi(lhcat) = 1. / pwmas(lhcat)
   pwen0(lhcat) = 1. / (pwmas(lhcat) + 1.)
   pwemb0(lhcat) = pwmas(lhcat) / (pwmas(lhcat) + 1.)
   c1 = 1.5 + .5 * pwvt(lhcat)

   glg = gammln(gnu(lcat))
   glg1 = gammln(gnu(lcat) + 1.)
   glg2 = gammln(gnu(lcat) + 2.)
   glgm = gammln(gnu(lcat) + pwmas(lhcat))
   glgc = gammln(gnu(lcat) + c1)
   glgmv = gammln(gnu(lcat) + pwmas(lhcat) + pwvt(lhcat))

   if (jnmb(lcat) .eq. 3) then
      cfemb0(lhcat) = cfmas(lhcat) * exp(glgm - glg)  &
         ** pwen0(lhcat) * (1. / parm(lcat)) ** pwemb0(lhcat)
      cfen0(lhcat) = parm(lcat) * (exp(glg - glgm) / parm(lcat))  &
         ** pwen0(lhcat)
   endif

   dnfac(lhcat) = (cfmasi * exp(glg - glgm)) ** pwmasi(lhcat)

   vtfac(lhcat) = cfvt(lhcat) * exp(glgmv - glgm)  &
      * (cfmasi * exp(glg - glgm)) ** (pwvt(lhcat) *pwmasi(lhcat))

   frefac1(lhcat) = var_shape(lhcat) * exp(glg1 - glg)  &
      * (cfmasi * exp(glg - glgm)) ** pwmasi(lhcat)

   frefac2(lhcat) = var_shape(lhcat) * 0.229 * sqrt(cfvt(lcat))  &
      * (cfmasi * exp(glg - glgm)) ** (pwmasi(lhcat) * c1)  &
      * exp(glgc - glg)

   sipfac(lhcat) = .785 * exp(glg2 - glg)  &
      * (cfmasi * exp(glg - glgm)) ** (2. * pwmasi(lhcat))

   cfmasft(lhcat) = cfmas(lhcat) * exp(gammln  &
      (gnu(lcat) + pwmas(lhcat)) - gammln(gnu(lcat)))

   dict(lcat) = float(nembc-1) / (emb1log(lcat) - emb0log(lcat))

   dpsmi(lhcat) = 1. / (cfmas(lhcat) * dps ** pwmas(lhcat))
   if (lhcat .le. 4) gamm(lhcat) = exp(glg)
   if (lhcat .le. 4) gamn1(lhcat) = exp(glg1)

! gam1   :  the integral of the pristine distribution from dps to infty
! gam2   :  the integral of the snow dist. from 0 to dps
! gam3   :  values of the exponential exp(-dps/dn)

enddo

!***********************************************************************
!************** Secondary Ice Production Arrays ************************
!***********************************************************************
flngi = 1. / float(ngam)
do ia=1,ngam
   dpsi = dps * 1.e6 / float(ia)

   gam(ia,1) = gammq(gnu(3) + 1., dpsi)
   gam(ia,2) = gammp(gnu(4) + 1., dpsi)
   gam(ia,3) = exp(-dpsi)

   GAMINC(IA,1)=GAMMQ(GNU(3),dpsi)
   GAMINC(IA,2)=GAMMP(GNU(4),dpsi)

   embsip = emb1(1) * float(ia) * flngi
   dnsip = dnfac(1) * embsip ** pwmasi(1)
   gamsip13(1,ia) = gammp(gnu(1),13.e-6/dnsip)
   gamsip24(1,ia) = gammq(gnu(1),24.e-6/dnsip)

   embsip = emb1(8) * float(ia) * flngi
   dnsip = dnfac(16) * embsip ** pwmasi(16)
   gamsip13(2,ia)= gammp(gnu(8),13.e-6/dnsip)
   gamsip24(2,ia)= gammq(gnu(8),24.e-6/dnsip)
enddo

!***********************************************************************
!************* CCN tracking and spectrum tables ************************
!***********************************************************************
!The Feingold parcel model calculates the log-normal CCN distribution
!with respect to the temperature and pressure of the parcel. This is
!adequate for the parcel model in deciding how many CCN get to nucleate
!depending upon the ambient conditions of T,W,P,#concen. For a simple
!standard distribution used to subtract CCN after nucleation, these
!dependencies should not be used. We simply want to know how much mass
!corressponds to a given number of CCN. The cldnuc lookup table already
!contains the environment influences in choosing what percent of available
!CCN to nucleate. The value ccncon is the minimum distribution for CCN
!which has a number concentration of 1/CC. The cldnuc routine calculates
!the necessary multiples of the chosen table value for nucleation.
!*****************************************************************************
!The largest value of the 1st print statement displays the bin in which the
!% change in number between adjacent bins is greatest over the % change in
!mass between adjacent bins. This is the point of the curve where % increase
!in number most closely equals the % increase in mass as one follows the curve.
!By eliminating CCN from the distribution from the largest mass bin up to this
!point, one can eliminate the tail overlaps between curves of different median
!radius while minimizing the total loss of mass and number. Since the cloud
!nucleation routine nucleates a certain number, the removal of the tail only
!removes the very few large CCN that are anomalous in nature anyway. Such
!giant CCN can be treated separately in RAMS. By doing this we prevent over
!estimating the mass of CCN that is nucleated. This is not done for GCCN!
!Perhaps since the number of GCCN is relatively much smaller, the depletion
!from the curves behaves a bit differently. In some cases either not enough
!mass is depleted, or the tendency equations influence the GCCN mass and
!number such that the median radius increases from one timestep to the next
!when only GCCN depletion is allowed. Typically this should shift the spectrum
!to smaller median radii.
!*****************************************************************************

if(iccnlev>=1) then
  do rgb = 1,9

   !Set up binned distribution mass and sizes
    rg = rg_ccn(rgb)
    if(rg<=0.02e-4) then
      rmsma = 1.0e-20
      rmlar = 1.0e-13
    elseif(rg>0.02e-4 .and. rg<=0.04e-4) then
      rmsma = 1.0e-19
      rmlar = 1.0e-12
    elseif(rg>0.04e-4 .and. rg<=0.08e-4) then
      rmsma = 1.0e-18
      rmlar = 1.0e-11
    elseif(rg>0.08e-4 .and. rg<=0.16e-4) then
      rmsma = 1.0e-17
      rmlar = 1.0e-10
    elseif(rg>0.16e-4 .and. rg<=0.36e-4) then
      rmsma = 5.0e-17
      rmlar = 5.0e-10
    elseif(rg>0.36e-4 .and. rg<=0.64e-4) then
      rmsma = 5.0e-16
      rmlar = 5.0e-09
    elseif(rg>0.64e-4) then
      rmsma = 1.0e-15
      rmlar = 1.0e-08
    endif

    power = alog10(rmsma/rmlar) / float(itbin-1)
    do ic=1,itbin
     !bin radius = [(3/4) x (1/pi) x (1/rho) x mass of bin] ^ (1/3)
     !rho(ammonium sulfate=1.769), rho(salt=2.165), Average=1.967g/cm3
     !The two rho values are averaged since we end up with aerosol
     !mixing and the 2 chemistries become indistinguishable
     smass(ic) = rmlar * 10.**(float(ic-1) * power) !solute masses (g)
     binrad(ic)=(0.1213688*smass(ic))**(.33333) !radius (cm)
    enddo

    !Set up binned distribution concentration
    rg = rg_ccn(rgb)
    do ic=1,itbin-1
     rcm=0.5*(binrad(ic)+binrad(ic+1))
     ccncon(ic,rgb) = 1.0/(1.47336267*rcm)*exp(-(alog(rcm/rg))**2/0.6909863)
     ccncon(ic,rgb) = ccncon(ic,rgb)*(binrad(ic)-binrad(ic+1))
     ccnmas(ic,rgb) = smass(ic) * ccncon(ic,rgb)
     if(ic > 1) then
        ccncon(ic,rgb) = ccncon(ic,rgb) + ccncon(ic-1,rgb)
        ccnmas(ic,rgb) = ccnmas(ic,rgb) + ccnmas(ic-1,rgb)
     endif
    enddo

    !Eliminate the largest CCN. This represents a very small #/mass ratio.
    pctccnmas = .085 * ccnmas(itbin-1,rgb)
    do ic=1,itbin-1
      if(ccnmas(ic,rgb) < pctccnmas) then
        ccnmass = ccnmas(ic,rgb)
        ccnnum  = ccncon(ic,rgb)
        ccnmas(ic,rgb) = 0.0
        ccncon(ic,rgb) = 0.0
        if(ccnmas(ic+1,rgb) >= pctccnmas) go to 101
      endif
101  continue
     enddo
     do ic=1,itbin-1
      if(ccnmas(ic,rgb) > 0.0) then
       ccnmas(ic,rgb) = ccnmas(ic,rgb) - ccnmass
       ccncon(ic,rgb) = ccncon(ic,rgb) - ccnnum
      endif
    enddo

  enddo !rgb loop
endif !ccnlev ifstatement

return
end subroutine micinit


!###########################################################################

subroutine haznuc()

!use micphys

implicit none

integer :: ithz,irhhz,k
real :: denccn,gnuccn,dnccn,ddccn,rhhz,c1hz,c2hz,c3hz,bhz,dm,sum  &
       ,dccn,y,dum,thz
real :: gammln

!  Haze nucleation table

denccn = 1.769    !Density of ammonium sulfate
gnuccn = 2.       !Saleeby(02-21-2007) Originally = 1.
dnccn = 2.*cnparm !Saleeby(02-21-2007) Originally = .075e-4
ddccn = .005e-4
do ithz = 1,nthz
   thz = -60. + dthz * float(ithz - 1)
   do irhhz = 1,nrhhz
      rhhz = 0.82 + drhhz * float(irhhz - 1)
      c1hz = (3.14159 * denccn / 6.) ** (-.333333)
      c2hz = -14.65 - 1.045 * thz
      c3hz = -492.35 - 8.34 * thz - 0.0608 * thz ** 2
      bhz = min(38., max(-38., c2hz + c3hz * (1. - rhhz)))
      dm = c1hz * 10 ** (-bhz/6.)

      sum = 0.
      dccn = 0.
      do k=1,200
         dccn = dccn + ddccn
         y=dccn / dnccn
         dum=min(50., (dccn / dm) ** 6)
         sum = sum + y ** (gnuccn-1.) * exp(-y) * (1. - exp(-dum))
      enddo
      frachz(irhhz,ithz) = sum*ddccn/(exp(gammln(gnuccn))*dnccn)
   enddo
enddo

return
end subroutine haznuc

!******************************************************************************

subroutine homfrzcl(dtlt,ngr)

!use micphys

implicit none

integer :: itc,ngr,k,idnc
real :: gnuc,ddc,ajlso,dnc,sum,dc,v1,tc,y,dtlt
real :: gammln

!  Make table for homogeneous freezing of cloud droplets
!  Uses characteristic diameter instead of true diameter

ddc = 0.5e-6
!ndnc = 11     | both preset in micphys.h
!ddnc = 2.e-6  |

do itc = 1,ntc
   tc = -50. + dtc * float(itc-1)
   y = -(606.3952+tc*(52.6611+tc*(1.7439+tc*(.0265+tc*1.536e-4))))
   ajlso = 1.e6 * 10. ** y

   do idnc = 1,ndnc
      dnc = ddnc * float(idnc)
      sum = 0.
      dc = 0.
      do k = 1,2000
         dc = dc + ddc
         v1 = 0.523599 * dc ** 3
         sum = sum + (dc / dnc) ** (gnu(1) - 1.) * exp(-dc / dnc)  &
            * (1. - exp(-ajlso * v1 * dtlt))
      enddo
      fracc(idnc,itc,ngr) = sum * ddc / (exp(gammln(gnu(1))) * dnc)
   enddo

enddo
return
end subroutine homfrzcl

!******************************************************************************

! As before, sedimentation is not yet designed to transfer hydrometeor mass
! between grids in the case where a nested grid does not reach the top and/or
! bottom of the model domain.  Thus, vertical nested grid boundaries should
! be avoided where sedimentation occurs.

subroutine mksedim_tab(m1,m2,m3,ngr,nembfall,maxkfall  &
   ,zm,dzt,pcpfillc,pcpfillr,sfcpcp,allpcp,dtsed)

!use micphys

implicit none

integer, parameter :: nbin=50
integer :: m1,m2,m3,iembs,lcat,lhcat,k,kkf,ibin,kk,nembfall,maxkfall  &
          ,jbin,ngr
integer, dimension(16) :: lcat0

real, dimension(16) :: dmbsed

real :: dmbodn,diam0,diam1,fac1,fac3,sumc,sumr,diam,fac2,fac4  &
       ,disp,ztopnew,zbotnew,fallin,delzsfc,dispemb,dispmax,dispmx
real :: gammln,gammp,dtsed,dispmax0,dispmax1
real, dimension(m1) :: zm,dzt
real, dimension(m1,maxkfall,nembfall,nhcat) :: pcpfillc,pcpfillr
real, dimension(maxkfall,nembfall,nhcat) :: sfcpcp
real, dimension(m1,nembfall,nhcat) :: allpcp
real, dimension(nbin) :: cbin,rbin,reldisp

data lcat0 /1,2,3,4,5,6,7,3,3,3,3,4,4,4,4,8/

data dmbsed/40.e-6,7.e-3,750.e-6,10.e-3,15.e-3,15.e-3,30.e-3,750.e-6 &
   ,750.e-6,750.e-6,750.e-6,10.e-3,10.e-3,10.e-3,10.e-3,150.e-6/

! Because timestep may now be variable in time, define sedtime0 and sedtime1
! here as 0.1 seconds and 3000 seconds.  The former is supposed to be
! less than 0.7 of the shortest timestep on any grid (sqrt(dn0i) never exceeds
! 0.7) and the latter is the longest timestep expected to ever be used (300
! seconds) times a factor of 2 for the maximum inverse of rtgt times a factor
! of 5 for the largest value of sqrt(dn0i).

sedtime0 = .1
sedtime1 = 3000.
dispmax0 = 500.

! Loop over hydrometeor categories

do lhcat = 1,nhcat
   lcat = lcat0(lhcat)

   !Displacement max based on largest size and timestep (Saleeby 8-17-05)
   dispmax1 = dtsed * cfvt(lhcat) * dmbsed(lhcat) ** pwvt(lhcat)
   dispmax = min(dispmax0,dispmax1)

   dispemb0(lhcat,ngr) = sedtime0 * cfvt(lhcat)  &
      * (emb0(lcat) / cfmas(lhcat)) ** (pwvt(lhcat) * pwmasi(lhcat))

   dispemb1(lhcat,ngr) = sedtime1 * cfvt(lhcat)  &
      * (emb1(lcat) / cfmas(lhcat)) ** (pwvt(lhcat) * pwmasi(lhcat))

!Bob (10/24/00):  Limit dispemb1 to a maximum of dispmax

   if (dispemb1(lhcat,ngr) .gt. dispmax) dispemb1(lhcat,ngr) = dispmax

! Loop over bins, filling them with fractional number, fractional mass,
! and displacement quotient relative to emb.

   dmbodn = (exp(gammln(gnu(lcat) + pwmas(lhcat))  &
      - gammln(gnu(lcat)))) ** pwmasi(lhcat)
   diam0 = 0.06 * dmbodn
   diam1 = 1.0 * dmbodn
   fac1 = gammp(gnu(lcat),diam0)
   fac3 = gammp(gnu(lcat) + pwmas(lhcat),diam0)
   sumc = 0.
   sumr = 0.

   do jbin = 1,nbin

      diam = diam0 * (diam1 / diam0) ** (float(jbin)/float(nbin))
      fac2 = gammp(gnu(lcat),diam)
      fac4 = gammp(gnu(lcat) + pwmas(lhcat),diam)
      cbin(jbin) = fac2 - fac1
      rbin(jbin) = fac4 - fac3
      fac1 = fac2
      fac3 = fac4
      sumc = sumc + cbin(jbin)
      sumr = sumr + rbin(jbin)
      reldisp(jbin) = diam ** pwvt(lhcat)

   enddo

   do jbin = 1,nbin
      cbin(jbin) = cbin(jbin) / sumc
      rbin(jbin) = rbin(jbin) / sumr
   enddo

! Loop over displacement distance for size emb.

   do iembs = 1,nembfall
      dispemb = dispemb0(lhcat,ngr)  &
         * (dispemb1(lhcat,ngr) / dispemb0(lhcat,ngr))  &
         ** (float(iembs-1) / float(nembfall-1))

! Zero out concentration and mass fill arrays and surface precip array
! before accumulation.

      do k = 1,m1
         do kkf = 1,maxkfall
            pcpfillc(k,kkf,iembs,lhcat) = 0.
            pcpfillr(k,kkf,iembs,lhcat) = 0.
         enddo
         allpcp(k,iembs,lhcat) = 0.
         if (k .le. maxkfall) sfcpcp(k,iembs,lhcat) = 0.
      enddo

! Loop over vertical grid index.

      do k = 2,m1-1

!Bob (10/24/00):  Limit disp distance to (maxkfall-1) levels

         dispmx = dispmax
         if (k .gt. maxkfall) then
            dispmx = min(dispmx,zm(k-1) - zm(k-maxkfall))
         endif

! Loop over bins

         do ibin = 1,nbin
            disp = dispemb * reldisp(ibin)
            if (disp .gt. dispmx) disp = dispmx
            ztopnew = zm(k) - disp
            zbotnew = zm(k-1) - disp

! Loop over grid cells that a parcel falls into.

            do kkf = 1,min(k-1,maxkfall)

               kk = k + 1 - kkf
               if (zbotnew .gt. zm(kk)) go to 50

               if (ztopnew .le. zm(kk-1)) then
                  fallin = 0.
               else
                  fallin = dzt(kk) *  &
                     (min(zm(kk),ztopnew) - max(zm(kk-1),zbotnew))
               endif

               pcpfillc(k,kkf,iembs,lhcat) = pcpfillc(k,kkf,iembs,lhcat)  &
                  + fallin * cbin(ibin)

               pcpfillr(k,kkf,iembs,lhcat) = pcpfillr(k,kkf,iembs,lhcat)  &
                  + fallin * rbin(ibin)

            enddo
50          continue

! Compute precip rate at levels above ground
            allpcp(k,iembs,lhcat) = allpcp(k,iembs,lhcat) + disp * rbin(ibin)

! Compute surface precipitation.
            if (zbotnew .lt. 0.) then
               delzsfc = min(0.,ztopnew) - zbotnew
               if (k .le. maxkfall) sfcpcp(k,iembs,lhcat)  &
                  = sfcpcp(k,iembs,lhcat) + delzsfc * rbin(ibin)
            endif

         enddo

      enddo
   enddo
enddo

return
end subroutine mksedim_tab

!******************************************************************************

subroutine tabmelt()

!!use micphys

implicit none

integer, parameter :: nbins=500

integer :: lhcat,lcat,ndns1,ibin,inc,iter,idns
integer, dimension(16) :: lcat0
real :: dn,gammaa,totfmg,totmass,vtx,fre,totqm,qmgoal,qmnow,totmdqdt,deltat  &
       ,pliqmass,picemass,critmass,vk,gammaa8,gammaa2,dn8,dn2

real, dimension(nbins) :: db,fmg,pmass,binmass,dqdt,q,db2,db8,fmg8,fmg2
real, dimension(ncat) :: dmean
real :: gammln
data lcat0/1,2,3,4,5,6,7,3,3,3,3,4,4,4,4,8/
data dmean/20.e-6,500.e-6,30.e-6,500.e-6,500.e-6,500.e-6,8000.e-6,60.e-6/
data vk/0.2123e-04/

do lhcat = 1,nhcat
   lcat = lcat0(lhcat)
   dn = dmean(lcat) / gnu(lcat)
   gammaa = exp(gammln(gnu(lcat)))

   rmlttab(1,lhcat) = 0.0
   rmlttab(ninc,lhcat) = 1.0
   enmlttab(1,lhcat) = 0.0
   enmlttab(ninc,lhcat) = 1.0

   ndns1 = 1
   if (lcat .eq. 7) ndns1 = ndns

   do idns = 1,ndns1
      shedtab(1,idns) = 0.0
      shedtab(ninc,idns) = 0.0

      if (ndns1 .gt. 1) dn = 1.e-3 * float(idns) / gnu(lcat)

      totfmg = 0.
      totmass = 0.
      do ibin = 1,nbins
         db(ibin) = 0.02 * dn * (float(ibin) - 0.5)
         fmg(ibin) = (db(ibin) / dn) ** (gnu(lcat) - 1.)  &
            / (dn * gammaa) * exp(-db(ibin) / dn)
         totfmg = totfmg + fmg(ibin)
         q(ibin) = 0.
         pmass(ibin) = cfmas(lhcat) * db(ibin) ** pwmas(lhcat)
         binmass(ibin) = pmass(ibin) * fmg(ibin)
         totmass = totmass + binmass(ibin)
         vtx = cfvt(lhcat) * db(ibin) ** pwvt(lhcat)
         fre = (1.0 + 0.229 * sqrt(vtx * db(ibin) / vk))  &
            * var_shape(lhcat)
         dqdt(ibin) = db(ibin) ** (1. - pwmas(lhcat)) * fre
      enddo
      totqm = totmass * 80.

      do inc = 2,ninc-1
         qmgoal = totqm * float(inc-1) / float(ninc-1)

         do iter = 1,2
            qmnow = 0.
            totmdqdt = 0.
            do ibin = 1,nbins
               if(q(ibin) .lt. 79.9999)then
                  totmdqdt = totmdqdt + binmass(ibin) * dqdt(ibin)
               endif
               qmnow = qmnow + q(ibin) * binmass(ibin)
            enddo
            deltat = max(0.,(qmgoal - qmnow) / totmdqdt)
            do ibin = 1,nbins
               q(ibin) = min(80.,q(ibin) + dqdt(ibin) * deltat)
            enddo
         enddo

!  For the current inc value (representing total liquid fraction), compute
!  melted mixing ratio (rmlttab) and number (enmlttab) from totally-melted
!  bins and compute shedded mixing ratio (shedtab) from partially-melted bins.

         if(idns .eq. 7 .or. ndns1 .eq. 1)then
            rmlttab(inc,lhcat) = 0.0
            do ibin = 1,nbins
               if(q(ibin) .gt. 79.9)then
                  rmlttab(inc,lhcat) = rmlttab(inc,lhcat) + binmass(ibin)
               endif
            enddo
            rmlttab(inc,lhcat) = rmlttab(inc,lhcat) / totmass
         endif

         if(idns .eq. 7 .or. ndns1 .eq. 1)then
            enmlttab(inc,lhcat) = 0.0
            do ibin = 1,nbins
               if(q(ibin) .gt. 79.9)then
                  enmlttab(inc,lhcat) = enmlttab(inc,lhcat)  &
                     + fmg(ibin)
               endif
            enddo
            enmlttab(inc,lhcat) = enmlttab(inc,lhcat) / totfmg
         endif

         if(lcat .eq. 7)then
            shedtab(inc,idns) = 0.0
!                  do ibin = kbin,nbins
            do ibin = 1,nbins
               if(q(ibin) .le. 79.9)then
                  pliqmass = pmass(ibin) * q(ibin) / 80.
                  picemass = pmass(ibin) - pliqmass
                  critmass = .268e-3 + .1389 * picemass
                  shedtab(inc,idns) = shedtab(inc,idns)  &
                     + max(0.0, pliqmass - critmass) * fmg(ibin)
               endif
            enddo
            shedtab(inc,idns) = shedtab(inc,idns) / totmass
         endif

      enddo
   enddo
enddo
return
end subroutine tabmelt

!******************************************************************************

subroutine mkcoltb

!!use micphys

implicit none

integer, parameter :: ndx=20
integer :: ihx,ix,ihy,iy,iemby,iembx,idx
integer, dimension(16) :: ix0,iy0
real :: gxm,dnminx,dnmaxx,dxlo,dxhi,gyn,gyn1,gyn2,gynp,gynp1,gynp2,gym  &
       ,dnminy,dnmaxy,dny,vny,dnx,ans
real :: gammln,xj

real, dimension(ndx) :: dx,fx,gx
data ix0/1,2,3,4,5,6,7,3,3,3,3,4,4,4,4,8/
data iy0/1,2,3,4,5,6,7,3,3,3,3,4,4,4,4,8/

do ihx = 1,nhcat

   ix = ix0(ihx)

   gxm = exp(gammln(gnu(ix)) - gammln(gnu(ix) + pwmas(ihx)))
   dnminx = ((emb0(ix) / cfmas(ihx)) * gxm) ** (1. / pwmas(ihx))
   dnmaxx = ((emb1(ix) / cfmas(ihx)) * gxm) ** (1. / pwmas(ihx))
   dxlo = .01 * dnminx
   dxhi = 10. * dnmaxx

do ihy = 1,nhcat

   print*, 'ihx,ihy,ipairc,ipairr',ihx,ihy,ipairc(ihx,ihy)  &
      ,ipairr(ihx,ihy)

   iy = iy0(ihy)

   if (ipairc(ihx,ihy) .gt. 0 .or. ipairr(ihx,ihy) .gt. 0) then
      gyn = exp(gammln(gnu(iy)))
      gyn1 = exp(gammln(gnu(iy) + 1.)) / gyn
      gyn2 = exp(gammln(gnu(iy) + 2.)) / gyn
      gynp = exp(gammln(gnu(iy) + pwvt(ihy))) / gyn
      gynp1 = exp(gammln(gnu(iy) + pwvt(ihy) + 1.)) / gyn
      gynp2 = exp(gammln(gnu(iy) + pwvt(ihy) + 2.)) / gyn

      gym = exp(gammln(gnu(iy)) - gammln(gnu(iy) + pwmas(ihy)))
      dnminy = ((emb0(iy) / cfmas(ihy)) * gym) ** (1. /pwmas(ihy))
      dnmaxy = ((emb1(iy) / cfmas(ihy)) * gym) ** (1. /pwmas(ihy))

      do iemby = 1,nembc
         dny = dnminy * (dnmaxy / dnminy) ** (float(iemby-1)  &
            / float(nembc-1))
         vny = cfvt(ihy) * dny ** pwvt(ihy)
         do iembx = 1,nembc

            dnx = dnminx * (dnmaxx / dnminx) ** (float(iembx-1)  &
               / float(nembc-1))
            do idx = 1,ndx
               dx(idx) = dxlo * (dxhi / dxlo)  &
                  ** (float(idx-1) / float(ndx-1))
               fx(idx) = xj(dx(idx),cfvt(ihx),pwvt(ihx),cfvt(ihy)  &
                  ,pwvt(ihy),vny,dnx,dny,gnu(ix),gnu(iy)  &
                  ,gyn1,gyn2,gynp,gynp1,gynp2)
               gx(idx) = fx(idx) * cfmas(ihx)  &
                  * dx(idx) ** pwmas(ihx)

            enddo
            if (ipairc(ihx,ihy) .gt. 0) then
               call avint(dx,fx,ndx,dxlo,dxhi,ans)
!nonlog10                     coltabc(iembx,iemby,ipairc(ihx,ihy))=max(0.,ans)
               coltabc(iembx,iemby,ipairc(ihx,ihy))=  &
                  -log10(max(1.e-30,ans))
            endif
            if (ipairr(ihx,ihy) .gt. 0) then
               call avint(dx,gx,ndx,dxlo,dxhi,ans)
!nonlog10                     coltabr(iembx,iemby,ipairr(ihx,ihy))=max(0.,ans)
               coltabr(iembx,iemby,ipairr(ihx,ihy))=  &
                  -log10(max(1.e-30,ans))
            endif
         enddo
      enddo
   endif
enddo
enddo
return
end subroutine mkcoltb

!******************************************************************************

real function xj(dx,cvx,pvx,cvy,pvy,vny,dnx,dny,xnu,ynu  &
                ,gyn1,gyn2,gynp,gynp1,gynp2)
implicit none
real :: dx,cvx,pvx,cvy,pvy,vny,dnx,dny,xnu,ynu,gyn1,gyn2,gynp,gynp1,gynp2  &
       ,dnxi,rdx,vx,dxy,ynup
real :: gammln,gammp,gammq
dnxi = 1. / dnx
rdx = dx * dnxi
vx = cvx * dx ** pvx
dxy = (vx / cvy) ** (1. / pvy) / dny
ynup = ynu + pvy

if (rdx .lt. 38.) then
   xj=exp(-rdx-gammln(xnu)-gammln(ynu))*rdx**(xnu-1.)*dnxi*(  &
       vx*(dx*dx*(gammp(ynu,dxy)-gammq(ynu,dxy))  &
         +2.*dx*dny*gyn1*(gammp(ynu+1.,dxy)-gammq(ynu+1.,dxy))  &
         +dny*dny*gyn2*(gammp(ynu+2.,dxy)-gammq(ynu+2.,dxy)))  &
     -vny*(dx*dx*gynp*(gammp(ynup,dxy)-gammq(ynup,dxy))  &
         +2.*dx*dny*gynp1*(gammp(ynup+1.,dxy)-gammq(ynup+1.,dxy))  &
         +dny*dny*gynp2*(gammp(ynup+2.,dxy)-gammq(ynup+2.,dxy))))
else
   xj = 0.
endif
return
end function

!******************************************************************************

subroutine make_autotab()

!!use micphys

implicit none

integer, parameter :: ithresh=14,ithresh1=17,ibins=36,icutoff=15
integer :: i,k,idcc,id1cd,id2cd,idccr,irrcr,idrcr,irrr,idrr,lcat,cld,maxcld
real :: r1,r2,r3,ri,en1,en2,en3,enice,en1i,en1i2,d1,d2,d3,di &
       ,sum1,sum10,sun1,sun10,sum2,sum20,sun2,sun20,sum3,sum30,sun3,sun30
real :: collectormass

real, dimension(ibins+1) :: x,diam,xcs,diamcs,xca,diamca,xcg,diamcg,xch,diamch
real, dimension(ibins+1,4) :: xi,diami
real, dimension(ibins) :: ank0,amk0,ank,amk,ank1,amk1,ank2,amk2,ank3,amk3
real, dimension(ibins,ibins,6) :: akbarx
real, dimension(ibins,ibins,4) :: akbarci,akbarxci
real, dimension(ibins,ibins) :: akbar,akbarcs,akbarca,akbarcg,akbarch

! This subroutine works in cgs units.
! read in mass grid x(k+1)=2*x(k), diameters (diam) and collection kernel kbar
! GNF: kernels for ice with cloud?

call mcphys_data(x,diam,akbar,ibins)
call data_cs(xcs,diamcs,akbarcs,ibins,cfmas(4),pwmas(4))
call data_ca(xca,diamca,akbarca,ibins,cfmas(5),pwmas(5))
call data_cg(xcg,diamcg,akbarcg,ibins,cfmas(6),pwmas(6))
call data_ch(xch,diamch,akbarch,ibins,cfmas(7),pwmas(7))

!call azero (ibins*ibins*6,akbarx)
!call azero (ibins*ibins*4,akbarci)
!call azero (ibins*ibins*4,akbarxci)
akbarx  =0.0
akbarci =0.0
akbarxci=0.0

do i=1,ibins+1
 xi(i,1)=xcs(i)
 xi(i,2)=xca(i)
 xi(i,3)=xcg(i)
 xi(i,4)=xch(i)
 diami(i,1)=diamcs(i)
 diami(i,2)=diamca(i)
 diami(i,3)=diamcg(i)
 diami(i,4)=diamch(i)
enddo
do i=1,ibins
 do k=1,ibins
  akbarci(i,k,1)=akbarcs(i,k)
  akbarci(i,k,2)=akbarca(i,k)
  akbarci(i,k,3)=akbarcg(i,k)
  akbarci(i,k,4)=akbarch(i,k)
 enddo
enddo

!Break up kernel in 6 different segments for cloud,drizzle,rain
do i=1,ibins
   do k=1,ibins
      !For CLOUD-CLOUD collection for 1 and 2 mode options
      if(i <  ithresh  .and. k <  ithresh  .and. jnmb(8).gt.0) &
         akbarx(i,k,1) = akbar(i,k)
      if(i <=  icutoff  .and. k <=  icutoff  .and. jnmb(8).eq.0) &
         akbarx(i,k,1) = akbar(i,k)

      !For RAIN-CLOUD collection
      if(i >  ithresh1 .and. k < ithresh  .and. jnmb(8).gt.0) &
         akbarx(i,k,3) = akbar(i,k)
      if(i >  icutoff  .and. k <= icutoff  .and. jnmb(8).eq.0) &
         akbarx(i,k,3) = akbar(i,k)

      if(jnmb(8).gt.0) then  !Dual-cloud mode option only
       !For DRIZZLE-CLOUD collection
       if(i >= ithresh  .and. i <= ithresh1 .and. k < ithresh) &
          akbarx(i,k,2) = akbar(i,k)

       !For DRIZZLE-DRIZZLE collection
       if(i<=ithresh1 .and. i>=ithresh .and. k<=ithresh1 .and. k>=ithresh) &
          akbarx(i,k,5) = akbar(i,k)

       !For RAIN-DRIZZLEcollection
       if(i >  ithresh1 .and. k <= ithresh1 .and. k >= ithresh) &
          akbarx(i,k,4) = akbar(i,k)
      endif

      !For RAIN-RAIN collection (not currently used)
      !if(i >  ithresh1 .and. k >  ithresh1) &
      !   akbarx(i,k,6) = akbar(i,k)
   enddo
enddo

!Break up kernel in 2 different segments for c1,c2,ice species
do lcat=4,7
 do i=1,ibins
   do k=1,ibins
      if(i > icutoff .and. k <= icutoff) &
         akbarxci(i,k,lcat-3) = akbarci(i,k,lcat-3)
   enddo
 enddo
enddo

!(d1min and d1max are equivalent to dmb0 and dmb1, but may have diff values)
!Diameters in cm (4.e-4 cm = 4 microns)
!Mixing ratios in g/cm3 (.01e-6g/cm3 ~ 1.e-11kg/kg for dn0=1)
d1min = 4.e-4
if(jnmb(8).eq.0) d1max = 50.e-4
if(jnmb(8).gt.0) d1max = 35.e-4
d2min = 65.e-4
d2max = 100.e-4
d3min = 1.e-2
d3max = 1.
r3min = .01e-6
r3max = 20.e-6
!For ICE SPECIES: SNOW, AGGREGATES, GRAUPEL, AND HAIL
dimin = 1.e-2
dimax = 1.
rimin = .01e-6
rimax = 20.e-6

d1ecr = log10 (d1max / d1min) / float(ndccr-1)
d2ecr = log10 (d2max / d2min) / float(ndccr-1)
r3ecr = log10 (r3max / r3min) / float(nrrcr-1)
r3err = log10 (r3max / r3min) / float(nrrr-1)
!For SNOW, AGGREGATES, GRAUPEL, AND HAIL
rieci = log10 (rimax / rimin) / float(nrrcr-1)

en1 = 100.
en2 = 1.
en3 = 1.e-6
!For SNOW, AGGREGATES, GRAUPEL, AND HAIL
enice = 1.e-6

en1i = 1. / en1
en1i2 = en1i ** 2.

!**************************************************************************
!**************************** CLOUD-CLOUD *********************************
!**************************************************************************
! Start 1 cc loop for r(1) and c(1,2) [r1 approx = (-r2)]
d2 = 65.e-4
r2 = en2 * .5236 * d2 ** 3
r3 = .01e-6

do idcc = 1,ndcc
   d1 = d1min + (d1max - d1min) * float(idcc-1) / float(ndcc-1)
   r1 = en1 * .5236 * d1 ** 3

   if(jnmb(8).gt.0) then
    call initg2mode(r1,r2,r3,en1,en2,en3,gnu(1),gnu(8),gnu(2),diam,x &
     ,amk0,ank0,ank1,amk1,ank2,amk2,ank3,amk3,ithresh,ithresh1,ibins)
    call sumn(ank0,amk0,1,ithresh-1,ibins,sun10,sum10)
    call sumn(ank0,amk0,ithresh,ithresh1,ibins,sun20,sum20)
    call sxy(x,amk0,ank0,amk,ank,akbarx(1,1,1),collectormass)
    call sumn(ank,amk,1,ithresh-1,ibins,sun1,sum1)
    call sumn(ank,amk,ithresh,ithresh1,ibins,sun2,sum2)
   endif
   if(jnmb(8).eq.0) then
    call initg1mode(r1,r3,en1,en3,gnu(1),gnu(2),diam,x &
     ,amk0,ank0,ank1,amk1,ank3,amk3,icutoff,ibins,2)
    call sumn(ank0,amk0,1,icutoff,ibins,sun10,sum10)
    call sumn(ank0,amk0,icutoff+1,ibins,ibins,sun20,sum20)
    call sxy(x,amk0,ank0,amk,ank,akbarx(1,1,1),collectormass)
    call sumn(ank,amk,1,icutoff,ibins,sun1,sum1)
    call sumn(ank,amk,icutoff+1,ibins,ibins,sun2,sum2)
   endif

   r1tabcc(idcc) = max(0.,max((sum10-sum1)*en1i2,(sum2-sum20)*en1i2))
   c1tabcc(idcc) = max(0.,(sun10-sun1)*en1i2)
   c2tabcc(idcc) = max(0.,(sun2-sun20)*en1i2)
   if(r1tabcc(idcc)==0.0 .or. c1tabcc(idcc)==0.0 .or. c2tabcc(idcc)==0.0) then
     r1tabcc(idcc) = 0.
     c1tabcc(idcc) = 0.
     c2tabcc(idcc) = 0.
   endif

!   write(*,'(a2,2X,f7.1,2X,i2,2X,e10.3,2X,e10.3,2X,e10.3)') &
!     'CC',d1*10000,idcc,r1tabcc(idcc),c1tabcc(idcc),c2tabcc(idcc)
enddo

!**************************************************************************
!**************************** RAIN-CLOUD **********************************
!**************************************************************************
! Start 3 cr loops for r(1) c(1,3)  [r1 approx = (-r3)]
d2 = 65.e-4
r2 = en2 * .5236 * d2 ** 3

do idccr = 1,ndccr
   d1 = d1min * 10. ** (d1ecr * float(idccr-1))
   r1 = en1 * .5236 * d1 ** 3

   do irrcr = 1,nrrcr
      r3 = r3min * 10. ** (r3ecr * float(irrcr-1))
      d3minx = max(d3min,(r3 / (.1 * .5236)) ** .333333)
      d3ecr = alog10(d3max / d3minx) / float(ndrcr-1)

      do idrcr = 1,ndrcr
       d3 = d3minx * 10. ** (d3ecr * float(idrcr-1))
       en3 = r3 / (.5236 * d3 ** 3)

       if(jnmb(8).gt.0) then
        call initg2mode(r1,r2,r3,en1,en2,en3,gnu(1),gnu(8),gnu(2),diam,x &
         ,amk0,ank0,ank1,amk1,ank2,amk2,ank3,amk3,ithresh,ithresh1,ibins)
        call sumn(ank0,amk0,1,ithresh-1,ibins,sun10,sum10)
        call sxy(x,amk0,ank0,amk,ank,akbarx(1,1,3),collectormass)
        call sumn(ank,amk,1,ithresh-1,ibins,sun1,sum1)
       endif
       if(jnmb(8).eq.0) then
        call initg1mode(r1,r3,en1,en3,gnu(1),gnu(2),diam,x &
         ,amk0,ank0,ank1,amk1,ank3,amk3,icutoff,ibins,2)
        call sumn(ank0,amk0,1,icutoff,ibins,sun10,sum10)
        call sxy(x,amk0,ank0,amk,ank,akbarx(1,1,3),collectormass)
        call sumn(ank,amk,1,icutoff,ibins,sun1,sum1)
       endif

       r1tabcr(idccr,irrcr,idrcr) = alog10(max(1.e-20,(sum10-sum1)*en1i))
       c1tabcr(idccr,irrcr,idrcr) = alog10(max(1.e-20,(sun10-sun1)*en1i))

!       write(*,'(a2,2X,i2,2X,i2,2X,i2,2X,f7.1,2X,f10.2,2X,f10.2)') &
!        'CR',idccr,irrcr,idrcr,d1*10000,r1tabcr(idccr,irrcr,idrcr) &
!        ,c1tabcr(idccr,irrcr,idrcr)
     enddo
   enddo
enddo

if(jnmb(8).gt.0) then
!**************************************************************************
!*************************** DRIZZLE-CLOUD ********************************
!**************************************************************************
r3 = .01e-6

do id1cd = 1,ndcd
   d1 = d1min + (d1max - d1min) * float(id1cd-1) / float(ndcd-1)
   r1 = en1 * .5236 * d1 ** 3

   do id2cd = 1,ndcd
      d2 = d2min + (d2max - d2min) * float(id2cd-1) / float(ndcd-1)
      r2 = en2 * .5236 * d2 ** 3

      call initg2mode(r1,r2,r3,en1,en2,en3,gnu(1),gnu(8),gnu(2),diam,x &
       ,amk0,ank0,ank1,amk1,ank2,amk2,ank3,amk3,ithresh,ithresh1,ibins)
      call sumn(ank0,amk0,1,ithresh-1,ibins,sun10,sum10)
      call sumn(ank0,amk0,ithresh,ithresh1,ibins,sun20,sum20)
      call sumn(ank0,amk0,ithresh1+1,ibins,ibins,sun30,sum30)
      call sxy(x,amk0,ank0,amk,ank,akbarx(1,1,2),collectormass)
      call sumn(ank,amk,1,ithresh-1,ibins,sun1,sum1)
      call sumn(ank,amk,ithresh,ithresh1,ibins,sun2,sum2)
      call sumn(ank,amk,ithresh1+1,ibins,ibins,sun3,sum3)

      r1tabcd(id1cd,id2cd) = (sum10-sum1)*en1i
      c1tabcd(id1cd,id2cd) = (sun10-sun1)*en1i
      r2tabcd(id1cd,id2cd) = (sum20-sum2)*en1i
      c2tabcd(id1cd,id2cd) = (sun20-sun2)*en1i

!      write(*,'(a3,2X,i2,2X,i2,2X,f6.1,2X,f6.1,6(2X,e10.3))') &
!       ,'CD',id1cd,id2cd,d1*10000,d2*10000 &
!       ,r1tabcd(id1cd,id2cd),r2tabcd(id1cd,id2cd),(sum30-sum3)*en1i &
!       ,c1tabcd(id1cd,id2cd),c2tabcd(id1cd,id2cd),(sun30-sun3)*en1i
   enddo
enddo

!**************************************************************************
!**************************** RAIN-DRIZZLE ********************************
!**************************************************************************
! Start 3 dr loops for r(2) c(2,3)   [r2 approx = (-r3)]
d1 = 4.e-4
r1 = en1 * .5236 * d1 ** 3

do idccr = 1,ndccr
   d2 = d2min * 10. ** (d2ecr * float(idccr-1))
   r2 = en2 * .5236 * d2 ** 3

   do irrcr = 1,nrrcr
      r3 = r3min * 10. ** (r3ecr * float(irrcr-1))
      d3minx = max(d3min,(r3 / (.1 * .5236)) ** .333333)
      d3ecr = alog10(d3max / d3minx) / float(ndrcr-1)

      do idrcr = 1,ndrcr
       d3 = d3minx * 10. ** (d3ecr * float(idrcr-1))
       en3 = r3 / (.5236 * d3 ** 3)

       call initg2mode(r1,r2,r3,en1,en2,en3,gnu(1),gnu(8),gnu(2),diam,x &
        ,amk0,ank0,ank1,amk1,ank2,amk2,ank3,amk3,ithresh,ithresh1,ibins)
       call sumn(ank0,amk0,ithresh,ithresh1,ibins,sun20,sum20)
       call sxy(x,amk0,ank0,amk,ank,akbarx(1,1,4),collectormass)
       call sumn(ank,amk,ithresh,ithresh1,ibins,sun2,sum2)

       r2tabcr(idccr,irrcr,idrcr) = alog10(max(1.e-20,sum20-sum2))
       c2tabcr(idccr,irrcr,idrcr) = alog10(max(1.e-20,sun20-sun2))

!       write(*,'(a3,2X,i2,2X,i2,2X,i2,2X,f10.2,2X,f10.2)') &
!        ,'DR',idccr,irrcr,idrcr,r2tabcr(idccr,irrcr,idrcr) &
!        ,c2tabcr(idccr,irrcr,idrcr)
      enddo
   enddo
enddo

!**************************************************************************
!**************************** DRIZZLE-DRIZZLE *****************************
!**************************************************************************
! Start 1 dd loop for r(2) and c(2,3) [r2=(-r3)]
d1 = 4.e-4
r1 = en1 * .5236 * d1 ** 3
r3 = .01e-6

do idcc = 1,ndcc
   d2 = d2min + (d2max - d2min) * float(idcc-1) / float(ndcc-1)
   r2 = en2 * .5236 * d2 ** 3

   call initg2mode(r1,r2,r3,en1,en2,en3,gnu(1),gnu(8),gnu(2),diam,x &
    ,amk0,ank0,ank1,amk1,ank2,amk2,ank3,amk3,ithresh,ithresh1,ibins)
   call sumn(ank0,amk0,ithresh,ithresh1,ibins,sun20,sum20)
   call sumn(ank0,amk0,ithresh1+1,ibins,ibins,sun30,sum30)
   call sxy(x,amk0,ank0,amk,ank,akbarx(1,1,5),collectormass)
   call sumn(ank,amk,ithresh,ithresh1,ibins,sun2,sum2)
   call sumn(ank,amk,ithresh1+1,ibins,ibins,sun3,sum3)

   r2tabdd(idcc) = max(0.,sum20-sum2)
   c2tabdd(idcc) = max(0.,sun20-sun2)
   c3tabdd(idcc) = max(0.,sun3-sun30)

!   write(*,'(a4,2X,f7.1,2X,i2,4(2X,e10.3))') ,'DD',d2*10000,idcc &
!    ,r2tabdd(idcc),sum30-sum3,c2tabdd(idcc),c3tabdd(idcc)
enddo

endif !IF CLOUD2 turned on

!**************************************************************************
!**************************** RAIN RAIN ***********************************
!**************************************************************************
! Start 2 rr loops for c(3)
!d1 = 4.e-4
!r1 = en1 * .5236 * d1 ** 3
!
!do irrr = 1,nrrr
!   r3 = r3min * 10. ** (r3err * float(irrr-1))
!   d3minx = max(d3min,(r3 / (.1 * .5236)) ** .333333)
!   d3err = alog10(d3max / d3minx) / float(ndrr-1)
!
!   do idrr = 1,ndrr
!      d3 = d3minx * 10. ** (d3err * float(idrr-1))
!      en3 = r3 / (.5236 * d3 ** 3)
!
!      call initg1mode(r1,r3,en1,en3,gnu(1),gnu(2),diam,x &
!       ,amk0,ank0,ank1,amk1,ank3,amk3,ithresh1,ibins,2)
!      call sumn(ank0,amk0,ithresh1+1,ibins,ibins,sun30,sum30)
!      call sxy(x,amk0,ank0,amk,ank,akbarx(1,1,6),collectormass)
!      call sumn(ank,amk,ithresh1+1,ibins,ibins,sun3,sum3)
!
!      c3tabrr(irrr,idrr) = alog10(max(1.e-25,sun30-sun3))
!      write(*,'(a2,2X,i2,2X,i2,2X,f10.2)') &
!         ,'RR',irrr,idrr,c3tabrr(irrr,idrr)
!   enddo
!enddo

!*************************************************************************
!*********************** CLOUD1/CLOUD2 ICE SPECIES ***********************
!*************************************************************************
maxcld=1
if(jnmb(8).gt.0) maxcld=2

do cld=1,maxcld
do lcat=4,7

do idccr = 1,ndccr
   if(cld==1)then
     d1 = d1min * 10. ** (d1ecr * float(idccr-1))
     r1 = en1 * .5236 * d1 ** 3
   elseif(cld==2)then
     d2 = d2min * 10. ** (d2ecr * float(idccr-1))
     r2 = en2 * .5236 * d2 ** 3
   endif

   do irrcr = 1,nrrcr
    ri = rimin * 10. ** (rieci * float(irrcr-1))
    !Convert ri (g to kg) & switch diminx (m to cm) after
    !Multiply by 200 instead of 100 to bump up the min diameter
    diminx=max(dimin,((ri/1000./cfmas(lcat))**(1./pwmas(lcat)))*200.)
    dieci = alog10(dimax / diminx) / float(ndrcr-1)

     do idrcr = 1,ndrcr
      di = diminx * 10. ** (dieci * float(idrcr-1))
      !Here ri has units of (g/cm3)
      enice = ri / (1000. * cfmas(lcat) * (di/100.) ** pwmas(lcat))

      if(cld==1) then
       call initg1mode(r1,ri,en1,enice,gnu(1),gnu(lcat) &
        ,diami(1,lcat-3),xi(1,lcat-3),amk0,ank0,ank1,amk1 &
        ,ank3,amk3,icutoff,ibins,lcat)
       call sumn(ank0,amk0,1,icutoff,ibins,sun10,sum10)
       call sxy(x,amk0,ank0,amk,ank,akbarxci(1,1,lcat-3),collectormass)
       call sumn(ank,amk,1,icutoff,ibins,sun1,sum1)
      endif
      if(cld==2) then
       call initg1mode(r2,ri,en2,enice,gnu(8),gnu(lcat) &
        ,diami(1,lcat-3),xi(1,lcat-3),amk0,ank0,ank2,amk2 &
        ,ank3,amk3,ithresh1,ibins,lcat)
       call sumn(ank0,amk0,1,ithresh1,ibins,sun10,sum10)
       call sxy(x,amk0,ank0,amk,ank,akbarxci(1,1,lcat-3),collectormass)
       call sumn(ank,amk,1,ithresh1,ibins,sun1,sum1)
      endif

   if(cld==1)then
    r1tabci(idccr,irrcr,idrcr,lcat-3) = alog10(max(1.e-20,(sum10-sum1)*en1i))
    c1tabci(idccr,irrcr,idrcr,lcat-3) = alog10(max(1.e-20,(sun10-sun1)*en1i))
    r1rimer(idccr,irrcr,idrcr,lcat-3) = alog10(max(1.e-20,collectormass*en1i))
!    write(*,'(a2,2X,a5,i2,3(2X,a6,i2),3(2X,f10.2))') &
!     'CI','lcat=',lcat,'idccr=',idccr,'irrcr=',irrcr,'idrcr=',idrcr &
!     ,r1tabci(idccr,irrcr,idrcr,lcat-3) &
!     ,c1tabci(idccr,irrcr,idrcr,lcat-3),r1rimer(idccr,irrcr,idrcr,lcat-3)
   endif

   if(cld==2)then
    r2tabci(idccr,irrcr,idrcr,lcat-3) = alog10(max(1.e-20,sum10-sum1))
    c2tabci(idccr,irrcr,idrcr,lcat-3) = alog10(max(1.e-20,sun10-sun1))
    r2rimer(idccr,irrcr,idrcr,lcat-3) = alog10(max(1.e-20,collectormass))
!    write(*,'(a3,2X,a5,i2,3(2X,a6,i2),3(2X,f10.2))') &
!     'DI','lcat=',lcat,'idccr=',idccr,'irrcr=',irrcr,'idrcr=',idrcr &
!     ,r2tabci(idccr,irrcr,idrcr,lcat-3) &
!     ,c2tabci(idccr,irrcr,idrcr,lcat-3),r2rimer(idccr,irrcr,idrcr,lcat-3)
   endif

   enddo
 enddo
enddo

!*********************************
enddo !for looping lcat=4,7
enddo !for looping maxcld=1,2

return
end subroutine make_autotab

!******************************************************************************

subroutine sxy(x,amkd,ankd,amk,ank,akbar,smcollector)
implicit none
integer, parameter :: ibins=36
integer :: i,ik,k,l
real, dimension(ibins+1) :: x
real, dimension(ibins,ibins) :: akbar
real, dimension(ibins) :: xave,ankd,ank,amkd,amk,am2,am3,am4,psi,f
real :: ap,pi,dm,dn,sm1,sm2,sm3,sm4,sm5,sn1,sn2,sn3,sn4,dm4,dm2
real :: dmcollector,smcollector

data ap/1.062500000/
data pi/3.141592654/

do l=1,ibins
   if(ankd(l).gt.0)then
      xave(l)=amkd(l)/ankd(l)
   else
      xave(l)=0.
   endif
enddo

smcollector=0.

do k=1,ibins

! calculation of the 2nd, 3rd, and 4th moments of the mass distribution
! based on equation 8 in reference.

   am2(k)=ap*xave(k)*amkd(k)
   am3(k)=ap*ap*xave(k)*am2(k)
   am4(k)=ap*ap*ap*xave(k)*am3(k)

! these functions come out of the linear approximations used to integrate
! over partial bins.  they are defined:
!      psi(k) = nk(k+1)
!        f(k) = nk(k)
! where nk is the distribution function.  see equation 13 in reference.

   psi(k)=2./x(k)*(amkd(k)/x(k)-ankd(k))
   f(k)=2./x(k)*(2.*ankd(k)-amkd(k)/x(k))

! zeroing the tendencies on the moments.

   sm1=0.
   sm2=0.
   sm3=0.
   sm4=0.
   sm5=0.
   sn1=0.
   sn2=0.
   sn3=0.
   sn4=0.

! calculation of tendencies on moments

   do i=k,ibins
      dm=akbar(i,k)*(am2(k)*ankd(i)+amkd(k)*amkd(i))
      dn=akbar(i,k)*(ankd(k)*amkd(i)+amkd(k)*ankd(i))

      sm5=sm5+dm
      sn4=sn4+dn
   enddo

   if(k.gt.1)then
      sm3=akbar(k-1,k-1)*(am2(k-1)*ankd(k-1)+amkd(k-1)**2)
      sn2=akbar(k-1,k-1)*ankd(k-1)*amkd(k-1)
      dn=sn2
      dm=sm3
   endif

   do i=1,k-1
      dm4=akbar(k,i)*(ankd(k)*am2(i)+amkd(k)*amkd(i))
      sm4=sm4+dm4


      if(xave(k).ge.x(k))then
         dm2=akbar(k,i)*(4.*x(k)**2*psi(k)*amkd(i)  &
            +0.5*x(k)*(4.*psi(k)+f(k))*am2(i)  &
            -(psi(k)-f(k))*am3(i)  &
            -0.5/(x(k))*(psi(k)-f(k))*am4(i))
         sm2=sm2+dm2
         dn=akbar(k,i)*(2.*x(k)*psi(k)*amkd(i)  &
            +0.5*f(k)*am2(i)  &
            -0.5/(x(k))*(psi(k)-f(k))*am3(i))
         sn3=sn3+dn

      endif
   enddo

   do i=1,k-2
      ik=k-1
      if(xave(ik).ge.x(ik))then

         dm=akbar(ik,i)*(4.*x(ik)**2*psi(ik)*amkd(i)  &
            +x(ik)/2.*(4.*psi(ik)+f(ik))*am2(i)  &
            -(psi(ik)-f(ik))*am3(i)  &
            -0.5/(x(ik))*(psi(ik)-f(ik))*am4(i))
         sm1=sm1+dm
         dn=akbar(ik,i)*(2.*x(ik)*psi(ik)*amkd(i)  &
            +0.5*f(ik)*am2(i)  &
            -0.5/(x(ik))*(psi(ik)-f(ik))*am3(i))
         sn1=sn1+dn

      endif
   enddo

   amk(k)=amkd(k)+sm1-sm2+sm3+sm4-sm5
   ank(k)=ankd(k)+sn1+sn2-sn3-sn4

   !This is not necessary for "sxy" but does provide the amount of
   !mixing ratio of the collector species that participates in collection
   !in the binned process. This pulls out the "dm2" and "sm2" terms.
   !This sums the total collector mass over all bins; for each bin, it
   !calculates the mass involved in collecting cloud water in all bins
   !smaller than the collector bin. Currently this is only applicable to
   !the ice species as collectors and is sensitive to the collection kernels.
   do i=1,k-1
      dmcollector = akbar(k,i)*(4.*x(k)**2*psi(k)*amkd(i)  &
            +0.5*x(k)*(4.*psi(k)+f(k))*am2(i)  &
            -(psi(k)-f(k))*am3(i)  &
            -0.5/(x(k))*(psi(k)-f(k))*am4(i))
      smcollector = smcollector + dmcollector
   enddo

enddo

return
end subroutine sxy

!******************************************************************************
subroutine mcphys_data(x,diam,akbar,ibins)

implicit none
integer :: l,i,j,ibins,n,kernel
real :: pi
real, dimension(ibins+1) :: x,diam
real, dimension(ibins,ibins) :: akbar
real, dimension(36,36) :: aabar,ahbar

kernel=1

if(kernel==1) then  !Long's collection kernel
data (aabar( 1,n),n=1, 1) /-.47757E-01  /
data (aabar( 2,n),n=1, 2) /-.26460E+00,-.47965E-01 /
data (aabar( 3,n),n=1, 3) /-.82258E+00,-.26760E+00,-.20453E-01 /
data (aabar( 4,n),n=1, 4) /-.19050E+01,-.82072E+00,-.11992E+00, .78909E-01 /
data (aabar( 5,n),n=1, 5) /-.39171E+01,-.18915E+01,-.33270E+00, .41936E+00  &
   ,.34801E+00 /
data (aabar( 6,n),n=1, 6) /-.76415E+01,-.38808E+01,-.73737E+00, .14121E+01  &
   ,.18851E+01,.99793E+00 /
data (aabar( 7,n),n=1, 7) /-.14595E+02,-.75638E+01,-.14861E+01, .33598E+01  &
   ,.61219E+01, .54314E+01, .24751E+01 /
data (aabar( 8,n),n=1, 8) /-.27720E+02,-.14442E+02,-.28741E+01, .69895E+01  &
   ,.14394E+02, .17479E+02, .13500E+02, .57110E+01 /
data (aabar( 9,n),n=1, 9) /-.52737E+02,-.27428E+02,-.54729E+01, .13703E+02  &
   ,.29792E+02, .40971E+02, .43267E+02, .31185E+02, .12630E+02 /
data (aabar(10,n),n=1,10) /-.10083E+03,-.52188E+02,-.10391E+02, .26218E+02  &
   ,.58283E+02, .84686E+02, .10128E+03, .99726E+02, .69014E+02, .27176E+02 /
data (aabar(11,n),n=1,11) /-.19396E+03,-.99799E+02,-.19790E+02, .49801E+02  &
   ,.11143E+03, .16558E+03, .20922E+03, .23326E+03, .22039E+03, .14858E+03  &
   ,.57396E+02 /
data (aabar(12,n),n=1,12) /-.37536E+03,-.19200E+03,-.37896E+02, .94692E+02  &
   ,.21165E+03, .31650E+03, .40896E+03, .48169E+03, .51524E+03, .47402E+03  &
   ,.31389E+03, .11962E+03 /
data (aabar(13,n),n=1,13) /-.73047E+03,-.37164E+03,-.73015E+02, .18089E+03  &
   ,.40253E+03, .60115E+03, .78166E+03, .94143E+03, .10638E+04, .11078E+04  &
   ,.10008E+04, .65436E+03, .24691E+03 /
data (aabar(14,n),n=1,14) /-.14285E+04,-.72333E+03,-.14152E+03, .34764E+03  &
   ,.76925E+03, .11434E+04, .14846E+04, .17993E+04, .20789E+04, .22870E+04  &
   ,.23385E+04, .20854E+04, .13509E+04, .50600E+03 /
data (aabar(15,n),n=1,15) /-.41365E+04,-.20869E+04,-.40697E+03, .99310E+03  &
   ,.21878E+04, .32394E+04, .41995E+04, .51084E+04, .59888E+04, .68297E+04  &
   ,.75528E+04, .79583E+04, .76785E+04, .62489E+04, .76776E+03 /
data (aabar(16,n),n=1,16) / .63760E+04, .64739E+04, .65970E+04, .67516E+04  &
   ,.69451E+04, .71861E+04, .74835E+04, .78448E+04, .82709E+04, .87453E+04  &
   ,.92111E+04, .95276E+04, .94079E+04, .83797E+04, .26045E+04, .89777E+03 /
data (aabar(17,n),n=1,17) / .62974E+04, .63746E+04, .64717E+04, .65934E+04  &
   ,.67457E+04, .69355E+04, .71702E+04, .74571E+04, .78005E+04, .81957E+04  &
   ,.86163E+04, .89879E+04, .91399E+04, .87394E+04, .46530E+04, .26045E+04  &
   ,.89777E+03 /
data (aabar(18,n),n=1,18) / .62353E+04, .62963E+04, .63729E+04, .64689E+04  &
   ,.65889E+04, .67383E+04, .69233E+04, .71502E+04, .74238E+04, .77446E+04  &
   ,.81009E+04, .84538E+04, .87067E+04, .86514E+04, .59471E+04, .46530E+04  &
   ,.26045E+04, .89777E+03 /
data (aabar(19,n),n=1,19) / .61862E+04, .62344E+04, .62949E+04, .63707E+04  &
   ,.64653E+04, .65831E+04, .67290E+04, .69080E+04, .71250E+04, .73819E+04  &
   ,.76742E+04, .79815E+04, .82491E+04, .83524E+04, .66125E+04, .59471E+04  &
   ,.46530E+04, .26045E+04, .89777E+03 /
data (aabar(20,n),n=1,20) / .61474E+04, .61855E+04, .62334E+04, .62932E+04  &
   ,.63679E+04, .64608E+04, .65759E+04, .67172E+04, .68887E+04, .70932E+04  &
   ,.73291E+04, .75856E+04, .78311E+04, .79911E+04, .68735E+04, .66125E+04  &
   ,.59471E+04, .46530E+04, .26045E+04, .89777E+03 /
data (aabar(21,n),n=1,21) / .61166E+04, .61468E+04, .61847E+04, .62320E+04  &
   ,.62910E+04, .63644E+04, .64552E+04, .65668E+04, .67023E+04, .68644E+04  &
   ,.70531E+04, .72625E+04, .74738E+04, .76415E+04, .69140E+04, .68735E+04  &
   ,.66125E+04, .59471E+04, .46530E+04, .26045E+04, .89777E+03 /
data (aabar(22,n),n=1,22) / .60923E+04, .61162E+04, .61462E+04, .61836E+04  &
   ,.62303E+04, .62883E+04, .63600E+04, .64481E+04, .65553E+04, .66836E+04  &
   ,.68338E+04, .70027E+04, .71786E+04, .73330E+04, .68498E+04, .69140E+04  &
   ,.68735E+04, .66125E+04, .59471E+04, .46530E+04, .26045E+04, .89777E+03 /
data (aabar(23,n),n=1,23) / .60730E+04, .60919E+04, .61157E+04, .61453E+04  &
   ,.61823E+04, .62281E+04, .62848E+04, .63545E+04, .64392E+04, .65408E+04  &
   ,.66601E+04, .67953E+04, .69391E+04, .70729E+04, .67447E+04, .68498E+04  &
   ,.69140E+04, .68735E+04, .66125E+04, .59471E+04, .46530E+04, .26045E+04  &
   ,.89777E+03 /
data (aabar(24,n),n=1,24) / .60577E+04, .60727E+04, .60915E+04, .61150E+04  &
   ,.61443E+04, .61806E+04, .62254E+04, .62805E+04, .63475E+04, .64279E+04  &
   ,.65225E+04, .66304E+04, .67467E+04, .68590E+04, .66311E+04, .67447E+04  &
   ,.68498E+04, .69140E+04, .68735E+04, .66125E+04, .59471E+04, .46530E+04  &
   ,.26045E+04, .89777E+03 /
data (aabar(25,n),n=1,25) / .77967E+04, .78122E+04, .78316E+04, .78560E+04  &
   ,.78863E+04, .79242E+04, .79713E+04, .80294E+04, .81008E+04, .81878E+04  &
   ,.82924E+04, .84158E+04, .85571E+04, .87104E+04, .86265E+04, .88325E+04  &
   ,.90719E+04, .93363E+04, .95996E+04, .98007E+04, .98157E+04, .94274E+04  &
   ,.83361E+04, .63023E+04, .57988E+03 /
data (aabar(26,n),n=1,26) / .69349E+04, .69458E+04, .69595E+04, .69766E+04  &
   ,.69979E+04, .70244E+04, .70573E+04, .70978E+04, .71473E+04, .72072E+04  &
   ,.72788E+04, .73623E+04, .74565E+04, .75566E+04, .74715E+04, .76064E+04  &
   ,.77647E+04, .79435E+04, .81311E+04, .82983E+04, .83827E+04, .82640E+04  &
   ,.77406E+04, .65488E+04, .15807E+04, .51662E+03 /
data (aabar(27,n),n=1,27) / .61704E+04, .61781E+04, .61877E+04, .61997E+04  &
   ,.62147E+04, .62333E+04, .62562E+04, .62843E+04, .63186E+04, .63598E+04  &
   ,.64086E+04, .64648E+04, .65271E+04, .65912E+04, .65100E+04, .65961E+04  &
   ,.66976E+04, .68135E+04, .69382E+04, .70574E+04, .71390E+04, .71194E+04  &
   ,.68816E+04, .62379E+04, .26526E+04, .14083E+04, .46025E+03 /
data (aabar(28,n),n=1,28) / .54916E+04, .54971E+04, .55038E+04, .55123E+04  &
   ,.55228E+04, .55357E+04, .55517E+04, .55712E+04, .55949E+04, .56232E+04  &
   ,.56562E+04, .56938E+04, .57345E+04, .57747E+04, .57001E+04, .57533E+04  &
   ,.58161E+04, .58880E+04, .59661E+04, .60426E+04, .61008E+04, .61062E+04  &
   ,.59940E+04, .56500E+04, .31742E+04, .23632E+04, .12546E+04, .41004E+03 /
data (aabar(29,n),n=1,29) / .48886E+04, .48924E+04, .48971E+04, .49031E+04  &
   ,.49104E+04, .49195E+04, .49306E+04, .49441E+04, .49604E+04, .49797E+04  &
   ,.50020E+04, .50269E+04, .50530E+04, .50774E+04, .50108E+04, .50422E+04  &
   ,.50792E+04, .51212E+04, .51667E+04, .52111E+04, .52447E+04, .52486E+04  &
   ,.51861E+04, .49913E+04, .32935E+04, .28279E+04, .21054E+04, .11177E+04  &
   ,.36530E+03 /
data (aabar(30,n),n=1,30) / .43524E+04, .43551E+04, .43585E+04, .43626E+04  &
   ,.43678E+04, .43741E+04, .43818E+04, .43912E+04, .44024E+04, .44155E+04  &
   ,.44304E+04, .44467E+04, .44631E+04, .44771E+04, .44188E+04, .44361E+04  &
   ,.44561E+04, .44786E+04, .45022E+04, .45241E+04, .45384E+04, .45339E+04  &
   ,.44893E+04, .43663E+04, .31847E+04, .29342E+04, .25193E+04, .18757E+04  &
   ,.99579E+03, .32545E+03 /
data (aabar(31,n),n=1,31) / .38756E+04, .38775E+04, .38799E+04, .38828E+04  &
   ,.38864E+04, .38908E+04, .38961E+04, .39026E+04, .39102E+04, .39191E+04  &
   ,.39290E+04, .39395E+04, .39494E+04, .39568E+04, .39066E+04, .39149E+04  &
   ,.39241E+04, .39340E+04, .39435E+04, .39507E+04, .39516E+04, .39392E+04  &
   ,.39006E+04, .38129E+04, .29707E+04, .28372E+04, .26141E+04, .22445E+04  &
   ,.16710E+04, .88715E+03, .28994E+03 /
data (aabar(32,n),n=1,32) / .30106E+04, .30118E+04, .30132E+04, .30149E+04  &
   ,.30171E+04, .30197E+04, .30229E+04, .30266E+04, .30309E+04, .30357E+04  &
   ,.30408E+04, .30456E+04, .30491E+04, .30494E+04, .30032E+04, .30013E+04  &
   ,.29981E+04, .29929E+04, .29844E+04, .29706E+04, .29480E+04, .29111E+04  &
   ,.28504E+04, .27503E+04, .20956E+04, .19717E+04, .17926E+04, .15245E+04  &
   ,.11243E+04, .55892E+03, .22135E+03, .00000E+00 /
data (aabar(33,n),n=1,33) / .23888E+04, .23895E+04, .23903E+04, .23914E+04  &
   ,.23927E+04, .23943E+04, .23962E+04, .23983E+04, .24007E+04, .24033E+04  &
   ,.24057E+04, .24075E+04, .24077E+04, .24049E+04, .23645E+04, .23582E+04  &
   ,.23497E+04, .23382E+04, .23225E+04, .23007E+04, .22699E+04, .22258E+04  &
   ,.21613E+04, .20655E+04, .15572E+04, .14495E+04, .13057E+04, .11057E+04  &
   ,.82055E+03, .41732E+03, .17636E+03, .00000E+00, .00000E+00 /
data (aabar(34,n),n=1,34) / .18955E+04, .18959E+04, .18964E+04, .18971E+04  &
   ,.18979E+04, .18988E+04, .18999E+04, .19011E+04, .19024E+04, .19036E+04  &
   ,.19045E+04, .19047E+04, .19033E+04, .18990E+04, .18647E+04, .18567E+04  &
   ,.18462E+04, .18326E+04, .18145E+04, .17905E+04, .17581E+04, .17138E+04  &
   ,.16525E+04, .15662E+04, .11695E+04, .10771E+04, .95995E+03, .80547E+03  &
   ,.59516E+03, .30433E+03, .13287E+03, .00000E+00, .00000E+00, .00000E+00 /
data (aabar(35,n),n=1,35) / .15041E+04, .15044E+04, .15047E+04, .15051E+04  &
   ,.15056E+04, .15061E+04, .15067E+04, .15074E+04, .15080E+04, .15084E+04  &
   ,.15085E+04, .15079E+04, .15058E+04, .15011E+04, .14725E+04, .14642E+04  &
   ,.14536E+04, .14399E+04, .14221E+04, .13988E+04, .13682E+04, .13274E+04  &
   ,.12725E+04, .11976E+04, .88682E+03, .80897E+03, .71340E+03, .59221E+03  &
   ,.43360E+03, .22074E+03, .97241E+02, .00000E+00, .00000E+00, .00000E+00  &
   ,.00000E+00 /
data (aabar(36,n),n=1,36) / .11936E+04, .11938E+04, .11940E+04, .11942E+04  &
   ,.11945E+04, .11948E+04, .11951E+04, .11954E+04, .11957E+04, .11957E+04  &
   ,.11954E+04, .11943E+04, .11921E+04, .11876E+04, .11640E+04, .11563E+04  &
   ,.11464E+04, .11337E+04, .11174E+04, .10963E+04, .10689E+04, .10330E+04  &
   ,.98554E+03, .92214E+03, .67808E+03, .61344E+03, .53582E+03, .44015E+03  &
   ,.31885E+03, .16089E+03, .70536E+02, .00000E+00, .00000E+00, .00000E+00  &
   ,.00000E+00, .00000E+00 /
endif

if(kernel==2) then  !Hall's collection kernel
data (ahbar( 1,n),n=1, 1)/ &
0.19826E+00/
data (ahbar( 2,n),n=1, 2)/ &
0.11115E+01, 0.31472E+00/
data (ahbar( 3,n),n=1, 3)/ &
0.29798E+01, 0.17644E+01, 0.49959E+00/
data (ahbar( 4,n),n=1, 4)/ &
0.48803E+01, 0.47301E+01, 0.28009E+01, 0.79305E+00/
data (ahbar( 5,n),n=1, 5)/ &
0.64832E+01, 0.77470E+01, 0.75086E+01, 0.44461E+01, 0.12589E+01/
data (ahbar( 6,n),n=1, 6)/ &
0.75069E+01, 0.10291E+02, 0.12298E+02, 0.11919E+02, 0.70577E+01, 0.19984E+01/
data (ahbar( 7,n),n=1, 7)/ &
0.97323E+01, 0.11917E+02, 0.16337E+02, 0.19521E+02, 0.18921E+02, 0.11203E+02,&
0.31722E+01/
data (ahbar( 8,n),n=1, 8)/ &
0.11705E+02, 0.15449E+02, 0.18916E+02, 0.25933E+02, 0.30988E+02, 0.30034E+02,&
0.17784E+02, 0.50356E+01/
data (ahbar( 9,n),n=1, 9)/ &
0.28432E+01, 0.18581E+02, 0.24524E+02, 0.30028E+02, 0.41166E+02, 0.49190E+02,&
0.47677E+02, 0.28231E+02, 0.79935E+01/
data (ahbar(10,n),n=1,10)/ &
0.11264E+01, 0.54577E+01, 0.24506E+02, 0.33433E+02, 0.44838E+02, 0.64398E+02,&
0.80924E+02, 0.77938E+02, 0.41342E+02, 0.10151E+02/
data (ahbar(11,n),n=1,11)/ &
0.19148E+00, 0.30048E+01, 0.99008E+01, 0.32365E+02, 0.45873E+02, 0.67472E+02,&
0.10098E+03, 0.13218E+03, 0.12667E+03, 0.61080E+02, 0.12790E+02/
data (ahbar(12,n),n=1,12)/ &
0.23333E+00, 0.24125E+00, 0.63028E+01, 0.17275E+02, 0.43141E+02, 0.63748E+02,&
0.10244E+03, 0.15873E+03, 0.21451E+03, 0.20480E+03, 0.91229E+02, 0.16114E+02/
data (ahbar(13,n),n=1,13)/ &
0.96499E+00, 0.18995E+01, 0.31461E+01, 0.21422E+02, 0.49047E+02, 0.95732E+02,&
0.19602E+03, 0.44815E+03, 0.86833E+03, 0.11774E+04, 0.10705E+04, 0.54659E+03,&
0.20565E+03/
data (ahbar(14,n),n=1,14)/ &
0.82031E+00, 0.11610E+02, 0.26319E+02, 0.45949E+02, 0.14169E+03, 0.27697E+03,&
0.47057E+03, 0.76235E+03, 0.13258E+04, 0.22194E+04, 0.28021E+04, 0.24852E+04,&
0.13356E+04, 0.61752E+03/
data (ahbar(15,n),n=1,15)/ &
0.49547E+01, 0.63623E+01, 0.11954E+03, 0.27411E+03, 0.48091E+03, 0.12745E+04,&
0.23941E+04, 0.39885E+04, 0.53365E+04, 0.64055E+04, 0.75273E+04, 0.78726E+04,&
0.71201E+04, 0.55775E+04, 0.12227E+04/
data (ahbar(16,n),n=1,16)/ &
0.19925E+02, 0.25489E+02, 0.32725E+02, 0.70048E+03, 0.16119E+04, 0.28297E+04,&
0.37418E+04, 0.48244E+04, 0.57636E+04, 0.68214E+04, 0.77280E+04, 0.83723E+04,&
0.84671E+04, 0.74915E+04, 0.23147E+04, 0.20493E+04/
data (ahbar(17,n),n=1,17)/ &
0.29658E+03, 0.37826E+03, 0.48383E+03, 0.62105E+03, 0.12813E+04, 0.21736E+04,&
0.33625E+04, 0.42800E+04, 0.53873E+04, 0.64521E+04, 0.73410E+04, 0.79541E+04,&
0.83280E+04, 0.80576E+04, 0.42900E+04, 0.25057E+04, 0.29399E+04/
data (ahbar(18,n),n=1,18)/ &
0.72103E+03, 0.91733E+03, 0.11698E+04, 0.14961E+04, 0.19199E+04, 0.25718E+04,&
0.34340E+04, 0.45765E+04, 0.53658E+04, 0.62564E+04, 0.68926E+04, 0.76242E+04,&
0.81438E+04, 0.82870E+04, 0.57426E+04, 0.44930E+04, 0.25977E+04, 0.27291E+04/
data (ahbar(19,n),n=1,19)/ &
0.96659E+03, 0.12273E+04, 0.15613E+04, 0.19908E+04, 0.25456E+04, 0.32656E+04,&
0.38523E+04, 0.46057E+04, 0.55962E+04, 0.62746E+04, 0.69491E+04, 0.70313E+04,&
0.78366E+04, 0.83524E+04, 0.66125E+04, 0.59471E+04, 0.46530E+04, 0.26045E+04,&
0.89099E+03/
data (ahbar(20,n),n=1,20)/ &
0.97638E+03, 0.12378E+04, 0.15716E+04, 0.19991E+04, 0.25486E+04, 0.32579E+04,&
0.41778E+04, 0.46749E+04, 0.52921E+04, 0.60950E+04, 0.66107E+04, 0.71113E+04,&
0.71824E+04, 0.77993E+04, 0.68735E+04, 0.66125E+04, 0.59471E+04, 0.46530E+04,&
0.26045E+04, 0.89099E+03/
data (ahbar(21,n),n=1,21)/ &
0.94823E+03, 0.12006E+04, 0.15220E+04, 0.19322E+04, 0.24575E+04, 0.31324E+04,&
0.40029E+04, 0.51305E+04, 0.55255E+04, 0.59887E+04, 0.65800E+04, 0.69248E+04,&
0.72481E+04, 0.73294E+04, 0.69140E+04, 0.68735E+04, 0.66125E+04, 0.59471E+04,&
0.46530E+04, 0.26045E+04, 0.89099E+03/
data (ahbar(22,n),n=1,22)/ &
0.82817E+03, 0.10475E+04, 0.13263E+04, 0.16812E+04, 0.21341E+04, 0.27139E+04,&
0.34583E+04, 0.44175E+04, 0.56582E+04, 0.59651E+04, 0.62990E+04, 0.67126E+04,&
0.69632E+04, 0.72083E+04, 0.68455E+04, 0.69140E+04, 0.68735E+04, 0.66125E+04,&
0.59471E+04, 0.46530E+04, 0.26045E+04, 0.89099E+03/
data (ahbar(23,n),n=1,23)/ &
0.69439E+03, 0.87760E+03, 0.11100E+04, 0.14053E+04, 0.17812E+04, 0.22609E+04,&
0.28745E+04, 0.36617E+04, 0.46749E+04, 0.59830E+04, 0.62384E+04, 0.64949E+04,&
0.67994E+04, 0.69710E+04, 0.66897E+04, 0.68477E+04, 0.69140E+04, 0.68735E+04,&
0.66125E+04, 0.59471E+04, 0.46530E+04, 0.26045E+04, 0.89099E+03/
data (ahbar(24,n),n=1,24)/ &
0.57838E+03, 0.73052E+03, 0.92325E+03, 0.11677E+04, 0.14782E+04, 0.18735E+04,&
0.23776E+04, 0.30221E+04, 0.38482E+04, 0.49098E+04, 0.62770E+04, 0.64812E+04,&
0.66606E+04, 0.68557E+04, 0.66311E+04, 0.67447E+04, 0.68498E+04, 0.69140E+04,&
0.68735E+04, 0.66125E+04, 0.59471E+04, 0.46530E+04, 0.26045E+04, 0.89099E+03/
data (ahbar(25,n),n=1,25)/ &
0.59084E+03, 0.74589E+03, 0.94211E+03, 0.11907E+04, 0.15060E+04, 0.19065E+04,&
0.24163E+04, 0.30665E+04, 0.38980E+04, 0.49639E+04, 0.63340E+04, 0.80991E+04,&
0.83645E+04, 0.85993E+04, 0.86224E+04, 0.88325E+04, 0.90719E+04, 0.93363E+04,&
0.95996E+04, 0.98007E+04, 0.98157E+04, 0.94274E+04, 0.83361E+04, 0.63023E+04,&
0.57550E+03/
data (ahbar(26,n),n=1,26)/ &
0.41712E+03, 0.52636E+03, 0.66448E+03, 0.83925E+03, 0.10606E+04, 0.13414E+04,&
0.16979E+04, 0.21515E+04, 0.27296E+04, 0.34680E+04, 0.44128E+04, 0.56235E+04,&
0.71759E+04, 0.73865E+04, 0.73763E+04, 0.76028E+04, 0.77647E+04, 0.79435E+04,&
0.81311E+04, 0.82983E+04, 0.83827E+04, 0.82640E+04, 0.77406E+04, 0.65488E+04,&
0.15807E+04, 0.51271E+03/
data (ahbar(27,n),n=1,27)/ &
0.29457E+03, 0.37160E+03, 0.46891E+03, 0.59194E+03, 0.74760E+03, 0.94473E+03,&
0.11947E+04, 0.15119E+04, 0.19153E+04, 0.24289E+04, 0.30837E+04, 0.39193E+04,&
0.49856E+04, 0.63431E+04, 0.63635E+04, 0.65119E+04, 0.66944E+04, 0.68135E+04,&
0.69382E+04, 0.70574E+04, 0.71390E+04, 0.71194E+04, 0.68816E+04, 0.62379E+04,&
0.26526E+04, 0.14083E+04, 0.45678E+03/
data (ahbar(28,n),n=1,28)/ &
0.20808E+03, 0.26243E+03, 0.33104E+03, 0.41773E+03, 0.52730E+03, 0.66592E+03,&
0.84143E+03, 0.10639E+04, 0.13461E+04, 0.17045E+04, 0.21602E+04, 0.27398E+04,&
0.34765E+04, 0.44108E+04, 0.54855E+04, 0.56238E+04, 0.57419E+04, 0.58852E+04,&
0.59661E+04, 0.60426E+04, 0.61008E+04, 0.61062E+04, 0.59940E+04, 0.56500E+04,&
0.31742E+04, 0.23632E+04, 0.12546E+04, 0.40694E+03/
data (ahbar(29,n),n=1,29)/ &
0.14702E+03, 0.18538E+03, 0.23379E+03, 0.29491E+03, 0.37212E+03, 0.46970E+03,&
0.59313E+03, 0.74934E+03, 0.94723E+03, 0.11981E+04, 0.15162E+04, 0.19199E+04,&
0.24314E+04, 0.30782E+04, 0.38273E+04, 0.48524E+04, 0.49649E+04, 0.50559E+04,&
0.51642E+04, 0.52111E+04, 0.52447E+04, 0.52486E+04, 0.51861E+04, 0.49913E+04,&
0.32935E+04, 0.28279E+04, 0.21054E+04, 0.11177E+04, 0.36254E+03/
data (ahbar(30,n),n=1,30)/ &
0.10389E+03, 0.13098E+03, 0.16515E+03, 0.20827E+03, 0.26271E+03, 0.33147E+03,&
0.41837E+03, 0.52824E+03, 0.66723E+03, 0.84316E+03, 0.10659E+04, 0.13479E+04,&
0.17045E+04, 0.21543E+04, 0.26789E+04, 0.33884E+04, 0.42884E+04, 0.43778E+04,&
0.44448E+04, 0.45220E+04, 0.45384E+04, 0.45339E+04, 0.44893E+04, 0.43663E+04,&
0.31847E+04, 0.29342E+04, 0.25193E+04, 0.18757E+04, 0.99579E+03, 0.32299E+03/
data (ahbar(31,n),n=1,31)/ &
0.73425E+02, 0.92555E+02, 0.11668E+03, 0.14712E+03, 0.18553E+03, 0.23402E+03,&
0.29525E+03, 0.37261E+03, 0.47038E+03, 0.59398E+03, 0.75026E+03, 0.94780E+03,&
0.11972E+04, 0.15111E+04, 0.18798E+04, 0.23734E+04, 0.29973E+04, 0.37859E+04,&
0.38548E+04, 0.39003E+04, 0.39497E+04, 0.39392E+04, 0.39006E+04, 0.38129E+04,&
0.29707E+04, 0.28372E+04, 0.26141E+04, 0.22445E+04, 0.16710E+04, 0.88715E+03,&
0.28775E+03/
data (ahbar(32,n),n=1,32)/ &
0.45270E+02, 0.57059E+02, 0.71923E+02, 0.90671E+02, 0.11432E+03, 0.14416E+03,&
0.18182E+03, 0.22936E+03, 0.28938E+03, 0.36518E+03, 0.46087E+03, 0.58157E+03,&
0.73358E+03, 0.92435E+03, 0.11470E+04, 0.14441E+04, 0.18176E+04, 0.22860E+04,&
0.28721E+04, 0.29038E+04, 0.29104E+04, 0.29097E+04, 0.28504E+04, 0.27503E+04,&
0.20956E+04, 0.19717E+04, 0.17926E+04, 0.15245E+04, 0.11243E+04, 0.55892E+03,&
0.21886E+03, 0.00000E+00/
data (ahbar(33,n),n=1,33)/ &
0.28509E+02, 0.35930E+02, 0.45286E+02, 0.57082E+02, 0.71958E+02, 0.90722E+02,&
0.11439E+03, 0.14425E+03, 0.18193E+03, 0.22946E+03, 0.28939E+03, 0.36488E+03,&
0.45978E+03, 0.57859E+03, 0.71675E+03, 0.90062E+03, 0.11306E+04, 0.14176E+04,&
0.17740E+04, 0.22141E+04, 0.22189E+04, 0.21974E+04, 0.21603E+04, 0.20655E+04,&
0.15572E+04, 0.14495E+04, 0.13057E+04, 0.11057E+04, 0.82055E+03, 0.41732E+03,&
0.17442E+03, 0.00000E+00, 0.00000E+00/
data (ahbar(34,n),n=1,34)/ &
0.17955E+02, 0.22627E+02, 0.28517E+02, 0.35941E+02, 0.45302E+02, 0.57105E+02,&
0.71988E+02, 0.90757E+02, 0.11442E+03, 0.14426E+03, 0.18184E+03, 0.22912E+03,&
0.28847E+03, 0.36263E+03, 0.44864E+03, 0.56281E+03, 0.70510E+03, 0.88179E+03,&
0.11001E+04, 0.13676E+04, 0.16919E+04, 0.16753E+04, 0.16315E+04, 0.15655E+04,&
0.11695E+04, 0.10771E+04, 0.95995E+03, 0.80547E+03, 0.59516E+03, 0.30433E+03,&
0.13143E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (ahbar(35,n),n=1,35)/ &
0.11309E+02, 0.14251E+02, 0.17959E+02, 0.22632E+02, 0.28524E+02, 0.35950E+02,&
0.45313E+02, 0.57115E+02, 0.71989E+02, 0.90728E+02, 0.11432E+03, 0.14397E+03,&
0.18114E+03, 0.22752E+03, 0.28119E+03, 0.35228E+03, 0.44063E+03, 0.54992E+03,&
0.68429E+03, 0.84805E+03, 0.10451E+04, 0.12774E+04, 0.12438E+04, 0.11823E+04,&
0.88640E+03, 0.80897E+03, 0.71340E+03, 0.59221E+03, 0.43360E+03, 0.22074E+03,&
0.96191E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (ahbar(36,n),n=1,36)/ &
0.71229E+01, 0.89755E+01, 0.11310E+02, 0.14253E+02, 0.17962E+02, 0.22636E+02,&
0.28527E+02, 0.35951E+02, 0.45304E+02, 0.57082E+02, 0.71899E+02, 0.90509E+02,&
0.11382E+03, 0.14286E+03, 0.17642E+03, 0.22080E+03, 0.27581E+03, 0.34365E+03,&
0.42675E+03, 0.52753E+03, 0.64803E+03, 0.78904E+03, 0.94844E+03, 0.90139E+03,&
0.66943E+03, 0.61315E+03, 0.53582E+03, 0.44015E+03, 0.31885E+03, 0.16089E+03,&
0.69774E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
endif

data pi/3.141592654/

!calculating the mass categories (c.g.s) with the lowest diameter
!of 3.125 microns and mass doubling every bin
diam(1)=1.5625*2.e-04
x(1)=pi/6.*diam(1)**3.
do l=2,ibins+1
   x(l)=2.*x(l-1)
   diam(l)=(6./pi*x(l))**(0.333333333333)
enddo

!Long's or Halls's collection kernel as calculated
!by Walko or Saleeby (1-36) with (1-36) weighted for (x+y)
do i=1,ibins
   do j=1,i
      if(kernel==1) akbar(i,j)=aabar(i,j)
      if(kernel==2) akbar(i,j)=ahbar(i,j)
      if(akbar(i,j).lt.0.) akbar(i,j)=0.
      akbar(36,i)=0.
   enddo
enddo

!Invert here for use in routine 'sxy'
do j=1,ibins
   do i=1,j
      akbar(i,j)=akbar(j,i)
   enddo
enddo

return
end subroutine mcphys_data

!******************************************************************************
subroutine data_cs(x,diam,akbarcs,ibins,cfmass,pwmass)

implicit none
integer :: l,i,j,ibins,n,skernel
real :: pi,cfmass,pwmass
real, dimension(ibins+1) :: x,diam
real, dimension(ibins,ibins) :: akbarcs
real, dimension(36,36) :: abarcs1,abarcs2,abarcs3

skernel=2

if(skernel==1) then  !Snow column kernel based on LENGTH (Wang & Ji 2000)
data (abarcs1( 1,n),n=1, 1)/ &
0.00000E+00/
data (abarcs1( 2,n),n=1, 2)/ &
0.00000E+00, 0.00000E+00/
data (abarcs1( 3,n),n=1, 3)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs1( 4,n),n=1, 4)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs1( 5,n),n=1, 5)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs1( 6,n),n=1, 6)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs1( 7,n),n=1, 7)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcs1( 8,n),n=1, 8)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcs1( 9,n),n=1, 9)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs1(10,n),n=1,10)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs1(11,n),n=1,11)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs1(12,n),n=1,12)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs1(13,n),n=1,13)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcs1(14,n),n=1,14)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcs1(15,n),n=1,15)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs1(16,n),n=1,16)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.31216E+02, 0.56711E+02, 0.88248E+02,&
0.18769E+03, 0.29148E+03, 0.34144E+03, 0.34193E+03, 0.21972E+03, 0.11915E+03,&
0.58562E+03, 0.11056E+04, 0.00000E+00, 0.00000E+00/
data (abarcs1(17,n),n=1,17)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.52550E+02, 0.10357E+03, 0.16736E+03,&
0.27698E+03, 0.39650E+03, 0.47465E+03, 0.51015E+03, 0.43415E+03, 0.24718E+03,&
0.21709E+03, 0.95866E+03, 0.17606E+04, 0.00000E+00, 0.00000E+00/
data (abarcs1(18,n),n=1,18)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.79937E+02, 0.15926E+03, 0.25901E+03,&
0.38224E+03, 0.51951E+03, 0.62226E+03, 0.68448E+03, 0.65184E+03, 0.53382E+03,&
0.25276E+03, 0.42129E+03, 0.16137E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs1(19,n),n=1,19)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.11801E+03, 0.23101E+03, 0.37356E+03,&
0.51499E+03, 0.67301E+03, 0.79486E+03, 0.87227E+03, 0.88162E+03, 0.81734E+03,&
0.62746E+03, 0.25786E+03, 0.78017E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcs1(20,n),n=1,20)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.16466E+03, 0.30611E+03, 0.48484E+03,&
0.64985E+03, 0.83343E+03, 0.97080E+03, 0.10702E+04, 0.11090E+04, 0.10888E+04,&
0.98557E+03, 0.71081E+03, 0.29253E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcs1(21,n),n=1,21)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.22895E+03, 0.39484E+03, 0.60459E+03,&
0.79997E+03, 0.10148E+04, 0.11613E+04, 0.12926E+04, 0.13472E+04, 0.13645E+04,&
0.13358E+04, 0.11816E+04, 0.78564E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs1(22,n),n=1,22)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.29257E+03, 0.48959E+03, 0.73871E+03,&
0.95937E+03, 0.12021E+04, 0.13695E+04, 0.15213E+04, 0.15911E+04, 0.16355E+04,&
0.16545E+04, 0.15832E+04, 0.13670E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs1(23,n),n=1,23)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.37772E+03, 0.61508E+03, 0.91517E+03,&
0.11586E+04, 0.14269E+04, 0.16203E+04, 0.17861E+04, 0.18685E+04, 0.19349E+04,&
0.19864E+04, 0.19684E+04, 0.18756E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs1(24,n),n=1,24)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.44759E+03, 0.72397E+03, 0.10733E+04,&
0.13497E+04, 0.16545E+04, 0.18767E+04, 0.20642E+04, 0.21600E+04, 0.22434E+04,&
0.23186E+04, 0.23352E+04, 0.23202E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs1(25,n),n=1,25)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.51497E+03, 0.83281E+03, 0.12344E+04,&
0.15518E+04, 0.19017E+04, 0.21564E+04, 0.23711E+04, 0.24811E+04, 0.25786E+04,&
0.26719E+04, 0.27101E+04, 0.27417E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcs1(26,n),n=1,26)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.59250E+03, 0.95807E+03, 0.14198E+04,&
0.17844E+04, 0.21859E+04, 0.24775E+04, 0.27227E+04, 0.28474E+04, 0.29582E+04,&
0.30664E+04, 0.31175E+04, 0.31775E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcs1(27,n),n=1,27)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.68175E+03, 0.11022E+04, 0.16331E+04,&
0.20519E+04, 0.25127E+04, 0.28465E+04, 0.31263E+04, 0.32670E+04, 0.33912E+04,&
0.35126E+04, 0.35710E+04, 0.36474E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs1(28,n),n=1,28)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.78447E+03, 0.12681E+04, 0.18786E+04,&
0.23598E+04, 0.28888E+04, 0.32710E+04, 0.35903E+04, 0.37486E+04, 0.38871E+04,&
0.40212E+04, 0.40831E+04, 0.41683E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs1(29,n),n=1,29)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.90270E+03, 0.14591E+04, 0.21612E+04,&
0.27143E+04, 0.33216E+04, 0.37595E+04, 0.41240E+04, 0.43024E+04, 0.44564E+04,&
0.46036E+04, 0.46666E+04, 0.47556E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs1(30,n),n=1,30)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.10388E+04, 0.16790E+04, 0.24865E+04,&
0.31222E+04, 0.38199E+04, 0.43218E+04, 0.47383E+04, 0.49395E+04, 0.51110E+04,&
0.52725E+04, 0.53348E+04, 0.54243E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs1(31,n),n=1,31)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.11955E+04, 0.19320E+04, 0.28609E+04,&
0.35918E+04, 0.43934E+04, 0.49691E+04, 0.54455E+04, 0.56730E+04, 0.58644E+04,&
0.60418E+04, 0.61024E+04, 0.61903E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcs1(32,n),n=1,32)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.13758E+04, 0.22233E+04, 0.32919E+04,&
0.41324E+04, 0.50537E+04, 0.57144E+04, 0.62596E+04, 0.65174E+04, 0.67318E+04,&
0.69274E+04, 0.69856E+04, 0.70705E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcs1(33,n),n=1,33)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.15833E+04, 0.25586E+04, 0.37881E+04,&
0.47548E+04, 0.58139E+04, 0.65723E+04, 0.71971E+04, 0.74898E+04, 0.77308E+04,&
0.79474E+04, 0.80027E+04, 0.80838E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs1(34,n),n=1,34)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.18222E+04, 0.29445E+04, 0.43592E+04,&
0.54711E+04, 0.66889E+04, 0.75601E+04, 0.82764E+04, 0.86095E+04, 0.88812E+04,&
0.91222E+04, 0.91745E+04, 0.92513E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs1(35,n),n=1,35)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.20972E+04, 0.33887E+04, 0.50166E+04,&
0.62958E+04, 0.76963E+04, 0.86973E+04, 0.95191E+04, 0.98988E+04, 0.10206E+05,&
0.10476E+05, 0.10525E+05, 0.10597E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs1(36,n),n=1,36)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.24137E+04, 0.39001E+04, 0.57734E+04,&
0.72450E+04, 0.88559E+04, 0.10006E+05, 0.10950E+05, 0.11383E+05, 0.11732E+05,&
0.12035E+05, 0.12081E+05, 0.12148E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
endif

if(skernel==2) then  !Snow column kernel based on WIDTH (Wang & Ji 2000)
data (abarcs2( 1,n),n=1, 1)/ &
0.00000E+00/
data (abarcs2( 2,n),n=1, 2)/ &
0.00000E+00, 0.00000E+00/
data (abarcs2( 3,n),n=1, 3)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs2( 4,n),n=1, 4)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs2( 5,n),n=1, 5)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs2( 6,n),n=1, 6)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs2( 7,n),n=1, 7)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcs2( 8,n),n=1, 8)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcs2( 9,n),n=1, 9)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs2(10,n),n=1,10)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs2(11,n),n=1,11)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs2(12,n),n=1,12)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs2(13,n),n=1,13)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcs2(14,n),n=1,14)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcs2(15,n),n=1,15)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs2(16,n),n=1,16)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.72704E+02, 0.14215E+03, 0.22811E+03,&
0.31388E+03, 0.40083E+03, 0.44514E+03, 0.42268E+03, 0.27852E+03, 0.15034E+03,&
0.80020E+03, 0.22832E+04, 0.49348E+04, 0.00000E+00/
data (abarcs2(17,n),n=1,17)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.12878E+03, 0.22237E+03, 0.33933E+03,&
0.44479E+03, 0.55366E+03, 0.61013E+03, 0.62603E+03, 0.54080E+03, 0.30859E+03,&
0.28557E+03, 0.14261E+04, 0.39526E+04, 0.00000E+00, 0.00000E+00/
data (abarcs2(18,n),n=1,18)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.19263E+03, 0.31128E+03, 0.46044E+03,&
0.57617E+03, 0.69915E+03, 0.77670E+03, 0.81805E+03, 0.77874E+03, 0.64207E+03,&
0.31390E+03, 0.53669E+03, 0.23503E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs2(19,n),n=1,19)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.22187E+03, 0.35881E+03, 0.53154E+03,&
0.66704E+03, 0.81385E+03, 0.91380E+03, 0.98361E+03, 0.98249E+03, 0.91835E+03,&
0.72268E+03, 0.30165E+03, 0.95257E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcs2(20,n),n=1,20)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.25540E+03, 0.41318E+03, 0.61253E+03,&
0.76982E+03, 0.94195E+03, 0.10637E+04, 0.11582E+04, 0.11851E+04, 0.11712E+04,&
0.10765E+04, 0.78283E+03, 0.32584E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcs2(21,n),n=1,21)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.29390E+03, 0.47552E+03, 0.70514E+03,&
0.88679E+03, 0.10866E+04, 0.12306E+04, 0.13479E+04, 0.13965E+04, 0.14180E+04,&
0.13911E+04, 0.12320E+04, 0.82728E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs2(22,n),n=1,22)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.33815E+03, 0.54709E+03, 0.81131E+03,&
0.10205E+04, 0.12511E+04, 0.14187E+04, 0.15584E+04, 0.16246E+04, 0.16720E+04,&
0.16914E+04, 0.16185E+04, 0.14065E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs2(23,n),n=1,23)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.38904E+03, 0.62936E+03, 0.93320E+03,&
0.11737E+04, 0.14391E+04, 0.16325E+04, 0.17953E+04, 0.18769E+04, 0.19440E+04,&
0.19958E+04, 0.19777E+04, 0.18870E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs2(24,n),n=1,24)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.44759E+03, 0.72397E+03, 0.10733E+04,&
0.13497E+04, 0.16545E+04, 0.18767E+04, 0.20642E+04, 0.21600E+04, 0.22434E+04,&
0.23186E+04, 0.23352E+04, 0.23202E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs2(25,n),n=1,25)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.51497E+03, 0.83281E+03, 0.12344E+04,&
0.15518E+04, 0.19017E+04, 0.21564E+04, 0.23711E+04, 0.24811E+04, 0.25786E+04,&
0.26719E+04, 0.27101E+04, 0.27417E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcs2(26,n),n=1,26)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.59250E+03, 0.95807E+03, 0.14198E+04,&
0.17844E+04, 0.21859E+04, 0.24775E+04, 0.27227E+04, 0.28474E+04, 0.29582E+04,&
0.30664E+04, 0.31175E+04, 0.31775E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcs2(27,n),n=1,27)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.68175E+03, 0.11022E+04, 0.16331E+04,&
0.20519E+04, 0.25127E+04, 0.28465E+04, 0.31263E+04, 0.32670E+04, 0.33912E+04,&
0.35126E+04, 0.35710E+04, 0.36474E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs2(28,n),n=1,28)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.78447E+03, 0.12681E+04, 0.18786E+04,&
0.23598E+04, 0.28888E+04, 0.32710E+04, 0.35903E+04, 0.37486E+04, 0.38871E+04,&
0.40212E+04, 0.40831E+04, 0.41683E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs2(29,n),n=1,29)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.90270E+03, 0.14591E+04, 0.21612E+04,&
0.27143E+04, 0.33216E+04, 0.37595E+04, 0.41240E+04, 0.43024E+04, 0.44564E+04,&
0.46036E+04, 0.46666E+04, 0.47556E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs2(30,n),n=1,30)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.10388E+04, 0.16790E+04, 0.24865E+04,&
0.31222E+04, 0.38199E+04, 0.43218E+04, 0.47383E+04, 0.49395E+04, 0.51110E+04,&
0.52725E+04, 0.53348E+04, 0.54243E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs2(31,n),n=1,31)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.11955E+04, 0.19320E+04, 0.28609E+04,&
0.35918E+04, 0.43934E+04, 0.49691E+04, 0.54455E+04, 0.56730E+04, 0.58644E+04,&
0.60418E+04, 0.61024E+04, 0.61903E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcs2(32,n),n=1,32)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.13758E+04, 0.22233E+04, 0.32919E+04,&
0.41324E+04, 0.50537E+04, 0.57144E+04, 0.62596E+04, 0.65174E+04, 0.67318E+04,&
0.69274E+04, 0.69856E+04, 0.70705E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcs2(33,n),n=1,33)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.15833E+04, 0.25586E+04, 0.37881E+04,&
0.47548E+04, 0.58139E+04, 0.65723E+04, 0.71971E+04, 0.74898E+04, 0.77308E+04,&
0.79474E+04, 0.80027E+04, 0.80838E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs2(34,n),n=1,34)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.18222E+04, 0.29445E+04, 0.43592E+04,&
0.54711E+04, 0.66889E+04, 0.75601E+04, 0.82764E+04, 0.86095E+04, 0.88812E+04,&
0.91222E+04, 0.91745E+04, 0.92513E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs2(35,n),n=1,35)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.20972E+04, 0.33887E+04, 0.50166E+04,&
0.62958E+04, 0.76963E+04, 0.86973E+04, 0.95191E+04, 0.98988E+04, 0.10206E+05,&
0.10476E+05, 0.10525E+05, 0.10597E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs2(36,n),n=1,36)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.24137E+04, 0.39001E+04, 0.57734E+04,&
0.72450E+04, 0.88559E+04, 0.10006E+05, 0.10950E+05, 0.11383E+05, 0.11732E+05,&
0.12035E+05, 0.12081E+05, 0.12148E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
endif

if(skernel==3) then  !Snow broad banched crystals (Wang & Ji 2000)
data (abarcs3( 1,n),n=1, 1)/ &
0.00000E+00/
data (abarcs3( 2,n),n=1, 2)/ &
0.00000E+00, 0.00000E+00/
data (abarcs3( 3,n),n=1, 3)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3( 4,n),n=1, 4)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3( 5,n),n=1, 5)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3( 6,n),n=1, 6)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3( 7,n),n=1, 7)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcs3( 8,n),n=1, 8)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcs3( 9,n),n=1, 9)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3(10,n),n=1,10)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3(11,n),n=1,11)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3(12,n),n=1,12)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3(13,n),n=1,13)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcs3(14,n),n=1,14)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcs3(15,n),n=1,15)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3(16,n),n=1,16)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.11795E+02, 0.95841E+02, 0.52095E+03, 0.10372E+04, 0.16533E+04, 0.19705E+04,&
0.35602E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3(17,n),n=1,17)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.41069E+02, 0.44906E+02, 0.50035E+02,&
0.11143E+03, 0.37952E+03, 0.15710E+04, 0.29705E+04, 0.46371E+04, 0.55880E+04,&
0.20941E+04, 0.12015E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3(18,n),n=1,18)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.13011E+03, 0.14183E+03, 0.15733E+03,&
0.27620E+03, 0.60142E+03, 0.17631E+04, 0.30710E+04, 0.46529E+04, 0.57245E+04,&
0.41760E+04, 0.35795E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3(19,n),n=1,19)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.26149E+03, 0.28431E+03, 0.31426E+03,&
0.51578E+03, 0.92432E+03, 0.20654E+04, 0.32879E+04, 0.48199E+04, 0.60845E+04,&
0.68823E+04, 0.66054E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcs3(20,n),n=1,20)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.36547E+03, 0.43104E+03, 0.51590E+03,&
0.91906E+03, 0.15345E+04, 0.26755E+04, 0.37939E+04, 0.51560E+04, 0.62769E+04,&
0.75561E+04, 0.87690E+04, 0.89202E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcs3(21,n),n=1,21)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.44249E+03, 0.52787E+03, 0.63774E+03,&
0.97834E+03, 0.15700E+04, 0.29526E+04, 0.39091E+04, 0.50700E+04, 0.60515E+04,&
0.68943E+04, 0.78737E+04, 0.88749E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3(22,n),n=1,22)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.61561E+03, 0.75660E+03, 0.93717E+03,&
0.14355E+04, 0.21229E+04, 0.31551E+04, 0.40377E+04, 0.50532E+04, 0.59179E+04,&
0.68270E+04, 0.79620E+04, 0.91102E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3(23,n),n=1,23)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.70516E+03, 0.79166E+03, 0.90233E+03,&
0.13386E+04, 0.19532E+04, 0.29246E+04, 0.37268E+04, 0.46359E+04, 0.55447E+04,&
0.62921E+04, 0.75094E+04, 0.85887E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3(24,n),n=1,24)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.70933E+03, 0.79568E+03, 0.90583E+03,&
0.13415E+04, 0.19525E+04, 0.29132E+04, 0.36938E+04, 0.45631E+04, 0.54062E+04,&
0.60593E+04, 0.71210E+04, 0.80045E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3(25,n),n=1,25)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.71379E+03, 0.80016E+03, 0.91006E+03,&
0.13459E+04, 0.19550E+04, 0.29087E+04, 0.36731E+04, 0.45117E+04, 0.53037E+04,&
0.58820E+04, 0.68196E+04, 0.75422E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcs3(26,n),n=1,26)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.71849E+03, 0.80501E+03, 0.91487E+03,&
0.13515E+04, 0.19601E+04, 0.29095E+04, 0.36622E+04, 0.44776E+04, 0.52297E+04,&
0.57492E+04, 0.65886E+04, 0.71811E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcs3(27,n),n=1,27)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.72339E+03, 0.81016E+03, 0.92017E+03,&
0.13581E+04, 0.19672E+04, 0.29147E+04, 0.36592E+04, 0.44572E+04, 0.51785E+04,&
0.56517E+04, 0.64137E+04, 0.69023E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3(28,n),n=1,28)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.72845E+03, 0.81556E+03, 0.92586E+03,&
0.13656E+04, 0.19760E+04, 0.29235E+04, 0.36624E+04, 0.44477E+04, 0.51455E+04,&
0.55824E+04, 0.62835E+04, 0.66895E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3(29,n),n=1,29)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.73367E+03, 0.82118E+03, 0.93188E+03,&
0.13737E+04, 0.19861E+04, 0.29350E+04, 0.36707E+04, 0.44470E+04, 0.51271E+04,&
0.55354E+04, 0.61890E+04, 0.65293E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3(30,n),n=1,30)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.73900E+03, 0.82698E+03, 0.93818E+03,&
0.13823E+04, 0.19974E+04, 0.29489E+04, 0.36831E+04, 0.44534E+04, 0.51202E+04,&
0.55064E+04, 0.61228E+04, 0.64110E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3(31,n),n=1,31)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.74445E+03, 0.83294E+03, 0.94471E+03,&
0.13915E+04, 0.20095E+04, 0.29647E+04, 0.36987E+04, 0.44654E+04, 0.51226E+04,&
0.54916E+04, 0.60792E+04, 0.63261E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcs3(32,n),n=1,32)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.74999E+03, 0.83902E+03, 0.95143E+03,&
0.14010E+04, 0.20224E+04, 0.29819E+04, 0.37171E+04, 0.44820E+04, 0.51326E+04,&
0.54882E+04, 0.60536E+04, 0.62678E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcs3(33,n),n=1,33)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.75561E+03, 0.84523E+03, 0.95832E+03,&
0.14108E+04, 0.20360E+04, 0.30005E+04, 0.37377E+04, 0.45023E+04, 0.51485E+04,&
0.54941E+04, 0.60425E+04, 0.62307E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3(34,n),n=1,34)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.76132E+03, 0.85154E+03, 0.96535E+03,&
0.14209E+04, 0.20500E+04, 0.30201E+04, 0.37600E+04, 0.45257E+04, 0.51693E+04,&
0.55073E+04, 0.60429E+04, 0.62105E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3(35,n),n=1,35)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.76709E+03, 0.85794E+03, 0.97252E+03,&
0.14313E+04, 0.20645E+04, 0.30405E+04, 0.37839E+04, 0.45515E+04, 0.51941E+04,&
0.55265E+04, 0.60527E+04, 0.62040E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs3(36,n),n=1,36)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.77294E+03, 0.86443E+03, 0.97980E+03,&
0.14418E+04, 0.20794E+04, 0.30617E+04, 0.38089E+04, 0.45794E+04, 0.52222E+04,&
0.55505E+04, 0.60699E+04, 0.62083E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
endif

data pi/3.141592654/

!calculating the mass categories (c.g.s) with the lowest diameter
!of 3.125 microns and mass doubling every bin
diam(1)=1.5625*2.e-04
x(1)=pi/6.*diam(1)**3.
do l=2,15
   x(l)=2.*x(l-1)
   diam(l)=(6./pi*x(l))**(0.333333333333)
enddo
do l=16,ibins+1
   x(l)=2.*x(l-1)
   !uses power law values in m.k.s. and converts back to c.g.s.
   diam(l)=100 * ((x(l) / 1000.0 / cfmass) ** (1./pwmass))
enddo

!(Wang & Ji 2000) collection kernel by Saleeby 6/23/04 (1-36) with (1-36)
! weighted for (x+y)
do i=1,ibins
   do j=1,i
      if(skernel==1) akbarcs(i,j)=abarcs1(i,j)
      if(skernel==2) akbarcs(i,j)=abarcs2(i,j)
      if(skernel==3) akbarcs(i,j)=abarcs3(i,j)
      if(akbarcs(i,j).lt.0.) akbarcs(i,j)=0.
      akbarcs(36,i)=0.
   enddo
enddo

!Invert here for use in routine 'sxy'
do j=1,ibins
   do i=1,j
      akbarcs(i,j)=akbarcs(j,i)
   enddo
enddo
return
end subroutine data_cs

!******************************************************************************
subroutine data_ca(x,diam,akbarca,ibins,cfmass,pwmass)

implicit none
integer :: l,i,j,ibins,n
real :: pi,cfmass,pwmass
real, dimension(ibins+1) :: x,diam
real, dimension(ibins,ibins) :: akbarca
real, dimension(36,36) :: abarca

data (abarca( 1,n),n=1, 1)/ &
0.00000E+00/
data (abarca( 2,n),n=1, 2)/ &
0.00000E+00, 0.00000E+00/
data (abarca( 3,n),n=1, 3)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca( 4,n),n=1, 4)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca( 5,n),n=1, 5)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca( 6,n),n=1, 6)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca( 7,n),n=1, 7)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarca( 8,n),n=1, 8)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarca( 9,n),n=1, 9)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(10,n),n=1,10)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(11,n),n=1,11)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(12,n),n=1,12)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(13,n),n=1,13)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarca(14,n),n=1,14)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarca(15,n),n=1,15)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(16,n),n=1,16)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.33459E+03, 0.23632E+04, 0.53117E+04, 0.62743E+04, 0.31416E+04,&
0.22027E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(17,n),n=1,17)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.16956E+02, 0.24102E+02, 0.33700E+02,&
0.10696E+03, 0.60496E+03, 0.30958E+04, 0.48229E+04, 0.59387E+04, 0.67221E+04,&
0.53597E+04, 0.67867E+03, 0.58085E+03, 0.00000E+00, 0.00000E+00/
data (abarca(18,n),n=1,18)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.84307E+02, 0.11916E+03, 0.16539E+03,&
0.52011E+03, 0.11898E+04, 0.28912E+04, 0.41915E+04, 0.50265E+04, 0.57334E+04,&
0.53183E+04, 0.30534E+04, 0.27060E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(19,n),n=1,19)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.14764E+03, 0.20771E+03, 0.28658E+03,&
0.89434E+03, 0.17192E+04, 0.27848E+04, 0.37669E+04, 0.44009E+04, 0.50567E+04,&
0.53169E+04, 0.47997E+04, 0.43101E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarca(20,n),n=1,20)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.20987E+03, 0.31676E+03, 0.45590E+03,&
0.11075E+04, 0.19046E+04, 0.25945E+04, 0.32667E+04, 0.38062E+04, 0.43542E+04,&
0.47991E+04, 0.50310E+04, 0.46364E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarca(21,n),n=1,21)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.26446E+03, 0.43289E+03, 0.65060E+03,&
0.11449E+04, 0.17435E+04, 0.22820E+04, 0.26740E+04, 0.31657E+04, 0.35495E+04,&
0.38183E+04, 0.40415E+04, 0.39004E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(22,n),n=1,22)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.28620E+03, 0.47498E+03, 0.71763E+03,&
0.11065E+04, 0.15672E+04, 0.19636E+04, 0.22605E+04, 0.25943E+04, 0.28461E+04,&
0.30174E+04, 0.31500E+04, 0.30876E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(23,n),n=1,23)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.27979E+03, 0.45329E+03, 0.67529E+03,&
0.97222E+03, 0.13240E+04, 0.16432E+04, 0.19107E+04, 0.21302E+04, 0.23112E+04,&
0.24208E+04, 0.24960E+04, 0.24245E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(24,n),n=1,24)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.26258E+03, 0.40737E+03, 0.59197E+03,&
0.83777E+03, 0.11270E+04, 0.13862E+04, 0.15987E+04, 0.17468E+04, 0.18788E+04,&
0.19486E+04, 0.19883E+04, 0.19155E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(25,n),n=1,25)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.21999E+03, 0.34088E+03, 0.49458E+03,&
0.69850E+03, 0.93713E+03, 0.11487E+04, 0.13190E+04, 0.14333E+04, 0.15312E+04,&
0.15756E+04, 0.15942E+04, 0.15253E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarca(26,n),n=1,26)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.18445E+03, 0.28554E+03, 0.41377E+03,&
0.58343E+03, 0.78109E+03, 0.95481E+03, 0.10925E+04, 0.11821E+04, 0.12561E+04,&
0.12845E+04, 0.12910E+04, 0.12287E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarca(27,n),n=1,27)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.15474E+03, 0.23938E+03, 0.34654E+03,&
0.48801E+03, 0.65226E+03, 0.79562E+03, 0.90792E+03, 0.97895E+03, 0.10360E+04,&
0.10542E+04, 0.10543E+04, 0.99952E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(28,n),n=1,28)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.12989E+03, 0.20081E+03, 0.29048E+03,&
0.40866E+03, 0.54549E+03, 0.66428E+03, 0.75644E+03, 0.81348E+03, 0.85811E+03,&
0.87007E+03, 0.86686E+03, 0.81977E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(29,n),n=1,29)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.10907E+03, 0.16854E+03, 0.24366E+03,&
0.34251E+03, 0.45674E+03, 0.55549E+03, 0.63153E+03, 0.67779E+03, 0.71327E+03,&
0.72127E+03, 0.71674E+03, 0.67687E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(30,n),n=1,30)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.91610E+02, 0.14152E+03, 0.20449E+03,&
0.28728E+03, 0.38279E+03, 0.46510E+03, 0.52812E+03, 0.56595E+03, 0.59453E+03,&
0.60006E+03, 0.59528E+03, 0.56191E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(31,n),n=1,31)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.76966E+02, 0.11886E+03, 0.17169E+03,&
0.24110E+03, 0.32106E+03, 0.38980E+03, 0.44221E+03, 0.47337E+03, 0.49665E+03,&
0.50063E+03, 0.49617E+03, 0.46850E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarca(32,n),n=1,32)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.64676E+02, 0.99860E+02, 0.14421E+03,&
0.20243E+03, 0.26944E+03, 0.32695E+03, 0.37066E+03, 0.39646E+03, 0.41561E+03,&
0.41861E+03, 0.41473E+03, 0.39195E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarca(33,n),n=1,33)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.54356E+02, 0.83913E+02, 0.12115E+03,&
0.17002E+03, 0.22623E+03, 0.27440E+03, 0.31093E+03, 0.33240E+03, 0.34827E+03,&
0.35065E+03, 0.34742E+03, 0.32878E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(34,n),n=1,34)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.45689E+02, 0.70524E+02, 0.10180E+03,&
0.14284E+03, 0.19002E+03, 0.23041E+03, 0.26100E+03, 0.27893E+03, 0.29215E+03,&
0.29412E+03, 0.29154E+03, 0.27636E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(35,n),n=1,35)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.38407E+02, 0.59278E+02, 0.85560E+02,&
0.12003E+03, 0.15965E+03, 0.19355E+03, 0.21920E+03, 0.23421E+03, 0.24529E+03,&
0.24697E+03, 0.24497E+03, 0.23265E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(36,n),n=1,36)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.32288E+02, 0.49831E+02, 0.71918E+02,&
0.10088E+03, 0.13416E+03, 0.16263E+03, 0.18416E+03, 0.19675E+03, 0.20607E+03,&
0.20755E+03, 0.20604E+03, 0.19608E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/

data pi/3.141592654/

!calculating the mass categories (c.g.s) with the lowest diameter
!of 3.125 microns and mass doubling every bin
diam(1)=1.5625*2.e-04
x(1)=pi/6.*diam(1)**3.
do l=2,15
   x(l)=2.*x(l-1)
   diam(l)=(6./pi*x(l))**(0.333333333333)
enddo
do l=16,ibins+1
   x(l)=2.*x(l-1)
   !uses power law values in m.k.s. and converts back to c.g.s.
   diam(l)=100 * ((x(l) / 1000.0 / cfmass) ** (1./pwmass))
enddo

!(Cober & List) collection kernel by Saleeby 6/23/04 for graupel/cloud droplet
! weighted for (x+y)
do i=1,ibins
   do j=1,i
      akbarca(i,j)=abarca(i,j)
      if(akbarca(i,j).lt.0.) akbarca(i,j)=0.
      akbarca(36,i)=0.
   enddo
enddo

!Invert here for use in routine 'sxy'
do j=1,ibins
   do i=1,j
      akbarca(i,j)=akbarca(j,i)
   enddo
enddo
return
end subroutine data_ca

!******************************************************************************
subroutine data_cg(x,diam,akbarcg,ibins,cfmass,pwmass)

implicit none
integer :: l,i,j,ibins,n
real :: pi,cfmass,pwmass
real, dimension(ibins+1) :: x,diam
real, dimension(ibins,ibins) :: akbarcg
real, dimension(36,36) :: abarcg

data (abarcg( 1,n),n=1, 1)/ &
0.00000E+00/
data (abarcg( 2,n),n=1, 2)/ &
0.00000E+00, 0.00000E+00/
data (abarcg( 3,n),n=1, 3)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg( 4,n),n=1, 4)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg( 5,n),n=1, 5)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg( 6,n),n=1, 6)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg( 7,n),n=1, 7)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcg( 8,n),n=1, 8)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcg( 9,n),n=1, 9)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(10,n),n=1,10)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(11,n),n=1,11)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(12,n),n=1,12)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(13,n),n=1,13)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcg(14,n),n=1,14)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcg(15,n),n=1,15)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(16,n),n=1,16)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.67855E+04, 0.92403E+04, 0.12604E+05,&
0.15715E+05, 0.19798E+05, 0.24668E+05, 0.29470E+05, 0.32895E+05, 0.36737E+05,&
0.40385E+05, 0.42478E+05, 0.40953E+05, 0.00000E+00/
data (abarcg(17,n),n=1,17)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.53655E+04, 0.73931E+04, 0.10134E+05,&
0.12602E+05, 0.15778E+05, 0.19424E+05, 0.23432E+05, 0.26028E+05, 0.28902E+05,&
0.32000E+05, 0.34669E+05, 0.35548E+05, 0.00000E+00, 0.00000E+00/
data (abarcg(18,n),n=1,18)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.42064E+04, 0.59049E+04, 0.81733E+04,&
0.10201E+05, 0.12744E+05, 0.15463E+05, 0.18807E+05, 0.20863E+05, 0.22905E+05,&
0.25241E+05, 0.27593E+05, 0.29240E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(19,n),n=1,19)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.33647E+04, 0.48074E+04, 0.67139E+04,&
0.83270E+04, 0.10339E+05, 0.12553E+05, 0.15029E+05, 0.16969E+05, 0.18385E+05,&
0.20050E+05, 0.21849E+05, 0.23431E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcg(20,n),n=1,20)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.26252E+04, 0.38617E+04, 0.54812E+04,&
0.68451E+04, 0.85109E+04, 0.10238E+05, 0.12227E+05, 0.14000E+05, 0.14972E+05,&
0.16126E+05, 0.17414E+05, 0.18663E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcg(21,n),n=1,21)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.20424E+04, 0.31113E+04, 0.45008E+04,&
0.56170E+04, 0.69781E+04, 0.84371E+04, 0.10059E+05, 0.11550E+05, 0.12366E+05,&
0.13155E+05, 0.14048E+05, 0.14951E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(22,n),n=1,22)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.15389E+04, 0.24686E+04, 0.36700E+04,&
0.46367E+04, 0.58016E+04, 0.70155E+04, 0.82899E+04, 0.95844E+04, 0.10347E+05,&
0.10883E+05, 0.11490E+05, 0.12115E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(23,n),n=1,23)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.11747E+04, 0.19882E+04, 0.30342E+04,&
0.38472E+04, 0.48214E+04, 0.58403E+04, 0.68960E+04, 0.79984E+04, 0.87531E+04,&
0.91183E+04, 0.95271E+04, 0.99473E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(24,n),n=1,24)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.87819E+03, 0.15930E+04, 0.25081E+04,&
0.32121E+04, 0.40498E+04, 0.49158E+04, 0.57581E+04, 0.67068E+04, 0.74672E+04,&
0.77254E+04, 0.79984E+04, 0.82757E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(25,n),n=1,25)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.65329E+03, 0.12834E+04, 0.20875E+04,&
0.27010E+04, 0.34268E+04, 0.41701E+04, 0.48833E+04, 0.56770E+04, 0.64231E+04,&
0.66079E+04, 0.67886E+04, 0.69684E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcg(26,n),n=1,26)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.57923E+03, 0.11365E+04, 0.18458E+04,&
0.23838E+04, 0.30174E+04, 0.36614E+04, 0.42729E+04, 0.49469E+04, 0.55694E+04,&
0.56970E+04, 0.58157E+04, 0.59298E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcg(27,n),n=1,27)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.51407E+03, 0.10077E+04, 0.16346E+04,&
0.21079E+04, 0.26632E+04, 0.32244E+04, 0.37525E+04, 0.43300E+04, 0.48556E+04,&
0.49441E+04, 0.50209E+04, 0.50915E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(28,n),n=1,28)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.45660E+03, 0.89433E+03, 0.14493E+04,&
0.18668E+04, 0.23552E+04, 0.28463E+04, 0.33052E+04, 0.38037E+04, 0.42521E+04,&
0.43137E+04, 0.43627E+04, 0.44045E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(29,n),n=1,29)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.40581E+03, 0.79436E+03, 0.12863E+04,&
0.16553E+04, 0.20859E+04, 0.25173E+04, 0.29181E+04, 0.33512E+04, 0.37369E+04,&
0.37800E+04, 0.38105E+04, 0.38338E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(30,n),n=1,30)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.36084E+03, 0.70600E+03, 0.11426E+04,&
0.14692E+04, 0.18497E+04, 0.22298E+04, 0.25812E+04, 0.29594E+04, 0.32935E+04,&
0.33240E+04, 0.33423E+04, 0.33538E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(31,n),n=1,31)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.32099E+03, 0.62779E+03, 0.10155E+04,&
0.13050E+04, 0.16419E+04, 0.19774E+04, 0.22866E+04, 0.26183E+04, 0.29094E+04,&
0.29311E+04, 0.29415E+04, 0.29457E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcg(32,n),n=1,32)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.28563E+03, 0.55846E+03, 0.90303E+03,&
0.11599E+04, 0.14585E+04, 0.17553E+04, 0.20281E+04, 0.23199E+04, 0.25748E+04,&
0.25905E+04, 0.25958E+04, 0.25956E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcg(33,n),n=1,33)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.25422E+03, 0.49694E+03, 0.80331E+03,&
0.10315E+04, 0.12964E+04, 0.15594E+04, 0.18005E+04, 0.20579E+04, 0.22820E+04,&
0.22935E+04, 0.22957E+04, 0.22930E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(34,n),n=1,34)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.22632E+03, 0.44231E+03, 0.71483E+03,&
0.91761E+03, 0.11528E+04, 0.13861E+04, 0.15996E+04, 0.18272E+04, 0.20247E+04,&
0.20334E+04, 0.20337E+04, 0.20298E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(35,n),n=1,35)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.20150E+03, 0.39376E+03, 0.63625E+03,&
0.81655E+03, 0.10256E+04, 0.12327E+04, 0.14220E+04, 0.16236E+04, 0.17981E+04,&
0.18048E+04, 0.18040E+04, 0.17997E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(36,n),n=1,36)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.17944E+03, 0.35059E+03, 0.56642E+03,&
0.72680E+03, 0.91267E+03, 0.10967E+04, 0.12647E+04, 0.14435E+04, 0.15980E+04,&
0.16033E+04, 0.16019E+04, 0.15976E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/

data pi/3.141592654/

!calculating the mass categories (c.g.s) with the lowest diameter
!of 3.125 microns and mass doubling every bin
diam(1)=1.5625*2.e-04
x(1)=pi/6.*diam(1)**3.
do l=2,15
   x(l)=2.*x(l-1)
   diam(l)=(6./pi*x(l))**(0.333333333333)
enddo
do l=16,ibins+1
   x(l)=2.*x(l-1)
   !uses power law values in m.k.s. and converts back to c.g.s.
   diam(l)=100 * ((x(l) / 1000.0 / cfmass) ** (1./pwmass))
enddo

!(Cober & List) collection kernel by Saleeby 6/23/04 for graupel/cloud droplet
! weighted for (x+y)
do i=1,ibins
   do j=1,i
      akbarcg(i,j)=abarcg(i,j)
      if(akbarcg(i,j).lt.0.) akbarcg(i,j)=0.
      akbarcg(36,i)=0.
   enddo
enddo

!Invert here for use in routine 'sxy'
do j=1,ibins
   do i=1,j
      akbarcg(i,j)=akbarcg(j,i)
   enddo
enddo
return
end subroutine data_cg

!******************************************************************************
subroutine data_ch(x,diam,akbarch,ibins,cfmass,pwmass)

implicit none
integer :: l,i,j,ibins,n
real :: pi,cfmass,pwmass
real, dimension(ibins+1) :: x,diam
real, dimension(ibins,ibins) :: akbarch
real, dimension(36,36) :: abarch

data (abarch( 1,n),n=1, 1)/ &
0.00000E+00/
data (abarch( 2,n),n=1, 2)/ &
0.00000E+00, 0.00000E+00/
data (abarch( 3,n),n=1, 3)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch( 4,n),n=1, 4)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch( 5,n),n=1, 5)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch( 6,n),n=1, 6)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch( 7,n),n=1, 7)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarch( 8,n),n=1, 8)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarch( 9,n),n=1, 9)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(10,n),n=1,10)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(11,n),n=1,11)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(12,n),n=1,12)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(13,n),n=1,13)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarch(14,n),n=1,14)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarch(15,n),n=1,15)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(16,n),n=1,16)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.18811E+05, 0.20803E+05, 0.23441E+05,&
0.26241E+05, 0.29685E+05, 0.33074E+05, 0.35718E+05, 0.38736E+05, 0.41990E+05,&
0.44862E+05, 0.46104E+05, 0.43975E+05, 0.00000E+00/
data (abarch(17,n),n=1,17)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.16123E+05, 0.17698E+05, 0.19768E+05,&
0.21963E+05, 0.24660E+05, 0.27382E+05, 0.29424E+05, 0.31632E+05, 0.34181E+05,&
0.36838E+05, 0.38983E+05, 0.39415E+05, 0.00000E+00, 0.00000E+00/
data (abarch(18,n),n=1,18)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.13837E+05, 0.15087E+05, 0.16717E+05,&
0.18462E+05, 0.20616E+05, 0.22884E+05, 0.24509E+05, 0.26075E+05, 0.27943E+05,&
0.30040E+05, 0.32103E+05, 0.33491E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(19,n),n=1,19)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.11886E+05, 0.12882E+05, 0.14175E+05,&
0.15582E+05, 0.17331E+05, 0.19271E+05, 0.20630E+05, 0.21732E+05, 0.23056E+05,&
0.24596E+05, 0.26246E+05, 0.27700E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarch(20,n),n=1,20)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.10286E+05, 0.11124E+05, 0.12204E+05,&
0.13341E+05, 0.14745E+05, 0.16296E+05, 0.17497E+05, 0.18309E+05, 0.19232E+05,&
0.20322E+05, 0.21537E+05, 0.22732E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarch(21,n),n=1,21)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.88853E+04, 0.95962E+04, 0.10510E+05,&
0.11433E+05, 0.12567E+05, 0.13827E+05, 0.14956E+05, 0.15578E+05, 0.16215E+05,&
0.16972E+05, 0.17830E+05, 0.18716E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(22,n),n=1,22)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.76635E+04, 0.82723E+04, 0.90515E+04,&
0.98024E+04, 0.10723E+05, 0.11761E+05, 0.12846E+05, 0.13366E+05, 0.13805E+05,&
0.14325E+05, 0.14917E+05, 0.15539E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(23,n),n=1,23)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.66687E+04, 0.71944E+04, 0.78652E+04,&
0.84747E+04, 0.92236E+04, 0.10090E+05, 0.10982E+05, 0.11552E+05, 0.11853E+05,&
0.12207E+05, 0.12610E+05, 0.13034E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(24,n),n=1,24)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.58010E+04, 0.62578E+04, 0.68392E+04,&
0.73428E+04, 0.79596E+04, 0.86767E+04, 0.94191E+04, 0.10029E+05, 0.10249E+05,&
0.10490E+05, 0.10760E+05, 0.11043E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(25,n),n=1,25)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.50520E+04, 0.54509E+04, 0.59576E+04,&
0.63920E+04, 0.69147E+04, 0.74869E+04, 0.81114E+04, 0.87268E+04, 0.89141E+04,&
0.90768E+04, 0.92573E+04, 0.94416E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarch(26,n),n=1,26)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.44066E+04, 0.47450E+04, 0.51740E+04,&
0.55505E+04, 0.59998E+04, 0.64767E+04, 0.70082E+04, 0.75638E+04, 0.77908E+04,&
0.79001E+04, 0.80190E+04, 0.81368E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarch(27,n),n=1,27)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.38541E+04, 0.41407E+04, 0.45036E+04,&
0.48316E+04, 0.52196E+04, 0.56163E+04, 0.60726E+04, 0.65482E+04, 0.68345E+04,&
0.69086E+04, 0.69859E+04, 0.70592E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(28,n),n=1,28)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.33767E+04, 0.36294E+04, 0.39491E+04,&
0.42324E+04, 0.45626E+04, 0.48801E+04, 0.52742E+04, 0.56933E+04, 0.60119E+04,&
0.60648E+04, 0.61143E+04, 0.61581E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(29,n),n=1,29)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.29476E+04, 0.31680E+04, 0.34464E+04,&
0.36778E+04, 0.39564E+04, 0.42702E+04, 0.46056E+04, 0.49581E+04, 0.52899E+04,&
0.53407E+04, 0.53715E+04, 0.53961E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(30,n),n=1,30)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.25662E+04, 0.27528E+04, 0.29885E+04,&
0.32049E+04, 0.34630E+04, 0.37386E+04, 0.39938E+04, 0.43087E+04, 0.46257E+04,&
0.47146E+04, 0.47332E+04, 0.47456E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(31,n),n=1,31)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.22451E+04, 0.24161E+04, 0.26320E+04,&
0.28083E+04, 0.30211E+04, 0.32637E+04, 0.34876E+04, 0.37542E+04, 0.40379E+04,&
0.41702E+04, 0.41808E+04, 0.41856E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarch(32,n),n=1,32)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.19632E+04, 0.21059E+04, 0.22858E+04,&
0.24414E+04, 0.26294E+04, 0.28437E+04, 0.30406E+04, 0.32718E+04, 0.35204E+04,&
0.36945E+04, 0.36999E+04, 0.37002E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarch(33,n),n=1,33)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.17259E+04, 0.18468E+04, 0.19993E+04,&
0.21369E+04, 0.23035E+04, 0.24931E+04, 0.26669E+04, 0.28529E+04, 0.30716E+04,&
0.32771E+04, 0.32793E+04, 0.32770E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(34,n),n=1,34)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.15113E+04, 0.16153E+04, 0.17464E+04,&
0.18684E+04, 0.20152E+04, 0.21782E+04, 0.23317E+04, 0.24914E+04, 0.26787E+04,&
0.28914E+04, 0.29100E+04, 0.29064E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(35,n),n=1,35)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.13164E+04, 0.14088E+04, 0.15253E+04,&
0.16336E+04, 0.17617E+04, 0.18942E+04, 0.20297E+04, 0.21858E+04, 0.23376E+04,&
0.25247E+04, 0.25847E+04, 0.25805E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(36,n),n=1,36)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.11508E+04, 0.12331E+04, 0.13367E+04,&
0.14329E+04, 0.15460E+04, 0.16595E+04, 0.17796E+04, 0.19113E+04, 0.20411E+04,&
0.22065E+04, 0.22975E+04, 0.22933E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/

data pi/3.141592654/

!calculating the mass categories (c.g.s) with the lowest diameter
!of 3.125 microns and mass doubling every bin
diam(1)=1.5625*2.e-04
x(1)=pi/6.*diam(1)**3.
do l=2,15
   x(l)=2.*x(l-1)
   diam(l)=(6./pi*x(l))**(0.333333333333)
enddo
do l=16,ibins+1
   x(l)=2.*x(l-1)
   !uses power law values in m.k.s. and converts back to c.g.s.
   diam(l)=100 * ((x(l) / 1000.0 / cfmass) ** (1./pwmass))
enddo

!(Greenan & List) collection kernel by Saleeby 12/07/04 for hail/cloud droplet
! weighted for (x+y)
do i=1,ibins
   do j=1,i
      akbarch(i,j)=abarch(i,j)
      if(akbarch(i,j).lt.0.) akbarch(i,j)=0.
      akbarch(36,i)=0.
   enddo
enddo

!Invert here for use in routine 'sxy'
do j=1,ibins
   do i=1,j
      akbarch(i,j)=akbarch(j,i)
   enddo
enddo
return
end subroutine data_ch

!******************************************************************************
subroutine initg2mode(r1,r2,r3,n1,n2,n3,gnu1,gnu2,gnu3,diam,x,amk,ank  &
       ,ank1,amk1,ank2,amk2,ank3,amk3,ithresh,ithresh1,ibins)
!!use micphys
implicit none
integer :: i,ibins,ithresh,ithresh1
real :: r1,r2,r3,n1,n2,n3,gnu1,gnu2,gnu3,dn1,dn2,dn3,trunc,fac1,fac2,pi  &
       ,dmean,sum,sumn,ex1
real, dimension(ibins+1) :: x,diam
real, dimension(ibins) :: amk,ank,amk1,ank1,amk2,ank2,amk3,ank3

real :: gammln,gammp

! * *
! * Initial triple gamma distribution: n(D) = n1(D) + n2(D) + n3(D)
! * *

data pi/3.141592654/
data ex1/0.333333333/

! * *
! * gamma spectrum
! * *

do i=1,ibins
  ank1(i)=0.
  amk1(i)=0.
  ank2(i)=0.
  amk2(i)=0.
  ank3(i)=0.
  amk3(i)=0.
enddo

!*******************************************
! CLOUD
!*******************************************
dmean = (6.*r1/(pi*n1))**ex1   !mass mean diam
dn1 = dmean * (exp(gammln(gnu1) - gammln(gnu1+3.))) ** ex1

do i=1,ithresh
 fac1=gammp(gnu1,diam(i)/dn1)
 fac2=gammp(gnu1,diam(i+1)/dn1)
 trunc=fac2-fac1
 ank1(i)=n1*trunc

 fac1=gammp(gnu1+3.,(diam(i)/dn1))
 fac2=gammp(gnu1+3.,(diam(i+1)/dn1))
 trunc=fac2-fac1
 amk1(i)=r1*trunc
enddo

!This is a single bin patch for cases of truncation of cloud1
ank1(ithresh-1) = ank1(ithresh-1) + ank1(ithresh)
amk1(ithresh-1) = amk1(ithresh-1) + amk1(ithresh)
ank1(ithresh) = 0.0
amk1(ithresh) = 0.0

!********************************************
! CLOUD 2
!********************************************
dmean = (6. * r2 / (pi * n2)) ** ex1   !mass mean diam
dn2 = dmean * (exp(gammln(gnu2) - gammln(gnu2+3.))) ** ex1

do i=ithresh,ithresh1+1
 fac1=gammp(gnu2,diam(i)/dn2)
 fac2=gammp(gnu2,diam(i+1)/dn2)
 trunc=fac2-fac1
 ank2(i)=n2*trunc

 fac1=gammp(gnu2+3.,(diam(i)/dn2))
 fac2=gammp(gnu2+3.,(diam(i+1)/dn2))
 trunc=fac2-fac1
 amk2(i)=r2*trunc
enddo

!This is a single bin patch for cases of truncation of cloud2
ank2(ithresh1) = ank2(ithresh1) + ank2(ithresh1+1)
amk2(ithresh1) = amk2(ithresh1) + amk2(ithresh1+1)
ank2(ithresh1+1) = 0.0
amk2(ithresh1+1) = 0.0

!********************************************
! RAIN
!********************************************
dmean = (6. * r3 / (pi * n3)) ** ex1  !mass mean diam
dn3 = dmean * (exp(gammln(gnu3) - gammln(gnu3+3.))) ** ex1

do i=ithresh1+1,ibins
 fac1=gammp(gnu3,diam(i)/dn3)
 fac2=gammp(gnu3,diam(i+1)/dn3)
 trunc=fac2-fac1
 ank3(i)=n3*trunc

 fac1=gammp(gnu3+3.,(diam(i)/dn3))
 fac2=gammp(gnu3+3.,(diam(i+1)/dn3))
 trunc=fac2-fac1
 amk3(i)=r3*trunc
enddo

!**********************************************
! SUM THE NUMBER AND MASS OF EACH DISTRIBUTION
!**********************************************
do i=1,ibins
   ank(i)=ank1(i)+ank2(i)+ank3(i)
   amk(i)=amk1(i)+amk2(i)+amk3(i)
enddo

return
end subroutine initg2mode

!******************************************************************************
subroutine sumn(ank,amk,imin,imax,ibins,sun,sum)
implicit none
integer :: imin,imax,ibins
real :: sun,sum
real, dimension(ibins) :: ank,amk

integer :: i

  sum=0.
  sun=0.

  do i=imin,imax
   sun=sun+ank(i)
   sum=sum+amk(i)
  enddo

return
end subroutine sumn

!******************************************************************************
subroutine initg1mode(r1,r3,n1,n3,gnu1,gnu3,diam,x,amk,ank  &
       ,ank1,amk1,ank3,amk3,ithresh,ibins,hyd)
!!use micphys
implicit none
integer :: i,ibins,ithresh,hyd
real :: r1,r3,n1,n3,gnu1,gnu3,dn1,dn3,trunc,fac1,fac2,pi  &
       ,dmean,sum,sumn,ex1,gammln,gammp
real, dimension(ibins+1) :: x,diam
real, dimension(ibins) :: amk,ank,amk1,ank1,amk3,ank3

data pi/3.141592654/
data ex1/0.333333333/

do i=1,ibins
  ank1(i)=0.
  amk1(i)=0.
  ank3(i)=0.
  amk3(i)=0.
enddo

!*******************************************
! CLOUD
!*******************************************
dmean = (6.*r1/(pi*n1))**ex1   !mass mean diam
dn1 = dmean * (exp(gammln(gnu1) - gammln(gnu1+3.))) ** ex1

do i=1,ithresh
 fac1=gammp(gnu1,diam(i)/dn1)
 fac2=gammp(gnu1,diam(i+1)/dn1)
 trunc=fac2-fac1
 ank1(i)=n1*trunc

 fac1=gammp(gnu1+3.,(diam(i)/dn1))
 fac2=gammp(gnu1+3.,(diam(i+1)/dn1))
 trunc=fac2-fac1
 amk1(i)=r1*trunc
enddo

!********************************************
! RAIN, SNOW, AGGREGATES, GRAUPEL, or HAIL
!********************************************
dmean = 100 * (r3 / n3 / 1000. / cfmas(hyd)) ** (1./pwmas(hyd))
dn3 = dmean * (exp(gammln(gnu3) - gammln(gnu3+pwmas(hyd)))) ** (1./pwmas(hyd))

do i=ithresh+1,ibins
 fac1=gammp(gnu3,diam(i)/dn3)
 fac2=gammp(gnu3,diam(i+1)/dn3)
 trunc=fac2-fac1
 ank3(i)=n3*trunc

 fac1=gammp(gnu3+pwmas(hyd),(diam(i)/dn3))
 fac2=gammp(gnu3+pwmas(hyd),(diam(i+1)/dn3))
 trunc=fac2-fac1
 amk3(i)=r3*trunc
enddo

!**********************************************
! SUM THE NUMBER AND MASS OF EACH DISTRIBUTION
!**********************************************
do i=1,ibins
   ank(i)=ank1(i)+ank3(i)
   amk(i)=amk1(i)+amk3(i)
enddo

return
end subroutine initg1mode



!******************************************************************************

subroutine tabhab()

!!use micphys

implicit none

integer, parameter :: nhab=0
integer :: it,is

if (nhab .eq.  0) print*,'VARIABLE HABIT PREDICTION'
if (nhab .eq.  3) print*,'ASSUMED HABIT IS COLUMNS'
if (nhab .eq.  8) print*,'ASSUMED HABIT IS HEX PLATES'
if (nhab .eq.  9) print*,'ASSUMED HABIT IS DENDRITES'
if (nhab .eq. 10) print*,'ASSUMED HABIT IS NEEDLES'
if (nhab .eq. 11) print*,'ASSUMED HABIT IS ROSETTES'
!c    if (nhab .eq.  x) print*,'ASSUMED HABIT IS SPHERES'

! nt is temp, ns = satur (liq)

do it = 1,31
   do is = 1,100
      if (nhab .eq. 0) then
         if (it .ge. 0 .and. it .le. 2) then
            if (is .le. 95) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 8
               jhabtab(it,is,2) = 12
            endif
         else if(it .gt. 2 .and. it .le. 4) then
            if (is .lt. 90) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 8
               jhabtab(it,is,2) = 12
            endif
         else if(it .gt. 4 .and. it .le. 6) then
            if (is .lt. 85) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 10
               jhabtab(it,is,2) = 14
            endif
         else if(it .gt. 6 .and. it .le. 9) then
            if (is .lt. 90) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 10
               jhabtab(it,is,2) = 14
            endif
         else if(it .gt. 9 .and. it .le. 22) then
            if (is .lt. 90) then
               jhabtab(it,is,1) = 8
               jhabtab(it,is,2) = 12
            else
               jhabtab(it,is,1) = 9
               jhabtab(it,is,2) = 13
            endif
         elseif(it .gt. 22 .and. it .le. 30) then
            if (is .lt. 80) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 10
               jhabtab(it,is,2) = 14
            endif
         elseif(it .gt. 30) then
            if (is .lt. 90) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 11
               jhabtab(it,is,2) = 15
            endif
         endif
      else
         jhabtab(it,is,1) = nhab
         jhabtab(it,is,2) = nhab + 4
         if (nhab .eq. 3) jhabtab(it,is,2) = 4
      endif
   enddo
enddo
return
end subroutine tabhab
!###########################################################################

subroutine thrmstr(m1,k1,k2,lpw,pp,thp,theta,pi0,rtp,rv,i,j)

use rconstants
!!use micphys

implicit none

integer :: m1,i,j,k,lcat,lpw
real :: fracliq,tcoal,tairstr
integer, dimension(11) :: k1,k2
real, dimension(m1) :: pp,thp,theta,pi0,rtp,rv

!OVER WHOLE COLUMN
do k = lpw,m1
   pitot(k) = pi0(k) + pp(k)
   press(k) = p00 * (pitot(k) * cpi) ** cpor
   tair(k) = theta(k) * pitot(k) * cpi
enddo
!LOWER VAPOR LAYER
do k = 1,k1(11)-1
   theta(k) = thp(k)
   rv(k) = rtp(k)
enddo
!UPPER VAPOR LAYER
do k = k2(11)+1,m1
   theta(k) = thp(k)
   rv(k) = rtp(k)
enddo
!ICE AND LIQUID LAYER
do k = k1(11),k2(11)
   til(k) = thp(k) * pitot(k) * cpi
   rliq(k) = 0.
   rice(k) = 0.
enddo
!CLOUD AND RAIN
do lcat = 1,2
   do k = k1(lcat),k2(lcat)
      rliq(k) = rliq(k) + rx(k,lcat)
   enddo
enddo
!DRIZZLE MODE
do lcat = 8,8
   do k = k1(lcat),k2(lcat)
      rliq(k) = rliq(k) + rx(k,lcat)
   enddo
enddo
!PRISTINE ICE, SNOW, AGGREGATES
do lcat = 3,5
   do k = k1(lcat),k2(lcat)
      rice(k) = rice(k) + rx(k,lcat)
   enddo
enddo
!GRAUPEL AND HAIL
!qtc gets temperature and fraction of liquid from q
do lcat = 6,7
   do k = k1(lcat),k2(lcat)
      call qtc(qx(k,lcat),tcoal,fracliq)
      rliq(k) = rliq(k) + rx(k,lcat) * fracliq
      rice(k) = rice(k) + rx(k,lcat) * (1. - fracliq)
   enddo
enddo
!FOR ICE/LIQUID LAYER get hydro heat and sat vap mix ratio
do k = k1(11),k2(11)
   qhydm(k) = alvl * rliq(k) + alvi * rice(k)
   rvstr(k) = rtp(k) - rliq(k) - rice(k)
   sa(k,1) = til(k) * qhydm(k) / (1.e-12 + rliq(k) + rice(k))
enddo
!FOR ICE/LIQUID LAYER get temp of sat air in Celcius
do k = k1(11),k2(11)
   if (tair(k) .gt. 253.) then
      tairstr = 0.5 * (til(k)  &
         + sqrt(til(k) * (til(k) + cpi4 * qhydm(k))))
      sa(k,1) = sa(k,1) * cpi / (2. * tairstr - til(k))
   else
      tairstr = til(k) * (1. + qhydm(k) * cp253i)
      sa(k,1) = sa(k,1) * cp253i
  endif
  tairstrc(k) = tairstr - 273.15
enddo

return
end subroutine thrmstr

!******************************************************************************

subroutine diffprep(m1,lcat,k1,k2,rv,dn0,i,j,mynum)

use rconstants
!!use micphys

implicit none

integer :: m1,lcat,k1,k2,i,j,k,mynum,if1,if4,if6,if8,lhcat
real :: fre,scdei
real, dimension(m1) :: rv,dn0

!CODE BASED ON WALKO ET AL 2000
!EFFICIENT COMPUTATION OF VAPOR AND HEAT DIFFUSION BETWEEN HYDROMETEORS
!IN A NUMERICAL MODEL

!DETERMINES WHETHER TO CALCULATE USING LIQUID OR ICE "SA" ARRAYS
if ((lcat .le. 2) .or. (lcat .eq. 8)) then
   if1 = 1
   if4 = 4
   if6 = 6
   if8 = 8
else
   if1 = 2
   if4 = 5
   if6 = 7
   if8 = 9
endif

do k = k1,k2
   lhcat = jhcat(k,lcat)

   if (rx(k,lcat) .lt. 1.e-12) go to 229

   fre = frefac1(lhcat) * emb(k,lcat) ** pwmasi(lhcat)  &
      + rdynvsci(k) * frefac2(lhcat) * emb(k,lcat) ** cdp1(lhcat)

   sb(k,lcat) = cx(k,lcat) * dn0(k) * fre * pi4dt
   su(k,lcat) = vapdif(k) * sb(k,lcat)

!zero for rain,graupel,hail
   sd(k,lcat) = sh(k,lcat) * rx(k,lcat)
   se(k,lcat) = su(k,lcat) * sa(k,if6) + sb(k,lcat) * thrmcon(k)
   sf(k,lcat) = su(k,lcat) * sl(if1) - sb(k,lcat) * sa(k,2)
!q term used for rain,graupel,hail
   sg(k,lcat) = su(k,lcat) * sa(k,if8) + sb(k,lcat) * sa(k,3)  &
              + sj(lcat) * qr(k,lcat)
!reduces to 1/se if not rain,graupel,hail
   scdei = 1. / (sc(if1) * sd(k,lcat) + se(k,lcat))
   ss(k,lcat) = sf(k,lcat) * scdei
!reduces to sg*scdei if not all-liquid species
   sw(k,lcat) = (sg(k,lcat) - sk(if1) * sd(k,lcat)) * scdei
   ttest(k,lcat) = ss(k,lcat) * rv(k) + sw(k,lcat)

229    continue

enddo

!FOR ALL ICE HYDROS, "SM" IS 1
!IF PRISTINE,SNOW,AGG ABOVE OR= ZERO, "SH" IS 1
if (lcat .ge. 3 .and. lcat .le. 5) then
   do k = k1,k2
      if (rx(k,lcat) .lt. 1.e-12) go to 228
      if (ttest(k,lcat) .ge. 0.) then
         sm(k,lcat) = 0.
         sh(k,lcat) = 1.
         sd(k,lcat) = sh(k,lcat) * rx(k,lcat)
         scdei = 1. / (sc(if1) * sd(k,lcat) + se(k,lcat))
         ss(k,lcat) = sf(k,lcat) * scdei
         sw(k,lcat) = (sg(k,lcat) - sk(if1) * sd(k,lcat)) * scdei
      else
         sm(k,lcat) = 1.
      endif
228        continue
   enddo
endif

!FOR MIXED-PHASE HYDROMETEORS, "SM" IS 0
if (lcat .ge. 6 .and. lcat .le. 7) then
   do k = k1,k2
      if (rx(k,lcat) .lt. 1.e-12) go to 227
      if (ttest(k,lcat) .ge. 0.) then
         sm(k,lcat) = 0.
      else
         sm(k,lcat) = 1.
      endif
227        continue
   enddo
endif

do k = k1,k2
   if (rx(k,lcat) .lt. 1.e-12) go to 226
   sy(k,lcat) = rvsrefp(k,if1) * sm(k,lcat) * sw(k,lcat) - sa(k,if4)
   sz(k,lcat) = 1. - rvsrefp(k,if1) * ss(k,lcat) * sm(k,lcat)
   sumuy(k) = sumuy(k) + su(k,lcat) * sy(k,lcat)
   sumuz(k) = sumuz(k) + su(k,lcat) * sz(k,lcat)

226      continue
enddo

return
end subroutine diffprep

!******************************************************************************

subroutine vapdiff (m1,kf1,kf2,rv,i,j,mynum)

!!use micphys

implicit none

integer :: m1,kf1,kf2,i,j,k,mynum
real, dimension(m1) :: rv

do k = kf1,kf2
   rv(k) = (rvstr(k) + sumuy(k)) / (1.0 + sumuz(k))
enddo

return
end subroutine vapdiff

!******************************************************************************

subroutine vapflux(m1,lcat,i,j,mynum,k1,k2,dn0,rv &
       ,vapliq,vapice &
       ,vapcld,vaprain,vappris,vapsnow,vapaggr,vapgrau,vaphail,vapdriz &
       ,vapliqt,vapicet &
       ,vapcldt,vapraint,vapprist,vapsnowt,vapaggrt,vapgraut,vaphailt,vapdrizt)

!!use micphys

implicit none

integer :: m1,lcat,i,j,k,mynum,k1,k2,if1,if4
real :: rxx
real, dimension(m1) :: dn0,rv
real, dimension(m1) :: vapliq,vapice,vapcld,vaprain,vappris,vapsnow &
                      ,vapaggr,vapgrau,vaphail,vapdriz &
                      ,vapliqt,vapicet,vapcldt,vapraint,vapprist,vapsnowt &
                      ,vapaggrt,vapgraut,vaphailt,vapdrizt

!UPDATES MIXING RATIO AND HEAT DUE TO FLUX OF VAPOR
!ALSO LINKED TO WALKO ET AL 2000

if ((lcat .le. 2) .or. (lcat .eq. 8)) then
   if1 = 1
   if4 = 4
else
   if1 = 2
   if4 = 5
endif

do k = k1,k2

   if (rx(k,lcat) .lt. 1.e-12) go to 229
   rxtemp = rx(k,lcat)
   tx(k,lcat) = (ss(k,lcat) * rv(k) + sw(k,lcat)) * sm(k,lcat)
   vap(k,lcat) = su(k,lcat) * (rv(k) + sa(k,if4) - rvsrefp(k,if1) * tx(k,lcat))

   if (vap(k,lcat) .gt. -rx(k,lcat)) then

      rxx = rx(k,lcat) + vap(k,lcat)

      if (sm(k,lcat) .gt. .5) then
         qx(k,lcat) = sc(if1) * tx(k,lcat) + sk(if1)
         qr(k,lcat) = qx(k,lcat) * rxx
      else
         qx(k,lcat) = (rv(k) * sf(k,lcat) + sg(k,lcat)  &
                    - tx(k,lcat) * se(k,lcat)) / sd(k,lcat)
         qx(k,lcat) = min(350000.,max(-100000.,qx(k,lcat)))
         qr(k,lcat) = qx(k,lcat) * rxx
      endif

   endif

!Do the following section if pristine ice totally melts: evaporate it too.

!Saleeby(10/5/06): CCN restoration source.
!If total evap occurs, transfer CCN mass within hydromet to CCN field and
!add CCN number that is equal to the number of totally evaporated droplets.
!If the evaporated hydromet is a liquid species, then # drops = # CCN added,
!except that CCN are in #/CC and drops are in #/KG, so convert. We calculate
!the proportion of evaporated mass and return the proportion of CCN mass
!contained with the hydrometeor species. We use the calculations of number
!concentration in "subroutine enemb" to determine the number to restore. We
!do not restore any CCN upon sublimation of ice species. We're assuming that
!when ice sublimes there is no way for previously dissolved CCN mass to fuse
!back together to form a full sized CCN that can activate and nucleate.
!Evaporated aerosol from cloud1 goes to CCN and cloud2 or rain goes to GCCN.

   rxferratio = 0.
   cxloss     = 0.

   !If we are doing FULL EVAPORATION
   if ((lcat .eq. 3 .and. qx(k,lcat) .gt. 330000.) .or.  &
      vap(k,lcat) .le. -rx(k,lcat)) then

      sumuy(k) = sumuy(k) - su(k,lcat) * sy(k,lcat)
      sumuz(k) = sumuz(k) - su(k,lcat) * sz(k,lcat)
      sumvr(k) = sumvr(k) + rx(k,lcat)
      rv(k) = (rvstr(k) + sumuy(k) + sumvr(k)) / (1.0 + sumuz(k))
      vap(k,lcat) = - rx(k,lcat)
      tx(k,lcat) = 0.
      rx(k,lcat) = 0.
      qx(k,lcat) = 0.
      qr(k,lcat) = 0.
      cxloss     = cx(k,lcat)
      cx(k,lcat) = 0.
      rxferratio = 1.

   !If we are doing VAPOR GROWTH or PARTIAL EVAPORATION
   else

      if(vap(k,lcat) .lt. 0.)then
        !Compute losses from partial evaporation
        fracmass = min(1.,-vap(k,lcat) / rx(k,lcat))
        cxloss = cx(k,lcat) * enmlttab( int(200.*fracmass)+1, jhcat(k,lcat) )
        cx(k,lcat) = cx(k,lcat) - cxloss
        rxferratio = rmlttab(int(200.*fracmass)+1,jhcat(k,lcat))
      endif

      rx(k,lcat) = rxx
   endif !if (full evaporation occurs) or (vapor growth or partial evap)

   !Attempt to restore aerosols if conditions are met
   if(iccnlev>=2 .and. vap(k,lcat)<0.0 .and. rxferratio > .0001) then
   if(lcat<=2.or.lcat==8) then
       ccnnum  = cxloss / 1.e6 * dn0(k)
       ccnmass = cnmhx(k,lcat) * rxferratio
       cnmhx(k,lcat) = cnmhx(k,lcat) - ccnmass
       !For approximation, test rg based on CCN sulfate constants
       rg=(0.02523*ccnmass/ccnnum)**.3333
       if(rg>0.96e-4 .and. jnmb(lcat)>=5) then
         gccmx(k) = gccmx(k) + ccnmass
         gccnx(k) = gccnx(k) + ccnnum
       else
         cccmx(k) = cccmx(k) + ccnmass
         cccnx(k) = cccnx(k) + ccnnum
       endif
       if(lcat==1.and.rxferratio>.01.and.ccnnum>1.0.and.(rg>96.e-6 .or. rg<1.0e-6))then
        print*,'Bad rg restore',k,rg,rxferratio
        print*,'      ',ccnnum,ccnmass,vap(k,lcat)
        print*,'      ',cccnx(k),cccmx(k)
        print*,'      ',cnmhx(k,lcat)
       endif
   endif
   endif

   !Vapor deposition and evaporation budgets for all species
   if(imbudget >= 1) then
     if(lcat.eq.1 .or. lcat.eq.2 .or. lcat.eq.8) &
       vapliq(k) = vapliq(k) + vap(k,lcat)*budget_scale
     if(lcat.ge.3 .and. lcat.le.7) &
       vapice(k) = vapice(k) + vap(k,lcat)*budget_scale
   endif
   if(imbudget == 2) then
     if(lcat==1) vapcld(k)  = vap(k,lcat)*budget_scale
     if(lcat==2) vaprain(k) = vap(k,lcat)*budget_scale
     if(lcat==3) vappris(k) = vap(k,lcat)*budget_scale
     if(lcat==4) vapsnow(k) = vap(k,lcat)*budget_scale
     if(lcat==5) vapaggr(k) = vap(k,lcat)*budget_scale
     if(lcat==6) vapgrau(k) = vap(k,lcat)*budget_scale
     if(lcat==7) vaphail(k) = vap(k,lcat)*budget_scale
     if(lcat==8) vapdriz(k) = vap(k,lcat)*budget_scale
   endif
   if(imbudtot >= 1) then
     if(lcat.eq.1 .or. lcat.eq.2 .or. lcat.eq.8) &
       vapliqt(k) = vapliqt(k) + vap(k,lcat)*budget_scalet
     if(lcat.ge.3 .and. lcat.le.7) &
       vapicet(k) = vapicet(k) + vap(k,lcat)*budget_scalet
   endif
   if(imbudtot == 2) then
     if(lcat==1) vapcldt(k)  = vapcldt(k)  + vap(k,lcat)*budget_scalet
     if(lcat==2) vapraint(k) = vapraint(k) + vap(k,lcat)*budget_scalet
     if(lcat==3) vapprist(k) = vapprist(k) + vap(k,lcat)*budget_scalet
     if(lcat==4) vapsnowt(k) = vapsnowt(k) + vap(k,lcat)*budget_scalet
     if(lcat==5) vapaggrt(k) = vapaggrt(k) + vap(k,lcat)*budget_scalet
     if(lcat==6) vapgraut(k) = vapgraut(k) + vap(k,lcat)*budget_scalet
     if(lcat==7) vaphailt(k) = vaphailt(k) + vap(k,lcat)*budget_scalet
     if(lcat==8) vapdrizt(k) = vapdrizt(k) + vap(k,lcat)*budget_scalet
   endif

229     continue

enddo
return
end subroutine vapflux

!******************************************************************************

subroutine psxfer(m1,k1,k2,dn0,i,j)

!!use micphys

implicit none

integer :: m1,k1,k2,i,j,k,lhcat,it
real :: embx,dn,xlim,dvap,dqr,dnum
real, dimension(m1) :: dn0
real :: old_m,old_c,old_r,prelim_m,delta_c, delta_r

do k = k1,k2

   if (vap(k,3) .gt. 0. .or. vap(k,4) .lt. 0.) then

      if (vap(k,3) .gt. 0.) then
         lhcat = jhcat(k,3)
         embx = max(1.e-12,rx(k,3)) / max(1.e-6,cx(k,3))
         dn = dnfac(lhcat) * embx ** pwmasi(lhcat)
         it = min(5000,max(1,nint(dn * 1.e6)))

         xlim = gam(it,3) * dps2 * (dps / dn) ** (gnu(3) - 1.)  &
            / (gamn1(3) * pwmas(lhcat) * dn ** 2)

         dvap = min(rx(k,3),  &
                    vap(k,3) * (xlim + gam(it,1) / gamn1(3)))
         dqr = dvap * qx(k,3)
         dnum = dvap * min(dpsmi(lhcat),1./embx)
      else
         lhcat = jhcat(k,4)
         embx = max(1.e-12,rx(k,4)) / max(1.e-6,cx(k,4))
         dn = dnfac(lhcat) * embx ** pwmasi(lhcat)
         it = min(5000,max(1,nint(dn * 1.e6)))

         xlim = gam(it,3) * dps2 * (dps / dn) ** (gnu(4) - 1.)  &
            / (gamn1(4) * pwmas(lhcat) * dn ** 2)

         dvap = max(-rx(k,4),vap(k,4) * xlim)
         dqr = dvap * qx(k,4)
         dnum = dvap * max(dpsmi(lhcat),1./embx)
      endif

      rx(k,3) = rx(k,3) - dvap
      cx(k,3) = cx(k,3) - dnum
      qr(k,3) = qr(k,3) - dqr
      rx(k,4) = rx(k,4) + dvap
      cx(k,4) = cx(k,4) + dnum
      qr(k,4) = qr(k,4) + dqr

      !Carrio 2003: Better xfer pristine to snow to prevent subroutine
      !enemb from artificially creating pristine number concentration
      !when a bounds problem appears with the pristine ice cutoff size.
      !Compare preliminary calculation pristine mean mass with maximum.
      delta_r = 0.0
      prelim_m = rx(k,3) / cx(k,3)
      if (prelim_m .gt. emb1(3)) then
          old_m = (rx(k,3) + dvap) / (cx(k,3) + dnum)
          old_c = cx(k,3) + dnum
          old_r = rx(k,3) + dvap
          delta_r = rx(k,3) - old_c * emb1(3)
          delta_c = delta_r / emb1(3) ! delta_c only adds to snow
          rx(k,3) = rx(k,3) - delta_r
          cx(k,3) = old_c
          rx(k,4) = rx(k,4) + delta_r
          cx(k,4) = cx(k,4) + delta_c
      endif

   endif
enddo
return
end subroutine psxfer

!******************************************************************************

subroutine newtemp(m1,kf1,kf2,rv,theta,i,j)

use rconstants
!!use micphys

implicit none

real rslf,rsif

integer :: m1,kf1,kf2,i,j,k
real, dimension(m1) :: rv,theta

do k = kf1,kf2
   tairc(k) = tairstrc(k) + sa(k,1) * (rvstr(k) - rv(k))
   tair(k)  = tairc(k) + 273.15
   theta(k) = tair(k) * cp / pitot(k)

   rvlsair(k) = rslf(press(k),tair(k))
   rvisair(k) = rsif (press(k),tair(k))
enddo

return
end subroutine newtemp
!******************************************************************************
end Module rams_microphysics_2M
