module rrtm_driv

     use node_mod, only: mynum ! INTENT(IN)
     use mem_rrtm, only: &
                i8,r8,aot_rrtm_lw,aot_rrtm_sw,sig,delsig,sigmid, &
                adjes,adjust,ch4ispresent,ch4pos,co2ispresent, &
                firsttime, iceflglw, iceflgsw,icld,idrv,inflglw, &
                liqflglw,nls,no2ispresent,no2pos,o2ispresent, &
                o2pos,co2pos,inflgsw,o3ispresent,o3pos,pi,rch4ar, &
                rco2ar,rh2oar,rno2ar,ro2ar,ro3ar,scon,bndmax_sw, &
                bndmax_lw,liqflgsw,bndmin_sw,bndmin_lw,initRRTM, &
                flip, irng, permuteseed

    use mem_tend, only: &
                tend

    use rrtmg_sw_rad, only: &
                rrtmg_sw

    use rrtmg_lw_rad, only: &
            rrtmg_lw

    use mcica_subcol_gen_sw, only: &
            mcica_subcol_sw

    use mcica_subcol_gen_lw, only: &
            mcica_subcol_lw
    
    use mem_grid    , only: &
         ngrid, time, dtlt, itime1, nzg, nzs, npatch, grid_g, & ! intent(in)
         nnzp, if_adap, zm, zt, naddsc, nzpmax, imonth1,      & ! intent(in)
         idate1, iyear1, centlat, centlon, ztop,dzm ,ngrids,  & ! intent(in)
         grid_vars, dzt

    use micphys     , only: &
         gnu, level, icloud, irain, ipris, & ! intent(in)
         isnow, iaggr, igraup, ihail,mcphys_type         ! intent(in)

    use mem_cuparm  , only: &
            cuparm_g, cuparm_vars, nnqparm  ! intent(in)

    use mem_basic   , only: &
        basic_g,basic_vars ! intent(inout)

    use mem_micro   , only: &
        micro_g,micro_vars ! intent(inout)

    use mem_radiate, only:        &
         ilwrtyp, iswrtyp,        & ! intent(in)
         radiate_g,radiate_vars,  & ! intent(inout)
         radfrq,                  & ! intent(in)
         ncall_i, prsnz, prsnzp,  & ! intent(inout)
         lonrad,radtun

    use parkind, only : &
                   im => kind_im, rb => kind_rb

    use rconstants  , only : &
                   cp,cpor,p00,stefan,cpi,pi180
    use parrrsw, only : &
                   nbndsw, ngptsw, naerec, jpband

    ! 2600-3250, 3250-4000, 4000-4650, 4650-5150, 5150-6150, 6150-7700,
    ! 7700-8050, 8050-12850, 12850-16000, 16000-22650, 22650-29000,
    !29000-38000, 38000-50000, and 820-2600

    use mem_leaf    , only: &
                   leaf_g, isfcl, albedo ! intent(in)

    use moddateutils, only: &
                   julday

    !use chem1_list , only: !nspecies,spc_name,o3,co2,no2,h2o,ch4,o2

    use mem_chem1, only: &
                   chem1_g, chemistry

    use ref_sounding, only: &
                   pi01dn ! intent(in)

    use parrrtm, only : &
                   nbndlw, ngptlw, maxxsec, mxmol

    use rrtmg_lw_cldprop, only: &
            cldprop

    use rrtmg_sw_cldprop, only: &
                   cldprop_sw

    use mem_tuv     , only: &
            carma_tuv, tuv2carma

    !use dump, only:

    use mem_grell_param, ONLY : maxiens ! INTENT(IN)

    use mem_scratch1_grell, only: &
                  ierr4d,xmb4d,zup5d,clwup5d,up_massdetr5d

    use leaf_coms, only: &
                  rlonga_a,rlonga_gs,rlonggs_v,rlonga_v,   &
                  rlongv_a,rlongv_gs,rlonggs_a,rshort_g,   &
                  rshort_v,snowfac,vf,tempk,slcpd,fracliq, &
                  emisv,emisg,slmsts,rshort_s

    use teb_spm_start, only: &
            teb_spm ! INTENT(IN)

    USE mem_carma, only: carma_aotMap

    implicit none

    real(kind=r8), allocatable, dimension(:,:,:) :: ozone
    real(kind=r8) :: ozsig(18)
    integer, parameter :: nlm_getoz=18
    logical :: first_getoz,inter_getoz
    integer :: year_getoz,mon_getoz
    integer, parameter :: nl=37
    integer, parameter :: ns=4
    include "aerosol_setup.f90"

    private
    public rrtm_driver

contains


    subroutine rrtm_driver(mzp, mxp, myp, ia, iz, ja, jz, mynum)

      integer, intent(in) :: mzp, mxp, myp, ia, iz, ja, jz, mynum
      !
      real,dimension(mzp,mxp,myp) :: lwl,iwl
      real,dimension(mxp,myp) :: rain
      real :: patch_area,can_temp,veg_temp,leaf_class,veg_height
      real :: veg_fracarea,veg_albedo,sfcwater_nlev,rshort,rlong
      real :: albedt,rlongup,cosz,hrAngleLocal,cdec
      real, dimension(nzg) :: soil_energy,soil_water,soil_text
      real, dimension(nzs) :: sfcwater_energy,sfcwater_depth
      real :: maxCloud_fraction
      real :: solfac
      integer :: ip,i,j,jday
      integer :: icount = 0
      integer, parameter :: ngpt = 141

      !- if not including radiation, return
      IF ((ilwrtyp + iswrtyp)==0) return

      icount = icount + 1
      if (icount>ngpt) icount = 0
       
!-srf moved the update to the end of the routine.
!!    !--- apply radiative tendencies to model tendencies
!!    call tend_accum_rtm(mzp, mxp, myp, ia, iz, ja, jz)

      !--- check if it is time to recompute radiative tendency and fluxes
      !
      !--- radiation calculation is updated only every radfrq seconds
      IF ( (mod(time+.001, radfrq) < dtlt .or. time<0.001)) THEN

        !--- set radiation tendency for theta to zero
        radiate_g(ngrid)%fthrd(1:mzp,1:mxp,1:myp) = 0.0

	!- call routine to calculate cloud properties for RRTM/CARMA
        !--1st, set zero to the local arrays
         radiate_g(ngrid)%cloud_fraction=0.0;  rain=0.0; lwl=0.0; iwl =0.0

	 call  cloud_prop_rrtm(mzp, mxp, myp, ia, iz, ja, jz &
		    !-- output
		    ,radiate_g(ngrid)%cloud_fraction	  &
		    ,rain		  &
		    ,lwl		  &
		    ,iwl		  &
		    )

!-srf tuning section for cloud fraction and other parameters for radiation
         if(radtun /= 1.0) then
           radiate_g(ngrid)%cloud_fraction=  min(1.,radtun* radiate_g(ngrid)%cloud_fraction)
	   rain          =  radtun*rain
	   lwl           =  radtun*lwl
	   iwl           =  radtun*iwl
         endif



        !- Compute solar zenith angle [cosz(i,j)] & solar constant factr [solfac].
         call zen_rtm(imonth1, idate1, iyear1, time, itime1, centlat, centlon, &
                     lonrad, pi180, ia, iz, ja, jz, jday, solfac, hranglelocal, cdec,&
		     mynum)

        !- compute patch-averaged surface albedo [albedt(i,j)] and up longwave
        !- radiative flux [rlongup(i,j)].
        !DSM ---- In case of JULES surface scheme, rlongup and albedt are already provided
        if (isfcl /= 5 .or. time<0.001) then
            radiate_g(ngrid)%albedt  = 0.
            radiate_g(ngrid)%rlongup = 0.

            do ip = 1,npatch
               do j = 1,jz
                  do i = 1,iz

               !By now just copy the memory to local variables
               !in order to send to sfcrad
                   soil_energy    (1:nzg)=leaf_g(ngrid)%soil_energy(1:nzg,i,j,ip)
                   soil_water     (1:nzg)=leaf_g(ngrid)%soil_water (1:nzg,i,j,ip)
                   soil_text      (1:nzg)=leaf_g(ngrid)%soil_text  (1:nzg,i,j,ip)
                   sfcwater_energy(1:nzs)=leaf_g(ngrid)%sfcwater_energy(1:nzs,i,j,ip)
                   sfcwater_depth (1:nzs)=leaf_g(ngrid)%sfcwater_depth (1:nzs,i,j,ip)
                   patch_area            =leaf_g(ngrid)%patch_area   (i,j,ip)
                   can_temp              =leaf_g(ngrid)%can_temp     (i,j,ip)
                   veg_temp              =leaf_g(ngrid)%veg_temp     (i,j,ip)
                   leaf_class            =leaf_g(ngrid)%leaf_class   (i,j,ip)
                   veg_height            =leaf_g(ngrid)%veg_height   (i,j,ip)
                   veg_fracarea          =leaf_g(ngrid)%veg_fracarea (i,j,ip)
                   veg_albedo            =leaf_g(ngrid)%veg_albedo   (i,j,ip)
                   sfcwater_nlev         =leaf_g(ngrid)%sfcwater_nlev(i,j,ip)

                   rshort                =radiate_g(ngrid)%rshort (i,j)
                   rlong                 =radiate_g(ngrid)%rlong  (i,j)
                   albedt                =radiate_g(ngrid)%albedt (i,j)
                   rlongup               =radiate_g(ngrid)%rlongup(i,j)
                   cosz                  =radiate_g(ngrid)%cosz   (i,j)

                   !
                   call sfcrad_rtm(nzg, nzs, ip,             &
                             soil_energy,    soil_water,     &
                             soil_text  ,    sfcwater_energy,&
                             sfcwater_depth, patch_area,     &
                             can_temp,       veg_temp,       &
                             leaf_class,     veg_height,     &
                             veg_fracarea,   veg_albedo,     &
                             sfcwater_nlev,                  &
                             rshort, rlong,  albedt,         &
                             rlongup, cosz                   &
                             )

                  !Copy back albedo and rlongUP to memory
                   radiate_g(ngrid)%albedt (i,j)=albedt
                   radiate_g(ngrid)%rlongup(i,j)=rlongup

                  end do
              end do
          end do
        endif



       !- RRTM Radiation
        call radrrtmdrv(ia,iz,ja,jz,mxp,myp,mzp,mynum&
                       , radiate_g(ngrid)%cloud_fraction         &
                       ,rain                   &
                       ,lwl                    &
                       ,iwl                    &
		       ,icount                 &
		       ,ngpt                   &
                                               )
      ENDIF
      !--- apply radiative tendencies to model tendencies
      call tend_accum_rtm(mzp, mxp, myp, ia, iz, ja, jz)

    end subroutine rrtm_driver

! ****************************************************************************

    subroutine tend_accum_rtm(m1,m2,m3,ia,iz,ja,jz)
          integer, intent(in) :: m1, m2, m3, ia, iz, ja, jz

          ! local variables:
          integer :: i, j, k,ipos

          ipos=0
          do j=1,m3
             do i=1,m2
                do k=1,m1
                  ipos=ipos+1
                  tend%tht(ipos) = tend%tht(ipos) + radiate_g(ngrid)%fthrd(k,i,j)
                end do
             end do
          end do

    end subroutine tend_accum_rtm
! ****************************************************************************

    subroutine zen_rtm(imonth1, idate1, iyear1, time, itime1, centlat, centlon, &
                lonrad, pi180, ia, iz, ja, jz, jday, solfac, hrangle, cdec,mynum)

                ! arguments:
                integer, intent(in)  :: imonth1, idate1, iyear1, itime1,mynum
                real, intent(in)     :: time
                real, intent(in)     :: centlat(:), centlon(:)
                integer, intent(in)  :: lonrad
                real, intent(in)     :: pi180
                integer, intent(in)  :: ia, iz, ja, jz
                integer, intent(out) :: jday
                real, intent(out)    :: solfac
                real, intent(out)    :: hrangle
                real, intent(out)    :: cdec
                ! local variables:
                integer :: i, j
                real    :: sdec, declin, d0, d02, dayhr, radlat, cslcsd, snlsnd, gglon, &
                           dayhrr, eqt

                jday   = julday(imonth1, idate1, iyear1)
                jday   = jday + nint(time/86400.)
                !      sdec - sine of declination, cdec - cosine of declination
                declin = -23.5*cos(6.283/365.*(jday + 9))*pi180
                sdec   = sin(declin)
                cdec   = cos(declin)

                ! find the factor, solfac, to multiply the solar constant to correct
                ! for earth's varying distance to the sun.

                d0     = 6.2831853*float(jday-1)/365.
                d02    = d0*2.
                solfac = 1.000110 + 0.034221*cos(d0) + 0.001280*sin(d0) + &
                         0.000719*cos(d02) + 0.000077*sin(d02)

                ! find the hour angle, then get cosine of zenith angle.

                !ner_i - including solar time equation ("eqt" must be defined, it is a new variable)
                eqt = (0.000075 + 0.001868*cos(d0) - 0.032077*sin(d0) - 0.014615*cos(d02) &
                       - 0.040849*sin(d02))*1440/(2*3.141593)
                !ner_f - including solar time equation

                !-ner dayhr  = time/3600. + float(itime1/100) + float(mod(itime1,100))/60.
                dayhr = (time / 3600. + float(itime1/100) + float(mod(itime1,100)) / 60.) !&
		      !+ (radfrq/(2.*3600.))
                
		!-ner (radfrq/(2*3600)) - rad transfer shift half of radfrq(improving rad tendency representativity)

                do j = ja,jz
                   do i = ia,iz
                      radlat = grid_g(ngrid)%glat(i,j)*pi180
                      if (lonrad==0)      radlat = centlat(1)*pi180
                      if (radlat==declin) radlat = radlat + 1.e-5
                      cslcsd = cos(radlat)*cdec
                      snlsnd = sin(radlat)*sdec
                      gglon  = grid_g(ngrid)%glon(i,j)
                      if (lonrad==0)      gglon = centlon(1)

                      !ner_i new hour angle calculation
                      hrangle=((dayhr+(gglon/15.)-(12.-eqt/60.))*15./1.)*3.141593/180.
                      !          dayhrr    = mod(dayhr+gglon/15.+24., 24.)
                      !          hrangl    = 15.*(dayhrr - 12.)*pi180

                      radiate_g(ngrid)%cosz(i,j) = snlsnd + cslcsd*cos(hrangle)

                      !radiate_g(ngrid)%cosz(i,j) = min(radiate_g(ngrid)%cosz(i,j), 1.0)
                      !radiate_g(ngrid)%cosz(i,j) = max(radiate_g(ngrid)%cosz(i,j),-1.0)

                      radiate_g(ngrid)%cosz(i,j) = max(radiate_g(ngrid)%cosz(i,j),1.0e-10)

                   end do
                end do
    end subroutine zen_rtm
! ****************************************************************************

    subroutine sfcrad_rtm(nzg, nzs, ip,                                            &
             soil_energy, soil_water, soil_text, sfcwater_energy, sfcwater_depth,  &
             patch_area, can_temp, veg_temp, leaf_class, veg_height, veg_fracarea, &
             veg_albedo, sfcwater_nlev, rshort, rlong, albedt, rlongup, cosz,      &
             ! For TEB
             G_URBAN, ETOWN, ALBTOWN, TSTOWN                                       &
             !
             )

          ! Arguments:
          integer, intent(IN) :: nzg, nzs, ip
          real, intent(IN)    :: soil_energy(nzg)
          real, intent(IN)    :: soil_water(nzg)
          real, intent(IN)    :: soil_text(nzg)
          real, intent(IN)    :: sfcwater_energy(nzs)
          real, intent(IN)    :: sfcwater_depth(nzs)
          real, intent(IN)    :: patch_area
          real, intent(IN)    :: can_temp
          real, intent(IN)    :: veg_temp
          real, intent(IN)    :: leaf_class
          real, intent(IN)    :: veg_height
          real, intent(IN)    :: veg_fracarea
          real, intent(IN)    :: veg_albedo
          real, intent(IN)    :: sfcwater_nlev
          real, intent(IN)    :: rshort
          real, intent(IN)    :: rlong
          real, intent(INOUT) :: albedt
          real, intent(INOUT) :: rlongup
          real, intent(IN)    :: cosz
          ! for TEB
          real, pointer, optional :: G_URBAN, ETOWN, ALBTOWN, TSTOWN
          !
          ! Local Variables:
          integer :: k, nsoil, nveg, ksn
          real :: alb, vfc, fcpct, alg, rad, als, fractrans, absg, algs, emv, emgs, &
               gslong, vlong, alv
          real :: vctr32(nint(sfcwater_nlev)+10)
          ! This routine is called by the radiation parameterization and by leaf.
          ! It computes net surface albedo plus radiative exchange between the
          ! atmosphere, vegetation, and the snow/ground given previously computed
          ! downward longwave and shortwave fluxes from the atmosphere.
          ! Also computed are functions of snowcover that are required for the above
          ! radiation calculations as well as other calculations in leaf.

          ! The shortwave parameterizations are only valid if the cosine of the zenith
          ! angle is greater than .03 .  Water albedo from Atwater and Bell (1981)

          ! alg, als, and alv are the albedos of the ground, snow, and vegetation
          ! (als needs a better formula based on age of the surface snow).

          ! absg and vctr32 are the actual fractions of shortwave incident on snow
          ! plus ground that get absorbed by the ground and each snow layer,
          ! respectively.  They currently use the variable fractrans, which is the
          ! fraction of light transmitted through each layer based on mass per square
          ! meter.  algs is the resultant albedo from snow plus ground.

          if (ip==1) then

             if (cosz>0.03) then
                alb    = min(max(-.0139 + .0467*tan(acos(cosz)), 0.03), 0.999)
                albedt = albedt + patch_area*alb
             endif

             call qtk(soil_energy(nzg), tempk(nzg), fracliq(nzg))
             rlongup = rlongup + patch_area*stefan*tempk(nzg)**4

          elseif (isfcl==0) then

             albedt  = albedt + patch_area*albedo
             rlongup = rlongup + patch_area*stefan*can_temp**4

          else

             ! Diagnose soil temperature and liquid fraction

             do k=1,nzg
                nsoil = nint(soil_text(k))
                call qwtk(soil_energy(k), soil_water(k)*1.e3,  &
                     slcpd(nsoil), tempk(k), fracliq(k))
             enddo

             ! Diagnose snow temperature and the influence of snow covering veg.

             nveg    = nint(leaf_class)
             ksn     = nint(sfcwater_nlev)
             snowfac = 0.
             do k=1,ksn
                snowfac = snowfac + sfcwater_depth(k)
                call qtk(sfcwater_energy(k), tempk(k+nzg), fracliq(k+nzg))
             enddo
             snowfac = min(0.99, snowfac/max(0.001, veg_height))

             vf  = veg_fracarea*(1. - snowfac)
             vfc = 1. - vf

             ! Shortwave radiation calculations

             !srf-25-02-2005
             nsoil = nint(soil_text(nzg))
             !srf

             fcpct = soil_water(nzg)/slmsts(nsoil)
             if (fcpct>0.5) then
                alg = 0.14
             else
                alg = 0.31 - 0.34*fcpct
             endif
             alv = veg_albedo

             rad = 1.
             if (ksn>0) then

                ! als = .14 (the wet soil value) for all-liquid

                als = 0.5 - 0.36*fracliq(ksn + nzg)
                rad = 1.  - als
             endif
             do k=ksn,1,-1
                fractrans = exp(-20.*sfcwater_depth(k))
                vctr32(k) = rad*(1. - fractrans)
                rad = rad * fractrans
             enddo
             absg = (1. - alg)*rad
             algs = 1. - absg
             do k=ksn,1,-1
                algs        = algs - vctr32(k)
                rshort_s(k) = rshort*vfc*vctr32(k)
             enddo
             rshort_g = rshort*vfc*absg
             rshort_v = rshort*vf*(1. - alv + vfc*algs)
             !  rshort_a = rshort*(vf*alv + vfc*vfc*algs)

             alb = vf*alv + vfc*vfc*algs

             ! TEB
             if (TEB_SPM==1 .and. present(G_URBAN) .and. present(ALBTOWN)) then
                if (NINT(G_URBAN)==0) then
                   albedt = albedt + patch_area*alb
                else
                   albedt = albedt + patch_area*ALBTOWN
                endif
             else
                albedt = albedt + patch_area*alb
             endif

             ! Longwave radiation calculations

             emv  = emisv(nveg)
             emgs = emisg(nsoil)
             if (ksn>0) emgs = 1.0
             gslong = emgs*stefan*tempk(ksn+nzg)**4
             vlong  = emv*stefan*veg_temp**4

             rlonga_v  = rlong*vf*(emv + vfc*(1. - emgs))
             rlonga_gs = rlong*vfc*emgs
             rlongv_gs = vlong*vf*emgs
             rlongv_a  = vlong*vf*(2. - emgs - vf + emgs*vf)
             rlonggs_v = gslong*vf*emv
             rlonggs_a = gslong*vfc
             rlonga_a  = rlong*(vf*(1. - emv) + vfc*vfc*(1. - emgs))

             ! TEB
             if (TEB_SPM==1 .and. present(G_URBAN) .and. present(ETOWN) .and. &
                  present(TSTOWN)) then
                IF (NINT(G_URBAN)==0) THEN
                   rlongup = rlongup + patch_area*(rlongv_a + rlonggs_a + rlonga_a)
                else
                   rlongup = rlongup + patch_area*ETOWN*STEFAN*TSTOWN**4
                ENDIF
             else
                rlongup = rlongup + patch_area*(rlongv_a + rlonggs_a + rlonga_a)
             endif

             !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
             ! In case rlong is not computed, zero out all longwave fluxes other
             ! than rlongup.  [On the first timestep, radiative fluxes may not
             ! be available until microphysics is called, and zeroing these fluxes
             ! avoids the imbalance of having upward longwave without downward longwave
             ! in LEAF-2.  Also, this allows LEAF-2 to run without radiation for all
             ! timesteps, if desired for testing.]

             if (rlong<0.1) then
                rlonga_v  = 0.
                rlonga_gs = 0.
                rlongv_gs = 0.
                rlongv_a  = 0.
                rlonggs_v = 0.
                rlonggs_a = 0.
                rlonga_a  = 0.
             endif

          endif

   end subroutine sfcrad_rtm

! ****************************************************************************
   subroutine radrrtmdrv(ia,iz,ja,jz,mxp,myp,mzp,mynum &
                        ,cloud_fraction           &
                        ,rain                     &
                        ,lwl                      &
                        ,iwl                      &
			,icount                   &
			,ngpt                     &
                                                  )
   use ccatt_start, only: &
                   ccatt
   use mem_basic, only: &
                   basic_g  ! INTENT(INOUT)
   use mem_grid, only: &
                   deltaxn, &
                   deltayn, &
                   zmn
   use optical, only: &
                   optProp

   use mem_carma, only: &
                  carma
   implicit none
  

   integer, intent(in) :: ia,iz,ja,jz
   integer, intent(in) :: icount
   integer, intent(in) :: ngpt
   
   integer, intent(in) :: mxp,myp,mzp,mynum
   real, intent(in), dimension(mzp,mxp,myp) :: cloud_fraction !cloud_fraction
   real, intent(in), dimension(    mxp,myp) :: rain !total rain water
   real, intent(in), dimension(mzp,mxp,myp) :: lwl !total cloud liquid water (kg/kg for carma and g/m2 for rrtm)
   real, intent(in), dimension(mzp,mxp,myp) :: iwl !total cloud ice water (kg/kg for carma and g/m2 for rrtm)

   integer,allocatable,dimension(:) :: ipos,jpos
   integer (kind=i8),allocatable,dimension(:)  :: imask
   logical, allocatable, dimension(:) :: aboveh

   integer :: nlay
   integer, parameter :: iplon=1
   
   real(kind=rb) :: taucloud_lw(mzp-1,nbndlw)

   real(kind=rb) :: taucloud_sw(mzp-1,jpband)
   real(kind=rb) :: tauc_lw(nbndlw,mzp-1)
   real(kind=rb) :: ssac_lw(nbndlw,mzp-1)
   real(kind=rb) :: asmc_lw(nbndlw,mzp-1)
   real(kind=rb) :: fsfc_lw(nbndlw,mzp-1)

   real(kind=rb) :: tauc_sw(nbndsw,mzp-1)
   real(kind=rb) :: ssac_sw(nbndsw,mzp-1)
   real(kind=rb) :: asmc_sw(nbndsw,mzp-1)
   real(kind=rb) :: fsfc_sw(nbndsw,mzp-1)
   
   

   real(kind=rb) :: taucldorig(mzp-1,jpband)

   real(kind=rb) :: ssacloud_sw(mzp-1,jpband)
   real(kind=rb) :: ssacloud_lw(mzp-1,jpband)

   real(kind=rb) :: asmcloud_sw(mzp-1,jpband)
   real(kind=rb) :: asmcloud_lw(mzp-1,jpband)

   real(kind=rb),parameter :: factorcor=1.0e-9_rb

   real(kind=rb), allocatable :: play(:,:)          ! layer pressures (hpa, mb)
                                                   !    dimensions: (ncol,nlay)
   real(kind=rb), allocatable :: plev(:,:)          ! interface pressures (hpa, mb)
                                                   !    dimensions: (ncol,nlay+1)
   real(kind=rb), allocatable :: tlay(:,:)          ! layer temperatures (k)
                                                   !    dimensions: (ncol,nlay)
   real(kind=rb), allocatable :: tlev(:,:)          ! interface temperatures (k)
                                                   !    dimensions: (ncol,nlay+1)
   real(kind=rb), allocatable :: tsfc(:)            ! surface temperature (k)
                                                   !    dimensions: (ncol)
   real(kind=rb), allocatable :: h2ovmr(:,:)        ! h2o volume mixing ratio
                                                   !    dimensions: (ncol,nlay)
   real(kind=rb), allocatable :: o3vmr(:,:)         ! o3 volume mixing ratio
                                                   !    dimensions: (ncol,nlay)
   real(kind=rb), allocatable :: co2vmr(:,:)        ! co2 volume mixing ratio
                                                   !    dimensions: (ncol,nlay)
   real(kind=rb), allocatable :: ch4vmr(:,:)        ! methane volume mixing ratio
                                                   !    dimensions: (ncol,nlay)
   real(kind=rb), allocatable :: n2ovmr(:,:)        ! nitrous oxide volume mixing ratio
                                                   !    dimensions: (ncol,nlay)
   real(kind=rb), allocatable :: o2vmr(:,:)         ! oxygen volume mixing ratio
                                                   !    dimensions: (ncol,nlay)
   real(kind=rb), allocatable :: cfc11vmr(:,:)      ! cfc11 volume mixing ratio
                                                   !    dimensions: (ncol,nlay)
   real(kind=rb), allocatable :: cfc12vmr(:,:)      ! cfc12 volume mixing ratio
                                                   !    dimensions: (ncol,nlay)
   real(kind=rb), allocatable :: cfc22vmr(:,:)      ! cfc22 volume mixing ratio
                                                   !    dimensions: (ncol,nlay)
   real(kind=rb), allocatable :: ccl4vmr(:,:)       ! ccl4 volume mixing ratio
                                                   !    dimensions: (ncol,nlay)
   real(kind=rb), allocatable :: emis(:,:)          ! surface emissivity
                                                   !    dimensions: (ncol,nbndlw)
   real(kind=rb), allocatable :: asdir(:)           ! uv/vis surface albedo direct rad
                                                   !    dimensions: (ncol)
   real(kind=rb), allocatable :: aldir(:)           ! near-ir surface albedo direct rad
                                                   !    dimensions: (ncol)
   real(kind=rb), allocatable :: asdif(:)           ! uv/vis surface albedo: diffuse rad
                                                   !    dimensions: (ncol)
   real(kind=rb), allocatable :: aldif(:)           ! near-ir surface albedo: diffuse rad
                                                   !    dimensions: (ncol)

   integer(kind=im) :: dyofyr          ! day of the year (used to get earth/sun
                                                   !  distance if adjflx not provided)
   real(kind=rb), allocatable :: coszen(:)         ! cosine of solar zenith angle
                                                   !    dimensions: (ncol)

                                                      !    in upward flux as a function of
                                                      !    surface temperature [0=off, 1=on]
                                                      !    0: normal forward calculation
                                                      !    1: normal forward calculation with
                                                      !       duflx_dt and duflxc_dt output

!!!Cloud properties from mcica_sw with dimensions: (ngptsw,ncol,nlay)

   !LFR For version 5.0 of rrtmg
   real(kind=rb), allocatable :: alpha_sw(:,:) ! cloud fraction decorrelation length (m)
   real(kind=rb), allocatable :: cldfmcl_sw(:,:,:)    ! cloud fraction from mcica_sw
   real(kind=rb), allocatable :: taucmcl_sw(:,:,:)    ! in-cloud optical depth 
   real(kind=rb), allocatable :: ssacmcl_sw(:,:,:)    ! in-cloud single scattering albedo 
   real(kind=rb), allocatable :: asmcmcl_sw(:,:,:)    ! in-cloud asymmetry parameter 
   real(kind=rb), allocatable :: fsfcmcl_sw(:,:,:)    ! in-cloud forward scattering fraction 
   real(kind=rb), allocatable :: ciwpmcl_sw(:,:,:)    ! in-cloud ice water path (g/m2) 
   real(kind=rb), allocatable :: clwpmcl_sw(:,:,:)    ! in-cloud liquid water path (g/m2) 
   real(kind=rb), allocatable :: reicmcl_sw(:,:)      ! cloud ice effective radius (microns) 
   real(kind=rb), allocatable :: relqmcl_sw(:,:)      ! cloud water drop effective radius (microns)

!!!Cloud properties from mcica_lw with dimensions: (ngptlw,ncol,nlay)

   real(kind=rb), allocatable :: cldfmcl_lw(:,:,:)    ! cloud fraction from mcica_lw
   real(kind=rb), allocatable :: taucmcl_lw(:,:,:)    ! in-cloud optical depth
   real(kind=rb), allocatable :: ssacmcl_lw(:,:,:)    ! in-cloud single scattering albedo
   real(kind=rb), allocatable :: asmcmcl_lw(:,:,:)    ! in-cloud asymmetry parameter
   real(kind=rb), allocatable :: fsfcmcl_lw(:,:,:)    ! in-cloud forward scattering fraction
   real(kind=rb), allocatable :: ciwpmcl_lw(:,:,:)    ! in-cloud ice water path (g/m2)
   real(kind=rb), allocatable :: clwpmcl_lw(:,:,:)    ! in-cloud liquid water path (g/m2)

!!!No subcolumns for these properties yet
   real(kind=rb), allocatable :: reicmcl_lw(:,:)      ! cloud ice effective radius (microns) 
   real(kind=rb), allocatable :: relqmcl_lw(:,:)      ! cloud water drop effective radius (microns)

!!!Cloud properties rrtmg with dimension dimensions: (ncol,nlay)
   real(kind=rb), allocatable :: cldf(:,:)            ! cloud fraction
   real(kind=rb), allocatable :: ciwp(:,:)            ! in-cloud ice water path (g/m2)
   real(kind=rb), allocatable :: clwp(:,:)            ! in-cloud liquid water path (g/m2)
   real(kind=rb), allocatable :: reic(:,:)            ! cloud ice effective radius (microns)
                                                      ! specific definition of reicmcl depends on setting of iceflgsw:
                                                      ! iceflgsw = 0: (inactive)
                                                      !
                                                      ! iceflgsw = 1: ice effective radius, r_ec, (ebert and curry, 1992),
                                                      !               r_ec range is limited to 13.0 to 130.0 microns
                                                      ! iceflgsw = 2: ice effective radius, r_k, (key, streamer ref. manual, 1996)
                                                      !               r_k range is limited to 5.0 to 131.0 microns
                                                      ! iceflgsw = 3: generalized effective size, dge, (fu, 1996),
                                                      !               dge range is limited to 5.0 to 140.0 microns
						      !               [dge = 0.9021 * r_ec - 7.0365] 
   real(kind=rb), allocatable :: relq(:,:)            ! cloud water drop effective radius (microns)

!!!Cloud properties rrtmg with dimensions : (nbndsw,ncol,nlay)

   real(kind=rb), allocatable :: taucl_sw(:,:,:)	     ! in-cloud optical depth
   real(kind=rb), allocatable :: ssacl_sw(:,:,:) 	   ! in-cloud single scattering albedo
   real(kind=rb), allocatable :: asmcl_sw(:,:,:) 	   ! in-cloud asymmetry parameter
   real(kind=rb), allocatable :: fsfcl_sw(:,:,:) 	   ! in-cloud forward scattering fraction

!!!Cloud properties rrtmg with dimensions: (nbndlw,ncol,nlay)
   !LFR For version 5.0 of rrtmg
   real(kind=rb), allocatable :: alpha_lw(:,:) ! cloud fraction decorrelation length (m)

   real(kind=rb), allocatable :: taucl_lw(:,:,:)	    ! in-cloud optical depth
   real(kind=rb), allocatable :: ssacl_lw(:,:,:)	    ! in-cloud single scattering albedo
   real(kind=rb), allocatable :: asmcl_lw(:,:,:)	    ! in-cloud asymmetry parameter
   real(kind=rb), allocatable :: fsfcl_lw(:,:,:)	    ! in-cloud forward scattering fraction

!!!Aerosol properties rrtmg_sw with dimensions: (ncol,nlay,nbndsw)

   real(kind=rb), allocatable :: tauaer_sw(:,:,:)     ! aerosol optical depth (iaer=10 only) (non-delta scaled)
   real(kind=rb), allocatable :: ssaaer_sw(:,:,:)     ! aerosol single scattering albedo (iaer=10 only) (non-delta scaled)
   real(kind=rb), allocatable :: asmaer_sw(:,:,:)     ! aerosol asymmetry parameter (iaer=10 only) (non-delta scaled)


!!!Aerosol properties rrtmg_sw with dimensions: (ncol,nlay,nbndsw)

   real(kind=rb), allocatable :: tauaer_lw(:,:,:)     ! aerosol optical depth (iaer=10 only) (non-delta scaled)
   real(kind=rb), allocatable :: ssaaer_lw(:,:,:)     ! aerosol single scattering albedo (iaer=10 only) (non-delta scaled)
   real(kind=rb), allocatable :: asmaer_lw(:,:,:)     ! aerosol asymmetry parameter (iaer=10 only) (non-delta scaled)

!!!Aerosol optical depth at 0.55 micron
   real(kind=rb), allocatable :: ecaer    (:,:,:)     ! aerosol optical depth at 0.55 micron (iaer=6 only) (non-delta scaled)

! -- output -----

   real(kind=rb), allocatable :: swuflx(:,:)       ! total sky shortwave upward flux (w/m2)
                                                   !    dimensions: (ncol,nlay+1)
   real(kind=rb), allocatable :: swdflx(:,:)       ! total sky shortwave downward flux (w/m2)
                                                   !    dimensions: (ncol,nlay+1)
   real(kind=rb), allocatable :: swhr(:,:)         ! total sky shortwave radiative heating rate (k/d)
                                                   !    dimensions: (ncol,nlay)
   real(kind=rb), allocatable :: swuflxc(:,:)      ! clear sky shortwave upward flux (w/m2)
                                                   !    dimensions: (ncol,nlay+1)
   real(kind=rb), allocatable :: swdflxc(:,:)      ! clear sky shortwave downward flux (w/m2)
                                                   !    dimensions: (ncol,nlay+1)
   real(kind=rb), allocatable :: swhrc(:,:)        ! clear sky shortwave radiative heating rate (k/d)
                                                      !    dimensions: (ncol,nlay)
   real(kind=rb), allocatable :: uflx(:,:)         ! total sky longwave upward flux (w/m2)
                                                      !        dimensions: (ncol,nlay+1)
   real(kind=rb), allocatable :: dflx(:,:)         ! total sky longwave downward flux (w/m2)
                                                      !        dimensions: (ncol,nlay+1)
   real(kind=rb), allocatable :: hr(:,:)           ! total sky longwave radiative heating rate (k/d)
                                                      !        dimensions: (ncol,nlay)
   real(kind=rb), allocatable :: uflxc(:,:)        ! clear sky longwave upward flux (w/m2)
                                                      !        dimensions: (ncol,nlay+1)
   real(kind=rb), allocatable :: dflxc(:,:)        ! clear sky longwave downward flux (w/m2)
                                                      !        dimensions: (ncol,nlay+1)
   real(kind=rb), allocatable :: hrc(:,:)           ! clear sky longwave radiative heating rate (k/d)
                                                      !        dimensions: (ncol,nlay)

   real(kind=rb), allocatable :: duflx_dt(:,:)
                                                    ! change in upward longwave flux (w/m2/k)
                                                    ! with respect to surface temperature
                                                    !    dimensions: (ncol,nlay+1)
   real(kind=rb), allocatable :: duflxc_dt(:,:)
                                                   ! change in clear sky upward longwave flux (w/m2/k)
                                                   ! with respect to surface temperature
                                                   !        dimensions: (ncol,nlay+1)
   integer :: nofcols,noc,np
   integer :: ilay,iSpc
   integer :: changeseed

   integer :: i,j,k,m
   real(kind=rb) :: picpi
   real(kind=r8) :: co2val,doy
   real(kind=rb),allocatable,dimension(:,:) :: p,t,q,fice,o3l,lmixr
   real(kind=rb),allocatable,dimension(:)   :: psurf,colrad
   character(len=2) :: cmzp
   integer :: itime(4)
   integer(kind=im) :: ncbandssw,ncbandslw

   integer(kind=im) :: iaer !LFR for the 5.0 version
   integer(kind=im), parameter :: isolvar=-1         ! Flag for solar variability method
   !#Solar variability scaling factors or indices (ISOLVAR=-1,1,2,3 only)
   !#            For ISOLVAR = 1:
   !#                    SOLVAR(1)    Facular (Mg) index amplitude scale factor
   !#                    SOLVAR(2)    Sunspot (SB) index amplitude scale factor
   !#
   !#            For ISOLVAR = 2:
   !#                    SOLVAR(1)    Facular (Mg) index as defined in the NRLSSI2 model;
   !#                                 used for modeling time-specific solar activity
   !#                    SOLVAR(2)    Sunspot (SB) index as defined in the NRLSSI2 model; 
   !#                                 used for modeling time-specific solar activity
   !#
   !#            For ISOLVAR = -1 or 3:
   !#                    SOLVAR(1:14) Band scale factors for modeling spectral variation of 
   !#                                 averaged solar cycle in each shortwave ban
   !real(kind=rb) :: indsolvar(2)
   !# Mg and SB index amplitude scale factors (isolvar=1), or
   !# Mg and SB indices as defined in the NRLSSI2 model (isolvar=2)
   !real(kind=rb) :: bndsolvar(nbndsw)
   !# Solar variability scale factors 
   !# for each shortwave band Dimensions: (nbndsw=14)
   !real(kind=rb) :: solcycfrac 
   !# Solar cycle fraction (0-1); fraction of the way through the mean 11-year
   !# cycle with 0.0 defined as the first day of year 1 and 1.0 defined as the
   !# last day of year 11 (ISOLVAR=1 only). Note that for the combined effect of
   !# the solar constant of 1360.85, and the mean facular brightening and sunspot 
   !# dimming components (without scaling), the minimum total solar irradiance of
   !# 1360.49 occurs at SOLCYCFRAC = 0.0265, and the maximum total solar irradiance 
   !# of 1361.34 occurs at SOLCYCFRAC = 0.3826. 

   
   iaer=10 !LFR for the 5.0 version

   nofcols=(iz-ia+1)*(jz-ja+1)
   nlay=mzp-1

   !- integer
   allocate(ipos  (nofcols))     ;ipos =0 !integer
   allocate(jpos  (nofcols))     ;jpos =0 !integer
   allocate(imask (nofcols))     ;imask=0 !integer
   !- real
   allocate(play  (nofcols,nlay)  ) ;play  =0.0
   allocate(plev  (nofcols,nlay+1)) ;plev  =0.0
   allocate(tlay  (nofcols,nlay)  ) ;tlay  =0.0
   allocate(tlev  (nofcols,nlay+1)) ;tlev  =0.0
   !
   allocate(h2ovmr(nofcols,nlay))    ;h2ovmr=0.0
   allocate(o3vmr (nofcols,nlay))    ;o3vmr =0.0
   allocate(co2vmr(nofcols,nlay))    ;co2vmr=0.0
   allocate(ch4vmr(nofcols,nlay))    ;ch4vmr=0.0
   allocate(n2ovmr(nofcols,nlay))    ;n2ovmr=0.0
   allocate(o2vmr (nofcols,nlay))    ;o2vmr =0.0
   !
   allocate(cfc11vmr(nofcols,nlay))  ; cfc11vmr=0.0
   allocate(cfc12vmr(nofcols,nlay))  ; cfc12vmr=0.0
   allocate(cfc22vmr(nofcols,nlay))  ; cfc22vmr=0.0
   allocate(ccl4vmr (nofcols,nlay))  ; ccl4vmr =0.0
   allocate(emis    (nofcols,nbndlw)); emis    =0.0
   !

   allocate(alpha_sw(nofcols,nlay)); alpha_sw=0.0
   allocate(cldfmcl_sw(ngptsw,nofcols,nlay)) ;cldfmcl_sw=0.0
   allocate(taucmcl_sw(ngptsw,nofcols,nlay)) ;taucmcl_sw=0.0
   allocate(ssacmcl_sw(ngptsw,nofcols,nlay)) ;ssacmcl_sw=0.0
   allocate(asmcmcl_sw(ngptsw,nofcols,nlay)) ;asmcmcl_sw=0.0
   allocate(fsfcmcl_sw(ngptsw,nofcols,nlay)) ;fsfcmcl_sw=0.0
   allocate(ciwpmcl_sw(ngptsw,nofcols,nlay)) ;ciwpmcl_sw=0.0
   allocate(clwpmcl_sw(ngptsw,nofcols,nlay)) ;clwpmcl_sw=0.0

   allocate(reicmcl_sw(nofcols,nlay)) ;reicmcl_sw=0.0
   allocate(relqmcl_sw(nofcols,nlay)) ;relqmcl_sw=0.0

   allocate(cldfmcl_lw(ngptlw,nofcols,nlay)) ;cldfmcl_lw=0.0
   allocate(taucmcl_lw(ngptlw,nofcols,nlay)) ;taucmcl_lw=0.0
   allocate(ssacmcl_lw(ngptlw,nofcols,nlay)) ;ssacmcl_lw=0.0
   allocate(asmcmcl_lw(ngptlw,nofcols,nlay)) ;asmcmcl_lw=0.0
   allocate(fsfcmcl_lw(ngptlw,nofcols,nlay)) ;fsfcmcl_lw=0.0
   allocate(ciwpmcl_lw(ngptlw,nofcols,nlay)) ;ciwpmcl_lw=0.0
   allocate(clwpmcl_lw(ngptlw,nofcols,nlay)) ;clwpmcl_lw=0.0

   allocate(reicmcl_lw(nofcols,nlay)) ;reicmcl_lw=0.0
   allocate(relqmcl_lw(nofcols,nlay)) ;relqmcl_lw=0.0

   allocate(cldf(nofcols,nlay)) 	     ;cldf=0.0
   allocate(ciwp(nofcols,nlay)) 	     ;ciwp=0.0
   allocate(clwp(nofcols,nlay)) 	     ;clwp=0.0
   allocate(reic(nofcols,nlay)) 	     ;reic=0.0
   allocate(relq(nofcols,nlay)) 	     ;relq=0.0

   allocate(taucl_sw(nbndsw,nofcols,nlay))    ;taucl_sw=0.0
   allocate(ssacl_sw(nbndsw,nofcols,nlay))    ;ssacl_sw=0.0
   allocate(asmcl_sw(nbndsw,nofcols,nlay))    ;asmcl_sw=0.0
   allocate(fsfcl_sw(nbndsw,nofcols,nlay))    ;fsfcl_sw=0.0

   allocate(alpha_lw(nofcols,nlay)); alpha_lw=0.0

   allocate(taucl_lw(nbndlw,nofcols,nlay))    ;taucl_lw=0.0
   allocate(ssacl_lw(nbndlw,nofcols,nlay))    ;ssacl_lw=0.0
   allocate(asmcl_lw(nbndlw,nofcols,nlay))    ;asmcl_lw=0.0
   allocate(fsfcl_lw(nbndlw,nofcols,nlay))    ;fsfcl_lw=0.0
					  
   allocate(tauaer_sw (nofcols,nlay,nbndsw)) ;tauaer_sw =0.0
   allocate(ssaaer_sw (nofcols,nlay,nbndsw)) ;ssaaer_sw =0.0
   allocate(asmaer_sw (nofcols,nlay,nbndsw)) ;asmaer_sw =0.0
   allocate(tauaer_lw (nofcols,nlay,nbndlw)) ;tauaer_lw =0.0
   allocate(ssaaer_lw (nofcols,nlay,nbndlw)) ;ssaaer_lw =0.0
   allocate(asmaer_lw (nofcols,nlay,nbndlw)) ;asmaer_lw =0.0

   allocate(ecaer     (nofcols,nlay,naerec)) ;ecaer     =0.0

   allocate(swhr      (nofcols,nlay))        ;swhr      =0.0
   allocate(swhrc     (nofcols,nlay))        ;swhrc     =0.0
   allocate(hr        (nofcols,nlay))        ;hr        =0.0
   allocate(hrc       (nofcols,nlay))        ;hrc       =0.0

   !
   allocate(swuflx   (nofcols,nlay+2));  swuflx   =0.0 !srf set "nlay+2" for
   allocate(swdflx   (nofcols,nlay+2));  swdflx   =0.0 !    the extra layer
   allocate(swuflxc  (nofcols,nlay+2));  swuflxc  =0.0 !    at the top of
   allocate(swdflxc  (nofcols,nlay+2));  swdflxc  =0.0 !    the model.
   allocate(uflx     (nofcols,nlay+2));  uflx	  =0.0 !    Otherwise, it can
   allocate(dflx     (nofcols,nlay+2));  dflx	  =0.0 !    be "nlay+1".
   allocate(uflxc    (nofcols,nlay+2));  uflxc    =0.0
   allocate(dflxc    (nofcols,nlay+2));  dflxc    =0.0
   allocate(duflx_dt (nofcols,nlay+2));  duflx_dt =0.0
   allocate(duflxc_dt(nofcols,nlay+2));  duflxc_dt=0.0

   allocate(p     (nofcols,mzp)) ;  p     =0.0
   allocate(t     (nofcols,mzp)) ;  t     =0.0
   allocate(q     (nofcols,mzp)) ;  q     =0.0
   allocate(o3l   (nofcols,mzp)) ;  o3l   =0.0
   allocate(fice  (nofcols,nlay));  fice  =0.0
   allocate(lmixr (nofcols,nlay));  lmixr =0.0
   allocate(psurf (nofcols))     ;  psurf =0.0
   allocate(colrad(nofcols))     ;  colrad=0.0
   allocate(tsfc  (nofcols))     ;  tsfc  =0.0
   allocate(asdir (nofcols))     ;  asdir =0.0
   allocate(aldir (nofcols))     ;  aldir =0.0
   allocate(asdif (nofcols))     ;  asdif =0.0
   allocate(aldif (nofcols))     ;  aldif =0.0
   allocate(coszen(nofcols))     ;  coszen=0.0

   
   tauc_sw      =0.0
   ssac_sw      =0.0
   asmc_sw      =0.0
   fsfc_sw      =0.0
   tauc_lw      =0.0
   ssac_lw      =0.0
   asmc_lw      =0.0
   fsfc_lw      =0.0
   taucloud_lw  =0.0
   taucloud_sw	=0.0
   taucldorig	=0.0
   ssacloud_sw	=0.0
   ssacloud_lw	=0.0
   asmcloud_sw	=0.0
   asmcloud_lw	=0.0

   !-initialization of rrtm memory
   if(firsttime) call initrrtm()
   !
   !- day of year
   dyofyr=julday(imonth1, idate1, iyear1)

   !prsnz  = (pi01dn(mzp,1)/cp)**cpor*p00
   aot_rrtm_sw=0.0
   aot_rrtm_lw=0.0

   if(firsttime .or. chemistry < 0) then

      tauaer_sw=0.0_rb
      ssaaer_sw=0.0_rb
      asmaer_sw=0.0_rb

      tauaer_lw=0.0_rb
      ssaaer_lw=0.0_rb
      asmaer_lw=0.0_rb
      
   end if

   !-
   !- column counter
   noc=0
   !- loop over the node domain
   do j=ja,jz  
      do i=ia,iz  
         noc=noc+1
         ipos(noc)=i
         jpos(noc)=j

         !- 2d input data
         coszen(noc)= radiate_g(ngrid)%cosz(i,j)
         tsfc  (noc)=(radiate_g(ngrid)%rlongup(i,j)/stefan)** 0.25e0_rb

         if(.not. firsttime) then
            asdir(noc)=radiate_g(ngrid)%albedt(i,j)
            aldir(noc)=radiate_g(ngrid)%albedt(i,j) !infrared
            asdif(noc)=radiate_g(ngrid)%albedt(i,j)
            aldif(noc)=radiate_g(ngrid)%albedt(i,j) !infrared
            emis(noc,:)=1.-radiate_g(ngrid)%albedt(i,j)
         end if
         !
         !- surface pressure
         psurf(noc)=(basic_g(ngrid)%pi0(1,i,j) + basic_g(ngrid)%pp(1,i,j)+ &
                     basic_g(ngrid)%pi0(2,i,j) + basic_g(ngrid)%pp(2,i,j))*cpi*0.50
         psurf(noc)=(p00*psurf(noc)**cpor)*adjust
         !land/ocean mask
         do np=1,npatch  
            if(leaf_g(ngrid)%leaf_class(i,j,np)==0) then
               imask(noc)=0
               cycle
            end if !ocean
            if(leaf_g(ngrid)%leaf_class(i,j,np)==2) then
               imask(noc)=13
               cycle
            end if!land ice
         end do   
         colrad(noc)=pi*(1.00e0_rb-(grid_g(ngrid)%glat(i,j)+90.00e0_rb)/180.00e0_rb)

         !- 3d input data
         do k=2,mzp   
            picpi = (basic_g(ngrid)%pi0(k,i,j) + basic_g(ngrid)%pp(k,i,j)) * cpi
            p(noc,k-1) = (p00 * picpi ** cpor)*adjust
            t(noc,k-1) =basic_g(ngrid)%theta(k,i,j) * picpi
            q(noc,k-1) =basic_g(ngrid)%rv   (k,i,j)
            !
            cldf(noc,k-1)=cloud_fraction(k,i,j) !cloud_fraction
            clwp(noc,k-1)=lwl(k,i,j)            !total cloud liquid water (kg/kg for carma and g/m2 for rrtm)
            ciwp(noc,k-1)=iwl(k,i,j)            !total cloud ice water (kg/kg for carma and g/m2 for rrtm)

            relq(noc,k-1)=micro_g(ngrid)%rel(k,i,j) ! effective radius liquid

            if(iceflgsw == 1) then
	    
	      reic(noc,k-1)=min(130.0,max(13.0,micro_g(ngrid)%rei(k,i,j))) !effective radius ice
            
	    elseif(iceflgsw == 2) then 
	    
	      reic(noc,k-1)=min(131.0,max(5.0,micro_g(ngrid)%rei(k,i,j))) !effective radius ice
	    
	    elseif(iceflgsw == 3) then 
	    
	      reic(noc,k-1)=min(140.0,max(5.0,0.9021*micro_g(ngrid)%rei(k,i,j)-7.0365)) ! generalized effective radius ice
            
	    endif

         end do  
	 q(noc,mzp) = q(noc,mzp-1)
	 p(noc,mzp)=p(noc,mzp-1)-2.*(p(noc,mzp-1)-p(noc,mzp-2))
         t(noc,mzp)=t(noc,mzp-1)+(t(noc,mzp-1)-t(noc,mzp-2))/&
	                         (p(noc,mzp-1)-p(noc,mzp-2))*&
				 (p(noc,mzp  )-p(noc,mzp-1))

	 !- these are values in the middle of the atmospheric layers
	 do k=1,nlay  
            play(noc,k)=p(noc,k)
            tlay(noc,k)=t(noc,k)
          end do   

	 !- these are values in the interfaces of the atmos layers
	 do k=2,nlay  
            plev(noc,k)=(p(noc,k)+p(noc,k-1))*0.5
            tlev(noc,k)=(t(noc,k)+t(noc,k-1))*0.5
         end do   
         !-special treatment for level 1 (interface layer)
         plev(noc,1)=psurf(noc)
         tlev(noc,1)=0.5* ( basic_g(ngrid)%theta(1,i,j)*&
	                   (basic_g(ngrid)%pi0(1,i,j) + basic_g(ngrid)%pp(1,i,j)) * cpi + &

			    basic_g(ngrid)%theta(2,i,j)*&
	                   (basic_g(ngrid)%pi0(2,i,j) + basic_g(ngrid)%pp(2,i,j)) * cpi   &
			  )

         !-special treatment for level mzp (interface layer)
	 plev(noc,mzp)=plev(noc,mzp-1)-2.*(plev(noc,mzp-1)-play(noc,nlay))

         tlev(noc,mzp)=tlev(noc,mzp-1)+(tlev(noc,mzp-1)-tlev(noc,mzp-2))/&
	                               (plev(noc,mzp-1)-plev(noc,mzp-2))*&
				       (plev(noc,mzp  )-plev(noc,mzp-1))




         if(.not. (firsttime .or. chemistry < 0)) then
            do k=2,mzp  
!KML 
               do np=1,nbndsw  
                  tauaer_sw(noc,k-1,np)=optProp(1,np)%tauaer(k,i,j)
                  aot_rrtm_sw(i,j,np) = aot_rrtm_sw(i,j,np)+tauaer_sw(noc,k-1,np)
!LFR *********************************************************************
                  !Total of AOT integrated in column for output
                  if(k==mzp) carma(ngrid)%aot(i,j,np)=aot_rrtm_sw(i,j,np)
		  !if(k==mzp) then
		  ! if(carma(ngrid)%aot(i,j,np)>0.1) print*,"aot=",mynum,np, carma(ngrid)%aot(i,j,np);call flush(6)
		  !endif
!LFR *********************************************************************
                  ssaaer_sw(noc,k-1,np)=optProp(ngrid,np)%ssa(k,i,j)
                  asmaer_sw(noc,k-1,np)=optProp(ngrid,np)%asp(k,i,j)
               end do  

               do np=1,nbndlw  
                  tauaer_lw(noc,k-1,np)=optProp(1,np+nbndsw)%tauaer(k,i,j)
                  aot_rrtm_lw(i,j,np) = aot_rrtm_lw(i,j,np)+tauaer_lw(noc,k-1,np)
!LFR *********************************************************************
                  !Total of AOT integrated in column for output
                  if(k==mzp) carma(ngrid)%aot(i,j,np+nbndsw)=aot_rrtm_lw(i,j,np)
		  !if(k==mzp) then
		  ! if(carma(ngrid)%aot(i,j,np)>0.1) print*,"aot=",mynum,np, carma(ngrid)%aot(i,j,np);call flush(6)
		  !endif
!LFR *********************************************************************
                  ssaaer_lw(noc,k-1,np)=optProp(ngrid,np+nbndsw)%ssa(k,i,j)
                  asmaer_lw(noc,k-1,np)=optProp(ngrid,np+nbndsw)%asp(k,i,j)
               end do  
!KML 
            end do  
          end if
      end do
   end do


   ecaer =0.0_rb

   do noc=1,nofcols  
      !calculating sig and delta sig
      do i=2,mzp  
!srf         sig(i)=plev(2,mzp-i+1)/ plev(2,1) !plev(2,mzp)
         sig(i)=plev(noc,i)/ plev(noc,1) !plev(2,mzp)
      end do   
      sig(1)=1.
      sig(mzp+1)=0.

      do i=1,mzp  
         if (sig(i)<0.1) exit
      end do      
      nls=mzp-i+1
      do k=1,mzp   
         delsig(k)=sig(k)-sig(k+1)
      end do       
      delsig(mzp+1)=delsig(mzp)
      do k=1,mzp   
          sigmid(k)=(sig(k)+sig(k+1))/2
          !print*,"sig=",sigmid(k),k,plev(2,mzp),sig(k)
      end do  
      sigmid(mzp+1)=sigmid(mzp)
      !print*,"sig=",sig
      call initgetoz_rrtm(365.2500_r8,mzp,sig)

    !kml calling getoz for both cases
    doy=real(dyofyr,r8)
    call getoz_rrtm (mzp,sigmid,colrad(noc),doy,o3l(noc,:))

    !call dumpBugRad(sig,mzp,ipos(noc),jpos(noc))
  enddo  

    if(o3ispresent) then

      !srf o3l(:,1)=o3l(:,2)    !kml eliminate o3 surface layer of the climatological profile

      do noc=1,nofcols  
         do k=2,mzp  
          !orig  o3vmr(noc,k-1)=co3fromtuv(k,ipos(noc),jpos(noc))*ro3ar/ &
          !orig                 basic_g(ngrid)%dn0(k,ipos(noc),jpos(noc))
          !orig  o3vmr(noc,k-1)=chem1_g(o3pos,ngrid)%sc_p(k,ipos(noc),jpos(noc))*ro3ar

          !kml sum o3 climatological and model profiles
              o3vmr(noc,k-1)=chem1_g(o3pos,ngrid)%sc_p(k,ipos(noc),jpos(noc))*       &
	                     ro3ar*basic_g(ngrid)%dn0(k,ipos(noc),jpos(noc))*factorcor +&
                             o3l(noc,flip(k-1))*ro3ar
         end do  
      end do  
   else
      do noc=1,nofcols  
         do k=1,mzp-1   
            o3vmr(noc,k)=o3l(noc,flip(k))*ro3ar !!! *1.e+4  << check!
            !print*,"o3vmr=",noc,k,real(o3vmr(noc,k),4)
         end do   
      end do   
   end if
!srf   o3vmr=o3vmr*factorcor

   if(co2ispresent) then
      do noc=1,nofcols
         do k=2,mzp
            co2vmr(noc,k-1)=chem1_g(co2pos,ngrid)%sc_p(k,ipos(noc),jpos(noc))*rco2ar*factorcor
         end do
      end do
   else
      co2vmr=379.e-6*rco2ar !fixed value from ipcc
   end if
!srf   co2vmr=co2vmr*factorcor

   if(ch4ispresent) then
      do noc=1,nofcols
         do k=2,mzp
            ch4vmr(noc,k-1)=chem1_g(ch4pos,ngrid)%sc_p(k,ipos(noc),jpos(noc))*rch4ar*factorcor
         end do
      end do
   else
      ch4vmr=1774.e-9*rch4ar !fixed value from ipcc
   end if
!srf   ch4vmr=ch4vmr*factorcor

   if(no2ispresent) then
      do noc=1,nofcols
         do k=2,mzp
            n2ovmr(noc,k-1)=chem1_g(no2pos,ngrid)%sc_p(k,ipos(noc),jpos(noc))*rno2ar*factorcor
         end do
      end do
   else
      n2ovmr=319.e-9 *rno2ar!fixed value from ipcc
   end if
!srf   n2ovmr=n2ovmr*factorcor

   if(o2ispresent) then
      do noc=1,nofcols
         do k=2,mzp
            o2vmr(noc,k-1)=chem1_g(o2pos,ngrid)%sc_p(k,ipos(noc),jpos(noc))*ro2ar*factorcor
         end do
      end do
   else
      o2vmr=0.209488*ro2ar !fixed value from ipcc
   end if
!srf   o2vmr=o2vmr*factorcor

   do noc=1,nofcols
      do k=2,mzp
         h2ovmr(noc,k-1)=q(noc,k)*rh2oar
      end do
   end do

   cfc11vmr=0.251e-9 !srf* factorcor
   cfc12vmr=0.538e-9 !srf* factorcor
   cfc22vmr=0.169e-9 !srf* factorcor
   ccl4vmr =0.093e-9 !srf* factorcor

!-----tmp
!print*,"1",maxval(h2ovmr),maxloc(h2ovmr(1,:))
!print*,"2",maxval(o3vmr ),maxloc(o3vmr (1,:))
!print*,"3",maxval(co2vmr)
!print*,"4",maxval(ch4vmr )
!print*,"5",maxval(n2ovmr )
!print*,"6",maxval(o2vmr   )
!print*,"7",maxval(cfc11vmr)
!print*,"8",maxval(cfc12vmr)
!print*,"9",maxval(cfc22vmr)
!print*,"10",maxval(ccl4vmr)
!call flush(6)
!-----tmp


   !- longwave section

     do noc=1,nofcols
        call cldprop(nlay, inflglw, iceflglw, liqflglw, cldf(noc,:), tauc_lw, &
                         ciwp(noc,:), clwp(noc,:), &
                         reic(noc,:), relq(noc,:), ncbandslw, taucloud_lw)
        do np=1,nbndlw
           taucl_lw(np,noc,:)=taucloud_lw(:,np)

	end do
     end do

      changeseed = permuteseed + icount 

      call mcica_subcol_lw(iplon, nofcols, nlay, icld, changeseed, irng, play, & 
                       cldf, ciwp, clwp, reic, relq, taucl_lw,  &
                       alpha_lw, &
                       cldfmcl_lw, ciwpmcl_lw, clwpmcl_lw, reicmcl_lw, relqmcl_lw, &
                       taucmcl_lw)


!-srf print for debug
!   print*,"max/min  taucmclc_lw",maxval(taucmcl_lw),minval(taucmcl_lw)
!   write(4 ,*) ,real(  nofcols    ,4)  ;call flush(4)
!   write(8 ,*) ,real(  nlay	  ,4) ;call flush(8)
!   write(7 ,*) ,real(  icld	  ,4) ;call flush(7)
!   write(9 ,*) ,real(  idrv	  ,4) ;call flush(9)
!   write(10,*) ,real(  play	  ,4) ;call flush(10)
!   write(11,*) ,real(  plev	  ,4) ;call flush(11)
!   write(12,*) ,real(  tlay	  ,4) ;call flush(12)
!   write(13,*) ,real(  tlev	  ,4) ;call flush(13)
!   write(14,*) ,real(  tsfc	  ,4) ;call flush(14)
!   write(15,*) ,real(  h2ovmr	  ,4) ;call flush(15)
!   write(16,*) ,real(  o3vmr	   ,4) ;call flush(16)
!   write(20,*) ,real(  co2vmr     ,4) ;call flush(20)
!   write(21,*) ,real(  ch4vmr     ,4) ;call flush(21)
!   write(62,*) ,real(  n2ovmr     ,4) ;call flush(62)
!   write(23,*) ,real(  o2vmr	   ,4) ;call flush(23)
!   write(24,*) ,real(  cfc11vmr    ,4) ;call flush(24)
!   write(25,*) ,real(  cfc12vmr   ,4) ;call flush(25)
!   write(26,*) ,real(  cfc22vmr   ,4) ;call flush(26)
!   write(27,*) ,real(  ccl4vmr    ,4) ;call flush(27)
!   write(28,*) ,real(  emis	   ,4)  ;call flush(28)
!   write(31,*) ,real(  inflglw    ,4)  ;call flush(31)
!   write(34,*) ,real(  iceflglw   ,4)  ;call flush(34)
!   write(40,*) ,real(  liqflglw   ,4)  ;call flush(40)
!   write(41,*) ,real(  cldfmcl_lw ,4)  ;call flush(41)
!   write(42,*) ,real(  taucmcl_lw ,4)  ;call flush(42)
!   write(43,*) ,real(  ciwpmcl_lw ,4)  ;call flush(43)
!   write(44,*) ,real(  clwpmcl_lw ,4)  ;call flush(44)
!   write(45,*) ,real(  reicmcl    ,4)  ;call flush(45)
!   write(46,*) ,real(  relqmcl    ,4)  ;call flush(46)
!   write(47,*) ,real(  tauaer_lw  ,4)  ;call flush(47)
!   write(48,*) ,real(  uflx	    ,4) ;call flush(48)
!   write(49,*) ,real(  dflx	    ,4) ;call flush(49)
!   write(50,*) ,real(  hr  ,4)         ;call flush(50)
!   write(51,*) ,real(  uflxc	   ,4) ;call flush(51)
!   write(52,*) ,real(  dflxc	   ,4) ;call flush(52)
!   write(54,*) ,real(  hrc  ,4)         ;call flush(54)
!   write(55,*) ,real(  duflx_dt   ,4)  ;call flush(55)
!   write(56,*) ,real(  duflxc_dt  ,4)  ;call flush(56)
!call flush(6)
!----------------------------------

   !- longwave calculations
   call rrtmg_lw (&
             nofcols   , &!01
             nlay      , &!02
             icld      , &!03
             idrv      , &!04
             play      , &!05
             plev      , &!06
             tlay      , &!07
             tlev      , &!08
             tsfc      , &!09
             h2ovmr    , &!10
             o3vmr     , &!11
             co2vmr    , &!12
             ch4vmr    , &!13
             n2ovmr    , &!14
             o2vmr     , &!15
             cfc11vmr  , &!16
             cfc12vmr  , &!17
             cfc22vmr  , &!18
             ccl4vmr   , &!19
             emis      , &!20
             inflglw   , &!21
             iceflglw  , &!22
             liqflglw  , &!23
             cldfmcl_lw, &!24
             taucmcl_lw, &!25
             ciwpmcl_lw, &!26
             clwpmcl_lw, &!27
             reicmcl_lw, &!28
             relqmcl_lw, &!29
             tauaer_lw , &!30
!             ssaaer_lw , &!31
!             asmaer_lw , &!32
             uflx      , &!33
             dflx      , &!34
             hr        , &!35
             uflxc     , &!36
             dflxc     , &!37
             hrc       , &!38
             duflx_dt  , &!39
             duflxc_dt   &!40
             )

   !--- shortwave section
   do noc=1,nofcols
          call cldprop_sw(nlay, inflgsw, iceflgsw, liqflgsw, cldf(noc,:) , &
                            tauc_sw, ssac_sw, asmc_sw, fsfc_sw, ciwp(noc,:), &
                            clwp(noc,:),reic(noc,:), relq(noc,:), &
                            taucldorig, taucloud_sw, ssacloud_sw, asmcloud_sw)
      do np=1,nbndsw

         taucl_sw(np,noc,:)= taucloud_sw(:,np+15)
         ssacl_sw(np,noc,:)= ssacloud_sw(:,np+15)
         asmcl_sw(np,noc,:)= asmcloud_sw(:,np+15)

      end do
   end do

      changeseed = changeseed + ngpt

   call mcica_subcol_sw(iplon, nofcols, nlay, icld, changeseed, irng, play, & 
                       cldf, ciwp, clwp, reic, relq, taucl_sw, ssacl_sw, asmcl_sw, fsfcl_sw, &
                       alpha_sw, &
                       cldfmcl_sw, ciwpmcl_sw, clwpmcl_sw, reicmcl_sw, relqmcl_sw, &
                       taucmcl_sw, ssacmcl_sw, asmcmcl_sw, fsfcmcl_sw)

      changeseed = changeseed - ngpt

   !examples of dump 3d
   !call dumpvar(real(reicmcl),'rei','-001',nofcols,nlay)
   !call dumpvar(real(relqmcl),'rel','-001',nofcols,nlay)

   !- shortwave calculations
   call rrtmg_sw &
            (nofcols , &!01                         
             nlay    , &!02
             icld    , &!03
             iaer    , &! LFR for the 5.0 version
             play    , &!04
             plev    , &!05
             tlay    , &!06
             tlev    , &!07
             tsfc    , &!08 *
             h2ovmr  , &!09
             o3vmr   , &!10
             co2vmr  , &!11
             ch4vmr  , &!12
             n2ovmr  , &!13
             o2vmr   , &!14 *
             asdir   , &!15
             asdif   , &!16
             aldir   , &!17
             aldif   , &!18 *
             coszen  , &!19
             adjes   , &!20
             dyofyr  , &!21
             scon    , &!22 *
             isolvar , & !LFR for the 5.0 version
             inflgsw , &!23
             iceflgsw, &!24
             liqflgsw, &!25
             cldfmcl_sw, &!26 *
             taucmcl_sw, &!27
             ssacmcl_sw, &!28
             asmcmcl_sw, &!29
             fsfcmcl_sw, &!30 *
             ciwpmcl_sw, &!31
             clwpmcl_sw, &!32
             reicmcl_sw, &!33
             relqmcl_sw, &!34 *
             tauaer_sw, &!35
             ssaaer_sw, &!36
             asmaer_sw, &!37
             ecaer   , &!38 *
             swuflx  , &!40
             swdflx  , &!41
             swhr    , &!42 *
             swuflxc , &!43
             swdflxc , &!44
             swhrc     &!45 *
             !bndsolvar, & !LFR for the 5.0 version
             !indsolvar, & !LFR for the 5.0 version
             !solcycfrac & !LFR for the 5.0 version
             )

   !-set initialization variable to false, so
   !- next timestep it will no be performed anymore.
   firsttime=.false.

   !- sending out theta tendency and surface radiative fluxes
   noc=0
   do j=ja,jz
      do i=ia,iz
         noc=noc+1
         radiate_g(ngrid)%rlong (i,j)=   dflx(noc,1)
         radiate_g(ngrid)%rshort(i,j)= swdflx(noc,1)
         if(swdflx(noc,1)<.5) radiate_g(ngrid)%rshort(i,j)=0.0
         do k=2,mzp-1
           !- radiative tendency of temperature (Kelvin/sec)
           radiate_g(ngrid)%fthrd(k,i,j)=(swhr(noc,k-1)+hr(noc,k-1))/86400.
         enddo
         radiate_g(ngrid)%fthrd(1  ,i,j)= radiate_g(ngrid)%fthrd(2    ,i,j)
         radiate_g(ngrid)%fthrd(mzp,i,j)= radiate_g(ngrid)%fthrd(mzp-1,i,j)
      end do
   end do

   deallocate(ipos  )
   deallocate(jpos  )
   deallocate(imask )
   deallocate(play  )
   deallocate(plev  )
   deallocate(tlay  )
   deallocate(tlev  )
   deallocate(h2ovmr)
   deallocate(o3vmr )
   deallocate(co2vmr)
   deallocate(ch4vmr)
   deallocate(n2ovmr)
   deallocate(o2vmr )
   deallocate(cfc11vmr)
   deallocate(cfc12vmr)
   deallocate(cfc22vmr)
   deallocate(ccl4vmr )
   deallocate(emis )

   deallocate(cldfmcl_sw)
   deallocate(taucmcl_sw)
   deallocate(ssacmcl_sw)
   deallocate(asmcmcl_sw)
   deallocate(fsfcmcl_sw)
   deallocate(ciwpmcl_sw)
   deallocate(clwpmcl_sw)

   deallocate(reicmcl_sw)
   deallocate(relqmcl_sw)

   deallocate(cldfmcl_lw)
   deallocate(taucmcl_lw)
   deallocate(ssacmcl_lw)
   deallocate(asmcmcl_lw)
   deallocate(fsfcmcl_lw)
   deallocate(ciwpmcl_lw)
   deallocate(clwpmcl_lw)

   deallocate(reicmcl_lw)
   deallocate(relqmcl_lw)

   deallocate(cldf)
   deallocate(ciwp)
   deallocate(clwp)
   deallocate(reic)
   deallocate(relq)

   deallocate(taucl_sw)
   deallocate(ssacl_sw)
   deallocate(asmcl_sw)
   deallocate(fsfcl_sw)
 
   deallocate(taucl_lw)
   deallocate(ssacl_lw)
   deallocate(asmcl_lw)
   deallocate(fsfcl_lw)
   
   deallocate(tauaer_sw) 
   deallocate(ssaaer_sw) 
   deallocate(asmaer_sw) 
   deallocate(tauaer_lw) 
   deallocate(ssaaer_lw) 
   deallocate(asmaer_lw) 
   
   deallocate(ecaer)

   deallocate(swhr      )
   deallocate(swhrc     )
   deallocate(hr        )
   deallocate(hrc       )

   deallocate(swuflx   )
   deallocate(swdflx   )
   deallocate(swuflxc  )
   deallocate(swdflxc  )
   deallocate(uflx     )
   deallocate(dflx     )
   deallocate(uflxc    )
   deallocate(dflxc    )
   deallocate(duflx_dt )
   deallocate(duflxc_dt)

   deallocate(p     )
   deallocate(t     )
   deallocate(q     )
   deallocate(o3l   )
   deallocate(fice  )
   deallocate(lmixr )
   deallocate(psurf )
   deallocate(colrad)
   deallocate(tsfc  )
   deallocate(asdir )
   deallocate(aldir )
   deallocate(asdif )
   deallocate(aldif )
   deallocate(coszen)
end subroutine radrrtmdrv

! ****************************************************************************

subroutine cloud_prop_rrtm(m1,m2,m3,ia, iz, ja, jz &
                     , cloud_fraction  &
                     , rain            &
                     , lwl             &
                     , iwl             &
                       )

  implicit none
  integer, intent(in) :: m1,m2,m3,ia,iz,ja,jz

  real, intent(out), dimension(m1,m2,m3) :: cloud_fraction !cloud_fraction
  real, intent(out), dimension(m2,m3   ) :: rain !total rain water
  real, intent(out), dimension(m1,m2,m3) :: lwl !total cloud liquid water (kg/kg for carma and g/m2 for rrtm)
  real, intent(out), dimension(m1,m2,m3) :: iwl !total cloud ice water (kg/kg for carma and g/m2 for rrtm)

  !-- local variables
  integer, parameter :: r8=4
  integer :: i,j,k,k700,kp1,kdthdp

  real, parameter :: rhminl = .90              ! minimum rh for low stable clouds
  real, parameter :: rhminh = .80              ! minimum rh for high stable clouds
  real, parameter :: sh1 = 0.07 ,sh2= 500.0   ! parameters for shallow convection cloud fraction
  real, parameter :: dp1 = 0.14 ,dp2= 500.0   ! parameters for deep convection cloud fraction
  real, parameter :: premit= 750.e2              ! top pressure bound for mid level cloud
  real, parameter :: pnot = 1.e5       ! reference pressure
  real, parameter :: lapse = 6.5e-3    ! u.s. standard atmsophere lapse rate
  real, parameter :: premib = 750.e2   ! bottom pressure bound of middle cloud
  real, parameter :: pretop = 1.0e2    ! pressure bounding high cloud
  real, parameter :: abeta = 0.07
  real, parameter :: bbeta = -0.14
  real, parameter :: pi = 3.14159265358979323846
  real, parameter :: bx = 100.* (3./(4.*pi))**(1./3.)
  real, parameter :: r13 = 1./3.
  real, parameter :: r13bbeta = 1./3. - 0.14

  real, dimension(m1) :: press,rh,cldst,concld,rhu00,rpdeli
  real, dimension(m1,m2,m3) :: dummy_vec
  real, dimension(m2,m3   ) :: upmf,upmfsh
  real, dimension(m1,m2,m3) :: zup,zupshallow,clwup,clwupsh
  real :: ocnfrac,picpi,temp,dztri,strat,rhpert,shallowcu,deepcu,upmfshx,upmfx&
             ,rhwght,rhdif,rhlim,ps,thetas,dthdpmn,dthdp,dummy,bb,factor
  logical cldbnd          ! region below high cloud boundary
  real, external :: rs
  character(len=3) :: cm1
  !- tuning parameters to include direct coupling between cupar and radiation
  real, parameter :: tun_rad_shall=0.02, tun_rad_deep =0.02
  integer,parameter :: coupl_rad_cupar=1 ! 0 -no, 1-yes
  integer :: n
! real, dimension(m1,m2,m3) :: clcn

  ! set defaults for rhu00
  rhu00(:) = 2.0
  ! define rh perturbation in order to estimate rhdfda
  rhpert = 0.01

  !--- initialization of cuparm parameters
  !--- convective cloud fraction
  !clcn = 0.0
  !
  if(nnqparm(ngrid) == 5 .or. nnqparm(ngrid) == 6 .or. nnqparm(ngrid) == 8) then
    upmf      (     1:m2,1:m3)= xmb4d(     1:m2,1:m3,1,ngrid) !- mass flux deep    convection
    upmfsh    (     1:m2,1:m3)= xmb4d(     1:m2,1:m3,2,ngrid) !- mass flux shallow convection
    zup       (1:m1,1:m2,1:m3)= zup5d(1:m1,1:m2,1:m3,1,ngrid) !- normalized mass flux
    zupshallow(1:m1,1:m2,1:m3)= zup5d(1:m1,1:m2,1:m3,2,ngrid) !- normalized mass flux

    if(coupl_rad_cupar == 1 ) then
      clwup     (1:m1,1:m2,1:m3)= tun_rad_deep *clwup5d(1:m1,1:m2,1:m3,1,ngrid)
      clwupsh   (1:m1,1:m2,1:m3)= tun_rad_shall*clwup5d(1:m1,1:m2,1:m3,2,ngrid)
      !if (nnqparm(ngrid) == 8) & ! includes mid convection
      !     clwup(1:m1,1:m2,1:m3) = clwup(1:m1,1:m2,1:m3) + tun_rad_deep *clwup5d(1:m1,1:m2,1:m3,3,ngrid)
    endif
    !--- convective cloud fraction (to implement this scheme, CLCN needs to be transported as a scalar)
    !do j=ja,jz
    ! do i=ia,iz
    !   do k=2,m1
    !      dztri= grid_g(ngrid)%rtgt(i,j)/dzt(k)
    !
    !      factor=dtlt/(dztri * basic_g(ngrid)%dn0(k,i,j))
    !
    ! 	  do n=1,maxiens
    !        clcn(k,i,j)=clcn(k,i,j)+(1.0 - clcn(k,i,j))* factor *&
    !                 xmb4d(i,j,n,ngrid) * up_massdetr5d(k-1,i,j,n,ngrid)
    !enddo;enddo;enddo;enddo
    !print*,"CL=",maxval(clcn),minval(clcn);call flush(6)
  else
    upmf      (     1:m2,1:m3)=0.0
    upmfsh    (     1:m2,1:m3)=0.0
    zup       (1:m1,1:m2,1:m3)=0.0
    zupshallow(1:m1,1:m2,1:m3)=0.0
    if(coupl_rad_cupar == 1 ) then
      clwup     (1:m1,1:m2,1:m3)=0.0
      clwupsh   (1:m1,1:m2,1:m3)=0.0
    endif
  endif


  ! print*,"=============radiation 0 ==================="
  ! print*,"max/min theta",maxval(basic_g(ngrid)%theta),minval(basic_g(ngrid)%theta)
  ! print*,"max/min    pi",maxval(basic_g(ngrid)%pi0),minval(basic_g(ngrid)%pp)
  ! print*,"=============radiation 0 ==================="
  ! call flush(6)

  do j=ja,jz
   do i=ia,iz

     ! evaluate potential temperature and relative humidity
     do k=1,m1
           picpi = (basic_g(ngrid)%pi0(k,i,j) + basic_g(ngrid)%pp(k,i,j)) * cpi
           press(k) = p00 * picpi ** cpor
           !LFR test
!           if(0.01*press(k) .gt.  700.) then
!            write(80,fmt='(3(I3.3,1X),3(E18.3,1X))') i,j,k,press(k),basic_g(ngrid)%pi0(k,i,j),basic_g(ngrid)%pp(k,i,j)
!           endif
           temp = basic_g(ngrid)%theta(k,i,j) * picpi
           rh(k) =min(1.,max(0.05,basic_g(ngrid)%rv(k,i,j)/rs(press(k),temp)))
           !------
           !
           cloud_fraction(k,i,j)  = 0.
           cldst (k)      = 0.
           concld(k)      = 0.
    enddo
    ps    =0.5*(press(1)    +press(2))
    thetas=0.5*(basic_g(ngrid)%theta(1,i,j)+basic_g(ngrid)%theta(2,i,j))
    !print*,"press=",ps,thetas;call flush(6)
    !
    ocnfrac=0.
    if(leaf_g(ngrid)%patch_area(i,j,1)>0.99) ocnfrac=1. !flag < 1 para land
                                          !flag  =1 para water

    do k=1,m1-1
          rpdeli(k) = 1./(press(k+1) - press(k))
    end do
    !
    ! cloud mass flux in si units of kg/m2/s; should produce typical numbers of 20%
    ! shallow and deep convective cloudiness are evaluated separately (since processes
    ! are evaluated separately) and summed
    !
    do k=2,m1
          shallowcu=0.0; upmfshx=0.0
          if(upmfsh(i,j)>1.e-8) then
	     upmfshx=zupshallow(k-1,i,j) * upmfsh(i,j)
             !print*,"upmfshx",k,upmfshx,upmfsh(i,j),zupshallow(k+1,i,j); call flush(6)
	     !-orig shallowcu = max(0.0,min(sh1*log(1.0+sh2*cmfmc2(i,k+1)),0.30))
             shallowcu = max(0.0,min(sh1*log(1.0+sh2*upmfshx),0.30))
          endif
	  deepcu = 0.
	  if(upmf(i,j) >1.e-8) then
	     upmfx     = zup(k-1,i,j) * upmf(i,j)

             !-orig deepcu = max(0.0,min(dp1*log(1.0+dp2*(cmfmc(i,k+1)-cmfmc2(i,k+1))),0.60))
             ! check possibility of (upmfx-upmfshx) <0. => log (x<0)
             deepcu = max(0.0,min(dp1*log(1.0+dp2*max(0.,(upmfx-upmfshx))),0.60))
          endif
          concld(k) = min(shallowcu + deepcu,0.80)
          rh(k) = (rh(k) - concld(k))/(1.0 - concld(k))
          !if(concld(k) > 0.) then
          !   print*,"1=",k, concld(k),rh(k) ,upmfshx,upmfx ;call flush(6)
          !endif
    end do

    do k=2,m1
        kp1 = min(k + 1,m1)
        !
        cldbnd = press(k).ge.pretop
          if ( press(k).ge.premib ) then
             !==============================================================
             ! this is the low cloud (below premib) block
             !==============================================================
             ! enhance low cloud activation over land with no snow cover
             if ( ocnfrac < 0.999 ) then !.and. (snowh(i) <= 0.000001)) then
                rhlim = rhminl - 0.10
             else
                rhlim = rhminl
             endif
             !
             rhdif = (rh(k) - rhlim)/(1.0-rhlim)
             cloud_fraction(k,i,j) = min(0.999,(max(rhdif,0.0))**2)
          else if ( press(k).lt.premit ) then
             !==============================================================
             ! this is the high cloud (above premit) block
             !==============================================================
             !
             rhlim = rhminh
             !
             rhdif = (rh(k) - rhlim)/(1.0-rhlim)
             cloud_fraction(k,i,j) = min(0.999,(max(rhdif,0.0))**2)
          else
             !==============================================================
             ! this is the middle cloud block
             !==============================================================
             !
             !       linear rh threshold transition between thresholds for low & high cloud
             !
             rhwght = (premib-(max(press(k),premit)))/(premib-premit)

             if ( ocnfrac < 0.999 ) then !if (land(i) .and. (snowh(i) <= 0.000001)) then
                rhlim = rhminh*rhwght + (rhminl - 0.10)*(1.0-rhwght)
             else
                rhlim = rhminh*rhwght + rhminl*(1.0-rhwght)
             endif
             rhdif = (rh(k) - rhlim)/(1.0-rhlim)
             cloud_fraction(k,i,j) = min(0.999,(max(rhdif,0.0))**2)
          end if
          !print*,"2=",k,i,j, cloud(k,i,j) ;call flush(6)
          !      ! save rhlim to rhu00, it handles well by itself for low/high cloud
          !      !
          rhu00(k)=rhlim
          !==================================================================================
      end do


      !--- stratus
      ! find most stable level below 750 mb for evaluating stratus regimes
      ! nothing triggers unless a stability greater than this minimum threshold is found
      dthdpmn = -0.125
      kdthdp  = 0

      do k=2,m1
          if (press(k) >= premib .and. ocnfrac.gt. 0.01) then
             ! i think this is done so that dtheta/dp is in units of dg/mb (jjh)
             dthdp = 100.0*(basic_g(ngrid)%theta(k,i,j) - basic_g(ngrid)%theta(k-1,i,j))*rpdeli(k-1)
             if (dthdp < dthdpmn) then
                dthdpmn = dthdp
                kdthdp  = k     ! index of interface of max inversion
             endif
          endif
      enddo

      ! also check between the bottom layer and the surface
      ! only perform this check if the criteria were not met above

      if ( kdthdp .eq. 0 .and. ocnfrac.gt.0.01) then
                 dthdp = 100.0 * (thetas - basic_g(ngrid)%theta(m1,i,j)) / (ps-press(m1))
                 if (dthdp < dthdpmn) then
                    dthdpmn = dthdp
                    kdthdp  = m1     ! index of interface of max inversion
                 endif
      endif
      do k=2,m1-1
           k700=-99
           !print*,"k70",k,press(k)
           if(0.01*press(k) .le.  700.) then
             k700=k
             exit
           endif
      enddo
      !print*,        k700
      if( k700 > m1 .or.   k700 < 1) then
         write (*,fmt='(a)') "wrong k700 at cloud_prop routine, see press:"
         write (cm1,fmt='(i3.3)') m1
         write (*,fmt='('//cm1//'(e9.3,1x))') (press(k),k=1,m1)
         write (*,fmt='(a,i3.3,a,i3.3,a,i3.3)') 'k700= ',k700, ' - for column ',i,' , ',j
         stop
      end if
      if (kdthdp /= 0) then
          k = kdthdp
          kp1 = min(k+1,m1)
          ! note: strat will be zero unless ocnfrac > 0.01
          strat = min(1.,max(0., ocnfrac * ((basic_g(ngrid)%theta(k700,i,j)-thetas)*0.057-0.5573)))
          !
          ! assign the stratus to the layer just below max inversion
          ! the relative humidity changes so rapidly across the inversion
          ! that it is not safe to just look immediately below the inversion
          ! so limit the stratus cloud by rh in both layers below the inversion
          !
           cldst(k) = min(strat,max(rh(k),rh(kp1)))
       endif

       !
       ! aggregate cloud contributions (cldst should be zero everywhere except at level kdthdp(i))
       !
       do k=1,m1
                 !
                 !       which is greater; standard layered cloud amount or stratocumulus diagnosis
                 !
                 cloud_fraction(k,i,j) = max(cloud_fraction(k,i,j),cldst(k))
                 !
                 !       add in the contributions of convective cloud (determined separately and accounted
                 !       for by modifications to the large-scale relative humidity.
                 !
                 cloud_fraction(k,i,j) = min(cloud_fraction(k,i,j)+concld(k), 1.0)
             !print*,"cloudfraction=",i,j,k,cloud(k,i,j),ocnfrac;call flush(6)
       end do
       !print*,"cloudfraction=",i,j,maxval(cloud(:,i,j)),ocnfrac;call flush(6)


  enddo;enddo


!-------start calculation of cloud ....


       dummy_vec = 0.0
       ! if level == 1 do nothing
       if (level==2) then
          lwl(1:m1,ia:iz,ja:jz) = micro_g(ngrid)%rcp(1:m1,ia:iz,ja:jz)

	  if (nnqparm(ngrid)/=0) then
             rain(ia:iz,ja:jz)= cuparm_g(ngrid)%conprr(ia:iz,ja:jz)* 3600.
          endif

       elseif (level>=3) then
        if (nnqparm(ngrid)/=0) then
          rain(ia:iz,ja:jz) = cuparm_g(ngrid)%conprr(ia:iz,ja:jz) + &
                          micro_g(ngrid)%pcpg(ia:iz,ja:jz)
        else
                rain(ia:iz,ja:jz) = micro_g(ngrid)%pcpg(ia:iz,ja:jz)
        endif
        rain(ia:iz,ja:jz) = rain(ia:iz,ja:jz)*3600.

        if (icloud>0) lwl(1:m1,ia:iz,ja:jz) = lwl(1:m1,ia:iz,ja:jz) + micro_g(ngrid)%rcp(1:m1,ia:iz,ja:jz)
        if (igraup>0) then
          if(mcphys_type <= 1) then
            do k=1,m1
              do i=ia,iz
                do j=ja,jz
                   call qtc(micro_g(ngrid)%q6(k,i,j), dummy,dummy_vec(k,i,j))
                enddo
              enddo
            enddo
          elseif(mcphys_type == 2 .or. mcphys_type == 3 .or. mcphys_type == 4 .or. &
                 mcphys_type == 5 .or. mcphys_type == 6 .or. mcphys_type == 7) then  !srf -gthompson/gfdl microphysics - graupel only in ice phase
            dummy_vec=0.0
          endif
          lwl(1:m1,ia:iz,ja:jz) = dummy_vec(1:m1,ia:iz,ja:jz)*micro_g(ngrid)%rgp(1:m1,ia:iz,ja:jz) &
                                + lwl(1:m1,ia:iz,ja:jz) !kg/kg
          dummy_vec(:,:,:)      = 1. - dummy_vec(:,:,:)
          iwl(1:m1,ia:iz,ja:jz) = dummy_vec(1:m1,ia:iz,ja:jz)*micro_g(ngrid)%rgp(1:m1,ia:iz,ja:jz) &
                                + iwl(1:m1,ia:iz,ja:jz) !kg/kg
        endif
        if (ihail>0) then
          dummy_vec = 0.
          if(mcphys_type <= 1) then
           do k=1,m1
             do i=ia,iz
                do j=ja,jz
                   call qtc(micro_g(ngrid)%q7(k,i,j), dummy, dummy_vec(k,i,j))
                enddo
             enddo
           enddo
          endif
          lwl(1: m1,ia:iz,ja:jz) = dummy_vec(1:m1,ia:iz,ja:jz)*micro_g(ngrid)%rhp(1:m1,ia:iz,ja:jz)&
                                + lwl(1:m1,ia:iz,ja:jz)

          dummy_vec(:,:,:) = 1.0 - dummy_vec(:,:,:)

          iwl(1:m1,ia:iz,ja:jz) = dummy_vec(1:m1,ia:iz,ja:jz)*micro_g(ngrid)%rhp(1:m1,ia:iz,ja:jz) &
                                + iwl(1:m1,ia:iz,ja:jz)   !kg/kg
        endif
        if (iaggr>0) &
            iwl(1:m1,ia:iz,ja:jz) = iwl(1:m1,ia:iz,ja:jz) + micro_g(ngrid)%rap(1:m1,ia:iz,ja:jz)   !kg/kg

        if (isnow>0) &
            iwl(1:m1,ia:iz,ja:jz) = iwl(1:m1,ia:iz,ja:jz) + micro_g(ngrid)%rsp(1:m1,ia:iz,ja:jz)   !kg/kg

        if (ipris>0) &
            iwl(1:m1,ia:iz,ja:jz) = iwl(1:m1,ia:iz,ja:jz) + micro_g(ngrid)%rpp(1:m1,ia:iz,ja:jz)  !kg/kg
       endif
       !- making direct couplig between liq/ice water from cupar to radiation
       if(coupl_rad_cupar == 1 ) then
         !print*,"curad=",maxval(lwl),maxval(iwl),(1./tun_rad_deep)*maxval(clwup)&
         !              ,(1./tun_rad_shall)*maxval(clwupsh);call flush(6)
         do j = ja,jz
          do i = ia,iz
           do k = 1,m1
            temp = basic_g(ngrid)%theta(k,i,j)*(basic_g(ngrid)%pp(k,i,j)+basic_g(ngrid)%pi0(k,i,j))*cpi ! air temp (kelvin)

            if(temp .gt. 253.)then
                 lwl(k,i,j)=lwl(k,i,j)+ clwup(k,i,j)+  clwupsh(k,i,j)
            else
                 iwl(k,i,j)=iwl(k,i,j)+ clwup(k,i,j)+  clwupsh(k,i,j)
            endif
         enddo; enddo; enddo
       endif

       do j=ja,jz
         do i=ia,iz
           do k=1,m1
                lwl(k,i,j) = max(0.,lwl(k,i,j))  !kg/kg
                iwl(k,i,j) = max(0.,iwl(k,i,j))  !kg/kg
           enddo
         enddo
       enddo
       do j=ja,jz
        do i=ia,iz
          rain(i,j) = max(0.,rain(i,j))
        enddo
       enddo

       !- rrtmg radiation  -  continue calculation of cloud effective radius

       if (ilwrtyp==6 .or. iswrtyp==6) then

          if(mcphys_type == 2 .or. mcphys_type == 4 .or.  mcphys_type == 5 .or. & 
             mcphys_type == 6 .or. mcphys_type == 7 ) then
          !- rei and rel are calculated by gt microphysics (no-aer option)

           do j=ja,jz
            do i=ia,iz
             do k=1,m1

                picpi = (basic_g(ngrid)%pi0(k,i,j) + basic_g(ngrid)%pp(k,i,j)) * cpi
                press(k) = p00 * picpi ** cpor
                temp = basic_g(ngrid)%theta(k,i,j) * picpi

                lwl(k,i,j) = basic_g(ngrid)%dn0(k,i,j)*lwl(k,i,j)* 1.e+3  !g/m3
                iwl(k,i,j) = basic_g(ngrid)%dn0(k,i,j)*iwl(k,i,j)* 1.e+3  !g/m3
                dztri= grid_g(ngrid)%rtgt(i,j)/dzt(k)

		          if(iwl(k,i,j)<1.0e-6 .or. temp>273.0) iwl(k,i,j)=0.0
                lwl(k,i,j) = lwl(k,i,j)* dztri  !g/m2
                iwl(k,i,j) = iwl(k,i,j)* dztri  !g/m2
             enddo
            enddo
           enddo

          else
          !- for all other options, rei and rel are calculated in the section below
          !
           do j=ja,jz
            do i=ia,iz
             do k=1,m1

              !- liquid cloud effective radius -----
              !- [liu&daum, 2000 and 2005. liu et al 2008]
              !- cloud drop number concentration
                dummy_vec(k,i,j) = micro_g(ngrid)%ccp(k,i,j) * basic_g(ngrid)%dn0(k,i,j) * 1.e-6  ! #/cm3
                dummy_vec(k,i,j) = max(150.,dummy_vec(k,i,j))

              !
              !- liquid cloud effective radius ----- [liu&daum, 2000 and 2005]

                lwl(k,i,j) = basic_g(ngrid)%dn0(k,i,j)*lwl(k,i,j)* 1.e+3  !g/m3
                iwl(k,i,j) = basic_g(ngrid)%dn0(k,i,j)*iwl(k,i,j)* 1.e+3  !g/m3

                !- u[lwl] = g/m3 , u[dummy_vec]= #/cm^3
                !- the cte bx with the units above, provides rel in 10^-6 meter
                !  micro_g(ngrid)%rel(k,i,j)= bx *  ( lwl(k,i,j) /dummy_vec(k,i,j))**r13 &
                !- the factor below is adimensional e provides correction
                !- for dispersion of the cloud spectrum:
                !- prefactor of liu et al (2008) - lwl must be in g/cm^3
                !          *abeta*(1.e-6*lwl(k,i,j) /dummy_vec(k,i,j))**bbeta
                !- expression to avoid nan when lwl = 0. in prefactor
                micro_g(ngrid)%rel(k,i,j)= bx *  ( lwl(k,i,j) /dummy_vec(k,i,j))**r13bbeta &
                            *abeta*6.92 !6.92=(1.e-6)**bbeta

                ! rel is limited between 2.5 and 60 micrometers as
                ! required by rrtm parameterization
                micro_g(ngrid)%rel(k,i,j) = max(2.5, min( 60.0, micro_g(ngrid)%rel(k,i,j) ) )

              !------ice cloud effective radius ----- [klaus wyser, 1998]

                picpi = (basic_g(ngrid)%pi0(k,i,j) + basic_g(ngrid)%pp(k,i,j)) * cpi
                press(k) = p00 * picpi ** cpor
                temp = basic_g(ngrid)%theta(k,i,j) * picpi
                dztri= grid_g(ngrid)%rtgt(i,j)/dzt(k)
                if(iwl(k,i,j)<1.0e-6 .or. temp>273.0) then
                  micro_g(ngrid)%rei(k,i,j)=5.0
                  iwl(k,i,j)=0.0
                else
                  bb = -2. + log10(iwl(k,i,j)/50.)*(1.e-3*(273.15-max(210.15,temp))**1.5)
                  micro_g(ngrid)%rei(k,i,j) =377.4 + 203.3 * bb+ 37.91 * bb **2 + 2.3696 * bb **3
                  !print*,"bb=",temp,micro_g(ngrid)%rei(k,i,j),bb,iwl(k,i,j);call flush(6)
                endif
!srf -tmp ---------------------
!if(lwl(k,i,j) > 1.e-6) then
!  write(4,1000) k,lwl(k,i,j),micro_g(ngrid)%rel(k,i,j), abeta*(1.e-6*lwl(k,i,j) /dummy_vec(k,i,j))**bbeta &
!                      ,micro_g(ngrid)%ccp(k,i,j)* basic_g(ngrid)%dn0(k,i,j) * 1.e-6, temp!&
!     1000 format(1x,i4,5e13.4)
!     call flush(4)
!endif
!if(iwl(k,i,j) > 1.e-6) then
!      write(14,1001) k,iwl(k,i,j),rei(k,i,j), temp!&
!     1001 format(1x,i4,3e13.4)
!     call flush(14)
!endif
!srf -tmp ---------------------

                lwl(k,i,j) = lwl(k,i,j)* dztri  !g/m2
                iwl(k,i,j) = iwl(k,i,j)* dztri  !g/m2
               enddo
             enddo
          enddo
         endif

       endif
 return
 print*,"=============mcphys-radiation coupling ==================="
 print*,"max-min  cl_frac:  ",maxval(cloud_fraction(2:m1,ia:iz,ja:jz)),minval(cloud_fraction(2:m1,ia:iz,ja:jz)) !cloud_fraction
 print*,"max-min     rain:  ",maxval(rain(ia:iz,ja:jz)),minval(rain(ia:iz,ja:jz)) !total rain water
 print*,"max-min lwl: ",maxval(lwl(2:m1,ia:iz,ja:jz)),minval(lwl(2:m1,ia:iz,ja:jz)) !total cloud liquid water (kg/kg for carma and g/m2 for rrtm)
 print*,"max-min iwl: ",maxval(iwl(2:m1,ia:iz,ja:jz)),minval(iwl(2:m1,ia:iz,ja:jz))     !total cloud ice water (kg/kg for carma and g/m2 for rrtm)
 print*,"max-min micro_g(ngrid)%rel: ",maxval(micro_g(ngrid)%rel(2:m1,ia:iz,ja:jz)),minval(micro_g(ngrid)%rel(2:m1,ia:iz,ja:jz))  !total cloud liquid water
 print*,"max-min rei: ",maxval(micro_g(ngrid)%rei(2:m1,ia:iz,ja:jz)),minval(micro_g(ngrid)%rei(2:m1,ia:iz,ja:jz))  !total cloud ice water

 !print*,"=============================================="
 call flush(6)

 end subroutine cloud_prop_rrtm

! ****************************************************************************


 subroutine InitGetoz_rrtm(yrl,kmax,sl)
    implicit none
    real(kind=r8),    intent(in   ) :: yrl
    integer, intent(in   ) :: kmax
    real(kind=r8),    intent(in   ) :: sl (kmax)
    integer                :: l, ll


    if(.not. allocated(ozone)) allocate(ozone(nlm_getoz,nl,ns))

    !
    !     four season climatological ozone data in nmc sigma layers
    !
    !     for seasonal variation
    !     season=1 - winter          season=2 - spring
    !     season=3 - summer          season=4 - fall
    !     unit of ozone mixing ratio is in ( 10**-4 g/g ).  the data is
    !     in 18 sigma layers from top to bottom.  for every layer, there
    !     are 37 latitudes at 5 degree interval from north pole to south
    !     pole.
    !     mrf86 18 layers
    !
    !
    !     1. winter
    !
    !     wint1(18,6)
    !
    ozone(1:18, 1:6, 1) = RESHAPE( (/ &
         .068467e0_r8,.052815e0_r8,.035175e0_r8,.022334e0_r8,.013676e0_r8,.007363e0_r8, &
         .003633e0_r8,.001582e0_r8,.001111e0_r8,.000713e0_r8,.000517e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .069523e0_r8,.052249e0_r8,.034255e0_r8,.021379e0_r8,.012306e0_r8,.006727e0_r8, &
         .003415e0_r8,.001578e0_r8,.001072e0_r8,.000681e0_r8,.000517e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .070579e0_r8,.051684e0_r8,.033335e0_r8,.020423e0_r8,.010935e0_r8,.006091e0_r8, &
         .003197e0_r8,.001573e0_r8,.001034e0_r8,.000650e0_r8,.000517e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .074885e0_r8,.049987e0_r8,.030140e0_r8,.017894e0_r8,.009881e0_r8,.005543e0_r8, &
         .002907e0_r8,.001379e0_r8,.000961e0_r8,.000644e0_r8,.000512e0_r8,.000463e0_r8, &
         .000451e0_r8,.000408e0_r8,.000385e0_r8,.000361e0_r8,.000351e0_r8,.000349e0_r8, &
         .079190e0_r8,.048290e0_r8,.026945e0_r8,.015366e0_r8,.008826e0_r8,.004995e0_r8, &
         .002616e0_r8,.001184e0_r8,.000887e0_r8,.000637e0_r8,.000508e0_r8,.000486e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .082443e0_r8,.047591e0_r8,.025358e0_r8,.014294e0_r8,.008233e0_r8,.004664e0_r8, &
         .002430e0_r8,.001068e0_r8,.000851e0_r8,.000644e0_r8,.000508e0_r8,.000474e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8/), &
         (/18,6/))
    !
    !     wint2(18,6)
    !
    ozone(1:18, 7:12, 1) = RESHAPE( (/ &
         .085695e0_r8,.046892e0_r8,.023772e0_r8,.013223e0_r8,.007640e0_r8,.004333e0_r8, &
         .002244e0_r8,.000951e0_r8,.000815e0_r8,.000650e0_r8,.000508e0_r8,.000463e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .089618e0_r8,.042869e0_r8,.019963e0_r8,.010502e0_r8,.005966e0_r8,.003525e0_r8, &
         .001936e0_r8,.000906e0_r8,.000769e0_r8,.000625e0_r8,.000508e0_r8,.000452e0_r8, &
         .000451e0_r8,.000408e0_r8,.000385e0_r8,.000361e0_r8,.000351e0_r8,.000349e0_r8, &
         .093540e0_r8,.038846e0_r8,.016155e0_r8,.007781e0_r8,.004292e0_r8,.002716e0_r8, &
         .001628e0_r8,.000862e0_r8,.000724e0_r8,.000600e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .097097e0_r8,.034916e0_r8,.012983e0_r8,.006240e0_r8,.003666e0_r8,.002259e0_r8, &
         .001336e0_r8,.000730e0_r8,.000629e0_r8,.000549e0_r8,.000499e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .100654e0_r8,.030986e0_r8,.009812e0_r8,.004698e0_r8,.003041e0_r8,.001803e0_r8, &
         .001044e0_r8,.000599e0_r8,.000533e0_r8,.000499e0_r8,.000491e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .101724e0_r8,.026500e0_r8,.007228e0_r8,.003391e0_r8,.002058e0_r8,.001285e0_r8, &
         .000811e0_r8,.000531e0_r8,.000478e0_r8,.000449e0_r8,.000440e0_r8,.000421e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    !
    !     wint3(18,6)
    !
    ozone(1:18, 13:18, 1) = RESHAPE( (/ &
         .102794e0_r8,.022015e0_r8,.004645e0_r8,.002084e0_r8,.001076e0_r8,.000767e0_r8, &
         .000577e0_r8,.000463e0_r8,.000423e0_r8,.000399e0_r8,.000389e0_r8,.000401e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .103456e0_r8,.018235e0_r8,.003195e0_r8,.001379e0_r8,.000771e0_r8,.000585e0_r8, &
         .000474e0_r8,.000411e0_r8,.000380e0_r8,.000362e0_r8,.000343e0_r8,.000348e0_r8, &
         .000346e0_r8,.000328e0_r8,.000317e0_r8,.000305e0_r8,.000302e0_r8,.000302e0_r8, &
         .104118e0_r8,.014455e0_r8,.001745e0_r8,.000674e0_r8,.000467e0_r8,.000403e0_r8, &
         .000370e0_r8,.000359e0_r8,.000337e0_r8,.000325e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8, &
         .104106e0_r8,.012997e0_r8,.001479e0_r8,.000639e0_r8,.000468e0_r8,.000422e0_r8, &
         .000392e0_r8,.000372e0_r8,.000342e0_r8,.000325e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8, &
         .104093e0_r8,.011539e0_r8,.001213e0_r8,.000604e0_r8,.000468e0_r8,.000442e0_r8, &
         .000414e0_r8,.000385e0_r8,.000347e0_r8,.000325e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8, &
         .104087e0_r8,.010726e0_r8,.000971e0_r8,.000538e0_r8,.000440e0_r8,.000434e0_r8, &
         .000418e0_r8,.000397e0_r8,.000375e0_r8,.000343e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    !
    !     wint4(18,6)
    !
    ozone(1:18, 19:24, 1) = RESHAPE( (/ &
         .102665e0_r8,.010977e0_r8,.001237e0_r8,.000590e0_r8,.000498e0_r8,.000479e0_r8, &
         .000458e0_r8,.000436e0_r8,.000421e0_r8,.000387e0_r8,.000326e0_r8,.000298e0_r8, &
         .000246e0_r8,.000227e0_r8,.000211e0_r8,.000200e0_r8,.000194e0_r8,.000186e0_r8, &
         .100892e0_r8,.012873e0_r8,.001886e0_r8,.000785e0_r8,.000643e0_r8,.000568e0_r8, &
         .000519e0_r8,.000487e0_r8,.000471e0_r8,.000437e0_r8,.000368e0_r8,.000305e0_r8, &
         .000201e0_r8,.000151e0_r8,.000117e0_r8,.000098e0_r8,.000090e0_r8,.000093e0_r8, &
         .100534e0_r8,.013704e0_r8,.002028e0_r8,.000861e0_r8,.000701e0_r8,.000604e0_r8, &
         .000546e0_r8,.000513e0_r8,.000504e0_r8,.000462e0_r8,.000381e0_r8,.000307e0_r8, &
         .000201e0_r8,.000151e0_r8,.000117e0_r8,.000098e0_r8,.000090e0_r8,.000093e0_r8, &
         .100218e0_r8,.015035e0_r8,.002537e0_r8,.001037e0_r8,.000790e0_r8,.000726e0_r8, &
         .000673e0_r8,.000628e0_r8,.000579e0_r8,.000512e0_r8,.000440e0_r8,.000374e0_r8, &
         .000307e0_r8,.000253e0_r8,.000227e0_r8,.000208e0_r8,.000194e0_r8,.000186e0_r8, &
         .099903e0_r8,.016365e0_r8,.003045e0_r8,.001214e0_r8,.000879e0_r8,.000848e0_r8, &
         .000801e0_r8,.000744e0_r8,.000654e0_r8,.000562e0_r8,.000499e0_r8,.000441e0_r8, &
         .000410e0_r8,.000358e0_r8,.000342e0_r8,.000322e0_r8,.000302e0_r8,.000302e0_r8, &
         .099547e0_r8,.017725e0_r8,.003693e0_r8,.001578e0_r8,.001125e0_r8,.000985e0_r8, &
         .000879e0_r8,.000795e0_r8,.000712e0_r8,.000643e0_r8,.000584e0_r8,.000521e0_r8, &
         .000482e0_r8,.000384e0_r8,.000351e0_r8,.000322e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    !
    !     wint5(18,6)
    !
    ozone(1:18, 25:30, 1) = RESHAPE( (/ &
         .099191e0_r8,.019085e0_r8,.004340e0_r8,.001943e0_r8,.001371e0_r8,.001122e0_r8, &
         .000957e0_r8,.000847e0_r8,.000770e0_r8,.000724e0_r8,.000669e0_r8,.000601e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .098107e0_r8,.020617e0_r8,.004758e0_r8,.002137e0_r8,.001516e0_r8,.001211e0_r8, &
         .000999e0_r8,.000848e0_r8,.000778e0_r8,.000730e0_r8,.000677e0_r8,.000603e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .097023e0_r8,.022148e0_r8,.005177e0_r8,.002332e0_r8,.001660e0_r8,.001300e0_r8, &
         .001041e0_r8,.000849e0_r8,.000786e0_r8,.000737e0_r8,.000686e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .093464e0_r8,.026177e0_r8,.008525e0_r8,.003892e0_r8,.002452e0_r8,.001609e0_r8, &
         .001116e0_r8,.000851e0_r8,.000809e0_r8,.000762e0_r8,.000690e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .089906e0_r8,.030206e0_r8,.011873e0_r8,.005453e0_r8,.003244e0_r8,.001918e0_r8, &
         .001192e0_r8,.000852e0_r8,.000832e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .080939e0_r8,.032414e0_r8,.014163e0_r8,.007241e0_r8,.004328e0_r8,.002522e0_r8, &
         .001481e0_r8,.000934e0_r8,.000861e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8/), &
         (/18,6/))
    !
    !     wint6(18,6)
    !
    ozone(1:18, 31:36, 1) = RESHAPE( (/ &
         .071972e0_r8,.034622e0_r8,.016453e0_r8,.009029e0_r8,.005413e0_r8,.003127e0_r8, &
         .001770e0_r8,.001015e0_r8,.000890e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .069820e0_r8,.035028e0_r8,.016929e0_r8,.009389e0_r8,.005645e0_r8,.003260e0_r8, &
         .001843e0_r8,.001055e0_r8,.000905e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .067669e0_r8,.035434e0_r8,.017406e0_r8,.009749e0_r8,.005876e0_r8,.003393e0_r8, &
         .001916e0_r8,.001094e0_r8,.000920e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .065518e0_r8,.035975e0_r8,.017854e0_r8,.010100e0_r8,.006534e0_r8,.003985e0_r8, &
         .002321e0_r8,.001240e0_r8,.000966e0_r8,.000774e0_r8,.000640e0_r8,.000548e0_r8, &
         .000479e0_r8,.000384e0_r8,.000346e0_r8,.000316e0_r8,.000302e0_r8,.000302e0_r8, &
         .063367e0_r8,.036516e0_r8,.018302e0_r8,.010452e0_r8,.007192e0_r8,.004577e0_r8, &
         .002727e0_r8,.001387e0_r8,.001012e0_r8,.000762e0_r8,.000585e0_r8,.000490e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .061216e0_r8,.037359e0_r8,.019151e0_r8,.010633e0_r8,.006845e0_r8,.004382e0_r8, &
         .002691e0_r8,.001511e0_r8,.001061e0_r8,.000749e0_r8,.000568e0_r8,.000465e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    !
    !     wint7(18)
    !
    ozone(1:18, 37, 1) = (/ &
         .059066e0_r8,.038201e0_r8,.019999e0_r8,.010813e0_r8,.006498e0_r8,.004188e0_r8, &
         .002656e0_r8,.001636e0_r8,.001110e0_r8,.000737e0_r8,.000551e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/)
    !
    !     2. spring
    !
    ozone(1:18, 1:6, 2) = RESHAPE( (/ &
         .074229e0_r8,.050084e0_r8,.030930e0_r8,.018676e0_r8,.011965e0_r8,.008165e0_r8, &
         .005428e0_r8,.003399e0_r8,.002098e0_r8,.001138e0_r8,.000780e0_r8,.000632e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .074927e0_r8,.049459e0_r8,.029215e0_r8,.018025e0_r8,.011754e0_r8,.007786e0_r8, &
         .004972e0_r8,.002926e0_r8,.001817e0_r8,.001025e0_r8,.000758e0_r8,.000632e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .075625e0_r8,.048835e0_r8,.027500e0_r8,.017375e0_r8,.011544e0_r8,.007407e0_r8, &
         .004516e0_r8,.002453e0_r8,.001536e0_r8,.000912e0_r8,.000737e0_r8,.000632e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .077409e0_r8,.048159e0_r8,.026661e0_r8,.016596e0_r8,.010962e0_r8,.006972e0_r8, &
         .004160e0_r8,.002132e0_r8,.001391e0_r8,.000868e0_r8,.000686e0_r8,.000601e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .079194e0_r8,.047483e0_r8,.025822e0_r8,.015818e0_r8,.010380e0_r8,.006537e0_r8, &
         .003804e0_r8,.001811e0_r8,.001245e0_r8,.000825e0_r8,.000635e0_r8,.000570e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .084591e0_r8,.046553e0_r8,.025037e0_r8,.015156e0_r8,.009841e0_r8,.006124e0_r8, &
         .003534e0_r8,.001693e0_r8,.001170e0_r8,.000793e0_r8,.000631e0_r8,.000537e0_r8, &
         .000551e0_r8,.000509e0_r8,.000486e0_r8,.000516e0_r8,.000548e0_r8,.000446e0_r8/), &
         (/18,6/))
    ozone(1:18, 7:12, 2) = RESHAPE( (/ &
         .089988e0_r8,.045622e0_r8,.024253e0_r8,.014495e0_r8,.009303e0_r8,.005711e0_r8, &
         .003264e0_r8,.001574e0_r8,.001096e0_r8,.000762e0_r8,.000627e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .092863e0_r8,.042419e0_r8,.020704e0_r8,.012034e0_r8,.007417e0_r8,.004504e0_r8, &
         .002590e0_r8,.001334e0_r8,.000977e0_r8,.000731e0_r8,.000622e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .095737e0_r8,.039215e0_r8,.017155e0_r8,.009572e0_r8,.005532e0_r8,.003296e0_r8, &
         .001916e0_r8,.001094e0_r8,.000858e0_r8,.000699e0_r8,.000618e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .097501e0_r8,.035382e0_r8,.014856e0_r8,.008207e0_r8,.004619e0_r8,.002720e0_r8, &
         .001610e0_r8,.001012e0_r8,.000829e0_r8,.000687e0_r8,.000610e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .099264e0_r8,.031548e0_r8,.012557e0_r8,.006841e0_r8,.003705e0_r8,.002144e0_r8, &
         .001304e0_r8,.000930e0_r8,.000799e0_r8,.000675e0_r8,.000601e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .101718e0_r8,.026523e0_r8,.008473e0_r8,.004382e0_r8,.002392e0_r8,.001505e0_r8, &
         .001036e0_r8,.000836e0_r8,.000727e0_r8,.000618e0_r8,.000550e0_r8,.000494e0_r8, &
         .000501e0_r8,.000479e0_r8,.000473e0_r8,.000509e0_r8,.000541e0_r8,.000445e0_r8/), &
         (/18,6/))
    ozone(1:18, 13:18, 2) = RESHAPE( (/ &
         .104172e0_r8,.021499e0_r8,.004389e0_r8,.001922e0_r8,.001078e0_r8,.000865e0_r8, &
         .000767e0_r8,.000743e0_r8,.000654e0_r8,.000562e0_r8,.000499e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .104145e0_r8,.018082e0_r8,.003274e0_r8,.001493e0_r8,.000919e0_r8,.000762e0_r8, &
         .000678e0_r8,.000641e0_r8,.000584e0_r8,.000531e0_r8,.000495e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .104118e0_r8,.014665e0_r8,.002159e0_r8,.001063e0_r8,.000759e0_r8,.000659e0_r8, &
         .000589e0_r8,.000539e0_r8,.000514e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .107719e0_r8,.013052e0_r8,.001822e0_r8,.000953e0_r8,.000701e0_r8,.000604e0_r8, &
         .000551e0_r8,.000525e0_r8,.000509e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .111320e0_r8,.011439e0_r8,.001485e0_r8,.000843e0_r8,.000642e0_r8,.000549e0_r8, &
         .000512e0_r8,.000512e0_r8,.000504e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .112375e0_r8,.011255e0_r8,.001357e0_r8,.000744e0_r8,.000585e0_r8,.000533e0_r8, &
         .000512e0_r8,.000512e0_r8,.000504e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8/), &
         (/18,6/))
    ozone(1:18, 19:24, 2) = RESHAPE( (/ &
         .109850e0_r8,.010424e0_r8,.001079e0_r8,.000567e0_r8,.000498e0_r8,.000479e0_r8, &
         .000463e0_r8,.000448e0_r8,.000418e0_r8,.000399e0_r8,.000389e0_r8,.000367e0_r8, &
         .000351e0_r8,.000328e0_r8,.000320e0_r8,.000337e0_r8,.000355e0_r8,.000304e0_r8, &
         .107002e0_r8,.009961e0_r8,.001025e0_r8,.000533e0_r8,.000497e0_r8,.000460e0_r8, &
         .000422e0_r8,.000385e0_r8,.000332e0_r8,.000300e0_r8,.000288e0_r8,.000249e0_r8, &
         .000202e0_r8,.000158e0_r8,.000132e0_r8,.000114e0_r8,.000104e0_r8,.000093e0_r8, &
         .107735e0_r8,.010146e0_r8,.001120e0_r8,.000576e0_r8,.000526e0_r8,.000477e0_r8, &
         .000430e0_r8,.000385e0_r8,.000332e0_r8,.000300e0_r8,.000288e0_r8,.000249e0_r8, &
         .000202e0_r8,.000158e0_r8,.000132e0_r8,.000114e0_r8,.000104e0_r8,.000093e0_r8, &
         .107021e0_r8,.012233e0_r8,.001533e0_r8,.000643e0_r8,.000556e0_r8,.000505e0_r8, &
         .000471e0_r8,.000448e0_r8,.000403e0_r8,.000362e0_r8,.000355e0_r8,.000296e0_r8, &
         .000251e0_r8,.000207e0_r8,.000180e0_r8,.000161e0_r8,.000152e0_r8,.000140e0_r8, &
         .106308e0_r8,.014320e0_r8,.001946e0_r8,.000709e0_r8,.000585e0_r8,.000533e0_r8, &
         .000512e0_r8,.000512e0_r8,.000473e0_r8,.000425e0_r8,.000423e0_r8,.000342e0_r8, &
         .000301e0_r8,.000257e0_r8,.000232e0_r8,.000212e0_r8,.000205e0_r8,.000209e0_r8, &
         .100592e0_r8,.015718e0_r8,.002411e0_r8,.001007e0_r8,.000802e0_r8,.000642e0_r8, &
         .000559e0_r8,.000526e0_r8,.000501e0_r8,.000474e0_r8,.000470e0_r8,.000439e0_r8, &
         .000430e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    ozone(1:18, 25:30, 2) = RESHAPE( (/ &
         .094877e0_r8,.017116e0_r8,.002877e0_r8,.001305e0_r8,.001018e0_r8,.000751e0_r8, &
         .000606e0_r8,.000539e0_r8,.000529e0_r8,.000524e0_r8,.000516e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .094163e0_r8,.020198e0_r8,.004594e0_r8,.001772e0_r8,.001077e0_r8,.000806e0_r8, &
         .000649e0_r8,.000565e0_r8,.000547e0_r8,.000537e0_r8,.000521e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .093449e0_r8,.023279e0_r8,.006312e0_r8,.002240e0_r8,.001135e0_r8,.000862e0_r8, &
         .000692e0_r8,.000591e0_r8,.000564e0_r8,.000549e0_r8,.000525e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .089886e0_r8,.026029e0_r8,.008558e0_r8,.003312e0_r8,.001655e0_r8,.001124e0_r8, &
         .000807e0_r8,.000631e0_r8,.000602e0_r8,.000568e0_r8,.000525e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .086323e0_r8,.028778e0_r8,.010805e0_r8,.004383e0_r8,.002175e0_r8,.001386e0_r8, &
         .000923e0_r8,.000671e0_r8,.000640e0_r8,.000587e0_r8,.000525e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .082715e0_r8,.031096e0_r8,.013350e0_r8,.006131e0_r8,.003205e0_r8,.002043e0_r8, &
         .001304e0_r8,.000842e0_r8,.000734e0_r8,.000631e0_r8,.000555e0_r8,.000494e0_r8, &
         .000480e0_r8,.000408e0_r8,.000385e0_r8,.000361e0_r8,.000351e0_r8,.000349e0_r8/), &
         (/18,6/))
    ozone(1:18, 31:36, 2) = RESHAPE( (/ &
         .079108e0_r8,.033415e0_r8,.015895e0_r8,.007878e0_r8,.004234e0_r8,.002700e0_r8, &
         .001686e0_r8,.001014e0_r8,.000829e0_r8,.000675e0_r8,.000584e0_r8,.000454e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .074807e0_r8,.034651e0_r8,.017056e0_r8,.008574e0_r8,.004769e0_r8,.002986e0_r8, &
         .001827e0_r8,.001079e0_r8,.000853e0_r8,.000675e0_r8,.000584e0_r8,.000454e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .070506e0_r8,.035887e0_r8,.018218e0_r8,.009270e0_r8,.005304e0_r8,.003271e0_r8, &
         .001969e0_r8,.001145e0_r8,.000878e0_r8,.000675e0_r8,.000584e0_r8,.000454e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .067669e0_r8,.037799e0_r8,.019680e0_r8,.009612e0_r8,.005481e0_r8,.003476e0_r8, &
         .002093e0_r8,.001123e0_r8,.000837e0_r8,.000631e0_r8,.000546e0_r8,.000447e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .064832e0_r8,.039712e0_r8,.021142e0_r8,.009954e0_r8,.005658e0_r8,.003681e0_r8, &
         .002218e0_r8,.001100e0_r8,.000796e0_r8,.000587e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .063734e0_r8,.039842e0_r8,.022004e0_r8,.010859e0_r8,.005712e0_r8,.003589e0_r8, &
         .002155e0_r8,.001174e0_r8,.000856e0_r8,.000612e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    ozone(1:18, 37, 2) = (/ &
         .062636e0_r8,.039972e0_r8,.022867e0_r8,.011765e0_r8,.005766e0_r8,.003498e0_r8, &
         .002092e0_r8,.001248e0_r8,.000917e0_r8,.000637e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/)
    !
    !     3. summer
    !
    ozone(1:18, 1:6, 3) = RESHAPE( (/ &
         .059066e0_r8,.038201e0_r8,.019999e0_r8,.010813e0_r8,.006498e0_r8,.004188e0_r8, &
         .002656e0_r8,.001636e0_r8,.001110e0_r8,.000737e0_r8,.000551e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .061216e0_r8,.037359e0_r8,.019151e0_r8,.010633e0_r8,.006845e0_r8,.004382e0_r8, &
         .002691e0_r8,.001511e0_r8,.001061e0_r8,.000749e0_r8,.000568e0_r8,.000465e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .063367e0_r8,.036516e0_r8,.018302e0_r8,.010452e0_r8,.007192e0_r8,.004577e0_r8, &
         .002727e0_r8,.001387e0_r8,.001012e0_r8,.000762e0_r8,.000585e0_r8,.000490e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .065518e0_r8,.035975e0_r8,.017854e0_r8,.010100e0_r8,.006534e0_r8,.003985e0_r8, &
         .002321e0_r8,.001240e0_r8,.000966e0_r8,.000774e0_r8,.000640e0_r8,.000548e0_r8, &
         .000479e0_r8,.000384e0_r8,.000346e0_r8,.000316e0_r8,.000302e0_r8,.000302e0_r8, &
         .067669e0_r8,.035434e0_r8,.017406e0_r8,.009749e0_r8,.005876e0_r8,.003393e0_r8, &
         .001916e0_r8,.001094e0_r8,.000920e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .069820e0_r8,.035028e0_r8,.016929e0_r8,.009389e0_r8,.005645e0_r8,.003260e0_r8, &
         .001843e0_r8,.001055e0_r8,.000905e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8/), &
         (/18,6/))
    ozone(1:18, 7:12, 3) = RESHAPE( (/ &
         .071972e0_r8,.034622e0_r8,.016453e0_r8,.009029e0_r8,.005413e0_r8,.003127e0_r8, &
         .001770e0_r8,.001015e0_r8,.000890e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .080939e0_r8,.032414e0_r8,.014163e0_r8,.007241e0_r8,.004328e0_r8,.002522e0_r8, &
         .001481e0_r8,.000934e0_r8,.000861e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .089906e0_r8,.030206e0_r8,.011873e0_r8,.005453e0_r8,.003244e0_r8,.001918e0_r8, &
         .001192e0_r8,.000852e0_r8,.000832e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .093464e0_r8,.026177e0_r8,.008525e0_r8,.003892e0_r8,.002452e0_r8,.001609e0_r8, &
         .001116e0_r8,.000851e0_r8,.000809e0_r8,.000762e0_r8,.000690e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .097023e0_r8,.022148e0_r8,.005177e0_r8,.002332e0_r8,.001660e0_r8,.001300e0_r8, &
         .001041e0_r8,.000849e0_r8,.000786e0_r8,.000737e0_r8,.000686e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .098107e0_r8,.020617e0_r8,.004758e0_r8,.002137e0_r8,.001516e0_r8,.001211e0_r8, &
         .000999e0_r8,.000848e0_r8,.000778e0_r8,.000730e0_r8,.000677e0_r8,.000603e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8/), &
         (/18,6/))
    ozone(1:18, 13:18, 3) = RESHAPE( (/ &
         .099191e0_r8,.019085e0_r8,.004340e0_r8,.001943e0_r8,.001371e0_r8,.001122e0_r8, &
         .000957e0_r8,.000847e0_r8,.000770e0_r8,.000724e0_r8,.000669e0_r8,.000601e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .099547e0_r8,.017725e0_r8,.003693e0_r8,.001578e0_r8,.001125e0_r8,.000985e0_r8, &
         .000879e0_r8,.000795e0_r8,.000712e0_r8,.000643e0_r8,.000584e0_r8,.000521e0_r8, &
         .000482e0_r8,.000384e0_r8,.000351e0_r8,.000322e0_r8,.000302e0_r8,.000302e0_r8, &
         .099903e0_r8,.016365e0_r8,.003045e0_r8,.001214e0_r8,.000879e0_r8,.000848e0_r8, &
         .000801e0_r8,.000744e0_r8,.000654e0_r8,.000562e0_r8,.000499e0_r8,.000441e0_r8, &
         .000410e0_r8,.000358e0_r8,.000342e0_r8,.000322e0_r8,.000302e0_r8,.000302e0_r8, &
         .100218e0_r8,.015035e0_r8,.002537e0_r8,.001037e0_r8,.000790e0_r8,.000726e0_r8, &
         .000673e0_r8,.000628e0_r8,.000579e0_r8,.000512e0_r8,.000440e0_r8,.000374e0_r8, &
         .000307e0_r8,.000253e0_r8,.000227e0_r8,.000208e0_r8,.000194e0_r8,.000186e0_r8, &
         .100534e0_r8,.013704e0_r8,.002028e0_r8,.000861e0_r8,.000701e0_r8,.000604e0_r8, &
         .000546e0_r8,.000513e0_r8,.000504e0_r8,.000462e0_r8,.000381e0_r8,.000307e0_r8, &
         .000201e0_r8,.000151e0_r8,.000117e0_r8,.000098e0_r8,.000090e0_r8,.000093e0_r8, &
         .100892e0_r8,.012873e0_r8,.001886e0_r8,.000785e0_r8,.000643e0_r8,.000568e0_r8, &
         .000519e0_r8,.000487e0_r8,.000471e0_r8,.000437e0_r8,.000368e0_r8,.000305e0_r8, &
         .000201e0_r8,.000151e0_r8,.000117e0_r8,.000098e0_r8,.000090e0_r8,.000093e0_r8/), &
         (/18,6/))
    ozone(1:18, 19:24, 3) = RESHAPE( (/ &
         .102665e0_r8,.010977e0_r8,.001237e0_r8,.000590e0_r8,.000498e0_r8,.000479e0_r8, &
         .000458e0_r8,.000436e0_r8,.000421e0_r8,.000387e0_r8,.000326e0_r8,.000298e0_r8, &
         .000246e0_r8,.000227e0_r8,.000211e0_r8,.000200e0_r8,.000194e0_r8,.000186e0_r8, &
         .104087e0_r8,.010726e0_r8,.000971e0_r8,.000538e0_r8,.000440e0_r8,.000434e0_r8, &
         .000418e0_r8,.000397e0_r8,.000375e0_r8,.000343e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8, &
         .104093e0_r8,.011539e0_r8,.001213e0_r8,.000604e0_r8,.000468e0_r8,.000442e0_r8, &
         .000414e0_r8,.000385e0_r8,.000347e0_r8,.000325e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8, &
         .104106e0_r8,.012997e0_r8,.001479e0_r8,.000639e0_r8,.000468e0_r8,.000422e0_r8, &
         .000392e0_r8,.000372e0_r8,.000342e0_r8,.000325e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8, &
         .104118e0_r8,.014455e0_r8,.001745e0_r8,.000674e0_r8,.000467e0_r8,.000403e0_r8, &
         .000370e0_r8,.000359e0_r8,.000337e0_r8,.000325e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8, &
         .103456e0_r8,.018235e0_r8,.003195e0_r8,.001379e0_r8,.000771e0_r8,.000585e0_r8, &
         .000474e0_r8,.000411e0_r8,.000380e0_r8,.000362e0_r8,.000343e0_r8,.000348e0_r8, &
         .000346e0_r8,.000328e0_r8,.000317e0_r8,.000305e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    ozone(1:18, 25:30, 3) = RESHAPE( (/ &
         .102794e0_r8,.022015e0_r8,.004645e0_r8,.002084e0_r8,.001076e0_r8,.000767e0_r8, &
         .000577e0_r8,.000463e0_r8,.000423e0_r8,.000399e0_r8,.000389e0_r8,.000401e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .101724e0_r8,.026500e0_r8,.007228e0_r8,.003391e0_r8,.002058e0_r8,.001285e0_r8, &
         .000811e0_r8,.000531e0_r8,.000478e0_r8,.000449e0_r8,.000440e0_r8,.000421e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .100654e0_r8,.030986e0_r8,.009812e0_r8,.004698e0_r8,.003041e0_r8,.001803e0_r8, &
         .001044e0_r8,.000599e0_r8,.000533e0_r8,.000499e0_r8,.000491e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .097097e0_r8,.034916e0_r8,.012983e0_r8,.006240e0_r8,.003666e0_r8,.002259e0_r8, &
         .001336e0_r8,.000730e0_r8,.000629e0_r8,.000549e0_r8,.000499e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .093540e0_r8,.038846e0_r8,.016155e0_r8,.007781e0_r8,.004292e0_r8,.002716e0_r8, &
         .001628e0_r8,.000862e0_r8,.000724e0_r8,.000600e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .089618e0_r8,.042869e0_r8,.019963e0_r8,.010502e0_r8,.005966e0_r8,.003525e0_r8, &
         .001936e0_r8,.000906e0_r8,.000769e0_r8,.000625e0_r8,.000508e0_r8,.000452e0_r8, &
         .000451e0_r8,.000408e0_r8,.000385e0_r8,.000361e0_r8,.000351e0_r8,.000349e0_r8/), &
         (/18,6/))
    ozone(1:18, 31:36, 3) = RESHAPE( (/ &
         .085695e0_r8,.046892e0_r8,.023772e0_r8,.013223e0_r8,.007640e0_r8,.004333e0_r8, &
         .002244e0_r8,.000951e0_r8,.000815e0_r8,.000650e0_r8,.000508e0_r8,.000463e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .082443e0_r8,.047591e0_r8,.025358e0_r8,.014294e0_r8,.008233e0_r8,.004664e0_r8, &
         .002430e0_r8,.001068e0_r8,.000851e0_r8,.000644e0_r8,.000508e0_r8,.000474e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .079190e0_r8,.048290e0_r8,.026945e0_r8,.015366e0_r8,.008826e0_r8,.004995e0_r8, &
         .002616e0_r8,.001184e0_r8,.000887e0_r8,.000637e0_r8,.000508e0_r8,.000486e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .074885e0_r8,.049987e0_r8,.030140e0_r8,.017894e0_r8,.009881e0_r8,.005543e0_r8, &
         .002907e0_r8,.001379e0_r8,.000961e0_r8,.000644e0_r8,.000512e0_r8,.000463e0_r8, &
         .000451e0_r8,.000408e0_r8,.000385e0_r8,.000361e0_r8,.000351e0_r8,.000349e0_r8, &
         .070579e0_r8,.051684e0_r8,.033335e0_r8,.020423e0_r8,.010935e0_r8,.006091e0_r8, &
         .003197e0_r8,.001573e0_r8,.001034e0_r8,.000650e0_r8,.000517e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .069523e0_r8,.052249e0_r8,.034255e0_r8,.021379e0_r8,.012306e0_r8,.006727e0_r8, &
         .003415e0_r8,.001578e0_r8,.001072e0_r8,.000681e0_r8,.000517e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    ozone(1:18, 37, 3) = (/ &
         .068467e0_r8,.052815e0_r8,.035175e0_r8,.022334e0_r8,.013676e0_r8,.007363e0_r8, &
         .003633e0_r8,.001582e0_r8,.001111e0_r8,.000713e0_r8,.000517e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/)
    !
    !     4. fall
    !
    ozone(1:18, 1:6, 4) = RESHAPE( (/ &
         .062636e0_r8,.039972e0_r8,.022867e0_r8,.011765e0_r8,.005766e0_r8,.003498e0_r8, &
         .002092e0_r8,.001248e0_r8,.000917e0_r8,.000637e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .063734e0_r8,.039842e0_r8,.022004e0_r8,.010859e0_r8,.005712e0_r8,.003589e0_r8, &
         .002155e0_r8,.001174e0_r8,.000856e0_r8,.000612e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .064832e0_r8,.039712e0_r8,.021142e0_r8,.009954e0_r8,.005658e0_r8,.003681e0_r8, &
         .002218e0_r8,.001100e0_r8,.000796e0_r8,.000587e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .067669e0_r8,.037799e0_r8,.019680e0_r8,.009612e0_r8,.005481e0_r8,.003476e0_r8, &
         .002093e0_r8,.001123e0_r8,.000837e0_r8,.000631e0_r8,.000546e0_r8,.000447e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .070506e0_r8,.035887e0_r8,.018218e0_r8,.009270e0_r8,.005304e0_r8,.003271e0_r8, &
         .001969e0_r8,.001145e0_r8,.000878e0_r8,.000675e0_r8,.000584e0_r8,.000454e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .074807e0_r8,.034651e0_r8,.017056e0_r8,.008574e0_r8,.004769e0_r8,.002986e0_r8, &
         .001827e0_r8,.001079e0_r8,.000853e0_r8,.000675e0_r8,.000584e0_r8,.000454e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    ozone(1:18, 7:12, 4) = RESHAPE( (/ &
         .079108e0_r8,.033415e0_r8,.015895e0_r8,.007878e0_r8,.004234e0_r8,.002700e0_r8, &
         .001686e0_r8,.001014e0_r8,.000829e0_r8,.000675e0_r8,.000584e0_r8,.000454e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .082715e0_r8,.031096e0_r8,.013350e0_r8,.006131e0_r8,.003205e0_r8,.002043e0_r8, &
         .001304e0_r8,.000842e0_r8,.000734e0_r8,.000631e0_r8,.000555e0_r8,.000494e0_r8, &
         .000480e0_r8,.000408e0_r8,.000385e0_r8,.000361e0_r8,.000351e0_r8,.000349e0_r8, &
         .086323e0_r8,.028778e0_r8,.010805e0_r8,.004383e0_r8,.002175e0_r8,.001386e0_r8, &
         .000923e0_r8,.000671e0_r8,.000640e0_r8,.000587e0_r8,.000525e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .089886e0_r8,.026029e0_r8,.008558e0_r8,.003312e0_r8,.001655e0_r8,.001124e0_r8, &
         .000807e0_r8,.000631e0_r8,.000602e0_r8,.000568e0_r8,.000525e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .093449e0_r8,.023279e0_r8,.006312e0_r8,.002240e0_r8,.001135e0_r8,.000862e0_r8, &
         .000692e0_r8,.000591e0_r8,.000564e0_r8,.000549e0_r8,.000525e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .094163e0_r8,.020198e0_r8,.004594e0_r8,.001772e0_r8,.001077e0_r8,.000806e0_r8, &
         .000649e0_r8,.000565e0_r8,.000547e0_r8,.000537e0_r8,.000521e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8/), &
         (/18,6/))
    ozone(1:18, 13:18, 4) = RESHAPE( (/ &
         .094877e0_r8,.017116e0_r8,.002877e0_r8,.001305e0_r8,.001018e0_r8,.000751e0_r8, &
         .000606e0_r8,.000539e0_r8,.000529e0_r8,.000524e0_r8,.000516e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .100592e0_r8,.015718e0_r8,.002411e0_r8,.001007e0_r8,.000802e0_r8,.000642e0_r8, &
         .000559e0_r8,.000526e0_r8,.000501e0_r8,.000474e0_r8,.000470e0_r8,.000439e0_r8, &
         .000430e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .106308e0_r8,.014320e0_r8,.001946e0_r8,.000709e0_r8,.000585e0_r8,.000533e0_r8, &
         .000512e0_r8,.000512e0_r8,.000473e0_r8,.000425e0_r8,.000423e0_r8,.000342e0_r8, &
         .000301e0_r8,.000257e0_r8,.000232e0_r8,.000212e0_r8,.000205e0_r8,.000209e0_r8, &
         .107021e0_r8,.012233e0_r8,.001533e0_r8,.000643e0_r8,.000556e0_r8,.000505e0_r8, &
         .000471e0_r8,.000448e0_r8,.000403e0_r8,.000362e0_r8,.000355e0_r8,.000296e0_r8, &
         .000251e0_r8,.000207e0_r8,.000180e0_r8,.000161e0_r8,.000152e0_r8,.000140e0_r8, &
         .107735e0_r8,.010146e0_r8,.001120e0_r8,.000576e0_r8,.000526e0_r8,.000477e0_r8, &
         .000430e0_r8,.000385e0_r8,.000332e0_r8,.000300e0_r8,.000288e0_r8,.000249e0_r8, &
         .000202e0_r8,.000158e0_r8,.000132e0_r8,.000114e0_r8,.000104e0_r8,.000093e0_r8, &
         .107002e0_r8,.009961e0_r8,.001025e0_r8,.000533e0_r8,.000497e0_r8,.000460e0_r8, &
         .000422e0_r8,.000385e0_r8,.000332e0_r8,.000300e0_r8,.000288e0_r8,.000249e0_r8, &
         .000202e0_r8,.000158e0_r8,.000132e0_r8,.000114e0_r8,.000104e0_r8,.000093e0_r8/), &
         (/18,6/))
    ozone(1:18, 19:24, 4) = RESHAPE( (/ &
         .109850e0_r8,.010424e0_r8,.001079e0_r8,.000567e0_r8,.000498e0_r8,.000479e0_r8, &
         .000463e0_r8,.000448e0_r8,.000418e0_r8,.000399e0_r8,.000389e0_r8,.000367e0_r8, &
         .000351e0_r8,.000328e0_r8,.000320e0_r8,.000337e0_r8,.000355e0_r8,.000304e0_r8, &
         .112375e0_r8,.011255e0_r8,.001357e0_r8,.000744e0_r8,.000585e0_r8,.000533e0_r8, &
         .000512e0_r8,.000512e0_r8,.000504e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .111320e0_r8,.011439e0_r8,.001485e0_r8,.000843e0_r8,.000642e0_r8,.000549e0_r8, &
         .000512e0_r8,.000512e0_r8,.000504e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .107719e0_r8,.013052e0_r8,.001822e0_r8,.000953e0_r8,.000701e0_r8,.000604e0_r8, &
         .000551e0_r8,.000525e0_r8,.000509e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .104118e0_r8,.014665e0_r8,.002159e0_r8,.001063e0_r8,.000759e0_r8,.000659e0_r8, &
         .000589e0_r8,.000539e0_r8,.000514e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .104145e0_r8,.018082e0_r8,.003274e0_r8,.001493e0_r8,.000919e0_r8,.000762e0_r8, &
         .000678e0_r8,.000641e0_r8,.000584e0_r8,.000531e0_r8,.000495e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8/), &
         (/18,6/))
    ozone(1:18, 25:30, 4) = RESHAPE( (/ &
         .104172e0_r8,.021499e0_r8,.004389e0_r8,.001922e0_r8,.001078e0_r8,.000865e0_r8, &
         .000767e0_r8,.000743e0_r8,.000654e0_r8,.000562e0_r8,.000499e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .101718e0_r8,.026523e0_r8,.008473e0_r8,.004382e0_r8,.002392e0_r8,.001505e0_r8, &
         .001036e0_r8,.000836e0_r8,.000727e0_r8,.000618e0_r8,.000550e0_r8,.000494e0_r8, &
         .000501e0_r8,.000479e0_r8,.000473e0_r8,.000509e0_r8,.000541e0_r8,.000445e0_r8, &
         .099264e0_r8,.031548e0_r8,.012557e0_r8,.006841e0_r8,.003705e0_r8,.002144e0_r8, &
         .001304e0_r8,.000930e0_r8,.000799e0_r8,.000675e0_r8,.000601e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .097501e0_r8,.035382e0_r8,.014856e0_r8,.008207e0_r8,.004619e0_r8,.002720e0_r8, &
         .001610e0_r8,.001012e0_r8,.000829e0_r8,.000687e0_r8,.000610e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .095737e0_r8,.039215e0_r8,.017155e0_r8,.009572e0_r8,.005532e0_r8,.003296e0_r8, &
         .001916e0_r8,.001094e0_r8,.000858e0_r8,.000699e0_r8,.000618e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .092863e0_r8,.042419e0_r8,.020704e0_r8,.012034e0_r8,.007417e0_r8,.004504e0_r8, &
         .002590e0_r8,.001334e0_r8,.000977e0_r8,.000731e0_r8,.000622e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8/), &
         (/18,6/))
    ozone(1:18, 31:36, 4) = RESHAPE( (/ &
         .089988e0_r8,.045622e0_r8,.024253e0_r8,.014495e0_r8,.009303e0_r8,.005711e0_r8, &
         .003264e0_r8,.001574e0_r8,.001096e0_r8,.000762e0_r8,.000627e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .084591e0_r8,.046553e0_r8,.025037e0_r8,.015156e0_r8,.009841e0_r8,.006124e0_r8, &
         .003534e0_r8,.001693e0_r8,.001170e0_r8,.000793e0_r8,.000631e0_r8,.000537e0_r8, &
         .000551e0_r8,.000509e0_r8,.000486e0_r8,.000516e0_r8,.000548e0_r8,.000446e0_r8, &
         .079194e0_r8,.047483e0_r8,.025822e0_r8,.015818e0_r8,.010380e0_r8,.006537e0_r8, &
         .003804e0_r8,.001811e0_r8,.001245e0_r8,.000825e0_r8,.000635e0_r8,.000570e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .077409e0_r8,.048159e0_r8,.026661e0_r8,.016596e0_r8,.010962e0_r8,.006972e0_r8, &
         .004160e0_r8,.002132e0_r8,.001391e0_r8,.000868e0_r8,.000686e0_r8,.000601e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .075625e0_r8,.048835e0_r8,.027500e0_r8,.017375e0_r8,.011544e0_r8,.007407e0_r8, &
         .004516e0_r8,.002453e0_r8,.001536e0_r8,.000912e0_r8,.000737e0_r8,.000632e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .074927e0_r8,.049459e0_r8,.029215e0_r8,.018025e0_r8,.011754e0_r8,.007786e0_r8, &
         .004972e0_r8,.002926e0_r8,.001817e0_r8,.001025e0_r8,.000758e0_r8,.000632e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8/), &
         (/18,6/))
    ozone(1:18, 37, 4) = (/ &
         .074229e0_r8,.050084e0_r8,.030930e0_r8,.018676e0_r8,.011965e0_r8,.008165e0_r8, &
         .005428e0_r8,.003399e0_r8,.002098e0_r8,.001138e0_r8,.000780e0_r8,.000632e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8/)

    ozsig(:) = (/ &
         .020747_r8,.073986_r8,.124402_r8,.174576_r8,.224668_r8,.274735_r8, &
         .324767_r8,.374806_r8,.424818_r8,.497450_r8,.593540_r8,.688125_r8, &
         .777224_r8,.856317_r8,.920400_r8,.960480_r8,.981488_r8,.995004_r8/)


    first_getoz=.TRUE.

    IF(first_getoz)THEN
       mon_getoz=yrl/12.0_r8
       year_getoz=yrl
       IF(nlm_getoz.NE.kmax)THEN
          inter_getoz=.TRUE.
       ELSE
          inter_getoz=.FALSE.
          DO l=1,nlm_getoz
             ll=nlm_getoz-l+1
             IF(ABS(ozsig(l)-sl(ll)).GT.0.0001_r8) inter_getoz=.TRUE.
          END DO
       ENDIF
       !print*,"inter_getoz=",inter_getoz
       first_getoz=.FALSE.
    ENDIF
  END SUBROUTINE InitGetoz_rrtm

  SUBROUTINE getoz_rrtm (kmax,sigmid,colrad,date,o3l)
    !
    ! input parameters and variables:
    !     ncols  =  number of atmospheric columns
    !     kmax   =  number of atmospheric layers
    !     sigmid =  sigma coordinate at middle of layer
    !     colrad =  colatitude of each column (0-3.14 from np to sp in radians)
    !     date   =  model julian date
    !
    ! tabulated data
    !     ozone  =  climatological ozone mixing ratio in 18 sigma layers
    !               and in 5 degree latitude interval
    !
    ! output variables:
    !     o3l   =  18 layers ozone mixing ratio in given lat and date
    !
    !==========================================================================
    ! :: kmax.....Number of grid points at vertical
    ! :: sigmid.......sigma coordinate at middle of layer
    ! :: pai......constant pi=3.1415926
    ! :: yrl......length of year in days
    !==========================================================================

    ! Input variables
    !INTEGER      ,    INTENT(IN   ) :: ncols,adjNCols
    INTEGER      ,    INTENT(IN   ) :: kmax
    REAL(KIND=r8),    INTENT(IN   ) :: sigmid (kmax)
    REAL(KIND=r8),    INTENT(IN   ) :: colrad
    REAL(KIND=r8),    INTENT(INOUT) :: date
    REAL(KIND=r8),PARAMETER :: pai=3.1415926e00_r8

    ! Output variable
    REAL(KIND=r8),    INTENT(OUT) :: o3l(kmax)

    ! Local variables
    REAL(KIND=r8) :: a1   (nlm_getoz)
    REAL(KIND=r8) :: a2   (nlm_getoz)
    REAL(KIND=r8) :: a3   (nlm_getoz)
    REAL(KIND=r8) :: a4   (nlm_getoz)
    REAL(KIND=r8) :: b1   (nlm_getoz)
    REAL(KIND=r8) :: b2   (nlm_getoz)
    REAL(KIND=r8) :: b3   (nlm_getoz)
    REAL(KIND=r8) :: b4   (nlm_getoz)
    REAL(KIND=r8) :: do3a (nlm_getoz)
    REAL(KIND=r8) :: do3b (nlm_getoz)
    REAL(KIND=r8) :: ozo3l(nlm_getoz)

    REAL(KIND=r8), PARAMETER :: rlag = 14.8125e0_r8

    INTEGER :: l
    INTEGER :: la
    INTEGER :: ll
    INTEGER :: kmx
    INTEGER :: imon
    INTEGER :: isea
    INTEGER :: k
    INTEGER :: i
    INTEGER :: kk    (kmax)
    REAL(KIND=r8)    :: theta
    REAL(KIND=r8)    :: flat
    REAL(KIND=r8)    :: rang
    REAL(KIND=r8)    :: rsin1
    REAL(KIND=r8)    :: rcos1
    REAL(KIND=r8)    :: rcos2
    REAL(KIND=r8)    :: rate
    REAL(KIND=r8)    :: aa
    REAL(KIND=r8)    :: bb
    LOGICAL :: notfound(kMax)

    kmx=nlm_getoz
    !
    !     find closest place in the data according to input slat.
    !
    IF(date.GT.year_getoz) date=date-year_getoz

    imon=date/mon_getoz + 1

    IF(imon.LT.1)imon=1

    isea=imon/3 + 1

    IF(isea.EQ.5) isea=1
    IF(isea.GT.5) THEN
       WRITE(12,"('0 ERROR IN ISEA - TERMINATION IN SUBROUTINE GETOZ')")
       WRITE(10,"('0 ERROR IN ISEA - TERMINATION IN SUBROUTINE GETOZ')")
       STOP 9954
    END IF
    !DO i=1,adjNCols
       theta = 90.0_r8-(180.0_r8/pai)*colrad ! colatitude -> latitude
! if(nmachs==1 .and. i==1318) &
!    write (80+nmachs,fmt='(A,2(F30.15,1X),I3.3)') 'theta,colrad(i),kmx=',theta,colrad(i),kmx
! if(nmachs==2 .and. mynum==2 .and. i==651) &
!    write (80+nmachs,fmt='(A,2(F30.15,1X),I3.3)') 'theta,colrad(i),kmx=',theta,colrad(i),kmx
       ! the 180 degrees are divided into 37 bands with 5deg each
       ! except for the first and last, which have 2.5 deg
       ! The centers of the bands are located at:
       !   90, 85, 80, ..., 5, 0, -5, ..., -85, -90 (37 latitudes)
       flat  = 0.2_r8*theta ! indexing the latitudes: goes from -18. to +18.
       ! find the latitude index before and after each latitude
       la    = 19.501e0_r8-flat !
       ll    = 19.001e0_r8-flat

       !
       !     find sin and cos coefficients for time interpolation.
       !
       rang=2.0e0_r8*pai*(date-rlag)/year_getoz
       rsin1=SIN(rang)
       rcos1=COS(rang)
       rcos2=COS(2.0e0_r8*rang)
       rate=REAL(19-ll,r8)-flat
       !
       !     ozone interpolation in latitude and time
       !
    !END DO
    DO k=1,kmx
       !DO i=1,adjNCols

          a1(k) =2.5e-1_r8*(ozone(k,la,1)+ozone(k,la,2)+ &
               ozone(k,la,3)+ozone(k,la,4))
          a2(k) =0.5e0_r8*(ozone(k,la,2)-ozone(k,la,4))
          a3(k) =0.5e0_r8*(ozone(k,la,1)-ozone(k,la,3))
          a4(k) =2.5e-1_r8*(ozone(k,la,1)+ozone(k,la,3)- &
               ozone(k,la,2)-ozone(k,la,4))
          b1(k) =2.5e-1_r8*(ozone(k,ll,1)+ozone(k,ll,2)+ &
               ozone(k,ll,3)+ozone(k,ll,4))
          b2(k) =0.5e0_r8*(ozone(k,ll,2)-ozone(k,ll,4))
          b3(k) =0.5e0_r8*(ozone(k,ll,1)-ozone(k,ll,3))
          b4(k) =2.5e-1_r8*(ozone(k,ll,1)+ozone(k,ll,3)- &
               ozone(k,ll,2)-ozone(k,ll,4))
          do3a(k)=a1(k)+rsin1*a2(k)+rcos1*a3(k)+rcos2*a4(k)
          do3b(k)=b1(k)+rsin1*b2(k)+rcos1*b3(k)+rcos2*b4(k)
          ozo3l(k)=do3a(k)+rate*(do3b(k)-do3a(k))
          ozo3l(k)=1.0e-04_r8*ozo3l(k)
       !END DO
    END DO
    !print *,'LFR->Oz: ', (ozo3l(1,k),k=1,kmx)
    IF(inter_getoz)THEN
       DO l=1,kmax
!          print *,'LFR->MQB: ',sigmid(l),ozsig(1),sigmid(l) > ozsig(1)
          notfound(l) = sigmid(l) > ozsig(1)
          IF (notfound(l)) THEN
             kk(l)=kmx
          ELSE
             kk(l)=2
          END IF
!	  print *,'lfr->1.kk(1,l) :',l,kk(1,l)
       END DO
       DO l=1,kmax
          IF (notfound(l)) THEN
             DO k=2,kmx
                IF(sigmid(l).GT.ozsig(k-1).AND.sigmid(l).LE.ozsig(k))THEN
                   kk(l)=k
                   EXIT
                END IF
             END DO
          END IF
!	  print *,'lfr->2.kk(1,l) :',l,kk(1,l)

       END DO
       !DO l = 1, kmax
          !DO i= 2, adjNCols
             !kk(l) = kk(l)
          !END DO
!	  print *,'lfr->3.kk(1,l) :',l,kk(1,l)
       !END DO
    END IF
!    print *,'LFR->inter_getoz=',inter_getoz
    IF(inter_getoz)THEN
       DO l=1,kmax
          !DO i=1,adjNCols
             aa=(ozo3l(kk(l))-ozo3l(kk(l)-1))/(ozsig(kk(l))-ozsig(kk(l)-1))
             bb=ozo3l(kk(l)-1)-aa*ozsig(kk(l)-1)
             o3l(kmax+1-l)=bb+aa*sigmid(l)
	     !print *,'LFR->',l,i,kmax+1-l,aa,bb,o3l(i,kmax+1-l)
          !END DO
       END DO
    END IF
    IF(.NOT.inter_getoz)THEN
       DO l=1,nlm_getoz
          !DO i=1,adjNCols
             o3l(l)=ozo3l(l)
          !END DO
       END DO
    ENDIF
!    stop 362
  END SUBROUTINE getoz_rrtm
END MODULE rrtm_driv

!.. subroutine dumpBugRad(var,mzp,i,j)
  !.. use node_mod, only:  &
       !.. nodei0, &
       !.. nodej0, mynum, nmachs
  !.. use mem_grid, only: time
  !.. integer, intent(in) :: i,j,mzp
  !.. real,intent(in) :: var(mzp)

  !.. if(nmachs==1) then
    !.. write (70,fmt='(2(I2.2,1X),35(E16.8,1X))') i+nodei0(mynum,1),j+nodej0(mynum,1),var
  !.. else
    !.. write (70+mynum,fmt='(2(I2.2,1X),35(E16.8,1X))') i+nodei0(mynum,1),j+nodej0(mynum,1),var
  !.. endif

!.. end subroutine dumpBugRad
