
MODULE carma_driv
    USE mem_tend, ONLY: tend ! INTENT(INOUT)
    USE mem_radiate, ONLY: &
         ilwrtyp, iswrtyp, &       ! INTENT(IN)
         radiate_g,        &       ! INTENT(INOUT)
         radfrq, &                 ! INTENT(IN)
         ncall_i, prsnz, prsnzp, & ! INTENT(INOUT)
         lonrad                    ! INTENT(IN)
    USE mem_micro, ONLY: micro_g,micro_vars ! INTENT(INOUT)
    USE mem_basic, ONLY: basic_g,basic_vars  ! INTENT(INOUT)
    USE micphys, ONLY: &
         gnu, level, icloud, irain, ipris, & ! INTENT(IN)
         isnow, iaggr, igraup, ihail,mcphys_type         ! INTENT(IN)

    USE mem_cuparm, ONLY: cuparm_g, cuparm_vars, nnqparm  ! INTENT(IN)
    USE mem_leaf, ONLY: leaf_g ! INTENT(IN)
    use rconstants  , only : &
                   cp,cpor,p00,stefan,cpi,pi180,solar

    USE mem_grid, ONLY: &
         ngrid, time, dtlt, itime1, nzg, nzs, npatch, grid_g, & ! INTENT(IN)
         nnzp, if_adap, zm, zt, naddsc, nzpmax, imonth1,      & ! INTENT(IN)
         idate1, iyear1, centlat, centlon, ztop,dzm ,ngrids     ! INTENT(IN)
    use mem_scratch1_grell, only: &
                  ierr4d,xmb4d,zup5d,clwup5d
    private
    public carma_driver

CONTAINS

  SUBROUTINE carma_driver(mzp, mxp, myp, ia, iz, ja, jz, mynum)



!    USE mem_leaf, ONLY: leaf_g ! INTENT(IN)


    USE rad_carma, ONLY: radcomp_carma 
    
    use teb_spm_start, only: TEB_SPM ! INTENT(IN)
    use mem_teb_common, only: tebc_g ! INTENT(INOUT)
                           
!    USE rconstants , ONLY : solar
    
    USE mem_carma, only: carma_aotMap,solfac
  
    implicit none
    ! arguments:
    integer, intent(in) :: mzp, mxp, myp, ia, iz, ja, jz, mynum
    ! local variables:
    real :: hranglelocal
    real :: solc
    real :: maxcloud_fraction
    real,allocatable,dimension(:,:,:) :: lwl,ice_frac,iwl,cloud_fraction
    real,allocatable,dimension(:,:) :: rain
    real :: dummy
    integer :: i, j, k, ncols
    ! teb_spm
    real, pointer :: emis_town(:,:), alb_town(:,:), ts_town(:,:), g_urban(:,:,:)
    !
    real      :: cdec
    integer   :: jday

    !srf -03/03/2013 cloud fraction
    real, parameter :: fxx = 0.2 

    !- if not including radiation, return
    IF ((ilwrtyp + iswrtyp)==0) RETURN
    !   
    !--- apply radiative tendencies to model tendencies
    CALL tend_accum(mzp, mxp, myp, ia, iz, ja, jz)

    !--- radiation is called on each radfrq seconds
    IF (.not. (mod(time+.001, radfrq) < dtlt .or. time<0.001)) return

      !- TEB_SPM
     	if (TEB_SPM==1) then
     	   EMIS_TOWN => tebc_g(ngrid)%EMIS_TOWN
     	   ALB_TOWN  => tebc_g(ngrid)%ALB_TOWN
     	   TS_TOWN   => tebc_g(ngrid)%TS_TOWN
     	   G_URBAN   => leaf_g(ngrid)%G_URBAN
     	else
     	   nullify(EMIS_TOWN)
     	   nullify(ALB_TOWN)
     	   nullify(TS_TOWN)
     	   nullify(G_URBAN)
     	endif

      ! Compute solar zenith angle, multiplier for solar constant, sfc albeDO,
      ! and surface upward longwave radiation.
 
       if (TEB_SPM==1) then 
          CALL radprep(TEB_SPM, imonth1, idate1, iyear1, time, itime1, &
               centlat, centlon, lonrad, pi180,                        &
               nzg, nzs, npatch, ia, iz, ja, jz, jday,       &
               leaf_g(ngrid)%soil_water,                               &
               leaf_g(ngrid)%soil_energy,                              &
               leaf_g(ngrid)%soil_text,                                &
               leaf_g(ngrid)%sfcwater_energy,                          &
               leaf_g(ngrid)%sfcwater_depth,                           &
               leaf_g(ngrid)%leaf_class,                               &
               leaf_g(ngrid)%veg_fracarea,                             &
               leaf_g(ngrid)%veg_height,                               &
               leaf_g(ngrid)%veg_albedo,                               &
               leaf_g(ngrid)%patch_area,                               &
               leaf_g(ngrid)%sfcwater_nlev,                            &
               leaf_g(ngrid)%veg_temp,                                 &
               leaf_g(ngrid)%can_temp,                                 &
               solfac,                                                 &
               grid_g(ngrid)%glat,                                     &
               grid_g(ngrid)%glon,                                     &
               radiate_g(ngrid)%rshort,                                &
               radiate_g(ngrid)%rlong,                                 &
               radiate_g(ngrid)%rlongup,                               &
               radiate_g(ngrid)%albedt,                                &
               radiate_g(ngrid)%cosz,                                  &
               hrAngleLocal,                                           &	       
               cdec,                                                   &
               ! TEB_SPM
               EMIS_TOWN, ALB_TOWN, TS_TOWN, G_URBAN                   )
       else
          
           CALL radprep(TEB_SPM, imonth1, idate1, iyear1, time, itime1, &
               centlat, centlon, lonrad, pi180,                        &
               nzg, nzs, npatch, ia, iz, ja, jz, jday,       &
               leaf_g(ngrid)%soil_water,                               &
               leaf_g(ngrid)%soil_energy,                              &
               leaf_g(ngrid)%soil_text,                                &
               leaf_g(ngrid)%sfcwater_energy,                          &
               leaf_g(ngrid)%sfcwater_depth,                           &
               leaf_g(ngrid)%leaf_class,                               &
               leaf_g(ngrid)%veg_fracarea,                             &
               leaf_g(ngrid)%veg_height,                               &
               leaf_g(ngrid)%veg_albedo,                               &
               leaf_g(ngrid)%patch_area,                               &
               leaf_g(ngrid)%sfcwater_nlev,                            &
               leaf_g(ngrid)%veg_temp,                                 &
               leaf_g(ngrid)%can_temp,                                 &
               solfac,                                                 &
               grid_g(ngrid)%glat,                                     &
               grid_g(ngrid)%glon,                                     &
               radiate_g(ngrid)%rshort,                                &
               radiate_g(ngrid)%rlong,                                 &
               radiate_g(ngrid)%rlongup,                               &
               radiate_g(ngrid)%albedt,                                &
               radiate_g(ngrid)%cosz,                                  &
               hrAngleLocal,                                           &	       
               cdec                                                    )
       endif


       !--- set radiation tendency for theta to zero
       radiate_g(ngrid)%fthrd(1:mzp,1:mxp,1:myp) = 0.0


       !--- get cloud properties         
       !--- alocate arrays for cloud/radiation interaction
	allocate(lwl(mzp,mxp,myp),ice_frac      (mzp,mxp,myp) &	   
                ,iwl(mzp,mxp,myp),cloud_fraction(mzp,mxp,myp))	
        allocate(rain(mxp,myp))						   	
        !- set zero to the local arrays:				   		
        cloud_fraction=0.0  ;ice_frac =0.0; rain=0.0 ; lwl =0.0 ;iwl =0.0 							   		

	call  cloud_prop_carma(mzp, mxp, myp, ia, iz, ja, jz &
		    !-- output  		  
		    ,cloud_fraction	  &
		    ,rain		  &
		    ,lwl		  &
		    ,iwl		  &
		    ,ice_frac		  &
		    )
       
       !- CARMA Radiation

       !-srf - including a 'cloud fraction" :
       !- included 01.03.2013 to improve near surf temp and precip
       DO i=1,mxp
          DO j=1,myp
            maxCloud_fraction=maxval(cloud_fraction(:,i,j))
            RAIN(i,j)= (1.0 - maxCloud_fraction) * RAIN(i,j)
         END DO
       END DO
       !
       LWL = (1.0 - cloud_fraction) * LWL
       IWL = (1.0 - cloud_fraction) * IWL
       !-srf - including a 'cloud fraction" :
       !- included 01.03.2013 to improve near surf temp and precip
       !rain= fxx * rain
       !lwl = fxx * lwl
       !iwl = fxx * iwl


       CALL radcomp_carma(mzp,mxp,myp,ia,iz,ja,jz,solfac  &
               ,basic_g(ngrid)%theta	   &
               ,basic_g(ngrid)%pi0	   &
               ,basic_g(ngrid)%pp	   &
               ,basic_g(ngrid)%rv	   &
               ,RAIN,LWL,IWL		   &
               ,basic_g(ngrid)%dn0	   &
               ,basic_g(ngrid)%rtp	   &
               ,radiate_g(ngrid)%fthrd     &
               ,grid_g(ngrid)%rtgt	   &
               ,grid_g(ngrid)%f13t	   &
               ,grid_g(ngrid)%f23t	   &
               ,grid_g(ngrid)%glat	   &
               ,grid_g(ngrid)%glon	   &
               ,radiate_g(ngrid)%rshort    &
               ,radiate_g(ngrid)%rlong     &
               ,radiate_g(ngrid)%albedt    &
               ,radiate_g(ngrid)%cosz	   &
               ,radiate_g(ngrid)%rlongup   &
               ,mynum			   &
               ,grid_g(ngrid)%fmapt	   &
               ,leaf_g(ngrid)%patch_area   &
               ,npatch  		   &
               ,hrAngleLocal		   &
	       ,carma_aotMap(ngrid)%aotMap &
               )

!    if(mynum== 5) then
!     write(mynum+1,*) "============= radiation-carma ==================="
!     write(mynum+1,*) "mynum=",mynum
!     write(mynum+1,*) "max/min rlong ",maxval(radiate_g(ngrid)%rlong),minval(radiate_g(ngrid)%rlong)
!     write(mynum+1,*) "max/min rshort",maxval(radiate_g(ngrid)%rshort),minval(radiate_g(ngrid)%rshort)
!     write(mynum+1,*) "max/min fthrd ",maxval(86400.*radiate_g(ngrid)%fthrd ),&
!                  		      minval(86400.*radiate_g(ngrid)%fthrd )
!     write(mynum+1,*) "max/min theta ",maxval(basic_g(ngrid)%theta),minval(basic_g(ngrid)%theta)
!     write(mynum+1,*) "============= radiation-carma ===================="
!     call flush(mynum+1)
!    endif

    
       deallocate(lwl,ice_frac,iwl,cloud_fraction,rain)
  
  END SUBROUTINE carma_driver

! ****************************************************************************

  SUBROUTINE radprep(TEB_SPM, imonth1, idate1, iyear1, time, itime1,        &
       centlat, centlon, lonrad, pi180,                                     &
       mzg, mzs, np, ia, iz, ja, jz, jday,                          &
       soil_water, soil_energy, soil_text, sfcwater_energy, sfcwater_depth, &
       leaf_class, veg_fracarea, veg_height, veg_albedo, patch_area,        &
       sfcwater_nlev, veg_temp, can_temp,                                   &
       solfac, glat, glon, rshort, rlong, rlongup, albedt, cosz,            &
       hrAngleLocal,                                                        &
       cdec,                                                                &
       EMIS_TOWN, ALB_TOWN, TS_TOWN, G_URBAN                                &
       !
       )

    use mem_leaf, only: isfcl  !DSM
  
    IMPLICIT NONE

    ! Arguments:
    INTEGER, INTENT(IN)      :: TEB_SPM
    INTEGER, INTENT(IN)      :: imonth1, idate1, iyear1, itime1
    REAL, INTENT(IN)         :: time
    REAL, INTENT(IN)         :: centlat(:), centlon(:)
    INTEGER, INTENT(IN)      :: lonrad
    REAL, INTENT(IN)         :: pi180
    INTEGER, INTENT(IN)      :: mzg, mzs, np, ia, iz, ja, jz
    INTEGER, INTENT(OUT)     :: jday
    !DIMENSION(mzg,m2,m3,np)
    REAL, INTENT(IN)         :: soil_water(:,:,:,:), soil_energy(:,:,:,:), &
         soil_text(:,:,:,:)
    !DIMENSION(mzs,m2,m3,np)
    REAL, INTENT(IN)         :: sfcwater_energy(:,:,:,:), &
         sfcwater_depth(:,:,:,:)
    !DIMENSION(m2,m3,np)
    REAL, INTENT(IN)         :: leaf_class(:,:,:), veg_fracarea(:,:,:), &
         veg_height(:,:,:), veg_albedo(:,:,:), patch_area(:,:,:),       &
         sfcwater_nlev(:,:,:), veg_temp(:,:,:), can_temp(:,:,:)
    !TEB_SPM
    !DIMENSION(m2,m3)
    REAL, pointer, optional :: EMIS_TOWN(:,:), ALB_TOWN(:,:), TS_TOWN(:,:)
    !DIMENSION(m2,m3,np)
    real, pointer, optional :: G_URBAN(:,:,:)
    !
    REAL, INTENT(OUT)        :: solfac
    !DIMENSION(m2,m3)
    REAL, INTENT(IN)         :: glat(:,:), glon(:,:), rshort(:,:), rlong(:,:)
    REAL, INTENT(OUT)        :: rlongup(:,:), albedt(:,:)
    REAL, INTENT(INOUT)      :: cosz(:,:)
    REAL, INTENT(OUT)        :: cdec
    
    REAL, INTENT(OUT)        :: hrAngleLocal
        
    ! Local Variables
    real, pointer :: L_EMIS_TOWN, L_ALB_TOWN, L_TS_TOWN, L_G_URBAN
    !
    INTEGER :: ip, i, j

    INTERFACE

       SUBROUTINE sfcrad(mzg, mzs, ip,                  &
            soil_energy, soil_water, soil_text,         &
            sfcwater_energy, sfcwater_depth,            &
            patch_area, can_temp, veg_temp, leaf_class, &
            veg_height, veg_fracarea, veg_albedo,       &
            sfcwater_nlev, rshort, rlong, albedt,       &
            rlongup, cosz,                              &
            G_URBAN, ETOWN, ALBTOWN, TSTOWN             )
         integer, intent(IN) :: mzg, mzs, ip
         real, intent(IN)    :: soil_energy(mzg)
         real, intent(IN)    :: soil_water(mzg)
         real, intent(IN)    :: soil_text(mzg)
         real, intent(IN)    :: sfcwater_energy(mzs)
         real, intent(IN)    :: sfcwater_depth(mzs)
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
         real, pointer, optional :: G_URBAN, ETOWN, ALBTOWN, TSTOWN
       END SUBROUTINE sfcrad
         
    END INTERFACE

    ! TEB_SPM
    nullify(L_EMIS_TOWN)
    nullify(L_ALB_TOWN)
    nullify(L_TS_TOWN)
    nullify(L_G_URBAN)

    ! Compute solar zenith angle [cosz(i,j)] & solar constant factr [solfac].

    CALL zen(imonth1, idate1, iyear1, time, itime1, centlat, centlon, &
         lonrad, pi180, ia, iz, ja, jz, jday, glat, glon, cosz, &
         solfac, hrAngleLocal, cdec)

    ! Compute patch-averaged surface albeDO [albedt(i,j)] and up longwave
    ! radiative flux [rlongup(i,j)].

  if (isfcl /= 5 .or. time==0.) then !DSM ---- o rlongup e o albedt deve vir do proprio jules
    albedt  = 0.
    rlongup = 0.

    DO ip = 1,np
       DO j = 1,jz
          DO i = 1,iz

             ! TEB_SPM
             if (TEB_SPM==1 .and. present(G_URBAN) .and.           &
                  present(EMIS_TOWN) .and. present(ALB_TOWN) .and. &
                  present(TS_TOWN))                                then
                L_G_URBAN   => G_URBAN(i,j,ip)
                L_EMIS_TOWN => EMIS_TOWN(i,j)
                L_ALB_TOWN  => ALB_TOWN(i,j)
                L_TS_TOWN   => TS_TOWN(i,j)
                CALL sfcrad(mzg, mzs, ip,                                     &
                     soil_energy(1:mzg,i,j,ip), soil_water(1:mzg,i,j,ip),     &
                     soil_text(1:mzg,i,j,ip),   sfcwater_energy(1:mzs,i,j,ip),&
                     sfcwater_depth(1:mzs,i,j,ip), patch_area(i,j,ip),        &
                     can_temp(i,j,ip),         veg_temp(i,j,ip),              &
                     leaf_class(i,j,ip),       veg_height(i,j,ip),            &
                     veg_fracarea(i,j,ip),     veg_albedo(i,j,ip),            &
                     sfcwater_nlev(i,j,ip),                                   &
                     rshort(i,j), rlong(i,j), albedt(i,j),                    &
                     rlongup(i,j), cosz(i,j),                                 &
                     ! TEB_SPM
                     L_G_URBAN, L_EMIS_TOWN, L_ALB_TOWN, L_TS_TOWN            &
                     !
                     )
             else
                CALL sfcrad(mzg, mzs, ip,                                     &
                     soil_energy(1:mzg,i,j,ip), soil_water(1:mzg,i,j,ip),     &
                     soil_text(1:mzg,i,j,ip),   sfcwater_energy(1:mzs,i,j,ip),&
                     sfcwater_depth(1:mzs,i,j,ip), patch_area(i,j,ip),        &
                     can_temp(i,j,ip),         veg_temp(i,j,ip),              &
                     leaf_class(i,j,ip),       veg_height(i,j,ip),            &
                     veg_fracarea(i,j,ip),     veg_albedo(i,j,ip),            &
                     sfcwater_nlev(i,j,ip),                                   &
                     rshort(i,j), rlong(i,j), albedt(i,j),                    &
                     rlongup(i,j), cosz(i,j)                                  &
                     )
             endif

          END DO
       END DO
    END DO

  endif  !DSM
  
  END SUBROUTINE radprep

  ! ***************************************************************************

  SUBROUTINE zen(imonth1, idate1, iyear1, time, itime1, centlat, centlon, &
       lonrad, pi180, &
      ia, iz, ja, jz, jday, glat, glon, cosz, solfac, hrangle,&
      cdec)

    USE ModDateUtils
    USE mem_radiate, only: radfrq
    use node_mod, only: mynum ! INTENT(IN)

    IMPLICIT NONE
    ! Arguments:
    INTEGER, INTENT(IN)  :: imonth1, idate1, iyear1, itime1
    REAL, INTENT(IN)     :: time
    REAL, INTENT(IN)     :: centlat(:), centlon(:)
    INTEGER, INTENT(IN)  :: lonrad
    REAL, INTENT(IN)     :: pi180
    INTEGER, INTENT(IN)  :: ia, iz, ja, jz
    INTEGER, INTENT(OUT) :: jday
    REAL, INTENT(IN)     :: glat(:,:) !(m2,m3)
    REAL, INTENT(IN)     :: glon(:,:) !(m2,m3)
    REAL, INTENT(INOUT)  :: cosz(:,:) !(m2,m3)
    REAL, INTENT(OUT)    :: solfac

    REAL, INTENT(OUT)    :: hrangle

    REAL, INTENT(OUT)    :: cdec
    ! Local Variables:
    INTEGER :: i, j! , julday
    REAL    :: sdec, declin, d0, d02, dayhr, radlat, cslcsd, snlsnd, gglon, &
         dayhrr !,tdec
 
    REAL :: eqt

    !common /radcom/ tdec,sdec,cdec,declin,rvr,rtr,dn0r,pird,prd,temprd  &
    !     ,fthrl,dzmr,dztr,fthrs

    jday   = julday(imonth1, idate1, iyear1)
    jday   = jday + nint(time/86400.)
    !      sdec - sine of declination, cdec - cosine of declination
    declin = -23.5*cos(6.283/365.*(jday + 9))*pi180
    sdec   = sin(declin)
    cdec   = cos(declin)

    ! Find the factor, solfac, to multiply the solar constant to correct
    ! for Earth's varying distance to the sun.

    d0     = 6.2831853*float(jday-1)/365.
    d02    = d0*2.
    solfac = 1.000110 + 0.034221*cos(d0) + 0.001280*sin(d0) + &
         0.000719*cos(d02) + 0.000077*sin(d02)

    ! Find the hour angle, THEN get cosine of zenith angle.

!--(DMK-CCATT-INI)-----------------------------------------------------
    !NER_i - including solar time equation ("eqt" must be defined, it is a new variable)
    eqt = (0.000075 + 0.001868*cos(d0) - 0.032077*sin(d0) - 0.014615*cos(d02) &
 	  - 0.040849*sin(d02))*1440/(2*3.141593)	   
    !NER_f - including solar time equation	
  
    dayhr = (time / 3600. + float(itime1/100) + float(mod(itime1,100)) / 60.)+ (radfrq/(2*3600))
    !NER (radfrq/(2*3600)) - rad transfer shift half of radfrq(improving rad tendency representativity)
!--(DMK-CCATT-OLD)-----------------------------------------------------
!    dayhr  = time/3600. + float(itime1/100) + float(mod(itime1,100))/60.
!--(DMK-CCATT-FIM)-----------------------------------------------------

    DO j = ja,jz
       DO i = ia,iz
          radlat = glat(i,j)*pi180
          IF (lonrad==0)      radlat = centlat(1)*pi180
          IF (radlat==declin) radlat = radlat + 1.e-5
          cslcsd = cos(radlat)*cdec
          snlsnd = sin(radlat)*sdec
          gglon  = glon(i,j)
          IF (lonrad==0)      gglon = centlon(1)

!--(DMK-CCATT-INI)-----------------------------------------------------
	!NER_i new hour angle calculation		
	hrangle=((dayhr+(gglon/15.)-(12.-eqt/60.))*15./1.)*3.141593/180.	
!--(DMK-CCATT-OLD)-----------------------------------------------------
!          dayhrr    = mod(dayhr+gglon/15.+24., 24.)
!          hrangl    = 15.*(dayhrr - 12.)*pi180
!--(DMK-CCATT-FIM)-----------------------------------------------------

          cosz(i,j) = snlsnd + cslcsd*cos(hrangle)
          !-srf - cosz > 1 no SX6

!--(DMK-CCATT-OLD)-----------------------------------------------------
!          cosz(i,j) = min(cosz(i,j)+1.0E-10, 1.0) !LFR: Prevent 90 degrees z angle
!--(DMK-CCATT-FIM)-----------------------------------------------------

          !cosz(i,j) = min(cosz(i,j), 1.0)
          !cosz(i,j) = max(cosz(i,j),-1.0)

!--(DMK-CCATT-INI)-----------------------------------------------------
          cosz(i,j) = max(cosz(i,j),1.0E-10)
!--(DMK-CCATT-FIM)-----------------------------------------------------

          !-srf        
       END DO
    END DO

  END SUBROUTINE zen
!--------------------------------------------------------------------------------

  SUBROUTINE tend_accum(m1,m2,m3,ia,iz,ja,jz)
      integer, intent(in) :: m1, m2, m3, ia, iz, ja, jz

      ! local variables:
      integer :: i, j, k,ipos
     
      ipos=0
      !do j=ja,jz
     	! do i=ia,iz
      do j=1,m3
     	 do i=1,m2
     	    do k=1,m1
     	      ipos=ipos+1	 
     	      tend%tht(ipos) = tend%tht(ipos) + radiate_g(ngrid)%fthrd(k,i,j)
     	    end do
     	 end do
      end do


 END SUBROUTINE tend_accum

!--------------------------------------------------------------------------------
subroutine cloud_prop_carma(m1,m2,m3,ia, iz, ja, jz &
                     ,cloud_fraction &
                     , rain             &
                     , lwl             &
                     , iwl             &
                     , ice_frac      &
                       )
  
  implicit none
  integer, intent(in) :: m1,m2,m3,ia,iz,ja,jz

  real, intent(out), dimension(m1,m2,m3) :: cloud_fraction !cloud_fraction
  real, intent(out), dimension(m2,m3   ) :: rain !total rain water 
  real, intent(out), dimension(m1,m2,m3) :: lwl !total cloud liquid water (kg/kg for carma and g/m2 for rrtm)
  real, intent(out), dimension(m1,m2,m3) :: iwl !total cloud ice water (kg/kg for carma and g/m2 for rrtm)
  real, intent(out), dimension(m1,m2,m3) :: ice_frac !total cloud ice water 

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
             ,rhwght,rhdif,rhlim,ps,thetas,dthdpmn,dthdp,dummy,bb
  logical cldbnd          ! region below high cloud boundary
  real, external :: rs
  character(len=3) :: cm1
  !- tuning parameters to include direct coupling between cupar and radiation
  real, parameter :: tun_rad_shall=0.02, tun_rad_deep =0.05
! real, parameter :: tun_rad_shall=0.05, tun_rad_deep =0.15
! real, parameter :: tun_rad_shall=0.15, tun_rad_deep =0.3
  integer,parameter :: coupl_rad_cupar=0 ! 0 -no, 1-yes
  
  ! set defaults for rhu00
  rhu00(:) = 2.0
  ! define rh perturbation in order to estimate rhdfda
  rhpert = 0.01 
  
  !- initialization of cuparm parameters
  if(nnqparm(ngrid) == 5 .or. nnqparm(ngrid) == 6) then
    upmf      (     1:m2,1:m3)= xmb4d(     1:m2,1:m3,1,ngrid) !- mass flux deep    convection
    upmfsh    (     1:m2,1:m3)= xmb4d(     1:m2,1:m3,2,ngrid) !- mass flux shallow convection
    zup       (1:m1,1:m2,1:m3)= zup5d(1:m1,1:m2,1:m3,1,ngrid) !- normalized mass flux
    zupshallow(1:m1,1:m2,1:m3)= zup5d(1:m1,1:m2,1:m3,2,ngrid) !- normalized mass flux
    if(coupl_rad_cupar == 1 ) then
      clwup     (1:m1,1:m2,1:m3)= tun_rad_deep *clwup5d(1:m1,1:m2,1:m3,1,ngrid)
      clwupsh   (1:m1,1:m2,1:m3)= tun_rad_shall*clwup5d(1:m1,1:m2,1:m3,2,ngrid)
    endif
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
    do k=1,m1-1
          upmfshx=zupshallow(k+1,i,j) * upmfsh(i,j)
          !-orig shallowcu = max(0.0,min(sh1*log(1.0+sh2*cmfmc2(i,k+1)),0.30))
          shallowcu = max(0.0,min(sh1*log(1.0+sh2*upmfshx),0.30))
          upmfx     = zup(k+1,i,j) * upmf(i,j)
          
          !-orig deepcu = max(0.0,min(dp1*log(1.0+dp2*(cmfmc(i,k+1)-cmfmc2(i,k+1))),0.60))
          ! check possibility of (upmfx-upmfshx) <0. => log (x<0)
          deepcu = max(0.0,min(dp1*log(1.0+dp2*(upmfx-upmfshx)),0.60))
          
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
          elseif(mcphys_type == 2 .or. mcphys_type ==3 .or. mcphys_type ==4 .or. mcphys_type ==6 .or. mcphys_type ==7) then  !srf -gthompson microphysics/gfdl - graupel only in ice phase
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
            lwl(1:m1,ia:iz,ja:jz) = dummy_vec(1:m1,ia:iz,ja:jz)*micro_g(ngrid)%rhp(1:m1,ia:iz,ja:jz)&
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
  
 return
 print*,"=============mcphys-carma radiation coupling ==================="
 print*,"max-min  cl_frac:  ",maxval(cloud_fraction(2:m1,ia:iz,ja:jz)),minval(cloud_fraction(2:m1,ia:iz,ja:jz)) !cloud_fraction
 print*,"max-min ice_frac:  ",maxval(ice_frac(2:m1,ia:iz,ja:jz)),minval(ice_frac(2:m1,ia:iz,ja:jz)) !total cloud ice water 
 print*,"max-min     rain:  ",maxval(rain(ia:iz,ja:jz)),minval(rain(ia:iz,ja:jz)) !total rain water 
 print*,"max-min lwl: ",maxval(lwl(2:m1,ia:iz,ja:jz)),minval(lwl(2:m1,ia:iz,ja:jz)) !total cloud liquid water (kg/kg for carma and g/m2 for rrtm)
 print*,"max-min iwl: ",maxval(iwl(2:m1,ia:iz,ja:jz)),minval(iwl(2:m1,ia:iz,ja:jz))     !total cloud ice water (kg/kg for carma and g/m2 for rrtm)
 print*,"max-min micro_g(ngrid)%rel: ",maxval(micro_g(ngrid)%rel(2:m1,ia:iz,ja:jz)),minval(micro_g(ngrid)%rel(2:m1,ia:iz,ja:jz))  !total cloud liquid water 
 print*,"max-min rei: ",maxval(micro_g(ngrid)%rei(2:m1,ia:iz,ja:jz)),minval(micro_g(ngrid)%rei(2:m1,ia:iz,ja:jz))  !total cloud ice water 

 !print*,"=============================================="
 call flush(6)
    
 end subroutine cloud_prop_carma

END MODULE carma_driv
