!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine datassim()

  use mem_tend, only: &
       tend

  use mem_basic, only: &
       basic_g

  use mem_grid, only: &
       jdim, ngrid, time

  use mem_varinit, only: &
       nud_cond, tcond_beg, tcond_end, varinit_g

  use node_mod, only: &
       ia, iz, ja, jz, mxp, myp, mzp, ibcon

!--(DMK-CCATT-INI)-----------------------------------------------------
  use chem1_list, only: &
       nspecies, fdda, on, spc_alloc

  use mem_chem1, only: &
       chem_assim,     &
       chem1_g
!--(DMK-CCATT-END)-----------------------------------------------------

  implicit none

  integer :: il,ir,jl,jr

!--(DMK-CCATT-INI)-----------------------------------------------------
  integer :: nspc
!--(DMK-CCATT-END)-----------------------------------------------------

  !     Set bounds for nudging this sub-domain

  il=ia  ;  ir=iz  ;  jl=ja  ;  jr=jz

  !     West and east boundaries.
  if(iand(ibcon,1) /= 0)  il=1
  if(iand(ibcon,2) /= 0)  ir=mxp


  !     South and north boundaries.
  if(jdim == 1) then
     if(iand(ibcon,4) /= 0) jl=1
     if(iand(ibcon,8) /= 0) jr=myp
  endif


  ! Basic boundary and analysis nudging scheme

  call nudge(mzp,mxp,myp,il,ir,jl,jr,varinit_g(ngrid)%varwts(1,1,1)  &

       ,varinit_g(ngrid)%varup(1,1,1),varinit_g(ngrid)%varvp(1,1,1)  &
       ,varinit_g(ngrid)%varpp(1,1,1),varinit_g(ngrid)%vartp(1,1,1)  &
       ,varinit_g(ngrid)%varrp(1,1,1) &

       ,varinit_g(ngrid)%varuf(1,1,1),varinit_g(ngrid)%varvf(1,1,1)  &
       ,varinit_g(ngrid)%varpf(1,1,1),varinit_g(ngrid)%vartf(1,1,1)  &
       ,varinit_g(ngrid)%varrf(1,1,1)  &

       ,basic_g(ngrid)%up(1,1,1)   ,basic_g(ngrid)%vp(1,1,1)  &
       ,basic_g(ngrid)%theta(1,1,1),basic_g(ngrid)%rtp(1,1,1)  &
       ,basic_g(ngrid)%pp(1,1,1)  &
       ,tend%ut(1),tend%vt(1),tend%tht(1),tend%rtt(1),tend%pt(1))

!--(DMK-CCATT-INI)-----------------------------------------------------
if (chem_assim == 1) then
   do nspc=1,nspecies
     if(spc_alloc(fdda,nspc) == on) then
       call chem_nudge(mzp,mxp,myp,il,ir,jl,jr, nspc  &
                  ,varinit_g(ngrid)%varwts_chem(1,1,1)  &
                  ,chem1_g(nspc,ngrid)%sc_pp(1,1,1) & !past
                  ,chem1_g(nspc,ngrid)%sc_pf(1,1,1) & !future
                  ,chem1_g(nspc,ngrid)%sc_p (1,1,1) & !current
 		  ,chem1_g(nspc,ngrid)%sc_t (1)     ) !tendency
     endif
   enddo
endif
!--(DMK-CCATT-END)-----------------------------------------------------

  ! Condensate nudging scheme

  if (nud_cond == 1 .and. time >= tcond_beg .and. time <= tcond_end) &
       call nudge_cond(mzp,mxp,myp,il,ir,jl,jr,varinit_g(ngrid)%varwts(1,1,1)  &
       ,varinit_g(ngrid)%varrph(1,1,1),varinit_g(ngrid)%varcph(1,1,1) &
       ,varinit_g(ngrid)%varrfh(1,1,1),varinit_g(ngrid)%varcfh(1,1,1)  &
       ,basic_g(ngrid)%rtp(1,1,1),tend%rtt(1))

  return
end subroutine datassim

!     ******************************************************************

subroutine nudge(m1,m2,m3,ia,iz,ja,jz,varwts  &
     ,varup,varvp,varpp,vartp,varrp  &
     ,varuf,varvf,varpf,vartf,varrf  &
     ,up,vp,theta,rtp,pp,ut,vt,tht,rtt,pt)

  use mem_grid, only: &
       time, ngrid, jdim, grid_g, nnzp, nnxp,nnyp

  use node_mod, only: &
       i0,            &
       j0,            &
       mchnum,        &
       master_num,    &
       nmachs,        &
       mxp,           &
       myp,           &
       mzp,           &
       nodei0,        &
       nodej0,        &
       nodemxp,       &
       nodemyp,       &
       mynum

  use mem_varinit, only: &
       nud_type, vtime1, vtime2, htime1, htime2, &
       wt_nudge_grid, wt_nudge_uv, wt_nudge_th, &
       wt_nudge_pi, wt_nudge_rt, timeWindowIAU, varinit_g

  use mem_scratch, only: &
       vctr1, vctr2, vctr3, vctr4, vctr5, &
       vctr10, vctr11, vctr12, vctr13, vctr14
       
  use mem_basic

  use ModEvaluation, only: &
    comunicateStatistic, &
    evaluate

  implicit none

  integer :: m1,m2,m3,ia,iz,ja,jz
  real, dimension(m1,m2,m3) :: varup,varvp,vartp,varrp,varpp  &
       ,varuf,varvf,vartf,varrf,varpf  &
       ,varwts,up,vp,theta,rtp,pp  &
       ,ut,vt,tht,rtt,pt

  integer :: i,j,k,iCount,v
  real :: tfact, wt_uv, wt_th, wt_pi, wt_rt
  real :: rmse(m1,5),bias(m1,5),soma(m1,5),somaQ(m1,5),centFact
  
  !srf - special weights for pressure and/or UV (only for operations) 
  real, dimension(m1,m2,m3) :: varwts_for_operations_only
  if(  wt_nudge_pi < 0. .or. wt_nudge_uv < 0. ) varwts_for_operations_only=0.

  !         Linearly interpolate values in time, then nudge.

  !-- initialize and update (if desired) the nudging weights (lat, top, center)
  call VariableWeight(nnzp(1), nodemxp(mynum,1), nodemyp(mynum,1), nnxp(1),&
                      nnyp(1), nodei0 (mynum,1), nodej0 (mynum,1),  &
                      grid_g(1)%topt(1,1), grid_g(1)%rtgt(1,1), varinit_g(1)%varwts(1,1,1),&
                      varwts_for_operations_only)

  if (nud_type == 2 .or. nud_type == 4) then
     tfact=(time-vtime1)/(vtime2-vtime1)
  elseif (nud_type == 1) then
     tfact=(time-htime1)/(htime2-htime1)
  endif

  wt_uv=wt_nudge_grid(ngrid)* wt_nudge_uv
  wt_th=wt_nudge_grid(ngrid)* wt_nudge_th
  wt_pi=wt_nudge_grid(ngrid)* wt_nudge_pi
  wt_rt=wt_nudge_grid(ngrid)* wt_nudge_rt
 
  if(evaluate==1) then
    soma=0
    somaQ=0
    !iCount=0
  endif
  do j=ja,jz
     do i=ia,iz
        !iCount=iCount+1
        do k=1,m1
           vctr1(k)=varup(k,i,j)+(varuf(k,i,j)-varup(k,i,j))*tfact
           vctr2(k)=varvp(k,i,j)+(varvf(k,i,j)-varvp(k,i,j))*tfact
           vctr3(k)=vartp(k,i,j)+(vartf(k,i,j)-vartp(k,i,j))*tfact
           vctr4(k)=varrp(k,i,j)+(varrf(k,i,j)-varrp(k,i,j))*tfact
           vctr5(k)=varpp(k,i,j)+(varpf(k,i,j)-varpp(k,i,j))*tfact

           if(evaluate==1) then
            !Sum the quadratic error to statistic RMSE
            somaQ(k,1)=somaQ(k,1)+(basic_g(1)%up(k,i,j)-vctr1(k))**2
            somaQ(k,2)=somaQ(k,2)+(basic_g(1)%vp(k,i,j)-vctr2(k))**2
            somaQ(k,3)=somaQ(k,3)+(basic_g(1)%theta(k,i,j)-vctr3(k))**2
            somaQ(k,4)=somaQ(k,4)+(basic_g(1)%rtp(k,i,j)-vctr4(k))**2
            somaQ(k,5)=somaQ(k,5)+(basic_g(1)%pp(k,i,j)-vctr5(k))**2

            !Sum the absolute error to statistic
            soma(k,1)=soma(k,1)+(basic_g(1)%up(k,i,j)-vctr1(k))
            soma(k,2)=soma(k,2)+(basic_g(1)%vp(k,i,j)-vctr2(k))
            soma(k,3)=soma(k,3)+(basic_g(1)%theta(k,i,j)-vctr3(k))
            soma(k,4)=soma(k,4)+(basic_g(1)%rtp(k,i,j)-vctr4(k))
            soma(k,5)=soma(k,5)+(basic_g(1)%pp(k,i,j)-vctr5(k))
           endif

           vctr10(k)=(varwts(k,i,j)+varwts(k,min(m2,i+1),j)   )*.5* wt_uv
           vctr11(k)=(varwts(k,i,j)+varwts(k,i,min(m3,j+jdim)))*.5* wt_uv
           vctr12(k)=varwts(k,i,j)* wt_th
           vctr13(k)=varwts(k,i,j)* wt_pi
           vctr14(k)=varwts(k,i,j)* wt_rt
        enddo

!-srf special weights for pressure (full domain)
        if(wt_nudge_pi < 0.) then   
           do k=1,m1
               vctr13(k)=varwts_for_operations_only(k,i,j)* abs (wt_pi)
           enddo
        endif
        if(wt_nudge_uv < 0.) then   
           do k=1,m1              
               vctr10(k)=(varwts_for_operations_only(k,i,j)+&
                          varwts_for_operations_only(k,min(m2,i+1),j)   )*.5* abs(wt_uv)
               vctr11(k)=(varwts_for_operations_only(k,i,j)+&
                          varwts_for_operations_only(k,i,min(m3,j+jdim)))*.5* abs(wt_uv)
           enddo
        endif
!-srf end

        do k=1,m1
           ut (k,i,j) = ut(k,i,j) + vctr10(k)*(vctr1(k)-up   (k,i,j))
           vt (k,i,j) = vt(k,i,j) + vctr11(k)*(vctr2(k)-vp   (k,i,j))
           tht(k,i,j) =tht(k,i,j) + vctr12(k)*(vctr3(k)-theta(k,i,j))
           pt (k,i,j) = pt(k,i,j) + vctr13(k)*(vctr5(k)-pp   (k,i,j))
           rtt(k,i,j) =rtt(k,i,j) + vctr14(k)*(vctr4(k)-rtp  (k,i,j))
        enddo

     enddo
  enddo

  !Comunicate the local errors to master and let the master compute RMSE and BIAS
  if(evaluate==1)  call comunicateStatistic(soma,somaQ)

  return
end subroutine nudge

!     ******************************************************************

subroutine VariableWeight(mzp, mxp, myp, nxp, nyp, i0, j0, &
     topt, rtgt, varwts,varwts_for_operations_only)

  use mem_grid, only: &
       ztop, zt, time,dtlt

  use mem_varinit, only: &
       nudlat, tnudcent, tnudtop, tnudlat, znudtop, timeWindowIAU, ramp
  
  use modIau, only:   applyIAU
  
  use node_mod, only: mynum,master_num,mchnum

  implicit none
  integer, intent(in)  :: mzp
  integer, intent(in)  :: mxp
  integer, intent(in)  :: myp
  integer, intent(in)  :: nxp
  integer, intent(in)  :: nyp
  integer, intent(in)  :: i0
  integer, intent(in)  :: j0
  real,    intent(in)    :: topt(mxp,myp)
  real,    intent(in)    :: rtgt(mxp,myp)
  real,    intent(inout) :: varwts(mzp,mxp,myp)
  !srf - special weights for pressure (only for operations) 
  real,    intent(inout) :: varwts_for_operations_only(mzp,mxp,myp)

  integer :: i,j,k
  integer :: iGlobal
  integer :: jGlobal
  real :: tnudcenti,tnudtopi,tnudlati
  real :: tnudtopi_x,tnudlati_x
  real :: rown,rows,rowe,roww
  real :: zloc,wttop,wtlat,delzi,centFact,delPBL,wtbot
  character(len=*), parameter :: h="**(VariableWeight)**"
  real, parameter :: znudbot=1000.  ! meters above local terrain

  !         Get weights for large scale and model tendencies
  if (nudlat .le. 0) return
  
  ! 1) time = 0. do it anyway (needs for initialization)
  if(applyIAU > 0)  then 
     if((timeWindowIAU==0.0 .and. time>0.0) .or. (timeWindowIAU>0.01 .and. time>timeWindowIAU+RAMP)) return 
  endif

  !--- t     => t+ iau => nudging at center with centFact =1 
  !--- t+iau => t+ iau + ramp => nudging decays to zero in this interval with centFact 1 -> 0 


  if(timeWindowIAU>0.01 .and. time <= timeWindowIAU+RAMP .and. RAMP > 0.) then
   !centFact=-1*( (time-timeWindowIAU) /(timeWindowIAU+RAMP)) + 1
    centFact=-1*( (time-timeWindowIAU) /(RAMP))     + 1
  else
    centFact=0.0
  endif
  
  !-- this is for the case of requesting nudging at center of model domain 
  !-- for the total time integration
  if(applyIAU == 0)  centFact=1.0
  
  !if(mynum == 1) print*,"ND=",time,timeWindowIAU,timeWindowIAU+RAMP,centFact;call flush(6)

  tnudcenti=0.
  if(tnudcent.gt. .01) tnudcenti=centFact*(1./tnudcent) 
  tnudtopi=0.
  if(tnudtop.gt. .01) then 
      tnudtopi  =1./tnudtop-tnudcenti
      tnudtopi_x=1./tnudtop
  endif
  tnudlati=0.
  if(tnudlat.gt. .01) then 
    tnudlati  =1./tnudlat-tnudcenti
    tnudlati_x=1./tnudlat
  endif
  
  if(ztop.gt.znudtop) then
     delzi=1./(ztop-znudtop)
  elseif(tnudtop.gt. .01) then
     print*,'Incorrect specification of znudtop ! ! !'
     print*,' znudtop = ',znudtop
     print*,'    ztop = ',ztop
     stop 'varwt-znud'
  endif
  
  if(applyIAU > 0)  then 
    if(time <= dtlt + timeWindowIAU .and. tnudcent > .01 .and. mchnum == master_num) then
      !print*,"=> doing nudging in center:",time,centFact,tnudcenti
       write(*,100) "===> doing nudging in center domain              :",time,centFact,tnudcenti
       call flush(6)
    endif
    if(time >  dtlt + timeWindowIAU .and. tnudcent > .01 .and. mchnum == master_num) then
      !print*,"=> ramping the nudging in center to zero:",time,centFact,tnudcenti
       write(*,100) "===> ramping the nudging in center domain to zero:",time,centFact,tnudcenti
       call flush(6)
    endif
  endif

  !-- turn off nudging from surface to znudbot
  delPBL= 0.
  if(znudbot>0.01) delPBL= 1./znudbot

  do j=1,myp
     jGlobal = j + j0
     do i=1,mxp
        iGlobal = i + i0

        ! quadratic weight function for lateral boundaries

        rown=max(0.,float(jGlobal+nudlat-nyp))
        rows=max(0.,float(nudlat+1-jGlobal))
        rowe=max(0.,float(iGlobal+nudlat-nxp))
        roww=max(0.,float(nudlat+1-iGlobal))
        wtlat=max(rown*rown,rows*rows,rowe*rowe,roww*roww)  &
             /float(nudlat*nudlat)

        ! linear weight function for top boundary

        do k=1,mzp
           zloc=zt(k)*rtgt(i,j)+topt(i,j)
           wttop=max(0.,(zloc-znudtop)*delzi)
	        !-- turn off nudging from surface to znudbot
	        if(znudbot>0.01) then 
             wtbot =min(1.,max(0.,(zloc-znudbot+topt(i,j))*delPBL))
           else
             wtbot =1.
	        endif

! if( mchnum == master_num) then
!    if(i==10 .and. j==10) print*,'wtbot=',k,wtbot,zloc-topt(i,j)
!    call flush(6)
! endif

           ! full 3-D weight function

!-srf orig varwts(k,i,j)=tnudcenti       + max(tnudlati*wtlat,tnudtopi*wttop)
!-srf excluding nudging in the PBL (wtbot)
           varwts(k,i,j)=tnudcenti*wtbot + max(tnudlati*wtlat,tnudtopi*wttop)

!-srf special weights for pressure (pressure will be nudged in the full domain)          
           varwts_for_operations_only(k,i,j) = tnudlati_x*wtbot + max(tnudlati_x*wtlat,tnudtopi_x*wttop)          
       end do
     end do
  end do

!PRINT*,"MX",maxval(varwts),maxval(varwts_for_operations_only)
!PRINT*,"MN",minval(varwts),minval(varwts_for_operations_only)

100 format(1x,A50,F10.1,F8.2,E12.3)

end subroutine VariableWeight



!subroutine evaluateVars(ia,iz,ja,jz,m1,uWind,vWind,rv,theta,)
!,basic_g(ngrid)%up(1,1,1)   ,basic_g(ngrid)%vp(1,1,1)  &
!,basic_g(ngrid)%theta(1,1,1),basic_g(ngrid)%rtp(1,1,1)  &
!,basic_g(ngrid)%pp(1,1,1)  &
!,tend%ut(1),tend%vt(1),tend%tht(1),tend%rtt(1),tend%pt(1))

!,up,vp,theta,rtp,pp,ut,vt,tht,rtt,pt)

!  use mem_grid, only: &
!       time, ngrid, jdim

!  use mem_varinit, only: &
!       nud_type, vtime1, vtime2, htime1, htime2, &
!       wt_nudge_grid, wt_nudge_uv, wt_nudge_th, &
!       wt_nudge_pi, wt_nudge_rt

!  ! Linearly interpolate values in time

!  if (nud_type == 2 .or. nud_type == 4) then
!     tfact=(time-vtime1)/(vtime2-vtime1)
!  elseif (nud_type == 1) then
!     tfact=(time-htime1)/(htime2-htime1)
!  endif

!  do j=ja,jz
!     do i=ia,iz

!        do k=1,m1
!           vctr1(k)=varup(k,i,j)+(varuf(k,i,j)-varup(k,i,j))*tfact !U
!           vctr2(k)=varvp(k,i,j)+(varvf(k,i,j)-varvp(k,i,j))*tfact !V
!           vctr3(k)=vartp(k,i,j)+(vartf(k,i,j)-vartp(k,i,j))*tfact !Theta
!           vctr4(k)=varrp(k,i,j)+(varrf(k,i,j)-varrp(k,i,j))*tfact !rtp
!           vctr5(k)=varpp(k,i,j)+(varpf(k,i,j)-varpp(k,i,j))*tfact !pp

!           vctr10(k)=(varwts(k,i,j)+varwts(k,min(m2,i+1),j))*.5* wt_uv
!           vctr11(k)=(varwts(k,i,j)+varwts(k,i,min(m3,j+jdim)))*.5* wt_uv
!           vctr12(k)=varwts(k,i,j)* wt_th
!           vctr13(k)=varwts(k,i,j)* wt_pi
!           vctr14(k)=varwts(k,i,j)* wt_rt
!        enddo

!        do k=1,m1
!           ut(k,i,j) = ut(k,i,j) + vctr10(k)*(vctr1(k)-up(k,i,j))
!           vt(k,i,j) = vt(k,i,j) + vctr11(k)*(vctr2(k)-vp(k,i,j))
!           tht(k,i,j)=tht(k,i,j) + vctr12(k)*(vctr3(k)-theta(k,i,j))
!           pt(k,i,j) = pt(k,i,j) + vctr13(k)*(vctr5(k)-pp(k,i,j))
!           rtt(k,i,j)=rtt(k,i,j) + vctr14(k)*(vctr4(k)-rtp(k,i,j))
!        enddo

!     enddo
!  enddo
  
! end subroutine evaluateVars



subroutine nudge_cond(m1,m2,m3,ia,iz,ja,jz,varwts  &
     ,varrph,varcph,varrfh,varcfh,rtp,rtt)

  use mem_grid, only: &
       time, ngrid

  use mem_varinit, only: &
       condtime1, condtime2, wt_nudgec_grid, t_nudge_rc

  use mem_scratch, only: &
       vctr1, vctr2, vctr3

  implicit none

  integer :: m1,m2,m3,ia,iz,ja,jz
  real, dimension(m1,m2,m3) :: varrph,varcph,varrfh,varcfh  &
       ,varwts,rtp,rtt

  integer :: i,j,k
  real :: tfact, wt_rc



  ! tfact is temporal interpolation weight,
  ! wt_rc is timescale and grid-dependent weight

  tfact=(time-condtime1)/(condtime2-condtime1)   &
       * wt_nudgec_grid(ngrid)/t_nudge_rc

  wt_rc=wt_nudgec_grid(ngrid)/t_nudge_rc

  do j=ja,jz
     do i=ia,iz

        do k=1,m1
           vctr1(k)=varrph(k,i,j)+(varrfh(k,i,j)-varrph(k,i,j))*tfact
           vctr2(k)=varcph(k,i,j)+(varcfh(k,i,j)-varcph(k,i,j))*tfact

           vctr3(k)=varwts(k,i,j)*wt_rc
        enddo

        ! Only nudging total water where condensate exists...

        do k=1,m1
           if (vctr2(k) > 0.)  &
                rtt(k,i,j)=rtt(k,i,j) + vctr3(k)*(vctr1(k)-rtp(k,i,j))
        enddo

     enddo
  enddo

  return
end subroutine nudge_cond

!     ****************************************************************

subroutine varweight(n1,n2,n3,varwts,topt,rtgt)

  use mem_grid, only: &
       ztop, nzp, zt, &
!--(DMK-CCATT-INI)-----------------------------------------------------
       ztn
!--(DMK-CCATT-FIM)-----------------------------------------------------

  use mem_varinit, only: &
       nudlat, tnudcent, tnudtop, tnudlat, znudtop

  implicit none
  integer :: n1,n2,n3
  real :: varwts(n1,n2,n3),topt(n2,n3),rtgt(n2,n3)

  integer :: i,j,k
  real :: tnudcenti,tnudtopi,tnudlati,rown,rows,rowe,roww,zloc,wttop &
       ,wtlat,delzi
  
  !         Get weights for large scale and model tendencies
  if (nudlat .le. 0) return

  tnudcenti=0.
  if(tnudcent.gt. .01) tnudcenti=1./tnudcent
  tnudtopi=0.
  if(tnudtop.gt. .01) tnudtopi=1./tnudtop-tnudcenti
  tnudlati=0.
  if(tnudlat.gt. .01) tnudlati=1./tnudlat-tnudcenti

  if(ztop.gt.znudtop) then
     delzi=1./(ztop-znudtop)
  elseif(tnudtop.gt. .01) then
     print*,'Incorrect specification of znudtop ! ! !'
     print*,' znudtop = ',znudtop
     print*,'    ztop = ',ztop
     stop 'varwt-znud'
  endif


  do j=1,n3
     do i=1,n2

        !                       quadratic weight function for lateral boundaries

        rown=max(0.,float(j+nudlat-n3))
        rows=max(0.,float(nudlat+1-j))
        rowe=max(0.,float(i+nudlat-n2))
        roww=max(0.,float(nudlat+1-i))
        wtlat=max(rown*rown,rows*rows,rowe*rowe,roww*roww)  &
             /float(nudlat*nudlat)

        !                       linear weight function for top boundary

!--(DMK-CCATT-INI)-----------------------------------------------------
      !srf-rams60 mod
      do k=1,n1
         zloc=ztn(k,1)*rtgt(i,j)+topt(i,j)
      !srf-rams60 mod
!--(DMK-original)------------------------------------------------------
!      do k=1,nzp
!         zloc=zt(k)*rtgt(i,j)+topt(i,j)
!--(DMK-CCATT-END)-----------------------------------------------------

           wttop=max(0.,(zloc-znudtop)*delzi)

           !                       full 3-D weight function

           varwts(k,i,j)=tnudcenti  &
                +max(tnudlati*wtlat,tnudtopi*wttop)

        enddo

     enddo
  enddo

  return
end subroutine varweight


!     *****************************************************************

subroutine vfintrpf(ifm,ifflag)

  use mem_grid, only: &
       nxtnest, nnxp, nnyp, maxnxp, maxnyp, grid_g

  use mem_scratch, only: &
       scratch

  use mem_varinit, only: &
       varinit_g

  use mem_basic, only: &
       basic_g

  implicit none
  real :: scr1(maxnxp*maxnyp)
  real :: scr2(maxnxp*maxnyp)
  integer :: ifm,ifflag,icm

  icm = nxtnest(ifm)
  if (icm == 0) return

  !    Temporarily fill VT2DA with interpolated topography from coarser grid

  call fillscr(1,maxnxp,maxnyp,1,nnxp(icm),nnyp(icm),1,1  &
       ,scr1(1),grid_g(icm)%topt(1,1))
  call eintp(scr1(1),scr2(1)  &
       ,1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
  call fillvar(1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1  &
       ,scr2(1),scratch%vt2da(1))

  if (ifflag == 1) then

     !     Interpolate varwts

     call fmint4(varinit_g(icm)%varwts(1,1,1)  &
          ,varinit_g(ifm)%varwts(1,1,1)  &
          ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
          ,scratch%vt2da(1),ifm,icm,'t',0)

  endif


  if (ifflag == 2) then

     !     Interpolate future level atmospheric variables

     call fmint4(varinit_g(icm)%varuf(1,1,1),varinit_g(ifm)%varuf(1,1,1)  &
          ,basic_g(icm)%dn0u(1,1,1),basic_g(ifm)%dn0u(1,1,1)  &
          ,scratch%vt2da(1),ifm,icm,'u',1)
     call fmint4(varinit_g(icm)%varvf(1,1,1),varinit_g(ifm)%varvf(1,1,1)  &
          ,basic_g(icm)%dn0v(1,1,1),basic_g(ifm)%dn0v(1,1,1)  &
          ,scratch%vt2da(1),ifm,icm,'v',1)
     call fmint4(varinit_g(icm)%varpf(1,1,1),varinit_g(ifm)%varpf(1,1,1)  &
          ,basic_g(icm)%dn0v(1,1,1),basic_g(ifm)%dn0v(1,1,1)  &
          ,scratch%vt2da(1),ifm,icm,'t',1)
     call fmint4(varinit_g(icm)%vartf(1,1,1),varinit_g(ifm)%vartf(1,1,1)  &
          ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
          ,scratch%vt2da(1),ifm,icm,'t',1)
     call fmint4(varinit_g(icm)%varrf(1,1,1),varinit_g(ifm)%varrf(1,1,1)  &
          ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
          ,scratch%vt2da(1),ifm,icm,'t',1)

  endif

  return
end subroutine vfintrpf




subroutine VarfIntrp(ifm,ifflag)

  use mem_grid, only: &
       nxtnest, nnxp, nnyp, maxnxp, maxnyp, grid_g

  use mem_scratch, only: &
       scratch

  use mem_varinit, only: &
       varinit_g

  use mem_basic, only: &
       basic_g

!--(DMK-CCATT-INI)-----------------------------------------------------
  use chem1_list, only: &
       nspecies, fdda, spc_alloc

  use mem_chem1, only: &
       chem1_g,        &
       chem_assim
!--(DMK-CCATT-FIM)-----------------------------------------------------

  use dump,only: &
      dumpMessage

  implicit none
  include "constants.f90"
  integer :: ifm,ifflag,icm
  real :: scr1(maxnxp*maxnyp)
  real :: scr2(maxnxp*maxnyp)
  character(len=*), parameter :: h="**(VarfIntrp)**"

!--(DMK-CCATT-INI)-----------------------------------------------------
  integer :: nspc
!--(DMK-CCATT-FIM)-----------------------------------------------------

  icm = nxtnest(ifm)
  if (icm == 0) return

  !    Temporarily fill VT2DA with interpolated topography from coarser grid
  !**(JP)** not converted

  !call fatal_error(h//"**(JP)** not converted")
  iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
    "**(JP)** not converted")

  call fillscr(1,maxnxp,maxnyp,1,nnxp(icm),nnyp(icm),1,1  &
       ,scr1(1),grid_g(icm)%topt(1,1))
  call eintp(scr1(1),scr2(1)  &
       ,1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
  call fillvar(1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1  &
       ,scr2(1),scratch%vt2da(1))

  if (ifflag == 1) then

     !     Interpolate varwts

     call fmint4(varinit_g(icm)%varwts(1,1,1)  &
          ,varinit_g(ifm)%varwts(1,1,1)  &
          ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
          ,scratch%vt2da(1),ifm,icm,'t',0)

!--(DMK-CCATT-INI)-----------------------------------------------------
!     Interpolate varwts_chem
      if(chem_assim == 1 ) call fmint4(varinit_g(icm)%varwts_chem(1,1,1)  &
                                       ,varinit_g(ifm)%varwts_chem(1,1,1)  &
                                       ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
                                       ,scratch%vt2da(1),ifm,icm,'t',0)
!--(DMK-CCATT-END)-----------------------------------------------------

  else if (ifflag == 2) then

     !     Interpolate future level atmospheric variables

     call fmint4(varinit_g(icm)%varuf(1,1,1),varinit_g(ifm)%varuf(1,1,1)  &
          ,basic_g(icm)%dn0u(1,1,1),basic_g(ifm)%dn0u(1,1,1)  &
          ,scratch%vt2da(1),ifm,icm,'u',1)
     call fmint4(varinit_g(icm)%varvf(1,1,1),varinit_g(ifm)%varvf(1,1,1)  &
          ,basic_g(icm)%dn0v(1,1,1),basic_g(ifm)%dn0v(1,1,1)  &
          ,scratch%vt2da(1),ifm,icm,'v',1)
     call fmint4(varinit_g(icm)%varpf(1,1,1),varinit_g(ifm)%varpf(1,1,1)  &
          ,basic_g(icm)%dn0v(1,1,1),basic_g(ifm)%dn0v(1,1,1)  &
          ,scratch%vt2da(1),ifm,icm,'t',1)
     call fmint4(varinit_g(icm)%vartf(1,1,1),varinit_g(ifm)%vartf(1,1,1)  &
          ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
          ,scratch%vt2da(1),ifm,icm,'t',1)
     call fmint4(varinit_g(icm)%varrf(1,1,1),varinit_g(ifm)%varrf(1,1,1)  &
          ,basic_g(icm)%dn0(1,1,1),basic_g(ifm)%dn0(1,1,1)  &
          ,scratch%vt2da(1),ifm,icm,'t',1)

!--(DMK-CCATT-INI)-----------------------------------------------------
!     Interpolate future level chemical variables
     if (chem_assim == 1) then
        do nspc=1,nspecies
           if(spc_alloc(fdda,nspc) == 1) then
              call fmint4(chem1_g(nspc,icm)%sc_pf(1,1,1) &! varinit_g(icm)%varrf(1,1,1)
       	                  ,chem1_g(nspc,ifm)%sc_pf(1,1,1) &! varinit_g(ifm)%varrf(1,1,1)  &
                          ,basic_g(icm)%dn0(1,1,1)  &
                          ,basic_g(ifm)%dn0(1,1,1)  &
                          ,scratch%vt2da(1),ifm,icm,'t',1)
           endif
        enddo
     endif
!--(DMK-CCATT-END)-----------------------------------------------------

  endif
end subroutine VarfIntrp

!--(DMK-CCATT-INI)-----------------------------------------------------------
subroutine chem_nudge(m1,m2,m3,ia,iz,ja,jz,nspc &
                     ,varwts_chem  &
                     ,sc_pp  &
                     ,sc_pf  &
                     ,scp    &
                     ,sct)

use mem_grid, only: &
     time

use mem_varinit, only: &
     vtime1, vtime2,   &
     htime1, htime2,   &
     nud_type

use mem_scratch, only: &
     vctr4, vctr14

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,nspc
real, dimension(m1,m2,m3) :: sc_pp,sc_pf,varwts_chem,scp,sct

integer :: i,j,k
real :: tfact,  wt_chem


!return

!         Linearly interpolate values in time, then nudge.

if (nud_type == 2 .or. nud_type == 4) then
   tfact=(time-vtime1)/(vtime2-vtime1)
elseif (nud_type == 1) then
   tfact=(time-htime1)/(htime2-htime1)
endif

! - specific nudging  weight for specie 'nspc'
!- for now the weight is 1 for all species.
wt_chem = 0.083333 !=1/12 !valor original= 1.
!- in the future, use different weights and pass them to the slaves via chem_mpass_init
!!wt_chem=wt_nudge_grid(ngrid)* wt_nudge_chem(nspc)

do j=ja,jz
   do i=ia,iz

      do k=1,m1
         vctr4(k)=sc_pp(k,i,j)+(sc_pf(k,i,j)-sc_pp(k,i,j))*tfact
         vctr14(k)=varwts_chem(k,i,j)* wt_chem
      enddo

      do k=1,m1
         sct(k,i,j)=sct(k,i,j) + vctr14(k)*(vctr4(k)-scp(k,i,j))
      enddo

   enddo
enddo

return
end subroutine chem_nudge



subroutine VariableWeightChem(mzp, mxp, myp, nxp, nyp, i0, j0, &
     topt, rtgt, varwts)

  use mem_grid, only: &
       ztop, zt, ztn

  use mem_varinit, only: &
       nudlat, tnudcent, tnudtop, tnudlat, znudtop

  implicit none
  integer, intent(in)  :: mzp
  integer, intent(in)  :: mxp
  integer, intent(in)  :: myp
  integer, intent(in)  :: nxp
  integer, intent(in)  :: nyp
  integer, intent(in)  :: i0
  integer, intent(in)  :: j0
  real,    intent(out) :: varwts(mzp,mxp,myp)
  real,    intent(in)  :: topt(mxp,myp)
  real,    intent(in)  :: rtgt(mxp,myp)

  integer :: i,j,k
  integer :: iGlobal
  integer :: jGlobal
  real :: tnudcenti,tnudtopi,tnudlati
  real :: rown,rows,rowe,roww
  real :: zloc,wttop,wtlat,delzi
  character(len=*), parameter :: h="**(VariableWeight)**"

  !         Get weights for large scale and model tendencies

  if (nudlat .le. 0) return

  tnudcenti=0.
  if(tnudcent.gt. .01) tnudcenti=1./tnudcent
  tnudtopi=0.
  if(tnudtop.gt. .01) tnudtopi=1./tnudtop-tnudcenti
  tnudlati=0.
  if(tnudlat.gt. .01) tnudlati=1./tnudlat-tnudcenti

  if(ztop.gt.znudtop) then
     delzi=1./(ztop-znudtop)
  elseif(tnudtop.gt. .01) then
     print*,'Incorrect specification of znudtop ! ! !'
     print*,' znudtop = ',znudtop
     print*,'    ztop = ',ztop
     stop 'varwt-znud'
  endif


  do j=1,myp
     jGlobal = j + j0
     do i=1,mxp
        iGlobal = i + i0

        ! quadratic weight function for lateral boundaries

        rown=max(0.,float(jGlobal+nudlat-nyp))
        rows=max(0.,float(nudlat+1-jGlobal))
        rowe=max(0.,float(iGlobal+nudlat-nxp))
        roww=max(0.,float(nudlat+1-iGlobal))
        wtlat=max(rown*rown,rows*rows,rowe*rowe,roww*roww)  &
             /float(nudlat*nudlat)

        ! linear weight function for top boundary

        do k=1,mzp
           zloc=ztn(k,1)*rtgt(i,j)+topt(i,j)

           wttop=max(0.,(zloc-znudtop)*delzi)

           ! full 3-D weight function

           varwts(k,i,j)=tnudcenti  &
                +max(tnudlati*wtlat,tnudtopi*wttop)
        end do
     end do
  end do
end subroutine VariableWeightChem



subroutine varweight_chem(n1,n2,n3,varwts_chem,topt,rtgt)

use mem_grid
use mem_varinit

implicit none
integer :: n1,n2,n3
real :: varwts_chem(n1,n2,n3),topt(n2,n3),rtgt(n2,n3)

integer :: i,j,k
real :: tnudcenti,tnudtopi,tnudlati,rown,rows,rowe,roww,zloc,wttop &
       ,wtlat,delzi

!         Get weights for large scale and model tendencies

if (nudlat .le. 0) return

tnudcenti=0.
!-srf
!- no chemistry nudging in the inner of domain is allowed
!if(tnudcent.gt. .01) tnudcenti=1./tnudcent

tnudtopi=0.
if(tnudtop.gt. .01) tnudtopi=1./tnudtop-tnudcenti
tnudlati=0.
if(tnudlat.gt. .01) tnudlati=1./tnudlat-tnudcenti



if(ztop.gt.znudtop) then
   delzi=1./(ztop-znudtop)
elseif(tnudtop.gt. .01) then
   print*,'Incorrect specification of znudtop ! ! !'
   print*,' znudtop = ',znudtop
   print*,'    ztop = ',ztop
   stop 'varwt-znud'
endif


do j=1,n3
   do i=1,n2

!                       quadratic weight function for lateral boundaries

      rown=max(0.,float(j+nudlat-n3))
      rows=max(0.,float(nudlat+1-j))
      rowe=max(0.,float(i+nudlat-n2))
      roww=max(0.,float(nudlat+1-i))
      wtlat=max(rown*rown,rows*rows,rowe*rowe,roww*roww)  &
           /float(nudlat*nudlat)

!                       linear weight function for top boundary

      do k=1,nzp

!srf-rams60
         zloc=ztn(k,1)*rtgt(i,j)+topt(i,j)
!         zloc=zt(k)*rtgt(i,j)+topt(i,j)

         wttop=max(0.,(zloc-znudtop)*delzi)

!                       full 3-D weight function

         varwts_chem(k,i,j)=tnudcenti  &
              +max(tnudlati*wtlat,tnudtopi*wttop)

      enddo

   enddo
enddo

return
end
!srf-chem-end
!--(DMK-CCATT-FIM)-----------------------------------------------------------
