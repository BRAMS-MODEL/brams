!------------------------------------------------------------------------------------------!
! Subroutine prep_advflx_to_stilt                                                          !
! Developed by Saulo R. Freitas (CPTEC/INPE)                                               !
!                                                                                          !
!   This subroutine prepares the advective fluxes to be used in STILT (or other            !
! Lagrangian models).                                                                      !
!------------------------------------------------------------------------------------------!
subroutine prep_advflx_to_stilt(mzp,mxp,myp,ia,iz,ja,jz,ng)
   
   use mem_grid   ,  only: dtlt
   use mem_scratch,  only: scratch
   use mem_stilt   ,  only: stilt_g                 & ! structure
                         , frqmassave             & ! intent(in)
                         , etime_adve             & ! intent(inout)
                         , zero_average_mass_adve ! ! subroutine

   implicit none

   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in) :: mzp
   integer, intent(in) :: mxp
   integer, intent(in) :: myp
   integer, intent(in) :: ia
   integer, intent(in) :: iz
   integer, intent(in) :: ja
   integer, intent(in) :: jz
   integer, intent(in) :: ng
   !----- Local variables. ----------------------------------------------------------------!
   real                :: dtlti
   real                :: frqmassi
   !---------------------------------------------------------------------------------------!


   !----- Find the inverse of the time step and the average time. -------------------------!
   dtlti    = 1./dtlt
   frqmassi = 1./frqmassave
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Update the step counter for advective fluxes.  If this exceeds the maximum number !
   ! of steps, it means it is time to reset it and start integrating again.                ! 
   !---------------------------------------------------------------------------------------!
   etime_adve(ng) = etime_adve(ng) + dtlt
   if (etime_adve(ng) > frqmassave + 0.1 * dtlt) then
      etime_adve(ng) = dtlt
      call zero_average_mass_adve(stilt_g(ng))
   end if
   !---------------------------------------------------------------------------------------!

  
   call get_adv_fluxes_for_stilt(mzp,mxp,myp,ia,iz,ja,jz)
 
   call compute_mass_flux(mzp,mxp,myp,ia,iz,ja,jz,dtlti,frqmassi                           &
           ,scratch%vt3da            , scratch%vt3db           , scratch%vt3dc             &
           ,stilt_g(ng)%afxu          , stilt_g(ng)%afxv         , stilt_g(ng)%afxw           &
           ,stilt_g(ng)%afxub         , stilt_g(ng)%afxvb        , stilt_g(ng)%afxwb          )


end subroutine prep_advflx_to_stilt
!------------------------------------------------------------------------------------------!

!==========================================================================================!
!==========================================================================================!
! Subroutine compute_mass_flux                                                             !
! Based on original Saulo R. Freitas (CPTEC/INPE) subroutine                               !
!                                                                                          !
! This subroutine compute the integrated mass flux from advection.                         !
!------------------------------------------------------------------------------------------!
subroutine compute_mass_flux(mzp,mxp,myp,ia,iz,ja,jz,dtlti,frqmassi,vt3da,vt3db,vt3dc,afxu &
                            ,afxv,afxw,afxub,afxvb,afxwb)
   implicit none
   integer                        , intent(in)    :: mzp,mxp,myp
   integer                        , intent(in)    :: ia,iz,ja,jz
   real                           , intent(in)    :: dtlti,frqmassi
   real   , dimension(mzp,mxp,myp), intent(in)    :: vt3da,vt3db,vt3dc
   real   , dimension(mzp,mxp,myp), intent(out)   :: afxu,afxv,afxw
   real   , dimension(mzp,mxp,myp), intent(inout) :: afxub,afxvb,afxwb
   integer                                        :: i,j,k

   do k=1,mzp
      do i=ia,iz
         do j=ja,jz
            afxu(k,i,j)  =                vt3da(k,i,j) * dtlti
            afxv(k,i,j)  =                vt3db(k,i,j) * dtlti
            afxw(k,i,j)  =                vt3dc(k,i,j) * dtlti
            afxub(k,i,j) = afxub(k,i,j) + vt3da(k,i,j) * frqmassi
            afxvb(k,i,j) = afxvb(k,i,j) + vt3db(k,i,j) * frqmassi
            afxwb(k,i,j) = afxwb(k,i,j) + vt3dc(k,i,j) * frqmassi
         end do
      end do
   end do

   return
end subroutine compute_mass_flux
!==========================================================================================!
!==========================================================================================!

!------------------------------------------------------------------------------------------!
! Subroutine prep_convflx_to_stilt                                                         !
! Developed by Saulo R. Freitas (CPTEC/INPE)                                               !
!                                                                                          !
!   This subroutine prepares the convective fluxes to be used in STILT (or other           !
! Lagrangian models).                                                                      !
!------------------------------------------------------------------------------------------!
subroutine prep_convflx_to_stilt(m1,m2,m3,ia,iz,ja,jz,mgmxp,mgmyp,mgmzp,maxiens,ngrid      &
                                ,ngrids_cp,ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d       &
                                ,kpbl4d,kstabi4d,kstabm4d,xmb4d,edt4d,zcup5d,pcup5d,enup5d &
                                ,endn5d,deup5d,dedn5d,zup5d,zdn5d,iens)

use mem_stilt
!-srf including mass fluxes from GF scheme
use mem_scratch1_grell, only:   up_massdetr5d, up_massentr5d,                       &
                                dd_massdetr5d ,dd_massentr5d
use mem_cuparm, only: NNQPARM ! INTENT(IN)
 

implicit none

integer, intent(in) :: mgmxp,mgmyp,mgmzp,ngrid,ngrids_cp,iens, maxiens,m1,m2,m3,ia,iz,ja,jz

integer, intent(in),dimension(mgmxp,mgmyp,maxiens,ngrids_cp) ::                            &
                         ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d,kpbl4d,kstabi4d,kstabm4d
                       
real, intent(in), dimension(mgmxp,mgmyp,maxiens,ngrids_cp) :: xmb4d,edt4d
                        
real, intent(in), dimension(mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp) ::                        &
                                      enup5d,endn5d,deup5d,dedn5d,zup5d,zdn5d,zcup5d,pcup5d

integer :: i, j


!zero out all convflx
if (iens == 1) then
!Deep conv

!--(DMK-BRAMS-5.0)---------------------------------------------
!  call azero(m1*m2*m3,stilt_g(ngrid)%cfxup1(1,1,1))
!  call azero(m1*m2*m3,stilt_g(ngrid)%cfxdn1(1,1,1))
!  call azero(m1*m2*m3,stilt_g(ngrid)%dfxup1(1,1,1))
!  call azero(m1*m2*m3,stilt_g(ngrid)%efxup1(1,1,1))
!  call azero(m1*m2*m3,stilt_g(ngrid)%dfxdn1(1,1,1))
!  call azero(m1*m2*m3,stilt_g(ngrid)%efxdn1(1,1,1))

  stilt_g(ngrid)%cfxup1 = 0.
  stilt_g(ngrid)%cfxdn1 = 0.
  stilt_g(ngrid)%dfxup1 = 0.
  stilt_g(ngrid)%efxup1 = 0.
  stilt_g(ngrid)%dfxdn1 = 0.
  stilt_g(ngrid)%efxdn1 = 0.
!--(DMK-BRAMS-5.0)---------------------------------------------

!shallow conv
elseif (iens == 2) then

!--(DMK-BRAMS-5.0)---------------------------------------------
! call azero(m1*m2*m3,stilt_g(ngrid)%cfxup2(1,1,1))
! call azero(m1*m2*m3,stilt_g(ngrid)%dfxup2(1,1,1))       
! call azero(m1*m2*m3,stilt_g(ngrid)%efxup2(1,1,1))

 stilt_g(ngrid)%cfxup2 = 0.
 stilt_g(ngrid)%dfxup2 = 0.
 stilt_g(ngrid)%efxup2 = 0.
!--(DMK-BRAMS-5.0)---------------------------------------------

end if

do j=ja,jz
  do i=ia,iz
!    if((iens == 1 .and. ierr4d(i,j,iens,ngrid) == 0) .or. iens == 2) then
    if(ierr4d(i,j,iens,ngrid) == 0) then
      call get_convflx(iens,i,j,mgmzp,m1,m2,m3,nnqparm(ngrid)                               &
                  ,   xmb4d(i,j,iens,ngrid),   edt4d(i,j,iens,ngrid)                       &
                  ,  jmin4d(i,j,iens,ngrid),  kdet4d(i,j,iens,ngrid)                       &
                  ,   k224d(i,j,iens,ngrid), kbcon4d(i,j,iens,ngrid)                       &
                  ,  ktop4d(i,j,iens,ngrid),  kpbl4d(i,j,iens,ngrid)                       &
                  ,kstabi4d(i,j,iens,ngrid),kstabm4d(i,j,iens,ngrid)                       &
                  ,zcup5d(1:mgmzp,i,j,iens,ngrid),pcup5d(1:mgmzp,i,j,iens,ngrid)           &
                  ,deup5d(1:mgmzp,i,j,iens,ngrid),enup5d(1:mgmzp,i,j,iens,ngrid)           &
                  ,dedn5d(1:mgmzp,i,j,iens,ngrid),endn5d(1:mgmzp,i,j,iens,ngrid)           &
                  , zup5d(1:mgmzp,i,j,iens,ngrid), zdn5d(1:mgmzp,i,j,iens,ngrid)           &
!-srf: including fluxes from GF scheme
                  ,up_massdetr5d(1:mgmzp,i,j,iens,ngrid)  &
                  ,up_massentr5d(1:mgmzp,i,j,iens,ngrid)  &
                  ,dd_massdetr5d(1:mgmzp,i,j,iens,ngrid)  & 
                  ,dd_massentr5d(1:mgmzp,i,j,iens,ngrid)  &
!-
                  ,stilt_g(ngrid)%cfxup1(1,1,1),stilt_g(ngrid)%cfxdn1(1,1,1)               &
                  ,stilt_g(ngrid)%dfxup1(1,1,1),stilt_g(ngrid)%efxup1(1,1,1)               &
                  ,stilt_g(ngrid)%dfxdn1(1,1,1),stilt_g(ngrid)%efxdn1(1,1,1)               &
                  ,stilt_g(ngrid)%cfxup2(1,1,1),stilt_g(ngrid)%dfxup2(1,1,1)               &
                  ,stilt_g(ngrid)%efxup2(1,1,1))
    endif 
  enddo 
enddo  

return
end subroutine prep_convflx_to_stilt
!------------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------------!
! Subroutine get_convflx                                                                   !
! Developed by Saulo R. Freitas (CPTEC/INPE)                                               !
!                                                                                          !
!   This subroutine aims at getting the convective fluxes from Grell shallow and           !
! deep convective parameterizations.                                                       !
!------------------------------------------------------------------------------------------!
subroutine get_convflx(iens,i,j,mgmzp,m1,m2,m3,nnqparm,xmb,edt,jmin,kdet,k22,kbcon,ktop,kpbl &
                      ,kstabi,kstabm,z_cup,p_cup,cd,entr,cdd,entrd,zu,zd&!,cfxup1,cfxdn1       &
!-srf for GF scheme
                      ,up_massdetro,up_massentro,dd_massdetro,dd_massentro                   &
                      ,cfxup1,cfxdn1,dfxup1,efxup1,dfxdn1,efxdn1,cfxup2,dfxup2,efxup2)

implicit none

integer, intent(in)                      ::  iens,i,j,mgmzp,m1,m2,m3,jmin,kdet,k22,kbcon   &
                                            ,ktop,kpbl,kstabi,kstabm,nnqparm
                                            
real, intent(in)                         ::  xmb,edt

real, intent(in),    dimension(mgmzp)    ::  z_cup,p_cup,cd,entr,cdd,entrd,zu,zd

real, intent(in),    dimension(mgmzp)    ::  up_massdetro,up_massentro,dd_massdetro,dd_massentro   

real, intent(inout), dimension(m1,m2,m3) ::  cfxup1,cfxdn1,dfxup1,efxup1,dfxdn1,efxdn1     &
                                            ,cfxup2,dfxup2,efxup2  
                                             
integer                                  ::  k,kr
real                                     ::  dz,totmas,entup,detup,entdoj,entupk,detupk    &
                                            ,detdo,entdo,subdown,subin,detdo1,detdo2

do k=1,ktop+1
       
  kr= k + 1   ! level K of conv grid  corresponds to level K + 1 of RAMS grid
  dz =  z_cup(kr) - z_cup(k)
   
  if(iens==2 .or. (iens==1 .and. nnqparm == 2)) then
   if(k==1) cycle ! k=1 only for GF scheme
   entup   = 0.
   detup   = 0.
   entdoj  = 0.
   entupk  = 0.
   detupk  = 0.
   detdo   = edt*cdd(k)*dz*zd(kr)
   entdo   = edt*entrd(k)*dz*zd(kr)
   subdown = zu(k ) - edt*zd(k )
   subin   = zu(kr) - edt*zd(kr)

   if (k >= kbcon .and. k < ktop) then
     entup  = entr(k)   *dz*zu(k)
     detup  =   cd(kr) *dz*zu(k)
   end if

   if(k == jmin)  entdoj = edt*zd(k)
   if(k == k22-1) entupk = zu(kpbl)
   if(k == ktop)  detupk = zu(ktop)
   if(k > kdet)   detdo  = 0.
   if(k == ktop)  subin  = 0.
   if(k < kbcon)  detup  = 0.
  endif

  if(iens==1 .and. nnqparm == 5) then ! GF scheme
      entdoj  = 0. 
      entupk  = 0. 
      detupk  = 0. 
     
      subdown = zu(k  ) - edt*zd(k  ) 
      subin   = zu(k+1) - edt*zd(k+1) 
   ! detrainment and entrainment for downdrafts
      detdo=edt*dd_massdetro(k)
      entdo=edt*dd_massentro(k)
   !
   ! entrainment/detrainment for updraft
      entup=up_massentro(k)
      detup=up_massdetro(k)
      
      if(k==jmin)  entdoj=edt*zd(k)

      if(k == ktop)then
               detupk=zu(ktop)
               detdo=0.
               entdo=0.
               entup=0.
               detup=0.
       endif 
  endif


  if(iens == 1) then ! Deep convection
      cfxup1(kr,i,j) =     xmb* zu(k)
      cfxdn1(kr,i,j) =-edt*xmb* zd(k)
      dfxup1(kr,i,j) =     xmb*(detup + detupk)
      efxup1(kr,i,j) =    -xmb*(entup + entupk)
      dfxdn1(kr,i,j) =     xmb*(detdo         ) !edt already is at detdo
      efxdn1(kr,i,j) =    -xmb*(entdo + entdoj) !edt already is at entdo,entdoj
  elseif(iens == 2)  then ! Shallow convection
      cfxup2(kr,i,j) =     xmb* zu(k)
      dfxup2(kr,i,j) =     xmb*(detup + detupk)
      efxup2(kr,i,j) =    -xmb*(entup + entupk)
  end if
!------------------------------------------------------------------------------------------!
! Checking the mass conservation                                                           !
!------------------------------------------------------------------------------------------!
  totmas=subin-subdown+detup-entup-entdo+detdo-entupk-entdoj+detupk
  if(abs(totmas) > 1.e-6) then
    write (unit=*,fmt='(a)')                 '----------- Subroutine Get_convflx ----------'
    write(unit=*, fmt='(4(a,1x,i3,1x))')     '  K= ',k,'   I=',i,'   J=',j,'   IENS=',iens
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  subin=  ',    subin,'subdown=',subdown
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  detup=  ',    detup,'entup=  ',entup
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  entdo=  ',    entdo,'detdo=  ',detdo
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  entupk= ',   entupk,'detupk= ',detupk
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  zu(k)=  ',    zu(k),'zd(k)=  ',zd(k)
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  zu(kr)= ',   zu(kr),'zd(kr)= ',zd(kr)
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  entdoj= ',   entdoj,'edt=    ',edt
    write(unit=*, fmt='(1(a,1x,es10.3,1x))') '  totmas= ',   totmas
    write(unit=*, fmt='(a)')                 '---------------------------------------------'
    stop 'The model will stop since it is not conserving mass...' 
  end if

end do

!------------------------------------------------------------------------------------------!
! Bottom layer                                                                             !
!------------------------------------------------------------------------------------------!
if(iens == 1 .and. nnqparm == 2) then  ! Deep convection
  k = 1
  kr= k + 1  ! the K-level of Grell is equivalent to the BRAMS K+1-level

  dz        =  z_cup(2)-z_cup(1)

  detdo1    = edt*zd(2)*  cdd(1)*dz
  detdo2    = edt*zd(1)
  entdo     = edt*zd(2)*entrd(1)*dz
  subin     =-edt*zd(2)           

  cfxup1(kr,i,j) = 0.
  cfxdn1(kr,i,j) =-edt*xmb* zd(1)
  dfxup1(kr,i,j) = 0.
  efxup1(kr,i,j) = 0.
  dfxdn1(kr,i,j) = xmb*(detdo1+detdo2) !edt already is at detdo1,2
  efxdn1(kr,i,j) =-xmb* entdo          !edt already is at entdo
 

!------------------------------------------------------------------------------------------!
! Checking the mass conservation                                                           !
!------------------------------------------------------------------------------------------!
  totmas = detdo1+detdo2-entdo+subin
  if(abs(totmas) > 1.e-6) then
    write (unit=*,fmt='(a)')                 '----------- Subroutine Get_convflx ----------'
    write(unit=*, fmt='(4(a,1x,i3,1x))')     '  K= ',k,'   I=',i,'   J=',j,'   IENS=',iens
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  subin=  ',    subin,'entdo=  ',entdo
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  detdo1= ',   detdo1,'detdo2= ',detdo2
    write(unit=*, fmt='(1(a,1x,es10.3,1x))') '  totmas= ',   totmas
    write(unit=*, fmt='(a)')                 '---------------------------------------------'
    stop 'The model will stop since it is not conserving mass...' 
  end if
end if

return
end subroutine get_convflx
!------------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------------!
! Subroutine output_mean_at_inst_analysis_advec                                            !
! Developed by Saulo Freitas (CPTEC/INPE)                                                  !
! Modified by Marcos Longo (EPS/Harvard U.) - Turbulence needs to be called separately,    !
!                                             because parameters vary accordinly to idiffk !
!                                                                                          !
!   The aim of this subroutine is simply save the mean values of some advection fluxes as  !
! at the regular and lite analysis.                                                        !
!------------------------------------------------------------------------------------------!
subroutine output_mean_at_inst_analysis_advec(m1,m2,m3,afxub,afxvb,afxwb,afxum,afxvm,afxwm)
implicit none
integer, intent(in)                    :: m1,m2,m3
real, dimension(m1,m2,m3), intent(in)  :: afxum, afxvm, afxwm
real, dimension(m1,m2,m3), intent(out) :: afxub, afxvb, afxwb
integer                                :: i, j, k

do k=1, m1
  do i=1, m2
    do j=1, m3
      afxub(k,i,j)=afxum(k,i,j)
      afxvb(k,i,j)=afxvm(k,i,j)
      afxwb(k,i,j)=afxwm(k,i,j)
    end do
  end do
end do

return
end subroutine output_mean_at_inst_analysis_advec
!------------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------------!
! Subroutine output_mean_at_inst_analysis_tke                                              !
! Developed by Saulo Freitas (CPTEC/INPE)                                                  !
! Modified by Marcos Longo (EPS/Harvard U.) - Turbulence needs to be called separately,    !
!                                             because just IDIFFK=1,4,5,6,7 produces TKE   !
!                                                                                          !
!   The aim of this subroutine is simply save the mean TKE value at the regular and lite   !
! analysis.                                                                                !
!------------------------------------------------------------------------------------------!
subroutine output_mean_at_inst_analysis_tke(m1,m2,m3,tkepb,tkepm)
implicit none
integer, intent(in)                    :: m1,m2,m3
real, dimension(m1,m2,m3), intent(in)  :: tkepm
real, dimension(m1,m2,m3), intent(out) :: tkepb
integer                                :: i, j, k

do k=1, m1
  do i=1, m2
    do j=1, m3
      tkepb(k,i,j)=tkepm(k,i,j)
    end do
  end do
end do

return
end subroutine output_mean_at_inst_analysis_tke
!------------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------------!
! Subroutine output_mean_at_inst_analysis_turb                                             !
! Developed by Saulo Freitas (CPTEC/INPE)                                                  !
! Modified by Marcos Longo (EPS/Harvard U.) - Turbulence needs to be called separately,    !
!                                             because just IDIFFK=7 produces sigmaw and Tl !
!                                                                                          !
!   The aim of this subroutine is simply save the mean values of Lagrangian time scale and !
! Sigma-w at the regular and lite analysis.                                                !
!------------------------------------------------------------------------------------------!
subroutine output_mean_at_inst_analysis_turb(m1,m2,m3,sigwb,tlb,sigwm,tlm)
implicit none
integer, intent(in)                    :: m1,m2,m3
real, dimension(m1,m2,m3), intent(in)  :: sigwm, tlm
real, dimension(m1,m2,m3), intent(out) :: sigwb, tlb
integer                                :: i, j, k

do k=1, m1
  do i=1, m2
    do j=1, m3
      sigwb(k,i,j)=sigwm(k,i,j)
      tlb  (k,i,j)=tlm  (k,i,j)
    end do
  end do
end do

return
end subroutine output_mean_at_inst_analysis_turb
!------------------------------------------------------------------------------------------!

subroutine get_adv_fluxes_for_stilt(mzp,mxp,myp,ia,iz,ja,jz)
  use mem_scratch, only: scratch, vctr1, vctr2
  use mem_grid, only: ngrid, nzpmax, grid_g, dtlt, if_adap, jdim, time, &
       zt, zm, dzm, dzt, hw4

  use mem_basic, only: basic_g

 
  use monotonic_adv, ONLY: advmnt 

  implicit none
  include "i8.h"
  integer :: mzp,mxp,myp,ia,iz,ja,jz!,izu,jzv,mynum,n
  integer(kind=i8) :: mxyzp
  

  real :: dtlto2
  integer :: i,j,k,ind

 
 
  mxyzp = mxp * myp * mzp

     dtlto2 = .5 * dtlt
     ind = 0
     do j = 1,myp
        do i = 1,mxp
           do k = 1,mzp
              ind = ind + 1
              scratch%vt3da(ind) = (basic_g(ngrid)%up(k,i,j)  &
                   + basic_g(ngrid)%uc(k,i,j)) * dtlto2
              scratch%vt3db(ind) = (basic_g(ngrid)%vp(k,i,j)  &
                   + basic_g(ngrid)%vc(k,i,j)) * dtlto2
              scratch%vt3dc(ind) = (basic_g(ngrid)%wp(k,i,j)  &
                   + basic_g(ngrid)%wc(k,i,j)) * dtlto2
           enddo
        enddo
     enddo

        call fa_preptc_for_stilt(mzp,mxp,myp        &
             ,scratch%vt3da        (1)     ,scratch%vt3db        (1)      &
             ,scratch%vt3dc        (1)      & 
             ,basic_g(ngrid)%dn0   (1,1,1) ,basic_g(ngrid)%dn0u  (1,1,1)  &
             ,basic_g(ngrid)%dn0v  (1,1,1) ,grid_g(ngrid)%rtgt   (1,1)    &
             ,grid_g(ngrid)%rtgu   (1,1)   ,grid_g(ngrid)%rtgv   (1,1)    &
             ,grid_g(ngrid)%fmapt  (1,1)   ,grid_g(ngrid)%fmapui (1,1)    &
             ,grid_g(ngrid)%fmapvi (1,1)   ,grid_g(ngrid)%f13t   (1,1)    &
             ,grid_g(ngrid)%f23t   (1,1)   ,grid_g(ngrid)%dxu    (1,1)    &
             ,grid_g(ngrid)%dyv    (1,1)   ,grid_g(ngrid)%dxt    (1,1)    &
             ,grid_g(ngrid)%dyt    (1,1)                            )



end subroutine get_adv_fluxes_for_stilt

!     *********************************************************************

subroutine fa_preptc_for_stilt(m1,m2,m3,vt3da,vt3db,vt3dc &
     !,vt3dd,vt3de,vt3df ,vt3dh,vt3di,vt3dj,vt3dk
     ,dn0,dn0u,dn0v  &
     ,rtgt,rtgu,rtgv,fmapt,fmapui,fmapvi,f13t,f23t  &
     ,dxu,dyv,dxt,dyt)

  use mem_grid, only: hw4,dzm,zm,zt
  use mem_scratch, only: vctr1,vctr2

  implicit none

  integer :: m1,m2,m3,j,i,k,im,ip,jm,jp

  real :: c1,c2,c3,c4,rtgti

  real, dimension(m1,m2,m3) :: vt3da,vt3db,vt3dc &!&,vt3dd,vt3de,vt3df  &
                                                  !,vt3dh,vt3di,vt3dj,vt3dk&
       ,dn0,dn0u,dn0v

  real, dimension(m2,m3) :: rtgt,rtgu,rtgv,fmapt,fmapui,fmapvi,f13t,f23t  &
       ,dxu,dyv,dxt,dyt

  ! VT3DA, VT3DB, and VT3DC are input as the velocity components (averaged
  ! between past and current time levels) times dtlt.

  ! Add contribution to VT3DC from horiz winds crossing sloping sigma surfaces,
  !    and include 1/rtgt factor in VT3DC
  ! Compute half Courant numbers: VT3DD, VT3DE, and VT3DF
  ! Compute weight at scalar point: VT3DH
  ! Compute advective weights for the linear term: VT3DI, VCTR1, and VCTR2

  do j = 1,m3
     jm = max(1,j-1)
     jp = min(m3,j+1)
     do i = 1,m2
        im = max(1,i-1)
        ip = min(m2,i+1)
        rtgti = 1. / rtgt(i,j)
        c1 = .5 * dxu(i,j)
        c2 = .5 * dyv(i,j)
        c3 = dxt(i,j) * fmapt(i,j) * rtgti
        c4 = dyt(i,j) * fmapt(i,j) * rtgti

        do k = 1,m1-1
           vt3dc(k,i,j) = ((vt3da(k,i,j) + vt3da(k+1,i,j)  &
                + vt3da(k,im,j) + vt3da(k+1,im,j)) * f13t(i,j)  &
                + (vt3db(k,i,j) + vt3db(k+1,i,j) + vt3db(k,i,jm)  &
                + vt3db(k+1,i,jm)) * f23t(i,j)) * hw4(k)  &
                + vt3dc(k,i,j) * rtgti
        enddo
    enddo
  enddo

  do k = 1,m1-1
     vctr1(k) = (zt(k+1) - zm(k)) * dzm(k)
     vctr2(k) =  (zm(k) - zt(k)) * dzm(k)
  enddo

  ! Convert velocity components * dtlt (VT3DA, VT3DB, VT3DC)
  ! into mass fluxes times dtlt.

  do j = 1,m3
     do i = 1,m2
        c1 = fmapui(i,j) * rtgu(i,j)
        c2 = fmapvi(i,j) * rtgv(i,j)
        do k = 1,m1-1
           vt3da(k,i,j) = vt3da(k,i,j) * c1 * dn0u(k,i,j)
           vt3db(k,i,j) = vt3db(k,i,j) * c2 * dn0v(k,i,j)
           vt3dc(k,i,j) = vt3dc(k,i,j) * .5  &
                * (dn0(k,i,j) + dn0(k+1,i,j))
        enddo
       vt3da(m1,i,j) = vt3da(m1-1,i,j) 
       vt3db(m1,i,j) = vt3db(m1-1,i,j) 
       vt3dc(m1,i,j) = vt3dc(m1-1,i,j) 
     enddo
  enddo

  return
end subroutine fa_preptc_for_stilt

!     *********************************************************************

subroutine prepare_timeavg_driver(mzp,mxp,myp,ia,iz,ja,jz,dtlt,ifm,idiffkk)
   use mem_stilt , only :   stilt_g                 & ! intent(inout) 
                          , etime_turb             & ! intent(inout)
                          , frqmassave             & ! intent(in)
                          , zero_average_mass_turb ! ! subroutine
   use mem_turb , only : turb_g   ! ! intent(in)
   implicit none
   !------ Arguments ----------------------------------------------------------------------!
   integer, intent(in) :: mzp, mxp, myp, ia, iz, ja, jz, ifm, idiffkk
   real   , intent(in) :: dtlt
   !---------------------------------------------------------------------------------------!
   

   !---------------------------------------------------------------------------------------!
   !     Update the step counter for turbulent variables.  If this exceeds the maximum     !
   ! number of steps, it means it is time to reset it and start integrating again.         ! 
   !---------------------------------------------------------------------------------------!
   etime_turb(ifm) = etime_turb(ifm) + dtlt
   if (etime_turb(ifm) > frqmassave + 0.1 * dtlt) then
      etime_turb(ifm) = dtlt
      call zero_average_mass_turb(stilt_g(ifm))
   end if
   !---------------------------------------------------------------------------------------!


   !----- Nakanishi-Niino closure, TKE, sig-W and Lagrangian time scale exist. ---!
      call prepare_timeavg_to_mass(mzp,mxp,myp,ia,iz,ja,jz,dtlt                            &
                                  ,turb_g(ifm)%tkep     ,stilt_g(ifm)%tkepb                 )

      call prepare_timeavg_to_mass(mzp,mxp,myp,ia,iz,ja,jz,dtlt                            &
                                  ,stilt_g(ifm)%sigw     ,stilt_g(ifm)%sigwb                 )

      call prepare_timeavg_to_mass(mzp,mxp,myp,ia,iz,ja,jz,dtlt                            &
                                  ,stilt_g(ifm)%ltscale  ,stilt_g(ifm)%ltscaleb              )

end subroutine prepare_timeavg_driver
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
! Subroutine prepare_timeavg_to_mass                                                       !
!                                                                                          !
!   The aim of this subroutine is simply save the mean value of any variable at the        !
! regular and lite analysis.                                                               !
!------------------------------------------------------------------------------------------!
subroutine prepare_timeavg_to_mass(m1,m2,m3,ia,iz,ja,jz,dtlt,var,avgvar)
   use mem_stilt, only : frqmassave
   implicit none
   integer, intent(in)                      :: m1,m2,m3
   integer, intent(in)                      :: ia,iz,ja,jz
   real                     , intent(in)    :: dtlt
   real, dimension(m1,m2,m3), intent(in)    :: var
   real, dimension(m1,m2,m3), intent(inout) :: avgvar
   integer                                  :: i, j, k
   real                                     :: timefac

   timefac = dtlt/frqmassave

   do k=1, m1
     do i= ia, iz
       do j= ja, jz
         avgvar(k,i,j)= avgvar(k,i,j) + var(k,i,j) * timefac
       end do
     end do
   end do

   return
end subroutine prepare_timeavg_to_mass
!==========================================================================================!
!==========================================================================================!

