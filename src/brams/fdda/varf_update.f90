!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine varf_adap(n1,n2,n3,varu,varv,varp,vart,varr,topta)

  use mem_scratch, only: &
       vctr1, vctr2, vctr2, vctr3, vctr4, vctr10, &
       vctr11, vctr12, vctr13, vctr14, vctr15, vctr16

  use mem_grid, only: &
       ztop, ztn, ngrid

  use rconstants, only: &
       g

  implicit none

  integer :: n1,n2,n3
  real, dimension(n1,n2,n3) :: varu,varv,varp,vart,varr
  real, dimension(n2,n3) :: topta

  integer :: i,j,k

  ! Interpolate from sigma-z varfile vertical coords to ADAP grid


  do j=1,n3
     do i=1,n2

        do k=1,n1
           vctr10(k)=topta(i,j) + (1.-topta(i,j)/ztop)*ztn(k,ngrid)
        enddo
        vctr1(1:n1)=varu(1:n1,i,j)
        vctr2(1:n1)=varv(1:n1,i,j)
        vctr3(1:n1)=vart(1:n1,i,j)
        vctr4(1:n1)=varr(1:n1,i,j)
        call htint2(n1,vctr1(1),vctr10,n1,vctr11(1),ztn(1,ngrid))
        call htint2(n1,vctr2(1),vctr10,n1,vctr12(1),ztn(1,ngrid))
        call htint2(n1,vctr3(1),vctr10,n1,vctr13(1),ztn(1,ngrid))
        call htint2(n1,vctr4(1),vctr10,n1,vctr14(1),ztn(1,ngrid))

        ! Do hydrostatic balance
        do k=1,n1
           vctr15(k) = vctr13(k)* (1.+.61*vctr14(k))
        enddo

        vctr16(n1)= varp(n1,i,j) + g * (ztn(n1,ngrid) - vctr10(n1) )  &
             / vctr15(n1)
        do k = n1-1,1,-1
           vctr16(k) = vctr16(k+1) + g * (ztn(k+1,ngrid)-ztn(k,ngrid))  &
                /((vctr15(k)+vctr15(k+1))*.5)
        enddo

        varu(1:n1,i,j)= vctr11(1:n1)
        varv(1:n1,i,j)= vctr12(1:n1)
        vart(1:n1,i,j)= vctr13(1:n1)
        varr(1:n1,i,j)= vctr14(1:n1)
        varp(1:n1,i,j)= vctr16(1:n1)

     enddo
  enddo

  return
end subroutine varf_adap

!     **************************************************************





subroutine PrtOpt()

  use ISO_FORTRAN_ENV, only: OUTPUT_UNIT
  use mem_grid, only: nnzp, ztn
  use rconstants, only: p00, cp, cpor
  use ref_sounding, only: &
       iref, jref, topref, &
       u01dn, v01dn, dn01dn, pi01dn, th01dn, rt01dn

  implicit none
  include "constants.f90"
  integer :: k
  real :: vctr1(nnzp(1))
  real :: vctr2(nnzp(1))
  character(len=*), parameter :: pfmt="(/,'--------REFERENCE STATE (AVERAGE OVER DOMAIN)'" &
       //",'   SFC ELEV (M)= ',F6.1,'-------------'"  &
       //",//,4X,'Z',6X,'U01D',4X,'V01D',4X,'DN01D',4X"  &
       //",'PI01D',4X,'PRESS',4X,'TH01D',4X,'THD',6X,'RT01D'"  &
       //",/,3X,'(m)',5X,'(m/s)',3X,'(m/s)',2X,'(kg/m3)',2X"  &
       //",'(J/kgK)',4X,'(Pa)',5X,'(K)',5X,'(K)',5X,'(kg/kg)'"  &
       //",//,(1X,F7.1,2F8.2,F8.3,F10.2,F10.1,2F8.2,F10.5))"

  do k=1,nnzp(1)
     vctr1(k)=p00*(pi01dn(k,1)/cp)**cpor
     vctr2(k)=th01dn(k,1)/(1.+.61*rt01dn(k,1))
  enddo
  open(unit=22,file='brams.log',position='append',action='write')
  write(unit=22,fmt=pfmt)topref,(ztn(k,1),u01dn(k,1),v01dn(k,1)  &
       ,dn01dn(k,1),pi01dn(k,1),vctr1(k),th01dn(k,1),vctr2(k)  &
       ,rt01dn(k,1),k=1,nnzp(1))
310 format(/,'--------REFERENCE STATE (AVERAGE OVER DOMAIN)' &
       ,'   SFC ELEV (M)= ',F6.1,'-------------'  &
       ,//,4X,'Z',6X,'U01D',4X,'V01D',4X,'DN01D',4X  &
       ,'PI01D',4X,'PRESS',4X,'TH01D',4X,'THD',6X,'RT01D'  &
       ,/,3X,'(m)',5X,'(m/s)',3X,'(m/s)',2X,'(kg/m3)',2X  &
       ,'(J/kgK)',4X,'(Pa)',5X,'(K)',5X,'(K)',5X,'(kg/kg)'  &
       ,//,(1X,F7.1,2F8.2,F8.3,F10.2,F10.1,2F8.2,F10.5))
   close(unit=22)
end subroutine PrtOpt

!--(DMK-CCATT-INI)-----------------------------------------------------
!-----------------------------

subroutine hi_interpInitial4(n1,n2,n3,vn,xm1,xt1,ym1,yt1,zm1,zt1,plat1,plon1,topt1,ztop1,m1,m2,m3,vm,ngm,ngr1,vname,idim)

 use mem_grid, only: grid_g, ztn
 use mem_scratch, only: vctr1, vctr2, vctr3, vctr10

 implicit none

 character (len=*), intent(in)			:: vname

 integer, intent(in) 				:: n1
 integer, intent(in) 				:: n2
 integer, intent(in) 				:: n3
 integer, intent(in) 				:: ngr1
 integer, intent(in) 				:: idim
 integer, intent(in) 				:: m1
 integer, intent(in)				:: m2
 integer, intent(in)				:: m3
 integer, intent(in)				:: ngm

 real, intent(in)				:: plat1
 real, intent(in) 				:: plon1
 real, intent(in) 				:: ztop1
 real, dimension(n2,n3,n1), intent(in)		:: vn
 real, dimension(n2,n3), intent(in)      	:: topt1
 real, dimension(n1), intent(in)         	:: zm1
 real, dimension(n1), intent(in)		:: zt1
 real, dimension(n2), intent(in)	 	:: xm1
 real, dimension(n2), intent(in)	 	:: xt1
 real, dimension(n3), intent(in)	 	:: ym1
 real, dimension(n3), intent(in)	 	:: yt1
 real, dimension(m1,m2,m3), intent(out)		:: vm



 integer :: i,j,k,np,ii,jj
 real :: xxm,yym,fixxm,fiyym,topoh,rtgth



  	!!print*, 'inside hi_interpInitil4', minval(vn), maxval(vn)


	do j=1,m3
		do i=1,m2


		        ! Find real grid point x,y relative to history file
     			call ll_xy(grid_g(ngm)%glat(i,j),grid_g(ngm)%glon(i,j), plat1, plon1, xxm, yym)

        		! See if point is on this grid.
        		if(xxm < xm1(1) .or. xxm > xm1(n2-1) .or. yym < ym1(1) .or. yym > ym1(n3-1) ) then

           			if (ngr1 == 1) then

              				! We are on the input coarse grid and point is not on this grid.
              				!    Stop immediately...
              				print*, grid_g(ngm)%glat(i,j),grid_g(ngm)%glon(i,j),plat1,plon1,xxm,yym
              				print*,xm1(1),xm1(n2-1),xxm
              				print*,ym1(1),ym1(n3-1),yym
              				print*, 'His_init: grid point not on history file grids'
              				stop 'his_init: point OB'
           			else
              				! Otherwise, go to next point
              				cycle
           			endif

        		endif
			! We are okay horizontally, now interpolate vertical column from
        		!     field

        		! Find x,y grid point locations on input field.
        		!     Assuming constant spacing and deal with stagger

        		if(vname == 'UP' .or. vname == 'UC') then
           			fixxm=1.+(xxm-xm1(1))/(xm1(2)-xm1(1))
        		else
           			fixxm=1.+(xxm-xt1(1))/(xt1(2)-xt1(1))
        		endif

        		if(vname == 'VP' .or. vname == 'VC') then
           			fiyym=1.+(yym-ym1(1))/(ym1(2)-ym1(1))
        		else
           			fiyym=1.+(yym-yt1(1))/(yt1(2)-yt1(1))
        		endif


           		do k=1,n1
              			call gdtost2(vn(1,1,k),n2,n3,fixxm,fiyym,vctr1(k))
           		enddo

           		if (idim == 3) then
              			! Interpolate this column vertically to actual grid if 3d variable
              			call gdtost2(topt1,n2,n3,fixxm,fiyym,topoh)

				rtgth=1.-topoh/ztop1

				do k=1,m1
                 			! Actual grid level heights
		 			vctr2(k) = grid_g(ngm)%topt(i,j) + ztn(k,1) * grid_g(ngm)%rtgt(i,j)
		      		enddo

				do k=1,n1
                 			! History grid level heights
                 			vctr3(k)= topoh + zt1(k) *rtgth
              			enddo

              			! Interpolate vertically

              			call htint(n1,vctr1(1),vctr3(1),m1,vctr10(1),vctr2(1))
              			vm(1:m1,i,j) = vctr10(1:m1)

			elseif (idim == 2) then
              			vm(1,i,j)=vctr1(1)
           		endif

     		end do
	end do

	if(vname == 'UP' .or. vname == 'UC') then
     		!CALL hi_avgu(m1,m2,m3,vm) ! Modif.by Alvaro L.Fazenda
     		call hi_avgu(m1,m2,m3,vm)
  	elseif (vname == 'VP' .or. vname == 'VC') then
     		call hi_avgv(m1,m2,m3,vm)
  	elseif (vname == 'WP' .or. vname == 'WC') then
     		call hi_avgw(m1,m2,m3,vm)
  	endif

end subroutine hi_interpInitial4
!--(DMK-CCATT-FIM)-----------------------------------------------------
