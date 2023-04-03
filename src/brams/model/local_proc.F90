module dtset

  use grid_dims, only : &
       maxgrds           ! INTENT(IN)

  implicit none

  real    :: ssodx(maxgrds)
  integer :: idelx(maxgrds)
  real    :: sscourn(maxgrds)
  logical :: dtSet_firstTime
  real,allocatable :: validDts(:)
  integer :: totValidDt
  real :: maxCflPercent

contains

  subroutine dtset_new(mynum, nndtflg, dxtmax_local)

    use mem_grid, only : &
         ideltat,        & ! INTENT(IN)
         ngrids,         & ! INTENT(IN)
         deltaxn,        & ! INTENT(IN)
         deltayn,        & ! intent(in)
         nnxp,           & ! INTENT(IN)
         nnyp,           & ! INTENT(IN)
         nnzp,           & ! INTENT(IN)
         grid_g,         & ! INTENT(IN)
         sspct,          & ! INTENT(INOUT) ! Calculada localmente nesta rotina
         nxtnest,        & ! INTENT(IN)
         dtlongn,        & ! INTENT(INOUT) ! Calculada localmente nesta rotina
         dtlong,         & ! INTENT(IN)    ! Falta passar para o escravo-ok
         nndtrat,        & ! INTENT(INOUT) ! Valor inicial passado-ok
         ! Calculado atualizacao nesta rot.
         ! Retirar da com.:Mest->escravo
         nnacoust,       & ! INTENT(INOUT) ! Calculada localmente nesta rotina
         nacoust,        & ! INTENT(IN)
         cflxy,          & ! INTENT(IN)    ! Calculado localmente em CFL
         cflz,           & ! INTENT(IN)    ! A ser calculado em modsched local
         iflag,          & ! INTENT(INOUT) ! definida localmente nesta rotina
         timmax,         & ! INTENT(IN)
         time,           & ! INTENT(IN)
	 dyncore_flag

    ! "cflxy" e "cflz" calculado localmente, precisa ser passado a todos os
    ! escravos para se determinar o maior valor.

    use rconstants, only : &
         cp,               & ! INTENT(IN)
         cv,               & ! INTENT(IN)
         rgas                ! INTENT(IN)

    use ref_sounding, only : &
         th01dn,             & ! INTENT(IN)
         pi01dn                ! INTENT(IN)
    ! Receber apos atualizacao-chamada de varf_read em RAMS_OUTPUT

    use grid_dims, only : &
         maxgrds,         & ! INTENT(IN)
         nzpmax             ! INTENT(IN)

    use io_params, only : &
         frqanl             ! INTENT(IN)

    use mem_stilt, only : &
         iexev             ! INTENT(IN)

    use ReadBcst, only: &
        Broadcast

    use node_mod, only:  &
       master_num

    use dump, only: &
      dumpMessage

    implicit none


    include "constants.f90"
    ! Arguments:
    integer, intent(in)  :: mynum
    integer, intent(out) :: nndtflg
    real, intent(in)     :: dxtmax_local(maxgrds)

    ! Local variables:
    integer, parameter  :: ndx=37,ndt=42
    integer, dimension(maxgrds) :: idelt,nndtrat1,nnacoust1
    real, dimension(maxgrds) :: dtlongn1,dtmax_rk3,dtlongn2
    real, dimension(nzpmax) :: vctr1
    real :: delx(ndx), delt(ndt)

    real,parameter :: cflnumh = .90
    real,parameter :: cflnumv = .90

    integer :: iabsdt,ifm,id,n2,n3,k,nn2,nn3,icm,ntf,ii,i,ierr
    real :: ssmax,tmax,dxtmax,sscnmax,sspct0,cflxyz,timeleft

    real :: newDeltaT

    real :: dxta, dxtb, dxtc, dxtd
    character(len=8) :: c0, c1, c2, c3, c4
    character(len=*), parameter :: h="**(dtset_new)**"
    character(len=*), parameter :: header="**(dtset_new)**"
    character(len=*), parameter :: version="local_proc.f90"

    delx(:) = (/ &
         200000.,150000.,100000.,80000.,70000.,60000.,40000.,  &
         030000., 20000., 10000., 6000., 4000., 3000., 2000.,  &
         001000.,   800.,   600.,  500.,  400.,  300.,  200.,  &
         000150.,   100.,    80.,   60.,   50.,   40.,   30.,  &
         000020.,    10.,     8.,    6.,    5.,    4.,    3.,  &
         000002.,     1. /)

    delt(:) = (/ &
         300.,  240.,   180.,  150.,  120.,   90.,   60.,  &
         050.,   40.,    30.,   20.,   15.,   12.,   10.,  &
         006.,    5.,     4.,    3.,   2.5,   2.0,   1.5,  &
         01.2,    1.0,     .8,    .6,   .5,    .4,    .3,  &
         00.2,     .1,     .08,   .06,  .05,   .04,   .03,  &
         00.02,    .01,    .008,  .006, .005,  .004,  .003/)

    idelx(1) = 0

    iabsdt = abs(ideltat)

    dtmax_rk3(:) = 0.

    ! On the first call to this subroutine, initialize idelx, ssodx, dtlongn,
    ! nnacoust, and if required, nndtrat.

!    if ( idelx(1)==0 ) then
    if(dtSet_firstTime) then
      dtSet_firstTime=.false.
      call createDts(frqanl)

      do ifm = 1,ngrids
        do id = ndx,1,-1
          if (delx(id)<=deltaxn(ifm))  idelx(ifm) = id
        enddo
        n2 = nnxp(ifm)
        n3 = nnyp(ifm)
        do k = 1,nnzp(ifm)
          vctr1(k) = th01dn(k,1) * pi01dn(k,1) / cp
        enddo
        tmax = maxval(vctr1(1:nnzp(ifm)))
        ssmax = sqrt(cp / cv * rgas * tmax)
        nn2 = nnxp(ifm)
        nn3 = nnyp(ifm)
        if (mynum==0) then
          ! Master process
          dxtmax = max(grid_g(ifm)%dxt(1,1), &
                  grid_g(ifm)%dxt(nn2,1),       &
                  grid_g(ifm)%dxt(nn2,nn3),     &
                  grid_g(ifm)%dxt(1,nn3)        )
        else
          ! Slave (node) process
          dxtmax = dxtmax_local(ifm)
        endif
        ssodx(ifm) = ssmax * dxtmax
      enddo

      !srf   if ( ideltat==0 .or. ideltat==3) then
!      if ( ideltat==0) then
        sspct = 1.
        do ifm = 1,ngrids
          icm = nxtnest(ifm)
          if ( icm==0 ) then
            dtlongn(ifm) = dtlong
          else
            dtlongn(ifm) = dtlongn(icm) / nndtrat(ifm)
          endif
          nnacoust(ifm) = nacoust

          sscourn(ifm) = 2. * ssodx(ifm) * dtlongn(ifm)

          sspct = min(sspct, .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))
          !-for RK3 only
          dtmax_rk3(ifm) = dtlongn(ifm)
        enddo
        !---srf: not allowing the artifical reduction of the sound speed for RK3
        if(dyncore_flag == 2) sspct = 1.
    !   else !ideltat/=0
    !     sscnmax = 0.
    !     do ifm = 1,ngrids
    !       icm = nxtnest(ifm)
    !       dtlongn(ifm) = delt(idelx(ifm)+iabsdt-1)
    !       !---srf: max dt allowed by RK3
    !       if(dyncore_flag >= 2) then
    !         if (mynum==0) then
    !           dxtmax = max(grid_g(ifm)%dxt(1,1), &
    !                      grid_g(ifm)%dxt(nn2,1),       &
    !                      grid_g(ifm)%dxt(nn2,nn3),     &
    !                      grid_g(ifm)%dxt(1,nn3)        )
    !         else
    !           ! Slave (node) process
    !           dxtmax = dxtmax_local(ifm)
    !         endif
    !         dtlongn  (ifm) = 6. * (1./dxtmax)/1000.
    !         ntf = nint(frqanl / dtlongn(1))
    !         dtlongn(ifm)   = frqanl / ntf
	   !        dtmax_rk3(ifm) = dtlongn(ifm)
    !       endif
    !       !print *, '1st time:',dtlongn  (ifm),dxtmax,1./dxtmax; call flush(6)
    !       ! For coarse grid(s) adjust dtlongn so that it is an integer
    !       ! divisor of FRQANL.  For nested grids, compute nndtrat(ifm) as the
    !       ! first integer greater than or equal to the timestep ratio between
    !       ! a nested grid's parent and the nested grid. Then adjust
    !       ! dtlongn(ifm) for the nested grid to be the parent grid timestep
    !       ! divided by nndtrat(ifm).
    !       if ( icm==0 ) then
    !         ntf = nint(frqanl / dtlongn(1))
    !         dtlongn(ifm)   = frqanl / ntf
    !         dtmax_rk3(ifm) = dtlongn(ifm)
    !       else
    !         nndtrat(ifm)   = min(10, nint(dtlongn(icm) / dtlongn(ifm)))
    !         dtlongn(ifm)   = dtlongn(icm) / nndtrat(ifm)
    !         dtmax_rk3(ifm) = dtlongn(ifm)
    !       endif
    !       ! Compute sst courant numbers (sstcourn(ifm)) for long timestep
    !       ! dtlongn.
    !       sscourn(ifm) = 2. * ssodx(ifm) * dtlongn(ifm)
    !       if ( sscourn(ifm)>sscnmax)  sscnmax = sscourn(ifm)
    !     enddo
    !     !---srf: not allowing the artifical reduction of the sound speed for RK3
    !     if(dyncore_flag >= 2) then
    !       sspct = 1.
    !       nnacoust(:) = nacoust
    !     else
    !       !Define trial sspct0 in terms of sscnmax using nonlinear formula
    !       ! intended to increase nnacoust as sspct decreases, but limit sspct0
    !       ! to a minimum of .2.
    !       sspct0 = min(1., (.95/sscnmax)**.5)
    !       if ( sspct0<.2 ) then
    !         print*, 'Sound speed percent is forced to be too low'
    !         stop 'low_sspct0'
    !       endif
    !       sspct = 1.
    !       do ifm = 1,ngrids
    !         nnacoust(ifm) = max(2, nint(sspct0 * sscourn(ifm) * 2. / .95))
    !         sspct = min(sspct, .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))
    !       enddo
    !     endif !dyncore_flag<2
    !     call adjustDt(dtlongn(1))
    !   endif !ideltat/=0
    endif !first time

    ! check Courant numbers
    nndtflg = 0
    if ( ideltat>=0) then
      do ifm = 1,ngrids
        cflxyz = max(cflxy(ifm)/cflnumh,cflz(ifm)/cflnumv)
        !print *,'1cflxyz:',cflxyz,cflxy(ifm),cflz(ifm)
        if ( cflxyz>1. ) then
          !print*,'cfl xy, z, cflnumh, cflnumv = ',cflxy(ifm),cflz(ifm),cflnumh, cflnumv,mynum,iflag; call flush(6)
          iflag = 1
          write(c0,"(f8.2)") cflxyz
          write(c3,"(f8.2)") cflxy(ifm)
          write(c4,"(f8.2)") cflz(ifm)
          write(c1,"(i8)") mynum
          write(c2,"(i8)") ifm
             !
          if ( dyncore_flag >= 2) then
            iflag=0
          ! don't break here, there exist a better suited criterion in subroutine 'cfl1'
          else
            call fatal_error(h//" Model will stop because CFL limit exceeded:"//&
                   &" clfxy="//trim(adjustl(c0))//" at proc "//trim(adjustl(c1))//&
                   &" on grid "//trim(adjustl(c2))//" cflxy="//trim(adjustl(c3))//&
	           &" cflz="//trim(adjustl(c4)))
          endif
	      endif
      enddo

      !       if(dtSet_firstTime) then
      !         !On the first time use a defined dtlong:
      !         if(ideltat==3 .and. dyncore_flag==2) dtlongn(1) = 6.*min(deltaxn(1), deltayn(1))/1000.
      !iErrNumber=dumpMessage(c_tty,c_yes,header,version,c_notice,'Using a adjustable DeltaT',dtlongn(1),'F6.1')
      !ATENCAO: ajustar para multiplos de ANL,POS,HIS
      !         dtSet_firstTime=.false.
      !         return
      !       endif

    !   !LFR ------------------------>
    !   if(ideltat==3) then
    !     sscnmax=0.
    !     dtlongn1(1) = dtlongn(1)
    !     if ( dyncore_flag >= 2 ) then
    !       ifm = 1
    !       if(mynum==0) then
    !         dxtmax = max(grid_g(ifm)%dxt(1,1), &
    !         grid_g(ifm)%dxt(nn2,1),       &
    !         grid_g(ifm)%dxt(nn2,nn3),     &
    !         grid_g(ifm)%dxt(1,nn3)        )
    !       else
    !         ! Slave (node) process
    !         dxtmax = dxtmax_local(ifm)
    !       endif
    !       dtlongn2  (ifm) = 6. * (1./dxtmax)/1000.
    !       ntf = nint(frqanl / dtlongn2(1))
    !       dtlongn2(1)   = frqanl / ntf
    !       dtmax_rk3(ifm) = dtlongn2(1)
    !       cflxyz = max(cflxy(1)/cflnumh,cflz(1)/cflnumv)
    !       if ( cflxyz > 0.7 * 1.62 ) then !Desestabilizando: decrementar
    !         !Adjust for 50% of initial time
    !         !print*,"3cflxyz > 0.7 * 1.62",mynum,cflxyz,dtlongn(1),dtmax_rk3  (1);call flush(6)
    !         dtlongn(1)=0.5*dtlongn(1)
    !         !   Verify if dtlong is too low (less then 10%)
    !         if( dtlongn(1)<= 0.1* dtmax_rk3  (1) ) then
    !           iErrNumber=dumpMessage(c_tty,c_yes,header,version,c_warning,'dtlong too low!')
    !           dtlongn(1) = dtlongn1(1) !set back
    !         endif
    !         !-
    !         ntf = nint(frqanl / dtlongn(1))
    !         dtlongn(1) = frqanl / ntf
    !       elseif (cflxyz <=0.7 * 1.62 ) then !Sobreestabilizado: incrementar
    !         !adjust for 125% of previous time
    !         !print*,"4cflxyz < 0.7 * 1.62",mynum,cflxyz,dtlongn(1),1.25*dtlongn(1),dtmax_rk3  (1);call flush(6)
    !         dtlongn(1)=min(1.25*dtlongn(1),dtmax_rk3  (1))
    !         !-
    !         ntf = nint(frqanl / dtlongn(1))
    !         dtlongn(1) = frqanl / ntf
    !         !print*,"5cflxyz < 0.7 * 1.62",mynum,cflxyz,dtlongn(1),1.25*dtlongn(1),dtmax_rk3  (1);call flush(6)
    !         !verify if dtlong is greater then maximum value
    !         if(dtlongn(1)> dtlongn2(1))  then
    !           iErrNumber=dumpMessage(c_tty,c_yes,header,version,c_warning,'dtlong too high!')
    !           dtlongn(1) = dtmax_rk3  (1)  !set original
    !           ntf = nint(frqanl / dtlongn(1))
    !           dtlongn(1) = frqanl / ntf
    !         endif
    !       endif
    !       call adjustDt(dtlongn(1))
    !     endif

    !     ! Compute sst courant numbers (sstcourn(ifm))for long timestep dtlongn.
    !     sscourn(1) = 2. * ssodx(1) * dtlongn(1)
    !     !print *,'sscourn(2)=',sscourn(1),dtlongn(1),ssodx(1),2. * ssodx(1) * dtlongn(1)

    !     if (sscourn(1) .gt. sscnmax) sscnmax = sscourn(1)
    !     !---srf: not allowing the artifical reduction of the sound speed for RK3
    !     if(dyncore_flag >= 2) then
    !       sspct = 1.
    !       nnacoust(ifm) = nacoust
    !     else
    !       ! Define trial sspct0 in terms of sscnmax using nonlinear formula
    !       ! intended to increase nnacoust as sspct decreases, but limit sspct0 to
    !       ! a minimum of .2.
    !       sspct0 = min(1., (.95/sscnmax) ** .5)
    !       if ( sspct0<.2 ) then
    !         iErrNumber=dumpMessage(c_tty,c_yes,header,version,c_warning,'Sound speed percent is forced to be too low!')
    !       endif
    !       sspct = 1.
    !       do ifm = 1,ngrids
    !         !srf - for�ando nacoust > 2 para iexev=2
    !         !nnacoust(ifm) = max(2, nint(sspct0 * sscourn(ifm) * 2. / .95))
    !         if(iexev==1)  nnacoust(ifm) = max(2, nint(sspct0 * sscourn(ifm) * 2. / .95))
    !         if(iexev==2)  nnacoust(ifm) = max(3, nint(sspct0 * sscourn(ifm) * 2. / .95))
    !         !srf
    !         if(nnacoust(ifm)>4) nnacoust=6
    !         if(nnacoust(ifm)<5) nnacoust=4

    !         sspct = min(sspct, .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))

    !         ! If there are any updates to dtlongn(ifm), print out new values.

    !         !if (abs(dtlongn(ifm) - dtlongn1(ifm)) > 1.e-3) then
    !         !   write(6,122) ifm,nndtrat(ifm),nnacoust(ifm),mynum,dtlongn(ifm),dtlongn1(ifm)
    !         !endif

    !         ! If there are any updates to nndtrat(ifm) or others, set
    !         ! nndtflg = 1 to flag new call to subroutine modsched and send new
    !         ! stuff to nodes.

    !         if (nndtrat(ifm) /= nndtrat1(ifm) .or. &
    !            nnacoust(ifm) /= nnacoust1(ifm) .or. &
    !            dtlongn(ifm) /= dtlongn1(ifm) ) nndtflg = 1

    !       enddo
    !     endif
    !   endif
    !   !LFR ------------------------<

    ! !if ( dyncore_flag == 2 .or. dyncore_flag == 3) then
	   ! ! don't break here, there exist a better suited criterion in subroutine 'cfl1'
    !   !endif
    ! else
    !   do ifm = 1,ngrids
    !     icm = nxtnest(ifm)
    !     cflxyz = max(cflxy(ifm)/cflnumh, cflz(ifm)/cflnumv)

    !     nndtrat1(ifm) = nndtrat(ifm)
    !     dtlongn1(ifm) = dtlongn(ifm)
    !     nnacoust1(ifm) = nnacoust(ifm)

    !     do id = ndt,1,-1
    !        if ( delt(id)*cflxyz<=dtlongn(ifm) ) idelt(ifm) = id
    !     enddo

    !     if ( idelt(ifm)>idelx(ifm)+4 ) then
    !       iflag = 1
    !       write(c1,"(i8)") mynum
    !       write(c2,"(i8)") ifm
    !       call fatal_error(h//" Model will stop because adjustable timestep"//&
    !               &" forced to become too small at proc "//trim(adjustl(c1))//&
    !               &" on grid "//trim(adjustl(c2)))
    !     else
    !       ii = max(idelx(ifm)+iabsdt-1,idelt(ifm))
    !       dtlongn(ifm) = delt(ii)

    !          ! For the coarse grid(s), adjust dtlongn(1) to not jump past an
    !          ! analysis write time or the end of the model run time.

    !          if ( icm==0 ) then
    !             timeleft =  &
    !                  min (timmax - time, frqanl - mod(time,frqanl))
    !             if ( dtlongn(1)>.95 * timeleft) then
    !                dtlongn(ifm) = timeleft
    !             endif
    !          else
    !             nndtrat(ifm) = min(10, nint(dtlongn(icm) / dtlongn(ifm)))
    !             dtlongn(ifm) = dtlongn(icm) / nndtrat(ifm)
    !          endif
    !       endif

    !       ! Compute sst courant numbers (sstcourn(ifm))for long timestep dtlongn.

    !       sscourn(ifm) = 2. * ssodx(ifm) * dtlongn(ifm)
    !       if (sscourn(ifm) .gt. sscnmax) sscnmax = sscourn(ifm)

    !    enddo

    !    !---srf: not allowing the artifical reduction of the sound speed for RK3
       if(dyncore_flag >= 2) then

         sspct = 1.

       else

       ! Define trial sspct0 in terms of sscnmax using nonlinear formula
       ! intended to increase nnacoust as sspct decreases, but limit sspct0 to
       ! a minimum of .2.

       sspct0 = min(1., (.95/sscnmax) ** .5)
       if ( sspct0<.2 ) then
          print*, 'Sound speed percent is forced to be too low'
          stop 'low_sspct0'
       endif

       sspct = 1.
       do ifm = 1,ngrids
          !srf - for�ando nacoust > 2 para iexev=2
          !nnacoust(ifm) = max(2, nint(sspct0 * sscourn(ifm) * 2. / .95))
          if(iexev==1)  nnacoust(ifm) = max(2, nint(sspct0 * sscourn(ifm) * 2. / .95))
          if(iexev==2)  nnacoust(ifm) = max(3, nint(sspct0 * sscourn(ifm) * 2. / .95))
          !srf

          sspct = min(sspct, .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))

          ! If there are any updates to dtlongn(ifm), print out new values.

          if (abs(dtlongn(ifm) - dtlongn1(ifm)) > 1.e-3) then
            write(6,122) ifm,nndtrat(ifm),nnacoust(ifm),mynum,dtlongn(ifm),dtlongn1(ifm)
122          format('Timestep update1: ngrid, nndtrat, nnacoust, mynum,'  &
                  ,'new dtlongn, old dtlongn = ',4i4,2f13.6)
          endif

          ! If there are any updates to nndtrat(ifm) or others, set
          ! nndtflg = 1 to flag new call to subroutine modsched and send new
          ! stuff to nodes.

          if (nndtrat(ifm) /= nndtrat1(ifm) .or. &
               nnacoust(ifm) /= nnacoust1(ifm) .or. &
               dtlongn(ifm) /= dtlongn1(ifm) ) nndtflg = 1

       enddo
        endif

    endif
!print *,'Check: ',dtlongn(1),nnacoust(1),nndtrat(1)
    return
  end subroutine dtset_new

  ! ******************************************************************

  subroutine dump_courn()

    use mem_grid, only : &
         ngrids,         & ! INTENT(IN)
         dtlongn,        & ! INTENT(IN)
         nndtrat,        & ! INTENT(IN)
         nnacoust,       & ! INTENT(IN)
         sspct,          & ! INTENT(IN)
         dyncore_flag
    use dump

    implicit none
    include "constants.f90"
    integer :: ifm,err

    ! Print out initial values of dtlongn, nndtrat, nnacoust, sscourn,
    ! and sspct
    open(unit=22,file='brams.log',position='append',action='write')
    write(unit=22,fmt="(a)") ' === Initial timestep info: ngrid, nndtrat, nnacoust,'//&
         '     dtlongn, sscourn,  sspct ==='

#ifdef color
    do ifm = 1,ngrids
      if(dyncore_flag>1) then
        if(sscourn(ifm)<=2.0) then
          write(unit=22,fmt="(29x,i3,1x,i8,1x,i9,f13.3,1x,A,f8.3,A,f8.3)") &
            ifm, nndtrat(ifm), nnacoust(ifm), dtlongn(ifm),c_green,sscourn(ifm),c_noColor,sspct
        else if(sscourn(ifm)>2.0 .and. sscourn(ifm)<=4.0) Then
          write(unit=22,fmt="(29x,i3,1x,i8,1x,i9,f13.3,1x,A,f8.3,A,f8.3)") &
            ifm, nndtrat(ifm), nnacoust(ifm), dtlongn(ifm),c_lightYellow,sscourn(ifm),c_noColor,sspct
        else
          write(unit=22,fmt="(29x,i3,1x,i8,1x,i9,f13.3,1x,A,f8.3,A,f8.3)") &
            ifm, nndtrat(ifm), nnacoust(ifm), dtlongn(ifm), c_red,sscourn(ifm),c_noColor ,sspct
        endif
      else
        if(sscourn(ifm)<=1.0) then
          write(unit=22,fmt="(29x,i3,1x,i8,1x,i9,f13.3,1x,A,f8.3,A,f8.3)") &
            ifm, nndtrat(ifm), nnacoust(ifm), dtlongn(ifm),c_green,sscourn(ifm),c_noColor,sspct
        else if(sscourn(ifm)>1.0 .and. sscourn(ifm)<=2.0) Then
            write(unit=22,fmt="(29x,i3,1x,i8,1x,i9,f13.3,1x,A,f8.3,A,f8.3)") &
              ifm, nndtrat(ifm), nnacoust(ifm), dtlongn(ifm),c_lightYellow,sscourn(ifm),c_noColor,sspct
        else
          write(unit=22,fmt="(29x,i3,1x,i8,1x,i9,f13.3,1x,A,f8.3,A,f8.3)") &
            ifm, nndtrat(ifm), nnacoust(ifm), dtlongn(ifm), c_red,sscourn(ifm),c_noColor ,sspct
        endif
      endif
    enddo
#else
  do ifm = 1,ngrids
    if(dyncore_flag>1) then
      if(sscourn(ifm)<=2.0) then
        write(unit=22,fmt="(29x,i3,1x,i8,1x,i9,f13.3,1x,A,f8.3,A,f8.3)") &
          ifm, nndtrat(ifm), nnacoust(ifm), dtlongn(ifm),'',sscourn(ifm),'',sspct
        else if(sscourn(ifm)>2.0 .and. sscourn(ifm)<=4.0) Then
          write(unit=22,fmt="(29x,i3,1x,i8,1x,i9,f13.3,1x,A,f8.3,A,f8.3)") &
          ifm, nndtrat(ifm), nnacoust(ifm), dtlongn(ifm),'',sscourn(ifm),'',sspct
        else
          write(unit=22,fmt="(29x,i3,1x,i8,1x,i9,f13.3,1x,A,f8.3,A,f8.3)") &
          ifm, nndtrat(ifm), nnacoust(ifm), dtlongn(ifm), '',sscourn(ifm),'' ,sspct
        endif
      else
        if(sscourn(ifm)<=1.0) then
          write(unit=22,fmt="(29x,i3,1x,i8,1x,i9,f13.3,1x,A,f8.3,A,f8.3)") &
          ifm, nndtrat(ifm), nnacoust(ifm), dtlongn(ifm),'',sscourn(ifm),'',sspct
        else if(sscourn(ifm)>1.0 .and. sscourn(ifm)<=2.0) Then
          write(unit=22,fmt="(29x,i3,1x,i8,1x,i9,f13.3,1x,A,f8.3,A,f8.3)") &
            ifm, nndtrat(ifm), nnacoust(ifm), dtlongn(ifm),'',sscourn(ifm),'',sspct
          else
            write(unit=22,fmt="(29x,i3,1x,i8,1x,i9,f13.3,1x,A,f8.3,A,f8.3)") &
            ifm, nndtrat(ifm), nnacoust(ifm), dtlongn(ifm), '',sscourn(ifm),'' ,sspct
          endif
        endif
      enddo
#endif
    close(unit=22)
    if(dyncore_flag == 2) then
       if(sspct<1.0) err=dumpMessage(c_tty,c_yes,'dumpCourn','local_proc',c_notice,'for RK3, sspct is not allowed to be less than 1:',sspct,'I2')
    endif


  end subroutine dump_courn

!=================================================================




  subroutine SetDt(mynum, nndtflg, dxtmax_local)

    use mem_grid, only : &
         ideltat,        & ! INTENT(IN)
         ngrids,         & ! INTENT(IN)
         deltaxn,        & ! INTENT(IN)
         nnxp,           & ! INTENT(IN)
         nnyp,           & ! INTENT(IN)
         nnzp,           & ! INTENT(IN)
         grid_g,         & ! INTENT(IN)
         sspct,          & ! INTENT(INOUT) ! Calculada localmente nesta rotina
         nxtnest,        & ! INTENT(IN)
         dtlongn,        & ! INTENT(INOUT) ! Calculada localmente nesta rotina
         dtlong,         & ! INTENT(IN)    ! Falta passar para o escravo-ok
         nndtrat,        & ! INTENT(INOUT) ! Valor inicial passado-ok
         ! Calculado atualizacao nesta rot.
         ! Retirar da com.:Mest->escravo
         nnacoust,       & ! INTENT(INOUT) ! Calculada localmente nesta rotina
         nacoust,        & ! INTENT(IN)
         cflxy,          & ! INTENT(IN)    ! Calculado localmente em CFL
         cflz,           & ! INTENT(IN)    ! A ser calculado em modsched local
         iflag,          & ! INTENT(INOUT) ! definida localmente nesta rotina
         timmax,         & ! INTENT(IN)
         time,           & ! INTENT(IN)
         dyncore_flag

    ! "cflxy" e "cflz" calculado localmente, precisa ser passado a todos os
    ! escravos para se determinar o maior valor.

    use rconstants, only : &
         cp,               & ! INTENT(IN)
         cv,               & ! INTENT(IN)
         rgas                ! INTENT(IN)

    use ref_sounding, only : &
         th01dn,             & ! INTENT(IN)
         pi01dn                ! INTENT(IN)
    ! Receber apos atualizacao-chamada de varf_read em RAMS_OUTPUT

    use grid_dims, only : &
         maxgrds,         & ! INTENT(IN)
         nzpmax             ! INTENT(IN)

    use io_params, only : &
         frqanl             ! INTENT(IN)


    use node_mod, only: &
         nodemxp, nodemyp

    use mem_stilt, only : &
         iexev             ! INTENT(IN)

    implicit none

    ! Arguments:
    integer, intent(in)  :: mynum
    integer, intent(out) :: nndtflg
    real, intent(in)     :: dxtmax_local(maxgrds)

    ! Local variables:
    integer, parameter  :: ndx=37,ndt=42
    integer, dimension(maxgrds) :: idelt,nndtrat1,nnacoust1
    real, dimension(maxgrds) :: dtlongn1
    real, dimension(nzpmax) :: vctr1
    real :: cflnumh, cflnumv, delx(ndx), delt(ndt)

    integer :: iabsdt,ifm,id,n2,n3,k,nn2,nn3,icm,ntf,ii
    real :: ssmax,tmax,dxtmax,sscnmax,sspct0,cflxyz,timeleft

    real :: dxta, dxtb, dxtc, dxtd
    character(len=8) :: c0, c1, c2
    character(len=*), parameter :: h="**(dtset_new)**"

    delx(:) = (/ &
         200000.,150000.,100000.,80000.,70000.,60000.,40000.,  &
         030000., 20000., 10000., 6000., 4000., 3000., 2000.,  &
         001000.,   800.,   600.,  500.,  400.,  300.,  200.,  &
         000150.,   100.,    80.,   60.,   50.,   40.,   30.,  &
         000020.,    10.,     8.,    6.,    5.,    4.,    3.,  &
         000002.,     1. /)

    delt(:) = (/ &
         300.,  240.,   180.,  150.,  120.,   90.,   60.,  &
         050.,   40.,    30.,   20.,   15.,   12.,   10.,  &
         006.,    5.,     4.,    3.,   2.5,   2.0,   1.5,  &
         01.2,    1.0,     .8,    .6,   .5,    .4,    .3,  &
         00.2,     .1,     .08,   .06,  .05,   .04,   .03,  &
         00.02,    .01,    .008,  .006, .005,  .004,  .003/)

    idelx(1) = 0

    iabsdt = abs(ideltat)

    ! initialize idelx, ssodx, dtlongn, nnacoust, and if required, nndtrat.


    cflnumh = .90
    cflnumv = .90

    do ifm = 1,ngrids

       do id = ndx,1,-1
          if (delx(id)<=deltaxn(ifm))  idelx(ifm) = id
       enddo

       n2 = nodemxp(mynum,ifm)
       n3 = nodemyp(mynum,ifm)
       do k = 1,nnzp(ifm)
          vctr1(k) = th01dn(k,1) * pi01dn(k,1) / cp
       enddo
       tmax = maxval(vctr1(1:nnzp(ifm)))
       ssmax = sqrt(cp / cv * rgas * tmax)

       nn2 = nodemxp(mynum,ifm)
       nn3 = nodemyp(mynum,ifm)

       dxtmax = max(grid_g(ifm)%dxt(1,1), &
            grid_g(ifm)%dxt(nn2,1),       &
            grid_g(ifm)%dxt(nn2,nn3),     &
            grid_g(ifm)%dxt(1,nn3)        )

       ssodx(ifm) = ssmax * dxtmax

    enddo

    if ( ideltat==0 ) then

       sspct = 1.
       do ifm = 1,ngrids
          icm = nxtnest(ifm)
          if ( icm==0 ) then
             dtlongn(ifm) = dtlong
          else
             dtlongn(ifm) = dtlongn(icm) / nndtrat(ifm)
          endif
          nnacoust(ifm) = nacoust
          sscourn(ifm) = 2. * ssodx(ifm) * dtlongn(ifm)
          sspct = min(sspct, .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))
       enddo
       !---srf: not allowing the artifical reduction of the sound speed for RK3
       if(dyncore_flag == 2) sspct = 1.
    else

       sscnmax = 0.
       do ifm = 1,ngrids
          icm = nxtnest(ifm)
          dtlongn(ifm) = delt(idelx(ifm)+iabsdt-1)

          !---srf: max dt allowed by RK3
          if(dyncore_flag == 2) then
	    dxtmax = max(grid_g(ifm)%dxt(1,1), &
                         grid_g(ifm)%dxt(nn2,1),       &
                         grid_g(ifm)%dxt(nn2,nn3),     &
                         grid_g(ifm)%dxt(1,nn3)        )
            dtlongn(ifm)= 6. * (1./dxtmax)/1000.
          endif

          !print*,"dt1=",dtlongn(ifm),icm

          ! For coarse grid(s) adjust dtlongn so that it is an integer
          ! divisor of FRQANL.  For nested grids, compute nndtrat(ifm) as the
          ! first integer greater than or equal to the timestep ratio between
          ! a nested grid's parent and the nested grid. Then adjust
          ! dtlongn(ifm) for the nested grid to be the parent grid timestep
          ! divided by nndtrat(ifm).

          if ( icm==0 ) then
             ntf = nint(frqanl / dtlongn(1))
             dtlongn(ifm) = frqanl / ntf
          else
             nndtrat(ifm) = min(10, nint(dtlongn(icm) / dtlongn(ifm)))
             dtlongn(ifm) = dtlongn(icm) / nndtrat(ifm)
          endif

          ! Compute sst courant numbers (sstcourn(ifm)) for long timestep
          ! dtlongn.

          sscourn(ifm) = 2. * ssodx(ifm) * dtlongn(ifm)
          if ( sscourn(ifm)>sscnmax)  sscnmax = sscourn(ifm)
          !print*,"dt2=",dtlongn(ifm),icm
       enddo

       !---srf: not allowing the artifical reduction of the sound speed for RK3
       if(dyncore_flag >= 2) then

         sspct = 1.
         nnacoust(:) = nacoust

       else

         ! Define trial sspct0 in terms of sscnmax using nonlinear formula
         ! intended to increase nnacoust as sspct decreases, but limit sspct0
         ! to a minimum of .2.

         sspct0 = min(1., (.95/sscnmax)**.5)

         if ( sspct0<.2 ) then
            print*, 'Sound speed percent is forced to be too low'
            stop 'low_sspct0'
         endif

         sspct = 1.
         do ifm = 1,ngrids
            nnacoust(ifm) = max(2, nint(sspct0 * sscourn(ifm) * 2. / .95))
            sspct = min(sspct, .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))
         enddo
       endif
    endif

    ! check Courant numbers

    nndtflg = 0

    if ( ideltat>=0 ) then

       do ifm = 1,ngrids
          cflxyz = max(cflxy(ifm)/cflnumh,cflz(ifm)/cflnumv)
          if ( cflxyz>1. ) then
             iflag = 1
             write(c0,"(f8.2)") cflxyz
             write(c1,"(i8)") mynum
             write(c2,"(i8)") ifm
             call fatal_error(h//" Model will stop because CFL limit exceeded:"//&
                  &" clfxy="//trim(adjustl(c0))//" at proc "//trim(adjustl(c1))//&
                  &" on grid "//trim(adjustl(c2)))
          end if
       enddo
       !print*,"dt4=",dtlongn(1),sspct,cflxyz

    else

       do ifm = 1,ngrids
          icm = nxtnest(ifm)
          cflxyz = max(cflxy(ifm)/cflnumh, cflz(ifm)/cflnumv)

          nndtrat1(ifm) = nndtrat(ifm)
          dtlongn1(ifm) = dtlongn(ifm)
          nnacoust1(ifm) = nnacoust(ifm)

          do id = ndt,1,-1
             if ( delt(id)*cflxyz<=dtlongn(ifm) ) idelt(ifm) = id
          enddo

          if ( idelt(ifm)>idelx(ifm)+4 ) then
             iflag = 1
             write(c1,"(i8)") mynum
             write(c2,"(i8)") ifm
             call fatal_error(h//" Model will stop because adjustable timestep"//&
                  &" forced to become too small at proc "//trim(adjustl(c1))//&
                  &" on grid "//trim(adjustl(c2)))
          else
             ii = max(idelx(ifm)+iabsdt-1,idelt(ifm))
             dtlongn(ifm) = delt(ii)

             ! For the coarse grid(s), adjust dtlongn(1) to not jump past an
             ! analysis write time or the end of the model run time.

             if ( icm==0 ) then
                timeleft =  &
                     min (timmax - time, frqanl - mod(time,frqanl))
                if ( dtlongn(1)>.95 * timeleft) then
                   dtlongn(ifm) = timeleft
                endif
             else
                nndtrat(ifm) = min(10, nint(dtlongn(icm) / dtlongn(ifm)))
                dtlongn(ifm) = dtlongn(icm) / nndtrat(ifm)
             endif
          endif

          ! Compute sst courant numbers (sstcourn(ifm))for long timestep dtlongn.

          sscourn(ifm) = 2. * ssodx(ifm) * dtlongn(ifm)
          if (sscourn(ifm) .gt. sscnmax) sscnmax = sscourn(ifm)
       enddo

       ! Define trial sspct0 in terms of sscnmax using nonlinear formula
       ! intended to increase nnacoust as sspct decreases, but limit sspct0 to
       ! a minimum of .2.

       sspct0 = min(1., (.95/sscnmax) ** .5)
       if ( sspct0<.2 ) then
          print*, 'Sound speed percent is forced to be too low'
          stop 'low_sspct0'
       endif
       !print*,"dt5=",dtlongn(1),sspct,cflxyz

       sspct = 1.
       do ifm = 1,ngrids
          !srf - for�ando nacoust > 2 para iexev=2
          !nnacoust(ifm) = max(2, nint(sspct0 * sscourn(ifm) * 2. / .95))
          if(iexev==1)  nnacoust(ifm) = max(2, nint(sspct0 * sscourn(ifm) * 2. / .95))
          if(iexev==2)  nnacoust(ifm) = max(3, nint(sspct0 * sscourn(ifm) * 2. / .95))
          !srf

          sspct = min(sspct, .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))

          ! If there are any updates to dtlongn(ifm), print out new values.

          if (abs(dtlongn(ifm) - dtlongn1(ifm)) > 1.e-3) then
             write(6,122) ifm,nndtrat(ifm),nnacoust(ifm),dtlongn(ifm)
122          format('Timestep update: ngrid, nndtrat, nnacoust,'  &
                  ,' dtlongn = ',3i3,f10.3)
          endif

          ! If there are any updates to nndtrat(ifm) or others, set
          ! nndtflg = 1 to flag new call to subroutine modsched and send new
          ! stuff to nodes.

          if (nndtrat(ifm) /= nndtrat1(ifm) .or. &
               nnacoust(ifm) /= nnacoust1(ifm) .or. &
               dtlongn(ifm) /= dtlongn1(ifm) ) nndtflg = 1

       enddo

    endif
    !print*,"dt6=",dtlongn(1),sspct,cflxyz
    return
  end subroutine SetDT

  subroutine createDts(frqanl)
    real, intent(in) :: frqanl
    integer :: i,ifrqanl

    allocate(validDts(int(frqanl)))

    ifrqanl=int(frqanl)
    totValidDt=0
    do i=1,ifrqanl
      if(mod(frqanl,real(i))==0.) then
        totValidDt=totValidDt+1
        validDts(totValidDt)=i
        !print *,'Dts: ',i,totValidDt,validDts(totValidDt); call flush(6)
      endif
    enddo

  end subroutine createDts

  subroutine adjustDt(dt)
    real, intent(inout) :: dt
    integer :: i
    real :: dtant

    do i=1,totValidDt
      if(validDts(i)>dt) Then
        dt=dtant
        exit
      endif
      dtant=validDts(i)
    enddo


  end subroutine adjustDt

end module dtset
