!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
module mem_tend

   type tend_vars

   real, pointer, dimension(:) ::  &
         ut, vt, wt, pt, tht, rtt  &
        ,rct, rrt, rpt, rst, rat, rgt ,rht  &
        ,cct, crt, cpt, cst, cat, cgt ,cht  &
        ,cccnt, cifnt, tket, epst           &
!Glauber Camponogara 2015
        ,rdt, cdt, gccnt, cccmt, gccmt            &
        ,cnm1t, cnm2t, cnm3t, cnm8t, md1nt, md2nt &
        ,salt_filmt, salt_jett, salt_spmt         &
!GC 2015
!-mb/srf for RK time integration
        ,ut_rk,vt_rk,wt_rk   &
        ,pt_rk               &
        ,tht_rk              &
!-srf for ABM3 time integration
        ,ut_past,vt_past     &
        ,wt_past,pt_past     &
	     ,tht_past            &
!srf for transport of cloud fraction
        ,cldfrt
   end type

   type (tend_vars) :: tend

contains
!---------------------------------------------------------------

   subroutine alloc_tend(nmzp,nmxp,nmyp,ngrs,naddsc,proc_type)

   use mem_basic, only: basic_g   ! Data Type INTENT(IN)
   use mem_grid,  only: dyncore_flag
   use mem_scalar, only: scalar_g ! Data Type INTENT(INOUT)
   use mem_micro, only: micro_g   ! Data Type INTENT(IN)
   use mem_turb, only: turb_g     ! Data Type INTENT(IN)
   use ModNamelistFile, only: namelistFile

   ! TEB_SPM
   use teb_spm_start, only: TEB_SPM ! INTENT(IN)
   use mem_gaspart, only:  gaspart_g !Data Type INTENT(INOUT)
   use mem_emiss, only: ichemi, isource ! INTENT(IN)

   implicit none

   ! Arguments:
   integer, intent(in) :: nmzp(:), nmxp(:), nmyp(:)
   integer, intent(in) :: ngrs, proc_type, naddsc

   ! Local Variables:
   integer :: ng, ntpts, nsc

   !         Find the maximum number of grid points needed for any grid.

   if(proc_type==1) then
      ntpts=1
   else
      ntpts=0
      do ng=1,ngrs
         ntpts=max( nmxp(ng)*nmyp(ng)*nmzp(ng),ntpts )
      enddo
   endif

   ! Allocate arrays based on options (if necessary).
   !   Do not need these arrays if it is a master process in a parallel run,
   !      so just allocate to 1 word.

!!!!!  WE ARE ONLY CHECKING GRID 1 !!!!!!!!!
!!!!!    All grids must have same scalars defined !!!!!!!

!--(DMK-LFR NEC-SX6)----------------------------------------------
!   if (associated(basic_g(1)%up))      allocate (tend%ut(ntpts))
!   if (associated(basic_g(1)%vp))      allocate (tend%vt(ntpts))
!   if (associated(basic_g(1)%wp))      allocate (tend%wt(ntpts))
!   if (associated(basic_g(1)%pp))      allocate (tend%pt(ntpts))
!
!   if (associated(basic_g(1)%thp))     allocate (tend%tht(ntpts))
!   if (associated(basic_g(1)%rtp))     allocate (tend%rtt(ntpts))
!   if (associated(micro_g(1)%rcp))     allocate (tend%rct(ntpts))
!   if (associated(micro_g(1)%rrp))     allocate (tend%rrt(ntpts))
!   if (associated(micro_g(1)%rpp))     allocate (tend%rpt(ntpts))
!   if (associated(micro_g(1)%rsp))     allocate (tend%rst(ntpts))
!   if (associated(micro_g(1)%rap))     allocate (tend%rat(ntpts))
!   if (associated(micro_g(1)%rgp))     allocate (tend%rgt(ntpts))
!   if (associated(micro_g(1)%rhp))     allocate (tend%rht(ntpts))
!   if (associated(micro_g(1)%ccp))     allocate (tend%cct(ntpts))
!   if (associated(micro_g(1)%crp))     allocate (tend%crt(ntpts))
!   if (associated(micro_g(1)%cpp))     allocate (tend%cpt(ntpts))
!   if (associated(micro_g(1)%csp))     allocate (tend%cst(ntpts))
!   if (associated(micro_g(1)%cap))     allocate (tend%cat(ntpts))
!   if (associated(micro_g(1)%cgp))     allocate (tend%cgt(ntpts))
!   if (associated(micro_g(1)%chp))     allocate (tend%cht(ntpts))
!   if (associated(micro_g(1)%cccnp))   allocate (tend%cccnt(ntpts))
!   if (associated(micro_g(1)%cifnp))   allocate (tend%cifnt(ntpts))
!   if (associated(turb_g(1)%tkep))     allocate (tend%tket(ntpts))
!   if (associated(turb_g(1)%epsp))     allocate (tend%epst(ntpts))

   if (associated(basic_g(1)%up))      then
   	allocate (tend%ut(ntpts))
	tend%ut = 0.
        if( dyncore_flag==2 .or. dyncore_flag==3) then
    	  allocate (tend%ut_rk(ntpts))
	  tend%ut_rk = 0.
        endif
        if( dyncore_flag==3 ) then
    	  allocate (tend%ut_past(ntpts))
	  tend%ut_past = 0.
        endif
   endif
   if (associated(basic_g(1)%vp))      then
   	allocate (tend%vt(ntpts))
	tend%vt = 0.
        if( dyncore_flag==2 .or. dyncore_flag==3) then
    	  allocate (tend%vt_rk(ntpts))
	  tend%vt_rk = 0.
        endif
        if( dyncore_flag==3 ) then
    	  allocate (tend%vt_past(ntpts))
	  tend%vt_past = 0.
        endif
   endif
   if (associated(basic_g(1)%wp))      then
   	allocate (tend%wt(ntpts))
	tend%wt = 0.
        if( dyncore_flag==2 .or. dyncore_flag==3) then
    	  allocate (tend%wt_rk(ntpts))
	  tend%wt_rk = 0.
        endif
        if( dyncore_flag==3 ) then
    	  allocate (tend%wt_past(ntpts))
	  tend%wt_past = 0.
        endif
   endif
   if (associated(basic_g(1)%pp))      then
   	allocate (tend%pt(ntpts))
	tend%pt = 0.
        if( dyncore_flag==2 .or. dyncore_flag==3) then
    	  allocate (tend%pt_rk(ntpts))
	  tend%pt_rk = 0.
        endif
        if( dyncore_flag==3 ) then
    	  allocate (tend%pt_past(ntpts))
	  tend%pt_past = 0.
        endif
   endif

   if (associated(basic_g(1)%thp))     then
        allocate (tend%tht(ntpts))
	tend%tht = 0.
        if( dyncore_flag==2 .or. dyncore_flag==3) then
    	  allocate (tend%tht_rk(ntpts))
	  tend%tht_rk = 0.
        endif
        if( dyncore_flag==3 ) then
    	  allocate (tend%tht_past(ntpts))
	  tend%tht_past = 0.
        endif
   endif
   if (associated(basic_g(1)%rtp))     then
        allocate (tend%rtt(ntpts))
	tend%rtt = 0.
   endif
   if (associated(micro_g(1)%rcp))     then
   	allocate (tend%rct(ntpts))
	tend%rct = 0.
   endif
   if (associated(micro_g(1)%rrp))     then
   	allocate (tend%rrt(ntpts))
	tend%rrt = 0.
   endif
   if (associated(micro_g(1)%rpp))     then
   	allocate (tend%rpt(ntpts))
	tend%rpt = 0.
   endif
   if (associated(micro_g(1)%rsp))     then
   	allocate (tend%rst(ntpts))
	tend%rst = 0.
   endif
   if (associated(micro_g(1)%rap))     then
   	allocate (tend%rat(ntpts))
	tend%rat = 0.
   endif
   if (associated(micro_g(1)%rgp))     then
   	allocate (tend%rgt(ntpts))
	tend%rgt = 0.
   endif
   if (associated(micro_g(1)%rhp))     then
   	allocate (tend%rht(ntpts))
	tend%rht = 0.
   endif
   if (associated(micro_g(1)%ccp))     then
   	allocate (tend%cct(ntpts))
	tend%cct = 0.
   endif
   if (associated(micro_g(1)%crp))     then
   	allocate (tend%crt(ntpts))
	tend%crt = 0.
   endif
   if (associated(micro_g(1)%cpp))     then
   	allocate (tend%cpt(ntpts))
	tend%cpt = 0.
   endif
   if (associated(micro_g(1)%csp))     then
   	allocate (tend%cst(ntpts))
	tend%cst = 0.
   endif
   if (associated(micro_g(1)%cap))     then
   	allocate (tend%cat(ntpts))
	tend%cat = 0.
   endif
   if (associated(micro_g(1)%cgp))     then
   	allocate (tend%cgt(ntpts))
	tend%cgt = 0.
   endif
   if (associated(micro_g(1)%chp))     then
   	allocate (tend%cht(ntpts))
	tend%cht = 0.
   endif
   if (associated(micro_g(1)%cccnp))   then
   	allocate (tend%cccnt(ntpts))
	tend%cccnt = 0.
   endif
   if (associated(micro_g(1)%cifnp))   then
   	allocate (tend%cifnt(ntpts))
	tend%cifnt = 0.
   endif
   if (associated(micro_g(1)%cldfr))   then
      allocate (tend%cldfrt(ntpts))
   tend%cldfrt = 0.
   endif
   if (associated(turb_g(1)%tkep))     then
   	allocate (tend%tket(ntpts))
	tend%tket = 0.
   endif
   if (associated(turb_g(1)%epsp))     then
   	allocate (tend%epst(ntpts))
	tend%epst = 0.
   endif
!
!-2015- for 2M microphysics (from G. Camponogara)
   if (associated(micro_g(1)%rdp))     then
       allocate (tend%rdt(ntpts))  ;tend%rdt=0.0
   endif

   if (associated(micro_g(1)%cdp))     then
       allocate (tend%cdt(ntpts))  ;tend%cdt=0.0
   endif

   if (associated(micro_g(1)%gccnp))   then
       allocate (tend%gccnt(ntpts))  ;tend%gccnt=0.0
   endif

   if (associated(micro_g(1)%cccmp))   then
       allocate (tend%cccmt(ntpts))  ;tend%cccmt=0.0
   endif

   if (associated(micro_g(1)%gccmp))   then
       allocate (tend%gccmt(ntpts))  ;tend%gccmt=0.0
   endif

   if (associated(micro_g(1)%cnm1p))   then
       allocate (tend%cnm1t(ntpts))  ;tend%cnm1t=0.0
   endif

   if (associated(micro_g(1)%cnm2p))   then
       allocate (tend%cnm2t(ntpts))  ;tend%cnm2t=0.0
   endif
   if (associated(micro_g(1)%cnm3p))   then
       allocate (tend%cnm3t(ntpts))  ;tend%cnm3t=0.0
   endif

   if (associated(micro_g(1)%cnm8p))   then
       allocate (tend%cnm8t(ntpts))  ;tend%cnm8t=0.0
   endif

   if (associated(micro_g(1)%md1np))   then
       allocate (tend%md1nt(ntpts))  ;tend%md1nt=0.0
   endif

   if (associated(micro_g(1)%md2np))   then
       allocate (tend%md2nt(ntpts))  ;tend%md2nt=0.0
   endif

   if (associated(micro_g(1)%salt_filmp)) then
       allocate (tend%salt_filmt(ntpts))  ;tend%salt_filmt=0.0
   endif
   if (associated(micro_g(1)%salt_jetp))  then
       allocate (tend%salt_jett(ntpts))  ;tend%salt_jett=0.0
   endif

   if (associated(micro_g(1)%salt_spmp))  then
       allocate (tend%salt_spmt(ntpts))  ;tend%salt_spmt=0.0
   endif
!GC - 2M microphysics

!--(DMK-LFR NEC-SX6)----------------------------------------------

   ! TEB_SPM
   if (TEB_SPM==1) then
      if(isource==1)then
         if (associated(gaspart_g(1)%pno).and.  &
              (.not.associated(gaspart_g(1)%pnot)))  &
	      allocate (gaspart_g(1)%pnot(ntpts))

         if (associated(gaspart_g(1)%pno2).and.  &
              (.not.associated(gaspart_g(1)%pno2t)))  &
	      allocate (gaspart_g(1)%pno2t(ntpts))

         if (associated(gaspart_g(1)%ppm25).and.  &
              (.not.associated(gaspart_g(1)%ppm25t)))  &
	      allocate (gaspart_g(1)%ppm25t(ntpts))

         if (associated(gaspart_g(1)%pco).and.  &
              (.not.associated(gaspart_g(1)%pcot)))  &
	      allocate (gaspart_g(1)%pcot(ntpts))

         if (associated(gaspart_g(1)%pso2).and.  &
              (.not.associated(gaspart_g(1)%pso2t)))  &
	      allocate (gaspart_g(1)%pso2t(ntpts))

         if (associated(gaspart_g(1)%pso4).and.  &
              (.not.associated(gaspart_g(1)%pso4t)))  &
	      allocate (gaspart_g(1)%pso4t(ntpts))

         if (associated(gaspart_g(1)%paer).and.  &
              (.not.associated(gaspart_g(1)%paert)))  &
	      allocate (gaspart_g(1)%paert(ntpts))

         if (associated(gaspart_g(1)%pvoc).and.  &
              (.not.associated(gaspart_g(1)%pvoct)))  &
	      allocate (gaspart_g(1)%pvoct(ntpts))

         if(ichemi==1)then

            if (associated(gaspart_g(1)%po3).and.  &
                 (.not.associated(gaspart_g(1)%po3t)))  &
                 allocate (gaspart_g(1)%po3t(ntpts))

            if (associated(gaspart_g(1)%prhco).and.  &
                 (.not.associated(gaspart_g(1)%prhcot)))  &
                 allocate (gaspart_g(1)%prhcot(ntpts))

            if (associated(gaspart_g(1)%pho2).and.  &
                 (.not.associated(gaspart_g(1)%pho2t)))  &
                 allocate (gaspart_g(1)%pho2t(ntpts))

            if (associated(gaspart_g(1)%po3p).and.  &
                 (.not.associated(gaspart_g(1)%po3pt)))  &
                 allocate (gaspart_g(1)%po3pt(ntpts))

            if (associated(gaspart_g(1)%po1d).and.  &
                 (.not.associated(gaspart_g(1)%po1dt)))  &
                 allocate (gaspart_g(1)%po1dt(ntpts))

            if (associated(gaspart_g(1)%pho).and.  &
                 (.not.associated(gaspart_g(1)%phot)))  &
                 allocate (gaspart_g(1)%phot(ntpts))

            if (associated(gaspart_g(1)%proo).and.  &
                 (.not.associated(gaspart_g(1)%proot)))  &
                 allocate (gaspart_g(1)%proot(ntpts))
         endif

         do ng=2,ngrs
            gaspart_g(ng)%pnot   => gaspart_g(1)%pnot
            gaspart_g(ng)%pno2t  => gaspart_g(1)%pno2t
            gaspart_g(ng)%ppm25t => gaspart_g(1)%ppm25t
            gaspart_g(ng)%pcot   => gaspart_g(1)%pcot
            gaspart_g(ng)%pso2t  => gaspart_g(1)%pso2t
            gaspart_g(ng)%pso4t  => gaspart_g(1)%pso4t
            gaspart_g(ng)%paert  => gaspart_g(1)%paert
            gaspart_g(ng)%pvoct  => gaspart_g(1)%pvoct
            if(ichemi==1)then
               gaspart_g(ng)%po3t   => gaspart_g(1)%po3t
               gaspart_g(ng)%prhcot => gaspart_g(1)%prhcot
               gaspart_g(ng)%pho2t  => gaspart_g(1)%pho2t
               gaspart_g(ng)%po3pt  => gaspart_g(1)%po3pt
               gaspart_g(ng)%po1dt  => gaspart_g(1)%po1dt
               gaspart_g(ng)%phot   => gaspart_g(1)%phot
               gaspart_g(ng)%proot  => gaspart_g(1)%proot
            endif
         enddo

      endif

   endif
   !

   do nsc=1,naddsc
      if (associated(scalar_g(nsc,1)%sclp).and.  &
           (.not.associated(scalar_g(nsc,1)%sclt)))  &
           allocate (scalar_g(nsc,1)%sclt(ntpts))
      do ng=2,ngrs
         scalar_g(nsc,ng)%sclt => scalar_g(nsc,1)%sclt
      enddo
   enddo

 end subroutine alloc_tend

!---------------------------------------------------------------

   subroutine nullify_tend(naddsc)

   use mem_scalar, only: scalar_g ! Data Type INTENT(INOUT)

   ! TEB_SPM
   use teb_spm_start, only: TEB_SPM ! INTENT(IN)
   use mem_gaspart, only:  gaspart_g !Data Type INTENT(INOUT)
   use mem_emiss, only: ichemi, isource ! INTENT(IN)

   implicit none

   ! Arguments:
   integer, intent(in) :: naddsc

   ! Local Variables:
   integer :: nsc

! Deallocate all tendency arrays

   if (associated(tend%ut))   nullify (tend%ut)
   if (associated(tend%vt))   nullify (tend%vt)
   if (associated(tend%wt))   nullify (tend%wt)
   if (associated(tend%pt))   nullify (tend%pt)
   if (associated(tend%tht))  nullify (tend%tht)
   if (associated(tend%rtt))  nullify (tend%rtt)
   if (associated(tend%rct))  nullify (tend%rct)
   if (associated(tend%rrt))  nullify (tend%rrt)
   if (associated(tend%rpt))  nullify (tend%rpt)
   if (associated(tend%rst))  nullify (tend%rst)
   if (associated(tend%rat))  nullify (tend%rat)
   if (associated(tend%rgt))  nullify (tend%rgt)
   if (associated(tend%rht))  nullify (tend%rht)
   if (associated(tend%cct))  nullify (tend%cct)
   if (associated(tend%crt))  nullify (tend%crt)
   if (associated(tend%cpt))  nullify (tend%cpt)
   if (associated(tend%cst))  nullify (tend%cst)
   if (associated(tend%cat))  nullify (tend%cat)
   if (associated(tend%cgt))  nullify (tend%cgt)
   if (associated(tend%cht))  nullify (tend%cht)
   if (associated(tend%cccnt))nullify (tend%cccnt)
   if (associated(tend%cifnt))nullify (tend%cifnt)
   if (associated(tend%cldfrt))nullify (tend%cldfrt)
   
   if (associated(tend%tket)) nullify (tend%tket)
   if (associated(tend%epst)) nullify (tend%epst)

!-2015- for 2M microphysics (from G. Camponogara)
   if (associated(tend%rdt  )) nullify(tend%rdt)
   if (associated(tend%cdt  )) nullify(tend%cdt)
   if (associated(tend%gccnt)) nullify(tend%gccnt)
   if (associated(tend%cccmt)) nullify(tend%cccmt)
   if (associated(tend%gccmt)) nullify(tend%gccmt)
   if (associated(tend%cnm1t)) nullify(tend%cnm1t)
   if (associated(tend%cnm2t)) nullify(tend%cnm2t)
   if (associated(tend%cnm3t)) nullify(tend%cnm3t)
   if (associated(tend%cnm8t)) nullify(tend%cnm8t)
   if (associated(tend%md1nt)) nullify(tend%md1nt)
   if (associated(tend%md2nt)) nullify(tend%md2nt)
   if (associated(tend%salt_filmt ))  nullify(tend%salt_filmt)
   if (associated(tend%salt_jett  ))  nullify(tend%salt_jett)
   if (associated(tend%salt_spmt  ))  nullify(tend%salt_spmt)
!-2015- for 2M microphysics (from G. Camponogara)
   !- for RK/ABM3 method
   if (associated(tend%ut_rk))   nullify (tend%ut_rk)
   if (associated(tend%vt_rk))   nullify (tend%vt_rk)
   if (associated(tend%wt_rk))   nullify (tend%wt_rk)
   if (associated(tend%pt_rk))   nullify (tend%pt_rk)
   if (associated(tend%tht_rk))  nullify (tend%tht_rk)
   !- for ABM3 method
   if (associated(tend%ut_past))   nullify (tend%ut_past)
   if (associated(tend%vt_past))   nullify (tend%vt_past)
   if (associated(tend%wt_past))   nullify (tend%wt_past)
   if (associated(tend%pt_past))   nullify (tend%pt_past)
   if (associated(tend%tht_past))  nullify (tend%tht_past)

   ! TEB_SPM
   if (TEB_SPM==1) then
      if(isource==1)then
         if (associated(gaspart_g(1)%pnot  )) nullify (gaspart_g(1)%pnot  )
         if (associated(gaspart_g(1)%pno2t )) nullify (gaspart_g(1)%pno2t )
         if (associated(gaspart_g(1)%ppm25t)) nullify (gaspart_g(1)%ppm25t)
         if (associated(gaspart_g(1)%pcot  )) nullify (gaspart_g(1)%pcot  )
         if (associated(gaspart_g(1)%pso2t )) nullify (gaspart_g(1)%pso2t )
         if (associated(gaspart_g(1)%pso4t )) nullify (gaspart_g(1)%pso4t )
         if (associated(gaspart_g(1)%paert )) nullify (gaspart_g(1)%paert )
         if (associated(gaspart_g(1)%pvoct )) nullify (gaspart_g(1)%pvoct )

         if(ichemi==1)then
            if (associated(gaspart_g(1)%po3t  )) nullify (gaspart_g(1)%po3t  )
            if (associated(gaspart_g(1)%prhcot)) nullify (gaspart_g(1)%prhcot)
            if (associated(gaspart_g(1)%pho2t )) nullify (gaspart_g(1)%pho2t )
            if (associated(gaspart_g(1)%po3pt )) nullify (gaspart_g(1)%po3pt )
            if (associated(gaspart_g(1)%po1dt )) nullify (gaspart_g(1)%po1dt )
            if (associated(gaspart_g(1)%phot  )) nullify (gaspart_g(1)%phot  )
            if (associated(gaspart_g(1)%proot )) nullify (gaspart_g(1)%proot )
         endif
      endif
   endif
   !

   do nsc=1,naddsc
      if (associated(scalar_g(nsc,1)%sclt)) nullify (scalar_g(nsc,1)%sclt)
   enddo

   return
   end subroutine
!---------------------------------------------------------------

   subroutine dealloc_tend(naddsc)

   use mem_scalar, only: scalar_g ! Data Type INTENT(INOUT)

   ! TEB_SPM
   use teb_spm_start, only: TEB_SPM ! INTENT(IN)
   use mem_gaspart, only:  gaspart_g ! Data Type INTENT(INOUT)
   use mem_emiss, only: ichemi, isource ! INTENT(IN)

   implicit none

   ! Arguments:
   integer, intent(in) :: naddsc

   ! Local Variables:
   integer :: nsc

! Deallocate all tendency arrays

   if (associated(tend%ut))   deallocate (tend%ut)
   if (associated(tend%vt))   deallocate (tend%vt)
   if (associated(tend%wt))   deallocate (tend%wt)
   if (associated(tend%pt))   deallocate (tend%pt)
   if (associated(tend%tht))  deallocate (tend%tht)
   if (associated(tend%rtt))  deallocate (tend%rtt)
   if (associated(tend%rct))  deallocate (tend%rct)
   if (associated(tend%rrt))  deallocate (tend%rrt)
   if (associated(tend%rpt))  deallocate (tend%rpt)
   if (associated(tend%rst))  deallocate (tend%rst)
   if (associated(tend%rat))  deallocate (tend%rat)
   if (associated(tend%rgt))  deallocate (tend%rgt)
   if (associated(tend%rht))  deallocate (tend%rht)
   if (associated(tend%cct))  deallocate (tend%cct)
   if (associated(tend%crt))  deallocate (tend%crt)
   if (associated(tend%cpt))  deallocate (tend%cpt)
   if (associated(tend%cst))  deallocate (tend%cst)
   if (associated(tend%cat))  deallocate (tend%cat)
   if (associated(tend%cgt))  deallocate (tend%cgt)
   if (associated(tend%cht))  deallocate (tend%cht)
   if (associated(tend%cccnt))deallocate (tend%cccnt)
   if (associated(tend%cifnt))deallocate (tend%cifnt)
   if (associated(tend%cldfrt))deallocate (tend%cldfrt)

   if (associated(tend%tket)) deallocate (tend%tket)
   if (associated(tend%epst)) deallocate (tend%epst)

!-2015- for 2M microphysics (from G. Camponogara)
   if (associated(tend%rdt  )) deallocate(tend%rdt)
   if (associated(tend%cdt  )) deallocate(tend%cdt)
   if (associated(tend%gccnt)) deallocate(tend%gccnt)
   if (associated(tend%cccmt)) deallocate(tend%cccmt)
   if (associated(tend%gccmt)) deallocate(tend%gccmt)
   if (associated(tend%cnm1t)) deallocate(tend%cnm1t)
   if (associated(tend%cnm2t)) deallocate(tend%cnm2t)
   if (associated(tend%cnm3t)) deallocate(tend%cnm3t)
   if (associated(tend%cnm8t)) deallocate(tend%cnm8t)
   if (associated(tend%md1nt)) deallocate(tend%md1nt)
   if (associated(tend%md2nt)) deallocate(tend%md2nt)
   if (associated(tend%salt_filmt ))  deallocate(tend%salt_filmt)
   if (associated(tend%salt_jett  ))  deallocate(tend%salt_jett)
   if (associated(tend%salt_spmt  ))  deallocate(tend%salt_spmt)
!-2015- for 2M microphysics (from G. Camponogara)

   ! TEB_SPM
   if (TEB_SPM==1) then
      if(isource==1)then
         if (associated(gaspart_g(1)%pnot  )) deallocate (gaspart_g(1)%pnot  )
         if (associated(gaspart_g(1)%pno2t )) deallocate (gaspart_g(1)%pno2t )
         if (associated(gaspart_g(1)%ppm25t)) deallocate (gaspart_g(1)%ppm25t)
         if (associated(gaspart_g(1)%pcot  )) deallocate (gaspart_g(1)%pcot  )
         if (associated(gaspart_g(1)%pso2t )) deallocate (gaspart_g(1)%pso2t )
         if (associated(gaspart_g(1)%pso4t )) deallocate (gaspart_g(1)%pso4t )
         if (associated(gaspart_g(1)%paert )) deallocate (gaspart_g(1)%paert )
         if (associated(gaspart_g(1)%pvoct )) deallocate (gaspart_g(1)%pvoct )

         if(ichemi==1)then
            if (associated(gaspart_g(1)%po3t  )) deallocate (gaspart_g(1)%po3t)
            if (associated(gaspart_g(1)%prhcot)) deallocate (gaspart_g(1)%prhcot)
            if (associated(gaspart_g(1)%pho2t )) deallocate (gaspart_g(1)%pho2t)
            if (associated(gaspart_g(1)%po3pt )) deallocate (gaspart_g(1)%po3pt)
            if (associated(gaspart_g(1)%po1dt )) deallocate (gaspart_g(1)%po1dt)
            if (associated(gaspart_g(1)%phot  )) deallocate (gaspart_g(1)%phot)
            if (associated(gaspart_g(1)%proot )) deallocate (gaspart_g(1)%proot)
         endif

      endif
   endif

   do nsc=1,naddsc
      if (associated(scalar_g(nsc,1)%sclt)) deallocate (scalar_g(nsc,1)%sclt)
   enddo
   if (associated(tend%ut_rk))   deallocate (tend%ut_rk)
   if (associated(tend%vt_rk))   deallocate (tend%vt_rk)
   if (associated(tend%wt_rk))   deallocate (tend%wt_rk)
   if (associated(tend%pt_rk))   deallocate (tend%pt_rk)
   if (associated(tend%tht_rk))  deallocate (tend%tht_rk)

   if (associated(tend%ut_past))   deallocate (tend%ut_past)
   if (associated(tend%vt_past))   deallocate (tend%vt_past)
   if (associated(tend%wt_past))   deallocate (tend%wt_past)
   if (associated(tend%pt_past))   deallocate (tend%pt_past)
   if (associated(tend%tht_past))  deallocate (tend%tht_past)

   return
   end subroutine

!---------------------------------------------------------------

   subroutine filltab_tend(basic,micro,turb,scalar, &
        ! TEB_SPM
        gaspart,                                    &
        !
        naddsc,ng)

   use mem_basic, only: basic_vars   !Type
   use mem_micro, only: micro_vars   !Type
   use mem_turb, only: turb_vars     !Type
   use mem_scalar, only: scalar_vars !Type
   use mem_grid,   only: dyncore_flag
   !use var_tables, only:  ! not used

   ! TEB_SPM
   use teb_spm_start, only: TEB_SPM ! INTENT(IN)
   use mem_gaspart, only: gaspart_vars !Type
   use mem_emiss, only: ichemi, isource ! INTENT(IN)

   implicit none

   ! Arguments:
   type (basic_vars), intent(in)   :: basic
   type (micro_vars), intent(in)   :: micro
   type (turb_vars), intent(in)    :: turb
   type (scalar_vars), intent(in)  :: scalar(:)
   type (gaspart_vars), pointer    :: gaspart ! TEB_SPM
   integer, intent(in)             :: naddsc, ng

   ! Local Variables:
   integer :: elements !ALF
   integer :: nsc
   character (len=7) :: sname

   ! Fill pointers to scalar arrays into scalar tables

   !if ( ( dyncore_flag==0) .or. (dyncore_flag==1) ) then
   !- not need to include THC in scalar table, including only THP
   !- will works better for leapfrog and RK schemes
     if (associated(tend%tht)) then
        call vtables_scalar (basic%thp(1,1,1),tend%tht(1),ng,'THP')
        elements = size(tend%tht)
        call vtables_scalar_new (basic%thp,tend%tht,ng,'THP',elements)
     endif
   !else if ( dyncore_flag==2 ) then
   !  if (associated(tend%tht)) then
   !     call vtables_scalar (basic%thc(1,1,1),tend%tht(1),ng,'THC')
   !     elements = size(tend%tht)
   !     call vtables_scalar_new (basic%thc,tend%tht,ng,'THC',elements)
   !  endif
   !end if

   if (associated(tend%rtt)) then
      call vtables_scalar (basic%rtp(1,1,1),tend%rtt(1),ng,'RTP')
      elements = size(tend%rtt)
      call vtables_scalar_new (basic%rtp(1,1,1),tend%rtt(1),ng,'RTP',elements)
   endif

   if (associated(tend%rct)) then
      call vtables_scalar (micro%rcp(1,1,1),tend%rct(1),ng,'RCP')
      elements = size(tend%rct)
      call vtables_scalar_new (micro%rcp(1,1,1),tend%rct(1),ng,'RCP',elements)
   endif

   if (associated(tend%rrt)) then
      call vtables_scalar (micro%rrp(1,1,1),tend%rrt(1),ng,'RRP')
      elements = size(tend%rrt)
      call vtables_scalar_new (micro%rrp(1,1,1),tend%rrt(1),ng,'RRP',elements)
   endif

   if (associated(tend%rpt)) then
      call vtables_scalar (micro%rpp(1,1,1),tend%rpt(1),ng,'RPP')
      elements = size(tend%rpt)
      call vtables_scalar_new (micro%rpp(1,1,1),tend%rpt(1),ng,'RPP',elements)
   endif

   if (associated(tend%rst)) then
      call vtables_scalar (micro%rsp(1,1,1),tend%rst(1),ng,'RSP')
      elements = size(tend%rst)
      call vtables_scalar_new (micro%rsp(1,1,1),tend%rst(1),ng,'RSP',elements)
   endif

   if (associated(tend%rat)) then
      call vtables_scalar (micro%rap(1,1,1),tend%rat(1),ng,'RAP')
      elements = size(tend%rat)
      call vtables_scalar_new (micro%rap(1,1,1),tend%rat(1),ng,'RAP',elements)
   endif

   if (associated(tend%rgt)) then
      call vtables_scalar (micro%rgp(1,1,1),tend%rgt(1),ng,'RGP')
      elements = size(tend%rgt)
      call vtables_scalar_new (micro%rgp(1,1,1),tend%rgt(1),ng,'RGP',elements)
   endif

   if (associated(tend%rht)) then
      call vtables_scalar (micro%rhp(1,1,1),tend%rht(1),ng,'RHP')
      elements = size(tend%rht)
      call vtables_scalar_new (micro%rhp(1,1,1),tend%rht(1),ng,'RHP',elements)
   endif

   if (associated(tend%cct)) then
      call vtables_scalar (micro%ccp(1,1,1),tend%cct(1),ng,'CCP')
      elements = size(tend%cct)
      call vtables_scalar_new (micro%ccp(1,1,1),tend%cct(1),ng,'CCP',elements)
   endif

   if (associated(tend%crt)) then
      call vtables_scalar (micro%crp(1,1,1),tend%crt(1),ng,'CRP')
      elements = size(tend%crt)
      call vtables_scalar_new (micro%crp(1,1,1),tend%crt(1),ng,'CRP',elements)
   endif

   if (associated(tend%cpt)) then
      call vtables_scalar (micro%cpp(1,1,1),tend%cpt(1),ng,'CPP')
      elements = size(tend%cpt)
      call vtables_scalar_new (micro%cpp(1,1,1),tend%cpt(1),ng,'CPP',elements)
   endif

   if (associated(tend%cst)) then
      call vtables_scalar (micro%csp(1,1,1),tend%cst(1),ng,'CSP')
      elements = size(tend%cst)
      call vtables_scalar_new (micro%csp(1,1,1),tend%cst(1),ng,'CSP',elements)
   endif

   if (associated(tend%cat)) then
      call vtables_scalar (micro%cap(1,1,1),tend%cat(1),ng,'CAP')
      elements = size(tend%cat)
      call vtables_scalar_new (micro%cap(1,1,1),tend%cat(1),ng,'CAP',elements)
   endif

   if (associated(tend%cgt)) then
      call vtables_scalar (micro%cgp(1,1,1),tend%cgt(1),ng,'CGP')
      elements = size(tend%cgt)
      call vtables_scalar_new (micro%cgp(1,1,1),tend%cgt(1),ng,'CGP',elements)
   endif

   if (associated(tend%cht)) then
      call vtables_scalar (micro%chp(1,1,1),tend%cht(1),ng,'CHP')
      elements = size(tend%cht)
      call vtables_scalar_new (micro%chp(1,1,1),tend%cht(1),ng,'CHP',elements)
   endif

   if (associated(tend%cccnt)) then
      call vtables_scalar (micro%cccnp(1,1,1),tend%cccnt(1),ng,'CCCNP')
      elements = size(tend%cccnt)
      call vtables_scalar_new (micro%cccnp(1,1,1),tend%cccnt(1),ng,'CCCNP',elements)
   endif

   if (associated(tend%cifnt)) then
      call vtables_scalar (micro%cifnp(1,1,1),tend%cifnt(1),ng,'CIFNP')
      elements = size(tend%cifnt)
      call vtables_scalar_new (micro%cifnp(1,1,1),tend%cifnt(1),ng,'CIFNP',elements)
   endif

  if (associated(tend%cldfrt)) then
      call vtables_scalar (micro%cldfr(1,1,1),tend%cldfrt(1),ng,'CLDFR')
      elements = size(tend%cldfrt)
      call vtables_scalar_new (micro%cldfr(1,1,1),tend%cldfrt(1),ng,'CLDFR',elements)
   endif

   if( associated(tend%tket)) then
      call vtables_scalar (turb%tkep(1,1,1),tend%tket(1),ng,'TKEP')
      elements = size(tend%tket)
      call vtables_scalar_new (turb%tkep(1,1,1),tend%tket(1),ng,'TKEP',elements)
   endif

   if( associated(tend%epst)) then
      call vtables_scalar     (turb%epsp(1,1,1),tend%epst(1),ng,'EPSP')
      elements = size(tend%epst)
      call vtables_scalar_new (turb%epsp(1,1,1),tend%epst(1),ng,'EPSP',elements)
   endif



!-2015- for 2M microphysics (from G. Camponogara)
   if (associated(tend%rdt  )) then
      call vtables_scalar (micro%rdp(1,1,1),tend%rdt(1),ng,'RDP')
      elements = size(tend%rdt)
      call vtables_scalar_new (micro%rdp(1,1,1),tend%rdt(1),ng,'RDP',elements)
   endif
   if (associated(tend%cdt  )) then
      call vtables_scalar (micro%cdp(1,1,1),tend%cdt(1),ng,'CDP')
      elements = size(tend%cdt)
      call vtables_scalar_new (micro%cdp(1,1,1),tend%cdt(1),ng,'CDP',elements)
   endif
   if (associated(tend%gccnt)) then
      call vtables_scalar (micro%gccnp(1,1,1),tend%gccnt(1),ng,'GCCNP')
      elements = size(tend%gccnt)
      call vtables_scalar_new (micro%gccnp(1,1,1),tend%gccnt(1),ng,'GCCNP',elements)
   endif
   if (associated(tend%cccmt))  then
      call vtables_scalar (micro%cccmp(1,1,1),tend%cccmt(1),ng,'CCCMP')
      elements = size(tend%cccmt)
      call vtables_scalar_new (micro%cccmp(1,1,1),tend%cccmt(1),ng,'CCCMP',elements)
   endif
   if (associated(tend%gccmt))  then
      call vtables_scalar (micro%gccmp(1,1,1),tend%gccmt(1),ng,'GCCMP')
      elements = size(tend%gccmt)
      call vtables_scalar_new (micro%gccmp(1,1,1),tend%gccmt(1),ng,'GCCMP',elements)
   endif
   if (associated(tend%cnm1t))  then
      call vtables_scalar (micro%cnm1p(1,1,1),tend%cnm1t(1),ng,'CNM1P')
      elements = size(tend%cnm1t)
      call vtables_scalar_new (micro%cnm1p(1,1,1),tend%cnm1t(1),ng,'CNM1P',elements)
   endif
   if (associated(tend%cnm2t))  then
      call vtables_scalar (micro%cnm2p(1,1,1),tend%cnm2t(1),ng,'CNM2P')
      elements = size(tend%cnm2t)
      call vtables_scalar_new (micro%cnm2p(1,1,1),tend%cnm2t(1),ng,'CNM2P',elements)
   endif
   if (associated(tend%cnm3t))  then
      call vtables_scalar (micro%cnm3p(1,1,1),tend%cnm3t(1),ng,'CNM3P')
      elements = size(tend%cnm3t)
      call vtables_scalar_new (micro%cnm3p(1,1,1),tend%cnm3t(1),ng,'CNM3P',elements)
   endif
   if (associated(tend%cnm8t))  then
      call vtables_scalar (micro%cnm8p(1,1,1),tend%cnm8t(1),ng,'CNM8P')
      elements = size(tend%cnm8t)
      call vtables_scalar_new (micro%cnm8p(1,1,1),tend%cnm8t(1),ng,'CNM8P',elements)
   endif
   if (associated(tend%md1nt))  then
      call vtables_scalar (micro%md1np(1,1,1),tend%md1nt(1),ng,'MD1NP')
      elements = size(tend%md1nt)
      call vtables_scalar_new (micro%md1np(1,1,1),tend%md1nt(1),ng,'MD1NP',elements)
   endif
   if (associated(tend%md2nt))  then
      call vtables_scalar (micro%md2np(1,1,1),tend%md2nt(1),ng,'MD2NP')
      elements = size(tend%md2nt)
      call vtables_scalar_new (micro%md2np(1,1,1),tend%md2nt(1),ng,'MD2NP',elements)
   endif
   if (associated(tend%salt_filmt))  then
      call vtables_scalar (micro%salt_filmp(1,1,1),tend%salt_filmt(1),ng,'SALT_FILMP')
      elements = size(tend%salt_filmt)
      call vtables_scalar_new (micro%salt_filmp(1,1,1),tend%salt_filmt(1),ng,'SALT_FILMP',elements)
   endif
   if (associated(tend%salt_jett))  then
      call vtables_scalar (micro%salt_jetp(1,1,1),tend%salt_jett(1),ng,'SALT_JETP')
      elements = size(tend%salt_jett)
      call vtables_scalar_new (micro%salt_jetp(1,1,1),tend%salt_jett(1),ng,'SALT_JETP',elements)
   endif
   if (associated(tend%salt_spmt))  then
      call vtables_scalar (micro%salt_spmp(1,1,1),tend%salt_spmt(1),ng,'SALT_SPMP')
      elements = size(tend%salt_spmt)
      call vtables_scalar_new (micro%salt_spmp(1,1,1),tend%salt_spmt(1),ng,'SALT_SPMP',elements)
   endif
!-2015- for 2M microphysics (from G. Camponogara)


   ! TEB_SPM
   if (TEB_SPM==1) then
      if(isource==1)then

         if (associated(gaspart%pnot)) then
            call vtables_scalar (gaspart%pno(1,1,1)  &
                 ,gaspart%pnot(1),ng,'PNO')
            elements = size(gaspart%pno)
            call vtables_scalar_new (gaspart%pno(1,1,1)  &
                 ,gaspart%pnot(1),ng,'PNO',elements)
         endif

         if (associated(gaspart%pno2t)) then
            call vtables_scalar (gaspart%pno2(1,1,1)  &
                 ,gaspart%pno2t(1),ng,'PNO2')
            call vtables_scalar_new (gaspart%pno2(1,1,1)  &
                 ,gaspart%pno2t(1),ng,'PNO2',elements)
         endif

         if (associated(gaspart%ppm25t)) then
            call vtables_scalar (gaspart%ppm25(1,1,1)  &
                 ,gaspart%ppm25t(1),ng,'PM25')
            call vtables_scalar_new (gaspart%ppm25(1,1,1)  &
                 ,gaspart%ppm25t(1),ng,'PM25',elements)
         endif

         if (associated(gaspart%pcot)) then
            call vtables_scalar (gaspart%pco(1,1,1)  &
                 ,gaspart%pcot(1),ng,'PCO')
            call vtables_scalar_new (gaspart%pco(1,1,1)  &
                 ,gaspart%pcot(1),ng,'PCO',elements)
         endif

         if (associated(gaspart%pso2t)) then
            call vtables_scalar (gaspart%pso2(1,1,1)  &
                 ,gaspart%pso2t(1),ng,'PSO2')
            call vtables_scalar_new (gaspart%pso2(1,1,1)  &
                 ,gaspart%pso2t(1),ng,'PSO2',elements)
         endif

         if (associated(gaspart%pso4t)) then
            call vtables_scalar (gaspart%pso4(1,1,1)  &
                 ,gaspart%pso4t(1),ng,'PSO4')
            call vtables_scalar_new (gaspart%pso4(1,1,1)  &
                 ,gaspart%pso4t(1),ng,'PSO4',elements)
         endif

         if (associated(gaspart%paert)) then
            call vtables_scalar (gaspart%paer(1,1,1)  &
                 ,gaspart%paert(1),ng,'PAER')
            call vtables_scalar_new (gaspart%paer(1,1,1)  &
                 ,gaspart%paert(1),ng,'PAER',elements)
         endif

         if (associated(gaspart%pvoct)) then
            call vtables_scalar (gaspart%pvoc(1,1,1)  &
                 ,gaspart%pvoct(1),ng,'PVOC')
            call vtables_scalar_new (gaspart%pvoc(1,1,1)  &
                 ,gaspart%pvoct(1),ng,'PVOC',elements)
         endif

         if(ichemi==1) then
            if (associated(gaspart%po3t)) then
               call vtables_scalar (gaspart%po3(1,1,1)  &
                    ,gaspart%po3t(1),ng,'PO3')
               call vtables_scalar_new (gaspart%po3(1,1,1)  &
                    ,gaspart%po3t(1),ng,'PO3',elements)
            endif

            if (associated(gaspart%prhcot)) then
               call vtables_scalar (gaspart%prhco(1,1,1)  &
                    ,gaspart%prhcot(1),ng,'PRHCO')
               call vtables_scalar_new (gaspart%prhco(1,1,1)  &
                    ,gaspart%prhcot(1),ng,'PRHCO',elements)
            endif

            if (associated(gaspart%pho2t)) then
               call vtables_scalar (gaspart%pho2(1,1,1)  &
                    ,gaspart%pho2t(1),ng,'PHO2')
               call vtables_scalar_new (gaspart%pho2(1,1,1)  &
                    ,gaspart%pho2t(1),ng,'PHO2',elements)
            endif

            if (associated(gaspart%po3pt)) then
               call vtables_scalar (gaspart%po3p(1,1,1)  &
                    ,gaspart%po3pt(1),ng,'PO3P')
               call vtables_scalar_new (gaspart%po3p(1,1,1)  &
                    ,gaspart%po3pt(1),ng,'PO3P',elements)
            endif

            if (associated(gaspart%po1dt)) then
               call vtables_scalar (gaspart%po1d(1,1,1)  &
                    ,gaspart%po1dt(1),ng,'PO1D')
               call vtables_scalar_new (gaspart%po1d(1,1,1)  &
                    ,gaspart%po1dt(1),ng,'PO1D',elements)
            endif

            if (associated(gaspart%phot)) then
               call vtables_scalar (gaspart%pho(1,1,1)  &
                    ,gaspart%phot(1),ng,'PHO')
               call vtables_scalar_new (gaspart%pho(1,1,1)  &
                    ,gaspart%phot(1),ng,'PHO',elements)
            endif

            if (associated(gaspart%proot)) then
               call vtables_scalar (gaspart%proo(1,1,1)  &
                    ,gaspart%proot(1),ng,'PROO')
               call vtables_scalar_new (gaspart%proo(1,1,1)  &
                    ,gaspart%proot(1),ng,'PROO',elements)
            endif
         endif

      endif

   endif
   !

   do nsc=1,naddsc
      write(sname,'(a4,i3.3)') 'SCLP',nsc
      if (associated(scalar(nsc)%sclt)) then
         call vtables_scalar (scalar(nsc)%sclp(1,1,1)  &
              ,scalar(nsc)%sclt(1),ng,sname)
         elements = size(scalar(nsc)%sclt)
         call vtables_scalar_new (scalar(nsc)%sclp(1,1,1)  &
              ,scalar(nsc)%sclt(1),ng,sname,elements)
      endif
   enddo

 end subroutine filltab_tend

end module mem_tend
