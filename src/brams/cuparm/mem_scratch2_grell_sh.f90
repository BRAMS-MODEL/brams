! Module necessary to Grell Cumulus param.
! Scratch variables - module 2

module mem_scratch2_grell_sh

  !TYPE scratch2_grell_vars

  !                   adapted in july-15-2002 for 5.x version
  !
  ! 3d dependence (mgmxp,mgmyp,ensdim)
  real, allocatable, dimension(:,:,:) :: massfln

  ! 2d dependence (mgmxp,mgmzp)
  real, allocatable, dimension(:,:) ::     &
       T,                                  &
       Q,                                  &
       P,                                  &
       PO,                                 &
       TN,                                 &
       QO,                                 &
       OUTT,                               &
       OUTQ,                               &
       outqc,                              &
       US,                                 & ! Substitui US original
       VS,                                 & ! Substitui VS original
       omeg,                               &
       dhdt                                  ! for  bound-layer-quasi-equil closure

  ! 2d dependence (mgmxp, mgmyp)
  integer, allocatable, dimension(:,:) :: KDT, &
       iact_gr,                                &
       iact_old_gr
  real, allocatable, dimension(:,:) :: xland,  &
       massflx,                                &
       tkeg,                                   &
       rcpg


  ! 1d dependence (mgmxp)
  real, allocatable, dimension(:) :: mconv, &
       umean,                               &
       vmean,                               &
       pmean,                               &
       direction
  real, allocatable, dimension(:) :: AA0, &
       PRET,                              &
       PSUR,                              &
       TER11,                             &
       z1
  integer, allocatable, dimension(:) :: KDET

  !END TYPE scratch2_grell_vars

contains

  subroutine alloc_scratch2_grell_sh !(scratch2_grell)

    use mem_grell_param, only : mgmxp,  & ! intent(in)
         mgmyp,                         & ! intent(in)
         mgmzp,                         & ! intent(in)
         ensdim                           ! intent(in)
    use node_mod, only : mynum, &   ! intent(in)
         mxp,                   &   ! intent(in)
         myp,                   &   ! intent(in)
         mzp,                   &   ! intent(in)
         ia,                    &   ! intent(in)
         iz,                    &   ! intent(in)
         ja,                    &   ! intent(in)
         jz,                    &   ! intent(in)
         i0,                    &   ! intent(in)
         j0                         ! intent(in)

    implicit none
    !type (scratch2_grell_vars) :: scratch2_grell

    integer :: i

    !write(6,'(a1,78a1)') ' ',('_',i=1,78)
    !print*,' '
    !print*,'---- In cuparth - alocando memoria--'
    !print*,'------------SCRATCH2_sh-------------'
    !print*,'---------Node configuration---------'
    !print*,'MYNUM_I0____J0____=',mynum,i0,j0
    !print*,'mgmzp_mgmxp_mgmyp_=',mgmzp,mgmxp,mgmyp
    !print*,'m1____m2___m3_____=', mzp, mxp, myp !,m1,m2,m3
    !print*,'ia__iz____ja__jz__=',ia, iz, ja, jz
    !print*,'ensdim_____ialloc_=',ensdim !,ialloc
    !print*,' '
    !write(6,'(a1,78a1)') ' ',('_',i=1,78)

    allocate (massfln(mgmxp, mgmyp, ensdim))

    allocate (mconv    (mgmxp))   ;mconv  = 0.0  
    allocate (umean    (mgmxp))   ;umean  = 0.0   
    allocate (vmean    (mgmxp))   ;vmean  = 0.0   
    allocate (pmean    (mgmxp))   ;pmean  = 0.0   
    allocate (direction(mgmxp))   ;direction= 0.0 
    allocate (aa0      (mgmxp))   ;aa0	 = 0.0    
    allocate (kdet     (mgmxp))   ;kdet  = 0.0    
    allocate (z1       (mgmxp))   ;z1	 = 0.0    
    allocate (pret     (mgmxp))   ;pret  = 0.0   
    allocate (psur     (mgmxp))   ;psur  = 0.0   
    allocate (ter11    (mgmxp))   ;ter11 = 0.0      

    allocate (t          (mgmxp, mgmzp))  ;t	   = 0.0 
    allocate (q          (mgmxp, mgmzp))  ;q	   = 0.0 
    allocate (p          (mgmxp, mgmzp))  ;p	   = 0.0 
    allocate (po         (mgmxp, mgmzp))  ;po	   = 0.0 
    allocate (tn         (mgmxp, mgmzp))  ;tn	   = 0.0 
    allocate (qo         (mgmxp, mgmzp))  ;qo	   = 0.0 
    allocate (outt       (mgmxp, mgmzp))  ;outt    = 0.0 
    allocate (outq       (mgmxp, mgmzp))  ;outq    = 0.0 
    allocate (outqc      (mgmxp, mgmzp))  ;outqc   = 0.0 
    allocate (us         (mgmxp, mgmzp))  ;us	   = 0.0 
    allocate (vs         (mgmxp, mgmzp))  ;vs	   = 0.0 
    allocate (omeg	 (mgmxp, mgmzp))  ;omeg    = 0.0 
    allocate (dhdt(mgmxp, mgmzp))	  ;dhdt    = 0.0 
    allocate (kdt        (mgmxp, mgmyp))  ;kdt     = 0.0 
    allocate (tkeg       (mgmxp, mgmzp))  ;tkeg    = 0.0 
    allocate (rcpg       (mgmxp, mgmzp))  ;rcpg    = 0.0 
    allocate (xland      (mgmxp, mgmyp))  ;xland   = 0.0    
    allocate (massflx    (mgmxp, mgmyp))  ;massflx = 0.0    
    allocate (iact_gr    (mgmxp, mgmyp))  ;iact_gr = 0.0   
    allocate (iact_old_gr(mgmxp, mgmyp))  ;iact_old_gr = 0.0  
					  
    return				  
  end subroutine alloc_scratch2_grell_sh

  subroutine dealloc_scratch2_grell_sh !(scratch2_grell)

    implicit none
    !type (scratch2_grell_vars) :: scratch2_grell

    deallocate (massfln)

    deallocate (mconv)
    deallocate (umean)
    deallocate (vmean)
    deallocate (pmean)
    deallocate (direction)

    deallocate (T)
    deallocate (Q)
    deallocate (P)
    deallocate (PO)
    deallocate (TN)
    deallocate (QO)
    deallocate (OUTT)
    deallocate (OUTQ)
    deallocate (outqc)
    deallocate (PRET)
    deallocate (PSUR)
    deallocate (TER11)
    deallocate (US)
    deallocate (VS)
    deallocate (omeg)
    deallocate (dhdt)
    deallocate (AA0)

    deallocate (KDET)
    deallocate (KDT)

    deallocate (xland)
    deallocate (massflx)
    deallocate (iact_gr)
    deallocate (iact_old_gr)

    return
  end subroutine dealloc_scratch2_grell_sh

!  subroutine zero_scratch2_grell_sh()
!    massfln=0.
!    t=0.
!    q=0.
!    p=0.
!    po=0.
!    tn=0.
!    qo=0.
!    outt=0.
!    outq=0.
!   outqc=0.
!    us=0.
!    vs=0.
!    omeg=0.
!    dhdt=0.
!    xland=0.
!    massflx=0.
!    tkeg=0.
!    rcpg=0.
!    mconv=0.
!    umean=0.
!    vmean=0.
!    pmean=0.
!    direction=0.
!    aa0=0.
!    pret=0.
!    psur=0.
!    ter11=0.
!    z1=0.
!
!    kdet=0
!    kdt=0
!    iact_gr=0
!    iact_old_gr=0
!
!  end subroutine zero_scratch2_grell_sh

!!$  subroutine filltab_scratch2_grell()
!!$
!!$    use var_tables
!!$
!!$    implicit none
!!$
!!$    ! can't think of anything to do here...
!!$
!!$    return
!!$  end subroutine filltab_scratch2_grell

end module mem_scratch2_grell_sh
