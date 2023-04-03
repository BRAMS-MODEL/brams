! Module necessary to Grell Cumulus param.
! Scratch variables - module 3

module mem_scratch3_grell_sh

  !TYPE scratch3_grell_vars
  ! adapted in july-15-2002 for 5.x version
  ! 3d dependence (mgmxp,mgmzp,maxens2)
  real, allocatable, dimension(:,:,:) :: dellat_ens,  &
       dellaq_ens,                                    &
       dellaqc_ens,                                   &
       pwo_ens
  ! 3d dependence (mgmxp,mgmyp,ensdim)
  real, allocatable, dimension(:,:,:) :: xf,  &
       xf_ens,                                &
       pr_ens,                                &
       outt_ens

  ! 2d dependence (mgmxp,mgmzp)
  real, allocatable, dimension(:,:) :: HE,      &
       HES,                                     &
       QES,                                     &
       Z,                                       &
       TV,                                      &
       DBY,                                     &
       QC,                                      &
       QRCD,                                    &
       PWD,                                     &
       PW
  real, allocatable, dimension(:,:) :: HEO,     &
       HESO,                                    &
       QESO,                                    &
       ZO,                                      &
       TVO,                                     &
       DBYO,                                    &
       QCO,                                     &
       QRCDO,                                   &
       PWDO,                                    &
       PWO
  real, allocatable, dimension(:,:) :: XHE,     &
       XHES,                                    &
       XQES,                                    &
       XZ,                                      &
       XTV,                                     &
       XT_Grell,                                & ! Substitui XT original
       XQ,                                      &
       XDBY,                                    &
       XQC,                                     &
       XQRCD,                                   &
       XPWD,                                    &
       XPW
  real, allocatable, dimension(:,:) :: hcd,     &
       hcdo,                                    &
       xhcd
  real, allocatable, dimension(:,:) :: qcd,     &
       qcdo,                                    &
       xqcd
  real, allocatable, dimension(:,:) :: dbyd,    &
       dbydo
  real, allocatable, dimension(:,:) :: hc,      &
       hco,                                     &
       xhc,                                     &
       qrc,                                     &
       qrco,                                    &
       xqrc,                                    &
       zu,                                      &
       zuo,                                     &
       xzu,                                     &
       zd,                                      &
       zdo,                                     &
       xzd
  real, allocatable, dimension(:,:) :: DELLAH,  &
       DELLAQ,                                  &
       DELLAT,                                  &
       DELLAQC
  real, allocatable, dimension(:,:) :: qes_cup,  &
       q_cup,                                    &
       he_cup,                                   &
       hes_cup,                                  &
       z_cup,                                    &
       p_cup,                                    &
       gamma_cup,                                &
       t_cup
  real, allocatable, dimension(:,:) :: qeso_cup,  &
       qo_cup,                                    &
       heo_cup,                                   &
       heso_cup,                                  &
       zo_cup,                                    &
       po_cup,                                    &
       gammao_cup,                                &
       tn_cup
  real, allocatable, dimension(:,:) :: xqes_cup,  &
       xq_cup,                                    &
       xhe_cup,                                   &
       xhes_cup,                                  &
       xz_cup,                                    &
       xt_cup, &
       xt
  real, allocatable, dimension(:,:) :: cd,  &
       cdd,                                 &
       scr1
  ! 2d dependence (mgmxp,maxens)
  real, allocatable, dimension(:,:) :: xaa0_ens
  ! 2d dependence (mgmxp,maxens2)
  real, allocatable, dimension(:,:) :: edtc
  ! 2d dependence (mgmxp,mgmyp)
  real, allocatable, dimension(:,:) :: cwf,  &
       pwf,                                  &
       pwdf,                                 &
       eddt,                                 &
       predb,                                &
       xmass

  ! 1d dependence (mgmxp)
  integer, allocatable, dimension(:) :: kzdown,  &
       KBMAX,                                    &
       IERR,                                     &
       K22,                                      &
       KBCON,                                    &
       KB,                                       &
       JMIN,                                     &
       KTOP,                                     &
       kstabi,                                   &
       kstabm,                                   &
       K22x,                                     &
       KBCONx,                                   &
       KBx,                                      &
       KTOPx, &
       kzi
  real, allocatable, dimension(:) :: EDT,  &
       EDTO,                               &
       EDTX,                               &
       AA1,                                &
       AA0,                                &
       XAA0,                               &
       HKB,                                & 
       HKBO,                               &
       aad,                                &
       XHKB,                               &
       QKB,                                &
       QKBO,                               &
       XMB,                                &
       !for CATT
       tkemax
  real, allocatable, dimension(:) :: XPWAV,  &
       XPWEV,                                &
       PWAV,                                 &
       PWEV,                                 &
       PWAVO,                                &
       PWEVO,                                &
       BU,                                   &
       BUO
  ! 1d dependence (maxens3)
  real, allocatable, dimension(:) :: xff_ens3
  ! 1d dependence (maxens)
  real, allocatable, dimension(:) :: xk
  ! 1d dependence (mgmxp)
  real, allocatable, dimension(:) :: xfac1
  real, allocatable, dimension(:) :: cap_max,cap_max_increment
  ! 1d dependence (maxens)
  real, allocatable, dimension(:) :: mbdt_ens
  ! 1d dependence (maxens2)
  real, allocatable, dimension(:) :: edt_ens
  ! 1d dependence (mgmxp)
  !srf variaveis para a rotine "cup_dd_edt"
  real, allocatable, dimension(:) :: vshear,   &
       sdp,                                    &
       vws

  !END TYPE scratch3_grell_vars

contains

  subroutine alloc_scratch3_grell_sh  !(scratch3_grell)

    use mem_grell_param, only : mgmxp,  & ! INTENT(IN)
         mgmyp,                         & ! INTENT(IN)
         mgmzp,                         & ! INTENT(IN)
         maxens,                        & ! INTENT(IN)
         maxens2,                       & ! INTENT(IN)
         maxens3,                       & ! INTENT(IN)
         ensdim                           ! INTENT(IN)

    implicit none
    !TYPE (scratch3_grell_vars) :: scratch3_grell

    integer :: i

    !write(6,'(a1,78a1)') ' ',('_',i=1,78)
    !print*,'In CUPENSS - Alocando memoria'
    !print*,'------SCRATCH3_GRELL---------'
    !tmp      print*,'mgmxp mgmyp mgmzp ensdim ialloc'
    !tmp      print*,mgmxp,mgmyp,mgmzp,ensdim,ialloc
    !tmp      print*,'istart iend mix mjx mkx'
    !tmp      print*,istart,iend,mix,mjx,mkx
    !write(6,'(a1,78a1)') ' ',('_',i=1,78)


    ! 3d dependence (mgmxp,mgmzp,maxens2)
    allocate (dellat_ens (mgmxp,mgmzp,maxens2))  ;dellat_ens  =0.0
    allocate (dellaq_ens (mgmxp,mgmzp,maxens2))  ;dellaq_ens  =0.0 
    allocate (dellaqc_ens(mgmxp,mgmzp,maxens2))  ;dellaqc_ens =0.0
    allocate (pwo_ens    (mgmxp,mgmzp,maxens2))  ;pwo_ens     =0.0   

    ! 3d dependence (mgmxp,mgmyp,ensdim)
    allocate (xf      (mgmxp,mgmyp,ensdim))  ;xf         =0.0
    allocate (xf_ens  (mgmxp,mgmyp,ensdim))  ;xf_ens  	 =0.0
    allocate (pr_ens  (mgmxp,mgmyp,ensdim))  ;pr_ens  	 =0.0
    allocate (outt_ens(mgmxp,mgmyp,ensdim))  ;outt_ens	 =0.0

    ! 2d dependence   (mgmxp,mgmzp)
    allocate (HE      (mgmxp,mgmzp))   ;HE	= 0.0
    allocate (HES     (mgmxp,mgmzp))   ;HES	= 0.0
    allocate (QES     (mgmxp,mgmzp))   ;QES	= 0.0
    allocate (Z       (mgmxp,mgmzp))   ;Z	= 0.0
    allocate (TV      (mgmxp,mgmzp))   ;TV	= 0.0
    allocate (DBY     (mgmxp,mgmzp))   ;DBY	= 0.0
    allocate (QC      (mgmxp,mgmzp))   ;QC	= 0.0
    allocate (QRCD    (mgmxp,mgmzp))   ;QRCD	= 0.0
    allocate (PWD     (mgmxp,mgmzp))   ;PWD	= 0.0
    allocate (PW      (mgmxp,mgmzp))   ;PW	= 0.0
    allocate (HEO     (mgmxp,mgmzp))   ;HEO	= 0.0
    allocate (HESO    (mgmxp,mgmzp))   ;HESO	= 0.0
    allocate (QESO    (mgmxp,mgmzp))   ;QESO	= 0.0
    allocate (ZO      (mgmxp,mgmzp))   ;ZO	= 0.0
    allocate (TVO     (mgmxp,mgmzp))   ;TVO	= 0.0
    allocate (DBYO    (mgmxp,mgmzp))   ;DBYO	= 0.0
    allocate (QCO     (mgmxp,mgmzp))   ;QCO	= 0.0
    allocate (QRCDO   (mgmxp,mgmzp))   ;QRCDO	= 0.0
    allocate (PWDO    (mgmxp,mgmzp))   ;PWDO	= 0.0
    allocate (PWO     (mgmxp,mgmzp))   ;PWO	= 0.0
    allocate (XHE     (mgmxp,mgmzp))   ;XHE	= 0.0
    allocate (XHES    (mgmxp,mgmzp))   ;XHES	= 0.0
    allocate (XQES    (mgmxp,mgmzp))   ;XQES	= 0.0
    allocate (XZ      (mgmxp,mgmzp))   ;XZ	= 0.0
    allocate (XTV     (mgmxp,mgmzp))   ;XTV	= 0.0
    allocate (XT_Grell(mgmxp,mgmzp))   ;XT_Grell= 0.0
    allocate (XQ      (mgmxp,mgmzp))   ;XQ	= 0.0
    allocate (XDBY    (mgmxp,mgmzp))   ;XDBY	= 0.0
    allocate (XQC     (mgmxp,mgmzp))   ;XQC	= 0.0
    allocate (XQRCD   (mgmxp,mgmzp))   ;XQRCD	= 0.0
    allocate (XPWD    (mgmxp,mgmzp))   ;XPWD	= 0.0
    allocate (XPW     (mgmxp,mgmzp))   ;XPW	= 0.0
    allocate (hcd     (mgmxp,mgmzp))   ;hcd	= 0.0
    allocate (hcdo    (mgmxp,mgmzp))   ;hcdo	= 0.0
    allocate (xhcd    (mgmxp,mgmzp))   ;xhcd	= 0.0
    allocate (qcd     (mgmxp,mgmzp))   ;qcd	= 0.0
    allocate (qcdo    (mgmxp,mgmzp))   ;qcdo	= 0.0
    allocate (xqcd    (mgmxp,mgmzp))   ;xqcd	= 0.0
    allocate (dbyd    (mgmxp,mgmzp))   ;dbyd	= 0.0
    allocate (dbydo   (mgmxp,mgmzp))   ;dbydo	= 0.0
    allocate (hc      (mgmxp,mgmzp))   ;hc	= 0.0
    allocate (hco     (mgmxp,mgmzp))   ;hco	= 0.0
    allocate (xhc     (mgmxp,mgmzp))   ;xhc	= 0.0
    allocate (qrc     (mgmxp,mgmzp))   ;qrc	= 0.0
    allocate (qrco    (mgmxp,mgmzp))   ;qrco	= 0.0
    allocate (xqrc    (mgmxp,mgmzp))   ;xqrc	= 0.0
    allocate (zu      (mgmxp,mgmzp))   ;zu	= 0.0
    allocate (zuo     (mgmxp,mgmzp))   ;zuo	= 0.0
    allocate (xzu     (mgmxp,mgmzp))   ;xzu	= 0.0
    allocate (zd      (mgmxp,mgmzp))   ;zd	= 0.0
    allocate (zdo     (mgmxp,mgmzp))   ;zdo	= 0.0
    allocate (xzd     (mgmxp,mgmzp))   ;xzd	= 0.0
    allocate (DELLAH  (mgmxp,mgmzp))   ;DELLAH  = 0.0
    allocate (DELLAQ  (mgmxp,mgmzp))   ;DELLAQ  = 0.0
    allocate (DELLAT  (mgmxp,mgmzp))   ;DELLAT  = 0.0
    allocate (DELLAQC (mgmxp,mgmzp))   ;DELLAQC = 0.0
    allocate (qes_cup (mgmxp,mgmzp))   ;qes_cup = 0.0
    allocate (q_cup   (mgmxp,mgmzp))   ;q_cup	  = 0.0
    allocate (he_cup  (mgmxp,mgmzp))   ;he_cup    = 0.0
    allocate (hes_cup (mgmxp,mgmzp))   ;hes_cup   = 0.0
    allocate (z_cup   (mgmxp,mgmzp))   ;z_cup	  = 0.0
    allocate (p_cup   (mgmxp,mgmzp))   ;p_cup	  = 0.0
    allocate (gamma_cup(mgmxp,mgmzp))  ;gamma_cup = 0.0
    allocate (t_cup   (mgmxp,mgmzp))   ;t_cup	  = 0.0
    allocate (qeso_cup(mgmxp,mgmzp))   ;qeso_cup  = 0.0
    allocate (qo_cup  (mgmxp,mgmzp))   ;qo_cup    = 0.0
    allocate (heo_cup (mgmxp,mgmzp))   ;heo_cup   = 0.0
    allocate (heso_cup(mgmxp,mgmzp))   ;heso_cup  = 0.0
    allocate (zo_cup  (mgmxp,mgmzp))   ;zo_cup    = 0.0
    allocate (po_cup  (mgmxp,mgmzp))   ;po_cup    = 0.0
    allocate (gammao_cup(mgmxp,mgmzp)) ;gammao_cup= 0.0
    allocate (tn_cup  (mgmxp,mgmzp))   ;tn_cup    = 0.0
    allocate (xqes_cup(mgmxp,mgmzp))   ;xqes_cup  = 0.0
    allocate (xq_cup  (mgmxp,mgmzp))   ;xq_cup    = 0.0
    allocate (xhe_cup (mgmxp,mgmzp))   ;xhe_cup   = 0.0
    allocate (xhes_cup(mgmxp,mgmzp))   ;xhes_cup  = 0.0
    allocate (xz_cup  (mgmxp,mgmzp))   ;xz_cup  = 0.0
    allocate (xt_cup  (mgmxp,mgmzp))   ;xt_cup  = 0.0
    allocate (xt      (mgmxp,mgmzp))   ;xt	= 0.0
    allocate (cd      (mgmxp,mgmzp))   ;cd	= 0.0
    allocate (cdd     (mgmxp,mgmzp))   ;cdd     = 0.0
    allocate (scr1    (mgmxp,mgmzp))   ;scr1    = 0.0
						
    ! 2d dependence   (mgmxp,maxens)		
    allocate (xaa0_ens(mgmxp,maxens))	;xaa0_ens = 0.0
						
    ! 2d dependence (mgmxp,maxens2)		
    allocate (edtc  (mgmxp,maxens2))	;edtc=0.0	
						
    ! 2d dependence(mgmxp,mgmyp)		
    allocate (cwf  (mgmxp,mgmyp))	;cwf	 = 0.0        
    allocate (pwf  (mgmxp,mgmyp))	;pwf	 = 0.0     
    allocate (pwdf (mgmxp,mgmyp))	;pwdf	 = 0.0        
    allocate (eddt (mgmxp,mgmyp))	;eddt	 = 0.0     
    allocate (predb(mgmxp,mgmyp))	;predb   = 0.0     
    allocate (xmass(mgmxp,mgmyp))	;xmass   = 0.0

    ! 1d dependence  (mgmxp)
    allocate (kzdown (mgmxp))     ;kzdown = 0.0
    allocate (KBMAX  (mgmxp))	  ;KBMAX  = 0.0
    allocate (IERR   (mgmxp))	  ;IERR   = 0.0
    allocate (K22    (mgmxp))	  ;K22    = 0.0
    allocate (KBCON  (mgmxp))	  ;KBCON  = 0.0
    allocate (KB     (mgmxp))	  ;KB	  = 0.0
    allocate (JMIN   (mgmxp))	  ;JMIN   = 0.0
    allocate (KTOP   (mgmxp))	  ;KTOP   = 0.0
    allocate (kstabi (mgmxp))	  ;kstabi = 0.0
    allocate (kstabm (mgmxp))	  ;kstabm = 0.0
    allocate (K22x   (mgmxp))	  ;K22x   = 0.0
    allocate (KBCONx (mgmxp))	  ;KBCONx = 0.0
    allocate (KBx    (mgmxp))	  ;KBx    = 0.0
    allocate (KTOPx  (mgmxp))	  ;KTOPx  = 0.0
    allocate (kzi    (mgmxp))	  ;kzi    = 0.0
    allocate (EDT    (mgmxp))	  ;EDT    = 0.0
    allocate (EDTO   (mgmxp))	  ;EDTO   = 0.0
    allocate (EDTX   (mgmxp))	  ;EDTX   = 0.0
    allocate (AA1    (mgmxp))	  ;AA1    = 0.0
    allocate (AA0    (mgmxp))	  ;AA0    = 0.0
    allocate (XAA0   (mgmxp))	  ;XAA0   = 0.0
    allocate (HKB    (mgmxp))	  ;HKB    = 0.0
    allocate (HKBO   (mgmxp))	  ;HKBO   = 0.0
    allocate (aad    (mgmxp))	  ;aad    = 0.0
    allocate (XHKB   (mgmxp))	  ;XHKB   = 0.0
    allocate (QKB    (mgmxp))	  ;QKB    = 0.0
    allocate (QKBO   (mgmxp))	  ;QKBO   = 0.0
    allocate (XMB    (mgmxp))	  ;XMB    = 0.0
    allocate (tkemax (mgmxp))	  ;tkemax = 0.0
    allocate (XPWAV  (mgmxp))	  ;XPWAV  = 0.0
    allocate (XPWEV  (mgmxp))	  ;XPWEV  = 0.0
    allocate (PWAV   (mgmxp))	  ;PWAV   = 0.0
    allocate (PWEV   (mgmxp))	  ;PWEV   = 0.0
    allocate (PWAVO  (mgmxp))	  ;PWAVO  = 0.0
    allocate (PWEVO  (mgmxp))	  ;PWEVO  = 0.0
    allocate (BU     (mgmxp))	  ;BU	  = 0.0
    allocate (BUO    (mgmxp))	  ;BUO    = 0.0
    allocate (xfac1  (mgmxp))	  ;xfac1  = 0.0
    allocate (cap_max(mgmxp))	  ;cap_max= 0.0
    allocate (cap_max_increment(mgmxp)) ;cap_max_increment   = 0.0
    allocate (vshear (mgmxp))		;vshear  = 0.0
    allocate (sdp    (mgmxp))		;sdp = 0.0
    allocate (vws    (mgmxp))		;vws = 0.0
					  
					  
    ! 1d dependence (maxens)		  
    allocate (xk      (maxens)) ;xk=0.0	  
    allocate (mbdt_ens(maxens)) ;mbdt_ens=0.0
					  
    ! 1d dependence (maxens2)		  
    allocate (edt_ens(maxens2)) ;edt_ens=0.0	  
					  
    ! 1d dependence (maxens3)		  
    allocate (xff_ens3(maxens3)) ;xff_ens3=0.0

    return
  end subroutine alloc_scratch3_grell_sh

  subroutine dealloc_scratch3_grell_sh !(scratch3_grell)

    implicit none
    !TYPE (scratch3_grell_vars) :: scratch3_grell

    ! 3d dependence (mgmxp,mgmzp,maxens2)
    deallocate (dellat_ens)
    deallocate (dellaq_ens)
    deallocate (dellaqc_ens)
    deallocate (pwo_ens)
    ! 3d dependence (mgmxp,mgmyp,ensdim)
    deallocate (xf)
    deallocate (xf_ens)
    deallocate (pr_ens)
    deallocate (outt_ens)
    deallocate (xf)

    ! 2d dependence (mgmxp,mgmzp)
    deallocate (HE)
    deallocate (HES)
    deallocate (QES)
    deallocate (Z)
    deallocate (TV)
    deallocate (DBY)
    deallocate (QC)
    deallocate (QRCD)
    deallocate (PWD)
    deallocate (PW)
    deallocate (HEO)
    deallocate (HESO)
    deallocate (QESO)
    deallocate (ZO)
    deallocate (TVO)
    deallocate (DBYO)
    deallocate (QCO)
    deallocate (QRCDO)
    deallocate (PWDO)
    deallocate (PWO)
    deallocate (XHE)
    deallocate (XHES)
    deallocate (XQES)
    deallocate (XZ)
    deallocate (XTV)
    deallocate (XT_Grell)
    deallocate (XQ)
    deallocate (XDBY)
    deallocate (XQC)
    deallocate (XQRCD)
    deallocate (XPWD)
    deallocate (XPW)
    deallocate (hcd)
    deallocate (hcdo)
    deallocate (xhcd)
    deallocate (qcd)
    deallocate (qcdo)
    deallocate (xqcd)
    deallocate (dbyd)
    deallocate (dbydo)
    deallocate (hc)
    deallocate (hco)
    deallocate (xhc)
    deallocate (qrc)
    deallocate (qrco)
    deallocate (xqrc)
    deallocate (zu)
    deallocate (zuo)
    deallocate (xzu)
    deallocate (zd)
    deallocate (zdo)
    deallocate (xzd)
    deallocate (DELLAH)
    deallocate (DELLAQ)
    deallocate (DELLAT)
    deallocate (DELLAQC)
    deallocate (qes_cup)
    deallocate (q_cup)
    deallocate (he_cup)
    deallocate (hes_cup)
    deallocate (z_cup)
    deallocate (p_cup)
    deallocate (gamma_cup)
    deallocate (t_cup)
    deallocate (qeso_cup)
    deallocate (qo_cup)
    deallocate (heo_cup)
    deallocate (heso_cup)
    deallocate (zo_cup)
    deallocate (po_cup)
    deallocate (gammao_cup)
    deallocate (tn_cup)
    deallocate (xqes_cup)
    deallocate (xq_cup)
    deallocate (xhe_cup)
    deallocate (xhes_cup)
    deallocate (xz_cup)
    deallocate (xt_cup)
    deallocate (xt)
    deallocate (cd)
    deallocate (cdd)
    deallocate (scr1)
    ! 2d dependence (mgmxp,maxens)
    deallocate (xaa0_ens)
    ! 2d dependence (mgmxp,maxens2)
    deallocate (edtc)
    ! 2d dependence (mgmxp,mgmyp)
    deallocate (cwf)
    deallocate (pwf)
    deallocate (pwdf)
    deallocate (eddt)
    deallocate (predb)
    deallocate (xmass)

    ! 1d dependence (mgmxp)
    deallocate (kzdown)
    deallocate (KBMAX)
    deallocate (IERR)
    deallocate (K22)
    deallocate (KBCON)
    deallocate (KB)
    deallocate (JMIN)
    deallocate (KTOP)
    deallocate (kstabi)
    deallocate (kstabm)
    deallocate (K22x)
    deallocate (KBCONx)
    deallocate (KBx)
    deallocate (KTOPx)
    deallocate (kzi)
    deallocate (EDT)
    deallocate (EDTO)
    deallocate (EDTX)
    deallocate (AA1)
    deallocate (AA0)
    deallocate (XAA0)
    deallocate (HKB)
    deallocate (HKBO)
    deallocate (aad)
    deallocate (XHKB)
    deallocate (QKB)
    deallocate (QKBO)
    deallocate (XMB)
    deallocate (tkemax)
    deallocate (XPWAV)
    deallocate (XPWEV)
    deallocate (PWAV)
    deallocate (PWEV)
    deallocate (PWAVO)
    deallocate (PWEVO)
    deallocate (BU)
    deallocate (BUO)
    ! 1d dependence (maxens3)
    deallocate (xff_ens3)
    ! 1d dependence (maxens)
    deallocate (xk)
    ! 1d dependence 
    deallocate (xfac1)
    deallocate (cap_max, cap_max_increment)
    ! 1d dependence (maxens)
    deallocate (mbdt_ens)
    ! 1d dependence (maxens2)
    deallocate (edt_ens)
    ! 1d dependence 
    !srf variaveis para a rotine "cup_dd_edt"
    deallocate (vshear)
    deallocate (sdp)
    deallocate (vws)

    return
  end subroutine dealloc_scratch3_grell_sh

!  subroutine zero_scratch3_grell_sh()
!
!    dellat_ens=0.
!    dellaq_ens=0.
!    dellaqc_ens=0.
!    pwo_ens=0.
!    xf=0.
!    xf_ens=0.
!    pr_ens=0.
!    outt_ens=0.
!    xf=0.
!    HE=0.
!    HES=0.
!    QES=0.
!    Z=0.
!    TV=0.
!    DBY=0.
!    QC=0.
!    QRCD=0.
!    PWD=0.
!    PW=0.
!    HEO=0.
!    HESO=0.
!    QESO=0.
!    ZO=0.
!    TVO=0.
!    DBYO=0.
!    QCO=0.
!    QRCDO=0.
!    PWDO=0.
!    PWO=0.
!    XHE=0.
!    XHES=0.
!    XQES=0.
!    XZ=0.
!    XTV=0.
!    XT_Grell=0.
!    XQ=0.
!    XDBY=0.
!    XQC=0.
!    XQRCD=0.
!    XPWD=0.
!    XPW=0.
!    hcd=0.
!    hcdo=0.
!    xhcd=0.
!    qcd=0.
!    qcdo=0.
!    xqcd=0.
!    dbyd=0.
!    dbydo=0.
!    hc=0.
!    hco=0.
!    xhc=0.
!    qrc=0.
!    qrco=0.
!    xqrc=0.
!    zu=0.
!    zuo=0.
!    xzu=0.
!    zd=0.
!    zdo=0.
!    xzd=0.
!    DELLAH=0.
!    DELLAQ=0.
!    DELLAT=0.
!    DELLAQC=0.
!    qes_cup=0.
!    q_cup=0.
!    he_cup=0.
!    hes_cup=0.
!    z_cup=0.
!    p_cup=0.
!    gamma_cup=0.
!    t_cup=0.
!    qeso_cup=0.
!    qo_cup=0.
!    heo_cup=0.
!    heso_cup=0.
!    zo_cup=0.
!    po_cup=0.
!    gammao_cup=0.
!    tn_cup=0.
!    xqes_cup=0.
!    xq_cup=0.
!    xhe_cup=0.
!    xhes_cup=0.
!    xz_cup=0.
!    xt_cup=0.
!    xt=0.
!    cd=0.
!    cdd=0.
!    scr1=0.
!    xaa0_ens=0.
!    edtc=0.
!    cwf=0.
!    pwf=0.
!    pwdf=0.
!    eddt=0.
!    predb=0.
!    xmass=0.
!
!    kzdown=0
!    KBMAX=0
!    IERR=0
!    K22=0
!    KBCON=0
!    KB=0
!    JMIN=0
!    KTOP=0
!    kstabi=0
!    kstabm=0
!    K22x=0
!    KBCONx=0
!    KBx=0
!    KTOPx=0
!    kzi=0
!
!    EDT=0.
!    EDTO=0.
!    EDTX=0.
!    AA1=0.
!    AA0=0.
!    XAA0=0.
!    HKB=0.
!    HKBO=0.
!    aad=0.
!    XHKB=0.
!    QKB=0.
!    QKBO=0.
!    XMB=0.
!    tkemax=0.
!    XPWAV=0.
!    XPWEV=0.
!    PWAV=0.
!    PWEV=0.
!    PWAVO=0.
!    PWEVO=0.
!    BU=0.
!    BUO=0.
!    xff_ens3=0.
!    xk=0.
!    xfac1=0.
!    cap_max=0.
!    cap_max_increment=0.
!
!    mbdt_ens=0.
!    edt_ens=0.
!    vshear=0.
!    sdp=0.
!    vws=0.
!
!  end subroutine zero_scratch3_grell_sh

!!$  SUBROUTINE filltab_scratch3_grell(=0.
!!$
!!$    USE var_tables
!!$
!!$    IMPLICIT NONE
!!$
!!$    ! Can't think of anything to do here...
!!$
!!$    RETURN
!!$  END SUBROUTINE filltab_scratch3_grell

end module mem_scratch3_grell_sh
