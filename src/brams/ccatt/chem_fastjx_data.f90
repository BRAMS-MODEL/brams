MODULE jx_data

  !    include 'parm_mie.f'  for fast-JX code v5.3 (prather 6/05)
  !     N_  = no. of levels in Mie scattering arrays
  !         = 2*NC+1 = 4*LPAR + 5 + 2*sum(JADDLV)
  !     M_  = no. of Gauss points used, must = 4 in fast_JX (no option)

  !-----------------------------------------------------------------------
  INTEGER,  PARAMETER ::   n_=501, m_=4
  !-----------------------------------------------------------------------
  !    include 'parm_CTM.f'  for fast-JX code v5.3 (prather 6/05)

  !     I_ = longitude dim of CTM grid
  !     J_ = latitude  dim of CTM grid
  !     L_ = altitude(levels) dim of CTM grid
  !     LWE_ = altitude(level) dim for trop processes (clouds, rain)
  !     JVL_ = vertical(levels) dim for J-values
  !     L2_  = 2*L1_ = 2*L_ + 2 = no. levels in the basic Fast-JX grid (mid-level_
  !     JVN_ =  no. of J-values
  !     W_   = dim = no. of Wavelength bins
  !     X_   = dim = no. of X-section data sets (input data)
  !     A_   = dim = no. of Aerosol/cloud Mie sets (input data)
  !     MX   = no. of aerosol/cloud types supplied from CTM
  !     NTR_ = no. of CTM tracers
  !     SZAMAX    Solar zenith angle cut-off, above which to skip calculation

  !-----------------------------------------------------------------------
  INTEGER, PARAMETER ::  i_=128, j_=64, l_=37, lwe_=37  !for EC T42L37
  INTEGER, PARAMETER ::  jvl_=37, jvn_=64, w_=18, x_=66, a_=40
 ! INTEGER, PARAMETER ::  jvl_=37, jvn_=62, w_=18, x_=64, a_=40
 
  INTEGER, PARAMETER ::  l1_=l_+1, l2_=2*l_+2
  INTEGER, PARAMETER ::  mx=4, ntr_=1
  DOUBLE PRECISION,  PARAMETER ::  szamax=98.0D0
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  ! Old 'cmn_JVdat.f'  for fast-JX code v5.3 (prather 6/05)

  ! NB - ALL of these common variables are set paramters,
  !    They are NOT to be used as variables for a local solution
  !    Thus this entire set is 'in' only after it is initialized
  !-----------------------------------------------------------------------
!srf  DOUBLE PRECISION :: rad,zzht,atau,atau0
  !	RAD	  Radius of Earth (cm)
  !	ZZHT	  Effective scale height above top of atmosphere (cm)
  DOUBLE PRECISION, parameter :: rad= 6375.d5,zzht= 5.d5
  DOUBLE PRECISION :: atau,atau0
  DOUBLE PRECISION :: wbin(w_+1),wl(w_),fl(w_),qo2(w_,3),qo3(w_,3),q1d(w_,3)
  DOUBLE PRECISION :: qqq(w_,2,x_),qrayl(w_+1),tqq(3,x_)
  DOUBLE PRECISION :: waa(5,a_),qaa(5,a_),paa(8,5,a_),raa(5,a_),ssa(5,a_)
  DOUBLE PRECISION :: jfacta(jvn_)
  INTEGER :: jind(jvn_),nratj,njval,nw1,nw2,naa,jtaumx
  CHARACTER (LEN=20) :: titlea(a_)
  CHARACTER (LEN=78) :: title0
  CHARACTER (LEN=7)  :: titlej(x_),titlej2,titlej3,jlabel(jvn_)
  !-----------------------------------------------------------------------
  !    Old'cmn_metdat.f'  for fast-JX code v5.3 (prather 6/05)

  !        needs 'parm_ctm.f' for dimensions
  !        delivers p, T, Surf Albedo, and Optical Depth from CTM to fastJX
  !        >>>>this is for standalone fast-JX ver 5.3   (6.05)
  !-----------------------------------------------------------------------

  DOUBLE PRECISION :: p(i_,j_)         !  Surface pressure
  DOUBLE PRECISION :: t(i_,j_,l_)      !  Temperature profile
  DOUBLE PRECISION :: od(i_,j_,lwe_)   !  Optical Depth profile

  DOUBLE PRECISION :: xgrd(i_)         !  Longitude (midpoint, radians)
  DOUBLE PRECISION :: xdgrd(i_)
  DOUBLE PRECISION :: ygrd(j_)         !  Latitude  (midpoint, radians)
  DOUBLE PRECISION :: ydgrd(j_)
  DOUBLE PRECISION :: etaa(l1_)       !  Eta(a) value for level boundaries
  DOUBLE PRECISION :: etab(l1_)       !  Eta(b) value for level boundaries
  DOUBLE PRECISION :: areaxy(i_,j_)    !  area (m^2)
  INTEGER :: month
  INTEGER :: nslat         ! Latitude(J) index of current column
  INTEGER :: nslon         ! Longitude(I) index of current column

  DOUBLE PRECISION, DIMENSION(i_,j_,l1_) :: tj, dm, do3, zh
  DOUBLE PRECISION, DIMENSION(i_,j_,l1_) :: daer1, daer2, daer3, odcld
  INTEGER,  DIMENSION(i_,j_,l1_) :: naer1, naer2, naer3, ncldx
  DOUBLE PRECISION, DIMENSION(i_,j_)     :: pmean, sa

  DOUBLE PRECISION ::  stt(i_,j_,l_,ntr_)
  DOUBLE PRECISION ::  tref(51,18,12),oref(51,18,12)

  CHARACTER(LEN=10) :: tcname(ntr_)


END MODULE jx_data







