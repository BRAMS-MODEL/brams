! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE fcdch_mod

USE phi_m_h_vol_mod, ONLY: phi_m_h_vol
USE phi_m_h_mod, ONLY: phi_m_h
USE sea_rough_int_mod, ONLY: sea_rough_int

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FCDCH_MOD'

CONTAINS

!   SUBROUTINE FCDCH-------------------------------------------------

!  Purpose: Calculate surface transfer coefficients at one or more
!           gridpoints.


!  Documentation: UM Documentation Paper No 24, section 8.
!  Code Owner: Please refer to ModuleLeaders.txt
!--------------------------------------------------------------------

SUBROUTINE fcdch (                                                            &
 cor_mo_iter,points,surft_pts,surft_index,pts_index,                          &
 db,vshr,z0m,z0h,zdt,zh,z1_uv,z1_uv_top,z1_tq,z1_tq_top,                      &
 wind_profile_factor,ddmfx,i_surfalg,charnock,                                &
 charnock_w,                                                                  &
 l_vegdrag,canht,lai,                                                         &
 cdv,chv,cdv_std,v_s,v_s_std,recip_l_mo,u_s_std                               &
)

USE atm_fields_bounds_mod, ONLY:tdims
USE theta_field_sizes, ONLY: t_i_length

USE planet_constants_mod, ONLY: g, vkman
USE jules_surface_mod, ONLY:                                                  &
             i_modiscopt, Limit_ObukhovL, ISrfExCnvGust,IP_SrfExWithCnv,      &
             Improve_Initial_Guess, third, beta, beta_cndd, beta_cnv_bl,      &
             cnst_cndd_0, cnst_cndd_1, cnst_cndd_2,                           &
             min_wind, min_ustar, Ri_m
USE jules_sea_seaice_mod, ONLY:                                               &
             ip_ss_surf_div_int,ip_ss_coare_mq,ip_ss_surf_div_int_coupled
USE bl_option_mod, ONLY: off, on

USE can_drag_mod, ONLY:                                                       &
  can_drag_z0, can_drag_phi_m_h, can_drag_phi_m_h_vol

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER ::                                                                    &
 cor_mo_iter          ! IN Switch for MO iteration correction

INTEGER ::                                                                    &
 points                                                                       &
                      ! IN Number of points.
,surft_pts                                                                    &
                      ! IN Number of tile points.
,surft_index(points)                                                          &
                      ! IN Index of tile points.
,pts_index(points) ! IN Index of land points.

INTEGER, INTENT(IN) ::                                                        &
 i_surfalg
                      ! Option for sea surface transfer.
                      ! Set to ip_ss_solid to mark solid surfaces

LOGICAL, INTENT(IN) ::                                                        &
 l_vegdrag
                      ! Option for vegetation canopy drag scheme.

REAL(KIND=real_jlslsm) ::                                                     &
 db(points)                                                                   &
               ! IN Buoyancy difference between surface and lowest
!                    !    temperature and humidity level in the
!                    !    atmosphere (m/s^2).
,vshr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                       ! IN Wind speed difference between the surf
!                    !    the lowest wind level in the atmosphere (m/s).
,z0m(points)                                                                  &
               ! INOUT Roughness length for momentum transport (m).
,z0h(points)                                                                  &
               ! INOUT Roughness length for heat and moisture (m).
,zdt(points)                                                                  &
               ! INOUT Difference between the canopy height and
!                    ! displacement height (m)
,zh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                      &
                       ! IN Depth of boundary layer (m).
,z1_uv(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                       ! IN Height of lowest wind level (m).
,z1_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                       ! IN Height of lowest temperature and
!                    !    humidity level (m).
,wind_profile_factor(points)
!                    ! IN for adjusting the surface transfer
!                    !    coefficients to remove form drag effects.
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  z1_uv_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                     ! Height of top of lowest uv-layer
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  z1_tq_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                     !  Height of top of lowest Tq-layer
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
                    ddmfx(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!                    !    Convective downdraught mass-flux
!                    !    at cloud-base.
REAL(KIND=real_jlslsm), INTENT(IN) :: charnock
!                    ! Prescribed value of Charnock's coefficient
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  charnock_w(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!                    ! Charnock's coefficient from the wave model
REAL(KIND=real_jlslsm), INTENT(IN) :: canht(points)
                     ! Canopy height (m)
REAL(KIND=real_jlslsm), INTENT(IN) :: lai(points)
                     ! Leaf area index
REAL(KIND=real_jlslsm) ::                                                     &
 cdv(points)                                                                  &
               ! OUT Surface transfer coefficient for momentum
!                    !     including orographic form drag (m/s).
,chv(points)                                                                  &
               ! OUT Surface transfer coefficient for
!                    !     heat, moisture & other scalars (m/s).
,cdv_std(points)                                                              &
!                    ! OUT Surface transfer coefficient for momentum
!                    !     excluding orographic form drag (m/s).
,v_s(points)                                                                  &
               ! OUT Surface layer scaling velocity
!                    !     including orographic form drag (m/s).
,v_s_std(points)                                                              &
!                    ! OUT Surface layer scaling velocity
!                    !     excluding orographic form drag (m/s).
,u_s_std(points)                                                              &
               ! OUT Scaling velocity from middle of MO iteration
!                    !     - picked up in error by dust code!
,recip_l_mo(points)
!                    ! OUT Reciprocal of the Monin-Obukhov length
!                    !     (m^-1).

!    Workspace usage----------------------------------------------------

!     Local work arrays.

REAL(KIND=real_jlslsm) ::                                                     &
 phi_m(points)                                                                &
!                    ! Monin-Obukhov stability function for momentum
!                    ! integrated to the model's lowest wind level.
,phi_h(points)
!                    ! Monin-Obukhov stability function for scalars
!                    ! integrated to the model's lowest temperature
!                    ! and humidity level.

REAL(KIND=real_jlslsm) ::                                                     &
        wstrcnvgust(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!                   ! Turbulent velocity scale for convective gustiness
REAL(KIND=real_jlslsm) :: wgst_tmp

! ----------------------------------------------------------------------

!  Define local variables

INTEGER :: i,j,k,l ! Loop counter; horizontal field index.
INTEGER :: it    ! Iteration loop counter.
INTEGER :: n_its ! Number of iterations for Monin-Obukhov length
!                   ! and stability functions.

REAL(KIND=real_jlslsm) ::                                                     &
 b_flux                                                                       &
              ! Surface buoyancy flux over air density.
,u_s2                                                                         &
              ! Iteration surface friction velocity squared
              ! (effective value)
,u_s_std2                                                                     &
              ! Non-effective version of U_S2
,w_s          ! Surface turbulent convective scaling velocity.

! Derived constants for the convective gustiness calculation
REAL(KIND=real_jlslsm) :: cnst_cndd_1a, cnst_cndd_2a

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FCDCH'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! 0. Initialise values
!-----------------------------------------------------------------------
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(l)                   &
!$OMP SHARED(points,cdv,chv,cdv_std,v_s,v_s_std,u_s_std,recip_l_mo,           &
!$OMP        phi_m,phi_h)
DO l = 1, points
  cdv(l)        = 0.0
  chv(l)        = 0.0
  cdv_std(l)    = 0.0
  v_s(l)        = 0.0
  v_s_std(l)    = 0.0
  u_s_std(l)    = 0.0
  recip_l_mo(l) = 0.0
  phi_m(l)      = 0.0
  phi_h(l)      = 0.0
END DO
!$OMP END PARALLEL DO
n_its         = 8  ! Found typically 0.2% from converged value
!-----------------------------------------------------------------------
! 1. Set initial values for the iteration.
!-----------------------------------------------------------------------
IF (cor_mo_iter == Improve_Initial_Guess) THEN

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(k,l,j,i,wgst_tmp,cnst_cndd_1a,cnst_cndd_2a) &
!$OMP& SHARED(surft_pts,surft_index,pts_index,t_i_length,ddmfx,wstrcnvgust,      &
!$OMP& beta_cnv_bl,db,vshr,recip_l_mo,zh,isrfexcnvgust,g) IF(surft_pts>1)
  IF (isrfexcnvgust == ip_srfexwithcnv) THEN
    cnst_cndd_1a = cnst_cndd_1 / g
    cnst_cndd_2a = cnst_cndd_2 / g**2
    !CDIR NODEP
!$OMP DO SCHEDULE(STATIC)
    DO k = 1,surft_pts
      l = surft_index(k)
      j =(pts_index(l) - 1) / t_i_length + 1
      i = pts_index(l) - (j-1) * t_i_length

      ! Compute convective gustiness, which is unchanged during the iteration.
      ! Redelsperger's adjustment is to U_10; here it is scaled to apply to
      ! the friction velocity.
      wgst_tmp = MIN( 0.6, MAX(0.0, ddmfx(i, j) ) )
      wstrcnvgust(i,j) = 0.067 *  LOG ( cnst_cndd_0 +                         &
                                      cnst_cndd_1a * wgst_tmp +               &
                                      cnst_cndd_2a * wgst_tmp ** 2 )
    END DO
!$OMP END DO
  END IF

  !CDIR NODEP
!$OMP DO SCHEDULE(STATIC)
  DO k = 1,surft_pts
    l = surft_index(k)
    j =(pts_index(l) - 1) / t_i_length + 1
    i = pts_index(l) - (j-1) * t_i_length
    IF (db(l)  <   0.0 .AND. vshr(i,j)  <   2.0) THEN
      !-----------------------------------------------------------------------
      !         Start the iteration from the convective limit.
      !-----------------------------------------------------------------------
      recip_l_mo(l) = -vkman /                                                &
                         (beta_cnv_bl * beta_cnv_bl * beta_cnv_bl * zh(i,j))
    ELSE
      !-----------------------------------------------------------------------
      !         Start the iteration from neutral values.
      !-----------------------------------------------------------------------
      recip_l_mo(l) = 0.0
    END IF
  END DO
!$OMP END DO
!$OMP END PARALLEL

  IF (l_vegdrag) THEN
    CALL can_drag_z0(points, surft_pts, surft_index,                          &
                     recip_l_mo, canht, lai,                                  &
                     z0m, z0h, zdt)
    SELECT CASE (i_modiscopt)
    CASE (off)
      CALL can_drag_phi_m_h(points, surft_pts, surft_index, pts_index,        &
                            recip_l_mo, z1_uv, z1_tq, canht, lai, z0m, z0h,   &
                            phi_m, phi_h)
    CASE (on)
      CALL can_drag_phi_m_h_vol(points, surft_pts, surft_index, pts_index,    &
                       recip_l_mo, z1_uv_top, z1_tq_top, canht, lai, z0m, z0h,&
                       phi_m, phi_h)
    END SELECT
  ELSE
    SELECT CASE (i_modiscopt)
    CASE (off)
      CALL phi_m_h ( points,surft_pts,surft_index,pts_index,                  &
                     recip_l_mo,z1_uv,z1_tq,z0m,z0h,                          &
                     phi_m,phi_h)
    CASE (on)
      CALL phi_m_h_vol ( points,surft_pts,surft_index,pts_index,              &
                     recip_l_mo,z1_uv_top,z1_tq_top,z0m,z0h,                  &
                     phi_m,phi_h)
    END SELECT
  END IF

  IF (l_vegdrag) THEN
    !CDIR NODEP
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k,l,j,i)             &
!$OMP& SHARED(surft_pts,surft_index,pts_index,t_i_length,db,vshr,             &
!$OMP& v_s_std,phi_h,zh,v_s,wind_profile_factor,chv,beta_cnv_bl,              &
!$OMP& cdv,cdv_std,phi_m) IF(surft_pts>1)
    DO k = 1,surft_pts
      l = surft_index(k)
      j=(pts_index(l) - 1) / t_i_length + 1
      i = pts_index(l) - (j-1) * t_i_length
      IF (db(l)  <   0.0 .AND. vshr(i,j)  <   2.0) THEN
        !---------------------------------------------------------------------
        !         Start the iteration from the convective limit
        !---------------------------------------------------------------------
        v_s_std(l) = beta_cnv_bl *                                            &
            SQRT( beta_cnv_bl * ( vkman / phi_h(l) ) * zh(i,j) * (-db(l)) )
        !---------------------------------------------------------------------
        !         Start from the neutral limit if only slightly unstable
        !---------------------------------------------------------------------
        v_s_std(l) = MAX(v_s_std(l), ( vkman / phi_m(l) ) * vshr(i,j))
        v_s(l) = v_s_std(l)
      ELSE
        !---------------------------------------------------------------------
        !         Start the iteration from neutral values.
        !---------------------------------------------------------------------
        v_s_std(l) = ( vkman / phi_m(l) ) * vshr(i,j)
        v_s(l) = v_s_std(l) / wind_profile_factor(l)
      END IF
      chv(l) = ( vkman / phi_h(l) ) * v_s_std(l)
      cdv_std(l) = ( vkman / phi_m(l) ) * v_s_std(l)
      cdv(l) = cdv_std(l) * ( v_s(l) / v_s_std(l) ) /                         &
                            wind_profile_factor(l)
    END DO
!$OMP END PARALLEL DO
  ELSE
    !CDIR NODEP
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k,l,j,i)             &
!$OMP& SHARED(surft_pts,surft_index,pts_index,t_i_length,db,vshr,             &
!$OMP& v_s_std,phi_h,zh,v_s,wind_profile_factor,chv,beta_cnv_bl,              &
!$OMP& cdv,cdv_std,phi_m) IF(surft_pts>1)
    DO k = 1,surft_pts
      l = surft_index(k)
      j=(pts_index(l) - 1) / t_i_length + 1
      i = pts_index(l) - (j-1) * t_i_length
      IF (db(l)  <   0.0 .AND. vshr(i,j)  <   2.0) THEN
        !---------------------------------------------------------------------
        !         Start the iteration from the convective limit.
        !---------------------------------------------------------------------
        v_s_std(l) = beta_cnv_bl *                                            &
            SQRT( beta_cnv_bl * ( vkman / phi_h(l) ) * zh(i,j) * (-db(l)) )
        !---------------------------------------------------------------------
        !         Start from the neutral limit if only slightly unstable
        !---------------------------------------------------------------------
        v_s_std(l) = MAX(v_s_std(l), ( vkman / phi_m(l) ) * vshr(i,j))
        v_s(l) = v_s_std(l)
      ELSE
        !---------------------------------------------------------------------
        !         Start the iteration from neutral values.
        !---------------------------------------------------------------------
        v_s(l) = ( vkman / phi_m(l) ) * vshr(i,j)
        v_s_std(l) = v_s(l) * wind_profile_factor(l)
      END IF
      chv(l) = ( vkman / phi_h(l) ) * v_s_std(l)
      cdv(l) = ( vkman / phi_m(l) ) * v_s(l)
      cdv_std(l) = cdv(l) * ( v_s_std(l) / v_s(l) ) *                         &
                          wind_profile_factor(l)
    END DO
!$OMP END PARALLEL DO
  END IF
  !-----------------------------------------------------------------------
  ! 2. Iterate to obtain sucessively better approximations for CD & CH.
  !-----------------------------------------------------------------------
  DO it = 1,n_its
    !
    ! --------------------------------------------------------------
    ! Sea Surface
    ! --------------------------------------------------------------
    ! Modify roughness lengths if required by the surface algorithm.
    SELECT CASE (i_surfalg)
      !
    CASE (ip_ss_surf_div_int, ip_ss_coare_mq, ip_ss_surf_div_int_coupled)
      CALL sea_rough_int (                                                    &
        points,surft_pts,surft_index,pts_index,                               &
        charnock,charnock_w,v_s,recip_l_mo,                                   &
        z0m,z0h                                                               &
        )
      !
    CASE DEFAULT
      !     Do not alter the roughness lengths at this point.
      !
    END SELECT
    !
    ! --------------------------------------------------------------
    !
    !CDIR NODEP
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(k,l,j,i,b_flux,u_s2,w_s,u_s_std2)                              &
!$OMP& SHARED(surft_pts,surft_index,pts_index,t_i_length,chv,db,cdv,          &
!$OMP& vshr,cdv_std,u_s_std,zh,v_s,v_s_std,recip_l_mo,isrfexcnvgust,          &
!$OMP& beta_cnv_bl,wstrcnvgust) IF(surft_pts>1)
    DO k = 1,surft_pts
      l = surft_index(k)
      j=(pts_index(l) - 1) / t_i_length + 1
      i = pts_index(l) - (j-1) * t_i_length
      b_flux = -chv(l) * db(l)

      IF ( vshr(i,j) > min_wind) THEN
        u_s2 = cdv(l) * vshr(i,j)
      ELSE
        u_s2 = min_ustar
      END IF
      u_s_std2 = cdv_std(l) * vshr(i,j)
      u_s_std(l) = SQRT( u_s_std2 )

      IF (db(l)  <   -EPSILON(0.0)) THEN
        w_s = MAX( (zh(i,j) * b_flux)**third, EPSILON(0.0) )
        ! Note that, during this iteration, CDV already includes
        ! this gust enhancement and so U_S2 is not simply the
        ! friction velocity arising from the mean wind, hence:
        SELECT CASE(isrfexcnvgust)
        CASE (off)
          v_s(l) = SQRT( 0.5 * ( beta_cnv_bl * beta_cnv_bl * w_s * w_s        &
                 + SQRT( (beta_cnv_bl * w_s)**4 + 4.0 * u_s2 * u_s2 )) )
          v_s_std(l) = SQRT( 0.5 * ( beta_cnv_bl * beta_cnv_bl * w_s * w_s    &
                     + SQRT( (beta_cnv_bl * w_s)**4                           &
                     + 4.0 * u_s_std2 * u_s_std2 )) )
        CASE (ip_srfexwithcnv)
          v_s(l) = SQRT( 0.5 * ( beta_cnv_bl * beta_cnv_bl * w_s * w_s        &
                 + wstrcnvgust(i,j)**2                                        &
                 + SQRT( ((beta_cnv_bl * w_s)**2                              &
                 + wstrcnvgust(i,j)**2)**2                                    &
                     + 4.0 * u_s2 * u_s2 )) )
          v_s_std(l) = SQRT( 0.5 * ( beta_cnv_bl * beta_cnv_bl * w_s * w_s    &
                 + wstrcnvgust(i,j)**2                                        &
                 + SQRT( ((beta_cnv_bl * w_s)**2                              &
                 + wstrcnvgust(i,j)**2)**2                                    &
                 + 4.0 * u_s_std2 * u_s_std2 )) )
        END SELECT
      ELSE
        v_s(l) = SQRT( u_s2 )
        v_s_std(l) = SQRT( u_s_std2 )
      END IF

      IF ( v_s(l) > min_wind) THEN
        recip_l_mo(l) = -vkman * b_flux                                       &
                      / (v_s(l) * v_s(l) * v_s(l))
      ELSE
        recip_l_mo(l) = 1.0 / SQRT(EPSILON(0.0))
      END IF
    END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k,l,j,i)             &
!$OMP& SHARED(surft_pts,surft_index,pts_index,t_i_length,db,z1_tq,            &
!$OMP&  z1_uv,recip_l_mo,vshr) IF(surft_pts>1)
    DO k = 1,surft_pts
      l = surft_index(k)
      j=(pts_index(l) - 1) / t_i_length + 1
      i = pts_index(l) - (j-1) * t_i_length

      IF ( db(l) > Ri_m * z1_tq(i,j) * (vshr(i,j) / z1_uv(i,j))**2.0 ) THEN
        recip_l_mo(l)= 1.5 * Ri_m**2.0 / z1_tq(i,j)
      END IF
    END DO
!$OMP END PARALLEL DO

    IF (l_vegdrag) THEN
      CALL can_drag_z0(points, surft_pts, surft_index,                        &
                       recip_l_mo, canht, lai,                                &
                       z0m, z0h, zdt)

      SELECT CASE (i_modiscopt)
      CASE (off)
        CALL can_drag_phi_m_h(points, surft_pts, surft_index, pts_index,      &
                              recip_l_mo, z1_uv, z1_tq, canht, lai, z0m, z0h, &
                              phi_m, phi_h)
      CASE (on)
        CALL can_drag_phi_m_h_vol(points, surft_pts, surft_index, pts_index,  &
                       recip_l_mo, z1_uv_top, z1_tq_top, canht, lai, z0m, z0h,&
                       phi_m, phi_h)
      END SELECT

      !CDIR NODEP
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k,l)                 &
!$OMP& SHARED(surft_pts,surft_index,chv,cdv,cdv_std,phi_h,                    &
!$OMP& phi_m,v_s_std,v_s,wind_profile_factor) IF(surft_pts>1)
      DO k = 1,surft_pts
        l = surft_index(k)
        chv(l) = ( vkman / phi_h(l) ) * v_s_std(l)
        cdv_std(l) = ( vkman / phi_m(l) ) * v_s_std(l)
        cdv(l) = cdv_std(l) * ( v_s(l) / v_s_std(l) ) /                       &
                              wind_profile_factor(l)
      END DO
!$OMP END PARALLEL DO
    ELSE
      SELECT CASE (i_modiscopt)
      CASE (off)
        CALL phi_m_h ( points,surft_pts,surft_index,pts_index,                &
                       recip_l_mo,z1_uv,z1_tq,z0m,z0h,                        &
                       phi_m,phi_h)
      CASE (on)
        CALL phi_m_h_vol ( points,surft_pts,surft_index,pts_index,            &
                       recip_l_mo,z1_uv_top,z1_tq_top,z0m,z0h,                &
                       phi_m,phi_h)
      END SELECT


      !CDIR NODEP
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k,l)                 &
!$OMP& SHARED(surft_pts,surft_index,chv,cdv,cdv_std,phi_h,                    &
!$OMP& phi_m,v_s_std,v_s,wind_profile_factor) IF(surft_pts>1)
      DO k = 1,surft_pts
        l = surft_index(k)
        chv(l) = ( vkman / phi_h(l) ) * v_s_std(l)
        cdv(l) = ( vkman / phi_m(l) ) * v_s(l)
        cdv_std(l) = cdv(l) * ( v_s_std(l) / v_s(l) ) *                       &
                              wind_profile_factor(l)
      END DO
!$OMP END PARALLEL DO
    END IF
  END DO ! Iteration loop

ELSE     ! cor_mo_iter earlier option than Improve_Initial_Guess

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(k,l,j,i,wgst_tmp,cnst_cndd_1a,cnst_cndd_2a) &
!$OMP& SHARED(surft_pts,surft_index,pts_index,t_i_length,ddmfx,wstrcnvgust,   &
!$OMP& db,vshr,recip_l_mo,zh,isrfexcnvgust,g) IF(surft_pts>1)
  IF (isrfexcnvgust == ip_srfexwithcnv) THEN
    cnst_cndd_1a = cnst_cndd_1 / g
    cnst_cndd_2a = cnst_cndd_2 / g**2
    !CDIR NODEP
!$OMP DO SCHEDULE(STATIC)
    DO k = 1,surft_pts
      l = surft_index(k)
      j =(pts_index(l) - 1) / t_i_length + 1
      i = pts_index(l) - (j-1) * t_i_length

      !         Compute convective gustiness, which is unchanged during
      !         the iteration. Redelsperger's adjustment is to U_10; here it
      !         is scaled to apply to the friction velocity.
      wgst_tmp = MIN( 0.6, MAX(0.0, ddmfx(i, j) ) )
      wstrcnvgust(i,j) = 0.067 *  LOG ( cnst_cndd_0 +                         &
                                     cnst_cndd_1a * wgst_tmp +                &
                                     cnst_cndd_2a * wgst_tmp ** 2 )
    END DO
!$OMP END DO
  END IF

  !CDIR NODEP
!$OMP DO SCHEDULE(STATIC)
  DO k = 1,surft_pts
    l = surft_index(k)
    j=(pts_index(l) - 1) / t_i_length + 1
    i = pts_index(l) - (j-1) * t_i_length
    IF (db(l)  <   0.0 .AND. vshr(i,j)  <   2.0) THEN
      !-----------------------------------------------------------------------
      !         Start the iteration from the convective limit.
      !-----------------------------------------------------------------------
      recip_l_mo(l) = -vkman / (beta * beta * beta * zh(i,j))
    ELSE
      !-----------------------------------------------------------------------
      !         Start the iteration from neutral values.
      !-----------------------------------------------------------------------
      recip_l_mo(l) = 0.0
    END IF
  END DO
!$OMP END DO
!$OMP END PARALLEL

  IF (l_vegdrag) THEN
    CALL can_drag_z0(points, surft_pts, surft_index,                          &
                     recip_l_mo, canht, lai,                                  &
                     z0m, z0h, zdt)
    SELECT CASE (i_modiscopt)
    CASE (off)
      CALL can_drag_phi_m_h(points, surft_pts, surft_index, pts_index,        &
                            recip_l_mo, z1_uv, z1_tq, canht, lai, z0m, z0h,   &
                            phi_m, phi_h)
    CASE (on)
      CALL can_drag_phi_m_h_vol(points, surft_pts, surft_index, pts_index,    &
                       recip_l_mo, z1_uv_top, z1_tq_top, canht, lai, z0m, z0h,&
                       phi_m, phi_h)
    END SELECT
  ELSE
    SELECT CASE (i_modiscopt)
    CASE (off)
      CALL phi_m_h ( points,surft_pts,surft_index,pts_index,                  &
                     recip_l_mo,z1_uv,z1_tq,z0m,z0h,                          &
                     phi_m,phi_h)
    CASE (on)
      CALL phi_m_h_vol ( points,surft_pts,surft_index,pts_index,              &
                     recip_l_mo,z1_uv_top,z1_tq_top,z0m,z0h,                  &
                     phi_m,phi_h)
    END SELECT
  END IF

  IF (l_vegdrag) THEN
    !CDIR NODEP
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k,l,j,i)             &
!$OMP& SHARED(surft_pts,surft_index,pts_index,t_i_length,db,vshr,             &
!$OMP& v_s_std,phi_h,zh,v_s,wind_profile_factor,chv,                          &
!$OMP& cdv,cdv_std,phi_m) IF(surft_pts>1)
    DO k = 1,surft_pts
      l = surft_index(k)
      j=(pts_index(l) - 1) / t_i_length + 1
      i = pts_index(l) - (j-1) * t_i_length
      IF (db(l)  <   0.0 .AND. vshr(i,j)  <   2.0) THEN
        !---------------------------------------------------------------------
        !         Start the iteration from the convective limit
        !---------------------------------------------------------------------
        v_s_std(l) = beta *                                                   &
            SQRT( beta * ( vkman / phi_h(l) ) * zh(i,j) * (-db(l)) )
        !---------------------------------------------------------------------
        !         Start from the neutral limit if only slightly unstable
        !---------------------------------------------------------------------
        v_s(l) = v_s_std(l)
      ELSE
        !---------------------------------------------------------------------
        !         Start the iteration from neutral values.
        !---------------------------------------------------------------------
        v_s_std(l) = ( vkman / phi_m(l) ) * vshr(i,j)
        v_s(l) = v_s_std(l) / wind_profile_factor(l)
      END IF
      chv(l) = ( vkman / phi_h(l) ) * v_s_std(l)
      cdv_std(l) = ( vkman / phi_m(l) ) * v_s_std(l)
      cdv(l) = cdv_std(l) * ( v_s(l) / v_s_std(l) ) /                         &
                            wind_profile_factor(l)
    END DO
!$OMP END PARALLEL DO
  ELSE
    !CDIR NODEP
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k,l,j,i)             &
!$OMP& SHARED(surft_pts,surft_index,pts_index,t_i_length,db,vshr,             &
!$OMP& v_s_std,phi_h,zh,v_s,wind_profile_factor,chv,                          &
!$OMP& cdv,cdv_std,phi_m) IF(surft_pts>1)
    DO k = 1,surft_pts
      l = surft_index(k)
      j=(pts_index(l) - 1) / t_i_length + 1
      i = pts_index(l) - (j-1) * t_i_length
      IF (db(l)  <   0.0 .AND. vshr(i,j)  <   2.0) THEN
        !---------------------------------------------------------------------
        !         Start the iteration from the convective limit.
        !---------------------------------------------------------------------
        v_s_std(l) = beta *                                                   &
            SQRT( beta * ( vkman / phi_h(l) ) * zh(i,j) * (-db(l)) )
        !---------------------------------------------------------------------
        !         Start from the neutral limit if only slightly unstable
        !---------------------------------------------------------------------
        v_s(l) = v_s_std(l)
      ELSE
        !---------------------------------------------------------------------
        !         Start the iteration from neutral values.
        !---------------------------------------------------------------------
        v_s(l) = ( vkman / phi_m(l) ) * vshr(i,j)
        v_s_std(l) = v_s(l) * wind_profile_factor(l)
      END IF
      chv(l) = ( vkman / phi_h(l) ) * v_s_std(l)
      cdv(l) = ( vkman / phi_m(l) ) * v_s(l)
      cdv_std(l) = cdv(l) * ( v_s_std(l) / v_s(l) ) *                         &
                          wind_profile_factor(l)
    END DO
!$OMP END PARALLEL DO
  END IF
  !-----------------------------------------------------------------------
  ! 2. Iterate to obtain sucessively better approximations for CD & CH.
  !-----------------------------------------------------------------------
  DO it = 1,n_its
    !
    ! --------------------------------------------------------------
    ! Sea Surface
    ! --------------------------------------------------------------
    ! Modify roughness lengths if required by the surface algorithm.
    SELECT CASE (i_surfalg)
      !
    CASE (ip_ss_surf_div_int, ip_ss_coare_mq, ip_ss_surf_div_int_coupled)
      CALL sea_rough_int (                                                    &
        points,surft_pts,surft_index,pts_index,                               &
        charnock,charnock_w,v_s,recip_l_mo,                                   &
        z0m,z0h                                                               &
        )
      !
    CASE DEFAULT
      !     Do not alter the roughness lengths at this point.
      !
    END SELECT
    !
    !-----------------------------------------------------------------------
    !  Corrected version
    !-----------------------------------------------------------------------
    !CDIR NODEP
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(k,l,j,i,b_flux,u_s2,w_s,u_s_std2)                              &
!$OMP& SHARED(surft_pts,surft_index,pts_index,t_i_length,chv,db,cdv,          &
!$OMP& vshr,cdv_std,u_s_std,zh,v_s,v_s_std,recip_l_mo,isrfexcnvgust,          &
!$OMP& wstrcnvgust,cor_mo_iter) IF(surft_pts>1)
    DO k = 1,surft_pts
      l = surft_index(k)
      j=(pts_index(l) - 1) / t_i_length + 1
      i = pts_index(l) - (j-1) * t_i_length
      b_flux = -chv(l) * db(l)

      IF ( vshr(i,j) > min_wind) THEN
        u_s2 = cdv(l) * vshr(i,j)
      ELSE
        u_s2 = min_ustar
      END IF
      u_s_std2 = cdv_std(l) * vshr(i,j)
      u_s_std(l) = SQRT( u_s_std2 )

      IF (db(l)  <   -EPSILON(0.0)) THEN
        w_s = MAX( (zh(i,j) * b_flux)**third, EPSILON(0.0) )
        !           ! Note that, during this iteration, CDV already includes
        !           ! this gust enhancement and so U_S2 is not simply the
        !           ! friction velocity arising from the mean wind, hence:
        SELECT CASE(isrfexcnvgust)
        CASE (off)
          v_s(l) = SQRT( 0.5 * ( beta * beta * w_s * w_s                      &
                 + SQRT( (beta * w_s)**4 + 4.0 * u_s2 * u_s2 )) )
          v_s_std(l) = SQRT( 0.5 * ( beta * beta * w_s * w_s                  &
                     + SQRT( (beta * w_s)**4                                  &
                     + 4.0 * u_s_std2 * u_s_std2 )) )
        CASE (ip_srfexwithcnv)
          v_s(l) = SQRT( 0.5 * ( beta_cndd * beta_cndd * w_s * w_s            &
                 + wstrcnvgust(i,j)**2                                        &
                 + SQRT( ((beta_cndd * w_s)**2                                &
                 + wstrcnvgust(i,j)**2)**2                                    &
                     + 4.0 * u_s2 * u_s2 )) )
          v_s_std(l) = SQRT( 0.5 * ( beta_cndd * beta_cndd * w_s * w_s        &
                 + wstrcnvgust(i,j)**2                                        &
                 + SQRT( ((beta_cndd * w_s)**2                                &
                 + wstrcnvgust(i,j)**2)**2                                    &
                 + 4.0 * u_s_std2 * u_s_std2 )) )
        END SELECT
      ELSE
        v_s(l) = SQRT( u_s2 )
        v_s_std(l) = SQRT( u_s_std2 )
      END IF

      IF ( v_s(l) > min_wind) THEN
        recip_l_mo(l) = -vkman * b_flux                                       &
                      / (v_s(l) * v_s(l) * v_s(l))
      ELSE
        IF (cor_mo_iter >= Limit_ObukhovL) THEN
          recip_l_mo(l) = 1.0 / SQRT(EPSILON(0.0))
        ELSE
          recip_l_mo(l) = SQRT(EPSILON(0.0))
        END IF
      END IF
    END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k,l,j,i)             &
!$OMP& SHARED(surft_pts,surft_index,pts_index,t_i_length,db,z1_tq,            &
!$OMP&  z1_uv,recip_l_mo,cor_mo_iter,vshr) IF(surft_pts>1)
    DO k = 1,surft_pts
      l = surft_index(k)
      j=(pts_index(l) - 1) / t_i_length + 1
      i = pts_index(l) - (j-1) * t_i_length
      !
      IF (cor_mo_iter >= Limit_ObukhovL) THEN
        IF ( db(l) > Ri_m  * z1_tq(i,j)                                       &
                       *(vshr(i,j) / z1_uv(i,j))**2.0 ) THEN
          recip_l_mo(l)= 1.5 * Ri_m**2.0 / z1_tq(i,j)
        END IF
      END IF
    END DO
!$OMP END PARALLEL DO

    IF (l_vegdrag) THEN
      CALL can_drag_z0(points, surft_pts, surft_index,                        &
                       recip_l_mo, canht, lai,                                &
                       z0m, z0h, zdt)

      SELECT CASE (i_modiscopt)
      CASE (off)
        CALL can_drag_phi_m_h(points, surft_pts, surft_index, pts_index,      &
                              recip_l_mo, z1_uv, z1_tq, canht, lai, z0m, z0h, &
                              phi_m, phi_h)
      CASE (on)
        CALL can_drag_phi_m_h_vol(points, surft_pts, surft_index, pts_index,  &
                       recip_l_mo, z1_uv_top, z1_tq_top, canht, lai, z0m, z0h,&
                       phi_m, phi_h)
      END SELECT

      !CDIR NODEP
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k,l)                 &
!$OMP& SHARED(surft_pts,surft_index,chv,cdv,cdv_std,phi_h,                    &
!$OMP& phi_m,v_s_std,v_s,wind_profile_factor) IF(surft_pts>1)
      DO k = 1,surft_pts
        l = surft_index(k)
        chv(l) = ( vkman / phi_h(l) ) * v_s_std(l)
        cdv_std(l) = ( vkman / phi_m(l) ) * v_s_std(l)
        cdv(l) = cdv_std(l) * ( v_s(l) / v_s_std(l) ) /                       &
                              wind_profile_factor(l)
      END DO
!$OMP END PARALLEL DO
    ELSE
      SELECT CASE (i_modiscopt)
      CASE (off)
        CALL phi_m_h ( points,surft_pts,surft_index,pts_index,                &
                       recip_l_mo,z1_uv,z1_tq,z0m,z0h,                        &
                       phi_m,phi_h)
      CASE (on)
        CALL phi_m_h_vol ( points,surft_pts,surft_index,pts_index,            &
                       recip_l_mo,z1_uv_top,z1_tq_top,z0m,z0h,                &
                       phi_m,phi_h)
      END SELECT


      !CDIR NODEP
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k,l)                 &
!$OMP& SHARED(surft_pts,surft_index,chv,cdv,cdv_std,phi_h,                    &
!$OMP& phi_m,v_s_std,v_s,wind_profile_factor) IF(surft_pts>1)
      DO k = 1,surft_pts
        l = surft_index(k)
        chv(l) = ( vkman / phi_h(l) ) * v_s_std(l)
        cdv(l) = ( vkman / phi_m(l) ) * v_s(l)
        cdv_std(l) = cdv(l) * ( v_s_std(l) / v_s(l) ) *                       &
                              wind_profile_factor(l)
      END DO
!$OMP END PARALLEL DO
    END IF
  END DO ! Iteration loop

END IF ! test on cor_mo_iter eq Improve_Initial_Guess
!-----------------------------------------------------------------------
! Set CD's and CH's to be dimensionless paremters
!-----------------------------------------------------------------------
!CDIR NODEP
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k,l,j,i)             &
!$OMP& SHARED(surft_pts,surft_index,pts_index,t_i_length,cdv,                 &
!$OMP& vshr,cdv_std,chv) IF(surft_pts>1)
DO k = 1,surft_pts
  l = surft_index(k)
  j=(pts_index(l) - 1) / t_i_length + 1
  i = pts_index(l) - (j-1) * t_i_length
  cdv(l) = cdv(l) / vshr(i,j)
  cdv_std(l) = cdv_std(l) / vshr(i,j)
  chv(l) = chv(l) / vshr(i,j)
END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fcdch
END MODULE fcdch_mod
