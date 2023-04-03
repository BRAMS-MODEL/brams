module ModPostOneField3d
   use micphys, only : &
         mcphys_type
   !
   use mem_grid   , only: time           ! INTENT(IN)
   use mem_varinit, only: vtime1, vtime2 ! INTENT(IN)
   use node_mod, only:  mchnum, mynum,master_num

   use ModOutputUtils, only : GetVarFromMemToOutput

   use ModPostGrid, only : OutputGradsField
   use ModBramsGrid, only : BramsGrid
   use ModPostTypes, only : PostGrid
   use ModPostOneFieldUtils, only : PrepareAndOutputGradsField
   !use ModPostOneFieldUtils, only : PostVarType
   use ModPostTypes, only: PostVarType

   use ModPostUtils, only : rams_comp_tempk
   use ModPostUtils, only : rams_comp_tempc
   use ModPostUtils, only : rams_comp_dewk
   use ModPostUtils, only : rams_comp_rh
   use ModPostUtils, only : rams_comp_press
   use ModPostUtils, only : rams_comp_dn0
   use ModPostUtils, only : rams_comp_z
   use ModPostUtils, only : rams_comp_rotate
   use ModPostUtils, only : rams_comp_avgu
   use ModPostUtils, only : rams_comp_avgv
   use ModPostUtils, only : calc_omeg

   use io_params, only : & ! 
       IPOS

   !LFR
   use modTimeLineFRN, only: writeTimeLineFRN

   use ModNamelistFile, only: namelistFile

   implicit none

   private

   !--(DMK-CCATT-INI)-------------------------------------------------------
   real, parameter :: PMAR = 28.96
   !--(DMK-CCATT-FIM)-------------------------------------------------------

   public :: Brams2Post_3d
   include "constants.f90"

contains


   subroutine Brams2Post_3d (one_post_variable, oneBramsGrid, onePostGrid, oneNamelistFile)
      type(NamelistFile), pointer :: oneNamelistFile
      type(PostVarType) :: one_post_variable
      type(BramsGrid), pointer :: oneBramsGrid
      type(PostGrid), pointer :: onePostGrid

      real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
      real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
      real :: ScrT3N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
      real :: ScrT3N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
      real :: ScrT3N04(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
      real :: ScrT3N05(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

      real :: ScrT2N01(oneBramsGrid%mxp, oneBramsGrid%myp)
      integer :: firstX, lastX, firstY, lastY, k
      real :: tfact
      tfact=(time-vtime1)/(vtime2-vtime1)
      
      select case (one_post_variable%fieldName)

      case ('U')
         call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, OutputField)
      case ('V')
         call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, OutputField)
      case ('REF_U')
         call GetVarFromMemToOutput ('VARUP', oneBramsGrid%currGrid, ScrT3N01)
         call GetVarFromMemToOutput ('VARUF', oneBramsGrid%currGrid, ScrT3N02)
         OutputField = ScrT3N01 +(ScrT3N02 -ScrT3N01 )*tfact

      case ('REF_V')
         call GetVarFromMemToOutput ('VARVP', oneBramsGrid%currGrid, ScrT3N01)
         call GetVarFromMemToOutput ('VARVF', oneBramsGrid%currGrid, ScrT3N02)
         OutputField = ScrT3N01 +(ScrT3N02 -ScrT3N01 )*tfact
	 !print*,"time=",tfact,time/3600.,vtime1/3600., vtime2/3600.!, mynum,master_num 
	 !call flush(6)

      case ('TEMPK')
         call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, OutputField)
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N01)
         call rams_comp_tempk (OutputField, ScrT3N01)
      case ('TKE')
         call GetVarFromMemToOutput ('TKEP', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('DEWPTC')
         call GetVarFromMemToOutput ('RV', oneBramsGrid%currGrid, OutputField)
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N01)
         call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N02)
         call rams_comp_dewk (OutputField, ScrT3N01, ScrT3N02)
         call rams_comp_tempc (OutputField)
      case ('CLOUD')
         call GetVarFromMemToOutput ('RCP', oneBramsGrid%currGrid, OutputField)
         OutputField = OutputField * 1.e3
         OutputField = max(OutputField, 0.0)
      case ('RV')
         call GetVarFromMemToOutput ('RV', oneBramsGrid%currGrid, OutputField)
         OutputField = OutputField * 1.e3
         OutputField = max(OutputField, 0.0)
      case ('REF_RV')
         call GetVarFromMemToOutput ('VARRP', oneBramsGrid%currGrid, ScrT3N01)
         call GetVarFromMemToOutput ('VARRF', oneBramsGrid%currGrid, ScrT3N02)
         OutputField = ScrT3N01 +(ScrT3N02-ScrT3N01 )*tfact
         OutputField = OutputField * 1.e3
         OutputField = max(OutputField, 0.0)

      case ('GEO')
         call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
         call rams_comp_z (OutputField, ScrT2N01, oneBramsGrid%ztn, oneBramsGrid%ztop)
      case ('UE_AVG')
         firstX = oneBramsGrid%nodei0(oneBramsGrid%mynum) + 1
         lastX = oneBramsGrid%nodei0(oneBramsGrid%mynum) + oneBramsGrid%mxp
         firstY = oneBramsGrid%nodej0(oneBramsGrid%mynum) + 1
         lastY = oneBramsGrid%nodej0(oneBramsGrid%mynum) + oneBramsGrid%myp
         call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, OutputField)
         call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, ScrT3N01)
         call rams_comp_rotate (OutputField, ScrT3N01, &
               oneBramsGrid%xtn(firstX : lastX), oneBramsGrid%ytn(firstY : lastY), &
               oneBramsGrid%polelat, oneBramsGrid%polelon)
         call rams_comp_avgu (OutputField)
      case ('VE_AVG')
         firstX = oneBramsGrid%nodei0(oneBramsGrid%mynum) + 1
         lastX = oneBramsGrid%nodei0(oneBramsGrid%mynum) + oneBramsGrid%mxp
         firstY = oneBramsGrid%nodej0(oneBramsGrid%mynum) + 1
         lastY = oneBramsGrid%nodej0(oneBramsGrid%mynum) + oneBramsGrid%myp
         call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, OutputField)
         call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, ScrT3N01)
         call rams_comp_rotate (ScrT3N01, OutputField, &
               oneBramsGrid%xtn(firstX : lastX), oneBramsGrid%ytn(firstY : lastY), &
               oneBramsGrid%polelat, oneBramsGrid%polelon)
         call rams_comp_avgv (OutputField)
!Introduzir o cálculo de Magnitude e Direção do vento
      case ('MAGUV')
         !Pega UE_AVG
         firstX = oneBramsGrid%nodei0(oneBramsGrid%mynum) + 1
         lastX = oneBramsGrid%nodei0(oneBramsGrid%mynum) + oneBramsGrid%mxp
         firstY = oneBramsGrid%nodej0(oneBramsGrid%mynum) + 1
         lastY = oneBramsGrid%nodej0(oneBramsGrid%mynum) + oneBramsGrid%myp
         call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, ScrT3N02)
         call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, ScrT3N01)
         !UE_AVG
         call rams_comp_rotate (ScrT3N02, ScrT3N01, &
               oneBramsGrid%xtn(firstX : lastX), oneBramsGrid%ytn(firstY : lastY), &
               oneBramsGrid%polelat, oneBramsGrid%polelon)
         call rams_comp_avgu (ScrT3N02) 
         !VE_AVG
         call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, ScrT3N03)
         call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, ScrT3N01)
         call rams_comp_rotate (ScrT3N01, ScrT3N03, &
               oneBramsGrid%xtn(firstX : lastX), oneBramsGrid%ytn(firstY : lastY), &
               oneBramsGrid%polelat, oneBramsGrid%polelon)
         call rams_comp_avgv (ScrT3N03)    
         OutputField = sqrt(ScrT3N02*ScrT3N02+ScrT3N03*ScrT3N03)

      case('DIRUV')
         !Pega UE_AVG
         firstX = oneBramsGrid%nodei0(oneBramsGrid%mynum) + 1
         lastX = oneBramsGrid%nodei0(oneBramsGrid%mynum) + oneBramsGrid%mxp
         firstY = oneBramsGrid%nodej0(oneBramsGrid%mynum) + 1
         lastY = oneBramsGrid%nodej0(oneBramsGrid%mynum) + oneBramsGrid%myp
         call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, ScrT3N02)
         call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, ScrT3N01)
         !UE_AVG
         call rams_comp_rotate (ScrT3N02, ScrT3N01, &
               oneBramsGrid%xtn(firstX : lastX), oneBramsGrid%ytn(firstY : lastY), &
               oneBramsGrid%polelat, oneBramsGrid%polelon)
         call rams_comp_avgu (ScrT3N02) 
         !VE_AVG
         call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, ScrT3N03)
         call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, ScrT3N01)
         call rams_comp_rotate (ScrT3N01, ScrT3N03, &
               oneBramsGrid%xtn(firstX : lastX), oneBramsGrid%ytn(firstY : lastY), &
               oneBramsGrid%polelat, oneBramsGrid%polelon)
         call rams_comp_avgv (ScrT3N03)    
         ScrT3N04 = atan(ScrT3N03/(ScrT3N02 + 1.E-20))
         OutputField = c_i_pi180*ScrT3N04


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      case ('TEMPC')
         call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, OutputField)
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N01)
         call rams_comp_tempk (OutputField, ScrT3N01)
         call rams_comp_tempc (OutputField)
      case ('REF_TEMPC')
         call GetVarFromMemToOutput ('VARTP', oneBramsGrid%currGrid, ScrT3N01)
         call GetVarFromMemToOutput ('VARPP', oneBramsGrid%currGrid, ScrT3N02)
         call GetVarFromMemToOutput ('VARTF', oneBramsGrid%currGrid, ScrT3N03)
         call GetVarFromMemToOutput ('VARPF', oneBramsGrid%currGrid, ScrT3N04)
         OutputField=ScrT3N01 +(ScrT3N03 -ScrT3N01 )*tfact ! theta
	 ScrT3N02=ScrT3N02 +(ScrT3N04 -ScrT3N02 )*tfact !pert exner function
	 call GetVarFromMemToOutput ('PI0', oneBramsGrid%currGrid, ScrT3N03)
	 ScrT3N01=ScrT3N02+ScrT3N03 !exner function
         call rams_comp_tempk (OutputField, ScrT3N01)
         call rams_comp_tempc (OutputField)
      case ('RH')
         call GetVarFromMemToOutput ('RV', oneBramsGrid%currGrid, OutputField)
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N01)
         call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N02)
         call rams_comp_rh (OutputField, ScrT3N01, ScrT3N02)
         OutputField = max(OutputField, 0.0)
      case ('W')
         call GetVarFromMemToOutput ('WP', oneBramsGrid%currGrid, OutputField)
      case ('OMEG')
         call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
         call rams_comp_dn0 (ScrT3N01, ScrT3N02, ScrT3N03, ScrT2N01, oneBramsGrid%pi01dn, &
               oneBramsGrid%th01dn, oneBramsGrid%ztn, oneBramsGrid%ztop, oneBramsGrid%dzmn)
         call GetVarFromMemToOutput ('WP', oneBramsGrid%currGrid, ScrT3N01)
         call calc_omeg (OutputField, ScrT3N01, ScrT3N03)
      case ('DFXUP1')
         call GetVarFromMemToOutput ('DFXUP1', oneBramsGrid%currGrid, OutputField)
      case ('EFXUP1')
         call GetVarFromMemToOutput ('EFXUP1', oneBramsGrid%currGrid, OutputField)
      case ('DFXDN1')
         call GetVarFromMemToOutput ('DFXDN1', oneBramsGrid%currGrid, OutputField)
      case ('EFXDN1')
         call GetVarFromMemToOutput ('EFXDN1', oneBramsGrid%currGrid, OutputField)
      case ('CFXUP1')
         call GetVarFromMemToOutput ('CFXUP1', oneBramsGrid%currGrid, OutputField)
      case ('CFXDN1')
         call GetVarFromMemToOutput ('CFXDN1', oneBramsGrid%currGrid, OutputField)
      case ('KHV')
         call GetVarFromMemToOutput ('VKH', oneBramsGrid%currGrid, OutputField)
      case ('CO2')
         call GetVarFromMemToOutput ('CO2P', oneBramsGrid%currGrid, OutputField)
         OutputField = OutputField * (28.96 / 44.) * 1.E-3  ! TRANSFORMACAO DE kg[CO2]/kg[AR] para ppm (PARTE POR MILHAO)
      case ('CO2_BURN')
         call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
         call GetVarFromMemToOutput ('CO2_bburn_SRC', oneBramsGrid%currGrid, OutputField)
         do k = 1, oneBramsGrid%mzp
            OutputField(:, :, k) = OutputField(:, :, k) * (1. - ScrT2N01(:, :) / oneBramsGrid%ztop) &
                  / oneBramsGrid%dztn(2) / 86400. * 1.e-3  ! convertendo de kg/m3/dia para mg/m2/s
         enddo
      case ('THETA')
         call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, OutputField)
      case ('REF_THETA')
         call GetVarFromMemToOutput ('VARTP', oneBramsGrid%currGrid, ScrT3N01)
         call GetVarFromMemToOutput ('VARTF', oneBramsGrid%currGrid, ScrT3N02)
         OutputField = ScrT3N01 +(ScrT3N02-ScrT3N01 )*tfact
      case ('AGGREGATES')
         !-For GThompson microphysics micphys_type>1 : RAP does not exist
         if(mcphys_type .le. 1) then
            call GetVarFromMemToOutput ('RAP', oneBramsGrid%currGrid, OutputField)
         else
            OutputField = 0.0
         endif
         OutputField = OutputField * 1.e3
         OutputField = max(OutputField, 0.0)
      case ('GRAUPEL')
         call GetVarFromMemToOutput ('RGP', oneBramsGrid%currGrid, OutputField)
         OutputField = OutputField * 1.e3
         OutputField = max(OutputField, 0.0)
      case ('PRISTINE')
         call GetVarFromMemToOutput ('RGP', oneBramsGrid%currGrid, OutputField)
         OutputField = OutputField * 1.e3
         OutputField = max(OutputField, 0.0)
      case ('HAIL')
         if(mcphys_type .le. 1) then
            call GetVarFromMemToOutput ('RHP', oneBramsGrid%currGrid, OutputField)
         else
            OutputField = 0.0
         endif
         OutputField = OutputField * 1.e3
         OutputField = max(OutputField, 0.0)
      case ('SNOW')
         call GetVarFromMemToOutput ('RSP', oneBramsGrid%currGrid, OutputField)
         OutputField = OutputField * 1.e3
         OutputField = max(OutputField, 0.0)
      case ('LIQUID')
         OutputField = 0.0
         call GetVarFromMemToOutput ('RCP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + ScrT3N01
         call GetVarFromMemToOutput ('RRP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + ScrT3N01

         if(mcphys_type .le. 1) then
            call GetVarFromMemToOutput ('RGP', oneBramsGrid%currGrid, ScrT3N01)
            call GetVarFromMemToOutput ('Q6', oneBramsGrid%currGrid, ScrT3N02)
            ScrT3N02 = min(ScrT3N02, 334000.)
            ScrT3N02 = max(ScrT3N02, 0.0)
            ScrT3N02 = ScrT3N02 / 334000.
            ScrT3N01 = ScrT3N01 * ScrT3N02
            OutputField = OutputField + ScrT3N01

            call GetVarFromMemToOutput ('RHP', oneBramsGrid%currGrid, ScrT3N01)
            call GetVarFromMemToOutput ('Q7', oneBramsGrid%currGrid, ScrT3N02)
            ScrT3N02 = min(ScrT3N02, 334000.)
            ScrT3N02 = max(ScrT3N02, 0.0)
            ScrT3N02 = ScrT3N02 / 334000.
            ScrT3N01 = ScrT3N01 * ScrT3N02
            OutputField = OutputField + ScrT3N01
         endif

         OutputField = OutputField * 1.e3
         OutputField = max(OutputField, 0.0)
      case ('RAIN')
         call GetVarFromMemToOutput ('RRP', oneBramsGrid%currGrid, OutputField)
         OutputField = OutputField * 1.e3
         OutputField = max(OutputField, 0.0)
      case ('ICE')
         OutputField = 0.0
         call GetVarFromMemToOutput ('RPP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + ScrT3N01
         call GetVarFromMemToOutput ('RSP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + ScrT3N01
         !-For GThompson microphysics micphys_type>1 : RAP does not exist
         if(mcphys_type .le. 1) then
            call GetVarFromMemToOutput ('RAP', oneBramsGrid%currGrid, ScrT3N01)
            OutputField = OutputField + ScrT3N01
         endif
         call GetVarFromMemToOutput ('RGP', oneBramsGrid%currGrid, ScrT3N01)

         if(mcphys_type .le. 1) then
            call GetVarFromMemToOutput ('Q6', oneBramsGrid%currGrid, ScrT3N02)

            ScrT3N02 = max(ScrT3N02, 0.0)
            ScrT3N02 = min(ScrT3N02, 334000.)
            ScrT3N02 = 1. - (ScrT3N02 / 334000.)
            ScrT3N01 = ScrT3N01 * ScrT3N02
         endif
         OutputField = OutputField + ScrT3N01

         if(mcphys_type .le. 1) then
            call GetVarFromMemToOutput ('RHP', oneBramsGrid%currGrid, ScrT3N01)
            call GetVarFromMemToOutput ('Q7', oneBramsGrid%currGrid, ScrT3N02)

            ScrT3N02 = max(ScrT3N02, 0.0)
            ScrT3N02 = min(ScrT3N02, 334000.)
            ScrT3N02 = 1. - (ScrT3N02 / 334000.)
            ScrT3N01 = ScrT3N01 * ScrT3N02
            OutputField = OutputField + ScrT3N01
         endif

         OutputField = OutputField * 1.e3
         OutputField = max(OutputField, 0.0)
      case ('DEWPTK')
         call GetVarFromMemToOutput ('RV', oneBramsGrid%currGrid, OutputField)
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N01)
         call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N02)
         call rams_comp_dewk(OutputField, ScrT3N01, ScrT3N02)
         call rams_comp_tempk(OutputField, ScrT3N01)
      case ('PRESS')
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, OutputField)
         call rams_comp_press(OutputField)
      case ('TOTAL_COND')
         OutputField = 0.0
         call GetVarFromMemToOutput ('RCP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + ScrT3N01
         call GetVarFromMemToOutput ('RRP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + ScrT3N01
         call GetVarFromMemToOutput ('RPP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + ScrT3N01
         call GetVarFromMemToOutput ('RSP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + ScrT3N01
         !-For GThompson microphysics micphys_type>1 : RAP does not exist
         if(mcphys_type .le. 1) then
            call GetVarFromMemToOutput ('RAP', oneBramsGrid%currGrid, ScrT3N01)
            OutputField = OutputField + ScrT3N01
         endif
         call GetVarFromMemToOutput ('RGP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + ScrT3N01
         if(mcphys_type .le. 1) then
            call GetVarFromMemToOutput ('RHP', oneBramsGrid%currGrid, ScrT3N01)
            OutputField = OutputField + ScrT3N01
         endif
         OutputField = OutputField * 1.e3
         OutputField = max(OutputField, 0.0)
         !--(DMK-CCATT-INI)-------------------------------------------------------
      case ('CO')
         ! ierr= RAMS_getvar('COP',idim_type,ngrd,a,b,flnm)
         call GetVarFromMemToOutput ('COP', oneBramsGrid%currGrid, OutputField)
         ! call RAMS_comp_noneg(n1,n2,n3,a)
         OutputField = max(OutputField, 0.0)
         ! call RAMS_transf_ppb(n1,n2,n3,a,28.)
         OutputField = OutputField * (PMAR / 28.)
      case ('NO')
         ! ierr= RAMS_getvar('NO2P',idim_type,ngrd,a,b,flnm)
         call GetVarFromMemToOutput ('NOP', oneBramsGrid%currGrid, OutputField)
         ! call RAMS_comp_noneg(n1,n2,n3,a)
         OutputField = max(OutputField, 0.0)
         ! call RAMS_transf_ppb(n1,n2,n3,a,30.)
         OutputField = OutputField * (PMAR / 30.)
      case ('HNO3')
         ! ierr= RAMS_getvar('HNO3P',idim_type,ngrd,a,b,flnm)
         call GetVarFromMemToOutput ('HNO3P', oneBramsGrid%currGrid, OutputField)
         ! call RAMS_comp_noneg(n1,n2,n3,a)
         OutputField = max(OutputField, 0.0)
         ! call RAMS_transf_ppb(n1,n2,n3,a,63.)
         OutputField = OutputField * (PMAR / 63.)
      case ('O3')
         ! ierr= RAMS_getvar('O3P',idim_type,ngrd,a,b,flnm)
         call GetVarFromMemToOutput ('O3P', oneBramsGrid%currGrid, OutputField)
         ! call RAMS_comp_noneg(n1,n2,n3,a)
         OutputField = max(OutputField, 0.0)
         ! call RAMS_transf_ppb(n1,n2,n3,a,48.)
         OutputField = OutputField * (PMAR / 48.)
      case ('NO2')
         ! ierr= RAMS_getvar('NO2P',idim_type,ngrd,a,b,flnm)
         call GetVarFromMemToOutput ('NO2P', oneBramsGrid%currGrid, OutputField)
         ! call RAMS_comp_noneg(n1,n2,n3,a)
         OutputField = max(OutputField, 0.0)
         ! call RAMS_transf_ppb(n1,n2,n3,a,46.)
         OutputField = OutputField * (PMAR / 46.)
      case ('PM25')
         ! ierr= RAMS_getvar('bburn2P',idim_type,ngrd,a,b,flnm)
         call GetVarFromMemToOutput ('bburn2P', oneBramsGrid%currGrid, OutputField)
         call GetVarFromMemToOutput ('urban2P', oneBramsGrid%currGrid, ScrT3N04)
         call GetVarFromMemToOutput ('urban3P', oneBramsGrid%currGrid, ScrT3N05)
         ! call RAMS_comp_noneg(n1,n2,n3,a)
         OutputField = max(OutputField + ScrT3N04 + ScrT3N05, 0.0)
         ! call RAMS_comp_mults(n1,n2,n3,a,1.e-9)  ! converte de 1e-9 kg/kg para 1 kg/kg
         OutputField = OutputField * 1.e-9
         ! ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
         call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
         ! call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)
         call rams_comp_dn0 (ScrT3N01, ScrT3N02, ScrT3N03, ScrT2N01, oneBramsGrid%pi01dn, &
               oneBramsGrid%th01dn, oneBramsGrid%ztn, oneBramsGrid%ztop, oneBramsGrid%dzmn)
         ! call RAMS_transf_ugm3(n1,n2,n3,a,d)
         OutputField = OutputField * scrT3N03 * 1.e+9
         ! call RAMS_comp_noneg(n1,n2,n3,a)
         !OutputField = max(OutputField, 0.0)
      case ('CO_SRC')
         call GetVarFromMemToOutput ('CO_bburn_SRC', oneBramsGrid%currGrid, OutputField)
         call GetVarFromMemToOutput ('CO_antro_SRC', oneBramsGrid%currGrid, ScrT3N01)
         call GetVarFromMemToOutput ('CO_bioge_SRC', oneBramsGrid%currGrid, ScrT3N02)
         OutputField(:, :, 2) = max(OutputField(:, :, 2) + ScrT3N01(:, :, 1) + ScrT3N02(:, :, 1), 0.0)
         OutputField = OutputField * 1.e-9 ! converte de mg/kg para kg/kg
      case ('NO_SRC')
         call GetVarFromMemToOutput ('NO_bburn_SRC', oneBramsGrid%currGrid, OutputField)
         call GetVarFromMemToOutput ('NO_antro_SRC', oneBramsGrid%currGrid, ScrT3N01)
         call GetVarFromMemToOutput ('NO_bioge_SRC', oneBramsGrid%currGrid, ScrT3N02)
         OutputField(:, :, 2) = max(OutputField(:, :, 2) + ScrT3N01(:, :, 1) + ScrT3N02(:, :, 1), 0.0)
         OutputField = OutputField * 1.e-9 ! converte de mg/kg para kg/kg
      case('NMVOCM')
         OutputField = 0.
         call GetVarFromMemToOutput ('ETHP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + max(ScrT3N01, 0.0)
         call GetVarFromMemToOutput ('ALKAP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + max(ScrT3N01, 0.0)
         call GetVarFromMemToOutput ('ALKEP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + max(ScrT3N01, 0.0)
         call GetVarFromMemToOutput ('BIOP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + max(ScrT3N01, 0.0)
         call GetVarFromMemToOutput ('AROP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + max(ScrT3N01, 0.0)
         call GetVarFromMemToOutput ('HCHOP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + max(ScrT3N01, 0.0)
         call GetVarFromMemToOutput ('ALDP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + max(ScrT3N01, 0.0)
         call GetVarFromMemToOutput ('KETP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + max(ScrT3N01, 0.0)
         call GetVarFromMemToOutput ('CRBOP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + max(ScrT3N01, 0.0)
         call GetVarFromMemToOutput ('ONITP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + max(ScrT3N01, 0.0)
         call GetVarFromMemToOutput ('PANP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + max(ScrT3N01, 0.0)
         call GetVarFromMemToOutput ('OP1P', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + max(ScrT3N01, 0.0)
         call GetVarFromMemToOutput ('OP2P', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + max(ScrT3N01, 0.0)
         call GetVarFromMemToOutput ('ORA1P', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + max(ScrT3N01, 0.0)
         call GetVarFromMemToOutput ('ORA2P', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + max(ScrT3N01, 0.0)
         ! call GetVarFromMemToOutput ('MO2P', oneBramsGrid%currGrid, ScrT3N01)
         ! OutputField = OutputField + max(ScrT3N01, 0.0)
         ! call GetVarFromMemToOutput ('AKAPP', oneBramsGrid%currGrid, ScrT3N01)
         ! OutputField = OutputField + max(ScrT3N01, 0.0)
         !  call GetVarFromMemToOutput ('AKEPP', oneBramsGrid%currGrid, ScrT3N01)
         !  OutputField = OutputField + max(ScrT3N01, 0.0)
         call GetVarFromMemToOutput ('BIOPP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + max(ScrT3N01, 0.0)
         !  call GetVarFromMemToOutput ('PHOP', oneBramsGrid%currGrid, ScrT3N01)
         !  OutputField = OutputField + max(ScrT3N01, 0.0)
         !  call GetVarFromMemToOutput ('ADDP', oneBramsGrid%currGrid, ScrT3N01)
         !  OutputField = OutputField + max(ScrT3N01, 0.0)
         call GetVarFromMemToOutput ('AROPP', oneBramsGrid%currGrid, ScrT3N01)
         OutputField = OutputField + max(ScrT3N01, 0.0)
         !call GetVarFromMemToOutput ('CBOPP', oneBramsGrid%currGrid, ScrT3N01)
         !OutputField = OutputField + max(ScrT3N01, 0.0)
         ! call GetVarFromMemToOutput ('OLNP', oneBramsGrid%currGrid, ScrT3N01)
         ! OutputField = OutputField + max(ScrT3N01, 0.0)
         !--(DMK-CCATT-FIM)-------------------------------------------------------
      case ('CUTHDP')
         call GetVarFromMemToOutput ('THSRC', oneBramsGrid%currGrid, OutputField)
         OutputField = OutputField * 86400.0
      case ('CURTDP')
         call GetVarFromMemToOutput ('RTSRC', oneBramsGrid%currGrid, OutputField)
         OutputField = OutputField * 86400.0 * 1000.0
      case ('CUUDP')
         call GetVarFromMemToOutput ('USRC', oneBramsGrid%currGrid, OutputField)
         OutputField = OutputField * 86400.0
      case ('CUVDP')
         call GetVarFromMemToOutput ('VSRC', oneBramsGrid%currGrid, OutputField)
         OutputField = OutputField * 86400.0
      case ('CUCLDP')
         call GetVarFromMemToOutput ('CLSRC', oneBramsGrid%currGrid, OutputField)

      case ('ZMFUP')
         call GetVarFromMemToOutput ('ZMFUP', oneBramsGrid%currGrid, OutputField)
      case ('ZMFDD')
         call GetVarFromMemToOutput ('ZMFDD', oneBramsGrid%currGrid, OutputField)
      case ('ZMFSH')
         call GetVarFromMemToOutput ('ZMFSH', oneBramsGrid%currGrid, OutputField)
      case ('ZUPMD')
         call GetVarFromMemToOutput ('ZUPMD', oneBramsGrid%currGrid, OutputField)
      case ('CUBUDP')
         call GetVarFromMemToOutput ('BUOYSRC', oneBramsGrid%currGrid, OutputField)
         OutputField = OutputField * 86400.0
      case ('BUX')
         call GetVarFromMemToOutput ('SCLP001', oneBramsGrid%currGrid, OutputField)

      case ('LSFTH')
         call GetVarFromMemToOutput ('LSFTH', oneBramsGrid%currGrid, OutputField)
         OutputField = OutputField * 86400.0
      case ('LSFRT')
         call GetVarFromMemToOutput ('LSFRT', oneBramsGrid%currGrid, OutputField)
         OutputField = OutputField * 86400.0* 1000.0
      case ('LSFTH_SH')
         call GetVarFromMemToOutput ('LSFTH_SH', oneBramsGrid%currGrid, OutputField)
         OutputField = OutputField * 86400.0
      case ('LSFRT_SH')
         call GetVarFromMemToOutput ('LSFRT_SH', oneBramsGrid%currGrid, OutputField)
         OutputField = OutputField * 86400.0* 1000.0

!Variáveis do Recycle
      case ('R_O3'     )
         call GetVarFromMemToOutput ('O3P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_H2O2'   )
         call GetVarFromMemToOutput ('H2O2P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_NO'     )
         call GetVarFromMemToOutput ('NOP', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_NO2'    )
         call GetVarFromMemToOutput ('NO2P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_NO3'    )
         call GetVarFromMemToOutput ('NO3P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_N2O5'   )
         call GetVarFromMemToOutput ('N2O5P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_HONO'   )
         call GetVarFromMemToOutput ('HONOP', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_HNO3'   )
         call GetVarFromMemToOutput ('HNO3P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_HNO4'   )
         call GetVarFromMemToOutput ('HNO4P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_SO2 '   )
         call GetVarFromMemToOutput ('SO2P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_SULF'   )
         call GetVarFromMemToOutput ('SULFP', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_CO'     )
         call GetVarFromMemToOutput ('COP', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_HO'     )
         call GetVarFromMemToOutput ('HOP', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_HO2'    )
         call GetVarFromMemToOutput ('HO2P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_CH4'    )
         call GetVarFromMemToOutput ('CH4P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_ETH'    )
         call GetVarFromMemToOutput ('ETHP', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_ALKA'   )
         call GetVarFromMemToOutput ('ALKAP', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_ALKE'   )
         call GetVarFromMemToOutput ('ALKEP', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_BIO'    )
         call GetVarFromMemToOutput ('BIOP', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_ARO'    )
         call GetVarFromMemToOutput ('AROP', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_HCHO'   )
         call GetVarFromMemToOutput ('HCHOP', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_ALD'    )
         call GetVarFromMemToOutput ('ALDP', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_KET'    )
         call GetVarFromMemToOutput ('KETP', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_CRBO'   )
         call GetVarFromMemToOutput ('CRBOP', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_ONIT'   )
         call GetVarFromMemToOutput ('ONITP', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_PAN'    )
         call GetVarFromMemToOutput ('PANP', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_OP1'    )
         call GetVarFromMemToOutput ('OP1P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_OP2'    )
         call GetVarFromMemToOutput ('OP2P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_ORA1'   )
         call GetVarFromMemToOutput ('ORA1P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_ORA2'   )
         call GetVarFromMemToOutput ('ORA2P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_BBURN2' )
         call GetVarFromMemToOutput ('bburn2P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_BBURN3' )
         call GetVarFromMemToOutput ('bburn3P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_URBAN2' )
         call GetVarFromMemToOutput ('urban2P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_URBAN3' )
         call GetVarFromMemToOutput ('urban3P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_MARIN1' )
         call GetVarFromMemToOutput ('marin1P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_MARIN2' )
         call GetVarFromMemToOutput ('marin2P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_MARIN3' )
         call GetVarFromMemToOutput ('marin3P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_V_ASH1' )
         call GetVarFromMemToOutput ('v_ash1P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_V_ASH2' )
         call GetVarFromMemToOutput ('v_ash2P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_V_ASH3' )
         call GetVarFromMemToOutput ('v_ash3P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_V_ASH4' )
         call GetVarFromMemToOutput ('v_ash4P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_V_ASH5' )
         call GetVarFromMemToOutput ('v_ash5P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_V_ASH6' )
         call GetVarFromMemToOutput ('v_ash6P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_V_ASH7' )
         call GetVarFromMemToOutput ('v_ash7P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_V_ASH8' )
         call GetVarFromMemToOutput ('v_ash8P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_V_ASH9' )
         call GetVarFromMemToOutput ('v_ash9P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('R_V_ASH10')
         call GetVarFromMemToOutput ('v_ash10P', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
!End vaviaveis do recycle


      case default
         write(*, "(a)") "**(OnePostField)** Post field 3d " // one_post_variable%fieldName // " not implemented!"
         stop
      end select

      if(IPOS == 11 .or. IPOS == 10) call writeTimeLineFRN(one_post_variable%fieldName,OutputField,oneBramsGrid%mxp &
      , oneBramsGrid%myp, oneBramsGrid%mzp,time,oneNamelistFile)

      if(IPOS == 11) return
      ! most of cases run this. Some just returns before
      call PrepareAndOutputGradsField(one_post_variable, oneBramsGrid, onePostGrid, OutputField)

   end subroutine Brams2Post_3d


end module ModPostOneField3d
