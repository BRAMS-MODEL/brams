module mem_tconv
  use chem1_list, only: nspecies_chem=> nspecies

  use aer1_list, only:  nspecies_aer=> nspecies,nmodes,ninorg

  implicit none
  integer, parameter :: ntotal=nspecies_chem+nspecies_aer*nmodes+nmodes+ninorg
  integer :: nspecies
  integer :: nchem_a,nchem_z
  integer, dimension(ntotal) :: ind_chem
  integer :: naer_a,naer_z
  integer :: naer_a_number,naer_z_number
  integer :: naer_a_inorg,naer_z_inorg
 
  integer, dimension(ntotal)   :: ind_aer
  integer, dimension(ntotal)   :: ind_mode
  integer, dimension(ntotal)   :: ind_mode_number
  integer, dimension(ntotal)   :: ind_aer_inorg

  logical           :: trans_conv_alloc = .false.
  real, allocatable ::&
       se     (:,:),    & ! environment scalar profile Z levels
       se_cup (:,:),    & ! environment scalar profile Z_cup levels
       sc_up  (:,:),    & ! updraft   gas-phase  scalar profile
       sc_dn  (:,:),    & ! downdraft gas-phase  scalar profile
       stcum1d(:,:),    & ! 1d convective tendency
       sc_up_c(:,:),    & ! updraft   aqueous-phase scalar profile
       sc_dn_c(:,:),    & ! DOwndraft aqueous-phase scalar profile
       pw_up  (:,:),    & ! updraft precitable gas/aer
       pw_dn  (:,:),	& ! downdraft precitable gas/aer
       henry_coef(:,:), & ! Henry's constant
       dn01d  (:)       ! 1d air density

contains

  subroutine alloc_trans_conv(NSPECIES_TRANSPORTED,mgmzp)

    use chem1_list, only: nspecies_chem=> nspecies &
                         ,spc_alloc_chem=> spc_alloc &
		         ,transport,on,off,spc_name
    use aer1_list, only:  nspecies_aer=> nspecies &
                         ,spc_alloc_aer=> spc_alloc &
			 ,mode_alloc  &
                         ,numb_alloc,inorg_mod_alloc&
			 ,inorg_alloc,numb_mod_alloc
			 
    USE mem_aer1, ONLY:  &
                        aerosol

    implicit none
    integer, intent (IN)::NSPECIES_TRANSPORTED,mgmzp
    integer :: ispc,nspecies,imode,ininorg


    if(trans_conv_alloc) then
       print *,'ERROR: trans_conv already allocated'
       print *,'Routine: trans_conv File: trans_conv.f90'
       stop
    end if

    nspecies = 0

    ! chemistry section : mapping

    nchem_a  = 1
    nchem_z  = 0
    do ispc = 1,nspecies_chem

     if(spc_alloc_chem(transport,ispc) == off) cycle
     nchem_z = nchem_z + 1
     ind_chem(nchem_z) = ispc
     !if(spc_name(ispc)=='CO') print*, 'CO INDEX= ',nchem_z, ind_chem(nchem_z)
    enddo
    nspecies = nchem_z ! number of chemical species with conv transport

    ! aerosol section : mapping
      naer_a = nchem_z + 1
      naer_z = nchem_z
!    naer_a = 1
!    naer_z = 0
    if(AEROSOL >= 1) then


      do ispc = 1,nspecies_aer
         do imode=1,nmodes

           if(mode_alloc   (imode,ispc          ) == on .and. &
              spc_alloc_aer(transport,imode,ispc) == on) then

	      naer_z = naer_z + 1
	      ind_aer (naer_z) = ispc
	      ind_mode(naer_z) = imode

            endif
         enddo
      enddo 
      !- total number of species (chem + aer mass) to be transported
      nspecies  =  naer_z
      
      naer_a_number = naer_z+1
      naer_z_number = naer_z
      naer_a_inorg  = naer_z+1
      naer_z_inorg  = naer_z
      if(AEROSOL == 2) then
      
      !- transport of aer number 
         ! naer_a_number = naer_z + 1
         ! naer_z_number = naer_z 
    
         do imode=1,nmodes
	    !print*,"aer2",imode,nmodes,numb_mod_alloc(imode),numb_alloc(transport,imode)
            if(numb_mod_alloc(imode)       == on .and. &
               numb_alloc(transport,imode) == on) then
	      
	      naer_z_number = naer_z_number + 1
	      ind_mode_number(naer_z_number) = imode
	      !print*,"index=",imode
            endif
         enddo
      !- transport of mass of inorg 
         naer_a_inorg  = naer_Z_number+1
         naer_z_inorg  = naer_Z_number
         do ininorg=1,ninorg
            
	    if(inorg_mod_alloc(ininorg)== on .and. &
               inorg_alloc(transport,ininorg) == on) then
	       
	       naer_z_inorg = naer_z_inorg + 1
	       ind_aer_inorg (naer_z_inorg) = ininorg

            endif
         enddo
      endif
      !- total number of species (chem + aer mass + aer number + aer mass inorg) to be transported
      nspecies  =  naer_z_inorg
    endif

    !print*,'xx=',nspecies,NSPECIES_TRANSPORTED,naer_z_inorg
    if(nspecies /=NSPECIES_TRANSPORTED) stop 'nspecies must not be different of nspecies_transported'
    if(nspecies == 0) return

    allocate( se    (nspecies,mgmzp), &
              se_cup(nspecies,mgmzp), &
	      sc_up (nspecies,mgmzp), &
	      sc_dn (nspecies,mgmzp), &
             stcum1d(nspecies,mgmzp), &
             sc_up_c(nspecies,mgmzp), &
	     sc_dn_c(nspecies,mgmzp), &
	     pw_up  (nspecies,mgmzp), &
	     pw_dn  (nspecies,mgmzp), &
             henry_coef(nspecies,mgmzp))


    allocate(dn01d(mgmzp) )

    trans_conv_alloc=.true.

  end subroutine alloc_trans_conv

  subroutine zero_tconv()

    implicit none

!!$    print *,'LFR->Zero tconv !'
    se=.0
    se_cup=.0
    sc_up=.0
    sc_dn=.0
    stcum1d=.0
    dn01d=.0
    sc_up_c=.0
    sc_dn_c=.0
    henry_coef=.0
    pw_up=.0
    pw_dn=0.

  end subroutine zero_tconv

end module mem_tconv
