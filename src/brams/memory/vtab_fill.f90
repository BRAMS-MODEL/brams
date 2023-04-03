!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


!-------------------------------------------------------------------------

subroutine lite_varset(proc_type)

  use var_tables, only: &
       nvgrids, & ! INTENT(IN)
       vtab_r,  & ! INTENT(INOUT)
       num_var    ! INTENT(IN)

  use io_params, only: &
       nlite_vars, & ! INTENT(IN)
       lite_vars     ! INTENT(IN)

  implicit none

  ! Arguments:
  integer, intent(in) :: proc_type

  ! Local variables:
  integer :: nv,ng,nvl,ifound


  ! Loop over each variable input in namelist "LITE_VARS" and set
  !   lite flag in var_tables

  do ng = 1,nvgrids   
     vtab_r(1:num_var(ng),ng)%ilite = 0
  enddo

  do nvl=1,nlite_vars
     ifound=0

     do ng=1,nvgrids

        do nv=1,num_var(ng)

           if (vtab_r(nv,ng)%name == lite_vars(nvl) ) then
              vtab_r(nv,ng)%ilite = 1
              ifound=1
           endif

        enddo

     enddo

     if (proc_type==0 .or. proc_type==1) then !Output only in Master Process
        if(ifound == 0) then
           print*,'!---------------------------------------------------------'
           print*,'! LITE_VARS variable does not exist in main variable table'
           print*,'!    variable name-->',lite_vars(nvl),'<--'
           print*,'!---------------------------------------------------------'
        else
           print*,'!---------------------------------------------------------'
           print*,'! LITE_VARS variable added--->',trim(lite_vars(nvl))
           print*,'!---------------------------------------------------------'
        endif
     endif

  enddo

  return
end subroutine lite_varset

!-------------------------------------------------------------------------

subroutine vtables_scalar(varp,vart,ng,tabstr)

  use var_tables

  implicit none
  real, target :: varp,vart
  integer, intent(in) :: ng
  character (len=*), intent(in) :: tabstr

  character (len=80) ::line
  character (len=1) ::toksep=':'
  character (len=32) ::tokens(10)
  character (len=16) :: cname,ctab

  integer :: ntok,nv,isnum,ns

  call tokenize1(tabstr,tokens,ntok,toksep)
  cname=tokens(1)
  !   ctab=tokens(2) 

  !    See if this scalar name is already in the table...

  isnum=0

  !   do ns=1,num_scalar(ng)
  !      if(cname == scalar_tab(ns,ng)%name) then
  !      	 isnum=ns
  !      	 exit
  !      endif
  !   enddo

  !    Fill in existing table slot or make new scalar slot

  !   if (isnum == 0) then
  num_scalar(ng)=num_scalar(ng)+1
  nv=num_scalar(ng)
  scalar_tab(nv,ng)%name = cname
  !   else
  !      nv=isnum
  !   endif

  scalar_tab(nv,ng)%var_p => varp
  scalar_tab(nv,ng)%var_t => vart

  !   if(ctab == 'sclp' ) then
  !      scalar_tab(nv,ng)%varp => varp
  !   elseif(ctab == 'sclt' ) then
  !      scalar_tab(nv,ng)%vart => vart
  !   else
  !      print*, 'Illegal scalar table specification for var:', cname
  !      stop 'bad scalar table'
  !   endif

  return
end subroutine vtables_scalar

!-------------------------------------------------------------------------

subroutine vtables_scalar_new(varp,vart,ng,tabstr,elements)

  use var_tables

  implicit none
  integer :: elements !ALF
  real, target :: varp(elements), vart(elements)
  integer, intent(in) :: ng
  character (len=*), intent(in) :: tabstr

  character (len=80) ::line
  character (len=1) ::toksep=':'
  character (len=32) ::tokens(10)
  character (len=16) :: cname,ctab

  integer :: ntok,nv,isnum,ns, i

  call tokenize1(tabstr,tokens,ntok,toksep)
  cname=tokens(1)
  !   ctab=tokens(2) 

  !    See if this scalar name is already in the table...

  isnum=0

  !   do ns=1,num_scalar(ng)
  !      if(cname == scalar_tab(ns,ng)%name) then
  !      	 isnum=ns
  !      	 exit
  !      endif
  !   enddo

  !    Fill in existing table slot or make new scalar slot

  !   if (isnum == 0) then
  !num_scalar(ng)=num_scalar(ng)+1
  nv=num_scalar(ng)
  !scalar_tab(nv,ng)%name = cname
  !   else
  !      nv=isnum
  !   endif

  !scalar_tab(nv,ng)%var_p => varp
  !scalar_tab(nv,ng)%var_t => vart

  ! ALF
  scalar_tab(nv,ng)%a_var_p => varp(:)
  scalar_tab(nv,ng)%a_var_t => vart(:)

  !
  !

  !   if(ctab == 'sclp' ) then
  !      scalar_tab(nv,ng)%varp => varp
  !   elseif(ctab == 'sclt' ) then
  !      scalar_tab(nv,ng)%vart => vart
  !   else
  !      print*, 'Illegal scalar table specification for var:', cname
  !      stop 'bad scalar table'
  !   endif

  return
end subroutine vtables_scalar_new






subroutine GetVarFromMem (nxp, nyp, nzp, nzg, nzs, npatch, &
     varName, itype, ngrd, arrayOut, sizeArray)
  use var_tables, only: vtab_r, num_var

  implicit none
  include "i8.h"
  integer,            intent(in)    :: nxp  ! as at vartable
  integer,            intent(in)    :: nyp  ! as at vartable
  integer,            intent(in)    :: nzp  ! as at vartable
  integer,            intent(in)    :: nzg  ! as at vartable
  integer,            intent(in)    :: nzs  ! as at vartable
  integer,            intent(in)    :: npatch  ! as at vartable
  character(LEN=*),   intent(in)    :: varName
  integer,            intent(in)    :: ngrd
  integer,            intent(out)   :: itype
  integer(kind=i8),   intent(in)    :: sizeArray
  real,	              intent(inout) :: arrayOut(sizeArray)

  character(len=16) :: c0, c1
  character(len=*), parameter :: h="**(GetVarFromMem)**"
  integer(kind=i8) :: ni
  integer(kind=i8) :: npts
  logical          :: found
  character(len=len(varName)) :: varnIn, varnOut
  real, pointer :: ptr ! points to a field at vartable
  real :: scr1(sizeArray)

  ! field name changes from vartable to analysis file
  ! in two cases (PI and HKH); 
  ! given output file field name, find vartable correspondent

  if (trim(varName) == 'PI') then
     varnIn = 'PP'
  else if (trim(varName) == 'HKH') then
     varnIn = 'HKM'
  else
     varnIn = varName
  end if

  ! search for vartable name at vartable
  ! store result at array 

  found = .false.
  do ni = 1, num_var(ngrd)
     if (trim(vtab_r(ni,ngrd)%name) == trim(varnIn)) then
        itype = vtab_r(ni,ngrd)%idim_type
        npts  = vtab_r(ni,ngrd)%npts
        ptr => vtab_r(ni,ngrd)%var_p
        if (npts > sizeArray) then
           write(c0,"(i16)") npts
           write(c1,"(i16)") sizeArray
           call fatal_error(h//&
                " array size for "//trim(varName)//&
                " is "//trim(adjustl(c1))//&
                ", smaller than "//trim(adjustl(c0))//" required")
        end if
        found = .true.
        exit
     end if
  end do

  ! halts if not there

  if (.not. found) then
     write(*,"(a)") h//" var "//trim(varnIn)//" not found in vtab_r; will dump vtab_r"
     call DumpVTab(ngrd)
     call fatal_error(h//" var "//trim(varnIn)//" not found in vtab_r")
  end if

  ! convert fields PP, HKM and VKH from vartables to analysis file
  ! or store field at scr1

  call PreProcForOutput(ngrd, varnIn, npts, ptr, scr1, varnOut)

  ! verify output name

  if (trim(varnOut) /= trim (varName)) then
     call fatal_error(h//" fails computing "//trim(varName))
  end if

  ! move verticals from first to third dimension, if required

  if (itype==3 .or. itype==4 .or. itype==5) then
     call RearrangeForOutput(nxp, nyp, nzp, nzg, nzs, npatch, &
          itype, scr1, arrayOut)
  else
     arrayOut = scr1
  end if
end subroutine GetVarFromMem

subroutine DumpVTab(ngrd)
  use var_tables, only: vtab_r, num_var

  implicit none
  include "i8.h"
  integer, intent(in) :: ngrd

  integer :: i
  character(len=*), parameter :: h="**(DumpVTab)**"

  write(*,"(a,i2)") h//" dump of vtab_r names for grid ",ngrd
  write(*,"(a)") h//" name            idim_type           npts"
  do i = 1, num_var(ngrd)
     write(*,"(1x,a16,1x,i8,1x,i16)") vtab_r(i,ngrd)%name, &
          vtab_r(i,ngrd)%idim_type, vtab_r(i,ngrd)%npts
  end do
end subroutine DumpVTab
  
  subroutine setInitial4Vtable(ng)
  
  use chem1_list, only: chem_name=>spc_name,    &
                        chem_alloc=>spc_alloc,  & 
			chem_on=>on,            &
			chem_fdda=>fdda,        &
			chem_transport=>transport, &
			chem_nspecies=>nspecies 
			
  use aer1_list, only: aer_name=>spc_name,     &
                       aer_nspecies=>nspecies, &
		       aer_alloc=>spc_alloc,     &
		       aer_fdda=>fdda, &
		       aer_transport=>transport, &
		       aer_on=>on,               &
		       aer_nmodes=>nmodes
			
  use mem_chem1, only: CHEM_ASSIM, CHEMISTRY
  use mem_aer1, only: AEROSOL, AER_ASSIM
  use var_tables, only: vtab_r, num_var
  
    integer,           intent(in) :: ng
    integer :: ni, nspc, imode
    character(len=2) :: cmode

    do ni = 1, num_var(ng)
    	
	vtab_r(ni,ng)%ianal=0	
	
	if (trim(vtab_r(ni,ng)%name) == 'TOPT'  .or. &
	    trim(vtab_r(ni,ng)%name) == 'UP'    .or. &
	    trim(vtab_r(ni,ng)%name) == 'VP'    .or. &
	    trim(vtab_r(ni,ng)%name) == 'THETA' .or. &
	    trim(vtab_r(ni,ng)%name) == 'PP'    .or. &
	    trim(vtab_r(ni,ng)%name) == 'RV') then
		vtab_r(ni,ng)%ianal=1
		cycle
        end if
	
	if(vtab_r(ni,ng)%irecycle == 1) vtab_r(ni,ng)%ianal=1
	
	
	if(CHEMISTRY >= 0) then 
		do nspc=1,chem_nspecies
			!print*, spc_alloc(fdda,nspc), on
			!print*, trim(vtab_r(ni,ng)%name), '>>', trim(spc_name(nspc))//'P'
        		if(chem_alloc(chem_fdda,nspc) == chem_on .and. trim(vtab_r(ni,ng)%name) == trim(chem_name(nspc))//'P') then 
				vtab_r(ni,ng)%ianal=1
				cycle
			end if
		end do
	end if
	if(AEROSOL == 1 .and. CHEMISTRY >= 0) then
		do nspc=1,aer_nspecies
			do imode = 1, aer_nmodes
			write(cmode, '(BN, I2)')imode
			cmode = adjustl(cmode)
				!print*, trim(vtab_r(ni,ng)%name), trim(aer_name(nspc))//trim(cmode)//'P'
				if(aer_alloc(aer_fdda,imode,nspc) == 1  .and. trim(vtab_r(ni,ng)%name) == trim(aer_name(nspc))//trim(cmode)//'P') then
					vtab_r(ni,ng)%ianal=1
					cycle
				end if
			end do
		end do
	end if
    end do

  end subroutine setInitial4Vtable
