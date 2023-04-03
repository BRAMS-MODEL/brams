!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine recycle()

  use dump, only: &
    dumpMessage
  use mem_grid, only: &
       ngrids, nnxp, nnyp, nzg, npatch, nzs, nnzp
  use var_tables, only : vtab_r, num_var
  use io_params, only: pastfn
  !srf -for carma AOT recycle
  use mem_aerad, only: nwave !INTENT(IN)
  use node_mod, only: &
       mzp, mxp, myp  ! INTENT(IN)
  use an_header,only: &
    anal_table,     &
    nvbtab
  use node_mod, only: &
       mzp, mxp, myp,  & ! INTENT(IN)
       izu, jzv,       & ! INTENT(IN)
       mynum,          & ! INTENT(IN)
       ibcon,          & ! INTENT(IN)
       nmachs,         &   ! INTENT(IN)
       mchnum, &
       master_num, &
       nodemxp, &
       nodemyp, &
       nodemzp, &
       nodei0, &
       nodej0
  USE ReadBcst, ONLY: &
       Broadcast
  use ParLib, only: &
       parf_bcast


  implicit none

  include "files.h"
  include "i8.h"
  include "constants.f90"
  integer, parameter :: i64 = selected_int_kind(14) !Kind for 64-bits Integer Numbers

  character(len=f_name_length) :: flnm
  !integer :: idtype,ng,nvars,lenf
  integer :: ng,nvars,lenf
  integer :: ierr
  real, allocatable :: scr1(:)
  real, allocatable :: scr2(:)
  character(len=f_name_length) :: flng
  integer(kind=i8)             :: npts
  integer(kind=i8)             :: fPosition
  character(len=*), parameter :: h="**(recycle)**"
  INTEGER :: npoints,i
  TYPE scridim
      real, allocatable :: scr(:,:,:,:)
  END TYPE scridim
  TYPE (scridim), DIMENSION(7) :: srcRead
  INTEGER :: ia,iz,ja,jz
  INTEGER :: m1,m2,m3
  INTEGER(kind=i64)  :: nzl,nxl,nyl,n4,j,k
  character(len=256) :: ctlFileName
  character(len=256) :: graFileName
  character :: cgrid
  character(len=256) :: lixo,dado
  integer :: xdef,ydef,zdef,recordLen,recn,tdef
  real :: undef
  real :: levels(45)


  !call readChemRecycleVars()
  !return

  allocate (srcRead(2)%scr(1,nnxp(1),nnyp(1),1))
  allocate (srcRead(3)%scr(nnzp(1),nnxp(1),nnyp(1),1))
  allocate (srcRead(4)%scr(nzg,nnxp(1),nnyp(1),npatch))
  allocate (srcRead(5)%scr(nzs,nnxp(1),nnyp(1),npatch))
  allocate (srcRead(6)%scr(1,nnxp(1),nnyp(1),npatch))
  allocate (srcRead(7)%scr(1,nnxp(1),nnyp(1),nwave))

  flnm=pastfn(1:len_trim(pastfn)-9)

  if(mchnum==master_num) then
    write (*,*) '+=====================================================================+'
    write (*,*) '|         Reading assimilation fields from past analysis file         |'
    write (*,*) '|Past File Name: ',trim(pastfn)
    write (*,*) '+=====================================================================+'
  endif

  call rams_read_header(flnm(1:len_trim(flnm)))

  ia = nodei0(mynum,1)+1
  iz = nodei0(mynum,1)+nodemxp(mynum,1)
  ja = nodej0(mynum,1)+1
  jz = nodej0(mynum,1)+nodemyp(mynum,1)
  m1=nodemzp(mynum,1)
  m2=nodemxp(mynum,1)
  m3=nodemyp(mynum,1)

  do ng=1,ngrids

     DO nvars=1,num_var(ng)

        if(vtab_r(nvars,ng)%irecycle == 1) then

           if(mchnum==master_num) &
           write (*,*) 'Reading assimilation field (from vtab):', vtab_r(nvars,ng)%name &
                , ' for grid:', ng                                     &
                , ' dim:', vtab_r(nvars,ng)%idim_type                  &
                , ' npts:', vtab_r(nvars,ng)%npts

           if(mchnum==master_num) then
              call FindFieldInAnalysisFile(vtab_r(nvars,ng)%name, ng,   &
                flnm(1:len_trim(flnm)), flng, npts, fPosition)
              npoints=npts
           endif

           CALL Broadcast(npoints, master_num, "npoints")
           npts=npoints

           allocate(scr1(npts), stat=ierr)
           if (ierr /= 0) & !call fatal_error(h//" allocating scr1")
            iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
              "error allocating scr1")
           allocate(scr2(npts), stat=ierr)
           if (ierr /= 0) &!call fatal_error(h//" allocating scr2")
           iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
             "error allocating scr2")
           IF(mchnum==master_num) call GetFieldInAnalysisFile(flng, npts, fPosition, scr1, scr2)

           CALL Broadcast(scr1, master_num, "scr1")

           if(vtab_r(nvars,ng)%idim_type == 4) then

              call unarrange_p(nnxp(ng),nnyp(ng),nzg,npatch  &
                   ,scr1(1),srcRead(4)%scr)

           elseif(vtab_r(nvars,ng)%idim_type == 5) then
              call unarrange_p(nnxp(ng),nnyp(ng),nzs,npatch  &
                   ,scr1(1),srcRead(5)%scr)
              !srf
              !use this for 3d atmospheric fields
           elseif(vtab_r(nvars,ng)%idim_type == 3) then
              call unarrange(nnzp(ng),nnxp(ng),nnyp(ng)  &
                   ,scr1(1),srcRead(3)%scr(:,:,:,1))

           elseif(vtab_r(nvars,ng)%idim_type == 6) then
              call rearrange_aot(npatch,nnxp(ng),nnyp(ng)  &
                   ,scr1(1),srcRead(6)%scr(1,:,:,:))
              !srf
              !use this for 3d (NX,NY,NWAVE) CARMA AOT fields
           elseif(vtab_r(nvars,ng)%idim_type == 7) then
              call rearrange_aot(nwave,nnxp(ng),nnyp(ng)  &
                   ,scr1(1),srcRead(7)%scr(1,:,:,:))

              !use this for 2 dim:
           else
              call atob(vtab_r(nvars,ng)%npts,  &
                   scr1(1),srcRead(2)%scr(1,:,:,1) )

           endif

           SELECT CASE (vtab_r(nvars,ng)%idim_type)
           CASE (2)
               nzl = 1
               nxl = nnxp(1)
               nyl = nnyp(1)
               n4 = 1
               call mk_2_buff(srcRead(vtab_r(nvars,ng)%idim_type)%scr(1,:,:,1), &
                              vtab_r(nvars,ng)%var_p, &
                              nnxp(ng), nnyp(ng), &
                              m2, m3, ia, iz, ja, jz)
            CASE (3)
               nzl = nnzp(1)
               nxl = nnxp(1)
               nyl = nnyp(1)
               n4 = 1
               call mk_3_buff(srcRead(vtab_r(nvars,ng)%idim_type)%scr(:,:,:,1), &
                              vtab_r(nvars,ng)%var_p, &
                              nnzp(ng),nnxp(ng), nnyp(ng), &
                              m1, m2, m3, ia, iz, ja, jz)
            CASE (4)
               nzl = nzg
               nxl = nnxp(1)
               nyl = nnyp(1)
               n4 = npatch
               call mk_4_buff(srcRead(vtab_r(nvars,ng)%idim_type)%scr, &
                              vtab_r(nvars,ng)%var_p, &
                              nzg,nnxp(ng), nnyp(ng),npatch, &
                              nzg, m2, m3,npatch, ia, iz, ja, jz)
            CASE (5)
               nzl = nzs
               nxl = nnxp(1)
               nyl = nnyp(1)
               n4 = npatch
               call mk_4_buff(srcRead(vtab_r(nvars,ng)%idim_type)%scr, &
                              vtab_r(nvars,ng)%var_p, &
                              nzs,nnxp(ng), nnyp(ng),npatch, &
                              nzs, m2, m3,npatch, ia, iz, ja, jz)
            CASE (6)
               nzl = 1
               nxl = nnxp(1)
               nyl = nnyp(1)
               n4 = npatch
               call mk_4_buff(srcRead(vtab_r(nvars,ng)%idim_type)%scr(1,:,:,:), &
                              vtab_r(nvars,ng)%var_p, &
                              1,nnxp(ng), nnyp(ng),npatch, &
                              1, m2, m3,npatch, ia, iz, ja, jz)
            CASE (7)
               nzl = 1
               nxl = nnxp(1)
               nyl = nnyp(1)
               n4 = nwave
               call mk_4_buff(srcRead(vtab_r(nvars,ng)%idim_type)%scr(1,:,:,:), &
                              vtab_r(nvars,ng)%var_p, &
                              1,nnxp(ng), nnyp(ng),nwave, &
                              1, m2, m3,nwave, ia, iz, ja, jz)
            CASE DEFAULT
               PRINT *, 'Wrong idim_type: ',vtab_r(nvars,ng)%idim_type
               STOP 'history_start'
            END SELECT

           deallocate(scr1)
           deallocate(scr2)

        endif

     enddo

  enddo
  if(mchnum==master_num) then
    write (*,*) '+=====================================================================+'
    write (*,*) '|                          End of recycle                             |'
    write (*,*) '+=====================================================================+'
  endif
  return
end subroutine recycle

!******************************************************************

subroutine rearrange_aot(nwave, nx, ny, array_i, array_o)
   implicit none
   integer, intent(in) :: nx, ny, nwave
   real, intent(in)    :: array_i(nx,ny,nwave)
   real, intent(out)   :: array_o(nx,ny,nwave)

   integer :: i, j, w

   do w=1, nwave
      do j=1, ny
         do i=1, nx
	    array_o(i,j,w) = array_i(i,j,w)
	    !if(w==12 .and. array_o(i,j,w) > 0.01) print*,i,j,array_o(i,j,w)
	 enddo
      enddo
   enddo
!       open(19,file='aot.gra',         &
!            form='unformatted',access='direct',status='unknown',  &
!            recl=4*nx*ny)
!       nrec=1
!       write(19,rec=nrec) array_o(:,:,12)
!       close (19)
!       stop 333

end subroutine rearrange_aot

!*******************************************************************************

subroutine unarrange_p(n2,n3,n4,n5,a,b)
implicit none

integer :: n2,n3,n4,n5
real :: a(n2,n3,n4,n5),b(n4,n2,n3,n5)

integer :: i,j,k,ip

do ip = 1,n5
   do k = 1,n4
      do j = 1,n3
         do i = 1,n2
            b(k,i,j,ip) = a(i,j,k,ip)
         enddo
      enddo
   enddo
enddo
return
end

! !=============================================================================================
! subroutine saveChemRecycleVars
!     !# salva as variaveis para recycle
!     !#
!     !# @note
!     !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
!     !#
!     !# **Brief**: salva as variaveis para recycle
!     !#
!     !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
!     !#
!     !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
!     !#
!     !# **Date**: 28 September 2020 (Monday)
!     !# @endnote
!     !#
!     !# @changes
!     !# &#9744; <br/>
!     !# @endchanges
!     !# @bug
!     !#
!     !#@endbug
!     !#
!     !#@todo
!     !#  &#9744; <br/>
!     !# @endtodo
!     !#
!     !# @warning
!     !# Now is under CC-GPL License, please see
!     !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
!     !# @endwarning
!     !#
    
!     !Use area
!     use dump
!     use var_tables, only : &
!         vtab_r, num_var
!     use node_mod, only: &
!        mzp, mxp, myp,  & ! INTENT(IN)
!        izu, jzv,       & ! INTENT(IN)
!        mynum,          & ! INTENT(IN)
!        ibcon,          & ! INTENT(IN)
!        nmachs,         &   ! INTENT(IN)
!        mchnum, &
!        master_num, &
!        nodemxp, &
!        nodemyp, &
!        nodemzp, &
!        nodei0, &
!        nodej0, &
!        ixb,ixe,iyb,iye

!     use mem_grid, only: &
!        ngrids, nnxp, nnyp, nzg, npatch, nzs, nnzp, &
!        oneGlobalGridData, &
!        iyear1,imonth1,idate1,ihour1,itime1 ! from RAMSIN
    
!     use chem1_list, only: &
!       nspecies_chem=>nspecies, &
!       spc_alloc_chem=>spc_alloc, &
!       spc_name_chem=>spc_name

!     use mem_chem1, only: &
!       chem1_g

!     use aer1_list, only: &
!       nspecies_aer=>nspecies, &
!       nmodes, &
!       mode_alloc_aer=>mode_alloc, &
!       aer_name

!     use mem_aer1, only: &
!       aer1_g

!     use ModDateUtils, only: &
!            date_add_to, date_add_to_dble

!     use ParLib

!     use mpi


!     implicit none

!     include "constants.f90"

!     character(len=*),parameter :: sourceName='recycle.f90' !Name of source code
!     character(len=*),parameter :: procedureName='**saveChemRecycleVars**' !Name of this procedure
!     !
!     !Local Parameters

!     !Input/Output variables

!     !Local variables
!     integer :: recordLen,irec
!     integer :: ng,nz,ispc,imode
!     real :: dlat,dlon
!     integer :: iyy,imm,idd,ihh,nvars,i,j,icnt,k,m
!     character(len=15) :: tDef
!     character(len=256) :: fileName
!     real, allocatable :: Globalvar(:,:,:),localVar(:)
!     integer :: gtag(nspecies_chem+nspecies_aer+nmodes,nmachs,nnzp(1))
!     type gv
!       integer :: ia
!       integer :: iz
!       integer :: ja
!       integer :: jz
!       real, pointer :: buffer(:)
!     end type gv
!     type(gv), allocatable :: getvar(:)

!     integer(kind=8) :: ncols
!     integer :: jnode,ia,iz,ja,jz,tag

!     integer(kind=8), allocatable :: getCols(:)
!     integer :: status(MPI_STATUS_SIZE)

!     !Code
!     ia = 2
!     iz = mxp-1
!     ja = 2
!     jz = myp-1


!     allocate(getCols(nmachs))
!     allocate(getvar(nmachs))
    
!     do jnode = 1,nmachs
!         !Alocando um pacote com todas as colunas

!         getVar(jnode)%ia=nodei0(jnode,1)+1
!         getVar(jnode)%iz=nodei0(jnode,1)+nodemxp(jnode,1)
!         getVar(jnode)%ja=nodej0(jnode,1)+1
!         getVar(jnode)%jz=nodej0(jnode,1)+nodemyp(jnode,1)
!         getcols(jnode)= (getVar(jnode)%iz-getVar(jnode)%ia+1)*(getVar(jnode)%jz-getVar(jnode)%ja+1)*nnzp(1)
!         allocate(getVar(jnode)%buffer(getCols(jNode)))

!         if(mchnum==master_num) write(*,"('  ',11i8)") jnode,ixb(jnode,1),ixe(jnode,1)  &
!              ,iyb(jnode,1),iye(jnode,1),getCols(jnode),getVar(jnode)%ia &
!              ,getVar(jnode)%iz,getVar(jnode)%ja,getVar(jnode)%jz

!      enddo
!      !print *,'Local: ',mynum,ia,iz,ja,jz,getCols(mynum)

!     if(mchnum==master_num) then

!       allocate(Globalvar(nnzp(1),nnxp(1),nnyp(1)))
!       GlobalVar=0.0

!       iErrNumber=dumpMessage(c_tty,c_yes,'','' &
!                ,c_notice,'Writing recycle files with chem and aerosols.')

!       dlon=oneGlobalGridData(1)%global_glon(2,1)-oneGlobalGridData(1)%global_glon(1,1)
!       dlat=oneGlobalGridData(1)%global_glat(1,2)-oneGlobalGridData(1)%global_glat(1,1)
  
!       call date_add_to_dble(iyear1,imonth1,idate1,0,dble(86400.0),'s' &
!                          ,iyy,imm,idd,ihh)
!       write(tDef,fmt='(I2.2,":00z",I2.2,A3,I4.4)')  ihh,idd,month_Name(imm),iyy
  
!       write(fileName,fmt='(A,I4.4,I2.2,I2.2,I2.2)') './recycle.',iyy,imm,idd,ihh
  
!       recordLen=4*nnxp(1)*nnyp(1)
  
!       open(unit=33,file=trim(fileName)//'.gra',action='WRITE',status='REPLACE' &
!                ,form='UNFORMATTED',access='DIRECT',recl=recordLen)

!       irec=1
!     endif

!     nvars=0

!     do ng=1,ngrids
!       do ispc=1,nspecies_chem
!         nvars=nvars+1
!         if(mchnum==master_num) then
!           Globalvar(:,getVar(mynum)%ia:getVar(mynum)%iz,getVar(mynum)%ja:getVar(mynum)%jz)=chem1_g(ispc,ng)%sc_p

!           do jnode = 1,nmachs
!             if(jnode/=mynum) then 
!               tag=jnode*1000
!               call parf_get_real(getVar(jnode)%buffer, getcols(jnode), jnode-1,tag)
!               call CopyBackLocal(getcols(jnode),mzp,getVar(jnode)%ia,getVar(jnode)%iz &
!                                                ,getVar(jnode)%ja,getVar(jnode)%jz &
!                                                ,getVar(jnode)%buffer &
!                                                ,Globalvar(:,getVar(jnode)%ia:getVar(jnode)%iz,getVar(jnode)%ja:getVar(jnode)%jz))
!             endif
!           enddo
!         else
!           tag=mynum*1000
!           call CopyLocal(getcols(mynum),mzp,1,mxp,1,myp,chem1_g(ispc,ng)%sc_p,getVar(mynum)%buffer)
!           print *,size(getVar(mynum)%buffer),getcols(mynum)
!           call parf_send_real(getVar(mynum)%buffer, getcols(mynum),master_num,tag)
!         endif
!         if(mchnum==master_num) then
!           do nz=1,nnzp(1)
!               write(33,rec=irec) Globalvar(nz,:,:)
!               irec=irec+1
!           enddo  
!         endif

!       enddo

!       do ispc=1,nspecies_aer
!         do imode=1,nmodes
!           if(mode_alloc_aer(imode,ispc) /= 1) cycle
!           nvars=nvars+1
!           if(mchnum==master_num) then
!             print *,imode,ispc,ng
!             Globalvar(:,getVar(mynum)%ia:getVar(mynum)%iz,getVar(mynum)%ja:getVar(mynum)%jz)=aer1_g(imode,ispc,ng)%sc_p
!             do jnode = 1,nmachs
!               if(jnode/=mynum) then 
!                 tag=jnode*1000
!                 call parf_get_real(getVar(jnode)%buffer, getcols(jnode), jnode-1,tag)
!                 call CopyBackLocal(getcols(jnode),mzp,getVar(jnode)%ia,getVar(jnode)%iz &
!                                                ,getVar(jnode)%ja,getVar(jnode)%jz &
!                                                ,getVar(jnode)%buffer &
!                                                ,Globalvar(:,getVar(jnode)%ia:getVar(jnode)%iz,getVar(jnode)%ja:getVar(jnode)%jz))
!               endif
!             enddo
!           else
!             tag=mynum*1000
!             call CopyLocal(getcols(mynum),mzp,1,mxp,1,myp,aer1_g(imode,ispc,ng)%sc_p,getVar(mynum)%buffer)
!             call parf_send_real(getVar(mynum)%buffer, getcols(mynum),master_num,tag)
!           endif
        
!           if(mchnum==master_num) then
!             do nz=1,nnzp(1)
!                 write(33,rec=irec) Globalvar(nz,:,:)
!                 irec=irec+1
!             enddo  
!           endif
  
!         enddo
!       enddo

!     enddo

!     if(mchnum==master_num) then  

!       close(unit=33)

!       open(unit=44, file=trim(fileName)//'.ctl', action='write', status='replace')
  
!       !writing the name of grads file
!       write(44,*) 'dset ^'//trim(fileName)//'.gra'
!       !writing others infos to ctl
!       write(44,*) 'undef -0.9990000E+34'
!       write(44,*) 'title Recycle data'
!       write(44,*) 'xdef ',nnxp(1),' linear ',oneGlobalGridData(1)%global_glon(1,1),dlon
!       write(44,*) 'ydef ',nnyp(1),' linear ',oneGlobalGridData(1)%global_glat(1,1),dlat
!       write(44,*) 'zdef ',nnzp(1),'levels',(i,i=1,nnzp(1))
!       write(44,*) 'tdef 1 linear '//tDef//' 1mo'
!       write(44,*) 'vars ',nvars
!       do ispc=1,nspecies_chem
!         write(44,*) trim(spc_name_chem(ispc)),nnzp(1),'99 ',trim(spc_name_chem(ispc))
!       enddo
!       do ispc=1,nspecies_aer
!         do imode=1,nmodes
!           if(mode_alloc_aer(imode,ispc) /= 1) cycle
!           write(44,*) trim(aer_name(imode,ispc)),nnzp(1),'99 ',trim(aer_name(imode,ispc))
!         enddo
!       enddo
  
!       write(44,*) 'endvars'

!       close(44)

!     endif

! end subroutine saveChemRecycleVars 

!=============================================================================================
subroutine CopyLocal(npts,nz,ia,iz,ja,jz,inVar,outVar)
    !# Copy the 2d dar to 1d var
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: copy 2d var to 1d var
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 29 September 2020 (Tuesday)
    !# @endnote
    !#
    !# @changes
    !# &#9744; <br/>
    !# @endchanges
    !# @bug
    !#
    !#@endbug
    !#
    !#@todo
    !#  &#9744; <br/>
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !# @endwarning
    !#
    
    !Use area
    use dump

    implicit none

    include "constants.f90"
    character(len=*),parameter :: sourceName='recycle.f90' !Name of source code
    character(len=*),parameter :: procedureName='**CopyLocal()**' !Name of this procedure
    !
    !Local Parameters

    !Input/Output variables
    integer, intent(in) :: ia,iz,ja,jz,nz
    integer(kind=8) :: npts
    real, intent(in) :: inVar(nz,ia:iz,ja:jz)
    real, intent(out):: outVar(npts)

    !Local variables
    integer :: Count
    integer :: i,j,iCount,k

    !Code
    count=(iz-ia+1)*(jz-ja+1)*nz

    if(npts/=count) then
      write (*,fmt='("(",5(I4.4,1X),") ",I8.8," /= ",I8.8)') nz,ia,iz,ja,jz,Count,npts
      iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
               ,c_fatal,'Size of two vars differ!')
    endif

    iCount=0
    do k=1,nz
      do i=ia,iz
        do j=ja,jz
          iCount=iCount+1
          outVar(iCount)=inVar(k,i,j)
        enddo
      enddo
    enddo

end subroutine CopyLocal

subroutine CopyBackLocal(npts,nz,ia,iz,ja,jz,inVar,outVar)
    !# Copy the 1d dar to 2d var
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: copy 1d var to 2d var
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 29 September 2020 (Tuesday)
    !# @endnote
    !#
    !# @changes
    !# &#9744; <br/>
    !# @endchanges
    !# @bug
    !#
    !#@endbug
    !#
    !#@todo
    !#  &#9744; <br/>
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !# @endwarning
    !#
    
    !Use area
    use dump

    implicit none

    include "constants.f90"
    character(len=*),parameter :: sourceName='recycle.f90' !Name of source code
    character(len=*),parameter :: procedureName='**CopyBackLocal()**' !Name of this procedure
    !
    !Local Parameters

    !Input/Output variables
    integer, intent(in) :: ia,iz,ja,jz,nz
    integer(kind=8) :: npts
    real, intent(out) :: OutVar(nz,ia:iz,ja:jz)
    real, intent(in):: inVar(npts)

    !Local variables
    integer :: Count
    integer :: i,j,iCount,k

    !Code
    count=(iz-ia+1)*(jz-ja+1)*nz

    if(npts/=count) then
      write (*,fmt='("(",5(I4.4,1X),") ",I8.8," /= ",I8.8)') nz,ia,iz,ja,jz,Count,npts
      iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
               ,c_fatal,'Size of two vars differ!')
    endif

    iCount=0
    do k=1,nz
      do i=ia,iz
        do j=ja,jz
          iCount=iCount+1
          outVar(k,i,j)=inVar(iCount)
        enddo
      enddo
    enddo

end subroutine CopyBackLocal


!=============================================================================================
subroutine readChemRecycleVars
    !# Le as variaveis para recycle
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: Le as variaveis para recycle
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 29 September 2020 (Monday)
    !# @endnote
    !#
    !# @changes
    !# &#9744; <br/>
    !# @endchanges
    !# @bug
    !#
    !#@endbug
    !#
    !#@todo
    !#  &#9744; <br/>
    !# @endtodo
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !# @endwarning
    !#
    
    !Use area
    use dump
    use var_tables, only : &
        vtab_r, num_var
    use node_mod, only: &
       mzp, mxp, myp,  & ! INTENT(IN)
       izu, jzv,       & ! INTENT(IN)
       mynum,          & ! INTENT(IN)
       ibcon,          & ! INTENT(IN)
       nmachs,         &   ! INTENT(IN)
       mchnum, &
       master_num, &
       nodemxp, &
       nodemyp, &
       nodemzp, &
       nodei0, &
       nodej0, &
       ixb,ixe,iyb,iye

    use mem_grid, only: &
       ngrids, nnxp, nnyp, nzg, npatch, nzs, nnzp, &
       oneGlobalGridData, &
       iyear1,imonth1,idate1,ihour1,itime1 ! from RAMSIN
    
    use chem1_list, only: &
      nspecies_chem=>nspecies, &
      spc_alloc_chem=>spc_alloc, &
      spc_name_chem=>spc_name

    use mem_chem1, only: &
      chem1_g

    use aer1_list, only: &
      nspecies_aer=>nspecies, &
      nmodes, &
      mode_alloc_aer=>mode_alloc, &
      aer_name

    use mem_aer1, only: &
      aer1_g

    use ModDateUtils, only: &
           date_add_to, date_add_to_dble

    use ParLib

    use ReadBcst, ONLY: &
       Broadcast

    implicit none

    include "constants.f90"
    character(len=*),parameter :: sourceName='recycle.f90' !Name of source code
    character(len=*),parameter :: procedureName='**leChemRecycleVars**' !Name of this procedure
    !
    !Local Parameters

    !Input/Output variables

    !Local variables
    integer :: recordLen,irec
    integer :: ng,nz,ispc,imode
    real :: dlat,dlon
    integer :: iyy,imm,idd,ihh,nvars,i,j,icnt,k,m
    character(len=15) :: tDef
    character(len=256) :: fileName
    real, allocatable :: Globalvar(:,:),localVar(:)
    integer :: gtag(nspecies_chem+nspecies_aer+nmodes,nmachs,nnzp(1))
    type gv
      integer :: ia
      integer :: iz
      integer :: ja
      integer :: jz
      real, pointer :: buffer(:)
    end type gv
    type(gv), allocatable :: getvar(:)

    integer(kind=8) :: ncols
    integer :: jnode,ia,iz,ja,jz,tag

    integer(kind=8), allocatable :: getCols(:)

    !Code
    ia = 2
    iz = mxp-1
    ja = 2
    jz = myp-1
    

    allocate(getCols(nmachs))
    allocate(getvar(nmachs))
    
    do jnode = 1,nmachs
        getcols(jnode)= (1+ixe(jnode,1)-ixb(jnode,1))  &
             *(1+iye(jnode,1)-iyb(jnode,1))

        allocate(getVar(jnode)%buffer(getCols(jNode)))
        getVar(jnode)%ia=nodei0(jnode,1)+2
        getVar(jnode)%iz=nodei0(jnode,1)+nodemxp(jnode,1)-1
        getVar(jnode)%ja=nodej0(jnode,1)+2
        getVar(jnode)%jz=nodej0(jnode,1)+nodemyp(jnode,1)-1

        !if(mchnum==master_num) write(*,"('Jnode=  ',11i8)") jnode,ixb(jnode,1),ixe(jnode,1)  &
        !     ,iyb(jnode,1),iye(jnode,1),getCols(jnode),getVar(jnode)%ia &
        !     ,getVar(jnode)%iz,getVar(jnode)%ja,getVar(jnode)%jz

     enddo
     !print *,'Local: ',mynum,ia,iz,ja,jz,getCols(mynum)

    allocate(Globalvar(nnxp(1),nnyp(1)))
    allocate(localVar(nnxp(1)*nnyp(1)))
    GlobalVar=0.0

    if(mchnum==master_num) then
 
      write(fileName,fmt='(A,I4.4,I2.2,I2.2,I2.2)') './recycle.',iyear1,imonth1,idate1,ihour1
  
      recordLen=4*nnxp(1)*nnyp(1)

      iErrNumber=dumpMessage(c_tty,c_yes,'','' &
               ,c_notice,'Reading recycle files with chem and aerosols: '//trim(fileName))
  
      open(unit=33,file=trim(fileName)//'.gra',action='read',status='old' &
               ,form='UNFORMATTED',access='DIRECT',recl=recordLen)

      irec=1
    endif

    nvars=0

    do ng=1,ngrids
      do ispc=1,nspecies_chem
        nvars=nvars+1
        
        do nz=1,nnzp(1)

          if(mchnum==master_num) then
            read(33,rec=irec) Globalvar
            irec=irec+1
            icnt=0
            do i=1,nnxp(1)
              do j=1,nnyp(1)
                icnt=icnt+1
                localVar(icnt)=Globalvar(i,j)
              enddo
            enddo
          endif

          CALL Broadcast(localVar, master_num, "recycle")

          if(mchnum/=master_num) then
            icnt=0
            do i=1,nnxp(1)
              do j=1,nnyp(1)
                icnt=icnt+1
                Globalvar(i,j)=localVar(icnt)
              enddo
            enddo
          endif 

          chem1_g(ispc,ng)%sc_p(nz,ia:iz,ja:jz)=Globalvar(getVar(mynum)%ia:getVar(mynum)%iz,getVar(mynum)%ja:getVar(mynum)%jz)

        enddo
      enddo

      do ispc=1,nspecies_aer
        do imode=1,nmodes
          if(mode_alloc_aer(imode,ispc) /= 1) cycle
          nvars=nvars+1
          do nz=1,nnzp(1)

            if(mchnum==master_num) then
              read(33,rec=irec) Globalvar
              irec=irec+1
              icnt=0
              do i=1,nnxp(1)
                do j=1,nnyp(1)
                  icnt=icnt+1
                  localVar(icnt)=Globalvar(i,j)
                enddo
              enddo
            endif

            CALL Broadcast(localVar, master_num, "recycle")

            if(mchnum/=master_num) then
              icnt=0
              do i=1,nnxp(1)
                do j=1,nnyp(1)
                  icnt=icnt+1
                  Globalvar(i,j)=localVar(icnt)
                enddo
              enddo
            endif 

            aer1_g(imode,ispc,ng)%sc_p(nz,ia:iz,ja:jz)=Globalvar(getVar(mynum)%ia:getVar(mynum)%iz,getVar(mynum)%ja:getVar(mynum)%jz)

          enddo
        enddo
      enddo
    enddo

end subroutine readChemRecycleVars 
   
   