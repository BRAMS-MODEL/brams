program leQueima
   use dump

   implicit none

   integer, parameter :: iunit=22

   include "constants.f90"
   integer :: nxp,nyp
   real :: dep_glon(2),dep_glat (2)
   character(len=9) :: date
   character(len=32) :: chemical_mechanism,aerosol_mechanism,section,read_spc_name
   character(len=32) :: read_src_name,read_units
   real, allocatable :: var(:,:,:)
   character(len=32) :: varName(1000),unitName(1000)
   character(len=256) :: cArg


   integer :: exdo,read_ident_chem_mec,read_ident_src,vCount,read_ident_aer,read_aer_mode
   integer :: nsp,recl,nrec,i
   
   i=1
   call get_command_argument(i, cArg)
   if(trim(cArg)=='-h' .or. trim(cArg)=='--help') then
      call help()
      return
      !
   endif

   open(unit=iunit,file=trim(cArg)//'.vfm',form='formatted',status='old') 
   
   read(iunit,*) nxp,(dep_glon(i),i=1,2)
   read(iunit,*) nyp,(dep_glat(i),i=1,2)
   read(iunit,*) date 
   write (*,fmt='(a)') '=== sources header (nxpoints,nypoint,lon,lat,date) ==='
   write (*,fmt='(4x,2(i4.4,1x),4(f8.3,1x),a)') &
   nxp,nyp,(dep_glon(i),i=1,2),(dep_glat(i),i=1,2),date
   read(iunit,*)  chemical_mechanism,aerosol_mechanism

   allocate(var(1000,nxp,nyp))
   
   print*,'   chem mechanism= ',trim(chemical_mechanism)
   print*,'   aer  mechanism= ',trim(aerosol_mechanism)

   vCount=0

   do while(.true.)

      read(iunit,*,iostat=exdo) section
      if (exdo<0) exit ! eof

      !- emission section  --------------------------------
      if(trim(section) == 'chemistry') then
         read(iunit,*)   read_spc_name &
                           , read_ident_chem_mec &
                           , read_src_name       &
                           , read_ident_src      &
                           , read_units

         print *,trim(read_spc_name) &
                           , read_ident_chem_mec &
                           , trim(read_src_name)       &
                           , read_ident_src      &
                           , trim(read_units)

         vCount=vCount+1
         varName(vCount)=trim(read_spc_name)//'_'//trim(read_src_name)
         unitName(vCount)=trim(read_units)

         CALL vfirec(iunit,var(vCount,:,:),nxp*nyp,'LIN')

         print *,vCount,maxval(var(vCount,:,:)),minVal(var(vCount,:,:))

      elseif(trim(section) == 'aerosol' ) then 
       
         read(iunit,*)   read_spc_name       &
                       , read_ident_aer      &
                       , read_aer_mode       &
                       , read_src_name       &
                       , read_ident_src      &
                       , read_units
  
         vCount=vCount+1
         varName(vCount)=trim(read_spc_name)//'_'//trim(read_src_name)
         unitName(vCount)=trim(read_units)
         
         print *, trim(read_spc_name)       &
                       , read_ident_aer      &
                       , read_aer_mode       &
                       , trim(read_src_name)       &
                       , read_ident_src      &
                       , trim(read_units)

         CALL vfirec(iunit,var(vCount,:,:),nxp*nyp,'LIN')

      endif

   enddo

   recl=nxp*nyp*4

   open(unit=33,file=trim(cArg)//'.ctl' &
            ,action='WRITE',status='replace',form='FORMATTED')

   !writing the name of grads file
   write(33,*) 'dset ^'//trim(cArg)//'.gra'
   !writing others infos to ctl
   write(33,*) 'undef -0.9990000E+34'
   write(33,*) 'title correlacao'
   write(33,*) 'xdef ',nxp,' linear ',(dep_glon(i),i=1,2)
   write(33,*) 'ydef ',nyp,' linear ',(dep_glat(i),i=1,2)
   write(33,*) 'zdef ',1,'levels',1000
   write(33,*) 'tdef 1 linear 00:00z'//date,' 1hr'
   write(33,*) 'vars ',vCount
   do nsp=1,vCount
      write(33,*) varName(nsp),1,'99 ',' '//trim(unitName(nsp))
   enddo
   write(33,*) 'endvars'

   close(unit=33)
   nrec=0

   open(unit=33,file=trim(cARg)//'.gra',access='direct',recl=recl,status='replace')
   do nsp=1,vCount
      nrec=nrec+1
      write(33,rec=nrec) var(nsp,:,:)
   enddo

   close(unit=33)

end program leQueima

subroutine help()

   write(*,*) ''
   write(*,*) ' lq_1.0 FILE, where: '
   write(*,*) ''
   write(*,*) '  FILE - Filename of Queima SRC file without extension'
   write(*,*) ''
   write(*,*) 'Example:'
   write(*,*) ''
   write(*,*) 'lq_1.0 Queima_source-T-2020-09-18-000000-g1'
   write(*,*) ''

end subroutine help
