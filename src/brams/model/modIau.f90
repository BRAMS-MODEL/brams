module modIau
    !######################################################################
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: implements the incremental analysis update - IAU 
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Authors**: Saulo R. Freitas & Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 25 February 2021
    !# @endnote
    !#
    !# @changes
    !#
    !# +
    !# @endchanges
    !# @bug
    !#
    !#@endbug
    !#
    !# @warning
    !# Now is under CC-GPL License, please see
    !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
    !# See the message in todo
    !# @endwarning
    !#######################################################################
    !---
    use dump, only: &
        dumpMessage     !subroutine
      
    use ModNamelistFile, only: namelistFile
        
    use mem_grid, only:                  &
        iyear1, imonth1, itime1,  idate1,   & ! all intent in
        time, jdim, grid_g, nnzp, nnxp,nnyp,&
        ztop, zt,dtlt, npatch

        
    use node_mod, only:     &
        i0,            &
        j0,            &
        mchnum,        &
        master_num,    &
        nmachs,        &
        mxp,           &
        myp,           &
        mzp,           &
        nodei0,        &
        nodej0,        &
        nodemxp,       &
        nodemyp,       &
        nodemzp,       &
        mynum
                                     
    use mem_varinit, only: &
        nudlat, tnudcent, tnudtop, tnudlat, znudtop
                 
    use ParLib, only: &
          parf_get_int, &
          parf_send_int, &
          parf_get_real, &
          parf_send_real


    implicit none
    include "files.h"
    include "i8.h"

    character(len=*),parameter :: sourceName='modIau.f90' !Name of source code

    type tend_vars_iau
        real   , pointer, dimension(:) :: ut,vt,tht,rtt,pt,iauwts
        real   , pointer, dimension(:) :: dfi
        integer                        :: istep
        integer                        :: nsteps
    end type
    type(tend_vars_iau), allocatable,dimension(:) :: tend_iau_g

    integer                      :: applyIAU      !-- indicates the Iau action (0,1,2 from RAMSIN)
    real                         :: timeWindowIAU !-- seconds 21600,from RAMSIN
    character(len=f_name_length) :: fileNameIAU   !-- prefix of file with the IAU tend, from RAMSIN
    real                         :: ramp          !-- seconds 3600, from RAMSIN (not in use)
   
    integer :: write_iau_tend_to_file  = 1
    integer :: read_iau_tend_from_file = 1
    integer :: get_iau_weights         = 1
    integer :: nsteps
    integer :: apply_digital_filter    = 1
    integer,parameter :: irecUT  = 1 &
        ,irecVT  = 2 &
        ,irecPT  = 3 &
        ,irecRTT = 4 &
        ,irecTHT = 5

    integer, parameter :: xFactor=10
    type cc 
      integer(kind=i8) :: size
      integer :: tag
      integer :: tia,tiz,tja,tjz
      real, pointer, dimension(:) :: varBuff
    end type cc
    type(cc), allocatable :: com(:)

    logical :: firstTimeIau
    real    :: real_byte_size = 1.0
    integer :: int_byte_size

contains


    !=======================================================================================
    subroutine allocIau(tend_iau_g,m1,m2,m3)
            !# Alocate iau variables
        !#
        !# @note
        !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
        !#
        !# **Brief**: allocate the variables for incremental analysis update
        !#
        !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
        !#
        !# **Author**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
        !#
        !# **Date**: 25 February 2021
        !# @endnote
        !#
        !# @changes
        !#
        !# +
        !# @endchanges
        !# @bug
        !#
        !#@endbug
        !#
        !#@todo
        !#  &#9744; Changes to work with another type of grads files<br/>
        !# @endtodo
        !#
        !# @warning
        !# Now is under CC-GPL License, please see
        !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
        !# See the message in todo
        !# @endwarning
        !#
        !#---
        character(len=*),parameter :: procedureName='**allocIau**' !Name of this procedure


        integer, intent(in) :: m2
        !# number of points in x direction
        integer, intent(in) :: m3
        !# number of points in y direction
        integer, intent(in) :: m1
        !# number of points in z direction
 
        type (tend_vars_iau) :: tend_iau_g
        integer :: ierr, ntpts
      
        ntpts = m1 * m2 * m3
        allocate (tend_iau_g%rtt(ntpts), STAT=ierr)
        if (ierr/=0) call fatal_error("Allocating tend_iau_g%rtt")
        tend_iau_g%rtt = 0.0
      
        allocate (tend_iau_g%ut(ntpts), STAT=ierr)
        if (ierr/=0) call fatal_error("Allocating tend_iau_g%rt")
        tend_iau_g%ut = 0.0
 
        allocate (tend_iau_g%vt(ntpts), STAT=ierr)
        if (ierr/=0) call fatal_error("Allocating tend_iau_g%rt")
        tend_iau_g%vt = 0.0

        allocate (tend_iau_g%pt(ntpts), STAT=ierr)
        if (ierr/=0) call fatal_error("Allocating tend_iau_g%rt")
        tend_iau_g%pt = 0.0

        allocate (tend_iau_g%tht(ntpts), STAT=ierr)
        if (ierr/=0) call fatal_error("Allocating tend_iau_g%rt")
        tend_iau_g%tht = 0.0

        allocate (tend_iau_g%iauwts(ntpts), STAT=ierr)
        if (ierr/=0) call fatal_error("Allocating tend_iau_g%iauwts")
        tend_iau_g%iauwts = 0.0
      
        firstTimeIau=.true.

    end subroutine allocIau

    !===========================================================================================

    subroutine nullifyIAU(tend_iau_g)
        implicit none
        type (tend_vars_iau) :: tend_iau_g
        
        if(associated(tend_iau_g%ut)     ) nullify (tend_iau_g%ut )
        if(associated(tend_iau_g%vt)     ) nullify (tend_iau_g%vt )
        if(associated(tend_iau_g%rtt)    ) nullify (tend_iau_g%rtt)
        if(associated(tend_iau_g%tht)    ) nullify (tend_iau_g%tht)
        if(associated(tend_iau_g%pt)     ) nullify (tend_iau_g%pt )
        if(associated(tend_iau_g%pt)     ) nullify (tend_iau_g%pt )
        if(associated(tend_iau_g%iauwts) ) nullify (tend_iau_g%iauwts )
        if(associated(tend_iau_g%dfi)    ) nullify (tend_iau_g%dfi )
        
    end subroutine nullifyIAU
    !===========================================================================================

    subroutine CreateIauTendency(ng,ntpts, m1, m2, m3,ia,iz,ja,jz&
        ,varup,varvp,varpp,vartp,varrp  &
        ,varuf,varvf,varpf,vartf,varrf  &
        ,up,vp,theta,rtp,pp)
        implicit none
        integer, intent(in) :: ng,ntpts, m1,m2,m3 ,ia,iz,ja,jz
        real,    intent(in), dimension(ntpts) ::         &
            varup,varvp,varpp,vartp,varrp  &
            ,varuf,varvf,varpf,vartf,varrf  &
            ,up,vp,theta,rtp,pp
        real :: dtinv
        character(len=f_name_length) :: sVarName
        character*1 cgrid
      
        if(write_iau_tend_to_file /= 1) return 

        if(firstTimeIau) call initComIau(ia,iz,ja,jz)
        
        ! --- time interval for the tendencies (must be in seconds)
        dtinv = 1./timeWindowIAU
        
        !- iau tend = ( ana - model)/timeWindowIAU 
        !- for an IAU procedure starting 21UTC and time widown of 6h,
        !- 'ana' is the analysis at time 00 UTC 
        !- (in the current configuration it is var*"F", and not var*"P".
        !-
        tend_iau_g(ng)%ut (1:ntpts) = (varuf (1:ntpts) - up   (1:ntpts))*dtinv
        tend_iau_g(ng)%vt (1:ntpts) = (varvf (1:ntpts) - vp   (1:ntpts))*dtinv
        tend_iau_g(ng)%pt (1:ntpts) = (varpf (1:ntpts) - pp   (1:ntpts))*dtinv
        tend_iau_g(ng)%rtt(1:ntpts) = (varrf (1:ntpts) - rtp  (1:ntpts))*dtinv
        tend_iau_g(ng)%tht(1:ntpts) = (vartf (1:ntpts) - theta(1:ntpts))*dtinv

        
        write(cgrid,'(i1)')ng
        call makefnam (sVarName, fileNameIAU, 0, iyear1, imonth1, idate1, itime1*100, 'V', 'g'//cgrid,'bin')        

        call writeReadIauTendency('write',ng, m1, m2, m3,ia,iz,ja,jz,sVarName,'iauUT' ,tend_iau_g(ng)%ut ,irecUT )
        call writeReadIauTendency('write',ng, m1, m2, m3,ia,iz,ja,jz,sVarName,'iauVT' ,tend_iau_g(ng)%vt ,irecVT )
        call writeReadIauTendency('write',ng, m1, m2, m3,ia,iz,ja,jz,sVarName,'iauPT' ,tend_iau_g(ng)%pt ,irecPT )
        call writeReadIauTendency('write',ng, m1, m2, m3,ia,iz,ja,jz,sVarName,'iauRTT',tend_iau_g(ng)%rtt,irecRTT)
        call writeReadIauTendency('write',ng, m1, m2, m3,ia,iz,ja,jz,sVarName,'iauTHT',tend_iau_g(ng)%tht,irecTHT)
          
        if(mchnum == master_num) then
                print*,"====================================================="
                print*," Dumping the IAU tendencies for the time period:"
                print*,  iyear1, imonth1, idate1, itime1/100, 'h to', nint(itime1/100+timeWindowIAU/3600.), 'h'
                print*," file:", sVarName(1:len_trim(sVarName))
                print*,"====================================================="
                call flush(6)
        endif
        write_iau_tend_to_file=2
        return

        !open(2, status='unknown', form='UNFORMATTED',file=sVarName(1:len_trim(sVarName)))
        !write(2) tend_iau_g(ng)%ut
        !write(2) tend_iau_g(ng)%vt
        !write(2) tend_iau_g(ng)%pt
        !write(2) tend_iau_g(ng)%rtT
        !write(2) tend_iau_g(ng)%tht
        !close(2)
          

    end subroutine CreateIauTendency
    !===========================================================================================
        
    subroutine GetIauTendency(ng,ntpts, m1, m2, m3,ia,iz,ja,jz,time,dtlt)
        use mem_tend, only: tend
        implicit none
        integer, intent(in) :: ng,ntpts, m1, m2, m3,ia,iz,ja,jz
        real   , intent(in) :: time,dtlt
        real :: centFact
        integer :: ierr
          
        if(time >  timeWindowIAU + RAMP) return
          
        !--- deallocate the arrays as they will not be longer needed after the IAU procedure
        !
        if( abs(time - (timeWindowIAU + RAMP)) < 0.01) then
          
            deallocate(tend_iau_g(ng)%ut,  tend_iau_g(ng)%vt ,tend_iau_g(ng)%pt     & 
                ,tend_iau_g(ng)%rtt, tend_iau_g(ng)%tht,tend_iau_g(ng)%iauwts )
            
            if(mchnum == master_num) then
                print*,"====================================================="
                print*,"===> IAU procedure ended"
                print*,"===> IAU tendencies deallocated at time=",time
                print*,"====================================================="
                call flush(6)
            endif
            
            return
        endif
          
        !---- initialize the IAU weights (lat, top, center)
        !
        if(get_iau_weights == 1) &
          
            call IAU_VariableWeight(nnzp(ng), nodemxp(mynum,ng), nodemyp(mynum,ng), nnxp(ng),&
            nnyp(ng), nodei0 (mynum,ng), nodej0 (mynum,ng)          ,&
            grid_g(ng)%topt, grid_g(ng)%rtgt, tend_iau_g(ng)%iauwts)
          
          
          
        !---- get the Digital Filter coefficients
        !
        tend_iau_g(ng)%nsteps = nint( timeWindowIAU/dtlt ) + 1
          
        if (.not.associated(tend_iau_g(ng)%dfi) .and. apply_digital_filter == 1) then
              
            allocate(tend_iau_g(ng)%dfi(tend_iau_g(ng)%nsteps), STAT=ierr)
            if (ierr/=0) call fatal_error("Allocating tend_iau_g%dfi")
              
            tend_iau_g(ng)%istep=0;  tend_iau_g(ng)%dfi(:)=0.0
              
            call get_dfi_coeffs(DTLT,timeWindowIAU,tend_iau_g(ng)%nsteps,tend_iau_g(ng)%dfi)
        endif
          
          
        !--- get the factor for the ramping after the IAU phase
                  !
        if(time < dtlt + timeWindowIAU) then
            centFact=1.0
        else
            centFact=-1*( (time-timeWindowIAU) /(RAMP)) + 1
        endif

        if(time <= dtlt + timeWindowIAU .and. mchnum == master_num) then
            write(*,100) "===> doing IAU in center domain              :",time,centFact
            call flush(6)
        endif
        if(time > dtlt + timeWindowIAU .and. mchnum == master_num) then
            write(*,100) "===> ramping the IAU in center domain to zero:",time,centFact
            call flush(6)
        endif
        
        
        !--- apply the digital filter
        !
        if( apply_digital_filter == 1 ) then
             
            tend_iau_g(ng)%istep = tend_iau_g(ng)%istep + 1
          
            centFact = centFact * tend_iau_g(ng)%dfi(tend_iau_g(ng)%istep)

            if(tend_iau_g(ng)%istep == tend_iau_g(ng)%nsteps -1 ) tend_iau_g(ng)%istep = 0
        endif

        !--- here the IAU happens ...
        !--- include the IAU tendencies for:

        !- U-wind
        tend%ut (1:ntpts) = tend%ut (1:ntpts) + &
            tend_iau_g(ng)%ut (1:ntpts)*tend_iau_g(ng)%iauwts(1:ntpts)*centFact

        !- V-wind
        tend%vt (1:ntpts) = tend%vt (1:ntpts) + &
            tend_iau_g(ng)%vt (1:ntpts)*tend_iau_g(ng)%iauwts(1:ntpts)*centFact

        !- Exner function perturbation
        tend%pt (1:ntpts) = tend%pt (1:ntpts) + &
            tend_iau_g(ng)%pt (1:ntpts)*tend_iau_g(ng)%iauwts(1:ntpts)*centFact

        !- Total water mixing ratio
        tend%rtt(1:ntpts) = tend%rtt(1:ntpts) + &
            tend_iau_g(ng)%rtt(1:ntpts)*tend_iau_g(ng)%iauwts(1:ntpts)*centFact

        !- Potential temperature
        tend%tht(1:ntpts) = tend%tht(1:ntpts) + &
            tend_iau_g(ng)%tht(1:ntpts)*tend_iau_g(ng)%iauwts(1:ntpts)*centFact


100     format(1x,A46,F10.1,'s',F8.3)

    end subroutine GetIauTendency
    !===========================================================================================

    subroutine readIauTendency(ng, ntpts, m1, m2, m3,ia,iz,ja,jz,time,dtlt)
        implicit none
        integer, intent(in) :: ng,ntpts, m1, m2, m3,ia,iz,ja,jz
        real   , intent(in) :: time,dtlt
        character(len=f_name_length) :: sVarName
        character*1 cgrid

        if (read_iau_tend_from_file == 0) return
        
        read_iau_tend_from_file = 0  ! do it just one time
         
        !--- reading the IAU tendencies from file
        write(cgrid,'(i1)')ng
        call makefnam (sVarName, fileNameIAU, 0, iyear1, imonth1, idate1, itime1*100, 'V', 'g'//cgrid,'bin')
          
        call writeReadIauTendency('read',ng, m1, m2, m3,ia,iz,ja,jz,sVarName,'iauUT' ,tend_iau_g(ng)%ut ,irecUT )
        call writeReadIauTendency('read',ng, m1, m2, m3,ia,iz,ja,jz,sVarName,'iauVT' ,tend_iau_g(ng)%vt ,irecVT )
        call writeReadIauTendency('read',ng, m1, m2, m3,ia,iz,ja,jz,sVarName,'iauPT' ,tend_iau_g(ng)%pt ,irecPT )
        call writeReadIauTendency('read',ng, m1, m2, m3,ia,iz,ja,jz,sVarName,'iauRTT',tend_iau_g(ng)%rtt,irecRTT)
        call writeReadIauTendency('read',ng, m1, m2, m3,ia,iz,ja,jz,sVarName,'iauTHT',tend_iau_g(ng)%tht,irecTHT)
          
        return
          
        !open(2, status='unknown', form='unformatted',file=sVarName(1:len_trim(sVarName)))
        !read(2) tend_iau_g(ng)%ut
        !read(2) tend_iau_g(ng)%vt
        !read(2) tend_iau_g(ng)%pt
        !read(2) tend_iau_g(ng)%rtt
        !read(2) tend_iau_g(ng)%tht
        !close(2)
       
    end subroutine readIauTendency
    !===============================================================================================
        
    subroutine StoreNamelistFileAtIAU(oneNamelistFile)
        implicit none
        type(namelistFile), pointer :: oneNamelistFile

        applyIAU      = oneNamelistFile%applyIAU
        timeWindowIAU = oneNamelistFile%timeWindowIAU
        fileNameIAU   = oneNamelistFile%fileNameIAU
        ramp          = oneNamelistFile%ramp
            
        !--- for IAU with DFI, no ramping is allowed
        if( apply_digital_filter == 1) RAMP = 0.
     
        inquire (iolength=int_byte_size) real_byte_size

    end subroutine StoreNamelistFileAtIAU
    !===============================================================================================
        
    subroutine IAU_VariableWeight(mzp, mxp, myp, nxp, nyp, i0, j0, &
        topt, rtgt, varwts)

        implicit none
        integer, intent(in)  :: mzp
        integer, intent(in)  :: mxp
        integer, intent(in)  :: myp
        integer, intent(in)  :: nxp
        integer, intent(in)  :: nyp
        integer, intent(in)  :: i0
        integer, intent(in)  :: j0
        real,    intent(in)    :: topt(mxp,myp)
        real,    intent(in)    :: rtgt(mxp,myp)
        real,    intent(inout) :: varwts(mzp,mxp,myp)
          
        integer :: i,j,k
        integer :: iGlobal
        integer :: jGlobal
        real :: tnudtopi,tnudlati
        real :: rown,rows,rowe,roww
        real :: zloc,wttop,wtlat,delzi
        character(len=*), parameter :: h="**(IAU_VariableWeight)**"
           
          
        tnudtopi=0.
        if(tnudtop.gt. .01) tnudtopi=1./tnudtop
        tnudlati=0.
        if(tnudlat.gt. .01) tnudlati=1./tnudlat
          
        if(ztop.gt.znudtop) then
            delzi=1./(ztop-znudtop)
        elseif(tnudtop.gt. .01) then
            print*,'Incorrect specification of znudtop ! ! !'
            print*,' znudtop = ',znudtop,' ztop = ',ztop
            stop 'varwt-znud'
        endif
           
        do j=1,myp
            jGlobal = j + j0
            do i=1,mxp
                iGlobal = i + i0

                ! quadratic weight function for lateral boundaries
                rown=max(0.,float(jGlobal+nudlat-nyp))
                rows=max(0.,float(nudlat+1-jGlobal))
                rowe=max(0.,float(iGlobal+nudlat-nxp))
                roww=max(0.,float(nudlat+1-iGlobal))
                wtlat=max(rown*rown,rows*rows,rowe*rowe,roww*roww)  &
                    /float(nudlat*nudlat)
                
                !-- for the IAU wtlat will be 1 in the inner of dommain
                !-- and zero at borders
                wtlat= 1.-wtlat
                
                               !if(j==20) print*,"i,j=",i,j,wtlat
               
                ! linear weight function for top boundary
                do k=1,mzp
                    zloc=zt(k)*rtgt(i,j)+topt(i,j)
                    wttop=max(0.,(zloc-znudtop)*delzi)


                    !-- for the IAU wttop will be 1 from the surface to znutop
                    wttop = 1.- wttop

                    ! full 3-D IAU weight function

                    varwts(k,i,j)= min(wtlat,wttop)
                    varwts(k,i,j)= max(0., min( varwts(k,i,j),1.) )

                !if(j==20.and.i==20)print*,"k=",k,wttop,wtlat,varwts(k,i,j)
                end do
            ! if(j==20) print*,"i=",i,wttop,wtlat,varwts(1,i,j)
            end do
        end do
        !
        !-- do it only one time
        get_iau_weights = 0
           
        return
           
        !
        !--- only for checking if the iau_weights was defined as planned
        print*,"nxp,nyp,nzp=",mxp,myp,mzp
        open(2, status='unknown', form='unformatted', access='direct', &
            recl=mxp*myp, file="IAU_weights.bin")
        do k=1,mzp
            write(2,rec=k) varwts(k,:,:)
        enddo
        close(2)
       ! example of CTL file for mxp=151 ,myp=167 ,mzp=45
       !dset ^IAU_weights.bin
       !undef  -0.9990000E+34
       !title IAU weights
       !xdef  151 linear     1     1
       !ydef  167 linear     1     1
       !zdef   45 linear     1     100
       !tdef   1 linear 21:00z01dec2020            6hr
       !vars   1
       !iauwts     45 99    - RAMS : IAU    [0-1        ]
       !endvars
       !
    end subroutine IAU_VariableWeight
    !=======================================================================================

    subroutine get_dfi_coeffs (DT,timeWindowDF,nsteps,dfi)
        implicit none

        real,   intent(in)  :: DT                ! model time step (sec)
        real,   intent(in)  :: timeWindowDF   ! IAU time scale  (sec)
        integer,intent(in)  :: nsteps         ! number of steps:  timeWindowDF/DT+1
        real,   intent(out) :: dfi(nsteps)

        real pi,arg,wc
        integer n,i,k,nhlf,np1

        !   Calculate DFI coefficients
        !   --------------------------
        pi   = 4.0*atan(1.)
        nhlf = (nsteps+1)/2
           
        do k = 1, nhlf-1
            n   = k-nhlf
            arg = n*pi/nhlf
            wc  = sin(arg)/arg ! Lanczos window
            dfi(k) = wc*sin(n*2.0*pi*DT/timeWindowDF)/(n*pi)
        end do
        dfi(nhlf) = 2*DT/timeWindowDF
        do i = nhlf+1, nsteps
            dfi(i) = dfi(nsteps-i+1)
        end do

        !   Normalize coefficients
        !   ----------------------
        where (dfi < 0.) dfi = 0.0
        dfi = dfi/sum(dfi)
        
        !!!!dfi = dfi/DT ! remember: dynamics multiplies by DT
        
        !--" You will have to look at the code for the subtle details. If done correctly, 
        !-- the applied increment time-averaged over the Corrector Duration must be equal to 
        !--the increment applied as a constant with no digital filter." Larry
        
        dfi = dfi * timeWindowDF/dt
        
        !return
        if(mchnum == master_num) then
            print*,"====================================================="
            print*, '===> DFI initialized for',nsteps,' steps'
            do i =1, nsteps-1
                print '(1x,a,i4,5x,a,g13.6,5x,a,i6)', 'i: ',i,'DFI Coeff: ',dfi(i)
            end do
            print*,"====================================================="
        endif
           
    end subroutine get_dfi_coeffs
    !=======================================================================================

    subroutine writeReadIauTendency(action,ifm,m1, m2, m3,ia,iz,ja,jz,sVarName,varn,localtend,ivar)
        use parlib  , only: parf_bcast ! subroutine
        use readbcst, only: gatherdata ! subroutine
        
        
        implicit none
        character(len=*), intent(in) :: action, sVarName
        character(len=*), intent(in) :: varn
        integer         , intent(in) :: ifm,m1, m2, m3,ia,iz,ja,jz,ivar
        real,  dimension(m1,m2,m3), intent(inout) :: localtend

        !---- local vars
       !real :: globalTend (nnzp(ifm), nnxp(ifm), nnyp(ifm))
        real, allocatable, dimension(:,:,:) :: globalTend 
        integer :: i0, j0, ia1, iz1, ja1, jz1,irec,nrec
        integer :: nunit= 2, pNum, i,j,k,ierr
        integer(kind=i8) :: msgSize
        integer(kind=i8) :: masterZ,masterX,masterY
        logical :: file_exists

 
        !--- allocate the array for the entire model dommain
        !
        allocate(globalTend(nnzp(ifm), nnxp(ifm), nnyp(ifm)), STAT=ierr)
        if (ierr/=0) call fatal_error("Allocating globalTend")
        globalTend = 0.0

        
        !-- initial grid definitions for grid ifm
        i0 = nodei0(mynum,ifm)
        j0 = nodej0(mynum,ifm)
        msgSize=m1*m2*m3
        
        masterZ=nnzp(ifm);  masterX=nnxp(ifm);   masterY=nnyp(ifm)
                

        if(action == 'write') then

            if (mchnum==master_num) then
              !print *,"========================================================================="     
              !print *,0,i0,j0,ia,iz,ja,jz        
              
              !-- Copia a parte do mestre para o global
              call  ex_3_buff(globalTend(:,:,:), localTend(:,:,:), &
                              nnzp(ifm), nnxp(ifm), nnyp(ifm), m1, m2, m3,i0, j0, ia, iz, ja, jz)
              
              do pNum=1,nmachs-1
                !-- Recebe o buffer do escravo pnum
                call parf_get_real(com(pNum)%varBuff,com(pNum)%size, pNUm,com(pNum)%tag)
                
                !-- Pega a posicao dos escravos no global
                i0 = nodei0(pNum+1,ifm)
                j0 = nodej0(pNum+1,ifm)

                !--Copia a parte dos escravos no array global
                call  ex_3_buff(globalTend(:,:,:), com(pNum)%varBuff                       &
                               ,nnzp   (ifm),      nnxp   (ifm),      nnyp   (ifm)         &
                               ,nodemzp(pnum+1,1), nodemxp(pnum+1,1), nodemyp(pnum+1,1)    &
                               ,i0, j0                                                     &
                               ,com(pNum)%tia, com(pNum)%tiz, com(pNum)%tja, com(pNum)%tjz)

               !subroutine ex_3_buff(a,b,n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2)
               !  a(1:n1,i1+i0:i2+i0,j1+j0:j2+j0) = b(1:m1,i1:i2,j1:j2)
               !write(88,*) com(pNum)%tia, com(pNum)%tiz, com(pNum)%tja, com(pNum)%tjz

              enddo
                   
              !print *,nnzp(ifm), nnxp(ifm), nnyp(ifm), m1, m2, m3
              !print *,"========================================================================="
              
              !--- Master opens the file and write the IAU tendencies
              !
              if(ivar==irecUT) open(unit=nunit, status='unknown', file=sVarName(1:len_trim(sVarName)) &
                                   ,form='unformatted', access='direct'                               &
                                   ,recl=int_byte_size*nnxp(ifm)*nnyp(ifm)*nnzp(ifm))
              write(unit=nunit,rec=ivar)  globalTend !- direct access
              
              if(ivar== irecTHT) close(unit=nunit)
            
            
              !--- print some information
              !
              if(ivar==irecUT) then 
                 print*,"====================================================="
                 print*," IAU tendencies for the time period:", iyear1, imonth1, idate1
                 print*," IAU tendencies for the time period:",itime1/100, 'to ', itime1/100+timeWindowIAU/3600.
                 print*," Current time is:",iyear1, imonth1, idate1,itime1/100+time/3600.
              endif
                 write(*,101) trim(varn)," max/min [X/day]:",86400.*maxval(globalTend),86400.*minval(globalTend)
                 call flush(6)
              
              if(ivar== irecTHT) then
                 write(*,100) " Model now should restart back to:", iyear1, imonth1, idate1, itime1/100, " UTC"
                 print*," Addding the IAU tendencies for", timeWindowIAU/3600., 'hours'
                 print*,"====================================================="
                 call flush(6)
               endif
               100     format(1x,A30,4I6,A5)
               101     format(1x,A7,A17,1x,2F13.5)
                          
            else  !(mchnum==master_num)
            
              !---Envia a parte desse escravo para o mestre
              call parf_send_real(localtend ,msgSize, master_num ,mchnum*xFactor)

            endif 
            !do i=1,nnxp(ifm)
            !    write(88,fmt='(I2.2,1X,60(I2,1X))') i,int(globalTend(1,i,:))
            ! enddo
                       
        endif

        
        if(action == 'read') then
            ia1 = nodei0(mynum,ifm) + 1
            iz1 = nodei0(mynum,ifm) + nodemxp(mynum,ifm)
            ja1 = nodej0(mynum,ifm) + 1
            jz1 = nodej0(mynum,ifm) + nodemyp(mynum,ifm)
         
            !if (mchnum==master_num) then
            !    open(unit=nunit, status='unknown', form='unformatted', access='direct', &
            !        recl=4*nnxp(ifm)*nnyp(ifm)*nnzp(ifm), file=sVarName(1:len_trim(sVarName)))
            !    read(unit=nunit,rec=ivar)  globalTend
            !    close(unit=nunit)
            !endif
            
            if (mchnum==master_num) then

 
              inquire(file=sVarName(1:len_trim(sVarName)), EXIST=file_exists)
              if(.not. file_exists) then 
                  print*,"====================================================="
                  print*,"IAU file does not exist:",sVarName(1:len_trim(sVarName))
                  print*,"model will stop"
                  print*,"====================================================="
                  stop "modIau"
              endif               
              
              !--- Master opens the file and read the IAU tendencies
              !
              if(ivar==irecUT) open(unit=nunit, status='old', file=sVarName(1:len_trim(sVarName)) &
                                   ,form='unformatted', access='direct'                               &
                                   ,recl=int_byte_size*nnxp(ifm)*nnyp(ifm)*nnzp(ifm))
              read(unit=nunit,rec=ivar)  globalTend !- direct access
              
              if(ivar== irecTHT) close(unit=nunit)

              !--- print some information
              !
              if(ivar==irecUT) then 
                  print*,"====================================================="
                  print*," Reading the IAU tendencies for the time period:"
                  print*,  iyear1, imonth1, idate1, itime1/100, 'to', nint(itime1/100+timeWindowIAU/3600.), 'UTC'
                  print*," file:", sVarName(1:len_trim(sVarName))
                  print*," Current time is:",iyear1, imonth1, idate1,itime1/100+time/3600.,'UTC'
              endif
                  write(*,101) trim(varn)," max/min [X/day]:",86400.*maxval(globalTend),86400.*minval(globalTend)
                  call flush(6)
              if(ivar== irecTHT) then
                 print*,"====================================================="
                 call flush(6)
               endif

            endif  ! (mchnum==master_num)
                        
            call parf_bcast(globalTend,masterZ,masterX,masterY,master_num)

            call mk_3_buff(globalTend(:,:,:), localTend(:,:,:), &
                           nnzp(ifm), nnxp(ifm), nnyp(ifm), m1, m2, m3,ia1, iz1, ja1, jz1)

        endif

        if(allocated(globalTend)) deallocate(globalTend)

    end subroutine writeReadIauTendency
    !=============================================================================================

    subroutine initComIau(ia,iz,ja,jz)
            !# Initializate IAU comunication
            !#
            !# @note
            !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
            !#
            !# **Brief**: The "master" get information about all slaves over size os comunication and
            !# each tag
            !#
            !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
            !#
            !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
            !#
            !# **Date**: 01 March 2021 (Monday)
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
    
            character(len=*),parameter :: procedureName='**initComIau**' !Name of this procedure
            !
            !Local Parameters

            integer(kind=i8), parameter :: sm=5
            !# Facto to multiplies pNum and get tag
    
            !Input/Output variables
            integer, intent(in) :: ia,iz,ja,jz
    
            !Local variables
            integer :: pNum
            !# slave processor number
            integer :: messSize
            !# Size of message
            integer :: ss(5),rs(5)
    
            !Code
            messSize=nnzp(1)*nnxp(1)*nnyp(1)
            if (mchnum==master_num) then
              if(.not. allocated(com)) allocate(com(nmachs))
              do pNum=1,nmachs-1
                call parf_get_int(rs, sm, pNUm, pNum*xFactor)
                !print *,'Mestre recebendo: ',rs
                com(pNum)%size=rs(1)
                !print *,'Sou o ',mchnum,' recebendo ',com(pNum)%size,' do ',pNum,pNum*xFactor
                com(pNum)%tag=pNum*xFactor
                allocate(com(pnum)%varBuff(rs(1)))
                com(pNum)%tia=rs(2)
                com(pNum)%tiz=rs(3)
                com(pNum)%tja=rs(4)
                com(pNum)%tjz=rs(5)
              enddo      
            else
              !print *, 'Sou o ',mchnum,":",ia,iz,ja,jz
              ss(1)=messSize
              ss(2)=ia
              ss(3)=iz
              ss(4)=ja
              ss(5)=jz
              call parf_send_int(ss, sm, master_num, mchnum*xFactor)
              !print *,'Sou o ',mchnum,' enviando ',messSize,' para o master ',mchnum*xFactor
            endif
    
            !if(mchnum==master_num) then
            !  print *,0,mxp*myp*mzp,ia,iz,ja,jz
            !  do pNum=1,nmachs-1
            !    print *,pNum,com(pNum)%size,com(pNum)%tia,com(pNum)%tiz,com(pNum)%tja,com(pNum)%tjz
            !
            !  enddo
            !endif

            firstTimeIau=.false.

    end subroutine initComIau 


end module modIau
