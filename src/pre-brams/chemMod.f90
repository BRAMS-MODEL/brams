!=============================================================================================
module chemMod
    !# Module for manipulate chem species
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: Module to manipulate chem Species
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 27 August 2020 (Thursday)
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
    character(len=*),parameter :: sourceName='chemMod.f90' !Name of source code
    character(len=*),parameter :: procedureName='**chemMod**' !Name of this procedure
    !
    !Local Parameters

    !Local variables
    character(len=32), allocatable :: chemVarNames(:)
    real, allocatable :: chemValues(:,:,:,:) !nvar,nx,ny,nz
    integer :: nSpc

    contains

    !=============================================================================================
    subroutine createChemEquivalence(onlyRead)
        !# Creates an equivalence table for chemical Species
        !#
        !# @note
        !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
        !#
        !# **Brief**: Creates an equivalence table for chemical species
        !#
        !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
        !#
        !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
        !#
        !# **Date**: 27 August 2020 (Thursday)
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

        use chem1_list, only: &
            nspecies,             &
            spc_alloc,            &
            spc_name,             &
            fdda,                 &
            on

         use utilsMod, only: &
            fileExist, &
            getUnit, &
            releaseUnit

         use memoryMod, only: &
            spcName, &
            spcCamsName, &
            whichCams, &
            factor

        implicit none
    
        include "constants.f90"
        character(len=*),parameter :: procedureName='**createEquivalence**' !Name of this procedure
        !
        !Local Parameters
        character(len=*),parameter :: fileName='./tables/equivalence.dat'
    
        !Input/Output variables
        logical, intent(in) :: onlyRead
    
        !Local variables
        integer :: lunit,i,j,k
        logical :: iFoundVar
    
        !Code

        if(.not. fileExist(fileName)) iErrNumber=dumpMessage(c_tty,c_yes &
            ,sourceName,procedureName,c_fatal,'File '//fileName &
              //' not found. Please, verify and solve it!')
         
         iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Reading species equivalence')
         lunit=getUnit()
         open(unit=lunit, file=trim(fileName), status='old')
         read(lunit,*) nSpc
         allocate(spcName(nSpc))
         allocate(spcCamsName(nSpc,5))
         allocate(whichCams(nSpc,5))
         allocate(factor(nSpc,5))
         whichCams=0
         do i=1,nSpc
            read(lunit,*) spcName(i),spcCamsName(i,1), factor(i,1) &
                                    ,spcCamsName(i,2), factor(i,2) &
                                    ,spcCamsName(i,3), factor(i,3) &
                                    ,spcCamsName(i,4), factor(i,4) &
                                    ,spcCamsName(i,5), factor(i,5)
         enddo
         ierrNumber=releaseUnit(lunit)

         !Verify if all species from chem1_list is present on table
         do i=1,nspecies
            if(spc_alloc(fdda,i)==on)then
               iFoundVar=.false.
               do j=1,nSpc
                  if(trim(spcName(j))==spc_name(i)) then
                     iFoundVar=.true.
                     exit
                  endif
               enddo
               if(.not. iFoundVar) iErrNumber=dumpMessage(c_tty,c_yes &
                  ,sourceName,procedureName,c_fatal,'Var '//trim(spc_Name(i)) &
                  //' from chem1_list not found in equivalence table. Please, verify and solve it!')
            endif
         enddo
         iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'OK!')
        
    end subroutine createChemEquivalence 

    !=============================================================================================
    subroutine fillChemArrays(iMonth,iHour)
        !# Fill chem arrays with correct values
        !#
        !# @note
        !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
        !#
        !# **Brief**: Fill chem arrays with correct values from 1 or 2 cams data
        !#
        !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
        !#
        !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
        !#
        !# **Date**: 27 August 2020 (Thursday)
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
        ! use dump

         use memoryMod, only: &
          spcName, &
          spcCamsName, &
          whichCams, &
          cams, &
          chemDate

         use utilsMod, only: &
            getUnit, &
            releaseUnit, &
            outRealSize

         implicit none
    
         include "constants.f90"
         character(len=*),parameter :: procedureName='**fillChemArrays**' !Name of this procedure
         !
         !Local Parameters
         logical, parameter :: localDump=.true.
    
        !Input/Output variables
        integer,intent(in) :: iMonth
        !# Count Time to be used
        integer,intent(in) :: iHour
        !# Hour to be used
    
        !Local variables
        integer :: i,j,k,iPart,z
        character(len=16) :: partName
        character(len=256) :: localFileName
        character(len=15) :: tDef
        integer :: lunit,recordLen,iRec
    
        ! !Code
         allocate(chemValues(nSpc,cams%nx,cams%ny,cams%nz))
         chemValues=0.0
         !Fill species name
         iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Filling chem values for time',(/imonth,iHour/),'I4.4')
!print *,'LFR->DGB 0',imonth,ihour
          do i=1,nSpc
             !chemVarNames(i)=spcName(i)
             iPart=0
             do j=1,5
                if(trim(spcCamsName(i,j))/='NONE') iPart=iPart+1
             enddo
             !Walking only in valid parts
             do j=1,iPart
                partName=trim(spcCamsName(i,j))
                do k=1,cams%nv
                    if(trim(partName)==trim(cams%VarName(k))) then
                         chemValues(i,:,:,:)=chemValues(i,:,:,:)+cams%Values(iMonth,iHour,k,:,:,:)*1.0e+9
                    endif
                enddo
             enddo
          enddo
!print *,'LFR->DGB 1'

        if(.not. localdump) return

        write(tDef,fmt='(I2.2,":00z",I2.2,A3,I4.4)')  chemDate(iMonth,iHour)%hour,1 &
                                                      ,month_Name(chemDate(iMonth,iHour)%month),2019
        write(localFileName,fmt='("Chem-",I2.2,"-",I2.2)') iMonth,iHour
        lunit=getUnit()
        open(unit=lunit,file=trim(localFileName)//'.ctl',action='write',status='replace')
        !writing the name of grads file
        write(lunit,*) 'dset ^'//trim(localFileName)//'.gra'
        !writing others infos to ctl
        write(lunit,*) 'undef -0.9990000E+34'
        write(lunit,*) 'title '//procedureName
        write(lunit,*) 'xdef ',cams%Nx,' linear ',cams%Lon(1),cams%dlon
        write(lunit,*) 'ydef ',cams%Ny,' linear ',cams%Lat(1),cams%dlat
        write(lunit,*) 'zdef ',cams%Nz,'levels',cams%Levels
        write(lunit,*) 'tdef 1 linear '//tDef//' 1mo'
        write(lunit,*) 'vars ',nSpc
        do i=1,nSpc
            write(lunit,*) SpcName(i),cams%Nz,'99 ',spcName(i)
        enddo
        write(lunit,*) 'endvars'
        ierrNumber=releaseUnit(lunit)

        recordLen=outRealSize()*cams%Nx*cams%Ny
        lunit=getUnit()
        open(unit=lunit,file=trim(localFileName)//'.gra',action='WRITE',status='REPLACE' &
             ,form='UNFORMATTED',access='DIRECT',recl=recordLen)

        irec=1
        do i=1,nSpc
            do z=1,cams%Nz
                write(lunit,rec=irec) chemValues(i,:,:,z)
                irec=irec+1
            enddo
        enddo
        ierrNumber=releaseUnit(lunit)

    end subroutine fillChemArrays

end module chemMod 