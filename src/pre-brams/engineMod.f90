!=============================================================================================
module engineMod
    !# Module to perform operations with atmos and chem data 
    !#
    !# @note
    !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
    !#
    !# **Brief**: Module to perform operations with atmos and chem data 
    !#
    !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
    !#
    !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
    !#
    !# **Date**: 28 August 2020 (Friday)
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
    character(len=*),parameter :: sourceName='engineMod.f90' !Name of source code
    character(len=*),parameter :: procedureName='**engineMod**' !Name of this procedure
    !
    !Local Parameters

    !Local variables

    contains

    !=============================================================================================
    subroutine concatenate(iTime)
        !# Concatenate data from Atmospheric and chemical
        !#
        !# @note
        !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
        !#
        !# **Brief**: concatenate data from atmospheric data and chemical data
        !#
        !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
        !#
        !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
        !#
        !# **Date**: 28 August 2020 (Friday)
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

        use memoryMod, only: &
            atmosValues, &
            atmosNz, &
            atmosNx, &
            atmosNy, &
            atmosDate, &
            atmosLat, &
            atmosLon, &
            atmosXs, &
            atmosYs, &
            atmosLevels, &
            chemDate, &
            cams

         use chemMod, only: &
            fillChemArrays, &
            chemValues, &
            chemVarNames, &
            nSpc

         use utilsMod, only: &
            interpolationBilinear, &
            sortZ

         use filesMod, only: &
            writeGradsCtlFile
    
        implicit none
    
        include "constants.f90"
        character(len=*),parameter :: procedureName='**concatenate**' !Name of this procedure
        !
        !Local Parameters
    
        !Input/Output variables
        integer, intent(in) :: iTime
        !# Time to be processed
    
        !Local variables
        integer, dimension(4) :: currentDate
        integer :: iMonth,icPos,i,k,iHour,iH
        real :: repr2Values(nSpc,atmosNx, atmosNy, cams%nz)
        real :: repr1Values(nSpc,cams%Nx,cams%Ny,atmosNz)
        real :: llat(2)
        real :: llon(2)
        integer :: zIndex(1),z
    
        !Code
        currentDate=(/atmosDate(iTime)%year,atmosDate(iTime)%month,atmosDate(iTime)%day,atmosDate(iTime)%hour/)
        iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Begin concatenation for data:',currentDate,'I4')

        !Get the current Month  Position in 
        iMonth=atmosDate(iTime)%month
        iHour=atmosDate(iTime)%hour
        iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Imonth, Ihour,nX,nY,nZ',(/iMonth,iHour,atmosnx,atmosny,atmosnz/),'I10')
        k=0
        do i=1,12
            k=k+1
            !print *,k,chemDate(k,1)%month,iMonth
            if(chemDate(k,1)%month==iMonth) then
               icPos=k
               exit
            endif
         enddo
         k=0
         do i=1,24
            k=k+1
            if(chemDate(icPos,k)%hour==iHour) then
               iH=k
               exit
            endif
         enddo

        !Copy correct data from input Cams to temporarily array
        call fillChemArrays(icPos,iH)

        !Vertical Nearest pressure interpolation
        do z =1,atmosNz
            call sortZ(cams%nz,cams%levels,atmosLevels(z),1,zIndex)
            !print *,z,atmosLevels(z),zIndex
            do i=1,nSpc
                repr1Values(i,:,:,z)=chemValues(i,:,:,zIndex(1))             
            enddo
        end do

        !Dealocando chemvalues para proxima interacao
        deallocate(chemValues)


        llat(1)=cams%Lat(1)
        llat(2)=cams%Lat(cams%ny)
        llon(1)=cams%Lon(1)
        llon(2)=cams%Lon(cams%nx)       

        !Horizontal Bilinear Interpolation
        iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Interpolating chem array for atmos grid!')
        do i=1,nSpc
            call interpolationBilinear(repr1Values(i,:,:,:), cams%Nx, cams%Ny, atmosNz, llat, llon,cams%dlon, cams%dlat, &
                                       repr2Values(i,:,:,:), atmosNx, atmosNy, atmosNz, atmosLat, atmosLon, atmosXs, atmosYs)
        enddo

        !call writeGradsCtlFile(iTime,nSpc,repr2Values)
        call writeGradsCtlFile(iTime,nSpc,repr2Values,atmosnx,atmosny,atmosnz,atmosLon(1),atmosLat(1),atmosXs,atmosYs,atmosLevels)
    
    end subroutine concatenate 

end module engineMod 