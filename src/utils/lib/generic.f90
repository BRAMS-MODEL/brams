module genericFunctions
   implicit none

   private

   public outRealSize,arcDistance

contains

   !=============================================================================================
   integer function outRealSize()
       !# get the output byte size accordingly the machine
       !#
       !# @note
       !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
       !#
       !# **Brief**: Return the output real size (1,4,8,etc) accordingly the machine
       !#
       !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
       !#
       !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
       !#
       !# **Date**: 28 July 2020 (Tuesday)
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
       character(len=*),parameter :: sourceName='generic.f90' !Name of source code
       character(len=*),parameter :: procedureName='**getOutputByteSize**' !Name of this procedure
       !
       !Local Parameters
   
       !Input/Output variables
   
       !Local variables
       integer :: output_byte_size
       !# output_byte_size
       real :: bytes_in_float
       !# bytes_in_float
   
       !Code
       inquire(iolength=output_byte_size) bytes_in_float
   
       outRealSize=output_byte_size
   
   end function outRealSize 

   !=============================================================================================
   real function arcDistance(latp1,lonp1,latp2,lonp2)
      !# Return the distance between 2 lat,lons points
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: Return the distance in meters between 2 (lon,lat) points
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 28 July 2020 (Tuesday)
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
      character(len=*),parameter :: sourceName='generic.f90' !Name of source code
      character(len=*),parameter :: procedureName='**arcDistance**' !Name of this procedure
      !
      !Local Parameters
   
      !Input/Output variables
      real,intent(in) :: latP1
      !# Latitude of point 1
      real,intent(in) :: lonP1
      !# Longitude of point 1
      real,intent(in) :: latp2
      !# Latitude of point 2
      real,intent(in) :: lonp2
      !# Longitude of point 2

      !Local variables
      real :: a1
      !# distance in lon direction
      real :: a2
      !# distance in lat direction
   
      !Code
      a1=(lonp1-lonp2)*c_arc
      a2=(latp1-latp2)*c_arc             
      arcDistance=sqrt(a1**2+a2**2)
   
   end function arcDistance 

end module genericFunctions