!=============================================================================================
program filterOut
   !# Program to filter the brams output
   !#
   !# @note
   !#
   !# **Documentation/Version Control**: <documents>
   !#
   !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
   !#
   !# **Date**: 05 October 2021 (Tuesday)
   !#
   !# **Full Description**: Program to filter the brams output
   !#
   !# @endnote
   !#
   !# @warning
   !# ![](https://licensebuttons.net/l/by-sa/3.0/88x31.png"")
   !# Now is under CC Attribution-ShareAlike 4.0 International, please see:
   !# &copy; <https://creativecommons.org/licenses/by-sa/4.0/>
   !# @endwarning
   !#
    
   !Use area
   use dump !Dump contains a lot of functions for debugs and formated printouts
   use utilsMod, only: fileExist,cutLine

   implicit none

   include "constants.f90"
   character(len=*),parameter :: sourceName='filterOut.f90' !Name of this source code
   character(len=*),parameter :: procedureName='**filterOut**' !Name of this procedure
   !
   !Local Parameters

   !Local variables
   character(len=256) :: num1char
   character(len=256) :: linha,linhap
   integer :: err, result, sizel, i

   !Code
   if(command_argument_count().ne.1)then
      iErrNumber=dumpMessage(c_tty,c_yes,'','',c_fatal,'The name of file to convert must be passed!')
   endif
   call get_command_argument(1,num1char)   !first, read in the file name
   if(.not. fileExist(num1char)) then
      iErrNumber=dumpMessage(c_tty,c_yes,'','',c_fatal,'File not exist. Please check it!')
   endif
   open(unit=22,file=num1char,status='old',action='read')
   open(unit=55,file='output.txt', status='replace',action='write')
   do
      read(unit=22, iostat=err, fmt='(A)') linha
      if(err/=0) exit
      do i=1,5
         linhap=cutLine(linha,c_darkGrey   )
         linha=cutLine(linhap,c_peach      )
         linhap=cutLine(linha,c_lightGreen )
         linha=cutLine(linhap,c_lightYellow)
         linhap=cutLine(linha,c_lightBlue  )
         linha=cutLine(linhap,c_pink       )
         linhap=cutLine(linha,c_lightAqua  )
         linha=cutLine(linhap,c_pearlWhite )
         linhap=cutLine(linha,c_black      )
         linha=cutLine(linhap,c_red        )
         linhap=cutLine(linha,c_green      )
         linha=cutLine(linhap,c_yellow     )
         linhap=cutLine(linha,c_blue       )
         linha=cutLine(linhap,c_purple     )
         linhap=cutLine(linha,c_aqua       )
         linha=cutLine(linhap,c_blink      )
         linhap=cutLine(linha,c_noColor    )
         linha=cutLine(linhap,c_Iblack     )
         linhap=cutLine(linha,c_Ired       )
         linha=cutLine(linhap,c_Igreen     )
         linhap=cutLine(linha,c_Iyellow    )
         linha=cutLine(linhap,c_Iblue      )
         linhap=cutLine(linha,c_IMagenta   )
         linha=cutLine(linhap,c_Icyan      )
         linhap=cutLine(linha,c_underline  )
         linha=cutLine(linhap,c_strike     )
         linhap=cutLine(linha,c_bold       )
         linha=cutLine(linhap,c_normal     )
         linhap=cutLine(linha,c_blinking   )
         linha=cutLine(linhap,c_reverse    )
         linhap=cutLine(linha,c_inverted   )
         linha=cutLine(linhap,c_bg_black   )
         linhap=cutLine(linha,c_bg_red     )
         linha=cutLine(linhap,c_bg_green   )
         linhap=cutLine(linha,c_bg_brown   )
         linha=cutLine(linhap,c_bg_blue    )
         linhap=cutLine(linha,c_bg_purple  )
         linha=cutLine(linhap,c_bg_cyan    )
         linhap=cutLine(linha,c_bg_lgray   )
         if (i==5) write(55,fmt='(A)') trim(linhap)
         linha=linhap
      enddo
   enddo
   close(unit=22)
   close(unit=55)
   open(unit=22,file='output.txt',status='old',action='read')
   open(unit=55,file=num1char, status='replace',action='write')
   do
      read(unit=22, iostat=err, fmt='(A)') linha
      if(err/=0) exit   
      write(55,fmt='(A)') trim(linha)
   enddo 
   close(unit=22)
   close(unit=55)   
end program filterOut 

  