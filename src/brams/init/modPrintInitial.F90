!=============================================================================================
module modPrintInitial
   !# SOme functions to print the initial information
   !#
   !# @note
   !#
   !# **Documentation/Version Control**: <documents>
   !#
   !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
   !#
   !# **Date**: 02 October 2021 (Saturday)
   !#
   !# **Full Description**: SOme functions to print the initial information
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

   implicit none

   include "constants.f90"
   character(len=*),parameter :: sourceName='modPrintInitial.F90' !Name of this source code
   character(len=*),parameter :: procedureName='**modPrintInitial**' !Name of this procedure

   !Parameters
   !0- leapfrog Robert-Asselin filter, 
   !1- leapfrog Robert-Asselin_Williams filter
   !2- Runge-Kutta 3rd order
   !3- Adams-Bashforth-Moulton 3rd order
   character(len=*), parameter, dimension(4) :: c_dyncore_flag=(/"leapfrog Robert-Asselin filter         " &
                                                  ,"leapfrog Robert-Asselin_Williams filter" &
                                                  ,"Runge-Kutta 3rd order                  " &
                                                  ,"Adams-Bashforth-Moulton 3rd order      "/)

   ! 0 = Forward 2nd order (non-monotonic)
   ! 1 = Walcek monotonic advection
   ! 2 = hybrid: forward 2nd order for thetail,microphysics, TKE and
   character(len=*), parameter,dimension(3) :: c_advmnt=(/'Forward 2nd order (non-monotonic)                      ' &
                                            ,'Walcek monotonic advection                             ' &
                                            ,'hybrid: forward 2nd order for thetail,microphysics, TKE'/)

   !0-none, 2-Mahrer/Pielke, 1-Chen, 3-Harrington
              !    4- CARMA 
   character(len=*), parameter ,dimension(6):: c_RadTYP=(/ 'off          ' &
                                             ,'Mahrer/Pielke' &
                                             ,'Che-Cotton   ' &
                                             ,'Harrington   ' &
                                             ,'UkMet        ' &
                                             ,'rrtmg        '/)

   ! 0- off,
   ! 1- Tremback formulation
   ! 2- Grell-Deveny scheme
   ! 3- Grell-3d formulation
   ! 4- Grell-Deveny scheme FIM/NOAA model
   ! 5- Grell-Freitas scheme 5
   ! 6- Grell-Freitas scheme 6
   ! 7- Grell-Freitas scheme 7
   ! 8- Grell-Freitas scheme 8
   character(len=*), parameter,dimension(9) :: c_nnqparm=(/'off                               ' &
                                           ,  'Tremback formulation              ' &
                                           ,  'Grell-Deveny scheme               ' &
                                           ,  'Grell-3d formulation              ' &
                                           ,  'Grell-Deveny scheme FIM/NOAA model' &
                                           ,  'Grell-Freitas scheme 5            ' &
                                           ,  'Grell-Freitas scheme 6            ' &
                                           ,  'Grell-Freitas scheme 7            ' &
                                           ,  'Grell-Freitas scheme 8            '/)

   !0-off, 1-Souza scheme, 2-Grell-Deveny scheme 3-?
   character(len=*), parameter,dimension(4) :: c_nnshcu=(/'off                ' &
                                             ,'Souza scheme       ' &
                                             ,'Grell-Deveny scheme' &
                                             ,'Nao sei scheme     '/)


   !
   integer,parameter :: maxGrids=5

   interface printOneVarGrid
      module procedure printOneVarGridReal
      module procedure printOneVarGridInt
      module procedure printOneVarGridLogical
      module procedure printOneVarGridchar
   end interface

   interface conv2String
      module procedure convReal2String
      module procedure convInt2String
      module procedure convLogical2String
      module procedure convString2String
   end interface

   type csv 
      !2, 'PODA', 'Partial Oxigen Density in the air index', 'g/m^3'
      integer :: nd 
      character(len=32) :: vname 
      character(len=60) :: varDescription
      character(len=32) :: unit
   end type csv
   type(csv),allocatable :: variables(:)

   integer :: totvars

   private

   public printGridHeader,printOneVarGrid,printGridTail,bramsHeader,printVarHeader,printVarTAil &
         ,conv2String,printOneLineVars,printOneFile,printFileHeader,printFileTail,csvTail,csvHeader &
         ,printOnecsv

   logical :: lastIsBold=.false.

   contains

   !=============================================================================================
   integer function bramsHeader(version,license,nmachs,runtype,namelist,mchnum,master_num,year &
                                ,month,day,hour,experiment,ngrids,timmax,timeunit)
      !# Print a header with some informations
      !#
      !# @note
      !#
      !# **Documentation/Version Control**: <documents>
      !#
      !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@gmail.com>
      !#
      !# **Date**: 02 October 2021 (Saturday)
      !#
      !# **Full Description**: Print a header with some informations
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
   
      implicit none

      character(len=*),parameter :: procedureName='**bramsHeader**' !Name of this procedure
      !
      !Local Parameters
   
      !Input/Output variables
      integer, intent(in) :: nMachs,mchnum,master_num
      integer, intent(in) :: year,month,day,hour,ngrids
      character(len=*),intent(in) :: experiment
      character(len=*), intent(in) :: version,license,runtype,namelist
      real, intent(in) :: timmax
      character(len=*), intent(in) :: timeunit
   
      !Local variables
      character(len=12) :: c0
      integer :: son
      include "modGitInfo.inc"
   
      !Code
          ! master prints initial banner and dumps namelist at stdio
    if (mchnum==master_num) then
      !Make the versioning visible for all
      print *,''
      print *,'----- Git versioning ----'
      print *,lastCommit
      print *,lastMerge
      print *,lastAuthor
      print *,lastGitDate
      print *,'--------------------------'
      print *,''
      write (*,fmt='(A)') c_empty
      write (*,fmt='(A)') c_empty
      write (*,fmt='(A)') c_empty
      write (*,fmt='(A)') c_green
      write (*,fmt='(A)')
      write (*,fmt='(A)')
      write (*,fmt='(A)') "BBBBBBBBBBBBBBBBB   RRRRRRRRRRRRRRRRR                  AAA               MMMMMMMM               MMMMMMMM   SSSSSSSSSSSSSSS "
      write (*,fmt='(A)') "B::::::::::::::::B  R::::::::::::::::R                A:::A              M:::::::M             M:::::::M SS:::::::::::::::S"
      write (*,fmt='(A)') "B::::::BBBBBB:::::B R::::::RRRRRR:::::R              A:::::A             M::::::::M           M::::::::MS:::::SSSSSS::::::S"
      write (*,fmt='(A)') "BB:::::B     B:::::BRR:::::R     R:::::R            A:::::::A            M:::::::::M         M:::::::::MS:::::S     SSSSSSS"
      write (*,fmt='(A)') "  B::::B     B:::::B  R::::R     R:::::R           A:::::::::A           M::::::::::M       M::::::::::MS:::::S            "
      write (*,fmt='(A)') "  B::::B     B:::::B  R::::R     R:::::R          A:::::A:::::A          M:::::::::::M     M:::::::::::MS:::::S            "
      write (*,fmt='(A)') "  B::::BBBBBB:::::B   R::::RRRRRR:::::R          A:::::A A:::::A         M:::::::M::::M   M::::M:::::::M S::::SSSS         "
      write (*,fmt='(A)') "  B:::::::::::::BB    R:::::::::::::RR          A:::::A   A:::::A        M::::::M M::::M M::::M M::::::M  SS::::::SSSSS    "
      write (*,fmt='(A)') "  B::::BBBBBB:::::B   R::::RRRRRR:::::R        A:::::A     A:::::A       M::::::M  M::::M::::M  M::::::M    SSS::::::::SS  "
      write (*,fmt='(A)') "  B::::B     B:::::B  R::::R     R:::::R      A:::::AAAAAAAAA:::::A      M::::::M   M:::::::M   M::::::M       SSSSSS::::S "
      write (*,fmt='(A)') "  B::::B     B:::::B  R::::R     R:::::R     A:::::::::::::::::::::A     M::::::M    M:::::M    M::::::M            S:::::S"
      write (*,fmt='(A)') "  B::::B     B:::::B  R::::R     R:::::R    A:::::AAAAAAAAAAAAA:::::A    M::::::M     MMMMM     M::::::M            S:::::S"
      write (*,fmt='(A)') "BB:::::BBBBBB::::::BRR:::::R     R:::::R   A:::::A             A:::::A   M::::::M               M::::::MSSSSSSS     S:::::S"
      write (*,fmt='(A)') "B:::::::::::::::::B R::::::R     R:::::R  A:::::A               A:::::A  M::::::M               M::::::MS::::::SSSSSS:::::S"
      write (*,fmt='(A)') "B::::::::::::::::B  R::::::R     R:::::R A:::::A                 A:::::A M::::::M               M::::::MS:::::::::::::::SS "
      write (*,fmt='(A)') "BBBBBBBBBBBBBBBBB   RRRRRRRR     RRRRRRRAAAAAAA                   AAAAAAAMMMMMMMM               MMMMMMMM SSSSSSSSSSSSSSS   "
      write (*,fmt='(A)') '--------------------- Brazilian developments on the Regional Atmospheric Modeling System ----------------------------------'
      write (*,fmt='(A)') '                                                '//c_pearlWhite//'Revision '//version//c_noColor
      write (*,fmt='(A)') '                             See more information >> http://brams.cptec.inpe.br'
      write (*,fmt='(A)')
      write (*,fmt='(A)') '                           *** Under license: '//license//' ***'
      write (*,fmt='(A)') c_noColor
      write (*,fmt='(A)')
      son=len(trim(experiment))
      write (*,fmt='(A)') c_pearlWhite//repeat(" ",80-(50-son/2))//trim(experiment)//c_noColor
      write (*,fmt='(A)')
      write (*,fmt='(A,I4,1X,A3,1X,I2.2," - ",I4.4,"h",A,F8.1,1X,A)') c_lightYellow//repeat(" ",25)//'Start at '&
                                                                           ,year,month_name(month),day,hour  &
                                                                           ,c_noColor//' - forecast for '//c_lightYellow &
                                                                           ,timmax,'s'//c_noColor
      write (*,fmt='(A)')
      write (*,fmt='(A)') c_noColor
      write(c0,"(i6.6)") nmachs
      write(*,"(a)") "+-------------------------------------------- Information about submission ---------------------------------------------------"
      write(*,"(a)") "                  Run Type: "//c_lightYellow//trim(runtype)//c_noColor//"            :            Input namelist filename: " &
                          //c_lightYellow//trim(namelist)//c_noColor
      write(*,"(a,I1,a)") "                    ******   Parallel execution using "//c_lightYellow//trim(adjustl(c0))//c_noColor//" processes and " &
                         //c_lightYellow,ngrids,c_noColor//' grids *****'
      write(*,"(a)") "+----------------------------------------------------------------------------------------------------------------------------+"   
      write(*,"(a)") "| "//c_lightYellow//"   More information about submission(non fatal errors, notices, warnings), please, see the file brams.log" &
                      //" and jules.log"//c_noColor//"    |"
      write(*,"(a)") "+----------------------------------------------------------------------------------------------------------------------------+"

   endif
   
      bramsHeader=0

   end function bramsHeader 

   !=============================================================================================
   integer function csvHeader(csvFile)
      !# Print a Header to select variables
      !#
      !# @note
      !#
      !# **Documentation/Version Control**: <documents>
      !#
      !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 05 October 2021 (Tuesday)
      !#
      !# **Full Description**: Print a Header to select variables
      !#
      !# @endnote
      !#
      !# @warning
      !# ![](https://licensebuttons.net/l/by-sa/3.0/88x31.png"")
      !# Now is under CC Attribution-ShareAlike 4.0 International, please see:
      !# &copy; <https://creativecommons.org/licenses/by-sa/4.0/>
      !# @endwarning
      !#
       
      character(len=*),parameter :: procedureName='**csvHeader**' !Name of this procedure
      !
      !Local Parameters
   
      !Input/Output variables
      character(len=*), intent(in) :: csvFile   
   
      !Local variables
      integer :: i,err,varfileUnit,nvar
      type(csv) :: postVar

      !Code
      varfileUnit = 135
      nvar=0
      if(.not. allocated(variables)) allocate(variables(300))
      open(varfileUnit, file = trim(csvFile), status = "old", action = "read")
      do
         read(varfileUnit, *, iostat = err) postVar%nd, postVar%vname, postVar%vardescription, postVar%unit
         if (err /= 0) then
            exit
         end if
         nvar=nvar+1
         variables(nvar)%nd=postVar%nd
         variables(nvar)%vname=postVar%vname
         variables(nvar)%vardescription=postVar%vardescription
         variables(nvar)%unit=postVar%unit
         totvars=nvar
      enddo
      close(varfileUnit)
      write(*,fmt='(A)') "+--------------+--+"//repeat("-",60)//"+"//repeat("-",20)//"+"
      write(*,fmt='("|",A14,"|",A2,"|",A60,"|",A20,"|")') 'Variable','#D','Description','Unit'
      write(*,fmt='(A)') "+--------------+--+"//repeat("-",60)//"+"//repeat("-",20)//"+"

      csvHeader=0
   
   end function csvHeader 

   !=============================================================================================
   integer function csvTail()
      !# Print a tail to select variables
      !#
      !# @note
      !#
      !# **Documentation/Version Control**: <documents>
      !#
      !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 05 October 2021 (Tuesday)
      !#
      !# **Full Description**: Print a tail to select variables
      !#
      !# @endnote
      !#
      !# @warning
      !# ![](https://licensebuttons.net/l/by-sa/3.0/88x31.png"")
      !# Now is under CC Attribution-ShareAlike 4.0 International, please see:
      !# &copy; <https://creativecommons.org/licenses/by-sa/4.0/>
      !# @endwarning
      !#
       
      character(len=*),parameter :: procedureName='**csvTail**' !Name of this procedure
      !
      !Local Parameters
   
      !Input/Output variables
   
      !Local variables
   
      !Code
      write(*,fmt='(A)') "+--------------+--+"//repeat("-",60)//"+"//repeat("-",20)//"+"

      csvTail=0
   
   end function csvTail

   !=============================================================================================
   integer function printOnecsv(varName)
      !# Print one csv line
      !#
      !# @note
      !#
      !# **Documentation/Version Control**: <documents>
      !#
      !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 05 October 2021 (Tuesday)
      !#
      !# **Full Description**: Print one csv line
      !#
      !# @endnote
      !#
      !# @warning
      !# ![](https://licensebuttons.net/l/by-sa/3.0/88x31.png"")
      !# Now is under CC Attribution-ShareAlike 4.0 International, please see:
      !# &copy; <https://creativecommons.org/licenses/by-sa/4.0/>
      !# @endwarning
      !#
       
      character(len=*),parameter :: procedureName='**printOnecsv**' !Name of this procedure
      !
      !Local Parameters
   
      !Input/Output variables
      character(len=*), intent(in) :: varName 

      !Local variables
      type(csv) :: one_post_variable
      integer :: all_size,i
      character(len=14) :: varNameUpper

      !Code
      if(lastIsBold) then 
          write(*,fmt='(A)',advance='no') c_noColor
          lastIsBold=.false.
     else
          write(*,fmt='(A)',advance='no') c_inverted
          lastIsBold=.True.
     endif
      varNameUpper = trim(csvUpperCase(varName))
      one_post_variable = getCsvVarible(varNameUpper)
      if (len(trim(one_post_variable%vName)) .eq. 0) then
         write(*,fmt='(A)',advance='no') c_red
         write(*,fmt='("|",A14,"|",I2.2,"|",A60,"|",A20,"|")',advance='no') varNameUpper,0 &
               ,'Variable not in file','Check the List'

      else
         write(*,fmt='("|",A14,"|",I2.2,"|",A60,"|",A20,"|")',advance='no') varNameUpper,one_post_variable%nd &
               ,one_post_variable%varDescription,one_post_variable%Unit
      endif

      write(*,fmt='(A)') c_noColor
      printOnecsv=0

   
   end function printOnecsv 

   function getCsvVarible(varName) result(one_post_variable)
      character(len = *), intent(in) :: varName
      type(csv) :: one_post_variable
      integer :: i

      do i = 1, totvars
         if(varName .eq. VARIABLES(i)%vname) then
            one_post_variable = variables(i)
            return
         end if
      end do
      one_post_variable%vname = ''
   end function getCsvVarible

   function csvUpperCase(strIn) result(strOut)

      ! converts lower case letters in "strIn" into upper case letters at strOut
      ! all remaining symbols in "strIn" are copied to strOut

      character(len = *), intent(in) :: strIn
      character(len = len(strIn)) :: strOut

      integer :: i
      integer :: charInInt
      integer, parameter :: firstLowCase = iachar("a")
      integer, parameter :: lastLowCase = iachar("z")
      integer, parameter :: firstUpperCase = iachar("A")
      character(len = *), parameter :: h = "**(UpperCase)**"

      strOut = ""
      do i = 1, len_trim(strIn)
         charInInt = iachar(strIn(i : i))
         if (charInInt>=firstLowCase .and. charInInt<=lastLowCase) then
            strOut(i : i) = achar(charInInt - firstLowCase + firstUpperCase)
         else
            strOut(i : i) = strIn(i : i)
         end if
      end do
   end function csvUpperCase



   !=============================================================================================
   integer function printFileHeader()
      !# Print a header to dump files
      !#
      !# @note
      !#
      !# **Documentation/Version Control**: <documents>
      !#
      !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 04 October 2021 (Monday)
      !#
      !# **Full Description**: Print a header to dump files
      !#
      !# @endnote
      !#
      !# @warning
      !# ![](https://licensebuttons.net/l/by-sa/3.0/88x31.png"")
      !# Now is under CC Attribution-ShareAlike 4.0 International, please see:
      !# &copy; <https://creativecommons.org/licenses/by-sa/4.0/>
      !# @endwarning
      !#
       
      character(len=*),parameter :: procedureName='**printFileHeader**' !Name of this procedure
      !
      !Local Parameters
   
      !Input/Output variables
   
      !Local variables
   
      !Code
      write(*,fmt='(A)') "+--------------+---+"//repeat("-",50)//"+"//repeat("-",50)//"+"
      write(*,fmt='("|",A14,"|",A3,"|",A50,"|",A50,"|")') 'Variable','I/O','Folder','File or Prefix'
      write(*,fmt='(A)') "+--------------+---+"//repeat("-",50)//"+"//repeat("-",50)//"+"

      printFileHeader=0
   
   end function printFileHeader 


   !=============================================================================================
   integer function printFileTail()
      !# Print a tail to dump files
      !#
      !# @note
      !#
      !# **Documentation/Version Control**: <documents>
      !#
      !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 04 October 2021 (Monday)
      !#
      !# **Full Description**: Print a tail to dump files
      !#
      !# @endnote
      !#
      !# @warning
      !# ![](https://licensebuttons.net/l/by-sa/3.0/88x31.png"")
      !# Now is under CC Attribution-ShareAlike 4.0 International, please see:
      !# &copy; <https://creativecommons.org/licenses/by-sa/4.0/>
      !# @endwarning
      !#
       
      character(len=*),parameter :: procedureName='**printFileTail**' !Name of this procedure
      !
      !Local Parameters
   
      !Input/Output variables
   
      !Local variables
   
      !Code
      write(*,fmt='(A)') "+--------------+---+"//repeat("-",50)//"+"//repeat("-",50)//"+"

      printFileTail=0
   
   end function printFileTail

   !=============================================================================================
   integer function printOneFile(varNam,compositFileName,io)
      !# Print a line with file information
      !#
      !# @note
      !#
      !# **Documentation/Version Control**: <documents>
      !#
      !# **Author(s)**: Luiz Flávio **e-mail:** <luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 04 October 2021 (Monday)
      !#
      !# **Full Description**: Print a line with file information
      !#
      !# @endnote
      !#
      !# @warning
      !# ![](https://licensebuttons.net/l/by-sa/3.0/88x31.png"")
      !# Now is under CC Attribution-ShareAlike 4.0 International, please see:
      !# &copy; <https://creativecommons.org/licenses/by-sa/4.0/>
      !# @endwarning
      !#
       

      character(len=*),parameter :: procedureName='**printOneFile**' !Name of this procedure
      !
      !Local Parameters
   
      !Input/Output variables
      character(len=*), intent(in) :: varNam,compositFileName,io
   
      !Local variables
      character(len=50) :: fnames(2)
      character(len=14) :: vname
      integer :: lvn
      character(len=3) :: cio
      logical :: dir_exist
  
     !Code
     lvn=len(trim(io))
     if(lvn==1) then
         cio=' '//trim(io)//' '
     else
         cio=trim(io)
     endif


     lvn=len(trim(varNam))
     if(lvn<14) then
         vname=trim(varNam)//repeat(" ",14-lvn)
     else
         vname=varNam(1:14)
     endif


     fnames=getFolder(compositFileName)
     inquire(file=trim(fnames(1))//'/.', exist=dir_exist)   
      if(dir_exist) then
         if(lastIsBold) then 
            write(*,fmt='(A)',advance='no') c_noColor
            lastIsBold=.false.
         else
            write(*,fmt='(A)',advance='no') c_inverted
            lastIsBold=.True.
         endif
     else
         if(lastIsBold) then 
            write(*,fmt='(A)',advance='no') c_noColor//c_blink
            lastIsBold=.false.
         else
            write(*,fmt='(A)',advance='no') c_inverted//c_blink
            lastIsBold=.True.
         endif
     endif         

     write(*,fmt='("|",A14,"|",A3,"|",A50,"|",A50,"|")',advance='no') vname,cio,fnames(1),fnames(2)
     
     printOneFile=0
     write(*,fmt='(A)') c_noColor
   end function printOneFile 

   !=============================================================================================
   function getFolder(compositName) result(fnames)
      !# Get the folder name
      !#
      !# @note
      !#
      !# **Documentation/Version Control**: <documents>
      !#
      !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 04 October 2021 (Monday)
      !#
      !# **Full Description**: Get the folder name
      !#
      !# @endnote
      !#
      !# @warning
      !# ![](https://licensebuttons.net/l/by-sa/3.0/88x31.png"")
      !# Now is under CC Attribution-ShareAlike 4.0 International, please see:
      !# &copy; <https://creativecommons.org/licenses/by-sa/4.0/>
      !# @endwarning
      !#
       

      character(len=*),parameter :: procedureName='**getFolder**' !Name of this procedure
      !
      !Local Parameters
   
      !Input/Output variables
      character(len=*), intent(in) :: compositName
   
      !Local variables
      integer :: slashPos
      character(len=50) :: fnames(2)
   
      !Code
      slashPos=index(compositName,'/',BACK=.true.)
      fnames(1)=compositName(1:slashPos)
      fnames(2)=compositName(slashPos+1:)

   
   end function getFolder 



  !=============================================================================================
  integer function printVarHeader(nvars)
     !# Print a header to dump
     !#
     !# @note
     !#
     !# **Documentation/Version Control**: <documents>
     !#
     !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
     !#
     !# **Date**: 02 October 2021 (Saturday)
     !#
     !# **Full Description**: Print a header to dump
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

     character(len=*),parameter :: procedureName='**printVarHeader**' !Name of this procedure
     !
     !Local Parameters
  
     !Input/Output variables
     integer, intent(in) :: nvars
  
     !Local variables
     integer :: i
     character(len=29) :: gName(nvars)
     character(len=1) :: cn

 
     !Code
     write(cn,fmt='(I1)') nvars

     do i=1,nvars
       write(gName(i),fmt='(A)') "   Variable   :       #      "
     enddo
     write(*,fmt='(A)') '+'//repeat('--------------:--------------+',nvars)
     write(*,fmt='(A,'//cn//'(A,"|"))') "|",gname
     write(*,fmt='(A)') '+'//repeat('--------------:--------------+',nvars)

     printVarHeader=0
  
  end function printVarHeader 

  !=============================================================================================
  integer function printVarTAil(nvars)
     !# Print a header to dump
     !#
     !# @note
     !#
     !# **Documentation/Version Control**: <documents>
     !#
     !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
     !#
     !# **Date**: 02 October 2021 (Saturday)
     !#
     !# **Full Description**: Print a header to dump
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

     character(len=*),parameter :: procedureName='**printVarTail**' !Name of this procedure
     !
     !Local Parameters
  
     !Input/Output variables
     integer, intent(in) :: nvars
  
     !Local variables
     character(len=1) :: cn

 
     !Code
     write(cn,fmt='(I1)') nvars

     write(*,fmt='(A)') '+'//repeat('--------------:--------------+',nvars)

     printVarTail=0
  
  end function printVarTail




  !=============================================================================================
  integer function printGridHeader(ngrids)
     !# Print a header to dump
     !#
     !# @note
     !#
     !# **Documentation/Version Control**: <documents>
     !#
     !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
     !#
     !# **Date**: 02 October 2021 (Saturday)
     !#
     !# **Full Description**: Print a header to dump
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

     character(len=*),parameter :: procedureName='**printGridHeader**' !Name of this procedure
     !
     !Local Parameters
  
     !Input/Output variables
     integer, intent(in) :: ngrids
  
     !Local variables
     integer :: i
     character(len=14) :: gName(ngrids)
     character(len=1) :: cn

 
     !Code
     write(cn,fmt='(I1)') ngrids

     do i=1,ngrids
       write(gName(i),fmt='("    grid ",I1,"    ")') i
     enddo
     write(*,fmt='(A)') '+'//repeat('--------------+',ngrids+1)
     write(*,fmt='("|  variable    |"'//cn//'(A14,"|"))') gName
     write(*,fmt='(A)') '+'//repeat('--------------+',ngrids+1)

     printGridHeader=0
  
  end function printGridHeader 

    !=============================================================================================
  integer function printGridTail(ngrids)
     !# Print a tail to dump
     !#
     !# @note
     !#
     !# **Documentation/Version Control**: <documents>
     !#
     !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
     !#
     !# **Date**: 02 October 2021 (Saturday)
     !#
     !# **Full Description**: Print a tail to dump
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

     character(len=*),parameter :: procedureName='**printGridTail**' !Name of this procedure
     !
     !Local Parameters
  
     !Input/Output variables
     integer, intent(in) :: ngrids
  
     !Local variables
  
     !Code
     write(*,fmt='(A)') '+'//repeat('--------------+',ngrids+1)

     printGridTail=0
  
  end function printGridTail 


  !=============================================================================================
  character(len=14) function convReal2String(value,format)
     !# Convert a real value to String
     !#
     !# @note
     !#
     !# **Documentation/Version Control**: <documents>
     !#
     !# **Author(s)**: Liz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
     !#
     !# **Date**: 04 October 2021 (Monday)
     !#
     !# **Full Description**: Convert a real value to String using format and size=14
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

     character(len=*),parameter :: procedureName='**convReal2String**' !Name of this procedure
     !
     !Local Parameters
  
     !Input/Output variables
     real, intent(in) :: value
     character(len=*),intent(in) :: format
  
     !Local variables
     character(len=14) :: cvar
  
     !Code
     write(cvar,fmt='('//format//')') value

     convReal2String=cvar
  
  end function convReal2String 

  !=============================================================================================
  character(len=14) function convInt2String(value,format)
     !# Convert a real value to String
     !#
     !# @note
     !#
     !# **Documentation/Version Control**: <documents>
     !#
     !# **Author(s)**: Liz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
     !#
     !# **Date**: 04 October 2021 (Monday)
     !#
     !# **Full Description**: Convert a real value to String using format and size=14
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

     character(len=*),parameter :: procedureName='**convReal2String**' !Name of this procedure
     !
     !Local Parameters
  
     !Input/Output variables
     integer, intent(in) :: value
     character(len=*),intent(in) :: format
  
     !Local variables
     character(len=14) :: cvar
  
     !Code
     write(cvar,fmt='('//format//')') value

     convInt2String=cvar
  
  end function convInt2String 

  !=========================================================================
  character(len=14) function convLogical2String(value,format)
     !# Convert a logical value to String
     !#
     !# @note
     !#
     !# **Documentation/Version Control**: <documents>
     !#
     !# **Author(s)**: Liz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
     !#
     !# **Date**: 04 October 2021 (Monday)
     !#
     !# **Full Description**: Convert a logical value to String using format and size=14
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

     character(len=*),parameter :: procedureName='**convReal2String**' !Name of this procedure
     !
     !Local Parameters
  
     !Input/Output variables
     logical, intent(in) :: value
     character(len=*),intent(in) :: format
  
     !Local variables
     character(len=14) :: cvar
     character(len=14) :: varl
  
     !Code
     cvar='     FALSE    '
     if(value) cvar='     TRUE     '

     convLogical2String=cvar
  
  end function convLogical2String 

!=========================================================================
  character(len=14) function convString2String(value,format)
     !# Convert a logical value to String
     !#
     !# @note
     !#
     !# **Documentation/Version Control**: <documents>
     !#
     !# **Author(s)**: Liz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
     !#
     !# **Date**: 04 October 2021 (Monday)
     !#
     !# **Full Description**: Convert a logical value to String using format and size=14
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

     character(len=*),parameter :: procedureName='**convString2String**' !Name of this procedure
     !
     !Local Parameters
  
     !Input/Output variables
     character(len=*), intent(in) :: value
     character(len=*),intent(in) :: format
  
     !Local variables
     character(len=14) :: vname
     integer :: lvn
  
     !Code
     lvn=len(trim(value))
     if(lvn<14) then
         vname=trim(value)//repeat(" ",14-lvn)
     else
         vname=value(1:14)
     endif

     convString2String=vname
  
  end function convString2String 


  !=============================================================================================
  integer function printOneLineVars(nvars,varNam,varVal)
     !# Print one line with four vars
     !#
     !# @note
     !#
     !# **Documentation/Version Control**: <documents>
     !#
     !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
     !#
     !# **Date**: 04 October 2021 (Monday)
     !#
     !# **Full Description**: Print one line with four vars
     !#
     !# @endnote
     !#
     !# @warning
     !# ![](https://licensebuttons.net/l/by-sa/3.0/88x31.png"")
     !# Now is under CC Attribution-ShareAlike 4.0 International, please see:
     !# &copy; <https://creativecommons.org/licenses/by-sa/4.0/>
     !# @endwarning
     !#
      
     character(len=*),parameter :: procedureName='**printOneLineVars**' !Name of this procedure
     !
     !Local Parameters
  
     !Input/Output variables
     integer, intent(in) :: nvars
     character(len=*), intent(in) :: varNam(nvars),varVal(nvars)
  
     !Local variables
     integer :: n


  
     !Code
     if(lastIsBold) then 
          write(*,fmt='(A)',advance='no') c_noColor
          lastIsBold=.false.
     else
          write(*,fmt='(A)',advance='no') c_inverted
          lastIsBold=.True.
     endif
     write(*,fmt='("|")',advance="no")
     do n=1,nvars
         write(*,fmt='(A14,":",A14,"|")',advance="no")  varNam(n),varVal(n)
     enddo 
     write(*,fmt='(A)') c_noColor

     printOneLineVars=0
  
  end function printOneLineVars 


  !=============================================================================================
  integer function printOneVarGridReal(ngrids,varname,var,fmt)
     !# Print one var for each grid
     !#
     !# @note
     !#
     !# **Documentation/Version Control**: <documents>
     !#
     !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
     !#
     !# **Date**: 02 October 2021 (Saturday)
     !#
     !# **Full Description**: Print one var for each grid
     !#
     !# @endnote
     !#
     !# @warning
     !# ![](https://licensebuttons.net/l/by-sa/3.0/88x31.png"")
     !# Now is under CC Attribution-ShareAlike 4.0 International, please see:
     !# &copy; <https://creativecommons.org/licenses/by-sa/4.0/>
     !# @endwarning
     !#
      
     character(len=*),parameter :: procedureName='**printOneVarGridReal**' !Name of this procedure
     !
     !Local Parameters
  
     !Input/Output variables
     integer, intent(in) :: ngrids
     real, intent(in) :: var(ngrids)
     character(len=*),intent(in) :: varname
     character(len=*), intent(in) :: fmt
  
     !Local variables
     character(len=14) :: vname
     integer :: lvn
     character(len=5) :: pc
     character(len=1) :: cn
  
     !Code
     write(cn,fmt='(I1)') ngrids
      if(lastIsBold) then 
          write(*,fmt='(A)',advance='no') c_noColor
          lastIsBold=.false.
      else
          write(*,fmt='(A)',advance='no') c_inverted
          lastIsBold=.True.
      endif

     lvn=len(trim(varname))
     if(lvn<14) then
         vname=trim(varname)//repeat(" ",14-lvn)
     else
         vname=varname(1:14)
     endif
     if(fmt=="") then
         write(*,fmt='("|",A14,"|",'//cn//'(E14.4,"|"),A)',advance='no') vname,var
     else
         write(*,fmt='("|",A14,"|",'//cn//'('//fmt//',"|"),A)',advance='no') vname,var
     endif
     printOneVarGridReal=0
     write(*,fmt='(A)') c_noColor
  
  end function printOneVarGridReal

    !=============================================================================================
  integer function printOneVarGridInt(ngrids,varname,var,fmt)
     !# Print one var for each grid
     !#
     !# @note
     !#
     !# **Documentation/Version Control**: <documents>
     !#
     !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
     !#
     !# **Date**: 02 October 2021 (Saturday)
     !#
     !# **Full Description**: Print one var for each grid
     !#
     !# @endnote
     !#
     !# @warning
     !# ![](https://licensebuttons.net/l/by-sa/3.0/88x31.png"")
     !# Now is under CC Attribution-ShareAlike 4.0 International, please see:
     !# &copy; <https://creativecommons.org/licenses/by-sa/4.0/>
     !# @endwarning
     !#
      
     character(len=*),parameter :: procedureName='**printOneVarGridInt**' !Name of this procedure
     !
     !Local Parameters
  
     !Input/Output variables
     integer, intent(in) :: ngrids
     Integer, intent(in) :: var(ngrids)
     character(len=*),intent(in) :: varname
     character(len=*), intent(in) :: fmt
  
     !Local variables
     character(len=14) :: vname
     integer :: lvn
     character(len=1) :: cn
  
     !Code
      write(cn,fmt='(I1)') ngrids
      if(lastIsBold) then 
          write(*,fmt='(A)',advance='no') c_noColor
          lastIsBold=.false.
      else
          write(*,fmt='(A)',advance='no') c_inverted
          lastIsBold=.True.
      endif
     lvn=len(trim(varname))
     if(lvn<14) then
         vname=trim(varname)//repeat(" ",14-lvn)
     else
         vname=varname(1:14)
     endif
     if(fmt=="") then
         write(*,fmt='("|",A14,"|",'//cn//'(I14,"|"))',advance='no') vname,var
     else
         write(*,fmt='("|",A14,"|",'//cn//'('//fmt//',"|"))',advance='no') vname,var
     endif
     printOneVarGridInt=0
     write(*,fmt='(A)') c_noColor
  
  end function printOneVarGridInt

    !=============================================================================================
  integer function printOneVarGridLogical(ngrids,varname,var)
     !# Print one var for each grid
     !#
     !# @note
     !#
     !# **Documentation/Version Control**: <documents>
     !#
     !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
     !#
     !# **Date**: 02 October 2021 (Saturday)
     !#
     !# **Full Description**: Print one var for each grid
     !#
     !# @endnote
     !#
     !# @warning
     !# ![](https://licensebuttons.net/l/by-sa/3.0/88x31.png"")
     !# Now is under CC Attribution-ShareAlike 4.0 International, please see:
     !# &copy; <https://creativecommons.org/licenses/by-sa/4.0/>
     !# @endwarning
     !#
      
     character(len=*),parameter :: procedureName='**printOneVarGridLogical**' !Name of this procedure
     !
     !Local Parameters
  
     !Input/Output variables
     integer, intent(in) :: ngrids
     logical, intent(in) :: var(ngrids)
     character(len=*),intent(in) :: varname
  
     !Local variables
     character(len=14) :: vname
     integer :: lvn,i
     character(len=14) :: varl(ngrids)
     character(len=1) :: cn

     !Code
      write(cn,fmt='(I1)') ngrids
      if(lastIsBold) then 
          write(*,fmt='(A)',advance='no') c_noColor
          lastIsBold=.false.
      else
          write(*,fmt='(A)',advance='no') c_inverted
          lastIsBold=.True.
      endif
     lvn=len(trim(varname))
     if(lvn<14) then
         vname=trim(varname)//repeat(" ",14-lvn)
     else
         vname=varname(1:14)
     endif
     varl='    FALSE     '
     do i=1,ngrids
         if(var(i)) varl(i)='     TRUE     '
     enddo
     write(*,fmt='("|",A14,"|",'//cn//'(A14,"|"))',advance='no') vname,varl

     printOneVarGridLogical=0
     write(*,fmt='(A)') c_noColor
  
  end function printOneVarGridLogical

    !=============================================================================================
  integer function printOneVarGridChar(ngrids,varname,var)
     !# Print one var for each grid
     !#
     !# @note
     !#
     !# **Documentation/Version Control**: <documents>
     !#
     !# **Author(s)**: Luiz Flávio Rodrigues **e-mail:** <luiz.rodrigues@inpe.br>
     !#
     !# **Date**: 02 October 2021 (Saturday)
     !#
     !# **Full Description**: Print one var for each grid
     !#
     !# @endnote
     !#
     !# @warning
     !# ![](https://licensebuttons.net/l/by-sa/3.0/88x31.png"")
     !# Now is under CC Attribution-ShareAlike 4.0 International, please see:
     !# &copy; <https://creativecommons.org/licenses/by-sa/4.0/>
     !# @endwarning
     !#
      
     character(len=*),parameter :: procedureName='**printOneVarGridChar**' !Name of this procedure
     !
     !Local Parameters
  
     !Input/Output variables
     integer, intent(in) :: ngrids
     character(len=*), intent(in) :: var(ngrids)
     character(len=*),intent(in) :: varname
  
     !Local variables
     character(len=14) :: vname
     integer :: lvn
     character(len=1) :: cn
  
     !Code
     write(cn,fmt='(I1)') ngrids
      if(lastIsBold) then 
          write(*,fmt='(A)',advance='no') c_noColor
          lastIsBold=.false.
      else
          write(*,fmt='(A)',advance='no') c_inverted
          lastIsBold=.True.
      endif
     lvn=len(trim(varname))
     if(lvn<14) then
         vname=trim(varname)//repeat(" ",14-lvn)
     else
         vname=varname(1:14)
     endif
     write(*,fmt='("|",A14,"|",'//cn//'(A14,"|"))',advance='no') vname,var

     printOneVarGridChar=0
     write(*,fmt='(A)') c_noColor
  
  end function printOneVarGridChar


  end module modPrintInitial 