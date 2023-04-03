Program geraCI
  use mpi
  use modMemory
  use modUtils, only: readNamelist,init,encerra
  use modAnalysis, only: analysisNcep,analysisNasa
  use ModDateUtils, only: date_add_to

  character(len=*), parameter :: cfn1='("'
  character(len=*), parameter :: cfn2='",I4.4,I2.2,' &
              //'I2.2,"_",I2.2,"+",I4.4,I2.2,I2.2,"_",I4.4,".V01.nc4")'
  integer :: i
  character(len=256) :: innpr
  integer :: iyy,imm,idd,ihh

  call readNamelist()

  call init(tIncrem,lastTime)

  if(source=='NCEP') then

    do i=ini(id),fim(id),tIncrem
       write(innpr,fmt='("'//trim(prefix)//'",I3.3)') i
       call analysisNcep(innpr,id,outFolder)
    enddo

  elseif(source=='NASA') then

    do i=ini(id),fim(id),tIncrem
      call date_add_to(iyear1,imonth1,idate1,itime1,real(i),'h' &
        ,iyy,imm,idd,ihh)
        ihh=ihh/100
      write(innpr,fmt=cfn1//trim(prefix)//cfn2) iYear1,imonth1,idate1 &
            ,itime1,iyy,imm,idd,ihh
      print *,trim(innpr)
      call analysisNasa(innpr,id,outFolder)

    enddo 
  endif 

  i=encerra()

end program geraCI



