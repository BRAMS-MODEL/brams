module MPI_IO_ENGINE

  implicit none

  private
  
  public :: OPEN_FILE_READ
  public :: OPEN_FILE_WRITE
  public :: CLOSE_FILE_READ
  public :: CLOSE_FILE_WRITE
  public :: READ_INITIAL_INFO
  public :: INIT_FILETYPES
  public :: READ_BRAMS_DATA
  public :: WRITE_BRAMS_DATA


  include 'mpif.h'

  INTEGER(KIND=MPI_OFFSET_KIND) :: DISPR, DISPW
  INTEGER                       :: H_INPUT, H_OUTPUT

  interface WRITE_BRAMS_DATA
     module procedure WRITE_BRAMS_DATA_2D
     module procedure WRITE_BRAMS_DATA_2D_I
     module procedure WRITE_BRAMS_DATA_3D
     module procedure WRITE_BRAMS_DATA_4D
  end interface

  interface READ_BRAMS_DATA
     module procedure READ_BRAMS_DATA_2D
     module procedure READ_BRAMS_DATA_2D_I
     module procedure READ_BRAMS_DATA_3D
     module procedure READ_BRAMS_DATA_4D
  end interface

CONTAINS

!!$$$$$$$$$$$$$$$$$$ INICIO ABRE ARQUIVO LEITURA $$$$$$$$$$$$$$$$$$
  SUBROUTINE OPEN_FILE_READ(fileName)

    implicit none

    CHARACTER*(*), intent(in) :: fileName

    INTEGER :: ierr

    DISPR = 0

    call MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_RDONLY, &
         MPI_INFO_NULL, H_INPUT, ierr);

    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro ao abrir arquivo.'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

  END SUBROUTINE OPEN_FILE_READ
!!$$$$$$$$$$$$$$$$$$ FIM ABRE ARQUIVO LEITURA $$$$$$$$$$$$$$$$$$$$$

!!$$$$$$$$$$$$$$$$$$ INICIO ABRE ARQUIVO ESCRITA $$$$$$$$$$$$$$$$$$
  SUBROUTINE OPEN_FILE_WRITE(fileName)

    implicit none

    CHARACTER*(*), intent(in) :: fileName

    INTEGER :: ierr

    DISPW = 0

    call MPI_File_open(MPI_COMM_WORLD, fileName, &
         MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, H_OUTPUT, ierr);

    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro ao abrir arquivo.'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

  END SUBROUTINE OPEN_FILE_WRITE
!!$$$$$$$$$$$$$$$$$$ FIM ABRE ARQUIVO ESCRITA $$$$$$$$$$$$$$$$$$$$$

!!$$$$$$$$$$$$$$$$$$ INICIO FECHA ARQUIVO LEITURA $$$$$$$$$$$$$$$$$$
  SUBROUTINE CLOSE_FILE_READ()

    implicit none

    INTEGER :: ierr

    call MPI_File_close(H_input, ierr);

    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro fechando arquivo.'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

  END SUBROUTINE CLOSE_FILE_READ
!!$$$$$$$$$$$$$$$$$$ FIM FECHA ARQUIVO LEITURA $$$$$$$$$$$$$$$$$$$$$

!!$$$$$$$$$$$$$$$$$$ INICIO FECHA ARQUIVO ESCRITA $$$$$$$$$$$$$$$$$$
  SUBROUTINE CLOSE_FILE_WRITE()

    implicit none

    INTEGER :: ierr

    call MPI_File_close(H_output, ierr);

    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro fechando arquivo de escrita.'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

  END SUBROUTINE CLOSE_FILE_WRITE
!!$$$$$$$$$$$$$$$$$$ FIM FECHA ARQUIVO ESCRITA $$$$$$$$$$$$$$$$$$$$$

!!$$$$$$$$$$$$$$$$$$ INICIO LEITURA BRAMS $$$$$$$$$$$$$$$$$$
  SUBROUTINE READ_INITIAL_INFO(rank, NNXP, NNYP, NNZP, NN4P, NN5P, &
       xbeg, xend, ybeg, yend, gZoneSize, fileName)

    implicit none

    INTEGER, intent(in)       :: rank
    INTEGER, intent(out)      :: nnxp, nnyp, nnzp, nn4p, nn5p, &
         xbeg, xend, ybeg, yend, gZoneSize
    CHARACTER*(*), intent(in) :: fileName

    INTEGER :: cSize, ierr, s1, s2
    INTEGER, POINTER :: idata(:,:)
    INTEGER :: iheader(6)
    real, pointer :: w(:,:)
    integer, pointer :: ixbeg(:), ixend(:), iybeg(:), iyend(:)

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, cSize, ierr)

    ALLOCATE(idata(4, cSize))

    OPEN(1, FILE = fileName)

    READ(1,*) iheader

    !READ(1,*) idata

    CLOSE(1)

    !print *,'iheader',iheader,'fimiheader'

    !print *,'idata',idata,'fimidata'

    NNXP = iheader(1)
    NNYP = iheader(2)
    NNZP = iheader(3)
    NN4P = iheader(4)
    NN5P = iheader(5)

    gZoneSize = iheader(6)

    !print *, rank, nnxp, nnyp, nnzp, gZoneSize

    allocate(w(nnxp,nnyp))
    allocate(ixbeg(cSize))
    allocate(ixend(cSize))
    allocate(iybeg(cSize))
    allocate(iyend(cSize))

    do s1=1,nnxp 
       do s2=1,nnyp 
          w(s1,s2) = 1
       enddo
    enddo

    call decomp_par(nnxp, nnyp, cSize, w, ixbeg, ixend, iybeg, iyend)

    do s1=1,cSize
       idata(1, s1) = ixbeg(s1)
       idata(2, s1) = ixend(s1)
       idata(3, s1) = iybeg(s1)
       idata(4, s1) = iyend(s1)
    enddo

    !! Soma ao rank
    !! 1 - Porque rank começa no indice 0
    xbeg = idata(1, rank + 1)
    xend = idata(2, rank + 1)
    ybeg = idata(3, rank + 1)
    yend = idata(4, rank + 1)

    !print *, rank, xbeg, xend, ybeg, yend

  END SUBROUTINE READ_INITIAL_INFO
!!$$$$$$$$$$$$$$$$$$ FIM LEITURA BRAMS $$$$$$$$$$$$$$$$$$$$$


!!$$$$$$$$$$$$$$$$$$ INICIO DEFINIÇÃO ESTRUTRAS UTILIZADAS $$$$$$$$$$$$$$$$$$
  SUBROUTINE INIT_FILETYPES(nnzp, nnxp, nnyp, nn4p, nn5p, &
       xbeg, xend, ybeg, yend, gZoneSize, gSizes, &
       subSizesRead, subSizesWrite, startRead, startWrite, startWriteGlobal, &
       fileTypeRead, fileTypeWrite, fileTypeWriteGlobal)

    implicit none

    integer, intent(in)  :: nnzp, nnxp, nnyp, nn4p, nn5p, xbeg, xend, ybeg, yend, gZoneSize
    integer, intent(out) :: gSizes(5), subSizesRead(5), subSizesWrite(5), &
         startRead(5), startWrite(5), startWriteGlobal(5)
    INTEGER, intent(out) :: fileTypeRead, fileTypeWrite, fileTypeWriteGlobal

    INTEGER :: ierr

    gSizes(1) = NNZP
    gSizes(2) = NNXP
    gSizes(3) = NNYP
    gSizes(4) = NN4P
    gSizes(5) = NN5P

!!$ Dimensões 1, 4 e 5 não possuem divisões de dominio e nem gosth zone
    startRead(1) = 0
    subSizesRead(1) = gSizes(1)

    startRead(4) = 0
    subSizesRead(4) = gSizes(4)

    startRead(5) = 0
    subSizesRead(5) = gSizes(5)

!!$ Subtrai um porque os indices do BRAMS começam em 1 e aqui eles começam em 0
!!$ Subtrai o gZone porque o BRAMS não o leva em consideração ao passar as coordenadas
    startRead(2) = xbeg -1 - gZoneSize
    startRead(3) = ybeg -1 - gZoneSize

!!$ Tamanho de leitura é o tamanho definido pelo BRAMS mais o GZone dos dois extremos
    subSizesRead(2) = (xend - xbeg + 1) + (2 * gZoneSize)
    subSizesRead(3) = (yend - ybeg + 1) + (2 * gZoneSize)

!!$ Para processos que não estão nas bordas
    subSizesWrite(1) = subSizesRead(1)
    subSizesWrite(2) = subSizesRead(2) - (2 * gZoneSize)
    subSizesWrite(3) = subSizesRead(3) - (2 * gZoneSize)
    subSizesWrite(4) = subSizesRead(4)
    subSizesWrite(5) = subSizesRead(5)

    startWrite(1) = 0
    startWrite(2) = gZoneSize
    startWrite(3) = gZoneSize
    startWrite(4) = 0
    startWrite(5) = 0

!!$ Se for um processo da borda superior
    if (startRead(2) == 0) then
       startWrite(2) = 0
       subSizesWrite(2) = subSizesWrite(2) + gZoneSize
    end if
!!$ Se for um processo da borda esquerda
    if (startRead(3) == 0) then
       startWrite(3) = 0
       subSizesWrite(3) = subSizesWrite(3) + gZoneSize
    end if
!!$ Se for um processo da borda inferior
    if ((gSizes(2) - subSizesRead(2)) == startRead(2)) then
       subSizesWrite(2) = subSizesWrite(2) + gZoneSize
    end if
!!$ Se for um processo da borda direita
    if ((gSizes(3) - subSizesRead(3)) == startRead(3)) then
       subSizesWrite(3) = subSizesWrite(3) + gZoneSize
    end if

!!$ Estrutura global do arquivo para leitura

    call MPI_Type_create_subarray(5, gSizes, subSizesRead, startRead, MPI_ORDER_FORTRAN, MPI_REAL, fileTypeRead, ierr);

    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro na criação da estrutura do arquivo leitura.'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

    call MPI_Type_commit(fileTypeRead, ierr)

    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro commit leitura'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

!!$ Estrutura global do arquivo para escrita

    startWriteGlobal(1) = startRead(1)
    startWriteGlobal(2) = startRead(2)
    startWriteGlobal(3) = startRead(3)
    startWriteGlobal(4) = startRead(4)
    startWriteGlobal(5) = startRead(5)

    if (startWriteGlobal(2) > 0) then
       startWriteGlobal(2) = startWriteGlobal(2) + 1
    end if
    if (startWriteGlobal(3) > 0) then
       startWriteGlobal(3) = startWriteGlobal(3) + 1
    end if

    call MPI_Type_create_subarray(5, gSizes, subSizesWrite, startWriteGlobal, &
         &MPI_ORDER_FORTRAN, MPI_REAL, fileTypeWriteGlobal, ierr);

    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro na criação da estrutura do arquivo escrita local.'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

    call MPI_Type_commit(fileTypeWriteGlobal, ierr)

    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro commit escrita local'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

!!$ Estrutura do processo para escrita no arquivo

    call MPI_Type_create_subarray(5, subSizesRead, subSizesWrite, startWrite, MPI_ORDER_FORTRAN, MPI_REAL, fileTypeWrite, ierr);

    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro na criação da estrutura do arquivo escrita.'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

    call MPI_Type_commit(fileTypeWrite, ierr)

    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro commit escrita global'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

  END SUBROUTINE INIT_FILETYPES
!!$$$$$$$$$$$$$$$$$$ FIM DEFINIÇÃO ESTRUTRAS UTILIZADAS $$$$$$$$$$$$$$$$$$


!!$$$$$$$$$$$$$$$$$$ INICIO DA LEITURA DO ARQUIVO DE ENTRADA $$$$$$$$$$$$$$$$$$
!!$  SUBROUTINE READ_BRAMS_DATA(nints, buf, fileTypeRead)
!!$
!!$    implicit none
!!$
!!$    INTEGER, intent(in)        :: nints
!!$    REAL, POINTER, intent(out) :: buf(:)
!!$    INTEGER, intent(in)        :: fileTypeRead
!!$    
!!$    INTEGER :: ierr, structureSize
!!$
!!$    call MPI_File_set_view(H_INPUT, DISPR, MPI_REAL, fileTypeRead, &
!!$         "native", MPI_INFO_NULL, ierr);
!!$    if (ierr .ne. MPI_SUCCESS) then
!!$       print *,'Erro set view 1.'
!!$       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
!!$    end if
!!$
!!$    call MPI_File_read_all(H_INPUT, buf, nints, MPI_REAL, &
!!$         MPI_STATUS_IGNORE, ierr);
!!$
!!$    if (ierr .ne. MPI_SUCCESS) then
!!$       print *,'Erro de leitura do arquivo.'
!!$       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
!!$    end if
!!$
!!$    call MPI_Type_extent(fileTypeRead, structureSize, ierr)
!!$
!!$    DISPR = DISPR + structureSize
!!$
!!$  END SUBROUTINE READ_BRAMS_DATA


  SUBROUTINE BRAMS_SET_VIEW_REAL_INPUT(fileTypeRead)

    implicit none

    integer, intent(in) :: fileTypeRead

    integer :: ierr

    call MPI_File_set_view(H_INPUT, DISPR, MPI_REAL, fileTypeRead, &
         "native", MPI_INFO_NULL, ierr);
    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro BRAMS_SET_VIEW_REAL_INPUT'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

  END SUBROUTINE BRAMS_SET_VIEW_REAL_INPUT



  SUBROUTINE BRAMS_SET_VIEW_INTEGER_INPUT(fileTypeRead)

    implicit none

    integer, intent(in) :: fileTypeRead

    integer :: ierr

    call MPI_File_set_view(H_INPUT, DISPR, MPI_INTEGER, fileTypeRead, &
         "native", MPI_INFO_NULL, ierr);
    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro BRAMS_SET_VIEW_REAL_INPUT'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

  END SUBROUTINE BRAMS_SET_VIEW_INTEGER_INPUT



  SUBROUTINE UPDATE_DISPR(fileTypeRead)

    implicit none

    integer, intent(in) :: fileTypeRead

    INTEGER :: ierr, structureSize

    call MPI_Type_extent(fileTypeRead, structureSize, ierr)

    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro UPDATE_DISPR'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

    DISPR = DISPR + structureSize

  END SUBROUTINE UPDATE_DISPR


  SUBROUTINE READ_BRAMS_DATA_2D(nints, bufr2d, fileTypeRead)

    implicit none

    INTEGER, intent(in) :: nints
    REAL, intent(out)   :: bufr2d(:,:)
    INTEGER, intent(in) :: fileTypeRead
    
    INTEGER :: ierr, structureSize

    call BRAMS_SET_VIEW_REAL_INPUT(fileTypeRead)

    call MPI_File_read_all(H_INPUT, bufr2d(1,1), nints, MPI_REAL, &
         MPI_STATUS_IGNORE, ierr);

    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro de leitura do arquivo.'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

    call UPDATE_DISPR(fileTypeRead)

  END SUBROUTINE READ_BRAMS_DATA_2D


  SUBROUTINE READ_BRAMS_DATA_2D_I(nints, bufi2d, fileTypeRead)

    implicit none

    INTEGER, intent(in)  :: nints
    INTEGER, intent(out) :: bufi2d(:,:)
    INTEGER, intent(in)  :: fileTypeRead
    
    INTEGER :: ierr, structureSize

    call BRAMS_SET_VIEW_INTEGER_INPUT(fileTypeRead)

    call MPI_File_read_all(H_INPUT, bufi2d(1,1), nints, MPI_INTEGER, &
         MPI_STATUS_IGNORE, ierr);

    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro de leitura do arquivo.'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

    call UPDATE_DISPR(fileTypeRead)

  END SUBROUTINE READ_BRAMS_DATA_2D_I


  SUBROUTINE READ_BRAMS_DATA_3D(nints, bufr3d, fileTypeRead)

    implicit none

    INTEGER, intent(in) :: nints
    REAL, intent(out)   :: bufr3d(:,:,:)
    INTEGER, intent(in) :: fileTypeRead
    
    INTEGER :: ierr, structureSize

    call BRAMS_SET_VIEW_REAL_INPUT(fileTypeRead)

    call MPI_File_read_all(H_INPUT, bufr3d(1,1,1), nints, MPI_REAL, &
         MPI_STATUS_IGNORE, ierr);

    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro de leitura do arquivo.'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

    call UPDATE_DISPR(fileTypeRead)

  END SUBROUTINE READ_BRAMS_DATA_3D


  SUBROUTINE READ_BRAMS_DATA_4D(nints, bufr4d, fileTypeRead)

    implicit none

    INTEGER, intent(in) :: nints
    REAL, intent(out)   :: bufr4d(:,:,:,:)
    INTEGER, intent(in) :: fileTypeRead
    
    INTEGER :: ierr, structureSize

    call BRAMS_SET_VIEW_REAL_INPUT(fileTypeRead)

    call MPI_File_read_all(H_INPUT, bufr4d(1,1,1,1), nints, MPI_REAL, &
         MPI_STATUS_IGNORE, ierr);

    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro de leitura do arquivo.'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

    call UPDATE_DISPR(fileTypeRead)

  END SUBROUTINE READ_BRAMS_DATA_4D

!!$$$$$$$$$$$$$$$$$$ FIM DA LEITURA DO ARQUIVO DE ENTRADA $$$$$$$$$$$$$$$$$$


!!$$$$$$$$$$$$$$$$$$ INICIO DA ESCRITA DO ARQUIVO DE SAIDA $$$$$$$$$$$$$$$$$$
!!$  SUBROUTINE WRITE_BRAMS_DATA(buf, fileTypeWrite, fileTypeWriteGlobal)
!!$
!!$    integer, intent(in) :: fileTypeWrite, fileTypeWriteGlobal
!!$    real, intent(in)    :: buf(:)
!!$
!!$    INTEGER :: nints, ierr, structureSize
!!$
!!$    call MPI_File_set_view(H_OUTPUT, DISPW, MPI_REAL, fileTypeWriteGlobal, "native", MPI_INFO_NULL, ierr);
!!$    if (ierr .ne. MPI_SUCCESS) then
!!$       print *,'Erro set_view output.'
!!$       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
!!$    end if
!!$
!!$    call MPI_File_write_all(H_OUTPUT, buf(1), 1, fileTypeWrite, MPI_STATUS_IGNORE, ierr);
!!$
!!$    if (ierr .ne. MPI_SUCCESS) then
!!$       print *,'Erro escrita.'
!!$       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
!!$    end if
!!$
!!$    call MPI_Type_extent(fileTypeWriteGlobal, structureSize, ierr)
!!$
!!$    DISPW = DISPW + structureSize
!!$
!!$  END SUBROUTINE WRITE_BRAMS_DATA
!!$
!!$  !---
!!$
!!$  SUBROUTINE WRITE_BRAMS_DATA_NEW(bufr2d, bufi2d, bufr3d, bufr4d, &
!!$       fileTypeWrite, fileTypeWriteGlobal)
!!$
!!$    integer, intent(in) :: fileTypeWrite, fileTypeWriteGlobal
!!$    real, optional, intent(in)    :: bufr2d(:,:)
!!$    integer, optional, intent(in) :: bufi2d(:,:)
!!$    real, optional, intent(in)    :: bufr3d(:,:,:)
!!$    real, optional, intent(in)    :: bufr4d(:,:,:,:)
!!$
!!$    INTEGER :: nints, ierr, structureSize
!!$
!!$    call MPI_File_set_view(H_OUTPUT, DISPW, MPI_REAL, fileTypeWriteGlobal, "native", MPI_INFO_NULL, ierr);
!!$    if (ierr .ne. MPI_SUCCESS) then
!!$       print *,'Erro set_view output.'
!!$       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
!!$    end if
!!$
!!$    if (present(bufr2d)) then
!!$       call MPI_File_write_all(H_OUTPUT, bufr2d(1,1), 1, fileTypeWrite, MPI_STATUS_IGNORE, ierr);
!!$    elseif (present(bufi2d)) then
!!$       call MPI_File_write_all(H_OUTPUT, bufi2d(1,1), 1, fileTypeWrite, MPI_STATUS_IGNORE, ierr);
!!$    elseif (present(bufr3d)) then
!!$       call MPI_File_write_all(H_OUTPUT, bufr3d(1,1,1), 1, fileTypeWrite, MPI_STATUS_IGNORE, ierr);
!!$    elseif (present(bufr3d)) then
!!$       call MPI_File_write_all(H_OUTPUT, bufr4d(1,1,1,1), 1, fileTypeWrite, MPI_STATUS_IGNORE, ierr);
!!$    else
!!$       print *,'Erro escrita: chamada sem tipo definido.'
!!$       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
!!$    endif
!!$
!!$    if (ierr .ne. MPI_SUCCESS) then
!!$       print *,'Erro escrita.'
!!$       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
!!$    end if
!!$
!!$    call MPI_Type_extent(fileTypeWriteGlobal, structureSize, ierr)
!!$
!!$    DISPW = DISPW + structureSize
!!$
!!$  END SUBROUTINE WRITE_BRAMS_DATA_NEW

  SUBROUTINE BRAMS_SET_VIEW_REAL(fileTypeWriteGlobal)

    implicit none

    integer, intent(in) :: fileTypeWriteGlobal

    integer :: ierr

    call MPI_File_set_view(H_OUTPUT, DISPW, MPI_REAL, fileTypeWriteGlobal, &
         "native", MPI_INFO_NULL, ierr);
    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro set_view output.'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

  END SUBROUTINE BRAMS_SET_VIEW_REAL


  SUBROUTINE BRAMS_SET_VIEW_INTEGER(fileTypeWriteGlobal)

    implicit none

    integer, intent(in) :: fileTypeWriteGlobal

    integer :: ierr

    call MPI_File_set_view(H_OUTPUT, DISPW, MPI_INTEGER, fileTypeWriteGlobal, &
         "native", MPI_INFO_NULL, ierr);
    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro set_view output.'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

  END SUBROUTINE BRAMS_SET_VIEW_INTEGER


  SUBROUTINE UPDATE_DISPW(fileTypeWriteGlobal)

    implicit none

    integer, intent(in) :: fileTypeWriteGlobal

    INTEGER :: ierr, structureSize

    call MPI_Type_extent(fileTypeWriteGlobal, structureSize, ierr)

    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro MPI_Type_extent.'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

    DISPW = DISPW + structureSize

  END SUBROUTINE UPDATE_DISPW


  SUBROUTINE WRITE_BRAMS_DATA_2D(bufr2d, fileTypeWrite, fileTypeWriteGlobal)

    implicit none

    integer, intent(in) :: fileTypeWrite, fileTypeWriteGlobal
    real, intent(in)    :: bufr2d(:,:)

    INTEGER :: ierr

    call BRAMS_SET_VIEW_REAL(fileTypeWriteGlobal)

    call MPI_File_write_all(H_OUTPUT, bufr2d(1,1), 1, fileTypeWrite, &
         MPI_STATUS_IGNORE, ierr);

    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro escrita.'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

    call UPDATE_DISPW(fileTypeWriteGlobal)

  END SUBROUTINE WRITE_BRAMS_DATA_2D


  SUBROUTINE WRITE_BRAMS_DATA_2D_I(bufi2d, fileTypeWrite, fileTypeWriteGlobal)

    implicit none

    integer, intent(in) :: fileTypeWrite, fileTypeWriteGlobal
    integer, intent(in) :: bufi2d(:,:)

    INTEGER :: ierr

    call BRAMS_SET_VIEW_INTEGER(fileTypeWriteGlobal)

    call MPI_File_write_all(H_OUTPUT, bufi2d(1,1), 1, &
         fileTypeWrite, MPI_STATUS_IGNORE, ierr);

    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro escrita.'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

    call UPDATE_DISPW(fileTypeWriteGlobal)

  END SUBROUTINE WRITE_BRAMS_DATA_2D_I


  SUBROUTINE WRITE_BRAMS_DATA_3D(bufr3d, fileTypeWrite, fileTypeWriteGlobal)

    implicit none

    integer, intent(in) :: fileTypeWrite, fileTypeWriteGlobal
    real, intent(in)    :: bufr3d(:,:,:)

    INTEGER :: ierr

    call BRAMS_SET_VIEW_REAL(fileTypeWriteGlobal)

    call MPI_File_write_all(H_OUTPUT, bufr3d(1,1,1), 1, fileTypeWrite, &
         MPI_STATUS_IGNORE, ierr);

    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro escrita.'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

    call UPDATE_DISPW(fileTypeWriteGlobal)

  END SUBROUTINE WRITE_BRAMS_DATA_3D


  SUBROUTINE WRITE_BRAMS_DATA_4D(bufr4d, fileTypeWrite, fileTypeWriteGlobal)

    implicit none

    integer, intent(in) :: fileTypeWrite, fileTypeWriteGlobal
    real, intent(in)    :: bufr4d(:,:,:,:)

    INTEGER :: ierr

    call BRAMS_SET_VIEW_REAL(fileTypeWriteGlobal)

    call MPI_File_write_all(H_OUTPUT, bufr4d(1,1,1,1), 1, fileTypeWrite, &
         MPI_STATUS_IGNORE, ierr);

    if (ierr .ne. MPI_SUCCESS) then
       print *,'Erro escrita.'
       call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr)
    end if

    call UPDATE_DISPW(fileTypeWriteGlobal)

  END SUBROUTINE WRITE_BRAMS_DATA_4D



!!$$$$$$$$$$$$$$$$$$ FIM DA ESCRITA DO ARQUIVO DE SAIDA $$$$$$$$$$$$$$$$$$

!!$$$$$$$$$$$$$$$$$$ INICIO ROTINA RETIRADA DO BRAMS PARA DIVISAO DE DOMINIO $$$$$$$$$$$$$$$$$$
!!$  subroutine decomp_par(nxp,nyp,nodes,work,ixb,ixe,iyb,iye)
!!$
!!$    implicit none
!!$    ! Arguments:
!!$    integer, intent(in)  :: nxp
!!$    integer, intent(in)  :: nyp
!!$    integer, intent(in)  :: nodes
!!$    real,    intent(in)  :: work(nxp,nyp)
!!$    integer, intent(out) :: ixb(nodes)
!!$    integer, intent(out) :: ixe(nodes)
!!$    integer, intent(out) :: iyb(nodes)
!!$    integer, intent(out) :: iye(nodes)
!!$    ! Local variables:
!!$    real :: workrow(nyp)
!!$    real :: workcol(nxp)
!!$    real :: workload(nodes)
!!$    real :: workblock(nodes)
!!$    integer :: jrows(nodes)
!!$    integer :: jrow(nodes)
!!$    integer :: nblocks(nodes)
!!$    real :: relspeed(nodes)
!!$
!!$    integer :: inode,i,j,islab,jnodes,nslabs,min_blocks,nbigslabs,iblock &
!!$         ,jnode,knode
!!$    real :: anodes,aslabs,totspeed,workdom,workaccum,worksofar &
!!$         ,slabspeed,workslab
!!$
!!$    ! default relspeed = 1.0 for nodes of uniform speed.
!!$    relspeed = 1.0
!!$
!!$    ! This routine decomposes grid domains of size (nnxp,nnyp) into a number,
!!$    ! specified by nodes, of rectangular subdomains.  The convention is followed
!!$    ! that any internal boundaries (between subdomains) that are parallel to
!!$    ! the x-axis run continuously across the full domain, while boundaries
!!$    ! parallel to the y-axis may or may not run the full distance across the
!!$    ! domain.  For convenience, regions of the domain bounded by adjacent
!!$    ! east-west internal boundaries are termed "slabs", while smaller divisions
!!$    ! within each slab are termed "blocks".  Each block is required to have
!!$    ! a minimum dimension of 6 by 6 grid cells.  If this cannot be satisfied
!!$    ! with the given input parameters, the subroutine stops.
!!$
!!$
!!$    ! Estimate the number of slabs to be used (aslabs), and compute a final
!!$    ! nearest integer value (nslabs) which is limited to allowable values.
!!$    ! Zero out array for accumulating number of columns for each node.
!!$
!!$    anodes = float(nodes)
!!$    aslabs = sqrt(anodes * float(nyp) / float(nxp))
!!$    nslabs = min(nodes,max(1,nint(aslabs)))
!!$
!!$    !          print*, 'nslabs',nslabs
!!$
!!$    totspeed = 0.
!!$    do inode = 1,nodes
!!$       ixe(inode) = 0
!!$       totspeed = totspeed + relspeed(inode)
!!$    enddo
!!$
!!$    !          print*, 'totspeed',totspeed
!!$
!!$    ! Compute total work load over each row and over entire domain.
!!$
!!$    workdom = 0.
!!$    do j = 1,nyp
!!$       workrow(j) = 0.
!!$       do i = 1,nxp
!!$          workrow(j) = workrow(j) + work(i,j)
!!$       enddo
!!$       workdom = workdom + workrow(j)
!!$
!!$       !          print*, 'j,workdom,workrow(j)',j,workdom,workrow(j)
!!$
!!$    enddo
!!$    workrow(2) = workrow(2) + workrow(1)
!!$    workrow(nyp-1) = workrow(nyp-1) + workrow(nyp)
!!$
!!$    ! Determine number of blocks and the average workload for each slab.
!!$
!!$    min_blocks = nodes / nslabs
!!$    nbigslabs = nodes - min_blocks * nslabs
!!$    inode = 0
!!$    do islab = 1,nslabs
!!$       workload(islab) = 0.
!!$       nblocks(islab) = min_blocks
!!$       if (islab .le. nbigslabs) nblocks(islab) = min_blocks + 1
!!$       do iblock = 1,nblocks(islab)
!!$          inode = inode + 1
!!$          workload(islab) = workload(islab)  &
!!$               + workdom * relspeed(inode) / totspeed
!!$
!!$          !           print*, 'islab,iblock,workload(islab),workdom,inode'
!!$          !           print*,  islab,iblock,workload(islab),workdom,inode
!!$
!!$       enddo
!!$    enddo
!!$
!!$    ! Assign all j-rows to their respective slabs in a way that balances the work
!!$    ! load among slabs according to their respective numbers of nodes (blocks).
!!$    ! The array jrows counts the number of rows in each slab, and the array
!!$    ! jrow is the index of the southernmost row in each slab.
!!$
!!$    do islab = 1,nslabs
!!$       jrows(islab) = 0
!!$    enddo
!!$
!!$    workaccum = 0.
!!$    worksofar = 0.
!!$    islab = 0
!!$
!!$    do j = 2,nyp-1
!!$       workaccum = workaccum + workrow(j)
!!$       if (workaccum - .5 * workrow(j) .gt. worksofar .and.  &
!!$            islab .lt. nslabs) then
!!$          islab = islab + 1
!!$          jrow(islab) = j
!!$          worksofar = worksofar + workload(islab)
!!$       endif
!!$       jrows(islab) = jrows(islab) + 1
!!$    enddo
!!$
!!$    inode = 0
!!$    jnode = 0
!!$    knode = 0
!!$    do islab = 1,nslabs
!!$
!!$       ! Compute the total work load for each slab and for each i-column in the
!!$       ! slab.
!!$
!!$       slabspeed = 0.
!!$       workslab = 0.
!!$       do i = 1,nxp
!!$          workcol(i) = 0.
!!$          do j = jrow(islab),jrow(islab)+jrows(islab)-1
!!$             workcol(i) = workcol(i) + work(i,j)
!!$          enddo
!!$          workslab = workslab + workcol(i)
!!$       enddo
!!$       workcol(2) = workcol(2) + workcol(1)
!!$       workcol(nxp-1) = workcol(nxp-1) + workcol(nxp)
!!$
!!$       ! Determine average workload for each block.
!!$
!!$       do iblock = 1,nblocks(islab)
!!$          jnode = jnode + 1
!!$          slabspeed = slabspeed + relspeed(jnode)
!!$
!!$          !           print*, 'r1:iblock,jnode,slabspeed,relspeed(jnode)'
!!$          !           print*,     iblock,jnode,slabspeed,relspeed(jnode)
!!$
!!$       enddo
!!$       do iblock = 1,nblocks(islab)
!!$          knode = knode + 1
!!$          workblock(iblock) = workslab  &
!!$               * relspeed(knode) / slabspeed
!!$
!!$          !       print*, 'islab,iblock,workblock,workslab,relspeed,slabspeed'
!!$          !       print*, islab,iblock,workblock(iblock),workslab,relspeed(knode)
!!$          !     +       ,slabspeed
!!$          !       print*, 'knode',knode
!!$
!!$       enddo
!!$
!!$       ! Assign the i-columns of each slab to their respective blocks in a way that
!!$       ! balances the work load among the blocks.  The array ncols counts the number
!!$       ! of i-columns on each node, and the array ncol is the index of the
!!$       ! westernmost i-column on each node.
!!$
!!$       workaccum = 0.
!!$       worksofar = 0.
!!$
!!$       iblock = 0
!!$       do i = 2,nxp-1
!!$          workaccum = workaccum + workcol(i)
!!$
!!$          !        print*, 'islab',islab
!!$          !        print*, 'i,workaccum,workcol(i),worksofar,iblock,nblocks'
!!$          !        print*, i,workaccum,workcol(i),worksofar,iblock,nblocks(islab)
!!$
!!$          if (workaccum - .5 * workcol(i) .gt. worksofar .and.  &
!!$               iblock .lt. nblocks(islab)) then
!!$             iblock = iblock + 1
!!$
!!$             !ccccccc defining node variables here ccccccccccccccccccccccccccccccccccc
!!$             inode = inode + 1
!!$             iyb(inode) = jrow(islab)
!!$             ixb(inode) = i
!!$             iye(inode) = iyb(inode) + jrows(islab) - 1
!!$
!!$             !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$
!!$             worksofar = worksofar + workblock(iblock)
!!$          endif
!!$          ixe(inode) = ixe(inode) + 1
!!$       enddo
!!$    enddo
!!$
!!$    !ccccccc defining node variable here ccccccccccccccccccccccccccccccccccc
!!$    do jnode = 1,nodes
!!$       ixe(jnode) = ixb(jnode) + ixe(jnode) - 1
!!$
!!$       !           print*, 'jnode,ixb,ixe',jnode,ixb(jnode),ixe(jnode)
!!$       !     +        ,(ixe(jnode)-ixb(jnode)+1),(iye(jnode)-iyb(jnode)+1)
!!$       !     +        ,(ixe(jnode)-ixb(jnode)+1)*(iye(jnode)-iyb(jnode)+1)
!!$
!!$    enddo
!!$    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$
!!$    ! Check to make sure that each subdomain has at least 2 interior 
!!$    ! rows and columns.
!!$
!!$    do jnode = 1,nodes
!!$       if (iye(jnode) - iyb(jnode) .lt. 1 .or.  &
!!$            ixe(jnode) - ixb(jnode) .lt. 1) then
!!$          print*, 'grid:',nxp,nyp,'  subdomain too small on node ',jnode
!!$          print*, '(ixb,ixe,iyb,iye) = '  &
!!$               ,ixb(jnode),ixe(jnode),iyb(jnode),iye(jnode)
!!$          !stop 'small_nodes'
!!$       endif
!!$    enddo
!!$
!!$    return
!!$  end subroutine decomp_par

!!$$$$$$$$$$$$$$$$$$ FIM ROTINA RETIRADA DO BRAMS PARA DIVISAO DE DOMINIO $$$$$$$$$$$$$$$$$$

end module MPI_IO_ENGINE
