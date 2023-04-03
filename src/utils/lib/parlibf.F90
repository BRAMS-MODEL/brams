module ParLib
  implicit none
  private
  public :: parf_init_mpi
  public :: parf_exit_mpi
  public :: parf_error
  public :: parf_send_int
  public :: parf_get_int
  public :: parf_send_real
  public :: parf_get_real
  public :: parf_send_block
  public :: parf_get_block
  public :: parf_get_block_any
  public :: parf_get_noblock
  public :: parf_get_noblock_real
  public :: parf_get_noblock_integer
  public :: parf_send_noblock
  public :: parf_send_noblock_real
  public :: parf_send_noblock_integer
  public :: parf_wait
  public :: parf_wait_nostatus
  public :: parf_wait_any_nostatus
  public :: parf_wait_all_nostatus
  public :: parf_pack
  private:: parf_pack_int_1d
  private:: parf_pack_int_2d
  private:: parf_pack_int_scalar
  private:: parf_pack_real_1d
  private:: parf_pack_real_scalar
  private:: parf_pack_char
  public :: parf_unpack
  private:: parf_unpack_int_1d
  private:: parf_unpack_int_1d_2d
  private:: parf_unpack_int_scalar
  private:: parf_unpack_real_1d
  private:: parf_unpack_real_scalar
  private:: parf_unpack_char
  public :: parf_barrier
  public :: parf_pack_max_size
  public :: parf_bcast
  private:: parf_bcast_logical_scalar
  private:: parf_bcast_real_scalar
  private:: parf_bcast_real_1d
  private:: parf_bcast_real_2d
  private:: parf_bcast_real_3d
  private:: parf_bcast_real_4d
  private:: parf_bcast_int_1d
  private:: parf_bcast_int_scalar
  private:: parf_bcast_char
  private:: parf_bcast_char_vec
  public :: parf_minloc
  public :: parf_reduce_max
  public :: parf_allreduce_max
  private:: parf_allreduce_sum_scalar
  private:: parf_allreduce_sum_vector
  public :: parf_allreduce_sum
  public :: parf_GatherAllChunks
  public :: parf_GatherPostSfc
  public :: parf_gather_real

  ! Interfaces
  interface parf_pack
     module procedure &
          parf_pack_int_1d, &
          parf_pack_int_2d, &
          parf_pack_int_scalar, &
          parf_pack_real_1d, &
          parf_pack_real_scalar, &
          parf_pack_char
  end interface

  interface parf_unpack
     module procedure &
          parf_unpack_int_1d, &
          parf_unpack_int_1d_2d, &
          parf_unpack_int_scalar, &
          parf_unpack_real_1d, &
          parf_unpack_real_scalar, &
          parf_unpack_char
  end interface

  interface parf_bcast
     module procedure &
          parf_bcast_logical_scalar, &
          parf_bcast_real_1d, &
          parf_bcast_real_2d, &
          parf_bcast_real_3d, &
          parf_bcast_real_4d, &
          parf_bcast_real_scalar, &
          parf_bcast_int_1d, &
          parf_bcast_int_scalar, &
          parf_bcast_char, &
          parf_bcast_char_vec
  end interface

  interface parf_allreduce_sum
     module procedure &
          parf_allreduce_sum_scalar, &
          parf_allreduce_sum_vector
  end interface

  include "i8.h"
#if defined (RAMS_MPI)
  include 'mpif.h'
#endif

contains


#if defined (RAMS_MPI)

  subroutine parf_init_mpi(mchnum, nmachs, master_num)
    integer, intent(out) :: mchnum
    integer, intent(out) :: nmachs
    integer, intent(out) :: master_num
    integer :: ierr

    master_num=0
    call MPI_Init(ierr)
    if(ierr /= MPI_SUCCESS) then
       call fatal_error("Error in MPI_Init")
    endif

    call MPI_Comm_rank(MPI_COMM_WORLD, mchnum, ierr)
    if(ierr /= MPI_SUCCESS) then
       call fatal_error("Error in MPI_Comm_rank")
    endif

    call MPI_Comm_size(MPI_COMM_WORLD, nmachs, ierr)
    if(ierr /= MPI_SUCCESS) then
       call fatal_error("Error in MPI_Comm_size")
    endif
  end subroutine parf_init_mpi

  !--------------------------------------------------------------------

  subroutine parf_exit_mpi()
    integer :: ierr

    call MPI_Finalize(ierr)
    if(ierr /= MPI_SUCCESS) then
       call fatal_error("Error in MPI_Finalize")
    endif
  end subroutine parf_exit_mpi

  ! -------------------------------------------------------------------

  subroutine parf_error()
    integer :: ierr

    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    stop
  end subroutine parf_error

  !--------------------------------------------------------------------

  subroutine parf_send_int(buff, buff_len, dest_host, tag)
    integer(kind=i8), intent(in) :: buff_len
    integer,          intent(in) :: buff(buff_len)
    integer,          intent(in) :: dest_host
    integer,          intent(in) :: tag
    ! Local Variables:
    integer              :: ierr, ierr_b, rank
    character(len=20)    :: string

    call MPI_send(buff, buff_len, MPI_Integer, dest_host, tag, &
         MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_send_int: rank, ierr="//trim(string))
    endif
  end subroutine parf_send_int

  !--------------------------------------------------------------------

  subroutine parf_get_int(buff, buff_len, source_host, tag)
    integer(kind=i8), intent(in)  :: buff_len
    integer,          intent(out) :: buff(buff_len)
!srf    real,         intent(out) :: buff(buff_len)
    integer,          intent(in)  :: source_host
    integer,          intent(in)  :: tag
    ! Local Variables:
    integer              :: ierr, ierr_b, rank
    character(len=20)    :: string
    integer              :: status(MPI_STATUS_SIZE)

    call MPI_recv(buff, buff_len, MPI_Integer, source_host, tag, &
         MPI_COMM_WORLD, status, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_get_int: rank, ierr="//trim(string))
    endif
  end subroutine parf_get_int

  !--------------------------------------------------------------------

  subroutine parf_send_real(buff, buff_len, dest_host, tag)
    integer(kind=i8), intent(in) :: buff_len
    real,             intent(in) :: buff(buff_len)
    integer,          intent(in) :: dest_host
    integer,          intent(in) :: tag
    ! Local Variables:
    integer              :: ierr, ierr_b, rank
    character(len=20)    :: string

    call MPI_send(buff, buff_len, MPI_Real, dest_host, tag, &
         MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_send_real: rank, ierr="//trim(string))
    endif
  end subroutine parf_send_real

  !--------------------------------------------------------------------

  subroutine parf_get_real(buff, buff_len, source_host, tag)
    integer(kind=i8), intent(in)  :: buff_len
    real,             intent(out) :: buff(buff_len)
    integer,          intent(in)  :: source_host
    integer,          intent(in)  :: tag
    ! Local Variables:
    integer              :: ierr, ierr_b, rank
    character(len=20)    :: string
    integer              :: status(MPI_STATUS_SIZE)

    call MPI_recv(buff, buff_len, MPI_Real, source_host, tag, &
         MPI_COMM_WORLD, status, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_get_real: rank, ierr="//trim(string))
    endif
  end subroutine parf_get_real

  !--------------------------------------------------------------------

  subroutine parf_send_block(buff, buff_len, dest_host, tag)
    integer(kind=i8), intent(in) :: buff_len
    real,             intent(in) :: buff(buff_len)
    integer,          intent(in) :: dest_host
    integer,          intent(in) :: tag
    ! Local Variables:
    integer              :: ierr, ierr_b, rank
    character(len=20)    :: string

    call MPI_send(buff, buff_len, MPI_Packed, dest_host, tag, &
         MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_send_block: rank, ierr="//trim(string))
    endif
  end subroutine parf_send_block

  !--------------------------------------------------------------------

  subroutine parf_get_block(buff, buff_len, source_host, tag)
    integer(kind=i8), intent(in)  :: buff_len
    real,             intent(out) :: buff(buff_len)
    integer,          intent(in)  :: source_host
    integer,          intent(in)  :: tag
    ! Local Variables:
    integer              :: ierr, ierr_b, rank
    character(len=20)    :: string
    integer              :: status(MPI_STATUS_SIZE)

    call MPI_recv(buff, buff_len, MPI_Packed, source_host, tag, &
         MPI_COMM_WORLD, status, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_get_block: rank, ierr="//trim(string))
    endif
  end subroutine parf_get_block

  !--------------------------------------------------------------------

  subroutine parf_get_block_any(buff, buff_len, tag, host_msg)
    integer(kind=i8), intent(in)  :: buff_len
    real,             intent(out) :: buff(buff_len)
    integer,          intent(in)  :: tag
    integer,          intent(out) :: host_msg
    ! Local Variables:
    integer              :: ierr, ierr_b, rank
    character(len=20)    :: string
    integer              :: status(MPI_STATUS_SIZE)

    call MPI_recv(buff, buff_len, MPI_Packed, MPI_ANY_SOURCE, tag, &
         MPI_COMM_WORLD, status, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_get_block: rank, ierr="//trim(string))
    endif

    host_msg = status(MPI_SOURCE)
  end subroutine parf_get_block_any

  !--------------------------------------------------------------------

  subroutine parf_get_noblock(buff, buff_len, source_host, tag, request)
    integer(kind=i8), intent(in)  :: buff_len
    real,             intent(out) :: buff(buff_len)
    integer,          intent(in)  :: source_host
    integer,          intent(in)  :: tag
    integer,          intent(out) :: request
    ! Local Variables:
    integer              :: ierr, ierr_b, rank
    character(len=20)    :: string


    call MPI_Irecv(buff, buff_len, MPI_PACKED, source_host, tag, &
         MPI_COMM_WORLD, request, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8,X,I8)') rank, request, ierr
       call fatal_error("Error in parf_get_noblock: rank, request, ierr="//trim(string))
    endif
  end subroutine parf_get_noblock


  !--(DMK-CCATT-INI)--------------------------------------------
  !--------------------------------------------------------------------
  subroutine parf_get_noblock_real(buff, buff_len, source_host, tag, request)

    ! Arguments:
    integer, intent(in)  :: buff_len, source_host, tag
    real,    intent(out) :: buff(buff_len)
    integer, intent(out) :: request
    ! Local Variables:
    integer              :: ierr, ierr_b, rank
    character(len=20)    :: string


    call MPI_Irecv(buff, buff_len, MPI_REAL, source_host, tag, &
         MPI_COMM_WORLD, request, ierr)

    if(ierr /= MPI_Success) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8,X,I8)') rank, request, ierr
       call fatal_error("Error in parf_get_noblock: rank, request, ierr="//trim(string))
    endif


  end subroutine parf_get_noblock_real

  !--------------------------------------

  subroutine parf_get_noblock_integer(buff, buff_len, source_host, tag, request)

    ! Arguments:
    integer, intent(in)  :: buff_len, source_host, tag
    integer, intent(out) :: buff(buff_len)
    integer, intent(out) :: request
    ! Local Variables:
    integer              :: ierr, ierr_b, rank
    character(len=20)    :: string


    call MPI_Irecv(buff, buff_len, MPI_INTEGER, source_host, tag, &
         MPI_COMM_WORLD, request, ierr)

    if(ierr /= MPI_Success) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8,X,I8)') rank, request, ierr
       call fatal_error("Error in parf_get_noblock_integer: rank, request, ierr="//trim(string))
    endif
  end subroutine parf_get_noblock_integer

  !--------------------------------------

  subroutine parf_send_noblock_real(buff, buff_len, dest_host, tag, request)

    ! Arguments:
    integer, intent(in)  :: buff_len, dest_host, tag
    real,    intent(in)  :: buff(buff_len)
    integer, intent(out) :: request
    ! Local Variables:
    integer              :: ierr, ierr_b, rank
    character(len=20)    :: string

    call MPI_Isend(buff, buff_len, MPI_REAL, dest_host, tag, &
         MPI_COMM_WORLD, request, ierr)

    if(ierr /= MPI_Success) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8,X,I8)') rank, request, ierr
       call fatal_error("Error in parf_send_noblock: rank, request, ierr="//&
            trim(string))
    endif

  end subroutine parf_send_noblock_real


  subroutine parf_send_noblock_integer(buff, buff_len, dest_host, tag, request)

    ! Arguments:
    integer, intent(in)  :: buff_len, dest_host, tag
    integer, intent(in)  :: buff(buff_len)
    integer, intent(out) :: request
    ! Local Variables:
    integer              :: ierr, ierr_b, rank
    character(len=20)    :: string

    call MPI_Isend(buff, buff_len, MPI_INTEGER, dest_host, tag, &
         MPI_COMM_WORLD, request, ierr)

    if(ierr /= MPI_Success) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8,X,I8)') rank, request, ierr
       call fatal_error("Error in parf_send_noblock_integer: rank, request, ierr="//&
            trim(string))
    endif
  end subroutine parf_send_noblock_integer

  !--------------------------------------

  subroutine parf_wait_any_nostatus(total,request,number)

    ! Arguments:
    integer, intent(inout) :: request(*)
    integer,intent(in) :: total

    integer, intent(out) :: number
    ! Local Variables:
    integer                :: status(MPI_STATUS_SIZE)
    integer                :: ierr, ierr_b, rank

    call MPI_Waitany(total, request, number, status, ierr)
    if(ierr /= MPI_Success) call fatal_error("Error in parf_wait_any")

  end subroutine parf_wait_any_nostatus

  !--------------------------------------------------------------------

  subroutine parf_wait_all_nostatus(total,request)

    ! Arguments:
    integer, intent(inout) :: request(*)
    integer,intent(in) :: total
    ! Local Variables:
    integer                :: status(MPI_STATUS_SIZE,total)
    integer                :: ierr, ierr_b, rank

    call MPI_Waitall(total, request, status, ierr)
    if(ierr /= MPI_Success) call fatal_error("Error in parf_wait_all")

  end subroutine parf_wait_all_nostatus

  !--------------------------------------------------------------------
  !--(DMK-CCATT-FIM)--------------------------------------------

  !--------------------------------------------------------------------

  subroutine parf_send_noblock(buff, buff_len, dest_host, tag, request)
    integer(kind=i8), intent(in)  :: buff_len
    real,             intent(in)  :: buff(buff_len)
    integer,          intent(in)  :: dest_host
    integer,          intent(in)  :: tag
    integer,          intent(out) :: request
    ! Local Variables:
    integer              :: ierr, ierr_b, rank
    character(len=20)    :: string

    call MPI_Isend(buff, buff_len, MPI_Packed, dest_host, tag, &
         MPI_COMM_WORLD, request, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8,X,I8)') rank, request, ierr
       call fatal_error("Error in parf_send_noblock: rank, request, ierr="//&
            trim(string))
    endif
  end subroutine parf_send_noblock

  !--------------------------------------------------------------------

  subroutine parf_wait(request, status)
    integer, intent(inout) :: request
    integer, intent(out)   :: status(MPI_STATUS_SIZE)
    ! Local Variables:
    integer                :: ierr, ierr_b, rank
    character(len=20)      :: string

    call MPI_Wait(request, status, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8,X,I8)') rank, request, ierr
       call fatal_error("Error in parf_wait: rank, request, ierr="//trim(string))
    endif
  end subroutine parf_wait

  !--------------------------------------------------------------------

  subroutine parf_wait_nostatus(request)
    integer, intent(inout) :: request
    ! Local Variables:
    integer                :: status(MPI_STATUS_SIZE)
    integer                :: ierr, ierr_b, rank
    character(len=20)      :: string

    call MPI_Wait(request, status, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8,X,I8)') rank, request, ierr
       call fatal_error("Error in parf_wait_nostatus: rank, request, ierr="//&
            trim(string))
    endif
  end subroutine parf_wait_nostatus

  !--------------------------------------------------------------------

  subroutine parf_pack_int_1d(in_buff, in_buff_len, &
       out_buff, out_buff_len, position)
    integer(kind=i8), intent(in)    :: in_buff_len
    integer,          intent(in)    :: in_buff(in_buff_len)
    integer(kind=i8), intent(in)    :: out_buff_len
    real,             intent(inout) :: out_buff(out_buff_len)
    integer(kind=i8), intent(inout) :: position
    ! Local Variables:
    integer                :: ierr, ierr_b, rank
    character(len=20)      :: string

    call MPI_Pack(in_buff, in_buff_len, MPI_INTEGER, out_buff, out_buff_len, &
         position, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_pack_int_1d: rank, ierr="//trim(string))
    endif
  end subroutine parf_pack_int_1d

  !--------------------------------------------------------------------

  subroutine parf_pack_int_2d(in_buff, in_buff_lenx, in_buff_leny, &
       out_buff, out_buff_len, position)
    integer(kind=i8), intent(in)    :: in_buff_lenx, in_buff_leny
    integer,          intent(in)    :: in_buff(in_buff_lenx, in_buff_leny)
    integer(kind=i8), intent(in)    :: out_buff_len
    real,             intent(inout) :: out_buff(out_buff_len)
    integer(kind=i8), intent(inout) :: position
    ! Local Variables:
    integer                :: ierr, ierr_b, rank
    character(len=20)      :: string
    integer(kind=i8)       :: in_buff_len

    in_buff_len = (in_buff_lenx*in_buff_leny)

    call MPI_Pack(in_buff, in_buff_len, MPI_INTEGER, out_buff, out_buff_len, &
         position, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_pack_int_2d: rank, ierr="//trim(string))
    endif
  end subroutine parf_pack_int_2d

  !--------------------------------------------------------------------

  subroutine parf_pack_int_scalar(in_buff, out_buff, out_buff_len, position)
    integer,          intent(in)    :: in_buff
    integer(kind=i8), intent(in)    :: out_buff_len
    real,             intent(inout) :: out_buff(out_buff_len)
    integer(kind=i8), intent(inout) :: position
    ! Local Variables:
    integer                :: ierr, ierr_b, rank
    character(len=20)      :: string

    call MPI_Pack(in_buff, 1, MPI_INTEGER, out_buff, out_buff_len, &
         position, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_pack_int_scalar: rank, ierr="//&
            trim(string))
    endif
  end subroutine parf_pack_int_scalar

  !--------------------------------------------------------------------

  subroutine parf_pack_real_1d(in_buff, in_buff_len, &
       out_buff, out_buff_len, position)
    integer(kind=i8), intent(in)    :: in_buff_len
    real,             intent(in)    :: in_buff(in_buff_len)
    integer(kind=i8), intent(in)    :: out_buff_len
    real,             intent(inout) :: out_buff(out_buff_len)
    integer(kind=i8), intent(inout) :: position
    ! Local Variables:
    integer                :: ierr, ierr_b, rank
    character(len=20)      :: string

    call MPI_Pack(in_buff, in_buff_len, MPI_REAL, out_buff, out_buff_len, &
         position, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_pack_real_1d: rank, ierr="//trim(string))
    endif
  end subroutine parf_pack_real_1d

  !--------------------------------------------------------------------

  subroutine parf_pack_real_scalar(in_buff, out_buff, out_buff_len, position)
    real,             intent(in)    :: in_buff
    integer(kind=i8), intent(in)    :: out_buff_len
    real,             intent(inout) :: out_buff(out_buff_len)
    integer(kind=i8), intent(inout) :: position
    ! Local Variables:
    integer                :: ierr, ierr_b, rank
    character(len=20)      :: string

    call MPI_Pack(in_buff, 1, MPI_REAL, out_buff, out_buff_len, &
         position, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_pack_real_scalar: rank, ierr="//&
            trim(string))
    endif
  end subroutine parf_pack_real_scalar

  !--------------------------------------------------------------------

  subroutine parf_pack_char(in_buff, in_buff_len, &
       out_buff, out_buff_len, position)
    integer(kind=i8),           intent(in)    :: in_buff_len
    character(len=in_buff_len), intent(in)    :: in_buff
    integer(kind=i8),           intent(in)    :: out_buff_len
    real,                       intent(inout) :: out_buff(out_buff_len)
    integer(kind=i8),           intent(inout) :: position
    ! Local Variables:
    integer                :: ierr, ierr_b, rank
    character(len=20)      :: string

    call MPI_Pack(in_buff, in_buff_len, MPI_CHARACTER, out_buff, out_buff_len, &
         position, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_pack_char: rank, ierr="//trim(string))
    endif
  end subroutine parf_pack_char

  !--------------------------------------------------------------------

  subroutine parf_unpack_int_1d(in_buff, in_buff_len, position, &
       out_buff, out_buff_len)
    integer(kind=i8), intent(in)    :: in_buff_len
    real,             intent(in)    :: in_buff(in_buff_len)
    integer(kind=i8), intent(inout) :: position
    integer(kind=i8), intent(in)    :: out_buff_len
    integer,          intent(inout) :: out_buff(out_buff_len)
    ! Local Variables:
    integer                :: ierr, ierr_b, rank
    character(len=20)      :: string

    call MPI_Unpack(in_buff, in_buff_len, position, out_buff, out_buff_len, &
         MPI_INTEGER, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_unpack_int_1d: rank, ierr="//&
            trim(string))
    endif
  end subroutine parf_unpack_int_1d

  !--------------------------------------------------------------------

  subroutine parf_unpack_int_1d_2d(in_buff, in_buff_len, position, &
       out_buff, out_buff_lenx, out_buff_leny)
    integer(kind=i8), intent(in)    :: in_buff_len
    real,             intent(in)    :: in_buff(in_buff_len)
    integer(kind=i8), intent(inout) :: position
    integer(kind=i8), intent(in)    :: out_buff_lenx, out_buff_leny
    integer,          intent(inout) :: out_buff(out_buff_lenx,out_buff_leny)
    ! Local Variables:
    integer                :: ierr, ierr_b, rank
    character(len=20)      :: string
    integer(kind=i8)       :: out_buff_len

    out_buff_len = (out_buff_lenx*out_buff_leny)

    call MPI_Unpack(in_buff, in_buff_len, position, out_buff, out_buff_len, &
         MPI_INTEGER, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_unpack_int_1d: rank, ierr="//&
            trim(string))
    endif
  end subroutine parf_unpack_int_1d_2d

  !--------------------------------------------------------------------

  subroutine parf_unpack_int_scalar(in_buff, position, out_buff, out_buff_len)
    real,             intent(in)    :: in_buff
    integer(kind=i8), intent(inout) :: position
    integer(kind=i8), intent(in)    :: out_buff_len
    integer,          intent(inout) :: out_buff(out_buff_len)
    ! Local Variables:
    integer                :: ierr, ierr_b, rank
    character(len=20)      :: string

    call MPI_Unpack(in_buff, 1, position, out_buff, out_buff_len, &
         MPI_INTEGER, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_unpack_int_scalar: rank, ierr="//&
            trim(string))
    endif
  end subroutine parf_unpack_int_scalar

  !--------------------------------------------------------------------

  subroutine parf_unpack_real_1d(in_buff, in_buff_len, position, &
       out_buff, out_buff_len)
    integer(kind=i8), intent(in)    :: in_buff_len
    real,             intent(in)    :: in_buff(in_buff_len)
    integer(kind=i8), intent(inout) :: position
    integer(kind=i8), intent(in)    :: out_buff_len
    real,             intent(inout) :: out_buff(out_buff_len)
    ! Local Variables:
    integer                :: ierr, ierr_b, rank
    character(len=20)      :: string

    call MPI_Unpack(in_buff, in_buff_len, position, out_buff, out_buff_len, &
         MPI_REAL, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_unpack_real_1d: rank, ierr="//&
            trim(string))
    endif
  end subroutine parf_unpack_real_1d

  !--------------------------------------------------------------------

  subroutine parf_unpack_real_scalar(in_buff, position, out_buff, out_buff_len)
    real,             intent(in)    :: in_buff
    integer(kind=i8), intent(inout) :: position
    integer(kind=i8), intent(in)    :: out_buff_len
    real,             intent(inout) :: out_buff(out_buff_len)
    ! Local Variables:
    integer                :: ierr, ierr_b, rank
    character(len=20)      :: string

    call MPI_Unpack(in_buff, 1, position, out_buff, out_buff_len, &
         MPI_REAL, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_unpack_real_scalar: rank, ierr="//&
            trim(string))
    endif
  end subroutine parf_unpack_real_scalar

  !--------------------------------------------------------------------

  subroutine parf_unpack_char(in_buff, in_buff_len, position, &
       out_buff, out_buff_len)
    ! Arguments
    integer, intent(in)                         :: in_buff_len, out_buff_len
    real,    intent(in)                         :: in_buff(in_buff_len)
    character (len=out_buff_len), intent(inout) :: out_buff
    integer, intent(inout)                      :: position
    ! Local Variables:
    integer                :: ierr, ierr_b, rank
    character(len=20)      :: string

    call MPI_Unpack(in_buff, in_buff_len, position, out_buff, out_buff_len, &
         MPI_CHARACTER, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_unpack_char: rank, ierr="//trim(string))
    endif
  end subroutine parf_unpack_char

  !--------------------------------------------------------------------

  subroutine parf_barrier(ibarrier)
    integer, intent(in) :: ibarrier
    ! Local Variables:
    integer             :: ierr, ierr_b, rank
    character(len=20)   :: string

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    if (ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string,FMT='(I6.6,X,I6,X,I8)') rank, ibarrier, ierr
       call fatal_error("Error in parf_barrier: "//trim(string))
    endif
  end subroutine parf_barrier

  !--------------------------------------------------------------------

  subroutine parf_pack_max_size(len, max_size)
    integer(kind=i8), intent(in)  :: len
    integer(kind=i8), intent(out) :: max_size
    ! Local Variables:
    integer              :: mpi_int_size, mpi_real_size
    integer              :: ierr, rank
    character(len=20)    :: string

    call MPI_Pack_size(len, MPI_INTEGER, MPI_COMM_WORLD, mpi_int_size, ierr)
    if (ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
       write(string,FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_pack_max_size:A: "//trim(string))
    endif
    call MPI_Pack_size(len, MPI_REAL, MPI_COMM_WORLD, mpi_real_size, ierr)
    if (ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
       write(string,FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_pack_max_size:B: "//trim(string))
    endif

    max_size = max(mpi_real_size, mpi_int_size)

    ! Extra lenght
    max_size = max_size + max_size/2
  end subroutine parf_pack_max_size

  !--------------------------------------------------------------------

  subroutine parf_bcast_logical_scalar(buff, source_host)
    integer, intent(in) :: source_host
    logical, intent(in) :: buff
    ! Local Variables:
    integer             :: ierr, ierr_b, rank
    character(len=20)   :: string

    call MPI_BCAST(buff, 1, MPI_LOGICAL, source_host, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string, FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_bcast_logical_scalar: rank, ierr="//&
            trim(string))
    endif
  end subroutine parf_bcast_logical_scalar

  !--------------------------------------------------------------------

  subroutine parf_bcast_real_3d(buff, buff_lenz, buff_lenx, buff_leny, &
       source_host)
    integer(kind=i8), intent(in)    :: buff_lenz, buff_lenx, buff_leny
    integer,          intent(in)    :: source_host
    real,             intent(inout) :: buff(buff_lenz, buff_lenx, buff_leny)
    ! Local Variables:
    integer(kind=i8)             :: buff_len
    integer                      :: ierr, ierr_b, rank
    character(len=20)            :: string

    buff_len = buff_lenz*buff_lenx*buff_leny

    call MPI_BCAST(buff, buff_len, MPI_REAL, source_host, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string, FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_bcast_real_3d: rank, ierr="//&
            trim(string))
    endif
  end subroutine parf_bcast_real_3d

  !--------------------------------------------------------------------

  subroutine parf_bcast_real_4d(buff, buff_lenz, buff_lenx, buff_leny, &
       buff_lenk, source_host)
    integer(kind=i8), intent(in)    :: buff_lenz, buff_lenx, buff_leny, &
         buff_lenk
    integer,          intent(in)    :: source_host
    real,             intent(inout) :: buff(buff_lenz, buff_lenx, buff_leny, &
         buff_lenk)
    ! Local Variables:
    integer(kind=i8)             :: buff_len
    integer                      :: ierr, ierr_b, rank
    character(len=20)            :: string

    buff_len = buff_lenz*buff_lenx*buff_leny*buff_lenk

    call MPI_BCAST(buff, buff_len, MPI_REAL, source_host, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string, FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_bcast_real_4d: rank, ierr="//&
            trim(string))
    endif
  end subroutine parf_bcast_real_4d

  !--------------------------------------------------------------------

  subroutine parf_bcast_real_2d(buff, buff_lenx, buff_leny, source_host)
    integer(kind=i8), intent(in)    :: buff_lenx, buff_leny
    integer,          intent(in)    :: source_host
    real,             intent(inout) :: buff(buff_lenx, buff_leny)
    ! Local Variables:
    integer(kind=i8)             :: buff_len
    integer                      :: ierr, ierr_b, rank
    character(len=20)            :: string

    buff_len = buff_lenx*buff_leny

    call MPI_BCAST(buff, buff_len, MPI_REAL, source_host, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string, FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_bcast_real_2d: rank, ierr="//&
            trim(string))
    endif
  end subroutine parf_bcast_real_2d

  !--------------------------------------------------------------------

  subroutine parf_bcast_real_1d(buff, buff_len, source_host)
    integer(kind=i8), intent(in)    :: buff_len
    integer,          intent(in)    :: source_host
    real,             intent(inout) :: buff(buff_len)
    ! Local Variables:
    integer                      :: ierr, ierr_b, rank
    character(len=20)            :: string

    call MPI_BCAST(buff, buff_len, MPI_REAL, source_host, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string, FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_bcast_real_1d: rank, ierr="//&
            trim(string))
    endif
  end subroutine parf_bcast_real_1d

  !--------------------------------------------------------------------

  subroutine parf_bcast_real_scalar(buff, source_host)
    integer,          intent(in) :: source_host
    real,             intent(in) :: buff
    ! Local Variables:
    integer                      :: ierr, ierr_b, rank
    character(len=20)            :: string

    call MPI_BCAST(buff, 1, MPI_REAL, source_host, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string, FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_bcast_real_scalar: rank, ierr="//&
            trim(string))
    endif
  end subroutine parf_bcast_real_scalar

  !--------------------------------------------------------------------

  subroutine parf_bcast_int_1d(buff, buff_len, source_host)
    integer(kind=i8), intent(in)    :: buff_len
    integer,          intent(inout) :: buff(buff_len)
    integer,          intent(in)    :: source_host
    ! Local Variables:
    integer                      :: ierr, ierr_b, rank
    character(len=20)            :: string

    call MPI_BCAST(buff, buff_len, MPI_INTEGER, source_host, &
         MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string, FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_bcast_int_1d: rank, ierr="//trim(string))
    endif
  end subroutine parf_bcast_int_1d

  !--------------------------------------------------------------------

  subroutine parf_bcast_int_scalar(buff, source_host)
    integer,          intent(in) :: buff
    integer,          intent(in) :: source_host
    ! Local Variables:
    integer                      :: ierr, ierr_b, rank
    character(len=20)            :: string

    call MPI_BCAST(buff, 1, MPI_INTEGER, source_host, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string, FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_bcast_int_scalar: rank, ierr="//&
            trim(string))
    endif
  end subroutine parf_bcast_int_scalar

  !--------------------------------------------------------------------

  subroutine parf_bcast_char(buff, buff_len, source_host)
    integer(kind=i8),        intent(in)    :: buff_len
    integer,                 intent(in)    :: source_host
    character(len=buff_len), intent(inout) :: buff
    ! Local Variables:
    integer                                :: ierr, ierr_b, rank
    character(len=20)                      :: string

    call MPI_BCAST(buff, buff_len, MPI_CHARACTER, source_host, &
         MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string, FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_bcast_char: rank, ierr="//trim(string))
    endif
  end subroutine parf_bcast_char

  !--------------------------------------------------------------------

  subroutine parf_bcast_char_vec(buff, buff_len, buff_size, source_host)
    integer(kind=i8),        intent(in)    :: buff_len
    integer(kind=i8),        intent(in)    :: buff_size
    integer,                 intent(in)    :: source_host
    character(len=buff_len), intent(inout) :: buff(buff_size)
    ! Local Variables:
    integer                                :: ierr, ierr_b, rank
    character(len=20)                      :: string

    call MPI_BCAST(buff, buff_len*buff_size, MPI_CHARACTER, source_host, &
         MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string, FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_bcast_char_vec: rank, ierr="//&
            trim(string))
    endif
  end subroutine parf_bcast_char_vec

  !--------------------------------------------------------------------

  subroutine parf_minloc(bufin, bufout)
    real, intent(in)  :: bufin(2)
    real, intent(out) :: bufout(2)
    ! Local Variables:
    integer           :: ierr, ierr_b, rank
    character(len=20) :: string

    call MPI_ALLREDUCE(bufin, bufout, 1, MPI_2REAL, &
         MPI_MINLOC, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string, FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_minloc: rank, ierr="//trim(string))
    endif
  end subroutine parf_minloc

  !--------------------------------------------------------------------

  subroutine parf_reduce_max(bufin, bufout, buff_len, master_num)
    integer(kind=i8), intent(in)  :: buff_len
    real,             intent(in)  :: bufin(buff_len)
    real,             intent(out) :: bufout(buff_len)
    integer,          intent(in)  :: master_num
    ! Local Variables:
    integer           :: ierr, ierr_b, rank
    character(len=20) :: string

    call MPI_REDUCE(bufin, bufout, buff_len, MPI_REAL, &
         MPI_MAX, master_num, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string, FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_reduce_max: rank, ierr="//trim(string))
    endif
  end subroutine parf_reduce_max

  !--------------------------------------------------------------------

  subroutine parf_allreduce_max(bufin, bufout, buff_len)
    integer(kind=i8), intent(in)  :: buff_len
    real,             intent(in)  :: bufin(buff_len)
    real,             intent(out) :: bufout(buff_len)
    ! Local Variables:
    integer           :: ierr, ierr_b, rank
    character(len=20) :: string

    call MPI_ALLREDUCE(bufin, bufout, buff_len, MPI_REAL, &
         MPI_MAX, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string, FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_allreduce_max: rank, ierr="//trim(string))
    endif
  end subroutine parf_allreduce_max

  !--------------------------------------------------------------------

  subroutine parf_allreduce_sum_vector(bufin, bufout, buff_len)
    integer(kind=i8), intent(in)  :: buff_len
    real,             intent(in)  :: bufin(buff_len)
    real,             intent(out) :: bufout(buff_len)
    ! Local Variables:
    integer           :: ierr, ierr_b, rank
    character(len=20) :: string

    call MPI_ALLREDUCE(bufin, bufout, buff_len, MPI_REAL, &
         MPI_SUM, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string, FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_allreduce_sum_vector: rank, ierr="//&
            trim(string))
    endif
  end subroutine parf_allreduce_sum_vector

  !--------------------------------------------------------------------

  subroutine parf_allreduce_sum_scalar(bufin, bufout)
    real,             intent(in)  :: bufin
    real,             intent(out) :: bufout
    ! Local Variables:
    integer           :: ierr, ierr_b, rank
    character(len=20) :: string

    call MPI_ALLREDUCE(bufin, bufout, 1, MPI_REAL, &
         MPI_SUM, MPI_COMM_WORLD, ierr)

    if(ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string, FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_allreduce_sum_scalar: rank, ierr="//&
            trim(string))
    endif
  end subroutine parf_allreduce_sum_scalar

  !--------------------------------------------------------------------

  subroutine parf_Gather_real(bufin, bufout, master_num, nmachs)
    real,    intent(in ) :: bufin
    real,    intent(out) :: bufout(0:nmachs)
    integer, intent(in ) :: master_num
    integer, intent(in ) :: nmachs
    ! Local Variables:
    integer :: ierr, ierr_b, rank
    character(len=20) :: string

    call MPI_GATHER (bufin, 1, MPI_REAL, bufout, 1, MPI_REAL, &
         master_num, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string, FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_gather_real:rank, ierr="//trim(string))
    end if
  end subroutine Parf_Gather_real

  !--------------------------------------------------------------------

  subroutine parf_GatherAllChunks(LocalChunk, thisChunkSize, idim_type, &
       localSize, disp, gathered, sizeGathered, master_num, nmachs)
    integer, intent(in ) :: thisChunkSize
    real,    intent(in ) :: LocalChunk(thisChunkSize)
    integer, intent(in ) :: idim_type
    integer, intent(in ) :: localSize(nmachs,2:7)
    integer, intent(in ) :: disp(nmachs,2:7)
    integer, intent(in ) :: sizeGathered
    real,    intent(out) :: gathered(sizeGathered)
    integer, intent(in ) :: master_num
    integer, intent(in ) :: nmachs
    ! Local Variables:
    integer :: ierr, ierr_b, rank
    character(len=20) :: string

    ! gather a field

    call MPI_GATHERV(LocalChunk, thisChunkSize, MPI_REAL, &
         gathered, localSize(1,idim_type), disp(1,idim_type), MPI_REAL, &
         master_num, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string, FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error("Error in parf_gatherv:rank, ierr="//trim(string))
    end if
  end subroutine Parf_GatherAllChunks


  subroutine parf_GatherPostSfc(LocalChunk, localSize, disp, &
       gathered, master_num)
    real,    intent(in ) :: LocalChunk(:)
    integer, intent(in ) :: localSize(:)
    integer, intent(in ) :: disp(:)
    real,    intent(out) :: gathered(:)
    integer, intent(in ) :: master_num

    ! gathers one surface of a packed post field

    ! Local Variables:
    integer :: ierr, ierr_b, rank
    character(len=20) :: string
    character(len=*), parameter :: h="**(parf_GatherPostSfc)**"

    ! gather a field

    call MPI_GATHERV(LocalChunk, size(LocalChunk), MPI_REAL, &
         gathered, localSize, disp, MPI_REAL, &
         master_num, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
       call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_b)
       write(string, FMT='(I6.6,X,I8)') rank, ierr
       call fatal_error(h//" rank, ierr="//trim(string))
    end if
  end subroutine parf_GatherPostSfc

  ! **************************************************************

#else
  ! No MPI - Serial Version

  !--------------------------------------------------------------------

  subroutine parf_init_mpi(mchnum, nmachs, master_num)
    integer, intent(out) :: mchnum
    integer, intent(out) :: nmachs
    integer, intent(out) :: master_num

    master_num=0
    mchnum  = 0
    nmachs = 1
  end subroutine parf_init_mpi

  !--------------------------------------------------------------------

  subroutine parf_exit_mpi()
  end subroutine parf_exit_mpi

  ! -------------------------------------------------------------------

  subroutine parf_error()
    stop
  end subroutine parf_error

  !--------------------------------------------------------------------

  subroutine parf_send_int(buff, buff_len, dest_host, tag)
    integer(kind=i8) :: buff_len
    integer          :: buff(buff_len)
    integer          :: dest_host
    integer          :: tag
  end subroutine parf_send_int

  !--------------------------------------------------------------------

  subroutine parf_get_int(buff, buff_len, source_host, tag)
    integer(kind=i8)  :: buff_len
    integer           :: source_host, tag
    real              :: buff(buff_len)
  end subroutine parf_get_int

  !--------------------------------------------------------------------

  subroutine parf_send_real(buff, buff_len, dest_host, tag)
    integer(kind=i8) :: buff_len
    real             :: buff(buff_len)
    integer          :: dest_host
    integer          :: tag
  end subroutine parf_send_real

  !--------------------------------------------------------------------

  subroutine parf_get_real(buff, buff_len, source_host, tag)
    integer(kind=i8) :: buff_len
    real             :: buff(buff_len)
    integer          :: source_host
    integer          :: tag
  end subroutine parf_get_real

  !--------------------------------------------------------------------

  subroutine parf_send_block(buff, buff_len, dest_host, tag)
    integer(kind=i8) :: buff_len
    real             :: buff(buff_len)
    integer          :: dest_host
    integer          :: tag
  end subroutine parf_send_block

  !--------------------------------------------------------------------

  subroutine parf_get_block(buff, buff_len, source_host, tag)
    integer(kind=i8)  :: buff_len
    real              :: buff(buff_len)
    integer           :: source_host
    integer           :: tag
  end subroutine parf_get_block

  !--------------------------------------------------------------------

  subroutine parf_get_block_any(buff, buff_len, tag, host_msg)
    integer(kind=i8)  :: buff_len
    real              :: buff(buff_len)
    integer           :: tag
    integer           :: host_msg
  end subroutine parf_get_block_any

  !--------------------------------------------------------------------

  subroutine parf_get_noblock(buff, buff_len, source_host, tag, request)
    integer(kind=i8)  :: buff_len
    real              :: buff(buff_len)
    integer           :: source_host
    integer           :: tag
    integer           :: request

    request = -1
  end subroutine parf_get_noblock

  !--------------------------------------------------------------------

  subroutine parf_send_noblock(buff, buff_len, dest_host, tag, request)
    integer(kind=i8) :: buff_len
    real             :: buff(buff_len)
    integer          :: dest_host
    integer          :: tag
    integer          :: request

    request = -1
  end subroutine parf_send_noblock

  !--------------------------------------------------------------------

  subroutine parf_wait(request, status)
    integer :: request
    integer :: status(*)
  end subroutine parf_wait

  !--------------------------------------------------------------------

  subroutine parf_wait_nostatus(request)
    integer :: request
  end subroutine parf_wait_nostatus

  !--------------------------------------------------------------------

  subroutine parf_pack_int_1d(in_buff, in_buff_len, &
       out_buff, out_buff_len, position)
    integer(kind=i8) :: in_buff_len
    integer          :: in_buff(in_buff_len)
    integer(kind=i8) :: out_buff_len
    real             :: out_buff(out_buff_len)
    integer(kind=i8) :: position
  end subroutine parf_pack_int_1d

  !--------------------------------------------------------------------

  subroutine parf_pack_int_2d(in_buff, in_buff_lenx, in_buff_leny, &
       out_buff, out_buff_len, position)
    integer(kind=i8), intent(in)    :: in_buff_lenx, in_buff_leny
    integer,          intent(in)    :: in_buff(in_buff_lenx, in_buff_leny)
    integer(kind=i8), intent(in)    :: out_buff_len
    real,             intent(inout) :: out_buff(out_buff_len)
    integer(kind=i8), intent(inout) :: position
  end subroutine parf_pack_int_2d

  !--------------------------------------------------------------------

  subroutine parf_pack_int_scalar(in_buff, out_buff, out_buff_len, position)
    integer          :: in_buff
    integer(kind=i8) :: out_buff_len
    real             :: out_buff(out_buff_len)
    integer(kind=i8) :: position
  end subroutine parf_pack_int_scalar

  !--------------------------------------------------------------------

  subroutine parf_pack_real_1d(in_buff, in_buff_len, &
       out_buff, out_buff_len, position)
    integer(kind=i8) :: in_buff_len
    real             :: in_buff(in_buff_len)
    integer(kind=i8) :: out_buff_len
    real             :: out_buff(out_buff_len)
    integer(kind=i8) :: position
  end subroutine parf_pack_real_1d

  !--------------------------------------------------------------------

  subroutine parf_pack_real_scalar(in_buff, out_buff, out_buff_len, position)
    real             :: in_buff
    integer(kind=i8) :: out_buff_len
    real             :: out_buff(out_buff_len)
    integer(kind=i8) :: position
  end subroutine parf_pack_real_scalar

  !--------------------------------------------------------------------

  subroutine parf_pack_char(in_buff, in_buff_len, &
       out_buff, out_buff_len, position)
    integer(kind=i8)           :: in_buff_len
    character(len=in_buff_len) :: in_buff
    integer(kind=i8)           :: out_buff_len
    real                       :: out_buff(out_buff_len)
    integer(kind=i8)           :: position
  end subroutine parf_pack_char

  !--------------------------------------------------------------------

  subroutine parf_unpack_int_1d(in_buff, in_buff_len, position, &
       out_buff, out_buff_len)
    integer(kind=i8) :: in_buff_len
    real             :: in_buff(in_buff_len)
    integer(kind=i8) :: position
    integer(kind=i8) :: out_buff_len
    integer          :: out_buff(out_buff_len)
  end subroutine parf_unpack_int_1d

  !--------------------------------------------------------------------

  subroutine parf_unpack_int_1d_2d(in_buff, in_buff_len, position, &
       out_buff, out_buff_lenx, out_buff_leny)
    integer(kind=i8), intent(in)    :: in_buff_len
    real,             intent(in)    :: in_buff(in_buff_len)
    integer(kind=i8), intent(inout) :: position
    integer(kind=i8), intent(in)    :: out_buff_lenx, out_buff_leny
    integer,          intent(inout) :: out_buff(out_buff_lenx,out_buff_leny)
  end subroutine parf_unpack_int_1d_2d

  !--------------------------------------------------------------------

  subroutine parf_unpack_int_scalar(in_buff, position, out_buff, out_buff_len)
    real,             intent(in)    :: in_buff
    integer(kind=i8), intent(inout) :: position
    integer(kind=i8), intent(in)    :: out_buff_len
    integer,          intent(inout) :: out_buff(out_buff_len)
  end subroutine parf_unpack_int_scalar

  !--------------------------------------------------------------------

  subroutine parf_unpack_real_1d(in_buff, in_buff_len, position, &
       out_buff, out_buff_len)
    integer(kind=i8) :: in_buff_len
    real             :: in_buff(in_buff_len)
    integer(kind=i8) :: position
    integer(kind=i8) :: out_buff_len
    real             :: out_buff(out_buff_len)
  end subroutine parf_unpack_real_1d

  !--------------------------------------------------------------------

  subroutine parf_unpack_real_scalar(in_buff, position, out_buff, out_buff_len)
    real,             intent(in)    :: in_buff
    integer(kind=i8), intent(inout) :: position
    integer(kind=i8), intent(in)    :: out_buff_len
    real,             intent(inout) :: out_buff(out_buff_len)
  end subroutine parf_unpack_real_scalar

  !--------------------------------------------------------------------

  subroutine parf_unpack_char(in_buff, in_buff_len, position, &
       out_buff, out_buff_len)
    integer(kind=i8)           :: in_buff_len
    character(len=in_buff_len) :: in_buff
    integer(kind=i8)           :: position
    integer(kind=i8)           :: out_buff_len
    real                       :: out_buff(out_buff_len)
  end subroutine parf_unpack_char

  !--------------------------------------------------------------------

  subroutine parf_barrier(ibarrier)
    integer :: ibarrier
  end subroutine parf_barrier

  !--------------------------------------------------------------------

  subroutine parf_pack_max_size(len, max_size)
    integer(kind=i8)  :: len
    integer(kind=i8) :: max_size

    max_size = 0
  end subroutine parf_pack_max_size

  !--------------------------------------------------------------------

  subroutine parf_bcast_real_3d(buff, buff_lenz, buff_lenx, buff_leny, &
       source_host)
    integer(kind=i8) :: buff_lenz, buff_lenx, buff_leny
    integer          :: source_host
    real             :: buff(buff_lenz,buff_lenx,buff_leny)
  end subroutine parf_bcast_real_3d

  !--------------------------------------------------------------------

  subroutine parf_bcast_real_4d(buff, buff_lenz, buff_lenx, buff_leny, &
       buff_lenk, source_host)
    integer(kind=i8), intent(in)    :: buff_lenz, buff_lenx, buff_leny, &
         buff_lenk
    integer,          intent(in)    :: source_host
    real,             intent(inout) :: buff(buff_lenz, buff_lenx, buff_leny, &
         buff_lenk)
  end subroutine parf_bcast_real_4d

  !--------------------------------------------------------------------

  subroutine parf_bcast_real_2d(buff, buff_lenx, buff_leny, source_host)
    integer(kind=i8) :: buff_lenx, buff_leny
    integer          :: source_host
    real             :: buff(buff_lenx,buff_leny)
  end subroutine parf_bcast_real_2d

  !--------------------------------------------------------------------

  subroutine parf_bcast_real_1d(buff, buff_len, source_host)
    integer(kind=i8) :: buff_len
    integer          :: source_host
    real             :: buff(buff_len)
  end subroutine parf_bcast_real_1d

  !--------------------------------------------------------------------

  subroutine parf_bcast_real_scalar(buff, source_host)
    integer          :: source_host
    real             :: buff
  end subroutine parf_bcast_real_scalar

  !--------------------------------------------------------------------

  subroutine parf_bcast_int_1d(buff, buff_len, source_host)
    integer(kind=i8) :: buff_len
    integer          :: buff(buff_len)
    integer          :: source_host
  end subroutine parf_bcast_int_1d

  !--------------------------------------------------------------------

  subroutine parf_bcast_int_scalar(buff, source_host)
    integer          :: buff
    integer          :: source_host
  end subroutine parf_bcast_int_scalar

  !--------------------------------------------------------------------

  subroutine parf_bcast_char(buff, buff_len, source_host)
    integer(kind=i8)        :: buff_len
    integer                 :: source_host
    character(len=buff_len) :: buff
  end subroutine parf_bcast_char

  !--------------------------------------------------------------------

  subroutine parf_bcast_char_vec(buff, buff_len, buff_size, source_host)
    integer(kind=i8),        intent(in)    :: buff_len
    integer(kind=i8),        intent(in)    :: buff_size
    integer,                 intent(in)    :: source_host
    character(len=buff_len), intent(inout) :: buff(buff_size)
  end subroutine parf_bcast_char_vec

  !--------------------------------------------------------------------

  subroutine parf_minloc(bufin, bufout)
    real  :: bufin(2)
    real :: bufout(2)
  end subroutine parf_minloc

  !--------------------------------------------------------------------

  subroutine parf_reduce_max(bufin, bufout, buff_len, machreduce)
    integer(kind=i8), intent(in)  :: buff_len
    real,             intent(in)  :: bufin(buff_len)
    real,             intent(out) :: bufout(buff_len)
    integer,          intent(in)  :: machreduce
  end subroutine parf_reduce_max

  !--------------------------------------------------------------------

  subroutine parf_allreduce_max(bufin, bufout, buff_len)
    integer(kind=i8), intent(in)  :: buff_len
    real,             intent(in)  :: bufin(buff_len)
    real,             intent(out) :: bufout(buff_len)
  end subroutine parf_allreduce_max

  !--------------------------------------------------------------------

  subroutine parf_allreduce_sum_vector(bufin, bufout, buff_len)
    integer(kind=i8), intent(in)  :: buff_len
    real,             intent(in)  :: bufin(buff_len)
    real,             intent(out) :: bufout(buff_len)
  end subroutine parf_allreduce_sum_vector

  !--------------------------------------------------------------------

  subroutine parf_allreduce_sum_scalar(bufin, bufout)
    real,             intent(in)  :: bufin
    real,             intent(out) :: bufout
  end subroutine parf_allreduce_sum_scalar

  !--------------------------------------------------------------------

  subroutine parf_Gather_real(bufin, bufout, master_num, nmachs)
    real,    intent(in ) :: bufin
    real,    intent(out) :: bufout(0:nmachs)
    integer, intent(in ) :: master_num
    integer, intent(in ) :: nmachs
  end subroutine Parf_Gather_real

  !--------------------------------------------------------------------

  subroutine parf_GatherAllChunks(LocalChunk, thisChunkSize, idim_type, &
       localSize, disp, gathered, sizeGathered, master_num, nmachs)
    integer :: thisChunkSize
    real    :: LocalChunk(thisChunkSize)
    integer :: idim_type
    integer :: localSize(nmachs,2:7)
    integer :: disp(nmachs,2:7)
    integer :: sizeGathered
    real    :: gathered(sizeGathered)
    integer :: master_num
    integer :: nmachs
  end subroutine Parf_GatherAllChunks

#endif

end module ParLib

subroutine fatal_error(msg)

  ! fatal_error: exception handling, to be invoked
  !   whenever a fatal error occurs.
  !   dumps a message at stdout and halts execution

  use ParLib, only: &
       parf_error
  use dump, only: &
    dumpMessage
  implicit none

  include "constants.f90"
  character(len=*), intent(in) :: msg

  iErrNumber=dumpMessage(c_tty,c_yes,'Undefined - Old',modelVersion,c_fatal,msg)

  call parf_error()
  stop

end subroutine fatal_error
