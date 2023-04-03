module hdf5_parallel_engine

  USE HDF5 ! This module contains all necessary modules 

  implicit none

  private

  integer(hid_t) :: hdf_file_id
!!$  integer(hid_t) :: hdf_file_id2
  integer, parameter :: io_maxdims = 4
  integer :: rank
  integer(hsize_t) :: localdims(io_maxdims)
  integer(hsize_t) :: offsetNoGhost(io_maxdims)
  integer(hsize_t) :: dimsNoGhost(io_maxdims)
  integer(hsize_t) :: globaldims(io_maxdims)
  integer(hsize_t) :: offsetThisProc(io_maxdims)

  public :: hdf_open
  public :: hdf_close
  public :: hdf_output
  public :: dims_output
  public :: hdf_parallel_output

contains

  subroutine hdf_open(hdf_name, parallel)
!!$    use hdf5
!!$    use an_header, only: hdf_fid
    implicit none
    include 'mpif.h'

    character(len=256), intent(in) :: hdf_name
!!$    integer(hid_t), intent(out) :: hdf_file_id
    logical, intent(in) :: parallel

    integer(hid_t) :: plist_id
    integer :: error

    character(len=*), parameter :: h="**(hdf_open)**"

    CALL h5open_f(error)
    if (error .ne. 0) then
       call fatal_error(h//" Error in h5open_f")
    endif

    if (parallel) then
       call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
       if (error .ne. 0) then
          call fatal_error(h//" Error in h5pcreate_f")
       endif

       call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
       if (error .ne. 0) then
          call fatal_error(h//" Error in h5pset_fapl_mpio_f")
       endif

       call h5fcreate_f(hdf_name, H5F_ACC_TRUNC_F, hdf_file_id, error, access_prp = plist_id)
       if (error .ne. 0) then
          call fatal_error(h//" Error in h5fcreate_f")
       endif

       call h5pclose_f(plist_id, error)
       if (error .ne. 0) then
          call fatal_error(h//" Error in h5pclose_f")
       endif
    else
       call h5fcreate_f(hdf_name, H5F_ACC_TRUNC_F, hdf_file_id, error)
       if (error .ne. 0) then
          call fatal_error(h//" Error in h5fcreate_f")
       endif
    endif
!!$    hdf_fid = hdf_file_id
  end subroutine hdf_open


  subroutine hdf_close()
!!$    use hdf5
    implicit none

!!$    integer(hid_t), intent(in) :: hdf_file_id
    integer :: error
    character(len=*), parameter :: h="**(hdf_close)**"


    call h5fclose_f(hdf_file_id, error)
    if (error .ne. 0) then
       call fatal_error(h//" Error in h5fclose_f")
    endif

  end subroutine hdf_close


  ! this routine assumes that the field was already rearranged (IJK, not KIJ)
  subroutine hdf_output(master_num, mchnum, nnxp, nnyp, nnzp, nzg, npatch, &
       nzs, nwave, varn, idim_type, fullField)
!!$    use node_mod, only: master_num, mchnum
!!$    use mem_grid, only: nnxp, nnyp, nnzp, nzg, npatch, nzs
!!$    use mem_aerad, only: nwave
!!$    use hdf5
    implicit none
    integer, intent(in) :: master_num, mchnum
    integer, intent(in) :: nnxp, nnyp, nnzp, nzg, npatch, nzs, nwave
!!$    integer(HID_T), intent(in) :: hdf_file_id
    character(len=16), intent(in) :: varn
    integer, intent(in) :: idim_type
!!$    integer, intent(in) :: ng
    integer :: sizeFullField ! I can compute that just like the GlobalSizes subroutine, it doesn't need to be a input parameter
!!$  real, intent(in) :: fullField(sizeFullField)
    real, intent(in) :: fullField(:)

    integer(HID_T) :: dspace_id
    integer(HID_T) :: dset_id
    integer :: error
    integer, parameter :: MAX_DIMS = 4
    integer(HSIZE_T), dimension(MAX_DIMS) :: dims
!!$    integer :: rank

    if (mchnum == master_num) then

       select case (idim_type)
       case (2)
          dims = (/ nnxp, nnyp, 0, 0 /)
          sizeFullField = nnxp*nnyp
          rank = 2
       case (3)
          !dims = (/ nnzp, nnxp, nnyp, 0 /) ! NOT REARRANGED
          dims = (/ nnxp, nnyp, nnzp, 0 /)
          sizeFullField = nnzp*nnxp*nnyp
          rank = 3
       case (4)
          !dims = (/ nzg, nnxp, nnyp, npatch /) ! NOT REARRANGED
          dims = (/ nnxp, nnyp, nzg, npatch /)
          sizeFullField = nzg*nnxp*nnyp*npatch
          rank = 4
       case (5)
          !dims = (/ nzs, nnxp, nnyp, npatch /) ! NOT REARRANGED
          dims = (/ nnxp, nnyp, nzs, npatch /)
          sizeFullField = nzs*nnxp*nnyp*npatch
          rank = 4
       case (6)
          dims = (/ nnxp, nnyp, npatch, 0 /)
          sizeFullField = nnxp*nnyp*npatch
          rank = 3
       case (7)
          dims = (/ nnxp, nnyp, nwave, 0 /)
          sizeFullField = nnxp*nnyp*nwave
          rank = 3
       end select

       call h5screate_simple_f(rank, dims, dspace_id, error)
       call h5dcreate_f(hdf_file_id, varn, H5T_NATIVE_REAL, dspace_id, &
            dset_id, error)
       call h5dwrite_f(dset_id, H5T_NATIVE_REAL, fullField, dims, error)
       call h5dclose_f(dset_id, error)
       call h5sclose_f(dspace_id, error)
    end if
  end subroutine hdf_output


  subroutine dims_output(nnxp, nnyp, nnzp, nzg, npatch, nzs, nwave, &
       nodemxp, nodemyp, nodeia, nodeiz, nodeja, nodejz, nodeibcon, ixb, iyb, &
       idim_type &
!!$       localdims, offsetNoGhost, dimsNoGhost, globaldims, offsetThisProc &
       )
!!$    use node_mod, only: nodemxp, nodemyp, nodeia, nodeiz, nodeja, nodejz, nodeibcon, ixb, iyb
!!$    use mem_grid, only: nnxp, nnyp, nnzp, nzg, npatch, nzs
!!$    use mem_aerad, only: nwave
!!$    use hdf5
    implicit none
!!$    integer, intent(in) :: mynum
!!$    integer, intent(in) :: ng
    integer, intent(in) :: nnxp, nnyp, nnzp, nzg, npatch, nzs, nwave
    integer, intent(in) :: nodemxp, nodemyp
    integer, intent(in) :: nodeia, nodeiz, nodeja, nodejz
    integer, intent(in) :: nodeibcon
    integer, intent(in) :: ixb, iyb
    integer, intent(in) :: idim_type
!!$    integer, intent(out) :: rank
!!$    integer, intent(in) :: io_maxdims

!!$  integer(hsize_t), intent(out) :: localdims(rank)
!!$  integer(hsize_t), intent(out) :: offsetNoGhost(rank)
!!$  integer(hsize_t), intent(out) :: dimsNoGhost(rank)
!!$  integer(hsize_t), intent(out) :: globaldims(rank)
!!$  integer(hsize_t), intent(out) :: offsetThisProc(rank)
!!$  integer(hsize_t), intent(inout) :: localdims(:)
!!$  integer(hsize_t), intent(inout) :: offsetNoGhost(:)
!!$  integer(hsize_t), intent(inout) :: dimsNoGhost(:)
!!$  integer(hsize_t), intent(inout) :: globaldims(:)
!!$  integer(hsize_t), intent(inout) :: offsetThisProc(:)
!!$    integer(hsize_t), intent(out) :: localdims(io_maxdims)
!!$    integer(hsize_t), intent(out) :: offsetNoGhost(io_maxdims)
!!$    integer(hsize_t), intent(out) :: dimsNoGhost(io_maxdims)
!!$    integer(hsize_t), intent(out) :: globaldims(io_maxdims)
!!$    integer(hsize_t), intent(out) :: offsetThisProc(io_maxdims)

    integer :: onoghost(2), dnoghost(2) ! offset and dimension without ghostzone and independent of idim_type
    integer :: othisproc(2) ! offset in the entire domain independent of idim_type

    logical, parameter :: dumpLocal=.false.

    character(len=7)            :: cProc
    character(len=16)           :: varn
    character(len=8)            :: c0, c1, c2
    character(len=*), parameter :: h="**(dims_output)**" 

    onoghost(1) = nodeia - 1
    onoghost(2) = nodeja - 1
    dnoghost(1) = nodemxp
    dnoghost(2) = nodemyp
    othisproc(1) = ixb - 1
    othisproc(2) = iyb - 1
    if (btest(nodeibcon, 0)) then
       onoghost(1) = onoghost(1) - 1
       othisproc(1) = othisproc(1) - 1
    else
       dnoghost(1) = dnoghost(1) - 1
    endif
    if (.not. btest(nodeibcon, 1)) then
       dnoghost(1) = dnoghost(1) - 1
    endif
    if (btest(nodeibcon, 2)) then
       onoghost(2) = onoghost(2) - 1
       othisproc(2) = othisproc(2) - 1
    else
       dnoghost(2) = dnoghost(2) - 1
    endif
    if (.not. btest(nodeibcon, 3)) then
       dnoghost(2) = dnoghost(2) - 1
    endif

    if (dumpLocal) then
       write(*,*) h//" checkpoint A: idim_type=", idim_type
       call flush(6)
    end if

    select case (idim_type)
    case (2) 
       rank = 2
       localdims(1) = nodemxp
       localdims(2) = nodemyp
       offsetNoGhost(1) = onoghost(1)
       offsetNoGhost(2) = onoghost(2)
       dimsNoGhost(1) = dnoghost(1)
       dimsNoGhost(2) = dnoghost(2)
       globaldims(1) = nnxp
       globaldims(2) = nnyp
       offsetThisProc(1) =  othisproc(1)
       offsetThisProc(2) = othisproc(2)
    case (3)
       rank = 3

       if (dumpLocal) then
          write(*,*) h//" checkpoint B0: idim_type, rank=", idim_type, rank
          write(*,*) h//" size(localdims)=", size(localdims)
          call flush(6)
       end if

       localdims(1) = nnzp
       localdims(2) = nodemxp
       localdims(3) = nodemyp

       if (dumpLocal) then
          write(*,*) h//" checkpoint B1-localdims: idim_type, rank=", idim_type, rank
          call flush(6)
       end if

       offsetNoGhost(1) = 0
       offsetNoGhost(2) = onoghost(1)
       offsetNoGhost(3) = onoghost(2)

       if (dumpLocal) then
          write(*,*) h//" checkpoint B2-offsetNoGhost: idim_type, rank=", idim_type, rank
          call flush(6)
       end if

       dimsNoGhost(1) = localdims(1)
       dimsNoGhost(2) = dnoghost(1)
       dimsNoGhost(3) = dnoghost(2)

       if (dumpLocal) then
          write(*,*) h//" checkpoint B3-dimsNoGhost: idim_type, rank=", idim_type, rank
          call flush(6)
       end if

       globaldims(1) = nnzp
       globaldims(2) = nnxp
       globaldims(3) = nnyp

       if (dumpLocal) then
          write(*,*) h//" checkpoint B4-globaldims: idim_type, rank=", idim_type, rank
          call flush(6)
       end if

       offsetThisProc(1) = 0
       offsetThisProc(2) = othisproc(1)
       offsetThisProc(3) = othisproc(2)

       if (dumpLocal) then
          write(*,*) h//" checkpoint B-fim: idim_type, rank=", idim_type, rank
          call flush(6)
       end if

    case (4)
       rank = 4
       localdims(1) = nzg
       localdims(2) = nodemxp
       localdims(3) = nodemyp
       localdims(4) = npatch
       offsetNoGhost(1) = 0
       offsetNoGhost(2) = onoghost(1)
       offsetNoGhost(3) = onoghost(2)
       offsetNoGhost(4) = 0
       dimsNoGhost(1) = localdims(1)
       dimsNoGhost(2) = dnoghost(1)
       dimsNoGhost(3) = dnoghost(2)
       dimsNoGhost(4) = localdims(4)
       globaldims(1) = nzg
       globaldims(2) = nnxp
       globaldims(3) = nnyp
       globaldims(4) = npatch
       offsetThisProc(1) = 0
       offsetThisProc(2) = othisproc(1)
       offsetThisProc(3) = othisproc(2)
       offsetThisProc(4) = 0
    case (5)
       rank = 4
       localdims(1) = nzs
       localdims(2) = nodemxp
       localdims(3) = nodemyp
       localdims(4) = npatch
       offsetNoGhost(1) = 0
       offsetNoGhost(2) = onoghost(1)
       offsetNoGhost(3) = onoghost(2)
       offsetNoGhost(4) = 0
       dimsNoGhost(1) = localdims(1)
       dimsNoGhost(2) = dnoghost(1)
       dimsNoGhost(3) = dnoghost(2)
       dimsNoGhost(4) = localdims(4)
       globaldims(1) = nzs
       globaldims(2) = nnxp
       globaldims(3) = nnyp
       globaldims(4) = npatch
       offsetThisProc(1) = 0
       offsetThisProc(2) = othisproc(1)
       offsetThisProc(3) = othisproc(2)
       offsetThisProc(4) = 0
    case (6)
       rank = 3
       localdims(1) = nodemxp
       localdims(2) = nodemyp
       localdims(3) = npatch
       offsetNoGhost(1) = onoghost(1)
       offsetNoGhost(2) = onoghost(2)
       offsetNoGhost(3) = 0
       dimsNoGhost(1) = dnoghost(1)
       dimsNoGhost(2) = dnoghost(2)
       dimsNoGhost(3) = localdims(3)
       globaldims(1) = nnxp
       globaldims(2) = nnyp
       globaldims(3) = npatch
       offsetThisProc(1) = othisproc(1)
       offsetThisProc(2) = othisproc(2)
       offsetThisProc(3) = 0
    case (7)
       rank = 3
       localdims(1) = nodemxp
       localdims(2) = nodemyp
       localdims(3) = nwave
       offsetNoGhost(1) = onoghost(1)
       offsetNoGhost(2) = onoghost(2)
       offsetNoGhost(3) = 0
       dimsNoGhost(1) = dnoghost(1)
       dimsNoGhost(2) = dnoghost(2)
       dimsNoGhost(3) = localdims(3)
       globaldims(1) = nnxp
       globaldims(2) = nnyp
       globaldims(3) = nwave
       offsetThisProc(1) = othisproc(1)
       offsetThisProc(2) = othisproc(2)
       offsetThisProc(3) = 0
    case default
       print*, "ERROR idim_type"
       call fatal_error(h//" ERROR idim_type")
    end select
  end subroutine dims_output


  subroutine hdf_parallel_output(varn, localchunk, localsize)
!!$    , &
!!$       localdims, offsetNoGhost, dimsNoGhost, &
!!$       globaldims, offsetThisProc)
!!$    use hdf5
    implicit none
    include 'mpif.h'
    include "i8.h"

!!$    integer(hid_t), intent(in) :: hdf_file_id
    character(len=16), intent(in) :: varn
!!$    integer, intent(in) :: rank
    integer(kind=i8), intent(in) :: localsize
    real, intent(in) :: localchunk(localsize)
!!$    integer(hsize_t), intent(in) :: localdims(rank)
!!$    integer(hsize_t), intent(in) :: offsetNoGhost(rank)
!!$    integer(hsize_t), intent(in) :: dimsNoGhost(rank)
!!$    integer(hsize_t), intent(in) :: globaldims(rank)
!!$    integer(hsize_t), intent(in) :: offsetThisProc(rank)

    integer :: error
    integer(hid_t) :: memspace, filespace
    integer(hid_t) :: dset_id
    integer(hid_t) :: plist_id

    integer :: mpi_rank

    character(len=*), parameter :: h="**(hdf_parallel_output)**"


    call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, error)

    !TODO A better approach is to create this hyperslab only once. Here it is created
    !     each time a field is written 

    call h5screate_simple_f(rank, localdims, memspace, error)
    if (error .ne. 0) then
       print*, "error h5screate_simple_f 1 - mpi rank: ", mpi_rank
       call fatal_error(h//" Error h5screate_simple_f")
    endif

    call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
         start=offsetNoGhost, count=dimsNoGhost, &
         hdferr=error)
    if (error .ne. 0) then
       print*, "error h5sselect_hyperslab_f - mpi rank: ", mpi_rank
       call fatal_error(h//" Error h5sselect_hyperslab_f")
    endif

    call h5screate_simple_f(rank, globaldims, filespace, error)
    if (error .ne. 0) then
       print*, "error h5screate_simple_f 2 - mpi rank: ", mpi_rank
       call fatal_error(h//" Error h5screate_simple_f")
    endif

    call h5dcreate_f(hdf_file_id, varn, H5T_NATIVE_REAL, filespace, &
         dset_id, error)
    if (error .ne. 0) then
       print*, "error h5dcreate_f - mpi rank: ", mpi_rank
       call fatal_error(h//" Error h5dcreate_f")
    endif

    call h5sclose_f(filespace, error) ! I don't know why I have to do that, really!
    if (error .ne. 0) then
       print*, "error h5sclose_f - mpi rank: ", mpi_rank
       call fatal_error(h//" Error h5sclose_f")
    endif

    call h5dget_space_f(dset_id, filespace, error) ! neither this
    if (error .ne. 0) then
       print*, "error h5dget_space_f - mpi rank: ", mpi_rank
       call fatal_error(h//" Error h5dget_space_f")
    endif

    call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offsetThisProc, &
         dimsNoGhost, error)
    if (error .ne. 0) then
       print*, "error h5sselect_hyperslab_f - mpi rank: ", mpi_rank
       call fatal_error(h//" Error h5sselect_hyperslab_f")
    endif

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (error .ne. 0) then
       print*, "error h5pcreate_f - mpi rank: ", mpi_rank
       call fatal_error(h//" Error h5pcreate_f")
    endif

    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    if (error .ne. 0) then
       print*, "error h5pset_dxpl_mpio_f - mpi rank: ", mpi_rank
       call fatal_error(h//" Error h5pset_dxpl_mpio_f")
    endif

    call h5dwrite_f(dset_id, H5T_NATIVE_REAL, localchunk, localdims, error, &
         file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = plist_id)
    if (error .ne. 0) then
       print*, "error h5dwrite_f - mpi rank: ", mpi_rank
       call fatal_error(h//" Error h5dwrite_f")
    endif

    call h5dclose_f(dset_id, error)
    if (error .ne. 0) then
       print*, "error h5dclose_f - mpi rank: ", mpi_rank
       call fatal_error(h//" Error h5dclose_f")
    endif

    call h5sclose_f(filespace, error)
    if (error .ne. 0) then
       print*, "error h5sclose_f - mpi rank: ", mpi_rank
       call fatal_error(h//" Error h5sclose_f")
    endif

    call h5sclose_f(memspace, error)
    if (error .ne. 0) then
       print*, "error h5sclose_f - mpi rank: ", mpi_rank
       call fatal_error(h//" Error h5sclose_f")
    endif

    call h5pclose_f(plist_id, error)
    if (error .ne. 0) then
       print*, "error h5pclose_f - mpi rank: ", mpi_rank
       call fatal_error(h//" Error h5pclose_f")
    endif

  end subroutine hdf_parallel_output

end module hdf5_parallel_engine
