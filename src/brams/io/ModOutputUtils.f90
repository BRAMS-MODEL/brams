module ModOutputUtils
  use var_tables, only: &
       var_tables_r, &
       GetVTabEntry

  use mem_basic, only: basic_g

  use mem_turb, only: &
       turb_g, &
       idiffk, &
       xkhkm

  use dump

  implicit none

  private
  public :: GetVarFromMemToOutput

  interface GetVarFromMemToOutput
     module procedure GetVarFromMemToOutput_2D
     module procedure GetVarFromMemToOutput_3D
     module procedure GetVarFromMemToOutput_4D
  end interface

  character(len=*),parameter :: sourceName='ModOutputUtils.f90' !Name of source code
  include "constants.f90"
contains


  subroutine GetVarFromMemToOutput_2D (varName, ngrd, arrayOut)
    character(LEN=*),   intent(in ) :: varName
    integer,            intent(in ) :: ngrd
    real,	        intent(out) :: arrayOut(:,:)

    type(var_tables_r), pointer   :: vtabPtr
    character(len=30), parameter :: h="**(GetVarFromMemToOutput_2D)**"


    ! get vtab_r entry that points to the field; stop if not there

    vtabPtr => null()
    call GetVTabEntry(varName, ngrd, vtabPtr)
    if (.not. associated(vtabPtr)) then
        iErrNumber=dumpMessage(c_tty,c_yes,sourceName,h &
              ,c_fatal,'var '//trim(varName)//' not found in vtab_r')
    end if

    ! copy the field

    arrayOut = vtabPtr%var_p_2D
  end subroutine GetVarFromMemToOutput_2D





  subroutine GetVarFromMemToOutput_3D (varName, ngrd, arrayOut)
    character(LEN=*),   intent(in ) :: varName
    integer,            intent(in ) :: ngrd
    real,	        intent(out) :: arrayOut(:,:,:)

    type(var_tables_r), pointer   :: vtabPtr
    real :: transposed(size(arrayOut,3),size(arrayOut,1),size(arrayOut,2))
    character(len=len(varName)) :: varnIn, varnOut
    character(len=8) :: c0, c1
    character(len=30), parameter :: h="**(GetVarFromMemToOutput_3D)**"
    integer :: i, j, k
    integer :: n1, n2, n3

    ! arrayOut has index order prepared for output, that is,
    ! horizontal planes for all verticals, patches or waves

    n1=size(arrayOut,1) ! x points
    n2=size(arrayOut,2) ! y points
    n3=size(arrayOut,3) ! z (verticals, patches or waves)

    ! output fields PI and HKH do not exist at vartables
    ! and have to be computed from fields PP and HKM
    ! get required fields from vartables

    if (trim(varName) == 'PI') then
       varnIn = 'PP'
    else if (trim(varName) == 'HKH') then
       varnIn = 'HKM'
    else
       varnIn = varName
    end if

    vtabPtr => null()
    call GetVTabEntry(varnIn, ngrd, vtabPtr)
    if (.not. associated(vtabPtr)) then
               iErrNumber=dumpMessage(c_tty,c_yes,sourceName,h &
              ,c_fatal,'var '//trim(varName)//' not found in vtab_r')
    end if

    ! copy the field observing that
    ! idim_type==3 has vtab indices (z,x,y), same as "transposed"
    ! idim_type==6 has vtab indices (x,y,patch), while "transposed" has (patch,x,y)
    ! idim_type==7 has vtab indices (x,y,wave), while "transposed" has (wave,x,y)
    ! while arrayOut is always (x, y, other dimension)
    !
    ! to allow ghost zone update, let transposed field as
    ! (other dimension, x, y)

    if (vtabPtr%idim_type == 3) then
       transposed = vtabPtr%var_p_3d
    else if (vtabPtr%idim_type == 6 .or. vtabPtr%idim_type == 7) then
       do j = 1, n2
          do i = 1, n1
             do k = 1, n3
                transposed(k,i,j) = vtabPtr%var_p_3d(i,j,k)
             end do
          end do
       end do
    else
       write(c0,"(i8)") vtabPtr%idim_type
       call fatal_error(h//" not prepared for idim_type=="//trim(adjustl(c0)))
    end if

    ! convert fields PP, HKM and VKH from vartables
    ! store field at arrayOut
    ! rearranging arrayOut from (k,i,j) to (i,j,k)


    select case (trim(vtabPtr%name))
    case ('PP')

       ! Output total Exner function

       do j = 1, n2
          do i = 1, n1
             do k = 1, n3
                arrayOut(i,j,k) = transposed(k,i,j) + &
                     basic_g(ngrd)%pi0(k,i,j)
             end do
          end do
       end do

    case ('HKM')

       ! Convert to HKM to HKH (note that VKH is HKH for Deardorff)
       ! and transpose

       if (idiffk(ngrd) <= 3) then

          do j = 1, n2
             do i = 1, n1
                do k = 1, n3
                   arrayOut(i,j,k) = transposed(k,i,j) * &
                        xkhkm(ngrd) / basic_g(ngrd)%dn0(k,i,j)
                end do
             end do
          end do

       else if (idiffk(ngrd) >= 4) then

          do j = 1, n2
             do i = 1, n1
                do k = 1, n3
                   arrayOut(i,j,k) = turb_g(ngrd)%vkh(k,i,j) / &
                        basic_g(ngrd)%dn0(k,i,j)
                end do
             end do
          end do
       endif

    case ('VKH')

       ! Un-density weight VKH and transpose

       do j = 1, n2
          do i = 1, n1
             do k = 1, n3
                arrayOut(i,j,k) = transposed(k,i,j) / &
                     basic_g(ngrd)%dn0(k,i,j)
             end do
          end do
       end do

    case default

       ! move verticals (or waves or patches) from first to third dimension

       do j = 1, n2
          do i = 1, n1
             do k = 1, n3
                arrayOut(i,j,k) = transposed(k,i,j)
             end do
          end do
       end do
    end select
  end subroutine GetVarFromMemToOutput_3D





  subroutine GetVarFromMemToOutput_4D (varName, ngrd, arrayOut)
    character(LEN=*),   intent(in ) :: varName
    integer,            intent(in ) :: ngrd
    real,	        intent(out) :: arrayOut(:,:,:,:)

    type(var_tables_r), pointer   :: vtabPtr
    integer :: i, j, k, l
    integer :: n1, n2, n3, n4
    real :: transposed(size(arrayOut,3),size(arrayOut,4),size(arrayOut,1),size(arrayOut,2))
    character(len=8) :: c0
    character(len=30), parameter :: h="**(GetVarFromMemToOutput_4D)**"

    ! arrayOut has index order prepared for output, that is,
    ! horizontal planes for all verticals, patches or waves

    n1=size(arrayOut,1) ! x points
    n2=size(arrayOut,2) ! y points
    n3=size(arrayOut,3) ! z (verticals)
    n4=size(arrayOut,4) ! patches or waves

    ! get vtab_r entry that points to the field; stop if not there

    vtabPtr => null()
    call GetVTabEntry(varName, ngrd, vtabPtr)
    if (.not. associated(vtabPtr)) then
               iErrNumber=dumpMessage(c_tty,c_yes,sourceName,h &
              ,c_fatal,'var '//trim(varName)//' not found in vtab_r')
    end if

    ! prepared for idim_type 4 or 5

    if (vtabPtr%idim_type /= 4 .and. vtabPtr%idim_type /= 5) then
       write(c0,"(i8)") vtabPtr%idim_type
       call fatal_error(h//" not prepared for idim_type="//trim(adjustl(c0)))
    end if

    ! copy the field observing that
    ! idim_type==4 or 5  has vtab indices (some_z,x,y,patch),
    ! while "transposed" has indices (some_z*patch, x, y) to
    ! simplify border updates


    do l = 1, n4
       do j = 1, n2
          do i = 1, n1
             do k = 1, n3
                transposed(k,l,i,j) = vtabPtr%var_p_4D(k,i,j,l)
             end do
          end do
       end do
    end do

    ! set verticals, patches and waves at desired indices

    do l = 1, n4
       do k = 1, n3
          do j = 1, n2
             do i = 1, n1
                arrayOut(i,j,k,l) = transposed(k,l,i,j)
             end do
          end do
       end do
    end do
  end subroutine GetVarFromMemToOutput_4D
end module ModOutputUtils
