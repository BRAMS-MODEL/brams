module ModBuffering

  use ModParallelEnvironment, only: &
       MsgDump


  ! ModBuffering: procedures to copy data 
  ! from/to array sections to/from 1D buffer

  implicit none
  private
  public :: FieldSection2Buffer
  public :: Buffer2FieldSection

  interface FieldSection2Buffer
     module procedure FieldSection2Buffer_2D
     module procedure FieldSection2Buffer_3D
     module procedure FieldSection2Buffer_4D
  end interface

  interface Buffer2FieldSection
     module procedure Buffer2FieldSection_2D
     module procedure Buffer2FieldSection_3D
     module procedure Buffer2FieldSection_4D
  end interface

  logical, parameter :: dumpLocal = .false.


contains


  ! FieldSection2Buffer_2D: copy field(iStart:iEnd,jStart:jEnd) into
  !                       buffer(k), starting at k=lastBuffer+1, in
  !                       array element order


  subroutine FieldSection2Buffer_2D(field, idim_type, &
       iStart, iEnd, jStart, jEnd, &
       buffer, lastBuffer)
    real, intent(in) :: field(:,:)
    integer, intent(in) :: idim_type
    integer, intent(in) :: iStart
    integer, intent(in) :: iEnd
    integer, intent(in) :: jStart
    integer, intent(in) :: jEnd
    real, intent(inout) :: buffer(:)
    integer, intent(inout) :: lastBuffer

    integer :: i, j, firstBuffer
    character(len=8) :: c0, c1, c2, c3, c4, c5
    character(len=*), parameter :: h="**(FieldSection2Buffer_2D)**"

    firstBuffer = lastBuffer + 1
    if(size(field,1)<iStart .or. size(field,2)<jend) print *,'fudeu 2d'

    select case (idim_type)
    case(2)
       do j=jStart, jEnd
          do i=iStart, iEnd
             lastBuffer=lastBuffer+1
             buffer(lastBuffer)=field(i,j)
          end do
       end do
       if (dumpLocal) then
          write(c0,"(i8)") firstBuffer
          write(c1,"(i8)") lastBuffer
          write(c2,"(i8)") iStart
          write(c3,"(i8)") iEnd
          write(c4,"(i8)") jStart
          write(c5,"(i8)") jEnd
          call MsgDump(h//" buffer"//&
               "["//trim(adjustl(c0))//":"//trim(adjustl(c1))//"]"//&
               "<-field"//&
               "["//trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
               trim(adjustl(c4))//":"//trim(adjustl(c5))//"]")
       end if
    case default
       write(c0,"(i8)") idim_type
       call fatal_error(h//" wrong idim_type="//trim(adjustl(c0)))
    end select
  end subroutine FieldSection2Buffer_2D


  ! FieldSection2Buffer_3D: according to idim_type, copy 
  !                       field(:,iStart:iEnd,jStart:jEnd) or
  !                       field(iStart:iEnd,jStart:jEnd,:) into
  !                       buffer(k), starting at k=lastBuffer+1, in
  !                       array element order


  subroutine FieldSection2Buffer_3D(field, idim_type,&
       iStart, iEnd, jStart, jEnd, &
       buffer, lastBuffer)
    real, intent(in) :: field(:,:,:)
    integer, intent(in) :: idim_type
    integer, intent(in) :: iStart
    integer, intent(in) :: iEnd
    integer, intent(in) :: jStart
    integer, intent(in) :: jEnd
    real, intent(inout) :: buffer(:)
    integer, intent(inout) :: lastBuffer

    integer :: i, j, k, kMax, firstBuffer
    character(len=8) :: c0, c1, c2, c3, c4, c5, c6
    character(len=*), parameter :: h="**(FieldSection2Buffer_3D)**"

    firstBuffer = lastBuffer + 1

    select case (idim_type)

    case(3)
       kMax = size(field,1)
       do j=jStart, jEnd
          do i=iStart, iEnd
             do k = 1, kMax
                lastBuffer=lastBuffer+1
                buffer(lastBuffer)=field(k,i,j)
             end do
          end do
       end do
       if (dumpLocal) then
          write(c0,"(i8)") firstBuffer
          write(c1,"(i8)") lastBuffer
          write(c2,"(i8)") iStart
          write(c3,"(i8)") iEnd
          write(c4,"(i8)") jStart
          write(c5,"(i8)") jEnd
          write(c6,"(i8)") kMax
          call MsgDump(h//" buffer"//&
               "["//trim(adjustl(c0))//":"//trim(adjustl(c1))//"]"//&
               "<-field"//&
               "[1:"//trim(adjustl(c6))//","//&
               trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
               trim(adjustl(c4))//":"//trim(adjustl(c5))//"]")
       end if

    case(6:7)
       do k = 1, size(field,3)
          do j=jStart, jEnd
             do i=iStart, iEnd
                lastBuffer=lastBuffer+1
                buffer(lastBuffer)=field(i,j,k)
             end do
          end do
       end do
       if (dumpLocal) then
          write(c0,"(i8)") firstBuffer
          write(c1,"(i8)") lastBuffer
          write(c2,"(i8)") iStart
          write(c3,"(i8)") iEnd
          write(c4,"(i8)") jStart
          write(c5,"(i8)") jEnd
          write(c6,"(i8)") kMax
          call MsgDump(h//" buffer"//&
               "["//trim(adjustl(c0))//":"//trim(adjustl(c1))//"]"//&
               "<-field"//&
               "["//trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
               trim(adjustl(c4))//":"//trim(adjustl(c5))//","//&
               "1:"//trim(adjustl(c6))//"]")
       end if
       
    case default
       write(c0,"(i8)") idim_type
       call fatal_error(h//" wrong idim_type="//trim(adjustl(c0)))
    end select
  end subroutine FieldSection2Buffer_3D


  ! FieldSection2Buffer_4D: according to idim_type, copy 
  !                       field(:,iStart:iEnd,jStart:jEnd,:) into
  !                       buffer(k), starting at k=lastBuffer+1, in
  !                       array element order


  subroutine FieldSection2Buffer_4D(field, idim_type, &
       iStart, iEnd, jStart, jEnd, &
       buffer, lastBuffer)
    real, intent(in) :: field(:,:,:,:)
    integer, intent(in) :: idim_type
    integer, intent(in) :: iStart
    integer, intent(in) :: iEnd
    integer, intent(in) :: jStart
    integer, intent(in) :: jEnd
    real, intent(inout) :: buffer(:)
    integer, intent(inout) :: lastBuffer

    integer :: i, j, k, l, kMax, lMax, firstBuffer
    character(len=8) :: c0, c1, c2, c3, c4, c5, c6, c7
    character(len=*), parameter :: h="**(FieldSection2Buffer_4D)**"

    firstBuffer = lastBuffer + 1

    select case (idim_type)

    case(4:5)
       kMax = size(field,1)
       lMax = size(field,4)
       do l = 1, lMax
          do j = jStart, jEnd
             do i = iStart, iEnd
                do k = 1, kMax
                   lastBuffer=lastBuffer+1
                   buffer(lastBuffer)=field(k,i,j,l)
                end do
             end do
          end do
       end do
       if (dumpLocal) then
          write(c0,"(i8)") firstBuffer
          write(c1,"(i8)") lastBuffer
          write(c2,"(i8)") iStart
          write(c3,"(i8)") iEnd
          write(c4,"(i8)") jStart
          write(c5,"(i8)") jEnd
          write(c6,"(i8)") kMax
          write(c7,"(i8)") lMax
          call MsgDump(h//" buffer"//&
               "["//trim(adjustl(c0))//":"//trim(adjustl(c1))//"]"//&
               "<-field"//&
               "[1:"//trim(adjustl(c6))//","//&
               trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
               trim(adjustl(c4))//":"//trim(adjustl(c5))//","//&
               "1:"//trim(adjustl(c7))//"]")
       end if

    case default
       write(c0,"(i8)") idim_type
       call fatal_error(h//" wrong idim_type="//trim(adjustl(c0)))
    end select
  end subroutine FieldSection2Buffer_4D


  ! Buffer2FieldSection_2D: copy buffer(k), starting at k=lastBuffer+1, into
  !                     field(iStart:iEnd,jStart:jEnd) in array element order


  subroutine Buffer2FieldSection_2D(field, idim_type, &
       iStart, iEnd, jStart, jEnd, &
       buffer, lastBuffer)
!    use CUPARM_GRELL3, only: g3d_g

    real, intent(inout) :: field(:,:)
    integer, intent(in) :: idim_type
    integer, intent(in) :: iStart
    integer, intent(in) :: iEnd
    integer, intent(in) :: jStart
    integer, intent(in) :: jEnd
    real, intent(in) :: buffer(:)
    integer, intent(inout) :: lastBuffer

    integer :: i, j, firstBuffer
    character(len=8) :: c0, c1, c2, c3, c4, c5, c6, c7
    character(len=*), parameter :: h="**(Buffer2FieldSection_2D)**"

    firstBuffer = lastBuffer+1

    select case (idim_type)
    case(2)
!print *,'LFR-DBG: inside BUF 1: ',size(g3d_g(1)%cugd_ttens,1), &
!                 size(g3d_g(1)%cugd_ttens,2),size(g3d_g(1)%cugd_ttens,3); call flush(6)
     ! idim_type codes dimensionality of var_p:
     ! idim_type == 2 means (nmxp, nmyp)
     ! idim_type == 3 means (nmzp, nmxp, nmyp)
     ! idim_type == 4 means (nzg, nmxp, nmyp, npatch)
     ! idim_type == 5 means (nzs, nmxp, nmyp, npatch)
     ! idim_type == 6 means (nmxp, nmyp, npatch)
     ! idim_type == 7 means (nmxp, nmyp, nwave)


       do j=jStart, jEnd
          do i=iStart, iEnd
             lastBuffer=lastBuffer+1
             field(i,j)=buffer(lastBuffer)
          end do
       end do
    
       if (dumpLocal) then
          write(c6,"(i8)") size(field,1)
          write(c7,"(i8)") size(field,2)
          write(c0,"(i8)") firstBuffer
          write(c1,"(i8)") lastBuffer
          write(c2,"(i8)") iStart
          write(c3,"(i8)") iEnd
          write(c4,"(i8)") jStart
          write(c5,"(i8)") jEnd
          call MsgDump(h//" field"//&
               "["//trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
               trim(adjustl(c4))//":"//trim(adjustl(c5))//"]"//&
               "<-buffer"//&
               "["//trim(adjustl(c0))//":"//trim(adjustl(c1))//"]"//"("//c6//","//c7//")")
       end if
! print *,'LFR-DBG: inside BUF 2: ',size(g3d_g(1)%cugd_ttens,1), &
 !                size(g3d_g(1)%cugd_ttens,2),size(g3d_g(1)%cugd_ttens,3); call flush(6)  
!LFR-DBG BEGIN
!LFR if(size(g3d_g(1)%cugd_ttens,3)==0) then
!LFR           write(c6,"(i8)") size(field,1)
!LFR           write(c7,"(i8)") size(field,2)
!LFR           write(c0,"(i8)") firstBuffer
!LFR           write(c1,"(i8)") lastBuffer
!LFR           write(c2,"(i8)") iStart
!LFR           write(c3,"(i8)") iEnd
!LFR           write(c4,"(i8)") jStart
!LFR           write(c5,"(i8)") jEnd
!LFR           print *,h//" field"//&
!LFR                "["//trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
!LFR                trim(adjustl(c4))//":"//trim(adjustl(c5))//"]"//&
!LFR                "<-buffer"//&
!LFR                "["//trim(adjustl(c0))//":"//trim(adjustl(c1))//"]"//"("//c6//","//c7//")"
!LFR end if
!LFR !LFR-DBG END    

    case default
       write(c0,"(i8)") idim_type
       call fatal_error(h//" wrong idim_type="//trim(adjustl(c0)))
    end select
  end subroutine Buffer2FieldSection_2D


  ! Buffer2FieldSection_3D: copy buffer(k), starting at k=lastBuffer+1, into
  !                     field(iStart:iEnd,jStart:jEnd) in array element order


  subroutine Buffer2FieldSection_3D(field, idim_type, &
       iStart, iEnd, jStart, jEnd, &
       buffer, lastBuffer)
    real, intent(inout) :: field(:,:,:)
    integer, intent(in) :: idim_type
    integer, intent(in) :: iStart
    integer, intent(in) :: iEnd
    integer, intent(in) :: jStart
    integer, intent(in) :: jEnd
    real, intent(in) :: buffer(:)
    integer, intent(inout) :: lastBuffer

    integer :: i, j, k, kMax, firstBuffer
    character(len=8) :: c0, c1, c2, c3, c4, c5, c6
    character(len=*), parameter :: h="**(Buffer2FieldSection_3D)**"

    firstBuffer = lastBuffer + 1

    select case (idim_type)

    case(3)

       kMax = size(field,1)
       do j=jStart, jEnd
          do i=iStart, iEnd
             do k = 1, kMax
                lastBuffer=lastBuffer+1
                field(k,i,j)=buffer(lastBuffer)
             end do
          end do
       end do
       if (dumpLocal) then
          write(c0,"(i8)") firstBuffer
          write(c1,"(i8)") lastBuffer
          write(c2,"(i8)") iStart
          write(c3,"(i8)") iEnd
          write(c4,"(i8)") jStart
          write(c5,"(i8)") jEnd
          write(c6,"(i8)") kMax
          call MsgDump(h//" field"//&
               "[1:"//trim(adjustl(c6))//","//&
               trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
               trim(adjustl(c4))//":"//trim(adjustl(c5))//"]"//&
               "<-buffer"//&
               "["//trim(adjustl(c0))//":"//trim(adjustl(c1))//"]")
       end if

    case(6:7)
!LFR          if(size(field,1)/=nodemxp(mynum,1) .or. size(field,2)/=nodemyp(mynum,1) &
!LFR             .or. size(field,3)/=npatch) print *, 'Fudeu 3d - 6'; call flush(6)
!LFR          if(size(field,1)/=nodemxp(mynum,1) .or. size(field,2)/=nodemyp(mynum,1) &
!LFR             .or. size(field,3)/=nwave) print *, 'Fudeu 3d - 7'; call flush(6)
       do k = 1, size(field,3)
          do j=jStart, jEnd
             do i=iStart, iEnd
                lastBuffer=lastBuffer+1
                field(i,j,k)=buffer(lastBuffer)
             end do
          end do
       end do
       if (dumpLocal) then
          write(c0,"(i8)") firstBuffer
          write(c1,"(i8)") lastBuffer
          write(c2,"(i8)") iStart
          write(c3,"(i8)") iEnd
          write(c4,"(i8)") jStart
          write(c5,"(i8)") jEnd
          write(c6,"(i8)") kMax
          call MsgDump(h//" field["//&
               trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
               trim(adjustl(c4))//":"//trim(adjustl(c5))//","//&
               "1:"//trim(adjustl(c6))//"]"//&
               "<-buffer"//&
               "["//trim(adjustl(c0))//":"//trim(adjustl(c1))//"]")
       end if
       
    case default
       write(c0,"(i8)") idim_type
       call fatal_error(h//" wrong idim_type="//trim(adjustl(c0)))
    end select
  end subroutine Buffer2FieldSection_3D


  ! Buffer2FieldSection_4D: copy buffer(k), starting at k=lastBuffer+1, into
  !                     field(:,iStart:iEnd,jStart:jEnd,:) in array element order


  subroutine Buffer2FieldSection_4D(field, idim_type, &
       iStart, iEnd, jStart, jEnd, &
       buffer, lastBuffer)
    real, intent(inout) :: field(:,:,:,:)
    integer, intent(in) :: idim_type
    integer, intent(in) :: iStart
    integer, intent(in) :: iEnd
    integer, intent(in) :: jStart
    integer, intent(in) :: jEnd
    real, intent(in) :: buffer(:)
    integer, intent(inout) :: lastBuffer

    integer :: i, j, k, l, kMax, lMax, firstBuffer
    character(len=8) :: c0, c1, c2, c3, c4, c5, c6, c7
    character(len=*), parameter :: h="**(Buffer2FieldSection_4D)**"

    firstBuffer = lastBuffer + 1

    select case (idim_type)

    case(4:5)
   
       kMax = size(field,1)
       lMax = size(field,4)
       do l = 1, lMax
          do j=jStart, jEnd
             do i=iStart, iEnd
                do k = 1, kMax
                   lastBuffer=lastBuffer+1
                   field(k,i,j,l)=buffer(lastBuffer)
                end do
             end do
          end do
       end do
       if (dumpLocal) then
          write(c0,"(i8)") firstBuffer
          write(c1,"(i8)") lastBuffer
          write(c2,"(i8)") iStart
          write(c3,"(i8)") iEnd
          write(c4,"(i8)") jStart
          write(c5,"(i8)") jEnd
          write(c6,"(i8)") kMax
          write(c7,"(i8)") lMax
          call MsgDump(h//" field"//&
               "[1:"//trim(adjustl(c6))//","//&
               trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
               trim(adjustl(c4))//":"//trim(adjustl(c5))//","//&
               "1:"//trim(adjustl(c5))//"]"//&
               "<-buffer"//&
               "["//trim(adjustl(c0))//":"//trim(adjustl(c1))//"]")
       end if

    case default
       write(c0,"(i8)") idim_type
       call fatal_error(h//" wrong idim_type="//trim(adjustl(c0)))
    end select
  end subroutine Buffer2FieldSection_4D
end module ModBuffering
