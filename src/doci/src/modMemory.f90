module modMemory
  implicit none

  integer,allocatable :: nFilesPerProc(:)
  integer,allocatable :: ini(:)
  integer,allocatable :: fim(:)
  integer ( kind = 4 ) id
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) p
  real ( kind = 8 ) wtime
  character(len=256) :: prefix
  character(len=256) :: outFolder
  integer :: imonth1
  integer :: idate1
  integer :: iyear1
  integer :: itime1
  integer :: ntimes
  integer :: tincrem
  character(len=4) :: source
  integer :: lastTime

end module modMemory