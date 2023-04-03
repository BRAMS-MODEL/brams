module ModPostTypes

  include "files.h"
   
  public :: PostVarType
  type PostVarType
    integer :: ivar_type        ! Dimension
    character(len = 16) :: fieldName        ! post field name
    character(len = 64) :: fieldDescription ! post field description
    character(len = 16) :: fieldUnits       ! post field units
    integer             :: netcdfId
  end type PostVarType

  type(PostVarType), allocatable :: all_post_variables(:)

  type fieldID
     ! post field dimensionality (encoded)

     ! ivar_type = 2 means 2D field (just surface)
     ! ivar_type = 3 means 3D field (surface and atmospheric verticals)
     ! ivar_type = 7 means 3D field (surface and patches)
     ! ivar_type = 8 means 4D field (surface and underground verticals and patches)
     ! ivar_type =10 means unknown

     integer           :: ivar_type
     integer           :: nLevCode         ! grads code for number of verticals
     character(len=12) :: fieldName        ! post field name
     character(len=64) :: fieldDescription ! post field description
     character(len=8)  :: fieldUnits       ! post field units
     type(fieldID), pointer :: next => null()
  end type fieldID

  real, parameter :: undef=-9.99e+33
  integer, parameter :: INT_UNDEF=0


  ! PostGrid:
  !
  ! Contains all info required for post grids
  !
  ! Post grids could be 2D, 3D or 4D arrays. In any case, they
  ! have surface components (X,Y), that ocupy leading first and
  ! second dimension. In some cases, there are atmospheric vertical levels,
  ! that ocupy the third dimension. In other cases, remaining dimensions
  ! represent sub-surface levels, or wave bins, or vegetation patches.
  ! Post grid dimensionality and the meaning of dimensions is encoded
  ! and stored into a fieldID type.
  !
  ! Surface components of post grids are regular lat-lon planes,
  ! dimensioned (nLon,nLat). Post grid points longitude (latitude) are stored
  ! at lon(i), i=1,nLon (lat(i), i=1,nLat). 
  !
  ! Post Grids are built from Brams Grids.
  !
  ! Brams Grid surface components are regular on a tangent plane to 
  ! the earth's surface. They are not regular in lat-lon.
  !
  ! Mapping of Brams Grid surface to Post Grid surface varies with user
  ! selection (field "proj" at POST namelist).
  !
  ! If user selects proj="no", mapping is 1 to 1 and brams grid is assumed
  ! to be regular lat-lon (which is not). Consequently, this option introduces
  ! distortions at the pictures.
  !
  ! If user selects proj="yes", post grid comes out of interpolating brams
  ! grid points, eliminating distortions.

  type PostGrid

     ! experimet title (for ctl files)

     character(len=f_name_length) :: title

     ! grid projection
     ! post grid may be a regular lat-lon grid (if project=true), and
     ! as so it reports a subset of brams grid, or it may be a irregular
     ! lat-lon grid (if project=false), and as so it reports the entire
     ! brams grid.
     ! this option is selected at namelist file.

     logical :: project

     ! list of saved field ids

     type(fieldID), pointer :: head => null()
     type(fieldID), pointer :: tail => null()

     ! whether field id list is being built or not.
     ! this field prevents list duplication whenever
     ! a second date is processed for the same grid.
     ! it allows building the list just once for each
     ! grid and use it for all dates.
     ! 
     ! before inserting a field on the list, the insertion
     ! procedure (saveCurrentFieldID) test this field
     ! and inserts iff true.
     !
     ! dinamics of this field values:
     ! the field is set true on variable creation 
     ! (at CreatePostGrid), and remains
     ! true until the list is written at control file,
     ! when is set false (at FillGradsControlFile). 

     logical :: buildingFieldIDList

     ! temporary storage of current field id

     integer           :: ivar_type
     character(len=16) :: fieldName        ! post field name
     character(len=64) :: fieldDescription ! post field description
     character(len=16) :: fieldUnits       ! post field units

     ! grads binary file for fields

     character(len=f_name_length) :: binFileName
     integer           :: unitBinFile
     integer           :: lastRec

     ! grads ascii file for grads control

     character(len=f_name_length) :: ctlFileName
     integer           :: unitCtlFile

     ! grads file options

     logical :: oneToOne     ! true iff one grads file for each date and grid
     logical :: template     ! if a template file should be built

     ! grads file details

     character(len=15) :: chstep    ! delta t among post outputs

     ! post grid surface in lat-lon

     real              :: firstLon  ! first longitude (degrees)
     real              :: delLon    ! delta longitude
     integer           :: nLon      ! number of longitudes
     real, allocatable :: lon(:)    ! longitudes

     real              :: firstLat  ! first latitude (degrees)
     real              :: delLat    ! delta latitude
     integer           :: nLat      ! number of latitudes
     real, allocatable :: lat(:)    ! latitudes

     ! map BRAMS grid surface onto post grid surface.
     ! In any case, post grid surface has size (nLon,nLat)
     ! and BRAMS grid surface has size (nnxp, nnyp)

     ! if no mapping is applied, post grid surface is
     ! the subset of BRAMS grid surface starting at BRAMS indices
     ! (xStart,yStart) with size (nLon,nLat)

     integer              :: xStart    ! first BRAMS x index
     integer              :: yStart    ! first BRAMS y index

     ! if mapping is applied, post grid surface stills has size
     ! (nLon, nLat), but does not start at BRAMS (xStart,yStart). 
     ! Instead, point (ix,iy) of the post grid is interpolated from BRAMS surface points
     ! (xmap(ix,iy),ymap(ix,iy)) + (0,0), + (0,1), + (1,0), + (1,1)
     ! with weights weight (ix,iy,1), (ix,iy,2), (ix,iy,3), (ix,iy,4)
     ! whenever the four points fall into BRAMS surface grid,
     ! or the field receives "undef" at this point and interpolation is made
     ! from (1,1).

     integer, allocatable :: xMap(:,:) ! BRAMS x index, indexed by post grid indices
     integer, allocatable :: yMap(:,:) ! BRAMS y index, indexed by post grid indices
     real,    allocatable :: weight(:,:,:) ! weight of four neighbour BRAMS points

     ! domain decomposition of post xy grid:
     ! the post xy grid does not have its own domain decomposition.
     ! instead, it uses BRAMS domain decomposition.
     ! each process computes its portion of the post xy grid that it
     ! owns on BRAMS domain decomposition.
     ! computed portion of the xy grid is packed into a 1D
     ! array that is sent to a single process; the sent arrays
     ! are gathered at a single array by MPI_Gatherv. Offset of first position of each
     ! sent array in the gathered array is defined by disp(proc), indexed by BRAMS process number.
     ! pack size at each process is localSize(:), indexed by BRAMS process number;
     ! pairs (X,Y) of local BRAMS indices are defined by (packXLocal(i),packYLocal(i))
     ! with i=1,localSizeThisProc. Observe that the packing takes only post grid points; 
     ! consequently, ghost zones are eliminated whenever unused at post grid.
     ! the packed array is unpacked at the receiving process at time of writing;
     ! the order of writing is unpackMap(i), i=1,nLat*nLon

     ! packWeight is the portion of weight that belongs to this process;
     ! it is indexed 1:localSizeThisProc, to simplify interpolation.

     integer, allocatable :: localSize(:)
     integer, allocatable :: disp(:)
     integer              :: localSizeThisProc
     integer, allocatable :: packXLocal(:)
     integer, allocatable :: packYLocal(:)
     real,    allocatable :: packWeight(:,:)

     integer, allocatable :: unpackMap(:)

     ! Summary of surface XY grids:
     ! There are three XY grids and two lat-lon. The XY grids are:
     ! (1) BRAMS global XY grid (not domain decomposed), dimensioned (1:nnxp, 1:nnyp);
     ! (2) BRAMS local XY grid (domain decomposed), dimensioned (1:mxp, 1:myp);
     ! (3) POST global XY grid (not domain decomposed), dimensioned (1:nLon, 1:nLat)
     !
     ! Mappings among grids:
     ! (2) from (1):
     !    BRAMS local (i,j) = BRAMS global (i+nodei0(proc,grid), j+nodej0(proc,grid))
     ! (1) from (2):
     !    BRAMS global(i,j) = BRAMS local  (i-nodei0(proc,grid), j-nodej0(proc,grid))
     ! (3) from (1):
     ! if no mapping (project==.false.)
     !    POST global(i,j)  = BRAMS global (i+xStart-1, j+yStart-1)
     ! if mapping (project==.true.)
     !    POST global(i,j)  = BRAMS global (xMap(i,j), yMap(i,j))
     ! (3) from (2):
     ! if no mapping (project==.false.)
     !    POST global(i,j)  = BRAMS local (i+xStart-1-nodei0, j+yStart-1-nodej0)
     ! if mapping (project==.true.)
     !    POST global(i,j)  = BRAMS local (xMap(i,j)-nodei0, yMap(i,j)-nodej0)
     !
     ! POST grid domain decomposition:
     ! Post grid values that correspond to BRAMS local values are computed and 
     ! stored at array localChunk by procedure OutputGradsField. This array has
     ! size localSizeThisProc. Its values are obtained from BRAMS local grid by
     ! localChunk(i) = BRAMS local (packXLocal(i),packYLocal(i))
     !
     ! The set of localChunk of all processes are gathered by MPI_GATHER at array
     ! gathered, dimensioned nLon*nLat, by procedure OutputGradFields. At this array,
     ! the localChunk of proc p starts at index disp(p)+1 and has size localSize(p).
     ! The gathered array is not ready for output, since element order is not the 
     ! required order. To output, use gathered(unpackMap(:)).

     ! verticals
     ! post will compute nVert vertical levels
     ! nVert is obtained from namelist values 
     ! zlevmax and inplevs;

     integer :: nVert ! # of output vertical levels 

     ! vertical scale unit is encoded into vertScaleCode:
     !  0 means BRAMS levels
     !  1 means pressure
     !  2 means height
     !  3 means user selected BRAMS levels
     ! vertScaleCode is obtained from namelist value ipresslev

     integer :: vertScaleCode ! encoded vertical scale 

     ! vertical scale values are stored at vertScaleValues,
     ! dimensioned nVert

     real, allocatable :: vertScaleValues(:)

     ! remaining fields could be used or not, depending on the
     ! selection of vertical scale unit (vertScaleCode):

     ! when vertScaleCode is 0 or 3 (that is, use original BRAMS verticals)
     ! then zMap maps Post verticals to BRAMS verticals, that is
     ! zMap(i), i=1,nVert, is the index of BRAMS field vertical

     ! when vertScaleCode is 1 or 2, zMap is undefined

     integer, allocatable :: zMap(:) 

     ! pi contains pressure, only when
     ! posting pressure verticals (vertScaleCode=1)

     real, allocatable :: pi(:,:,:)

     ! topo contains topography only when
     ! posting pressure verticals (vertScaleCode=1) or
     ! posting height verticals (vertScaleCode=2)

     real, allocatable :: topo(:,:)
  end type PostGrid

  contains

   function getPostVarible(varName) result(one_post_variable)
      character(len = *), intent(in) :: varName
      type(PostVarType) :: one_post_variable
      integer :: i

      do i = 1, size(all_post_variables)
         if(varName .eq. all_post_variables(i)%fieldName) then
            one_post_variable = all_post_variables(i)
            return
         end if
      end do
      one_post_variable%fieldName = ''
   end function getPostVarible

end module ModPostTypes
