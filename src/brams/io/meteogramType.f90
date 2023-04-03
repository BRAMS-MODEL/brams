module meteogramType


!!$ type ModelPixelPointer
!!$  real,              dimension(:), pointer :: cellPointer !$ pointer to a cell or vector of cells in vtables. (only real)
!!$  type(ModelPixelPointer),         pointer :: prev
!!$ end type ModelPixelPointer
!!$
!!$ type MeteoVars
!!$  character(len=50)                        :: name        !$ variable output identification.
!!$  real                                     :: localSum    !$ sum over all pixels that includes current polygon .
!!$  real                                     :: localMin    !$ min value of all pixels that includes current polygon.
!!$  real                                     :: localMax    !$ max value of all pixels that includes current polygon.
!!$  character(len=50), dimension(:), pointer :: compName    
!!$  integer                                  :: nVars       !$ how many variable are necessary to create requested one.
!!$  type(ModelPixelPointer), pointer         :: vTablesVars !$ pixel list to vtables.
!!$ end type MeteoVars


 type ModelPixelPointer
  real, pointer :: cellPointer !$ pointer to a cell or vector of cells in vtables. (only real)
 end type ModelPixelPointer

 type PolygonVertices
  character(len=250),      dimension(:), allocatable :: vTableName
  type(ModelPixelPointer), dimension(:), allocatable :: var
 end type PolygonVertices

 type PolygonContainer
  character(len=250)                            :: name      !$ polygon identification.
  real                                          :: lat
  real                                          :: lon
  integer                                       :: grid
  integer                                       :: totalVertices
  type(PolygonVertices),  pointer, dimension(:) :: vertices  !$ meteogram variables.
  type(PolygonContainer), pointer               :: prev
 end type PolygonContainer

 integer :: meteogramOutUnit !$ local variable.
 integer :: nCities          !$ local variable.
 
 end module meteogramType
