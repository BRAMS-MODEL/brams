module satPolyColision

 implicit none

 type tPoint
  real :: x
  real :: y
 end type tPoint

 type tPolygon
  character(50)                           :: name
  real                                    :: lat
  real                                    :: lon
  type(tPoint), allocatable, dimension(:) :: points
 end type tPolygon

 contains

 !#
 !# Collision check by using Separating Axis Theorem
 !#
 function check2dConvexPolyCollision(poly0, poly1) result (ret)
  type(tPolygon), intent(in) :: poly0
  type(tPolygon), intent(in) :: poly1
  logical                    :: ret

  integer            :: idx
  real, dimension(2) :: mmPoly0
  real, dimension(2) :: mmPoly1
  type(tPoint)       :: axis

   ret = .true.

   !@checking for polygon 1
   do idx = 1, size(poly0%points)-1
    axis = Get2dNormal(poly0%points(idx), poly0%points(idx+1), .true.)
    mmPoly0 = Project1d(axis, poly0%points)
    mmPoly1 = Project1d(axis, poly1%points)
    if( ((mmPoly0(1) - mmPoly1(2)) .gt. 0.0) .or. ((mmPoly1(1) - mmPoly0(2)) .gt. 0.0) )then
     ret = .false.
     return
    end if
   end do ! # do idx = 1, size(poly0%point)-1

   !@ polygon 2
   do idx = 1, size(poly1%points)-1
    axis = Get2dNormal(poly1%points(idx), poly1%points(idx+1), .true.)
    mmPoly0 = Project1d(axis, poly0%points)
    mmPoly1 = Project1d(axis, poly1%points)

    if( ((mmPoly0(1) - mmPoly1(2)) .gt. 0.0) .or. ((mmPoly1(1) - mmPoly0(2)) .gt. 0.0) )then
     ret = .false.
     return
    end if

   end do ! # do idx = 1, size(poly1%point)-1


 end function check2dConvexPolyCollision

 !#
 !# Normalize vector size.
 !#
 subroutine Normalize(p0)
  type(tPoint), intent(inout) :: p0

  real :: len

   len  = sqrt((p0%x*p0%x)+(p0%y*p0%y))
   p0%x = real(p0%x / len)
   p0%y = real(p0%y / len)

 end subroutine Normalize

 !#
 !# returns normalized vector.
 !#
 function Get2dNormal(p0, p1,opt) result(pN)
  type(tPoint), intent(in) :: p0
  type(tPoint), intent(in) :: p1
  logical, intent(in)      :: opt
  type(tPoint)             :: pN

   !# right or left normals.
   if(opt)then
    pN%x = -(p1%y - p0%y)
    pN%y =  (p1%x - p0%x)
   else
    pN%x =  (p1%y - p0%y)
    pN%y = -(p1%x - p0%x)
   end if

   call Normalize(pN)

 end function Get2dNormal

 !#
 !# point dot product.
 !#
 function dot(p0, p1) result (r)
  type(tPoint), intent(in) :: p0
  type(tPoint), intent(in) :: p1

  real :: r

   r = (p0%x * p1%x) + (p0%y * p1%y)

 end function dot

 !#
 !# Project polygons vertices into a normalized vector.
 !#
 function Project1d(axis, verts) result(mm)
  type(tPoint), intent(in)               :: axis
  type(tPoint), dimension(:), intent(in) :: verts

  real, dimension(2)                     :: mm

  integer :: idx
  real    :: curr
  real    :: min
  real    :: max

   curr = dot(verts(1), axis)
   mm   = curr

   do idx = 2, size(verts)-1
    curr = dot(verts(idx), axis)
    if(curr .lt. mm(1)) mm(1) = curr
    if(curr .gt. mm(2)) mm(2) = curr
   end do

 end function Project1d

end module satPolyColision
