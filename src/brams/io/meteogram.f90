! RMF begin.
!####################################### 
!# This module generates meteogram outpus based on political areas of cities, or in user defined areas ** only at surface level **
!#---------------------------
!# autor : Rafael Fonseca
!# date  : 14/03/13
!# last  : 24/04/13
!#######################################

module meteogram

 use satPolyColision, only: &
 tPolygon,                  &
 check2dConvexPolyCollision

 use meteogramType
 
 use ModNamelistFile, only: namelistFile

 implicit none
 
 include "files.h"
 
 private
 
 !RAMSIN
 logical, public                      :: applyMeteogram
 real, public                         :: meteogramFreq 
 character(len=f_name_length), public :: meteogramMap  
 character(len=f_name_length), public :: meteogramDir 
 

 public :: initMeteogram
 public :: ProcessLocalMeteogram
 public :: PolygonContainer  
 public :: StoreNamelistFileAtmeteogram
 
 contains
 
 
 subroutine StoreNamelistFileAtmeteogram(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    
       	applyMeteogram = oneNamelistFile%applyMeteogram
     	meteogramFreq  = oneNamelistFile%meteogramFreq
     	meteogramMap   = oneNamelistFile%meteogramMap
	meteogramDir   = oneNamelistFile%meteogramDir
       
 end subroutine StoreNamelistFileAtmeteogram

 !# 
 !# this subroutine uses fixed values of polygons representing cities to create pixels lists. 
 !# Each proccess execute it and some can have no list of pixels to work in. 
 subroutine InitMeteogram(meteoPolys, nGrid, filename)
   use node_mod, only: &
        mchnum,        &
        master_num,    &
        ia, iz, ja, jz

   use mem_grid, only: &
        grid_g

   use var_tables, only: &
        var_tables_r,    &
        GetVTabEntry

   type(PolygonContainer), pointer :: meteoPolys
   integer, intent(in)             :: nGrid
   character(len=*), intent(in)    :: filename

   integer, external :: AvailableFileUnit

   ! pixel list 
   type t_pixelList
      integer                    :: i
      integer                    :: j
      type(t_pixelList), pointer :: prev
   end type t_pixelList
   type(t_pixelList), pointer                :: localList
   type(t_pixelList), pointer                :: auxPixel

   type(tPolygon)                            :: modelGrid
   type(tPolygon), dimension(:), allocatable :: polygons

   type(PolygonContainer), pointer           :: auxMeteoPoly

   type(ModelPixelPointer), pointer          :: auxPixelPtr

   integer                                   :: allocStat
   integer                                   :: i
   integer                                   :: j
   integer                                   :: v
   integer                                   :: np
   integer                                   :: nv
   integer                                   :: globalPts
   integer                                   :: localPts
 
   real                                      :: dx
   real                                      :: dy
   type t_vtabPtrContainer
      type(var_tables_r), pointer :: vtabPtr
   end type t_vtabPtrContainer   
   type(t_vtabPtrContainer), dimension(:), allocatable :: vtabPointers

   ! variable configuration
   integer, parameter :: nVarTables = 15
   character(len=10), parameter, dimension(nVarTables) :: varNamesInput = &
	      (/'U10MJ     ', 'V10MJ     ', 'T2MJ      ', 'RV2MJ     ', &
       'TOPT      ', 'PI0       ', 'THETA     ', 'ACCPR     ', &
       'ACCPP     ', 'ACCPS     ', 'ACCPA     ', 'ACCPG     ', &
        'ACCPH     ', 'ACONPR    ', 'PP        '/)  
				
   nCities = 0          
   ! saving variable address from vtables.
   allocate(vtabPointers(nVarTables), stat=allocStat)  
   do nv = 1, nVarTables
      call GetVTabEntry(varNamesInput(nv), nGrid, vtabPointers(nv)%vtabPtr)
      if(.not. associated(vtabPointers(nv)%vtabPtr)) then
      	if(mchnum .eq. master_num) print*, '**meteogram warning** variable:', varNamesInput(nv), 'not found!'
      endif
   end do

   ! initializing user polygons.
   meteogramOutUnit = AvailableFileUnit() + mchnum  
   call LoadPolygonFromFile(meteogramOutUnit, filename, polygons)

   ! initializing meteo data type.
   nullify(meteoPolys)
   
   allocate(modelGrid%points(5))

   do np = 1, size(polygons)
   
      nullify(localList)
      localPts = 0

      do j = ja, jz
         do i = ia, iz
            ! checking all points of this processor against current polygon.
            dx = (grid_g(nGrid)%glon(i,j) - grid_g(nGrid)%glon(i-1,j))/2.
            dy = (grid_g(nGrid)%glat(i,j) - grid_g(nGrid)%glat(i,j-1))/2.

            modelGrid%points(1)%x = grid_g(nGrid)%glon(i,j) - dx
            modelGrid%points(1)%y = grid_g(nGrid)%glat(i,j) - dy

            dx = (grid_g(nGrid)%glon(i+1,j) - grid_g(nGrid)%glon(i,j))/2.

            modelGrid%points(2)%x = grid_g(nGrid)%glon(i,j) + dx
            modelGrid%points(2)%y = grid_g(nGrid)%glat(i,j) - dy

            dy = (grid_g(nGrid)%glat(i,j+1) - grid_g(nGrid)%glat(i,j))/2.

            modelGrid%points(3)%x = grid_g(nGrid)%glon(i,j) + dx
            modelGrid%points(3)%y = grid_g(nGrid)%glat(i,j) + dy

            dx = (grid_g(nGrid)%glon(i,j) - grid_g(nGrid)%glon(i-1,j))/2.

            modelGrid%points(4)%x = grid_g(nGrid)%glon(i,j) - dx
            modelGrid%points(4)%y = grid_g(nGrid)%glat(i,j) + dy

            modelGrid%points(5)%x =  modelGrid%points(1)%x
            modelGrid%points(5)%y =  modelGrid%points(1)%y

            if(check2dConvexPolyCollision(modelGrid, polygons(np)))then
               localPts = localPts + 1
               allocate(auxPixel)
               auxPixel%i    = i
               auxPixel%j    = j
               auxPixel%prev => localList
               localList     => auxPixel      
            end if

         end do !$ do i = ia, iz
      end do  !$ do j = ja, jz  

      ! checking whether any process has points on current polygon.
      call AllReduceMeteogramNpts(localPts, globalPts) 
      
      if(globalPts .gt. 0)then
	 		nCities = nCities + 1
      
         allocate(auxMeteoPoly, stat=allocStat)
         if(localPts .gt. 0)then
	 	allocate(auxMeteoPoly%vertices(localPts), stat=allocStat)
	 else 
	 	nullify(auxMeteoPoly%vertices)
	 endif

         auxMeteoPoly%name          = polygons(np)%name  
	 auxMeteoPoly%lat           = polygons(np)%lat
	 auxMeteoPoly%lon           = polygons(np)%lon
         auxMeteoPoly%totalVertices = globalPts
         auxMeteoPoly%grid          = nGrid
         v = 1
         auxPixel => localList
	 
         do while(associated(auxPixel))
	 
            allocate(auxMeteoPoly%vertices(v)%vTableName(nVarTables), stat=allocStat)
            allocate(auxMeteoPoly%vertices(v)%var(nVarTables), stat=allocStat)

	    do nv = 1, nVarTables
               auxMeteoPoly%vertices(v)%vTableName(nv) = varNamesInput(nv)
	       
	       !rmf- very ugly workaround to allow new/old mcphys_type
	       if(.not. associated(vtabPointers(nv)%vtabPtr))then
	        nullify(auxMeteoPoly%vertices(v)%var(nv)%cellPointer)
	       else
               	select case(vtabPointers(nv)%vtabPtr%idim_type)
               	case(2)
                  auxMeteoPoly%vertices(v)%var(nv)%cellPointer => vtabPointers(nv)%vtabPtr%var_p_2D(auxPixel%i,auxPixel%j)
               	case(3)
                  auxMeteoPoly%vertices(v)%var(nv)%cellPointer => vtabPointers(nv)%vtabPtr%var_p_3D(2, auxPixel%i,auxPixel%j) ! only at surface
               	case default
                  stop 'fatal error: idim_type not allowed - meteogram.f90'
               	end select
	       end if
	       

            end do ! nv = 1, nVarTables
            v = v + 1
            auxPixel => auxPixel%prev
         end do ! do while(associated(auxPixel))   

         auxMeteoPoly%prev => meteoPolys
         meteoPolys        => auxMeteoPoly

      end if ! if(globalPts .gt. 0)
      
      !! call MeteogramBarrier()	
   end do ! np = 1, size(polygons)

   deallocate(vtabPointers)

   if(mchnum .eq. master_num)then
      meteogramOutUnit = AvailableFileUnit()
      open(meteogramOutUnit, file=trim(meteogramDir)//'-MeteogramASC.out', status='replace')
   end if

   !call  ProcessLocalMeteogram(meteoPolys)
   !stop 'RMF initMeteogram'
   
   !write ctl down.

 end subroutine initMeteogram
 

  !
 ! calculates points average.
 subroutine ProcessLocalMeteogram(meteoPolys)

   use node_mod, only: &
        mynum,              &
        mchnum,             &
        master_num,         &
        nmachs

   use mem_grid, only: &
        ztn
	
   use mem_grid,   only: &
       ngrids,          &
       time,            &
       iyear1,          &
       imonth1,         &
       idate1,          &
       itime1	
   use ModPostUtils, only: RAMS_comp_dewK_2m

   type(PolygonContainer), pointer :: meteoPolys

   ! local variables
   integer                          :: nv
   integer                          :: np
   integer                          :: vertice
   type(PolygonContainer), pointer  :: auxMeteoPoly

   real, parameter                    :: gMetar=9.8
   real, parameter                    :: RMetar=287.04
   real                               :: alturaMetar
   real                               :: tmMetar
   real                               :: aux0
   real                               :: aux1
   real                               :: aux2

   real, dimension(1,1) :: td2m, pi2d, t2mj2d


   real, dimension(:,:), allocatable    :: localSum
   real, dimension(:,:), allocatable    :: localMin
   real, dimension(:,:), allocatable    :: localMax   
   real, dimension(:,:), allocatable    :: totalSum
   real, dimension(:,:), allocatable    :: totalMin
   real, dimension(:,:), allocatable    :: totalMax  
   
   ! temporary fix
   real :: accpr, accpp, accps, accpa, accpg, accph, aconpr             
    
   !new stuff - not optimized!! RMF
   character(len=250), dimension(:), allocatable :: cityName
   integer, dimension(:), allocatable            :: nVertices
   real, dimension(:), allocatable               :: cityLat
   real, dimension(:), allocatable               :: cityLon
   integer :: city
   !new stuff - not optimized !! RMF

   integer, parameter :: nvarOutput = 6
   character(len=10), parameter, dimension(nvarOutput) :: varNamesOutput = &
	    (/'U10MJ     ', 'V10MJ     ', 'T2MJ      ', 'TD2MJ     ', &
        'SLPMETAR  ', 'PRECIP    '/)


   include "post_rconstants.h"

   if(.not. associated(meteoPolys)) return

   auxMeteoPoly => meteoPolys

   allocate(localSum(nCities, nvarOutput))
   allocate(localMin(nCities, nvarOutput))
   allocate(localMax(nCities, nvarOutput))

   allocate(totalSum(nCities, nvarOutput))
   allocate(totalMin(nCities, nvarOutput))
   allocate(totalMax(nCities, nvarOutput))
   
   allocate(nVertices(nCities))
   allocate(cityName(nCities))
   allocate(cityLat(nCities))
   allocate(cityLon(nCities))
   
   totalSum  = 0.0
   totalMin  = 0.0
   totalMax  = 0.0
   localSum  = 0.0
   localMin  =  huge(1.0)
   localMax  = -huge(1.0)
   nVertices = 1
   city      = 0
  
   do while(associated(auxMeteoPoly))
   
      city = city + 1

      nVertices(city) = auxMeteoPoly%totalVertices
      cityName(city)  = trim(auxMeteoPoly%name)
      cityLat(city)   = auxMeteoPoly%lat
      cityLon(city)   = auxMeteoPoly%lon

      !print*, trim(auxMeteoPoly%name), mchnum, size(auxMeteoPoly%vertices)
      if(associated(auxMeteoPoly%vertices))then
	
         do vertice = 1, size(auxMeteoPoly%vertices)

            do nv = 1, 3
               !U10MJ, V10MJ, T2MJ
               localSum(city, nv) = auxMeteoPoly%vertices(vertice)%var(nv)%cellPointer + localSum(city, nv)
               localMin(city, nv) = min(localMin(city, nv), auxMeteoPoly%vertices(vertice)%var(nv)%cellPointer)
               localMax(city, nv) = max(localMax(city, nv), auxMeteoPoly%vertices(vertice)%var(nv)%cellPointer)
            end do
	    !RAMS_comp_dewK_2m (Brams2Post_td2mj) - 4-RV2MJ


	    td2m(1,1)   = auxMeteoPoly%vertices(vertice)%var(4)%cellPointer * 1.e3
            pi2d(1,1)   = (auxMeteoPoly%vertices(vertice)%var(6)%cellPointer+auxMeteoPoly%vertices(vertice)%var(15)%cellPointer)
	    t2mj2d(1,1) = auxMeteoPoly%vertices(vertice)%var(3)%cellPointer

            call RAMS_comp_dewK_2m(td2m, pi2d, t2mj2d)

            localSum(city,4) = localSum(city,4) + td2m(1,1)
	    localMin(city,4) = min(localMin(city,4), td2m(1,1))
	    localMax(city,4) = max(localMax(city,4), td2m(1,1))

            !SLP METAR - 5-TOPT, 6-PI0 7-THETA 15-PP    
            ! RAMS_comp_press
            aux0 = ((auxMeteoPoly%vertices(vertice)%var(6)%cellPointer+auxMeteoPoly%vertices(vertice)%var(15)%cellPointer)&
						       /cp)**cpor*p00*.01
            ! RAMS_comp_tempk
            !aux1 = (auxMeteoPoly%vertices(vertice)%var(7)%cellPointer * aux0 /cp
            ! comp_slp_metar
            alturaMetar = ztn(2,auxMeteoPoly%grid) + auxMeteoPoly%vertices(vertice)%var(5)%cellPointer
            tmMetar     = ( 288.15 - 0.0065 * alturaMetar / 2. )

            aux1 = aux0/exp( -(gMetar*alturaMetar)/(rMetar*tmMetar) )

            localSum(city, 5) = aux1 + localSum(city, 5)
            localMin(city, 5) = min(localMin(city, 5), aux1)
            localMax(city, 5) = max(localMax(city, 5), aux1)

            !PRECIP -8) ACCPR 9)ACCPP 10)ACCPS 11)ACCPA 12)ACCPG 13)ACCPH 14)ACONPR

            !aux1 = auxMeteoPoly%vertices(vertice)%var(8)%cellPointer  + &
            !       auxMeteoPoly%vertices(vertice)%var(9)%cellPointer  + &
            !       auxMeteoPoly%vertices(vertice)%var(10)%cellPointer + &
            !       auxMeteoPoly%vertices(vertice)%var(11)%cellPointer + &
            !       auxMeteoPoly%vertices(vertice)%var(12)%cellPointer + &
            !       auxMeteoPoly%vertices(vertice)%var(13)%cellPointer + & 
            !       auxMeteoPoly%vertices(vertice)%var(14)%cellPointer 

	    
	    !rmf- very ugly workaround to allow new/old mcphys_type
	    if(associated(auxMeteoPoly%vertices(vertice)%var(8)%cellPointer))then
	     accpr = auxMeteoPoly%vertices(vertice)%var(8)%cellPointer 
	    else
	     accpr = 0.0
	    endif 

	    if(associated(auxMeteoPoly%vertices(vertice)%var(9)%cellPointer))then
	     accpp = auxMeteoPoly%vertices(vertice)%var(9)%cellPointer 
	    else
	     accpp = 0.0
	    endif 
	    if(associated(auxMeteoPoly%vertices(vertice)%var(10)%cellPointer))then
	     accps = auxMeteoPoly%vertices(vertice)%var(10)%cellPointer 
	    else
	     accps = 0.0
	    endif 
	    if(associated(auxMeteoPoly%vertices(vertice)%var(11)%cellPointer))then
	     accpa = auxMeteoPoly%vertices(vertice)%var(11)%cellPointer 
	    else
	     accpa = 0.0
	    endif 
	    if(associated(auxMeteoPoly%vertices(vertice)%var(12)%cellPointer))then
	     accpg = auxMeteoPoly%vertices(vertice)%var(12)%cellPointer 
	    else
	     accpg = 0.0
	    endif 
	    if(associated(auxMeteoPoly%vertices(vertice)%var(13)%cellPointer))then
	     accph = auxMeteoPoly%vertices(vertice)%var(13)%cellPointer 
	    else
	     accph = 0.0
	    endif 
	    if(associated(auxMeteoPoly%vertices(vertice)%var(14)%cellPointer))then
	     aconpr = auxMeteoPoly%vertices(vertice)%var(14)%cellPointer 
	    else
	     aconpr = 0.0
	    endif 	    	    	    	    	    
	    
    	    aux1 = accpr+accpp+accps+accpa+accpg+accph+aconpr

            localSum(city, 6) = aux1 + localSum(city, 6)
            localMin(city, 6) = min(localMin(city, 6), aux1)
            localMax(city, 6) = max(localMax(city, 6), aux1)


         end do ! do vertice = 1, size(auxMeteoPoly%vertices)
      end if ! associated(auxMeteoPoly%vertices))

      auxMeteoPoly => auxMeteoPoly%prev      
   end do ! do while(associated(auxMeteoPoly)

   
   call GatherMeteogram(master_num,                  &
                        mchnum,                      &
                        nCities,                     &
                        nvarOutput,                  &
                        localSum,                    &
                        localMin,                    &
                        localMax,                    &
                        meteogramOutUnit,            &
                        totalSum,                    &
                        totalMin,                    &
                        totalMax)

   if(mchnum .eq. master_num)then
   
      write(meteogramOutUnit,*)'nCities: ',nCities
   	 	write(meteogramOutUnit, *)time, iyear1, imonth1, idate1, itime1	
      write(meteogramOutUnit,*)&
			'polygon_name - st points    Lat        Lon       Uave    Umax    Umin    '//  &
			'Vave    Vmax    Vmin Tave   Tmax  Tmin TDave TDmax TDmin SLPave SLPmax SLPmin PRECave PRECmax PRECmin'
   	do city = 1, nCities

   	 write(meteogramOutUnit,'(A,1x,i6,1x,f9.5,1x,f10.5,1x,6(f7.2,1x),6(f5.1,1x),3(f6.1,1x),3(f6.2,1x))') trim(cityName(city))  &
            , nVertices(city), cityLat(city), cityLon(city) &
            ,totalSum(city, 1)/nVertices(city),totalMax(city, 1),totalMin(city, 1) &
            ,totalSum(city, 2)/nVertices(city),totalMax(city, 2),totalMin(city, 2) &
            ,totalSum(city, 3)/nVertices(city),totalMax(city, 3),totalMin(city, 3) &
            ,totalSum(city, 4)/nVertices(city),totalMax(city, 4),totalMin(city, 4) &
            ,totalSum(city, 5)/nVertices(city),totalMax(city, 5),totalMin(city, 5) &
            ,totalSum(city, 6)/nVertices(city),totalMax(city, 6),totalMin(city, 6)
   	end do
   end if

!   if(mchnum .eq. master_num)then
!   	do city = 1, nCities
!        	write(meteogramOutUnit,*)'---------------------------BEGIN METEOGRAM -----------------------------'
!          	write(meteogramOutUnit,*)'polygon name: '//trim(cityName(city))//' used points: ', nVertices(city), cityLat(city), cityLon(city)
! 	 	write(meteogramOutUnit, *)time, iyear1, imonth1, idate1, itime1	
!          	do nv = 1, nvarOutput
!             		write(meteogramOutUnit,*) 'var: ', varNamesOutput(nv)
!             		write(meteogramOutUnit,*) 'sum-avg: ', totalSum(city, nv), totalSum(city, nv)/nVertices(city)
!             		write(meteogramOutUnit,*) 'max: ', totalMax(city, nv)
!             		write(meteogramOutUnit,*) 'min: ', totalMin(city, nv)
!          	end do
!	end do
!   end if

   
   deallocate(nVertices)
   deallocate(cityName)
   deallocate(cityLat)
   deallocate(cityLon)

   deallocate(localSum)
   deallocate(localMin)
   deallocate(localMax)

   deallocate(totalSum)
   deallocate(totalMin)
   deallocate(totalMax)    


 end subroutine ProcessLocalMeteogram
 
  subroutine LoadPolygonFromFile(fileUnit, filename, userMap)
    integer, intent(in)                                    :: fileUnit
    character(len=*), intent(in)                           :: filename
    type(tPolygon), dimension(:), allocatable, intent(out) :: userMap
   
    integer :: poly
    integer :: nPolygons
    integer :: allocStat
    integer :: pop
    real    :: loni
    real    :: lonf
    real    :: lati
    real    :: latf
    integer :: vertice
    integer :: nVertices
 
     open(unit=fileUnit, file=trim(filename), status='old')
     read(fileUnit,*) nPolygons
 
     allocate(userMap(nPolygons), stat=allocStat)
 
     do poly = 1, nPolygons
      read(fileUnit,'(A)')userMap(poly)%name
      read(fileUnit,*)pop
      read(fileUnit,*)loni, lati
      read(fileUnit,*)lonf, latf
      read(fileUnit,*)nVertices
    
      userMap(poly)%lon = (lonf + loni)/2.
      userMap(poly)%lat = (latf + lati)/2.
      
      allocate(userMap(poly)%points(nVertices), stat=allocStat)
 
      do vertice = 1, nVertices
       read(fileUnit,*)userMap(poly)%points(vertice)%x, userMap(poly)%points(vertice)%y
      end do
 
     end do
 
    close(fileUnit)
  
  end subroutine LoadPolygonFromFile

end module meteogram
! RMF end.
