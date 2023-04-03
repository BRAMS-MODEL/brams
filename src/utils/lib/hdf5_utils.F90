module hdf5_utils
#ifdef HDF5
   use HDF5
#endif 
   implicit none 

   private

   public :: h5_readNDVI

   contains

   subroutine h5_readNDVI(nxll, nyll, data, file_name, h5name)
     integer, intent(IN)                      :: nxll
     integer, intent(IN)                      :: nyll
     real, dimension(nxll, nyll), intent(OUT) :: data(:,:)
     character(*), intent(IN)                 :: file_name
     character(*), intent(IN)                 :: h5name

#ifndef HDF5
    write(*,fmt='(A)') "To use HDF5 the model must be compiled with HDF5."
    stop
#endif


#ifdef HDF5
     integer(HID_T)     :: file_id         ! File identifier
     integer(HID_T)     :: dset_id         ! Dataset identifier
     integer            :: nmembers        ! Number of group members
     character(LEN=256) :: name_buffer     ! Buffer to hold object's name
     integer            :: otype           ! Type of the object
     integer(HID_T)     :: dspace_id       ! dataspace identifier
     integer            :: ndims           ! number of dimensions
     integer(Hsize_t), dimension(2)    :: datadims, maxdatadims ! size of each dimension
     integer            :: i, error

     ! Initialize variables ---------------------------------------
     nmembers    = 0
     name_buffer = ""
     otype       = 0
     dspace_id   = -1
     ndims       = 0
     datadims    = 0
     maxdatadims = 0
     i           = 0
     error       = 0 
     ! ------------------------------------------------------------
 
     ! Open file "file_name" read only 
     CALL h5open_f(error)
     CALL h5fopen_f(file_name, H5F_ACC_RDONLY_F, file_id, error)

     CALL h5gn_members_f(file_id, "/", nmembers, error)

     if( nmembers > 1) then
       print *, "hdf5_utils: h5_readNDVI error!"
       print *, "Number of data set are ", nmembers
       print *, "Expected 1 data set!"
       print *, "List of data sets in root group"
       do i = 0, nmembers - 1
         CALL h5gget_obj_info_idx_f(file_id, "/", i, name_buffer, otype, error)
	 write(*,*) name_buffer
       end do
       call exit()
     end if 

     ! Get 1st and only dataset name
     CALL h5gget_obj_info_idx_f(file_id, "/", 0, name_buffer, otype, error)

     ! load dataset identifier
     CALL h5dopen_f(file_id, name_buffer, dset_id, error)
     if(trim(name_buffer) /= trim(h5name)) then
       print *, "hdf5_utils: h5_readNDVI error!"
       print *, "Dasaset name differ to the HEADER"
       print *, "HEADER      = ", h5name
       print *, "Inside file = ", name_buffer
       call exit()
     end if

     ! get dataspace id
     CALL h5dget_space_f(dset_id,dspace_id,error)

     ! get dimensions of dataset
     CALL h5sget_simple_extent_ndims_f(dspace_id,ndims,error)
     if( ndims /= 2) then
       print *, "hdf5_utils: h5_readNDVI error!"
       print *, "Expected array 2 dimensions"
       print *, "Number of dimensions = ", ndims
       call exit()
     end if

     CALL h5sget_simple_extent_dims_f(dspace_id,datadims,maxdatadims,error)
     if(nxll /= datadims(1) .OR. nyll /= datadims(2)) then
       print *, "hdf5_utils: h5_readNDVI error!"
       print *, "Dimension(s) of dataset differ to the HEADER"
       print *, "HEADER      dim1=",nxll," dim2=",nyll
       print *, "Inside file dim1=",datadims(1)," dim2=",datadims(2)
       call exit()
     end if
    
     ! read data
     CALL h5dread_f(dset_id, H5T_NATIVE_REAL, data, maxdatadims, error)
    
     CALL h5dclose_f(dset_id, error)
     CALL h5fclose_f(file_id, error)
     CALL h5close_f(error)
   
#endif
   end subroutine h5_readNDVI


end module hdf5_utils
