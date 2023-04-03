MODULE HDF5UTILS
       INTERFACE hdf5_attr
              MODULE PROCEDURE hdf5_attr_r, hdf5_attr_i
       END INTERFACE hdf5_attr

       INTERFACE hdf5_read
              MODULE PROCEDURE hdf5_read_r
       END INTERFACE hdf5_read
CONTAINS

      SUBROUTINE hdf5_attr_r(attribute, dset, value)

        USE HDF5

        IMPLICIT NONE

        INTEGER :: hdferr, i
        INTEGER(HID_T), INTENT(IN)  :: dset
        INTEGER(HID_T)              :: space, attr
        CHARACTER(*), INTENT(IN) :: attribute
        INTEGER(HSIZE_T), DIMENSION(1) :: maxdims
        INTEGER(hsize_t), DIMENSION(1) :: dims = (/1/)
        REAL, DIMENSION(:), ALLOCATABLE, TARGET :: attr_data_r
        REAL, INTENT(OUT) :: value

        CALL h5aopen_f(dset, trim(attribute), attr, hdferr)
        CALL checkError(hdferr, 'Erro ao abrir atributo')     
        CALL h5aget_space_f(attr, space, hdferr)
        CALL checkError(hdferr, 'Erro ao determinar o espaço do atributo')   
        CALL h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)
        CALL checkError(hdferr, 'Erro ao determinar tamanho da memória para alocar o valor do atributo')   

        ALLOCATE(attr_data_r(1:dims(1)))
       
        CALL h5aread_f( attr, H5T_NATIVE_REAL, attr_data_r, dims, hdferr)
        CALL checkError(hdferr, 'Erro ao ler o valor do atributo')
        
        value = attr_data_r(1)

        DEALLOCATE(attr_data_r)

        CALL h5aclose_f(attr , hdferr)
        CALL h5sclose_f(space, hdferr)

      END SUBROUTINE hdf5_attr_r

      SUBROUTINE hdf5_attr_i(attribute, dset, value)

        USE HDF5

        IMPLICIT NONE

        INTEGER :: hdferr, i
        INTEGER(HID_T), INTENT(IN)  :: dset
        INTEGER(HID_T)              :: space, attr
        CHARACTER(*), INTENT(IN) :: attribute
        INTEGER(HSIZE_T), DIMENSION(1) :: maxdims
        INTEGER(hsize_t),   DIMENSION(1) :: dims = (/1/)
        INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: attr_data_i
        INTEGER, INTENT(OUT) :: value

        CALL h5aopen_f(dset, trim(attribute), attr, hdferr)
        CALL checkError(hdferr, 'Erro ao abrir atributo')     
        CALL h5aget_space_f(attr, space, hdferr)
        CALL checkError(hdferr, 'Erro ao determinar o espaço do atributo')   
        CALL h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)
        CALL checkError(hdferr, 'Erro ao determinar tamanho da memória para alocar o valor do atributo')   

        ALLOCATE(attr_data_i(1:dims(1)))
       
        CALL h5aread_f( attr, H5T_NATIVE_INTEGER, attr_data_i, dims, hdferr)
        CALL checkError(hdferr, 'Erro ao ler o valor do atributo')           

        value = attr_data_i(1)

        DEALLOCATE(attr_data_i)

        CALL h5aclose_f(attr , hdferr)
        CALL h5sclose_f(space, hdferr)

      END SUBROUTINE hdf5_attr_i

      SUBROUTINE hdf5_read_r(file, dataset, data)
  
        USE HDF5

        IMPLICIT NONE

        INTEGER(HSIZE_T), DIMENSION(1) :: data_dims 
        INTEGER(HID_T), INTENT(IN)  :: file
        INTEGER(HID_T)  :: dset
        CHARACTER(*), INTENT(IN)  :: dataset
        INTEGER, DIMENSION(:), ALLOCATABLE :: idata
        REAL, DIMENSION(:), POINTER :: data
        REAL :: scaling_factor
        INTEGER ::  nx, ny, hdferr

        CALL h5dopen_f (file, trim(dataset), dset, hdferr)
        CALL checkError(hdferr, 'Erro ao abrir dataset')   

        CALL hdf5_attr("SCALING_FACTOR", dset, scaling_factor)
        CALL hdf5_attr("N_COLS", dset, nx)
        CALL hdf5_attr("N_LINES", dset, ny)

        data_dims(1) = ny
        
!ALTEREI PARA O UBUNTU
!        IF ( associated(data) ) THEN
!          print *, 'Data is currently associated'
!          nullify(data)
!          DEALLOCATE(data)
!        ENDIF

        ALLOCATE(idata(ny), data(ny))
        
        CALL h5dread_f(dset, H5T_NATIVE_INTEGER, idata, data_dims, hdferr)
        CALL checkError(hdferr, 'Erro ao ler o dado do dataset')           
        
        data(:) = float(idata(:))/(scaling_factor)

        DEALLOCATE(idata)

        CALL h5dclose_f(dset , hdferr)

        IF(hdferr .LT. 0) THEN              
          CALL h5fclose_f(file , hdferr)
          CALL h5close_f(hdferr)
          STOP
        ENDIF

      END SUBROUTINE hdf5_read_r

END MODULE HDF5UTILS
