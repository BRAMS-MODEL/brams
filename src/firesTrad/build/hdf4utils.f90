MODULE HDF4UTILS

	INTERFACE hdf4_read
		MODULE PROCEDURE hdf4_read_r
		MODULE PROCEDURE hdf4_read_i
		MODULE PROCEDURE hdf4_read_uint8
		MODULE PROCEDURE hdf4_read_int16
	END INTERFACE hdf4_read

CONTAINS

	SUBROUTINE hdf4_read_r(sd_id, sds_name, output)
    IMPLICIT NONE

		! I/O:
		INTEGER, INTENT(IN)  :: sd_id
		CHARACTER(*), INTENT(IN)  	:: sds_name
		REAL, DIMENSION(:), POINTER :: output

		! Funções:
		INTEGER		:: sfn2index, sfselect, sfrdata, sfendacc
		INTEGER		:: sffattr, sfrnatt, sfginfo, sfgdinfo

		! Variáveis auxiliares:
		INTEGER		:: n_dims, dim_id, dim_size, op_size, attr_index, i, firepix
		INTEGER		:: sds_index, sds_id, ierr, n_attrs, data_type, rank
    CHARACTER	(len=64):: ds_name

		INTEGER, DIMENSION(:), ALLOCATABLE :: edges, start, stride, dim_sizes

		! Teste de existência de dados no arquivo.
		attr_index = sffattr(sd_id, 'FirePix')
		CALL checkError(attr_index, "Erro ao obter o índice do atributo FirePix")
		ierr = sfrnatt(sd_id, attr_index, firepix)
		CALL checkError(ierr, "Erro ao ler o atributo FirePix.")
	
		IF (firepix .EQ. 0) THEN
			PRINT *, "Arquivo vazio."
			STOP
		ENDIF

		! Obtenção o índice do sds a partir do nome.
		! Pode haver erro se existir mais de um sds com mesmo nome,
		! ver seções 3.7.4 e 3.7.5 da documentação do HDF4.
		sds_index = sfn2index(sd_id, sds_name)

		! Seleção do SDS a partir do índice.
		sds_id = sfselect(sd_id, sds_index)
		CALL checkError(sds_id, "Erro ao abrir o data set.")
		
		! Determinação de start, stride, edges.
		ALLOCATE(dim_sizes(1))
		ALLOCATE(start(firepix), stride(firepix), edges(firepix))

		start(:) = 0
		stride(:) = 1
		edges(:) = firepix

		ierr = sfginfo(sds_id, ds_name, rank, dim_sizes, data_type, n_attrs)

		! Teste de ocupação do array de saída.
! ALTEREI PARA O UBUNTU
!		IF ( ASSOCIATED(output) ) THEN
!			PRINT *, 'Ponteiro para o array de saída está associado. Desalocando...'
!			DEALLOCATE(output)
!		ENDIF

		! Alocação do vetor de saída.
		ALLOCATE(output(firepix))

		! Leitura dos dados para o array 'output'.
		ierr = sfrdata (sds_id, start, stride, edges, output)
    CALL checkError(ierr, 'Erro na leitura do data set.')
		
		! Desalocação dos arrays.
		DEALLOCATE(start, stride, edges)

		! Término do acesso ao data set.
		ierr = sfendacc(sds_id)
		CALL checkError(ierr, 'Erro ao fechar data set.')
	END SUBROUTINE hdf4_read_r

	SUBROUTINE hdf4_read_i(sd_id, sds_name, output)
    IMPLICIT NONE

		! I/O:
		INTEGER, INTENT(IN)  :: sd_id
		CHARACTER(*), INTENT(IN)  	:: sds_name
		INTEGER, DIMENSION(:), POINTER :: output

		! Funções:
		INTEGER		:: sfn2index, sfselect, sfrdata, sfendacc
		INTEGER		:: sffattr, sfrnatt, sfginfo, sfgdinfo

		! Variáveis auxiliares:
		INTEGER		:: n_dims, dim_id, dim_size, op_size, attr_index, i, firepix
		INTEGER		:: sds_index, sds_id, ierr, n_attrs, data_type, rank
    CHARACTER	(len=64):: ds_name

		INTEGER, DIMENSION(:), ALLOCATABLE :: edges, start, stride, dim_sizes

		! Teste de existência de dados no arquivo.
		attr_index = sffattr(sd_id, 'FirePix')
		CALL checkError(attr_index, "Erro ao obter o índice do atributo FirePix")
		ierr = sfrnatt(sd_id, attr_index, firepix)
		CALL checkError(ierr, "Erro ao ler o atributo FirePix.")
	
		IF (firepix .EQ. 0) THEN
			PRINT *, "Arquivo vazio."
			STOP
		ENDIF

		! Obtenção o índice do sds a partir do nome.
		! Pode haver erro se existir mais de um sds com mesmo nome,
		! ver seções 3.7.4 e 3.7.5 da documentação do HDF4.
		sds_index = sfn2index(sd_id, sds_name)

		! Seleção do SDS a partir do índice.
		sds_id = sfselect(sd_id, sds_index)
		CALL checkError(sds_id, "Erro ao abrir o data set.")
		
		! Determinação de start, stride, edges.
		ALLOCATE(dim_sizes(firepix))
		ALLOCATE(start(firepix), stride(firepix), edges(firepix))

		dim_sizes(1) = firepix
		start(:) = 0
		stride(:) = 1
		edges(:) = firepix

		ierr = sfginfo(sds_id, ds_name, rank, dim_sizes, data_type, n_attrs)

		! Teste de ocupação do array de saída.

!ALTEREI PARA O UBUNTU
!		IF ( ASSOCIATED(output) ) THEN
!			PRINT *, 'Ponteiro para o array de saída está associado. Desalocando...'
!			DEALLOCATE(output)
!		ENDIF

		! Alocação do vetor de saída.
		ALLOCATE(output(firepix))

		! Leitura dos dados para o array 'output'.
		ierr = sfrdata (sds_id, start, stride, edges, output)
    CALL checkError(ierr, 'Erro na leitura do data set.')
		
		! Desalocação dos arrays.
		DEALLOCATE(start, stride, edges)

		! Término do acesso ao data set.
		ierr = sfendacc(sds_id)
		CALL checkError(ierr, 'Erro ao fechar data set.')
	END SUBROUTINE hdf4_read_i

	SUBROUTINE hdf4_read_int16(sd_id, sds_name, output)
    IMPLICIT NONE

		! I/O:
		INTEGER, INTENT(IN)  :: sd_id
		CHARACTER(*), INTENT(IN)  	:: sds_name
		INTEGER*2, DIMENSION(:), POINTER :: output

		! Funções:
		INTEGER		:: sfn2index, sfselect, sfrdata, sfendacc
		INTEGER		:: sffattr, sfrnatt, sfginfo, sfgdinfo

		! Variáveis auxiliares:
		INTEGER		:: n_dims, dim_id, dim_size, op_size, attr_index, i, firepix
		INTEGER		:: sds_index, sds_id, ierr, n_attrs, data_type, rank
    CHARACTER	(len=64):: ds_name

		INTEGER, DIMENSION(:), ALLOCATABLE :: edges, start, stride, dim_sizes

		! Teste de existência de dados no arquivo.
		attr_index = sffattr(sd_id, 'FirePix')
		CALL checkError(attr_index, "Erro ao obter o índice do atributo FirePix")
		ierr = sfrnatt(sd_id, attr_index, firepix)
		CALL checkError(ierr, "Erro ao ler o atributo FirePix.")
	
		IF (firepix .EQ. 0) THEN
			PRINT *, "Arquivo vazio."
			STOP
		ENDIF

		! Obtenção o índice do sds a partir do nome.
		! Pode haver erro se existir mais de um sds com mesmo nome,
		! ver seções 3.7.4 e 3.7.5 da documentação do HDF4.
		sds_index = sfn2index(sd_id, sds_name)

		! Seleção do SDS a partir do índice.
		sds_id = sfselect(sd_id, sds_index)
		CALL checkError(sds_id, "Erro ao abrir o data set.")
		
		! Determinação de start, stride, edges.
		ALLOCATE(dim_sizes(firepix))
		ALLOCATE(start(firepix), stride(firepix), edges(firepix))

		dim_sizes(1) = firepix
		start(:) = 0
		stride(:) = 1
		edges(:) = firepix

		ierr = sfginfo(sds_id, ds_name, rank, dim_sizes, data_type, n_attrs)

		! Teste de ocupação do array de saída.
!UBUNTU
!		IF ( ASSOCIATED(output) ) THEN
!			PRINT *, 'Ponteiro para o array de saída está associado. Desalocando...'
!			DEALLOCATE(output)
!		ENDIF

		! Alocação do vetor de saída.
		ALLOCATE(output(firepix))

		! Leitura dos dados para o array 'output'.
		ierr = sfrdata (sds_id, start, stride, edges, output)
    CALL checkError(ierr, 'Erro na leitura do data set.')
		
		! Desalocação dos arrays.
		DEALLOCATE(start, stride, edges)

		! Término do acesso ao data set.
		ierr = sfendacc(sds_id)
		CALL checkError(ierr, 'Erro ao fechar data set.')
	END SUBROUTINE hdf4_read_int16

	SUBROUTINE hdf4_read_uint8(sd_id, sds_name, output)
    IMPLICIT NONE

		! I/O:
		INTEGER, INTENT(IN)  :: sd_id
		CHARACTER(*), INTENT(IN)  	:: sds_name
		INTEGER*1, DIMENSION(:), POINTER :: output

		! Funções:
		INTEGER		:: sfn2index, sfselect, sfrdata, sfendacc
		INTEGER		:: sffattr, sfrnatt, sfginfo, sfgdinfo

		! Variáveis auxiliares:
		INTEGER		:: n_dims, dim_id, dim_size, op_size, attr_index, i, firepix
		INTEGER		:: sds_index, sds_id, ierr, n_attrs, data_type, rank
    CHARACTER	(len=64):: ds_name

		INTEGER, DIMENSION(:), ALLOCATABLE :: edges, start, stride, dim_sizes

		! Teste de existência de dados no arquivo.
		attr_index = sffattr(sd_id, 'FirePix')
		CALL checkError(attr_index, "Erro ao obter o índice do atributo FirePix")
		ierr = sfrnatt(sd_id, attr_index, firepix)
		CALL checkError(ierr, "Erro ao ler o atributo FirePix.")
	
		IF (firepix .EQ. 0) THEN
			PRINT *, "Arquivo vazio."
			STOP
		ENDIF

		! Obtenção o índice do sds a partir do nome.
		! Pode haver erro se existir mais de um sds com mesmo nome,
		! ver seções 3.7.4 e 3.7.5 da documentação do HDF4.
		sds_index = sfn2index(sd_id, sds_name)

		! Seleção do SDS a partir do índice.
		sds_id = sfselect(sd_id, sds_index)
		CALL checkError(sds_id, "Erro ao abrir o data set.")
		
		! Determinação de start, stride, edges.
		ALLOCATE(dim_sizes(firepix))
		ALLOCATE(start(firepix), stride(firepix), edges(firepix))

		dim_sizes(1) = firepix
		start(:) = 0
		stride(:) = 1
		edges(:) = firepix

		ierr = sfginfo(sds_id, ds_name, rank, dim_sizes, data_type, n_attrs)

		! Teste de ocupação do array de saída.
!UBUNTU
!		IF ( ASSOCIATED(output) ) THEN
!			PRINT *, 'Ponteiro para o array de saída está associado. Desalocando...'
!			DEALLOCATE(output)
!		ENDIF

		! Alocação do vetor de saída.
		ALLOCATE(output(firepix))

		! Leitura dos dados para o array 'output'.
		ierr = sfrdata (sds_id, start, stride, edges, output)
    CALL checkError(ierr, 'Erro na leitura do data set.')
		
		! Desalocação dos arrays.
		DEALLOCATE(start, stride, edges)

		! Término do acesso ao data set.
		ierr = sfendacc(sds_id)
		CALL checkError(ierr, 'Erro ao fechar data set.')
	END SUBROUTINE hdf4_read_uint8

	SUBROUTINE checkError(hdferr, msg)
  	IMPLICIT NONE

 		INTEGER 			:: hdferr
 		CHARACTER(*)	:: msg

 		IF(hdferr .LT. 0) THEN
  	 	WRITE(*,'("ERROR: ",A)') trim(msg)
 		ENDIF
	END SUBROUTINE checkError

END MODULE HDF4UTILS

