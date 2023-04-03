#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


PROGRAM BRAMS

!DSM USE netcdf, ONLY: nf90_open, NF90_NOWRITE
!DSM implicit none
!DSM integer :: STATUS,ncid,ncid1,ncid2,nf_INQ_DIMLEN,nx,ny,nt,ntime,ntype,var_id,nf_inq_varid,nf_get_var
!DSM character :: filename*80, filename1*80, filename2*80
!DSM real, allocatable :: glat(:,:), glon(:,:), land_frac(:,:), totsnow(:,:,:), swdown(:,:,:)
!DSM real, allocatable :: lwdown(:,:,:), totrain(:,:,:), temp(:,:,:), wind(:,:,:), pstar(:,:,:)
!DSM real, allocatable :: q(:,:,:), frac(:,:,:)
!DSM 
!DSM filename='/home/jonathan_alves/IC/JULES/from_BRAMS/grid_info_BRAMS.nc'
!DSM 
!DSM STATUS = nf90_open(TRIM(filename), NF90_NOWRITE, ncid); print*,'nf90_open:',STATUS
!DSM STATUS = nf_INQ_DIMLEN(ncid, 1, nx); print*,'nf_INQ_DIMLEN-nx:',STATUS                
!DSM STATUS = nf_INQ_DIMLEN(ncid, 2, ny); print*,'nf_INQ_DIMLEN-ny:',STATUS
!DSM STATUS = nf_INQ_DIMLEN(ncid, 3, nt); print*,'nf_INQ_DIMLEN-nt:',STATUS
!DSM 
!DSM if (nx>1000 .or. nx == 0) then
!DSM         print*, "Erro na leitura de grid_info_BRAMS.nc, nx=", nx
!DSM         stop
!DSM endif
!DSM 
!DSM print*,ncid,nx,ny,nt
!DSM      
!DSM allocate (glat(nx,ny), glon(nx,ny), land_frac(nx,ny))                                                         !declarando a dimensao da matriz de numeros reais
!DSM 
!DSM STATUS = nf_inq_varid(ncid, 'latitude', var_id) ; print*, 'nf_inq_varid-latitude:', STATUS
!DSM STATUS = nf_get_var(ncid, var_id, glat) ; print*, 'nf_inq_var-latitude:', STATUS
!DSM STATUS = nf_inq_varid(ncid, 'longitude', var_id) ; print*, 'nf_inq_varid-longitude:', STATUS
!DSM STATUS = nf_get_var(ncid, var_id, glon) ; print*, 'nf_inq_var-longitude:', STATUS
!DSM STATUS = nf_inq_varid(ncid, 'land_frac', var_id) ; print*, 'nf_inq_varid-land_frac:', STATUS
!DSM STATUS = nf_get_var(ncid, var_id, land_frac) ; print*, 'nf_inq_var-land_frac:', STATUS
!DSM 
!DSM !filename1='/home/jonathan_alves/IC/JULES/from_BRAMS/drive_from_BRAMS_2010010300.nc'
!DSM !
!DSM !STATUS = nf90_open(TRIM(filename1), NF90_NOWRITE, ncid1); print*,'nf90_open:',STATUS
!DSM !STATUS = nf_INQ_DIMLEN(ncid1, 1, nx); print*,'nf_INQ_DIMLEN-nx:',STATUS                
!DSM !STATUS = nf_INQ_DIMLEN(ncid1, 2, ny); print*,'nf_INQ_DIMLEN-ny:',STATUS
!DSM !STATUS = nf_INQ_DIMLEN(ncid1, 3, nt); print*,'nf_INQ_DIMLEN-nt:',STATUS
!DSM !
!DSM !if (nx>1000 .or. nx == 0) then
!DSM !        print*, "Erro na leitura de drive_from_BRAMS.nc, nx=", nx
!DSM !        stop
!DSM !endif
!DSM !
!DSM !print*,ncid1,nx,ny,nt
!DSM !
!DSM allocate (totsnow(nx,ny,nt), swdown(nx,ny,nt), lwdown(nx,ny,nt), totrain(nx,ny,nt))
!DSM allocate (wind(nx,ny,nt), pstar(nx,ny,nt), q(nx,ny,nt), temp(nx,ny,nt))
!DSM !
!DSM !STATUS = nf_inq_varid(ncid1, 'swdown', var_id) ; print*, 'nf_inq_varid-swdown:', STATUS
!DSM !STATUS = nf_get_var(ncid1, var_id, swdown) ; print*, 'nf_inq_var-swdown:', STATUS
!DSM !STATUS = nf_inq_varid(ncid1, 'lwdown', var_id) ; print*, 'nf_inq_varid-lwdown:', STATUS
!DSM !STATUS = nf_get_var(ncid1, var_id, lwdown) ; print*, 'nf_inq_var-lwdown:', STATUS
!DSM !STATUS = nf_inq_varid(ncid1, 'totrain', var_id) ; print*, 'nf_inq_varid-totrain:', STATUS
!DSM !STATUS = nf_get_var(ncid1, var_id, totrain) ; print*, 'nf_inq_var-totrain:', STATUS
!DSM !STATUS = nf_inq_varid(ncid1, 'totsnow', var_id) ; print*, 'nf_inq_varid-totsnow:', STATUS
!DSM !STATUS = nf_get_var(ncid1, var_id, totsnow) ; print*, 'nf_inq_var-totsnow:', STATUS
!DSM !STATUS = nf_inq_varid(ncid1, 'temp', var_id) ; print*, 'nf_inq_varid-temp:', STATUS
!DSM !STATUS = nf_get_var(ncid1, var_id, temp) ; print*, 'nf_inq_var-temp:', STATUS
!DSM !STATUS = nf_inq_varid(ncid1, 'wind', var_id) ; print*, 'nf_inq_varid-wind:', STATUS
!DSM !STATUS = nf_get_var(ncid1, var_id, wind) ; print*, 'nf_inq_var-wind:', STATUS
!DSM !STATUS = nf_inq_varid(ncid1, 'pstar', var_id) ; print*, 'nf_inq_varid-pstar:', STATUS
!DSM !STATUS = nf_get_var(ncid1, var_id, pstar) ; print*, 'nf_inq_var-pstar:', STATUS
!DSM !STATUS = nf_inq_varid(ncid1, 'q', var_id) ; print*, 'nf_inq_varid-q:', STATUS
!DSM !STATUS = nf_get_var(ncid1, var_id, q) ; print*, 'nf_inq_var-q:', STATUS
!DSM 
!DSM filename2='/home/jonathan_alves/IC/JULES/from_BRAMS/frac.nc'
!DSM 
!DSM STATUS = nf90_open(TRIM(filename2), NF90_NOWRITE, ncid2); print*,'nf90_open:',STATUS
!DSM STATUS = nf_INQ_DIMLEN(ncid2, 1, ntype); print*,'nf_INQ_DIMLEN-ntype:',STATUS                
!DSM STATUS = nf_INQ_DIMLEN(ncid2, 2, ntime); print*,'nf_INQ_DIMLEN-ntime:',STATUS
!DSM STATUS = nf_INQ_DIMLEN(ncid2, 3, ny); print*,'nf_INQ_DIMLEN-ny:',STATUS
!DSM STATUS = nf_INQ_DIMLEN(ncid2, 4, nx); print*,'nf_INQ_DIMLEN-nx:',STATUS
!DSM 
!DSM if (nx>1000 .or. nx == 0) then
!DSM         print*, "Erro na leitura de frac.nc, nx=", nx
!DSM         stop
!DSM endif
!DSM 
!DSM print*,ntype,ntime,ny,nx
!DSM 
!DSM allocate (frac(nx,ny,ntype))
!DSM 
!DSM STATUS = nf_inq_varid(ncid2, 'frac', var_id) ; print*, 'nf_inq_varid-frac:', STATUS
!DSM STATUS = nf_get_var(ncid2, var_id, frac) ; print*, 'nf_inq_var-frac:', STATUS
!DSM 
!DSM !glat,glon,land_frac,totsnow,swdown,lwdown,totrain,temp,wind,pstar,q,frac
!DSM call sfclyr_jules(glat,glon,nx,ny,nt,ntype,land_frac,totsnow,swdown,lwdown,totrain,temp,wind,pstar,q,frac)

END PROGRAM BRAMS


!DSM SUBROUTINE read_drive_BRAMS(FILE,varID,var)
!DSM USE netcdf, ONLY: nf90_open, NF90_NOWRITE
!DSM USE mem_brams_jules, ONLY: nxB,nyB,ntB,swdownB,lwdownB,totrainB,totsnowB,tempB,windB,pstarB,qB
!DSM USE file_mod
!DSM USE datetime_mod, ONLY: datetime_str_len, datetime, datetime_from_string
!DSM USE model_time_mod, ONLY:                                                     &
!DSM   current_time, is_spinup, spinup_cycle
!DSM USE io_constants, ONLY:                                                       &
!DSM   max_file_name_len, max_dim_var, mode_write
!DSM USE string_utils_mod, ONLY:                                                   &
!DSM   to_string
!DSM implicit none
!DSM TYPE(datetime) :: data_dt  ! The datetime for data that applies for
!DSM TYPE(file_handle), INTENT(IN) :: FILE  ! The file to read from
!DSM character :: filename*180
!DSM integer :: STATUS,ncid,nf_INQ_DIMLEN,nx,ny,var_id,nf_inq_varid,nf_get_var
!DSM integer :: varID, i, j, k
!DSM real, allocatable :: swdown(:,:), lwdown(:,:), totrain(:,:), totsnow(:,:)
!DSM real, allocatable :: temp(:,:), wind(:,:), pstar(:,:), q(:,:)
!DSM real :: var(nxB*nyB)
!DSM CHARACTER(LEN=10) :: dt_string
!DSM 
!DSM 
!DSM WRITE(dt_string, '(I4.4,I2.2,I2.2,I2.2)') current_time%year,                     &
!DSM                                        current_time%month,                    &
!DSM                                        current_time%day,current_time%time/3600
!DSM !  dt_string = TRIM(dt_string) // "." // TRIM(to_string(current_time%time))
!DSM !Print*,current_time%year,current_time%month,current_time%day,current_time%time
!DSM    filename = '/home/jonathan_alves/IC/JULES/run_espacial/from_BRAMS/drive_from_BRAMS_'//TRIM(dt_string)//'.nc'
!DSM 
!DSM    Print*, 'file name=', trim(filename)
!DSM    STATUS = nf90_open(TRIM(filename), NF90_NOWRITE, ncid)!; print*,'nf90_open:',STATUS
!DSM    STATUS = nf_INQ_DIMLEN(ncid, 1, nx)!; print*,'nf_INQ_DIMLEN-nx:',STATUS                
!DSM    STATUS = nf_INQ_DIMLEN(ncid, 2, ny)!; print*,'nf_INQ_DIMLEN-ny:',STATUS
!DSM 
!DSM    allocate(swdown(nx,ny),lwdown(nx,ny),totrain(nx,ny))
!DSM    allocate(totsnow(nx,ny),temp(nx,ny),wind(nx,ny),pstar(nx,ny),q(nx,ny))
!DSM 
!DSM if (varID == 1) then
!DSM    
!DSM    STATUS = nf_inq_varid(ncid, 'swdown', var_id)!; print*, 'nf_inq_varid-swdown:', STATUS
!DSM    STATUS = nf_get_var(ncid, var_id, swdown)!; print*, 'nf_inq_var-swdown:', STATUS
!DSM    k=1
!DSM    do j = 1 , ny
!DSM       do i = 1 , nx
!DSM          var(k) = swdown(i,j)
!DSM          k=k+1
!DSM       enddo
!DSM    enddo
!DSM 
!DSM else if (varID == 2) then
!DSM    
!DSM    STATUS = nf_inq_varid(ncid, 'lwdown', var_id)!; print*, 'nf_inq_varid-lwdown:', STATUS
!DSM    STATUS = nf_get_var(ncid, var_id, lwdown)!; print*, 'nf_inq_var-lwdown:', STATUS
!DSM    k=1
!DSM    do j = 1 , ny
!DSM       do i = 1 , nx
!DSM          var(k) = lwdown(i,j)
!DSM          k=k+1
!DSM       enddo
!DSM    enddo
!DSM 
!DSM else if (varID == 3) then
!DSM 
!DSM    STATUS = nf_inq_varid(ncid, 'totrain', var_id)!; print*, 'nf_inq_varid-totrain:', STATUS
!DSM    STATUS = nf_get_var(ncid, var_id, totrain)!; print*, 'nf_inq_var-totrain:', STATUS
!DSM    k=1
!DSM    do j = 1 , ny
!DSM       do i = 1 , nx
!DSM          var(k) = totrain(i,j)
!DSM          k=k+1
!DSM       enddo
!DSM    enddo
!DSM 
!DSM else if (varID == 4) then
!DSM 
!DSM    STATUS = nf_inq_varid(ncid, 'totsnow', var_id)!; print*, 'nf_inq_varid-totsnow:', STATUS
!DSM    STATUS = nf_get_var(ncid, var_id, totsnow)!; print*, 'nf_inq_var-totsnow:', STATUS
!DSM    k=1
!DSM    do j = 1 , ny
!DSM       do i = 1 , nx
!DSM          var(k) = totsnow(i,j)
!DSM          k=k+1
!DSM       enddo
!DSM    enddo
!DSM 
!DSM else if (varID == 5) then
!DSM 
!DSM    STATUS = nf_inq_varid(ncid, 'temp', var_id)!; print*, 'nf_inq_varid-temp:', STATUS
!DSM    STATUS = nf_get_var(ncid, var_id, temp)!; print*, 'nf_inq_var-temp:', STATUS
!DSM    k=1
!DSM    do j = 1 , ny
!DSM       do i = 1 , nx
!DSM          var(k) = temp(i,j)
!DSM          k=k+1
!DSM       enddo
!DSM    enddo
!DSM 
!DSM else if (varID == 6) then
!DSM 
!DSM    STATUS = nf_inq_varid(ncid, 'wind', var_id)!; print*, 'nf_inq_varid-wind:', STATUS
!DSM    STATUS = nf_get_var(ncid, var_id, wind)!; print*, 'nf_inq_var-wind:', STATUS
!DSM    k=1
!DSM    do j = 1 , ny
!DSM       do i = 1 , nx
!DSM          var(k) = wind(i,j)
!DSM          k=k+1
!DSM       enddo
!DSM    enddo
!DSM 
!DSM else if (varID == 7) then
!DSM    
!DSM    STATUS = nf_inq_varid(ncid, 'pstar', var_id)!; print*, 'nf_inq_varid-pstar:', STATUS
!DSM    STATUS = nf_get_var(ncid, var_id, pstar)!; print*, 'nf_inq_var-pstar:', STATUS
!DSM    k=1
!DSM    do j = 1 , ny
!DSM       do i = 1 , nx
!DSM          var(k) = pstar(i,j)
!DSM          k=k+1
!DSM       enddo
!DSM    enddo
!DSM 
!DSM else if (varID == 8) then
!DSM    
!DSM    STATUS = nf_inq_varid(ncid, 'q', var_id)!; print*, 'nf_inq_varid-q:', STATUS
!DSM    STATUS = nf_get_var(ncid, var_id, q)!; print*, 'nf_inq_var-q:', STATUS
!DSM    k=1
!DSM    do j = 1 , ny
!DSM       do i = 1 , nx
!DSM          var(k) = q(i,j)
!DSM          k=k+1
!DSM       enddo
!DSM    enddo
!DSM 
!DSM end if
!DSM 
!DSM  
!DSM END SUBROUTINE read_drive_BRAMS

#endif
