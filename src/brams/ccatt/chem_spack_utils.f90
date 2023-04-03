MODULE Spack_utils

  IMPLICIT NONE

  TYPE, PUBLIC ::  index_dim
    INTEGER,ALLOCATABLE,DIMENSION(:)     :: block_end
    INTEGER,ALLOCATABLE,DIMENSION(:,:)   :: indexi,indexj,indexk !Colapsed 3d whilst 2d
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: ijkindex
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: bindex
    INTEGER,ALLOCATABLE,DIMENSION(:,:)   :: kij_index ! index of tendency array
    DOUBLE PRECISION,POINTER,DIMENSION(:):: last_accepted_dt
  END TYPE

  PRIVATE

  INTEGER        , PUBLIC                           :: maxblock_size
  INTEGER        , PUBLIC, ALLOCATABLE,DIMENSION(:) :: nob
  TYPE(index_dim), PUBLIC, ALLOCATABLE,DIMENSION(:) :: index_g

  PUBLIC :: AllocIndex ! Subroutine

CONTAINS
 
  !Note from Lufla
  !For chem we can't have column dependencies, then we use collapse all atmosphere array (3d=k,i,j) 
  !to only 1d arrays (ijk) divided in blocks. Each block have only "block_size" size or few.
  !The arrays indexi, indexed for position and number_of_block return a true position x of cell in 
  !athmosphere. At same way the other 2 dimensions y and k we obtain the cell positions.
  !The array ijkindex and bindex mapped the inverse, i.e., we obtain the value of x,y and z, if we put
  !the value of block and ijk in it.

  !For example: the point x=2, y=5, z=6 for blocks with size=30 have nb=2 and ijk=30
  ! The array indexi(30,2)=2
  ! The array indexj(30,2)=5 and
  ! The array indexk(30,2)=6 while
  ! The array ijkindex(6,2,5)=30 and
  ! The array bindex(6,2,5)=2

  SUBROUTINE AllocIndex(block_size,mmxp,mmyp,mmzp,mia,miz,mja,mjz, &
!--(DMK-BRAMS-5.0-INI)------------------------------------------------------------    
! not used
!  	mibcon,&
!--(DMK-BRAMS-5.0-FIM)------------------------------------------------------------    
       ngrids,dtlongn,N_DYN_CHEM)

    INTEGER, INTENT(IN) :: block_size,N_DYN_CHEM
    INTEGER, INTENT(IN) :: mmxp(ngrids)
    INTEGER, INTENT(IN) :: mmyp(ngrids)
    INTEGER, INTENT(IN) :: mmzp(ngrids)
    INTEGER, INTENT(IN) :: mia(ngrids)
    INTEGER, INTENT(IN) :: miz(ngrids)
    INTEGER, INTENT(IN) :: mja(ngrids)
    INTEGER, INTENT(IN) :: mjz(ngrids)
!--(DMK-BRAMS-5.0-INI)------------------------------------------------------------    
! not used
!   INTEGER, INTENT(IN) :: mibcon(ngrids)
!--(DMK-BRAMS-5.0-FIM)------------------------------------------------------------    
    INTEGER, INTENT(IN) :: ngrids
    REAL   , INTENT(IN) :: dtlongn(ngrids)

    INTEGER :: m1,m2,m3,ia,iz,ja,jz
    INTEGER :: i1,j1,k1,cc,inb
    
!--(DMK-BRAMS-5.0-INI)------------------------------------------------------------    
! not used
!    INTEGER :: ibcon
!--(DMK-BRAMS-5.0-FIM)------------------------------------------------------------    

    INTEGER :: resto,i,ijk,ngr   


    !calculate the total os athmosphere points

    maxblock_size = 0
    ALLOCATE(index_g(ngrids),nob(ngrids))
    nob(:) = 0

    DO ngr=1,ngrids
       m2=mmxp(ngr)
       m3=mmyp(ngr)
       m1=mmzp(ngr)
       ia =mia(ngr)
       iz =miz(ngr)
       ja =mja(ngr)
       jz =mjz(ngr)
!--(DMK-BRAMS-5.0-INI)------------------------------------------------------------    
! not used
!       ibcon=mibcon(ngr)
!--(DMK-BRAMS-5.0-FIM)------------------------------------------------------------    
       ! write(6,'10I5') ngr,mynum,m1,m2,m3,ia+mi0(ngr),iz+mi0(ngr)&
            !                 ,ja+mj0(ngr),jz+mj0(ngr),ibcon

       ijk=0;cc=0
       !-srf  DO k1=1,m1
       !DO k1=2,m1-1
       DO j1=ja,jz
          DO i1=ia,iz
             DO k1=2,m1-1
                ijk=ijk+1
             END DO
          END DO
       END DO
       !PRINT *,'LFR->Total points:',ijk, mynum 
       !call flush(6)
       !Calculate the end point of each block. Depends on max block size.
       IF(block_size>=ijk) THEN
          nob(ngr)=1
          ALLOCATE(index_g(ngr)%block_end(1))
          index_g(ngr)%block_end(1)=ijk
       ELSE
          nob(ngr)=ijk/block_size
          ALLOCATE(index_g(ngr)%block_end(nob(ngr)))
          index_g(ngr)%block_end=block_size
       END IF
       !Calculate the remaining points in the athmosphere
       resto=ijk-(nob(ngr)*block_size)
       !PRINT *,'LFR->Resto,nob(ngr):',resto,nob(ngr), mynum 
       !call flush(6)
       DO WHILE (resto/=0)
          DO i=1,nob(ngr)
             index_g(ngr)%block_end(i)=index_g(ngr)%block_end(i)+1
             resto=resto-1
             IF(resto==0) EXIT
          END DO
       END DO


       maxblock_size=MAXVAL(index_g(ngr)%block_end)

       !Alocando somente a quantidade necessaria de colunas
       ALLOCATE(index_g(ngr)%indexi   (maxblock_size,nob(ngr)))
       ALLOCATE(index_g(ngr)%indexj   (maxblock_size,nob(ngr)))
       ALLOCATE(index_g(ngr)%indexk   (maxblock_size,nob(ngr)))
       ALLOCATE(index_g(ngr)%ijkindex (1:m1,ia:iz,ja:jz)      )
       ALLOCATE(index_g(ngr)%bindex   (1:m1,ia:iz,ja:jz)      )
       ALLOCATE(index_g(ngr)%kij_index(maxblock_size,nob(ngr)))
       ALLOCATE(index_g(ngr)%last_accepted_dt(nob(ngr)))

       index_g(ngr)%ijkindex=0
       ijk=0
       inb=1
       !srf DO k1=1,m1  
       !DO k1=2,m1-1
       DO j1=ja,jz
          DO i1=ia,iz
             DO k1=2,m1-1
                ijk=ijk+1
                IF(ijk>index_g(ngr)%block_end(inb)) THEN
                   ijk=1
                   inb=inb+1
                END IF
                !indexing (mapping) one position in others
                index_g(ngr)%indexi(ijk,inb)=i1
                index_g(ngr)%indexj(ijk,inb)=j1
                index_g(ngr)%indexk(ijk,inb)=k1
                index_g(ngr)%ijkindex(k1,i1,j1)=ijk
                index_g(ngr)%bindex(k1,i1,j1)=inb
             END DO
          END DO
       END DO
       ! index of tendency array
       !  a(k,i,j) = a(k +     m1*(i-1) +  	   (m1*(m2-1)+m1)*(j-1))
       DO i=1,nob(ngr)
          DO ijk=1,maxblock_size
             index_g(ngr)%kij_index(ijk,i) =                 index_g(ngr)%indexk(ijk,i)   &! k
                  +  m1*           (index_g(ngr)%indexi(ijk,i)-1)&! m1*(i-1)
                  + (m1*(m2-1)+m1)*(index_g(ngr)%indexj(ijk,i)-1) ! m1*(m2-1)+m1)*(j-1)
          ENDDO
       ENDDO
       !- timestep initialization for the timestep control procedure
       !- no futuro incluir no history file para reproducibilidade quanto reiniciando a partir do modo "history"
       index_g(ngr)%last_accepted_dt(1:nob(ngr)) = DBLE(dtlongn(ngr)*N_DYN_CHEM) 

    END DO ! ngrids

    !- calculates the maximum block size (including all grids) to allocate memory 
    !- for spack arrays (see mem_spack.f90)
    maxblock_size = 0
    DO ngr=1,ngrids
       maxblock_size=MAX(maxblock_size, MAXVAL(index_g(ngr)%block_end))
    END DO

  END SUBROUTINE AllocIndex


END MODULE Spack_utils
