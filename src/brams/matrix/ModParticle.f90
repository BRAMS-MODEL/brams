MODULE ModParticle
   private
   ! Particles data
   TYPE part
      REAL,POINTER,DIMENSION(:,:,:,:) :: radius !meter
      REAL,POINTER,DIMENSION(:,:,:,:) :: dens !kg/m^3
   END TYPE part
   TYPE(part), ALLOCATABLE,DIMENSION(:) :: particle
   public :: particle,  copyDensRadius,  allocateParticle
   CONTAIns
   
      SUBROUTINE copyDensRadius(ia,iz,ja,jz,m1,ngrid,ac,dens_mode_dry,nmodes)
          INTEGER,INTENT(IN) :: nmodes,m1,ia,iz,ja,jz,ngrid             
          DOUBLE PRECISION,INTENT(IN) :: dens_mode_dry(m1,ia:iz,ja:jz,nmodes) ! average mode density calculated from component concentrations [g/cm^3]
          DOUBLE PRECISION,INTENT(IN) :: ac      (m1,ia:iz,ja:jz,nmodes) ! minimum dry radius to activate for each mode [um]. 
          INTEGER :: i,j,k,nm
   
          DO nm=1,nmodes
             DO j=ja,jz
               DO i=ia,iz
                  DO k=1,m1
                     particle(ngrid)%radius(k,i,j,nm)=real(ac(k,i,j,nm)*1.0e-6)
                     particle(ngrid)%dens  (k,i,j,nm)=real(dens_mode_dry(k,i,j,nm)*1.0e+3)
                  END DO
               END DO
             END DO
          END DO

      END SUBROUTINE copyDensRadius


      SUBROUTINE allocateParticle(m1,m2,m3,ngrids)
         INTEGER,INTENT(IN) :: ngrids
         INTEGER,INTENT(IN) :: m1(ngrids)
         INTEGER,INTENT(IN) :: m2(ngrids)
         INTEGER,INTENT(IN) :: m3(ngrids)
      
         INTEGER :: n
      
         ALLOCATE(particle(ngrids))
         DO n=1,ngrids
            ALLOCATE(particle(n)%radius(m1(n),m2(n),m3(n),nmodes))
            ALLOCATE(particle(n)%dens(m1(n),m2(n),m3(n),nmodes))
         END DO
      
      END SUBROUTINE allocateParticle
   
end module ModParticle