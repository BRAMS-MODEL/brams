MODULE uv_atten

  IMPLICIT NONE

  INTEGER , PARAMETER :: yes=1,no=0

  INTEGER :: uv_atten_initialized=0

  TYPE uv_atten_vars 
     REAL,POINTER,DIMENSION  (:,:)  :: att 
  END TYPE uv_atten_vars

  TYPE(uv_atten_vars), ALLOCATABLE,DIMENSION(:)   :: uv_atten_g

  INTEGER, PARAMETER :: Length=10
  INTEGER, PARAMETER :: MaxCity=81
  CHARACTER(Length), DIMENSION (MaxCity) :: CPrefix 
  DATA CPrefix/ &
       'BUENOSAI', &
       'santiago', &
       'RIOGRAND', &
       'BAGEMERS', &
       'CAXIASRS', &
       'SANTOANG', &
       'URUGUAIA', &
       'PORTOALE', &
       'SANTAMAR', &
       'LAGESERS', &
       'CHAPECRS', &
       'FLORIPSC', &
       'POMERODE', &
       'JOIVILLE', &
       'MARINGPR', &
       'CASCAVEL', &
       'FOZIGUAC', &
       'GUARAPVA', &
       'CURITIBA', &
       'RIBEIRAP', &
       'SERTAOZI', &
       'BIRIGEPR', &
       'CAMPINAS', &
       'SANTOSSP', &
       'BARRETOS', &
       'SAOJOSER', &
       'SAOJOSEC', &
       'SAOPAULO', &
       'ITAQQTUB', &
       'SOROCABA', &
       'CABOFRRJ', &
       'NOFBURGO', &
       'TEREZOPO', &
       'RIODEJAN', &
       'PETROPOL', &
       'UBERABAA', &
       'UBERLAND', &
       'ARAXAERJ', &
       'MOCLAROS', &
       'BELOHORZ', &
       'VARGINHA', &
       'GOVVALAD', &
       'TEOFOTON', &
       'CAMPOGMS', &
       'DOURADOS', &
       'RONDONOP', &
       'CUIABAMT', &
       'VARZEAGR', &
       'SERRAEES', &
       'VILAVELH', &
       'VITORIAA', &
       'GUARAPAR', &
       'LINHARES', &
       'BRASILDF', &
       'GOIANIGO', &
       'CATALAGO', &
       'RIOVERDE', &
       'PORTOVEL', &
       'VILHENAA', &
       'BARREIRA', &
       'JEQUIEBA', &
       'SALVADOR', &
       'ILHEUSBA', &
       'VITORIAC', &
       'ARACAJSE', &
       'MACEIOAL', &
       'PETROLIN', &
       'RECIFEPE', &
       'OLINDAPE', &
       'JAUARAPS', &
       'JOAOPESS', &
       'NATALERN', &
       'FORTALEZ', &
       'TERESINA', &
       'SAOLUIZZ', &
       'SANTAREM', &
       'BELEMEPA', &
       'MACAPAAP', &
       'MANAUSAM', &
       'BOAVISTA', &
       'RIOBRANC'/


  REAL,DIMENSION(MaxCity) ::alat
  DATA alat/ &
       -34.58,&  !Buenos Aires      
       -33.47,&  !santiago		
       -32.00,&  !Rio Grande        
       -31.00,&  !Bage	       
       -29.00,&  !Caxias do Sul     
       -28.00,&  !Santo Angelo    
       -30.00,&  !Uruguaiana        
       -30.03,&  !Porto Alegre      
       -30.00,&  !Santa Maria       
       -28.00,&  !Lages	       
       -27.00,&  !Chapeco	       
       -28.00,&  !Florianopolis     
       -27.00,&  !Pomerode	       
       -26.00,&  !Joinville         
       -23.00,&  !Maringa	       
       -25.00,&  !Cascavel	       
       -26.00,&  !Foz do Iguacu     
       -25.00,&  !Guarapuava        
       -25.43,&  !Curitiba	       
       -21.00,&  !Ribeirao Preto    
       -21.00,&  !Sertaozinho       
       -21.00,&  !Birigui	       
       -22.54,&  !Campinas	       
       -24.00,&  !Santos	       
       -21.00,&  !Barretos	       
       -21.00,&  !SaoJosedoRioPreto 
       -23.10,&  !SaoJosedosCampos  
       -23.62,&  !Sao Paulo         
       -23.00,&  !Itaquaquecetuba   
       -23.30,&  !Sorocaba	       
       -23.00,&  !Cabo Frio         
       -22.00,&  !Nova Friburgo     
       -22.00,&  !Teresopolis       
       -22.90,&  !Rio de Janeiro    
       -23.00,&  !Petropolis        
       -20.00,&  !Uberaba	       
       -19.00,&  !Uberlandia        
       -20.00,&  !Araxa	       
       -17.00,&  !Montes Claros     
       -19.91,&  !Belo Horizonte    
       -22.00,&  !Varginha	       
       -19.00,&  !Governador Valad. 
       -18.00,&  !Teofilo Otoni     
       -20.44,&  !Campo Grande      
       -22.00,&  !Dourados	       
       -16.00,&  !Rondonopolis      
       -15.59,&  !Cuiaba	       
       -16.00,&  !Varzea Grande     
       -20.00,&  !Serra	       
       -20.00,&  !Vila Velha        
       -20.00,&  !Vitoria	       
       -21.00,&  !Guarapari         
       -19.00,&  !Linhares	       
       -15.78,&  !Brasilia	       
       -16.68,&  !Goiania	       
       -18.00,&  !Catalao	       
       -18.00,&  !Rio Verde         
       -08.71,&  !Porto Velho       
       -13.00,&  !Vilhena	       
       -12.00,&  !Barreiras         
       -14.00,&  !Jequiao	       
       -13.00,&  !Salvador	       
       -15.00,&  !Ilheus	       
       -15.00,&  !VitoriadaConquista
       -11.00,&  !Aracaju	       
       -10.00,&  !Maceio	       
       -09.00,&  !Petrolina         
       -08.00,&  !Recife	       
       -08.00,&  !Olinda	       
       -08.00,&  !Jaboatao dos Guar 
       -07.00,&  !Joao Pessoa       
       -06.00,&  !Natal	       
       -04.00,&  !Fortaleza         
       -05.00,&  !Teresina	       
       -03.00,&  !Sao Luiz	       
       -02.00,&  !Santarem	       
       -01.00,&  !Belem	       
       00.00,&   !Macapa	       
       -03.06,&  !Manaus	       
       -03.00,&  !Boa Vista         
       -09.90/   !Rio Branco        


  REAL,DIMENSION(MaxCity) ::alon
  DATA alon/ &
       -58.67,&    !Buenos Aires	 
       -70.75,&    !santiago 	 
       -52.00,&    !Rio Grande	 
       -54.00,&    !Bage		 
       -51.00,&    !Caxias do Sul	 
       -54.00,&    !Santo Angelo	     
       -57.00,&    !Uruguaiana	 
       -51.23,&    !Porto Alegre	 
       -54.00,&    !Santa Maria	 
       -50.00,&    !Lages		 
       -53.00,&    !Chapeco  	 
       -49.00,&    !Florianopolis	 
       -49.00,&    !Pomerode 	 
       -59.00,&    !Joinville	 
       -52.00,&    !Maringa  	 
       -54.00,&    !Cascavel 	 
       -55.00,&    !Foz do Iguacu	 
       -51.00,&    !Guarapuava	 
       -49.27,&    !Curitiba 	 
       -48.00,&    !Ribeirao Preto	 
       -48.00,&    !Sertaozinho	 
       -50.00,&    !Birigui  	 
       -47.03,&    !Campinas 	 
       -46.00,&    !Santos		 
       -49.00,&    !Barretos 	 
       -49.00,&    !SaoJosedoRioPreto 
       -45.53,&    !SaoJosedosCampos  
       -46.65,&    !Sao Paulo	 
       -46.00,&    !Itaquaquecetuba   
       -47.27,&    !Sorocaba 	 
       -42.00,&    !Cabo Frio	 
       -43.00,&    !Nova Friburgo	 
       -43.00,&    !Teresopolis	 
       -43.21,&    !Rio de Janeiro	 
       -43.00,&    !Petropolis	 
       -48.00,&    !Uberaba  	 
       -48.00,&    !Uberlandia	 
       -47.00,&    !Araxa		 
       -44.00,&    !Montes Claros	 
       -43.93,&    !Belo Horizonte	 
       -45.00,&    !Varginha 	 
       -42.00,&    !Governador Valad. 
       -42.00,&    !Teofilo Otoni	 
       -55.65,&    !Campo Grande	 
       -55.00,&    !Dourados 	 
       -55.00,&    !Rondonopolis	 
       -56.09,&    !Cuiaba		 
       -56.00,&    !Varzea Grande	 
       -40.00,&    !Serra		 
       -40.00,&    !Vila Velha	 
       -40.00,&    !Vitoria  	 
       -40.00,&    !Guarapari	 
       -40.00,&    !Linhares 	 
       -47.93,&    !Brasilia 	 
       -49.26,&    !Goiania  	 
       -48.00,&    !Catalao  	 
       -51.00,&    !Rio Verde	 
       -63.90,&    !Porto Velho	 
       -60.00,&    !Vilhena  	 
       -45.00,&    !Barreiras	 
       -40.00,&    !Jequiao  	 
       -39.00,&    !Salvador 	 
       -39.00,&    !Ilheus		 
       -41.00,&    !VitoriadaConquista
       -37.00,&    !Aracaju  	 
       -36.00,&    !Maceio		 
       -41.00,&    !Petrolina	 
       -35.00,&    !Recife		 
       -35.00,&    !Olinda		 
       -35.00,&    !Jaboatao dos Guar 
       -35.00,&    !Joao Pessoa	 
       -35.00,&    !Natal		 
       -39.00,&    !Fortaleza	 
       -43.00,&    !Teresina 	 
       -44.00,&    !Sao Luiz 	 
       -55.00,&    !Santarem 	 
       -49.00,&    !Belem		 
       -51.00,&    !Macapa		 
       -59.98,&    !Manaus		 
       -61.00,&    !Boa Vista	 
       -67.80/     !Rio Branco	 

  REAL,    DIMENSION(MaxCity) :: city_size
  DATA city_size/ &
                                ! area (m^2)  
                                !----------------------
       203.    ,&	  !Buenos Aires      '203.0,& 
       -99.    ,&	  !santiago	      '1536.1,& 
       -99.    ,&	  !Rio Grande	     '2719.4,& 
       -99.    ,&	  !Bage 	     '4655.0,& 
       -99.    ,&	  !Caxias do Sul     !1470.0,&
       -99.    ,&	  !Santo Angelo      !572.3 ,&
       -99.    ,&	  !Uruguaiana	     !5339.5,&
       470.    ,&	  !Porto Alegre      !481.8,&
       -99.    ,&	  !Santa Maria       !1729.4,&
       -99.    ,&	  !Lages	     !2603.1,&
       -99.    ,&	  !Chapeco	     !574.4	  ,&
       -99.    ,&	  !Florianopolis     !437.0	  ,&
       -99.    ,&	  !Pomerode	     !185.1	  ,&
       -99.    ,&	  !Joinville	     !1046.3,&
       -99.    ,&	  !Maringa	     !483.3	  ,&
       -99.    ,&	  !Cascavel	     !1928.1,&
       -99.    ,&	  !Foz do Iguacu     !586.1	  ,&
       -99.    ,&	  !Guarapuava	     !2891.7,&
       -99.    ,&	  !Curitiba	     !431.0	  ,&
       274.1   ,&	  !Ribeirao Preto    !649.3	  ,&
       -99.    ,&	  !Sertaozinho       !386.4	  ,&
       -99.    ,&	  !Birigui	     !513.4	  ,&
       785.5   ,&	  !Campinas	     !785.5	  ,&
       -99.    ,&	  !Santos	     !279.5	  ,&
       -99.    ,&	  !Barretos	     !1493.4,&
       83.46   ,&	  !SaoJosedoRioPreto !409.1	  ,&
       288.    ,&	  !SaoJosedosCampos  !1089.6,&
       1441.4  ,&	  !Sao Paulo	     !1441.4,&   2.209 km ** 2 == Regiao metropolitana de Sao Paulo
       -99.    ,&	  !Itaquaquecetuba   !82.2	  ,&
       249.2   ,&	  !Sorocaba	     !443.0	  ,&
       -99.    ,&	  !Cabo Frio	     !338.2	  ,&
       -99.    ,&	  !Nova Friburgo     !819.0	  ,&
       -99.    ,&	  !Teresopolis       !646.8	  ,&
       1265.8  ,&	  !Rio de Janeiro    !1265.8,&
       -99.    ,&	  !Petropolis	     !734.1	  ,&
       -99.    ,&	  !Uberaba	     !4407.4,&
       135.    ,&	  !Uberlandia	     !4024.5,&
       -99.    ,&	  !Araxa	     !1150.0,&
       -99.    ,&	  !Montes Claros     !3394.2,&
       333.2   ,&	  !Belo Horizonte    !333.2	  ,&
       -99.    ,&	  !Varginha	     !381.6	  ,&
       -99.    ,&	  !Governador Valad. !2252.8,&
       -99.    ,&	  !Teofilo Otoni     !2589.7,&
       -99.    ,&	  !Campo Grande      !8038.2,&
       -99.    ,&	  !Dourados	     !3729.6,&
       -99.    ,&	  !Rondonopolis      !3950.9,&
       -99.    ,&	  !Cuiaba	     !3935.0,&
       -99.    ,&	  !Varzea Grande     !889.7	  ,&
       -99.    ,&	  !Serra	     !547.8	  ,&
       -99.    ,&	  !Vila Velha	     !218.6	  ,&
       -99.    ,&	  !Vitoria	     !89.1	  ,&
       -99.    ,&	  !Guarapari	     !544.8	  ,&
       -99.    ,&	  !Linhares	     !2859.0,&
       -99.    ,&	  !Brasilia	     !5589.9,&
       -99.    ,&	  !Goiania	     !739.8	  ,&
       -99.    ,&	  !Catalao	     !3388.6,&
       -99.    ,&	  !Rio Verde	     !7631.6,&
       -99.    ,&	  !Porto Velho       !27929.5,&
       -99.    ,&	  !Vilhena	     !10766.2,&
       -99.    ,&	  !Barreiras	     !10525.8,&
       -99.    ,&	  !Jequiao	     !2697.6,&
       -99.    ,&	  !Salvador	     !709.9 ,&
       -99.    ,&	  !Ilheus	     !1349.9,&
       -99.    ,&	  !VitoriadaConquista!2764.0,&
       -99.    ,&	  !Aracaju	     !182.1	  ,&
       -99.    ,&	  !Maceio	     !512.1	  ,&
       -99.    ,&	  !Petrolina	     !3622.6,&
       -99.    ,&	  !Recife	     !218.8	  ,&
       -99.    ,&	  !Olinda	     !37.3	  ,&
       -99.    ,&	  !Jaboatao dos Guar !252.0	  ,&
       -99.    ,&	  !Joao Pessoa       !211.9	  ,&
       -99.    ,&	  !Natal	     !170.6	  ,&
       -99.    ,&	  !Fortaleza	     !314.3	  ,&
       -99.    ,&	  !Teresina	     !1593.3,&
       -99.    ,&	  !Sao Luiz	     !984.9	  ,&
       -99.    ,&	  !Santarem	     !17249.7,&
       -99.    ,&	  !Belem	     !1063.8,&
       -99.    ,&	  !Macapa	     !6279.1,&
       377.    ,&	  !Manaus	     !11402.2,&
       -99.    ,&	  !Boa Vista	     !23.1	  ,&
       -99.    /	  !Rio Branco	     !8909.4/


  REAL,    DIMENSION(MaxCity) :: city_local_att
  DATA city_local_att/ & !   Local attenuation due local pollution
       0.8	 ,&	 !Buenos Aires    
       0.8	 ,&	 !SANTIAGO	  
       1.0	 ,&	 !Rio Grande	  
       1.0	 ,&	 !Bage  	  
       1.0	 ,&	  !Caxias do Sul     
       1.0	 ,&	  !Santo Angelo      
       1.0	 ,&	  !Uruguaiana	     
       0.9	 ,&	  !Porto Alegre      
       1.0	 ,&	  !Santa Maria       
       1.0	 ,&	  !Lages	     
       1.0	 ,&	  !Chapeco	     
       1.0	 ,&	  !Florianopolis     
       1.0	 ,&	  !Pomerode	     
       1.0	 ,&	  !Joinville	     
       1.0	 ,&	  !Maringa	     
       1.0	 ,&	  !Cascavel	     
       1.0	  ,&	  !Foz do Iguacu     
       1.0	  ,&	  !Guarapuava	     
       1.0	  ,&	  !Curitiba	     
       1.0	  ,&	  !Ribeirao Preto    
       1.0	  ,&	  !Sertaozinho       
       1.0	  ,&	  !Birigui	     
       0.95    ,&	  !Campinas	     
       1.0	  ,&	  !Santos	     
       1.0	  ,&	  !Barretos	     
       1.0	  ,&	  !SaoJosedoRioPreto 
       0.98    ,&	  !SaoJosedosCampos  
       0.8	  ,&	  !Sao Paulo	     
       1.0	  ,&	  !Itaquaquecetuba   
       1.0	  ,&	  !Sorocaba	     
       1.0	  ,&	  !Cabo Frio	     
       1.0	  ,&	  !Nova Friburgo     
       1.0	  ,&	  !Teresopolis       
       0.9	  ,&	  !Rio de Janeiro    
       1.0	  ,&	  !Petropolis	     
       1.0	  ,&	  !Uberaba	     
       1.	  ,&	  !Uberlandia	     
       1.0	  ,&	  !Araxa	     
       1.0	  ,&	  !Montes Claros     
       0.9	  ,&	  !Belo Horizonte    
       1.0	  ,&	  !Varginha	     
       1.0	  ,&	  !Governador Valad. 
       1.0	  ,&	  !Teofilo Otoni     
       1.0	  ,&	  !Campo Grande      
       1.0	  ,&	  !Dourados	     
       1.0	  ,&	  !Rondonopolis      
       1.0	  ,&	  !Cuiaba	     
       1.0	  ,&	  !Varzea Grande     
       1.0	  ,&	  !Serra	     
       1.0	  ,&	  !Vila Velha	     
       1.0	  ,&	  !Vitoria	     
       1.0	  ,&	  !Guarapari	     
       1.0	  ,&	  !Linhares	     
       1.0	  ,&	  !Brasilia	     
       1.0	  ,&	  !Goiania	     
       1.0	  ,&	  !Catalao	     
       1.0	  ,&	  !Rio Verde	     
       1.0	  ,&	  !Porto Velho       
       1.0	  ,&	  !Vilhena	     
       1.0	  ,&	  !Barreiras	     
       1.0	  ,&	  !Jequiao	     
       1.0	  ,&	  !Salvador	     
       1.0	  ,&	  !Ilheus	     
       1.0	  ,&	  !VitoriadaConquista
       1.0	  ,&	  !Aracaju	     
       1.0	  ,&	  !Maceio	     
       1.0	  ,&	  !Petrolina	     
       1.0	  ,&	  !Recife	     
       1.0	  ,&	  !Olinda	     
       1.0	  ,&	  !Jaboatao dos Guar 
       1.0	  ,&	  !Joao Pessoa       
       1.0	  ,&	  !Natal	     
       1.0	  ,&	  !Fortaleza	     
       1.0	  ,&	  !Teresina	     
       1.0	  ,&	  !Sao Luiz	     
       1.0	  ,&	  !Santarem	     
       1.0	  ,&	  !Belem	     
       1.0	  ,&	  !Macapa	     
       1.0	  ,&	  !Manaus	     
       1.0	  ,&	  !Boa Vista	     
       1.0	  /	  !Rio Branco	     

  PUBLIC :: uv_attenuation ! Subroutine

CONTAINS
  !----------------------------------------------------

  SUBROUTINE uv_attenuation(m2,m3,ia,iz,ja,jz,i0,j0,   &
                            platn,plonn,deltaxn,   &
                            deltayn,nxpmax,nypmax,xt,yt,  &
                            nwave,aot, ilwrtyp, &
                            iswrtyp,att,i550)

    INTEGER , INTENT(IN) :: m2
    INTEGER , INTENT(IN) :: m3
    INTEGER , INTENT(IN) :: ia
    INTEGER , INTENT(IN) :: iz
    INTEGER , INTENT(IN) :: ja
    INTEGER , INTENT(IN) :: jz
    INTEGER , INTENT(IN) :: i0
    INTEGER , INTENT(IN) :: j0
    integer , intent(in) :: i550
    REAL    , INTENT(IN) :: platn !(ngrids)
    REAL    , INTENT(IN) :: plonn !(ngrids)
    REAL    , INTENT(IN) :: deltaxn !(ngrids)
    REAL    , INTENT(IN) :: deltayn !(ngrids)
    INTEGER , INTENT(IN) :: nxpmax
    INTEGER , INTENT(IN) :: nypmax
    REAL    , INTENT(IN) :: xt(nxpmax)
    REAL    , INTENT(IN) :: yt(nypmax)
    INTEGER , INTENT(IN) :: nwave
    REAL    , INTENT(IN) :: aot(m2,m3,nwave)
    INTEGER , INTENT(IN) :: ilwrtyp
    INTEGER , INTENT(IN) :: iswrtyp
    REAL    , INTENT(INOUT) :: att(m2,m3)

    INTEGER i,j,idaot2,ic,icity,jcity
!!$    INTEGER, DIMENSION (ngrids) :: isp,jsp
    REAL qx,qy,qlon,qlat

  !print*,'821'; call flush(6)	  

!!$    IF(uv_atten_initialized==no) THEN
!!$  !print*,'822'; call flush(6)	  
!!$       CALL initialize_uv_atten(ngrids, mmxp,mmyp)
!!$       uv_atten_initialized=yes
!!$    ENDIF
!!$ !print*,'823'; call flush(6)	  

    !- if not using, CARMA, attenuation will be always 1 (see allocation routine)
    IF (ilwrtyp .NE. 4 .OR. iswrtyp .NE. 4) RETURN
 !print*,'824'; call flush(6)	  

    DO j=ja,jz
       DO i=ia,iz

          !- get attenuation (MFA)

!          idaot2 = MAX(MIN(INT(10*((ANINT(10.*carma(ng)%aot(i,j,11))/10.)+0.2)/3.),7),1)
          idaot2 = MAX(MIN(INT(10*((ANINT(10.*aot(i,j,i550))/10.)+0.2)/3.),7),1)

          IF (idaot2.EQ.1) att(i,j) =0.898
          IF (idaot2.EQ.2) att(i,j) =0.865
          IF (idaot2.EQ.3) att(i,j) =0.832
          IF (idaot2.EQ.4) att(i,j) =0.799
          IF (idaot2.EQ.5) att(i,j) =0.766
          IF (idaot2.EQ.6) att(i,j) =0.733
          IF (idaot2.EQ.7) att(i,j) =0.7
          IF(att(i,j) <  0.7) STOP 'wrong att < 0.7'

       ENDDO
    ENDDO
    !- srf
    !- special treatment for megacities, reduction to account the local poluttion
    !- not included in the aot sources (only bioburn, for now)

    DO ic=1,MaxCity
       !!print*,'city-',cprefix(ic) ; call flush(6)
       qlon=alon(ic)
       qlat=alat(ic)

       CALL ge_to_xy(platn,plonn,qlon,qlat,qx,qy)
       icity =	(NINT((qx-xt(1+i0))/deltaxn)) + 1	     
       jcity =	(NINT((qy-yt(1+j0))/deltayn)) + 1

       IF(icity .LT. 1 .OR. icity .GT. iz)  icity = -999
       IF(jcity .LT. 1 .OR. jcity .GT. jz)  jcity = -999

       IF(jcity == -999 .OR.  icity == -999) CYCLE
       att(icity,jcity) = city_local_att(ic)*att(icity,jcity)
    ENDDO

  END SUBROUTINE uv_attenuation
  !----------------------------------------------------

  SUBROUTINE initialize_uv_atten(ngrids, mmxp,mmyp)

    INTEGER , INTENT(IN) :: ngrids
    INTEGER , INTENT(IN) :: mmxp(ngrids)
    INTEGER , INTENT(IN) :: mmyp(ngrids)

    INTEGER :: ng

    ALLOCATE (uv_atten_g(ngrids))
    DO ng=1,ngrids
       ALLOCATE(uv_atten_g(ng)%att (mmxp(ng),mmyp(ng)) )
       uv_atten_g(ng)%att= 1.0
    ENDDO

  END SUBROUTINE initialize_uv_atten
  !----------------------------------------------------

  SUBROUTINE ge_to_xy(polelat,polelon,xlon,xlat,x,y)

    REAL , INTENT(IN)  :: polelat
    REAL , INTENT(IN)  :: polelon
    REAL , INTENT(IN)  :: xlon
    REAL , INTENT(IN)  :: xlat
    REAL , INTENT(OUT) :: x
    REAL , INTENT(OUT) :: y
    
    REAL, PARAMETER :: rt=6367000.,p=3.14159265360/180.

    REAL f,b

    !transformacao horizontal:
    b = 1.0+SIN(p*xlat)*SIN(p*polelat)+  	       &
         COS(p*xlat)*COS(p*polelat)*COS(p*(xlon-polelon))

    f = 2.00*rt/b

    y = f*(COS(p*polelat)*SIN(p*xlat) -  	       &
         SIN(p*polelat)*COS(p*xlat)*COS(p*(xlon-polelon)))

    x = f*(COS(p*xlat)*SIN(p*(xlon - polelon)))

  END SUBROUTINE ge_to_xy
  !----------------------------------------------------

END MODULE uv_atten
