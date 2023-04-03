MODULE mod_chem_trans_liq

  USE mem_chem1   , ONLY: &
       chem1_vars            ! Type

  USE mem_chem1aq , ONLY: &
       chem1aq_vars          ! Type

           
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: trans_liq ! Subroutine
       
CONTAINS

  !======================================================================
  SUBROUTINE trans_liq(m1,m2,m3,ia,iz,ja,jz,chem1_g,chem1aq_g, &
                       nspeciesaq,ind_gas,coll,sedimr,rcp,rrp)

    !====================================================================
    !
    !  author M. Pirre 18/08/2008 
    !
    !==================================================================

    INTEGER , INTENT(IN) :: m1
    INTEGER , INTENT(IN) :: m2
    INTEGER , INTENT(IN) :: m3
    INTEGER , INTENT(IN) :: ia
    INTEGER , INTENT(IN) :: iz
    INTEGER , INTENT(IN) :: ja
    INTEGER , INTENT(IN) :: jz

    ! mem_chem1
    TYPE (chem1_vars)   , INTENT(INOUT) :: chem1_g(nspeciesaq)

    ! mem_chem1aq
    TYPE (chem1aq_vars) , INTENT(INOUT) :: chem1aq_g(nspeciesaq)

    ! chem1aq_list
    INTEGER , INTENT(IN) :: nspeciesaq
    INTEGER , INTENT(IN) :: ind_gas(nspeciesaq)

    ! mem_chemic

    REAL    , INTENT(IN) :: coll  (m1,m2,m3)
    REAL    , INTENT(IN) :: sedimr(m1,m2,m3)

    ! mem_micro
    REAL    , INTENT(IN) :: rcp(m1,m2,m3)
    REAL    , INTENT(IN) :: rrp(m1,m2,m3)

    INTEGER iaq
    !
    !
    DO iaq = 1, nspeciesaq
       CALL trans_l (m1,m2,m3,iaq,ia,iz,ja,jz	 &
            ,chem1aq_g(iaq)%sc_pc	 &
            ,chem1aq_g(iaq)%sc_pr	 &
            ,chem1_g(ind_gas(iaq))%sc_p  &
            ,coll		 &
            ,sedimr		 &
            ,rcp	         &
            ,rrp )	 

    ENDDO

  END SUBROUTINE trans_liq

  !--------------------------------------------------------------
  !
  SUBROUTINE trans_l(m1,m2,m3,iaq,ia,iz,ja,jz,sc_pc,sc_pr,sc_p &
       ,coll,sedimr,rcp,rrp)

    INTEGER , INTENT(IN)    :: m1
    INTEGER , INTENT(IN)    :: m2
    INTEGER , INTENT(IN)    :: m3
    INTEGER , INTENT(IN)    :: iaq
    INTEGER , INTENT(IN)    :: ia
    INTEGER , INTENT(IN)    :: iz
    INTEGER , INTENT(IN)    :: ja
    INTEGER , INTENT(IN)    :: jz
    REAL    , INTENT(INOUT) :: sc_pc (m1,m2,m3)
    REAL    , INTENT(INOUT) :: sc_pr (m1,m2,m3)
    REAL    , INTENT(INOUT) :: sc_p  (m1,m2,m3)
    REAL    , INTENT(IN)    :: coll  (m1,m2,m3)
    REAL    , INTENT(IN)    :: sedimr(m1,m2,m3)
    REAL    , INTENT(IN)    :: rrp   (m1,m2,m3)
    REAL    , INTENT(IN)    :: rcp   (m1,m2,m3)

    INTEGER i,j,k

    DO j = ja , jz
       DO i = ia , iz
          DO k = 2 , m1

             sc_pr(k,i,j)=sc_pr(k,i,j)-sedimr(k,i,j)*sc_pr(k,i,j)   &
                  +coll  (k,i,j)*sc_pc(k,i,j)
             sc_pc(k,i,j)=sc_pc(k,i,j)-coll(k,i,j)*sc_pc(k,i,j)


             IF(rcp(k,i,j) .LE. 1.e-9)THEN
                sc_p (k,i,j)=sc_p(k,i,j)+sc_pc(k,i,j)
                sc_pc(k,i,j)=0.
             ENDIF

             IF(rrp(k,i,j) .LE. 1.e-9)THEN
                sc_p (k,i,j)= sc_p(k,i,j)+sc_pr(k,i,j)
                sc_pr(k,i,j)= 0.
             ENDIF

          ENDDO
       ENDDO
    ENDDO
    !
    !
    ! DO j = ja , jz
    !  DO i = ia , iz
    !    DO k = 2 , m1
    !       if(rcp(k,i,j) .le. 1.e-9)then
    !            sc_p (k,i,j)=sc_p(k,i,j)+sc_pc(k,i,j)
    !            sc_pc(k,i,j)=0.
    !       endif
    !     enddo
    !  enddo
    ! enddo
    !
    ! DO j = ja , jz
    !   DO i = ia , iz
    !     DO k = 2 , m1
    !      if(rrp(k,i,j) .le. 1.e-9)then
    !	sc_p (k,i,j)= sc_p(k,i,j)+sc_pr(k,i,j)
    !	sc_pr(k,i,j)= 0.
    !      endif
    !    enddo
    !  enddo
    ! enddo

    !     DO j = ja , jz
    !       DO i = ia , iz
    !         sc_pr(1,i,j)=0.
    !       enddo
    !     enddo
    !
  END SUBROUTINE trans_l

  !
  !-----------------------------------------------------------------
END MODULE mod_chem_trans_liq
