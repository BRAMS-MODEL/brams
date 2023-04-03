MODULE lu
  !This set of routines use a LU decomposition an Lu solver to solver a diferencial equation System
  !
  !Author: Luiz Flavio Rodrigues - CPTEC/INPE - 2007
  !This subroutine do a LU decomposition using a partial pivoting from a quadratic system of equation. 
  !Adapted and based in a sample routine from Numerical Recipes.
  !References:
  !Numerical Recipes in Fortran 77, The Art of Scientific Computing. Cambridge University Press. 1986-1992
  !Ralston, A., and Rabinowitz, P. 1978, A First Course in Numerical Analysis, 2nd ed. (New York:
  !McGraw-Hill), #9.3-1.
  !Isaacson, E., and Keller, H.B. 1966, Analysis of Numerical Methods (New York: Wiley), #2.1.
  !Johnson, L.W., and Riess, R.D. 1982, Numerical Analysis, 2nd ed. (Reading, MA: Addison-
  !Wesley), #2.2.1.
  !Westlake, J.R. 1968, A Handbook of Numerical Matrix Inversion and Solution of Linear Equations
  !(New York: Wiley).

  INTEGER,ALLOCATABLE,DIMENSION(:) :: indx


  CONTAINS
  
  !==================================
  SUBROUTINE Solve(coef,indep,n,keep)
  !==================================
  !Suroutine Solve. Solve a system o differential equations
  IMPLICIT NONE
  
  INTEGER,INTENT(IN)              :: n
  DOUBLE PRECISION,DIMENSION(n,n),INTENT(IN)  :: coef  !Matrix coeficients to solve
  DOUBLE PRECISION,DIMENSION(n),INTENT(INOUT) :: indep !Independent terms at input and solutions at output
  LOGICAL,INTENT(IN),OPTIONAL     :: keep  !To mantain the allocate of indx
  
  DOUBLE PRECISION,DIMENSION(n,n) :: a
  
  a=coef !Copy a input coeficients to manipulate array
  
  IF(.NOT. allocated(indx)) CALL Prepare_Lu(n) ! Allocated the indexarray
  
  CALL ludecomp(a,n)
  
  CALL lusolver(a,n,indep)
  
  IF(PRESENT(keep)) THEN
    IF(.NOT. keep) CALL Destroy_lu()
  ELSE
    CALL Destroy_lu()
  END IF
  
  END SUBROUTINE Solve
  
  
  !========================
  SUBROUTINE Prepare_Lu(n)
  !========================
  !Subroutine Prepare_LU - Allocate a transfer index between decomp and solver
  IMPLICIT NONE
  
  INTEGER,INTENT(IN) :: n
  
  IF(.NOT. ALLOCATED(indx)) ALLOCATE(indx(n))
  !PRINT *, 'Passei'
  
  END SUBROUTINE Prepare_LU
  
  
  !========================
  SUBROUTINE Destroy_Lu()
  !========================
  !Subroutine to destroy a index
  
  IF(allocated(indx)) DEALLOCATE(indx)
  
  END SUBROUTINE Destroy_Lu
  
  !=========================
  SUBROUTINE ludecomp (a, n)
  !=========================  
  !Subroutine ludecomp -  LU decomposition from a
  !                         L.U = A      
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: n                    !Size of quadratic system
  DOUBLE PRECISION,DIMENSION(n,n),INTENT(INOUT) :: a      !Coeficients Terms of equation (input)
					      !LU decomposition of matrix (output)
  !INTEGER,DIMENSION(N),INTENT(OUT) :: indx    !Vector with index of row permutation
  DOUBLE PRECISION,PARAMETER :: zero=1.0D-30

  DOUBLE PRECISION,DIMENSION(n) :: vv
  INTEGER :: i, imax, j, k  
  DOUBLE PRECISION :: aamax, dum, temp

  do i = 1, n  
    aamax = 0.  
    do j = 1, n  
      if (abs(a(i,j))>aamax) aamax=abs(a(i,j))  
    end do  
    if (aamax==0.) THEN
          Write (*,*) 'stopping: singular matrix in ludecomp'
	  stop
    END IF
    vv(i) = 1./aamax  
  end do  

  do j = 1, n  
    do i = 1, j - 1  
      temp = a (i, j)  
      do k = 1, i - 1  
        temp = temp - a (i, k) * a (k, j)  
      end do  
      a (i, j) = temp  
    end do  

    aamax = 0.  

    do i = j, n  
      temp = a (i, j)  
      do k = 1, j - 1  
       temp = temp - a (i, k) * a (k, j)  
      end do  

      a (i, j) = temp  

      dum = vv (i) * abs (temp)  
      if (dum>=aamax) then  
        imax = i  
        aamax = dum  
      end if  

    end do  

    if (j/=imax) then  
      do k = 1, n  
        dum = a (imax, k)  
        a (imax, k) = a (j, k)  
        a (j, k) = dum  
      end do  
      vv (imax) = vv (j)  
    end if  
    indx (j) = imax  

    if (a (j, j)==0.) a (j, j) = zero  
    if (j/=n) then  
      dum = 1. / a (j, j)  
      do i = j + 1, n  
        a (i, j) = a (i, j) * dum  
      end do  
    end if  
  end do  

  end subroutine ludecomp

  !===========================
  SUBROUTINE lusolver (a,n,b) 
  !===========================  
  !Subroutine lusolver -  forward substitution and backsubstitution
  !                   A . x = (L . U) . x = L . (U . x) = b
  IMPLICIT NONE
  INTEGER,INTENT(IN)              :: n           !Size of quadratic system
  !INTEGER,DIMENSION(n),INTENT(IN) :: indx        !Index for pivot obtained from LU decomposition

  DOUBLE PRECISION,DIMENSION(n,n),INTENT(IN)  :: a           !Array obtained from LU decomposition
  DOUBLE PRECISION,DIMENSION(n),INTENT(INOUT) :: b           !b - Independent terms of equations (in input)
                                                 !    Solution of system (in output)
  
  INTEGER :: i, ii, j, ll  

  DOUBLE PRECISION :: temp 

  ii = 0  

  do i = 1, n  

    ll = indx (i)  
    temp= b (ll)  
    b (ll) = b (i)  
    if (ii/=0) then  
       do j = ii, i - 1  
         temp= temp- a (i, j) * b (j)  
       end do  
    else if (temp/=0.) then  
       ii = i  
    end if  
    b (i) = temp 
  end do  
  do i = n, 1, - 1  
    temp= b (i)  
    do j = i + 1, n  
      temp= temp- a (i, j) * b (j)  
    end do  
    b (i) = temp/ a (i, i)  
  end do  

end subroutine lusolver


END MODULE lu
