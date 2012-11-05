SUBROUTINE RESHAPE_MAT2(MV,n,m,V)
USE global ! module for global variables

INTEGER :: i, j, k, n, m
DOUBLE PRECISION, DIMENSION(n*m) :: V
DOUBLE PRECISION, DIMENSION(n,m) :: MV

DO j=1, m
   DO i=1, n
      k = (j-1)*n + i
      V(k) = MV(i,j)
   END DO
END DO

 
END SUBROUTINE RESHAPE_MAT2
