SUBROUTINE TRAPZ(X_T,X_P,n,INTE)
USE global ! module for global variables

INTEGER :: n, i

DOUBLE PRECISION, DIMENSION(n) :: X_T, X_P

DOUBLE PRECISION :: INTE

INTE = 0.0D0

DO i=1, n-1
   INTE = INTE + X_P(i)*(X_T(i+1)-X_T(i))
END DO

END SUBROUTINE TRAPZ
