SUBROUTINE ANGLES
USE global ! module for global variables
IMPLICIT NONE

INTEGER :: i, j
DOUBLE PRECISION, DIMENSION(3) :: V1, V2, V3, V4
DOUBLE PRECISION :: ALPHAM
DOUBLE PRECISION :: AMPV1, AMPV2, AMPV3, AMPV4
DOUBLE PRECISION :: ALPHA1, ALPHA3

ALPHA_MAX = 0.0D0

DO i=1, NR
   DO j=1, NC
      ! 1st vector
      V1(1) = X0(i,j+1)-X0(i,j)
      V1(2) = Y0(i,j+1)-Y0(i,j) 
      V1(3) = Z0(i,j+1)-Z0(i,j)
      ! 2nd vector
      V2(1) = X0(i+1,j+1)-X0(i,j+1)
      V2(2) = Y0(i+1,j+1)-Y0(i,j+1) 
      V2(3) = Z0(i+1,j+1)-Z0(i,j+1)
	  ! 3rd vector
      V3(1) = X0(i+1,j)-X0(i+1,j+1)
      V3(2) = Y0(i+1,j)-Y0(i+1,j+1) 
      V3(3) = Z0(i+1,j)-Z0(i+1,j+1)
	  ! 4th vector
      V4(1) = X0(i,j)-X0(i+1,j)
      V4(2) = Y0(i,j)-Y0(i+1,j) 
      V4(3) = Z0(i,j)-Z0(i+1,j)

      AMPV1 = SQRT(V1(1)**2+V1(2)**2+V1(3)**2)
	  AMPV2 = SQRT(V2(1)**2+V2(2)**2+V2(3)**2)
	  AMPV3 = SQRT(V3(1)**2+V3(2)**2+V3(3)**2)
	  AMPV4 = SQRT(V4(1)**2+V4(2)**2+V4(3)**2)
	  
	  ALPHA1 = ACOS((V1(1)*V4(1)+V1(2)*V4(2)+V1(3)*V4(3))/(AMPV1*AMPV4))
	  ALPHA3 = ACOS((V2(1)*V3(1)+V2(2)*V3(2)+V2(3)*V3(3))/(AMPV2*AMPV3))
	  !WRITE(*,*) ALPHA1, ALPHA3
	  !WRITE(*,*)
	  ALPHAM = MAX(ABS(ALPHA1-PI/2.0D0),ABS(ALPHA3-PI/2.0D0))
	  
	  IF (ALPHAM .GT. ALPHA_MAX) THEN
	  ALPHA_MAX = ALPHAM
	  END IF

   END DO
END DO

END SUBROUTINE ANGLES
