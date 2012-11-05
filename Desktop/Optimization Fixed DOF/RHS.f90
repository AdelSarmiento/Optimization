SUBROUTINE RHS
USE global ! module for global variables
IMPLICIT NONE

INTEGER :: i,j,k
k=0
DO i=1, NR
   DO j=1, NC
      k = k + 1
      NORM_FLOW(i,j,i_time) = Outward_X(i,j,i_time)*U_KIN(i,j,i_time) & 
	  + Outward_Y(i,j,i_time)*V_KIN(i,j,i_time) + Outward_Z(i,j,i_time)& 
	  *W_KIN(i,j,i_time)
      LL(k,i_time) = - NORM_FLOW(i,j,i_time)
   END DO
END DO

END SUBROUTINE RHS
