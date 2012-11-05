SUBROUTINE VEL_KIN
USE global ! module for global variables

INTEGER :: i, j
! Compute velocities due to wing motion (finite difference)

DO i=1, NR
   DO j=1, NC
      IF(i_time .EQ. 1) THEN
         U_KIN(i,j,i_time) = 0.0D0
         V_KIN(i,j,i_time) = 0.0D0
         W_KIN(i,j,i_time) = 0.0D0
      ELSE
         U_KIN(i,j,i_time) = -( XC(i,j,i_time)-XC(i,j,i_time-1) )/t_step
         V_KIN(i,j,i_time) = -( YC(i,j,i_time)-YC(i,j,i_time-1) )/t_step
         W_KIN(i,j,i_time) = -( ZC(i,j,i_time)-ZC(i,j,i_time-1) )/t_step
      END IF
   END DO
END DO

END SUBROUTINE VEL_KIN
