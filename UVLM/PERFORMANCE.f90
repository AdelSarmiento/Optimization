SUBROUTINE PERFORMANCE
USE global ! module for global variables



INTEGER :: i



  ! Average lift
  CALL TRAPZ(t(N_Avg_Begin:N_steps),CL(N_Avg_Begin:N_steps),N_steps-N_Avg_Begin+1,CL_AVG)
  CL_AVG = CL_AVG/(t(N_steps)-t(N_Avg_Begin)) 
  ! Average thrust
  CALL TRAPZ(t(N_Avg_Begin:N_steps),-CD(N_Avg_Begin:N_steps),N_steps-N_Avg_Begin+1,CT_AVG)
  CT_AVG = CT_AVG/(t(N_steps)-t(N_Avg_Begin)) 
  ! Power input
  CALL TRAPZ(t(N_Avg_Begin:N_steps),POWER(N_Avg_Begin:N_steps)/(0.5D0*rho*V_ref**3),N_steps-N_Avg_Begin+1,CP_IN_AVG)
  CP_IN_AVG = CP_IN_AVG/(t(N_steps)-t(N_Avg_Begin)) 
 
  
 ! Efficiency
 EFF = ABS(CT_AVG/CP_IN_AVG)


END SUBROUTINE PERFORMANCE
