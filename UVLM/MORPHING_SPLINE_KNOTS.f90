SUBROUTINE MORPHING_SPLINE_KNOTS(KNOTS_VAL)

USE global ! module for global variables

INTEGER :: i
DOUBLE PRECISION :: yp1, ypn, OUT_spline
DOUBLE PRECISION :: OUT_spline2, OUT_spline3, OUT_spline4
DOUBLE PRECISION, DIMENSION(14) :: KNOTS_VAL

 DO i= 1, 3*N_K
    t_spline(i) = REAL(i-1)*2.0D0*PI/Omega/N_K + t(N_Avg_Begin) - 2.0D0*PI/Omega! time spline vector
    
 END DO

 TW_spline(1) = KNOTS_VAL(1)
 TW_spline(2) = KNOTS_VAL(2)
 TW_spline(3) = KNOTS_VAL(3)
 TW_spline(4) = KNOTS_VAL(4)
 TW_spline(5) = KNOTS_VAL(5)
 TW_spline(6) = KNOTS_VAL(6)
 TW_spline(7) = KNOTS_VAL(7)

 TW2_spline(1) = 0.0D0
 TW2_spline(2) = 0.0D0
 TW2_spline(3) = 0.0D0
 TW2_spline(4) = 0.0D0
 TW2_spline(5) = 0.0D0
 TW2_spline(6) = 0.0D0
 TW2_spline(7) = 0.0D0

 

 B_spline(1) = KNOTS_VAL(8)
 B_spline(2) = KNOTS_VAL(9)
 B_spline(3) = KNOTS_VAL(10)
 B_spline(4) = KNOTS_VAL(11)
 B_spline(5) = KNOTS_VAL(12)
 B_spline(6) = KNOTS_VAL(13)
 B_spline(7) = KNOTS_VAL(14)

 B2_spline(1) = 0.0D0
 B2_spline(2) = 0.0D0
 B2_spline(3) = 0.0D0
 B2_spline(4) = 0.0D0
 B2_spline(5) = 0.0D0
 B2_spline(6) = 0.0D0
 B2_spline(7) = 0.0D0

 DO i= N_k+1, 3*N_K
    TW_spline(i) = TW_spline(i-N_K)
    TW2_spline(i) = TW2_spline(i-N_K)
    B_spline(i) = B_spline(i-N_K)
    B2_spline(i) = B2_spline(i-N_K)
 END DO
 

 
 yp1 = 1.e+30
 ypn = 1.e+30
 
 CALL SPLINE(t_spline,TW_spline,3*N_K,yp1,ypn,y2it1)
 CALL SPLINE(t_spline,TW2_spline,3*N_K,yp1,ypn,y2it2)
 CALL SPLINE(t_spline,B_spline,3*N_K,yp1,ypn,y2ib1)
 CALL SPLINE(t_spline,B2_spline,3*N_K,yp1,ypn,y2ib2)

 DO i=1, N_steps
    CALL SPLINT(t_spline,TW_spline,y2it1,3*N_K,t(i),OUT_spline)
    CALL SPLINT(t_spline,TW2_spline,y2it2,3*N_K,t(i),OUT_spline2)
    CALL SPLINT(t_spline,B_spline,y2ib1,3*N_K,t(i),OUT_spline3)
    CALL SPLINT(t_spline,B2_spline,y2ib2,3*N_K,t(i),OUT_spline4)
    Y_spline(i) = OUT_spline
    Y2_spline(i) = OUT_spline2
    Y3_spline(i) = OUT_spline3
    Y4_spline(i) = OUT_spline4
    
 END DO


END SUBROUTINE MORPHING_SPLINE_KNOTS
