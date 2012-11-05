SUBROUTINE INF_MATRIXWaWi
USE global ! module for global variables
IMPLICIT NONE

! Procedure interfaces 
INTERFACE
  FUNCTION  BSL_SYM(P1,P2,PC) RESULT(OUTP)
    DOUBLE PRECISION, DIMENSION(3) :: P1, P2, PC
    DOUBLE PRECISION, DIMENSION(3) :: OUTP
  END
END INTERFACE

INTEGER :: K, L, m, n
INTEGER :: m1, n1

DOUBLE PRECISION, DIMENSION(3) :: P1, P2, PC

DOUBLE PRECISION, DIMENSION(3) :: Vel1, Vel2, Vel3, Vel4
DOUBLE PRECISION, DIMENSION(3) :: UVW, UVW_stream, UVW_stream_and_TE

IF(i_time .GT. 1) THEN

K = 0 

DO m=1, NR
   DO n=1, NC
      K = K + 1
      L = 0
      DO m1=1, i_time-1
         DO n1=1, NC
            L = L + 1
            ! control point
            PC(1) = XC(m,n,i_time)
            PC(2) = YC(m,n,i_time)
            PC(3) = ZC(m,n,i_time) 

            ! 1st segment
            P2(1) = X_Wake(m1+1,n1+1)
            P2(2) = Y_Wake(m1+1,n1+1)
            P2(3) = Z_Wake(m1+1,n1+1)

            P1(1) = X_Wake(m1+1,n1)
            P1(2) = Y_Wake(m1+1,n1)
            P1(3) = Z_Wake(m1+1,n1)
           
            Vel1(:) = BSL_SYM(P1,P2,PC)

            ! 2nd segement
            P2(1) = X_Wake(m1,n1+1)
            P2(2) = Y_Wake(m1,n1+1)
            P2(3) = Z_Wake(m1,n1+1)

            P1(1) = X_Wake(m1+1,n1+1)
            P1(2) = Y_Wake(m1+1,n1+1)
            P1(3) = Z_Wake(m1+1,n1+1)
           
            Vel2(:) = BSL_SYM(P1,P2,PC)

            ! 3rd segement
            P2(1) = X_Wake(m1,n1)
            P2(2) = Y_Wake(m1,n1)
            P2(3) = Z_Wake(m1,n1)

            P1(1) = X_Wake(m1,n1+1)
            P1(2) = Y_Wake(m1,n1+1)
            P1(3) = Z_Wake(m1,n1+1)
           
            Vel3(:) = BSL_SYM(P1,P2,PC)

            ! left segement
            P2(1) = X_Wake(m1+1,n1)
            P2(2) = Y_Wake(m1+1,n1)
            P2(3) = Z_Wake(m1+1,n1)

            P1(1) = X_Wake(m1,n1)
            P1(2) = Y_Wake(m1,n1)
            P1(3) = Z_Wake(m1,n1)
           
            Vel4(:) = BSL_SYM(P1,P2,PC)
 
            UVW(:) = Vel1(:) + Vel2(:) + Vel3(:) + Vel4(:)
            UVW_stream(:) = Vel2(:) + Vel4(:)
            UVW_stream_and_TE(:) = Vel2(:) + Vel3(:) + Vel4(:) 
            ! Wake on wing influence matrix (total contribution)
            C_mat_Wake(K,L) = UVW(1)*Outward_X(m,n,i_time) + UVW(2)*Outward_Y(m,n,i_time) + UVW(3)*Outward_Z(m,n,i_time)

            ! Wake on wing influence matrix (just induced downwash)
            B_mat_Wake(K,L) = UVW(1)*N_lift_X(m,n,i_time) + UVW(2)*N_lift_Y(m,n,i_time) + UVW(3)*N_lift_Z(m,n,i_time)
            
            END DO
    END DO
 END DO
END DO


END IF 

END SUBROUTINE INF_MATRIXWaWi
