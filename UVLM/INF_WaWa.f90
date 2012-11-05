SUBROUTINE INF_WaWa
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

DO m=wake_start, i_time-1
   DO n=1, NC+1
      K = K + 1
      L = 0
      DO m1=1, i_time-1
         DO n1=1, NC
            L = L + 1
            ! wake point
            PC(1) = X_Wake(m,n)
            PC(2) = Y_Wake(m,n)
            PC(3) = Z_Wake(m,n) 

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
            
            ! wake on wake influence matrix
            C2_rollup(K,L,:) = UVW(:)
            
            END DO
    END DO
 END DO
END DO


END IF 

END SUBROUTINE INF_WaWa
