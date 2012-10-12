SUBROUTINE INF_WiWa
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
  IF(i_time .LT. N_W) THEN
     wake_start = 1
  ELSE
     wake_start = i_time - N_W
  END IF

  K=0
  DO m = wake_start, i_time-1
     DO n = 1, NC+1
        K = K + 1
        L = 0
        DO m1=1, NR
           DO n1=1, NC
              L = L + 1
              ! wake point
              PC(1) = X_Wake(m,n)
              PC(2) = Y_Wake(m,n)
              PC(3) = Z_Wake(m,n) 

              ! upper segment
              P2(1) = XR(m1,n1+1,i_time)
              P2(2) = YR(m1,n1+1,i_time)
              P2(3) = ZR(m1,n1+1,i_time)

              P1(1) = XR(m1,n1,i_time)
              P1(2) = YR(m1,n1,i_time)
              P1(3) = ZR(m1,n1,i_time)
           
              Vel1(:) = BSL_SYM(P1,P2,PC)

              ! right segement
              P2(1) = XR(m1+1,n1+1,i_time)
              P2(2) = YR(m1+1,n1+1,i_time)
              P2(3) = ZR(m1+1,n1+1,i_time)

              P1(1) = XR(m1,n1+1,i_time)
              P1(2) = YR(m1,n1+1,i_time)
              P1(3) = ZR(m1,n1+1,i_time)
           
              Vel2(:) = BSL_SYM(P1,P2,PC)

              ! lower segement
              P2(1) = XR(m1+1,n1,i_time)
              P2(2) = YR(m1+1,n1,i_time)
              P2(3) = ZR(m1+1,n1,i_time)

              P1(1) = XR(m1+1,n1+1,i_time)
              P1(2) = YR(m1+1,n1+1,i_time)
              P1(3) = ZR(m1+1,n1+1,i_time)
           
              Vel3(:) = BSL_SYM(P1,P2,PC)

              ! left segement
              P2(1) = XR(m1,n1,i_time)
              P2(2) = YR(m1,n1,i_time)
              P2(3) = ZR(m1,n1,i_time)

              P1(1) = XR(m1+1,n1,i_time)
              P1(2) = YR(m1+1,n1,i_time)
              P1(3) = ZR(m1+1,n1,i_time)
           
              Vel4(:) = BSL_SYM(P1,P2,PC)
 
              UVW(:) = Vel1(:) + Vel2(:) + Vel3(:) + Vel4(:)
              
              ! Wing on wake influence matrix

              C1_rollup(K,L,:) = UVW(:)
              
          END DO
       END DO
    END DO
  END DO

END IF
END SUBROUTINE INF_WiWa
