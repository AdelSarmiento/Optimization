SUBROUTINE INF_MATRIXWiWi
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

K = 0 

DO m=1, NR
   DO n=1, NC
      K = K + 1
      L = 0
      DO m1=1, NR
         DO n1=1, NC
            L = L + 1
            ! control point
            PC(1) = XC(m,n,i_time)
            PC(2) = YC(m,n,i_time)
            PC(3) = ZC(m,n,i_time) 

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
            UVW_stream(:) = Vel2(:) + Vel4(:)
            UVW_stream_and_TE(:) = Vel2(:) + Vel3(:) + Vel4(:) 
            
            C_mat(K,L,i_time) = UVW(1)*Outward_X(m,n,i_time) & 
			+ UVW(2)*Outward_Y(m,n,i_time) + UVW(3)*Outward_Z(m,n,i_time)

            IF(m1 .LT. NR) THEN
               B_mat(K,L,i_time) = UVW_stream(1)*N_lift_X(m,n,i_time) & 
			   + UVW_stream(2)*N_lift_Y(m,n,i_time) + UVW_stream(3)*N_lift_Z(m,n,i_time)
            ELSE 
               B_mat(K,L,i_time) = UVW_stream_and_TE(1)& 
			   *N_lift_X(m,n,i_time) + UVW_stream_and_TE(2)*N_lift_Y(m,n,i_time) & 
			   + UVW_stream_and_TE(3)*N_lift_Z(m,n,i_time)
            END IF
            END DO
    END DO
 END DO
END DO 

END SUBROUTINE INF_MATRIXWiWi
