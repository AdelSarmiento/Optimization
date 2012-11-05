SUBROUTINE VECTS
USE global ! module for global variables
IMPLICIT NONE

! cross product of 2 vectors
INTERFACE
  FUNCTION  CROSS_PDT(INP1,INP2) RESULT(OUTP)
    DOUBLE PRECISION, DIMENSION(3) :: INP1, INP2
    DOUBLE PRECISION, DIMENSION(3) :: OUTP
  END
END INTERFACE

INTEGER :: i, j
DOUBLE PRECISION, DIMENSION(3) :: V1, V2
DOUBLE PRECISION, DIMENSION(3) :: VN, VT
DOUBLE PRECISION :: AMPVN, AMPVT

Area = 0.0D0

DO i=1, NR
   DO j=1, NC
      ! 1st vector
      V1(1) = XR(i,j,i_time)-XR(i+1,j+1,i_time)
      V1(2) = YR(i,j,i_time)-YR(i+1,j+1,i_time) 
      V1(3) = ZR(i,j,i_time)-ZR(i+1,j+1,i_time)
      ! 2nd vector
      V2(1) = XR(i+1,j,i_time)-XR(i,j+1,i_time)
      V2(2) = YR(i+1,j,i_time)-YR(i,j+1,i_time) 
      V2(3) = ZR(i+1,j,i_time)-ZR(i,j+1,i_time)

      ! normal vector (unit vector)
      VN(:) = CROSS_PDT(V1(:),V2(:))
      AMPVN = SQRT(VN(1)**2+VN(2)**2+VN(3)**2)
	  Area = Area + AMPVN/2.0D0
      AreaElt(i,j) = AMPVN/2.0D0	  
      VN(:) = VN(:)/AMPVN
 
      Outward_X(i,j,i_time) = VN(1)
      Outward_Y(i,j,i_time) = VN(2)
      Outward_Z(i,j,i_time) = VN(3)

      ! tangent vector 
      VT(1) = 0.5D0*(XR(i+1,j,i_time)+XR(i+1,j+1,i_time)) & 
	                         - 0.5D0*(XR(i,j,i_time)+XR(i,j+1,i_time))
      VT(2) = 0.5D0*(YR(i+1,j,i_time)+YR(i+1,j+1,i_time)) & 
	                         - 0.5D0*(YR(i,j,i_time)+YR(i,j+1,i_time))
      VT(3) = 0.5D0*(ZR(i+1,j,i_time)+ZR(i+1,j+1,i_time)) & 
	                         - 0.5D0*(ZR(i,j,i_time)+ZR(i,j+1,i_time))
      AMPVT = SQRT(VT(1)**2+VT(2)**2+VT(3)**2)
      VT(:) = VT(:)/AMPVT

      Tau_X(i,j,i_time) = VT(1)
      Tau_Y(i,j,i_time) = VT(2)
      Tau_Z(i,j,i_time) = VT(3)
   END DO
END DO

 CALL RESHAPE_MAT2(AreaElt,NR,NC,Varea)

 !WRITE(*,*) MINVAL(Varea(:))

END SUBROUTINE VECTS
