SUBROUTINE LD_VECTS
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
 
DOUBLE PRECISION, DIMENSION(3) :: VN, VT, VQ
DOUBLE PRECISION, DIMENSION(3) :: V_drag, V_lift
DOUBLE PRECISION :: V_H, V_V

DOUBLE PRECISION, DIMENSION(3,3) :: T_foo, T_alpha

! Compute angle of attack for each panel
DO i=1, NR
   DO j=1, NC
      IF(i_time == 1) THEN
      ALPHA(i,j,i_time) = 0.0D0
      ELSE
      V_H = U_KIN(i,j,i_time)*Tau_X(i,j,i_time) + V_KIN(i,j,i_time)& 
	  *Tau_Y(i,j,i_time) + W_KIN(i,j,i_time)*Tau_Z(i,j,i_time) ! horizontal component
      V_V = U_KIN(i,j,i_time)*Outward_X(i,j,i_time) + V_KIN(i,j,i_time)& 
	  *Outward_Y(i,j,i_time) + W_KIN(i,j,i_time)*Outward_Z(i,j,i_time) ! vertical component
      ALPHA(i,j,i_time) = ATAN(V_V/V_H)
      END IF
      VN(1) = Outward_X(i,j,i_time)
      VN(2) = Outward_Y(i,j,i_time)
      VN(3) = Outward_Z(i,j,i_time)

      VT(1) = Tau_X(i,j,i_time)
      VT(2) = Tau_Y(i,j,i_time)
      VT(3) = Tau_Z(i,j,i_time)

      VQ(:) = CROSS_PDT(VT(:),VN(:)) 
 
      T_foo(:,1) = VT(:)
      T_foo(:,2) = VN(:)
      T_foo(:,3) = VQ(:)

      T_foo = TRANSPOSE(T_foo)

      T_alpha(:,:) =   0.0D0

      T_alpha(1,1) =   COS(ALPHA(i,j,i_time))
      T_alpha(2,1) =   SIN(ALPHA(i,j,i_time))
      T_alpha(1,2) = - SIN(ALPHA(i,j,i_time))
      T_alpha(2,2) =   COS(ALPHA(i,j,i_time))
      T_alpha(3,3) =   1.0D0

      V_drag(:)  =   MATMUL( VT(:), MATMUL( TRANSPOSE(T_foo), & 
	  MATMUL(TRANSPOSE(T_alpha), T_foo )  ) )
      V_lift(:)  =   MATMUL( VN(:), MATMUL( TRANSPOSE(T_foo), & 
	  MATMUL(TRANSPOSE(T_alpha), T_foo )  ) )

      N_drag_X(i,j,i_time) = V_drag(1)
      N_drag_Y(i,j,i_time) = V_drag(2)
      N_drag_Z(i,j,i_time) = V_drag(3)   
      
      N_lift_X(i,j,i_time) = V_lift(1)
      N_lift_Y(i,j,i_time) = V_lift(2)
      N_lift_Z(i,j,i_time) = V_lift(3)  
     
   END DO
END DO

END SUBROUTINE LD_VECTS
