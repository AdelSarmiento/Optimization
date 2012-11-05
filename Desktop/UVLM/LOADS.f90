SUBROUTINE LOADS
USE global ! module for global variables

INTEGER :: i, j
DOUBLE PRECISION :: DL, DD
DOUBLE PRECISION, DIMENSION(NR,NC) :: CHORD, SPAN
DOUBLE PRECISION, DIMENSION(NR,NC,N_steps) :: F_X, F_Y, F_Z

Sigma1(1,:,i_time) = 0.5D0 * Circ(1,:,i_time)
Sigma1(2:NR,:,i_time) = 0.5D0 * (Circ(1:NR-1,:,i_time)+Circ(2:NR,:,i_time))


IF(i_time .EQ. 1) THEN 
   DFDT(:,:,i_time) = SIGMA1(:,:,i_time)/t_step
ELSE
   DFDT(:,:,i_time) = (SIGMA1(:,:,i_time)-SIGMA1(:,:,i_time-1))/t_step
END IF
 


 POWER(i_time) = 0.0D0
 CD(i_time) = 0.0D0
 CL(i_time) = 0.0D0

DO i=1, NR
   DO j=1, NC

      CHORD(i,j) = SQRT((0.5D0*(XX(i+1,j,i_time)+XX(i+1,j+1,i_time))- & 
	  0.5D0*(XX(i,j,i_time)+XX(i,j+1,i_time)))**2 &
                        + (0.5D0*(YY(i+1,j,i_time)+YY(i+1,j+1,i_time))- & 
						0.5D0*(YY(i,j,i_time)+YY(i,j+1,i_time)))**2 &
                        + (0.5D0*(ZZ(i+1,j,i_time)+ZZ(i+1,j+1,i_time))- & 
						0.5D0*(ZZ(i,j,i_time)+ZZ(i,j+1,i_time)))**2    )

      SPAN(i,j) = SQRT( (0.5D0*(XX(i,j+1,i_time)+XX(i+1,j+1,i_time))- & 
	  0.5D0*(XX(i,j,i_time)+XX(i+1,j,i_time)))**2 &
                        + (0.5D0*(YY(i,j+1,i_time)+YY(i+1,j+1,i_time))- & 
						0.5D0*(YY(i,j,i_time)+YY(i+1,j,i_time)))**2 &
                        + (0.5D0*(ZZ(i,j+1,i_time)+ZZ(i+1,j+1,i_time))- & 
						0.5D0*(ZZ(i,j,i_time)+ZZ(i+1,j,i_time)))**2     )

      ! lift on the panel
      DL = rho * SPAN(i,j) * ( SQRT( U_KIN(i,j,i_time)**2 + & 
	  V_KIN(i,j,i_time)**2 + W_KIN(i,j,i_time)**2 ) & 
	  * Circ_local(i,j,i_time) & 
                  + CHORD(i,j)*DFDT(i,j,i_time) )*COS(ALPHA(i,j,i_time))                                       
      ! drag on the panel
      DD = rho * SPAN(i,j) & 
	  * ( - W_ind(i,j,i_time)*Circ_local(i,j,i_time) & 
                 + CHORD(i,j)*DFDT(i,j,i_time)*SIN(ALPHA(i,j,i_time)) )
      !WRITE(*,*) DL
      !WRITE(*,*)
      !WRITE(*,*) DD
      !WRITE(*,*)
      
      ! total force
      F_X(i,j,i_time) = DL*N_lift_X(i,j,i_time) + DD*N_drag_X(i,j,i_time) 
      F_Y(i,j,i_time) = DL*N_lift_Y(i,j,i_time) + DD*N_drag_Y(i,j,i_time) 
      F_Z(i,j,i_time) = DL*N_lift_Z(i,j,i_time) + DD*N_drag_Z(i,j,i_time)
      
      ! pressure 
      DCP(i,j,i_time) = (F_X(i,j,i_time)*Outward_X(i,j,i_time) & 
	  + F_Y(i,j,i_time)*Outward_Y(i,j,i_time) & 
	  + F_Z(i,j,i_time)*Outward_Z(i,j,i_time))/(CHORD(i,j)*SPAN(i,j))
      
      ! power
      POWER(i_time) = POWER(i_time) + 2.0D0*DCP(i,j,i_time)& 
	  *(  U_KIN(i,j,i_time)*Outward_X(i,j,i_time) & 
	  + V_KIN(i,j,i_time)*Outward_Y(i,j,i_time) &
    + W_KIN(i,j,i_time)*Outward_Z(i,j,i_time)   ) * CHORD(i,j)*SPAN(i,j)
      ! drag coefficient
      CD(i_time) = CD(i_time) + F_X(i,j,i_time)*2.0D0& 
	  /(rho*0.5D0*V_ref**2)
      ! lift coefficient
      CL(i_time) = CL(i_time) + F_Z(i,j,i_time)*2.0D0&
	  /(spanL*rho*0.5D0*V_ref**2)

   END DO
END DO
 !WRITE(*,*) Area  
  WRITE(1,1000) t(i_time)/(2.0D0*PI/omega)-0.25D0, CD(i_time), & 
                CL(i_time), POWER(i_time)
1000     FORMAT(1x, 4(1pe12.5,2x))  
 
END SUBROUTINE LOADS
