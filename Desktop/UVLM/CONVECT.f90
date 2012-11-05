SUBROUTINE CONVECT
USE global ! module for global variables

INTEGER :: i, j
IF(i_time .GT. 1) THEN

 MVV4 = TRANSPOSE(Circ(:,:,i_time))
 CALL RESHAPE_MAT2(MVV4,NC,NR,VV4)
  
 MVV5 = TRANSPOSE(Circ_Wake(:,:))
 CALL RESHAPE_MAT2(MVV5,NC,i_time-1,VV5)
  
 DEL_X = ( MATMUL(VV4, TRANSPOSE(C1_rollup(:,:,1)) ) + MATMUL(VV5, TRANSPOSE(C2_rollup(:,:,1)) ) ) * t_step
 DEL_Y = ( MATMUL(VV4, TRANSPOSE(C1_rollup(:,:,2)) ) + MATMUL(VV5, TRANSPOSE(C2_rollup(:,:,2)) ) ) * t_step
 DEL_Z = ( MATMUL(VV4, TRANSPOSE(C1_rollup(:,:,3)) ) + MATMUL(VV5, TRANSPOSE(C2_rollup(:,:,3)) ) ) * t_step
 
 CALL RESHAPE_MAT(DEL_X,NC+1,(N_steps-1),MDEL_X)
 CALL RESHAPE_MAT(DEL_Y,NC+1,(N_steps-1),MDEL_Y)
 CALL RESHAPE_MAT(DEL_Z,NC+1,(N_steps-1),MDEL_Z)

 X_Wake(wake_start:i_time-1,:) = X_Wake(wake_start:i_time-1,:) + TRANSPOSE(MDEL_X(:,:))
 Y_Wake(wake_start:i_time-1,:) = Y_Wake(wake_start:i_time-1,:) + TRANSPOSE(MDEL_Y(:,:))
 Z_Wake(wake_start:i_time-1,:) = Z_Wake(wake_start:i_time-1,:) + TRANSPOSE(MDEL_Z(:,:))
 
 !WRITE(*,*) X_Wake(:,:)
 !WRITE(*,*) 
 !WRITE(*,*) Y_Wake(:,:)
 !WRITE(*,*) 
 !WRITE(*,*) Z_Wake(:,:)
 
END IF
END SUBROUTINE CONVECT
