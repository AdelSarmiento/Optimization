SUBROUTINE OUTPUT_DATA(GEOMVAL3, ELTL)
USE global ! module for global variables 
USE iga_global
DOUBLE PRECISION, DIMENSION(18) :: GEOMVAL3
DOUBLE PRECISION, DIMENSION(5) :: ELTL

 CALL GRID_GENERATION(GEOMVAL3)

 CALL KINEMATICS_PARAMETER_VALUES 

 CALL TIME_VARYING_KINEMATICS
 
 CALL MORPHING_WING_SHAPE 

 CALL INITIALIZE

 !!!!!!!!!!!!!!!!!!!!! TIME MARCHING !!!!!!!!!!!!!!!!!!!

 DO i_time=1, N_steps
    
    CALL VEL_KIN

    CALL VECTS

    CALL LD_VECTS

    CALL CONNECT_WAKE

    CALL INF_MATRIXWiWi

    CALL RHS
    
    CALL INF_MATRIXWaWi
    
    ! Solve linear system A*G = RHS
    Mat(1:NR*NC,1:NR*NC) = C_mat(1:NR*NC,1:NR*NC,i_time)
    IF(i_time == 1) THEN
        R(:) = LL(:,i_time)
        CALL DGESV(NEL, 1, Mat, NEL, IPIV, R, nel, INFO) 
        G(1:NEL) = R(1:NEL) 
        
        ! induced downwash
        VV = MATMUL(G(:),TRANSPOSE(B_mat(:,:,i_time)))
        CALL RESHAPE_MAT(VV,NC,NR,MVV)
        W_ind(:,:,i_time) = TRANSPOSE(MVV)
        
     ELSE  
        MVV2 = TRANSPOSE(Circ_Wake(:,:))
        CALL RESHAPE_MAT2(MVV2,NC,i_time-1,VV2) 
        ! right-hand side (including wake effect)
        R(:) = LL(:,i_time) - MATMUL(VV2,TRANSPOSE(C_mat_Wake))
        CALL DGESV(NEL, 1, Mat, NEL, IPIV, R, nel, INFO) 
        G(1:NEL) = R(1:NEL) 
        
        ! induced downwash
        VV = MATMUL(G(:),TRANSPOSE(B_mat(:,:,i_time)))
        
        CALL RESHAPE_MAT(VV,NC,NR,MVV)
        VV3 = MATMUL(VV2,TRANSPOSE(B_mat_Wake)) 
        CALL RESHAPE_MAT(VV3,NC,NR,MVV3)
        W_ind(:,:,i_time) = TRANSPOSE(MVV) + TRANSPOSE(MVV3) 
        
     END IF
     
     CALL RESHAPE_MAT(G(:),NC,NR,MVV)
     Circ(:,:,i_time) = TRANSPOSE(MVV)
     
     Circ_foo(:,:,i_time) = Circ(:,:,i_time)
     Circ_foo(2:NR,:,i_time) = Circ(2:NR,:,i_time) - Circ(1:NR-1,:,i_time)
     Circ_local(:,:,i_time) = Circ_foo(:,:,i_time)

     
     CALL LOADS

     CALL INF_WiWa 

     CALL INF_WaWa

     CALL CONVECT
     
     
 END DO

     CALL PERFORMANCE

      ELTL(1) = -EFF
      ELTL(2) = -CL_AVG
      ELTL(3) = -CT_AVG
	  ELTL(4) = Area
	  ELTL(5) = ALPHA_MAX
	  
	  CLOSE(1)
     
END SUBROUTINE OUTPUT_DATA
