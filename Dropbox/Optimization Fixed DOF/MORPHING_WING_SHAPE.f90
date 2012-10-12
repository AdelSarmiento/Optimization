SUBROUTINE MORPHING_WING_SHAPE
USE global ! module for global variables

INTEGER :: i, j, k
DOUBLE PRECISION, DIMENSION(3) :: XYZ, FOO

! Compute morphed wing shape in body-attached ref frame
DO i=1, N_steps
   DO j=1, NR+1
      DO k=1, NC+1
         twist_1 = twist1(i) * (Y0(j,k)/(chordL/2.0D0)) 
         twist_2 = twist2(i) * ( -192.0D0*(Y0(j,k)/chordL)**4 & 
		 + 224.0D0*(Y0(j,k)/chordL)**3 -60.0D0*(Y0(j,k)/chordL)**2)
         bend = bend1(i) * (Y0(j,k)/(chordL/2.0D0))**2 + bend2(i) & 
		 * ( -192.0D0*(Y0(j,k)/chordL)**4 + 224.0D0*(Y0(j,k)/chordL)**3 & 
		 -60.0D0*(Y0(j,k)/chordL)**2)
         XX(j,k,i) = X0(j,k) - X0(j,k) * (1.0D0 - COS( twist_1  ) ) & 
		 - X0(j,k) * (1.0D0 - COS( twist_2  ) )
         YY(j,k,i) = Y0(j,k) + fore(i)* (Y0(j,k)/(chordL/2.0D0))**2
         ZZ(j,k,i) = Z0(j,k) - X0(j,k) * SIN( twist_1 ) - X0(j,k) & 
		 * SIN( twist_2 ) + bend
         !WRITE(*,*) twist
      END DO
   END DO
   !WRITE(*,*) X(6,9,i), Y(6,9,i), Z(6,9,i) 
END DO
 
! Compute ring coordinates in body-attached ref frame
DO i=1, N_steps
   DO k=1, NC+1   
      DO j=1, NR
         XRb(j,k,i) = (XX(j+1,k,i)-XX(j,k,i))/4.0D0 + XX(j,k,i)
         YRb(j,k,i) =  YY(j,k,i)
         ZRb(j,k,i) = (ZZ(j+1,k,i)-ZZ(j,k,i))/4.0D0 + ZZ(j,k,i)
      END DO
         XRb(NR+1,k,i) = (XX(NR+1,k,i)-XX(NR,k,i))/4.0D0 + XX(NR+1,k,i)
         YRb(NR+1,k,i) =  YY(NR,k,i)
         ZRb(NR+1,k,i) =  ZZ(NR+1,k,i)
   END DO
END DO 

! Compute control points in body-attached ref frame
DO i=1, N_steps
   DO k=1, NC
      DO j=1, NR
         XCb(j,k,i) = (XRb(j,k,i)+XRb(j+1,k,i)+XRb(j,k+1,i)& 
		 +XRb(j+1,k+1,i))/4.0D0
         YCb(j,k,i) = (YRb(j,k,i)+YRb(j+1,k,i)+YRb(j,k+1,i)& 
		 +YRb(j+1,k+1,i))/4.0D0
         ZCb(j,k,i) = (ZRb(j,k,i)+ZRb(j+1,k,i)+ZRb(j,k+1,i)& 
		 +ZRb(j+1,k+1,i))/4.0D0
      END DO
   END DO
END DO


DO i=1, N_steps

   ! transformation matrix between inertial and body-attached frames
   
   ! roll/flap (rotation about x-axis)
   T_roll(:,:) = 0.0D0

   T_roll(1,1) =   1.0D0
   T_roll(2,2) =   COS(Roll(i))
   T_roll(3,2) = - SIN(Roll(i))
   T_roll(2,3) =   SIN(Roll(i))
   T_roll(3,3) =   COS(Roll(i))

   ! pitch (rotation about y-axis)
   T_pitch(:,:) = 0.0D0

   T_pitch(1,1) =   COS(Pitch(i))
   T_pitch(3,1) =   SIN(Pitch(i))
   T_pitch(2,2) =   1.0D0
   T_pitch(1,3) = - SIN(Pitch(i))
   T_pitch(3,3) =   COS(Pitch(i))

   ! transformation matrix
   T_matrix = MATMUL( T_roll, T_pitch )

  DO j=1, NR+1
     DO k=1, NC+1
        XYZ(1) = XRb(j,k,i)
        XYZ(2) = YRb(j,k,i)
        XYZ(3) = ZRb(j,k,i)
        FOO(:) = MATMUL( XYZ(:), T_matrix(:,:) )  
        IF (XYZ(2) .EQ. 0.0D0) THEN
           FOO(2) = 0.0D0   ! prevents root of cambered wing from crossing to negative Y values
        END IF
       ! ring element nodes in inertial ref frame 
       XR(j,k,i) = FOO(1) + Path_X(i)
       YR(j,k,i) = FOO(2) + Path_Y(i)
       ZR(j,k,i) = FOO(3) + Path_Z(i)
        
       ! control points in inertial ref frame 
       IF( (j .LT. NR+1) .AND. (k .LT. NC+1) ) THEN
           XYZ(1) = XCb(j,k,i)
           XYZ(2) = YCb(j,k,i)
           XYZ(3) = ZCb(j,k,i)
           FOO(:) = MATMUL( XYZ(:), T_matrix(:,:) ) 
           XC(j,k,i) = FOO(1) + Path_X(i)
           YC(j,k,i) = FOO(2) + Path_Y(i)
           ZC(j,k,i) = FOO(3) + Path_Z(i)
       END IF
     END DO ! wing columns
  END DO ! wing rows

END DO ! time

END SUBROUTINE MORPHING_WING_SHAPE
