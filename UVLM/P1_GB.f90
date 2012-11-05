!********+*********+*********+*********+*********+*********+*********+
!
      SUBROUTINE INITI(M,N,GEPS,XVAL,XMIN,XMAX,FMAX,A,C)
      
	  USE global
	  USE iga_global
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XVAL(1),XMIN(1),XMAX(1),FMAX(4),A(4),C(4)

      M = 4 ! number of constraints
      N = 12 ! number of design variables 
	  
      GEPS=0.0000001
      !OPEN(7777, FILE = 'DATA_OPT3.dat')
      !READ(7777,*) (INVAL(j), j=1,12)
      DO J=1,N
      XMIN(J)= -5.0D0
      XMAX(J)= 5.0D0
      XVAL(J)= INVAL(J)
      END DO
	  
	  XMIN(N-1)= -0.5D0
      XMAX(N-1)= 0.5D0
      !XVAL(N-1)= 0.0D0
	  
	  XMIN(N)= -0.5D0
      XMAX(N)= 0.5D0
      !XVAL(N)= 0.0D0

      
       
      FMAX(1)=-4.484D0
      A(1)=0.0D0
      C(1)=1000.0D0

      FMAX(2)=-0.1852D0
      A(2)=0.0D0
      C(2)=1000.0D0
	  
	  FMAX(3)=3.03D0
      A(3)=0.0D0
      C(3)=100000.0D0
	  
	  FMAX(4)=0.25D0
      A(4)=0.0D0
      C(4)=100000.0D0
!
      RETURN
      END
!
!********+*********+*********+*********+*********+*********+*********+
!
      SUBROUTINE FUNC1(M,N,X,F0,F)
      USE global
!  !-----------------------------------------------------!
!  !  The author of this subroutine is Krister Svanberg  !
!  !-- --------------------------------------------------!
!
!     Calculation of function values but no gradients for the
!     three bar truss example.
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(1),DF0DX(1),F(4),DFDX(1)
!
      CALL OUTPUT_DATA(X,ELT)
      
      F0 = ELT(1)
      F(1) = ELT(2)
      F(2) = ELT(3)
	  F(3) = ELT(4)
	  F(4) = ELT(5)
                  
      RETURN
      END
!
!********+*********+*********+*********+*********+*********+*********+
!
      SUBROUTINE FUNC2(M,N,X,F0,DF0DX,F,DFDX)
      USE global
!  !-----------------------------------------------------!
!  !  The author of this subroutine is Krister Svanberg  !
!  !-- --------------------------------------------------!
!
!     Calculation of function values and gradients for the
!     three bar truss example.
!     Since this problem is so small, both the function values
!     and the gradients are calculated analytically.
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	  DIMENSION X(1),DF0DX(1),F(4),DFDX(1)
      DIMENSION FPER(4)

      CALL OUTPUT_DATA(X,ELT)
!
      F0 = ELT(1)
      F(1) = ELT(2)
      F(2) = ELT(3)
	  F(3) = ELT(4)
	  F(4) = ELT(5)
      
      DO J=1,N
      X(J) = X(J) + 0.000001D0
      CALL OUTPUT_DATA(X,ELT)
      F0PER = ELT(1)
      FPER(1) = ELT(2)
      FPER(2) = ELT(3)
	  FPER(3) = ELT(4)
	  FPER(4) = ELT(5)
      DF0DX(J) = (F0PER-F0)/0.000001D0
      DO I=1, M
         IJ=M*(J-1)+I
         DFDX(IJ) = (FPER(I)-F(I))/0.000001D0
      END DO
      X(J) = X(J) - 0.000001D0
      END DO

      
      RETURN
      END
!
!********+*********+*********+*********+*********+*********+*********+
!
      SUBROUTINE OUTXOF(ITER,INNER,M,N,X,F0VAL,FVAL,FMAX)
      USE global
!     OUTXOF writes the current solution
!     for the three bar truss problem.
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(28),FVAL(2),FMAX(1)
	  
!
      WRITE(*,*)
      WRITE(*,90) ITER
 90   FORMAT(' Iteration number:',I4)
      IF(ITER.GE.1) WRITE(*,91) INNER
 91   FORMAT(' Inner iterations =',I3)
      WRITE(*,93) -F0VAL
 93   FORMAT('        Efficiency:',F12.5)
      WRITE(*,94)  -FVAL(1)
 94   FORMAT('         Lift:',F12.5)
      WRITE(*,95)  -FVAL(2)
 95   FORMAT('         Thrust:',F12.5)
  WRITE(*,96)  Area
 96   FORMAT('         Area:',F12.5)
 WRITE(*,97)  CP_IN_AVG
 97   FORMAT('         Power:',F12.5)
 WRITE(*,98)  ALPHA_MAX*180.0D0/PI
 98   FORMAT('         Angle:',F12.5)
!
      WRITE(1111,1000) REAL(ITER), -F0VAL, -FVAL(1),&
	                   -FVAL(2), Area, CP_IN_AVG, X(1:12)

1000     FORMAT(1x, 18F12.5)
      RETURN
      END
!
!********+*********+*********+*********+*********+*********+*********+
