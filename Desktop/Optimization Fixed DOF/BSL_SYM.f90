FUNCTION  BSL_SYM(P1,P2,PC) RESULT(OUTP)   ! Function to implement the Biot-Savart Law for symmetric flow
USE global
IMPLICIT NONE

! cross product of 2 vectors
INTERFACE
  FUNCTION  CROSS_PDT(INP1,INP2) RESULT(OUTP)
    DOUBLE PRECISION, DIMENSION(3) :: INP1, INP2
    DOUBLE PRECISION, DIMENSION(3) :: OUTP
  END
END INTERFACE

! variables & parameters 

! Dummy variables.
DOUBLE PRECISION, DIMENSION(3) :: P1, P2, PC, PCi
DOUBLE PRECISION, DIMENSION(3) :: OUTP

DOUBLE PRECISION, DIMENSION(3) :: R0, R1, R2

DOUBLE PRECISION, DIMENSION(3) :: AA, CC

DOUBLE PRECISION, DIMENSION(3) :: UVW1, UVW2

! Local variables.
DOUBLE PRECISION :: BB

 R0(:) = P2(:) - P1(:)
 R1(:) = PC(:) - P1(:)
 R2(:) = PC(:) - P2(:)
 
 AA(:) = CROSS_PDT(R1(:),R2(:))
 
 BB = AA(1)**2+AA(2)**2+AA(3)**2
 IF(BB .GT. Cut_off) THEN 
   CC(:) = R1(:)/(SQRT(R1(1)**2+R1(2)**2+R1(3)**2))-R2(:)/(SQRT(R2(1)**2+R2(2)**2+R2(3)**2))
   UVW1(:) = (R0(1)*CC(1)+R0(2)*CC(2)+R0(3)*CC(3))/BB * AA(:)/4.0D0/PI
 ELSE
   UVW1(:) = 0.0D0
 END IF
 
! the contribution from the image
 
 PCi(1) =   PC(1)
 PCi(2) = - PC(2)
 PCi(3) =   PC(3)

 R0(:) = P2(:) - P1(:)
 R1(:) = PCi(:) - P1(:)
 R2(:) = PCi(:) - P2(:)

 AA(:) = CROSS_PDT(R1(:),R2(:))
 BB = AA(1)**2+AA(2)**2+AA(3)**2
 IF(BB .GT. Cut_off) THEN 
   CC(:) = R1(:)/(SQRT(R1(1)**2+R1(2)**2+R1(3)**2))-R2(:)/(SQRT(R2(1)**2+R2(2)**2+R2(3)**2))
   UVW2(:) = (R0(1)*CC(1)+R0(2)*CC(2)+R0(3)*CC(3))/BB * AA(:)/4.0D0/PI
 ELSE
   UVW2(:) = 0.0D0
 END IF

 UVW2(2) = - UVW2(2)

 
! The velocity components

OUTP(:) = UVW1(:) + UVW2(:)


END FUNCTION  BSL_SYM
