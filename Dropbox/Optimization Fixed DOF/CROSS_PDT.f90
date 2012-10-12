FUNCTION  CROSS_PDT(V1,V2) RESULT(OUTP)   ! Function to compute the cross product of 2 vectors
IMPLICIT NONE

! variables & parameters 

! Dummy variables.
DOUBLE PRECISION, DIMENSION(3) :: V1, V2, OUTP

OUTP(1) = V1(2)*V2(3) -V1(3)*V2(2)
OUTP(2) = -(V1(1)*V2(3) -V1(3)*V2(1))
OUTP(3) = V1(1)*V2(2) -V1(2)*V2(1)

END FUNCTION  CROSS_PDT
