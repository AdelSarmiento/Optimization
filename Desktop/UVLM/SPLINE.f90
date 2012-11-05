SUBROUTINE SPLINE(xi,yi,n,yp1,ypn,y2)
USE global ! module for global variables

INTEGER :: n
DOUBLE PRECISION :: yp1, ypn
DOUBLE PRECISION, DIMENSION(n) :: xi, yi, y2

INTEGER, PARAMETER :: NMAX = 500

INTEGER :: i,k
DOUBLE PRECISION :: p, qn, sig, un
DOUBLE PRECISION, DIMENSION(NMAX) :: u

IF(yp1 .GT. 0.99E30) THEN
  y2(1) = 0.0D0
  u(1) = 0.0D0
ELSE
  y2(1) = -0.5D0
  u(1) = (3.0D0/( xi(2)-xi(1) ))*( (yi(2) - yi(1))/(xi(2)-xi(1)) - yp1 )
END IF

DO i=2, n-1
   sig = (xi(i) - xi(i-1))/(xi(i+1) - xi(i-1))
   p = sig*y2(i-1) + 2.0D0
   y2(i) = (sig - 1.0D0)/p
   u(i) = ( 6.0D0*((yi(i+1)-yi(i))/(xi(i+1)-xi(i))-(yi(i)-yi(i-1))/(xi(i)-xi(i-1)))/(xi(i+1)-xi(i-1))-sig*u(i-1) )/p
END DO   

IF(ypn .GT. 0.99E30) THEN
 qn = 0.0D0
 un = 0.0D0
ELSE
 qn = 0.5D0
 un = ( 3.0D0/(xi(n)-xi(n-1)))*(ypn-(yi(n)-yi(n-1))/(xi(n)-xi(n-1)))
END IF

y2(n) = (un-qn*u(n-1))/(qn*y2(n-1) + 1.0D0)

DO k=n-1,1,-1
   y2(k) = y2(k)*y2(k+1) + u(k)
END DO

RETURN

END SUBROUTINE SPLINE
