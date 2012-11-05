SUBROUTINE SPLINT(xa,ya,y2a,n,xv,yv)
USE global ! module for global variables

INTEGER :: n
DOUBLE PRECISION :: xv, yv
DOUBLE PRECISION, DIMENSION(n) :: xa, ya, y2a

INTEGER :: k, khi, klo
DOUBLE PRECISION :: a, bb2, h

klo = 1
khi = n
1 IF((khi-klo) .GT. 1) THEN
    k = (khi+klo)/2
    IF(xa(k) .GT. xv) THEN
      khi = k
    ELSE
      klo = k
    END IF

GOTO 1

END IF

h = xa(khi)-xa(klo)
!IF(h .EQ. 0.0D0) pause 'bad xa input in splint

a = (xa(khi) - xv)/h
bb2 = (xv - xa(klo))/h

yv = a*ya(klo) + bb2*ya(khi) + ( (a**3-a)*y2a(klo)+(bb2**3-bb2)*y2a(khi) )*(h**2)/6.0D0

RETURN 

END SUBROUTINE SPLINT
