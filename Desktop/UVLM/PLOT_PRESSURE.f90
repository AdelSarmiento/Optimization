SUBROUTINE PLOT_PRESSURE
USE global
IMPLICIT NONE

! variables & parameters 

! Local variables.
 INTEGER :: i, j, nrec, nS
 CHARACTER*50 :: fname

 nS = 10000 + i_time


    WRITE(fname,341)nS 
341 FORMAT('SnapPr_',i5,'.dat')


OPEN(unit=97,file=fname,form='FORMATTED')

WRITE(97,*)'VARIABLE "X","Y","Z","P"'

WRITE(97,*)'ZONE T="ZONE1", I=',nr,' ,J=',2*nc

DO j=nc, 1, -1
   DO i=1, nr
      WRITE(97,7) XC(i,j,i_time), YC(i,j,i_time), ZC(i,j,i_time), DCP(i,j,i_time)/(rho*V_ref**2)
   END DO
END DO

DO j=1, nc
   DO i=1, nr
      WRITE(97,7) XC(i,j,i_time), - YC(i,j,i_time), ZC(i,j,i_time), DCP(i,j,i_time)/(rho*V_ref**2)
   END DO
END DO

IF (i_time .GT. 10) THEN 
   DO j=1, nc
      DO i=1, nr
         WRITE(123,77) DCP(i,j,i_time)/(rho*V_ref**2)
      END DO
   END DO
ENDIF


7     FORMAT(1x, 4(1pe12.5,2x))
      CLOSE(97)
77    FORMAT(1x, 41(1pe12.5,2x))

END SUBROUTINE  PLOT_PRESSURE
