SUBROUTINE CONNECT_WAKE
USE global ! module for global variables
IMPLICIT NONE

! Add to wake geometry
X_Wake(i_time,:) = XR(NR+1,:,i_time)
Y_Wake(i_time,:) = YR(NR+1,:,i_time)
Z_Wake(i_time,:) = ZR(NR+1,:,i_time)

! Set circulation of newest wake rings to be the old TE circulation

IF(i_time .GT. 1) THEN
  Circ_Wake(i_time-1,:) = Circ(NR,:,i_time-1)
END IF

END SUBROUTINE CONNECT_WAKE
