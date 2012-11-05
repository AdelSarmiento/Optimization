SUBROUTINE TIME_VARYING_KINEMATICS
USE global ! module for global variables


Ones(:) = 1.0D0 

Path_X(:) = -V_ref * t(:)
Path_Y(:) = 0.0D0
Path_Z(:) = PLamp * SIN(Omega*t(:)+PLphase*Ones(:))

Roll(:) = FLamp * SIN(Omega*t(:)+FLphase*Ones(:))
Pitch(:) = Pitchamp

END SUBROUTINE TIME_VARYING_KINEMATICS
