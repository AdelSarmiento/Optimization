SUBROUTINE KINEMATICS_PARAMETER_VALUES
USE global ! module for global variables

! flapping
FLamp = 45.0D0 
FLphase = 0.0D0 

! plunging
PLamp = 0.0D0 
PLphase = 0.0D0 

! pitch angle
Pitchamp = 5.0D0 
 

! Conversion to radians
FLamp = FLamp * (PI/180.0D0)
FLphase = FLphase * (PI/180.0D0)

PLphase = PLphase * (PI/180.0D0)
Pitchamp = Pitchamp * (PI/180.0D0)


END SUBROUTINE KINEMATICS_PARAMETER_VALUES
