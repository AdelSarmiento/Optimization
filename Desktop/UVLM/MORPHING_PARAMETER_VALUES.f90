SUBROUTINE MORPHING_PARAMETER_VALUES
USE global ! module for global variables

! 1st twist
TW1amp = 0.0D0!8.94316D0 
TW1phase = 90.0D0!-0.49469D0 

! 2nd twist
TW2amp = 0.0D0 
TW2phase = -44.4D0 

! 1st bending
TB1amp = 0.0D0!0.69383D0 
TB1phase = -40.0D0!-89.011D0 

! 2nd bending
TB2amp = 0.0D0 
TB2phase = -156.0D0

! Conversion to radians
TW1amp = TW1amp * (PI/180.0D0)
TW2amp = TW2amp * (PI/180.0D0)
TW1phase = TW1phase * (PI/180.0D0)
TW2phase = TW2phase * (PI/180.0D0)
TB1phase = TB1phase * (PI/180.0D0)
TB2phase = TB2phase * (PI/180.0D0)


END SUBROUTINE MORPHING_PARAMETER_VALUES
