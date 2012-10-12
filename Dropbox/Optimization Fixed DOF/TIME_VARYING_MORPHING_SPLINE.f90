SUBROUTINE TIME_VARYING_MORPHING_SPLINE
USE global ! module for global variables

INTEGER :: i

Ones(:) = 1.0D0 

twist1(:) = Y_spline(:)*PI/180.0D0

bend1(:) = Y3_spline(:)

twist2(:) = Y2_spline(:)*PI/180.0D0

bend2(:) = Y4_spline(:)

DO i=1, N_steps
   fore(i) = -(1.0D0/chordL)*(bend1(i)**2+bend2(i)**2)
   !WRITE(2,50) t(i)/(2.0D0*PI/omega)-0.25D0, Y_spline(i), Y2_spline(i), Y3_spline(i), Y4_spline(i)
END DO

50     FORMAT(1x, 5(1pe12.5,2x))

END SUBROUTINE TIME_VARYING_MORPHING_SPLINE
