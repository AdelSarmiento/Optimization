SUBROUTINE TIME_VARYING_MORPHING
USE global ! module for global variables

INTEGER :: i

Ones(:) = 1.0D0 

twist1(:) = TW1amp * SIN(Omega*t(:)+TW1phase*Ones(:))

bend1(:) = TB1amp * SIN(Omega*t(:)+TB1phase*Ones(:))

twist2(:) = TW2amp * SIN(Omega*t(:)+TW2phase*Ones(:))

bend2(:) = TB2amp * SIN(Omega*t(:)+TB2phase*Ones(:))

DO i=1, N_steps
   fore(i) = -(1.0D0/chordL)*(bend1(i)+bend2(i))**2
END DO


END SUBROUTINE TIME_VARYING_MORPHING
