SUBROUTINE INPUT
USE global ! module for global variables
USE iga_global ! module for global variables (igalib)

INTEGER :: i


 NR = 24 ! number of chordwise panels
 NC = 20 ! number of spanwise panels

 NEL = NR*NC ! number of panels/elements

 PI = ACOS(-1.0D0)

 spanL = 1.0D0 ! chordL (m)
 chordL = 6.0D0 ! spanL (m)
 
 AR = chordL/spanL ! aspect ratio

 V_ref = 10.0D0 ! freestream velocity (m/s)
 
 Omega = 2.0D0*1.0D0! frequency of flapping (rad/s)

 Red_Freq = Omega*(spanL/2.0D0)/V_ref ! reduced frequency

 NST = 120 ! number of time step per cycle

 t_step = (2.0D0*PI/Omega)/REAL(NST) ! time step

 N_steps = 151 ! number of time steps

 N_W = 10 ! number of deforming wake rows (only update the newest N_W wake rings)

 N_K = 7 ! number of knots
 
 !!!!!!!!!!! bisplines geometry !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !PX = 6
 PY = 3 ! (3->cubic)
 
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 CALL MEMORY_ALLOCATION

 DO i= 1, N_steps
    t(i) = REAL(i-1)*t_step ! time vector
 END DO

 N_Avg_Begin = 31
 
 Cut_off = 0.000001D0 ! cutoff radius (Biot-Savart law)

 rho = 1.225D0 ! fluid density 

 WRITE (*,*) "Input data:"
 WRITE (*,*) "Number of time steps:", N_steps
 WRITE (*,*) "Chord:", spanL
 WRITE (*,*) "Aspect ratio:", AR
 WRITE (*,*) "Number of rows:", NR
 WRITE (*,*) "Number of columns:", NC
 WRITE (*,*) "Reduced frequency:", Red_freq
 WRITE (*,*)
 
END SUBROUTINE INPUT

 

