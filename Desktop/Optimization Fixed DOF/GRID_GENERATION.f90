SUBROUTINE GRID_GENERATION(GEOMVAL2)
USE REAL_PRECISION
USE IGA
USE bspline
USE Srf
USE global ! module for global variables
USE iga_global ! module for global variables (igalib)
INTEGER :: i, j
REAL(KIND = R8) :: uux, vvy, DEL
REAL(KIND = R8), DIMENSION(18) :: GEOMVAL2
  
REAL(KIND = R8), DIMENSION(:), ALLOCATABLE :: UX0, UY0
REAL(kind=8), DIMENSION(:,:,:), ALLOCATABLE :: GEOM0

REAL(KIND = R8), DIMENSION(:), ALLOCATABLE :: UX1, UY1
REAL(kind=8), DIMENSION(:,:,:), ALLOCATABLE :: GEOM1

REAL(KIND = R8), DIMENSION(:), ALLOCATABLE :: UX2, UY2
REAL(kind=8), DIMENSION(:,:,:), ALLOCATABLE :: GEOM2
REAL(kind=8), DIMENSION(:,:,:), ALLOCATABLE :: GEOM2i

REAL(KIND = R8), DIMENSION(:), ALLOCATABLE :: UYP, Xp

INTEGER :: MX0, NDX0, PX0
INTEGER :: MY0, NDY0, PY0

INTEGER :: MX1, NDX1, PX1
INTEGER :: MY1, NDY1, PY1

INTEGER :: MX2, NDX2, PX2
INTEGER :: MY2, NDY2, PY2

INTEGER :: kp, ep, cp

 CHARACTER*50 :: fname1, fname2, fname3
 CHARACTER*50 :: fname4, fname5, fname6
 CHARACTER*50 :: fname7, fname8, fname9
 CHARACTER*50 :: fname0

      WRITE(fname0,340)PY 
340 FORMAT('ConvergenceWR_P',i1,'.dat')

OPEN(unit=1111,file=fname0,form='FORMATTED')

      WRITE(fname1,341)PY 
341 FORMAT('GridxWR_P',i1,'.dat')

     WRITE(fname2,342)PY 
342 FORMAT('GridyWR_P',i1,'.dat')

     WRITE(fname3,343)PY 
343 FORMAT('GridzWR_P',i1,'.dat')


OPEN(unit=111,file=fname1,form='FORMATTED')
OPEN(unit=112,file=fname2,form='FORMATTED')
OPEN(unit=113,file=fname3,form='FORMATTED')

      WRITE(fname4,344)PY 
344 FORMAT('CTR_PTSxWR_P',i1,'.dat')

     WRITE(fname5,345)PY 
345 FORMAT('CTR_PTSyWR_P',i1,'.dat')

     WRITE(fname6,346)PY 
346 FORMAT('CTR_PTSzWR_P',i1,'.dat')

OPEN(unit=211,file=fname4,form='FORMATTED')
OPEN(unit=212,file=fname5,form='FORMATTED')
OPEN(unit=213,file=fname6,form='FORMATTED')


      WRITE(fname7,347)PY 
347 FORMAT('KNT_VECTXWR_P',i1,'.dat')

     WRITE(fname8,348)PY 
348 FORMAT('KNT_VECTYWR_P',i1,'.dat')


OPEN(unit=311,file=fname7,form='FORMATTED')
OPEN(unit=312,file=fname8,form='FORMATTED')

      WRITE(fname9,349)PY 
349 FORMAT('AEROLOADSWR_P',i1,'.dat')

OPEN(unit=1,file=fname9,form='FORMATTED')

! Grid generation for NACA 83XX (for an uncambered wing set Z0 to 0) 
PX0 = 6
PY0 = 1
MX0 = 2*(PX0+1) - 1
MY0 = 2*(PY0+1) - 1
NDX0 = PX0
NDY0 = PY0
ALLOCATE(UX0(0:MX0))
ALLOCATE(UY0(0:MY0))
ALLOCATE(GEOM0(1:3,0:NDY0,0:NDX0))
 CALL KnotVector(1,PX0,0,0.0_R8,spanL,.FALSE.,MX0,UX0)
 CALL KnotVector(1,PY0,0,0.0_R8,chordL/2.0_R8,.FALSE.,MY0,UY0)
 

 
DO i=0, 6
   GEOM0(1,:,i) = i*spanL/6.0_R8
END DO
DO j=0, 1
   GEOM0(2,j,:) = j*chordL/2.0_R8
END DO
GEOM0(3,:,0) = 0.0_R8
GEOM0(3,:,1) = 0.089250000000001564_R8
GEOM0(3,:,2) = 0.12646666666666181_R8
GEOM0(3,:,3) = 0.034380000000007863_R8
GEOM0(3,:,4) = 0.091453333333326003_R8
GEOM0(3,:,5) = 0.031200000000003506_R8
GEOM0(3,:,6) = 0.0_R8


! Degree elevate in the y direction
PX1 = 6
PY1 = PY
MX1 = MX0
MY1 = 2*(PY1+1) - 1
NDX1 = PX1
NDY1 = PY1 
ALLOCATE(UX1(0:MX1))
ALLOCATE(UY1(0:MY1))
ALLOCATE(GEOM1(1:3,0:NDY1,0:NDX1))
 CALL DegreeElevate2(3,NDX0,PX0,UX0,NDY0,PY0,UY0,GEOM0,0,PY1-PY0,& 
		NDX1,UX1,NDY1,UY1,GEOM1)
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PX2 = PX1
PY2 = PY1
kp = 5 -PY2
ep = kp + 1
 cp = PY2-1
MX2 = 2*(PX2+1) - 1
MY2 = 2*(PY2+1) + kp -1
NDX2 = PX2
NDY2 = PY2 + kp 
		

ALLOCATE(UYP(1:kp))		

DO i=1, kp
   UYP(i) = 3.0_R8*REAL(i, KIND = R8)/REAL(ep, KIND = R8)
END DO


ALLOCATE(UX2(0:MX2))
ALLOCATE(UY2(0:MY2))
ALLOCATE(GEOM2(1:3,0:NDY2,0:NDX2))
ALLOCATE(GEOM2i(1:3,0:NDY2,0:NDX2))
ALLOCATE(Xp(0:0))
DO i=0, 0
   Xp(i) = 0.0_R8
END DO

 call RefineKnotVector2(3,NDX1,PX1,UX1,NDY1,PY1,UY1,GEOM1,-1,Xp,kp-1,& 
                        UYP,UX2,UY2,GEOM2)


DO j=1, NDY2
   GEOM2(1,j,0) = GEOM2(1,j,0) + GEOMVAL2(j)
   GEOM2(1,j,NDX2) = GEOM2(1,j,NDX2) + GEOMVAL2(NDY2+j)   
END DO

GEOM2(2,NDY2,:) = GEOM2(2,NDY2,:) + GEOMVAL2(2*NDY2+1)

GEOM2i(:,:,:) = GEOM2(:,:,:)

GEOM2(1,0,0) = GEOM2(1,0,0) + GEOMVAL2(2*NDY2+2)

DO i=1, NDX2
   DO j=0, NDY2
      GEOM2(1,j,i) = (GEOM2i(1,0,i)-GEOM2i(1,0,0))/spanL & 
	                   *(GEOM2(1,j,NDX1)-GEOM2(1,j,0))+GEOM2(1,j,0)
   END DO
END DO

! Sample the surface

DX = spanL/REAL(NR, KIND = R8) ! element choordwise length
DY = (chordL/2.0_R8)/REAL(NC, KIND = R8) ! element spanwise length

DO i=1, NR+1
   DO j=1, NC+1
      uux = REAL(i-1, KIND = R8)*DX
      vvy = REAL(j-1, KIND = R8)*DY

      CALL SurfacePoint(3,NDX2,PX2,UX2,NDY2,PY2,UY2,GEOM2,uux,vvy,S)
      X0(i,j) = S(1) 
      Y0(i,j) = S(2)
      Z0(i,j) = S(3)
      !WRITE(*,*) X0(i,j), Y0(i,j), Z0(i,j)

   END DO
END DO

 CALL ANGLES


DO i=1, NR+1
   WRITE(111,1001) (X0(i,j), j=1,NC+1)
END DO

DO i=1, NR+1
   WRITE(112,1001) (Y0(i,j), j=1,NC+1)
END DO

DO i=1, NR+1
   WRITE(113,1001) (Z0(i,j), j=1,NC+1)
END DO

DO i=0, NDX2
   WRITE(211,1001) (GEOM2(1,j,i), j=0,NDY2)
END DO

DO i=0, NDX2
   WRITE(212,1001) (GEOM2(2,j,i), j=0,NDY2)
END DO

DO i=0, NDX2
   WRITE(213,1001) (GEOM2(3,j,i), j=0,NDY2)
END DO

WRITE(311,1001) (UX2(j), j=0,MX2)

WRITE(312,1001) (UY2(j), j=0,MY2)

 1001     FORMAT(1x, 101(1pe12.5,2x))
 
 
 CLOSE(111)
 CLOSE(112)
 CLOSE(113)
 
 CLOSE(211)
 CLOSE(212)
 CLOSE(213)
 
 CLOSE(311)
 CLOSE(312)

END SUBROUTINE GRID_GENERATION
