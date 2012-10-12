PROGRAM MAIN_UVLM
USE REAL_PRECISION
USE global ! module for global variables
USE iga_global ! module for global variables (igalib)

IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION XVAL(100),XOLD1(100),XOLD2(100),XMMA(100), &
                XMIN(100), XMAX(100), XLOW(100),XUPP(100), &
                ALFA(100), BETA(100),DF0DX(100), &
                A(60),B(60),C(60),Y(60),RAA(60),ULAM(60), &
                FVAL(60),FAPP(60),FNEW(60),FMAX(60), &
                DFDX(6000),P(6000),Q(6000),P0(100),Q0(100), &
                UU(60),GRADF(60),DSRCH(60),HESSF(1830) 
      INTEGER IYFREE(60)
  

 CALL INPUT
 
 
 CALL INITI(M,N,GEPS,XVAL,XMIN,XMAX,FMAX,A,C)
 
      IF(M.EQ.0) GOTO 100
      IF(N.EQ.0) GOTO 100
      INNMAX=15
      MAXITE=1
      ITER=0
      ITE=0
      INNER=0
!
!  The USER should now calculate function values at XVAL.
!  The result should be put in F0VAL and FVAL.
!
      
      CALL FUNC1(M,N,XVAL,F0VAL,FVAL)
!
!  The USER may now write the current (starting) solution.
! 
      CALL OUTXOF(ITER,INNER,M,N,XVAL,F0VAL,FVAL,FMAX)

 30   CONTINUE
!
      ITER=ITER+1
      ITE=ITE+1
!
!  The USER should now calculate function values and gradients
!  at XVAL. The result should be put in F0VAL,DF0DX,FVAL,DFDX.
!
      
      CALL FUNC2(M,N,XVAL,F0VAL,DF0DX,FVAL,DFDX)
!
!
!  RAA0,RAA,XLOW,XUPP,ALFA and BETA are calculated.
!
      CALL RAASTA(M,N,RAA0,RAA,XMIN,XMAX,DF0DX,DFDX)
      CALL ASYMPG(ITER,M,N,XVAL,XMIN,XMAX,XOLD1,XOLD2, &
                  XLOW,XUPP,ALFA,BETA)
!      CALL ASYMPG(ITER,M,N,XVAL,XMIN,XMAX,XOLD1,XOLD2,
!     1            XLOW,XUPP,ALFA,BETA,RAA0,RAA)
!
!  The inner iterative process starts.
!
      INNER=0
 40   CONTINUE
!
!  The subproblem is generated and solved.
!
      CALL MMASUG(ITER,M,N,GEPS,IYFREE,XVAL,XMMA, &
                  XMIN,XMAX,XLOW,XUPP,ALFA,BETA, &
                  A,B,C,Y,Z,RAA0,RAA,ULAM, &
                  F0VAL,FVAL,F0APP,FAPP,FMAX,DF0DX,DFDX, &
                  P,Q,P0,Q0,UU,GRADF,DSRCH,HESSF)
!
!  The USER should now calculate function values at XMMA.
!  The result should be put in F0NEW and FNEW.
!
      CALL FUNC1(M,N,XMMA,F0NEW,FNEW)
!
      IF(INNER.GE.INNMAX) GOTO 60
!
!  It is checked if the approximations were conservative.
!
      CALL CONSER(M,ICONSE,GEPS,F0NEW,F0APP,FNEW,FAPP)
      IF(ICONSE.EQ.1) GOTO 60
!
!  The approximations were not conservative, so RAA0 and RAA
!  are updated and one more inner iteration is started.
!
      INNER=INNER+1
      CALL RAAUPD(M,N,GEPS,XMMA,XVAL,XMIN,XMAX,XLOW,XUPP, &
                 F0NEW,FNEW,F0APP,FAPP,RAA0,RAA)
      GOTO 40
!
 60   CONTINUE
!
!  The inner iterative process has terminated, which means
!  that an outer iteration has been completed.
!  The variables are updated so that XVAL stands for the new
!  outer iteration point. The fuction values are also updated.
!
      CALL XUPDAT(N,ITER,XMMA,XVAL,XOLD1,XOLD2)
      CALL FUPDAT(M,F0NEW,FNEW,F0VAL,FVAL)
!
!  The USER may now write the current solution.
! 
      CALL OUTXOF(ITER,INNER,M,N,XVAL,F0VAL,FVAL,FMAX)
!
!  One more outer iteration is started as long as
!  ITE is less than MAXITE:
!
      IF(ITE.LT.MAXITE) GOTO 30
!
      WRITE(*,90)
 90   FORMAT(' How many more iterations? (0 to stop)')
      READ(*,*) MAXITE
      IF(MAXITE.EQ.0) GOTO 100
      ITE=0
      GOTO 30
!
 100  CONTINUE
!
      STOP



END PROGRAM MAIN_UVLM
