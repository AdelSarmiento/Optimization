!********+*********+*********+*********+*********+*********+*********+
!
      SUBROUTINE MAXIM(M,N,GEPS,IYFREE,GRADF,DSRCH,HESSF,X,Y,Z,ULAM, &
                      UU,XLOW,XUPP,ALFA,BETA,A,B,C,P,Q,P0,Q0)
!
!       Version "December 2006".
!    !-----------------------------------------!
!    !  The author of this subroutine is       !
!    !  Krister Svanberg <krille@math.kth.se>  !
!    !-----------------------------------------!
!
!     MAXIM solves the dual MMA subproblem.
!     The dual variables are ulam(i), i=1,..,m,
!     which are required to be non-negative.
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GRADF(1),DSRCH(1),HESSF(1),X(1),Y(1),ULAM(1), &
                UU(1),XLOW(1),XUPP(1),ALFA(1),BETA(1), &
                A(1),B(1),C(1),P(1),Q(1),P0(1),Q0(1)
      INTEGER IYFREE(1)
!
      ITR=0
      M3=3*M+30
!
      DO 10 I=1,M
      ULAM(I)=0.
      IYFREE(I)=1
 10   CONTINUE
!
      CALL GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA, &
                 A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
!
      GMX=0.
      DO 20 I=1,M
      IYFREE(I)=0
      IF(GRADF(I).GT.GEPS) IYFREE(I)=1
      IF(GRADF(I).GT.GMX) GMX=GRADF(I)
 20   CONTINUE
!
      IF(GMX.LE.GEPS) GOTO 100
!     Vi avbryter optimeringen, ulam=0 ar optimal losning.
!
 30   CONTINUE
      ITR=ITR+1
      IF(ITR.GT.M3) GOTO 100
!     Vi avbryter optimeringen pga for manga subspa-anrop.
!
      CALL SUBSPA(ITR,M,N,GEPS,F,IYFREE,GRADF,DSRCH,HESSF, &
                  X,Y,ULAM,UU,XLOW,XUPP,ALFA,BETA,A,B,C, &
                  P,Q,P0,Q0,IHITY)
!
      IF(IHITY.EQ.0) GOTO 40
!     Om ihity = 0 sa ar ulam optimal pa aktuellt underrum.
!     Om ihity > 0 sa har vi slagit i ett nytt bivillkor.
      IYFREE(IHITY)=0
      ULAM(IHITY)=0.
      GOTO 30
!
 40   CONTINUE
      CALL GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA, &
                 A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
!
      GMX=0.
      IGMX=0
      DO 50 I=1,M
      IF(IYFREE(I).EQ.1) GOTO 50
      IF(GRADF(I).LE.GMX) GOTO 50
      GMX=GRADF(I)
      IGMX=I
 50   CONTINUE
!
      IF(GMX.LE.GEPS) GOTO 100
!     Om gmx =< geps sa ar ulam optimal losning.
!     Om gmx > geps sa tar vi bort kravet att ulam(igmx)=0.
      IYFREE(IGMX)=1
      GOTO 30
!
 100  CONTINUE
!     Nu ar antingen ulam optimal losning eller itr>m3.
      CALL XYZLAM(M,N,X,Y,Z,ULAM,XLOW,XUPP,ALFA,BETA, &
                  A,B,C,P,Q,P0,Q0,IYFREE)
      IF(ITR.GT.M3) WRITE(*,911)
 911  FORMAT(' ITR GT M3 IN MAXIM')
!
      RETURN
      END
!
!********+*********+*********+*********+*********+*********+*********+
!
      SUBROUTINE SUBSPA(ITR,M,N,GEPS,F,IYFREE,GRADF,DSRCH,HESSF, &
                        X,Y,ULAM,UU,XLOW,XUPP,ALFA,BETA, &
                        A,B,C,P,Q,P0,Q0,IHITY)
!
!       Version "December 2006".
!    !-----------------------------------------!
!    !  The author of this subroutine is       !
!    !  Krister Svanberg <krille@math.kth.se>  !
!    !-----------------------------------------!
!
!    SUBSPA maximizes the dual objective function on the subspace
!    defined by ulam(i) = 0 for every i such that iyfree(i) = 0.
!    The first three iterations a steepest ascent method is used,
!    and after that a Newton method is used.
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GRADF(1),DSRCH(1),HESSF(1),X(1),Y(1), &
                ULAM(1),UU(1),XLOW(1),XUPP(1),ALFA(1),BETA(1), &
                A(1),B(1),C(1),P(1),Q(1),P0(1),Q0(1) 
      INTEGER IYFREE(1)
!
      IHITY=0
      ITESUB=0
      NYDIM=0
      DSRTOL=-0.0000001*GEPS
      CALL GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA, &
                 A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
!
      DO 10 I=1,M
      DSRCH(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 10
      NYDIM=NYDIM+1
      DSRCH(I)=GRADF(I)
 10   CONTINUE
!
      IF(NYDIM.EQ.0) GOTO 100
!     Vi avbryter med ihity = 0, ty inga variabler ulam(i) ar fria.
      ITEMAX=50+5*NYDIM
!
 15   ITESUB=ITESUB+1
!     Har startar en ny iteration.
!
      TMAX0=1.0D8
      TMAX=TMAX0
      IHITY=0
      GTD=0.
      DO 20 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 20
      GTD=GTD+GRADF(I)*DSRCH(I)
      IF(DSRCH(I).GE.0.) GOTO 20
      IF((-DSRCH(I)).LE.(ULAM(I)/TMAX0)) GOTO 20
      IF(DSRCH(I).GT.DSRTOL) GOTO 20
      T=ULAM(I)/(-DSRCH(I))
      IF(T.GE.TMAX) GOTO 20
      TMAX=T
      IHITY=I
 20   CONTINUE
      IF(TMAX.LT.0.) TMAX=0.
      IF(GTD.GT.0.) GOTO 25
      IHITY=0
      WRITE(*,912)
      GOTO 100
!     Vi avbryter med ihity = 0, ty dsrch ar ej en ascentriktning.
!
 25   CONTINUE
      CALL LINSE(M,N,ITESUB,IHITMX,IYFREE,TMAX,TOPT,ULAM,DSRCH, &
                 X,Y,UU,XLOW,XUPP,ALFA,BETA,A,B,C,P,Q,P0,Q0)
!
      DO 30 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 30
      ULAM(I)=ULAM(I)+TOPT*DSRCH(I)
      IF(ULAM(I).LT.0.) ULAM(I)=0.
 30   CONTINUE
!
      CALL GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA, &
                 A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
!
      IF(IHITMX.EQ.1.AND.IHITY.GT.0) GOTO 100
!     Vi avbryter med ihity > 0, ty vi har slagit i det tidigare
!     inaktiva bivillkoret ulam(ihity) >= 0.
      IHITY=0
      IOPT=1
      DO 40 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 40
      IF(DABS(GRADF(I)).GT.GEPS) IOPT=0
 40   CONTINUE
!
      IF(IOPT.EQ.1) GOTO 100
!     Vi avbryter med ihity = 0, ty optimal losning hittad.
      IF(ITESUB.GT.ITEMAX) GOTO 97
!     Vi avbryter med ihity = 0, ty for manga iterationer.
      IF(ITESUB.GE.3) GOTO 55
!     Om itesub>=3 sa byter vi fran steepest ascent till Newton.
!     Om itesub=<2 sa fortsatter vi med steepest ascent.
      DO 50 I=1,M
      DSRCH(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 50
      DSRCH(I)=GRADF(I)
 50   CONTINUE
      GOTO 15
!
 55   CONTINUE
!
      CALL HESSI(M,N,ULAM,HESSF,X,Y,ALFA,BETA, &
                 A,B,C,P,Q,P0,Q0,XLOW,XUPP,IYFREE)
!
      IK=0
      IKRED=0
      DO 70 K=1,M
      DO 65 I=K,M
      IK=IK+1
      IF(IYFREE(K).EQ.0) GOTO 65
      IF(IYFREE(I).EQ.0) GOTO 65
      IKRED=IKRED+1
      HESSF(IKRED)=HESSF(IK)
 65   CONTINUE
 70   CONTINUE
!
      HTRACE=0.
      IKRED=0
      ZZZZ=0.
      DO 73 K=1,NYDIM
      DO 72 I=K,NYDIM
      IKRED=IKRED+1
      IF(I.EQ.K) HTRACE=HTRACE+HESSF(IKRED)
      IF(I.EQ.K) ZZZZ=ZZZZ+1.
 72   CONTINUE
 73   CONTINUE
!
      HESMOD=0.0001*HTRACE/ZZZZ
      IF(HESMOD.LT.GEPS) HESMOD=GEPS
      IKRED=0
      DO 77 K=1,NYDIM
      DO 76 I=K,NYDIM
      IKRED=IKRED+1
      IF(I.EQ.K) HESSF(IKRED)=HESSF(IKRED)+HESMOD
 76   CONTINUE
 77   CONTINUE
!
      CALL LDLFAC(NYDIM,GEPS,HESSF,UU)
!
      IRED=0
      DO 79 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 79
      IRED=IRED+1
      UU(IRED)=GRADF(I)
 79   CONTINUE
!
      CALL LDLSOL(NYDIM,UU,HESSF,DSRCH)
!
      DO 80 I=1,M
      UU(I)=DSRCH(I)
 80   CONTINUE
!
      IRED=0
      DO 85 I=1,M
      DSRCH(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 85
      IRED=IRED+1
      DSRCH(I)=UU(IRED)
 85   CONTINUE
!
      GOTO 15
!
 97   CONTINUE
      WRITE(*,911)
!
 100  CONTINUE
!
      DO 110 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 110
      IF(ULAM(I).LT.0.) ULAM(I)=0.
 110  CONTINUE
!
 911  FORMAT(' ITESUB GT ITEMAX IN SUBSPA')
 912  FORMAT(' GTD LE 0 IN SUBSPA')
!
      RETURN
      END
!
!********+*********+*********+*********+*********+*********+*********+
!
      SUBROUTINE LINSE(M,N,ITESUB,IHITMX,IYFREE,TMAX,TOPT,ULAM,DSRCH, &
                       X,Y,UU,XLOW,XUPP,ALFA,BETA,A,B,C,P,Q,P0,Q0)
!
!       Version "December 2006".
!    !-----------------------------------------!
!    !  The author of this subroutine is       !
!    !  Krister Svanberg <krille@math.kth.se>  !
!    !-----------------------------------------!
!
!     LINSE makes an approximate line search (maximization) in the
!     direction DSRCH from the point ULAM.
!     Main input:  ULAM, DSRCH, TMAX.
!     Main output: TOPT, IHITMX.
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ULAM(1),DSRCH(1),X(1),Y(1),UU(1),XLOW(1),XUPP(1), &
                ALFA(1),BETA(1),A(1),B(1),C(1),P(1),Q(1),P0(1),Q0(1)
      INTEGER IYFREE(1)
!
      ITT1=0
      ITT2=0
      ITT3=0
!
      CALL LINDER(M,N,TMAX,DFDTMX,ULAM,DSRCH,X,Y,UU,XLOW,XUPP, &
                  ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(DFDTMX.GE.0.) GOTO 80
!     Linjesokningen klar. Optimalt steg ar T=TMAX. IHITMX=1.
      IF(TMAX.GT.1.) GOTO 40
      T2=TMAX
!
 30   CONTINUE
!     Nu sker en upprepad minskning av steget.
      ITT1=ITT1+1
      IF(ITT1.GT.13) GOTO 90
      T1=T2/2.
      IF(ITESUB.LE.3) T1=T2/16.
      CALL LINDER(M,N,T1,DFDT1,ULAM,DSRCH,X,Y,UU,XLOW,XUPP, &
                  ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(DFDT1.GT.0.) GOTO 60
      T2=T1
      GOTO 30
!
 40   CONTINUE
!     Nu testas enhetssteget, dvs T=1.
      T1=1.
      T2=T1
      CALL LINDER(M,N,T1,DFDT1,ULAM,DSRCH,X,Y,UU,XLOW,XUPP, &
                  ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(ITESUB.GE.6.AND.DFDT1.GE.0.) GOTO 90
!     Linjesokningen klar. Enhetssteget duger. T=1, IHITMX=0.
      IF(ITESUB.LE.5.AND.DFDT1.GT.0.) GOTO 50
!     Enhetssteget ar for kort.
      GOTO 30
!     Enhetssteget ar for langt.
!
 50   ITT2=ITT2+1
!     Nu sker en upprepad okning av steget.
      IF(ITT2.GT.10) GOTO 90
      T2=2.*T1
      IF(ITESUB.LE.3) T2=16.*T1
      IF(T2.LT.TMAX) GOTO 55
      T2=TMAX
      GOTO 60
 55   CALL LINDER(M,N,T2,DFDT2,ULAM,DSRCH,X,Y,UU,XLOW,XUPP, &
                  ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(DFDT2.LE.0.) GOTO 60
      T1=T2
      GOTO 50
!
 60   CONTINUE
!     Nu sker en upprepad krympning av intervallet T1,T2.
      SQT1=DSQRT(T1)
      SQT2=DSQRT(T2)
 62   ITT3=ITT3+1
      IF(ITT3.GT.10) GOTO 90
      TM=SQT1*SQT2
      CALL LINDER(M,N,TM,DFDTM,ULAM,DSRCH,X,Y,UU,XLOW,XUPP, &
                  ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(DFDTM.GT.0.) GOTO 65
      T2=TM
      TKVOT=T1/T2
      IF(TKVOT.GT.0.97) GOTO 90
!     Linjesokningen klar. T1 ar approx optimal. IHITMX=0.
      SQT2=DSQRT(T2)
      GOTO 62
 65   T1=TM
      TKVOT=T1/T2
      IF(TKVOT.GT.0.97) GOTO 90
!     Linjesokningen klar. T1 ar approx optimal. IHITMX=0.
      SQT1=DSQRT(T1)
      GOTO 62
!
 80   TOPT=TMAX
      IHITMX=1
      GOTO 100
 90   TOPT=T1
      IHITMX=0
      IF(ITT1.GT.13) WRITE(*,911)
      IF(ITT2.GT.10) WRITE(*,912)
      IF(ITT3.GT.10) WRITE(*,913)
 911  FORMAT(' ITT1 GT 13 in LINSE')
 912  FORMAT(' ITT2 GT 10 in LINSE')
 913  FORMAT(' ITT3 GT 10 in LINSE')
 100  CONTINUE
!
      RETURN
      END
!
!********+*********+*********+*********+*********+*********+*********+
