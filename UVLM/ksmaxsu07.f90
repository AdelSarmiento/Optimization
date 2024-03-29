!********+*********+*********+*********+*********+*********+*********+
!
      SUBROUTINE XUPDAT(N,ITER,XMMA,XVAL,XOLD1,XOLD2)
!
!       Version "December 2006".
!    !-----------------------------------------!
!    !  The author of this subroutine is       !
!    !  Krister Svanberg <krille@math.kth.se>  !
!    !-----------------------------------------!
!
!     XUPDAT updates XVAL, XOLD1 and if possible XOLD2.
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XMMA(1),XVAL(1),XOLD1(1),XOLD2(1)
!
      DO 70 J=1,N
      IF(ITER.GE.2) XOLD2(J)=XOLD1(J)
      XOLD1(J)=XVAL(J)
      XVAL(J)=XMMA(J)
 70   CONTINUE
!
      RETURN
      END
!
!********+*********+*********+*********+*********+*********+*********+
!
      SUBROUTINE XYZLAM(M,N,X,Y,Z,ULAM,XLOW,XUPP,ALFA,BETA, &
                       A,B,C,P,Q,P0,Q0,IYFREE)
!
!       Version "December 2006".
!    !-----------------------------------------!
!    !  The author of this subroutine is       !
!    !  Krister Svanberg <krille@math.kth.se>  !
!    !-----------------------------------------!
!
!     XYZLAM calculates the X,Y,Z that minimize the Lagrange
!     function, given the vector ULAM of Lagrange multipliers.
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(1),Y(1),ULAM(1),XLOW(1),XUPP(1),ALFA(1),BETA(1), &
                A(1),B(1),C(1),P(1),Q(1),P0(1),Q0(1)
      INTEGER IYFREE(1)
!
      DO 10 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 10
      IF(ULAM(I).LT.0.) ULAM(I)=0.
 10   CONTINUE
!
      DO 30 J=1,N
      PJ=P0(J)
      QJ=Q0(J)
      MJ1=M*(J-1)
!
      DO 20 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 20
      PIJ=P(MJ1+I)
      QIJ=Q(MJ1+I)
      PJ=PJ+ULAM(I)*PIJ
      QJ=QJ+ULAM(I)*QIJ
 20   CONTINUE
!
      SRPJ=DSQRT(PJ)
      SRQJ=DSQRT(QJ)
      XJ=(SRPJ*XLOW(J)+SRQJ*XUPP(J))/(SRPJ+SRQJ)
      IF(XJ.LT.ALFA(J)) XJ=ALFA(J)
      IF(XJ.GT.BETA(J)) XJ=BETA(J)
      X(J)=XJ
 30   CONTINUE
!
      UA=0.
      DO 40 I=1,M
      Y(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 40
      UA=UA+ULAM(I)*A(I)
      YI=ULAM(I)-C(I)
      IF(YI.GT.0.) Y(I)=YI
 40   CONTINUE
!
      Z=0.
      UA1=UA-1.
      IF(UA1.GT.0.) Z=10.*UA1
!
      RETURN
      END
!
!********+*********+*********+*********+*********+*********+*********+
!
      SUBROUTINE GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA, &
                       A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
!
!       Version "December 2006".
!    !-----------------------------------------!
!    !  The author of this subroutine is       !
!    !  Krister Svanberg <krille@math.kth.se>  !
!    !-----------------------------------------!
!
!     GRADI calculates the gradient GRADF of the dual
!     objective function, given the vector ULAM of dual variables.
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(1),Y(1),ULAM(1),XLOW(1),XUPP(1),ALFA(1),BETA(1), &
                A(1),B(1),C(1),P(1),Q(1),P0(1),Q0(1),GRADF(1)
      INTEGER IYFREE(1)
!
      CALL XYZLAM(M,N,X,Y,Z,ULAM,XLOW,XUPP,ALFA,BETA, &
                  A,B,C,P,Q,P0,Q0,IYFREE)
!
      DO 10 I=1,M
      GRADF(I)=-B(I)-Y(I)-A(I)*Z
 10   CONTINUE
!
      DO 30 J=1,N
      MJ1=M*(J-1)
      UJXJ=XUPP(J)-X(J)
      XJLJ=X(J)-XLOW(J)
!
      DO 20 I=1,M
      PIJ=P(MJ1+I)
      QIJ=Q(MJ1+I)
      GRADF(I)=GRADF(I)+PIJ/UJXJ+QIJ/XJLJ
 20   CONTINUE
 30   CONTINUE
!
      RETURN
      END
!
!********+*********+*********+*********+*********+*********+*********+
!
      SUBROUTINE LINDER(M,N,T,DFDT,ULAM,DSRCH,X,Y,UU,XLOW,XUPP, &
                        ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
!
!       Version "December 2006".
!    !-----------------------------------------!
!    !  The author of this subroutine is       !
!    !  Krister Svanberg <krille@math.kth.se>  !
!    !-----------------------------------------!
!
!     LINDER calculates the scalar product DFDT of GRADF and DSRCH
!     (= the directional derivative) at the point ULAM + T*DSRCH.
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ULAM(1),DSRCH(1),X(1),Y(1),UU(1), &
                XLOW(1),XUPP(1),ALFA(1),BETA(1), &
                A(1),B(1),C(1),P(1),Q(1),P0(1),Q0(1)
      INTEGER IYFREE(1)
!
      DO 10 I=1,M
      UU(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 10
      UU(I)=ULAM(I)+T*DSRCH(I)
      IF(UU(I).LT.0.) UU(I)=0.
 10   CONTINUE
!
      CALL XYZLAM(M,N,X,Y,Z,UU,XLOW,XUPP,ALFA,BETA, &
                  A,B,C,P,Q,P0,Q0,IYFREE)
!
      DO 20 I=1,M
      UU(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 20
      UU(I)=-B(I)-Y(I)-A(I)*Z
   20 CONTINUE
!
      DO 40 J=1,N
      MJ1=M*(J-1)
      UJXJ=XUPP(J)-X(J)
      XJLJ=X(J)-XLOW(J)
      DO 30 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 30
      PIJ=P(MJ1+I)
      QIJ=Q(MJ1+I)
      UU(I)=UU(I)+PIJ/UJXJ+QIJ/XJLJ
   30 CONTINUE
   40 CONTINUE
!
      DFDT=0.
      DO 50 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 50
      DFDT=DFDT+UU(I)*DSRCH(I)
   50 CONTINUE
!
      RETURN
      END
!
!********+*********+*********+*********+*********+*********+*********+
!
      SUBROUTINE HESSI(M,N,ULAM,HESSF,X,Y,ALFA,BETA, &
                       A,B,C,P,Q,P0,Q0,XLOW,XUPP,IYFREE)
!
!       Version "December 2006".
!    !-----------------------------------------!
!    !  The author of this subroutine is       !
!    !  Krister Svanberg <krille@math.kth.se>  !
!    !-----------------------------------------!
!
!   HESSI calculates HESSF = minus the reduced Hessian matrix of the
!   dual objective function, given the vector ULAM of dual variables.
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ULAM(1),HESSF(1),X(1),Y(1),ALFA(1),BETA(1),A(1), &
                B(1),C(1),P(1),Q(1),P0(1),Q0(1),XLOW(1),XUPP(1)
      INTEGER IYFREE(1)
!
      IK=0
      DO 12 K=1,M
      DO 11 I=K,M
      IK=IK+1
      HESSF(IK)=0.
 11   CONTINUE
 12   CONTINUE
!
      ULAMTA=0.
      II=1
      DO 15 I=1,M
      IF(I.GT.1) II=II+M+2-I
      IF(IYFREE(I).EQ.0) GOTO 15
      IF(ULAM(I).LT.0.) ULAM(I)=0.
      IF(ULAM(I).GT.C(I)) HESSF(II)=1.
      ULAMTA=ULAMTA+ULAM(I)*A(I)
 15   CONTINUE
!
      IF(ULAMTA.LE.1.) GOTO 40
!
      KK=1
      DO 30 K=1,M
      IF(K.GT.1) KK=KK+M+2-K
      IF(IYFREE(K).EQ.0) GOTO 30
      DO 20 I=K,M
      IF(IYFREE(I).EQ.0) GOTO 20
      IK=KK+I-K
      HESSF(IK)=HESSF(IK)+10.*A(I)*A(K)
 20   CONTINUE
 30   CONTINUE
!
 40   CONTINUE
!
      DO 100 J=1,N
      PJ=P0(J)
      QJ=Q0(J)
      MJ1=M*(J-1)
      DO 50 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 50
      PIJ=P(MJ1+I)
      QIJ=Q(MJ1+I)
      PJ=PJ+ULAM(I)*PIJ
      QJ=QJ+ULAM(I)*QIJ
 50   CONTINUE
!
      SRPJ=DSQRT(PJ)
      SRQJ=DSQRT(QJ)
      XJ=(SRPJ*XLOW(J)+SRQJ*XUPP(J))/(SRPJ+SRQJ)
      IF(XJ.GE.BETA(J)) GOTO 100
      IF(XJ.LE.ALFA(J)) GOTO 100
!
      UJXJ=XUPP(J)-XJ
      XJLJ=XJ-XLOW(J)
      UJXJ2=UJXJ**2
      XJLJ2=XJLJ**2
      RR=2.*PJ/UJXJ**3+2.*QJ/XJLJ**3
!
      KK=1
      DO 80 K=1,M
      IF(K.GT.1) KK=KK+M+2-K
      IF(IYFREE(K).EQ.0) GOTO 80
      PKJ=P(MJ1+K)
      QKJ=Q(MJ1+K)
      TTK=PKJ/UJXJ2-QKJ/XJLJ2
      DO 70 I=K,M
      IF(IYFREE(I).EQ.0) GOTO 70
      IK=KK+I-K
      PIJ=P(MJ1+I)
      QIJ=Q(MJ1+I)
      TTI=PIJ/UJXJ2-QIJ/XJLJ2
      HESSF(IK)=HESSF(IK)+TTI*TTK/RR
 70   CONTINUE
 80   CONTINUE
!
 100  CONTINUE
!
      RETURN
      END
!
!********+*********+*********+*********+*********+*********+*********+
!
      SUBROUTINE LDLFAC(N,EPS,ADL,E)
!
!       Version "December 2006".
!    !-----------------------------------------!
!    !  The author of this subroutine is       !
!    !  Krister Svanberg <krille@math.kth.se>  !
!    !-----------------------------------------!
!
!    LDLFAC makes a factorization of a given symmetric matrix A.
!    If A is positive definite, then A = L*D*LT.
!    If A is not positive definite, then A + E = L*D*LT,
!    where E is a positive semidefinite diagonal matrix such that
!    A + E is positive definite.
!    On entry, ADL defines the given matrix A.
!    On leave, ADL defines the calculated matrices D and L.
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ADL(1),E(1)
!
      JJ=1
!
      DO 100 J=1,N
      E(J)=0.
      IF(J.GT.1) JJ=JJ+N+2-J
      IF(J.EQ.1) GOTO 25
      KK=JJ
      JK=JJ
      DO 20 L=1,J-1
      K=J-L
      JK=JK-N+K
      KK=KK-N+K-1
      ADL(JJ)=ADL(JJ)-ADL(KK)*ADL(JK)*ADL(JK)
 20   CONTINUE
 25   IF(ADL(JJ).GE.EPS) GOTO 35
      E(J)=EPS-ADL(JJ)
      ADL(JJ)=EPS
 35   IF(J.EQ.N) GOTO 100
      IJ=JJ
      DO 50 I=J+1,N
      IJ=IJ+1
      ADLIJ=ADL(IJ)
      IF(J.EQ.1) GOTO 45
      IK=IJ
      JK=JJ
      KK=JJ
      DO 40 L=1,J-1
      K=J-L
      IK=IK-N+K
      JK=JK-N+K
      KK=KK-N+K-1
      ADLIJ=ADLIJ-ADL(KK)*ADL(IK)*ADL(JK)
 40   CONTINUE
 45   ADL(IJ)=ADLIJ/ADL(JJ)
 50   CONTINUE
 100  CONTINUE
!
      RETURN
      END
!
!********+*********+*********+*********+*********+*********+*********+
!
      SUBROUTINE LDLSOL(N,B,DL,X)
!
!       Version "December 2006".
!    !-----------------------------------------!
!    !  The author of this subroutine is       !
!    !  Krister Svanberg <krille@math.kth.se>  !
!    !-----------------------------------------!
!
!     LDLSOL solves a system of linear equations: A*X = B,
!     where A has already been factorized as L*D*Ltranspose.
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(1),DL(1),X(1)
!
      JJ=1
!
      DO 30 J=1,N
      X(J)=B(J)
      IF(J.EQ.1) GOTO 30
      JJ=JJ+N+2-J
      JK=JJ
      DO 20 L=1,J-1
      K=J-L
      JK=JK-N+K
      X(J)=X(J)-DL(JK)*X(K)
 20   CONTINUE
 30   CONTINUE
!
      JJ=1
      DO 40 J=1,N
      IF(J.GT.1) JJ=JJ+N+2-J
      X(J)=X(J)/DL(JJ)
 40   CONTINUE
!
      DO 60 L=1,N-1
      J=N-L
      JJ=JJ-N+J-1
      KJ=JJ
      DO 50 K=J+1,N
      KJ=KJ+1
      X(J)=X(J)-DL(KJ)*X(K)
 50   CONTINUE
 60   CONTINUE
!
      RETURN
      END
!
!********+*********+*********+*********+*********+*********+*********+
