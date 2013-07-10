C-----------------------------------------------------------------------
C
C#######################################################################
      SUBROUTINE INVERR(ALPHA,BETA,V1,V2,G,MEV,VS,MP,ERR,NDIS)
C#######################################################################
C
C calculate error estimates for computed isolated 'good' eigenvalues of
C a tridiagonal Lanczos matrix (with respect to eigenvalues of the
C original matrix) using inverse iteration
C
C this routine is copied from C&W with the following changes:
C omission of the final part of the routine and elimination of the NG
C  variable
C omission of all C&W's comment lines and WRITE statements (and the
C  IWRITE variable)
C addition of my comment lines
C minor formatting and FORTRAN style changes
C IMPLICIT typing (the variable NORM is renamed XNORM and G becomes
C  DOUBLE PRECISION)
C explicit dimensions for all arrays
C reordering of call parameters
C replacement of top section of G array by the ERR array, added to the
C  call parameters
C storage of computed errors in ERR (reset internally) in a manner
C  comensurate with VS and MP instead of according to NISO, which is
C  eliminated
C elimination of (useless) call parameters MMB and N
C replacement of call parameters EPS, IKL and IT by internal settings
C replacement of FINPRO by SCALAR for scalar product
C addition of explanatory writing/timing and the QWRITE parameter
C
C meaning of call parameters:
C ALPHA: diagonal matrix elements, unchanged on exit
C BETA: subdiagonal matrix elements, unchanged on exit
C V1,V2,G: work spaces
C MEV: order of matrix, unchanged on exit
C VS: computed distinct eigenvalues in increasing order, unchanged on
C  exit
C MP: corresponding multiplicities (except that 0 means spurious and -1
C  means good but close to a spurious - see ISOEV), unchanged on exit
C ERR: calculated error estimates (0 for levels for which MP<>1 and
C  negative if the iteration limit was reached
C NDIS: no. of computed distinct eigenvalues, unchanged on exit
C QWRITE: if true, do timing and write out
C
      IMPLICIT DOUBLE PRECISION (A-H,O-P,R-Z),LOGICAL (Q)
      DIMENSION ALPHA(MEV),BETA(MEV+1),V1(MEV),V2(MEV+1),G(MEV),
     1          VS(NDIS),MP(NDIS),ERR(NDIS)
C
c      IF(QWRITE) THEN
c       WRITE(6,'(1X,''calculate Lanczos eigenvalue errors...'')')
c       TZERO=TIMER(0D0)
c      ENDIF
C
C reset ERR
      DO 1 I=1,NDIS
       ERR(I)=0D0
1     CONTINUE
C
C set EPS,IKL,IT
C
C according to C&W's comments EPS is twice machine epsilon (see my
C comments to BISEC)
CMD
C NAG routine
C      EPS=2D0*X02AJF()
      EPS=2.D0*dlamch('e')
c for ibm 43p260
c      EPS=2.D0*2.2204460492503131d-16

C FPS FORTRAN runtime library routine
C     EPS=2D0*DFFRAC()
C
C IKL is the seed for the GENRAN routine; set to arbitrary value
      IKL=1
C
C IT is the maximum no. of inverse iteration steps for each eigenvalue;
C set to arbitrary small number
      IT=10
C
      ZERO=0D0
      ONE=1D0
      ITER=IT
      MP1=MEV+1
      MM1=MEV-1
      BETAM=BETA(MP1)
      BETA(MP1)=ZERO
      TSUM=ABS(ALPHA(1))
      DO 30 I=2,MEV
30    TSUM=TSUM+ABS(ALPHA(I))+BETA(I)
      EPS3=EPS*TSUM
      EPS4=DBLE(MEV)*EPS3
      ILL=IKL
      CALL GENRAN(ILL,G,MEV)
      GSUM=ZERO
      DO 40 I=1,MEV
40    GSUM=GSUM+ABS(G(I))
      GSUM=EPS4/GSUM
      DO 50 I=1,MEV
50    G(I)=GSUM*G(I)
      DO 180 JEV=1,NDIS
      IF(MP(JEV).NE.1) GO TO 180
      IT=1
      X1=VS(JEV)
      DO 60 I=1,MEV
60    V2(I)=G(I)
70    CONTINUE
      U=ALPHA(1)-X1
      Z=BETA(2)
      DO 90 I=2,MEV
      IF(BETA(I).GT.ABS(U)) GO TO 80
      V1(I-1)=Z/U
      V2(I-1)=V2(I-1)/U
      V2(I)=V2(I)-BETA(I)*V2(I-1)
      RATIO=BETA(I)/U
      U=ALPHA(I)-X1-Z*RATIO
      Z=BETA(I+1)
      GO TO 90
80    CONTINUE
      RATIO=U/BETA(I)
      BETA(I)=-BETA(I)
      V1(I-1)=ALPHA(I)-X1
      U=Z-RATIO*V1(I-1)
      Z=-RATIO*BETA(I+1)
      TEMP=V2(I-1)
      V2(I-1)=V2(I)
      V2(I)=TEMP-RATIO*V2(I)
90    CONTINUE
      IF(U.EQ.ZERO) U=EPS3
      V2(MEV)=V2(MEV)/U
      DO 110 II=1,MM1
      I=MEV-II
      IF(BETA(I+1).LT.ZERO) GO TO 100
      V2(I)=V2(I)-V1(I)*V2(I+1)
      GO TO 110
100   BETA(I+1)=-BETA(I+1)
      V2(I)=(V2(I)-V1(I)*V2(I+1)-BETA(I+2)*V2(I+2))/BETA(I+1)
110   CONTINUE
      XNORM=ABS(V2(MEV))
      DO 120 II=1,MM1
      I=MEV-II
120   XNORM=XNORM+ABS(V2(I))
      IF(XNORM.GE.ONE) GO TO 140
      IT=IT+1
      IF(IT.GT.ITER) GO TO 140
      XU=EPS4/XNORM
      DO 130 I=1,MEV
130   V2(I)=V2(I)*XU
      GO TO 70
140   CONTINUE
      SUM=SCALAR(V2(1),V2(1),MEV)
      SUM=ONE/SQRT(SUM)
      DO 150 II=1,MEV
150   V2(II)=SUM*V2(II)
      EST=BETAM*ABS(V2(MEV))
      IF(IT.GT.ITER) EST=-EST
      ERR(JEV)=EST
180   CONTINUE
      BETA(MP1)=BETAM
C
c      IF(QWRITE) THEN
c       TIME=TIMER(TZERO)
c       WRITE(6,'(1X,''time taken = '',F14.6,'' seconds'',/)')TIME
c      ENDIF
C
      RETURN
      END
