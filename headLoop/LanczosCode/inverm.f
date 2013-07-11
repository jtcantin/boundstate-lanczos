C-----------------------------------------------------------------------
C
C#######################################################################
      SUBROUTINE INVERM(ALPHA,BETA,V1,V2,G,MEV,X1,ERROR,QCONV,QWRITE)
C#######################################################################
C
C calculate the eigenvector of a tridiagonal Lanczos matrix
C corresponding to a given 'good' eigenvalue using inverse iteration
C
C this routine is copied from C&W with the following changes:
C omission of all C&W's comment lines and WRITE statements (and the
C  IWRITE variable)
C addition of my comment lines
C minor formatting and FORTRAN style changes
C IMPLICIT typing (the variable NORM is renamed XNORM and G becomes
C  DOUBLE PRECISION)
C explicit dimensions for all arrays
C reordering of call parameters
C replacement of call parameters EPS and IT by internal settings
C replacement of FINPRO by SCALAR for scalar product
C elimination of ERRORV parameter and its calculation
C addition of QCONV parameter
C addition of explanatory writing/timing and the QWRITE parameter
C
C meaning of call parameters:
C ALPHA: diagonal matrix elements, unchanged on exit
C BETA: subdiagonal matrix elements, unchanged on exit
C V1: work space
C V2: the calculated eigenvector (normalised)
C G: random vector, unchanged on exit
C MEV: order of matrix, unchanged on exit
C X1: the 'good' eigenvalue, unchanged on exit
C ERROR: error estimate for corresponding Ritz vector
C QCONV: false if iteration limit was reached, otherwise true
C QWRITE: if true, do timing and write out
C
      IMPLICIT DOUBLE PRECISION (A-H,O-P,R-Z),LOGICAL (Q)
      DIMENSION ALPHA(MEV),BETA(MEV+1),V1(MEV),V2(MEV),G(MEV)
C
c      IF(QWRITE) THEN
c       WRITE(6,'(1X,''calculate eigenvector of Lanczos matrix...'')')
c       TZERO=TIMER(0D0)
c      ENDIF
C
C EPS is twice machine epsilon, as specified C&W vol.II p.22
CMD
C NAG routine
c      EPS=2D0*X02AJF()
      EPS=2.D0*dlamch('e')
C IT is the maximum no. of inverse iteration steps; set to arbitrary
C small number
      IT=10
C
      ONE=1D0
      ZERO=0D0
      ITER=IT
      MP1=MEV+1
      MM1=MEV-1
      BETAM=BETA(MP1)
      BETA(MP1)=ZERO
      TSUM=ABS(ALPHA(1))
      DO 10 I=2,MEV
10    TSUM=TSUM+ABS(ALPHA(I))+BETA(I)
      EPS3=EPS*TSUM
      EPS4=DBLE(MEV)*EPS3
      GSUM=ZERO
      DO 20 I=1,MEV
20    GSUM=GSUM+ABS(G(I))
      GSUM=EPS4/GSUM
      DO 30 I=1,MEV
30    V2(I)=GSUM*G(I)
      IT=1
40    CONTINUE
      U=ALPHA(1)-X1
      Z=BETA(2)
      DO 60 I=2,MEV
      IF(BETA(I).GT.ABS(U)) GOTO 50
      V1(I-1)=Z/U
      V2(I-1)=V2(I-1)/U
      V2(I)=V2(I)-BETA(I)*V2(I-1)
      RATIO=BETA(I)/U
      U=ALPHA(I)-X1-Z*RATIO
      Z=BETA(I+1)
      GOTO 60
50    CONTINUE
      RATIO=U/BETA(I)
      BETA(I)=-BETA(I)
      V1(I-1)=ALPHA(I)-X1
      U=Z-RATIO*V1(I-1)
      Z=-RATIO*BETA(I+1)
      TEMP=V2(I-1)
      V2(I-1)=V2(I)
      V2(I)=TEMP-RATIO*V2(I)
60    CONTINUE
      IF(U.EQ.ZERO) U=EPS3
      V2(MEV)=V2(MEV)/U
      DO 80 II=1,MM1
      I=MEV-II
      IF(BETA(I+1).LT.ZERO) GOTO 70
      V2(I)=V2(I)-V1(I)*V2(I+1)
      GOTO 80
70    BETA(I+1)=-BETA(I+1)
      V2(I)=(V2(I)-V1(I)*V2(I+1)-BETA(I+2)*V2(I+2))/BETA(I+1)
80    CONTINUE
      XNORM=ABS(V2(MEV))
      DO 90 II=1,MM1
      I=MEV-II
90    XNORM=XNORM+ABS(V2(I))
      IF(XNORM.GE.ONE) GOTO 110
      IT=IT+1
      IF(IT.GT.ITER) GOTO 110
      XU=EPS4/XNORM
      DO 100 I=1,MEV
100   V2(I)=V2(I)*XU
      GOTO 40
110   CONTINUE
      SUM=SCALAR(V2(1),V2(1),MEV)
      SUM=ONE/SQRT(SUM)
      DO 120 II=1,MEV
120   V2(II)=SUM*V2(II)
      ERROR=ABS(V2(MEV))
      V1(MEV)=ALPHA(MEV)*V2(MEV)+BETA(MEV)*V2(MEV-1)-X1*V2(MEV)
      DO 130 J=2,MM1
      JM=MP1-J
      V1(JM)=ALPHA(JM)*V2(JM)+BETA(JM)*V2(JM-1)+BETA(JM+1)*V2(JM+1)
     1 -X1*V2(JM)
130   CONTINUE
      V1(1)=ALPHA(1)*V2(1)+BETA(2)*V2(2)-X1*V2(1)
      IF(IT.GT.ITER) THEN
       QCONV=.FALSE.
      ELSE
       QCONV=.TRUE.
      ENDIF
      BETA(MP1)=BETAM
C
c      IF(QWRITE) THEN
c       TIME=TIMER(TZERO)
c       WRITE(6,'(1X,''time taken = '',F14.6,'' seconds'',/)')TIME
c      ENDIF
C
      RETURN
      END
