C-----------------------------------------------------------------------
C
C#######################################################################
      SUBROUTINE BISEC(ALPHA,BETA,BETA2,MEV,VS,MP,VB,NMAX,BB,UB,NDIS)
C#######################################################################
C
C calculate eigenvalues of a tridiagonal Lanczos matrix by the method
C of bisection and detect those that are 'spurious'
C
C this routine is copied from C&W with the following changes:
C omission of all C&W's comment lines and WRITE statements (and the
C  IWRITE variable)
C addition of my comment lines
C minor formatting and FORTRAN style changes
C IMPLICIT typing (the variable LB is renamed BB)
C restriction to just one interval, removing the arrays LBD,UBD,IDEF,
C  the 430 loop and the variables JIND,NINT,ISKIP
C insertion of STOP to replace the C&W WRITEs
C replacement of some GO TO's with IF blocks
C replacement of MP1 variable with MEV+1 everywhere
C explicit dimensions for all arrays
C reordering of call parameters
C addition of NMAX call parameter (half the MMAX of C&W's comments)
C replacement of call parameters EPS, TTOL by internally-determined
C  values
C replacement of IC call parameter by an internal setting for MXSTUR,
C  and elimination of related ICT and IC0 counters
C addition of explanatory writing/timing and the QWRITE parameter
C
C meaning of call parameters:
C ALPHA: diagonal matrix elements, unchanged on exit
C BETA: subdiagonal matrix elements, unchanged on exit
C BETA2: computed squared subdiagonal matrix elements
C MEV: order of matrix, unchanged on exit
C VS: computed distinct eigenvalues in increasing order
C MP: corresponding multiplicities (except that 0 means spurious)
C VB: workspace
C NMAX: max. no. of eigenvalues to be computed, unchanged on exit
C BB: (open) lower bound on eigenvalues to be found, unchanged on exit
C UB: (closed) upper bound on eigenvalues to be found, unchanged on exit
C NDIS: no. of computed distinct eigenvalues
C QWRITE: if true, do timing and write out
C
C note BETA(1) is assumed to be zero
C
      IMPLICIT DOUBLE PRECISION (A-H,O-P,R-Z),LOGICAL (Q)
      DIMENSION ALPHA(MEV),BETA(MEV+1),BETA2(MEV+1),VS(NMAX),MP(NMAX),
     1          VB(NMAX*2)
C
C I could insert a STOP to check NDIS<NMAX, but I have not
C
c      QWRITE = .FALSE.
c      IF(QWRITE) THEN
c       WRITE(6,'(1X,''calculate eigenvalues of Lanczos matrix...'')')
c       TZERO=TIMER(0D0)
c      ENDIF
C
C according to C&W's comments:
C
C (i) EPS is twice machine epsilon MACHEP; I assume that what they mean
C by MACHEP is the smallest positive real number x for which 1+x>1 (x
C is, confusingly, often called EPS)
CMD
C NAG routine
C     EPS=2D0*X02AJF()
      EPS=2.D0*dlamch('e')
c for ibm 43p260
c      EPS=2.D0*2.2204460492503131d-16
c      EPS=2.2204460492503131d-16

C FPS FORTRAN runtime library routine
C     EPS=2D0*DFFRAC()
C
C (ii) TTOL=EPS*TKMAX, with
      TKMAX=0.D0
      DO 1 I=1,MEV
       TKMAX=MAX(ABS(ALPHA(I)),BETA(I),TKMAX)
1     CONTINUE
      TTOL=EPS*TKMAX
C
      ZERO=0.D0
      ONE=1.D0
      HALF=0.5D0
C should be large enough even if ALL the evalues are calculated; I used
C to have 30*MEV, but quite often ran into the STOP below
c      MXSTUR=100*MEV
      MXSTUR=30*MEV
      NDIS=0
      BETAM=BETA(MEV+1)
      BETA(MEV+1)=ZERO
      DO 10 I=1,MEV+1
10    BETA2(I)=BETA(I)**2
      TEMP=DBLE(MEV+1000)
      EP0=TEMP*TTOL
      EP1=SQRT(TEMP)*TTOL
      NA=0
      MD=0
      NG=0
      X1=UB
      ISTURM=1
      GO TO 330
50    NA=NEV
      X1=BB
      ISTURM=2
      GO TO 330
60    CONTINUE
      MT=NEV
      IEST=30*MT
      IF(IEST.GE.MXSTUR)
     1 STOP'estimated number of Sturms exceeds limit in BISEC'
      IF(MT.EQ.0) GOTO 410
      DO 140 I=1,MT
       VB(I)=BB
       MTI=MT+I
       VB(MTI)=UB
140   CONTINUE
      K=MT
150   CONTINUE
      XL=VB(K)
      MTK=MT+K
      XU=VB(MTK)
      ISTURM=3
      X1=XU
      GO TO 330
160   NU=NEV
      ISTURM=4
170   CONTINUE
      X1=(XL+XU)*HALF
      XS=ABS(XL)+ABS(XU)
      X0=XU-XL
      EPT=EPS*XS+EP1
      IF(X0.LE.EPT) GO TO 230
      GO TO 330
180   CONTINUE
      IF(NEV.LT.K) GO TO 190
      XL=X1
      GO TO 170
190   CONTINUE
      XU=X1
      NU=NEV
      IF(NEV.EQ.0) GO TO 210
      DO 200 I=1,NEV
200   VB(I)=MAX(X1,VB(I))
210   NEV1=NEV+1
      DO 220 II=NEV1,K
       I=MT+II
       VB(I)=MIN(X1,VB(I))
220   CONTINUE
      GO TO 170
230   CONTINUE
      NDIS=NDIS+1
      MD=MD+1
      VS(NDIS)=X1
      JSTURM=1
      X1=XL-EP0
      GO TO 370
240   KL=KEV
      JL=JEV
      JSTURM=2
      X1=XU+EP0
      GO TO 370
250   JU=JEV
      KU=KEV
      IF(KL-KU-1.EQ.0) GO TO 290
      IF(KU.EQ.NU) GO TO 280
260   CONTINUE
      ISTURM=5
      X1=X1+EP0
      GO TO 330
270   KNE=KU-NEV
      KU=NEV
      IF(KNE.NE.0) GO TO 260
280   MPEV=KL-KU
      KNEW=KU
      GO TO 300
290   CONTINUE
      MPEV=1
      IF(JU.LT.JL) MPEV=0
      KNEW=K-1
300   K=KNEW
      MP(NDIS)=MPEV
      IF(MPEV.GE.1) NG=NG+1
      IF(K.LE.0) GO TO 410
      DO 320 I=1,KNEW
320   VB(I)=MAX(X1,VB(I))
      GO TO 150
330   NEV=-NA
      YU=ONE
      DO 360 I=1,MEV
       IF(YU.EQ.ZERO) THEN
        YV=BETA(I)/EPS
       ELSE
        YV=BETA2(I)/YU
       ENDIF
       YU=X1-ALPHA(I)-YV
       IF(YU.LT.ZERO) NEV=NEV+1
360   CONTINUE
      GO TO (50,60,160,180,270), ISTURM
370   KEV=-NA
      YU=ONE
      DO 400 I=MEV,1,-1
       IF(YU.EQ.ZERO) THEN
        YV=BETA(I+1)/EPS
       ELSE
        YV=BETA2(I+1)/YU
       ENDIF
       YU=X1-ALPHA(I)-YV
       JEV=0
       IF(YU.LT.ZERO) THEN
        KEV=KEV+1
        JEV=1
       ENDIF
400   CONTINUE
      JEV=KEV-JEV
      GO TO (240,250), JSTURM
410   CONTINUE
      BETA(MEV+1)=BETAM
C
c      IF(QWRITE) THEN
c       TIME=TIMER(TZERO)
c       WRITE(6,'(1X,''time taken = '',F14.6,'' seconds'',/)')TIME
c      ENDIF
C
      RETURN
      END
