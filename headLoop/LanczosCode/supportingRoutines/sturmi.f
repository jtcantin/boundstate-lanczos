C-----------------------------------------------------------------------
C
C#######################################################################
      SUBROUTINE STURMI(ALPHA,BETA,BETA2,MMAX,X1,MK1,MK2,QWRITE)
C#######################################################################
C
C calculate the sizes of Lanczos matrix for which the specified
C eigenvalue (i) first converges and (ii) doubles
C
C this routine is copied from C&W with the following changes:
C omission of all C&W's comment lines and WRITE statements (and the
C  IWRITE variable)
C addition of my comment lines
C minor formatting and FORTRAN style changes
C IMPLICIT typing
C replacement of TOLN call parameter by internally-determined value
C addition of the BETA2 array as a call parameter
C explicit dimensions for all arrays
C reordering of call parameters
C replacement of call parameter EPSM by internally-determined value
C elimination of IC call parameter and replacement by STOP
C addition of explanatory writing/timing and the QWRITE parameter
C
C note the internal determination of the tolerances each time this
C routine is called might be wasteful
C
C meaning of call parameters:
C ALPHA: diagonal matrix elements, unchanged on exit
C BETA: subdiagonal matrix elements, unchanged on exit
C BETA2: squared subdiagonal matrix elements, unchanged on exit
C MMAX: maximum order of matrix, unchanged on exit
C X1: eigenvalue, unchanged on exit
C MK1: size for 1 eigenvalue
C MK2: size for 2 eigenvalues or MMAX, whichever is smaller
C QWRITE: if true, do timing and write out
C
      IMPLICIT DOUBLE PRECISION (A-H,O-P,R-Z),LOGICAL (Q)
      DIMENSION ALPHA(MMAX),BETA(MMAX),BETA2(MMAX)
C
      IF(QWRITE) THEN
       WRITE(6,'(1X,
     1 ''calculate Lanczos matrix sizes for eigenvectors...'')')
C       TZERO=TIMER(0D0)
      ENDIF
C
C EPSM is twice machine epsilon, as specified C&W vol.II p.22
CMD
C NAG routine
c      EPSM=2D0*X02AJF()
      EPSM=2.D0*dlamch('e')

C
C TTOL defined C&W vol.II p.22
      TKMAX=0D0
      DO 1 I=1,MMAX
       TKMAX=MAX(ABS(ALPHA(I)),BETA(I),TKMAX)
1     CONTINUE
      TTOL=EPSM*TKMAX
C MULTOL defined C&W vol.II p.22, renamed XMLTOL
      XMLTOL=DBLE(MMAX+1000)*TTOL
C SCALE0 assigned C&W vol.II p.49
      SCALE0=5D0
C STUTOL assigned C&W vol.II p.54
      STUTOL=SCALE0*XMLTOL
C RELTOL specified C&W vol.II p.22
      RELTOL=1D-10
      TEMP=ABS(X1)*RELTOL
      TOLN=MAX(TEMP,STUTOL)
C
      MK1=0
      MK2=0
      ZERO=0D0
      ONE=1D0
      BETA(1)=ZERO
      EVL=X1-TOLN
      EVU=X1+TOLN
      U1=ONE
      U2=ONE
      IC0=0
      IC1=0
      IC2=0
      DO 60 I=1,MMAX
      IF(U1.NE.ZERO) GOTO 10
      V1=BETA(I)/EPSM
      GOTO 20
10    V1=BETA2(I)/U1
20    U1=EVL-ALPHA(I)-V1
      IF(U1.LT.ZERO) IC1=IC1+1
      IF(U2.NE.ZERO) GOTO 30
      V2=BETA(I)/EPSM
      GOTO 40
30    V2=BETA2(I)/U2
40    U2=EVU-ALPHA(I)-V2
      IF(U2.LT.ZERO) IC2=IC2+1
      ICD=IC1-IC2
      IC=ICD-IC0
      IF(IC.GE.1) GOTO 50
      GOTO 60
50    CONTINUE
      IF(IC0.EQ.0) MK1=I
      IC0=IC0+1
      IF(IC0.GT.1) GOTO 70
60    CONTINUE
      I=I-1
      IF(IC0.EQ.0) MK1=MMAX
70    MK2=I
      IC=ICD
      IF(IC.EQ.0)STOP'eigenvalue not found in STURMI'
C
      IF(QWRITE) THEN
C       TIME=TIMER(TZERO)
       WRITE(6,'(1X,''time taken = '',F14.6,'' seconds'',/)')TIME
      ENDIF
C
      RETURN
      END
