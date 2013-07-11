C-----------------------------------------------------------------------
C
C#######################################################################
      FUNCTION SCALAR(A,B,N)
C#######################################################################
C
Copyright Matthew J. Bramley 1992
C
C calculate scalar product A.B
C
      IMPLICIT DOUBLE PRECISION (A-H,O-P,R-Z),LOGICAL (Q)
      DIMENSION A(N),B(N)
C
      SUM=0D0
      DO 1 I=1,N
       SUM=SUM+A(I)*B(I)
1     CONTINUE
      SCALAR=SUM
C
      RETURN
      END
