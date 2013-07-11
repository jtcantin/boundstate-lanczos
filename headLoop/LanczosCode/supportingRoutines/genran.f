C-----------------------------------------------------------------------
C
C#######################################################################
      SUBROUTINE GENRAN(ISEED,R,N)
C#######################################################################
C
Copyright Matthew J. Bramley 1992
C
C generate a sequence of N random numbers between 0 and 1 in R from the
C seed ISEED
C
C beware the comments in NR on the poor quality of library random
C number generators
C
      IMPLICIT DOUBLE PRECISION (A-H,O-P,R-Z),LOGICAL (Q)
      DIMENSION R(N)
C
CMD
C library routines SRAND and RAND on SG
      CALL SRAND(ISEED)
      DO 1 I=1,N
       R(I)=RAND()
1     CONTINUE
C
C library routine DRAND on FPS
C     X=DRAND(ISEED)
C     DO 2 I=1,N
C      R(I)=DRAND(0)
2     CONTINUE
C
      RETURN
      END
