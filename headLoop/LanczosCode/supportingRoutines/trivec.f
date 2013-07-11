C-----------------------------------------------------------------------
C
C#######################################################################
      SUBROUTINE TRIVEC(ALPHA,BETA,BET2,WRK1,WRK2,NIT,EVAL,NGOOD,EVTR,
     1                  MAMIN)
C#######################################################################
C
C calculate many eigenvectors of a tridiagonal Lanczos matrix for sizes
C of the latter appropriate for later calculating Ritz vectors
C
      IMPLICIT DOUBLE PRECISION (A-H,O-P,R-Z),LOGICAL (Q)
      DIMENSION ALPHA(NIT),BETA(NIT+1),BET2(NIT+1),WRK1(NIT),WRK2(NIT),
     1 EVAL(NGOOD),EVTR(NIT,NGOOD),MAMIN(NGOOD)
C
      WRITE(6,'(1X,
     1 ''calculate eigenvectors of tridiagonal matrix...'')')
C      TZERO=TIMER(0D0)
C
C obtain random vector WRK2 for INVERM (see also comments in INVERR):
C set seed to arbitrary value of 1
      CALL GENRAN(1,WRK2,NIT)
C
C loop good eigenvalues
      DO 1 IGOOD=1,NGOOD
C
C calculate number of Lanczos iterations needed (i) to obtain the
C eigenvalue (even if it is unconverged, since the unconverged value
C will be searched for) MK1 and (ii) to obtain the eigenvalue doubled
C MK2
       CALL STURMI(ALPHA,BETA,BET2,NIT,EVAL(IGOOD),MK1,MK2,.FALSE.)
C
C choose endpoints MA,ML and increment IDELTA for searching between MK1
C and MK2 for an appropriate size of Lanczos matrix for use in
C calculating the corresponding eigenvector
       IF(MK2.EQ.NIT) THEN
C change made from C&W: MA must not be larger than NIT
        MA=MIN((5*MK1)/4+1,NIT)
        ML=NIT
C change made from C&W in the line below because NIT, not (11*NIT)/8+12,
C is the largest matrix size; also 10 search intervals replaced by 5
        IDELTA=(MIN(NIT,(13*MK1)/8+1)-MA)/10+1
C        IDELTA=(MIN(NIT,(13*MK1)/8+1)-MA)/5+1
       ELSE
        MA=(3*MK1+MK2)/4+1
        ML=(MK1+3*MK2)/4
        IDELTA=((3*MK1+5*MK2)/8+1-MA)/10+1
C        IDELTA=((3*MK1+5*MK2)/8+1-MA)/5+1
       ENDIF
C
C loop on different Lanczos matrix orders to find an order for which the
C eigenvector error is minimal
       ICOUNT=0
       QFIND=.FALSE.
C initial large value for minimum error
       ERRMIN=1.D0
2      CONTINUE
       ICOUNT=ICOUNT+1
C
C obtain corresponding eigenvector of Lanczos matrix
       CALL INVERM(ALPHA,BETA,WRK1,EVTR(1,IGOOD),WRK2,MA,EVAL(IGOOD),
     1             ERROR,QCONV,.FALSE.)
C
C QFIND becomes true if any of the INVERM calls converges
       QFIND=QFIND.OR.QCONV
C
       IF(ERROR.LE.ERRMIN) THEN
        ERRMIN=ERROR
        MAMIN(IGOOD)=MA
       ENDIF
C
       MA=MA+IDELTA
C 10 search intervals replaced by 5
C       IF((MA.LT.ML).AND.(ICOUNT.LE.5)) GOTO 2
       IF((MA.LT.ML).AND.(ICOUNT.LE.10)) GOTO 2
C
       IF(.NOT.QFIND) THEN
        WRITE(6,'(1X,''IGOOD='',I5)')IGOOD
        STOP'INVERM never converged in TRIVEC'
       ENDIF
C
C call INVERM again for the optimum order
       CALL INVERM(ALPHA,BETA,WRK1,EVTR(1,IGOOD),WRK2,MAMIN(IGOOD),
     1             EVAL(IGOOD),ERROR,QCONV,.FALSE.)
C
       IF(ERRMIN.GT.1D-10)
     1 WRITE(6,'(1X,''vector no. '',I4,'' has ERRMIN value > 1D-10: '',
     2 D8.2,/)')IGOOD,ERRMIN
C
1     CONTINUE
C
C      TIME=TIMER(TZERO)
c      WRITE(6,'(1X,''time taken = '',F14.6,'' seconds'',/)')TIME
C
      RETURN
      END
