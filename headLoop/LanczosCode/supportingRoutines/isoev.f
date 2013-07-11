C-----------------------------------------------------------------------
C
C#######################################################################
      SUBROUTINE ISOEV(VS,MP,NDIS)
C#######################################################################
C
C flag 'good' eigenvalues of a tridiagonal Lanczos matrix that are not
C isolated (i.e. ones with a 'spurious' eigenvalue as closest neighbour)
C
C this routine is a rewritten and simplified version of that of C&W;
C note that the middle 1 in the case 0--too close--1--even closer--1 is
C said to be isolated (see C&W's comments to INVERR)
C
C meaning of call parameters:
C VS: computed distinct eigenvalues in increasing order
C MP: corresponding multiplicities (except that 0 means spurious); this
C  routine replaces multiplicities 1 by -1 if the corresponding
C  eigenvalue is not isolated
C NDIS: no. of computed distinct eigenvalues
C
      IMPLICIT DOUBLE PRECISION (A-H,O-P,R-Z),LOGICAL (Q)
      DIMENSION VS(NDIS),MP(NDIS)
C
C see comments on p.21-22 of C&W vol.II
      GAPTOL=1.D-8
C
      DO 1 J=1,NDIS
       IF(MP(J).EQ.1) THEN
        IF(J.NE.NDIS) T1=VS(J+1)-VS(J)
        IF(J.NE.1) T0=VS(J)-VS(J-1)
        IF((J.NE.1).AND.(J.NE.NDIS)) THEN
         IF(T0.LT.T1) THEN
          G=-T0
         ELSE
          G=T1
         ENDIF
        ELSE
         IF(J.EQ.1) G=VS(2)-VS(1)
         IF(J.EQ.NDIS) G=VS(NDIS-1)-VS(NDIS)
        ENDIF
        IF(G.GE.0D0) THEN
         I=J+1
        ELSE
         I=J-1
        ENDIF
        IF(MP(I).EQ.0) THEN
         GAP=ABS(G)
         TEMP=MAX(GAPTOL,GAPTOL*ABS(VS(J)))
         IF(GAP.LE.TEMP) MP(J)=-MP(J)
        ENDIF
       ENDIF
1     CONTINUE
C
      RETURN
      END











