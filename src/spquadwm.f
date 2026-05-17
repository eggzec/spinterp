C     *******************************************************************
C
C       SPQUADW_M - Quadrature weights, Maximum-norm grid.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPQUADW_M(LEVELSEQ, NLEVELS, D, W, NW)
C     *******************************************************************
C
C       SPQUADW_M returns the quadrature weights for the Maximum-norm
C       sparse grid.
C
C       1-D weight vectors:
C            lev=0: [1/4, 1/2, 1/4]  (3 nodes: 0, 0.5, 1)
C            lev>=1: 1/2^(lev+1) for all 2^lev nodes
C
C       The D-dimensional weight is the product of 1-D weights.
C
C       Parameters:
C
C         Input,  INTEGER LEVELSEQ(NLEVELS,D) - multi-index array
C         Input,  INTEGER NLEVELS             - rows in LEVELSEQ
C         Input,  INTEGER D                   - dimension
C         Output, DOUBLE PRECISION W(NW)      - quadrature weights
C         Input,  INTEGER NW                  - length of W
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER NLEVELS, D, NW
      INTEGER LEVELSEQ(NLEVELS, D)
      DOUBLE PRECISION W(NW)

      INTEGER KL, K, L, LVAL, NPTS, INDEX
      INTEGER NPTS_DIM(50), REP(50)
      DOUBLE PRECISION W1D(3, 50), WVAL
      INTEGER J

      INDEX = 1

      DO KL = 1, NLEVELS
            NPTS = 1
            DO L = 1, D
                  LVAL = LEVELSEQ(KL, L)
                  IF (LVAL .EQ. 0) THEN
                        NPTS_DIM(L) = 3
                        W1D(1, L) = 0.25D+00
                        W1D(2, L) = 0.50D+00
                        W1D(3, L) = 0.25D+00
                  ELSE
                        NPTS_DIM(L) = 2**LVAL
                        DO J = 1, NPTS_DIM(L)
                              W1D(J, L) = 1.0D+00/DBLE(2**(LVAL+1))
                        END DO
                  END IF
                  NPTS = NPTS * NPTS_DIM(L)
            END DO

            REP(1) = 1
            DO L = 2, D
                  REP(L) = REP(L-1) * NPTS_DIM(L-1)
            END DO

            DO K = 0, NPTS-1
                  WVAL = 1.0D+00
                  DO L = 1, D
                        IF (D .EQ. 1) THEN
                              J = MOD(K, NPTS_DIM(L)) + 1
                        ELSE
                              J = MOD(K/REP(L), NPTS_DIM(L)) + 1
                        END IF
                        WVAL = WVAL * W1D(J, L)
                  END DO
                  W(INDEX + K) = WVAL
            END DO

            INDEX = INDEX + NPTS
      END DO

      RETURN
      END
