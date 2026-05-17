C     *******************************************************************
C
C       SPQUADW_CC - Quadrature weights, Clenshaw-Curtis / Chebyshev grid.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPQUADW_CC(LEVELSEQ, NLEVELS, D, W, NW)
C     *******************************************************************
C
C       SPQUADW_CC returns the quadrature weights for the CC / Chebyshev
C       sparse grid.
C
C       The weight for each node in a subgrid is:
C            w = 1 / (product over dims of nw_k)
C       where:
C            lev=0: nw=1, np=1
C            lev=1,2: nw=4, np=2
C            lev>=3: nw=2^lev, np=2^(lev-1)
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

      INTEGER KL, K, L, LVAL, NP, INDEX
      DOUBLE PRECISION WVAL

      INDEX = 1

      DO KL = 1, NLEVELS
            NP = 1
            WVAL = 1.0D+00
            DO L = 1, D
                  LVAL = LEVELSEQ(KL, L)
                  IF (LVAL .EQ. 0) THEN
                        NP = NP * 1
                        WVAL = WVAL * 1.0D+00
                  ELSE IF (LVAL .LE. 2) THEN
                        NP = NP * 2
                        WVAL = WVAL * 4.0D+00
                  ELSE
                        NP = NP * 2**(LVAL-1)
                        WVAL = WVAL * DBLE(2**LVAL)
                  END IF
            END DO
            DO K = INDEX, INDEX + NP - 1
                  W(K) = 1.0D+00 / WVAL
            END DO
            INDEX = INDEX + NP
      END DO

      RETURN
      END
