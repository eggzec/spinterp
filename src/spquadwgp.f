C     *******************************************************************
C
C       SPQUADW_GP - Quadrature weights, Gauss-Patterson grid.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPQUADW_GP(LEVELSEQ, NLEVELS, D, W, NW,
     &                       W1D, NW1D, STARTID)
C     *******************************************************************
C
C       SPQUADW_GP returns the quadrature weights for the Gauss-Patterson
C       sparse grid.  Caller provides 1-D weights W1D and start indices
C       STARTID computed by GP_WEIGHTS.
C
C       Parameters:
C
C         Input,  INTEGER LEVELSEQ(NLEVELS,D) - multi-index array
C         Input,  INTEGER NLEVELS             - rows in LEVELSEQ
C         Input,  INTEGER D                   - dimension
C         Output, DOUBLE PRECISION W(NW)      - quadrature weights
C         Input,  INTEGER NW                  - length of W
C         Input,  DOUBLE PRECISION W1D(NW1D)  - 1-D weight table
C         Input,  INTEGER NW1D               - length of W1D
C         Input,  INTEGER STARTID(*)          - start index per level
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER NLEVELS, D, NW, NW1D
      INTEGER LEVELSEQ(NLEVELS, D)
      DOUBLE PRECISION W(NW)
      DOUBLE PRECISION W1D(NW1D)
      INTEGER STARTID(*)

      INTEGER KL, K, L, LVAL, NPTS, INDEX
      INTEGER NPTS_DIM(50), REP(50)
      INTEGER WIDSTART(50), WIDEND(50), WID(50), WID1
      DOUBLE PRECISION WVAL

      INDEX = 1

      DO KL = 1, NLEVELS
            NPTS = 1
            DO L = 1, D
                  LVAL = LEVELSEQ(KL, L)
                  NPTS_DIM(L) = 2**LVAL
                  NPTS = NPTS * NPTS_DIM(L)
                  WIDSTART(L) = STARTID(LVAL+1)
                  WID(L) = WIDSTART(L)
                  WIDEND(L) = WID(L) + NPTS_DIM(L) - 1
            END DO

            WID1 = WID(1)
            DO K = INDEX, INDEX + NPTS - 1
                  WVAL = W1D(WID1)
                  DO L = 2, D
                        WVAL = WVAL * W1D(WID(L))
                  END DO
                  W(K) = WVAL

                  IF (WID1 .LT. WIDEND(1)) THEN
                        WID1 = WID1 + 1
                  ELSE
                        WID1 = WIDSTART(1)
                        DO L = 2, D
                              IF (WID(L) .LT. WIDEND(L)) THEN
                                    WID(L) = WID(L) + 1
                                    GOTO 10
                              ELSE
                                    WID(L) = WIDSTART(L)
                              END IF
                        END DO
                  END IF
 10               CONTINUE
            END DO

            INDEX = INDEX + NPTS
      END DO

      RETURN
      END
