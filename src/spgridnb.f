C     *******************************************************************
C
C       SPGRID_NB - NoBoundary sparse grid point coordinates.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPGRID_NB(LEVELSEQ, NLEVELS, D, X, TOTALPOINTS)
C     *******************************************************************
C
C       SPGRID_NB computes NoBoundary sparse grid points on [0,1]^D.
C       All nodes are interior midpoints; no boundary nodes included.
C
C       1-D node formula (all levels):
C            lev >= 0 ->  {(2i-1)/2^(lev+1), i=1,...,2^lev}
C
C       Parameters:
C
C         Input,  INTEGER LEVELSEQ(NLEVELS,D)    - multi-index array
C         Input,  INTEGER NLEVELS                - rows in LEVELSEQ
C         Input,  INTEGER D                      - dimension
C         Output, DOUBLE PRECISION X(TOTALPTS,D) - grid coordinates
C         Input,  INTEGER TOTALPOINTS            - total grid points
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER NLEVELS, D, TOTALPOINTS
      INTEGER LEVELSEQ(NLEVELS, D)
      DOUBLE PRECISION X(TOTALPOINTS, D)

      INTEGER KL, K, I, J, LEV, IDX, NPTS_KL
      INTEGER NPTS_DIM(50), REP(50)

      IDX = 1

      DO KL = 1, NLEVELS

            NPTS_KL = 1
            DO K = 1, D
                  NPTS_DIM(K) = 2**LEVELSEQ(KL, K)
                  NPTS_KL = NPTS_KL * NPTS_DIM(K)
            END DO

            REP(1) = 1
            DO K = 2, D
                  REP(K) = REP(K-1) * NPTS_DIM(K-1)
            END DO

            DO I = 0, NPTS_KL - 1
                  DO K = 1, D
                        LEV = LEVELSEQ(KL, K)
                        IF (D .EQ. 1) THEN
                              J = MOD(I, NPTS_DIM(K))
                        ELSE
                              J = MOD(I / REP(K), NPTS_DIM(K))
                        END IF
                        X(IDX+I, K) = DBLE(2*J+1) /
     &                                DBLE(2**(LEV+1))
                  END DO
            END DO

            IDX = IDX + NPTS_KL

      END DO

      RETURN
      END
