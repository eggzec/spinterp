C     *******************************************************************
C
C       SPGRID_GP - Gauss-Patterson sparse grid point coordinates.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPGRID_GP(LEVELSEQ, NLEVELS, D, X, TOTALPOINTS)
C     *******************************************************************
C
C       SPGRID_GP computes Gauss-Patterson sparse grid points on [0,1]^D.
C       Uses the nested Gauss-Patterson nodes (max level 6).
C
C       Node count per 1-D level:
C            lev = 0  ->  1  (midpoint)
C            lev >= 1 ->  2^lev  (new GP nodes at this level)
C       Full node set at level lev has 2^(lev+1)-1 nodes.
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
      INTEGER DIM_ACTIVE(50), NDIMS
      INTEGER REP(50), NPTS_DIM(50)
      INTEGER NX
      DOUBLE PRECISION ABSC(127)

      IDX = 1

      DO K = 1, D
            DO I = 1, TOTALPOINTS
                  X(I, K) = 0.5D+00
            END DO
      END DO

      DO KL = 1, NLEVELS

            NDIMS = 0
            NPTS_KL = 1
            DO K = 1, D
                  LEV = LEVELSEQ(KL, K)
                  IF (LEV .EQ. 0) THEN
                        NPTS_DIM(K) = 1
                  ELSE
                        NPTS_DIM(K) = 2**LEV
                        NDIMS = NDIMS + 1
                        DIM_ACTIVE(NDIMS) = K
                  END IF
                  NPTS_KL = NPTS_KL * NPTS_DIM(K)
            END DO

            REP(1) = 1
            DO K = 2, NDIMS
                  REP(K) = REP(K-1) * NPTS_DIM(DIM_ACTIVE(K-1))
            END DO

            DO I = 0, NPTS_KL - 1
                  DO K = 1, NDIMS
                        LEV = LEVELSEQ(KL, DIM_ACTIVE(K))
                        NX = 2**(LEV+1) - 1
                        CALL GP_ABSC(LEV, ABSC, NX)
                        IF (NDIMS .EQ. 1) THEN
                              J = MOD(I, NPTS_DIM(DIM_ACTIVE(K)))
                        ELSE
                              J = MOD(I/REP(K),
     &                              NPTS_DIM(DIM_ACTIVE(K)))
                        END IF
C                       New GP nodes occupy the even indices (0,2,4,...)
C                       of the full node set. Skip by 2 to get new nodes.
                        X(IDX+I, DIM_ACTIVE(K)) = ABSC(J*2 + 1)
                  END DO
            END DO

            IDX = IDX + NPTS_KL

      END DO

      RETURN
      END
