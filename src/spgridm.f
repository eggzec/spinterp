C     *******************************************************************
C
C       SPGRID_M - Maximum-norm sparse grid point coordinates.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPGRID_M(LEVELSEQ, NLEVELS, D, X, TOTALPOINTS)
C     *******************************************************************
C
C       SPGRID_M computes Maximum-norm sparse grid points on [0,1]^D.
C
C       1-D node formula:
C            lev = 0  ->  {0, 0.5, 1}  (3 nodes)
C            lev >= 1 ->  {(2i-1)/2^(lev+1), i=1,...,2^lev}
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

      INTEGER KL, K, I, J, LEV, IDX, NPTS_KL, NPTS_DIM(50)
      INTEGER REP(50)
      DOUBLE PRECISION COORD

      IDX = 1

      DO KL = 1, NLEVELS

            NPTS_KL = 1
            DO K = 1, D
                  LEV = LEVELSEQ(KL, K)
                  IF (LEV .EQ. 0) THEN
                        NPTS_DIM(K) = 3
                  ELSE
                        NPTS_DIM(K) = 2**LEV
                  END IF
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

                        IF (LEV .EQ. 0) THEN
                              IF (J .EQ. 0) THEN
                                    COORD = 0.0D+00
                              ELSE IF (J .EQ. 1) THEN
                                    COORD = 0.5D+00
                              ELSE
                                    COORD = 1.0D+00
                              END IF
                        ELSE
                              COORD = DBLE(2*J+1) /
     &                                DBLE(2**(LEV+1))
                        END IF
                        X(IDX+I, K) = COORD
                  END DO
            END DO

            IDX = IDX + NPTS_KL

      END DO

      RETURN
      END
