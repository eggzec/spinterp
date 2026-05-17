C     *******************************************************************
C
C       SPGRID_CB - Chebyshev sparse grid point coordinates.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPGRID_CB(LEVELSEQ, NLEVELS, D, X, TOTALPOINTS)
C     *******************************************************************
C
C       SPGRID_CB computes Chebyshev sparse grid points on [0,1]^D.
C       Nodes are the extrema of Chebyshev polynomials (Chebyshev-Lobatto).
C
C       1-D node formula (same point count as CC):
C            lev = 0        ->  0.5
C            lev = 1        ->  {0, 1}
C            lev >= 2       ->  {0.5 - cos((2i-1)*pi/2^lev)/2,
C                                i=1,...,2^(lev-1)}
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

      DOUBLE PRECISION PI
      PARAMETER (PI = 3.14159265358979323846D+00)

      INTEGER NLEVELS, D, TOTALPOINTS
      INTEGER LEVELSEQ(NLEVELS, D)
      DOUBLE PRECISION X(TOTALPOINTS, D)

      INTEGER KL, K, I, J, LEV, IDX, NPTS_KL
      INTEGER DIM_ACTIVE(50), NDIMS
      INTEGER REP(50), NPTS_DIM(50)

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
                  ELSE IF (LEV .LE. 2) THEN
                        NPTS_DIM(K) = 2
                        NDIMS = NDIMS + 1
                        DIM_ACTIVE(NDIMS) = K
                  ELSE
                        NPTS_DIM(K) = 2**(LEV-1)
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
                        IF (NDIMS .EQ. 1) THEN
                              J = MOD(I, NPTS_DIM(DIM_ACTIVE(K)))
                        ELSE
                              J = MOD(I/REP(K),
     &                              NPTS_DIM(DIM_ACTIVE(K)))
                        END IF

                        IF (LEV .EQ. 1) THEN
                              IF (J .EQ. 0) THEN
                                    X(IDX+I, DIM_ACTIVE(K)) = 0.0D+00
                              ELSE
                                    X(IDX+I, DIM_ACTIVE(K)) = 1.0D+00
                              END IF
                        ELSE
C                             Chebyshev-Lobatto new nodes at level lev:
C                             0.5 - cos((2j+1)*pi/2^lev) / 2
                              X(IDX+I, DIM_ACTIVE(K)) = 0.5D+00
     &                              - DCOS(DBLE(2*J+1)*PI
     &                              / DBLE(2**LEV)) / 2.0D+00
                        END IF
                  END DO
            END DO

            IDX = IDX + NPTS_KL

      END DO

      RETURN
      END
