C     *******************************************************************
C
C       SPINTERP_CB - Polynomial interpolation, Chebyshev grid.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPINTERP_CB(D, Z, NZ, Y, NINTERP,
     &                        LEVELSEQ, NLEVELS, IP)
C     *******************************************************************
C
C       SPINTERP_CB evaluates the Chebyshev sparse grid interpolant at
C       NINTERP query points Y using hierarchical surpluses Z.
C
C       For each subgrid with >1 node, calls BARY_PD_STEP_CB with
C       dimensions sorted descending by level so the largest node count
C       is the innermost loop.
C
C       Parameters:
C
C         Input,  INTEGER D                     - dimension
C         Input,  DOUBLE PRECISION Z(NZ)        - hierarchical surpluses
C         Input,  INTEGER NZ                    - length of Z
C         Input,  DOUBLE PRECISION Y(NINTERP,D) - query points in [0,1]^D
C         Input,  INTEGER NINTERP               - number of query points
C         Input,  INTEGER LEVELSEQ(NLEVELS,D)   - multi-index set
C         Input,  INTEGER NLEVELS               - rows in LEVELSEQ
C         Output, DOUBLE PRECISION IP(NINTERP)  - interpolated values
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER D, NZ, NINTERP, NLEVELS
      DOUBLE PRECISION Z(NZ)
      DOUBLE PRECISION Y(NINTERP, D)
      INTEGER LEVELSEQ(NLEVELS, D)
      DOUBLE PRECISION IP(NINTERP)

      INTEGER KL, K, L, LVAL, NPTS, NDIMS, INDEX
      INTEGER ALLNX(50), DIMS(50), NORD(50), LEVEL(50)
      INTEGER ORDERARR(50), XTOT
      DOUBLE PRECISION XBUF(16384)
      DOUBLE PRECISION IPTEMP(65536)

      DO K = 1, NINTERP
            IP(K) = 0.0D+00
      END DO

      INDEX = 1

      DO KL = 1, NLEVELS

            NPTS = 1
            NDIMS = 0
            DO L = 1, D
                  LVAL = LEVELSEQ(KL, L)
                  LEVEL(L) = LVAL
                  IF (LVAL .EQ. 0) THEN
                  ELSE IF (LVAL .LE. 2) THEN
                        NPTS = NPTS * 2
                  ELSE
                        NPTS = NPTS * 2**(LVAL-1)
                  END IF
            END DO

            IF (NPTS .EQ. 1) THEN
C                 Only one (constant) node, no interpolation needed
                  DO K = 1, NINTERP
                        IP(K) = IP(K) + Z(INDEX)
                  END DO
                  INDEX = INDEX + 1

            ELSE
C                 Sort dimensions descending by level
                  DO L = 1, D
                        ORDERARR(L) = L
                  END DO
                  DO L = 1, D-1
                        DO K = L+1, D
                              IF (LEVEL(ORDERARR(K)) .GT.
     &                              LEVEL(ORDERARR(L))) THEN
                                    LVAL = ORDERARR(L)
                                    ORDERARR(L) = ORDERARR(K)
                                    ORDERARR(K) = LVAL
                              END IF
                        END DO
                  END DO

C                 Build ALLNX and DIMS for active (non-zero level) dims
                  NDIMS = 0
                  XTOT = 0
                  DO L = 1, D
                        LVAL = LEVEL(ORDERARR(L))
                        IF (LVAL .GT. 0) THEN
                              NDIMS = NDIMS + 1
                              DIMS(NDIMS) = ORDERARR(L)
                              ALLNX(NDIMS) = 2**LVAL + 1
                              XTOT = XTOT + ALLNX(NDIMS)
                        END IF
                  END DO

C                 Build node array
                  CALL GET_CHEB_NODES(ALLNX, NDIMS, XBUF, XTOT)

C                 Call barycentric step
                  DO K = 1, NINTERP
                        IPTEMP(K) = 0.0D+00
                  END DO
                  CALL BARY_PD_STEP_CB(Z(INDEX), NPTS,
     &                  ALLNX, DIMS, NDIMS,
     &                  XBUF, XTOT, Y, NINTERP, D, IPTEMP)

                  DO K = 1, NINTERP
                        IP(K) = IP(K) + IPTEMP(K)
                  END DO

                  INDEX = INDEX + NPTS
            END IF

      END DO

      RETURN
      END
