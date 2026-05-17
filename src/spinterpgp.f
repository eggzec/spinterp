C     *******************************************************************
C
C       SPINTERP_GP - Polynomial interpolation, Gauss-Patterson grid.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPINTERP_GP(D, Z, NZ, Y, NINTERP,
     &                        LEVELSEQ, NLEVELS, IP)
C     *******************************************************************
C
C       SPINTERP_GP evaluates the Gauss-Patterson sparse grid interpolant
C       at NINTERP query points Y using hierarchical surpluses Z.
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
      INTEGER ALLNX(50), DIMS(50), XTOT
      DOUBLE PRECISION XBUF(16384), WBUF(16384)
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
                  IF (LVAL .GT. 0) THEN
                        NPTS = NPTS * 2**LVAL
                        NDIMS = NDIMS + 1
                        ALLNX(NDIMS) = 2**(LVAL+1) - 1
                        DIMS(NDIMS) = L
                  END IF
            END DO

            IF (NPTS .EQ. 1) THEN
                  DO K = 1, NINTERP
                        IP(K) = IP(K) + Z(INDEX)
                  END DO
                  INDEX = INDEX + 1

            ELSE
                  XTOT = 0
                  DO L = 1, NDIMS
                        XTOT = XTOT + ALLNX(L)
                  END DO

                  CALL GET_GP_NODES(ALLNX, NDIMS, XBUF, XTOT)
                  CALL GET_GP_BARY_W(ALLNX, NDIMS, WBUF, XTOT)

                  DO K = 1, NINTERP
                        IPTEMP(K) = 0.0D+00
                  END DO
                  CALL BARY_PD_STEP_GP(Z(INDEX), NPTS,
     &                  ALLNX, DIMS, NDIMS,
     &                  XBUF, XTOT, Y, NINTERP, D, WBUF, IPTEMP)

                  DO K = 1, NINTERP
                        IP(K) = IP(K) + IPTEMP(K)
                  END DO

                  INDEX = INDEX + NPTS
            END IF

      END DO

      RETURN
      END
