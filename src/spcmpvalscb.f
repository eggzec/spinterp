C     *******************************************************************
C
C       SPCMPVALS_CB - Compute hierarchical surpluses, Chebyshev grid.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPCMPVALS_CB(D, Z, NZ, Y, NY,
     &                         NEWLEVELSEQ, NNEWLEVELS,
     &                         LEVELSEQ, NLEVELS, IP)
C     *******************************************************************
C
C       SPCMPVALS_CB computes the interpolant at new Chebyshev grid
C       points Y using old surpluses Z from previous levels.
C
C       Parameters:
C
C         Input,  INTEGER D                          - dimension
C         Input,  DOUBLE PRECISION Z(NZ)             - old surpluses
C         Input,  INTEGER NZ                         - length of Z
C         Input,  DOUBLE PRECISION Y(NY,D)           - new grid points
C         Input,  INTEGER NY                         - rows in Y
C         Input,  INTEGER NEWLEVELSEQ(NNEWLEVELS,D)  - new multi-indices
C         Input,  INTEGER NNEWLEVELS                 - rows in NEWLEVELSEQ
C         Input,  INTEGER LEVELSEQ(NLEVELS,D)        - old multi-indices
C         Input,  INTEGER NLEVELS                    - rows in LEVELSEQ
C         Output, DOUBLE PRECISION IP(NY)            - interpolant at Y
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER D, NZ, NY, NNEWLEVELS, NLEVELS
      DOUBLE PRECISION Z(NZ)
      DOUBLE PRECISION Y(NY, D)
      INTEGER NEWLEVELSEQ(NNEWLEVELS, D)
      INTEGER LEVELSEQ(NLEVELS, D)
      DOUBLE PRECISION IP(NY)

      INTEGER KL, NKL, K, L, LVAL, NPTS, NDIMS, INDEX
      INTEGER ALLNX(50), DIMS(50), LEVEL(50), ORDERARR(50)
      INTEGER NNEWPTS(512), SKIPLEVEL, XTOT
      INTEGER KSTART, KEND
      DOUBLE PRECISION XBUF(16384)
      DOUBLE PRECISION IPTEMP(65536)

      DO K = 1, NY
            IP(K) = 0.0D+00
      END DO

C     Pre-count new subgrid points
      DO NKL = 1, NNEWLEVELS
            NPTS = 1
            DO L = 1, D
                  LVAL = NEWLEVELSEQ(NKL, L)
                  IF (LVAL .EQ. 0) THEN
                  ELSE IF (LVAL .LE. 2) THEN
                        NPTS = NPTS * 2
                  ELSE
                        NPTS = NPTS * 2**(LVAL-1)
                  END IF
            END DO
            NNEWPTS(NKL) = NPTS
      END DO

      INDEX = 1

      DO KL = 1, NLEVELS

            NPTS = 1
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

            IF (NPTS .GT. 1) THEN
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

                  CALL GET_CHEB_NODES(ALLNX, NDIMS, XBUF, XTOT)

                  KSTART = 1
                  DO NKL = 1, NNEWLEVELS
                        KEND = KSTART + NNEWPTS(NKL) - 1

                        SKIPLEVEL = 0
                        DO L = 1, D
                              IF (LEVELSEQ(KL,L) .GT.
     &                              NEWLEVELSEQ(NKL,L)) THEN
                                    SKIPLEVEL = 1
                                    GOTO 80
                              END IF
                        END DO

 80                     IF (SKIPLEVEL .EQ. 0) THEN
                              DO K = 1, KEND - KSTART + 1
                                    IPTEMP(K) = 0.0D+00
                              END DO
                              CALL BARY_PD_STEP_CB(
     &                              Z(INDEX), NPTS,
     &                              ALLNX, DIMS, NDIMS,
     &                              XBUF, XTOT,
     &                              Y(KSTART,1),
     &                              KEND-KSTART+1, D,
     &                              IPTEMP)
                              DO K = KSTART, KEND
                                    IP(K) = IP(K) +
     &                                    IPTEMP(K-KSTART+1)
                              END DO
                        END IF

                        KSTART = KEND + 1
                  END DO

            ELSE
C                 Single-node subgrid: constant contribution
                  DO K = 1, NY
                        IP(K) = IP(K) + Z(INDEX)
                  END DO
            END IF

            INDEX = INDEX + NPTS

      END DO

      RETURN
      END
