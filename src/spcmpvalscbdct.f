C     *******************************************************************
C
C       SPCMPVALS_CB_DCT - CB surplus computation using DCT.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPCMPVALS_CB_DCT(D, Z, NZ, Y, NY,
     &    NEWLEVELSEQ, NNEWLEVELS,
     &    LEVELSEQ, NLEVELS, IP)
C     *******************************************************************
C
C       SPCMPVALS_CB_DCT computes the hierarchical surplus increments for
C       Chebyshev grid using DCT-based upsampling instead of barycentric
C       interpolation.
C
C       For each old subgrid (LEVELSEQ row kl), for each new subgrid
C       (NEWLEVELSEQ row nkl) that it contributes to, calls SP_DCT_UP_STEP.
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

      INTEGER KL, NKL, K, L, LVAL, NPTS, NDIMS
      INTEGER INDEX, KSTART, KEND
      INTEGER NNEWPTS(512), SKIPLEVEL
      INTEGER OLDLEV(50), NEWLEV(50), DIMS(50)
      DOUBLE PRECISION IPBUF(16384)

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

C           Count old subgrid points and collect active dims/levels
            NPTS  = 1
            NDIMS = 0
            DO L = 1, D
                  LVAL = LEVELSEQ(KL, L)
                  IF (LVAL .GT. 0) THEN
                        NDIMS = NDIMS + 1
                        OLDLEV(NDIMS) = LVAL
                        DIMS(NDIMS)   = L
                        IF (LVAL .LE. 2) THEN
                              NPTS = NPTS * 2
                        ELSE
                              NPTS = NPTS * 2**(LVAL-1)
                        END IF
                  END IF
            END DO

            IF (NPTS .GT. 1) THEN

                  KSTART = 1
                  DO NKL = 1, NNEWLEVELS
                        KEND = KSTART + NNEWPTS(NKL) - 1

                        SKIPLEVEL = 0
                        DO L = 1, D
                              IF (LEVELSEQ(KL,L) .GT.
     &                            NEWLEVELSEQ(NKL,L)) THEN
                                    SKIPLEVEL = 1
                                    GOTO 80
                              END IF
                        END DO

 80                     IF (SKIPLEVEL .EQ. 0) THEN
C                             Collect new levels for this NKL
                              DO L = 1, NDIMS
                                    NEWLEV(L) = NEWLEVELSEQ(NKL,DIMS(L))
                              END DO

                              DO K = 1, KEND - KSTART + 1
                                    IPBUF(K) = 0.0D+00
                              END DO

                              CALL SP_DCT_UP_STEP(Z(INDEX), NPTS,
     &                              OLDLEV, NDIMS,
     &                              NEWLEV, NDIMS, DIMS,
     &                              IPBUF, KEND-KSTART+1)

                              DO K = KSTART, KEND
                                    IP(K) = IP(K) +
     &                                    IPBUF(K-KSTART+1)
                              END DO
                        END IF

                        KSTART = KEND + 1
                  END DO

            ELSE
C                 Single-node constant subgrid
                  DO K = 1, NY
                        IP(K) = IP(K) + Z(INDEX)
                  END DO
            END IF

            INDEX = INDEX + NPTS

      END DO

      RETURN
      END
