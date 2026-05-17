C     *******************************************************************
C
C       SPCMPVALS_CC - Compute hierarchical surpluses,
C                      Clenshaw-Curtis sparse grid.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPCMPVALS_CC(D, Z, NZ, Y, NY,
     &                         NEWLEVELSEQ, NNEWLEVELS,
     &                         LEVELSEQ, NLEVELS, IP)
C     *******************************************************************
C
C       SPCMPVALS_CC computes the hierarchical surplus increments for
C       the NEW subgrid levels (NEWLEVELSEQ) using the already-computed
C       surpluses from previous levels (Z, LEVELSEQ).
C
C       Caller subtracts IP from the actual function values at Y to
C       obtain the true surpluses for the new subgrid.
C
C       Algorithm mirrors MATLAB spcmpvalscc:
C         For each old subgrid (LEVELSEQ row kl):
C           For each new subgrid (NEWLEVELSEQ row nkl):
C             Skip if any level(l) > newlevel(l) (no contribution).
C             Otherwise accumulate hat-function weights times surplus.
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

      INTEGER KL, NKL, K, L
      INTEGER LVAL, XP
      INTEGER INDEX, INDEX2(D), INDEX3
      INTEGER REPVEC(D), NPTS, SKIPLEVEL
      INTEGER NNEWPTS(NNEWLEVELS)
      DOUBLE PRECISION TEMP, SCALE, YT
      INTEGER KSTART, KEND

      DO K = 1, NY
            IP(K) = 0.0D+00
      END DO

C     Pre-compute point counts for each new subgrid.
      DO NKL = 1, NNEWLEVELS
            NPTS = 1
            DO L = 1, D
                  LVAL = NEWLEVELSEQ(NKL, L)
                  IF (LVAL .EQ. 0) THEN
                        CONTINUE
                  ELSE IF (LVAL .LT. 3) THEN
                        NPTS = NPTS * 2
                  ELSE
                        NPTS = NPTS * 2**(LVAL-1)
                  END IF
            END DO
            NNEWPTS(NKL) = NPTS
      END DO

      INDEX = 1

      DO KL = 1, NLEVELS

C           Build repvec (cumulative strides) for old subgrid kl.
            NPTS = 1
            DO L = 1, D
                  LVAL = LEVELSEQ(KL, L)
                  IF (LVAL .EQ. 0) THEN
                        REPVEC(L) = 1
                  ELSE IF (LVAL .LT. 3) THEN
                        REPVEC(L) = 2
                  ELSE
                        REPVEC(L) = 2**(LVAL-1)
                  END IF
                  NPTS = NPTS * REPVEC(L)
                  IF (L .GT. 1) THEN
                        REPVEC(L) = REPVEC(L) * REPVEC(L-1)
                  END IF
            END DO

C           Walk through new subgrids to determine which Y rows this
C           old subgrid contributes to.
            KSTART = 1
            DO NKL = 1, NNEWLEVELS
                  KEND = KSTART + NNEWPTS(NKL) - 1

C                 Skip if old level exceeds new level in any dim.
                  SKIPLEVEL = 0
                  DO L = 1, D
                        IF (LEVELSEQ(KL,L) .GT.
     &                      NEWLEVELSEQ(NKL,L)) THEN
                              SKIPLEVEL = 1
                              GOTO 80
                        END IF
                  END DO

 80               IF (SKIPLEVEL .EQ. 0) THEN
                        DO K = KSTART, KEND
                              TEMP = 1.0D+00
                              L = 1

 10                           IF (L .LE. D) THEN
                                    LVAL = LEVELSEQ(KL, L)
                                    YT   = Y(K, L)

                                    IF (LVAL .EQ. 0) THEN
                                          INDEX2(L) = 0

                                    ELSE IF (LVAL .EQ. 1) THEN
                                          IF (YT .EQ. 1.0D+00) THEN
                                                INDEX2(L) = 1
                                          ELSE
                                                XP = INT(YT*2.0D+00)
                                                IF (XP .EQ. 0) THEN
                                                      TEMP = TEMP *
     &                                               2.0D+00 *
     &                                               (0.5D+00 - YT)
                                                ELSE
                                                      TEMP = TEMP *
     &                                               2.0D+00 *
     &                                               (YT - 0.5D+00)
                                                END IF
                                                INDEX2(L) = XP
                                          END IF

                                    ELSE
                                          SCALE = DBLE(2**LVAL)
                                          XP = INT(YT*SCALE/2.0D+00)
                                          TEMP = TEMP * (1.0D+00
     &                                         - SCALE * DABS(YT
     &                                         - DBLE(XP*2+1)/SCALE))
                                          INDEX2(L) = XP
                                    END IF

                                    L = L + 1
                                    GOTO 10
                              END IF

                              INDEX3 = INDEX + INDEX2(1)
                              DO L = 2, D
                                    INDEX3 = INDEX3 +
     &                                    REPVEC(L-1) * INDEX2(L)
                              END DO
                              IP(K) = IP(K) + TEMP * Z(INDEX3)

                        END DO
                  END IF

                  KSTART = KEND + 1
            END DO

            INDEX = INDEX + NPTS

      END DO

      RETURN
      END
