C     *******************************************************************
C
C       SPINTERP_M - Multi-linear interpolation, Maximum-norm grid.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPINTERP_M(D, Z, NZ, Y, NINTERP,
     &                       LEVELSEQ, NLEVELS, IP)
C     *******************************************************************
C
C       SPINTERP_M evaluates the Maximum-norm sparse grid interpolant.
C
C       For lev=0 each dimension contributes 3 nodes (at 0, 0.5, 1).
C       All 2^nlevelzero combinations of midpoint vs boundary nodes are
C       accumulated for each query point by an inner repeat loop.
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

      INTEGER KL, K, L, LVAL, NPTS, NLEVELZERO, NREPEATS
      INTEGER INDEX, INDEX2(50), INDEX3, REPVEC(50)
      INTEGER REPEAT(50), REPSTEP, XP
      DOUBLE PRECISION YT, TEMP, SCALE

      DO K = 1, NINTERP
            IP(K) = 0.0D+00
      END DO

      INDEX = 1

      DO KL = 1, NLEVELS

            NPTS = 1
            NLEVELZERO = 0
            DO L = 1, D
                  LVAL = LEVELSEQ(KL, L)
                  IF (LVAL .EQ. 0) THEN
                        REPVEC(L) = 3
                        NLEVELZERO = NLEVELZERO + 1
                  ELSE
                        REPVEC(L) = 2**LVAL
                  END IF
                  REPEAT(L) = 0
                  NPTS = NPTS * REPVEC(L)
                  IF (L .GT. 1) REPVEC(L) = REPVEC(L) * REPVEC(L-1)
            END DO

            NREPEATS = 2**NLEVELZERO - 1

            DO K = 1, NINTERP
                  REPSTEP = 0

 30               IF (REPSTEP .LE. NREPEATS) THEN
                        TEMP = 1.0D+00
                        L = 1

 10                     IF (L .LE. D) THEN
                              LVAL = LEVELSEQ(KL, L)
                              YT   = Y(K, L)
                              IF (YT .LT. 0.0D+00) YT = 0.0D+00
                              IF (YT .GT. 1.0D+00) YT = 1.0D+00

                              IF (LVAL .EQ. 0) THEN
                                    IF (REPEAT(L) .EQ. 0) THEN
                                          INDEX2(L) = 1
                                          TEMP = TEMP*(1.0D+00
     &                                          - 2.0D+00
     &                                          *DABS(YT-0.5D+00))
                                    ELSE
                                          IF (YT .EQ. 1.0D+00) THEN
                                                INDEX2(L) = 2
                                          ELSE
                                                XP = INT(YT*2.0D+00)*2
                                                IF (XP .EQ. 0) THEN
                                                      TEMP = TEMP *
     &                                               2.0D+00*
     &                                               (0.5D+00-YT)
                                                ELSE
                                                      TEMP = TEMP *
     &                                               2.0D+00*
     &                                               (YT-0.5D+00)
                                                END IF
                                                INDEX2(L) = XP
                                          END IF
                                    END IF
                              ELSE IF (YT .EQ. 1.0D+00) THEN
                                    TEMP = 0.0D+00
                              ELSE
                                    SCALE = DBLE(2**LVAL)
                                    XP = INT(YT * SCALE)
                                    TEMP = TEMP*(1.0D+00
     &                                    - 2.0D+00*SCALE
     &                                    *DABS(YT-(DBLE(XP)
     &                                    +0.5D+00)/SCALE))
                                    INDEX2(L) = XP
                              END IF

                              IF (TEMP .EQ. 0.0D+00) GOTO 20
                              L = L + 1
                              GOTO 10
                        END IF

                        IF (TEMP .GT. 0.0D+00) THEN
                              INDEX3 = INDEX + INDEX2(1)
                              DO L = 2, D
                                    INDEX3 = INDEX3
     &                                    + REPVEC(L-1)*INDEX2(L)
                              END DO
                              IP(K) = IP(K) + TEMP * Z(INDEX3)
                        END IF

 20                     REPSTEP = REPSTEP + 1
                        IF (REPSTEP .LE. NREPEATS) THEN
                              DO L = 1, D
                                    IF (LEVELSEQ(KL,L) .EQ. 0) THEN
                                          REPEAT(L) = REPEAT(L) + 1
                                          IF (REPEAT(L) .GT. 1) THEN
                                                REPEAT(L) = 0
                                          ELSE
                                                GOTO 30
                                          END IF
                                    END IF
                              END DO
                        END IF
                        GOTO 30
                  END IF

                  DO L = 1, D
                        REPEAT(L) = 0
                  END DO

            END DO

            INDEX = INDEX + NPTS

      END DO

      RETURN
      END
