C     *******************************************************************
C
C       SPCMPVALS_M - Compute hierarchical surpluses, Maximum-norm grid.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPCMPVALS_M(D, Z, NZ, Y, NY,
     &                        NEWLEVELSEQ, NNEWLEVELS,
     &                        LEVELSEQ, NLEVELS, IP)
C     *******************************************************************
C
C       SPCMPVALS_M computes interpolant at new Maximum-norm grid points Y
C       using surpluses Z from previous levels, for surplus subtraction.
C
C       For lev=0 dimensions the repeat mechanism iterates over all
C       2^nlevelzero combinations of midpoint vs boundary nodes.
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
      INTEGER INDEX, INDEX2(50), INDEX3
      INTEGER REPVEC(50), REPEAT(50), NPTS, NLEVELZERO, NREPEATS
      INTEGER SKIPLEVEL, NNEWPTS(512)
      INTEGER REPSTEP
      DOUBLE PRECISION TEMP, SCALE, YT
      INTEGER KSTART, KEND

      DO K = 1, NY
            IP(K) = 0.0D+00
      END DO

      DO NKL = 1, NNEWLEVELS
            NPTS = 1
            DO L = 1, D
                  LVAL = NEWLEVELSEQ(NKL, L)
                  IF (LVAL .EQ. 0) THEN
                        NPTS = NPTS * 3
                  ELSE
                        NPTS = NPTS * 2**LVAL
                  END IF
            END DO
            NNEWPTS(NKL) = NPTS
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

            KSTART = 1
            DO NKL = 1, NNEWLEVELS
                  KEND = KSTART + NNEWPTS(NKL) - 1

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
                              REPSTEP = 0

 30                           IF (REPSTEP .LE. NREPEATS) THEN
                                    TEMP = 1.0D+00
                                    L = 1

 10                                 IF (L .LE. D) THEN
                                          LVAL = LEVELSEQ(KL, L)
                                          YT = Y(K, L)
                                          IF (YT .LT. 0.0D+00)
     &                                          YT = 0.0D+00
                                          IF (YT .GT. 1.0D+00)
     &                                          YT = 1.0D+00

                                          IF (LVAL .EQ. 0) THEN
                                                IF (REPEAT(L).EQ.0)THEN
                                                      INDEX2(L) = 1
                                                      TEMP = TEMP*
     &                                               (1.0D+00
     &                                               - 2.0D+00
     &                                               *DABS(YT
     &                                               -0.5D+00))
                                                ELSE
                                                      IF (YT.EQ.
     &                                               1.0D+00) THEN
                                                            INDEX2(L)=2
                                                      ELSE
                                                            XP=INT(YT
     &                                                    *2.0D+00)*2
                                                            IF(XP.EQ.0)
     &                                                      THEN
                                                            TEMP=TEMP*
     &                                                      2.0D+00*
     &                                                      (0.5D+00
     &                                                      -YT)
                                                            ELSE
                                                            TEMP=TEMP*
     &                                                      2.0D+00*
     &                                                      (YT
     &                                                      -0.5D+00)
                                                            END IF
                                                            INDEX2(L)
     &                                                      =XP
                                                      END IF
                                                END IF
                                          ELSE IF (YT.EQ.1.0D+00) THEN
                                                TEMP = 0.0D+00
                                          ELSE
                                                SCALE = DBLE(2**LVAL)
                                                XP = INT(YT*SCALE)
                                                TEMP = TEMP*(1.0D+00
     &                                               - 2.0D+00*SCALE
     &                                               *DABS(YT-(DBLE(XP)
     &                                               +0.5D+00)/SCALE))
                                                INDEX2(L) = XP
                                          END IF

                                          IF (TEMP.EQ.0.0D+00)
     &                                          GOTO 20
                                          L = L + 1
                                          GOTO 10
                                    END IF

                                    INDEX3 = INDEX + INDEX2(1)
                                    DO L = 2, D
                                          INDEX3 = INDEX3 +
     &                                          REPVEC(L-1)*INDEX2(L)
                                    END DO
                                    IP(K) = IP(K) + TEMP*Z(INDEX3)

 20                                 REPSTEP = REPSTEP + 1
                                    IF (REPSTEP .LE. NREPEATS) THEN
                                          DO L = 1, D
                                                IF (LEVELSEQ(KL,L)
     &                                               .EQ. 0) THEN
                                                      REPEAT(L) =
     &                                               REPEAT(L) + 1
                                                      IF (REPEAT(L)
     &                                               .GT. 1) THEN
                                                            REPEAT(L)=0
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
                  END IF

                  KSTART = KEND + 1
            END DO

            INDEX = INDEX + NPTS

      END DO

      RETURN
      END
