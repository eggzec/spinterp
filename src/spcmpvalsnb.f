C     *******************************************************************
C
C       SPCMPVALS_NB - Compute hierarchical surpluses, NoBoundary grid.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPCMPVALS_NB(D, Z, NZ, Y, NY,
     &                         NEWLEVELSEQ, NNEWLEVELS,
     &                         LEVELSEQ, NLEVELS, IP)
C     *******************************************************************
C
C       SPCMPVALS_NB computes interpolant at new NoBoundary grid points Y
C       using surpluses Z from previous levels, for surplus subtraction.
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
      INTEGER REPVEC(50), NPTS, SKIPLEVEL
      INTEGER NNEWPTS(512)
      DOUBLE PRECISION TEMP, SCALE, YT
      INTEGER KSTART, KEND

      DO K = 1, NY
            IP(K) = 0.0D+00
      END DO

      DO NKL = 1, NNEWLEVELS
            NPTS = 1
            DO L = 1, D
                  NPTS = NPTS * 2**NEWLEVELSEQ(NKL, L)
            END DO
            NNEWPTS(NKL) = NPTS
      END DO

      INDEX = 1

      DO KL = 1, NLEVELS

            NPTS = 1
            DO L = 1, D
                  LVAL = LEVELSEQ(KL, L)
                  REPVEC(L) = 2**LVAL
                  NPTS = NPTS * REPVEC(L)
                  IF (L .GT. 1) REPVEC(L) = REPVEC(L) * REPVEC(L-1)
            END DO

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
                              TEMP = 1.0D+00
                              L = 1

 10                           IF (L .LE. D) THEN
                                    LVAL = LEVELSEQ(KL, L)
                                    YT   = Y(K, L)

                                    IF (LVAL .EQ. 0) THEN
                                          INDEX2(L) = 0
                                    ELSE
                                          SCALE = DBLE(2**LVAL)
                                          XP = INT(YT * SCALE)
                                          IF (XP .EQ. 0) THEN
                                                TEMP = TEMP*(1.0D+00
     &                                               - 2.0D+00*SCALE
     &                                               *(YT-0.5D+00
     &                                               /SCALE))
                                          ELSE IF (XP .EQ.
     &                                          INT(SCALE)-1) THEN
                                                TEMP = TEMP*(1.0D+00
     &                                               + 2.0D+00*SCALE
     &                                               *(YT-(SCALE
     &                                               -0.5D+00)/SCALE))
                                          ELSE
                                                TEMP = TEMP*(1.0D+00
     &                                               - 2.0D+00*SCALE
     &                                               *DABS(YT-(DBLE(XP)
     &                                               +0.5D+00)/SCALE))
                                          END IF
                                          INDEX2(L) = XP
                                    END IF

                                    L = L + 1
                                    GOTO 10
                              END IF

                              INDEX3 = INDEX + INDEX2(1)
                              DO L = 2, D
                                    INDEX3 = INDEX3 +
     &                                    REPVEC(L-1)*INDEX2(L)
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
