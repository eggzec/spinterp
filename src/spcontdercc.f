C     *******************************************************************
C
C       SPCONT_DERIV_CC - Continuous derivatives, CC grid (full levelseq).
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPCONT_DERIV_CC(D, Z, NZ, Y, NINTERP,
     &    LEVELSEQ, NLEVELS, MAXLEV, IP, IPDER, IPDER2)
C     *******************************************************************
C
C       SPCONT_DERIV_CC computes interpolated values, gradient, and
C       augmented derivative values at neighboring grid cell boundaries
C       for post-processing to produce continuous derivatives.
C
C       IPDER(k,l)  = derivative at Y(k,l) using hat functions
C       IPDER2(k,l) = derivative at the neighboring cell boundary
C                     (used for continuity processing by PP_DERIV)
C
C       Augmented points:
C            if maxlev==1: ytd1=0.25, ytd2=0.75
C            else:         snap to nearest cell boundary
C
C       Parameters:
C
C         Input,  INTEGER D                         - dimension
C         Input,  DOUBLE PRECISION Z(NZ)            - hierarchical surpluses
C         Input,  INTEGER NZ                        - length of Z
C         Input,  DOUBLE PRECISION Y(NINTERP,D)     - query points
C         Input,  INTEGER NINTERP                   - number of query points
C         Input,  INTEGER LEVELSEQ(NLEVELS,D)       - multi-index set
C         Input,  INTEGER NLEVELS                   - rows in LEVELSEQ
C         Input,  INTEGER MAXLEV                    - maximum level
C         Output, DOUBLE PRECISION IP(NINTERP)      - interpolated values
C         Output, DOUBLE PRECISION IPDER(NINTERP,D) - gradient at Y
C         Output, DOUBLE PRECISION IPDER2(NINTERP,D)- gradient at aug pts
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER D, NZ, NINTERP, NLEVELS, MAXLEV
      DOUBLE PRECISION Z(NZ)
      DOUBLE PRECISION Y(NINTERP, D)
      INTEGER LEVELSEQ(NLEVELS, D)
      DOUBLE PRECISION IP(NINTERP)
      DOUBLE PRECISION IPDER(NINTERP, D)
      DOUBLE PRECISION IPDER2(NINTERP, D)

      INTEGER KL, K, L, L2, LVAL, NPTS, INDEX, INDEX2(50), INDEX3
      INTEGER REPVEC(50), XP
      DOUBLE PRECISION YT, YTD, TEMP, TEMP2, SCALE, DIST
      DOUBLE PRECISION TEMPVEC(50), DERVEC(50)
      DOUBLE PRECISION TEMPVEC2(50), DERVEC2(50)
      DOUBLE PRECISION STEPSIZE, HALFSTEP

      DO K = 1, NINTERP
            IP(K) = 0.0D+00
            DO L = 1, D
                  IPDER(K, L)  = 0.0D+00
                  IPDER2(K, L) = 0.0D+00
            END DO
      END DO

      INDEX = 1

      DO KL = 1, NLEVELS

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
                  IF (L .GT. 1) REPVEC(L) = REPVEC(L) * REPVEC(L-1)
            END DO

            DO K = 1, NINTERP

                  DO L = 1, D
                        LVAL = LEVELSEQ(KL, L)
                        YT   = Y(K, L)
                        IF (YT .LT. 0.0D+00) YT = 0.0D+00
                        IF (YT .GT. 1.0D+00) YT = 1.0D+00

C                       Compute augmented point YTD
                        IF (MAXLEV .EQ. 1) THEN
                              IF (YT .LE. 0.5D+00) THEN
                                    YTD = 0.25D+00
                              ELSE
                                    YTD = 0.75D+00
                              END IF
                        ELSE
                              IF (LVAL .EQ. 0) THEN
                                    YTD = YT
                              ELSE
                                    STEPSIZE = 1.0D+00 / DBLE(2**LVAL)
                                    HALFSTEP = STEPSIZE / 2.0D+00
                                    XP = INT(YT / STEPSIZE)
                                    IF (XP .GE. 2**LVAL)
     &                                    XP = 2**LVAL - 1
                                    YTD = DBLE(XP)*STEPSIZE + HALFSTEP
                                    IF (YT .GT. YTD) THEN
                                          YTD = YTD + STEPSIZE
                                    END IF
                                    IF (YTD .GT. 1.0D+00)
     &                                    YTD = 1.0D+00
                              END IF
                        END IF

                        IF (LVAL .EQ. 0) THEN
                              INDEX2(L)   = 0
                              TEMPVEC(L)  = 1.0D+00
                              DERVEC(L)   = 0.0D+00
                              TEMPVEC2(L) = 1.0D+00
                              DERVEC2(L)  = 0.0D+00

                        ELSE IF (LVAL .EQ. 1) THEN
                              IF (YT .EQ. 1.0D+00) THEN
                                    INDEX2(L)   = 1
                                    TEMPVEC(L)  = 1.0D+00
                                    DERVEC(L)   = 2.0D+00
                                    TEMPVEC2(L) = 1.0D+00
                                    DERVEC2(L)  = 2.0D+00
                              ELSE
                                    XP = INT(YT * 2.0D+00)
                                    IF (XP .EQ. 0) THEN
                                          TEMPVEC(L) = 2.0D+00*
     &                                          (0.5D+00 - YT)
                                          DERVEC(L)  = -2.0D+00
                                    ELSE
                                          TEMPVEC(L) = 2.0D+00*
     &                                          (YT - 0.5D+00)
                                          DERVEC(L)  = 2.0D+00
                                    END IF
                                    INDEX2(L) = XP
                                    IF (YTD .EQ. 1.0D+00) THEN
                                          TEMPVEC2(L) = 1.0D+00
                                          DERVEC2(L)  = 2.0D+00
                                    ELSE
                                          XP = INT(YTD * 2.0D+00)
                                          IF (XP .EQ. 0) THEN
                                                TEMPVEC2(L) = 2.0D+00*
     &                                               (0.5D+00-YTD)
                                                DERVEC2(L)  = -2.0D+00
                                          ELSE
                                                TEMPVEC2(L) = 2.0D+00*
     &                                               (YTD-0.5D+00)
                                                DERVEC2(L)  =  2.0D+00
                                          END IF
                                    END IF
                              END IF

                        ELSE
                              SCALE = DBLE(2**LVAL)
                              IF (YT .EQ. 1.0D+00) THEN
                                    INDEX2(L)   = INT(SCALE/2.0D+00)-1
                                    TEMPVEC(L)  = 0.0D+00
                                    DERVEC(L)   = -SCALE
                                    TEMPVEC2(L) = 0.0D+00
                                    DERVEC2(L)  = -SCALE
                              ELSE
                                    XP = INT(YT * SCALE / 2.0D+00)
                                    INDEX2(L) = XP
                                    DIST = YT - DBLE(XP*2+1)/SCALE
                                    TEMPVEC(L) = 1.0D+00
     &                                    - SCALE*DABS(DIST)
                                    IF (DIST .GE. 0.0D+00) THEN
                                          DERVEC(L) = -SCALE
                                    ELSE
                                          DERVEC(L) =  SCALE
                                    END IF
                                    DIST = YTD - DBLE(XP*2+1)/SCALE
                                    TEMPVEC2(L) = 1.0D+00
     &                                    - SCALE*DABS(DIST)
                                    IF (TEMPVEC2(L) .LT. 0.0D+00)
     &                                    TEMPVEC2(L) = 0.0D+00
                                    IF (DIST .GE. 0.0D+00) THEN
                                          DERVEC2(L) = -SCALE
                                    ELSE
                                          DERVEC2(L) =  SCALE
                                    END IF
                              END IF
                        END IF
                  END DO

                  INDEX3 = INDEX + INDEX2(1)
                  DO L = 2, D
                        INDEX3 = INDEX3 + REPVEC(L-1)*INDEX2(L)
                  END DO

                  TEMP = TEMPVEC(1)
                  DO L = 2, D
                        TEMP = TEMP * TEMPVEC(L)
                  END DO
                  IP(K) = IP(K) + TEMP * Z(INDEX3)

                  DO L = 1, D
                        TEMP  = 1.0D+00
                        TEMP2 = 1.0D+00
                        DO L2 = 1, D
                              IF (L .EQ. L2) THEN
                                    TEMP  = TEMP  * DERVEC(L2)
                                    TEMP2 = TEMP2 * DERVEC2(L2)
                              ELSE
                                    TEMP  = TEMP  * TEMPVEC(L2)
                                    TEMP2 = TEMP2 * TEMPVEC2(L2)
                              END IF
                        END DO
                        IPDER(K, L)  = IPDER(K, L)  + TEMP  * Z(INDEX3)
                        IPDER2(K, L) = IPDER2(K, L) + TEMP2 * Z(INDEX3)
                  END DO

            END DO

            INDEX = INDEX + NPTS

      END DO

      RETURN
      END
