C     *******************************************************************
C
C       SPINTERP_CC - Multi-linear sparse grid interpolation,
C                     Clenshaw-Curtis grid.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPINTERP_CC(D, Z, NZ, Y, NINTERP, LEVELSEQ, NLEVELS,
     &                        IP)
C     *******************************************************************
C
C       SPINTERP_CC evaluates the sparse grid interpolant at NINTERP
C       query points Y using hierarchical surpluses Z and the
C       multi-index set LEVELSEQ.
C
C       The grid is on the unit cube [0,1]^D.  If your domain differs,
C       normalise Y before calling.
C
C       Basis function for dimension k at level LEV:
C            LEV=0 : phi = 1  (constant, no dependency on y)
C            LEV=1 : hat functions centred at 0 and 1 covering [0,1]
C            LEV>=2: hat function of width 2/2^LEV centred at
C                    (2*xp+1)/2^LEV where xp = floor(y*2^(LEV-1))
C
C       Parameters:
C
C         Input,  INTEGER D                    - dimension
C         Input,  DOUBLE PRECISION Z(NZ)       - hierarchical surpluses
C         Input,  INTEGER NZ                   - length of Z
C         Input,  DOUBLE PRECISION Y(NINTERP,D)- query points in [0,1]^D
C         Input,  INTEGER NINTERP              - number of query points
C         Input,  INTEGER LEVELSEQ(NLEVELS,D)  - multi-index set
C         Input,  INTEGER NLEVELS              - rows in LEVELSEQ
C         Output, DOUBLE PRECISION IP(NINTERP) - interpolated values
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER D, NZ, NINTERP, NLEVELS
      DOUBLE PRECISION Z(NZ)
      DOUBLE PRECISION Y(NINTERP, D)
      INTEGER LEVELSEQ(NLEVELS, D)
      DOUBLE PRECISION IP(NINTERP)

      INTEGER KL, K, L, LVAL, NPTS
      INTEGER INDEX, INDEX2(D), INDEX3
      INTEGER REPVEC(D), XP
      DOUBLE PRECISION YT, TEMP, SCALE

      DO K = 1, NINTERP
            IP(K) = 0.0D+00
      END DO

      INDEX = 1

      DO KL = 1, NLEVELS

C           Compute npoints and cumulative repvec (column strides).
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

            DO K = 1, NINTERP
                  TEMP = 1.0D+00
                  L = 1

 10               IF (L .LE. D) THEN
                        LVAL = LEVELSEQ(KL, L)
                        YT   = Y(K, L)

C                       Clamp to [0,1].
                        IF (YT .LT. 0.0D+00) YT = 0.0D+00
                        IF (YT .GT. 1.0D+00) YT = 1.0D+00

                        IF (LVAL .EQ. 0) THEN
                              INDEX2(L) = 0

                        ELSE IF (LVAL .EQ. 1) THEN
                              IF (YT .EQ. 1.0D+00) THEN
                                    INDEX2(L) = 1
                              ELSE
                                    XP = INT(YT * 2.0D+00)
                                    IF (XP .EQ. 0) THEN
                                          TEMP = TEMP *
     &                                          2.0D+00*(0.5D+00 - YT)
                                    ELSE
                                          TEMP = TEMP *
     &                                          2.0D+00*(YT - 0.5D+00)
                                    END IF
                                    INDEX2(L) = XP
                              END IF

                        ELSE
                              IF (YT .EQ. 1.0D+00) THEN
                                    TEMP = 0.0D+00
                              ELSE
                                    SCALE = DBLE(2**LVAL)
                                    XP    = INT(YT * SCALE / 2.0D+00)
                                    TEMP  = TEMP * (1.0D+00 - SCALE *
     &                                      DABS(YT -
     &                                      DBLE(XP*2+1)/SCALE))
                                    INDEX2(L) = XP
                              END IF
                        END IF

                        IF (TEMP .EQ. 0.0D+00) GOTO 20

                        L = L + 1
                        GOTO 10
                  END IF

C                 Accumulate contribution: IP(K) += TEMP * Z(index3).
                  INDEX3 = INDEX + INDEX2(1)
                  DO L = 2, D
                        INDEX3 = INDEX3 + REPVEC(L-1) * INDEX2(L)
                  END DO
                  IP(K) = IP(K) + TEMP * Z(INDEX3)

 20               CONTINUE

            END DO

            INDEX = INDEX + NPTS

      END DO

      RETURN
      END
