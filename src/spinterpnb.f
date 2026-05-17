C     *******************************************************************
C
C       SPINTERP_NB - Multi-linear interpolation, NoBoundary grid.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPINTERP_NB(D, Z, NZ, Y, NINTERP,
     &                        LEVELSEQ, NLEVELS, IP)
C     *******************************************************************
C
C       SPINTERP_NB evaluates the NoBoundary sparse grid interpolant at
C       NINTERP query points Y using hierarchical surpluses Z.
C
C       Basis function for dimension k at level LEV (scale = 2^lev):
C            lev = 0  : phi = 1  (constant)
C            xp = floor(yt * scale)
C            xp = 0         : 1 - 2*scale*(yt - 0.5/scale)  (left-ext)
C            xp = scale-1   : 1 + 2*scale*(yt - (scale-0.5)/scale)
C            otherwise      : 1 - 2*scale*|yt - (xp+0.5)/scale|
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

      INTEGER KL, K, L, LVAL, NPTS
      INTEGER INDEX, INDEX2(50), INDEX3
      INTEGER REPVEC(50), XP
      DOUBLE PRECISION YT, TEMP, SCALE

      DO K = 1, NINTERP
            IP(K) = 0.0D+00
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

            DO K = 1, NINTERP
                  TEMP = 1.0D+00
                  L = 1

 10               IF (L .LE. D) THEN
                        LVAL = LEVELSEQ(KL, L)
                        YT   = Y(K, L)
                        IF (YT .LT. 0.0D+00) YT = 0.0D+00
                        IF (YT .GT. 1.0D+00) YT = 1.0D+00

                        IF (LVAL .EQ. 0) THEN
                              INDEX2(L) = 0
                        ELSE
                              SCALE = DBLE(2**LVAL)
                              IF (YT .EQ. 1.0D+00) THEN
                                    XP = INT(SCALE) - 1
                              ELSE
                                    XP = INT(YT * SCALE)
                              END IF
                              IF (XP .EQ. 0) THEN
                                    TEMP = TEMP * (1.0D+00
     &                                    - 2.0D+00*SCALE*(YT
     &                                    - 0.5D+00/SCALE))
                              ELSE IF (XP .EQ. INT(SCALE)-1) THEN
                                    TEMP = TEMP * (1.0D+00
     &                                    + 2.0D+00*SCALE*(YT
     &                                    - (SCALE-0.5D+00)/SCALE))
                              ELSE
                                    TEMP = TEMP * (1.0D+00
     &                                    - 2.0D+00*SCALE
     &                                    *DABS(YT-(DBLE(XP)
     &                                    +0.5D+00)/SCALE))
                              END IF
                              INDEX2(L) = XP
                        END IF

                        IF (TEMP .EQ. 0.0D+00) GOTO 20
                        L = L + 1
                        GOTO 10
                  END IF

                  INDEX3 = INDEX + INDEX2(1)
                  DO L = 2, D
                        INDEX3 = INDEX3 + REPVEC(L-1)*INDEX2(L)
                  END DO
                  IP(K) = IP(K) + TEMP * Z(INDEX3)

 20               CONTINUE

            END DO

            INDEX = INDEX + NPTS

      END DO

      RETURN
      END
