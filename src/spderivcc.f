C     *******************************************************************
C
C       SPDERIV_CC - Interpolated values and gradient, Clenshaw-Curtis.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPDERIV_CC(D, Z, NZ, Y, NINTERP,
     &                       LEVELSEQ, NLEVELS, IP, IPDER)
C     *******************************************************************
C
C       SPDERIV_CC computes both the interpolated values and the gradient
C       vectors of the Clenshaw-Curtis sparse grid interpolant.
C
C       The gradient uses the piecewise-linear (hat function) derivative:
C            lev=0: d/dy = 0
C            lev=1: d/dy = -2 (left hat) or +2 (right hat)
C            lev>=2: d/dy = -scale (right side) or +scale (left side)
C                    where scale = 2^lev
C
C       Parameters:
C
C         Input,  INTEGER D                        - dimension
C         Input,  DOUBLE PRECISION Z(NZ)           - hierarchical surpluses
C         Input,  INTEGER NZ                       - length of Z
C         Input,  DOUBLE PRECISION Y(NINTERP,D)    - query points
C         Input,  INTEGER NINTERP                  - number of query points
C         Input,  INTEGER LEVELSEQ(NLEVELS,D)      - multi-index set
C         Input,  INTEGER NLEVELS                  - rows in LEVELSEQ
C         Output, DOUBLE PRECISION IP(NINTERP)     - interpolated values
C         Output, DOUBLE PRECISION IPDER(NINTERP,D)- gradient vectors
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER D, NZ, NINTERP, NLEVELS
      DOUBLE PRECISION Z(NZ)
      DOUBLE PRECISION Y(NINTERP, D)
      INTEGER LEVELSEQ(NLEVELS, D)
      DOUBLE PRECISION IP(NINTERP)
      DOUBLE PRECISION IPDER(NINTERP, D)

      INTEGER KL, K, L, L2, LVAL, NPTS, INDEX, INDEX2(50), INDEX3
      INTEGER REPVEC(50), XP
      DOUBLE PRECISION YT, TEMP, SCALE, DIST
      DOUBLE PRECISION TEMPVEC(50), DERVEC(50)

      DO K = 1, NINTERP
            IP(K) = 0.0D+00
            DO L = 1, D
                  IPDER(K, L) = 0.0D+00
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

                        IF (LVAL .EQ. 0) THEN
                              INDEX2(L)  = 0
                              TEMPVEC(L) = 1.0D+00
                              DERVEC(L)  = 0.0D+00
                        ELSE IF (LVAL .EQ. 1) THEN
                              IF (YT .EQ. 1.0D+00) THEN
                                    INDEX2(L)  = 1
                                    TEMPVEC(L) = 1.0D+00
                                    DERVEC(L)  = 2.0D+00
                              ELSE
                                    XP = INT(YT * 2.0D+00)
                                    IF (XP .EQ. 0) THEN
                                          TEMPVEC(L) = 2.0D+00*
     &                                          (0.5D+00-YT)
                                          DERVEC(L)  =-2.0D+00
                                    ELSE
                                          TEMPVEC(L) = 2.0D+00*
     &                                          (YT-0.5D+00)
                                          DERVEC(L)  = 2.0D+00
                                    END IF
                                    INDEX2(L) = XP
                              END IF
                        ELSE
                              SCALE = DBLE(2**LVAL)
                              IF (YT .EQ. 1.0D+00) THEN
                                    INDEX2(L)  = INT(SCALE/2.0D+00)-1
                                    TEMPVEC(L) = 0.0D+00
                                    DERVEC(L)  = -SCALE
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
                        TEMP = 1.0D+00
                        DO L2 = 1, D
                              IF (L .EQ. L2) THEN
                                    TEMP = TEMP * DERVEC(L2)
                              ELSE
                                    TEMP = TEMP * TEMPVEC(L2)
                              END IF
                        END DO
                        IPDER(K, L) = IPDER(K, L) + TEMP * Z(INDEX3)
                  END DO

            END DO

            INDEX = INDEX + NPTS

      END DO

      RETURN
      END
