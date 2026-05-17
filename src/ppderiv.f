C     *******************************************************************
C
C       PP_DERIV - Post-processing for continuous derivatives.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE PP_DERIV(IPDER, IPDER2, NINTERP, D,
     &    MAXLEVVEC, Y)
C     *******************************************************************
C
C       PP_DERIV applies post-processing to produce continuously
C       differentiable derivatives by blending IPDER and IPDER2.
C
C       For each point k and dimension l:
C         - Compute ytd1 (cell center) and stepsize from maxlev.
C         - If ipd1 * ipd2 >= 0: linear interpolation between boundaries.
C         - Else if yt <= ytd1 + halfstep: decay from ipd1 toward zero.
C         - Else: ramp from zero toward ipd2.
C
C       Parameters:
C
C         In/out, DOUBLE PRECISION IPDER(NINTERP,D)  - gradient (modified)
C         Input,  DOUBLE PRECISION IPDER2(NINTERP,D) - augmented gradient
C         Input,  INTEGER NINTERP                     - number of points
C         Input,  INTEGER D                           - dimension
C         Input,  INTEGER MAXLEVVEC(D)                - max level per dim
C         Input,  DOUBLE PRECISION Y(NINTERP,D)       - query points
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER NINTERP, D
      DOUBLE PRECISION IPDER(NINTERP, D)
      DOUBLE PRECISION IPDER2(NINTERP, D)
      INTEGER MAXLEVVEC(D)
      DOUBLE PRECISION Y(NINTERP, D)

      INTEGER K, L, MAXLEV, XP
      DOUBLE PRECISION YT, YTD1, STEPSIZE, HALFSTEP
      DOUBLE PRECISION IPD1, IPD2, ALPHA

      DO L = 1, D
            MAXLEV = MAXLEVVEC(L)
            IF (MAXLEV .EQ. 0) GOTO 200

            STEPSIZE = 1.0D+00 / DBLE(2**MAXLEV)
            HALFSTEP = STEPSIZE / 2.0D+00

            DO K = 1, NINTERP
                  YT   = Y(K, L)
                  IPD1 = IPDER(K, L)
                  IPD2 = IPDER2(K, L)

C                 Find the cell center ytd1
                  XP = INT(YT / STEPSIZE)
                  IF (XP .GE. 2**MAXLEV) XP = 2**MAXLEV - 1
                  YTD1 = DBLE(XP)*STEPSIZE + HALFSTEP

                  IF (IPD1 * IPD2 .GE. 0.0D+00) THEN
C                       Same sign: linear blend
                        IF (STEPSIZE .GT. 0.0D+00) THEN
                              ALPHA = (YT - YTD1) / STEPSIZE
                              IF (ALPHA .LT. 0.0D+00) ALPHA = 0.0D+00
                              IF (ALPHA .GT. 1.0D+00) ALPHA = 1.0D+00
                              IPDER(K, L) = (1.0D+00 - ALPHA)*IPD1
     &                                    + ALPHA*IPD2
                        END IF

                  ELSE IF (YT .LE. YTD1 + HALFSTEP) THEN
C                       Left half of cell: taper ipd1 toward zero
                        IF (HALFSTEP .GT. 0.0D+00) THEN
                              ALPHA = (YT - YTD1) / HALFSTEP
                              IF (ALPHA .LT. 0.0D+00) ALPHA = 0.0D+00
                              IF (ALPHA .GT. 1.0D+00) ALPHA = 1.0D+00
                              IPDER(K, L) = IPD1 * (1.0D+00 - ALPHA)
                        END IF

                  ELSE
C                       Right half of cell: ramp from zero toward ipd2
                        IF (HALFSTEP .GT. 0.0D+00) THEN
                              ALPHA = (YT - YTD1 - HALFSTEP) / HALFSTEP
                              IF (ALPHA .LT. 0.0D+00) ALPHA = 0.0D+00
                              IF (ALPHA .GT. 1.0D+00) ALPHA = 1.0D+00
                              IPDER(K, L) = IPD2 * ALPHA
                        END IF
                  END IF

            END DO

 200        CONTINUE
      END DO

      RETURN
      END
