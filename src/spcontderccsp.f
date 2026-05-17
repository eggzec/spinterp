C     *******************************************************************
C
C       SPCONT_DERIV_CC_SP - Continuous derivatives, CC, sparse index.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPCONT_DERIV_CC_SP(D, Z, NZ, Y, NINTERP,
     &    INDICESNDIIMS, NSUBGRIDS,
     &    INDICESDIMS, INDICESLEVS, INDICESADDR, NADDR,
     &    SUBGRIDPOINTS, MAXLEV, IP, IPDER, IPDER2)
C     *******************************************************************
C
C       SPCONT_DERIV_CC_SP computes interpolated values, gradient, and
C       augmented derivative values using the sparse index structure.
C       Same as SPCONT_DERIV_CC but uses sparse format for subgrid data.
C
C       Parameters:
C
C         Input,  INTEGER D                          - dimension
C         Input,  DOUBLE PRECISION Z(NZ)             - hierarchical surpluses
C         Input,  INTEGER NZ                         - length of Z
C         Input,  DOUBLE PRECISION Y(NINTERP,D)      - query points
C         Input,  INTEGER NINTERP                    - number of query points
C         Input,  INTEGER INDICESNDIIMS(NSUBGRIDS)   - ndims per subgrid
C         Input,  INTEGER NSUBGRIDS                  - number of subgrids
C         Input,  INTEGER INDICESDIMS(NADDR)          - packed dim indices
C         Input,  INTEGER INDICESLEVS(NADDR)          - packed levels
C         Input,  INTEGER INDICESADDR(NSUBGRIDS)      - start addr per subgrid
C         Input,  INTEGER NADDR                       - length of packed arrays
C         Input,  INTEGER SUBGRIDPOINTS(NSUBGRIDS)    - points per subgrid
C         Input,  INTEGER MAXLEV                      - maximum level
C         Output, DOUBLE PRECISION IP(NINTERP)        - interpolated values
C         Output, DOUBLE PRECISION IPDER(NINTERP,D)   - gradient at Y
C         Output, DOUBLE PRECISION IPDER2(NINTERP,D)  - gradient at aug pts
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER D, NZ, NINTERP, NSUBGRIDS, NADDR, MAXLEV
      DOUBLE PRECISION Z(NZ)
      DOUBLE PRECISION Y(NINTERP, D)
      INTEGER INDICESNDIIMS(NSUBGRIDS)
      INTEGER INDICESDIMS(NADDR)
      INTEGER INDICESLEVS(NADDR)
      INTEGER INDICESADDR(NSUBGRIDS)
      INTEGER SUBGRIDPOINTS(NSUBGRIDS)
      DOUBLE PRECISION IP(NINTERP)
      DOUBLE PRECISION IPDER(NINTERP, D)
      DOUBLE PRECISION IPDER2(NINTERP, D)

      INTEGER CI, K, DID, L2, NDIMS, ADDR
      INTEGER LVAL, NPTS
      INTEGER ACTIVEDIMS(50), ACTIVELEV(50)
      INTEGER REPVEC(50), NPTS_DIM(50)
      INTEGER INDEX, INDEX2(50), INDEX3
      INTEGER XP
      DOUBLE PRECISION YT, YTD, TEMP, TEMP2, SCALE, DIST
      DOUBLE PRECISION TEMPVEC(50), DERVEC(50)
      DOUBLE PRECISION TEMPVEC2(50), DERVEC2(50)
      DOUBLE PRECISION STEPSIZE

      DO K = 1, NINTERP
            IP(K) = 0.0D+00
            DO L2 = 1, D
                  IPDER(K, L2)  = 0.0D+00
                  IPDER2(K, L2) = 0.0D+00
            END DO
      END DO

      INDEX = 1

      DO CI = 1, NSUBGRIDS

            NDIMS = INDICESNDIIMS(CI)
            ADDR  = INDICESADDR(CI)
            NPTS  = SUBGRIDPOINTS(CI)

            IF (NDIMS .EQ. 0) THEN
                  DO K = 1, NINTERP
                        IP(K) = IP(K) + Z(INDEX)
                  END DO
                  INDEX = INDEX + NPTS
                  GOTO 200
            END IF

            DO DID = 1, NDIMS
                  ACTIVEDIMS(DID) = INDICESDIMS(ADDR + DID - 1)
                  ACTIVELEV(DID)  = INDICESLEVS(ADDR + DID - 1)
            END DO

            DO DID = 1, NDIMS
                  LVAL = ACTIVELEV(DID)
                  IF (LVAL .EQ. 0) THEN
                        NPTS_DIM(DID) = 1
                  ELSE IF (LVAL .LT. 3) THEN
                        NPTS_DIM(DID) = 2
                  ELSE
                        NPTS_DIM(DID) = 2**(LVAL-1)
                  END IF
            END DO
            REPVEC(1) = 1
            DO DID = 2, NDIMS
                  REPVEC(DID) = REPVEC(DID-1) * NPTS_DIM(DID-1)
            END DO

            DO K = 1, NINTERP

                  DO DID = 1, NDIMS
                        LVAL = ACTIVELEV(DID)
                        YT   = Y(K, ACTIVEDIMS(DID))
                        IF (YT .LT. 0.0D+00) YT = 0.0D+00
                        IF (YT .GT. 1.0D+00) YT = 1.0D+00

                        IF (MAXLEV .EQ. 1) THEN
                              IF (YT .LE. 0.5D+00) THEN
                                    YTD = 0.25D+00
                              ELSE
                                    YTD = 0.75D+00
                              END IF
                        ELSE IF (LVAL .EQ. 0) THEN
                              YTD = YT
                        ELSE
                              STEPSIZE = 1.0D+00 / DBLE(2**LVAL)
                              XP = INT(YT / STEPSIZE)
                              IF (XP .GE. 2**LVAL) XP = 2**LVAL - 1
                              YTD = DBLE(XP)*STEPSIZE + STEPSIZE/2.0D+00
                              IF (YT .GT. YTD) YTD = YTD + STEPSIZE
                              IF (YTD .GT. 1.0D+00) YTD = 1.0D+00
                        END IF

                        IF (LVAL .EQ. 0) THEN
                              INDEX2(DID)   = 0
                              TEMPVEC(DID)  = 1.0D+00
                              DERVEC(DID)   = 0.0D+00
                              TEMPVEC2(DID) = 1.0D+00
                              DERVEC2(DID)  = 0.0D+00

                        ELSE IF (LVAL .EQ. 1) THEN
                              IF (YT .EQ. 1.0D+00) THEN
                                    INDEX2(DID)   = 1
                                    TEMPVEC(DID)  = 1.0D+00
                                    DERVEC(DID)   = 2.0D+00
                                    TEMPVEC2(DID) = 1.0D+00
                                    DERVEC2(DID)  = 2.0D+00
                              ELSE
                                    XP = INT(YT * 2.0D+00)
                                    IF (XP .EQ. 0) THEN
                                          TEMPVEC(DID) = 2.0D+00*
     &                                          (0.5D+00-YT)
                                          DERVEC(DID)  = -2.0D+00
                                    ELSE
                                          TEMPVEC(DID) = 2.0D+00*
     &                                          (YT-0.5D+00)
                                          DERVEC(DID)  =  2.0D+00
                                    END IF
                                    INDEX2(DID) = XP
                                    XP = INT(YTD * 2.0D+00)
                                    IF (XP .EQ. 0) THEN
                                          TEMPVEC2(DID) = 2.0D+00*
     &                                          (0.5D+00-YTD)
                                          DERVEC2(DID)  = -2.0D+00
                                    ELSE
                                          TEMPVEC2(DID) = 2.0D+00*
     &                                          (YTD-0.5D+00)
                                          DERVEC2(DID)  =  2.0D+00
                                    END IF
                              END IF
                        ELSE
                              SCALE = DBLE(2**LVAL)
                              IF (YT .EQ. 1.0D+00) THEN
                                    INDEX2(DID)   = INT(SCALE/2.0D+00)-1
                                    TEMPVEC(DID)  = 0.0D+00
                                    DERVEC(DID)   = -SCALE
                                    TEMPVEC2(DID) = 0.0D+00
                                    DERVEC2(DID)  = -SCALE
                              ELSE
                                    XP = INT(YT * SCALE / 2.0D+00)
                                    INDEX2(DID) = XP
                                    DIST = YT - DBLE(XP*2+1)/SCALE
                                    TEMPVEC(DID) = 1.0D+00
     &                                    - SCALE*DABS(DIST)
                                    IF (DIST .GE. 0.0D+00) THEN
                                          DERVEC(DID) = -SCALE
                                    ELSE
                                          DERVEC(DID) =  SCALE
                                    END IF
                                    DIST = YTD - DBLE(XP*2+1)/SCALE
                                    TEMPVEC2(DID) = 1.0D+00
     &                                    - SCALE*DABS(DIST)
                                    IF (TEMPVEC2(DID) .LT. 0.0D+00)
     &                                    TEMPVEC2(DID) = 0.0D+00
                                    IF (DIST .GE. 0.0D+00) THEN
                                          DERVEC2(DID) = -SCALE
                                    ELSE
                                          DERVEC2(DID) =  SCALE
                                    END IF
                              END IF
                        END IF
                  END DO

                  INDEX3 = INDEX + INDEX2(1)
                  DO DID = 2, NDIMS
                        INDEX3 = INDEX3 + REPVEC(DID)*INDEX2(DID)
                  END DO

                  TEMP = TEMPVEC(1)
                  DO DID = 2, NDIMS
                        TEMP = TEMP * TEMPVEC(DID)
                  END DO
                  IP(K) = IP(K) + TEMP * Z(INDEX3)

                  DO DID = 1, NDIMS
                        TEMP  = 1.0D+00
                        TEMP2 = 1.0D+00
                        DO L2 = 1, NDIMS
                              IF (L2 .EQ. DID) THEN
                                    TEMP  = TEMP  * DERVEC(L2)
                                    TEMP2 = TEMP2 * DERVEC2(L2)
                              ELSE
                                    TEMP  = TEMP  * TEMPVEC(L2)
                                    TEMP2 = TEMP2 * TEMPVEC2(L2)
                              END IF
                        END DO
                        IPDER(K, ACTIVEDIMS(DID)) =
     &                        IPDER(K, ACTIVEDIMS(DID)) +
     &                        TEMP * Z(INDEX3)
                        IPDER2(K, ACTIVEDIMS(DID)) =
     &                        IPDER2(K, ACTIVEDIMS(DID)) +
     &                        TEMP2 * Z(INDEX3)
                  END DO

            END DO

            INDEX = INDEX + NPTS

 200        CONTINUE
      END DO

      RETURN
      END
