C     *******************************************************************
C
C       SPINTERP_CC_SP - Multi-linear interpolation, CC grid, sparse idx.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPINTERP_CC_SP(D, Z, NZ, Y, NINTERP,
     &    INDICESNDIIMS, NSUBGRIDS,
     &    INDICESDIMS, INDICESLEVS, INDICESADDR, NADDR,
     &    SUBGRIDPOINTS, IP)
C     *******************************************************************
C
C       SPINTERP_CC_SP evaluates the CC sparse grid interpolant using
C       the sparse index structure (packed active dimensions/levels).
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
C         Output, DOUBLE PRECISION IP(NINTERP)        - interpolated values
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER D, NZ, NINTERP, NSUBGRIDS, NADDR
      DOUBLE PRECISION Z(NZ)
      DOUBLE PRECISION Y(NINTERP, D)
      INTEGER INDICESNDIIMS(NSUBGRIDS)
      INTEGER INDICESDIMS(NADDR)
      INTEGER INDICESLEVS(NADDR)
      INTEGER INDICESADDR(NSUBGRIDS)
      INTEGER SUBGRIDPOINTS(NSUBGRIDS)
      DOUBLE PRECISION IP(NINTERP)

      INTEGER CI, K, DID, NDIMS, ADDR
      INTEGER LVAL, NPTS
      INTEGER ACTIVEDIMS(50), ACTIVELEV(50)
      INTEGER REPVEC(50), NPTS_DIM(50)
      INTEGER INDEX, INDEX2(50), INDEX3
      INTEGER XP
      DOUBLE PRECISION YT, TEMP, SCALE

      DO K = 1, NINTERP
            IP(K) = 0.0D+00
      END DO

      INDEX = 1

      DO CI = 1, NSUBGRIDS

            NDIMS = INDICESNDIIMS(CI)
            ADDR  = INDICESADDR(CI)
            NPTS  = SUBGRIDPOINTS(CI)

            IF (NDIMS .EQ. 0) THEN
C                 Constant subgrid: add Z(INDEX) to all points
                  DO K = 1, NINTERP
                        IP(K) = IP(K) + Z(INDEX)
                  END DO
                  INDEX = INDEX + NPTS
                  GOTO 200
            END IF

C           Read active dims and levels
            DO DID = 1, NDIMS
                  ACTIVEDIMS(DID) = INDICESDIMS(ADDR + DID - 1)
                  ACTIVELEV(DID)  = INDICESLEVS(ADDR + DID - 1)
            END DO

C           Compute npts_dim and repvec for active dims
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
                  TEMP = 1.0D+00
                  DID = 1

 10               IF (DID .LE. NDIMS) THEN
                        LVAL = ACTIVELEV(DID)
                        YT   = Y(K, ACTIVEDIMS(DID))

                        IF (YT .LT. 0.0D+00) YT = 0.0D+00
                        IF (YT .GT. 1.0D+00) YT = 1.0D+00

                        IF (LVAL .EQ. 0) THEN
                              INDEX2(DID) = 0

                        ELSE IF (LVAL .EQ. 1) THEN
                              IF (YT .EQ. 1.0D+00) THEN
                                    INDEX2(DID) = 1
                              ELSE
                                    XP = INT(YT * 2.0D+00)
                                    IF (XP .EQ. 0) THEN
                                          TEMP = TEMP *
     &                                          2.0D+00*(0.5D+00 - YT)
                                    ELSE
                                          TEMP = TEMP *
     &                                          2.0D+00*(YT - 0.5D+00)
                                    END IF
                                    INDEX2(DID) = XP
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
                                    INDEX2(DID) = XP
                              END IF
                        END IF

                        IF (TEMP .EQ. 0.0D+00) GOTO 100

                        DID = DID + 1
                        GOTO 10
                  END IF

C                 Accumulate: INDEX3 = INDEX + sum_k(REPVEC(k)*INDEX2(k))
                  INDEX3 = INDEX + INDEX2(1)
                  DO DID = 2, NDIMS
                        INDEX3 = INDEX3 + REPVEC(DID) * INDEX2(DID)
                  END DO
                  IP(K) = IP(K) + TEMP * Z(INDEX3)

 100              CONTINUE

            END DO

            INDEX = INDEX + NPTS

 200        CONTINUE
      END DO

      RETURN
      END
