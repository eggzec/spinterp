C     *******************************************************************
C
C       SPCMPVALS_CC_SP - CC hierarchical surpluses, sparse index struct.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPCMPVALS_CC_SP(D, Z, NZ, Y, NY,
     &    INDICESNDIIMS, NSUBGRIDS,
     &    INDICESDIMS, INDICESLEVS, INDICESADDR, NADDR,
     &    BACKWARDNEIGHBORS, FORWARDNEIGHBORS, NFWD,
     &    SUBGRIDPOINTS, SUBGRIDADDR,
     &    FROMINDEX, TOINDEX, IP)
C     *******************************************************************
C
C       SPCMPVALS_CC_SP computes the hierarchical surplus increments for
C       new subgrids FROMINDEX..TOINDEX using a sparse index structure
C       and neighbor traversal (ancestor walk).
C
C       For each new subgrid, walks back through ancestor subgrids using
C       BACKWARDNEIGHBORS / FORWARDNEIGHBORS, applies hat function
C       accumulation from each ancestor.
C
C       Parameters:
C
C         Input,  INTEGER D                          - dimension
C         Input,  DOUBLE PRECISION Z(NZ)             - old surpluses
C         Input,  INTEGER NZ                         - length of Z
C         Input,  DOUBLE PRECISION Y(NY,D)           - new grid points
C         Input,  INTEGER NY                         - rows in Y
C         Input,  INTEGER INDICESNDIIMS(NSUBGRIDS)   - ndims per subgrid
C         Input,  INTEGER NSUBGRIDS                  - number of subgrids
C         Input,  INTEGER INDICESDIMS(NADDR)          - packed dim indices
C         Input,  INTEGER INDICESLEVS(NADDR)          - packed levels
C         Input,  INTEGER INDICESADDR(NSUBGRIDS)      - start addr per subgrid
C         Input,  INTEGER NADDR                       - length of packed arrays
C         Input,  INTEGER BACKWARDNEIGHBORS(NADDR)    - back neighbor per dim
C         Input,  INTEGER FORWARDNEIGHBORS(NFWD)      - fwd neighbor flat array
C         Input,  INTEGER NFWD                        - length of fwd array
C         Input,  INTEGER SUBGRIDPOINTS(NSUBGRIDS)    - points per subgrid
C         Input,  INTEGER SUBGRIDADDR(NSUBGRIDS)      - start Z addr per subgrid
C         Input,  INTEGER FROMINDEX                   - first new subgrid
C         Input,  INTEGER TOINDEX                     - last new subgrid
C         Output, DOUBLE PRECISION IP(NY)             - interpolant at Y
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER D, NZ, NY, NSUBGRIDS, NADDR, NFWD
      INTEGER FROMINDEX, TOINDEX
      DOUBLE PRECISION Z(NZ)
      DOUBLE PRECISION Y(NY, D)
      INTEGER INDICESNDIIMS(NSUBGRIDS)
      INTEGER INDICESDIMS(NADDR)
      INTEGER INDICESLEVS(NADDR)
      INTEGER INDICESADDR(NSUBGRIDS)
      INTEGER BACKWARDNEIGHBORS(NADDR)
      INTEGER FORWARDNEIGHBORS(NFWD)
      INTEGER SUBGRIDPOINTS(NSUBGRIDS)
      INTEGER SUBGRIDADDR(NSUBGRIDS)
      DOUBLE PRECISION IP(NY)

C     Local variables
      INTEGER NKL, CI, DID, K, L
      INTEGER NDIMS, ADDR
      INTEGER NEWLEV(50), DIMVEC(50)
      INTEGER CURLEV(50), CURINDEX
      INTEGER LVAL, NPTS
      INTEGER ACTIVEDIMS(50), ACTIVELEV(50)
      INTEGER REPVEC(50), NPTS_DIM(50)
      INTEGER INDEX2(50), INDEX3, INDEX
      INTEGER XP, KSTART, KEND
      DOUBLE PRECISION YT, TEMP, SCALE
      INTEGER STACKIDX(100), STACKDID(100), NSTACK
      INTEGER MAXANCESTORS, IANC
      INTEGER ANCINDEX, ANCADDR, ANCNDIMS, ANCNPTS
      INTEGER ANCACTIVEDIMS(50), ANCACTIVELEV(50)
      INTEGER SKIPLEVEL

      DO K = 1, NY
            IP(K) = 0.0D+00
      END DO

C     Compute kstart/kend offsets for new subgrids
      KSTART = 1
      DO CI = 1, FROMINDEX - 1
            KSTART = KSTART + SUBGRIDPOINTS(CI)
      END DO

      DO NKL = FROMINDEX, TOINDEX

            NDIMS = INDICESNDIIMS(NKL)
            ADDR  = INDICESADDR(NKL)
            NPTS  = SUBGRIDPOINTS(NKL)
            KEND  = KSTART + NPTS - 1

C           Read new subgrid's active dims and levels
            DO DID = 1, NDIMS
                  DIMVEC(DID)  = INDICESDIMS(ADDR + DID - 1)
                  NEWLEV(DID)  = INDICESLEVS(ADDR + DID - 1)
                  CURLEV(DID)  = 0
            END DO

C           Walk all ancestors using nested loop descent
C           Ancestor subgrid starts at the (0,0,...,0) level subgrid
C           and climbs up through all levels up to NEWLEV.
C           Use a simple exhaustive approach: iterate over all old subgrids
C           and check if they are ancestors of NKL.

            DO CI = 1, NSUBGRIDS
C                 Skip new subgrids themselves
                  IF (CI .GE. FROMINDEX) GOTO 300

                  ANCNDIMS = INDICESNDIIMS(CI)
                  ANCADDR  = INDICESADDR(CI)
                  ANCNPTS  = SUBGRIDPOINTS(CI)

C                 Check if CI is an ancestor: all its active dims must have
C                 levels <= NEWLEV in the same dim, and only dims in DIMVEC
                  SKIPLEVEL = 0
                  DO L = 1, ANCNDIMS
                        ANCACTIVEDIMS(L) = INDICESDIMS(ANCADDR + L - 1)
                        ANCACTIVELEV(L)  = INDICESLEVS(ANCADDR + L - 1)
                  END DO

C                 Each ancestor dim must match a dim in DIMVEC and have lev<=
                  DO L = 1, ANCNDIMS
                        DID = 0
                        DO K = 1, NDIMS
                              IF (ANCACTIVEDIMS(L) .EQ. DIMVEC(K)) THEN
                                    IF (ANCACTIVELEV(L) .LE.
     &                                  NEWLEV(K)) THEN
                                          DID = K
                                    ELSE
                                          SKIPLEVEL = 1
                                    END IF
                              END IF
                        END DO
                        IF (DID .EQ. 0 .AND.
     &                      ANCACTIVELEV(L) .GT. 0) THEN
                              SKIPLEVEL = 1
                        END IF
                        IF (SKIPLEVEL .EQ. 1) GOTO 300
                  END DO

C                 CI is an ancestor. Apply hat function from CI to Y(KSTART..KEND)
                  INDEX = SUBGRIDADDR(CI)

C                 Build repvec for ancestor
                  DO L = 1, ANCNDIMS
                        LVAL = ANCACTIVELEV(L)
                        IF (LVAL .EQ. 0) THEN
                              NPTS_DIM(L) = 1
                        ELSE IF (LVAL .LT. 3) THEN
                              NPTS_DIM(L) = 2
                        ELSE
                              NPTS_DIM(L) = 2**(LVAL-1)
                        END IF
                  END DO
                  REPVEC(1) = 1
                  DO L = 2, ANCNDIMS
                        REPVEC(L) = REPVEC(L-1) * NPTS_DIM(L-1)
                  END DO

                  DO K = KSTART, KEND
                        TEMP = 1.0D+00
                        L = 1

 10                     IF (L .LE. ANCNDIMS) THEN
                              LVAL = ANCACTIVELEV(L)
                              YT   = Y(K, ANCACTIVEDIMS(L))

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
     &                                            2.0D+00*(0.5D+00-YT)
                                          ELSE
                                                TEMP = TEMP *
     &                                            2.0D+00*(YT-0.5D+00)
                                          END IF
                                          INDEX2(L) = XP
                                    END IF

                              ELSE
                                    IF (YT .EQ. 1.0D+00) THEN
                                          TEMP = 0.0D+00
                                    ELSE
                                          SCALE = DBLE(2**LVAL)
                                          XP = INT(YT*SCALE/2.0D+00)
                                          TEMP = TEMP * (1.0D+00
     &                                         - SCALE * DABS(YT
     &                                         - DBLE(XP*2+1)/SCALE))
                                          INDEX2(L) = XP
                                    END IF
                              END IF

                              IF (TEMP .EQ. 0.0D+00) GOTO 100

                              L = L + 1
                              GOTO 10
                        END IF

                        INDEX3 = INDEX + INDEX2(1)
                        DO L = 2, ANCNDIMS
                              INDEX3 = INDEX3 + REPVEC(L)*INDEX2(L)
                        END DO
                        IP(K) = IP(K) + TEMP * Z(INDEX3)

 100                    CONTINUE
                  END DO

 300              CONTINUE
            END DO

            KSTART = KEND + 1

      END DO

      RETURN
      END
