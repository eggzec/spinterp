C     *******************************************************************
C
C       SPCMPVALS_CBGP_SP - CB/GP hierarchical surpluses, sparse struct.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPCMPVALS_CBGP_SP(D, Z, NZ, Y, NY,
     &    INDICESNDIIMS, NSUBGRIDS,
     &    INDICESDIMS, INDICESLEVS, INDICESADDR, NADDR,
     &    BACKWARDNEIGHBORS, FORWARDNEIGHBORS, NFWD,
     &    SUBGRIDPOINTS, SUBGRIDADDR,
     &    FROMINDEX, TOINDEX, ISGP, IP)
C     *******************************************************************
C
C       SPCMPVALS_CBGP_SP computes hierarchical surplus increments for
C       new subgrids using Chebyshev (ISGP=0) or GP (ISGP=1) barycentric
C       interpolation and the sparse index structure.
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
C         Input,  INTEGER ISGP                        - 0=Cheb, 1=GP
C         Output, DOUBLE PRECISION IP(NY)             - interpolant at Y
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER D, NZ, NY, NSUBGRIDS, NADDR, NFWD
      INTEGER FROMINDEX, TOINDEX, ISGP
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

      INTEGER NKL, CI, DID, K, L
      INTEGER NDIMS, ADDR
      INTEGER NEWLEV(50), DIMVEC(50)
      INTEGER LVAL, NPTS
      INTEGER ANCACTIVEDIMS(50), ANCACTIVELEV(50)
      INTEGER ALLNX(50), DIMS(50), XTOT
      INTEGER SKIPLEVEL
      INTEGER ANCNDIMS, ANCADDR, ANCNPTS
      INTEGER KSTART, KEND, INDEX
      INTEGER ORDERARR(50), TMP, I, J
      DOUBLE PRECISION XBUF(16384), WBUF(16384)
      DOUBLE PRECISION IPTEMP(65536)

      DO K = 1, NY
            IP(K) = 0.0D+00
      END DO

      KSTART = 1
      DO CI = 1, FROMINDEX - 1
            KSTART = KSTART + SUBGRIDPOINTS(CI)
      END DO

      DO NKL = FROMINDEX, TOINDEX

            NDIMS = INDICESNDIIMS(NKL)
            ADDR  = INDICESADDR(NKL)
            NPTS  = SUBGRIDPOINTS(NKL)
            KEND  = KSTART + NPTS - 1

            DO DID = 1, NDIMS
                  DIMVEC(DID) = INDICESDIMS(ADDR + DID - 1)
                  NEWLEV(DID) = INDICESLEVS(ADDR + DID - 1)
            END DO

C           Walk all ancestor subgrids
            DO CI = 1, NSUBGRIDS
                  IF (CI .GE. FROMINDEX) GOTO 300

                  ANCNDIMS = INDICESNDIIMS(CI)
                  ANCADDR  = INDICESADDR(CI)
                  ANCNPTS  = SUBGRIDPOINTS(CI)
                  INDEX    = SUBGRIDADDR(CI)

                  SKIPLEVEL = 0
                  DO L = 1, ANCNDIMS
                        ANCACTIVEDIMS(L) = INDICESDIMS(ANCADDR + L - 1)
                        ANCACTIVELEV(L)  = INDICESLEVS(ANCADDR + L - 1)
                  END DO

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

                  IF (ANCNPTS .EQ. 1) THEN
C                       Constant subgrid: add Z(INDEX) to all new points
                        DO K = KSTART, KEND
                              IP(K) = IP(K) + Z(INDEX)
                        END DO
                        GOTO 300
                  END IF

C                 Build ALLNX/DIMS sorted by level descending
                  DO L = 1, ANCNDIMS
                        ORDERARR(L) = L
                  END DO
                  DO I = 2, ANCNDIMS
                        TMP = ORDERARR(I)
                        J = I - 1
 10                     IF (J .GE. 1 .AND.
     &                      ANCACTIVELEV(ORDERARR(J)) .LT.
     &                      ANCACTIVELEV(TMP)) THEN
                              ORDERARR(J+1) = ORDERARR(J)
                              J = J - 1
                              GOTO 10
                        END IF
                        ORDERARR(J+1) = TMP
                  END DO

                  XTOT = 0
                  DO I = 1, ANCNDIMS
                        L = ORDERARR(I)
                        LVAL = ANCACTIVELEV(L)
                        IF (ISGP .EQ. 0) THEN
                              ALLNX(I) = 2**LVAL + 1
                        ELSE
                              ALLNX(I) = 2**(LVAL+1) - 1
                        END IF
                        DIMS(I) = ANCACTIVEDIMS(L)
                        XTOT = XTOT + ALLNX(I)
                  END DO

                  DO K = 1, KEND - KSTART + 1
                        IPTEMP(K) = 0.0D+00
                  END DO

                  IF (ISGP .EQ. 0) THEN
                        CALL GET_CHEB_NODES(ALLNX, ANCNDIMS,
     &                        XBUF, XTOT)
                        CALL BARY_PD_STEP_CB(Z(INDEX), ANCNPTS,
     &                        ALLNX, DIMS, ANCNDIMS,
     &                        XBUF, XTOT,
     &                        Y(KSTART,1), KEND-KSTART+1, D,
     &                        IPTEMP)
                  ELSE
                        CALL GET_GP_NODES(ALLNX, ANCNDIMS, XBUF, XTOT)
                        CALL GET_GP_BARY_W(ALLNX, ANCNDIMS, WBUF, XTOT)
                        CALL BARY_PD_STEP_GP(Z(INDEX), ANCNPTS,
     &                        ALLNX, DIMS, ANCNDIMS,
     &                        XBUF, XTOT,
     &                        Y(KSTART,1), KEND-KSTART+1, D,
     &                        WBUF, IPTEMP)
                  END IF

                  DO K = KSTART, KEND
                        IP(K) = IP(K) + IPTEMP(K - KSTART + 1)
                  END DO

 300              CONTINUE
            END DO

            KSTART = KEND + 1

      END DO

      RETURN
      END
