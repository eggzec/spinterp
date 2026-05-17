C     *******************************************************************
C
C       SPQUADW_GP_SP - GP quadrature weights, sparse index structure.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPQUADW_GP_SP(INDICESNDIIMS, NSUBGRIDS,
     &    INDICESLEVS, INDICESADDR, NADDR,
     &    SUBGRIDPOINTS, W, NW, W1D, NW1D, STARTID)
C     *******************************************************************
C
C       SPQUADW_GP_SP returns the GP quadrature weights using the sparse
C       index structure. npoints per dim = 2^lev (not 2^(lev-1)).
C
C       Parameters:
C
C         Input,  INTEGER INDICESNDIIMS(NSUBGRIDS)   - ndims per subgrid
C         Input,  INTEGER NSUBGRIDS                  - number of subgrids
C         Input,  INTEGER INDICESLEVS(NADDR)          - packed levels
C         Input,  INTEGER INDICESADDR(NSUBGRIDS)      - start addr per subgrid
C         Input,  INTEGER NADDR                       - length of packed arrays
C         Input,  INTEGER SUBGRIDPOINTS(NSUBGRIDS)    - points per subgrid
C         Output, DOUBLE PRECISION W(NW)              - quadrature weights
C         Input,  INTEGER NW                          - length of W
C         Input,  DOUBLE PRECISION W1D(NW1D)          - 1-D weight table
C         Input,  INTEGER NW1D                        - length of W1D
C         Input,  INTEGER STARTID(*)                  - start index per level
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER NSUBGRIDS, NADDR, NW, NW1D
      INTEGER INDICESNDIIMS(NSUBGRIDS)
      INTEGER INDICESLEVS(NADDR)
      INTEGER INDICESADDR(NSUBGRIDS)
      INTEGER SUBGRIDPOINTS(NSUBGRIDS)
      DOUBLE PRECISION W(NW)
      DOUBLE PRECISION W1D(NW1D)
      INTEGER STARTID(*)

      INTEGER CI, DID, NDIMS, ADDR, LVAL, NPTS
      INTEGER K, INDEX
      INTEGER NPTS_DIM(50), WIDSTART(50), WIDEND(50), WID(50), WID1
      DOUBLE PRECISION WVAL

      INDEX = 1

      DO CI = 1, NSUBGRIDS

            NDIMS = INDICESNDIIMS(CI)
            ADDR  = INDICESADDR(CI)
            NPTS  = SUBGRIDPOINTS(CI)

            DO DID = 1, NDIMS
                  LVAL = INDICESLEVS(ADDR + DID - 1)
                  NPTS_DIM(DID) = 2**LVAL
                  WIDSTART(DID) = STARTID(LVAL+1)
                  WID(DID) = WIDSTART(DID)
                  WIDEND(DID) = WID(DID) + NPTS_DIM(DID) - 1
            END DO

            IF (NDIMS .EQ. 0) THEN
                  DO K = INDEX, INDEX + NPTS - 1
                        W(K) = 1.0D+00
                  END DO
                  INDEX = INDEX + NPTS
                  GOTO 200
            END IF

            WID1 = WID(1)
            DO K = INDEX, INDEX + NPTS - 1
                  WVAL = W1D(WID1)
                  DO DID = 2, NDIMS
                        WVAL = WVAL * W1D(WID(DID))
                  END DO
                  W(K) = WVAL

                  IF (WID1 .LT. WIDEND(1)) THEN
                        WID1 = WID1 + 1
                  ELSE
                        WID1 = WIDSTART(1)
                        DO DID = 2, NDIMS
                              IF (WID(DID) .LT. WIDEND(DID)) THEN
                                    WID(DID) = WID(DID) + 1
                                    GOTO 100
                              ELSE
                                    WID(DID) = WIDSTART(DID)
                              END IF
                        END DO
                  END IF
 100              CONTINUE
            END DO

            INDEX = INDEX + NPTS

 200        CONTINUE
      END DO

      RETURN
      END
