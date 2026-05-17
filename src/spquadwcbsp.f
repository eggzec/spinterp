C     *******************************************************************
C
C       SPQUADW_CB_SP - CB quadrature weights, sparse index structure.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPQUADW_CB_SP(INDICESNDIIMS, NSUBGRIDS,
     &    INDICESLEVS, INDICESADDR, NADDR,
     &    SUBGRIDPOINTS, W, NW, W1D, NW1D, STARTID)
C     *******************************************************************
C
C       SPQUADW_CB_SP returns the Chebyshev quadrature weights using
C       the sparse index structure. Levels are sorted descending so
C       innermost loop iterates over the largest dimension.
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
      INTEGER ACTIVELEV(50), ORDER(50), TMP, I, J
      DOUBLE PRECISION WVAL

      INDEX = 1

      DO CI = 1, NSUBGRIDS

            NDIMS = INDICESNDIIMS(CI)
            ADDR  = INDICESADDR(CI)
            NPTS  = SUBGRIDPOINTS(CI)

C           Read active levels
            DO DID = 1, NDIMS
                  ACTIVELEV(DID) = INDICESLEVS(ADDR + DID - 1)
                  ORDER(DID) = DID
            END DO

C           Sort by level descending
            DO I = 2, NDIMS
                  TMP = ORDER(I)
                  J = I - 1
 10               IF (J .GE. 1 .AND.
     &                ACTIVELEV(ORDER(J)) .LT. ACTIVELEV(TMP)) THEN
                        ORDER(J+1) = ORDER(J)
                        J = J - 1
                        GOTO 10
                  END IF
                  ORDER(J+1) = TMP
            END DO

            DO DID = 1, NDIMS
                  LVAL = ACTIVELEV(ORDER(DID))
                  IF (LVAL .EQ. 0) THEN
                        NPTS_DIM(DID) = 1
                  ELSE IF (LVAL .LE. 2) THEN
                        NPTS_DIM(DID) = 2
                  ELSE
                        NPTS_DIM(DID) = 2**(LVAL-1)
                  END IF
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
