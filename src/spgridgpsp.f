C     *******************************************************************
C
C       SPGRID_GP_SP - GP grid points, sparse index structure.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPGRID_GP_SP(INDICESNDIIMS, NSUBGRIDS,
     &    INDICESDIMS, INDICESLEVS, INDICESADDR, NADDR,
     &    SUBGRIDPOINTS, D, X, TOTALPOINTS, FROMINDEX, TOINDEX)
C     *******************************************************************
C
C       SPGRID_GP_SP computes Gauss-Patterson sparse grid points using
C       the sparse index structure. For each lev, NX=2^(lev+1)-1, calls
C       GP_ABSC and takes every other node (odd 1-based indices).
C
C       Parameters:
C
C         Input,  INTEGER INDICESNDIIMS(NSUBGRIDS)   - ndims per subgrid
C         Input,  INTEGER NSUBGRIDS                  - number of subgrids
C         Input,  INTEGER INDICESDIMS(NADDR)          - packed dim indices
C         Input,  INTEGER INDICESLEVS(NADDR)          - packed levels
C         Input,  INTEGER INDICESADDR(NSUBGRIDS)      - start addr per subgrid
C         Input,  INTEGER NADDR                       - length of packed arrays
C         Input,  INTEGER SUBGRIDPOINTS(NSUBGRIDS)    - points per subgrid
C         Input,  INTEGER D                           - dimension
C         Output, DOUBLE PRECISION X(TOTALPOINTS,D)  - grid coordinates
C         Input,  INTEGER TOTALPOINTS                 - total points
C         Input,  INTEGER FROMINDEX                   - first subgrid (1-based)
C         Input,  INTEGER TOINDEX                     - last subgrid (1-based)
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER NSUBGRIDS, NADDR, D, TOTALPOINTS, FROMINDEX, TOINDEX
      INTEGER INDICESNDIIMS(NSUBGRIDS)
      INTEGER INDICESDIMS(NADDR)
      INTEGER INDICESLEVS(NADDR)
      INTEGER INDICESADDR(NSUBGRIDS)
      INTEGER SUBGRIDPOINTS(NSUBGRIDS)
      DOUBLE PRECISION X(TOTALPOINTS, D)

      INTEGER CI, DID, NDIMS, ADDR, LVAL, NPTS
      INTEGER ACTIVEDIMS(50), ACTIVELEV(50)
      INTEGER NPTS_DIM(50), REP(50)
      INTEGER I, K, J, IDX, NX
      DOUBLE PRECISION ABSC(127)

      DO K = 1, D
            DO I = 1, TOTALPOINTS
                  X(I, K) = 0.5D+00
            END DO
      END DO

      IDX = 1
      DO CI = 1, FROMINDEX - 1
            IDX = IDX + SUBGRIDPOINTS(CI)
      END DO

      DO CI = FROMINDEX, TOINDEX

            NDIMS = INDICESNDIIMS(CI)
            ADDR  = INDICESADDR(CI)
            NPTS  = SUBGRIDPOINTS(CI)

            IF (NDIMS .EQ. 0) THEN
                  IDX = IDX + NPTS
                  GOTO 200
            END IF

            DO DID = 1, NDIMS
                  ACTIVEDIMS(DID) = INDICESDIMS(ADDR + DID - 1)
                  ACTIVELEV(DID)  = INDICESLEVS(ADDR + DID - 1)
            END DO

            DO DID = 1, NDIMS
                  LVAL = ACTIVELEV(DID)
                  NPTS_DIM(DID) = 2**LVAL
            END DO

            REP(1) = 1
            DO DID = 2, NDIMS
                  REP(DID) = REP(DID-1) * NPTS_DIM(DID-1)
            END DO

            DO I = 0, NPTS - 1
                  DO DID = 1, NDIMS
                        LVAL = ACTIVELEV(DID)
                        NX   = 2**(LVAL+1) - 1
                        CALL GP_ABSC(LVAL, ABSC, NX)
                        J = MOD(I / REP(DID), NPTS_DIM(DID))
C                       New GP nodes occupy odd 1-based positions (1,3,5,...).
C                       Skip by 2: index = J*2+1 (1-based).
                        X(IDX+I, ACTIVEDIMS(DID)) = ABSC(J*2 + 1)
                  END DO
            END DO

            IDX = IDX + NPTS

 200        CONTINUE
      END DO

      RETURN
      END
