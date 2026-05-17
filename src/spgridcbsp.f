C     *******************************************************************
C
C       SPGRID_CB_SP - Chebyshev grid points, sparse index structure.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPGRID_CB_SP(INDICESNDIIMS, NSUBGRIDS,
     &    INDICESDIMS, INDICESLEVS, INDICESADDR, NADDR,
     &    SUBGRIDPOINTS, D, X, TOTALPOINTS, FROMINDEX, TOINDEX)
C     *******************************************************************
C
C       SPGRID_CB_SP computes Chebyshev sparse grid points using sparse
C       index structure.
C
C       1-D node formula (new nodes only):
C            lev = 1   ->  {0, 1}
C            lev >= 2  ->  {0.5 - cos((2i-1)*pi/2^lev)/2,
C                           i=1,...,2^(lev-1)}
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

      DOUBLE PRECISION PI
      PARAMETER (PI = 3.14159265358979323846D+00)

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
      INTEGER I, K, J, IDX
      DOUBLE PRECISION COORD

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
                  IF (LVAL .EQ. 0) THEN
                        NPTS_DIM(DID) = 1
                  ELSE IF (LVAL .LE. 2) THEN
                        NPTS_DIM(DID) = 2
                  ELSE
                        NPTS_DIM(DID) = 2**(LVAL-1)
                  END IF
            END DO

            REP(1) = 1
            DO DID = 2, NDIMS
                  REP(DID) = REP(DID-1) * NPTS_DIM(DID-1)
            END DO

            DO I = 0, NPTS - 1
                  DO DID = 1, NDIMS
                        LVAL = ACTIVELEV(DID)
                        J = MOD(I / REP(DID), NPTS_DIM(DID))

                        IF (LVAL .EQ. 0) THEN
                              COORD = 0.5D+00
                        ELSE IF (LVAL .EQ. 1) THEN
                              IF (J .EQ. 0) THEN
                                    COORD = 0.0D+00
                              ELSE
                                    COORD = 1.0D+00
                              END IF
                        ELSE
                              COORD = 0.5D+00 - DCOS(DBLE(2*J+1)*PI
     &                                / DBLE(2**LVAL)) / 2.0D+00
                        END IF

                        X(IDX+I, ACTIVEDIMS(DID)) = COORD
                  END DO
            END DO

            IDX = IDX + NPTS

 200        CONTINUE
      END DO

      RETURN
      END
