C     *******************************************************************
C
C       SPQUADW_CC_SP - CC quadrature weights, sparse index structure.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPQUADW_CC_SP(INDICESNDIIMS, NSUBGRIDS,
     &    INDICESLEVS, INDICESADDR, NADDR,
     &    SUBGRIDPOINTS, W, NW)
C     *******************************************************************
C
C       SPQUADW_CC_SP returns the CC quadrature weights using the sparse
C       index structure.
C
C       Weight per node:
C            lev = 0,1,2  ->  1/4 contribution per dim
C            lev >= 3     ->  1/2^lev contribution per dim
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
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER NSUBGRIDS, NADDR, NW
      INTEGER INDICESNDIIMS(NSUBGRIDS)
      INTEGER INDICESLEVS(NADDR)
      INTEGER INDICESADDR(NSUBGRIDS)
      INTEGER SUBGRIDPOINTS(NSUBGRIDS)
      DOUBLE PRECISION W(NW)

      INTEGER CI, DID, NDIMS, ADDR, LVAL, NPTS
      INTEGER K, INDEX
      DOUBLE PRECISION WVAL

      INDEX = 1

      DO CI = 1, NSUBGRIDS

            NDIMS = INDICESNDIIMS(CI)
            ADDR  = INDICESADDR(CI)
            NPTS  = SUBGRIDPOINTS(CI)

            WVAL = 1.0D+00
            DO DID = 1, NDIMS
                  LVAL = INDICESLEVS(ADDR + DID - 1)
                  IF (LVAL .LT. 3) THEN
                        WVAL = WVAL * 0.25D+00
                  ELSE
                        WVAL = WVAL / DBLE(2**LVAL)
                  END IF
            END DO

            DO K = INDEX, INDEX + NPTS - 1
                  W(K) = WVAL
            END DO

            INDEX = INDEX + NPTS

      END DO

      RETURN
      END
