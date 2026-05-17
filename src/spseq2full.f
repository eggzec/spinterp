C     *******************************************************************
C
C       SPSEQ2FULL - Convert sparse index structure to full dense levelseq.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPSEQ2FULL(INDICESNDIIMS, NSUBGRIDS,
     &                       INDICESDIMS, INDICESLEVS, INDICESADDR,
     &                       NADDR, D, FULLSEQ)
C     *******************************************************************
C
C       SPSEQ2FULL expands a sparse index representation (packed arrays
C       of active dimension indices and levels) into a full dense level
C       matrix FULLSEQ of size NSUBGRIDS x D.
C
C       For each subgrid k, the active dims and their levels are stored
C       starting at address INDICESADDR(k) in INDICESDIMS / INDICESLEVS.
C       Inactive dimensions default to zero.
C
C       Parameters:
C
C         Input,  INTEGER INDICESNDIIMS(NSUBGRIDS) - ndims per subgrid
C         Input,  INTEGER NSUBGRIDS               - number of subgrids
C         Input,  INTEGER INDICESDIMS(NADDR)       - packed dim indices
C         Input,  INTEGER INDICESLEVS(NADDR)       - packed levels
C         Input,  INTEGER INDICESADDR(NSUBGRIDS)   - start addr per subgrid
C         Input,  INTEGER NADDR                    - length of packed arrays
C         Input,  INTEGER D                        - full dimension
C         Output, INTEGER FULLSEQ(NSUBGRIDS,D)     - full level matrix
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER NSUBGRIDS, NADDR, D
      INTEGER INDICESNDIIMS(NSUBGRIDS)
      INTEGER INDICESDIMS(NADDR)
      INTEGER INDICESLEVS(NADDR)
      INTEGER INDICESADDR(NSUBGRIDS)
      INTEGER FULLSEQ(NSUBGRIDS, D)

      INTEGER K, DID, ADDR, DIM

      DO K = 1, NSUBGRIDS
C           Zero out all dimensions
            DO DID = 1, D
                  FULLSEQ(K, DID) = 0
            END DO

C           Fill in active dimensions
            ADDR = INDICESADDR(K)
            DO DID = 1, INDICESNDIIMS(K)
                  DIM = INDICESDIMS(ADDR + DID - 1)
                  FULLSEQ(K, DIM) = INDICESLEVS(ADDR + DID - 1)
            END DO
      END DO

      RETURN
      END
