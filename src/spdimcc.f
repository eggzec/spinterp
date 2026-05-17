C     *******************************************************************
C
C       SPDIM_CC - Count grid points, Clenshaw-Curtis/Chebyshev grid.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPDIM_CC(LEVELSEQ, NLEVELS, D, TOTALPOINTS)
C     *******************************************************************
C
C       SPDIM_CC counts the total grid points for a Clenshaw-Curtis
C       or Chebyshev sparse grid defined by LEVELSEQ.
C
C       Node count per 1-D level:
C            lev = 0        ->  1  (midpoint 0.5)
C            lev = 1 or 2   ->  2  (boundary or first interior pair)
C            lev >= 3       ->  2^(lev-1)  (interior nodes)
C
C       Parameters:
C
C         Input,  INTEGER LEVELSEQ(NLEVELS,D) - multi-index array
C         Input,  INTEGER NLEVELS             - rows in LEVELSEQ
C         Input,  INTEGER D                   - dimension
C         Output, INTEGER TOTALPOINTS         - total grid points
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER NLEVELS, D
      INTEGER LEVELSEQ(NLEVELS, D)
      INTEGER TOTALPOINTS

      INTEGER KL, K, LEV, NPTS

      TOTALPOINTS = 0
      DO KL = 1, NLEVELS
            NPTS = 1
            DO K = 1, D
                  LEV = LEVELSEQ(KL, K)
                  IF (LEV .EQ. 0) THEN
                        CONTINUE
                  ELSE IF (LEV .LE. 2) THEN
                        NPTS = NPTS * 2
                  ELSE
                        NPTS = NPTS * 2**(LEV-1)
                  END IF
            END DO
            TOTALPOINTS = TOTALPOINTS + NPTS
      END DO

      RETURN
      END
