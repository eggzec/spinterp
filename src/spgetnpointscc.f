C     *******************************************************************
C
C       SPGET_NPOINTS_CC - Count grid points per subgrid, CC grid.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPGET_NPOINTS_CC(LEVELSEQ, NLEVELS, D,
     &                             TOTALPOINTS, NPOINTS)
C     *******************************************************************
C
C       SPGET_NPOINTS_CC counts the number of grid points in each subgrid
C       row of LEVELSEQ for the Clenshaw-Curtis grid, and returns the
C       total.
C
C       Point count per dimension:
C            lev = 0      ->  1
C            lev = 1,2   ->  2
C            lev >= 3    ->  2^(lev-1)
C
C       Parameters:
C
C         Input,  INTEGER LEVELSEQ(NLEVELS,D) - multi-index array
C         Input,  INTEGER NLEVELS             - rows in LEVELSEQ
C         Input,  INTEGER D                   - dimension
C         Output, INTEGER TOTALPOINTS         - total number of points
C         Output, INTEGER NPOINTS(NLEVELS)    - points per subgrid
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER NLEVELS, D
      INTEGER LEVELSEQ(NLEVELS, D)
      INTEGER TOTALPOINTS
      INTEGER NPOINTS(NLEVELS)

      INTEGER KL, L, LVAL, NP

      TOTALPOINTS = 0

      DO KL = 1, NLEVELS
            NP = 1
            DO L = 1, D
                  LVAL = LEVELSEQ(KL, L)
                  IF (LVAL .EQ. 0) THEN
                        NP = NP * 1
                  ELSE IF (LVAL .LT. 3) THEN
                        NP = NP * 2
                  ELSE
                        NP = NP * 2**(LVAL-1)
                  END IF
            END DO
            NPOINTS(KL) = NP
            TOTALPOINTS = TOTALPOINTS + NP
      END DO

      RETURN
      END
