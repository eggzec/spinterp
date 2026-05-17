C     *******************************************************************
C
C       SPGET_NPOINTS_NB - Count grid points per subgrid, NoBoundary grid.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPGET_NPOINTS_NB(LEVELSEQ, NLEVELS, D,
     &                             TOTALPOINTS, NPOINTS)
C     *******************************************************************
C
C       SPGET_NPOINTS_NB counts the number of grid points in each subgrid
C       row of LEVELSEQ for the NoBoundary grid, and returns the total.
C
C       Point count per dimension (all levels):
C            lev >= 0 ->  2^lev
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

      INTEGER KL, L, NP

      TOTALPOINTS = 0

      DO KL = 1, NLEVELS
            NP = 1
            DO L = 1, D
                  NP = NP * 2**LEVELSEQ(KL, L)
            END DO
            NPOINTS(KL) = NP
            TOTALPOINTS = TOTALPOINTS + NP
      END DO

      RETURN
      END
