C     *******************************************************************
C
C       SPDIM_M - Count grid points, Maximum / NoBoundary grid.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPDIM_M(LEVELSEQ, NLEVELS, D, TOTALPOINTS, BOUNDARY)
C     *******************************************************************
C
C       SPDIM_M counts total grid points for Maximum or NoBoundary grids.
C
C       Node count per 1-D level:
C            Maximum (BOUNDARY=1):
C               lev = 0  ->  3  (nodes at 0, 0.5, 1)
C               lev >= 1 ->  2^lev  (interior midpoints)
C            NoBoundary (BOUNDARY=0):
C               lev >= 0 ->  2^lev  (interior midpoints)
C
C       Parameters:
C
C         Input,  INTEGER LEVELSEQ(NLEVELS,D) - multi-index array
C         Input,  INTEGER NLEVELS             - rows in LEVELSEQ
C         Input,  INTEGER D                   - dimension
C         Output, INTEGER TOTALPOINTS         - total grid points
C         Input,  INTEGER BOUNDARY            - 1=Maximum, 0=NoBoundary
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER NLEVELS, D, BOUNDARY
      INTEGER LEVELSEQ(NLEVELS, D)
      INTEGER TOTALPOINTS

      INTEGER KL, K, LEV, NPTS

      TOTALPOINTS = 0
      DO KL = 1, NLEVELS
            NPTS = 1
            DO K = 1, D
                  LEV = LEVELSEQ(KL, K)
                  IF (BOUNDARY .EQ. 1 .AND. LEV .EQ. 0) THEN
                        NPTS = NPTS * 3
                  ELSE
                        NPTS = NPTS * 2**LEV
                  END IF
            END DO
            TOTALPOINTS = TOTALPOINTS + NPTS
      END DO

      RETURN
      END
