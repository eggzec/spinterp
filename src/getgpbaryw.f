C     *******************************************************************
C
C       GET_GP_BARY_W - Gauss-Patterson barycentric weights array.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE GET_GP_BARY_W(ALLNX, D, W, NW)
C     *******************************************************************
C
C       GET_GP_BARY_W builds a flat array of GP barycentric weights for
C       each of the D dimensions.  ALLNX(k) = 2^(lev+1)-1 for dim k.
C
C       Parameters:
C
C         Input,  INTEGER ALLNX(D)        - full node counts per dim
C         Input,  INTEGER D               - number of dimensions
C         Output, DOUBLE PRECISION W(NW)  - concatenated weight array
C         Input,  INTEGER NW              - total weights = sum(ALLNX)
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER D, NW
      INTEGER ALLNX(D)
      DOUBLE PRECISION W(NW)

      INTEGER K, AID, N, LEV

      AID = 0
      DO K = 1, D
            N = ALLNX(K)
            LEV = 0
 10         IF (2**(LEV+1) - 1 .LT. N) THEN
                  LEV = LEV + 1
                  GOTO 10
            END IF
            CALL GP_BARY_W(LEV, W(AID+1), N)
            AID = AID + N
      END DO

      RETURN
      END
