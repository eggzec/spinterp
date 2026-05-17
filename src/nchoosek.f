C     *******************************************************************
C
C       NCHOOSEK - Binomial coefficient C(N, K).
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      INTEGER FUNCTION NCHOOSEK(N, K)
C     *******************************************************************
C
C       NCHOOSEK computes the binomial coefficient C(N, K).
C
C       Parameters:
C
C         Input, INTEGER N - total elements
C         Input, INTEGER K - chosen elements
C
C         Return, INTEGER NCHOOSEK - binomial coefficient
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER N, K, I, RESULT

      RESULT = 1
      DO I = 0, K - 1
            RESULT = RESULT * (N - I) / (I + 1)
      END DO
      NCHOOSEK = RESULT

      RETURN
      END
