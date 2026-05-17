C     *******************************************************************
C
C       SPNLEVELS - Count multi-index tuples summing to N in D dims.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPNLEVELS(N, D, NLEVELS)
C     *******************************************************************
C
C       SPNLEVELS computes C(N+D-1, D-1): the number of multi-index
C       tuples (l_1,...,l_D) with l_1 + ... + l_D = N, l_i >= 0.
C
C       Parameters:
C
C         Input,  INTEGER N       - level
C         Input,  INTEGER D       - dimension
C         Output, INTEGER NLEVELS - number of multi-indices
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER N, D, NLEVELS
      INTEGER NCHOOSEK

      NLEVELS = NCHOOSEK(N + D - 1, D - 1)

      RETURN
      END
