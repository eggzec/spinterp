C     *******************************************************************
C
C       GET_CHEB_NODES - Chebyshev nodes for barycentric interpolation.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE GET_CHEB_NODES(ALLNX, D, X, NX)
C     *******************************************************************
C
C       GET_CHEB_NODES builds a flat array of Chebyshev-Lobatto nodes
C       for each of the D dimensions with ALLNX(k) nodes in dim k.
C
C       Node formula for dimension k with n = ALLNX(k) nodes:
C            x_i = 0.5 - cos(i*pi/(n-1)) / 2,  i = 0,...,n-1
C
C       Parameters:
C
C         Input,  INTEGER ALLNX(D)        - node counts per dimension
C         Input,  INTEGER D               - number of dimensions
C         Output, DOUBLE PRECISION X(NX)  - concatenated node array
C         Input,  INTEGER NX              - total nodes = sum(ALLNX)
C
C     *******************************************************************

      IMPLICIT NONE

      DOUBLE PRECISION PI
      PARAMETER (PI = 3.14159265358979323846D+00)

      INTEGER D, NX
      INTEGER ALLNX(D)
      DOUBLE PRECISION X(NX)

      INTEGER K, I, AID, N

      AID = 0
      DO K = 1, D
            N = ALLNX(K)
            DO I = 0, N-1
                  X(AID+I+1) = 0.5D+00
     &                  - DCOS(DBLE(I)*PI/DBLE(N-1)) / 2.0D+00
            END DO
            AID = AID + N
      END DO

      RETURN
      END
