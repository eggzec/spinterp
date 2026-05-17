C     *******************************************************************
C
C       POP_HEAP - Pop max from max-heap and return its index.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE POP_HEAP(A, NA, G, NG, INDEX)
C     *******************************************************************
C
C       POP_HEAP removes the max element from the heap (A(1)), replacing
C       it with A(NA), then sifts down to restore the heap property.
C       Returns the popped index in INDEX.
C
C       Caller must set A(NA) to the replacement value before calling,
C       and decrement NA after the call.
C
C       Parameters:
C
C         In/out, INTEGER A(NA)          - heap indices (A(NA) = replacement)
C         Input,  INTEGER NA             - current heap size
C         Input,  DOUBLE PRECISION G(NG) - priority values
C         Input,  INTEGER NG             - length of G
C         Output, INTEGER INDEX          - popped index (old A(1))
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER NA, NG
      INTEGER A(NA)
      DOUBLE PRECISION G(NG)
      INTEGER INDEX

      INTEGER K, LEFT, RIGHT, LARGEST, TMP

      INDEX  = A(1)
      A(1)   = A(NA)

      K = 1

 10   CONTINUE
            LEFT  = 2 * K
            RIGHT = 2 * K + 1

C           Find largest among k, left, right (within NA-1 elements)
            LARGEST = K
            IF (LEFT .LE. NA - 1) THEN
                  IF (G(A(LEFT)) .GT. G(A(LARGEST))) THEN
                        LARGEST = LEFT
                  END IF
            END IF
            IF (RIGHT .LE. NA - 1) THEN
                  IF (G(A(RIGHT)) .GT. G(A(LARGEST))) THEN
                        LARGEST = RIGHT
                  END IF
            END IF

            IF (LARGEST .NE. K) THEN
                  TMP        = A(K)
                  A(K)       = A(LARGEST)
                  A(LARGEST) = TMP
                  K = LARGEST
                  GOTO 10
            END IF

      RETURN
      END
