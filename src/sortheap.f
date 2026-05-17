C     *******************************************************************
C
C       SORT_HEAP - Max-heap sift-up after inserting at position NA.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SORT_HEAP(A, NA, G, NG)
C     *******************************************************************
C
C       SORT_HEAP sifts up the element at position NA in the heap A,
C       restoring the max-heap property with respect to priorities G.
C
C       The heap stores indices into G; A(k) is an index such that
C       G(A(k)) is the priority at heap position k.
C
C       Parameters:
C
C         In/out, INTEGER A(NA)          - heap indices (1-based)
C         Input,  INTEGER NA             - current heap size (element added)
C         Input,  DOUBLE PRECISION G(NG) - priority values
C         Input,  INTEGER NG             - length of G
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER NA, NG
      INTEGER A(NA)
      DOUBLE PRECISION G(NG)

      INTEGER K, PREV, TMP

      K = NA

 10   IF (K .GT. 1) THEN
            PREV = K / 2
            IF (G(A(PREV)) .LT. G(A(K))) THEN
                  TMP    = A(PREV)
                  A(PREV) = A(K)
                  A(K)   = TMP
                  K = PREV
                  GOTO 10
            END IF
      END IF

      RETURN
      END
