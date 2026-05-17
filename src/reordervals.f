C     *******************************************************************
C
C       REORDER_VALS - Reorder hierarchical surplus array for CB grid.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE REORDER_VALS(Z, NZ, LEVELSEQ, NLEVELS, D)
C     *******************************************************************
C
C       REORDER_VALS permutes the hierarchical surplus array Z so that
C       dimensions are sorted by level descending within each subgrid.
C       This matches the layout expected by BARY_PD_STEP_CB.
C
C       Uses a temporary buffer ZTMP of size MAXA=16384.
C       Permutation is implemented by building old->new index mappings
C       and copying in and out of the temp buffer.
C
C       Parameters:
C
C         In/out, DOUBLE PRECISION Z(NZ)       - surpluses to reorder
C         Input,  INTEGER NZ                   - length of Z
C         Input,  INTEGER LEVELSEQ(NLEVELS,D)  - multi-index set
C         Input,  INTEGER NLEVELS              - rows in LEVELSEQ
C         Input,  INTEGER D                    - dimension
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER MAXA
      PARAMETER (MAXA = 16384)

      INTEGER NZ, NLEVELS, D
      DOUBLE PRECISION Z(NZ)
      INTEGER LEVELSEQ(NLEVELS, D)

      INTEGER KL, K, L, LVAL, NPTS
      INTEGER NPTS_DIM(50), ORDER(50), ORDERARR(50), LEVEL(50)
      INTEGER NEWPTS_DIM(50)
      INTEGER STRIDE_OLD(50), STRIDE_NEW(50)
      INTEGER FLAT_OLD, FLAT_NEW
      INTEGER INDEX
      INTEGER COORD(50), COORD_PERM(50)
      INTEGER DIM_PERM
      DOUBLE PRECISION ZTMP(MAXA)
      INTEGER I, J, TMP

C     Working index into Z
      INDEX = 1

      DO KL = 1, NLEVELS

C           Compute npts per dim and total
            NPTS = 1
            DO L = 1, D
                  LVAL = LEVELSEQ(KL, L)
                  LEVEL(L) = LVAL
                  IF (LVAL .EQ. 0) THEN
                        NPTS_DIM(L) = 1
                  ELSE IF (LVAL .LE. 2) THEN
                        NPTS_DIM(L) = 2
                  ELSE
                        NPTS_DIM(L) = 2**(LVAL-1)
                  END IF
                  NPTS = NPTS * NPTS_DIM(L)
            END DO

            IF (NPTS .EQ. 1) THEN
                  INDEX = INDEX + 1
                  GOTO 100
            END IF

C           Build sort order: insertion sort by level descending
            DO L = 1, D
                  ORDER(L) = L
            END DO
            DO L = 2, D
                  TMP = ORDER(L)
                  J = L - 1
 10               IF (J .GE. 1 .AND.
     &                LEVEL(ORDER(J)) .LT. LEVEL(TMP)) THEN
                        ORDER(J+1) = ORDER(J)
                        J = J - 1
                        GOTO 10
                  END IF
                  ORDER(J+1) = TMP
            END DO

C           Check if order is identity - if so, skip
            DO L = 1, D
                  IF (ORDER(L) .NE. L) GOTO 20
            END DO
            INDEX = INDEX + NPTS
            GOTO 100

 20         CONTINUE

C           Build inverse permutation: ORDERARR(dim_new_pos) = dim_old_pos
C           ORDER(l) = which old dim goes to new position l
C           So new dim l has points from old dim ORDER(l)
            DO L = 1, D
                  NEWPTS_DIM(L) = NPTS_DIM(ORDER(L))
            END DO

C           Compute strides for old layout (Fortran col-major over dims 1..D)
            STRIDE_OLD(1) = 1
            DO L = 2, D
                  STRIDE_OLD(L) = STRIDE_OLD(L-1) * NPTS_DIM(L-1)
            END DO

C           Compute strides for new layout (after permutation)
            STRIDE_NEW(1) = 1
            DO L = 2, D
                  STRIDE_NEW(L) = STRIDE_NEW(L-1) * NEWPTS_DIM(L-1)
            END DO

C           Copy permuted data into ZTMP
C           For each flat new index, compute old flat index
            DO I = 0, NPTS - 1
C                 Decode new flat index into coordinates
                  K = I
                  DO L = D, 1, -1
                        COORD(L) = K / STRIDE_NEW(L)
                        K = K - COORD(L) * STRIDE_NEW(L)
                  END DO

C                 COORD(L) is coordinate in new-dim L = old dim ORDER(L)
C                 Remap to old coordinates: old dim ORDER(L) gets COORD(L)
                  DO L = 1, D
                        COORD_PERM(ORDER(L)) = COORD(L)
                  END DO

C                 Compute old flat index
                  FLAT_OLD = 0
                  DO L = 1, D
                        FLAT_OLD = FLAT_OLD +
     &                            COORD_PERM(L)*STRIDE_OLD(L)
                  END DO

                  ZTMP(I+1) = Z(INDEX + FLAT_OLD)
            END DO

C           Copy back
            DO I = 0, NPTS - 1
                  Z(INDEX + I) = ZTMP(I+1)
            END DO

            INDEX = INDEX + NPTS

 100        CONTINUE
      END DO

      RETURN
      END
