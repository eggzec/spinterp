C     *******************************************************************
C
C       SPGETSEQ_SP - Build sparse multi-index structure for sparse grid.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPGETSEQ_SP(N, D, GRIDCODE,
     &    INDICESNDIIMS, INDICESDIMS, INDICESLEVS, INDICESADDR,
     &    BACKWARDNEIGHBORS, FORWARDNEIGHBORS,
     &    SUBGRIDPOINTS, SUBGRIDADDR, ACTIVEINDICES,
     &    NSUBGRIDS, MAXADDR, CURRENTINDEX)
C     *******************************************************************
C
C       SPGETSEQ_SP builds the complete sparse multi-index structure
C       for a sparse grid of level N in D dimensions.
C
C       The structure stores multi-indices (l_1,...,l_D) with
C       l_1 + ... + l_D <= N using a compact representation:
C       only non-zero dimensions are stored per subgrid.
C
C       GRIDCODE:
C            0 = Clenshaw-Curtis or Chebyshev
C            2 = NoBoundary or Gauss-Patterson
C
C       Parameters:
C
C         Input,  INTEGER N              - max level
C         Input,  INTEGER D              - dimension
C         Input,  INTEGER GRIDCODE       - grid type code
C         Output, INTEGER INDICESNDIIMS(NSUBGRIDS) - active dims per subgrid
C         Output, INTEGER INDICESDIMS(MAXADDR)     - packed dim indices (1-based)
C         Output, INTEGER INDICESLEVS(MAXADDR)     - packed levels
C         Output, INTEGER INDICESADDR(NSUBGRIDS)   - start addr per subgrid
C         Output, INTEGER BACKWARDNEIGHBORS(MAXADDR) - backward neighbor idx
C         Output, INTEGER FORWARDNEIGHBORS(NSUBGRIDS*D) - forward neighbors
C         Output, INTEGER SUBGRIDPOINTS(NSUBGRIDS) - point count per subgrid
C         Output, INTEGER SUBGRIDADDR(NSUBGRIDS)   - start addr in point array
C         Output, INTEGER ACTIVEINDICES(NSUBGRIDS) - 1=active, 0=passive
C         Input,  INTEGER NSUBGRIDS      - C(N+D, D) total subgrids
C         Input,  INTEGER MAXADDR        - size of packed arrays
C         Output, INTEGER CURRENTINDEX   - next index to process
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER N, D, GRIDCODE, NSUBGRIDS, MAXADDR
      INTEGER INDICESNDIIMS(NSUBGRIDS)
      INTEGER INDICESDIMS(MAXADDR)
      INTEGER INDICESLEVS(MAXADDR)
      INTEGER INDICESADDR(NSUBGRIDS)
      INTEGER BACKWARDNEIGHBORS(MAXADDR)
      INTEGER FORWARDNEIGHBORS(NSUBGRIDS * D)
      INTEGER SUBGRIDPOINTS(NSUBGRIDS)
      INTEGER SUBGRIDADDR(NSUBGRIDS)
      INTEGER ACTIVEINDICES(NSUBGRIDS)
      INTEGER CURRENTINDEX

      INTEGER I, J, K, L, LK, KK
      INTEGER NI, ADDR, NBID, NBAD
      INTEGER DIM, LEV, LVAL, INSERT
      INTEGER NNEWDIMS, NPOINTS, NNEWLEVELS
      INTEGER NCHOOSEK
      INTEGER BN(50), BDIM(50), BID(50), LEVEL(50)
      INTEGER ISADMISSIBLE, BACKOFNEW
      INTEGER LOOPK, LOOPNL, SUBNL

C     Initialize first subgrid: all-zero multi-index (0,...,0)
      DO I = 1, NSUBGRIDS
            INDICESNDIIMS(I) = 0
            INDICESADDR(I) = 0
            SUBGRIDPOINTS(I) = 0
            SUBGRIDADDR(I) = 0
            ACTIVEINDICES(I) = 0
      END DO
      DO I = 1, NSUBGRIDS * D
            FORWARDNEIGHBORS(I) = 0
      END DO
      DO I = 1, MAXADDR
            INDICESDIMS(I) = 0
            INDICESLEVS(I) = 0
            BACKWARDNEIGHBORS(I) = 0
      END DO

C     First subgrid: index 1, no active dims (all levels = 0)
      INDICESNDIIMS(1) = 0
      INDICESADDR(1) = 1
      IF (GRIDCODE .EQ. 1) THEN
C           Maximum grid: 3^D points for level-0 subgrid
            NPOINTS = 1
            DO I = 1, D
                  NPOINTS = NPOINTS * 3
            END DO
      ELSE
            NPOINTS = 1
      END IF
      SUBGRIDPOINTS(1) = NPOINTS
      SUBGRIDADDR(1) = 1

      NI = 1
      ADDR = 1
      CURRENTINDEX = 1

C     Loop over levels 0..N-1, generating new subgrids
      DO LK = 0, N-1

            NNEWLEVELS = NCHOOSEK(D + LK - 1, D - 1)

            DO LOOPK = 1, NNEWLEVELS

                  ACTIVEINDICES(CURRENTINDEX) = 0

C                 Read backward neighbors of current index
                  NBID = INDICESNDIIMS(CURRENTINDEX)
                  NBAD = INDICESADDR(CURRENTINDEX)

                  DO J = 1, NBID
                        BDIM(J) = INDICESDIMS(NBAD + J - 1)
                        BID(J) = BACKWARDNEIGHBORS(NBAD + J - 1)
                        LEVEL(BDIM(J)) =
     &                        INDICESLEVS(NBAD + J - 1)
                  END DO

C                 Try each dimension as a new direction
                  DO I = 1, D

                        ISADMISSIBLE = 1

C                       Check admissibility
                        DO J = 1, NBID
                              IF (I .EQ. BDIM(J)) GOTO 20
                              BACKOFNEW = FORWARDNEIGHBORS(
     &                              (BID(J)-1)*D + I)
                              IF (BACKOFNEW .GT. 0) THEN
                                    IF (ACTIVEINDICES(BACKOFNEW)
     &                                    .EQ. 1) THEN
                                          ISADMISSIBLE = 0
                                          GOTO 30
                                    END IF
                                    BN(BDIM(J)) = BACKOFNEW
                              ELSE
                                    ISADMISSIBLE = 0
                                    GOTO 30
                              END IF
 20                           CONTINUE
                        END DO

 30                     IF (ISADMISSIBLE .EQ. 1) THEN
                              BN(I) = CURRENTINDEX

C                             Advance NI for the new subgrid
                              IF (NI .GT. 0) THEN
                                    ADDR = INDICESADDR(NI) +
     &                                    INDICESNDIIMS(NI)
                              END IF
                              NI = NI + 1
                              IF (NI .GT. NSUBGRIDS) GOTO 999
                              INDICESADDR(NI) = ADDR

                              IF (NBID .EQ. 0) THEN
C                                   First non-zero dimension
                                    NNEWDIMS = 1
                                    INDICESDIMS(ADDR) = I
                                    INDICESLEVS(ADDR) = 1
                                    BACKWARDNEIGHBORS(ADDR) = BN(I)
                                    FORWARDNEIGHBORS((BN(I)-1)*D+I)
     &                                    = NI
                                    NPOINTS = 2
                                    ADDR = ADDR + 1
                              ELSE
                                    INSERT = 1
                                    NNEWDIMS = 0
                                    NPOINTS = 1
                                    J = 1

 40                                 IF (J .LE. NBID) THEN
                                          NNEWDIMS = NNEWDIMS + 1
                                          DIM = BDIM(J)
                                          IF (I .GT. DIM) THEN
                                                LEV = LEVEL(DIM)
                                                J = J + 1
                                          ELSE IF (I .EQ. DIM) THEN
                                                INSERT = 0
                                                LEV = LEVEL(DIM)+1
                                                J = J + 1
                                          ELSE IF (INSERT .EQ. 1) THEN
                                                INSERT = 0
                                                LEV = 1
                                                DIM = I
                                          ELSE
                                                LEV = LEVEL(DIM)
                                                J = J + 1
                                          END IF

                                          FORWARDNEIGHBORS(
     &                                          (BN(DIM)-1)*D+DIM)
     &                                          = NI
                                          INDICESDIMS(ADDR) = DIM
                                          INDICESLEVS(ADDR) = LEV
                                          BACKWARDNEIGHBORS(ADDR) =
     &                                          BN(DIM)

C                                         Compute point count for dim
                                          IF (GRIDCODE .EQ. 0) THEN
                                                IF (LEV .LT. 3) THEN
                                                      NPOINTS=NPOINTS*2
                                                ELSE
                                                      NPOINTS=NPOINTS*
     &                                               2**(LEV-1)
                                                END IF
                                          ELSE IF (GRIDCODE .EQ. 1)
     &                                          THEN
                                                IF (LEV .EQ. 0) THEN
                                                      NPOINTS=NPOINTS*3
                                                ELSE
                                                      NPOINTS=NPOINTS*
     &                                               2**LEV
                                                END IF
                                          ELSE
                                                NPOINTS = NPOINTS *
     &                                                2**LEV
                                          END IF

                                          ADDR = ADDR + 1
                                          GOTO 40
                                    END IF

C                                   Handle remaining new dim if not inserted
                                    IF (INSERT .EQ. 1) THEN
                                          NNEWDIMS = NNEWDIMS + 1
                                          INDICESDIMS(ADDR) = I
                                          INDICESLEVS(ADDR) = 1
                                          BACKWARDNEIGHBORS(ADDR)=BN(I)
                                          FORWARDNEIGHBORS(
     &                                         (BN(I)-1)*D+I) = NI
                                          IF (GRIDCODE .EQ. 0) THEN
                                                NPOINTS = NPOINTS * 2
                                          ELSE IF (GRIDCODE .EQ. 1)
     &                                          THEN
                                                NPOINTS = NPOINTS * 2
                                          ELSE
                                                NPOINTS = NPOINTS * 2
                                          END IF
                                          ADDR = ADDR + 1
                                    END IF
                              END IF

                              INDICESNDIIMS(NI) = NNEWDIMS
                              SUBGRIDPOINTS(NI) = NPOINTS
                              IF (NI .GT. 1) THEN
                                    SUBGRIDADDR(NI) = SUBGRIDADDR(NI-1)
     &                                    + SUBGRIDPOINTS(NI-1)
                              ELSE
                                    SUBGRIDADDR(NI) = 1
                              END IF
                              ACTIVEINDICES(NI) = 1
                        END IF

                  END DO

                  CURRENTINDEX = CURRENTINDEX + 1

            END DO
      END DO

 999  CONTINUE

      RETURN
      END
