C     *******************************************************************
C
C       SPINTERP_CB_SP - Chebyshev interpolation, sparse index structure.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPINTERP_CB_SP(D, Z, NZ, Y, NINTERP,
     &    INDICESNDIIMS, NSUBGRIDS,
     &    INDICESDIMS, INDICESLEVS, INDICESADDR, NADDR,
     &    SUBGRIDPOINTS, IP)
C     *******************************************************************
C
C       SPINTERP_CB_SP evaluates the Chebyshev sparse grid interpolant
C       using the sparse index structure (packed active dims/levels).
C       For each subgrid with active dims, sort by level descending and
C       call BARY_PD_STEP_CB.
C
C       Parameters:
C
C         Input,  INTEGER D                          - dimension
C         Input,  DOUBLE PRECISION Z(NZ)             - hierarchical surpluses
C         Input,  INTEGER NZ                         - length of Z
C         Input,  DOUBLE PRECISION Y(NINTERP,D)      - query points
C         Input,  INTEGER NINTERP                    - number of query points
C         Input,  INTEGER INDICESNDIIMS(NSUBGRIDS)   - ndims per subgrid
C         Input,  INTEGER NSUBGRIDS                  - number of subgrids
C         Input,  INTEGER INDICESDIMS(NADDR)          - packed dim indices
C         Input,  INTEGER INDICESLEVS(NADDR)          - packed levels
C         Input,  INTEGER INDICESADDR(NSUBGRIDS)      - start addr per subgrid
C         Input,  INTEGER NADDR                       - length of packed arrays
C         Input,  INTEGER SUBGRIDPOINTS(NSUBGRIDS)    - points per subgrid
C         Output, DOUBLE PRECISION IP(NINTERP)        - interpolated values
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER D, NZ, NINTERP, NSUBGRIDS, NADDR
      DOUBLE PRECISION Z(NZ)
      DOUBLE PRECISION Y(NINTERP, D)
      INTEGER INDICESNDIIMS(NSUBGRIDS)
      INTEGER INDICESDIMS(NADDR)
      INTEGER INDICESLEVS(NADDR)
      INTEGER INDICESADDR(NSUBGRIDS)
      INTEGER SUBGRIDPOINTS(NSUBGRIDS)
      DOUBLE PRECISION IP(NINTERP)

      INTEGER CI, K, DID, NDIMS, ADDR
      INTEGER LVAL, NPTS, XTOT
      INTEGER ACTIVEDIMS(50), ACTIVELEV(50)
      INTEGER ALLNX(50), DIMS(50), ORDER(50)
      INTEGER TMP, I, J
      DOUBLE PRECISION XBUF(16384)
      DOUBLE PRECISION IPTEMP(65536)
      INTEGER INDEX

      DO K = 1, NINTERP
            IP(K) = 0.0D+00
      END DO

      INDEX = 1

      DO CI = 1, NSUBGRIDS

            NDIMS = INDICESNDIIMS(CI)
            ADDR  = INDICESADDR(CI)
            NPTS  = SUBGRIDPOINTS(CI)

            IF (NDIMS .EQ. 0) THEN
                  DO K = 1, NINTERP
                        IP(K) = IP(K) + Z(INDEX)
                  END DO
                  INDEX = INDEX + NPTS
                  GOTO 200
            END IF

C           Read active dims and levels
            DO DID = 1, NDIMS
                  ACTIVEDIMS(DID) = INDICESDIMS(ADDR + DID - 1)
                  ACTIVELEV(DID)  = INDICESLEVS(ADDR + DID - 1)
                  ORDER(DID) = DID
            END DO

C           Sort by level descending (insertion sort)
            DO I = 2, NDIMS
                  TMP = ORDER(I)
                  J = I - 1
 10               IF (J .GE. 1 .AND.
     &                ACTIVELEV(ORDER(J)) .LT. ACTIVELEV(TMP)) THEN
                        ORDER(J+1) = ORDER(J)
                        J = J - 1
                        GOTO 10
                  END IF
                  ORDER(J+1) = TMP
            END DO

C           Build ALLNX and DIMS in sorted order
            XTOT = 0
            DO I = 1, NDIMS
                  DID = ORDER(I)
                  LVAL = ACTIVELEV(DID)
                  ALLNX(I) = 2**LVAL + 1
                  DIMS(I)  = ACTIVEDIMS(DID)
                  XTOT = XTOT + ALLNX(I)
            END DO

            CALL GET_CHEB_NODES(ALLNX, NDIMS, XBUF, XTOT)

            DO K = 1, NINTERP
                  IPTEMP(K) = 0.0D+00
            END DO
            CALL BARY_PD_STEP_CB(Z(INDEX), NPTS,
     &            ALLNX, DIMS, NDIMS,
     &            XBUF, XTOT, Y, NINTERP, D, IPTEMP)

            DO K = 1, NINTERP
                  IP(K) = IP(K) + IPTEMP(K)
            END DO

            INDEX = INDEX + NPTS

 200        CONTINUE
      END DO

      RETURN
      END
