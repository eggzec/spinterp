C     *******************************************************************
C
C       SPINTERP_GP_SP - GP interpolation, sparse index structure.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SPINTERP_GP_SP(D, Z, NZ, Y, NINTERP,
     &    INDICESNDIIMS, NSUBGRIDS,
     &    INDICESDIMS, INDICESLEVS, INDICESADDR, NADDR,
     &    SUBGRIDPOINTS, IP)
C     *******************************************************************
C
C       SPINTERP_GP_SP evaluates the Gauss-Patterson sparse grid
C       interpolant using the sparse index structure.
C       For each subgrid, ALLNX(k) = 2^(lev+1)-1, calls BARY_PD_STEP_GP.
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
      INTEGER ALLNX(50), DIMS(50)
      DOUBLE PRECISION XBUF(16384), WBUF(16384)
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
            END DO

C           Build ALLNX and DIMS
            XTOT = 0
            DO DID = 1, NDIMS
                  LVAL = ACTIVELEV(DID)
                  ALLNX(DID) = 2**(LVAL+1) - 1
                  DIMS(DID)  = ACTIVEDIMS(DID)
                  XTOT = XTOT + ALLNX(DID)
            END DO

            CALL GET_GP_NODES(ALLNX, NDIMS, XBUF, XTOT)
            CALL GET_GP_BARY_W(ALLNX, NDIMS, WBUF, XTOT)

            DO K = 1, NINTERP
                  IPTEMP(K) = 0.0D+00
            END DO
            CALL BARY_PD_STEP_GP(Z(INDEX), NPTS,
     &            ALLNX, DIMS, NDIMS,
     &            XBUF, XTOT, Y, NINTERP, D, WBUF, IPTEMP)

            DO K = 1, NINTERP
                  IP(K) = IP(K) + IPTEMP(K)
            END DO

            INDEX = INDEX + NPTS

 200        CONTINUE
      END DO

      RETURN
      END
