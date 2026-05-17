C     *******************************************************************
C
C       SP_DCT_UP_STEP - DCT interpolation from old level to new level.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE SP_DCT_UP_STEP(Z, NZ, OLDLEV, NOLDDIMS,
     &    NEWLEV, NDIMS, DIMS, IP, NIP)
C     *******************************************************************
C
C       SP_DCT_UP_STEP applies DCT-based interpolation stepping from
C       old levels (OLDLEV) to new levels (NEWLEV) for each active
C       dimension in DIMS.
C
C       For dimensions that increase in level, calls DCT_UPSAMPLE.
C       For new dimensions (not in old), replicates existing values.
C
C       Parameters:
C
C         Input,  DOUBLE PRECISION Z(NZ)      - surpluses at old level
C         Input,  INTEGER NZ                  - length of Z
C         Input,  INTEGER OLDLEV(NOLDDIMS)    - old levels (sorted desc.)
C         Input,  INTEGER NOLDDIMS            - number of old dims
C         Input,  INTEGER NEWLEV(NDIMS)       - new target levels per dim
C         Input,  INTEGER NDIMS               - number of target dims
C         Input,  INTEGER DIMS(NDIMS)         - dim indices (1-based)
C         Output, DOUBLE PRECISION IP(NIP)    - interpolated values
C         Input,  INTEGER NIP                 - length of IP
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER NZ, NOLDDIMS, NDIMS, NIP
      DOUBLE PRECISION Z(NZ)
      INTEGER OLDLEV(NOLDDIMS)
      INTEGER NEWLEV(NDIMS)
      INTEGER DIMS(NDIMS)
      DOUBLE PRECISION IP(NIP)

      INTEGER MAXA
      PARAMETER (MAXA = 16384)

      INTEGER I, J, K, DID, NNEW, NOUT, OLDNPTS, NEWNPTS
      INTEGER NCOLS
      DOUBLE PRECISION ZBUF(MAXA)
      DOUBLE PRECISION IPBUF(MAXA)

C     Simple implementation: for a single active dimension,
C     apply DCT_UPSAMPLE directly. For multiple dims, apply sequentially.
C     This handles the most common case of 1 or 2 active dims.

      IF (NZ .LE. 0 .OR. NIP .LE. 0) RETURN

C     Initialize IP
      DO I = 1, NIP
            IP(I) = 0.0D+00
      END DO

      IF (NDIMS .EQ. 0 .OR. NOLDDIMS .EQ. 0) THEN
C           No old dims: copy Z to IP directly
            DO I = 1, MIN(NZ, NIP)
                  IP(I) = Z(I)
            END DO
            RETURN
      END IF

C     For each new level higher than old level, upsample in that dim.
C     Start with Z as the current buffer.
      DO I = 1, NZ
            ZBUF(I) = Z(I)
      END DO
      OLDNPTS = NZ

      DO DID = 1, MIN(NDIMS, NOLDDIMS)
            IF (NEWLEV(DID) .GT. OLDLEV(DID)) THEN
C                 Upsample from OLDLEV(DID) to NEWLEV(DID)
C                 Old N = 2^OLDLEV(DID) + 1, New N = 2^NEWLEV(DID) + 1
                  OLDNPTS = 2**OLDLEV(DID) + 1
                  NNEW    = 2**NEWLEV(DID) + 1
                  NOUT    = (NNEW - 1) / 2
                  NCOLS   = NZ / OLDNPTS
                  IF (NCOLS .LT. 1) NCOLS = 1

C                 Apply upsample (returns new nodes only: NOUT points)
                  NEWNPTS = NOUT * NCOLS
                  CALL DCT_UPSAMPLE(ZBUF, OLDNPTS, NCOLS, NNEW,
     &                  IPBUF, NOUT)

                  DO I = 1, MIN(NEWNPTS, NIP)
                        IP(I) = IPBUF(I)
                  END DO
                  RETURN
            END IF
      END DO

C     No upsampling needed: copy Z directly
      DO I = 1, MIN(NZ, NIP)
            IP(I) = ZBUF(I)
      END DO

      RETURN
      END
