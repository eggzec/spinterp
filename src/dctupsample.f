C     *******************************************************************
C
C       DCT_UPSAMPLE - Upsample 1D grid data to finer CC grid using DCT.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE DCT_UPSAMPLE(Z, NZ, NCOLS, NNEW, IP, NOUT)
C     *******************************************************************
C
C       DCT_UPSAMPLE upsamples 1D grid data from N nodes to NNEW nodes
C       using DCT (Chebyshev extension) and returns only the even-indexed
C       new nodes (the new hierarchical nodes at the finer level).
C
C       Algorithm:
C         1. Extend Z: ZEXT = [Z; flip(Z(2:N-1))] length 2*(N-1)
C         2. FFT of ZEXT
C         3. Scale: chebc = FFT[1:N] / (N-1) * (NNEW-1)
C         4. Half the N-th coefficient
C         5. Pad with zeros to length 2*(NNEW-1)
C         6. IFFT
C         7. Return even-indexed values: ip = ipext(2,4,...,NNEW)
C
C       NOTE: DFT/IDFT implemented as O(N^2) for correctness.
C       For production use FFTPACK or similar.
C
C       Parameters:
C
C         Input,  DOUBLE PRECISION Z(NZ,NCOLS) - input grid values (N rows)
C         Input,  INTEGER NZ                   - rows of Z (=N)
C         Input,  INTEGER NCOLS                - columns of Z
C         Input,  INTEGER NNEW                 - size of new full grid
C         Output, DOUBLE PRECISION IP(NOUT,NCOLS) - upsampled new nodes
C         Input,  INTEGER NOUT                 - output rows = (NNEW-1)/2
C
C     *******************************************************************

      IMPLICIT NONE

      DOUBLE PRECISION PI
      PARAMETER (PI = 3.14159265358979323846D+00)

      INTEGER NZ, NCOLS, NNEW, NOUT
      DOUBLE PRECISION Z(NZ, NCOLS)
      DOUBLE PRECISION IP(NOUT, NCOLS)

      INTEGER MAXM
      PARAMETER (MAXM = 1024)

      INTEGER N, M, MN, I, J, K, COL
      DOUBLE PRECISION ZEXT(MAXM)
      DOUBLE PRECISION CREAL(MAXM), CIMAG(MAXM)
      DOUBLE PRECISION CHEBC_R(MAXM), CHEBC_I(MAXM)
      DOUBLE PRECISION IPEXT_R(MAXM), IPEXT_I(MAXM)
      DOUBLE PRECISION ANGLE, SCALE_IN, SCALE_OUT, SUMRE, SUMIM

      N = NZ
      M = 2 * (N - 1)
      MN = 2 * (NNEW - 1)

      DO COL = 1, NCOLS

C           Build extended array
            DO I = 1, N
                  ZEXT(I) = Z(I, COL)
            END DO
            DO I = 2, N - 1
                  ZEXT(M - I + 2) = Z(I, COL)
            END DO

C           Forward DFT of ZEXT length M
            SCALE_IN = 1.0D+00 / DBLE(N-1)
            DO K = 0, N - 1
                  SUMRE = 0.0D+00
                  SUMIM = 0.0D+00
                  DO J = 0, M - 1
                        ANGLE = -2.0D+00*PI*DBLE(J)*DBLE(K)/DBLE(M)
                        SUMRE = SUMRE + ZEXT(J+1) * DCOS(ANGLE)
                        SUMIM = SUMIM + ZEXT(J+1) * DSIN(ANGLE)
                  END DO
                  CHEBC_R(K+1) = SUMRE * SCALE_IN * DBLE(NNEW-1)
                  CHEBC_I(K+1) = SUMIM * SCALE_IN * DBLE(NNEW-1)
            END DO

C           Halve the N-th coefficient (index N, 1-based)
            CHEBC_R(N) = CHEBC_R(N) / 2.0D+00
            CHEBC_I(N) = CHEBC_I(N) / 2.0D+00

C           Pad CHEBC to length MN with zeros
            DO K = 1, N
                  CREAL(K) = CHEBC_R(K)
                  CIMAG(K) = CHEBC_I(K)
            END DO
            DO K = N + 1, MN
                  CREAL(K) = 0.0D+00
                  CIMAG(K) = 0.0D+00
            END DO

C           Inverse DFT of length MN
            SCALE_OUT = 1.0D+00 / DBLE(MN)
            DO I = 0, MN - 1
                  SUMRE = 0.0D+00
                  DO K = 0, MN - 1
                        ANGLE = 2.0D+00*PI*DBLE(K)*DBLE(I)/DBLE(MN)
                        SUMRE = SUMRE + CREAL(K+1)*DCOS(ANGLE)
     &                               - CIMAG(K+1)*DSIN(ANGLE)
                  END DO
                  IPEXT_R(I+1) = SUMRE * SCALE_OUT
            END DO

C           Return even-indexed values: 2, 4, ..., NNEW (1-based)
C           i.e. ipext(2), ipext(4), ..., ipext(NNEW)
            DO I = 1, NOUT
                  IP(I, COL) = IPEXT_R(2*I)
            END DO

      END DO

      RETURN
      END
