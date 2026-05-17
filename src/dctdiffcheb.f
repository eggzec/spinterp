C     *******************************************************************
C
C       DCT_DIFF_CHEB - DCT-based derivative of Chebyshev polynomial.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE DCT_DIFF_CHEB(N, Z, NPTS, NCOLS, Y, NY, IPDER)
C     *******************************************************************
C
C       DCT_DIFF_CHEB computes the derivative of the Chebyshev polynomial
C       at the new nodes using DCT-based coefficient extraction.
C
C       N    = full Chebyshev node count  (2^lev + 1)
C       NPTS = new surplus values only:  lev=1 -> 2,  lev>1 -> (N-1)/2
C       NCOLS = 1 (univariate) or NY (one column per query point)
C
C       Algorithm (MATLAB dctdiffcheb):
C         n==3:  ZEXT = [Z(2), 0, Z(1), 0]
C         n>3:   ZEXT(2i)   = Z(NPTS+1-i)  i=1..NPTS  (reversed)
C                ZEXT(2NPTS+2i) = Z(i)      i=1..NPTS  (direct)
C         Then DFT -> Chebyshev coeffs -> differentiate -> Clenshaw eval.
C
C       Parameters:
C
C         Input,  INTEGER N           - full Chebyshev node count
C         Input,  DOUBLE PRECISION Z(NPTS,NCOLS) - new surplus values
C         Input,  INTEGER NPTS        - number of new surplus values
C         Input,  INTEGER NCOLS       - 1 (univariate) or NY
C         Input,  DOUBLE PRECISION Y(NY) - eval points in [0,1]
C         Input,  INTEGER NY          - number of eval points
C         Output, DOUBLE PRECISION IPDER(NY) - accumulated derivative
C
C     *******************************************************************

      IMPLICIT NONE

      DOUBLE PRECISION PI
      PARAMETER (PI = 3.14159265358979323846D+00)

      INTEGER N, NPTS, NCOLS, NY
      DOUBLE PRECISION Z(NPTS, NCOLS)
      DOUBLE PRECISION Y(NY)
      DOUBLE PRECISION IPDER(NY)

      INTEGER MAXM
      PARAMETER (MAXM = 2050)

      INTEGER M, I, J, K, COL
      DOUBLE PRECISION ZEXT(MAXM)
      DOUBLE PRECISION CHEBCOEF(MAXM)
      DOUBLE PRECISION DCHEB(MAXM)
      DOUBLE PRECISION ANGLE, T, B0, B1, B2, SUMRE
      LOGICAL UNIVARIATE

      DO K = 1, NY
            IPDER(K) = 0.0D+00
      END DO

      IF (N .LE. 1 .OR. NPTS .LE. 0) RETURN

      UNIVARIATE = (NCOLS .EQ. 1)
      M = 2 * (N - 1)

      DO COL = 1, NCOLS

            DO I = 1, M
                  ZEXT(I) = 0.0D+00
            END DO

            IF (N .EQ. 3) THEN
C                 Special: new nodes at x=0 (Z(1)) and x=1 (Z(2))
C                 ZEXT = [Z(2), 0, Z(1), 0]
                  IF (NPTS .GE. 2) THEN
                        ZEXT(1) = Z(2, COL)
                        ZEXT(3) = Z(1, COL)
                  ELSE
                        ZEXT(1) = Z(1, COL)
                  END IF
            ELSE
C                 ZEXT(2i) = Z(NPTS+1-i) for i=1..NPTS
C                 ZEXT(2*NPTS + 2i) = Z(i) for i=1..NPTS
                  DO I = 1, NPTS
                        ZEXT(2*I) = Z(NPTS + 1 - I, COL)
                  END DO
                  DO I = 1, NPTS
                        IF (2*NPTS + 2*I .LE. M)
     &                      ZEXT(2*NPTS + 2*I) = Z(I, COL)
                  END DO
            END IF

C           DFT of ZEXT -> Chebyshev coefficients
            DO K = 0, N - 1
                  SUMRE = 0.0D+00
                  DO J = 0, M - 1
                        ANGLE = -2.0D+00*PI*DBLE(J)*DBLE(K)/DBLE(M)
                        SUMRE = SUMRE + ZEXT(J+1) * DCOS(ANGLE)
                  END DO
                  CHEBCOEF(K+1) = SUMRE / DBLE(N-1)
            END DO
            CHEBCOEF(1) = CHEBCOEF(1) / 2.0D+00
            CHEBCOEF(N) = CHEBCOEF(N) / 2.0D+00

C           Differentiate Chebyshev coefficients (MATLAB backward recurrence)
            DO K = 1, N
                  DCHEB(K) = 0.0D+00
            END DO
            IF (N .GE. 3) THEN
                  DCHEB(N-1) = 2.0D+00*DBLE(N-1)*CHEBCOEF(N)
                  DO K = N - 2, 2, -1
                        DCHEB(K) = DCHEB(K+2)
     &                        + 2.0D+00*DBLE(K)*CHEBCOEF(K+1)
                  END DO
                  DCHEB(1) = DCHEB(3)/2.0D+00 + CHEBCOEF(2)
            ELSE IF (N .EQ. 2) THEN
                  DCHEB(1) = CHEBCOEF(2)
            END IF

C           Scale by 2 for [0,1] domain
            DO K = 1, N - 1
                  DCHEB(K) = DCHEB(K) * 2.0D+00
            END DO

C           Clenshaw evaluation
            IF (UNIVARIATE) THEN
                  DO I = 1, NY
                        T = 2.0D+00 * Y(I) - 1.0D+00
                        B0 = 0.0D+00
                        B1 = 0.0D+00
                        B2 = 0.0D+00
                        DO K = N - 1, 2, -1
                              B2 = B1
                              B1 = B0
                              B0 = DCHEB(K) + 2.0D+00*T*B1 - B2
                        END DO
                        IPDER(I) = IPDER(I) + DCHEB(1) + T*B0 - B1
                  END DO
            ELSE
                  T = 2.0D+00 * Y(COL) - 1.0D+00
                  B0 = 0.0D+00
                  B1 = 0.0D+00
                  B2 = 0.0D+00
                  DO K = N - 1, 2, -1
                        B2 = B1
                        B1 = B0
                        B0 = DCHEB(K) + 2.0D+00*T*B1 - B2
                  END DO
                  IPDER(COL) = IPDER(COL) + DCHEB(1) + T*B0 - B1
            END IF

      END DO

      RETURN
      END
