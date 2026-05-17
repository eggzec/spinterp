C     *******************************************************************
C
C       CHEB_WEIGHTS - 1-D Chebyshev-Lobatto quadrature weights.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE CHEB_WEIGHTS(MAXLEV, WEIGHTS, NW, STARTID)
C     *******************************************************************
C
C       CHEB_WEIGHTS builds the flat 1-D weight array and level start
C       indices for SPQUADW_CB, mirroring MATLAB chebweights.m.
C
C       Parameters:
C
C         Input,  INTEGER MAXLEV               - max Chebyshev level
C         Output, DOUBLE PRECISION WEIGHTS(NW) - flat weight array
C         Input,  INTEGER NW                   - size of WEIGHTS
C         Output, INTEGER STARTID(MAXLEV+1)    - level start indices
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER MAXLEV, NW
      DOUBLE PRECISION WEIGHTS(NW)
      INTEGER STARTID(*)

      INTEGER LEV, I, WID, NW1
      DOUBLE PRECISION W(512)

      WID = 0

      DO LEV = 0, MAXLEV

            IF (LEV .EQ. 0) THEN
                  W(1) = 1.0D+00
                  NW1 = 1
            ELSE IF (LEV .EQ. 1) THEN
                  W(1) = 1.0D+00 / 6.0D+00
                  NW1 = 1
            ELSE IF (LEV .EQ. 2) THEN
                  W(1) = 4.0D+00 / 15.0D+00
                  NW1 = 1
            ELSE IF (LEV .EQ. 3) THEN
                  W(1) = 0.7310932460800910D-01
                  W(2) = 0.1808589293602449D+00
                  NW1 = 2
            ELSE IF (LEV .EQ. 4) THEN
                  W(1) = 0.1868435141860280D-01
                  W(2) = 0.5445277629094545D-01
                  W(3) = 0.8158633214085165D-01
                  W(4) = 0.9625693230646280D-01
                  NW1 = 4
            ELSE IF (LEV .EQ. 5) THEN
                  W(1) = 0.4696598981477508D-02
                  W(2) = 0.1422895833861684D-01
                  W(3) = 0.2313138141887588D-01
                  W(4) = 0.3113605477264700D-01
                  W(5) = 0.3794190022069424D-01
                  W(6) = 0.4328876922091372D-01
                  W(7) = 0.4697162221938437D-01
                  W(8) = 0.4884909410402779D-01
                  NW1 = 8
            ELSE IF (LEV .EQ. 6) THEN
                  W( 1) = 0.1175745337655852D-02
                  W( 2) = 0.3596346580868057D-02
                  W( 3) = 0.5961697357106385D-02
                  W( 4) = 0.8267493828644795D-02
                  W( 5) = 0.1049313721486872D-01
                  W( 6) = 0.1261753249087738D-01
                  W( 7) = 0.1462032659873417D-01
                  W( 8) = 0.1648227328498816D-01
                  W( 9) = 0.1818546014331959D-01
                  W(10) = 0.1971349435647805D-01
                  W(11) = 0.2105166555570905D-01
                  W(12) = 0.2218708961962866D-01
                  W(13) = 0.2310883375546279D-01
                  W(14) = 0.2380802229262510D-01
                  W(15) = 0.2427792242857052D-01
                  W(16) = 0.2451400921551278D-01
                  NW1 = 16
            ELSE
                  NW1 = 2**(LEV-2)
                  DO I = 1, NW1
                        W(I) = 0.0D+00
                  END DO
            END IF

            STARTID(LEV+1) = WID + 1
            DO I = 1, NW1
                  WEIGHTS(WID + I) = W(I)
            END DO
            WID = WID + NW1

            IF (LEV .GT. 0) THEN
                  DO I = 1, NW1
                        WEIGHTS(WID + I) = W(NW1 - I + 1)
                  END DO
                  WID = WID + NW1
            END IF

      END DO

      RETURN
      END
