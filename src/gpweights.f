C     *******************************************************************
C
C       GP_WEIGHTS - 1-D Gauss-Patterson quadrature weights.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE GP_WEIGHTS(MAXLEV, WEIGHTS, NW, STARTID)
C     *******************************************************************
C
C       GP_WEIGHTS builds the flat 1-D weight array and level start
C       indices for SPQUADW_GP, mirroring MATLAB gpweights.m.
C
C       For each level: weights include both new and mirrored half.
C       All weights are divided by 2 (range [0,1] vs [-1,1]).
C
C       Parameters:
C
C         Input,  INTEGER MAXLEV          - max GP level (0..6)
C         Output, DOUBLE PRECISION WEIGHTS(NW) - flat weight array
C         Input,  INTEGER NW              - size of WEIGHTS
C         Output, INTEGER STARTID(MAXLEV+1)    - level start indices
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER MAXLEV, NW
      DOUBLE PRECISION WEIGHTS(NW)
      INTEGER STARTID(*)

      INTEGER LEV, I, WID, NW1
      DOUBLE PRECISION W(64)

      WID = 0

      DO LEV = 0, MAXLEV

            IF (LEV .EQ. 0) THEN
                  W(1) = 2.0D+00
                  NW1 = 1
            ELSE IF (LEV .EQ. 1) THEN
                  W(1) = 0.555555555555555555556D+00
                  NW1 = 1
            ELSE IF (LEV .EQ. 2) THEN
                  W(1) = 0.104656226026467265194D+00
                  W(2) = 0.401397414775962222905D+00
                  NW1 = 2
            ELSE IF (LEV .EQ. 3) THEN
                  W(1) = 0.170017196299402603390D-01
                  W(2) = 0.929271953151245376859D-01
                  W(3) = 0.171511909136391380787D+00
                  W(4) = 0.219156858401587496404D+00
                  NW1 = 4
            ELSE IF (LEV .EQ. 4) THEN
                  W(1) = 0.254478079156187441540D-02
                  W(2) = 0.164460498543878109338D-01
                  W(3) = 0.359571033071293220968D-01
                  W(4) = 0.569795094941233574122D-01
                  W(5) = 0.768796204990035310427D-01
                  W(6) = 0.936271099812644736167D-01
                  W(7) = 0.105669893580234809744D+00
                  W(8) = 0.111956873020953456880D+00
                  NW1 = 8
            ELSE IF (LEV .EQ. 5) THEN
                  W(1)  = 0.363221481845530659694D-03
                  W(2)  = 0.257904979468568827243D-02
                  W(3)  = 0.611550682211724633968D-02
                  W(4)  = 0.104982469096213218983D-01
                  W(5)  = 0.154067504665594978021D-01
                  W(6)  = 0.205942339159127111492D-01
                  W(7)  = 0.258696793272147469108D-01
                  W(8)  = 0.310735511116879648799D-01
                  W(9)  = 0.360644327807825726401D-01
                  W(10) = 0.407155101169443189339D-01
                  W(11) = 0.449145316536321974143D-01
                  W(12) = 0.485643304066731987159D-01
                  W(13) = 0.515832539520484587768D-01
                  W(14) = 0.539054993352660639269D-01
                  W(15) = 0.554814043565593639878D-01
                  W(16) = 0.562776998312543012726D-01
                  NW1 = 16
            ELSE IF (LEV .EQ. 6) THEN
                  W(1)  = 0.505360952078625176247D-04
                  W(2)  = 0.377746646326984660274D-03
                  W(3)  = 0.938369848542381500794D-03
                  W(4)  = 0.168114286542146990631D-02
                  W(5)  = 0.256876494379402037313D-02
                  W(6)  = 0.357289278351729964938D-02
                  W(7)  = 0.467105037211432174741D-02
                  W(8)  = 0.584344987583563950756D-02
                  W(9)  = 0.707248999543355546805D-02
                  W(10) = 0.834283875396815770558D-02
                  W(11) = 0.964117772970253669530D-02
                  W(12) = 0.109557333878379016480D-01
                  W(13) = 0.122758305600827700870D-01
                  W(14) = 0.135915710097655467896D-01
                  W(15) = 0.148936416648151820348D-01
                  W(16) = 0.161732187295777199419D-01
                  W(17) = 0.174219301594641737472D-01
                  W(18) = 0.186318482561387901863D-01
                  W(19) = 0.197954950480974994880D-01
                  W(20) = 0.209058514458120238522D-01
                  W(21) = 0.219563663053178249393D-01
                  W(22) = 0.229409642293877487608D-01
                  W(23) = 0.238540521060385400804D-01
                  W(24) = 0.246905247444876769091D-01
                  W(25) = 0.254457699654647658126D-01
                  W(26) = 0.261156733767060976805D-01
                  W(27) = 0.266966229274503599062D-01
                  W(28) = 0.271855132296247918192D-01
                  W(29) = 0.275797495664818730349D-01
                  W(30) = 0.278772514766137016085D-01
                  W(31) = 0.280764557938172466068D-01
                  W(32) = 0.281763190330166021307D-01
                  NW1 = 32
            ELSE
                  NW1 = 2**(LEV-2)
                  DO I = 1, NW1
                        W(I) = 0.0D+00
                  END DO
            END IF

C           Divide by 2 for [0,1] range
            DO I = 1, NW1
                  W(I) = W(I) / 2.0D+00
            END DO

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
