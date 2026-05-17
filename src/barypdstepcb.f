C     *******************************************************************
C
C       BARY_PD_STEP_CB - Barycentric d-dim Lagrange interpolation step,
C                         Chebyshev grid.
C
C       License:
C            Sparse Grid Interpolation Toolbox
C            Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart
C            Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
C            Copyright (c) 2026 eggzec. All Rights Reserved.
C            See LICENSE for details.
C
C     *******************************************************************

      SUBROUTINE BARY_PD_STEP_CB(Z, NZ, ALLNX, DIMS, D,
     &                            X, NX, Y, NINTERP, DYFULL, IP)
C     *******************************************************************
C
C       BARY_PD_STEP_CB evaluates the multi-dimensional polynomial
C       interpolant using barycentric Lagrange interpolation on the
C       Chebyshev-Lobatto grid.  Only the NZ new (hierarchical) nodes
C       contribute; the node structure has NX(k) total nodes per dim,
C       with new nodes at even positions (0,2,4,...).
C
C       Parameters:
C
C         Input,  DOUBLE PRECISION Z(NZ)       - hierarchical surpluses
C         Input,  INTEGER NZ                   - number of new nodes
C         Input,  INTEGER ALLNX(D)             - full node count per dim
C         Input,  INTEGER DIMS(D)              - active dimension indices
C         Input,  INTEGER D                    - active dimensions
C         Input,  DOUBLE PRECISION X(NX)       - concatenated node arrays
C         Input,  INTEGER NX                   - total nodes in X
C         Input,  DOUBLE PRECISION Y(NINTERP,DYFULL) - query points
C         Input,  INTEGER NINTERP              - number of query points
C         Input,  INTEGER DYFULL               - full dimension of Y
C         Output, DOUBLE PRECISION IP(NINTERP) - interpolated values
C
C     *******************************************************************

      IMPLICIT NONE

      INTEGER NZ, D, NX, NINTERP, DYFULL
      DOUBLE PRECISION Z(NZ)
      INTEGER ALLNX(D), DIMS(D)
      DOUBLE PRECISION X(NX)
      DOUBLE PRECISION Y(NINTERP, DYFULL)
      DOUBLE PRECISION IP(NINTERP)

      INTEGER K, L, K2, KXI, AID, XID, ACTNX, ACTNA
      INTEGER NNA(50), IDADDFIX(50), AIDSTART(50), AIDMAX(50)
      INTEGER AIDVEC(50), IDVEC(50), IDADD(50)
      INTEGER NID, NITER
      INTEGER ZIDBASE, ZID, ID(50)
      INTEGER AID1, AID2, AIDMAX1, AIDMAX2, ADDID1, ADDID2
      INTEGER ID1, ID2, IDL, AIDL, AIDLMAX
      DOUBLE PRECISION A(16384)
      DOUBLE PRECISION C1(50), C2(50), S, FRAC, T
      DOUBLE PRECISION C1D, C2D, C11, C12
      INTEGER ISZERONODE

      DO K = 1, NINTERP
            IP(K) = 0.0D+00
      END DO

C     NNA(k) = number of new (non-zero) nodes per dim k
C     For level 1 (ALLNX=3): NNA = 2; else NNA = (ALLNX-1)/2
      DO L = 1, D
            IF (ALLNX(L) .EQ. 3) THEN
                  NNA(L) = 2
            ELSE
                  NNA(L) = (ALLNX(L) - 1) / 2
            END IF
      END DO

C     Column strides in Z (new-node layout)
      IDADDFIX(1) = 1
      AIDSTART(1) = 1
      DO L = 2, D
            IDADDFIX(L) = IDADDFIX(L-1) * NNA(L-1)
            AIDSTART(L) = AIDSTART(L-1) + NNA(L-1)
      END DO
      DO L = 1, D
            AIDMAX(L) = AIDSTART(L) + NNA(L) - 1
      END DO

      DO K = 1, NINTERP
            XID = 0
            AID = 0
            ISZERONODE = 0

            DO L = 1, D
                  T = Y(K, DIMS(L))
                  S = 0.0D+00
                  ID(L) = 0
                  ACTNX = ALLNX(L)
                  ACTNA = NNA(L)

                  IF (ACTNX .EQ. 3) THEN
C                       Level 1: 3 nodes (x[1]=0, x[2]=0.5, x[3]=1)
C                       x[2] (index 2, zero node) sets ip=0 if hit
                        IF (T .EQ. X(XID+1)) THEN
                              ID(L) = 1
                              C2(L) = 1.0D+00
                              AID = AID + ACTNA
                              XID = XID + ACTNX
                              GOTO 30
                        ELSE IF (T .EQ. X(XID+2)) THEN
                              ISZERONODE = 1
                              GOTO 50
                        ELSE IF (T .EQ. X(XID+3)) THEN
                              ID(L) = 2
                              C2(L) = 1.0D+00
                              AID = AID + ACTNA
                              XID = XID + ACTNX
                              GOTO 30
                        END IF
                        FRAC = 0.5D+00 / (T - X(XID+1))
                        A(AID+1) = FRAC
                        S = S + FRAC
                        S = S - 1.0D+00 / (T - X(XID+2))
                        FRAC = 0.5D+00 / (T - X(XID+3))
                        A(AID+2) = FRAC
                        S = S + FRAC

                  ELSE
C                       General level >= 2
C                       Check if t hits a support node (even pos)
C                       or a zero node (odd pos) using acos approximation
                        IF (T .GE. 0.0D+00 .AND. T .LE. 1.0D+00) THEN
                              K2 = XID + NINT(DBLE(ACTNX) - DACOS(
     &                              (T-0.5D+00)*2.0D+00)
     &                              *DBLE(ACTNX-1)/
     &                              3.14159265358979323846D+00)
                              IF (K2 .LT. XID+1) K2 = XID+1
                              IF (K2 .GT. XID+ACTNX) K2 = XID+ACTNX
                              IF (T .EQ. X(K2)) THEN
                                    IF (MOD(K2-XID, 2) .EQ. 1) THEN
                                          ISZERONODE = 1
                                          GOTO 50
                                    ELSE
                                          ID(L) = (K2-XID)/2
                                          C2(L) = 1.0D+00
                                          AID = AID + ACTNA
                                          XID = XID + ACTNX
                                          GOTO 30
                                    END IF
                              END IF
                        END IF

                        S = S + 0.5D+00 / (T - X(XID+1))
                        KXI = XID + 2
                        DO K2 = AID+1, AID+ACTNA-1
                              FRAC = -1.0D+00 / (T - X(KXI))
                              A(K2) = FRAC
                              S = S + FRAC + 1.0D+00/(T - X(KXI+1))
                              KXI = KXI + 2
                        END DO
                        S = S - 0.5D+00 / (T - X(XID+ACTNX))
                  END IF

                  C2(L) = S
                  AID = AID + ACTNA
                  XID = XID + ACTNX
 30               CONTINUE
            END DO

            IF (ISZERONODE .NE. 0) GOTO 40

C           Determine free (non-hit) dimensions
            NITER = 1
            NID = 0
            DO L = 1, D
                  IF (ID(L) .EQ. 0) THEN
                        NITER = NITER * NNA(L)
                        NID = NID + 1
                        IDVEC(NID) = L
                  ELSE
                        ID(L) = ID(L) - 1
                  END IF
            END DO

            IF (NID .EQ. 0) THEN
                  ZIDBASE = 1
                  DO L = 1, D
                        ZIDBASE = ZIDBASE + ID(L)*IDADDFIX(L)
                  END DO
                  IP(K) = Z(ZIDBASE)
                  GOTO 40
            END IF

            ZIDBASE = 1
            DO L = 1, D
                  C1(L) = 0.0D+00
                  ZIDBASE = ZIDBASE + ID(L)*IDADDFIX(L)
                  AIDVEC(L) = AIDSTART(L)
            END DO

            DO L = 1, NID
                  IDL = IDVEC(L)
                  IF (L .GT. 1) THEN
                        IDADD(IDL) = IDADDFIX(IDL)
     &                        - IDADDFIX(IDVEC(L-1)+1)
                  ELSE
                        IDADD(IDL) = IDADDFIX(IDL)
                  END IF
            END DO

            C11 = 0.0D+00
            C12 = 0.0D+00
            ID1 = IDVEC(1)
            ADDID1 = IDADD(ID1)
            AID1 = AIDVEC(ID1)
            AIDMAX1 = AIDMAX(ID1)
            ZID = ZIDBASE

            IF (NID .EQ. 1) THEN
                  DO K2 = AID1, AIDMAX1
                        C11 = C11 + Z(ZID) * A(K2)
                        ZID = ZID + ADDID1
                  END DO
                  C1D = C11

            ELSE
                  NITER = NITER / NNA(ID1)
                  ID2 = IDVEC(2)
                  AID2 = AIDVEC(ID2)
                  AIDMAX2 = AIDMAX(ID2)
                  ADDID2 = IDADD(ID2)

                  DO K2 = 1, NITER
                        DO AID = AID1, AIDMAX1
                              C11 = C11 + Z(ZID) * A(AID)
                              ZID = ZID + ADDID1
                        END DO
                        ZID = ZID + ADDID2
                        C12 = C12 + A(AID2) * C11
                        C11 = 0.0D+00
                        AID2 = AID2 + 1

                        IF (NID .GT. 2 .AND. AID2 .GT. AIDMAX2) THEN
                              ZID = ZID + IDADD(IDVEC(3))
                              C1(3) = C1(3) + A(AIDVEC(IDVEC(3)))*C12
                              C12 = 0.0D+00
                              AID2 = AIDSTART(ID2)
                              DO L = 3, NID
                                    IDL = IDVEC(L)
                                    AIDL = AIDVEC(IDL)
                                    AIDLMAX = AIDMAX(IDL)
                                    IF (AIDL .EQ. AIDLMAX) THEN
                                          IF (L .LT. NID) THEN
                                                ZID = ZID
     &                                          + IDADD(IDVEC(L+1))
                                                C1(L+1) = C1(L+1)
     &                                          + A(AIDVEC(IDVEC(L+1)))
     &                                          * C1(L)
                                                C1(L) = 0.0D+00
                                                AIDVEC(IDL) =
     &                                               AIDSTART(IDL)
                                          END IF
                                    ELSE
                                          AIDVEC(IDL) = AIDL + 1
                                          GOTO 60
                                    END IF
                              END DO
 60                           CONTINUE
                        END IF
                  END DO

                  IF (NID .EQ. 2) THEN
                        C1D = C12
                  ELSE
                        C1D = C1(NID)
                  END IF
            END IF

            C2D = 1.0D+00
            DO L = 1, D
                  C2D = C2D * C2(L)
            END DO
            IP(K) = C1D / C2D

 40         CONTINUE
      END DO

 50   RETURN
      END
