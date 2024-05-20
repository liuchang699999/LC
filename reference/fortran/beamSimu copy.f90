PROGRAM SFSAP
   ! ANALYSIS PROGRAM FOR SPACE FRAME
   REAL :: K(100,100), KE(12,12), AKE(12,12), X(50), Y(50), Z(50), AL(50), EAI(6,20), PJ(20), PF(2,20), R(12,12), P(50), D_D(20), E, G, D_E, DMAX, FE(12), D(50), ADE(12), DE(12), RT(12,12), F(6), AFE(12), AF(12)
   INTEGER :: JE(3,50), JN(6,50), JPJ(20), JPF(3,20), M(12), JEAI(40)
   OPEN (6, FILE='SFSAP.IN')
   OPEN (8, FILE='SFSAP.OUT')
   CALL READ(NJ, NJJ, N, NE, NM, NPJ, NPF, JN, X, Y, Z, JE, JEAI, JPJ, PJ, JPF, PF)
   DO IE = 1, NE
      D_D(IE) = 0.01 ! 初始化杆直径。
   END DO
10 CONTINUE
   DO I = 1, N
      P(I) = 0.
      DO J = 1, N
         K(I,J) = 0.
      END DO
   END DO
   E = 2.0E+11
   PI = 3.14159
   G = 7.672311E+10

   DO IE = 1, NE
      CALL MKE(KE, IE, JE, JEAI, EAI, X, Y, Z, AL, D_D)
      CALL MR(R, IE, JE, X, Y, Z)
      CALL MAKE(KE, R, AKE)
      CALL CALM(M, IE, JN, JE)
      CALL MK(K, AKE, M) ! 形成总刚
   END DO
   DO IP = 1, NPF
      CALL MR(R, JPF(1,IP), JE, X, Y, Z)
      CALL TRAN(R, RT)
      CALL PE(FE, IP, JPF, PF, AL)
      CALL MULV12(RT, FE, AFE)
      CALL CALM(M, JPF(1,IP), JN, JE)
      CALL MF(P, AFE, M)
   END DO
   DO I = 1, NPJ
      P(JPJ(I)) = P(JPJ(I)) + PJ(I) ! 形成节点载荷。
   END DO
   CALL SLOV(K, P, D, N)
   WRITE(8, '(// 2(22(1H*),A))') 'RESULTS OF CALCULATION'
   WRITE(8, '(// 27X,A)') 'DISPLACEMENT'
   WRITE(8, 40)
40 FORMAT(/ 'NO.E', 5X, 'DX', 10X, 'DY', 10X, 'DZ', 10X, 'RX', 10X, 'RY', 10X, 'RZ')
   DO KK = 1, NJ
      DO II = 1, 6
         F(II) = 0.
         I1 = JN(II,KK)
         IF (I1.GT.0) F(II) = D(I1)
      END DO
      WRITE(8, 70) KK, F(1), F(2), F(3), F(4), F(5), F(6) ! 读出节点位移。
70    FORMAT(I2, 2X, 6G12.4)
   END DO
   WRITE(8, 80)
80 FORMAT(/ '单元号', 2X, '内力DMAX', 2X, '临界应力D_E', 2X, '直径D_D', 2X)
   DO IE = 1, NE
      CALL MADE(IE, JN, JE, D, ADE)
      CALL MULV12(R, ADE, DE) ! 求局部坐标系下的单元杆端力。
      II = JE(1,IE)
      JJ = JE(2,IE)
      MT = JEAI(IE)
      barL = SQRT((X(JJ)-X(II))**2 + (Y(JJ)-Y(II))**2 + (Z(JJ)-Z(II))**2)
      DE7 = DE(7)
      DE1 = DE(1)
      DE2 = DE(2)
      DE3 = DE(3)
      DE5 = DE(5)
      DE6 = DE(6)
      DE4 = DE(4)
      DE8 = DE(8)
      DE12 = DE(12)
      DE10 = DE(10)
      DE9 = DE(9)
      DE11 = DE(11)
      D1MAX = (DE7-DE1)*E/barL
      D2MAX = 3.*D_D(IE)*E*(DE8-DE2)/barL**2
      tmd11 = D_D(IE)*E*(DE12+2.*DE6)/barL
      D2MAX = D2MAX + tmd11
      D3MAX = 3.*D_D(IE)*E*(DE9-DE3)/barL**2
      tmd11 = D_D(IE)*E*(DE11+2.*DE5)/barL
      D3MAX = D3MAX + tmd11
      DWMAX = SQRT(D2MAX**2 + D3MAX**2)
      DTMAX = 0.5*D_D(IE)*G*(DE10-DE4)/barL
      DMAX = SQRT((D1MAX+DWMAX)**2 + 3*DTMAX**2) ! 求杆端最大内应力。
      D_E = 41*PI**2*E*D_D(IE)**2/(1600*barL**2) ! 求临界应力。
      IF (DMAX.GT.D_E) THEN
         D_D(IE) = D_D(IE) + 0.0001 ! 线性搜索最优直径。
         WRITE(8, 90) IE, DMAX, D_E, D_D(IE)
90       FORMAT(I2, 2X, 3G12.4)
      END IF
   END DO
   ! pause
   GOTO 10 ! 循环语句形成线性搜索。
END PROGRAM SFSAP


SUBROUTINE READ(NJ, NJJ, N, NE, NM, NPJ, NPF, JN, X, Y, Z, JE, JEAI, JPJ, PJ, JPF, PF)
   REAL :: X(50), Y(50), Z(50), PJ(20), PF(2,20)
   INTEGER :: JE(3,50), JN(6,50), JPJ(20), JPF(3,20), TITLE(20), JEAI(20)
   READ(6, '(20A4)') (TITLE(I), I=1, 20)
   WRITE(8, '(7X,20A4)') TITLE
   READ(6, *) NJ, NJJ, N, NE, NM, NPJ, NPF
   WRITE(8, '(// 4(5X,A4,1H:,I2))') 'NJ=', NJ, 'NJJ=', NJJ, &
      'N=', N, 'NE=', NE, 'NPJ=', NPJ, 'NPF=', NPF
   WRITE(8, '(A75)') 'NO(1) (2) (3) (4) (5) (6) X Y Z'
   DO I = 1, NJ
      READ(6, *) JN(1,I), JN(2,I), JN(3,I), JN(4,I), JN(5,I), JN(6,I), &
         X(I), Y(I), Z(I)
   END DO
   DO I = 1, NJ
      WRITE(8, '(2X,1H(,I2,1H),6I6,4X,3F10.3)') I, JN(1,I), JN(2,I), &
         JN(3,I), JN(4,I), JN(5,I), JN(6,I), X(I), Y(I), Z(I)
   END DO
   WRITE(8, '(//,A38)') 'V-NODE NO.XYZ'
   READ(6, 30) (X(I), Y(I), Z(I), I = NJ+1, NJ+NJJ)
30 FORMAT (3(3G16.4))
   DO I = NJ+1, NJ+NJJ
      WRITE(8, '(2X,1H(,I2,1H),4X,3F10.3)') I, X(I), Y(I), Z(I)
   END DO
   WRITE(8, '(//)')
   WRITE(8, *) 'ELEMENT NO.NODE-1 NODE-2 NODE-VMATERIALS'
   READ(6, 50) (JE(1,I), JE(2,I), JE(3,I), JEAI(I), I = 1, NE)
50 FORMAT (3(4I5))
   DO I = 1, NE
      WRITE(8, '(4X,I2,4(7X,I3))') I, JE(1,I), JE(2,I), JE(3,I), JEAI(I)
   END DO
   IF (NPJ .NE. 0) THEN
      WRITE(8, '(//20X,16HNODEL LOADS)')
      WRITE(8, '(16XA)') ' NO.DISPVALUE'
      READ(6, 90) (JPJ(I), PJ(I), I = 1, NPJ)
90    FORMAT (3(I5,G16.4))
      DO I = 1, NPJ
         WRITE(8, '(14X,I7,F16.3)') JPJ(I), PJ(I)
      END DO
   END IF
   IF (NPF .NE. 0) THEN
      WRITE(8, '(//20X,16HNON-NODEL LOADS)')
      WRITE(8, '(A,9X,A,10X,A)') 'NO.ENO.LOAD.MODEL SURFACE', 'A', 'C'
      READ(6, 120) (JPF(1,I), JPF(2,I), JPF(3,I), PF(1,I), &
         PF(2,I), I = 1, NPF)
120   FORMAT (2(3I5,2G16.4))
      DO I = 1, NPF
         WRITE(8, '(3(I3,7X),2F10.3)') (JPF(J,I), J = 1, 3), PF(1,I), PF(2,I)
      END DO
   END IF
   RETURN
END

SUBROUTINE MKE(KE, IE, JE, JEAI, EAI, X, Y, Z, AL, D_D)
   REAL :: KE(12,12), X(50), Y(50), Z(50), EAI(6,20), AL(20), D_D(50), PI
   INTEGER :: JE(3,50), JEAI(20)
   PI = 3.14159
   II = JE(1,IE)
   JJ = JE(2,IE)
   MT = JEAI(IE)
   barL = SQRT((X(JJ)-X(II))**2 + (Y(JJ)-Y(II))**2 + (Z(JJ)-Z(II))**2)
   AL(IE) = barL
   EAI(1,MT) = 2.0E+11
   EAI(2,MT) = 7.672311E+10
   EAI(3,MT) = 9.*PI*D_D(IE)**2/25.
   EAI(4,MT) = 369.*PI*D_D(IE)**4/500.
   EAI(5,MT) = 1.5*EAI(2,MT)
   EAI(6,MT) = 2.625*PI*EAI(3,MT)
   EAI(3,MT) = EAI(3,MT)/barL
   EAI(4,MT) = EAI(4,MT)/barL**3
   EAI(5,MT) = EAI(5,MT)/barL
   EAI(6,MT) = EAI(6,MT)/barL**2
   DO I = 1, 12
      DO J = 1, 12
         KE(I,J) = 0.
      END DO
   END DO
   KE(1,1) = EAI(1,MT)*AL(IE)/barL + EAI(3,MT)*barL/3.
   KE(1,7) = KE(1,1) - EAI(3,MT)*barL/6.
   KE(1,8) = EAI(5,MT)*AL(IE) + EAI(6,MT)*barL/2.
   KE(1,9) = -KE(1,7)
   KE(1,10) = KE(1,8)
   KE(1,11) = 0.
   KE(1,12) = -KE(1,8)/2.
   KE(2,2) = KE(1,1)
   KE(2,4) = KE(1,8)
   KE(2,5) = -KE(1,1)
   KE(2,6) = KE(1,8)/2.
   KE(2,8) = EAI(1,MT)*AL(IE)/barL + EAI(3,MT)*barL/3.
   KE(2,9) = -KE(2,8)/2.
   KE(2,10) = -EAI(5,MT)*AL(IE) - EAI(6,MT)*barL/2.
   KE(2,11) = KE(2,8)
   KE(2,12) = -KE(2,10)
   KE(3,3) = EAI(2,MT)*AL(IE)/barL + EAI(4,MT)*barL/3.
   KE(3,7) = -KE(3,3)
   KE(3,8) = -KE(2,8)
   KE(3,9) = KE(2,8)
   KE(3,10) = 0.
   KE(3,11) = KE(3,3) - EAI(4,MT)*barL/6.
   KE(3,12) = -KE(2,8)/2.
   KE(4,4) = KE(3,3)
   KE(4,5) = -KE(3,3)
   KE(4,6) = -KE(2,8)/2.
   KE(4,8) = EAI(2,MT)*AL(IE)/barL + EAI(4,MT)*barL/3.
   KE(4,9) = KE(2,8)/2.
   KE(4,10) = KE(2,8)
   KE(4,11) = KE(3,3) - EAI(4,MT)*barL/6.
   KE(4,12) = -KE(2,8)
   KE(5,5) = EAI(2,MT)*AL(IE)/barL + EAI(4,MT)*barL/3.
   KE(5,7) = KE(3,3) - EAI(4,MT)*barL/6.
   KE(5,8) = KE(2,8)/2.
   KE(5,9) = -KE(2,8)/2.
   KE(5,10) = KE(3,3)
   KE(5,11) = -KE(2,8)
   KE(5,12) = KE(2,8)
   KE(6,6) = 2.*EAI(2,MT)*AL(IE)/barL + EAI(4,MT)*barL/3.
   KE(6,8) = KE(2,8)
   KE(6,9) = -KE(2,8)
   KE(6,10) = -2.*KE(2,8)
   KE(6,11) = -EAI(3,MT)*barL/3.
   KE(6,12) = 2.*EAI(6,MT)*barL/3.
   KE(7,7) = EAI(3,MT)*barL/3.
   KE(7,8) = 0.
   KE(7,9) = -KE(7,7)
   KE(7,10) = 0.
   KE(7,11) = -EAI(3,MT)*barL/6.
   KE(7,12) = 0.
   KE(8,8) = EAI(5,MT)*AL(IE) + EAI(6,MT)*barL/2.
   KE(8,9) = 0.
   KE(8,10) = 0.
   KE(8,11) = KE(8,8)
   KE(8,12) = -KE(8,8)/2.
   KE(9,9) = KE(8,8)
   KE(9,10) = -EAI(5,MT)*AL(IE) - EAI(6,MT)*barL/2.
   KE(9,11) = -KE(8,8)
   KE(9,12) = KE(9,10)
   KE(10,10) = EAI(5,MT)*AL(IE) + EAI(6,MT)*barL/2.
   KE(10,11) = KE(8,8)
   KE(10,12) = -KE(10,10)/2.
   KE(11,11) = EAI(3,MT)*barL/3.
   KE(11,12) = 0.
   KE(12,12) = EAI(6,MT)*barL/3.
   RETURN
END

SUBROUTINE MR(R, IE, JE, X, Y, Z)
   REAL :: R(12,12), X(50), Y(50), Z(50)
   INTEGER :: JE(3,50)
   II = JE(1,IE)
   JJ = JE(2,IE)
   XM = (X(II)+X(JJ))/2.
   YM = (Y(II)+Y(JJ))/2.
   ZM = (Z(II)+Z(JJ))/2.
   L = SQRT((X(JJ)-X(II))**2 + (Y(JJ)-Y(II))**2 + (Z(JJ)-Z(II))**2)
   A = (X(JJ)-X(II))/L
   B = (Y(JJ)-Y(II))/L
   C = (Z(JJ)-Z(II))/L
   AN = SQRT(B**2 + C**2)
   BN = -B/AN
   CN = -C/AN
   DO I = 1, 3
      DO J = 1, 3
         R(I,J) = 0.
      END DO
   END DO
   R(1,1) = A
   R(1,2) = B
   R(1,3) = C
   R(2,1) = 0.
   R(2,2) = AN
   R(2,3) = BN
   R(3,1) = 0.
   R(3,2) = CN
   R(3,3) = AN
   R(4,4) = 1.
   R(5,5) = A
   R(5,6) = B
   R(5,7) = C
   R(6,5) = 0.
   R(6,6) = AN
   R(6,7) = BN
   R(7,5) = 0.
   R(7,6) = CN
   R(7,7) = AN
   R(8,8) = 1.
   R(9,9) = 1.
   R(10,10) = A
   R(10,11) = B
   R(10,12) = C
   R(11,10) = 0.
   R(11,11) = AN
   R(11,12) = BN
   R(12,10) = 0.
   R(12,11) = CN
   R(12,12) = AN
   RETURN
END

SUBROUTINE MAKE(KE, R, AKE)
   REAL :: KE(12,12), R(12,12), AKE(12,12)
   DO I = 1, 12
      DO J = 1, 12
         AKE(I,J) = 0.
         DO K = 1, 12
            DO L = 1, 12
               AKE(I,J) = AKE(I,J) + R(I,K)*KE(K,L)*R(J,L)
            END DO
         END DO
      END DO
   END DO
   RETURN
END

SUBROUTINE CALM(M, IE, JN, JE)
   INTEGER :: M(12), JE(3,50), JN(6,50)
   DO I = 1, 12
      M(I) = 0
   END DO
   DO I = 1, 6
      JJ = JE(I,IE)
      IF (JJ.GT.0) THEN
         DO J = 1, 6
            II = JN(J,JJ)
            IF (II.GT.0) M(II) = 1
         END DO
      END IF
   END DO
   RETURN
END

SUBROUTINE MK(K, AKE, M)
   REAL :: K(100,100), AKE(12,12)
   INTEGER :: M(12)
   DO I = 1, 12
      IF (M(I).EQ.0) THEN
         DO J = 1, 12
            DO K = 1, 12
               DO L = 1, 12
                  K(12*(I-1)+J,12*(I-1)+K) = K(12*(I-1)+J,12*(I-1)+K) + AKE(J,K)
               END DO
            END DO
         END DO
      END IF
   END DO
   RETURN
END

SUBROUTINE PE(FE, IP, JPF, PF, AL)
   REAL :: FE(12), PF(2,20), AL(20)
   INTEGER :: JPF(3,20)
   DO I = 1, 12
      FE(I) = 0.
   END DO
   DO I = 1, 20
      IF (JPF(1,I).EQ.IP) THEN
         DO J = 1, 12
            FE(J) = FE(J) + PF(1,I)*AL(I)
         END DO
      ELSE IF (JPF(2,I).EQ.IP) THEN
         DO J = 1, 12
            FE(J) = FE(J) + PF(2,I)*AL(I)
         END DO
      END IF
   END DO
   RETURN
END

SUBROUTINE MF(P, AFE, M)
   REAL :: P(50), AFE(12)
   INTEGER :: M(12)
   DO I = 1, 12
      IF (M(I).EQ.0) THEN
         P(I) = P(I) + AFE(I)
      END IF
   END DO
   RETURN
END

SUBROUTINE MADE(IE, JN, JE, D, ADE)
   REAL :: D(50), ADE(12)
   INTEGER :: JE(3,50), JN(6,50)
   DO I = 1, 12
      ADE(I) = 0.
   END DO
   II = JE(1,IE)
   JJ = JE(2,IE)
   DO I = 1, 6
      JI = JN(I,II)
      JJ = JN(I,JJ)
      IF (JI.GT.0 .AND. JJ.GT.0) THEN
         DO J = 1, 12
            ADE(J) = ADE(J) + D(JI)*(JJ-1) + D(JJ)*JJ
         END DO
      END IF
   END DO
   RETURN
END

SUBROUTINE MULV12(R, AFE, DE)
   REAL :: R(12,12), AFE(12), DE(12)
   DO I = 1, 12
      DE(I) = 0.
   END DO
   DO I = 1, 12
      DO J = 1, 12
         DE(I) = DE(I) + R(I,J)*AFE(J)
      END DO
   END DO
   RETURN
END

SUBROUTINE SLOV(K, P, D, N)
   REAL :: K(100,100), P(50), D(50)
   INTEGER :: N
   CALL GAUSS(K, P, D, N)
   RETURN
END

SUBROUTINE GAUSS(K, P, D, N)
   REAL :: K(100,100), P(50), D(50)
   INTEGER :: N, I, J, K1
   DO I = 1, N-1
      DO K1 = I+1, N
         P(K1) = P(K1) - K(K1,I)/K(I,I)*P(I)
         DO J = I+1, N
            K(K1,J) = K(K1,J) - K(K1,I)/K(I,I)*K(I,J)
         END DO
      END DO
   END DO
   D(N) = P(N)/K(N,N)
   DO I = N-1, 1, -1
      D(I) = P(I)
      DO J = I+1, N
         D(I) = D(I) - K(I,J)*D(J)
      END DO
      D(I) = D(I)/K(I,I)
   END DO
   RETURN
END

