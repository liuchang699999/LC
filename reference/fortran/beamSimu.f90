PROGRAM SFSAP
!    C ANALYSIS PROGRAM FOR SPACE FRAME
   REAL :: K(100,100), KE(12,12), AKE(12,12), X(50), Y(50), Z(50), AL(50), EAI(6,20), PJ(20), PF(2,20), &
      R(12,12), P(50), D_D(20), E, G, D_E, DMAX, FE(12), D(50), ADE(12), DE(12), RT(12,12), F(6), &
      AFE(12), AF(12)

   OPEN (6,FILE='SFSAP.IN')
   OPEN (8,FILE='SFSAP.OUT')
   CALL READ(NJ,NJJ,N,NE,NM,NPJ,NPF,JN,X,Y,Z,JE,JEAI,JPJ,PJ,JPF,PF)
   DO IE=1,NE
      D_D(IE)=0.01 !初化杆直径?
   END DO

10 DO I=1,N
      P(I)=0.
      DO J=1,N
         K(I,J)=0.
      END DO
   END DO
   E=2.0E+11
   PI=3.14159
   G=7.672311E+10

   DO IE=1,NE
      CALL MKE(KE,IE,JE,JEAI,EAI,X,Y,Z,AL,D_D)
      CALL MR(R,IE,JE,X,Y,Z)
      CALL MAKE(KE,R,AKE)
      CALL CALM(M,IE,JN,JE)
      CALL MK(K,AKE,M) !形成总刚
   END DO
   DO IP=1,NPF
      CALL MR(R,JPF(1,IP),JE,X,Y,Z)
      CALL TRAN(R,RT)
      CALL PE(FE,IP,JPF,PF,AL)
      CALL MULV12 (RT,FE,AFE)
      CALL CALM(M,JPF(1,IP),JN,JE)
      CALL MF(P,AFE,M)
   END DO
   DO I=1,NPJ
      P(JPJ(I))=P(JPJ(I))+PJ(I)  !形成节点载荷?
   END DO
   CALL SLOV(K,P,D,N)
   WRITE(8,'(/2(22(1H*),A))') 'RESULTS OF CALCULATION'
   WRITE(8,'(/27X,A)')'DISPLACEMENT'
   WRITE(8,40)
40 FORMAT(/'NO.E',5X,'DX',10X,'DY',10X,'DZ',10X,'RX',10X,'RY',10X,'RZ')
   DO KK=1,NJ
      DO II=1,6
         F(II)=0.
         I1=JN(II,KK)
         IF(I1.GT.0)F(II)=D(I1)
      END DO
      WRITE(8,70)KK,F(1),F(2),F(3),F(4),F(5),F(6)!读出节点位移?
70    FORMAT(I2,2X,6G12.4)
   END DO
   WRITE(8,80)
80 FORMAT(/'单元?',2X,'内力DMAX',2X,'临界应力D_E',2X,'直径D_D',2X)
   DO IE=1,NE
      CALL MADE(IE,JN,JE,D,ADE)
      CALL MULV12(R,ADE,DE)!求局部坐标系下的单元杆�??力�?
      II=JE(1,IE)
      JJ=JE(2,IE)
      MT=JEAI(IE)
      barL=SQRT((X(JJ)-X(II))**2+(Y(JJ)-Y(II))**2+(Z(JJ)-Z(II))**2)
      DE7=DE(7)
      DE1=DE(1)
      DE2=DE(2)
      DE3=DE(3)
      DE5=DE(5)
      DE6=DE(6)
      DE4=DE(4)
      DE8=DE(8)
      DE12=DE(12)
      DE10=DE(10)
      DE9=DE(9)
      DE11=DE(11)
      D1MAX=(DE7-DE1)*E/barL
      D2MAX=3.*D_D(IE)*E*(DE8-DE2)/barL**2
      tmd11=D_D(IE)*E*(DE12+2.*DE6)/barL
      D2MAX=D2MAX+tmd11
      D3MAX=3.*D_D(IE)*E*(DE9-DE3)/barL**2
      tmd11=D_D(IE)*E*(DE11+2.*DE5)/barL
      D3MAX=D3MAX+tmd11
      DWMAX=SQRT(D2MAX**2+D3MAX**2)
      DTMAX=0.5*D_D(IE)*G*(DE10-DE4)/barL
      DMAX=SQRT((D1MAX+DWMAX)**2+3*DTMAX**2)  !求杆�?最大内应力�?
      D_E=41*pI**2*E*D_D(IE)**2/(1600*barL**2) !求临界应力�?
      if(DMAX.GT.D_E) THEN
         D_D(IE)=D_D(IE)+0.0001 !线性搜索最优直径。
         WRITE(8,90)ie, dmax,d_e,d_d(ie)
90       FORMAT(I2,2X,3G12.4)

         !WRITE(*,*) ie, dmax,d_e,d_d(ie)
      END If
   END DO
   !pause
   GOTO 10 !循环语句形成线性搜索。

END PROGRAM SFSAP


SUBROUTINE READ(NJ,NJJ,N,NE,NM,NPJ,NPF,JN,X,Y,Z,JE, JEAI,JPJ,PJ,JPF,PF)
   REAL :: X(50),Y(50),Z(50),PJ(20),PF(2,20)
   INTEGER :: JE(3,50),JN(6,50),JPJ(20),JPF(3,20),TITLE(20),JEAI(20)
   READ(6,'(20A4)') (TITLE(I),I=1,20)
   WRITE(8,'(7X,20A4)') TITLE
   READ (6,*) NJ,NJJ,N,NE,NPJ,NPF
   WRITE(8,'(/4(5X,A4,1H:,I2))') 'NJ=',NJ,'NJJ=',NJJ,'N=',N,'NE=',NE,'NPJ=',NPJ,'NPF=',NPF
   WRITE (8,'(/A75)') 'NO(1) (2) (3) (4) (5) (6) X Y Z'
   DO I=1,NJ
      READ (6,*) JN(1,I),JN(2,I),JN(3,I),JN(4,I),JN(5,I),JN(6,I),X(I),Y(I),Z(I)
   END DO
   DO I=1,NJ
      WRITE (8,'(2X,1H(,I2,1H),6I6,4X,3F10.3)') I,JN(1,I),JN(2,I),JN(3,I),JN(4,I),JN(5,I),JN(6,I),X(I),Y(I),Z(I)
   END DO
   WRITE (8,'(/,A38)')'V-NODE NO.XYZ'
   READ (6,30) (X(I),Y(I),Z(I),I=NJ+1,NJ+NJJ)
30 FORMAT (3(3G16.4))
   DO I=NJ+1,NJ+NJJ
      WRITE (8,'(2X,1H(,I2,1H),4X,3F10.3)') I,X(I),Y(I),Z(I)
   END DO
   WRITE (8,'(/)')
   WRITE (8,*)'ELEMENT NO.NODE-1 NODE-2 NODE-VMATERIALS'
   READ (6,50) (JE(1,I),JE(2,I),JE(3,I),JEAI(I),I=1,NE)
50 FORMAT (3(4I5))
   DO I=1,NE
      WRITE (8,'(4X,I2,4(7X,I3))') I,JE(1,I),JE(2,I),JE(3,I),JEAI(I)
   END DO
   IF(NPJ.NE.0)   THEN
      IF(NPJ.NE.0)   THEN
         !  WRITE (8,"(/20X,16HNODEL LOADS)")
         !  WRITE (8,'(16XA)') ' NO.DISPVALUE'
         READ (6,90)(JPJ(I),PJ(I),I=1,NPJ)
90       FORMAT (3(I5,G16.4))
         DO I=1,NPJ
            WRITE(8,'(14X,I7,F16.3)') JPJ(I),PJ(I)
         END DO
      ELSE
      END IF
      IF(NPF.NE.0) THEN
         !  WRITE (8,'(/20X,16HNON-NODEL LOADS)')
         !  WRITE (8,'(A,9X,A,10X,A)') 'NO.ENO.LOAD.MODEL SURFACE','A','C'
         READ (6,120) (JPF(1,I),JPF(2,I),JPF(3,I),PF(1,I),PF(2,I),I=1,NPF)
120      FORMAT (2(3I5,2G16.4))
         DO I=1,NPF
            WRITE(8,'(3(I3,7X),2F10.3)')(JPF(J,I),J=1,3),PF(1,I),PF(2,I)
         END DO
      ELSE
      END IF
      RETURN
   END

   SUBROUTINE MKE(KE,IE,JE,JEAI,EAI,X,Y,Z,AL,D_D)!子程序求单刚
      REAL :: KE(12,12),X(50),Y(50),Z(50),EAI(6,20),AL(20),barL,D_D(50),PI
      INTEGER :: JE(3,50),JEAI(20)
      PI=3.14159
      II=JE(1,IE)
      JJ=JE(2,IE)
      MT=JEAI(IE)
      barL=SQRT((X(JJ)-X(II))**2+(Y(JJ)-Y(II))**2+(Z(JJ)-Z(II))**2)
      AL(IE)=barL
      EAI(1,MT)=2.0E+11
      EAI(2,MT)=7.672311E+10
      EAI(3,MT)=9.*PI*D_D(IE)**2/25.
      EAI(4,MT)=369.*PI*D_D(IE)**4/20000.
      EAI(5,MT)=369.*PI*D_D(IE)**4/40000.
      EAI(6,MT)=369.*PI*D_D(IE)**4/40000.
      A1=EAI(1,MT)*EAI(3,MT)/barL
      A2=EAI(1,MT)*EAI(6,MT)/barL**3
      A3=EAI(1,MT)*EAI(6,MT)/barL**2
      A4=EAI(1,MT)*EAI(5,MT)/barL**3
      A5=EAI(1,MT)*EAI(5,MT)/barL**2
      A6=EAI(2,MT)*EAI(4,MT)/barL
      A7=EAI(1,MT)*EAI(5,MT)/barL
      A8=EAI(1,MT)*EAI(6,MT)/barL
      KE(1,1)=A1
      KE(1,7)=-A1
      KE(2,2)=12*A2
      KE(2,6)=6*A3
      KE(2,8)=-12*A2
      KE(2,12)=6*A3
      KE(3,3)=12*A4
      KE(3,5)=-6*A5
      KE(3,9)=-12*A4
      KE(3,11)=-6*A5
      KE(4,4)=A6
      KE(4,10)=-A6
      KE(5,5)=4*A7
      KE(5,9)=6*A5
      KE(5,11)=2*A7
      KE(6,6)=4*A8
      KE(6,8)=-6*A3
      KE(6,12)=2*A8
      KE(7,7)=A1
      KE(8,8)=12*A2
      KE(8,12)=-6*A3
      KE(9,9)=12*A4
      KE(9,11)=6*A5
      KE(10,10)=A6
      KE(11,11)=4*A7
      KE(12,12)=4*A8
      DO I=1 ,12
         DO K=I ,12
            KE(K,I)=KE(I,K)
         END DO
      END DO
      RETURN
   END


   SUBROUTINE MR(R,IE,JE,X,Y,Z)!子程序求坐标矩阵
      REAL R(12,12),X(50),Y(50),Z(50),barL,L2,L3,LXX,LXY,LXZ,LYX,LYY,LYZ,LZX,LZY,LZZ
      INTEGER JE(3,50)
      I=JE(1,IE)
      J=JE(2,IE)
      K=JE(3,IE)
      barL=SQRT((X(J)-X(I))**2+(Y(J)-Y(I))**2+(Z(J)-Z(I))**2)
      LXX=(X(J)-X(I))/barL
      LXY=(Y(J)-Y(I))/barL
      LXZ=(Z(J)-Z(I))/barL
      YZ=((Y(K)-Y(I))*(Z(K)-Z(J))-(Z(K)-Z(I))*(Y(K)-Y(J)))
      ZX=-((X(K)-X(I))*(Z(K)-Z(J))-(Z(K)-Z(I))*(X(K)-X(J)))
      XY=((X(K)-X(I))*(Y(K)-Y(J))-(Y(K)-Y(I))*(X(K)-X(J)))
      L2=SQRT(YZ**2+ZX**2+XY**2)
      LZX=YZ/L2
      LZY=ZX/L2
      LZZ=XY/L2
      S1=(1-LXX**2)*(X(K)-X(I))-LXX*LXY*(Y(K)-Y(I))-LXX*LXZ*(Z(K)-Z(I))
      S2=(1-LXY**2)*(Y(K)-Y(I))-LXY*LXX*(X(K)-X(I))-LXY*LXZ*(Z(K)-Z(I))
      S3=(1-LXZ**2)*(Z(K)-Z(I))-LXZ*LXX*(X(K)-X(I))-LXZ*LXY*(Y(K)-Y(I))
      L3=SQRT(S1**2+S2**2+S3**2)
      LYX=S1/L3
      LYY=S2/L3
      LYZ=S3/L3
      DO I=1,12
         DO J=1,12
            R(I,J)=0
         END DO
      END DO
      DO II=1,10,3
         R(II,II)=LXX
         R(II,II+1)=LXY
         R(II,II+2)=LXZ
         R(II+1,II)=LYX
         R(II+1,II+1)=LYY
         R(II+1,II+2)=LYZ
         R(II+2,II)=LZX
         R(II+2,II+1)=LZY
         R(II+2,II+2)=LZZ
      END DO
      RETURN
   END
   SUBROUTINE MAKE(KE,R,AKE) !子程序求整体坐标系下单刚矩阵
      REAL KE(12,12),R(12,12),RT(12,12),TMP(12,12),AKE(12,12)
      CALL TRAN(R,RT)
      CALL MULV(RT,KE,TMP)
      CALL MULV(TMP,R,AKE)
      RETURN
   END


   SUBROUTINE TRAN(R,RT) !子程序求�?�?
      REAL R(12,12),RT(12,12)
      DO I=1,12
         DO J=1,12
            RT(I,J)=R(J,I)
         END DO
      END DO
      RETURN
   END

   SUBROUTINE MULV(A,B,C)
      REAL A(12,12),B(12,12),C(12,12)
      DO I=1,12
         DO J=1,12
            C(I,J)=0.
            DO K=1,12
               C(I,J)=C(I,J)+A(I,K)*B(K,J)
            END DO
         END DO
      END DO
      RETURN
   END

   SUBROUTINE CALM(M,IE,JN,JE)
      INTEGER M(12),JN(6,50),JE(3,50),IE
      DO I=1,6
         M(I)=JN(I,JE(1,IE))
         M(I+6)=JN(I,JE(2,IE))
      END DO
      RETURN
   END

   SUBROUTINE MK(K,AKE,M)
      REAL K(100,100),AKE(12,12)
      INTEGER M(12)
      DO I=1,12
         DO J=1,12
            IF(M(I).NE.0.AND.M(J).NE.0)
            &K(M(I),M(J))=K(M(I),M(J))+AKE(I,J)
         END DO
      END DO
      RETURN
   END


   SUBROUTINE PE(FE,IP,JPF,PF,AL)!子程序�?�算单元等效结点���荷
      REAL :: FE(12),PF(2,20),AL(20),barL
      INTEGER :: JPF(3,20)
      A=PF(1,IP)
      C=PF(2,IP)
      barL=AL(JPF(1,IP))
      IND=JPF(2,IP)
      DO I=1, 12
         FE(I)=0
      END DO
      GOTO (10,20,30,40,50,60),IND
10    FE(2)=(7*A/20+3*C/20)*barL
      FE(6)=(A/20+C/30)*barL**2
      FE(8)=(3*A/20+7*C/20)*barL
      FE(12)=-(A/30+C/20)*barL**2
      GOTO 70
20    FE(8)=A*C**3*(2*barL-C)/2/barL**3
      FE(2)=A*C-FE(8)
      FE(6)=A*C**2*(6*barL*barL-8*C*barL+3*C*C)/12/barL/barL
      FE(12)=-A*C**3*(4*barL-3*C)/12/barL/barL
      GOTO 70
30    FE(2)=A*(barL-C)**2*(barL+2*C)/barL**3
      FE(6)=A*C*(barL-C)**2/barL**2
      FE(8)=A-FE(2)
      FE(12)=-A*C**2*(barL-C)/barL**2
      GOTO 70
40    FE(2)=-6*A*C*(barL-C)/barL**3
      FE(6)=A*(barL-C)*(barL-3*C)/barL**2
      FE(8)=-FE(2)
      FE(12)=A*C*(3*C-2*barL)/barL**2
      GOTO 70
50    FE(1)=A*(1-C/barL)
      FE(7)=A*C/barL
      GOTO 70
60    FE(1)=C*barL/2.
      FE(7)=FE(1)
70    CONTINUE
      IF (JPF(3,IP).EQ.2) THEN
         P=FE(5)
         FE(5)=-FE(6)
         FE(6)=-P
         P=FE(2)
         FE(2)=FE(3)
         FE(3)=P
         P=FE(11)
         FE(11)=-FE(12)
         FE(12)=-P
         P=FE(8)
         FE(8)=FE(9)
         FE(9)=P
      ENDIF
      RETURN
   END

   SUBROUTINE MULV12(A,B,C)!子程序矩阵相�?
      REAL C(12),A(12,12),B(12)
      DO I=1,12
         C(I)=0.
         DO J=1,12
            C(I)=C(I)+A(I,J)*B(J)
         END DO
      END DO
      RETURN
   END

   SUBROUTINE MF(P,AFE,M)!子程序组集整体坐标系下载荷列�?
      REAL P(50),AFE(12)
      INTEGER M(12)
      DO I=1,12
         IF(M(I).NE.0)P(M(I))=AFE(I)+P(M(I))
      END DO
      RETURN
   END

   SUBROUTINE SLOV(AK,P,D,N)!子程序解方程求自由结点位�?
      REAL AK(100,100),P(50),D(50)
      DO I=1,N
         D(I)=P(I)
      END DO
      DO K=1,N-1
         DO I=K+1,N
            C=-AK(K,I)/AK(K,K)
            DO J=I,N
               AK(I,J)=AK(I,J)+C*AK(K,J)
            END DO
            D(I)=D(I)+C*D(K)
         END DO
      END DO
      D(N)=D(N)/AK(N,N)
      DO I=N-1,1,-1
         DO J=I+1,N
            D(I)=D(I)-AK(I,J)*D(J)
         END DO
         D(I)=D(I)/AK(I,I)
      END DO
      RETURN
   END

   SUBROUTINE MADE(IE,JN,JE,D,ADE) !子程序求整体坐标系下单元杆�??位移
      REAL ADE(12),D(50)
      INTEGER IE,JN(6,50),JE(3,20)
      DO I=1,12
         ADE(I)=0
      END DO
      DO I=1,6
         IF (JN(I,JE(1,IE)).NE.0) ADE(I)=D(JN(I,JE(1,IE)))
         IF (JN(I,JE(2,IE)).NE.0) ADE(I+6)=D(JN(I,JE(2,IE)))
      END DO
      RETURN
   END
