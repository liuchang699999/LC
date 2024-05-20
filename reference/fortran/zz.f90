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
          WRITE (8,"(/20X,16HNODEL LOADS)")
          WRITE (8,'(16XA)') ' NO.DISPVALUE'
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
   END IF
END

