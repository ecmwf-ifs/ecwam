      SUBROUTINE NLWEIGT (ML, KL)

! ----------------------------------------------------------------------

!**** *NLWEIGT* - COMPUTATION OF INDEX ARRAYS AND WEIGHTS FOR THE
!                 COMPUTATION OF THE NONLINEAR TRANSFER RATE.

!     SUSANNE HASSELMANN JUNE 86.

!     H. GUNTHER   ECMWF/GKSS  DECEMBER 90 - CYCLE_4 MODIFICATIONS.
!                                            4 FREQUENCIES ADDED.

!*    PURPOSE.
!     --------

!       COMPUTATION OF PARAMETERS USED IN DISCRETE INTERACTION
!       PARAMETERIZATION OF NONLINEAR TRANSFER.

!**   INTERFACE.
!     ----------

!       *CALL* *NLWEIGT (ML, KL)*
!          *ML*     INTEGER   NUMBER OF FREQUENCIES.
!          *KL*     INTEGER   NUMBER OF DIRECTIONS.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       *JAFU*      - FUNCTION FOR COMPUTATION OF ANGULAR INDICES
!                     OF K(F,THET).

!     REFERENCE.
!     ----------
!       S. HASSELMANN AND K. HASSELMANN, JPO, 1985 B.


! ----------------------------------------------------------------------

      USE YOWFRED  , ONLY : FR       ,DELTH    ,FRATIO
      USE YOWINDN  , ONLY : IKP      ,IKP1     ,IKM      ,IKM1     ,
     &            K1W      ,K2W      ,K11W     ,K21W     ,AF11     ,
     &            FKLAP    ,FKLAP1   ,FKLAM    ,FKLAM1   ,ACL1     ,
     &            ACL2     ,CL11     ,CL21     ,DAL1     ,DAL2     ,
     &            FRH      ,KFRH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : DEG
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

!*    *PARAMETER*  FOR DISCRETE APPROXIMATION OF NONLINEAR TRANSFER

      PARAMETER (ALAMD=0.25, CON=3000., DELPHI1=-11.48, DELPHI2=33.56)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *ALAMD*     REAL      LAMBDA
!      *CON*       REAL      WEIGHT FOR DISCRETE APPROXIMATION OF
!                            NONLINEAR TRANSFER
!      *DELPHI1*   REAL
!      *DELPHI2*   REAL

! ----------------------------------------------------------------------
!     LOCAL ARRAYS
      INTEGER, ALLOCATABLE :: JA1(:,:)
      INTEGER, ALLOCATABLE :: JA2(:,:)
      REAL, ALLOCATABLE :: FRLON(:)

! ----------------------------------------------------------------------

!     0. ALLOCATE ARRAYS
!        ---------------

      ALLOCATE(JA1(NANG,2))
      ALLOCATE(JA2(NANG,2))
      ALLOCATE(FRLON(2*NFRE+2))

      ALLOCATE(IKP(NFRE+4))
      ALLOCATE(IKP1(NFRE+4))
      ALLOCATE(IKM(NFRE+4))
      ALLOCATE(IKM1(NFRE+4))
      ALLOCATE(K1W(NANG,2))
      ALLOCATE(K2W(NANG,2))
      ALLOCATE(K11W(NANG,2))
      ALLOCATE(K21W(NANG,2))
      ALLOCATE(AF11(NFRE+4))
      ALLOCATE(FKLAP(NFRE+4))
      ALLOCATE(FKLAP1(NFRE+4))
      ALLOCATE(FKLAM(NFRE+4))
      ALLOCATE(FKLAM1(NFRE+4))
      ALLOCATE(FRH(KFRH))

!*    1. COMPUTATION FOR ANGULAR GRID.
!        -----------------------------

      DELTHA = DELTH*DEG
      CL1 = DELPHI1/DELTHA
      CL2 = DELPHI2/DELTHA

!*    1.1 COMPUTATION OF INDICES OF ANGULAR CELL.
!         ---------------------------------------

      KLP1 = KL+1
      IC = 1
      DO KH=1,2
        KLH = KL
        IF (KH.EQ.2) KLH=KLP1
        DO K=1,KLH
          KS = K
          IF (KH.GT.1) KS=KLP1-K+1
          IF (KS.GT.KL) GO TO 1002
          CH = IC*CL1
          JA1(KS,KH) = JAFU(CH,K,KLP1)
          CH = IC*CL2
          JA2(KS,KH) = JAFU(CH,K,KLP1)
 1002     CONTINUE
        ENDDO
        IC = -1
      ENDDO

!*    1.2 COMPUTATION OF ANGULAR WEIGHTS.
!         -------------------------------

      ICL1 = CL1
      CL1  = CL1-ICL1
      ICL2 = CL2
      CL2  = CL2-ICL2
      ACL1 = ABS(CL1)
      ACL2 = ABS(CL2)
      CL11 = 1.-ACL1
      CL21 = 1.-ACL2
      AL11 = (1.+ALAMD)**4
      AL12 = (1.-ALAMD)**4
      DAL1 = 1./AL11
      DAL2 = 1./AL12

!*    1.3 COMPUTATION OF ANGULAR INDICES.
!         -------------------------------

      ISG = 1
      DO KH=1,2
        CL1H = ISG*CL1
        CL2H = ISG*CL2
        DO K=1,KL
          KS = K
          IF (KH.EQ.2) KS = KL-K+2
          IF(K.EQ.1) KS = 1
          K1 = JA1(K,KH)
          K1W(KS,KH) = K1
          IF (CL1H.LT.0.) THEN
            K11 = K1-1
            IF (K11.LT.1) K11 = KL
          ELSE
            K11 = K1+1
            IF (K11.GT.KL) K11 = 1
          ENDIF
          K11W(KS,KH) = K11
          K2 = JA2(K,KH)
          K2W(KS,KH) = K2
          IF (CL2H.LT.0) THEN
            K21 = K2-1
            IF(K21.LT.1) K21 = KL
          ELSE
            K21 = K2+1
            IF (K21.GT.KL) K21 = 1
          ENDIF
          K21W(KS,KH) = K21
        ENDDO
        ISG = -1
      ENDDO

!*    2. COMPUTATION FOR FREQUENCY GRID.
!        -------------------------------

      DO M=1,ML
        FRLON(M) = FR(M)
      ENDDO
      DO M=ML+1,2*ML+2
        FRLON(M) = FRATIO*FRLON(M-1)
      ENDDO
      F1P1 = LOG10(FRATIO)
      DO M=1,ML+4
        FRG = FRLON(M)
        AF11(M) = CON * FRG**11
        FLP = FRG*(1.+ALAMD)
        FLM = FRG*(1.-ALAMD)
        IKN = INT(LOG10(1.+ALAMD)/F1P1+.000001)
        IKN = M+IKN
        IKP(M) = IKN
        FKP = FRLON(IKP(M))
        IKP1(M) = IKP(M)+1
        FKLAP(M) = (FLP-FKP)/(FRLON(IKP1(M))-FKP)
        FKLAP1(M) = 1.-FKLAP(M)
        IF (FRLON(1).GE.FLM) THEN
          IKM(M) = 1
          IKM1(M) = 1
          FKLAM(M) = 0.
          FKLAM1(M) = 0.
        ELSE
          IKN = INT(LOG10(1.-ALAMD)/F1P1+.0000001)
          IKN = M+IKN-1
          IF (IKN.LT.1) IKN = 1
          IKM(M) = IKN
          FKM = FRLON(IKM(M))
          IKM1(M) = IKM(M)+1
          FKLAM(M) = (FLM-FKM)/(FRLON(IKM1(M))-FKM)
          FKLAM1(M) = 1.-FKLAM(M)
        ENDIF
      ENDDO

!*    3. COMPUTE TAIL FREQUENCY RATIOS.
!        ------------------------------

      IE = MIN(30,ML+3)
      DO I=1,IE
        M = ML+I-1
        FRH(I) = (FRLON(ML)/FRLON(M))**5
      ENDDO

!*    4. PRINTER PROTOCOL.
!        -----------------

      WRITE(IU06,'(1H1,'' NON LINEAR INTERACTION PARAMETERS:'')')
      WRITE(IU06,'(1H0,'' COMMON INDNL: CONSTANTS'')')
      WRITE(IU06,'(1X,''    ACL1       ACL2   '',
     &             ''    CL11       CL21   '',
     &             ''    DAL1       DAL2'')')
      WRITE(IU06,'(1X,6F11.8)') ACL1, ACL2, CL11, CL21, DAL1, DAL2

      WRITE(IU06,'(1H0,'' COMMON INDNL: FREQUENCY ARRAYS'')')
      WRITE(IU06,'(1X,'' M   IKP IKP1  IKM IKM1'',
     &          ''   FKLAP       FKLAP1 '',
     &          ''   FKLAM       FKLAM1     AF11'')')
      DO M=1,ML+4
        WRITE(IU06,'(1X,I2,4I5,4F11.8,E11.3)')
     &   M, IKP(M), IKP1(M), IKM(M), IKM1(M),
     &   FKLAP(M), FKLAP1(M), FKLAM(M), FKLAM1(M), AF11(M)
      ENDDO

      WRITE(IU06,'(1H0,'' COMMON INDNL: ANGULAR ARRAYS'')')
      WRITE(IU06,'(1X,''  |--------KH = 1----------|'',
     &              ''|--------KH = 2----------|'')')
      WRITE(IU06,'(1X,'' K   K1W   K2W  K11W  K21W'',
     &              ''   K1W   K2W  K11W  K21W'')')
      DO K=1,KL
        WRITE(IU06,'(1X,I2,8I6)') K,(K1W(K,KH), K2W(K,KH), K11W(K,KH),
     &   K21W(K,KH),KH=1,2)
      ENDDO
      WRITE(IU06,'(1H0,'' COMMON INDNL: TAIL ARRAY FRH'')')
      WRITE(IU06,'(1X,8F10.7)') (FRH(M),M=1,30)

!     5. DEALLOCATE LOCAL ARRAYS
!        -----------------------

      DEALLOCATE(JA1)
      DEALLOCATE(JA2)
      DEALLOCATE(FRLON)

      RETURN
      END SUBROUTINE NLWEIGT
