      SUBROUTINE SDISSIP_ARD (GFL, FLD, SL, IJS, IJL, KIJS, KIJL, &
     &                        USNEW, THWNEW, ROAIRN)
! ----------------------------------------------------------------------

!**** *SDISSIP_ARD* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.

!     LOTFI AOUF       METEO FRANCE 2013
!     FABRICE ARDHUIN  IFREMER  2013


!*    PURPOSE.
!     --------
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SDISSIP_ARD (GFL, FLD, IJS, IJL, KIJS, KIJL, SL,*
!                            USNEW, THWNEW,ROAIRN)*
!          *GFL*    - SPECTRUM.
!          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *IJS:IJL - 1st DIMENSION of GFL
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *SL*     - TOTAL SOURCE FUNCTION ARRAY
!          *USNEW*  - NEW FRICTION VELOCITY IN M/S.
!          *ROAIRN* - AIR DENSITY IN KG/M3
!          *THWNEW* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC.


!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       ARDHUIN et AL. JPO DOI:10.1175/20110JPO4324.1


! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWFRED  , ONLY : FR      , TH     ,ZPIFR
      USE YOWPCONS , ONLY : G        ,ZPI    ,ROWATER
      USE YOWPARAM , ONLY : NANG    ,NFRE
      USE YOWPHYS  , ONLY : SDSBR   ,ISDSDTH ,ISB     ,IPSAT    ,      &
&                  SSDSC2  , SSDSC4, SSDSC6,  MICHE, SSDSC3, SSDSBRF1, &
&                  BRKPBCOEF ,SSDSC5, NSDSNTH, NDIKCUMUL,              &
&                  INDICESSAT, SATWEIGHTS, CUMULW
      USE YOWSHAL  , ONLY : TFAK    ,INDEP, TCGOND
      USE YOWSTAT  , ONLY : ISHALLO
      USE YOMHOOK   ,ONLY : LHOOK   ,DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: GFL

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT) :: FLD, SL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: USNEW, THWNEW, ROAIRN 


      INTEGER(KIND=JWIM) :: IJ, K, M, I, J, M2, K2, KK, NANGD
      INTEGER(KIND=JWIM), DIMENSION(NANG) :: KKD

      REAL(KIND=JWRB) :: TPIINV, TPIINVH, TMP01, TMP03
      REAL(KIND=JWRB) :: EPSR, ROG, SSDSC6M1
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      REAL(KIND=JWRB), DIMENSION(NFRE) :: SIG, SSDSC2_SIG 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: FACTURB
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: XK
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: FACSAT, FACWTRB, TEMP1 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: BTH0 !saturation spectrum
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: BTH !saturation spectrum 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: TEMP2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: D
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: SCUMUL 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: RENEWALFREQ
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,0:NANG/2,NFRE) :: WCUMUL

      LOGICAL :: LLSSDSC3,  LLSSDSC5

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SDISSIP_ARD',0,ZHOOK_HANDLE)

      ! INITIALISATION

      EPSR=SQRT(SDSBR)

      TPIINV = 1.0_JWRB/ZPI
      TPIINVH= 0.5_JWRB*TPIINV

      NANGD=NANG/2

      ROG = ROWATER*G

      LLSSDSC3=(SSDSC3.NE.0.0_JWRB)

      TMP03 = 1.0_JWRB/(SDSBR*MICHE)

      LLSSDSC5=(SSDSC5.NE.0.0_JWRB)

      SSDSC6M1=1._JWRB-SSDSC6

      DO M=1,NFRE
        SIG(M) = ZPIFR(M)
        SSDSC2_SIG(M)=SSDSC2*SIG(M)
      END DO

      IF (ISHALLO.EQ.0) THEN
        DO M=1, NFRE
          DO IJ=KIJS,KIJL
            XK(IJ,M) = TFAK(INDEP(IJ),M)
            FACSAT(IJ,M) = XK(IJ,M)**3*TPIINV*TCGOND(INDEP(IJ),M)
          ENDDO
        ENDDO
      ELSE
        DO M=1, NFRE
          DO IJ=KIJS,KIJL
            XK(IJ,M) = (SIG(M)**2)/G
            FACSAT(IJ,M) = XK(IJ,M)**3*TPIINVH*G/SIG(M)
          ENDDO
        ENDDO
      ENDIF


      ! COMPUTE SATURATION SPECTRUM
      DO M=1,NFRE
        DO IJ=KIJS,KIJL
          BTH0(IJ,M) = 0.0_JWRB
        ENDDO
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            BTH(IJ,K,M)=0.0_JWRB
          ENDDO
        ENDDO
      ENDDO

      DO M=1,NFRE
        DO K=1,NANG
          ! integrates in directional sector
          DO K2=1,NSDSNTH*2+1
            KK=INDICESSAT(K,K2)
            DO IJ=KIJS,KIJL
              BTH(IJ,K,M) = BTH(IJ,K,M) + SATWEIGHTS(K,K2)*GFL(IJ,KK,M)
            ENDDO
          ENDDO
          DO IJ=KIJS,KIJL
            BTH(IJ,K,M)=BTH(IJ,K,M)*FACSAT(IJ,M)
            BTH0(IJ,M)=MAX(BTH0(IJ,M),BTH(IJ,K,M))
          ENDDO
        ENDDO
      ENDDO


      ! SATURATION TERM

      DO  M=1,NFRE
        DO IJ=KIJS,KIJL
          TEMP1(IJ,M)=SSDSC6*(MAX(0._JWRB,BTH0(IJ,M)*TMP03-SSDSC4))**IPSAT
        ENDDO
      ENDDO
      DO  M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            D(IJ,K,M)= SSDSC2_SIG(M)*(TEMP1(IJ,M)+SSDSC6M1*(MAX(0._JWRB,BTH(IJ,K,M)*TMP03-SSDSC4))**IPSAT)
          ENDDO
        ENDDO
      ENDDO


      ! CUMULATIVE TERM
      IF (LLSSDSC3) THEN

        DO M2=1,NFRE-NDIKCUMUL
          DO IJ=KIJS,KIJL
            IF(BTH0(IJ,M2).GT.SDSBR) THEN
              TEMP1(IJ,M2)=1.0_JWRB
            ELSE
              TEMP1(IJ,M2)=0.0_JWRB
            ENDIF
          ENDDO
        ENDDO
        DO M2=1,NFRE-NDIKCUMUL
          DO K2=1,NANG
            DO IJ=KIJS,KIJL
              SCUMUL(IJ,K2,M2)=TEMP1(IJ,M2)*(MAX(SQRT(BTH(IJ,K2,M2))-EPSR,0.0_JWRB))**2
            ENDDO
          ENDDO
        ENDDO

        DO M=NDIKCUMUL+1,NFRE
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              RENEWALFREQ(IJ,K,M)=0.0_JWRB
            ENDDO
          ENDDO
        ENDDO


        DO M=NDIKCUMUL+1,NFRE

          DO M2=1,M-NDIKCUMUL
            DO KK=0,NANGD
              DO IJ=KIJS,KIJL
                WCUMUL(IJ,KK,M2)=CUMULW(INDEP(IJ),KK,M2,M)
              ENDDO
            ENDDO
          ENDDO

          DO K=1,NANG
          ! Correction of saturation level for shallow-water kinematics
          ! Cumulative effect based on lambda   (breaking probability is
          ! the expected rate of sweeping by larger breaking waves)

            DO K2=1,NANG
              KKD(K2)=ABS(K2-K)
              IF(KKD(K2).GT.NANGD) KKD(K2)=KKD(K2)-NANGD
            ENDDO

            DO M2=1,M-NDIKCUMUL
              DO K2=1,NANG
                KK=KKD(K2)
                DO IJ=KIJS,KIJL
                  ! Integrates over frequencies M2 and directions K2 to 
                  ! Integration is performed from M2=1 to a frequency lower than M: IK-NDIKCUMUL
                  RENEWALFREQ(IJ,K,M)=RENEWALFREQ(IJ,K,M)+ WCUMUL(IJ,KK,M2)*SCUMUL(IJ,K2,M2)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO

        DO M=NDIKCUMUL+1,NFRE
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              D(IJ,K,M)= D(IJ,K,M) + RENEWALFREQ(IJ,K,M)
            ENDDO
          ENDDO
        ENDDO
      ENDIF


!     WAVE-TURBULENCE INTERACTION TERM
      IF (LLSSDSC5) THEN
        TMP01 = 2._JWRB*SSDSC5/ROG
        DO IJ=KIJS,KIJL
          FACTURB(IJ) = TMP01*ROAIRN(IJ)*USNEW(IJ)*USNEW(IJ)
        ENDDO
        DO M=1, NFRE
          DO IJ=KIJS,KIJL
            FACWTRB(IJ,M) = SIG(M)*XK(IJ,M)*FACTURB(IJ)
          ENDDO
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              D(IJ,K,M)= D(IJ,K,M)- FACWTRB(IJ,M)*COS(THWNEW(IJ)-TH(K))
            ENDDO
          ENDDO
        ENDDO
      ENDIF


      ! ADD ALL CONTRIBUTIONS TO SOURCE TERM
      DO  M=1, NFRE
        DO K=1, NANG
          DO IJ=KIJS,KIJL
            SL(IJ,K,M) = SL(IJ,K,M)+D(IJ,K,M)*GFL(IJ,K,M)
            FLD(IJ,K,M) = FLD(IJ,K,M)+D(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('SDISSIP_ARD',1,ZHOOK_HANDLE)

      END SUBROUTINE SDISSIP_ARD
