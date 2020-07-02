      SUBROUTINE SDISSIP_ARD (F, FL, SL, IJS, IJL, &
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

!       *CALL* *SDISSIP_ARD (F, FL, IJS, IJL, SL,*
!                            USNEW, THWNEW,ROAIRN)*
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY
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

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM) :: IJ, K, M, I, J, M2, K2, KK, NANGD
      INTEGER(KIND=JWIM), DIMENSION(NANG) :: KKD

      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: USNEW, THWNEW, ROAIRN 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: FL,SL

      REAL(KIND=JWRB) :: XK(IJS:IJL,NFRE)
      REAL(KIND=JWRB) :: TPIINV, TPIINVH, TMP01, TMP03
      REAL(KIND=JWRB) :: EPSR
      REAL(KIND=JWRB) :: ROG
      REAL(KIND=JWRB) :: SSDSC6M1
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: SIG(NFRE)
      REAL(KIND=JWRB) :: SSDSC2_SIG(NFRE)
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: FACTURB
      REAL(KIND=JWRB) :: FACSAT(IJS:IJL,NFRE)
      REAL(KIND=JWRB) :: FACWTRB(IJS:IJL,NFRE)
      REAL(KIND=JWRB) :: TEMP1(IJS:IJL,NFRE)
      REAL(KIND=JWRB) :: BTH0(IJS:IJL,NFRE)  !saturation spectrum 
      REAL(KIND=JWRB) :: BTH(IJS:IJL,NANG,NFRE)  !saturation spectrum 
      REAL(KIND=JWRB) :: TEMP2(IJS:IJL,NANG,NFRE)
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: D
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: SCUMUL 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: RENEWALFREQ
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,0:NANG/2,NFRE) :: WCUMUL

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
          DO IJ=IJS,IJL
            XK(IJ,M) = TFAK(INDEP(IJ),M)
            FACSAT(IJ,M) = XK(IJ,M)**3*TPIINV*TCGOND(INDEP(IJ),M)
          ENDDO
        ENDDO
      ELSE
        DO M=1, NFRE
          DO IJ=IJS,IJL
            XK(IJ,M) = (SIG(M)**2)/G
            FACSAT(IJ,M) = XK(IJ,M)**3*TPIINVH*G/SIG(M)
          ENDDO
        ENDDO
      ENDIF


      ! COMPUTE SATURATION SPECTRUM
      DO M=1,NFRE
        DO IJ=IJS,IJL
          BTH0(IJ,M) = 0.0_JWRB
        ENDDO
        DO K=1,NANG
          DO IJ=IJS,IJL
            BTH(IJ,K,M)=0.0_JWRB
          ENDDO
        ENDDO
      ENDDO

      DO M=1,NFRE
        DO K=1,NANG
          ! integrates in directional sector
          DO K2=1,NSDSNTH*2+1
            KK=INDICESSAT(K,K2)
            DO IJ=IJS,IJL
              BTH(IJ,K,M) = BTH(IJ,K,M) + SATWEIGHTS(K,K2)*F(IJ,KK,M)
            ENDDO
          ENDDO
          DO IJ=IJS,IJL
            BTH(IJ,K,M)=BTH(IJ,K,M)*FACSAT(IJ,M)
            BTH0(IJ,M)=MAX(BTH0(IJ,M),BTH(IJ,K,M))
          ENDDO
        ENDDO
      ENDDO


      ! SATURATION TERM

      DO  M=1,NFRE
        DO IJ=IJS,IJL
          TEMP1(IJ,M)=SSDSC6*(MAX(0._JWRB,BTH0(IJ,M)*TMP03-SSDSC4))**IPSAT
        ENDDO
      ENDDO
      DO  M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            D(IJ,K,M)= SSDSC2_SIG(M)*(TEMP1(IJ,M)+SSDSC6M1*(MAX(0._JWRB,BTH(IJ,K,M)*TMP03-SSDSC4))**IPSAT)
          ENDDO
        ENDDO
      ENDDO


      ! CUMULATIVE TERM
      IF (LLSSDSC3) THEN

        DO M2=1,NFRE-NDIKCUMUL
          DO IJ=IJS,IJL
            IF(BTH0(IJ,M2).GT.SDSBR) THEN
              TEMP1(IJ,M2)=1.0_JWRB
            ELSE
              TEMP1(IJ,M2)=0.0_JWRB
            ENDIF
          ENDDO
        ENDDO
        DO M2=1,NFRE-NDIKCUMUL
          DO K2=1,NANG
            DO IJ=IJS,IJL
              SCUMUL(IJ,K2,M2)=TEMP1(IJ,M2)*(MAX(SQRT(BTH(IJ,K2,M2))-EPSR,0.0_JWRB))**2
            ENDDO
          ENDDO
        ENDDO

        DO M=NDIKCUMUL+1,NFRE
          DO K=1,NANG
            DO IJ=IJS,IJL
              RENEWALFREQ(IJ,K,M)=0.0_JWRB
            ENDDO
          ENDDO
        ENDDO


        DO M=NDIKCUMUL+1,NFRE

          DO M2=1,M-NDIKCUMUL
            DO KK=0,NANGD
              DO IJ=IJS,IJL
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
                DO IJ=IJS,IJL
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
            DO IJ=IJS,IJL
              D(IJ,K,M)= D(IJ,K,M) + RENEWALFREQ(IJ,K,M)
            ENDDO
          ENDDO
        ENDDO
      ENDIF


!     WAVE-TURBULENCE INTERACTION TERM
      IF (LLSSDSC5) THEN
        TMP01 = 2._JWRB*SSDSC5/ROG
        DO IJ=IJS,IJL
          FACTURB(IJ) = TMP01*ROAIRN(IJ)*USNEW(IJ)*USNEW(IJ)
        ENDDO
        DO M=1, NFRE
          DO IJ=IJS,IJL
            FACWTRB(IJ,M) = SIG(M)*XK(IJ,M)*FACTURB(IJ)
          ENDDO
          DO K=1,NANG
            DO IJ=IJS,IJL
              D(IJ,K,M)= D(IJ,K,M)- FACWTRB(IJ,M)*COS(THWNEW(IJ)-TH(K))
            ENDDO
          ENDDO
        ENDDO
      ENDIF


      ! ADD ALL CONTRIBUTIONS TO SOURCE TERM
      DO  M=1, NFRE
        DO K=1, NANG
          DO IJ=IJS,IJL
            SL(IJ,K,M) = SL(IJ,K,M)+D(IJ,K,M)*F(IJ,K,M)
            FL(IJ,K,M) = FL(IJ,K,M)+D(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('SDISSIP_ARD',1,ZHOOK_HANDLE)

      END SUBROUTINE SDISSIP_ARD
