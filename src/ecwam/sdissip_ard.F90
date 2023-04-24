! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SDISSIP_ARD (KIJS, KIJL, FL1, FLD, SL,          &
     &                        INDEP, WAVNUM, XK2CG,              &
     &                        UFRIC, COSWDIF, RAORW)
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

!       *CALL* *SDISSIP_ARD (KIJS, KIJL, FL1, FLD,SL,*
!                            INDEP, WAVNUM, XK2CG,
!                            UFRIC, COSWDIF, RAORW)*
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *FL1*    - SPECTRUM.
!          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*     - TOTAL SOURCE FUNCTION ARRAY
!          *INDEP*  - DEPTH INDEX
!          *WAVNUM* - WAVE NUMBER
!          *XK2CG*  - (WAVE NUMBER)**2 * GROUP SPEED
!          *UFRIC*  - FRICTION VELOCITY IN M/S.
!          *RAORW*  - RATIO AIR DENSITY TO WATER DENSITY
!          *COSWDIF*-  COS(TH(K)-WDWAVE(IJ))


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
      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOWPARAM , ONLY : NANG    ,NFRE, NANG_PARAM
      USE YOWPHYS  , ONLY : SDSBR   ,ISDSDTH ,ISB     ,IPSAT    ,      &
&                  SSDSC2  , SSDSC4, SSDSC6,  MICHE, SSDSC3, SSDSBRF1, &
&                  BRKPBCOEF ,SSDSC5, NSDSNTH, NDIKCUMUL,              &
&                  INDICESSAT, SATWEIGHTS, CUMULW

      USE YOMHOOK  , ONLY : LHOOK   ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FLD, SL
      INTEGER(KIND=JWIM), DIMENSION(KIJL), INTENT(IN) :: INDEP
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: WAVNUM, XK2CG 
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: UFRIC, RAORW 
      REAL(KIND=JWRB), DIMENSION(KIJL, NANG), INTENT(IN) :: COSWDIF 


      INTEGER(KIND=JWIM) :: IJ, K, M, I, J, M2, K2, KK, NANGD
      INTEGER(KIND=JWIM), DIMENSION(NANG_PARAM) :: KKD

      REAL(KIND=JWRB) :: TPIINV, TPIINVH, TMP01, TMP03
      REAL(KIND=JWRB) :: EPSR, SSDSC6M1, ZCOEF, ZCOEFM1
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


      REAL(KIND=JWRB) :: SSDSC2_SIG 
      REAL(KIND=JWRB), DIMENSION(KIJL) :: FACTURB, BTH, BTH0
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG_PARAM) :: SCUMUL, D 

      REAL(KIND=JWRB), DIMENSION(KIJL) :: RENEWALFREQ

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SDISSIP_ARD',0,ZHOOK_HANDLE)

      ! INITIALISATION

      EPSR = SQRT(SDSBR)

      TPIINV = 1.0_JWRB/ZPI
      TPIINVH= 0.5_JWRB*TPIINV
      TMP03 = 1.0_JWRB/(SDSBR*MICHE)
      SSDSC6M1=1._JWRB-SSDSC6
      NANGD=NANG/2

      DO M=1,NFRE

      ! COMPUTE SATURATION SPECTRUM
        DO IJ=KIJS,KIJL
          BTH0(IJ) = 0.0_JWRB
        ENDDO

        ! SATURATION TERM
        SSDSC2_SIG = SSDSC2*ZPIFR(M)
        ZCOEF = SSDSC2_SIG*SSDSC6
        ZCOEFM1 = SSDSC2_SIG*SSDSC6M1

        DO K=1,NANG
          DO IJ=KIJS,KIJL
            BTH(IJ) = 0.0_JWRB
          ENDDO

          ! integrates in directional sector
          DO K2=1,NSDSNTH*2+1
            KK=INDICESSAT(K,K2)
            DO IJ=KIJS,KIJL
              BTH(IJ) = BTH(IJ) + SATWEIGHTS(K,K2)*FL1(IJ,KK,M)
            ENDDO
          ENDDO
          DO IJ=KIJS,KIJL
            BTH(IJ)=BTH(IJ)*WAVNUM(IJ,M)*TPIINV*XK2CG(IJ,M)
            BTH0(IJ)=MAX(BTH0(IJ), BTH(IJ))
          ENDDO
        ENDDO

        DO K=1,NANG
          DO IJ=KIJS,KIJL
            BTH(IJ) = 0.0_JWRB
          ENDDO
          ! integrates in directional sector
          DO K2=1,NSDSNTH*2+1
            KK=INDICESSAT(K,K2)
            DO IJ=KIJS,KIJL
              BTH(IJ) = BTH(IJ) + SATWEIGHTS(K,K2)*FL1(IJ,KK,M)
            ENDDO
          ENDDO
          DO IJ=KIJS,KIJL
            BTH(IJ)=BTH(IJ)*WAVNUM(IJ,M)*TPIINV*XK2CG(IJ,M)
          ENDDO

          DO IJ=KIJS,KIJL
            D(IJ,K) = ZCOEF * (MAX(0._JWRB, BTH0(IJ)*TMP03-SSDSC4))**IPSAT + &
 &                      ZCOEFM1 * (MAX(0._JWRB, BTH(IJ)*TMP03-SSDSC4))**IPSAT

            ! cumulative term
            SCUMUL(IJ,K) = 0._JWRB
            IF(BTH0(IJ) > SDSBR)THEN
              SCUMUL(IJ,K) = MAX(SQRT(BTH(IJ)) - EPSR, 0._JWRB)**2
            ENDIF
          ENDDO

        ENDDO

        IF(M .GT. NDIKCUMUL)THEN
          ! CUMULATIVE TERM
          IF (SSDSC3 /= 0.0_JWRB) THEN
    
            DO K=1,NANG
            ! Correction of saturation level for shallow-water kinematics
            ! Cumulative effect based on lambda   (breaking probability is
            ! the expected rate of sweeping by larger breaking waves)
    
              DO K2=1,NANG
                KKD(K2)=ABS(K2-K)
                IF ( KKD(K2) > NANGD) KKD(K2)=KKD(K2)-NANGD
              ENDDO
    
              DO IJ=KIJS,KIJL
                RENEWALFREQ(IJ)=0.0_JWRB
              ENDDO
    
              DO M2=1,M-NDIKCUMUL
                DO K2=1,NANG
                  KK=KKD(K2)
                  DO IJ=KIJS,KIJL
                    ! Integrates over frequencies M2 and directions K2 to 
                    ! Integration is performed from M2=1 to a frequency lower than M: IK-NDIKCUMUL
                    RENEWALFREQ(IJ)=RENEWALFREQ(IJ)+ CUMULW(INDEP(IJ),KK,M2,M)*SCUMUL(IJ,K2)
                  ENDDO
                ENDDO
              ENDDO
    
              DO IJ=KIJS,KIJL
                D(IJ,K)= D(IJ,K) + RENEWALFREQ(IJ)
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!       WAVE-TURBULENCE INTERACTION TERM
        IF (SSDSC5 /= 0.0_JWRB) THEN
          TMP01 = 2._JWRB*SSDSC5/G
          DO IJ=KIJS,KIJL
            FACTURB(IJ) = TMP01*RAORW(IJ)*UFRIC(IJ)*UFRIC(IJ)
          ENDDO
            DO K=1,NANG
              DO IJ=KIJS,KIJL
                D(IJ,K)= D(IJ,K)- ZPIFR(M)*WAVNUM(IJ,M)*FACTURB(IJ)*COSWDIF(IJ,K)
              ENDDO
            ENDDO
        ENDIF


      ! ADD ALL CONTRIBUTIONS TO SOURCE TERM
        DO K=1, NANG
          DO IJ=KIJS,KIJL
            SL(IJ,K,M) = SL(IJ,K,M)+D(IJ,K)*FL1(IJ,K,M)
            FLD(IJ,K,M) = FLD(IJ,K,M)+D(IJ,K)
          ENDDO
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('SDISSIP_ARD',1,ZHOOK_HANDLE)

      END SUBROUTINE SDISSIP_ARD
