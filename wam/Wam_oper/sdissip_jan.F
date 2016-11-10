      SUBROUTINE SDISSIP_JAN (F, FL, IJS, IJL, SL,
     &                        EMEAN, F1MEAN, XKMEAN,
     &                        PHIEPS, TAUWD, MIJ, LCFLX)

! ----------------------------------------------------------------------

!**** *SDISSIP* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.

!     S.D.HASSELMANN.
!     MODIFIED TO SHALLOW WATER : G. KOMEN , P. JANSSEN
!     OPTIMIZATION : L. ZAMBRESKY
!     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
!     J. BIDLOT   ECMWF  NOVEMBER 2004  REFORMULATION BASED ON XKMEAN
!                                       AND F1MEAN.

!*    PURPOSE.
!     --------
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SDISSIP (F, FL, IJS, IJL, SL,*
!                        EMEAN, F1MEAN, XKMEAN,*
!                        PHIEPS, TAUWD, MIJ, LCFLX)
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY
!          *EMEAN* - MEAN ENERGY DENSITY 
!          *F1MEAN* - MEAN FREQUENCY BASED ON 1st MOMENT.
!          *XKMEAN* - MEAN WAVE NUMBER BASED ON 1st MOMENT.
!          *PHIEPS* - ENERGY FLUX FROM WAVES TO OCEAN INTEGRATED OVER 
!                     THE PROGNOSTIC RANGE.
!          *TAUWD*  - DISSIPATION STRESS INTEGRATED OVER
!                     THE PROGNOSTIC RANGE.
!          *MIJ*    - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *LCFLX*  - TRUE IF THE CALCULATION FOR THE FLUXES ARE 
!                     PERFORMED.


!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       G.KOMEN, S. HASSELMANN AND K. HASSELMANN, ON THE EXISTENCE
!          OF A FULLY DEVELOPED WINDSEA SPECTRUM, JGR, 1984.

! ----------------------------------------------------------------------

      USE YOWFRED  , ONLY : FR       ,DELTH    ,DFIM     ,FRATIO
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER
      USE YOWSHAL  , ONLY : DEPTH    ,TFAK     ,INDEP
      USE YOWSTAT  , ONLY : ISHALLO  ,LBIWBK
      USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL
      INTEGER, DIMENSION(IJS:IJL), INTENT(IN) :: MIJ
      INTEGER :: IJ, K, M

      REAL,DIMENSION(IJS:IJL), INTENT(IN) :: EMEAN, F1MEAN, XKMEAN
      REAL,DIMENSION(IJS:IJL), INTENT(INOUT) :: PHIEPS, TAUWD
      REAL,DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F
      REAL,DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: FL, SL

      REAL :: SCDFM, ROG, CONSD, X, SDISS
      REAL :: ZHOOK_HANDLE
      REAL :: CONSS
      REAL,DIMENSION(IJS:IJL) :: CM, TEMP1, SDS
      REAL, DIMENSION(NFRE) :: FAC
      REAL, DIMENSION(IJS:IJL,NFRE) :: CONSTFM

      LOGICAL, INTENT(IN) :: LCFLX

      REAL, PARAMETER :: CDIS = 1.33
      REAL, PARAMETER :: DELTA = 0.5

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SDISSIP_JAN',0,ZHOOK_HANDLE)
#endif

!*    1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE
!*       FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.
!        --------------------------------------------------------------

      ROG = ROWATER*G

      IF (LCFLX) THEN
        PHIEPS(IJS:IJL) = 0.
        TAUWD(IJS:IJL) = 0.
!       !!!! CONSTFM is only defined up to M=MIJ(IJ)
        SCDFM = 0.5*DELTH*(1.-1./FRATIO)
        DO IJ=IJS,IJL
          DO M=1,MIJ(IJ)-1
            CONSTFM(IJ,M) = ROG*DFIM(M)
          ENDDO
          CONSTFM(IJ,MIJ(IJ)) = ROG*SCDFM*FR(MIJ(IJ))
        ENDDO
      ENDIF

      IF (ISHALLO.EQ.1) THEN
        CONSD = -CDIS*ZPI**9/G**4
        DO IJ=IJS,IJL
          SDS(IJ)=CONSD*F1MEAN(IJ)*EMEAN(IJ)**2*F1MEAN(IJ)**8
        ENDDO
        DO M=1,NFRE
          FAC(M) = ZPI*FR(M)
          CM(IJS)  = FAC(M)/G
          DO IJ=IJS,IJL
            X         = (FR(M)/F1MEAN(IJ))**2
            TEMP1(IJ) = SDS(IJ)*( (1.-DELTA)*X + DELTA*X**2)
          ENDDO
          DO K=1,NANG
            DO IJ=IJS,IJL
              SDISS = TEMP1(IJ)*F(IJ,K,M)
              SL(IJ,K,M) = SL(IJ,K,M)+TEMP1(IJ)*F(IJ,K,M)
              FL(IJ,K,M) = FL(IJ,K,M)+TEMP1(IJ)

              IF (LCFLX.AND.M.LE.MIJ(IJ)) THEN
                PHIEPS(IJ) = PHIEPS(IJ)+SDISS*CONSTFM(IJ,M)
                TAUWD(IJ)  = TAUWD(IJ)+CM(IJS)*SDISS*CONSTFM(IJ,M)
              ENDIF

            ENDDO
          ENDDO     
        ENDDO

      ELSE

        CONSS = -CDIS*ZPI
        DO IJ=IJS,IJL
          SDS(IJ)=CONSS*F1MEAN(IJ)*EMEAN(IJ)**2*XKMEAN(IJ)**4
        ENDDO

        DO M=1,NFRE
          FAC(M) = ZPI*FR(M)
          DO IJ=IJS,IJL
            X         = TFAK(INDEP(IJ),M)/XKMEAN(IJ)
            TEMP1(IJ) = SDS(IJ)*( (1.-DELTA)*X + DELTA*X**2)
            CM(IJ)    = TFAK(INDEP(IJ),M)/FAC(M)
          ENDDO
          DO K=1,NANG
            DO IJ=IJS,IJL
              SDISS = TEMP1(IJ)*F(IJ,K,M)
              SL(IJ,K,M) = SL(IJ,K,M)+TEMP1(IJ)*F(IJ,K,M)
              FL(IJ,K,M) = FL(IJ,K,M)+TEMP1(IJ)

              IF (LCFLX.AND.M.LE.MIJ(IJ)) THEN
                PHIEPS(IJ) = PHIEPS(IJ)+SDISS*CONSTFM(IJ,M)
                TAUWD(IJ)  = TAUWD(IJ)+CM(IJ)*SDISS*CONSTFM(IJ,M)
              ENDIF

            ENDDO
          ENDDO
        ENDDO

      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SDISSIP_JAN',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE SDISSIP_JAN
