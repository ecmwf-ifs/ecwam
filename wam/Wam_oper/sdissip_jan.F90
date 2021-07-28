      SUBROUTINE SDISSIP_JAN (GFL, FLD, SL, IJS, KIJL, KIJS, KIJL,  &
     &                        EMEAN, F1MEAN, XKMEAN)

! ----------------------------------------------------------------------

!**** *SDISSIP_JAN* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.

!     S.D.HASSELMANN.
!     MODIFIED TO SHALLOW WATER : G. KOMEN , P. JANSSEN
!     OPTIMIZATION : L. ZAMBRESKY
!     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
!     J. BIDLOT   ECMWF  NOVEMBER 2004  REFORMULATION BASED ON XKMEAN
!                                       AND F1MEAN.
!                        AUGUST 2020 Added small viscous dissipation term

!*    PURPOSE.
!     --------
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SDISSIP_JAN (GFL, FLD, IJS, IJL, KIJS, KIJL, SL,
!                            EMEAN,F1MEAN, XKMEAN,)*
!          *GFL*   - SPECTRUM.
!          *FLD*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY
!          IJS:IJL - 1st DIMENSION OF GFL
!          *KIJS* - INDEX OF FIRST GRIDPOINT
!          *KIJL* - INDEX OF LAST GRIDPOINT
!          *EMEAN* - MEAN ENERGU DENSITY 
!          *F1MEAN* - MEAN FREQUENCY BASED ON 1st MOMENT.
!          *XKMEAN* - MEAN WAVE NUMBER BASED ON 1st MOMENT.


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

! ---------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DELTH    ,DFIM     ,FRATIO
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,ZPI4GM2
      USE YOWPHYS  , ONLY : CDIS     ,DELTA_SDIS, RNU    ,CDISVIS
      USE YOWSHAL  , ONLY : TFAK     ,INDEP
      USE YOWSTAT  , ONLY : ISHALLO
      USE YOWTEST  , ONLY : IU06     ,ITEST

      USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: GFL

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT):: FLD, SL

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN):: EMEAN, F1MEAN, XKMEAN

      INTEGER(KIND=JWIM) :: IJ, K, M

      REAL(KIND=JWRB) :: SCDFM, CONSD, CONSS, DELTA_SDISM1, CVIS
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL) :: CM, TEMP1, SDS, X
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL) :: XK2

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SDISSIP_JAN',0,ZHOOK_HANDLE)

!*    1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE
!*       FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.
!        --------------------------------------------------------------

      DELTA_SDISM1=1.0_JWRB-DELTA_SDIS

        IF (ITEST.GE.2) THEN
          WRITE(IU06,*) '   SUB. SDISSIP_JAN: START DO-LOOP (ISHALLO=0)'
          CALL FLUSH (IU06)
        ENDIF

      IF (ISHALLO.EQ.1) THEN
!       DEEP
        CONSD = CDIS*ZPI**9/G**4
        DO IJ=KIJS,KIJL
          SDS(IJ)=CONSD*F1MEAN(IJ)*EMEAN(IJ)**2*F1MEAN(IJ)**8
        ENDDO
      ELSE
!       SHALLOW
        CONSS = CDIS*ZPI
       DO IJ=KIJS,KIJL
          SDS(IJ)=CONSS*F1MEAN(IJ)*EMEAN(IJ)**2*XKMEAN(IJ)**4
        ENDDO
      ENDIF

      DO M=1,NFRE
!       DEEP
        IF (ISHALLO.EQ.1) THEN
          DO IJ=KIJS,KIJL
            X(IJ) = (FR(M)/F1MEAN(IJ))**2
            XK2(IJ) = ZPI4GM2 * FR(M)**4 
          ENDDO
        ELSE
!         SHALLOW
          DO IJ=KIJS,KIJL
            X(IJ) = TFAK(INDEP(IJ),M)/XKMEAN(IJ)
            XK2(IJ) = TFAK(INDEP(IJ),M)**2
          ENDDO
        ENDIF

        CVIS=RNU*CDISVIS
        DO IJ=KIJS,KIJL
          TEMP1(IJ) = SDS(IJ)*X(IJ)*(DELTA_SDISM1 + DELTA_SDIS*X(IJ)) + CVIS*XK2(IJ)
        ENDDO

        DO K=1,NANG
          DO IJ=KIJS,KIJL
            FLD(IJ,K,M) = FLD(IJ,K,M) + TEMP1(IJ)
            SL(IJ,K,M) = SL(IJ,K,M) + TEMP1(IJ)*GFL(IJ,K,M)
          ENDDO
        ENDDO

      ENDDO

      IF (LHOOK) CALL DR_HOOK('SDISSIP_JAN',1,ZHOOK_HANDLE)

      END SUBROUTINE SDISSIP_JAN
