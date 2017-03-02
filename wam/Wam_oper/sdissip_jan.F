      SUBROUTINE SDISSIP_JAN (F, FL, SL, IJS, IJL,
     &                    EMEAN, F1MEAN, XKMEAN)

! ----------------------------------------------------------------------

!**** *SDISSIP_JAN* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.

!     S.D.HASSELMANN.
!     MODIFIED TO SHALLOW WATER : G. KOMEN , P. JANSSEN
!     OPTIMIZATION : L. ZAMBRESKY
!     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
!     J. BIDLOT   ECMWF  NOVEMBER 2004  REFORMULATION BASED ON XKMEAN
!                                       AND F1MEAN.
!     P. JANSSEN  ECMWF  JANUARY 2006   ADD BOTTOM-INDUCED DISSIPATION.

!*    PURPOSE.
!     --------
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SDISSIP_JAN (F, FL, IJS, IJL, SL, EMEAN,F1MEAN, XKMEAN,)*
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
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

      USE YOWFRED  , ONLY : FR       ,DELTH    ,DFIM     ,FRATIO
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOWSHAL  , ONLY : DEPTH    ,CINV     ,TFAK     ,INDEP
      USE YOWSTAT  , ONLY : ISHALLO
      USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL

      REAL,DIMENSION(IJS:IJL), INTENT(IN) :: EMEAN, F1MEAN, XKMEAN
      REAL,DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F
      REAL,DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: FL, SL

      INTEGER :: IJ, K, M

      REAL :: SCDFM, CONSD, CONSS, DELTAM1
      REAL :: ZHOOK_HANDLE
      REAL,DIMENSION(IJS:IJL) :: CM, TEMP1, SDS, X


      REAL, PARAMETER :: CDIS = 1.33
      REAL, PARAMETER :: DELTA = 0.5

! ----------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SDISSIP_JAN',0,ZHOOK_HANDLE)
#endif

!*    1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE
!*       FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.
!        --------------------------------------------------------------

      DELTAM1=1.-DELTA

        IF (ITEST.GE.2) THEN
          WRITE(IU06,*) '   SUB. SDISSIP_JAN: START DO-LOOP (ISHALLO=0)'
          CALL FLUSH (IU06)
        ENDIF

      IF (ISHALLO.EQ.1) THEN
!       DEEP
        CONSD = -CDIS*ZPI**9/G**4
        DO IJ=IJS,IJL
          SDS(IJ)=CONSD*F1MEAN(IJ)*EMEAN(IJ)**2*F1MEAN(IJ)**8
        ENDDO
      ELSE
!       SHALLOW
        CONSS = -CDIS*ZPI
       DO IJ=IJS,IJL
          SDS(IJ)=CONSS*F1MEAN(IJ)*EMEAN(IJ)**2*XKMEAN(IJ)**4
        ENDDO
      ENDIF

      DO M=1,NFRE
!       DEEP
        IF (ISHALLO.EQ.1) THEN
          DO IJ=IJS,IJL
            X(IJ) = (FR(M)/F1MEAN(IJ))**2
          ENDDO
        ELSE
!         SHALLOW
          DO IJ=IJS,IJL
            X(IJ) = TFAK(INDEP(IJ),M)/XKMEAN(IJ)
          ENDDO
        ENDIF

        DO IJ=IJS,IJL
          TEMP1(IJ) = SDS(IJ)*X(IJ)*(DELTAM1 + DELTA*X(IJ))
        ENDDO

        DO K=1,NANG
          DO IJ=IJS,IJL
            FL(IJ,K,M) = FL(IJ,K,M) + TEMP1(IJ)
            SL(IJ,K,M) = SL(IJ,K,M) + TEMP1(IJ)*F(IJ,K,M)
          ENDDO
        ENDDO

      ENDDO


#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SDISSIP_JAN',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE SDISSIP_JAN
