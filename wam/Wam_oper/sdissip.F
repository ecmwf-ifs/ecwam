      SUBROUTINE SDISSIP (F, FL, IJS, IJL, SL, F1MEAN, XKMEAN)

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

!       *CALL* *SDISSIP (F, FL, IJS, IJL, SL, F1MEAN, XKMEAN)*
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY
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

! ----------------------------------------------------------------------

      USE YOWFRED  , ONLY : FR
      USE YOWMEAN  , ONLY : EMEAN
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOWSHAL  , ONLY : TFAK     ,INDEP
      USE YOWSTAT  , ONLY : ISHALLO

! ----------------------------------------------------------------------

      DIMENSION TEMP1(IJS:IJL), SDS(IJS:IJL)

! ----------------------------------------------------------------------

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F,FL,SL
      REAL,DIMENSION(IJS:IJL) :: F1MEAN, XKMEAN

      PARAMETER (CDIS = 2.1)
      PARAMETER (DELTA = 0.6)

! ----------------------------------------------------------------------

!*    1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE
!*       FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.
!        --------------------------------------------------------------

      IF (ISHALLO.EQ.1) THEN
        CONSD = -CDIS*ZPI**9/G**4
        DO IJ=IJS,IJL
          SDS(IJ) = CONSD*F1MEAN(IJ)*EMEAN(IJ)**2*F1MEAN(IJ)**8
        ENDDO
        DO M=1,NFRE
          DO IJ=IJS,IJL
            X         = (FR(M)/F1MEAN(IJ))**2
            TEMP1(IJ) = SDS(IJ)*( (1.-DELTA)*X + DELTA*X**2)
          ENDDO
          DO K=1,NANG
            DO IJ=IJS,IJL
              SL(IJ,K,M) = SL(IJ,K,M)+TEMP1(IJ)*F(IJ,K,M)
              FL(IJ,K,M) = FL(IJ,K,M)+TEMP1(IJ)
            ENDDO
          ENDDO
        ENDDO
      ELSE
!SHALLOW
        CONSS = -CDIS*ZPI
        DO IJ=IJS,IJL
          SDS(IJ) = CONSS*F1MEAN(IJ)*EMEAN(IJ)**2*XKMEAN(IJ)**4
        ENDDO

        DO M=1,NFRE
          DO IJ=IJS,IJL
            X         = TFAK(INDEP(IJ),M)/XKMEAN(IJ)
            TEMP1(IJ) = SDS(IJ)*( (1.-DELTA)*X + DELTA*X**2)
          ENDDO
          DO K=1,NANG
            DO IJ=IJS,IJL
              SL(IJ,K,M) = SL(IJ,K,M)+TEMP1(IJ)*F(IJ,K,M)
              FL(IJ,K,M) = FL(IJ,K,M)+TEMP1(IJ)
            ENDDO
          ENDDO
        ENDDO
!SHALLOW
      ENDIF

      RETURN
      END SUBROUTINE SDISSIP
