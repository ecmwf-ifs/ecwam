      SUBROUTINE SDISSIP (F, FL, IJS, IJL, SL)

! ----------------------------------------------------------------------

!**** *SDISSIP* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.

!     S.D.HASSELMANN.
!     MODIFIED TO SHALLOW WATER : G. KOMEN , P. JANSSEN
!     OPTIMIZATION : L. ZAMBRESKY
!     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL

!*    PURPOSE.
!     --------
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SDISSIP (F, FL, IJS, IJL, SL)*
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY

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
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,AKMEAN
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOWSHAL  , ONLY : TFAK     ,INDEP
      USE YOWSTAT  , ONLY : ISHALLO

! ----------------------------------------------------------------------

      DIMENSION TEMP1(IJS:IJL), SDS(IJS:IJL)

! ----------------------------------------------------------------------

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F,FL,SL
      PARAMETER (CDIS = 4.5)

! ----------------------------------------------------------------------

!*    1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE
!*       FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.
!        --------------------------------------------------------------

      CONSD = -0.5*CDIS*ZPI**9/G**4
      CONSS = -0.5*CDIS*ZPI

      IF (ISHALLO.EQ.1) THEN
        DO IJ=IJS,IJL
          SDS(IJ) = CONSD*EMEAN(IJ)**2*FMEAN(IJ)**9
        ENDDO

        DO M=1,NFRE
          DO IJ=IJS,IJL
            X         = (FR(M)/FMEAN(IJ))**2
            TEMP1(IJ) = SDS(IJ)*( X + X**2)
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
        DO IJ=IJS,IJL
          SDS(IJ) = CONSS*FMEAN(IJ)*EMEAN(IJ)**2*AKMEAN(IJ)**4
        ENDDO

        DO M=1,NFRE
          DO IJ=IJS,IJL
            X         = TFAK(INDEP(IJ),M)/AKMEAN(IJ)
            TEMP1(IJ) = SDS(IJ)*( X + X**2)
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
