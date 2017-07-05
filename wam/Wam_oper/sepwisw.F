      SUBROUTINE SEPWISW (IJS, IJL, MIJ, FL3, XLLWS,
     &                    USNEW, U10NEW, THWNEW,
     &                    ESWELL   ,FSWELL   ,THSWELL  ,
     &                    P1SWELL  ,P2SWELL  ,SPRDSWELL,
     &                    ESEA     ,FSEA     ,THWISEA  ,
     &                    P1SEA    ,P2SEA    ,SPRDSEA  ,
     &                    EMTRAIN  ,THTRAIN  ,PMTRAIN,
     &                    LLPARTITION)

! ----------------------------------------------------------------------

!**** *SEPWISW* - COMPUTES THE SWELL ENERGY, THE MEAN SWELL DIRECTION,
!****             THE MEAN SWELL FREQUENCY AND THE MEAN WINDSEA DIR.

!     P.LIONELLO     FEBRUARY 87

!     L.ZAMBRESKY    NOVEMBER 87   GKSS/ECMWF   OPTIMIZED SUB.
!     J. BIDLOT      FEBRARY 1996       ECMWF   MESSAGE PASSING
!     J. BIDLOT      MAY     1999       ECMWF   ADD HIGH FREQUENCY TAIL
!     J. BIDLOT      APRIL   2000       ECMWF   ADD  EXTRA PARAMETERS
!     J. BIDLOT      DECEMBER2003       ECMWF   MOVE ALL ALLOCATION TO
!                                               *OUTBS*
!                                               TO MAKE IT CALLABLE IN AN
!                                               OPENMP LOOP.
!     J. BIDLOT     JUNE 213           ECMWF    ADD PARTITIONING SCHEME

!*    PURPOSE.
!     --------

!       TO SEPARATE THE SWELL FROM THE WIND INTERACTING SEA

!**   INTERFACE.
!     ----------

!       *CALL* SEPWISW (IJS, IJL, MIJ, FL3, XLLWS,
!    &                  USNEW, U10NEW, THWNEW,
!    &                  ESWELL   ,FSWELL   ,THSWELL  ,
!    &                  P1SWELL  ,P2SWELL  ,SPRDSWELL,
!    &                  ESEA     ,FSEA     ,THWISEA  ,
!    &                  P1SEA    ,P2SEA    ,SPRDSEA  ,
!    &                  EMTRAIN  ,THTRAIN  ,PMTRAIN,
!    &                  LLPARTITION)
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *MIJ* - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *FL3* - BLOCK OF SPECTRA
!          *XLLWS* - WINDSEA MASK FROM INPUT SOURCE TERM
!          *USNEW* - NEW FRICTION VELOCITY IN M/S.
!          *U10NEW* - LATESt WIND SPEED.
!          *THWNEW* - LATEST WIND DIRECTION.
!          *        - SWELL and WINDSEA PARAMETERS.
!          *LLPARTITION* TRUE IF SPECTRAL PARTITIONING IS USED.

!     METHOD.
!     -------

!       THE WAVES WHICH DO NOT INTERACT WITH THE WIND ARE
!       CONSIDERED SWELL.
!       IF LLPARTITION THEN A SPECTRAl PARTITIONING SCHEME WILL BE USED
!       TO SPLIT THE SWELL SPECTRUM INTO ITS TWO MAIN COMPONENTS.

!     EXTERNALS.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCOUT  , ONLY : IPFGTBL  ,NTRAIN
      USE YOWFRED  , ONLY : C        ,TH       ,FRIC     ,OLDWSFC
      USE YOWPCONS , ONLY : EPSMIN
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! -----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS,IJL
      INTEGER, DIMENSION(IJS:IJL), INTENT(IN) :: MIJ
      REAL, DIMENSION(IJS:IJL), INTENT(IN) :: USNEW, U10NEW, THWNEW
      REAL, DIMENSION(IJS:IJL), INTENT(OUT) :: ESWELL ,FSWELL ,THSWELL
      REAL, DIMENSION(IJS:IJL), INTENT(OUT) :: P1SWELL,P2SWELL,SPRDSWELL
      REAL, DIMENSION(IJS:IJL), INTENT(OUT) :: ESEA   ,FSEA   ,THWISEA
      REAL, DIMENSION(IJS:IJL), INTENT(OUT) :: P1SEA  ,P2SEA  ,SPRDSEA
      REAL, DIMENSION(IJS:IJL,NTRAIN), INTENT(OUT) :: EMTRAIN
      REAL, DIMENSION(IJS:IJL,NTRAIN), INTENT(OUT) :: THTRAIN, PMTRAIN
      REAL, DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: FL3
      REAL, DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: XLLWS 
      LOGICAL, INTENT(IN) :: LLPARTITION

      INTEGER :: IJ,K,M

      REAL :: COEF
      REAL :: CHECKTA
      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION(NFRE) :: CM
      REAL, DIMENSION(IJS:IJL) :: ETNP, R
      REAL, DIMENSION(IJS:IJL,NFRE) :: XINVWVAGE
      REAL, DIMENSION(IJS:IJL,NANG) :: DIRCOEF
      REAL, DIMENSION(IJS:IJL,NANG,NFRE) :: SWM, FL1

      LOGICAL :: LLPEAKF

! ----------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SEPWISH',0,ZHOOK_HANDLE)
#endif

      ESWELL(:)=0.0
      ESEA(:)=0.0
      EMTRAIN(:,:)=0.0

!*    1. THE SWELL DISTRIBUTION IS COMPUTED.
!        -----------------------------------

      COEF = OLDWSFC*FRIC

      DO M=1,NFRE
        CM(M)=1.0/C(M)
        DO IJ=IJS,IJL
          XINVWVAGE(IJ,M)=USNEW(IJ)*CM(M)
        ENDDO
      ENDDO

      DO K=1,NANG
        DO IJ=IJS,IJL
          DIRCOEF(IJ,K)=COEF*COS(TH(K)-THWNEW(IJ))
        ENDDO
      ENDDO

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            IF (XLLWS(IJ,K,M).NE.0.) THEN
              ! this is windsea 
              SWM(IJ,K,M)=0.
            ELSE
              CHECKTA=XINVWVAGE(IJ,M)*DIRCOEF(IJ,K)
              IF (CHECKTA.GE.1.) THEN
                ! this is extra windsea 
                SWM(IJ,K,M)=0.
              ELSE
                ! this is swell
                SWM(IJ,K,M)=1.
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

!     CHECK THAT TOTAL SWELL MEAN FREQUENCY IS LOWER THAN WINDSEA ONE
!     OTHERWISE RESET WIND SECTOR TO WINDSEA

!     SWELL:
      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            FL1(IJ,K,M)=FL3(IJ,K,M)*SWM(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO
      CALL FEMEAN(FL1, IJS, IJL, ESWELL, FSWELL)

!     WINDSEA:
      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            FL1(IJ,K,M)=MAX(FL3(IJ,K,M)-FL1(IJ,K,M),0.)
          ENDDO
        ENDDO
      ENDDO
      CALL FEMEAN(FL1, IJS, IJL, ESEA, FSEA)

!     CHECK:
      DO IJ=IJS,IJL
         IF(FSWELL(IJ).GT.FSEA(IJ)) THEN
           R(IJ)=1.
         ELSE
           R(IJ)=0.
         ENDIF
      ENDDO
      DO K=1,NANG
        DO IJ=IJS,IJL
          ! add 0.25 factor to extend windsea area
          DIRCOEF(IJ,K)=R(IJ)*COEF*SIGN(1.,0.4+COS(TH(K)-THWNEW(IJ)))
        ENDDO
      ENDDO
      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            CHECKTA = XINVWVAGE(IJ,M)*DIRCOEF(IJ,K)
            IF (CHECKTA.GE.1.) THEN
              ! this is additional windsea 
              SWM(IJ,K,M)=0.
            ENDIF
          ENDDO
        ENDDO
      ENDDO

!     ADJUST THE LOW FREQUENCY BOUNDARY OF THE WINDSEAS AREA 
!     TO TOPOLOGICALLY CONNECT IT TO THE WINDSEA
      DO IJ=IJS,IJL
        DO K=1,NANG
          DO M=NFRE,2,-1
            IF(SWM(IJ,K,M).EQ.1. .AND. SWM(IJ,K,M-1).EQ.1.) THEN
              EXIT
            ELSEIF(SWM(IJ,K,M).EQ.0. .AND. SWM(IJ,K,M-1).EQ.1.) THEN
               IF(FL3(IJ,K,M).GE.FL3(IJ,K,M-1)) SWM(IJ,K,M-1)=0. 
            ENDIF
          ENDDO
        ENDDO
      ENDDO

!*    2.1 COMPUTATION OF TOTAL SWELL OUTPUT PARAMETERS
!         --------------------------------------------

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            FL1(IJ,K,M)=MAX(FL3(IJ,K,M),EPSMIN)*SWM(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO

      IF(IPFGTBL(12).NE.0 .OR. IPFGTBL(16).NE.0 .OR.
     &   IPFGTBL(29).NE.0 .OR. IPFGTBL(31).NE.0 .OR.
     &   LLPARTITION ) THEN
        CALL FEMEAN(FL1, IJS, IJL, ESWELL, FSWELL)
      ENDIF
      
      IF(IPFGTBL(14).NE.0 .OR. LLPARTITION) THEN
        CALL STHQ(FL1, IJS, IJL, THSWELL)
      ENDIF

      IF(IPFGTBL(29).NE.0) THEN
        CALL MWP1(FL1, IJS, IJL, ESWELL, P1SWELL)
      ENDIF

      IF(IPFGTBL(31).NE.0) THEN
        CALL MWP2(FL1, IJS, IJL, ESWELL, P2SWELL)
      ENDIF

      IF(IPFGTBL(33).NE.0) THEN
        LLPEAKF = .TRUE.
        CALL WDIRSPREAD (FL1, IJS, IJL, LLPEAKF, SPRDSWELL)
      ENDIF


!*    2.2 COMPUTATION OF THE PARTITIONED SWELL OUTPUT PARAMETERS
!         ------------------------------------------------------

      IF(LLPARTITION) THEN
        CALL SEP3TR (FL1, IJS, IJL, MIJ, U10NEW, THWNEW,
     &               ESWELL, FSWELL, THSWELL, FSEA,
     &               ETNP,
     &               EMTRAIN  ,THTRAIN  ,PMTRAIN)
      ENDIF


!*    3. COMPUTATION OF WIND SEA OUTPUT PARAMETERS
!        -----------------------------------------

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            FL1(IJ,K,M)=MAX(FL3(IJ,K,M)-FL1(IJ,K,M)+EPSMIN,0.)
          ENDDO
        ENDDO
      ENDDO

      IF(IPFGTBL(15).NE.0 .OR. IPFGTBL(11).NE.0 .OR.
     &   IPFGTBL(13).NE.0 .OR. IPFGTBL(28).NE.0 .OR.
     &   IPFGTBL(30).NE.0) THEN
        CALL FEMEAN(FL1, IJS, IJL, ESEA, FSEA)
      ENDIF

      IF(IPFGTBL(13).NE.0) THEN
        CALL STHQ(FL1, IJS, IJL, THWISEA)
!       if there isn't any windsea energy, set the windsea mean direction
!       to the wind direction.
        DO IJ=IJS,IJL
          IF(ESEA(IJ).LE.0.) THEN
            THWISEA(IJ)=THWNEW(IJ)
          ENDIF
        ENDDO
      ENDIF

      IF(IPFGTBL(28).NE.0) THEN
        CALL MWP1(FL1, IJS, IJL, ESEA, P1SEA)
      ENDIF

      IF(IPFGTBL(30).NE.0) THEN
        CALL MWP2(FL1, IJS, IJL, ESEA, P2SEA)
      ENDIF

      IF(IPFGTBL(32).NE.0) THEN
        LLPEAKF = .TRUE.
        CALL WDIRSPREAD(FL1, IJS, IJL, LLPEAKF, SPRDSEA)
      ENDIF


!     RE-ASSIGN THE NON PARTITONED ENERGY TO WINDSEA
!!!!!! it does not work as it creates holes in the total swell field
!!!!!! as feature that cannot be handled by our users
!     AND REMOVE IT FROM THE SWELL
!      IF(LLPARTITION) THEN
!        IF(IPFGTBL(11).NE.0 ) THEN
!          DO IJ=IJS,IJL
!            ESEA(IJ) = ESEA(IJ)+ETNP(IJ)
!          ENDDO
!        ENDIF
!        IF(IPFGTBL(12).NE.0) THEN
!          DO IJ=IJS,IJL
!            ESWELL(IJ) = MAX(ESWELL(IJ)-ETNP(IJ),0.0)
!          ENDDO
!        ENDIF
!      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SEPWISH',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE SEPWISW
