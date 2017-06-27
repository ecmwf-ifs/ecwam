      SUBROUTINE SEPWISW (FL3, FL1, IJS, IJL, THWOLD, USOLD, USNEW,
     &                    U10NEW, THWNEW, LLPARTITION)

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

!       *CALL* *SEPWISW (FL3, FL1, IJS, IJL, THWOLD,USOLD,USNEW,
!                        U10NEW, THWNEW, LLPARTITION)*
!          *FL3* - BLOCK OF SPECTRA
!          *FL1* - SWELL FLAG ARRAY
!                  = 1 IF COMPONENT IS SWELL
!                  = 0 IF COMPONENT IS SEA
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *THWOLD*    INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
!                      WIND VELOCITY.
!          *USNEW*     NEW FRICTION VELOCITY IN M/S.
!          *USOLD*     INTERMEDIATE STORAGE OF MODULUS OF FRICTION
!                      VELOCITY.
!          *U10NEW* - LATESt WIND SPEED.
!          *THWNEW* - LATEST WIND DIRECTION.
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

      USE YOWCOUT  , ONLY : IPFGTBL
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFIMOFR  ,C        ,
     &            TH       ,FRIC     ,OLDWSFC
      USE YOWPCONS , ONLY : EPSMIN
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWMEAN  , ONLY : EMEAN
      USE YOWSWEL  , ONLY : ESWELL   ,FSWELL   ,THSWELL  ,ESEA     ,
     &            FSEA     ,THWISEA  ,P1SWELL  ,P2SWELL  ,SPRDSWELL,
     &            P1SEA    ,P2SEA    ,SPRDSEA
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! -----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS,IJL
      INTEGER :: IJ,K,M
      INTEGER, DIMENSION(IJS:IJL) :: MIJ

      REAL, DIMENSION(IJS:IJL), INTENT(IN) :: THWOLD, USOLD
      REAL, DIMENSION(IJS:IJL), INTENT(IN) :: USNEW, U10NEW, THWNEW
      REAL, DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: FL3
      REAL, DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: FL1
      REAL, DIMENSION(NFRE) :: CM
      REAL :: CHECKTA
      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION(IJS:IJL) :: THRESHOLD, ETNP
      REAL, DIMENSION(IJS:IJL) :: EMWS, FMWS, FPWS, EMAXWS, DUM
      REAL, DIMENSION(IJS:IJL,NFRE) :: RHOWGDFTH


      LOGICAL, INTENT(IN) :: LLPARTITION
      LOGICAL :: LLEPSMIN, LLAK, LLPEAKF

! ----------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SEPWISH',0,ZHOOK_HANDLE)
#endif

!*    1. THE SWELL DISTRIBUTION IS COMPUTED.
!        -----------------------------------

      DO IJ=IJS,IJL
        EMWS(IJ) = EPSMIN
        FMWS(IJ) = EPSMIN
        FPWS(IJ) = FR(NFRE) 
        EMAXWS(IJ) = EPSMIN
        DUM(IJ) = 0.
      ENDDO

      DO M=1,NFRE
        CM(M) = OLDWSFC*FRIC/C(M)
      ENDDO

      DO K=1,NANG
        DO IJ=IJS,IJL
          THRESHOLD(IJ) = USOLD(IJ)*COS(TH(K)-THWOLD(IJ))
        ENDDO
        DO M=1,NFRE
          DO IJ=IJS,IJL
            CHECKTA = CM(M)*THRESHOLD(IJ)
            IF (CHECKTA.LT.1.) THEN
              FL1(IJ,K,M)=1.
            ELSE
              FL1(IJ,K,M)=0.
              EMWS(IJ) = EMWS(IJ)+FL3(IJ,K,M)*DFIM(M)
              FMWS(IJ) = FMWS(IJ)+FL3(IJ,K,M)*DFIMOFR(M)
              IF(FL3(IJ,K,M).GT.EMAXWS(IJ)) THEN
                 EMAXWS(IJ) = FL3(IJ,K,M)
                 FPWS(IJ) = FR(M)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        FMWS(IJ) = EMWS(IJ)/FMWS(IJ)
      ENDDO


!*    2.1 COMPUTATION OF TOTAL SWELL OUTPUT PARAMETERS
!         --------------------------------------------

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            FL1(IJ,K,M) = MAX(FL3(IJ,K,M),EPSMIN)*FL1(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO

      IF(IPFGTBL(16).NE.0) THEN
        LLAK=.FALSE.
        CALL FEMEAN(FL1, IJS, IJL, ESWELL(IJS), FSWELL(IJS), LLAK)
      ENDIF

!     femean must be called before semean since LLEPSMIN is false
      IF(IPFGTBL(12).NE.0) THEN
        LLEPSMIN=.FALSE.
        CALL SEMEAN(FL1, IJS, IJL, ESWELL(IJS), LLEPSMIN)
        DO IJ=IJS,IJL
          IF (USNEW(IJ).EQ.0.) ESWELL(IJ) = EMEAN(IJ)
        ENDDO
      ENDIF
      
      IF(IPFGTBL(14).NE.0) THEN
        CALL STHQ   (FL1, IJS, IJL, THSWELL(IJS))
      ENDIF

      IF(IPFGTBL(29).NE.0) THEN
        IF(IPFGTBL(12).EQ.0) THEN
          LLEPSMIN=.FALSE.
          CALL SEMEAN(FL1, IJS, IJL, ESWELL(IJS), LLEPSMIN)
        ENDIF
        CALL MWP1 (FL1, IJS, IJL, ESWELL(IJS), P1SWELL(IJS))
      ENDIF

      IF(IPFGTBL(31).NE.0) THEN
        IF(IPFGTBL(12).EQ.0) THEN
          LLEPSMIN=.FALSE.
          CALL SEMEAN(FL1, IJS, IJL, ESWELL(IJS), LLEPSMIN)
        ENDIF
        CALL MWP2 (FL1, IJS, IJL, ESWELL(IJS), P2SWELL(IJS))
      ENDIF

      IF(IPFGTBL(33).NE.0) THEN
        LLPEAKF = .TRUE.
        CALL WDIRSPREAD (FL1, IJS, IJL, LLPEAKF, SPRDSWELL(IJS))
      ENDIF


!*    2.2 COMPUTATION OF THE PARTITIONED SWELL OUTPUT PARAMETERS
!         ------------------------------------------------------

      IF(LLPARTITION) THEN
        CALL FRCUTINDEX(IJS, IJL, FPWS, FMWS, USNEW, DUM, MIJ,RHOWGDFTH)
        CALL SEP3TR (FL1, IJS, IJL, MIJ, U10NEW, THWNEW, ETNP)
      ENDIF


!*    3. COMPUTATION OF WIND SEA OUTPUT PARAMETERS
!        -----------------------------------------

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            FL1(IJ,K,M) = MAX(FL3(IJ,K,M)-FL1(IJ,K,M)+EPSMIN,0.)
          ENDDO
        ENDDO
      ENDDO

      IF(IPFGTBL(15).NE.0) THEN
        LLAK=.FALSE.
        CALL FEMEAN(FL1, IJS, IJL, ESEA(IJS), FSEA(IJS), LLAK)
      ENDIF

!     femean must be called before semean since LLEPSMIN is false
      IF(IPFGTBL(11).NE.0 .OR. IPFGTBL(13).NE.0 ) THEN
        LLEPSMIN=.FALSE.
        CALL SEMEAN(FL1, IJS, IJL, ESEA(IJS), LLEPSMIN)
      ENDIF

      IF(IPFGTBL(13).NE.0) THEN
        CALL STHQ   (FL1, IJS, IJL, THWISEA(IJS))
!       if there isn't any windsea energy, set the windsea mean direction
!       to the wind direction.
        DO IJ=IJS,IJL
          IF(ESEA(IJ).LE.0.) THEN
            THWISEA(IJ)=THWOLD(IJ)
          ENDIF
        ENDDO
      ENDIF

      IF(IPFGTBL(28).NE.0) THEN
        IF(IPFGTBL(11).EQ.0) THEN
          LLEPSMIN=.FALSE.
          CALL SEMEAN(FL1, IJS, IJL, ESEA(IJS), LLEPSMIN)
        ENDIF
        CALL MWP1 (FL1, IJS, IJL, ESEA(IJS), P1SEA(IJS))
      ENDIF

      IF(IPFGTBL(30).NE.0) THEN
        IF(IPFGTBL(11).EQ.0) THEN
          LLEPSMIN=.FALSE.
          CALL SEMEAN(FL1, IJS, IJL, ESEA(IJS), LLEPSMIN)
        ENDIF
        CALL MWP2 (FL1, IJS, IJL, ESEA(IJS), P2SEA(IJS))
      ENDIF

      IF(IPFGTBL(32).NE.0) THEN
        LLPEAKF = .TRUE.
        CALL WDIRSPREAD (FL1, IJS, IJL, LLPEAKF, SPRDSEA(IJS))
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

      RETURN
      END SUBROUTINE SEPWISW
