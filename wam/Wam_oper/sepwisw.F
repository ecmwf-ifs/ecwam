      SUBROUTINE SEPWISW (FL3, FL1, IJS, IJL, IG,THWOLD,USOLD,USNEW)

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

!*    PURPOSE.
!     --------

!       TO SEPARATE THE SWELL FROM THE WIND INTERACTING SEA

!**   INTERFACE.
!     ----------

!       *CALL* *SEPWISW (FL3, FL1, IJS, IJL, IG,THWOLD,USOLD,USNEW)*
!          *FL3* - BLOCK OF SPECTRA
!          *FL1* - SWELL FLAG ARRAY
!                  = 1 IF COMPONENT IS SWELL
!                  = 0 IF COMPONENT IS SEA
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *IG * - BLOCK NUMBER
!          *THWOLD*    INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
!                      WIND VELOCITY.
!          *USNEW*     NEW FRICTION VELOCITY IN M/S.
!          *USOLD*     INTERMEDIATE STORAGE OF MODULUS OF FRICTION
!                      VELOCITY.

!     METHOD.
!     -------

!       THE WAVES WHICH DO NOT INTERACT WITH THE WIND ARE
!       CONSIDERED SWELL.

!     EXTERNALS.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCOUT  , ONLY : IPFGTBL
      USE YOWFRED  , ONLY : FR       ,DFIM     ,C        ,
     &            TH       ,COSTH    ,SINTH    ,FRIC
      USE YOWMEAN  , ONLY : EMEAN
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO
      USE YOWSWEL  , ONLY : ESWELL   ,FSWELL   ,THSWELL  ,ESEA     ,
     &            FSEA     ,THWISEA  ,P1SWELL  ,P2SWELL  ,SPRDSWELL,
     &            P1SEA    ,P2SEA    ,SPRDSEA

! ----------------------------------------------------------------------

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE)  :: FL1,FL3
      REAL,DIMENSION(NINF:NSUP,NBLO) :: THWOLD,USOLD
      REAL,DIMENSION(NINF:NSUP) :: USNEW
      REAL, ALLOCATABLE :: THRESHOLD(:)

      LOGICAL :: LLEPSMIN, LLAK, LLPEAKF

! ----------------------------------------------------------------------


!*    1. THE SWELL DISTRIBUTION IS COMPUTED.
!        -----------------------------------

      ALLOCATE(THRESHOLD(IJS:IJL))
      DO K=1,NANG
        DO IJ=IJS,IJL
          THRESHOLD(IJ) = 1.2*USOLD(IJ,IG)*COS(TH(K)-THWOLD(IJ,IG))
        ENDDO
        DO M=1,NFRE
          CM = FRIC/C(M)
          DO IJ=IJS,IJL
            CHECKTA = CM*THRESHOLD(IJ)
            IF (CHECKTA.LT.1.) THEN
              FL1(IJ,K,M)=1.
            ELSE
              FL1(IJ,K,M)=0.
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(THRESHOLD)


!*    2. COMPUTATION OF SWELL OUTPUT PARAMETERS
!        ---------------------------------------

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            FL1(IJ,K,M) = FL3(IJ,K,M)*FL1(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO
  
      IF(IPFGTBL(16).NE.0) THEN
        LLAK=.FALSE.
        CALL FEMEAN(FL1, IJS, IJL, ESWELL(IJS), FSWELL(IJS),
     &              FSWELL(IJS), LLAK)
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



!*    3. COMPUTATION OF WIND SEA OUTPUT PARAMETERS
!        -----------------------------------------

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            FL1(IJ,K,M) = FL3(IJ,K,M)-FL1(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO
  
      IF(IPFGTBL(15).NE.0) THEN
        LLAK=.FALSE.
        CALL FEMEAN(FL1, IJS, IJL, ESEA(IJS), FSEA(IJS), FSEA(IJS),LLAK)
      ENDIF

!     femean must be called before semean since LLEPSMIN is false
      IF(IPFGTBL(11).NE.0 .OR. IPFGTBL(13).NE.0) THEN
        LLEPSMIN=.FALSE.
        CALL SEMEAN(FL1, IJS, IJL, ESEA(IJS), LLEPSMIN)
      ENDIF
      
      IF(IPFGTBL(13).NE.0) THEN
        CALL STHQ   (FL1, IJS, IJL, THWISEA(IJS))
!       if there isn't any windsea energy, set the windsea mean direction
!       to the wind direction.
        DO IJ=IJS,IJL
          IF(ESEA(IJ).LE.0.) THEN
            THWISEA(IJ)=THWOLD(IJ,IG)
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

      RETURN
      END SUBROUTINE SEPWISW
