      SUBROUTINE SEPWISW (FL3, FL1, IJS, IJL, IG,THWOLD,USOLD,USNEW)

! ----------------------------------------------------------------------

!**** *SEPWISW* - COMPUTES THE SWELL ENERGY, THE MEAN SWELL DIRECTION,
!****             THE MEAN SWELL FREQUENCY AND THE MEAN WINDSEA DIR.

!     P.LIONELLO     FEBRUARY 87

!     L.ZAMBRESKY    NOVEMBER 87   GKSS/ECMWF   OPTIMIZED SUB.
!     J. BIDLOT      FEBRARY 1996       ECMWF   MESSAGE PASSING
!     J. BIDLOT      MAY     1999       ECMWF   ADD HIGH FREQUENCY TAIL
!     J. BIDLOT      APRIL   200        ECMWF   ADD  EXTRA PARAMETERS

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
     &            TH       ,COSTH    ,SINTH
      USE YOWMEAN  , ONLY : EMEAN
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO
      USE YOWSWEL  , ONLY : ESWELL   ,FSWELL   ,THSWELL  ,ESEA     ,
     &            FSEA     ,THWISEA  ,P1SWELL  ,P2SWELL  ,SPRDSWELL,
     &            P1SEA    ,P2SEA    ,SPRDSEA

! ----------------------------------------------------------------------

      REAL, PARAMETER :: FRIC=28.

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE)  :: FL1,FL3
      REAL,DIMENSION(NINF:NSUP,NBLO) :: THWOLD,USOLD
      REAL,DIMENSION(NINF:NSUP) :: USNEW
      REAL, ALLOCATABLE :: THRESHOLD(:)

      LOGICAL :: LLEPSMIN, LLAK, LLPEAKF

! ----------------------------------------------------------------------


!*    1. THE SWELL DISTRIBUTION IS COMPUTED.
!        -----------------------------------

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            FL1(IJ,K,M) = 0.
          ENDDO
        ENDDO
      ENDDO

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
        IF(.NOT.ALLOCATED(ESWELL)) ALLOCATE(ESWELL(IJS:IJL))
        IF(.NOT.ALLOCATED(FSWELL)) ALLOCATE(FSWELL(IJS:IJL))
        LLAK=.FALSE.
        CALL FEMEAN(FL1, IJS, IJL, ESWELL, FSWELL, FSWELL, LLAK)
      ENDIF

!     femean must be called before semean since LLEPSMIN is false
      IF(IPFGTBL(12).NE.0) THEN
        IF(.NOT.ALLOCATED(ESWELL)) ALLOCATE(ESWELL(IJS:IJL))
        LLEPSMIN=.FALSE.
        CALL SEMEAN(FL1, IJS, IJL, ESWELL, LLEPSMIN)
        DO IJ=IJS,IJL
          IF (USNEW(IJ).EQ.0.) ESWELL(IJ) = EMEAN(IJ)
        ENDDO
      ENDIF
      
      IF(IPFGTBL(14).NE.0) THEN
        IF(.NOT.ALLOCATED(THSWELL)) ALLOCATE(THSWELL(IJS:IJL))
        CALL STHQ   (FL1, IJS, IJL, THSWELL)
      ENDIF

      IF(IPFGTBL(29).NE.0) THEN
        IF(.NOT.ALLOCATED(P1SWELL)) ALLOCATE(P1SWELL(IJS:IJL))
        IF(.NOT.ALLOCATED(ESWELL)) THEN
          ALLOCATE(ESWELL(IJS:IJL))
          LLEPSMIN=.FALSE.
          CALL SEMEAN(FL1, IJS, IJL, ESWELL, LLEPSMIN)
        ENDIF
        CALL MWP1 (FL1, IJS, IJL, ESWELL, P1SWELL)
      ENDIF

      IF(IPFGTBL(31).NE.0) THEN
        IF(.NOT.ALLOCATED(P2SWELL)) ALLOCATE(P2SWELL(IJS:IJL))
        IF(.NOT.ALLOCATED(ESWELL)) THEN
          ALLOCATE(ESWELL(IJS:IJL))
          LLEPSMIN=.FALSE.
          CALL SEMEAN(FL1, IJS, IJL, ESWELL, LLEPSMIN)
        ENDIF
        CALL MWP2 (FL1, IJS, IJL, ESWELL, P2SWELL)
      ENDIF

      IF(IPFGTBL(33).NE.0) THEN
        LLPEAKF = .TRUE.
        IF(.NOT.ALLOCATED(SPRDSWELL)) ALLOCATE(SPRDSWELL(IJS:IJL))
        CALL WDIRSPREAD (FL1, IJS, IJL, LLPEAKF, SPRDSWELL)
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
        IF(.NOT.ALLOCATED(ESEA)) ALLOCATE(ESEA(IJS:IJL))
        IF(.NOT.ALLOCATED(FSEA)) ALLOCATE(FSEA(IJS:IJL))
        LLAK=.FALSE.
        CALL FEMEAN(FL1, IJS, IJL, ESEA, FSEA, FSEA, LLAK)
      ENDIF

!     femean must be called before semean since LLEPSMIN is false
      IF(IPFGTBL(11).NE.0) THEN
        IF(.NOT.ALLOCATED(ESEA)) ALLOCATE(ESEA(IJS:IJL))
        LLEPSMIN=.FALSE.
        CALL SEMEAN(FL1, IJS, IJL, ESEA, LLEPSMIN)
      ENDIF
      
      IF(IPFGTBL(13).NE.0) THEN
        IF(.NOT.ALLOCATED(THWISEA)) ALLOCATE(THWISEA(IJS:IJL))
        CALL STHQ   (FL1, IJS, IJL, THWISEA)
      ENDIF

      IF(IPFGTBL(28).NE.0) THEN
        IF(.NOT.ALLOCATED(P1SEA)) ALLOCATE(P1SEA(IJS:IJL))
        IF(.NOT.ALLOCATED(ESEA)) THEN
          ALLOCATE(ESEA(IJS:IJL))
          LLEPSMIN=.FALSE.
          CALL SEMEAN(FL1, IJS, IJL, ESEA, LLEPSMIN)
        ENDIF
        CALL MWP1 (FL1, IJS, IJL, ESEA, P1SEA)
      ENDIF

      IF(IPFGTBL(30).NE.0) THEN
        IF(.NOT.ALLOCATED(P2SEA)) ALLOCATE(P2SEA(IJS:IJL))
        IF(.NOT.ALLOCATED(ESEA)) THEN
          ALLOCATE(ESEA(IJS:IJL))
          LLEPSMIN=.FALSE.
          CALL SEMEAN(FL1, IJS, IJL, ESEA, LLEPSMIN)
        ENDIF
        CALL MWP2 (FL1, IJS, IJL, ESEA, P2SEA)
      ENDIF

      IF(IPFGTBL(32).NE.0) THEN
        LLPEAKF = .TRUE.
        IF(.NOT.ALLOCATED(SPRDSEA)) ALLOCATE(SPRDSEA(IJS:IJL))
        CALL WDIRSPREAD (FL1, IJS, IJL, LLPEAKF, SPRDSEA)
      ENDIF


      RETURN
      END SUBROUTINE SEPWISW
