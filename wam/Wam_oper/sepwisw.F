      SUBROUTINE SEPWISW (FL3, FL1, IJS, IJL, IG,THWOLD,USOLD,USNEW)

! ----------------------------------------------------------------------

!**** *SEPWISW* - COMPUTES THE SWELL ENERGY, THE MEAN SWELL DIRECTION,
!****             THE MEAN SWELL FREQUENCY AND THE MEAN WINDSEA DIR.

!     P.LIONELLO     FEBRUARY 87

!     L.ZAMBRESKY    NOVEMBER 87   GKSS/ECMWF   OPTIMIZED SUB.
!     J. BIDLOT      FEBRARY 1996       ECMWF   MESSAGE PASSING

!*    PURPOSE.
!     --------

!       TO SEPARATE THE SWELL FROM THE WIND INTERACTING SEE

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

      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO
      USE YOWPCONS , ONLY : EPSMIN
      USE YOWFRED  , ONLY : FR       ,DFIM     ,C        ,TH       ,
     &            COSTH    ,SINTH
      USE YOWMEAN  , ONLY : EMEAN
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWSWEL  , ONLY : ESWELL   ,FSWELL   ,THSWELL  ,ESEA     ,
     &            FSEA     ,THWISEA

! ----------------------------------------------------------------------

!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE)  :: FL1,FL3

      REAL,DIMENSION(NINF:NSUP,NBLO) :: THWOLD,USOLD
      REAL,DIMENSION(NINF:NSUP) :: USNEW

! ----------------------------------------------------------------------
      PARAMETER (FRIC = 28.)
      REAL,DIMENSION(IJS:IJL) :: CIS, SIS, CIWS, SIWS, SF3SW, SF3WIS


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

      DO K=1,NANG
        DO IJ=IJS,IJL
          SIS(IJ) = 1.2*USOLD(IJ,IG)*COS(TH(K)-THWOLD(IJ,IG))
        ENDDO
        DO M=1,NFRE
          CM = FRIC/C(M)
          DO IJ=IJS,IJL
            CHECKTA = CM*SIS(IJ)
            IF (CHECKTA.LT.1.) THEN
              FL1(IJ,K,M)=1.
            ENDIF
          ENDDO
        ENDDO
      ENDDO

!*    2. COMPUTATION OF SWELL TOTAL ENERGY,MEAN FREQUENCY,MEAN DIRECTION
!*       AND OF WIND-SEA MEAN DIRECTION.
!        ---------------------------------------------------------------

      DO IJ=IJS,IJL
        ESWELL(IJ) = 0.
        FSWELL(IJ) = 0.
        SIS(IJ)    = 0.
        CIS(IJ)    = 0.
        ESEA(IJ)   = 0.
        FSEA(IJ)   = 0.
        CIWS(IJ)   = 0.
        SIWS(IJ)   = 0.
      ENDDO

      DO K=1,NANG
        CIK = COSTH(K)
        SIK = SINTH(K)
        DO IJ=IJS,IJL
          SF3SW(IJ)  = 0.
          SF3WIS(IJ) = 0.
        ENDDO
        DO M=1,NFRE
          DF = DFIM(M)
          FD = DF/FR(M)
          DO IJ=IJS,IJL
            F3SW  = FL3(IJ,K,M)*FL1(IJ,K,M)
            F3WIS = FL3(IJ,K,M)-F3SW
            SF3SW(IJ) = SF3SW(IJ)+F3SW*DF
            SF3WIS(IJ) = SF3WIS(IJ)+F3WIS*DF
            ESWELL(IJ) = ESWELL(IJ)+DF*F3SW
            FSWELL(IJ) = FSWELL(IJ)+FD*F3SW
            ESEA(IJ)   = ESEA(IJ)+DF*F3WIS
            FSEA(IJ)   = FSEA(IJ)+FD*F3WIS
          ENDDO
        ENDDO
        DO IJ=IJS,IJL
          SIS(IJ)  = SIS(IJ)+SF3SW(IJ)*SIK
          CIS(IJ)  = CIS(IJ)+SF3SW(IJ)*CIK
          SIWS(IJ) = SIWS(IJ)+SF3WIS(IJ)*SIK
          CIWS(IJ) = CIWS(IJ)+SF3WIS(IJ)*CIK
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        IF (CIWS(IJ).EQ.0.) CIWS(IJ) = EPSMIN
        IF (CIS(IJ).EQ.0.) CIS(IJ) = EPSMIN
        IF (FSWELL(IJ).EQ.0.) FSWELL(IJ) = EPSMIN
        IF (FSEA(IJ).EQ.0.) FSEA(IJ) = EPSMIN
        IF (USNEW(IJ).EQ.0.) ESWELL(IJ) = EMEAN(IJ)
      ENDDO

      DO IJ=IJS,IJL
        THWISEA(IJ) = ATAN2(SIWS(IJ),CIWS(IJ))
        THSWELL(IJ) = ATAN2(SIS(IJ),CIS(IJ))
        FSWELL(IJ)  = ESWELL(IJ)/(FSWELL(IJ))
        FSEA(IJ)    = ESEA(IJ)/(FSEA(IJ))
      ENDDO

      RETURN
      END SUBROUTINE SEPWISW
