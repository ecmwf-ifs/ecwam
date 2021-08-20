!-------------------------------------------------------------------

      SUBROUTINE MICEP (IPARAM, CICVR, CITH, IJS, IJL)

!-------------------------------------------------------------------

!**** *MICEP* - CLEAN UP SEA FRACTION FOR ALL SEA POINTS.
!               DETERNINE THE SEA ICE THICKNESS IF NOT SUPPLIED
!               AND CHECk THE CONSISTENCY BETWEEN SEA ICE COVER AND
!               THICKNESS. IMPOSE MINIMUM THICKNESS.

!     R. PORTZ     MPI HAMBURG   OCTOBER 1992
!     J. BIDLOT    ECMWF         JUNE 1996    MESSAGE PASSING
!     J. BIDLOT    ECMWF         JANUARY 1998 : CORRECT THRESHOLD FOR
!                                               ICE TO BE LT 271.50
!     J. BIDLOT    ECMWF         FEBRUARY 2000 : USE OF SEA ICE FRACTION
!     J. BIDLOT    ECMWF         AUG 2006 :  


!     PURPOSE.
!     --------
!            SELECTS POINTS (TEMP < 271.5 KELVIN) IN THE TEMPERATURE
!            GRID AND DEFINE SEA COVER AS 1., OTHERWISE 0. 
!            BUT IF SEA ICE FRACTION IS GIVEN DEFINE MODEL SEA ICE
!            COVER MAKING SURE IT IS ALWAYS >=0. AND <=1. 

!**   INTERFACE
!     ---------
!             *CALL MICEP* *(IPARAM, CICVR, CITH, IJS, IJL)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *IPARAM*    INTEGER   GRIB PARAMETER OF FIELDG*CICOVER
!      *CICVR*     REAL      SEA ICE COVER.
!      *CITH*      REAL      SEA ICE THICKNESS.
!      *IJS*       INDEX OF FIRST GRIDPOINT
!      *IJL*       INDEX OF LAST GRIDPOINT


!     METHOD.
!     -------
!             NONE

!    EXTERNAL.
!    ---------
!             *NONE* 

!---------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWICE   , ONLY : CITHRSH  ,LICERUN ,LMASKICE   ,LICETH     , &
     &               HICMIN
      USE YOWMAP   , ONLY : IFROMIJ  ,JFROMIJ
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : NGX      ,NGY     ,CLDOMAIN   ,SWAMPCITH
      USE YOWPCONS , ONLY : ZMISS
      USE YOWTEST  , ONLY : IU06
      USE YOWWIND  , ONLY : FIELDG 
      USE YOWCOUP  , ONLY : LWCOU    ,LWNEMOCOUCIC, LWNEMOCOUCIT
      USE YOWNEMOFLDS, ONLY : NEMOCICOVER, NEMOCITHICK, LNEMOICEREST

      USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IPARAM, IJS, IJL

      REAL(KIND=JWRB), DIMENSION (IJS:IJL), INTENT(INOUT) :: CICVR, CITH


      INTEGER(KIND=JWIM) :: IJ, IX, IY
      INTEGER(KIND=JWIM) :: NICE

!     CONSTANTS FOR PARAMETRISATION OF SEA ICE THICKNESS:
      REAL(KIND=JWRB), PARAMETER :: C1=0.2_JWRB
      REAL(KIND=JWRB), PARAMETER :: C2=0.4_JWRB

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

!     1. INITIALIZATION
!     ------------------

      IF (LHOOK) CALL DR_HOOK('MICEP',0,ZHOOK_HANDLE)

!*    2. WE DEFINE ICE FOR THOSE POINTS WHERE SST < 271.5 KELVIN). 
!        OR CLEAN SEA ICE FRACTION TO DEFINE MODEL SEA ICE COVER.
!        -------------------------------------------------------

      IF (LWNEMOCOUCIC) THEN
        IF (LWCOU) THEN
          DO IJ=IJS,IJL
            IX = IFROMIJ(IJ)
            IY = JFROMIJ(IJ)
            IF (FIELDG(IX,IY)%LKFR <= 0.0_JWRB ) THEN
!            if lake cover = 0, we assume open ocean point, then get sea ice directly from NEMO 
              CICVR(IJ) = NEMOCICOVER(IJ)
            ELSE
!            get ice information from atmopsheric model
              IF (FIELDG(IX,IY)%CICOVER == ZMISS .OR.                   &
     &            FIELDG(IX,IY)%CICOVER < 0.01_JWRB .OR.                &
     &            FIELDG(IX,IY)%CICOVER > 1.01_JWRB ) THEN 
                CICVR(IJ) = 0.0_JWRB
              ELSEIF (FIELDG(IX,IY)%CICOVER > 0.95_JWRB) THEN 
                CICVR(IJ) = 1.0_JWRB
              ELSE
                CICVR(IJ) = FIELDG(IX,IY)%CICOVER 
              ENDIF
            ENDIF
          ENDDO
        ELSE
          DO IJ=IJS,IJL
            IX = IFROMIJ(IJ)
            IY = JFROMIJ(IJ)
            CICVR(IJ) = NEMOCICOVER(IJ)
          ENDDO
        ENDIF
      ELSEIF (IPARAM == 31) THEN
        DO IJ=IJS,IJL
          IX = IFROMIJ(IJ)
          IY = JFROMIJ(IJ)
          IF (FIELDG(IX,IY)%CICOVER == ZMISS .OR.                       &
     &        FIELDG(IX,IY)%CICOVER < 0.01_JWRB .OR.                    &
     &        FIELDG(IX,IY)%CICOVER > 1.01_JWRB ) THEN 
            CICVR(IJ) = 0.0_JWRB
          ELSEIF (FIELDG(IX,IY)%CICOVER .GT. 0.95_JWRB) THEN 
            CICVR(IJ) = 1.0_JWRB
          ELSE
            CICVR(IJ) = FIELDG(IX,IY)%CICOVER 
          ENDIF
        ENDDO
      ELSEIF (IPARAM == 139) THEN
        DO IJ=IJS,IJL
          IX = IFROMIJ(IJ)
          IY = JFROMIJ(IJ)
          IF (FIELDG(IX,IY)%CICOVER < 271.5_JWRB) THEN
            CICVR(IJ) = 1.0_JWRB
          ELSE
            CICVR(IJ) = 0.0_JWRB
          ENDIF
        ENDDO
      ENDIF 


      IF (.NOT. LICERUN .OR. LMASKICE) THEN
!       SEA ICE THICKNESS IN CASE IT IS NOT SUPPLIED AS INPUT:
        DO IJ=IJS,IJL
          CITH(IJ)=0.0_JWRB
        ENDDO
      ELSEIF (CLDOMAIN == 's') THEN
        DO IJ=IJS,IJL
          IF (CICVR(IJ) > 0.0_JWRB) THEN
            CITH(IJ)=SWAMPCITH
          ELSE
            CITH(IJ)=0.0_JWRB
          ENDIF
        ENDDO
      ELSEIF ((.NOT.LICETH) .AND. (.NOT.LWNEMOCOUCIT)) THEN
!       SEA ICE THICKNESS IS PARAMETERISED:
        DO IJ=IJS,IJL
          IF (CICVR(IJ) > 0.0_JWRB) THEN
            CITH(IJ)=MAX(C1+C2*CICVR(IJ),0.0_JWRB)
          ELSE
            CITH(IJ)=0.0_JWRB
          ENDIF
        ENDDO

      ELSEIF (LWNEMOCOUCIT) THEN
        IF (LWCOU) THEN
          DO IJ=IJS,IJL
            IX = IFROMIJ(IJ)
            IY = JFROMIJ(IJ)
            IF (FIELDG(IX,IY)%LKFR <= 0.0_JWRB ) THEN
!             if lake cover = 0, we assume open ocean point, then get sea ice thickness directly from NEMO 
              IF (LNEMOICEREST) THEN
                CITH(IJ)=NEMOCITHICK(IJ)
              ELSE
                CITH(IJ)=CICVR(IJ)*NEMOCITHICK(IJ)
              ENDIF
            ELSE
!           We should get ice thickness information from atmopsheric model
!           but it is not yet coded. For now, parameterise it from the cover...
              IF (CICVR(IJ) > 0.0_JWRB) THEN
                CITH(IJ)=MAX(C1+C2*CICVR(IJ),0.0_JWRB)
              ELSE
               CITH(IJ)=0.0_JWRB
              ENDIF
            ENDIF
          ENDDO

        ELSE
          IF (LNEMOICEREST) THEN
            DO IJ=IJS,IJL
               CITH(IJ)=NEMOCITHICK(IJ)
            ENDDO
          ELSE
            DO IJ=IJS,IJL
               CITH(IJ)=CICVR(IJ)*NEMOCITHICK(IJ)
            ENDDO
          ENDIF
        ENDIF

!       CONSISTENCY CHECK:
!       no ice if thickness < 0.5*HICMIN
        DO IJ=IJS,IJL
          IF (CICVR(IJ) > 0.0_JWRB .AND. CITH(IJ) < 0.5_JWRB*HICMIN) THEN
            CICVR(IJ)=0.0_JWRB
            CITH(IJ)=0.0_JWRB
          ENDIF
        ENDDO

      ELSE

!!!!   define a representative sea ice thickness that account for sea ice coverage
        DO IJ=IJS,IJL
          CITH(IJ)=CICVR(IJ)*CITH(IJ)
        ENDDO

!       CONSISTENCY CHECK:
!       no ice if thickness < 0.5*HICMIN
        DO IJ=IJS,IJL
          IF (CICVR(IJ) > 0.0_JWRB .AND. CITH(IJ) < 0.5_JWRB*HICMIN) THEN
            CICVR(IJ)=0.0_JWRB
            CITH(IJ)=0.0_JWRB
          ENDIF
        ENDDO

      ENDIF

      IF (LHOOK) CALL DR_HOOK('MICEP',1,ZHOOK_HANDLE)

      END SUBROUTINE MICEP
