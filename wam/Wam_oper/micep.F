!-------------------------------------------------------------------

      SUBROUTINE MICEP (IPARAM, CICVR, CITH, IJS, IJL, IG)

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
!             *CALL MICEP* *(IPARAM, CICVR, CITH, IJS, IJL, IG)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *IPARAM*    INTEGER   GRIB PARAMETER OF FIELDG*CIFR
!      *CICVR*     REAL      SEA ICE COVER.
!      *CITH*      REAL      SEA ICE THICKNESS.
!      *IJS*       INDEX OF FIRST GRIDPOINT
!      *IJL*       INDEX OF LAST GRIDPOINT
!      *IG*        BLOCK NUMBER


!     METHOD.
!     -------
!             NONE

!    EXTERNAL.
!    ---------
!             *NONE* 

!---------------------------------------------------------------------

      USE YOWICE   , ONLY : CITHRSH  ,LICERUN ,LMASKICE   ,LICETH     ,
     &               HICMIN
      USE YOWMAP   , ONLY : IXLG     ,KXLT
      USE YOWMPP   , ONLY : IRANK    ,NPROC   ,NINF       ,NSUP
      USE YOWPARAM , ONLY : NGX      ,NGY     ,NBLO       ,CLDOMAIN   ,
     &            SWAMPCITH
      USE YOWPCONS , ONLY : ZMISS
      USE YOWTEST  , ONLY : IU06     ,ITEST 
      USE YOWUNPOOL ,ONLY : LLUNSTR
      USE YOWWIND  , ONLY : FIELDG 
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWCOUP  , ONLY : LWNEMOCOUCIC, LWNEMOCOUCIT
      USE YOWNEMOFLDS, ONLY : NEMOCICOVER, NEMOCITHICK, LNEMOICEREST

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IG
      INTEGER, INTENT(IN) :: IPARAM, IJS, IJL

      REAL, DIMENSION (NINF:NSUP,NBLO), INTENT(INOUT) :: CICVR, CITH 


      INTEGER :: IJ, IX, IY
      INTEGER :: NICE

!     CONSTANTS FOR PARAMETRISATION OF SEA ICE THICKNESS:
      REAL, PARAMETER :: C1=0.2
      REAL, PARAMETER :: C2=0.4

      REAL :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

!     1. INITIALIZATION
!     ------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('MICEP',0,ZHOOK_HANDLE)
#endif

!-----------------------------------------------------------------

!*    2. WE DEFINE ICE FOR THOSE POINTS WHERE SST < 271.5 KELVIN). 
!        OR CLEAN SEA ICE FRACTION TO DEFINE MODEL SEA ICE COVER.
!        -------------------------------------------------------


      IF(LWNEMOCOUCIC) THEN
        DO IJ=IJS,IJL
          CICVR(IJ,IG) = NEMOCICOVER(IJ)
        ENDDO
      ELSEIF (IPARAM.EQ.31) THEN
        DO IJ=IJS,IJL
          IF (LLUNSTR) THEN
            IX = IJ
            IY = 1
          ELSE
            IX = IXLG(IJ,IG)
            IY = NGY-KXLT(IJ,IG) + 1
          ENDIF
          IF (FIELDG(IX,IY)%CIFR .EQ. ZMISS .OR. 
     &        FIELDG(IX,IY)%CIFR .LT. 0.01 .OR.
     &        FIELDG(IX,IY)%CIFR .GT. 1.01 ) THEN 
            CICVR(IJ,IG) = 0.0
          ELSEIF (FIELDG(IX,IY)%CIFR .GT. 0.95) THEN 
            CICVR(IJ,IG) = 1.0
          ELSE
            CICVR(IJ,IG) = FIELDG(IX,IY)%CIFR 
          ENDIF
        ENDDO
      ELSE IF(IPARAM.EQ.139) THEN
        DO IJ=IJS,IJL
          IF (LLUNSTR) THEN
            IX = IJ
            IY = 1
          ELSE
            IX = IXLG(IJ,IG)
            IY = NGY-KXLT(IJ,IG) + 1
          ENDIF
          IF (FIELDG(IX,IY)%CIFR .LT. 271.5) THEN
            CICVR(IJ,IG) = 1.0
          ELSE
            CICVR(IJ,IG) = 0.0
          ENDIF
        ENDDO
      ENDIF 


      IF( .NOT. LICERUN .OR. LMASKICE ) THEN
!      SEA ICE THICKNESS IN CASE IT IS NOT SUPPLIED AS INPUT:
        DO IJ=IJS,IJL
          CITH(IJ,IG)=0.0
        ENDDO
      ELSE IF(CLDOMAIN == 's' ) THEN
        DO IJ=IJS,IJL
          IF(CICVR(IJ,IG).GT.0.0) THEN
            CITH(IJ,IG)=SWAMPCITH
          ELSE
            CITH(IJ,IG)=0.0
          ENDIF
        ENDDO
      ELSE IF((.NOT.LICETH).AND.(.NOT.LWNEMOCOUCIT)) THEN
!       SEA ICE THICKNESS IS PARAMETERISED:
        DO IJ=IJS,IJL
          IF(CICVR(IJ,IG).GT.0.0) THEN
            CITH(IJ,IG)=MAX(C1+C2*CICVR(IJ,IG),0.)
          ELSE
            CITH(IJ,IG)=0.0
          ENDIF
        ENDDO
      ELSE IF(LWNEMOCOUCIT) THEN

         IF(LNEMOICEREST) THEN
           DO IJ=IJS,IJL
              CITH(IJ,IG)=NEMOCITHICK(IJ)
           ENDDO
         ELSE
           DO IJ=IJS,IJL
              CITH(IJ,IG)=CICVR(IJ,IG)*NEMOCITHICK(IJ)
           ENDDO
         ENDIF

!       CONSISTENCY CHECK:
!       no ice if thickness < 0.5*HICMIN
        DO IJ=IJS,IJL
          IF(CICVR(IJ,IG).GT.0.0 .AND. CITH(IJ,IG).LT.0.5*HICMIN ) THEN
            CICVR(IJ,IG)=0.0
            CITH(IJ,IG)=0.0
          ENDIF
        ENDDO

      ELSE

!!!!   define a representative sea ice thickness that account for sea ice coverage
        DO IJ=IJS,IJL
          CITH(IJ,IG)=CICVR(IJ,IG)*CITH(IJ,IG)
        ENDDO

!       CONSISTENCY CHECK:
!       no ice if thickness < 0.5*HICMIN
        DO IJ=IJS,IJL
          IF(CICVR(IJ,IG).GT.0.0 .AND. CITH(IJ,IG).LT.0.5*HICMIN ) THEN
            CICVR(IJ,IG)=0.0
            CITH(IJ,IG)=0.0
          ENDIF
        ENDDO

      ENDIF


      IF(ITEST.GE.1) THEN
        NICE = 0
        DO IJ=IJS,IJL
          IF(CICVR(IJ,IG) .GT. CITHRSH) NICE = NICE + 1 
        ENDDO
        IF(IPARAM.EQ.139) THEN
          WRITE(IU06,*)' SUB. MICEP: POINTS WITH TEMP < 271.5 K :',
     &     NICE,' ON PE ',IRANK, ' OUT OF A TOTAL OF ',NPROC, 'PEs'
        ELSEIF (IPARAM.EQ.31) THEN
          WRITE(IU06,*)' SUB. MICEP: ',
     &           'POINTS WITH SEA ICE FRACTION > ', CITHRSH, ': ',
     &     NICE,' ON PE ',IRANK, ' OUT OF A TOTAL OF ',NPROC, 'PEs'
        ELSE
          WRITE(IU06,*)' ********************************************'
          WRITE(IU06,*)' SUB. MICEP: PARAMETER FOR ICE MASK NOT KNOWN'
          WRITE(IU06,*)' SUB. MICEP: NO ICE MASK WAS SET UP'
          WRITE(IU06,*)' ********************************************'
        ENDIF
      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('MICEP',1,ZHOOK_HANDLE)
#endif
      RETURN
      END SUBROUTINE MICEP
