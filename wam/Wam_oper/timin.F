      SUBROUTINE TIMIN (CDTWIS, CDTWIE,
     &                  MIJS, MIJL,
     &                  U10OLD, THWOLD, USOLD, TAUW, Z0OLD,
     &                  ROAIRO, ZIDLOLD, CICOVER, CITHICK, IREAD,
     &                  LWCUR)

! ----------------------------------------------------------------------
!     MODIFIED  J. BIDLOT  FEBRARY 1996  MESSAGE PASSING
!               S. ABDALLA OCTOBER 2001  INCLUSION OF AIR DENSITY & Zi/L
!               J BIDLOT AUGUST 2006 INCLUSION OF ICE MASK.
!               J BIDLOT AUGUST 2008 REMOVE INFORMATION FROM COUPLED RUNS
!                                    AND CLEAN UP OBSOLETE OPTIONS.



!**** *TIMIN* - STEERING MODULE IF TIME INTERPOLATION IS WANTED.

!     H. GUNTHER    ECMWF    MAY 1990     MODIFIED FOR SUB VERSION.
!     H. GUNTHER    ECMWF    DECEMBER 90  MODIFIED FOR CYCLE_4.

!*    PURPOSE.
!     --------

!       TIME INTERPOLATION: PROCESS WINDFIELDS.

!**   INTERFACE.
!     ----------

!       *CALL* *TIMIN (CDTWIS, CDTWIE,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
!    &                 ROAIRO, ZIDLOLD, CICOVER, CITHICK, IREAD,
!    &                 LWCUR)*
!          *CDTWIS* - DATE OF FIRST WIND FIELD.
!          *CDTWIE* - DATE OF LAST WIND FIELD.
!          *U10OLD* - WIND SPEED.
!          *THWOLD* - WIND DIRECTION (RADIANS).
!          *USOLD*  - FRICTION VELOCITY.
!          *TAUW*   - WAVE STRESS.
!          *Z0OLD*  - ROUGHNESS LENGTH IN M.
!          *ROAIRO* - AIR DENSITY IN KG/M3.
!          *ZIDLOLD*- Zi/L
!                     (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!          *CICOVER*- SEA ICE COVER.
!          *CITHICK*- SEA ICE THICKNESS.
!          *IREAD*  - PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
!          *LWCUR*  - LOGICAL INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS


!     METHOD.
!     -------

!       WINDFIELDS ARE READ IN EVERY IDELWI SECONDS (U,V),
!       INTERPOLATED IN SPACE, BLOCKED AND SAVED ON SCRATCH UNITS.
!       MAGNITUDE AND DIRECTION ARE INTERPOLATED LINEARLY IN TIME AND
!       WRITTEN TO THE OUTPUT FILES.

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.
!       *AIRSEA*    - TOTAL STRESS IN SURFACE LAYER.
!       *GETWND*    - READ A WINDFIELD AND COMPUTE WIND
!                     FOR ALL BLOCKS (US,DS).
!       *INCDATE*   - INCREMENT DATE.
!       *LOCINT*    - SPACE INTERPOLATION.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPCONS , ONLY : PI       ,ZPI
      USE YOWSTAT  , ONLY : IDELPRO  ,IDELWI   ,IDELWO   ,NPROMA_WAM
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWIND  , ONLY : CDA      ,CDTNEXT  ,NSTORE   ,FF_NEXT
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "airsea.intfb.h"
#include "getwnd.intfb.h"
#include "incdate.intfb.h"

      INTEGER, INTENT(IN) :: MIJS, MIJL
      INTEGER, INTENT(IN) :: IREAD
      
      REAL,DIMENSION(MIJS:MIJL), INTENT(INOUT) ::
     &               U10OLD, THWOLD, USOLD, Z0OLD, TAUW,
     &               ROAIRO, ZIDLOLD, CICOVER, CITHICK

      CHARACTER(LEN=14), INTENT(IN) :: CDTWIS, CDTWIE

      LOGICAL, INTENT(IN) :: LWCUR


      INTEGER :: JKGLO,KIJS,KIJL,NPROMA
      INTEGER :: ILEV, IJ
      INTEGER :: ICODE_WND
      INTEGER :: NTS, N

      REAL :: ZHOOK_HANDLE
      REAL :: DEL, D
      REAL,DIMENSION(MIJS:MIJL) :: U10,US,DS,ADS,ZIDL,CICR,CITH
      REAL,DIMENSION(MIJS:MIJL) :: U102,US2,DS2,ADS2,ZIDL2,CICR2,CITH2
      REAL,DIMENSION(MIJS:MIJL) :: U103,US3, DS3,ADS3,ZIDL3,CICR3,CITH3

      CHARACTER(LEN=14) :: CDTH, CDT1, CDT2, ZERO
      CHARACTER(LEN=80) :: FILENAME

      LOGICAL :: LWNDFILE, LCLOSEWND

! ----------------------------------------------------------------------

!*    1. INITIALIZE TIMECOUNTER.
!        -----------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('TIMIN',0,ZHOOK_HANDLE)
#endif

      ZERO= ' '

      NPROMA=NPROMA_WAM

      LCLOSEWND=.FALSE.
      IF (LWCOU) THEN
        LWNDFILE=.FALSE.
      ELSE
        LWNDFILE=.TRUE.
      ENDIF
      ILEV=1

      IF (CDA.EQ.ZERO) THEN
        CDT1 = CDTWIS
        CALL GETWND (MIJS, MIJL,
     &               U10(MIJS), US(MIJS),
     &               DS(MIJS),
     &               ADS(MIJS), ZIDL(MIJS),
     &               CICR(MIJS), CITH(MIJS),
     &               CDT1, LWNDFILE, LCLOSEWND, IREAD,
     &               LWCUR, ICODE_WND)

        CALL GSTATS(1494,0)
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ)
        DO JKGLO=MIJS,MIJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,MIJL)

          DO IJ=KIJS,KIJL
            U10OLD(IJ) = U10(IJ)
            THWOLD(IJ) = DS(IJ)
            ROAIRO(IJ) = ADS(IJ)
            ZIDLOLD(IJ)= ZIDL(IJ)
            CICOVER(IJ)= CICR(IJ)
            CITHICK(IJ)= CITH(IJ)
          ENDDO

          CALL AIRSEA (U10OLD(KIJS), TAUW(KIJS), USOLD(KIJS),
     &                 Z0OLD(KIJS), KIJS, KIJL, ILEV, ICODE_WND)

        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1494,1)

        IF (ITEST.GE.3) THEN
          WRITE(IU06,'(''       SUB. TIMIN: FIRST WIND FIELD '',
     &     ''SAVED IN COMMON WIND'')')
          CALL FLUSH(IU06)
        ENDIF

      ELSE
        CDT1 = CDA
        CALL GSTATS(1494,0)
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ)
        DO JKGLO=MIJS,MIJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,MIJL)
          DO IJ=KIJS,KIJL
            U10(IJ) = U10OLD(IJ)
            US(IJ) = USOLD(IJ)
            DS(IJ) = THWOLD(IJ)
            ADS(IJ) = ROAIRO(IJ)
            ZIDL(IJ)= ZIDLOLD(IJ)
            CICR(IJ)= CICOVER(IJ)
            CITH(IJ)= CITHICK(IJ)
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1494,1)
        CDTH = CDT1
        CALL INCDATE (CDTH,IDELWO)
        IF (CDTWIS.NE.CDTH) THEN
          WRITE(IU06,*) ' *******************************************'
          WRITE(IU06,*) ' *                                         *'
          WRITE(IU06,*) ' *        FATAL ERROR IN --TIMIN--         *'
          WRITE(IU06,*) ' *        =========================        *'
          WRITE(IU06,*) ' * DATES DO NOT MATCH.                     *'
          WRITE(IU06,*) ' * START DATE FOR WIND IS       CDTWIS = ',
     &     CDTWIS
          WRITE(IU06,*) ' * LAST DATE SAVED IN COM WIND IS CDT1 = ',
     &     CDT1
          WRITE(IU06,*) ' * PROCESSING WILL BE ABORTED              *'
          WRITE(IU06,*) ' *                                         *'
          WRITE(IU06,*) ' *******************************************'
          CALL ABORT1
        ENDIF
      ENDIF
      CDT2 = CDT1
      CALL INCDATE(CDT2,IDELWI)

      NTS = IDELWI/IDELWO
      DEL = REAL(IDELWO)/REAL(IDELWI)

!*    2. LOOP OVER INPUT WINDFIELDS.
!        ---------------------------

 2000 CONTINUE

!*    2.1 READ ONE WINDFIELD AND TRANSFORM ALL BLOCKS.
!         -------------------------------------------

      CDT2 = CDT1
      CALL INCDATE(CDT2,IDELWI)
      CALL GETWND (MIJS, MIJL,
     &             U102(MIJS), US2(MIJS),
     &             DS2(MIJS),
     &             ADS2(MIJS), ZIDL2(MIJS),
     &             CICR2(MIJS), CITH2(MIJS),
     &             CDT2, LWNDFILE, LCLOSEWND, IREAD,
     &             LWCUR, ICODE_WND)

!*    2.2 SAVE BLOCKED WIND FIELD ON SCRATCH UNITS.
!         -----------------------------------------

      NSTORE=NTS

      IF(.NOT.ALLOCATED(CDTNEXT)) ALLOCATE(CDTNEXT(NSTORE))
        IF(.NOT.ALLOCATED(FF_NEXT))
     &          ALLOCATE(FF_NEXT(MIJS:MIJL,NSTORE))

        CDTH = CDT1
        DO N=1,NTS
          CALL INCDATE(CDTH,IDELWO)
          CDTNEXT(N)=CDTH
          CALL GSTATS(1494,0)
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ,D)
          DO JKGLO=MIJS,MIJL,NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,MIJL)
            DO IJ=KIJS,KIJL
              U103(IJ) = U10(IJ)+REAL(N)*DEL*(U102(IJ)-U10(IJ))
              US3(IJ) = US(IJ)+REAL(N)*DEL*(US2(IJ)-US(IJ))
              D = DS2(IJ)-DS(IJ)
              IF (ABS(D).GT.PI) D = D-ZPI*SIGN(1.,D)
              D = DS(IJ)+REAL(N)*DEL*D
              DS3(IJ) = MOD(D+ZPI,ZPI)
              ADS3(IJ) = ADS(IJ)+
     &                        REAL(N)*DEL*(ADS2(IJ)-ADS(IJ))
              ZIDL3(IJ) = ZIDL(IJ)+
     &                         REAL(N)*DEL*(ZIDL2(IJ)-ZIDL(IJ))
              CICR3(IJ) = CICR(IJ)+
     &                         REAL(N)*DEL*(CICR2(IJ)-CICR(IJ))
              CITH3(IJ) = CITH(IJ)+
     &                         REAL(N)*DEL*(CITH2(IJ)-CITH(IJ))

              FF_NEXT(IJ,N)%WSWAVE=U103(IJ)
              FF_NEXT(IJ,N)%USTAR=US3(IJ)
              FF_NEXT(IJ,N)%WDWAVE=DS3(IJ)
              FF_NEXT(IJ,N)%AIRD=ADS3(IJ)
              FF_NEXT(IJ,N)%ZIDL=ZIDL3(IJ)
              FF_NEXT(IJ,N)%CIFR=CICR3(IJ)
              FF_NEXT(IJ,N)%CITH=CITH3(IJ)
            ENDDO
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1494,1)
        ENDDO

        DO IJ=MIJS,MIJL
          U10(IJ) = U102(IJ)
          US(IJ) = US2(IJ)
          DS(IJ) = DS2(IJ)
          ADS(IJ) = ADS2(IJ)
          ZIDL(IJ) = ZIDL2(IJ)
          CICR(IJ) = CICR2(IJ)
          CITH(IJ) = CITH2(IJ)
        ENDDO

      CDT1 = CDT2
      IF (ITEST.GE.3) THEN
        WRITE(IU06,*) '       SUB. TIMIN: LAST WIND FIELD AT ',
     &   'CDTH = ', CDTH,' WRITTEN TO SCRATCH UNITS'
          CALL FLUSH (IU06)
      ENDIF

!*    2.3 UPDATE WIND FIELD REQUEST TIME AND READ NEXT IF REQUESTED.
!         ----------------------------------------------------------

      CALL INCDATE (CDTH,IDELWI)
!!! we don't allow to loop anymore
!!!!      IF (CDTH.LE.CDTWIE) GOTO 2000


      IF (ITEST.GE.3) THEN
        WRITE(IU06,*) '       SUB. TIMIN: LAST WIND FIELD AT ',
     &   'CDTH = ', CDTH,' WRITTEN TO OUTPUT UNIT'
          CALL FLUSH (IU06)
      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('TIMIN',1,ZHOOK_HANDLE)
#endif
      RETURN
      END SUBROUTINE TIMIN
