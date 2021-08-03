      SUBROUTINE SEPWISW (IJS, IJL, KIJS, KIJL, MIJ, GFL, GXLLWS, CINV, &
     &                    USNEW, U10NEW, THWNEW,                        &
     &                    ESWELL   ,FSWELL   ,THSWELL  ,                &
     &                    P1SWELL  ,P2SWELL  ,SPRDSWELL,                &
     &                    ESEA     ,FSEA     ,THWISEA  ,                &
     &                    P1SEA    ,P2SEA    ,SPRDSEA  ,                &
     &                    EMTRAIN  ,THTRAIN  ,PMTRAIN)

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

!       *CALL* SEPWISW (IJS, IJL, KIJS, KIJL, MIJ, GFL, GXLLWS,
!    &                  USNEW, U10NEW, THWNEW,
!    &                  ESWELL   ,FSWELL   ,THSWELL  ,
!    &                  P1SWELL  ,P2SWELL  ,SPRDSWELL,
!    &                  ESEA     ,FSEA     ,THWISEA  ,
!    &                  P1SEA    ,P2SEA    ,SPRDSEA  ,
!    &                  EMTRAIN  ,THTRAIN  ,PMTRAIN)
!          *IJS:IJL* - 1st DIMENSION of GFL and GXLLWS
!          *KIJS* - INDEX OF FIRST GRIDPOINT
!          *KIJL* - INDEX OF LAST GRIDPOINT
!          *MIJ* - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *GFL* - BLOCK OF SPECTRA
!          *GXLLWS* - WINDSEA MASK FROM INPUT SOURCE TERM
!          *USNEW* - NEW FRICTION VELOCITY IN M/S.
!          *U10NEW* - LATESt WIND SPEED.
!          *THWNEW* - LATEST WIND DIRECTION.
!          *        - SWELL and WINDSEA PARAMETERS.

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

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : NTRAIN   ,LLPARTITION
      USE YOWFRED  , ONLY : FR       ,TH       ,FRIC     ,OLDWSFC, ZPIFR
      USE YOWPCONS , ONLY : G        ,EPSMIN
      USE YOWSTAT  , ONLY : ISHALLO
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,CLDOMAIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! -----------------------------------------------------------------------

      IMPLICIT NONE
#include "femean.intfb.h"
#include "mwp1.intfb.h"
#include "mwp2.intfb.h"
#include "sep3tr.intfb.h"
#include "sthq.intfb.h"
#include "wdirspread.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, KIJS, KIJL
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN) :: MIJ

      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: GFL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: GXLLWS
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE), INTENT(IN) :: CINV 

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: USNEW, U10NEW, THWNEW
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: ESWELL ,FSWELL ,THSWELL, P1SWELL, P2SWELL, SPRDSWELL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: ESEA   ,FSEA   ,THWISEA, P1SEA  , P2SEA  , SPRDSEA
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NTRAIN), INTENT(OUT) :: EMTRAIN, THTRAIN, PMTRAIN


      INTEGER(KIND=JWIM) :: IJ, K, M

      REAL(KIND=JWRB) :: COEF
      REAL(KIND=JWRB) :: CHECKTA
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: R
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: XINVWVAGE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG) :: DIRCOEF
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: SWM, F1

      LOGICAL :: LLPEAKF

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SEPWISH',0,ZHOOK_HANDLE)

!*    1. THE SWELL DISTRIBUTION IS COMPUTED.
!        -----------------------------------

      COEF = OLDWSFC*FRIC

      DO M=1,NFRE
        DO IJ=KIJS,KIJL
          XINVWVAGE(IJ,M)=USNEW(IJ)*CINV(IJ,M)
        ENDDO
      ENDDO

      DO K=1,NANG
        DO IJ=KIJS,KIJL
          DIRCOEF(IJ,K)=COEF*COS(TH(K)-THWNEW(IJ))
        ENDDO
      ENDDO

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            IF (GXLLWS(IJ,K,M).NE.0.0_JWRB) THEN
              ! this is windsea 
              SWM(IJ,K,M)=0.0_JWRB
            ELSE
              CHECKTA=XINVWVAGE(IJ,M)*DIRCOEF(IJ,K)
              IF (CHECKTA.GE.1.0_JWRB) THEN
                ! this is extra windsea 
                SWM(IJ,K,M)=0.0_JWRB
              ELSE
                ! this is swell
                SWM(IJ,K,M)=1.0_JWRB
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      IF(.NOT.(CLDOMAIN.EQ.'s')) THEN
!     CHECK THAT TOTAL SWELL MEAN FREQUENCY IS LOWER THAN WINDSEA ONE
!     OTHERWISE RESET WIND SECTOR TO WINDSEA

!     SWELL:
      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            F1(IJ,K,M)=GFL(IJ,K,M)*SWM(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO
      CALL FEMEAN(F1, KIJS, KIJL, KIJS, KIJL, ESWELL, FSWELL)

!     WINDSEA:
      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            F1(IJ,K,M)=MAX(GFL(IJ,K,M)-F1(IJ,K,M),0.0_JWRB)
          ENDDO
        ENDDO
      ENDDO
      CALL FEMEAN(F1, KIJS, KIJL, KIJS, KIJL, ESEA, FSEA)

!     CHECK:
      DO IJ=KIJS,KIJL
         IF(FSWELL(IJ).GT.0.96_JWRB*FSEA(IJ)) THEN
           R(IJ)=1.0_JWRB
         ELSE
           R(IJ)=0.0_JWRB
         ENDIF
      ENDDO
      DO K=1,NANG
        DO IJ=KIJS,KIJL
          ! add factor to extend windsea area
          DIRCOEF(IJ,K)=R(IJ)*COEF*SIGN(1.0_JWRB,0.4_JWRB+COS(TH(K)-THWNEW(IJ)))
        ENDDO
      ENDDO
      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            CHECKTA = XINVWVAGE(IJ,M)*DIRCOEF(IJ,K)
            IF (CHECKTA.GE.1.0_JWRB) THEN
              ! this is additional windsea 
              SWM(IJ,K,M)=0.0_JWRB
            ENDIF
          ENDDO
        ENDDO
      ENDDO

!     ADJUST THE LOW FREQUENCY BOUNDARY OF THE WINDSEAS AREA 
!     TO TOPOLOGICALLY CONNECT IT TO THE WINDSEA
      DO IJ=KIJS,KIJL
        DO K=1,NANG
          DO M=NFRE,2,-1
            IF(SWM(IJ,K,M).EQ.1.0_JWRB .AND. SWM(IJ,K,M-1).EQ.1.0_JWRB) THEN
              EXIT
            ELSEIF(SWM(IJ,K,M).EQ.0.0_JWRB .AND. SWM(IJ,K,M-1).EQ.1.0_JWRB) THEN
               IF(GFL(IJ,K,M).GE.GFL(IJ,K,M-1)) SWM(IJ,K,M-1)=0.0_JWRB
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      ENDIF ! CLDOMAIN

!     SWELL SPECTRUM
      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            F1(IJ,K,M)=MAX(GFL(IJ,K,M),EPSMIN)*SWM(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO

!*    2.1 COMPUTATION OF THE PARTITIONED SWELL OUTPUT PARAMETERS
!         ------------------------------------------------------

      IF(LLPARTITION) THEN
        CALL FEMEAN(F1, KIJS, KIJL, KIJS, KIJL, ESWELL, FSWELL)
        CALL STHQ(F1, KIJS, KIJL, KIJS, KIJL, THSWELL)
        CALL SEP3TR (GFL, IJS, IJL, KIJS, KIJL, MIJ, U10NEW, THWNEW,   &
     &               ESWELL, FSWELL, THSWELL, FSEA,                    &
     &               F1, SWM,                                          &
     &               EMTRAIN  ,THTRAIN  ,PMTRAIN)
      ELSE
        EMTRAIN(:,:) = 0.0_JWRB
        THTRAIN(:,:) = 0.0_JWRB
        PMTRAIN(:,:) = 0.0_JWRB
      ENDIF

!*    2.2 COMPUTATION OF TOTAL SWELL OUTPUT PARAMETERS
!         --------------------------------------------

      CALL FEMEAN(F1, KIJS, KIJL, KIJS, KIJL, ESWELL, FSWELL)
      
      CALL STHQ(F1, KIJS, KIJL, KIJS, KIJL, THSWELL)

      CALL MWP1(F1, KIJS, KIJL, KIJS, KIJL, P1SWELL)

      CALL MWP2(F1, KIJS, KIJL, KIJS, KIJL, P2SWELL)

      LLPEAKF = .TRUE.
      CALL WDIRSPREAD (F1, KIJS, KIJL, KIJS, KIJL, ESWELL, LLPEAKF, SPRDSWELL)


!*    3. COMPUTATION OF WIND SEA OUTPUT PARAMETERS
!        -----------------------------------------

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            F1(IJ,K,M)=MAX(GFL(IJ,K,M)-F1(IJ,K,M)+EPSMIN,0.0_JWRB)
          ENDDO
        ENDDO
      ENDDO

      CALL FEMEAN(F1, KIJS, KIJL, KIJS, KIJL, ESEA, FSEA)

      CALL STHQ(F1, KIJS, KIJL, KIJS, KIJL, THWISEA)
!     if there isn't any windsea energy, set the windsea mean direction
!     to the wind direction.
      DO IJ=KIJS,KIJL
        IF(ESEA(IJ).LE.0.0_JWRB) THEN
          THWISEA(IJ)=THWNEW(IJ)
        ENDIF
      ENDDO

      CALL MWP1(F1, KIJS, KIJL, KIJS, KIJL, P1SEA)

      CALL MWP2(F1, KIJS, KIJL, KIJS, KIJL, P2SEA)

      LLPEAKF = .TRUE.
      CALL WDIRSPREAD(F1, KIJS, KIJL, KIJS, KIJL, ESEA, LLPEAKF, SPRDSEA)

      IF (LHOOK) CALL DR_HOOK('SEPWISH',1,ZHOOK_HANDLE)

      END SUBROUTINE SEPWISW
