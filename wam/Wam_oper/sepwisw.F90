      SUBROUTINE SEPWISW (IJS, IJL, MIJ, FL3, XLLWS,                    &
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

!       *CALL* SEPWISW (IJS, IJL, MIJ, FL3, XLLWS,
!    &                  USNEW, U10NEW, THWNEW,
!    &                  ESWELL   ,FSWELL   ,THSWELL  ,
!    &                  P1SWELL  ,P2SWELL  ,SPRDSWELL,
!    &                  ESEA     ,FSEA     ,THWISEA  ,
!    &                  P1SEA    ,P2SEA    ,SPRDSEA  ,
!    &                  EMTRAIN  ,THTRAIN  ,PMTRAIN)
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *MIJ* - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *FL3* - BLOCK OF SPECTRA
!          *XLLWS* - WINDSEA MASK FROM INPUT SOURCE TERM
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
      USE YOWFRED  , ONLY : C        ,TH       ,FRIC     ,OLDWSFC
      USE YOWPCONS , ONLY : EPSMIN
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,CLDOMAIN
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! -----------------------------------------------------------------------

      IMPLICIT NONE
#include "femean.intfb.h"
#include "mwp1.intfb.h"
#include "mwp2.intfb.h"
#include "sep3tr.intfb.h"
#include "sthq.intfb.h"
#include "wdirspread.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS,IJL
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL), INTENT(IN) :: MIJ

      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: USNEW, U10NEW, THWNEW
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: ESWELL ,FSWELL ,THSWELL, P1SWELL, P2SWELL, SPRDSWELL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: ESEA   ,FSEA   ,THWISEA, P1SEA  , P2SEA  , SPRDSEA
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NTRAIN), INTENT(OUT) :: EMTRAIN, THTRAIN, PMTRAIN
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: FL3
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: XLLWS

      INTEGER(KIND=JWIM) :: IJ, K, M

      REAL(KIND=JWRB) :: COEF
      REAL(KIND=JWRB) :: CHECKTA
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NFRE) :: CM
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: R
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE) :: XINVWVAGE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG) :: DIRCOEF
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: SWM, FL1

      LOGICAL :: LLPEAKF

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SEPWISH',0,ZHOOK_HANDLE)

!*    1. THE SWELL DISTRIBUTION IS COMPUTED.
!        -----------------------------------

      COEF = OLDWSFC*FRIC

      DO M=1,NFRE
        CM(M)=1.0_JWRB/C(M)
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
            IF (XLLWS(IJ,K,M).NE.0.0_JWRB) THEN
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
            FL1(IJ,K,M)=MAX(FL3(IJ,K,M)-FL1(IJ,K,M),0.0_JWRB)
          ENDDO
        ENDDO
      ENDDO
      CALL FEMEAN(FL1, IJS, IJL, ESEA, FSEA)

!     CHECK:
      DO IJ=IJS,IJL
         IF(FSWELL(IJ).GT.FSEA(IJ)) THEN
           R(IJ)=1.0_JWRB
         ELSE
           R(IJ)=0.0_JWRB
         ENDIF
      ENDDO
      DO K=1,NANG
        DO IJ=IJS,IJL
          ! add factor to extend windsea area
          DIRCOEF(IJ,K)=R(IJ)*COEF*SIGN(1.0_JWRB,0.4_JWRB+COS(TH(K)-THWNEW(IJ)))
        ENDDO
      ENDDO
      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
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
      DO IJ=IJS,IJL
        DO K=1,NANG
          DO M=NFRE,2,-1
            IF(SWM(IJ,K,M).EQ.1.0_JWRB .AND. SWM(IJ,K,M-1).EQ.1.0_JWRB) THEN
              EXIT
            ELSEIF(SWM(IJ,K,M).EQ.0.0_JWRB .AND. SWM(IJ,K,M-1).EQ.1.0_JWRB) THEN
               IF(FL3(IJ,K,M).GE.FL3(IJ,K,M-1)) SWM(IJ,K,M-1)=0.0_JWRB
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      ENDIF ! CLDOMAIN

!     SWELL SPECTRUM
      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            FL1(IJ,K,M)=MAX(FL3(IJ,K,M),EPSMIN)*SWM(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO

!*    2.1 COMPUTATION OF THE PARTITIONED SWELL OUTPUT PARAMETERS
!         ------------------------------------------------------

      IF(LLPARTITION) THEN
        CALL FEMEAN(FL1, IJS, IJL, ESWELL, FSWELL)
        CALL STHQ(FL1, IJS, IJL, THSWELL)
        CALL SEP3TR (FL3, IJS, IJL, MIJ, U10NEW, THWNEW,                &
     &               ESWELL, FSWELL, THSWELL, FSEA,                     &
     &               FL1, SWM,                                          &
     &               EMTRAIN  ,THTRAIN  ,PMTRAIN)
      ELSE
        EMTRAIN(:,:)=0.0_JWRB
        THTRAIN(:,:)=0.0_JWRB
        PMTRAIN(:,:)=0.0_JWRB
      ENDIF

!*    2.2 COMPUTATION OF TOTAL SWELL OUTPUT PARAMETERS
!         --------------------------------------------

      CALL FEMEAN(FL1, IJS, IJL, ESWELL, FSWELL)
      
      CALL STHQ(FL1, IJS, IJL, THSWELL)

      CALL MWP1(FL1, IJS, IJL, P1SWELL)

      CALL MWP2(FL1, IJS, IJL, P2SWELL)

      LLPEAKF = .TRUE.
      CALL WDIRSPREAD (FL1, IJS, IJL, ESWELL, LLPEAKF, SPRDSWELL)


!*    3. COMPUTATION OF WIND SEA OUTPUT PARAMETERS
!        -----------------------------------------

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            FL1(IJ,K,M)=MAX(FL3(IJ,K,M)-FL1(IJ,K,M)+EPSMIN,0.0_JWRB)
          ENDDO
        ENDDO
      ENDDO

      CALL FEMEAN(FL1, IJS, IJL, ESEA, FSEA)

      CALL STHQ(FL1, IJS, IJL, THWISEA)
!     if there isn't any windsea energy, set the windsea mean direction
!     to the wind direction.
      DO IJ=IJS,IJL
        IF(ESEA(IJ).LE.0.0_JWRB) THEN
          THWISEA(IJ)=THWNEW(IJ)
        ENDIF
      ENDDO

      CALL MWP1(FL1, IJS, IJL, P1SEA)

      CALL MWP2(FL1, IJS, IJL, P2SEA)

      LLPEAKF = .TRUE.
      CALL WDIRSPREAD(FL1, IJS, IJL, ESEA, LLPEAKF, SPRDSEA)

      IF (LHOOK) CALL DR_HOOK('SEPWISH',1,ZHOOK_HANDLE)

      END SUBROUTINE SEPWISW
