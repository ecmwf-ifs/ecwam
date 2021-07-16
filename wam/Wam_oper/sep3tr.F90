      SUBROUTINE SEP3TR (F3, IJS, IJL, MIJ, U10NEW, THWNEW ,           &
     &                   ESWELL, FSWELL, THSWELL, FSEA,                 &
     &                   F1, SWM,                                      &
     &                   EMTRAIN  ,THTRAIN  ,PMTRAIN)

! ----------------------------------------------------------------------

!**** *SEP3TR* -  COMPUTE ENERGY, DIRECTION AND PERIOD FOR THE 3 WAVE
!****             TRAINS: SWELL1, SWELL2, SWELL3

!     D.PETTENUZZO     MAY 2012
!     JEAN BIDLOT      JUNE 2013  SIMPLIFIED BY REMOVING THE WIND SEA PART
!                                 AND ONLY PARTITIONED THE SWELL SPECTRUM
!                                 AND ONLY LOOKING IN THE PROGNOSTIC RANGE


!*    PURPOSE.
!     --------

!     CREATED TO TEST WAVE SWELL SPECTRA PARTITIONING INTO 3 WAVE TRAINS 

!**   INTERFACE.
!     ----------

!       *CALL* *SEP3TR (F3, IJS, IJL, MIJ, U10NEW, THWNEW,
!                       ESWELL, FSWELL, THSWELL, FSEA,
!                       F1, SWM,
!                       EMTRAIN  ,THTRAIN  ,PMTRAIN)
!          *F3*    - BLOCK OF FULL SPECTRA
!          *IJS*    - INDEX OF FIRST GRIDPOINT
!          *IJL*    - INDEX OF LAST GRIDPOINT
!          *MIJ*    - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *U10NEW* - LATESt WIND SPEED.
!          *THWNEW* - LATEST WIND DIRECTION.
!          *ESWELL* - TOTAL SWELL ENERGY
!          *FSWELL* - TOTAL SWELL MEAN FREQUENCY
!          *THSWELL*- TOTAL SWELL MEAN DIRECTION
!          *FSEA*   - WINDSEA MEAN FREQUENCY
!          *F1*    - SWELL SPECTRA
!          *SWM*    - ORIGINAL SWELL MASK
!                     IT MIGHT NEED TO BE ADJUSTED
!                     THIS IS POSSIBLE BECAUSE WE SUPPLY SWELL SPECTRA
!                     IN WHICH THE WIND SEA HAS BEEN REMOVED. SOME OF
!                     THE ENERGY LEAKING OUT FROM THE WINDSEA SPECTRUM
!                     MIGHT NOT BE PARTITIONED.
!          *EMTRAIN* - PARTITIONED ENERGY
!          *THTRAIN* - PARTITIONED MEAN DIRECTION
!          *PMTRAIN* - PARTITIONED MEAN PERIOD

!     METHOD.
!     -------

!       HANSON AND PHILLIPS 2001

!     EXTERNALS.
!     ----------

!       *FNDPRT*  - COMPUTE PARTITION MASKS

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : NTRAIN
      USE YOWFRED  , ONLY : FR       ,DFIM     ,C        ,              &
     &            TH       ,COSTH    ,SINTH    ,FRIC
      USE YOWICE   , ONLY : FLMIN
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : ZPI      ,G        ,EPSMIN
      USE YOWTEST  , ONLY : IU06
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "fndprt.intfb.h"
#include "semean.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL), INTENT(IN) :: MIJ
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: U10NEW, THWNEW 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: ESWELL,  FSWELL, THSWELL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: FSEA
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F3
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: F1, SWM
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NTRAIN), INTENT(OUT) :: EMTRAIN, THTRAIN, PMTRAIN


      INTEGER(KIND=JWIM), PARAMETER :: NPMAX=20
      INTEGER(KIND=JWIM) :: NPMAX_LOC
      INTEGER(KIND=JWIM) :: IJ, M, K, IP
      INTEGER(KIND=JWIM) :: ISORT, I, IPLOC
      INTEGER(KIND=JWIM) :: IFL, IFH, ITHL, ITHH
      INTEGER(KIND=JWIM) :: KM, KP
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL) :: NPEAK, NPK
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL) :: FRINVMIJ
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL) :: MMIN, MMAX 
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL) :: IPNOW 
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL,NTRAIN) :: IENERGY
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL,NPMAX) :: NFRP, NTHP

      ! relative value above max swell value that is considered above noise level
      REAL(KIND=JWRB), PARAMETER :: XNOISELEVEL=0.005_JWRB
      ! minimum Hs (m) value to keep a partition HSMIN=HSMIN_INTER+HSMIN_SLOPE*PERIOD
      REAL(KIND=JWRB), PARAMETER :: HSMIN_INTER=0.05_JWRB
      REAL(KIND=JWRB), PARAMETER :: HSMIN_SLOPE=-0.0017_JWRB
      REAL(KIND=JWRB) :: THRS
      REAL(KIND=JWRB) :: HSMIN

      REAL(KIND=JWRB) :: COSDIFF
      REAL(KIND=JWRB) :: DELDW
      REAL(KIND=JWRB) :: COSDIR, FRLIMIT 
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: FLLOWEST
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: ENEX, SUMETRAIN
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: ETT
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: ENMAX, FLNOISE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG) :: SPRD
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,0:NPMAX) :: DIR, PER, ENE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NTRAIN) :: TEMPDIR, TEMPPER, TEMPENE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: FL, FLLOW

      LOGICAL :: LLEPSMIN
      LOGICAL :: LLADDPART
      LOGICAL, DIMENSION(IJS:IJL,NTRAIN) :: LPWSECTOR
      LOGICAL, DIMENSION(IJS:IJL,NANG) :: LLCOSDIFF

! ----------------------------------------------------------------------

!*    LOOP THROUGH THE GRID POINTS
!     ----------------------------

      IF (LHOOK) CALL DR_HOOK('SEP3TR',0,ZHOOK_HANDLE)

      DO IJ=IJS,IJL
        FRINVMIJ(IJ)=1.0_JWRB/FR(MIJ(IJ))
      ENDDO

      DO IJ=IJS,IJL
        ENMAX(IJ)=0.0_JWRB
      ENDDO
!     SMOOTH INPUT SPECTRA (in direction only)
      DO M=1,NFRE
        DO K=1,NANG
          KM=K-1
          IF (KM.LT.1) KM=NANG
          KP=K+1
          IF (KP.GT.NANG) KP=1
          DO IJ=IJS,IJL
!           RE-IMPOSE THE WINDSEA MASK
            IF(F1(IJ,K,M).LE.0.0_JWRB) THEN
              FL(IJ,K,M) = 0.0_JWRB
            ELSE
              FL(IJ,K,M) = 0.10_JWRB*(F1(IJ,KM,M)+F1(IJ,KP,M)) + 0.80_JWRB*F1(IJ,K,M) 
              ENMAX(IJ)=MAX(ENMAX(IJ),FL(IJ,K,M))
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        FLNOISE(IJ) = XNOISELEVEL*ENMAX(IJ)
      ENDDO

      DO ISORT=1,NTRAIN
        DO IJ=IJS,IJL
          EMTRAIN(IJ,ISORT) = 0.0_JWRB
          PMTRAIN(IJ,ISORT) = 0.0_JWRB
          IF(ESWELL(IJ).GT.0.0_JWRB) THEN
            THTRAIN(IJ,ISORT) = THSWELL(IJ)
          ELSE
            THTRAIN(IJ,ISORT) = THWNEW(IJ)
          ENDIF
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        NPEAK(IJ)=0
      ENDDO

!*    1. DETERMINATES MAXIMA POSITIONS AND NUMBER
!        ----------------------------------------

      DO K=1,NANG
        DO IJ=IJS,IJL
          COSDIFF=COS(TH(K)-THWNEW(IJ))
          LLCOSDIFF(IJ,K)=(COSDIFF.LT.-0.4_JWRB)
          SPRD(IJ,K)=MAX(0.0_JWRB,COSDIFF)**2
        ENDDO
      ENDDO

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            FLLOW(IJ,K,M)=FLMIN
          ENDDO
        ENDDO
      ENDDO


      DO IJ=IJS,IJL
        OUT: DO M= 2,MIJ(IJ)-1
          IFL    = MAX ( 1 , M-1 )
          IFH    = MIN ( NFRE , M+1 )
          DO K=1,NANG

            FLLOWEST = MAX(FLLOW(IJ,K,M),FLNOISE(IJ))

            IF(FL(IJ,K,M) .GT. FLLOWEST) THEN
              ITHL   = 1 + MOD(NANG+K-2,NANG)
              ITHH   = 1 + MOD(K,NANG)
              IF ( FL(IJ,ITHL,M   ) .GT. 0.0_JWRB .AND.                 &
     &             FL(IJ,ITHH,M   ) .GT. 0.0_JWRB .AND.                 &
     &             FL(IJ,K   ,IFL ) .GT. 0.0_JWRB .AND.                 &
     &             FL(IJ,K   ,IFH ) .GT. 0.0_JWRB .AND.                 &
     &             FL(IJ,ITHL,IFL ) .GT. 0.0_JWRB .AND.                 &
     &             FL(IJ,ITHL,IFH ) .GT. 0.0_JWRB .AND.                 &
     &             FL(IJ,ITHH,IFL ) .GT. 0.0_JWRB .AND.                 &
     &             FL(IJ,ITHH,IFH ) .GT. 0.0_JWRB ) THEN


                IF ( FL(IJ,K,M) .GE. FL(IJ,K   ,IFL ) .AND.             &
     &               FL(IJ,K,M) .GE. FL(IJ,K   ,IFH ) .AND.             &
     &               FL(IJ,K,M) .GE. FL(IJ,ITHL,IFL ) .AND.             &
     &               FL(IJ,K,M) .GE. FL(IJ,ITHL,M   ) .AND.             &
     &               FL(IJ,K,M) .GE. FL(IJ,ITHL,IFH ) .AND.             &
     &               FL(IJ,K,M) .GE. FL(IJ,ITHH,IFL ) .AND.             &
     &               FL(IJ,K,M) .GE. FL(IJ,ITHH,M   ) .AND.             &
     &               FL(IJ,K,M) .GE. FL(IJ,ITHH,IFH ) ) THEN
                    NPEAK(IJ) = NPEAK(IJ) + 1
                    IF(NPEAK(IJ).GT.NPMAX) EXIT OUT
                    NFRP(IJ,NPEAK(IJ)) = M
                    NTHP(IJ,NPEAK(IJ)) = K
                ENDIF
              ENDIF
            ENDIF

          ENDDO

        ENDDO OUT
      ENDDO

      DO IJ=IJS,IJL
        NPEAK(IJ)=MIN(NPEAK(IJ),NPMAX)
      ENDDO

!*    2. GENERATE MASK FOR EACH PARTITION AND COMPUTE STATISTICS
!        -------------------------------------------------------

      CALL FNDPRT(IJS, IJL, NPMAX,                                      &
     &            NPEAK, MIJ, NTHP, NFRP,                               &
     &            FLLOW, LLCOSDIFF, FLNOISE,                            &
     &            FL, SWM,                                              &
     &            ENE, DIR, PER)

      ! UPDATE SWELL SPECTRUM ??????
      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            F1(IJ,K,M)=MAX(F3(IJ,K,M),EPSMIN)*SWM(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO

!     TOTAL ENERGY IN THE UPDATED SWELL SPECTRA
      LLEPSMIN=.FALSE.
      CALL SEMEAN (F1, IJS, IJL, IJS, IJL, ETT, LLEPSMIN)

!     REMOVE PARTITION WITH Hs < HSMIN 
!     MEAN PERIOD > 1/FR(MIJ)
      DO IJ=IJS,IJL
        NPK(IJ)=NPEAK(IJ)
        DO IP=1,NPEAK(IJ)
          HSMIN=HSMIN_INTER+HSMIN_SLOPE*PER(IJ,IP)
          THRS=0.0625_JWRB*HSMIN**2
          IF(ENE(IJ,IP).LT.THRS .OR. PER(IJ,IP).LE.FRINVMIJ(IJ)) THEN
            ENE(IJ,IP)=0.
            DIR(IJ,IP)=0.
            PER(IJ,IP)=0.
            NPK(IJ)=NPK(IJ)-1
          ENDIF
        ENDDO
      ENDDO

!     IF NO SWELL PARTITION BUT TOTAL SWELL EXIST THEN
!     ASSIGN FIRST PARTITION TO TOTAL SWELL
!     IF IT HAS SWELL CHARACTERSITICS
      DO IJ=IJS,IJL
        IF(NPK(IJ).LE.0 .AND. ESWELL(IJ).GT.0.) THEN
          IF(FSWELL(IJ).LT.FSEA(IJ)) THEN
            NPEAK(IJ)=1
            ENE(IJ,1)=ESWELL(IJ)
            DIR(IJ,1)=THSWELL(IJ)
            PER(IJ,1)=1.0_JWRB/FSWELL(IJ)
          ENDIF 
        ENDIF 
      ENDDO

!*    5. SORT PARTITIONS ACCORDING TO ENERGY AND ASSIGN THE FIRST NTRAIN TRAINS
!        WE ONLY NEED THE FIRST NTRAIN TO BE SORTED
!        -----------------------------------------------------------------

      NPMAX_LOC=MAXVAL(NPEAK(:))

      DO ISORT=1,NTRAIN
        DO IJ=IJS,IJL
          IPNOW(IJ)=0
          ENMAX(IJ)=0.0_JWRB
        ENDDO

        DO IP=1,NPMAX_LOC
          DO IJ=IJS,IJL
            IF (ENE(IJ,IP).GT.ENMAX(IJ)) THEN
              IPNOW(IJ) = IP
              ENMAX(IJ) = ENE(IJ,IP)
            ENDIF
          ENDDO
        ENDDO

        DO IJ=IJS,IJL
          EMTRAIN(IJ,ISORT)=ENE(IJ,IPNOW(IJ))
          THTRAIN(IJ,ISORT)=DIR(IJ,IPNOW(IJ))
          PMTRAIN(IJ,ISORT)=PER(IJ,IPNOW(IJ))
          ENE(IJ,IPNOW(IJ))=0.0_JWRB
          IENERGY(IJ,ISORT)=MIN(IPNOW(IJ),1)
        ENDDO
      ENDDO


!*    6. PRESERVE TOTAL ENERGY
!        ---------------------

!     6.1 Distribute extra energy proportionally to swell trains

      DO IJ=IJS,IJL
        SUMETRAIN(IJ)=MAX(EMTRAIN(IJ,1),EPSMIN)
      ENDDO
      DO ISORT=2,NTRAIN
        DO IJ=IJS,IJL
          SUMETRAIN(IJ)=SUMETRAIN(IJ)+EMTRAIN(IJ,ISORT)
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        ENEX(IJ)=MAX((ETT(IJ)-SUMETRAIN(IJ)),0.0_JWRB)/SUMETRAIN(IJ)
      ENDDO

      DO ISORT=1,NTRAIN
        DO IJ=IJS,IJL
          EMTRAIN(IJ,ISORT)=EMTRAIN(IJ,ISORT)+ENEX(IJ)*EMTRAIN(IJ,ISORT)
        ENDDO
      ENDDO

!     RE-ORDER SWELL SYSTEM ACCORDING TO TOTAL ENERGY

      DO ISORT=1,NTRAIN
        DO IJ=IJS,IJL
          TEMPENE(IJ,ISORT)=EMTRAIN(IJ,ISORT)
          TEMPDIR(IJ,ISORT)=THTRAIN(IJ,ISORT)
          TEMPPER(IJ,ISORT)=PMTRAIN(IJ,ISORT)
        ENDDO
      ENDDO

      DO ISORT=1,NTRAIN
        DO IJ=IJS,IJL
          IPNOW(IJ)=0
          ENMAX(IJ)=0.0_JWRB
        ENDDO

        DO IP=1,NTRAIN
          DO IJ=IJS,IJL
            IF (TEMPENE(IJ,IP).GT.ENMAX(IJ)) THEN
              IPNOW(IJ) = IP
              ENMAX(IJ) = TEMPENE(IJ,IP)
            ENDIF
          ENDDO
        ENDDO

        DO IJ=IJS,IJL
          IPLOC=MAX(IPNOW(IJ),1)
          EMTRAIN(IJ,ISORT)=TEMPENE(IJ,IPLOC)
          THTRAIN(IJ,ISORT)=TEMPDIR(IJ,IPLOC)
          PMTRAIN(IJ,ISORT)=TEMPPER(IJ,IPLOC)
          TEMPENE(IJ,IPLOC)=0.0_JWRB
        ENDDO
      ENDDO

!     PREPARE OUTPUT

      DO ISORT=1,NTRAIN
        DO IJ=IJS,IJL
          IF(IENERGY(IJ,ISORT).EQ.0) THEN
            EMTRAIN(IJ,ISORT) = 0.0_JWRB
            THTRAIN(IJ,ISORT) = THWNEW(IJ)
            PMTRAIN(IJ,ISORT) = 0.0_JWRB
          ENDIF
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('SEP3TR',1,ZHOOK_HANDLE)

      END SUBROUTINE SEP3TR
