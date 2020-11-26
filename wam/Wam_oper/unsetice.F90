      SUBROUTINE UNSETICE(FL,IJS,IJL,U10OLD,THWOLD)
! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      AUGUST 2006 

!*    PURPOSE.
!     --------
!     TRY TO CREATE A JONSWAP SPECTRUM AT POINTS WHICH WERE ICE IN THE
!     PREVIOUS RUN AND WHICH ARE NOW OPEN SEA.
!     ALSO REINSTATE THE MIMIMUM ENERGY LEVEL.

!**   INTERFACE.
!     ----------
!     *CALL* *UNSETICE(FL,IJS,IJL,U10OLD,THWOLD)
!     *FL*        ARRAY CONTAINING THE SPECTRA CONTRIBUTION ON EACH PE
!     *IJS*       INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!     *IJL*       INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!     *U10OLD*    WIND SPEED. (used with fetch law to fill empty 
!                 sea points)
!     *THWOLD*    WIND DIRECTION (RADIANS).


!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!     NONE

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRID  , ONLY : DELPHI
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE YOWJONS  , ONLY : AJONS    ,BJONS    ,DJONS    ,EJONS
      USE YOWICE   , ONLY : FLMIN    ,LICERUN  ,LMASKICE ,CICOVER  ,CITHRSH
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,CLDOMAIN
      USE YOWPCONS , ONLY : PI       ,G        ,R        ,EPSMIN
      USE YOWSHAl  , ONLY : DEPTH    ,EMAXDPT
      USE YOWSTAT  , ONLY : ISHALLO  ,LBIWBK
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWTEXT  , ONLY : LRESTARTED
      USE MPL_MODULE

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "sdepthlim.intfb.h"
#include "jonswap.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM) :: IJ, K, M

      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: U10OLD,THWOLD
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: FL


      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: CILIMIT
      REAL(KIND=JWRB) :: GAMMA, SA, SB, FETCH, GX, FPMAX, ZDP
      REAL(KIND=JWRB) :: U10, GXU, UG, FLLOWEST 
      REAL(KIND=JWRB) :: STK(NANG)
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: FPK, ALPHAV
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE) :: ET
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG) :: SPRD, COS2NOISE

      LOGICAL :: LICE2SEA(IJS:IJL)

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('UNSETICE',0,ZHOOK_HANDLE)

      IF(.NOT.LRESTARTED .AND. CLDOMAIN.NE.'s') THEN

        ZDP=2.0_JWRB/PI
        DO K=1,NANG
          DO IJ=IJS,IJL
            COS2NOISE(IJ,K)=MAX(0.0_JWRB,COS(TH(K)-THWOLD(IJ)))**2
            SPRD(IJ,K)=ZDP*COS2NOISE(IJ,K)
          ENDDO
        ENDDO

        IF(LICERUN .AND. LMASKICE) THEN
          CILIMIT=CITHRSH
        ELSE
!         WHEN WE ARE NOT IMPOSING A FIXED ICE BOUNDARY THEN
!         ONLY POINTS WITH NO ICE WILL BE RESET
          CILIMIT=0.0_JWRB
        ENDIF

!       TRY TO CREATE A JONSWAP SPECTRUM AT POINTS WHICH WERE ICE IN THE
!       PREVIOUS RUN AND WHICH ARE NOW OPEN SEA.
!       THEY ARE CHARCTERISED BY HAVING NO WAVE ENERGY BESIDE NOISE
!       AND HAVE AN CICOVER > CITHRSH .
        DO IJ=IJS,IJL
          LICE2SEA(IJ)=.TRUE.
        ENDDO
        DO K=1,NANG
          DO M=1,NFRE
            DO IJ=IJS,IJL
               IF(FL(IJ,K,M) .GT. EPSMIN) LICE2SEA(IJ) = .FALSE. 
            ENDDO
          ENDDO
        ENDDO

        GAMMA=3.0_JWRB
        SA=7.0E-02_JWRB
        SB=9.0E-02_JWRB
        FPMAX=FR(NFRE-1)

        DO IJ=IJS,IJL
          IF(LICE2SEA(IJ) .AND. CICOVER(IJ,1).LE.CILIMIT) THEN
            IF(DEPTH(IJ,1).LE.10.0_JWRB) THEN
              FETCH=MIN(0.5_JWRB*DELPHI,10000.0_JWRB)
            ELSEIF(DEPTH(IJ,1).LE.50.0_JWRB) THEN
              FETCH=MIN(DELPHI,50000.0_JWRB)
            ELSEIF(DEPTH(IJ,1).LE.150.0_JWRB) THEN
              FETCH=MIN(2*DELPHI,100000.0_JWRB)
            ELSEIF(DEPTH(IJ,1).LE.250.0_JWRB) THEN
              FETCH=MIN(3*DELPHI,200000.0_JWRB)
            ELSE
              FETCH=MIN(5*DELPHI,250000.0_JWRB)
            ENDIF

!           FIND PEAK PERIOD AND ENERGY LEVEL FROM FETCH LAW
!           THE SAME FORMULATION AS IN SUBROUTINE PEAK IS USED.
            IF (U10OLD(IJ) .GT. 0.1E-08_JWRb ) THEN
              GX = G * FETCH
              U10 = U10OLD(IJ)
              GXU = GX/(U10*U10)
              UG = G/U10
              FPK(IJ) = AJONS * GXU ** DJONS
              FPK(IJ) = MAX(0.13_JWRB, FPK(IJ))
              FPK(IJ) = MIN(FPK(IJ), FPMAX/UG)
              ALPHAV(IJ) = BJONS * FPK(IJ)** EJONS
              ALPHAV(IJ) = MAX(ALPHAV(IJ), 0.0081_JWRB)
              FPK(IJ) = FPK(IJ)*UG
            ELSE
              FPK(IJ)=0.
              ALPHAV(IJ)=0.0_JWRB
            ENDIF
          ELSE
            FPK(IJ)=FPMAX 
            ALPHAV(IJ)=0.0_JWRB
          ENDIF
        ENDDO

        CALL JONSWAP(ALPHAV, GAMMA, SA, SB, FPK, IJS, IJL, ET)

        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=IJS,IJL
              FL(IJ,K,M) = MAX(FL(IJ,K,M),ET(IJ,M)*SPRD(IJ,K))
              FLLOWEST = FLMIN*COS2NOISE(IJ,K)
              FL(IJ,K,M) = MAX(FL(IJ,K,M),FLLOWEST)
            ENDDO
          ENDDO
        ENDDO

        IF(ISHALLO.NE.1 .AND. LBIWBK) THEN
          CALL SDEPTHLIM(IJS,IJL,EMAXDPT(IJS),FL)
        ENDIF

        IF (ITEST.GE.2) THEN
         IF(LICERUN) THEN
         WRITE(IU06,*) ' UNSETICE: SPECTRA INITIALISED OVER OLD SEA ICE'
         ENDIF
         WRITE(IU06,*) ' UNSETICE: NOISE LEVEL REIMPOSED'
         WRITE(IU06,*) ' '
        ENDIF

      ENDIF

      IF (LHOOK) CALL DR_HOOK('UNSETICE',1,ZHOOK_HANDLE)

      END SUBROUTINE UNSETICE
