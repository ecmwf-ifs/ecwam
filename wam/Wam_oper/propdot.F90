      SUBROUTINE PROPDOT(IJS, IJL, THDC, THDD, SDOT)

! ----------------------------------------------------------------------

!**** *PROPDOT* - PROPAGATION DOT TERMS FROM DEPTH AND CURRENT GRADIENT.

!     H. GUNTHER   GKSS/ECMWF   17/02/91
!     J. BIDLOT    ECMWF        FEBRUARY 97 MESSAGE PASSING
!     J. BIDLOT    ECMWF        JANUARY 2004 REMOVE MULTI-BLOCK OPTION.

!*    PURPOSE.
!     --------

!       COMPUTATION OF COMMON REFDOT FOR PROPAGATION.

!**   INTERFACE.
!     ----------

!       *CALL* *PROPDOT*(IJS, IJL, THDC, THDD, SDOT)
!         THDC and THDD ARRAYS TO KEEP DEPTH AND CURRENT REFRACTION FOR THETA DOT
!         AND SDOT FOR SIGMA DOT

!     METHOD.
!     -------

!       THE DEPTH AND CURRENT GRADIENTS ARE COMPUTED

!     EXTERNALS.
!     ----------

!       *GRADI*     - COMPUTES DEPTH AND CURRENT GRADIENTS.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCURR  , ONLY : U        ,V
      USE YOWFRED  , ONLY : COSTH    ,SINTH
      USE YOWGRID  , ONLY : COSPHM1
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSHAL  , ONLY : TCGOND   ,TFAK     ,TSIHKD   ,INDEP 
      USE YOWSTAT  , ONLY : ICASE    ,ISHALLO  ,IREFRA
      USE YOWTEST  , ONLY : IU06
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "gradi.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL, NANG), INTENT(OUT) :: THDC, THDD
      REAL(KIND=JWRB), DIMENSION(IJS:IJL, NANG, NFRE), INTENT(OUT) :: SDOT


      INTEGER(KIND=JWIM) :: IG
      INTEGER(KIND=JWIM) :: IJ, K, M

      REAL(KIND=JWRB) :: CD, SD, SS, SC, CC
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: DDPHI, DDLAM, DUPHI, DULAM
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: DVPHI, DVLAM, DCO, OMDD

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('PROPDOT',0,ZHOOK_HANDLE)

      IG = 1

!*    2.2 DEPTH AND CURRENT GRADIENTS.
!         ----------------------------

        CALL GRADI (IJS, IJL, IREFRA, DDPHI, DDLAM, DUPHI,    &
     &              DULAM, DVPHI, DVLAM)

!*    2.3 COSINE OF LATITUDES IF SPHERICAL PROPAGATION.
!         ---------------------------------------------

        IF (ICASE.EQ.1) THEN
          DO IJ = IJS,IJL
            DCO(IJ) = COSPHM1(IJ,IG)
          ENDDO
        ELSE
          DO IJ = IJS,IJL
            DCO(IJ) = 1.0_JWRB
          ENDDO
        ENDIF

!*    2.4 DEPTH GRADIENT PART OF SIGMA DOT.
!         ---------------------------------

        IF (ISHALLO.NE.1 ) THEN
          IF (IREFRA.EQ.3) THEN
            DO IJ = IJS,IJL
              OMDD(IJ) = V(IJ,IG)*DDPHI(IJ) + U(IJ,IG)*DDLAM(IJ)*DCO(IJ)
            ENDDO
          ELSEIF (IREFRA.EQ.2) THEN
            DO IJ = IJS,IJL
              OMDD(IJ) = 0.0_JWRB
            ENDDO
          ENDIF
        ENDIF

!*    2.5. LOOP OVER DIRECTIONS.
!          ---------------------

        DO K=1,NANG
          SD = SINTH(K)
          CD = COSTH(K)

!*    2.5.1. DEPTH GRADIENT OF THETA DOT.
!            ----------------------------

          IF (ISHALLO.NE.1) THEN
            IF (IREFRA.EQ.1 .OR. IREFRA.EQ.3) THEN
              DO IJ = IJS,IJL
                THDD(IJ,K) = SD*DDPHI(IJ) - CD*DDLAM(IJ)*DCO(IJ)
              ENDDO
            ELSE
              DO IJ = IJS,IJL
                THDD(IJ,K) = 0.0_JWRB
              ENDDO
            ENDIF
          ENDIF

!*    2.5.2 SIGMA DOT AND THETA DOT PART FROM CURRENT GRADIENT.
!           ---------------------------------------------------

          IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN
            SS  = SD**2
            SC  = SD*CD
            CC  = CD**2
            DO IJ = IJS,IJL
              SDOT(IJ,K,NFRE) = -SC*DUPHI(IJ) - CC*DVPHI(IJ)            &
     &                        - (SS*DULAM(IJ) + SC*DVLAM(IJ))*DCO(IJ)
              THDC(IJ,K) =  SS*DUPHI(IJ) + SC*DVPHI(IJ)                 &
     &                    - (SC*DULAM(IJ) + CC*DVLAM(IJ))*DCO(IJ)
            ENDDO

!*    2.5.3 LOOP OVER FREQUENCIES.
!           ----------------------

            IF (ISHALLO.NE.1) THEN
              DO M=1,NFRE
                DO IJ=IJS,IJL
                  SDOT(IJ,K,M) = (SDOT(IJ,K,NFRE)*TCGOND(INDEP(IJ),M)   &
     &             + OMDD(IJ)*TSIHKD(INDEP(IJ),M))                      &
     &             * TFAK(INDEP(IJ),M)
                ENDDO

!*    BRANCH BACK TO 2.5.3 FOR NEXT FREQUENCY.

              ENDDO
            ENDIF
          ENDIF

!*    BRANCH BACK TO 2.5 FOR NEXT DIRECTION.

        ENDDO

      IF (LHOOK) CALL DR_HOOK('PROPDOT',1,ZHOOK_HANDLE)

      END SUBROUTINE PROPDOT
