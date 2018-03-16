      SUBROUTINE SPECTRA (IJS, IJL, IG, FL1)

! ----------------------------------------------------------------------

!**** *SPECTRA*  - COMPUTATION OF 2-D SPECTRA FOR ONE BLOCK.

!     S. HASSELMANN  - JULY 1990
!     H. GUNTHER     - DECEMBER 1990   MODIFIED FOR CYCLE_4.
!     J. BIDLOT      - FEBRUARY 1996   MESSAGE PASSING

!*    PURPOSE.
!     --------

!       INITIALISATION OF A BLOCK BY 2-D SPECTRA.

!**   INTERFACE.
!     ----------

!       *CALL* *SPECTRA (IJS, IJL, IG, FL1)*
!          *IJS*     INTEGER  FIRST POINT IN BLOCK.
!          *IJL*     INTEGER  LAST  POINT IN BLOCK.
!          *IG*      INTEGER  BLOCK NUMBER.
!          *FL1*      REAL      2-D SPECTRUM FOR EACH GRID POINT

!     METHOD.
!     -------
!       1-D JONSWAP SPECTRA AND COSINE**2 SPREADING FUNCTIONS ARE
!       COMPUTED FROM GIVEN WINDS AND PARAMETERS AT EACH GRID POINT.
!       THE 1-D SPECTRA ARE SPREAD OVER THE DIRECTIONS BY MULTIPLICATION
!       WITH THE SPREADING FUNCTION.

!     EXTERNALS.
!     ----------

!       *SPT*       - COMPUTATION OF COS**2 SPREADING FUNCTION.


!     REFERENCES.
!     -----------

!       K.HASSELMAN,D.B.ROSS,P.MUELLER AND W.SWELL
!          A PARAMETRIC WAVE PREDICTION MODEL
!          JOURNAL OF PHYSICAL OCEANOGRAPHY, VOL. 6, NO. 2, MARCH 1976

! ----------------------------------------------------------------------

      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWJONS  , ONLY : FP       ,ALPHJ    ,THES     ,GAMMA    ,
     &            SA       ,SB
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "jonswap.intfb.h"
#include "spr.intfb.h"

      INTEGER, INTENT(IN) :: IJS, IJL
      INTEGER, INTENT(IN) :: IG
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE), INTENT(OUT) :: FL1

      INTEGER :: IJ, K, M

      REAL, DIMENSION(NANG) :: STH
      REAL, DIMENSION(IJS:IJL,NANG) :: ST
      REAL, DIMENSION(IJS:IJL,NFRE) :: ET


! ----------------------------------------------------------------------


!*    1. COMPUTE JONSWAP SPECTRUM.
!        -------------------------

      CALL JONSWAP (ALPHJ(IJS), GAMMA, SA, SB, FP(IJS), IJS, IJL, ET) 

! ----------------------------------------------------------------------

!*    2. COMPUTATION OF SPREADING FUNCTION.
!        ----------------------------------

      DO IJ=IJS,IJL
        CALL SPR (NANG, THES(IJ), TH, STH)
        DO K=1,NANG
          ST(IJ,K) = STH(K)
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    3. COMPUTATION OF 2-D SPECTRUM.
!        ----------------------------

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            FL1(IJ,K,M) = ET(IJ,M)*ST(IJ,K)
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE SPECTRA
