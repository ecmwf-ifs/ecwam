      SUBROUTINE OUTBC (FL1, IJS, IJL, IU19)

! ----------------------------------------------------------------------

!**** *OUTBC* - OUTPUT OF THE COARSE GRID BOUNDARY VALUES.

!     R. PORTZ     MPI          JANUARY 1991
!     J. BIDLOT    ECMWF        FEBRARY 1996  MESSAGE PASSING

!*    PURPOSE.
!     --------

!        WRITE THE BOUNDARY VALUE OUTPUT FILE(s).

!**   INTERFACE.
!     ----------

!    *CALL* *OUTBC (FL1, IJS, IJL, IU19)*
!      *FL1*     - BLOCK OF SPECTRA.
!      *IJS*     - INDEX OF FIRST GRIDPOINT.
!      *IJL*     - INDEX OF LAST GRIDPOINT.
!      *IU19*    - OUTPUT UNIT OF BOUNDARY VALUES.


!     METHOD.
!     -------

!       SEQUENCIAL UNFORMATED WRITE TO UNIT(s).

!     EXTERNALS.
!     ----------

!       *FEMEAN*    - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT.
!       *SEMEAN*    - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.
!       *STHQ*      - COMPUTATION OF MEAN WAVE DIRECTION AT EACH
!                     GRID POINT.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE     ,LL1D
      USE YOWCPBO  , ONLY : NBOUNC   ,IJARC    ,GBOUNC   ,IPOGBO
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMAP   , ONLY : BLK2GLO  ,AMOWEP   ,AMOSOP   ,    &
     &            XDELLA   ,ZDELLO
      USE YOWMPP   , ONLY : NPROC  , IRANK
      USE YOWSTAT  , ONLY : CDTPRO   ,NPROMA_WAM
      USE YOWSPEC  , ONLY : IJ2NEWIJ

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "femean.intfb.h"
#include "mpgatherbc.intfb.h"
#include "sthq.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM),DIMENSION(GBOUNC), INTENT(IN) :: IU19

      REAL(KIND=JWRB), INTENT(IN) :: FL1(IJS:IJL, NANG, NFRE)

      INTEGER(KIND=JWIM), PARAMETER :: NSCFLD=3

      INTEGER(KIND=JWIM):: IJ, II, NGOU, IX, KX, K ,M 
      INTEGER(KIND=JWIM):: JKGLO, KIJS, KIJL, NPROMA
      INTEGER(KIND=JWIM) :: IRECV

      REAL(KIND=JWRB) :: XLON, XLAT
      REAL(KIND=JWRB),DIMENSION(IJS:IJL) :: EM, FM, TQ
      REAL(KIND=JWRB),DIMENSION(NBOUNC) :: EMPTS, FMPTS, TQPTS
      REAL(KIND=JWRB),DIMENSION(NBOUNC,NANG,NFRE) :: FLPTS


! ----------------------------------------------------------------------

      NPROMA=NPROMA_WAM


!*    1. COMPUTE MEAN PARAMETERS.
!        ------------------------

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
      DO JKGLO=IJS,IJL,NPROMA
        KIJS=JKGLO
        KIJL=MIN(KIJS+NPROMA-1,IJL)

        CALL FEMEAN (KIJS, KIJL, FL1(KIJS:KIJL,:,:), EM(KIJS), FM(KIJS) )
        CALL STHQ (KIJS, KIJL, FL1(KIJS:KIJL,:,:), TQ(KIJS))
      ENDDO
!$OMP END PARALLEL DO

      IRECV=1
      CALL MPGATHERBC(IRECV, IJS, IJL, NSCFLD,                          &
     &                FL1, EM, TQ, FM,                                  &
     &                FLPTS, EMPTS, TQPTS, FMPTS)


!*    2. WRITE INFORMATION TO FILE(s) IU19.
!        ----------------------------------

      IF (IRANK == IRECV) THEN
        DO II=1,GBOUNC
          DO NGOU=IPOGBO(II-1)+1,IPOGBO(II)
            IF (LL1D .OR. NPROC == 1) THEN
              IJ = IJARC(NGOU)
            ELSE
              IJ = IJ2NEWIJ(IJARC(NGOU))
            ENDIF
            IX = BLK2GLO(IJ)%IXLG 
            KX = BLK2GLO(IJ)%KXLT 
            XLON = AMOWEP + REAL(IX-1)*ZDELLO(KX)
            XLAT = AMOSOP + REAL(KX-1)*XDELLA

            WRITE(IU19(II)) XLON, XLAT, CDTPRO,                         &
     &          EMPTS(NGOU), TQPTS(NGOU), FMPTS(NGOU) 
            WRITE(IU19(II)) ((FLPTS(NGOU,K,M),K=1,NANG),M=1,NFRE)
          ENDDO
        ENDDO
      ENDIF

      END SUBROUTINE OUTBC
