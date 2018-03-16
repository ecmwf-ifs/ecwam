      SUBROUTINE OUTBC (FL3, IJS, IJL, IG, IU19)

! ----------------------------------------------------------------------

!**** *OUTBC* - OUTPUT OF THE COARSE GRID BOUNDARY VALUES.

!     R. PORTZ     MPI          JANUARY 1991
!     J. BIDLOT    ECMWF        FEBRARY 1996  MESSAGE PASSING

!*    PURPOSE.
!     --------

!        WRITE THE BOUNDARY VALUE OUTPUT FILE(s).

!**   INTERFACE.
!     ----------

!    *CALL* *OUTBC (FL3, IJS, IJL, IG, IU19)*
!      *FL3*     - BLOCK OF SPECTRA.
!      *IJS*     - INDEX OF FIRST GRIDPOINT.
!      *IJL*     - INDEX OF LAST GRIDPOINT.
!      *IG*      - BLOCK NUMBER.
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

      USE YOWPARAM , ONLY : NANG     ,NFRE     ,LL1D
      USE YOWCPBO  , ONLY : NBOUNC   ,IJARC    ,IGARC    ,GBOUNC   ,
     &            IPOGBO
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,AMOWEP   ,AMOSOP   ,
     &            XDELLA   ,ZDELLO
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWMPP   , ONLY : NINF     ,NSUP     , NPROC  , IRANK
      USE YOWSTAT  , ONLY : CDTPRO   ,NPROMA_WAM
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWSPEC  , ONLY : IJ2NEWIJ

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "femean.intfb.h"
#include "mpgatherbc.intfb.h"
#include "sthq.intfb.h"

      INTEGER, INTENT(IN) :: IJS, IJL, IG 
      INTEGER,DIMENSION(GBOUNC), INTENT(IN) :: IU19

      REAL, INTENT(IN) :: FL3(NINF-1:NSUP,NANG,NFRE)

      INTEGER, PARAMETER :: NSCFLD=3

      INTEGER:: IJ, II, NGOU, IX, KX, K ,M 
      INTEGER:: JKGLO, KIJS, KIJL, NPROMA
      INTEGER :: IRECV

      REAL :: XLON, XLAT
      REAL,DIMENSION(IJS:IJL) :: EM, FM, TQ
      REAL,DIMENSION(NBOUNC) :: EMPTS, FMPTS, TQPTS
      REAL,DIMENSION(NBOUNC,NANG,NFRE) :: FLPTS


! ----------------------------------------------------------------------

      NPROMA=NPROMA_WAM


!*    1. COMPUTE MEAN PARAMETERS.
!        ------------------------

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
      DO JKGLO=IJS,IJL,NPROMA
        KIJS=JKGLO
        KIJL=MIN(KIJS+NPROMA-1,IJL)

        CALL FEMEAN (FL3(KIJS:KIJL,:,:), KIJS, KIJL, EM(KIJS),
     &               FM(KIJS))
        CALL STHQ (FL3(KIJS:KIJL,:,:), KIJS, KIJL, TQ(KIJS))
      ENDDO
!$OMP END PARALLEL DO

      IF (ITEST.GE.3) THEN
        WRITE(IU06,*) '      SUB. OUTBC: INTEGRATED',
     &   ' PARAMETERS COMPUTED FOR BOUNDARY POINTS OUTPUT'
        CALL FLUSH(IU06)
      ENDIF

      IRECV=1
      CALL MPGATHERBC(IRECV, IJS, IJL, NSCFLD,
     &                FL3(IJS:IJL,:,:), EM, TQ, FM,
     &                FLPTS, EMPTS, TQPTS, FMPTS)


!*    2. WRITE INFORMATION TO FILE(s) IU19.
!        ----------------------------------

      IF(IRANK.EQ.IRECV) THEN
        DO II=1,GBOUNC
          DO NGOU=IPOGBO(II-1)+1,IPOGBO(II)
            IF(LL1D .OR. NPROC.EQ.1) THEN
              IJ = IJARC(NGOU)
            ELSE
              IJ = IJ2NEWIJ(IJARC(NGOU))
            ENDIF
            IX = IXLG(IJ,IG)
            KX = KXLT(IJ,IG)
            XLON = AMOWEP + REAL(IX-1)*ZDELLO(KX)
            XLAT = AMOSOP + REAL(KX-1)*XDELLA

            WRITE(IU19(II)) XLON, XLAT, CDTPRO,
     &          EMPTS(NGOU), TQPTS(NGOU), FMPTS(NGOU) 
            WRITE(IU19(II)) ((FLPTS(NGOU,K,M),K=1,NANG),M=1,NFRE)
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE OUTBC
