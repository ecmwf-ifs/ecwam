SUBROUTINE OUTBS (IJS, IJL, MIJ, FL1, XLLWS)
! ----------------------------------------------------------------------

!**** *OUTBS* - MODEL OUTPUT FROM BLOCK TO FILE, PRINTER AND COMMON.

!*    PURPOSE.
!     --------

!       CONTROL OUTPUT OF WAVE AND WIND FIELDS (except spectrum).

!**   INTERFACE.
!     ----------
!      *CALL*OUTBS (IJS, IJL, MIJ, FL1, XLLWS)
!      *IJS:IJL - FIRST DIMENSION OF ARRAYS MIJ, FL1, XLLWS.
!      *MIJ*    - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!      *FL1*    - INPUT SPECTRUM.
!      *XLLWS*  - WINDSEA MASK FROM INPUT SOURCE TERM

!     EXTERNALS.
!     ----------

!       *OUTBLOCK*  - GET ALL OUTPUT PARAMETERS
!   
!     METHOD.
!     -------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : JPPFLAG  ,NIPRMOUT    ,BOUT
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWGRID  , ONLY : IJSLOC   ,IJLLOC
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSHAL  , ONLY : DEPTH    ,WAVNUM      ,CINV        ,CGROUP
      USE YOWSTAT  , ONLY : NPROMA_WAM
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "outblock.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL), INTENT(IN) :: MIJ
           
      REAL(KIND=JWRB), DIMENSION(IJS:IJL, NANG, NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(IJS:IJL, NANG, NFRE), INTENT(IN) :: XLLWS

      INTEGER(KIND=JWIM) :: M, IJ, JKGLO, KIJS, KIJL, NPROMA

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('OUTBS',0,ZHOOK_HANDLE)

!*    1. COMPUTE MEAN PARAMETERS.
!        ------------------------

      NPROMA=NPROMA_WAM

!     COMPUTE MEAN PARAMETERS

      CALL GSTATS(1502,0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
      DO JKGLO = IJSLOC, IJLLOC, NPROMA
        KIJS=JKGLO
        KIJL=MIN(KIJS+NPROMA-1,IJLLOC)
        CALL OUTBLOCK(IJS, IJL, KIJS, KIJL, MIJ(KIJS),                          &
     &                FL1(IJS,1,1), XLLWS(IJS,1,1),                             &
     &                WAVNUM(IJ,1), CINV(IJS,1), CGROUP(IJS,1),                 &
     &                DEPTH(KIJS),                                              &
     &                BOUT(IJS,1))
      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1502,1)

      IF (ITEST.GE.3) THEN
          WRITE(IU06,*) '      SUB. OUTBS: INTEGRATED PARAMETERS COMPUTED FOR OUTPUT'
      ENDIF

      IF (LHOOK) CALL DR_HOOK('OUTBS',1,ZHOOK_HANDLE)

END SUBROUTINE OUTBS
