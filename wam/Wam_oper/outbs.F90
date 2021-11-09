SUBROUTINE OUTBS (IJS, IJL, MIJ, FL1, XLLWS,                  &
 &                WVPRPT, WVENVI, FF_NOW, INTFLDS, NEMO2WAM,  &
 &                BOUT)
! ----------------------------------------------------------------------

!**** *OUTBS* - MODEL OUTPUT FROM BLOCK TO FILE, PRINTER AND COMMON.

!*    PURPOSE.
!     --------

!       CONTROL OUTPUT OF WAVE AND WIND FIELDS (except spectrum).

!**   INTERFACE.
!     ----------
!      *CALL*OUTBS (IJS, IJL, MIJ, FL1, XLLWS,
!                   WVPRPT, WVENVI, FF_NOW, INTFLDS, NEMO2WAM, BOUT)
!      *IJS:IJL  - FIRST DIMENSION OF ARRAYS MIJ, FL1, XLLWS, BOUT.
!      *MIJ*     - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!      *FL1*     - INPUT SPECTRUM.
!      *XLLWS*   - WINDSEA MASK FROM INPUT SOURCE TERM
!      *WVENVI*  - WAVE ENVIRONMENT (depth, currents,...)
!      *FF_NOW*  - FORCING FIELDS
!      *INTFLDS* - INTEGRATED/DERIVED PARAMETERS
!      *NEMO2WAM*- FIELDS FRON OCEAN MODEL to WAM
!      *BOUT*    - OUTPUT PARAMETERS BUFFER 



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
      USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FREQUENCY, FORCING_FIELDS,  &
     &                         INTGT_PARAM_FIELDS, OCEAN2WAVE

      USE YOWCOUT  , ONLY : JPPFLAG  ,NIPRMOUT
      USE YOWCOUP  , ONLY : LLNORMWAMOUT
      USE YOWGRID  , ONLY : IJSLOC   ,IJLLOC
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS  , ONLY : ZMISS
      USE YOWSTAT  , ONLY : NPROMA_WAM
      USE YOWUNPOOL ,ONLY : LLUNSTR

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "outblock.intfb.h"
#include "outwnorm.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL), INTENT(IN) :: MIJ
      REAL(KIND=JWRB), DIMENSION(IJS:IJL, NANG, NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(IJS:IJL, NANG, NFRE), INTENT(IN) :: XLLWS
      TYPE(FREQUENCY), DIMENSION(IJS:IJL,NFRE), INTENT(IN) :: WVPRPT
      TYPE(ENVIRONMENT), DIMENSION(IJS:IJL), INTENT(IN) :: WVENVI
      TYPE(FORCING_FIELDS), DIMENSION(IJS:IJL), INTENT(IN) :: FF_NOW
      TYPE(INTGT_PARAM_FIELDS), DIMENSION(IJS:IJL), INTENT(IN) :: INTFLDS
      TYPE(OCEAN2WAVE), DIMENSION(IJS:IJL), INTENT(IN) :: NEMO2WAM
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NIPRMOUT), INTENT(OUT) :: BOUT

      INTEGER(KIND=JWIM) :: M, IJ, JKGLO, KIJS, KIJL, NPROMA

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      LOGICAL :: LDREPROD

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('OUTBS',0,ZHOOK_HANDLE)

!*    1. COMPUTE MEAN PARAMETERS.
!        ------------------------

      NPROMA=NPROMA_WAM

!     COMPUTE MEAN PARAMETERS

      CALL GSTATS(1502,0)

      IF (LLUNSTR) THEN
!!!!     the openMP option with unstructured grid will need to be coded better
!!!!     but currently not used operationally....
        BOUT(:,:) = ZMISS

!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO = IJSLOC, IJLLOC, NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJLLOC)
          CALL OUTBLOCK(KIJS, KIJL, MIJ(KIJS),                      &
     &                  FL1(KIJS:KIJL,:,:), XLLWS(KIJS:KIJL,:,:),   &
     &                  WVPRPT(KIJS:KIJL,:),                        &
     &                  WVENVI(KIJS), FF_NOW(KIJS), INTFLDS(KIJS),  &
     &                  NEMO2WAM(KIJS),                             &
     &                  BOUT(KIJS:KIJL,:))
        ENDDO
!$OMP   END PARALLEL DO

      ELSE

!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO = IJS, IJL, NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL)
          CALL OUTBLOCK(KIJS, KIJL, MIJ(KIJS),                      &
     &                  FL1(KIJS:KIJL,:,:), XLLWS(KIJS:KIJL,:,:),   &
     &                  WVPRPT(KIJS:KIJL,:),                        &
     &                  WVENVI(KIJS), FF_NOW(KIJS), INTFLDS(KIJS),  &
     &                  NEMO2WAM(KIJS),                             &
     &                  BOUT(KIJS:KIJL,:))
        ENDDO
!$OMP   END PARALLEL DO

      ENDIF
      CALL GSTATS(1502,1)

!     PRINT OUT NORMS
!!!1 to do: decide if there are cases where we might want LDREPROD false
      LDREPROD=.TRUE.
      IF (LLNORMWAMOUT) CALL OUTWNORM(LDREPROD, IJS, IJL, BOUT)


      IF (LHOOK) CALL DR_HOOK('OUTBS',1,ZHOOK_HANDLE)

END SUBROUTINE OUTBS
