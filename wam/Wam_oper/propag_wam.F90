      SUBROUTINE PROPAG_WAM (IJS, IJL, FL1)

! ----------------------------------------------------------------------

!**** *PROPAG_WAM* - WAVE PROPGATION

!*    PURPOSE.
!     --------

!     PROPAGATION

!**   INTERFACE.
!     ----------

!     *CALL* *PROPAG_WAM (IJS, IJL, FL1)*
!          *IJS* - INDEX OF FIRST GRIDPOINT.
!          *IJL* - INDEX OF LAST GRIDPOINT.
!          *FL1*  - SPECTRUM



!     METHOD.
!     -------

!     EXTERNALS.
!     ----------


! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCURR  , ONLY : LLCHKCFL ,LLCHKCFLA
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_RED ,NIBLO
      USE YOWSTAT  , ONLY : IPROPAGS ,NPROMA_WAM
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWUBUF  , ONLY : LUPDTWGHT
      USE UNWAM    , ONLY : PROPAG_UNWAM
      USE YOWUNPOOL ,ONLY : LLUNSTR
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "ctuwupdt.intfb.h"

#include "mpexchng.intfb.h"
#include "propags.intfb.h"
#include "propags1.intfb.h"
#include "propags2.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL, NANG, NFRE), INTENT(INOUT) :: FL1

      INTEGER(KIND=JWIM) :: IJ, K, M, J
      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA, MTHREADS
!$    INTEGER,EXTERNAL :: OMP_GET_MAX_THREADS

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
!     Spectra extended with the hallo exchange for the propagation
!     But limited to NFRE_RED frequencies
      REAL(KIND=JWRB), DIMENSION(NINF-1:NSUP,NANG,NFRE_RED) :: FLEXT

      LOGICAL :: L1STCALL

      DATA L1STCALL / .TRUE. /

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('PROPAG_WAM',0,ZHOOK_HANDLE)


      IF(NIBLO.GT.1) THEN

        NPROMA=NPROMA_WAM
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,K,M,IJ) 
        DO JKGLO=IJS,IJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL)
          DO M=1,NFRE_RED
            DO K=1,NANG
              DO IJ=KIJS,KIJL
                FLEXT(IJ,K,M) = FL1(IJ,K,M)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO


        IF (LLUNSTR) THEN 

           CALL PROPAG_UNWAM(FLEXT, FL1)

        ELSE

!          SET THE DUMMY LAND POINTS TO 0.
           FLEXT(NINF-1,:,:) = 0.0_JWRB 

!          OBTAIN INFORMATION AT NEIGHBORING GRID POINTS (HALLO)
           CALL MPEXCHNG(FLEXT,NANG,NFRE_RED)
           IF (ITEST.GE.2) THEN
             WRITE(IU06,*) '   SUB. PROPAG_WAM: MPEXCHNG CALLED' 
             CALL FLUSH (IU06)
           ENDIF


           CALL GSTATS(1430,0)

!          IPROPAGS is the default option (the other options are kept but usually not used)
           IF(IPROPAGS.EQ.2) THEN

             IF(LUPDTWGHT) THEN
               CALL CTUWUPDT(IJS, IJL)
               LUPDTWGHT=.FALSE.
             ENDIF

             MTHREADS=1
!$           MTHREADS=OMP_GET_MAX_THREADS()
             NPROMA=(IJL-IJS+1)/MTHREADS + 1
!$OMP        PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JKGLO,KIJS,KIJL)
             DO JKGLO=IJS,IJL,NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1,IJL)
               CALL PROPAGS2(FLEXT, FL1, IJS, IJL, KIJS, KIJL)
             ENDDO
!$OMP        END PARALLEL DO

             IF (ITEST.GE.2) THEN
               WRITE(IU06,*) '   SUB. PROPAG_WAM: PROPAGS2 CALLED'
               CALL FLUSH (IU06)
             ENDIF

           ELSE IF(IPROPAGS.EQ.1) THEN
             IF(L1STCALL .OR. LLCHKCFLA) LLCHKCFL=.TRUE.
             NPROMA=NPROMA_WAM
!$OMP        PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
             DO JKGLO=IJS,IJL,NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1,IJL)
               CALL PROPAGS1(FLEXT, FL1, IJS, IJL, KIJS, KIJL, L1STCALL)
             ENDDO
!$OMP       END PARALLEL DO
             IF (ITEST.GE.2) THEN
               WRITE(IU06,*) '   SUB. PROPAG_WAM: PROPAGS1 CALLED'
               CALL FLUSH (IU06)
             ENDIF

           ELSE
             IF(L1STCALL .OR. LLCHKCFLA) LLCHKCFL=.TRUE.
             NPROMA=NPROMA_WAM
!$OMP        PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
             DO JKGLO=IJS,IJL,NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1,IJL)
               CALL PROPAGS(FLEXT, FL1, IJS, IJL, KIJS, KIJL, L1STCALL)
             ENDDO
!$OMP        END PARALLEL DO
             IF (ITEST.GE.2) THEN
               WRITE(IU06,*) '   SUB. PROPAG_WAM: PROPAGS CALLED'
               CALL FLUSH (IU06)
             ENDIF
           ENDIF

           CALL GSTATS(1430,1)

        ENDIF  ! end propagation

      ENDIF ! more than one grid point

      L1STCALL=.FALSE.
      LLCHKCFL=.FALSE.

      IF (LHOOK) CALL DR_HOOK('PROPAG_WAM',1,ZHOOK_HANDLE)

      END SUBROUTINE PROPAG_WAM
