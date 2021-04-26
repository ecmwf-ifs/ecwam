      SUBROUTINE PROPAG_WAM (FL1)

! ----------------------------------------------------------------------

!**** *PROPAG_WAM* - PROPGATION ON STRUCTURED GRID. 

!*    PURPOSE.
!     --------

!     PROPAGATION

!**   INTERFACE.
!     ----------

!     *CALL* *PROPAG_WAM (FL1)*
!          *FL1*  - SPECTRUM


!     METHOD.
!     -------

!     EXTERNALS.
!     ----------


! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCURR  , ONLY : LLCHKCFL ,LLCHKCFLA
      USE YOWFRED  , ONLY : FR5      ,FRM5
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NIBLO
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

      REAL(KIND=JWRB), DIMENSION(NINF-1:NSUP,NANG,NFRE), INTENT(INOUT) :: FL1

      INTEGER(KIND=JWIM), PARAMETER :: NFRE_PRO=29
 
      INTEGER(KIND=JWIM) :: IG
      INTEGER(KIND=JWIM) :: IJ, K, M, J
      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA, MTHREADS
!$    INTEGER,EXTERNAL :: OMP_GET_MAX_THREADS

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS(1):IJL(1),NANG,NFRE_PRO) :: FLNEW

      LOGICAL :: L1STCALL

      DATA L1STCALL / .TRUE. /

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('PROPAG_WAM',0,ZHOOK_HANDLE)

      IG=1

      IF(NIBLO.GT.1) THEN

        IF (LLUNSTR) THEN 

           CALL PROPAG_UNWAM(FL1, FLNEW)

        ELSE

!          SET THE DUMMY LAND POINTS TO 0.
           FL1(NINF-1,:,:) = 0.0_JWRB 

!          OBTAIN INFORMATION AT NEIGHBORING GRID POINTS
           CALL MPEXCHNG(FL1,NANG,NFRE)
           IF (ITEST.GE.2) THEN
             WRITE(IU06,*) '   SUB. PROPAG_WAM: MPEXCHNG CALLED' 
             CALL FLUSH (IU06)
           ENDIF


           CALL GSTATS(1430,0)

           IF(IPROPAGS.EQ.1) THEN
             IF(L1STCALL .OR. LLCHKCFLA) LLCHKCFL=.TRUE.
             NPROMA=NPROMA_WAM
!$OMP        PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
             DO JKGLO=IJS(IG),IJL(IG),NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
               CALL PROPAGS1(FL1,FLNEW(KIJS:KIJL,:,:),KIJS,KIJL,L1STCALL)
             ENDDO
!$OMP       END PARALLEL DO
             IF (ITEST.GE.2) THEN
               WRITE(IU06,*) '   SUB. PROPAG_WAM: PROPAGS1 CALLED'
               CALL FLUSH (IU06)
             ENDIF

           ELSEIF(IPROPAGS.EQ.2) THEN

             IF(LUPDTWGHT) THEN
               CALL CTUWUPDT(IJS(IG), IJL(IG))
               LUPDTWGHT=.FALSE.
             ENDIF

             MTHREADS=1
!$           MTHREADS=OMP_GET_MAX_THREADS()
             NPROMA=(IJL(IG)-IJS(IG)+1)/MTHREADS + 1
!$OMP        PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JKGLO,KIJS,KIJL)
             DO JKGLO=IJS(IG),IJL(IG),NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
               CALL PROPAGS2(FL1,FLNEW,NFRE_PRO,IJS(1),IJL(1),KIJS,KIJL)
             ENDDO
!$OMP        END PARALLEL DO
             IF (ITEST.GE.2) THEN
               WRITE(IU06,*) '   SUB. PROPAG_WAM: PROPAGS2 CALLED'
               CALL FLUSH (IU06)
             ENDIF
           ELSE
             IF(L1STCALL .OR. LLCHKCFLA) LLCHKCFL=.TRUE.
             NPROMA=NPROMA_WAM
!$OMP        PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
             DO JKGLO=IJS(IG),IJL(IG),NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
               CALL PROPAGS(FL1,FLNEW(KIJS:KIJL,:,:),KIJS,KIJL,L1STCALL)
             ENDDO
!$OMP        END PARALLEL DO
             IF (ITEST.GE.2) THEN
               WRITE(IU06,*) '   SUB. PROPAG_WAM: PROPAGS CALLED'
               CALL FLUSH (IU06)
             ENDIF
           ENDIF
           CALL GSTATS(1430,1)

        ENDIF  ! end propagation



!       UPDATING SPECTRA

!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,K,M,IJ) 
        DO JKGLO=IJS(IG),IJL(IG),NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
          DO M=1,NFRE_PRO
            DO K=1,NANG
              DO IJ=KIJS,KIJL
                FL1(IJ,K,M) = FLNEW(IJ,K,M)
              ENDDO
            ENDDO
          ENDDO
          DO M=NFRE_PRO+1,NFRE
            DO K=1,NANG
              DO IJ=KIJS,KIJL
                FL1(IJ,K,M) = FL1(IJ,K,NFRE_PRO)*FR5(NFRE_PRO)*FRM5(M)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO


      ENDIF ! more than one grid point

      L1STCALL=.FALSE.
      LLCHKCFL=.FALSE.

      IF (LHOOK) CALL DR_HOOK('PROPAG_WAM',1,ZHOOK_HANDLE)

      END SUBROUTINE PROPAG_WAM
