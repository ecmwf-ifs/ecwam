SUBROUTINE CTUWUPDT (IJS, IJL)

! ----------------------------------------------------------------------

!**** *CTUWUPDT*

!*    PURPOSE.
!     --------

!     COMPUTES THE WEIGHTS FOR THE CTU ADVECTION SCHEME (IPROPAGS=2).

!**   INTERFACE.
!     ----------

!     *CALL* *CTUWUPDT (IJS, IJL)*

!     METHOD.
!     -------

! -------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWCURR  , ONLY : LLCFLCUROFF
USE YOWREFD  , ONLY : THDD     ,THDC     ,SDOT
USE YOWMPP   , ONLY : IRANK    ,NPROC
USE YOWPARAM , ONLY : NANG     ,NFRE_RED
USE YOWSTAT  , ONLY : IREFRA   ,NPROMA_WAM 
USE YOWTEST  , ONLY : IU06     ,ITEST
USE YOWUBUF  , ONLY : SUMWN    ,                                            &
&                     WLATN    ,WLONN    ,WCORN    ,WKPMN    ,WMPMN    ,    &
&                     LLWLATN  ,LLWLONN  ,LLWCORN  ,LLWKPMN  ,LLWMPMN

USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      
! ----------------------------------------------------------------------
      IMPLICIT NONE

#include "abort1.intfb.h"
#include "ctuw.intfb.h"

INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL

INTEGER(KIND=JWIM) :: IJ, K, M, J, ICALL
INTEGER(KIND=JWIM) :: IC, ICL, ICR
INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA

REAL(KIND=JWRB) :: ZHOOK_HANDLE

LOGICAL :: LL2NDCALL
LOGICAL, DIMENSION(IJS:IJL) :: LCFLFAIL

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CTUWUPDT',0,ZHOOK_HANDLE)

! THE CTU SCHEME IS USED, COMPUTE THE WEIGHTS

IF(.NOT. ALLOCATED(SUMWN)) ALLOCATE(SUMWN(IJS:IJL,NANG,NFRE_RED))
IF(.NOT. ALLOCATED(WLATN)) ALLOCATE(WLATN(IJS:IJL,NANG,NFRE_RED,2,2))
IF(.NOT. ALLOCATED(LLWLATN)) ALLOCATE(LLWLATN(NANG,NFRE_RED,2,2))

IF(.NOT. ALLOCATED(WLONN)) ALLOCATE(WLONN(IJS:IJL,NANG,NFRE_RED,2))
IF(.NOT. ALLOCATED(LLWLONN))  ALLOCATE(LLWLONN(NANG,NFRE_RED,2))

IF(.NOT. ALLOCATED(WCORN)) ALLOCATE(WCORN(IJS:IJL,NANG,NFRE_RED,4,2))
IF(.NOT. ALLOCATED(LLWCORN)) ALLOCATE(LLWCORN(NANG,NFRE_RED,4,2))

IF(.NOT. ALLOCATED(WKPMN)) ALLOCATE(WKPMN(IJS:IJL,NANG,NFRE_RED,-1:1))
IF(.NOT. ALLOCATED(LLWKPMN)) ALLOCATE(LLWKPMN(NANG,NFRE_RED,-1:1))

IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN
  IF(.NOT. ALLOCATED(WMPMN)) ALLOCATE(WMPMN(IJS:IJL,NANG,NFRE_RED,-1:1))
  IF(.NOT. ALLOCATED(LLWMPMN)) ALLOCATE(LLWMPMN(NANG,NFRE_RED,-1:1))
ENDIF

NPROMA=NPROMA_WAM

!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL,ICALL,IJ,LL2NDCALL)
DO JKGLO=IJS,IJL,NPROMA
  KIJS=JKGLO
  KIJL=MIN(KIJS+NPROMA-1,IJL)
  ICALL = 1
  CALL CTUW(KIJS,KIJL,LCFLFAIL(KIJS),ICALL)

! WHEN SURFACE CURRENTS ARE USED AND LLCFLCUROFF IS TRUE
! THEN TRY TO SATISFY THE CFL CONDITION WITHOUT THE CURRENTS
! IF IT WAS VIOLATE IN THE FIRST PLACE
  IF (LLCFLCUROFF .AND. (IREFRA.EQ.2 .OR. IREFRA.EQ.3)) THEN
    LL2NDCALL=.FALSE.
    DO IJ=KIJS,KIJL
       IF(LCFLFAIL(IJ)) THEN
         LL2NDCALL=.TRUE.
         EXIT
       ENDIF
    ENDDO
    IF(LL2NDCALL) THEN
       ICALL = 2
       CALL CTUW(KIJS,KIJL,LCFLFAIL(KIJS),ICALL)
    ENDIF
  ENDIF
ENDDO
!$OMP  END PARALLEL DO

IF (ITEST.GE.2) THEN
  WRITE(IU06,*) '   SUB. CTUWUPDT: CTUW CALLED'
  CALL FLUSH (IU06)
ENDIF

DO IJ=IJS,IJL
  IF(LCFLFAIL(IJ)) THEN
    CALL FLUSH (IU06)
    WRITE(0,*) '!!! ********************************* !!'
    WRITE(0,*) '!!! WAVE MODEL HAS ABORTED !!!'
    WRITE(0,*) '!!! FOLLOWING CFL CRITERION VIOLATION !!'
    WRITE(0,*) '!!! ON PE ',IRANK
    WRITE(0,*) '!!! ********************************* !!'
    WRITE(IU06,*) '!!! ********************************* !!'
    WRITE(IU06,*) '!!! WAVE MODEL HAS ABORTED !!!'
    WRITE(IU06,*) '!!! FOLLOWING CFL CRITERION VIOLATION !!'
    WRITE(IU06,*) '!!! ON PE ',IRANK
    WRITE(IU06,*) '!!! ********************************* !!'
    CALL ABORT1
  ENDIF
ENDDO

! FIND THE LOGICAL FLAGS THAT WILL LIMIT THE EXTEND OF THE CALCULATION IN PROPAGS2

              DO IC=1,2
                DO ICL=1,2
                  DO K=1,NANG
                    DO M=1,NFRE_RED
                      LLWLATN(K,M,IC,ICL)=.FALSE.
                      DO IJ=IJS,IJL
                        IF(WLATN(IJ,K,M,IC,ICL).GT.0.0_JWRB) THEN
                          LLWLATN(K,M,IC,ICL)=.TRUE.
                          EXIT
                        ENDIF
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO

              DO IC=1,2
                DO M=1,NFRE_RED
                  DO K=1,NANG
                    LLWLONN(K,M,IC)=.FALSE.
                    DO IJ=IJS,IJL
                      IF(WLONN(IJ,K,M,IC).GT.0.0_JWRB) THEN
                        LLWLONN(K,M,IC)=.TRUE.
                        EXIT
                      ENDIF
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO

              DO ICL=1,2
                DO ICR=1,4
                  DO M=1,NFRE_RED
                    DO K=1,NANG
                      LLWCORN(K,M,ICR,ICL)=.FALSE.
                      DO IJ=IJS,IJL
                        IF(WCORN(IJ,K,M,ICR,ICL).GT.0.0_JWRB) THEN
                          LLWCORN(K,M,ICR,ICL)=.TRUE.
                          EXIT
                        ENDIF
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO

              DO IC=-1,1
                DO M=1,NFRE_RED
                  DO K=1,NANG
                   LLWKPMN(K,M,IC)=.FALSE.
                    DO IJ=IJS,IJL
                      IF(WKPMN(IJ,K,M,IC).GT.0.0_JWRB) THEN
                        LLWKPMN(K,M,IC)=.TRUE.
                        EXIT
                      ENDIF
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO

              IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN
                DO IC=-1,1
                  DO M=1,NFRE_RED
                    DO K=1,NANG
                     LLWMPMN(K,M,IC)=.FALSE.
                      DO IJ=IJS,IJL
                        IF(WMPMN(IJ,K,M,IC).GT.0.0_JWRB) THEN
                          LLWMPMN(K,M,IC)=.TRUE.
                          EXIT
                        ENDIF
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF


IF(ALLOCATED(THDD)) DEALLOCATE(THDD)
IF(ALLOCATED(THDC)) DEALLOCATE(THDC)
IF(ALLOCATED(SDOT)) DEALLOCATE(SDOT)

IF (LHOOK) CALL DR_HOOK('CTUWUPDT',1,ZHOOK_HANDLE)

END SUBROUTINE CTUWUPDT
