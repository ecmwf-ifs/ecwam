! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE PROPAGS2 (F1, F3, NINF, NSUP, KIJS, KIJL, NANG, ND3SF1, ND3EF1, ND3S, ND3E)

! ----------------------------------------------------------------------

!**** *PROPAGS2* -  ADVECTION USING THE CORNER TRANSPORT SCHEME IN SPACE

!*    PURPOSE.
!     --------

!       COMPUTATION OF A PROPAGATION TIME STEP.

!**   INTERFACE.
!     ----------

!       *CALL* *PROPAGS2(F1, F3, NINF, NSUP, KIJS, KIJL, NANG, ND3SF1, ND3EF1, ND3S, ND3E)*
!          *F1*          - SPECTRUM AT TIME T (with exchange halo).
!          *F3*          - SPECTRUM AT TIME T+DELT
!          *NINF:NSUP+1* - 1st DIMENSION OF F1 and F3
!          *KIJS*        - ACTIVE INDEX OF FIRST POINT
!          *KIJL*        - ACTIVE INDEX OF LAST POINT
!          *NANG*        - NUMBER OF DIRECTIONS
!          *ND3SF1*      - LOWER 3rd DIMENSION OF F1 
!          *ND3EF1*      - UPPER 3d DIMENSION OF F1 
!          *ND3S*        - FREQUENCY INDEX SOLVED BY THIS CALL ND3S:ND3E
!          *ND3E*        - FREQUENCY INDEX SOLVED BY THIS CALL ND3S:ND3E

!     METHOD.
!     -------

!!  need text !!!!!!!!

!     EXTERNALS.
!     ----------


!     REFERENCE.
!     ----------

!      See IFS Documentation, part VII 

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : COSTH    ,SINTH
      USE YOWSTAT  , ONLY : ICASE    ,IREFRA
      USE YOWTEST  , ONLY : IU06
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,KCOR      ,             &
     &            WLATN    ,WLONN    ,WCORN    ,WKPMN    ,WMPMN     ,   &
     &            LLWLATN  ,LLWLONN  ,LLWCORN  ,LLWKPMN  ,LLWMPMN   ,   &
     &            SUMWN    ,                                            &
     &            JXO      ,JYO      ,KCR      ,KPM      ,MPM

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE YOWABORT , ONLY : WAM_ABORT

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"

      REAL(KIND=JWRB),DIMENSION(NINF:NSUP+1, NANG, ND3SF1:ND3EF1), INTENT(IN) :: F1
      REAL(KIND=JWRB),DIMENSION(NINF:NSUP+1, NANG, ND3S:ND3E), INTENT(OUT) :: F3
      INTEGER(KIND=JWIM), INTENT(IN) :: NINF, NSUP, KIJS, KIJL
      INTEGER(KIND=JWIM), INTENT(IN) :: NANG, ND3SF1, ND3EF1, ND3S, ND3E


      INTEGER(KIND=JWIM) :: K, M, IJ
      INTEGER(KIND=JWIM) :: IC, ICR, ICL 
      INTEGER(KIND=JWIM) :: KP1, KM1, MM1, MP1, KNS, KEW

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PROPAGS2',0,ZHOOK_HANDLE)

!*    SPHERICAL OR CARTESIAN GRID?
!     ----------------------------
      IF (ICASE == 1) THEN

!*      SPHERICAL GRID.
!       ---------------

        IF (IREFRA /= 2 .AND. IREFRA /= 3 ) THEN
!*      WITHOUT DEPTH OR/AND CURRENT REFRACTION.
!       ----------------------------------------

          !$acc parallel loop independent collapse(3) &
          !$acc & present(F1,F3,KLON,KLAT,KCOR,SUMWN,WLONN,WLATN,WCORN,JXO,JYO,KCR,WKPMN,KPM)
          DO K = 1, NANG
            DO M = ND3S, ND3E

!DIR$ IVDEP
!DIR$ PREFERVECTOR
              DO IJ = KIJS, KIJL
                F3(IJ,K,M) =                                            &
     &                (1.0_JWRB-SUMWN(IJ,K,M))* F1(IJ           ,K  ,M) &

     &         + WLONN(IJ,K,M,JXO(K,1))   * F1(KLON(IJ,JXO(K,1))  ,K  ,M) &
     &         + WLATN(IJ,K,M,JYO(K,1),1) * F1(KLAT(IJ,JYO(K,1),1),K  ,M) &
     &         + WLATN(IJ,K,M,JYO(K,1),2) * F1(KLAT(IJ,JYO(K,1),2),K  ,M) &
     &         + WCORN(IJ,K,M,1,1)        * F1(KCOR(IJ,KCR(K,1),1),K  ,M) &
     &         + WCORN(IJ,K,M,1,2)        * F1(KCOR(IJ,KCR(K,1),2),K  ,M) &
     &         + WKPMN(IJ,K,M,-1)         * F1(IJ,KPM(K,-1),M)            &
     &         + WKPMN(IJ,K,M, 1)         * F1(IJ,KPM(K, 1),M)
              ENDDO

            ENDDO
          ENDDO
          !$acc end parallel loop 

        ELSE
!*      DEPTH AND CURRENT REFRACTION.
!       -----------------------------

          !$acc parallel loop independent collapse(3) &
          !$acc & present(F1,F3,SUMWN,WLONN,KLON,WLATN,KLAT,WCORN,KCOR, &
          !$acc &         WKPMN,KPM,WMPMN,MPM)
          DO M = ND3S, ND3E
            DO K = 1, NANG

              !$loki loop-fusion
              DO IJ = KIJS, KIJL
                F3(IJ,K,M) = (1.0_JWRB-SUMWN(IJ,K,M))* F1(IJ,K,M)
              ENDDO

              !$loki loop-unroll
              DO IC=1,2
                IF (LLWLONN(K,M,IC)) THEN
                  !$loki loop-fusion
                  DO IJ = KIJS, KIJL
                    F3(IJ,K,M) = F3(IJ,K,M) + WLONN(IJ,K,M,IC)*F1(KLON(IJ,IC),K,M)
                  ENDDO
                ENDIF
              ENDDO

              !$loki loop-unroll
              DO ICL=1,2
                !$loki loop-unroll
                DO IC=1,2
                  IF (LLWLATN(K,M,IC,ICL)) THEN
                    !$loki loop-fusion
                    DO IJ = KIJS, KIJL
                      F3(IJ,K,M) = F3(IJ,K,M) + WLATN(IJ,K,M,IC,ICL)*F1(KLAT(IJ,IC,ICL),K,M)
                    ENDDO
                  ENDIF
                ENDDO

                !$loki loop-unroll
                DO ICR=1,4
                  IF (LLWCORN(K,M,ICR,ICL)) THEN
                    !$loki loop-fusion
                    DO IJ = KIJS, KIJL
                      F3(IJ,K,M) = F3(IJ,K,M) + WCORN(IJ,K,M,ICR,ICL)*F1(KCOR(IJ,KCR(K,ICR),ICL),K,M)
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO

              !$loki loop-unroll
              DO IC=-1,1,2

                IF (LLWKPMN(K,M,IC)) THEN
                  !$loki loop-fusion
                  DO IJ = KIJS, KIJL
                    F3(IJ,K,M) = F3(IJ,K,M) + WKPMN(IJ,K,M,IC)* F1(IJ,KPM(K,IC),M)
                  ENDDO
                ENDIF

                IF (LLWMPMN(K,M,IC)) THEN
                  !$loki loop-fusion
                  DO IJ = KIJS, KIJL
                    F3(IJ,K,M) = F3(IJ,K,M) + WMPMN(IJ,K,M,IC)* F1(IJ,K,MPM(M,IC))
                  ENDDO
                ENDIF

              ENDDO

            ENDDO
          ENDDO
          !$acc end parallel loop

        ENDIF

      ELSE
!*    CARTESIAN GRID.
!     ---------------
        IF (IREFRA == 2 .OR. IREFRA == 3 ) THEN
!*      WITHOUT DEPTH OR/AND CURRENT REFRACTION.
!       ----------------------------------------
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '* PROPAGS2:                              *'
          WRITE (IU06,*) '* CORNER TRANSPORT SCHEME NOT YET READY  *' 
          WRITE (IU06,*) '* FOR  CARTESIAN GRID !                  *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.      *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          CALL ABORT1
        ELSE
!*      DEPTH AND CURRENT REFRACTION.
!       ----------------------------
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '* PROPAGS2:                              *'
          WRITE (IU06,*) '* CORNER TRANSPORT SCHEME NOT YET READY  *' 
          WRITE (IU06,*) '* FOR  CARTESIAN GRID !                  *'
          WRITE (IU06,*) '* FOR DEPTH OR/AND CURRENT REFRACTION !  *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.      *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          CALL ABORT1
        ENDIF

      ENDIF

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PROPAGS2',1,ZHOOK_HANDLE)

END SUBROUTINE PROPAGS2 
