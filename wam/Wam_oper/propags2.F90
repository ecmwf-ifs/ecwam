      SUBROUTINE PROPAGS2 (F1, F3, IJS, IJL, MIJS, MIJL)

! ----------------------------------------------------------------------

!**** *PROPAGS2* - COMPUTATION OF A PROPAGATION TIME STEP.

! ----------------------------------------------------------------------

!**** *PROPAGS2* -  ADVECTION USING THE CORNER TRANSPORT SCHEME IN SPACE


!*    PURPOSE.
!     --------

!       COMPUTATION OF A PROPAGATION TIME STEP.

!**   INTERFACE.
!     ----------

!       *CALL* *PROPAGS2(F1, F3, IJS, IJL, MIJS, MIJL)*
!          *F1*   - BLOCK SPECTRUM AT TIME T.
!          *F3(IJS:IJL,:,:)*   - SPECTRUM AT TIME T+DELT.
!          *MIJS* - INDEX OF FIRST POINT
!          *MIJL* - INDEX OF LAST POINT


!     METHOD.
!     -------

!!  need text !!!!!!!!


!     EXTERNALS.
!     ----------


!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : COSTH    ,SINTH
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSTAT  , ONLY : ICASE    ,ISHALLO  ,IREFRA
      USE YOWTEST  , ONLY : IU06
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,KCOR      ,             &
     &            WLATN    ,WLONN    ,WCORN    ,WKPMN    ,WMPMN     ,   &
     &            LLWLATN  ,LLWLONN  ,LLWCORN  ,LLWKPMN  ,LLWMPMN   ,   &
     &            SUMWN    ,                                            &
     &            JXO      ,JYO      ,KCR      ,KPM      ,MPM
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: MIJS, MIJL, IJS, IJL

      REAL(KIND=JWRB),DIMENSION(NINF-1:NSUP,NANG,NFRE), INTENT(IN) :: F1
      REAL(KIND=JWRB),DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT)  :: F3

      INTEGER(KIND=JWIM) :: K, M, IJ
      INTEGER(KIND=JWIM) :: IC, ICR, ICL 
      INTEGER(KIND=JWIM) :: KP1, KM1, MM1, MP1, KNS, KEW
      INTEGER(KIND=JWIM) :: JJK, JJY, JJX

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB),DIMENSION(MIJS:MIJL) :: FJ1, FJ2, FJ3, FJ4, FJ5

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('PROPAGS2',0,ZHOOK_HANDLE)

!*    SPHERICAL OR CARTESIAN GRID?
!     ----------------------------
      IF (ICASE.EQ.1) THEN

!*      SPHERICAL GRID.
!       ---------------

        IF (IREFRA.NE.2 .AND. IREFRA.NE.3 ) THEN
!*      WITHOUT DEPTH OR/AND CURRENT REFRACTION.
!       ----------------------------------------

          DO K=1,NANG
            JJX=JXO(K,1)
            JJY=JYO(K,1)
            JJY=JYO(K,1)
            JJK=KCR(K,1)
            DO M=1,NFRE
              DO IJ=MIJS,MIJL
                FJ1(IJ)= F1(KLON(IJ,JJX)  ,K  ,M)
                FJ2(IJ)= F1(KLAT(IJ,JJY,1),K  ,M)
                FJ3(IJ)= F1(KLAT(IJ,JJY,2),K  ,M)
                FJ4(IJ)= F1(KCOR(IJ,JJK,1),K  ,M)
                FJ5(IJ)= F1(KCOR(IJ,JJK,2),K  ,M)
              ENDDO
!JFH Loop split to enhance vectorisation
!DIR$ IVDEP
!DIR$ PREFERVECTOR
              DO IJ=MIJS,MIJL
                F3(IJ,K,M) =                                            &
     &                (1.0_JWRB-SUMWN(IJ,K,M))* F1(IJ           ,K  ,M) &
!    &         + WLONN(IJ,K,M,JXO(K,1)) * F1(KLON(IJ,JXO(K,1))  ,K  ,M) &
!    &         +WLATN(IJ,K,M,JYO(K,1),1)* F1(KLAT(IJ,JYO(K,1),1),K  ,M) &
!    &         +WLATN(IJ,K,M,JYO(K,1),2)* F1(KLAT(IJ,JYO(K,1),2),K  ,M) &
!    &         +       WCORN(IJ,K,M,1,1)* F1(KCOR(IJ,KCR(K,1),1),K  ,M) &
!    &         +       WCORN(IJ,K,M,1,2)* F1(KCOR(IJ,KCR(K,1),2),K  ,M) &
!    &         + WLONN(IJ,K,M,JJX) * F1(KLON(IJ,JJX)  ,K  ,M)           &
!    &         +WLATN(IJ,K,M,JJY,1)* F1(KLAT(IJ,JJY,1),K  ,M)           &
!    &         +WLATN(IJ,K,M,JJY,2)* F1(KLAT(IJ,JJY,2),K  ,M)           &
!    &         +       WCORN(IJ,K,M,1,1)* F1(KCOR(IJ,JJK,1),K  ,M)      &
!    &         +       WCORN(IJ,K,M,1,2)* F1(KCOR(IJ,JJK,2),K  ,M)      &
     &         + WLONN(IJ,K,M,JJX) * FJ1(IJ)                            &
     &         +WLATN(IJ,K,M,JJY,1)* FJ2(IJ)                            &
     &         +WLATN(IJ,K,M,JJY,2)* FJ3(IJ)                            &
     &         +       WCORN(IJ,K,M,1,1)* FJ4(IJ)                       &
     &         +       WCORN(IJ,K,M,1,2)* FJ5(IJ)
              ENDDO

              DO IC=-1,1,2
                IF(LLWKPMN(K,M,IC)) THEN
                  DO IJ=MIJS,MIJL
                    F3(IJ,K,M) = F3(IJ,K,M)                             &
     &         +      WKPMN(IJ,K,M,IC)* F1(IJ,KPM(K,IC),M)
                  ENDDO
                ENDIF
              ENDDO

            ENDDO
          ENDDO

        ELSE
!*      DEPTH AND CURRENT REFRACTION.
!       -----------------------------

          DO M=1,NFRE
            DO K=1,NANG

              DO IJ=MIJS,MIJL
                F3(IJ,K,M) = (1.0_JWRB-SUMWN(IJ,K,M))* F1(IJ,K,M)
              ENDDO
              
              DO IC=1,2
                IF(LLWLONN(K,M,IC)) THEN
                  DO IJ=MIJS,MIJL
                    F3(IJ,K,M) = F3(IJ,K,M)                             &
     &           +      WLONN(IJ,K,M,IC)*F1(KLON(IJ,IC) ,K,M)
                  ENDDO
                ENDIF
              ENDDO

              DO ICL=1,2
                DO IC=1,2
                  IF(LLWLATN(K,M,IC,ICL)) THEN
                    DO IJ=MIJS,MIJL
                      F3(IJ,K,M) = F3(IJ,K,M)                           &
     &         +        WLATN(IJ,K,M,IC,ICL)*F1(KLAT(IJ,IC,ICL) ,K,M)
                    ENDDO
                 ENDIF
                ENDDO
              ENDDO

              DO ICL=1,2
                DO ICR=1,4
                  IF(LLWCORN(K,M,ICR,ICL)) THEN
                    DO IJ=MIJS,MIJL
                      F3(IJ,K,M) = F3(IJ,K,M)                           &
     &         +   WCORN(IJ,K,M,ICR,ICL)*F1(KCOR(IJ,KCR(K,ICR),ICL),K,M)
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO

              DO IC=-1,1,2
                IF(LLWKPMN(K,M,IC)) THEN
                  DO IJ=MIJS,MIJL
                    F3(IJ,K,M) = F3(IJ,K,M)                             &
     &         +      WKPMN(IJ,K,M,IC)* F1(IJ,KPM(K,IC),M)
                  ENDDO
                ENDIF
              ENDDO


              DO IC=-1,1,2
                IF(LLWMPMN(K,M,IC)) THEN
                  DO IJ=MIJS,MIJL
                    F3(IJ,K,M) = F3(IJ,K,M)                             &
     &         +      WMPMN(IJ,K,M,IC)* F1(IJ,K  ,MPM(M,IC))
                  ENDDO
                ENDIF
              ENDDO


            ENDDO

          ENDDO

        ENDIF

      ELSE
!*    CARTESIAN GRID.
!     ---------------
        IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3 ) THEN
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
