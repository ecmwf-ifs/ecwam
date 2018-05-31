      SUBROUTINE SDIWBK (IJS, IJL, F, FL, SL,
     &                   DEPTH, EMAXDPT, EMEAN, F1MEAN)

! ----------------------------------------------------------------------

!**** *SDIWBK* - COMPUTATION OF BOTTOM-INDUCED WAVE BREAKING DISSIPATION


!*    PURPOSE.
!     --------
!       COMPUTE BOTTOM-INDUCED DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SDIWBK (IJS, IJL, F, FL, SL,
!    &                   DEPTH, EMAXDPT, EMEAN, F1MEAN)
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY
!          *DEPTH* WATER DEPTH
!          *EMAXDPT* MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH
!          *EMEAN* - MEAN ENERGY DENSITY
!          *F1MEAN* - MEAN FREQUENCY BASED ON 1st MOMENT.

!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSTAT  , ONLY : LBIWBK
      USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL

      INTEGER :: IJ, K, M, IC

      REAL,DIMENSION(IJS:IJL), INTENT(IN) :: DEPTH,EMAXDPT,EMEAN,F1MEAN
      REAL,DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F
      REAL,DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: FL, SL

      REAL :: ZHOOK_HANDLE
      REAL :: ALPH, ARG, Q, Q_OLD, REL_ERR, EXPQ
      REAL,DIMENSION(IJS:IJL) :: SDS

      REAL, PARAMETER :: ALPH_B_J = 1.0
      REAL, PARAMETER :: COEF_B_J=2*ALPH_B_J
      REAL, PARAMETER :: DEPTHTRS=50. 

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SDIWBK',0,ZHOOK_HANDLE)
#endif

!*    1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE
!*       FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.
!        --------------------------------------------------------------

      IF (LBIWBK) THEN
!       (FOLLOWING BATTJES-JANSSEN AND BEJI)
        DO IJ=IJS,IJL
           IF(DEPTH(IJ).LT.DEPTHTRS) THEN
             ALPH = 2.*EMAXDPT(IJ)/EMEAN(IJ)
             ARG  = MIN(ALPH,50.)
             Q_OLD = EXP(-ARG)
!            USE NEWTON-RAPHSON METHOD
             DO IC=1,15
               EXPQ = EXP(-ARG*(1.-Q_OLD))
               Q = Q_OLD - (EXPQ-Q_OLD)/(ARG*EXPQ-1.)
               REL_ERR=ABS(Q-Q_OLD)/Q_OLD
               IF(REL_ERR.LT.0.00001) EXIT
               Q_OLD = Q
             ENDDO
             Q=MIN(Q,1.)
             SDS(IJ) = COEF_B_J*ALPH*Q*F1MEAN(IJ)
           ENDIF
        ENDDO 
      
        DO M=1,NFRE
           DO K=1,NANG
              DO IJ=IJS,IJL
                IF(DEPTH(IJ).LT.DEPTHTRS) THEN
                  SL(IJ,K,M) = SL(IJ,K,M)-SDS(IJ)*F(IJ,K,M)
                  FL(IJ,K,M) = FL(IJ,K,M)-SDS(IJ)
                ENDIF
              ENDDO
           ENDDO
        ENDDO
     
      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SDIWBK',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE SDIWBK
