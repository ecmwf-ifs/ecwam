      SUBROUTINE SDISSIP (F, FL, SL, IJS, IJL, &
     &                    EMEAN, F1MEAN, XKMEAN, &
     &                    USNEW, THWNEW, ROAIRN)
! ----------------------------------------------------------------------

!**** *SDISSIP* - COMPUTATION OF DEEP WATER DISSIPATION SOURCE FUNCTION.


!*    PURPOSE.
!     --------
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SDISSIP (F, FL, IJS, IJL, SL,*
!                        EMEAN, F1MEAN, XKMEAN,*
!                        USNEW, THWNEW,ROAIRN)*
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *EMEAN* - MEAN ENERGY DENSITY 
!          *F1MEAN* - MEAN FREQUENCY BASED ON 1st MOMENT.
!          *XKMEAN* - MEAN WAVE NUMBER BASED ON 1st MOMENT.
!          *USNEW*  - NEW FRICTION VELOCITY IN M/S.
!          *ROAIRN* - AIR DENSITY IN KG/M3
!          *THWNEW* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC.

!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       ARDHUIN et AL. JPO DOI:10.1175/20110JPO4324.1


! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSTAT  , ONLY : IPHYS
      USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL

      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: EMEAN, F1MEAN, XKMEAN
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: USNEW, THWNEW, ROAIRN 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: FL, SL

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SDISSIP',0,ZHOOK_HANDLE)
#endif

      SELECT CASE (IPHYS)
      CASE(0)
         CALL SDISSIP_JAN (F ,FL, SL, IJS, IJL, &
     &                     EMEAN, F1MEAN, XKMEAN)

      CASE(1) 
         CALL SDISSIP_ARD (F ,FL, SL, IJS, IJL, &
     &                     USNEW, THWNEW, ROAIRN)
      END SELECT 


#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SDISSIP',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE SDISSIP
