      SUBROUTINE SBOTTOM (F, FL, IJS, IJL, IG, SL)

!SHALLOW
! ----------------------------------------------------------------------

!**** *SBOTTOM* - COMPUTATION OF BOTTOM FRICTION.

!     G.J.KOMEN AND Q.D.GAO
!     OPTIMIZED BY L.F. ZAMBRESKY
!     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL

!*    PURPOSE.
!     --------

!       COMPUTATION OF BOTTOM FRICTION DISSIPATION

!**   INTERFACE.
!     ----------

!       *CALL* *SBOTTOM (F, FL, IJS, IJL, IG, SL)*
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *IG*  - BLOCK NUMBER
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY

!     METHOD.
!     -------

!       SEE REFERENCES.

!     REFERENCES.
!     -----------

!       HASSELMANN ET AL, D. HYDR. Z SUPPL A12(1973) (JONSWAP)
!       BOUWS AND KOMEN, JPO 13(1983)1653-1658

! ----------------------------------------------------------------------

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G
      USE YOWSHAL  , ONLY : DEPTH    ,TFAK     ,INDEP
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      PARAMETER (CONST = -2.0*0.038/G)
      REAL,DIMENSION(IJS:IJL,NANG,NFRE) :: F,FL,SL
      DIMENSION SBO(IJS:IJL)

      REAL ZHOOK_HANDLE

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('SBOTTOM',0,ZHOOK_HANDLE)

      DO M=1,NFRE
        DO IJ=IJS,IJL
          IF(DEPTH(IJ,IG).LT.999) THEN
            ARG = 2.* DEPTH(IJ,IG)*TFAK(INDEP(IJ),M)
            ARG = MIN(ARG,50.)
            SBO(IJ) = CONST*TFAK(INDEP(IJ),M)/SINH(ARG)
          ENDIF
        ENDDO

        DO K=1,NANG
          DO IJ=IJS,IJL
            IF(DEPTH(IJ,IG).LT.999) THEN
              SL(IJ,K,M) = SL(IJ,K,M)+SBO(IJ)*F(IJ,K,M)
              FL(IJ,K,M) = FL(IJ,K,M)+SBO(IJ)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('SBOTTOM',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE SBOTTOM
