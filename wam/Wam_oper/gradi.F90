      SUBROUTINE GRADI (MIJS, MIJL, IREFRA, DDPHI, DDLAM, DUPHI,       &
     &                  DULAM, DVPHI, DVLAM)

! ----------------------------------------------------------------------

!**** *GRADI* - CALCULATES DEPTH AND CURRENT VELOCITY GRADIENTS.

!     K.P. HUBBERT              AUGUST   1988
!     H. GUNTHER    ECMWF/GKSS  DECEMBER 1990  MODIFIED FOR CYCLE_4.

!*    PURPOSE.
!     --------

!       CALCULATES DEPTH AND CURRENT VELOCITY GRADIENTS OF A BLOCK.

!**   INTERFACE.
!     ----------

!       *CALL* *GRADI (MIJS, MIJL,IREFRA, DDPHI, DDLAM, DUPHI,
!                      DULAM, DVPHI, DVLAM)*
!          *MIJS*   - STARTING INDEX
!          *MIJL*   - ENDING INDEX
!          *IREFRA* - REFRACTION OPTION.
!          *DDPHI*  - LATITUDE DEPTH GRADIENT.
!          *DDLAM*  - LONGITUDE DEPTH GRADIENT.
!          *DUPHI*  - LATITUDE  U-COMPONENT GRADIENT.
!          *DULAM*  - LONGITUDE U-COMPONENT GRADIENT.
!          *DVPHI*  - LATITUDE  V-COMPONENT GRADIENT.
!          *DVLAM*  - LONGITUDE V-COMPONENT GRADIENT.

!     METHOD.
!     ------

!       CENTRAL DIFFERENCING FOR DEPTH AND CURRENT VELOCITY GRADIENTS.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCURR  , ONLY : U        ,V        ,CURRENT_GRADIENT_MAX
!                       !!! debile
      USE YOWGRID  , ONLY : DELPHI   ,DELLAM   ,COSPH
      USE YOWMAP   , ONLY : KXLT
      USE YOWMPP   , ONLY : NSUP 
      USE YOWSHAL  , ONLY : DEPTH
!                      !!!!!!!!!debile
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,WLAT

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: MIJS
      INTEGER(KIND=JWIM), INTENT(IN) :: MIJL
      INTEGER(KIND=JWIM), INTENT(IN) :: IREFRA
      REAL(KIND=JWRB),DIMENSION(MIJS:MIJL), INTENT(OUT) :: DDPHI
      REAL(KIND=JWRB),DIMENSION(MIJS:MIJL), INTENT(OUT) :: DDLAM
      REAL(KIND=JWRB),DIMENSION(MIJS:MIJL), INTENT(OUT) :: DUPHI
      REAL(KIND=JWRB),DIMENSION(MIJS:MIJL), INTENT(OUT) :: DULAM
      REAL(KIND=JWRB),DIMENSION(MIJS:MIJL), INTENT(OUT) :: DVPHI
      REAL(KIND=JWRB),DIMENSION(MIJS:MIJL), INTENT(OUT) :: DVLAM


      INTEGER(KIND=JWIM) :: NLAND, IJ, IPP, IPM, IPP2, IPM2, ILP, ILM, KX
      REAL(KIND=JWRB) :: DPTP, DPTM, UP, UM, VP, VM 
      REAL(KIND=JWRB) :: CGMAX
      REAL(KIND=JWRB) :: ONEO2DELPHI
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('GRADI',0,ZHOOK_HANDLE)

!*    1. INITIALISE.
!        -----------

      NLAND=NSUP+1
      ONEO2DELPHI = 0.5_JWRB/DELPHI

! ----------------------------------------------------------------------

!*    2. CALCULATE DEPTH GRADIENTS.
!        --------------------------

      IF (IREFRA == 1 .OR. IREFRA == 3) THEN
        DO IJ=MIJS,MIJL
          IPP = KLAT(IJ,2,1)
          IPM = KLAT(IJ,1,1)
          IPP2 = KLAT(IJ,2,2)
          IPM2 = KLAT(IJ,1,2)
          IF (IPP /= NLAND  .AND. IPM /= NLAND .AND.                    &
     &        IPP2 /= NLAND .AND. IPM2 /= NLAND ) THEN
            DPTP=WLAT(IJ,2)*DEPTH(IPP)+(1.0_JWRB-WLAT(IJ,2))*DEPTH(IPP2)
            DPTM=WLAT(IJ,1)*DEPTH(IPM)+(1.0_JWRB-WLAT(IJ,1))*DEPTH(IPM2)
            DDPHI(IJ) = (DPTP-DPTM)*ONEO2DELPHI
          ELSEIF (IPP /= NLAND .AND. IPM /= NLAND) THEN
            DPTP=DEPTH(IPP)
            DPTM=DEPTH(IPM)
            DDPHI(IJ) = (DPTP-DPTM)*ONEO2DELPHI
          ELSEIF (IPP2 /= NLAND .AND. IPM2 /= NLAND) THEN
            DPTP=DEPTH(IPP2)
            DPTM=DEPTH(IPM2)
            DDPHI(IJ) = (DPTP-DPTM)*ONEO2DELPHI
          ELSE
            DDPHI(IJ) = 0.0_JWRB
          ENDIF
          ILP = KLON(IJ,2)
          ILM = KLON(IJ,1)
          KX  = KXLT(IJ)
          IF (ILP /= NLAND .AND. ILM /= NLAND) THEN
            DDLAM(IJ)=(DEPTH(ILP)-DEPTH(ILM))/(2._JWRB*DELLAM(KX))
          ELSE
            DDLAM(IJ) = 0.0_JWRB 
          ENDIF
        ENDDO
      ELSE
        DO IJ=MIJS,MIJL
          DDPHI(IJ) = 0.0_JWRB
          DDLAM(IJ) = 0.0_JWRB
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    3. CALCULATE CURRENT VELOCITY GRADIENTS.
!        -------------------------------------

      IF (IREFRA == 2 .OR. IREFRA == 3) THEN
        DO IJ=MIJS,MIJL
          IPP = KLAT(IJ,2,1)
!         exact 0 means that the current field was not defined, hence
!         no gradient should be extrapolated
          IF (U(IPP) == 0.0_JWRB .AND. V(IPP) == 0.0_JWRB) IPP = NLAND
          IPM = KLAT(IJ,1,1)
          IF (U(IPM) == 0.0_JWRB .AND. V(IPM) == 0.0_JWRB) IPM = NLAND
          IPP2 = KLAT(IJ,2,2)
          IF (U(IPP2) == 0.0_JWRB .AND. V(IPP2) == 0.0_JWRB) IPP2 = NLAND
          IPM2 = KLAT(IJ,1,2)
          IF (U(IPM2) == 0.0_JWRB .AND. V(IPM2) == 0.0_JWRB) IPM2 = NLAND

          IF (IPP /= NLAND .AND. IPM /= NLAND .AND.                     &
     &        IPP2 /= NLAND .AND. IPM2 /= NLAND) THEN
            UP = WLAT(IJ,2)*U(IPP)+(1.0_JWRB-WLAT(IJ,2))*U(IPP2)
            VP = WLAT(IJ,2)*V(IPP)+(1.0_JWRB-WLAT(IJ,2))*V(IPP2)
            UM = WLAT(IJ,1)*U(IPM)+(1.0_JWRB-WLAT(IJ,1))*U(IPM2)
            VM = WLAT(IJ,1)*V(IPM)+(1.0_JWRB-WLAT(IJ,1))*V(IPM2)
            DUPHI(IJ) = (UP-UM)*ONEO2DELPHI
            DVPHI(IJ) = (VP-VM)*ONEO2DELPHI
          ELSEIF (IPP /= NLAND .AND. IPM /= NLAND) THEN
            UP = U(IPP)
            VP = V(IPP)
            UM = U(IPM)
            VM = V(IPM)
            DUPHI(IJ) = (UP-UM)*ONEO2DELPHI
            DVPHI(IJ) = (VP-VM)*ONEO2DELPHI
          ELSE
            DUPHI(IJ) = 0.0_JWRB
            DVPHI(IJ) = 0.0_JWRB
          ENDIF

          ILP = KLON(IJ,2)
!         exact 0 means that the current field was not defined, hence
!         no gradient should be extrapolated
          IF (U(ILP) == 0.0_JWRB .AND. V(ILP) == 0.0_JWRB) ILP=NLAND
          ILM = KLON(IJ,1)
          IF (U(ILM) == 0.0_JWRB .AND. V(ILM) == 0.0_JWRB) ILM=NLAND
          KX  = KXLT(IJ)
          IF (ILP /= NLAND .AND. ILM /= NLAND) THEN
            DULAM(IJ) = (U(ILP)-U(ILM))/(2.0_JWRB*DELLAM(KX))
            DVLAM(IJ) = (V(ILP)-V(ILM))/(2.0_JWRB*DELLAM(KX))
          ELSE
            DULAM(IJ) = 0.0_JWRB
            DVLAM(IJ) = 0.0_JWRB
          ENDIF
        ENDDO

        DO IJ=MIJS,MIJL
          KX  = KXLT(IJ)
          CGMAX = CURRENT_GRADIENT_MAX*COSPH(KX)
          DUPHI(IJ) = SIGN(MIN(ABS(DUPHI(IJ)),CGMAX),DUPHI(IJ))
          DVPHI(IJ) = SIGN(MIN(ABS(DVPHI(IJ)),CGMAX),DVPHI(IJ))
          DULAM(IJ) = SIGN(MIN(ABS(DULAM(IJ)),CGMAX),DULAM(IJ))
          DVLAM(IJ) = SIGN(MIN(ABS(DVLAM(IJ)),CGMAX),DVLAM(IJ))
        ENDDO

      ELSE
        DO IJ=MIJS,MIJL
          DUPHI(IJ) = 0.0_JWRB
          DVPHI(IJ) = 0.0_JWRB
          DULAM(IJ) = 0.0_JWRB
          DVLAM(IJ) = 0.0_JWRB
        ENDDO
      ENDIF

      IF (LHOOK) CALL DR_HOOK('GRADI',1,ZHOOK_HANDLE)

      END SUBROUTINE GRADI
