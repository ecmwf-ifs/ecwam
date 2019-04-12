      SUBROUTINE GRADI (IG, MIJS, MIJL, IREFRA, DDPHI, DDLAM, DUPHI,    &
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

!       *CALL* *GRADI (IG, MIJS, MIJL,IREFRA, DDPHI, DDLAM, DUPHI,
!                      DULAM, DVPHI, DVLAM)*
!          *IG*     - BLOCK NUMBER.
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
      USE YOWGRID  , ONLY : DELPHI   ,DELLAM
      USE YOWMAP   , ONLY : KXLT
      USE YOWMPP   , ONLY : NINF
      USE YOWSHAL  , ONLY : DEPTH
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,WLAT

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IG
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
      REAL(KIND=JWRB) :: ONEO2DELPHI
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('GRADI',0,ZHOOK_HANDLE)

!*    1. INITIALISE.
!        -----------

      NLAND=NINF-1
      ONEO2DELPHI = 0.5_JWRB/DELPHI

! ----------------------------------------------------------------------

!*    2. CALCULATE DEPTH GRADIENTS.
!        --------------------------

      IF (IREFRA.EQ.1 .OR. IREFRA.EQ.3) THEN
        DO IJ=MIJS,MIJL
          IPP = KLAT(IJ,2,1)
          IPM = KLAT(IJ,1,1)
          IPP2 = KLAT(IJ,2,2)
          IPM2 = KLAT(IJ,1,2)
          IF (IPP.NE.NLAND  .AND. IPM.NE.NLAND .AND.                    &
     &        IPP2.NE.NLAND .AND. IPM2.NE.NLAND      ) THEN
            DPTP=WLAT(IJ,2)*DEPTH(IPP,IG)+(1.0_JWRB-WLAT(IJ,2))*DEPTH(IPP2,IG)
            DPTM=WLAT(IJ,1)*DEPTH(IPM,IG)+(1.0_JWRB-WLAT(IJ,1))*DEPTH(IPM2,IG)
            DDPHI(IJ) = (DPTP-DPTM)*ONEO2DELPHI
          ELSEIF (IPP.NE.NLAND .AND. IPM.NE.NLAND) THEN
            DPTP=DEPTH(IPP,IG)
            DPTM=DEPTH(IPM,IG)
            DDPHI(IJ) = (DPTP-DPTM)*ONEO2DELPHI
          ELSEIF (IPP2.NE.NLAND .AND. IPM2.NE.NLAND) THEN
            DPTP=DEPTH(IPP2,IG)
            DPTM=DEPTH(IPM2,IG)
            DDPHI(IJ) = (DPTP-DPTM)*ONEO2DELPHI
          ELSE
            DDPHI(IJ) = 0.0_JWRB
          ENDIF
          ILP = KLON(IJ,2)
          ILM = KLON(IJ,1)
          KX  = KXLT(IJ,IG)
          IF (ILP.NE.NLAND .AND. ILM.NE.NLAND) THEN
            DDLAM(IJ)=(DEPTH(ILP,IG)-DEPTH(ILM,IG))/(2._JWRB*DELLAM(KX))
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

      IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN
        DO IJ=MIJS,MIJL
          IPP = KLAT(IJ,2,1)
!         exact 0 means that the current field was not defined, hence
!         no gradient should be extrapolated
          IF (U(IPP,IG).EQ.0.0_JWRB .AND. V(IPP,IG).EQ.0.0_JWRB) IPP = NLAND
          IPM = KLAT(IJ,1,1)
          IF (U(IPM,IG).EQ.0.0_JWRB .AND. V(IPM,IG).EQ.0.0_JWRB) IPM = NLAND
          IPP2 = KLAT(IJ,2,2)
          IF (U(IPP2,IG).EQ.0.0_JWRB .AND. V(IPP2,IG).EQ.0.0_JWRB) IPP2 = NLAND
          IPM2 = KLAT(IJ,1,2)
          IF (U(IPM2,IG).EQ.0.0_JWRB .AND. V(IPM2,IG).EQ.0.0_JWRB) IPM2 = NLAND

          IF (IPP.NE.NLAND .AND. IPM.NE.NLAND .AND.                     &
     &        IPP2.NE.NLAND .AND. IPM2.NE.NLAND) THEN
            UP = WLAT(IJ,2)*U(IPP,IG)+(1.0_JWRB-WLAT(IJ,2))*U(IPP2,IG)
            VP = WLAT(IJ,2)*V(IPP,IG)+(1.0_JWRB-WLAT(IJ,2))*V(IPP2,IG)
            UM = WLAT(IJ,1)*U(IPM,IG)+(1.0_JWRB-WLAT(IJ,1))*U(IPM2,IG)
            VM = WLAT(IJ,1)*V(IPM,IG)+(1.0_JWRB-WLAT(IJ,1))*V(IPM2,IG)
            DUPHI(IJ) = (UP-UM)*ONEO2DELPHI
            DVPHI(IJ) = (VP-VM)*ONEO2DELPHI
          ELSEIF (IPP.NE.NLAND .AND. IPM.NE.NLAND) THEN
            UP = U(IPP,IG)
            VP = V(IPP,IG)
            UM = U(IPM,IG)
            VM = V(IPM,IG)
            DUPHI(IJ) = (UP-UM)*ONEO2DELPHI
            DVPHI(IJ) = (VP-VM)*ONEO2DELPHI
          ELSEIF (IPP2.NE.NLAND .AND. IPM2.NE.NLAND ) THEN
            UP = U(IPP2,IG)
            VP = V(IPP2,IG)
            UM = U(IPM2,IG)
            VM = V(IPM2,IG)
            DUPHI(IJ) = (UP-UM)*ONEO2DELPHI
            DVPHI(IJ) = (VP-VM)*ONEO2DELPHI
          ELSE
            DUPHI(IJ) = 0.0_JWRB
            DVPHI(IJ) = 0.0_JWRB
          ENDIF

          ILP = KLON(IJ,2)
!         exact 0 means that the current field was not defined, hence
!         no gradient should be extrapolated
          IF (U(ILP,IG).EQ.0.0_JWRB .AND. V(ILP,IG).EQ.0.0_JWRB) ILP=NLAND
          ILM = KLON(IJ,1)
          IF (U(ILM,IG).EQ.0.0_JWRB .AND. V(ILM,IG).EQ.0.0_JWRB) ILM=NLAND
          KX  = KXLT(IJ,IG)
          IF (ILP.NE.NLAND .AND. ILM.NE.NLAND) THEN
            DULAM(IJ) = (U(ILP,IG)-U(ILM,IG))/(2.0_JWRB*DELLAM(KX))
            DVLAM(IJ) = (V(ILP,IG)-V(ILM,IG))/(2.0_JWRB*DELLAM(KX))
          ELSE
            DULAM(IJ) = 0.0_JWRB
            DVLAM(IJ) = 0.0_JWRB
          ENDIF
        ENDDO

        DO IJ=MIJS,MIJL
          DUPHI(IJ) = SIGN(MIN(ABS(DUPHI(IJ)),CURRENT_GRADIENT_MAX),DUPHI(IJ))
          DVPHI(IJ) = SIGN(MIN(ABS(DVPHI(IJ)),CURRENT_GRADIENT_MAX),DVPHI(IJ))
          DULAM(IJ) = SIGN(MIN(ABS(DULAM(IJ)),CURRENT_GRADIENT_MAX),DULAM(IJ))
          DVLAM(IJ) = SIGN(MIN(ABS(DVLAM(IJ)),CURRENT_GRADIENT_MAX),DVLAM(IJ))
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
