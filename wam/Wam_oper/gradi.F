      SUBROUTINE GRADI (IG, IREFRA, DDPHI, DDLAM, DUPHI, DULAM,
     &                  DVPHI, DVLAM)

! ----------------------------------------------------------------------

!**** *GRADI* - CALCULATES DEPTH AND CURRENT VELOCITY GRADIENTS.

!     K.P. HUBBERT              AUGUST   1988
!     H. GUNTHER    ECMWF/GKSS  DECEMBER 1990  MODIFIED FOR CYCLE_4.

!*    PURPOSE.
!     --------

!       CALCULATES DEPTH AND CURRENT VELOCITY GRADIENTS OF A BLOCK.

!**   INTERFACE.
!     ----------

!       *CALL* *GRADI (IG, IREFRA, DDPHI, DDLAM, DUPHI, DULAM,
!                      DVPHI, DVLAM)*
!          *IG*     - BLOCK NUMBER.
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

      USE YOWCURR  , ONLY : U        ,V
      USE YOWMAP   , ONLY : KXLT
      USE YOWGRID  , ONLY : DELPHI   ,DELLAM   ,IJS      ,IJL
      USE YOWPARAM , ONLY : NIBLD    ,NIBLC
      USE YOWSHAL  , ONLY : DEPTH
      USE YOWUBUF  , ONLY : KLAT     ,KLON

! ----------------------------------------------------------------------

      DIMENSION DDPHI(NIBLD), DDLAM(NIBLD), DUPHI(NIBLC), DULAM(NIBLC),
     &          DVPHI(NIBLC), DVLAM(NIBLC)

! ----------------------------------------------------------------------

!*    1. INITIALISE.
!        -----------

! ----------------------------------------------------------------------

!*    2. CALCULATE DEPTH GRADIENTS.
!        --------------------------

      DO IJ=1,NIBLD
        DDPHI(IJ) = 0.0
        DDLAM(IJ) = 0.0
      ENDDO
      DO IJ=IJS(IG),IJL(IG)
        IPP = KLAT(IJ,2)
        IPM = KLAT(IJ,1)
        IF (IPP.GT.0 .AND. IPM.GT.0) THEN
          DDPHI(IJ) = (DEPTH(IPP,IG)-DEPTH(IPM,IG))/(2.*DELPHI)
        ENDIF
        ILP = KLON(IJ,2)
        ILM = KLON(IJ,1)
        KX  = KXLT(IJ,IG)
        IF (ILP.GT.0 .AND. ILM.GT.0) THEN
          DDLAM(IJ) = (DEPTH(ILP,IG)-DEPTH(ILM,IG))/(2.*DELLAM(KX))
        ENDIF
      ENDDO

! ----------------------------------------------------------------------

!*    3. CALCULATE CURRENT VELOCITY GRADIENTS.
!        -------------------------------------

      IF (IREFRA.EQ.2) THEN
        DO IJ=1,NIBLC
          DUPHI(IJ) = 0.0
          DULAM(IJ) = 0.0
          DVPHI(IJ) = 0.0
          DVLAM(IJ) = 0.0
        ENDDO
        DO IJ=IJS(IG),IJL(IG)
          IPP = KLAT(IJ,2)
          IPM = KLAT(IJ,1)
          IF (IPP.GT.0 .AND. IPM.GT.0) THEN
            DUPHI(IJ) = (U(IPP,IG)-U(IPM,IG))/(2.*DELPHI)
            DVPHI(IJ) = (V(IPP,IG)-V(IPM,IG))/(2.*DELPHI)
          ENDIF
          ILP = KLON(IJ,2)
          ILM = KLON(IJ,1)
          KX  = KXLT(IJ,IG)
          IF (ILP.GT.0 .AND. ILM.GT.0) THEN
            DULAM(IJ) = (U(ILP,IG)-U(ILM,IG))/(2.*DELLAM(KX))
            DVLAM(IJ) = (V(ILP,IG)-V(ILM,IG))/(2.*DELLAM(KX))
          ENDIF
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE GRADI
