      SUBROUTINE GRADI (IG, MIJS, MIJL,IREFRA, DDPHI, DDLAM, DUPHI,
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

      USE YOWCURR  , ONLY : U        ,V
      USE YOWGRID  , ONLY : DELPHI   ,DELLAM
      USE YOWMAP   , ONLY : KXLT
      USE YOWMPP   , ONLY : NINF
      USE YOWSHAL  , ONLY : DEPTH
      USE YOWSPEC  , ONLY : NSTART   ,NEND
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,WLAT

! ----------------------------------------------------------------------

      REAL, DIMENSION(MIJS:MIJL) :: DDPHI,DDLAM,DUPHI,DULAM,DVPHI,DVLAM

! ----------------------------------------------------------------------

!*    1. INITIALISE.
!        -----------

      NLAND=NINF-1

! ----------------------------------------------------------------------

!*    2. CALCULATE DEPTH GRADIENTS.
!        --------------------------

      IF (IREFRA.EQ.1 .OR. IREFRA.EQ.3) THEN
        DO IJ=MIJS,MIJL
          IPP = KLAT(IJ,2,1)
          IPM = KLAT(IJ,1,1)
          IPP2 = KLAT(IJ,2,2)
          IPM2 = KLAT(IJ,1,2)
          IF (IPP.NE.NLAND .AND. IPM.NE.NLAND .AND.
     &        IPP2.NE.NLAND .AND. IPM2.NE.NLAND) THEN
            DPTP=WLAT(IJ,2)*DEPTH(IPP,IG)+(1.-WLAT(IJ,2))*DEPTH(IPP2,IG)
            DPTM=WLAT(IJ,1)*DEPTH(IPM,IG)+(1.-WLAT(IJ,1))*DEPTH(IPM2,IG)
            DDPHI(IJ) = (DPTP-DPTM)/(2.*DELPHI)
          ELSEIF (IPP.NE.NLAND .AND. IPM.NE.NLAND ) THEN
            DPTP=DEPTH(IPP,IG)
            DPTM=DEPTH(IPM,IG)
            DDPHI(IJ) = (DPTP-DPTM)/(2.*DELPHI)
          ELSEIF (IPP2.NE.NLAND .AND. IPM2.NE.NLAND ) THEN
            DPTP=DEPTH(IPP2,IG)
            DPTM=DEPTH(IPM2,IG)
            DDPHI(IJ) = (DPTP-DPTM)/(2.*DELPHI)
          ELSE
            DDPHI(IJ) = 0.0
          ENDIF
          ILP = KLON(IJ,2)
          ILM = KLON(IJ,1)
          KX  = KXLT(IJ,IG)
          IF (ILP.NE.NLAND .AND. ILM.NE.NLAND) THEN
            DDLAM(IJ) = (DEPTH(ILP,IG)-DEPTH(ILM,IG))/(2.*DELLAM(KX))
          ELSE
            DDLAM(IJ) = 0.0 
          ENDIF
        ENDDO
      ELSE
        DO IJ=MIJS,MIJL
          DDPHI(IJ) = 0.0
          DDLAM(IJ) = 0.0
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    3. CALCULATE CURRENT VELOCITY GRADIENTS.
!        -------------------------------------

      IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN
        DO IJ=MIJS,MIJL
          IPP = KLAT(IJ,2,1)
          IPM = KLAT(IJ,1,1)
          IPP2 = KLAT(IJ,2,2)
          IPM2 = KLAT(IJ,1,2)
          IF (IPP.NE.NLAND .AND. IPM.NE.NLAND .AND.
     &        IPP2.NE.NLAND .AND. IPM2.NE.NLAND) THEN
            UP=WLAT(IJ,2)*U(IPP,IG)+(1.-WLAT(IJ,2))*U(IPP2,IG)
            VP=WLAT(IJ,2)*V(IPP,IG)+(1.-WLAT(IJ,2))*V(IPP2,IG)
            UM=WLAT(IJ,1)*U(IPM,IG)+(1.-WLAT(IJ,1))*U(IPM2,IG)
            VM=WLAT(IJ,1)*V(IPM,IG)+(1.-WLAT(IJ,1))*V(IPM2,IG)
            DUPHI(IJ) = (UP-UM)/(2.*DELPHI)
            DVPHI(IJ) = (VP-VM)/(2.*DELPHI)
          ELSEIF (IPP.NE.NLAND .AND. IPM.NE.NLAND ) THEN
            UP=U(IPP,IG)
            VP=V(IPP,IG)
            UM=U(IPM,IG)
            VM=V(IPM,IG)
            DUPHI(IJ) = (UP-UM)/(2.*DELPHI)
            DVPHI(IJ) = (VP-VM)/(2.*DELPHI)
          ELSEIF (IPP2.NE.NLAND .AND. IPM2.NE.NLAND ) THEN
            UP=U(IPP2,IG)
            VP=V(IPP2,IG)
            UM=U(IPM2,IG)
            VM=V(IPM2,IG)
            DUPHI(IJ) = (UP-UM)/(2.*DELPHI)
            DVPHI(IJ) = (VP-VM)/(2.*DELPHI)
          ELSE
            DUPHI(IJ) = 0.0
            DVPHI(IJ) = 0.0
          ENDIF
          ILP = KLON(IJ,2)
          ILM = KLON(IJ,1)
          KX  = KXLT(IJ,IG)
          IF (ILP.NE.NLAND .AND. ILM.NE.NLAND) THEN
            DULAM(IJ) = (U(ILP,IG)-U(ILM,IG))/(2.*DELLAM(KX))
            DVLAM(IJ) = (V(ILP,IG)-V(ILM,IG))/(2.*DELLAM(KX))
          ELSE
            DULAM(IJ) = 0.0
            DVLAM(IJ) = 0.0
          ENDIF
        ENDDO
      ELSE
        DO IJ=MIJS,MIJL
          DUPHI(IJ) = 0.0
          DVPHI(IJ) = 0.0
          DULAM(IJ) = 0.0
          DVLAM(IJ) = 0.0
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE GRADI
