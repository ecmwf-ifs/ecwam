      SUBROUTINE STRESSO (F, IJS, IJL, THWNEW, USNEW, Z0NEW,
     &                    ROAIRN, TAUW, TAUWLF, PHIWA,
     &                    PHIAWDIAG, PHIAWUNR, SL, 
     &                    MIJ, LCFLX)

! ----------------------------------------------------------------------

!**** *STRESSO* - COMPUTATION OF WAVE STRESS.

!     H. GUNTHER      GKSS/ECMWF NOVEMBER   1989 CODE MOVED FROM SINPUT.
!     P.A.E.M. JANSSEN     KNMI  AUGUST     1990
!     J. BIDLOT            ECMWF FEBRUARY   1996-97
!     S. ABDALLA           ECMWF OCTOBER    2001 INTRODUCTION OF VARIABLE
!                                                AIR DENSITY
!     J. BIDLOT            ECMWF            2007  ADD MIJ
!     P.A.E.M. JANSSEN     ECMWF            2011  ADD FLUX CALULATIONS

!*    PURPOSE.
!     --------

!       COMPUTE NORMALIZED WAVE STRESS FROM INPUT SOURCE FUNCTION

!**   INTERFACE.
!     ----------

!       *CALL* *STRESSO (F, IJS, IJL, THWNEW, USNEW, Z0NEW,
!                        ROAIRN, TAUW, TAUWLF, PHIWA,
!                        PHIAWDIAG, PHIAWUNR, SL, 
!                        MIJ, LCFLX)*
!         *F*           - WAVE SPECTRUM.
!         *IJS*         - INDEX OF FIRST GRIDPOINT.
!         *IJL*         - INDEX OF LAST GRIDPOINT.
!         *THWNEW*      - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                         NOTATION (POINTING ANGLE OF WIND VECTOR,
!                         CLOCKWISE FROM NORTH).
!         *USNEW*       - NEW FRICTION VELOCITY IN M/S.
!         *Z0NEW*       - ROUGHNESS LENGTH IN M.
!         *ROAIRN*      - AIR DENSITY IN KG/M3.
!         *TAUW*        - WAVE STRESS IN (M/S)**2
!         *TAUWLF*      - LOW-FREQUENCY PART OF WAVE STRESS
!                         (i.e. INTEGRATED OVER  THE PROGNOSTIC RANGE)
!         *PHIWA*       - ENERGY FLUX FROM WIND INTO WAVES INTEGRATED
!                         OVER THE FULL DISCRETISED FREQUENCY RANGE.
!         *PHIAWDIAG*   - DIAGNOSTIC PART OF ENERGY FLUX INTO WAVES
!                         (i.e. INTEGRATED FROM THE END OF THE PROGNOSTIC RANGE
!                          TO THE LAST DISCRETISED FREQUENCY) 
!         *PHIAWUNR*    - UNRESOLVED PART OF ENERGY FLUX
!                         (i.e. INTEGRATED FROM THE LAST DISCRETISED FREQUENCY
!                          TO INFINITY).
!         *SL*          - TOTAL SOURCE FUNCTION ARRAY.
!         *MIJ*         - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!         *LCFLX*       - TRUE IF THE CALCULATION FOR THE FLUXES ARE 
!                         PERFORMED.

!     METHOD.
!     -------

!       THE INPUT SOURCE FUNCTION IS INTEGRATED OVER FREQUENCY
!       AND DIRECTIONS.
!       BECAUSE ARRAY *SL* IS USED, ONLY THE INPUT SOURCE
!       HAS TO BE STORED IN *SL* (CALL FIRST SINPUT, THEN
!       STRESSO, AND THEN THE REST OF THE SOURCE FUNCTIONS)

!     EXTERNALS.
!     -----------

!       NONE.

!     REFERENCE.
!     ----------
!       R SNYDER ET AL,1981.
!       G. KOMEN, S. HASSELMANN AND K. HASSELMANN, JPO, 1984.
!       P. JANSSEN, JPO, 1985

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : ALPHA    ,TAUWSHELTER, ITSHELT
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,TH       ,
     &            COSTH    ,SINTH    ,FR5      ,FRATIO
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER  ,YINVEPS
      USE YOWTABL  , ONLY : IUSTAR   ,IALPHA   ,EPS1     ,TAUHFT   ,
     &            TAUHFT2  ,ILEVTAIL ,DELUST   ,DELALP   ,DELTAIL
      USE YOWSTAT  , ONLY : IPHYS

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL
      INTEGER, DIMENSION(IJS:IJL), INTENT(IN) :: MIJ
      INTEGER :: IJ, M, K, I, J, II

      REAL,DIMENSION(IJS:IJL), INTENT(IN) :: THWNEW, USNEW, Z0NEW
      REAL,DIMENSION(IJS:IJL), INTENT(IN) :: ROAIRN
      REAL,DIMENSION(IJS:IJL), INTENT(OUT) :: TAUW, TAUWLF
      REAL,DIMENSION(IJS:IJL), INTENT(OUT) :: PHIWA, PHIAWDIAG, PHIAWUNR
      REAL,DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F
      REAL,DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: SL

      REAL :: CONST, ROG, DFIMLOC, CNST
      REAL :: XI, XJ, DELI1, DELI2, DELJ1, DELJ2, TAU1, XK, DELK1, DELK2
      REAL :: COSW, UST, SCDFM, SCDFP, BETAM, CONSTPHI
      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION(NFRE) :: ZPIROFR, CONSTFT
      REAL, DIMENSION(IJS:IJL) :: TAUHF, TEMP1, TEMP2, XSTRESS, YSTRESS
      REAL, DIMENSION(IJS:IJL) :: UST2, TAUX, TAUY, TAUPX, TAUPY
      REAL, DIMENSION(IJS:IJL) :: USDIRP, CONST0
      REAL, DIMENSION(IJS:IJL,NFRE) :: CONSTFM, CONSTFE

      LOGICAL, INTENT(IN) :: LCFLX

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('STRESSO',0,ZHOOK_HANDLE)
#endif

!*    1. PRECOMPUTE FREQUENCY SCALING.
!        -----------------------------

      CONST = DELTH*(ZPI)**4/G**2
      ROG   = ROWATER*G

      DO IJ=IJS,IJL
        CONST0(IJ)  = CONST*FR5(MIJ(IJ))
      ENDDO

      DO M=1,NFRE
        ZPIROFR(M) = ZPI*ROWATER*FR(M)
      ENDDO


!     !!!! CONSTFM is only defined up to M=MIJ(IJ)
      SCDFM = 0.5*DELTH*(1.-1./FRATIO)
      DO IJ=IJS,IJL
!       !!!! CONSTFM is only defined up to M=MIJ(IJ)
        DO M=1,MIJ(IJ)-1
          CONSTFM(IJ,M) = ZPIROFR(M)*DFIM(M)
        ENDDO
        CONSTFM(IJ,MIJ(IJ)) = ZPIROFR(MIJ(IJ))*SCDFM*FR(MIJ(IJ))
      ENDDO


!!!!!!!!! CONSTFT and CONSTFE are only used if LCFLX is true
      IF (LCFLX) THEN
        DO M=1,NFRE
          CONSTFT(M) = ROG*DFIM(M)
        ENDDO
!     !!!! CONSTFE is only defined from M=MIJ(IJ)
        SCDFP = 0.5*DELTH*(FRATIO-1.)
        DO IJ=IJS,IJL
          IF(MIJ(IJ).LT.NFRE) THEN
            CONSTFE(IJ,MIJ(IJ)) = ROG*SCDFP*FR(MIJ(IJ))
          ELSE
            CONSTFE(IJ,MIJ(IJ)) = 0. 
          ENDIF
          DO M=MIJ(IJ)+1,NFRE
            CONSTFE(IJ,M) = ROG*DFIM(M)
          ENDDO
        ENDDO
      ENDIF

!*    2. COMPUTE WAVE STRESS OF ACTUEL BLOCK.
!        ------------------------------------

!*    2.2 INTEGRATE INPUT SOURCE FUNCTION OVER FREQUENCY AND DIRECTIONS.
!         --------------------------------------------------------------

      DO IJ=IJS,IJL
        XSTRESS(IJ) = 0.
        YSTRESS(IJ) = 0.
        DO M=1,MIJ(IJ)
          DO K=1,NANG
            CNST = SL(IJ,K,M)*CONSTFM(IJ,M)
            XSTRESS(IJ) = XSTRESS(IJ)+CNST*SINTH(K)
            YSTRESS(IJ) = YSTRESS(IJ)+CNST*COSTH(K)
          ENDDO
        ENDDO
      ENDDO

      IF (LCFLX) THEN
        DO IJ=IJS,IJL
           TAUWLF(IJ) = SQRT(XSTRESS(IJ)**2+YSTRESS(IJ)**2)
        ENDDO

        PHIWA(IJS:IJL) = 0.
        PHIAWDIAG(IJS:IJL) = 0.
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=IJS,IJL
              PHIWA(IJ) = PHIWA(IJ)+SL(IJ,K,M)*CONSTFT(M)
              IF (M.GE.MIJ(IJ)) THEN
                PHIAWDIAG(IJ) = PHIAWDIAG(IJ)+SL(IJ,K,M)*CONSTFE(IJ,M)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      DO IJ=IJS,IJL
        XSTRESS(IJ) = XSTRESS(IJ)/MAX(ROAIRN(IJ),1.)
        YSTRESS(IJ) = YSTRESS(IJ)/MAX(ROAIRN(IJ),1.)
      ENDDO

!*    2.3 CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
!     ----------------------------------------------------

      IF (IPHYS.EQ.0 .OR. ITSHELT.EQ.0) THEN
        DO IJ=IJS,IJL
          USDIRP(IJ)=THWNEW(IJ)
        ENDDO
      ELSE
        DO IJ=IJS,IJL
          TAUX(IJ)=(USNEW(IJ)**2)*SIN(THWNEW(IJ))
          TAUY(IJ)=(USNEW(IJ)**2)*COS(THWNEW(IJ))
          TAUPX(IJ)=TAUX(IJ)-ABS(TAUWSHELTER)*XSTRESS(IJ)
          TAUPY(IJ)=TAUY(IJ)-ABS(TAUWSHELTER)*YSTRESS(IJ)
          USDIRP(IJ)=ATAN2(TAUPX(IJ),TAUPY(IJ))
        ENDDO
      ENDIF

      K=1
      DO IJ=IJS,IJL
        COSW     = MAX(COS(TH(K)-THWNEW(IJ)),0.)
        TEMP1(IJ) = F(IJ,K,MIJ(IJ))*COSW**3
        TEMP2(IJ) = F(IJ,K,NFRE)*COSW**2
      ENDDO

      DO K=2,NANG
        DO IJ=IJS,IJL
          COSW     = MAX(COS(TH(K)-THWNEW(IJ)),0.)
          TEMP1(IJ) = TEMP1(IJ)+F(IJ,K,MIJ(IJ))*COSW**3
          TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,NFRE)*COSW**2
        ENDDO
      ENDDO

      IF (ITSHELT.EQ.0) THEN
        DO IJ=IJS,IJL
          UST   = MAX(USNEW(IJ),0.000001)
          UST2(IJ) = UST**2
          XI    = UST / DELUST
          XI    = MIN(REAL(IUSTAR),XI)
          I     = MIN (IUSTAR-1, INT(XI))
          I     = MAX (0, I)
          DELI1 = MIN (1. ,XI-REAL(I))
          DELI2   = 1. - DELI1

          XJ    = (G*Z0NEW(IJ)/UST2(IJ)-ALPHA) / DELALP
          XJ    = MIN(REAL(IALPHA),XJ)
          J     = MIN (IALPHA-1, INT(XJ))
          J     = MAX (0, J)
          DELJ1 = MAX(MIN (1. ,XJ-REAL(J)),0.)
          DELJ2   = 1. - DELJ1

          TAU1 = ( TAUHFT(I  ,J  ,MIJ(IJ))*DELI2 +
     &             TAUHFT(I+1,J  ,MIJ(IJ))*DELI1 )*DELJ2
     &         + ( TAUHFT(I  ,J+1,MIJ(IJ))*DELI2 +
     &             TAUHFT(I+1,J+1,MIJ(IJ))*DELI1 )*DELJ1

          TAUHF(IJ) = CONST0(IJ)*TEMP1(IJ)*UST2(IJ)*TAU1
        ENDDO
      ELSE
        DO IJ=IJS,IJL
          UST   = MAX(USNEW(IJ),0.000001)
          UST2(IJ) = UST**2
          XI    = UST / DELUST
          XI    = MIN(REAL(IUSTAR),XI)
          I     = MIN (IUSTAR-1, INT(XI))
          I     = MAX (0, I)
          DELI1 = MIN (1. ,XI-REAL(I))
          DELI2   = 1. - DELI1

          XJ    = (G*Z0NEW(IJ)/UST2(IJ)-ALPHA) / DELALP
          XJ    = MIN(REAL(IALPHA),XJ)
          J     = MIN (IALPHA-1, INT(XJ))
          J     = MAX (0, J)
          DELJ1 = MAX(MIN (1. ,XJ-REAL(J)),0.)
          DELJ2   = 1. - DELJ1

          XK=CONST0(IJ)*TEMP1(IJ)/DELTAIL
          II=MIN(ILEVTAIL-1,INT(XK))
          DELK1= MIN (1. ,XK-FLOAT(II))
          DELK2=1. -DELK1

          TAU1 = ( (TAUHFT2(I  ,J  ,II  )*DELI2 +
     &              TAUHFT2(I+1,J  ,II  )*DELI1 )*DELJ2 +
     &             (TAUHFT2(I  ,J+1,II  )*DELI2 +
     &              TAUHFT2(I+1,J+1,II  )*DELI1)*DELJ1 ) * DELK2
     &           +
     &           ( (TAUHFT2(I  ,J  ,II+1)*DELI2 +
     &              TAUHFT2(I+1,J  ,II+1)*DELI1 )*DELJ2+
     &             (TAUHFT2(I  ,J+1,II+1)*DELI2 +
     &              TAUHFT2(I+1,J+1,II+1)*DELI1)*DELJ1 ) * DELK1

          TAUHF(IJ) = CONST0(IJ)*TEMP1(IJ)*UST2(IJ)*TAU1
        ENDDO
      ENDIF

      DO IJ=IJS,IJL
        XSTRESS(IJ) = XSTRESS(IJ)+TAUHF(IJ)*SIN(USDIRP(IJ))
        YSTRESS(IJ) = YSTRESS(IJ)+TAUHF(IJ)*COS(USDIRP(IJ))
        TAUW(IJ) = SQRT(XSTRESS(IJ)**2+YSTRESS(IJ)**2)

        TAUW(IJ) = MIN(TAUW(IJ),UST2(IJ)-EPS1)
        TAUW(IJ) = MAX(TAUW(IJ),0.)
      ENDDO
!
!*    4. UNRESOLVED PART ENERGY FLUX.
!        ----------------------------
!
      IF (LCFLX) THEN
        BETAM = 26.
        CONST = BETAM*ZPI**3*FR(NFRE)**4/G*DELTH

        DO IJ=IJS,IJL
          CONSTPHI     = CONST*ROAIRN(IJ) 
          PHIAWUNR(IJ) = CONSTPHI*USNEW(IJ)**2*TEMP2(IJ)
        ENDDO
      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('STRESSO',1,ZHOOK_HANDLE)
#endif

      RETURN
      END SUBROUTINE STRESSO
