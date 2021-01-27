      SUBROUTINE SINPUT_ARD (NGST,F,FL,IJS,IJL,THWNEW,USNEW,Z0NEW,&
     &                       ROAIRN,WSTAR,RNFAC,SL,SPOS,XLLWS)
! ----------------------------------------------------------------------

!**** *SINPUT_ARD* - COMPUTATION OF INPUT SOURCE FUNCTION.


!*    PURPOSE.
!     ---------

!       COMPUTE THE WIND INPUT SOURCE TRERM BASED ON ARDHUIN ET AL. 2010.

!       COMPUTE INPUT SOURCE FUNCTION AND STORE ADDITIVELY INTO NET
!       SOURCE FUNCTION ARRAY, ALSO COMPUTE FUNCTIONAL DERIVATIVE OF
!       INPUT SOURCE FUNCTION.
!
!       GUSTINESS IS INTRODUCED FOLL0WING THE APPROACH OF JANSSEN(1986),
!       USING A GAUSS-HERMITE APPROXIMATION SUGGESTED BY MILES(1997).
!       IN THE PRESENT VERSION ONLY TWO HERMITE POLYNOMIALS ARE UTILISED
!       IN THE EVALUATION OF THE PROBABILITY INTEGRAL. EXPLICITELY ONE THEN
!       FINDS:
!
!             <GAMMA(X)> = 0.5*( GAMMA(X(1+SIG)) + GAMMA(X(1-SIG)) )
!
!       WHERE X IS THE FRICTION VELOCITY AND SIG IS THE RELATIVE GUSTINESS
!       LEVEL.

!**   INTERFACE.
!     ----------

!     *CALL* *SINPUT_ARD (NGST, F, FL, IJS, IJL, U10NEW, THWNEW, USNEW, Z0NEW,
!    &                   ROAIRN,WSTAR, RNFAC, SL, SPOS, XLLWS)
!            *NGST* - IF = 1 THEN NO GUSTINESS PARAMETERISATION
!                   - IF = 2 THEN GUSTINESS PARAMETERISATION
!            *F* - SPECTRUM.
!           *FL* - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!          *IJS* - INDEX OF FIRST GRIDPOINT.
!          *IJL* - INDEX OF LAST GRIDPOINT.
!       *THWNEW* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!        *USNEW* - NEW FRICTION VELOCITY IN M/S.
!        *Z0NEW* - ROUGHNESS LENGTH IN M.
!       *ROAIRN* - AIR DENSITY IN KG/M3
!        *WSTAR* - FREE CONVECTION VELOCITY SCALE (M/S).
!        *RNFAC* - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
!           *SL* - TOTAL SOURCE FUNCTION ARRAY.
!         *SPOS* - POSITIVE SOURCE FUNCTION ARRAY.
!         *XLLWS*- 1 WHERE SINPUT IS POSITIVE

!     METHOD.
!     -------

!       SEE REFERENCE.

!     EXTERNALS.
!     ----------

!       WSIGSTAR.

!     MODIFICATIONS
!     -------------

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LLCAPCHNK,LLNORMAGAM
      USE YOWFRED  , ONLY : FR       ,TH       ,DFIM     ,COSTH  ,SINTH, ZPIFR
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO
      USE YOWPCONS , ONLY : G        ,GM1      ,ROWATER  ,EPSMIN, EPSUS
      USE YOWPHYS  , ONLY : ZALP     ,TAUWSHELTER, XKAPPA, BETAMAXOXKAPPA2,    &
     &                      BMAXOKAPDTH, RN1_RN, &
     &                      RNU      ,RNUM, &
     &                      SWELLF   ,SWELLF2  ,SWELLF3  , SWELLF4  , SWELLF5, &
     &                      SWELLF6  ,SWELLF7  ,Z0RAT    , Z0TUBMAX , ABMIN  ,ABMAX
      USE YOWSHAL  , ONLY : TFAK     ,CINV     ,INDEP    ,TCGOND
      USE YOWSTAT  , ONLY : ISHALLO
!debile debug
      USE YOWSTAT  , ONLY : cdtpro 

      USE YOWTEST  , ONLY : IU06
      USE YOWTABL  , ONLY : IAB      ,SWELLFT
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "wsigstar.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: NGST
      INTEGER(KIND=JWIM), INTENT(IN) :: IJS,IJL

      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: THWNEW, USNEW, Z0NEW
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: ROAIRN, WSTAR, RNFAC
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(OUT) :: FL, SL, SPOS
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(OUT) :: XLLWS


      INTEGER(KIND=JWIM) :: IJ,K,M,IND,IGST

      REAL(KIND=JWRB) :: ROG
      REAL(KIND=JWRB) :: AVG_GST, ABS_TAUWSHELTER 
      REAL(KIND=JWRB) :: CONST1
      REAL(KIND=JWRB) :: UFAC0, ZN 
      REAL(KIND=JWRB) :: X1,X2,ZLOG,ZLOG1,ZLOG2,ZLOG2X,XV1,XV2,ZBETA1,ZBETA2
      REAL(KIND=JWRB) :: XI,X,DELI1,DELI2
      REAL(KIND=JWRB) :: FU,FUD,NU_AIR,SMOOTH, HFTSWELLF6,Z0TUB
      REAL(KIND=JWRB) :: FAC_NU_AIR,FACM1_NU_AIR
      REAL(KIND=JWRB) :: ARG, DELABM1
      REAL(KIND=JWRB) :: TAUPX,TAUPY
      REAL(KIND=JWRB) :: DSTAB2
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NFRE) :: CONST, SIG, SIGM1, SIG2, COEF, COEF5, DFIM_SIG2
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: CONSTF, CONST11, CONST22
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: Z0VIS, Z0NOZ, FWW
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: PVISC, PTURB
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: ZCN
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: SIG_N, UORBT, AORB, TEMP, RE, RE_C, ZORB
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: CNSN, SUMF
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: FLP_AVG, SLP_AVG
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: GZ0, ROGOROAIR, ROAIRN_PVISC
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: BMAXFAC
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NGST) :: XSTRESS, YSTRESS, FLP, SLP
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NGST) :: USG2, TAUX, TAUY, USTP, USTPM1, USDIRP, UCN
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NGST) :: UCNZALPD
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE) :: CM, XK
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE) :: XNGAMCONST
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NGST) :: GAMNORMA ! ! RENORMALISATION FACTOR OF THE GROWTH RATE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: DSTAB1, TEMP1, TEMP2
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NGST) :: COSLP, UFAC, DSTAB

      LOGICAL :: LTAUWSHELTER
! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SINPUT_ARD',0,ZHOOK_HANDLE)

      ROG = ROWATER*G
      AVG_GST = 1.0_JWRB/NGST
      CONST1  = BETAMAXOXKAPPA2/ROWATER
      NU_AIR = RNU
      FAC_NU_AIR= RNUM
      FACM1_NU_AIR=4.0_JWRB/NU_AIR

      FU=ABS(SWELLF3)
      FUD=SWELLF2
      DELABM1= REAL(IAB)/(ABMAX-ABMIN)

      ABS_TAUWSHELTER=ABS(TAUWSHELTER)
      IF(ABS_TAUWSHELTER .EQ. 0.0_JWRB ) THEN
        LTAUWSHELTER = .FALSE.
      ELSE
        LTAUWSHELTER = .TRUE.
      ENDIF

      IF(LLNORMAGAM) THEN
        BMAXFAC(:) = BMAXOKAPDTH * RNFAC(:)
        IF (ISHALLO.EQ.1) THEN
          DO M=1,NFRE
            DO IJ=IJS,IJL
              XNGAMCONST(IJ,M) = BMAXFAC(IJ)*0.5_JWRB*ZPIFR(M)**3*FR(M)*GM1
            ENDDO
          ENDDO
        ELSE
          DO M=1,NFRE
            DO IJ=IJS,IJL
              XNGAMCONST(IJ,M) = BMAXFAC(IJ)*FR(M)*TFAK(INDEP(IJ),M)**2*TCGOND(INDEP(IJ),M)
            ENDDO
          ENDDO
        ENDIF
      ELSE
        XNGAMCONST(:,:) = 0.0_JWRB
      ENDIF


!     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
      IF(NGST.GT.1) CALL WSIGSTAR (IJS, IJL, USNEW, Z0NEW, WSTAR, SIG_N)

! ----------------------------------------------------------------------
! computation of Uorb and Aorb

      DO IJ=IJS,IJL
        UORBT(IJ) = EPSMIN 
        AORB(IJ) = EPSMIN
      ENDDO

      DO M=1,NFRE
        SIG(M) = ZPIFR(M)
        SIGM1(M) = 1.0_JWRB/SIG(M)
        SIG2(M) = SIG(M)**2
        DFIM_SIG2(M)=DFIM(M)*SIG2(M)
        CONST(M)=SIG(M)*CONST1
        COEF(M) =-SWELLF*16._JWRB*SIG2(M)/(G*ROWATER)
        COEF5(M)=-SWELLF5*2._JWRB*SQRT(2._JWRB*NU_AIR*SIG(M))/ROWATER

        K=1
        DO IJ=IJS,IJL
          TEMP(IJ) = F(IJ,K,M)
        ENDDO
        DO K=2,NANG
          DO IJ=IJS,IJL
            TEMP(IJ) = TEMP(IJ)+F(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=IJS,IJL
          UORBT(IJ) = UORBT(IJ)+DFIM_SIG2(M)*TEMP(IJ)
          AORB(IJ) = AORB(IJ)+DFIM(M)*TEMP(IJ)
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        UORBT(IJ) = 2.0_JWRB*SQRT(UORBT(IJ))  ! this is the significant orbital amplitude
        AORB(IJ)  = 2.0_JWRB*SQRT(AORB(IJ))   ! this 1/2 Hs
        RE(IJ)    = FACM1_NU_AIR*UORBT(IJ)*AORB(IJ) ! this is the Reynolds number 
        Z0VIS(IJ) = FAC_NU_AIR/MAX(USNEW(IJ),0.0001_JWRB)
        Z0TUB = Z0RAT*MIN(Z0TUBMAX,Z0NEW(IJ))
        Z0NOZ(IJ) = MAX(Z0VIS(IJ),Z0TUB)
        ZORB(IJ)  = AORB(IJ)/Z0NOZ(IJ)

! compute fww
        XI=(LOG10(MAX(ZORB(IJ),3.0_JWRB))-ABMIN)*DELABM1
        IND  = MIN (IAB-1, INT(XI))
        DELI1= MIN (1.0_JWRB ,XI-REAL(IND,JWRB))
        DELI2= 1.0_JWRB - DELI1
        FWW(IJ)= SWELLFT(IND)*DELI2+SWELLFT(IND+1)*DELI1
        TEMP2(IJ) = FWW(IJ)*UORBT(IJ)
      ENDDO

! Define the critical Reynolds number
      IF( SWELLF6 .EQ. 1.0_JWRB) THEN
        DO IJ=IJS,IJL
          RE_C(IJ) = SWELLF4
        ENDDO
      ELSE
        HFTSWELLF6=1.0_JWRB-SWELLF6
        DO IJ=IJS,IJL
          RE_C(IJ) = SWELLF4*(2.0_JWRB/AORB(IJ))**HFTSWELLF6
        ENDDO
      ENDIF

! Swell damping weight between viscous and turbulent boundary layer
      IF (SWELLF7.GT.0.0_JWRB) THEN
        DO IJ=IJS,IJL
          SMOOTH=0.5_JWRB*TANH((RE(IJ)-RE_C(IJ))/SWELLF7)
          PTURB(IJ)=0.5_JWRB+SMOOTH
          PVISC(IJ)=0.5_JWRB-SMOOTH
        ENDDO
      ELSE
        DO IJ=IJS,IJL
          IF (RE(IJ).LE.RE_C(IJ)) THEN
            PTURB(IJ)=0.0_JWRB
            PVISC(IJ)=0.5_JWRB
          ELSE
            PTURB(IJ)=0.5_JWRB
            PVISC(IJ)=0.0_JWRB
          ENDIF
        ENDDO
      ENDIF


! Initialisation

      IF(NGST.EQ.1) THEN
        DO IJ=IJS,IJL
          USTP(IJ,1) = USNEW(IJ)
        ENDDO
      ELSE IF (NGST.EQ.2) THEN
        DO IJ=IJS,IJL
          USTP(IJ,1) = USNEW(IJ)*(1.0_JWRB+SIG_N(IJ))
          USTP(IJ,2) = USNEW(IJ)*(1.0_JWRB-SIG_N(IJ))
        ENDDO
      ELSE
         WRITE (IU06,*) '**************************************'
         WRITE (IU06,*) '*    FATAL ERROR                     *'
         WRITE (IU06,*) '*    ===========                     *'
         WRITE (IU06,*) '* IN SINPUT_ARD: NGST > 2            *'
         WRITE (IU06,*) '* NGST = ', NGST
         WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.  *'
         WRITE (IU06,*) '*                                    *'
         WRITE (IU06,*) '**************************************'
         CALL ABORT1
      ENDIF

      DO IGST=1,NGST
        DO IJ=IJS,IJL
          USTPM1(IJ,IGST) = 1.0_JWRB/MAX(USTP(IJ,IGST),EPSUS)
        ENDDO
      ENDDO

      IF(LTAUWSHELTER) THEN
        DO IGST=1,NGST
          DO IJ=IJS,IJL
            XSTRESS(IJ,IGST)=0.0_JWRB
            YSTRESS(IJ,IGST)=0.0_JWRB
            USG2(IJ,IGST)=USTP(IJ,IGST)**2
            TAUX(IJ,IGST)=USG2(IJ,IGST)*SIN(THWNEW(IJ))
            TAUY(IJ,IGST)=USG2(IJ,IGST)*COS(THWNEW(IJ))
          ENDDO
        ENDDO
      ELSE
        DO IGST=1,NGST
          DO K=1,NANG
            DO IJ=IJS,IJL
              COSLP(IJ,K,IGST) = COS(TH(K)-THWNEW(IJ))
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      DO IJ=IJS,IJL
        GZ0(IJ) = G*Z0NEW(IJ)
      ENDDO
      DO IJ=IJS,IJL
        ROGOROAIR(IJ) = ROG/MAX(ROAIRN(IJ),1.0_JWRB)
      ENDDO

      DO IJ=IJS,IJL
        ROAIRN_PVISC(IJ) = ROAIRN(IJ)*PVISC(IJ)
      ENDDO


!     INVERSE OF PHASE VELOCITIES AND WAVE NUMBER.
      IF (ISHALLO.EQ.1) THEN
        DO M=1,NFRE
          DO IJ=IJS,IJL
            CM(IJ,M) = SIG(M)/G
            XK(IJ,M) = SIG2(M)/G
          ENDDO
        ENDDO
      ELSE
        DO M=1,NFRE
          DO IJ=IJS,IJL
            XK(IJ,M) = TFAK(INDEP(IJ),M)
            CM(IJ,M) = XK(IJ,M)*SIGM1(M)
          ENDDO
        ENDDO
      ENDIF

!*    2. MAIN LOOP OVER FREQUENCIES.
!        ---------------------------
      DO M=1,NFRE

        IF(LTAUWSHELTER) THEN
          DO IGST=1,NGST
            DO IJ=IJS,IJL
              TAUPX=TAUX(IJ,IGST)-ABS_TAUWSHELTER*XSTRESS(IJ,IGST)
              TAUPY=TAUY(IJ,IGST)-ABS_TAUWSHELTER*YSTRESS(IJ,IGST)
              USDIRP(IJ,IGST)=ATAN2(TAUPX,TAUPY)
              USTP(IJ,IGST)=(TAUPX**2+TAUPY**2)**0.25_JWRB
              USTPM1(IJ,IGST)=1.0_JWRB/MAX(USTP(IJ,IGST),EPSUS)
            ENDDO
          ENDDO

          DO IGST=1,NGST
            DO K=1,NANG
              DO IJ=IJS,IJL
                COSLP(IJ,K,IGST) = COS(TH(K)-USDIRP(IJ,IGST))
              ENDDO
            ENDDO
          ENDDO

          DO IJ=IJS,IJL
            CONSTF(IJ) = ROGOROAIR(IJ)*CINV(INDEP(IJ),M)*DFIM(M)
          ENDDO
        ENDIF


!*      PRECALCULATE FREQUENCY DEPENDENCE.
!       ----------------------------------

        DO IGST=1,NGST
          DO IJ=IJS,IJL
            UCN(IJ,IGST) = USTP(IJ,IGST)*CM(IJ,M)
            UCNZALPD(IJ,IGST) = XKAPPA/(UCN(IJ,IGST) + ZALP)
          ENDDO
        ENDDO
        DO IJ=IJS,IJL
          ZCN(IJ) = LOG(XK(IJ,M)*Z0NEW(IJ))
          CNSN(IJ) = CONST(M)*ROAIRN(IJ)
        ENDDO

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------

        DO K=1,NANG
          DO IJ=IJS,IJL
            XLLWS(IJ,K,M)=0.0_JWRB
          ENDDO
        ENDDO

        DO IGST=1,NGST
          DO K=1,NANG
            DO IJ=IJS,IJL
              IF (COSLP(IJ,K,IGST).GT.0.01_JWRB) THEN
                X    = COSLP(IJ,K,IGST)*UCN(IJ,IGST)
                ZLOG = ZCN(IJ) + UCNZALPD(IJ,IGST)/COSLP(IJ,K,IGST)
                IF (ZLOG.LT.0.0_JWRB) THEN
                  ZLOG2X=ZLOG*ZLOG*X
                  UFAC(IJ,K,IGST) = EXP(ZLOG)*ZLOG2X*ZLOG2X
                  XLLWS(IJ,K,M) = 1.0_JWRB
                ELSE
                  UFAC(IJ,K,IGST) = 0.0_JWRB
                ENDIF
              ELSE
                UFAC(IJ,K,IGST) = 0.0_JWRB
              ENDIF
            ENDDO
          ENDDO
        ENDDO

      IF(LLNORMAGAM) THEN

!       windsea part of the spectrum
        SUMF(:) = 0.0_JWRB
        DO K=1,NANG
          DO IJ=IJS,IJL
            SUMF(IJ) = SUMF(IJ) + XLLWS(IJ,K,M)*F(IJ,K,M)*MAX(COS(TH(K)-THWNEW(IJ)),0.0_JWRB)**2
          ENDDO
        ENDDO

!       Computes the growth rate in the wind direction
        DO IGST=1,NGST
            DO IJ=IJS,IJL
              X    = UCN(IJ,IGST)
              ZLOG = ZCN(IJ) + UCNZALPD(IJ,IGST)
              ZLOG = MIN(ZLOG,0.0_JWRB)
              ZLOG2X = ZLOG*ZLOG*X
              UFAC0 = ZLOG2X*ZLOG2X*EXP(ZLOG)
              ZN = XNGAMCONST(IJ,M)*UFAC0*USTPM1(IJ,IGST)*SUMF(IJ)
              GAMNORMA(IJ,IGST) = (1.0_JWRB + RN1_RN*ZN)/(1.0_JWRB + ZN)
            ENDDO
        ENDDO
      ELSE
        GAMNORMA(:,:) = 1.0_JWRB
      ENDIF


!       SWELL DAMPING:

        DO IJ=IJS,IJL
          DSTAB1(IJ) = COEF5(M)*ROAIRN_PVISC(IJ)*XK(IJ,M)
          TEMP1(IJ) = COEF(M)*ROAIRN(IJ)
        ENDDO

        DO IGST=1,NGST
          DO K=1,NANG
            DO IJ=IJS,IJL
              DSTAB2 = TEMP1(IJ)*(TEMP2(IJ)+(FU+FUD*COSLP(IJ,K,IGST))*USTP(IJ,IGST))
              DSTAB(IJ,K,IGST) = DSTAB1(IJ)+PTURB(IJ)*DSTAB2
            ENDDO
          ENDDO
        ENDDO


!*    2.2 UPDATE THE SHELTERING STRESS (in any),
!         AND THEN ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
!         ---------------------------------------------------------

        DO K=1,NANG

          DO IGST=1,NGST
            DO IJ=IJS,IJL
              ! SLP: only the positive contributions
              SLP(IJ,IGST) =  UFAC(IJ,K,IGST)*CNSN(IJ) * GAMNORMA(IJ,IGST)
              FLP(IJ,IGST) = SLP(IJ,IGST)+DSTAB(IJ,K,IGST)
            ENDDO
          ENDDO

          DO IGST=1,NGST
            DO IJ=IJS,IJL
!!debile debug
       if (IGST == 1 .and. K == 1 ) then
       write(iu06,'(a9,1x,a12,1x,2(f14.8,1x))') 'debile_lf',cdtpro,TFAK(INDEP(IJ),M),SLP(IJ,IGST)/ZPIFR(M)
       endif

              SLP(IJ,IGST) = SLP(IJ,IGST)*F(IJ,K,M)
            ENDDO
          ENDDO

          IF(LTAUWSHELTER) THEN
            DO IJ=IJS,IJL
              CONST11(IJ)=CONSTF(IJ)*SINTH(K)
              CONST22(IJ)=CONSTF(IJ)*COSTH(K)
            ENDDO
            DO IGST=1,NGST
              DO IJ=IJS,IJL
                XSTRESS(IJ,IGST)=XSTRESS(IJ,IGST)+SLP(IJ,IGST)*CONST11(IJ)
                YSTRESS(IJ,IGST)=YSTRESS(IJ,IGST)+SLP(IJ,IGST)*CONST22(IJ)
              ENDDO
            ENDDO
          ENDIF

          IGST=1
            DO IJ=IJS,IJL
              SLP_AVG(IJ) = SLP(IJ,IGST)
              FLP_AVG(IJ) = FLP(IJ,IGST)
            ENDDO
          DO IGST=2,NGST
            DO IJ=IJS,IJL
              SLP_AVG(IJ) = SLP_AVG(IJ)+SLP(IJ,IGST)
              FLP_AVG(IJ) = FLP_AVG(IJ)+FLP(IJ,IGST)
            ENDDO
          ENDDO

          DO IJ=IJS,IJL
            SPOS(IJ,K,M) = AVG_GST*SLP_AVG(IJ)
          ENDDO

          DO IJ=IJS,IJL
            FL(IJ,K,M) = AVG_GST*FLP_AVG(IJ)
            SL(IJ,K,M) = FL(IJ,K,M)*F(IJ,K,M)
          ENDDO

        ENDDO

      ENDDO ! END LOOP OVER FREQUENCIES

      IF (LHOOK) CALL DR_HOOK('SINPUT_ARD',1,ZHOOK_HANDLE)

      END SUBROUTINE SINPUT_ARD
