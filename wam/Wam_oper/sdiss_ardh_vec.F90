      SUBROUTINE SDISS_ARDH_VEC (F, FL, IJS, IJL, SL, F1MEAN, XKMEAN,&
     &                           PHIEPS, TAUWD, MIJ, LCFLX)
! ----------------------------------------------------------------------

!**** *SDISSIP_ARDH_VEC* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.

!     LOTFI AOUF       METEO FRANCE 2013
!     FABRICE ARDHUIN  IFREMER  2013


!*    PURPOSE.
!     --------
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SDISSIP (F, FL, IJS, IJL, SL, F1MEAN, XKMEAN,*
!                        PHIEPS, TAUWD, MIJ, LCFLX)
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY
!          *F1MEAN* - MEAN FREQUENCY BASED ON 1st MOMENT.
!          *XKMEAN* - MEAN WAVE NUMBER BASED ON 1st MOMENT.
!          *PHIEPS* - ENERGY FLUX FROM WAVES TO OCEAN INTEGRATED OVER 
!                     THE PROGNOSTIC RANGE.
!          *TAUWD*  - DISSIPATION STRESS INTEGRATED OVER
!                     THE PROGNOSTIC RANGE.
!          *MIJ*    - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *LCFLX*  - TRUE IF THE CALCULATION FOR THE FLUXES ARE 
!                     PERFORMED.


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

      USE YOWFRED  , ONLY : FR, TH, FRATIO, DELTH, GOM, COSTH, SINTH, DFIM
      USE YOWPCONS , ONLY : RAD     ,G        ,ZPI      ,ROWATER  ,YEPS
      USE YOWMEAN  , ONLY : EMEAN, FMEAN
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSHAL  , ONLY : TFAK     ,INDEP, TCGOND
      USE YOWSPEC  , ONLY : U10NEW, THWNEW, USNEW
      USE YOWSTAT  , ONLY : ISHALLO
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM), INTENT(IN) :: MIJ(IJS:IJL)
      INTEGER(KIND=JWIM) :: IJ, K , M
      INTEGER(KIND=JWIM) :: I, J, I1, J1
      INTEGER(KIND=JWIM) :: IS, SDSNTH, DIKCUMUL
      INTEGER(KIND=JWIM) :: I_INT, J_INT, M2, K2
      INTEGER(KIND=JWIM) :: NSMOOTH(IJS:IJL,NFRE)
      INTEGER(KIND=JWIM) :: ISDSDTH
      INTEGER(KIND=JWIM) :: ISSDSBRFDF
      INTEGER(KIND=JWIM), ALLOCATABLE :: SATINDICES(:,:)

      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: F1MEAN, XKMEAN
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: PHIEPS, TAUWD
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: FL,SL
      REAL(KIND=JWRB) :: XK(IJS:IJL,NFRE), CG(IJS:IJL,NFRE)
      REAL(KIND=JWRB) :: ALFAMEAN
      REAL(KIND=JWRB) :: TPIINV, TMP00, TMP01, TMP02, TMP03, TMP04  
      REAL(KIND=JWRB) :: COSWIND
      REAL(KIND=JWRB) :: DTURB
      REAL(KIND=JWRB) :: EPSR
      REAL(KIND=JWRB) :: W, P0, MICHE, DELTA1, DELTA2
      REAL(KIND=JWRB) :: SCDFM
      REAL(KIND=JWRB) :: SDISS
      REAL(KIND=JWRB) :: SDSBR, SDSBR2
      REAL(KIND=JWRB) :: SATURATION2,FACSAT
      REAL(KIND=JWRB) :: SSDSC3, SSDSC4, SSDSC5, SSDSC6, SDSCOS
      REAL(KIND=JWRB) :: SSDSHF, SSDSLF, X, DTEMP, TEMP
      REAL(KIND=JWRB) :: SSDSBRF1, XFR, SXFR,SSDSBR,SDSC1,SSDSP,SSDSC2
      REAL(KIND=JWRB) :: DSIP_
      REAL(KIND=JWRB) :: ROG
      REAL(KIND=JWRB) :: SSDSC1, SSDSC , SSDSISO 
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: SIG(NFRE)
      REAL(KIND=JWRB) :: BSIGBAJ(IJS:IJL)
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: WNMEAN2, CM
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: FACTOR, FACTURB
      REAL(KIND=JWRB) :: RENEWALFREQ(IJS:IJL,NANG)
      REAL(KIND=JWRB) :: BTH0(IJS:IJL,NFRE)  !saturation spectrum 
      REAL(KIND=JWRB) :: BTH0S(IJS:IJL,NFRE)    !smoothed saturation spectrum 
      REAL(KIND=JWRB) :: BTH(IJS:IJL,NANG,NFRE)  !saturation spectrum 
      REAL(KIND=JWRB) :: BTHS(IJS:IJL,NANG,NFRE)  !smoothed saturation spectrum 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE) :: CONSTFM
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: A, D
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:) :: C_,C_C,C2_,C2_C2,DSIP_05_C2
      REAL(KIND=JWRB), ALLOCATABLE :: SATWEIGHTS(:,:)
      REAL(KIND=JWRB), ALLOCATABLE :: CUMULW(:,:,:,:)

      LOGICAL, INTENT(IN) :: LCFLX
      LOGICAL :: LLTEST

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SDISSIP_ARDH_VEC',0,ZHOOK_HANDLE)
#endif

      W=26.0_JWRB
      TPIINV = 1.0_JWRB/ZPI

      ROG = ROWATER*G

      IF (LCFLX) THEN
        PHIEPS(IJS:IJL) = 0.0_JWRB
        TAUWD(IJS:IJL) = 0.0_JWRB
!       !!!! CONSTFM is only defined up to M=MIJ(IJ)
        SCDFM = 0.5_JWRB*DELTH*(1.0_JWRB-1.0_JWRB/FRATIO)
        DO IJ=IJS,IJL
          DO M=1,MIJ(IJ)-1
            CONSTFM(IJ,M) = ROG*DFIM(M)
          ENDDO
          CONSTFM(IJ,MIJ(IJ)) = ROG*SCDFM*FR(MIJ(IJ))
        ENDDO
      ENDIF

      XFR = 1.1_JWRB !! ??? JEAN BIDLOT: what is XFR is it FRATIO ????
      SSDSBRF1   = 0.5_JWRB
      SXFR = 0.5_JWRB*(FRATIO-1/FRATIO)
!      SDSBR     = 1.40E-3_JWRB       
!      SDSBR     = 1.20E-3_JWRB       ! from Babanin (personal communication)
      SDSBR     = 9.0E-4_JWRB
      SDSBR2 = 0.8_JWRB
      SSDSBR = SDSBR

!      SDSC1  = -2.1_JWRB !! This is Bidlot et al. 2005,  Otherwise WAM4 uses -4.5
      SDSC1  = 0.0_JWRB
      SSDSC1  = SDSC1
      SSDSC5  = 0.0_JWRB
      SSDSISO  = 1.0_JWRB
      SSDSP = 2.0_JWRB
!      SSDSC2 = -3.0E-5_JWRB ! Retuned by Fabrice from VdW
      SSDSC2 = -2.8E-5_JWRB
      SSDSC3 = 0.8_JWRB
      SSDSC4 = 1.0_JWRB
      SSDSC6 = 0.3_JWRB
      ISDSDTH = 80
      SDSCOS = 2.0_JWRB
      SSDSHF = 1.0_JWRB
      SSDSLF = 1.0_JWRB
      EPSR=SQRT(SSDSBR)

      TMP00 = SSDSC1*ZPI
      TMP01 = SSDSC5/G*YEPS
      DO M=1, NFRE
        SIG(M) = ZPI*FR(M)
      END DO
      DO IJ=IJS,IJL
         ALFAMEAN = (XKMEAN(IJ)**2)*EMEAN(IJ)
         FACTOR(IJ) = TMP00     *F1MEAN(IJ)*(ALFAMEAN*ALFAMEAN)
         FACTURB(IJ) = TMP01*USNEW(IJ)*USNEW(IJ)
         WNMEAN2(IJ) = MAX( 1.E-10_JWRB , XKMEAN(IJ)  )
      END DO
      IF (ISHALLO.EQ.0) THEN
        DO M=1, NFRE
          DO IJ=IJS,IJL
            XK(IJ,M) = TFAK(INDEP(IJ),M)
            CG(IJ,M) = TCGOND(INDEP(IJ),M)
          END DO
        END DO
      ELSE
        DO M=1, NFRE
          DO IJ=IJS,IJL
            XK(IJ,M) = (SIG(M)**2)/G
            CG(IJ,M) = G/(2*ZPI*FR(M))
          END DO
        END DO
      ENDIF
      DO M=1, NFRE
!cdir outerunroll=4
        DO K=1, NANG
          DO IJ=IJS,IJL
            A(IJ,K,M) = TPIINV*CG(IJ,M)*F(IJ,K,M)
          END DO
        ENDDO
      END DO


      IF (ISDSDTH.NE.180) THEN
        SDSNTH  = MIN(NINT(ISDSDTH*RAD/(DELTH)),NANG/2-1)
        ALLOCATE(SATINDICES(NANG,SDSNTH*2+1))
        ALLOCATE(SATWEIGHTS(NANG,SDSNTH*2+1))
        DO K=1,NANG
          DO I_INT=K-SDSNTH, K+SDSNTH  
            J_INT=I_INT
            IF (I_INT.LT.1)  J_INT=I_INT+NANG
            IF (I_INT.GT.NANG) J_INT=I_INT-NANG
            SATINDICES(K,I_INT-(K-SDSNTH)+1)=J_INT
            SATWEIGHTS(K,I_INT-(K-SDSNTH)+1)=COS(TH(K)-TH(J_INT))**SDSCOS
          END DO
        END DO
      END IF

!     calcul de cumulw
      SSDSBRF1   = 0.5_JWRB
      SXFR = 0.5_JWRB*(FRATIO-1/FRATIO)
! initialise CUMULW
      ALLOCATE(CUMULW(NFRE,NANG,NFRE,NANG))
      DO I=1,NFRE
       DO J=1,NANG
         DO I1=1,NFRE
           DO J1=1,NANG
        CUMULW(I,J,I1,J1)=0.0_JWRB
           ENDDO
         ENDDO
       ENDDO
      END DO
      IF (SSDSC3.NE.0.0_JWRB) THEN

!        DIKCUMUL is the  integer difference in frequency bands
!        between the "large breakers" and short "wiped-out waves"
        DIKCUMUL = NINT(SSDSBRF1/(XFR-1.))
        ALLOCATE(C_(NFRE))
        ALLOCATE(C_C(NFRE))
        ALLOCATE(C2_(NFRE-DIKCUMUL))
        ALLOCATE(C2_C2(NFRE-DIKCUMUL))
        ALLOCATE(DSIP_05_C2(NFRE-DIKCUMUL))
        DO M=1,NFRE  
          C_(M)=G/SIG(M)  ! Valid in deep water only
          C_C(M)=C_(M)*C_(M)
        END DO
        DO M=1,NFRE  
          DO K=1,NANG
        DO M2=1,M-DIKCUMUL
          C2_(M2)=G/SIG(M2)
          C2_C2(M2)=C2_(M2)*C2_(M2)
          DSIP_ = SIG(M2)*SXFR
          DSIP_05_C2(M2)=DSIP_/(0.5_JWRB*C2_(M2))
            DO K2=1,NANG
                CUMULW(M,K,M2,K2)=SQRT(C_C(M)+C2_C2(M2)    &
     &          -2*C_(M)*C2_(M2)*COSTH(1+ABS(K2-K)))*DSIP_05_C2(M2) 
              END DO
            END DO 
          END DO
        END DO
        DEALLOCATE(C_)
        DEALLOCATE(C_C)
        DEALLOCATE(C2_)
        DEALLOCATE(C2_C2)
        DEALLOCATE(DSIP_05_C2)

! Multiplies by lambda(k,theta)=1/(2*pi**2) and 
! and the coefficient that transforms  SQRT(B) to Banner et al. (2000)'s epsilon
! 2.26 is equal to 5.55 (Banner & al. 2000) times 1.6**2 / 2pi where
! 1.6 is the ratio between Banner's epsilon and SQRT(B)

        TMP02 = 2*TPIINV*2.26
        DO I=1,NFRE
         DO J=1,NANG 
           DO I1=1,NFRE
             DO J1=1,NANG
               CUMULW(I,J,I1,J1)=CUMULW(I,J,I1,J1)*TMP02
             END DO
           ENDDO
         END DO
        END DO

      END IF

        DO  M=1, NFRE
          DO IJ=IJS,IJL
            FACSAT=(XK(IJ,M)**3)*DELTH
            BTH0(IJ,M)=SUM(A(IJ,1:NANG,M))*FACSAT
          END DO
        END DO
!      DO K=1,NANG
      DO  M=1, NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            FACSAT=(XK(IJ,M)**3)*DELTH
            ! integrates around full circle
            BTH(IJ,K,M) = SUM(SATWEIGHTS(K,1:SDSNTH*2+1)*A(IJ,SATINDICES(K,1:SDSNTH*2+1),M))*FACSAT
          END DO
        END DO
          DO IJ=IJS,IJL
            BTH0(IJ,M) = MAXVAL(BTH(IJ,1:NANG,M))
          END DO
      END DO

!/ST3      SDSBR     = 1.20E-3 ! Babanin (personnal communication)
      ISSDSBRFDF  = 22    ! test pour DC
      ISSDSBRFDF  = 0
!/ST3      SDSBRF1   = 0.5
!/ST3      SDSBRF2   = 0.
      IF (ISSDSBRFDF.GT.0.AND.ISSDSBRFDF.LT.NFRE/2) THEN 
!cdir collapse
        BTH0S(:,:)=BTH0(:,:)
!cdir collapse
        NSMOOTH(:,:)=1
!cdir collapse
        BTHS(:,:,:)=BTH(:,:,:)
!cdir outerunroll=4
        DO M=1, ISSDSBRFDF
          DO IJ=IJS,IJL
            BTH0S  (IJ,1+ISSDSBRFDF)=BTH0S  (IJ,1+ISSDSBRFDF)+BTH0(IJ,M)
            NSMOOTH(IJ,1+ISSDSBRFDF)=NSMOOTH(IJ,1+ISSDSBRFDF)+1
          END DO 
!cdir collapse
          DO K=1,NANG       
            DO IJ=IJS,IJL
              BTHS(IJ,K,M)=BTHS(IJ,K,M)+BTH(IJ,K,M)
            END DO
          END DO 
        ENDDO
        DO M=2+ISSDSBRFDF,1+2*ISSDSBRFDF
!cdir nodep
          DO IJ=IJS,IJL
            BTH0S  (IJ,1+ISSDSBRFDF)=BTH0S  (IJ,1+ISSDSBRFDF)+BTH0(IJ,M)
            NSMOOTH(IJ,1+ISSDSBRFDF)=NSMOOTH(IJ,1+ISSDSBRFDF)+1
          END DO 
!cdir collapse
          DO K=1,NANG       
            DO IJ=IJS,IJL
              BTHS(IJ,K,M)=BTHS(IJ,K,M)+BTH(IJ,K,M)
            END DO
          END DO
        ENDDO
        DO M=ISSDSBRFDF,1,-1
!cdir nodep
          DO IJ=IJS,IJL
            BTH0S  (IJ,M)=BTH0S  (IJ,M+1)-BTH0(IJ,M+ISSDSBRFDF+1)
            NSMOOTH(IJ,M)=NSMOOTH(IJ,M+1)-1
          END DO
!cdir collapse
          DO K=1,NANG
            DO IJ=IJS,IJL
              BTHS(IJ,K,M)=BTHS(IJ,K,M)-BTH(IJ,K,M)
            ENDDO
          END DO
        ENDDO
        DO M=2+ISSDSBRFDF,NFRE-ISSDSBRFDF
!cdir nodep
          DO IJ=IJS,IJL
            BTH0S  (IJ,M)=BTH0S  (IJ,M-1)-BTH0(IJ,M-ISSDSBRFDF-1)+BTH0(IJ,M+ISSDSBRFDF)
            NSMOOTH(IJ,M)=NSMOOTH(IJ,M-1)
          END DO
!cdir collapse
          DO K=1,NANG       
            DO IJ=IJS,IJL
              BTHS(IJ,K,M)=BTHS(IJ,K,M)-BTH(IJ,K,M)+BTH(IJ,K,M)
            END DO
          END DO
        ENDDO
!cdir novector
        DO M=NFRE-ISSDSBRFDF+1,NFRE
!cdir nodep
          DO IJ=IJS,IJL
            BTH0S  (IJ,M)=BTH0S  (IJ,M-1)-BTH0(IJ,M-ISSDSBRFDF)
            NSMOOTH(IJ,M)=NSMOOTH(IJ,M-1)-1
          END DO
!cdir collapse
          DO K=1,NANG       
            DO IJ=IJS,IJL
              BTHS(IJ,K,M)=BTHS(IJ,K,M)-BTH(IJ,K,M)
            END DO
          END DO
        END DO

!  final division by NSMOOTH

!cdir collapse
        DO M=1,NFRE
          DO IJ=IJS,IJL
            BTH0(IJ,M)=MAX(0.0_JWRB,BTH0S(IJ,M)/NSMOOTH(IJ,M))
          END DO
        END DO 
        DO M=1,NFRE
!cdir outerunroll=4
          DO K=1,NANG
            DO IJ=IJS,IJL
              BTH(IJ,K,M)=MAX(0.0_JWRB,BTHS(IJ,K,M)/NSMOOTH(IJ,M))
            END DO
          END DO
        END DO 
           
      END IF

!      DELTA1 = 0.4_JWRB
!      DELTA2 = 0.6_JWRB
      DELTA1 = 0.0_JWRB
      DELTA2 = 0.0_JWRB
      MICHE = 1.0_JWRB
      TMP03 = 1.0_JWRB/(SSDSBR*MICHE)

      DO  M=1, NFRE
        LLTEST = (SSDSC3.NE.0.AND.M.GT.DIKCUMUL)

        DO IJ=IJS,IJL
          IF (XKMEAN(IJ).NE.0) THEN
            X         = TFAK(INDEP(IJ),M)/XKMEAN(IJ)
            BSIGBAJ(IJ) = FACTOR(IJ)*( (1.0_JWRB-DELTA2)*X + DELTA2*X**2)
          ELSE
            BSIGBAJ(IJ) = 0.0_JWRB
          ENDIF
        END DO

        IF (ISHALLO.EQ.0) THEN
          DO IJ=IJS,IJL
            CM(IJ)=TFAK(INDEP(IJ),M)/SIG(M)
          ENDDO
        ELSE
          DO IJ=IJS,IJL
            CM(IJS)=SIG(M)/G
          ENDDO
        ENDIF

        DO K=1,NANG

              DO IJ=IJS,IJL
                  RENEWALFREQ(IJ,K)=0.0_JWRB
              ENDDO
          ! Correction of saturation level for shallow-water kinematics
          ! Cumulative effect based on lambda   (breaking probability is
          ! the expected rate of sweeping by larger breaking waves)
          IF (LLTEST) THEN
            DO M2=1,M-DIKCUMUL  
              IF (BTH0(IJ,M2).GT.SSDSBR) THEN
                DO K2=1,NANG
                  DO IJ=IJS,IJL
                  ! Integrates over frequencies M2 and directions K2 to 
                  ! Integration is performed from M2=1 to a frequency lower than M: IK-DIKCUMUL
                 RENEWALFREQ(IJ,K)=RENEWALFREQ(IJ,K)+ CUMULW(M,K,M2,K2) &
     &              *(MAX(SQRT(BTH(IJ,K2,M2))-EPSR,0.0_JWRB))**2
                  END DO
                END DO
              ENDIF
            END DO
          ENDIF

          DO IJ=IJS,IJL
            SATURATION2=TANH(10*(((BTH(IJ,K,M)/SSDSBR)**0.5_JWRB)-SDSBR2))
            COSWIND=(COSTH(K)*COS(THWNEW(IJ))+SINTH(K)*SIN(THWNEW(IJ)))   ! ÃvÃ©rifier K ?
            DTURB=-2.*SIG(M)*XK(IJ,M)*FACTURB(IJ)*COSWIND  ! Theory -> stress direction
            P0=SSDSP ! -0.5*SSDSC3*(1-TANH(W*USTAR*XK(IJ,M)/SIG(M)-0.1))  ! for SDSC3=1 this is vdW et al. 

            TMP04 = SSDSC3*RENEWALFREQ(IJ,K)
!            DTEMP=SSDSC2 * SIG(M) &
            DTEMP=SSDSC2 * SIG(M) &
     &    * (  SSDSC6 *(MAX(0._JWRB,BTH0(IJ,    M)*TMP03-SSDSC4))**P0 &
     &    + (1-SSDSC6)*(MAX(0._JWRB,BTH (IJ,K,M)*TMP03-SSDSC4))**P0)&
     &    - (TMP04+DTURB)  !terme cumulatif
 
            D(IJ,K,M) = DTEMP + BSIGBAJ(IJ)*SSDSLF *0.5_JWRB*(1-SATURATION2) &
     &                    + BSIGBAJ(IJ)*SSDSHF *0.5_JWRB*(SATURATION2+1)
 
          END DO
        END DO


!cdir outerunroll=4
        DO K=1, NANG
          DO IJ=IJS,IJL
            SL(IJ,K,M) = SL(IJ,K,M)+D(IJ,K,M)*F(IJ,K,M)
            FL(IJ,K,M) = FL(IJ,K,M)+D(IJ,K,M)
            IF (LCFLX.AND.M.LE.MIJ(IJ)) THEN
              SDISS = D(IJ,K,M)*F(IJ,K,M)
              PHIEPS(IJ) = PHIEPS(IJ)+SDISS*CONSTFM(IJ,M)
              TAUWD(IJ)  = TAUWD(IJ)+CM(IJS)*SDISS*CONSTFM(IJ,M)
            ENDIF
          END DO
        END DO
      END DO

      IF (ALLOCATED(CUMULW)) DEALLOCATE (CUMULW)
      IF (ALLOCATED(SATWEIGHTS)) DEALLOCATE (SATWEIGHTS)
      IF (ALLOCATED(SATINDICES)) DEALLOCATE (SATINDICES)

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SDISSIP_ARDH_VEC',1,ZHOOK_HANDLE)
#endif

      RETURN
      END SUBROUTINE SDISS_ARDH_VEC
