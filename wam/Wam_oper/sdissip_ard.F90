      SUBROUTINE SDISSIP_ARD (F, FL, IJS, IJL, SL,&
     &                        USNEW, THWNEW, ROAIRN,&
     &                        PHIEPS, TAUWD, MIJ, LCFLX)
! ----------------------------------------------------------------------

!**** *SDISSIP_ARD* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.

!     LOTFI AOUF       METEO FRANCE 2013
!     FABRICE ARDHUIN  IFREMER  2013


!*    PURPOSE.
!     --------
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SDISSIP (F, FL, IJS, IJL, SL,*
!                        USNEW, THWNEW,ROAIRN,*
!                        PHIEPS, TAUWD, MIJ, LCFLX)*
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY
!          *USNEW*  - NEW FRICTION VELOCITY IN M/S.
!          *ROAIRN* - AIR DENSITY IN KG/M3
!          *THWNEW* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC.

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
      USE YOWPCONS , ONLY : RAD     ,G        ,ZPI      ,ROWATER
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSHAL  , ONLY : TFAK     ,INDEP, TCGOND
      USE YOWSTAT  , ONLY : ISHALLO
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM), INTENT(IN) :: MIJ(IJS:IJL)
      INTEGER(KIND=JWIM) :: IJ, K , M
      INTEGER(KIND=JWIM) :: I, J, I1, J1
      INTEGER(KIND=JWIM) :: IS, SDSNTH, NDIKCUMUL
      INTEGER(KIND=JWIM) :: I_INT, J_INT, M2, K2, KK
      INTEGER(KIND=JWIM) :: NSMOOTH(IJS:IJL,NFRE)
      INTEGER(KIND=JWIM) :: ISDSDTH
      INTEGER(KIND=JWIM) :: ISSDSBRFDF
      INTEGER(KIND=JWIM), ALLOCATABLE :: SATINDICES(:,:)

      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: USNEW, THWNEW, ROAIRN 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: PHIEPS, TAUWD
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: FL,SL
      REAL(KIND=JWRB) :: XK(IJS:IJL,NFRE), CGOZPI(IJS:IJL,NFRE)
      REAL(KIND=JWRB) :: TPIINV, TMP01, TMP02, TMP03
      REAL(KIND=JWRB) :: DTURB
      REAL(KIND=JWRB) :: EPSR
      REAL(KIND=JWRB) :: MICHE
      REAL(KIND=JWRB) :: SCDFM
      REAL(KIND=JWRB) :: SDISS
      REAL(KIND=JWRB) :: SDSBR
      REAL(KIND=JWRB) :: SSDSC3, SSDSC4, SSDSC5, SSDSC6
      REAL(KIND=JWRB) :: SSDSBRF1, SXFR,SSDSBR,SSDSC2
      REAL(KIND=JWRB) :: DSIP_
      REAL(KIND=JWRB) :: ROG
      REAL(KIND=JWRB) :: SSDSC
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: SIG(NFRE)
      REAL(KIND=JWRB) :: SSDSC2_SIG(NFRE)
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: CM
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: FACTURB
      REAL(KIND=JWRB) :: RENEWALFREQ(IJS:IJL)
      REAL(KIND=JWRB) :: FACSAT(IJS:IJL,NFRE)
      REAL(KIND=JWRB) :: BTH0(IJS:IJL,NFRE)  !saturation spectrum 
      REAL(KIND=JWRB) :: BTH0S(IJS:IJL,NFRE)    !smoothed saturation spectrum 
      REAL(KIND=JWRB) :: BTH(IJS:IJL,NANG,NFRE)  !saturation spectrum 
      REAL(KIND=JWRB) :: BTHS(IJS:IJL,NANG,NFRE)  !smoothed saturation spectrum 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE) :: CONSTFM
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: A, D
      REAL(KIND=JWRB), DIMENSION(NFRE,NANG,NFRE,NANG) :: CUMULW
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:) :: C_,C_C,C2_,C2_C2,DSIP_05_C2
      REAL(KIND=JWRB), ALLOCATABLE :: SATWEIGHTS(:,:)

      LOGICAL, INTENT(IN) :: LCFLX
      LOGICAL :: LLSSDSC3,  LLSSDSC5

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SDISSIP_ARD',0,ZHOOK_HANDLE)
#endif

!     Br:
      SDSBR  = 9.0E-4_JWRB
      SSDSBR = SDSBR

!     Saturation dissipation coefficient
      SSDSC2 = -2.8E-5_JWRB
      SSDSC4 = 1.0_JWRB
      SSDSC6 = 0.3_JWRB
      ISDSDTH = 80

!     Cumulative dissipation coefficient
      SSDSC3 = 0.8_JWRB

!     Wave-turbulence interaction coefficient 
      SSDSC5  = 0.0_JWRB

!-----------------------------------------------------------------------------------

      EPSR=SQRT(SSDSBR)
      TPIINV = 1.0_JWRB/ZPI
      ROG = ROWATER*G

      LLSSDSC3=(SSDSC3.NE.0.0_JWRB)
      LLSSDSC5=(SSDSC5.NE.0.0_JWRB)


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

      DO M=1, NFRE
        SIG(M) = ZPI*FR(M)
      END DO

      IF (LLSSDSC5) THEN
        TMP01 = SSDSC5/ROG
        DO IJ=IJS,IJL
           FACTURB(IJ) = TMP01*ROAIRN(IJ)*USNEW(IJ)*USNEW(IJ)
        END DO
      ENDIF

      IF (ISHALLO.EQ.0) THEN
        DO M=1, NFRE
          DO IJ=IJS,IJL
            XK(IJ,M) = TFAK(INDEP(IJ),M)
            CGOZPI(IJ,M) = TPIINV*TCGOND(INDEP(IJ),M)
          ENDDO
        ENDDO
      ELSE
        DO M=1, NFRE
          DO IJ=IJS,IJL
            XK(IJ,M) = (SIG(M)**2)/G
            CGOZPI(IJ,M) = TPIINV*G/(2*ZPI*FR(M))
          ENDDO
        ENDDO
      ENDIF

      DO M=1, NFRE
        DO K=1, NANG
          DO IJ=IJS,IJL
            A(IJ,K,M) = CGOZPI(IJ,M)*F(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO


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
            SATWEIGHTS(K,I_INT-(K-SDSNTH)+1)=COS(TH(K)-TH(J_INT))**2
          END DO
        END DO
      END IF

!     compute cumulw
      SSDSBRF1   = 0.5_JWRB
      SXFR = 0.5_JWRB*(FRATIO-1/FRATIO)

      IF (LLSSDSC3) THEN
        DO J1=1,NANG
          DO I1=1,NFRE
            DO J=1,NANG
              DO I=1,NFRE
                CUMULW(I,J,I1,J1)=0.0_JWRB
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!        NDIKCUMUL is the  integer difference in frequency bands
!        between the "large breakers" and short "wiped-out waves"
        NDIKCUMUL = NINT(SSDSBRF1/(FRATIO-1.))
        ALLOCATE(C_(NFRE))
        ALLOCATE(C_C(NFRE))
        ALLOCATE(C2_(NFRE-NDIKCUMUL))
        ALLOCATE(C2_C2(NFRE-NDIKCUMUL))
        ALLOCATE(DSIP_05_C2(NFRE-NDIKCUMUL))
        DO M=1,NFRE  
          C_(M)=G/SIG(M)  ! Valid in deep water only
          C_C(M)=C_(M)*C_(M)
        END DO
        DO M=1,NFRE  
          DO K=1,NANG
            DO M2=1,M-NDIKCUMUL
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
             ENDDO
           ENDDO
         ENDDO
        ENDDO

      ENDIF

      DO  M=1, NFRE
        DO IJ=IJS,IJL
          FACSAT(IJ,M)=(XK(IJ,M)**3)*DELTH
        END DO
      ENDDO

      DO M=1,NFRE
        DO IJ=IJS,IJL
          BTH0(IJ,M) = 0.0_JWRB
        ENDDO
        DO K=1,NANG
          ! integrates around full circle
          DO IJ=IJS,IJL
            BTH(IJ,K,M)=0.0_JWRB
          ENDDO
          DO K2=1,SDSNTH*2+1
            KK=SATINDICES(K,K2)
            DO IJ=IJS,IJL
            BTH(IJ,K,M) = BTH(IJ,K,M) + SATWEIGHTS(K,K2)*A(IJ,KK,M)
            ENDDO
          ENDDO
          DO IJ=IJS,IJL
            BTH(IJ,K,M)=BTH(IJ,K,M)*FACSAT(IJ,M)
            BTH0(IJ,M)=MAX(BTH0(IJ,M),BTH(IJ,K,M))
          ENDDO
        ENDDO
      ENDDO

!/ST3      SDSBR     = 1.20E-3 ! Babanin (personnal communication)
!!!      ISSDSBRFDF  = 22    ! test pour DC
      ISSDSBRFDF  = 0

      IF (ISSDSBRFDF.GT.0 .AND. ISSDSBRFDF.LT.NFRE/2) THEN 
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

        DO M=1,NFRE
          DO IJ=IJS,IJL
            BTH0(IJ,M)=MAX(0.0_JWRB,BTH0S(IJ,M)/NSMOOTH(IJ,M))
          ENDDO
        ENDDO 
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=IJS,IJL
              BTH(IJ,K,M)=MAX(0.0_JWRB,BTHS(IJ,K,M)/NSMOOTH(IJ,M))
            ENDDO
          ENDDO
        ENDDO 
           
      ENDIF ! (ISSDSBRFDF.GT.0.AND.ISSDSBRFDF.LT.NFRE/2)

      MICHE = 1.0_JWRB
      TMP03 = 1.0_JWRB/(SSDSBR*MICHE)

      DO  M=1, NFRE
        SSDSC2_SIG(M)=SSDSC2*SIG(M)
      ENDDO

      ! Saturation term
      DO  M=1, NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            D(IJ,K,M)= &
     &      SSDSC2_SIG(M) &
     &    * (         SSDSC6 *(MAX(0._JWRB,BTH0(IJ,  M)*TMP03-SSDSC4))**2 &
     &    +  (1._JWRB-SSDSC6)*(MAX(0._JWRB,BTH (IJ,K,M)*TMP03-SSDSC4))**2)
          ENDDO
        ENDDO
      ENDDO

      ! Cumulative term
      IF (LLSSDSC3) THEN
        DO  M=NDIKCUMUL+1,NFRE
          DO K=1,NANG
            DO IJ=IJS,IJL
              RENEWALFREQ(IJ)=0.0_JWRB
            ENDDO
          ! Correction of saturation level for shallow-water kinematics
          ! Cumulative effect based on lambda   (breaking probability is
          ! the expected rate of sweeping by larger breaking waves)
            DO M2=1,M-NDIKCUMUL  
              DO K2=1,NANG
                DO IJ=IJS,IJL
                  IF (BTH0(IJ,M2).GT.SSDSBR) THEN
                  ! Integrates over frequencies M2 and directions K2 to 
                  ! Integration is performed from M2=1 to a frequency lower than M: IK-NDIKCUMUL
                    RENEWALFREQ(IJ)=RENEWALFREQ(IJ)+ CUMULW(M,K,M2,K2)*(MAX(SQRT(BTH(IJ,K2,M2))-EPSR,0.0_JWRB))**2
                  ENDIF
                ENDDO
              ENDDO
            ENDDO

            DO IJ=IJS,IJL
              D(IJ,K,M)= D(IJ,K,M)- SSDSC3*RENEWALFREQ(IJ)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!     Wave-turbulence interaction term
      IF (LLSSDSC5) THEN
        DO  M=1, NFRE
          DO K=1,NANG
            DO IJ=IJS,IJL
              D(IJ,K,M)= D(IJ,K,M)- 2._JWRB*SIG(M)*XK(IJ,M)*FACTURB(IJ)*COS(THWNEW(IJ)-COSTH(K))
            ENDDO
          ENDDO
        ENDDO
      ENDIF


      ! Add all contributions to sources term
      DO  M=1, NFRE
        DO K=1, NANG
          DO IJ=IJS,IJL
            SL(IJ,K,M) = SL(IJ,K,M)+D(IJ,K,M)*F(IJ,K,M)
            FL(IJ,K,M) = FL(IJ,K,M)+D(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO

      IF (LCFLX) THEN
        DO  M=1, NFRE
          IF (ISHALLO.EQ.0) THEN
            DO IJ=IJS,IJL
              CM(IJ)=TFAK(INDEP(IJ),M)/SIG(M)
            ENDDO
          ELSE
            DO IJ=IJS,IJL
              CM(IJ)=SIG(M)/G
            ENDDO
          ENDIF
          DO K=1, NANG
            DO IJ=IJS,IJL
              IF (M.LE.MIJ(IJ)) THEN
                SDISS = D(IJ,K,M)*F(IJ,K,M)
                PHIEPS(IJ) = PHIEPS(IJ)+SDISS*CONSTFM(IJ,M)
                TAUWD(IJ)  = TAUWD(IJ)+CM(IJ)*SDISS*CONSTFM(IJ,M)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF (ALLOCATED(SATWEIGHTS)) DEALLOCATE (SATWEIGHTS)
      IF (ALLOCATED(SATINDICES)) DEALLOCATE (SATINDICES)

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SDISSIP_ARD',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE SDISSIP_ARD
