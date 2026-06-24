SUBROUTINE FUSTAR (KIJS, KIJL, DEPTH, USMO, USA, EWFG, EWOI, EWA, T, FMWA)

! ----------------------------------------------------------------------

!**** *FUSTAR* - ESTIMATE THE ANALYSED FRICTION VELOCITY FROM THE 
!                DURATION AND THE "MEASURED" WINDSEA ENERGY. 
!                ESTIMATE THE ANALYSED WINDSEA MEAN FREQUENCY 
!                FROM THE "MEASURED" WINDSEA ENERGY 

!     P.LIONELLO     ECMWF       FEBRUARY 1989 
!         (FROM MODIFICATION OF A CODE BY P.JANSSEN 
!           AND P.LIONELLO - SUMMER '87 - ) 

!     PURPOSE. 
!     -------- 

!        ESTIMATE THE ANALYSED FRICTION VELOCITY FROM THE DURATION AND 
!        THE "MEASURED" WINDSEA ENERGY. ESTIMATE THE ANALYSED WINDSEA 
!        MEAN FREQUENCY FROM THE "MEASURED" WINDSEA ENERGY. 

!**   INTERFACE. 
!     ---------- 

!        *CALL* *FUSTAR (KIJS, KIJL, DEPTH, USMO, USA, EWFG, EWOI, EWA, T, FMWA)*

!        *KIJS*   FIRST INDEX IN BLOCK.                        
!        *KIJL*   LAST  INDEX IN BLOCK.                        
!        *DEPTH*  WATER DEPTH
!        *USMO*   FIRST GUESS FRICTION VELOCITY. 
!        *USA*    ANALYSED FRICTION VELOCITY. 
!        *EWFG*   WINDSEA ENERGY (FIRST GUESS). 
!        *EWOI*   ESTIMATED WIND-SEA ENERGY FROM MEASUREMENTS. 
!        *EWA*    ANALYSED WIND-SEA ENERGY. 
!        *T*      DURATION. 
!        *FMWA*   ANALYSED WINDSEA MEAN FREQUENCY. 

!      METHOD. 
!      -------                                                          

!       THE DURATION OF THE WIND SEA IS USED TO COMPUTE THE FRICTION  
!       VELOCITY BY NEWTONS METHOD FROM THE MODEL DURATION CURVE AND 
!       THE MEASURED ENERGY.  (FIVE ITERATIONS ARE DONE) 
!       THIS VALUE IS USED TO UPDATE THE MODEL IF IT IS 
!       IN REASONABLE AGREEMENT WITH THE MEASUREMENT FRICTION VELOCITY. 
!       IF THERE IS NO AGREEMENT THE MEASURED WIND SPEED IS USED  
!       TO COMPUTE THE WIND SEA. 

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWALTAS , ONLY : EGRCRV   ,AGRCRV   ,BGRCRV   ,              &
     &                      AFCRV   ,BFCRV   ,ESH      ,ASH     ,BSH
      USE YOWPCONS , ONLY : G       ,EPSUS
      USE YOWTABL  , ONLY : USTARM
      USE YOWTEST  , ONLY : IU06

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: USMO, EWFG, EWOI, T
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: USA, EWA, FMWA


      INTEGER(KIND=JWIM), PARAMETER :: NITER=5
      INTEGER(KIND=JWIM) :: IJ, ITER

      REAL(KIND=JWRB) :: YNU, EMAX, X, DELF
      REAL(KIND=JWRB) :: G2, BGRCRVM1, EWFGMIN
      REAL(KIND=JWRB) :: TSTAR, ESTAR, DSTAR, XX, YY, ZZ, ERATIO
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) ::  USTANAL, USTANAL0

      LOGICAL, DIMENSION(KIJS:KIJL) :: LLUPDT  

! ----------------------------------------------------------------------

!     INLINE FUNCTION.
!     ----------------

!     DIMENSIONLESS FREQUENCY AS FUNCTION OF DIMENSIONLESS ENERGY

      YNU(X) = (X/AFCRV)**(1.0_JWRB/BFCRV)

!     DIMENSIONLESS ENERGY LIMIT AS FUNCTION OF DIMENSIONLESS DEPTH 

      EMAX(X) = ESH*TANH(ASH*X**BSH)


! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FUSTAR',0,ZHOOK_HANDLE)

!*    1. INITIALIZE OUTPUT FIELDS.                                      
!        -------------------------                                      

      DO IJ=KIJS,KIJL
        USA(IJ)  =   -99.0_JWRB
        EWA(IJ)  =  -999.0_JWRB
        FMWA(IJ) = -9999.0_JWRB
      ENDDO

!*    2. LOOP OVER GRID POINTS.                                         
!        ----------------------                                         

      G2=G*G
      BGRCRVM1=BGRCRV-1.0_JWRB
      EWFGMIN=0.0001_JWRB

!*    2.1  DETERMINE FRICTION VELOCITY FROM DURATION GROWTH CURVE       
!          IF WINDSEA IS PRESENT.                                       
!           ------------------------------------------------------      
      DO IJ=KIJS,KIJL                                                
        IF (EWFG(IJ) > EWFGMIN .AND. T(IJ) > 0.0_JWRB) THEN
          USTANAL0(IJ) = MAX(USMO(IJ), EPSUS)
          USTANAL(IJ) = USTANAL0(IJ)
          LLUPDT(IJ) = .TRUE.  
        ELSE
          LLUPDT(IJ) = .FALSE.  
        ENDIF                                                          
      ENDDO

!     SOLVE ITERATIVELY THE RELATION FOR USTAR
!     ESTAR=EGRCRV*(TSTAR/(AGRCRV+TSTAR))**BGRCRV
!     USING THE 1 ORDER NEWTON METHOD Xn+1 = Xn - f(Xn)/f'(Xn)
!     WHERE f  = 0 IS A FUNCTION OBTAINED FROM THE PREVIOUS
!     RELATION IN TERM OF USTAR ONLY.

      DO ITER=1,NITER
        DO IJ=KIJS,KIJL
          IF (LLUPDT(IJ)) THEN
            TSTAR = G*T(IJ)/USTANAL(IJ)
            XX = AGRCRV/TSTAR
            YY = 1.0_JWRB+XX
            ERATIO = (EWOI(IJ)*G2/USTANAL(IJ)**4)/EGRCRV
            ZZ = (1.0_JWRB/YY)**BGRCRV

            DELF = 4.0_JWRB*ZZ - ERATIO*BGRCRV*XX/YY

            IF (DELF /= 0.0_JWRB) USTANAL(IJ) = USTANAL(IJ)*( 1.0_JWRB - (ZZ-ERATIO)/DELF )

            IF (USTANAL(IJ) <= 0.0_JWRB) THEN
              USTANAL0(IJ) = 2.0_JWRB*USTANAL0(IJ)
              USTANAL(IJ) = USTANAL0(IJ)
            ENDIF
            USTANAL(IJ) = MAX(USTANAL(IJ), EPSUS)
          ENDIF
        ENDDO
      ENDDO

!          -----------------------------------------------------

      DO IJ=KIJS,KIJL
        IF (LLUPDT(IJ)) THEN
          IF (ABS(USTANAL(IJ)-USMO(IJ)) < 0.5_JWRB*USMO(IJ)) THEN
            USA(IJ) = MIN(USTANAL(IJ),USTARM)
            EWA(IJ) = EWOI(IJ)
            ESTAR=EWOI(IJ)*G2/USTANAL(IJ)**4
            DSTAR = DEPTH(IJ)*G/USTANAL(IJ)**2
            ESTAR = MIN(EMAX(DSTAR), ESTAR)
            FMWA(IJ) = YNU(ESTAR)*G/USTANAL(IJ)
          ENDIF
        ENDIF
      ENDDO

IF (LHOOK) CALL DR_HOOK('FUSTAR',1,ZHOOK_HANDLE)

END SUBROUTINE FUSTAR
