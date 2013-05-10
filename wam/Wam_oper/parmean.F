      SUBROUTINE PARMEAN (F,EM,THQ,MEANWP1,FP)

! ----------------------------------------------------------------------

!**** *PARMEAN* - COMPUTATION OF TOTAL ENERGY, MEAN DIRECTION,
!                  AND MEAN PERIOD BASED ON THE FIRST MOMENT

!      J-R BIDLOT    ECMWF     MARCH 2000
!      D PETTENUZZO  MAY 2012 MERGED 3 ROUTINES FOR PAR COMPUTATION

!*    PURPOSE.
!     --------

!       COMPUTE TOTAL ENERGY, MEAN DIRECTION AND MEAN PERIOD.

!**   INTERFACE.
!     ----------

!       *CALL* *PARMEAN (F,EM,THQ,MEANWP1,FP)*
!              *F*       - SPECTRUM.
!              *EM*      - MEAN WAVE ENERGY FOR THE WAVE SYSTEM 
!              *THQ*     - MEAN WAVE DIRECTION
!              *MEANWP1* - MEAN PERIOD BASED ON THE FIRST MOMENT.
!              *FP*      - PEAK FREQUENCY OF THE 1D SPECTRUM.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWFRED  ,ONLY : FR ,DFIM ,DFFR ,DELTH,WP1TAIL,COSTH,
     &                     SINTH,WETAIL
      USE YOWPARAM ,ONLY : NANG,NFRE
      USE YOWPCONS ,ONLY : EPSMIN,ZPI

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: K,M,MMAX
      REAL    :: F(NANG,NFRE)
      REAL    :: TEMP2(NFRE)
      REAL    :: TEMP, EMEAN, MEANWP1, THQ, FP, SI, CI, EM
      REAL    :: XMAX, XP1, XM1, F10, F1P1, F1M1
      REAL    :: A, B, C
      REAL    :: COEF_FR, DELT25

! ----------------------------------------------------------------------


!!!!! need to re-write to combine similar operations

!     1D FREQUENCy SPECTRUM
      DO M=1,NFRE
        K=1
        TEMP2(M) = F(K,M)
        DO K=2,NANG
          TEMP2(M) = TEMP2(M)+F(K,M)
        ENDDO
      ENDDO

!*    TOTAL ENERGY
!     ------------

!*    1. Initialise energy array.
!        ------------------------

      EM = EPSMIN

! ----------------------------------------------------------------------

!*    2. Integrate over frequencies and direction.
!        -----------------------------------------

      DO M=1,NFRE
        EM = EM+DFIM(M)*TEMP2(M)
      ENDDO

! ----------------------------------------------------------------------

!*    3. Add tail energy.
!        ----------------

      DELT25 = WETAIL*FR(NFRE)*DELTH
      EM = EM+DELT25*TEMP2(NFRE)

! ----------------------------------------------------------------------

!*    MEAN WAVE DIRECTION
!     -------------------

!*    1. Initialize sin and cos arrays
!        ------------------------------

      SI = 0.
      CI = 0.

! ----------------------------------------------------------------------

!*    2. Integrate over frequencies and directions.
!        ------------------------------------------

      DO K=1,NANG
        TEMP = F(K,1)*DFIM(1)
        DO M=2,NFRE
          TEMP = TEMP+F(K,M)*DFIM(M)
        ENDDO
        SI = SI+SINTH(K)*TEMP
        CI = CI+COSTH(K)*TEMP
      ENDDO

! ----------------------------------------------------------------------

!*    3. Compute mean direction.
!        -----------------------

      IF (CI.EQ.0.) CI = EPSMIN
      THQ = ATAN2(SI,CI)
      IF (THQ.LT.0.) THQ = THQ + ZPI

! ----------------------------------------------------------------------

!*    MEAN PERIOD BASED ON FIRST MOMENT
!     ---------------------------------

      MEANWP1 = EPSMIN 

      DO M=1,NFRE
        MEANWP1 = MEANWP1+DFFR(M)*TEMP2(M)
      ENDDO

!*    ADD TAIL CORRECTION TO MEAN PERIOD AND
!     AND TRANSFORM TO PERIOD.

      COEF_FR = WP1TAIL*DELTH*FR(NFRE)**2

      MEANWP1 = MEANWP1+COEF_FR*TEMP2(NFRE)
      IF(EM.GT.EPSMIN .AND. MEANWP1.GT.EPSMIN ) THEN
        MEANWP1 = EM/MEANWP1
      ELSE
        MEANWP1 = 0.
      ENDIF

!*    PEAK FREQUENCY
!     --------------
      XMAX = 0.
      MMAX = 1

      DO M=1,NFRE
        IF (TEMP2(M).GT.XMAX) THEN
          MMAX = M
          XMAX = TEMP2(M)
        ENDIF
      ENDDO

      FP= FR(MMAX)

      IF (XMAX.GT.0..AND.MMAX.LT.NFRE.AND.MMAX.GT.1) THEN
        XP1 = FR(MMAX+1)-FR(MMAX)
        XM1 = FR(MMAX-1)-FR(MMAX)
        F10  = TEMP2(MMAX) 
        F1P1 = TEMP2(MMAX+1)-F10
        F1M1 = TEMP2(MMAX-1)-F10

        A = F10
        B = 1./(XM1-XP1)*(XM1*F1P1/XP1-XP1*F1M1/XM1)
        C = 1./(XM1-XP1)*(F1M1/XM1-F1P1/XP1)

        IF (C.LT.0.) THEN
          FP = FR(MMAX)-B/(2.*C)
        ENDIF
      ENDIF


! ----------------------------------------------------------------------

      RETURN
      END SUBROUTINE PARMEAN
