      SUBROUTINE WEFLUX (F, IJS, IJL,
     &                   NFRE, NANG, DFIM, DELTH,
     &                   COSTH, SINTH, NDEPTH, INDEP, TCGOND,
     &                   WEFMAG, WEFDIR)
!
! ----------------------------------------------------------------------
!
!**** *WEFLUX* - COMPUTATION OF WAVE ENERGY FLUX AT EACH GRID POINT


!*    PURPOSE.
!     --------

!       COMPUTE WAVE ENERGY FLUX MAGNITUDE AND MEAN DIRECTION AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *WEFLUX* (F, IJS, IJL,
!                        NFRE, NANG, DFIM, DELTH,
!                        COSTH, SINTH, NDEPTH, INDEP, TCGOND,
!                        WEFMAG, WEFDIR)
!              *F*      - SPECTRUM.
!              *IJS*    - INDEX OF FIRST GRIDPOINT
!              *IJL*    - INDEX OF LAST GRIDPOINT
!              *NFRE*   - NUMBER OF FREQUENCIES
!              *NANG*   - NUMBER OF DIRECTIONS
!              *DFIM*   - FREQUENCY-DIRECTION INCREMENT
!              *DELTH*  - DIRECTIONAL INCREMENT
!              *COSTH*  - COS(THETA)
!              *SINTH*  - SIN(THETA)
!              *NDEPTH* - NUMBER OF DEPTHS
!              *INDEP*  - DEPTH INDEX ARRAY
!              *TCGOND* - GROUP VELOCITY TABLE
!              *WEFMAG* - WAVE ENERGY FLUX MAGNITUDE (W/m) 
!              *WEFDIR* - WAVE ENERGY FLUX MEAN DIRECTION (radian oceanographic convention) 
!

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.
!
!     REFERENCE.
!     ----------
!
!       NONE.
!
! ----------------------------------------------------------------------
 
      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER   ,EPSMIN
      USE YOWFRED  , ONLY : FRTAIL
 
! ----------------------------------------------------------------------
 
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL, NFRE, NANG, NDEPTH
      INTEGER,DIMENSION(IJS:IJL), INTENT(IN) :: INDEP

      REAL,INTENT(IN) :: DELTH
      REAL,DIMENSION(NFRE), INTENT(IN) :: DFIM
      REAL,DIMENSION(NANG), INTENT(IN) :: COSTH, SINTH
      REAL,DIMENSION(NDEPTH,NFRE), INTENT(IN) :: TCGOND 
      REAL,DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F

      REAL,DIMENSION(IJS:IJL), INTENT(OUT) :: WEFMAG, WEFDIR

      INTEGER :: IJ,M,K

      REAL :: ROG
      REAL :: FCG
      REAL :: DELT
      REAL,DIMENSION(IJS:IJL) :: TEMP, TEMPX, TEMPY, CGROUP
      REAL,DIMENSION(IJS:IJL) :: WEFX, WEFY 
!
! ----------------------------------------------------------------------

!*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
!        ------------------------------------------------

      DO IJ=IJS,IJL
         WEFMAG(IJ)= 0.
         WEFX(IJ)  = 0.
         WEFY(IJ)  = 0.
      ENDDO

      ROG = ROWATER*G  
      DELT = FRTAIL*DELTH*G/(2.*ZPI)
!
!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      DO M=1,NFRE
         DO IJ=IJS,IJL
            CGROUP(IJ)=TCGOND(INDEP(IJ),M)
         ENDDO
         K=1
         DO IJ=IJS,IJL
            FCG = F(IJ,K,M)*CGROUP(IJ)
            TEMP(IJ) = FCG
            TEMPX(IJ) = FCG*SINTH(K) 
            TEMPY(IJ) = FCG*COSTH(K) 
         ENDDO
         DO K=2,NANG
            DO IJ=IJS,IJL
              FCG = F(IJ,K,M)*CGROUP(IJ)
              TEMP(IJ) = TEMP(IJ)+FCG
              TEMPX(IJ) = TEMPX(IJ)+FCG*SINTH(K)
              TEMPY(IJ) = TEMPY(IJ)+FCG*COSTH(K)
            ENDDO
         ENDDO
         DO IJ=IJS,IJL
            WEFMAG(IJ)  = WEFMAG(IJ)+DFIM(M)*TEMP(IJ)
            WEFX(IJ)  = WEFX(IJ)+DFIM(M)*TEMPX(IJ)
            WEFY(IJ)  = WEFY(IJ)+DFIM(M)*TEMPY(IJ)
         ENDDO
      ENDDO


!*    3. ADD TAIL CORRECTION
!        -------------------

       K=1
       DO IJ=IJS,IJL
          FCG = F(IJ,K,NFRE)
          TEMP(IJ) = FCG
          TEMPX(IJ) = FCG*SINTH(K) 
          TEMPY(IJ) = FCG*COSTH(K) 
       ENDDO
       DO K=2,NANG
          DO IJ=IJS,IJL
            FCG = F(IJ,K,NFRE)
            TEMP(IJ) = TEMP(IJ)+FCG
            TEMPX(IJ) = TEMPX(IJ)+FCG*SINTH(K)
            TEMPY(IJ) = TEMPY(IJ)+FCG*COSTH(K)
          ENDDO
       ENDDO
       DO IJ=IJS,IJL
          WEFMAG(IJ) = WEFMAG(IJ)+DELT*TEMP(IJ)
          WEFX(IJ)  = WEFX(IJ)+DELT*TEMPX(IJ)
          WEFY(IJ)  = WEFY(IJ)+DELT*TEMPY(IJ)
       ENDDO

!*    4. MULTIPLY MAGNITUDE BY ROG:
!        -------------------------
      DO IJ=IJS,IJL
        WEFMAG(IJ)  = ROG*WEFMAG(IJ)
      ENDDO

!*    5. FLUX MEAN DIRECTION
!        -------------------

      DO IJ=IJS,IJL
        IF (WEFY(IJ).EQ.0.) WEFY(IJ) = EPSMIN
      ENDDO

      DO IJ=IJS,IJL
        WEFDIR(IJ)  = ATAN2(WEFX(IJ),WEFY(IJ)) 
      ENDDO

      DO IJ=IJS,IJL
        IF (WEFDIR(IJ).LT.0.) WEFDIR(IJ) = WEFDIR(IJ) + ZPI
      ENDDO

      END SUBROUTINE WEFLUX
