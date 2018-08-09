      SUBROUTINE WEFLUX (F, IJS, IJL,                                   &
     &                   NFRE, NANG, DFIM, DELTH,                       &
     &                   COSTH, SINTH, CGROUP,                          &
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
!                        COSTH, SINTH, CGROUP,
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
!              *CGROUP* - GROUP VELOCITIES
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

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
 
      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER   ,EPSMIN
      USE YOWFRED  , ONLY : FRTAIL
 
! ----------------------------------------------------------------------
 
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, NFRE, NANG

      REAL(KIND=JWRB), INTENT(IN) :: DELTH
      REAL(KIND=JWRB), DIMENSION(NFRE), INTENT(IN) :: DFIM
      REAL(KIND=JWRB), DIMENSION(NANG), INTENT(IN) :: COSTH, SINTH
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE), INTENT(IN) :: CGROUP
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F

      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: WEFMAG, WEFDIR

      INTEGER(KIND=JWIM) :: IJ, M, K

      REAL(KIND=JWRB) :: ROG, FCG, DELT
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: TEMP, TEMPX, TEMPY
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: WEFX, WEFY 
!
! ----------------------------------------------------------------------

!*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
!        ------------------------------------------------

      DO IJ=IJS,IJL
        WEFMAG(IJ)= 0.0_JWRB
        WEFX(IJ)  = 0.0_JWRB
        WEFY(IJ)  = 0.0_JWRB
      ENDDO

      ROG = ROWATER*G  
      DELT = FRTAIL*DELTH*G/(2.0_JWRB*ZPI)
!
!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      DO M=1,NFRE
         K=1
         DO IJ=IJS,IJL
           FCG = F(IJ,K,M)*CGROUP(IJ,M)
           TEMP(IJ) = FCG
           TEMPX(IJ) = FCG*SINTH(K) 
           TEMPY(IJ) = FCG*COSTH(K) 
         ENDDO
         DO K=2,NANG
            DO IJ=IJS,IJL
              FCG = F(IJ,K,M)*CGROUP(IJ,M)
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
        IF (WEFY(IJ).EQ.0.0_JWRB) WEFY(IJ) = EPSMIN
      ENDDO

      DO IJ=IJS,IJL
        WEFDIR(IJ)  = ATAN2(WEFX(IJ),WEFY(IJ)) 
      ENDDO

      DO IJ=IJS,IJL
        IF (WEFDIR(IJ).LT.0.0_JWRB) WEFDIR(IJ) = WEFDIR(IJ) + ZPI
      ENDDO

      END SUBROUTINE WEFLUX
