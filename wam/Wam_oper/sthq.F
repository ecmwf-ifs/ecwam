      SUBROUTINE STHQ (F3, IJS, IJL, THQ)

! ----------------------------------------------------------------------

!**** *STHQ* - COMPUTATION OF MEAN WAVE DIRECTION AT EACH GRID POINT.

!     S.D. HASSELMANN
!     OPTIMIZED BY L. ZAMBRESKY

!*    PURPOSE.
!     --------

!       TO COMPUTE MEAN WAVE DIRECTION AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *STHQ (F3, IJS, IJL, THQ)*
!          *F3*  - SPECTRUM.
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *THQ* - MEAN WAVE DIRECTION

!     METHOD.
!     -------

!       INTEGRATION OF SPECTRUM TIMES SIN AND COS OVER DIRECTION.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : ZPI      ,EPSMIN
      USE YOWFRED  , ONLY : DFIM     ,COSTH    ,SINTH
      USE YOWMPP   , ONLY : NINF     ,NSUP

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IJ,M,K,IJS,IJL
      REAL :: F3(NINF-1:NSUP,NANG,NFRE)
      REAL, DIMENSION(IJS:IJL) :: TEMP, THQ, SI, CI

! ----------------------------------------------------------------------

!*    1. INITIALISE SIN AND COS ARRAYS.
!        ------------------------------

      DO IJ=IJS,IJL
        SI(IJ) = 0.
        CI(IJ) = 0.
      ENDDO

! ----------------------------------------------------------------------

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      DO K=1,NANG
        DO IJ=IJS,IJL
          TEMP(IJ) = 0.
        ENDDO
        DO M=1,NFRE
          DO IJ=IJS,IJL
            TEMP(IJ) = TEMP(IJ)+F3(IJ,K,M)*DFIM(M)
          ENDDO
        ENDDO
        DO IJ=IJS,IJL
          SI(IJ) = SI(IJ)+SINTH(K)*TEMP(IJ)
          CI(IJ) = CI(IJ)+COSTH(K)*TEMP(IJ)
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    3. COMPUTE MEAN DIRECTION.
!        -----------------------

      DO IJ=IJS,IJL
        IF (CI(IJ).EQ.0.) CI(IJ) = EPSMIN
      ENDDO
      DO IJ=IJS,IJL
        THQ(IJ) = ATAN2(SI(IJ),CI(IJ))
      ENDDO
      DO IJ=IJS,IJL
        IF (THQ(IJ).LT.0.) THQ(IJ) = THQ(IJ) + ZPI
      ENDDO

      RETURN
      END SUBROUTINE STHQ
