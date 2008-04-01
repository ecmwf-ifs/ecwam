      SUBROUTINE FEMEAN (F, IJS, IJL, EM, FM, LLAK)

! ----------------------------------------------------------------------

!**** *FEMEAN* - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT
!                AND MEAN WAVE NUMBER .
!                THE COMPUTATION OF THE MEAN WAVE ENERGU WAS ALSO
!                ADDED SUCH THAT A CALL TO FEMEAN DOES NOT NEED
!                TO BE PRECEDED BY A CALL TO SEMEAN.

!     S.D. HASSELMANN
!     MODIFIED : P.JANSSEN (INTEGRATION OF F**-4 TAIL)
!     OPTIMIZED BY : L. ZAMBRESKY AND H. GUENTHER


!*    PURPOSE.
!     --------

!       COMPUTE MEAN FREQUENCY AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *FEMEAN (F, IJS, IJL, EM, FM, LLAK)*
!              *F*   - SPECTRUM.
!              *IJS* - INDEX OF FIRST GRIDPOINT
!              *IJL* - INDEX OF LAST GRIDPOINT
!              *EM*  - MEAN WAVE ENERGY (INPUT)
!              *FM*  - MEAN WAVE FREQUENCY (OUTPUT)
!              *LLAK*- TRUE IF MEAN WAVE NUMBER IS COMPUTED

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

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFIMOFR  ,DELTH    ,
     &                WETAIL    ,FRTAIL
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN
      USE YOWSTAT  , ONLY : ISHALLO
      USE YOWSHAL  , ONLY : TFAK     ,INDEP

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IJ,M,K,IJS,IJL
      REAL :: DELT25, DELT2, DEL2
      REAL :: F(NINF-1:NSUP,NANG,NFRE)
      REAL,DIMENSION(IJS:IJL) :: TEMP1, TEMP2, EM, FM
      LOGICAL :: LLAK

! ----------------------------------------------------------------------

!*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
!        ------------------------------------------------

      DO IJ=IJS,IJL
        EM(IJ) = EPSMIN
        FM(IJ) = EPSMIN
      ENDDO

      DELT25 = WETAIL*FR(NFRE)*DELTH
      DELT2 = FRTAIL*DELTH

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      IF (ISHALLO.EQ.1 .OR. .NOT.LLAK) THEN

!*    2.1 DEEP WATER INTEGRATION.
!         -----------------------

        DO M=1,NFRE
          K=1
          DO IJ=IJS,IJL
            TEMP2(IJ) = F(IJ,K,M)
          ENDDO
          DO K=2,NANG
            DO IJ=IJS,IJL
              TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
            ENDDO
          ENDDO
          DO IJ=IJS,IJL
            EM(IJ) = EM(IJ)+TEMP2(IJ)*DFIM(M)
            FM(IJ) = FM(IJ)+DFIMOFR(M)*TEMP2(IJ)
          ENDDO
        ENDDO

      ELSE

!*    2.2 SHALLOW WATER INTEGRATION.
!         --------------------------

        DO M=1,NFRE
          K=1
          DO IJ=IJS,IJL
            TEMP1(IJ) = DFIM(M)/SQRT(TFAK(INDEP(IJ),M))
            TEMP2(IJ) = F(IJ,K,M) 
          ENDDO
          DO K=2,NANG
            DO IJ=IJS,IJL
              TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
            ENDDO
          ENDDO
          DO IJ=IJS,IJL
            EM(IJ) = EM(IJ)+TEMP2(IJ)*DFIM(M)
            FM(IJ) = FM(IJ)+DFIMOFR(M)*TEMP2(IJ)
          ENDDO
        ENDDO

      ENDIF

!*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
!*       NORMALIZE WITH TOTAL ENERGY.
!        ------------------------------------------

      IF (ISHALLO.EQ.1 .OR. .NOT.LLAK) THEN

        DO IJ=IJS,IJL
          EM(IJ) = EM(IJ)+DELT25*TEMP2(IJ)
          FM(IJ) = FM(IJ)+DELT2*TEMP2(IJ)
          FM(IJ) = EM(IJ)/FM(IJ)
        ENDDO

      ELSE

        DEL2 = DELT2*SQRT(G)/ZPI
        DO IJ=IJS,IJL
          EM(IJ) = EM(IJ)+DELT25*TEMP2(IJ)
          FM(IJ) = FM(IJ)+DELT2*TEMP2(IJ)
          FM(IJ) = EM(IJ)/FM(IJ)
        ENDDO

      ENDIF

      RETURN
      END SUBROUTINE FEMEAN
