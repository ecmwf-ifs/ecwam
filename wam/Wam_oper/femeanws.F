      SUBROUTINE FEMEANWS (F, IJS, IJL, EM, FM, LLWS)

! ----------------------------------------------------------------------

!**** *FEMEANWS* - COMPUTATION OF MEAN ENERGY, MEAN FREQUENCY 
!              !!! FOR PART OF THE SPECTRUM WITH LLWS TRUE.


!*    PURPOSE.
!     --------

!       COMPUTE MEAN FREQUENCY AT EACH GRID POINT FOR PART OF THE
!       SPECTRUM WHERE LLWS IS TRUE.

!**   INTERFACE.
!     ----------

!       *CALL* *FEMEANWS (F, IJS, IJL, EM, FM, LLWS)*
!              *F*   - SPECTRUM.
!              *IJS* - INDEX OF FIRST GRIDPOINT
!              *IJL* - INDEX OF LAST GRIDPOINT
!              *EM*  - MEAN WAVE ENERGY (INPUT)
!              *FM*  - MEAN WAVE FREQUENCY (OUTPUT)
!              *LLWS*- TRUE WHERE SINPUT IS POSITIVE.

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

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IJ,M,K,IJS,IJL
      REAL :: DELT25, DELT2
      REAL :: F(NINF-1:NSUP,NANG,NFRE)
      REAL,DIMENSION(IJS:IJL) :: TEMP2, EM, FM
      LOGICAL, DIMENSION(IJS:IJL,NANG,NFRE) :: LLWS

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
      
      DO M=1,NFRE
        K=1
        DO IJ=IJS,IJL
          IF(LLWS(IJ,K,M)) THEN
            TEMP2(IJ) = F(IJ,K,M)
          ELSE
            TEMP2(IJ) = 0.
          ENDIF
        ENDDO
        DO K=2,NANG
          DO IJ=IJS,IJL
            IF(LLWS(IJ,K,M)) TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=IJS,IJL
          EM(IJ) = EM(IJ)+TEMP2(IJ)*DFIM(M)
          FM(IJ) = FM(IJ)+DFIMOFR(M)*TEMP2(IJ)
        ENDDO
      ENDDO

!*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
!*       NORMALIZE WITH TOTAL ENERGY.
!        ------------------------------------------

      DO IJ=IJS,IJL
        EM(IJ) = EM(IJ)+DELT25*TEMP2(IJ)
        FM(IJ) = FM(IJ)+DELT2*TEMP2(IJ)
        FM(IJ) = EM(IJ)/FM(IJ)
      ENDDO

    

      RETURN
      END SUBROUTINE FEMEANWS
