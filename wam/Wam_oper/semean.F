      SUBROUTINE SEMEAN (F3, IJS, IJL)

! ----------------------------------------------------------------------

!**** *SEMEAN* - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.

!     S.D. HASSELMANN.
!     OPTIMIZED BY: L. ZAMBRESKY AND H. GUENTHER

!*    PURPOSE.
!     --------

!       TO COMPUTE TOTAL ENERGY AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *SEMEAN(F3, IJS, IJL)*
!          *F3*  - SPECTRUM.
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT

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

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,WETAIL
      USE YOWMEAN  , ONLY : EMEAN
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSMIN

! ----------------------------------------------------------------------

      DIMENSION F3(NINF-1:NSUP,NANG,NFRE), TEMP(IJS:IJL)

! ----------------------------------------------------------------------

!*    1. INITIALISE ENERGY ARRAY.
!        ------------------------

      DO IJ=IJS,IJL
        EMEAN(IJ) = EPSMIN
      ENDDO

! ----------------------------------------------------------------------

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTION.
!        -----------------------------------------

      DO M=1,NFRE
        DO IJ=IJS,IJL
          TEMP(IJ) = 0.
        ENDDO
        DO K=1,NANG
          DO IJ=IJS,IJL
            TEMP(IJ) = TEMP(IJ)+F3(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=IJS,IJL
          EMEAN(IJ) = EMEAN(IJ)+DFIM(M)*TEMP(IJ)
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    3. ADD TAIL ENERGY.
!        ----------------

      DELT25 = WETAIL*FR(NFRE)*DELTH
      DO IJ=IJS,IJL
        EMEAN(IJ) = EMEAN(IJ)+DELT25*TEMP(IJ)
      ENDDO

      RETURN
      END SUBROUTINE SEMEAN
