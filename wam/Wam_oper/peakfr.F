      SUBROUTINE PEAKFR (F, IJS, IJL, PEAKF)

! ----------------------------------------------------------------------

!**** *PEAKFR* - COMPUTATION OF PEAK FREQUENCY AT EACH GRID POINT.

!     H. GUNTHER      ECMWF            DECEMBER 1989
!     (CODE REMOVED FROM SUB. FEMEAN)

!*    PURPOSE.
!     --------

!       COMPUTE PEAK FREQUENCY AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *PEAKFR (F, IJS, IJL, PEAKF)*
!          *F*     - SPECTRUM.
!          *IJS*   - INDEX OF FIRST GRIDPOINT
!          *IJL*   - INDEX OF LAST GRIDPOINT
!          *PEAKF* - PEAKFREQUENCY

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

      USE YOWFRED  , ONLY : FR
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NIBLO

! ----------------------------------------------------------------------

      DIMENSION TEMP1(NINF:NSUP), TEMP2(NINF:NSUP)
      DIMENSION F(NINF-1:NSUP,NANG,NFRE), PEAKF(NIBLO)

!*    1. INITIALIZE ARRAYS
!        -----------------

      DO IJ = IJS,IJL
        TEMP1(IJ) = 0.
        PEAKF(IJ) = 0.
      ENDDO

!*    2. LOOP OVER FREQUENCIES
!        ---------------------

      DO M = 1,NFRE

!*    2.1 COMPUTE 1-D SPECTRUM (WITHOUT DELTA THETA)
!        -------------------------------------------

        DO IJ = IJS,IJL
          TEMP2(IJ) = 0.
        ENDDO
        DO K = 1,NANG
          DO IJ = IJS,IJL
            TEMP2(IJ) = TEMP2(IJ) + F(IJ,K,M)
          ENDDO
        ENDDO

!*    2.2 DEFINE PEAK FREQUENCY
!         ---------------------

        DO IJ = IJS,IJL
          IF (TEMP1(IJ).LT.TEMP2(IJ)) THEN
            TEMP1(IJ) = TEMP2(IJ)
            PEAKF(IJ) = FR(M)
          ENDIF
        ENDDO

      ENDDO

      RETURN
      END SUBROUTINE PEAKFR
