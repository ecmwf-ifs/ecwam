      SUBROUTINE DOTDC (IJS, IJL, SHLFAC)

! ----------------------------------------------------------------------

!**** *DOTDC* - SIGMA DOT FOR CURRENT REFRACTION.

!     H. GUNTHER   GKSS/ECMWF   17/02/91
!     J. BIDLOT    ECMWF FEBRUARY 1997 : MESSAGE PASSING
!     J. BIDLOT    ECMWF JANUARY 2004 : REMOVE READING FROM FILE.
!                                       THE SIGMA DOT ARE NOW IN SDOT.

!*    PURPOSE.
!     --------

!       SCATTER SIGMA/ SINH (2*K*D) TABLE.

!**   INTERFACE.
!     ----------

!       *CALL* *DOTDC (IJS, IJL, SHLFAC)*
!          *IJS*     - INDEX OF FIRST GRID POINT.
!          *IJL*     - INDEX OF LAST  GRID POINT.
!          *SHLFAC * - SIGMA/SINH(2KD)

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NFRE_RED
      USE YOWSHAL  , ONLY : TSIHKD   ,INDEP

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM) :: IJ, M
      REAL(KIND=JWRB),DIMENSION(IJS:IJL,NFRE_RED), INTENT(OUT):: SHLFAC

! ----------------------------------------------------------------------

!*    3. GATHER SIGMA /SINH(2KD) FROM TABLE.
!        -----------------------------------

      DO M=1,NFRE_RED
        DO IJ=IJS,IJL
          SHLFAC(IJ,M) = TSIHKD(INDEP(IJ),M)
        ENDDO
      ENDDO

      END SUBROUTINE DOTDC
