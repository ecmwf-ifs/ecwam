      SUBROUTINE SPR (NANG, THETAQ, THETA, ST)

! ----------------------------------------------------------------------

!**** *SPR* - ROUTINE TO COMPUTE SPREADING FACTOR.

!     SUSANNE HASSELMANN  JULY 1986.

!*    PURPOSE.
!     --------

!       COMPUTATION OF COS**2 SPREADING FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SPR (NANG, THETAQ, THETA, ST)*
!          *NANG*    INTEGER  DIMENSION FOR ANGULAR INTERVALS
!          *THETAQ*  REAL     MEAN WAVE DIRECTION
!          *THETA*   REAL     ANGULAR DIRECTIONS
!          *ST*      REAL     SPREADING FUNCTION

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.


!     REFERENCES.
!     -----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWPCONS , ONLY : PI

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NANG
      REAL, INTENT(IN) :: THETAQ
      REAL, DIMENSION(NANG), INTENT(IN) :: THETA
      REAL, DIMENSION(NANG), INTENT(OUT) :: ST

      INTEGER :: K

      REAL :: ZDP, THE

! ----------------------------------------------------------------------

      ZDP=2./PI

!     SPREAD FCT. WITH COS**2.

      DO K=1,NANG
        THE = COS(THETA(K)-THETAQ)
        IF (THE.GT.0.) THEN
          ST(K) = ZDP*THE**2
          IF(ST(K).LT.0.1E-08) ST(K)=0.
        ELSE
          ST(K) = 0.
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE SPR
