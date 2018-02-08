      SUBROUTINE WAMNORM(FIELD,NX,NY,WNORM)

! ----------------------------------------------------------------------

!**** *WAMNORM* - COMPUTES THE AVERAGE, MINIMUM AND MAXIMUM OF FIELD 

!     J. BIDLOT       ECMWF   MAY 1998 

!*    PURPOSE.
!     --------

!     COMPUTES THE AVERAGE, MINIMUM AND MAXIMUM OF FIELD FOR ALL VALUES
!     DIFFERENT THAN ZMISS 

!**   INTERFACE.
!     ----------
!     *CALL*  WAMNORM(FIELD,NX,NY,WNORM)
!          *FIELD*     REAL    INPUT FIELD
!          *NX*        INTEGER FIRST DIMENSION OF FIELD
!          *NY*        INTEGER SECOND DIMENSION OF FIELD
!          *WNORM*     REAL    WAM NORM

!     METHOD.
!     -------

!     COMPUTES THE AVERAGE, MINIMUM AND MAXIMUM OF FIELD FOR ALL VALUES
!     DIFFERENT THAN ZMISS.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWPCONS , ONLY : ZMISS

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NX, NY

      REAL, INTENT(IN) :: FIELD(NX,NY)
      REAL, INTENT(OUT):: WNORM(4) ! average, minmum, maximum, number of points

      INTEGER :: I, J
      REAL :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('WAMNORM',0,ZHOOK_HANDLE)

      WNORM(1)=0.
      WNORM(2)=FIELD(1,1)
      WNORM(3)=FIELD(1,1)
      WNORM(4)=0.
      DO J=1,NY
        DO I=1,NX
          IF(FIELD(I,J).NE.ZMISS) THEN
            WNORM(1)=WNORM(1)+FIELD(I,J)
            WNORM(2)=MIN(WNORM(2),FIELD(I,J))
            WNORM(3)=MAX(WNORM(3),FIELD(I,J))
            WNORM(4)=WNORM(4)+1
          ENDIF
        ENDDO
      ENDDO
      WNORM(1)=WNORM(1)/MAX(WNORM(4),1.)

      IF (LHOOK) CALL DR_HOOK('WAMNORM',1,ZHOOK_HANDLE)

      END SUBROUTINE WAMNORM
