! ======================================================================

      SUBROUTINE WAM_SORTI (IARRAY, INDEX, N)
!
! ======================================================================
 
!**** *WAM_SORTI* - SORTS AN INTEGER ARRAY.
 
!     HEINZ GUNTHER    ECMWF/GKSS    JANUARY 1991
 
!     PURPOSE
!     -------
 
!        SORT AN INTEGER ARRAY BY AN INDEX ARRAY.
 
!     METHOD
!     ------
 
!       THE DATA IN THE INPUT ARRAY ARE SORTED AS DEFINED IN THE INPUT
!       INDEX ARRAY.
!       THIS ARRAY MAY BE GENEGATED BY SUB. WAM_SORTINI.
 
!**   INTERFACE
!     ---------
 
!        *CALL* *WAM_SORTI (IARRAY, INDEX, N)
 
!          IARRAY  INTEGER   INPUT/ OUTPUT DATA ARRAY.
!          INDEX   INTEGER   INPUT INDEX ARRAY.
!          N       INTEGER   NUMBER OF DATA TO BE SORTED.
 
!     EXTERNALS
!     ---------
 
!          NONE
 
!     REFERENCES
!     ----------
 
!          NONE
 
! ----------------------------------------------------------------------

     USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: N
      INTEGER(KIND=JWIM), DIMENSION(N), INTENT(INOUT) :: IARRAY 
      INTEGER(KIND=JWIM), DIMENSION(N), INTENT(IN) :: INDEX

      INTEGER(KIND=JWIM) :: I
      INTEGER(KIND=JWIM), DIMENSION(N) :: IWORK

! ----------------------------------------------------------------------
 
!*   1. SORT DATA INTO WORK ARRAY.
!       --------------------------
 
      DO I=1,N
        IWORK(I) = IARRAY(INDEX(I))
      ENDDO

! ----------------------------------------------------------------------
 
!*   2. COPY BACK TO DATA ARRAY.
!       ------------------------
 
      DO I=1,N
        IARRAY(I) = IWORK(I)
      ENDDO

      END SUBROUTINE WAM_SORTI
