! ======================================================================
 
      SUBROUTINE WAM_SORTINI (IARRAY, INDEX, N)
 
! ======================================================================
 
!**** *WAM_SORTINI* - SORTS AN INDEX ARRAY FORM AN INTEGER DATA ARRAY.
 
!     HEINZ GUNTHER    ECMWF/GKSS    JANUARY 1991
 
!     PURPOSE
!     -------
 
!        FIND THE INDEX ARRAY WHICH SORTS THE DATA.
 
!     METHOD
!     ------
 
!       THE INPUT ARRAY IS SCANNED AND AN INDEX ARRAY IS GENERATED,
!       WHICH CONTAINS THE INDEX OF THE INPUT ARRAY TO SORT THE DATA
!       IN INCREASING ORDER. THE INPUT ARRAY IS NOT CHANGED.
 
!**   INTERFACE
!     ---------
 
!        *CALL* *WAM_SORTINI (IARRAY, INDEX, N)
 
!          IARRAY  INTEGER   INPUT DATA TO BE SCANNED
!          INDEX   INTEGER   OUTPUT ARRAY OF INDICES.
!          N       INTEGER   NUMBER OF DATA TO BE SCANNED
 
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
      INTEGER(KIND=JWIM), DIMENSION(N), INTENT(IN) :: IARRAY 
      INTEGER(KIND=JWIM), DIMENSION(N), INTENT(OUT) :: INDEX

      INTEGER(KIND=JWIM) :: I, J, IH, IHM

      LOGICAL :: LLEXIT

! ----------------------------------------------------------------------
 
!*   1. INITIAL INDEX ARRAY.
!       -------------------
 
      DO I=1,N
        INDEX(I) = I
      ENDDO
 
! ----------------------------------------------------------------------
 
!*   2. CHECK FIRST TWO DATA.
!       ---------------------
 
      IF (N.LT.2) RETURN
      IF (IARRAY(2).LT.IARRAY(1)) THEN
         INDEX(2)=1
         INDEX(1)=2
      ENDIF
      IF (N.EQ.2) RETURN
 
! ----------------------------------------------------------------------
 
!*   3. DO REST OF DATA.
!       ----------------
 
      DO I=3,N
         IH = INDEX(I)
         IHM = INDEX(I-1)
         IF (IARRAY(IH).LT.IARRAY(IHM)) THEN
            LLEXIT=.FALSE. 
            DO J=I-1,2,-1
              INDEX(J+1) = INDEX(J)
              IF (IARRAY(INDEX(J-1)).LE.IARRAY(IH)) THEN
                 INDEX(J) = IH
                 LLEXIT=.TRUE. 
                 EXIT
              ENDIF
            ENDDO
            IF(.NOT.LLEXIT) THEN
              INDEX(2) = INDEX(1)
              INDEX(1) = IH
            ENDIF
         ENDIF
      ENDDO

      END SUBROUTINE WAM_SORTINI
