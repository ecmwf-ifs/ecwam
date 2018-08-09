      SUBROUTINE PACKR (NIN, NOUT, NDIM, IFLAG, RPACK)

! ----------------------------------------------------------------------

!**** *PACKR* - PACKS A REAL ARRAY.

!     R. PORTZ     MPI         15/01/1991

!*    PURPOSE.
!     -------

!       TO REMOVE FLAGGED POINTS FROM AN REAL ARRAY.

!**   INTERFACE.
!     ----------

!       *CALL* *PACKR (NIN, NOUT, NDIM, IFLAG, RPACK)*
!          *NIN*   -  NUMBER OP POINTS IN ARRAYS BEFORE PACKING.
!          *NOUT*  -  NUMBER OF POINTS IN ARRAYES AFTER PACKING.
!          *NDIM*  -  DIMENSION OF ARRAYS.
!          *IFLAG* -  FLAG ARRAY.
!          *RPACK* -  ARRAY TO BE PACKED / PACKED ARRAY AT OUTPUT.

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

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: NIN, NDIM
      INTEGER(KIND=JWIM), INTENT(OUT) :: NOUT
      INTEGER(KIND=JWIM), DIMENSION(NDIM), INTENT(INOUT) :: IFLAG
      REAL(KIND=JWRB), DIMENSION(NDIM), INTENT(INOUT) :: RPACK


      INTEGER(KIND=JWIM) :: K, I

! ----------------------------------------------------------------------


!*    1. LOOP OVER POINTS.
!        -----------------

      K = 0
      DO I = 1, NIN
        IF (IFLAG(I).EQ.0) THEN
          K = K + 1
        ELSE
          RPACK(I-K) = RPACK(I)
        ENDIF
      ENDDO

!*    2. NUMBER OF POINTS LEFT.
!        ----------------------

      NOUT = NIN - K

      END SUBROUTINE PACKR
