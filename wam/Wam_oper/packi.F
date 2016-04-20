      SUBROUTINE PACKI (NIN, NOUT, NDIM, IFLAG, IPACK)

! ----------------------------------------------------------------------

!**** *PACKI* - PACKS A INTEGER ARRAY.

!     R. PORTZ     MPI         15/01/1991

!*    PURPOSE.
!     --------

!       TO REMOVE FLAGGED POINTS FROM AN INTEGER ARRAY.

!**   INTERFACE.
!     ----------

!       *CALL* *PACKI (NIN, NOUT, NDIM, IFLAG, IPACK)*
!          *NIN*   -  NUMBER OP POINTS IN ARRAYS BEFORE PACKING.
!          *NOUT*  -  NUMBER OF POINTS IN ARRAYES AFTER PACKING.
!          *NDIM*  -  DIMENSION OF ARRAYS.
!          *IFLAG* -  FLAG ARRAY.
!          *IPACK* -  ARRAY TO BE PACKED / PACKED ARRAY AT OUTPUT.

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

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NIN, NDIM
      INTEGER, INTENT(OUT) :: NOUT
      INTEGER, DIMENSION(NDIM), INTENT(INOUT) :: IFLAG, IPACK


      INTEGER :: K, I

! ----------------------------------------------------------------------


!*    1. LOOP OVER POINTS.
!        -----------------

      K = 0
      DO I = 1, NIN
        IF (IFLAG(I).EQ.0) THEN
          K = K + 1
        ELSE
          IPACK(I-K) = IPACK(I)
        ENDIF
      ENDDO

!*    2. NUMBER OF POINTS LEFT.
!        ----------------------

      NOUT = NIN - K

      RETURN
      END SUBROUTINE PACKI
