      SUBROUTINE PACKR (NIN, NOUT, NDIM, IFLAG, PACK)

! ----------------------------------------------------------------------

!**** *PACKR* - PACKS A REAL ARRAY.

!     R. PORTZ     MPI         15/01/1991

!*    PURPOSE.
!     -------

!       TO REMOVE FLAGGED POINTS FROM AN REAL ARRAY.

!**   INTERFACE.
!     ----------

!       *CALL* *PACKR (NIN, NOUT, NDIM, IFLAG, PACK)*
!          *NIN*   -  NUMBER OP POINTS IN ARRAYS BEFORE PACKING.
!          *NOUT*  -  NUMBER OF POINTS IN ARRAYES AFTER PACKING.
!          *NDIM*  -  DIMENSION OF ARRAYS.
!          *IFLAG* -  FLAG ARRAY.
!          *PACK*  -  ARRAY TO BE PACKED / PACKED ARRAY AT OUTPUT.

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

      DIMENSION IFLAG(NDIM), PACK(NDIM)

! ----------------------------------------------------------------------


!*    1. LOOP OVER POINTS.
!        -----------------

      K = 0
      DO I = 1, NIN
        IF (IFLAG(I).EQ.0) THEN
          K = K + 1
        ELSE
          PACK(I-K) = PACK(I)
        ENDIF
      ENDDO

!*    2. NUMBER OF POINTS LEFT.
!        ----------------------

      NOUT = NIN - K

      RETURN
      END SUBROUTINE PACKR




