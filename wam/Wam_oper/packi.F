      SUBROUTINE PACKI (NIN, NOUT, NDIM, IFLAG, IPACK)

C ----------------------------------------------------------------------
C
C**** *PACKI* - PACKS A INTEGER ARRAY.
C
C     R. PORTZ     MPI         15/01/1991
C
C*    PURPOSE.
C     --------
C
C       TO REMOVE FLAGGED POINTS FROM AN INTEGER ARRAY.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *PACKI (NIN, NOUT, NDIM, IFLAG, IPACK)*
C          *NIN*   -  NUMBER OP POINTS IN ARRAYS BEFORE PACKING.
C          *NOUT*  -  NUMBER OF POINTS IN ARRAYES AFTER PACKING.
C          *NDIM*  -  DIMENSION OF ARRAYS.
C          *IFLAG* -  FLAG ARRAY.
C          *IPACK* -  ARRAY TO BE PACKED / PACKED ARRAY AT OUTPUT.
C
C     METHOD.
C     -------
C
C       NONE.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
      DIMENSION IFLAG(NDIM), IPACK(NDIM)
C
C ----------------------------------------------------------------------
C
C
C*    1. LOOP OVER POINTS.
C        -----------------
C
 1000 CONTINUE
      K = 0
      DO 1001 I = 1, NIN
         IF (IFLAG(I).EQ.0) THEN
            K = K + 1
         ELSE
            IPACK(I-K) = IPACK(I)
         ENDIF
 1001 CONTINUE
C
C*    2. NUMBER OF POINTS LEFT.
C        ----------------------
C
      NOUT = NIN - K

      RETURN
      END
