      SUBROUTINE FINDBICE (NDIM, NBOUN, BLATB, BLNGB, IGARB, IJARB,
     1                     NSTART,NEND)
C ----------------------------------------------------------------------
C
C**** *FINDBICE* - FIND BLOCK AND GRID POINT NUMBER.
C
C     R. PORTZ     MPI         15/01/1991
C     J. BIDLOT    ECMWF       JUNE 1996  MESSAGE PASSING
C     ORIGINALLY CALLED FINDB BUT MODIFIED TO ACCOMODATE MPP VERSION
C
C*    PURPOSE.
C     -------
C
C       FIND BLOCK AND GRID POINT NUMBER FOR A GIVEN ARRAY
C       OF LONGITUDES AND LATITUDES.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *FINDBICE (NDIM, NBOUN, BLATB, BLNGB, IGARB, IJARB)*
C                         NSTART,NEND)
C          *NDIM*    - DIMENSION OF ARRAYS.
C          *NBOUN*   - NUMBER OF POINTS IN ARRAYS.
C          *BLATB*   - INPUT LATITUDES.
C          *BLNGB*   - INPUT LONGITUDES.
C          *IGARB*   - OUTPUT BLOCK NUMBERS.
C          *IJARB*   - OUTPUT SEA POINT NUMBERS.
C  MESSAGE PASSING SPECIFIC
C      *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
C      *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
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
C     THIS IS NOW A DUMMY ROUTINE
C
      RETURN
      END
