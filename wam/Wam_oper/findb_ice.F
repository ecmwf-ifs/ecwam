      SUBROUTINE FINDBICE (NDIM, NBOUN, BLATB, BLNGB, IGARB, IJARB,
     &                     NSTART,NEND)
! ----------------------------------------------------------------------

!**** *FINDBICE* - FIND BLOCK AND GRID POINT NUMBER.

!     R. PORTZ     MPI         15/01/1991
!     J. BIDLOT    ECMWF       JUNE 1996  MESSAGE PASSING
!     ORIGINALLY CALLED FINDB BUT MODIFIED TO ACCOMODATE MPP VERSION

!*    PURPOSE.
!     -------

!       FIND BLOCK AND GRID POINT NUMBER FOR A GIVEN ARRAY
!       OF LONGITUDES AND LATITUDES.

!**   INTERFACE.
!     ----------

!       *CALL* *FINDBICE (NDIM, NBOUN, BLATB, BLNGB, IGARB, IJARB)*
!                         NSTART,NEND)
!          *NDIM*    - DIMENSION OF ARRAYS.
!          *NBOUN*   - NUMBER OF POINTS IN ARRAYS.
!          *BLATB*   - INPUT LATITUDES.
!          *BLNGB*   - INPUT LONGITUDES.
!          *IGARB*   - OUTPUT BLOCK NUMBERS.
!          *IJARB*   - OUTPUT SEA POINT NUMBERS.
!  MESSAGE PASSING SPECIFIC
!      *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!      *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN

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

!     THIS IS NOW A DUMMY ROUTINE

      RETURN
      END SUBROUTINE FINDBICE






