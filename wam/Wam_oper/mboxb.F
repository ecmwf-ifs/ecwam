      SUBROUTINE MBOXB (NBOUN, AMOWEB, AMOSOB, AMOEAB,
     &                  AMONOB, BLATB, BLNGB)

! ----------------------------------------------------------------------

!**** *MBOXB* - MAKE BOX OF FINE GRID IN COARSE GRID.

!     R. PORTZ     MPI         15/01/1991

!*    PURPOSE.
!     -------

!       COMPUTE LATITUDES AND LONGITUDES OF COARSE GRID POINTS
!       AT NEST BOUNDARY.

!**   INTERFACE.
!     ----------

!       *CALL* *MBOXB (NBOUN, AMOWEB, AMOSOB, AMOEAB, AMONOB,
!                      BLATB, BLNGB)*
!          *NBOUN*  - NUMBER OF BOUNDARY POINTS.
!          *AMOWEB* - WESTERN  LONGITUDE OF COARSE GRID.
!          *AMOSOB* - SOUTHERN LATITUDE  OF COARSE GRID.
!          *AMOEAB* - EASTERN  LONGITUDE OF COARSE GRID.
!          *AMONOB* - NORTHERN LATITUDE  OF COARSE GRID.
!          *BLATB*  - LATITUDE  OF BOUNDARY POINTS.
!          *BLNGB*  - LONGITUDE OF BOUNDARY POINTS.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWMAP   , ONLY : XDELLA   ,XDELLO
      USE YOWPRPROC, ONLY : NBMAX
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER, INTENT(OUT) :: NBOUN
      REAL, INTENT(IN) :: AMOWEB, AMOSOB, AMOEAB, AMONOB
      REAL, DIMENSION(NBMAX), INTENT(INOUT) :: BLNGB, BLATB

      INTEGER :: I, K
      INTEGER :: NLNGB, NLATB

! ----------------------------------------------------------------------

!*    1. COMPUTED THE SQUARE BOX
!     --------------------------

      NLNGB = NINT((AMOEAB - AMOWEB) / XDELLO) + 1

      NLATB = (NINT((AMONOB - AMOSOB) / XDELLA) - 1) * 2

      NBOUN = (NLNGB * 2) + NLATB

!     DIMENSION CHECK: NBMAX > NBOUN

      IF (NBMAX .LT. NBOUN) THEN
        WRITE(IU06,*) ' **********************************'
        WRITE(IU06,*) ' *                                *'
        WRITE(IU06,*) ' *   FATAL ERROR IN SUB. MBOX     *'
        WRITE(IU06,*) ' *   ========================     *'
        WRITE(IU06,*) ' *  NUMBER OF BOUNDARY POINTS     *'
        WRITE(IU06,*) ' *  EXCEEDS DIMENSION.            *'
        WRITE(IU06,*) ' *  DIMENSION IS        NBMAX = ', NBMAX
        WRITE(IU06,*) ' *  NUMBER OF POINTS IS NBOUN = ', NBOUN
        WRITE(IU06,*) ' *                                *'
        WRITE(IU06,*) ' **********************************'
        CALL ABORT1
      ENDIF

!     COMPUTES THE BOUNDARY POINTS FOR THE FIRST AND THE LAST
!     LATITUDE

      K = NLATB + NLNGB
      DO I = 1, NLNGB
        BLATB(I) = AMOSOB
        BLNGB(I) = AMOWEB + REAL(I-1) * XDELLO
        BLATB(K+I) = AMONOB
        BLNGB(K+I) = BLNGB(I)
      ENDDO

!     COMPUTED THE EAST AND THE WEST BOUNDARY POINT FOR
!     EACH LATITUDE

      K = 0
      DO I = 2,NLATB,2
        K = K + 1
        BLATB(NLNGB+I-1) = AMOSOB + REAL(K) * XDELLA
        BLNGB(NLNGB+I-1) = AMOWEB
        BLATB(NLNGB+I)   = BLATB(NLNGB+I-1)
        BLNGB(NLNGB+I)   = AMOEAB
      ENDDO

      RETURN
      END SUBROUTINE MBOXB
