      SUBROUTINE  MNINTW (                                              &
!                     DIMENSIONS
     &                  NSPECW, MPARTSW,                                &
!                     INPUT
     &                  NPARTW, NIPARTINFW,                             &
!                     OUTPUT
     &                  NINTW)

!---------------------------------------------------------------------

!**** *MNINTW*  READS GRID OUTPUT FROM PREPROC.

!     R.BROKOPF MAY 94.

!*    PURPOSE.
!     --------
!       COMPUTES INDEX OF WIND-SEA PARTITIONING FOR EACH SPECTRUM.

!**   INTERFACE.
!     ----------
!       *CALL* *MNINTW(
!                     DIMENSIONS
!     &                  NSPECW , NSPECW , MPARTSW ,
!                     INPUT
!     &                  NPARTW, NIPARTINFW,
!                     OUTPUT
!     &                  NINTW)*

!     METHOD.
!     -------

!       NONE

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: NSPECW , MPARTSW

      INTEGER(KIND=JWIM) :: NPARTW(NSPECW) , NINTW(NSPECW),             &
     &                      NIPARTINFW(NSPECW,MPARTSW)

      INTEGER(KIND=JWIM) :: IPART, ISPEC

! ----------------------------------------------------------------------

      DO ISPEC = 1, NSPECW
         NINTW(ISPEC) = 0
      END DO

      DO ISPEC = 1, NSPECW
        DO IPART = 1,NPARTW(ISPEC)
          IF(NIPARTINFW(ISPEC,IPART).NE.0)                              &
     &       NINTW(ISPEC)=IPART
        ENDDO
      ENDDO

      END SUBROUTINE MNINTW
