      SUBROUTINE OUTUBUF (IU08, IU18, IFORM)

! ----------------------------------------------------------------------

!**** *OUTUBUF* - ROUTINE TO WRITE COMMON UBUF TO DISK.

!     H.GUNTHER            ECMWF       04/04/1990

!*    PURPOSE.
!     -------

!       TO WRITE OUT THE COMPUTED COMMON UBUF.

!**   INTERFACE.
!     ----------

!       *CALL* *OUTUBUF (IU08, IU18, IFORM)*
!         *IU08*   - LOGICAL UNIT FOR  UNFORMATED WRITE
!         *IU18*   - LOGICAL UNIT FOR    FORMATED WRITE
!         *IFORM*  - FORMAT OPTION  = 1  UNFORMATED WRITE
!                                   = 2  FORMATED WRITE
!                                   OTHERWISE BOTH
!     METHOD.
!     -------

!       FORMATED AND OR UNFORMATED WRITE TO UNITS.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,KLAT2    ,WADV

! ----------------------------------------------------------------------

  998 FORMAT(10I8)
  999 FORMAT(10E6.3)

! ----------------------------------------------------------------------

!*    1. WRITE COMMON UBUF UNFORMATED.
!        -----------------------------

      IF (IFORM.NE.2) WRITE (IU08) KLAT, KLON
      IF (IFORM.NE.2) WRITE (IU08) KLAT2
      IF (IFORM.NE.2) WRITE (IU08) WADV

! ----------------------------------------------------------------------

!*    2. WRITE COMMON UBUF FORMATED.
!        ---------------------------

      IF (IFORM.NE.1) WRITE (IU18,998) KLAT, KLON
      IF (IFORM.NE.1) WRITE (IU18,998) KLAT2
      IF (IFORM.NE.1) WRITE (IU18,999) WADV
      RETURN
      END SUBROUTINE OUTUBUF
