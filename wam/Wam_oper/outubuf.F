      SUBROUTINE OUTUBUF (IU08, IU18, IFORM,NIBLO,ML,KOBSLAT,KOBSLON)

! ----------------------------------------------------------------------

!**** *OUTUBUF* - ROUTINE TO WRITE COMMON UBUF TO DISK.
!     H.GUNTHER            ECMWF       04/04/1990

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     J. BIDLOT            ECMWF       11/2003
!                          IF YOUR ARE RUNNING AT ECMWF:
!                          BE AWARE THAT IF YOU CHANGE ANYTHING TO THE
!                          STRUCTURE OF THE OUTPUT FILE YOU WILL HAVE TO
!                          MAKE SURE THAT IT IS CREATED FOR YOUR RUN,
!                          OTHERWISE IT MIGHT PICK UP THE DEFAULT ONE
!                          THAT IS ALREADY ON DISK.
!                          YOU ALSO HAVE TO CHANGE MPDECOMP.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

      INTEGER, DIMENSION(NIBLO,ML,2) :: KOBSLON, KOBSLAT

  998 FORMAT(10I8)
  999 FORMAT(10E6.3)

! ----------------------------------------------------------------------

!*    1. WRITE COMMON UBUF UNFORMATED.
!        -----------------------------

      IF (IFORM.NE.2) WRITE (IU08) KLAT, KLON
      IF (IFORM.NE.2) WRITE (IU08) KLAT2
      IF (IFORM.NE.2) WRITE (IU08) WADV
      IF (IFORM.NE.2) WRITE (IU08) KOBSLAT
      IF (IFORM.NE.2) WRITE (IU08) KOBSLON

! ----------------------------------------------------------------------

!*    2. WRITE COMMON UBUF FORMATED.
!        ---------------------------

      IF (IFORM.NE.1) WRITE (IU18,998) KLAT, KLON
      IF (IFORM.NE.1) WRITE (IU18,998) KLAT2
      IF (IFORM.NE.1) WRITE (IU18,999) WADV
      IF (IFORM.NE.1) WRITE (IU18,998) KOBSLAT
      IF (IFORM.NE.1) WRITE (IU18,998) KOBSLON

      RETURN
      END SUBROUTINE OUTUBUF
