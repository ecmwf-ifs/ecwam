      SUBROUTINE OUTUBUF (IU08, NIBLO, ML,
     &                    KOBSLAT,KOBSLON,KOBSCOR,KOBSRLAT,KOBSRLON)

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

!       *CALL* *OUTUBUF (IU08)*
!         *IU08*   - LOGICAL UNITS FOR  UNFORMATED WRITE
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

      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,KCOR     ,
     &                      KRLAT    ,KRLON    ,
     &                      WLAT     ,WCOR     ,WRLAT    ,WRLON

! ----------------------------------------------------------------------

      INTEGER :: IU08(0:2)
      INTEGER, DIMENSION(NIBLO,ML,2) :: KOBSLON, KOBSLAT
      INTEGER, DIMENSION(NIBLO,ML,2) :: KOBSRLON, KOBSRLAT
      INTEGER, DIMENSION(NIBLO,ML,4) :: KOBSCOR

! ----------------------------------------------------------------------

!*    1. WRITE COMMON UBUF UNFORMATED.
!        -----------------------------

      DO IP=0,2
        WRITE (IU08(IP)) KLAT
        WRITE (IU08(IP)) KLON
      ENDDO

      WRITE (IU08(1)) KRLAT
      WRITE (IU08(1)) KRLON

      WRITE (IU08(2)) KCOR

      DO IP=0,2
        WRITE (IU08(IP)) WLAT 
      ENDDO

      WRITE (IU08(1)) WRLAT
      WRITE (IU08(1)) WRLON

      WRITE (IU08(2)) WCOR 

      DO IP=0,2
        WRITE (IU08(IP)) KOBSLAT
        WRITE (IU08(IP)) KOBSLON
      ENDDO

      WRITE (IU08(1)) KOBSRLAT 
      WRITE (IU08(1)) KOBSRLON

      WRITE (IU08(2)) KOBSCOR

      RETURN
      END SUBROUTINE OUTUBUF
