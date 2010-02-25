!-----------------------------------------------------------------------

      SUBROUTINE SETICE (FL3,MIJS,MIJL,IG)

!-----------------------------------------------------------------------

!**** *SETICE* ROUTINE TO SET SPECTRA ON ICE TO ZERO.

!     R.PORTZ      MPI         OKT.1992
!     J. BIDLOT    ECMWF       JUNE 1996  MESSAGE PASSING
!     P. JANSSEN   ECMWF       OCTOBER 1997 SIMPLIFY BY INTRODUCING
!                                           ICEMASK

!     PURPOSE.
!     -------

!          *SETICE* SET ICE SPECTRA (FL3) TO ZERO

!**   INTERFACE.
!     ----------

!         *CALL* *SETICE*

!     METHOD.
!     -------

!          NONE.


!     EXTERNALS.
!     ----------

!          .NONE.

!     REFERENCE.
!     ----------

!          NONE.

! ----------------------------------------------------------------------

      USE YOWICE   , ONLY : ICEMASK
      USE YOWPARAM , ONLY : NANG     ,NFRE

      REAL :: FL3(MIJS:MIJL,NANG,NFRE)

! ----------------------------------------------------------------------

!*    1. SET SPECTRA TO ZERO AT ICE POINTS.
!     ------------------------------------

      DO M = 1, NFRE
        DO K = 1, NANG
          DO IJ = MIJS,MIJL
            IF(ICEMASK(IJ,IG).EQ.0) THEN
              FL3(IJ,K,M) = 0.
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      RETURN

      END SUBROUTINE SETICE
