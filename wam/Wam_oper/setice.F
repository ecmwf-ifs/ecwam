!ICE
!-----------------------------------------------------------------------

      SUBROUTINE SETICE (FL3,IJS,IJL,IG)

!-----------------------------------------------------------------------

!**** *SETICE* ROUTINE TO SET SPECTRA ON ICE TO ZERO.

!     R.PORTZ      MPI         OKT.1992
!     J. BIDLOT    ECMWF       JUNE 1996  MESSAGE PASSING
!     P. JANSSEN   ECMWF       OCTOBER 1997 SIMPLIFY BY INTRODUCING
!                                           ICEMASK

!     PURPOSE.
!     -------

!          *READICE* SET ICE SPECTRA (FL3) TO ZERO

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
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE

      DIMENSION FL3(NINF-1:NSUP,NANG,NFRE)

! ----------------------------------------------------------------------

!*    1. SET SPECTRA TO ZERO AT ICE POINTS.
!     ------------------------------------

      DO M = 1, NFRE
        DO K = 1, NANG
          DO IJ = IJS,IJL
            FL3(IJ,K,M) = ICEMASK(IJ,IG)*FL3(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO

      RETURN

      END SUBROUTINE SETICE
!ICE
