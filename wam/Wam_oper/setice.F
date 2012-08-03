!-----------------------------------------------------------------------

      SUBROUTINE SETICE (FL3,MIJS,MIJL,IG)

!-----------------------------------------------------------------------

!**** *SETICE* ROUTINE TO SET SPECTRA ON ICE TO ZERO.

!     R.PORTZ      MPI         OKT.1992
!     J. BIDLOT    ECMWF       JUNE 1996  MESSAGE PASSING

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

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWICE   , ONLY : CICOVER  ,CITHRSH
      USE YOWPARAM , ONLY : NANG     ,NFRE

      REAL :: FL3(MIJS:MIJL,NANG,NFRE)
      REAL ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('SETICE',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

!*    1. SET SPECTRA TO ZERO AT ICE POINTS.
!     ------------------------------------

      DO M = 1, NFRE
        DO K = 1, NANG
          DO IJ = MIJS,MIJL
            IF(CICOVER(IJ,IG).GT.CITHRSH) THEN
              FL3(IJ,K,M) = 0.
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('SETICE',1,ZHOOK_HANDLE)

      RETURN

      END SUBROUTINE SETICE
