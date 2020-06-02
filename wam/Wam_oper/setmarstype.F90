      SUBROUTINE SETMARSTYPE

! ----------------------------------------------------------------------

!**** *SETMARSTYPE* -

!     J. BIDLOT     ECMWF   JANUARY 1998  SETS MARSTYPE

!*    PURPOSE.
!     --------
!      SETS VARIABLE MARSTYPE

!**   INTERFACE.
!     ----------
!     *CALL SETMARSTYPE*

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWSTAT  , ONLY : MARSTYPE ,CDATEA   ,CDATEF   ,CDTPRO   ,    &
     &            IASSI    ,NENSFNB  ,NTOTENS  ,NSYSNB   ,NMETNB   ,    &
     &            LANAONLY ,L4VTYPE  ,ISTREAM  ,IASSI_ORIG
      USE YOWTEST  , ONLY : IU06
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
      USE ALGORITHM_STATE_MOD, ONLY : GET_NUPTRA, GET_MUPTRA, &
     &                                GET_ALGOR_TYPE
      USE YOWCOUT  , ONLY : LWAMANOUT, LWAMANOUT_ORIG

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "wstream_strg.intfb.h"

      INTEGER(KIND=JWIM) :: KSTREAM
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      CHARACTER(LEN=4):: CSTREAM
      LOGICAL :: LASTREAM
      LOGICAL, SAVE :: LFRST_OOPS

      DATA LFRST_OOPS /.TRUE./

! ---------------------------------------------------------------------

!*    1.  THE FIRST CALL TO WAVEMDL PERFORMS INITIALIZATION.

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SETMARSTYPE',0,ZHOOK_HANDLE)

      IF (GET_ALGOR_TYPE() == 'OOPS') THEN
        ! OOPS-IFS may do wave assimilation only in the final outer loop
        IF (LFRST_OOPS) THEN
          LFRST_OOPS = .FALSE.
          IASSI_ORIG = IASSI
          LWAMANOUT_ORIG = LWAMANOUT
        ENDIF
        IF (GET_NUPTRA() /= GET_MUPTRA() - 1) THEN
          IASSI = 0
          LWAMANOUT = .FALSE.
        ELSE
          IASSI = IASSI_ORIG
          LWAMANOUT = LWAMANOUT_ORIG
        ENDIF
        WRITE(IU06,*) ' SUB. WAVEMDL CALLED FROM ', GET_ALGOR_TYPE(), &
          &           ' FOR NUPTRA: ', GET_NUPTRA(), ' AND MUPTRA: ', &
          &           GET_MUPTRA(), ' --> IASSI reset to', IASSI
      ENDIF

      IF (IASSI .EQ. 0 ) THEN
        MARSTYPE = 'an'
      ELSE
        IF(L4VTYPE) THEN
          MARSTYPE = '4v'
        ELSE
          MARSTYPE = 'fg'
        ENDIF
      ENDIF

      IF (CDTPRO.GE.CDATEF .OR. CDATEA.EQ.CDATEF) THEN
        CALL WSTREAM_STRG(ISTREAM,CSTREAM,NENSFNB,NTOTENS,MARSTYPE,     &
     &                    KSTREAM, LASTREAM)
        IF(CSTREAM.EQ.'****') THEN
          WRITE(IU06,*) '*****************************************'
          WRITE(IU06,*) ''
          WRITE(IU06,*) ' ERROR IN SETMARSTYPE !!!!'
          WRITE(IU06,*) ' FORECAST STREAM UNKNOWN '
          WRITE(IU06,*) ' INPUT ISTREAM = ', ISTREAM
          WRITE(IU06,*) ' BUT NOT DEFINED IN WSTREAM_STRG !!!!'
          WRITE(IU06,*) ''
          WRITE(IU06,*) '*****************************************'
          CALL ABORT1
        ENDIF
      ENDIF

      IF(LANAONLY) THEN
        IF (IASSI .EQ. 0 ) THEN
          MARSTYPE = 'an'
        ELSE
          MARSTYPE = 'fg'
        ENDIF
      ENDIF

      IF (LHOOK) CALL DR_HOOK('SETMARSTYPE',1,ZHOOK_HANDLE)

      END SUBROUTINE SETMARSTYPE
