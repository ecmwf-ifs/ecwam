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
      USE YOWSTAT  , ONLY : MARSTYPE ,CDATEA   ,CDATEF   ,CDTPRO   ,
     &            IASSI    ,NENSFNB  ,NTOTENS  ,NSYSNB   ,NMETNB   ,
     &            LANAONLY ,L4VTYPE  ,ISTREAM
      USE YOWTEST  , ONLY : IU06 
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "wstream_strg.intfb.h"

      INTEGER :: KSTREAM
      REAL ::ZHOOK_HANDLE
      CHARACTER(LEN=4):: CSTREAM
      LOGICAL :: LASTREAM

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SETMARSTYPE',0,ZHOOK_HANDLE)
#endif

      IF (IASSI .EQ. 0 ) THEN
        MARSTYPE = 'an'
      ELSE
        IF(L4VTYPE) THEN
          MARSTYPE = '4v'
        ELSE
          MARSTYPE = 'fg'
        ENDIF
      ENDIF

      IF (CDTPRO.GE.CDATEF.OR.CDATEA.EQ.CDATEF) THEN
        CALL WSTREAM_STRG(ISTREAM,CSTREAM,NENSFNB,NTOTENS,MARSTYPE,
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

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SETMARSTYPE',1,ZHOOK_HANDLE)
#endif

      RETURN
      END SUBROUTINE SETMARSTYPE
