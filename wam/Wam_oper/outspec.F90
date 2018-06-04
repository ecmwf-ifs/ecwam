      SUBROUTINE OUTSPEC (SPEC) 

!----------------------------------------------------------------------

!**** *OUTSPEC*  ENCODES SPECTRA AS PARAMETER 251 USING GRIB API
!                AND WRITES TO FILE OR TO FDB.

!     J. BIDLOT   ECMWF  APRIL 2010 

!*    PURPOSE.
!     --------

!       ENCODES SPECTRA AND WRITES TO FDB OR TO A FILE.

!**   INTERFACE.
!     ----------

!     SUBROUTINE OUTSPEC (SPEC)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *SPEC*     REAL      LOCAL SPECTRA OF CURRENT PE.

!     METHOD.
!     -------

!           ENCODE SPECTRA PER FREQUENCY AND DIRECTION
!      INTO GRIB AND WRITE TO FDB OR TO A SINGLE FILE.

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!       NONE.

!-------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWCOU, LIFS_IO_SERV_ENABLED
      USE YOWCOUT  , ONLY : LWAM_USE_IO_SERV
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSTAT  , ONLY : CDATEE   ,CDATEF   ,CDTPRO   ,CDATEA   ,
     &            MARSTYPE
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOMHOOK  , ONLY : LHOOK, DR_HOOK

!-----------------------------------------------------------------------
      IMPLICIT NONE
#include "outwspec_io_serv.intfb.h"
#include "outwspec.intfb.h"
#include "difdate.intfb.h"

      REAL(KIND=JWRB), DIMENSION(NINF-1:NSUP, NANG, NFRE), INTENT(IN) :: SPEC

      INTEGER(KIND=JWIM) :: IFCST

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      CHARACTER(LEN=14) :: CDATE 

!-----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('OUTSPEC',0,ZHOOK_HANDLE)
      
      CALL GSTATS(2080,0)

      IF (ITEST.GT.1) THEN
        WRITE(IU06,*) '*      THIS IS OUTSPEC         *'
        CALL FLUSH (IU06)
      ENDIF

      IF(CDTPRO.LE.CDATEF) THEN
!*    0.1.  THIS IS AN ANALYSIS DATE.
        IF (LWCOU .AND. MARSTYPE .EQ. 'fg') THEN
          CDATE=CDATEA
          CALL DIFDATE (CDATEA, CDTPRO, IFCST)
          IFCST = IFCST/3600
        ELSEIF (LWCOU .AND. MARSTYPE .EQ. '4v') THEN
          CDATE=CDATEA
          CALL DIFDATE (CDATEA, CDTPRO, IFCST)
          IFCST = IFCST/3600
        ELSE
          CDATE=CDTPRO
          IFCST = 0
        ENDIF
      ELSE
!*    0.2.  THIS IS A  FORECAST DATE.
        CDATE=CDATEF
        CALL DIFDATE (CDATEF, CDTPRO, IFCST)
        IFCST = IFCST/3600
      ENDIF
!-----------------------------------------------------------------------

!*    3. OUTPUT GRIB DATA
!     -------------------
      
      ! Use IFS IO server?
      IF (LWCOU .AND. LIFS_IO_SERV_ENABLED .AND. LWAM_USE_IO_SERV) THEN
        CALL OUTWSPEC_IO_SERV(SPEC, MARSTYPE, CDATE, IFCST)
      ELSE
        CALL OUTWSPEC(SPEC, MARSTYPE, CDATE, IFCST)
      ENDIF

      CALL GSTATS(2080,1)

      IF (LHOOK) CALL DR_HOOK('OUTWSPEC',1,ZHOOK_HANDLE)

      END SUBROUTINE OUTSPEC
