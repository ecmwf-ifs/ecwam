SUBROUTINE OUTSPEC (IJS, IJL, FL, CICOVER)

!----------------------------------------------------------------------

!**** *OUTSPEC*  ENCODES SPECTRA AS PARAMETER 251 USING GRIB API
!                AND WRITES TO FILE OR TO FDB.

!     J. BIDLOT   ECMWF  APRIL 2010 

!*    PURPOSE.
!     --------

!       ENCODES SPECTRA AND WRITES TO FDB OR TO A FILE.

!**   INTERFACE.
!     ----------

!     SUBROUTINE OUTSPEC (IJS, IJL, FL, CICOVER)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *IJS:IJL - FIRST DIMENSION OF FL
!      *FL*        REAL      LOCAL SPECTRA OF CURRENT PE.
!      *CICOVER*   REAL      SEA ICE COVER.

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

      USE YOWCOUP  , ONLY : LWCOU, LIFS_IO_SERV_ENABLED,                &
                            OUTWSPEC_IO_SERV_HANDLER
      USE YOWCOUT  , ONLY : LWAM_USE_IO_SERV
      USE YOWICE   , ONLY : LICERUN  ,CITHRSH  ,FLMIN
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSTAT  , ONLY : CDATEE   ,CDATEF   ,CDTPRO   ,CDATEA   ,    &
     &            MARSTYPE ,LLSOURCE
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOMHOOK  , ONLY : LHOOK, DR_HOOK

!-----------------------------------------------------------------------
      IMPLICIT NONE
#include "outwspec_io_serv.intfb.h"
#include "outwspec.intfb.h"
#include "difdate.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL, NANG, NFRE), INTENT(IN) :: FL 
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP), INTENT(IN) :: CICOVER

      INTEGER(KIND=JWIM) :: IJ, K, M
      INTEGER(KIND=JWIM) :: IFCST

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL, NANG, NFRE) :: SPEC

      CHARACTER(LEN=14) :: CDATE 

!-----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('OUTSPEC',0,ZHOOK_HANDLE)
      
      CALL GSTATS(2080,0)

      IF (ITEST.GT.1) THEN
        WRITE(IU06,*) '*      THIS IS OUTSPEC         *'
        CALL FLUSH (IU06)
      ENDIF

!*    APPLY SEA ICE MASK TO THE OUTPUT SPECTRA (IF NEEDED)
      IF (LICERUN .AND. LLSOURCE) THEN
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ = IJS, IJL
              IF (CICOVER(IJ).GT.CITHRSH) THEN
                SPEC(IJ,K,M) = FLMIN 
              ELSE
                SPEC(IJ,K,M) = FL(IJ,K,M)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ELSE
        SPEC(:,:,:) = FL(:,:,:)
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
        IF (.NOT.ASSOCIATED(OUTWSPEC_IO_SERV_HANDLER)) CALL ABOR1('OUTWSPEC_IO_SERV_HANDLER IS NOT INITIALIZED')
        CALL OUTWSPEC_IO_SERV_HANDLER(IJS, IJL, SPEC, MARSTYPE, CDATE, IFCST)
      ELSE
        CALL OUTWSPEC(IJS, IJL, SPEC, MARSTYPE, CDATE, IFCST)
      ENDIF

      CALL GSTATS(2080,1)

      IF (LHOOK) CALL DR_HOOK('OUTWSPEC',1,ZHOOK_HANDLE)

END SUBROUTINE OUTSPEC
