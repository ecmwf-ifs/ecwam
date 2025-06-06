! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE WVWAMINIT (LLCOUPLED, IULOG, LLRNL,                    &
     &                      NGAUSSW, NLON, NLAT, RSOUTW, RNORTW)

!****  *WVWAMINIT* - INITIALISES A FEW VARIABLE AND READ WAVE MODEL
!                    NAMELIST (if LLRNL is true) and
!                    CONFIGURATION PARAMETERS AND TABLES.


!     J. BIDLOT    ECMWF   FEBRUARY 2008: MOVE THIS PART OUT OF MPDECOMP 

!*    INTERFACE.
!     ----------
!     CALL *WVWAMINIT* (LLCOUPLED,IULOG,LLRNL,
!                       NGAUSSW, NLON, NLAT, RSOUTW, RNORTW)

!      *LLCOUPLED*  LOGICAL  TRUE IF ATMOSPHERIC COUPLED RUN.
!      *IULOG*      INTEGER  STANDARD OUTPUT UNIT USED BY IFS.
!      *LLRNL*      LOGICAL  IF TRUE READS WAM NAMELIST.
!      *NGAUSSW*    INTEGER  1: GAUSSIAN GRID, 0: LAT-LON GRID
!      *NLON*       INTEGER  MAXIMUM NUMBER OF LONGITUDES USED BY WAM. 
!      *NLAT*       INTEGER  MAXIMUM NUMBER OF LATITUDES USED BY WAM. 
!      *RSOUTW*     REAL     MOST SOUTHERN LATITUDE IN WAM.
!      *RNORTW*     REAL     MOST NORTHERn LATITUDE IN WAM.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMAP   , ONLY : AMOSOP   ,AMONOP   ,IQGAUSS  ,NGX      ,NGY
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : KWAMVER  ,LLUNSTR
      
      USE YOWTEST  , ONLY : IU06
      USE YOWUNIT  , ONLY : LWVWAMINIT
      USE YOWSTAT  , ONLY : IPROPAGS

      USE MPL_MODULE, ONLY : MPL_MYRANK, MPL_NPROC
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "expand_string.intfb.h"
#include "iniwcst.intfb.h"
#include "iwam_get_unit.intfb.h"
#include "mfredir.intfb.h"
#include "mpuserin.intfb.h"
#include "setwavphys.intfb.h"
#include "readmdlconf.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IULOG
      INTEGER(KIND=JWIM), INTENT(OUT) :: NGAUSSW, NLON, NLAT
      REAL(KIND=JWRB), INTENT(OUT) :: RSOUTW, RNORTW
      LOGICAL, INTENT(IN) :: LLRNL
      LOGICAL, INTENT(IN) :: LLCOUPLED


      REAL(KIND=JWRB) :: PRPLRADI
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      CHARACTER(LEN=80) :: LOGFILENAME

      LOGICAL :: LLEXIST
      LOGICAL, SAVE :: LFRST

      DATA LFRST /.TRUE./

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WVWAMINIT',0,ZHOOK_HANDLE)

      LWVWAMINIT=.FALSE.

      IRANK = MPL_MYRANK()
      NPROC = MPL_NPROC()

!     STANDARD OUTPUT UNIT
!     --------------------

      IF (LLCOUPLED) THEN
        IU06=IULOG
      ELSE
        IU06=65
        LLEXIST=.TRUE.
        DO WHILE (LLEXIST)
          IU06=IU06+1
          CALL GET_ENVIRONMENT_VARIABLE("ECWAM_LOG_FILE",LOGFILENAME)
          IF (LOGFILENAME == '') LOGFILENAME='wam.log'
          IF (NPROC > 1 .AND. INDEX(LOGFILENAME, '%p') == 0) LOGFILENAME=TRIM(LOGFILENAME)//".%p"
          CALL EXPAND_STRING(IRANK,NPROC,0,0,LOGFILENAME,1)
          INQUIRE(UNIT=IU06, OPENED=LLEXIST)
          IF (.NOT. LLEXIST) THEN
            OPEN(IU06,FILE=LOGFILENAME,STATUS='UNKNOWN')
            WRITE(IU06,'(A,I0,A,I0)') '  # WAM LOG FILE OF PE ', IRANK, ' / ', NPROC
            WRITE(IU06,*) ' '
          ENDIF
        ENDDO
      ENDIF

      WRITE(IU06,*) ' WAM SOFTWARE VERSION: ', KWAMVER
      WRITE(IU06,*) ' '

!     GET CONTROLLING FLAG FROM INPUT NAMELIST 
!     ----------------------------------------
      IF (LFRST .AND. LLRNL) THEN
        CALL MPUSERIN
        LFRST=.FALSE.
      ENDIF

!     INITIALISE THE FREQUENCY AND DIRECTION ARRAYS
!     ---------------------------------------------
      PRPLRADI = 1.0_JWRB
      CALL INIWCST(PRPLRADI)  !! it will be call again in INITMDL once PRPLRADI is known
      CALL MFREDIR

!     SET WAVE PHYSICS PACKAGE
!     ------------------------

      CALL SETWAVPHYS

!     READ MODEL GRID CONFIGURATION
!     -----------------------------

      CALL READMDLCONF

      WRITE(IU06,*) ' WVWAMINT: WAVE MODEL GRID CONFIGURATION READ IN'
      CALL FLUSH (IU06)


!     RETURN NECESSARY PARAMETERS
      NGAUSSW = IQGAUSS
      NLON=NGX
      NLAT=NGY
      RSOUTW=AMOSOP
      RNORTW=AMONOP

      IF (LHOOK) CALL DR_HOOK('WVWAMINIT',1,ZHOOK_HANDLE)

      END SUBROUTINE WVWAMINIT