      SUBROUTINE WRITEFL(FL,IJINF,IJSUP,MINF,MSUP,FILENAME,IUNIT,
     1                   LPBOPEN,LPBCLOSE,CDMARS,CDCLASS,CDEXPVER,
     2                   CDFDB2DSP)
C ----------------------------------------------------------------------
C     J. BIDLOT    ECMWF      MARCH 1997
C
C*    PURPOSE.
C     --------
C     WRITES ARRAY FL TO FILE.
C
C**   INTERFACE.
C     ----------
C     CALL *WRITEFL*(FL,IJINF,IJSUP,MINF,MSUP,FILENAME,IUNIT,
C    1               LPBOPEN,LPBCLOSE)
C     *FL*        ARRAY TO BE WRITTEN TO FILE
C     *IJINF*     FIRST LOWER DIMENSION BOUND OF FL
C     *IJSUP*     FIRST UPPER DIMENSION BOUND OF FL
C     note: NANG  is the SECOND DIMENSION
C     *MINF*      THIRD LOWER DIMENSION BOUND OF FL
C     *MSUP*      THIRD UPPER DIMENSION BOUND OF FL
C     *FILENAME*  FILENAME (INCLUDING PATH) OF TARGET FILE
C     *IUNIT*     PBIO UNIT (ONLY ACTIVE IF PBIO OUTPUT)
C     *LPBOPEN*   LOGICAL, TRUE IF PBOPEN HAS TO BE CALLED
C     *LPBCLOSE*  LOGICAL, TRUE IF PBCLOSE HAS TO BE CALLED
C     *CDMARS*    MARS TYPE OF FIELD.
C     *CDCLASS*   MARS CLASS OF FIELD.
C     *CDEXPVER*  EXPERIMENT VERSION IDENTIFIER.
C     *CDFDB2DSP* ROOT DIRECTORY FOR 2D SPECTRA NON GRIB DATA.
C
C     METHOD.
C     -------
C     WRITES ARRAY FL TO FILE. A FORTRAN WRITE, FDB WRITE OR A PBIO
C     WRITE IS USED DEPENDING ON THE CHOICE OF OUTPUT MODE AS SELECTED
C     BY LPBIOOUT
C
C     IF PBIO OR FDB IS USED, THEN FL WILL BE SAVED IN FILENAME USING
C     PBWRITE ON CONDITION THAT FILENAME WAS NOT OPENED BEFORE AND
C     CONNECTED TO IUNIT.
C
C     IT WILL BE APPENDED TO THE BOTTOM OF THE CONTENT OF FILENAME IF
C     FILENAME WAS ALREADY OPENED AND CONNECTED TO UNIT IUNIT.
C
C     IF FORTRAN WRITE IS USED THEN IN THE CASE OF MESSAGE PASSING
C     FL IS WRITTEN AS UNFORMATTED BINARY TO FILENAME CONNECTED
C     TO UNIT IU12. IN THE NON MESSAGE PASSING CASE, FL IS WRITTEN AS
C     UNFORMATTED BINARY TO FORTRAN UNIT SPECIFIED BY IU12.
C
C
C     EXTERNALS.
C     ----------
C     *PBOPEN*
C     *PBWRITE*
C     *PBCLOSE*
C     *ABORT1*
C     *IINITFDB*    - INITIALISE DATA BASE ELSEWHERE(!).
C     *IOPENFDB*    - OPEN       DATA BASE.
C     *ISETVALFDB*  - SET PARAMETERS ( CLASS, TYPE, DATE ETC )
C     *IWRITEFDB*   - WRITE FIELD INTO DATA BASE..
C     *ICLOSEFDB*   - CLOSE      DATA BASE.
C
C     REFERENCE.
C     ----------
C     NONE
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "comcout.h"
C
#include "comgrid.h"
C
#include "comtest.h"
C
#include "comtext.h"
C
#include "comunit.h"
C
#include "txtunit.h"
C
#include "comwind.h"
C
#include "commesspass.h"
C
#include "commpp.h"
C
C     HELPS TO CALL IINITFDB ONLY ONCE.
C     ---------------------------------
C
      LOGICAL       LINITFDB
      CHARACTER*  1 CLMODE
      DATA          LINITFDB /.TRUE./,
     1              CLMODE   / "w" /
C
C ----------------------------------------------------------------------
C
      CHARACTER*10 CDTPRO,CDATEF,CDATER,CDATES
      CHARACTER*80 FILENAME
      LOGICAL LPBOPEN,LPBCLOSE
      CHARACTER *1 CLDOMAIN
      CHARACTER *2 CDMARS
      CHARACTER *2 CDCLASS
      CHARACTER *4 CDEXPVER
      CHARACTER *6 CLDATE
      CHARACTER *2 CLTIME
      CHARACTER *4 CLSTEP
      CHARACTER *4 CLFREQ
      CHARACTER *4 CLSTREAM
      CHARACTER *(*) CDFDB2DSP
C
      DATA CLDOMAIN  / region /
C
C ----------------------------------------------------------------------
C
      REAL,DIMENSION(IJINF:IJSUP,NANG,MINF:MSUP) :: FL
C
C ----------------------------------------------------------------------
      IF(LPBIOOUT .AND. .NOT. LFDBIOOUT) THEN
        IF(LPBOPEN) THEN
          CALL PBOPEN(IUNIT,FILENAME,CLMODE,KRET)
          IF(KRET.LT.0) THEN
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBOPEN       *'
          WRITE (IU06,*) '*   IN WRITEFL                           *'
          IF(KRET.EQ.-1)
     1        WRITE (IU06,*) 'COULD NOT OPEN FILE ',FILENAME
          IF(KRET.EQ.-2) WRITE (IU06,*) 'INVALID FILENAME ',FILENAME
          IF(KRET.EQ.-3) WRITE (IU06,*) 'INVALID OPEN MODE SPECIFIED'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          CALL ABORT1
          ENDIF
        ENDIF
        KOUNT = (IJSUP-IJINF+1)*NANG*(MSUP-MINF+1)*NPRECI
        CALL PBWRITE(IUNIT,FL,KOUNT,KRET)
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '*****************************************'
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*    ERROR FOLLOWING CALL TO PBWRITE    *'
          WRITE (IU06,*) '*    IN WRITEFL                         *'
          WRITE (IU06,*) '*    FILE ',FILENAME
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*****************************************'
          CALL ABORT1
        ENDIF
        IF(LPBCLOSE) THEN
          CALL PBCLOSE(IUNIT,KRET)
          IF(KRET.LT.0) THEN
             WRITE (IU06,*) '************************************'
             WRITE (IU06,*) '*                                  *'
             WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
             WRITE (IU06,*) '* IN WRITEFL                       *'
             WRITE (IU06,*) '* FILE ',FILENAME
             WRITE (IU06,*) '*                                  *'
             WRITE (IU06,*) '************************************'
             CALL ABORT1
          ENDIF
        ENDIF
      ELSEIF(LFDBIOOUT) THEN
        IF (LINITFDB) THEN
          ISTAT=IINITFDB()
          LINITFDB = .FALSE.
        ENDIF
C
        IF(LPBOPEN) THEN
          ICDFDB2DSP = LEN_TRIM(CDFDB2DSP)
          ISTAT = IOPENFDB ('fdb', IUNIT, CLMODE   )
          IF (ISTAT .NE. 0) THEN
            WRITE(IU06,'("Error\ /writefl/ Action specified: ",a)')
     .      "open"
            WRITE(IU06,'("Error\ /writefl/ iopenfdb return status:",
     .      I3)') ISTAT
            CALL ABORT1
          ELSE
            WRITE(IU06,'(" OPEN FDB ON ROOT ",a)')
     1            CDFDB2DSP(1:ICDFDB2DSP)
            CALL FLUSH(IU06)
          ENDIF
C
          ISTAT = ISET_FDB_ROOT(IUNIT,    CDFDB2DSP(1:ICDFDB2DSP) )
C
          IF (ISTAT .NE. 0) THEN
            WRITE (IU06,'("ERROR FOLLOWING CALL TO SET_FDB_ROOT status",
     1                    I3)') ISTAT
            CALL ABORT1
          ENDIF
          ISTAT = ISETVALFDB (IUNIT, 'levt',   's'        )
          ISTAT = ISETVALFDB (IUNIT, 'level',  '0000'     )
          ISTAT = ISETVALFDB (IUNIT, 'param',  '250'      )
          ISTAT = ISETVALFDB (IUNIT, 'repres', 'll'       )
        ENDIF
C
        IL1 = INDEX(FILENAME,'BLS',.TRUE.) + 3
        IL2 = IL1 + 5
        CLDATE = FILENAME(IL1:IL2)
        CLTIME = FILENAME(IL2+1: IL2+2)
        CLSTEP = FILENAME(IL2+3: IL2+6)
        WRITE(CLFREQ,'(I4.4)') MINF

        IF (ITEST .GE. 2)  THEN
          WRITE(IU06,'(" WRITEFL TO FDB: date=", a6, ". time=", a2, 
     .                 ". step=", a4, ". type=", a2, ". freq=", a4,
     .                 ". il1=", i2, ". il2=", i2, ".")')
     .    CLDATE, CLTIME, CLSTEP, CDMARS, CLFREQ, IL1, IL2
          CALL FLUSH (IU06)
        ENDIF
C
        IF (CDMARS .EQ. 'cf' .OR. CDMARS .EQ. 'pf') THEN
          CLSTREAM = 'waef'
        ELSE
          CLSTREAM = 'wave'
        ENDIF

        ISTAT = ISETVALFDB (IUNIT, 'STREAM', CLSTREAM )
        ISTAT = ISETVALFDB (IUNIT, 'class',  CDCLASS  )
        ISTAT = ISETVALFDB (IUNIT, 'stream', 'wave'   )
        ISTAT = ISETVALFDB (IUNIT, 'domain', 'g'      )
        ISTAT = ISETVALFDB (IUNIT, 'expver', CDEXPVER )
        ISTAT = ISETVALFDB (IUNIT, 'date',   CLDATE   )
        ISTAT = ISETVALFDB (IUNIT, 'type',   CDMARS   )
        ISTAT = ISETVALFDB (IUNIT, 'time',   CLTIME   )
        ISTAT = ISETVALFDB (IUNIT, 'fstep',  'spectr' )
        ISTAT = ISETVALFDB (IUNIT, 'step',   CLSTEP   )
        ISTAT = ISETVALFDB (IUNIT, 'leve',   CLFREQ   )
C
        KOUNT = (IJSUP-IJINF+1)*NANG*(MSUP-MINF+1) * NPRECI / 4
C
C       KOUNT should be number of values;
C
C       ISTAT = IWRITEFDB (IUNIT, FL(IOFFSET),KOUNT)
        ISTAT = IWRITEFDB (IUNIT, FL, KOUNT)
C
        IF (ISTAT .NE. 0) THEN
          WRITE (IU06,'("ERROR FOLLOWING CALL TO IWRITEFDB status ",
     1                   I4)') ISTAT
          CALL ABORT1
        ENDIF
C
        IF(LPBCLOSE) THEN
          ISTAT=ICLOSEFDB(IUNIT)
          IF (ISTAT .NE. 0) THEN
            WRITE (IU06,'("ERROR FOLLOWING CALL TO ICLOSEFDB status "
     1                    , I4)') ISTAT
            CALL ABORT1
          ENDIF
        ENDIF
C
      ELSE
        IF(LMESSPASS) THEN
          OPEN(IU12,FILE=FILENAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
C
          WRITE(IU12) (((FL(J1,J2,J3),
     .                      J1=IJINF,IJSUP),
     .                      J2=1,NANG),
     .                      J3=MINF,MSUP)
          IF (ITEST.GE.2)
     1      WRITE(IU06,*) ' SUB. WRITEFL: OUTPUT OF FL DONE'
          CLOSE(IU12)
        ELSE
          IF (IGL.EQ.1) THEN
            REWIND IU12
            WRITE(IU12) FL
            IF (ITEST.GE.2)
     1        WRITE(IU06,*) ' SUB. WRITEFL: OUTPUT OF FL DONE'
          ENDIF
        ENDIF
      ENDIF
C
      RETURN
      END
