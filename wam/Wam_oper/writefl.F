      SUBROUTINE WRITEFL(FL,IJINF,IJSUP,MINF,MSUP,FILENAME,IUNIT,
     &                   LPBOPEN,LPBCLOSE,CDMARS,CDCLASS,CDEXPVER,
     &                   CDFDB2DSP)
! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      MARCH 1997

!*    PURPOSE.
!     --------
!     WRITES ARRAY FL TO FILE.

!**   INTERFACE.
!     ----------
!     CALL *WRITEFL*(FL,IJINF,IJSUP,MINF,MSUP,FILENAME,IUNIT,
!    1               LPBOPEN,LPBCLOSE)
!     *FL*        ARRAY TO BE WRITTEN TO FILE
!     *IJINF*     FIRST LOWER DIMENSION BOUND OF FL
!     *IJSUP*     FIRST UPPER DIMENSION BOUND OF FL
!     note: NANG  is the SECOND DIMENSION
!     *MINF*      THIRD LOWER DIMENSION BOUND OF FL
!     *MSUP*      THIRD UPPER DIMENSION BOUND OF FL
!     *FILENAME*  FILENAME (INCLUDING PATH) OF TARGET FILE
!     *IUNIT*     PBIO UNIT (ONLY ACTIVE IF PBIO OUTPUT)
!     *LPBOPEN*   LOGICAL, TRUE IF PBOPEN HAS TO BE CALLED
!     *LPBCLOSE*  LOGICAL, TRUE IF PBCLOSE HAS TO BE CALLED
!     *CDMARS*    MARS TYPE OF FIELD.
!     *CDCLASS*   MARS CLASS OF FIELD.
!     *CDEXPVER*  EXPERIMENT VERSION IDENTIFIER.
!     *CDFDB2DSP* ROOT DIRECTORY FOR 2D SPECTRA NON GRIB DATA.

!     METHOD.
!     -------
!     WRITES ARRAY FL TO FILE. A FORTRAN WRITE, FDB WRITE OR A PBIO
!     WRITE IS USED DEPENDING ON THE CHOICE OF OUTPUT MODE AS SELECTED
!     BY LPBIOOUT

!     IF PBIO OR FDB IS USED, THEN FL WILL BE SAVED IN FILENAME USING
!     PBWRITE ON CONDITION THAT FILENAME WAS NOT OPENED BEFORE AND
!     CONNECTED TO IUNIT.

!     IT WILL BE APPENDED TO THE BOTTOM OF THE CONTENT OF FILENAME IF
!     FILENAME WAS ALREADY OPENED AND CONNECTED TO UNIT IUNIT.

!     IF FORTRAN WRITE IS USED THEN IN THE CASE OF MESSAGE PASSING
!     FL IS WRITTEN AS UNFORMATTED BINARY TO FILENAME CONNECTED
!     TO UNIT IU12. IN THE NON MESSAGE PASSING CASE, FL IS WRITTEN AS
!     UNFORMATTED BINARY TO FORTRAN UNIT SPECIFIED BY IU12.


!     EXTERNALS.
!     ----------
!     *PBOPEN*
!     *PBWRITE*
!     *PBCLOSE*
!     *ABORT1*
!     *IOPENFDB*    - OPEN       DATA BASE.
!     *ISETVALFDB*  - SET PARAMETERS ( CLASS, TYPE, DATE ETC )
!     *IWRITEFDB*   - WRITE FIELD INTO DATA BASE..
!     *ICLOSEFDB*   - CLOSE      DATA BASE.

!     REFERENCE.
!     ----------
!     NONE

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NANG     ,CLDOMAIN
      USE YOWGRID  , ONLY : IGL
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWUNIT  , ONLY : IU12
      USE YOWMESPAS, ONLY : LMESSPASS,LPBIOOUT ,LFDBIOOUT
      USE YOWMPP   , ONLY : NPRECR
      USE YOMGSTATS, ONLY : LSYNCSTATS

! ----------------------------------------------------------------------

      CHARACTER*  1 CLMODE
      DATA          CLMODE   / "w" /
      CHARACTER*120 FILENAME
      LOGICAL LPBOPEN,LPBCLOSE
      CHARACTER *2 CDMARS
      CHARACTER *2 CDCLASS
      CHARACTER *4 CDEXPVER
      CHARACTER *6 CLDATE
      CHARACTER *2 CLTIME
      CHARACTER *6 CLSTEP
      CHARACTER *4 CLFREQ
      CHARACTER *4 CLSTREAM
      CHARACTER *(*) CDFDB2DSP

! ----------------------------------------------------------------------

      REAL,DIMENSION(IJINF:IJSUP,NANG,MINF:MSUP) :: FL

! ----------------------------------------------------------------------
      IF(LPBIOOUT .AND. .NOT. LFDBIOOUT) THEN
        IF(LPBOPEN) THEN
          CALL PBOPEN(IUNIT,FILENAME,CLMODE,KRET)
          IF(KRET.LT.0) THEN
            WRITE (IU06,*) '******************************************'
            WRITE (IU06,*) '*                                        *'
            WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBOPEN       *'
            WRITE (IU06,*) '*   IN WRITEFL                           *'
            IF(KRET.EQ.-1)
     &       WRITE (IU06,*) 'COULD NOT OPEN FILE ',FILENAME
            IF(KRET.EQ.-2) WRITE (IU06,*) 'INVALID FILENAME ',FILENAME
            IF(KRET.EQ.-3) WRITE (IU06,*) 'INVALID OPEN MODE SPECIFIED'
            WRITE (IU06,*) '*                                        *'
            WRITE (IU06,*) '******************************************'
            CALL ABORT1
          ENDIF
        ENDIF
        KOUNT = (IJSUP-IJINF+1)*NANG*(MSUP-MINF+1)*NPRECR
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

        WRITE (IU06,*) '************************************'
        WRITE (IU06,*) '*                                  *'
        WRITE (IU06,*) '* IN WRITEFL                       *'
        WRITE (IU06,*) '* THE OPTION TO WRITE BINARY FILE  *'
        WRITE (IU06,*) '* TO FDB IS OBSOLETE               *'
        WRITE (IU06,*) '* PROGRAM WILL ABORT               *'
        WRITE (IU06,*) '*                                  *'
        WRITE (IU06,*) '************************************'
        CALL ABORT1

        IF(LPBOPEN) THEN
          ICDFDB2DSP = LEN_TRIM(CDFDB2DSP)
          ISTAT = IOPENFDB ('fdb', IUNIT, CLMODE   )
          IF (ISTAT .NE. 0) THEN
            WRITE(IU06,'("Error\ /writefl/ Action specified: ",a)')
     &      "open"
            WRITE(IU06,'("Error\ /writefl/ iopenfdb return status:",
     &      I3)') ISTAT
            CALL ABORT1
          ELSEIF (ITEST .GE. 2 ) THEN
            WRITE(IU06,'(" OPEN FDB ON ROOT ",a)')
     &            CDFDB2DSP(1:ICDFDB2DSP)
            CALL FLUSH(IU06)
          ENDIF

          ISTAT = ISET_FDB_ROOT(IUNIT,    CDFDB2DSP(1:ICDFDB2DSP) )

          IF (ISTAT .NE. 0) THEN
            WRITE (IU06,'("ERROR FOLLOWING CALL TO SET_FDB_ROOT status",
     &                    I3)') ISTAT
            CALL ABORT1
          ENDIF
          ISTAT = ISETVALFDB (IUNIT, 'levt',   's'        )
          ISTAT = ISETVALFDB (IUNIT, 'level',  '0000'     )
          ISTAT = ISETVALFDB (IUNIT, 'param',  '250'      )
          ISTAT = ISETVALFDB (IUNIT, 'repres', 'll'       )
        ENDIF

        IL1 = INDEX(FILENAME,'BLS',.TRUE.) + 3
        IL2 = IL1 + 5
        CLDATE = FILENAME(IL1:IL2)
        CLTIME = FILENAME(IL2+1: IL2+2)
        CLSTEP = '00'//FILENAME(IL2+3: IL2+6)
        WRITE(CLFREQ,'(I4.4)') MINF

        IF (ITEST .GE. 2)  THEN
          WRITE(IU06,'(" WRITEFL TO FDB: date=", a6, ". time=", a2, 
     &                 ". step=", a6, ". type=", a2, ". freq=", a4,
     &                 ". il1=", i2, ". il2=", i2, ".")')
     &    CLDATE, CLTIME, CLSTEP, CDMARS, CLFREQ, IL1, IL2
          CALL FLUSH (IU06)
        ENDIF

        IF (CDMARS .EQ. 'cf' .OR. CDMARS .EQ. 'pf') THEN
          CLSTREAM = 'waef'
        ELSE
          CLSTREAM = 'wave'
        ENDIF

        ISTAT = ISETVALFDB (IUNIT, 'STREAM', CLSTREAM )
        ISTAT = ISETVALFDB (IUNIT, 'class',  CDCLASS  )
        ISTAT = ISETVALFDB (IUNIT, 'stream', 'wave'   )
        ISTAT = ISETVALFDB (IUNIT, 'domain', CLDOMAIN )
        ISTAT = ISETVALFDB (IUNIT, 'expver', CDEXPVER )
        ISTAT = ISETVALFDB (IUNIT, 'date',   CLDATE   )
        ISTAT = ISETVALFDB (IUNIT, 'type',   CDMARS   )
        ISTAT = ISETVALFDB (IUNIT, 'time',   CLTIME   )
        ISTAT = ISETVALFDB (IUNIT, 'fstep',  'spectr' )
        ISTAT = ISETVALFDB (IUNIT, 'step',   CLSTEP   )
        ISTAT = ISETVALFDB (IUNIT, 'leve',   CLFREQ   )

        KOUNT = (IJSUP-IJINF+1)*NANG*(MSUP-MINF+1) * NPRECR / 4

!       KOUNT should be number of values;

!       ISTAT = IWRITEFDB (IUNIT, FL(IOFFSET),KOUNT)
        ISTAT = IWRITEFDB (IUNIT, FL, KOUNT)

        IF (ISTAT .NE. 0) THEN
          WRITE (IU06,'("ERROR FOLLOWING CALL TO IWRITEFDB status ",
     &                   I4)') ISTAT
          CALL ABORT1
        ENDIF

        IF(LPBCLOSE) THEN
          ISTAT=ICLOSEFDB(IUNIT)
          IF (ISTAT .NE. 0) THEN
            WRITE (IU06,'("ERROR FOLLOWING CALL TO ICLOSEFDB status "
     &                    , I4)') ISTAT
            CALL ABORT1
          ENDIF
        ENDIF

      ELSE
        IF(LMESSPASS) THEN
          OPEN(IU12,FILE=FILENAME,STATUS='UNKNOWN',FORM='UNFORMATTED')

          WRITE(IU12) (((FL(J1,J2,J3),
     &                      J1=IJINF,IJSUP),
     &                      J2=1,NANG),
     &                      J3=MINF,MSUP)
          IF (ITEST.GE.2)
     &      WRITE(IU06,*) ' SUB. WRITEFL: OUTPUT OF FL DONE'
          CLOSE(IU12)
        ELSE
          IF (IGL.EQ.1) THEN
            REWIND IU12
            WRITE(IU12) FL
            IF (ITEST.GE.2)
     &        WRITE(IU06,*) ' SUB. WRITEFL: OUTPUT OF FL DONE'
          ENDIF
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE WRITEFL
