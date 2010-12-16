      SUBROUTINE WGRIBOUT (IU06, ITEST, I1, I2, FIELD,
     &                     ITABLE, IPARAM, KLEV, IK, IM,
     &                     CDATE, IFCST, MARSTYPE,
     &                     LFDB, CDFDBSF, KFDB, LFDBOPEN, IU)

! ----------------------------------------------------------------------

!****  *WGRIBOUT*  ENCODES WAM MODEL FIELD INTO GRIB CODE AND OUTPUT

!       J. BIDLOT    ECMWF JULY 2009: USE GRIB API 

!       PURPOSE.
!       --------
!         SUBROUTINE PACKS WAVE FIELDS INTO THE GRIB CODE
!         AND ARCHIVE INTO FDB OR INTO FILE.
!         !!!! PRESET_WGRIB_TEMPLATE NEEDS to BE CALLED BEFORE TO
!         !!!! INITIALISE THE GRIB TEMPLATES.

!**    INTERFACE.
!      ----------
!        *CALL* *WGRIBOUT (IU06, ITEST, I1, I2, FIELD,
!                          ITABLE, IPARAM, KLEV, IK, IM,
!                          CDATE, IFCST,
!                          LFDB, CDFDBSF, KFDB, LFDBOPEN, IU,
!          *IU06*    LOGFILE OUTPUT UNIT.
!          *ITEST*   TEST OUTPUT GIVEN IF ITEST GT 2.
!          *I1*      FIRST DIMENSION OF FIELD.
!          *I2*      SECOND DIMENSION OF FIELD.
!          *FIELD*   FIELD TO BE PACKED.
!          *ITABLE*  GRIB TABLE NUMBER.
!          *IPARAM*  PARAMETER IDENTIFIER.
!          *KLEV*    REFERENCE LEVEL IN full METER
!                    (SHOULD BE 0 EXCEPT FOR 233 AND 245)
!          *IK*      DIRECTION INDEX,
!                    ONLY MEANINGFUL FOR SPECTRAL PARAMETERS.
!          *IM*      FREQUENCY INDEX,
!                    ONLY MEANINGFUL FOR SPECTRAL PARAMETERS.
!          *CDATE*   DATE AND TIME.
!          *IFCST*   FORECAST STEP IN HOURS.
!          *MARSTYPE* TYPE OF CURRENT FIELD
!           *LFDB*    LOGICAL SWITCH TO ACTIVATE FDB PROCESSING
!                    .T. TO ACTIVATE, .F. TO IGNORE
!          *IU*      LOGICAL UNIT FOR PACKED FIELD DATA.
!                    !!! IT WILL NEED TO BE OPEN AN
!                    !!! WITH GRIB_API SOFTWARE !!!

!      METHOD.
!      -------

!      EXTERNALS.
!      ----------

!      REFERENCES.
!      -----------

! ----------------------------------------------------------------------

      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRIB_HANDLES , ONLY :NGRIB_HANDLE_WAM_I,NGRIB_HANDLE_WAM_S
      USE YOWGRIBHD, ONLY : PPMISS   , PPEPS    ,PPREC    ,NTENCODE ,
     &            NGRBRESS ,PPRESOL  ,LGRHDIFS ,LNEWLVTP  ,PPMIN_RESET,
     &            LPADPOLES
      USE YOWGRID  , ONLY : NLONRGG
      USE YOWMAP   , ONLY : IRGG     ,AMONOP   ,AMOSOP   ,XDELLA
      USE YOWMPP   , ONLY : NPRECI   ,IRANK
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,CLDOMAIN
      USE YOWCOUP  , ONLY : KCOUSTEP
      USE YOWCOUT  , ONLY : LRSTST0
      USE YOWPCONS , ONLY : G        ,ZMISS    ,DEG
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH
      USE GRIB_API_INTERFACE
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: IU06, I1, I2, ITABLE, IPARAM, KLEV, IK, IM, IU, ITEST
      INTEGER :: KFDB
      INTEGER :: ICLASS, IFCST, ISTEP, ISTEP_HRS 
      INTEGER :: IC, JC, ITABPAR, IDATE, ITIME
      INTEGER :: ICOUNT, NN, I, J, JSN, KK, MM
      INTEGER :: IERR
      INTEGER :: IGRIB_HANDLE
      INTEGER :: KBYTES, ISIZE
      INTEGER, ALLOCATABLE :: KGRIB_BUFR(:)

      REAL :: TEMP
      REAL :: PMISS 
      REAL :: PPMAX, PPMIN, DELTAPP, ABSPPREC
      REAL :: ZHOOK_HANDLE
      REAL :: FIELD(I1,I2)
      REAL, ALLOCATABLE :: VALUES(:)

      CHARACTER(LEN=*) :: CDFDBSF
      CHARACTER(LEN=2) :: MARSTYPE
      CHARACTER(LEN=4) :: CSTREAM
      CHARACTER(LEN=4) :: CEXPVER
      CHARACTER(LEN=12) :: C12
      CHARACTER(LEN=14) :: CDATE
      CHARACTER(LEN=24) :: FILNM

      LOGICAL :: LFDB
      LOGICAL :: LFDBOPEN

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WGRIBOUT',0,ZHOOK_HANDLE)

      CALL GSTATS(1709,0)

      IF(ITEST.GT.0) THEN
        WRITE(IU06,*) '   SUB. WGRIBOUT CALLED FOR ',IPARAM
        CALL FLUSH(IU06)
      ENDIF

!     CLONE GRIB TEMPLATE:

      IF(IPARAM.EQ.251) THEN
        IF(LGRHDIFS) THEN
          CALL PRESET_WGRIB_TEMPLATE("S",IGRIB_HANDLE)
        ELSE
          IGRIB_HANDLE=-99
          CALL IGRIB_CLONE(NGRIB_HANDLE_WAM_S,IGRIB_HANDLE)
        ENDIF
      ELSE
        IF(LGRHDIFS) THEN
          CALL PRESET_WGRIB_TEMPLATE("I",IGRIB_HANDLE)
        ELSE
          IGRIB_HANDLE=-99
          CALL IGRIB_CLONE(NGRIB_HANDLE_WAM_I,IGRIB_HANDLE)
        ENDIF
      ENDIF

      ALLOCATE(VALUES(NTENCODE))

!*    0. PUT FIELD INTO GLOBAL MATRIX VALUES.
!        -----------------------------------
      IF(IPARAM.EQ.251) THEN
        PMISS=0.
      ELSE
        PMISS=ZMISS
      ENDIF

      IF(IRGG.EQ.1 .OR. CLDOMAIN == 'm' ) THEN
        ICOUNT=1
      ELSEIF(CLDOMAIN == 's' ) THEN
        ICOUNT=1
      ELSE
        ICOUNT = (NINT((90. - AMONOP ) / XDELLA))*I1 + 1
        VALUES=PMISS
      ENDIF

!     PAD THE POLES IF INCLUDED IN THE GRID
      IF(LPADPOLES) THEN
        IF((NINT((90. - AMONOP ) / XDELLA)).EQ.0) THEN
          TEMP=0.
          NN=0
          J=2
          JSN=I2-J+1
          DO I=1,NLONRGG(JSN)
            IF(FIELD(I,J).NE.ZMISS) THEN
              TEMP=TEMP+FIELD(I,J)
              NN=NN+1
            ENDIF
          ENDDO
          IF(NN.GT.0) THEN
            TEMP=TEMP/NN
          ELSE
            TEMP=ZMISS
          ENDIF
          J=1
          JSN=I2-J+1
          DO I=1,NLONRGG(JSN)
            FIELD(I,J)=TEMP
          ENDDO
        ENDIF
        IF((NINT((-90. - AMOSOP ) / XDELLA)).EQ.0) THEN
          TEMP=0.
          NN=0
          J=I2-1
          JSN=I2-J+1
          DO I=1,NLONRGG(JSN)
            IF(FIELD(I,J).NE.ZMISS) THEN
              TEMP=TEMP+FIELD(I,J)
              NN=NN+1
            ENDIF
          ENDDO
          IF(NN.GT.0) THEN
            TEMP=TEMP/NN
          ELSE
            TEMP=ZMISS
          ENDIF
          J=I2
          JSN=I2-J+1
          DO I=1,NLONRGG(JSN)
            FIELD(I,J)=TEMP
          ENDDO
        ENDIF
      ENDIF

!     FILL THE ENCODING ARRAY
!     -----------------------

      IF(IPARAM.EQ.251) THEN
        DELTAPP=(2**NGRBRESS-1)*PPRESOL
        ABSPPREC=ABS(PPREC)
        DO J=1,I2
          JSN=I2-J+1
          DO I=1,NLONRGG(JSN)
            IF(FIELD(I,J).NE.ZMISS) THEN
              VALUES(ICOUNT)=FIELD(I,J)
            ELSE
              VALUES(ICOUNT)=PMISS
            ENDIF
            ICOUNT=ICOUNT+1
          ENDDO
        ENDDO

        CALL GSTATS(1433,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JC)
        DO JC=1,NTENCODE
          VALUES(JC) = LOG10(VALUES(JC)+PPEPS)+ABSPPREC
        ENDDO
!$OMP END PARALLEL DO
        PPMAX=VALUES(1)
        DO JC=2,NTENCODE
          PPMAX=MAX(PPMAX,VALUES(JC))
        ENDDO
        PPMIN=MIN(PPMISS,PPMAX-DELTAPP)
        PPMIN=MIN(PPMIN,PPMIN_RESET)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JC)
        DO JC=1,NTENCODE
          IF ( VALUES(JC) .LE. PPMIN ) VALUES(JC)=ZMISS
        ENDDO
!$OMP END PARALLEL DO

      CALL GSTATS(1433,1)

      ELSE
        DO J=1,I2
          JSN=I2-J+1
          DO I=1,NLONRGG(JSN)
            VALUES(ICOUNT)=FIELD(I,J)
            ICOUNT=ICOUNT+1
          ENDDO
        ENDDO
      ENDIF

!*    1. FIX PARAMETERS AND PACK DATA.
!        -----------------------------

!     SPECIFIC VALUES :

!     MISSING DATA:
      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'missingValue',ZMISS)

!     GRIB TABLE AND PARAMETER NUMBER
      ITABPAR=ITABLE*1000+IPARAM
      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'paramId',ITABPAR)

!     LEVEL DEFINITION
      IF(.NOT.LNEWLVTP) THEN
        IF(KLEV.NE.0) THEN
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'levtype',105)
        ELSE
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'levtype',102)
        ENDIF
      ENDIF
      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'level',KLEV)


      IF(.NOT. LGRHDIFS .OR. 
     &   (MARSTYPE .EQ. 'an' .AND. IFCST .EQ. 0) .OR.
     &   (MARSTYPE .EQ. 'fg' .AND. IFCST .EQ. 0) ) THEN

        READ(CDATE(1:8),'(I8)') IDATE
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'date',IDATE)
        READ(CDATE(9:12),'(I4)') ITIME

        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'time',ITIME)

        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'stepUnits','h')
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'stepType','instant')
!!!   for compatibility with previous coding, impose:
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'timeRangeIndicator',10)
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'endStep',IFCST)

        IF ( MARSTYPE .EQ. 'fg' .AND. IFCST .EQ. 0 ) THEN
          ICLASS = 1
        ELSEIF ( MARSTYPE .EQ. '4v' ) THEN
          ICLASS = 6
        ELSEIF ( MARSTYPE .EQ. 'an' .AND. IFCST .EQ. 0 ) THEN
          ICLASS = 2
        ELSEIF ( MARSTYPE .EQ. 'cf' ) THEN
          ICLASS = 10 
        ELSEIF ( MARSTYPE .EQ. 'pf' ) THEN
          ICLASS = 11 
        ELSE
          ICLASS = 9
        ENDIF

        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'type',ICLASS)

        IF (ICLASS.NE.9.AND.ICLASS.NE.10.AND.ICLASS.NE.11
     &     .AND.ICLASS.NE.6.AND.IFCST.GT.0) THEN
          WRITE(IU06,*)' SUB: WGRIBOUT: THIS IS A FORECAST'
          WRITE(IU06,*)' BUT MARSTYPE DOES NOT KNOW ABOUT IT'
          WRITE(IU06,*)'  '
          WRITE(IU06,*)' CALL ABORT1 '
          WRITE(IU06,*)'  '
          CALL ABORT1
        ENDIF

      ELSE
!       TAKE THE COMMON INFORMATION FROM THE IFS
!       MOST OF IT WAS ALREADY OBTAINED WHEN CLONE OF IFS HANDLE WAS TAKEN
!       IN *PRESET_WGRIB_TEMPLATE*

        IF(LRSTST0) THEN
!         NEED TO TEMPORARLY RESET THE IFS FORECAST STEP
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'stepUnits','s')
          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'endStep',ISTEP)

          ISTEP_HRS=ISTEP
          ISTEP=ISTEP-KCOUSTEP
          IF(ISTEP.EQ.0) THEN 
            ISTEP_HRS=0
          ELSE
            ISTEP_HRS=ISTEP_HRS/3600
          ENDIF
        ELSE
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'stepUnits','h')
          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'endStep',ISTEP_HRS)
        ENDIF

!!!   for compatibility with previous coding, impose:
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'timeRangeIndicator',10)

        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'stepUnits','h')
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'endStep',ISTEP_HRS)

      ENDIF

      IF(IPARAM.EQ.251) THEN
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'directionNumber',IK)
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'frequencyNumber',IM)
      ENDIF

!     ENCODE DATA:
      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'values',VALUES)

      DEALLOCATE(VALUES)

! ----------------------------------------------------------------------

!*    2. SAVE ENCODED RESULT
!        -------------------

      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'stream',C12)
      CSTREAM=C12(1:4)
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'expver',C12)
      CEXPVER=C12(1:4)
      WRITE(IU06,*)'  '
      WRITE(IU06,*)'   SUB. WGRIBOUT : PARAMETER ',IPARAM,
     &                              ' EXPVER=',CEXPVER,
     &                              ' STREAM= ', CSTREAM,
     &                              ' TYPE = ',MARSTYPE
      CALl FLUSH(IU06)


      IF(LFDB) THEN
!*     2.1 WRITE DATA TO FDB
!          -----------------

        CALL WGRIB2FDB (IU06, ITEST, IGRIB_HANDLE, 
     &                  CDFDBSF,KFDB, LFDBOPEN,
     &                  IERR)
        IF(ITEST.GT.0) THEN
          WRITE(IU06,*) '   SUB. WGRIBOUT: GRIB DATA WRITTEN TO FDB',
     &     '(REFERENCE=',KFDB,')'
           CALl FLUSH(IU06)
        ENDIF

        IF(IERR.NE.0)THEN
          WRITE(IU06,*) ' ------------------------'
          WRITE(IU06,*) ' ERROR ACCESSING FDB '
          WRITE(IU06,*) ' FDB ERROR CODE IS ',IERR
          WRITE(IU06,*) ' ------------------------'
          CALL ABORT1
        ENDIF

      ELSE

!*      2.2 WRITE PACKED DATA TO FILE CONNECTED TO UNIT IU.
!           ----------------------------------------------

        CALL IGRIB_GET_MESSAGE_SIZE(IGRIB_HANDLE,KBYTES)
        ISIZE=(KBYTES+NPRECI-1)/NPRECI
        ALLOCATE(KGRIB_BUFR(ISIZE))
        CALL IGRIB_GET_MESSAGE(IGRIB_HANDLE,KGRIB_BUFR)
        CALL IGRIB_WRITE_BYTES(IU,KGRIB_BUFR,KBYTES)

        DEALLOCATE(KGRIB_BUFR)

        IF (ITEST.GT.2) THEN
            WRITE(IU06,*) '   SUB. WGRIBOUT : DATA WRITTEN TO FILE'
        ENDIF
      ENDIF

!     RELEASE GRIB HANDLE

      CALL IGRIB_RELEASE(IGRIB_HANDLE)

      CALL GSTATS(1709,1)

      IF (LHOOK) CALL DR_HOOK('WGRIBOUT',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE WGRIBOUT
