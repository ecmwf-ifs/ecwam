      SUBROUTINE GRIBPAC (IU06,IU,NTGRIB,FIELD,I1,I2,IPARAM,IY,IM,ID,IH,
     &                    IFCST,KLEV,LFDB,ITEST,MARSTYPE,YCLASS,YEXPVER,
     &                    HOPERI, CDFDBSF,NENSFNB,NTOTENS,Kfdb)

! ----------------------------------------------------------------------

!****  *GRIBPAC*  PACKS WAM MODEL INTEGRATED FIELDS INTO GRIB CODE.

!       L. ZAMBRESKY (GKSS/ECMWF)  JUNE 1987
!       H. GUNTHER   (GKSS/ECMWF)  JUNE 1991   GRIB EDITION 1.
!       A. GUILLAUME (ECMWF)  JUNE 1992 (FDB)
!       J. BIDLOT    ECMWF  JULY 1996 MESSAGE-PASSING
!                           SEPTEMBER 1996 REDUCED GRID
!                           NTGRIB is now used instead of NGYFUL !!!!
!       IU IS LIMITED TO BE 30 OR 32 !!!
!       B. HANSEN    ECMWF  NOVEMBER 1996 REDUCED GRID ALSO FOR
!                                         MEDITEARRANEAN FIELDS.
!       J. BIDLOT    ECMWF NOVEMBER 1997 ADD NENSFNB,NTOTENS FOR
!                                        GRIBBING OF ENSEMBLE WAVE DATA.
!       J. BIDLOT    ECMWF NOVEMBER 1998 USE DEFAULT GRIB HEADER VALUES

!       PURPOSE.
!       --------
!         SUBROUTINE PACKS WAVE FIELDS INTO THE WMO GRIB CODE
!         FORMAT AND ARCHIVE INTO FDB, GLOBAL VERSION, PUT MATRIX
!         INTO GLOBAL MATRIX

!**    INTERFACE.
!      ----------
!        *CALL* *GRIBPAC (IU06,IU,NTGRIB,FIELD,I1,I2,IPARAM,IY,IM,ID,IH,
!                         IFCST,LFDB,ITEST,MARSTYPE,YCLASS,YEXPVER,
!                         NENSFNB,NTOTENS,Kfdb)
!          *IU06*    PRINTER OUTPUT UNIT.
!          *IU*      LOGICAL UNIT FOR PACKED FIELD DATA.
!          *NTGRIB*  TOTAL NUMBER OF GRID POINTS PACKED INTO GRIB 
!          *FIELD*   FIELD TO BE PACKED.
!          *I1*      FIRST DIMENSION OF FIELD.
!          *I2*      SECOND DIMENSION OF FIELD.
!          *IPARAM*  PARAMETER IDENTIFIER.
!          *IY*      YEAR OF LAST ANALYSIS (format yyyy).
!          *IM*      MONTH OF LAST ANALYSIS.
!          *ID*      DAY OF LAST ANALYSIS.
!          *IH*      HOUR OF LAST ANALYSIS.
!          *IFCST*   FORECAST TIME FROM BASE DATE-TIME-GROUP.
!                    SPECIFIED IN HOURS.
!          *KLEV*   - REFERENCE LEVEL IN full METER
!                     (SHOULD BE 0 EXCEPT FOR 233 AND 245)
!          *LFDB*   - LOGICAL SWITCH TO ACTIVATE FDB PROCESSING
!                     .T. TO ACTIVATE, .F. TO IGNORE
!          *ITEST*  - TEST OUTPUT GIVEN IF ITEST GT 2.
!          *MARSTYPE* TYPE OF CURRENT FIELD
!          *YCLASS* - CLASS OF THE CURRENT RUN.
!          *YEXPVER*- EXPERIMENT VERSION NUMBER OF THE CURRENT RUN.
!          *NENSFNB*- ENSEMBLE FORECAST NUMBER
!          *NTOTENS*- TOTAL NUMBER OF ENSEMBLE FORECASTS.
!          *Kfdb*   - Reference to the FDB database.

!      METHOD.
!      -------
!        GRIBEX IS CALLED TO PACK THE FIELDS INTO GRIB CODE.

!      EXTERNALS.
!      ----------
!       *ABORT1*         ABORTS ON ERROR FROM GRIBEX
!       *PBOPEN*         OPEN FILE TO BE READ BY PBGRIB
!       *PBGRIB*         READ GRIB FILE
!       *PBCLOSE*        CLOSE FILE
!        GRIBEX
!        fld2fdb (provided by operation section)
!      REFERENCES.
!      -----------
!        SEE MARS USER GUIDE FOR DOCUMENTATION OF MARS FIELDS.

! ----------------------------------------------------------------------

      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3   ,
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,KSEK0    ,KSEK1    ,
     &            KSEK2    ,KSEK3    ,KSEK4    ,PSEK2    ,PSEK3
      USE YOWGRID  , ONLY : NLONRGG
      USE YOWMAP   , ONLY : IRGG     ,AMONOP   ,XDELLA
      USE YOWPARAM , ONLY : CLDOMAIN
      USE YOWPCONS , ONLY : G        ,ZMISS
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH

! ----------------------------------------------------------------------

      INTEGER KSEC0(NKSEK0), KSEC1(NKSEK1), KSEC2(NKSEK2),
     &        KSEC3(NKSEK3), KSEC4(NKSEK4)
      INTEGER, ALLOCATABLE :: KGRIB(:)

      REAL PSEC2(NPSEK2), PSEC3(NPSEK3), FIELD(I1,I2), ZWORK(NTGRIB)
      LOGICAL LFDB, FRSTIME30, FRSTIME32
      CHARACTER*2 MARSTYPE
      CHARACTER*2 YCLASS
      CHARACTER*4 YEXPVER, CSTREAM
      CHARACTER*24 FILNM
      CHARACTER*(*) CDFDBSF, HOPERI
      DATA FRSTIME30 / .TRUE. /
      DATA FRSTIME32 / .TRUE. /

! ----------------------------------------------------------------------
      CALL GSTATS(1709,0)

!*    0. PUT FIELD INTO GLOBAL MATRIX ZWORK.
!        -----------------------------------

      IF(IRGG.EQ.1 .OR. CLDOMAIN == 'm' ) THEN
        INDEX=1
      ELSEIF(CLDOMAIN == 's' ) THEN
        INDEX=1
      ELSE
        INDEX = (NINT((90. - AMONOP ) / XDELLA))*I1 + 1
        ZWORK=ZMISS
      ENDIF
      DO J=1,I2
        JSN=I2-J+1
        DO I=1,NLONRGG(JSN)
          ZWORK(INDEX)=FIELD(I,J)
          INDEX=INDEX+1
        ENDDO
      ENDDO

!*    1. FIX PARAMETERS AND PACK DATA.
!        -----------------------------


!     DEFAULT VALUES :

      KSEC0 = KSEK0
      KSEC1 = KSEK1
      KSEC2 = KSEK2
      KSEC3 = KSEK3
      KSEC4 = KSEK4

      PSEC2 = PSEK2
      PSEC3 = PSEK3

!     SPECIFIC VALUES :

      KSEC1(6)  = IPARAM
      IF(KLEV.NE.0) THEN
        KSEC1(7)=105
      ELSE
        KSEC1(7)=102
      ENDIF
      KSEC1( 8) = KLEV
      KSEC1(10) = MOD(IY,100) 
      IF(KSEC1(10).EQ.0) THEN
        KSEC1(21) = IY/100
        KSEC1(10) = 100
      ELSE
        KSEC1(21) = (IY/100)+1
      ENDIF
      KSEC1(11) = IM
      KSEC1(12) = ID
      KSEC1(13) = IH
      IF(IFCST >= (2**16-1)) THEN
!     the number of hours is too large, use coding in days or
!     multiple of hours instead
        IF(MOD(IFCST,24) == 0 ) THEN
          KSEC1(15)= 2 
          KSEC1(16) = IFCST/24
        ELSEIF(MOD(IFCST,12) == 0 ) THEN
          KSEC1(15)= 12 
          KSEC1(16) = IFCST/12
        ELSEIF(MOD(IFCST,6) == 0 ) THEN
          KSEC1(15)= 11 
          KSEC1(16) = IFCST/6
        ELSEIF(MOD(IFCST,3) == 0 ) THEN
          KSEC1(15)= 10 
          KSEC1(16) = IFCST/3
        ELSE
          WRITE(IU06,*)' SUB: GRIBPAC: THE FORECAST IS TOO LONG'
          WRITE(IU06,*)' IT CANNOT BE CODED IN GRIB'
          WRITE(IU06,*)' IFCST = ',IFCST
          WRITE(IU06,*)' CALL ABORT1 '
          WRITE(IU06,*)'  '
        ENDIF 
      ELSE
        KSEC1(15)= 1 
        KSEC1(16) = IFCST
      ENDIF


      IF ( MARSTYPE .EQ. 'fg' .AND. IFCST .EQ. 0 ) THEN
        KSEC1(39) = 1
      ELSEIF ( MARSTYPE .EQ. '4v' ) THEN
        KSEC1(39) = 6
      ELSEIF ( MARSTYPE .EQ. 'an' .AND. IFCST .EQ. 0 ) THEN
        KSEC1(39) = 2
      ELSEIF ( MARSTYPE .EQ. 'sf' ) THEN
        IF(NTOTENS.EQ.0) THEN
          KSEC1(39) = 9
        ELSE
          IF(NENSFNB.EQ.0) THEN
            KSEC1(39) = 10 
          ELSE
            KSEC1(39) = 9 
          ENDIF
        ENDIF
      ELSE
        IF(NTOTENS.EQ.0) THEN
          KSEC1(39) = 9
        ELSE
          IF(NENSFNB.EQ.0) THEN
            KSEC1(39) = 10
          ELSE
            KSEC1(39) = 11
          ENDIF
        ENDIF
      ENDIF
      IF (KSEC1(39).NE.9.AND.KSEC1(39).NE.10.AND.KSEC1(39).NE.11
     &   .AND.KSEC1(39).NE.6.AND.IFCST.GT.0) THEN
        WRITE(IU06,*)' SUB: GRIBPAC: THIS IS A FORECAST'
        WRITE(IU06,*)' BUT MARSTYPE DOES NOT KNOW ABOUT IT'
        WRITE(IU06,*)'  '
        WRITE(IU06,*)' CALL ABORT1 '
        WRITE(IU06,*)'  '
        CALL ABORT1
      ENDIF

      KLENP = NTGRIB 
      KSEC4(1) = KLENP

      KLENG =NINT(NTGRIB*(KSEC4(2)+1.)/(KIND(IU06)*8))+NKSEK0+NKSEK1+
     &  NKSEK2+NKSEK3+NKSEK4+NPSEK2+NPSEK3 

      ALLOCATE(KGRIB(KLENG))

      KRET  = 1
      IF (ITEST .GT. 8 ) CALL GRSDBG (1)
      CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &             ZWORK, KLENP, KGRIB, KLENG, KWORD, HOPERI, KRET)

! ----------------------------------------------------------------------

!*    2. CHECK ERROR CODE.
!        ----------------

      CALL WSTREAM_STRG(KSEC1(40),CSTREAM)

      WRITE(IU06,*)'  '
      WRITE(IU06,*)'   SUB. GRIBPAC : PARAMETER ',IPARAM,
     &                              ' STREAM= ', CSTREAM,
     &                              ' EXPVER=',YEXPVER,
     &                              ' TYPE OF MARS DATA = ',MARSTYPE

      IF(KRET.EQ.0) THEN
        IF(LFDB) THEN

!*       2.1 WRITE Field to FDB
!            ------------------

          IDIR=0
          IFRE=0
          CALL FLD2FDB (IU06, ITEST, KSEC1, KGRIB, KWORD, 'WRITE',
     &     IERR,
     &     MARSTYPE, YCLASS, YEXPVER, CDFDBSF, CLDOMAIN,
     &     IDIR, IFRE , KFDB)
          IF(ITEST.GT.0) THEN
            WRITE(IU06,*) '   SUB. GRIBPAC:  GRIB FIELD LENGTH ',
     &       KWORD, ' WRITTEN TO FDB ',
     &       '(REFERENCE=',Kfdb,')'
          ENDIF

!*       2.2 PRINT ERROR CODE.
!            -----------------

          IF(IERR.NE.0)THEN
            WRITE(IU06,*) ' ------------------------'
            WRITE(IU06,*) ' ERROR ACCESSING FDB '
            WRITE(IU06,*) ' FDB ERROR CODE IS ',IERR
            WRITE(IU06,*) ' ------------------------'
            CALL ABORT1
          ENDIF
        ELSE

!*      2.1 OPEN FILE UNIT
!           --------------
          IF (FRSTIME30.AND.(IU.EQ.30)) THEN
            WRITE(FILNM,'(A,I2)') 'fort.',IU
            CALL PBOPEN(IUNIT,CPATH(1:ICPLEN)//'/'//FILNM,'w',KRET)
            FRSTIME30=.FALSE.
            IF(KRET.LT.0) WRITE(IU06,*)
     &       'PBOPEN (30) IN GRIBPAC FAILED, UNIT ',IU, 'KRET =',KRET
          ELSEIF (FRSTIME32.AND.(IU.EQ.32)) THEN
            WRITE(FILNM,'(A,I2)') 'fort.',IU
            CALL PBOPEN(IUNIT,CPATH(1:ICPLEN)//'/'//FILNM,'w',KRET)
            FRSTIME32=.FALSE.
            IF(KRET.LT.0) WRITE(IU06,*)
     &       'PBOPEN (32) IN GRIBPAC FAILED, UNIT ',IU, 'KRET =',KRET
          ELSE
            WRITE(FILNM,'(A,I2)') 'fort.',IU
            CALL PBOPEN(IUNIT,CPATH(1:ICPLEN)//'/'//FILNM,'a',KRET)
            IF(KRET.LT.0) WRITE(IU06,*)
     &       'PBOPEN (a) IN GRIBPAC FAILED, UNIT ',IU, 'KRET =',KRET
          ENDIF


!*       2.2 WRITE PACKED DATA ON UNIT IU.
!            -----------------------------

          KOUNT=KWORD*4
          CALL PBWRITE(IUNIT,KGRIB,KOUNT,KRET)
          IF(KRET.LT.0) WRITE(IU06,*) 'PBWRITE FAILED, UNIT ',IU,
     &     ' CORRUPTED'

          IF (ITEST.GT.2) THEN
            IF (IU.GT.111) THEN
              WRITE(IU06,2200)'   SUB. GRIBPAC : PACKED FIELD LENGTH ',
     &         KWORD, '    WRITTEN TO UNIT IU = ',IU
 2200         FORMAT(' ',a,i6,a,a3)
            ELSE
              WRITE(IU06,*) '   SUB. GRIBPAC : PACKED FIELD LENGTH ',
     &         KWORD, '    WRITTEN TO UNIT IU = ',IU
            ENDIF
          ENDIF
          CALL PBCLOSE(IUNIT,KRET)
          IF(KRET.LT.0) WRITE(IU06,*)
     &     'PBCLOSE IN GRIBPAC FAILED, UNIT ',IU, 'KRET =',KRET
        ENDIF
      ELSE

!*    2.3 PRINT ERROR CODE.
!         -----------------

        WRITE(IU06,*) ' -----------------------------------------'
        WRITE(IU06,*) ' ERROR PACKING FIELD BY SUB GRIBEX. '
        WRITE(IU06,*) ' PARAMETER : ',IPARAM 
        WRITE(IU06,*) ' GRIBEX ERROR CODE IS ',KRET
        WRITE(IU06,*) ' -----------------------------------------'
        WRITE(*,*) ' -----------------------------------------'
        WRITE(*,*) ' ERROR PACKING FIELD BY SUB GRIBEX. '
        WRITE(*,*) ' PARAMETER : ',IPARAM 
        WRITE(*,*) ' GRIBEX ERROR CODE IS ',KRET
        WRITE(*,*) ' -----------------------------------------'
        CALL ABORT1
      ENDIF

      DEALLOCATE(KGRIB)

! ----------------------------------------------------------------------

!*    3. PRINT GRIB SECTION 1 AND 2.
!         --------------------------

      IF (ITEST.GT.2) THEN
        CALL GRPRS1 (KSEC0, KSEC1)
        CALL GRPRS2 (KSEC0, KSEC2, PSEC2)
      ENDIF

      CALL GSTATS(1709,1)

      RETURN
      END SUBROUTINE GRIBPAC
