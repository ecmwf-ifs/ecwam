      SUBROUTINE GRIBPAC (IU06,IU,NTGRIB,PSEC4,I1,I2,IPARAM,IY,IM,ID,IH,
     &                    IFCST,KLEV,LFDB,ITEST,MARSTYPE,YCLASS,YEXPVER,
     &                    CDFDBSF,NENSFNB,NTOTENS,Kfdb)
C
C ----------------------------------------------------------------------
C
C****  *GRIBPAC*  PACKS WAM MODEL INTEGRATED FIELDS INTO GRIB CODE.
C
C       L. ZAMBRESKY (GKSS/ECMWF)  JUNE 1987
C       H. GUNTHER   (GKSS/ECMWF)  JUNE 1991   GRIB EDITION 1.
C       A. GUILLAUME (ECMWF)  JUNE 1992 (FDB)
C       J. BIDLOT    ECMWF  JULY 1996 MESSAGE-PASSING
C                           SEPTEMBER 1996 REDUCED GRID
C                           NTGRIB is now used instead of NGYFUL !!!!
C       IU IS LIMITED TO BE 30 OR 32 !!!
C       B. HANSEN    ECMWF  NOVEMBER 1996 REDUCED GRID ALSO FOR
C                                         MEDITEARRANEAN FIELDS.
C          ATTN!! NEW CHARACTER VARIABLE CLDOMAIN DERIVED FROM
C                 PREPROCESSOR DEFINITION region.
C                 SPECIFY 'region' AS FOLLOWS:
C                 -Dregion=\'g\' IF THE DOMAIN SHOULD BE GLOBAL.
C                 -Dregion=\'m\' IF THE DOMAIN SHOULD BE MEDITE.
C
C       J. BIDLOT    ECMWF NOVEMBER 1997 ADD NENSFNB,NTOTENS FOR
C                                        GRIBBING OF ENSEMBLE WAVE DATA.
C       PURPOSE.
C       --------
C         SUBROUTINE PACKS WAVE FIELDS INTO THE WMO GRIB CODE
C         FORMAT AND ARCHIVE INTO FDB, GLOBAL VERSION, PUT MATRIX
C         INTO GLOBAL MATRIX
C
C**    INTERFACE.
C      ----------
C        *CALL* *GRIBPAC (IU06,IU,NTGRIB,PSEC4,I1,I2,IPARAM,IY,IM,ID,IH,
C                         IFCST,LFDB,ITEST,MARSTYPE,YCLASS,YEXPVER,
C                         NENSFNB,NTOTENS,Kfdb)
C          *IU06*    PRINTER OUTPUT UNIT.
C          *IU*      LOGICAL UNIT FOR PACKED FIELD DATA.
C          *NTGRIB*  TOTAL NUMBER OF GRID POINTS PACKED INTO GRIB 
C          *PSEC4*   FIELD TO BE PACKED.
C          *I1*      FIRST DIMENSION OF FIELD.
C          *I2*      SECOND DIMENSION OF FIELD.
C          *IPARAM*  PARAMETER IDENTIFIER.
C          *IY*      YEAR OF LAST ANALYSIS (format yyyy).
C          *IM*      MONTH OF LAST ANALYSIS.
C          *ID*      DAY OF LAST ANALYSIS.
C          *IH*      HOUR OF LAST ANALYSIS.
C          *IFCST*   FORECAST TIME FROM BASE DATE-TIME-GROUP.
C          *KLEV*   - REFERENCE LEVEL IN full METER
C                     (SHOULD BE 0 EXCEPT FOR 233 AND 245)
C          *LFDB*   - LOGICAL SWITCH TO ACTIVATE FDB PROCESSING
C                     .T. TO ACTIVATE, .F. TO IGNORE
C          *ITEST*  - TEST OUTPUT GIVEN IF ITEST GT 2.
C          *MARSTYPE* TYPE OF CURRENT FILED ( FG, AN, OR FC )
C          *YCLASS* - CLASS OF THE CURRENT RUN.
C          *YEXPVER*- EXPERIMENT VERSION NUMBER OF THE CURRENT RUN.
C          *NENSFNB*- ENSEMBLE FORECAST NUMBER
C          *NTOTENS*- TOTAL NUMBER OF ENSEMBLE FORECASTS.
C          *Kfdb*   - Reference to the FDB database.
C
C      METHOD.
C      -------
C        GRIBEX IS CALLED TO PACK THE FIELDS INTO GRIB CODE.
C
C      EXTERNALS.
C      ----------
C       *PBOPEN*         OPEN FILE TO BE READ BY PBGRIB
C       *PBGRIB*         READ GRIB FILE
C       *PBCLOSE*        CLOSE FILE
C        GRIBEX
C        fld2fdb (provided by operation section)
C      REFERENCES.
C      -----------
C        SEE MARS USER GUIDE FOR DOCUMENTATION OF MARS FIELDS.
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "comgrid.h"
C
#include "commap.h"
C
      LOGICAL FRSTIME30,FRSTIME32
C
      CHARACTER*1 CLDOMAIN
      CHARACTER*2 MARSTYPE
      CHARACTER*2 YCLASS
      CHARACTER*4 YEXPVER
      CHARACTER*8 YEXPVERHLP
      CHARACTER*24 FILNM
      CHARACTER*(*) CDFDBSF
      DATA FRSTIME30 / .TRUE. /
      DATA FRSTIME32 / .TRUE. /
      DATA CLDOMAIN  / region /
C
C ----------------------------------------------------------------------
C
#if region=='s'
#else
#include "param_pack.h"
C
#include "parcons.h"
C
C ----------------------------------------------------------------------
C
C*    0. PUT FIELD INTO GLOBAL MATRIX ZWORK.
C        -----------------------------------
C
      ZWORK=ZMISS
      IF(IRGG.EQ.1 .OR. CLDOMAIN == 'm' ) THEN
       INDEX=1
      ELSE
       INDEX = (NINT((90. - AMONOP ) / XDELLA))*NGX + 1
      ENDIF
      DO J=1,NGY
           JSN=NGY-J+1
        DO I=1,NLONRGG(JSN)
           ZWORK(INDEX)=PSEC4(I,J)
           INDEX=INDEX+1
        ENDDO
      ENDDO
C
C*    1. FIX PARAMETERS AND PACK DATA.
C        -----------------------------
C
 1000 CONTINUE
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
      KSEC1(16) = IFCST
C
C FOR REDUCED GRID OUTPUT
C
      IF(IRGG.EQ.1) THEN
        KSEC2(1) = 0
        KSEC2(6) = 128
        KSEC2(9) = 0
        KSEC2(17) = 1
        ISTART=0
        IF (CLDOMAIN == 'g' ) ISTART = NINT((90. - AMONOP ) / XDELLA)
        DO J=1,NGY
           JSN=NGY-J+1
           KSEC2(22+J+ISTART) = NLONRGG(JSN)
        ENDDO
      ENDIF

      KLENP = NTGRIB 
      KSEC4(1) = KLENP
C #elif region=='m'
C  ---------------------------------------------------------------------
C
C*    1. FIX PARAMETERS AND PACK DATA.
C        -----------------------------
C
C 1000 CONTINUE
C      DO 1001 J12=1,I1*I2
C      ZWORK(J12)=PSEC4(J12)
C 1001  CONTINUE
C      KSEC1(6)  = IPARAM
C      KSEC1(10) = IY
C      KSEC1(11) = IM
C      KSEC1(12) = ID
C      KSEC1(13) = IH
C      KSEC1(16) = IFCST
C      KLENP = I1*I2
C      KSEC4(1) = I1*I2
C#endif
C
      KSEC1(24) = 1
      DO 1002 j=25,36
        KSEC1(j)=0
 1002 CONTINUE
      KSEC1(37) = 1
      IF ( YCLASS .EQ. 'od' ) THEN
        KSEC1(38) = 1
      ELSEIF ( YCLASS .EQ. 'rd' ) THEN
        KSEC1(38) = 2
      ELSEIF ( YCLASS .EQ. 'er' ) THEN
        KSEC1(38) = 3
      ELSEIF ( YCLASS .EQ. 'e4' ) THEN
        KSEC1(38) = 5
      ELSE
        WRITE(IU06,*) ' -------------------------'
        WRITE(IU06,*) ' ERROR DETERMINING CLASS  '
        WRITE(IU06,*) ' THE CODE FOUND IS : ',  YCLASS
        WRITE(IU06,*) ' -------------------------'
        CALL ABORT1
      ENDIF
      IF ( MARSTYPE .EQ. 'fg' .AND. IFCST .EQ. 0 ) THEN
        KSEC1(39) = 1
      ELSEIF ( MARSTYPE .EQ. 'an' .AND. IFCST .EQ. 0 ) THEN
        KSEC1(39) = 2
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
      IF ((KSEC1(39).NE.9.AND.KSEC1(39).NE.10.AND.KSEC1(39).NE.11)
     &   .AND.IFCST.GT.0) THEN
        WRITE(IU06,*)' SUB: GRIBPAC: THIS IS A FORECAST'
        WRITE(IU06,*)' BUT MARSTYPE DOES NOT KNOW ABOUT IT'
        WRITE(IU06,*)'  '
        WRITE(IU06,*)' CALL ABORT1 '
        WRITE(IU06,*)'  '
        CALL ABORT1
      ENDIF
      IF(MARSTYPE.EQ.'cf'.OR.MARSTYPE.EQ.'pf') THEN 
        KSEC1(40)=1081
      ELSE
        KSEC1(40)=1045
      ENDIF
      YEXPVERHLP = ' '
C
Cbh   YEXPVERHLP(5:8) = YEXPVER 
Cbh   READ(YEXPVERHLP,'(A8)') KSEC1(41)
C
      KSEC1(41) = 0
      DO ji = 1, 4
        KSEC1(41) = KSEC1(41)*256 + ICHAR(YEXPVER(ji:ji))
      ENDDO
C
      IF(MARSTYPE.EQ.'cf'.OR.MARSTYPE.EQ.'pf') THEN 
        KSEC1(42) = NENSFNB 
        KSEC1(43) = NTOTENS 
      ELSE
        KSEC1(42) = 0
        KSEC1(43) = 0 
      ENDIF

      WRITE(IU06,*)'  '
      WRITE(IU06,*)'    SUB: GRIBPAC - EXPERIMENT VERSION=',YEXPVER,
     .                           '   TYPE OF MARS DATA = ',MARSTYPE
C      
      KRET  = 1
      IF (ITEST .GT. 8 ) CALL GRSDBG (1)
      CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     1             ZWORK, KLENP, KGRIB, KLENG, KWORD, HOPER, KRET)
C
C ----------------------------------------------------------------------
C
C*    2. CHECK ERROR CODE.
C        ----------------
C
      WRITE(IU06,*) '    SUB. GRIBPAC :  GRIBEX STATUS CODE ',KRET,
     .              ' FOR PARAMETER ', IPARAM
 2000 CONTINUE
C 
      IF(KRET.EQ.0) THEN
         IF(LFDB) THEN
C
C*       2.1 WRITE Field to FDB
C            ------------------
C
            IDIR=0
            IFRE=0
            CALL FLD2FDB (IU06, ITEST, KSEC1, KGRIB, KWORD, 'WRITE',
     .                    IERR,
     .                    MARSTYPE, YCLASS, YEXPVER, CDFDBSF, CLDOMAIN,
     .                    IDIR, IFRE , KFDB)
            WRITE(IU06,*) '    SUB. GRIBPAC:  PACKED FIELD LENGTH ',
     .                    KWORD, ' WRITTEN TO FDB ',
     .                    '(REFERENCE=',Kfdb,')'
C
C*       2.2 PRINT ERROR CODE.
C            -----------------
C
            IF(IERR.NE.0)THEN
               WRITE(IU06,*) ' ------------------------'
               WRITE(IU06,*) ' ERROR ACCESSING FDB '
               WRITE(IU06,*) ' FDB ERROR CODE IS ',IERR
               WRITE(IU06,*) ' ------------------------'
               CALL ABORT1
            ENDIF
         ELSE
C
C*      2.1 OPEN FILE UNIT
C           --------------
           IF (FRSTIME30.AND.(IU.EQ.30)) THEN
               WRITE(FILNM,'(A,I2)') 'fort.',IU
               CALL PBOPEN(IUNIT,FILNM,'w',KRET)
               FRSTIME30=.FALSE.
               IF(KRET.LT.0) WRITE(IU06,*)
     1       'PBOPEN (30) IN GRIBPAC FAILED, UNIT ',IU, 'KRET =',KRET
           ELSEIF (FRSTIME32.AND.(IU.EQ.32)) THEN
               WRITE(FILNM,'(A,I2)') 'fort.',IU
               CALL PBOPEN(IUNIT,FILNM,'w',KRET)
               FRSTIME32=.FALSE.
               IF(KRET.LT.0) WRITE(IU06,*)
     1       'PBOPEN (32) IN GRIBPAC FAILED, UNIT ',IU, 'KRET =',KRET
           ELSE
               WRITE(FILNM,'(A,I2)') 'fort.',IU
               CALL PBOPEN(IUNIT,FILNM,'a',KRET)
               IF(KRET.LT.0) WRITE(IU06,*)
     1       'PBOPEN (a) IN GRIBPAC FAILED, UNIT ',IU, 'KRET =',KRET
           ENDIF


C*       2.2 WRITE PACKED DATA ON UNIT IU.
C            -----------------------------
C
            KOUNT=KWORD*4
            CALL PBWRITE(IUNIT,KGRIB,KOUNT,KRET)
            IF(KRET.LT.0) WRITE(IU06,*) 'PBWRITE FAILED, UNIT ',IU,
     & ' CORRUPTED'

            IF (IU.GT.111) THEN
              WRITE(IU06,2200)'    SUB. GRIBPAC:  PACKED FIELD LENGTH ',
     .                        KWORD, '    WRITTEN TO UNIT IU = ',IU
 2200         FORMAT(' ',a,i6,a,a3)
            ELSE
              WRITE(IU06,*) '    SUB. GRIBPAC:  PACKED FIELD LENGTH ',
     .                      KWORD, '    WRITTEN TO UNIT IU = ',IU
            ENDIF
               CALL PBCLOSE(IUNIT,KRET)
               IF(KRET.LT.0) WRITE(IU06,*)
     1       'PBCLOSE IN GRIBPAC FAILED, UNIT ',IU, 'KRET =',KRET
         ENDIF
      ELSE
C
C*    2.3 PRINT ERROR CODE.
C         -----------------
C
         WRITE(IU06,*) ' -----------------------------------------'
         WRITE(IU06,*) ' ERROR PACKING FIELD BY SUB GRIBEX. '
         WRITE(IU06,*) ' GRIBEX ERROR CODE IS ',KRET
         WRITE(IU06,*) ' -----------------------------------------'
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    3. PRINT GRIB SECTION 1 AND 2.
C         --------------------------
C
 3000 CONTINUE
      IF (ITEST.GT.2) THEN
        CALL GRPRS1 (KSEC0, KSEC1)
        CALL GRPRS2 (KSEC0, KSEC2, PSEC2)
      ENDIF
#endif
      RETURN
      END
