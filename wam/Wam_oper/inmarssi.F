      SUBROUTINE INMARSSI (KUSO, IUNIT, IDY, ISX, ISY, PMISS,
     .                     CDATE, IFORP, IPARAM, IRGG, IDOMRGG,
     .                     NGX, NGY, AMOWEP, AMOSOP, AMOEAP, AMONOP,
     .                     NSX, NSY, PX, PY, IPLENG, IGRIB, KRET, MSG)
C
C ======================================================================
C
C***  *INMARSSI* - READS IN GRIB HEADERS FOR  GRIB FORMATTED SPECTRA.
C                  AND RETURNS THE UNPACKED GRIB FIELD.
C
C      J. BIDLOT    ECMWF   JANUARY 1998. 
C
C     PURPOSE.
C     --------
C
C         READS GRIB FORMATTED FILES CONTAINING SPECTRA
C         AND INFERS FROM THIS
C             A) THE DATE OF THE FIELD.
C             B) THE DEFINITION OF THE GRID.
C             C) THE DEFINITION OF THE SPECTRA.
C             D) THE TYPE OF DATA.
C
C**   INTERFACE.
C     ----------
C
C      *CALL INMARSSI* *(KUSO, IUNIT, ISX, ISY,
C                       CDATE, IFORP, IPARAM,
C                       NGX, NGY, AMOWEP, AMOSOP, AMOEAP, AMONOP,
C                       NSX, NSY, PX, PY, IPLENG, IGRIB, KRET, MSG)*
C   I/   *KUSO*   - STANDARD OUTPUT UNIT.
C   I/   *IUNIT*  - DATA INPUT UNIT.
C   I/   *IDY  *  - THE LENGTH OF THE DATA ARRAY 'IDOMRGG'.
C   I/   *ISX  *  - THE LENGTH OF THE MATRIX DESCRIPTION ARRAY PX.
C   I/   *ISY  *  - THE LENGTH OF THE MATRIX DESCRIPTION ARRAY PY.
C   I/   *PMISS*  - REAL VALUE USED TO INDICATE MISSING DATA IN FIELD.
C    /O  *CDATE*  - DATE/TIME OF THE DATA READ (ANALYSIS DATE).
C    /O  *IFORP*  - FORCAST PERIOD IN SECONDS.
C    /O  *IPARAM* - DATA CODE: 250 = SPECTRAL DENSITY (METRES**2*SEC/RAD.)
C    /O  *IRGG*   - GRID CODE: 0 = REGULAR, 1 = IRREGULAR.
C    /O  *IDOMRGG*- NUMBER OF GRID POINTS FOR EACH LATITUDE FOR THE DOMAIN. 
C    /O  *NGX  *  - NUMBER OF COLUMNS IN DATA ARRAY AS FROM MARS.
C    /O  *NGY  *  - NUMBER OF ROWS    IN DATA ARRAY AS FROM MARS.
C    /O  *AMOWEP* - MOST WESTERN LONGITUDE IN GRID (  1, ? ).
C    /O  *AMOSOP* - MOST SOUTHERN LATITUDE IN GRID.( ? ,NGY).
C    /O  *AMOEAP* - MOST EASTERN LONGITUDE IN GRID (NGX, ? ).
C    /O  *AMONOP* - MOST NORTHERN LATITUDE IN GRID ( ? , 1 ).
C    /O  *NSX*    - LENGTH OF MATRIX DESCRIPTION ARRAY PX AS FROM MARS.
C    /O  *NSY*    - LENGTH OF MATRIX DESCRIPTION ARRAY PY AS FROM MARS.
C    /O  *PX*     - MATRIX DESCRIPTION ARRAY PX FIRST INDEX OF MATRIX.
C    /O  *PY*     - MATRIX DESCRIPTION ARRAY PY SEC. INDEX OF MATRIX.
C    /O  *IPLENG* - LENGTH OF THE UNPACKED GRIB DATA ARRAY IGRIB.
C    /O  *IGRIB*  - UNPACKED GRIB DATA ARRAY.
C    /O  *KRET*   - ERROR NUMBER.
C    /O  *MSG*    _ ERROR MESSAGE.
C
C     EXTERNALS.
C     ----------
C
C     EMOSLIB
C
C     *GRIBEX*         CODES/DECODES GRIB FORMATTED DATA.
C     *GRPRS0*         PRINTS SECTION 0 OF GRIB.
C     *GRPRS1*         PRINTS PRODUCT DEFINITION BLOCK.
C     *GRPRS2*         PRINTS GRID  DEFINITION BLOCK.
C
C ----------------------------------------------------------------------
C
      INTEGER   KSEC0(2),
     .          KSEC1(43),
     .          KSEC2(384),
     .          KSEC3(2),
     .          KSEC4(128)
C
      DIMENSION PSEC2(96),
     .          PSEC3(2),
     .          PSEC4(1),
     .          PX(ISX),
     .          PY(ISY),
     .          ZKSEC4(128),
     .          IDOMRGG(IDY)
      EQUIVALENCE (KSEC4(1),ZKSEC4(1))
C
      INTEGER :: IGRIB(IPLENG)

      CHARACTER HOPER*1
      CHARACTER CLMSG*40
      DATA HOPER / 'D' /
      DATA KSEC3 /0,     0 /
      DATA PSEC3 /0., -999./
      CHARACTER*40 ERRMSG(-4:6), MSG
      CHARACTER*10 CDATE
C
      DATA ERRMSG/
     4            ' DECODED WITH BIT MAP                  *',
     3            '                                       *',
     2            '                                       *',
     1            '                                       *',
     1            ' NO ERROR                              *',
     1            ' END OF FILE ENCOUNTED                 *',
     2            ' DECODING ERROR SEE GRIBEX DESCRIPTION *',
     3            ' SUSPICIOUS TIME UNIT IN BLOCK 1       *',
     4            ' DIMENSION 1   IS TOO SMALL            *',
     5            ' NO 2D SPECTRA IN INPUT USE "INMARSB"  *',
     6            ' DIMENSIONS ISX OR ISY ARE TO SMALL    *'/
C
C ----------------------------------------------------------------------
C
#include <parcons.h>
C
C ---------------------------------------------------------------------
C
#include <outctrl.h>
      DATA csubna /"inmarssi"/

C
C ---------------------------------------------------------------------
C
C*    1. INITIALIZATIONS.
C        ----------------
 1000 CONTINUE
      msub = msub + 1
      CALL  uraopc (csubna, msub, llsuvi, lltest)
      IF (llsuvi .OR. lltest) WRITE(kuso,*)' ',csubna,': '
C
      KSEC4=-99
      PSEC3(2) = PMISS
      IF (llsuvi .OR. lltest) WRITE(kuso,*)' ',csubna,': '
C
C ----------------------------------------------------------------------
C
C*    1.0 READ MARS DATA IN GRIB CODE FORMAT.
C         -----------------------------------
C
      CALL PBGRIB(IUNIT,IGRIB,IPLENG,ILENG,IRET)
      IF     (iret.EQ.-1) THEN
        WRITE (kuso,*) ' SUB. INMARSSI : ' 
        WRITE (kuso,*) ' REACHED EOF ON UNIT ',iunit
      ELSEIF (iret .EQ. -2) THEN
        WRITE (kuso,*) ' SUB. INMARSSI : ' 
        WRITE (kuso,*) ' ERROR IN FILE HANDLING ON UNIT ', iunit
      ELSEIF (iret .EQ. -3) THEN
        WRITE (kuso,*) ' SUB. INMARSSI : ' 
        WRITE (kuso,*) ' IGRIB IS NOT BIG ENOUGH. CURRENT SIZE: ',IPLENG
        WRITE (kuso,*) ' REQUIRED SIZE SHOULD BE: ',ILENG
      ENDIF
      IF(iret .LT. -1) CALL abort
C
C ----------------------------------------------------------------------
C
C*    2.0 UNPACK MARS FIELDS.
C         -------------------
C
      KRET  = 1
      ILENP = 1
      IF (lltest) CALL GRSDBG (1)
      CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     1             PSEC4, ILENP,   IGRIB, ileng , KWORD, HOPER, KRET)
      IF (llsuvi) WRITE(kuso,*)' GRIBEX DONE'


      IF(KRET.GT.0) THEN
        MSG = ERRMSG(2)
        CALL GRPRS0 (KSEC0)
        CALL GRPRS1 (KSEC0, KSEC1)
        msub = msub - 1
        WRITE (kuso,*) MSG 
        CALL abort
      ENDIF
      IF(KSEC4(8).NE.64) THEN
        MSG = ERRMSG(5)
        CALL GRPRS0 (KSEC0)
        CALL GRPRS1 (KSEC0, KSEC1)
        CALL GRPRS2 (KSEC0, KSEC1, PSEC2)
        CALL GRPRS3 (KSEC0, KSEC3, PSEC3)
        msub = msub - 1
        RETURN
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    3.0 PRODUCT DEFINITION BLOCK.
C         -------------------------
C
C*    3.1 PRINT BLOCK.
C         ------------
C
      IF (lltest) THEN
        CALL GRPRS0 (KSEC0)
        CALL GRPRS1 (KSEC0, KSEC1)
      ENDIF
C
C*    3.2 DETERMINE DATE FROM KSEC1(10)-KSEC1(16).
C         ----------------------------------------
C
C         START DATE.
C
      WRITE (CDATE,'(5I2.2)') KSEC1(10),KSEC1(11),KSEC1(12),KSEC1(13),
     .                        KSEC1(14)
C
C*         DETERMINE TIME UNIT.
C
      IF (KSEC1(15).EQ.0) THEN
         JCONS = 60
      ELSE IF (KSEC1(15).EQ.1) THEN
         JCONS = 3600
      ELSE IF (KSEC1(15).EQ.2) THEN
         JCONS = 86400
      ELSE
         KRET = 3
         MSG = ERRMSG(KRET)
      ENDIF
C
C         DETERMINE FORECAST PERIOD.
C
      IFORP = JCONS*KSEC1(16)
C
C*    3.3 DETERMINE CODE FOR DATA FIELD TYPE.
C         -----------------------------------
C
      IPARAM = KSEC1(6)
C
C ----------------------------------------------------------------------
C
C*    4.0 GRID DEFINITION BLOCK.
C         ----------------------
C
C*    4.1 PRINT BLOCK.
C         ------------
C
      IF (lltest) CALL GRPRS2 (KSEC0, KSEC2, PSEC2)
C
C*    4.2 DETERMINE GRID AND MATRIX PARAMETERS.
C         -------------------------------------
C
      AMONOP = FLOAT(KSEC2(4)/1000)+0.1*(MOD(KSEC2(4),1000)/100)+
     &         0.01*(MOD(KSEC2(4),100)/10)+0.001*MOD(KSEC2(4),10)

      AMOWEP = FLOAT(KSEC2(5)/1000)+0.1*(MOD(KSEC2(5),1000)/100)+
     &         0.01*(MOD(KSEC2(5),100)/10)+0.001*MOD(KSEC2(5),10)

      AMOSOP = FLOAT(KSEC2(7)/1000)+0.1*(MOD(KSEC2(7),1000)/100)+
     &         0.01*(MOD(KSEC2(7),100)/10)+0.001*MOD(KSEC2(7),10)

      AMOEAP = FLOAT(KSEC2(8)/1000)+0.1*(MOD(KSEC2(8),1000)/100)+
     &         0.01*(MOD(KSEC2(8),100)/10)+0.001*MOD(KSEC2(8),10)

      NSX = KSEC4(50)
      NSY = KSEC4(51)

      NGY = KSEC2(3)

      IRGG = KSEC2(17)

      IF(IRGG.EQ.1) THEN

        NGX = 0
        NTOT=0
        DO J=1,NGY
           JSN=NGY-J+1
           IDOMRGG(JSN) = KSEC2(22+J)
           NGX = MAX(NGX,IDOMRGG(JSN))
           NTOT=NTOT+KSEC2(22+J)
        ENDDO

      ELSE IF(IRGG.EQ.0) THEN
        NGX = KSEC2(2)
        NTOT=NGX*NGY
        IDOMRGG=NGX
      ELSE
        WRITE(kuso,*) 'INMARSSI: REPRESENTATION OF THE FIELD NOT KNOWN'
        CALL abort
      ENDIF

C
      NSXL=NSX
      IF (NSX.GT.ISX) THEN
         KRET = 6
         MSG = ERRMSG(KRET)
         WRITE (kuso,*) MSG 
         msub = msub - 1
         NSXL=ISX
      ENDIF

      NSYL=NSY
      IF (NSY.GT.ISY) THEN
         KRET = 6
         MSG = ERRMSG(KRET)
         WRITE (kuso,*) MSG 
         msub = msub - 1
         NSYL=ISY
      ENDIF
C
      IPX=0
C     ZKSEC4 = TRANSFER(KSEC4, 1, 128)
      DO 4305 J1 = 60, 59+NSXL
        IPX=IPX+1
        PX(IPX) = ZKSEC4(J1)/DEG
 4305 CONTINUE

      IPY=0
      DO 4310 J1 = 60+NSX, 59+NSX+NSYL
        IPY=IPY+1
        PY(IPY) = ZKSEC4(J1)
 4310 CONTINUE

      IF (lltest) THEN
        WRITE(kuso,*) '  J1  KSEC4: '
        DO 4320 J1=1,128
          WRITE(kuso,'(I5,x,I19)') J1, KSEC4(J1)
 4320   CONTINUE
        WRITE(kuso,*) '  J1  PSEC4: '
        DO 4330 J1=1,300
          WRITE(kuso,'(I5,x,f10.4)') J1, PSEC4(J1)
 4330   CONTINUE
        WRITE(kuso,*) '  J1  ZKSEC4: '
        DO 4340 J1=1,100
          WRITE(kuso,'(I5,x,f10.4)') J1, ZKSEC4(J1)
 4340   CONTINUE
      ENDIF
C
C
      MSG = ERRMSG(KRET)
      msub = msub - 1

      RETURN
C
C ----------------------------------------------------------------------
C
C*    6. END OF FILE BUFFER IN.
C        ----------------------
C
 6000 CONTINUE
      KRET = 1
      MSG = ERRMSG(KRET)
      msub = msub - 1
      WRITE (kuso,*) MSG 
      CALL abort
      END
