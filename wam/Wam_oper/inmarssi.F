      SUBROUTINE INMARSSI (KUSO, IUNIT, IDY, ISX, ISY, PMISS,
     .                     CDATE, IFORP, IPARAM, IRGG, IDOMRGG, KK, KM,
     .                     NGX, NGY, AMOWEP, AMOSOP, AMOEAP, AMONOP,
     .                     NSX, NSY, PX, PY, IPLENG, ISIZE, IGRIB,
     .                     KRET, MSG)
C
C ======================================================================
C
C***  *INMARSSI* - READS IN GRIB HEADERS FOR GRIB FORMATTED SPECTRA.
C                  AND RETURNS THE NECESSARY INFORMATION AND
C                  THE UNPACKED GRIB FIELD.
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
C      *CALL *INMARSSI* *(KUSO, IUNIT, ISX, ISY,
C                         CDATE, IFORP, IPARAM,
C                         NGX, NGY, AMOWEP, AMOSOP, AMOEAP, AMONOP,
C                         NSX, NSY, PX, PY, ISIZE, IGRIB, KRET, MSG)*
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
C    /O  *KK*     - 1ST INDEX OF THE READ SPECTRA IF PARAM 251, ELSE = NSX
C    /O  *KM*     - 2ND INDEX OF THE READ SPECTRA IF PARAM 251, ELSE = NSY
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
C   I/   *IPLENG* - LENGTH OF THE UNPACKED GRIB DATA ARRAY IGRIB.
C   I/   *ISIZE*  - SIZE IN BYTES OF THE UNPACKED GRIB DATA.
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
     .          KSEC1(49+ISX+ISY),
     .          KSEC2(384),
     .          KSEC3(2),
     .          KSEC4(128)
C
      DIMENSION PSEC2(96),
     .          PSEC3(2),
     .          PSEC4(1),
     .          PX(ISX),
     .          PY(ISY),
     .          IDOMRGG(IDY)
C
      INTEGER :: IGRIB(ISIZE)
      REAL :: ONE,ZTHETA,ZFRE
      CHARACTER HOPER*1
      CHARACTER CLMSG*40
      DATA HOPER / 'J' /
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
      IF (llsuvi .OR. lltest) WRITE(kuso,*)' ',csubna,': isize=', isize
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
      ELSE 
        WRITE (kuso,*) ' SUB. INMARSSI : '
        WRITE (kuso,*) ' Allocated ', ISIZE, ' WORDS, SUFFICIANT FOR ',
     .                 IPLENG, ' bytes for IGRIB. Read ',
     .                 ILENG, ' bytes from file. Status= ', IRET
      ENDIF
      IF (LLTEST.OR.LLSUVI) CALL FLUSH (KUSO)
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
      IF (llsuvi) THEN
         WRITE(kuso,*)' GRIBEX DONE status=' , KRET
         CALL flush (kuso)
      ENDIF
      IF(KRET.GT.0) THEN
        MSG = ERRMSG(2)
        CALL GRPRS0 (KSEC0)
        CALL GRPRS1 (KSEC0, KSEC1)
        msub = msub - 1
        WRITE (kuso,*) MSG 
        CALL abort
      ENDIF
      IF((KSEC4(8).NE.64.AND.KSEC1(6).EQ.250).OR.
     &   (KSEC4(8).NE.0.AND.KSEC1(6).EQ.251)) THEN
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
C
      NGY = KSEC2(3)
C
      IRGG = KSEC2(17)
C
      IF(IRGG.EQ.1) THEN
C
        IF(NGY.GT.IDY) THEN
          WRITE(kuso,*) 'INMARSSI: IDY is not big enough to accomodate' 
          WRITE(kuso,*) '          the value found in grib header (NGY)' 
          call abort
        ENDIF
C
        ISTART=0
        DO WHILE(KSEC2(23+ISTART).EQ.0)
           ISTART=ISTART+1
        ENDDO
C
        ISTOP=0
        DO WHILE(KSEC2(23+NGY-1-ISTOP).EQ.0)
           ISTOP=ISTOP+1
        ENDDO
C
        IC=0
        DO J=1,NGY
          IF(KSEC2(22+J).NE.0) IC=IC+1
        ENDDO
        NGY=IC
C
        NGX = 0
        DO J=1,NGY
           JSN=NGY-J+1
           IDOMRGG(JSN) = KSEC2(22+J+ISTART)
           NGX = MAX(NGX,IDOMRGG(JSN))
        ENDDO
C
        IDUM = KSEC2(4)-ISTART*KSEC2(10)
        AMONOP = FLOAT(IDUM/1000)+0.1*(MOD(IDUM,1000)/100)+
     &           0.01*(MOD(IDUM,100)/10)+0.001*MOD(IDUM,10)
C
        IDUM = KSEC2(7)+ISTOP*KSEC2(10)
        AMOSOP = FLOAT(IDUM/1000)+0.1*(MOD(IDUM,1000)/100)+
     &         0.01*(MOD(IDUM,100)/10)+0.001*MOD(IDUM,10)
C
      ELSE IF(IRGG.EQ.0) THEN
        NGX = KSEC2(2)
        IDOMRGG=NGX
      ELSE
        WRITE(kuso,*) 'INMARSSI: REPRESENTATION OF THE FIELD NOT KNOWN'
        CALL abort
      ENDIF
C
      IF(IPARAM.EQ.250) THEN
        NSX = KSEC4(53)
        NSY = KSEC4(55)
C
        NSXL=NSX
        IF (NSX.GT.ISX) THEN
           KRET = 6
           MSG = ERRMSG(KRET)
           WRITE (kuso,*) MSG 
           msub = msub - 1
           NSXL=ISX
        ENDIF
C
        NSYL=NSY
        IF (NSY.GT.ISY) THEN
           KRET = 6
           MSG = ERRMSG(KRET)
           WRITE (kuso,*) MSG 
           msub = msub - 1
           NSYL=ISY
        ENDIF
C
        KK=NSX
        KM=NSY
        IPX=0
        DO 4305 J1 = 60, 59+NSXL
          IPX=IPX+1
          ZTHETA = TRANSFER(KSEC4(J1),ONE)
          PX(IPX) = ZTHETA/DEG
 4305   CONTINUE
C
        IPY=0
        DO 4310 J1 = 60+NSX, 59+NSX+NSYL
          IPY=IPY+1
          ZFRE = TRANSFER(KSEC4(J1),ONE)
          PY(IPY) = ZFRE 
 4310   CONTINUE
C
      ELSE
        NSX = KSEC1(46)
        NSXL=NSX
        IF (NSX.GT.ISX) THEN
           KRET = 6
           MSG = ERRMSG(KRET)
           WRITE (kuso,*) MSG 
           msub = msub - 1
           NSXL=ISX
        ENDIF
C
        NSY = KSEC1(47)
        NSYL=NSY
        IF (NSY.GT.ISY) THEN
           KRET = 6
           MSG = ERRMSG(KRET)
           WRITE (kuso,*) MSG 
           msub = msub - 1
           NSYL=ISY
        ENDIF
C
        KK=KSEC1(44)
        KM=KSEC1(45)
C
        DO J1=1,NSXL
          PX(J1) = (FLOAT(KSEC1(49+J1))/KSEC1(48))/DEG
        ENDDO
        DO J1=1,NSYL
          PY(J1) = FLOAT(KSEC1(49+NSX+J1))/KSEC1(49) 
        ENDDO
      ENDIF

      IF (lltest) THEN
        WRITE(kuso,*) '  J1  KSEC4: '
        DO 4320 J1=1,128
          WRITE(kuso,'(I5,x,I19)') J1, KSEC4(J1)
 4320   CONTINUE
        WRITE(kuso,*) '  J1  PSEC4: '
        DO 4330 J1=1,300
          WRITE(kuso,'(I5,x,f10.4)') J1, PSEC4(J1)
 4330   CONTINUE
        WRITE(kuso,*) '  J1  KSEC4: '
        DO 4340 J1=1,128
          WRITE(kuso,'(I5,x,f10.4)') J1, KSEC4(J1)
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
