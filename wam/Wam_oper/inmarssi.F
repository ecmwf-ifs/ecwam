      SUBROUTINE INMARSSI (KUSO, IUNIT, IDY, ISX, ISY, PMISS,
     &                     CDATE, IFORP, IPARAM, IRGG, IDOMRGG, KK, KM,
     &                     NGX, NGY, AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                     NSX, NSY, PX, PY, IPLENG, ISIZE, IGRIB,
     &                     KRET, MSG)

! ======================================================================

!***  *INMARSSI* - READS IN GRIB HEADERS FOR GRIB FORMATTED SPECTRA.
!                  AND RETURNS THE NECESSARY INFORMATION AND
!                  THE UNPACKED GRIB FIELD.

!      J. BIDLOT    ECMWF   JANUARY 1998. 

!     PURPOSE.
!     --------

!         READS GRIB FORMATTED FILES CONTAINING SPECTRA
!         AND INFERS FROM THIS
!             A) THE DATE OF THE FIELD.
!             B) THE DEFINITION OF THE GRID.
!             C) THE DEFINITION OF THE SPECTRA.
!             D) THE TYPE OF DATA.

!**   INTERFACE.
!     ----------

!      *CALL *INMARSSI* *(KUSO, IUNIT, ISX, ISY,
!                         CDATE, IFORP, IPARAM,
!                         NGX, NGY, AMOWEP, AMOSOP, AMOEAP, AMONOP,
!                         NSX, NSY, PX, PY, ISIZE, IGRIB, KRET, MSG)*
!   I/   *KUSO*   - STANDARD OUTPUT UNIT.
!   I/   *IUNIT*  - DATA INPUT UNIT.
!   I/   *IDY  *  - THE LENGTH OF THE DATA ARRAY 'IDOMRGG'.
!   I/   *ISX  *  - THE LENGTH OF THE MATRIX DESCRIPTION ARRAY PX.
!   I/   *ISY  *  - THE LENGTH OF THE MATRIX DESCRIPTION ARRAY PY.
!   I/   *PMISS*  - REAL VALUE USED TO INDICATE MISSING DATA IN FIELD.
!    /O  *CDATE*  - DATE/TIME OF THE DATA READ (ANALYSIS DATE).
!    /O  *IFORP*  - FORCAST PERIOD IN SECONDS.
!    /O  *IPARAM* - DATA CODE: 250 = SPECTRAL DENSITY (METRES**2*SEC/RAD.)
!    /O  *IRGG*   - GRID CODE: 0 = REGULAR, 1 = IRREGULAR.
!    /O  *IDOMRGG*- NUMBER OF GRID POINTS FOR EACH LATITUDE FOR THE DOMAIN. 
!    /O  *KK*     - 1ST INDEX OF THE READ SPECTRA IF PARAM 251, ELSE = NSX
!    /O  *KM*     - 2ND INDEX OF THE READ SPECTRA IF PARAM 251, ELSE = NSY
!    /O  *NGX  *  - NUMBER OF COLUMNS IN DATA ARRAY AS FROM MARS.
!    /O  *NGY  *  - NUMBER OF ROWS    IN DATA ARRAY AS FROM MARS.
!    /O  *AMOWEP* - MOST WESTERN LONGITUDE IN GRID (  1, ? ).
!    /O  *AMOSOP* - MOST SOUTHERN LATITUDE IN GRID.( ? ,NGY).
!    /O  *AMOEAP* - MOST EASTERN LONGITUDE IN GRID (NGX, ? ).
!    /O  *AMONOP* - MOST NORTHERN LATITUDE IN GRID ( ? , 1 ).
!    /O  *NSX*    - LENGTH OF MATRIX DESCRIPTION ARRAY PX AS FROM MARS.
!    /O  *NSY*    - LENGTH OF MATRIX DESCRIPTION ARRAY PY AS FROM MARS.
!    /O  *PX*     - MATRIX DESCRIPTION ARRAY PX FIRST INDEX OF MATRIX.
!    /O  *PY*     - MATRIX DESCRIPTION ARRAY PY SEC. INDEX OF MATRIX.
!   I/   *IPLENG* - LENGTH OF THE UNPACKED GRIB DATA ARRAY IGRIB.
!   I/   *ISIZE*  - SIZE IN BYTES OF THE UNPACKED GRIB DATA.
!    /O  *IGRIB*  - UNPACKED GRIB DATA ARRAY.
!    /O  *KRET*   - ERROR NUMBER.
!    /O  *MSG*    _ ERROR MESSAGE.

!     EXTERNALS.
!     ----------

!     EMOSLIB

!     *GRIBEX*         CODES/DECODES GRIB FORMATTED DATA.
!     *GRPRS0*         PRINTS SECTION 0 OF GRIB.
!     *GRPRS1*         PRINTS PRODUCT DEFINITION BLOCK.
!     *GRPRS2*         PRINTS GRID  DEFINITION BLOCK.

! ----------------------------------------------------------------------

      USE YOWPCONS , ONLY : DEG

! ----------------------------------------------------------------------

      INTEGER   KSEC0(2),
     &          KSEC1(49+ISX+ISY),
     &          KSEC2(384),
     &          KSEC3(2),
     &          KSEC4(128)

      DIMENSION PSEC2(96),
     &          PSEC3(2),
     &          PSEC4(1),
     &          PX(ISX),
     &          PY(ISY),
     &          IDOMRGG(IDY)

      INTEGER :: IGRIB(ISIZE)
      REAL :: ONE,ZTHETA,ZFRE
      CHARACTER HOPER*1
      CHARACTER CLMSG*40
      DATA HOPER / 'J' /
      DATA KSEC3 /0,     0 /
      DATA PSEC3 /0., -999./
      CHARACTER*40 ERRMSG(-4:6), MSG
      CHARACTER*12 CDATE

      DATA ERRMSG/
     &            ' DECODED WITH BIT MAP                  *',
     &            '                                       *',
     &            '                                       *',
     &            '                                       *',
     &            ' NO ERROR                              *',
     &            ' END OF FILE ENCOUNTED                 *',
     &            ' DECODING ERROR SEE GRIBEX DESCRIPTION *',
     &            ' SUSPICIOUS TIME UNIT IN BLOCK 1       *',
     &            ' DIMENSION 1   IS TOO SMALL            *',
     &            ' NO 2D SPECTRA IN INPUT USE "INMARSB"  *',
     &            ' DIMENSIONS ISX OR ISY ARE TO SMALL    *'/

! ---------------------------------------------------------------------

#include <outctrl.h>
      DATA csubna /"inmarssi"/


! ---------------------------------------------------------------------

!*    1. INITIALIZATIONS.
!        ----------------
      msub = msub + 1
      CALL  uraopc (csubna, msub, llsuvi, lltest)
      IF (llsuvi .OR. lltest) WRITE(kuso,*)' ',csubna,': isize=', isize

      KSEC4=-99
      PSEC3(2) = PMISS
      IF (llsuvi .OR. lltest) WRITE(kuso,*)' ',csubna,': '

! ----------------------------------------------------------------------

!*    1.0 READ MARS DATA IN GRIB CODE FORMAT.
!         -----------------------------------

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
     &                 IPLENG, ' bytes for IGRIB. Read ',
     &                 ILENG, ' bytes from file. Status= ', IRET
      ENDIF
      IF (LLTEST.OR.LLSUVI) CALL FLUSH (KUSO)
      IF(iret .LT. -1) CALL abort

! ----------------------------------------------------------------------

!*    2.0 UNPACK MARS FIELDS.
!         -------------------

      KRET  = 1
      ILENP = 1
      IF (lltest) CALL GRSDBG (1)
      CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &             PSEC4, ILENP,   IGRIB, ileng , KWORD, HOPER, KRET)
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

! ----------------------------------------------------------------------

!*    3.0 PRODUCT DEFINITION BLOCK.
!         -------------------------

!*    3.1 PRINT BLOCK.
!         ------------

      IF (lltest) THEN
        CALL GRPRS0 (KSEC0)
        CALL GRPRS1 (KSEC0, KSEC1)
      ENDIF

!*    3.2 DETERMINE DATE FROM KSEC1(10)-KSEC1(16).
!         ----------------------------------------

!         START DATE.

      IYYYY=(KSEC1(21)-1)*100+KSEC1(10)
      WRITE (CDATE,'(I4.4,4I2.2)') IYYYY,KSEC1(11),KSEC1(12),
     &                             KSEC1(13),KSEC1(14)

!*         DETERMINE TIME UNIT.

      IF (KSEC1(15).EQ.0) THEN
        JCONS = 60
      ELSEIF (KSEC1(15).EQ.1) THEN
        JCONS = 3600
      ELSEIF (KSEC1(15).EQ.2) THEN
        JCONS = 86400
      ELSE
        KRET = 3
        MSG = ERRMSG(KRET)
      ENDIF

!         DETERMINE FORECAST PERIOD.

      IFORP = JCONS*KSEC1(16)

!*    3.3 DETERMINE CODE FOR DATA FIELD TYPE.
!         -----------------------------------

      IPARAM = KSEC1(6)

! ----------------------------------------------------------------------

!*    4.0 GRID DEFINITION BLOCK.
!         ----------------------

!*    4.1 PRINT BLOCK.
!         ------------

      IF (lltest) CALL GRPRS2 (KSEC0, KSEC2, PSEC2)

!*    4.2 DETERMINE GRID AND MATRIX PARAMETERS.
!         -------------------------------------

      AMONOP = REAL(KSEC2(4)/1000)+0.1*(MOD(KSEC2(4),1000)/100)+
     &         0.01*(MOD(KSEC2(4),100)/10)+0.001*MOD(KSEC2(4),10)

      AMOWEP = REAL(KSEC2(5)/1000)+0.1*(MOD(KSEC2(5),1000)/100)+
     &         0.01*(MOD(KSEC2(5),100)/10)+0.001*MOD(KSEC2(5),10)

      AMOSOP = REAL(KSEC2(7)/1000)+0.1*(MOD(KSEC2(7),1000)/100)+
     &         0.01*(MOD(KSEC2(7),100)/10)+0.001*MOD(KSEC2(7),10)

      AMOEAP = REAL(KSEC2(8)/1000)+0.1*(MOD(KSEC2(8),1000)/100)+
     &         0.01*(MOD(KSEC2(8),100)/10)+0.001*MOD(KSEC2(8),10)

      NGY = KSEC2(3)

      IRGG = KSEC2(17)

      IF(IRGG.EQ.1) THEN

        IF(NGY.GT.IDY) THEN
          WRITE(kuso,*) 'INMARSSI: IDY is not big enough to accomodate'
          WRITE(kuso,*) '          the value found in grib header (NGY)'
          call abort
        ENDIF

        ISTART=0
        DO WHILE(KSEC2(23+ISTART).EQ.0)
          ISTART=ISTART+1
        ENDDO

        ISTOP=0
        DO WHILE(KSEC2(23+NGY-1-ISTOP).EQ.0)
          ISTOP=ISTOP+1
        ENDDO

        IC=0
        DO J=1,NGY
          IF(KSEC2(22+J).NE.0) IC=IC+1
        ENDDO
        NGY=IC

        NGX = 0
        DO J=1,NGY
          JSN=NGY-J+1
          IDOMRGG(JSN) = KSEC2(22+J+ISTART)
          NGX = MAX(NGX,IDOMRGG(JSN))
        ENDDO

        IDUM = KSEC2(4)-ISTART*KSEC2(10)
        AMONOP = REAL(IDUM/1000)+0.1*(MOD(IDUM,1000)/100)+
     &           0.01*(MOD(IDUM,100)/10)+0.001*MOD(IDUM,10)

        IDUM = KSEC2(7)+ISTOP*KSEC2(10)
        AMOSOP = REAL(IDUM/1000)+0.1*(MOD(IDUM,1000)/100)+
     &         0.01*(MOD(IDUM,100)/10)+0.001*MOD(IDUM,10)

      ELSEIF(IRGG.EQ.0) THEN
        NGX = KSEC2(2)
        IDOMRGG=NGX
      ELSE
        WRITE(kuso,*) 'INMARSSI: REPRESENTATION OF THE FIELD NOT KNOWN'
        CALL abort
      ENDIF

      IF(IPARAM.EQ.250) THEN
        NSX = KSEC4(53)
        NSY = KSEC4(55)

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

        KK=NSX
        KM=NSY
        IPX=0
        DO J1 = 60, 59+NSXL
          IPX=IPX+1
          ZTHETA = TRANSFER(KSEC4(J1),ONE)
          PX(IPX) = ZTHETA/DEG
        ENDDO

        IPY=0
        DO J1 = 60+NSX, 59+NSX+NSYL
          IPY=IPY+1
          ZFRE = TRANSFER(KSEC4(J1),ONE)
          PY(IPY) = ZFRE 
        ENDDO

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

        NSY = KSEC1(47)
        NSYL=NSY
        IF (NSY.GT.ISY) THEN
          KRET = 6
          MSG = ERRMSG(KRET)
          WRITE (kuso,*) MSG 
          msub = msub - 1
          NSYL=ISY
        ENDIF

        KK=KSEC1(44)
        KM=KSEC1(45)

        DO J1=1,NSXL
          PX(J1) = (REAL(KSEC1(49+J1))/KSEC1(48))/DEG
        ENDDO
        DO J1=1,NSYL
          PY(J1) = REAL(KSEC1(49+NSX+J1))/KSEC1(49) 
        ENDDO
      ENDIF

      IF (lltest) THEN
        WRITE(kuso,*) '  J1  KSEC4: '
        DO J1=1,128
          WRITE(kuso,'(I5,x,I19)') J1, KSEC4(J1)
        ENDDO
        WRITE(kuso,*) '  J1  PSEC4: '
        DO J1=1,300
          WRITE(kuso,'(I5,x,f10.4)') J1, PSEC4(J1)
        ENDDO
        WRITE(kuso,*) '  J1  KSEC4: '
        DO J1=1,128
          WRITE(kuso,'(I5,x,f10.4)') J1, KSEC4(J1)
        ENDDO
      ENDIF


      MSG = ERRMSG(KRET)
      msub = msub - 1

      RETURN

! ----------------------------------------------------------------------

!*    6. END OF FILE BUFFER IN.
!        ----------------------

      KRET = 1
      MSG = ERRMSG(KRET)
      msub = msub - 1
      WRITE (kuso,*) MSG 
      CALL abort
      END SUBROUTINE INMARSSI


