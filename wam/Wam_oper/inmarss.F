      SUBROUTINE INMARSS (KUSO, IUNIT, IDY, IDL, ISX, ISY, PMISS,
     &                    CDATE, IFORP, IPARAM, IRGG, IDOMRGG,
     &                    NGX, NGY, AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                    NSX, NSY, PSEC4, PX, PY, KRET, MSG)

! ======================================================================

!***  *INMARSS* - READS IN AND UNPACKS GRIB FORMATTED SPECTRA.

!      B. HANSEN    MPIfM/ECMWF   OCTOBER 1992.
!      J. BIDLOT    ECMWF         OCTOBER 1997. add irregular grid features

!     PURPOSE.
!     --------

!         READS GRIB FORMATTED FILES CONTAINING SPECTRA
!         AND INFERS FROM THIS
!             A) THE DATE OF THE FIELD.
!             B) THE DEFINITION OF THE GRID.
!             C) THE DEFINITION OF THE SPECTRA.
!             C) THE SPECTRA IN A ONE DIMENSIONAL ARRAY.
!             D) THE TYPE OF DATA.

!**   INTERFACE.
!     ----------

!      *CALL INMARSS* *(KUSO, IUNIT, IDL, ISX, ISY,
!                       CDATE, IFORP, IPARAM,
!                       NGX, NGY, AMOWEP, AMOSOP, AMOEAP, AMONOP,
!                       NSX, NSY, PSEC4, PX, PY, KRET, MSG)*
!   I/   *KUSO*   - STANDARD OUTPUT UNIT.
!   I/   *IUNIT*  - DATA INPUT UNIT.
!   I/   *IDY  *  - THE LENGTH OF THE DATA ARRAY 'IDOMRGG'.
!   I/   *IDL  *  - THE LENGTH OF THE DATA ARRAY 'PSEC4'.
!   I/   *ISX  *  - THE LENGTH OF THE MATRIX DESCRIPTION ARRAY PX.
!   I/   *ISY  *  - THE LENGTH OF THE MATRIX DESCRIPTION ARRAY PY.
!   I/   *PMISS*  - REAL VALUE USED TO INDICATE MISSING DATA IN FIELD.
!    /O  *CDATE*  - DATE/TIME OF THE DATA READ (ANALYSIS DATE).
!    /O  *IFORP*  - FORCAST PERIOD IN SECONDS.
!    /O  *IPARAM* - DATA CODE: 250 = SPECTRAL DENSITY (METRES**2*SEC/RAD.)
!    /O  *IRGG*   - GRID CODE: 0 = REGULAR, 1 = IRREGULAR.
!    /O  *IDOMRGG*- NUMBER OF GRID POINTS FOR EACH LATITUDE FOR THE DOMAIN. 
!    /O  *NGX  *  - NUMBER OF COLUMNS IN DATA ARRAY 'PSEC4' AS FROM MARS.
!    /O  *NGY  *  - NUMBER OF ROWS    IN DATA ARRAY 'PSEC4' AS FROM MARS.
!    /O  *AMOWEP* - MOST WESTERN LONGITUDE IN GRID (  1, ? ).
!    /O  *AMOSOP* - MOST SOUTHERN LATITUDE IN GRID.( ? ,NGY).
!    /O  *AMOEAP* - MOST EASTERN LONGITUDE IN GRID (NGX, ? ).
!    /O  *AMONOP* - MOST NORTHERN LATITUDE IN GRID ( ? , 1 ).
!    /O  *NSX*    - LENGTH OF MATRIX DESCRIPTION ARRAY PX AS FROM MARS.
!    /O  *NSY*    - LENGTH OF MATRIX DESCRIPTION ARRAY PY AS FROM MARS.
!    /O  *PSEC4*  - UNPACKED DATA.
!    /O  *PX*     - MATRIX DESCRIPTION ARRAY PX FIRST INDEX OF MATRIX.
!    /O  *PY*     - MATRIX DESCRIPTION ARRAY PY SEC. INDEX OF MATRIX.
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
     &          KSEC1(43),
     &          KSEC2(384),
     &          KSEC3(2),
     &          KSEC4(128)

      DIMENSION PSEC2(96),
     &          PSEC3(2),
     &          PSEC4(IDL),
     &          PX(ISX),
     &          PY(ISY),
     &          ZKSEC4(128),
     &          IDOMRGG(IDY)
      EQUIVALENCE (KSEC4(1),ZKSEC4(1))

      INTEGER, ALLOCATABLE :: KGRIB(:)
      INTEGER I4(2)
      CHARACTER HOPER*1
      CHARACTER CLMSG*40
      DATA HOPER / 'D' /
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
     &            ' DIMENSION IDL IS TOO SMALL            *',
     &            ' NO 2D SPECTRA IN INPUT USE "INMARSB"  *',
     &            ' DIMENSIONS ISX OR ISY ARE TO SMALL    *'/

! ---------------------------------------------------------------------

#include <outctrl.h>
      DATA csubna /"inmarss"/


! ---------------------------------------------------------------------

!*    1. INITIALIZATIONS.
!        ----------------
      msub = msub + 1
      CALL  uraopc (csubna, msub, llsuvi, lltest)
      IF (llsuvi .OR. lltest) WRITE(kuso,*)' ',csubna,': '

      KSEC4=-99
      PSEC3(2) = PMISS
      IF (llsuvi .OR. lltest) WRITE(kuso,*)' ',csubna,': '

      I4=1
      NPRECI = KIND(I4)

!     GET SIZE OF KGRIB

      CALL PBSIZE (IUNIT, IPLENG)
      ISIZE=(IPLENG+NPRECI-1)/NPRECI 
      ALLOCATE (KGRIB(ISIZE))


! ----------------------------------------------------------------------

!*    1.0 READ MARS DATA IN GRIB CODE FORMAT.
!         -----------------------------------

      CALL PBGRIB(IUNIT,KGRIB,IPLENG,ILENG,IRET)
      IF     (iret.EQ.-1) THEN
        WRITE (kuso,*) ' SUB. INMARSS : ' 
        WRITE (kuso,*) ' REACHED EOF ON UNIT ',iunit
      ELSEIF (iret .EQ. -2) THEN
        WRITE (kuso,*) ' SUB. INMARSS : ' 
        WRITE (kuso,*) ' ERROR IN FILE HANDLING ON UNIT ', iunit
      ELSEIF (iret .EQ. -3) THEN
        WRITE (kuso,*) ' SUB. INMARSS : ' 
        WRITE (kuso,*) ' KGRIB IS NOT BIG ENOUGH. CURRENT SIZE: ',IPLENG
        WRITE (kuso,*) ' REQUIRED SIZE SHOULD BE: ',ILENG
      ENDIF
      IF(iret .LT. -1) CALL abort

! ----------------------------------------------------------------------

!*    2.0 UNPACK MARS FIELDS.
!         -------------------

      KRET  = 1
      IF (lltest) CALL GRSDBG (1)
      CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &             PSEC4, IDL,   KGRIB, ISIZE , KWORD, HOPER, KRET)
      IF (llsuvi) WRITE(kuso,*)' GRIBEX DONE'

      DEALLOCATE(KGRIB)

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
        CALL GRPRS4 (KSEC0, KSEC4, PSEC4)
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

      ELSEIF(IRGG.EQ.0) THEN
        NGX = KSEC2(2)
        NTOT=NGX*NGY
        IDOMRGG=NGX
      ELSE
        WRITE(kuso,*) 'INMARSS: REPRESENTATION OF THE FIELD NOT KNOWN'
        CALL abort
      ENDIF


      IF ((NSX.GT.ISX).OR.(NSY.GT.ISY)) THEN
        KRET = 6
        MSG = ERRMSG(KRET)
        msub = msub - 1
        RETURN
      ENDIF

      IPX=0
!     ZKSEC4 = TRANSFER(KSEC4, 1, 128)
      DO J1 = 60, 59+NSX
        IPX=IPX+1
        PX(IPX) = ZKSEC4(J1)/DEG
      ENDDO
      IPY=0
      DO J1 = 60+NSX, 59+NSX+NSY
        IPY=IPY+1
        PY(IPY) = ZKSEC4(J1)
      ENDDO
      IF (lltest) THEN
        WRITE(kuso,*) '  J1  KSEC4: '
        DO J1=1,128
          WRITE(kuso,'(I5,x,I19)') J1, KSEC4(J1)
        ENDDO
        WRITE(kuso,*) '  J1  PSEC4: '
        DO J1=1,300
          WRITE(kuso,'(I5,x,f10.4)') J1, PSEC4(J1)
        ENDDO
        WRITE(kuso,*) '  J1  ZKSEC4: '
        DO J1=1,100
          WRITE(kuso,'(I5,x,f10.4)') J1, ZKSEC4(J1)
        ENDDO
      ENDIF

      IF ((NTOT*NSX*NSY .GT. IDL)) THEN
        KRET = 4
        MSG = ERRMSG(KRET)
        msub = msub - 1
        WRITE (kuso,*) ' SUB. INMARSS : ' 
        WRITE (kuso,*) MSG 
        WRITE (kuso,*) 'REQUIRED DIMENSION = ',NTOT*NSX*NSY  
        WRITE (kuso,*) 'IDL                = ',IDL
        CALL abort
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
      END SUBROUTINE INMARSS


