C
      SUBROUTINE SPEC2FDB (KU06, KPROC, KRANK, LDWCOU, PTHETA, PFREQ,
     .                     PSPEC, KINF, KSUP, KANG, KFRE, KGX, KGY,
     .                     KRGG, KLONRGG, KXLT, KXLG, KIBLO, KBLO,
     .                     PMISS, PMONOP, PMOWEP, PMOEAP, PMOSOP,
     .                     PDELLA, PDELLO, KSTART, KEND,
     .                     CDTPRO, CDATEF, CDATER, CDATES, CDATEA,
     .                     CDTYPE, CDCLASS, CDEXPVER, CDFDB2DSP,
     .                     KENSFNB, KTOTENS, LFDBIOOUT, KTEST, KFDB )
C
C----------------------------------------------------------------------
C
C**** *SPEC2FDB*  ENCODES SPECTRA AS PARAMETER 251
C                 AND WRITES TO FILE OR TO FDB.
C
C     B. HANSEN   ECMWF   FEBRUARY 1998 
C
C*    PURPOSE.
C     --------
C
C       ENCODES SPECTRA AND WRITES TO FDB OR TO A FILE.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *SPEC2FDB* (KU06, PSPEC, KINF, KSUP, KANG, KFRE,
C                          KSTART, KEND,
C                          CDTPRO, CDATEF, CDATER, CDATES,
C                          CDTYPE, CDCLASS, CDEXPVER, CDFDB2DSP,
C                          KENSFNB, KTOTENS, LFDBIOOUT, KTEST, KFDB )
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C      *KU06*      INTEGER   OUTPUT UNIT FOR STANDARD OUTPUT.
C      *KPROC*     INTEGER   TOTAL NUMBER OF PES.
C      *KRANK*     INTEGER   CURRENT PE NUMBER.
C      *LDWCOU*    LOGICAL   TRUE: THIS IS A COUPLED RUN.
C      *PTHETA*    REAL      DIRECTION VECTOR.
C      *PFREQ*     REAL      FREQUENCY VECTOR.
C      *PSPEC*     REAL      LOCAL SPECTRA OF CURRENT PE.
C      *KINF*      INTEGER   SMALLEST GRID POINT INDEX USED BY THE PE.
C      *KSUP*      INTEGER   LARGEST GRID POINT INDEX USED BY THE PE.
C      *KANG*      INTEGER   NUMBER OF ANGLES.
C      *KFRE*      INTEGER   NUMBER OF FREQUENCIES.
C      *KGX*       INTEGER   NUMBER OF LONGITUDES IN GRID.
C      *KGY*       INTEGER   NUMBER OF LATITUDES IN GRID.
C      *KRGG*      INTEGER   GRID CODE: 0 = REGULAR, 1 = IRREGULAR.
C      *KLONRGG*   INTEGER   NUMBER OF POINTS ALONG A MERIDIAN.
C      *KXLT*      INTEGER   LAT. GRID INDEX.
C      *KXLG*      INTEGER   LON. GRID INDEX.
C      *KIBLO*     INTEGER   NUMBER OF POINTS IN BLOCK.
C      *KBLO*      INTEGER   NUMBER OF BLOCKS.
C      *PMISS*     REAL      MISSING DATA VALUE.
C      *PMONOP*    REAL      MOST NORTHERN LATITUDE IN GRID (DEGREE).
C      *PMOWEP*    REAL      MOST WESTERN LATITUDE IN GRID (DEGREE).
C      *PMOEAP*    REAL      MOST EASTERN LATITUDE IN GRID (DEGREE).
C      *PMOSOP*    REAL      MOST SOUTHERN LATITUDE IN GRID (DEGREE).
C      *PDELLA*    REAL      DISTANCE BETWEEN LATITUDES.
C      *PDELLO*    REAL      DISTANCE BETWEEN LONGITUDES.
C      *KSTART*    INTEGER   FIRST INDICES OF LOCAL BLOCKS OF SPECTRA.
C      *KEND*      INTEGER   LAST INDICES OF LOCAL BLOCKS OF SPECTRA.
C      *CDTPRO*    CHAR*10   END DATE OF PROPAGATION.
C      *CDATEF*    CHAR*10   END DATE OF ANALYSIS RUN (YYMMDDHHMM).
C      *CDATER*    CHAR*10   DATE FOR OUTPUT OF BOTH RESTART FILES
C      *CDATES*    CHAR*10   LAST DATE FOR OUTPUT OF RESTART FILES
C                            IF NOT SET IT WILL BE SET BY DEFAULT TO
C                            CDATEE
C      *CDATEA*    CHAR*10   BEGIN DATE OF ANALYSIS RUN (YYMMDDHHMM).
C      *CDTYPE*    CHAR*     MARS TYPE OF FIELD.
C      *CDCLASS*   CHAR*     MARS CLASS OF FIELD.
C      *CDEXPVER*  CHAR*     EXPERIMENT VERSION IDENTIFIER.
C      *CDFDB2DSP*   CHAR*     ROOT DIRECTORY FOR 2D GRIB SPECTRA.
C      *KENSFNB*   INTEGER   ENSEMBLE FORECAST NUMBER
C      *KTOTENS*   INTEGER   TOTAL NUMBER OF ENSEMBLE FORECASTS.
C      *LFDBIOOUT* LOGICAL   IF TRUE OUTPUT TO FDB OTHERWISE TO FILE
C      *KTEST*     INTEGER   PRINT DIAGNOSTICS TO STDOUT IF GT 1.
C
C     METHOD.
C     -------
C
C	   ENCODE SPECTRA PER FREQUENCY AND DIRECTION
C      INTO GRIB AND WRITE TO FDB OR TO A SINGLE FILE.
C
C     EXTERNALS.
C     ----------
C
C      *ABORT1*
C      *GRIBEX* 
C      *FLD2FDB*
C      *MPGATHERSCFLD*
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C-------------------------------------------------------------------
C
      IMPLICIT NONE
C
C*  DUMMY ARGUMENTS:
C   ----------------
C
      LOGICAL LDWCOU, LFDBIOOUT, LOUTALL
C
      INTEGER KU06, KINF, KSUP, KANG, KFRE, KGX, KGY, KRGG, KIBLO, KBLO,
     .        KRANK, KPROC, KSEND, IK, IM, KENSFNB, KTOTENS,
     .        KLONRGG(KGY), LFILE, IUOUT, KOUNT, KTEST, KFDB, KST,
     .        KXLT(KIBLO, KBLO), KXLG(KIBLO, KBLO), KDEXN, KDEXS
      INTEGER, DIMENSION(KPROC) :: KSTART, KEND 
C
      REAL    PMISS, PMONOP, PDELLA, PDELLO, PMOEAP, PMOSOP, PMOWEP,
     .        PSPEC (KINF-1:KSUP, KANG, KFRE),
     .        PTHETA(KANG),
     .        PFREQ(KFRE)
C
      CHARACTER*1 CMODE
      CHARACTER*(*) CDTPRO, CDATEF, CDATER, CDATES, CDATEA,
     .              CDTYPE, CDCLASS, CDEXPVER, CDFDB2DSP
C
C-----------------------------------------------------------------------
C
      INTEGER, ALLOCATABLE :: IOUTGRIB(:), IPOINTER(:)
C
      INTEGER ISEC0(2), ISEC1(49+KANG+KFRE), ISEC2(384), ISEC3(2),
     .        ISEC4(128), JSN, JK, ITOT, JC, IX, ISTEP, IST,
     .        IY, IMO, ID, IH, IMN, IFCST, itheta(KANG), ifreq(KFRE), 
     .        J, IJ, JANG, J1, J2, J3, JI, ILENGOUT, IC,
     .        IRET, IWORD, IERR, ILENP, ITAG, ISTAT

C
      REAL, ALLOCATABLE :: ZSEC4OUT(:),WORK(:)
      REAL ZSEC2(96), ZSEC3(2), ZMONOP, ZMOSOP
      REAL , DIMENSION(KANG) :: ZTHETA
      REAL , DIMENSION(KFRE) :: ZFREQ
C
      REAL, PARAMETER :: PPDEG = 180.0 / 3.1415927
      REAL, PARAMETER :: PPMISS = -3.0
      REAL, PARAMETER :: PPEPS  = 1.0e-10
      REAL, PARAMETER :: PPREC  = 0.0
C
      CHARACTER*  1 CLDOMAIN
      CHARACTER* 40 CLMSG
      CHARACTER* 80 OUTFILEN
C
      DATA CLDOMAIN  / region /
C
      LOGICAL LLFIRST
      DATA LLFIRST  / .TRUE. /
C
C-----------------------------------------------------------------------
C
C*    0. FIXED SETTINGS FOR THE GRIB HEADERS.
C     ---------------------------------------
C
      IF (KTEST.GT.1) THEN
        WRITE(KU06,*) '*********************************'
        WRITE(KU06,*) '*      THIS IS SPEC2FDB         *'
        WRITE(KU06,*) '*********************************'
        CALL FLUSH (KU06)
      ENDIF
      ISEC1 = 0
      ISEC2 = 0
      ISEC4 = 0
      ZSEC2 = 0.0
C
      IF(CDTPRO.LE.CDATEF) THEN
C
C*    0.1.  THIS IS AN ANALYSIS DATE.
C           -------------------------
C
        IF (LDWCOU .AND. CDTYPE .EQ. 'fg') THEN
          READ (CDATEA, '(5I2)') IY, IMO, ID, IH, IMN
          CALL DIFDATE (CDATEA, CDTPRO, IFCST)
          IFCST = IFCST/3600
        ELSE
          READ (CDTPRO, '(5I2)') IY, IMO, ID, IH, IMN
          IFCST = 0
        ENDIF
      ELSE
C
C*    0.2.  THIS IS A  FORECAST DATE.
C           -------------------------
C
        READ (CDATEF, '(5I2)') IY, IMO, ID, IH, IMN
        CALL DIFDATE (CDATEF, CDTPRO, IFCST)
        IFCST = IFCST/3600
      ENDIF
C
C    *SECTION 1*  INTEGER  PRODUCT DEFINITION BLOCK.
C    -----------------------------------------------
C       WORD        VALUE   CONTENTS.
C       ----        -----   ---------
      ISEC1( 1) =   140     ! VERSION NUMBER OF CODE TABLE 2.
      ISEC1( 2) =    98     ! ORIGINATING CENTRE IDENTIFIER.
      ISEC1( 3) =   107     ! MODEL IDENTIFICATION.
      ISEC1( 4) =   255     ! GRID DEFINITION.
      ISEC1( 5) =   192     ! FLAG (WMO CODE TABLE 1)
      ISEC1( 6) =   251     ! PARAMETER IDENTIFIER (WMO CODE TABLE 2)
      ISEC1( 7) =   102     ! TYPE OF LEVEL (WMO CODE TABLE 3)
      ISEC1( 8) =     0     ! TYPE OF LEVEL (WMO CODE TABLE 3)
      ISEC1( 9) =     0     ! VALUE(S) OF LEVEL (WMO CODE TABLE 3)
      ISEC1(10) =   IY      ! YEAR OF DATA
      ISEC1(11) =   IMO     ! MONTH OF DATA
      ISEC1(12) =   ID      ! DAY OF DATA
      ISEC1(13) =   IH      ! HOUR OF DATA
      ISEC1(14) =   IMN     ! MINUTE OF DATA
      ISEC1(15) =     1     ! TIME UNIT (WMO CODE TABLE 4)
      ISEC1(16) =  IFCST    ! TIME RANGE ONE (FCST TIME FROM BASE DTG)
      ISEC1(17) =     0     ! TIME RANGE TWO
      ISEC1(18) =     0     ! TIME RANGE FLAG (WMO CODE TABLE 5)
      ISEC1(19) =     0     ! NUMBER AVERAGED.
      ISEC1(20) =     0     ! NUMBER MISSING FROM AVERAGE.
      ISEC1(21) =    20     ! CENTURY OF REFERENCE TIME.
      ISEC1(22) =     0     ! 
      ISEC1(23) =     0     !
      ISEC1(24) =     1     ! INDICATES LOCAL USE OF SECTION 1.
C       25-36                                                  0
      ISEC1(37) =    13     ! INDICATES LABELLING WITH MARS PARAMETERS 
C                           ! FOR WAVE SPECTRA.
C        38     INDICATES CLASS 1=OD 2=RD
C        39     INDICATES TYPE OF FIELD 1=FG 2=AN 9=FC
C        40     INDICATES STREAM  1045=WAVE
C        41     INDICATES EXPERIMENT VERSION NUMBER
C               (4 ASCII characters, righ justified)
C       42,43                                                  0
C
      IF ( CDCLASS .EQ. 'OD' .OR. CDCLASS .EQ. 'od' ) THEN
        ISEC1(38) = 1
      ELSEIF (CDCLASS .EQ. 'ER' .OR. CDCLASS .EQ. 'er' ) THEN
        ISEC1(38) = 3
      ELSEIF (CDCLASS .EQ. 'RD' .OR. CDCLASS .EQ. 'rd' ) THEN
        ISEC1(38) = 2
      ELSE
        WRITE(KU06,*) ' -----------------------------------'
        WRITE(KU06,*) ' SPEC2FDB: ERROR DETERMINING CLASS  '
        WRITE(KU06,*) ' THE CODE FOUND IS : ',  CDCLASS
        WRITE(KU06,*) ' -----------------------------------'
        CALL ABORT1
      ENDIF
      IF ( CDTYPE .EQ. 'fg' .AND. IFCST .EQ. 0 ) THEN
        ISEC1(39) = 1
      ELSEIF ( CDTYPE .EQ. 'an' .AND. IFCST .EQ. 0 ) THEN
        ISEC1(39) = 2
      ELSE
        IF (KTOTENS .EQ. 0) THEN
          ISEC1(39) = 9
        ELSE
          IF(KENSFNB .EQ. 0) THEN
            ISEC1(39) = 10
          ELSE
            ISEC1(39) = 11
          ENDIF
        ENDIF
      ENDIF
C* STREAM.
C  -------
      IF(CDTYPE .EQ. 'cf' .OR. CDTYPE .EQ. 'pf') THEN
        ISEC1(40)=1081
      ELSE
        ISEC1(40)=1045
      ENDIF
C* EXPERIMENT VERSION.
C  -------------------
      ISEC1(41) = 0
      DO JI = 1, 4
        ISEC1(41) = ISEC1(41)*256 + ICHAR(CDEXPVER(JI:JI))
      ENDDO
C* ENSEMBLE LABELS.
C  ----------------
      IF(CDTYPE.EQ.'cf'.OR.CDTYPE.EQ.'pf') THEN
        ISEC1(42) = KENSFNB
        ISEC1(43) = KTOTENS
      ELSE
        ISEC1(42) = 0
        ISEC1(43) = 0
      ENDIF
C* DIRECTION AND FREQUENCY INDEX
C  -----------------------------
C
C     SEE BELOW in LOOP 
C* TOTAL NUMBER OF DIRECTIONS AND FREQUENCIES
C  ------------------------------------------
C
      ISEC1(46) = KANG
      ISEC1(47) = KFRE
C
C* SCALING FACTORS FOR DIRECTION AND FREQUENCY
C  -------------------------------------------
      ISEC1(48) = 1000
      ISEC1(49) = 1000000
C
C* DIRECTIONS
C  ----------
      DO IC=1,KANG
         ISEC1(49+IC) = NINT(ISEC1(48)*PTHETA(IC)*PPDEG)
      ENDDO
C
C* FREQUENCIES
C  -----------
      DO IC=1,KFRE
         ISEC1(49+KANG+IC) = NINT(ISEC1(49)*PFREQ(IC))
      ENDDO
C
C    *SECTION 2*  INTEGER  GRID DEFINITION BLOCK.
C
C       WORD        VALUE   CONTENTS.
C       ----        -----   ---------
C
      ISEC2( 1) =       0  ! DATA REPRESENTATION TYPE: WMO CODE TABLE 6
      ISEC2( 5) =  NINT(PMOWEP*1000) ! LONGITUDE OF ORIGIN (WEST -)
      ISEC2( 6) =     128  ! RESOLUTION FLAG (WMO CODE TABLE 7)
      ISEC2( 8) =  NINT(PMOEAP*1000) ! LONGITUDE OF EXTREME POINT (WEST)
      ISEC2( 9) =  NINT(PDELLO*1000) ! LONGITUDE INCREMENT
      ISEC2(11) =       0  ! SCANNING MODE FLAGS (WMO CODE TABLE 8)
      ISEC2(17) =    KRGG  ! GRID TYPE.
C       12[17]-68                                               56*0
C
C    *ISEC3*  INTEGER  BIT MAP SECTION.
C    ----------------------------------
C
      ISEC3(1) = 0
      ZSEC3(2) = PMISS
C
C    *SECTION 4*  INTEGER  DATA SECTION.
C    -----------------------------------
C
C       WORD    CONTENTS                                  VALUE.
C       ----    --------                                  ------
C         1     NUMBER OF DATA VALUES IN PSEC4             I1*I2
C         2     NUMBER OF BITS FOR EACH PACKED DATA           16
C       3-42                                                40*0
C
C
C ----------------------------------------------------------------------
C
      ISEC4=0
 1000 CONTINUE
      DO 1001 JANG=1,KANG
        ZTHETA(JANG) = PTHETA(JANG)*PPDEG
 1001 CONTINUE
      ZFREQ=PFREQ
      ilenp = kang*kfre
C
      isec4( 1) = ilenp
      isec4( 2) =  9 ! 9 or 16 or 32 
                     ! Number of bits used for each encoded value
      isec4( 3) =  0
      isec4( 4) =  0
      isec4( 5) =  0
      isec4( 6) = 16
      isec4( 7) =  0
C SET HEADER TO SINGLE DATUM AT EACH GRID POINT AND NO SECONDARY BITMAPS
      isec4( 8) =  0
      isec4( 9) =  0
      isec4(10) =  0
      isec4(11) =  0
      isec4(12) =  0
      isec4(13) =  0
      isec4(14) =  0
      isec4(15) =  0
      isec4(16) =  0
      isec4(17) =  0
      isec4(18) =  0
      isec4(19) =  0
      isec4(50) =  1
      isec4(51) =  1
      isec4(52) =  0
      isec4(53) =  KANG 
      isec4(54) =  0
      isec4(55) =  KFRE 
      isec4(56) =  1
      isec4(57) =  2
C
C
      ij=0
      DO j = 60, 60-1+kang
        ij=ij+1
        itheta(ij) = TRANSFER (ztheta(ij), 1)
        isec4(j) = itheta(ij)
      ENDDO
C
      ij=0
      DO j = 60+kang, 60+kang-1+kfre
        ij=ij+1
        ifreq(ij)  = TRANSFER (zfreq(ij), 1)
        isec4(j) = ifreq(ij)
      ENDDO
C-----------------------------------------------------------------------
C
C*    2. CALCULATE OFFSETS.
C     ---------------------
      ALLOCATE(IPOINTER(1:KEND(KPROC)))
C
      IF(KRGG.EQ.1 .OR. CLDOMAIN == 'm' ) THEN
       KDEXN = 0
       KDEXS = 0 
      ELSE
       KDEXN = (NINT((90. - PMONOP ) / PDELLA))*KGX
       KDEXS = (NINT((PMOSOP + 90. ) / PDELLA))*KGX
      ENDIF
C
      IF (KRGG .EQ. 1 ) THEN
        DO IJ = 1, KEND(KPROC)
          IX = KXLG(IJ,1)
          IY = KGY- KXLT(IJ,1) +1
          IPOINTER(IJ)=0
          DO JK = 1,IY-1
            JSN = KGY - JK + 1
            IPOINTER(IJ) = KLONRGG(JSN) + IPOINTER(IJ)
          ENDDO
          IPOINTER(IJ) = IPOINTER(IJ) + IX
        ENDDO
        ITOT=0
        DO JK=1,KGY
          JSN=KGY-JK+1
          ITOT = KLONRGG(JSN) + ITOT 
        ENDDO
      ELSE
        DO IJ = 1, KEND(KPROC)
          IX = KXLG(IJ,1)
          IY = KGY- KXLT(IJ,1) +1
          IPOINTER(IJ) = (IY-1)*KGX + IX + KDEXN
        ENDDO
        ITOT = KGY * KGX + KDEXN + KDEXS
      ENDIF
C
C
C-----------------------------------------------------------------------
C
C*    3. OUTPUT GRIB DATA 
C     -------------------
C
      ISTEP=KPROC
      IF (.NOT.LFDBIOOUT) ISTEP=1
C
      ALLOCATE(ZSEC4OUT(ITOT))
C
      DO IC=1,KFRE*KANG,ISTEP
C
        ALLOCATE(WORK(KIBLO))
        IST=ISTEP
        IF(IC+ISTEP.GT.KFRE*KANG ) IST=KFRE*KANG-IC+1
        ZSEC4OUT = 0.0
C
C       COLLECT CONTRIBUTIONS ON ALL PE's (or a subset if less fields
C       ---------------------------------- than PE's)
        DO KSEND=1,IST
          IM=(((IC-1)+KSEND-1)/KANG)+1
          IK=(IC-1)+KSEND-(IM-1)*KANG
C
          DO J1=KSTART(KRANK),KEND(KRANK)
            WORK(J1) = PSPEC(J1,IK,IM)
          ENDDO
C
          ITAG=(IM-1)*KANG+IK
          CALL MPGATHERSCFLD(KSEND,ITAG,KSTART,KEND,WORK)
          IF (KTEST.GT.1)
     &        WRITE(KU06,*) 'SUB SPEC2FDB: MPGATHERSCFLD DONE' 
C
          IF(ISTEP.EQ.1.AND.KPROC.GT.1) THEN
            CALL MPE_BARRIER(IRET)
            IF(IRET.LT.0)
     &      CALL MPEI_ABORT('MPE_BARRIER  ERROR IN  SPEC2FDB')
          ENDIF
C
          IF(KSEND.EQ.KRANK) THEN
            DO J1=1,KEND(KPROC)
              ZSEC4OUT(IPOINTER(J1)) = WORK(J1)
            ENDDO
          ENDIF
        ENDDO 
        IF(ALLOCATED(WORK)) DEALLOCATE(WORK)
C
C-----------------------------------------------------------------------
C
C       OUTPUT FROM ALL PE's WHICH CONTAIN A CONTRIBUTION
C       -------------------------------------------------
        IF(KRANK.LE.IST) THEN
          IM=(((IC-1)+KRANK-1)/KANG)+1
          IK=(IC-1)+KRANK-(IM-1)*KANG
C
C*        PREPARE FOR ENCODING.
C         ----------------------
C
          DO JC=1,ITOT
            ZSEC4OUT(JC) = LOG10(ZSEC4OUT(JC)+PPEPS)+ABS(PPREC)
            IF ( ZSEC4OUT(JC) .LE. PPMISS ) ZSEC4OUT(JC)=PMISS
          ENDDO
C
          ILENGOUT=KGX*KGY/KIND(KU06)
          ALLOCATE(IOUTGRIB(ILENGOUT))
C
          ISEC1(44) = IK
          ISEC1(45) = IM
          ISEC4( 1) = ITOT
          ISEC2( 2) = KGX
          IF( CLDOMAIN == 'm' ) THEN
            ISEC2( 3) =  KGY
            ISEC2( 4) =  NINT(PMONOP*1000)
            ISEC2( 7) =  NINT(PMOSOP*1000)
          ELSE
            ISEC2( 4) =  90000
            ISEC2( 7) = -90000
            ISEC2( 3) =  NINT(180./PDELLA) + 1
          ENDIF
C
          ISEC2(10) = NINT(PDELLA*1000)
C
          IF (KRGG .EQ. 1) THEN
            IF( CLDOMAIN == 'm' ) THEN
              KST = 0
            ELSE
              KST = NINT((90. - PMONOP ) / PDELLA)
            ENDIF
              DO JC = 1, KGY
                JSN=KGY-JC+1
                ISEC2(22+JC+KST) = KLONRGG(JSN)
              ENDDO
          ENDIF
C
          IF (KTEST.GT.1) THEN
            WRITE(KU06,'(" I KIND(ztheta) :", I10 )') KIND(ztheta)
            WRITE(KU06,'(" I SECTION 0:", /,80(8x, 8I10,/))') ISEC0
            WRITE(KU06,'(" I SECTION 1:", /,80(8x, 8I10,/))') ISEC1
            WRITE(KU06,'(" I SECTION 2:", /,80(8x, 8I10,/))') ISEC2
            WRITE(KU06,'(" I SECTION 3:", /,80(8x, 8I10,/))') ISEC3
            WRITE(KU06,'(" I SECTION 4:", /,80(8x, 8I10,/))') ISEC4
            WRITE(KU06,'(" R SECTION 2:", /,80(8x, 8F10.5,/))') ZSEC2
            WRITE(KU06,'(" R SECTION 3:", /,80(8x, 8F10.5,/))') ZSEC3
            WRITE(KU06,'(" ZTHETA:     ", /,80(8x, 8F10.5,/))') ZTHETA
            WRITE(KU06,'(" ZFREQ:      ", /,80(8x, 8F10.5,/))') ZFREQ
            CALL FLUSH (KU06)
          ENDIF
C
C         CODING
C         ------
          CALL GRIBEX (ISEC0, ISEC1, ISEC2, ZSEC2, ISEC3, ZSEC3, ISEC4,
     &               ZSEC4OUT, ITOT, IOUTGRIB, ILENGOUT ,IWORD,'C',IRET)
          CLMSG='ENCODING ZSEC4OUT'
          IF(IRET.GT.0) GOTO 3500
C
C         OUTPUT
C         ------
          IF (LFDBIOOUT) THEN
            CALL FLD2FDB (ku06, ktest, isec1, ioutgrib, iword, 'write',
     .                    ierr,cdtype, cdclass,
     .                    cdexpver, CDFDB2DSP, cldomain,
     .                    IK, IM, kfdb)
          ELSE
C         output to file should only take place on one PE, thus istep=1
            CMODE='a'
            IF(IM.EQ.1.AND.IK.EQ.1) CMODE='w'
            CALL GRSTNAME(KU06,CDTPRO,CDATEF,1,'SGS',KRANK,KPROC,
     &                    OUTFILEN)
            LFILE=LEN_TRIM(OUTFILEN)
            CALL PBOPEN(IUOUT,OUTFILEN(1:LFILE),CMODE,IRET)
            IF(IRET.LT.0) GOTO 3200
            KOUNT=IWORD*KIND(IWORD)
            CALL PBWRITE(IUOUT,IOUTGRIB,KOUNT,IRET)
            IF(IRET.LT.0) GOTO 3300
            CALL PBCLOSE(IUOUT,IRET)
            IF(IRET.LT.0) GOTO 3400
          ENDIF

        ENDIF
C
C-----------------------------------------------------------------------
C
        IF(ALLOCATED(IOUTGRIB)) DEALLOCATE(IOUTGRIB)
C
      ENDDO
C
C*      5. DEALLOCATION.
C       ----------------
C
      IF(ALLOCATED(ZSEC4OUT)) DEALLOCATE(ZSEC4OUT)
      IF(ALLOCATED(IPOINTER)) DEALLOCATE(IPOINTER)

      RETURN
C
C-----------------------------------------------------------------------
C
C*    6. ERROR MESSAGES.
C     ------------------
C
3200  CONTINUE
      WRITE (KU06,*) '*************************************'
      WRITE (KU06,*) '*                                   *'
      WRITE (KU06,*) '*   ERROR FOLLOWING CALL TO PBOPEN  *'
      WRITE (KU06,*) '*   IN SPEC2FDB                     *'
      IF(IRET.EQ.-1)
     &        WRITE (KU06,*) ' COULD NOT OPEN FILE ',OUTFILEN
      IF(IRET.EQ.-2)
     &        WRITE (KU06,*) ' INVALID FILENAME ',OUTFILEN
      IF(IRET.EQ.-3)
     &        WRITE (KU06,*)' INVALID OPEN MODE SPECIFIED'
      WRITE (KU06,*) '*                                   *'
      WRITE (KU06,*) '*************************************'
      CALL ABORT1
3300  CONTINUE
      WRITE (KU06,*) '**************************************'
      WRITE (KU06,*) '*                                    *'
      WRITE (KU06,*) '*  ERROR FOLLOWING CALL TO PBWRITE   *'
      WRITE (KU06,*) '*  IN SPEC2FDB                       *'
      WRITE (KU06,*) '*  FILE ',OUTFILEN
      WRITE (KU06,*) '*                                    *'
      WRITE (KU06,*) '**************************************'
      CALL ABORT1
3400  CONTINUE
      WRITE (KU06,*) '************************************'
      WRITE (KU06,*) '*                                  *'
      WRITE (KU06,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
      WRITE (KU06,*) '* IN SPEC2FDB                      *'
      WRITE (KU06,*) '* FILE ',OUTFILEN
      WRITE (KU06,*) '*                                  *'
      WRITE (KU06,*) '************************************'
      CALL ABORT1
3500  CONTINUE
      WRITE(KU06,*) '*********************************************'
      WRITE(KU06,*) '*                                           *'
      WRITE(KU06,*) '*         FATAL ERROR IN SUB. SPEC2FDB      *'
      WRITE(KU06,*) '*         ============================      *'
      WRITE(KU06,*) '* GRIBEX ERROR - IRET = ',IRET
      WRITE(KU06,*) '* WHILE ',CLMSG
      WRITE(KU06,*) '*                                           *'
      WRITE(KU06,*) '*    PROGRAM ABORTS     PROGRAM ABORTS      *'
      WRITE(KU06,*) '*                                           *'
      WRITE(KU06,*) '*********************************************'
      CALL ABORT1
C
      END
