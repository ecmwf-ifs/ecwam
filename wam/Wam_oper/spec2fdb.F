
      SUBROUTINE SPEC2FDB (KU06, KPROC, KRANK, LDWCOU, PTHETA, PFREQ,
     &                     PSPEC, KINF, KSUP, KANG, KFRE, KGX, KGY,
     &                     KLONRGG, KXLT, KXLG, KIBLO, KBLO,
     &                     PMONOP, PDELLA, KSTART, KEND,
     &                     NKSEK0, NKSEK1, NKSEK2, NKSEK3, NKSEK4,
     &                     NPSEK2, NPSEK3, KSEK0, KSEK1, KSEK2, KSEK3,
     &                     KSEK4, PSEK2, PSEK3, NTENCODE,
     &                     PPDEG, PPMISS, PPEPS, PPREC,
     &                     CDTPRO, CDATEF, CDATEA,
     &                     CDDOMAIN, CDTYPE, CDCLASS, CDEXPVER,
     &                     CDFDB2DSP, KENSFNB, KTOTENS, KSYSNB, KMETNB,
     &                     LFDBIOOUT, KCPLEN, CPAD, KTEST, KFDB )

!----------------------------------------------------------------------

!**** *SPEC2FDB*  ENCODES SPECTRA AS PARAMETER 251
!                 AND WRITES TO FILE OR TO FDB.

!     B. HANSEN   ECMWF   FEBRUARY 1998 
!     J. BIDLOT   ECMWF   NOVEMBER 1998

!*    PURPOSE.
!     --------

!       ENCODES SPECTRA AND WRITES TO FDB OR TO A FILE.

!**   INTERFACE.
!     ----------

!     SUBROUTINE SPEC2FDB (KU06, KPROC, KRANK, LDWCOU, PTHETA, PFREQ,
!    &                     PSPEC, KINF, KSUP, KANG, KFRE, KGX, KGY,
!    &                     KLONRGG, KXLT, KXLG, KIBLO, KBLO,
!    &                     PMONOP, PDELLA, KSTART, KEND, 
!    &                     NKSEK0, NKSEK1, NKSEK2, NKSEK3, NKSEK4,
!    &                     NPSEK2, NPSEK3, KSEK0, KSEK1, KSEK2, KSEK3,
!    &                     KSEK4, PSEK2, PSEK3, NTENCODE,
!    &                     PPDEG, PPMISS, PPEPS, PPREC,
!    &                     CDTPRO, CDATEF, CDATEA,
!    &                     CDDOMAIN, CDTYPE, CDCLASS, CDEXPVER,
!    &                     CDFDB2DSP, KENSFNB, KTOTENS, LFDBIOOUT,
!    &                     KCPLEN, CPAD, KTEST, KFDB )

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *KU06*      INTEGER   OUTPUT UNIT FOR STANDARD OUTPUT.
!      *KPROC*     INTEGER   TOTAL NUMBER OF PES.
!      *KRANK*     INTEGER   CURRENT PE NUMBER.
!      *LDWCOU*    LOGICAL   TRUE: THIS IS A COUPLED RUN.
!      *PTHETA*    REAL      DIRECTION VECTOR.
!      *PFREQ*     REAL      FREQUENCY VECTOR.
!      *PSPEC*     REAL      LOCAL SPECTRA OF CURRENT PE.
!      *KINF*      INTEGER   SMALLEST GRID POINT INDEX USED BY THE PE.
!      *KSUP*      INTEGER   LARGEST GRID POINT INDEX USED BY THE PE.
!      *KANG*      INTEGER   NUMBER OF ANGLES.
!      *KFRE*      INTEGER   NUMBER OF FREQUENCIES.
!      *KGX*       INTEGER   NUMBER OF LONGITUDES IN GRID.
!      *KGY*       INTEGER   NUMBER OF LATITUDES IN GRID.
!      *KLONRGG*   INTEGER   NUMBER OF POINTS ALONG A MERIDIAN.
!      *KXLT*      INTEGER   LAT. GRID INDEX.
!      *KXLG*      INTEGER   LON. GRID INDEX.
!      *KIBLO*     INTEGER   NUMBER OF POINTS IN BLOCK.
!      *KBLO*      INTEGER   NUMBER OF BLOCKS.
!      *PMONOP*    REAL      MOST NORTHERN LATITUDE IN GRID (DEGREE).
!      *PDELLA*    REAL      DISTANCE BETWEEN LATITUDES.
!      *KSTART*    INTEGER   FIRST INDICES OF LOCAL BLOCKS OF SPECTRA.
!      *KEND*      INTEGER   LAST INDICES OF LOCAL BLOCKS OF SPECTRA.
!      *NKSEK0*    INTEGER   SIZE OF KSEK0
!      *NKSEK1*    INTEGER   SIZE OF KSEK1
!      *NKSEK2*    INTEGER   SIZE OF KSEK2
!      *NKSEK3*    INTEGER   SIZE OF KSEK3
!      *NKSEK4*    INTEGER   SIZE OF KSEK4
!      *NPSEK2*    INTEGER   SIZE OF PSEK2
!      *NPSEK3*    INTEGER   SIZE OF PSEK3
!      *KSEK0*     INTEGER   DEFAULT VALUES FOR GRIB HEADER 0
!      *KSEK1*     INTEGER   DEFAULT VALUES FOR GRIB HEADER 1
!      *KSEK2*     INTEGER   DEFAULT VALUES FOR GRIB HEADER 2
!      *KSEK3*     INTEGER   DEFAULT VALUES FOR GRIB HEADER 3
!      *KSEK4*     INTEGER   DEFAULT VALUES FOR GRIB HEADER 4
!      *PSEK2*     REAL      DEFAULT VALUES FOR REAL GRIB HEADER 2
!      *PSEK3*     REAL      DEFAULT VALUES FOR REAL GRIB HEADER 3
!      *NTENCODE*  INTEGER   TOTAL NUMBER OF GRID DATA FOR ENCODING
!      *PPDEG*     REAL      PI/180
!      *PPMISS*    REAL      MINIMUM VALUE FOR CODING OF LOGARITHMIC
!                            SPECTRA.
!      *PPEPS*     REAL      SMALL NUMBER FOR CODING OF SPECTRA
!      *PPREC*     REAL      REFERENCE VALUE FOR CODING OF SPECTRA
!      *CDTPRO*    CHAR*12   END DATE OF PROPAGATION.
!      *CDATEF*    CHAR*12   END DATE OF ANALYSIS RUN (YYMMDDHHMM).
!      *CDATEA*    CHAR*12   BEGIN DATE OF ANALYSIS RUN (YYMMDDHHMM).
!      *CDDOMAIN   CHAR*1    GRIB DOMAIN DEFINITION FOR THE FDB
!      *CDTYPE*    CHAR*     MARS TYPE OF FIELD.
!      *CDCLASS*   CHAR*     MARS CLASS OF FIELD.
!      *CDEXPVER*  CHAR*     EXPERIMENT VERSION IDENTIFIER.
!      *CDFDB2DSP* CHAR*     ROOT DIRECTORY FOR 2D GRIB SPECTRA.
!      *KENSFNB*   INTEGER   ENSEMBLE FORECAST NUMBER
!      *KTOTENS*   INTEGER   TOTAL NUMBER OF ENSEMBLE FORECASTS.
!      *LFDBIOOUT* LOGICAL   IF TRUE OUTPUT TO FDB OTHERWISE TO FILE
!      *KCPLEN*    INTEGER    LENGTH OF CPAD
!      *CPAD*      CHARACTER  PATH FOR OUTPUT TO DISK
!      *KTEST*     INTEGER   PRINT DIAGNOSTICS TO STDOUT IF GT 1.
!      *KFDB*      INTEGER    DATA BASE REFERENCE

!     METHOD.
!     -------

!           ENCODE SPECTRA PER FREQUENCY AND DIRECTION
!      INTO GRIB AND WRITE TO FDB OR TO A SINGLE FILE.

!     EXTERNALS.
!     ----------

!      *ABORT1*
!      *GRIBEX* 
!      *FLD2FDB*
!      *MPGATHERSCFLD*

!     REFERENCE.
!     ----------

!       NONE.

!-------------------------------------------------------------------

      IMPLICIT NONE

!*  DUMMY ARGUMENTS:
!   ----------------

      LOGICAL LDWCOU, LFDBIOOUT

      INTEGER KU06, KINF, KSUP, KANG, KFRE, KGX, KGY, KIBLO, KBLO,
     &        KRANK, KPROC, KSEND, IK, IM, KENSFNB, KTOTENS, KSYSNB,
     &        KMETNB, KLONRGG(KGY), LFILE, IUOUT, KOUNT, KTEST, KFDB,
     &        KXLT(KIBLO, KBLO), KXLG(KIBLO, KBLO), KDEXN,
     &        NKSEK0, NKSEK1, NKSEK2, NKSEK3, NKSEK4, NPSEK2, NPSEK3,
     &        NTENCODE, KCPLEN
      INTEGER, DIMENSION(KPROC) :: KSTART, KEND 

      INTEGER KSEK0(NKSEK0), KSEK1(NKSEK1), KSEK2(NKSEK2),
     &        KSEK3(NKSEK3), KSEK4(NKSEK4)

      REAL    PMONOP, PDELLA,
     &        PSPEC (KINF-1:KSUP, KANG, KFRE),
     &        PPDEG, PPMISS, PPEPS, PPREC,
     &        PTHETA(KANG),
     &        PFREQ(KFRE)
      REAL PSEK2(NPSEK2), PSEK3(NPSEK3)

      CHARACTER*(*) CDTPRO, CDATEF, CDATEA, CPAD,
     &              CDTYPE, CDCLASS, CDEXPVER, CDFDB2DSP

!-----------------------------------------------------------------------

      INTEGER, ALLOCATABLE :: IOUTGRIB(:), IPOINTER(:)

      INTEGER ISEC0(NKSEK0), ISEC1(NKSEK1), ISEC2(NKSEK2),
     &        ISEC3(NKSEK3), ISEC4(NKSEK4)
      INTEGER JSN, JK, JC, IX, ISTEP, IST, IYYYY,
     &        IY, IMO, ID, IH, IMN, IFCST, ITHETA(KANG), IFREQ(KFRE),
     &        IJ, JANG, J1, ILENGOUT, IC,
     &        IRET, IWORD, IERR, ILENP, ITAG

      REAL, ALLOCATABLE :: ZSEC4OUT(:),WORK(:)
      REAL ZSEC2(NPSEK2), ZSEC3(NPSEK3)
      REAL , DIMENSION(KANG) :: ZTHETA
      REAL , DIMENSION(KFRE) :: ZFREQ

      CHARACTER*  1 CDDOMAIN
      CHARACTER* 40 CLMSG
      CHARACTER* 120 OUTFILEN

!-----------------------------------------------------------------------

!*    0. FIXED SETTINGS FOR THE GRIB HEADERS.
!     ---------------------------------------

      ISEC0 = KSEK0

      IF (KTEST.GT.1) THEN
        WRITE(KU06,*) '*      THIS IS SPEC2FDB         *'
        CALL FLUSH (KU06)
      ENDIF

      IF(CDTPRO.LE.CDATEF) THEN

!*    0.1.  THIS IS AN ANALYSIS DATE.
!           -------------------------

        IF (LDWCOU .AND. CDTYPE .EQ. 'fg') THEN
          READ (CDATEA, '(I4,4I2)') IYYYY, IMO, ID, IH, IMN
          CALL DIFDATE (CDATEA, CDTPRO, IFCST)
          IFCST = IFCST/3600
        ELSE
          READ (CDTPRO, '(I4,4I2)') IYYYY, IMO, ID, IH, IMN
          IFCST = 0
        ENDIF
      ELSE

!*    0.2.  THIS IS A  FORECAST DATE.
!           -------------------------

        READ (CDATEF, '(I4,4I2)') IYYYY, IMO, ID, IH, IMN
        CALL DIFDATE (CDATEF, CDTPRO, IFCST)
        IFCST = IFCST/3600
      ENDIF

!    *SECTION 1*  INTEGER  PRODUCT DEFINITION BLOCK.
!    -----------------------------------------------
!       WORD        VALUE   CONTENTS.
!       ----        -----   ---------

      ISEC1 = KSEK1

      ISEC1( 6) =   251     ! PARAMETER IDENTIFIER (WMO CODE TABLE 2)
      ISEC1(10) = MOD(IYYYY,100)      ! YEAR OF DATA
      IF(ISEC1(10).EQ.0) ISEC1(10) = 100
      ISEC1(11) =   IMO     ! MONTH OF DATA
      ISEC1(12) =   ID      ! DAY OF DATA
      ISEC1(13) =   IH      ! HOUR OF DATA
      ISEC1(14) =   IMN     ! MINUTE OF DATA
      ISEC1(16) =  IFCST    ! TIME RANGE ONE (FCST TIME FROM BASE DTG)
      ISEC1(18) =     0     ! TIME RANGE FLAG (WMO CODE TABLE 5)
      IF(ISEC1(10).EQ.100) THEN
        ISEC1(21) = IYYYY/100 ! CENTURY OF REFERENCE TIME.
      ELSE
        ISEC1(21) = (IYYYY/100)+1 ! CENTURY OF REFERENCE TIME.
      ENDIF
      ISEC1(37) =    13     ! INDICATES LABELLING WITH MARS PARAMETERS 
!                           ! FOR WAVE SPECTRA.
!        39     INDICATES TYPE OF FIELD 1=FG 2=AN 9=FC
      IF ( CDTYPE .EQ. 'fg' .AND. IFCST .EQ. 0 ) THEN
        ISEC1(39) = 1
      ELSEIF ( CDTYPE .EQ. 'an' .AND. IFCST .EQ. 0 ) THEN
        ISEC1(39) = 2
      ELSEIF ( CDTYPE .EQ. 'sf' ) THEN
        ISEC1(39) = 9
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
!* DIRECTION AND FREQUENCY INDEX
!  -----------------------------

!     SEE BELOW in LOOP 
!* TOTAL NUMBER OF DIRECTIONS AND FREQUENCIES
!  ------------------------------------------

      ISEC1(46) = KANG
      ISEC1(47) = KFRE

!* SCALING FACTORS FOR DIRECTION AND FREQUENCY
!  -------------------------------------------
      ISEC1(48) = 1000
      ISEC1(49) = 1000000

!* DIRECTIONS
!  ----------
      DO IC=1,KANG
        ISEC1(49+IC) = NINT(ISEC1(48)*PTHETA(IC)*PPDEG)
      ENDDO

!* FREQUENCIES
!  -----------
      DO IC=1,KFRE
        ISEC1(49+KANG+IC) = NINT(ISEC1(49)*PFREQ(IC))
      ENDDO

!* SYSTEM AND METHOD (for seasonal forecast only)
!  -----------------
      IF ( CDTYPE .EQ. 'sf' ) THEN
        ISEC1(50+KANG+KFRE) = KSYSNB 
        ISEC1(51+KANG+KFRE) = KMETNB 
      ELSE
        ISEC1(50+KANG+KFRE) = 65535
        ISEC1(51+KANG+KFRE) = 65535
      ENDIF

!    *SECTION 2* GRID DEFINITION BLOCK. 
!    -----------------------------------

      ISEC2 = KSEK2

      ZSEC2 = PSEK2 

!    *SECTION 3* BIT MAP SECTION.
!    ----------------------------

      ISEC3 = KSEK3 

      ZSEC3 = PSEK3 

!    *SECTION 4*  INTEGER  DATA SECTION.
!    -----------------------------------

      ISEC4=KSEK4

      DO JANG=1,KANG
        ZTHETA(JANG) = PTHETA(JANG)*PPDEG
      ENDDO
      ZFREQ=PFREQ
      ILENP = KANG*KFRE

      ISEC4( 1) = ILENP
      ISEC4( 2) =  9 ! 9 or 16 or 32 
      ! Number of bits used for each encoded value
      ISEC4( 6) = 16 ! additional flags
      ISEC4(50) =  1
      ISEC4(51) =  1
      ISEC4(53) =  KANG 
      ISEC4(55) =  KFRE 
      ISEC4(56) =  1
      ISEC4(57) =  2

      DO IJ = 1, KANG
        ITHETA(IJ) = TRANSFER (ZTHETA(IJ), 1)
        ISEC4(60+IJ-1) = ITHETA(IJ)
      ENDDO

      DO IJ = 1, KFRE
        IFREQ(IJ)  = TRANSFER (ZFREQ(IJ), 1)
        ISEC4(60+KANG+IJ-1) = IFREQ(IJ)
      ENDDO
!-----------------------------------------------------------------------

!*    2. CALCULATE OFFSETS.
!     ---------------------
      ALLOCATE(IPOINTER(1:KEND(KPROC)))

      IF(KSEK2(17).EQ.1 .OR. CDDOMAIN == 'm' ) THEN
        KDEXN = 0
      ELSEIF(CDDOMAIN == 's' ) THEN
        KDEXN = 0
      ELSE
        KDEXN = (NINT((90. - PMONOP ) / PDELLA))*KGX
      ENDIF

      IF (KSEK2(17) .EQ. 1 ) THEN
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
      ELSE
        DO IJ = 1, KEND(KPROC)
          IX = KXLG(IJ,1)
          IY = KGY- KXLT(IJ,1) +1
          IPOINTER(IJ) = (IY-1)*KGX + IX + KDEXN
        ENDDO
      ENDIF


!-----------------------------------------------------------------------

!*    3. OUTPUT GRIB DATA 
!     -------------------

      ISTEP=KPROC
      IF (.NOT.LFDBIOOUT) ISTEP=1

      ALLOCATE(ZSEC4OUT(NTENCODE))

      DO IC=1,KFRE*KANG,ISTEP

        IST=ISTEP
        IF(IC+ISTEP.GT.KFRE*KANG ) IST=KFRE*KANG-IC+1
        ZSEC4OUT = 0.0

!       COLLECT CONTRIBUTIONS ON ALL PE's (or a subset if less fields
!       ---------------------------------- than PE's)
        IF(KPROC.GT.1) THEN
          ALLOCATE(WORK(KIBLO))

          DO KSEND=1,IST
            IM=(((IC-1)+KSEND-1)/KANG)+1
            IK=(IC-1)+KSEND-(IM-1)*KANG

            DO J1=KSTART(KRANK),KEND(KRANK)
              WORK(J1) = PSPEC(J1,IK,IM)
            ENDDO

            ITAG=(IM-1)*KANG+IK
            CALL MPGATHERSCFLD(KSEND,ITAG,KSTART,KEND,WORK)
            IF (KTEST.GT.1)
     &          WRITE(KU06,*) 'SUB SPEC2FDB: MPGATHERSCFLD DONE' 

            IF(ISTEP.EQ.1.AND.KPROC.GT.1) THEN
              CALL MPE_BARRIER(IRET)
              IF(IRET.LT.0)
     &        CALL MPEI_ABORT('MPE_BARRIER  ERROR IN  SPEC2FDB')
            ENDIF

            IF(KSEND.EQ.KRANK) THEN
              DO J1=1,KEND(KPROC)
                ZSEC4OUT(IPOINTER(J1)) = WORK(J1)
              ENDDO
            ENDIF
          ENDDO 
          IF(ALLOCATED(WORK)) DEALLOCATE(WORK)
        ELSE
          IM=((IC-1)/KANG)+1
          IK=IC-(IM-1)*KANG
 
          DO J1=1,KEND(KPROC)
            ZSEC4OUT(IPOINTER(J1)) = PSPEC(J1,IK,IM)
          ENDDO
        ENDIF

!-----------------------------------------------------------------------

!       OUTPUT FROM ALL PE's WHICH CONTAIN A CONTRIBUTION
!       -------------------------------------------------
        IF(KRANK.LE.IST) THEN
          IM=(((IC-1)+KRANK-1)/KANG)+1
          IK=(IC-1)+KRANK-(IM-1)*KANG

!*        PREPARE FOR ENCODING.
!         ----------------------

          DO JC=1,NTENCODE
            ZSEC4OUT(JC) = LOG10(ZSEC4OUT(JC)+PPEPS)+ABS(PPREC)
            IF ( ZSEC4OUT(JC) .LE. PPMISS ) ZSEC4OUT(JC)=ZSEC3(2)
          ENDDO

          ILENGOUT = NINT(NTENCODE*(ISEC4(2)+1.)/(KIND(KU06)*8))+
     &      NKSEK0+NKSEK1+NKSEK2+NKSEK3+NKSEK4+NPSEK2+NPSEK3

          ALLOCATE(IOUTGRIB(ILENGOUT))

          ISEC1(44) = IK
          ISEC1(45) = IM
          ISEC4( 1) = NTENCODE 

          IF (KTEST.GT.3) THEN
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

!         CODING
!         ------
          CALL GRIBEX (ISEC0, ISEC1, ISEC2, ZSEC2, ISEC3, ZSEC3, ISEC4,
     &          ZSEC4OUT, NTENCODE, IOUTGRIB, ILENGOUT ,IWORD,'C',IRET)
          CLMSG='ENCODING ZSEC4OUT'
          IF(IRET.GT.0) GOTO 3500

!         OUTPUT
!         ------
          IF (LFDBIOOUT) THEN
            CALL FLD2FDB (ku06, ktest, isec1, ioutgrib, iword, 'write',
     &                    ierr,cdtype, cdclass,
     &                    cdexpver, CDFDB2DSP, cddomain,
     &                    IK, IM, kfdb)
          ELSE
!         output to file should only take place on one PE, thus istep=1
            IF(IM.EQ.1.AND.IK.EQ.1) THEN
              CALL GRSTNAME(CDTPRO,CDATEF,'SGS',KCPLEN,CPAD,OUTFILEN)
              LFILE=LEN_TRIM(OUTFILEN)
              CALL PBOPEN(IUOUT,OUTFILEN(1:LFILE),'w',IRET)
              IF(IRET.LT.0) GOTO 3200
            ENDIF

            KOUNT=IWORD*KIND(IWORD)
            CALL PBWRITE(IUOUT,IOUTGRIB,KOUNT,IRET)
            IF(IRET.LT.0) GOTO 3300

            IF(IM.EQ.KFRE.AND.IK.EQ.KANG) THEN
              CALL PBCLOSE(IUOUT,IRET)
              IF(IRET.LT.0) GOTO 3400
            ENDIF
          ENDIF

        ENDIF

!-----------------------------------------------------------------------

        IF(ALLOCATED(IOUTGRIB)) DEALLOCATE(IOUTGRIB)

      ENDDO

!*      5. DEALLOCATION.
!       ----------------

      IF(ALLOCATED(ZSEC4OUT)) DEALLOCATE(ZSEC4OUT)
      IF(ALLOCATED(IPOINTER)) DEALLOCATE(IPOINTER)

      RETURN

!-----------------------------------------------------------------------

!*    6. ERROR MESSAGES.
!     ------------------

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

      END SUBROUTINE SPEC2FDB




