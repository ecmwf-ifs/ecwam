      SUBROUTINE gribpacs
     .          ( kuso, kugrs, psec4, ktgrib, kgx, kgy, kang, kfre,
     .            krgg, ptheta, pfreq, kparam, ky, km, kd, kh, kfcst,
     .            pmowep, pmosop, pmoeap, pmonop, klonrgg, kst, knd,
     .            pdella, pdello, pmiss, cdclass, cdexpver, cddomain,
     .            ldfdb, ksect1_3, krank, jpidl, kfdb)
C
C ======================================================================
C
C****  *GRIBPACS*  PACKS WAM MODEL 2D SPECTRA INTO GRIB CODE . 
C                  if it's of the old type, i.e param=250
C
C       L. ZAMBRESKY (GKSS/ECMWF)  JUNE 1987
C       H. GUNTHER   (GKSS/ECMWF)  JUNE 1991   GRIB EDITION 1.
C       B. HANSEN    (MPIfM/ECMWF)  OCTOBER 1992
C                   USED SUBROUTINE GRIBPAC AS A BASE FOR GRIBPACS.
C                   WHICH IS AN INTERFACE TO GRIBEX IN ORDER TO
C                   CODE 2D SPECTRA.
C       J BIDLOT    jpidl will be passed as an input parametre 
C                   jpmap1 is computed locally (common block removed)
C
C       PURPOSE.
C       --------
C         SUBROUTINE PACKS WAVE 2D SPECTRA INTO THE WMO GRIB CODE
C         FORMAT FOR ARCHIVING INTO MARS.
C
C**    INTERFACE.
C      ----------
C        *CALL* *GRIBPAC
C               ( KUSO, KUGRS, PSEC4, KTGRIB, KGX, KGY, KANG, KFRE, 
C                 KRGG, PTHETA, PFREQ, KPARAM, KY, KM, KD, KH, KFCST
C                 PMOWEP, PMOSOP, PMOEAP, PMONOP, KLONRGG, KST, KND,
C                 PDELLA, PDELLO, PMISS, CDCLASS, CDEXPVER, LDFDB,
C                 KSECT1_3 )*
C
C          *KUSO*      LOGICAL UNIT FOR STANDARD OUTPUT.
C          *KUGRS*        LOGICAL UNIT FOR PACKED FIELD DATA.
C          *PSEC4*     FIELD TO BE PACKED.
C          *KTGRIB*    ACTUAL SIZE OF PSEC4 FOR THE FIRST TWO DIM.
C          *KGX*       NUMBER OF POINTS ALONG A LATITUDE (1st DIM).
C          *KGY*       NUMBER OF POINTS ALONG A MERIDIAN (2nd DIM).
C          *KANG*      NUMBER OF DIRECTIONS IN THE MATRIX (3rd DIM).
C          *KFRE*      NUMBER OF FREQUENCIES IN THE MATRIX (4th DIM).
C          *KRGG*      GRID REPESENTATION (0=REGULAR, 1=REDUCED).
C          *PTHETA*    DIRECTION VECTOR.
C          *PFREQ*     FREQUENCY VECTOR.
C          *KPARAM*    PARAMETER IDENTIFIER.
C          *KY*        YEAR OF LAST ANALYSIS (format yyyy).
C          *KM*        MONTH OF LAST ANALYSIS.
C          *KD*        DAY OF LAST ANALYSIS.
C          *KH*        HOUR OF LAST ANALYSIS.
C          *KFCST*     FORECAST TIME FROM BASE DATE-TIME-GROUP.
C          *PMOWEP*    MOST WESTERN LONGITUDE OF THE FIELD.
C          *PMOSOP*    MOST SOUTHERN LATITUDE OF THE FIELD.
C          *PMOEAP*    MOST EASTERN LONGITUDE OF THE FIELD.
C          *PMONOP*    MOST NORTHERN LATITUDE OF THE FIELD.
C          *KLONRGG*   VECTOR OF NUMBER OF POINTS FOR A GIVEN LATITUDE.
C          *KST*       FIRST Y-INDEX FOR PARTIAL CODING.
C          *KND*       LAST Y-INDEX FOR PARTIAL CODING.
C          *PDELLA*    DISTANCE BETWEEN LATITUDES.
C          *PDELLO*    DISTANCE BETWEEN LONGITUDES.
C          *PMISS*     MISSING VALUE INDICATOR.
C          *CDCLASS*   CLASS TO BE USED.
C          *CDEXPVER*  EXPERIMENT VERSION.
C          *LDFDB*     IF TRUE WRITE GRIB TO FDB. OTHERWISE TO FILE.
C          *KSECT1_3*  MODEL VERSION NUMBER.
C          *KFDB*      FDB REFERENCE.
C          *JPIDL*     NUMBER OF BYTES USED TO REPRESENT REAL NUMBERS
C
C      METHOD.
C      -------
C        GRIBEX IS CALLED TO PACK THE FIELDS INTO GRIB CODE.
C
C      EXTERNALS.
C      ----------
C        GRIBEX
C        FLD2FDB
C
C      REFERENCES.
C      -----------
C        SEE MARS USER GUIDE FOR DOCUMENTATION OF MARS FIELDS.
C
C ----------------------------------------------------------------------
C
#include <outctrl.h>
      DATA csubna /"GRIBPAC"/
C
C ----------------------------------------------------------------------
C
#include <parcons.h>
C
C ----------------------------------------------------------------------
C
C ----------------------------------------------------------------------
C
      INTEGER   isec0(2), isec1(43), isec2(384), isec3(2)
      INTEGER   isec4(59+kang+kfre),
     1          igrib(kgx*kgy*kang*kfre/jpidl)
C
      DIMENSION zsec2(96), zsec3(2), psec4(ktgrib*kang*kfre)
C
      CHARACTER cloper*1
C
      DATA cloper / 'C' /
C
C     DATA isec1 /128, 98, 1, 255, 192, 0, 102, 7*0, 1, 5*0, 20, 4*0/
C     DATA isec2 /0, 120, 61, 90000, 0, 128, -90000, 357000,
C    *            3000, 3000, 58*0 /
C     DATA zsec2  /  56*0  /
C     DATA isec3  / 0,  0  /
C     DATA zsec3  / 0., -999. /
C
C    *ISEC1*  INTEGER  PRODUCT DEFINITION BLOCK.
C
C       WORD    CONTENTS                                  VALUE.
C       ----    --------                                  ------
C         1     VERSION NUMBER OF CODE TABLE 2.              128
C         2     ORIGINATING CENTRE IDENTIFIER.                98
C         3     MODEL IDENTIFICATION.                          1
C         4     GRID DEFINITION.                             255
C         5     FLAG (WMO CODE TABLE 1)                      128
C         6     PARAMETER IDENTIFIER (WMO CODE TABLE 2)   KPARAM
C        7-8    TYPE OF LEVEL (WMO CODE TABLE 3)             2*0
C         9     VALUE(S) OF LEVEL (WMO CODE TABLE 3)           0
C        10     YEAR OF DATA                                  KY
C        11     MONTH OF DATA                                 KM
C        12     DAY OF DATA                                   KD
C        13     HOUR OF DATA                                  KH
C        14     MINUTE OF DATA                                 0
C        15     TIME UNIT (WMO CODE TABLE 4)                   1
C        16     TIME RANGE ONE (FCST TIME FROM BASE DTG)   KFCST
C        17     TIME RANGE TWO                                 0
C        18     TIME RANGE FLAG (WMO CODE TABLE 5)             0
C        19     NUMBER AVERAGED.                               0
C        20     NUMBER MISSING FROM AVERAGE.                   0
C        21     CENTURY OF REFERENCE TIME                     20
C        22                                                    0
C        23                                                    0
C        24     0=NO LOCAL USE OF SECTION 1; 1=LOCAL USE       1
C        25-36                                                 0
C        37     1=LABELLING WITH MARS PARAMETERS; 0=OTHERWISE  1
C        38     1=CLASS IS OP 2=CLASS IS RD                 1 or 2
C        39     1=TYPE IS FIRST GUESS; 2=TYPE IS ANALYSIS   1,2 or 9
C               9=TYPE IS FORCAST
C        40     1045=STREAM IS WAVE                         1045
C        41     4 ASCII character righ justified version or experiment
C        42,43                                                 0
C
C
C
C    *ISEC2*  INTEGER  GRID DEFINITION BLOCK.
C
C       WORD    CONTENTS                                  VALUE.
C       ----    --------                                  ------
C         1     DATA REPRESENTATION TYPE (WMO CODE TABLE 6)    0
C         2     NUMBER OF POINTS ALONG A LATITUDE            120
C         3     NUMBER OF POINTS ALONG A MERIDIAN             46
C         4     LATITUDE OF ORIGIN (SOUTH -)               72000
C         5     LONGITUDE OF ORIGIN (WEST -)                   0
C         6     RESOLUTION FLAG (WMO CODE TABLE 7)           128
C         7     LATITUDE OF EXTREME POINT (SOUTH -)       -63000
C         8     LONGITUDE OF EXTREME POINT (WEST -)       357000
C         9     LATITUDE INCREMENT                          3000
C        10     LONGITUDE INCREMENT                         3000
C        11     SCANNING MODE FLAGS (WMO CODE TABLE 8)         0
C       12-68                                               57*0
C
C    *ZSEC2*  REAL     GRID DEFINITON BLOCK   NOT USED.
C
C    *ISEC3*  INTEGER  BIT MAP SECTION        NOT USED.
C
C    *ZSEC3*  REAL     BIT MAP SECTION        NOT USED.
C
C    *ISEC4*  INTEGER  DATA SECTION.
C
C       WORD    CONTENTS                                  VALUE.
C       ----    --------                                  ------
C         1     NUMBER OF DATA VALUES IN PSEC4       KTGRIB*kang*kfre
C         2     NUMBER OF BITS FOR EACH PACKED DATA            9
C       3-42                                                40*0
C
C    *PSEC4*  REAL      UNPACKED DATA FIELD.
C    *ILENP*  INTEGER   LENGTH OF THIS PSEC4.
C    *IGRIB*  INTEGER   ARRAY TO RECEIVE DATA IN GRIB CODE.
C    *IWORD*  INTEGER   NUMBER OF WORDS OCCUPIED BY CODED DATA IN IGRIB.
C    *IRET*   INTEGER   ERROR RETURN FLAG.
C    *CLOPER* CHARACTER REQUESTED FUNCTION
C
C ----------------------------------------------------------------------
C
      INTEGER   klonrgg(kgy)
C
      DIMENSION ptheta(kang), pfreq(kfre)
C
      REAL :: ONE 
      REAL , DIMENSION(KANG) :: ZTHETA
      REAL , DIMENSION(KFRE) :: ZFREQ
C
      DIMENSION itheta(kang), ifreq(kfre)
C
      CHARACTER*1 cddomain, clfdbsf
      CHARACTER*2 cltype
      CHARACTER*2 cdclass
      CHARACTER*4 cdexpver
      CHARACTER*8 clexpverhlp
      LOGICAL ldfdb
C
      DIMENSION ipr(2)
      DATA ipr / 2*1 /
C
C     FOR TIMING:
      REAL wam_user_clock,   ! the clock function.
     .     time0,        ! the reference time.
     .     time_g,       ! time spend in gribex.
     .     time_w        ! time used for writing.
C
C ----------------------------------------------------------------------
C
C*    1. FIX PARAMETERS AND PACK DATA.
C        -----------------------------
C
      ONE=1.0
      msub = msub + 1
      CALL  uraopc (csubna, msub, llsuvi, lltest)
      IF (llsuvi .OR. lltest) THEN
        WRITE(kuso,*)' ',csubna,': '
        itest = 2
      ELSE
        itest = 0
      ENDIF
C
      ibts = n_precision(ipr)  ! ECMWF UTILITY RETURN NUMBER OF BYTES
C                              ! PER INTEGER IN THIS CASE.
C
      IF (lltest ) THEN
        WRITE(kuso,*) ' SUBROUTINE: GRIBPACKS'
        WRITE(kuso,*) ' KANG  =', kang
        WRITE(kuso,*) ' KFRE  =', kfre
        WRITE(kuso,*) ' KGX   =', kgx
        WRITE(kuso,*) ' KGY   =', kgy
        WRITE(kuso,*) ' KTGRIB=', ktgrib
      ENDIF
C
C     Parameters for section 1.
C
C     DATA isec1 /128, 98, 1, 255, 192, 0, 102, 7*0, 1, 5*0, 20, 4*0/
C
      isec1( 1) = 140   !128
      isec1( 2) = 98
      isec1( 3) = ksect1_3
      isec1( 4) = 255
      isec1( 5) = 192
      isec1( 6) = kparam
      isec1( 7) = 102
      isec1( 8) = 0
      isec1( 9) = 0
      isec1(10) = MOD(ky,100)
      if(isec1(10).eq.0) isec1(10) = 100
      isec1(11) = km
      isec1(12) = kd
      isec1(13) = kh
      isec1(14) = 0
      isec1(15) = 1
      isec1(16) = kfcst
      isec1(17) = 0
      isec1(18) = 0
      isec1(19) = 0
      isec1(20) = 0
      if(isec1(10).eq.100) then
        isec1(21) = (ky/100)
      else
        isec1(21) = (ky/100)+1
      endif
      isec1(22) = 0
      isec1(23) = 0
      isec1(24) = 1
      DO i=25,36
        isec1(i)= 0
      ENDDO
      isec1(37) = 1
      IF ( cdclass .EQ. 'OD' ) THEN
        isec1(38) = 1
      ELSEIF ( cdclass .EQ. 'RD' ) THEN
        isec1(38) = 2
      ELSEIF ( cdclass .EQ. 'ER' ) THEN
        isec1(38) = 3
      ELSEIF ( cdclass .EQ. 'E4' ) THEN
        isec1(38) = 5
      ELSE
        WRITE(kuso,*)' CLASS ', cdclass ,' NOT SUPPORTED '
        WRITE(kuso,*)' PROGRAM ABORTS '
        CALL abort1
      ENDIF
      IF ( kfcst .EQ. 0 ) THEN
        isec1(39) = 2
        cltype='an'
      ELSE
        isec1(39) = 9
        cltype='fc'
      ENDIF
      isec1(40) = 1045
C
C     clexpverhlp(5:8)=cdexpver
C     READ(clexpverhlp,'(A8)') isec1(41)
C
      isec1(41) = 0
      DO ji = 1, 4
        isec1(41) = isec1(41)*256 + ICHAR(cdexpver(ji:ji))
      ENDDO
C
      isec1(42) = 0
      isec1(43) = 1
      IF (llsuvi) THEN
        WRITE(kuso,*)' SUB: GRIBPACS - EXPERIMENT VERSION = ', cdexpver
        WRITE(kuso,*)' SUB: GRIBPACS - ISEC1(41) = ', isec1(41)
      ENDIF
C
C
C     Set parameters for section 2.
C
      isec2     =  0
      isec2( 1) =  0
      isec2( 2) =  kgx
      isec2( 3) =  knd - kst +1 ! kgy
      isec2( 4) =  NINT(pmonop*1000)    !   90000
      isec2( 5) =  NINT(pmowep*1000)    !       0
      isec2( 6) =  128
      isec2( 7) =  NINT(pmosop*1000)    !  -90000
      isec2( 8) =  NINT(pmoeap*1000)    !  357000
      isec2( 9) =  NINT(pdella*1000)    !    3000
      isec2(10) =  NINT(pdello*1000)    !    3000
      isec2(11) =  0
      isec2(12) =  0
      isec2(13) =  0
      isec2(14) =  0
      isec2(15) =  0
      isec2(16) =  0
      isec2(17) =  0
      isec2(18) =  0
      isec2(19) =  0
      isec2(20) =  0
      isec2(21) =  0
      isec2(22) =  0
C
      IF (krgg .EQ. 1) THEN
        isec2(1) = 0
        isec2(6) = 128
        isec2(9) = 0
        isec2(17) = 1
C
        ipadding = 0  ! IPADDING IS JUST FOR CONVIENIENCE TO PRETEND
C                     ! THAT WE DEAL WITH FULL GLOBAL FIELDS.
        istart = 0    ! ISTART IS AN OFFSET USED WHEN THERE ARE 
C                     ! A NUMBER (i.e. istart) OF LATITUDES IN THE
C                     ! BEGINNING OF THE FIELD NOT USED.
                      ! istart = NINT((90. - pmonop ) / pdella)
CC  SECTION COMMENTED OUT : IF CHANGED NEED TO CHANGE INMARSS AS WELL !!!
C       IF ( isec2(7) .GE. 0 ) THEN
C         istart = NINT((90. - pmonop ) / pdella)
C         isec2( 4) =   90000
C       ELSEIF ( isec2(4) .LT. 0 ) THEN
C         ipadding = NINT((-90. - pmosop ) / (pdella * (-1.))) -1
C         isec2( 7) =  -90000
C       ENDIF
CC
C
        isec2( 3) =  knd - kst +1 +istart +ipadding
C
        IF (lltest) WRITE(kuso,*)'  jy, isn, isec2(22-kst+1+jy+istart),'
     .                 , '  22-kst+1+jy+istart '
        DO jy=kst,knd  !kgy
           isn=kgy-jy+1
           isec2(22-kst+1+jy+istart) = klonrgg(isn)
           IF (lltest) WRITE(kuso,*) jy, isn, isec2(22-kst+1+jy+istart),
     .                 22-kst+1+jy+istart
        ENDDO
      ENDIF
C
C
C     parameters for section 3.
C
      isec3(1) = 0
CJH   isec3(2) = 0
CJH   zsec3(1) = 0.0
      zsec3(2) = pmiss
C
C     parameters for section 4.
C
      isec4=0
 1000 CONTINUE
      DO 1001 jang=1,kang
        ztheta(jang) = ptheta(jang)*deg
 1001 CONTINUE
      zfreq=pfreq
      ilenp = ktgrib*kang*kfre
      isec4( 1) = ilenp
      isec4( 2) =  9 ! 9 or 16 or 32 !Number of bits used for each encoded value
      isec4( 3) =  0
      isec4( 4) =  0
      isec4( 5) =  0
      isec4( 6) = 16
      isec4( 7) =  0
      isec4( 8) = 64
      isec4( 9) = 32
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
      isec4(50) = kang
      isec4(51) = kfre
      isec4(52) =  0
      isec4(53) = kang
      isec4(54) =  0
      isec4(55) = kfre
      isec4(56) =  1
      isec4(57) =  2
      iret  = 1
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
        ifreq(ij)  = TRANSFER (zfreq(ij),  1)
        isec4(j) = ifreq(ij)
      ENDDO
C
      IF (lltest) THEN
        WRITE(kuso,*) ' ISEC4'
        DO 1332 j1=1,59+kang+kfre
          WRITE(kuso,'(i4,1X,I19)')j1,isec4(J1)
 1332   CONTINUE
C
        WRITE(kuso,*) ' ANGLES '
        WRITE(kuso,'(8F8.0)') ztheta
        WRITE(kuso,*) ' FREQUENCIES '
        WRITE(kuso,'(8F8.5)') zfreq
C
        i1=knd-kst/2
        i2=i1+1
        DO 1333 j1=i1, i2
          WRITE(kuso,*) ' SPECTRUM ',j1
          js = ((j1-1)*kang*kfre)+1
          je =  (j1-1)*kang*kfre +kang*kfre
          WRITE(kuso,'(12f8.4)')(psec4(j2),j2=js,je)
 1333   CONTINUE
      ENDIF
C
      kkk = 1
      IF (lltest) THEN
        CALL grsdbg ( kkk )
C
C                  CALL GRSVCK (I) where
C                              I = 0, to switch off checking.
C                                  Non-zero to switch on checking.
C
        CALL grsvck (   1 )
      ENDIF
C
        CALL grsdbg ( kkk )
        CALL grsvck (   1 )
C
      IF (lltest) WRITE(kuso,*) (isec2(jsec2),jsec2=1,50)
      time0=-wam_user_clock()
      IF (lltest.OR.llsuvi) CALL FLUSH(kuso)
      CALL gribex (isec0, isec1, isec2, zsec2, isec3, zsec3, isec4,
     .             psec4, ilenp, igrib, kang*kfre*kgx*kgy/jpidl, iword,
     .             cloper, iret)
C
      IF (llsuvi) THEN
        time_g=time0+wam_user_clock()
        time_g=time_g*1E-06
        WRITE (kuso,*) ' ++++++++++++++++++++++++++++++'
        WRITE (kuso,*) ' +   TIME SPEND IN GRIBEX     +'
        WRITE (kuso,*) ' +                            +'
        WRITE (kuso,*) ' + ', time_g
        WRITE (kuso,*) ' +                            +'
        WRITE (kuso,*) ' ++++++++++++++++++++++++++++++'
      ENDIF
C
      IF (lltest) WRITE(kuso,*)' ', kkk
C
C ----------------------------------------------------------------------
C
C*    2. CHECK ERROR CODE.
C        ----------------
C
 2000 CONTINUE
C
      IF (llsuvi)
     &WRITE (kuso,*) '  SUB. GRIBPACS:  GRIBEX RETURN CODE ',iret,
     .              ' FOR PARAMETER ', kparam
C
      IF (iret .EQ. 0) THEN
C
C*    2.1 WRITE OUT PACKED DATA.
C        -----------------------
C
        IF ( ldfdb ) THEN
C
C*       2.1.1 WRITE FIELD TO FDB.
C              -------------------
C
          IF (llsuvi) time0=-wam_user_clock()
          idir = 0
          ifre = 0
          clfdbsf = ' '
          CALL FLD2FDB (kuso, itest, isec1, igrib, iword, 'WRITE', ierr,
     .                  cltype, cdclass, cdexpver, clfdbsf, cddomain,
     .                  idir, ifre, kfdb)
          IF (llsuvi) THEN
            time_w=time0+wam_user_clock()
            time_w=time_w*1E-06
            WRITE (kuso,*) ' ++++++++++++++++++++++++++++++'
            WRITE (kuso,*) ' + TIME USED FOR FLD2FDB      +'
            WRITE (kuso,*) ' +                            +'
            WRITE (kuso,*) ' + ', time_w
            WRITE (kuso,*) ' +                            +'
            WRITE (kuso,*) ' ++++++++++++++++++++++++++++++'
          ENDIF
C
          WRITE (kuso,*) '    SUB. GRIBPACS:  PACKED FIELD LENGTH ',
     .                  iword, '    WRITTEN TO FIELDS DATA BASE ',
     .                  '(REFERENCE=',kfdb,')'
C
C*       2.1.2 PRINT FDB ERROR CODE.
C              ---------------------
C
          IF (ierr .NE. 0) THEN
            WRITE(kuso,*) ' ------------------------'
            WRITE(kuso,*) ' ERROR ACCESSING FDB '
            WRITE(kuso,*) ' FDB ERROR CODE IS ',ierr
            WRITE(kuso,*) ' ------------------------'
            CALL abort
          ENDIF
        ELSE
C
C*       2.1.3 WRITE FIELD TO FILE.
C              --------------------
C
          CALL PBWRITE(kugrs, igrib, iword*ibts, iret)
          IF (iret .LT. 0) THEN
            WRITE(kuso,*) 'PBWRITE: IERR = ', ierr
            WRITE(kuso,*) 'Error in writing to grib file'
            CALL abort
          ENDIF
C           
          IF (llsuvi) THEN
            IF (kugrs .GT. 111) THEN
              WRITE (kuso,2200)'  SUB. GRIBPACS:  PACKED FIELD LENGTH ',
     .                     iword, '    WRITTEN TO UNIT KUGRS = ', kugrs
 2200         FORMAT(1x,a,i6,a,i8)
            ELSE
              WRITE (kuso,*) '  SUB. GRIBPACS:  PACKED FIELD LENGTH ',
     .                     iword, '    WRITTEN TO UNIT KUGRS = ', kugrs
            ENDIF
          ENDIF
 100      CONTINUE
        ENDIF
      ELSE
C
C*    2.3 PRINT GRIBEX ERROR CODE.
C         ------------------------
C
        WRITE(kuso,*) ' -----------------------------------------'
        WRITE(kuso,*) ' ERROR PACKING FIELD BY SUB GRIBEX. '
        WRITE(kuso,*) ' GRIBEX ERROR CODE IS ',iret
        WRITE(kuso,*) ' -----------------------------------------'
        CALL abort
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    3. PRINT GRIB SECTION 1 AND 2.
C         --------------------------
C
 3000 CONTINUE
      IF (llsuvi) CALL GRPRS1 (isec0, isec1)
      IF (llsuvi) CALL GRPRS2 (isec0, isec2, zsec2)
      IF (llsuvi) WRITE(kuso,*) ' PACKED FIELD LENGTH ',iword
      msub = msub - 1
      RETURN
      END
