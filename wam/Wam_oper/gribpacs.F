      SUBROUTINE gribpacs
     &          ( kuso, kugrs, psec4, ktgrib, kgx, kgy, kang, kfre,
     &            krgg, ptheta, pfreq, kparam, ky, km, kd, kh, kfcst,
     &            pmowep, pmosop, pmoeap, pmonop, klonrgg, kst, knd,
     &            pdella, pdello, pmiss, cdclass, cdexpver, cddomain,
     &            ldfdb, ksect1_3, krank, jpidl, kfdb)

! ======================================================================

!****  *GRIBPACS*  PACKS WAM MODEL 2D SPECTRA INTO GRIB CODE . 
!                  if it's of the old type, i.e param=250

!       L. ZAMBRESKY (GKSS/ECMWF)  JUNE 1987
!       H. GUNTHER   (GKSS/ECMWF)  JUNE 1991   GRIB EDITION 1.
!       B. HANSEN    (MPIfM/ECMWF)  OCTOBER 1992
!                   USED SUBROUTINE GRIBPAC AS A BASE FOR GRIBPACS.
!                   WHICH IS AN INTERFACE TO GRIBEX IN ORDER TO
!                   CODE 2D SPECTRA.
!       J BIDLOT    jpidl will be passed as an input parametre 
!                   jpmap1 is computed locally (common block removed)

!       PURPOSE.
!       --------
!         SUBROUTINE PACKS WAVE 2D SPECTRA INTO THE WMO GRIB CODE
!         FORMAT FOR ARCHIVING INTO MARS.

!**    INTERFACE.
!      ----------
!        *CALL* *GRIBPAC
!               ( KUSO, KUGRS, PSEC4, KTGRIB, KGX, KGY, KANG, KFRE, 
!                 KRGG, PTHETA, PFREQ, KPARAM, KY, KM, KD, KH, KFCST
!                 PMOWEP, PMOSOP, PMOEAP, PMONOP, KLONRGG, KST, KND,
!                 PDELLA, PDELLO, PMISS, CDCLASS, CDEXPVER, LDFDB,
!                 KSECT1_3 )*

!          *KUSO*      LOGICAL UNIT FOR STANDARD OUTPUT.
!          *KUGRS*        LOGICAL UNIT FOR PACKED FIELD DATA.
!          *PSEC4*     FIELD TO BE PACKED.
!          *KTGRIB*    ACTUAL SIZE OF PSEC4 FOR THE FIRST TWO DIM.
!          *KGX*       NUMBER OF POINTS ALONG A LATITUDE (1st DIM).
!          *KGY*       NUMBER OF POINTS ALONG A MERIDIAN (2nd DIM).
!          *KANG*      NUMBER OF DIRECTIONS IN THE MATRIX (3rd DIM).
!          *KFRE*      NUMBER OF FREQUENCIES IN THE MATRIX (4th DIM).
!          *KRGG*      GRID REPESENTATION (0=REGULAR, 1=REDUCED).
!          *PTHETA*    DIRECTION VECTOR.
!          *PFREQ*     FREQUENCY VECTOR.
!          *KPARAM*    PARAMETER IDENTIFIER.
!          *KY*        YEAR OF LAST ANALYSIS (format yyyy).
!          *KM*        MONTH OF LAST ANALYSIS.
!          *KD*        DAY OF LAST ANALYSIS.
!          *KH*        HOUR OF LAST ANALYSIS.
!          *KFCST*     FORECAST TIME FROM BASE DATE-TIME-GROUP.
!          *PMOWEP*    MOST WESTERN LONGITUDE OF THE FIELD.
!          *PMOSOP*    MOST SOUTHERN LATITUDE OF THE FIELD.
!          *PMOEAP*    MOST EASTERN LONGITUDE OF THE FIELD.
!          *PMONOP*    MOST NORTHERN LATITUDE OF THE FIELD.
!          *KLONRGG*   VECTOR OF NUMBER OF POINTS FOR A GIVEN LATITUDE.
!          *KST*       FIRST Y-INDEX FOR PARTIAL CODING.
!          *KND*       LAST Y-INDEX FOR PARTIAL CODING.
!          *PDELLA*    DISTANCE BETWEEN LATITUDES.
!          *PDELLO*    DISTANCE BETWEEN LONGITUDES.
!          *PMISS*     MISSING VALUE INDICATOR.
!          *CDCLASS*   CLASS TO BE USED.
!          *CDEXPVER*  EXPERIMENT VERSION.
!          *LDFDB*     IF TRUE WRITE GRIB TO FDB. OTHERWISE TO FILE.
!          *KSECT1_3*  MODEL VERSION NUMBER.
!          *KFDB*      FDB REFERENCE.
!          *JPIDL*     NUMBER OF BYTES USED TO REPRESENT REAL NUMBERS

!      METHOD.
!      -------
!        GRIBEX IS CALLED TO PACK THE FIELDS INTO GRIB CODE.

!      EXTERNALS.
!      ----------
!        GRIBEX
!        FLD2FDB

!      REFERENCES.
!      -----------
!        SEE MARS USER GUIDE FOR DOCUMENTATION OF MARS FIELDS.

! ----------------------------------------------------------------------

      USE YOWPCONS , ONLY : DEG 
      USE YOMGSTATS, ONLY : LSYNCSTATS

#include <outctrl.h>
      DATA csubna /"GRIBPAC"/

! ----------------------------------------------------------------------

! ----------------------------------------------------------------------

      INTEGER   isec0(2), isec1(43), isec2(384), isec3(2)
      INTEGER   isec4(59+kang+kfre),
     &          igrib(kgx*kgy*kang*kfre/jpidl)

      DIMENSION zsec2(96), zsec3(2), psec4(ktgrib*kang*kfre)

      CHARACTER cloper*1

      DATA cloper / 'C' /

!     DATA isec1 /128, 98, 1, 255, 192, 0, 102, 7*0, 1, 5*0, 20, 4*0/
!     DATA isec2 /0, 120, 61, 90000, 0, 128, -90000, 357000,
!    *            3000, 3000, 58*0 /
!     DATA zsec2  /  56*0  /
!     DATA isec3  / 0,  0  /
!     DATA zsec3  / 0., -999. /

!    *ISEC1*  INTEGER  PRODUCT DEFINITION BLOCK.

!       WORD    CONTENTS                                  VALUE.
!       ----    --------                                  ------
!         1     VERSION NUMBER OF CODE TABLE 2.              128
!         2     ORIGINATING CENTRE IDENTIFIER.                98
!         3     MODEL IDENTIFICATION.                          1
!         4     GRID DEFINITION.                             255
!         5     FLAG (WMO CODE TABLE 1)                      128
!         6     PARAMETER IDENTIFIER (WMO CODE TABLE 2)   KPARAM
!        7-8    TYPE OF LEVEL (WMO CODE TABLE 3)             2*0
!         9     VALUE(S) OF LEVEL (WMO CODE TABLE 3)           0
!        10     YEAR OF DATA                                  KY
!        11     MONTH OF DATA                                 KM
!        12     DAY OF DATA                                   KD
!        13     HOUR OF DATA                                  KH
!        14     MINUTE OF DATA                                 0
!        15     TIME UNIT (WMO CODE TABLE 4)                   1
!        16     TIME RANGE ONE (FCST TIME FROM BASE DTG)   KFCST
!        17     TIME RANGE TWO                                 0
!        18     TIME RANGE FLAG (WMO CODE TABLE 5)             0
!        19     NUMBER AVERAGED.                               0
!        20     NUMBER MISSING FROM AVERAGE.                   0
!        21     CENTURY OF REFERENCE TIME                     20
!        22                                                    0
!        23                                                    0
!        24     0=NO LOCAL USE OF SECTION 1; 1=LOCAL USE       1
!        25-36                                                 0
!        37     1=LABELLING WITH MARS PARAMETERS; 0=OTHERWISE  1
!        38     1=CLASS IS OP 2=CLASS IS RD                 1 or 2
!        39     1=TYPE IS FIRST GUESS; 2=TYPE IS ANALYSIS   1,2 or 9
!               9=TYPE IS FORCAST
!        40     1045=STREAM IS WAVE                         1045
!        41     4 ASCII character righ justified version or experiment
!        42,43                                                 0



!    *ISEC2*  INTEGER  GRID DEFINITION BLOCK.

!       WORD    CONTENTS                                  VALUE.
!       ----    --------                                  ------
!         1     DATA REPRESENTATION TYPE (WMO CODE TABLE 6)    0
!         2     NUMBER OF POINTS ALONG A LATITUDE            120
!         3     NUMBER OF POINTS ALONG A MERIDIAN             46
!         4     LATITUDE OF ORIGIN (SOUTH -)               72000
!         5     LONGITUDE OF ORIGIN (WEST -)                   0
!         6     RESOLUTION FLAG (WMO CODE TABLE 7)           128
!         7     LATITUDE OF EXTREME POINT (SOUTH -)       -63000
!         8     LONGITUDE OF EXTREME POINT (WEST -)       357000
!         9     LATITUDE INCREMENT                          3000
!        10     LONGITUDE INCREMENT                         3000
!        11     SCANNING MODE FLAGS (WMO CODE TABLE 8)         0
!       12-68                                               57*0

!    *ZSEC2*  REAL     GRID DEFINITON BLOCK   NOT USED.

!    *ISEC3*  INTEGER  BIT MAP SECTION        NOT USED.

!    *ZSEC3*  REAL     BIT MAP SECTION        NOT USED.

!    *ISEC4*  INTEGER  DATA SECTION.

!       WORD    CONTENTS                                  VALUE.
!       ----    --------                                  ------
!         1     NUMBER OF DATA VALUES IN PSEC4       KTGRIB*kang*kfre
!         2     NUMBER OF BITS FOR EACH PACKED DATA            9
!       3-42                                                40*0

!    *PSEC4*  REAL      UNPACKED DATA FIELD.
!    *ILENP*  INTEGER   LENGTH OF THIS PSEC4.
!    *IGRIB*  INTEGER   ARRAY TO RECEIVE DATA IN GRIB CODE.
!    *IWORD*  INTEGER   NUMBER OF WORDS OCCUPIED BY CODED DATA IN IGRIB.
!    *IRET*   INTEGER   ERROR RETURN FLAG.
!    *CLOPER* CHARACTER REQUESTED FUNCTION

! ----------------------------------------------------------------------

      INTEGER   klonrgg(kgy)

      DIMENSION ptheta(kang), pfreq(kfre)

      REAL :: ONE 
      REAL , DIMENSION(KANG) :: ZTHETA
      REAL , DIMENSION(KFRE) :: ZFREQ

      DIMENSION itheta(kang), ifreq(kfre)

      CHARACTER*1 cddomain, clfdbsf
      CHARACTER*2 cltype
      CHARACTER*2 cdclass
      CHARACTER*4 cdexpver
      CHARACTER*8 clexpverhlp
      LOGICAL ldfdb

      DIMENSION ipr(2)
      DATA ipr / 2*1 /

!     FOR TIMING:
      REAL wam_user_clock,   ! the clock function.
     &     time0,        ! the reference time.
     &     time_g,       ! time spend in gribex.
     &     time_w        ! time used for writing.

! ----------------------------------------------------------------------

!*    1. FIX PARAMETERS AND PACK DATA.
!        -----------------------------

      ONE=1.0
      msub = msub + 1
      CALL  uraopc (csubna, msub, llsuvi, lltest)
      IF (llsuvi .OR. lltest) THEN
        WRITE(kuso,*)' ',csubna,': '
        itest = 2
      ELSE
        itest = 0
      ENDIF

      ibts = KIND(ipr)

      IF (lltest ) THEN
        WRITE(kuso,*) ' SUBROUTINE: GRIBPACKS'
        WRITE(kuso,*) ' KANG  =', kang
        WRITE(kuso,*) ' KFRE  =', kfre
        WRITE(kuso,*) ' KGX   =', kgx
        WRITE(kuso,*) ' KGY   =', kgy
        WRITE(kuso,*) ' KTGRIB=', ktgrib
      ENDIF

!     Parameters for section 1.

!     DATA isec1 /128, 98, 1, 255, 192, 0, 102, 7*0, 1, 5*0, 20, 4*0/

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

!     clexpverhlp(5:8)=cdexpver
!     READ(clexpverhlp,'(A8)') isec1(41)

      isec1(41) = 0
      DO ji = 1, 4
        isec1(41) = isec1(41)*256 + ICHAR(cdexpver(ji:ji))
      ENDDO

      isec1(42) = 0
      isec1(43) = 1
      IF (llsuvi) THEN
        WRITE(kuso,*)' SUB: GRIBPACS - EXPERIMENT VERSION = ', cdexpver
        WRITE(kuso,*)' SUB: GRIBPACS - ISEC1(41) = ', isec1(41)
      ENDIF


!     Set parameters for section 2.

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

      IF (krgg .EQ. 1) THEN
        isec2(1) = 0
        isec2(6) = 128
        isec2(9) = 0
        isec2(17) = 1

        ipadding = 0  ! IPADDING IS JUST FOR CONVIENIENCE TO PRETEND
!                     ! THAT WE DEAL WITH FULL GLOBAL FIELDS.
        istart = 0    ! ISTART IS AN OFFSET USED WHEN THERE ARE 
!                     ! A NUMBER (i.e. istart) OF LATITUDES IN THE
!                     ! BEGINNING OF THE FIELD NOT USED.
        ! istart = NINT((90. - pmonop ) / pdella)
!C  SECTION COMMENTED OUT : IF CHANGED NEED TO CHANGE INMARSS AS WELL !!!
!       IF ( isec2(7) .GE. 0 ) THEN
!         istart = NINT((90. - pmonop ) / pdella)
!         isec2( 4) =   90000
!       ELSEIF ( isec2(4) .LT. 0 ) THEN
!         ipadding = NINT((-90. - pmosop ) / (pdella * (-1.))) -1
!         isec2( 7) =  -90000
!       ENDIF
!C

        isec2( 3) =  knd - kst +1 +istart +ipadding

        IF (lltest) WRITE(kuso,*)'  jy, isn, isec2(22-kst+1+jy+istart),'
     &                 , '  22-kst+1+jy+istart '
        DO jy=kst,knd  !kgy
          isn=kgy-jy+1
          isec2(22-kst+1+jy+istart) = klonrgg(isn)
          IF (lltest) WRITE(kuso,*) jy, isn, isec2(22-kst+1+jy+istart),
     &     22-kst+1+jy+istart
        ENDDO
      ENDIF


!     parameters for section 3.

      isec3(1) = 0
!JH   isec3(2) = 0
!JH   zsec3(1) = 0.0
      zsec3(2) = pmiss

!     parameters for section 4.

      isec4=0
      DO jang=1,kang
        ztheta(jang) = ptheta(jang)*deg
      ENDDO
      zfreq=pfreq
      ilenp = ktgrib*kang*kfre
      isec4( 1) = ilenp
      isec4( 2) =  9 ! 9,16,32 !Number of bits used to encode each value
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

      ij=0
      DO j = 60, 60-1+kang
        ij=ij+1
        itheta(ij) = TRANSFER (ztheta(ij), 1)
        isec4(j) = itheta(ij)
      ENDDO

      ij=0
      DO j = 60+kang, 60+kang-1+kfre
        ij=ij+1
        ifreq(ij)  = TRANSFER (zfreq(ij),  1)
        isec4(j) = ifreq(ij)
      ENDDO

      IF (lltest) THEN
        WRITE(kuso,*) ' ISEC4'
        DO j1=1,59+kang+kfre
          WRITE(kuso,'(i4,1X,I19)')j1,isec4(J1)
        ENDDO

        WRITE(kuso,*) ' ANGLES '
        WRITE(kuso,'(8F8.0)') ztheta
        WRITE(kuso,*) ' FREQUENCIES '
        WRITE(kuso,'(8F8.5)') zfreq

        i1=knd-kst/2
        i2=i1+1
        DO j1=i1, i2
          WRITE(kuso,*) ' SPECTRUM ',j1
          js = ((j1-1)*kang*kfre)+1
          je =  (j1-1)*kang*kfre +kang*kfre
          WRITE(kuso,'(12f8.4)')(psec4(j2),j2=js,je)
        ENDDO
      ENDIF

      kkk = 1
      IF (lltest) THEN
        CALL grsdbg ( kkk )

!                  CALL GRSVCK (I) where
!                              I = 0, to switch off checking.
!                                  Non-zero to switch on checking.

        CALL grsvck (   1 )
      ENDIF

      CALL grsdbg ( kkk )
      CALL grsvck (   1 )

      IF (lltest) WRITE(kuso,*) (isec2(jsec2),jsec2=1,50)
      time0=-wam_user_clock()
      IF (lltest.OR.llsuvi) CALL FLUSH(kuso)
      CALL gribex (isec0, isec1, isec2, zsec2, isec3, zsec3, isec4,
     &             psec4, ilenp, igrib, kang*kfre*kgx*kgy/jpidl, iword,
     &             cloper, iret)

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

      IF (lltest) WRITE(kuso,*)' ', kkk

! ----------------------------------------------------------------------

!*    2. CHECK ERROR CODE.
!        ----------------


      IF (llsuvi)
     &WRITE (kuso,*) '  SUB. GRIBPACS:  GRIBEX RETURN CODE ',iret,
     &              ' FOR PARAMETER ', kparam

      IF (iret .EQ. 0) THEN

!*    2.1 WRITE OUT PACKED DATA.
!        -----------------------

        IF ( ldfdb ) THEN

!*       2.1.1 WRITE FIELD TO FDB.
!              -------------------

          IF (llsuvi) time0=-wam_user_clock()
          idir = 0
          ifre = 0
          clfdbsf = ' '
          CALL FLD2FDB (kuso, itest, isec1, igrib, iword, 'WRITE', ierr,
     &                  cltype, cdclass, cdexpver, clfdbsf, cddomain,
     &                  idir, ifre, kfdb)
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

          WRITE (kuso,*) '    SUB. GRIBPACS:  PACKED FIELD LENGTH ',
     &                  iword, '    WRITTEN TO FIELDS DATA BASE ',
     &                  '(REFERENCE=',kfdb,')'

!*       2.1.2 PRINT FDB ERROR CODE.
!              ---------------------

          IF (ierr .NE. 0) THEN
            WRITE(kuso,*) ' ------------------------'
            WRITE(kuso,*) ' ERROR ACCESSING FDB '
            WRITE(kuso,*) ' FDB ERROR CODE IS ',ierr
            WRITE(kuso,*) ' ------------------------'
            CALL abort
          ENDIF
        ELSE

!*       2.1.3 WRITE FIELD TO FILE.
!              --------------------

          CALL PBWRITE(kugrs, igrib, iword*ibts, iret)
          IF (iret .LT. 0) THEN
            WRITE(kuso,*) 'PBWRITE: IERR = ', ierr
            WRITE(kuso,*) 'Error in writing to grib file'
            CALL abort
          ENDIF

          IF (llsuvi) THEN
            IF (kugrs .GT. 111) THEN
              WRITE (kuso,2200)'  SUB. GRIBPACS:  PACKED FIELD LENGTH ',
     &                     iword, '    WRITTEN TO UNIT KUGRS = ', kugrs
 2200         FORMAT(1x,a,i6,a,i8)
            ELSE
              WRITE (kuso,*) '  SUB. GRIBPACS:  PACKED FIELD LENGTH ',
     &                     iword, '    WRITTEN TO UNIT KUGRS = ', kugrs
            ENDIF
          ENDIF
        ENDIF
      ELSE

!*    2.3 PRINT GRIBEX ERROR CODE.
!         ------------------------

        WRITE(kuso,*) ' -----------------------------------------'
        WRITE(kuso,*) ' ERROR PACKING FIELD BY SUB GRIBEX. '
        WRITE(kuso,*) ' GRIBEX ERROR CODE IS ',iret
        WRITE(kuso,*) ' -----------------------------------------'
        CALL abort
      ENDIF

! ----------------------------------------------------------------------

!*    3. PRINT GRIB SECTION 1 AND 2.
!         --------------------------

      IF (llsuvi) CALL GRPRS1 (isec0, isec1)
      IF (llsuvi) CALL GRPRS2 (isec0, isec2, zsec2)
      IF (llsuvi) WRITE(kuso,*) ' PACKED FIELD LENGTH ',iword
      msub = msub - 1
      RETURN
      END SUBROUTINE gribpacs
