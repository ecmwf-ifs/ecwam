      PROGRAM IGS2GRS
C
C ----------------------------------------------------------------------
C
C**** IGS2GRS*  - RECOMBINE SPLIT GRIB SPECTRA INTO ONE SINGLE GRIB FILE
C                 THE spectra were split in as many files as there were PE's
C                 
C
C     JEAN  BIDLOT   ECMWF    FEBRUARY 1998
C
C*    PURPOSE.
C     --------
C
C**   INTERFACE.
C     ----------
C
C     METHOD.
C     -------
C
C     EXTERNALS.
C     ----------
C     ABORT1
C     FLD2FDB
C     GRSTNAME
C     MPDECOMP
C     READPRE
C     U2L1CR
C     FBDLIB :
C     IREADFDB
C     IWRITEFDB
C     ICLOSEFDB
C     IOPENFDB
C     EMOSLIB:
C     PBOPEN
C     PBSIZE
C     PBGRIB
C     PBCLOSE
C     PBWRITE
C     GRSDBG
C
C     REFERENCE
C     ---------
C     NONE
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "comgrid.h"
C
#include "txtgrid.h"
C
#include "commap.h"
C
#include "txtmap.h"
C
#include "comtext.h"
C
#include "parcons.h"
C
C ----------------------------------------------------------------------
C
      INTEGER I4(2)
      INTEGER,ALLOCATABLE :: NSTARTWAM(:),NENDWAM(:)
      INTEGER,ALLOCATABLE :: KLENTOPWAM(:),KLENBOTWAM(:)
      INTEGER, ALLOCATABLE :: INGRIB(:)
      INTEGER :: KSEC0(2),KSEC1(43),KSEC2(384),KSEC3(2),KSEC4(128)
      INTEGER :: ISEC0(2),ISEC1(43),ISEC2(384),ISEC3(2),ISEC4(128)

C
      REAL X4(2)
      REAL :: PSEC2(96),PSEC3(2)
      REAL :: ZSEC2(96),ZSEC3(2)
      REAL, ALLOCATABLE :: PSEC4(:),PSEC4SUB(:)

C

      CHARACTER*  1 CDDOMAIN
      CHARACTER*2 CDTYPE, CDCLASS
      CHARACTER*3 INPREFIX, OUTPREFIX 
      CHARACTER*4 CLEVEL, CDEXPVER
      CHARACTER*6 CLSTEP
      CHARACTER*10 ZERO, CDTPRO, CDATEF, CDATE, ZDATE 
      CHARACTER*80 CLMSG
      CHARACTER*256 FILENAME, INPATH, OUTPATH, CDFDB2DSP
C
      LOGICAL LLEXIST, LFDBIOIN, LFDBIOOUT
C
      PARAMETER(ZERO=' ')
C
C ----------------------------------------------------------------------
C
      NAMELIST /NALINE/ 
     &          ITEST, NPRWAM, CDTPRO, CDATEF, ISPLIT,
     1          INPREFIX, INPATH, OUTPREFIX, OUTPATH,
     1          LFDBIOIN, LFDBIOOUT,
     1          CDCLASS, CDEXPVER, CDFDB2DSP, CDDOMAIN
C ----------------------------------------------------------------------
C
C*    0. SET DEFAULT VALUES FOR THE NAMELIST ELEMENTS.
C        ---------------------------------------------
C
      ITEST  =   -9
      NPRWAM =    0
      CDTPRO = ZERO
      CDATEF = ZERO
      ISPLIT = 2
      INPREFIX = 'IGS'
      INPATH   = ZERO
      OUTPREFIX = 'GRS'
      OUTPATH  = ZERO
C
      LFDBIOIN  = .FALSE.
      LFDBIOOUT = .FALSE.
C
      CDCLASS  = 'od'
      CDEXPVER = '0001' ! cdexpver is needed as it will be used to
C                         access  the fdb BUT also to set isec1(41) 
      CDFDB2DSP = ' '
      CDDOMAIN = 'g' ! cddomain is used for input. it will be reset for
C                       output depending on isplit.
C
C     1. INITIALISATION
C        --------------

      IU05 = 5
      IU06 = 6
      IU07 = 7

      X4=1.
      NPRECR = N_PRECISION(X4)
      I4=1
      NPRECI = N_PRECISION(I4)
C
C*    2. READ NAMELIST NALINE.
C        ---------------------
C
      READ (IU05, NALINE)
C
C*    3. READ PREPROC OUTPUT.
C        --------------------
C
      CALL READPRE (IU07, 0)
C
C
C*    4. RECOMBINE SPLIT SPECTRA
C        -----------------------
C
C
C     DETERMINE GRID DOMAIN DECOMPOSITION OF THE NPRWAM SUB DOMAINS
C
      ALLOCATE (NSTARTWAM(NPRWAM))
      ALLOCATE (NENDWAM(NPRWAM))
      ALLOCATE (KLENBOTWAM(NPRWAM))
      ALLOCATE (KLENTOPWAM(NPRWAM))
C
      CALL MPDECOMP(NPRWAM,NSTARTWAM,NENDWAM,KLENTOPWAM,
     &              KLENBOTWAM,MAXLEN)
      IF (ITEST.GT.0) THEN
        WRITE (IU06,*) ' SUB. MPDECOMP DONE'
        WRITE(IU06,*)
        WRITE(IU06,*)  ' SUB MODEL DOMAIN DECOMPOSITION  : '
        WRITE(IU06,*)  ' ==============================='
         DO IP=1,NPRWAM
           WRITE(IU06,*)
           WRITE(IU06,*) ' SUB DOMAIN NUMBER : ',IP
           WRITE(IU06,*) ' NSTARTWAM : ',NSTARTWAM(IP)
           WRITE(IU06,*) ' NENDWAM   : ',NENDWAM(IP)
           WRITE(IU06,*) ' ----------------------- '
         END DO
      ENDIF
C
      IF (IRGG.EQ.1) THEN
        KLENTOT=0
        DO K=1,NGY
          JSN=NGY-K+1
          KLENTOT = NLONRGG(JSN) + KLENTOT
        ENDDO
        KLENTOT = KLENTOT * NANG * NFRE
      ELSE
        KLENTOT = NGY * NGX * NANG * NFRE
      ENDIF
      ALLOCATE(PSEC4(KLENTOT))
      PSEC4=ZMISS
C
      CPATH=INPATH
      ICPLEN=LEN_TRIM(CPATH)

      ISIZEN=0
      ISIZES=0
C
      IF(ISPLIT.EQ.1) THEN
        SPLITLAT=AMOSOP
      ELSE
        SPLITLAT=0.0
      ENDIF
C
C     INITIALIZE AND OPEN FDB.
C
      IF(LFDBIOIN.OR.LFDBIOOUT) THEN
        istat = iinitfdb()
        WRITE(iu06,'("  FDB INITIALIZED status=", I3)') istat
      ENDIF
C
C
C     LOOP ON EACH SUB REGION
C
      DO IPR=NPRWAM,1,-1
C
C       GET DATA FROM FILE OR FDB
C
        IF (LFDBIOIN) THEN
C***********************************************************************
C
C          SET FIXED PARAMETERS FOR FDB.

          istat = iopenfdb ('fdb', ifdb, 'r')
          WRITE(iu06,'("  FDB OPENED status=", I3)') istat
          istat = isetvalfdb (ifdb, 'levt',   'm'      )
          istat = isetvalfdb (ifdb, 'param',  '250'    )
          istat = isetvalfdb (ifdb, 'repres', 'll'     )
          istat = isetvalfdb (ifdb, 'stream', 'wave'   )
C---<     istat = isetvalfdb (ifdb, 'fstep',  '000000' )
C
          IF (CDTPRO .LE. CDATEF) THEN
C
C           THIS IS AN ANALYSIS DATE.
C           -------------------------
C
            kfcst = 0
          ELSE
C
C           THIS IS A  FORECAST DATE.
C           -------------------------
C
            CALL DIFDATE (CDATEF, CDTPRO, kfcst)
            kfcst=kfcst/3600
          ENDIF
C
          IF ( kfcst .EQ. 0 ) THEN
            cdtype='an'
            istat = isetvalfdb (ifdb, 'date', '19'//cdtpro(1:6))
            istat = isetvalfdb (ifdb, 'time',   cdtpro(7:8)  )
          ELSE
            cdtype='fc'
            istat = isetvalfdb (ifdb, 'date', '19'//cdatef(1:6))
            istat = isetvalfdb (ifdb, 'time',   cdatef(7:8)  )
          ENDIF
C
          CALL u2l1cr( cdtype )
          CALL u2l1cr( cdclass )
          CALL u2l1cr( cdexpver )
C
          ipe = ipr
          WRITE( clevel, '(i4.4)' ) ipe
          WRITE( clstep, '(i6.6)' ) kfcst
          istat = isetvalfdb (ifdb, 'step',   clstep        )
          istat = isetvalfdb (ifdb, 'level',  clevel        )
          istat = isetvalfdb (ifdb, 'class',  cdclass       )
          istat = isetvalfdb (ifdb, 'expver', cdexpver      )
          istat = isetvalfdb (ifdb, 'type',   cdtype       )
          istat = isetvalfdb (ifdb, 'domain', cddomain )
C         istat = isetrankfdb (ifdb, ipe)

          WRITE(IU06,'("\ /igs2grs/ PARAM: ", a8)') '250'
          WRITE(IU06,'("\ /igs2grs/ STEP:  ", a8)') clstep
          WRITE(IU06,'("\ /igs2grs/ CLASS: ", a8)') cdclass
          WRITE(IU06,'("\ /igs2grs/ EXPVER:", a8)') cdexpver
          WRITE(IU06,'("\ /igs2grs/ DOMAIN:", a8)') cddomain
          WRITE(IU06,'("\ /igs2grs/ TY:    ", a8)') cdtype
          WRITE(IU06,'("\ /igs2grs/ DATE:  ", a8)') cdtpro(1:6)
          WRITE(IU06,'("\ /igs2grs/ TIME:  ", a8)') cdtpro(7:8)
          WRITE(IU06,'("\ /igs2grs/ LEVEL: ", a8)') clevel
C
C* try to estimate the size of the next product:
C* sea points times size spectrum + 10% for land points time bits used
C* for one word.
C
          iblo = ( nendwam(ipr) - nstartwam(ipr) ) * (nang*nfre)
          iblo = iblo + (iblo / 10)
          isize = (iblo * 9)  / 4
          isize = 1129120
          ALLOCATE (ingrib(isize))
          WRITE(IU06,'(" igs2grs: allocate", i10, " words for ingrib",
     .      " nstart=", i8, " nend=", i8, " spectrum=" , i4)
     .               ') isize, nstartwam(ipr), nendwam(ipr), nang*nfre

          irstat = ireadfdb (ifdb, ingrib, isize      )

          WRITE(IU06,'(" igs2grs: readfdb read ",i10,
     .    " words (status:", i3, ")" )')isize, irstat
          IF (irstat.NE.0 .OR. isize.LE.0) THEN
            WRITE (IU06,'("***************************************")')
            WRITE (IU06,'("*                                     *")')
            WRITE (IU06,'("*   ERROR FOLLOWING CALL TO ireadfdb  *")')
            WRITE (IU06,'("*   IN igs2grs                        *")')
            WRITE (IU06,'("*   STAT CODE irstat =", I10)') irstat
            WRITE (IU06,'("*   LENGTH OF MESSAGE=", I10)') isize
            WRITE (IU06,'("*                                     *")')
            WRITE (IU06,'("***************************************")')
            CALL ABORT1
          ENDIF
C
          iclose = iclosefdb (ifdb)
C
C***********************************************************************
        ELSE
          CALL GRSTNAME(IU06,CDTPRO,CDATEF,2,INPREFIX,IPR,NPRWAM,
     &                  FILENAME)
          LFILE=IECF_LEN(FILENAME)
          INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)
          IF(LLEXIST) THEN
            CALL PBOPEN(IUGRS,FILENAME(1:LFILE),'r',KRET)
            IF(KRET.LT.0) GOTO 6100
          ELSE
            WRITE(IU06,*)'*******************************************'
            WRITE(IU06,*)'*                                         *'
            WRITE(IU06,*)'*IN IGS2GRS : GRIB SPECTRA NOT FOUND IN   *'
            WRITE(IU06,*)  FILENAME
            WRITE(IU06,*)'*PROGRAM WILL ABORT                       *'
            WRITE(IU06,*)'*                                         *'
            WRITE(IU06,*)'*******************************************'
            CALL ABORT1
          ENDIF
C
          CALL PBSIZE (IUGRS, IPLENG)
          ISIZE=(IPLENG+NPRECI-1)/NPRECI
          ALLOCATE(INGRIB(ISIZE))
C
          CALL PBGRIB(IUGRS,INGRIB,IPLENG,ILENG,KRET)
          IF (KRET.EQ.-1) THEN
            WRITE (IU06,*) ' PROGRAM IGS2GRS : '
            WRITE (IU06,*) ' REACHED EOF ON UNIT ',IUGRS
            CALL ABORT1
          ELSEIF (KRET .EQ. -2) THEN
            WRITE (IU06,*) ' PROGRAM IGS2GRS : '
            WRITE (IU06,*) ' ERROR IN FILE HANDLING ON UNIT ',IUGRS
            CALL ABORT1
          ELSEIF (KRET .EQ. -3) THEN
            WRITE (IU06,*) ' PROGRAM IGS2GRS : '
            WRITE (IU06,*) ' INGRIB IS NOT BIG ENOUGH.'
            CALL ABORT1
          ENDIF
C
          CALL PBCLOSE(IUGRS,KRET)
          IF(KRET.LT.0) GOTO 6400
        ENDIF
C
C       FIND THE NORTH AND SOUTH BOUNDARY OF EACH SUB REGION
C
        IKN = NGY-KXLT(NENDWAM(IPR),1)+1
        IKS = NGY-KXLT(NSTARTWAM(IPR),1)+1

        SMONOP = AMONOP-(IKN-1)*XDELLA
        SMOSOP = AMONOP-(IKS-1)*XDELLA
        KNGY = IKS-IKN+1
C
        IF(SMONOP.GE.SPLITLAT) ISIZEN=ISIZEN+ISIZE
        IF(SMOSOP.LT.SPLITLAT) ISIZES=ISIZES+ISIZE
C
C       DETERMINE THE RELATIVE POSSITION OF THE SUBREGION DATA
C       IN GLOBAL ARRAY, AND THE NUMBER OF DATA IN THAT SUBREGION
C
        IF (IRGG.EQ.1) THEN
          IOFFSET=0
          DO K=1,IKN-1
            JSN=NGY-K+1
            IOFFSET = NLONRGG(JSN) + IOFFSET
          ENDDO
          IOFFSET = IOFFSET * NANG * NFRE
          NTOT=0
          DO K=IKN,IKS
            JSN=NGY-K+1
            NTOT = NLONRGG(JSN) + NTOT
          ENDDO
          NTOT = NTOT * NANG * NFRE
        ELSE
          IOFFSET = (NGY-IKN) * NGX * NFRE * NANG
          NTOT = KNGY * NGX * NANG * NFRE
        ENDIF

        ALLOCATE(PSEC4SUB(NTOT))
        PSEC4SUB=ZMISS
C
C       DECODE THE SUBREGION GRIB DATA
C
        KSEC3(2)=0
        PSEC3(2)=ZMISS
C
        CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &               PSEC4SUB, NTOT, INGRIB, ISIZE, KWORD,
     &               'D', KRET)
        IF(KRET.GT.0) THEN
          WRITE(IU06,*) '************************************'
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '* FATAL ERROR IN PGRM IGS2GRS      *'
          WRITE(IU06,*) '* GRIBEX ERROR - KRET = ',KRET
          WRITE(IU06,*) '* WHILE DECODING INGRB             *'
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '*    PROGRAM ABORTS                *'
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '************************************'
          CALL ABORT1
        ENDIF
C
        ISEC0 = KSEC0
        ISEC1 = KSEC1
        ISEC2 = KSEC2
        ISEC3 = KSEC3
        ISEC4 = KSEC4
        ZSEC2 = PSEC2
        ZSEC3 = PSEC3
C
C       CHECK THAT INPUT DESCRIPTION IS CONSISTENT
C
        IF (KSEC1(6).NE.250) THEN
           WRITE(IU06,*)'********************************'
           WRITE(IU06,*)'*                              *'
           WRITE(IU06,*)'* FATAL ERROR IN PRGM IGS2GRS  *'
           WRITE(IU06,*)'* ===========================  *'
           WRITE(IU06,*)'*                              *'
           WRITE(IU06,*)'* GRIB PARAMETER  ',KSEC1(6)
           WRITE(IU06,*)'* WAS READ INSTEAD OF 250' 
           WRITE(IU06,*)'*                              *'
           WRITE(IU06,*)'********************************'
           CALL ABORT1
        ENDIF

        WRITE (ZDATE,'(5I2.2)') KSEC1(10),KSEC1(11),KSEC1(12),
     .                          KSEC1(13),KSEC1(14)
        IF(CDTPRO.GT.CDATEF) THEN 
          CDATE=CDATEF
        ELSE
          CDATE=CDTPRO
        ENDIF
        IF (ZDATE.NE.CDATE) THEN
           WRITE(IU06,*)'**********************************'
           WRITE(IU06,*)'*                                *'
           WRITE(IU06,*)'* FATAL ERROR IN PRGM IGS2GRS    *'
           WRITE(IU06,*)'* ===========================    *'
           WRITE(IU06,*)'*                                *'
           WRITE(IU06,*)'* REQUESTED DATE IS NOT EQUAL TO *'
           WRITE(IU06,*)'* RETRIEVED DATE.                *'
           WRITE(IU06,*)'* ZDATE = ',ZDATE
           WRITE(IU06,*)'* CDATE = ',CDATE
           WRITE(IU06,*)'*                                *'
           WRITE(IU06,*)'**********************************'
           CALL ABORT1
        ENDIF
C
        ZMONOP = FLOAT(KSEC2(4)/1000)+0.1*(MOD(KSEC2(4),1000)/100)+
     &           0.01*(MOD(KSEC2(4),100)/10)+0.001*MOD(KSEC2(4),10)

        ZMOWEP = FLOAT(KSEC2(5)/1000)+0.1*(MOD(KSEC2(5),1000)/100)+
     &           0.01*(MOD(KSEC2(5),100)/10)+0.001*MOD(KSEC2(5),10)

        ZMOSOP = FLOAT(KSEC2(7)/1000)+0.1*(MOD(KSEC2(7),1000)/100)+
     &           0.01*(MOD(KSEC2(7),100)/10)+0.001*MOD(KSEC2(7),10)

        ZMOEAP = FLOAT(KSEC2(8)/1000)+0.1*(MOD(KSEC2(8),1000)/100)+
     &           0.01*(MOD(KSEC2(8),100)/10)+0.001*MOD(KSEC2(8),10)

        IF (AMOWEP.NE.ZMOWEP .OR. AMOEAP.NE.ZMOEAP .OR.
     +      SMONOP.NE.ZMONOP .OR. SMOSOP.NE.ZMOSOP) THEN
            WRITE(IU06,*) '***********************************'
            WRITE(IU06,*) '*                                 *'
            WRITE(IU06,*) '*  SUB REGION NOT RECOGNISED      *'
            WRITE(IU06,*) '*  ==========================     *'
            WRITE(IU06,*) '*                                 *'
            WRITE(IU06,*) '*  THIS IS NOT THE SAME GRID AREA *'
            WRITE(IU06,*) '*  AS DELARED                     *'
            WRITE(IU06,*) '*                                 *'
            WRITE(IU06,*) '*  AMOSOP: ', SMOSOP, 'ZMOSOP : ',ZMOSOP
            WRITE(IU06,*) '*  AMONOP: ', SMONOP, 'ZMONOP : ',ZMONOP
            WRITE(IU06,*) '*  AMOWEP: ', AMOWEP, 'ZMOWEP : ',ZMOWEP
            WRITE(IU06,*) '*  AMOEAP: ', AMOEAP, 'ZMOEAP : ',ZMOEAP
            WRITE(IU06,*) '*                                 *'
            WRITE(IU06,*) '***********************************'
            CALL ABORT1
        ENDIF
        INANG = KSEC4(50)
        IF (INANG .NE. NANG) THEN
           WRITE(IU06,*)'**********************************'
           WRITE(IU06,*)'*                                *'
           WRITE(IU06,*)'* FATAL ERROR IN  IGS2GRS        *'
           WRITE(IU06,*)'* ===========================    *'
           WRITE(IU06,*)'*                                *'
           WRITE(IU06,*)'* INANG SHOULD BE EQUAL TO NANG  *'
           WRITE(IU06,*)'* NANG  = ',NANG
           WRITE(IU06,*)'* INANG  = ',INANG
           WRITE(IU06,*)'*                                *'
           WRITE(IU06,*)'**********************************'
           CALL ABORT1
        ENDIF
        INFRE =  KSEC4(51)
        IF (INFRE .NE. NFRE) THEN
           WRITE(IU06,*)'**********************************'
           WRITE(IU06,*)'*                                *'
           WRITE(IU06,*)'* FATAL ERROR IN IGS2GRS         *'
           WRITE(IU06,*)'* ===========================    *'
           WRITE(IU06,*)'*                                *'
           WRITE(IU06,*)'* INFRE SHOULD BE EQUAL TO NFRE  *'
           WRITE(IU06,*)'* NFRE  = ',NFRE
           WRITE(IU06,*)'* INFRE  = ',INFRE
           WRITE(IU06,*)'*                                *'
           WRITE(IU06,*)'**********************************'
           CALL ABORT1
        ENDIF
        INGY = KSEC2(3)
        IF (INGY .NE. KNGY) THEN
           WRITE(IU06,*)'********************************'
           WRITE(IU06,*)'*                              *'
           WRITE(IU06,*)'* FATAL ERROR IN IGS2GRS       *'
           WRITE(IU06,*)'* ERROR IN GRID REPRESENTATION *' 
           WRITE(IU06,*)'* ===========================  *'
           WRITE(IU06,*)'*                              *'
           WRITE(IU06,*)'* INGY IS NOT EQUAL TO KNGY    *'
           WRITE(IU06,*)'* KNGY  = ',KNGY
           WRITE(IU06,*)'* INGY  = ',INGY
           WRITE(IU06,*)'*                              *'
           WRITE(IU06,*)'********************************'
           CALL ABORT1
        ENDIF
        KRGG = KSEC2(17)
        IF (KRGG.NE.IRGG) THEN
           WRITE(IU06,*)'********************************'
           WRITE(IU06,*)'*                              *'
           WRITE(IU06,*)'* FATAL ERROR IN IGS2GRS       *'
           WRITE(IU06,*)'* ERROR IN GRID REPRESENTATION *' 
           WRITE(IU06,*)'* ============================ *'
           WRITE(IU06,*)'*                              *'
           WRITE(IU06,*)'* KRGG IS NOT EQUAL TO IRGG    *'
           WRITE(IU06,*)'* KRGG  = ',KRGG
           WRITE(IU06,*)'* IRGG  = ',IRGG
           WRITE(IU06,*)'*                              *'
           WRITE(IU06,*)'********************************'
           CALL ABORT1
        ENDIF
        IF(IRGG.EQ.1) THEN
          IC=0
          DO J=IKN,IKS
            JSN=NGY-J+1
            IC=IC+1
            IF( KSEC2(22+IC).NE.NLONRGG(JSN)) THEN
               WRITE(IU06,*)'***********************************'
               WRITE(IU06,*)'*                                 *'
               WRITE(IU06,*)'* FATAL ERROR IN IGS2GRS          *'
               WRITE(IU06,*)'* INCONSISTENT IRREGULAR GRID     *'
               WRITE(IU06,*)'* ===========================     *'
               WRITE(IU06,*)'*                                 *'
               WRITE(IU06,*)'* THE IRREGULAR GRID DEFINITIONS  *'
               WRITE(IU06,*)'* OF THE INPUT FILE DO NOT MATCH  *'
               WRITE(IU06,*)'* WITH THE MODEL ONE.             *'
               WRITE(IU06,*)'* KDOMRGG IS NOT EQUAL TO NLONRGG *'
               IIC=0
               DO JJ=IKN,IKS
                 JSN=NGY-JJ+1
                 IIC=IIC+1
                 WRITE(IU06,*)'* JSN = ',JSN
                 WRITE(IU06,*)'* KDOMRGG  = ',KSEC2(22+IIC)
                 WRITE(IU06,*)'* NLONRGG  = ',NLONRGG(JSN)
               ENDDO
               WRITE(IU06,*)'*                                 *'
               WRITE(IU06,*)'***********************************'
               CALL ABORT1
            ENDIF
          ENDDO
        ELSE 
          INGX = KSEC2(2)
          IF (INGX .NE. NGX) THEN
             WRITE(IU06,*)'********************************'
             WRITE(IU06,*)'*                              *'
             WRITE(IU06,*)'* FATAL ERROR IN IGS2GRS       *'
             WRITE(IU06,*)'* INCONSISTENT REGULAR GRID    *'
             WRITE(IU06,*)'* ===========================  *'
             WRITE(IU06,*)'*                              *'
             WRITE(IU06,*)'* INGX IS NOT EQUAL TO NGX     *'
             WRITE(IU06,*)'* NGX  = ',NGX
             WRITE(IU06,*)'* INGX  = ',INGX
             WRITE(IU06,*)'*                              *'
             WRITE(IU06,*)'********************************'
             CALL ABORT1
          ENDIF
        ENDIF
C
C       TRANSFER PSEC4SUB INTO PSEC4
C
C       MOST NORTHERN LATITUDE OF SUB REGION
        IXN = IXLG(NENDWAM(IPR),1)

        ILEN = IXN*NANG*NFRE
        DO IC=1,ILEN
           PSEC4(IOFFSET+IC)=PSEC4SUB(IC)
        ENDDO
C
C       SUB REGION MID SECTION
        IF (IRGG.EQ.1) THEN
          ILEN = 0 
          DO K=IKN+1,IKS-1
            JSN=NGY-K+1
            ILEN = NLONRGG(JSN) + ILEN
          ENDDO
          ILEN = ILEN * NANG * NFRE
          IOFFSUB = NLONRGG(NGY-IKN+1) * NANG * NFRE
        ELSE
          ILEN = (IKS-IKN-1) * NANG * NFRE
          IOFFSUB = NGX * NANG * NFRE
        ENDIF
C
        DO IC=IOFFSUB+1,IOFFSUB+ILEN
           PSEC4(IOFFSET+IC)=PSEC4SUB(IC)
        ENDDO
C
C       MOST SOUTHERN LATITUDE OF SUB REGION
        IXS = IXLG(NSTARTWAM(IPR),1)
        IOFFSUB = IOFFSUB + ILEN + (IXS-1) * NANG * NFRE
        IF (IRGG.EQ.1) THEN
          ILEN = (NLONRGG(NGY-IKS+1)-IXS+1) * NANG * NFRE
        ELSE
          ILEN = (NGX-IXS+1) * NANG * NFRE
        ENDIF
C
        DO IC=IOFFSUB+1,IOFFSUB+ILEN
           PSEC4(IOFFSET+IC)=PSEC4SUB(IC)
        ENDDO
C
C
        IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)
        IF(ALLOCATED(PSEC4SUB)) DEALLOCATE(PSEC4SUB)

      ENDDO

      IF(ALLOCATED(NSTARTWAM)) DEALLOCATE(NSTARTWAM)
      IF(ALLOCATED(NENDWAM)) DEALLOCATE(NENDWAM)
      IF(ALLOCATED(KLENBOTWAM)) DEALLOCATE (KLENBOTWAM)
      IF(ALLOCATED(KLENTOPWAM)) DEALLOCATE (KLENTOPWAM)
C
C
C*    5. PACK AND SAVE TOTAL SPECTRA
C        ---------------------------
C
      CPATH=OUTPATH
      ICPLEN=LEN_TRIM(CPATH)
      CALL GRSTNAME(IU06,CDTPRO,CDATEF,1,OUTPREFIX,0,0,FILENAME)
      LFILE=IECF_LEN(FILENAME)

      IF (LFDBIOOUT) THEN
        CONTINUE
      ELSE
        WRITE(IU06,*)'RECOMBINED SPECTRA WILL BE SAVED IN '
        WRITE(IU06,*) FILENAME(1:LFILE)
        WRITE(IU06,*) ' ' 
        CALL PBOPEN(IUOUT,FILENAME(1:LFILE),'w',KRET)
        IF(KRET.LT.0) GOTO 6100
      ENDIF

      K=1
      DO WHILE (AMONOP-(K-1)*XDELLA.GT.SPLITLAT)
        K=K+1
      ENDDO
      IKSPLIT=MIN(K,NGY)

      IF (IRGG.EQ.1) THEN
        KLENPN=0
        DO K=1,IKSPLIT
          JSN=NGY-K+1
          KLENPN=KLENPN+NLONRGG(JSN)
        ENDDO
        KLENPN = KLENPN * NANG * NFRE

        KLENPS=0
        DO K=IKSPLIT+1,NGY
          JSN=NGY-K+1
          KLENPS = NLONRGG(JSN) + KLENPS
        ENDDO
        KLENPS = KLENPS * NANG * NFRE
      ELSE
        KLENPN= IKSPLIT *  NGX * NFRE * NANG
        KLENPS = (NGY-IKSPLIT) * NGX * NFRE * NANG
      ENDIF
C
C     FIND PSEC4 MAXIMUM 
      PMAX=-3
      DO IC=1,KLENTOT
         PMAX=MAX(PMAX,PSEC4(IC))
      ENDDO

      IF(ITEST.GT.0) CALL GRSDBG(1)
C
C     BOTH HEMISPHERES OR NORTHERN HEMISPHERE
C
      CALL u2l1cr( cdclass )
      if(cdclass.eq.'od') then
        isec1(38)=1
      elseif(cdclass.eq.'rd') then
        isec1(38)=2
      elseif(cdclass.eq.'er') then
        isec1(38)=3
      else
        write(iu06,*) ' PROGRAM IGS2GRS '
        write(iu06,*) ' class not supported ',isec1(38)
        write(iu06,*) ' PROGRAM ABORTS '
        call abort1
      endif
c
      READ(cdexpver,'(A4)') ISEC1(41)
c
      SMONOP = AMONOP
      SMOSOP = SPLITLAT
      ISEC2( 3) =  IKSPLIT 
      ISEC2( 4) =  NINT(SMONOP*1000)
      ISEC2( 5) =  NINT(AMOWEP*1000)
      ISEC2( 7) =  NINT(SMOSOP*1000)
      ISEC2( 8) =  NINT(AMOEAP*1000)
      IF (IRGG.EQ.1) THEN
        DO K=1,IKSPLIT
          JSN=NGY-K+1
          ISEC2(22+K) = NLONRGG(JSN)
        ENDDO
      ELSE
        DO K=1,IKSPLIT
          ISEC2(22+K) = NGX 
        ENDDO
      ENDIF
      ISEC4(1) = KLENPN 
C
      ALLOCATE(INGRIB(ISIZEN))
C
C     SET MAXIMUM FOR CODING
      CALL GRSMAX(PMAX)
C
      IF (ITEST.GT.0) CALL GRSDBG (1) 
      CALL GRIBEX (ISEC0, ISEC1, ISEC2, ZSEC2, ISEC3, ZSEC3, ISEC4,
     &             PSEC4, KLENPN, INGRIB, ISIZEN , KWORD, 'C', KRET)
      CLMSG='CODING THE FULL FIELD OR THE NORTHERN HEMISPHERE'
      IF(KRET.GT.0) GOTO 6200
C
C     OUTPUT
C
      IF (LFDBIOOUT) THEN
        if(isec1(39).eq.1) then
          cdtype='fg' 
        elseif(isec1(39).eq.2) then
          cdtype='an' 
        elseif(isec1(39).eq.9) then
          cdtype='fc'
        elseif(isec1(39).eq.10) then
          cdtype='cf'
        elseif(isec1(39).eq.11) then
          cdtype='pf'
        else
          write(iu06,*) ' PROGRAM IGS2GRS '
          write(iu06,*) ' type not supported ',isec1(39)
          write(iu06,*) ' PROGRAM ABORTS '
          call abort1
        endif
c
        IF(ISPLIT.EQ.2) then
          cddomain='n' 
        ELSE
          cddomain='g'
        ENDIF
C
        cdfdb2dsp = ' '
C
        CALL FLD2FDB (IU06, itest, isec1, ingrib, kword, 'write', ierr,
     .                cdtype, cdclass, cdexpver, cdfdb2dsp, cddomain,
     .                0, 0, kfdb)
        IF (ierr .NE. 0) THEN
          WRITE(IU06,*) ' ------------------------'
          write(iu06,*) ' PROGRAM IGS2GRS '
          WRITE(IU06,*) ' ERROR ACCESSING FDB '
          WRITE(IU06,*) ' FDB ERROR CODE IS ',ierr
          WRITE(IU06,*) ' ------------------------'
          CALL abort1
        ENDIF
        IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. FLD2FDB FOR NH DONE'
        ISTAT = iflushfdb (kfdb)
        CALL FLD2FDB (IU06, itest, isec1, outgrib, kount, 'close', kerr,
     .                cdtype, cdclass, cdexpver, cdfdb2dsp, cddomain,
     .                0 ,0 , kfdb)
        IF (ierr .NE. 0) THEN
          WRITE(IU06,*) ' ------------------------'
          write(iu06,*) ' PROGRAM IGS2GRS '
          WRITE(IU06,*) ' ERROR ACCESSING FDB '
          WRITE(IU06,*) ' FDB ERROR CODE IS ',ierr
          WRITE(IU06,*) ' ------------------------'
          CALL abort1
        ENDIF
      ELSE
        KOUNT=KWORD*NPRECI
        CALL PBWRITE(IUOUT,INGRIB,KOUNT,KRET)
        IF(KRET.LT.0) GOTO 6300
      ENDIF

      IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)
C
C     SOUTHERN HEMISPHERE
C
      IF(ISPLIT.EQ.2) THEN
        CALL U2L1CR( CDCLASS )
        if(cdclass.eq.'od') then
          isec1(38)=1
        elseif(cdclass.eq.'rd') then
          isec1(38)=2
        elseif(cdclass.eq.'er') then
          isec1(38)=3
        else
          write(iu06,*) ' PROGRAM IGS2GRS '
          write(iu06,*) ' class not supported ',isec1(38)
          WRITE(IU06,*) ' PROGRAM ABORTS '
          CALL ABORT1
        ENDIF
        READ(cdexpver,'(A4)') ISEC1(41)
        SMONOP = SPLITLAT - XDELLA
        SMOSOP = AMOSOP
        ISEC2( 3) =  NGY - IKSPLIT 
        ISEC2( 4) =  NINT(SMONOP*1000)
        ISEC2( 5) =  NINT(AMOWEP*1000)
        ISEC2( 7) =  NINT(SMOSOP*1000)
        ISEC2( 8) =  NINT(AMOEAP*1000)
        IC=1
        IF (IRGG.EQ.1) THEN
          DO K=IKSPLIT+1,NGY
            JSN=NGY-K+1
            ISEC2(22+IC) = NLONRGG(JSN)
            IC=IC+1
          ENDDO
        ELSE
          DO K=IKSPLIT+1,NGY
            ISEC2(22+IC) = NGX 
            IC=IC+1
          ENDDO
        ENDIF
        ISEC4(1) = KLENPS 
C
        ALLOCATE(INGRIB(ISIZES))
C
C       SET MAXIMUM FOR CODING
        CALL GRSMAX(PMAX)
C
        CALL GRIBEX (ISEC0, ISEC1, ISEC2, ZSEC2, ISEC3, ZSEC3, ISEC4,
     &               PSEC4(KLENPN+1), KLENPS, INGRIB, ISIZES , KWORD,
     &               'C', KRET)
        CLMSG='CODING THE SOUTHERN HEMISPHERE'
        IF(KRET.GT.0) GOTO 6200
C
C       OUTPUT
C
        IF (LFDBIOOUT) THEN
          if(isec1(39).eq.1) then
            cdtype='fg' 
          elseif(isec1(39).eq.2) then
            cdtype='an' 
          elseif(isec1(39).eq.9) then
            cdtype='fc'
          elseif(isec1(39).eq.10) then
            cdtype='cf'
          elseif(isec1(39).eq.11) then
            cdtype='pf'
          else
            write(iu06,*) ' PROGRAM IGS2GRS '
            write(iu06,*) ' type not supported ',isec1(39)
            write(iu06,*) ' PROGRAM ABORTS '
           call abort1
          endif
          if(isec1(38).eq.1) then
            cdclass='OD'
          elseif(isec1(38).eq.2) then
            cdclass='RD'
          elseif(isec1(38).eq.3) then
            cdclass='ER'
          else
            write(iu06,*) ' PROGRAM IGS2GRS '
            write(iu06,*) ' class not supported ',isec1(38)
            write(iu06,*) ' PROGRAM ABORTS '
           call abort1
          endif
          cddomain='s' 
          cdfdb2dsp = ' '
C
          CALL FLD2FDB (IU06, itest, isec1, ingrib, kword, 'write',
     .                  ierr,
     .                  cdtype, cdclass, cdexpver, cdfdb2dsp, cddomain,
     .                  0, 0, kfdb)
          IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. FLD2FDB FOR SH DONE'
          ISTAT = iflushfdb (kfdb)
          CALL FLD2FDB (IU06, itest, isec1, idum, idum, 'close',
     .                  ierr,
     .                  cdtype, cdclass, cdexpver, cdfdb2dsp, cddomain,
     .                  0, 0, kfdb)
          IF (ierr .NE. 0) THEN
            WRITE(IU06,*) ' ------------------------'
            write(iu06,*) ' PROGRAM IGS2GRS '
            WRITE(IU06,*) ' ERROR ACCESSING FDB '
            WRITE(IU06,*) ' FDB ERROR CODE IS ',ierr
            WRITE(IU06,*) ' ------------------------'
            CALL abort1
          ENDIF

        IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. FLD2FDB close DONE'

        ELSE
          KOUNT=KWORD*NPRECI
          CALL PBWRITE(IUOUT,INGRIB,KOUNT,KRET)
          IF(KRET.LT.0) GOTO 6300
C
          CALL PBCLOSE(IUOUT,KRET)
          IF(KRET.LT.0) GOTO 6400
        ENDIF
C
        IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)
      ENDIF
C
      IF(ALLOCATED(PSEC4)) DEALLOCATE(PSEC4)
C
      STOP 
C
C     6. ERROR MESSAGES
C
6100  CONTINUE
      WRITE (IU06,*) '*************************************'
      WRITE (IU06,*) '*                                   *'
      WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBOPEN  *'
      WRITE (IU06,*) '*   IN IGS2GRS                      *'
      IF(KRET.EQ.-1)
     &        WRITE (IU06,*) ' COULD NOT OPEN FILE ',FILENAME
      IF(KRET.EQ.-2)
     &        WRITE (IU06,*) ' INVALID FILENAME ',FILENAME
      IF(KRET.EQ.-3)
     &        WRITE (IU06,*)' INVALID OPEN MODE SPECIFIED'
      WRITE (IU06,*) '*                                   *'
      WRITE (IU06,*) '*************************************'
      CALL ABORT1
6200  CONTINUE
      WRITE(IU06,*) '*********************************************'
      WRITE(IU06,*) '*                                           *'
      WRITE(IU06,*) '*         FATAL ERROR IN IGS2GRS            *'
      WRITE(IU06,*) '*         ============================      *'
      WRITE(IU06,*) '* GRIBEX ERROR - KRET = ',KRET
      WRITE(IU06,*) '* WHILE ',CLMSG
      WRITE(IU06,*) '*                                           *'
      WRITE(IU06,*) '*    PROGRAM ABORTS     PROGRAM ABORTS      *'
      WRITE(IU06,*) '*                                           *'
      WRITE(IU06,*) '*********************************************'
      CALL ABORT1
6300  CONTINUE
      WRITE (IU06,*) '**************************************'
      WRITE (IU06,*) '*                                    *'
      WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO PBWRITE   *'
      WRITE (IU06,*) '*  IN IGS2GRS                        *'
      WRITE (IU06,*) '*  FILE ',FILENAME
      WRITE (IU06,*) '*                                    *'
      WRITE (IU06,*) '**************************************'
      CALL ABORT1
6400  CONTINUE
      WRITE (IU06,*) '************************************'
      WRITE (IU06,*) '*                                  *'
      WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
      WRITE (IU06,*) '* IN IGS2GRS                       *'
      WRITE (IU06,*) '* FILE ',FILENAME
      WRITE (IU06,*) '*                                  *'
      WRITE (IU06,*) '************************************'
      CALL ABORT1
C
      END
