      SUBROUTINE GETSPEC(FL,NSTART,NEND,IREAD)
C ----------------------------------------------------------------------
C     J. BIDLOT    ECMWF      SEPTEMBER 1997 
C     B. HANSEN    ECMWF      FEBRUARY  1998  UTILIS FDB TO INPUT 2DSP.
C          ATTN!! NEW CHARACTER VARIABLE CLDOMAIN DERIVED FROM
C                 PREPROCESSOR DEFINITION region.
C                 SPECIFY 'region' AS FOLLOWS:
C                 -Dregion=\'g\' IF THE DOMAIN SHOULD BE GLOBAL.
C                 -Dregion=\'m\' IF THE DOMAIN SHOULD BE MEDITE.
C
C*    PURPOSE.
C     --------
C     READ THE SPECTRA FROM DISK.
C
C**   INTERFACE.
C     ----------
C     *CALL* *GETSPEC(FL,NSTART,NEND,IREAD)
C     *FL*        ARRAY CONTAINING THE SPECTRA CONTRIBUTION ON EACH PE
C     *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
C     *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
C     *IREAD*     PROCESSOR WHICH WILL ACCESS THE FILE ON DISK 
C
C     METHOD.
C     -------
C     USE READFL TO READ IN THE SPECTRA DEPENDING ON THE USE OF THE
C     PBIO SOFTWARE OR NOT. PBIO WILL LIMIT THE SIZE OF THE ARRAY
C     NECESSARY TO READ THE INPUT SPECTRA. THE READING IS ONLY DONE ON
C     PE 1, THEREFORE THE RELEVANT INFORMATION IS SENT TO THE OTHER
C     PE'S USING MPDISTRIBFL
C
C     EXTERNALS.
C     ----------
C     GETENV
C     GRIBEX
C     GRSTNAME
C     MPDISTRIBFL
C     MPDISTRIBSCFLD
C     MPE_BARRIER
C     PBOPEN
C     PBGRIB
C     PBCLOSE 
C     READFL
C
C     REFERENCE.
C     ----------
C     NONE
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "parcons.h"
C
#include "comcout.h"
C
#include "comfred.h"
C
#include "comgrid.h"
C
#include "commap.h"
C
#include "commesspass.h"
C
#include "commpp.h"
C
#include "comstat.h"
C
#include "comtest.h"
C
#include "comtext.h"
C
      CHARACTER*  1 MODE
      CHARACTER*  9 CMBXSIZE
      CHARACTER* 80 FILENAME
C
C* fdb variables.
C
      CHARACTER*  4 CLEVEL, CLEXPVER
      CHARACTER*  2 CLTYPE, CLCLASS, CLDIR, CLFRE
      CHARACTER*  1 CLDOMAIN
      DATA CLDOMAIN  / region /
C
      LOGICAL LFRSTIME, LPBOPEN, LPBCLOSE, LLEXIST
C
      INTEGER ifdb
      SAVE ifdb
C
C ----------------------------------------------------------------------
C
      INTEGER,DIMENSION(NPROC) :: NSTART, NEND
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL
C
C ----------------------------------------------------------------------
C     define local allocatable arrays
C
      INTEGER,ALLOCATABLE :: INGRIB(:), IPOINTER(:)
      INTEGER,ALLOCATABLE :: KSEC0(:),KSEC1(:),KSEC2(:),KSEC3(:)
      INTEGER,ALLOCATABLE :: KSEC4(:)
      REAL, ALLOCATABLE,DIMENSION(:) :: PSEC2,PSEC3,PSEC4,WORK,WORKRCV
C
      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: RFL
      REAL :: ONE,ZTHETA,ZFRE
C
      DATA MDEL, LFRSTIME, ONE / 1, .TRUE., 1.0/
C
      LPBOPEN = .TRUE.
      LPBCLOSE = .TRUE.
      ONE=1.0
C
      IF(LMESSPASS) THEN
C
        IF(LFRSTIME.AND.(.NOT.LGRIBIN.OR.LRESTARTED)) THEN
C
C       DETERMINE MAXIMUM NUMBER OF FREQUENCIES WHICH CAN BE USED
C       TO SEND SPECTRA FROM ONE PE  TO THE OTHERS
C       WITHOUT EXCEEDING THE MAXIMUM MAILBOX SIZE 
C
          CALL GETENV('VPP_MBX_SIZE',CMBXSIZE)
          READ(CMBXSIZE,'(I9)') MBXSIZE
          MINSIZE=NINT(FLOAT(NPROC-1)*NIBLO*NANG*NPRECR/NPROC)
          IF(MINSIZE.GT.MBXSIZE) THEN
            WRITE (IU06,*) '***************************************'
            WRITE (IU06,*) '*                                     *'
            WRITE (IU06,*) '* MESSAGE MAIL BOX SIZE IS TOO SMALL  *'
            WRITE (IU06,*) '* TO DISTRIBUTE THE RESTART SPECTRUM  *'
            WRITE (IU06,*) '* FILE                                *'
            WRITE (IU06,*) '* THE PRESCRIBED SIZE SIZE IS         *'
            WRITE (IU06,*) '* ', MBXSIZE,' BYTES                  *'
            WRITE (IU06,*) '* PLEASE INCREASE SIZE TO AT LEAST    *'
            WRITE (IU06,*) '* ', MINSIZE, 'BYTES                  *'
            WRITE (IU06,*) '* SEE ENVIRONMENT VARIABLE            *'
            WRITE (IU06,*) '* VPP_MBX_SIZE                        *' 
            WRITE (IU06,*) '*                                     *'
            WRITE (IU06,*) '***************************************'
            CALL ABORT1
          ENDIF
C
          MDEL=1
CCC          MDEL = NFRE 
CCC          DO WHILE (MINSIZE*MDEL.GE.MBXSIZE.AND.MDEL.GT.1)
CCC             MDEL=MDEL-1
CCC          ENDDO

          IF (ITEST.GE.2) THEN
            WRITE(IU06,*)''
            WRITE(IU06,*)' IN GETSPEC :'
            WRITE(IU06,*)'  MESSAGE MAIL BOX SIZE OF ',MBXSIZE,' BYTES'
            IF(MDEL.EQ.1) THEN
              WRITE(IU06,*)'  ',MDEL,
     .                     ' FREQUENCY IS DISTRIBUTED AT ONE TIME'
            ELSE
              WRITE(IU06,*)'  ',MDEL,
     .                     ' FREQUENCIES ARE DISTRIBUTED EACH TIME'
            ENDIF
            WRITE(IU06,*)'  WHICH SHOULD RESULT IN A MAXIMUM USAGE OF'
            WRITE(IU06,*)'  ',MINSIZE*MDEL,
     .                   ' BYTES FOR EACH MPDISTRIBFL CALL'
            WRITE(IU06,*)''
          ENDIF
        ENDIF
        ISEND=IREAD
C
        IF(LGRIBIN.AND..NOT.LRESTARTED) THEN
C
C       INPUT SPECTRA ARE IN GRIB FORMAT
C       --------------------------------
C
          IF(IRGG.EQ.1 .OR. CLDOMAIN == 'm') THEN
            KDEXN = 0
            KDEXS = 0
          ELSE
            KDEXN = (NINT((90. - AMONOP ) / XDELLA))*NGX
            KDEXS = (NINT((AMOSOP + 90. ) / XDELLA))*NGX
          ENDIF
C
          IF (IRGG.EQ.1) THEN
            KLENP=0
            DO K=1,NGY
              JSN=NGY-K+1
              KLENP = NLONRGG(JSN) + KLENP
            ENDDO
          ELSE
            KLENP = NGY * NGX + KDEXN + KDEXS
          ENDIF
C
          ALLOCATE(KSEC0(2),KSEC1(49+NANG+NFRE))
          ALLOCATE(KSEC2(384),KSEC3(2),KSEC4(128))
          ALLOCATE(PSEC2(96),PSEC3(2))
          ALLOCATE(IPOINTER(1:NEND(NPROC)))
          ALLOCATE(WORK(NIBLO))
C
          IF (IRGG.EQ.1) THEN
            DO IJ = 1,NEND(NPROC) 
              IX = IXLG(IJ,1)
              IY = NGY- KXLT(IJ,1) +1
              IPOINTER(IJ)=0
              DO IC=1,IY-1
                JSN=NGY-IC+1
                IPOINTER(IJ) = NLONRGG(JSN) + IPOINTER(IJ)
              ENDDO
              IPOINTER(IJ) = IPOINTER(IJ) + IX 
            ENDDO
          ELSE
            DO IJ = 1,NEND(NPROC)
              IX = IXLG(IJ,1)
              IY = NGY- KXLT(IJ,1) +1
              IPOINTER(IJ) = (IY-1)*NGX + IX + KDEXN
            ENDDO
          ENDIF
C
          IF (IRANK.EQ.IREAD) THEN
              IF (.NOT.LFDBIOIN) THEN
C               CONNECT INPUT PE (IREAD) WITH INPUT FILE
C
                CALL GRSTNAME(IU06,CDTPRO,CDATEF,0,'SGS',IRANK,NPROC,
     &                        FILENAME)
                LFILE=IECF_LEN(FILENAME)
                INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)
                IF(LLEXIST) THEN
C
                  CALL PBOPEN(IUGRS,FILENAME(1:LFILE),'r',KRET)
                  IF(KRET.LT.0) THEN 
                    WRITE (IU06,*) '***********************************'
                    WRITE (IU06,*) '*                                 *'
                    WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBOPEN  *'
                    WRITE (IU06,*) '* IN GETSPEC                      *'
                    IF(KRET.EQ.-1)
     &                WRITE (IU06,*) ' COULD NOT OPEN FILE ',FILENAME
                    IF(KRET.EQ.-2)
     &                WRITE (IU06,*) ' INVALID FILENAME ',FILENAME
                    IF(KRET.EQ.-3)
     &                WRITE (IU06,*)' INVALID OPEN MODE SPECIFIED'
                    WRITE (IU06,*) '*                                 *'
                    WRITE (IU06,*) '***********************************'
                    CALL ABORT1
                  ENDIF
                ELSE
                  WRITE(IU06,*)'**************************************'
                  WRITE(IU06,*)'*                                    *'
                  WRITE(IU06,*)'*GETSPEC : GRIB SPECTRA NOT FOUND IN *'
                  WRITE(IU06,*)  FILENAME
                  WRITE(IU06,*)'*PROGRAM WILL ABORT                  *'
                  WRITE(IU06,*)'*                                    *'
                  WRITE(IU06,*)'**************************************'
                  CALL ABORT1
                ENDIF
              ENDIF
          ENDIF
C
C
C         GET GRIB DATA FROM SPLIT (NFRE*NANG) FIELDS
C
C
C         READ NPROC-1 FIELDS AND SEND THEM SUCCESIVELY TO ALL OTHER PE'S
C         FOR DECODING
C
          ISTEP=MAX(NPROC-1,1)
          DO IC=1,NFRE*NANG,ISTEP 
            IST=ISTEP
            IF(IC+ISTEP.GT.NFRE*NANG ) IST=NFRE*NANG-IC+1
C
            DO IDUM=1,IST
              IF(NPROC.EQ.1) THEN
                KSEND=1
              ELSE IF(IDUM.LT.IREAD) THEN
                KSEND=IDUM
              ELSE
                KSEND=IDUM+1
              ENDIF
              M=(((IC-1)+IDUM-1)/NANG)+1
              K=(IC-1)+IDUM-(M-1)*NANG
C
C           DATA ARE READ IN ON PE IREAD
C           ----------------------------
              IF (LFDBIOIN.AND.IRANK.EQ.IREAD) THEN
C
C             READ GRIB DATA FROM FDB (IF REQUESTED)
C
                IF(M.EQ.1.AND.K.EQ.1) THEN
C
C                 INITIALIZE AND OPEN FDB.
C
                  istat = iopenfdb ('fdb', ifdb, 'r')
                  IF (ITEST.GE.2)WRITE(iu06,'("  FDB OPENED status=",I3,
     .            " Reference=", I2)') istat, ifdb
C
                  ICFDB2DSP = LEN_TRIM(CFDB2DSP)
                  IF (ICFDB2DSP .GT. 0 ) THEN
                    ISTAT = ISET_FDB_ROOT(IFDB, CFDB2DSP(1:ICFDB2DSP) )
                    IF ( ISTAT .NE. 0 ) THEN
                      WRITE(IU06,'("Error\ /Getspec/ Root dir specified:
     .                       ",A)') CFDB2DSP(1:ICFDB2DSP)
                      WRITE(IU06,'("Error\ /Getspec/
     .                        iopenfdb return status:",I3)') istat
                      kerr = 1
                      CALL ABORT1
                    ELSEIF (ITEST.GE.2) THEN
                      WRITE(IU06,'("\ /Getspec/ iset_fdb_root status:",
     .                      i3,
     .                " Root: ", a25)') istat, CFDB2DSP(1:ICFDB2DSP)
                    ENDIF
                  ENDIF
C
C                 SET FIXED PARAMETERS FOR FDB.
C
                  istat = isetvalfdb (ifdb, 'levty',  's'      )
                  istat = isetvalfdb (ifdb, 'param',  '251'    )
                  istat = isetvalfdb (ifdb, 'repres', 'll'     )
                  istat = isetvalfdb (ifdb, 'stream', 'wave'   )
                  istat = isetvalfdb (ifdb, 'step',   '000000' )
C---<             istat = isetvalfdb (ifdb, 'fstep',  '000000' )
C
                  kfcst = 0
                  IF ( kfcst .EQ. 0 ) THEN
                    cltype='an'
                  ELSE
                    cltype='fc'
                  ENDIF
C
                  clclass = yclass
                  clexpver = yexpver
C
                  CALL u2l1cr( cltype )
                  CALL u2l1cr( clclass )
                  CALL u2l1cr( clexpver )
C

                  istat = isetvalfdb (ifdb, 'expver', clexpver    )
                  istat = isetvalfdb (ifdb, 'date',   '19'//cdtpro(1:6))
                  istat = isetvalfdb (ifdb, 'type',   cltype     )
                  istat = isetvalfdb (ifdb, 'time',   cdtpro(7:8))
                  istat = isetvalfdb (ifdb, 'domain', cldomain )
C---<             istat = isetrankfdb (ifdb, ipe)
C
                  IF (ITEST.GE.2) THEN
                    WRITE(IU06,'("\ /getspec/ PARAM: ", a8)') '251'
                    WRITE(IU06,'("\ /getspec/ STEP:  ", a8)') '000000'
                    WRITE(IU06,'("\ /getspec/ CLASS: ", a8)') clclass
                    WRITE(IU06,'("\ /getspec/ EXPVER:", a8)') clexpver
                    WRITE(IU06,'("\ /getspec/ DOMAIN:", a8)') cldomain
                    WRITE(IU06,'("\ /getspec/ TY:   ", a8)') cltype
                    WRITE(IU06,'("\ /getspec/ DATE: ", a8)') cdtpro(1:6)
                    WRITE(IU06,'("\ /getspec/ TIME: ", a8)') cdtpro(7:8)
                  ENDIF
                ENDIF
C
C                 SET VARIABLE PARAMETERS FOR FDB.
C
                WRITE( cldir, '(i2.2)' ) K 
                istat = isetvalfdb( ifdb, 'direction', cldir )
                WRITE( clfre, '(i2.2)' ) M
                istat = isetvalfdb( ifdb, 'frequency', clfre )
                istat = isetvalfdb (ifdb, 'class',  clclass       )
                IF (ITEST.GE.2) THEN
                  WRITE(IU06,'("\ /getspec/ DIR:   ", a8)') cldir
                  WRITE(IU06,'("\ /getspec/ FREQ:  ", a8)') clfre
                ENDIF
C
C* try to estimate the size of the next product:
C
                icount = 1129120
                ALLOCATE (ingrib(icount))
                IF (ITEST.GE.2) WRITE(IU06,'(" GETSPEC: allocate", I10,
     .            " words for ingrib", " nstart=", i8, " nend=", i8,
     .            " spectrum=" , i4) ')icount,nstart(irank),nend(irank),
     .            nang*nfre
                irstat = ireadfdb (ifdb, ingrib, icount      )
                IF (irstat.NE.0 .OR. icount.LE.0) THEN
                  WRITE (IU06,'("***********************************")')
                  WRITE (IU06,'("*                                 *")')
                  WRITE (IU06,'("*ERROR FOLLOWING CALL TO ireadfdb *")')
                  WRITE (IU06,'("*IN GETSPEC                       *")')
                  WRITE (IU06,'("*STAT CODE irstat =", I10)') irstat
                  WRITE (IU06,'("*LENGTH OF MESSAGE=", I10)') icount
                  WRITE (IU06,'("*                                 *")')
                  WRITE (IU06,'("***********************************")')
                  CALL ABORT1
                ELSE IF (ITEST.GE.2) THEN
                  WRITE(IU06,'(" GETSPEC: readfdb read ",i10,
     .            " words (status:", i3, ")" )')icount, irstat
                ENDIF
C
                IF(M.EQ.NFRE.AND.K.EQ.NANG) THEN 
                  istat = iclosefdb(ifdb)
                  IF (ISTAT .NE. 0) CALL ABORT1
                ENDIF
C
                isize = icount
                IF (ITEST.GE.2) CALL FLUSH (IU06)
C
C             READ GRIB DATA FROM FILE (IF REQUESTED)
C
              ELSE IF (.NOT.LFDBIOIN.AND.IRANK.EQ.IREAD) THEN
                CALL PBSIZE (IUGRS, IPLENG)
                ISIZE=(IPLENG+NPRECI-1)/NPRECI
                ALLOCATE(INGRIB(ISIZE))
C
                CALL PBGRIB(IUGRS,INGRIB,IPLENG,ILENG,KRET)
                IF (KRET.EQ.-1) THEN
                  WRITE (IU06,*) ' SUB. GETSPEC : '
                  WRITE (IU06,*) ' REACHED EOF ON UNIT ',IUGRS
                  CALL ABORT1
                ELSEIF (KRET .EQ. -2) THEN
                  WRITE (IU06,*) ' SUB. GETSPEC : '
                  WRITE (IU06,*) ' FILE HANDLING ERROR ON UNIT ',IUGRS
                  CALL ABORT1
                ELSEIF (KRET .EQ. -3) THEN
                  WRITE (IU06,*) ' SUB. GETSPEC : '
                  WRITE (IU06,*) ' INGRIB IS NOT BIG ENOUGH.'
                  CALL ABORT1
                ENDIF
C
              ENDIF
C
              IF(IRANK.EQ.IREAD.AND.NPROC.NE.1) THEN
C 
C               SEND GRIB DATA TO PE KSEND
C
                ITAG=(M-1)*NANG+K
                CALL MPE_SEND(ISIZE,1,1,KSEND,ITAG,0,0,0,IERR)
                IF(IERR.LT.0) CALL MPEI_ABORT
     1                      ('MPE_SEND ERROR AT 1 in GETSPEC' )

                ITAG=((M-1)*NANG+K)*2
                CALL MPE_SEND(INGRIB,ISIZE,1,KSEND,ITAG,0,0,0,IERR)
                IF(IERR.LT.0) CALL MPEI_ABORT
     1                      ('MPE_SEND ERROR AT 2 in GETSPEC' )
C
                IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)
              ENDIF
C
              IF (IRANK.EQ.KSEND.AND.NPROC.NE.1) THEN
C
C               RECEIVED DATA FROM PE IREAD
C
                ITAG=(M-1)*NANG+K
                CALL MPE_RECV(ISIZE,1,1,IREAD,ITAG,0,0,0,
     $                        KRCOUNT,KRFROM,KRTAG,IERR)
                IF(IERR.LT.0) CALL MPEI_ABORT
     1                           ('MPE_RECV ERROR AT 1 in GETSPEC' )
                IF(KRTAG.NE.ITAG) CALL MPEI_ABORT
     1          ('MPE_RECV ERROR AT 1 in GETSPEC:  MISMATCHED TAGS' )
C
                ALLOCATE(INGRIB(ISIZE))
                ITAG=((M-1)*NANG+K)*2
                CALL MPE_RECV(INGRIB,ISIZE,1,IREAD,ITAG,0,0,0,
     $                        KRCOUNT,KRFROM,KRTAG,IERR)
                IF(IERR.LT.0) CALL MPEI_ABORT
     1                       ('MPE_RECV ERROR AT 2 in GETSPEC ' )
                IF(KRCOUNT.NE.ISIZE) CALL MPEI_ABORT
     1          ('MPE_RECV ERROR in 2 in GETSPEC:MISMATCHED MSG LENGTH')
                IF(KRTAG.NE.ITAG) CALL MPEI_ABORT
     1          ('MPE_RECV ERROR in 2 in GETSPEC:MISMATCHED TAGS' )
C
              ENDIF


              IF (IRANK.EQ.KSEND.OR.NPROC.EQ.1) THEN
C
C               DECODE THE GRIB DATA ON PE KSEND 
C
                KSEC3(2)=0
                PSEC3(2)=ZMISS
                ALLOCATE(PSEC4(KLENP))
C
                CALL GRIBEX (KSEC0,KSEC1,KSEC2,PSEC2,KSEC3,PSEC3,KSEC4,
     &                       PSEC4, KLENP,INGRIB,ISIZE,KWORD,'D',KRET)
                IF(KRET.GT.0) THEN 
                  WRITE(IU06,*) '************************************'
                  WRITE(IU06,*) '* FATAL ERROR IN SUB. GETSPEC      *'
                  WRITE(IU06,*) '* GRIBEX ERROR - KRET = ',KRET
                  WRITE(IU06,*) '* WHILE DECODING INGRB             *'
                  WRITE(IU06,*) '*                                  *'
                  WRITE(IU06,*) '*    PROGRAM ABORTS                *'
                  WRITE(IU06,*) '*                                  *'
                  WRITE(IU06,*) '************************************'
                  CALL ABORT1
                ENDIF
C
                IF(K.NE.KSEC1(44)) THEN
                  WRITE(IU06,*) '************************************'
                  WRITE(IU06,*) '* FATAL ERROR IN SUB. GETSPEC      *'
                  WRITE(IU06,*) '* REQUESTED AND DECODED DIRECTIONAL*'
                  WRITE(IU06,*) '* INDEX ARE DIFFERENT :            *'
                  WRITE(IU06,*) '* REQUESTED : ',K 
                  WRITE(IU06,*) '* DECODED   : ',KSEC1(44) 
                  WRITE(IU06,*) '*                                  *'
                  WRITE(IU06,*) '************************************'
                  CALL ABORT1
                ENDIF
C
                ZTHETA = FLOAT(KSEC1(49+K))/KSEC1(48) 
                IF(ABS(ZTHETA-(TH(K)*DEG)).GT.1.) THEN
                  WRITE(IU06,*) '************************************'
                  WRITE(IU06,*) '* FATAL ERROR IN SUB. GETSPEC      *'
                  WRITE(IU06,*) '* REQUESTED AND DECODED DIRECTIONS *'
                  WRITE(IU06,*) '* ARE DIFFERENT :                  *'
                  WRITE(IU06,*) '* ZTHETA = ',ZTHETA
                  WRITE(IU06,*) '* TH(K)*DEG = ',TH(K)*DEG
                  WRITE(IU06,*) '*                                  *'
                  WRITE(IU06,*) '*    PROGRAM ABORTS                *'
                  WRITE(IU06,*) '*                                  *'
                  WRITE(IU06,*) '************************************'
                  CALL ABORT1
                ENDIF
C
                IF(M.NE.KSEC1(45)) THEN
                  WRITE(IU06,*) '************************************'
                  WRITE(IU06,*) '* FATAL ERROR IN SUB. GETSPEC      *'
                  WRITE(IU06,*) '* REQUESTED AND DECODED FREQUENCY  *'
                  WRITE(IU06,*) '* INDEX ARE DIFFERENT :            *'
                  WRITE(IU06,*) '* REQUESTED : ',M 
                  WRITE(IU06,*) '* DECODED   : ',KSEC1(45) 
                  WRITE(IU06,*) '*                                  *'
                  WRITE(IU06,*) '************************************'
                  CALL ABORT1
                ENDIF
C
                ZFRE = FLOAT(KSEC1(49+NANG+M))/KSEC1(49) 
                IF(ABS(ZFRE-FR(M)).GT.0.01*FR(M)) THEN
                  WRITE(IU06,*) '************************************'
                  WRITE(IU06,*) '* FATAL ERROR IN SUB. GETSPEC      *'
                  WRITE(IU06,*) '* REQUESTED AND DECODED FREQUENCIES*'
                  WRITE(IU06,*) '* ARE DIFFERENT :                  *'
                  WRITE(IU06,*) '* ZFRE = ',ZFRE
                  WRITE(IU06,*) '* FR(M) = ',FR(M)
                  WRITE(IU06,*) '*                                  *'
                  WRITE(IU06,*) '*    PROGRAM ABORTS                *'
                  WRITE(IU06,*) '*                                  *'
                  WRITE(IU06,*) '************************************'
                  CALL ABORT1
                ENDIF
C
                IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)
C
                DO IJ = 1,NEND(NPROC) 
                  WORK(IJ) = PSEC4(IPOINTER(IJ)) 
                ENDDO
                IF(ALLOCATED(PSEC4)) DEALLOCATE(PSEC4)
C
                PPREC=0.0
                XEPS=1.0e-10
                DO IJ = 1,NEND(NPROC) 
                  IF (WORK(IJ) .NE. ZMISS) THEN
                    WORK(IJ) = 10.**(WORK(IJ)-ABS(PPREC))- XEPS
                  ELSE
                    WORK(IJ) = 0.0
                  ENDIF
                ENDDO
C
              ENDIF
C
            END DO
C
C           DISTRIBUTE THE RESPECTIVE CONTRIBUTION OF WORK TO EACH PE.
C
            ALLOCATE(WORKRCV(NIBLO))
C
            DO IDUM=1,IST
              IF(NPROC.EQ.1) THEN
                KSEND=1
              ELSE IF(IDUM.LT.IREAD) THEN
                KSEND=IDUM
              ELSE
                KSEND=IDUM+1
              ENDIF
              M=(((IC-1)+IDUM-1)/NANG)+1
              K=(IC-1)+IDUM-(M-1)*NANG

              ITAG=((M-1)*NANG+K)*3
              IF(IRANK.EQ.KSEND) THEN
                CALL MPDISTRIBSCFLD(KSEND,ITAG,NSTART,NEND,WORK)
                DO IJ = NSTART(IRANK),NEND(IRANK) 
                  FL(IJ,K,M) = WORK(IJ)
                ENDDO
              ELSE
                CALL MPDISTRIBSCFLD(KSEND,ITAG,NSTART,NEND,WORKRCV)
                DO IJ = NSTART(IRANK),NEND(IRANK) 
                  FL(IJ,K,M) = WORKRCV(IJ)
                ENDDO
              ENDIF
C
              FL(NINF-1,K,M) = 0.
            END DO
C
            IF(ALLOCATED(WORKRCV)) DEALLOCATE(WORKRCV)
C
          END DO
C
          IF (.NOT.LFDBIOIN.AND.IRANK.EQ.IREAD) THEN
            CALL PBCLOSE(IUGRS,KRET)
            IF(KRET.LT.0) THEN 
              WRITE (IU06,*) '**********************************'
              WRITE (IU06,*) '*                                *'
              WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBCLOSE '
              WRITE (IU06,*) '* IN GETSPEC                     *'
              WRITE (IU06,*) '* FILE ',FILENAME
              WRITE (IU06,*) '*                                *'
              WRITE (IU06,*) '**********************************'
              CALL ABORT1
            ENDIF
          ENDIF
C
          IF(ALLOCATED(KSEC0)) DEALLOCATE(KSEC0)
          IF(ALLOCATED(KSEC1)) DEALLOCATE(KSEC1)
          IF(ALLOCATED(KSEC2)) DEALLOCATE(KSEC2)
          IF(ALLOCATED(KSEC3)) DEALLOCATE(KSEC3)
          IF(ALLOCATED(KSEC4)) DEALLOCATE(KSEC4)
          IF(ALLOCATED(PSEC2)) DEALLOCATE(PSEC2)
          IF(ALLOCATED(PSEC3)) DEALLOCATE(PSEC3)
          IF(ALLOCATED(IPOINTER)) DEALLOCATE(IPOINTER)
          IF(ALLOCATED(WORK)) DEALLOCATE(WORK)
C
C
        ELSE IF(LPBIOIN.AND.(.NOT.LGRIBIN.OR.LRESTARTED)) THEN
C
C       INPUT SPECTRA ARE IN A BINARY FILE
C
          IF (IRANK.EQ.ISEND)
     1        CALL GRSTNAME(IU06,CDTPRO,CDATEF,0,'BLS',IRANK,NPROC,
     &                      FILENAME)

C use the ability of the PBIO software to read input by frequency band

          DO MLOOP=1,NFRE,MDEL

             MINF=MLOOP
             MSUP=MIN(MLOOP+MDEL-1,NFRE)


             ALLOCATE(RFL(0:NIBLO,NANG,MINF:MSUP))

C            READ RESTART SPECTRA FROM PE ISEND (IREAD) 
             IF (IRANK.EQ.ISEND) THEN

               LPBOPEN = .FALSE.
               LPBCLOSE = .FALSE.
               IF(MINF.EQ.1) LPBOPEN = .TRUE.
               IF(MSUP.EQ.NFRE) LPBCLOSE = .TRUE.

               CALL READFL(RFL,MINF,MSUP,FILENAME,IUNIT,
     1                     LPBOPEN,LPBCLOSE)
             ENDIF

C            PUT A BARRIER TO AVOID FLOODING THE CROSSBAR 
             CALL MPE_BARRIER(IERR)
             IF(IERR.LT.0) 
     1       CALL MPEI_ABORT('MPE_BARRIER ERROR in GETSPEC')

             CALL MPDISTRIBFL(ISEND,KTAG,NSTART,NEND,MINF,MSUP,
     1                        MINF,MSUP,RFL)
             KTAG=KTAG+1
             IF (ITEST.GE.2)
     1       WRITE(IU06,*)
     2       'SUB. GETSPEC: RESTART SPECTRUM COLLECTED, PART:',MLOOP

C            KEEP CORRESPONDING CONTRIBUTION TO FL
             DO M=MINF,MSUP
               DO K=1,NANG
                  DO IJ=NSTART(IRANK),NEND(IRANK)
                    FL(IJ,K,M) = RFL(IJ,K,M)
                  END DO
               END DO
             END DO

             DO M=MINF,MSUP
               DO K=1,NANG
                  FL(NINF-1,K,M) = 0. 
               END DO
             END DO

             IF(ALLOCATED(RFL)) DEALLOCATE(RFL)

          END DO

        ELSE

           ALLOCATE(RFL(0:NIBLO,NANG,NFRE))

C          READ RESTART SPECTRA FROM PE ISEND (IREAD)
           IF (IRANK.EQ.ISEND) THEN
             CALL READFL(RFL,1,NFRE,FILENAME,IUNIT,LPBOPEN,LPBCLOSE)
           ENDIF

           DO MLOOP=1,NFRE,MDEL
              MINF=MLOOP
              MSUP=MIN(MLOOP+MDEL-1,NFRE)

C            PUT A BARRIER TO AVOID FLOODING THE CROSSBAR 
             CALL MPE_BARRIER(IERR)
             IF(IERR.LT.0) 
     1       CALL MPEI_ABORT('MPE_BARRIER ERROR in GETSPEC 2')

             CALL MPDISTRIBFL(ISEND,KTAG,NSTART,NEND,MINF,MSUP,
     1                         1,NFRE,RFL)
             KTAG=KTAG+1
           ENDDO

C            KEEP CORRESPONDING CONTRIBUTION TO FL
             DO M=1,NFRE
               DO K=1,NANG
                  DO IJ=NSTART(IRANK),NEND(IRANK)
                    FL(IJ,K,M) = RFL(IJ,K,M)
                  END DO
               END DO
             END DO

             DO M=1,NFRE
               DO K=1,NANG
                  FL(NINF-1,K,M) = 0. 
               END DO
             END DO

           IF (ITEST.GE.2)
     1     WRITE(IU06,*)
     2     'SUB. GETSPEC: RESTART SPECTRUM COLLECTED'


           IF(ALLOCATED(RFL)) DEALLOCATE(RFL)

        ENDIF
C
C
      ELSE

        IF(LGRIBIN) THEN
          WRITE(IU06,*) ' starting from grib spectra in the non message'
          WRITE(IU06,*) ' passing configuaration is not available yet'
          WRITE(IU06,*) ' program will abort'
          CALL ABORT1
        ENDIF
        IF(LPBIOIN) THEN
          CALL GRSTNAME(IU06,CDTPRO,CDATEF,0,'BLS',IRANK,NPROC,FILENAME)
          LPBOPEN = .TRUE.
          LPBCLOSE = .TRUE.
        ENDIF

        CALL READFL(FL,1,NFRE,FILENAME,IUNIT,LPBOPEN,LPBCLOSE)

      END IF
C
        WRITE(IU06,*) ' SPECTRUM FILE READ IN........',
     1                ' CDTPRO  = ', CDTPRO
C
      RETURN
      END
