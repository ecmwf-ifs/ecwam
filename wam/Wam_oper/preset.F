      PROGRAM preset 
C
C ----------------------------------------------------------------------
C
C**** *PRESET* - GENERATES ALL FILES REQUIRED FOR A WAMODEL START.
C
C     SUSANNE HASSELMANN  MPI     JULY 1986.
C
C     ANNEGRET SPEIDEL    MPI      MAY 1988 PARAMETER STATEMENTS.
C
C     ANNEGRET SPEIDEL    MPI NOVEMBER 1988 CRAY-2 VERSION.
C
C     CYCLE_3 MODICIFATIONS:
C     ----------------------
C
C     RENATE PORTZ       MPI      JUNE 1990 COMPUTATION OF INITIAL
C                                           JONSWAP SPECTRA FROM
C                                           INITIAL WIND FIELD.
C
C     CYCLE_4 MODICIFATIONS:
C     ----------------------
C
C     H. GUNTHER  GKSS/ECMWF  DECEMBER 1990
C
C     J. BIDLOT    ECMWF   FEBRUARY 1996  MESSAGE PASSING
C     J. BIDLOT    ECMWF   MARCH 1997  MODIFY ROUTINES FOR OUTPUT OF
C                          RESTART FILE 
C     B. HANSEN    ECMWF   APRIL    1997  RESTART FACILITY.
C     B. HANSEN    ECMWF   JANUARY  1998  NAMELIST INPUT.
C     B. HANSEN    ECMWF   FEBRUARY 1998  WRITE 2DSP TO FDB.
C     J. BIDLOT    ECMWF   MARCH 1998 PRODUCE SPLIT GRIB SPECTRA AND
C                                     INTERPOLATE TO OUTPUT GRID.
C                                     NO LAW FILE IS PRODUCED THEN.
C
C*    PURPOSE.
C     --------
C
C       TO INITIALISE ALL FILES REQUESTED BY THE WAMODEL.
C
C**   INTERFACE.
C     ----------
C
C       *IU01*   INTEGER    INPUT  UNIT UNBLOCKED WIND FILE.
C                           (SEE SUB READWND).
C       *IU05*   INTEGER    USER INPUT UNIT.
C       *IU06*   INTEGER    PRINTER OUTPUT.
C       *IU07*   INTEGER    INPUT  UNIT PREPROC GRID OUTPUT.
C       *IU12*   INTEGER    OUTPUT UNIT BLOCKS OF SPECTRA.
C       *IU14*   INTEGER    OUTPUT UNIT SECOND LAT OF BLOCKS.
C
C       *IU15*   INTEGER    OUTPUT UNIT LAST WINDFIELDS.
C                 !OR!      OUTPUT UNIT LAST DRAG COEFFICIENT.
C
C     METHOD.
C     -------
C
C       JONSWAP PARAMETERS ARE DEFINED EITHER BY USER INPUT OR
C       BY FETCH LAWS. THE 2-D SPECTRA ARE COMPUTED FOR EACH
C       POINT IN A BLOCK, THE WAMODEL BLOCKS ARE INITIALISED BY THESE
C       SPECTRA AND ALL BLOCKS AND OVERLAPPING LATITUDES ARE SAVED.
C       IF FETCH LAWS ARE USED TO DEFINE PARAMETERS THE FIRST WIND
C       FIELD IS GENERATED OTHERWISE THE WIND FIELD IS INITIALISED
C       WITH ZEROS. THE MAIN MODEL WILL RECONSTRUCT THE WIND ANY HOW.
C
C       THE FILE HANDLING OF THE RESTART FILES IS COMPUTER DEPENDENT.
C       SUB GSFILE HAS TO BE MODIFIED, TO COPY THE UNIT ALIAS FILES
C       (UNITS IU12,IU14, AND IU15) TO PERMANENT FILES.
C
C     EXTERNALS.
C     ----------
C
C       *ABORT1*     - TERMINATES PROCESSING.
C       *CREWFN*    - CREATES A WIND FILE NAME.
C       *GRSTNAME*  - BUILD NAME FOR INPUT/OUTPUT FILES 
C       *GSFILE*    - GETS OR SAVES FILES (COMPUTER DEPENDENT).
C       *iecf_len*  - LENGTH OF A CHARACTER ARRAY.
C       *INCDATE*   - UPDATES A DATE TINE GROUP.
C       *INTGRS*    - PERFORMS INTERPOLATIONS OF INPUT SPECTRA 
C       *LOCINT*    - INTERPOLATES TO MODEL GRID.
C       *MSTART*    - GENERATES THE RESTART FILES.
C       *NOTIM*     - CONTROLS WIND GENERATION (NO TIME INTERPOLATION).
C       *OPENFIL*   - OPENS A FILE.
C       *PEAK*      - COMPUTES PARAMETERS BY FETCH LAWS.
C       *PREWIND*   - PREPARES WINDS.
C       *SAVSPEC*   - SAVES THE RESTART SPECTRA.
C       *SAVSTRESS* - SAVES THE RESTART WIND AND STRESS FIELDS.
C       *SPECTRA*   - COMPUTES SPECTRA FROM PARAMETERS.
C       *SPLITBL*   - SAVES SECOND TO LAST LATITUDES.
C       *SPLITGRS*  - PERFORMS SPLITTING OF SPECTRA ON SUB DOMAINS 
C       *SPR*       - DIRECTIONAL DISTRIBUTION.
C       *READGRS*   - READS GRIBBED SPECTRUM DATA
C       *READPRE*   - READS PREPROC OUTPUT.
C       *READWND*   - READS A WIND FIELD.
C       *TIMIN*     - CONTROLS WIND GENERATION (TIME INTERPOLATION).
C
C     REFERENCES
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------

      USE WAVE_MODULE

C -------------------------------------------------------------------

C
#include "param.h"
C
#include "txtpara.h"
C
#include "parameter_wind.h"
C
#include "txtpara_wind.h"
C
#include "comcoup.h"
C
#include "txtcoup.h"
C
#include "comcout.h"
C
#include "txtcout.h"
C
#include "comcurr.h"
C
#include "txtcurr.h"
C
#include "comfred.h"
C
#include "txtfred.h"
C
#include "comgrid.h"
C
#include "txtgrid.h"
C
#include "comindn.h"
C
#include "txtindn.h"
C
#include "comjons.h"
C
#include "txtjons.h"
C
#include "commap.h"
C
#include "txtmap.h"
C
#include "comshal.h"
C
#include "txtshal.h"
C
#include "comstat.h"
C
#include "txtstat.h"
C
#include "comtabl.h"
C
#include "txttabl.h"
C
#include "comtest.h"
C
#include "txttest.h"
C
#include "comtext.h"
C
#include "txttext.h"
C
#include "comunit.h"
C
#include "txtunit.h"
C
C#include "comwami.h"
C
#include "comwind.h"
C
#include "txtwind.h"
C
#include "comwndg.h"
C
#include "txtwndg.h"
C
#include "parcons.h"
C
C
C
C THIS COMMON BLOCK IS NECESSARY SINCE SOME ROUTINE ARE USED IN WAMODEL
C WHERE MESSAGE PASSING CAN BE ACTIVITED. IN PRESET THE PARAMETERS ARE
C SET TO THEIR DEFAULT VALUES FOR RUN ON ONE PROCESSOR.
C
#include "commesspass.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
C* 0.    COMMON OUTPUT CONTROL.
C  ----------------------------
#include <outctrl.h>
      DATA csubna /"PRESET "/
C
C ----------------------------------------------------------------------
      REAL X4(2)
      INTEGER I4(2)
C
      CHARACTER*  1 CLTUNIT, C_DOMAIN, CMODE
      CHARACTER* 10 ZERO, CDUM, FILENAME
      CHARACTER* 70 HEADER
      CHARACTER* 72 LINE
      CHARACTER* 80 SFILENAME, ISFILENAME
C
      PARAMETER (ZERO=' ', CDUM='0000000000')
C
      LOGICAL LDSTOP, LDWRRE, LWAVEWIND, LINTGRS
      LOGICAL LSPLITIN, LSPLITOUT
      DATA LDSTOP, LDWRRE / 2*.FALSE. /
C
C DIMENSION DUMMY COUPLED VARIABLES
C
      INTEGER IUGB(1),IVGB(1)
C
C
C ----------------------------------------------------------------------
C
      NAMELIST /NALINE/ HEADER,
     1          IOPTI, ITEST, ITESTB,
     1          ALFA, FM, GAMMA, SA, SB, THETA,
     1          FETCH,
     1          USERID, RUNID, PATH, CPATH,
     1          CDATEA, IDELWI, CLTUNIT,
     1          LPBIOIN, LPBIOOUT, LFDBIOIN, LFDBIOOUT,
     1          LSPLITIN,LSPLITOUT,
     1          MARSTYPE, YCLASS, YEXPVER, CFDB2DSP, C_DOMAIN
C
C ----------------------------------------------------------------------
C
C
C*    0. SET DEFAULT VALUES FOR THE NAMELIST ELEMENTS.
C        ---------------------------------------------
C
      msub = 0
      csuvi = ' '
      ctest = ' '
      mtest = 0
      msuvi = 0
C
      HEADER = ZERO
      IOPTI  =   -9
      ITEST  =   -9
      ITESTB =   -9
      ALFA   =    0.0
      FM     =    0.0
      GAMMA  =    0.0
      SA     =    0.0
      SB     =    0.0
      THETA  =    0.0
      FETCH  =    0.0
      USERID = ZERO
      RUNID  = ZERO
      PATH   = ZERO
      CPATH  = ZERO
      CDATEA = ZERO
      CLTUNIT= ZERO
      IDELWI =    0
C
      LMESSPASS = .FALSE.
C
      LPBIOIN   = .FALSE.
      LPBIOOUT  = .FALSE.
      LGRIBIN   = .FALSE.
      LFDBIOIN  = .FALSE. 
      LFDBIOOUT = .FALSE. 
      LSPLITIN  = .FALSE.
      LSPLITOUT = .FALSE.
C
      MARSTYPE = 'an'
      YCLASS   = 'rd'
      YEXPVER  = USERID//'a'
      C_DOMAIN = 'g'
C
C*    1. DEFINE UNIT NAMES.
C        ------------------
C
 1000 CONTINUE
      IU04 = 4
      IU05 = 5
      IU06 = 6
      IU07 = 7
      IU12 = 12
      IU14 = 14
      IU15 = 15
C
C     1.1 INITIALISATION OF MPP VARIABLES WITH DEFAULT VALUES
C         ---------------------------------------------------
C
       IRANK=1
       NPROC=1
       NINF=1
       NSUP=NIBLO
       X4=1.
       NPRECR = N_PRECISION(X4)
       I4=1
       NPRECI = N_PRECISION(I4)
       KTAG=100
C
C
C ----------------------------------------------------------------------
C
C ALLOCATE NECESSARY ARRAYS
C
      ALLOCATE (NSTART(NPROC),NEND(NPROC))
      NSTART=1
      NEND=NIBLO
C
C ----------------------------------------------------------------------
C
C*    2. READ NAMELIST NALINE.
C        ---------------------
C
      CDATEWO = ' '
      CDAWIFL = ' '
      CDATEFL = ' '

      READ (IU05, NALINE)
C
      IF (CLTUNIT .EQ. 'H') IDELWI = IDELWI*3600
      CDATEF = CDATEA
      CDATER = '0000000000'
      ICPLEN=LEN_TRIM(CPATH)
C
C ----------------------------------------------------------------------
C
C*    3. READ PREPROC OUTPUT.
C        --------------------
C
 3000 CONTINUE
C
      CALL READPRE (IU07, 0)
C
      NEND=IJL(1)
C
C ----------------------------------------------------------------------
C
C*    4. PRINTER PROTOCOL OF INPUT.
C        --------------------------
C
 4000 CONTINUE
      WRITE (IU06,'(" PROGRAM PRESET:", / , "
     . WAM Cycle4 Release8 CI 4 (common with IFS CY18R6)
     .    CC=wam_sth_CY18R5_oam_for_18r6", /,
     . "  USER INPUT:" )')

      WRITE (IU06,'(A70)') HEADER
      WRITE (IU06,'('' MODEL OPTIONS  :'',/)')
      IF (IOPTI.EQ.0) THEN
         WRITE (IU06,'('' INITIAL VALUES ARE COMPUTED FROM'',
     1                 '' INPUT PARAMETERS.'')')
      ELSE IF (IOPTI.EQ.1) THEN
         WRITE (IU06,'('' INITIAL VALUES ARE COMPUTED FROM'',
     1                 '' LOCAL WIND.'')')
         WRITE (IU06,'('' WAVE ENERGY IS ZERO IN CALM WIND AREAS.'')')
      ELSE IF (IOPTI.EQ.2) THEN
         WRITE (IU06,'('' INITIAL VALUES ARE COMPUTED FROM'',
     1                 '' LOCAL WIND.'')')
         WRITE (IU06,'('' PARAMETERS USED IN CALM WIND AREAS.'')')
      ELSE IF (IOPTI.EQ.-1.OR.IOPTI.EQ.-2) THEN
         WRITE (IU06,'('' INITIAL SPECTRAL VALUES ARE COMPUTED FROM'',
     1                 '' ARCHIVED SPECTRUM AND DRAG COEF. FIELDS.'')')
         WRITE (IU06,'('' AND ARCHIVED WAVE MODEL ANALYSED WIND SPEED'',
     1                 '' IF AVAILABLE.'')')
      END IF
C
      WRITE(IU06,*) '  '
      IF(LPBIOIN) THEN
         WRITE (IU06,*)
     1   ' PBIO SOFTWARE IS USED TO READ SELECTED INPUT FILES' 
      ENDIF
      IF(LPBIOOUT) THEN
         WRITE (IU06,*)
     1   ' PBIO SOFTWARE IS USED TO WRITE TO SELECTED OUTPUT FILES' 
      ENDIF
      WRITE(IU06,*) '  '
C
      IF(LFDBIOIN) THEN
         WRITE (IU06,*)
     1   ' FDB SOFTWARE IS USED TO READ SELECTED INPUT FILES' 
      ENDIF
      IF(LFDBIOOUT) THEN
         WRITE (IU06,*)
     1   ' FDB SOFTWARE IS USED TO WRITE TO SELECTED OUTPUT FILES' 
         WRITE(IU06,'("   ITS ROOT DIRECTORY IS ", A)')
     1        CFDB2DSP(1:LEN_TRIM(CFDB2DSP))
      ENDIF
C
      IF (IOPTI.EQ.-1.OR.IOPTI.EQ.-2) THEN
        IF(LSPLITIN) THEN
          WRITE (IU06,*)'  '
          WRITE (IU06,*)' ******************************************'
          WRITE (IU06,*)' THE INPUT SPECTRA ARE GRIB SPLIT FILES '
          WRITE (IU06,*)' ******************************************'
          WRITE (IU06,*)'  '
        ELSE
          WRITE (IU06,*)'  '
          WRITE (IU06,*)' ******************************************'
          WRITE (IU06,*)' INPUT SPECTRA ARE OLDFASHION GRIB FILES '
          WRITE (IU06,*)' ******************************************'
          WRITE (IU06,*)'  '
        ENDIF
        WRITE (IU06,*)'  '
        WRITE (IU06,*)' ******************************************'
        IF(LSPLITOUT) THEN
          WRITE (IU06,*)' THE OUTPUT OF SPECTRA WILL BE GRIBBED '
          WRITE (IU06,*)' AND SPLIT INTO ',NFRE*NANG,' FILES'
        ELSE
          WRITE (IU06,*)' THE OUTPUT OF SPECTRA WILL BE BINARY  '
        ENDIF
        WRITE (IU06,*)' ******************************************'
        WRITE (IU06,*)'  '
      ENDIF
      WRITE (IU06,*) ' TEST OUTPUT LEVEL IS .......... ITEST = ', ITEST
      WRITE (IU06,*) ' TEST OUTPUT IN BLOCK LOOP UPTO ITESTB = ', ITESTB
C
      IF (IOPTI.NE.-1.AND.IOPTI.NE.-2) THEN
        WRITE (IU06,'('' JONSWAP PARAMETERS  :'',/)')
        WRITE (IU06,'('' ALFA : '',F10.5,'' FM : '',F10.5,'' GAMMA : '',
     1                F10.5,'' SA : '',F10.5,'' SB : '',F10.5)')
     2                ALFA, FM, GAMMA, SA, SB
        WRITE (IU06,'('' MEAN WAVE DIRECTION :  THETA = '',F10.5,
     1                '' DEGREE'')')  THETA
      ENDIF
      WRITE (IU06,*) '  '
      WRITE (IU06,*) ' MODEL FILE NAMES:'
      WRITE (IU06,*) ' USERID IS ........: ', USERID
      WRITE (IU06,*) ' RUN IDENTIFIER IS : ', RUNID
      WRITE (IU06,*) ' PATH NAME IS .....: ', PATH
      WRITE (IU06,*) '  '
      WRITE (IU06,*) ' WIND INPUT TIMESTEP (SECONDS)      : ',IDELWI
      WRITE (IU06,*) '  '
      WRITE (IU06,*) ' END OF USER INPUT PROTOCOLL'
      WRITE (IU06,'(''  NUMBER OF DIRECTION BINS  NANG = '',I4)') NANG
      WRITE (IU06,'(''  NUMBER OF FREQUENCY BINS  NFRE = '',I4)') NFRE
C
      IF (LFDBIOOUT .OR. LFDBIOIN) THEN
        WRITE(IU06,*) ' ACCESS TO FIELDS DATA BASE REQUESTED '
        WRITE(IU06,*) ' DATA BASE WILL BE INITIALIZED '
        INITFDB_STAT = IINITFDB()
        IF (INITFDB_STAT .NE. 0) THEN
           WRITE(IU06,'("****************************************")')
           WRITE(IU06,'("*                                      *")')
           WRITE(IU06,'("*    FATAL ERROR IN SUB. USERIN        *")')
           WRITE(IU06,'("*    ==========================        *")')
           WRITE(IU06,'("*    INITIALIZATIO OF THE              *")')
           WRITE(IU06,'("*    FIELDS DATA BASE                  *")')
           WRITE(IU06,'("*     F A I L E D .                    *")')
           WRITE(IU06,'("*                                      *")')
           WRITE(IU06,'("*  ERROR CODE (ISTAT) =",i7, "         *")')
     .     INITFDB_STAT
           WRITE(IU06,'("*                                      *")')
           WRITE(IU06,'("****************************************")')
           CALL ABORT1
        ENDIF
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    5. PREPARE WINDFIELD.
C        ------------------
C
 5000 CONTINUE
      CDA = ZERO
      CDTPRO  = CDATEA
      CDATEE  = CDATEA
      CDAWIFL = CDATEA
      IDELPRO = IDELWI
      IDELWO  = IDELWI
C
      LWCOU=.FALSE.
      ILEN=1
C
      IF (IOPTI.NE.-1.AND.IOPTI.NE.-2) THEN
        IF(.NOT.ALLOCATED(FL1)) ALLOCATE (FL1(0:NIBLO,NANG,NFRE))
        IF(.NOT.ALLOCATED(U10OLD)) ALLOCATE(U10OLD(NIBLO,NBLO))
        IF(.NOT.ALLOCATED(THWOLD)) ALLOCATE(THWOLD(NIBLO,NBLO))
        IF(.NOT.ALLOCATED(USOLD)) ALLOCATE(USOLD(NIBLO,NBLO))
        IF(.NOT.ALLOCATED(Z0OLD)) ALLOCATE(Z0OLD(NIBLO,NBLO))
        IF(.NOT.ALLOCATED(TAUW)) ALLOCATE(TAUW(NIBLO,NBLO))
        U10OLD = 0.
        THWOLD = 0.
        USOLD = 0.
        TAUW = 0.
        Z0OLD = 0.
      ENDIF

      IF (IOPTI.GT.0) THEN
         CALL PREWIND (NC,NR,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     1              IUGB, IVGB, ILEN)
         IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. PREWIND DONE'
C
      ELSE IF (IOPTI.EQ.-1.OR.IOPTI.EQ.-2) THEN
C
C       GET SPECTRA FROM GRIBBED DATA.
C
        IF(IOPTI.EQ.-1) THEN
          ISPLIT=1
        ELSE
          ISPLIT=2
        ENDIF
C
        CALL GRSTNAME(IU06,CDTPRO,CDATEF,0,'GRS',IRANK,NPROC,
     &                SFILENAME)
        CALL GRSTNAME(IU06,CDTPRO,CDATEF,0,'IGS',IRANK,NPROC,
     &                ISFILENAME)
C
C       INTERPOLATE THE GLOBAL FIELD IF IT IS NECESSARY
C
        LINTGRS=.FALSE. 
        CALL INTGRS (IU06,SFILENAME,ISPLIT,CDTPRO,LINTGRS,ISFILENAME)

        IF(LINTGRS) SFILENAME=ISFILENAME 

        IF(.NOT.LSPLITOUT) THEN
C         THE SPECTRA WILL BE SAVED AS A BINARY FILE
CCCC      note !! FL1 is allocated in readgrs !!!!
          CALL READGRS(IU06,SFILENAME,ISPLIT,CDTPRO)
          IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. READGRS DONE'
        ELSE
C
C         SPLIT SPECTRA INTO NFRE*NANG GRIB FILES
C
          CALL SPLITGRS (IU06,SFILENAME,ISPLIT,CDTPRO,ITEST,LFDBIOOUT,
     &                   MARSTYPE, YCLASS, YEXPVER, CFDB2DSP,C_DOMAIN)
          IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. SPLITGRS DONE'

        ENDIF
      ELSE
         IF (ITEST.GT.0) WRITE (IU06,*) ' WIND SET TO ZERO'
      END IF
C
C ----------------------------------------------------------------------
C
C*    6. DEFINE FETCH AND MAXIMUM PEAK FREQUENCY.
C        ----------------------------------------
C
 6000 CONTINUE
C
      IF (IOPTI.NE.-1.AND.IOPTI.NE.-2) THEN
        IF (FETCH.LT.0.1E-5) FETCH = 0.5*DELPHI
        FRMAX = FM
        IF (IOPTI.NE.0) THEN
           WRITE (IU06,*) ' FETCH USED (METRES)       : ', FETCH
           WRITE (IU06,*) ' MAXIMUM PEAK FREQUENCY IS : ', FRMAX
        ENDIF
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    7. GENERATE AND WRITE START FILES.
C        -------------------------------
C
 7000 CONTINUE
      IF (IOPTI.NE.-1.AND.IOPTI.NE.-2) THEN
        THETAQ = THETA * RAD
        CALL MSTART (IU12, IU14, IU15, IOPTI, FETCH, FRMAX,
     1              FL1,U10OLD,THWOLD,USOLD,Z0OLD,TAUW)
        IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. MSTART DONE'
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    8. DISPOSE START FILES.
C        --------------------
C
 8000 CONTINUE
C
      IF (IOPTI .NE. -1 .AND. IOPTI .NE. -2  ) THEN
        CALL SAVSTRESS(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,NSTART,NEND,
     1                 LDSTOP,LDWRRE,1,CDATEA,CDATEA)
        IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. SAVSTRESS DONE'
        IF (LSPLITOUT.AND..NOT.LPBIOOUT)
     &     CALL GSFILE (IU06, IU15, 0, CDATEA, CDATEA, 'LAW', 'S')
      ENDIF


      IF (IOPTI .NE. -1 .AND. IOPTI .NE. -2 .AND.
     .    LSPLITOUT) THEN

C       THE SPECTRA WILL BE SAVED AS GRIB FILES IN THE FDB.

        NTOTENS = 0
        CALL SPEC2FDB (IU06, NPROC, IRANK, LWCOU, TH, FR,
     1                 FL1, NINF, NSUP, NANG, NFRE, NGX, NGY,IRGG,
     2                 NLONRGG, KXLT, IXLG, NIBLO, NBLO, ZMISS,
     3                 AMONOP, AMOWEP, AMOEAP, AMOSOP, XDELLA,
     4                 XDELLO, NSTART,NEND,
     5                 CDATEA, CDATEF, CDATER, CDATES, CDATEA,
     6                 MARSTYPE, YCLASS, YEXPVER, CFDB2DSP,
     7                 NENSFNB, NTOTENS, LFDBIOOUT, ITEST, IFDB)
      ELSEIF (.NOT.LSPLITOUT) THEN
        CALL SAVSPEC(FL1,NSTART,NEND,LDSTOP,LDWRRE,1,
     1               CDATEA,CDATEA,CDUM,CDUM, MARSTYPE,
     2               YCLASS, YEXPVER, CFDB2DSP)

      ENDIF
        IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. SAVSPEC DONE'
C
C
C ----------------------------------------------------------------------
C
C*    9. END OF JOB: DELETE WORK FILES.
C        ------------------------------
C
 9000 CONTINUE
C DEALLOCATE ARRAYS
C
      IF(ALLOCATED(U10OLD)) DEALLOCATE(U10OLD)
      IF(ALLOCATED(THWOLD)) DEALLOCATE(THWOLD)
      IF(ALLOCATED(USOLD)) DEALLOCATE(USOLD)
      IF(ALLOCATED(Z0OLD)) DEALLOCATE(Z0OLD)
      IF(ALLOCATED(TAUW)) DEALLOCATE(TAUW)
      IF(ALLOCATED(FL1)) DEALLOCATE(FL1)
      IF(ALLOCATED(NSTART)) DEALLOCATE(NSTART)
      IF(ALLOCATED(NEND)) DEALLOCATE(NEND)

C
      WRITE (IU06,*) ' '
      WRITE (IU06,*) ' PROGRAM PRESET: ALL DONE'
C
      END
