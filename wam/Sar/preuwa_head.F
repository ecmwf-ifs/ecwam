      program preuwa_head
C
C**** *preuwa_head*
C
C     PURPOSE.
C     --------
C          generates a list with output positions for the WAModel   
C          in order to produce a collocation file of wave spectra 
C          to be used for the SAR-INVERSION system.
C
C**   INTERFACE.
C     ----------
C      INPUT:
C        DATE     -       an environment variable (YYMMDDhhmm) indi-
C                         cating the time of which a collocation file
C                         is to be produced.
C        GRID     -       the output file on unit 7 of PREPROC.
C        PARWAM   -       a file with a namelist of the main.
C                         arry dimensions used by the WAModel.
C        UWAyymmddhhmm  - the FDP SAR DATA.
C                 
C      OUTPUT:
C        COLyymmddhhmm  - the list with output positions.
C        SARyymmddhhmm  - the total content of FDP UWAX decoded.
C        SERyymmddhhmm  - Control output and warnings from ers1sar.
C
C      SCRATCH FILE:
C        VALyymmddhhmm  - intermedite data file.
C
C     METHOD.
C     -------
C          driven by DATE the file UWAXyymmddhhmm and COLyymmddhhmm
C          are opened. Using the lat - long information of SAR spectra
C          FINDB is used to get the corresponding block and block index 
C          of the nearest WAModel grid point. This information is 
C          written to the output file.
C
C     EXTERNALS.
C     ----------
C          DECOUWA  - generates a list of output points
C          ERS1SAR  - generates decoded and reduced FDP UWA data for
C                     SAR inversion..
C
C     REFERENCE.
C     ----------
C          NONE.
C
C     AUTHOR.
C     -------
C          B. HANSEN        *ECMWF*       10/02/91.
C
C     MODIFICATIONS.
C     --------------
C          B. HANSEN        *ECMWF*       FEB 93    
C            Top routine preuwa_head added in order to be able
C            to read the array dimensions of the WAModel arrays
C            from a namelist. 
C            ECFILE call inside the program dismissed.
C            Per call only one product will be processed. Loops
C            over time have to be done on jcl-level.
C
C----------------------------------------------------------------------
C
      CHARACTER* 20 cpuwam, cpugrd
      PARAMETER (
     .   jpuso  =         6, ! LOGICAL UNIT PRINTER OUTPUT
     .   cpuwam = 'PARWAM' , ! HOLORITH UNIT FOR FILE GRID
     .   cpugrd =   'GRID' , ! HOLORITH UNIT FOR FILE GRID
     .   jpuwax =        10, ! LOGICAL UNIT FOR PRODUCT UWAyymmddhhmm
     .   jpucol =        11, ! LOGICAL UNIT FOR FILE COLyymmddhhmm
     .   jpuval =        12, ! LOGICAL UNIT FOR FILE VALyymmddhhmm
     .   jpusar =        13, ! LOGICAL UNIT FOR FILE SARyymmddhhmm
     .   jpuerro=        14  ! LOGICAL UNIT FOR FILE SERyymmddhhmm
     .          )
C
C               
      CHARACTER* 3 youwa             ! FILE ID FOR THE SAR INPUT FILE.
      CHARACTER* 3 yocol             ! FILE ID FOR THE COLLOCATION FILE.
      CHARACTER* 3 yoval             ! FILE ID FOR THE SCRATCH SAR DATA.
      CHARACTER* 3 yosar             ! FILE ID FOR THE DECODED SAR DATA.
      CHARACTER* 3 yoerro            ! FILE ID FOR THE CONTROL OUTPUT.  
      CHARACTER*10 yohi              ! FOR THE yymmddhhss part of file
C                                      names.
      DATA youwa, yocol, yoval, yosar, yoerro 
     .   / 'UWA', 'COL', 'VAL', 'SAR', 'SER' /

      INTEGER f_getenv
      CHARACTER*256 YOISCOM
      CHARACTER*100 PATH
C
C ----------------------------------------------------------------------
C
C*       PARAMETER NAMELIST AS IN COMMON PARAM OF THE WAM-MODEL.
C        -------------------------------------------------------
C
      LOGICAL llper          !  the land sea mask is periodic (T)
      NAMELIST /PARWAM/
     .  iang, ifre, igx, igy, iblo, iiblo, iover, ioutp, ioutt,
     .  imaxc, imaxf, ibinp, iibl1, iibld, ibld, iiblc, iblc ,
     .  itaumax, iumax, iustar, ialpha, idepth, irefra, llper
C
C ----------------------------------------------------------------------
C
#include <outctrl.h>
      DATA csubna /"uwapreh"/
C
C ----------------------------------------------------------------------
C
C* 1.    INITIALIZATIONS.
C  ----------------------
 1000 CONTINUE
      msub = 1     ! DO NOT CHANGE THIS VALUE !
C
C   THERE WILL BE STANDARD OUTPUT OF EACH SUB-ROUTINE DOWN
C   TO THE LEVEL msuvi
C   THERE WILL BE TEST OUTPUT OF EACH SUB-ROUTINE DOWN
C   TO THE LEVEL mtest
      msuvi = 2
      mtest = 0
C
C   THERE WILL BE STANDARD OUTPUT OF  EACH
C   SUB-ROUTINE MENTIONED IN THE LIST "csuvi"
      csuvi = ' '
C
C   THERE WILL BE TEST OUTPUT OF  EACH
C   SUB-ROUTINE MENTIONED IN THE LIST "ctest"
      ctest = ' '
C
      CALL uraopc (csubna, msub, llsuvi, lltest)
C
      IF (llsuvi) THEN
        WRITE(jpuso,'(a)')' '
        WRITE(jpuso,'(1h1,a)')
     .  'UWAPRE----------------------------------------'//
     .  '--------------------------------'
        WRITE(jpuso,'(a)')' '
      ENDIF

 1001 FORMAT(A50)
 1002 FORMAT(A50,I10)
 1003 FORMAT(A6,A13,A31)
C
C-----------------------------------------------------------------------
C
C* 2.    GET DATE FROM ENVIRONMENT VARIABLE DATE.
C  ----------------------------------------------
C
      yoiscom='DATE'
      ir=f_getenv(yoiscom,yohi)
      IF (llsuvi) WRITE(jpuso,*)'  CURRENT DATE IS >',yohi,'<'
C
C* 3.    READ WAMODEL-PARAMETERS FROM FILE PARWAM.
C  -----------------------------------------------
C
      iuparwam=i_get_unit(jpuso, cpuwam, 's', 'f', 0)
      rewind iuparwam  
      READ(iuparwam, PARWAM)
C
C* 4.    OPEN THE REQUIRED FILES FOR LATER USE.
C  --------------------------------------------
C
      OPEN(jpuwax,FILE=youwa//yohi,status='old',
     .     iostat=ios,form='unformatted')
      IF (ios .NE. 0)THEN
        IF (llsuvi) THEN
          WRITE(jpuso,1003)' FILE ',youwa//yohi,' not opened    '
          WRITE(jpuso,1002)' iostat: ',ios
          WRITE(jpuso,1001)' this is a fatal error              '
          WRITE(jpuso,1001)' PROGRAM STOPS   NOTHING DONE       '
          STOP 'NO UWAX'
        ENDIF
      ENDIF
C
      OPEN(jpucol,FILE=yocol//yohi,status='new',
     .     iostat=ios,form='unformatted')
      IF (ios .NE. 0)THEN
        IF (llsuvi) THEN
          WRITE(jpuso,1003)' FILE ',yocol//yohi,' not opened    '
          WRITE(jpuso,1002)' iostat: ',ios
          WRITE(jpuso,1001)' this is a fatal error              '
          WRITE(jpuso,1001)' PROGRAM STOPS   NOTHING DONE       '
          STOP 'COL ? '
        ENDIF
      ELSE
        WRITE(jpucol) yohi
      ENDIF
C    
      OPEN(jpuval,FILE=yoval//yohi,status='new',
     .     iostat=ios,form='unformatted')
      IF (ios .NE. 0)THEN
        IF (llsuvi) THEN
          WRITE(jpuso,1003)' FILE ',yoval//yohi,' not opened    '
          WRITE(jpuso,1002)' iostat: ',ios
          WRITE(jpuso,1001)' this is a fatal error              '
          WRITE(jpuso,1001)' PROGRAM STOPS   NOTHING DONE       '
          STOP 'VAL ? '
        ENDIF
      ELSE
        WRITE(jpuval) yohi
      ENDIF
C    
      OPEN(jpusar,FILE=yosar//yohi,status='new',
     .     iostat=ios,form='unformatted')
      IF (ios .NE. 0)THEN
        IF (llsuvi) THEN
          WRITE(jpuso,1003)' FILE ',yosar//yohi,' not opened    '
          WRITE(jpuso,1002)' iostat: ',ios
          WRITE(jpuso,1001)' this is a fatal error              '
          WRITE(jpuso,1001)' PROGRAM STOPS   NOTHING DONE       '
          STOP 'SAR ? '
        ENDIF
      ELSE
        WRITE(jpusar) yohi
      ENDIF
C    
      OPEN(jpuerro,FILE=yoerro//yohi,status='new',
     .     iostat=ios,form='formatted')
      IF (ios .NE. 0)THEN
        IF (llsuvi) THEN
          WRITE(jpuso,1003)' FILE ',yoerro//yohi,' not opened   '
          WRITE(jpuso,1002)' iostat: ',ios
          WRITE(jpuso,1001)' this is a fatal error              '
          WRITE(jpuso,1001)' PROGRAM STOPS   NOTHING DONE       '
          STOP 'SER ? '
        ENDIF
      ELSE
        WRITE(jpuerro,'(a11)') yohi
      ENDIF
C
C* 5.    GENERATE THE LIST.
C  ------------------------
C
      iugrd=i_get_unit(jpuso, cpugrd, 's', 'u', 0)
      CALL uwapre (jpuso, iugrd, jpuwax, jpucol, jpuval, jpusar,
     .     jpuerro,
     .     iang, ifre, igx, igy, iblo, iiblo, iover, ioutp, ioutt,
     .     imaxc, imaxf, ibinp, iibl1, iibld, ibld, iiblc, iblc ,
     .     itaumax, iumax, iustar, ialpha, idepth, irefra,
     .     kfail)
      END
