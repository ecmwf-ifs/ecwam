      SUBROUTINE READWND (CDTWIR, UWND, VWND, NC, NR,
     1                    IUGB, IVGB, ILEN)

C ----------------------------------------------------------------------
C
C**** *READWND* - ROUTINE TO READ WIND INPUT.
C
C     H. GUNTHER      GKSS/ECMWF    JULY 1990
C     J. DOYLE        NRL/ECMWF     OCT. 1996
C
C*    PURPOSE.
C     --------
C
C        READ WIND DATA.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *READWND (CDTWIR, UWND, VWND, NC, NR)*
C         *CDTWIR* - DATE/TIME OF THE DATA READ.
C         *UWND  * - HORIZONTAL WIND COMPONENS.
C         *VWND  * - VERTICAL WIND COMPONENTS.
C         *NC    * - NUMBER OF COLUMNS IN INPUT WIND ARRAY (DIMENSION).
C         *NR    * - NUMBER OF ROWS    IN INPUT WIND ARRAY (DIMENSION).
C
C     METHOD.
C     -------
C
C       THE SUB. EXPECTS A SEQUENTIAL, UNFORMATED WIND INPUT FILE
C       WHICH IS ASSIGNED TO UNIT IU01 = 1.
C       IN ITS PRESENT SET UP A HEADER IS ONLY READ WHEN THE
C       SUB IS CALLED THE FIRST TIME. INFORMATION ABOUT THE WINDS
C       AND THE WIND GRID IS STORED IN COMMON WNDGRD.
C       FOR EACH CALL A WINDFIELD IS READ TOGETHER WITH THE
C       DATE/TIME GROUP.
C       THE OUTPUT ARRAYS MUST CONTAIN THE WIND COMPONENTS AND MUST
C       BE ORGANISED AS FOLLOWS:
C
C                  ARRAY INDICES |   LONGITUDE  LATITUDE
C                  --------------|----------------------
C                      1     1   |     RLONL      RLATS
C                   KCOL     1   |     RLONR      RLATS
C                      1  KROW   |     RLONL      RLATN
C                   KCLO  KROW   |     RLONR      RLATN
C
C       NOTE: THIS SUBROUTINE HAS TO BE MODIFIED BY THE USER
C             FOR HIS SPECIAL WIND FORMATS.
C
C     EXTERNALS.
C     ----------
C
C       *ABORT1*     - TERMINATES PROCESSING.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
#include "comcoup.h"
C
#include "comwndg.h"
C
#include "comtest.h"
C
#include "txtwndg.h"
C
      DIMENSION UWND(NC,NR), VWND(NC,NR)
      INTEGER   KSEC1(25), KSEC2(22+NR), KSEC3(2), KSEC4(42), KSEC0(2),
C ---------------------------------------------------------------------
      INTEGER IUGB(ILEN), IVGB(ILEN)
C ---------------------------------------------------------------------
      DIMENSION PSEC2(10), PSEC3(2), PSEC4(NC*NR)
      CHARACTER*12 CDTWIR
      CHARACTER HOPER*1, FILNM*24
C
      DATA IFIL, HOPER / 63, 'R' /
C
      DATA IU01 /1/, ITEL /0/
C
C       *IU01*  -  WIND INPUT UNIT.
C       *ITEL*   - NUMBER OF WIND FIELDS READ.
C
C*    1. HEADER OF WIND DATA FILE AND CALCULATE STEPSIZES IN WINDGRID.
C        -------------------------------------------------------------
C
      IF(.NOT.LWCOU) THEN 
 1000 CONTINUE
      IF (ITEL.EQ.0) THEN
C
C*    1.1 READ HEADER OF WIND DATA FILE.
C         ------------------------------
C
         READ (IU01, ERR=3000, END=3000, IOSTAT=IOS)
     1        KCOL, KROW, RLATS, RLATN, RLONL, RLONR, ICOORD, IWPER
         DLAM  = (RLONR-RLONL)/REAL(KCOL-1)
         DPHI  = (RLATN-RLATS)/REAL(KROW-1)
         IF (ITEST.GT.0) THEN
            WRITE (IU06,*) '  '
            WRITE (IU06,*) ' SUB. READWND: WIND FILE HEADER READ'
            WRITE (IU06,*) ' WIND GRID SPECIFICATION ARE:'
            WRITE (IU06,*) ' NUMBER OF COLUMNS IN GRID KCOL = ', KCOL
            WRITE (IU06,*) ' NUMBER OF ROWS    IN GRID KROW = ', KROW
            WRITE (IU06,*) ' SOUTHERN MOST LATITUDE   RLATS = ', RLATS
            WRITE (IU06,*) ' NORTHERN MOST LATITUDE   RLATN = ', RLATN
            WRITE (IU06,*) ' WESTERN MOST LONGITUDE   RLONL = ', RLONL
            WRITE (IU06,*) ' EASTERN MOST LONGITUDE   RLONR = ', RLONR
            WRITE (IU06,*) ' COORDINATE SYSTEM CODE  ICOORD = ', ICOORD
            WRITE (IU06,*) ' PERIODIC GRID INDICATOR  IWPER = ', IWPER
            WRITE (IU06,*) ' LATITUDE  STEP IS         DPHI = ', DPHI
            WRITE (IU06,*) ' LONGITUDE STEP IS         DLAM = ', DLAM
         ENDIF
C
C*    1.2 CHECK DIMENSIONS.
C         -----------------
C                            
         NCNR = NC*NR
C
         IF (KCOL*KROW.GT.NCNR) THEN
            WRITE (IU06,*) ' *****************************************'
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READWND     *'
            WRITE (IU06,*) ' *       ===========================     *'
            WRITE (IU06,*) ' * DIMENSION OF WIND ARRAY IS TO SMALL   *'
            WRITE (IU06,*) ' * DIMENSION WIND ARRAY IS     NCNR = ',
     1                     NCNR
            WRITE (IU06,*) ' * DIMENSION REQUESTED IS KCOL*KROW = ',
     1                     KCOL*KROW
            WRITE (IU06,*) ' * CHANGE PARAMETER STATEMENT IN WAVMDL  *'
            WRITE (IU06,*) ' * VALUES HAVE TO BE:                    *'
            WRITE (IU06,*) ' *      NC = ', KCOL
            WRITE (IU06,*) ' *      NR = ', KROW
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS      *'
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' *****************************************'
            CALL ABORT1
         ENDIF
         IF (KROW.LT.2 .OR. KCOL.LT.2) THEN
            WRITE (IU06,*) ' *****************************************'
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READWND     *'
            WRITE (IU06,*) ' *       ===========================     *'
            WRITE (IU06,*) ' * DIMENSIONS OF WIND INPUT GRID ARE     *'
            WRITE (IU06,*) ' * LESS THAN 2.                          *'
            WRITE (IU06,*) ' * NUMBER OF COLUMNS   KCOL = ', KCOL
            WRITE (IU06,*) ' * NUMBER OF ROWS      KROW = ', KROW
            WRITE (IU06,*) ' * CHECK WIND INPUT AND MAKE SURE, THAT  *'
            WRITE (IU06,*) ' * PARAMETER STATEMENT IN WAVMDL IS OK   *'
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS      *'
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' *****************************************'
            CALL ABORT1
         ENDIF
C
C*    1.3 CHECK PERIODIC INDICATOR.
C         -------------------------
C
         IF (IWPER.NE.0 .AND. IWPER.NE.1) THEN
            WRITE (IU06,*) ' *****************************************'
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READWND     *'
            WRITE (IU06,*) ' *       ===========================     *'
            WRITE (IU06,*) ' * PERIODIC GRID INDICATOR HAS TO BE:    *'
            WRITE (IU06,*) ' *   0  IF GRID IS NOT PERIODIC          *'
            WRITE (IU06,*) ' *   1  IF GRID IS     PERIODIC          *'
            WRITE (IU06,*) ' * VALUE READ FROM WIND FILE IS IWPER = ',
     1                     IWPER
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS      *'
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' *****************************************'
            CALL ABORT1
         ENDIF
C
C*    1.4 COORDINATE SYSTEM INDICATOR.
C         ----------------------------
C
         IF (ICOORD.NE.1) THEN
            WRITE (IU06,*) ' *****************************************'
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READWND     *'
            WRITE (IU06,*) ' *       ===========================     *'
            WRITE (IU06,*) ' * COORDINATE SYSTEM INDICATOR HAS TO BE *'
            WRITE (IU06,*) ' *   1 FOR EQUIDISTANT LAT/LONG GRID     *'
            WRITE (IU06,*) ' *   OTHER GRIDS ARS NOT INSTALLED       *'
            WRITE (IU06,*) ' * VALUE READ FROM WIND FILE IS ICOORD = ',
     1                     ICOORD
            WRITE (IU06,*) ' * CHECK WIND INPUT OR UPDATE SUBROUTINE *'
            WRITE (IU06,*) ' * WAMWND AND THIS CHECK                 *'
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS      *'
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' *****************************************'
            CALL ABORT1
         ENDIF
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    2. READ WIND DATA.
C        ---------------
C
 2000 CONTINUE
      ITEL = ITEL + 1
C
C*    2.1 DATE / TIME OF WIND.
C         --------------------
C
      READ (IU01, ERR=3100, END=3100, IOSTAT=IOS) CDTWIR
      IF (ITEST.GT.0) THEN
         WRITE (IU06,*) '  '
         WRITE (IU06,*) ' SUB. READWND: NEW WIND DATE READ CDTWIR = ',
     1                   CDTWIR
      ENDIF
C
C*    2.2 READ HORIZONTAL COMPONENTS.
C         ---------------------------
C
      READ (IU01, ERR=3200, END=3200, IOSTAT=IOS)
     1     ICODE, ((UWND(I,J),I=1,KCOL),J=1,KROW)
      IF (ITEST.GT.0) WRITE (IU06,*) ' U: ',(UWND(I,1),I=1,10)
C
C*    2.3 READ VERTICAL COMPONENTS.
C         -------------------------
C
      READ (IU01, ERR=3300, END=3300, IOSTAT=IOS)
     1     KCODE, ((VWND(I,J),I=1,KCOL),J=1,KROW)
      IF (ITEST.GT.0) WRITE (IU06,*) ' V: ',(VWND(I,1),I=1,10)
C
C*    2.4 CHECK WIND CODES.
C         -----------------
C
      IF (ICODE.NE.KCODE) THEN
         WRITE (IU06,*) ' *****************************************'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READWND     *'
         WRITE (IU06,*) ' *       ===========================     *'
         WRITE (IU06,*) ' * CODES OF U AND V WIND COMPONENT ARE   *'
         WRITE (IU06,*) ' * DIFFERENT.                            *'
         WRITE (IU06,*) ' * CODE OF U COMPONENT IS ICODE = ', ICODE
         WRITE (IU06,*) ' * CODE OF V COMPONENT IS KCODE = ', KCODE
         WRITE (IU06,*) ' * DATE OF WIND INPUT IS CDTWIR = ', CDTWIR
         WRITE (IU06,*) ' * CODE HAS TO BE :                      *'
         WRITE (IU06,*) ' *    1  FOR FRICTION VELOCITIES         *'
         WRITE (IU06,*) ' *    2  FOR SURFACE STRESSES            *'
         WRITE (IU06,*) ' *    3  FOR WIND IN 10M HEIGHT          *'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS      *'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *****************************************'
         CALL ABORT1
      ENDIF
      IF (ICODE.NE.1 .AND. ICODE.NE.2 .AND. ICODE.NE.3) THEN
         WRITE (IU06,*) ' *****************************************'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READWND     *'
         WRITE (IU06,*) ' *       ===========================     *'
         WRITE (IU06,*) ' * WIND FIELD CODE IS NOT ALLOWED        *'
         WRITE (IU06,*) ' * CODE OF U COMPONENT IS ICODE = ', ICODE
         WRITE (IU06,*) ' * CODE OF V COMPONENT IS KCODE = ', KCODE
         WRITE (IU06,*) ' * DATE OF WIND INPUT IS CDTWIR = ', CDTWIR
         WRITE (IU06,*) ' * CODE HAS TO BE :                      *'
         WRITE (IU06,*) ' *    1  FOR FRICTION VELOCITIES         *'
         WRITE (IU06,*) ' *    2  FOR SURFACE STRESSES            *'
         WRITE (IU06,*) ' *    3  FOR WIND IN 10M HEIGHT          *'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS      *'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *****************************************'
         CALL ABORT1
      ENDIF

      RETURN
C
C ----------------------------------------------------------------------
C
C*    3. ERROR MESSAGES.
C        ---------------
C
 3000 CONTINUE
         WRITE (IU06,*) ' *****************************************'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READWND     *'
         WRITE (IU06,*) ' *       ===========================     *'
         WRITE (IU06,*) ' * READ ERROR OR EOF ON WIND FILE:       *'
         WRITE (IU06,*) ' * FILE HEADER EXPECTED                  *'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS      *'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *****************************************'
         CALL ABORT1
 3100 CONTINUE
         WRITE (IU06,*) ' *****************************************'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READWND     *'
         WRITE (IU06,*) ' *       ===========================     *'
         WRITE (IU06,*) ' * READ ERROR OR EOF ON WIND FILE:       *'
         WRITE (IU06,*) ' * WIND DATE EXPECTED                    *'
         WRITE (IU06,*) ' * PROGRAM TRIES TO READ THE ', ITEL,
     1                  ' WINDFIELD'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS      *'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *****************************************'
         CALL ABORT1
 3200 CONTINUE
         WRITE (IU06,*) ' *****************************************'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READWND     *'
         WRITE (IU06,*) ' *       ===========================     *'
         WRITE (IU06,*) ' * READ ERROR OR EOF ON WIND FILE:       *'
         WRITE (IU06,*) ' * U COMPONENT EXPECTED                  *'
         WRITE (IU06,*) ' * DATE OF WIND INPUT IS CDTWIR = ', CDTWIR
         WRITE (IU06,*) ' * PROGRAM TRIES TO READ THE ', ITEL,
     1                  ' WINDFIELD'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS      *'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *****************************************'
         CALL ABORT1
 3300 CONTINUE
         WRITE (IU06,*) ' *****************************************'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READWND     *'
         WRITE (IU06,*) ' *       ===========================     *'
         WRITE (IU06,*) ' * READ ERROR OR EOF ON WIND FILE:       *'
         WRITE (IU06,*) ' * V COMPONENT EXPECTED                  *'
         WRITE (IU06,*) ' * DATE OF WIND INPUT IS CDTWIR = ', CDTWIR
         WRITE (IU06,*) ' * PROGRAM TRIES TO READ THE ', ITEL,
     1                  ' WINDFIELD'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS      *'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *****************************************'
         CALL ABORT1

      ELSE
C
C*    3.0 READ ATMOS. MODEL GRIB FILES FOR COUPLED MODEL
C         -------------------
C
      DO 4000 IVAR=1,2
C         -------------------
C
      KWORD = LENGTH(IFIL)
      KLENP = NC*NR
      KRET  = 1
C     CALL GRSDBG(1)
      IF(IVAR.EQ.1) THEN
        CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     1               PSEC4, KLENP, IUGB,  ILEN , KWORD, HOPER, KRET)
      ELSE
        CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     1               PSEC4, KLENP, IVGB,  ILEN , KWORD, HOPER, KRET)
      ENDIF
      IF(KRET.NE.0) GOTO 5200
C         
C*    3.1 PRODUCT DEFINITION BLOCK.
C         -------------------------
C
      IF (ITEST.GT.1)  CALL GRPRS1 (KSEC0, KSEC1)
C
C*    3.1.1 DETERMINE DATE OF WIND FIELD FROM KSEC1(10)-KSEC1(14).
C
C     START DATE
C
      IYYYY=(KSEC1(21)-1)*100+KSEC1(10)
      WRITE(CDTWIR,'(I4.4,4I2.2)') IYYYY,KSEC1(11),KSEC1(12),KSEC1(13),
     1                             KSEC1(14)
C
C     DETERMINE TIME UNIT
C
      IF (KSEC1(15).EQ.0) THEN
         JCONS = 60
      ELSE IF (KSEC1(15).EQ.1) THEN
         JCONS = 3600
      ELSE IF (KSEC1(15).EQ.2) THEN
         JCONS = 86400
      ELSE
         WRITE(IU06,*) ' *********************************************'
         WRITE(IU06,*) ' *                                           *'
         WRITE(IU06,*) ' *         FATAL ERROR IN SUB. READWND       *'
         WRITE(IU06,*) ' *         ===========================       *'
         WRITE(IU06,*) ' * SUSPICIOUS TIME UNIT IN MARS BLOCK 1      *'
         WRITE(IU06,*) ' * TIME UNIT IS KSEC1(15) = ', KSEC1(15)
         WRITE(IU06,*) ' *                                           *'
         WRITE(IU06,*) ' *    PROGRAM ABORTS     PROGRAM ABORTS      *'
         WRITE(IU06,*) ' *                                           *'
         WRITE(IU06,*) ' *********************************************'
         CALL ABORT
      ENDIF
C         
C     DETERMINE WIND FIELD DATE.
C         
      CALL INCDATE(CDTWIR,JCONS*KSEC1(16))
C         
C         
C*    3.1.2 DETERMINE CODE FOR WIND FIELD TYPE:
C           165 = U10 FIELD , 180 = STRESS FIELD.
C           ------------------------------------
C         
      IF (KSEC1(6).EQ.165 .OR. KSEC1(6).EQ.166) THEN
         ICODE = 3
      ELSE IF (KSEC1(6).EQ.180 .OR. KSEC1(6).EQ.181) THEN
         ICODE = 2
      ELSE
         WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++'
         WRITE(IU06,*) ' +                                           +'
         WRITE(IU06,*) ' +       WARNING ERROR IN SUB. READWND       +'
         WRITE(IU06,*) ' +       =============================       +'
         WRITE(IU06,*) ' + SUSPICIOUS WIND FIELD TYPE IN MARS BLOCK 1+'
         WRITE(IU06,*) ' + TYPE IS KSEC1(6) = ', KSEC1(6)
         WRITE(IU06,*) ' +                                           +'
         WRITE(IU06,*) ' +    PROGRAM CONTINUES                      +'
         WRITE(IU06,*) ' +                                           +'
         WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++'
         GOTO 2000
      ENDIF
C         
C*    3.2 GRID DEFINITION BLOCK.
C         ---------------------
C         
      IF (ITEST.GT.1) CALL GRPRS2 (KSEC0, KSEC2, PSEC2)
C         
C*    3.2.1 DETERMINE GRID PARAMETERS.
C           --------------------------
C         
      KCOL = KSEC2(2)
      KROW = KSEC2(3)
      IF ((KCOL.GT.NC).OR.(KROW.GT.NR)) THEN
         WRITE(IU06,*) ' *********************************************'
         WRITE(IU06,*) ' *                                           *'
         WRITE(IU06,*) ' *         FATAL ERROR IN SUB. READWND       *'
         WRITE(IU06,*) ' *         ===========================       *'
         WRITE(IU06,*) ' * DIMENSION DO NOT FIT                      *'
         WRITE(IU06,*) ' * REQUIRED IS     KCOL = ', KCOL
         WRITE(IU06,*) ' * DIMENSION IS      NC = ', NC
         WRITE(IU06,*) ' * REQUIRED IS     KROW = ', KROW
         WRITE(IU06,*) ' * DIMENSION IS      NR = ', NR
         WRITE(IU06,*) ' *                                           *'
         WRITE(IU06,*) ' *    PROGRAM ABORTS     PROGRAM ABORTS      *'
         WRITE(IU06,*) ' *                                           *'
         WRITE(IU06,*) ' *********************************************'
         CALL ABORT
      END IF
      RLATN=.001*FLOAT(KSEC2(4))
      RLONL=.001*FLOAT(KSEC2(5))
      RLATS=.001*FLOAT(KSEC2(7))
      RLONR=.001*FLOAT(KSEC2(8))
C
      DLAM  = (RLONR-RLONL)/REAL(KCOL-1)
      DPHI  = (RLATN-RLATS)/REAL(KROW-1)
C         
C --------------------------------------------------------------------
C         
C*    4. REARRANGE WIND DATA.
C        --------------------
C         
      NR1 = KROW+1
      IF (KSEC1(6).EQ.165 .OR. KSEC1(6).EQ.180) THEN
         L = 0
         DO 3100 J=1,KROW
         DO 3100 I=1,KCOL
            L = L+1
            UWND(I,NR1-J)=PSEC4(L)
 3100    CONTINUE
      ELSE IF (KSEC1(6).EQ.166 .OR. KSEC1(6).EQ.181) THEN
         L = 0
         DO 3200 J=1,KROW
         DO 3200 I=1,KCOL
            L = L+1
            VWND(I,NR1-J)=PSEC4(L)
 3200    CONTINUE
      ENDIF
C
 4000 CONTINUE
C
      RETURN
C
 5200 CONTINUE
         WRITE(IU06,*) ' *********************************************'
         WRITE(IU06,*) ' *                                           *'
         WRITE(IU06,*) ' *         FATAL ERROR IN SUB. READWND       *'
         WRITE(IU06,*) ' *         ===========================       *'
         WRITE(IU06,*) ' * GRIBEX ERROR - KRET = ',KRET
         WRITE(IU06,*) ' *                                           *'
         WRITE(IU06,*) ' *    PROGRAM ABORTS     PROGRAM ABORTS      *'
         WRITE(IU06,*) ' *                                           *'
         WRITE(IU06,*) ' *********************************************'
         CALL ABORT
C
      ENDIF
      RETURN
C
      END
