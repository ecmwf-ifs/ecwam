      SUBROUTINE READCUR (IU02, CDTCR, UCUR, VCUR, NCC, NRC)

C ----------------------------------------------------------------------
C
C**** *READCUR* - ROUTINE TO READ A CURRENT FIELD.
C
C     SUSANNE HASSELMANN        JUNE     1990.
C     H. GUNTHER    ECMWF/GKSS  DECEMBER 1990  MODIFIED FOR CYCLE_4.
C
C*    PURPOSE.
C     --------
C       INPUT OF A CURRENT FIELD.
C
C**   INTERFACE.
C     ----------
C       *CALL* *READCUR (IU02, CDTCR, UCUR, VCUR, NCC, NRC)*
C          *IU02*   - INPUT UNIT FOR CURRENTS.
C          *CDTCR*  - DATE/TIME OF THE DATA READ.
C          *UCUR*   - HORIZONTAL CURRENT COMPONENT.
C          *VCUR*   - VERTICAL CURRENT COMPONENT.
C          *NCC*    - NUMBER OF COLUMNS IN INPUT ARRAYS (DIMENSION).
C          *NRC*    - NUMBER OF ROWS    IN INPUT ARRAYS (DIMENSION).
C
C     METHOD.
C     -------
C
C       FORMATED READ; CURRENT GRID IS LAT / LONG GRID.
C
C     EXTERNALS.
C     ----------
C
C       *ABORT1*     - TERMINATES PROCESSING.
C       *ADJUST*    - CORRECTS LONGITUDE INPUT.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
#include "comcurg.h"
C
#include "comtest.h"
C
C ----------------------------------------------------------------------
C
      DIMENSION UCUR(NCC,NRC), VCUR(NCC,NRC)
      CHARACTER*12 CDTCR, ZERO
      PARAMETER(ZERO=' ') 
C
C ----------------------------------------------------------------------
C
C*    1. READ CURRENT FILE HEADER.
C        -------------------------
C
 1000 CONTINUE
C
      READ (IU02, '(6F10.5,2I3)',END=3001,ERR=3001)
     1    UCLO, UCLH, UCLW, UCLE, DUCLO, DUCLA, KRCI, KCCI
      IF (ITEST.GT.0) THEN
         WRITE (IU06,*) '  '
         WRITE (IU06,*) ' SUB. READCUR: FILE HEADER READ'
         WRITE (IU06,*) ' GRID SPECIFICATION ARE:'
         WRITE (IU06,*) ' NUMBER OF COLUMNS IN GRID KCCI = ', KCCI
         WRITE (IU06,*) ' NUMBER OF ROWS    IN GRID KRCI = ', KRCI
         WRITE (IU06,*) ' SOUTHERN MOST LATITUDE    UCLO = ', UCLO
         WRITE (IU06,*) ' NORTHERN MOST LATITUDE    UCLH = ', UCLH
         WRITE (IU06,*) ' WESTERN MOST LONGITUDE    UCLW = ', UCLW
         WRITE (IU06,*) ' EASTERN MOST LONGITUDE    UCLE = ', UCLE
         WRITE (IU06,*) ' LATITUDE  STEP IS        DUCLA = ', DUCLA
         WRITE (IU06,*) ' LONGITUDE STEP IS        DUCLO = ', DUCLO
      ENDIF
C
C*    1.2  CHECK DIMENSIONS
C     ----------------------
C
      IF (KCCI.GT.NCC .OR. KRCI.GT.NRC) THEN
         WRITE (IU06,*) ' *****************************************'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *       FATAL ERROR IN SUB. READCUR     *'
         WRITE (IU06,*) ' *       ===========================     *'
         WRITE (IU06,*) ' * DIMENSION OF INPUT ARRAYS ARE TO SMALL*'
         WRITE (IU06,*) ' * DIMENSION OF COLUMNS IS  NCC = ', NCC
         WRITE (IU06,*) ' * DIMENSION REQUESTED  IS KCCI = ', KCCI
         WRITE (IU06,*) ' * DIMENSION OF ROWS    IS  NRC = ', NRC
         WRITE (IU06,*) ' * DIMENSION REQUESTED  IS KRCI = ', KRCI
         WRITE (IU06,*) ' * CHANGE PARAMETER STATEMENT            *'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS      *'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *****************************************'
         CALL ABORT1
      ENDIF
      IF (KRCI.LT.2 .OR. KCCI.LT.2) THEN
         WRITE (IU06,*) ' *****************************************'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *     FATAL ERROR IN SUB. READCUR       *'
         WRITE (IU06,*) ' *     ===========================       *'
         WRITE (IU06,*) ' * LESS THAN 2 LATITUDES OR LONGITUDES   *'
         WRITE (IU06,*) ' * IN CURRENT INPUT.                     *'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS      *'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *****************************************'
         CALL ABORT1
      ENDIF

      CALL ADJUST (UCLW, UCLE)
      IPCUR=0
      IF (MOD(UCLE+DUCLO-UCLW+720.,360.).LT.0.0000001) IPCUR=1
C
C ----------------------------------------------------------------------
C
C*    2. READ CURRENT FIELD.
C        -------------------
C
 2000 CONTINUE
      CDTCR = ZERO
      DO 2001 J=1,KRCI
         READ(IU02,'(8E9.3)',END=3002,ERR=3002) (UCUR(I,J),I=1,KCCI)
 2001 CONTINUE

      DO 2002 J=1,KRCI
         READ(IU02,'(8E9.3)',END=3003,ERR=3003) (VCUR(I,J),I=1,KCCI)
 2002 CONTINUE
C
      RETURN
C
C ----------------------------------------------------------------------
C
C*    3. ERROR HANDLING.
C        ---------------
C
 3001 CONTINUE
         WRITE (IU06,*) ' *****************************************'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *      FATAL ERROR IN SUB. READCUR      *'
         WRITE (IU06,*) ' *      ===========================      *'
         WRITE (IU06,*) ' * READ ERROR OR EOF ON CURRENT FILE.    *'
         WRITE (IU06,*) ' * FILE HEADER EXPECTED                  *'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS      *'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *****************************************'
         CALL ABORT1
 3002 CONTINUE
         WRITE (IU06,*) ' *****************************************'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *      FATAL ERROR IN SUB. READCUR      *'
         WRITE (IU06,*) ' *      ===========================      *'
         WRITE (IU06,*) ' * READ ERROR OR EOF ON CURRENT FILE.    *'
         WRITE (IU06,*) ' * U - COMPONETS EXPECTED                *'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS      *'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *****************************************'
         CALL ABORT1
 3003 CONTINUE
         WRITE (IU06,*) ' *****************************************'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *      FATAL ERROR IN SUB. READCUR      *'
         WRITE (IU06,*) ' *      ===========================      *'
         WRITE (IU06,*) ' * READ ERROR OR EOF ON CURRENT FILE.    *'
         WRITE (IU06,*) ' * V - COMPONETS EXPECTED                *'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS      *'
         WRITE (IU06,*) ' *                                       *'
         WRITE (IU06,*) ' *****************************************'
         CALL ABORT1

      END
