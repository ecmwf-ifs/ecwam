      SUBROUTINE GETWND (US, DS, CDTWIS, UWND, VWND, NC, NR,
     1                   IUGB, IVGB, ILEN)

C ----------------------------------------------------------------------
C
C**** *GETWND* - ROUTINE TO READ AND PROCESS ONE WINDFIELD.
C
C*    PURPOSE.
C     --------
C
C        READ A WINDFIELD FROM THE WINDFILE (SEARCH FOR IT)
C        AND CALCULATES THE WIND VELOCITY  AND DIRECTION
C        FOR ALL WAM BLOCKS.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *GETWND (US, DS, CDTWIS, UWND, VWND, NC, NR)*
C         *US*     - MAGNITUDE OF USTAR IN EACH POINT AND BLOCK.
C         *DS*     - DIRECTION OF USTAR IN EACH POINT AND BLOCK.
C         *CDTWIS* - DATE OF WIND FIELD TO BE LOOKED FOR.
C         *UWND*   - INPUT WIND FIELD U COMPONENT.
C         *VWND*   - INPUT WIND FIELD V COMPONENT.
C         *NC*     - DIMENSION OF WIND INPUT ARRAYS NUMBER OF COLUMNES.
C         *NR*     - DIMENSION OF WIND INPUT ARRAYS NUMBER OF ROWS.
C
C     METHOD.
C     -------
C
C       NONE.
C
C     EXTERNALS.
C     ----------
C
C       *ABORT1*     - TERMINATES PROCESSING.
C       *READWND*   - READING WINDS.
C       *WAMWND*    - CALCULATE WIND IN WAM POINTS.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C
C    MODIFIED BY:
C    ------------
C    B. HANSEN    ECMWF 1997
C                 RESTRUCTURE CALL TO READWND.
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
C
#include "commpp.h"
C
#include "comcoup.h"
C
#include "txtmpp.h"
C
#include "comtest.h"
C
C ----------------------------------------------------------------------
      INTEGER IUGB(ILEN), IVGB(ILEN)
C ----------------------------------------------------------------------
C
      DIMENSION UWND(NC,NR), VWND(NC,NR)
      REAL, DIMENSION (NINF:NSUP,NBLO) :: US, DS
C
      CHARACTER*10 CDTWIS, CDTWIR
C
      SAVE CDTWIR
C
C ----------------------------------------------------------------------
C
C*    1. WIND DATA ARE READ
C        ------------------
C
 1000 CONTINUE
#if region=='s'
      CALL READWND_RSWAMP2 (CDTWIR, UWND, VWND, NC, NR)
#elif region=='g'  || region=='m'
      CALL READWND_MARS (CDTWIR, UWND, VWND, NC, NR,
     1                   IUGB, IVGB, ILEN)
#else
      WRITE(0,*) ' GETWND: NO CALL TO READWND FOUND. '
      CALL ABORT1
#endif
C
C*    2. CHECK WIND FIELD DATE
C        ----------------------
C
      IF (CDTWIR.LT.CDTWIS) THEN
C
C     2.1 DATE OF INPUT FIELD IS BEFORE REQUESTED DATE
C         TRY AGAIN
C         --------------------------------------------
         IF(.NOT.LWCOU) THEN
           GOTO 1000
C         --------------------------------------------
         ELSE
           WRITE (IU06,*) ' ******************************************'
           WRITE (IU06,*) ' *                                         *'
           WRITE (IU06,*) ' *      FATAL ERROR SUB. GETWND            *'
           WRITE (IU06,*) ' *      =======================            *'
           WRITE (IU06,*) ' * WIND DATE READ IS EARLIER THAN EXPECTED *'
           WRITE (IU06,*) ' * DATE READ IS     CDTWIR = ', CDTWIR
           WRITE (IU06,*) ' * DATE EXPECTED IS CDTWIS = ', CDTWIS
           WRITE (IU06,*) ' *                                         *'
           WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS        *'
           WRITE (IU06,*) ' *                                         *'
           WRITE (IU06,*) ' ******************************************'
           CALL ABORT1
         ENDIF
      ELSE IF (CDTWIR.GT.CDTWIS) THEN
C
C     2.2 DATE OF INPUT FIELD IS LATER THAN REQUESTED DATE
C         ------------------------------------------------
         WRITE (IU06,*) ' ******************************************'
         WRITE (IU06,*) ' *                                        *'
         WRITE (IU06,*) ' *      FATAL ERROR SUB. GETWND           *'
         WRITE (IU06,*) ' *      =======================           *'
         WRITE (IU06,*) ' * WIND DATE READ IS LATER THAN EXPECTED  *'
         WRITE (IU06,*) ' * DATE READ IS     CDTWIR = ', CDTWIR
         WRITE (IU06,*) ' * DATE EXPECTED IS CDTWIS = ', CDTWIS
         WRITE (IU06,*) ' *                                        *'
         WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS       *'
         WRITE (IU06,*) ' *                                        *'
         WRITE (IU06,*) ' ******************************************'
         CALL ABORT1
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    3. INTERPOLATE AND BLOCK WINDFIELD
C        -------------------------------
C
      CALL WAMWND (US, DS, UWND, VWND, NC, NR)

      RETURN
      END
