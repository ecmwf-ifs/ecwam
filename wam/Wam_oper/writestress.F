      SUBROUTINE WRITESTRESS(IJINF,IJSUP,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     1                       FILENAME,LPBOPEN,LPBCLOSE) 

C ----------------------------------------------------------------------
C     J. BIDLOT    ECMWF      APRIL 1997
C
C*    PURPOSE.
C     --------
C     WRITES THE RESTART WIND AND STRESS FIELDS TO FILE.
C
C**   INTERFACE.
C     ----------
C
C     CALL *WRITESTRESS*(IJINF,IJSUP,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
C                        FILENAME,LPBOPEN,LPBCLOSE) 
C     *IJINF*    FIRST LOWER DIMENSION BOUND OF U10OLD,THWOLD,USOLD,TAUW,Z0OLD 
C     *IJSUP*    FIRST UPPER DIMENSION BOUND OF U10OLD,THWOLD,USOLD,TAUW,Z0OLD 
C     *U10OLD*   WIND SPEED.
C     *THWOLD*   WIND DIRECTION (RADIANS).
C     *USOLD*    FRICTION VELOCITY.
C     *TAUW*     WAVE STRESS.
C     *Z0OLD*    ROUGHNESS LENGTH IN.
C     *FILENAME* FILENAME (INCLUDING PATH) OF TARGET FILE
C     *LPBOPEN*  LOGICAL, TRUE IF PBOPEN HAS TO BE CALLED
C     *LPBCLOSE* LOGICAL, TRUE IF PBCLOSE HAS TO BE CALLED
C
C     METHOD.
C     -------
C     WRITES ARRAYS TO FILE. FORTRAN WRITES OR PBIO WRITES ARE
C     USED DEPENDING ON THE CHOICE OF OUTPUT MODE AS SELECTED BY LPBIOOUT
C
C     IF PBIO IS USED, THEN THE ARRAYS  WILL BE SAVED IN FILENAME 
C     USING PBWRITE ON CONDITION THAT FILENAME WAS NOT OPENED BEFORE
C     AND CONNECTED TO IUNIT.
C     IT WILL BE APPENDED TO THE BOTTOM OF THE CONTENT OF FILENAME IF
C     FILENAME WAS ALREADY OPENED AND CONNECTED TO UNIT IUNIT.
C
C     IF FORTRAN WRITES ARE USED THEN IN THE CASE OF MESSAGE PASSING
C     THE OUTPUT IS WRITTEN AS UNFORMATTED BINARY TO FILENAME CONNECTED
C     TO UNIT IU15. IN THE NON MESSAGE PASSING CASE, IT IS WRITTEN AS
C     UNFORMATTED BINARY TO FORTRAN UNIT SPECIFIED BY IU15.
C
C
C     EXTERNALS.
C     ----------
C     *PBOPEN*
C     *PBWRITE*
C     *PBCLOSE*
C     *ABORT1*
C
C     REFERENCE.
C     ----------
C     NONE
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "comstat.h"
C
#include "comtest.h"
C
#include "comunit.h"
C
#include "comwind.h"
C
#include "commesspass.h"
C
#include "commpp.h"

      CHARACTER*1 MODE
      CHARACTER*80 FILENAME
      LOGICAL LPBOPEN,LPBCLOSE
C
C ----------------------------------------------------------------------
      REAL,DIMENSION(IJINF:IJSUP,NBLO) :: U10OLD,THWOLD,USOLD,
     1                                    Z0OLD,TAUW
C ----------------------------------------------------------------------

      IF(LPBIOOUT) THEN

        IF(LPBOPEN) THEN
          MODE = 'w'
          CALL PBOPEN(IUNIT,FILENAME,MODE,KRET)
          IF(KRET.LT.0) THEN
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBOPEN       *'
          WRITE (IU06,*) '*   IN WRITESTRESS                       *' 
          IF(KRET.EQ.-1)
     1        WRITE (IU06,*) 'COULD NOT OPEN FILE ',FILENAME
          IF(KRET.EQ.-2) WRITE (IU06,*) 'INVALID FILENAME ',FILENAME
          IF(KRET.EQ.-3) WRITE (IU06,*) 'INVALID OPEN MODE SPECIFIED'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          CALL ABORT1
          ENDIF
        ENDIF

        CALL PBWRITE(IUNIT,CDTPRO,10,KRET)
        IF(KRET.LT.0) GOTO 1099
        CALL PBWRITE(IUNIT,CDATEWO,10,KRET)
        IF(KRET.LT.0) GOTO 1099
        CALL PBWRITE(IUNIT,CDAWIFL,10,KRET)
        IF(KRET.LT.0) GOTO 1099
        CALL PBWRITE(IUNIT,CDATEFL,10,KRET)
        IF(KRET.LT.0) GOTO 1099

        KOUNT = (IJSUP-IJINF+1)*NBLO*NPRECI
        CALL PBWRITE(IUNIT,U10OLD,KOUNT,KRET)
        IF(KRET.LT.0) GOTO 1099
        CALL PBWRITE(IUNIT,THWOLD,KOUNT,KRET)
        IF(KRET.LT.0) GOTO 1099
        CALL PBWRITE(IUNIT,USOLD,KOUNT,KRET)
        IF(KRET.LT.0) GOTO 1099
        CALL PBWRITE(IUNIT,TAUW,KOUNT,KRET)
        IF(KRET.LT.0) GOTO 1099
        CALL PBWRITE(IUNIT,Z0OLD,KOUNT,KRET)
        IF(KRET.LT.0) GOTO 1099

1099    CONTINUE 
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '*****************************************'
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*    ERROR FOLLOWING CALL TO PBWRITE    *'
          WRITE (IU06,*) '*    IN WRITESTRESS                     *' 
          WRITE (IU06,*) '*    FILE ',FILENAME 
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*****************************************'
          CALL ABORT1
        ENDIF


        IF(LPBCLOSE) THEN
          CALL PBCLOSE(IUNIT,KRET)
          IF(KRET.LT.0) THEN
             WRITE (IU06,*) '************************************'
             WRITE (IU06,*) '*                                  *'
             WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
             WRITE (IU06,*) '* IN WRITESTRESS                   *'
             WRITE (IU06,*) '* FILE ',FILENAME
             WRITE (IU06,*) '*                                  *'
             WRITE (IU06,*) '************************************'
             CALL ABORT1
          ENDIF
        ENDIF


      ELSE
         IF(LMESSPASS) THEN

            OPEN(IU15,FILE=FILENAME,STATUS='UNKNOWN',FORM='UNFORMATTED')

            WRITE(IU15) CDTPRO,CDATEWO,CDAWIFL,CDATEFL
            WRITE(IU15)((U10OLD(J1,IG),J1=IJINF,IJSUP),IG=1,NBLO)
            WRITE(IU15)((THWOLD(J1,IG),J1=IJINF,IJSUP),IG=1,NBLO)
            WRITE(IU15)((USOLD(J1,IG),J1=IJINF,IJSUP),IG=1,NBLO)
            WRITE(IU15)((TAUW(J1,IG),J1=IJINF,IJSUP),IG=1,NBLO)
            WRITE(IU15)((Z0OLD(J1,IG),J1=IJINF,IJSUP),IG=1,NBLO)

            CLOSE(IU15)

         ELSE

            REWIND IU15
            WRITE(IU15) CDTPRO,CDATEWO,CDAWIFL,CDATEFL
            WRITE(IU15) U10OLD
            WRITE(IU15) THWOLD
            WRITE(IU15) USOLD
            WRITE(IU15) TAUW
            WRITE(IU15) Z0OLD

         ENDIF
      ENDIF
      IF (ITEST.GE.2)
     1    WRITE(IU06,*) ' SUB. WRITESTRESS: OUTPUT DONE'
C
      RETURN
      END
