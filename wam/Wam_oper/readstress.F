      SUBROUTINE READSTRESS(RU10OLD,RTHWOLD,RUSOLD,RTAUW,RZ0OLD,
     1                      FILENAME)

C ----------------------------------------------------------------------
C     J. BIDLOT    ECMWF      MARCH 1997
C
C*    PURPOSE.
C     --------
C     READS THE RESTART WIND AND STRESS FIELDS
C
C**   INTERFACE.
C     ----------
C     CALL *READSTRESS*(RU10OLD,RTHWOLD,RUSOLD,RTAUW,RZ0OLD,FILENAME)
C     *RU10OLD*   WIND SPEED.
C     *RTHWOLD*   WIND DIRECTION (RADIANS).
C     *RUSOLD*    FRICTION VELOCITY.
C     *RTAUW*     WAVE STRESS.
C     *RZ0OLD*    ROUGHNESS LENGTH IN.
C     *FILENAME* FILENAME (INCLUDING PATH) OF INPUT FILE (ONLY IF MESSAGE 
C                PASSING IS USED)
C
C     METHOD.
C     -------
C     READS ARRAYS FROM FILE. FORTRAN READS OR PBIO READ ARE
C     USED DEPENDING ON THE CHOICE OF OUTPUT MODE AS SELECTED BY LPBIOOUT
C
C     IF PBIO IS USED, THEN PBREAD IS CALLED TO READ EACH ARRAY.
C     IF FORTRAN READ IS USED THEN FILNAME IS READ AS UNFORMATTED BINARY 
C     ONLY THE NECESSARY CONTRIBUTION FROM FL IS KEPT ON PE IN MESSAGE.
C
C     EXTERNALS.
C     ----------
C     *PBOPEN*
C     *PBREAD*
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
C
C ----------------------------------------------------------------------
      REAL,DIMENSION(NIBLO,NBLO) :: RU10OLD,RTHWOLD,RUSOLD,RZ0OLD,RTAUW
C ----------------------------------------------------------------------

      IF(LPBIOIN) THEN

        MODE = 'r'
        CALL PBOPEN(IUNIT,FILENAME,MODE,KRET)
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBOPEN       *'
          WRITE (IU06,*) '*   IN READSTRESS                        *' 
          IF(KRET.EQ.-1)
     1      WRITE (IU06,*) 'COULD NOT OPEN FILE ',FILENAME
          IF(KRET.EQ.-2) WRITE (IU06,*) 'INVALID FILENAME ',FILENAME
          IF(KRET.EQ.-3) WRITE (IU06,*) 'INVALID OPEN MODE SPECIFIED'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          CALL ABORT1
        ENDIF


        KOUNT = 10 
        CALL PBREAD(IUNIT,CDTPRO,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

        CALL PBREAD(IUNIT,CDATEWO,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

        CALL PBREAD(IUNIT,CDAWIFL,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

        CALL PBREAD(IUNIT,CDATEFL,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

        KOUNT = NIBLO*NBLO*NPRECI
        CALL PBREAD(IUNIT,RU10OLD,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

        CALL PBREAD(IUNIT,RTHWOLD,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

        CALL PBREAD(IUNIT,RUSOLD,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

        CALL PBREAD(IUNIT,RTAUW,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

        CALL PBREAD(IUNIT,RZ0OLD,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

1099    CONTINUE
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '*****************************************'
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*    ERROR FOLLOWING CALL TO PBREAD     *'
          WRITE (IU06,*) '*    IN READSTRESS                      *' 
          WRITE (IU06,*) '*    FILE ',FILENAME 
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*    KRET = ', KRET 
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*****************************************'
          CALL ABORT1
        ELSE IF(KRET.NE.KOUNT) THEN
          WRITE (IU06,*) '*****************************************'
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*    ERROR FOLLOWING CALL TO PBREAD     *'
          WRITE (IU06,*) '*    IN READSTRESS                      *' 
          WRITE (IU06,*) '*    FILE ',FILENAME 
          WRITE (IU06,*) '*    NUMBER OF BYTES READ : ',KRET
          WRITE (IU06,*) '*    IS NOT EQUAL TO THE NUMBER OF BYTES*'
          WRITE (IU06,*) '*    REQUESTED : ',KOUNT 
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*****************************************'
          CALL ABORT1
        ENDIF


        CALL PBCLOSE(IUNIT,KRET)
        IF(KRET.LT.0) THEN
           WRITE (IU06,*) '************************************'
           WRITE (IU06,*) '*                                  *'
           WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
           WRITE (IU06,*) '* IN READSTRESS                    *'
           WRITE (IU06,*) '* FILE ',FILENAME
           WRITE (IU06,*) '*                                  *'
           WRITE (IU06,*) '************************************'
           CALL ABORT1
        ENDIF

      ELSE

        READ(IU15) CDTPRO,CDATEWO,CDAWIFL,CDATEFL
        READ(IU15) RU10OLD
        READ(IU15) RTHWOLD
        READ(IU15) RUSOLD
        READ(IU15) RTAUW
        READ(IU15) RZ0OLD

        CLOSE(IU15)


      ENDIF

      IF (ITEST.GE.2)
     1    WRITE(IU06,*) ' SUB. READSTRESS: INPUT OF RESTART FILES DONE'
C
      RETURN
      END
