      SUBROUTINE READFL(FL,MINF,MSUP,FILENAME,IUNIT,LPBOPEN,LPBCLOSE)

C ----------------------------------------------------------------------
C     J. BIDLOT    ECMWF      SEPTEMBER 1997
C
C*    PURPOSE.
C     --------
C     READS ARRAY FL
C
C**   INTERFACE.
C     ----------
C     CALL *READFL*(FL,MINF,MSUP,FILENAME,IUNIT,LPBOPEN,LPBCLOSE)
C     *FL*       ARRAY TO BE WRITTEN TO FILE
C     *MINF*     THIRD LOWER DIMENSION BOUND OF FL
C     *MSUP*     THIRD UPPER DIMENSION BOUND OF FL
C     note the first and second dimension of FL are given by 0:NIBLO and NANG
C     *FILENAME* FILENAME (INCLUDING PATH) OF INPUT FILE
C     *IUNIT*    PBIO UNIT (ONLY ACTIVE IF PBIO OUTPUT)
C     *LPBOPEN*  LOGICAL, TRUE IF PBOPEN HAS TO BE CALLED
C     *LPBCLOSE* LOGICAL, TRUE IF PBCLOSE HAS TO BE CALLED
C     *FILENAME* FILENAME (INCLUDING PATH) OF INPUT FILE
C
C     METHOD.
C     -------
C     READS ARRAY FL FROM FILE. A FORTRAN READ OR PBIO READ ARE
C     USED DEPENDING ON THE CHOICE OF OUTPUT MODE AS SELECTED BY LPBIOIN
C
C     IF PBIO IS USED, THEN PBREAD IS CALLED TO READ (MSUP-MINF+1) 
C     FREQUENCY CONTRIBUTION TO FL.
C     IF FORTRAN READ IS USED THEN FILNAME IS READ AS UNFORMATTED BINARY 
C     MINF SHOULD BE SET TO 1 AND MSUP TO NFRE
C  
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
#include "comgrid.h"
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
      REAL,DIMENSION(0:NIBLO,NANG,MINF:MSUP) :: FL
C ----------------------------------------------------------------------

      IF(LPBIOIN) THEN

        IF(LPBOPEN) THEN
          MODE = 'r'
          CALL PBOPEN(IUNIT,FILENAME,MODE,KRET)
          IF(KRET.LT.0) THEN
             WRITE (IU06,*) '******************************************'
             WRITE (IU06,*) '*                                        *'
             WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBOPEN       *'
             WRITE (IU06,*) '*   IN READFL                            *' 
             IF(KRET.EQ.-1)
     1           WRITE (IU06,*) 'COULD NOT OPEN FILE ',FILENAME
             IF(KRET.EQ.-2) WRITE (IU06,*) 'INVALID FILENAME ',FILENAME
             IF(KRET.EQ.-3) WRITE (IU06,*) 'INVALID OPEN MODE SPECIFIED'
             WRITE (IU06,*) '*                                        *'
             WRITE (IU06,*) '******************************************'
             CALL ABORT1
          ENDIF
        ENDIF

        KOUNT = (NIBLO+1)*NANG*(MSUP-MINF+1)*NPRECR
        CALL PBREAD(IUNIT,FL,KOUNT,KRET)
        IF(KRET.LT.0) THEN
           WRITE (IU06,*) '*****************************************'
           WRITE (IU06,*) '*                                       *'
           WRITE (IU06,*) '*    ERROR FOLLOWING CALL TO PBREAD     *'
           WRITE (IU06,*) '*    IN READFL                          *' 
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
           WRITE (IU06,*) '*    IN READFL                          *' 
           WRITE (IU06,*) '*    FILE ',FILENAME 
           WRITE (IU06,*) '*    NUMBER OF BYTES READ : ',KRET
           WRITE (IU06,*) '*    IS NOT EQUAL TO THE NUMBER OF BYTES*'
           WRITE (IU06,*) '*    REQUESTED : ',KOUNT 
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
              WRITE (IU06,*) '* IN READFL                        *'
              WRITE (IU06,*) '* FILE ',FILENAME
              WRITE (IU06,*) '*                                  *'
              WRITE (IU06,*) '************************************'
              CALL ABORT1
           ENDIF
        ENDIF


      ELSE

         IF(MINF.NE.1.OR.MSUP.NE.NFRE) THEN
           WRITE (IU06,*) '************************************'
           WRITE (IU06,*) '*                                  *'
           WRITE (IU06,*) '* ERROR IN READFL                  *' 
           WRITE (IU06,*) '* MINF MUST BE = 1                 *'
           WRITE (IU06,*) '* MSUP MUST BE = NFRE              *'
           WRITE (IU06,*) '* MINF WAS ', MINF
           WRITE (IU06,*) '* MSUP WAS ', MSUP
           WRITE (IU06,*) '* NFRE WAS ', NFRE 
           WRITE (IU06,*) '*                                  *'
           WRITE (IU06,*) '************************************'
           CALL ABORT1
         ENDIF

         READ(IU12) FL
         CLOSE(IU12)

      END IF

      IF (ITEST.GE.2)
     1            WRITE(IU06,*) ' SUB. READFL: INPUT OF FL DONE'
C
      RETURN
      END
