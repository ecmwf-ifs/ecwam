C ----------------------------------------------------------------------
C
      SUBROUTINE FILE_TRANSFER (IU06, IUNIT, FILENA, CDOPT) 

C ----------------------------------------------------------------------
C
C**** *FILE_TRANSFER* - COPY FORTRAN UNIT FILE TO A FILE 
C
C     J. BIDLOT     ECMWF JUNE  1996  MESSAGE PASSING
C     B. HANSEN     ECMWF APRIL 1996  MOVE CYMPPATH TO COMTEXT AND ADD 
C                                     CDOPT.
C                                                   
C*    PURPOSE.
C     --------
C
C     COPY FORTRAN UNIT FILE TO A FILE 
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *FILE_TRANSFER*(IU06, IUNIT, FILENA) 
C
C         *IU06*   INTEGER    UNIT FOR PRINTER MESSAGES.
C         *IUNIT*  INTEGER    FORTRAN UNIT.
C         *FILENA* CHARACTER  DESTINATION FILENAME
C         *CDOPT*  CHARACTER  OPTION
C                             = 'S' SAVE FILE CONNECTED TO IUNIT
C                             = 'G' GET FILE AND CONNECT FILE TO IUNIT
C
C     EXTERNALS.
C     ----------
C
C       *PBOPEN*
C       *PBREAD2*
C       *PBWRITE*
C       *PBCLOSE* 
C       *IECF_LEN*  -  LENGTH OF A CHARACTER VARIABLE.
C       *EXPRAND_STRING - ADD EXTENSION TO FILENAME FUNCTION OF THE PE 
C
C     METHOD.
C     -------
C
C     REFERENCES.
C      -----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
C
#include "param.h"
C
!
#include "commesspass.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "comstat.h"
C
#include "comtext.h"
C
      INTEGER IU06, IUNIT
      CHARACTER*  1 MODE,
     1              CDOPT
      CHARACTER* 16 FILENA
      CHARACTER* 80 FILENAM,
     1              FILENAMEF, 
     2              FILENAMET 
C
      INTEGER, ALLOCATABLE :: KARRAY(:) 
      LOGICAL LMORE
C
C ----------------------------------------------------------------------
        LMORE=.TRUE.
C
        IF(LMESSPASS.AND.IUNIT.EQ.12.AND.(IREST.EQ.2.OR.IREST.EQ.4))THEN
          FILENAMEF='BLS_%n.%p'
          CALL EXPAND_STRING(IRANK,NPROC,0,0,FILENAMEF,1)
        ELSEIF(LMESSPASS.AND.IUNIT.EQ.15.AND.
     1         (IREST.EQ.2.OR.IREST.EQ.4)) THEN
          FILENAMEF='LAW_%n.%p'
          CALL EXPAND_STRING(IRANK,NPROC,0,0,FILENAMEF,1)
        ELSE
          IF (IUNIT.LT.10) THEN
            WRITE(FILENAMEF,'(A,I1)') 'fort.',IUNIT
          ELSE
            WRITE(FILENAMEF,'(A,I2)') 'fort.',IUNIT
          ENDIF
        ENDIF
C
        IF(LMESSPASS.AND.IUNIT.EQ.12.AND.(IREST.EQ.2.OR.IREST.EQ.4))THEN
          FILENAM=FILENA//'_%n.%p'
          CALL EXPAND_STRING(IRANK,NPROC,0,0,FILENAM,1)
        ELSEIF(LMESSPASS.AND.IUNIT.EQ.15.AND.
     1         (IREST.EQ.2.OR.IREST.EQ.4)) THEN
          FILENAM=FILENA//'_%n.%p'
          CALL EXPAND_STRING(IRANK,NPROC,0,0,FILENAM,1)
        ELSE
           FILENAM=FILENA
        ENDIF

        LNAME = IECF_LEN(FILENAM)

        IF(ICPLEN.EQ.0) THEN
          FILENAMET=FILENAM(1:LNAME)
        ELSE
          FILENAMET=CPATH(1:ICPLEN)//'/'//FILENAM(1:LNAME)
        ENDIF
C
        IF (CDOPT.EQ.'G') THEN
          FILENAM = FILENAMET
          FILENAMET = FILENAMEF
          FILENAMEF = FILENAM
        ENDIF
C
        WRITE(IU06,*) 'FILE_TRANSFER: COPY FILE ',
     1                 FILENAMEF(1:IECF_LEN(FILENAMEF)), ' TO ',
     2                 FILENAMET(1:IECF_LEN(FILENAMET))

C allocate KARRAY, could tune this in function of what file has to be
C moved !!
#if resolution==50
        KLENG=((NSUP-NINF+1)*NANG*NFRE*NPRECI)/9+8+1
#else
        KLENG=(NSUP-NINF+1)*NANG*NFRE*NPRECI+8+1
#endif
        ALLOCATE(KARRAY(KLENG)) 
C
C       1. OPEN INPUT AND OUTPUT FILES
C          ---------------------------
        MODE='r'
        CALL PBOPEN(KUNITF,FILENAMEF,MODE,KRET)
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBOPEN       *'
          WRITE (IU06,*) '*   IN FILE_TRANSFER                     *'
          IF(KRET.EQ.-1)
     1       WRITE (IU06,*) ' COULD NOT OPEN FILE ',FILENAMEF
          IF(KRET.EQ.-2) WRITE (IU06,*) ' INVALID FILENAME ',FILENAMEF
          IF(KRET.EQ.-3) WRITE (IU06,*) ' INVALID OPEN MODE SPECIFIED'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          CALL ABORT1
        ENDIF
C
        MODE='w'
        CALL PBOPEN(KUNITO,FILENAMET,MODE,KRET)
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBOPEN       *'
          WRITE (IU06,*) '*   IN FILE_TRANSFER                     *'
          IF(KRET.EQ.-1)
     1        WRITE (IU06,*) 'COULD NOT OPEN FILE ',FILENAMET
          IF(KRET.EQ.-2) WRITE (IU06,*) 'INVALID FILENAME ',FILENAMET
          IF(KRET.EQ.-3) WRITE (IU06,*) 'INVALID OPEN MODE SPECIFIED'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          CALL ABORT1
        ENDIF
C
C       2. COPY FILES
C          ----------
        DO WHILE (LMORE)
C
          CALL PBREAD2(KUNITF,KARRAY,KLENG,KRET)
          IF(KRET.LT.0) THEN
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBREAD2      *'
          WRITE (IU06,*) '*   IN FILE_TRANSFER                     *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          CALL ABORT1
          ENDIF
C
          IF(KRET.LT.KLENG) THEN 
            KOUNT=KRET 
            LMORE=.FALSE.
          ELSE
            KOUNT=KLENG 
            LMORE=.TRUE.
          ENDIF
C
          CALL PBWRITE(KUNITO,KARRAY,KOUNT,KRET)
          IF(KRET.LT.0) THEN
          WRITE (IU06,*) '*****************************************'
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*    ERROR FOLLOWING CALL TO PBWRITE    *'
          WRITE (IU06,*) '*    IN FILE_TRANSFER                   *'
          WRITE (IU06,*) '*      FILE ',FILENAMET
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*****************************************'
          CALL ABORT1
          ENDIF
C
        ENDDO
C
C       3. CLOSE FILES
C          -----------
        CALL PBCLOSE(KUNITF,KRET) 
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '************************************'
          WRITE (IU06,*) '*                                  *'
          WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
          WRITE (IU06,*) '* IN FILE_TRANSFER                 *'
          WRITE (IU06,*) '* FILE ',FILENAMEF
          WRITE (IU06,*) '*                                  *'
          WRITE (IU06,*) '************************************'
          CALL ABORT1
        ENDIF
C
        CALL PBCLOSE(KUNITO,KRET) 
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '*******************************************'
          WRITE (IU06,*) '*                                         *'
          WRITE (IU06,*) '*     ERROR FOLLOWING CALL TO PBCLOSE     *'
          WRITE (IU06,*) '*     IN FILE_TRANSFER                    *'
          WRITE (IU06,*) '*     FILE ',FILENAMET
          WRITE (IU06,*) '*                                         *'
          WRITE (IU06,*) '*******************************************'
          CALL ABORT1
        ENDIF
C
      DEALLOCATE(KARRAY)
C
      RETURN
      END
