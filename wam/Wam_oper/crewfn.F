      SUBROUTINE CREWFN (IUNIT, WDATE)

C ----------------------------------------------------------------------
C
C**** *CREWFN* - CREATES A WIND FILE NAME.
C
C     H. GUNTHER       ECMWF      29/5/90.
!     J. BIDLOT        ECMWF      MARCH 1996  MESSAGE PASSING
C
C*    PURPOSE.
C     --------
C
C       TO CREATE A FILE NAME AND ASSIGN IT TO A UNIT.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *CREWFN (IUNIT, WDATE)*
C
C         *IUNIT*  - LOCIAL UNIT TO WHICH A FILE NAME HAS
C                    TO BE ASSIGNED.
C         *WDATE*  - DATE TIME GROUP FOR FILE NAME.
C
C     METHOD.
C     -------
C
C       FROM THE DATE TIME GROUP YYYYMMDDHHMM THE FILE NAME
C       WDDHHMM IS GENERATED AND ASSIGNED TO IUNIT.
!       IF A MPP COMPUTER IS USED THEN THE NAME BECOMES
!       WYYYYMMDDHHMMxxx WHERE xxx IS THE PROCESS NUMBER IN CHARACTER*3
!       FORMAT
C
C     EXTERNALS.
C     ---------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
#include "commesspass.h"

C
      CHARACTER WDATE*12, WFN*16, CTAG*3

      ILEN=13
      WFN(1:1) = 'W'
      WFN(2:13) = WDATE(1:12)
      IF(LMESSPASS) THEN
        ILEN=16
        IRANK=MPE_MYRANK()
        WRITE(CTAG,'(I3.3)') IRANK
        WFN(14:16) = CTAG 
      ENDIF
      OPEN (UNIT=IUNIT, FILE=WFN(1:ILEN), FORM='UNFORMATTED',
     .      STATUS='UNKNOWN')

      RETURN
      END
