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
C       *CALL* *CREWFN (IUNIT, IDATE)*
C
C         *IUNIT*  - LOCIAL UNIT TO WHICH A FILE NAME HAS
C                    TO BE ASSIGNED.
C         *IDATE*  - DATE TIME GROUP FOR FILE NAME.
C
C     METHOD.
C     -------
C
C       FROM THE DATE TIME GROUP YYMMDDHHMM THE FILE NAME
C       WDDHHMM IS GENERATED AND ASSIGNED TO IUNIT.
!       IF A MPP COMPUTER IS USED THEN THE NAME BECOMES
!       WDDHHMMxxx WHERE xxx IS THE PROCESS NUMBER IN CHARACTER*3
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
      CHARACTER WDATE*10, WFN*10, CTAG*3

      ILEN=7
      WFN(1:1) = 'W'
      WFN(2:7) = WDATE(5:10)
      IF(LMESSPASS) THEN
        ILEN=10
        IRANK=MPE_MYRANK()
        WRITE(CTAG,'(I3.3)') IRANK
        WFN(8:10) = CTAG 
      ENDIF
      OPEN (UNIT=IUNIT, FILE=WFN(1:ILEN), FORM='UNFORMATTED',
     .      STATUS='UNKNOWN')

      RETURN
      END
