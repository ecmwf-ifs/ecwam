      SUBROUTINE CREWFN (IUNIT, WDATE)

! ----------------------------------------------------------------------

!**** *CREWFN* - CREATES A WIND FILE NAME.

!     H. GUNTHER       ECMWF      29/5/90.
!     J. BIDLOT        ECMWF      MARCH 1996  MESSAGE PASSING

!*    PURPOSE.
!     --------

!       TO CREATE A FILE NAME AND ASSIGN IT TO A UNIT.

!**   INTERFACE.
!     ----------

!       *CALL* *CREWFN (IUNIT, WDATE)*

!         *IUNIT*  - LOCIAL UNIT TO WHICH A FILE NAME HAS
!                    TO BE ASSIGNED.
!         *WDATE*  - DATE TIME GROUP FOR FILE NAME.

!     METHOD.
!     -------

!       FROM THE DATE TIME GROUP YYYYMMDDHHMM THE FILE NAME
!       WDDHHMM IS GENERATED AND ASSIGNED TO IUNIT.
!       IF A MPP COMPUTER IS USED THEN THE NAME BECOMES
!       WYYYYMMDDHHMMxxx WHERE xxx IS THE PROCESS NUMBER IN CHARACTER*3
!       FORMAT

!     EXTERNALS.
!     ---------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------


      USE YOWMESPAS, ONLY : LMESSPASS
      USE MPL_MODULE


      CHARACTER WDATE*14, WFN*18, CTAG*3

      ILEN=15
      WFN(1:1) = 'W'
      WFN(2:15) = WDATE(1:14)
      IF(LMESSPASS) THEN
        ILEN=18
        IRANK=MPL_MYRANK()
        WRITE(CTAG,'(I3.3)') IRANK
        WFN(16:18) = CTAG 
      ENDIF
      OPEN (UNIT=IUNIT, FILE=WFN(1:ILEN), FORM='UNFORMATTED',
     &      STATUS='UNKNOWN')

      RETURN
      END SUBROUTINE CREWFN
