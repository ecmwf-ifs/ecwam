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
     &      STATUS='UNKNOWN')

      RETURN
      END SUBROUTINE CREWFN
