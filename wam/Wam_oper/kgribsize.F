      SUBROUTINE KGRIBSIZE(IU06, KLEN, ISIZE, CALLINGSUB) 

! ----------------------------------------------------------------------

!**** *KGRIBSIZE* -

!     J. BIDLOT     ECMWF  JULY 2000 


!*    PURPOSE.
!     --------
!     UTILITY USED TO DETERMINE APOSTERIORI THE MINIMUM SIZE OF GRIB
!     INPUT DATA AND TO RESET THE FILE POINTER TO THE BEGINNING OF
!     THE FILE.

!**   INTERFACE.
!     ----------
!     *CALL KGRIBSIZE(IU06, KLEN, ISIZE, CALLINGSUB)*

!     METHOD.
!     -------
!     EXTERNALS.
!     ----------
!     REFERENCE.
!     ----------
! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER ::  IU06
      INTEGER(KIND=8) ::  KLEN
      INTEGER ::  ISIZE
      CHARACTER*(*) :: CALLINGSUB

      WRITE(IU06,*) ' '
      WRITE(IU06,*) ' ***** WARNING IN ',CALLINGSUB
      WRITE(IU06,*) ' SIZE  OF THE GRIB BUFFER IS NOT BIG ENOUGH.'
      WRITE(IU06,*) ' IT WAS ', ISIZE

      ISIZE=(KLEN+KIND(ISIZE)-1)/KIND(ISIZE)+1

      WRITE(IU06,*) ' IT SHOULD AT LEAST BE ', ISIZE
      WRITE(IU06,*) ' THE SIZE WAS RESET AUTOMATICALLY'
      WRITE(IU06,*) ' AND THE FIELD READ WITH THE NEW SIZE'
      WRITE(IU06,*) ' IF THIS PROBLEM OCCURS TOO OFTEN'
      WRITE(IU06,*) ' MODIFY THE VALUE OF NBIT IN SOURCE'
      WRITE(IU06,*) ' ***** WARNING ****** WARNING *****'
      WRITE(IU06,*) ' '
      CALL FLUSH(IU06)

      RETURN
      END SUBROUTINE KGRIBSIZE
