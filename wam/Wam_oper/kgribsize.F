      SUBROUTINE KGRIBSIZE(IU06, IUIN, KLEN, ISIZE, CALLINGSUB) 

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
!     *CALL KGRIBSIZE(IU06, IUIN, KLEN, ISIZE, CALLINGSUB)*

!     METHOD.
!     -------
!     EXTERNALS.
!     ----------
!     REFERENCE.
!     ----------
! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER ::  IU06, IUIN, KRET, ISIZE, KLEN, KOFFSET
      CHARACTER*(*) :: CALLINGSUB

      WRITE(IU06,*) ' '
      WRITE(IU06,*) ' ***** WARNING IN ',CALLINGSUB
      WRITE(IU06,*) ' SIZE  OF KGRIB IS NOT BIG ENOUGH.'
      WRITE(IU06,*) ' IT WAS ', ISIZE

      KRET=0

      ISIZE=(KLEN+KIND(ISIZE)-1)/KIND(ISIZE)
      WRITE(IU06,*) ' IT SHOULD AT LEAST BE ', ISIZE
      WRITE(IU06,*) ' THE SIZE WAS RESET AUTOMATICALLY'
      WRITE(IU06,*) ' AND THE FIELD READ WITH THE NEW SIZE'
      WRITE(IU06,*) ' IF THIS PROBLEM OCCURS TOO OFTEN'
      WRITE(IU06,*) ' MODIFY THE VALUE OF NBIT IN SOURCE'
      WRITE(IU06,*) ' ***** WARNING ****** WARNING *****'
      WRITE(IU06,*) ' '
      CALL FLUSH(IU06)

!     RESET THE FILE POINTER TO READ FIELD AGAIN
      KOFFSET=-KLEN
      CALL PBSEEK(IUIN,KOFFSET,1,KRET)
      IF(KRET.EQ.-1) THEN
        WRITE(IU06,*) '***********************************'
        WRITE(IU06,*) '*  SUB. KGRIBSIZE : ' 
        WRITE(IU06,*) '*  PBSEEK : END OF FILE ENCOUNTED'
        WRITE(IU06,*) '***********************************'
        CALL ABORT1
      ENDIF
      IF(KRET.EQ.-2) THEN
        WRITE(IU06,*) '***********************************'
        WRITE(IU06,*) '*  SUB. KGRIBSIZE : ' 
        WRITE(IU06,*) '*  PBSEEK : FILE HANDLING ERROR'
        WRITE(IU06,*) '***********************************'
        CALL ABORT1
      ENDIF

      RETURN
      END SUBROUTINE KGRIBSIZE
