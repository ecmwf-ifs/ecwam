! THE PURPOSE OF THIS DUMMY SUBROUTINE PACKAGE IS TO SUPPLEMEMT
! THE COMPILER WITH A FICTIVE DEFINITION OF SUBROUTINES
! FOUND IN THE WAM CODE (IE: THERE IS A CALL TO IT) BUT THEIR 
! DEFINITION RELY ON DIFFERENT LIBRARIES THAT ARE IN USE AT ECMWF
! OR COULD BE REPLACED BY A SIMILAR ROUTINE.

! JEAN BIDLOT  ECMWF NOVEMBER 1995

! -----------------------------------------------------------------
!

       SUBROUTINE ABOR1(CDTEXT)
!      This is a substitute routine for ABOR1 that is normally found
!      in the IFS library !
       CHARACTER(LEN=*) :: CDTEXT 

       WRITE(*,*) CDTEXT

       CALL ABORT1

       RETURN
       END SUBROUTINE ABOR1

       SUBROUTINE ABOR1FL(CDFILE, KLINENUM, CDTEXT)
!      This is a substitute routine for ABOR1 that is normally found
!      in the IFS library !
       USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
       CHARACTER(LEN=*) :: CDFILE, CDTEXT
       INTEGER(KIND=JWIM), INTENT(IN) :: KLINENUM

       WRITE(*,*) CDFILE
       WRITE(*,*) CDTEXT

       CALL ABORT1

       RETURN
       END SUBROUTINE ABOR1FL

       SUBROUTINE IFSSIG
       WRITE(*,*) 'THE ROUTINE IFSSIG WAS REPLACED BY A DUMMY'
       WRITE(*,*) 'IFSSIG is specific to signal handling for runs'
       WRITE(*,*) 'at  ECMWF, it can be commented out.'
       RETURN
       END SUBROUTINE IFSSIG

       SUBROUTINE CHESIG 
       WRITE(*,*) 'THE ROUTINE CHESIG WAS REPLACED BY A DUMMY'
       WRITE(*,*) 'CHESIG is specific to signal handling for runs'
       WRITE(*,*) 'at  ECMWF, it can be commented out.'
       RETURN
       END SUBROUTINE CHESIG

       SUBROUTINE SIGMASTER 
       WRITE(*,*) 'THE ROUTINE SIGMASTER WAS REPLACED BY A DUMMY'
       WRITE(*,*) 'SIGMASTER is specific to signal handling for runs'
       WRITE(*,*) 'at  ECMWF, it can be commented out.'
       RETURN
       END SUBROUTINE SIGMASTER

       SUBROUTINE GSTATS(KNUM,KSWITCH)

!     *GSTATS*  - The original GSTATS gathers timing statistics
!                 it was replaced here by a dummy to limit the amount
!                 of routines that should otherwise be obtained from
!                 the ifsaux library.
       USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
       INTEGER(KIND=JWIM) :: KNUM
       INTEGER(KIND=JWIM) :: KSWITCH

       RETURN
       END SUBROUTINE GSTATS

       SUBROUTINE JFH_BIND 
       WRITE(*,*) 'THE ROUTINE JFH_BIND WAS REPLACED BY A DUMMY'
       WRITE(*,*) 'It is specific to signal handling for runs'
       WRITE(*,*) 'at  ECMWF, it can be commented out.'
       RETURN
       END SUBROUTINE JFH_BIND

       SUBROUTINE EC_BIND 
       WRITE(*,*) 'THE ROUTINE EC_BIND WAS REPLACED BY A DUMMY'
       WRITE(*,*) 'It is specific to signal handling for runs'
       WRITE(*,*) 'at  ECMWF, it can be commented out.'
       RETURN
       END SUBROUTINE EC_BIND
