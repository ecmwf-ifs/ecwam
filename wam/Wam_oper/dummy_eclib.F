! THE PURPOSE OF THIS DUMMY SUBROUTINE PACKAGE IS TO SUPPLEMEMT
! THE COMPILER WITH A FICTIVE DEFINITION OF SUBROUTINES
! FOUND IN THE WAM CODE (IE: THERE IS A CALL TO IT) BUT THEIR 
! DEFINITION RELY ON DIFFERENT LIBRARIES THAT ARE IN USE AT ECMWF
! BUT ARE NOT NEEDED TO RUN WAM WITHOUT DATA ASSIMILATION
! OR COULD BE REPLACED BY A SIMILAR ROUTINE.

! JEAN BIDLOT  ECMWF NOVEMBER 1995

! -----------------------------------------------------------------
!
       SUBROUTINE WAMASSI 
       WRITE(*,*) 'FATAL ERROR : ATTEMPT TO USE THE DUMMY DEFINITION OF'
       WRITE(*,*) '              SUBROUTINE WAMASSI !!!'
       WRITE(*,*) 'THE ROUTINE WAS REPLACED BY A DUMMY'
       WRITE(*,*) 'DEFINITION AS IT SHOULD NOT BE USED IF NO DATA'
       WRITE(*,*) 'ASSIMILATION IS DONE !!!'
       CALL ABORT1
       RETURN
       END
       SUBROUTINE SARINVERT 
       WRITE(*,*) 'FATAL ERROR : ATTEMPT TO USE THE DUMMY DEFINITION OF'
       WRITE(*,*) '              SUBROUTINE SARINVERT !!!'
       WRITE(*,*) 'THE ROUTINE WAS REPLACED BY A DUMMY'
       WRITE(*,*) 'DEFINITION AS IT SHOULD NOT BE USED BECAUSE'
       WRITE(*,*) 'THE SAR INVERSION SOFTWARE WAS NOT PROVIDED !!!'
       CALL ABORT1
       RETURN
       END
       SUBROUTINE SYMINV 
       WRITE(*,*) 'FATAL ERROR : ATTEMPT TO USE THE DUMMY DEFINITION OF'
       WRITE(*,*) '              SUBROUTINE SYMINV !!!'
       WRITE(*,*) 'THE ROUTINE WAS REPLACED BY A DUMMY'
       WRITE(*,*) 'DEFINITION AS IT SHOULD NOT BE USED IF NO DATA'
       WRITE(*,*) 'ASSIMILATION IS DONE !!!'
       WRITE(*,*) 'THIS ROUTINE IS PART OF THE ECMWF LIBRARY :'
       WRITE(*,*) 'EMOSLIB OR ECLIB '
       CALL ABORT1
       RETURN
       END
       END
       SUBROUTINE SYMINV 
       WRITE(*,*) 'FATAL ERROR : ATTEMPT TO USE THE DUMMY DEFINITION OF'
       WRITE(*,*) '              SUBROUTINE SYMINV !!!'
       WRITE(*,*) 'THE ROUTINE WAS REPLACED BY A DUMMY'
       WRITE(*,*) 'DEFINITION AS IT SHOULD NOT BE USED IF NO DATA'
       WRITE(*,*) 'ASSIMILATION IS DONE !!!'
       WRITE(*,*) 'THIS ROUTINE IS PART OF THE ECMWF LIBRARY :'
       WRITE(*,*) 'EMOSLIB OR ECLIB '
       CALL ABORT1
       RETURN
       END

       SUBROUTINE DGELSS
       WRITE(*,*) 'FATAL ERROR : ATTEMPT TO USE THE DUMMY DEFINITION OF'
       WRITE(*,*) '              SUBROUTINE DGELSS !!!'
       WRITE(*,*) 'THE ROUTINE WAS REPLACED BY A DUMMY'
       WRITE(*,*) 'DEFINITION AS IT SHOULD NOT BE USED IF NO DATA'
       WRITE(*,*) 'ASSIMILATION IS DONE !!!'
       WRITE(*,*) 'THIS ROUTINE IS PART OF THE LAPACK LIBRARY'
       WRITE(*,*) 'WHICH IS NOT PROVIDED HERE !'
       CALL ABORT1
       RETURN
       END

       SUBROUTINE C06FUF
       WRITE(*,*) 'FATAL ERROR : ATTEMPT TO USE THE DUMMY DEFINITION OF'
       WRITE(*,*) '              SUBROUTINE C06FUF !!!'
       WRITE(*,*) 'THE ROUTINE WAS REPLACED BY A DUMMY'
       WRITE(*,*) 'DEFINITION AS IT SHOULD NOT BE USED IF NO DATA'
       WRITE(*,*) 'ASSIMILATION IS DONE !!!'
       WRITE(*,*) 'THIS ROUTINE IS PART OF THE NAG LIBRARY'
       WRITE(*,*) 'WHICH IS NOT PROVIDED HERE !'
       CALL ABORT1
       RETURN
       END

       SUBROUTINE ABOR1(CDTEXT)
!      This is a substitute routine for ABOR1 that is normally found
!      in the IFS library !
       CHARACTER(LEN=*) :: CDTEXT 

       WRITE(*,*) CDTEXT

       CALL ABORT1

       RETURN
       END

       SUBROUTINE IFSSIG
       WRITE(*,*) 'THE ROUTINE IFSSIG WAS REPLACED BY A DUMMY'
       WRITE(*,*) 'IFSSIG is specific to signal handling for runs'
       WRITE(*,*) 'at  ECMWF, it can be commented out.'
       RETURN
       END

       SUBROUTINE GSTATS(KNUM,KSWITCH)

!     *GSTATS*  - The original GSTATS gathers timing statistics
!                 it was replaced here by a dummy to limit the amount
!                 of routines that should otherwise be obtained from
!                 the ifsaux library.

       INTEGER :: KNUM
       INTEGER :: KSWITCH

       RETURN
       END

       SUBROUTINE JFH_BIND 
       WRITE(*,*) 'THE ROUTINE JFH_BIND WAS REPLACED BY A DUMMY'
       WRITE(*,*) 'IFSSIG is specific to signal handling for runs'
       WRITE(*,*) 'at  ECMWF, it can be commented out.'
       RETURN
       END

       SUBROUTINE C_DRHOOK_PRINT 
       WRITE(*,*) 'THE ROUTINE C_DRHOOK_PRINT WAS REPLACED BY A DUMMY'
       WRITE(*,*) 'IFSSIG is specific to signal handling for runs'
       WRITE(*,*) 'at  ECMWF, it can be commented out.'
       RETURN
       END
