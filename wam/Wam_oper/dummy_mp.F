! THE PURPOSE OF THIS DUMMY SUBROUTINE PACKAGE IS TO SUPPLEMEMT
! THE FUJITSU COMPILER WITH A FICTIVE DEFINITION OF SUBROUTINES
! FOUND IN THE WAM CODE (IE: THERE IS A CALL TO IT) BUT THEIR 
! DEFINITION RELIES ON THE MESSAGE PASSING LIBRARY TO BE SUPPLEMENTED
! WHEN AVAILABLE
!
! BJORN HANSEN  ECMWF NOVEMBER 1996 
!
! -----------------------------------------------------------------
!
       SUBROUTINE dummy_mp()

       ENTRY mpl_abort
       ENTRY mpl_myrank
       ENTRY mpl_init
       ENTRY mpl_nproc
       ENTRY mpl_end
       ENTRY mpl_barrier
       ENTRY mpl_send
       ENTRY mpl_recv
       ENTRY mpl_broadcast
       ENTRY mpabort

       WRITE(*,*) ' NO MESSAGE PASSING'
       WRITE(*,*) ' NO MESSAGE PASSING'
       WRITE(*,*) ' NO MESSAGE PASSING'
       WRITE(*,*) ' NO MESSAGE PASSING'
       WRITE(*,*) ' NO MESSAGE PASSING'
       WRITE(*,*) ' NO MESSAGE PASSING'
       WRITE(*,*) ' NO MESSAGE PASSING'
       WRITE(*,*) ' NO MESSAGE PASSING'
       WRITE(*,*) ' NO MESSAGE PASSING'
       WRITE(*,*) ' NO MESSAGE PASSING'
       WRITE(*,*) ' NO MESSAGE PASSING'
       CALL ABORT
       END
