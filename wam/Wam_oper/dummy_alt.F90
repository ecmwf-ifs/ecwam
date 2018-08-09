! THE PURPOSE OF THIS DUMMY SUBROUTINE PACKAGE IS TO SUPPLEMEMT
! THE FUJITSU COMPILER WITH A FICTIVE DEFINITION OF SUBROUTINES
! FOUND IN THE WAM CODE (IE: THERE IS A CALL TO IT) BUT THEIR 
! REFERENCES TO EXTERNALS CANNOT BE RESOLVED ON THE FUJITSU.
!
! BJORN HANSEN  ECMWF DECEMBER 1997
!
! -----------------------------------------------------------------
!
       SUBROUTINE dummy_alt()
!
       ENTRY popen
       ENTRY pclose
       ENTRY pcoast
       ENTRY pcont
       ENTRY pwind
       ENTRY ptext
       ENTRY pops
       ENTRY pnew
       ENTRY pact
       ENTRY pgrib
       ENTRY preset
       ENTRY pset
!
       ENTRY psetc
       ENTRY pseti
       ENTRY psetr
!
       ENTRY pset1c
       ENTRY pset1i
       ENTRY pset1r
!
       ENTRY pset2i
       ENTRY pset2r
!
       ENTRY pset3i
       ENTRY pset3r
!
       ENTRY penqc
       ENTRY penqi
       ENTRY penqr
!
       ENTRY penq1c
       ENTRY penq1i
       ENTRY penq1r
!
       ENTRY penq2i
       ENTRY penq2r
!
       ENTRY penq3i
       ENTRY penq3r
!
       ENTRY psymb
       ENTRY playout

!
       WRITE(0,*) ' ************************************************** '
       WRITE(0,*) ' *                                                * '
       WRITE(0,*) ' *     SUB: DUMMY_ALT                             * '
       WRITE(0,*) ' *     REFERENCE TO UNRESOLVED EXTERNAL MADE      * '
       WRITE(0,*) ' *                                                * '
       WRITE(0,*) ' * PROGRAM ABORTS  PROGRAM ABORTS  PROGRAM ABORTS * '
       WRITE(0,*) ' *                                                * '
       WRITE(0,*) ' ************************************************** '
       CALL ABORT1
      END SUBROUTINE dummy_alt
