C THE PURPOSE OF THIS DUMMY SUBROUTINE PACKAGE IS TO SUPPLEMEMT
C THE FUJITSU COMPILER WITH A FICTIVE DEFINITION OF SUBROUTINES
C FOUND IN THE WAM CODE (IE: THERE IS A CALL TO IT) BUT THEIR 
C REFERENCES TO EXTERNALS CANNOT BE RESOLVED ON THE FUJITSU.
C
C BJORN HANSEN  ECMWF DECEMBER 1997
C
C -----------------------------------------------------------------
C
       SUBROUTINE dummy_alt()
C
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
C
       ENTRY psetc
       ENTRY pseti
       ENTRY psetr
C
       ENTRY pset1c
       ENTRY pset1i
       ENTRY pset1r
C
       ENTRY pset2i
       ENTRY pset2r
C
       ENTRY pset3i
       ENTRY pset3r
C
       ENTRY penqc
       ENTRY penqi
       ENTRY penqr
C
       ENTRY penq1c
       ENTRY penq1i
       ENTRY penq1r
C
       ENTRY penq2i
       ENTRY penq2r
C
       ENTRY penq3i
       ENTRY penq3r
C
       ENTRY psymb
       ENTRY playout

C
       WRITE(0,*) ' ************************************************** '
       WRITE(0,*) ' *                                                * '
       WRITE(0,*) ' *     SUB: DUMMY_ALT                             * '
       WRITE(0,*) ' *     REFERENCE TO UNRESOLVED EXTERNAL MADE      * '
       WRITE(0,*) ' *                                                * '
       WRITE(0,*) ' * PROGRAM ABORTS  PROGRAM ABORTS  PROGRAM ABORTS * '
       WRITE(0,*) ' *                                                * '
       WRITE(0,*) ' ************************************************** '
       CALL ABORT
       END
