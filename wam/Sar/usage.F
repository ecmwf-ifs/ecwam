C=======================================================================
C
      SUBROUTINE USAGE
C
C=======================================================================
C
C   SEND MESSAGE TO INFORM ABOUT CORRECT USE OF PROGRAM
C                      AND ABORTS
C
      WRITE(*,*) '      SELECTED OPTION IS INVALID '
      WRITE(*,*) '      CORRECT OPTIONS  AND ARGUMENTS ARE :'
      WRITE(*,*) '      '
      WRITE(*,*) '        -i    filename of INV file (required) '
      WRITE(*,*) '        -b    filename of BNV file (required) '
      WRITE(*,*) '      '
      WRITE(*,*) '      ENDEINV bufferizes data of INV files and '
      WRITE(*,*) '              writes the buferd data to BNV files '
      WRITE(*,*) '              ore vice versa depending on which '
      WRITE(*,*) '              file is an existing file.  '
      WRITE(*,*) '              Therefor one file should exist and '
      WRITE(*,*) '              one file should not exist '
      WRITE(*,*) '      '
      CALL ABORT
      END
