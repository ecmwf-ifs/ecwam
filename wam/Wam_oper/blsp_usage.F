      SUBROUTINE BLSP_USAGE (ldhelp)

! ======================================================================

!   SEND MESSAGE TO INFORM ABOUT CORRECT USE OF PROGRAM
!                      AND ABORTS

      LOGICAL ldhelp

      IF ( ldhelp ) THEN
        WRITE(*,*)'   BLS2GRS converts a blsp file into a grsp file '
        WRITE(*,*)'           ore vice versa depending on which     '
        WRITE(*,*)'           file is an existing file.             '
        WRITE(*,*)'           Therefor one file should exist and    '
        WRITE(*,*)'           one file should not exist             '
        WRITE(*,*)' '
        WRITE(*,*)'   blsp2grs -b blsp file  -g grsp file  -C class '
        WRITE(*,*)'            -V version  -D domain  -f  -m iii -h '
        WRITE(*,*)' '
      ELSE
        WRITE(*,*)'   SELECTED OPTION IS INVALID                    '
        WRITE(*,*)'   CORRECT OPTIONS  AND ARGUMENTS ARE :          '
      ENDIF
      WRITE(*,*)' '
      WRITE(*,*)'     -b    filename of blsp file (required)      '
      WRITE(*,*)'     -g    filename of grsp file (required)      '
      WRITE(*,*)'     -m    iii  (iii is the model version number,'
      WRITE(*,*)'                 required)                       '
      WRITE(*,*)' '
      WRITE(*,*)'    Convention on filenames:                     '
      WRITE(*,*)' '
      WRITE(*,*)'    FIDYYMMDDHHFFFF   where:                     '
      WRITE(*,*)' '
      WRITE(*,*)'    FID   =  3 CHARACTER FILE ID                 '
      WRITE(*,*)'    YY    =  2 DIGIT YEAR INDICATOR              '
      WRITE(*,*)'    MM    =  2 DIGIT MONTH OF YEAR INDICATOR     '
      WRITE(*,*)'    DD    =  2 DIGIT DAY OF MONTH INDICATOR      '
      WRITE(*,*)'    HH    =  2 DIGIT HOUR OF DAY INDICATOR       '
      WRITE(*,*)'    FFFF  =  4 DIGIT FORECAST TIME IN HOURS      '
      WRITE(*,*)'             FROM BASE DATE-TIME-GROUP           '
      WRITE(*,*)' '
      WRITE(*,*)'   Further options are:                          '
      WRITE(*,*)'   -C class (default is RD)                      '
      WRITE(*,*)'   -V version (default is xxxx)                  '
      WRITE(*,*)'   -D domain ( medite or global, later is default'
      WRITE(*,*)'               "north", "south" or "both" as arguments'
      WRITE(*,*)'               will split global fields respectively.'
      WRITE(*,*)'   -f fdbaccess (if this option is provided 2d   ',
     &            ' spectra are written into fdb and not into file)'
      WRITE (*,*) '      -h help (prints this page) '
      IF ( .NOT. ldhelp ) CALL ABORT
      STOP ' helped? '
      END SUBROUTINE BLSP_USAGE


