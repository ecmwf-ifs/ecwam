C ======================================================================
C
      SUBROUTINE USERIN (kuso, kuin, kpar, kareas, pba, pwo, pwa, ppo,
     .                   cdreg, cdoption, cdmode, ldmesy)
C
C ======================================================================
C
C**** *USERIN*
C
C      B. HANSEN         ECMWF      SEPTEMBER    1991
C
C     PURPOSE.
C     --------
C     READS THE USERINPUT, EVALUATS IT AND SETS SEVERAL SWICHTES.
C
C**   INTERFACE.
C     ----------
C     I/   kuso     PRINTER OUTPUT UNIT.
C     I/   kuin     USER INPUT UNIT.
C     I/   kpar     DIMENSION OF ARRAYS pba, pwo, pwa AND ppo
C      /O  kareas   NUMBER OF AEREAS TO BE EXAMINED.
C      /O  pba      NORTH BORDER (ARRAY).
C      /O  pwo      WEST BORDER  (ARRAY).
C      /O  pwa      SOUTH BORDER (ARRAY).
C      /O  ppo      EAST BORDER  (ARRAY).
C      /O  cdreg    LINE 4 OF THE HEADLINES; THE REGION NAME FOR THAT
C                   AREA (ARRAY ).
C      /O  cdoption IF SET TO 'C' PLOT WILL BE IN COLOUR; ELSE B&W
C      /O  cdmode   IF SET TO 'C' PLOT IS DONE BY COUNTOURING; 
C                   ELSE SQARED BOXES
C
C**   LIBRARIES.
C     ----------
C          NONE.
C
C**   EXTERNALS.
C     ----------
C          ABORT      TERMINATION ROUTINE IN ERROR CASE.
C          INCDATE    INCREMENTS DATE TIME GROUP.
C
C     METHOD.
C     -------
C          NONE
C
C     REFERENCE.
C     ----------
C          NONE
C
C ---------------------------------------------------------------------
C
C
C*    0. DEFINITIONS.
C        ------------
      REAL      pwa(kpar),pwo(kpar),pba(kpar),ppo(kpar)
C
      CHARACTER*  1 cdoption, cdmode
      CHARACTER* 60 cdreg(kpar)
      CHARACTER*120 cluser
C
      LOGICAL ldmesy
C
C ----------------------------------------------------------------------
C
#include <outctrl.h>
      DATA csubna /"userin "/
C
C ---------------------------------------------------------------------
C
C*    1. INITIALIZATIONS.
C        ----------------
 1000 CONTINUE
C
      msub = msub + 1
      CALL  uraopc (csubna, msub, llsuvi, lltest)
      iregion=  0
C
C ---------------------------------------------------------------------
C
C*    2. LOOP OVER USER INPUT.
C        ---------------------
 2000 CONTINUE
C
      DO 2100 juser = 1,900
        READ(kuin,'(a112)',END=2200) cluser
C
        IF (cluser(1:2) .EQ. 'C ' .OR. cluser(1:2) .EQ. 'C-' .OR.
     .      cluser(1:5) .EQ. 'LSUVI' .OR. 
     .      cluser(1:5) .EQ. 'LTEST' .OR.
     .      cluser(1:5) .EQ. 'CSUVI' .OR.
     .      cluser(1:5) .EQ. 'CTEST' ) THEN
          CONTINUE
        ELSEIF (cluser(1:6) .EQ. 'COLOUR') THEN
          cdoption=cluser(13:13)
          WRITE(kuso,*) '  COLOUR CHOICE IS ',cdoption
          cdmode=cluser(22:22)
          WRITE(kuso,*) '  PLOT MODE IS ', cdmode
          IF (cluser(30:30) .EQ. 'N' ) ldmesy = .FALSE.
        ELSEIF (cluser(1:5) .EQ. 'AREAS') THEN
          READ  (cluser(8:9),'(i2)') kareas
          WRITE(kuso,*) ' NO. OF AREAS: ',kareas
          IF (kareas.GT.kpar) THEN
            WRITE(kuso,*) ' THERE IS ONLY SPACE FOR ', kpar,' REGIONS!'
            CALL abort
          ENDIF
        ELSEIF (cluser(1:6) .EQ. 'REGION') THEN
          iregion = iregion + 1
          READ (cluser(13:20),'(f8.3)')pwa(iregion)
          READ (cluser(23:30),'(f8.3)')pwo(iregion)
          READ (cluser(33:40),'(f8.3)')pba(iregion)
          READ (cluser(43:50),'(f8.3)')ppo(iregion)
          cdreg(iregion)=cluser(53:112)
          IF (lltest) THEN
            WRITE(kuso,*) ' CORNER COORDINATES OF AREA <',iregion,'>'
            WRITE(kuso,*) ' LOW-LEFT-LAT  : ', pwa(iregion)
            WRITE(kuso,*) ' LOW-LEFT-LON  : ', pwo(iregion)
            WRITE(kuso,*) ' UPP-RIGHT-LAT : ', pba(iregion)
            WRITE(kuso,*) ' UPP-RIGHT-LON : ', ppo(iregion)
            WRITE(kuso,*) ' REGION        : ', cdreg(iregion)
          ENDIF
        ELSEIF (cluser(1:3) .EQ. 'END') THEN
          IF (iregion .NE. kareas) THEN
            WRITE(kuso,*)' no. of areas and no. of given corner'//
     .                    ' coordinates do not match!!'
            WRITE(kuso,*) ' no. of areas: ',kareas
            WRITE(kuso,*) ' no. of corners: ',iregion
            WRITE(kuso,*) ' PROGRAM  TERMINATES NOW   '
            WRITE(kuso,*)'                                 '
            CALL abort
          ELSE
            WRITE(kuso,*)'                                 '
            msub = msub - 1
            RETURN
          ENDIF
        ELSE
          WRITE(kuso,*)' ',cluser
          WRITE(kuso,*)'  userinput not correct   '
          WRITE(kuso,*) ' PROGRAM  TERMINATES NOW   '
          CALL abort
        ENDIF
 2100 CONTINUE
 2200 WRITE(kuso,*)'  userinput not correct   '
      WRITE(kuso,*) ' PROGRAM  TERMINATES NOW   '
      WRITE(kuso,*)' ',cluser
      CALL abort
      END
