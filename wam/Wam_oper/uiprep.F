      SUBROUTINE UIPREP (IU05, IFORM, IREFRA, ML, KL)

C ----------------------------------------------------------------------
C
C**** *UIPREP* - ROUTINE TO READ NAMELIST INPUT FOR PREPROC.
C
C     H.GUNTHER            ECMWF       04/04/1990
C
C*    PURPOSE.
C     -------
C
C       TO READ USER INPUT OF PROGRAM PREPROC AND CHECKS CONSISTENCY.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *UIPREP (IU05, IFORM, IREFRA, ML, KL, NX, NY)*
C          *IU05*    - LOGICAL INPUT UNIT.
C          *IFORM*   - OUTPUT FORMAT OPTION = 1 UNFORMATED
C                                           = 2 FORMATED
C                                           OTHERWISE BOTH
C          *IREFRA*  - REFRACTION OPTION = 2 CURRENT PROCESSING.
C                                        OTHERWISE NOTHING.
C          *ML*      - NUMBER OF FREQUENCIES.
C          *KL*      - NUMBER OF DIRECTIONS.
C
C     METHOD.
C     -------
C
C       FORMATED READ FROM LOGICAL UNIT IU05. A PRINTER PROTOCOL
C       IS DONE AND THE DATA ARE CHECKED FOR CONSISTENCY.
C
C     EXTERNALS.
C     ----------
C
C       *ABORT1*     - TERMINATES PROCESSING.
C       *ADJUST*    - CORRECTS LONGITUDE INPUT.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C     MODIFICATIONS.                                                    
C     --------------                                                    
C                                                                       
C        B. HANSEN    *ECMWF*      JAN 98
C           REPLACE FORMATED READ OF USER INPUT BY NAMELIST READ.
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "parameter_preproc.h"
C
#include "parameter_topo.h"
C
#include "comcpbo.h"
C
#include "comcinp.h"
C
#include "comcout.h"
C
#include "comfpbo.h"
C
#include "comfred.h"
C
#include "comgrid.h"
C
#include "commap.h"
C
#include "comtest.h"
C
C ----------------------------------------------------------------------
C
      CHARACTER* 70 CLINE
      NAMELIST /NALINE/ CLINE, ML, KL, FR1, IRGG, XDELLA, XDELLO,
     .                  AMOSOP, AMONOP, AMOWEP, AMOEAP,
     .                  IFORM, IREFRA, ITEST, ITESTB,
     .                  IBOUNC, IBOUNF, AMOSOC, AMONOC, AMOWEC, AMOEAC
      NAMELIST /NACORR/ ZOUTS, ZOUTN, ZOUTW, ZOUTE, IOUTD
      NAMELIST /NAOUTP/ ZOUTLAT, ZOUTLONG 
C
C
C ----------------------------------------------------------------------
C
C*    0. SET DEFAULT VALUES FOR THE NAMELIST ELEMENTS.
C        ---------------------------------------------
C
      CLINE  = ' '
      ML     =   0
      KL     =   0
      FR1    =   0.0
      IRGG   =  -1
      XDELLA =   0.0
      XDELLO =   0.0
      AMOSOP =-100.0
      AMONOP =-100.0
      AMOWEP =   0.0
      AMOEAP =   0.0
      IFORM  =  -1
      IREFRA =  -1
      ITEST  =  -1
      ITESTB =  -1
      IBOUNC =  -1
      IBOUNF =  -1
      AMOSOC =-100.0
      AMONOC =-100.0
      AMOWEC =   0.0
      AMOEAC =   0.0
      XOUTS  =-100.0
      XOUTN  =-100.0
      XOUTW  =   0.0
      XOUTE  =   0.0
      NOUTD  =   0
      OUTLAT =-100.0
      OUTLONG=   0.0
C
C ----------------------------------------------------------------------
C
C*    1. READ NAMELIST NALINE.
C        ---------------------
C
      READ (IU05, NALINE)
      FR(1) = FR1
C
      WRITE (IU06,'(" PROGRAM PREPROC:",//,"  USER INPUT:",/)')
      WRITE (IU06,'(2X, A70,/)') CLINE(1:70)
C
C ----------------------------------------------------------------------
C
C
C*    2. PRINT CONTENTS OF NAMELIST AND PERFORM CONSISTENCY CHECKS.
C        ----------------------------------------------------------
C
C*    2.1 FREQUENCY AND DIRECTION GRID DEFINITIONS.
C
      WRITE (IU06,'("   FREQUENCY / DIRECTION GRID"/)')
      WRITE (IU06,'("   NUMBER OF FREQUENCIES IS ML = ",I6)') ML
      WRITE (IU06,'("   MINIMUM FREQUENCY IS  FR(1) = ",F10.6)') FR(1)
      WRITE (IU06,'("   NUMBER OF DIRECTIONS  IS KL = ",I6)') KL
C
C     CHECK DIMENSION.
C
      IF (ML.GT.NFRE .OR. KL.GT.NANG) THEN
         WRITE (IU06,*) '**********************************************'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '*       FATAL ERROR IN SUB. UIPREP           *'
         WRITE (IU06,*) '*       ==========================           *'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '*  NUMBER OF FREQUENCIES OR DIRECTIONS       *'
         WRITE (IU06,*) '*  EXCEEDS DIMENSION.                        *'
         WRITE (IU06,*) '*  NUMBER OF FREQUENCIES IS      ML = ', ML
         WRITE (IU06,*) '*  DIMENSION FOR FREQUENIES IS NFRE = ', NFRE
         WRITE (IU06,*) '*  NUMBER OF DIRECTIONS  IS      KL = ', KL
         WRITE (IU06,*) '*  DIMENSION FOR DIRECTIONS IS NANG = ', NANG
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '**********************************************'
         CALL ABORT1
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    2.2 OUTPUT GRID DEFINITIONS.
C
      WRITE (IU06,'(/"  OUTPUT GRID"/)')
      IF (IRGG.EQ.1) THEN
         WRITE(IU06,'("   USE A REDUCED GAUSSIAN GRID")')
      ELSE
         WRITE(IU06,'("   USE A REGULAR LON LAT GRID")')
      END IF      
C
C      
      CALL ADJUST (AMOWEP, AMOEAP)
      IPER = 0
      IF (MOD(AMOEAP-AMOWEP+XDELLA+720., 360.).EQ.0.) IPER = 1

      WRITE (IU06,'("   RESOLUTION LAT-LON ",2F8.3)') XDELLA, XDELLO
      WRITE (IU06,'("    SOUTHERN LAT  NORTHERN LAT ",
     1  " WESTERN LONG "," EASTERN LONG",
     2  /,2X,4F14.3)') AMOSOP, AMONOP, AMOWEP, AMOEAP
      IF (IPER.EQ.1) WRITE (IU06,*) '   THE GRID IS EAST-WEST PERIODIC'
C
C*    CHECK DIMENSION.
C
      NX = NINT((AMOEAP-AMOWEP)/XDELLO) + 1
      NY = NINT((AMONOP-AMOSOP)/XDELLA) + 1
      IF (NX.GT.NGX .OR. NY.GT.NGY) THEN
         WRITE (IU06,*) '**********************************************'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '*       FATAL ERROR IN SUB. UIPREP           *'
         WRITE (IU06,*) '*       ==========================           *'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '*  NUMBER OF GRID POINTS IN OUTPUT GRID      *'
         WRITE (IU06,*) '*  EXCEEDS DIMENSION.                        *'
         WRITE (IU06,*) '*  NUMBER OF LONGITUDES IS       NX = ', NX
         WRITE (IU06,*) '*  DIMENSION FOR LONGITUDES IS  NGX = ', NGX
         WRITE (IU06,*) '*  NUMBER OF LATITUDES IS        NY = ', NY
         WRITE (IU06,*) '*  DIMENSION FOR LONGITUDES IS  NGY = ', NGY
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '**********************************************'
         CALL ABORT1
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    2.3 OUTPUT GRID CORRECTIONS READ FROM NAMELIST NACORR.
C         --------------------------------------------------
C
      REWIND (IU05)
      NOUT=0
      CORR: DO
         READ(IU05, NACORR, END=2300)
         NOUT=NOUT+1
         IF (NOUT.GT.IOUTA) THEN
            WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++'
            WRITE (IU06,*) '+                                      +'
            WRITE (IU06,*) '+     WARINING ERROR IN SUB. UIPREP    +'
            WRITE (IU06,*) '+     =============================    +'
            WRITE (IU06,*) '+                                      +'
            WRITE (IU06,*) '+  NUMBER OF AREAS TO BE CORRECTED     +'
            WRITE (IU06,*) '+  EXCEEDS  DIMENSION   IOUTA = ', IOUTA
            WRITE (IU06,*) '+  THE FIRST IOUTA AREAS ARE ONLY USED.+'
            WRITE (IU06,*) '+                                      +'
            WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++'
            NOUT = NOUT - 1
            EXIT CORR
         ELSE
            CALL ADJUST (ZOUTW, ZOUTE)
            XOUTS(NOUT) = ZOUTS
            XOUTN(NOUT) = ZOUTN
            XOUTW(NOUT) = ZOUTW
            XOUTE(NOUT) = ZOUTE
            NOUTD(NOUT) = IOUTD
         ENDIF
      ENDDO CORR
 2300 IF (NOUT.GT.0) THEN
         WRITE (IU06,'(/4X," AREAS TO BE CORRECTED IN OUTPUT GRID",
     1                 /,4X,"  NO.   SOUTHERN LAT ",
     2                 " NORTHERN LAT  WESTERN LONG ",
     3                 " EASTERN LONG  DEPTH")')
         DO I=1,NOUT
            CALL ADJUST (XOUTW(I), XOUTE(I))
            WRITE (IU06,'(4X,I5,1X,4F14.3,I7 )') I, XOUTS(I),
     1               XOUTN(I), XOUTW(I), XOUTE(I), NOUTD(I)
         ENDDO
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    2.4 OUTPUT POINTS READ FROM NAMELIST NAOUTP.
C         ----------------------------------------
C
      REWIND (IU05)
      NGOUT = 0
      OUTPP: DO
         READ(IU05, NAOUTP, END=2400)
         NGOUT = NGOUT + 1
         IF (NGOUT.GT.MOUTP) THEN
            WRITE (IU06,*) '+++++++++++++++++++++++++++++++++++++++++++'
            WRITE (IU06,*) '+                                         +'
            WRITE (IU06,*) '+     WARINING ERROR IN SUB. UIPREP       +'
            WRITE (IU06,*) '+     =============================       +'
            WRITE (IU06,*) '+                                         +'
            WRITE (IU06,*) '+  NUMBER OF OUTPUT POINTS EXCEEDS        +'
            WRITE (IU06,*) '+  EXCEEDS  DIMENSION      MOUTP = ', MOUTP
            WRITE (IU06,*) '+  THE FIRST MOUTP POINTS ARE ONLY USED.  +'
            WRITE (IU06,*) '+                                         +'
            WRITE (IU06,*) '+++++++++++++++++++++++++++++++++++++++++++'
            NGOUT = NGOUT - 1
            EXIT OUTPP
         ELSE
            OUTLAT(NGOUT) = ZOUTLAT
            OUTLONG(NGOUT) = ZOUTLONG
         ENDIF
      ENDDO OUTPP
 2400 IF (NGOUT.GT.0) THEN
         WRITE (IU06,'(" OUTPUT POINTS FOR SPECTRA AS DEFINED",
     1                 " BY USER INPUT",/,
     2                 "     NO.    LAT.   LONG.")')
         DO I=1,NGOUT
            WRITE (IU06,'(3X,I5,2F8.2)') I,OUTLAT(I),OUTLONG(I)
         ENDDO
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    2.5 MODEL OPTIONS.
C         --------------
C
      WRITE (IU06,'(" OUTPUT OPTION IS       IFORM =",I3,
     1         " (1: UNFORM.  2: FORM.  3: BOTH)")') IFORM
      WRITE (IU06,'(" REFRACTION OPTION IS  IREFRA =",I3,
     1         " (2: CURRENTS  OTHERWISE: NOTHING)")') IREFRA
      WRITE (IU06,'(" TEST OUTPUT OPTION IS  ITEST =",I3,
     1         " (0: NO  >0: UP TO LEVEL ITEST)")') ITEST
      WRITE (IU06,'(" BLOCK TEST OPTION IS  ITESTB =",I3,
     1         " (0: NO  >0: UP TO BLOCK ITESTB)")') ITESTB
C
C ----------------------------------------------------------------------
C
C*    2.6 NESTED GRID INFORMATION.
C         ------------------------
C
C*    2.6.1 COARSE GRID OPTION.
C           -------------------
C
      WRITE (IU06,'(" COARSE GRID OPTION IS IBOUNC = ",I3)') IBOUNC
      IF (IBOUNC .EQ. 1) THEN
         CALL ADJUST (AMOWEC, AMOEAC)
         WRITE (IU06,*) '   THIS IS A COARSE GRID RUN  INFORMATION',
     1                  ' FOR A FOLLOW UP FINE GRID WILL BE GENERATED'
         WRITE (IU06,'(/4X," NEST AREA IN COARSE GIRD IS",
     1               /,4X,"  SOUTHERN LAT  NORTHERN LAT ",
     2                    " WESTERN LONG  EASTERN LONG ")')
         WRITE (IU06,'(4X,4F14.3)') AMOSOC, AMONOC, AMOWEC, AMOEAC
C
C*    2.6.1.1 ARE ALL CORNER POINTS OF THE NEST GRID POINTS?
C             ----------------------------------------------
C
         WEST = MOD(AMOWEC - AMOWEP + 720., 360.)
         EAST = MOD(AMOEAC - AMOWEP + 720., 360.)
         IF ((MOD(WEST, XDELLA) .NE. 0) .OR.
     1       (MOD(EAST, XDELLA) .NE. 0) .OR.
     2       (MOD(AMOSOC - AMOSOP, XDELLO) .NE. 0) .OR.
     3       (MOD(AMONOC - AMONOP, XDELLO) .NE. 0)) THEN
            WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++++'
            WRITE (IU06,*) '+                                        +'
            WRITE (IU06,*) '+    WARNING ERROR IN SUB. UIPREP        +'
            WRITE (IU06,*) '+    ============================        +'
            WRITE (IU06,*) '+ ERROR IN NEST SPECIFICATIONS.          +'
            WRITE (IU06,*) '+ ONE OR MORE CORNER POINTS ARE NOT      +'
            WRITE (IU06,*) '+ COARSE GRID POINTS.                    +'
            WRITE (IU06,*) '+ NEST INFORMATION WILL NOT BE GENERATED +'
            WRITE (IU06,*) '+                                        +'
            WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++++'
            IBOUNC = 0
         ENDIF
C
C*    2.6.1.2 INCLUDES THE COURSE GRID THE NEST GRID?
C             ---------------------------------------
C
         IF ((IPER.NE.1) .AND.
     1       (AMOWEP.GT.AMOWEC      .OR. AMOEAC.GT.AMOEAP     ) .AND.
     2       (AMOWEP.GT.AMOWEC+360. .OR. AMOEAC+360..GT.AMOEAP) .AND.
     3       (AMOWEP.GT.AMOWEC-360. .OR. AMOEAC-360..GT.AMOEAP) ) THEN
            WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++++'
            WRITE (IU06,*) '+                                        +'
            WRITE (IU06,*) '+     WARNING ERROR IN SUB. UIPREP       +'
            WRITE (IU06,*) '+     ============================       +'
            WRITE (IU06,*) '+ ERROR IN NEST SPECIFICATIONS.          +'
            WRITE (IU06,*) '+ WEST OR EAST BOUNDARY IS NOT IN COARSE +'
            WRITE (IU06,*) '+ GRID AREA.                             +'
            WRITE (IU06,*) '+                                        +'
            WRITE (IU06,*) '+     AMOWEC AND/OR AMOEAC ARE WRONG     +'
            WRITE (IU06,*) '+                                        +'
            WRITE (IU06,*) '+ NEST INFORMATION WILL NOT BE GENERATED +'
            WRITE (IU06,*) '+                                        +'
            WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++++'
            IBOUNC = 0
         ENDIF
         IF (AMOSOP .GT. AMOSOC .OR. AMONOC .GT. AMONOP .OR.
     1       AMOSOC .GE. AMONOC) THEN
            WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++++'
            WRITE (IU06,*) '+                                        +'
            WRITE (IU06,*) '+     WARNING ERROR IN SUB. UIPREP       +'
            WRITE (IU06,*) '+     ============================       +'
            WRITE (IU06,*) '+ ERROR IN NEST SPECIFICATIONS.          +'
            WRITE (IU06,*) '+ NORTH OR SOUTH BOUNDARY IS NOT IN      +'
            WRITE (IU06,*) '+ COARSE GRID AREA, OR SOUTH IS GE NORTH +'
            WRITE (IU06,*) '+                                        +'
            WRITE (IU06,*) '+     AMOSOC AND/OR AMONOC ARE WRONG     +'
            WRITE (IU06,*) '+                                        +'
            WRITE (IU06,*) '+ NEST INFORMATION WILL NOT BE GENERATED +'
            WRITE (IU06,*) '+                                        +'
            WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++++'
            IBOUNC = 0
         ENDIF
      ELSE
         WRITE (IU06,*) '   A NEST IS NOT INCLUDED IN THIS GRID'
      ENDIF
C
C*    2.6.2 FINE GRID OPTION.
C           -----------------
C
      WRITE (IU06,'(" FINE GRID OPTION IS   IBOUNF = ",I3)') IBOUNF
      IF (IBOUNF .EQ. 1) THEN
         WRITE (IU06,*) '   THIS IS A FINE GRID RUN, INPUT FROM',
     1                  ' A COARSE GRID IS EXPECTED'
      ELSE
         WRITE (IU06,*) '   BOUNDARY VALUES FROM A COARSE GRID',
     1                  ' ARE NOT EXPECTED'
      ENDIF

      WRITE (IU06,'(" END OF USER INPUT PROTOCOL",/)')
C
      RETURN
      END
