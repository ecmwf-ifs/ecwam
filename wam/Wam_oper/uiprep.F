      SUBROUTINE UIPREP (IU05, IFORM, IREFRA, ML, KL)

! ----------------------------------------------------------------------

!**** *UIPREP* - ROUTINE TO READ NAMELIST INPUT FOR PREPROC.

!     H.GUNTHER            ECMWF       04/04/1990

!*    PURPOSE.
!     -------

!       TO READ USER INPUT OF PROGRAM PREPROC AND CHECKS CONSISTENCY.

!**   INTERFACE.
!     ----------

!       *CALL* *UIPREP (IU05, IFORM, IREFRA, ML, KL)*
!          *IU05*    - LOGICAL INPUT UNIT.
!          *IFORM*   - OUTPUT FORMAT OPTION = 1 UNFORMATED
!                                           = 2 FORMATED
!                                           OTHERWISE BOTH
!          *IREFRA*  - REFRACTION OPTION = 2 CURRENT PROCESSING.
!                                        OTHERWISE NOTHING.
!          *ML*      - NUMBER OF FREQUENCIES.
!          *KL*      - NUMBER OF DIRECTIONS.

!     METHOD.
!     -------
!         NAMELIST READ.

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.
!       *ADJUST*    - CORRECTS LONGITUDE INPUT.

!     REFERENCE.
!     ----------

!       NONE.

!     MODIFICATIONS.                                                    
!     --------------                                                    

!        B. HANSEN    *ECMWF*      JAN 98
!           REPLACE FORMATED READ OF USER INPUT BY NAMELIST READ.

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      ,
     &            NBLO     ,NIBLO    ,CLDOMAIN
      USE YOWCPBO  , ONLY : IBOUNC   ,AMOSOC   ,AMONOC   ,AMOEAC   ,
     &            AMOWEC
      USE YOWCINP  , ONLY : NOUT     ,XOUTW    ,XOUTS    ,XOUTE    ,
     &            XOUTN    ,NOUTD    ,OUTLONG  ,OUTLAT
      USE YOWCOUT  , ONLY : NGOUT
      USE YOWFPBO  , ONLY : IBOUNF
      USE YOWFRED  , ONLY : FR       ,FRATIO
      USE YOWMAP   , ONLY : NX       ,NY       ,IPER     ,IRGG     ,
     &            AMOWEP   ,AMOSOP   ,AMOEAP   ,AMONOP   ,XDELLA   ,
     &            XDELLO
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB

! ----------------------------------------------------------------------

      CHARACTER* 70 CLINE
      REAL, ALLOCATABLE :: XDUMP(:)
      INTEGER, ALLOCATABLE :: NDUMP(:)
      NAMELIST /NALINE/ CLINE, ML, KL, FR1, IRGG, XDELLA, XDELLO,
     &                  AMOSOP, AMONOP, AMOWEP, AMOEAP,
     &                  IFORM, IREFRA, ITEST, ITESTB,
     &                  IBOUNC, IBOUNF, AMOSOC, AMONOC, AMOWEC, AMOEAC,
     &                  NBLO, NIBLO, CLDOMAIN
      NAMELIST /NACORR/ ZOUTS, ZOUTN, ZOUTW, ZOUTE, IOUTD
      NAMELIST /NAOUTP/ ZOUTLAT, ZOUTLONG 


! ----------------------------------------------------------------------

!*    0. SET DEFAULT VALUES FOR THE NAMELIST ELEMENTS.
!        ---------------------------------------------

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

!     THE FOLLOWING DEFAULT VALUES INDICATES THAT THEY ARE DETERMINED
!     INTERNALLY BY PREPROC EXCEPT IF PROVIDED IN THE NAMELIST
      NBLO   = -1
      NIBLO  = -1 ! 
      CLDOMAIN= '-'

! ----------------------------------------------------------------------

!*    1. READ NAMELIST NALINE.
!        ---------------------

      READ (IU05, NALINE)

      ALLOCATE(FR(ML))

      IF(ML.EQ.30) THEN
        FR(1) = (FR1/FRATIO)/FRATIO
        FR1 = FR(1)
      ELSE
        FR(1) = FR1
      ENDIF

      WRITE (IU06,'(2X, A70,/)') CLINE(1:70)

! ----------------------------------------------------------------------


!*    2. PRINT CONTENTS OF NAMELIST AND PERFORM CONSISTENCY CHECKS.
!        ----------------------------------------------------------

!*    2.1 FREQUENCY AND DIRECTION GRID DEFINITIONS.

      WRITE (IU06,'("   FREQUENCY / DIRECTION GRID"/)')
      WRITE (IU06,'("   NUMBER OF FREQUENCIES IS ML = ",I6)') ML
      IF(ML.EQ.30) THEN
        WRITE (IU06,'("!!SINCE THE NUMBER OF FREQUENCIES IS",I6)')ML
        WRITE (IU06,'("!!MINIMUM FREQUENCY WAS RESET TO",F10.6)') FR(1)
      ELSE
        WRITE (IU06,'("   MINIMUM FREQUENCY IS  FR(1) = ",F10.6)') FR(1)
      ENDIF
      WRITE (IU06,'("   NUMBER OF DIRECTIONS  IS KL = ",I6)') KL

!     SET DIMENSIONS.

      NFRE = ML
      NANG = KL

! ----------------------------------------------------------------------

!*    2.2 OUTPUT GRID DEFINITIONS.

      WRITE (IU06,'(/"  OUTPUT GRID"/)')
      IF (IRGG.EQ.1) THEN
        WRITE(IU06,'("   USE A REDUCED GAUSSIAN GRID")')
      ELSE
        WRITE(IU06,'("   USE A REGULAR LON LAT GRID")')
      ENDIF      


      CALL ADJUST (AMOWEP, AMOEAP)
      IPER = 0
      IF (MOD(AMOEAP-AMOWEP+XDELLO+720., 360.).EQ.0.) IPER = 1

      WRITE (IU06,'("   RESOLUTION LAT-LON ",2F8.3)') XDELLA, XDELLO
      WRITE (IU06,'("    SOUTHERN LAT  NORTHERN LAT ",
     &  " WESTERN LONG "," EASTERN LONG",
     &  /,2X,4F14.3)') AMOSOP, AMONOP, AMOWEP, AMOEAP
      IF (IPER.EQ.1) WRITE (IU06,*) '   THE GRID IS EAST-WEST PERIODIC'

      IF(CLDOMAIN == '-' ) THEN
        IF (IPER.EQ.1) THEN
          CLDOMAIN = 'g'
        ELSE
          CLDOMAIN = 'm'
        ENDIF 
      ENDIF 

!*    SET DIMENSIONS.

      NX = NINT((AMOEAP-AMOWEP)/XDELLO) + 1
      NY = NINT((AMONOP-AMOSOP)/XDELLA) + 1

      NGX = NX
      NGY = NY

! ----------------------------------------------------------------------

!*    2.3 OUTPUT GRID CORRECTIONS READ FROM NAMELIST NACORR.
!         --------------------------------------------------

      REWIND (IU05)
      NOUT=0
      IOS=-9
      IOUTA=500
      ALLOCATE(XOUTW(IOUTA))
      ALLOCATE(XOUTS(IOUTA))
      ALLOCATE(XOUTE(IOUTA))
      ALLOCATE(XOUTN(IOUTA))
      ALLOCATE(NOUTD(IOUTA))
      XOUTS  =-100.0
      XOUTN  =-100.0
      XOUTW  =   0.0
      XOUTE  =   0.0
      NOUTD  =   0
      CORR: DO
        READ(IU05, NACORR, ERR=3000, IOSTAT=IOS, END=2300)
        NOUT=NOUT+1
        IF (NOUT.GT.IOUTA) THEN
          ALLOCATE(XDUMP(IOUTA))
          IOUTANEW=IOUTA+100
          WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++'
          WRITE (IU06,*) '+                                      +'
          WRITE (IU06,*) '+     WARNING IN SUB. UIPREP           +'
          WRITE (IU06,*) '+     =============================    +'
          WRITE (IU06,*) '+                                      +'
          WRITE (IU06,*) '+  NUMBER OF AREAS TO BE CORRECTED     +'
          WRITE (IU06,*) '+  EXCEEDS  DIMENSION   IOUTA = ', IOUTA
          WRITE (IU06,*) '+  THE DIMEMSION WAS RESET TO ',IOUTANEW
          WRITE (IU06,*) '+                                      +'
          WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++'

          XDUMP=XOUTW
          DEALLOCATE(XOUTW)
          ALLOCATE(XOUTW(IOUTANEW))
          DO IDUM=1,IOUTA
            XOUTW(IDUM)=XDUMP(IDUM)
          ENDDO

          XDUMP=XOUTS
          DEALLOCATE(XOUTS)
          ALLOCATE(XOUTS(IOUTANEW))
          DO IDUM=1,IOUTA
            XOUTS(IDUM)=XDUMP(IDUM)
          ENDDO

          XDUMP=XOUTE
          DEALLOCATE(XOUTE)
          ALLOCATE(XOUTE(IOUTANEW))
          DO IDUM=1,IOUTA
            XOUTE(IDUM)=XDUMP(IDUM)
          ENDDO

          XDUMP=XOUTN
          DEALLOCATE(XOUTN)
          ALLOCATE(XOUTN(IOUTANEW))
          DO IDUM=1,IOUTA
            XOUTN(IDUM)=XDUMP(IDUM)
          ENDDO

          DEALLOCATE(XDUMP)

          ALLOCATE(NDUMP(IOUTA))
          NDUMP=NOUTD
          DEALLOCATE(NOUTD)
          ALLOCATE(NOUTD(IOUTANEW))
          DO IDUM=1,IOUTA
            NOUTD(IDUM)=NDUMP(IDUM)
          ENDDO
          DEALLOCATE(NDUMP)

          IOUTA=IOUTANEW

          CALL ADJUST (ZOUTW, ZOUTE)
          XOUTS(NOUT) = ZOUTS
          XOUTN(NOUT) = ZOUTN
          XOUTW(NOUT) = ZOUTW
          XOUTE(NOUT) = ZOUTE
          NOUTD(NOUT) = IOUTD
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
     &   /,4X,"  NO.   SOUTHERN LAT ",
     &   " NORTHERN LAT  WESTERN LONG ",
     &   " EASTERN LONG  DEPTH")')
        DO I=1,NOUT
          CALL ADJUST (XOUTW(I), XOUTE(I))
          WRITE (IU06,'(4X,I5,1X,4F14.3,I7 )') I, XOUTS(I),
     &     XOUTN(I), XOUTW(I), XOUTE(I), NOUTD(I)
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    2.4 OUTPUT POINTS READ FROM NAMELIST NAOUTP.
!         ----------------------------------------

      REWIND (IU05)
      MOUTP = 100 
      ALLOCATE(OUTLAT(MOUTP))
      ALLOCATE(OUTLONG(MOUTP))
      OUTLAT =-100.0
      OUTLONG=   0.0
      NGOUT = 0
      OUTPP: DO
        READ(IU05, NAOUTP, ERR=3000, IOSTAT=IOS, END=2400)
        NGOUT = NGOUT + 1
        IF (NGOUT.GT.MOUTP) THEN
          ALLOCATE(XDUMP(MOUTP))
          MOUTPNEW=MOUTP+100
          WRITE (IU06,*) '+++++++++++++++++++++++++++++++++++++++++++'
          WRITE (IU06,*) '+                                         +'
          WRITE (IU06,*) '+     WARINING IN SUB. UIPREP       +'
          WRITE (IU06,*) '+     =============================       +'
          WRITE (IU06,*) '+                                         +'
          WRITE (IU06,*) '+  NUMBER OF OUTPUT POINTS EXCEEDS        +'
          WRITE (IU06,*) '+  EXCEEDS  DIMENSION      MOUTP = ', MOUTP
          WRITE (IU06,*) '+  THE DIMEMSION WAS RESET TO ', MOUTPNEW
          WRITE (IU06,*) '+                                         +'
          WRITE (IU06,*) '+++++++++++++++++++++++++++++++++++++++++++'

          XDUMP=OUTLAT
          DEALLOCATE(OUTLAT)
          ALLOCATE(OUTLAT(MOUTPNEW))
          DO IDUM=1,MOUTP
            OUTLAT(IDUM)=XDUMP(IDUM)
          ENDDO

          XDUMP=OUTLONG
          DEALLOCATE(OUTLONG)
          ALLOCATE(OUTLONG(MOUTPNEW))
          DO IDUM=1,MOUTP
            OUTLONG(IDUM)=XDUMP(IDUM)
          ENDDO

          DEALLOCATE(XDUMP)
          MOUTP=MOUTPNEW
          OUTLAT(NGOUT) = ZOUTLAT
          OUTLONG(NGOUT) = ZOUTLONG
        ELSE
          OUTLAT(NGOUT) = ZOUTLAT
          OUTLONG(NGOUT) = ZOUTLONG
        ENDIF
      ENDDO OUTPP
 2400 IF (NGOUT.GT.0) THEN
        WRITE (IU06,'(" OUTPUT POINTS FOR SPECTRA AS DEFINED",
     &   " BY USER INPUT",/,
     &   "     NO.    LAT.   LONG.")')
        DO I=1,NGOUT
          WRITE (IU06,'(3X,I5,2F8.2)') I,OUTLAT(I),OUTLONG(I)
        ENDDO

        ALLOCATE(XDUMP(NGOUT))
        MOUTPNEW=NGOUT
        DO IDUM=1,MOUTPNEW
          XDUMP(IDUM)=OUTLAT(IDUM)
        ENDDO
        DEALLOCATE(OUTLAT)
        ALLOCATE(OUTLAT(MOUTPNEW))
        DO IDUM=1,MOUTPNEW
          OUTLAT(IDUM)=XDUMP(IDUM)
        ENDDO

        DO IDUM=1,MOUTPNEW
          XDUMP(IDUM)=OUTLONG(IDUM)
        ENDDO
        DEALLOCATE(OUTLONG)
        ALLOCATE(OUTLONG(MOUTPNEW))
        DO IDUM=1,MOUTPNEW
          OUTLONG(IDUM)=XDUMP(IDUM)
        ENDDO
        DEALLOCATE(XDUMP)
      ENDIF

! ----------------------------------------------------------------------

!*    2.5 MODEL OPTIONS.
!         --------------

      WRITE (IU06,'(" OUTPUT OPTION IS       IFORM =",I3,
     &         " (1: UNFORM.  2: FORM.  3: BOTH)")') IFORM
      WRITE (IU06,'(" REFRACTION OPTION IS  IREFRA =",I3,
     &         " (2: CURRENTS  OTHERWISE: NOTHING)")') IREFRA
      WRITE (IU06,'(" TEST OUTPUT OPTION IS  ITEST =",I3,
     &         " (0: NO  >0: UP TO LEVEL ITEST)")') ITEST
      WRITE (IU06,'(" BLOCK TEST OPTION IS  ITESTB =",I3,
     &         " (0: NO  >0: UP TO BLOCK ITESTB)")') ITESTB

! ----------------------------------------------------------------------

!*    2.6 NESTED GRID INFORMATION.
!         ------------------------

!*    2.6.1 COARSE GRID OPTION.
!           -------------------

      WRITE (IU06,'(" COARSE GRID OPTION IS IBOUNC = ",I3)') IBOUNC
      IF (IBOUNC .EQ. 1) THEN
        CALL ADJUST (AMOWEC, AMOEAC)
        WRITE (IU06,*) '   THIS IS A COARSE GRID RUN  INFORMATION',
     &   ' FOR A FOLLOW UP FINE GRID WILL BE GENERATED'
        WRITE (IU06,'(/4X," NEST AREA IN COARSE GIRD IS",
     &   /,4X,"  SOUTHERN LAT  NORTHERN LAT ",
     &   " WESTERN LONG  EASTERN LONG ")')
        WRITE (IU06,'(4X,4F14.3)') AMOSOC, AMONOC, AMOWEC, AMOEAC

!*    2.6.1.1 ARE ALL CORNER POINTS OF THE NEST GRID POINTS?
!             ----------------------------------------------

        WEST = MOD(AMOWEC - AMOWEP + 720., 360.)
        EAST = MOD(AMOEAC - AMOWEP + 720., 360.)
        IF ((MOD(WEST, XDELLA) .NE. 0) .OR.
     &   (MOD(EAST, XDELLA) .NE. 0) .OR.
     &   (MOD(AMOSOC - AMOSOP, XDELLO) .NE. 0) .OR.
     &   (MOD(AMONOC - AMONOP, XDELLO) .NE. 0)) THEN
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

!*    2.6.1.2 INCLUDES THE COARSE GRID THE NEST GRID?
!             ---------------------------------------

        IF ((IPER.NE.1) .AND.
     &   (AMOWEP.GT.AMOWEC      .OR. AMOEAC.GT.AMOEAP     ) .AND.
     &   (AMOWEP.GT.AMOWEC+360. .OR. AMOEAC+360..GT.AMOEAP) .AND.
     &   (AMOWEP.GT.AMOWEC-360. .OR. AMOEAC-360..GT.AMOEAP) ) THEN
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
     &   AMOSOC .GE. AMONOC) THEN
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

!*    2.6.2 FINE GRID OPTION.
!           -----------------

      WRITE (IU06,'(" FINE GRID OPTION IS   IBOUNF = ",I3)') IBOUNF
      IF (IBOUNF .EQ. 1) THEN
        WRITE (IU06,*) '   THIS IS A FINE GRID RUN, INPUT FROM',
     &   ' A COARSE GRID IS EXPECTED'
      ELSE
        WRITE (IU06,*) '   BOUNDARY VALUES FROM A COARSE GRID',
     &   ' ARE NOT EXPECTED'
      ENDIF

      WRITE (IU06,'(" END OF USER INPUT PROTOCOL",/)')

      RETURN

! ----------------------------------------------------------------------


!*    3. ERROR HANDLING.
!        ---------------

 3000 CONTINUE
      WRITE (IU06,*) '**********************************************'
      WRITE (IU06,*) '*                                            *'
      WRITE (IU06,*) '*       FATAL ERROR IN SUB. UIPREP           *'
      WRITE (IU06,*) '*       ==========================           *'
      WRITE (IU06,*) '*                                            *'
      WRITE (IU06,*) '*  ERROR WHILE READING NAMELIST  STATUS=', IOS
      WRITE (IU06,*) '*                                            *'
      WRITE (IU06,*) '**********************************************'
      CALL ABORT1 
      END SUBROUTINE UIPREP
