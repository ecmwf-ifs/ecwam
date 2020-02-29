      SUBROUTINE UIPREP (IFORM, ML, KL, LLGRID)

! ----------------------------------------------------------------------

!**** *UIPREP* - ROUTINE TO READ NAMELIST INPUT FOR PREPROC.

!     H.GUNTHER            ECMWF       04/04/1990

!*    PURPOSE.
!     -------

!       TO READ USER INPUT OF PROGRAM PREPROC AND CHECKS CONSISTENCY.

!**   INTERFACE.
!     ----------

!       *CALL* *UIPREP (IFORM, ML, KL, LLGRID)*
!          *IFORM*   - OUTPUT FORMAT OPTION = 1 UNFORMATED
!                                           = 2 FORMATED
!                                           OTHERWISE BOTH
!          *ML*      - NUMBER OF FREQUENCIES.
!          *KL*      - NUMBER OF DIRECTIONS.
!          *LLGRID*  - TRUE IF THE GRID DEFINITION HAS BEEN SPECIFIED
!                      IN INPUT FILE grid_description

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

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      ,    &
     &            NBLO     ,NIBLO    ,CLDOMAIN
      USE YOWCPBO  , ONLY : IBOUNC   ,GBOUNC_MAX, GBOUNC ,              &
     &            AMOSOC   ,AMONOC   ,AMOEAC   ,AMOWEC
      USE YOWCINP  , ONLY : NOUT     ,XOUTW    ,XOUTS    ,XOUTE    ,    &
     &            XOUTN    ,NOUTD    ,OUTLONG  ,OUTLAT
      USE YOWCOUT  , ONLY : NGOUT
      USE YOWFPBO  , ONLY : IBOUNF
      USE YOWFRED  , ONLY : FR       ,FRATIO
      USE YOWGRID  , ONLY : NLONRGG

      USE YOWMAP   , ONLY : NX       ,NY       ,IPER     ,IRGG     ,    &
     &            AMOWEP   ,AMOSOP   ,AMOEAP   ,AMONOP   ,              &
     &            XDELLA   ,XDELLO   ,LLOBSTRCT,LAQUA
      USE YOWSHAL  , ONLY : NDEPTH   ,DEPTHA   ,DEPTHD 
      USE YOWSTAT  , ONLY : IPHYS
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWUNPOOL, ONLY : LLUNSTR, LPREPROC, LVECTOR, IVECTOR

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "adjust.intfb.h"

      INTEGER(KIND=JWIM), INTENT(OUT) :: IFORM, ML, KL
      LOGICAL, INTENT(OUT) :: LLGRID

      INTEGER(KIND=JWIM) :: IWAM_GET_UNIT
      INTEGER(KIND=JWIM) :: K, M, I, II, KSN
      INTEGER(KIND=JWIM) :: IU05, IU
      INTEGER(KIND=JWIM) :: IFRE1, ISPECTRUNC
      INTEGER(KIND=JWIM) :: MOUTP, MOUTPNEW 
      INTEGER(KIND=JWIM) :: IOS, IOUTA, IOUTANEW, IDUM
      INTEGER(KIND=JWIM), ALLOCATABLE :: NDUMP(:)

      REAL(KIND=JWRB) :: FR1 
      REAL(KIND=JWRB) :: DEPTHMAX
      REAL(KIND=JWRB) :: ZOUTS, ZOUTN, ZOUTW, ZOUTE, IOUTD
      REAL(KIND=JWRB) :: ZOUTLAT, ZOUTLONG 
      REAL(KIND=JWRB) :: WEST, EAST, DW, DE, DS, DN
      REAL(KIND=JWRB), ALLOCATABLE :: XDUMP(:)

      CHARACTER(LEN=70) :: CLINE, FILENAME

! ----------------------------------------------------------------------

      NAMELIST /NALINE/ CLINE, ML, KL, FR1, IRGG, XDELLA, XDELLO,       &
     &                  AMOSOP, AMONOP, AMOWEP, AMOEAP,                 &
     &                  IFORM, ITEST, ITESTB,                           &
     &                  IPHYS,                                          &
     &                  IBOUNC, IBOUNF, AMOSOC, AMONOC, AMOWEC, AMOEAC, &
     &                  NBLO, NIBLO, CLDOMAIN,LLOBSTRCT,                &
     &                  NDEPTH   ,DEPTHA   ,DEPTHD,                     &
     &                  IFRE1, LAQUA, LLUNSTR, LPREPROC

      NAMELIST /NACORR/ ZOUTS, ZOUTN, ZOUTW, ZOUTE, IOUTD
      NAMELIST /NAOUTP/ ZOUTLAT, ZOUTLONG 


! ----------------------------------------------------------------------

!*    0. SET DEFAULT VALUES FOR THE NAMELIST ELEMENTS.
!        ---------------------------------------------

      CLINE  = ' '

      LLUNSTR = .FALSE. 
      LPREPROC = .FALSE. 
      LVECTOR   =.FALSE.
      IVECTOR   = 1
      ML     =   0
      KL     =   0
      FR1    =   0.0_JWRB
      IRGG   =  -1
      XDELLA =   0.0_JWRB
      XDELLO =   0.0_JWRB
      AMOSOP =-100.0_JWRB
      AMONOP =-100.0_JWRB
      AMOWEP =   0.0_JWRB
      AMOEAP =   0.0_JWRB
      IFORM  =  -1
      ITEST  =  -1
      ITESTB =  -1
!      IPHYS  =   0  ! ECMWF PHYSICS
      IPHYS  =   1  ! ARDHUIN PHYSICS
      IBOUNC =  -1
      IBOUNF =  -1
      AMOSOC =-100.0_JWRB
      AMONOC =-100.0_JWRB
      AMOWEC =   0.0_JWRB
      AMOEAC =   0.0_JWRB
      LLOBSTRCT = .TRUE.
      NDEPTH = 74
      DEPTHA = 1.0_JWRB
      DEPTHD = 1.1_JWRB
      IFRE1 = 1
      LAQUA =.FALSE. ! to force an aqua planet on irregular grid with xdella

!     THE FOLLOWING DEFAULT VALUES INDICATES THAT THEY ARE DETERMINED
!     INTERNALLY BY PREPROC EXCEPT IF PROVIDED IN THE NAMELIST
      NBLO   = -1
      NIBLO  = -1 ! 
      CLDOMAIN= '-'

! ----------------------------------------------------------------------

!*    1. READ NAMELIST NALINE.
!        ---------------------


      IU05 =  IWAM_GET_UNIT (IU06, 'procin', 'r', 'f', 0)

      READ (IU05, NALINE)

      ALLOCATE(FR(ML))

      FR(IFRE1) = FR1
      DO M=IFRE1-1,1,-1
        FR(M) = (FR(M+1)/FRATIO)
      ENDDO

      WRITE (IU06,'(2X, A70,/)') CLINE(1:70)


!     1.1 CHECK IF FILE grid_description IS PRESENT
!         IF IT IS THERE IT WILL SUPERSEDE THE NAMLIST DEFINITION

      FILENAME='grid_description'
      INQUIRE(FILE=FILENAME,EXIST=LLGRID)
      IF(LLGRID) THEN
        IU=IWAM_GET_UNIT(IU06,FILENAME,'S','F',0)
        OPEN(IU,FILE=FILENAME,STATUS='OLD', FORM='FORMATTED')
        READ (IU,*) ISPECTRUNC
        READ (IU,*) AMONOP
        READ (IU,*) AMOSOP
        READ (IU,*) AMOWEP
        READ (IU,*) AMOEAP
        READ (IU,*) IPER
        READ (IU,*) IRGG
        READ (IU,*) NY
        WRITE(IU06,*) "grid_description read in "
      ENDIF

! ----------------------------------------------------------------------


!*    2. PRINT CONTENTS OF NAMELIST AND PERFORM CONSISTENCY CHECKS.
!        ----------------------------------------------------------

!*    2.1 FREQUENCY AND DIRECTION GRID DEFINITIONS.

      WRITE (IU06,*) "   PHYSICS BASED ON IPHYS = ",IPHYS
      WRITE (IU06,'("   FREQUENCY / DIRECTION GRID"/)')
      WRITE (IU06,'("   NUMBER OF FREQUENCIES IS ML = ",I6)') ML
      IF(IFRE1.NE.1) THEN
        WRITE (IU06,'("!!MINIMUM FREQUENCY WAS RESET TO",F10.6)') FR(1)
      ELSE
        WRITE (IU06,'("  MINIMUM FREQUENCY IS  FR(1) = ",F10.6)') FR(1)
      ENDIF
      WRITE (IU06,'("   NUMBER OF DIRECTIONS  IS KL = ",I6)') KL

!     SET DIMENSIONS.

      NFRE = ML
      NANG = KL

      GBOUNC = 0
      DO II=1,GBOUNC_MAX
        IF(AMOSOC(II).EQ.-100.0_JWRB .AND. AMONOC(II).EQ.-100.0_JWRB .AND. &
     &    AMOWEC(II).EQ.0.0_JWRB .AND. AMOEAC(II).EQ.0.0_JWRB)EXIT
        GBOUNC=II
      ENDDO
      IF(IBOUNC .EQ. 1 .AND. GBOUNC.LE.0 ) THEN
        WRITE (IU06,*) '**********************************************'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '*       FATAL ERROR IN SUB. UIPREP           *'
        WRITE (IU06,*) '*       ==========================           *'
        WRITE (IU06,*) '* THE PRODUCTION OF FINE GRID BOUNDARY WAS   *'
        WRITE (IU06,*) '* REQUESTED BUT NO BOUNDING BOX VALUES WERE  *'
        WRITE (IU06,*) '* SUPPLIED IN THE NAMELIST!                  *'
        WRITE (IU06,*) '**********************************************'
        CALL ABORT1 
      ENDIF


! ----------------------------------------------------------------------

!*    2.2 OUTPUT GRID DEFINITIONS.

      WRITE (IU06,'(/"  OUTPUT GRID"/)')
      IF (IRGG.EQ.1) THEN
        WRITE(IU06,'("   USE A REDUCED GAUSSIAN GRID")')
      ELSE
        WRITE(IU06,'("   USE A REGULAR LON LAT GRID")')
      ENDIF      

      IF(LAQUA) THEN
        LLOBSTRCT=.FALSE.
        WRITE(IU06,*) ' '
        WRITE(IU06,*) '   AN AQUA PLANET WAS REQUESTED '
        WRITE(IU06,*) '   THE BATHYMETRY WILL BE DEEP EVERYWHERE !!!'
        WRITE(IU06,*) ' '
      ENDIF

      IF(LLOBSTRCT) THEN
        WRITE(IU06,*) ' '
        WRITE(IU06,*) '   THE NEW FORM OF BATHYMETRY INPUT IS USED IN'
        WRITE(IU06,*) '   CONJUNCTION WITH THE OBSTRUCTION COEFFICIENT'
        WRITE(IU06,*) ' '
      ENDIF

      CALL ADJUST (AMOWEP, AMOEAP)

!*    SET DIMENSIONS.

      IF(LLGRID) THEN
        XDELLA = (AMONOP-AMOSOP)/(NY-1)
        ALLOCATE(NLONRGG(NY))

        NX = 0
        DO K=1,NY
          KSN=NY-K+1
          READ(IU,*) NLONRGG(KSN)
          NX = MAX(NX,NLONRGG(KSN))
        ENDDO

        IF(IPER.EQ.1) THEN
          XDELLO  = 360._JWRB/REAL(NX)
          AMOEAP = AMOWEP + 360._JWRB - XDELLO
        ELSE
          XDELLO = (AMOEAP-AMOWEP)/(NX-1)
        ENDIF


        NGX = NX
        NGY = NY

        CLOSE(IU)

      ELSE
!       RESET FOR AQUA PLANET IF SELECTED
        IF(LAQUA) THEN
          NY = NINT((AMONOP-AMOSOP)/XDELLA) + 1
          WRITE (IU06,*) ' !! RESETING TO AQUA PLANET CONFIGURATION !!'
          AMONOP=90.0_JWRB
          AMOSOP=-90.0_JWRB
          NY=INT((AMONOP-AMOSOP)/XDELLA)+1
          IRGG=1
        ENDIF

        IPER = 0
        IF (ABS(AMOEAP-AMOWEP+1.5_JWRB*XDELLO) .GE. 360.0_JWRB ) IPER = 1
        NX = NINT((AMOEAP-AMOWEP)/XDELLO) + 1
        NY = NINT((AMONOP-AMOSOP)/XDELLA) + 1

        NGX = NX
        NGY = NY

        ALLOCATE(NLONRGG(NY))

      ENDIF

      WRITE (IU06,'("   RESOLUTION LAT-LON ",2F12.7)') XDELLA, XDELLO
      WRITE (IU06,'("    SOUTHERN LAT  NORTHERN LAT ",                  &
     &  " WESTERN LONG "," EASTERN LONG",                               &
     &  /,2X,4F14.3)') AMOSOP, AMONOP, AMOWEP, AMOEAP
      IF (IPER.EQ.1) WRITE (IU06,*) '   THE GRID IS EAST-WEST PERIODIC'

      IF(CLDOMAIN == '-' ) THEN
        IF (IPER.EQ.1) THEN
          CLDOMAIN = 'g'
        ELSE
          CLDOMAIN = 'm'
        ENDIF 
      ENDIF 

      WRITE (IU06,*) ' '
      WRITE (IU06,*) '   MINIMUM DEPTH IN TABLES IS ', DEPTHA
      WRITE (IU06,*) '   WITH A STEP ', DEPTHD 
      DEPTHMAX=DEPTHA*DEPTHD**(NDEPTH-1)
      WRITE (IU06,*) '   MAXIMUM DEPTH IN TABLES IS ', DEPTHMAX
      WRITE (IU06,*) ' '
 

! ----------------------------------------------------------------------

!*    2.3 OUTPUT GRID CORRECTIONS READ FROM NAMELIST NACORR.
!         --------------------------------------------------

      IF(LLOBSTRCT) THEN
!     when the new bathymetry used, it is assumed that there is not
!     any need for manual correction of the data since it should have
!     been taken care when creating the file. 
        NOUT=0
      ELSE
        REWIND (IU05)
        NOUT=0
        IOS=-9
        IOUTA=500
        ALLOCATE(XOUTW(IOUTA))
        ALLOCATE(XOUTS(IOUTA))
        ALLOCATE(XOUTE(IOUTA))
        ALLOCATE(XOUTN(IOUTA))
        ALLOCATE(NOUTD(IOUTA))
        XOUTS  =-100.0_JWRB
        XOUTN  =-100.0_JWRB
        XOUTW  =   0.0_JWRB
        XOUTE  =   0.0_JWRB
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
      ENDIF
 2300 IF (NOUT.GT.0) THEN
        WRITE (IU06,'(/4X," AREAS TO BE CORRECTED IN OUTPUT GRID",      &
     &   /,4X,"  NO.   SOUTHERN LAT ",                                  &
     &   " NORTHERN LAT  WESTERN LONG ",                                &
     &   " EASTERN LONG  DEPTH")')
        DO I=1,NOUT
          CALL ADJUST (XOUTW(I), XOUTE(I))
          WRITE (IU06,'(4X,I5,1X,4F14.3,I7 )') I, XOUTS(I),             &
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
      OUTLAT =-100.0_JWRB
      OUTLONG=   0.0_JWRB
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
        WRITE (IU06,'(" OUTPUT POINTS FOR SPECTRA AS DEFINED",          &
     &   " BY USER INPUT",/,                                            &
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

      WRITE (IU06,'(" OUTPUT OPTION IS       IFORM =",I3,               &
     &         " (1: UNFORM.  2: FORM.  3: BOTH)")') IFORM
      WRITE (IU06,'(" TEST OUTPUT OPTION IS  ITEST =",I3,               &
     &         " (0: NO  >0: UP TO LEVEL ITEST)")') ITEST
      WRITE (IU06,'(" BLOCK TEST OPTION IS  ITESTB =",I3,               &
     &         " (0: NO  >0: UP TO BLOCK ITESTB)")') ITESTB

! ----------------------------------------------------------------------

!*    2.6 NESTED GRID INFORMATION.
!         ------------------------

      IF ((IBOUNC.EQ.1.OR.IBOUNF.EQ.1) .AND. IRGG.EQ.1) THEN
        WRITE (IU06,*) '**********************************************'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '*       FATAL ERROR IN SUB. UIPREP           *'
        WRITE (IU06,*) '*       ==========================           *'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '* THE NESTING OPTION DOES NOT WORK FOR       *'
        WRITE (IU06,*) '* IRREGULAR ALT-LON GRID                     *' 
        WRITE (IU06,*) '* IRGG = ',IRGG
        WRITE (IU06,*) '* IBOUNC = ',IBOUNC
        WRITE (IU06,*) '* IBOUNF = ',IBOUNF
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '**********************************************'
        CALL ABORT1 
      ENDIF


!*    2.6.1 COARSE GRID OPTION.
!           -------------------

      WRITE (IU06,'(" COARSE GRID OPTION IS IBOUNC = ",I3)') IBOUNC

      IF (IBOUNC .EQ. 1) THEN

        DO I=1,GBOUNC
          CALL ADJUST (AMOWEC(I), AMOEAC(I))
          WRITE (IU06,*) '   THIS IS A COARSE GRID RUN  INFORMATION',   &
     &     ' FOR A FOLLOW UP FINE GRID WILL BE GENERATED'
          WRITE (IU06,'(/4X," NEST AREA IN COARSE GIRD IS",             &
     &     /,4X,"  SOUTHERN LAT  NORTHERN LAT ",                        &
     &     " WESTERN LONG  EASTERN LONG ")')
          WRITE (IU06,'(4X,4F14.3)') AMOSOC(I), AMONOC(I), AMOWEC(I),   &
     &     AMOEAC(I)

!*    2.6.1.1 ARE ALL CORNER POINTS OF THE NEST GRID POINTS?
!             ----------------------------------------------

          WEST = MOD(AMOWEC(I) - AMOWEP + 720.0_JWRB, 360.0_JWRB)
          EAST = MOD(AMOEAC(I) - AMOEAP + 720.0_JWRB, 360.0_JWRB)
          DW=ABS(NINT(WEST/ XDELLO)-(WEST/ XDELLO)) 
          DE=ABS(NINT(EAST/ XDELLO)-(EAST/ XDELLO))
          DS=ABS(NINT((AMOSOC(I) - AMOSOP)/ XDELLA)-((AMOSOC(I)-AMOSOP)/ XDELLA))
          DN=ABS(NINT((AMONOC(I) - AMONOP)/ XDELLA)-((AMONOC(I)-AMONOP)/ XDELLA))

          IF ((DW .GT. 1.E-10_JWRB) .OR.                                &
     &        (DE .GT. 1.E-10_JWRB) .OR.                                &
     &        (DS .GT. 1.E-10_JWRB) .OR.                                &
     &        (DN .GT. 1.E-10_JWRB)) THEN
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
            WRITE (IU06,*) DW, AMOWEC(I), AMOWEP
            WRITE (IU06,*) DE, AMOWEC(I), AMOEAP
            WRITE (IU06,*) DS, AMOSOC(I), AMOSOP
            WRITE (IU06,*) DN, AMONOC(I), AMONOP
            IBOUNC = 0
            CALL ABORT1 
          ENDIF


!*    2.6.1.2 INCLUDES THE COARSE GRID THE NEST GRID?
!             ---------------------------------------

          IF ((IPER.NE.1) .AND.                                          &
     &     (AMOWEP.GT.AMOWEC(I)     .OR. AMOEAC(I).GT.AMOEAP     ) .AND. &
     &     (AMOWEP.GT.AMOWEC(I)+360..OR. AMOEAC(I)+360..GT.AMOEAP) .AND. &
     &     (AMOWEP.GT.AMOWEC(I)-360..OR. AMOEAC(I)-360..GT.AMOEAP)) THEN

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
          IF (AMOSOP .GT. AMOSOC(I) .OR. AMONOC(I) .GT. AMONOP .OR.     &
     &     AMOSOC(I) .GE. AMONOC(I)) THEN
            WRITE (IU06,*) '++++++++++++++++++++++++++++++++++++++++++'
            WRITE (IU06,*) '+                                        +'
            WRITE (IU06,*) '+     WARNING ERROR IN SUB. UIPREP       +'
            WRITE (IU06,*) '+     ============================       +'
            WRITE (IU06,*) '+ ERROR IN NEST SPECIFICATIONS.          +'
            WRITE (IU06,*) '+ OF FINE GRID NUMBER : ', I
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
        ENDDO
        WRITE (IU06,*) '   NUMBERS OF FINE GRIDS ',GBOUNC
      ELSE
        WRITE (IU06,*) '   A NEST IS NOT INCLUDED IN THIS GRID'
      ENDIF

!*    2.6.2 FINE GRID OPTION.
!           -----------------

      WRITE (IU06,'(" FINE GRID OPTION IS   IBOUNF = ",I3)') IBOUNF
      IF (IBOUNF .EQ. 1) THEN
        WRITE (IU06,*) '   THIS IS A FINE GRID RUN, INPUT FROM',        &
     &   ' A COARSE GRID IS EXPECTED'
      ELSE
        WRITE (IU06,*) '   BOUNDARY VALUES FROM A COARSE GRID',         &
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
