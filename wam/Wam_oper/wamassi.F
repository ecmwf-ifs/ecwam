      SUBROUTINE WAMASSI(BETAG,NLONW,NLATW,LDSTOP,LDWRRE)

! ----------------------------------------------------------------------

!****  *WAMASSI* - SUPERVISES EXECUTION OF MAIN MODULES                 
!****              OF THE WAVE DATA ASSIMILATION.                       

!      P. LIONELLO     ECMWF     APRIL   1990
!      J. BIDLOT       ECMWF     FEBRARY 1996 MESSAGE PASSING.
!      J. BIDLOT       ECMWF     FEBRARY 1997 MODULE and ATMOSPHERIC
!                                COUPLING.
!      J. BIDLOT       ECMWF     MARCH   1997  ADD SAVSTRESS AND SAVSPEC
!      B. HANSEN       ECMWF     APRIL   1997  SIGNAL HANDLING.
!      B. HANSEN       ECMWF     JANUARY 1998  1.1 OUTSOURCED.

!     PURPOSE.                                                          
!     --------                                                          

!         TO ANALYSE THE WAVE FIELD SUBSTITUTING THE FIRST GUESS        
!         SPECTRA WITH ANALYSED SPECTRA AND FIRST GUESS VALUES          
!         WITH ANALYSED VALUES IN THE GLOBAL GRID                       

!*    INTERFACE.                                                        
!     ----------                                                        

!     *CALL* *WAMASSI(BETAG,NLONW,NLATW,LDSTOP,LDWRRE)
!        *BETAG*     REAL      CHARNOCK PARAMETER TO BE PASSED BACK TO
!                              THE ATMOSPHERIC MODEL
!        *NLONW*     INTEGER   FIRST DIMENSION OF BETAG
!        *NLATW*     INTEGER   SECOND DIMENSION OF BETAG
!        *LDSTOP*    LOGICAL   SET .TRUE. IF STOP SIGNAL RECEIVED.
!        *LDWRRE*    LOGICAL   SET .TRUE. IF RESTART SIGNAL RECEIVED.

!     METHOD.                                                           
!     -------                                                           

!     EXTERNALS.                                                        
!     ----------                                                        

!        *ANALYSE*  - UTILITY OF OIFIELD . IT PRODUCES THE OPTIMALLY    
!                     INTERPOLATED VALUE AT A SINGLE GRIDPOINT.         
!        *FDUR*     - COMPUTES WINDSEA DURATION.                        
!        *FEMEAN*   - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT. 
!        *FUSTAR*   - COMPUTES ANALYSED USTAR ACCORDING TO ANALYSED     
!                     WINDSEA ENERGY.                                   
!        *FWSEA*    - DETERMINES IF WINDSEA IS PRESENT IN THE SPECTRUM,,
!                     COMPUTING EVENTUALLY WINDSEA ENERGY AND MEAN FREQ.
!        *GRFIELD*  - READS GRIDDED MEASUREMENTS AND DISTRIBUTES THEM.
!        *GRDATA*   - TRANSFERS MEASUREMENTS TO THE MODEL GRID.         
!        *GSFILE*   - ROUTINE TO DYNAMICALLY FETCH OR DISPOSE FILES.    
!        *INCDATE*  - UPDATE DATE TIME GROUP.                           
!        *MAKEBLO*  - MAKES BLOCKS FROM GRIDDED FIEDS.                  
!        *MAKEGRID* - MAKES GRIDDED FIELDS.                             
!        *OIFIELD*  -  MERGES MODEL AND MEASUREMENTS BY O.I.            
!        *OUTAN*    - TRANSFERS INTEGRATED ANALYSED QUANTITIES TO FILE. 
!        *PEAKFR*   - COMUTE PEAK FREQUENCY.                            
!        *SAVSTRESS*- DISPOSE STRESS/WIND RESTART FILES.
!        *SAVSPEC   - DISPOSE SPECTRUM RESTART FILES.
!        *SEMEAN*   - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.   
!        *SERMAT*   - UTILITY OF OIFIELD TO SUM TWO MATRICES.           
!        *SPLITBL*  - STORING OVERLAPPING LATITUDES OF BLOCKS.          
!        *STHQ*     - COMPUTATION OF MEAN WAVE DIRECTION AT EACH        
!                     GRID POINT.                                       
!        *UPDATE*   - ANALYSES THE WAM MODEL SPECTRA.                   
!        *UPWSPEC*  - PRODUCES THE UPDATED SPECTRUM.                    
!        *WSMFEN*   - UTILITY OF FWSEA , COMPUTING  WINDSEA ENERGY AND  
!                     DURATION IN EACH GRID POINT.                      
!        *OUTWNORM* - COMPUTES A FEW NORMS OF GRIDDED FIELDS

!     REFERENCES.                                                       
!     -----------                                                       

!          NONE                                                         

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWCOUT  , ONLY : FFLAG20  ,FFLAG21  ,FFLAG25  ,FFLAG26  ,
     &            GFLAG20  ,LFDB     ,IPFGTBL
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3   ,
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,KSEK0    ,KSEK1    ,
     &            KSEK2    ,KSEK3    ,KSEK4    ,PSEK2    ,PSEK3    ,
     &            PPMISS   ,PPEPS    ,PPREC    ,NTENCODE
      USE YOWGRID  , ONLY : NLONRGG  ,IGL      ,IJS      ,IJL2     ,IJL
      USE YOWINTP  , ONLY : WHGTTG   ,WDIRTG   ,WPKFTG   ,WMNFTG   ,
     &            USTARG   ,UDIRG    ,TAUWG    ,CDG      ,SMEANG   ,
     &            U10G     ,WHGTAG   ,CWHTAG   ,RANGAG   ,
     &            MWP1G    ,MWP2G    ,WSPRDG
      USE YOWINTS  , ONLY : WHGTSG   ,WDIRSG   ,WMNFSG   ,WHGTWG   ,
     &            WDIRWG   ,WMNFWG
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IRGG     ,AMOWEP   ,
     &            AMOSOP   ,AMOEAP   ,AMONOP   ,XDELLA   ,XDELLO
      USE YOWMESPAS, ONLY : LMESSPASS,LFDBIOOUT,LGRIBOUT
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      ,
     &            NBLO     ,NIBLO    ,CLDOMAIN
      USE YOWPCONS , ONLY : ZMISS    ,DEG
      USE YOWSPEC, ONLY : NSTART   ,NEND     ,U10NEW   ,U10OLD   ,
     &            THWNEW   ,THWOLD   ,USNEW    ,USOLD    ,Z0NEW    ,
     &            Z0OLD    ,TAUW     ,FL1      ,FL2      ,FL3
      USE YOWSTAT  , ONLY : CDATEE   ,CDATEF   ,CDTPRO   ,CDTRES   ,
     &            CDATER   ,CDATES   ,CDTINTT  ,CDTINTS  ,CDTSPT   ,
     &            CDTSPS   ,CFDBSF   ,CFDB2DSP ,MFDBSF   ,IDELWI   ,
     &            IREST    ,IASSI    ,NTASKS   ,NSIZE    ,NENSFNB  ,
     &            NTOTENS  ,NTOTENS  ,NSYSNB   ,CDATEA   ,MARSTYPE ,
     &            YCLASS   ,YEXPVER  ,LALTAS   ,LSARAS
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH    ,CWI
      USE YOWUNIT  , ONLY : IU11     ,IU12     ,IU13     ,IU14     ,
     &            IU22     ,IU23     ,IU27     ,IU28     ,IU04     ,IU32
      USE YOWWAMI  , ONLY : CBPLTDT  ,CEPLTDT  ,IANALPD  ,IFOREPD  ,
     &            IDELWIN  ,NFCST    ,ISTAT

! ----------------------------------------------------------------------
      REAL BETAG(NLONW,NLATW)
! ----------------------------------------------------------------------

      LOGICAL LDSTOP, LDWRRE, LLFLUSH

      WRITE (IU06,*) '  '
      WRITE (IU06,*)
     &   ' START OF DATA ASSIMILATION: DATE IS CDTPRO: ', CDTPRO
      LLFLUSH = .TRUE.

      IF(.NOT.ALLOCATED(WHGTAG)) ALLOCATE(WHGTAG(NGX,NGY))
      IF(.NOT.ALLOCATED(CWHTAG)) ALLOCATE(CWHTAG(NGX,NGY))
      IF(.NOT.ALLOCATED(RANGAG)) ALLOCATE(RANGAG(NGX,NGY))

!     1. CALLING THE MAIN DATA ASSIMILATION ROUTINE
!        ------------------------------------------

      IF(LALTAS) THEN
        CALL ALTAS
        IF (ITEST.GE.1) WRITE(IU06,*) '   SUB. WAMASSI: ALTAS DONE'
      ENDIF

      IF(LSARAS) THEN
        CALL SARAS
        IF (ITEST.GE.1) WRITE(IU06,*) '   SUB. WAMASSI: SARAS DONE'
      ENDIF


!*    2. PREPARING OUTPUT AND UPDATE TO CHARNOCK PARAMETER
!        -------------------------------------------------

      IF (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO .OR.
     &    CDTSPT.EQ.CDTPRO .OR. CDTSPS.EQ.CDTPRO ) THEN
        IF(.NOT.LMESSPASS) THEN
          REWIND IU12
          IF (IGL.NE.1) THEN
            REWIND IU11
            REWIND IU12
            REWIND IU13
            REWIND IU14
            IFIL=IU11
            IU11=IU12
            IU12=IFIL
            IFIL=IU13
            IU13=IU14
            IU14=IFIL
            READ(IU11) FL1
            IF (ITEST.GE.2)
     &       WRITE(IU06,*) '   SUB. WAMASSI: START OF READING',
     &       ' FIRST BLOCK OF SPECTRA'
          ENDIF
        ENDIF
      ENDIF

!*    2.1 LOOP OVER THE BLOCKS.
!         ---------------------
      BLOCK : DO IG=1,IGL

!       UPDATE U10NEW,THWNEW,USNEW
!       --------------------------
        DO IJ=IJS(IG),IJL(IG)
          U10NEW(IJ) = U10OLD(IJ,IG)
          THWNEW(IJ) = THWOLD(IJ,IG)
          USNEW(IJ) = USOLD(IJ,IG)
        ENDDO


!       MODEL OUTPUT INTEGRATED DATA ARE SAVED IN MODULES.
!       --------------------------------------------------
        IF (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO .OR.
     &      CDTSPT.EQ.CDTPRO .OR. CDTSPS.EQ.CDTPRO ) THEN

          IF (MOD(IG,2).EQ.1) THEN
            IF (.NOT.LMESSPASS.AND.IG.NE.IGL) READ(IU11) FL2
            CALL OUTBS (FL1, FL3, IJS(IG), IJL(IG), IG, IGL,
     &        IU27, IU28,THWOLD,USOLD,U10NEW,THWNEW,USNEW,TAUW,
     &        NSTART,NEND)
            IF (ITEST.GE.2.AND.ITESTB.GE.IG) THEN
              WRITE(IU06,*) '   SUB. WAMASSI: OUTPUT DONE IN SUB OUTBS'
            ENDIF
          ELSE
            IF (IG.NE.IGL) READ(IU11) FL1
            CALL OUTBS (FL2, FL3, IJS(IG), IJL(IG), IG, IGL,
     &        IU27, IU28,THWOLD,USOLD,U10NEW,THWNEW,USNEW,TAUW,
     &        NSTART,NEND)
            IF (ITEST.GE.2.AND.ITESTB.GE.IG) THEN
              WRITE(IU06,*) '   SUB. WAMASSI: OUTPUT DONE IN SUB OUTBS'
            ENDIF
          ENDIF
        ENDIF

!       UPDATING THE CHARNOCK PARAMETER FIELD.
!       -------------------------------------
        IF(LWCOU) THEN
          CALL OUTBETA ( IJS(IG), IJL(IG), IG, IGL,USNEW,Z0NEW,
     &     NSTART, NEND, BETAG, NLONW, NLATW)
        ENDIF

!     BRANCHING BACK TO 2.1 FOR NEXT BLOCK OF LATITUDES
      ENDDO BLOCK

!*    3.1 SAVE FILES.
!         -----------

!     3.1.1 WRITE INTEGRATED DATA TO FILE AND/OR PRINTER
!           DATA WERE COLLECTED INSIDE THE BLOCK LOOP.
!           --------------------------------------------

      IF (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO) THEN

!       PRINTS OUT NORMS
        CALL OUTWNORM

        CALL OUTINT (IU06, IU22, IU23, IU32, ITEST, MFDBSF)

        IF(LMESSPASS) THEN
          CALL MPE_BARRIER(IERR)
          IF(IERR.LT.0) CALL MPEI_ABORT('MPE_BARRIER ERROR in WAMASSI')
        ENDIF
        IF (ITEST.GE.2) WRITE(IU06,*) '   SUB. WAMASSI: OUTINT DONE'
      ENDIF

      IF (CDATEE.EQ.CDTPRO.OR.LDWRRE) THEN

!*    3.1.2 STORE 2-D SPECTRA AT SELECTED GRID POINTS,
!           INTEGRATED PARAMETERS OF ENTIRE GRID.
!           -------------------------------------

!           SAVE INTEGRATED DATA FILE.

        IF (GFLAG20 .AND. .NOT. LFDB) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            CALL GSFILE (IU06, IU32, 0, CDTPRO, CDATEF, 'APP', 'S')
          ENDIF
        ENDIF
        IF (FFLAG20) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.NPROC)).OR.(.NOT.LMESSPASS)) THEN
            CALL GSFILE (IU06, IU22, 0, CDTPRO, CDATEF, 'AMP', 'S')
          ENDIF
        ENDIF

!           SAVE INTEGRATED SWELL AND SEA DATA FILE.
!           ----------------------------------------

        IF (FFLAG21) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.MIN(2,NPROC)))
     &     .OR.(.NOT.LMESSPASS)) THEN
            CALL GSFILE (IU06, IU23, 0, CDTPRO, CDATEF, 'ASS', 'S')
          ENDIF
        ENDIF

!           SAVE SPECTRA FILE.

        IF (FFLAG25) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            CALL GSFILE (IU06, IU27, 0, CDTPRO, CDATEF, 'AUT', 'S')
          ENDIF
        ENDIF

!           SAVE SWELL SPECTRA FILE. 

        IF (FFLAG26) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            CALL GSFILE (IU06, IU28, 0, CDTPRO, CDATEF, 'ASW', 'S')
          ENDIF
        ENDIF
      ENDIF

      IF (CDTRES.EQ.CDTPRO.OR.CDATEE.EQ.CDTPRO.OR.CDTPRO.EQ.CDATER.OR.
     &    LDWRRE) THEN


!*    3.1.4 SAVE RESTART FIELDS.
!           --------------------

        IF ( (IREST.EQ.1.AND.(CDTPRO.EQ.CDATER .OR. CDTPRO.LE.CDATES))
     &       .OR.
     &        LDWRRE ) THEN 

!*          STRESS FIELD IS DUMPED FOR RESTART.
!           -----------------------------------
!           NOTE :
!           STRESS FIELD IS DUMPED FOR RESTART
!           ONLY IF CURRENT TIME CORRESPONDS TO CDATER

          IF((CDTPRO .EQ. CDATER  ).OR.
     &       (CDATER .EQ. '000000000000' ).OR.
     &        LDWRRE) THEN

            IF ( .NOT.LGRIBOUT ) THEN
              CALL SAVSTRESS(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,NSTART,NEND,
     &                       CDTPRO,CDATEF)
            ENDIF

          ENDIF

!           SAVE SPECTRUM
!           -------------

          IF ( LGRIBOUT ) THEN
            CALL SPEC2FDB (IU06, NPROC, IRANK, LWCOU, TH, FR,
     &                     FL1, NINF, NSUP, NANG, NFRE, NGX, NGY,
     &                     NLONRGG, KXLT, IXLG, NIBLO, NBLO,
     &                     AMONOP, XDELLA, NSTART, NEND,
     &                     NKSEK0, NKSEK1, NKSEK2, NKSEK3, NKSEK4,
     &                     NPSEK2, NPSEK3, KSEK0, KSEK1, KSEK2, KSEK3,
     &                     KSEK4, PSEK2, PSEK3, NTENCODE,
     &                     DEG, PPMISS, PPEPS, PPREC,
     &                     CDTPRO, CDATEF, CDATEA,
     &                     CLDOMAIN, MARSTYPE, YCLASS, YEXPVER,
     &                     CFDB2DSP, NENSFNB, NTOTENS, NSYSNB, NMETNB,
     &                     LFDBIOOUT,ICPLEN, CPATH, ITEST, MFDBSF)
            ISTAT = IFLUSHFDB (MFDBSF)
            LLFLUSH = .FALSE.
          ELSE
            CALL SAVSPEC(FL1,NSTART,NEND,CDTPRO,CDATEF,CDATER,
     &                   MARSTYPE,YCLASS, YEXPVER, CFDB2DSP)
          ENDIF
          WRITE(IU06,*) ' '
          WRITE(IU06,*) '  WAVE SPECTRA FIELDS DISPOSED AT........',
     &                  ' CDTPRO  = ', CDTPRO
          WRITE(IU06,*) ' '


!*    3.1.5 UPDATE, WRITE AND SAVE WAMINFO FILE.
!           ------------------------------------

          IF (LDSTOP .AND. 
     &        (LMESSPASS .AND. IRANK .EQ. 1 .OR. .NOT. LMESSPASS )) THEN
            ICH = 7
            CALL DIFDATE (CDATEF,CDATEE,IFOREPD)
            CALL DIFDATE (CDTPRO,CDATEF,IANALPD)
            IF (CDTPRO.LE.CDATEF) THEN
              NFCST = 1
            ELSE
              NFCST = 0
              CALL DIFDATE (CDTPRO,CDATEE,IFOREPD)
            ENDIF
            ISTAT = 0
            IF (CDTPRO.EQ.CDATEE) ISTAT(1) = 1
            IDELWIN = IDELWI
            CBPLTDT = CDATEF
            CEPLTDT = CDATEE

            OPEN (IU04,FILE=CWI(1:ICPLEN+8),
     &            STATUS='UNKNOWN',FORM='FORMATTED')
            CALL WRITSTA (IU04, CDTPRO, CDATEE, IANALPD, IFOREPD,
     &                    IDELWIN, CDATER, CDATES, CBPLTDT, CEPLTDT,
     &                    IASSI, NFCST, ISTAT)

            CLOSE (IU04)
            WRITE(IU06,*) ' WAMINFO FILE WRITTEN FOR RESTART...',
     &                      ' CDTPRO  = ', CDTPRO
            WRITE(IU06,*) '                                    ',
     &                      ' CDATEF  = ', CDATEF
            WRITE(IU06,*) ' TO ', CWI(1:ICPLEN+8)
          ENDIF
        ENDIF

        IF(LMESSPASS) THEN
          CALL MPE_BARRIER(IERR)
          IF(IERR.LT.0) CALL MPEI_ABORT('MPE_BARRIER ERROR IN  WAMASSI')
        ENDIF
      ENDIF
      IF (LFDB .AND. LLFLUSH .AND. MFDBSF.NE.-5 .AND. (
     &  (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO))) THEN
        ISTAT = IFLUSHFDB (MFDBSF)
        WRITE(IU06,*) ' DB ', MFDBSF , ' FLUSHED AT ',
     &    CDTPRO, ' FROM WAMASSI. '
        CALL FLUSH (IU06)
      ENDIF

!     DEALLOCATE OUTPUT ARRAYS

      IF(ALLOCATED(WHGTTG)) DEALLOCATE(WHGTTG)
      IF(ALLOCATED(WDIRTG)) DEALLOCATE(WDIRTG)
      IF(ALLOCATED(WPKFTG)) DEALLOCATE(WPKFTG)
      IF(ALLOCATED(WMNFTG)) DEALLOCATE(WMNFTG)
      IF(ALLOCATED(USTARG)) DEALLOCATE(USTARG)
      IF(ALLOCATED(UDIRG)) DEALLOCATE(UDIRG)
      IF(ALLOCATED(TAUWG)) DEALLOCATE(TAUWG)
      IF(ALLOCATED(CDG)) DEALLOCATE(CDG)
      IF(ALLOCATED(SMEANG)) DEALLOCATE(SMEANG)
      IF(ALLOCATED(U10G)) DEALLOCATE(U10G)
      IF(ALLOCATED(WHGTAG)) DEALLOCATE(WHGTAG)
      IF(ALLOCATED(CWHTAG)) DEALLOCATE(CWHTAG)
      IF(ALLOCATED(RANGAG)) DEALLOCATE(RANGAG)
      IF(ALLOCATED(MWP1G)) DEALLOCATE(MWP1G)
      IF(ALLOCATED(MWP2G)) DEALLOCATE(MWP2G)
      IF(ALLOCATED(WSPRDG)) DEALLOCATE(WSPRDG)

      IF(ALLOCATED(WHGTSG)) DEALLOCATE(WHGTSG)
      IF(ALLOCATED(WDIRSG)) DEALLOCATE(WDIRSG)
      IF(ALLOCATED(WMNFSG)) DEALLOCATE(WMNFSG)
      IF(ALLOCATED(WHGTWG)) DEALLOCATE(WHGTWG)
      IF(ALLOCATED(WDIRWG)) DEALLOCATE(WDIRWG)
      IF(ALLOCATED(WMNFWG)) DEALLOCATE(WMNFWG)

      RETURN
      END SUBROUTINE WAMASSI
