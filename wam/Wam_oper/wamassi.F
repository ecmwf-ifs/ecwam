      SUBROUTINE WAMASSI(LDSTOP,LDWRRE)

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
!      S. ABDALLA      ECMWF     OCTOBER 2000 INCLUDE AIR DENSITY & Zi/L
!                                             CALL OUTBETA MOVED TO WAVEMDL

!     PURPOSE.                                                          
!     --------                                                          

!         TO ANALYSE THE WAVE FIELD SUBSTITUTING THE FIRST GUESS        
!         SPECTRA WITH ANALYSED SPECTRA AND FIRST GUESS VALUES          
!         WITH ANALYSED VALUES IN THE GLOBAL GRID                       

!*    INTERFACE.                                                        
!     ----------                                                        

!     *CALL* *WAMASSI(LDSTOP,LDWRRE)*
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
     &            GFLAG20  ,LFDB     ,LOUTINT  ,LOUTSPP  ,LWAMANOUT,
     &            LRSTPARALW
      USE YOWCURR  , ONLY : CDTCUR
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRIBHD, ONLY : PPMISS   ,PPEPS    ,PPREC    ,NTENCODE ,
     &            NGRBRESS ,HOPERS   ,PPRESOL  ,LGRHDIFS ,LNEWLVTP
      USE YOWGRID  , ONLY : NLONRGG  ,IGL      ,IJS      ,IJL2     ,
     &            IJL     , IJSLOC   ,IJLLOC   ,IJGLOBAL_OFFSET
      USE YOWICE   , ONLY : CICOVER  ,CITHICK
      USE YOWINTP  , ONLY : WHGTTG   ,WDIRTG   ,WPKFTG   ,WMNFTG   ,
     &            USTARG   ,UDIRG    ,TAUWG    ,CDG      ,SMEANG   ,
     &            U10G     ,WHGTAG   ,CWHTAG   ,RANGAG   ,MWP1G    ,
     &            MWP2G    ,WSPRDG   ,C4G      ,BFG      ,QPG      ,
     &            DEPTHG   ,HMAXG    ,TMAXG    ,USTOKESG ,VSTOKESG ,
     &            UCURG    ,VCURG    ,PHIEPSG  ,PHIAWG   ,TAUOCG   ,
     &            STRNMSG  ,RHOAG    ,WSTARG   ,CICG     ,CIHG     ,
     &            C3G      ,
     &            WX1G     ,WX2G     ,WX3G     ,WX4G     ,WX5G     ,
     &            NSSTG    ,NICCG    ,NICTG    ,NUCUG    ,NVCUG
      USE YOWINTS  , ONLY : WHGTSG   ,WDIRSG   ,WMNFSG   ,WHGTWG   ,
     &            WDIRWG   ,WMNFWG
      USE YOWINTT  , ONLY : WHGTTRG  ,WDIRTRG  ,MWPMTRG 
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IRGG     ,AMOWEP   ,
     &            AMOSOP   ,AMOEAP   ,AMONOP   ,XDELLA   ,XDELLO
      USE YOWMESPAS, ONLY : LMESSPASS,LFDBIOOUT,LGRIBOUT
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : ZMISS    ,DEG
      USE YOWSPEC, ONLY : NBLKS   ,NBLKE     ,U10NEW   ,U10OLD   ,
     &            THWNEW   ,THWOLD   ,USNEW    ,USOLD    ,Z0NEW    ,
     &            Z0OLD    ,TAUW     ,
     &            ROAIRN   ,ROAIRO   ,ZIDLNEW  ,ZIDLOLD  ,
     &            FL1
      USE YOWSTAT  , ONLY : CDATEE   ,CDATEF   ,CDTPRO   ,CDTRES   ,
     &            CDATER   ,CDATES   ,CDTINTT  ,CDTINTS  ,CDTSPT   ,
     &            CDTSPS   ,CFDBSF   ,CFDB2DSP ,NWFDBREF ,IDELWI   ,
     &            IREST    ,IREFRA   ,IASSI    ,NTASKS   ,NSIZE    ,
     &            NENSFNB  ,NTOTENS  ,NTOTENS  ,NSYSNB   ,CDATEA   ,
     &            MARSTYPE ,YCLASS   ,YEXPVER  ,LALTAS   ,LSARAS   ,
     &            LSARINV  ,IREFDATE ,LFDBOPEN
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH    ,CWI
      USE YOWUNIT  , ONLY : IU11     ,IU12     ,IU13     ,IU14     ,
     &            IU22     ,IU23     ,IU25     ,IU26     ,IU27     ,
     &            IU28     ,IU04     ,IU32
      USE YOWWAMI  , ONLY : CBPLTDT  ,CEPLTDT  ,IANALPD  ,IFOREPD  ,
     &            IDELWIN  ,NFCST    ,ISTAT
      USE MPL_MODULE
      USE FDBSUBS_MOD
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: IFIL, IJ, K, M
      INTEGER :: IG
      INTEGER :: I_GET_UNIT
      INTEGER :: JKGLO,KIJS,KIJL,NPROMA
      INTEGER :: IERR
      INTEGER, DIMENSION(IJSLOC:IJLLOC) :: MIJ

      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION(IJSLOC:IJLLOC,NANG,NFRE) :: XLLWS

      LOGICAL :: LDSTOP, LDWRRE, LLFLUSH
      LOGICAL :: LSV, LOUT
      LOGICAL :: LLOUTBS

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WAMASSI',0,ZHOOK_HANDLE)
#endif

      WRITE (IU06,*) '  '
      WRITE (IU06,*) '  '
      WRITE (IU06,*)
     &   ' START OF DATA ASSIMILATION: DATE IS CDTPRO: ', CDTPRO
      WRITE (IU06,*) '  '
      LLFLUSH = .TRUE.
      LLOUTBS=.FALSE.

!     1. CALLING THE MAIN DATA ASSIMILATION ROUTINE
!        ------------------------------------------

      IF(LSARINV) THEN
        CALL SARINVERT
        IF (ITEST.GE.1) WRITE(IU06,*) '   SUB. WAMASSI: SARINVERT DONE'
      ENDIF

      IF(LSARAS) THEN
        CALL SARAS
        IF (ITEST.GE.1) WRITE(IU06,*) '   SUB. WAMASSI: SARAS DONE'
      ENDIF

      IF(LALTAS) THEN
        CALL ALTAS(IJSLOC, IJLLOC, NANG, NFRE, MIJ, XLLWS)
        IF (ITEST.GE.1) WRITE(IU06,*) '   SUB. WAMASSI: ALTAS DONE'
      ENDIF


!*    2. PREPARING OUTPUT AND UPDATE TO CHARNOCK PARAMETER
!        -------------------------------------------------


!       UPDATE U10NEW,THWNEW,USNEW
!       --------------------------

        IG=1

        DO IJ=IJSLOC, IJLLOC
          U10NEW(IJ) = U10OLD(IJ,IG)
          THWNEW(IJ) = THWOLD(IJ,IG)
          USNEW(IJ) = USOLD(IJ,IG)
          ROAIRN(IJ) = ROAIRO(IJ,IG)
          ZIDLNEW(IJ) = ZIDLOLD(IJ,IG)
        ENDDO

!       MODEL OUTPUT INTEGRATED DATA ARE SAVED IN MODULES.
!       --------------------------------------------------
        IF (CDTINTT.EQ.CDTPRO .OR. CDTINTS.EQ.CDTPRO .OR.
     &       CDTSPT .EQ.CDTPRO .OR. CDTSPS .EQ.CDTPRO    ) THEN

          CALL OUTBS (IJSLOC, IJLLOC, IJGLOBAL_OFFSET,
     &                MIJ(IJSLOC),
     &                FL1(IJSLOC:IJLLOC,:,:), XLLWS(IJSLOC:IJLLOC,:,:),
     &                IG, IU25, IU26, LLOUTBS)


          IF (ITEST.GE.2) THEN
            WRITE(IU06,*) '   SUB. WAMASSI: OUTPUT DONE IN SUB OUTBS'
          ENDIF
        ENDIF


!*    3.1 SAVE FILES.
!         -----------

!     3.1.1 WRITE INTEGRATED DATA TO FILE AND/OR PRINTER
!           DATA WERE COLLECTED INSIDE THE BLOCK LOOP.
!           --------------------------------------------

      IF (LWAMANOUT .AND. 
     &    (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO) ) THEN

!       PRINTS OUT NORMS
        CALL OUTWNORM

        CALL OUTINT (IU06, IU22, IU23, IU32, ITEST, NWFDBREF, LFDBOPEN)

        IF(LMESSPASS) THEN
          CALL GSTATS(752,0)
          CALL MPL_BARRIER(CDSTRING='WAMASSI 1:')
          CALL GSTATS(752,1)
        ENDIF
        IF (ITEST.GE.2) WRITE(IU06,*) '   SUB. WAMASSI: OUTINT DONE'
      ENDIF

      IF (CDATEE.EQ.CDTPRO.OR.LDWRRE) THEN

!*    3.1.2 STORE 2-D SPECTRA AT SELECTED GRID POINTS,
!           INTEGRATED PARAMETERS OF ENTIRE GRID.
!           -------------------------------------

!           SAVE INTEGRATED DATA FILE.

        IF (GFLAG20 .AND. .NOT. LFDB .AND. LOUTINT) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            CALL GSFILE (IU06, IU32, 0, CDTPRO, CDATEF, 'APP', 'S')
          ENDIF
        ENDIF
        IF (FFLAG20 .AND. LOUTINT) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.NPROC)).OR.(.NOT.LMESSPASS)) THEN
            CALL GSFILE (IU06, IU22, 0, CDTPRO, CDATEF, 'AMP', 'S')
          ENDIF
        ENDIF

!           SAVE INTEGRATED SWELL AND SEA DATA FILE.
!           ----------------------------------------

        IF (FFLAG21 .AND. LOUTINT) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.MIN(2,NPROC)))
     &     .OR.(.NOT.LMESSPASS)) THEN
            CALL GSFILE (IU06, IU23, 0, CDTPRO, CDATEF, 'ASS', 'S')
          ENDIF
        ENDIF

!           SAVE SPECTRA FILE.

        IF (FFLAG25 .AND. LOUTSPP) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            CALL GSFILE (IU06, IU27, 0, CDTPRO, CDATEF, 'AUT', 'S')
          ENDIF
        ENDIF

!           SAVE SWELL SPECTRA FILE. 

        IF (FFLAG26 .AND. LOUTSPP) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            CALL GSFILE (IU06, IU28, 0, CDTPRO, CDATEF, 'ASW', 'S')
          ENDIF
        ENDIF
      ENDIF

      LSV=(CDTRES.EQ.CDTPRO .OR. CDATEE.EQ.CDTPRO 
     &                      .OR. CDTPRO.EQ.CDATER)

      IF (LSV.OR.LDWRRE) THEN

!*    3.1.4 SAVE RESTART FIELDS.
!           --------------------
        LOUT=(IREST.EQ.1 
     &        .AND. (CDTPRO.EQ.CDATER .OR. CDTPRO.LE.CDATES)
     &        .AND. LWAMANOUT )
        IF ( LOUT .OR. LDWRRE ) THEN 

!           SAVE SPECTRUM
!           -------------

          IF ( LOUT.AND.LGRIBOUT ) THEN
!           SAVE SPECTRUM IN GRIB
            CALL OUTWSPEC(FL1)
            CALL GSTATS(1976,0)
            IF (NWFDBREF.NE.-5 ) IERR = IFLUSHFDBSUBS (NWFDBREF)
            CALL GSTATS(1976,1)
            CALL GSTATS(752,0)
            CALL MPL_BARRIER(CDSTRING='WAMASSI 2:')
            CALL GSTATS(752,1)
            LLFLUSH = .FALSE.
            WRITE(IU06,*) ' '
            WRITE(IU06,*) '  GRIB WAVE SPECTRA DISPOSED AT........',
     &       ' CDTPRO  = ', CDTPRO
            WRITE(IU06,*) ' '
            CALL FLUSH(IU06)
          ENDIF

          IF ( .NOT.LGRIBOUT .OR. LDWRRE ) THEN
            CALL SAVSTRESS(U10OLD, THWOLD, USOLD, TAUW, Z0OLD,
     &                     ROAIRO, ZIDLOLD, CICOVER, CITHICK,
     &                     NBLKS,NBLKE,CDTPRO,CDATEF)
            WRITE(IU06,*) '  BINARY STRESS FILE DISPOSED AT........',
     &       ' CDTPRO  = ', CDTPRO
            WRITE(IU06,*) ' '

            CALL SAVSPEC(FL1,NBLKS,NBLKE,CDTPRO,CDATEF,CDATER)
            WRITE(IU06,*) '  BINARY WAVE SPECTRA DISPOSED AT........',
     &       ' CDTPRO  = ', CDTPRO
            WRITE(IU06,*) ' '
            CALL FLUSH(IU06)
          ENDIF


!*    3.1.5 UPDATE, WRITE AND SAVE WAMINFO FILE.
!           ------------------------------------

          IF ((LDSTOP.OR.LDWRRE) .AND. 
     &        (LMESSPASS .AND. IRANK .EQ. 1 .OR. .NOT. LMESSPASS )) THEN
            CALL DIFDATE (CDATEF,CDATEE,IFOREPD)
            CALL DIFDATE (CDTPRO,CDATEF,IANALPD)
            IF (CDTPRO.LE.CDATEF) THEN
              NFCST = 1
            ELSE
              NFCST = 0
              CALL DIFDATE (CDTPRO,CDATEE,IFOREPD)
            ENDIF
            ISTAT(:) = 0
            IF (CDTPRO.EQ.CDATEE) ISTAT(1) = 1
            IDELWIN = IDELWI
            CBPLTDT = CDATEF
            CEPLTDT = CDATEE

            IU04 =  I_GET_UNIT (IU06,CWI(1:ICPLEN+8) , 'w', 'f', 0)

            CALL WRITSTA (IU04, CDTPRO, CDATEE, IANALPD, IFOREPD,
     &                    IDELWIN, CDATER, CDATES, CBPLTDT, CEPLTDT,
     &                    IASSI, NFCST, ISTAT, CDTCUR,
     &                    LRSTPARALW, NPROC)

            CLOSE (IU04)
            WRITE(IU06,*) ' WAMINFO FILE WRITTEN FOR RESTART...',
     &                      ' CDTPRO  = ', CDTPRO
            WRITE(IU06,*) '                                    ',
     &                      ' CDATEF  = ', CDATEF
            WRITE(IU06,*) ' TO ', CWI(1:ICPLEN+8)
          ENDIF
        ENDIF

        IF(LMESSPASS) THEN
          CALL GSTATS(752,0)
          CALL MPL_BARRIER(CDSTRING='WAMASSI 3:')
          CALL GSTATS(752,1)
        ENDIF
      ENDIF

      IF (LFDB .AND. LLFLUSH .AND. NWFDBREF.NE.-5 .AND. (
     &  (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO))) THEN
        CALL GSTATS(1976,0)
        IERR = IFLUSHFDBSUBS (NWFDBREF)
        CALL GSTATS(1976,1)
        WRITE(IU06,*) '   DB ', NWFDBREF , ' FLUSHED AT ',
     &    CDTPRO, ' FROM WAMASSI. '
        CALL FLUSH (IU06)
      ENDIF
      CALL GSTATS(752,0)
      CALL MPL_BARRIER(CDSTRING='WAMASSI 4:')
      CALL GSTATS(752,1)

!     DEALLOCATE OUTPUT ARRAYS
      CALL CLEAN_OUTBS(LLOUTBS)

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WAMASSI',1,ZHOOK_HANDLE)
#endif

      RETURN
      END SUBROUTINE WAMASSI
