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
     &            GFLAG20  ,LFDB     ,IPFGTBL  ,LOUTINT  ,LOUTSPP
      USE YOWCURR  , ONLY : CDTCUR
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRIBHD, ONLY : PPMISS   ,PPEPS    ,PPREC    ,NTENCODE ,
     &            NGRBRESS ,HOPERS   ,PPRESOL  ,LGRHDIFS ,LNEWLVTP
      USE YOWGRID  , ONLY : NLONRGG  ,IGL      ,IJS      ,IJL2     ,IJL
      USE YOWICE   , ONLY : CICOVER  ,CITHICK
      USE YOWINTP  , ONLY : WHGTTG   ,WDIRTG   ,WPKFTG   ,WMNFTG   ,
     &            USTARG   ,UDIRG    ,TAUWG    ,CDG      ,SMEANG   ,
     &            U10G     ,WHGTAG   ,CWHTAG   ,RANGAG   ,MWP1G    ,
     &            MWP2G    ,WSPRDG   ,C4G      ,BFG      ,QPG      ,
     &            DEPTHG   ,HMAXG    ,TMAXG    ,USTOKESG ,VSTOKESG ,
     &            UCURG    ,VCURG    ,PHIEPSG  ,PHIAWG   ,TAUOCG   ,
     &            STRNMSG  ,RHOAG    ,WSTARG   ,CICG     ,CIHG     ,
     &            WX1G     ,WX2G     ,WX3G     ,WX4G     ,WX5G
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
     &            Z0OLD    ,TAUW     ,
     &            ROAIRN   ,ROAIRO   ,ZIDLNEW  ,ZIDLOLD  ,
     &            FL1      ,FL2      ,FL3      ,SL
      USE YOWSTAT  , ONLY : CDATEE   ,CDATEF   ,CDTPRO   ,CDTRES   ,
     &            CDATER   ,CDATES   ,CDTINTT  ,CDTINTS  ,CDTSPT   ,
     &            CDTSPS   ,CFDBSF   ,CFDB2DSP ,NWFDBREF ,IDELWI   ,
     &            IREST    ,IREFRA   ,IASSI    ,NTASKS   ,NSIZE    ,
     &            NENSFNB  ,NTOTENS  ,NTOTENS  ,NSYSNB   ,CDATEA   ,
     &            MARSTYPE ,YCLASS   ,YEXPVER  ,LALTAS   ,LSARAS   ,
     &            LSARINV  ,IREFDATE ,LFDBOPEN ,NWAM_BLKS
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH    ,CWI
      USE YOWUNIT  , ONLY : IU11     ,IU12     ,IU13     ,IU14     ,
     &            IU22     ,IU23     ,IU27     ,IU28     ,IU04     ,IU32
      USE YOWWAMI  , ONLY : CBPLTDT  ,CEPLTDT  ,IANALPD  ,IFOREPD  ,
     &            IDELWIN  ,NFCST    ,ISTAT
      USE MPL_MODULE
      USE FDBSUBS_MOD
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      LOGICAL LDSTOP, LDWRRE, LLFLUSH
      LOGICAL :: LSV, LOUT, LLINTPOL
      REAL ZHOOK_HANDLE
! Mod for OPENMP
      INTEGER JKGLO,KIJS,KIJL,NPROMA
! End Mod for OPENMP


      IF (LHOOK) CALL DR_HOOK('WAMASSI',0,ZHOOK_HANDLE)

      WRITE (IU06,*) '  '
      WRITE (IU06,*) '  '
      WRITE (IU06,*)
     &   ' START OF DATA ASSIMILATION: DATE IS CDTPRO: ', CDTPRO
      WRITE (IU06,*) '  '
      LLFLUSH = .TRUE.

      IF(.NOT.ALLOCATED(WHGTAG) .AND. IRANK.EQ.IPFGTBL(22))
     &  ALLOCATE(WHGTAG(NGX,NGY))
      IF(.NOT.ALLOCATED(CWHTAG) .AND. IRANK.EQ.IPFGTBL(23))
     &  ALLOCATE(CWHTAG(NGX,NGY))
      IF(.NOT.ALLOCATED(RANGAG) .AND. IRANK.EQ.IPFGTBL(24))
     &  ALLOCATE(RANGAG(NGX,NGY))

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
        CALL ALTAS
        IF (ITEST.GE.1) WRITE(IU06,*) '   SUB. WAMASSI: ALTAS DONE'
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
          ROAIRN(IJ) = ROAIRO(IJ,IG)
          ZIDLNEW(IJ) = ZIDLOLD(IJ,IG)
        ENDDO


!       MODEL OUTPUT INTEGRATED DATA ARE SAVED IN MODULES.
!       --------------------------------------------------
        LLINTPOL=.FALSE.
        IF (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO .OR.
     &      CDTSPT.EQ.CDTPRO .OR. CDTSPS.EQ.CDTPRO ) THEN

          IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN
!!! it is assumed the FL1 still contains the spectra in the relative
!!! frame of reference. It needs to be transformed to the absolute frame 
            CALL GSTATS(1236,0)
            IRA=1
            NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1
!$OMP       PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
            DO JKGLO=IJS(IG),IJL(IG),NPROMA
              KIJS=JKGLO
              KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
              CALL INTPOL (FL1, SL, KIJS, KIJL, IG, IRA)
            ENDDO
!$OMP       END PARALLEL DO
            LLINTPOL=.TRUE.
            CALL GSTATS(1236,1)

            CALL OUTBS (SL, FL3, IJS(IG), IJL(IG), IG, IGL, IU25, IU26)
          ELSE
            NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1
!$OMP       PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
            DO JKGLO=IJS(IG),IJL(IG),NPROMA
              KIJS=JKGLO
              KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
              DO M=1,NFRE
                DO K=1,NANG
                  DO IJ=KIJS,KIJL
                    SL(IJ,K,M) = FL1(IJ,K,M)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
!$OMP       END PARALLEL DO

            CALL OUTBS (SL, FL3, IJS(IG), IJL(IG), IG, IGL, IU25, IU26)
          ENDIF

          IF (ITEST.GE.2.AND.ITESTB.GE.IG) THEN
            WRITE(IU06,*) '   SUB. WAMASSI: OUTPUT DONE IN SUB OUTBS'
          ENDIF
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

        CALL OUTINT (IU06, IU22, IU23, IU32, ITEST, NWFDBREF, LFDBOPEN)

        IF(LMESSPASS) THEN
          CALL GSTATS(752,0)
          CALL MPL_BARRIER(CDSTRING='WAMASSI:')
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

      LSV=(CDTRES.EQ.CDTPRO.OR.CDATEE.EQ.CDTPRO.OR.CDTPRO.EQ.CDATER)

      IF (LSV.OR.LDWRRE) THEN

!*    3.1.4 SAVE RESTART FIELDS.
!           --------------------
        LOUT=(IREST.EQ.1.AND.(CDTPRO.EQ.CDATER .OR. CDTPRO.LE.CDATES))
        IF ( LOUT .OR. LDWRRE ) THEN 

!*          STRESS FIELD IS DUMPED FOR RESTART.
!           -----------------------------------
!           NOTE :
!           STRESS FIELD IS DUMPED FOR RESTART
!           ONLY IF CURRENT TIME CORRESPONDS TO CDATER

          IF((CDTPRO .EQ. CDATER  ).OR.
     &       (CDATER .EQ. '000000000000' ).OR.
     &        LDWRRE) THEN

            IF ( .NOT.LGRIBOUT .OR. LDWRRE ) THEN
              CALL SAVSTRESS(U10OLD, THWOLD, USOLD, TAUW, Z0OLD,
     &                       ROAIRO, ZIDLOLD, CICOVER, CITHICK,
     &                       NSTART,NEND,CDTPRO,CDATEF)
            ENDIF

          ENDIF

!           SAVE SPECTRUM
!           -------------

          IF ( LOUT.AND.LGRIBOUT ) THEN
            CALL OUTWSPEC(FL1)
            CALL GSTATS(1976,0)
            ISTAT = IFLUSHFDBSUBS (NWFDBREF)
            CALL MPL_BARRIER(CDSTRING='WAMASSI:')
            CALL GSTATS(1976,1)
            LLFLUSH = .FALSE.
          ENDIF

          IF ( .NOT.LGRIBOUT .OR. LDWRRE ) THEN
            CALL SAVSPEC(FL1,NSTART,NEND,CDTPRO,CDATEF,CDATER)
          ENDIF

          WRITE(IU06,*) ' '
          WRITE(IU06,*) '  WAVE SPECTRA FIELDS DISPOSED AT........',
     &                  ' CDTPRO  = ', CDTPRO
          WRITE(IU06,*) ' '
          CALL FLUSH (IU06)


!*    3.1.5 UPDATE, WRITE AND SAVE WAMINFO FILE.
!           ------------------------------------

          IF ((LDSTOP.OR.LDWRRE) .AND. 
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

            IU04 =  I_GET_UNIT (IU06,CWI(1:ICPLEN+8) , 'w', 'f', 0)

            CALL WRITSTA (IU04, CDTPRO, CDATEE, IANALPD, IFOREPD,
     &                    IDELWIN, CDATER, CDATES, CBPLTDT, CEPLTDT,
     &                    IASSI, NFCST, ISTAT, CDTCUR)

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
          CALL MPL_BARRIER(CDSTRING='WAMASSI:')
          CALL GSTATS(752,1)
        ENDIF
      ENDIF
      IF (LFDB .AND. LLFLUSH .AND. NWFDBREF.NE.-5 .AND. (
     &  (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO))) THEN
        CALL GSTATS(1976,0)
        ISTAT = IFLUSHFDBSUBS (NWFDBREF)
        CALL MPL_BARRIER(CDSTRING='WAMASSI:')
        CALL GSTATS(1976,1)
        WRITE(IU06,*) '   DB ', NWFDBREF , ' FLUSHED AT ',
     &    CDTPRO, ' FROM WAMASSI. '
        CALL FLUSH (IU06)
      ELSE
        CALL GSTATS(752,0)
        CALL MPL_BARRIER(CDSTRING='WAMASSI:')
        CALL GSTATS(752,1)
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
      IF(ALLOCATED(C4G)) DEALLOCATE(C4G)
      IF(ALLOCATED(BFG)) DEALLOCATE(BFG)
      IF(ALLOCATED(QPG)) DEALLOCATE(QPG)
      IF(ALLOCATED(DEPTHG)) DEALLOCATE(DEPTHG)
      IF(ALLOCATED(HMAXG)) DEALLOCATE(HMAXG)
      IF(ALLOCATED(TMAXG)) DEALLOCATE(TMAXG)
      IF(ALLOCATED(USTOKESG)) DEALLOCATE(USTOKESG)
      IF(ALLOCATED(VSTOKESG)) DEALLOCATE(VSTOKESG)
      IF(ALLOCATED(UCURG)) DEALLOCATE(UCURG)
      IF(ALLOCATED(VCURG)) DEALLOCATE(VCURG)
      IF(ALLOCATED(PHIEPSG)) DEALLOCATE(PHIEPSG)
      IF(ALLOCATED(PHIAWG)) DEALLOCATE(PHIAWG)
      IF(ALLOCATED(TAUOCG)) DEALLOCATE(TAUOCG)
      IF(ALLOCATED(STRNMSG)) DEALLOCATE(STRNMSG)
      IF(ALLOCATED(RHOAG)) DEALLOCATE(RHOAG)
      IF(ALLOCATED(WSTARG)) DEALLOCATE(WSTARG)
      IF(ALLOCATED(CICG)) DEALLOCATE(CICG)
      IF(ALLOCATED(CIHG)) DEALLOCATE(CIHG)
      IF(ALLOCATED(WX1G)) DEALLOCATE(WX1G)
      IF(ALLOCATED(WX2G)) DEALLOCATE(WX2G)
      IF(ALLOCATED(WX3G)) DEALLOCATE(WX3G)
      IF(ALLOCATED(WX4G)) DEALLOCATE(WX4G)
      IF(ALLOCATED(WX5G)) DEALLOCATE(WX5G)

      IF(ALLOCATED(WHGTSG)) DEALLOCATE(WHGTSG)
      IF(ALLOCATED(WDIRSG)) DEALLOCATE(WDIRSG)
      IF(ALLOCATED(WMNFSG)) DEALLOCATE(WMNFSG)
      IF(ALLOCATED(WHGTWG)) DEALLOCATE(WHGTWG)
      IF(ALLOCATED(WDIRWG)) DEALLOCATE(WDIRWG)
      IF(ALLOCATED(WMNFWG)) DEALLOCATE(WMNFWG)

      IF (LHOOK) CALL DR_HOOK('WAMASSI',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE WAMASSI
