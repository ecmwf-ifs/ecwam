      SUBROUTINE WAVEMDL (CBEGDAT, KDELWI, KDURAT,
     &                    NFIELDS, NGPTOTG, NC, NR, FCRANGE,
     &                    ISEC1, ISEC2, RMISS, FIELDS,
     &                    NWVFIELDS, WVFLDG,
     &                    NLONW, NLATW, LDSTOP, LDWRRE,
     &                    LDRESTARTED, ZDELATM, KQGAUSS,
     &                    LDWCOUNORMS, MASK_IN, MASK_OUT,
     &                    NFDBREF,
     &                    FRSTIME, NADV)

!****  *WAVEMDL* - SUPERVISES EXECUTION OF MAIN MODULES                 
!****              OF THE WAVE MODEL                                    

!      LIANA ZAMBRESKY    GKSS/ECMWF    OCTOBER 1988                    

!      MODIFICATION.
!      -------------
!         H. GUNTHER   ECMWF  MARCH    1990                      
!         P. LIONELLO  ECMWF  APRIL    1990  DATA ASSIMILATION   
!                                            MODULE WAMASSI ADDED.
!         J. BIDLOT    ECMWF  FEBRUARY 1996  MESSAGE PASSING.
!         J. DOYLE     ECMWF  OCTOBER  1996  ATMOSPHERIC COUPLING  .
!         J. BIDLOT    ECMWF  FEBRUARY 1997  MESSAGE PASSING.
!         B. HANSEN    ECMWF  MARCH    1997  SIGNAL HANDLING.
!          LDSTOP* - SET .TRUE. IF STOP SIGNAL RECEIVED.
!          LDWRRE* - SET .TRUE. IF RESTART FILE SIGNAL RECEIVED.
!          KDELWI* - WIND INPUT TIME STEP WHEN COUPLED TO THE IFS.
!         J. BIDLOT    ECMWF  April    1997  ADD ZDELATM IN PARAM. LIST 
!         S. ABDALLA   ECMWF  OCTOBER  2001  INCLUSION OF AIR DENSITY & Zi/L
!                                            GENERALIZE THE INTERFACE WITH
!                                            THE ATMOSPHERIC MODEL
!         G.Mozdzynski ECMWF  January  2005  OPTIMISE COUPLING COMMS
!         J BIDLOT     ECMWF  August 2006  PASS RMISS TO SPECIFY MISSING DATA.
!         M. Drusch    ECMWF  Sep 2007  Re-initialize through FRSTIME and NADV
!                                            

!     PURPOSE.                                                          
!     --------                                                          

!          THIS SUBROUTINE SUPERVISES THE EXECUTION OF                  
!          MAIN MODULES FOR WAM MODEL INITIALIZATION,                   
!          WIND FIELD PREPROCESSING, WAM MODEL EXECUTION,               
!          AND WAVE DATA ASSIMILATION.                                  

!*    INTERFACE.                                                        
!     ----------                                                        

!          SEE MAIN MODULES SUB INITMDL, PREWIND, WAMODEL, WAMASSI.     

!     METHOD.                                                           
!     -------                                                           

!          THE FIRST TIME WAVEMDL IS CALLED, THE WAM MODEL IS           
!          INITIALIZED. THIS INITIALIZATION INCLUDES GETTING            
!          FROM ECFILE THE INITIAL SEA STATE FILES, FILLING             
!          COMMON BLOCKS DEFINING THE GRID AND SETTING GENERAL          
!          PARAMETERS. IN THE FIRST AND ALL SUBSEQUENT CALLS TO         
!          WAVEMDL  PREWIND REFORMATS THE WINDS INTO THE WAM            
!          MODEL BLOCKED STRUCTURE AND THE WAM MODEL IS EXECUTED.       
!          EACH CALL TO WAMODEL INTEGRATES THE WAVE SPECTRA FORWARD     
!          IN TIME BY ONE INPUT WIND TIME STEP OR PROPAGATION TIME      
!          STEP, WHAT EVER IS GREATER.                                  

!     EXTERNALS.                                                        
!     ----------                                                        

!          INITMDL  -  INITIALIZES THE WAM MODEL.                       
!                      GETS RECOVERY FILES OUT OF ECFILE,               
!                      SETS COMMON BLOCKS NECESSARY TO DEFINE           
!                      THE GRID AND BLOCKING STRUCTURE.                 
!                      DEFINES GENERAL PARAMETERS.                      

!          PREWIND  -  REFORMATS WINDS ON THE GAUSSIAN GRID             
!                      INTO THE WAM MODEL BLOCKED STRUCTURE.            

!          WAMODEL  -  INTEGRATES THE WAVE SPECTRA FORWARD IN TIME BY   
!                      ONE WIND INPUT TIME STEP OR ONE PROPAGATION      
!                      TIME STEP, WHATEVER IS GREATER.                  

!          WAMASSI  -  SUPERVISES DATA ASSIMILATION:                    
!                      PREPROCESSES DATA; PRODUCES ANALYSED             
!                      INTEGRATED QUANTITIES BY OPTIMAL INTERPOLATION;  
!                      ANALYSES WAVE SPECTRA ; SAVE ANALYSIS FOR        
!                      OUTPUT AND NEXT TIME STEP MODEL COMPUTATION.     

!     REFERENCES.                                                       
!     -----------                                                       

!          NONE                                                         

! -------------------------------------------------------------------   

      USE YOWCOUT  , ONLY : CASS     ,NASS
      USE YOWCOUP  , ONLY : LWCOU    ,LWCOU2W  ,LWCOUNORMS,KCOUSTEP,
     &         LMASK_OUT_NOT_SET, LMASK_TASK_STR
      USE YOWGRIBHD, ONLY : NKSEK1   ,IFS_KSEC1,LNEWLVTP
      USE YOWCURR  , ONLY : U        ,V        ,CDTCUR   ,IDELCUR  ,
     &             LLCHKCFL,LLCHKCFLA
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWICE   , ONLY : ICEMASK
      USE YOWMAP   , ONLY : ZDELLO   ,IQGAUSS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NBLO     ,LL1D
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSTAT  , ONLY : MARSTYPE ,CDATEA   ,CDATEE   ,CDATEF   ,
     &            CDTPRO   ,IDELPRO  ,IDELWI   ,IDELWO   ,IREFRA   ,
     &            IASSI    ,NWFDBREF ,LFDBOPEN ,LSARINV
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWTEXT  , ONLY : LRESTARTED
      USE YOWSPEC, ONLY : U10OLD   ,THWOLD   ,USOLD    ,Z0OLD    ,
     &            TAUW   ,ROAIRO   ,ZIDLOLD  ,USNEW    ,Z0NEW    ,
     &            NSTART ,NEND
      USE YOWWIND  , ONLY : CDAWIFL  ,IUNITW
      USE FDBSUBS_MOD
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! --------------------------------------------------------------------- 


!     COUPLING ARRAYS
      REAL WVFLDG(NLONW,NLATW,NWVFIELDS), ZDELATM(NLATW)
      REAL, DIMENSION(NWVFIELDS) :: FAVG,FMIN,FMAX
      REAL, DIMENSION(NWVFIELDS) :: DEFVAL
      REAL  ZHOOK_HANDLE

      INTEGER MASK_IN(NGPTOTG)
      INTEGER MASK_OUT(NLONW,NLATW)

! Mod for OPENMP
      INTEGER JKGLO,KIJS,KIJL,NPROMA
! End Mod for OPENMP

      INTEGER NFIELDS, NGPTOTG, NC, NR
      INTEGER, DIMENSION(NWVFIELDS) :: NTOT 
      REAL FCRANGE
      INTEGER ISEC1(NKSEK1), ISEC2(22+NR)
      REAL FIELDS(NGPTOTG,NFIELDS)
      REAL, ALLOCATABLE :: WVBLOCK(:,:)
      
      CHARACTER CBEGDAT*14, CDTASS*14
      CHARACTER*14 :: CDATEIN, CDTNEWCUR
      CHARACTER*12 :: FLABEL(NWVFIELDS)
      CHARACTER*24 FILNM

      LOGICAL :: LDSTOP, LDWRRE, LDRESTARTED, LDWCOUNORMS
      LOGICAL :: LLCURRENT
      LOGICAL :: FRSTIME
      LOGICAL :: L1STCALL
      LOGICAL :: LLGLOBAL_WVFLDG
      INTEGER :: NADV
      INTEGER :: NFDBREF

!     NADV  NUMBER OF PROPAGATION TIME STEPS IN ONE CALL OF WAMODEL.    
! --------------------------------------------------------------------- 

!*    1.  THE FIRST CALL TO WAVEMDL PERFORMS INITIALIZATION.            
!         --------------------------------------------------            

 1010 FORMAT ('  **************************************************')
 1011 FORMAT ('  *                                                *')
 1012 FORMAT ('  *    WAVEMDL: INITMDL STARTS FOR (RE-)START.     *')
 1013 FORMAT ('  *    ========================================    *')
 1014 FORMAT ('  *    TWO-WAY INTERACTION WIND AND WAVES          *')
 1015 FORMAT ('  *    ONE-WAY INTERACTION WIND AND WAVES          *')
 1016 FORMAT ('  * START DATE OF RUN      ', A14,  ' (CDATEA)   *')
 1018 FORMAT ('  * END DATE OF RUN        ', A14,  ' (CDATEE)   *')
 1019 FORMAT ('  * TOTAL RUN IN HOURS     ', I10,  ' (KDURAT)     *')
 1020 FORMAT ('  * INTERACTION TIME STEP  ', I10,  ' (KDELWI)     *')

!     KEEP GRIB HEADER KSEC1 FROM THE IFS

      IF (LHOOK) CALL DR_HOOK('WAVEMDL',0,ZHOOK_HANDLE)

      ZMISS=RMISS

      IF(LWCOU) THEN
        IF (ITEST.GE.1) THEN
            WRITE (IU06,*) ' WAVEMDL: KSEC1 FROM IFS'
        ENDIF
        DO IC=1,NKSEK1
          IFS_KSEC1(IC)=ISEC1(IC)
          IF (ITEST.GE.1) THEN
            WRITE (IU06,*) 'KSEC1(',IC,')', ISEC1(IC)
          ENDIF
        ENDDO
      ENDIF

      KCOUSTEP=KDELWI


      IF (FRSTIME) THEN                                                 
        WRITE (IU06,*) ' WAVEMDL:'
        WRITE (IU06,*) ' FRSTIME', FRSTIME
        CALL FLUSH(IU06)

        LMASK_OUT_NOT_SET=.TRUE.
        LMASK_TASK_STR=.TRUE.

        IF(LWCOU) THEN
          IU06=20
          IL = LEN_TRIM(CBEGDAT)
          IF(IL.NE.14) THEN
            WRITE (IU06,*) ' NON-Y2K COMPLIANT DATE IN CALLING WAVEMDL'
            WRITE (IU06,*) ' CBEGDAT =  ',CBEGDAT
            WRITE (IU06,*) ' PROGRAM WILL ABORT '
            CALL FLUSH(IU06)
            CALL ABORT1
          ENDIF
          CDATEA = CBEGDAT
          CDATEF = CDATEA
          CDATEE = CDATEA
          IDURAT = KDURAT*60
          CALL INCDATE (CDATEE,IDURAT)
          WRITE (IU06,1010)
          WRITE (IU06,1011)
          WRITE (IU06,1012)
          WRITE (IU06,1013)
          IF (LWCOU2W) THEN
            WRITE (IU06,1014)
          ELSE
            WRITE (IU06,1015)
          ENDIF
          WRITE (IU06,1011)
          WRITE (IU06,1016) CDATEA
          WRITE (IU06,1018) CDATEE
          WRITE (IU06,1019) KDURAT/60
          WRITE (IU06,1020) KDELWI
          WRITE (IU06,1011)
          WRITE (IU06,1010)

        ENDIF

!       INQUIRE IF IUNITW IS ALREADY OPEN THEN CLOSE IT
        IF(IUNITW.NE.0) THEN
          CALL PBCLOSE(IUNITW,KRET)
          IUNITW=0
        ENDIF

        CALL INITMDL (NADV)

        LLCHKCFL=.FALSE.

        KQGAUSS=IQGAUSS
!       if new levtype is used in coupled run than it assumed
!       that WAM and IFS will be sharing the same stream, therefore
!       the FDB should be open on the IFS side.
        IF(LWCOU .AND. LNEWLVTP) THEN
          IF(NFDBREF .GT. 0) THEN
            NWFDBREF=NFDBREF 
            LFDBOPEN=.TRUE.
          ELSE
            NWFDBREF = -5
            LFDBOPEN=.FALSE.
          ENDIF          
        ELSE
!       if you change the -5 value, do the same in wamodel
          NWFDBREF = -5
          LFDBOPEN=.FALSE.
        ENDIF

        FRSTIME = .FALSE.                                              
        L1STCALL = .TRUE.

        IF (ITEST.GE.1) THEN
          WRITE(IU06,*) ' SUB. WAVEMDL: INITMDL DONE'                 
        ENDIF

        IF(LWCOU) THEN
          IF(NLONW.NE.NGX.OR.NLATW.NE.NGY) THEN
            WRITE (IU06,*) ' *********************************'
            WRITE (IU06,*) ' *                               *'
            WRITE (IU06,*) ' * PROBLEM IN WAVEMDL..........  *'
            WRITE (IU06,*) ' * PROBLEM WITH NLONW AND NLATW  *'
            WRITE (IU06,*) ' * NOT EQUAL TO NGX   AND NGY  : *'
            WRITE (IU06,*) ' * ============================= *'
            WRITE (IU06,*) ' *                               *'
            WRITE (IU06,*) ' * NLONW=',NLONW
            WRITE (IU06,*) ' * NLATW=',NLATW
            WRITE (IU06,*) ' * NGX=',NGX
            WRITE (IU06,*) ' * NGY=',NGY
            WRITE (IU06,*) ' *                               *'
            WRITE (IU06,*) ' *                               *'
            WRITE (IU06,*) ' *********************************'
            CALL ABORT1
          ENDIF
          IF (LRESTARTED)  THEN
            LDRESTARTED = .TRUE.
!            RETURN
          ELSE
            LDRESTARTED = .FALSE.
          ENDIF
        ENDIF

      ELSE
!       CHECK IF THE INTERACTION TIME STEP HAS NOT CHANGED
!       OTHERWISE RESET NADV
        IF(IDELWI.NE.KDELWI .AND. LWCOU) THEN
          IDELWO = KDELWI
          IDELWI = KDELWI
          NADV = IDELWI/IDELPRO
          CDAWIFL=CDTPRO
          CALL INCDATE(CDAWIFL,IDELWI)
          WRITE (IU06,*) ' ***********************************'
          WRITE (IU06,*) ' *                                 *'
          WRITE (IU06,*) ' * IN WAVEMDL :                    *'
          WRITE (IU06,*) ' * INTERACTION TIME STEP WAS RESET *' 
          WRITE (IU06,*) ' * KDELWI = ', KDELWI 
          WRITE (IU06,*) ' *                                 *'
          WRITE (IU06,*) ' ***********************************'
          CALL FLUSH(IU06)
          IF(IDELPRO*NADV.NE.IDELWI) THEN
            WRITE (IU06,*) ' ****************************************'
            WRITE (IU06,*) ' *                                      *'
            WRITE (IU06,*) ' * PROBLEM IN WAVEMDL :                 *'
            WRITE (IU06,*) ' * THE NEW INTERACTION TIME STEP IS NOT *'
            WRITE (IU06,*) ' * A MULTIPLE OF IDELPRO  !!!           *'
            WRITE (IU06,*) ' * KDELWI = ', KDELWI 
            WRITE (IU06,*) ' * IDELPRO = ', IDELPRO 
            WRITE (IU06,*) ' *                                      *'
            WRITE (IU06,*) ' ****************************************'
            CALL ABORT1
          ENDIF
        ENDIF

!       KEEP ATMOSPHERIC MODEL INFORMED ABOUT OUR GRID
        KQGAUSS=IQGAUSS

        IF(LWCOU .AND. LNEWLVTP .AND. NFDBREF .GT. 0) THEN
          NWFDBREF=NFDBREF 
          LFDBOPEN=.TRUE.
        ENDIF

        L1STCALL = .FALSE.
      ENDIF

! --------------------------------------------------------------------  

!*    2.   INTEGRATE THE WAVE SPECTRA FORWARD IN TIME.                  
!          -------------------------------------------                  


!     2.0  GET NEW CURRENTS IF IT IS REQUIRED.
!          ----------------------------------

      CALL GSTATS(1962,0)

      IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN

        CDTNEWCUR=CDTCUR
        CALL INCDATE(CDTNEWCUR,1800*IDELCUR)

        IF(CDTPRO.GE.CDTNEWCUR) THEN
          IF(.NOT.ALLOCATED(U)) ALLOCATE(U(NINF-1:NSUP,NBLO))
          IF(.NOT.ALLOCATED(V)) ALLOCATE(V(NINF-1:NSUP,NBLO))

          LLCURRENT=.FALSE.

          CALL INCDATE(CDTCUR,3600*IDELCUR)

          FILNM = 'currents'
          LIU   = LEN_TRIM(FILNM)
          FILNM=FILNM(1:LIU)
          INQUIRE(FILE=FILNM,EXIST=LLCURRENT)

          IREAD=1

          IF (LLCURRENT) THEN
            WRITE(IU06,*) ' '
            WRITE(IU06,*) '   SUB. WAVEMDL: GETTING OCEAN CURRENTS'
            WRITE(IU06,*) '   FOR DATE ',CDTCUR
            CALL FLUSH(IU06)

            CALL CURRENT2WAM (FILNM,IREAD,CDATEIN)

            IF(CDATEIN.NE.CDTCUR) THEN
              WRITE (IU06,*) ' **************************************'
              WRITE (IU06,*) ' *                                    *'
              WRITE (IU06,*) ' * PROBLEM IN WAVEMDL :               *'
              WRITE (IU06,*) ' * THE REQUESTED DATE FOR THE CURRENTS*'
              WRITE (IU06,*) ' * DOES NOT CORRESPOND TO THE DECODED *'
              WRITE (IU06,*) ' * DATE !!!!                          *'
              WRITE (IU06,*) ' * CDTCUR =',CDTCUR 
              WRITE (IU06,*) ' * CDATEIN=',CDATEIN
              WRITE (IU06,*) ' *                                    *'
              WRITE (IU06,*) ' **************************************'
              CALL ABORT1
            ENDIF
          ELSE
            U=0.
            V=0.
            WRITE(IU06,*)' '
            WRITE(IU06,*)'    ****************************'
            WRITE(IU06,*)'     FILE ',FILNM,' NOT FOUND '
            WRITE(IU06,*)'     CURRENTS ARE SET TO 0. '
            WRITE(IU06,*)'    ****************************'
            WRITE(IU06,*)' '
            CALL FLUSH(IU06)
          ENDIF


!         COMPUTE REFRACTION TERMS
          CALL PROPDOT
          IF (ITEST.GE.2)
     &     WRITE(IU06,*) '    SUB. WAVEMDL: REFRACTION TERMS',
     &     'INITIALIZED '
           CALL FLUSH(IU06)

          LLCHKCFLA=.TRUE.
        ELSE
          LLCHKCFLA=.FALSE.
        ENDIF
      ELSE
        LLCHKCFLA=.FALSE.
      ENDIF

      CALL GSTATS(1962,1)


!*    2.1  REFORMAT WINDS FROM GAUSSIAN TO BLOCKED.                     
!          ----------------------------------------                     

      CALL PREWIND (U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     &              ROAIRO, ZIDLOLD, ICEMASK,
     &              NFIELDS, NGPTOTG, NC, NR, FCRANGE,
     &              ISEC1, ISEC2, FIELDS, MASK_IN)

!     WRITE(IU06,*) ' WAVEMDL: AFTER PREWIND, NGPTOTG=', NGPTOTG,
!    & ' SUM(MASK_IN)=', SUM(MASK_IN)

      IF (ITEST.GE.1) THEN
        WRITE(IU06,*) ' SUB. WAVEMDL: PREWIND DONE'                   
      ENDIF


!*    2.2  INTEGRATE THE WAVE SPECTRA FORWARD IN TIME.                  
!          -------------------------------------------                  

      CALL SETMARSTYPE

      CALL WAMODEL (NADV, LDSTOP, LDWRRE, L1STCALL)

      IF (ITEST.GE.1) THEN
        WRITE(IU06,*) ' SUB. WAVEMDL: WAMODEL DONE'
      ENDIF

!*    2.3  DATA ASSIMILATION AND/OR SAR INVERSION.
!*         IF REQUESTED AND MODEL IS IN ANALYSIS PERIOD.
!          THE DATA ASSIMILATION SOFTWARE IS NOT AVAILABLE
!          FOR GENERAL DISSIMINATION !
!          ---------------------------------------------

      IF (IASSI.EQ.1) THEN

        MARSTYPE = 'an'

!       UPDATE ANALYSIS TIME
        DO J=1,NASS
          IF(CDTPRO.EQ.CASS(J)) THEN
            CDTASS=CDTPRO
            EXIT
          ENDIF
        ENDDO

        IF (NASS.GT.0 ) THEN
          IF ( CDTPRO.EQ.CDTASS ) THEN
            CALL WAMASSI (LDSTOP, LDWRRE)
            IF (ITEST.GE.1) THEN
              WRITE(IU06,*) ' SUB. WAVEMDL : WAMASSI DONE'
            ENDIF
          ENDIF
        ELSEIF ((.NOT.LWCOU .AND. CDTPRO.LE.CDATEF ) 
     &      .OR. 
     &      (LWCOU .AND. CDTPRO.EQ.CDATEF)) THEN
          CALL WAMASSI (LDSTOP, LDWRRE)
          IF (ITEST.GE.1) THEN
            WRITE(IU06,*) ' SUB. WAVEMDL: WAMASSI DONE'
          ENDIF
        ENDIF
      ELSE IF (IASSI.NE.1 .AND. LSARINV) THEN

!!!!
!!!!! we may have to introduce a list for sar inversion time !!!!
!       UPDATE ANALYSIS TIME (which is now the inversion time)
        DO J=1,NASS
          IF(CDTPRO.EQ.CASS(J)) THEN
            CDTASS=CDTPRO
            EXIT
          ENDIF
        ENDDO

        IF (NASS.GT.0 ) THEN
          IF ( CDTPRO.EQ.CDTASS ) THEN
            CALL SARINVERT
            IF(ITEST.GE.1)
     &         WRITE(IU06,*) '   SUB. WAVEMDL: SARINVERT DONE'
          ENDIF
        ELSEIF ((.NOT.LWCOU .AND. CDTPRO.LE.CDATEF ) 
     &      .OR. 
     &      (LWCOU .AND. CDTPRO.EQ.CDATEF)) THEN
          CALL SARINVERT
          IF(ITEST.GE.1) WRITE(IU06,*) '   SUB. WAVEMDL: SARINVERT done'
        ENDIF
      ENDIF

!----------------------------------------------------------------------

!     3. PREPARE FIELDS THAT ARE RETURNED TO IFS (if coupled).
!        -----------------------------------------------------

      IF(LWCOU) THEN

        IG=1

        ALLOCATE(WVBLOCK(IJS(IG):IJL(IG),NWVFIELDS))

!       GRID LAYOUT
        DO K=1,NLATW
          ZDELATM(K) = ZDELLO(K)
        ENDDO

!       FIELDS TO BE PASSED TO THE ATMOSPHERIC MODEL ARE:

!       1. CHARNOCK FIELD
        FLABEL(1)='Charnock'
        DEFVAL(1)=0.018  ! DEFAULT VALUE FOR GRID POINTS NOT COVERED BY
                         ! THE WAVE MODEL ICE FREE SEA POINTS.

        NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1
        CALL GSTATS(1443,0)
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=IJS(IG),IJL(IG),NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
          CALL OUTBETA (KIJS, KIJL, IG, USOLD, Z0OLD,
     &                  DEFVAL(1),WVBLOCK(KIJS,1))
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1443,1)

!!!!!!
!       if other fields are needed for IFS add  them here
!!!!!!


!       GRIDDED FIELDS ARE NEEDED FOR IFS, GATHERING OF THE NECESSARY 
!       INFORMATION

        CALL MPFLDTOIFS(IG, IJS(IG), IJL(IG), NWVFIELDS, WVBLOCK,
     &                  WVFLDG, DEFVAL, MASK_OUT, LLGLOBAL_WVFLDG)


        DEALLOCATE(WVBLOCK)


!       COMPUTATION OF THE NORMS OF OUTPUT FIELDS
        WRITE(IU06,*) ' '
        IF(LLGLOBAL_WVFLDG)THEN
          WRITE(IU06,*) ' GLOBAL NORM OF FIELDS RETURNED TO IFS :'

          CALL GSTATS(1443,0)
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(JF,J,I)
          DO IFLD=1,NWVFIELDS
            JF=IFLD
            NTOT(JF)=NLATW*NLONW
            FAVG(JF)=0.
            FMIN(JF)=WVFLDG(1,1,JF)
            FMAX(JF)=WVFLDG(1,1,JF)
            DO J=1,NLATW
              DO I=1,NLONW
                FAVG(JF)=FAVG(JF)+WVFLDG(I,J,JF)
                FMIN(JF)=MIN(FMIN(JF),WVFLDG(I,J,JF))
                FMAX(JF)=MAX(FMAX(JF),WVFLDG(I,J,JF))
              ENDDO
            ENDDO
            FAVG(JF)=FAVG(JF)/NTOT(JF)
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1443,1)

          DO IFLD=1,NWVFIELDS
            WRITE(IU06,*) FLABEL(IFLD),
     &                    FAVG(IFLD),FMIN(IFLD),FMAX(IFLD),NTOT(IFLD)
          ENDDO

        ELSE
          WRITE(IU06,*) ' Local NORM OF FIELDS RETURNED TO IFS :'

          CALL GSTATS(1443,0)
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(JF,J,I)
          DO IFLD=1,NWVFIELDS
            JF=IFLD
            NTOT(JF)=0
            FAVG(JF)=0.
            FMIN(JF)=0.
            FMAX(JF)=0.
            LOOP_J: DO J=1,NLATW
              DO I=1,NLONW
                IF(MASK_OUT(I,J).EQ.1) THEN
                  FMIN(JF)=WVFLDG(I,J,JF)
                  FMAX(JF)=WVFLDG(I,J,JF)
                  EXIT LOOP_J
                ENDIF
              ENDDO
            ENDDO LOOP_J
            DO J=1,NLATW
              DO I=1,NLONW
                IF(MASK_OUT(I,J).EQ.1) THEN
                  NTOT(JF)=NTOT(JF)+1
                  FAVG(JF)=FAVG(JF)+WVFLDG(I,J,JF)
                  FMIN(JF)=MIN(FMIN(JF),WVFLDG(I,J,JF))
                  FMAX(JF)=MAX(FMAX(JF),WVFLDG(I,J,JF))
                ENDIF
              ENDDO
            ENDDO
            FAVG(JF)=FAVG(JF)/MAX(NTOT(JF),1)
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1443,1)

          DO IFLD=1,NWVFIELDS
            WRITE(IU06,*) FLABEL(IFLD),
     &                    FAVG(IFLD),FMIN(IFLD),FMAX(IFLD),NTOT(IFLD),
     &                    IRANK, NPROC, LL1D
          ENDDO
        ENDIF
        WRITE(IU06,*) ' '
        CALL FLUSH(IU06)




!       TELL ATMOS MODEL THE WAM REQUIREMENT FOR GLOBAL NORMS
        LDWCOUNORMS=LWCOUNORMS
      ENDIF


!     4. END OF RUN ?
!        -----------
      IF ( NWFDBREF .GT. 0 .AND. CDATEE .EQ. CDTPRO) THEN
        IF(LWCOU .AND. LNEWLVTP) THEN
!         the fdb should be closed by IFS
          IF (ITEST.GE.1) THEN
            WRITE(IU06,*) ' SUB. WAVEMDL: END OF WAVE MODEL RUN'
            CALL FLUSH(IU06)
          ENDIF
        ELSE
        CALL GSTATS(1787,0)
        ISTAT = ICLOSEFDBSUBS (NWFDBREF)
        CALL GSTATS(1787,1)
          IF (ITEST.GE.1) THEN
            WRITE(IU06,*) ' SUB. WAVEMDL: END OF RUN: CLOSE FDB'
            CALL FLUSH(IU06)
          ENDIF
        ENDIF
      ENDIF

      IF (LHOOK) CALL DR_HOOK('WAVEMDL',1,ZHOOK_HANDLE)

      RETURN                                                            
      END SUBROUTINE WAVEMDL
