      SUBROUTINE WAVEMDL (CBEGDAT, KDELWI, KDURAT, IUGB, IVGB, ILEN,
     &                    BETAG, NLONW, NLATW, LDSTOP, LDWRRE,
     &                    LDRESTARTED, ZDELATM)

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

      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWCOUT  , ONLY : CASS     ,NASS
      USE YOWCOUP  , ONLY : LWCOU    ,LWCOU2W
      USE YOWSTAT  , ONLY : MARSTYPE ,CDATEA   ,CDATEE   ,CDATEF   ,
     &            CDTPRO   ,IASSI    ,MFDBSF
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWMAP   , ONLY : ZDELLO
      USE YOWTEXT  , ONLY : LRESTARTED
      USE YOWSPEC, ONLY : U10OLD   ,THWOLD   ,USOLD    ,Z0OLD    ,
     &            TAUW

! --------------------------------------------------------------------- 

!     COUPLING ARRAYS
      REAL BETAG(NLONW,NLATW), ZDELATM(NLATW)
      INTEGER IUGB(ILEN),IVGB(ILEN)
      CHARACTER CBEGDAT*12, CDTASS*12

      LOGICAL FRSTIME, LDSTOP, LDWRRE, LDRESTARTED

      DATA FRSTIME, NADV /.TRUE., 0 /

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
 1016 FORMAT ('  * START DATE OF RUN      ', A12,  ' (CDATEA)   *')
 1018 FORMAT ('  * END DATE OF RUN        ', A12,  ' (CDATEE)   *')
 1019 FORMAT ('  * TOTAL RUN IN HOURS     ', I10,  ' (KDURAT)     *')
 1020 FORMAT ('  * INTERACTION TIME STEP  ', I10,  ' (KDELWI)     *')

      IF (FRSTIME) THEN                                                 
!       if you change the -5 value, do the same in wamodel
        MFDBSF = -5
        IF(LWCOU) THEN
          IU06=20
          IL = LEN_TRIM(CBEGDAT)
          IF(IL.NE.12) THEN
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
        CALL INITMDL (KDELWI,NADV)
        FRSTIME = .FALSE.                                              
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

      ENDIF

! --------------------------------------------------------------------  

!*    2.0  INTEGRATE THE WAVE SPECTRA FORWARD IN TIME.                  
!          -------------------------------------------                  


!*    2.1  REFORMAT WINDS FROM GAUSSIAN TO BLOCKED.                     
!          ----------------------------------------                     

      CALL PREWIND (U10OLD,THWOLD,USOLD,TAUW,Z0OLD,IUGB, IVGB, ILEN)

      IF (ITEST.GE.1) THEN
        WRITE(IU06,*) ' SUB. WAVEMDL: PREWIND DONE'                   
      ENDIF

!*    2.2  INTEGRATE THE WAVE SPECTRA FORWARD IN TIME.                  
!          -------------------------------------------                  

      CALL SETMARSTYPE

      CALL WAMODEL (NADV, BETAG, NLONW, NLATW, LDSTOP, LDWRRE)

      IF (ITEST.GE.1) THEN
        WRITE(IU06,*) ' SUB. WAVEMDL: WAMODEL DONE'                   
      ENDIF

!*    2.3  DATA ASSIMILATION.                                           
!*         IF REQUESTED AND MODEL IS IN ANALYSIS PERIOD.                
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
            CALL WAMASSI (BETAG, NLONW, NLATW, LDSTOP, LDWRRE)
            IF (ITEST.GE.1) THEN
              WRITE(IU06,*) ' SUB. WAVEMDL : WAMASSI DONE'
            ENDIF
          ENDIF
        ELSEIF ((.NOT.LWCOU .AND. CDTPRO.LE.CDATEF ) 
     &      .OR. 
     &      (LWCOU .AND. CDTPRO.EQ.CDATEF)) THEN
          CALL WAMASSI (BETAG, NLONW, NLATW, LDSTOP, LDWRRE)
          IF (ITEST.GE.1) THEN
            WRITE(IU06,*) ' SUB. WAVEMDL: WAMASSI DONE'                 
          ENDIF
        ENDIF

      ENDIF

! PUT ZDELLO IN ZDELATM FOR THE ATMOSPHERIC COUPLING
      IF(LWCOU) THEN
        DO K=1,NLATW
          ZDELATM(K) = ZDELLO(K)
        ENDDO
      ENDIF
      IF ( MFDBSF .GT. 0 .AND. CDATEE .EQ. CDTPRO) THEN
        ISTAT = ICLOSEFDB (MFDBSF)
        IF (ITEST.GE.1) THEN
          WRITE(IU06,*) ' SUB. WAVEMDL: END OF RUN: CLOSE FDB'
          CALL FLUSH(IU06)
        ENDIF
      ENDIF

      RETURN                                                            
      END SUBROUTINE WAVEMDL
