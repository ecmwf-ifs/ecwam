      SUBROUTINE WAVEMDL (CBEGDAT, KDELWI, KDURAT, IUGB, IVGB, ILEN,
     1                    BETAG, NLONW, NLATW, LDSTOP, LDWRRE,
     2                    LDRESTARTED, ZDELATM)
C                                                                       
C****  *WAVEMDL* - SUPERVISES EXECUTION OF MAIN MODULES                 
C****              OF THE WAVE MODEL                                    
C                                                                       
C      LIANA ZAMBRESKY    GKSS/ECMWF    OCTOBER 1988                    
C                                                                       
C      MODIFICATION.
C      -------------
C         H. GUNTHER   ECMWF  MARCH    1990                      
C         P. LIONELLO  ECMWF  APRIL    1990  DATA ASSIMILATION   
C                                            MODULE WAMASSI ADDED.
C         J. BIDLOT    ECMWF  FEBRUARY 1996  MESSAGE PASSING.
C         J. DOYLE     ECMWF  OCTOBER  1996  ATMOSPHERIC COUPLING  .
C         J. BIDLOT    ECMWF  FEBRUARY 1997  MESSAGE PASSING.
C         B. HANSEN    ECMWF  MARCH    1997  SIGNAL HANDLING.
C          LDSTOP* - SET .TRUE. IF STOP SIGNAL RECEIVED.
C          LDWRRE* - SET .TRUE. IF RESTART FILE SIGNAL RECEIVED.
C          KDELWI* - WIND INPUT TIME STEP WHEN COUPLED TO THE IFS.
C         J. BIDLOT    ECMWF  April    1997  ADD ZDELATM IN PARAM. LIST 
C                                                                       
C     PURPOSE.                                                          
C     --------                                                          
C                                                                       
C          THIS SUBROUTINE SUPERVISES THE EXECUTION OF                  
C          MAIN MODULES FOR WAM MODEL INITIALIZATION,                   
C          WIND FIELD PREPROCESSING, WAM MODEL EXECUTION,               
C          AND WAVE DATA ASSIMILATION.                                  
C                                                                       
C*    INTERFACE.                                                        
C     ----------                                                        
C                                                                       
C          SEE MAIN MODULES SUB INITMDL, PREWIND, WAMODEL, WAMASSI.     
C                                                                       
C     METHOD.                                                           
C     -------                                                           
C                                                                       
C          THE FIRST TIME WAVEMDL IS CALLED, THE WAM MODEL IS           
C          INITIALIZED. THIS INITIALIZATION INCLUDES GETTING            
C          FROM ECFILE THE INITIAL SEA STATE FILES, FILLING             
C          COMMON BLOCKS DEFINING THE GRID AND SETTING GENERAL          
C          PARAMETERS. IN THE FIRST AND ALL SUBSEQUENT CALLS TO         
C          WAVEMDL  PREWIND REFORMATS THE WINDS INTO THE WAM            
C          MODEL BLOCKED STRUCTURE AND THE WAM MODEL IS EXECUTED.       
C          EACH CALL TO WAMODEL INTEGRATES THE WAVE SPECTRA FORWARD     
C          IN TIME BY ONE INPUT WIND TIME STEP OR PROPAGATION TIME      
C          STEP, WHAT EVER IS GREATER.                                  
C                                                                       
C     EXTERNALS.                                                        
C     ----------                                                        
C                                                                       
C          INITMDL  -  INITIALIZES THE WAM MODEL.                       
C                      GETS RECOVERY FILES OUT OF ECFILE,               
C                      SETS COMMON BLOCKS NECESSARY TO DEFINE           
C                      THE GRID AND BLOCKING STRUCTURE.                 
C                      DEFINES GENERAL PARAMETERS.                      
C                                                                       
C          PREWIND  -  REFORMATS WINDS ON THE GAUSSIAN GRID             
C                      INTO THE WAM MODEL BLOCKED STRUCTURE.            
C                                                                       
C          WAMODEL  -  INTEGRATES THE WAVE SPECTRA FORWARD IN TIME BY   
C                      ONE WIND INPUT TIME STEP OR ONE PROPAGATION      
C                      TIME STEP, WHATEVER IS GREATER.                  
C                                                                       
C          WAMASSI  -  SUPERVISES DATA ASSIMILATION:                    
C                      PREPROCESSES DATA; PRODUCES ANALYSED             
C                      INTEGRATED QUANTITIES BY OPTIMAL INTERPOLATION;  
C                      ANALYSES WAVE SPECTRA ; SAVE ANALYSIS FOR        
C                      OUTPUT AND NEXT TIME STEP MODEL COMPUTATION.     
C                                                                       
C     REFERENCES.                                                       
C     -----------                                                       
C                                                                       
C          NONE                                                         
C                                                                       
C -------------------------------------------------------------------   
      USE WAVE_MODULE
C                                                                       
#include "param.h"
C                                                                       
#include "txtpara.h"
C                                                                       
#include "parameter_wind.h"
C                                                                       
#include "txtpara_wind.h"
C                                                                       
#include "comcoup.h"
C                                                                       
#include "comstat.h"
C                                                                       
#include "txtstat.h"
C                                                                       
#include "comtest.h"
C                                                                       
#include "txttest.h"
C                                                                       
#include "commpp.h"
C                                                                       
#include "txtmpp.h"
C
#include "commap.h"
C
C
#include "comtext.h"
C
C     COUPLING ARRAYS
      REAL BETAG(NLONW,NLATW), ZDELATM(NLATW)
      INTEGER IUGB(ILEN),IVGB(ILEN)
      CHARACTER CBEGDAT*10
C
      LOGICAL FRSTIME, LDSTOP, LDWRRE, LDRESTARTED
C
      DATA FRSTIME, NADV /.TRUE., 0 /
C                                                                       
C     NADV  NUMBER OF PROPAGATION TIME STEPS IN ONE CALL OF WAMODEL.    
C                                                                       
C --------------------------------------------------------------------- 
C                                                                       
C*    1.  THE FIRST CALL TO WAVEMDL PERFORMS INITIALIZATION.            
C         --------------------------------------------------            
C                                                                       
 1000 CONTINUE                                                          
 1010 FORMAT ('  **************************************************')
 1011 FORMAT ('  *                                                *')
 1012 FORMAT ('  *    WAVEMDL: INITMDL STARTS FOR (RE-)START.     *')
 1013 FORMAT ('  *    ========================================    *')
 1014 FORMAT ('  *    TWO-WAY INTERACTION WIND AND WAVES          *')
 1015 FORMAT ('  *    ONE-WAY INTERACTION WIND AND WAVES          *')
 1016 FORMAT ('  * START DATE OF RUN      ', A10,  ' (CDATEA)     *')
 1017 FORMAT ('  * START DATE OF FORECAST ', A10,  ' (CDATEF)     *')
 1018 FORMAT ('  * END DATE OF RUN        ', A10,  ' (CDATEE)     *')
 1019 FORMAT ('  * TOTAL RUN IN HOURS     ', I10,  ' (KDURAT)     *')
 1020 FORMAT ('  * INTERACTION TIME STEP  ', I10,  ' (KDELWI)     *')
C
      IF (FRSTIME) THEN                                                 
         IF(LWCOU) THEN
           IU06=20
           CDATEA(1:10) = CBEGDAT(1:10)
           CDATEF(1:10) = CDATEA(1:10)
           CDATEE(1:10) = CDATEA(1:10)
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
           WRITE (IU06,1017) CDATEF
           WRITE (IU06,1018) CDATEE
           WRITE (IU06,1019) KDURAT/60
           WRITE (IU06,1020) KDELWI
           WRITE (IU06,1011)
           WRITE (IU06,1010)
C
         ENDIF
         CALL INITMDL (KDELWI,NADV)
         FRSTIME = .FALSE.                                              
         IF (ITEST.GE.1) THEN                                           
            WRITE(IU06,*) ' SUB. WAVEMDL: INITMDL DONE'                 
         ENDIF                                                          
C
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
C            RETURN
           ELSE
             LDRESTARTED = .FALSE.
           ENDIF
         ENDIF
C
      ENDIF                                                             
C                                                                       
C --------------------------------------------------------------------  
C                                                                       
C*    2.0  INTEGRATE THE WAVE SPECTRA FORWARD IN TIME.                  
C          -------------------------------------------                  
C                                                                       
 2000 CONTINUE                                                          
C                                                                       
C*    2.1  REFORMAT WINDS FROM GAUSSIAN TO BLOCKED.                     
C          ----------------------------------------                     
C                                                                       
      CALL PREWIND (NC,NR,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     1              IUGB, IVGB, ILEN)
      IF (ITEST.GE.1) THEN                                              
          WRITE(IU06,*) ' SUB. WAVEMDL: PREWIND DONE'                   
      ENDIF                                                             
C                                                                       
C*    2.2  INTEGRATE THE WAVE SPECTRA FORWARD IN TIME.                  
C          -------------------------------------------                  
C                                                                       
      IF (IASSI .EQ. 0 ) THEN                                           
        marstype = 'an'                                                 
      ELSE                                                              
        marstype = 'fg'                                                 
      ENDIF                                                             
                                                                        
      IF (CDTPRO.GT.CDATEF.OR.CDATEA.EQ.CDATEF) THEN

        IF(NTOTENS.EQ.0) THEN
          marstype = 'fc'
        ELSE
          IF(NENSFNB.EQ.0) THEN
            marstype = 'cf'
          ELSE
            marstype = 'pf'
          ENDIF
        ENDIF

      ENDIF                                                             
                                                                        
      CALL WAMODEL (NADV, BETAG, NLONW, NLATW, LDSTOP, LDWRRE)

      IF (ITEST.GE.1) THEN                                              
          WRITE(IU06,*) ' SUB. WAVEMDL: WAMODEL DONE'                   
      ENDIF                                                             
C                                                                       
C*    2.3  DATA ASSIMILATION.                                           
C*         IF REQUESTED AND MODEL IS IN ANALYSIS PERIOD.                
C          ---------------------------------------------                
C                                                                       
      IF ((.NOT.LWCOU .AND. CDTPRO.LE.CDATEF .AND. IASSI.EQ.1) 
     1    .OR. 
     2    (LWCOU .AND. CDTPRO.EQ.CDATEF .AND. IASSI.EQ.1)) THEN                       
         marstype = 'an'                                                

         CALL WAMASSI (BETAG, NLONW, NLATW, LDSTOP, LDWRRE)
C
         IF (ITEST.GE.1) THEN                                           
            WRITE(IU06,*) ' SUB. WAVEMDL: WAMASSI DONE'                 
         ENDIF                                                          
      ENDIF                                                             
C
C PUT ZDELLO IN ZDELATM FOR THE ATMOSPHERIC COUPLING
      IF(LWCOU) THEN
         DO K=1,NLATW
            ZDELATM(K) = ZDELLO(K)
         ENDDO
      ENDIF
C
      RETURN                                                            
      END                                                               
