      SUBROUTINE ERSFILE (IU91, IU92, IG, IGL)                          
                                                                        
C -------------------------------------------------------------------   
C                                                                       
C**** *ERSFILE* - FIILE MANAGEMENT OF ERS_1 FILES.                      
C                                                                       
C     H. GUNTHER               ECMWF         JULY 1991                  
C                                                                       
C     PURPOSE.                                                          
C     --------                                                          
C                                                                       
C       CONTROL INPUT AND OUTPUT FILES OF ERS_1 SPECTRAL DATA.          
C                                                                       
C                                                                       
C*** INTERFACE.                                                         
C    ----------                                                         
C      *CALL* *ERSFILE (IU91, IU92 IG, IGL)*                            
C         *IU91*   - INPUT UNIT OF LCATIOB FILE.                        
C         *IU92*   - OUTPUT UNIT FOR SPECTRA.                           
C         *IG*     - BLOCK NUMBER.                                      
C         *IGL*    - NUMBER OF BLOCKS.                                  
C                                                                       
C     EXTERNALS.                                                        
C     ----------                                                        
C                                                                       
C        NONE.                                                          
C                                                                       
C                                                                       
C     METHOD.                                                           
C     -------                                                           
C                                                                       
C       NONE.                                                           
C                                                                       
C     REFERENCES.                                                       
C     -----------                                                       
C                                                                       
C        NONE.                                                          
C                                                                       
C     MODIFICATIONS.                                                    
C     --------------                                                    
C                                                                       
C        B. HANSEN    *ECMWF*      FEB 93'                              
C           ACCESS TO ECFILE DELETED. INPUT FILES COLYYMMDDhhmm HAVE    
C           TO BE ONLINE ALREADY OTHERWISE NO OUTPUT OF WAVE SPECTRA    
C           CLOSE TO POSITIONS OF ERS1-SAR-SPECTRA.                     
C           WAVE SPECTRA WRITTEN TO UNIT 71. FILES ARE OPEND FOLLOWING  
C           THE KNOWN NAMING KONVENTION COSYYMMDDhhmm. THESE FILES      
C           ARE CLOSED (STATUS=KEEP) AND ARE TO BE SAVED AFTER THE      
C           MODEL RUN.                                                  
C
C     J. BIDLOT       *ECMWF*      JUNE 1996 MIGRATION TO FUJITSU
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
#include "comcoer.h"
C                                                                       
#include "comstat.h"
C                                                                       
#include "comtest.h"
C
#include "commesspass.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
      CHARACTER * 2 CLH       ! HELPER STRING TO STORE HOUR.            
      CHARACTER *13 FNAME                                               
      CHARACTER *10 CL_DATE
      CHARACTER *14 CL_TIME
      LOGICAL FRSTIME                                                   
      DATA FRSTIME / .TRUE./                                            
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    1.  SAVE SPECTRA FILE.                                            
C         ------------------                                            
C
 1000 CONTINUE
      IF (IG.EQ.-9 .AND.IGL.EQ.-9) THEN
C        DO NOT DEAL WITH IU91 AND IU92 IF NOT ON PE 1
         IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
         CLOSE (IU92,STATUS='KEEP')
           IF(ITEST.GT.0) THEN
             WRITE (IU06,*) '      SUB. ERSFILE: UNIT ',
     1                          IU92, '  CLOSED    '
             WRITE (IU06,*) '                                      '
             CALL FLUSH (IU06)
           ENDIF
         ENDIF
         IERS = 0
         RETURN
      ENDIF
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    2. INITIAL DATES IF FIRST CALL.                                   
C        ----------------------------                                   
C                                                                       
 2000 CONTINUE                                                          
      IF (FRSTIME) THEN                                                 
         IERS = 0                                                       
         IDELERS = 21600                                                
         CDTERS  = CDATEA                                               
         CLH=CDATEA(7:8)
         IF (CLH.EQ.'09'.OR.CLH.EQ.'15'.OR.
     1       CLH.EQ.'21'.OR.CLH.EQ.'03') THEN
            CALL INCDATE (CDTERS, -10800)
         ENDIF
         FRSTIME = .FALSE.                                              
         IF (ITEST.GE.3) THEN
           WRITE(IU06,*) " SUB: ERSFILE (for the first time) "//
     1       " CDATEA=", CDATEA, " CDTERS=", CDTERS, " CDATEF=", CDATEF
           CALL FLUSH (IU06)
         ENDIF
      ENDIF                                                             
      IF (CDTERS.LT.CDTPRO) THEN                                        
 1001   CONTINUE                                                        
        CALL INCDATE (CDTERS, IDELERS)                                  
        IF (CDTERS.LT.CDTPRO) GOTO 1001                                 
      ENDIF                                                             
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    3. FETCH LOCATION FILE.                                           
C        --------------------                                           
C                                                                       
 3000 CONTINUE                                                          

C     DO NOT DEAL WITH IU91 AND IU92 IF NOT ON PE 1
      IF(LMESSPASS.AND.(IRANK.NE.1)) RETURN

      IERS = 0                                                          
      IF (CDTERS.LE.CDATEF.OR.MARSTYPE.EQ.'an'.OR.MARSTYPE.EQ.'fg') THEN
         WRITE(FNAME,'(''COL'',A10)') CDTERS                            
         OPEN (UNIT=IU91, FILE=FNAME, FORM='UNFORMATTED', STATUS='OLD', 
     1         IOSTAT=IERR)                                             
         IF (IERR.NE.0) THEN                                            
           IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
             IF(ITEST.GT.0) THEN
               WRITE (IU06,*) '                                     '
               WRITE (IU06,*) '      SUB. ERSFILE: FILE ',
     1                       FNAME, ' NOT OPENED'
               WRITE (IU06,*) '                                     ' 
             ENDIF
           ENDIF
           RETURN                                                      
         ELSE                                                           
            READ(IU91)CL_DATE
           IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
             IF(ITEST.GT.0) THEN
               WRITE (IU06,*) '                                     '
               WRITE (IU06,*) '      SUB. ERSFILE: FILE ',
     1                         FNAME, ' OPENED    '
             ENDIF
           ENDIF
C                                                                       
           WRITE(FNAME,'(''COS'',A10)') CDTERS                         
           OPEN (UNIT=IU92, FILE=FNAME(1:IECF_LEN(FNAME))//MARSTYPE,   
     1           FORM='UNFORMATTED', STATUS='UNKNOWN')                 
             IF(ITEST.GT.0) THEN
               WRITE (IU06,*) '      SUB. ERSFILE: FILE ',
     1                         FNAME, ' OPENED'
               CALL FLUSH (IU06)
             ENDIF
           WRITE(IU92) CDTPRO                                    
         ENDIF
 3001    CONTINUE                                                       
         IERS = IERS + 1                                                
         IF (IERS.GT.NERS) THEN                                         
           IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++'   
            WRITE (IU06,*) ' +                                     +'   
            WRITE (IU06,*) ' +     WARNING ERROR SUB ERSFILE.      +'   
            WRITE (IU06,*) ' +     ==========================      +'   
            WRITE (IU06,*) ' + NUMBER OF OUTPUT POINTS             +'   
            WRITE (IU06,*) ' + EXCEEDS DIMENSION NERS = ', NERS         
            WRITE (IU06,*) ' +                                     +'   
            WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++'   
           ENDIF
            IERS = NERS                                                 
            RETURN                                                      
         ENDIF                                                          
         READ (IU91, END=3002, ERR=3002) CL_TIME, DUMLE, DUMPE,
     1        DUMLM, DUMPM, IGERS(IERS), IJERS(IERS)                    
         IF (IGERS(IERS).EQ.0) THEN                                     
            IERS = IERS-1                                               
            GOTO 3001                                                   
         ENDIF                                                          
         IF (IERS.GT.1) THEN                                            
            DO 3003 I=1,IERS-1                                          
               IF (IGERS(I).EQ.IGERS(IERS) .AND.                        
     1             IJERS(I).EQ.IJERS(IERS)) THEN                        
                  IERS = IERS-1                                         
                  GOTO 3001                                             
               ENDIF                                                    
 3003       CONTINUE                                                    
         ENDIF                                                          
         GOTO 3001                                                      
 3002    CONTINUE                                                       
         IERS = IERS - 1                                                
         CLOSE (IU91, STATUS='KEEP')                                    
         IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            IF(ITEST.GT.0) THEN
              WRITE(IU06,*) '      SUB. ERSFILE: NEW OUTPUT',
     1                   ' POINTS FOR ERS_1 IERS = ', IERS
              CALL FLUSH (IU06)
            ENDIF
         ENDIF
      ENDIF                                                             
C
      RETURN                                                            
      END                                                               
