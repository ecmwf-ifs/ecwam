      SUBROUTINE UPDATE (F, ETOI, USOI, USMO, THMO, IJS, IJL)           
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C**** *UPDATE* - ANALYSE THE WAVE SPECTRUM, PRODUCING A CONSISTENT      
C****            UPDATE OF WAVE AND WIND FIELD.                         
C                                                                       
C     PIERO LIONELLO      ECMWF     JUNE 1990                           
C                                                                       
C     PURPOSE.                                                          
C     --------                                                          
C                                                                       
C       TO ANALYSE THE WAVE SPECTRUM, PRODUCING A CONSISTENT            
C       UPDATE OF WAVE AND WIND FIELD.                                  
C                                                                       
C**   INTERFACE.                                                        
C     ----------                                                        
C                                                                       
C       *CALL* *UPDATE (F, ETOI, USOI, USMO, THMO, IJS, IJL)*           
C                                                                       
C         *F*       REAL      SPECTRUM ( FIRST GUESS IN INPUT ;         
C                             ANALYSIS IN OUTPUT ).                     
C         *ETOI*    REAL      WAVE ENERGY (FROM O.I.).                  
C         *USOI*    REAL      USTAR ( FROM O.I.).                       
C         *USMO*    REAL      USTAR ( FIRST GUESS IN INPUT ;            
C                             ANALYSIS IN OUTPUT ).                     
C         *THMO*    REAL      WIND DIRECTION.                           
C         *IJS*     INTEGER   FIRST INDEX IN BLOCK.                     
C         *IJL*     INTEGER   LAST  INDEX IN BLOCK.                     
C                                                                       
C     METHOD.                                                           
C     -------                                                           
C                                                                       
C     EXTERNALS.
C     ----------
C     *SEMEAN*
C     *FWSEA*
C     *FDUR*
C     *FUSTAR*
C     *UPWSPEC*
C
C ----------------------------------------------------------------------
C                                                                       
#include "param.h"
C                                                                       
#include "commean.h"
C                                                                       
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "comtest.h"
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
      DIMENSION F(NINF-1:NSUP,NANG,NFRE)
      REAL,DIMENSION(NINF:NSUP) :: ETOI, USOI, USMO, THMO         
      REAL,DIMENSION(NINF:NSUP) :: EWFG, EWOI, EWA,FMWFG, FMWA,
     1                             THEW, TDUR, USA                       
C                                                                       
C       *EWFG*     : WINDSEA ENERGY (FIRST GUESS).                      
C       *EWOI*     : WINDSEA ENERGY (MEASUREMENT).                      
C       *EWA*      : WINDSEA ENERGY (ANALYSIS).                         
C       *ETOI*     : TOTAL ENERGY (FROM O.I ,IDENTICAL WITH ANALYSIS IN 
C                    THIS METHOD ).                                     
C       *FMWFG*    : MEAN FREQUENCY OF THE WINDSEA (FIRST GUESS).       
C       *FMWA*     : MEAN FREQUENCY OF THE WINDSEA (ANALYSIS).          
C       *THEW*     : DIRECTION OF THE WINDSEA PEAK (FIRST GUESS).       
C       *TDUR*     : WINDSEA DURATION.                                  
C       *USA*      : ESTIMATE OF USTAR FROM THE FIRST GUESS ENERGY.     
C                    AND DURATION (I.E. ANALYSED USTAR).                
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    1. FIND THE WINDSEA.                                              
C        -----------------                                              
C                                                                       
      CALL SEMEAN (F, IJS, IJL)
C
      CALL FWSEA (F, EWFG, THEW, FMWFG, ETOI, USMO, THMO, IJS, IJL)     
C                                                                       
C*    2. CORRECT OVERESTIMATE OF FIRST GUESS WINDSEA ENERGY.            
C        ---------------------------------------------------            
C                                                                       
      DO 2001 IJ = IJS,IJL                                              
         EWFG(IJ) = MIN(EWFG(IJ),EMEAN(IJ))                             
 2001 CONTINUE                                                          
                                                                        
      IF (ITEST.GE.3) WRITE (IU06,*) '      SUB. UPDATE: ',             
     1   ' WINDSEA ENERGY FIXED IN SUB. FWSEA'                          
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    3. FIND THE DURATION.                                             
C        ------------------                                             
C                                                                       
      CALL FDUR (EWFG, TDUR, USMO, IJS, IJL)                               
                                                                        
      IF (ITEST.GE.3) WRITE (IU06,*) '      SUB. UPDATE: ',             
     1   ' WINDSEA DURATION FIXED IN SUB. FDUR'                         
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    4. THE RATIO WINDSEA/SWELL IS ASSUMED CORRECT.                    
C        -------------------------------------------                    
C                                                                       
      DO 4001 IJ=IJS,IJL                                                
         IF (EMEAN(IJ).GT.0.) EWOI(IJ) = EWFG(IJ)/EMEAN(IJ)*ETOI(IJ)    
 4001 CONTINUE                                                          
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    5. FIND THE NEW USTAR AND THE NEW WINDSEA MEAN FREQUENCY.         
C        ------------------------------------------------------         
C                                                                       
      CALL FUSTAR (USMO, USA, USOI, EWFG, EWOI, EWA, TDUR, FMWA,
     1             IJS, IJL) 
C                                                                       
      DO 5001 IJ=IJS,IJL                                                
         IF (USA(IJ).GT.0.) USMO(IJ) = USA(IJ)                          
 5001 CONTINUE                                                          
                                                                        
      IF (ITEST.GE.3) WRITE (IU06,*) '      SUB. UPDATE: ',             
     1   ' USTAR COMPUTED IN SUB. FUSTAR'                               
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    6. THE ASSIMILATION IS ADJUSTED ON THE DOMINANT PART              
C*       OF THE SPECTRUM.                                               
C        -------------------------------------------------              
C                                                                       
C        THE FOLLOWING RETURN CODE ARE POSSIBLE :                       
C          THERE IS NO WINDSEA ->                                       
C                 EWFG  =  -999.                                        
C                 FMWFG = -9999.                                        
C          THE RATIO WINDSEA/SWELL IS NOT CONSISTENT WITH USOI ->       
C                 EWA   =  -999.                                        
C                 FMWA  = -9999.                                        
C                                                                       
      DO 6001 IJ=IJS,IJL                                                
         IF (EWA(IJ).GT.0. .AND. EWA(IJ).LT.0.5*ETOI(IJ)) THEN          
            EWFG(IJ)  =  -999.                                          
            FMWFG(IJ) = -9999.                                          
         END IF                                                         
 6001 CONTINUE                                                          
C                                                                       
C     IF THERE WERE WINDSEA BUT THE RATIO WINDSEA SWELL                 
C     WAS WRONG, DO NOT ASSIMILATE.                                     
C                                                                       
      DO 6002 IJ=IJS,IJL                                                
         IF (EWFG(IJ).GT.0. .AND. EWA(IJ).LT.0.) THEN                   
            ETOI(IJ)=-999.                                              
         ENDIF                                                          
 6002 CONTINUE                                                          
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    7. UPDATE THE SPECTRA.                                            
C        -------------------                                            
C                                                                       
      CALL UPWSPEC (F, EMEAN(IJS), ETOI, FMWFG, FMWA, IJS, IJL)               
                                                                        
      IF (ITEST.GE.3) WRITE (IU06,*) '      SUB. UPDATE: ',             
     1   ' NEW SPECTRA COMPUTED IN SUB. UPWSPEC'                         
                                                                        
      RETURN                                                            
      END                                                               
