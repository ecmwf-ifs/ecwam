      SUBROUTINE UPWSPEC (F, ETFG, ETA, FMWFG, FMWA, IJS, IJL)           
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C**** *UPWSPEC* - MODIFY THE SPECTRUM BY STRETCHING AND SCALING.         
C                                                                       
C     P.LIONELLO      ECMWF         APRIL 1990                          
C     J.BIDLOT        ECMWF         JUNE 1997 UPSPEC >> UPWSPEC 
C                                                                       
C     PURPOSE.                                                          
C     --------                                                          
C                                                                       
C         TO MODIFY THE SPECTRUM BY STRETCHING AND SCALING.             
C                                                                       
C**   INTERFACE.                                                        
C     ----------                                                        
C                                                                       
C        *CALL* *UPWSPEC (F, ETFG, ETA, FMWFG, FMWA, IJS, IJL)*          
C                                                                       
C          *F*      REAL     INPUT IS THE OLD SWELL SPECTRUM.           
C                            OUTPUT IS  THE ANALISED SPECTRUM.          
C          *ETFG*   REAL     FIRST GUESS TOTAL ENERGY.                  
C          *ETA*    REAL     ANALISED TOTAL ENERGY.                     
C          *FMWFG   REAL     MEAN FREQUENCY OF THE WINDSEA FROM THE     
C                            FIRST GUESS SPECTRUM.                      
C          *FMWA*   REAL     ANALYSED MEAN FREQUENCY OF THE WINDSEA.    
C          *IJS*    INTEGER  FIRST INDEX IN BLOCK.                      
C          *IJL*    INTEGER  LAST  INDEX IN BLOCK.                      
C                                                                       
C     METHOD.                                                           
C     -------                                                           
C       A NEW SPECTRUM IN THE FORM                                      
C          F   (IPOINT,F,K) =  A F   (IPOINT,BF,K)                      
C           NEW                   OLD                                   
C       IS BUILD.                                                       
C                                                                       
C       IF THERE IS MAINLY SWELL THE SPECTRUM IS UPDATED USING THE      
C       AVERAGE STEEPNESS CRITERIUM. TO CONSERVE EXACTLY THE            
C       STEEPNESS AND CHANGE THE ENERGY THE CONSTANT A AND B            
C       MUST BE GIVEN BY                                                
C                                                                       
C                A = (ETA/ETFG)**1.25  B = (ETA/ETFG)**.25              
C                                                                       
C       A SMALL CORRECTION , WHICH IS SUGGESTED BY THE MODEL            
C       DECAY CURVE ,IS ACTUALLY INTRODUCED ACCORDING TO                
C                                                                       
C            DELTA = 1 - .006 * ( HNEW - HOLD )                         
C                                                                       
C                A = DELTA * (ETA/ETFG)**1.25                           
C                B = DELTA * (ETA/ETFG)**.25                            
C                                                                       
C       IF THERE IS MAINLY WINDSEA, THE STRETCHING CONSTANT B IS        
C       COMPUTED TO PRODUCE IN THE ANALYSED  SPECTRUM THE ANALYSED      
C       MEAN FREQUENCY DERIVED BY THE MODEL GROWTH CURVE                
C       ( AS CARRIED OUT IN THE SUBROUTINE FUSTAR ). IN THIS CASE       
C                                                                       
C                  B = FMWFG/FMWA                                       
C                  A = (ETA/ETFG)*B                                     
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
#include "param.h"
C                                                                       
#include "comfred.h"
C                                                                       
#include "commpp.h"
C
#include "txtmpp.h"
C
C ----------------------------------------------------------------------
C                                                                       
      DIMENSION F(NINF-1:NSUP,NANG,NFRE)
      REAL,DIMENSION(NINF:NSUP) :: ETA, FMWFG, FMWA
ccc note ETFG corresponds to EMEAN which is still in a common block in UPDATE
ccc therefore it has to have a different definition as the other arrays
      REAL,DIMENSION(IJS:IJL) :: ETFG
C
      PARAMETER (XL11 = 0.0953101)                                      
      DIMENSION FTEMP(NANG,NFRE)                                        
C                                                                       
C          *XL11*   REAL     ALOG(1.1)                                  
C          *FTEMP*  REAL     TEMPORAR STORAGE OF THE SPECTRUM.          
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
      FR1=FR(1)                                                         
C                                                                       
C*    1. LOOP OVER GRID POINTS.                                         
C        ---------------------                                          
C                                                                       
 1000 CONTINUE                                                          
      DO 1001 IJ=IJS,IJL                                                
C                                                                       
C*    1.1 SKIP LAND POINTS AND POINTS WHERE THERE ARE NO RELIABLE DATA. 
C         ------------------------------------------------------------- 
C                                                                       
         IF (ETFG(IJ).LE.0.001.OR.ETA(IJ).LE.0.001) GOTO 1001           
         FNEW = FMWA(IJ)                                                
         HNEW = 4.*SQRT(ETA(IJ))                                        
         HOLD = 4.*SQRT(ETFG(IJ))                                       
         FOLD = FMWFG(IJ)                                               
C                                                                       
C*    1.2 COMPUTE SCALING AND STRETCHING FACTORS.                       
C         ---------------------------------------                       
C                                                                       
         IF (FOLD.LE.0.) THEN                                           
C                                                                       
C*    1.2.1 THE SPECTRUM IS MAINLY SWELL (THE WINDSEA MEAN FREQUENCY    
C*          OF THE FIRST GUESS SPECTRUM IS NEGATIVE).                   
C           --------------------------------------------------------    
C                                                                       
            XDELTA = 1.-0.006*(HNEW-HOLD)                               
            XR = XDELTA*(HNEW/HOLD)**2.5                                
            XB = XDELTA* SQRT(HNEW/HOLD)                                
         ELSE                                                           
C                                                                       
C*    1.2.2 THE SPECTRUM IS MAINLY WINDSEA.                             
C           -------------------------------                             
C                                                                       
            XR = (HNEW/HOLD)**2*FOLD/FNEW                               
            XB = FOLD/FNEW                                              
         END IF                                                         
C                                                                       
C*    1.3 LOOP OVER NEW FREQUENCIES.                                    
C         --------------------------                                    
C                                                                       
         DO 1301 M=1,NFRE                                               
            FU = FR(M)*XB                                               
            M1 = IFIX(LOG(FU/FR1)/XL11)+1                               
            IF (M1.GE.NFRE.OR.M1.LT.1) THEN                             
C                                                                       
C*    1.3.1 LOOP OVER DIRECTIONS FOR FREQUENCIES GETTING ENERGY FORM    
C*          FREQUENCIES OUT OF RANGE.                                   
C           ---------------------------------------------------------   
C                                                                       
               DO 1311 K=1,NANG                                         
                  FTEMP(K,M) = 0.                                       
 1311          CONTINUE                                                 
            ELSE                                                        
C                                                                       
C*    1.3.2 LOOP OVER DIRECTIONS FOR FREQUENCIES GETTING ENERGY.        
C           ----------------------------------------------------        
               M2 = M1+1                                                
               DFRE = FR(M1)*.1                                         
               DO 1321 K=1,NANG                                         
                  F1 = F(IJ,K,M1)                                       
                  F2 = F(IJ,K,M2)                                       
                  DE = (F2-F1)/DFRE * (FU-FR(M1))                       
                  FTEMP(K,M) = XR*MAX(0.,F1+DE)                         
 1321          CONTINUE                                                 
            END IF                                                      
C                                                                       
C*    BRANCH BACK TO 1.3 FOR NEXT FREQUENCY.                            
C                                                                       
 1301    CONTINUE                                                       
C                                                                       
C*    1.4 THE UPDATED SPECTRUM IS STORED.                               
C         -------------------------------                               
C                                                                       
         DO 1401 K=1,NANG                                               
            DO 1402 M=1,NFRE                                            
               F(IJ,K,M) = FTEMP(K,M)                                   
 1402       CONTINUE                                                    
 1401    CONTINUE                                                       
C                                                                       
C*    BRANCH BACK TO 1. FOR NEXT GRIDPOINT.                             
C                                                                       
 1001 CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
