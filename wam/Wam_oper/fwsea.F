      SUBROUTINE FWSEA (F, EWFG, THEW, FMWFG, ETOI, USMO, THMO,         
     1                  IJS, IJL)                                       
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C**** *FWSEA* - ANALYSE THE MODEL SPECTRA TO PROVIDE THE WIND-SEA PART  
C               OF THE SPECTRUM, ITS ENERGY, ITS MEAN FREQUENCY.        
C                                                                       
C      P.LIONELLO     ECMWF         FEBRUARY 1989                       
C        (MODIFICATION OF A CODE BY P.JANSSEN                           
C         AND P.LIONELLO - SUMMER '87 - )                               
C      J.BIDLOT       ECMWF         JANUARY 1997 : CORRECTION FOR INITIAL
C                                                  WAVE DIRECTION NOT BEING 0
C                                                                       
C      METHOD.                                                          
C      -------                                                          
C                                                                       
C      A PEAK IN THE SPECTRUM IS SOUGHT AROUND THE WIND DIRECTION.      
C      IF A PEAK IS FOUND ,ITS FREQUENCY IS COMPARED TO THE P.M.        
C      FREQUENCY TO ESTABLISH IF IT IS WINDSEA.                         
C      THE WHOLE MODEL SPECTRUM IN THE REGION AROUND THE PEAK IS ASSUMED
C      TO BE WINDSEA. IN THE REMAINIG PART OF THE SPECTRUM THE MINIMUM  
C      OF THE MODEL SPECTRUM AND THE CORRESPONDING JONSWAP SPECTRUM IS  
C      TAKEN AS WINDSEA                                                 
C      WIND SEA ENERGU AND MEAN FREQUENCY ARE COMPUTED BY THE EXTERNAL  
C      WSMFEN                                                           
C                                                                       
C     F        :     MODEL SPECTRUM (FIRST GUESS)                       
C     FSEA     :     WIND-SEA PART OF THE MODEL SPECTRUM  (FIRST GUESS) 
C     EWFG     :     FIRST GUESS WIND-SEA ENERGY                        
C     THEW     :     DIRECTION OF THE WIND-SEA PEAK (FIRST GUESS)       
C     THMO     :     WIND DIRECTION                                     
C     USMO     :     USTAR (FIRST GUESS)                                
C     FMWFG    :     WINDSEA MEAN FREQUENCY (FIRST GUESS)               
C     ETOI     :     TOTAL ENERGY FROM O.I. (ANALYSIS)                  
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
#include "param.h"
C                                                                       
#include "comfred.h"
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
      REAL,DIMENSION(NINF:NSUP) :: EWFG, THEW, FMWFG, ETOI, THMO, USMO
      DIMENSION FJONS(NANG,NFRE), FSEA(NANG,NFRE)                       
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
#include "parcons.h"
C                                                                       
      PARAMETER (CGAM = 0.00608, ESTARM = 955.,                         
     1           A = 6.02*10.**(-5), B = 0.695,                         
     2           AF = 0.000168, BF = - 3.27, EGAM = 1.439, ESAO = 2.338)
      PARAMETER (XL11 = 0.0953101 )                                     
C                                                                       
C     CGAM       : P.M. FREQUENCY IN THE WAM MODEL                      
C     ESTARM,A,B : PARAMETERS OF THE ENERGY GROWTH CURVE                
C     AF,BF      : PARAMETERS IN THE RELATION BETWEEN ENERGY            
C                  AND MEAN FREQUENCY                                   
C     EGAM,ESAO  : PARAMETERS IN THE WINDSEA SPECTRUM SHAPE             
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C     INLINE FUNCTIONS.                                                 
C     -----------------                                                 
C                                                                       
C     OVERSHOOT AS FUNCTION OF DIMENSIONLESS FREQUENCY.                 
C                                                                       
      GAM(X) = MAX (1.,1+2.70*(1-(CGAM/X)**EGAM))                       
C                                                                       
C     PEAK EXTENSION AS FUNCTION OF DIMENSIONLESS ENERGY.               
C                                                                       
      SAO(X) = MAX (.02,.02+.16*(1-(CGAM/X)**ESAO))                     
C                                                                       
C     DIMENSIONLESS ENERGY AS FUNCTION OF DIMENSIONLESS FREQUENCY.      
C                                                                       
      EN(X) = AF*X**BF                                                  
C                                                                       
C     DIMENSIONLESS FREQUENCY AS FUNCTION OF DIMENSIONLESS ENERGY.      
C                                                                       
      YNU(X) = (X/AF)**(1./BF)                                          
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C*   1. PARAMETERS AND INITIALIZATION.                                  
C       ------------------------------                                  
C                                                                       
 1000 CONTINUE                                                          
      FR1=FR(1)                                                         
      NWP=0                                                             
      DO 1001 IJ=IJS,IJL                                                
         THEW(IJ)  = 999999999.                                         
         EWFG(IJ)  = -999.                                              
         FMWFG(IJ) = -9999.                                             
 1001 CONTINUE                                                          
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C                                                                       
C*    2. LOOP ON THE POINTS INSIDE THE BLOCK.                           
C        ------------------------------------                           
C                                                                       
 2000 CONTINUE                                                          
      DO 2001 IJ=IJS,IJL                                                
C                                                                       
C*    2.1 SKIP THE LAND POINTS.                                         
C         ---------------------                                         
C                                                                       
         IF (EMEAN(IJ).LE.0.01.OR.ETOI(IJ).LE.0.01) GOTO 2001           
C                                                                       
C*    2.2 INITIALIZE SPECTRUM BY ZERO.                                  
C         ----------------------------                                  
C                                                                       
         DO 2201 M=1,NFRE                                               
            DO 2202 K=1,NANG                                            
               FSEA(K,M) = 0.                                           
 2202       CONTINUE                                                    
 2201    CONTINUE                                                       
C                                                                       
C     2.3 SEARCH A PEAK OF THE SPECTRUM AROUND THE WIND DIRECTION.      
C         --------------------------------------------------------      
C                                                                       
         THEWI = THMO(IJ)                                               
         IDW = NINT(MOD(THEWI-TH(1),ZPI)/DELTH)+1
         FSMAX = 0.                                                     
         MP=1                                                           
         KP=1                                                           
         DO 2301 II=IDW-1,IDW+1                                         
            K = II                                                      
            IF (K.GT.NANG) K = K-NANG                                   
            IF (K.LE.0) K = K+NANG                                      
            DO 2302 M=1,NFRE                                            
               FMODXX = F(IJ,K,M)                                       
               IF (FMODXX.GT.FSMAX) THEN                                
                  FSMAX = FMODXX                                        
                  MP = M                                                
                  KP = K                                                
               ENDIF                                                    
 2302       CONTINUE                                                    
 2301    CONTINUE                                                       
C                                                                       
C*    2.4 IF A PEAK HAS BEEN FOUND THE PEAK FREQUENCY IS EVALUATED.     
C         ---------------------------------------------------------     
C                                                                       
         IF (MP.GT.1 .AND. MP.LT.NFRE) THEN                             
            FMOD=0.                                                     
            FMOD3=0.                                                    
            FMOD1=0.                                                    
C                                                                       
C*    2.4.1 THE SPECTRUM IS INTEGRATED AROUND THE PEAK TO OBTAIN        
C           A BETTER ESTIMATE OF THE PEAK IN THE 1-DIM SPECTRUM.        
C           ----------------------------------------------------        
C                                                                       
            DO 2411 I=-1,1                                              
               KPP=KP+I                                                 
               IF (KPP.EQ.NANG+1) KPP = 1                               
               IF (KPP.EQ.0) KPP = NANG                                 
               FMOD  = F(IJ,KPP,MP  )+FMOD                              
               FMOD3 = F(IJ,KPP,MP+1)+FMOD3                             
               FMOD1 = F(IJ,KPP,MP-1)+FMOD1                             
 2411       CONTINUE                                                    
C                                                                       
C*    2.4.2 CHECK IF THERE IS A PEAK IN THE 1-DIM SPECTRUM.             
C*          IF NOT THEN THE 2D SPECTRUM IS USED TO COMPUTE              
C*          THE PEAK FREQUENCY.                                         
C           -----------------------------------------------             
C                                                                       
            IF (FMOD.LE.FMOD3 .OR. FMOD.LE.FMOD1) THEN                  
               FMOD  = F(IJ,KP,MP  )                                    
               FMOD3 = F(IJ,KP,MP+1)                                    
               FMOD1 = F(IJ,KP,MP-1)                                    
            END IF                                                      
            X1 = FR(MP)**2-FR(MP+1)**2                                  
            X2 = FR(MP)-FR(MP+1)                                        
            X3 = FMOD-FMOD3                                             
            X4 = FR(MP-1)**2-FR(MP+1)**2                                
            X5 = FR(MP-1)-FR(MP+1)                                      
            X6 = FMOD1-FMOD3                                            
            FPWW = 0.5*(X3*X4-X1*X6)/(X3*X5-X2*X6)                      
            XNU = USMO(IJ)*FPWW/G                                       
C                                                                       
C*    2.4.3 EVALUATE PEAK DIRECTION                                     
C*          IF PEAK FREQUENCY IS GE PM FREQUENCY.                       
C           -------------------------------------                       
C                                                                       
            IF (XNU.GE.CGAM) THEN                                       
               THE = TH(KP)                                             
               KP1 = KP+1                                               
               KM1 = KP-1                                               
               IF (KP1.GT.NANG) KP1=KP1-NANG                            
               IF (KM1.LT.1) KM1=KM1+NANG                               
               FMOD  = F(IJ,KP ,MP)                                     
               FMOD3 = F(IJ,KP1,MP)                                     
               FMOD1 = F(IJ,KM1,MP)                                     
               DDIR = DELTH*(FMOD3-FMOD1)/(2.*FMOD-FMOD1-FMOD3)         
               THE = THE+DDIR/2.                                        
               DTHE = MOD((THE-THEWI),ZPI)                              
               IF (ABS(DTHE).GT.DELTH) GOTO 2001                        
               NWP = NWP+1                                              
C                                                                       
C*    2.4.4 BUILD THE CORRESPONDING WIND-SEA SPECTRUM                   
C*          A PARAMETRIZATION OF THE MODEL SPECTRUM HAVING A "TOBA"     
C*          F**(-4) TAIL IS BUILT .                                     
C                                                                       
               SPINT = 0.                                               
C                                                                       
C             TOBA CONSTANT ACCORDING TO THE WAM MODEL                  
C                                                                       
               ALPHATB = .145                                             
               GAMMA = GAM(XNU)                                         
               SA = SAO(XNU)                                            
               CALL F4SPEC (ALPHATB, FPWW, GAMMA, SA, THE, USMO(IJ),      
     1                      SPINT, FJONS)                               
               DO 2441 K=1,NANG                                         
                  DO 2442 M=1,NFRE                                      
                     FSEA(K,M) = MIN(F(IJ,K,M),FJONS(K,M))              
 2442             CONTINUE                                              
 2441          CONTINUE                                                 
C                                                                       
C*   2.4.5 AROUND THE WIND SEA PEAK THE WHOLE MODEL SPECTRUM IS ASSUMED 
C*         TO BE WINDSEA.                                               
C          ------------------------------------------------------------ 
C                                                                       
               IDW = NINT(MOD(THE-TH(1),ZPI)/DELTH)+1
               DO 2451 II=IDW-2,IDW+2                                   
                  K = MOD(II,NANG)                                      
                  IF (K.LE.0)  K=K+NANG                                 
                  FSH1 = COS(THE-TH(K))                                 
                  IF (FSH1.GT.0.001) THEN                               
                     XLL = 0.7*FPWW/(FSH1*FR1)                          
                     JP = NINT(1.+ALOG(XLL)/XL11)                       
                     IF (JP.LE.NFRE) THEN                               
                        JP = MAX (JP,1)                                 
                        DO 2452 M=JP,NFRE                               
                           FSEA(K,M) = F(IJ,K,M)                        
 2452                   CONTINUE                                        
                     ENDIF                                              
                  ENDIF                                                 
 2451          CONTINUE                                                 
C                                                                       
C*    2.4.6 THE WIND DIRECTION IS STORED.                               
C           -----------------------------                               
C                                                                       
               THEW(IJ) = THE                                           
C                                                                       
C*    2.4.7 COMPUTATION OF THE WINDSEA ENERGY WINDSEA MEAN FREQUENCY.   
C           ---------------------------------------------------------   
C                                                                       
               CALL WSMFEN (FSEA, EWFG(IJ), FMWFG(IJ), USMO(IJ))        
            ENDIF                                                       
         ENDIF                                                          
 2001 CONTINUE                                                          
      IF (ITEST.GT.3)  WRITE(IU06,*)                                    
     1   '          SUB. FUSTAR: WIND SEA FOUND AT ',NWP, ' POINTS'     
                                                                        
      RETURN                                                            
      END                                                               
