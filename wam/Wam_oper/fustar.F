      SUBROUTINE FUSTAR (USMO, USA, USOI,                               
     &                  EWFG, EWOI, EWA, T, FMWA, IJS, IJL)             
                                                                        
! ----------------------------------------------------------------------

!**** *FUSTAR* - ESTIMATE THE ANALYSED FRICTION VELOCITY FROM THE       
!                DURATION AND THE "MEASURED" WINDSEA ENERGY.            
!                ESTIMATE THE ANALYSED WINDSEA MEAN FREQUENCY           
!                FROM THE "MEASURED" WINDSEA ENERGY                     

!     P.LIONELLO     ECMWF       FEBRUARY 1989                          
!         (FROM MODIFICATION OF A CODE BY P.JANSSEN                     
!           AND P.LIONELLO - SUMMER '87 - )                             

!     PURPOSE.                                                          
!     --------                                                          

!        ESTIMATE THE ANALYSED FRICTION VELOCITY FROM THE URATION AND   
!        THE "MEASURED" WINDSEA ENERGY. ESTIMATE THE ANALYSED WINDSEA   
!        MEAN FREQUENCY FROM THE "MEASURED" WINDSEA ENERGY.             

!**   INTERFACE.                                                        
!     ----------                                                        

!        *CALL* *FUSTAR (USMO, USA, USOI,                               
!                        EWFG, EWOI, EWA, T, FMWA, IJS, IJL)*           

!        *USMO*   REAL     FIRST GUESS FRICTION VELOCITY.               
!        *USA*    REAL     ANALISED FRICTION VELOCITY.                  
!        *USOI*   REAL     MEASURED FRICTION VELOCITY.                  
!        *EWFG*   REAL     WINDSEA ENERGY (FIRST GUESS).                
!        *EWOI*   REAL     ESTIMATED WIND-SEA ENERGY FROM MEASUREMENTS. 
!        *EWA*    REAL     ANALISED WIND-SEA ENERGY.                    
!        *T*      REAL     DURATION.                                    
!        *FMWA*   REAL     ANALYSED WINDSEA MEAN FREQUENCY.             
!        *IJS*    INTEGER  FIRST INDEX IN BLOCK.                        
!        *IJL*    INTEGER  LAST  INDEX IN BLOCK.                        

!      METHOD.                                                          
!      -------                                                          

!       THE DURATION OF THE WIND SEA IS USED TO COMPUTE THE FRICTION    
!       VELOCITY BY NEWTONS METHOD FROM THE MODEL DURATION CURVE AND    
!       THE MEASURED ENERGY.  (FIVE ITERATIONS ARE DONE)                
!       THIS VALUE IS USED TO UPDATE THE MODEL IF IT IS                 
!       IN REASONABLE AGREEMENT WITH THE MEASUREMENT FRICTION VELOCITY. 
!       IF THERE IS NO AGREEMENT THE MEASURED WIND SPEED IS USED        
!       TO COMPUTE THE WIND SEA.                                        

! ----------------------------------------------------------------------

      USE YOWPCONS , ONLY : G       ,EGRCRV   ,AGRCRV   ,BGRCRV   ,
     &             AFCRV   ,BFCRV
      USE YOWMPP   , ONLY : NINF    ,NSUP
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      REAL,DIMENSION(NINF:NSUP) :: USMO, USA, USOI, EWOI,      
     &                             EWA, T, FMWA, EWFG          

! ----------------------------------------------------------------------

!     INLINE FUNCTION.                                                  
!     ----------------                                                  

!     DIMENSIONLESS FREQUENCY AS FUNCTION OF DIMENSIONLESS ENERGY       

      YNU(X) = (X/AFCRV)**(1./BFCRV)                                          

! ----------------------------------------------------------------------

!*    1. INITIALIZE OUTPUT FIELDS.                                      
!        -------------------------                                      

      NNW=0                                                             
      DO IJ=IJS,IJL                                                
        USA(IJ)  =   -99.                                              
        EWA(IJ)  =  -999.                                              
        FMWA(IJ) = -9999.                                              
      ENDDO

!*    2. LOOP OVER GRID POINTS.                                         
!        ----------------------                                         

      G2=G*G
      DO IJ=IJS,IJL                                                
        IF (EWFG(IJ).GT.0.01) THEN                                     

!*    2.1  DETERMINE FRICTION VELOCITY FROM DURATION GROWTH CURVE       
!          IF WINDSEA IS PRESENT.                                       
!           ------------------------------------------------------      

          UMODEL = USMO(IJ)                                           
          USTANAL = UMODEL                                            
          EWSR = EWOI(IJ)                                             
          TX = T(IJ)                                                  
          DO ITER=1,5                                            
            TSTAR = G*TX/USTANAL                                     
            XX = AGRCRV*TSTAR**BGRCRV
            IF (XX.GE.2.) THEN                                       
              USTANAL = (EWSR*G2/EGRCRV)**(0.25)                      
            ELSE                                                     
              TT = TANH(XX)                                         
              RKSI = EWSR*G2/USTANAL**4/EGRCRV
              USTANAL = USTANAL*(1.-(RKSI-TT)                       
     &         /(BGRCRV*XX*(1.-TT**2)-4.*TT))               
            ENDIF                                                   
          ENDDO

!*    2.2  IF SIGNS OF (USA-USMO) AND (USOI-USMO) DIFFER, THEN          
!          THE MODEL IS WRONG ABOUT THE CORRECT RATIO BETWEEN           
!          WINDSEA AND SWELL. THEREFORE NO ASSIMILATION IS DONE.        
!          -----------------------------------------------------        

          IF (USOI(IJ).GT.0.) THEN                                    
            SIGNR  = SIGN(1.,(USTANAL-UMODEL))                       
            SIGNM2 = SIGN(1.,(0.8*USOI(IJ)-UMODEL))                  
            SIGNM1 = SIGN(1.,(1.2*USOI(IJ)-UMODEL))                  
            IF (SIGNR.NE.SIGNM1.AND.SIGNR.NE.SIGNM2) THEN            
              NNW = NNW+1                                           
            ELSE                                                     
              USA(IJ) = USTANAL                                     
              EWA(IJ) = EWSR                                        
              FMWA(IJ) = YNU(EWSR*G2/USTANAL**4)*G/USTANAL        
            ENDIF                                                   
          ELSE                                                        
            USA(IJ) = USTANAL                                        
            EWA(IJ) = EWSR                                           
            FMWA(IJ) = YNU(EWSR*G2/USTANAL**4)*G/USTANAL           
          ENDIF                                                       
        ENDIF                                                          
      ENDDO

      WRITE(IU06,*) ' WINDSEA / USTAR ASSIMILATION FAILED AT ',
     &             NNW, ' POINTS'                                       

      RETURN                                                            
      END SUBROUTINE FUSTAR
