      SUBROUTINE FWSEA (F, EWFG, THEW, FMWFG, ETOI, USMO, THMO,         
     &                  IJS, IJL)                                       

! ----------------------------------------------------------------------

!**** *FWSEA* - ANALYSE THE MODEL SPECTRA TO PROVIDE THE WIND-SEA PART  
!               OF THE SPECTRUM, ITS ENERGY, ITS MEAN FREQUENCY.        

!      P.LIONELLO     ECMWF    FEBRUARY 1989                       
!        (MODIFICATION OF A CODE BY P.JANSSEN                           
!         AND P.LIONELLO - SUMMER '87 - )                               
!      J.BIDLOT       ECMWF    JANUARY 1997 : CORRECTION FOR INITIAL
!                                             WAVE DIRECTION NOT BEING 0

!      METHOD.                                                          
!      -------                                                          

!      A PEAK IN THE SPECTRUM IS SOUGHT AROUND THE WIND DIRECTION.      
!      IF A PEAK IS FOUND ,ITS FREQUENCY IS COMPARED TO THE P.M.        
!      FREQUENCY TO ESTABLISH IF IT IS WINDSEA.                         
!      THE WHOLE MODEL SPECTRUM IN THE REGION AROUND THE PEAK IS ASSUMED
!      TO BE WINDSEA. IN THE REMAINIG PART OF THE SPECTRUM THE MINIMUM  
!      OF THE MODEL SPECTRUM AND THE CORRESPONDING JONSWAP SPECTRUM IS  
!      TAKEN AS WINDSEA                                                 
!      WIND SEA ENERGU AND MEAN FREQUENCY ARE YOWPUTED BY THE EXTERNAL  
!      WSMFEN                                                           

!     F        :     MODEL SPECTRUM (FIRST GUESS)                       
!     FSEA     :     WIND-SEA PART OF THE MODEL SPECTRUM  (FIRST GUESS) 
!     EWFG     :     FIRST GUESS WIND-SEA ENERGY                        
!     THEW     :     DIRECTION OF THE WIND-SEA PEAK (FIRST GUESS)       
!     THMO     :     WIND DIRECTION                                     
!     USMO     :     USTAR (FIRST GUESS)                                
!     FMWFG    :     WINDSEA MEAN FREQUENCY (FIRST GUESS)               
!     ETOI     :     TOTAL ENERGY FROM O.I. (ANALYSIS)                  

! ----------------------------------------------------------------------

      USE YOWFRED  , ONLY : FR       ,DELTH    ,TH
      USE YOWMEAN  , ONLY : EMEAN
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------

      DIMENSION F(NINF-1:NSUP,NANG,NFRE)
      REAL,DIMENSION(NINF:NSUP) :: EWFG, THEW, FMWFG, ETOI, THMO, USMO
      DIMENSION FJONS(NANG,NFRE), FSEA(NANG,NFRE)                       

! ----------------------------------------------------------------------

      PARAMETER (CGAM = 0.00608,EGAM = 1.439, ESAO = 2.338)
      PARAMETER (XL11 = 0.0953101 )                                     

!     CGAM       : P.M. FREQUENCY IN THE WAM MODEL                      
!     EGAM,ESAO  : PARAMETERS IN THE WINDSEA SPECTRUM SHAPE             

! ----------------------------------------------------------------------

!     INLINE FUNCTIONS.                                                 
!     -----------------                                                 

!     OVERSHOOT AS FUNCTION OF DIMENSIONLESS FREQUENCY.                 

      GAM(X) = MAX (1.,1+2.70*(1-(CGAM/X)**EGAM))                       

!     PEAK EXTENSION AS FUNCTION OF DIMENSIONLESS ENERGY.               

      SAO(X) = MAX (.02,.02+.16*(1-(CGAM/X)**ESAO))                     

!-----------------------------------------------------------------------

!*   1. PARAMETERS AND INITIALIZATION.                                  
!       ------------------------------                                  

      FR1=FR(1)                                                         
      NWP=0                                                             
      DO IJ=IJS,IJL                                                
        THEW(IJ)  = 999999999.                                         
        EWFG(IJ)  = -999.                                              
        FMWFG(IJ) = -9999.                                             
      ENDDO

! ----------------------------------------------------------------------


!*    2. LOOP ON THE POINTS INSIDE THE BLOCK.                           
!        ------------------------------------                           

      DO IJ=IJS,IJL                                                

!*    2.1 SKIP THE LAND POINTS.                                         
!         ---------------------                                         

        IF (EMEAN(IJ).LE.0.01.OR.ETOI(IJ).LE.0.01) GOTO 2001           

!*    2.2 INITIALIZE SPECTRUM BY ZERO.                                  
!         ----------------------------                                  

        DO M=1,NFRE                                               
          DO K=1,NANG                                            
            FSEA(K,M) = 0.                                           
          ENDDO
        ENDDO

!     2.3 SEARCH A PEAK OF THE SPECTRUM AROUND THE WIND DIRECTION.      
!         --------------------------------------------------------      

        THEWI = THMO(IJ)                                               
        IDW = NINT(MOD(THEWI-TH(1),ZPI)/DELTH)+1
        FSMAX = 0.                                                     
        MP=1                                                           
        KP=1                                                           
        DO II=IDW-1,IDW+1                                         
          K = II                                                      
          IF (K.GT.NANG) K = K-NANG                                   
          IF (K.LE.0) K = K+NANG                                      
          DO M=1,NFRE                                            
            FMODXX = F(IJ,K,M)                                       
            IF (FMODXX.GT.FSMAX) THEN                                
              FSMAX = FMODXX                                        
              MP = M                                                
              KP = K                                                
            ENDIF                                                    
          ENDDO
        ENDDO

!*    2.4 IF A PEAK HAS BEEN FOUND THE PEAK FREQUENCY IS EVALUATED.     
!         ---------------------------------------------------------     

        IF (MP.GT.1 .AND. MP.LT.NFRE) THEN                             
          FMOD=0.                                                     
          FMOD3=0.                                                    
          FMOD1=0.                                                    

!*    2.4.1 THE SPECTRUM IS INTEGRATED AROUND THE PEAK TO OBTAIN        
!           A BETTER ESTIMATE OF THE PEAK IN THE 1-DIM SPECTRUM.        
!           ----------------------------------------------------        

          DO I=-1,1                                              
            KPP=KP+I                                                 
            IF (KPP.EQ.NANG+1) KPP = 1                               
            IF (KPP.EQ.0) KPP = NANG                                 
            FMOD  = F(IJ,KPP,MP  )+FMOD                              
            FMOD3 = F(IJ,KPP,MP+1)+FMOD3                             
            FMOD1 = F(IJ,KPP,MP-1)+FMOD1                             
          ENDDO

!*    2.4.2 CHECK IF THERE IS A PEAK IN THE 1-DIM SPECTRUM.             
!*          IF NOT THEN THE 2D SPECTRUM IS USED TO COMPUTE              
!*          THE PEAK FREQUENCY.                                         
!           -----------------------------------------------             

          IF (FMOD.LE.FMOD3 .OR. FMOD.LE.FMOD1) THEN                  
            FMOD  = F(IJ,KP,MP  )                                    
            FMOD3 = F(IJ,KP,MP+1)                                    
            FMOD1 = F(IJ,KP,MP-1)                                    
          ENDIF                                                      
          X1 = FR(MP)**2-FR(MP+1)**2                                  
          X2 = FR(MP)-FR(MP+1)                                        
          X3 = FMOD-FMOD3                                             
          X4 = FR(MP-1)**2-FR(MP+1)**2                                
          X5 = FR(MP-1)-FR(MP+1)                                      
          X6 = FMOD1-FMOD3                                            
          FPWW = 0.5*(X3*X4-X1*X6)/(X3*X5-X2*X6)                      
          XNU = USMO(IJ)*FPWW/G                                       

!*    2.4.3 EVALUATE PEAK DIRECTION                                     
!*          IF PEAK FREQUENCY IS GE PM FREQUENCY.                       
!           -------------------------------------                       

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

!*    2.4.4 BUILD THE CORRESPONDING WIND-SEA SPECTRUM                   
!*          A PARAMETRIZATION OF THE MODEL SPECTRUM HAVING A "TOBA"     
!*          F**(-4) TAIL IS BUILT .                                     

            SPINT = 0.                                               

!             TOBA CONSTANT ACCORDING TO THE WAM MODEL                  

            ALPHATB = .145                                             
            GAMMA = GAM(XNU)                                         
            SA = SAO(XNU)                                            
            CALL F4SPEC (ALPHATB, FPWW, GAMMA, SA, THE, USMO(IJ),      
     &       SPINT, FJONS)                               
            DO K=1,NANG                                         
              DO M=1,NFRE                                      
                FSEA(K,M) = MIN(F(IJ,K,M),FJONS(K,M))              
              ENDDO
            ENDDO

!*   2.4.5 AROUND THE WIND SEA PEAK THE WHOLE MODEL SPECTRUM IS ASSUMED 
!*         TO BE WINDSEA.                                               
!          ------------------------------------------------------------ 

            IDW = NINT(MOD(THE-TH(1),ZPI)/DELTH)+1
            DO II=IDW-2,IDW+2                                   
              K = MOD(II,NANG)                                      
              IF (K.LE.0)  K=K+NANG                                 
              FSH1 = COS(THE-TH(K))                                 
              IF (FSH1.GT.0.001) THEN                               
                XLL = 0.7*FPWW/(FSH1*FR1)                          
                JP = NINT(1.+LOG(XLL)/XL11)                       
                IF (JP.LE.NFRE) THEN                               
                  JP = MAX (JP,1)                                 
                  DO M=JP,NFRE                               
                    FSEA(K,M) = F(IJ,K,M)                        
                  ENDDO
                ENDIF                                              
              ENDIF                                                 
            ENDDO

!*    2.4.6 THE WIND DIRECTION IS STORED.                               
!           -----------------------------                               

            THEW(IJ) = THE                                           

!*    2.4.7 COMPUTATION OF THE WINDSEA ENERGY WINDSEA MEAN FREQUENCY.   
!           ---------------------------------------------------------   

            CALL WSMFEN (FSEA, EWFG(IJ), FMWFG(IJ), USMO(IJ))        
          ENDIF                                                       
        ENDIF                                                          
 2001   continue
      ENDDO
      IF (ITEST.GT.3)  WRITE(IU06,*)                                    
     &   '          SUB. FUSTAR: WIND SEA FOUND AT ',NWP, ' POINTS'     
                                                                        
      RETURN                                                            
      END SUBROUTINE FWSEA
