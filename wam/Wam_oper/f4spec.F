      SUBROUTINE F4SPEC (ALFA, FM, GAMMA, SA, THETAQ, USTAR, SPINT, SPT)

!---------------------------------------------------------------------- 

!**** *F4SPEC* - PRODUCE A WINDSEA SPECTRUM FOR ASSIMILATION.           

!      P. LIONELLO    ECMWF        FEBRUARY 1990                        
!      (DERIVED FROM JSPEC OF PETER JANSSEN  1987.)                     

!**   INTERFACE.                                                        
!     ----------                                                        

!       *CALL* *F4SPEC (ALFA, FM, GAMMA, SA, THETAQ, USTAR, SPINT,SPT)* 

!          *ALFA*  - TOBAS' CONSTANT                                    
!          *FM*    - PEAK FREQUENCY                                     
!          *GAMMA* - OVERSHOOT PARAMETER                                
!          *SA*    - WIDTH PARAMETER                                    
!          *THETAQ*- MEAN WINDSEA DIRECTION                             
!          *USTAR* - FRICTION VELOCITY                                  
!          *SPINT* - TOTAL ENERGY                                       
!          *SPT*   - 2-D SPECTRUM                                       

!      METHOD.                                                          
!      -------                                                          

!      THE ONE DIMENSIONAL PARAMETRIC SPECTRUM IS EVALUATED AND         
!      DISTRIBUTED AROUND THE MEAN SEA DIRECTION.                       
!      THE TOTAL ENERGY IS COMPUTED.                                    

!      EXTERNALS.                                                       
!      ----------                                                       

!         *SPR*  -  SUB TO COMPUTE COS**2 SPREADING FUNCTION.           

! --------------------------------------------------------------------- 

      USE YOWFRED  , ONLY : FR       ,DELTH    ,TH       ,FRATIO
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI

! ----------------------------------------------------------------------

      DIMENSION SPT(NANG,NFRE)                                          
      DIMENSION DELFR(NFRE), ET(NFRE), ST(NANG)                         

!      *SPT*        REAL    2-D SPECTRUM.                               
!      *DELFR*      REAL    FREQUENCY INCREMENTS.                       
!      *ET*         REAL    1-D JONSWAP SPECTRUM.                       
!      *ST*         REAL    SPREADING FUNCTION.                         

!---------------------------------------------------------------------- 

      ML = NFRE
      KL = NANG

!*    1. COMPUTE FREQUENCY INTERVALLS.                                  
!        -----------------------------                                  

      DO M=1,ML                                                    
        DELFR(M) = FR(M)*(FRATIO-1.)                                           
      ENDDO

!---------------------------------------------------------------------- 

!*    2. COMPUTE TOBA SPECTRUM.                                         
!        ----------------------                                         

      FAKAL = ALFA*G*USTAR/(ZPI**3)                                       
      DO M=1,ML                                                    
        FRH = FR(M)                                                    
        EARG = (.5*((FRH-FM)/(SA*FM))**2)                              
        FJON = GAMMA**EXP(-EARG)                                       
        FMPF = 4./6.*(FM/FRH)**6                                       
        FJONH = EXP(-FMPF)                                             
        ET(M) = FAKAL*FJONH*FJON/FRH**4                                  
      ENDDO

!---------------------------------------------------------------------- 

!*  3. COMPUTATION OF SPREADING FUNCTION.                               
!      ----------------------------------                               

      CALL SPR (KL, NANG, THETAQ, TH, ST)                               

!---------------------------------------------------------------------- 

!*  4. COMPUTATION OF 2-D SPECTRUM.                                     
!      ---------------------------                                      

      DO M=1,ML                                                    
        DO K=1,KL                                                 
          SPT(K,M)=ET(M)*ST(K)                                        
        ENDDO
      ENDDO

! --------------------------------------------------------------------- 

!*  5. COMPUTE TOTAL ENERGY.                                            
!      --------------------                                             


!*  5.1 INTEGRATE OVER DIRECTIONS.                                      
!       --------------------------                                      

      KL2 = KL-2                                                        
      SPINT = 0.                                                        
      DO M=1,ML                                                    
        ET(M) = SPT(1,M)                                               
        DO K=2,KL2,2                                              
          ET(M) = ET(M)+4.*SPT(K,M)+2.*SPT(K+1,M)                     
        ENDDO
        ET(M) = (ET(M)+4.*SPT(KL,M)+SPT(1,M))*DELTH/3.                 
      ENDDO

!*  5.2 INTEGRATE 1-D SPECTRUM OVER FREQUENCIES.                        
!       ----------------------------------------                        

      SPINT = ET(1)*DELFR(1)                                            
      ML2 = ML-2                                                        
      DO M=2,ML2,2                                                 
        SPINT = SPINT+(4.*ET(M)*DELFR(M)+2.*ET(M+1)*DELFR(M+1))        
      ENDDO
      IF (2*(ML/2).EQ.ML) THEN                                          
        SPINT = (SPINT+4.*ET(ML)*DELFR(ML))/3.                         
      ELSE                                                              
        SPINT = (SPINT+4.*ET(ML-1)*DELFR(ML-1)                         
     &   + ET(ML)*DELFR(ML))/3.                                   
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END SUBROUTINE F4SPEC
