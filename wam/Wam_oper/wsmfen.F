      SUBROUTINE WSMFEN (FSEA, EW, FM, USTT, DPT)

!---------------------------------------------------------------------- 

!**** *WSMFEN* - COMPUTES MEAN FREQUENCY AND ENERGY OF THE WINDSEA      
!****            SPECTRUM.                                              

!     P.LIONELLO     ECMWF       APRIL 1990                             

!     PURPOSE.                                                          
!     --------                                                          

!       COMPUTES MEAN FREQUENCY AND ENERGY OF THE WINDSEA SPECTRUM.     

!**   INTERFACE.                                                        
!     ----------                                                        

!       *CALL* *WSMFEN (FSEA, EW, FM, USTT, DPT)*                            

!        *FSEA*   REAL     WINDSEA PART OF THE WAVE SPECTRUM.           
!        *EW*     REAL     WINDSEA ENERGY .                             
!        *FM*     REAL     MEAN FREQUENCY OF THE WINDSEA.               
!        *USTT*   REAL     FRICTION VELOCITY.                           
!        *DPT*    REAL     DEPTH.                           

!     METHOD.                                                           
!     -------                                                           

!       ENERGY AND MEAN FREQUENCY ARE FIRST COMPUTED BY INTEGRATION     
!       OF THE WINDSEA PART OF THE SPECTRUM:                            
!       THIS IMPLIES AN UNDERESTIMATION BOTH OF ENERGY AND MEAN         
!       FREQUENCY. THE UNDERESTIMATED VALUES ARE USED IN THE GROWTH     
!       CURVE PRODUCING OVERSTIMATES OF BOTH ENERGY AND MEAN FREQUENCY. 
!       THE AVERAGE OF THE TWO ESTIMATES IS TAKEN TO PROVIDE A BEST     
!       ESTIMATE.                                                       

! ----------------------------------------------------------------------

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH   ,FRATIO
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,EPSMIN   ,AFCRV   ,BFCRV     ,
     &            ESH      ,ASH     ,BSH
      USE YOWSTAT  , ONLY : ISHALLO


! ----------------------------------------------------------------------

      REAL :: FSEA(NANG,NFRE)

! ----------------------------------------------------------------------

!     INLINE FUNCTIONS.
!     -----------------

!     ENERGY AS FUNCTION OF MEAN FREQUENCY.

      EN(X)=AFCRV*X**BFCRV

!     MEAN FREQUENCY AS FUNCTION OF ENERGY.

      YNU(X)=(X/AFCRV)**(1./BFCRV)

!     DIMENSIONLESS ENERGY LIMIT AS FUNCTION OF DIMENSIONLESS DEPTH 

      EMAX(X) = ESH*TANH(ASH*X**BSH)

! ----------------------------------------------------------------------

!*    1. INTEGRATING THE WINDSEA PART OF THE SPECTRUM.                  
!*       (UNDERESTIMATE IS IMPLIED)                                     
!        ---------------------------------------------                  


!*    1.1 COMPUTATION OF THE WIND-SEA ENERGY.                           
!         -----------------------------------                           

      EW = EPSMIN                                                      
      DO M=1,NFRE                                                  
        TEMP = 0.                                                      
        DO K=1,NANG                                               
          TEMP = TEMP+FSEA(K,M)                                       
        ENDDO
        EW = EW+TEMP*DFIM(M)                                           
      ENDDO
      DELT25 = FR(NFRE)/3.*DELTH                                        
      EW = EW+DELT25*TEMP                                               

!*    1.2 COMPUTATION OF THE MEAN FREQUENCY.                            
!         ----------------------------------                            

      FM = EPSMIN                                                     
      DO M=1,NFRE                                                  
        SPFB = 0.                                                      
        DO K=1,NANG                                               
          SPFB = SPFB+FSEA(K,M)                                       
        ENDDO
        SPFB = SPFB*DELTH                                              
        FM = FM+SPFB*FR(M)**2.*(FRATIO-1./FRATIO)/2.
      ENDDO
      FM = FM+SPFB*FR(NFRE)**2./3.                                      
      FM = FM/EW                                                        

! ----------------------------------------------------------------------

!*    2. ESTIMATES FROM MODEL RELATIONS.                                
!        -------------------------------                                

!*    2.1 ENERGY IS DERIVED FROM UNDERESTIMATED MEAN FREQUENCY.         
!         -----------------------------------------------------         

      XTEMP = FM*USTT/G                                                 
      IF (XTEMP.LE.0.) WRITE(*,*) FM,USTT,EW                            

      ESTAR=EN(FM*USTT/G)
      IF (ISHALLO.NE.1) THEN
         DSTAR = DPT*G/USTT**2
         ESTAR = MIN(EMAX(DSTAR),ESTAR)
      ENDIF
      SPINTDI = USTT**4/G**2 * ESTAR 

!*    2.2 MEAN FREQUENCY IS DERIVED FROM THE UNDERESTIMATED ENERGY.     
!         ---------------------------------------------------------     

      ESTAR=EW*G**2/USTT**4                            
      IF (ISHALLO.NE.1) THEN
         DSTAR = DPT*G/USTT**2
         ESTAR = MIN(EMAX(DSTAR),ESTAR)
      ENDIF

      FREQDI = G/USTT * YNU(ESTAR)                            

! ----------------------------------------------------------------------

!*    3. FINAL ESTIMATE.                                                
!        ---------------                                                

!*    3.1 AVERAGING THE ENERGY ESTIMATES.                               
!         -------------------------------                               

      EW = (EW+SPINTDI)*0.5                                             

!*    3.2 AVERAGING THE MEAN FREQUENCY ESTIMATES.                       
!         ---------------------------------------                       

      FM = (FM+FREQDI)*0.5                                              
                                                                        
      RETURN                                                            
      END SUBROUTINE WSMFEN
