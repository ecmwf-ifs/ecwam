      SUBROUTINE UPWSPEC (F, ETFG, ETA, FMWFG, FMWA, IJS, IJL)

!-----------------------------------------------------------------------

!**** *UPWSPEC* - MODIFY THE SPECTRUM BY STRETCHING AND SCALING.         

!     P.LIONELLO      ECMWF         APRIL 1990                          
!     J.BIDLOT        ECMWF         JUNE 1997 UPSPEC >> UPWSPEC 

!     PURPOSE.                                                          
!     --------                                                          

!         TO MODIFY THE SPECTRUM BY STRETCHING AND SCALING.             

!**   INTERFACE.                                                        
!     ----------                                                        

!        *CALL* *UPWSPEC (F, ETFG, ETA, FMWFG, FMWA, IJS, IJL)*          

!          *F*      REAL     INPUT IS THE OLD SWELL SPECTRUM.           
!                            OUTPUT IS  THE ANALISED SPECTRUM.          
!          *ETFG*   REAL     FIRST GUESS TOTAL ENERGY.                  
!          *ETA*    REAL     ANALISED TOTAL ENERGY.                     
!          *FMWFG   REAL     MEAN FREQUENCY OF THE WINDSEA FROM THE     
!                            FIRST GUESS SPECTRUM.                      
!          *FMWA*   REAL     ANALYSED MEAN FREQUENCY OF THE WINDSEA.    
!          *IJS*    INTEGER  FIRST INDEX IN BLOCK.                      
!          *IJL*    INTEGER  LAST  INDEX IN BLOCK.                      

!     METHOD.                                                           
!     -------                                                           
!       A NEW SPECTRUM IN THE FORM                                      
!          F   (IPOINT,F,K) =  A F   (IPOINT,BF,K)                      
!           NEW                   OLD                                   
!       IS BUILD.                                                       

!       IF THERE IS MAINLY SWELL THE SPECTRUM IS UPDATED USING THE      
!       AVERAGE STEEPNESS CRITERIUM. TO CONSERVE EXACTLY THE            
!       STEEPNESS AND CHANGE THE ENERGY THE CONSTANT A AND B            
!       MUST BE GIVEN BY                                                

!                A = (ETA/ETFG)**1.25  B = (ETA/ETFG)**.25              

!       A SMALL CORRECTION , WHICH IS SUGGESTED BY THE MODEL            
!       DECAY CURVE ,IS ACTUALLY INTRODUCED ACCORDING TO                

!            DELTA = 1 - .006 * ( HNEW - HOLD )                         

!                A = DELTA * (ETA/ETFG)**1.25                           
!                B = DELTA * (ETA/ETFG)**.25                            

!       IF THERE IS MAINLY WINDSEA, THE STRETCHING CONSTANT B IS        
!       COMPUTED TO PRODUCE IN THE ANALYSED  SPECTRUM THE ANALYSED      
!       MEAN FREQUENCY DERIVED BY THE MODEL GROWTH CURVE                
!       ( AS CARRIED OUT IN THE SUBROUTINE FUSTAR ). IN THIS CASE       

!                  B = FMWFG/FMWA                                       
!                  A = (ETA/ETFG)*B                                     

!-----------------------------------------------------------------------

      USE YOWFRED  , ONLY : FR
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE

! ----------------------------------------------------------------------

      DIMENSION F(NINF-1:NSUP,NANG,NFRE)
      REAL,DIMENSION(NINF:NSUP) :: ETA, FMWFG, FMWA
!cc note ETFG corresponds to EMEAN which is module yowmean in UPDATE
!cc therefore it has to have a different definition as the other arrays
      REAL,DIMENSION(IJS:IJL) :: ETFG

! ----------------------------------------------------------------------

      REAL, PARAMETER :: XL11 = 0.0953101
      REAL, ALLOCATABLE :: FTEMP(:,:)

!          *XL11*   REAL     ALOG(1.1)                                  
!          *FTEMP*  REAL     TEMPORAR STORAGE OF THE SPECTRUM.          

! ----------------------------------------------------------------------

      ALLOCATE(FTEMP(NANG,NFRE))                                        

      FR1=FR(1)                                                         

!*    1. LOOP OVER GRID POINTS.                                         
!        ---------------------                                          

      DO IJ=IJS,IJL                                                

!*    1.1 SKIP LAND POINTS AND POINTS WHERE THERE ARE NO RELIABLE DATA. 
!         ------------------------------------------------------------- 

        IF (ETFG(IJ).LE.0.001.OR.ETA(IJ).LE.0.001) GOTO 1001           
        FNEW = FMWA(IJ)                                                
        HNEW = 4.*SQRT(ETA(IJ))                                        
        HOLD = 4.*SQRT(ETFG(IJ))                                       
        FOLD = FMWFG(IJ)                                               

!*    1.2 COMPUTE SCALING AND STRETCHING FACTORS.                       
!         ---------------------------------------                       

        IF (FOLD.LE.0.) THEN                                           

!*    1.2.1 THE SPECTRUM IS MAINLY SWELL (THE WINDSEA MEAN FREQUENCY    
!*          OF THE FIRST GUESS SPECTRUM IS NEGATIVE).                   
!           --------------------------------------------------------    

          XDELTA = 1.-0.006*(HNEW-HOLD)                               
          XR = XDELTA*(HNEW/HOLD)**2.5                                
          XB = XDELTA* SQRT(HNEW/HOLD)                                
        ELSE                                                           

!*    1.2.2 THE SPECTRUM IS MAINLY WINDSEA.                             
!           -------------------------------                             

          XR = (HNEW/HOLD)**2*FOLD/FNEW                               
          XB = FOLD/FNEW                                              
        ENDIF                                                         

!*    1.3 LOOP OVER NEW FREQUENCIES.                                    
!         --------------------------                                    

        DO M=1,NFRE                                               
          FU = FR(M)*XB                                               
          M1 = INT(LOG(FU/FR1)/XL11)+1                               
          IF (M1.GE.NFRE.OR.M1.LT.1) THEN                             

!*    1.3.1 LOOP OVER DIRECTIONS FOR FREQUENCIES GETTING ENERGY FORM    
!*          FREQUENCIES OUT OF RANGE.                                   
!           ---------------------------------------------------------   

            DO K=1,NANG                                         
              FTEMP(K,M) = 0.                                       
            ENDDO
          ELSE                                                        

!*    1.3.2 LOOP OVER DIRECTIONS FOR FREQUENCIES GETTING ENERGY.        
!           ----------------------------------------------------        
            M2 = M1+1                                                
            DFRE = FR(M1)*.1                                         
            DO K=1,NANG                                         
              F1 = F(IJ,K,M1)                                       
              F2 = F(IJ,K,M2)                                       
              DE = (F2-F1)/DFRE * (FU-FR(M1))                       
              FTEMP(K,M) = XR*MAX(0.,F1+DE)                         
            ENDDO
          ENDIF                                                      

!*    BRANCH BACK TO 1.3 FOR NEXT FREQUENCY.                            

        ENDDO

!*    1.4 THE UPDATED SPECTRUM IS STORED.                               
!         -------------------------------                               

        DO K=1,NANG                                               
          DO M=1,NFRE                                            
            F(IJ,K,M) = FTEMP(K,M)                                   
          ENDDO
        ENDDO

!*    BRANCH BACK TO 1. FOR NEXT GRIDPOINT.                             

 1001   continue
      ENDDO
                                                                        
      DEALLOCATE(FTEMP)                                        

      RETURN                                                            
      END SUBROUTINE UPWSPEC
