      SUBROUTINE UPDATE (F, ETOI, USOI, USMO, THMO, IJS, IJL)           

! ----------------------------------------------------------------------

!**** *UPDATE* - ANALYSE THE WAVE SPECTRUM, PRODUCING A CONSISTENT      
!****            UPDATE OF WAVE AND WIND FIELD.                         

!     PIERO LIONELLO      ECMWF     JUNE 1990                           

!     PURPOSE.                                                          
!     --------                                                          

!       TO ANALYSE THE WAVE SPECTRUM, PRODUCING A CONSISTENT            
!       UPDATE OF WAVE AND WIND FIELD.                                  

!**   INTERFACE.                                                        
!     ----------                                                        

!       *CALL* *UPDATE (F, ETOI, USOI, USMO, THMO, IJS, IJL)*           

!         *F*       REAL      SPECTRUM ( FIRST GUESS IN INPUT ;         
!                             ANALYSIS IN OUTPUT ).                     
!         *ETOI*    REAL      WAVE ENERGY (FROM O.I.).                  
!         *USOI*    REAL      USTAR ( FROM O.I.).                       
!         *USMO*    REAL      USTAR ( FIRST GUESS IN INPUT ;            
!                             ANALYSIS IN OUTPUT ).                     
!         *THMO*    REAL      WIND DIRECTION.                           
!         *IJS*     INTEGER   FIRST INDEX IN BLOCK.                     
!         *IJL*     INTEGER   LAST  INDEX IN BLOCK.                     

!     METHOD.                                                           
!     -------                                                           

!     EXTERNALS.
!     ----------
!     *SEMEAN*
!     *FWSEA*
!     *FDUR*
!     *FUSTAR*
!     *UPWSPEC*

! ----------------------------------------------------------------------

      USE YOWMEAN  , ONLY : EMEAN
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------

      LOGICAL :: LLEPSMIN
      DIMENSION F(NINF-1:NSUP,NANG,NFRE)
      REAL,DIMENSION(NINF:NSUP) :: ETOI, USOI, USMO, THMO         
      REAL,DIMENSION(NINF:NSUP) :: EWFG, EWOI, EWA,FMWFG, FMWA,
     &                             THEW, TDUR, USA

!       *EWFG*     : WINDSEA ENERGY (FIRST GUESS).                      
!       *EWOI*     : WINDSEA ENERGY (MEASUREMENT).                      
!       *EWA*      : WINDSEA ENERGY (ANALYSIS).                         
!       *ETOI*     : TOTAL ENERGY (FROM O.I ,IDENTICAL WITH ANALYSIS IN 
!                    THIS METHOD ).                                     
!       *FMWFG*    : MEAN FREQUENCY OF THE WINDSEA (FIRST GUESS).       
!       *FMWA*     : MEAN FREQUENCY OF THE WINDSEA (ANALYSIS).          
!       *THEW*     : DIRECTION OF THE WINDSEA PEAK (FIRST GUESS).       
!       *TDUR*     : WINDSEA DURATION.                                  
!       *USA*      : ESTIMATE OF USTAR FROM THE FIRST GUESS ENERGY.     
!                    AND DURATION (I.E. ANALYSED USTAR).                

! ----------------------------------------------------------------------

!*    1. FIND THE WINDSEA.                                              
!        -----------------                                              

      LLEPSMIN=.TRUE.
      CALL SEMEAN (F, IJS, IJL, EMEAN(IJS), LLEPSMIN)

      CALL FWSEA (F, EWFG, THEW, FMWFG, ETOI, USMO, THMO, IJS, IJL)     

!*    2. CORRECT OVERESTIMATE OF FIRST GUESS WINDSEA ENERGY.            
!        ---------------------------------------------------            

      DO IJ = IJS,IJL                                              
        EWFG(IJ) = MIN(EWFG(IJ),EMEAN(IJ))                             
      ENDDO
                                                                        
      IF (ITEST.GE.3) WRITE (IU06,*) '      SUB. UPDATE: ',             
     &   ' WINDSEA ENERGY FIXED IN SUB. FWSEA'                          

! ----------------------------------------------------------------------

!*    3. FIND THE DURATION.                                             
!        ------------------                                             

      CALL FDUR (EWFG, TDUR, USMO, IJS, IJL)
                                                                        
      IF (ITEST.GE.3) WRITE (IU06,*) '      SUB. UPDATE: ',             
     &   ' WINDSEA DURATION FIXED IN SUB. FDUR'                         

! ----------------------------------------------------------------------

!*    4. THE RATIO WINDSEA/SWELL IS ASSUMED CORRECT.                    
!        -------------------------------------------                    

      DO IJ=IJS,IJL                                                
        IF (EMEAN(IJ).GT.0.) EWOI(IJ) = EWFG(IJ)/EMEAN(IJ)*ETOI(IJ)    
      ENDDO

! ----------------------------------------------------------------------

!*    5. FIND THE NEW USTAR AND THE NEW WINDSEA MEAN FREQUENCY.         
!        ------------------------------------------------------         

      CALL FUSTAR (USMO, USA, USOI, EWFG, EWOI, EWA, TDUR, FMWA,
     &             IJS, IJL) 

      DO IJ=IJS,IJL                                                
        IF (USA(IJ).GT.0.) USMO(IJ) = USA(IJ)                          
      ENDDO
                                                                        
      IF (ITEST.GE.3) WRITE (IU06,*) '      SUB. UPDATE: ',             
     &   ' USTAR COMPUTED IN SUB. FUSTAR'                               

! ----------------------------------------------------------------------

!*    6. THE ASSIMILATION IS ADJUSTED ON THE DOMINANT PART              
!*       OF THE SPECTRUM.                                               
!        -------------------------------------------------              

!        THE FOLLOWING RETURN CODE ARE POSSIBLE :                       
!          THERE IS NO WINDSEA ->                                       
!                 EWFG  =  -999.                                        
!                 FMWFG = -9999.                                        
!          THE RATIO WINDSEA/SWELL IS NOT CONSISTENT WITH USOI ->       
!                 EWA   =  -999.                                        
!                 FMWA  = -9999.                                        

      DO IJ=IJS,IJL                                                
        IF (EWA(IJ).GT.0. .AND. EWA(IJ).LT.0.5*ETOI(IJ)) THEN          
          EWFG(IJ)  =  -999.                                          
          FMWFG(IJ) = -9999.                                          
        ENDIF                                                         
      ENDDO

!     IF THERE WERE WINDSEA BUT THE RATIO WINDSEA SWELL     
!     WAS WRONG, DO NOT ASSIMILATE.                                     

      DO IJ=IJS,IJL                                                
        IF (EWFG(IJ).GT.0. .AND. EWA(IJ).LT.0.) THEN        
          ETOI(IJ)=-999.                                              
        ENDIF                                                          
      ENDDO

! ----------------------------------------------------------------------

!*    7. UPDATE THE SPECTRA.                                            
!        -------------------                                            

      CALL UPWSPEC (F, EMEAN(IJS), ETOI, FMWFG, FMWA, IJS, IJL)
                                                                        
      IF (ITEST.GE.3) WRITE (IU06,*) '      SUB. UPDATE: ',             
     &   ' NEW SPECTRA COMPUTED IN SUB. UPWSPEC'
                                                                        
      RETURN                                                            
      END SUBROUTINE UPDATE
