      SUBROUTINE UPDATE (F, ETOI, USMO, THMO, IJS, IJL)           

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

!       *CALL* *UPDATE (F, ETOI, USMO, THMO, IJS, IJL)*           

!         *F*       REAL      SPECTRUM ( FIRST GUESS IN INPUT ;         
!                             ANALYSIS IN OUTPUT ).                     
!         *ETOI*    REAL      WAVE ENERGY (FROM O.I.).                  
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

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------

      LOGICAL :: LLEPSMIN
      DIMENSION F(IJS:IJL,NANG,NFRE)
      REAL,DIMENSION(IJS:IJL) :: EMN
      REAL,DIMENSION(IJS:IJL) :: ETOI, USMO, THMO         
      REAL,DIMENSION(IJS:IJL) :: EWFG, EWOI, EWA,FMWFG, FMWA,
     &                           TDUR, USA

!       *EWFG*     : WINDSEA ENERGY (FIRST GUESS).                      
!       *EWOI*     : WINDSEA ENERGY (MEASUREMENT).                      
!       *EWA*      : WINDSEA ENERGY (ANALYSIS).                         
!       *ETOI*     : TOTAL ENERGY (FROM O.I ,IDENTICAL WITH ANALYSIS IN 
!                    THIS METHOD ).                                     
!       *FMWFG*    : MEAN FREQUENCY OF THE WINDSEA (FIRST GUESS).       
!       *FMWA*     : MEAN FREQUENCY OF THE WINDSEA (ANALYSIS).          
!       *TDUR*     : WINDSEA DURATION.                                  
!       *USA*      : ESTIMATE OF USTAR FROM THE FIRST GUESS ENERGY.     
!                    AND DURATION (I.E. ANALYSED USTAR).                

! ----------------------------------------------------------------------

!*    1. FIND THE WINDSEA.                                              
!        -----------------                                              

      LLEPSMIN=.TRUE.
      CALL SEMEAN (F, IJS, IJL, EMN, LLEPSMIN)

      CALL FWSEA (F, EMN, EWFG, FMWFG, ETOI, USMO, THMO, IJS, IJL)     

!*    2. CORRECT OVERESTIMATE OF FIRST GUESS WINDSEA ENERGY.            
!        ---------------------------------------------------            

      DO IJ = IJS,IJL                                              
        EWFG(IJ) = MIN(EWFG(IJ),EMN(IJ))                             
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
        IF (EMN(IJ).GT.0.) EWOI(IJ) = EWFG(IJ)/EMN(IJ)*ETOI(IJ)    
      ENDDO

! ----------------------------------------------------------------------

!*    5. FIND THE NEW USTAR AND THE NEW WINDSEA MEAN FREQUENCY.         
!        ------------------------------------------------------         

      CALL FUSTAR (USMO, USA, EWFG, EWOI, EWA, TDUR, FMWA, IJS, IJL)

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

      DO IJ=IJS,IJL                                                
        IF (EWA(IJ).GT.0. .AND. EWA(IJ).LT.0.6*ETOI(IJ)) THEN          
          EWFG(IJ)  =  -999.                                          
          FMWFG(IJ) = -9999.                                          
        ENDIF                                                         
      ENDDO


!     IF THERE WAS A PROBLEM IN DETERMINING USTAR THEN
!     NO WINDSEA UPDATE SHOULD BE MADE
      DO IJ=IJS,IJL                                                
        IF (USA(IJ).LT.0.) THEN        
          EWFG(IJ)  =  -999.                                          
          FMWFG(IJ) = -9999.                                          
        ENDIF                                                          
      ENDDO

!     UPDATE USTAR ONLY WHEN THE MEAN WINDSEA FREQUENCY AND
!     ANALYSED WAVE ENERGY ARE POSITIVE
!     I.E. WHEN THE UPDATE ON SPECTRA WILL ASSUME IT IS MAINLY WINDSEA

      DO IJ=IJS,IJL
        IF (USA(IJ).GT.0. .AND. FMWFG(IJ) .GT. 0. .AND.
     &      ETOI(IJ) .GT. 0.) USMO(IJ) = USA(IJ)
      ENDDO
                                                                        
! ----------------------------------------------------------------------

!*    7. UPDATE THE SPECTRA.                                            
!        -------------------                                            

      CALL UPWSPEC (F, EMN, ETOI, FMWFG, FMWA, IJS, IJL)
                                                                        
      IF (ITEST.GE.3) WRITE (IU06,*) '      SUB. UPDATE: ',             
     &   ' NEW SPECTRA COMPUTED IN SUB. UPWSPEC'
                                                                        
      RETURN                                                            
      END SUBROUTINE UPDATE
