      SUBROUTINE UPWSPEC (F, ETFG, ETA, FMWFG, FMWA, IJS, IJL)

!-----------------------------------------------------------------------

!**** *UPWSPEC* - MODIFY THE SPECTRUM BY STRETCHING AND SCALING.         

!     P.LIONELLO      ECMWF         APRIL 1990                          
!     J.BIDLOT        ECMWF         JUNE 1997 UPSPEC >> UPWSPEC 
!     J.BIDLOT        ECMWF         APRIL 2000 ADD f-5 TAIL WHEN
!                                   UPDATED SPECTRUM SWIFTS TO LOWER
!                                   FREQUENCIES
!     J.BIDLOT        ECMWF         CHECK IF UPDATED SPECTRUM IS 
!                                   IN AGREEMENT WITH ENERGY INCREMENT.
!                                   DO NOT ALLOW ANY UPDATES FOR THE
!                                   VERY LOW FREQUENCIES.

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

      USE YOWFRED  , ONLY : FR       ,DFIM     ,FRATIO
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE

! ----------------------------------------------------------------------

      INTEGER :: M1(NFRE)
      REAL :: F(NINF-1:NSUP,NANG,NFRE)
      REAL,DIMENSION(NINF:NSUP) :: ETA, FMWFG, FMWA
!cc note ETFG corresponds to EMEAN which is module yowmean in UPDATE
!cc therefore it has to have a different definition as the other arrays
      REAL,DIMENSION(IJS:IJL) :: ETFG, XR, XB 

      REAL :: XL11
      REAL :: FTEMP(NANG,NFRE), DFRE(NFRE), FU(NFRE)
      REAL :: ESPFG(IJS:IJL), ESPFG_1D(IJS:IJL) 
      REAL :: DELLOW(IJS:IJL), DELHIGH(IJS:IJL) 

! ----------------------------------------------------------------------

      XL11 = LOG(FRATIO)

      IF(NFRE.EQ.30 .OR. NFRE.EQ.36) THEN
         NFREINF=4
      ELSE
         NFREINF=2
      ENDIF

      DO M=1,NFRE                                               
        DFRE(M) = FR(M)*(FRATIO-1.)
      ENDDO

      DO IJ=IJS,IJL                                                
!*      SKIP LAND POINTS AND POINTS WHERE THERE ARE NO RELIABLE DATA.
        IF (ETFG(IJ).GT.0.001.AND.ETA(IJ).GT.0.001) THEN
          HNEW = 4.*SQRT(ETA(IJ))
          FNEW = FMWA(IJ)
          HOLD = 4.*SQRT(ETFG(IJ))
          FOLD = FMWFG(IJ)

!*        COMPUTE SCALING AND STRETCHING FACTORS. 

          IF (FOLD.LE.0.) THEN
!*          THE SPECTRUM IS MAINLY SWELL (THE WINDSEA MEAN FREQUENCY
!*          OF THE FIRST GUESS SPECTRUM IS NEGATIVE).
            XDELTA = 1.-0.006*(HNEW-HOLD)
            XR(IJ) = XDELTA*(HNEW/HOLD)**2.5
            XB(IJ) = XDELTA* SQRT(HNEW/HOLD)
          ELSE 
!*          THE SPECTRUM IS MAINLY WINDSEA.
            XR(IJ) = (HNEW/HOLD)**2*FOLD/FNEW
            XB(IJ) = FOLD/FNEW
          ENDIF

          DELET=ABS(ETA(IJ)-ETFG(IJ))
          DELLOW(IJ)=0.8*DELET
          DELHIGH(IJ)=1.2*DELET 
        ENDIF
      ENDDO

      DO IJ=IJS,IJL                                                
        ESPFG(IJ)=0.
      ENDDO

      DO M=NFREINF,NFRE
        DO IJ=IJS,IJL                                                
          ESPFG_1D(IJ)=0.
        ENDDO
        DO K=1,NANG
          DO IJ=IJS,IJL                                                
            ESPFG_1D(IJ)=ESPFG_1D(IJ)+F(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=IJS,IJL                                                
          ESPFG(IJ)=ESPFG(IJ)+DFIM(M)*ESPFG_1D(IJ)
        ENDDO
      ENDDO

      DO IJ=IJS,IJL                                                

!*      SKIP LAND POINTS AND POINTS WHERE THERE ARE NO RELIABLE DATA.

        IF (ETFG(IJ).GT.0.001.AND.ETA(IJ).GT.0.001) THEN

!*        LOOP OVER NEW FREQUENCIES.

          DO M=NFREINF,NFRE
            FU(M) = FR(M)*XB(IJ)
            M1(M) = INT(LOG(FU(M)/FR(1))/XL11)+1
          ENDDO

          DO M=NFREINF,NFRE

            IF (M1(M).GE.NFRE) THEN
!             A F**-5 LAW IS ASSUMED FOR FREQUENCIES ABOVE FR(NFRE)
              DO K=1,NANG
                FTEMP(K,M) = XR(IJ)*F(IJ,K,NFRE)*(FR(NFRE)/FU(M))**5 
              ENDDO

            ELSEIF (M1(M).LT.1) THEN
              DO K=1,NANG
                FTEMP(K,M) = 0.
              ENDDO

            ELSE
              WEIGHT=(FU(M)-FR(M1(M)))/DFRE(M1(M))
              DO K=1,NANG
                F1 = F(IJ,K,M1(M))
                F2 = F(IJ,K,M1(M)+1)
                DE = (F2-F1)*WEIGHT
                FTEMP(K,M) = MAX(XR(IJ)*(F1+DE),0.)
              ENDDO
            ENDIF

          ENDDO

!*        THE UPDATED SPECTRUM IS STORED. (provided it agrees with
!                                          the energy increment)
          ESPAN=0.
          DO M=NFREINF,NFRE
            ESPAN_1D=0.
            DO K=1,NANG
              ESPAN_1D=ESPAN_1D+FTEMP(K,M) 
            ENDDO
            ESPAN=ESPAN+DFIM(M)*ESPAN_1D
          ENDDO

          IF(ABS(ESPAN-ESPFG(IJ)).LT.DELHIGH(IJ) .AND.
     &       ABS(ESPAN-ESPFG(IJ)).GT.DELLOW(IJ) ) THEN
            DO K=1,NANG
              DO M=NFREINF,NFRE
                F(IJ,K,M) = FTEMP(K,M)
              ENDDO
            ENDDO
          ENDIF
   
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE UPWSPEC
