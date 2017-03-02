      SUBROUTINE UPWSPEC (F, ETFG, AKFG, ETA, FMWFG, FMWA, IJS, IJL)

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

!        *CALL* *UPWSPEC (F, ETFG, AKFG, ETA, FMWFG, FMWA, IJS, IJL)*          

!          *F*      REAL     INPUT IS THE OLD SWELL SPECTRUM.           
!                            OUTPUT IS  THE ANALISED SPECTRUM.          
!          *ETFG*   REAL     FIRST GUESS TOTAL ENERGY.                  
!          *AKFG*   REAL     FIRST GUESS MEAN WAVE NUMBER.                  
!          *ETA*    REAL     ANALYSED TOTAL ENERGY.                     
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
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : ZPI      ,EPSMIN   ,ASWKM    ,BSWKM 
      USE YOWSHAL  , ONLY : DEPTH
      USE YOWSTAT  , ONLY : ISHALLO
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: IJ, IJS, IJL, M, K, ICOUNT 
      INTEGER :: NFREINF
      INTEGER :: IG 
      INTEGER :: M1(NFRE)

      REAL :: XL11
      REAL :: HNEW, FNEW, HOLD, FOLD
      REAL :: XKMEAN_MIN, XDELTA, XKMEAN_NEW
      REAL :: OMEGA, XK, ESPFG_1D
      REAL :: DELET, DELLOW, DELHIGH, ESPFG 
      REAL :: WEIGHT, F1, F2, DE 
      REAL :: ESPAN, ESPAN_1D    
      REAL :: AKI 
      REAL :: ZHOOK_HANDLE
      REAL :: F(IJS:IJL,NANG,NFRE)
      REAL, DIMENSION(IJS:IJL) :: ETA, FMWFG, FMWA
      REAL, DIMENSION(IJS:IJL) :: ETFG, AKFG, XR, XB 
      REAL :: FTEMP(NANG,NFRE), DFRE(NFRE), FU(NFRE)

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('UPWSPEC',0,ZHOOK_HANDLE)
#endif

      IG=1

      XL11 = LOG(FRATIO)

      IF(FR(1).LE.0.04) THEN
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

!           FOR RELATIVELY SHALLOW CASES,
!           THE ASSUMPTION OF CONSTANT STEEPNESS IS NOT APPROPRIATE.
!           RELAX THE SCHEME TO ONLY ADAPT ENERGY LEVEL  
!           (I.E. THE WAVES WILL GET STEEPER)
!           SUCH THAT THE MINIMUM MEAN WAVE NUMBER IS INFORCED.

            XKMEAN_MIN=ASWKM/DEPTH(IJ,IG)**BSWKM
            IF(ISHALLO.NE.1 .AND. AKFG(IJ).LT.XKMEAN_MIN) THEN
              XR(IJ) = (HNEW/HOLD)**2
              XB(IJ) = 1. 
            ELSE
              XDELTA = 1.-0.006*(HNEW-HOLD)
              XB(IJ) = XDELTA* SQRT(HNEW/HOLD)
!!!testing that new xkmean does not violate the condition on minimum xkmean???
              IF(ISHALLO.NE.1) THEN
                ICOUNT=0
1111            ICOUNT=ICOUNT+1
                IF(XB(IJ).GT.1.AND.ICOUNT.LE.10) THEN
                  XKMEAN_NEW=EPSMIN
                  DO M=1,NFRE
                     OMEGA=ZPI*FR(M)/XB(IJ)
                     XK=AKI(OMEGA,DEPTH(IJ,IG))
                     ESPFG_1D=0.
                     DO K=1,NANG
                       ESPFG_1D=ESPFG_1D+F(IJ,K,M)
                     ENDDO
                     XKMEAN_NEW=XKMEAN_NEW+DFIM(M)*ESPFG_1D/SQRT(XK)
                  ENDDO
                  XKMEAN_NEW=(HOLD**2/(16*XKMEAN_NEW))**2
                  IF(XKMEAN_NEW.LT.XKMEAN_MIN) THEN
                    XB(IJ)=0.5*(1.+XB(IJ))
!                   iterate
                    GOTO 1111
                  ENDIF
                ENDIF
              ENDIF
!!!!
              XR(IJ) = XB(IJ)*(HNEW/HOLD)**2
            ENDIF
          ELSE 
!*          THE SPECTRUM IS MAINLY WINDSEA.
            XR(IJ) = (HNEW/HOLD)**2*FOLD/FNEW
            XB(IJ) = FOLD/FNEW
          ENDIF

          DELET=ABS(ETA(IJ)-ETFG(IJ))
          DELLOW=0.8*DELET
          DELHIGH=1.2*DELET 

          ESPFG=0.

          DO M=NFREINF,NFRE
            ESPFG_1D=0.
            DO K=1,NANG
              ESPFG_1D=ESPFG_1D+F(IJ,K,M)
            ENDDO
            ESPFG=ESPFG+DFIM(M)*ESPFG_1D
          ENDDO

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

          IF(ABS(ESPAN-ESPFG).LT.DELHIGH .AND.
     &       ABS(ESPAN-ESPFG).GT.DELLOW ) THEN
            DO K=1,NANG
              DO M=NFREINF,NFRE
                F(IJ,K,M) = FTEMP(K,M)
              ENDDO
            ENDDO
          ENDIF
   
        ENDIF
      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('UPWSPEC',1,ZHOOK_HANDLE)
#endif
      RETURN
      END SUBROUTINE UPWSPEC
