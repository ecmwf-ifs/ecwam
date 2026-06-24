      SUBROUTINE UPWSPEC (KIJS, KIJL, FL1, DEPTH, ETFG, AKFG, ETA, FMWFG, FMWA)

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

!        *CALL* *UPWSPEC (KIJS, KIJL, FL1, DEPTH, ETFG, AKFG, ETA, FMWFG, FMWA, KIJS, KIJL)*          

!          *KIJS*   INTEGER   FIRST INDEX IN BLOCK.
!          *KIJL*   INTEGER   LAST  INDEX IN BLOCK.
!          *FL1*    REAL     INPUT IS THE OLD SWELL SPECTRUM.           
!                            OUTPUT IS  THE ANALYSED SPECTRUM.          
!          *DEPTH*  REAL     WATER DEPTH.
!          *ETFG*   REAL     FIRST GUESS TOTAL ENERGY.                  
!          *AKFG*   REAL     FIRST GUESS MEAN WAVE NUMBER.                  
!          *ETA*    REAL     ANALYSED TOTAL ENERGY.                     
!          *FMWFG   REAL     MEAN FREQUENCY OF THE WINDSEA FROM THE     
!                            FIRST GUESS SPECTRUM.                      
!          *FMWA*   REAL     ANALYSED MEAN FREQUENCY OF THE WINDSEA.    

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

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWALTAS , ONLY : ASWKM    ,BSWKM 
      USE YOWFRED  , ONLY : FR       ,ZPIFR   ,DFIM     ,FRATIO
      USE YOWMAP   , ONLY : XDELLA
      USE YOWSHAL  , ONLY : DEPTHA
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "aki.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: ETFG, AKFG
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: ETA, FMWFG, FMWA

      INTEGER(KIND=JWIM) :: IJ, M, K, ICOUNT 
      INTEGER(KIND=JWIM) :: NFREINF
      INTEGER(KIND=JWIM), DIMENSION(NFRE) :: M1

      REAL(KIND=JWRB) :: ZMINDEPTH
      REAL(KIND=JWRB) :: XL11
      REAL(KIND=JWRB) :: HNEW, FNEW, HOLD, FOLD
      REAL(KIND=JWRB) :: XKMEAN_MIN, XDELTA, XKMEAN_NEW
      REAL(KIND=JWRB) :: OMEGA, XK, ESPFG_1D
      REAL(KIND=JWRB) :: DELET, DELLOW, DELHIGH, ESPFG 
      REAL(KIND=JWRB) :: WEIGHT, F1, F2, DE 
      REAL(KIND=JWRB) :: ESPAN, ESPAN_1D    
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      REAL(KIND=JWRB), DIMENSION(NFRE) :: DFRE, FU
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE) :: FTEMP
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: XR, XB 

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('UPWSPEC',0,ZHOOK_HANDLE)

      XL11 = LOG(FRATIO)

      IF (FR(1) <= 0.04_JWRB) THEN
         NFREINF=4
      ELSE
         NFREINF=2
      ENDIF

      DO M=1,NFRE                                               
        DFRE(M) = FR(M)*(FRATIO-1.0_JWRB)
      ENDDO

      IF(XDELLA < 0.125_JWRB) THEN
!       This is a hack to avoid cases very near coast in very shallow water at high resolution
!       which retuned very weird mean wave period increments
!       This issue will need to get solved more robustly.
        ZMINDEPTH = 1.5_JWRB * DEPTHA
      ELSE
        ZMINDEPTH = 0.0_JWRB
      ENDIF

      DO IJ=KIJS,KIJL                                                
!*      SKIP LAND POINTS AND POINTS WHERE THERE ARE NO RELIABLE DATA.
        IF (DEPTH(IJ) > ZMINDEPTH .AND.  ETFG(IJ) > 0.001_JWRB .AND. ETA(IJ) > 0.001_JWRB) THEN
          HNEW = 4.0_JWRB*SQRT(ETA(IJ))
          FNEW = FMWA(IJ)
          HOLD = 4.0_JWRB*SQRT(ETFG(IJ))
          FOLD = FMWFG(IJ)

!*        COMPUTE SCALING AND STRETCHING FACTORS. 

          IF (FOLD <= 0.0_JWRB) THEN
!*          THE SPECTRUM IS MAINLY SWELL (THE WINDSEA MEAN FREQUENCY
!*          OF THE FIRST GUESS SPECTRUM IS NEGATIVE).

!           FOR RELATIVELY SHALLOW CASES,
!           THE ASSUMPTION OF CONSTANT STEEPNESS IS NOT APPROPRIATE.
!           RELAX THE SCHEME TO ONLY ADAPT ENERGY LEVEL  
!           (I.E. THE WAVES WILL GET STEEPER)
!           SUCH THAT THE MINIMUM MEAN WAVE NUMBER IS INFORCED.

            XKMEAN_MIN=ASWKM/DEPTH(IJ)**BSWKM
            IF (AKFG(IJ) < XKMEAN_MIN) THEN
              XR(IJ) = (HNEW/HOLD)**2
              XB(IJ) = 1.0_JWRB
            ELSE
              XDELTA = 1.0_JWRB-0.006_JWRB*(HNEW-HOLD)
              XB(IJ) = XDELTA* SQRT(HNEW/HOLD)
!!!testing that new xkmean does not violate the condition on minimum xkmean???
                ICOUNT=0
1111            ICOUNT=ICOUNT+1
                IF (XB(IJ) > 1.0_JWRB .AND. ICOUNT <= 10) THEN
                  XKMEAN_NEW=EPSMIN
                  DO M=1,NFRE
                     OMEGA=ZPIFR(M)/XB(IJ)
                     XK=AKI(OMEGA,DEPTH(IJ))
                     ESPFG_1D=0.0_JWRB
                     DO K=1,NANG
                       ESPFG_1D=ESPFG_1D+FL1(IJ,K,M)
                     ENDDO
                     XKMEAN_NEW=XKMEAN_NEW+DFIM(M)*ESPFG_1D/SQRT(XK)
                  ENDDO
                  XKMEAN_NEW=(HOLD**2/(16._JWRB*XKMEAN_NEW))**2
                  IF (XKMEAN_NEW < XKMEAN_MIN) THEN
                    XB(IJ)=0.5_JWRB*(1.0_JWRB+XB(IJ))
!                   iterate
                    GOTO 1111
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
          DELLOW=0.8_JWRB*DELET
          DELHIGH=1.2_JWRB*DELET 

          ESPFG=0.0_JWRB

          DO M=NFREINF,NFRE
            ESPFG_1D=0.0_JWRB
            DO K=1,NANG
              ESPFG_1D=ESPFG_1D+FL1(IJ,K,M)
            ENDDO
            ESPFG=ESPFG+DFIM(M)*ESPFG_1D
          ENDDO

!*        LOOP OVER NEW FREQUENCIES.

          DO M=NFREINF,NFRE
            FU(M) = FR(M)*XB(IJ)
            M1(M) = INT(LOG(FU(M)/FR(1))/XL11)+1
          ENDDO

          DO M=NFREINF,NFRE

            IF (M1(M) >= NFRE) THEN
!             A F**-5 LAW IS ASSUMED FOR FREQUENCIES ABOVE FR(NFRE)
              DO K=1,NANG
                FTEMP(K,M) = XR(IJ)*FL1(IJ,K,NFRE)*(FR(NFRE)/FU(M))**5 
              ENDDO

            ELSEIF (M1(M) < 1) THEN
              DO K=1,NANG
                FTEMP(K,M) = 0.0_JWRB
              ENDDO

            ELSE

              WEIGHT=(FU(M)-FR(M1(M)))/DFRE(M1(M))
              DO K=1,NANG
                F1 = FL1(IJ,K,M1(M))
                F2 = FL1(IJ,K,M1(M)+1)
                DE = (F2-F1)*WEIGHT
                FTEMP(K,M) = MAX(XR(IJ)*(F1+DE),0.0_JWRB)
              ENDDO
            ENDIF

          ENDDO

!*        THE UPDATED SPECTRUM IS STORED. (provided it agrees with
!                                          the energy increment)
          ESPAN=0.0_JWRB
          DO M=NFREINF,NFRE
            ESPAN_1D=0.0_JWRB
            DO K=1,NANG
              ESPAN_1D=ESPAN_1D+FTEMP(K,M) 
            ENDDO
            ESPAN=ESPAN+DFIM(M)*ESPAN_1D
          ENDDO

          IF (ABS(ESPAN-ESPFG) < DELHIGH .AND. ABS(ESPAN-ESPFG) > DELLOW ) THEN
            DO K=1,NANG
              DO M=NFREINF,NFRE
                FL1(IJ,K,M) = FTEMP(K,M)
              ENDDO
            ENDDO
          ENDIF
   
        ENDIF
      ENDDO

      IF (LHOOK) CALL DR_HOOK('UPWSPEC',1,ZHOOK_HANDLE)

      END SUBROUTINE UPWSPEC
