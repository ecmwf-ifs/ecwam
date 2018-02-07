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

      USE YOWALTAS , ONLY : AFCRV   ,BFCRV     ,
     &            ESH      ,ASH     ,BSH
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH   ,DFFR      ,
     &            WETAIL   ,WP1TAIL
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,EPSMIN
      USE YOWSTAT  , ONLY : ISHALLO

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL, INTENT(IN) :: USTT, DPT
      REAL, INTENT(IN) :: FSEA(NANG,NFRE)
      REAL, INTENT(OUT) :: EW, FM

      INTEGER :: M, K
      REAL :: EN, YNU, EMAX, X
      REAL :: DELT25, COEF1, XTEMP, TEMP, SPFB
      REAL :: ESTAR, DSTAR, SPINTDI, FREQDI 

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
!*    1.2 COMPUTATION OF THE MEAN FREQUENCY.                            
!         -----------------------------------                           

      EW = EPSMIN
      FM = EPSMIN
      DO M=1,NFRE
        TEMP = 0.
        DO K=1,NANG
          TEMP = TEMP+FSEA(K,M)
        ENDDO
        EW = EW+TEMP*DFIM(M)
        FM = FM+TEMP*DFFR(M)
      ENDDO
      DELT25 = WETAIL*FR(NFRE)*DELTH
      EW = EW+DELT25*TEMP                                               
      COEF1 = WP1TAIL*DELTH*FR(NFRE)**2
      FM = FM+COEF1*TEMP
      FM = FM/EW

! ----------------------------------------------------------------------

!*    2. ESTIMATES FROM MODEL RELATIONS.                                
!        -------------------------------                                

!*    2.1 ENERGY IS DERIVED FROM UNDERESTIMATED MEAN FREQUENCY.         
!         -----------------------------------------------------         

      XTEMP = FM*USTT/G

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

!*    3. FINAL ESTIMATE
!        ---------------

!*    3.1 AVERAGING THE ENERGY ESTIMATES
!         -------------------------------

      EW = (EW+SPINTDI)*0.5

!*    3.2 AVERAGING THE MEAN FREQUENCY ESTIMATES.                       
!         ---------------------------------------                       

      FM = (FM+FREQDI)*0.5
 
      END SUBROUTINE WSMFEN
