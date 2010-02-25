      SUBROUTINE FUSTAR (USMO, USA, EWFG, EWOI, EWA, T, FMWA, IJS, IJL)

! ----------------------------------------------------------------------

!**** *FUSTAR* - ESTIMATE THE ANALYSED FRICTION VELOCITY FROM THE       
!                DURATION AND THE "MEASURED" WINDSEA ENERGY.            
!                ESTIMATE THE ANALYSED WINDSEA MEAN FREQUENCY           
!                FROM THE "MEASURED" WINDSEA ENERGY                     

!     P.LIONELLO     ECMWF       FEBRUARY 1989                          
!         (FROM MODIFICATION OF A CODE BY P.JANSSEN                     
!           AND P.LIONELLO - SUMMER '87 - )                             

!     PURPOSE.                                                          
!     --------                                                          

!        ESTIMATE THE ANALYSED FRICTION VELOCITY FROM THE DURATION AND   
!        THE "MEASURED" WINDSEA ENERGY. ESTIMATE THE ANALYSED WINDSEA   
!        MEAN FREQUENCY FROM THE "MEASURED" WINDSEA ENERGY.             

!**   INTERFACE.                                                        
!     ----------                                                        

!        *CALL* *FUSTAR (USMO, USA, EWFG, EWOI, EWA, T, FMWA, IJS, IJL)*

!        *USMO*   REAL     FIRST GUESS FRICTION VELOCITY.               
!        *USA*    REAL     ANALISED FRICTION VELOCITY.                  
!        *EWFG*   REAL     WINDSEA ENERGY (FIRST GUESS).                
!        *EWOI*   REAL     ESTIMATED WIND-SEA ENERGY FROM MEASUREMENTS. 
!        *EWA*    REAL     ANALYSED WIND-SEA ENERGY.                    
!        *T*      REAL     DURATION.                                    
!        *FMWA*   REAL     ANALYSED WINDSEA MEAN FREQUENCY.             
!        *IJS*    INTEGER  FIRST INDEX IN BLOCK.                        
!        *IJL*    INTEGER  LAST  INDEX IN BLOCK.                        

!      METHOD.                                                          
!      -------                                                          

!       THE DURATION OF THE WIND SEA IS USED TO COMPUTE THE FRICTION    
!       VELOCITY BY NEWTONS METHOD FROM THE MODEL DURATION CURVE AND    
!       THE MEASURED ENERGY.  (FIVE ITERATIONS ARE DONE)                
!       THIS VALUE IS USED TO UPDATE THE MODEL IF IT IS                 
!       IN REASONABLE AGREEMENT WITH THE MEASUREMENT FRICTION VELOCITY. 
!       IF THERE IS NO AGREEMENT THE MEASURED WIND SPEED IS USED        
!       TO COMPUTE THE WIND SEA.                                        

! ----------------------------------------------------------------------

      USE YOWPCONS , ONLY : G       ,EGRCRV   ,AGRCRV   ,BGRCRV   ,
     &             AFCRV   ,BFCRV
      USE YOWTABL  , ONLY : USTARM
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      INTEGER, PARAMETER :: NITER=5
      REAL,DIMENSION(IJS:IJL) :: USMO, USA, EWOI, EWA, T, FMWA,
     &                           EWFG, USTANAL , USTANAL0

! ----------------------------------------------------------------------

!     INLINE FUNCTION.                                                  
!     ----------------                                                  

!     DIMENSIONLESS FREQUENCY AS FUNCTION OF DIMENSIONLESS ENERGY       

      YNU(X) = (X/AFCRV)**(1./BFCRV)                                          

! ----------------------------------------------------------------------

!*    1. INITIALIZE OUTPUT FIELDS.                                      
!        -------------------------                                      

      DO IJ=IJS,IJL                                                
        USA(IJ)  =   -99.                                              
        EWA(IJ)  =  -999.                                              
        FMWA(IJ) = -9999.                                              
      ENDDO

!*    2. LOOP OVER GRID POINTS.                                         
!        ----------------------                                         

      G2=G*G
      BGRCRVM1=BGRCRV-1.
      EWFGMIN=0.0001

!*    2.1  DETERMINE FRICTION VELOCITY FROM DURATION GROWTH CURVE       
!          IF WINDSEA IS PRESENT.                                       
!           ------------------------------------------------------      
      DO IJ=IJS,IJL                                                
        IF (EWFG(IJ).GT.EWFGMIN.AND.T(IJ).GT.0.) THEN
          USTANAL0(IJ) = ABS(USMO(IJ))
          USTANAL(IJ) = USTANAL0(IJ) 
        ENDIF                                                          
      ENDDO

!     SOLVE ITERATIVELY THE RELATION FOR USTAR
!        ESTAR=EGRCRV*(TSTAR/(AGRCRV+TSTAR))**BGRCRV
!     USING THE 1 ORDER NEWTON METHOD Xn+1 = Xn - f(Xn)/f'(Xn)
!     WHERE f  = 0 IS A FUNCTION OBTAINED FROM THE PREVIOUS
!     RELATION IN TERM OF USTAR ONLY.

      DO ITER=1,NITER
        DO IJ=IJS,IJL
          IF (EWFG(IJ).GT.EWFGMIN.AND.T(IJ).GT.0.) THEN
            TSTAR = G*T(IJ)/USTANAL(IJ)
            XX = AGRCRV/TSTAR
            YY = 1.+XX
            ERATIO = (EWOI(IJ)*G2/USTANAL(IJ)**4)/EGRCRV

            USTANAL(IJ) = USTANAL(IJ)*(1.-(1.-ERATIO*YY**BGRCRV)
     &       /(4. - ERATIO*BGRCRV*XX*YY**BGRCRVM1))

            IF(USTANAL(IJ).LE.0.) THEN
              USTANAL0(IJ)=2.*USTANAL0(IJ)
              USTANAL(IJ)=USTANAL0(IJ)
            ENDIF
          ENDIF
        ENDDO
      ENDDO

!          -----------------------------------------------------
      DO IJ=IJS,IJL
        IF (EWFG(IJ).GT.EWFGMIN.AND.T(IJ).GT.0.) THEN
          IF(ABS(USTANAL(IJ)-USMO(IJ))/USMO(IJ).LT.0.5) THEN
            USA(IJ) = MIN(USTANAL(IJ),USTARM)
            EWA(IJ) = EWOI(IJ)
            FMWA(IJ) = YNU(EWOI(IJ)*G2/USTANAL(IJ)**4)*G/USTANAL(IJ)
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE FUSTAR
