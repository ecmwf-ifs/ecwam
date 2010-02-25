      SUBROUTINE FDUR (EWFG, T, USMO, IJS, IJL)                         

!--------------------------------------------------------------------   

!**** *FDUR* - EVALUATE THE WINDSEA DURATION FROM THE MODEL GROTH CURVE.

!     P.LIONELLO     ECMWF       APRIL 1990                             

!     PURPOSE.                                                          
!     --------                                                          

!       EVALUATE THE WINDSEA DURATION FROM THE MODEL GROTH CURVE.       

!**   INTERFACE.                                                        
!     ----------                                                        

!       *CALL* *FDUR (EWFG, T, USMO, IJS, IJL)*                         

!        *EWFG*   REAL     WINDSEA ENERGY (FIRST GUESS).                
!        *T*      REAL     WINDSEA DURATION (OUTPUT).                   
!        *USMO*   REAL     FIRST GUESS USTAR.                           
!        *IJS*    INTEGER  FIRST INDEX IN BLOCK.                        
!        *IJL*    INTEGER  LAST  INDEX IN BLOCK.                        

!     METHOD.                                                           
!     -------                                                           

!        THE DURATION IS DERIVED BY THE NON DIMENSIONAL GROWTH CURVE

!        ESTAR=EGRCRV*(TSTAR/(AGRCRV+TSTAR))**BGRCRV

! ----------------------------------------------------------------------

      USE YOWPCONS , ONLY : G       ,EGRCRV   ,AGRCRV   ,BGRCRV

! ----------------------------------------------------------------------

      REAL,DIMENSION(IJS:IJL) :: EWFG, T, USMO 


!*    1. INITIALIZE DURATION.                                           
!        -------------------                                            

      DO IJ=IJS,IJL
        T(IJ) = -9.
      ENDDO

      G2=G*G
      BGRCRV_INV=1./BGRCRV

!*    2. LOOP OVER POINTS IN BLOCK.
!        --------------------------

      DO IJ=IJS,IJL
        IF (EWFG(IJ).GT.0.) THEN

!         DETERMINE THE EFFECTIVE DURATION OF THE WINDSEA ACCORDING TO
!         THE WAM MODEL GROWTH CURVE, IF WINDSEA IS PRESENT.

          ESTAR = EWFG(IJ)*G2/USMO(IJ)**4
          ERATIO= ESTAR/EGRCRV
          ERATIO= MIN (ERATIO, .93)
          ERATIO= ERATIO**BGRCRV_INV
          TSTAR = ERATIO*AGRCRV/(1.-ERATIO) 
          T(IJ) = USMO(IJ)/G*TSTAR
        ENDIF
      ENDDO
                                                                        
      RETURN
      END SUBROUTINE FDUR
