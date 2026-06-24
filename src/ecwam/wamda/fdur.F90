      SUBROUTINE FDUR (KIJS, KIJL, DEPTH, EWFG, T, USMO)

!--------------------------------------------------------------------   

!**** *FDUR* - EVALUATE THE WINDSEA DURATION FROM THE MODEL GROTH CURVE.

!     P.LIONELLO     ECMWF       APRIL 1990                             

!     PURPOSE.                                                          
!     --------                                                          

!       EVALUATE THE WINDSEA DURATION FROM THE MODEL GROTH CURVE.       

!**   INTERFACE.                                                        
!     ----------                                                        

!       *CALL* *FDUR (KIJS, KIJL, DEPTH, EWFG, T, USMO)*                         

!        *KIJS*   FIRST INDEX IN BLOCK.                        
!        *KIJL*   LAST  INDEX IN BLOCK.                        
!        *DEPTH*  WATER DEPTH
!        *EWFG*   WINDSEA ENERGY (FIRST GUESS).                
!        *T*      WINDSEA DURATION (OUTPUT).                   
!        *USMO*   FIRST GUESS USTAR.                           

!     METHOD.                                                           
!     -------                                                           

!        THE DURATION IS DERIVED BY THE NON DIMENSIONAL GROWTH CURVE

!        ESTAR=EGRCRV*(TSTAR/(AGRCRV+TSTAR))**BGRCRV

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWALTAS , ONLY : EGRCRV   ,AGRCRV   ,BGRCRV    ,             &
     &                      ESH      ,ASH     ,BSH
      USE YOWPCONS , ONLY : G

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: EWFG, USMO 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: T


      INTEGER(KIND=JWIM) :: IJ

      REAL(KIND=JWRB) :: G2, BGRCRV_INV, ESTAR, DSTAR, TSTAR, ERATIO
      REAL(KIND=JWRB) :: EMAX, X
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


!     INLINE FUNCTION.                                                  
!     ----------------                                                  

!     DIMENSIONLESS ENERGY LIMIT AS FUNCTION OF DIMENSIONLESS DEPTH 

      EMAX(X) = ESH*TANH(ASH*X**BSH)

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('FDUR',0,ZHOOK_HANDLE)

      G2=G*G
      BGRCRV_INV=1.0_JWRB/BGRCRV

!*    2. LOOP OVER POINTS IN BLOCK.
!        --------------------------

      DO IJ=KIJS,KIJL
        IF (EWFG(IJ) > 0._JWRB) THEN

!         DETERMINE THE EFFECTIVE DURATION OF THE WINDSEA ACCORDING TO
!         THE WAM MODEL GROWTH CURVE, IF WINDSEA IS PRESENT.

          ESTAR = EWFG(IJ)*G2/USMO(IJ)**4

          DSTAR = DEPTH(IJ)*G/USMO(IJ)**2
          ESTAR = MIN(EMAX(DSTAR),ESTAR)

          ERATIO= ESTAR/EGRCRV
!!!!why? May not be necessary once we use emax
          ERATIO= MIN (ERATIO, 0.93_JWRB)

          ERATIO= ERATIO**BGRCRV_INV
          TSTAR = ERATIO*AGRCRV/(1.0_JWRB-ERATIO) 
          T(IJ) = USMO(IJ)/G*TSTAR
        ELSE
          T(IJ) = -9.0_JWRB
        ENDIF
      ENDDO
                                                                        
      IF (LHOOK) CALL DR_HOOK('FDUR',1,ZHOOK_HANDLE)

      END SUBROUTINE FDUR
