      SUBROUTINE UPDATE (KIJS, KIJL, FL1,   &
     &                   WAVNUM, DEPTH,     &
     &                   ETOI, USMO, THMO)

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

!       *CALL* *UPDATE (KIJS, KIJL, FL1, WAVNUM, DEPTH, ETOI, USMO, THMO)*

!         *KIJS*    INTEGER   FIRST INDEX IN BLOCK.
!         *KIJL*    INTEGER   LAST  INDEX IN BLOCK.
!         *FL1*     REAL      SPECTRUM ( FIRST GUESS IN INPUT, ANALYSIS IN OUTPUT ).
!         *WAVNUM*  REAL      WAVE NUMBER.
!         *DEPTH*   REAL      WATER DEPTH.
!         *ETOI*    REAL      WAVE ENERGY (FROM O.I.).
!         *USMO*    REAL      USTAR : FIRST GUESS IN INPUT BUT ANALYSIS IN OUTPUT  !!!! 
!         *THMO*    REAL      WIND DIRECTION.

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

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "fdur.intfb.h"
#include "fkmean.intfb.h"
#include "fustar.intfb.h"
#include "fwsea.intfb.h"
#include "upwspec.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: ETOI, THMO

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: USMO


      INTEGER(KIND=JWIM) :: IJ

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: EMN, FMN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: F1MN, AKMN, XKMN 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: EWFG, EWOI, EWA, FMWFG
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: FMWA, TDUR, USA

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

      IF (LHOOK) CALL DR_HOOK('UPDATE',0,ZHOOK_HANDLE)

!*    1. FIND THE WINDSEA.                                              
!        -----------------                                              

      CALL FKMEAN(KIJS, KIJL, FL1, WAVNUM, EMN, FMN, F1MN, AKMN, XKMN)

      CALL FWSEA (KIJS, KIJL, FL1, DEPTH,  &
     &            EMN, ETOI, USMO, THMO,   &
     &            EWFG, FMWFG)

!*    2. CORRECT OVERESTIMATE OF FIRST GUESS WINDSEA ENERGY.            
!        ---------------------------------------------------            

      DO IJ = KIJS,KIJL
        EWFG(IJ) = MIN(EWFG(IJ),EMN(IJ))
      ENDDO
                                                                        
! ----------------------------------------------------------------------

!*    3. FIND THE DURATION.                                             
!        ------------------                                             

      CALL FDUR (KIJS, KIJL, DEPTH, EWFG, TDUR, USMO)
                                                                        
! ----------------------------------------------------------------------

!*    4. THE RATIO WINDSEA/SWELL IS ASSUMED CORRECT.                    
!        -------------------------------------------                    

      DO IJ=KIJS,KIJL
        IF (EMN(IJ) > 0.0_JWRB) EWOI(IJ) = EWFG(IJ)/EMN(IJ)*ETOI(IJ)
      ENDDO

! ----------------------------------------------------------------------

!*    5. FIND THE NEW USTAR AND THE NEW WINDSEA MEAN FREQUENCY.         
!        ------------------------------------------------------         

      CALL FUSTAR (KIJS, KIJL, DEPTH, USMO, USA, EWFG, EWOI, EWA, TDUR, FMWA)

! ----------------------------------------------------------------------

!*    6. THE ASSIMILATION IS ADJUSTED ON THE DOMINANT PART              
!*       OF THE SPECTRUM.                                               
!        -------------------------------------------------              

!        THE FOLLOWING RETURN CODE ARE POSSIBLE :                       
!          THERE IS NO WINDSEA ->                                       
!                 EWFG  =  -999.                                        
!                 FMWFG = -9999.                                        

      DO IJ=KIJS,KIJL
        IF (EWA(IJ) > 0.0_JWRB .AND. EWA(IJ) < 0.6_JWRB*ETOI(IJ)) THEN
          EWFG(IJ)  =  -999.0_JWRB
          FMWFG(IJ) = -9999.0_JWRB
        ENDIF
      ENDDO


!     IF THERE WAS A PROBLEM IN DETERMINING USTAR THEN
!     NO WINDSEA UPDATE SHOULD BE MADE
      DO IJ=KIJS,KIJL
        IF (USA(IJ) < 0.0_JWRB) THEN
          EWFG(IJ)  =  -999.0_JWRB
          FMWFG(IJ) = -9999.0_JWRB
        ENDIF
      ENDDO

!     UPDATE USTAR ONLY WHEN THE MEAN WINDSEA FREQUENCY AND
!     ANALYSED WAVE ENERGY ARE POSITIVE
!     I.E. WHEN THE UPDATE ON SPECTRA WILL ASSUME IT IS MAINLY WINDSEA

      DO IJ=KIJS,KIJL
        IF (USA(IJ) > 0.0_JWRB .AND. FMWFG(IJ) > 0.0_JWRB .AND. ETOI(IJ) > 0.0_JWRB) USMO(IJ) = USA(IJ)
      ENDDO
                                                                        
! ----------------------------------------------------------------------

!*    7. UPDATE THE SPECTRA.                                            
!        -------------------                                            

      CALL UPWSPEC (KIJS, KIJL, FL1, DEPTH, EMN, AKMN, ETOI, FMWFG, FMWA)
                                                                        
      IF (LHOOK) CALL DR_HOOK('UPDATE',1,ZHOOK_HANDLE)

      END SUBROUTINE UPDATE
