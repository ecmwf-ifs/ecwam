FUNCTION IRANGE(X0,X1,DX) RESULT(IX)

  ! ----------------------------------------------------------------------------
  !
  !  1. Purpose :
  !
  !         Generate a sequence of linear-spaced integer numbers.
  !         Used for instance array addressing (indexing).
  ! ----------------------------------------------------------------------------
  !
  !     INTERFACE VARIABLES.
  !     --------------------
  
  !     ORIGIN.
  !     ----------
  !     Adapted from Babanin Young Donelan & Banner (BYDB) physics 
  !     as implemented as ST6 in WAVEWATCH-III 
  !     WW3 module:       W3SRC6MD    
  !     WW3 subroutine:   IRANGE
  !     Implementation into ECWAM DECEMBER 2021 by J. Kousal 
  
  ! ----------------------------------------------------------------------------
  
        USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
        USE YOMHOOK  , ONLY : LHOOK   ,DR_HOOK, JPHOOK
  
  ! ----------------------------------------------------------------------------
  
        IMPLICIT NONE
  
        INTEGER(KIND=JWIM),              INTENT(IN)     :: X0, X1, DX
  
        INTEGER(KIND=JWIM), ALLOCATABLE                 :: IX(:)
        INTEGER(KIND=JWIM)                              :: N, I
  
        REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
  ! ----------------------------------------------------------------------------
  !
  
        IF (LHOOK) CALL DR_HOOK('IRANGE',0,ZHOOK_HANDLE)
  
        N = INT(REAL(X1-X0)/REAL(DX))+1
        ALLOCATE(IX(N))
        DO I = 1, N
           IX(I) = X0+ (I-1)*DX
        END DO
  
        IF (LHOOK) CALL DR_HOOK('IRANGE',1,ZHOOK_HANDLE)
  
        END FUNCTION IRANGE
  