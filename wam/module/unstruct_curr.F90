MODULE UNSTRUCT_CURR
  !
  ! Rule of conversion of angle in WAM.
  ! theta_{WAM} = 90 - theta_{trigonometric}
  !
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      IMPLICIT NONE
      REAL(KIND=JWRU), ALLOCATABLE :: CURTXY (:,:)

END MODULE
