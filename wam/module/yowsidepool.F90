module yowSidepool
  USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
  implicit none

  !> Number of sides, global
  !> Number of nodes neighbors?.
  integer(KIND=JWIM) :: ns_global = 0

  !> Number of sides, local
  integer(KIND=JWIM) :: ns = 0
end module yowSidepool
