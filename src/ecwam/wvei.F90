! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

      FUNCTION WVEI(x) RESULT(res)

! ----------------------------------------------------------------------------
!
!  1. Purpose :
!
      ! Real-valued exponential integral Ei(x) approximation for real x
      ! Uses power series for small/negative x and an asymptotic expansion for large positive x.
      ! Reasonable accuracy for typical geophysical ranges; replace with a library routine if you need higher precision.

!----------------------------------------------------------------------
!
!     INTERFACE VARIABLES.
!     --------------------

!     ORIGIN.
!     ----------
!     Implementation into ECWAM DECEMBER 2025 by J. Kousal based on available online solvers

! ----------------------------------------------------------------------------
!

        USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
        USE YOWPCONS    , ONLY : GAMMA_E
        USE YOMHOOK  , ONLY : LHOOK   ,DR_HOOK, JPHOOK

!----------------------------------------------------------------------

        IMPLICIT NONE
        REAL(KIND=JWRB), INTENT(IN)  :: x
        REAL(KIND=JWRB) :: term, sum, res
        INTEGER(KIND=JWIM) :: k, kmax
        REAL(KIND=JWRB) :: eps

        REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------------
!

      IF (LHOOK) CALL DR_HOOK('WVEI',0,ZHOOK_HANDLE)

            
      eps = 1.0E-12_JWRB
      if (x < 0.0_JWRB .or. abs(x) <= 6.0_JWRB) then
            ! series: Ei(x) = gamma + ln(|x|) + sum_{k=1..inf} x^k/(k*k!)
            if (x == 0.0_JWRB) then
                  res = -huge(1.0_JWRB)  ! singular; caller should not hit exact zero normally
                  return
            end if
            sum = 0.0_JWRB
            term = x
            k = 1
            kmax = 200
            do while (k <= kmax)
                  sum = sum + term/real(k,kind=JWRB)
                  term = term * x/real(k+1,kind=JWRB)
                  if (abs(term/real(k+1,kind=JWRB)) < abs(sum)*eps) exit
                  k = k + 1
            end do
            res = GAMMA_E + log(abs(x)) + sum
      else
            ! asymptotic for large positive x: Ei(x) ~ exp(x)/x * (1 + 1/x + 2!/x^2 + 6/x^3 + ...)
            kmax = 50
            sum = 1.0_JWRB
            term = 1.0_JWRB
            do k = 1, kmax
                  term = term * real(k,kind=JWRB) / x
                  sum = sum + term
                  if (abs(term) < abs(sum)*eps) exit
            end do
            res = exp(x) / x * sum
      end if

      IF (LHOOK) CALL DR_HOOK('WVEI',1,ZHOOK_HANDLE)

      END FUNCTION WVEI            