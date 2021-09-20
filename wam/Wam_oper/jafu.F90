INTEGER(KIND=JWIM) FUNCTION JAFU (CL, J, IAN)

! ----------------------------------------------------------------------

!**** *JAFU* - FUNCTION TO COMPUTE THE INDEX ARRAY FOR THE
!              ANGLES OF THE INTERACTING WAVENUMBERS.

!     S. HASSELMANN        MPIFM        01/12/1985.

!*    PURPOSE.
!     --------

!       INDICES DEFINING BINS IN FREQUENCY AND DIRECTION PLANE INTO
!       WHICH NONLINEAR ENERGY TRANSFER INCREMENTS ARE STORED. NEEDED
!       FOR COMPUTATION OF THE NONLINEAR ENERGY TRANSFER.

!**   INTERFACE.
!     ----------

!       *FUNCTION* *JAFU (CL, J, IAN)*
!          *CL*  - WEIGHTS.
!          *J*   - INDEX IN ANGULAR ARRAY.
!          *IAN* - NUMBER OF ANGLES IN ARRAY.

!     METHOD.
!     -------

!       SEE REFERENCE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!        S. HASSELMANN AND K. HASSELMANN,JPO, 1985 B.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: J, IAN
      INTEGER(KIND=JWIM) :: IDPH, JA
      REAL(KIND=JWRB):: CL

      IDPH = CL
      JA = J+IDPH
      IF (JA <= 0)   JA = IAN+JA-1
      IF (JA >= IAN) JA = JA-IAN+1
      JAFU = JA

END FUNCTION JAFU
