      SUBROUTINE LOC2GLO (ALOC, GLO, IJS, IJL, IG)
! ----------------------------------------------------------------------

!**** *LOC2GLO* - TRANSFERS A LOCAL BLOCK ARRAY TO A GLOBAL BLOCK ARRAY 

!     J. BIDLOT   ECMWF       MARCH 2000 

!*    PURPOSE.
!     --------

!    TRANSFERS A LOCAL BLOCK ARRAY TO A GLOBAL BLOCK ARRAY. 
!    IT WILL ALSO IMPOSE THE ICE MASK IF PRESENT


!**   INTERFACE.
!     ----------

!       *CALL* *LOC2GLO (ALOC, GLO, IJS, IJL, IG)
!         *ALOC*   - LOCAL BLOCK ARRAY.
!         *GLO*    - LOCAL BLOCK ARRAY.
!         *IJS*    - INDEX OF FIRST GRIDPOINT.
!         *IJL*    - INDEX OF LAST GRIDPOINT.
!         *IG*     - BLOCK NUMBER.


! ----------------------------------------------------------------------

      USE YOWICE   , ONLY : LICERUN  ,ICEMASK
      USE YOWPARAM , ONLY : NIBLO
      USE YOWPCONS , ONLY : ZMISS

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IJ, IJS, IJL, IG
      REAL :: ALOC(IJS:IJL), GLO(NIBLO)

      IF (LICERUN) THEN
        DO IJ = IJS,IJL
          IF (ICEMASK(IJ,IG).EQ.0) THEN
             GLO(IJ) = ZMISS
          ELSE
             GLO(IJ) = ALOC(IJ)
          ENDIF
        ENDDO
      ELSE
        DO IJ = IJS,IJL
          GLO(IJ) = ALOC(IJ)
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE LOC2GLO
