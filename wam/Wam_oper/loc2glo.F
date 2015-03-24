      SUBROUTINE LOC2GLO (ALOC, GLO, IJS, IJL, IJ_OFFSET, IG, THRSCI)
! ----------------------------------------------------------------------

!**** *LOC2GLO* - TRANSFERS A LOCAL BLOCK ARRAY TO A GLOBAL BLOCK ARRAY 

!     J. BIDLOT   ECMWF       MARCH 2000 

!*    PURPOSE.
!     --------

!    TRANSFERS A LOCAL BLOCK ARRAY TO A GLOBAL BLOCK ARRAY. 
!    IT WILL ALSO IMPOSE THE ICE MASK IF PRESENT


!**   INTERFACE.
!     ----------

!       *CALL* *LOC2GLO (ALOC, GLO, IJS, IJL, IJS_OFFSET, IG, THRSCI)
!         *ALOC*   - LOCAL BLOCK ARRAY.
!         *GLO*    - LOCAL BLOCK ARRAY.
!         *IJS*    - INDEX OF FIRST GRIDPOINT.
!         *IJL*    - INDEX OF LAST GRIDPOINT.
!         *IJ_OFFSET* OFFSET TO POINT TO GLOBAL DATA BLOCK
!         *IG*     - BLOCK NUMBER.
!         *THRSCI* - SEA ICE THRESHOLD FOR GIVEN PARAMETER.

! ----------------------------------------------------------------------

      USE YOWICE   , ONLY : LICERUN  , LMASKICE, CICOVER, CITHRSH
      USE YOWMPP   , ONLY : IRANK
      USE YOWPARAM , ONLY : NIBLO
      USE YOWPCONS , ONLY : ZMISS

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IJ, IJS, IJL, IJ_OFFSET, IG

      REAL :: THRSCI
      REAL :: ALOC(IJS:IJL), GLO(NIBLO)

! ----------------------------------------------------------------------

      IF (LICERUN) THEN
        DO IJ = IJS,IJL
          IF (CICOVER(IJ,IG).GT.THRSCI) THEN
             GLO(IJ+IJ_OFFSET) = ZMISS
          ELSE
             GLO(IJ+IJ_OFFSET) = ALOC(IJ)
          ENDIF
        ENDDO
      ELSE
        DO IJ = IJS,IJL
          GLO(IJ+IJ_OFFSET) = ALOC(IJ)
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE LOC2GLO
