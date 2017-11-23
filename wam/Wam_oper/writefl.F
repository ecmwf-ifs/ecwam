      SUBROUTINE WRITEFL(FL, IJINF, IJSUP, KINF, KSUP, MINF, MSUP,
     &                   FILENAME,IUNIT,LOUNIT,LRSTPARAL)

! ---------------------------------------------------------------------
!     J. BIDLOT    ECMWF      MARCH 1997

!*    PURPOSE.
!     --------
!     WRITES ARRAY FL TO FILE.

!**   INTERFACE.
!     ----------
!     CALL *WRITEFL(FL, IJINF, IJSUP, KINF, KSUP, MINF, MSUP,
!    &              FILENAME,IUNIT,LOUNIT,LRSTPARAL)
!     *FL*        ARRAY TO BE WRITTEN TO FILE
!     *IJINF*     FIRST LOWER DIMENSION OF FL
!     *IJSUP*     FIRST UPPER DIMENSION OF FL
!     *KINF*      SECOND LOWER DIMENSION BOUND OF FL
!     *KSUP*      SECOND UPPER DIMENSION BOUND OF FL
!     *MINF*      THIRD LOWER DIMENSION BOUND OF FL
!     *MSUP*      THIRD UPPER DIMENSION BOUND OF FL
!     *FILENAME*  FILENAME (INCLUDING PATH) OF TARGET FILE
!     *IUNIT*     PBIO UNIT (ONLY ACTIVE IF PBIO OUTPUT)
!     *LOUNIT*    LOGICAL, TRUE IF FREE UNIT HAS TO BE FOUND
!     *LRSTPARAL* LOGICAL, TRUE THEN WRITE IN PARALLEL

!     METHOD.
!     -------
!     WRITES ARRAY FL TO FILE (A FORTRAN WRITE)

!     FL IS WRITTEN AS UNFORMATTED BINARY TO FILENAME CONNECTED
!     TO UNIT IUNIT.
!     IT IS IN A FORM THAT IS INDEPENDENT OF THE MODEL DECOMPOSITION.



!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!     NONE

! ----------------------------------------------------------------------

      USE YOWMPP   , ONLY : NPROC
      USE YOWPARAM , ONLY : LL1D
      USE YOWSPEC  , ONLY : IJ2NEWIJ
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWUNPOOL, ONLY : LLUNSTR
      USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "unblkrord.intfb.h"

      INTEGER, INTENT(IN) :: IJINF, IJSUP, KINF, KSUP, MINF, MSUP
      INTEGER, INTENT(INOUT) :: IUNIT

      REAL, DIMENSION(IJINF:IJSUP,KINF:KSUP,MINF:MSUP), 
     &     INTENT(INOUT) :: FL

      CHARACTER(LEN=296), INTENT(IN) :: FILENAME

      LOGICAL, INTENT(IN) :: LOUNIT, LRSTPARAL

      INTEGER :: IJ, K, M, J2, J3
      INTEGER :: KRET, KOUNT, LFILE
      INTEGER :: I_GET_UNIT

      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION(IJINF:IJSUP,KINF:KSUP,MINF:MSUP) :: FL_G
      REAL, ALLOCATABLE :: RFL(:,:,:)

! ----------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WRITEFL',0,ZHOOK_HANDLE)
#endif

      LFILE=0
      IF (FILENAME.NE. ' ') LFILE=LEN_TRIM(FILENAME)
      IF(LOUNIT) THEN
        IUNIT=I_GET_UNIT(IU06, FILENAME(1:LFILE) , 'w', 'u', 0)
      ELSE
        IUNIT=I_GET_UNIT(IU06, FILENAME(1:LFILE) , 'a', 'u', 0)
      ENDIF

      IF(LLUNSTR .AND. .NOT.LRSTPARAL) THEN
        FL_G(0,:,:)=0.
        CALL UNBLKRORD(1,IJINF,IJSUP,KINF,KSUP,MINF,MSUP,
     &                 FL(IJINF:IJSUP,KINF:KSUP,MINF:MSUP),
     &               FL_G(IJINF:IJSUP,KINF:KSUP,MINF:MSUP))
        WRITE(IUNIT) (((FL_G(IJ,J2,J3),
     &                  IJ=IJINF,IJSUP),
     &                  J2=KINF,KSUP),
     &                  J3=MINF,MSUP)

      ELSE IF(LRSTPARAL .OR. LL1D .OR. NPROC.EQ.1) THEN
        WRITE(IUNIT) (((FL(IJ,J2,J3),
     &                  IJ=IJINF,IJSUP),
     &                  J2=KINF,KSUP),
     &                  J3=MINF,MSUP)
      ELSE
!     WHEN 2-D DECOMPOSITION IS USED THEN THE INDEXES IJ ARE RE-LABELLED
!     BUT THE SINGLE BINARY INPUT FILES SHOULD BE IN THE OLD MAPPING
        WRITE(IUNIT) (((FL(IJ2NEWIJ(IJ),J2,J3),
     &                  IJ=IJINF,IJSUP),
     &                  J2=KINF,KSUP),
     &                  J3=MINF,MSUP)

      ENDIF

      IF (ITEST.GE.2)
     &  WRITE(IU06,*) ' SUB. WRITEFL: OUTPUT OF FL DONE',
     &                 KINF,KSUP,MINF,MSUP

      CLOSE(IUNIT)

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WRITEFL',1,ZHOOK_HANDLE)
#endif
      RETURN
      END SUBROUTINE WRITEFL
