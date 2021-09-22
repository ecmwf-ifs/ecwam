SUBROUTINE NEWWIND (IJS, IJL, CDATE, CDATEWH,             &
 &                  LLNEWREAD, LLNEWFILE,                 &
 &                  FF_NOW, FF_NEXT,                      &
 &                  CGROUP,                               &
 &                  CIWA)
! ----------------------------------------------------------------------

!**** *NEWWIND* - HANDLING OF WINDS OUTSIDE MULTITASKED AREA    

!     P.A.E.M. JANSSEN  KNMI/ECMWF  SEPTEMBER 1994
!     J. BIDLOT         ECMWF       FEBRUARY 1996  MESSAGE PASSING 
!     S. ABDALLA        ECMWF       OCTOBER 2001   INCLUSION OF AIR
!                                                  DENSITY AND Zi/L
!     J. BIDLOT         ECMWF       AUGUST 2008: CLEAN-UP 

!*    PURPOSE.
!     --------

!       GETS FORCING FIELDS WHEN NEEDED.                                   

!**   INTERFACE.
!     ----------

!       *CALL* *NEWWIND (IJS, IJL, CDATE, LLNEWREAD, LLNEWFILE,
!                        FF_NOW, FF_NEXT,
!                        CGROUP,
!                        CICOVER, CITHICK, CIWA)
!      *IJS*     - INDEX OF FIRST GRIDPOINT
!      *IJL*     - INDEX OF LAST GRIDPOINT
!      *CDATE*   - START DATE OF SOURCE FUNCTION INTEGRATION
!      *LLNEWREAD* TRUE IF NEW WINDS HAVE BEEN READ
!      *LLNEWFILE* TRUE IF NEW WIND FILE HAS BEEN OPENED
!      *FF_NOW*  - DATA STRUCTURE WITH THE CURRENT FORCING FIELDS
!      *FF_NEXT* - DATA STRUCTURE WITH THE NEXT FORCING FIELDS
!      *CGROUP*  - GROUP SPEED.
!      *CIWA*    - SEA ICE WAVE ATTENUATION FACTOR.

!     METHOD.
!     -------

!       READ NEW WINDS AND ASSOCIATED FIELDS WHEN NEEDED.

!     EXTERNALS.
!     ---------

!       *INCDATE*   - UPDATE DATE TIME GROUP.

!     REFERENCE.
!     ----------

!       NONE

! ----------------------------------------------------------------------

!*    *PARAMETER*  FOR ARRAY DIMENSIONS.

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : FORCING_FIELDS

      USE YOWPCONS , ONLY : ACD      ,BCD      ,EPSMIN
      USE YOWCOUP  , ONLY : LWCOU
      USE YOWPARAM , ONLY : NFRE
      USE YOWPHYS  , ONLY : ALPHA
      USE YOWSTAT  , ONLY : IDELWO   ,NPROMA_WAM
      USE YOWTEST  , ONLY : IU06
      USE YOWWIND  , ONLY : CDATEWL  ,CDAWIFL  ,CDATEFL  ,CDTNEXT  ,    &
     &            NSTORE   ,WSPMIN_RESET_TAUW  ,USTMIN_RESET_TAUW
      USE YOWWNDG  , ONLY : ICODE    ,ICODE_CPL

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "cireduce.intfb.h"
#include "incdate.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      CHARACTER(LEN=14), INTENT(IN) :: CDATE
      CHARACTER(LEN=14), INTENT(INOUT) :: CDATEWH
      LOGICAL, INTENT(INOUT) :: LLNEWREAD, LLNEWFILE
      TYPE(FORCING_FIELDS), DIMENSION(IJS:IJL), INTENT(INOUT) :: FF_NOW
      TYPE(FORCING_FIELDS), DIMENSION(IJS:IJL), INTENT(IN) :: FF_NEXT
      REAL(KIND=JWRB),DIMENSION(IJS:IJL,NFRE), INTENT(IN) :: CGROUP
      REAL(KIND=JWRB),DIMENSION(IJS:IJL,NFRE), INTENT(INOUT) :: CIWA


      INTEGER(KIND=JWIM) :: ICODE_WND
      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA, IJ

      REAL(KIND=JWRB) :: WGHT, TLWMAX
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('NEWWIND',0,ZHOOK_HANDLE)

ASSOCIATE(WSWAVE => FF_NOW%WSWAVE, &
 &        WDWAVE => FF_NOW%WDWAVE, &
 &        UFRIC => FF_NOW%UFRIC, &
 &        TAUW => FF_NOW%TAUW, &
 &        CHNK => FF_NOW%CHNK, &
 &        AIRD => FF_NOW%AIRD, &
 &        WSTAR => FF_NOW%WSTAR, &
 &        CICOVER => FF_NOW%CICOVER, &
 &        CITHICK => FF_NOW%CITHICK)

      NPROMA=NPROMA_WAM

      IF (LWCOU) THEN
        ICODE_WND = ICODE_CPL
      ELSE
        ICODE_WND = ICODE
      ENDIF

      WGHT=1.0_JWRB/MAX(WSPMIN_RESET_TAUW,EPSMIN)

      IF (CDATE >= CDATEWH) THEN

!*    2. NEW WIND INPUT.
!        ---------------
        IF (CDATE >= CDATEFL) THEN
            LLNEWFILE = .TRUE.
        ENDIF

!*    2.2 NEW WINDS ARE READ IN.
!         ----------------------

        CDATEWL = CDTNEXT

        CALL GSTATS(1492,0)
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ,TLWMAX)
        DO JKGLO=IJS,IJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL)
          IF (ICODE_WND == 3 ) THEN
            DO IJ = KIJS, KIJL
              WSWAVE(IJ)=FF_NEXT(IJ)%WSWAVE
! adapt first estimate of wave induced stress for low winds
! to a fraction of the simple relation u*^2 = Cd(U10) * U10^2
! where this fraction varies from 0 for U10=0 to 1 for U10=WSPMIN_RESET_TAUW
              IF (WSWAVE(IJ) < WSPMIN_RESET_TAUW) THEN
                TLWMAX=WGHT*(ACD+BCD*WSWAVE(IJ))*WSWAVE(IJ)**3
                TAUW(IJ)=MIN(TAUW(IJ),TLWMAX)
              ENDIF
            ENDDO
          ELSE
            DO IJ = KIJS, KIJL
              UFRIC(IJ)=FF_NEXT(IJ)%UFRIC
! update the estimate of TAUW
              TAUW(IJ)=UFRIC(IJ)**2*(1.0_JWRB-(ALPHA/CHNK(IJ))**2)
! adapt first estimate of wave induced stress for low winds
              IF (UFRIC(IJ) < USTMIN_RESET_TAUW) TAUW(IJ)=0.0_JWRB
            ENDDO
          ENDIF

          DO IJ = KIJS, KIJL
            WDWAVE(IJ)=FF_NEXT(IJ)%WDWAVE
            AIRD(IJ)=FF_NEXT(IJ)%AIRD
            WSTAR(IJ)=FF_NEXT(IJ)%WSTAR
            CICOVER(IJ)=FF_NEXT(IJ)%CICOVER
            CITHICK(IJ)=FF_NEXT(IJ)%CITHICK
          ENDDO

        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1492,1)


        CALL INCDATE(CDATEWH, IDELWO)
        LLNEWREAD = .TRUE.   

!       UPDATE THE SEA ICE REDUCTION FACTOR
        CALL CIREDUCE (IJS, IJL, CGROUP, CICOVER, CITHICK, CIWA)

      ENDIF

END ASSOCIATE
      IF (LHOOK) CALL DR_HOOK('NEWWIND',1,ZHOOK_HANDLE)

END SUBROUTINE NEWWIND
