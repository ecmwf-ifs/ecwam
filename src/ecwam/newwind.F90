SUBROUTINE NEWWIND (CDATE, CDATEWH, LLNEWFILE,           &
 &                  WVPRPT, FF_NOW, FF_NEXT)
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

!       *CALL* *NEWWIND (CDATE, CDATEWH, LLNEWFILE,
!                        WVPRPT, FF_NOW, FF_NEXT,
!      *CDATE*    - START DATE OF SOURCE FUNCTION INTEGRATION
!      *CDATEWH*  - DATE OF THE NEXT FORCING FIELDS
!      *LLNEWFILE*- TRUE IF NEW WIND FILE HAS BEEN OPENED
!      *WVPRPT*   - WAVE PROPERTIES FIELDS
!      *FF_NOW*   - DATA STRUCTURE WITH THE CURRENT FORCING FIELDS
!      *FF_NEXT*  - DATA STRUCTURE WITH THE NEXT FORCING FIELDS

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
      USE YOWDRVTYPE  , ONLY : FREQUENCY, FORCING_FIELDS

      USE YOWPCONS , ONLY : ACD      ,BCD      ,EPSMIN
      USE YOWCOUP  , ONLY : LWCOU
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWPARAM , ONLY : NFRE
      USE YOWPHYS  , ONLY : ALPHA
      USE YOWSTAT  , ONLY : IDELWO
      USE YOWTEST  , ONLY : IU06
      USE YOWWIND  , ONLY : CDATEWL  ,CDAWIFL  ,CDATEFL  ,CDTNEXT  ,         &
     &                      NSTORE   ,WSPMIN_RESET_TAUW  ,USTMIN_RESET_TAUW
      USE YOWWNDG  , ONLY : ICODE    ,ICODE_CPL

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "cireduce.intfb.h"
#include "incdate.intfb.h"

      CHARACTER(LEN=14), INTENT(IN)    :: CDATE
      CHARACTER(LEN=14), INTENT(INOUT) :: CDATEWH
      LOGICAL, INTENT(INOUT)           :: LLNEWFILE
      TYPE(FREQUENCY), DIMENSION(NPROMA_WAM, NFRE, NCHNK), INTENT(INOUT) :: WVPRPT
      TYPE(FORCING_FIELDS), DIMENSION(NPROMA_WAM, NCHNK), INTENT(INOUT)  :: FF_NOW
      TYPE(FORCING_FIELDS), DIMENSION(NPROMA_WAM, NCHNK), INTENT(IN)     :: FF_NEXT


      INTEGER(KIND=JWIM) :: ICODE_WND
      INTEGER(KIND=JWIM) :: ICHNK, KIJS, KIJL, IJ

      REAL(KIND=JWRB) :: WGHT, TLWMAX
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('NEWWIND',0,ZHOOK_HANDLE)

ASSOCIATE(CGROUP => WVPRPT%CGROUP, &
 &        CIWA => WVPRPT%CIWA, &
 &        WSWAVE => FF_NOW%WSWAVE, &
 &        WDWAVE => FF_NOW%WDWAVE, &
 &        UFRIC => FF_NOW%UFRIC, &
 &        TAUW => FF_NOW%TAUW, &
 &        CHRNCK => FF_NOW%CHRNCK, &
 &        AIRD => FF_NOW%AIRD, &
 &        WSTAR => FF_NOW%WSTAR, &
 &        CICOVER => FF_NOW%CICOVER, &
 &        CITHICK => FF_NOW%CITHICK)

      IF (LWCOU) THEN
        ICODE_WND = ICODE_CPL
      ELSE
        ICODE_WND = ICODE
      ENDIF

      WGHT = 1.0_JWRB/MAX(WSPMIN_RESET_TAUW,EPSMIN)

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
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, KIJL, IJ, TLWMAX)
        DO ICHNK = 1, NCHNK
          KIJS = 1
          KIJL = NPROMA_WAM
          IF (ICODE_WND == 3 ) THEN
            DO IJ = KIJS, KIJL
              WSWAVE(IJ, ICHNK) = FF_NEXT(IJ, ICHNK)%WSWAVE
! adapt first estimate of wave induced stress for low winds
! to a fraction of the simple relation u*^2 = Cd(U10) * U10^2
! where this fraction varies from 0 for U10=0 to 1 for U10=WSPMIN_RESET_TAUW
              IF (WSWAVE(IJ, ICHNK) < WSPMIN_RESET_TAUW) THEN
                TLWMAX = WGHT * (ACD+BCD*WSWAVE(IJ, ICHNK)) * WSWAVE(IJ, ICHNK)**3
                TAUW(IJ, ICHNK) = MIN(TAUW(IJ, ICHNK), TLWMAX)
              ENDIF
            ENDDO
          ELSE
            DO IJ = KIJS, KIJL
              UFRIC(IJ, ICHNK) = FF_NEXT(IJ, ICHNK)%UFRIC
! update the estimate of TAUW
              TAUW(IJ, ICHNK) = UFRIC(IJ, ICHNK)**2 * (1.0_JWRB-(ALPHA/CHRNCK(IJ, ICHNK))**2)
! adapt first estimate of wave induced stress for low winds
              IF (UFRIC(IJ, ICHNK) < USTMIN_RESET_TAUW) TAUW(IJ, ICHNK) = 0.0_JWRB
            ENDDO
          ENDIF

          WDWAVE(KIJS:KIJL, ICHNK)  = FF_NEXT(KIJS:KIJL, ICHNK)%WDWAVE
          AIRD(KIJS:KIJL, ICHNK)    = FF_NEXT(KIJS:KIJL, ICHNK)%AIRD
          WSTAR(KIJS:KIJL, ICHNK)   = FF_NEXT(KIJS:KIJL, ICHNK)%WSTAR
          CICOVER(KIJS:KIJL, ICHNK) = FF_NEXT(KIJS:KIJL, ICHNK)%CICOVER
          CITHICK(KIJS:KIJL, ICHNK) = FF_NEXT(KIJS:KIJL, ICHNK)%CITHICK

        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1492,1)


        CALL INCDATE(CDATEWH, IDELWO)

!       UPDATE THE SEA ICE REDUCTION FACTOR
        CALL CIREDUCE (CGROUP, CICOVER, CITHICK, CIWA)

      ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('NEWWIND',1,ZHOOK_HANDLE)

END SUBROUTINE NEWWIND