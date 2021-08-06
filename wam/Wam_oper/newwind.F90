      SUBROUTINE NEWWIND (IJS, IJL, CDATE, CDATEWH,             &
     &                    NEWREAD, NEWFILE,                     &
     &                    U10OLD, THWOLD, U10NEW, THWNEW,       &
     &                    USOLD, USNEW,                         &
     &                    ROAIRO, ROAIRN, ZIDLOLD, ZIDLNEW,     &
     &                    CGROUP,                               &
     &                    CICOVER, CITHICK, CIWA,               &
     &                    TAUW, BETAOLD)
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

!       *CALL* *NEWWIND (IJS, IJL, CDATE, NEWREAD, NEWFILE,
!                        U10OLD,THWOLD,U10NEW,THWNEW, 
!                        USOLD, USNEW,
!                        ROAIRO, ROAIRN, ZIDLOLD,ZIDLNEW,
!                        CGROUP,
!                        CICOVER, CITHICK, CIWA,
!                        TAUW, BETAOLD)
!      *IJS*     - INDEX OF FIRST GRIDPOINT
!      *IJL*     - INDEX OF LAST GRIDPOINT
!      *CDATE*   - START DATE OF SOURCE FUNCTION INTEGRATION
!      *NEWREAD* - TRUE IF NEW WINDS HAVE BEEN READ
!      *NEWFILE* - TRUE IF NEW WIND FILE HAS BEEN OPENED
!      *U10NEW*  - NEW WIND SPEED IN M/S.
!      *U10OLD*  - INTERMEDIATE STORAGE OF MODULUS OF WIND
!                  VELOCITY.
!      *THWNEW*  - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!      *THWOLD*  - INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
!                  WIND VELOCITY.
!      *USOLD*   - INTERMEDIATE STORAGE OF FRICTIOn VELOCITY
!      *USNEW*   - NEW FRICTIOn VELOCITY
!      *ROAIRN*  - AIR DENSITY IN KG/M3.
!      *ROAIRO*  - INTERMEDIATE STORAGE OF AIR DENSITY.
!      *ZIDLNEW* - Zi/L (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!      *ZIDLOLD* - INTERMEDIATE STORAGE OF Zi/L.
!      *CGROUP*  - GROUP SPEED.
!      *CICOVER* - SEA ICE COVER. 
!      *CITHICK* - SEA ICE THICKNESS. 
!      *CIWA*    - SEA ICE WAVE ATTENUATION FACTOR.
!      *TAUW*    - WAVE STRESS
!      *BETAOLD* - CHARNOCK PARAMETER 

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

      USE YOWPCONS , ONLY : ACD      ,BCD      ,EPSMIN
      USE YOWCOUP  , ONLY : LWCOU
      USE YOWPARAM , ONLY : NFRE
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPHYS  , ONLY : ALPHA
      USE YOWSTAT  , ONLY : IDELWO   ,NPROMA_WAM
      USE YOWTEST  , ONLY : IU06
      USE YOWWIND  , ONLY : CDA      ,CDAWIFL  ,CDATEFL  ,CDTNEXT  ,    &
     &            NSTORE   ,FF_NEXT  ,WSPMIN_RESET_TAUW  ,              &
     &            USTMIN_RESET_TAUW
      USE YOWWNDG  , ONLY : ICODE    ,ICODE_CPL

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "cireduce.intfb.h"
#include "incdate.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      CHARACTER(LEN=14), INTENT(INOUT) :: CDATEWH
      LOGICAL, INTENT(INOUT) :: NEWREAD, NEWFILE
      REAL(KIND=JWRB),DIMENSION(NINF:NSUP), INTENT(IN) :: U10OLD, THWOLD
      REAL(KIND=JWRB),DIMENSION(NINF:NSUP), INTENT(IN) :: USOLD
      REAL(KIND=JWRB),DIMENSION(NINF:NSUP), INTENT(IN) :: ROAIRO, ZIDLOLD
      REAL(KIND=JWRB),DIMENSION(NINF:NSUP), INTENT(INOUT) :: TAUW
      REAL(KIND=JWRB),DIMENSION(NINF:NSUP), INTENT(IN)    :: BETAOLD 
      REAL(KIND=JWRB),DIMENSION(NINF:NSUP), INTENT(INOUT) :: U10NEW, THWNEW
      REAL(KIND=JWRB),DIMENSION(NINF:NSUP), INTENT(INOUT) :: USNEW
      REAL(KIND=JWRB),DIMENSION(NINF:NSUP), INTENT(INOUT) :: ROAIRN, ZIDLNEW
      REAL(KIND=JWRB),DIMENSION(IJS:IJL,NFRE), INTENT(IN) :: CGROUP
      REAL(KIND=JWRB),DIMENSION(NINF:NSUP), INTENT(INOUT) :: CICOVER, CITHICK 
      REAL(KIND=JWRB),DIMENSION(IJS:IJL,NFRE), INTENT(INOUT) :: CIWA
      CHARACTER(LEN=14), INTENT(IN) :: CDATE


      INTEGER(KIND=JWIM), SAVE :: ISTORE = 0

      INTEGER(KIND=JWIM) :: ICODE_WND
      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA, IJ

      REAL(KIND=JWRB) :: WGHT, TLWMAX
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('NEWWIND',0,ZHOOK_HANDLE)

      NPROMA=NPROMA_WAM

      IF (LWCOU) THEN
        ICODE_WND = ICODE_CPL
      ELSE
        ICODE_WND = ICODE
      ENDIF

      WGHT=1.0_JWRB/MAX(WSPMIN_RESET_TAUW,EPSMIN)

!*    1. WINDS ARE TAKEN FROM INTERMEDIATE STORAGE.
!        ------------------------------------------
      IF (CDATE < CDATEWH) THEN
        CALL GSTATS(1492,0)
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ)
        DO JKGLO=IJS,IJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL)
          DO IJ = KIJS, KIJL
            U10NEW(IJ) = U10OLD(IJ)
            USNEW(IJ) = USOLD(IJ)
            THWNEW(IJ) = THWOLD(IJ)
            ROAIRN(IJ) = ROAIRO(IJ)
            ZIDLNEW(IJ) = ZIDLOLD(IJ)
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1492,1)
      ELSE
!*    2. NEW WIND INPUT.
!        ---------------
        IF (CDATE >= CDATEFL) THEN
            NEWFILE = .TRUE.
            ISTORE=1
        ENDIF

!*    2.2 NEW WINDS ARE READ IN.
!         ----------------------
        IF (ISTORE > NSTORE) THEN
          WRITE(IU06,*) ' ********************************************'
          WRITE(IU06,*) '      FATAL ERROR IN SUB. NEWWIND:           '
          WRITE(IU06,*) '      ISTORE > NSTORE !!!                    '
          WRITE(IU06,*) '      ISTORE = ',ISTORE 
          WRITE(IU06,*) '      NSTORE = ',NSTORE 
          WRITE(IU06,*) ' ********************************************'
          CALL ABORT1
        ENDIF

        CDA=CDTNEXT(ISTORE)

        CALL GSTATS(1492,0)
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ,TLWMAX)
        DO JKGLO=IJS,IJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL)
          IF (ICODE_WND == 3 ) THEN
            DO IJ = KIJS, KIJL
              U10NEW(IJ)=FF_NEXT(IJ,ISTORE)%WSWAVE
! adapt first estimate of wave induced stress for low winds
! to a fraction of the simple relation u*^2 = Cd(U10) * U10^2
! where this fraction varies from 0 for U10=0 to 1 for U10=WSPMIN_RESET_TAUW
              IF (U10NEW(IJ) < WSPMIN_RESET_TAUW) THEN
                TLWMAX=WGHT*(ACD+BCD*U10NEW(IJ))*U10NEW(IJ)**3
                TAUW(IJ)=MIN(TAUW(IJ),TLWMAX)
              ENDIF
            ENDDO
          ELSE
            DO IJ = KIJS, KIJL
              USNEW(IJ)=FF_NEXT(IJ,ISTORE)%USTAR
! update the estimate of TAUW
              TAUW(IJ)=USNEW(IJ)**2*(1.0_JWRB-(ALPHA/BETAOLD(IJ))**2)
! adapt first estimate of wave induced stress for low winds
              IF (USNEW(IJ) < USTMIN_RESET_TAUW) TAUW(IJ)=0.0_JWRB
            ENDDO
          ENDIF
          DO IJ = KIJS, KIJL
            THWNEW(IJ)=FF_NEXT(IJ,ISTORE)%WDWAVE
            ROAIRN(IJ)=FF_NEXT(IJ,ISTORE)%AIRD
            ZIDLNEW(IJ)=FF_NEXT(IJ,ISTORE)%ZIDL
            CICOVER(IJ)=FF_NEXT(IJ,ISTORE)%CIFR
            CITHICK(IJ)=FF_NEXT(IJ,ISTORE)%CITH
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1492,1)

        ISTORE=ISTORE+1

        CALL INCDATE(CDATEWH, IDELWO)
        NEWREAD = .TRUE.   

!       UPDATE THE SEA ICE REDUCTION FACTOR
        CALL CIREDUCE (IJS, IJL, CGROUP, CICOVER, CITHICK, CIWA)

      ENDIF

      IF (LHOOK) CALL DR_HOOK('NEWWIND',1,ZHOOK_HANDLE)

      END SUBROUTINE NEWWIND
