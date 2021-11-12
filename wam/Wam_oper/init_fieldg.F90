SUBROUTINE INIT_FIELDG(IJS, IJL, BLK2LOC,                      &
 &                     LLALLOC_ONLY, LLINIALL, LLOCAL)

! ----------------------------------------------------------------------

!**** *INIT_FIELDG* - INITIALISES DATA STRUCTURE FIELDG 

!*    PURPOSE.
!     --------

!     ALLOCATES AND INITIALISES DATA STRUCTURE FIELDG 

!**   INTERFACE.
!     ----------

!     *CALL* *INIT_FIELDG(LLALLOC_ONLY,LLINIALL,LLOCAL)*

!          LLALLOC_ONLY LOGICAL : IF TRUE THEN ONLY ALLOCATION OF FIELDG IS DONE.
!          LLINIALL  LOGICAL : IF TRUE ALL STRUCTURE IS INITIALISED
!                               OTHERWISE ONLY XLON and YLAT
!          LLOCAL    LOGICAL : IF TRUE ONLY THE LOCAL SEA POINTS ARE GIVEN
!                              VALID COORDINATES (XLON, YLAT) 
!                              PROVIDED THE LOCAL INDEX ARRAYS EXIST.

!     REFERENCE.
!     -----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDLOC

      USE YOWGRID  , ONLY : NPROMA_WAM, NBLOC
      USE YOWMAP   , ONLY : AMOWEP   ,AMOSOP   ,XDELLA   ,ZDELLO, NLONRGG
      USE YOWPCONS , ONLY : ZMISS    ,ROAIR    ,WSTAR0
      USE YOWWIND  , ONLY : FIELDG   ,NXFF     ,NYFF    ,WSPMIN
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
      USE YOWUNPOOL ,ONLY : LLUNSTR
      USE YOWPD    , ONLY : MNP=>npa , XP=>x, YP=>y

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      TYPE(WVGRIDLOC), DIMENSION(IJS:IJL), INTENT(IN) :: BLK2LOC
      LOGICAL, INTENT(IN) :: LLALLOC_ONLY, LLINIALL, LLOCAL


      INTEGER(KIND=JWIM) :: IJ, IX, JY, JSN
      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      LOGICAL :: LLOCAL_EXIST
! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('INIT_FIELDG',0,ZHOOK_HANDLE)

      IF (ALLOCATED(FIELDG)) DEALLOCATE(FIELDG)
      ALLOCATE(FIELDG(NXFF,NYFF))

      IF (.NOT.LLALLOC_ONLY) THEN

      CALL GSTATS(1501,0)
      IF (LLINIALL) THEN
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JY,IX)
        DO JY=1,NYFF
          DO IX=1,NXFF
            FIELDG(IX,JY)%XLON    = ZMISS 
            FIELDG(IX,JY)%YLAT    = ZMISS
            FIELDG(IX,JY)%UWND    = 0.0_JWRB
            FIELDG(IX,JY)%VWND    = 0.0_JWRB
            FIELDG(IX,JY)%WSWAVE  = WSPMIN 
            FIELDG(IX,JY)%WDWAVE  = 0.0_JWRB
            FIELDG(IX,JY)%UFRIC   = 0.0_JWRB
            FIELDG(IX,JY)%CICOVER = 0.0_JWRB
            FIELDG(IX,JY)%CITHICK = 0.0_JWRB
            FIELDG(IX,JY)%LKFR    = 0.0_JWRB
            FIELDG(IX,JY)%AIRD    = ROAIR
            FIELDG(IX,JY)%WSTAR   = WSTAR0
            FIELDG(IX,JY)%UCUR    = 0.0_JWRB
            FIELDG(IX,JY)%VCUR    = 0.0_JWRB
            FIELDG(IX,JY)%TAUW    = 0.0_JWRB
            FIELDG(IX,JY)%TAUWDIR = 0.0_JWRB
            FIELDG(IX,JY)%Z0M     = 0.0_JWRB
            FIELDG(IX,JY)%Z0B     = 0.0_JWRB
            FIELDG(IX,JY)%CHNK    = 0.0_JWRB
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO
      ELSE
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JY,IX)
        DO JY=1,NYFF
          DO IX=1,NXFF
            FIELDG(IX,JY)%XLON   = ZMISS 
            FIELDG(IX,JY)%YLAT   = ZMISS
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO
      ENDIF

!     COORDINATES OF THE POINTS IN FIELDG THAT ARE NEEDED
      LLOCAL_EXIST=LLOCAL

      IF (LLOCAL_EXIST) THEN
        NPROMA=NPROMA_WAM
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ,IX,JY,JSN)
        DO JKGLO=IJS,IJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL)
          DO IJ=KIJS,KIJL
            IX = BLK2LOC(IJ)%IFROMIJ
            JY = BLK2LOC(IJ)%JFROMIJ
            JSN= BLK2LOC(IJ)%KFROMIJ
            IF (LLUNSTR) THEN
              FIELDG(IX,JY)%XLON = XP(IJ)
              FIELDG(IX,JY)%YLAT = YP(IJ)
            ELSE
              FIELDG(IX,JY)%XLON = AMOWEP + (IX-1)*ZDELLO(JSN)
              FIELDG(IX,JY)%YLAT = AMOSOP + (JSN-1)*XDELLA
            ENDIF
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO

      ELSE
        IF (LLUNSTR) THEN
!!!!
          write(*,*) 'In INIT_FIELDG : not yet ready for unstructured grid '
          call abort1

        ELSE
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(JY,IX,JSN)
          DO JY=1,NYFF
            JSN=NYFF-JY+1
            DO IX=1,NLONRGG(JSN)
              FIELDG(IX,JY)%XLON = AMOWEP + (IX-1)*ZDELLO(JSN)
              FIELDG(IX,JY)%YLAT = AMOSOP + (JSN-1)*XDELLA
            ENDDO
          ENDDO
!$OMP     END PARALLEL DO
        ENDIF

      ENDIF
      CALL GSTATS(1501,1)

      ENDIF ! LLALLOC_ONLY

IF (LHOOK) CALL DR_HOOK('INIT_FIELDG',1,ZHOOK_HANDLE)

END SUBROUTINE INIT_FIELDG
