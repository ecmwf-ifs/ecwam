      SUBROUTINE MAKEGRID (BLOCK, GRID, IG, LSQRT, PMISS)

! ----------------------------------------------------------------------

!**** *MAKEGRID* - MAKE GRIDDED WAM MODEL FIELDS FROM BLOCKED FIELDS.

!     H. GUNTHER       ECMWF    NOVEMBER 1989

!*    PURPOSE.
!     --------

!       GRIDDED WAVE FIELDS ARE CREATED FROM BLOCKED WAVE FIELDS.

!**   INTERFACE.
!     ----------

!        *CALL MAKEGRID*
!           *BLOCK*   REAL   DATA IN BLOCKED FORMAT
!           *GRID*    REAL   DATA IN GRID FORMAT
!           *IG*      INTEGER BLOCK NUMBER
!           *LSQRT*   LOGICAL IF TRUE THE 4*SQUARE ROOT OF BLOCK IS
!                     COMPUTED WHEN TRANSFERRING TO GRID EXCEPT IF
!                     BLOCK EQUALS PMISS.
!           *PMISS*   REAL  VALUE GIVEN FOR ALL NON SEA POINTS. 

!     METHOD.
!     -------

!       THE PARAMETER, WHICH IS GIVEN IN BLOCKED FORMAT, IS
!       DISTRIBUTED IN GRID FORMAT. ONE BLOCK IS DONE IN ONE
!       CALL. BEFORE THE FIRST BLOCK IS TRANSFORMED THE GRID ARRAY IS
!       INITIALISED WITH PMISS.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,AMOWEP   ,AMONOP   ,
     &            XDELLA   ,ZDELLO   ,IPER
      USE YOWMPP   , ONLY : NPROC
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NIBLO
      USE YOWUNPOOL, ONLY : LLUNSTR  ,IE_OUTPTS
      USE YOWUNPOOL, ONLY : NIBLO_OUT, OUT_METHOD
      USE OUTPUT_STRUCT, ONLY : IXarr, IYarr
      USE YOWPD,     ONLY : NODES=>nodes_global,INE_GLOBAL,rank, np
      USE OUTPUT_STRUCT, ONLY : INTELEMENT_IPOL
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "unblkrord.intfb.h"

      INTEGER :: IG
      INTEGER :: IY, IX, J, I, IJ, IE, KI, IR, IP
      INTEGER :: NI(3)

      REAL    :: PMISS
      REAL    :: ZHOOK_HANDLE
      REAL    :: GRID(NGX,NGY)
      REAL, DIMENSION(NIBLO_OUT) :: BLOCK
      REAL, allocatable :: BLOCK_G(:)
      REAL*8  :: XLO, XLA
      REAL*8  :: BLOCK8, GRID8, PMS8
      REAL*8  :: XYELE(2,3), SKALAR(3)

      LOGICAL :: LSQRT
      LOGICAL :: LLCLST
      integer iNode
! ----------------------------------------------------------------------
      
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('MAKEGRID',0,ZHOOK_HANDLE)
#endif

      CALL GSTATS(1996,0)

!*    1. INITIALIZATION
!     -----------------

      PMS8=PMISS

      DO J = 1,NGY
        DO I = 1,NGX
          GRID(I,J) = PMISS
        ENDDO
      ENDDO

      IF(LSQRT) THEN
        DO IJ = 1,NIBLO_OUT
          IF(BLOCK(IJ).NE.PMISS) THEN
            BLOCK(IJ) = 4.*SQRT(BLOCK(IJ))
          ELSE
            BLOCK(IJ) = PMISS
          ENDIF
        ENDDO
      ENDIF

!*    2. MAKE GRIDDED FIELD
!     ---------------------

      IF (LLUNSTR) THEN
        IF (OUT_METHOD .eq. 1) THEN
          allocate(BLOCK_G(NIBLO_OUT))
          CALL UNBLKRORD(1,1,NIBLO,1,1,1,1,BLOCK(1),BLOCK_G(1))
!!! I believe one could have an openmp loop here !!!!
          DO IY=1,NGY
            XLA = AMONOP-REAL(IY-1)*XDELLA
            DO IX=1,NGX
              XLO = AMOWEP+REAL(IX-1)*ZDELLO(IY)
!!debile
!!! need to find out unstructured grid left longitude (I assume -180 for now)
              IF(XLO.GE.180.d0) XLO=max(XLO-IPER*360.d0,-180.d0)

              IE=IE_OUTPTS(IX,IY)
              IF(IE.NE.-1) THEN
!!! this is not very efficient
                LLCLST=.FALSE.
                NI = INE_GLOBAL(:,IE)
                XYELE(1,:)=NODES(NI)%X
                XYELE(2,:)=NODES(NI)%Y
                SKALAR(:)=BLOCK_G(NI)
                IF(ANY(BLOCK_G(NI).EQ.PMISS)) LLCLST=.TRUE.
                CALL INTELEMENT_IPOL(LLCLST,XYELE,SKALAR,XLO,XLA,IE,
     &                               PMS8,GRID8)
                GRID(IX,IY)=GRID8
              ENDIF
            ENDDO
          ENDDO
          deallocate(BLOCK_G)
        ENDIF
        IF (OUT_METHOD .eq. 2) THEN
          DO iNode=1,NIBLO_OUT
            IX=IXarr(iNode)
            IY=IYarr(iNode)
            GRID(IX,IY) = BLOCK(iNode)
          END DO
        END IF
      ELSE
        DO IJ = 1, NIBLO 
          IX = IXLG(IJ,IG)
          IY = NGY- KXLT(IJ,IG) +1
          GRID(IX,IY) = BLOCK(IJ)
        ENDDO
      ENDIF ! LLUNSTR

      CALL GSTATS(1996,1)

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('MAKEGRID',1,ZHOOK_HANDLE)
#endif

      RETURN
      END SUBROUTINE MAKEGRID
