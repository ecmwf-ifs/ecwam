      SUBROUTINE PROPDOT

! ----------------------------------------------------------------------

!**** *PROPDOT* - PROPAGATION DOT TERMS FROM DEPTH AND CURRENT GRADIENT.

!     H. GUNTHER   GKSS/ECMWF   17/02/91
!     J. BIDLOT    ECMWF        FEBRUARY 97 MESSAGE PASSING
!     J. BIDLOT    ECMWF        JANUARY 2004 REMOVE MULTI-BLOCK OPTION.

!*    PURPOSE.
!     --------

!       COMPUTATION OF COMMON REFDOT FOR PROPAGATION.

!**   INTERFACE.
!     ----------

!       *CALL* *PROPDOT*

!     METHOD.
!     -------

!       IN A LOOP OVER THE BLOCKS THE COMMON UBUF IS READ,
!       THE DEPTH AND CURRENT GRADIENTS ARE COMPUTED,
!       COMMON REFDOT (DEPTH AND CURRENT REFRACTION FOR THETA DOT)
!       IS COMPUTED AND WRITTEN TO MASS STORAGE (IU16) (if multiblock).
!       IN CASE OF CURRENT REFRACTION THE COMPLETE SIGMA DOT TERM
!       IS COMPUTED AND WRITTEN TO IU16 ADDITIONALLY (if multiblock).
!       WRITE OPERATIONS ARE NOT DONE FOR COMMON UBUF AND REFDOT
!       IF THIS IS A ONE BLOCK MODEL.

!     EXTERNALS.
!     ----------

!       *GRADI*     - COMPUTES DEPTH AND CURRENT GRADIENTS.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCURR  , ONLY : U        ,V
      USE YOWFRED  , ONLY : COSTH    ,SINTH
      USE YOWGRID  , ONLY : COSPHM1   ,IGL      ,IJS      ,IJL
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWREFD  , ONLY : THDD     ,THDC     ,SDOT
      USE YOWSHAL  , ONLY : TCGOND   ,TFAK     ,TSIHKD   ,INDEP 
      USE YOWSTAT  , ONLY : ICASE    ,ISHALLO  ,IREFRA
      USE YOWSPEC  , ONLY : NSTART   ,NEND
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "gradi.intfb.h"

      INTEGER :: IG
      INTEGER :: IJ, K, M

      REAL :: CD, SD, SS, SC, CC

      REAL, ALLOCATABLE, DIMENSION(:) :: DDPHI, DDLAM, DUPHI, DULAM
      REAL, ALLOCATABLE, DIMENSION(:) :: DVPHI, DVLAM, DCO, OMDD

! ----------------------------------------------------------------------

!     ARRAY TO KEEP DEPTH AND CURRENT REFRACTION FOR THETA DOT
!     AND SIGMA DOT
      IF(.NOT.ALLOCATED(THDC))
     &     ALLOCATE(THDC(NSTART(IRANK):NEND(IRANK),NANG))
      IF(.NOT.ALLOCATED(THDD))
     &    ALLOCATE(THDD(NSTART(IRANK):NEND(IRANK),NANG))
      IF(.NOT.ALLOCATED(SDOT))
     &    ALLOCATE(SDOT(NSTART(IRANK):NEND(IRANK),NANG,NFRE))



!*    2. LOOP OVER BLOCKS.
!        -----------------

      IF (IGL.NE.1) THEN
        write(*,*) 'PROPDOT DOES NOT WORK FOR THIS CONFIGURATION !!!'
        WRITE(IU06,*) 'PROPDOT DOES NOT WORK FOR THIS CONFIGURATION !!!'
        CALL ABORT1
      ENDIF

      DO IG = 1,IGL

!       WORK ARRAY
        ALLOCATE(DDPHI(IJS(IG):IJL(IG)))
        ALLOCATE(DDLAM(IJS(IG):IJL(IG)))
        ALLOCATE(DUPHI(IJS(IG):IJL(IG)))
        ALLOCATE(DULAM(IJS(IG):IJL(IG)))
        ALLOCATE(DVPHI(IJS(IG):IJL(IG)))
        ALLOCATE(DVLAM(IJS(IG):IJL(IG)))

        ALLOCATE(DCO(IJS(IG):IJL(IG)))
        ALLOCATE(OMDD(IJS(IG):IJL(IG)))


!*    2.2 DEPTH AND CURRENT GRADIENTS.
!         ----------------------------

        CALL GRADI (IG, IJS(IG),IJL(IG),IREFRA, DDPHI, DDLAM, DUPHI,
     &              DULAM, DVPHI, DVLAM)

!*    2.3 COSINE OF LATITUDES IF SPHERICAL PROPAGATION.
!         ---------------------------------------------

        IF (ICASE.EQ.1) THEN
          DO IJ = IJS(IG),IJL(IG)
            DCO(IJ) = COSPHM1(IJ,IG)
          ENDDO
        ELSE
          DO IJ = IJS(IG),IJL(IG)
            DCO(IJ) = 1.
          ENDDO
        ENDIF

!*    2.4 DEPTH GRADIENT PART OF SIGMA DOT.
!         ---------------------------------

        IF (ISHALLO.NE.1 ) THEN
          IF (IREFRA.EQ.3) THEN
            DO IJ = IJS(IG),IJL(IG)
              OMDD(IJ) = V(IJ,IG)*DDPHI(IJ) + U(IJ,IG)*DDLAM(IJ)*DCO(IJ)
            ENDDO
          ELSEIF (IREFRA.EQ.2) THEN
            DO IJ = IJS(IG),IJL(IG)
              OMDD(IJ) = 0.0
            ENDDO
          ENDIF
        ENDIF

!*    2.5. LOOP OVER DIRECTIONS.
!          ---------------------

        DO K=1,NANG
          SD = SINTH(K)
          CD = COSTH(K)

!*    2.5.1. DEPTH GRADIENT OF THETA DOT.
!            ----------------------------

          IF (ISHALLO.NE.1) THEN
            IF (IREFRA.EQ.1 .OR. IREFRA.EQ.3) THEN
              DO IJ = IJS(IG),IJL(IG)
                THDD(IJ,K) = SD*DDPHI(IJ) - CD*DDLAM(IJ)*DCO(IJ)
              ENDDO
            ELSE
              DO IJ = IJS(IG),IJL(IG)
                THDD(IJ,K) = 0.0
              ENDDO
            ENDIF
          ENDIF

!*    2.5.2 SIGMA DOT AND THETA DOT PART FROM CURRENT GRADIENT.
!           ---------------------------------------------------

          IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN
            SS  = SD**2
            SC  = SD*CD
            CC  = CD**2
            DO IJ = IJS(IG),IJL(IG)
              SDOT(IJ,K,NFRE) = -SC*DUPHI(IJ) - CC*DVPHI(IJ)
     &                        - (SS*DULAM(IJ) + SC*DVLAM(IJ))*DCO(IJ)
              THDC(IJ,K) =  SS*DUPHI(IJ) + SC*DVPHI(IJ)
     &                    - (SC*DULAM(IJ) + CC*DVLAM(IJ))*DCO(IJ)
            ENDDO

!*    2.5.3 LOOP OVER FREQUENCIES.
!           ----------------------

            IF (ISHALLO.NE.1) THEN
              DO M=1,NFRE
                DO IJ=IJS(IG),IJL(IG)
                  SDOT(IJ,K,M) = (SDOT(IJ,K,NFRE)*TCGOND(INDEP(IJ),M)
     &             + OMDD(IJ)*TSIHKD(INDEP(IJ),M))
     &             * TFAK(INDEP(IJ),M)
                ENDDO

!*    BRANCH BACK TO 2.5.3 FOR NEXT FREQUENCY.

              ENDDO
            ENDIF
          ENDIF

!*    BRANCH BACK TO 2.5 FOR NEXT DIRECTION.

        ENDDO

        DEALLOCATE(DDPHI, DDLAM, DUPHI, DULAM)
        DEALLOCATE(DVPHI, DVLAM, DCO, OMDD)

!*    BRANCH BACK TO 2. FOR NEXT BLOCK.
      ENDDO


      RETURN
      END SUBROUTINE PROPDOT
