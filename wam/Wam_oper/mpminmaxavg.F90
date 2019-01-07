      SUBROUTINE MPMINMAXAVG(LLGLOBAL, IRECV, LDREPROD, WNORM)
! ----------------------------------------------------------------------

!****  *MPMINMAXAVG* - FIND GLOBAL MIN, MAX AND AVERAGE OF BOUT

!           LLGLOBAL -  IF TRUE NORMS FOR SELECTED OUPTUT FIELDS WILL BE GLOBAL
!                       ON PROCESSOR IRECV
!           LDREPROD -  only if LLGLOBAL is false
!                       Reproducibility flag for SUMmation-operator.
!                       Meaningful only for REAL-numbers.
!                       Three modes (applicable for REAL-number only):
!                       1) Not provided at all (the default) ==> MPL_ABORT
!                       2) Provided and .TRUE. ==> Use home-written binary tree
!                          No MPI_ALLREDUCE used.
!                       3) Provided, but .FALSE. ==> let MPI_ALLREDUCE do the summation.

!           *WNORM(1,:) : MINIMUM
!           *WNORM(2,:) : MAXIMUM
!           *WNORM(3,:) : AVERAGE
!           *WNORM(4,:) : NUMBER OF VALUES USED TO PRODUCE THE AVERAGE


! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT   , ONLY : NFLAG    ,JPPFLAG , ITOBOUT ,NIPRMOUT, BOUT
      USE YOWGRID   , ONLY : IJSLOC   ,IJLLOC
      USE YOWMPP    , ONLY : IRANK    ,NPROC
      USE YOWPARAM  , ONLY : NIBLO
      USE YOWSPEC   , ONLY : NBLKS    ,NBLKE
      USE YOWPCONS  , ONLY : ZMISS
      USE MPL_MODULE, ONLY : MPL_ALLREDUCE
      USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "abort1.intfb.h"
#include "mpgatherscfld.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IRECV
      REAL(KIND=JWRB),DIMENSION(4,NIPRMOUT), INTENT(OUT) :: WNORM
      LOGICAL, INTENT(IN) :: LLGLOBAL
      LOGICAL, INTENT(IN) :: LDREPROD

      INTEGER(KIND=JWIM) :: IJ, I, IP, ITG, IT

      REAL(KIND=JWRB),DIMENSION(NIPRMOUT) :: ZMIN, ZMAX
      REAL(KIND=JWRB),DIMENSION(2*NIPRMOUT) :: ZSUM
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), ALLOCATABLE :: ZGLOBAL(:)
!     -------------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('MPMINMAXAVG',0,ZHOOK_HANDLE)

      ZSUM(:)=0.0_JWRB
      ZMIN(:)=HUGE(ZMIN(:))
      ZMAX(:)=-HUGE(ZMAX(:))

      IF(LLGLOBAL) THEN

        WNORM(:,:)=0.0_JWRB

        IF(IJSLOC.NE.NBLKS(IRANK) .OR. IJLLOC.NE.NBLKE(IRANK) ) THEN
          WRITE(IU06,*) '************************************'
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '*  FATAL ERROR IN SUB. MPMINMAXAVG *'
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '* IJSLOC.NE.NBLKS(IRANK) .OR.      *'
          WRITE(IU06,*) '* IJLLOC.NE.NBLKE(IRANK)           *'
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '* THIS SHOULD NOT HAPPEN  !        *'
          WRITE(IU06,*) '************************************'
          CALL ABORT1
        ENDIF

        ALLOCATE(ZGLOBAL(NIBLO))
        ZGLOBAL(:)=ZMISS

        DO ITG=1,JPPFLAG
          IF(NFLAG(ITG)) THEN
            IT = ITOBOUT(ITG)
            DO IJ=IJSLOC,IJLLOC
              ZGLOBAL(IJ)= BOUT(IJ,IT)
            ENDDO
            CALL MPGATHERSCFLD(IRECV, NBLKS, NBLKE, ZGLOBAL, NIBLO)

            IF(IRANK.EQ.IRECV) THEN
              DO IJ=1,NIBLO
                IF(ZGLOBAL(IJ) /= ZMISS) THEN
                  ZSUM(IT) = ZSUM(IT) + ZGLOBAL(IJ)
                  ZMIN(IT) = MIN(ZMIN(IT), ZGLOBAL(IJ))
                  ZMAX(IT) = MAX(ZMAX(IT), ZGLOBAL(IJ))
                ENDIF
              ENDDO
              WNORM(1,IT)=ZSUM(IT)/NIBLO
              WNORM(2,IT)=ZMIN(IT)
              WNORM(3,IT)=ZMAX(IT)
              WNORM(4,IT)=NIBLO
            ENDIF

          ENDIF
        ENDDO

        DEALLOCATE(ZGLOBAL)
      ELSE
        DO I=1,NIPRMOUT
          DO IJ=IJSLOC,IJLLOC
            IF(BOUT(IJ,I) /= ZMISS) THEN
              ZSUM(I) = ZSUM(I) + BOUT(IJ,I)
              ZSUM(NIPRMOUT+I) = ZSUM(NIPRMOUT+I) + 1.0_JWRB
              ZMIN(I) = MIN(ZMIN(I), BOUT(IJ,I))
              ZMAX(I) = MAX(ZMAX(I), BOUT(IJ,I))
            ENDIF
          ENDDO
        ENDDO

        IF (NPROC > 1) THEN
          CALL MPL_ALLREDUCE(ZSUM,'SUM',LDREPROD=LDREPROD,                &
     &                     CDSTRING='MPMINMAXAVG:')
          CALL MPL_ALLREDUCE(ZMIN,'MIN',CDSTRING='MPMINMAXAVG VALMIN:')
          CALL MPL_ALLREDUCE(ZMAX,'MAX',CDSTRING='MPMINMAXAVG VALMAX:')
        ENDIF

        DO I=1,NIPRMOUT
          WNORM(2,I)=ZMIN(I)
          WNORM(3,I)=ZMAX(I)
          WNORM(4,I)=ZSUM(NIPRMOUT+I)
          IF (WNORM(4,I)<1.0_JWRB) THEN
            WNORM(1,I)=-HUGE(WNORM(3,I))
          ELSE
            WNORM(1,I)=ZSUM(I)/WNORM(4,I)
          ENDIF
        ENDDO

      ENDIF
     

      IF (LHOOK) CALL DR_HOOK('MPMINMAXAVG',1,ZHOOK_HANDLE)

      END SUBROUTINE MPMINMAXAVG
