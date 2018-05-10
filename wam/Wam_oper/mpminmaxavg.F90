      SUBROUTINE MPMINMAXAVG(IJS, IJL, NDIM, BLOCK, ZMISS, LDREPROD, WNORM)
! ----------------------------------------------------------------------

!****  *MPMINMAXAVG* - FIND GLOBAL MIN, MAX AND AVERAGE OF BLOCK(:,I) 
!                      FOR EACH I=1,NDIM

!           ZMISS : VALUE TO EXCLUDE FROM THE SEARCH

!           LDREPROD -  Reproducibility flag for SUMmation-operator.
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

      USE YOWMPP    , ONLY : NPROC
      USE MPL_MODULE, ONLY : MPL_ALLREDUCE
      USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, NDIM
      REAL(KIND=JWRB) :: ZMISS 
      REAL(KIND=JWRB),DIMENSION(IJS:IJL,NDIM), INTENT(IN) :: BLOCK
      REAL(KIND=JWRB),DIMENSION(4,NDIM), INTENT(OUT) :: WNORM
      LOGICAL, INTENT(IN) :: LDREPROD

      INTEGER(KIND=JWIM) :: IJ, I

      REAL(KIND=JWRB),DIMENSION(2*NDIM) :: ZMIN, ZMAX
      REAL(KIND=JWRB),DIMENSION(2*NDIM) :: ZSUM
      REAL(KIND=JWRB) :: ZHOOK_HANDLE


!     -------------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('MPMINMAXAVG',0,ZHOOK_HANDLE)

      ZSUM(:)=0.0_JWRB
      ZMIN(:)=HUGE(ZMIN(:))
      ZMAX(:)=-HUGE(ZMAX(:))

      DO I=1,NDIM
        DO IJ=IJS,IJL
          IF(BLOCK(IJ,I) /= ZMISS) THEN
            ZSUM(I) = ZSUM(I) + BLOCK(IJ,I)
            ZSUM(2*I) = ZSUM(2*I) + 1.0_JWRB
            ZMIN(I) = MIN(ZMIN(I), BLOCK(IJ,I))
            ZMAX(I) = MAX(ZMAX(I), BLOCK(IJ,I))
          ENDIF
        ENDDO
      ENDDO

      IF (NPROC > 1) THEN
        CALL MPL_ALLREDUCE(ZSUM,'SUM',LDREPROD=LDREPROD,                &
     &                     CDSTRING='MPMINMAXAVG:')
        CALL MPL_ALLREDUCE(ZMIN,'MIN',CDSTRING='MPMINMAXAVG VALMIN:')
        CALL MPL_ALLREDUCE(ZMAX,'MAX',CDSTRING='MPMINMAXAVG VALMAX:')
      ENDIF

      DO I=1,NDIM
        WNORM(1,I)=ZMIN(I)
        WNORM(2,I)=ZMAX(I)
        WNORM(4,I)=ZSUM(2*I)
        IF (WNORM(4,I)<1.0_JWRB) THEN
          WNORM(3,I)=-HUGE(WNORM(3,I))
        ELSE
          WNORM(3,I)=ZSUM(I)/WNORM(4,I)
        ENDIF
      ENDDO

      IF (LHOOK) CALL DR_HOOK('MPMINMAXAVG',1,ZHOOK_HANDLE)

      END SUBROUTINE MPMINMAXAVG
