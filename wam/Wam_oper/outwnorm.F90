      SUBROUTINE OUTWNORM(IJS, IJL, BOUT, LDREPROD)

! ----------------------------------------------------------------------

!**** *WAMNORM* - PRINTS AVERAGE, MINIMUM AND MAXIMUM VALUES OF
!                 INTEGRATED PARAMETER FIELDS 

!*    PURPOSE.
!     --------


!**   INTERFACE.
!     ----------

!     LDREPROD -  Reproducibility flag for SUMmation-operator.
!                 .TRUE. ==> Use home-written binary tree, No MPI_ALLREDUCE used.
!                 .FALSE. ==> let MPI_ALLREDUCE do the summation.

!     METHOD.
!     -------

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : NFLAG, NFLAGALL, JPPFLAG,                   &
     &                      COUTNAME, NIPRMOUT, ITOBOUT
      USE YOWMPP   , ONLY : IRANK
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWTEST  , ONLY : IU06
      USE MPL_MODULE
      USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NIPRMOUT), INTENT(IN) :: BOUT
      LOGICAL, INTENT(IN) :: LDREPROD

      INTEGER(KIND=JWIM) :: ITG, IT, I, IRECV

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(4,NIPRMOUT) :: WNORM

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTWNORM',0,ZHOOK_HANDLE)
#endif
      IF(NFLAGALL .AND. NIPRMOUT > 0) THEN

       CALL MPMINMAXAVG(IJS, IJL, NIPRMOUT, BOUT, ZMISS, LDREPROD, WNORM)

        IRECV=1
!       WRITE NORM TO LOGFILE
        IF(IRANK.EQ.IRECV) THEN
          WRITE(IU06,*) ' ' 
          WRITE(IU06,*) '   WAMNORM ON ',CDTPRO
          WRITE(IU06,*) '          AVERAGE            MINIMUM  ',       &
     &     '           MAXIMUM        SEA POINTS' 
          DO ITG=1,JPPFLAG
            IF(NFLAG(ITG)) THEN
              IT = ITOBOUT(ITG)
              WRITE(IU06,*) COUTNAME(ITG)
              WRITE(IU06,*) '  ',(WNORM(I,IT),I=1,4)
              WRITE(IU06,111) (WNORM(I,IT),I=1,4)
            ENDIF
          ENDDO
        ELSE
          WRITE(IU06,*) ' ' 
          WRITE(IU06,*) '   WAMNORM ON ',CDTPRO,' SEE OUTPUT PE',IRECV
        ENDIF
111     FORMAT(4x,'HEX: ',4(Z16.16,2x))
        CALL FLUSH(IU06)

      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTWNORM',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE OUTWNORM
