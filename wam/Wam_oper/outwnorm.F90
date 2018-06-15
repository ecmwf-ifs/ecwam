      SUBROUTINE OUTWNORM(LDREPROD)

! ----------------------------------------------------------------------

!**** *OUTWNORM* - PRINTS AVERAGE, MINIMUM AND MAXIMUM VALUES OF
!                  INTEGRATED PARAMETER FIELDS 

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
     &                      COUTNAME, ITOBOUT ,NIPRMOUT, BOUT
      USE YOWGRID  , ONLY : IJSLOC   ,IJLLOC
      USE YOWMPP   , ONLY : IRANK   ,NPROC
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWTEST  , ONLY : IU06
      USE MPL_MODULE
      USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "mpminmaxavg.intfb.h"

      LOGICAL, INTENT(IN) :: LDREPROD

      INTEGER(KIND=JWIM) :: ITG, IT, I, IRECV

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(4,NIPRMOUT) :: WNORM

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('OUTWNORM',0,ZHOOK_HANDLE)

      IF(NFLAGALL .AND. NIPRMOUT > 0) THEN

       CALL MPMINMAXAVG(IJSLOC, IJLLOC, NIPRMOUT, BOUT, ZMISS, LDREPROD, WNORM)

        IRECV=1
!       WRITE NORM TO LOGFILE
          WRITE(IU06,*) ' ' 
          IF(LDREPROD) THEN
            WRITE(IU06,*) '  WAMNORM ON ',CDTPRO
            WRITE(IU06,*) '  !!!!!!!!! REPRODUCEABLE ONLY IF SAME NPROC'
          ELSE
            WRITE(IU06,*) '  !!!!!!!!!!!! NON reproduceable WAMNORM ON ',CDTPRO
          ENDIF
          WRITE(IU06,*) '          AVERAGE               MINIMUM  ',    &
     &     '           MAXIMUM    NON MISSING POINTS  NPROC'
          DO ITG=1,JPPFLAG
            IF(NFLAG(ITG)) THEN
              IT = ITOBOUT(ITG)
              WRITE(IU06,*) '  WAMNORM FOR ',COUTNAME(ITG)
              WRITE(IU06,*) '  ',(WNORM(I,IT),I=1,3),INT(WNORM(I,IT)),NPROC
              WRITE(IU06,111) (WNORM(I,IT),I=1,3),INT(WNORM(I,IT)),NPROC
            ENDIF
          ENDDO
111     FORMAT(4x,'HEX: ',3(Z16.16,2x),1x,i8,1x,i6)
        CALL FLUSH(IU06)

      ENDIF

      IF (LHOOK) CALL DR_HOOK('OUTWNORM',1,ZHOOK_HANDLE)

      END SUBROUTINE OUTWNORM
