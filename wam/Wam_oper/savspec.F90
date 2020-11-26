      SUBROUTINE SAVSPEC(FL, NBLKS, NBLKE, CDTPRO, CDATEF, CDATER)

! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      MARCH 1997
!     J. BIDLOT    ECMWF      SEPTEMBER 1997 : use PBIO SOFTWARE
!     B. HANSEN    ECMWF      FEBRUARY  1998 : use FDB SOFTWARE

!*    PURPOSE.
!     --------
!     WRITE BINARY SPECTRA TO DISK.

!**   INTERFACE.
!     ----------
!     *CALL* *SAVSPEC(FL,NBLKS,NBLKE,CDTPRO,CDATEF,CDATER)
!     *FL*        ARRAY CONTAINING THE SPECTRA CONTRIBUTION ON EACH PE
!     *NBLKS*     INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!     *NBLKE*     INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!     *CDTPRO*    CHAR*14   END DATE OF PROPAGATION.
!     *CDATEF*    CHAR*14   END DATE OF ANALYSIS RUN (YYMMDDHHMM).
!     *CDATER*    CHAR*14   DATE FOR OUTPUT OF BOTH RESTART FILES

!     METHOD.
!     -------
!     IN CASE OF MESSAGE PASSING, THE OUTPUT IS DIRECTED STRAIGHT TO ITS
!     PLACE ON DISK, THE CORRESPONDING NAME AND DIRECTORY IS DETERMINED
!     BY GRSTNAME OR BY THE FDB SOFTWARE IF FDB IS USED.
!     THE CONTRIBUTION FROM ALL PE'S NEED TO BE GATHERED ON THE OUTPUT
!     PE BY CALLING MPGATHERFL.
!     IN ORDER TO USE MPL_GATHERV IN MPGATHERFL THE SPECTRA WILL BE
!     GATHERED MDEL FREQUENCIES AT A TIME
!     AND KDEL DIRECTIONS AT A TIME.

!     EXTERNALS.
!     ----------
!     GETENV
!     GRSTNAME
!     MPGATHERFL
!     WRITEFL

!     REFERENCE.
!     ----------
!     NONE
! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : JPPFLAG  ,IPFGTBL  ,KDEL    ,MDEL      ,    &
     &           LRSTPARALW
      USE YOWGRID  , ONLY : IJSLOC   ,IJLLOC   ,IJGLOBAL_OFFSET
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP     ,    &
     &            NPRECR
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NIBLO
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH
      USE YOWUNIT  , ONLY : IU12     ,IU14     ,IU15

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "expand_string.intfb.h"
#include "grstname.intfb.h"
#include "mpgatherfl.intfb.h"
#include "writefl.intfb.h"

      INTEGER(KIND=JWIM), DIMENSION(NPROC), INTENT(IN) :: NBLKS, NBLKE
      REAL(KIND=JWRB), DIMENSION(NINF-1:NSUP,NANG,NFRE), INTENT(INOUT) :: FL
      CHARACTER(LEN=14), INTENT(IN) :: CDTPRO, CDATEF, CDATER

      INTEGER(KIND=JWIM) :: IJ, K, M
      INTEGER(KIND=JWIM) :: IUNIT
      INTEGER(KIND=JWIM) :: IRECV, MLOOP, KLOOP, MINF, MSUP, KINF, KSUP 
      INTEGER(KIND=JWIM) :: LNAME

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:,:) :: RFL

      CHARACTER(LEN=296) :: FILENAME

      LOGICAL :: LOUNIT

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SAVSPEC',0,ZHOOK_HANDLE)

      LOUNIT = .TRUE.

      CALL GRSTNAME(CDTPRO,CDATEF,'BLS',ICPLEN,CPATH,FILENAME)

      IF(LRSTPARALW) THEN
!        RESTART FILES FROM ALL PS's
         LNAME = LEN_TRIM(FILENAME)
         FILENAME=FILENAME(1:LNAME)//'.%p_%n'
         CALL EXPAND_STRING(IRANK,NPROC,0,0,FILENAME,1)

         CALL WRITEFL(FL, NINF-1, NSUP, 1, NANG, 1, NFRE,               &
     &                FILENAME, IUNIT, LOUNIT, LRSTPARALW)

      ELSE

!       COLLECT THE DIFFERENT CONTRIBUTIONS AND PROCEED TO THE OUTPUT OF
!       SPECTRUM FROM PE IPFGTBL(JPPFLAG+1)
        IRECV=IPFGTBL(JPPFLAG+1)

        DO MLOOP=1,NFRE,MDEL
          MINF=MLOOP
          MSUP=MIN(MLOOP+MDEL-1,NFRE)
          DO KLOOP=1,NANG,KDEL
            KINF=KLOOP
            KSUP=MIN(KLOOP+KDEL-1,NANG)

            LOUNIT = .FALSE.
            IF(MINF.EQ.1 .AND. KINF.EQ.1) LOUNIT = .TRUE.

            ALLOCATE(RFL(0:NIBLO,KINF:KSUP,MINF:MSUP))
            RFL(0,:,:)=0.0_JWRB
            DO M=MINF,MSUP
              DO K=KINF,KSUP
                DO IJ = IJSLOC, IJLLOC
                  RFL(IJ+IJGLOBAL_OFFSET,K,M) = FL(IJ,K,M)
                ENDDO
              ENDDO
            ENDDO

            CALL MPGATHERFL(IRECV,NBLKS,NBLKE,KINF,KSUP,MINF,MSUP,RFL)

            IF (ITEST.GE.2)                                             &
     &       WRITE(IU06,*)                                              &
     &       'SUB. SAVSPEC: RESTART SPECTRUM COLLECTED, :',MLOOP,KLOOP

            IF (IRANK.EQ.IPFGTBL(JPPFLAG+1) .OR. NPROC.EQ.1) THEN
              CALL WRITEFL(RFL, 0, NIBLO, KINF, KSUP, MINF, MSUP,       &
     &                     FILENAME, IUNIT, LOUNIT, LRSTPARALW)
            ENDIF
            DEALLOCATE(RFL)
          ENDDO
        ENDDO
      ENDIF

      WRITE(IU06,*) ' SPECTRUM FILE DISPOSED AT........',               &
     &                ' CDTPRO  = ', CDTPRO

      IF (LHOOK) CALL DR_HOOK('SAVSPEC',1,ZHOOK_HANDLE)

      END SUBROUTINE SAVSPEC
