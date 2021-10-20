      SUBROUTINE OUTINT(CDATE, IFCST)

! ----------------------------------------------------------------------

!**** *OUTINT* - OUTPUT OF INTEGRATED FIELDS WITHOUT THE I/O SERVER

!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!       *CALL* *OUTINT(CDATE, IFCST)
!          *CDATE*   DATE AND TIME.
!          *IFCST*   FORECAST STEP IN HOURS.

!     EXTERNALS.
!     ----------
!       *OUTGRID*
!       *WGRIBENOUT*

!     METHOD.
!     -------

!       NONE.

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : FFLAG    ,FFLAG20  ,GFLAG    ,GFLAG20  ,    &
     &            JPPFLAG  ,LFDB     ,IPFGTBL  ,ITOBOUT  ,INFOBOUT ,    &
     &            LOUTINT
      USE YOWGRID  , ONLY : NLONRGG  ,DELPHI
      USE YOWINTP  , ONLY : GOUT
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,    &
     &            AMONOP   ,ZDELLO
      USE YOWMPP   , ONLY : IRANK
      USE YOWPARAM , ONLY : NGX      ,NGY      ,CLDOMAIN
      USE YOWSTAT  , ONLY : CDATEF   ,CDTPRO   ,CDTINTT  ,MARSTYPE
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH
      USE YOWUNIT  , ONLY : IU20

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE GRIB_API_INTERFACE

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "abort1.intfb.h"
#include "grstname.intfb.h"
#include "outgrid.intfb.h"
#include "out_onegrdpt.intfb.h"
#include "wgribenout.intfb.h"

      CHARACTER(LEN=14), INTENT(IN) :: CDATE
      INTEGER(KIND=JWIM), INTENT(IN) :: IFCST

      INTEGER(KIND=JWIM) :: I, J, ITG, IFLAG, IT
      INTEGER(KIND=JWIM) :: ICOUNT      ! Field counter
      INTEGER(KIND=JWIM) :: IPARAM, ITABLE, IZLEV
      INTEGER(KIND=JWIM) :: IERR
      INTEGER(KIND=JWIM) :: LFILE, IUOUT 
  
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      CHARACTER(LEN=3) :: FILEID
      CHARACTER(LEN=296) :: OUTFILEN

      LOGICAL :: LLIUOUTOPEN

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('OUTINT',0,ZHOOK_HANDLE)

      LOUTINT=.TRUE.

!     1. COLLECT INTEGRATED PARAMETERS FOR OUTPUT ON SELECTED PE's 
!        i.e BOUT => GOUT
!        ---------------------------------------------------------
      CALL OUTGRID


!     2. ONE GRID POINT OUTPUT
!        ---------------------
      IF(CLDOMAIN.EQ.'s') CALL OUT_ONEGRDPT(IU06)


!*    3. OUTPUT IN PURE BINARY FORM (obsolete)
!        --------------------------
      IF(FFLAG20 .AND. CDTINTT.EQ.CDTPRO ) THEN
        DO IFLAG=1,JPPFLAG
          IF (FFLAG(IFLAG) .AND. IRANK.EQ.IPFGTBL(IFLAG)) THEN
            WRITE (IU20) CDTPRO, NGX, NGY, IRGG
            WRITE (IU20) AMOWEP,AMOSOP,AMOEAP,AMONOP
            WRITE (IU20) ZDELLO, NLONRGG, DELPHI
            WRITE (IU20) GOUT(IFLAG,:,:)
          ENDIF
        ENDDO
      ENDIF

!*    4.  PACK INTEGRATED PARAMETERS INTO GRIB AND OUTPUT
!         -----------------------------------------------
      IF (GFLAG20 .AND. CDTINTT.EQ.CDTPRO ) THEN

!       DEFINE OUTPUT FILE FOR GRIB DATA (if not written to FDB)
        IF(.NOT.LFDB .AND.(IRANK.EQ.1)) THEN
       !  output to file should only take place on one PE
         IF(MARSTYPE == 'fg') THEN
           ! fg only exist in uncoupled anlysis experiments
           FILEID = 'MFG'
         ELSE
           FILEID = 'MPP'
         ENDIF
         CALL GRSTNAME(CDTPRO,CDATEF,IFCST,FILEID,ICPLEN,CPATH,OUTFILEN)
         LFILE=LEN_TRIM(OUTFILEN)
         CALL IGRIB_OPEN_FILE(IUOUT,OUTFILEN(1:LFILE),'w')
         WRITE(IU06,*) '  INTEGRATED PARAMETERS WRITTEN TO FILE ',OUTFILEN(1:LFILE)
         LLIUOUTOPEN=.TRUE.
        ELSE
          LLIUOUTOPEN=.FALSE.
          IUOUT=0
        ENDIF

!       OUTPUT:
        ICOUNT=0
        DO IFLAG=1,JPPFLAG
          IF (GFLAG(IFLAG)) THEN
            IF (IRANK.EQ.IPFGTBL(IFLAG)) THEN
              ICOUNT=ICOUNT+1
              IF(ICOUNT.GT. SIZE(GOUT,1)) THEN
                WRITE(*,*) ' -------------------------------------'
                WRITE(*,*) ' ERROR in OUTINT '
                WRITE(*,*) ' ACCESSING MORE FIELDS THAN AVAILABLE'
                WRITE(*,*) ' SIZE(GOUT,1) = ',SIZE(GOUT,1) 
                WRITE(*,*) ' -------------------------------------'
                CALL ABORT1
              ENDIF

              IT=ITOBOUT(IFLAG)
              ITABLE=INFOBOUT(IT,1)
              IPARAM=INFOBOUT(IT,2)
              IZLEV=INFOBOUT(IT,3)

              CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, GOUT(ICOUNT,:,:),  &
     &                        ITABLE, IPARAM, IZLEV, 0 , 0,             &
     &                        CDATE, IFCST, MARSTYPE, LFDB, IUOUT)
            ENDIF
          ENDIF
        ENDDO

        IF(LLIUOUTOPEN) THEN
          CALL IGRIB_CLOSE_FILE(IUOUT)
        ENDIF

      ENDIF ! end grib output

!     CLEAN-UP
      IF(ALLOCATED(GOUT)) DEALLOCATE(GOUT)


! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('OUTINT',1,ZHOOK_HANDLE)

      END SUBROUTINE OUTINT
