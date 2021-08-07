      SUBROUTINE GETSTRESS(U10OLD, THWOLD, USOLD, TAUW, TAUWDIR, Z0OLD, &
     &                     ROAIRO, ZIDLOLD, CICOVER, CITHICK,           &
     &                     NBLKS, NBLKE, IREAD)
! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      SEPTEMBER 1997 

!     S. ABDALLA   ECMWF      OCTOBER 2001  AIR DENSITY & Zi/L

!*    PURPOSE.
!     --------
!     READS RESTART WIND AND STRESS FIELDS FROM DISK OR USE GRIB WIND
!     AND DRAG COEFFICIENT TO RECONSTRUCT ALL WIND AND STRESS FIELDS.

!**   INTERFACE.
!     ----------
!     *CALL**GETSTRESS(U10OLD,THWOLD,USOLD,TAUW,TAUWDIR,Z0OLD,NBLKS,NBLKE,IREAD)
!     *U10OLD*    WIND SPEED.
!     *THWOLD*    WIND DIRECTION (RADIANS).
!     *USOLD*     FRICTION VELOCITY.
!     *TAUW*      WAVE STRESS.
!     *TAUWDIR*   WAVE STRESS DIRECTION.
!     *Z0OLD*     ROUGHNESS LENGTH IN M.
!     *ROAIRO*    AIR DENSITY IN KG/M3.
!     *ZIDLOLD*   Zi/L (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!     *CICOVER*   SEA ICE COVER. 
!     *CITHICK*   SEA ICE THICKNESS. 
!     *NBLKS*     INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!     *NBLKE*     INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!     *IREAD*     PROCESSOR WHICH WILL ACCESS THE FILE ON DISK

!     METHOD.
!     -------
!     THE READING IS ONLY DONE ON PE 1, THEREFORE
!     THE RELEVANT INFORMATION IS SENT TO THE OTHER PE'S USING 
!     MPDISTRIBSCFLD

!     EXTERNALS.
!     ----------
!     BUILDSTRESS
!     GRSTNAME
!     MPDISTRIBSCFLD
!     MPL_BARRIER
!     READSTRESS

!     REFERENCE.
!     ----------
!     NONE
! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : NREAL    ,LRSTPARALR

      USE YOWCURR  , ONLY : U        ,V
!!!                !!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE YOWGRID  , ONLY : IJS      ,IJL
      USE YOWMESPAS, ONLY : LGRIBIN  ,LWAVEWIND
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NSUP     ,KTAG
      USE YOWPARAM , ONLY : NANG     ,NFRE_RED ,NIBLO

      USE YOWREFD  , ONLY : THDD     ,THDC     ,SDOT

      USE YOWSTAT  , ONLY : CDATEA   ,CDATEF   ,CDTPRO   ,IREFRA   ,    &
     &            NPROMA_WAM,LNSESTART
      USE YOWTEST  , ONLY : IU06
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH    ,LRESTARTED
      USE YOWUBUF  , ONLY : LUPDTWGHT
      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
      USE MPL_MODULE

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "buildstress.intfb.h"
#include "expand_string.intfb.h"
#include "grstname.intfb.h"
#include "mpdistribscfld.intfb.h"
#include "mpexchng.intfb.h"
#include "propdot.intfb.h"
#include "readstress.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      INTEGER(KIND=JWIM),DIMENSION(NPROC), INTENT(IN) :: NBLKS, NBLKE

      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(INOUT) :: U10OLD, THWOLD, USOLD, TAUW, TAUWDIR, Z0OLD
      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(INOUT) :: ROAIRO, ZIDLOLD, CICOVER, CITHICK


      INTEGER(KIND=JWIM) :: IBUFLENGTH
      INTEGER(KIND=JWIM) :: IFLD, IJ, KCOUNT, IC
      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA
      INTEGER(KIND=JWIM) :: IJINF, IJSUP
      INTEGER(KIND=JWIM) :: LNAME 
      INTEGER(KIND=JWIM),ALLOCATABLE,DIMENSION(:) :: IBUF

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:,:) :: RFIELD

      CHARACTER(LEN= 14) :: ZERO
      CHARACTER(LEN=296) :: FILENAME

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('GETSTRESS',0,ZHOOK_HANDLE)

      ZERO = ' '
      NPROMA=NPROMA_WAM

!     READ RESTART FILE FROM PE IREAD

      IF (LGRIBIN .AND. .NOT.LRESTARTED) THEN
!       GRIB RESTART
!       CREATES WIND AND STRESS FIELDS FROM GRIB WINDS AND DRAG COEFFICIENT.
        CALL BUILDSTRESS(IJS, IJL,                                      &
     &                   U10OLD(IJS), THWOLD(IJS),                      &
     &                   USOLD(IJS), TAUW(IJS), TAUWDIR(IJS),           &
     &                   Z0OLD(IJS),                                    &
     &                   ROAIRO(IJS), ZIDLOLD(IJS),                     &
     &                   CICOVER(IJS), CITHICK(IJS),                    &
     &                   IREAD)


      ELSE

!       BINARY RESTART

        IF (LRSTPARALR) THEN
          IJINF=IJS
          IJSUP=IJL
        ELSE
          IJINF=1
          IJSUP=NIBLO
        ENDIF

        ALLOCATE(RFIELD(IJINF:IJSUP,NREAL))

        CALL GRSTNAME(CDATEA,CDATEF,'LAW',ICPLEN,CPATH,FILENAME)
        IF (LRSTPARALR) THEN
!          RESTART FILES FROM ALL PE's
           LNAME = LEN_TRIM(FILENAME)
           FILENAME=FILENAME(1:LNAME)//'.%p_%n'
           CALL EXPAND_STRING(IRANK, NPROC, 0, 0, FILENAME, 1)
        ENDIF

        IF (LRSTPARALR .OR. IRANK == IREAD) THEN
          CALL READSTRESS(IJINF, IJSUP, NREAL, RFIELD, FILENAME, LRSTPARALR)
        ENDIF

        IF (.NOT.LRSTPARALR) THEN
!       DISTRIBUTE THE DIFFERENT CONTRIBUTIONS TO THE OTHER PE'S
          DO IFLD=1,NREAL
            CALL MPDISTRIBSCFLD(IREAD, KTAG, NBLKS, NBLKE,RFIELD(:,IFLD))
            KTAG=KTAG+1
          ENDDO
        ENDIF

!       KEEP CORRESPONDING CONTRIBUTION 
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ)
        DO JKGLO=IJS,IJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL)
          DO IJ=KIJS,KIJL
            U10OLD(IJ)=RFIELD(IJ,1)
            THWOLD(IJ)=RFIELD(IJ,2)
            USOLD(IJ)=RFIELD(IJ,3)
            TAUW(IJ)=RFIELD(IJ,4)
            TAUWDIR(IJ)=RFIELD(IJ,5)
            Z0OLD(IJ)=RFIELD(IJ,6)
            ROAIRO(IJ)=RFIELD(IJ,7)
            ZIDLOLD(IJ)=RFIELD(IJ,8)
            CICOVER(IJ)=RFIELD(IJ,9)
            CITHICK(IJ)=RFIELD(IJ,10)
!           for U and V see below
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO

        IF (ALLOCATED(U) .AND. ALLOCATED(V)) THEN
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ)
          DO JKGLO=IJS,IJL,NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,IJL)
            DO IJ=KIJS,KIJL
              U(IJ)=RFIELD(IJ,11)
              V(IJ)=RFIELD(IJ,12)
            ENDDO
          ENDDO
!$OMP     END PARALLEL DO
          U(NSUP+1)=0.0_JWRB
          V(NSUP+1)=0.0_JWRB

!!!       U AND V MUST ALSO BE DEFINED OVER THE HALO
          CALL MPEXCHNG(U, 1, 1)
          CALL MPEXCHNG(V, 1, 1)

          IF (IREFRA .NE. 0) THEN
!         RE-COMPUTE REFRACTION TERMS
            IF (.NOT.ALLOCATED(THDC)) ALLOCATE(THDC(IJS:IJL,NANG))
            IF (.NOT.ALLOCATED(THDD)) ALLOCATE(THDD(IJS:IJL,NANG))
            IF (.NOT.ALLOCATED(SDOT)) ALLOCATE(SDOT(IJS:IJL,NANG,NFRE_RED))

            NPROMA=NPROMA_WAM
!$OMP       PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
            DO JKGLO = IJS, IJL, NPROMA
              KIJS=JKGLO
              KIJL=MIN(KIJS+NPROMA-1,IJL)
              CALL PROPDOT(KIJS, KIJL, THDC(KIJS:KIJL,:), THDD(KIJS:KIJL,:), SDOT(KIJS:KIJL,:,:))
            ENDDO
!$OMP       END PARALLEL DO

          END IF
  
!         SET LOGICAL TO RECOMPUTE THE WEIGHTS IN CTUW.
          LUPDTWGHT=.TRUE.

        ENDIF

        DEALLOCATE(RFIELD)

      ENDIF

!     BROADCAST THE 4 DATES FROM RESTART FILE TO THE OTHER PE'S
!     AS WELL AS LWAVEWIND
      IF (.NOT. (IREAD == 0 .OR. NPROC == 1)) THEN
        IBUFLENGTH=57+1
        ALLOCATE(IBUF(IBUFLENGTH))
        IF (IRANK == IREAD) THEN
          IF (LWAVEWIND) THEN
            IBUF(1)=1
          ELSE
            IBUF(1)=0
          ENDIF
          KCOUNT=1
          DO IC=1,14
            KCOUNT=KCOUNT+1
            READ(CDTPRO(IC:IC),'(I1)') IBUF(KCOUNT)
          ENDDO
          DO IC=1,14
            KCOUNT=KCOUNT+1
            READ(CDATEWO(IC:IC),'(I1)') IBUF(KCOUNT)
          ENDDO
          DO IC=1,14
            KCOUNT=KCOUNT+1
            READ(CDAWIFL(IC:IC),'(I1)') IBUF(KCOUNT)
          ENDDO
          DO IC=1,14
            KCOUNT=KCOUNT+1
            READ(CDATEFL(IC:IC),'(I1)') IBUF(KCOUNT)
            IF (KCOUNT > IBUFLENGTH) THEN
              WRITE (IU06,*) ' '
              WRITE (IU06,*) ' **************************************'
              WRITE (IU06,*) ' *                                    *'
              WRITE (IU06,*) ' *   FATAL ERROR IN SUB. GETSTRESS:   *'
              WRITE (IU06,*) ' *   ==============================   *'
              WRITE (IU06,*) ' *                                    *'
              WRITE (IU06,*) ' * IBUFLENGTH TOO SMALL !!!           *'
              WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS   *'
              WRITE (IU06,*) ' *                                    *'
              WRITE (IU06,*) ' **************************************'
              CALL ABORT1
            ENDIF
          ENDDO
        ENDIF

        CALL GSTATS(621,0)
        CALL MPL_BROADCAST(IBUF(1:IBUFLENGTH), KROOT=IREAD,             &
     &                     KTAG=KTAG,CDSTRING='GETSTRESS:')
        CALL GSTATS(621,1)

        KTAG=KTAG+1

        IF (IRANK /= IREAD) THEN
          IF (IBUF(1) == 1) THEN
            LWAVEWIND=.TRUE. 
          ELSE
            LWAVEWIND=.FALSE. 
          ENDIF
          KCOUNT=1
          DO IC=1,14
            KCOUNT=KCOUNT+1
            WRITE(CDTPRO(IC:IC),'(I1)') IBUF(KCOUNT)
          ENDDO
          DO IC=1,14
            KCOUNT=KCOUNT+1
            WRITE(CDATEWO(IC:IC),'(I1)') IBUF(KCOUNT)
          ENDDO
          DO IC=1,14
            KCOUNT=KCOUNT+1
            WRITE(CDAWIFL(IC:IC),'(I1)') IBUF(KCOUNT)
          ENDDO
          DO IC=1,14
            KCOUNT=KCOUNT+1
            WRITE(CDATEFL(IC:IC),'(I1)') IBUF(KCOUNT)
            IF (KCOUNT > IBUFLENGTH) THEN
              WRITE (IU06,*) ' '
              WRITE (IU06,*) ' **************************************'
              WRITE (IU06,*) ' *                                    *'
              WRITE (IU06,*) ' *   FATAL ERROR IN SUB. GETSTRESS:   *'
              WRITE (IU06,*) ' *   ==============================   *'
              WRITE (IU06,*) ' *                                    *'
              WRITE (IU06,*) ' * IBUFLENGTH TOO SMALL !!!           *'
              WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS   *'
              WRITE (IU06,*) ' *                                    *'
              WRITE (IU06,*) ' **************************************'
              CALL ABORT1
            ENDIF
          ENDDO
        ENDIF

        DEALLOCATE(IBUF)

      ENDIF

      IF (CDTPRO == '00000000000000') CDTPRO = ZERO
      IF (CDATEWO == '00000000000000')CDATEWO = ZERO
      IF (CDAWIFL == '00000000000000')CDAWIFL = ZERO
      IF (CDATEFL == '00000000000000') CDATEFL= ZERO


      IF (LNSESTART .AND. .NOT.LRESTARTED) THEN
!       WHEN INITAL SPECTRA SET TO NOISE LEVEL,
!       RESET WAVE INDUCED STRESS TO ZERO 
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ)
        DO JKGLO=IJS, IJL, NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL)
          DO IJ=KIJS,KIJL
            TAUW(IJ)=0.0_JWRB
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO
      ENDIF


      WRITE(IU06,*) ''
      WRITE(IU06,*) ' WIND AND STRESS FILES READ IN....... CDTPRO  = ', CDTPRO
      CALL FLUSH (IU06)

      IF (LHOOK) CALL DR_HOOK('GETSTRESS',1,ZHOOK_HANDLE)

      END SUBROUTINE GETSTRESS
