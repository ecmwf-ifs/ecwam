      SUBROUTINE GETSTRESS(IJS, IJL, BLK2LOC,                &
     &                     WVENVI, FF_NOW, NEMO2WAM,         &
     &                     NBLKS, NBLKE, IREAD)

! ----------------------------------------------------------------------

!*    PURPOSE.
!     --------
!     READS RESTART WIND AND STRESS FIELDS FROM DISK OR USE GRIB WIND
!     AND DRAG COEFFICIENT TO RECONSTRUCT ALL WIND AND STRESS FIELDS.

!**   INTERFACE.
!     ----------
!     *CALL**GETSTRESS(IJS, IJL, BLK2LOC, WVENVI, FF_NOW, NBLKS, NBLKE, IREAD)
!     *IJS*       INDEX OF FIRST GRIDPOINT.
!     *IJL*       INDEX OF LAST GRIDPOINT.
!     *BLK2LOC*   POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!     *WVENVI*    WAVE ENVIRONMENT
!     *FF_NOW*    FORCING FIELDS AT CURRENT TIME
!     *NEMO2WAM*  FIELDS FRON OCEAN MODEL to WAM
!     *NBLKS*     INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!     *NBLKE*     INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!     *IREAD*     PROCESSOR WHICH WILL ACCESS THE FILE ON DISK

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDLOC, ENVIRONMENT, FORCING_FIELDS, OCEAN2WAVE

      USE YOWCOUT  , ONLY : NREAL    ,LRSTPARALR
      USE YOWMESPAS, ONLY : LGRIBIN  ,LWAVEWIND
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NSUP     ,KTAG
      USE YOWPARAM , ONLY : NIBLO
      USE YOWREFD  , ONLY : LLUPDTTD
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
#include "readstress.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      TYPE(WVGRIDLOC), DIMENSION(IJS:IJL), INTENT(IN) :: BLK2LOC
      TYPE(ENVIRONMENT), DIMENSION(IJS:IJL), INTENT(INOUT) :: WVENVI
      TYPE(FORCING_FIELDS), DIMENSION(IJS:IJL), INTENT(INOUT) :: FF_NOW
      TYPE(OCEAN2WAVE), DIMENSION(IJS:IJL), INTENT(INOUT) :: NEMO2WAM
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      INTEGER(KIND=JWIM),DIMENSION(NPROC), INTENT(IN) :: NBLKS, NBLKE


      INTEGER(KIND=JWIM) :: IBUFLENGTH
      INTEGER(KIND=JWIM) :: IFLD, IJ, KCOUNT, IC
      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA
      INTEGER(KIND=JWIM) :: IJINF, IJSUP
      INTEGER(KIND=JWIM) :: LNAME 
      INTEGER(KIND=JWIM) :: IFCST
      INTEGER(KIND=JWIM),ALLOCATABLE,DIMENSION(:) :: IBUF

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:,:) :: RFIELD

      CHARACTER(LEN= 14) :: ZERO
      CHARACTER(LEN=296) :: FILENAME

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('GETSTRESS',0,ZHOOK_HANDLE)

ASSOCIATE(UCUR => WVENVI%UCUR, &
 &        VCUR => WVENVI%VCUR, &
 &        WSWAVE => FF_NOW%WSWAVE, &
 &        WDWAVE => FF_NOW%WDWAVE, &
 &        UFRIC => FF_NOW%UFRIC, &
 &        Z0M => FF_NOW%Z0M, &
 &        Z0B => FF_NOW%Z0B, &
 &        TAUW => FF_NOW%TAUW, &
 &        TAUWDIR => FF_NOW%TAUWDIR, &
 &        AIRD => FF_NOW%AIRD, &
 &        WSTAR => FF_NOW%WSTAR, &
 &        CICOVER => FF_NOW%CICOVER, &
 &        CITHICK => FF_NOW%CITHICK )


      ZERO = ' '
      NPROMA=NPROMA_WAM

      Z0B(:) = 0.0_JWRB

!     READ RESTART FILE FROM PE IREAD

      IF (LGRIBIN .AND. .NOT.LRESTARTED) THEN
!       GRIB RESTART
!       CREATES WIND AND STRESS FIELDS FROM GRIB WINDS AND DRAG COEFFICIENT.
        CALL BUILDSTRESS(IJS, IJL, BLK2LOC,                &
     &                   WVENVI, FF_NOW, NEMO2WAM,         &
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

        IFCST = 0
        CALL GRSTNAME(CDATEA,CDATEF,IFCST,'LAW',ICPLEN,CPATH,FILENAME)
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
            WSWAVE(IJ)=RFIELD(IJ,1)
            WDWAVE(IJ)=RFIELD(IJ,2)
            UFRIC(IJ)=RFIELD(IJ,3)
            TAUW(IJ)=RFIELD(IJ,4)
            TAUWDIR(IJ)=RFIELD(IJ,5)
            Z0M(IJ)=RFIELD(IJ,6)
            AIRD(IJ)=RFIELD(IJ,7)
            WSTAR(IJ)=RFIELD(IJ,8)
            CICOVER(IJ)=RFIELD(IJ,9)
            CITHICK(IJ)=RFIELD(IJ,10)
            UCUR(IJ)=RFIELD(IJ,11)
            VCUR(IJ)=RFIELD(IJ,12)
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO

        IF (IREFRA /= 0) THEN
          LLUPDTTD = .TRUE.
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

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('GETSTRESS',1,ZHOOK_HANDLE)

END SUBROUTINE GETSTRESS
