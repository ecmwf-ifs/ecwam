      SUBROUTINE BUILDSTRESS(IJS, IJL,                                &
     &                       WSWAVE, WDWAVE,                          &
     &                       UFRIC, TAUW, TAUWDIR, Z0M,               &
     &                       AIRD, WSTAR,                        &
     &                       CICOVER, CITHICK,                        &
     &                       IREAD)

! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF   APRIL 1998 

!     J. BIDLOT    ECMWF   FEBRUARY 1999 TAUT --> SQRT(TAUT)

!     S. ABDALLA   ECMWF   OCTOBER 1999 MODIFICATION THE CALL TO GETWND
 
!     J. BIDLOT    ECMWF   AUGUST 2008 : MAKE IT MORE PARALLEL.

!*    PURPOSE.
!     --------
!     CREATES WIND AND STRESS FIELDS FROM GRIB WINDS AND CD.

!**   INTERFACE.
!     ----------
!     CALL *BUILDSTRESS*(IJS, IJL,
!    &                   WSWAVE,WDWAVE,UFRIC,TAUW,TAUWDIR,Z0M,AIRD,
!    &                   AIRD, WSTAR, CICOVER, CITHICK,
!    &                   IREAD)*
!     *IJS*        INDEX OF FIRST GRIDPOINT
!     *IJL*        INDEX OF LAST GRIDPOINT
!     *WSWAVE *    WIND SPEED.
!     *WDWAVE *    WIND DIRECTION (RADIANS).
!     *UFRIC*      FRICTION VELOCITY.
!     *TAUW*       WAVE STRESS MAGNITUDE.
!     *TAUWDIR*    WAVE STRESS DIRECTION.
!     *Z0M*      ROUGHNESS LENGTH IN M.
!     *RAD0OLD*    AIR DENSITY IN KG/M3.
!     *RWSTAR0OLD* CONVECTIVE VELOCITY.
!     *CICOVER*    SEA ICE COVER.
!     *CITHICK*    SEA ICE THICKNESS.
!     *IREAD*      PROCESSOR WHICH WILL ACCESS THE FILE ON DISK

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------
!     *ABORT1*
!     *GETWND*
!     *READWGRIB*

!     REFERENCE.
!     ----------
!     NONE
! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWCOU    , LLCAPCHNK, LLGCBZ0  ,            &
     &                      LWNEMOCOUCIC, LWNEMOCOUCIT
      USE YOWMPP   , ONLY : NPROC    ,IRANK    ,KTAG
      USE YOWMESPAS, ONLY : LMESSPASS,LNOCDIN  ,LWAVEWIND
      USE YOWPCONS , ONLY : G        ,GM1      ,ROAIR    ,EPSUS    ,EPSU10
      USE YOWPHYS  , ONLY : ALPHA    ,XKAPPA   ,XNLEV    ,RNUM
      USE YOWSTAT  , ONLY : CDATEA   ,CDTPRO   ,NPROMA_WAM
      USE YOWTEST  , ONLY : IU06
      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL  ,FIELDG   ,    &
     &                      NXFF     ,NYFF
      USE YOWNEMOFLDS,ONLY: NEMOCICOVER, NEMOCITHICK

      USE MPL_MODULE
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "cdustarz0.intfb.h"
#include "chnkmin.intfb.h"
#include "getwnd.intfb.h"
#include "init_fieldg.intfb.h"
#include "readwgrib.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD

      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: WSWAVE, WDWAVE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: UFRIC, Z0M, TAUW, TAUWDIR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: AIRD, WSTAR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: CICOVER,CITHICK


      INTEGER(KIND=JWIM) :: ICODE_WND
      INTEGER(KIND=JWIM) :: ILEN, LIU, IPARAM, KZLEVUWAVE, KZLEVCD
      INTEGER(KIND=JWIM) :: IJ
      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA
      INTEGER(KIND=JWIM) :: NWAVEWIND(1)

      REAL(KIND=JWRB) :: RUSE
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: TEMPXNLEV, CDSQRTINV, Z0TOT, USTAR, CHARNOCKOG
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: CD
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: ALPHAOG

      CHARACTER(LEN=24) :: FILNM

      LOGICAL :: LLONLYPOS
      LOGICAL :: LWNDFILE, LCLOSEWND
      LOGICAL :: LCR
      LOGICAL :: LLALLOC_ONLY, LLINIALL, LLOCAL

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('BUILDSTRESS',0,ZHOOK_HANDLE)

      CDATEWO = ' '
      CDAWIFL = ' '
      CDATEFL = ' '
      CDTPRO = CDATEA
      LCR=.FALSE.

!     GETWND AND READWGRIB REQUIRES FIELDG TO BE ALLOCATED !
      LLALLOC_ONLY=.FALSE.
      LLINIALL=.TRUE.
      LLOCAL=.TRUE.
      CALL INIT_FIELDG(LLALLOC_ONLY,LLINIALL,LLOCAL)

!     1.1 GET ATMOSPHERIC MODEL FORCINGS FIELDS 
!         -------------------------------------

      ILEN=1
      LWNDFILE=.TRUE.
      LCLOSEWND=.TRUE.

      IF (LWNEMOCOUCIC .OR. LWNEMOCOUCIT) THEN
        ALLOCATE(NEMOCICOVER(IJS:IJL),NEMOCITHICK(IJS:IJL))
        NEMOCICOVER(IJS:IJL)=0.0_JWRB
        NEMOCITHICK(IJS:IJL)=0.0_JWRB
      ENDIF

      CALL GETWND (IJS, IJL,                                &
     &             WSWAVE, UFRIC,                           &
     &             WDWAVE,                                  &
     &             AIRD, WSTAR,                             &
     &             CICOVER, CITHICK,                        &
     &             CDTPRO, LWNDFILE, LCLOSEWND, IREAD,      &
     &             LCR, ICODE_WND)

      IF (LWNEMOCOUCIC .OR. LWNEMOCOUCIT) THEN
        DEALLOCATE(NEMOCICOVER,NEMOCITHICK)
      ENDIF


!     1.2 USE DATA FROM A FILE CONTAINING WIND SPEED MODIFIED BY
!         ----------------------------------------------------
!     A PREVIOUS WAVE MODEL RUN ON WHICH THIS RESTART IS BASED
!     IF IT EXISTS, OTHERWISE IT WILL USE THE INFORMATION FROM
!     THE ATMOSPHERIC MODEL WIND SPEED ARCHIVE.


      FILNM='uwavein'
      LIU = LEN_TRIM(FILNM)
      FILNM=FILNM(1:LIU)
      IF (IREAD == IRANK) THEN
        INQUIRE(FILE=FILNM,EXIST=LWAVEWIND)
        IF (LWAVEWIND) THEN
          NWAVEWIND(1)=1
        ELSE
          NWAVEWIND(1)=0
        ENDIF
      ENDIF

!     USE MESSAGE PASSING TO SEND FILE STATUS TO THE OTHER PE'S
      CALL GSTATS(696,0)
      IF (NPROC > 1) THEN
        CALL MPL_BROADCAST(NWAVEWIND,KROOT=IREAD,KTAG=KTAG,             &
     &    CDSTRING='BUILDSTRESS :')
        KTAG=KTAG+1
        IF (NWAVEWIND(1) == 1) THEN
          LWAVEWIND=.TRUE.
        ELSE
          LWAVEWIND=.FALSE.
        ENDIF
      ENDIF
      CALL GSTATS(696,1)

      IF (LWAVEWIND) THEN
        IPARAM=245
        LLONLYPOS=.FALSE.
        CALL READWGRIB(IU06, FILNM, IPARAM, CDTPRO, IJS, IJL,         &
     &                 WSWAVE(IJS), KZLEVUWAVE, LLONLYPOS, IREAD)

        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' A DATA FILE CONTAINING WIND SPEED INFORMATION'
        WRITE(IU06,*) ' AS PROVIDED BY A PREVIOUS WAVE MODEL RUN'
        WRITE(IU06,*) ' WAS USED TO UPDATE THE INPUT ATMOSPHERIC WINDS'
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' THE INPUT WINDS AND DRAG COEFFICIENT ARE FOUND'
        WRITE(IU06,*) ' TO HAVE BEEN DETERMINED FOR HEIGHT AT ',        &
     &                  KZLEVUWAVE,' m'


      ELSE

        KZLEVUWAVE=10
        WRITE(IU06,*) ' '
        WRITE(IU06,*) '          !!!! NOTE !!!!'
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' NO INFORMATION ON WIND SPEEDS FROM THE'
        WRITE(IU06,*) ' WAVE MODEL WAS PROVIDED'
        WRITE(IU06,*) ' NO UPDATE OF THE INPUT ATMOSPHERIC WINDS'
        WRITE(IU06,*) ' WAS POSSIBLE'
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' THE INPUT WIND AND DRAG COEFFICIENT ARE ASSUMED'
        WRITE(IU06,*) ' TO HAVE BEEN DETERMINED FOR HEIGHT AT ',        &
     &                  KZLEVUWAVE,' m'

      ENDIF
 
!     TEST WHETHER THE HEIGHT OF THE INPUT WINDS IS THE SAME AS DEFINED
!     FOR THE REST OF THE RUN (SEE CALL TO STRESS IN PREPROC). IF NOT
!     RECOMPUTE THE TABLE. THE ORIGINAL TABLE WILL BE SWAP BACK AT THE
!     END OF THE ROUTINE.

      IF (KZLEVUWAVE /= NINT(XNLEV)) THEN
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' THE REFERENCE HEIGHT TO BE USED IN WAMODEL'
        WRITE(IU06,*) ' ',XNLEV 
        WRITE(IU06,*) ' IS DIFFERENT THAN THE INPUT FIELDS HEIGHT'
        WRITE(IU06,*) ' ',FLOAT(KZLEVUWAVE)
        WRITE(IU06,*) ' THE NECESSARY ADJUSTMENTS WILL BE MADE'
        WRITE(IU06,*) ' TO DETERMINE THE INITIAL FIELDS.'
        WRITE(IU06,*) ' '

        TEMPXNLEV=KZLEVUWAVE
      ELSE
        TEMPXNLEV=XNLEV
      ENDIF

!     1.3 INITIALISE CD USING THE FRICTION VELOCITY FOR TAUW=0.
!         ----------------------------------------------------

! Mod for OPENMP
      NPROMA=NPROMA_WAM
      CALL GSTATS(1444,0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL,IJ) 
      DO JKGLO=IJS,IJL,NPROMA
        KIJS=JKGLO
        KIJL=MIN(KIJS+NPROMA-1,IJL)

        CALL CDUSTARZ0 (KIJS, KIJL, WSWAVE(KIJS), TEMPXNLEV, CD(KIJS), UFRIC(KIJS), Z0M(KIJS))

        DO IJ=KIJS,KIJL
          TAUW(IJ) = 0.1_JWRB * UFRIC(IJ)**2
          TAUWDIR(IJ) = WDWAVE(IJ)
        ENDDO

      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1444,1)

!     1.4  GET DRAG COEFFICIENT
!          --------------------
      IF (.NOT.LNOCDIN) THEN
        IPARAM=233
        LLONLYPOS=.TRUE.
        FILNM='cdwavein'
!       !!!! CD was initialised above !!!!
        CALL READWGRIB(IU06, FILNM, IPARAM, CDTPRO, IJS, IJL,         &
     &                 CD(IJS), KZLEVCD, LLONLYPOS, IREAD)

!       TEST REFERENCE LEVEL FOR UWAVE AND CD

        IF (KZLEVUWAVE /= 0 .AND. KZLEVCD /= 0) THEN
          IF (KZLEVUWAVE /= KZLEVCD) THEN
            WRITE(IU06,*)'************************************'
            WRITE(IU06,*)'*                                  *'
            WRITE(IU06,*)'* FATAL ERROR IN SUB BUILDSTRESS   *'
            WRITE(IU06,*)'* ==============================   *'
            WRITE(IU06,*)'* REFERENCE LEVEL FOR CD AND UWAVE *'
            WRITE(IU06,*)'* DIFFER. PROGRAM WILL ABORT       *'
            WRITE(IU06,*)'* REF. LEVEL CD = ',KZLEVCD
            WRITE(IU06,*)'* REF. LEVEL UWAVE = ',KZLEVUWAVE
            WRITE(IU06,*)'*                                  *'
            WRITE(IU06,*)'************************************'
            CALL ABORT1
          ENDIF
        ENDIF

!       1.5 COMPUTE TAUW, UFRIC AND Z0M
!           ----------------------------

        IF (LLCAPCHNK) THEN
          DO IJ=IJS,IJL
            ALPHAOG(IJ)= CHNKMIN(WSWAVE(IJ))*GM1
          ENDDO
        ELSE
          DO IJ=IJS,IJL
            ALPHAOG(IJ)= ALPHA*GM1
          ENDDO
        ENDIF

        IF (LLGCBZ0) THEN
          RUSE = 0.0_JWRB
        ELSE 
          RUSE = 1.0_JWRB
        ENDIF

! Mod for OPENMP
          NPROMA=NPROMA_WAM
          CALL GSTATS(1444,0)
!$OMP     PARALLEL DO SCHEDULE(DYNAMIC,1) & 
!$OMP&    PRIVATE(JKGLO,KIJS,KIJL,IJ,CDSQRTINV,Z0TOT,USTAR,CHARNOCKOG)
          DO JKGLO=IJS,IJL,NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,IJL)
            DO IJ=KIJS,KIJL
!!            UFRIC WILL FIRST CONTAIN ITS SQUARE
!!            THE NUMERICAL RELATION BETWEEN UFRIC AND WSWAVE SHOULD
!!            ALWAYS BE AS IN OUTGRID
              UFRIC(IJ) = CD(IJ)*MAX(WSWAVE(IJ)**2,EPSU10)
              UFRIC(IJ) = MAX(UFRIC(IJ),EPSUS)
              USTAR = SQRT(UFRIC(IJ))
              CDSQRTINV = MIN(1._JWRB/SQRT(CD(IJ)),100.0_JWRB)
              Z0TOT = TEMPXNLEV*EXP(-XKAPPA*CDSQRTINV)
!             Z0M ONLY CONTAINS CHARNOCK CONTRIBUTION (see taut_z0)
              Z0M(IJ) = MAX(Z0TOT - RUSE*RNUM/USTAR,ALPHAOG(IJ)*UFRIC(IJ))
              CHARNOCKOG = Z0M(IJ)/UFRIC(IJ)
              CHARNOCKOG = MAX(CHARNOCKOG,ALPHAOG(IJ))
              TAUW(IJ) = MAX(UFRIC(IJ)*(1._JWRB-(ALPHAOG(IJ)/CHARNOCKOG)**2),0._JWRB)
              TAUWDIR(IJ) = WDWAVE(IJ)
              UFRIC(IJ) = USTAR 
            ENDDO
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1444,1)
      ENDIF

      DEALLOCATE(FIELDG)

      WRITE(IU06,*) ' SUB. BUILDSTRESS: INPUT OF RESTART FILES DONE'

      IF (LHOOK) CALL DR_HOOK('BUILDSTRESS',1,ZHOOK_HANDLE)

      END SUBROUTINE BUILDSTRESS
