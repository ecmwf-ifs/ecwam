      SUBROUTINE BUILDSTRESS(MIJS, MIJL,                                &
     &                       U10OLD, THWOLD,                            &
     &                       USOLD, TAUW, Z0OLD,                        &
     &                       ROAIRO, ZIDLOLD,                           &
     &                       CICOVER, CITHICK,                          &
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
!     CALL *BUILDSTRESS*(MIJS, MIJL,
!    &                   U10OLD,THWOLD,USOLD,TAUW,Z0OLD,ROAIRO,
!    &                   ROAIRO, ZIDLOLD, CICOVER, CITHICK,
!    &                   IREAD)*
!     *MIJS*      INDEX OF FIRST GRIDPOINT
!     *MIJL*      INDEX OF LAST GRIDPOINT
!     *U10OLD*   WIND SPEED.
!     *THWOLD*   WIND DIRECTION (RADIANS).
!     *USOLD*    FRICTION VELOCITY.
!     *TAUW*     WAVE STRESS.
!     *Z0OLD*    ROUGHNESS LENGTH IN M.
!     *RAD0OLD*   AIR DENSITY IN KG/M3.
!     *RZIDL0OLD* Zi/L (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!     *CICOVER*   SEA ICE COVER.
!     *CITHICK*   SEA ICE THICKNESS.
!     *IREAD*     PROCESSOR WHICH WILL ACCESS THE FILE ON DISK

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------
!     *ABORT1*
!     *AIRSEA*
!     *GETWND*
!     *READWGRIB*

!     REFERENCE.
!     ----------
!     NONE
! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWCOU    ,ALPHA    ,XKAPPA   ,XNLEV,        &
     &                      RNUM     ,LLCAPCHNK,                        &
     &                      LWNEMOCOUCIC, LWNEMOCOUCIT
      USE YOWMPP   , ONLY : NPROC    ,IRANK    ,KTAG
      USE YOWMESPAS, ONLY : LMESSPASS,LNOCDIN  ,LWAVEWIND
      USE YOWPCONS , ONLY : G        ,ROAIR    ,EPSUS    ,EPSU10
      USE YOWSTAT  , ONLY : CDATEA   ,CDTPRO   ,NPROMA_WAM
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL  ,FIELDG   ,    &
     &                      NXFF     ,NYFF
      USE YOWNEMOFLDS,ONLY: NEMOCICOVER, NEMOCITHICK

      USE MPL_MODULE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "airsea.intfb.h"
#include "getwnd.intfb.h"
#include "init_fieldg.intfb.h"
#include "readwgrib.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: MIJS, MIJL
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD

      REAL(KIND=JWRB), DIMENSION(MIJS:MIJL), INTENT(OUT) :: U10OLD, THWOLD
      REAL(KIND=JWRB), DIMENSION(MIJS:MIJL), INTENT(OUT) :: USOLD, Z0OLD, TAUW
      REAL(KIND=JWRB), DIMENSION(MIJS:MIJL), INTENT(OUT) :: ROAIRO, ZIDLOLD
      REAL(KIND=JWRB), DIMENSION(MIJS:MIJL), INTENT(OUT) :: CICOVER,CITHICK

      INTEGER(KIND=JWIM) :: ICODE_WND
      INTEGER(KIND=JWIM) :: ILEN, LIU, IPARAM, ILEV, KZLEVUWAVE, KZLEVCD
      INTEGER(KIND=JWIM) :: IJ
      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA
      INTEGER(KIND=JWIM) :: NWAVEWIND(1)

      REAL(KIND=JWRB) :: CHNKMIN
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: TEMPXNLEV, CDINV, CDSQRTINV, Z0TOT, USTAR, CHARNOCKOG
      REAL(KIND=JWRB) :: GM1 
      REAL(KIND=JWRB), DIMENSION(MIJS:MIJL) :: CD
      REAL(KIND=JWRB), DIMENSION(MIJS:MIJL) :: ALPHAOG

      CHARACTER(LEN=24) :: FILNM

      LOGICAL :: LLONLYPOS
      LOGICAL :: LWNDFILE, LCLOSEWND
      LOGICAL :: LCR
      LOGICAL :: LLALLOC_ONLY, LLINIALL, LLOCAL
      LOGICAL :: LLADJHGT

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('BUILDSTRESS',0,ZHOOK_HANDLE)

      CDATEWO = ' '
      CDAWIFL = ' '
      CDATEFL = ' '
      CDTPRO = CDATEA
      LCR=.FALSE.
      GM1= 1.0_JWRB/G

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

      IF(LWNEMOCOUCIC.OR.LWNEMOCOUCIT) THEN
        ALLOCATE(NEMOCICOVER(MIJS:MIJL),NEMOCITHICK(MIJS:MIJL))
        NEMOCICOVER(MIJS:MIJL)=0.0_JWRB
        NEMOCITHICK(MIJS:MIJL)=0.0_JWRB
      ENDIF

      CALL GETWND (MIJS, MIJL,                                          &
     &             U10OLD(MIJS), USOLD(MIJS),                           &
     &             THWOLD(MIJS),                                        &
     &             ROAIRO(MIJS), ZIDLOLD(MIJS),                         &
     &             CICOVER(MIJS), CITHICK(MIJS),                        &
     &             CDTPRO, LWNDFILE, LCLOSEWND, IREAD,                  &
     &             LCR, ICODE_WND)

      IF(LWNEMOCOUCIC.OR.LWNEMOCOUCIT) THEN
        DEALLOCATE(NEMOCICOVER,NEMOCITHICK)
      ENDIF

      IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. GETWND DONE'


!     1.2 USE DATA FROM A FILE CONTAINING WIND SPEED MODIFIED BY
!         ----------------------------------------------------
!     A PREVIOUS WAVE MODEL RUN ON WHICH THIS RESTART IS BASED
!     IF IT EXISTS, OTHERWISE IT WILL USE THE INFORMATION FROM
!     THE ATMOSPHERIC MODEL WIND SPEED ARCHIVE.


      FILNM='uwavein'
      LIU = LEN_TRIM(FILNM)
      FILNM=FILNM(1:LIU)
      IF(IREAD.EQ.IRANK) THEN
        INQUIRE(FILE=FILNM,EXIST=LWAVEWIND)
        IF(LWAVEWIND) THEN
          NWAVEWIND(1)=1
        ELSE
          NWAVEWIND(1)=0
        ENDIF
      ENDIF

!     USE MESSAGE PASSING TO SEND FILE STATUS TO THE OTHER PE'S
      CALL GSTATS(696,0)
      IF(NPROC.GT.1) THEN
        CALL MPL_BROADCAST(NWAVEWIND,KROOT=IREAD,KTAG=KTAG,             &
     &    CDSTRING='BUILDSTRESS :')
        KTAG=KTAG+1
        IF(NWAVEWIND(1).EQ.1) THEN
          LWAVEWIND=.TRUE.
        ELSE
          LWAVEWIND=.FALSE.
        ENDIF
      ENDIF
      CALL GSTATS(696,1)

      IF(LWAVEWIND) THEN
        IPARAM=245
        LLONLYPOS=.FALSE.
        CALL READWGRIB(IU06, FILNM, IPARAM, CDTPRO, MIJS, MIJL,         &
     &                 U10OLD(MIJS), KZLEVUWAVE, LLONLYPOS, IREAD)
        IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. READWGRIB DONE FOR ',FILNM

        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' A DATA FILE CONTAINING WIND SPEED INFORMATION'
        WRITE(IU06,*) ' AS PROVIDED BY A PREVIOUS WAVE MODEL RUN'
        WRITE(IU06,*) ' WAS USED TO UPDATE THE INPUT ATMOSPHERIC WINDS'
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' THE INPUT WINDS AND DRAG COEFFICIENT ARE FOUND'
        WRITE(IU06,*) ' TO HAVE BEEN DETERMINED FOR HEIGHT AT ',        &
     &                  KZLEVUWAVE,' m'
        IF (ITEST.GT.0) CALL FLUSH(IU06) 

        ILEV=1

      ELSE

        ILEV=1
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

        IF (ITEST.GT.0) CALL FLUSH(IU06) 
      ENDIF
 
!     TEST WHETHER THE HEIGHT OF THE INPUT WINDS IS THE SAME AS DEFINED
!     FOR THE REST OF THE RUN (SEE CALL TO STRESS IN PREPROC). IF NOT
!     RECOMPUTE THE TABLE. THE ORIGINAL TABLE WILL BE SWAP BACK AT THE
!     END OF THE ROUTINE.

      IF (KZLEVUWAVE.NE.NINT(XNLEV(1))) THEN
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' THE REFERENCE HEIGHT TO BE USED IN WAMODEL'
        WRITE(IU06,*) ' ',XNLEV(1) 
        WRITE(IU06,*) ' IS DIFFERENT THAN THE INPUT FIELDS HEIGHT'
        WRITE(IU06,*) ' ',FLOAT(KZLEVUWAVE)
        WRITE(IU06,*) ' THE NECESSARY ADJUSTMENTS WILL BE MADE'
        WRITE(IU06,*) ' TO DETERMINE THE INITIAL FIELDS.'
        WRITE(IU06,*) ' '
        IF (ITEST.GT.0) CALL FLUSH(IU06) 

        LLADJHGT=.TRUE.
        TEMPXNLEV=XNLEV(1)
        XNLEV(1)=KZLEVUWAVE
      ELSE
        LLADJHGT=.FALSE.
      ENDIF

!     1.3 INITIALISE CD USING THE FRICTION VELOCITY FOR TAUW=0.
!         ----------------------------------------------------

! Mod for OPENMP
      NPROMA=NPROMA_WAM
      CALL GSTATS(1444,0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL,IJ,CDINV) 
      DO JKGLO=MIJS,MIJL,NPROMA
        KIJS=JKGLO
        KIJL=MIN(KIJS+NPROMA-1,MIJL)

        DO IJ=KIJS,KIJL
          TAUW(IJ)=0._JWRB
        ENDDO

        CALL AIRSEA (U10OLD(KIJS),TAUW(KIJS),USOLD(KIJS),               &
     &               Z0OLD(KIJS), KIJS, KIJL, ILEV, ICODE_WND)

        DO IJ=KIJS,KIJL
!!        THE NUMERICAL RELATION BETWEEN USOLD AND U10OLD SHOULD
!!        ALWAYS BE AS IN OUTGRID
          CDINV = MAX(U10OLD(IJ)**2,EPSU10)/MAX(USOLD(IJ)**2,EPSUS)
          CDINV = MIN(CDINV,10000.0_JWRB) 
          CD(IJ) = 1._JWRB/CDINV
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1444,1)
      IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. AIRSEA DONE AT 1'

!     1.4  GET DRAG COEFFICIENT
!          --------------------
      IF(.NOT.LNOCDIN) THEN
        IPARAM=233
        LLONLYPOS=.TRUE.
        FILNM='cdwavein'
!       !!!! CD was initialised above !!!!
        CALL READWGRIB(IU06, FILNM, IPARAM, CDTPRO, MIJS, MIJL,         &
     &                 CD(MIJS), KZLEVCD, LLONLYPOS, IREAD)
        IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. READWGRIB DONE FOR ',FILNM

!       TEST REFERENCE LEVEL FOR UWAVE AND CD

        IF(KZLEVUWAVE.NE.0.AND.KZLEVCD.NE.0) THEN
          IF(KZLEVUWAVE.NE.KZLEVCD) THEN
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

!       1.5 COMPUTE TAUW,USOLD AND Z0OLD
!           ----------------------------

        ILEV=1

        IF(LLCAPCHNK) THEN
          DO IJ=MIJS,MIJL
            ALPHAOG(IJ)= CHNKMIN(U10OLD(IJ))*GM1
          ENDDO
        ELSE
          DO IJ=MIJS,MIJL
            ALPHAOG(IJ)= ALPHA*GM1
          ENDDO
        ENDIF


! Mod for OPENMP
          NPROMA=NPROMA_WAM
          CALL GSTATS(1444,0)
!$OMP     PARALLEL DO SCHEDULE(DYNAMIC,1) & 
!$OMP&    PRIVATE(JKGLO,KIJS,KIJL,IJ,CDSQRTINV,Z0TOT,USTAR,CHARNOCKOG)
          DO JKGLO=MIJS,MIJL,NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,MIJL)
            DO IJ=KIJS,KIJL
!!            USOLD WILL FIRST CONTAIN ITS SQUARE
!!            THE NUMERICAL RELATION BETWEEN USOLD AND U10OLD SHOULD
!!            ALWAYS BE AS IN OUTGRID
              USOLD(IJ) = CD(IJ)*MAX(U10OLD(IJ)**2,EPSU10)
              USOLD(IJ) = MAX(USOLD(IJ),EPSUS)
              USTAR = SQRT(USOLD(IJ))
              CDSQRTINV = MIN(1._JWRB/SQRT(CD(IJ)),100.0_JWRB)
              Z0TOT = XNLEV(ILEV)*EXP(-XKAPPA*CDSQRTINV)
!             Z0OLD ONLY CONTAINS CHARNOCK CONTRIBUTION (see taut_z0)
              Z0OLD(IJ) = MAX(Z0TOT - RNUM/USTAR,ALPHAOG(IJ)*USOLD(IJ))
              CHARNOCKOG = Z0OLD(IJ)/USOLD(IJ)
              CHARNOCKOG = MAX(CHARNOCKOG,ALPHAOG(IJ))
              TAUW(IJ) = MAX(USOLD(IJ)*(1._JWRB-(ALPHAOG(IJ)/CHARNOCKOG)**2),0._JWRB)
              USOLD(IJ) = USTAR 
            ENDDO
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1444,1)
      ENDIF

      IF(LLADJHGT) THEN
        XNLEV(1)=TEMPXNLEV
      ENDIF

      DEALLOCATE(FIELDG)

      WRITE(IU06,*) ' SUB. BUILDSTRESS: INPUT OF RESTART FILES DONE'

      IF (LHOOK) CALL DR_HOOK('BUILDSTRESS',1,ZHOOK_HANDLE)

      END SUBROUTINE BUILDSTRESS
