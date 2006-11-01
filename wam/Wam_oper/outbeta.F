      SUBROUTINE OUTBETA (IJS, IJL, IG, IGL,USNEW,Z0NEW,NSTART,NEND,
     &                    BETAG, NLONW, NLATW)

! ----------------------------------------------------------------------

!**** *OUTBETA* - OUTPUT OF CHARNOCK CONSTANT FOR ATMOSPHERIC
!****             MODEL IN GRIDDED FORMAT.

!     P.JANSSEN      KNMI/ECMWF  JANUARY 1992
!     J.BIDLOT       ECMWF       FEBRUARY 1996  MESSAGE PASSING


!*    PURPOSE.
!     --------

!       OUTPUT OF CHARNOCK CONSTANT FOR ATMOSPHERIC MODEL


!**   INTERFACE.
!     ----------

!       *CALL* *OUTBETA (IJS, IJL, IG, IGL,USNEW,Z0NEW,NSTART,NEND,
!    &                   BETAG, NLONW, NLATW)*
!         *IJS*    - INDEX OF FIRST GRIDPOINT.
!         *IJL*    - INDEX OF LAST GRIDPOINT.
!         *IG*     - BLOCK NUMBER.
!         *USNEW*     NEW FRICTION VELOCITY IN M/S.
!         *Z0NEW*     ROUGHNESS LENGTH IN M.
!         *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!         *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!         *BETAG*   GRIDDED CHARNOCK FIELD
!         *NLONW*   FIRST DIMENSION OF BETAG
!         *NLATW*   SECOND DIMENSION OF BETAG
!         

!     EXTERNALS.
!     ----------

!       *MAKEGRID*  - GENERATES GRIDDED FIELDS.
!       *MPBCASTSCFLD DISTRIBUTE TO ALL PE's THE CONTRIBUTION OF A BLOCK
!                     FILE

!     METHOD.
!     -------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

!*    *PARAMETER*  FOR ARRAY DIMENSIONS.


      USE YOWPARAM , ONLY : NGX      ,NGY      ,NIBLO
      USE YOWPCONS , ONLY : G        ,EPSUS
      USE YOWICE   , ONLY : LICERUN  ,ICEMASK
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : NPROC    ,NINF     ,NSUP
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
      REAL,DIMENSION(NINF:NSUP) :: USNEW,Z0NEW
      REAL ZHOOK_HANDLE

! ----------------------------------------------------------------------
      DIMENSION BETAG(NLONW,NLATW)
! ----------------------------------------------------------------------

      LOGICAL :: LSQRT
      DIMENSION BETA(NIBLO)

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('OUTBETA',0,ZHOOK_HANDLE)

!*    1. COMPUTE CHARNOCK 'CONSTANT' BETA.
!        ---------------------------------

      RCHAR=0.018

      DO IJ = IJS,IJL
        IF (LICERUN.AND.ICEMASK(IJ,IG).EQ.0) THEN
          BETA(IJ) = RCHAR
        ELSE
          BETA(IJ) = G*Z0NEW(IJ)/MAX(USNEW(IJ)**2,EPSUS)
        ENDIF
      ENDDO

!*    1.1 STORE PARAMETERS FOR ACTUAL BLOCK IN GRID ARRAYS.
!         -------------------------------------------------

      LSQRT=.FALSE.

      IF(LMESSPASS) THEN
        CALL MPBCASTSCFLD(20,NSTART,NEND,BETA)
        CALL MAKEGRID (BETA , BETAG  , 1, NEND(NPROC), IG, LSQRT, RCHAR)
      ELSE
        CALL MAKEGRID (BETA , BETAG  , IJS, IJL, IG, LSQRT, RCHAR)
      ENDIF

      IF (LHOOK) CALL DR_HOOK('OUTBETA',1,ZHOOK_HANDLE)

! ----------------------------------------------------------------------
      RETURN
      END SUBROUTINE OUTBETA
