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
      USE YOWPCONS , ONLY : G
      USE YOWICE   , ONLY : LICERUN  ,ICEMASK
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : NPROC    ,NINF     ,NSUP

! ----------------------------------------------------------------------

!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
      REAL,DIMENSION(NINF:NSUP) :: USNEW,Z0NEW

! ----------------------------------------------------------------------
      DIMENSION BETAG(NLONW,NLATW)
! ----------------------------------------------------------------------

      DIMENSION BETA(NIBLO)

! ----------------------------------------------------------------------

!*    1. COMPUTE CHARNOCK 'CONSTANT' BETA.
!        ---------------------------------

      DO IJ = IJS,IJL
        BETA(IJ) = G*Z0NEW(IJ)/(USNEW(IJ)**2 +0.0001)
      ENDDO

      IF (LICERUN) THEN
        DO IJ = IJS,IJL
          ZICE     = REAL(ICEMASK(IJ,IG))
          BETA(IJ) = ZICE*BETA(IJ)
        ENDDO
      ENDIF   

!*    1.1 STORE PARAMETERS FOR ACTUAL BLOCK IN GRID ARRAYS.
!         -------------------------------------------------

      IF(LMESSPASS) THEN

        CALL MPBCASTSCFLD(20,NSTART,NEND,BETA)

        CALL MAKEGRID (BETA , BETAG  , 1, NEND(NPROC), IG)

      ELSE
        CALL MAKEGRID (BETA , BETAG  , IJS, IJL, IG)
      ENDIF


!*    1.2 SET LANDVALUES TO CHARNOCK CONSTANT.
!         ------------------------------------

      IF (IG.EQ.IGL) THEN
        DO I=1,NGX
          DO J=1,NGY
            IF (BETAG(I,J).LE.0.0) THEN
              BETAG(I,J) = 0.018
            ENDIF
          ENDDO
        ENDDO
      ENDIF

! ----------------------------------------------------------------------
      RETURN
      END SUBROUTINE OUTBETA
