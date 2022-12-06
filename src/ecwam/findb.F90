      SUBROUTINE FINDB (NDIM, NBOUN, BLATB, BLNGB, IJARB,        &
     &                  NPROC, NSTART, NEND, IRANK)

! ----------------------------------------------------------------------

!**** *FINDB* - FIND BLOCK AND GRID POINT NUMBER.

!     R. PORTZ     MPI         15/01/1991
!     J. BIDLOT    ECMWF       JUNE 1996  MESSAGE PASSING

!*    PURPOSE.
!     -------

!       FIND BLOCK AND GRID POINT NUMBER FOR A GIVEN ARRAY
!       OF LONGITUDES AND LATITUDES.
!       IN CASE OF MESSAGE PASSING ARCHITECHTURE, THE SEARCH IS ONLY
!       DONE ON THE SUBREGION WHICH IS DEFINED ON THAT PE.

!**   INTERFACE.
!     ----------

!       *CALL* *FINDB (NDIM, NBOUN, BLATB, BLNGB, IJARB)*
!          *NDIM*    - DIMENSION OF ARRAYS.
!          *NBOUN*   - NUMBER OF POINTS IN ARRAYS.
!          *BLATB*   - INPUT LATITUDES.
!          *BLNGB*   - INPUT LONGITUDES.
!          *IJARB*   - OUTPUT SEA POINT NUMBERS.
!          *NPROC*     DIMENSION OF NSTART AND NEND
!          *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!          *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!          *IRANK*     PE INDEX


!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWGRID  , ONLY : IJS      ,IJL
      USE YOWMAP   , ONLY : BLK2GLO  ,KXLTMIN  ,KXLTMAX   ,IPER    ,   &
     &            AMOWEP   ,AMOSOP   ,AMONOP   ,XDELLA   ,XDELLO   ,   &
     &            ZDELLO   ,NLONRGG
      USE YOWPARAM , ONLY : NGY

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: NDIM, NBOUN, NPROC, IRANK

      INTEGER(KIND=JWIM), DIMENSION(NPROC), INTENT(IN) :: NSTART,NEND
      INTEGER(KIND=JWIM), DIMENSION(NDIM), INTENT(INOUT) :: IJARB

      REAL(KIND=JWRB), DIMENSION(NDIM),INTENT(IN) :: BLATB, BLNGB


      INTEGER(KIND=JWIM) :: IJ
      INTEGER(KIND=JWIM) :: ILATS, ILATN, IO, IOLT, IOLG

      REAL(KIND=JWRB) :: ZDEL, ALONG 
! ----------------------------------------------------------------------

!*    1. LOOP OVER INPUT LATITUDES, LONGITUDES.
!        --------------------------------------

      ILATS=KXLTMIN(IRANK)-1
      ILATN=KXLTMAX(IRANK)+1

      DO IO = 1,NBOUN

!*    1.1 COMPUTE GRID MATRIX INDICES.
!         ----------------------------

        IOLT = NINT((BLATB(IO)-AMOSOP)/XDELLA+1.0_JWRB)
        ALONG = MOD(BLNGB(IO)-AMOWEP+720.0_JWRB,360.0_JWRB)
        IF (IOLT < ILATS .OR. IOLT > ILATN) THEN
          ZDEL=XDELLO
        ELSE
          ZDEL=ZDELLO(IOLT)
        ENDIF
        IOLG = NINT(ALONG/ZDEL+1.0_JWRB)
        IF (IOLG == (NLONRGG(IOLT)+1) .AND. IPER == 1) IOLG = 1

        IF (IOLG < 1     .OR. IOLG > NLONRGG(IOLT).OR.                 &
     &      IOLT < ILATS .OR. IOLT > ILATN             ) THEN
          IJARB(IO) = 0
        ELSE
!*    1.2 SEARCH BLOCK NUMBER AND SEA POINT NUMBER.
!         -----------------------------------------
          IJARB(IO) = 0
          DO IJ = NSTART(IRANK),NEND(IRANK)
            IF (BLK2GLO(IJ)%IXLG == IOLG .AND. BLK2GLO(IJ)%KXLT == IOLT) THEN
              IJARB(IO) = IJ
              EXIT
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      END SUBROUTINE FINDB