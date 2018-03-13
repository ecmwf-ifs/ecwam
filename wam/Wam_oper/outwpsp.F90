      SUBROUTINE OUTWPSP (IJSLOC, IJLLOC, IJ_OFFSET, FL1, IG, IU25, IU26)
! ----------------------------------------------------------------------

!**** *OUTWPSP* - MODEL OUTPUT OF SPECTRA AT GIVEN LOCATIONS 

!*    PURPOSE.
!     --------

!       CONTROL OUTPUT OF WAVE AND WIND FIELDS.


!**   INTERFACE.
!     ----------
!      *CALL*OUTWPSP (IJSLOC, IJLLOC, IJ_OFFSET, FL1, IG, IU25, IU26)
!      *IJSLOC* - INDEX OF FIRST LOCAL GRIDPOINT
!      *IJLLOC* - INDEX OF LAST LOCAL GRIDPOINT
!      *IJ_OFFSET* OFFSET to point IJSLOC and IJLLOC to the global block of data
!                   only meaningful if unstructured grid
!      *FL1*    - INPUT SPECTRUM.
!      *IG*     - BLOCK NUMBER
!      *IU25*   - OUTPUT UNIT FOR SPECTRA.
!      *IU26*   - OUTPUT UNIT FOR SEA AND SWELL SPECTRA.

!     EXTERNALS.
!     ----------

!       *OUTERS*    - OUTPUT OF SATELLITE COLOCATION SPECTRA.
!       *OUTSPP*    - OUTPUT OF SPECTRA AT SELECTED POINTS.
!   
!     METHOD.
!     -------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSTAT  , ONLY : CDATEA   ,CDTPRO   ,CDTSPT   , CDTSPS  , MARSTYPE
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJSLOC, IJLLOC
      INTEGER(KIND=JWIM), INTENT(IN) :: IJ_OFFSET
      INTEGER(KIND=JWIM), INTENT(IN) :: IG
      INTEGER(KIND=JWIM), INTENT(IN) :: IU25, IU26 
      REAL(KIND=JWRB), DIMENSION(IJSLOC:IJLLOC,NANG,NFRE), INTENT(IN) :: FL1

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTWPSP',0,ZHOOK_HANDLE)
#endif

!*    1. OUTPUT OF SPECTRA AT SELECTED GRID POINTS.
!        ------------------------------------------

      IF (CDTSPT.EQ.CDTPRO .OR. CDTSPS.EQ.CDTPRO) THEN
        CALL OUTSPP (FL1, IJSLOC, IJLLOC, IJ_OFFSET, IG, IU25, IU26)
        IF (ITEST.GE.3) THEN
            WRITE(IU06,*) '      SUB. OUTWPSP: OUTPUT OF SPECTRA',      &
     &       ' AT SELECTED POINTS DONE'
        ENDIF
      ENDIF


!*    2. OUTPUT OF SPECTRA FOR SATELLITE COLLOCATION.
!        --------------------------------------------

      IF (MARSTYPE.EQ.'an'.OR.MARSTYPE.EQ.'fg'.OR.MARSTYPE.EQ.'4v') THEN
        IF (CDTPRO.NE.CDATEA) THEN
          CALL OUTERS (FL1, IJSLOC, IJLLOC, CDTPRO)
          IF (ITEST.GE.3) THEN
            WRITE(IU06,*) '      SUB. OUTWPSP: OUTPUT OF SPECTRA',      &
     &       ' FOR SATELLITE COLLOCATION DONE FOR ', CDTPRO
          ENDIF
        ENDIF
      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTWPSP',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE OUTWPSP
