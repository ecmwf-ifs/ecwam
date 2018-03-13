      SUBROUTINE OUTBS (IJSLOC, IJLLOC, MIJ, IG, FL1, XLLWS, BOUT)
! ----------------------------------------------------------------------

!**** *OUTBS* - MODEL OUTPUT FROM BLOCK TO FILE, PRINTER AND COMMON.

!*    PURPOSE.
!     --------

!       CONTROL OUTPUT OF WAVE AND WIND FIELDS (except spectrum).

!**   INTERFACE.
!     ----------
!      *CALL*OUTBS (IJSLOC, IJLLOC, MIJ, IG, FL1, XLLWS, BOUT) 
!      *IJSLOC* - INDEX OF FIRST LOCAL GRIDPOINT
!      *IJLLOC* - INDEX OF LAST LOCAL GRIDPOINT
!      *IJ_OFFSET* OFFSET to point IJSLOC and IJLLOC to the global block of data
!                   only meaningful if unstructured grid
!      *MIJ*    - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!      *IG*     - BLOCK NUMBER
!      *FL1*    - INPUT SPECTRUM.
!      *XLLWS*  - WINDSEA MASK FROM INPUT SOURCE TERM
!      *BOUT*   - BLOCK OF SELECTED OUTPUT PARAMETERS.

!     EXTERNALS.
!     ----------

!       *OUTERS*    - OUTPUT OF SATELLITE COLOCATION SPECTRA.
!       *OUTSPP*    - OUTPUT OF SPECTRA AT SELECTED POINTS.
!       *OUTBLOCK*  - GET ALL OUTPUT PARAMETERS
!   
!     METHOD.
!     -------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : JPPFLAG  ,NIPRMOUT
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSHAL  , ONLY : DEPTH       ,INDEP    , TCGOND
      USE YOWSTAT  , ONLY : NPROMA_WAM
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJSLOC, IJLLOC
      INTEGER(KIND=JWIM), INTENT(IN) :: IG
      INTEGER(KIND=JWIM), DIMENSION(IJSLOC:IJLLOC), INTENT(IN) :: MIJ
           
      REAL(KIND=JWRB), DIMENSION(IJSLOC:IJLLOC,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(IJSLOC:IJLLOC,NANG,NFRE), INTENT(IN) :: XLLWS 
      REAL(KIND=JWRB), DIMENSION(IJSLOC:IJLLOC,NIPRMOUT), INTENT(OUT) :: BOUT

      INTEGER(KIND=JWIM) :: M, IJ, JKGLO, KIJS, KIJL, NPROMA

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJSLOC:IJLLOC,NFRE) :: CGGROUP 

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTBS',0,ZHOOK_HANDLE)
#endif
!
!  cgroup should be passed and used elsewhere
!  to fix later !!!
      DO M=1,NFRE
        DO IJ=IJSLOC,IJLLOC
          CGROUP(IJ,M)=TCGOND(INDEP(IJ),M)
        ENDDO
      ENDDO

!*    1. COMPUTE MEAN PARAMETERS.
!        ------------------------

      NPROMA=NPROMA_WAM

!     COMPUTE MEAN PARAMETERS

      CALL GSTATS(1502,0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
      DO JKGLO=IJSLOC,IJLLOC,NPROMA
        KIJS=JKGLO
        KIJL=MIN(KIJS+NPROMA-1,IJLLOC)
        CALL OUTBLOCK(KIJS, KIJL, MIJ(KIJS), IG,                        &
     &                FL1(KIJS:KIJL,:,:), XLLWS(KIJS:KIJL,:,:),         &
     &                DEPTH(KIJS:KIJL,IG),CGROUP(KIJS:KIJL,:),          &
     &                BOUT(KIJS,:))
      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1502,1)

      IF (ITEST.GE.3) THEN
          WRITE(IU06,*) '      SUB. OUTBS: INTEGRATED',                 &
     &     ' PARAMETERS COMPUTED FOR OUTPUT'
      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTBS',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE OUTBS
