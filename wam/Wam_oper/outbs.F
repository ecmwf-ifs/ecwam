      SUBROUTINE OUTBS (IJSLOC, IJLLOC, IJ_OFFSET, MIJ,                 &
     &                  FL1, XLLWS,                                     &
     &                  BOUTINFO, BOUT,                                 &
     &                  IG, IU25, IU26, LLOUTBS)
! ----------------------------------------------------------------------

!**** *OUTBS* - MODEL OUTPUT FROM BLOCK TO FILE, PRINTER AND COMMON.

!     H. GUNTHER         GKSS/ECMWF         JUNE 1990
!     J. BIDLOT          ECMWF   FEBRUARY 1996 MESSAGE PASSING
!     D. PETTENUZZO      MAY 2012 ADDED 3 TRAIN DECOMPOSITION CALL

!*    PURPOSE.
!     --------

!       CONTROL OUTPUT OF WAVE AND WIND FIELDS.


!**   INTERFACE.
!     ----------
!      *CALL*OUTBS (IJSLOC, IJLLOC, IJ_OFFSET,
!     &             MIJ, FL1, XLLWS,
!     &             IG, IU25, IU26, LLOUTBS)
!      *IJSLOC* - INDEX OF FIRST LOCAL GRIDPOINT
!      *IJLLOC* - INDEX OF LAST LOCAL GRIDPOINT
!      *IJ_OFFSET* OFFSET to point IJSLOC and IJLLOC to the global block of data
!                   only meaningful if unstructured grid
!      *MIJ*    - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!      *FL1*    - INPUT SPECTRUM.
!      *XLLWS*  - WINDSEA MASK FROM INPUT SOURCE TERM
!      *IG*     - BLOCK NUMBER
!      *IU25*   - OUTPUT UNIT FOR SPECTRA.
!      *IU26*   - OUTPUT UNIT FOR SEA AND SWELL SPECTRA.
!      *LLOUTBS - TRUE WHEN ACCESSING OUTBS

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
      USE YOWCOUT  , ONLY : NTRAIN   ,JPPFLAG  ,IPFGTBL  ,LSECONDORDER, &
     &            NIPRMOUT, ITOBOUT
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSHAL, ONLY: NDEPTH      ,DEPTH       ,INDEP    ,IODP  ,    &
     &            TCGOND
      USE YOWSTAT  , ONLY : CDATEA   ,CDTPRO   ,CDTSPT   , CDTSPS  ,    &
     &            MARSTYPE ,NPROMA_WAM
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJSLOC, IJLLOC
      INTEGER, INTENT(IN) :: IJ_OFFSET
      INTEGER, INTENT(IN) :: IG
      INTEGER, INTENT(IN) :: IU25, IU26 
      INTEGER, DIMENSION(IJSLOC:IJLLOC), INTENT(IN) :: MIJ

      REAL, DIMENSION(JPPFLAG), INTENT(OUT) :: BOUTINFO
      REAL, DIMENSION(IJSLOC:IJLLOC,NANG,NFRE), INTENT(IN) :: FL1
      REAL, DIMENSION(IJSLOC:IJLLOC,NANG,NFRE), INTENT(IN) :: XLLWS 
      REAL, DIMENSION(IJSLOC:IJLLOC,JPPFLAG), INTENT(OUT) :: BOUT

      LOGICAL, INTENT(OUT) :: LLOUTBS

      INTEGER :: ITT, ITG, ICT, IR, ITR
      INTEGER :: IJ, M 
      INTEGER :: JKGLO, KIJS, KIJL, NPROMA

      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION(IJSLOC:IJLLOC,NFRE) :: CGGROUP 

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

      LLOUTBS=.TRUE.

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


!*    4. OUTPUT OF SPECTRA AT SELECTED GRID POINTS.
!        ------------------------------------------

      IF (CDTSPT.EQ.CDTPRO .OR. CDTSPS.EQ.CDTPRO) THEN
        CALL OUTSPP (FL1, IJSLOC, IJLLOC, IJ_OFFSET, IG, IU25, IU26)
        IF (ITEST.GE.3) THEN
            WRITE(IU06,*) '      SUB. OUTBS: OUTPUT OF SPECTRA',        &
     &       ' AT SELECTED POINTS DONE'
        ENDIF
      ENDIF


!*    6. OUTPUT OF SPECTRA FOR SATELLITE COLLOCATION.
!        --------------------------------------------

      IF (MARSTYPE.EQ.'an'.OR.MARSTYPE.EQ.'fg'.OR.MARSTYPE.EQ.'4v') THEN
        IF (CDTPRO.NE.CDATEA) THEN
          CALL OUTERS (FL1, IJSLOC, IJLLOC, CDTPRO)
          IF (ITEST.GE.3) THEN
            IF (ITESTB.GE.IG) THEN
              WRITE(IU06,*) '      SUB. OUTBS: OUTPUT OF SPECTRA',      &
     &         ' FOR SATELLITE COLLOCATION DONE FOR ', CDTPRO
            ENDIF
          ENDIF
        ENDIF
      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTBS',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE OUTBS
