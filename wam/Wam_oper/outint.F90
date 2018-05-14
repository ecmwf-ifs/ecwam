      SUBROUTINE OUTINT(CDATE, IFCST)

! ----------------------------------------------------------------------

!**** *OUTINT* - OUTPUT OF INTEGRATED FIELDS WITHOUT THE I/O SERVER

!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!       *CALL* *OUTINT(CDATE, IFCST)
!          *CDATE*   DATE AND TIME.
!          *IFCST*   FORECAST STEP IN HOURS.

!     EXTERNALS.
!     ----------
!       *OUTGRID*
!       *WGRIBENOUT*

!     METHOD.
!     -------

!       NONE.

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : FFLAG    ,FFLAG20  ,GFLAG    ,GFLAG20  ,    &
     &            JPPFLAG  ,LFDB     ,IPFGTBL  ,ITOBOUT  ,INFOBOUT ,    &
     &            LOUTINT
      USE YOWGRID  , ONLY : NLONRGG  ,DELPHI
      USE YOWINTP  , ONLY : GOUT
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,    &
     &            AMONOP   ,ZDELLO
      USE YOWMPP   , ONLY : IRANK
      USE YOWPARAM , ONLY : NGX      ,NGY      ,CLDOMAIN
      USE YOWSTAT  , ONLY : CDTPRO   ,CDTINTT  ,                        &
     &            CFDBSF   ,MARSTYPE ,NWFDBREF ,LFDBOPEN
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWUNIT  , ONLY : IU20     ,IU30

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE FDBSUBS_MOD
      USE GRIB_API_INTERFACE

! ----------------------------------------------------------------------
      IMPLICIT NONE

      CHARACTER(LEN=14), INTENT(IN) :: CDATE
      INTEGER(KIND=JWIM), INTENT(IN) :: IFCST

      INTEGER(KIND=JWIM) :: I, J, ITG, IFLAG, IT
      INTEGER(KIND=JWIM) :: IGLOBAL, ILOCAL      ! FDB field counters
      INTEGER(KIND=JWIM) :: IPARAM, ITABLE, IZLEV
      INTEGER(KIND=JWIM) :: IERR
      INTEGER(KIND=JWIM) :: IUOUT 
  
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      CHARACTER(LEN=24) :: OFILENAME

      LOGICAL, SAVE :: FRSTIME30

      DATA FRSTIME30 / .TRUE. /

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('OUTINT',0,ZHOOK_HANDLE)

      LOUTINT=.TRUE.

!     1. COLLECT INTEGRATED PARAMETERS FOR OUTPUT ON SELECTED PE's 
!        i.e BOUT => GOUT
!        ---------------------------------------------------------
!!!debile
        write(*,*) 'before OUTGRID'

      CALL OUTGRID

!!!debile
        write(*,*) 'after OUTGRID'


!     2. ONE GRID POINT OUTPUT
!        ---------------------
      IF(CLDOMAIN.EQ.'s') CALL OUT_ONEGRDPT(IU06)


!*    3. OUTPUT IN PURE BINARY FORM
!        --------------------------
      IF(FFLAG20 .AND. CDTINTT.EQ.CDTPRO ) THEN
        DO IFLAG=1,JPPFLAG
          IF (FFLAG(IFLAG) .AND. IRANK.EQ.IPFGTBL(IFLAG)) THEN
            WRITE (IU20) CDTPRO, NGX, NGY, IRGG
            WRITE (IU20) AMOWEP,AMOSOP,AMOEAP,AMONOP
            WRITE (IU20) ZDELLO, NLONRGG, DELPHI
            WRITE (IU20) GOUT(IFLAG,:,:)
          ENDIF
        ENDDO
      ENDIF

!*    4.  PACK INTEGRATED PARAMETERS INTO GRIB AND OUTPUT
!         -----------------------------------------------
      IF (GFLAG20 .AND. CDTINTT.EQ.CDTPRO ) THEN

!       DEFINE OUTPUT FILE FOR GRIB DATA (if not written to FDB)
        IF(.NOT.LFDB) THEN
            WRITE(OFILENAME,'(A,I2)') 'fort.',IU30
            IF (FRSTIME30) THEN
              CALL IGRIB_OPEN_FILE(IUOUT,OFILENAME,'w')
              FRSTIME30=.FALSE.
            ELSE
              CALL IGRIB_OPEN_FILE(IUOUT,OFILENAME,'a')
            ENDIF
        ELSE
          IUOUT=0
        ENDIF

!       OUTPUT:
        IGLOBAL=0
        ILOCAL=0
        DO IFLAG=1,JPPFLAG
          IF (GFLAG(IFLAG)) THEN
            IGLOBAL=IGLOBAL+1
            IF (IRANK.EQ.IPFGTBL(IFLAG)) THEN
              ILOCAL=ILOCAL+1
              IF(ILOCAL.GT. SIZE(GOUT,1)) THEN
                WRITE(*,*) ' -------------------------------------'
                WRITE(*,*) ' ERROR in OUTINT '
                WRITE(*,*) ' ACCESSING MORE FIELDS THAN AVAILABLE'
                WRITE(*,*) ' SIZE(GOUT,1) = ',SIZE(GOUT,1) 
                WRITE(*,*) ' -------------------------------------'
                CALL ABORT1
              ENDIF

              IT=ITOBOUT(IFLAG)
              IPARAM=INFOBOUT(IT,2)
              ITABLE=INFOBOUT(IT,1)
              IZLEV=INFOBOUT(IT,3)

              CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, GOUT(ILOCAL,:,:),  &
     &                        ITABLE, IPARAM, IZLEV, 0 , 0,             &
     &                        CDATE, IFCST, MARSTYPE,                   &
     &                        LFDB, CFDBSF, NWFDBREF, LFDBOPEN, IUOUT)
            ENDIF
          ENDIF
        ENDDO

        IF(LFDB.AND.ILOCAL.GT.0) THEN
          IERR = ISETFIELDCOUNTFDBSUBS(NWFDBREF,IGLOBAL,ILOCAL)
          IF(IERR.NE.0)THEN
            WRITE(IU06,*) ' ------------------------'
            WRITE(IU06,*) ' ERROR setting fdb field count '
            WRITE(IU06,*) ' in routine OUTINT '
            WRITE(IU06,*) ' FDB ERROR CODE IS ',IERR
            WRITE(IU06,*) ' ------------------------'
            CALL ABORT1
          ENDIF
        ELSEIF(.NOT.LFDB) THEN
          CALL IGRIB_CLOSE_FILE(IUOUT)
        ENDIF

      ENDIF ! end grib output

!     CLEAN-UP
      IF(ALLOCATED(GOUT)) DEALLOCATE(GOUT)


! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('OUTINT',1,ZHOOK_HANDLE)

      END SUBROUTINE OUTINT
