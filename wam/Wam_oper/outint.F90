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

!!1      CALL OUTGRID

!!!debile
        write(*,*) 'after OUTGRID'


!     CLEAN-UP
      IF(ALLOCATED(GOUT)) DEALLOCATE(GOUT)


! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('OUTINT',1,ZHOOK_HANDLE)

      END SUBROUTINE OUTINT
