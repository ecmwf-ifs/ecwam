      SUBROUTINE OUTINT

! ----------------------------------------------------------------------

!**** *OUTINT* - OUTPUT OF INTEGRATED FIELDS WITHOUT THE I/O SERVER

!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!       *CALL* *OUTINT

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
      USE YOWSTAT  , ONLY : CDATEA   ,CDATEF   ,CDTPRO   ,CDTINTT  ,    &
     &            CFDBSF   ,MARSTYPE ,NWFDBREF ,LFDBOPEN
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWUNIT  , ONLY : IU20     ,IU30
      USE FDBSUBS_MOD
      USE GRIB_API_INTERFACE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
  
      REAL(KIND=JWRB) :: ZHOOK_HANDLE


! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('OUTINT',0,ZHOOK_HANDLE)

      write(*,*) ' hello world !'

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('OUTINT',1,ZHOOK_HANDLE)

      END SUBROUTINE OUTINT
