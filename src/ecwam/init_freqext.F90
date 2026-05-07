! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE INIT_FREQEXT

! ----------------------------------------------------------------------

!**** *INIT_FREQEXT* -

!*    PURPOSE.
!     ---------

!     INITIALISATION FOR EXTENDED FREQUENCY-SPACE


!**   INTERFACE.
!     ----------

!       *CALL* *INIT_FREQEXT*

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB

      USE YOWFRED  , ONLY : DELTH    ,FR       ,DFIM     ,FRATIO   ,    &
     &                      SIG      ,DSII     ,SIGM1    ,DF       ,    &
     &                      SIG_EXT  ,DSII_EXT ,IFRE_EXT ,NFRE_EXT ,DDEN
      USE YOWPARAM , ONLY : NFRE
      USE YOWPHYS  , ONLY : FRQMAX
      USE YOWPCONS , ONLY : ZPI

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: M

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('INIT_FREQEXT',0,ZHOOK_HANDLE)

      IF (.NOT.ALLOCATED(DF))    ALLOCATE(DF(NFRE))
      IF (.NOT.ALLOCATED(SIG))   ALLOCATE(SIG(NFRE))
      IF (.NOT.ALLOCATED(DDEN))  ALLOCATE(DDEN(NFRE))
      IF (.NOT.ALLOCATED(DSII))  ALLOCATE(DSII(NFRE))
      IF (.NOT.ALLOCATED(SIGM1)) ALLOCATE(SIGM1(NFRE))
      DO M=1,NFRE
            DF(M)    = DFIM(M)/DELTH
            SIG(M)   = ZPI*FR(M)
            DSII(M)  = ZPI*DF(M)
            DDEN(M)  = ZPI*DFIM(M)*SIG(M)
            SIGM1(M) = 1.0_JWRB/SIG(M)
      ENDDO

      ! DETERMINE THE NUMBER OF FREQUENCIES TO EXTEND TO
      NFRE_EXT = CEILING(LOG(FRQMAX/FR(1))/LOG(FRATIO))+1
      NFRE_EXT = MAX(NFRE,NFRE_EXT)
      IF (ALLOCATED(IFRE_EXT)) THEN
            IF (SIZE(IFRE_EXT) /= NFRE_EXT) DEALLOCATE(IFRE_EXT)
      END IF
      IF (.NOT.ALLOCATED(IFRE_EXT)) ALLOCATE(IFRE_EXT(NFRE_EXT))
      IFRE_EXT = (/ (REAL(M, KIND=JWRB), M=1,NFRE_EXT) /)
      IF (.NOT.ALLOCATED(SIG_EXT))  ALLOCATE(SIG_EXT(NFRE_EXT))
      IF (.NOT.ALLOCATED(DSII_EXT)) ALLOCATE(DSII_EXT(NFRE_EXT))
      IF (NFRE .LT. NFRE_EXT) THEN
            SIG_EXT                   = SIG(1)*FRATIO**(IFRE_EXT-1.0_JWRB)
            DSII_EXT                  = 0.5_JWRB * SIG_EXT * (FRATIO-1.0_JWRB/FRATIO)
            ! The first and last frequency bin:
            DSII_EXT(1)               = 0.5_JWRB * SIG_EXT(1) * (FRATIO-1.0_JWRB)
            DSII_EXT(NFRE_EXT)        = 0.5_JWRB * SIG_EXT(NFRE_EXT) * (FRATIO-1.0_JWRB) / FRATIO
      ELSE
            SIG_EXT        = SIG
            DSII_EXT       = DSII
      END IF


      IF (LHOOK) CALL DR_HOOK('INIT_FREQEXT',1,ZHOOK_HANDLE)

      END SUBROUTINE INIT_FREQEXT
