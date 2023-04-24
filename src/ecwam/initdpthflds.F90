! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE INITDPTHFLDS(WVENVI, WVPRPT, WVPRPT_LAND)
! ----------------------------------------------------------------------

!**** *INITDPTHFLDS* - 

!*    PURPOSE.
!     --------

!     INITIALISE GRID POINT FIELDS DEPENDENT ON WATER DEPTH:
!     EMAXDPT

!     INITIALISE GRID POINT FIELDS DEPENDENT ON WATER DEPTH AND FREQUENCY:

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FREQUENCY, FREQUENCY_LAND

      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWFRED  , ONLY : ZPIFR
      USE YOWPARAM , ONLY : NFRE
      USE YOWSHAL  , ONLY : NDEPTH,                                &
     &                      TFAK, TCGOND, TSIHKD, TFAC_ST, GAM_B_J

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      TYPE(ENVIRONMENT), INTENT(INOUT) :: WVENVI
      TYPE(FREQUENCY), INTENT(INOUT) :: WVPRPT
      TYPE(FREQUENCY_LAND), INTENT(INOUT) :: WVPRPT_LAND


      INTEGER(KIND=JWIM) :: M, IJ, ICHNK, KIJS, KIJL

      REAL(KIND=JWRB) :: GAM
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INITDPTHFLDS',0,ZHOOK_HANDLE)


      CALL GSTATS(1504,0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ICHNK, KIJS, KIJL, IJ, M, GAM)
      DO ICHNK = 1, NCHNK
        KIJS = 1
        KIJL = NPROMA_WAM

!       COMPUTE THE MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH
        DO IJ = KIJS, KIJL
!         REDUCE GAMMA FOR SMALL DEPTH ( < 4m)
!         (might need to be revisted when grid is fine resolution)
          IF (WVENVI%DEPTH(IJ,ICHNK) < 4.0_JWRB) THEN
            GAM=GAM_B_J*WVENVI%DEPTH(IJ,ICHNK)/4.0_JWRB
          ELSE
            GAM=GAM_B_J
          ENDIF
          WVENVI%EMAXDPT(IJ,ICHNK) = 0.0625_JWRB*(GAM*WVENVI%DEPTH(IJ,ICHNK))**2
        ENDDO

        DO M= 1, NFRE
          DO IJ = KIJS, KIJL
            WVPRPT%WAVNUM(IJ, M, ICHNK) = TFAK(WVENVI%INDEP(IJ,ICHNK), M)
            WVPRPT%CINV(IJ, M, ICHNK)   = WVPRPT%WAVNUM(IJ, M, ICHNK)/ZPIFR(M)
            WVPRPT%CGROUP(IJ, M, ICHNK) = TCGOND(WVENVI%INDEP(IJ,ICHNK), M)
            WVPRPT%XK2CG(IJ, M, ICHNK) = WVPRPT%WAVNUM(IJ, M, ICHNK)**2 * WVPRPT%CGROUP(IJ, M, ICHNK)
            WVPRPT%OMOSNH2KD(IJ, M, ICHNK) = TSIHKD(WVENVI%INDEP(IJ,ICHNK), M)
            WVPRPT%STOKFAC(IJ, M, ICHNK) = TFAC_ST(WVENVI%INDEP(IJ,ICHNK), M)

            WVPRPT%CIWA(IJ,M, ICHNK) = 1.0_JWRB
          ENDDO
        ENDDO

      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1504,1)


      ! Fictitious values for land point (NSUP+1)
      DO M=1,NFRE
        WVPRPT_LAND%WAVNUM(M) = TFAK(NDEPTH,M)
        WVPRPT_LAND%CINV(M)   = WVPRPT_LAND%WAVNUM(M)/ZPIFR(M)
        WVPRPT_LAND%CGROUP(M) = TCGOND(NDEPTH,M)
        WVPRPT_LAND%XK2CG(M) = WVPRPT_LAND%WAVNUM(M)**2 * WVPRPT_LAND%CGROUP(M)
        WVPRPT_LAND%OMOSNH2KD(M) = TSIHKD(NDEPTH,M)
        WVPRPT_LAND%STOKFAC(M) = TFAC_ST(NDEPTH,M)

        WVPRPT_LAND%CIWA(M) = 1.0_JWRB
      ENDDO

IF (LHOOK) CALL DR_HOOK('INITDPTHFLDS',1,ZHOOK_HANDLE)

END SUBROUTINE INITDPTHFLDS
