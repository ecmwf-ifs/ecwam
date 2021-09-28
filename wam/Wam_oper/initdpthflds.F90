SUBROUTINE INITDPTHFLDS
! ----------------------------------------------------------------------

!**** *INITDPTHFLDS* - 

!*    PURPOSE.
!     --------

!     INITIALISE GRID POINT FIELDS DEPENDENT ON WATER DEPTH:
!     EMAXDPT

!     INITIALISE GRID POINT FIELDS DEPENDENT ON WATER DEPTH AND FREQUENCY:

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE,   ONLY : FREQUENCY

      USE YOWFRED  , ONLY : ZPIFR
      USE YOWGRID  , ONLY : IJS, IJL
      USE YOWPARAM , ONLY : NFRE
      USE YOWSHAL  , ONLY : DEPTH, NDEPTH,                              &
     &                      INDEP, TFAK, TCGOND, TSIHKD, TFAC_ST,       &
     &                      GAM_B_J, EMAXDPT,                           &
     &                      WVPRPT, WVPRPT_LAND
      USE YOWSTAT  , ONLY : NPROMA_WAM

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: M, IJ, JKGLO, KIJS, KIJL, NPROMA

      REAL(KIND=JWRB) :: GAM
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INITDPTHFLDS',0,ZHOOK_HANDLE)

      IF (.NOT.ALLOCATED(EMAXDPT)) ALLOCATE(EMAXDPT(IJS:IJL))

      IF (.NOT.ALLOCATED(WVPRPT)) ALLOCATE(WVPRPT(IJS:IJL,NFRE))
      IF (.NOT.ALLOCATED(WVPRPT_LAND)) ALLOCATE(WVPRPT_LAND(NFRE))


      CALL GSTATS(1502,0)
      NPROMA=NPROMA_WAM

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL,M,IJ,GAM)
      DO JKGLO = IJS, IJL, NPROMA
        KIJS=JKGLO
        KIJL=MIN(KIJS+NPROMA-1,IJL)

!       COMPUTE THE MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH
        DO IJ = KIJS, KIJL
!         REDUCE GAMMA FOR SMALL DEPTH ( < 4m)
!         (might need to be revisted when grid is fine resolution)
          IF (DEPTH(IJ) < 4.0_JWRB) THEN
            GAM=GAM_B_J*DEPTH(IJ)/4.0_JWRB
          ELSE
            GAM=GAM_B_J
          ENDIF
          EMAXDPT(IJ)=0.0625_JWRB*(GAM*DEPTH(IJ))**2
        ENDDO

        DO M=1,NFRE
          DO IJ=KIJS,KIJL
            WVPRPT(IJ,M)%WAVNUM = TFAK(INDEP(IJ),M)
            WVPRPT(IJ,M)%CINV   = WVPRPT(IJ,M)%WAVNUM/ZPIFR(M)
            WVPRPT(IJ,M)%CGROUP = TCGOND(INDEP(IJ),M)
            WVPRPT(IJ,M)%OMOSNH2KD = TSIHKD(INDEP(IJ),M)
            WVPRPT(IJ,M)%STOKFAC = TFAC_ST(INDEP(IJ),M)

            WVPRPT(IJ,M)%CIWA = 1.0_JWRB
          ENDDO
        ENDDO

        ! Fictitious values for land point (NSUP+1)
        DO M=1,NFRE
          WVPRPT_LAND(M)%WAVNUM = TFAK(NDEPTH,M)
          WVPRPT_LAND(M)%CINV   = WVPRPT_LAND(M)%WAVNUM/ZPIFR(M)
          WVPRPT_LAND(M)%CGROUP = TCGOND(NDEPTH,M)
          WVPRPT_LAND(M)%OMOSNH2KD = TSIHKD(NDEPTH,M)
          WVPRPT_LAND(M)%STOKFAC = TFAC_ST(NDEPTH,M)

          WVPRPT_LAND(M)%CIWA = 1.0_JWRB
        ENDDO

      ENDDO
!$OMP END PARALLEL DO

      CALL GSTATS(1502,1)

IF (LHOOK) CALL DR_HOOK('INITDPTHFLDS',1,ZHOOK_HANDLE)

END SUBROUTINE INITDPTHFLDS
