SUBROUTINE INITDPTHFLDS(IJS, IJL, WVENVI, WVPRPT, WVPRPT_LAND)
! ----------------------------------------------------------------------

!**** *INITDPTHFLDS* - 

!*    PURPOSE.
!     --------

!     INITIALISE GRID POINT FIELDS DEPENDENT ON WATER DEPTH:
!     EMAXDPT

!     INITIALISE GRID POINT FIELDS DEPENDENT ON WATER DEPTH AND FREQUENCY:

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FREQUENCY

      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWFRED  , ONLY : ZPIFR
      USE YOWPARAM , ONLY : NFRE
      USE YOWSHAL  , ONLY : NDEPTH,                             &
     &                      TFAK, TCGOND, TSIHKD, TFAC_ST,      &
     &                      GAM_B_J

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      TYPE(ENVIRONMENT), DIMENSION(IJS:IJL), INTENT(INOUT) :: WVENVI
      TYPE(FREQUENCY), DIMENSION(IJS:IJL,NFRE), INTENT(OUT) :: WVPRPT
      TYPE(FREQUENCY), DIMENSION(NFRE), INTENT(OUT) :: WVPRPT_LAND


      INTEGER(KIND=JWIM) :: M, IJ, JKGLO, KIJS, KIJL, NPROMA

      REAL(KIND=JWRB) :: GAM
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INITDPTHFLDS',0,ZHOOK_HANDLE)

ASSOCIATE(DEPTH => WVENVI%DEPTH, &
 &        INDEP => WVENVI%INDEP, &
 &        EMAXDPT => WVENVI%EMAXDPT, &
 &        WAVNUM => WVPRPT%WAVNUM, &
 &        CINV => WVPRPT%CINV, &
 &        CGROUP => WVPRPT%CGROUP, &
 &        OMOSNH2KD => WVPRPT%OMOSNH2KD, &
 &        STOKFAC => WVPRPT%STOKFAC, &
 &        CIWA => WVPRPT%CIWA )


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
          EMAXDPT(IJ) = 0.0625_JWRB*(GAM*DEPTH(IJ))**2
        ENDDO

        DO M=1,NFRE
          DO IJ=KIJS,KIJL
            INDEP(IJ) = INDEP(IJ)
            WAVNUM(IJ,M) = TFAK(INDEP(IJ),M)
            CINV(IJ,M)   = WAVNUM(IJ,M)/ZPIFR(M)
            CGROUP(IJ,M) = TCGOND(INDEP(IJ),M)
            OMOSNH2KD(IJ,M) = TSIHKD(INDEP(IJ),M)
            STOKFAC(IJ,M) = TFAC_ST(INDEP(IJ),M)

            CIWA(IJ,M) = 1.0_JWRB
          ENDDO
        ENDDO

      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1502,1)


      ! Fictitious values for land point (NSUP+1)
      DO M=1,NFRE
        WVPRPT_LAND(M)%WAVNUM = TFAK(NDEPTH,M)
        WVPRPT_LAND(M)%CINV   = WVPRPT_LAND(M)%WAVNUM/ZPIFR(M)
        WVPRPT_LAND(M)%CGROUP = TCGOND(NDEPTH,M)
        WVPRPT_LAND(M)%OMOSNH2KD = TSIHKD(NDEPTH,M)
        WVPRPT_LAND(M)%STOKFAC = TFAC_ST(NDEPTH,M)

        WVPRPT_LAND(M)%CIWA = 1.0_JWRB
      ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('INITDPTHFLDS',1,ZHOOK_HANDLE)

END SUBROUTINE INITDPTHFLDS
