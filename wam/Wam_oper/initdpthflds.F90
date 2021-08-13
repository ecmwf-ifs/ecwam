SUBROUTINE INITDPTHFLDS
! ----------------------------------------------------------------------

!**** *INITDPTHFLDS* - 

!*    PURPOSE.
!     --------

!     INITIALISE GRID POINT FIELDS DEPENDENT ON WATER DEPTH:
!     EMAXDPT

!     INITIALISE GRID POINT FIELDS DEPENDENT ON WATER DEPTH AND FREQUENCY:
!     WAVNUM, CINV, CGROUP

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : ZPIFR
      USE YOWGRID  , ONLY : IJS, IJL
      USE YOWPARAM , ONLY : NFRE
      USE YOWSHAL  , ONLY : DEPTH, NDEPTH,                              &
     &                      INDEP, TFAK, TCGOND, TSIHKD, TFAC_ST,       &
     &                      GAM_B_J, EMAXDPT,                           &
     &                      WAVNUM, CINV, CGROUP, OMOSNH2KD, STOKFAC,   &
     &                      WAVNUM_LAND, CGROUP_LAND, OMOSNH2KD_LAND
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

      IF (.NOT.ALLOCATED(WAVNUM)) ALLOCATE(WAVNUM(IJS:IJL,NFRE))
      IF (.NOT.ALLOCATED(WAVNUM_LAND)) ALLOCATE(WAVNUM_LAND(NFRE))
      IF (.NOT.ALLOCATED(CINV)) ALLOCATE(CINV(IJS:IJL,NFRE))
      IF (.NOT.ALLOCATED(CGROUP)) ALLOCATE(CGROUP(IJS:IJL,NFRE))
      IF (.NOT.ALLOCATED(CGROUP_LAND)) ALLOCATE(CGROUP_LAND(NFRE))
      IF (.NOT.ALLOCATED(OMOSNH2KD)) ALLOCATE(OMOSNH2KD(IJS:IJL,NFRE))
      IF (.NOT.ALLOCATED(OMOSNH2KD_LAND)) ALLOCATE(OMOSNH2KD_LAND(NFRE))
      IF (.NOT.ALLOCATED(STOKFAC)) ALLOCATE(STOKFAC(IJS:IJL,NFRE))


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
            WAVNUM(IJ,M) = TFAK(INDEP(IJ),M) 
            CINV(IJ,M)   = WAVNUM(IJ,M)/ZPIFR(M) 
            CGROUP(IJ,M) = TCGOND(INDEP(IJ),M)
            OMOSNH2KD(IJ,M) = TSIHKD(INDEP(IJ),M) 
            STOKFAC(IJ,M) = TFAC_ST(INDEP(IJ),M)
          ENDDO
        ENDDO

        ! Fictitious values for land point (NSUP+1)
        DO M=1,NFRE
          WAVNUM_LAND(M) = TFAK(NDEPTH,M) 
          CGROUP_LAND(M) = TCGOND(NDEPTH,M) 
          OMOSNH2KD_LAND(M) = TSIHKD(NDEPTH,M) 
        ENDDO

      ENDDO
!$OMP END PARALLEL DO

      CALL GSTATS(1502,1)

IF (LHOOK) CALL DR_HOOK('INITDPTHFLDS',1,ZHOOK_HANDLE)

END SUBROUTINE INITDPTHFLDS
