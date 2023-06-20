! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
SUBROUTINE WAMINTGR_CUDA_GPU(CDTPRA, CDATE, CDATEWH, CDTIMP, CDTIMPNEXT, BLK2GLO, WVENVI, WVPRPT, FF_NOW, FF_NEXT, INTFLDS,  &
& WAM2NEMO, MIJ, FL1, XLLWS, TIME1)
  
  
  ! ----------------------------------------------------------------------
  
  !**** *WAMINTGR* - 3-G WAM MODEL - TIME INTEGRATION OF WAVE FIELDS.
  
  !*    PURPOSE.
  !     --------
  
  !       COMPUTATION OF THE 2-D FREQUENCY-DIRECTION WAVE SPECTRUM AT ALL
  !       GRID POINTS FOR A GIVEN INITIAL SPECTRUM AND FORCING SURFACE
  !       STRESS FIELD.
  
  !     REFERENCE.
  !     ----------
  
  !         IFS DOCUMENTATION, part VII
  
  ! -------------------------------------------------------------------
  
  USE IMPLSCH_CUF_PARAMETRISE_MOD, ONLY: IMPLSCH_CUF_PARAMETRISE
  USE ISO_C_BINDING, ONLY: C_SIZEOF
  USE OPENACC
  USE CUDAFOR
  USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
  USE YOWDRVTYPE, ONLY: WVGRIDGLO, ENVIRONMENT, FREQUENCY, FORCING_FIELDS, INTGT_PARAM_FIELDS, WAVE2OCEAN
  
  USE YOWCOUP, ONLY: LWNEMOCOU, NEMONTAU
  USE YOWGRID, ONLY: NPROMA_WAM, NCHNK
  USE YOWPARAM, ONLY: NIBLO, NANG, NFRE
  USE YOWPCONS, ONLY: EPSMIN
  USE YOWINDN, ONLY: MLSTHG
  USE YOWSTAT, ONLY: CDTPRO, IDELPRO, IDELT, IDELWI, LLSOURCE
  USE YOWWIND, ONLY: CDAWIFL, CDATEWO, CDATEFL
  USE YOWFIELD_MOD, ONLY: FREQUENCY_FIELD, ENVIRONMENT_FIELD, FORCING_FIELDS_FIELD, WAVE2OCEAN_FIELD, INTGT_PARAM_FIELDS_FIELD,  &
  & SOURCE_CONTRIBS_FIELD
  ![Loki::GlobalVarOffload].....Adding global variables to driver symbol table for offload instructions
  USE YOWFRED, ONLY: fr, dfimfr2, flogsprdm1, xkm_gc, omega_gc, delkcc_gc_ns, flmax, c2osqrtvg_gc, sinth, cofrm4, dfimfr,  &
  & dfim_sim, rhowg_dfim, costh, delth, xkmsqrtvgoc2_gc, xk_gc, dfimofr, omxkm3_gc, th, zpifr, dfim, delkcc_omxkm3_gc, fr5,  &
  & cm_gc, om3gmkm_gc, nwav_gc, fr_d, dfimfr2_d, flogsprdm1_d, xkm_gc_d, omega_gc_d, delkcc_gc_ns_d, flmax_d, c2osqrtvg_gc_d,  &
  & sinth_d, cofrm4_d, dfimfr_d, dfim_sim_d, rhowg_dfim_d, costh_d, delth_d, xkmsqrtvgoc2_gc_d, xk_gc_d, dfimofr_d, omxkm3_gc_d,  &
  & th_d, zpifr_d, dfim_d, delkcc_omxkm3_gc_d, fr5_d, cm_gc_d, om3gmkm_gc_d, nwav_gc_d
  USE YOWINDN, ONLY: k11w, dal1, k1w, k2w, mfrstlw, ikm1, rnlcoef, fklap, dal2, fklam1, ikp, inlcoef, k21w, fklam, af11, ikp1,  &
  & fklap1, ikm, kfrh, k11w_d, dal1_d, k1w_d, k2w_d, mfrstlw_d, ikm1_d, rnlcoef_d, fklap_d, dal2_d, fklam1_d, ikp_d, inlcoef_d,  &
  & k21w_d, fklam_d, af11_d, ikp1_d, fklap1_d, ikm_d, kfrh_d, mlsthg_d
  USE YOWPHYS, ONLY: dthrn_a, ndikcumul, rnu, bmaxokap, chnkmin_u, tailfactor, indicessat, rn1_rn, satweights, alpha, z0rat,  &
  & tailfactor_pm, cdisvis, cumulw, alphapmax, dthrn_u, nsdsnth, delta_sdis, cdis, tauwshelter, alphamin, z0tubmax, swellf5,  &
  & ang_gc_c, zalp, ang_gc_a, gamnconst, betamaxoxkappa2, ang_gc_b, rnum, dthrn_a_d, ndikcumul_d, rnu_d, bmaxokap_d,  &
  & chnkmin_u_d, tailfactor_d, indicessat_d, rn1_rn_d, satweights_d, alpha_d, z0rat_d, tailfactor_pm_d, cdisvis_d, cumulw_d,  &
  & alphapmax_d, dthrn_u_d, nsdsnth_d, delta_sdis_d, cdis_d, tauwshelter_d, alphamin_d, z0tubmax_d, swellf5_d, ang_gc_c_d,  &
  & zalp_d, ang_gc_a_d, gamnconst_d, betamaxoxkappa2_d, ang_gc_b_d, rnum_d
  USE YOWPCONS, ONLY: zpi, sqrtgosurft, zpi4gm1, zpi4gm2, gm1, g, zpi_d, sqrtgosurft_d, zpi4gm1_d, zpi4gm2_d, gm1_d, g_d
  USE YOWCOUP, ONLY: lwnemocousend, lwvflx_snl, lwnemocoustk, wtauhf, lwnemotauoc, llnormagam, lwflux, lwcou, llgcbz0,  &
  & llcapchnk, lwnemocoustrn, x0tauhf, lwnemocousend_d, lwvflx_snl_d, lwnemocoustk_d, wtauhf_d, lwnemotauoc_d, llnormagam_d,  &
  & lwflux_d, lwcou_d, llgcbz0_d, lwnemocou_d, llcapchnk_d, lwnemocoustrn_d, x0tauhf_d
  USE YOWWIND, ONLY: wspmin, wspmin_d
  USE YOWWNDG, ONLY: icode, icode_cpl, icode_d, icode_cpl_d
  USE YOWSTAT, ONLY: idamping, lbiwbk, iphys, isnonlin, idamping_d, lbiwbk_d, iphys_d, idelt_d, isnonlin_d
  USE YOWALTAS, ONLY: afcrv, bfcrv, egrcrv, afcrv_d, bfcrv_d, egrcrv_d
  USE YOWICE, ONLY: ciblock, lmaskice, cithrsh_tail, cdicwa, cithrsh, lciwabr, licerun, lwamrsetci, ciblock_d, lmaskice_d,  &
  & cithrsh_tail_d, cdicwa_d, cithrsh_d, lciwabr_d, licerun_d, lwamrsetci_d
  USE YOWCOUT, ONLY: lwfluxout, lwfluxout_d
  USE YOWPARAM, ONLY: nfre_odd, nfre_red, llunstr, nfre_d, nfre_odd_d, nang_d, nfre_red_d, llunstr_d
  USE YOWTABL, ONLY: swellft, swellft_d
  
  
  ! ----------------------------------------------------------------------
  
  IMPLICIT NONE
  INTERFACE
    SUBROUTINE INCDATE (CDATE, ISHIFT)
      USE parkind_wave, ONLY: jwim
      INTEGER(KIND=JWIM), INTENT(IN) :: ISHIFT
      CHARACTER(LEN=*), INTENT(INOUT) :: CDATE
    END SUBROUTINE INCDATE
  END INTERFACE
  INTERFACE
    SUBROUTINE NEWWIND (CDATE, CDATEWH, LLNEWFILE, WVPRPT, FF_NOW, FF_NEXT)
      USE YOWDRVTYPE, ONLY: FREQUENCY, FORCING_FIELDS
      CHARACTER(LEN=14), INTENT(IN) :: CDATE
      CHARACTER(LEN=14), INTENT(INOUT) :: CDATEWH
      LOGICAL, INTENT(INOUT) :: LLNEWFILE
      TYPE(FREQUENCY), INTENT(INOUT) :: WVPRPT
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FF_NOW
      TYPE(FORCING_FIELDS), INTENT(IN) :: FF_NEXT
    END SUBROUTINE NEWWIND
  END INTERFACE
  INTERFACE
    SUBROUTINE PROPAG_WAM (BLK2GLO, WVENVI, WVPRPT, FL1)
      USE parkind_wave, ONLY: jwrb
      USE YOWDRVTYPE, ONLY: WVGRIDGLO, ENVIRONMENT, FREQUENCY
      USE yowgrid, ONLY: nproma_wam, nchnk
      USE yowparam, ONLY: nang, nfre
      TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO
      TYPE(ENVIRONMENT), INTENT(IN) :: WVENVI
      TYPE(FREQUENCY), INTENT(IN) :: WVPRPT
      REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK) :: FL1
    END SUBROUTINE PROPAG_WAM
  END INTERFACE
  INTERFACE
    FUNCTION WAM_USER_CLOCK ()
      USE parkind_wave, ONLY: jwru
      REAL(KIND=JWRU) :: WAM_USER_CLOCK
    END FUNCTION WAM_USER_CLOCK
  END INTERFACE
  CHARACTER(LEN=14), INTENT(IN) :: CDTPRA  ! DATE FOR CALL PROPAGATION
  CHARACTER(LEN=14), INTENT(INOUT) :: CDATE  ! CURRENT DATE
  CHARACTER(LEN=14), INTENT(INOUT) :: CDATEWH  ! DATE OF THE NEXT FORCING FIELDS
  CHARACTER(LEN=14), INTENT(INOUT) :: CDTIMP  ! START DATE OF SOURCE FUNCTION INTEGRATION
  CHARACTER(LEN=14), INTENT(INOUT) :: CDTIMPNEXT  ! NEXT START DATE OF SOURCE FUNCTION INTEGRATION
  TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO  ! BLOCK TO GRID TRANSFORMATION
  TYPE(ENVIRONMENT), INTENT(INOUT) :: WVENVI  !  WAVE ENVIRONMENT FIELDS
  TYPE(FREQUENCY), INTENT(INOUT) :: WVPRPT  ! WAVE PROPERTIES FIELDS
  TYPE(FORCING_FIELDS), INTENT(INOUT) :: FF_NOW  ! FORCING FIELDS AT CURRENT TIME
  TYPE(FORCING_FIELDS), INTENT(IN) :: FF_NEXT  !  DATA STRUCTURE WITH THE NEXT FORCING FIELDS
  TYPE(INTGT_PARAM_FIELDS), INTENT(INOUT) :: INTFLDS  ! INTEGRATED/DERIVED PARAMETERS
  TYPE(WAVE2OCEAN), INTENT(INOUT) :: WAM2NEMO  ! WAVE FIELDS PASSED TO NEMO
  INTEGER(KIND=JWIM), INTENT(INOUT), DIMENSION(NPROMA_WAM, NCHNK) :: MIJ  ! LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE
  REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK) :: FL1
  REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK) :: XLLWS  ! TOTAL WINDSEA MASK FROM INPUT SOURCE TERM
  
  REAL(KIND=JWRB), INTENT(INOUT) :: TIME1(3)
  REAL(KIND=JWRB) :: TIME0, TIME2
  
  
  INTEGER(KIND=JWIM) :: IJ, K, M
  INTEGER(KIND=JWIM) :: ICHNK
  INTEGER(KIND=JWIM) :: IDELWH
  
  ! Objects to store fields
  TYPE(FREQUENCY_FIELD) :: WVPRPT_FIELD
  TYPE(ENVIRONMENT_FIELD) :: WVENVI_FIELD
  TYPE(FORCING_FIELDS_FIELD) :: FF_NOW_FIELD
  TYPE(WAVE2OCEAN_FIELD) :: WAM2NEMO_FIELD
  TYPE(INTGT_PARAM_FIELDS_FIELD) :: INTFLDS_FIELD
  TYPE(SOURCE_CONTRIBS_FIELD) :: SRC_CONTRIBS
  
  ! DEVICE POINTERS
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: FL1_DPTR(:, :, :, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: XLLWS_DPTR(:, :, :, :) => NULL()
  INTEGER(KIND=JWIM), POINTER, CONTIGUOUS :: MIJ_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WAVNUM_DPTR(:, :, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CGROUP_DPTR(:, :, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CIWA_DPTR(:, :, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CINV_DPTR(:, :, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: XK2CG_DPTR(:, :, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: STOKFAC_DPTR(:, :, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: EMAXDPT_DPTR(:, :) => NULL()
  INTEGER(KIND=JWIM), POINTER, CONTIGUOUS :: INDEP_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: DEPTH_DPTR(:, :) => NULL()
  INTEGER(KIND=JWIM), POINTER, CONTIGUOUS :: IOBND_DPTR(:, :) => NULL()
  INTEGER(KIND=JWIM), POINTER, CONTIGUOUS :: IODP_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CICOVER_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WSWAVE_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WDWAVE_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: AIRD_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WSTAR_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: UFRIC_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUW_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUWDIR_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: Z0M_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: Z0B_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CHRNCK_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CITHICK_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NEMOUSTOKES_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NEMOVSTOKES_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NEMOSTRN_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NPHIEPS_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NTAUOC_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NSWH_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NMWP_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NEMOTAUX_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NEMOTAUY_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NEMOWSWAVE_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NEMOPHIF_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WSEMEAN_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WSFMEAN_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: USTOKES_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: VSTOKES_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: STRNMS_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUXD_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUYD_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUOCXD_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUOCYD_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUOC_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: PHIOCD_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: PHIEPS_DPTR(:, :) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: PHIAW_DPTR(:, :) => NULL()
  
  REAL(KIND=JWRB), DEVICE, ALLOCATABLE :: ENH_D(:,:,:)
!  real(kind=jwrb), device, allocatable :: ciwa_d(:,:,:)
  
  
  LOGICAL :: LLNEWFILE
  
  DATA LLNEWFILE / .false. /
  INTEGER :: istat
  TYPE(DIM3) :: GRIDDIM, BLOCKDIM

!  allocate(ciwa_d(nproma_wam, nfre, nchnk))
  ALLOCATE(ENH_D(NPROMA_WAM, MLSTHG, NCHNK))

!  ciwa_d = wvprpt%ciwa

  ALLOCATE (indicessat_d, SOURCE=indicessat)
  CALL ACC_MAP_DATA(indicessat, indicessat_d, SIZE(indicessat)*C_SIZEOF(INT(1, kind=JWIM)))
  CALL ACC_MAP_DATA(isnonlin, isnonlin_d, C_SIZEOF(INT(1, kind=JWIM)))
  CALL ACC_MAP_DATA(ciblock, ciblock_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(kfrh, kfrh_d, C_SIZEOF(INT(1, kind=JWIM)))
  CALL ACC_MAP_DATA(tailfactor, tailfactor_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(cdis, cdis_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(lwnemocousend, lwnemocousend_d, C_SIZEOF(LOGICAL(.true.)))
  CALL ACC_MAP_DATA(rnum, rnum_d, C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (ikm1_d, SOURCE=ikm1)
  CALL ACC_MAP_DATA(ikm1, ikm1_d, SIZE(ikm1)*C_SIZEOF(INT(1, kind=JWIM)))
  CALL ACC_MAP_DATA(zpi4gm1, zpi4gm1_d, C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (xk_gc_d, SOURCE=xk_gc)
  CALL ACC_MAP_DATA(xk_gc, xk_gc_d, SIZE(xk_gc)*C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (zpifr_d, SOURCE=zpifr)
  CALL ACC_MAP_DATA(zpifr, zpifr_d, SIZE(zpifr)*C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(lwflux, lwflux_d, C_SIZEOF(LOGICAL(.true.)))
  CALL ACC_MAP_DATA(alphamin, alphamin_d, C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (dfim_sim_d, SOURCE=dfim_sim)
  CALL ACC_MAP_DATA(dfim_sim, dfim_sim_d, SIZE(dfim_sim)*C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(dthrn_u, dthrn_u_d, C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (af11_d, SOURCE=af11)
  CALL ACC_MAP_DATA(af11, af11_d, SIZE(af11)*C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(MLSTHG, mlsthg_d, C_SIZEOF(INT(1, kind=JWIM)))
  CALL ACC_MAP_DATA(cithrsh_tail, cithrsh_tail_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(lwnemocoustrn, lwnemocoustrn_d, C_SIZEOF(LOGICAL(.true.)))
  ALLOCATE (omega_gc_d, SOURCE=omega_gc)
  CALL ACC_MAP_DATA(omega_gc, omega_gc_d, SIZE(omega_gc)*C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(swellf5, swellf5_d, C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (k21w_d, SOURCE=k21w)
  CALL ACC_MAP_DATA(k21w, k21w_d, SIZE(k21w)*C_SIZEOF(INT(1, kind=JWIM)))
  CALL ACC_MAP_DATA(gm1, gm1_d, C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (om3gmkm_gc_d, SOURCE=om3gmkm_gc)
  CALL ACC_MAP_DATA(om3gmkm_gc, om3gmkm_gc_d, SIZE(om3gmkm_gc)*C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(lwamrsetci, lwamrsetci_d, C_SIZEOF(LOGICAL(.true.)))
  CALL ACC_MAP_DATA(dal1, dal1_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(tailfactor_pm, tailfactor_pm_d, C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (k2w_d, SOURCE=k2w)
  CALL ACC_MAP_DATA(k2w, k2w_d, SIZE(k2w)*C_SIZEOF(INT(1, kind=JWIM)))
  ALLOCATE (dfim_d, SOURCE=dfim)
  CALL ACC_MAP_DATA(dfim, dfim_d, SIZE(dfim)*C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (fr_d, SOURCE=fr)
  CALL ACC_MAP_DATA(fr, fr_d, SIZE(fr)*C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(lwcou, lwcou_d, C_SIZEOF(LOGICAL(.true.)))
  ALLOCATE (delkcc_omxkm3_gc_d, SOURCE=delkcc_omxkm3_gc)
  CALL ACC_MAP_DATA(delkcc_omxkm3_gc, delkcc_omxkm3_gc_d, SIZE(delkcc_omxkm3_gc)*C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (k11w_d, SOURCE=k11w)
  CALL ACC_MAP_DATA(k11w, k11w_d, SIZE(k11w)*C_SIZEOF(INT(1, kind=JWIM)))
  ALLOCATE (rhowg_dfim_d, SOURCE=rhowg_dfim)
  CALL ACC_MAP_DATA(rhowg_dfim, rhowg_dfim_d, SIZE(rhowg_dfim)*C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (ikp1_d, SOURCE=ikp1)
  CALL ACC_MAP_DATA(ikp1, ikp1_d, SIZE(ikp1)*C_SIZEOF(INT(1, kind=JWIM)))
  ALLOCATE (fklap_d, SOURCE=fklap)
  CALL ACC_MAP_DATA(fklap, fklap_d, SIZE(fklap)*C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (fklap1_d, SOURCE=fklap1)
  CALL ACC_MAP_DATA(fklap1, fklap1_d, SIZE(fklap1)*C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(z0tubmax, z0tubmax_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(lmaskice, lmaskice_d, C_SIZEOF(LOGICAL(.true.)))
  ALLOCATE (inlcoef_d, SOURCE=inlcoef)
  CALL ACC_MAP_DATA(inlcoef, inlcoef_d, SIZE(inlcoef)*C_SIZEOF(INT(1, kind=JWIM)))
  ALLOCATE (c2osqrtvg_gc_d, SOURCE=c2osqrtvg_gc)
  CALL ACC_MAP_DATA(c2osqrtvg_gc, c2osqrtvg_gc_d, SIZE(c2osqrtvg_gc)*C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(idamping, idamping_d, C_SIZEOF(INT(1, kind=JWIM)))
  CALL ACC_MAP_DATA(sqrtgosurft, sqrtgosurft_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(g, g_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(zpi4gm2, zpi4gm2_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(chnkmin_u, chnkmin_u_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(wspmin, wspmin_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(cdicwa, cdicwa_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(nfre_red, nfre_red_d, C_SIZEOF(INT(1, kind=JWIM)))
  ALLOCATE (ikp_d, SOURCE=ikp)
  CALL ACC_MAP_DATA(ikp, ikp_d, SIZE(ikp)*C_SIZEOF(INT(1, kind=JWIM)))
  CALL ACC_MAP_DATA(egrcrv, egrcrv_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(icode, icode_d, C_SIZEOF(INT(1, kind=JWIM)))
  ALLOCATE (th_d, SOURCE=th)
  CALL ACC_MAP_DATA(th, th_d, SIZE(th)*C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(betamaxoxkappa2, betamaxoxkappa2_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(x0tauhf, x0tauhf_d, C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (fr5_d, SOURCE=fr5)
  CALL ACC_MAP_DATA(fr5, fr5_d, SIZE(fr5)*C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (xkmsqrtvgoc2_gc_d, SOURCE=xkmsqrtvgoc2_gc)
  CALL ACC_MAP_DATA(xkmsqrtvgoc2_gc, xkmsqrtvgoc2_gc_d, SIZE(xkmsqrtvgoc2_gc)*C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(zpi, zpi_d, C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (k1w_d, SOURCE=k1w)
  CALL ACC_MAP_DATA(k1w, k1w_d, SIZE(k1w)*C_SIZEOF(INT(1, kind=JWIM)))
  CALL ACC_MAP_DATA(delth, delth_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(llnormagam, llnormagam_d, C_SIZEOF(LOGICAL(.true.)))
  CALL ACC_MAP_DATA(dal2, dal2_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(wtauhf, wtauhf_d, SIZE(wtauhf)*C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (rnlcoef_d, SOURCE=rnlcoef)
  CALL ACC_MAP_DATA(rnlcoef, rnlcoef_d, SIZE(rnlcoef)*C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(alpha, alpha_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(dthrn_a, dthrn_a_d, C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (cumulw_d, SOURCE=cumulw)
  CALL ACC_MAP_DATA(cumulw, cumulw_d, SIZE(cumulw)*C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (omxkm3_gc_d, SOURCE=omxkm3_gc)
  CALL ACC_MAP_DATA(omxkm3_gc, omxkm3_gc_d, SIZE(omxkm3_gc)*C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (dfimfr_d, SOURCE=dfimfr)
  CALL ACC_MAP_DATA(dfimfr, dfimfr_d, SIZE(dfimfr)*C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(ang_gc_c, ang_gc_c_d, C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (sinth_d, SOURCE=sinth)
  CALL ACC_MAP_DATA(sinth, sinth_d, SIZE(sinth)*C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(afcrv, afcrv_d, C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (cofrm4_d, SOURCE=cofrm4)
  CALL ACC_MAP_DATA(cofrm4, cofrm4_d, SIZE(cofrm4)*C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (ikm_d, SOURCE=ikm)
  CALL ACC_MAP_DATA(ikm, ikm_d, SIZE(ikm)*C_SIZEOF(INT(1, kind=JWIM)))
  CALL ACC_MAP_DATA(mfrstlw, mfrstlw_d, C_SIZEOF(INT(1, kind=JWIM)))
  CALL ACC_MAP_DATA(zalp, zalp_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(nfre_odd, nfre_odd_d, C_SIZEOF(INT(1, kind=JWIM)))
  CALL ACC_MAP_DATA(lwvflx_snl, lwvflx_snl_d, C_SIZEOF(LOGICAL(.true.)))
  CALL ACC_MAP_DATA(nsdsnth, nsdsnth_d, C_SIZEOF(INT(1, kind=JWIM)))
  CALL ACC_MAP_DATA(lciwabr, lciwabr_d, C_SIZEOF(LOGICAL(.true.)))
  CALL ACC_MAP_DATA(ang_gc_b, ang_gc_b_d, C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (fklam1_d, SOURCE=fklam1)
  CALL ACC_MAP_DATA(fklam1, fklam1_d, SIZE(fklam1)*C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(NFRE, nfre_d, C_SIZEOF(INT(1, kind=JWIM)))
  ALLOCATE (xkm_gc_d, SOURCE=xkm_gc)
  CALL ACC_MAP_DATA(xkm_gc, xkm_gc_d, SIZE(xkm_gc)*C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (fklam_d, SOURCE=fklam)
  CALL ACC_MAP_DATA(fklam, fklam_d, SIZE(fklam)*C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(gamnconst, gamnconst_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(ang_gc_a, ang_gc_a_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(iphys, iphys_d, C_SIZEOF(INT(1, kind=JWIM)))
  CALL ACC_MAP_DATA(icode_cpl, icode_cpl_d, C_SIZEOF(INT(1, kind=JWIM)))
  CALL ACC_MAP_DATA(LWNEMOCOU, lwnemocou_d, C_SIZEOF(LOGICAL(.true.)))
  ALLOCATE (delkcc_gc_ns_d, SOURCE=delkcc_gc_ns)
  CALL ACC_MAP_DATA(delkcc_gc_ns, delkcc_gc_ns_d, SIZE(delkcc_gc_ns)*C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (flmax_d, SOURCE=flmax)
  CALL ACC_MAP_DATA(flmax, flmax_d, SIZE(flmax)*C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (dfimofr_d, SOURCE=dfimofr)
  CALL ACC_MAP_DATA(dfimofr, dfimofr_d, SIZE(dfimofr)*C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(alphapmax, alphapmax_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(z0rat, z0rat_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(licerun, licerun_d, C_SIZEOF(LOGICAL(.true.)))
  CALL ACC_MAP_DATA(ndikcumul, ndikcumul_d, C_SIZEOF(INT(1, kind=JWIM)))
  CALL ACC_MAP_DATA(lwfluxout, lwfluxout_d, C_SIZEOF(LOGICAL(.true.)))
  CALL ACC_MAP_DATA(delta_sdis, delta_sdis_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(llunstr, llunstr_d, C_SIZEOF(LOGICAL(.true.)))
  CALL ACC_MAP_DATA(bmaxokap, bmaxokap_d, C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (costh_d, SOURCE=costh)
  CALL ACC_MAP_DATA(costh, costh_d, SIZE(costh)*C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(lwnemocoustk, lwnemocoustk_d, C_SIZEOF(LOGICAL(.true.)))
  CALL ACC_MAP_DATA(tauwshelter, tauwshelter_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(swellft, swellft_d, SIZE(swellft)*C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(llcapchnk, llcapchnk_d, C_SIZEOF(LOGICAL(.true.)))
  CALL ACC_MAP_DATA(flogsprdm1, flogsprdm1_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(nwav_gc, nwav_gc_d, C_SIZEOF(INT(1, kind=JWIM)))
  CALL ACC_MAP_DATA(rnu, rnu_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(lbiwbk, lbiwbk_d, C_SIZEOF(LOGICAL(.true.)))
  CALL ACC_MAP_DATA(lwnemotauoc, lwnemotauoc_d, C_SIZEOF(LOGICAL(.true.)))
  ALLOCATE (cm_gc_d, SOURCE=cm_gc)
  CALL ACC_MAP_DATA(cm_gc, cm_gc_d, SIZE(cm_gc)*C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(IDELT, idelt_d, C_SIZEOF(INT(1, kind=JWIM)))
  CALL ACC_MAP_DATA(cithrsh, cithrsh_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(cdisvis, cdisvis_d, C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(rn1_rn, rn1_rn_d, C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (satweights_d, SOURCE=satweights)
  CALL ACC_MAP_DATA(satweights, satweights_d, SIZE(satweights)*C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(bfcrv, bfcrv_d, C_SIZEOF(REAL(1, kind=JWRB)))
  ALLOCATE (dfimfr2_d, SOURCE=dfimfr2)
  CALL ACC_MAP_DATA(dfimfr2, dfimfr2_d, SIZE(dfimfr2)*C_SIZEOF(REAL(1, kind=JWRB)))
  CALL ACC_MAP_DATA(NANG, nang_d, C_SIZEOF(INT(1, kind=JWIM)))
  CALL ACC_MAP_DATA(llgcbz0, llgcbz0_d, C_SIZEOF(LOGICAL(.true.)))
  
  ! ----------------------------------------------------------------------
  
  
  !*     PROPAGATION TIME
  !      ----------------
  
  IF (CDATE == CDTPRA) THEN
    TIME0 = -WAM_USER_CLOCK()
    CALL PROPAG_WAM(BLK2GLO, WVENVI, WVPRPT, FL1)
    TIME1(1) = TIME1(1) + (TIME0 + WAM_USER_CLOCK())*1.E-06
    CDATE = CDTPRO
  END IF
  
  
  !* RETRIEVING NEW FORCING FIELDS IF NEEDED.
  !  ----------------------------------------
  CALL NEWWIND(CDTIMP, CDATEWH, LLNEWFILE, WVPRPT, FF_NOW, FF_NEXT)
  
  ! IT IS TIME TO INTEGRATE THE SOURCE TERMS
  ! ----------------------------------------
  IF (CDATE >= CDTIMPNEXT) THEN
    ! COMPUTE UPDATE DUE TO SOURCE TERMS
    CALL GSTATS(1431, 0)
    IF (LLSOURCE) THEN
      
      TIME2 = -WAM_USER_CLOCK()
      CALL WVPRPT_FIELD%INIT(WAVNUM=WVPRPT%WAVNUM, CGROUP=WVPRPT%CGROUP, CIWA=WVPRPT%CIWA, CINV=WVPRPT%CINV, XK2CG=WVPRPT%XK2CG,  &
      & STOKFAC=WVPRPT%STOKFAC)
      CALL WVENVI_FIELD%INIT(EMAXDPT=WVENVI%EMAXDPT, INDEP=WVENVI%INDEP, DEPTH=WVENVI%DEPTH, IOBND=WVENVI%IOBND, IODP=WVENVI%IODP &
      & )
      CALL FF_NOW_FIELD%INIT(AIRD=FF_NOW%AIRD, WDWAVE=FF_NOW%WDWAVE, CICOVER=FF_NOW%CICOVER, WSWAVE=FF_NOW%WSWAVE,  &
      & WSTAR=FF_NOW%WSTAR, UFRIC=FF_NOW%UFRIC, TAUW=FF_NOW%TAUW, TAUWDIR=FF_NOW%TAUWDIR, Z0M=FF_NOW%Z0M, Z0B=FF_NOW%Z0B,  &
      & CHRNCK=FF_NOW%CHRNCK, CITHICK=FF_NOW%CITHICK)
      CALL WAM2NEMO_FIELD%INIT(NEMOUSTOKES=WAM2NEMO%NEMOUSTOKES, NEMOVSTOKES=WAM2NEMO%NEMOVSTOKES, NEMOSTRN=WAM2NEMO%NEMOSTRN,  &
      & NPHIEPS=WAM2NEMO%NPHIEPS, NTAUOC=WAM2NEMO%NTAUOC, NSWH=WAM2NEMO%NSWH, NMWP=WAM2NEMO%NMWP, NEMOTAUX=WAM2NEMO%NEMOTAUX,  &
      & NEMOTAUY=WAM2NEMO%NEMOTAUY, NEMOWSWAVE=WAM2NEMO%NEMOWSWAVE, NEMOPHIF=WAM2NEMO%NEMOPHIF)
      CALL INTFLDS_FIELD%INIT(WSEMEAN=INTFLDS%WSEMEAN, WSFMEAN=INTFLDS%WSFMEAN, USTOKES=INTFLDS%USTOKES,  &
      & VSTOKES=INTFLDS%VSTOKES, STRNMS=INTFLDS%STRNMS, TAUXD=INTFLDS%TAUXD, TAUYD=INTFLDS%TAUYD, TAUOCXD=INTFLDS%TAUOCXD,  &
      & TAUOCYD=INTFLDS%TAUOCYD, TAUOC=INTFLDS%TAUOC, PHIOCD=INTFLDS%PHIOCD, PHIEPS=INTFLDS%PHIEPS, PHIAW=INTFLDS%PHIAW)
      CALL SRC_CONTRIBS%INIT(FL1=FL1, XLLWS=XLLWS, MIJ=MIJ)
      
!$acc update device(  &
!$acc & fr,dfimfr2,k11w,dthrn_a,dal1,zpi,ndikcumul,nfre,k1w,k2w,mfrstlw,lwnemocousend,rnu,ikm1,flogsprdm1,xkm_gc,lwvflx_snl,sqrtgosurft,omega_gc,bmaxokap,delkcc_gc_ns,chnkmin_u,lwnemocoustk,zpi4gm1,lwnemotauoc,rnlcoef,fklap,tailfactor,wspmin,icode,llnormagam,idamping,indicessat,rn1_rn,satweights,flmax,alpha,c2osqrtvg_gc,sinth,z0rat,cofrm4,dfimfr,afcrv,tailfactor_pm,lwflux,ciblock,dal2,dfim_sim,lwfluxout,rhowg_dfim,lwcou,zpi4gm2,cdisvis,fklam1,costh,cumulw,nfre_odd,gm1,lbiwbk,delth,swellft,lmaskice,xkmsqrtvgoc2_gc,xk_gc,bfcrv,alphapmax,dthrn_u,nsdsnth,dfimofr,delta_sdis,cdis,iphys,cithrsh_tail,tauwshelter,ikp,nang,alphamin,omxkm3_gc,egrcrv,cdicwa,z0tubmax,idelt,swellf5,nfre_red,ang_gc_c,zalp,cithrsh,ang_gc_a,inlcoef,icode_cpl,k21w,llgcbz0,lwnemocou,g,th,gamnconst,fklam,zpifr,dfim,af11,ikp1,fklap1,llcapchnk,delkcc_omxkm3_gc,ikm,kfrh,fr5,cm_gc,x0tauhf,lciwabr,licerun,betamaxoxkappa2,mlsthg,llunstr,om3gmkm_gc,ang_gc_b,lwamrsetci,lwnemocoustrn,nwav_gc,wtauhf,isnonlin,rnum &
!$acc &  )
      
      CALL WVPRPT_FIELD%UPDATE_DEVICE(WAVNUM=WAVNUM_DPTR, CGROUP=CGROUP_DPTR, CIWA=CIWA_DPTR, CINV=CINV_DPTR, XK2CG=XK2CG_DPTR,  &
      & STOKFAC=STOKFAC_DPTR)
      CALL WVENVI_FIELD%UPDATE_DEVICE(EMAXDPT=EMAXDPT_DPTR, INDEP=INDEP_DPTR, DEPTH=DEPTH_DPTR, IOBND=IOBND_DPTR, IODP=IODP_DPTR)
      CALL FF_NOW_FIELD%UPDATE_DEVICE(AIRD=AIRD_DPTR, WDWAVE=WDWAVE_DPTR, CICOVER=CICOVER_DPTR, WSWAVE=WSWAVE_DPTR,  &
      & WSTAR=WSTAR_DPTR, UFRIC=UFRIC_DPTR, TAUW=TAUW_DPTR, TAUWDIR=TAUWDIR_DPTR, Z0M=Z0M_DPTR, Z0B=Z0B_DPTR,  &
      & CHRNCK=CHRNCK_DPTR, CITHICK=CITHICK_DPTR)
      CALL WAM2NEMO_FIELD%UPDATE_DEVICE(NEMOUSTOKES=NEMOUSTOKES_DPTR, NEMOVSTOKES=NEMOVSTOKES_DPTR, NEMOSTRN=NEMOSTRN_DPTR,  &
      & NPHIEPS=NPHIEPS_DPTR, NTAUOC=NTAUOC_DPTR, NSWH=NSWH_DPTR, NMWP=NMWP_DPTR, NEMOTAUX=NEMOTAUX_DPTR,  &
      & NEMOTAUY=NEMOTAUY_DPTR, NEMOWSWAVE=NEMOWSWAVE_DPTR, NEMOPHIF=NEMOPHIF_DPTR)
      CALL INTFLDS_FIELD%UPDATE_DEVICE(WSEMEAN=WSEMEAN_DPTR, WSFMEAN=WSFMEAN_DPTR, USTOKES=USTOKES_DPTR, VSTOKES=VSTOKES_DPTR,  &
      & STRNMS=STRNMS_DPTR, TAUXD=TAUXD_DPTR, TAUYD=TAUYD_DPTR, TAUOCXD=TAUOCXD_DPTR, TAUOCYD=TAUOCYD_DPTR, TAUOC=TAUOC_DPTR,  &
      & PHIOCD=PHIOCD_DPTR, PHIEPS=PHIEPS_DPTR, PHIAW=PHIAW_DPTR)
      CALL SRC_CONTRIBS%UPDATE_DEVICE(FL1=FL1_DPTR, XLLWS=XLLWS_DPTR, MIJ=MIJ_DPTR)
      
!$acc data present(FL1_DPTR,XLLWS_DPTR,MIJ_DPTR,WAVNUM_DPTR,CGROUP_DPTR,CIWA_DPTR,CINV_DPTR,XK2CG_DPTR,STOKFAC_DPTR,&
!$acc &            EMAXDPT_DPTR,INDEP_DPTR,DEPTH_DPTR,IOBND_DPTR,IODP_DPTR,CICOVER_DPTR,WSWAVE_DPTR,WDWAVE_DPTR,AIRD_DPTR,&
!$acc &            WSTAR_DPTR,UFRIC_DPTR,TAUW_DPTR,TAUWDIR_DPTR,Z0M_DPTR,Z0B_DPTR,CHRNCK_DPTR,CITHICK_DPTR,NEMOUSTOKES_DPTR,&
!$acc &            NEMOVSTOKES_DPTR,NEMOSTRN_DPTR,NPHIEPS_DPTR,NTAUOC_DPTR,NSWH_DPTR,NMWP_DPTR,NEMOTAUX_DPTR,NEMOTAUY_DPTR,&
!$acc &            NEMOWSWAVE_DPTR,NEMOPHIF_DPTR,WSEMEAN_DPTR,WSFMEAN_DPTR,USTOKES_DPTR,VSTOKES_DPTR,STRNMS_DPTR,TAUXD_DPTR,&
!$acc &            TAUYD_DPTR,TAUOCXD_DPTR,TAUOCYD_DPTR,TAUOC_DPTR,PHIOCD_DPTR,PHIEPS_DPTR,PHIAW_DPTR)
      TIME0 = -WAM_USER_CLOCK()
      GRIDDIM = DIM3(1, 1, CEILING(REAL(NCHNK) / REAL(NPROMA_WAM)))
      BLOCKDIM = DIM3(NPROMA_WAM, 1, 1)
!$acc host_data use_device(FL1_DPTR,XLLWS_DPTR,MIJ_DPTR,WAVNUM_DPTR,CGROUP_DPTR,CIWA_DPTR,CINV_DPTR,XK2CG_DPTR,STOKFAC_DPTR,&
!$acc &            EMAXDPT_DPTR,INDEP_DPTR,DEPTH_DPTR,IOBND_DPTR,IODP_DPTR,CICOVER_DPTR,WSWAVE_DPTR,WDWAVE_DPTR,AIRD_DPTR,&
!$acc &            WSTAR_DPTR,UFRIC_DPTR,TAUW_DPTR,TAUWDIR_DPTR,Z0M_DPTR,Z0B_DPTR,CHRNCK_DPTR,CITHICK_DPTR,NEMOUSTOKES_DPTR,&
!$acc &            NEMOVSTOKES_DPTR,NEMOSTRN_DPTR,NPHIEPS_DPTR,NTAUOC_DPTR,NSWH_DPTR,NMWP_DPTR,NEMOTAUX_DPTR,NEMOTAUY_DPTR,&
!$acc &            NEMOWSWAVE_DPTR,NEMOPHIF_DPTR,WSEMEAN_DPTR,WSFMEAN_DPTR,USTOKES_DPTR,VSTOKES_DPTR,STRNMS_DPTR,TAUXD_DPTR,&
!$acc &            TAUYD_DPTR,TAUOCXD_DPTR,TAUOCYD_DPTR,TAUOC_DPTR,PHIOCD_DPTR,PHIEPS_DPTR,PHIAW_DPTR)
      CALL IMPLSCH_CUF_PARAMETRISE<<<GRIDDIM,BLOCKDIM>>>(1, NPROMA_WAM, FL1_DPTR, WAVNUM_DPTR, CGROUP_DPTR, CIWA_DPTR,  &
      & CINV_DPTR, XK2CG_DPTR, STOKFAC_DPTR, EMAXDPT_DPTR, INDEP_DPTR, DEPTH_DPTR, IOBND_DPTR, IODP_DPTR, AIRD_DPTR,  &
      & WDWAVE_DPTR, CICOVER_DPTR, WSWAVE_DPTR, WSTAR_DPTR, UFRIC_DPTR, TAUW_DPTR, TAUWDIR_DPTR, Z0M_DPTR, Z0B_DPTR,  &
      & CHRNCK_DPTR, CITHICK_DPTR, NEMOUSTOKES_DPTR, NEMOVSTOKES_DPTR, NEMOSTRN_DPTR, NPHIEPS_DPTR, NTAUOC_DPTR, NSWH_DPTR,  &
      & NMWP_DPTR, NEMOTAUX_DPTR, NEMOTAUY_DPTR, NEMOWSWAVE_DPTR, NEMOPHIF_DPTR, WSEMEAN_DPTR, WSFMEAN_DPTR, USTOKES_DPTR,  &
      & VSTOKES_DPTR, STRNMS_DPTR, TAUXD_DPTR, TAUYD_DPTR, TAUOCXD_DPTR, TAUOCYD_DPTR, TAUOC_DPTR, PHIOCD_DPTR, PHIEPS_DPTR,  &
      & PHIAW_DPTR, MIJ_DPTR, XLLWS_DPTR, ENH_D, NCHNK)
!$acc end host_data
      istat = cudaDeviceSynchronize()
      TIME1(2) = TIME1(2) + (TIME0 + WAM_USER_CLOCK())*1.E-06
      CALL WVPRPT_FIELD%ENSURE_HOST()
      CALL WVENVI_FIELD%ENSURE_HOST()
      CALL FF_NOW_FIELD%ENSURE_HOST()
      CALL WAM2NEMO_FIELD%ENSURE_HOST()
      CALL INTFLDS_FIELD%ENSURE_HOST()
      CALL SRC_CONTRIBS%ENSURE_HOST()
!      wvprpt%ciwa = ciwa_d
      
      TIME1(3) = TIME1(3) + (TIME2 + WAM_USER_CLOCK())*1.E-06
!$acc end data
      
      IF (LWNEMOCOU) NEMONTAU = NEMONTAU + 1
      
    ELSE
      !   NO SOURCE TERM CONTRIBUTION
      !$OMP      PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK)
          DO ICHNK = 1, NCHNK
            MIJ(:,ICHNK) = NFRE
            FL1(:,:,:,ICHNK) = MAX(FL1(:,:,:,ICHNK), EPSMIN)
            XLLWS(:,:,:,ICHNK) = 0.0_JWRB
          ENDDO
      !$OMP      END PARALLEL DO
    END IF
    CALL GSTATS(1431, 1)
    
    
    !*       UPDATE FORCING FIELDS TIME COUNTER
    !        ----------------------------------
    IF (LLNEWFILE) THEN
      LLNEWFILE = .false.
      IDELWH = MAX(IDELWI, IDELPRO)
      CALL INCDATE(CDAWIFL, IDELWH)
      CALL INCDATE(CDATEFL, IDELWH)
    END IF
    
    CDATEWO = CDATEWH
    CDTIMP = CDTIMPNEXT
    CALL INCDATE(CDTIMPNEXT, IDELT)
    
  END IF
  
!  deallocate(ciwa_d)
  DEALLOCATE(ENH_D)
  
  CALL ACC_UNMAP_DATA(fr)
  DEALLOCATE (fr_d)
  CALL ACC_UNMAP_DATA(dfimfr2)
  DEALLOCATE (dfimfr2_d)
  CALL ACC_UNMAP_DATA(k11w)
  DEALLOCATE (k11w_d)
  CALL ACC_UNMAP_DATA(dthrn_a)
  CALL ACC_UNMAP_DATA(dal1)
  CALL ACC_UNMAP_DATA(zpi)
  CALL ACC_UNMAP_DATA(ndikcumul)
  CALL ACC_UNMAP_DATA(NFRE)
  CALL ACC_UNMAP_DATA(k1w)
  DEALLOCATE (k1w_d)
  CALL ACC_UNMAP_DATA(k2w)
  DEALLOCATE (k2w_d)
  CALL ACC_UNMAP_DATA(mfrstlw)
  CALL ACC_UNMAP_DATA(lwnemocousend)
  CALL ACC_UNMAP_DATA(rnu)
  CALL ACC_UNMAP_DATA(ikm1)
  DEALLOCATE (ikm1_d)
  CALL ACC_UNMAP_DATA(flogsprdm1)
  CALL ACC_UNMAP_DATA(xkm_gc)
  DEALLOCATE (xkm_gc_d)
  CALL ACC_UNMAP_DATA(lwvflx_snl)
  CALL ACC_UNMAP_DATA(sqrtgosurft)
  CALL ACC_UNMAP_DATA(omega_gc)
  DEALLOCATE (omega_gc_d)
  CALL ACC_UNMAP_DATA(bmaxokap)
  CALL ACC_UNMAP_DATA(delkcc_gc_ns)
  DEALLOCATE (delkcc_gc_ns_d)
  CALL ACC_UNMAP_DATA(chnkmin_u)
  CALL ACC_UNMAP_DATA(lwnemocoustk)
  CALL ACC_UNMAP_DATA(zpi4gm1)
  CALL ACC_UNMAP_DATA(wtauhf)
  CALL ACC_UNMAP_DATA(lwnemotauoc)
  CALL ACC_UNMAP_DATA(rnlcoef)
  DEALLOCATE (rnlcoef_d)
  CALL ACC_UNMAP_DATA(fklap)
  DEALLOCATE (fklap_d)
  CALL ACC_UNMAP_DATA(tailfactor)
  CALL ACC_UNMAP_DATA(wspmin)
  CALL ACC_UNMAP_DATA(icode)
  CALL ACC_UNMAP_DATA(llnormagam)
  CALL ACC_UNMAP_DATA(idamping)
  CALL ACC_UNMAP_DATA(indicessat)
  DEALLOCATE (indicessat_d)
  CALL ACC_UNMAP_DATA(rn1_rn)
  CALL ACC_UNMAP_DATA(satweights)
  DEALLOCATE (satweights_d)
  CALL ACC_UNMAP_DATA(flmax)
  DEALLOCATE (flmax_d)
  CALL ACC_UNMAP_DATA(alpha)
  CALL ACC_UNMAP_DATA(c2osqrtvg_gc)
  DEALLOCATE (c2osqrtvg_gc_d)
  CALL ACC_UNMAP_DATA(sinth)
  DEALLOCATE (sinth_d)
  CALL ACC_UNMAP_DATA(z0rat)
  CALL ACC_UNMAP_DATA(cofrm4)
  DEALLOCATE (cofrm4_d)
  CALL ACC_UNMAP_DATA(dfimfr)
  DEALLOCATE (dfimfr_d)
  CALL ACC_UNMAP_DATA(afcrv)
  CALL ACC_UNMAP_DATA(tailfactor_pm)
  CALL ACC_UNMAP_DATA(lwflux)
  CALL ACC_UNMAP_DATA(ciblock)
  CALL ACC_UNMAP_DATA(dal2)
  CALL ACC_UNMAP_DATA(dfim_sim)
  DEALLOCATE (dfim_sim_d)
  CALL ACC_UNMAP_DATA(lwfluxout)
  CALL ACC_UNMAP_DATA(rhowg_dfim)
  DEALLOCATE (rhowg_dfim_d)
  CALL ACC_UNMAP_DATA(lwcou)
  CALL ACC_UNMAP_DATA(zpi4gm2)
  CALL ACC_UNMAP_DATA(cdisvis)
  CALL ACC_UNMAP_DATA(fklam1)
  DEALLOCATE (fklam1_d)
  CALL ACC_UNMAP_DATA(costh)
  DEALLOCATE (costh_d)
  CALL ACC_UNMAP_DATA(cumulw)
  DEALLOCATE (cumulw_d)
  CALL ACC_UNMAP_DATA(nfre_odd)
  CALL ACC_UNMAP_DATA(gm1)
  CALL ACC_UNMAP_DATA(lbiwbk)
  CALL ACC_UNMAP_DATA(delth)
  CALL ACC_UNMAP_DATA(swellft)
  CALL ACC_UNMAP_DATA(lmaskice)
  CALL ACC_UNMAP_DATA(xkmsqrtvgoc2_gc)
  DEALLOCATE (xkmsqrtvgoc2_gc_d)
  CALL ACC_UNMAP_DATA(xk_gc)
  DEALLOCATE (xk_gc_d)
  CALL ACC_UNMAP_DATA(bfcrv)
  CALL ACC_UNMAP_DATA(alphapmax)
  CALL ACC_UNMAP_DATA(dthrn_u)
  CALL ACC_UNMAP_DATA(nsdsnth)
  CALL ACC_UNMAP_DATA(dfimofr)
  DEALLOCATE (dfimofr_d)
  CALL ACC_UNMAP_DATA(delta_sdis)
  CALL ACC_UNMAP_DATA(cdis)
  CALL ACC_UNMAP_DATA(iphys)
  CALL ACC_UNMAP_DATA(cithrsh_tail)
  CALL ACC_UNMAP_DATA(tauwshelter)
  CALL ACC_UNMAP_DATA(ikp)
  DEALLOCATE (ikp_d)
  CALL ACC_UNMAP_DATA(NANG)
  CALL ACC_UNMAP_DATA(alphamin)
  CALL ACC_UNMAP_DATA(omxkm3_gc)
  DEALLOCATE (omxkm3_gc_d)
  CALL ACC_UNMAP_DATA(egrcrv)
  CALL ACC_UNMAP_DATA(cdicwa)
  CALL ACC_UNMAP_DATA(z0tubmax)
  CALL ACC_UNMAP_DATA(IDELT)
  CALL ACC_UNMAP_DATA(swellf5)
  CALL ACC_UNMAP_DATA(nfre_red)
  CALL ACC_UNMAP_DATA(ang_gc_c)
  CALL ACC_UNMAP_DATA(zalp)
  CALL ACC_UNMAP_DATA(cithrsh)
  CALL ACC_UNMAP_DATA(ang_gc_a)
  CALL ACC_UNMAP_DATA(inlcoef)
  DEALLOCATE (inlcoef_d)
  CALL ACC_UNMAP_DATA(icode_cpl)
  CALL ACC_UNMAP_DATA(k21w)
  DEALLOCATE (k21w_d)
  CALL ACC_UNMAP_DATA(llgcbz0)
  CALL ACC_UNMAP_DATA(LWNEMOCOU)
  CALL ACC_UNMAP_DATA(g)
  CALL ACC_UNMAP_DATA(th)
  DEALLOCATE (th_d)
  CALL ACC_UNMAP_DATA(gamnconst)
  CALL ACC_UNMAP_DATA(fklam)
  DEALLOCATE (fklam_d)
  CALL ACC_UNMAP_DATA(zpifr)
  DEALLOCATE (zpifr_d)
  CALL ACC_UNMAP_DATA(dfim)
  DEALLOCATE (dfim_d)
  CALL ACC_UNMAP_DATA(af11)
  DEALLOCATE (af11_d)
  CALL ACC_UNMAP_DATA(ikp1)
  DEALLOCATE (ikp1_d)
  CALL ACC_UNMAP_DATA(fklap1)
  DEALLOCATE (fklap1_d)
  CALL ACC_UNMAP_DATA(llcapchnk)
  CALL ACC_UNMAP_DATA(delkcc_omxkm3_gc)
  DEALLOCATE (delkcc_omxkm3_gc_d)
  CALL ACC_UNMAP_DATA(ikm)
  DEALLOCATE (ikm_d)
  CALL ACC_UNMAP_DATA(kfrh)
  CALL ACC_UNMAP_DATA(fr5)
  DEALLOCATE (fr5_d)
  CALL ACC_UNMAP_DATA(cm_gc)
  DEALLOCATE (cm_gc_d)
  CALL ACC_UNMAP_DATA(lciwabr)
  CALL ACC_UNMAP_DATA(licerun)
  CALL ACC_UNMAP_DATA(betamaxoxkappa2)
  CALL ACC_UNMAP_DATA(MLSTHG)
  CALL ACC_UNMAP_DATA(llunstr)
  CALL ACC_UNMAP_DATA(om3gmkm_gc)
  DEALLOCATE (om3gmkm_gc_d)
  CALL ACC_UNMAP_DATA(ang_gc_b)
  CALL ACC_UNMAP_DATA(lwamrsetci)
  CALL ACC_UNMAP_DATA(lwnemocoustrn)
  CALL ACC_UNMAP_DATA(nwav_gc)
  CALL ACC_UNMAP_DATA(x0tauhf)
  CALL ACC_UNMAP_DATA(isnonlin)
  CALL ACC_UNMAP_DATA(rnum)
END SUBROUTINE WAMINTGR_CUDA_GPU
