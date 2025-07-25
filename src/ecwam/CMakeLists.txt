# (C) Copyright 2022- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

set( ecwam ${PROJECT_NAME} )

list( APPEND unwam_srcs
  unwam/output_struct.F90
  unwam/unstruct_bound.F90
  unwam/unstruct_curr.F90
  unwam/unwam.F90
  unwam/yow_rank_gloloc.F90
  unwam/yowchecksmodule.F90
  unwam/yowdatapool.F90
  unwam/yowelementpool.F90
  unwam/yowerror.F90
  unwam/yowexchangeModule.F90
  unwam/yowincludes.h
  unwam/yowmpiModule.F90
  unwam/yownodepool.F90
  unwam/yowpd.F90
  unwam/yowpdlibmain.F90
  unwam/yowrankModule.F90
  unwam/yowsidepool.F90
  unwam/yowunblkrord.F90
  unwam/yowunpool.F90
)

if( NOT TARGET ecwam_intfb )
   ecbuild_generate_fortran_interfaces(
       TARGET ecwam_intfb
       DIRECTORIES ecwam
       SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..
       DESTINATION intfb
       INCLUDE_DIRS ecwam_intfb_includes
       PARALLEL ${FCM_PARALLEL}
   )
   install( DIRECTORY ${ecwam_intfb_includes}/
            DESTINATION include/ecwam
            COMPONENT interfaces )
endif()

list( APPEND ecwam_srcs
  abort1.F90
  adjust.F90
  airsea.F90
  aki.F90
  aki_ice.F90
  alphap_tail.F90
  bouinpt.F90
  buildstress.F90
  cal_second_order_spec.F90
  cdustarz0.F90
  check.F90
  checkcfl.F90
  chkoops.F90
  chnkmin.F90
  cigetdeac.F90
  cimsstrn.F90
  ctcor.F90
  ctuw.F90
  ctuwdrv.F90
  ctuwini.F90
  ctuwupdt.F90
  current2wam.F90
  depthprpt.F90
  difdate.F90
  dominant_period.F90
  expand_string.F90
  femean.F90
  femeanws.F90
  file_transfer.F90
  findb.F90
  fkmean.F90
  fldinter.F90
  fndprt.F90
  frcutindex.F90
  gc_dispersion.h
  get_preset_wgrib_template.F90
  getbobstrct.F90
  getgrbobstrct.F90
  getcurr.F90
  getfrstwnd.F90
  getspec.F90
  getstress.F90
  getwnd.F90
  gradi.F90
  grib2wgrid.F90
  grstname.F90
  gsfile_new.F90
  h_max.F90
  halphap.F90
  headbc.F90
  imphftail.F90
  implsch.F90
  incdate.F90
  inisnonlin.F90
  init_fieldg.F90
  init_sdiss_ardh.F90
  init_x0tauhf.F90
  initdpthflds.F90
  initgc.F90
  initialint.F90
  initmdl.F90
  initnemocpl.F90
  iniwcst.F90
  intpol.F90
  intrpolchk.F90
  intspec.F90
  inwgrib.F90
  iwam_get_unit.F90
  jafu.F90
  jonswap.F90
  kerkei.F90
  kgribsize.F90
  ktoobs.F90
  kurtosis.F90
  kzeone.F90
  makegrid.F90
  mblock.F90
  mbounc.F90
  mbounf.F90
  mboxb.F90
  mchunk.F90
  mcout.F90
  means.F90
  meansqs.F90
  meansqs_gc.F90
  meansqs_lf.F90
  mfr.F90
  mfredir.F90
  mgrid.F90
  micep.F90
  mintf.F90
  mnintw.F90
  mpabort.F90
  mpbcastgrid.F90
  mpbcastintfld.F90
  mpclose_unit.F90
  mpcrtbl.F90
  mpdecomp.F90
  mpdistribfl.F90
  mpdistribscfld.F90
  mpexchng.F90
  mpfldtoifs.F90
  mpgatherbc.F90
  mpgatherfl.F90
  mpgatherscfld.F90
  mpminmaxavg.F90
  mpuserin.F90
  mstart.F90
  mswell.F90
  mubuf.F90
  mwp1.F90
  mwp2.F90
  newwind.F90
  nlweigt.F90
  notim.F90
  ns_gc.F90
  omegagc.F90
  out_onegrdpt.F90
  out_onegrdpt_sp.F90
  outbc.F90
  outbeta.F90
  outblock.F90
  outbs.F90
  outcom.F90
  outgrid.F90
  outint.F90
  outmdldcp.F90
  outnam.F90
  outpp.F90
  outsetwmask.F90
  outspec.F90
  outstep0.F90
  outwint.F90
  outwnorm.F90
  outwpsp.F90
  outwspec.F90
  packi.F90
  packr.F90
  parkind_wave.F90
  parmean.F90
  peak.F90
  peak_ang.F90
  peak_freq.F90
  peakfri.F90
  preset_wgrib_template.F90
  prewind.F90
  proenvhalo.F90
  propag_wam.F90
  propags.F90
  propags1.F90
  propags2.F90
  propconnect.F90
  propdot.F90
  readbou.F90
  readfl.F90
  readmdlconf.F90
  readpre.F90
  readsta.F90
  readstress.F90
  readwgrib.F90
  readwind.F90
  recvnemofields.F90
  rotspec.F90
  runwam.F90
  savspec.F90
  savstress.F90
  sbottom.F90
  scosfl.F90
  sdepthlim.F90
  sdissip.F90
  sdissip_ard.F90
  sdissip_jan.F90
  sdiwbk.F90
  sdice.F90
  sdice1.F90
  sdice2.F90
  sdice3.F90
  icebreak_modify_attenuation.F90
  ibrmemout.F90
  se10mean.F90
  sebtmean.F90
  second_order_lib.F90
  secondhh.F90
  secondhh_gen.F90
  secspom.F90
  semean.F90
  sep3tr.F90
  sepwisw.F90
  set_wflags.F90
  setice.F90
  setmarstype.F90
  setwavphys.F90
  sinflx.F90
  sinput.F90
  sinput_ard.F90
  sinput_jan.F90
  skewness.F90
  snonlin.F90
  spectra.F90
  spr.F90
  stat_nl.F90
  sthq.F90
  stokesdrift.F90
  stokestrn.F90
  stress_gc.F90
  stresso.F90
  strspec.F90
  tables_2nd.F90
  tabu_swellft.F90
  tau_phi_hf.F90
  taut_z0.F90
  topoar.F90
  transf.F90
  transf_bfi.F90
  transf_r.F90
  transf_snl.F90
  uibou.F90
  uiprep.F90
  unsetice.F90
  updnemofields.F90
  updnemostress.F90
  userin.F90
  vmin.F90
  vmin_d.F90
  vplus.F90
  vplus_d.F90
  w_cdf.F90
  w_cdfmax.F90
  w_maxh.F90
  w_mode_st.F90
  w_pdf.F90
  w_pmax.F90
  wam_multio_mod.F90
  wam_nproma.F90
  wam_sorti.F90
  wam_sortini.F90
  wam_u2l1cr.F90
  wam_user_clock.F90
  wamadswstar.F90
  wamcur.F90
  wamintgr.F90
  wamodel.F90
  wamwnd.F90
  wavemdl.F90
  wdfluxes.F90
  wdirspread.F90
  weflux.F90
  wgrib2fdb.F90
  wgrib_edition.F90
  wgribencode.F90
  wgribencode_model.F90
  wgribencode_values.F90
  wgribenout.F90
  wgribout.F90
  wnfluxes.F90
  wposnam.F90
  writefl.F90
  writestress.F90
  writsta.F90
  wsigstar.F90
  wsmfen.F90
  wstream_strg.F90
  wvalloc.F90
  wvchkmid.F90
  wvdealloc.F90
  wvfricvelo.F90
  wvopenbathy.F90
  wvopensubbathy.F90
  wvwamdecomp.F90
  wvwaminit.F90
  wvwaminit1.F90
  yowabort.F90
  yowaltas.F90
  yowassi.F90
  yowcard.F90
  yowcinp.F90
  yowcoer.F90
  yowconst_2nd.F90
  yowcoup.F90
  yowcout.F90
  yowcpbo.F90
  yowcurg.F90
  yowcurr.F90
  yowdes.F90
  yowfpbo.F90
  yowfred.F90
  yowgrib.F90
  yowgrib_handles.F90
  yowgribhd.F90
  yowgribinfo.F90
  yowgrid.F90
  yowgstats.F90
  yowice.F90
  yowindn.F90
  yowintp.F90
  yowjons.F90
  yowmap.F90
  yowmean.F90
  yowmespas.F90
  yowmpp.F90
  yownemoflds.F90
  yownemoio.F90
  yowparam.F90
  yowpcons.F90
  yowphys.F90
  yowprproc.F90
  yowrefd.F90
  yowshal.F90
  yowspec.F90
  yowsphere.F90
  yowstat.F90
  yowtabl.F90
  yowtemp.F90
  yowtest.F90
  yowtext.F90
  yowtrains.F90
  yowubuf.F90
  yowunit.F90
  yowwami.F90
  yowwind.F90
  yowwndg.F90
  z0wave.F90
)

# Loki GPU transformations must use the alternate driver-layer with the FIELD_API
# offload instructions
if(HAVE_LOKI AND NOT LOKI_MODE MATCHES "idem|idem-stack")
  list(APPEND ecwam_srcs wamintgr_loki_gpu.F90)
  list(REMOVE_ITEM ecwam_srcs wamintgr.F90)
  list(APPEND ecwam_srcs cireduce_loki_gpu.F90)
  list(APPEND ecwam_srcs outbs_loki_gpu.F90)
endif()

# expand derived-types using src/ecwam/yowdrvtype_config.yml
ecwam_expand_drv_types()

unset( MPI_Fortran_LIBRARIES )
if( HAVE_GPU_AWARE_MPI )
  set( MPI_Fortran_LIBRARIES MPI::MPI_Fortran )
endif()

if( HAVE_UNWAM )
  list(APPEND ecwam_srcs ${unwam_srcs})
  list(APPEND ECWAM_DEFINITIONS WAM_HAVE_UNWAM)
  if( HAVE_MPI )
      if( NOT DEFINED HAVE_MPI_F08 )
         set( HAVE_MPI_F08 1 )
      endif()
      if( DEFINED ENV{OPENMPI_VERSION} )
        set( OPENMPI_VERSION $ENV{OPENMPI_VERSION} )
        if( OPENMPI_VERSION VERSION_LESS 4.0 )
           ecbuild_warn("There are compilation problems with F08 MPI API with OpenMPI < 4.0. Defaulting to F77 API.")
           set( HAVE_MPI_F08 0 )
        endif()
      endif()
      if( HAVE_MPI_F08 )
        list(APPEND ECWAM_DEFINITIONS WAM_HAVE_MPI_F08 )
      endif()
      set( MPI_Fortran_LIBRARIES MPI::MPI_Fortran )
  else()
      set( MPI_Fortran_LIBRARIES ${MPI_SERIAL_LIBRARIES} )
  endif()
endif()

set(MULTIO_LIBRARIES)
if( HAVE_MULTIO )
  list(APPEND MULTIO_LIBRARIES multio-fapi)
  list(APPEND ECWAM_DEFINITIONS WAM_HAVE_MULTIO)
endif()

if( HAVE_OCEAN_COUPLING )
  list(APPEND ECWAM_DEFINITIONS WITH_NEMO )
endif()

if( HAVE_ECFLOW )
  list(APPEND ECWAM_PRIVATE_DEFINITIONS WAM_HAVE_ECFLOW)
endif()

if( HAVE_GPU_AWARE_MPI )
  list(APPEND ECWAM_DEFINITIONS WITH_GPU_AWARE_MPI )
endif()

foreach( prec sp dp )
  if( HAVE_${prec} )
    set( ${PNAME}_OCEANMODEL_LIBRARIES "" )
    if( HAVE_OCEAN_COUPLING )
       set( _OCEAN_PREC "")
       if( DEFINED OCEAN_PREC )
          set( _OCEAN_PREC "${OCEAN_PREC}" )
       else()
          string( TOUPPER "${prec}" _OCEAN_PREC )
       endif()
       set( ${PNAME}_OCEANMODEL_LIBRARIES nemogcmcoup.${_OCEAN_PREC} )
    endif()

    ecbuild_add_library(
        TARGET           ecwam_${prec}
        SOURCES          ${ecwam_srcs}
        PUBLIC_LIBS      fiat parkind_${prec} ecwam_intfb
                         ${MPI_Fortran_LIBRARIES}
                         ${${PNAME}_OCEANMODEL_LIBRARIES}
                         field_api_${prec}
        PRIVATE_LIBS     eccodes_f90
                         ${MULTIO_LIBRARIES}
                         ${OpenMP_Fortran_LIBRARIES}
                         $<${HAVE_ECFLOW}:ecflow_lightf>
                         $<${HAVE_ACC}:OpenACC::OpenACC_Fortran>
        PUBLIC_INCLUDES  $<INSTALL_INTERFACE:include>
        PRIVATE_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
        PRIVATE_DEFINITIONS
           ${ECWAM_PRIVATE_DEFINITIONS}
           $<${HAVE_CUDA}:WAM_HAVE_CUDA>
           $<$<STREQUAL:${_OCEAN_PREC},SP>:PARKIND1_SINGLE_NEMO>
           $<$<STREQUAL:${prec},sp>:WAM_HAVE_SINGLE_PRECISION>
        PUBLIC_DEFINITIONS ${ECWAM_DEFINITIONS}
    )

    ecwam_target_fortran_module_directory(
        TARGET ecwam_${prec}
        MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module/ecwam_${prec}
        INSTALL_DIRECTORY module/ecwam_${prec}
    )

    if( HAVE_ACC AND CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC")
       target_link_options( ecwam_${prec} INTERFACE SHELL:${OpenACC_Fortran_FLAGS} )
       if( CMAKE_CUDA_ARCHITECTURES MATCHES 80 )
         target_compile_options( ecwam_${prec} PRIVATE "-gpu=maxregcount:168,fastmath" )
       else()
         target_compile_options( ecwam_${prec} PRIVATE "-gpu=fastmath" )
       endif()
    endif()

    # add autopromotion of reals if we are compiling in double precision
    if( prec MATCHES "dp" )
      target_compile_options( ecwam_${prec} PRIVATE ${autopromote_flags} )
    endif()

    if( HAVE_CUDA )
      target_link_options( ecwam_${prec} PUBLIC "-cuda;-gpu=pinned" )
    endif()

    ecwam_target_compile_definitions_FILENAME( ecwam_${prec} )

  endif()

endforeach()

### The file mubuf.F90, which is only used for "preproc" is sensitive to optimisations
#   possibly leading to different wam_grid_<1,2,3> files.
#   This in turn leads to non-neglibible differences
#   of average 'swh' when running "chief".

if( CMAKE_Fortran_COMPILER_ID MATCHES Intel )
  set_source_files_properties( mubuf.F90 PROPERTIES COMPILE_OPTIONS "-fp-model;strict" )
  set_source_files_properties( propconnect.F90 PROPERTIES COMPILE_OPTIONS "-fp-model;strict" )
  set_source_files_properties( stresso.F90 PROPERTIES COMPILE_OPTIONS "-fp-model;strict" )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES GNU )
  set_source_files_properties( mubuf.F90 PROPERTIES COMPILE_OPTIONS "-ffp-contract=off" )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC" AND CMAKE_BUILD_TYPE MATCHES "Bit")
  set_source_files_properties(
      w_maxh.F90 sbottom.F90 PROPERTIES COMPILE_FLAGS " -g -O1 -Mflushz -Mno-signed-zeros "
  )
  set_source_files_properties( mubuf.F90 PROPERTIES COMPILE_OPTIONS "-Mnofma" )
  set_source_files_properties( initgc.F90 PROPERTIES COMPILE_FLAGS " -g -O1 -Mflushz -Mno-signed-zeros " )
  set_source_files_properties( iniwcst.F90 PROPERTIES COMPILE_FLAGS " -g -O1 -Mflushz -Mno-signed-zeros " )
  set_source_files_properties( depthprpt.F90 PROPERTIES COMPILE_FLAGS " -g -O1 -Mflushz -Mno-signed-zeros " )
  if( HAVE_SINGLE_PRECISION )
     set_source_files_properties( aki.F90 PROPERTIES COMPILE_FLAGS " -g -O1 -Mflushz -Mno-signed-zeros " )
     set_source_files_properties( kurtosis.F90 PROPERTIES COMPILE_FLAGS " -g -O1 -Mflushz -Mno-signed-zeros " )
     set_source_files_properties( stat_nl.F90 PROPERTIES COMPILE_FLAGS " -g -O1 -Mflushz -Mno-signed-zeros " )
     set_source_files_properties( transf_bfi.F90 PROPERTIES COMPILE_FLAGS " -g -O1 -Mflushz -Mno-signed-zeros " )
  endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC" AND CMAKE_BUILD_TYPE MATCHES "Debug")
  string(REPLACE "-Ktrap=fp" "" ${PNAME}_Fortran_FLAGS_DEBUG ${${PNAME}_Fortran_FLAGS_DEBUG})
  set_source_files_properties( outbeta.F90 PROPERTIES COMPILE_OPTIONS "${${PNAME}_Fortran_FLAGS_DEBUG} -Ktrap=divz")
  set_source_files_properties( secondhh.F90 PROPERTIES COMPILE_OPTIONS "${${PNAME}_Fortran_FLAGS_DEBUG} -Ktrap=inv,ovf")
endif()

### The file grib2wgrid.F90 is sensitive to optimizations in single precision builds.
#   This leads to non-neglibible differences
#   of average 'swh' when running "chief".

if( CMAKE_Fortran_COMPILER_ID MATCHES Intel AND HAVE_SINGLE_PRECISION )
  set_source_files_properties( grib2wgrid.F90 PROPERTIES COMPILE_OPTIONS "-fp-model;strict" )
endif()

if( HAVE_LOKI )

  set( LOKI_FRONTEND "fp" CACHE STRING "Frontend parser for Loki source transformations" )

  if( LOKI_MODE MATCHES "idem|idem-stack" )
     set( LOKI_CONFIG_FILE ${CMAKE_CURRENT_SOURCE_DIR}/ecwam_loki.config )
  else()
     set( LOKI_CONFIG_FILE ${CMAKE_CURRENT_SOURCE_DIR}/ecwam_loki_gpu.config )
  endif()

  foreach(prec sp dp )
     if( HAVE_${prec} )

        if( NOT LOKI_MODE MATCHES "idem|idem-stack" )
           target_compile_definitions( ecwam_${prec} PRIVATE WAM_GPU )
        endif()

        if( NOT LOKI_MODE MATCHES "idem|idem-stack")
          # Preprocess propags2.F90 for GPU enabled runs using Loki
          loki_transform_target( TARGET ecwam_${prec}
            MODE "loop-unroll-fuse"
            FRONTEND ${LOKI_FRONTEND}
            CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/ecwam_propags2_loki.config
            PLAN ${CMAKE_CURRENT_BINARY_DIR}/loki_propags2_plan_ecwam.cmake
            SOURCES
              ${CMAKE_CURRENT_SOURCE_DIR}/
          )
        endif()

        # Apply Loki source file transformation to lib target
        loki_transform_target( TARGET ecwam_${prec}
          MODE ${LOKI_MODE}
          FRONTEND ${LOKI_FRONTEND}
          CONFIG ${LOKI_CONFIG_FILE}
          PLAN ${CMAKE_CURRENT_BINARY_DIR}/loki_plan_ecwam.cmake
          CPP
          DEFINITIONS WAM_GPU
          INCLUDES ${ecwam_intfb_includes}
          SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/
        )
     endif()
  endforeach()

endif()
