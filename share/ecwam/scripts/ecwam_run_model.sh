#!/usr/bin/env bash

# (C) Copyright 2021- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

set -o pipefail
set -e
set +v

SCRIPTS_DIR="$( cd $( dirname "${BASH_SOURCE[0]}" ) && pwd -P )"
if [[ $(basename ${SCRIPTS_DIR}) == "bin" ]]; then
  SCRIPTS_DIR="$( cd $( dirname "${SCRIPTS_DIR}/$(readlink "${BASH_SOURCE[0]}")" ) && pwd -P )"
fi

ECWAM_CONTEXT=model
source ${SCRIPTS_DIR}/ecwam_runtime.sh
source ${SCRIPTS_DIR}/ecwam_parse_commandline.sh
source ${SCRIPTS_DIR}/ecwam_helper_functions.sh

assert_executable_is_available ${MODEL}-${prec} || abort 4

begofrn=$(read_config begin               --format="%Y%m%d%H%M%S")
endofrn=$(read_config end                 --format="%Y%m%d%H%M%S")
begoffo=$(read_config forcing.at[1].begin --format="%Y%m%d%H%M%S" --default=${endofrn} )

find_preproc_files

lgribin=$(read_config lgribin --default=F)
lgribout=$(read_config lgribout --default=F)
# Note that lnocdin=T is only meaningful for grib restart (lgribin=T) (rather than pure binary)
# for which the model drag coefficient and the corresponding 10m wind speed are not available 
# For cycling the runs, one will need to extract the relevant parameters from the grib output MPP*
# and place them in the relevant CDWAVEIN* and UWAVEIN* files 
lnocdin=T


if [[ ${lgribin} = T ]]; then
  if [[ ${lnocdin} = T ]]; then
    initial_file_header_list="SGS"
  else
    initial_file_header_list="SGS CDWAVEIN UWAVEIN"
  fi
else
  initial_file_header_list="LAW BLS"
fi
find_preset_files ${begofrn} ${initial_file_header_list}

output_fields=$(read_config output.fields.name[:] --default="")
output_statistics=$(read_config output.statistics.name[:] --default="")

output_parameters_index=""
[[ -z "${output_fields}" ]] || output_parameters_index=$(${ECWAM_PYTHON_INTERP} ${SCRIPTS_DIR}/ecwam_parameters.py --request=index ${output_fields})

OUTPUT_FLAGS=""
if [ $(read_config output.fields.format --default=grib) = grib ] ; then
  for index in ${output_parameters_index}; do
    [[ -z "${OUTPUT_FLAGS}" ]] || OUTPUT_FLAGS+=$'\n  '
    OUTPUT_FLAGS+="GFLAG(${index})=T,"
  done
else
  for index in ${output_parameters_index}; do
    [[ -z "${OUTPUT_FLAGS}" ]] || OUTPUT_FLAGS+=$'\n  '
    OUTPUT_FLAGS+="FFLAG(${index})=T,"
  done
fi

output_parameters_index=""
[[ -z "${output_statistics}" ]] || output_parameters_index=$(${ECWAM_PYTHON_INTERP} ${SCRIPTS_DIR}/ecwam_parameters.py --request=index ${output_statistics})
for index in ${output_parameters_index}; do
    [[ -z "${OUTPUT_FLAGS}" ]] || OUTPUT_FLAGS+=$'\n  '
    OUTPUT_FLAGS+="NFLAG(${index})=T,"
done

NAWI=""
for i in $(seq 1 $(read_config forcings.at.size)); do
  [[ -z "${NAWI}" ]] || NAWI+=$'\n'
  forcing_end=$(read_config forcings.at[${i}-1].end --format=%Y%m%d%H%M%S)
  forcing_step=$(read_config forcings.at[${i}-1].timestep --format=seconds)
  NAWI+="&NAWI IDWI=${forcing_step}, IDWO=${forcing_step}, CLWOUT=\"${forcing_end}\" /"
done

NAOS=""
for restart_time in $(read_config output.restart.at[:].time --format=%Y%m%d%H%M%S --default=" "); do
  [[ -z "${NAOS}" ]] || NAOS+=$'\n'
  NAOS+="&NAOS CLSOUT=\"${restart_time}\" /"
done

# MODEL SETUP:
##############
cldomain=$(read_config cldomain --default=g)
wamnang=$(read_config directions)
wamnfre=$(read_config frequencies)
fr1=$(read_config fr1 --default=4.177248E-02)
ifre1=$(read_config ifre1 --default=1)
lsubgrid=$(read_config lsubgrid --default=T)
irefra=$(read_config advection.irefra --default=0)

iphys=$(read_config physics.iphys --default=1)
llgcbz0=$(read_config physics.llgcbz0 --default=F)
llnormagam=$(read_config physics.llnormagam --default=F)
lmaskice=$(read_config physics.lmaskice --default=T)
lciwa1=$(read_config physics.lciwa1 --default=F)
lciwa2=$(read_config physics.lciwa2 --default=F)
lciwa3=$(read_config physics.lciwa3 --default=F)
lciscal=$(read_config physics.lciscal --default=F)

nproma=$(read_config nproma --default=24)

# read timesteps
phys_tstp=$(read_config physics.timestep --format=seconds --default=900)
adv_base_tstp=$(read_config advection.timestep --format=seconds --default=900)
adv_fast_tstp=$(read_config advection.fast_waves.timestep --format=seconds --default=$adv_base_tstp)
ifrelfmax=$(read_config advection.fast_waves.max_frequency --default=0)
idelcur=$(read_config currents.input_step --default=86400)

# verify timesteps
if [ $(( $adv_base_tstp%$adv_fast_tstp )) -ne 0 ] ; then
   echo "ERROR: Base advection timestep should be a multiple of fast-wave advection timestep"
   exit 4
fi
if [ $(( $phys_tstp%$adv_base_tstp )) -ne 0 ] ; then
   echo "ERROR: Physics timestep should be a multiple of base advection timestep"
   exit 4
fi

if [[ $(read_config forcings.sea_ice --default=True) == "True" ]] ; then
  licerun=T
else
  licerun=F
fi

ppfreq=$(read_config output.fields.at[0].timestep --format=hours --default=02:00)

forcings_file=$(read_config forcings.file)
currents_file=$(read_config currents.file --default=NOTFOUND )

function cleanup() {
  shopt -s nullglob
  LOG_DIR=${RUN_DIR}/logs/model
  mkdir -p ${LOG_DIR}
  rm -rf ${LOG_DIR}/*
  mv ${WORK_DIR}/stdout.log ${LOG_DIR}/
  mv ${WORK_DIR}/stderr.log ${LOG_DIR}/
  [[ -f ${WORK_DIR}/meminfo.txt    ]] && mv ${WORK_DIR}/meminfo.txt    ${LOG_DIR}/
  [[ -f ${WORK_DIR}/gstats.xml     ]] && mv ${WORK_DIR}/gstats.xml     ${LOG_DIR}/
  [[ -f ${WORK_DIR}/wam.log        ]] && mv ${WORK_DIR}/wam.log        ${LOG_DIR}/
  [[ -f ${WORK_DIR}/wam.log.1      ]] && mv ${WORK_DIR}/wam.log.*      ${LOG_DIR}/
  [[ -f ${WORK_DIR}/drhook.prof.1  ]] && mv ${WORK_DIR}/drhook*        ${LOG_DIR}/
  [[ -f ${WORK_DIR}/gstats.xml     ]] && mv ${WORK_DIR}/gstats.xml     ${LOG_DIR}/
  [[ -f ${WORK_DIR}/gstats.1       ]] && mv ${WORK_DIR}/gstats*        ${LOG_DIR}/
  [[ -f ${WORK_DIR}/gstats         ]] && mv ${WORK_DIR}/gstats         ${LOG_DIR}/
  [[ -f ${WORK_DIR}/statistics.log ]] && mv ${WORK_DIR}/statistics.log ${LOG_DIR}/

  # Delete symbolic links
  find ${WORK_DIR} -maxdepth 1 -type l -delete

  mkdir -p ${RUN_DIR}/output/
  for outfile in ${WORK_DIR}/MPP*; do
    mv ${outfile} ${RUN_DIR}/output/
  done

  mkdir -p ${RUN_DIR}/restart/
  for outfile in ${WORK_DIR}/LAW* ${WORK_DIR}/BLS* ${WORK_DIR}/SGS* ; do
    mv ${outfile} ${RUN_DIR}/restart/
  done

  mv ${WORK_DIR}/wam_namelist ${RUN_DIR}/

  unset -f echo
  echo
  echo "Log     files have been copied to    ${LOG_DIR}"
  echo "Output  files have been copied to    ${RUN_DIR}/output/"
  echo "Restart files have been copied to    ${RUN_DIR}/restart/"

  cd ${RUN_DIR}
  rm -rf ${WORK_DIR}
  rm workdir # symlink

  shopt -u nullglob
}

workdir model
ln -sf ../config.yml config.yml
echo "*******************************************************************************"
log ${ECWAM_PROJECT_NAME} --info
echo "*******************************************************************************"

# General grid information and model constants
##############################################

for preproc_file in wam_grid_tables wam_subgrid_0 wam_subgrid_1 wam_subgrid_2; do
  if [[ -r ${RUN_DIR}/${preproc_file} ]] ; then
    ln -s  ${RUN_DIR}/${preproc_file} ${preproc_file}
  else
    echo "\n\n\t\tPREPROC file ${preproc_file} was not found in ${RUN_DIR}\n\n"
    trace_ls ${RUN_DIR}
    exit 1
  fi
done

# Initial stress and spectrum
#############################

rm -rf specwavein cdwavein uwavein

for header in ${initial_file_header_list}; do
  preset_file=${header}${begofrn}_000000000000

  if [[ -r ${RUN_DIR}/restart/${preset_file} ]] ; then
    if [[ ${lgribin} = T ]]; then
      # Grib restart files
      if [[ ${header} = SGS ]] ; then
        ln -sf ${RUN_DIR}/restart/${preset_file} specwavein
      elif [[ ${header} = CDWAVEIN ]] ; then
        ln -sf ${RUN_DIR}/restart/${preset_file} cdwavein 
      elif [[ ${header} = UWAVEIN ]] ; then
        ln -sf ${RUN_DIR}/restart/${preset_file} uwavein 
      else
        ln -sf ${RUN_DIR}/restart/${preset_file} ${preset_file}
      fi
    else
      # Pure binary restart files
      ln -sf ${RUN_DIR}/restart/${preset_file} ${preset_file}
    fi
  else
    echo "\n\n\tPRESET FILE ${preset_file} was not found in ${RUN_DIR}/restart"
    trace_ls ${RUN_DIR}/restart
    exit 1
  fi

done

# Forcing
#########

log ${SCRIPTS_DIR}/ecwam_retrieve.sh ${forcings_file} ${DATA_DIR}/${forcings_file} || {
  echo "ERROR: Could not retrieve forcings ${forcings_file}"
  exit 4
}

ln -s ${DATA_DIR}/${forcings_file} sfcwindin


# Surface currents (if needed)
##############################
if [ "${currents_file}" != NOTFOUND ]; then
  log ${SCRIPTS_DIR}/ecwam_retrieve.sh ${currents_file} ${DATA_DIR}/${currents_file} || {
    echo "ERROR: Could not retrieve surface currents ${currents_file}"
    exit 4
}

  ln -s ${DATA_DIR}/${currents_file} currents 
fi

# NAMELIST
##########

cat > wam_namelist << EOF
&NALINE
  CLHEADER              = " WAVE MODEL ",
  CLDOMAIN              = "${cldomain}",
  NPROMA_WAM            = ${nproma},
  LL1D                  = F,
  NANG                  = ${wamnang},
  NFRE                  = 36,
  NFRE_RED              = ${wamnfre},
  FR1                   = ${fr1},
  IFRE1                 = ${ifre1},
  CBPLTDT               = "${begofrn}",
  CEPLTDT               = "${endofrn}",
  CDATEF                = "${begoffo}",
  CDATECURA             = "${begofrn}",
  DELPRO_LF             = ${adv_fast_tstp},
  IFRELFMAX             = ${ifrelfmax},
  IDELPRO               = ${adv_base_tstp},
  IDELT                 = ${phys_tstp},
  IDELCUR               = ${idelcur}
  IPHYS                 = ${iphys},
  LLGCBZ0               = ${llgcbz0},
  LLNORMAGAM            = ${llnormagam},
  ISNONLIN              = 0,
  LICERUN               = ${licerun},
  LMASKICE              = ${lmaskice},
  LCIWA1                = ${lciwa1},
  LCIWA2                = ${lciwa2},
  LCIWA3                = ${lciwa3},
  LICETH                = F,
  ZALPFACB              = 1.0,
  LCISCAL               = ${lciscal},
  IPROPAGS              = 2,
  LSUBGRID              = ${lsubgrid},
  IREFRA                = ${irefra},
  LRELWIND              = T,
  RWFAC                 = 0.5,
  LFDBIOOUT             = F,
  LFDB                  = F,
  NGRIB_VERSION         = 2,
  LLRSTGRIBPARAM        = T,
  YCLASS                = "rd",
  YEXPVER               = "wave",
  ISTREAM               = 1045,
  CPATH                 = "${WORK_DIR}",
  LGRIBIN               = ${lgribin},
  LNOCDIN               = ${lnocdin},
  LGRIBOUT              = ${lgribout},
  LFRSTFLD              = T,
  IDELRES               = 0,  ! regular output for restart spectra, ignored if NAOS sections exists
  IDELINT               = ${ppfreq},
  LSECONDORDER          = F,
  IREST                 = 1,
  LRSTPARALR            = F,
  LRSTPARALW            = F,
  LLNORMWAMOUT          = T,
  LLNORMWAMOUT_GLOBAL   = T,
  CNORMWAMOUT_FILE      = "statistics.log",
  ${OUTPUT_FLAGS}
/
${NAWI}
${NAOS}
EOF

trace_ls $(pwd)

echo "*******************************************************************************"
echo "NAMELIST INPUT: wam_namelist"
echo "*******************************************************************************"
log cat wam_namelist
echo

precision="double"
if [ "${prec}" = "sp" ]; then
    precision="single"
fi;
if ${dryrun}; then
  echo "*******************************************************************************"
  echo "WAVE MODEL DRYRUN"
  echo "*******************************************************************************"
  echo "#1 Go to workdir\n"
  echo "  cd ${RUN_DIR}/workdir"
  echo
  echo "#2 Execute wave model\n"
  echo "  $(abs_path ${MODEL}-${prec})"
  echo 
  echo "#3 Validate results\n"
  echo "  ${SCRIPTS_DIR}/ecwam_validation.py config.yml statistics.log --section=validation.${precision}_precision"
  echo
  echo "#4 Clean up workdir when ready"
  echo "  ${RUN_DIR}/workdir/cleanup.sh"
  exit 0
fi

echo "*******************************************************************************"
echo "WAVE MODEL START"
echo "\n+ ${LAUNCH_PARALLEL} $(which ${MODEL}-${prec})\n"
echo "*******************************************************************************"


START=$(date +%s)
log ${LAUNCH_PARALLEL} $(which ${MODEL}-${prec}) || {
  sleep 1
  echo
  echo "*******************************************************************************"
  echo "WAVE MODEL FAILED EXECUTION"
  echo "*******************************************************************************"
  trace_ls $(pwd)
  cleanup
  exit 14
}

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "\n\n\t Running ${LAUNCH_PARALLEL} $(which ${MODEL}-${prec}) took $DIFF seconds\n"

trace_ls $(pwd)

if [ "$(read_config validation.${precision}_precision --default=NOTFOUND)" != NOTFOUND ]; then
  echo "\n\n"
  echo "*******************************************************************************"
  echo "Validation"
  echo "*******************************************************************************"
  mkdir -p ${RUN_DIR}/logs/model
  cp statistics.log ${RUN_DIR}/logs/model/statistics.log
  echo " + ${ECWAM_PYTHON_INTERP} ${SCRIPTS_DIR}/ecwam_validation.py ${RUN_DIR}/config.yml ${RUN_DIR}/logs/model/statistics.log --section=validation.${precision}_precision"
  log ${ECWAM_PYTHON_INTERP} ${SCRIPTS_DIR}/ecwam_validation.py ${RUN_DIR}/config.yml ${RUN_DIR}/logs/model/statistics.log --section=validation.${precision}_precision || {
    cleanup
    exit 1
  }
fi

cleanup

touch ${RUN_DIR}/model.completed
