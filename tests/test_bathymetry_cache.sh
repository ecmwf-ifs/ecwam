#!/usr/bin/env bash

# Test cache behavior, not scientific bathymetry calculations.
# Usage: bash test_bathymetry_cache.sh /path/to/share/ecwam/scripts
#
# The real bathymetry wrapper and command-line parser run in an isolated
# temporary directory. Runtime settings, configuration, helpers, retrieval,
# and generation are replaced with lightweight stubs. Fake cached ETOPO1,
# ETOPO2, reference levels, and generated products are created locally.
# A test-only cache version is reused for the runtime stub, cache checksum,
# and product directory; the production version does not matter for this test.
# The real runtime is also sourced in a subshell to check dataset and reference-level
# defaults and overrides; it only sets environment variables, not downloads.
#
# For both ETOPO1 and ETOPO2, this test verifies that:
# - Normal runs reuse cached products without invoking the generator.
# - --regenerate invokes the generator and replaces bathymetry and all three
#   subgrid files without retrieving generated-product caches.
# - ECWAM_REGENERATE_BATHYMETRY=true/1 also regenerates products; false/0
#   reuses caches, invalid values fail, and --regenerate overrides false.
# - Cached ETOPO datasets and reference levels remain unchanged.
# - Disabling obstruction output creates empty subgrid placeholders.
# - Generator failure returns an error and preserves the old bathymetry.
# - Absolute dataset and reference-level paths, including spaces, are linked
#   directly without retrieval; missing absolute files fail without downloads.
#
# NO DOWNLOADS OR NETWORK REQUESTS ARE MADE BY THIS TEST, on any server.
# The real retrieval script is never invoked. Its stub only logs requests
# and checks for an existing local file; a missing file causes failure,
# never a download. The generator stub writes small text fixtures instead
# of processing ETOPO data. No real ETOPO datasets or user caches are touched.
# All fixtures and outputs are confined to the temporary directory, which
# is removed automatically when the test exits.
#
# This isolation applies to the test only. The real wrapper can download
# missing ETOPO inputs: --regenerate bypasses generated-product caches,
# not ETOPO input caching or downloads.

set -euo pipefail

scripts_dir=$1
unset ECWAM_REGENERATE_BATHYMETRY ECWAM_BATHYMETRY_VERSION
test_dir=$(mktemp -d)
trap 'status=$?; if [[ ${status} -ne 0 && -f ${test_dir}/output.log ]]; then cat "${test_dir}/output.log"; fi; rm -rf "${test_dir}"' EXIT
export TEST_DIR=${test_dir}
mkdir -p "${test_dir}/scripts" "${test_dir}/data/inputs"
default_bathymetry_version=$(
  source "${scripts_dir}/ecwam_runtime.sh"
  printf '%s' "${ECWAM_BATHYMETRY_VERSION}"
)
[[ -n ${default_bathymetry_version} ]]
test_versions=("<version-does-not-matter-for-this-test>" "<another-version-does-not-matter-for-this-test>")
(
  unset ETOPO1 ETOPO2 ETOPO1_reference_levels ETOPO2_reference_levels
  source "${scripts_dir}/ecwam_runtime.sh"
  [[ ${ECWAM_BATHYMETRY_VERSION} == "${default_bathymetry_version}" ]]
  [[ ${ETOPO1} == data/bathymetry/ETOPO1/ETOPO1_Ice_g_int.bin ]]
  [[ ${ETOPO2} == data/bathymetry/ETOPO2/etopo2_2006apr.dat ]]
  [[ ${ETOPO1_reference_levels} == data/bathymetry/ETOPO1/reference_levels.v2 ]]
  [[ ${ETOPO2_reference_levels} == data/bathymetry/ETOPO2/reference_levels.v1 ]]
  export ETOPO1=inputs/custom_etopo1
  export ETOPO2=inputs/custom_etopo2
  export ETOPO1_reference_levels=inputs/custom_etopo1_reference_levels
  export ETOPO2_reference_levels=inputs/custom_etopo2_reference_levels
  source "${scripts_dir}/ecwam_runtime.sh"
  source "${scripts_dir}/ecwam_runtime.sh"
  [[ ${ETOPO1} == inputs/custom_etopo1 ]]
  [[ ${ETOPO2} == inputs/custom_etopo2 ]]
  [[ ${ETOPO1_reference_levels} == inputs/custom_etopo1_reference_levels ]]
  [[ ${ETOPO2_reference_levels} == inputs/custom_etopo2_reference_levels ]]
)
for version in '' "${test_versions[@]}"; do
  (
    export ECWAM_BATHYMETRY_VERSION=${version}
    source "${scripts_dir}/ecwam_runtime.sh"
    source "${scripts_dir}/ecwam_runtime.sh"
    [[ ${ECWAM_BATHYMETRY_VERSION} == "${version:-${default_bathymetry_version}}" ]]
  )
done
export ECWAM_BATHYMETRY_VERSION=${test_versions[0]}
export ETOPO1=inputs/etopo1
export ETOPO2=inputs/etopo2
export ETOPO1_reference_levels=inputs/reference_levels
export ETOPO2_reference_levels=inputs/reference_levels
for script in ecwam_run_create_bathymetry.sh ecwam_parse_commandline.sh; do
  ln -s "${scripts_dir}/${script}" "${test_dir}/scripts/${script}"
done

cat > "${test_dir}/scripts/ecwam_runtime.sh" <<'EOF'
export RUN_DIR=${TEST_DIR}/run
export DATA_DIR=${TEST_DIR}/data
export ECWAM_CACHE_PATH=${DATA_DIR}
export ECWAM_PROJECT_NAME=ecwam
export ETOPO1=${ETOPO1:-inputs/etopo1}
export ETOPO2=${ETOPO2:-inputs/etopo2}
export ETOPO1_reference_levels=${ETOPO1_reference_levels:-inputs/reference_levels}
export ETOPO2_reference_levels=${ETOPO2_reference_levels:-inputs/reference_levels}
export CREATE_WAM_BATHYMETRY_ETOPO1=${TEST_DIR}/generator
export CREATE_WAM_BATHYMETRY_ETOPO2=${TEST_DIR}/generator
export LAUNCH_SERIAL=""
EOF

cat > "${test_dir}/scripts/ecwam_helper_functions.sh" <<'EOF'
function log() { ("$@" | tee -a stdout.log) 3>&1 1>&2 2>&3 | tee -a stderr.log; }
function echo() { log builtin echo -e "$@"; }
function read_config() {
  case "$1" in
    grid) builtin echo O80 ;;
    frequencies) builtin echo 29 ;;
    bathymetry) builtin echo "${TEST_BATHYMETRY:-ETOPO1}" ;;
  esac
}
function assert_executable_is_available() { command -v "$1" >/dev/null; }
function trace_ls() { log ls "$1"; }
EOF

cat > "${test_dir}/scripts/ecwam_configure.sh" <<'EOF'
function ecwam_configure() {
  xdella=1; amosop=-89; amonop=89; amowep=0; amoeap=359
  iper=1; irgg=1; fr1=0.04177248; ifre1=3; cwamresol=00080
  llobstrout=${TEST_OBSTRUCTIONS:-T}; llgribout=T
  printf 'O80 test grid\n' > grid_description
}
EOF

cat > "${test_dir}/scripts/ecwam_retrieve.sh" <<'EOF'
#!/usr/bin/env bash
printf '%s\n' "$1" >> "${TEST_DIR}/retrieved"
test -f "$2"
EOF

cat > "${test_dir}/generator-dp" <<'EOF'
#!/usr/bin/env bash
set -eu
printf 'generated\n' >> "${TEST_DIR}/generated"
test -r reference_levels
test -r ETOPO1_Ice_g_int.xyz || test -r etopo2_2006apr.dat
if [[ -n ${TEST_EXPECTED_REFERENCE_LEVELS:-} ]]; then
  [[ $(readlink -f reference_levels) == "${TEST_EXPECTED_REFERENCE_LEVELS}" ]]
  if [[ ${TEST_BATHYMETRY} == ETOPO1 ]]; then
    [[ $(readlink -f ETOPO1_Ice_g_int.xyz) == "${ETOPO1}" ]]
  else
    [[ $(readlink -f etopo2_2006apr.dat) == "${ETOPO2}" ]]
  fi
fi
if [[ ${TEST_FAIL_GENERATOR:-false} == true ]]; then
  exit 1
fi
printf 'fresh bathymetry\n' > wam_topo_00080
for index in 0 1 2; do
  printf 'fresh subgrid\n' > "wam_grib_subgrid_${index}"
done
EOF
chmod +x "${test_dir}/scripts/ecwam_retrieve.sh" "${test_dir}/generator-dp"

printf 'cached ETOPO1\n' > "${test_dir}/data/inputs/etopo1"
printf 'cached ETOPO2\n' > "${test_dir}/data/inputs/etopo2"
printf 'cached reference levels\n' > "${test_dir}/data/inputs/reference_levels"
input_checksum=$(cksum "${test_dir}"/data/inputs/*)
checksum=$(printf '%s\nO80 test grid\ncached reference levels\n' "${ECWAM_BATHYMETRY_VERSION}" | md5sum | cut -d ' ' -f 1)
products_dir=${test_dir}/data/data/bathymetry/v${ECWAM_BATHYMETRY_VERSION}
mkdir -p "${products_dir}"

function seed_products() {
  printf 'old bathymetry\n' > "${products_dir}/bathymetry_dp_O80_nfre29_${TEST_BATHYMETRY}_${checksum}"
  for index in 0 1 2; do
    printf 'old subgrid\n' > "${products_dir}/wam_grib_subgrid_dp_O80_${index}_29_${TEST_BATHYMETRY}_${checksum}"
  done
}

function run_wrapper() {
  rm -rf "${test_dir}/run"
  mkdir -p "${test_dir}/run"
  printf 'test configuration\n' > "${test_dir}/run/config.yml"
  : > "${test_dir}/retrieved"
  rm -f "${test_dir}/generated"
  bash "${test_dir}/scripts/ecwam_run_create_bathymetry.sh" "$@" > "${test_dir}/output.log" 2>&1 || {
    cat "${test_dir}/output.log"
    return 1
  }
}

for bathymetry in ETOPO1 ETOPO2; do
  export TEST_BATHYMETRY=${bathymetry}
  seed_products # Provide old cached products for the reuse checks.
  run_wrapper # Default behavior must reuse products without generating them.
  [[ ! -e ${test_dir}/generated ]]
  [[ $(grep -c '^data/bathymetry/' "${test_dir}/retrieved") -eq 4 ]]
  [[ $(cat "${test_dir}/run/wam_topo") == 'old bathymetry' ]]

  for setting in false 0; do
    export ECWAM_REGENERATE_BATHYMETRY=${setting}
    run_wrapper # Explicitly disabling regeneration must still reuse the cache.
    [[ ! -e ${test_dir}/generated ]]
    [[ $(cat "${test_dir}/run/wam_topo") == 'old bathymetry' ]]
  done
  for setting in true 1; do
    seed_products # Each enabled value must replace old products independently.
    export ECWAM_REGENERATE_BATHYMETRY=${setting}
    run_wrapper # The environment alone must force generation without cache retrieval.
    [[ $(wc -l < "${test_dir}/generated") -eq 1 ]]
    ! grep -q '^data/bathymetry/' "${test_dir}/retrieved"
    [[ $(cat "${test_dir}/run/wam_topo") == 'fresh bathymetry' ]]
    [[ $(cksum "${test_dir}"/data/inputs/*) == "${input_checksum}" ]]
  done
  export ECWAM_REGENERATE_BATHYMETRY=invalid
  # Reject invalid settings before retrieving inputs or launching the generator.
  if run_wrapper; then
    echo 'Expected invalid regeneration setting to fail' >&2
    exit 1
  fi
  [[ ! -e ${test_dir}/generated && ! -s ${test_dir}/retrieved ]]
  seed_products # Restore old products to verify that the CLI forces replacement.
  export ECWAM_REGENERATE_BATHYMETRY=false
  run_wrapper --regenerate # The CLI must override the disabled environment setting.
  [[ $(wc -l < "${test_dir}/generated") -eq 1 ]]
  ! grep -q '^data/bathymetry/' "${test_dir}/retrieved"
  [[ $(wc -l < "${test_dir}/retrieved") -eq 2 ]]
  [[ $(cat "${test_dir}/run/wam_topo") == 'fresh bathymetry' ]]
  for index in 0 1 2; do
    [[ $(cat "${test_dir}/run/wam_subgrid_${index}") == 'fresh subgrid' ]]
  done
  [[ $(cksum "${test_dir}"/data/inputs/*) == "${input_checksum}" ]]
  unset ECWAM_REGENERATE_BATHYMETRY

  export TEST_OBSTRUCTIONS=F
  run_wrapper --regenerate # Disabled obstructions must leave empty subgrid placeholders.
  for index in 0 1 2; do
    [[ -f ${test_dir}/run/wam_subgrid_${index} && ! -s ${test_dir}/run/wam_subgrid_${index} ]]
  done
  unset TEST_OBSTRUCTIONS

  seed_products # Keep known old products to check preservation after generator failure.
  export TEST_FAIL_GENERATOR=true
  # Forced regeneration must report failure without replacing the old bathymetry.
  if run_wrapper --regenerate; then
    echo 'Expected generator failure' >&2
    exit 1
  fi
  [[ $(cat "${products_dir}/bathymetry_dp_O80_nfre29_${bathymetry}_${checksum}") == 'old bathymetry' ]]
  unset TEST_FAIL_GENERATOR
done

mkdir -p "${test_dir}/absolute inputs"
cp "${test_dir}/data/inputs/etopo1" "${test_dir}/absolute inputs/ETOPO1 input"
cp "${test_dir}/data/inputs/etopo2" "${test_dir}/absolute inputs/ETOPO2 input"
printf 'custom reference levels\n' > "${test_dir}/absolute inputs/reference levels"
absolute_input_checksum=$(cksum "${test_dir}"/absolute\ inputs/*)
export ETOPO1="${test_dir}/absolute inputs/ETOPO1 input"
export ETOPO2="${test_dir}/absolute inputs/ETOPO2 input"
export ETOPO1_reference_levels="${test_dir}/absolute inputs/reference levels"
export ETOPO2_reference_levels=${ETOPO1_reference_levels}
export TEST_EXPECTED_REFERENCE_LEVELS=${ETOPO1_reference_levels}
absolute_checksum=$(printf '%s\nO80 test grid\ncustom reference levels\n' "${ECWAM_BATHYMETRY_VERSION}" | md5sum | cut -d ' ' -f 1)

for bathymetry in ETOPO1 ETOPO2; do
  export TEST_BATHYMETRY=${bathymetry}
  run_wrapper --regenerate # Absolute inputs must be linked directly without retrieval.
  [[ ! -s ${test_dir}/retrieved ]]
  [[ $(wc -l < "${test_dir}/generated") -eq 1 ]]
  [[ $(cat "${test_dir}/run/wam_topo") == 'fresh bathymetry' ]]
  [[ $(readlink "${test_dir}/run/wam_topo") == "${products_dir}/bathymetry_dp_O80_nfre29_${bathymetry}_${absolute_checksum}" ]]
  [[ $(cksum "${test_dir}"/absolute\ inputs/*) == "${absolute_input_checksum}" ]]

  if [[ ${bathymetry} == ETOPO1 ]]; then
    export ETOPO1_reference_levels=${test_dir}/missing/reference_levels
  else
    export ETOPO2_reference_levels=${test_dir}/missing/reference_levels
  fi
  # A missing absolute reference file must fail without cache or network fallback.
  if run_wrapper --regenerate; then
    echo 'Expected missing absolute reference-level file to fail' >&2
    exit 1
  fi
  [[ ! -e ${test_dir}/generated && ! -s ${test_dir}/retrieved ]]
  export ETOPO1_reference_levels=${TEST_EXPECTED_REFERENCE_LEVELS}
  export ETOPO2_reference_levels=${TEST_EXPECTED_REFERENCE_LEVELS}

  if [[ ${bathymetry} == ETOPO1 ]]; then
    export ETOPO1=${test_dir}/missing/etopo1
  else
    export ETOPO2=${test_dir}/missing/etopo2
  fi
  # A missing absolute dataset must fail before retrieval or generator execution.
  if run_wrapper --regenerate; then
    echo 'Expected missing absolute ETOPO file to fail' >&2
    exit 1
  fi
  [[ ! -e ${test_dir}/generated && ! -s ${test_dir}/retrieved ]]
  export ETOPO1="${test_dir}/absolute inputs/ETOPO1 input"
  export ETOPO2="${test_dir}/absolute inputs/ETOPO2 input"
done

for version in "${test_versions[@]}"; do
  export ECWAM_BATHYMETRY_VERSION=${version}
  export TEST_BATHYMETRY=ETOPO1
  version_checksum=$(printf '%s\nO80 test grid\ncustom reference levels\n' "${version}" | md5sum | cut -d ' ' -f 1)
  version_dir=${test_dir}/data/data/bathymetry/v${version}
  run_wrapper --regenerate
  [[ $(readlink "${test_dir}/run/wam_topo") == "${version_dir}/bathymetry_dp_O80_nfre29_ETOPO1_${version_checksum}" ]]
  for index in 0 1 2; do
    [[ $(readlink "${test_dir}/run/wam_subgrid_${index}") == "${version_dir}/wam_grib_subgrid_dp_O80_${index}_29_ETOPO1_${version_checksum}" ]]
  done
  run_wrapper
  [[ ! -e ${test_dir}/generated ]]
  [[ $(grep -c "^data/bathymetry/v${version}/" "${test_dir}/retrieved") -eq 4 ]]
done

printf 'Bathymetry cache tests passed\n'