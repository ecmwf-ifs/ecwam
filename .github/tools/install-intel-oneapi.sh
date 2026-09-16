#!/usr/bin/env bash

set -euo pipefail

compiler=intel-classic
classic_version=2023.2.0

while [ "$#" -gt 0 ]; do
    case "$1" in
        --compiler)
            compiler="$2"
            shift
            ;;
        *)
            echo "Unrecognized argument '$1'" >&2
            exit 1
            ;;
    esac
    shift
done

case "${compiler}" in
    intel-classic|intel-modern)
        ;;
    *)
        echo "Unknown Intel compiler '${compiler}'" >&2
        exit 1
        ;;
esac

if [ "$(id -u)" -eq 0 ]; then
    sudo=()
else
    sudo=(sudo)
fi

keyring=/usr/share/keyrings/intel-oneapi-keyring.asc
wget -qO intel-oneapi-keyring.asc \
    https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
"${sudo[@]}" install -m 0644 intel-oneapi-keyring.asc "${keyring}"
rm intel-oneapi-keyring.asc
echo "deb [signed-by=${keyring}] https://apt.repos.intel.com/oneapi all main" | \
    "${sudo[@]}" tee /etc/apt/sources.list.d/oneAPI.list >/dev/null

"${sudo[@]}" apt-get update

if [ "${compiler}" = intel-modern ]; then
    "${sudo[@]}" apt-get install -y intel-oneapi-toolkit
else
    "${sudo[@]}" apt-get install -y \
        "intel-oneapi-compiler-fortran-${classic_version}" \
        "intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic-${classic_version}" \
        intel-oneapi-mpi-devel-2021.10.0 \
        "intel-oneapi-mkl-${classic_version}"
fi
