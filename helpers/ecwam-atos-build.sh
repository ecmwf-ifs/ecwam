#!/usr/bin/env bash

# Check out dependency packages (usually only has to be done once)
./package/bundle/ecwam-bundle create --update  --bundle package/bundle/bundle.yml 

# Build the project
./package/bundle/ecwam-bundle build --arch=package/bundle/arch/ecmwf/hpc2020/intel/2021.4.0/hpcx-openmpi/2.9.0 --with-fckit