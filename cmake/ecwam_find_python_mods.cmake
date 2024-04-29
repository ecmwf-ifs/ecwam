# (C) Copyright 2020- ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

macro( ecwam_find_python_mods )
   # Look for fypp pre-processor
   find_program( FYPP fypp HINTS ${fypp_ROOT} )
   if( FYPP )
     ecbuild_info( "${ECWAM_PROJECT_NAME} FOUND fypp" )
   else()
     ecbuild_critical( "${ECWAM_PROJECT_NAME} FAILED to find required package fypp" )
   endif()
   # We do a QUIET ecbuild_find_package to update the ecbuild project summary
   ecbuild_find_package( fypp QUIET )

   # Look for python interpreter and ruamel package
   find_package( Python3 COMPONENTS Interpreter )
   execute_process(
       COMMAND python3 -c "import ruamel.yaml"
       RESULT_VARIABLE EXIT_CODE
       OUTPUT_QUIET
   )
   if( EXIT_CODE EQUAL 0 )
     ecbuild_info("${ECWAM_PROJECT_NAME} FOUND ruamel.yaml")
   else()
     ecbuild_critical( "${ECWAM_PROJECT_NAME} FAILED to find required package ruamel.yaml" )
   endif()
   # We do a QUIET ecbuild_find_package to update the ecbuild project summary
   ecbuild_find_package( ruamel QUIET)
endmacro()
