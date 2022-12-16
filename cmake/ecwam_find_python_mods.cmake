# (C) Copyright 2020- ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

macro( ecwam_find_python_mods )
   ecbuild_find_python()
   find_program( FYPP_LOC fypp REQUIRED )
   if( FYPP_LOC )
     ecbuild_info( "${ECWAM_PROJECT_NAME} FOUND fypp" )
   endif()
   execute_process(
       COMMAND python3 -c "import yaml"
       RESULT_VARIABLE EXIT_CODE
       OUTPUT_QUIET
   )
   if( EXIT_CODE EQUAL 0 )
       ecbuild_info("${ECWAM_PROJECT_NAME} FOUND Python interpreter and required modules")
   else()
       ecbuild_critical("PyYAML is needed to build ${ECWAM_PROJECT_NAME} with field_api")
   endif()
endmacro()
