# (C) Copyright 2020- ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

macro( ecwam_find_python_mods )

   # Look for fypp pre-processor
   find_program( FYPP fypp )
   if( fypp_FOUND )
     ecbuild_info( "${ECWAM_PROJECT_NAME} FOUND fypp" )

     # We do a QUIET ecbuild_find_package to update the ecbuild project summary
     ecbuild_find_package( fypp QUIET )
   endif()

   set( yaml_FOUND OFF )

   # Look for python interpreter and pyyaml
   find_package( Python3 COMPONENTS Interpreter REQUIRED)
   execute_process(
       COMMAND python3 -c "import yaml"
       RESULT_VARIABLE EXIT_CODE
       OUTPUT_QUIET ERROR_QUIET
   )
   if( EXIT_CODE EQUAL 0 )
     ecbuild_info("${ECWAM_PROJECT_NAME} FOUND pyyaml")
     set( yaml_FOUND ON)

     # We do a QUIET ecbuild_find_package to update the ecbuild project summary
     ecbuild_find_package( pyyaml QUIET)
   endif()

   if( NOT yaml_FOUND OR NOT fypp_FOUND )
     ecbuild_find_package( fckit REQUIRED )
   endif()

endmacro()
