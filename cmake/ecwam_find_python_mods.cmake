# (C) Copyright 2020- ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

macro( ecwam_find_python_mods )

   # Find the system python install
   set( Python3_FIND_VIRTUALENV STANDARD )

   # Look for python interpreter and pyyaml
   find_package( Python3 COMPONENTS Interpreter REQUIRED)
   set(ECWAM_PYTHON_INTERP ${Python3_EXECUTABLE})

   set( pyyaml_FOUND OFF )
   # Look for pyyaml
   execute_process(
       COMMAND ${ECWAM_PYTHON_INTERP} -c "import yaml"
       RESULT_VARIABLE EXIT_CODE
       OUTPUT_QUIET ERROR_QUIET
   )
   if( EXIT_CODE EQUAL 0 )
     ecbuild_info("${ECWAM_PROJECT_NAME} FOUND pyyaml")

     set( pyyaml_FOUND ON )
   endif()

   set( fypp_FOUND OFF )
   ecbuild_find_package( fckit )

   # Look for (non-fckit) fypp pre-processor
   execute_process(
       COMMAND ${ECWAM_PYTHON_INTERP} -c "import fypp"
       RESULT_VARIABLE EXIT_CODE
       OUTPUT_QUIET ERROR_QUIET
   )
   if( EXIT_CODE EQUAL 0 )
     ecbuild_info("${ECWAM_PROJECT_NAME} FOUND fypp")

     set( fypp_FOUND ON )
     set( FYPP ${ECWAM_PYTHON_INTERP} -m fypp )
   endif()

   if( NOT pyyaml_FOUND OR NOT fypp_FOUND )
     if( NOT fckit_HAVE_FCKIT_VENV )
        ecbuild_critical("${ECWAM_PROJECT_NAME}: fckit python virtual environment was not installed")
     endif()

     ecbuild_info("${ECWAM_PROJECT_NAME} did not find either fypp or pyyaml, therefore using fckit python virtual environment")

     set( FYPP ${FCKIT_VENV_EXE} -m fypp )
     set(ECWAM_PYTHON_INTERP ${FCKIT_VENV_EXE})
   endif()

endmacro()
