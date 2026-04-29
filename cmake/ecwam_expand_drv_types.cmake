# (C) Copyright 2024- ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


# Read src/ecwam/yowdrvtype_config.yml and expand the derived-types accordingly.
# We create one sourcefile per derived-type because of *very* slow compilation times
# with nvfortran if they are all placed in the same module. Unfortunately the problem
# still persists even if compiler optimisations are disabled.

macro( ecwam_expand_drv_types )

   if( HAVE_LOKI AND NOT LOKI_MODE MATCHES "idem|idem-stack" )
      list(APPEND FYPP_ARGS -DWAM_GPU)
   endif()

   execute_process(
       COMMAND ${ECWAM_PYTHON_INTERP} -c
       "import sys; sys.path.append('${CMAKE_CURRENT_SOURCE_DIR}/../../share/ecwam/scripts'); \
       from ecwam_yaml_reader import yaml; f = open('${CMAKE_CURRENT_SOURCE_DIR}/yowdrvtype_config.yml'); \
       yml = f.read(); f.close(); objtypes = yaml.safe_load(yml)['objtypes']; print(list(objtypes))"
       RESULT_VARIABLE EXIT_CODE
       OUTPUT_VARIABLE TYPE_NAMES
       ERROR_QUIET
   )

   if( NOT EXIT_CODE EQUAL 0 )
      ecbuild_critical("${ECWAM_PROJECT_NAME} FAILED TO READ yowdrvtype_config.yml")
   endif()

   string(REGEX MATCHALL "\'[A-Za-z0-9_]+\'" type_names "${TYPE_NAMES}")

   foreach(type IN LISTS type_names)
      string(REPLACE "'" "" _type ${type})
      add_custom_command(
          OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_type}_type_mod.F90
          COMMAND ${FYPP} ${FYPP_ARGS} -m io -m os -DTYPE_NAME='${_type}'
                  -M ${CMAKE_CURRENT_SOURCE_DIR}/../../share/ecwam/scripts -m ecwam_yaml_reader
                  ${CMAKE_CURRENT_SOURCE_DIR}/drvtype_mod.fypp > ${_type}_type_mod.F90
          DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/drvtype_mod.fypp
          VERBATIM)
      list( APPEND ecwam_srcs ${CMAKE_CURRENT_BINARY_DIR}/${_type}_type_mod.F90)
   endforeach()
   
   add_custom_command(
       OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/yowdrvtype.F90
       COMMAND ${FYPP} -m io -m os -M ${CMAKE_CURRENT_SOURCE_DIR}/../../share/ecwam/scripts -m ecwam_yaml_reader
               ${CMAKE_CURRENT_SOURCE_DIR}/yowdrvtype.fypp > yowdrvtype.F90
       DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/yowdrvtype.fypp
       VERBATIM)
   list( APPEND ecwam_srcs ${CMAKE_CURRENT_BINARY_DIR}/yowdrvtype.F90)

endmacro()
