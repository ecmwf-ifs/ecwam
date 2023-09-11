# (C) Copyright 2020- ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


ecbuild_find_package(NAME field_api QUIET)
set(clone_field_api TRUE)

## Test if field_api was cloned by ecWAM
if(field_api_FOUND)
   cmake_path(RELATIVE_PATH field_api_DIR BASE_DIRECTORY ${CMAKE_BINARY_DIR} OUTPUT_VARIABLE path_var)
   cmake_path(GET path_var PARENT_PATH parent_path)
   string(FIND ${parent_path} "../" result_var)

#  If field_api is found but was not cloned by ecWAM, clone_field_api is set to FALSE
   string(COMPARE EQUAL ${result_var} "-1" clone_field_api)
endif()

if( clone_field_api )
   include(FetchContent)
   FetchContent_Declare(
      field_api
      GIT_REPOSITORY git@github.com:ecmwf-ifs/field_api.git
      GIT_TAG v0.2.1
      OVERRIDE_FIND_PACKAGE
   )
   
   set(FIELD_API_ENABLE_TESTS OFF)
   set(FIELD_API_ENABLE_ACC OFF)
   FetchContent_MakeAvailable(field_api)
endif()
