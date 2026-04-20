# (C) Copyright 2022- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

macro( ecwam_find_nemo )
  # We use this slighlty convoluted way to find nemo to ensure backwards compatibility with 50R1.
  # This can be dropped once the 50R1 maintenance window has lapsed.

  set( nemo_found 0 )
  ecbuild_find_package("nemo" QUIET)
  if( DEFINED OCEAN_PREC )
    set( _ocean_prec "" )
    string( TOLOWER "${OCEAN_PREC}" _ocean_prec )
    ecbuild_find_package("nemo_${_ocean_prec}" QUIET)
    if( nemo_${_ocean_prec}_FOUND )
        set( nemo_found 1 )
    elseif( nemo_FOUND AND nemo_HAVE_${OCEAN_PREC})
        set( nemo_found 1 )
    endif()
  else()
    if( HAVE_sp )
       ecbuild_find_package("nemo_sp" QUIET)
       if( nemo_sp_FOUND )
         set( nemo_found 1 )
       elseif( nemo_FOUND AND nemo_HAVE_SP)
         set( nemo_found 1 )
       else()
         set( nemo_found 0 )
       endif()
    endif()
    if( HAVE_dp )
       ecbuild_find_package("nemo_dp" QUIET)
       if( nemo_dp_FOUND )
         set( nemo_found 1 )
       elseif( nemo_FOUND AND nemo_HAVE_DP)
         set( nemo_found 1 )
       else()
         set( nemo_found 0 )
       endif()
    endif()
  endif()

endmacro()
