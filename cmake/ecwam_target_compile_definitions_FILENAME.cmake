function( ecwam_target_compile_definitions_FILENAME target )
    get_target_property( sources "${target}" SOURCES)
    foreach( src ${sources} )
        get_filename_component( filename ${src} NAME )
        set_property(
            SOURCE ${src}
            APPEND
            PROPERTY COMPILE_DEFINITIONS __FILENAME__="${filename}" )
    endforeach()
endfunction()

