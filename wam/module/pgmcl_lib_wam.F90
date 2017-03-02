#if defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
#else
MODULE PGMCL_LIB_WAM
  implicit none
  real, dimension(1) :: US_coupl, Z0_coupl
  logical :: HAVE_NEW_COUPLING_FIELDS

  contains
     SUBROUTINE WAV_initialize_all_coupling()
        write(*,*) ' calling WAV_initialize_all_coupling' 
        write(*,*) '!!!!!!!!!!!! this is a dummy version of it !!!'
        call abort
     END SUBROUTINE WAV_initialize_all_coupling
     SUBROUTINE WAV_assign_var()
        write(*,*) ' calling WAV_assign_var' 
        write(*,*) '!!!!!!!!!!!! this is a dummy version of it !!!'
        call abort
     END SUBROUTINE WAV_assign_var
     SUBROUTINE WAV_all_import_export()
        write(*,*) ' calling  WAV_all_import_export'
        write(*,*) '!!!!!!!!!!!! this is a dummy version of it !!!'
        call abort
     END SUBROUTINE  WAV_all_import_export
     SUBROUTINE WAV_coupl_prewind()
        write(*,*) ' calling  WAV_coupl_prewind'
        write(*,*) '!!!!!!!!!!!! this is a dummy version of it !!!'
        call abort
     END SUBROUTINE WAV_coupl_prewind
END MODULE PGMCL_LIB_WAM
#endif
