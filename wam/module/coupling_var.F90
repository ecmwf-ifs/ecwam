#if defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
#else
MODULE COUPLING_VAR 
  implicit none
  integer :: WAV_COMM_WORLD 
  real*8 :: WAV_PresTime, CouplingTime_ATM_WAV, ATM_EXCHANGE, TEST_DIV
END MODULE COUPLING_VAR
#endif
