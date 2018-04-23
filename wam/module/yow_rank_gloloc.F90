MODULE YOW_RANK_GLOLOC
#if defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
  USE coupling_var, ONLY : MyRankGlobal, MyRankLocal
#else
  integer :: MyRankGlobal, MyRankLocal
#endif
END MODULE
