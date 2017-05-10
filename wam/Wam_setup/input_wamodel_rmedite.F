 &NALINE
 CLHEADER=" MEDITERRANEAN WAVE MODEL RUN                      MODEL VERSION 4.",
 CBPLTDT="${begofrn}",
 CEPLTDT="${endofrn}",
 CDATEF="${begoffo}",
 IDELPRO= 1200
 IDELT= 1200,
 IDELWO= 21600,
 IDELWI= 21600
 CLMTSU="S", "S", "S", "S"
 IDELINT= 10800,
 IDELINS= 10800,
 IDELSPT= 10800
 IDELSPS= 10800,
 IDELRES= 86400
 CLOTSU="S", "S", "S", "S", "S", "S", "S"
 CDATER="${outofrf}",
 CDATES="${outof2d}"
 PFLAG= F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
F, F, F, F, F, F, F, F, F, F, F, F
 FFLAG= F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
F, F, F, F, F, F, F, F, F, F, F, F
 GFLAG= T, T, T, F, F, T, T, F, F, T, F, F, F, F, F, F, F, F, F, F, F, F, F, F,
F, F, F, F, F, F, F, F, F, F, F, F
 LFDB= ${lfdb},
 CFDBSF   ="${FDB_ROOT}",
 CFDB2DSP ="${FDB_ROOT}",
 LGRIBIN= ${lgribin},
 LPBIOIN= T,
 LGRIBOUT= ${lgribout},
 LPBIOOUT= T,
 LFDBIOOUT=${lfdbioout},
 LNOCDIN=${lnocdin},
 ICASE= 1,
 ISHALLO= 1,
 ITEST= 0,
 ITESTB= 1,
 IREST= 1,
 IASSI= ${iassi},
 LALTAS=${laltas},
 LSARINV=${lsarinv},
 LSARAS=${lsaras},
 YEXPVER="$expver",
 NWAM_BLKS=32,
 LL1D=F,
 CPATH="$STORAGE_PATH"
 /
C  &NAOT CLOUT="" /
