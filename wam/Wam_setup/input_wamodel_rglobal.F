 &NALINE
 CLHEADER=" GLOBAL WAVE MODEL RUN                             MODEL VERSION 4.0 
 ",CBPLTDT="${begofrn}",CEPLTDT="${endofrn}",CDATEF="${begoffo}",
 IDELPRO= 3600, IDELT= 1200,IDELWO= 21600,IDELWI= 21600
 CLMTSU="S","S","S","S"
 IDELALT= 21600,IDELINT= 21600,IDELINS= 21600,IDELAGI= 21600,IDELSPT= 21600
 IDELSPS= 21600,IDELAGS= 21600,IDELRES= 86400
 CLOTSU="S","S","S","S","S","S","S"
 CDATER="${outofrf}",CDATES="${outof2d}"
 PFLAG= F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F
 FFLAG= F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F
 GFLAG= T, T, T, F, F, T, T, F, T, T, T, T, T, T, T, T, F, F, F, F, F
 LFDB= T,LPBIOIN= T,LGRIBIN= T,LPBIOOUT= ${lgribin},
 LFDBIOIN=F,
 LFDBIOOUT=T,
 ICASE= 1,ISHALLO= 1,IREFRA= 0,ITEST= 0,ITESTB= 1,IREST= 1,
 IASSI= $AS,
 IFORCA= 1,IBOUNC= 0,IBOUNF= 0,USERID="   ",RUNID="GLO"
 PATH="                                                            ",
 YCLASS="$CLASS"
 YEXPVER="$XID",
 NENSFNB= $NENS,
 NTOTENS= $TNE,
 LMESSPASS= T,LWCOU= F,NTASKS= 0,NSIZE= 0,
 CPATH="$STORAGE_PATH",
 CFDBSF   = "$cfdbsf",
 CFDB2DSP = "$cfdb2dsp"
 /
