cat >POSTPRO.f <<\EOFPRO
#include "climate.h"
#ifdef 
#include "pospro_global_oper.h"
#include "fldplt_global.h"
#include "plotint_global.h"
#endif
#ifdef 
#include "pospro_medite_oper.h"
#include "fldplt_medite.h"
#include "plotint_medite.h"
#endif
EOFPRO

