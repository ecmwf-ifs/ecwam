#
#  THIS JOB COMPILES, LOADS AND SAVES THE WAMICE ABSOLUTE ON THE C90
#  AND ECFILE
#
#  ANNE GUILLAUME      ECMWF    NOVEMBER 1992
#
# QSUB -r wamice_tc90
# QSUB -q high
# QSUB -lT 200
# QSUB -eo
# QSUB

set +v
set -a
PATH=$PATH:/tmp/emos_sms/bin:/ec/emos_sms/bin:.
cd $TMPDIR

#======================================================================
#   Define variables
#======================================================================

EXPER=glo3ice
USEID=sth

#======================================================================

EXPPATH=/$USEID/$EXPER
YMPPATH=$TEMP/$EXPER

# ecfile -p /$USEID/$EXPER add
# mkdir $TEMP/$EXPER

#=======================================================================
#   Copy the source code to the compiler input.
#=======================================================================

cat > wamice.n <<\EOFPRO
#include "readice.h"
#include "readt.h"
#include "micep.h"
#include "msort.h"
#include "adjust.h"
#include "inmarsb.h"
#include "findb.h"
#include "outpp.h"
#include "cormars.h"
#include "readpre.h"
EOFPRO

   cft77 -eszf -Rabc wamice.f || { asa wamice.l ; exit ; }
#  cft77 -eszf       wamice.f || { asa wamice.l ; exit ; }
   segldr wamice.o              \
     -l /usr/lib/libsort.a      \
        $EMOSLIB                \
        $ECLIB                  \
     -V -g -o $YMPPATH/wamicea

# ===================================================================
#  Check return code.
# ===================================================================

   if [ $? != 0 ] ; then
      echo "\n\n\t\t LOADING OF WAMICE  FAILED\n\n"
   fi
   ecfile -p $EXPPATH/wamice.o -hO -uD save wamice.o
   asa wamice.l

#=======================================================================
# END OF JOB.
#=======================================================================


