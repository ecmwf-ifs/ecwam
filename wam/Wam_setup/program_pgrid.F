#
#   PRINTS FIELDS OF INTEGRATED DATA.
#
#   HEINZ GUNTHER     ECMWF   JANUARY 1990
#   HEINZ GUNTHER     ECMWF   OCTOBER 1990    (MIGRATION TO UNICOS)
#
# QSUB -r PRMAP
# QSUB -lT 10
# QSUB -eo
# QSUB
#
PATH=$PATH:.
cd $TMPDIR
set +v
#
#=======================================================================
#   Copy the source code to the compiler input.
#=======================================================================
#
cat >PROGRAM.f <<\EOFPRO
#include "pgrid.h"
#include "uipgrid.h"
EOFPRO
#
#=======================================================================
#  COMPILE AND LOAD.
#=======================================================================
#
cf77 -l $HOME/wamoper/wamcrlib2 PROGRAM.f
#
#  Check return code.
#
if [ $? != 0 ]
then
   echo "\n\n\t\tCOMPILE FAILED\n\n"
   echo "\n\n\t\tPROGRAM WILL TERMINATE\n\n"
   exit 3
fi
#
#=======================================================================
#  RUN THE PROGRAM.
#=======================================================================
#
a.out <<\EOD > tempout
#include "input_pgrid.h"
EOD
#
#=======================================================================
#  FORMAT OUTPUT and copy to stdout.
#=======================================================================
#
asa < tempout
#
#=======================================================================
#  END OF JOB
#=======================================================================
#
