#
#   PRINTS SPECTRA FROM WAMODEL OUTPUT FILES.
#
#   HEINZ GUNTHER     ECMWF   JANUARY 1990
#   HEINZ GUNTHER     ECMWF   OCTOBER 1990    (MIGRATION TO UNICOS)
#
# QSUB -r PRSPEC
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
#include "pspec.h"
#include "iecf_len.h"
#include "incdate.h"
#include "inspp.h"
#include "openfil_unicos.h"
#include "prspp.h"
#include "uipspec.h"
EOFPRO
#
#=======================================================================
#  COMPILE AND LOAD.
#=======================================================================
#
cf77 PROGRAM.f
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
#include "input_pspec.h"
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
