#
#  PROGRAM BOUNINT INTERPOLATES ARRAYS OF
#  SPECTRA IN TIME. THE PROGRAM IS SET UP TO
#  FETCH AND STORE ALL DATA DYNAMICALLY.
#
# QSUB -r BOUNINT
# QSUB -q high
# QSUB -lT 40
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
cat >bounint.f <<\EOFPRO
#include "bouint.h"
#include "difdate.h"
#include "iecf_len.h"
#include "incdate.h"
#include "intspec.h"
#include "openfil_unicos.h"
#include "rotspec.h"
#include "strspec.h"
#include "uibou.h"
EOFPRO
#
#===============================================================
#  COMPILE AND LOAD.
#===============================================================
#
cf77 bounint.f
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
#===============================================================
#  RUN THE PROGRAM AND FORMAT THE OUTPUT
#===============================================================
#
a.out <<\EOD > tempout
#include "input_bouint.h"
EOD
#
#===============================================================
#  FORMAT PROGRAM OUTPUT.
#===============================================================
#
asa < tempout
#
#===============================================================
#  END OF JOB
#===============================================================
#
