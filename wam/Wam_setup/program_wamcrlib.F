#ifdef JOB
#.
#. CREATES THE WAMCRLIB
#QSUB -q vpp700.normal
#QSUB -eo
#QSUB -r wamcrlib_pre_release
#QSUB -lt 200
#QSUB -lP 0
#QSUB -lV 20mb
#QSUB -lM 64mb
#QSUB -mb -me -mi
#QSUB

# MAKE USE OF MEMORY RESIDENT FILE SYSTEM for TMPDIR
TMPDIR=$MRFSDIR

#.
PATH=$PATH:.
cd $TMPDIR
set -v
set -e
#.
#include "path.h"
#.
cat >wamcrlib.src <<\EOFPRO
#endif
#include "abort1.F"
#include "adjust.F"
#include "aki.F"
#include "chesig.F"
#include "cormars.F"
#include "crewfn.F"
#include "difdate.F"
#include "distc.F"
#include "distn.F"
#include "distpi.F"
#include "distpr.F"
#include "dists.F"
#include "ersfile.F"
#include "expand_string.F"
#include "expoint.F"
#include "extra.F"
#include "fld2fdb.F"
#include "getrest.F"
#include "gsfile_new.F"
#include "iecf_len.F"
#include "incdate.F"
#include "ingrid.F"
#include "ingrids.F"
#include "inmars.F"
#include "inmarss.F"
#include "inspp.F"
#include "inspps.F"
#include "intspec.F"
#include "ispoff.F"
#include "maskl.F"
#include "moment.F"
#include "mpclose_unit.F"
#include "outpp.F"
#include "outnam.F"
#include "packd.F"
#include "prinsta.F"
#include "prspp.F"
#include "prspps.F"
#include "readsta.F"
#include "rotate.F"
#include "rotspec.F"
#include "savrest.F"
#include "sortc.F"
#include "sorti.F"
#include "sortini.F"
#include "sortinr.F"
#include "sortr.F"
#include "stats.F"
#include "statsi.F"
#include "statsc.F"
#include "statse.F"
#include "statsp.F"
#include "strspec.F"
#include "ustru10.F"
#include "wam_user_clock.F"
#include "writsta.F"
#include "myindex.F"
#include "dummy_cray.F"
#include "uraopc.F"
#include "makegrids.F"
#include "gribpacs.F"
#include "i_get_unit.F"
#ifdef JOB
EOFPRO
#.
#. COMPILE SOURCE CODE
#.
/home/rd/sth/bin/f90split < wamcrlib.src
if [ -f WAMINFOT.f -a -f WRITSTA.f ] ; then
  cat WAMINFOT.f >> WRITSTA.f
fi
#.
[ ! -d tmp ] && mkdir tmp
[ -f ${VPPPATH}/lib$WAMCRLIB.f.tar ] && {
  dd if=${VPPPATH}/lib$WAMCRLIB.f.tar of=lib$WAMCRLIB.f.tar bs=400k
  cd tmp
  tar xf ../lib$WAMCRLIB.f.tar
  cd ..
}
#.
[ -f ${VPPPATH}/lib$WAMCRLIB.a ] && {
  dd if=${VPPPATH}/lib$WAMCRLIB.a of=lib$WAMCRLIB.a bs=400k
}
#.
i=0
export file WAMCRLIB
for file in $(ls -1 *.f) ; do
  (
  if [ "$REMAKE" != YES ] ; then
    if [ -f tmp/$file ] ; then
      if cmp $file tmp/$file ; then
        if ar tv lib$WAMCRLIB.a | grep ${file%.f}.o ; then
          exit
        fi
      fi
    fi
  fi
  frt $rp -Of -M. -Wv,-Of -Pos -Z ${file%.f}.l -X9 \
      -c -o ${file%.f}.o ${file}
  )&
  (( i += 1 ))
  [[ $i -gt 5 ]] && {
    i=0
    wait
  }
done
wait
print - "compile done"
#.
#. BUILD LIBRARY
#.
if [[ $(ls -1 *.o | wc -l ) -gt 0 ]] ; then
  ar rv lib$WAMCRLIB.a *.o 
fi
#.
tar cf lib$WAMCRLIB.l.tar *.l
tar cf lib$WAMCRLIB.f.tar *.f
#.
#. SAVE LIBRARY
#.
ecp -o lib$WAMCRLIB.* ec:$DHSPATH/.
#.
ls -la
#.
dd if=lib$WAMCRLIB.a     of=$VPPPATH/lib$WAMCRLIB.a       bs=400k
dd if=lib$WAMCRLIB.f.tar of=${VPPPATH}/lib$WAMCRLIB.f.tar bs=400k
dd if=lib$WAMCRLIB.l.tar of=${VPPPATH}/lib$WAMCRLIB.l.tar bs=400k
#.
chmod go+rx $VPPPATH/lib$WAMCRLIB.*
#.
cat *.l
#.
#.END OF JOB
#.
#endif
