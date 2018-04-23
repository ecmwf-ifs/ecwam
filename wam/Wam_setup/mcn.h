cat > MAKE_CFS_NODE <<\EOF
#! /bin/ksh
set +xv
[ $DEBUG ] && set -xv

#
#  Create a new node in CFS even if the preceding node dos not exist.
#  ------------------------------------------------------------------

if [ $# -ne 1 ] ; then
 echo '\n\t'This script alows only one argument but there were $#
 exit 1
fi

err=$(pwd)/err

NODE=$1

ecfile -p $NODE -o,o -da list 1> ${err} 2> /dev/null
grep -i directory ${err} > /dev/null && exit

MAKE_CFS_NODE `dirname $NODE` || exit 0

ecfile -p $NODE -o,o -da list 1> ${err} 2> /dev/null
grep -i tape ${err} > /dev/null && \
  {
   [ $TALK ] || echo '\n'
   echo '\t'"$1 exists but is not a node"
   exit 1
  }

grep -i disk ${err} > /dev/null && \
  { 
   [ $TALK ] || echo '\n'
   echo '\t'"$1 exists but is not a node"
   exit  1
  }

getlog=/dev/null
[ $TALK ] && getlog=$(pwd)/getlog

ecfile -p $NODE -da add > ${getlog} || exit 1
[ $TALK ] && { cat $getlog ; \rm $getlog ; } || true
EOF
chmod u+x MAKE_CFS_NODE
